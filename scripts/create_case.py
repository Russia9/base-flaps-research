#!/usr/bin/env python3
"""Create one OpenFOAM case by cloning an existing case and retargeting Mach.

Example:
    python3 scripts/create_case.py --force --case openfoam/arc_M2p0 \
        --N 4 --xi 90 --L 140 --Mach 2.0

The source case (--template, default openfoam/arc) supplies the mesh
dictionaries and boundary conditions; this script rewrites
constant/caseProperties and the Mach-dependent entries of
constant/freestreamProperties. Pass --N 0 for the clean-body baseline.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TEMPLATE = ROOT / "openfoam" / "arc"


def parse_scalar(text: str, name: str) -> float:
    match = re.search(rf"^\s*{name}\s+([-+0-9.eE]+)\s*;", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find scalar '{name}'")
    return float(match.group(1))


def format_number(value: float) -> str:
    return f"{value:.8g}"


def replace_scalar(text: str, name: str, value: float) -> str:
    pattern = rf"(^\s*{name}\s+)([^;]+)(;.*$)"
    replacement = rf"\g<1>{format_number(value)}\3"
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"Could not replace scalar '{name}'")
    return new_text


def replace_vector(text: str, name: str, values: tuple[float, float, float]) -> str:
    vector = f"({format_number(values[0])} {format_number(values[1])} {format_number(values[2])})"
    pattern = rf"(^\s*{name}\s+)\([^;]+\)(;.*$)"
    replacement = rf"\g<1>{vector}\2"
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"Could not replace vector '{name}'")
    return new_text


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def write_case_properties(
    case: Path,
    *,
    diameter_mm: float,
    n_fins: int,
    xi: float,
    chord_mm: float,
    mach: float,
    gamma: float,
) -> None:
    fin_note = "clean-body baseline: no stabilizers" if n_fins == 0 else f"L/D = {format_number(chord_mm / diameter_mm)}"
    text = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2512                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      caseProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

D           {format_number(diameter_mm)};   // mm
N           {n_fins};
xi          {format_number(xi)};
L           {format_number(chord_mm)};  // mm, stabilizer chord ({fin_note})
Mach        {format_number(mach)};
gamma       {format_number(gamma)};


// ************************************************************************* //
"""
    (case / "constant" / "caseProperties").write_text(text)


def update_freestream(case: Path, mach: float, gamma: float) -> float:
    path = case / "constant" / "freestreamProperties"
    text = path.read_text()
    t_inf = parse_scalar(text, "TInf")
    r_gas = parse_scalar(text, "RGas")
    p_inf = parse_scalar(text, "pInf")
    intensity = parse_scalar(text, "I")
    mu_ratio = parse_scalar(text, "muRatio")
    u_mag = mach * math.sqrt(gamma * r_gas * t_inf)

    # Freestream turbulence held at constant intensity and eddy-viscosity ratio
    # across cases, so k and omega scale with Mach. Sutherland viscosity is
    # read from thermophysicalProperties to keep a single source of truth.
    thermo = (case / "constant" / "thermophysicalProperties").read_text()
    a_s = parse_scalar(thermo, "As")
    t_s = parse_scalar(thermo, "Ts")
    mu_inf = a_s * t_inf**1.5 / (t_inf + t_s)
    rho_inf = p_inf / (r_gas * t_inf)
    nu_inf = mu_inf / rho_inf
    k_inf = 1.5 * (intensity * u_mag) ** 2
    omega_inf = k_inf / (mu_ratio * nu_inf)

    text = replace_scalar(text, "UInfMag", u_mag)
    text = replace_vector(text, "UInf", (u_mag, 0.0, 0.0))
    text = replace_scalar(text, "kInf", k_inf)
    text = replace_scalar(text, "omegaInf", omega_inf)
    path.write_text(text)
    return u_mag


def create_case(
    case: Path,
    *,
    template: Path = DEFAULT_TEMPLATE,
    diameter_mm: float,
    n_fins: int,
    xi: float,
    chord_mm: float,
    mach: float,
    gamma: float,
    force: bool,
) -> float:
    case = case.resolve()
    template = template.resolve()
    openfoam_root = (ROOT / "openfoam").resolve()

    if case == template:
        raise ValueError("Refusing to overwrite the source case")
    if n_fins < 0:
        raise ValueError("N must be >= 0")
    if n_fins > 0 and not 0 < xi < 360:
        raise ValueError("xi must be in (0, 360) degrees")

    if case.exists():
        if not force:
            raise FileExistsError(f"{case} already exists; pass --force to replace it")
        if not is_relative_to(case, openfoam_root):
            raise ValueError(f"Refusing to overwrite path outside {openfoam_root}: {case}")
        shutil.rmtree(case)

    case.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(template, case)

    write_case_properties(
        case,
        diameter_mm=diameter_mm,
        n_fins=n_fins,
        xi=xi,
        chord_mm=chord_mm,
        mach=mach,
        gamma=gamma,
    )
    u_mag = update_freestream(case, mach, gamma)
    (case / "case.foam").touch()
    return u_mag


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", required=True, type=Path, help="case directory to create")
    p.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE, help="source case directory to clone")
    p.add_argument("--force", action="store_true", help="replace an existing case under openfoam/")
    p.add_argument("--D", type=float, default=80.0, help="body diameter in mm")
    p.add_argument("--N", type=int, required=True, help="number of fins; 0 for the clean-body baseline")
    p.add_argument("--xi", type=float, default=90.0, help="fin arc angle in degrees")
    p.add_argument("--L", type=float, default=140.0, help="fin chord in mm")
    p.add_argument("--Mach", type=float, required=True, help="freestream Mach number")
    p.add_argument("--gamma", type=float, default=1.4, help="specific heat ratio for UInf calculation")
    return p


def main(argv: list[str]) -> int:
    args = parser().parse_args(argv)
    try:
        u_mag = create_case(
            args.case,
            template=args.template,
            diameter_mm=args.D,
            n_fins=args.N,
            xi=args.xi,
            chord_mm=args.L,
            mach=args.Mach,
            gamma=args.gamma,
            force=args.force,
        )
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"case      : {args.case}")
    print(f"geometry  : D={format_number(args.D)}mm N={args.N} xi={format_number(args.xi)} L={format_number(args.L)}mm")
    print(f"freestream: Mach={format_number(args.Mach)} UInfMag={format_number(u_mag)} m/s")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
