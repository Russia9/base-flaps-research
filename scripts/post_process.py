#!/usr/bin/env python3
"""Convert OpenFOAM force/moment logs into body-axis coefficient CSV.

Usage:
    python scripts/post_process.py <case-dir>

Reads:
    <case>/constant/freestreamProperties     (pInf, TInf, UInfMag, RGas)
    newest valid <case>/postProcessing/forces/<time>/force.dat
    newest valid <case>/postProcessing/forces/<time>/moment.dat

Writes:
    results/<case-name>/coefficients.csv   (normalized by q_inf*S and q_inf*S*D)
    results/<case-name>/forces.csv          (raw force/moment, N and N*m)

Coefficients use D = 0.08 m as reference length and S = pi*D^2/4 as
reference area; moments are referenced to the nose tip (origin).
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

D_REF = 0.08
S_REF = math.pi * D_REF * D_REF / 4.0


def read_scalar(dict_text: str, name: str) -> float:
    m = re.search(rf"^\s*{name}\s+([-+0-9.eE]+)\s*;", dict_text, re.MULTILINE)
    if not m:
        raise ValueError(f"Could not find scalar '{name}' in freestreamProperties")
    return float(m.group(1))


def parse_freestream(case: Path) -> dict[str, float]:
    text = (case / "constant" / "freestreamProperties").read_text()
    p_inf = read_scalar(text, "pInf")
    t_inf = read_scalar(text, "TInf")
    u_mag = read_scalar(text, "UInfMag")
    r_gas = read_scalar(text, "RGas")
    rho = p_inf / (r_gas * t_inf)
    q = 0.5 * rho * u_mag * u_mag
    return {"pInf": p_inf, "TInf": t_inf, "UInfMag": u_mag, "RGas": r_gas,
            "rhoInf": rho, "qInf": q}


def parse_force_file(path: Path) -> list[dict[str, float]]:
    """Parse OpenFOAM v2512 force.dat / moment.dat.

    Each non-comment line is whitespace-separated columns (no parentheses):
        <time>  total(x y z)  pressure(x y z)  viscous(x y z)  [porous(x y z)]
    where total = pressure + viscous (+ porous). We keep the reported total
    and the pressure/viscous split; the template never enables porosity.
    """
    rows: list[dict[str, float]] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        nums = [float(tok) for tok in line.split()]
        # time + at least total/pressure/viscous triplets.
        if len(nums) < 10:
            raise ValueError(f"Unexpected line in {path}: {raw!r}")
        rows.append({
            "time": nums[0],
            "x": nums[1], "y": nums[2], "z": nums[3],
            "px": nums[4], "py": nums[5], "pz": nums[6],
            "vx": nums[7], "vy": nums[8], "vz": nums[9],
        })
    return rows


def find_forces_dir(case: Path) -> Path:
    base = case / "postProcessing" / "forces"
    if not base.is_dir():
        raise FileNotFoundError(f"No postProcessing/forces in {case}")
    valid: list[tuple[float, Path]] = []
    for sub in base.iterdir():
        if not sub.is_dir():
            continue
        try:
            start_time = float(sub.name)
        except ValueError:
            continue
        if (sub / "force.dat").is_file() and (sub / "moment.dat").is_file():
            valid.append((start_time, sub))
    if not valid:
        raise FileNotFoundError(f"No complete force.dat/moment.dat pair in {base}")
    return max(valid, key=lambda item: item[0])[1]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2

    case = Path(argv[1]).resolve()
    fs = parse_freestream(case)
    q_s = fs["qInf"] * S_REF
    q_s_d = q_s * D_REF

    forces_dir = find_forces_dir(case)
    forces = parse_force_file(forces_dir / "force.dat")
    moments = parse_force_file(forces_dir / "moment.dat")
    if not forces:
        raise ValueError(f"No force samples found in {forces_dir / 'force.dat'}")
    if not moments:
        raise ValueError(f"No moment samples found in {forces_dir / 'moment.dat'}")
    if len(forces) != len(moments):
        raise ValueError("force.dat and moment.dat sample counts differ")

    out_dir = Path("results") / case.name
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "coefficients.csv"

    cols = ["time",
            "Cx", "Cy", "Cz", "Mx", "My", "Mz",
            "Cx_p", "Cy_p", "Cz_p", "Cx_v", "Cy_v", "Cz_v",
            "Mx_p", "My_p", "Mz_p", "Mx_v", "My_v", "Mz_v"]
    with out.open("w") as fh:
        fh.write(",".join(cols) + "\n")
        for f, m in zip(forces, moments):
            row = [
                f["time"],
                f["x"] / q_s, f["y"] / q_s, f["z"] / q_s,
                m["x"] / q_s_d, m["y"] / q_s_d, m["z"] / q_s_d,
                f["px"] / q_s, f["py"] / q_s, f["pz"] / q_s,
                f["vx"] / q_s, f["vy"] / q_s, f["vz"] / q_s,
                m["px"] / q_s_d, m["py"] / q_s_d, m["pz"] / q_s_d,
                m["vx"] / q_s_d, m["vy"] / q_s_d, m["vz"] / q_s_d,
            ]
            fh.write(",".join(f"{v:.8e}" for v in row) + "\n")

    # Raw dimensional forces/moments, straight from force.dat/moment.dat
    # (total plus pressure/viscous split), N and N*m. Same layout as the
    # coefficient CSV, just undivided.
    out_forces = out_dir / "forces.csv"
    force_cols = ["time",
                  "Fx", "Fy", "Fz", "Mx", "My", "Mz",
                  "Fx_p", "Fy_p", "Fz_p", "Fx_v", "Fy_v", "Fz_v",
                  "Mx_p", "My_p", "Mz_p", "Mx_v", "My_v", "Mz_v"]
    with out_forces.open("w") as fh:
        fh.write(",".join(force_cols) + "\n")
        for f, m in zip(forces, moments):
            row = [
                f["time"],
                f["x"], f["y"], f["z"],
                m["x"], m["y"], m["z"],
                f["px"], f["py"], f["pz"],
                f["vx"], f["vy"], f["vz"],
                m["px"], m["py"], m["pz"],
                m["vx"], m["vy"], m["vz"],
            ]
            fh.write(",".join(f"{v:.8e}" for v in row) + "\n")

    # Mean over the last 10% of samples for a quick eyeball.
    n = len(forces)
    tail = max(1, n // 10)
    def mean(rows, key, denom):
        return sum(r[key] for r in rows[-tail:]) / tail / denom

    print(f"case          : {case}")
    print(f"samples       : {n} (mean over last {tail})")
    print(f"freestream    : pInf={fs['pInf']:.0f} Pa  TInf={fs['TInf']:.2f} K  "
          f"UInfMag={fs['UInfMag']:.1f} m/s  rhoInf={fs['rhoInf']:.4f} kg/m^3")
    print(f"qInf          : {fs['qInf']:.1f} Pa")
    print(f"D, S          : {D_REF} m, {S_REF:.4e} m^2")
    print(f"Cx, Cy, Cz    : {mean(forces, 'x', q_s):+.4f}, "
          f"{mean(forces, 'y', q_s):+.4f}, {mean(forces, 'z', q_s):+.4f}")
    print(f"Mx, My, Mz    : {mean(moments, 'x', q_s_d):+.4f}, "
          f"{mean(moments, 'y', q_s_d):+.4f}, {mean(moments, 'z', q_s_d):+.4f}")
    print(f"output        : {out}")
    print(f"                {out_forces}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
