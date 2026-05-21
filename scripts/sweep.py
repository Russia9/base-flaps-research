#!/usr/bin/env python3
"""Enumerate or create the 108-case aft-base flap sweep."""

from __future__ import annotations

import argparse
import itertools
import shlex
import sys
from pathlib import Path

from create_case import ROOT, create_case, format_number

N_VALUES = [1, 2, 3, 4]
XI_VALUES = [30.0, 45.0, 90.0]
LD_VALUES = [0.5, 1.0, 1.5]
MACH_VALUES = [1.5, 2.0, 2.5]


def token(value: float) -> str:
    return format_number(value).replace(".", "p")


def case_name(n_fins: int, xi: float, ld: float, mach: float) -> str:
    return f"N{n_fins}_xi{token(xi)}_LD{token(ld)}_M{token(mach)}"


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("openfoam/cases"), help="directory for generated cases")
    p.add_argument("--TD", type=float, default=0.02, help="fin thickness divided by D")
    p.add_argument("--D", type=float, default=80.0, help="body diameter in mm")
    p.add_argument("--gamma", type=float, default=1.4, help="specific heat ratio for UInf calculation")
    p.add_argument("--dry-run", action="store_true", help="print commands without creating cases")
    p.add_argument("--force", action="store_true", help="replace existing cases when creating")
    return p


def shell_join(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def main(argv: list[str]) -> int:
    args = parser().parse_args(argv)
    cases = list(itertools.product(N_VALUES, XI_VALUES, LD_VALUES, MACH_VALUES))

    if args.dry_run:
        for n_fins, xi, ld, mach in cases:
            case = args.root / case_name(n_fins, xi, ld, mach)
            create_cmd = [
                "python3", "scripts/create_case.py", "--force",
                "--case", str(case),
                "--N", str(n_fins),
                "--xi", format_number(xi),
                "--LD", format_number(ld),
                "--TD", format_number(args.TD),
                "--Mach", format_number(mach),
            ]
            print(shell_join(create_cmd))
            print(shell_join(["./rebuild-mesh.sh", str(case)]))
            print(shell_join(["./run-simulation.sh", "--dry-run", str(case)]))
            print(shell_join(["./run-simulation.sh", str(case)]))
        print(f"# total cases: {len(cases)}", file=sys.stderr)
        return 0

    for n_fins, xi, ld, mach in cases:
        case = args.root / case_name(n_fins, xi, ld, mach)
        u_mag = create_case(
            ROOT / case,
            diameter_mm=args.D,
            n_fins=n_fins,
            xi=xi,
            ld=ld,
            td=args.TD,
            mach=mach,
            gamma=args.gamma,
            force=args.force,
        )
        print(f"{case}: UInfMag={format_number(u_mag)}")

    print(f"created cases: {len(cases)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
