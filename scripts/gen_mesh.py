#!/usr/bin/env python3
# /// script
# dependencies = ["gmsh"]
# ///
"""Zoned structured Gmsh mesh of the ogive-cylinder body, bare or with arc fins.

Usage:
    uv run scripts/gen_mesh.py openfoam/<case> [--geometry-only] [--gui]
    uv run scripts/gen_mesh.py openfoam/<case> --fins

Reads D, N, xi and L from <case>/constant/caseProperties and the layout
from <case>/system/gmshParams (see scripts/gmesh). --geometry-only builds
the block layout, runs the checks and writes a coarse preview mesh to
<case>/preview (case.foam for ParaView). --fins builds only the fins and
compares them with the scad's STL.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gmsh
from create_case import parse_scalar
from gmesh import checks, export, layout, occ, params
from gmesh.shapes import Body, FinSpec


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--gui", action="store_true", help="open the result in the Gmsh GUI")
    ap.add_argument("--geometry-only", action="store_true",
                    help="stop after the block layout and write a coarse preview mesh to <case>/preview")
    ap.add_argument("--preview-cells", type=int, default=4, metavar="N",
                    help="cells along every block edge in the preview mesh (default 4)")
    ap.add_argument("--fins", action="store_true", help="build only the fins and compare them with the scad")
    args = ap.parse_args(argv)
    case = args.case
    try:
        props = (case / "constant" / "caseProperties").read_text()
        D = parse_scalar(props, "D")
        N = int(parse_scalar(props, "N"))
        P = params.read(case)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    body = Body(D)
    spec = FinSpec(body, N, parse_scalar(props, "xi"), parse_scalar(props, "L")) if N else None

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        if args.fins:
            if not spec:
                print("error: N = 0, no fins", file=sys.stderr)
                return 1
            fins = occ.fin_solids(spec)
            print("\n".join(checks.fin_check(spec, fins, case, D)))
            return 0
        info = layout.build(body, P, spec)
        gmsh.write(str(case / "geometry.brep"))
        t = info["tip"]
        print(f"nose core: tip half-angle {t['beta']:.2f} deg, patch edge x {t['x_a']:.3f}-{t['x_45']:.3f} mm, "
              f"half-width {t['h_in']:.2f} mm at the zone start")
        if spec:
            print("\n".join(info["layout"].fins.summary()))
        print("face groups:")
        for key, surfs in sorted(info["groups"].items()):
            if not key.startswith("shared_"):
                print(f"  {key:28s} {len(surfs)} faces")
        print("checks:")
        print("\n".join(checks.check(body, P, info)))
        if args.geometry_only:
            print(f"preview mesh: {export.preview_mesh(info, case / 'preview', args.preview_cells)}")
        if args.gui:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
