#!/usr/bin/env python3
# /// script
# dependencies = ["gmsh"]
# ///
"""Build the zoned geometry for the bare ogive-cylinder body (OpenCASCADE).

Usage:
    uv run scripts/gen_gmsh.py openfoam/<case> [--gui]

Step 1 of the zoned-mesh plan: geometry only, no mesh. Reads D and N from
<case>/constant/caseProperties and the zone layout from <case>/system/gmshParams,
builds the body and the domain, splits the fluid into zones, prints checks
against the analytic geometry, and writes <case>/geometry.brep.

Zones (each meshed on its own later and stitched in OpenFOAM, so a zone's
cell counts never leak into another): nose, body and aft inside the cylinder
r <= rZone, far outside it. Seams are the planes x = xNoseEnd, x = xAftStart
and the boundary of the inner cylinder.

Geometry follows geometry/arc_stabilizers.scad (secant ogive, 10 D body)
except the nose, which is sharp as in the report.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import gmsh
from create_case import parse_scalar

PARAMS = {
    "xInlet": -800.0, "xOutlet": 5600.0, "rFar": 600.0,
    "xZoneStart": -40.0, "xNoseEnd": 40.0, "xAftStart": 620.0, "xZoneEnd": 1200.0,
    "rZone": 120.0,
}


def read_params(case: Path) -> dict[str, float]:
    text = (case / "system" / "gmshParams").read_text()
    params = {}
    for name, default in PARAMS.items():
        try:
            params[name] = parse_scalar(text, name)
        except ValueError:
            params[name] = default
    return params


class Body:
    """Secant ogive through (0, 0) and (2D, D/2), radius 8.5 D; cylinder to 10 D."""

    def __init__(self, D: float):
        self.R, self.rho = D / 2, 8.5 * D
        self.l_ogive, self.total = 2.0 * D, 10.0 * D
        k = (self.l_ogive**2 + self.R**2) / (2 * self.l_ogive)
        m = self.R / self.l_ogive
        self.yc = (k * m - math.sqrt(k**2 * m**2 - (1 + m**2) * (k**2 - self.rho**2))) / (1 + m**2)
        self.xc = k - m * self.yc

    def r(self, x: float) -> float:
        if x <= 0:
            return 0.0
        if x >= self.l_ogive:
            return self.R
        return self.yc + math.sqrt(self.rho**2 - (x - self.xc) ** 2)

    def slope_deg(self, x: float) -> float:
        return math.degrees(math.atan((self.xc - x) / (self.r(x) - self.yc)))

    def integrate(self, x0: float, x1: float, n: int = 20000):
        """Volume and lateral area of the body between x0 and x1 (Simpson)."""
        h = (x1 - x0) / n
        vol = area = 0.0
        for i in range(n + 1):
            x = x0 + i * h
            w = 1 if i in (0, n) else (4 if i % 2 else 2)
            r = self.r(x)
            drdx = 0.0 if x >= self.l_ogive else math.tan(math.radians(self.slope_deg(x)))
            vol += w * math.pi * r * r
            area += w * 2 * math.pi * r * math.sqrt(1 + drdx * drdx)
        return vol * h / 3, area * h / 3


def build(body: Body, P: dict) -> dict:
    occ = gmsh.model.occ
    R, total = body.R, body.total

    # Body: revolve the closed profile (apex, ogive arc, cylinder, base, axis).
    p_apex = occ.addPoint(0, 0, 0)
    p_sh = occ.addPoint(body.l_ogive, R, 0)
    p_rim = occ.addPoint(total, R, 0)
    p_base = occ.addPoint(total, 0, 0)
    p_c = occ.addPoint(body.xc, body.yc, 0)
    loop = occ.addCurveLoop([occ.addCircleArc(p_apex, p_c, p_sh), occ.addLine(p_sh, p_rim),
                             occ.addLine(p_rim, p_base), occ.addLine(p_base, p_apex)])
    profile = occ.addPlaneSurface([loop])
    solid = [t for d, t in occ.revolve([(2, profile)], 0, 0, 0, 1, 0, 0, 2 * math.pi) if d == 3]
    occ.remove([(0, p_c)])

    # Domain and the three inner zone cylinders.
    x_in, x_out, r_ff, r_z = P["xInlet"], P["xOutlet"], P["rFar"], P["rZone"]
    xs = [P["xZoneStart"], P["xNoseEnd"], P["xAftStart"], P["xZoneEnd"]]
    domain = occ.addCylinder(x_in, 0, 0, x_out - x_in, 0, 0, r_ff)
    zones = [occ.addCylinder(xs[i], 0, 0, xs[i + 1] - xs[i], 0, 0, r_z) for i in range(3)]
    # Cut the body out of the domain and out of every zone, then fragment so
    # neighbouring zones share their seam faces exactly.
    tool = [(3, s) for s in solid]
    fluid, _ = occ.cut([(3, domain)], tool, removeTool=False)
    inner, _ = occ.cut([(3, z) for z in zones], tool, removeTool=False)
    occ.fragment(fluid, inner)
    occ.remove(tool, recursive=True)
    occ.synchronize()

    # Name the zones by position (OCC pads bounding boxes slightly).
    names = {}
    pad = 1e-3
    for _, v in gmsh.model.getEntities(3):
        bb = gmsh.model.getBoundingBox(3, v)
        if bb[4] < r_z + pad and bb[0] > xs[0] - pad and bb[3] < xs[3] + pad:
            mid = 0.5 * (bb[0] + bb[3])
            names[v] = "nose" if mid < xs[1] else ("body" if mid < xs[2] else "aft")
        else:
            names[v] = "far"
    if sorted(names.values()) != ["aft", "body", "far", "nose"]:
        raise RuntimeError(f"unexpected zone split: {names}")

    # Classify boundary faces: wall, outer patches, seams between zones.
    owners: dict[int, list[str]] = {}
    for v, name in names.items():
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            owners.setdefault(s, []).append(name)
    groups: dict[str, list[int]] = {}
    for s, own in owners.items():
        bb = gmsh.model.getBoundingBox(2, s)
        if len(own) == 2:
            key = "seam_" + "_".join(sorted(own))
        elif bb[3] < x_in + pad:
            key = "inlet"
        elif bb[0] > x_out - pad:
            key = "outlet"
        elif bb[4] > r_ff - pad and bb[1] < -r_ff + pad:
            key = "farfield"
        else:
            key = "fuselage"
        groups.setdefault(key, []).append(s)

    for v, name in names.items():
        gmsh.model.addPhysicalGroup(3, [v], name=name)
    for key, surfs in groups.items():
        gmsh.model.addPhysicalGroup(2, surfs, name=key)
    return {"zones": names, "groups": groups, "xs": xs}


def check(body: Body, P: dict, info: dict) -> list[str]:
    """Compare the OCC model against the analytic body and zone volumes."""
    occ = gmsh.model.occ
    lines, ok = [], True
    R, total, r_z = body.R, body.total, P["rZone"]
    xs = info["xs"]

    def row(label, got, want, tol_rel=1e-5):
        nonlocal ok
        err = abs(got - want) / abs(want)
        ok &= err < tol_rel
        lines.append(f"  {label:34s} {got:16.6f} {want:16.6f}  {err:.1e} {'ok' if err < tol_rel else 'FAIL'}")

    lines.append(f"  {'quantity':34s} {'OCC':>16s} {'analytic':>16s}  rel.err")
    v_body, a_lat = body.integrate(0.0, total)
    a_wall = sum(occ.getMass(2, s) for s in info["groups"]["fuselage"])
    row("body wall area (mm^2)", a_wall, a_lat + math.pi * R * R)
    v_dom = math.pi * P["rFar"] ** 2 * (P["xOutlet"] - P["xInlet"])
    v_fluid = sum(occ.getMass(3, v) for v in info["zones"])
    row("fluid volume (mm^3)", v_fluid, v_dom - v_body)
    for name, (x0, x1) in (("nose", xs[0:2]), ("body", xs[1:3]), ("aft", xs[2:4])):
        vb, _ = body.integrate(max(0.0, x0), min(total, x1))
        v = next(v for v, n in info["zones"].items() if n == name)
        row(f"zone {name} volume (mm^3)", occ.getMass(3, v), math.pi * r_z**2 * (x1 - x0) - vb)

    # Wall radius against the analytic profile, sampled on each lateral wall
    # face over its own parameter range (the base disk is checked by area).
    worst, n = 0.0, 0
    for s in info["groups"]["fuselage"]:
        if gmsh.model.getType(2, s) == "Plane":
            continue
        (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, s)
        uv = [c for i in range(21) for j in range(21)
              for c in (u0 + (u1 - u0) * i / 20, v0 + (v1 - v0) * j / 20)]
        pts = gmsh.model.getValue(2, s, uv)
        for k in range(0, len(pts), 3):
            x, y, z = pts[k:k + 3]
            worst = max(worst, abs(math.hypot(y, z) - body.r(x)))
            n += 1
    ok &= worst < 1e-6
    lines.append(f"  wall radius vs profile, {n} points: max |OCC - analytic| = {worst:.2e} mm  "
                 f"{'ok' if worst < 1e-6 else 'FAIL'}")
    lines.append(f"  tip half-angle {body.slope_deg(0.0):.3f} deg, shoulder kink "
                 f"{body.slope_deg(body.l_ogive - 1e-9):.3f} deg, ogive centre "
                 f"({body.xc:.4f}, {body.yc:.4f}) mm")
    lines.append("  " + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    return lines


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--gui", action="store_true", help="open the geometry in the Gmsh GUI")
    args = ap.parse_args(argv)
    case = args.case
    try:
        props = (case / "constant" / "caseProperties").read_text()
        D = parse_scalar(props, "D")
        N = int(parse_scalar(props, "N"))
        if N != 0:
            raise ValueError(f"N = {N}: only the bare body (N = 0) is supported so far")
        params = read_params(case)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    body = Body(D)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        info = build(body, params)
        gmsh.write(str(case / "geometry.brep"))
        print("zones:")
        for v, name in sorted(info["zones"].items(), key=lambda t: t[1]):
            print(f"  {name:5s} volume {v}, {len(gmsh.model.getBoundary([(3, v)], oriented=False))} faces")
        print("face groups:")
        for key, surfs in sorted(info["groups"].items()):
            print(f"  {key:18s} {len(surfs)} faces")
        print("checks:")
        print("\n".join(check(body, params, info)))
        if args.gui:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
