#!/usr/bin/env python3
# /// script
# dependencies = ["gmsh"]
# ///
"""Build the zoned geometry for the bare ogive-cylinder body (OpenCASCADE).

Usage:
    uv run scripts/gen_gmsh.py openfoam/<case> [--gui]

Step 1 of the zoned-mesh plan: geometry and block layout, no mesh. Reads D and N from
<case>/constant/caseProperties and the zone layout from <case>/system/gmshParams,
builds the body and the domain, splits the fluid into zones, prints checks
against the analytic geometry, and writes <case>/geometry.brep.

Zones (each meshed on its own later and stitched in OpenFOAM, so a zone's
cell counts never leak into another): nose, body and aft inside the cylinder
r <= rZone, far outside it. Seams are the planes x = xNoseEnd, x = xAftStart
and the boundary of the inner cylinder.

Every zone is split into hexahedral blocks (8 sectors of 45 degrees; a 2x2
core where the axis is in the fluid) so it can carry a structured mesh. At
the sharp tip the nose core is a square frustum whose downstream faces lie on
the ogive; the apex is an ordinary block corner.

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
    "rZone": 120.0, "coreTip": 0.5, "xTipInterface": -30.0, "coreWake": 20.0, "coreFar": 60.0,
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


SQ2 = math.sqrt(2.0)
PAD = 1e-3   # OCC pads bounding boxes slightly


def half_plane(theta: float, x0: float, x1: float, r0: float, r1: float):
    """Rectangle x0..x1, r0..r1 in the meridional half-plane at angle theta."""
    occ = gmsh.model.occ
    s = occ.addRectangle(x0, r0, 0, x1 - x0, r1 - r0)
    occ.rotate([(2, s)], 0, 0, 0, 1, 0, 0, theta)
    return (2, s)


def annulus(x: float, r0: float, r1: float):
    """Disk (r0 = 0) or annulus in the plane x = const."""
    occ = gmsh.model.occ
    d = occ.addDisk(0, 0, 0, r1, r1)
    if r0 > 0:
        (_, d), = occ.cut([(2, d)], [(2, occ.addDisk(0, 0, 0, r0, r0))])[0]
    occ.rotate([(2, d)], 0, 0, 0, 0, 1, 0, math.pi / 2)
    occ.translate([(2, d)], x, 0, 0)
    return (2, d)


def square_prism(x0: float, h0: float, x1: float, h1: float):
    """Square-section solid along x, half-width h0 at x0 and h1 at x1, with planar faces
    (a frustum if h0 != h1)."""
    occ = gmsh.model.occ
    corners = ((1, 1), (-1, 1), (-1, -1), (1, -1))
    p = [[occ.addPoint(x, sy * h, sz * h) for sy, sz in corners] for x, h in ((x0, h0), (x1, h1))]
    ring = [[occ.addLine(p[e][i], p[e][(i + 1) % 4]) for i in range(4)] for e in range(2)]
    axial = [occ.addLine(p[0][i], p[1][i]) for i in range(4)]
    faces = [occ.addPlaneSurface([occ.addCurveLoop(ring[e])]) for e in range(2)]
    faces += [occ.addPlaneSurface([occ.addCurveLoop([ring[0][i], axial[(i + 1) % 4], -ring[1][i], -axial[i]])])
              for i in range(4)]
    return (3, occ.addVolume([occ.addSurfaceLoop(faces)]))


def polygon_half_plane(theta: float, pts: list[tuple[float, float]]):
    """Planar polygon given as (x, r) points in the meridional half-plane at angle theta."""
    occ = gmsh.model.occ
    p = [occ.addPoint(x, r, 0) for x, r in pts]
    loop = occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))])
    s = occ.addPlaneSurface([loop])
    occ.rotate([(2, s)], 0, 0, 0, 1, 0, 0, theta)
    return (2, s)


def fragment(zones: dict[int, str], tools: list) -> dict[int, str]:
    """Fragment zone volumes with tools; return {new volume: zone}, dropping tool-only pieces."""
    occ = gmsh.model.occ
    objects = list(zones)
    _, children = occ.fragment([(3, v) for v in objects], tools)
    new = {}
    for i, v in enumerate(objects):
        for d, t in children[i]:
            if d == 3:
                new[t] = zones[v]
    occ.synchronize()
    stray = sorted({(3, t) for kids in children[len(objects):] for d, t in kids if d == 3 and t not in new})
    if stray:
        occ.remove(stray)
    occ.synchronize()
    for dim in (2, 1, 0):   # drop leftover tool pieces that bound nothing
        free = [(dim, t) for _, t in gmsh.model.getEntities(dim)
                if len(gmsh.model.getAdjacencies(dim, t)[0]) == 0]
        if free:
            occ.remove(free)
        occ.synchronize()
    return new


def build(body: Body, P: dict) -> dict:
    occ = gmsh.model.occ
    R, total = body.R, body.total
    x_in, x_out, r_ff, r_z = P["xInlet"], P["xOutlet"], P["rFar"], P["rZone"]
    xs = [P["xZoneStart"], P["xNoseEnd"], P["xAftStart"], P["xZoneEnd"]]
    a, s_w, s_f = P["coreTip"], P["coreWake"], P["coreFar"]
    thetas = [k * math.pi / 4 for k in range(8)]

    # Body: revolve the closed profile (apex, ogive arc, cylinder, base, axis).
    p_apex = occ.addPoint(0, 0, 0)
    p_sh = occ.addPoint(body.l_ogive, R, 0)
    p_rim = occ.addPoint(total, R, 0)
    p_base = occ.addPoint(total, 0, 0)
    p_c = occ.addPoint(body.xc, body.yc, 0)
    loop = occ.addCurveLoop([occ.addCircleArc(p_apex, p_c, p_sh), occ.addLine(p_sh, p_rim),
                             occ.addLine(p_rim, p_base), occ.addLine(p_base, p_apex)])
    profile = occ.addPlaneSurface([loop])
    solid = [(3, t) for d, t in occ.revolve([(2, profile)], 0, 0, 0, 1, 0, 0, 2 * math.pi) if d == 3]
    occ.remove([(0, p_c)])

    # Zones. The inner zones (nose, body, aft) and the far zone are built and
    # split separately, so the seams between them are two independent copies
    # of the same surface: each side keeps its own block layout and the
    # meshes are stitched non-conformally in OpenFOAM.
    zones = {}
    for name, (x0, x1) in (("nose", xs[0:2]), ("body", xs[1:3]), ("aft", xs[2:4])):
        cyl = occ.addCylinder(x0, 0, 0, x1 - x0, 0, 0, r_z)
        (_, v), = occ.cut([(3, cyl)], solid, removeTool=False)[0]
        zones[v] = name
    domain = occ.addCylinder(x_in, 0, 0, x_out - x_in, 0, 0, r_ff)
    hole = occ.addCylinder(xs[0], 0, 0, xs[3] - xs[0], 0, 0, r_z)
    (_, far), = occ.cut([(3, domain)], [(3, hole)])[0]
    occ.remove(solid, recursive=True)
    occ.synchronize()

    # Nose core: a square frustum whose sides leave the tip patch at 45 - beta
    # from the axis, halfway between the wall normal and the wall itself, so
    # neither the core blocks nor the ring blocks around them get a sharp
    # wedge along the patch edge. Its half-width is `a` where it meets the
    # ogive at phi = 0 and grows upstream to the zone start.
    beta = math.radians(body.slope_deg(0.0))
    k = math.tan(math.pi / 4 - beta)
    x_a = bisect(lambda x: body.r(x) - a, 0.0, body.l_ogive)
    h = lambda x: a - k * (x - x_a)
    x_45 = bisect(lambda x: body.r(x) - SQ2 * h(x), x_a, x_a + a / k)   # corner meets the wall
    x_end = x_45 + 0.5 * (x_a + a / k - x_45)                           # still inside the body

    # Nose interface, built before the sector cuts. Cutting the nose zone
    # with the core frustum alone leaves the tip patch bounded by four whole
    # edges, one per core side. A closed ruled funnel from that loop out to
    # the circle r = r_z at x = xTipInterface then separates the ring around
    # the core (upstream) from the ring along the ogive. The sector cuts
    # applied afterwards cross the funnel transversally; nothing in the
    # boolean is tangent or coplanar.
    nose = fragment({v: n for v, n in zones.items() if n == "nose"},
                    [square_prism(xs[0], h(xs[0]), x_end, h(x_end))])
    zones = {v: n for v, n in zones.items() if n != "nose"} | nose
    funnel = tip_funnel(a, x_a, x_45, P["xTipInterface"], r_z)

    # Inner zones: 8 sectors; the 45-degree cuts stop at the core corners.
    tools = [half_plane(t, xs[0], xs[3], 0.0, r_z) for t in thetas[0::2]]
    for t in thetas[1::2]:
        # Upstream of the base the cut follows the core's corner line to
        # x_end, then runs inside the body and meets the base plane only at
        # the rim, so the base disk is not cut; behind the base it stops at
        # the wake core's corner.
        tools += [polygon_half_plane(t, [(xs[0], SQ2 * h(xs[0])), (x_end, SQ2 * h(x_end)),
                                         (x_end, 0.0), (total - 0.5 * R, 0.0), (total, R),
                                         (total, r_z), (xs[0], r_z)]),
                  half_plane(t, total, xs[3], SQ2 * s_w, r_z)]
    tools += funnel + [square_prism(total, s_w, xs[3], s_w),
                       (3, occ.addCylinder(total, 0, 0, xs[3] - total, 0, 0, R)),
                       annulus(body.l_ogive, 0.0, r_z),          # shoulder
                       annulus(total, R, r_z)]                    # base plane, outside the base disk
    zones = fragment(zones, tools)

    # Far zone: core + two rings upstream and downstream, one ring around
    # the inner zones.
    far_tools = [half_plane(t, x_in, x_out, 0.0, r_ff) for t in thetas[0::2]]
    for t in thetas[1::2]:
        far_tools += [half_plane(t, x_in, xs[0], SQ2 * s_f, r_ff),
                      half_plane(t, xs[0], xs[3], r_z, r_ff),
                      half_plane(t, xs[3], x_out, SQ2 * s_f, r_ff)]
    far_tools += [square_prism(x_in, s_f, xs[0], s_f), square_prism(xs[3], s_f, x_out, s_f),
                  (3, occ.addCylinder(x_in, 0, 0, xs[0] - x_in, 0, 0, r_z)),
                  (3, occ.addCylinder(xs[3], 0, 0, x_out - xs[3], 0, 0, r_z)),
                  annulus(xs[0], r_z, r_ff), annulus(xs[3], r_z, r_ff)]
    zones = zones | fragment({far: "far"}, far_tools)

    # Faces: wall, outer patches, seams. A face of one zone that lies on the
    # inner/far interface is a seam; inner|inner seams are shared faces.
    owners: dict[int, list[int]] = {}
    for v in zones:
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            owners.setdefault(s, []).append(v)
    groups: dict[str, list[int]] = {}
    for s, own in owners.items():
        names = sorted({zones[v] for v in own})
        if len(own) == 2 and len(names) == 1:
            continue                                   # internal block face
        bb = gmsh.model.getBoundingBox(2, s)
        rmax = max(abs(v) for v in bb[1:3] + bb[4:6])
        on_end = abs(bb[0] - bb[3]) < PAD and min(abs(bb[0] - xs[0]), abs(bb[0] - xs[3])) < PAD
        on_interface = (rmax < r_z + PAD and bb[0] > xs[0] - PAD and bb[3] < xs[3] + PAD
                        and (on_end or _on_cylinder(s, r_z)))
        if len(names) == 2:
            key = "seam_" + "_".join(names)
        elif on_interface:
            key = f"seam_{names[0]}_far" if names[0] != "far" else "seam_far_inner"
        elif bb[3] < x_in + PAD:
            key = "inlet"
        elif bb[0] > x_out - PAD:
            key = "outlet"
        elif rmax > r_ff - PAD and _on_cylinder(s, r_ff):
            key = "farfield"
        else:
            key = "fuselage"
        groups.setdefault(key, []).append(s)

    for name in ("nose", "body", "aft", "far"):
        gmsh.model.addPhysicalGroup(3, [v for v, n in zones.items() if n == name], name=name)
    for key, surfs in groups.items():
        gmsh.model.addPhysicalGroup(2, surfs, name=key)
    return {"zones": zones, "groups": groups, "xs": xs,
            "tip": {"beta": math.degrees(beta), "x_a": x_a, "x_45": x_45, "x_end": x_end,
                    "h_in": h(xs[0])}}


def tip_funnel(a: float, x_a: float, x_45: float, x_f: float, r_z: float) -> list:
    """Closed ruled surface from the tip-patch loop to the circle r = r_z at x = x_f."""
    occ = gmsh.model.occ
    tol = 0.01 * a   # OCC's spline edges and padded boxes miss the analytic x_45 slightly
    edges = []
    for _, c in gmsh.model.getEntities(1):
        bb = gmsh.model.getBoundingBox(1, c)
        if (bb[0] > x_a - tol and bb[3] < x_45 + tol
                and max(abs(v) for v in bb[1:3] + bb[4:6]) < SQ2 * a + tol):
            lo, hi = gmsh.model.getParametrizationBounds(1, c)
            mid = gmsh.model.getValue(1, c, [0.5 * (lo[0] + hi[0])])
            edges.append((math.atan2(mid[2], mid[1]) % (2 * math.pi), c))
    # One per core side, plus a split where the ogive's parameter seam crosses.
    if not 4 <= len(edges) <= 5:
        raise RuntimeError(f"expected 4-5 tip-patch edges, found {len(edges)}")
    # Chain the edges head to tail (flipping where needed) into a closed loop.
    def ends(c):
        lo, hi = gmsh.model.getParametrizationBounds(1, c)
        return [gmsh.model.getValue(1, c, [u]) for u in (lo[0], hi[0])]

    def same(p, q):
        return sum((p[i] - q[i]) ** 2 for i in range(3)) < (1e-3 * a) ** 2

    todo = [c for _, c in sorted(edges)]
    chain = [todo.pop(0)]
    tip = ends(chain[0])[1]
    while todo:
        c = next(c for c in todo if any(same(tip, e) for e in ends(c)))
        todo.remove(c)
        head, tail = ends(c)
        chain.append(c if same(tip, head) else -c)
        tip = tail if same(tip, head) else head
    loop = occ.addWire(chain)

    # Outer wire: one arc per loop edge, between the same angles, so the
    # ruled surface pairs edge with edge and its rulings run outward from
    # each patch vertex.
    centre = occ.addPoint(x_f, 0, 0)

    def on_circle(p):
        r = math.hypot(p[1], p[2])
        return occ.addPoint(x_f, r_z * p[1] / r, r_z * p[2] / r)

    first = ends(abs(chain[0]))[0 if chain[0] > 0 else 1]
    p_prev = p0 = on_circle(first)
    arcs = []
    for i, c in enumerate(chain):
        head, tail = ends(abs(c))
        p_next = p0 if i == len(chain) - 1 else on_circle(tail if c > 0 else head)
        arcs.append(occ.addCircleArc(p_prev, centre, p_next))
        p_prev = p_next
    surf = [dt for dt in occ.addThruSections([loop, occ.addWire(arcs)], makeSolid=False,
                                             makeRuled=True) if dt[0] == 2]
    occ.synchronize()
    twist = 0.0
    for _, sf in surf:
        (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, sf)
        for i in range(9):
            u = u0 + (u1 - u0) * i / 8
            p, q = (gmsh.model.getValue(2, sf, [u, v]) for v in (v0, v1))
            d = math.atan2(p[2], p[1]) - math.atan2(q[2], q[1])
            twist = max(twist, abs(math.atan2(math.sin(d), math.cos(d))))
    best = (twist, surf)
    if best[0] > math.radians(10):
        raise RuntimeError(f"tip funnel rulings twisted by {math.degrees(best[0]):.1f} deg")
    return best[1]


def _on_cylinder(s: int, r: float) -> bool:
    """True if face s lies on the cylinder of radius r about the x axis."""
    (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, s)
    pts = gmsh.model.getValue(2, s, [u0, v0, u1, v1, 0.5 * (u0 + u1), 0.5 * (v0 + v1)])
    return all(abs(math.hypot(pts[i + 1], pts[i + 2]) - r) < PAD for i in range(0, 9, 3))


def bisect(f, lo: float, hi: float) -> float:
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def real_edges(face):
    """Edges of a face, without OCC's zero-length pole edges (surface-of-revolution apex)."""
    return [e for e in gmsh.model.getBoundary([face], oriented=False)
            if gmsh.model.occ.getMass(1, e[1]) > 1e-9]


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
        got = sum(occ.getMass(3, v) for v, n in info["zones"].items() if n == name)
        row(f"zone {name} volume (mm^3)", got, math.pi * r_z**2 * (x1 - x0) - vb)

    # Block topology: every block must be a hexahedron (6 faces of 4 edges,
    # 12 edges, 8 corners) so it can carry a structured mesh.
    expected = {"nose": 20, "body": 16, "aft": 28, "far": 48}
    for name, want in expected.items():
        blocks = [v for v, n in info["zones"].items() if n == name]
        bad = []
        for v in blocks:
            faces = gmsh.model.getBoundary([(3, v)], oriented=False)
            edges = {c for f in faces for _, c in real_edges(f)}
            corners = {p for c in edges for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)}
            per_face = [len(real_edges(f)) for f in faces]
            if len(faces) != 6 or len(edges) != 12 or len(corners) != 8 or set(per_face) != {4}:
                bad.append(f"{v}: {len(faces)}f/{len(edges)}e/{len(corners)}c {per_face}")
        good = len(blocks) == want and not bad
        ok &= good
        lines.append(f"  zone {name:4s}: {len(blocks):3d} blocks (expected {want}), "
                     f"{len(blocks) - len(bad)} hexahedral  {'ok' if good else 'FAIL'}")
        for b in bad[:6]:
            lines.append(f"      not hex: block {b}")

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
        t = info["tip"]
        print(f"nose core: tip half-angle {t['beta']:.2f} deg, patch edge x {t['x_a']:.3f}-{t['x_45']:.3f} mm, "
              f"half-width {t['h_in']:.2f} mm at the zone start")
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
