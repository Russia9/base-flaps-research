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
import struct
import sys
from pathlib import Path

import gmsh
from create_case import parse_scalar

PARAMS = {
    "xInlet": -800.0, "xOutlet": 5600.0, "rFar": 600.0,
    "xZoneStart": -40.0, "xNoseEnd": 40.0, "xAftStart": 620.0, "xZoneEnd": 1200.0,
    "rZone": 120.0, "coreTip": 0.5, "xTipInterface": -30.0, "coreWake": 20.0, "coreFar": 60.0,
    # spacing (mm)
    "nTheta": 16, "firstLayer": 0.003, "layerRatio": 1.2, "nLayers": 20, "growth": 1.1,
    "hUpstream": 2.0, "hOgive": 1.0, "hShoulder": 0.5, "hWall": 2.0, "hRingOut": 4.0,
    "hWake": 4.0, "hSeam": 4.0, "hFarMid": 8.0, "hFar": 30.0,
    "xWakeSplit": 820.0, "hRelax": 2.0, "hWakeRadial": 1.0, "firstLayerTip": 0.02,
}
INTEGER = ("nTheta", "nLayers")


def read_params(case: Path) -> dict[str, float]:
    text = (case / "system" / "gmshParams").read_text()
    params = {}
    for name, default in PARAMS.items():
        try:
            params[name] = parse_scalar(text, name)
        except ValueError:
            params[name] = default
    for name in INTEGER:
        params[name] = int(params[name])
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


def x_cylinder(x0: float, dx: float, r: float) -> tuple[int, int]:
    """Cylinder along x whose curved surface is eight 45-degree patches.

    A single occ.addCylinder has one periodic face with its parameter seam on
    a sector boundary, and structured faces that end on the seam fold (Gmsh
    picks the wrong parameter, 0 or 2 pi, at the seam corner). Here each
    sector's patch is a separate revolved, non-periodic face; the patches and
    the two end disks are sewn into one solid.
    """
    occ = gmsh.model.occ
    faces = []
    for k in range(8):
        line = occ.addLine(occ.addPoint(x0, r, 0), occ.addPoint(x0 + dx, r, 0))
        occ.rotate([(1, line)], 0, 0, 0, 1, 0, 0, k * math.pi / 4)
        faces += [t for d, t in occ.revolve([(1, line)], 0, 0, 0, 1, 0, 0, math.pi / 4) if d == 2]
    for x in (x0, x0 + dx):
        c = occ.addPoint(x, 0, 0)
        p = [occ.addPoint(x, r * math.cos(k * math.pi / 4), r * math.sin(k * math.pi / 4)) for k in range(8)]
        arcs = [occ.addCircleArc(p[k], c, p[(k + 1) % 8]) for k in range(8)]
        faces.append(occ.addPlaneSurface([occ.addCurveLoop(arcs)]))
        occ.remove([(0, c)])
    return (3, occ.addVolume([occ.addSurfaceLoop(faces, sewing=True)]))


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
        pieces, _ = occ.cut([x_cylinder(x0, x1 - x0, r_z)], solid, removeTool=False)
        zones |= {v: name for d, v in pieces if d == 3}
    pieces, _ = occ.cut([x_cylinder(x_in, x_out - x_in, r_ff)], [x_cylinder(xs[0], xs[3] - xs[0], r_z)])
    far = {v: "far" for d, v in pieces if d == 3}
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
                       x_cylinder(total, xs[3] - total, R),
                       annulus(body.l_ogive, 0.0, r_z),          # shoulder
                       annulus(total, R, r_z),                    # base plane, outside the base disk
                       annulus(P["xWakeSplit"], 0.0, r_z)]        # wake: near part | relaxing part
    zones = fragment(zones, tools)

    # Far zone: core + two rings upstream and downstream, one ring around
    # the inner zones.
    far_tools = [half_plane(t, x_in, x_out, 0.0, r_ff) for t in thetas[0::2]]
    for t in thetas[1::2]:
        far_tools += [half_plane(t, x_in, xs[0], SQ2 * s_f, r_ff),
                      half_plane(t, xs[0], xs[3], r_z, r_ff),
                      half_plane(t, xs[3], x_out, SQ2 * s_f, r_ff)]
    far_tools += [square_prism(x_in, s_f, xs[0], s_f), square_prism(xs[3], s_f, x_out, s_f),
                  x_cylinder(x_in, xs[0] - x_in, r_z),
                  x_cylinder(xs[3], x_out - xs[3], r_z),
                  annulus(xs[0], r_z, r_ff), annulus(xs[3], r_z, r_ff)]
    zones = zones | fragment(far, far_tools)

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
            # Named by side and position: stitchMesh couples each smooth piece
            # (end disk or cylinder) separately; across the 90-degree edges
            # between them its projection fails.
            where = ("up" if abs(bb[0] - xs[0]) < PAD else "down") if on_end else "side"
            key = f"seam_{'far' if names[0] == 'far' else 'inner'}_{where}"
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
    # One per core side, or one per 45-degree sector once the zone is split.
    if not 4 <= len(edges) <= 8:
        raise RuntimeError(f"expected 4-8 tip-patch edges, found {len(edges)}")
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
    p0 = on_circle(first)
    ring = [p0]
    for c in chain[:-1]:
        head, tail = ends(abs(c))
        ring.append(on_circle(tail if c > 0 else head))
    # OCC orients and starts the loop wire itself; try every start and both
    # directions for the arc wire and keep the untwisted surface.
    best = None
    n = len(ring)
    rings = [ring[k:] + ring[:k] for k in range(n)]
    rings += [r[::-1] for r in rings]
    for r in rings:
        pts = r + [r[0]]
        arcs = [occ.addCircleArc(pts[i], centre, pts[i + 1]) for i in range(len(pts) - 1)]
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
        if best is None or twist < best[0]:
            if best:
                occ.remove(best[1], recursive=True)
            best = (twist, surf)
        else:
            occ.remove(surf, recursive=True)
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


# ── cell spacing along an edge ───────────────────────────────────────────────

class End:
    """Cell sizes growing away from one end of an edge.

    Cell k is h * qLayer^min(k, nLayers) * q^max(0, k - nLayers).
    """

    def __init__(self, h: float, q: float = 1.0, n_layers: int = 0, q_layer: float = 1.0):
        self.h, self.q, self.n_layers, self.q_layer = h, q, n_layers, q_layer

    def size(self, k: int) -> float:
        return self.h * self.q_layer ** min(k, self.n_layers) * self.q ** max(0, k - self.n_layers)


def size_driven(length: float, a: End, b: End, hmax: float) -> list[float]:
    """Grow cells from both ends (capped at hmax) until the edge is full."""
    A: list[float] = []
    B: list[float] = []
    total = 0.0
    while True:
        na = min(a.size(len(A)), hmax)
        nb = min(b.size(len(B)), hmax)
        side, nxt = (A, na) if na <= nb else (B, nb)
        if total + nxt > length:
            break
        side.append(nxt)
        total += nxt
    gap = length - total
    if gap > 0.5 * min(na, nb):
        A.append(gap)
        gap = 0.0
    sizes = A + B[::-1]
    lo = min(a.n_layers, len(A))
    hi = len(sizes) - min(b.n_layers, len(B))
    free = sum(sizes[lo:hi])
    if free <= 0:
        return [x * length / sum(sizes) for x in sizes]
    scale = 1 + gap / free
    return sizes[:lo] + [x * scale for x in sizes[lo:hi]] + sizes[hi:]


def fixed_n(length: float, n: int, start: End | None = None, end: End | None = None) -> list[float]:
    """n cells over length, clustered at one end (or uniform)."""
    if start is None and end is None:
        return [length / n] * n
    if start is None:
        return fixed_n(length, n, end)[::-1]
    nl = min(start.n_layers, n - 1)
    layers = [start.h * start.q_layer**k for k in range(nl)]
    rem = length - sum(layers)
    m = n - nl
    base = start.h * start.q_layer**nl
    if rem <= 0 or base * m * 1e-3 > rem:
        return [length / n] * n

    def total(r: float) -> float:
        return base * (m if abs(r - 1) < 1e-12 else (r**m - 1) / (r - 1))

    lo, hi = 1e-3, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < rem:
            lo = mid
        else:
            hi = mid
    r = 0.5 * (lo + hi)
    return layers + [base * r**i for i in range(m)]


# ── step 2: structured surface mesh ──────────────────────────────────────────

class Curve:
    """Endpoints, length and an arc-length table of one model curve."""

    def __init__(self, c: int):
        self.tag = c
        (lo,), (hi,) = gmsh.model.getParametrizationBounds(1, c)
        self.u = [lo + (hi - lo) * i / 400 for i in range(401)]
        pts = gmsh.model.getValue(1, c, self.u)
        self.pts = [pts[3 * i:3 * i + 3] for i in range(401)]
        self.s = [0.0]
        for i in range(400):
            self.s.append(self.s[-1] + math.dist(self.pts[i], self.pts[i + 1]))
        self.length = self.s[-1]
        self.p, self.q = self.pts[0], self.pts[-1]
        self.rp, self.rq = math.hypot(*self.p[1:]), math.hypot(*self.q[1:])

    def u_at(self, s: float) -> float:
        i = min(max(1, next((k for k, v in enumerate(self.s) if v >= s), 400)), 400)
        f = (s - self.s[i - 1]) / max(self.s[i] - self.s[i - 1], 1e-300)
        return self.u[i - 1] + f * (self.u[i] - self.u[i - 1])


def classify(cv: Curve, inner: bool, body: Body, P: dict, tip: dict) -> str:
    """Edge family: circ, core (axis or apex to a core-square side), radial or axial,
    the last two with the zone segment they belong to."""
    tol = 1e-3
    r_z, r_ff, R = P["rZone"], P["rFar"], body.R
    xs = [P["xZoneStart"], P["xNoseEnd"], P["xAftStart"], P["xZoneEnd"]]
    if (cv.rp < tol) != (cv.rq < tol):
        return "core"
    if cv.rp > tol and cv.rq > tol:
        d = math.atan2(cv.p[2], cv.p[1]) - math.atan2(cv.q[2], cv.q[1])
        if abs(math.atan2(math.sin(d), math.cos(d))) > math.radians(30):
            return "circ"
    x0, x1 = sorted((cv.p[0], cv.q[0]))
    r0, r1 = sorted((cv.rp, cv.rq))
    if abs(r1 - r0) > abs(x1 - x0):           # radial (the tip funnel's lines included)
        if inner:
            if abs(r1 - r_z) < tol:
                return "radial:ring"
            if abs(r1 - R) < tol and x0 > body.total - tol:
                return "radial:wake"
        else:
            if abs(r1 - r_z) < tol:
                return "radial:farcore"
            if abs(r1 - r_ff) < tol:
                return "radial:far"
    else:
        def at(a, b):
            return abs(a - b) < tol
        if inner:
            if at(x0, xs[0]) and x1 < tip["x_45"] + 0.01:
                return "axial:za"
            if at(x1, xs[1]) and x0 < tip["x_45"] + 0.01:
                return "axial:zb"
            for name, (a, b) in (("b1", (xs[1], body.l_ogive)), ("b2", (body.l_ogive, xs[2])),
                                 ("a1", (xs[2], body.total)), ("wk", (body.total, P["xWakeSplit"])),
                                 ("wk2", (P["xWakeSplit"], xs[3]))):
                if at(x0, a) and at(x1, b):
                    return "axial:" + name
        else:
            for name, (a, b) in (("fu", (P["xInlet"], xs[0])), ("fm", (xs[0], xs[3])),
                                 ("fd", (xs[3], P["xOutlet"]))):
                if at(x0, a) and at(x1, b):
                    return "axial:" + name
    raise RuntimeError(f"unclassified edge {cv.tag}: {cv.p} -> {cv.q}")


def surface_mesh(body: Body, P: dict, info: dict) -> dict:
    """Assign structured spacing to every block edge, mesh curves and faces (quads)."""
    zones, R, r_z = info["zones"], body.R, P["rZone"]
    q, nt = P["growth"], P["nTheta"]
    wall = End(P["firstLayer"], q, P["nLayers"], P["layerRatio"])
    tol = 1e-3

    inner_vols = {v for v, n in zones.items() if n != "far"}
    curves, family, inner_of = {}, {}, {}
    for _, c in gmsh.model.getEntities(1):
        cv = Curve(c)
        if cv.length < 1e-9:
            continue                                   # OCC pole edge at the apex
        faces = gmsh.model.getAdjacencies(1, c)[0]
        vols = {v for f in faces for v in gmsh.model.getAdjacencies(2, f)[0]}
        inner_of[c] = bool(vols & inner_vols)
        curves[c] = cv
        family[c] = classify(cv, inner_of[c], body, P, info["tip"])

    def first(name, pred=lambda cv: True):
        return next(cv for c, cv in curves.items() if family[c] == name and pred(cv))

    # Controlling edges set each family's cell count; every other edge of the
    # family gets the same count with its own distribution.
    ctrl = {}
    ctrl["radial:ring"] = size_driven(r_z - R, wall, End(P["hRingOut"], q), P["hRingOut"])
    ctrl["radial:wake"] = size_driven(R - P["coreWake"], End(P["coreWake"] / nt, q), wall, R)
    ctrl["radial:farcore"] = fixed_n(r_z - P["coreFar"], max(1, round((r_z - P["coreFar"])
                                                                     / (P["coreFar"] / nt))))
    ctrl["radial:far"] = size_driven(P["rFar"] - r_z, End(P["hSeam"], q), End(P["hFar"], q), P["hFar"])
    # The core blocks meet the tip patch at a shallow angle; 3 um layers
    # interpolated across them cross the wall there. The patch (the first
    # ~1.5 mm of a sharp tip, where the boundary layer starts from nothing)
    # gets a thicker first cell instead.
    tip_wall = End(P["firstLayerTip"], q)
    ctrl["axial:za"] = size_driven(-P["xZoneStart"], End(P["hUpstream"], q), tip_wall, P["hUpstream"])
    apex_mer = first("core", lambda cv: min(cv.p[0], cv.q[0]) > -tol and max(cv.rp, cv.rq) < 1.0)
    h_tip = apex_mer.length / nt
    zb = first("axial:zb", lambda cv: min(cv.rp, cv.rq) < r_z - tol)
    ctrl["axial:zb"] = size_driven(zb.length, End(h_tip, q), End(P["hOgive"], q), P["hOgive"])
    b1 = first("axial:b1", lambda cv: min(cv.rp, cv.rq) < r_z - tol)
    ctrl["axial:b1"] = size_driven(b1.length, End(P["hOgive"], q), End(P["hShoulder"], q), P["hOgive"])
    ctrl["axial:b2"] = size_driven(P["xAftStart"] - body.l_ogive, End(P["hShoulder"], q),
                                   End(P["hWall"], q), P["hWall"])
    ctrl["axial:a1"] = size_driven(body.total - P["xAftStart"], End(P["hWall"], q), wall, P["hWall"])
    # Wake: the near part (base -> xWakeSplit) keeps the rim's radial wall
    # clustering at both ends and only grows axially to hRelax; past it all
    # axial edges share one spacing and the radial spacing relaxes towards
    # the zone end. One block varying in both directions gets sheared, and
    # 3 um radial cells against 4 mm axial ones break the aspect ratio.
    ctrl["axial:wk"] = size_driven(P["xWakeSplit"] - body.total, wall, End(P["hRelax"], q), P["hRelax"])
    ctrl["axial:wk2"] = size_driven(P["xZoneEnd"] - P["xWakeSplit"], End(P["hRelax"], q),
                                    End(P["hWake"], q), P["hWake"])
    ctrl["axial:fu"] = size_driven(P["xZoneStart"] - P["xInlet"], End(P["hFar"], q),
                                   End(P["hSeam"], q), P["hFar"])
    ctrl["axial:fm"] = fixed_n(1.0, max(1, round((P["xZoneEnd"] - P["xZoneStart"]) / P["hFarMid"])))
    ctrl["axial:fd"] = size_driven(P["xOutlet"] - P["xZoneEnd"], End(P["hSeam"], q), End(P["hFar"], q),
                                   P["hFar"])
    counts = {k: len(v) for k, v in ctrl.items()} | {"circ": nt, "core": nt}

    def sizes_for(c: int) -> list[float]:
        """Cell sizes from the curve's parametric start to its end."""
        cv, fam = curves[c], family[c]
        n, L = counts[fam], cv.length
        at_r_lo = cv.rp <= cv.rq                          # parametric start is the inner end
        at_x_lo = cv.p[0] <= cv.q[0]                      # parametric start is upstream
        outer = min(cv.rp, cv.rq) > r_z - tol             # an axial line on r = r_z
        at_wake_end = fam.startswith("radial:") and abs(cv.p[0] - P["xZoneEnd"]) < tol
        if fam in ("circ", "core") or (fam.startswith("axial:") and outer and inner_of[c]
                                       and fam != "axial:wk2"):
            sz = fixed_n(L, n)
        elif fam == "radial:ring":
            sz = fixed_n(L, n, End(P["hWakeRadial"], q) if at_wake_end else wall)
            sz = sz if at_r_lo else sz[::-1]
        elif fam == "radial:wake":
            sz = fixed_n(L, n, None, End(P["hWakeRadial"], q) if at_wake_end else wall)
            sz = sz if at_r_lo else sz[::-1]
        elif fam == "radial:farcore":
            sz = fixed_n(L, n)
        elif fam == "radial:far":
            sz = fixed_n(L, n, End(P["hSeam"]))           # single ratio, from the seam
            sz = sz if at_r_lo else sz[::-1]
        elif fam == "axial:fu":
            sz = fixed_n(L, n, None, End(P["hSeam"]))     # single ratio, towards the seam
            sz = sz if at_x_lo else sz[::-1]
        elif fam == "axial:fd":
            sz = fixed_n(L, n, End(P["hSeam"]))           # single ratio, from the seam
            sz = sz if at_x_lo else sz[::-1]
        elif fam == "axial:za":                           # inner: wall (tip patch) at the end
            sz = fixed_n(L, n, None, tip_wall)
            sz = sz if at_x_lo else sz[::-1]
        elif fam == "axial:zb":                           # inner: starts at the tip patch
            sz = fixed_n(L, n, End(h_tip, q))
            sz = sz if at_x_lo else sz[::-1]
        elif fam == "axial:fm":
            sz = fixed_n(L, n)
        else:                                             # same length on every edge of the family
            ref = ctrl[fam]
            sz = [x * L / sum(ref) for x in ref]
            sz = sz if at_x_lo else sz[::-1]
        return sz

    sizes = {c: sizes_for(c) for c in curves}

    # Transfinite curves, faces (quads) and, for step 3, volumes. Far-zone
    # edges carry one geometric ratio each and use Gmsh's own Progression:
    # moving nodes on the graded far edges next to a cylinder's parameter
    # seam folds the faces there. Inner-zone edges need multi-stage spacing
    # (the wall stack), so their nodes are moved after meshing the curves.
    relocate = {}
    for c, sz in sizes.items():
        if inner_of[c]:
            gmsh.model.mesh.setTransfiniteCurve(c, len(sz) + 1)
            relocate[c] = sz
        else:
            ratio = sz[1] / sz[0] if len(sz) > 1 else 1.0
            gmsh.model.mesh.setTransfiniteCurve(c, len(sz) + 1, "Progression", ratio)
    bad_faces = []
    for _, f in gmsh.model.getEntities(2):
        signed = [t for _, t in gmsh.model.getBoundary([(2, f)], oriented=True)]
        loop = [abs(t) for t in signed if abs(t) in sizes]
        # Faces touching the apex also carry OCC's zero-length pole edge;
        # give those their four corners in loop order. Gmsh finds the rest.
        corners = []
        if len(loop) != len(signed):
            ends = {c: [abs(p) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)] for c in loop}
            todo = list(loop)
            here = ends[todo[0]][0]
            while todo:                        # walk the loop edge to edge
                c = next(c for c in todo if here in ends[c])
                todo.remove(c)
                corners.append(here)
                here = ends[c][1] if ends[c][0] == here else ends[c][0]
        ends = {c: {p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)} for c in loop}
        opposite = [(a, b) for i, a in enumerate(loop) for b in loop[i + 1:] if not ends[a] & ends[b]]
        if len(loop) != 4 or len(opposite) != 2 or any(len(sizes[a]) != len(sizes[b]) for a, b in opposite):
            bad_faces.append((f, loop, [len(sizes[c]) for c in loop]))
        gmsh.model.mesh.setTransfiniteSurface(f, cornerTags=corners)
        gmsh.model.mesh.setRecombine(2, f)

    # Mesh the curves, move their nodes to the exact distributions, then mesh
    # the faces from those nodes (existing curve meshes are kept).
    gmsh.model.mesh.generate(1)
    for c, sz in relocate.items():
        cv = curves[c]
        tags, _, u = gmsh.model.mesh.getNodes(1, c, includeBoundary=False, returnParametricCoord=True)
        order = sorted(range(len(tags)), key=lambda i: u[i])
        if cv.u[-1] < cv.u[0]:
            order = order[::-1]
        s, total = 0.0, sum(sz)
        for k, i in enumerate(order):
            s += sz[k]
            uu = cv.u_at(s * cv.length / total)
            gmsh.model.mesh.setNode(tags[i], gmsh.model.getValue(1, c, [uu]), [uu])
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.model.mesh.generate(2)

    # Cells per block: product of the counts on the three edges at one corner.
    cells = {}
    for v, name in zones.items():
        edges = {abs(t) for f in gmsh.model.getBoundary([(3, v)], oriented=False)
                 for _, t in gmsh.model.getBoundary([f], oriented=False) if abs(t) in sizes}
        ends = {c: {p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)} for c in edges}
        corner = min(next(iter(ends.values())))
        at_corner = [c for c in edges if corner in ends[c]]
        cells[name] = cells.get(name, 0) + math.prod(len(sizes[c]) for c in at_corner[:3])
    return {"counts": counts, "sizes": sizes, "curves": curves, "family": family,
            "bad_faces": bad_faces, "h_tip": h_tip, "cells": cells, "wall": wall}


def surface_checks(body: Body, P: dict, info: dict, sm: dict) -> list[str]:
    lines, ok = [], True
    lines.append("  edge families (cells):")
    for k in sorted(sm["counts"]):
        lines.append(f"    {k:16s} {sm['counts'][k]}")
    good = not sm["bad_faces"]
    ok &= good
    lines.append(f"  faces with opposite-edge count mismatch: {len(sm['bad_faces'])}  {'ok' if good else 'FAIL'}")
    for f, loop, n in sm["bad_faces"][:5]:
        lines.append(f"      face {f}: edges {loop} cells {n}")

    # Element types on all faces: transfinite + recombine must give quads only.
    tri = quad = 0
    for _, f in gmsh.model.getEntities(2):
        for et, tags, _ in zip(*gmsh.model.mesh.getElements(2, f)):
            name = gmsh.model.mesh.getElementProperties(et)[0]
            tri += len(tags) if name.startswith("Triangle") else 0
            quad += len(tags) if name.startswith("Quadrilateral") else 0
    ok &= tri == 0
    lines.append(f"  surface elements: {quad} quads, {tri} triangles  {'ok' if tri == 0 else 'FAIL'}")

    # Quad shape: minimum scaled Jacobian (1 = rectangle, <= 0 = folded).
    worst = (2.0, None)
    for _, f in gmsh.model.getEntities(2):
        for et, tags, _ in zip(*gmsh.model.mesh.getElements(2, f)):
            qual = gmsh.model.mesh.getElementQualities(tags, "minSJ")
            i = min(range(len(qual)), key=qual.__getitem__)
            if qual[i] < worst[0]:
                worst = (qual[i], f)
    good = worst[0] > 0.1
    ok &= good
    lines.append(f"  worst quad scaled Jacobian {worst[0]:.3f} (face {worst[1]})  {'ok' if good else 'FAIL'}")

    # First cell off the wall, measured on the meshed wall-normal edges.
    first = P["firstLayer"]
    worst = 0.0
    for c, cv in sm["curves"].items():
        fam = sm["family"][c]
        if fam == "radial:ring" and cv.p[0] < body.total + 1e-3 \
                and abs(min(cv.rp, cv.rq) - body.r(cv.p[0])) < 1e-3 \
                or fam == "axial:wk" and min(cv.rp, cv.rq) < body.R + 1e-3:
            sz = sm["sizes"][c]
            h = sz[0] if (cv.rp <= cv.rq if fam.startswith("radial") else cv.p[0] <= cv.q[0]) else sz[-1]
            worst = max(worst, abs(h - first) / first)
    ok &= worst < 0.05
    lines.append(f"  first wall cell {first * 1000:.1f} um, max deviation {worst * 100:.2f} %  "
                 f"{'ok' if worst < 0.05 else 'FAIL'}")

    # Wall-cell aspect ratio: longest wall-parallel edge over the first cell.
    longest = 0.0
    for s in info["groups"]["fuselage"]:
        for et, tags, nodes in zip(*gmsh.model.mesh.getElements(2, s)):
            k = gmsh.model.mesh.getElementProperties(et)[3]
            for e in range(len(tags)):
                ids = nodes[k * e:k * e + k]
                xyz = [gmsh.model.mesh.getNode(n)[0] for n in ids]
                longest = max(longest, max(math.dist(xyz[i], xyz[(i + 1) % k]) for i in range(k)))
    ar = longest / first
    ok &= ar < 1000
    lines.append(f"  wall cell aspect ratio: longest wall edge {longest:.3f} mm -> {ar:.0f}  "
                 f"{'ok' if ar < 1000 else 'FAIL'}")
    lines.append(f"  tip patch cell {sm['h_tip'] * 1000:.0f} um; estimated cells per zone: "
                 + ", ".join(f"{k} {v / 1e6:.2f} M" for k, v in sm["cells"].items())
                 + f"; total {sum(sm['cells'].values()) / 1e6:.2f} M")
    lines.append("  " + ("ALL SURFACE CHECKS PASSED" if ok else "SOME SURFACE CHECKS FAILED"))
    return lines


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
    expected = {"nose": 20, "body": 16, "aft": 48, "far": 48}
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


def hex_corners(v: int) -> list[int]:
    """The 8 corners of a hexahedral block in Gmsh's transfinite order (one face
    in loop order, then the opposite face, each corner above its partner), or
    [] to let Gmsh choose. Needed for blocks touching the apex, whose zero-length
    pole edge otherwise makes Gmsh collapse cells there."""
    faces = gmsh.model.getBoundary([(3, v)], oriented=False)
    if all(len(real_edges(f)) == len(gmsh.model.getBoundary([f], oriented=False)) for f in faces):
        return []
    ends = {}
    for f in faces:
        for _, c in real_edges(f):
            ends[c] = [abs(p) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)]
    face_pts = []
    for f in faces:
        loop = [c for _, c in real_edges(f)]
        here, ring, todo = ends[loop[0]][0], [], list(loop)
        while todo:
            c = next(c for c in todo if here in ends[c])
            todo.remove(c)
            ring.append(here)
            here = ends[c][1] if ends[c][0] == here else ends[c][0]
        face_pts.append(ring)
    bottom = face_pts[0]
    top_set = next(set(r) for r in face_pts[1:] if not set(r) & set(bottom))
    top = [next(b if a == p else a for a, b in ends.values() if p in (a, b) and ({a, b} - {p}) <= top_set)
           for p in bottom]
    return bottom + top


def volume_mesh(info: dict, out: Path) -> dict:
    """Structured hexahedra in every block, from the step-2 surface mesh; write
    MSH 2.2 for gmshToFoam with one cell zone per mesh zone."""
    for v in info["zones"]:
        gmsh.model.mesh.setTransfiniteVolume(v, cornerTags=hex_corners(v))
        gmsh.model.mesh.setRecombine(3, v)
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.model.mesh.generate(3)

    counts, worst, low = {}, (2.0, None), 0
    for v in info["zones"]:
        for et, tags, _ in zip(*gmsh.model.mesh.getElements(3, v)):
            name = gmsh.model.mesh.getElementProperties(et)[0]
            counts[name] = counts.get(name, 0) + len(tags)
            q = gmsh.model.mesh.getElementQualities(tags, "minSJ")
            low += sum(1 for x in q if x < 0.3)
            i = min(range(len(q)), key=q.__getitem__)
            if q[i] < worst[0]:
                c = gmsh.model.mesh.getElement(tags[i])[1]
                xyz = [gmsh.model.mesh.getNode(n)[0] for n in c]
                worst = (q[i], [round(sum(p[k] for p in xyz) / len(xyz), 3) for k in range(3)])

    # Shared inner seams are internal faces of one conformal mesh: drop them
    # from the patches before writing.
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        name = gmsh.model.getPhysicalName(dim, tag)
        if name in ("seam_body_nose", "seam_aft_body"):
            gmsh.model.removePhysicalGroups([(dim, tag)])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(str(out))
    return {"counts": counts, "worst": worst, "low": low}


def write_vtk(path: Path, info: dict) -> dict[str, int]:
    """Write every meshed face as a binary legacy VTK file of quads for ParaView.

    Cell data: group (wall, seams, outer patches, internal block faces; the
    legend is returned), zone (0 nose, 1 body, 2 aft, 3 far), face (model tag)
    and minSJ (quad scaled Jacobian).
    """
    names = sorted(info["groups"]) + ["block_face"]
    group_of = {f: names.index(k) for k, fs in info["groups"].items() for f in fs}
    zone_ids = {"nose": 0, "body": 1, "aft": 2, "far": 3}
    node_ids, points, quads, cell = {}, [], [], {"group": [], "zone": [], "face": [], "minSJ": []}
    for _, f in gmsh.model.getEntities(2):
        owners = gmsh.model.getAdjacencies(2, f)[0]
        zone = min(zone_ids[info["zones"][v]] for v in owners) if len(owners) else -1
        for et, tags, nodes in zip(*gmsh.model.mesh.getElements(2, f)):
            if gmsh.model.mesh.getElementProperties(et)[3] != 4:
                continue
            q = gmsh.model.mesh.getElementQualities(tags, "minSJ")
            for e in range(len(tags)):
                ids = []
                for n in nodes[4 * e:4 * e + 4]:
                    if n not in node_ids:
                        node_ids[n] = len(points)
                        points.append(gmsh.model.mesh.getNode(n)[0])
                    ids.append(node_ids[n])
                quads.append(ids)
                cell["group"].append(group_of.get(f, len(names) - 1))
                cell["zone"].append(zone)
                cell["face"].append(f)
                cell["minSJ"].append(q[e])
    with open(path, "wb") as out:
        out.write(b"# vtk DataFile Version 3.0\nzoned surface mesh\nBINARY\nDATASET UNSTRUCTURED_GRID\n")
        out.write(f"POINTS {len(points)} double\n".encode())
        out.write(struct.pack(f">{3 * len(points)}d", *[c for p in points for c in p]))
        out.write(f"\nCELLS {len(quads)} {5 * len(quads)}\n".encode())
        out.write(struct.pack(f">{5 * len(quads)}i", *[v for qd in quads for v in (4, *qd)]))
        out.write(f"\nCELL_TYPES {len(quads)}\n".encode())
        out.write(struct.pack(f">{len(quads)}i", *([9] * len(quads))))
        out.write(f"\nCELL_DATA {len(quads)}\n".encode())
        for key in ("group", "zone", "face"):
            out.write(f"SCALARS {key} int 1\nLOOKUP_TABLE default\n".encode())
            out.write(struct.pack(f">{len(quads)}i", *cell[key]))
            out.write(b"\n")
        out.write(b"SCALARS minSJ double 1\nLOOKUP_TABLE default\n")
        out.write(struct.pack(f">{len(quads)}d", *cell["minSJ"]))
        out.write(b"\n")
    return {k: i for i, k in enumerate(names)}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--gui", action="store_true", help="open the result in the Gmsh GUI")
    ap.add_argument("--geometry-only", action="store_true", help="stop after the block layout (step 1)")
    ap.add_argument("--volume", action="store_true", help="also mesh the volume and write mesh.msh (step 3)")
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
        if not args.geometry_only:
            sm = surface_mesh(body, params, info)
            print("surface mesh:")
            print("\n".join(surface_checks(body, params, info, sm)))
            legend = write_vtk(case / "surface.vtk", info)
            print(f"wrote {case / 'surface.vtk'}; cell data group: "
                  + ", ".join(f"{i} {k}" for k, i in legend.items())
                  + "; zone: 0 nose, 1 body, 2 aft, 3 far")
            if args.volume:
                vm = volume_mesh(info, case / "mesh.msh")
                total = sum(vm["counts"].values())
                print(f"volume mesh: {total} cells {vm['counts']}; hex minSJ worst {vm['worst'][0]:.3f} "
                      f"at {vm['worst'][1]}, {vm['low']} cells below 0.3; wrote {case / 'mesh.msh'}")
        if args.gui:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
