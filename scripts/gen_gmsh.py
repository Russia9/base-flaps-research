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
import shutil
import struct
import subprocess
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
    "xSlabEnd": 900.0, "finWrap": 6.0,
    # fin slab spacing (mm)
    "hFinChord": 2.0, "hFinSpan": 1.5, "hSlabOuter": 3.0, "finSpanCut": -8.0, "hLE": 0.2,
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


def x_cylinder(x0: float, dx: float, r: float, angles: list[float] | None = None) -> tuple[int, int]:
    """Cylinder along x whose curved surface is eight 45-degree patches.

    A single occ.addCylinder has one periodic face with its parameter seam on
    a sector boundary, and structured faces that end on the seam fold (Gmsh
    picks the wrong parameter, 0 or 2 pi, at the seam corner). Here each
    sector's patch is a separate revolved, non-periodic face; the patches and
    the two end disks are sewn into one solid. `angles` (rad, increasing)
    moves the patch edges, e.g. off the fin roots.
    """
    occ = gmsh.model.occ
    a = angles or [k * math.pi / 4 for k in range(8)]
    a = a + [a[0] + 2 * math.pi]
    faces = []
    for k in range(len(a) - 1):
        line = occ.addLine(occ.addPoint(x0, r, 0), occ.addPoint(x0 + dx, r, 0))
        occ.rotate([(1, line)], 0, 0, 0, 1, 0, 0, a[k])
        faces += [t for d, t in occ.revolve([(1, line)], 0, 0, 0, 1, 0, 0, a[k + 1] - a[k]) if d == 2]
    for x in (x0, x0 + dx):
        c = occ.addPoint(x, 0, 0)
        p = [occ.addPoint(x, r * math.cos(t), r * math.sin(t)) for t in a[:-1]]
        arcs = [occ.addCircleArc(p[k], c, p[(k + 1) % len(p)]) for k in range(len(p))]
        faces.append(occ.addPlaneSurface([occ.addCurveLoop(arcs)]))
        occ.remove([(0, c)])
    return (3, occ.addVolume([occ.addSurfaceLoop(faces, sewing=True)]))


def half_plane(theta: float, x0: float, x1: float, r0: float, r1: float):
    """Rectangle x0..x1, r0..r1 in the meridional half-plane at angle theta."""
    occ = gmsh.model.occ
    s = occ.addRectangle(x0, r0, 0, x1 - x0, r1 - r0)
    occ.rotate([(2, s)], 0, 0, 0, 1, 0, 0, theta)
    return (2, s)


def annulus(x: float, r0: float, r1: float, seam: float = -math.pi / 2, seam_in: float | None = None):
    """Disk (r0 = 0) or annulus in the plane x = const, the outer circle's
    seam vertex at phi = seam and the inner one's at seam_in (rad, default
    seam)."""
    occ = gmsh.model.occ

    def disk(r, phi):
        d = occ.addDisk(0, 0, 0, r, r)
        occ.rotate([(2, d)], 0, 0, 0, 0, 1, 0, math.pi / 2)
        occ.rotate([(2, d)], 0, 0, 0, 1, 0, 0, phi + math.pi / 2)
        return d

    d = disk(r1, seam)
    if r0 > 0:
        (_, d), = occ.cut([(2, d)], [(2, disk(r0, seam if seam_in is None else seam_in))])[0]
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


def slab_tools(spec: FinSpec, P: dict, slab: tuple[float, float], body: Body) -> list:
    """Cutting surfaces for the fin slab, swept along x (see the F3 layout).

    Per fin sector (phi = -45..+45 deg about the fin root), in the fin's own
    polar frame about its arc centre: the three fin arcs and the two arcs of a
    wrap delta outside the faces; the tip-face plane, the wrap's corner lines
    and its tip offset; the tip extension, straight up the root's radial line
    to r = rZone; and the diagonals from the wrap's tip corners to
    phi = -45/+45 deg on r = rZone. Plus the sector lines between fins and the
    axial cuts where the fin walls start and end. Surfaces that would coincide
    with fin walls are swept only over the x-ranges where the fin is absent.
    """
    occ = gmsh.model.occ
    x0, x1 = slab
    total, R, r_z = body.total, body.R, P["rZone"]
    delta = P["finWrap"]
    z0, z1, z2, z3 = spec.z
    Ri, Re, Ro = spec.R_in, spec.R_edge, spec.R_out
    tip, dA = spec.th_tip, delta / Re

    def in_body(rho):
        """A fin-frame angle on radius rho a little inside the body (r < R)."""
        lo, hi = -math.pi / 2, tip
        for _ in range(200):
            m = 0.5 * (lo + hi)
            if math.hypot(*spec.point(rho, m, 0.0)[1:]) < R:
                lo = m
            else:
                hi = m
        return lo - 0.05

    def sweep(curve_pts, xa, xb, arc_centre=None):
        """Surface swept along x from xa to xb by a line (2 points) or an arc."""
        pts = [occ.addPoint(xa, y, z) for y, z in curve_pts]
        if arc_centre:
            c = occ.addPoint(xa, *arc_centre)
            curve = occ.addCircleArc(pts[0], c, pts[1])
            occ.remove([(0, c)])
        else:
            curve = occ.addLine(pts[0], pts[1])
        return [dt for dt in occ.extrude([(1, curve)], xb - xa, 0, 0) if dt[0] == 2]

    def yz(rho, a):
        return spec.point(rho, a, 0.0)[1:]

    centre = (spec.Yc, spec.Zc)
    E = (r_z, 0.0)
    Cm = (r_z * math.cos(math.pi / 4), -r_z * math.sin(math.pi / 4))
    Cp = (r_z * math.cos(math.pi / 4), r_z * math.sin(math.pi / 4))
    X = (r_z * math.cos(math.pi / 8), -r_z * math.sin(math.pi / 8))
    # Span station: every wrap is split at the same arc angle, and across the
    # gap fin 0's outer-wrap point Ws joins the next fin's inner-wrap point Vs.
    # Without the cut, the blocks between fins would carry the full span
    # count (wall stacks at the body and the tip) in both directions.
    phi_w = spec.wrap_root(delta, R)                                      # seams on r = R sit on this line
    a_s = math.radians(P["finSpanCut"])
    turn = 2 * math.pi / spec.N
    Ws = yz(Ro + delta, a_s)
    v = yz(Ri - delta, a_s)
    Vs = (v[0] * math.cos(turn) - v[1] * math.sin(turn), v[0] * math.sin(turn) + v[1] * math.cos(turn))
    Pm, Pp, Q = yz(Ri - delta, tip + dA), yz(Ro + delta, tip + dA), yz(Re, tip + dA)
    fin0 = []
    for rho in (Ri - delta, Ro + delta):                                  # the wrap
        fin0 += sweep([yz(rho, in_body(rho)), yz(rho, tip + dA)], x0, x1, centre)
    for rho, ranges in ((Ri, ((x0, z1), (z2, x1))), (Ro, ((x0, z1), (z2, x1))),
                        (Re, ((x0, z0), (z3, x1)))):                      # fin arcs off the walls
        for xa, xb in ranges:
            fin0 += sweep([yz(rho, in_body(rho)), yz(rho, tip)], xa, xb, centre)
    fin0 += sweep([yz(Re, tip), Q], x0, x1, centre)                       # mid arc past the tip
    for rho in (Ri, Ro):                                                  # tip-face plane ahead of the wall
        quad = [(x0, *yz(rho, tip)), (x0, *yz(Re, tip)), (z0, *yz(Re, tip)), (z1, *yz(rho, tip))]
        pts = [occ.addPoint(*q) for q in quad]
        loop = occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)])
        fin0.append((2, occ.addPlaneSurface([loop])))
    fin0 += sweep([yz(Ri, tip), Pm], x0, x1) + sweep([yz(Ro, tip), Pp], x0, x1)  # wrap corner lines
    fin0 += sweep([Pm, Pp], x0, x1)                                       # wrap tip offset
    # Outer lines from the cap top to r = rZone: Pm to -45 deg, the mid arc's
    # end Q to -22.5 deg, Pp to 0 deg. (Pp to +45 deg would give the block
    # above the cap's R_out half a 237-degree corner at Pp: the fin's bow
    # tilts the tip towards -Z.)
    fin0 += sweep([Q, X], x0, x1)                                         # tip extension
    fin0 += sweep([Pm, Cm], x0, x1) + sweep([Pp, E], x0, x1)              # diagonals
    fin0 += sweep([yz(Ri - delta, a_s), Ws], x0, x1)                      # span station through the wrap
    fin0 += sweep([Ws, Vs], x0, x1) + sweep([Ws, Cp], x0, x1)             # across the gap, and up to +45 deg
    # Leading-edge funnel (see the "Fin Tip Funnel" layout). The cap is not
    # cut at z1, so the tip wall from the LE to z2 is one face. The split in
    # front of it leaves the wall at the LE edge: ruled through the cap from
    # the mid arc at z0 to the wrap corner line at z1, then bent by a second
    # ruled surface so it reaches r = rZone on the x = z1 circle.
    for rho, Pc, Cc in ((Ri, Pm, Cm), (Ro, Pp, E)):
        fin0 += ruled(("arc", (z0, *yz(Re, tip)), (z0, *Q), (z0, *centre)),
                      ("line", (z1, *yz(rho, tip)), (z1, *Pc)))
        fin0 += ruled(("line", (z0, *Q), (z1, *Pc)), ("arc", (z1, *X), (z1, *Cc), (z1, 0.0, 0.0)))
    tools = list(fin0)
    for k in range(1, spec.N):
        copy = occ.copy(fin0)
        occ.rotate(copy, 0, 0, 0, 1, 0, 0, k * 2 * math.pi / spec.N)
        tools += copy
    # The z1 plane stops at the cap and the blocks above it, where the funnel
    # takes over.
    # Circle seams: on r = rZone at 0 deg (a block corner), on r = R at the
    # outer wrap's root.
    plane = annulus(z1, 0.0, r_z, 0.0)
    out = 1.3
    # Its corners on the fin side run on past the tip corners into the fin,
    # so its edge across the tip lies inside the solid and does not imprint
    # x = z1 on the tip wall.
    def past(a, b):
        return (b[0] + 0.25 * (b[0] - a[0]), b[1] + 0.25 * (b[1] - a[1]))

    notch_yz = [past(Pm, yz(Ri, tip)), Pm, Cm, (out * Cm[0], out * Cm[1]), (out * X[0], out * X[1]),
                (out * r_z, 0.0), E, Pp, past(Pp, yz(Ro, tip))]
    notches = []
    for k in range(spec.N):
        c, sn = math.cos(k * 2 * math.pi / spec.N), math.sin(k * 2 * math.pi / spec.N)
        pts = [occ.addPoint(z1, y * c - z * sn, y * sn + z * c) for y, z in notch_yz]
        loop = occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))])
        notches.append((2, occ.addPlaneSurface([loop])))
    plane = occ.cut([plane], notches)[0]
    tools += plane + [annulus(z2, 0.0, r_z, 0.0), annulus(total, R, r_z, 0.0, phi_w)]
    return tools


def ruled(w0: tuple, w1: tuple) -> list:
    """Ruled surface between two curves, each ("line", A, B) or ("arc", A, B,
    centre) with 3-D points, ruled A0-A1 and B0-B1. OCC's ThruSections picks
    the pairing itself and can return the crossed (twisted) surface; this
    checks the side edges and rebuilds with the second curve reversed."""
    occ = gmsh.model.occ

    def curve(w, rev):
        kind, a, b = w[0], w[1], w[2]
        if rev:
            a, b = b, a
        pa, pb = occ.addPoint(*a), occ.addPoint(*b)
        if kind == "line":
            return occ.addLine(pa, pb)
        c = occ.addPoint(*w[3])
        t = occ.addCircleArc(pa, c, pb)
        occ.remove([(0, c)])
        return t

    def near(p, q):
        return math.dist(p, q) < 1e-6

    for rev in (False, True):
        faces = [dt for dt in occ.addThruSections([occ.addWire([curve(w0, False)]), occ.addWire([curve(w1, rev)])],
                                                  makeSolid=False, makeRuled=True) if dt[0] == 2]
        occ.synchronize()
        ends = []
        for _, c in gmsh.model.getBoundary(faces, oriented=False):
            pts = [gmsh.model.getValue(0, p, []) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)]
            if len(pts) == 2:
                ends.append(pts)
        crossed = any((near(p, w0[1]) and near(q, w1[2])) or (near(q, w0[1]) and near(p, w1[2]))
                      for p, q in ends)
        if not crossed:
            return faces
        occ.remove(faces, recursive=True)
    raise RuntimeError(f"ruled surface between {w0} and {w1} stays twisted")


def cap_boxes(spec: FinSpec, P: dict, x1: float) -> list:
    """The tip caps behind z2 (one per fin): their own zone, so the prism
    chain from the trailing-edge root triangle ends on a flat seam."""
    occ = gmsh.model.occ
    z2, c = spec.z[2], spec.cap(P["finWrap"])
    pts = [occ.addPoint(z2, *c[k]) for k in ("bi", "bo", "Pp", "Pm")]
    quad = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)])])
    box = [dt for dt in occ.extrude([(2, quad)], x1 - z2, 0, 0) if dt[0] == 3]
    out = list(box)
    for k in range(1, spec.N):
        copy = occ.copy(box)
        occ.rotate(copy, 0, 0, 0, 1, 0, 0, k * 2 * math.pi / spec.N)
        out += copy
    return out


def cap_volume(spec: FinSpec, P: dict, x1: float) -> float:
    """Analytic volume of the N tip caps (straight-sided quadrilateral section)."""
    q = [spec.cap(P["finWrap"])[k] for k in ("bi", "bo", "Pp", "Pm")]
    area = 0.5 * abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(q, q[1:] + q[:1])))
    return spec.N * area * (x1 - spec.z[2])


def cap_tools(spec: FinSpec, P: dict, x1: float) -> list:
    """Split each cap at the mid arc and along the trailing-edge wedge line:
    per half, a prism over the wall triangle and a hexahedron behind it."""
    occ = gmsh.model.occ
    z2, z3, c = spec.z[2], spec.z[3], spec.cap(P["finWrap"])
    centre = (spec.Yc, spec.Zc)

    def arc(x, a, b):
        cc = occ.addPoint(x, *centre)
        t = occ.addCircleArc(occ.addPoint(x, *a), cc, occ.addPoint(x, *b))
        occ.remove([(0, cc)])
        return t

    mid = arc(z2, c["e"], c["Q"])
    fin0 = [dt for dt in occ.extrude([(1, mid)], x1 - z2, 0, 0) if dt[0] == 2]
    for b, Pc in (("bi", "Pm"), ("bo", "Pp")):
        fin0 += ruled(("line", (z2, *c[b]), (z2, *c[Pc])), ("arc", (z3, *c["e"]), (z3, *c["Q"]), (z3, *centre)))
    # The base plane through the cap, as through the wrap layers beside it:
    # its bottom, corner and front faces then match the slab's face for face
    # (stitched as perfect matches, clear of the 3 um cells at the fin walls)
    # and the trailing-edge prism chain runs on to the cap's top seam.
    pts = [occ.addPoint(z3, *c[k]) for k in ("bi", "bo", "Pp", "Pm")]
    fin0.append((2, occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % 4])
                                                           for i in range(4)])])))
    tools = list(fin0)
    for k in range(1, spec.N):
        copy = occ.copy(fin0)
        occ.rotate(copy, 0, 0, 0, 1, 0, 0, k * 2 * math.pi / spec.N)
        tools += copy
    return tools


def cap_piece(s: int, spec: FinSpec, P: dict, x1: float):
    """(fin, side) if face s lies on a tip cap's boundary behind z2, else None.
    Sides: bottom (tip plane), top (wrap offset), ci/co (corner planes),
    front (x = z2)."""
    tol = 1e-3
    bb = gmsh.model.getBoundingBox(2, s)
    z2 = spec.z[2]
    if bb[0] < z2 - tol or bb[3] > x1 + tol:
        return None
    x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
    k = round(math.atan2(z, y) / (2 * math.pi / spec.N)) % spec.N
    a = -k * 2 * math.pi / spec.N
    y, z = y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)
    c = spec.cap(P["finWrap"])

    def on(p, q):
        (px, py), (qx, qy) = p, q
        L2 = (qx - px) ** 2 + (qy - py) ** 2
        t = max(0.0, min(1.0, ((y - px) * (qx - px) + (z - py) * (qy - py)) / L2))
        return math.hypot(y - px - t * (qx - px), z - py - t * (qy - py)) < tol

    for side, (p, q) in (("bottom", ("bi", "bo")), ("top", ("Pm", "Pp")), ("ci", ("bi", "Pm")),
                         ("co", ("bo", "Pp"))):
        if bb[3] - bb[0] > tol and on(c[p], c[q]):
            return k, side
    if bb[3] - bb[0] < tol and abs(x - z2) < tol:
        quad = [c[n] for n in ("bi", "bo", "Pp", "Pm")]
        cross = [(q[0] - p[0]) * (z - p[1]) - (q[1] - p[1]) * (y - p[0])
                 for p, q in zip(quad, quad[1:] + quad[:1])]
        if all(v > 0 for v in cross) or all(v < 0 for v in cross):
            return k, "front"
    return None


def build(body: Body, P: dict, spec: FinSpec | None = None) -> dict:
    """Zones, split into hexahedral blocks, and the named faces.

    Zones sharing a group are fragmented together and share their faces;
    different groups are built separately, so the seams between them are two
    independent copies of the same surface, stitched non-conformally later.
    The bare body (no fins) uses groups {nose, body, aft} and {far}. With fins,
    a slab x = xAftStart..xSlabEnd holds the fins and the base, from the axis
    to the far field, with its own circumferential count: {nose, body},
    {slab}, {wake}, {far}, {far_slab}, {far_aft}. The count changes only across
    the slab's x-planes, so every cylindrical seam keeps the same count on
    both sides.
    """
    occ = gmsh.model.occ
    R, total = body.R, body.total
    x_in, x_out, r_ff, r_z = P["xInlet"], P["xOutlet"], P["rFar"], P["rZone"]
    xs = [P["xZoneStart"], P["xNoseEnd"], P["xAftStart"], P["xZoneEnd"]]
    a, s_w, s_f = P["coreTip"], P["coreWake"], P["coreFar"]
    thetas = [k * math.pi / 4 for k in range(8)]
    slab = (xs[2], P["xSlabEnd"]) if spec else None
    if slab and not total < slab[1] < xs[3]:
        raise ValueError("xSlabEnd must lie between the base and xZoneEnd")

    # Body: revolve the closed profile (apex, ogive arc, cylinder, base, axis).
    p_apex = occ.addPoint(0, 0, 0)
    p_sh = occ.addPoint(body.l_ogive, R, 0)
    p_rim = occ.addPoint(total, R, 0)
    p_base = occ.addPoint(total, 0, 0)
    p_c = occ.addPoint(body.xc, body.yc, 0)
    p_sh0 = occ.addPoint(body.l_ogive, 0, 0)
    cut = occ.addLine(p_sh, p_sh0)
    nose = occ.addPlaneSurface([occ.addCurveLoop([occ.addCircleArc(p_apex, p_c, p_sh), cut,
                                                  occ.addLine(p_sh0, p_apex)])])
    # A revolution's parameter seam starts at phi = 0, where the nose and
    # body blocks have an edge. In the fin slab phi = 0 is fin 0's root,
    # where the seam would split the faces around the root, so the body
    # behind the slab start is a separate revolution turned onto the outer
    # wrap's root line (a block edge on r = R there). The bare body keeps one
    # cylinder, turned onto 45 deg.
    x_cut = slab[0] if slab else None
    profiles = [nose]
    if x_cut:
        p_a, p_a0 = occ.addPoint(x_cut, R, 0), occ.addPoint(x_cut, 0, 0)
        split = occ.addLine(p_a, p_a0)
        profiles.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_sh, p_a), split,
                                                               occ.addLine(p_a0, p_sh0), cut])]))
        profiles.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_a, p_rim), occ.addLine(p_rim, p_base),
                                                               occ.addLine(p_base, p_a0), split])]))
        turns = [0.0, 0.0, spec.wrap_root(P["finWrap"], R)]
    else:
        profiles.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_sh, p_rim), occ.addLine(p_rim, p_base),
                                                               occ.addLine(p_base, p_sh0), cut])]))
        turns = [0.0, math.pi / 4]
    parts = []
    for f, turn in zip(profiles, turns):
        vol = [(3, t) for d, t in occ.revolve([(2, f)], 0, 0, 0, 1, 0, 0, 2 * math.pi) if d == 3]
        if turn:
            occ.rotate(vol, 0, 0, 0, 1, 0, 0, turn)
        parts += vol
    occ.remove([(0, p_c)])
    v_body = sum(occ.getMass(*v) for v in parts)
    fin_info = {}
    if spec:
        # The pieces stay separate solids (a fuse merges the two cylinders'
        # faces and with them the seams); the cuts below only use them as
        # tools. The fused copy is just for the wall checks.
        fins = build_fins(spec)
        whole, _ = occ.fuse(occ.copy(parts), occ.copy(fins))
        occ.synchronize()
        fin_info = {"exposed_volume": occ.getMass(*whole[0]) - v_body,
                    "wall_area": sum(occ.getMass(*f) for f in gmsh.model.getBoundary(whole, oriented=False))}
        occ.remove(whole, recursive=True)
        solid = parts + fins
    else:
        solid, _ = occ.fuse(parts[:1], parts[1:])

    # Zone layout: (name, x0, x1, inner?, group).
    if slab:
        layout = [("nose", xs[0], xs[1], True, 0), ("body", xs[1], slab[0], True, 0),
                  ("slab", slab[0], slab[1], True, 1), ("slab_core", total, slab[1], True, 6),
                  ("cap", spec.z[2], slab[1], True, 7),
                  ("wake", slab[1], xs[3], True, 2),
                  ("far", x_in, slab[0], False, 3), ("far_slab", slab[0], slab[1], False, 4),
                  ("far_aft", slab[1], x_out, False, 5)]
    else:
        layout = [("nose", xs[0], xs[1], True, 0), ("body", xs[1], xs[2], True, 0),
                  ("aft", xs[2], xs[3], True, 0), ("far", x_in, x_out, False, 1)]
    def overlapping(x0, x1):
        """The body and fin solids reaching into x0..x1 (a solid that only
        touches the range would imprint its seam vertices on the zone)."""
        out = []
        for v in solid:
            bb = occ.getBoundingBox(*v)
            if bb[0] < x1 - PAD and bb[3] > x0 + PAD:
                out.append(v)
        return out

    zones = {}
    for name, x0, x1, inner, _ in layout:
        if name == "slab_core":
            # Behind the base, r < R: its own zone, so the fin sectors' many
            # block corners on r = R stay out of the core (joined with
            # matching nodes on r = R).
            pieces = [x_cylinder(x0, x1 - x0, R, [math.pi / 4 + k * math.pi / 2 for k in range(4)])]
        elif name == "slab":
            phi_w = spec.wrap_root(P["finWrap"], R)
            core = x_cylinder(total, x1 - total, R,
                              [phi_w + k * 2 * math.pi / spec.N for k in range(spec.N)])
            caps = cap_boxes(spec, P, x1)
            pieces, _ = occ.cut([x_cylinder(x0, x1 - x0, r_z)], overlapping(x0, x1) + [core] + caps,
                                removeTool=False)
            occ.remove([core], recursive=True)
        elif name == "cap":
            pieces = caps
        elif inner:
            cyl, tools = x_cylinder(x0, x1 - x0, r_z), overlapping(x0, x1)
            pieces = occ.cut([cyl], tools, removeTool=False)[0] if tools else [cyl]
        else:
            hole = max(x0, xs[0]), min(x1, xs[3])
            tools = [x_cylinder(hole[0], hole[1] - hole[0], r_z)] if hole[1] > hole[0] else []
            pieces, _ = occ.cut([x_cylinder(x0, x1 - x0, r_ff)], tools) if tools else ([x_cylinder(
                x0, x1 - x0, r_ff)], None)
        zones |= {v: name for d, v in pieces if d == 3}
    occ.remove(solid, recursive=True)
    occ.synchronize()
    group_of = {name: g for name, *_, g in layout}

    def of_group(g):
        return {v: n for v, n in zones.items() if group_of[n] == g}

    def others(g):
        return {v: n for v, n in zones.items() if group_of[n] != g}

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

    # Group 0 (nose, body and, without fins, aft): 8 sectors; the 45-degree
    # cuts stop at the core corners. Upstream of the base the cut follows the
    # core's corner line to x_end, then runs inside the body; it meets the
    # base plane only at the rim, so the base disk is not cut. Behind the base
    # it stops at the wake core's corner.
    x_g0 = slab[0] if slab else xs[3]
    tools = [half_plane(t, xs[0], x_g0, 0.0, r_z) for t in thetas[0::2]]
    for t in thetas[1::2]:
        if slab:
            tools.append(polygon_half_plane(t, [(xs[0], SQ2 * h(xs[0])), (x_end, SQ2 * h(x_end)),
                                                (x_end, 0.0), (x_g0, 0.0), (x_g0, r_z), (xs[0], r_z)]))
        else:
            tools += [polygon_half_plane(t, [(xs[0], SQ2 * h(xs[0])), (x_end, SQ2 * h(x_end)),
                                             (x_end, 0.0), (total - 0.5 * R, 0.0), (total, R),
                                             (total, r_z), (xs[0], r_z)]),
                      half_plane(t, total, xs[3], SQ2 * s_w, r_z)]
    tools += funnel + [annulus(body.l_ogive, 0.0, r_z)]          # shoulder
    if not slab:
        tools += [square_prism(total, s_w, xs[3], s_w),
                  x_cylinder(total, xs[3] - total, R),
                  annulus(total, R, r_z),                         # base plane, outside the base disk
                  annulus(P["xWakeSplit"], 0.0, r_z)]             # wake: near part | relaxing part
    zones = others(0) | fragment(of_group(0), tools)

    def wake_tools(x0, x1):
        """Core square, inner ring (to r = R), outer ring (to r_z) in 8 sectors."""
        return ([half_plane(t, x0, x1, 0.0, r_z) for t in thetas[0::2]]
                + [half_plane(t, x0, x1, SQ2 * s_w, r_z) for t in thetas[1::2]]
                + [square_prism(x0, s_w, x1, s_w), x_cylinder(x0, x1 - x0, R)])

    def far_tools(x0, x1):
        """Far zone between x0 and x1: core + two rings outside the inner zones'
        x-range, one ring around them."""
        t_ = [half_plane(t, x0, x1, 0.0, r_ff) for t in thetas[0::2]]
        for t in thetas[1::2]:
            for a0, a1, r0 in ((x0, min(x1, xs[0]), SQ2 * s_f), (max(x0, xs[0]), min(x1, xs[3]), r_z),
                               (max(x0, xs[3]), x1, SQ2 * s_f)):
                if a1 > a0:
                    t_.append(half_plane(t, a0, a1, r0, r_ff))
        for a0, a1 in ((x0, min(x1, xs[0])), (max(x0, xs[3]), x1)):
            if a1 > a0:
                t_ += [square_prism(a0, s_f, a1, s_f), x_cylinder(a0, a1 - a0, r_z)]
        for xp in (xs[0], xs[3]):
            if x0 < xp < x1:
                t_.append(annulus(xp, r_z, r_ff))
        return t_

    if slab:
        zones = others(1) | fragment(of_group(1), slab_tools(spec, P, slab, body))
        # One square and four ring blocks (corners at 45 + 90k deg): the
        # square's opposite sides then match under the fins' 90-degree
        # symmetry, and the ring arcs copy the slab's nodes on r = R.
        core = ([half_plane(t, total, slab[1], SQ2 * s_w, R) for t in thetas[1::2]]
                + [square_prism(total, s_w, slab[1], s_w)])
        zones = others(6) | fragment(of_group(6), core)
        zones = others(7) | fragment(of_group(7), cap_tools(spec, P, slab[1]))
        zones = others(2) | fragment(of_group(2), wake_tools(slab[1], xs[3]))
        for g, (x0, x1) in ((3, (x_in, slab[0])), (4, slab), (5, (slab[1], x_out))):
            tools_g = far_tools(x0, x1)
            if g == 4:
                # The slab's x-stations, so the far slab's axial lines on
                # r = rZone match the slab's segments and the seam's nodes
                # coincide (stitched as a perfect match).
                tools_g += [annulus(xp, r_z, r_ff) for xp in (spec.z[1], spec.z[2], total)]
            zones = others(g) | fragment(of_group(g), tools_g)
    else:
        zones = others(1) | fragment(of_group(1), far_tools(x_in, x_out))

    # Faces: wall (body, fins), outer patches, shared faces between zones of
    # one group, and stitched seams between groups. A seam is named by its
    # pair and side; stitchMesh couples each smooth piece separately (across
    # the 90-degree edges between a cylinder and a disk its projection fails),
    # with the coarser side as master.
    owners: dict[int, list[int]] = {}
    for v in zones:
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            owners.setdefault(s, []).append(v)
    planes = {xs[0]: "up", xs[3]: "down"}
    if slab:
        planes |= {slab[0]: "slab_in", slab[1]: "slab_out"}
    fine = {"slab", "far_slab", "slab_core", "cap"}
    groups: dict[str, list[int]] = {}
    shared = set()
    for s, own in owners.items():
        names = sorted({zones[v] for v in own})
        if len(own) == 2 and len(names) == 1:
            continue                                   # internal block face
        if len(names) == 2:
            key = "shared_" + "_".join(names)          # conformal face between zones of a group
            shared.add(key)
            groups.setdefault(key, []).append(s)
            continue
        zone = names[0]
        bb = gmsh.model.getBoundingBox(2, s)
        rmax = max(abs(v) for v in bb[1:3] + bb[4:6])
        flat_x = abs(bb[0] - bb[3]) < PAD
        cp = cap_piece(s, spec, P, slab[1]) if slab else None
        if cp:
            # Tip cap behind z2: flat seams to the slab, apart from the
            # trailing-edge wall triangles on its bottom.
            k, side = cp
            x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
            a = -k * 2 * math.pi / spec.N
            yr, zr = y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)
            rho = math.hypot(yr - spec.Yc, zr - spec.Zc)
            z2, z3 = spec.z[2], spec.z[3]
            if side == "bottom" and x < z3 and abs(rho - spec.R_edge) < (
                    (spec.R_edge - spec.R_in) * (z3 - x) / (z3 - z2)):
                key = "stabilizers"
            else:
                key = f"seam_cap{k}_{side}_{'slave' if zone == 'cap' else 'master'}"
            groups.setdefault(key, []).append(s)
            continue
        plane = next((n for xp, n in planes.items() if flat_x and abs(bb[0] - xp) < PAD), None)
        far_side = zone.startswith("far")
        if plane and (rmax < r_z + PAD or plane.startswith("slab")):
            if plane.startswith("slab") and rmax > r_z + PAD:
                plane = "f" + plane                    # far-zone part of a slab plane
            if plane in ("up", "down"):
                role = "master" if far_side else "slave"
            else:
                role = "slave" if zone in fine else "master"
            key = f"seam_{plane}_{role}"
        elif slab and bb[0] > total - PAD and bb[3] < slab[1] + PAD and _on_cylinder(s, R):
            key = f"seam_core_side_{'master' if zone == 'slab_core' else 'slave'}"
        elif rmax < r_z + PAD and bb[0] > xs[0] - PAD and bb[3] < xs[3] + PAD and _on_cylinder(s, r_z):
            part = "side"
            if slab:
                mid = 0.5 * (bb[0] + bb[3])
                part = "side_fore" if mid < slab[0] else ("side_slab" if mid < slab[1] else "side_aft")
            key = f"seam_{part}_{'master' if far_side else 'slave'}"
        elif bb[3] < x_in + PAD:
            key = "inlet"
        elif bb[0] > x_out - PAD:
            key = "outlet"
        elif rmax > r_ff - PAD and _on_cylinder(s, r_ff):
            key = "farfield"
        elif spec and not _on_body(s, body):
            key = "stabilizers"
        else:
            key = "fuselage"
        groups.setdefault(key, []).append(s)

    for name, *_ in layout:
        gmsh.model.addPhysicalGroup(3, [v for v, n in zones.items() if n == name], name=name)
    for key, surfs in groups.items():
        gmsh.model.addPhysicalGroup(2, surfs, name=key)
    return {"zones": zones, "groups": groups, "xs": xs, "layout": layout, "shared": shared,
            "slab": slab, "fins": fin_info, "body_volume": v_body, "fins_n": spec.N if spec else 0,
            "cap_volume": cap_volume(spec, P, slab[1]) if slab else 0.0,
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


def _on_body(s: int, body: Body) -> bool:
    """True if face s lies on the body (ogive, cylinder or base disk)."""
    (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, s)
    uv = [c for f in (0.21, 0.5, 0.79) for g in (0.23, 0.5, 0.77)
          for c in (u0 + (u1 - u0) * f, v0 + (v1 - v0) * g)]
    pts = gmsh.model.getValue(2, s, uv)
    seen = 0
    for i in range(0, len(pts), 3):
        p = pts[i:i + 3]
        if not gmsh.model.isInside(2, s, p, parametric=False):
            continue                    # a trimmed face's parameter box reaches past the face
        seen += 1
        x, r = p[0], math.hypot(p[1], p[2])
        on_base = abs(x - body.total) < 1e-6 and r < body.R + 1e-6
        if not on_base and abs(r - body.r(x)) > 1e-4:
            return False
    return seen > 0


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
        try:
            return self.h * self.q_layer ** min(k, self.n_layers) * self.q ** max(0, k - self.n_layers)
        except OverflowError:
            return math.inf


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


# ── generic spacing: edge families and cell counts (F5) ──────────────────────

WALLS = ("fuselage", "stabilizers")


class Topo:
    """The block topology of the built model: every block edge with its end
    vertices (parametric start first), the zones it borders, its direction
    class, and the wall and seam faces at each vertex."""

    def __init__(self, info: dict):
        self.info = info
        zones = info["zones"]
        self.face_zones: dict[int, set] = {}
        for v, name in zones.items():
            for _, f in gmsh.model.getBoundary([(3, v)], oriented=False):
                self.face_zones.setdefault(f, set()).add(name)
        self.group_of = {s: k for k, ss in info["groups"].items() for s in ss}
        self.loops: dict[int, list[int]] = {}
        self.curves: dict[int, Curve] = {}
        self.ends: dict[int, tuple[int, int]] = {}
        self.edge_zones: dict[int, set] = {}
        self.xyz: dict[int, list[float]] = {}
        for f, names in self.face_zones.items():
            loop = [c for _, c in real_edges((2, f))]
            self.loops[f] = loop
            for c in loop:
                self.edge_zones.setdefault(c, set()).update(names)
                if c in self.curves:
                    continue
                cv = Curve(c)
                pts = [p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)]
                for p in pts:
                    self.xyz.setdefault(p, gmsh.model.getValue(0, p, []))
                start = min(pts, key=lambda p: math.dist(self.xyz[p], cv.p))
                end = next((p for p in pts if p != start), start)
                self.curves[c], self.ends[c] = cv, (start, end)
        # wall and seam faces at each vertex
        self.at_vertex: dict[int, list[int]] = {}
        for f in self.face_zones:
            key = self.group_of.get(f, "")
            if key in WALLS or key.startswith("seam_"):
                for c in self.loops[f]:
                    for p in self.ends[c]:
                        self.at_vertex.setdefault(p, [])
                        if f not in self.at_vertex[p]:
                            self.at_vertex[p].append(f)

    def cycle(self, f: int) -> tuple[list[int], list[int], list[int]]:
        """Face f's corners in loop order, its edges between them, and each
        edge's orientation along the loop (+1: parametric start first)."""
        loop = list(self.loops[f])
        here = self.ends[loop[0]][0]
        corners, edges, signs = [], [], []
        while loop:
            c = next(c for c in loop if here in self.ends[c])
            loop.remove(c)
            corners.append(here)
            edges.append(c)
            a, b = self.ends[c]
            signs.append(1 if a == here else -1)
            here = b if a == here else a
        return corners, edges, signs

    def direction(self, c: int) -> str:
        cv = self.curves[c]
        L = max(cv.length, 1e-12)
        dx, dr = abs(cv.q[0] - cv.p[0]), abs(cv.rq - cv.rp)
        if dx > 0.9 * L:
            return "axial"
        if dx < 1e-6 * max(1.0, L):
            if dr > 0.9 * L:
                return "radial"
            if dr < 0.05 * L and min(cv.rp, cv.rq) > 1e-3:
                return "circ"
            return "cross"
        return "radial" if dr > dx else "axial"     # dominant direction, as classify()

    def tangent(self, c: int, p: int) -> list[float]:
        """Unit tangent of edge c at its end p, pointing into the edge."""
        pts = self.curves[c].pts
        a, b = (pts[0], pts[3]) if p == self.ends[c][0] else (pts[-1], pts[-4])
        d = [b[i] - a[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in d)) or 1.0
        return [x / n for x in d]

    def normal_faces(self, c: int, p: int, keys) -> list[int]:
        """Faces in the groups `keys` at vertex p that edge c leaves at a steep
        angle (|cos| > 0.8 between the edge and the face normal, within 37
        deg of the normal), excluding faces the edge lies in. A 45-degree
        wedge apex (the fin LE/TE) is not steep: stacks along the wedge
        into its apex fold the cells there."""
        out, t = [], self.tangent(c, p)
        for f in self.at_vertex.get(p, []):
            key = self.group_of.get(f, "")
            if not any(key == k or key.startswith(k) for k in keys) or c in self.loops[f]:
                continue
            uv = gmsh.model.getParametrization(2, f, self.xyz[p])
            n = gmsh.model.getNormal(f, uv)
            if abs(sum(n[i] * t[i] for i in range(3))) > 0.8:
                out.append(f)
        return out


class Families:
    """Union-find over block edges with orientation: edges in one family share
    a cell count, and their parametric directions are related by a sign."""

    def __init__(self):
        self.parent: dict[int, int] = {}
        self.sign: dict[int, int] = {}
        self.conflicts: list[tuple[int, int]] = []
        self.links: list[tuple[int, int, int]] = []    # coincident (master, slave, sign) across seams

    def find(self, c: int) -> tuple[int, int]:
        if self.parent.setdefault(c, c) == c:
            self.sign.setdefault(c, 1)
            return c, 1
        root, s = self.find(self.parent[c])
        self.sign[c] *= s
        self.parent[c] = root
        return root, self.sign[c]

    def union(self, a: int, b: int, s: int):
        """b runs in direction s (+1 same, -1 opposite) relative to a."""
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            if pa * pb != s:
                self.conflicts.append((a, b))
            return
        self.parent[rb], self.sign[rb] = ra, s * pa * pb

    def groups(self) -> dict[int, list[int]]:
        out: dict[int, list[int]] = {}
        for c in list(self.parent):
            out.setdefault(self.find(c)[0], []).append(c)
        return out


def corner_angle(p, a, b) -> float:
    """Angle at p between the directions to a and b (rad)."""
    u, w = [a[j] - p[j] for j in range(3)], [b[j] - p[j] for j in range(3)]
    return math.acos(max(-1.0, min(1.0, sum(u[j] * w[j] for j in range(3)) / math.dist(a, p) / math.dist(b, p))))


def edge_families(topo: Topo) -> Families:
    """Opposite edges of every quad face; the two edges at a triangle's
    collapsed corner (the widest one, as prism_corners collapses it; the apex
    for faces with OCC's pole edge); coincident edges across a seam."""
    fam = Families()
    for c in topo.curves:
        fam.find(c)
    for f, loop in topo.loops.items():
        corners, edges, signs = topo.cycle(f)
        if len(edges) == 4:
            for i in (0, 1):
                fam.union(edges[i], edges[i + 2], -signs[i] * signs[i + 2])
        elif len(edges) == 3:
            all_edges = gmsh.model.getBoundary([(2, f)], oriented=False)
            if len(all_edges) != len(loop):            # pole edge: collapsed at the apex
                apex = min(corners, key=lambda p: math.hypot(*topo.xyz[p][1:]) + abs(topo.xyz[p][0]))
                k = corners.index(apex)
            else:
                pts = [topo.xyz[p] for p in corners]
                k = max(range(3), key=lambda i: corner_angle(pts[i], pts[i - 1], pts[(i + 1) % 3]))
            out_e, in_e = edges[k], edges[k - 1]       # leaves corner k / arrives at it
            fam.union(out_e, in_e, -signs[k] * signs[k - 1])
    # seams: a slave edge coinciding with a master edge (same ends and midpoint)
    def key(p):
        return tuple(round(x, 6) for x in p)
    by_pair: dict[str, dict[str, set]] = {}
    for f in topo.face_zones:
        g = topo.group_of.get(f, "")
        if g.startswith("seam_") and g.rsplit("_", 1)[1] in ("master", "slave"):
            pair, role = g[5:].rsplit("_", 1)
            by_pair.setdefault(pair, {}).setdefault(role, set()).update(topo.loops[f])
    # Cylindrical seams: axial edges on the two sides spanning the same x
    # range share a family (same axial nodes), also where the sides' block
    # corners differ in phi (the core seam behind the base).
    for pair, sides in by_pair.items():
        if not (pair.startswith("side") or pair == "core_side"):
            continue
        spans: dict[str, dict[tuple, int]] = {}
        for role, edges in sides.items():
            for c in edges:
                cv = topo.curves[c]
                if topo.direction(c) != "axial":
                    continue
                x0, x1 = sorted((cv.p[0], cv.q[0]))
                spans.setdefault(role, {}).setdefault((round(x0, 6), round(x1, 6), round(cv.rp, 6)), c)
        for k, cm in spans.get("master", {}).items():
            cs = spans.get("slave", {}).get(k)
            if cs:
                same = (topo.curves[cm].p[0] < topo.curves[cm].q[0]) == (topo.curves[cs].p[0] < topo.curves[cs].q[0])
                fam.union(cm, cs, 1 if same else -1)
    for pair, sides in by_pair.items():
        # Only seams that need matching nodes link families: the cylindrical
        # ones (same circumferential nodes on both sides) and the tip caps
        # (node-conformal flat seams). The other flat seams take any mesh.
        if not (pair.startswith(("side", "cap")) or pair == "core_side"):
            continue
        index = {}
        for c in sides.get("master", ()):
            cv = topo.curves[c]
            index[(key(cv.p), key(cv.q), key(cv.pts[200]))] = (c, 1)
            index[(key(cv.q), key(cv.p), key(cv.pts[200]))] = (c, -1)
        for c in sides.get("slave", ()):
            cv = topo.curves[c]
            hit = index.get((key(cv.p), key(cv.q), key(cv.pts[200])))
            if hit:
                fam.union(hit[0], c, hit[1])
                fam.links.append((hit[0], c, hit[1]))
    return fam


def edge_needs(topo: Topo, body: Body, P: dict, spec) -> dict[int, tuple]:
    """Per edge: (end condition at its start, at its end, cap on the cell size).
    An end condition is the wall stack where the edge leaves a wall steeply
    (the 20 um tip-patch cell near the apex), the seam size where it leaves a
    seam steeply (far zones: hSeam, inner zones on r = rZone: hRingOut), the
    shoulder size on axial edges at the ogive-cylinder kink, or None (free)."""
    q = P["growth"]
    wall = End(P["firstLayer"], q, P["nLayers"], P["layerRatio"])
    tip = End(P["firstLayerTip"], q)
    x45 = topo.info["tip"]["x_45"]
    r_z, R, total = P["rZone"], body.R, body.total
    r_tip = spec.tip_radius() + 2 * P["finWrap"] if spec else 0.0         # axial fin sizes inside this
    fin_x = (spec.z[0] - P["finWrap"], total) if spec else None
    def in_wrap(pt):
        """Inside a fin's wrap (fin frame: between the wrap arcs, up to the
        wrap's tip offset), where the fin sizes apply."""
        if not spec:
            return False
        _, y, z = pt
        k = round(math.atan2(z, y) / (2 * math.pi / spec.N))
        a = -k * 2 * math.pi / spec.N
        y, z = y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)
        rho = math.hypot(y - spec.Yc, z - spec.Zc)
        ang = math.atan2(y - spec.Yc, z - spec.Zc)               # spec.point's angle convention
        d, tol = P["finWrap"], 1e-3
        return (spec.R_in - d - tol <= rho <= spec.R_out + d + tol
                and ang <= spec.th_tip + d / spec.R_edge + tol)

    def on_fin_edge(p, c):
        """Vertex p on the fin's LE or TE edge (x = z0 or z3) and on a fin wall
        face that edge c does not lie in (the wedges are too shallow for the
        wall rule)."""
        x = topo.xyz[p][0]
        if not (abs(x - spec.z[0]) < 1e-6 or abs(x - spec.z[3]) < 1e-6):
            return False
        return any(topo.group_of.get(f) == "stabilizers" and c not in topo.loops[f]
                   for f in topo.at_vertex.get(p, []))

    def on_te_line(p):
        """Vertex p on a fin trailing edge's downstream line (fin frame:
        R_edge, between the root and the tip)."""
        _, y, z = topo.xyz[p]
        k = round(math.atan2(z, y) / (2 * math.pi / spec.N))
        a = -k * 2 * math.pi / spec.N
        y, z = y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)
        rho = math.hypot(y - spec.Yc, z - spec.Zc)
        return abs(rho - spec.R_edge) < 1e-6 and math.atan2(y - spec.Yc, z - spec.Zc) <= spec.th_tip + 1e-6

    def on_te_surface(pt):
        _, y, z = pt
        k = round(math.atan2(z, y) / (2 * math.pi / spec.N))
        a = -k * 2 * math.pi / spec.N
        y, z = y * math.cos(a) - z * math.sin(a), y * math.sin(a) + z * math.cos(a)
        return abs(math.hypot(y - spec.Yc, z - spec.Zc) - spec.R_edge) < 1e-4

    # Graded faces: walls, the base rim's shear layer (r = R behind the base)
    # and a fin trailing edge's wake surface (R_edge behind the TE). In a
    # block with such a face, the edges running from it to the opposite face
    # get its grading at that end, at any angle, so every edge across the
    # block's layer is graded alike.
    def face_kind(f):
        if topo.group_of.get(f) in WALLS:
            return "wall"
        bb = gmsh.model.getBoundingBox(2, f)
        if bb[0] < total - PAD:
            return None
        if _on_cylinder(f, R):
            return "shear"
        if spec and bb[0] > spec.z[3] - PAD:
            (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, f)
            pts = gmsh.model.getValue(2, f, [u0, v0, u1, v1, 0.5 * (u0 + u1), 0.5 * (v0 + v1)])
            if all(on_te_surface(pts[i:i + 3]) for i in range(0, 9, 3)):
                return "te"
        return None

    from_face: dict[tuple[int, int], str] = {}
    kinds = {f: face_kind(f) for f in topo.face_zones}
    for v in topo.info["zones"]:
        faces = [f for _, f in gmsh.model.getBoundary([(3, v)], oriented=False)]
        block_edges = {c for f in faces for c in topo.loops[f]}
        for f in faces:
            if not kinds[f]:
                continue
            corners = {p for c in topo.loops[f] for p in topo.ends[c]}
            for c in block_edges - set(topo.loops[f]):
                on = [p for p in topo.ends[c] if p in corners]
                if len(on) == 1 and (c, on[0]) not in from_face:
                    from_face[(c, on[0])] = kinds[f]

    needs = {}
    for c, cv in topo.curves.items():
        zones = topo.edge_zones[c]
        far = all(z.startswith("far") for z in zones)
        kind = topo.direction(c)
        xm = 0.5 * (cv.p[0] + cv.q[0])
        rm = max(0.5 * (cv.rp + cv.rq), 1e-3)
        slabby = bool(zones & {"slab", "cap", "slab_core"})
        nt = P["nTheta"]

        if far:
            h_ax = P["hFarMid"] if P["xZoneStart"] - PAD < xm < P["xZoneEnd"] + PAD else P["hFar"]
            h_cr = P["hFar"]
        else:
            if xm < 0:
                h_ax = P["hUpstream"]
            elif xm < body.l_ogive:
                h_ax = P["hOgive"]
            elif fin_x and fin_x[0] <= xm <= fin_x[1] and rm < r_tip:
                h_ax = P["hFinChord"]
            elif slabby and max(cv.rp, cv.rq) > R + 1e-6:
                h_ax = P["hSlabOuter"]                      # slab, away from the fins and the body
            elif xm < total:
                h_ax = P["hWall"]
            else:
                h_ax = P["hRelax"] if xm < P["xWakeSplit"] or slabby else P["hWake"]
            if slabby:
                h_cr = P["hFinSpan"] if all(in_wrap(pt) for pt in (cv.p, cv.pts[200], cv.q)) else P["hSlabOuter"]
            else:
                h_cr = P["hRingOut"]
        if kind == "axial":
            h = h_ax
        elif kind == "circ":
            h = rm * (math.pi / 4) / nt
            if slabby:
                h = min(h, h_cr)
        elif kind == "radial":
            h = h_cr if far or slabby else (R if xm > total else P["hRingOut"])
        else:
            h = h_cr

        conds = []
        for p in topo.ends[c]:
            x = topo.xyz[p][0]
            cond = None
            via = from_face.get((c, p))
            if via == "wall" or topo.normal_faces(c, p, WALLS):
                cond = tip if x < x45 + 1e-6 else wall
            elif via == "shear":
                cond = wall if x <= P["xWakeSplit"] + 1e-6 else End(P["hWakeRadial"], q)
            elif via == "te":
                cond = End(P["hLE"], q)
            elif topo.normal_faces(c, p, ("seam_",)):
                cond = End(P["hSeam"], q) if far else (End(P["hRingOut"], q) if abs(
                    math.hypot(*topo.xyz[p][1:]) - r_z) < PAD else None)
            elif kind == "axial" and abs(x - body.l_ogive) < 1e-6 and not far:
                cond = End(P["hShoulder"], q)
            elif kind == "axial" and spec and on_fin_edge(p, c):
                cond = End(P["hLE"], q)                     # chordwise at the LE/TE wedges
            elif kind != "axial" and x > total + 1e-6 and abs(math.hypot(*topo.xyz[p][1:]) - R) < 1e-6:
                # the base rim's shear layer, carried downstream on r = R
                cond = wall if x <= P["xWakeSplit"] + 1e-6 else End(P["hWakeRadial"], q)
            elif kind != "axial" and spec and x > total + 1e-6 and on_te_line(p):
                cond = End(P["hLE"], q)                     # the fin trailing edge's wake
            conds.append(cond)
        needs[c] = (conds[0], conds[1], h)

    # Face-consistent wall stacks: if one edge running from a block face has
    # the wall stack at that face (e.g. the fin tip reached along a fin
    # wall), the other edges from that face in the block take it too. Only
    # across faces that border a wall (an edge in a wall face), so the stacks
    # do not spread out from wall corners into the free flow.
    wall_edges = {c for f, k in kinds.items() if k == "wall" for c in topo.loops[f]}
    connecting: list[list[tuple[int, int]]] = []
    for v in topo.info["zones"]:
        faces = [f for _, f in gmsh.model.getBoundary([(3, v)], oriented=False)]
        block_edges = {c for f in faces for c in topo.loops[f]}
        for f in faces:
            if not set(topo.loops[f]) & wall_edges:
                continue
            corners = {p for c in topo.loops[f] for p in topo.ends[c]}
            group = []
            for c in block_edges - set(topo.loops[f]):
                on = [p for p in topo.ends[c] if p in corners]
                if len(on) == 1:
                    group.append((c, on[0]))
            connecting.append(group)

    def at(c, p):
        return needs[c][0] if p == topo.ends[c][0] else needs[c][1]

    def put(c, p, cond):
        a, b, _ = needs[c]
        needs[c] = (cond, b, h) if p == topo.ends[c][0] else (a, cond, h)

    changed = True
    while changed:
        changed = False
        for group in connecting:
            stacked = [at(c, p) for c, p in group if at(c, p) is not None and at(c, p).n_layers]
            if not stacked:
                continue
            best = min(stacked, key=lambda e: e.h)
            for c, p in group:
                cur = at(c, p)
                if cur is None or not cur.n_layers:
                    put(c, p, best)
                    changed = True
    return needs


def family_counts(topo: Topo, fam: Families, needs: dict) -> dict[int, int]:
    """Cells per family: the most any member needs (wall stacks and seam sizes
    at its ends, grown at `growth` up to its size cap)."""
    counts = {}
    for root, members in fam.groups().items():
        n = 1
        for c in members:
            a, b, h = needs[c]
            L = topo.curves[c].length
            a = a or End(h)
            b = b or End(h)
            n = max(n, len(size_driven(L, a, b, h)))
        counts[root] = n
    return counts


def block_cells(topo: Topo, fam: Families, counts: dict) -> dict[str, int]:
    """Cells per zone: per block, the product of the counts on three edges
    meeting at one corner (a prism block counts as half its hexahedron)."""
    cells: dict[str, int] = {}
    for v, name in topo.info["zones"].items():
        faces = [f for _, f in gmsh.model.getBoundary([(3, v)], oriented=False)]
        edges = {c for f in faces for c in topo.loops[f]}
        corner_edges: dict[int, list[int]] = {}
        for c in edges:
            for p in topo.ends[c]:
                corner_edges.setdefault(p, []).append(c)
        three = next(es for es in corner_edges.values() if len(es) == 3)
        n = math.prod(counts[fam.find(c)[0]] for c in three)
        if len(faces) == 5:
            n //= 2
        cells[name] = cells.get(name, 0) + n
    return cells


def stronger(a, b):
    """The end condition with the smaller first cell (a wall stack wins ties)."""
    if a is None or b is None:
        return a or b
    if abs(a.h - b.h) > 1e-12:
        return a if a.h < b.h else b
    return a if a.n_layers >= b.n_layers else b


def two_sided(L: float, n: int, a, b) -> tuple[list[float], bool]:
    """n cell sizes over L from the start condition a to the end condition b
    (End or None for free). Returns the sizes and whether both conditions
    were met (False: too few cells, wall stacks were merged)."""
    if a is None and b is None:
        return [L / n] * n, True
    if b is None:
        return fixed_n(L, n, a), True
    if a is None:
        return fixed_n(L, n, None, b), True
    lo, hi = min(a.h, b.h, L), L                   # a cap below the first cells means nothing
    if len(size_driven(L, a, b, lo)) < n:          # very short edge: split its largest cells
        sz = size_driven(L, a, b, lo)
        while len(sz) < n:
            i = max(range(len(sz)), key=sz.__getitem__)
            sz[i:i + 1] = [sz[i] / 2] * 2
        return sz, True
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if len(size_driven(L, a, b, mid)) > n:
            lo = mid
        else:
            hi = mid
    sz = size_driven(L, a, b, hi)
    if len(sz) == n:
        return sz, True
    sz = size_driven(L, a, b, lo)                   # one or two cells too many: merge in the middle
    ok = True
    while len(sz) > n:
        la, lb = min(a.n_layers, len(sz) // 2), min(b.n_layers, len(sz) // 2)
        rng = range(la, len(sz) - lb - 1)
        if not rng:
            rng, ok = range(len(sz) - 1), False
        i = min(rng, key=lambda k: sz[k] + sz[k + 1])
        sz[i:i + 2] = [sz[i] + sz[i + 1]]
    return sz, ok


def seam_stations(topo: Topo) -> list[dict]:
    """The circumferential arcs on each cylindrical seam, grouped by station
    (x, r) and seam side, with their phi ranges."""
    by: dict[tuple, dict] = {}
    for f, g in ((f, topo.group_of.get(f, "")) for f in topo.face_zones):
        if not g.startswith(("seam_side", "seam_core_side")):
            continue
        pair, role = g[5:].rsplit("_", 1)
        for c in topo.loops[f]:
            if topo.direction(c) != "circ":
                continue
            cv = topo.curves[c]
            key = (pair, round(cv.p[0], 6), round(cv.rp, 6))
            p0 = math.atan2(cv.p[2], cv.p[1])
            p1 = p0 + math.atan2(math.sin(math.atan2(cv.q[2], cv.q[1]) - p0),
                                 math.cos(math.atan2(cv.q[2], cv.q[1]) - p0))
            by.setdefault(key, {"master": {}, "slave": {}})[role][c] = (p0, p1)
    return [{"key": k, **v} for k, v in by.items()]


def generic_spacing(body: Body, P: dict, info: dict, spec) -> dict:
    """Cell sizes on every block edge (F5b): family counts, per-edge
    distributions from the end conditions (axial families share theirs),
    and node positions copied across the cylindrical seams."""
    topo = Topo(info)
    fam = edge_families(topo)
    needs = edge_needs(topo, body, P, spec)
    groups = fam.groups()
    ends: dict[int, tuple] = {}
    for root, members in groups.items():
        if all(topo.direction(c) == "axial" for c in members):
            A = B = None
            for c in members:
                a, b, _ = needs[c]
                s = fam.find(c)[1]
                ra, rb = (a, b) if s == 1 else (b, a)
                A, B = stronger(A, ra), stronger(B, rb)
            for c in members:
                ends[c] = (A, B) if fam.find(c)[1] == 1 else (B, A)
        else:
            for c in members:
                ends[c] = needs[c][:2]
    counts = {}
    for root, members in groups.items():
        n = 1
        for c in members:
            a, b = ends[c]
            h = needs[c][2]
            n = max(n, len(size_driven(topo.curves[c].length, a or End(h), b or End(h), h)))
        counts[root] = n

    # Cylindrical seams whose arcs do not coincide: the side with more arcs
    # at a station is the source; a target arc spanning several source arcs
    # takes their cell count and node positions, and a target corner inside
    # a source arc forces a node there.
    targets: dict[int, list[tuple]] = {}       # target arc -> [(source arc, phi range)] in phi order
    forced: dict[int, list[float]] = {}        # source arc -> phis that must be nodes
    for st in seam_stations(topo):
        sides = sorted((st["master"], st["slave"]), key=len)
        tgt, src = sides
        if not tgt or not src:
            continue
        for c, (t0, t1) in tgt.items():
            lo, hi = min(t0, t1), max(t0, t1)
            inside = []
            for e, (s0, s1) in src.items():
                a0, a1 = min(s0, s1), max(s0, s1)
                for shift in (0.0, 2 * math.pi, -2 * math.pi):
                    b0, b1 = a0 + shift, a1 + shift
                    if b1 > lo + 1e-9 and b0 < hi - 1e-9:
                        inside.append((e, b0, b1))
                        for phi in (lo, hi):
                            if b0 + 1e-9 < phi < b1 - 1e-9:
                                forced.setdefault(e, []).append(phi - shift)
                        break
            if len(inside) == 1 and abs(inside[0][1] - lo) < 1e-9 and abs(inside[0][2] - hi) < 1e-9:
                continue                                    # coincident: same family already
            targets[c] = (lo, hi, sorted(inside, key=lambda t: t[1]))

    def dist(c, n):
        """A free end takes the edge's size cap, so growth from a stack or a
        seam size stops at the cap instead of running on geometrically."""
        a, b = ends[c]
        h = needs[c][2]
        if a is None and b is None:
            return two_sided(topo.curves[c].length, n, None, None)
        return two_sided(topo.curves[c].length, n, a or End(h), b or End(h))

    sizes, short = {}, []
    for root, members in groups.items():
        for c in members:
            sizes[c], ok = dist(c, counts[root])
            if not ok:
                short.append(c)

    def phis(c):
        """Node angles along arc c, in its parametric order, from its sizes."""
        cv = topo.curves[c]
        p0 = math.atan2(cv.p[2], cv.p[1])
        p1 = p0 + math.atan2(math.sin(math.atan2(cv.q[2], cv.q[1]) - p0), math.cos(math.atan2(cv.q[2], cv.q[1]) - p0))
        out, acc = [p0], 0.0
        for x in sizes[c]:
            acc += x
            out.append(p0 + (p1 - p0) * acc / cv.length)
        return out

    def set_phis(c, ph):
        cv = topo.curves[c]
        sz = [abs(ph[i + 1] - ph[i]) for i in range(len(ph) - 1)]
        tot = sum(sz)
        sizes[c] = [x * cv.length / tot for x in sz]

    for e, fs in forced.items():                        # snap the nearest nodes onto forced corners
        ph = phis(e)
        for phi in fs:
            i = min(range(1, len(ph) - 1), key=lambda k: abs(ph[k] - phi))
            old = ph[i]
            for k in range(1, len(ph) - 1):             # stretch each side linearly
                if k <= i:
                    ph[k] = ph[0] + (ph[k] - ph[0]) * (phi - ph[0]) / (old - ph[0])
                else:
                    ph[k] = phi + (ph[k] - old) * (ph[-1] - phi) / (ph[-1] - old)
            ph[i] = phi
        set_phis(e, ph)
    for c, (lo, hi, chain) in targets.items():
        cv = topo.curves[c]
        p0 = math.atan2(cv.p[2], cv.p[1])
        nodes = set()                                   # source nodes within the target's own range
        for e, b0, b1 in chain:
            shift = b0 - min(phis(e)[0], phis(e)[-1])
            for ph in phis(e):
                x = ph + shift
                if lo - 1e-9 <= x <= hi + 1e-9:
                    nodes.add(round(x, 12))
        ph = sorted(nodes)
        def gap(a, b):
            return abs(math.atan2(math.sin(a - b), math.cos(a - b)))
        if gap(ph[0], p0) > gap(ph[-1], p0):            # the target runs from its high end
            ph = ph[::-1]
        root = fam.find(c)[0]
        n = len(ph) - 1
        if counts[root] != n:                           # the target's whole family takes the count
            counts[root] = n
            for m in groups[root]:
                sizes[m], _ = dist(m, n)
        set_phis(c, ph)
    return {"topo": topo, "fam": fam, "counts": counts, "sizes": sizes, "ends": ends,
            "targets": targets, "forced": forced, "short": short}


def generic_mesh(body: Body, P: dict, info: dict, spec, out: Path, convert: bool) -> list[str]:
    """F5b/F6: mesh every block edge with its generic spacing, then the faces
    and blocks, and write the volume mesh (see write_case)."""
    gs = generic_spacing(body, P, info, spec)
    topo, sizes, counts, fam = gs["topo"], gs["sizes"], gs["counts"], gs["fam"]
    n_forced = sum(len(v) for v in gs["forced"].values())
    lines = [(f"  {len(topo.curves)} block edges, {len(counts)} families, {len(gs['targets'])} seam arcs copied, "
              f"{n_forced} forced nodes, {len(gs['short'])} edges too short for both wall stacks")]
    bad = 0
    for f in topo.loops:
        _, edges, _ = topo.cycle(f)
        if len(edges) == 4 and (len(sizes[edges[0]]) != len(sizes[edges[2]])
                                or len(sizes[edges[1]]) != len(sizes[edges[3]])):
            bad += 1
    lines.append(f"  faces with unequal opposite counts: {bad}  {'ok' if not bad else 'FAIL'}")
    cells = block_cells(topo, fam, counts)
    lines.append(f"  estimated cells: {sum(cells.values()) / 1e6:.3f} M")
    wall_first = [sz[0] if e[0] and e[0].n_layers else sz[-1]
                  for c, sz in sizes.items() for e in [gs["ends"][c]] if (e[0] and e[0].n_layers) or (e[1] and e[1].n_layers)]
    if wall_first:
        lines.append(f"  first wall cell: min {min(wall_first) * 1e3:.2f} um, max {max(wall_first) * 1e3:.2f} um "
                     f"(target {P['firstLayer'] * 1e3:.1f} um)")
    for c, sz in sizes.items():
        gmsh.model.mesh.setTransfiniteCurve(c, len(sz) + 1)
    set_transfinite_blocks(info)
    gmsh.model.mesh.generate(1)
    for c, sz in sizes.items():                         # move every curve's nodes to its sizes
        cv = topo.curves[c]
        tags, _, u = gmsh.model.mesh.getNodes(1, c, includeBoundary=False, returnParametricCoord=True)
        order = sorted(range(len(tags)), key=lambda i: u[i])
        if cv.u[-1] < cv.u[0]:
            order = order[::-1]
        s, total = 0.0, sum(sz)
        for k, i in enumerate(order):
            s += sz[k]
            uu = cv.u_at(s * cv.length / total)
            gmsh.model.mesh.setNode(tags[i], gmsh.model.getValue(1, c, [uu]), [uu])
    # Coincident seam edges: the slave copies the master's node coordinates.
    # They are separate curves, and arc-length placement on two B-spline
    # parametrisations differs by ~1e-6 mm, more than stitchMesh -perfect allows.
    for cm, cs, sign in fam.links:
        mt, mc, mu = gmsh.model.mesh.getNodes(1, cm, includeBoundary=False, returnParametricCoord=True)
        st, _, su = gmsh.model.mesh.getNodes(1, cs, includeBoundary=False, returnParametricCoord=True)
        if len(mt) != len(st):
            continue
        mo = sorted(range(len(mt)), key=lambda i: mu[i])
        so = sorted(range(len(st)), key=lambda i: su[i])
        cvm, cvs = topo.curves[cm], topo.curves[cs]
        if cvm.u[-1] < cvm.u[0]:
            mo = mo[::-1]
        if cvs.u[-1] < cvs.u[0]:
            so = so[::-1]
        if sign < 0:
            so = so[::-1]
        for i, j in zip(mo, so):
            xyz = list(mc[3 * i:3 * i + 3])
            gmsh.model.mesh.setNode(st[j], xyz, gmsh.model.getParametrization(1, cs, xyz))
    # seam check: the copied arcs' nodes against their sources'
    worst = 0.0
    for c, (_, _, chain) in gs["targets"].items():
        tn = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)[1]
        tp = [tn[3 * i:3 * i + 3] for i in range(len(tn) // 3)]
        sp = []
        for e, _, _ in chain:
            en = gmsh.model.mesh.getNodes(1, e, includeBoundary=True)[1]
            sp += [en[3 * i:3 * i + 3] for i in range(len(en) // 3)]
        for p in tp:
            worst = max(worst, min(math.dist(p, q) for q in sp))
    lines.append(f"  copied seam nodes: max distance to a source node {worst:.2e} mm  "
                 f"{'ok' if worst < 1e-6 else 'FAIL'}")
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.model.mesh.generate(3)
    lines.append("  " + write_case(info, out, convert))
    return lines


def budget(body: Body, P: dict, info: dict, spec) -> list[str]:
    """F5a: edge families, their counts and the cells per zone, no meshing."""
    gs = generic_spacing(body, P, info, spec)
    topo, fam, counts = gs["topo"], gs["fam"], gs["counts"]
    needs = edge_needs(topo, body, P, spec)
    cells = block_cells(topo, fam, counts)
    lines = [f"  {len(topo.curves)} block edges in {len(counts)} families; {len(fam.conflicts)} orientation conflicts"]
    for a, b in fam.conflicts[:5]:
        lines.append(f"    conflict: edges {a} {b}")
    total = sum(cells.values())
    for name, n in sorted(cells.items(), key=lambda kv: -kv[1]):
        lines.append(f"  zone {name:9s} {n / 1e6:7.3f} M cells")
    lines.append(f"  total          {total / 1e6:7.3f} M cells")
    top = sorted(counts.items(), key=lambda kv: -kv[1])[:30]
    lines.append("  largest families (count, members, an edge: kind, zones, from -> to):")
    groups = fam.groups()
    for root, n in top:
        c = max(groups[root], key=lambda e: topo.curves[e].length)
        cv = topo.curves[c]
        a, b, _ = needs[c]
        lab = "".join("W" if e is not None and e.n_layers else ("h" if e is not None else "-") for e in (a, b))
        lines.append(f"    {n:5d} x{len(groups[root]):4d} {lab} {topo.direction(c):7s} "
                     f"{','.join(sorted(topo.edge_zones[c]))[:28]:28s} "
                     f"({cv.p[0]:.1f}, r {cv.rp:.1f}) -> ({cv.q[0]:.1f}, r {cv.rq:.1f})")
    return lines


# ── fins (benchmark configuration, report F1) ──────────────────────────────

class FinSpec:
    """Arc fin exactly as geometry/arc_stabilizers.scad builds it (report F1).

    Radii are absolute (Table 4: R3/R1/R2); the scad asserts R_out == R.
    In OpenFOAM axes the cross-section of fin 0 is (Y, Z) = (Yc + r sin a,
    Zc + r cos a) about the arc centre (Yc, Zc); the scad's mirror makes the
    bow bulge towards +Z, and fin k is fin 0 rotated by k * 360 / N about +X.
    """

    R_in, R_edge, R_out = 36.0, 38.0, 40.0
    root_embed, delta = 2.0, 45.0

    def __init__(self, body: Body, N: int, xi_deg: float, L: float):
        if abs(self.R_out - body.R) > 1e-9:
            raise ValueError("R_out must equal the body radius (scad assert)")
        self.N, self.L = N, L
        self.xi = math.radians(xi_deg)
        self.root_y = body.R - self.root_embed
        self.Yc = self.root_y + self.R_edge * math.sin(self.xi / 2)
        self.Zc = -self.R_edge * math.cos(self.xi / 2)
        self.chamfer = (self.R_out - self.R_in) / 2 / math.tan(math.radians(self.delta) / 2)
        self.z = [body.total - L, body.total - L + self.chamfer, body.total - self.chamfer, body.total]
        self.th = {r: math.atan2(self.root_y - self.Yc, self.Zc_at(r) - self.Zc)
                   for r in (self.R_in, self.R_edge, self.R_out)}
        self.th_tip = self.th[self.R_edge] + self.xi

    def Zc_at(self, r: float) -> float:
        """Z where the circle of radius r about the centre meets Y = root_y (root side)."""
        return self.Zc + math.sqrt(r * r - (self.root_y - self.Yc) ** 2)

    def point(self, r: float, a: float, x: float) -> tuple[float, float, float]:
        return (x, self.Yc + r * math.sin(a), self.Zc + r * math.cos(a))

    def cap(self, delta: float) -> dict[str, tuple[float, float]]:
        """Fin 0's tip cap in (Y, Z): the tip-face corners bi/bo, the wrap
        corners Pm/Pp a wrap `delta` outside the faces and past the tip, and
        Q on the mid arc between them. All four sides are straight: the tip
        face and the wrap offset lie on radial lines of the arc centre."""
        yz = lambda r, a: self.point(r, a, 0.0)[1:]
        t, dA = self.th_tip, delta / self.R_edge
        return {"bi": yz(self.R_in, t), "bo": yz(self.R_out, t), "e": yz(self.R_edge, t),
                "Pm": yz(self.R_in - delta, t + dA), "Pp": yz(self.R_out + delta, t + dA),
                "Q": yz(self.R_edge, t + dA)}

    def wrap_root(self, delta: float, R: float) -> float:
        """phi (rad) where fin 0's outer wrap arc (R_out + delta) meets r = R."""
        rho = self.R_out + delta
        lo, hi = -math.pi / 2, self.th_tip
        for _ in range(200):
            m = 0.5 * (lo + hi)
            if math.hypot(*self.point(rho, m, 0.0)[1:]) < R:
                lo = m
            else:
                hi = m
        y, z = self.point(rho, lo, 0.0)[1:]
        return math.atan2(z, y)

    def tip_radius(self) -> float:
        return math.hypot(*self.point(self.R_edge, self.th_tip, 0.0)[1:])


def build_fins(spec: FinSpec) -> list[tuple[int, int]]:
    """The N fin solids. Each is sewn from the scad's faces: ruled LE/TE wedges
    between the edge arc and the inner/outer arcs (proportional parametrisation,
    as the scad's index-matched polyhedron), cylindrical barrels, and planar root
    and tip caps."""
    occ = gmsh.model.occ
    z0, z1, z2, z3 = spec.z
    Ri, Re, Ro = spec.R_in, spec.R_edge, spec.R_out

    def arc(r, x):
        c = occ.addPoint(x, spec.Yc, spec.Zc)
        p0 = occ.addPoint(*spec.point(r, spec.th[r], x))
        p1 = occ.addPoint(*spec.point(r, spec.th_tip, x))
        a = occ.addCircleArc(p0, c, p1)
        occ.remove([(0, c)])
        return occ.addWire([a])

    faces = []
    for (ra, xa), (rb, xb) in (((Re, z0), (Ro, z1)), ((Re, z0), (Ri, z1)), ((Ro, z1), (Ro, z2)),
                               ((Ri, z1), (Ri, z2)), ((Ro, z2), (Re, z3)), ((Ri, z2), (Re, z3))):
        faces += [t for d, t in occ.addThruSections([arc(ra, xa), arc(rb, xb)], makeSolid=False,
                                                    makeRuled=True) if d == 2]
    for a_of in (lambda r: spec.th[r], lambda r: spec.th_tip):        # root cap, tip cap
        ring = [(Re, z0), (Ro, z1), (Ro, z2), (Re, z3), (Ri, z2), (Ri, z1)]
        pts = [occ.addPoint(*spec.point(r, a_of(r), x)) for r, x in ring]
        lines = [occ.addLine(pts[i], pts[(i + 1) % 6]) for i in range(6)]
        faces.append(occ.addPlaneSurface([occ.addCurveLoop(lines)]))
    fin = (3, occ.addVolume([occ.addSurfaceLoop(faces, sewing=True)]))
    out = [fin]
    for k in range(1, spec.N):
        copy = occ.copy([fin])
        occ.rotate(copy, 0, 0, 0, 1, 0, 0, k * 2 * math.pi / spec.N)
        out += copy
    occ.synchronize()
    return out


def stl_mass(path: Path) -> tuple[float, float, float]:
    """Volume, area and max radius about the x axis of an ASCII or binary STL."""
    data = path.read_bytes()
    tris = []
    if data[:5] == b"solid" and b"facet" in data[:400]:
        v = [tuple(map(float, ln.split()[1:4])) for ln in data.decode().splitlines()
             if ln.strip().startswith("vertex")]
        tris = [v[i:i + 3] for i in range(0, len(v), 3)]
    else:
        n = struct.unpack("<I", data[80:84])[0]
        for i in range(n):
            f = struct.unpack("<12f", data[84 + 50 * i:84 + 50 * i + 48])
            tris.append([f[3:6], f[6:9], f[9:12]])
    vol = area = rmax = 0.0
    for a, b, c in tris:
        cr = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        vol += (cr[0] * c[0] + cr[1] * c[1] + cr[2] * c[2]) / 6
        u = [b[i] - a[i] for i in range(3)]
        w = [c[i] - a[i] for i in range(3)]
        n_ = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0])
        area += math.sqrt(sum(x * x for x in n_)) / 2
        rmax = max(rmax, *(math.hypot(p[1], p[2]) for p in (a, b, c)))
    return abs(vol), area, rmax


def fin_check(spec: FinSpec, fins: list, case: Path, props: dict) -> list[str]:
    """Compare the OCC fins with the scad's own stabilizers.stl and the report."""
    occ = gmsh.model.occ
    vol = sum(occ.getMass(*f) for f in fins)
    area = sum(occ.getMass(*s) for f in fins for s in gmsh.model.getBoundary([f], oriented=False))
    lines = [(f"  fins: N = {spec.N}, xi = {math.degrees(spec.xi):.2f} deg, chord {spec.L:g} mm "
              f"(x {spec.z[0]:g}-{spec.z[3]:g}), wedge run {spec.chamfer:.3f} mm"),
             (f"  tip on the root's radial line at r = {spec.tip_radius():.3f} mm "
              f"(scad: R - root_embed + 2 R_edge sin(xi/2))"),
             f"  arc centre (Y, Z) = ({spec.Yc:.4f}, {spec.Zc:.4f}) mm; root angles inner/edge/outer "
             + "/".join(f"{math.degrees(spec.th[r]):.2f}" for r in (spec.R_in, spec.R_edge, spec.R_out))
             + " deg"]
    scad = Path(__file__).resolve().parents[1] / "geometry" / "arc_stabilizers.scad"
    stl = case / "fins_scad.stl"
    openscad = next((e for e in ("openscad-nightly", "openscad") if shutil.which(e)), None)
    if openscad:
        defs = f"D={props['D']:g}; N={spec.N}; xi={math.degrees(spec.xi):.6g}; L={spec.L:g}; EXPORT=\"stabilizers\";"
        subprocess.run([openscad, "-o", str(stl), "-D", defs, str(scad)], check=True, capture_output=True)
        sv, sa, sr = stl_mass(stl)
        # farthest fin point from the axis: the tip's outer corner, on both sides
        r_occ = max(math.hypot(*gmsh.model.getValue(0, p, [])[1:])
                    for f in fins for _, p in gmsh.model.getBoundary([f], recursive=True))
        ok = True
        for label, got, want, tol in (("fin volume (mm^3)", vol, sv, 1e-3), ("fin area (mm^2)", area, sa, 1e-3),
                                      ("max radius (mm)", r_occ, sr, 1e-4)):
            err = abs(got - want) / want
            ok &= err < tol
            lines.append(f"  {label:20s} OCC {got:14.4f}  scad STL {want:14.4f}  rel.err {err:.1e} "
                         f"{'ok' if err < tol else 'FAIL'}")
        lines.append("  " + ("FIN CHECKS PASSED" if ok else "FIN CHECKS FAILED"))
    else:
        lines.append(f"  OCC fin volume {vol:.4f} mm^3, area {area:.4f} mm^2 (no OpenSCAD to compare)")
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
    fins = info["fins"]
    a_wall = sum(occ.getMass(2, s) for k in ("fuselage", "stabilizers") for s in info["groups"].get(k, []))
    if fins:
        row("body + fin wall area (mm^2)", a_wall, fins["wall_area"])
    else:
        row("body wall area (mm^2)", a_wall, a_lat + math.pi * R * R)
    v_dom = math.pi * P["rFar"] ** 2 * (P["xOutlet"] - P["xInlet"])
    v_fluid = sum(occ.getMass(3, v) for v in info["zones"])
    v_fins = fins.get("exposed_volume", 0.0)
    row("fluid volume (mm^3)", v_fluid, v_dom - v_body - v_fins)
    for name, x0, x1, inner, _ in info["layout"]:
        got = sum(occ.getMass(3, v) for v, n in info["zones"].items() if n == name)
        if inner:
            lo, hi = max(0.0, x0), min(total, x1)
            vb = body.integrate(lo, hi)[0] if hi > lo else 0.0
            want = math.pi * r_z**2 * (x1 - x0) - vb
            if name == "slab":
                want -= v_fins + math.pi * R**2 * (x1 - total) + info["cap_volume"]
            elif name == "cap":
                want = info["cap_volume"]
            elif name == "slab_core":
                want = math.pi * R**2 * (x1 - x0)
        else:
            overlap = max(0.0, min(x1, xs[3]) - max(x0, xs[0]))
            want = math.pi * (P["rFar"] ** 2 * (x1 - x0) - r_z**2 * overlap)
        row(f"zone {name} volume (mm^3)", got, want)

    # Block topology: every block must be a hexahedron (6 faces of 4 edges,
    # 12 edges, 8 corners) so it can carry a structured mesh. The fin slab
    # is not split into blocks yet.
    if info["slab"]:
        n = info["fins_n"]
        expected = {"nose": 20, "body": 16, "slab": 52 * n, "slab_core": 5, "cap": 6 * n, "wake": 20,
                    "far": 28, "far_slab": 32, "far_aft": 28}
    else:
        expected = {"nose": 20, "body": 16, "aft": 48, "far": 48}
    for name, want in expected.items():
        blocks = [v for v, n in info["zones"].items() if n == name]
        if want is None:
            lines.append(f"  zone {name:8s}: {len(blocks):3d} volume(s), block split pending")
            continue
        bad, prisms = [], 0
        for v in blocks:
            faces = gmsh.model.getBoundary([(3, v)], oriented=False)
            edges = {c for f in faces for _, c in real_edges(f)}
            corners = {p for c in edges for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)}
            per_face = [len(real_edges(f)) for f in faces]
            hexa = len(faces) == 6 and len(edges) == 12 and len(corners) == 8 and set(per_face) == {4}
            prism = (len(faces) == 5 and len(edges) == 9 and len(corners) == 6
                     and sorted(per_face) == [3, 3, 4, 4, 4])
            if prism:
                prisms += 1
            if not (hexa or prism):
                bad.append(f"{v}: {len(faces)}f/{len(edges)}e/{len(corners)}c {per_face}")
        good = len(blocks) == want and not bad
        ok &= good
        lines.append(f"  zone {name:9s}: {len(blocks):3d} blocks (expected {want}), "
                     f"{len(blocks) - len(bad) - prisms} hexahedral, {prisms} prisms  {'ok' if good else 'FAIL'}")
        for b in bad[:6]:
            lines.append(f"      not hex: block {b}")

    # Block shapes: the corner angles of every planar block face, measured
    # inside the face. Topology alone passes a block with a reflex corner,
    # whose cells fold.
    worst_corner = {}
    for v, name in info["zones"].items():
        for f in gmsh.model.getBoundary([(3, v)], oriented=False):
            if gmsh.model.getType(*f) != "Plane":
                continue
            for ang, where in face_corner_angles(f[1]):
                if ang > worst_corner.get(name, (0.0, None))[0]:
                    worst_corner[name] = (ang, where)
    for name, (ang, where) in sorted(worst_corner.items(), key=lambda kv: -kv[1][0]):
        good = ang < 178.0
        ok &= good
        lines.append(f"  zone {name:9s}: largest block-face corner {ang:6.1f} deg at "
                     f"({where[0]:.1f}, {where[1]:.1f}, {where[2]:.1f})  {'ok' if good else 'FAIL'}")

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


def face_corner_angles(f: int) -> list[tuple[float, tuple]]:
    """(angle in degrees, corner point) at each corner of planar face f,
    measured inside the face from the edge tangents (over 180 = reflex)."""
    edges = [c for _, c in real_edges((2, f))]
    ends = {c: [p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)] for c in edges}
    out = []
    for p in {p for c in edges for p in ends[c]}:
        cs = [c for c in edges if p in ends[c]]
        if len(cs) != 2:
            continue
        P0 = gmsh.model.getValue(0, p, [])
        dirs = []
        for c in cs:
            (t0,), (t1,) = gmsh.model.getParametrizationBounds(1, c)
            a, b = gmsh.model.getValue(1, c, [t0]), gmsh.model.getValue(1, c, [t1])
            start = math.dist(a, P0) < math.dist(b, P0)
            t = t0 + (t1 - t0) * (1e-4 if start else 1 - 1e-4)
            q = gmsh.model.getValue(1, c, [t])
            d = [q[i] - P0[i] for i in range(3)]
            n = math.sqrt(sum(x * x for x in d))
            dirs.append([x / n for x in d])
        cosang = max(-1.0, min(1.0, sum(dirs[0][i] * dirs[1][i] for i in range(3))))
        ang = math.degrees(math.acos(cosang))
        bis = [dirs[0][i] + dirs[1][i] for i in range(3)]
        nb = math.sqrt(sum(x * x for x in bis))
        if nb > 1e-9:
            probe = [P0[i] + 0.02 * bis[i] / nb for i in range(3)]
            if not gmsh.model.isInside(2, f, probe):
                ang = 360.0 - ang
        out.append((ang, tuple(P0)))
    return out


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
        if gmsh.model.getPhysicalName(dim, tag) in info["shared"]:
            gmsh.model.removePhysicalGroups([(dim, tag)])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(str(out))
    return {"counts": counts, "worst": worst, "low": low}


PREVIEW_CONTROLDICT = """FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application none; startFrom startTime; startTime 0; stopAt endTime; endTime 0; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 8; timeFormat general;
"""


def prism_corners(v: int) -> list[int]:
    """For a prism block: its 6 corners (one triangle, then the other, each
    corner above its partner), with both triangles' transfinite corners set
    to match, so their collapsed corners lie on the same lateral edge. [] for
    any other block."""
    faces = gmsh.model.getBoundary([(3, v)], oriented=False)
    tris = [f for f in faces if len(real_edges(f)) == 3]
    if len(faces) != 5 or len(tris) != 2:
        return []
    pts = [[abs(p) for _, p in gmsh.model.getBoundary([f], combined=True, oriented=False, recursive=True)]
           for f in tris]
    # Gmsh collapses a transfinite triangle at its first corner, fanning one
    # row of prisms around it: use the widest corner. A fan at the sharp
    # (22.5-degree) trailing-edge corner gives sliver prisms whose centres
    # fall outside their faces.
    xyz = {p: gmsh.model.getValue(0, p, []) for p in pts[0]}

    def angle(p):
        a, b = [q for q in pts[0] if q != p]
        u = [xyz[a][i] - xyz[p][i] for i in range(3)]
        w = [xyz[b][i] - xyz[p][i] for i in range(3)]
        return math.acos(sum(u[i] * w[i] for i in range(3)) / math.dist(xyz[a], xyz[p]) / math.dist(xyz[b], xyz[p]))

    first = max(pts[0], key=angle)
    i = pts[0].index(first)
    bottom, top_set = pts[0][i:] + pts[0][:i], set(pts[1])
    # right-handed: the bottom triangle's normal points to the top one
    top_c = [sum(gmsh.model.getValue(0, p, [])[k] for p in pts[1]) / 3 for k in range(3)]
    a, b, c = (xyz[p] for p in bottom)
    u, w = [b[k] - a[k] for k in range(3)], [c[k] - a[k] for k in range(3)]
    nrm = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    if sum(nrm[k] * (top_c[k] - a[k]) for k in range(3)) < 0:
        bottom = [bottom[0], bottom[2], bottom[1]]
    ends = []
    for f in faces:
        for _, c in real_edges(f):
            ends.append([abs(p) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)])
    top = [next(b if a == p else a for a, b in ends if p in (a, b) and ({a, b} - {p}) <= top_set) for p in bottom]
    gmsh.model.mesh.setTransfiniteSurface(tris[0][1], cornerTags=bottom)
    gmsh.model.mesh.setTransfiniteSurface(tris[1][1], cornerTags=top)
    return bottom + top


def set_transfinite_blocks(info: dict):
    """Transfinite quads on every face and hexahedra/prisms in every block
    (the curves must be set already)."""
    for _, f in gmsh.model.getEntities(2):
        loop = [c for _, c in real_edges((2, f))]
        corners = []
        if len(loop) != len(gmsh.model.getBoundary([(2, f)], oriented=False)):
            # a face at the apex: its zero-length pole edge hides a corner
            ends = {c: [abs(p) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)] for c in loop}
            here, todo = ends[loop[0]][0], list(loop)
            while todo:
                c = next(c for c in todo if here in ends[c])
                todo.remove(c)
                corners.append(here)
                here = ends[c][1] if ends[c][0] == here else ends[c][0]
        gmsh.model.mesh.setTransfiniteSurface(f, cornerTags=corners)
        gmsh.model.mesh.setRecombine(2, f)
    for v in info["zones"]:
        corners = prism_corners(v) or hex_corners(v)
        gmsh.model.mesh.setTransfiniteVolume(v, cornerTags=corners)
        gmsh.model.mesh.setRecombine(3, v)


def write_case(info: dict, out: Path, convert: bool = True) -> str:
    """Write the volume mesh as out/mesh.msh (MSH 2.2, one cell zone per mesh
    zone, shared faces dropped) and, with convert and gmshToFoam on PATH,
    turn out/ into an OpenFOAM case (case.foam for ParaView, every seam its
    own unstitched patch) and run checkMesh there."""
    counts = {}
    for v in info["zones"]:
        for et, tags, _ in zip(*gmsh.model.mesh.getElements(3, v)):
            name = gmsh.model.mesh.getElementProperties(et)[0]
            counts[name] = counts.get(name, 0) + len(tags)
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        if gmsh.model.getPhysicalName(dim, tag) in info["shared"]:
            gmsh.model.removePhysicalGroups([(dim, tag)])
    out.mkdir(parents=True, exist_ok=True)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(str(out / "mesh.msh"))
    summary = f"{sum(counts.values())} cells {counts}, wrote {out / 'mesh.msh'}"
    if not convert:
        return summary
    if not shutil.which("gmshToFoam"):
        return summary + "; gmshToFoam not on PATH, so no OpenFOAM case"
    (out / "system").mkdir(exist_ok=True)
    (out / "system" / "controlDict").write_text(PREVIEW_CONTROLDICT)
    for name in ("fvSchemes", "fvSolution"):                   # checkMesh reads them
        shutil.copy(out.parent / "system" / name, out / "system" / name)
    (out / "case.foam").touch()
    shutil.rmtree(out / "constant" / "polyMesh", ignore_errors=True)
    for cmd in (["gmshToFoam", "mesh.msh"], ["checkMesh", "-constant", "-noZero"]):
        log = out / f"log.{cmd[0]}"
        with log.open("w") as fh:
            subprocess.run(cmd, cwd=out, stdout=fh, stderr=subprocess.STDOUT, check=False)
    text = (out / "log.checkMesh").read_text()
    flags = [ln.strip() for ln in text.splitlines() if ln.strip().startswith("***")]
    verdict = "Mesh OK" if "Mesh OK" in text else "; ".join(flags) or "checkMesh failed, see log.checkMesh"
    return summary + f"; OpenFOAM case {out / 'case.foam'}; checkMesh: {verdict}"


def preview_mesh(info: dict, out: Path, n: int) -> str:
    """The block layout as a coarse mesh: n cells along every block edge, so
    each block shows as an n x n x n lattice (prisms where a block is one).
    Equal counts satisfy every transfinite constraint, so no spacing is
    needed."""
    for _, c in gmsh.model.getEntities(1):
        if Curve(c).length > 1e-9:                     # not OCC's pole edge at the apex
            gmsh.model.mesh.setTransfiniteCurve(c, n + 1)
    set_transfinite_blocks(info)
    gmsh.model.mesh.generate(3)
    return write_case(info, out)


def write_vtk(path: Path, info: dict) -> tuple[dict[str, int], dict[str, int]]:
    """Write every meshed face as a binary legacy VTK file for ParaView.

    Quads and triangles. Cell data: group (wall, seams, outer patches,
    internal block faces), zone (layout order), face (model tag) and minSJ
    (element scaled Jacobian). Returns the group and zone legends.
    """
    names = sorted(info["groups"]) + ["block_face"]
    group_of = {f: names.index(k) for k, fs in info["groups"].items() for f in fs}
    zone_ids = {name: i for i, (name, *_) in enumerate(info["layout"])}
    node_ids, points, cells, cell = {}, [], [], {"group": [], "zone": [], "face": [], "minSJ": []}
    for _, f in gmsh.model.getEntities(2):
        owners = gmsh.model.getAdjacencies(2, f)[0]
        zone = min(zone_ids[info["zones"][v]] for v in owners) if len(owners) else -1
        for et, tags, nodes in zip(*gmsh.model.mesh.getElements(2, f)):
            k = gmsh.model.mesh.getElementProperties(et)[3]
            if k not in (3, 4):
                continue
            q = gmsh.model.mesh.getElementQualities(tags, "minSJ")
            for e in range(len(tags)):
                ids = []
                for n in nodes[k * e:k * e + k]:
                    if n not in node_ids:
                        node_ids[n] = len(points)
                        points.append(gmsh.model.mesh.getNode(n)[0])
                    ids.append(node_ids[n])
                cells.append(ids)
                cell["group"].append(group_of.get(f, len(names) - 1))
                cell["zone"].append(zone)
                cell["face"].append(f)
                cell["minSJ"].append(q[e])
    size = sum(len(c) + 1 for c in cells)
    with open(path, "wb") as out:
        out.write(b"# vtk DataFile Version 3.0\nzoned mesh faces\nBINARY\nDATASET UNSTRUCTURED_GRID\n")
        out.write(f"POINTS {len(points)} double\n".encode())
        out.write(struct.pack(f">{3 * len(points)}d", *[c for p in points for c in p]))
        out.write(f"\nCELLS {len(cells)} {size}\n".encode())
        out.write(struct.pack(f">{size}i", *[v for c in cells for v in (len(c), *c)]))
        out.write(f"\nCELL_TYPES {len(cells)}\n".encode())
        out.write(struct.pack(f">{len(cells)}i", *[9 if len(c) == 4 else 5 for c in cells]))
        out.write(f"\nCELL_DATA {len(cells)}\n".encode())
        for key in ("group", "zone", "face"):
            out.write(f"SCALARS {key} int 1\nLOOKUP_TABLE default\n".encode())
            out.write(struct.pack(f">{len(cells)}i", *cell[key]))
            out.write(b"\n")
        out.write(b"SCALARS minSJ double 1\nLOOKUP_TABLE default\n")
        out.write(struct.pack(f">{len(cells)}d", *cell["minSJ"]))
        out.write(b"\n")
    return {k: i for i, k in enumerate(names)}, zone_ids


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--gui", action="store_true", help="open the result in the Gmsh GUI")
    ap.add_argument("--geometry-only", action="store_true",
                    help="stop after the block layout and write a coarse preview mesh to <case>/preview")
    ap.add_argument("--generic", action="store_true",
                    help="mesh with the generic spacing engine and write the volume mesh to <case>/preview "
                         "(the finned case always uses it)")
    ap.add_argument("--budget", action="store_true",
                    help="after the block checks, print edge families and cells per zone (step F5a)")
    ap.add_argument("--preview-cells", type=int, default=4, metavar="N",
                    help="cells along every block edge in the preview mesh (default 4)")
    ap.add_argument("--volume", action="store_true", help="also mesh the volume and write mesh.msh (step 3)")
    ap.add_argument("--fins", action="store_true",
                    help="build only the fins, check them against the scad, write fins.stl (fin step F1)")
    args = ap.parse_args(argv)
    case = args.case
    try:
        props = (case / "constant" / "caseProperties").read_text()
        D = parse_scalar(props, "D")
        N = int(parse_scalar(props, "N"))
        if N != 0 and not (args.fins or args.geometry_only or args.budget or args.generic or args.volume):
            raise ValueError(f"N = {N}: with fins only --fins and --geometry-only work so far")
        params = read_params(case)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    body = Body(D)
    if args.fins:
        spec = FinSpec(body, N, parse_scalar(props, "xi"), parse_scalar(props, "L"))
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        try:
            fins = build_fins(spec)
            print("\n".join(fin_check(spec, fins, case, {"D": D})))
            gmsh.write(str(case / "fins.brep"))
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
            gmsh.model.mesh.generate(2)
            gmsh.write(str(case / "fins.stl"))
            print(f"wrote {case / 'fins.stl'} (and the scad's {case / 'fins_scad.stl'}) for ParaView")
            if args.gui:
                gmsh.fltk.run()
        finally:
            gmsh.finalize()
        return 0
    spec = FinSpec(body, N, parse_scalar(props, "xi"), parse_scalar(props, "L")) if N else None
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        info = build(body, params, spec)
        gmsh.write(str(case / "geometry.brep"))
        t = info["tip"]
        print(f"nose core: tip half-angle {t['beta']:.2f} deg, patch edge x {t['x_a']:.3f}-{t['x_45']:.3f} mm, "
              f"half-width {t['h_in']:.2f} mm at the zone start")
        print("face groups:")
        for key, surfs in sorted(info["groups"].items()):
            print(f"  {key:18s} {len(surfs)} faces")
        print("checks:")
        print("\n".join(check(body, params, info)))
        if args.budget:
            print("budget:")
            print("\n".join(budget(body, params, info, spec)))
            return 0
        if args.generic or (args.volume and spec):
            # The generic engine: always for fins; --volume writes <case>/mesh.msh
            # for rebuild-mesh.sh, otherwise a converted preview case.
            print("generic mesh:")
            out = case if args.volume else case / "preview"
            print("\n".join(generic_mesh(body, params, info, spec, out, convert=not args.volume)))
            return 0
        if args.geometry_only:
            out = preview_mesh(info, case / "preview", args.preview_cells)
            print(f"preview mesh: {out}")
        else:
            sm = surface_mesh(body, params, info)
            print("surface mesh:")
            print("\n".join(surface_checks(body, params, info, sm)))
            legend, zone_ids = write_vtk(case / "surface.vtk", info)
            print(f"wrote {case / 'surface.vtk'}; group: " + ", ".join(f"{i} {k}" for k, i in legend.items())
                  + "; zone: " + ", ".join(f"{i} {k}" for k, i in zone_ids.items()))
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
