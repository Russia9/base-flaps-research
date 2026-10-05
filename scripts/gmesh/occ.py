"""OpenCASCADE primitives (Gmsh occ kernel). Lengths in mm, x along the axis.

Cross-section points are (Y, Z) pairs; 3-D points are (x, Y, Z).
"""

from __future__ import annotations

import math

import gmsh

from .shapes import Body, FinSpec

PAD = 1e-3          # OCC pads bounding boxes slightly
Dt = tuple[int, int]


def occ():
    return gmsh.model.occ


def point(p) -> int:
    return occ().addPoint(*p)


def loop_face(pts3: list) -> Dt:
    """Planar face bounded by the closed polygon pts3 (3-D points)."""
    p = [point(q) for q in pts3]
    lines = [occ().addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))]
    return (2, occ().addPlaneSurface([occ().addCurveLoop(lines)]))


def section_face(poly: list, x: float) -> Dt:
    """The cross-section polygon poly as a face in the plane x."""
    return loop_face([(x, *q) for q in poly])


def x_prism(poly: list, x0: float, x1: float) -> Dt:
    """Solid with cross-section poly between x0 and x1 (planar sides)."""
    vol = [dt for dt in occ().extrude([section_face(poly, x0)], x1 - x0, 0, 0) if dt[0] == 3]
    return vol[0]


def x_strip(p, q, x0: float, x1: float) -> Dt:
    """Planar quad: the cross-section segment p-q swept from x0 to x1."""
    return loop_face([(x0, *p), (x0, *q), (x1, *q), (x1, *p)])


def meridional(phi: float, pts: list[tuple[float, float]]) -> Dt:
    """Planar polygon of (x, r) points in the meridional half-plane at phi."""
    c, s = math.cos(phi), math.sin(phi)
    return loop_face([(x, r * c, r * s) for x, r in pts])


def arc_curve(centre, a, b, x: float) -> int:
    """Circle arc in the plane x from a to b about centre (short way)."""
    c = point((x, *centre))
    t = occ().addCircleArc(point((x, *a)), c, point((x, *b)))
    occ().remove([(0, c)])
    return t


def x_arc_sweep(centre, a, b, x0: float, x1: float) -> list[Dt]:
    """Cylindrical face: the arc a-b about centre swept from x0 to x1."""
    return [dt for dt in occ().extrude([(1, arc_curve(centre, a, b, x0))], x1 - x0, 0, 0) if dt[0] == 2]


def fin_revolve(spec: FinSpec, xr: list[tuple[float, float]], a0: float, a1: float, k: int = 0) -> list[Dt]:
    """The (x, rho) polyline of fin k's frame at angle a0, revolved about the
    fin's arc axis to angle a1 (surfaces of revolution: cones, cylinders,
    planes)."""
    pts = [point(spec.point(rho, a0, x, k)) for x, rho in xr]
    lines = [(1, occ().addLine(pts[i], pts[i + 1])) for i in range(len(pts) - 1)]
    yc, zc = spec.centre(k)
    # increasing fin-frame angle turns from +Z towards +Y: negative about +x
    return [dt for dt in occ().revolve(lines, 0, yc, zc, 1, 0, 0, -(a1 - a0)) if dt[0] == 2]


def x_cylinder(x0: float, dx: float, r: float, angles: list[float]) -> Dt:
    """Cylinder along x whose curved surface is separate patches between
    the given angles (rad, increasing). A single occ.addCylinder has one
    periodic face, and structured faces that end on its parameter seam fold;
    every patch edge must sit on a block edge."""
    a = list(angles) + [angles[0] + 2 * math.pi]
    faces = []
    for k in range(len(a) - 1):
        line = occ().addLine(point((x0, r, 0)), point((x0 + dx, r, 0)))
        occ().rotate([(1, line)], 0, 0, 0, 1, 0, 0, a[k])
        faces += [t for d, t in occ().revolve([(1, line)], 0, 0, 0, 1, 0, 0, a[k + 1] - a[k]) if d == 2]
    for x in (x0, x0 + dx):
        c = point((x, 0, 0))
        p = [point((x, r * math.cos(t), r * math.sin(t))) for t in a[:-1]]
        arcs = [occ().addCircleArc(p[k], c, p[(k + 1) % len(p)]) for k in range(len(p))]
        faces.append(occ().addPlaneSurface([occ().addCurveLoop(arcs)]))
        occ().remove([(0, c)])
    return (3, occ().addVolume([occ().addSurfaceLoop(faces, sewing=True)]))


def disk_hole(poly: list, x: float, r: float, angles: list[float]) -> Dt:
    """The polygon at x with the disk of radius r (arcs between angles) removed."""
    outer = section_face(poly, x)
    c = point((x, 0, 0))
    a = list(angles) + [angles[0] + 2 * math.pi]
    p = [point((x, r * math.cos(t), r * math.sin(t))) for t in a[:-1]]
    arcs = [occ().addCircleArc(p[k], c, p[(k + 1) % len(p)]) for k in range(len(p))]
    occ().remove([(0, c)])
    hole = (2, occ().addPlaneSurface([occ().addCurveLoop(arcs)]))
    (face,), _ = occ().cut([outer], [hole])
    return face


def ruled(w0: tuple, w1: tuple) -> list[Dt]:
    """Ruled surface between two curves, each ("line", A, B) or ("arc", A, B,
    centre) with 3-D points, ruled A0-A1 and B0-B1. OCC's ThruSections picks
    the pairing itself and can return the crossed (twisted) surface; this
    checks the side edges and rebuilds with the second curve reversed."""

    def curve(w, rev):
        kind, a, b = w[0], w[1], w[2]
        if rev:
            a, b = b, a
        pa, pb = point(a), point(b)
        if kind == "line":
            return occ().addLine(pa, pb)
        c = point(w[3])
        t = occ().addCircleArc(pa, c, pb)
        occ().remove([(0, c)])
        return t

    def near(p, q):
        return math.dist(p, q) < 1e-6

    for rev in (False, True):
        faces = [dt for dt in occ().addThruSections([occ().addWire([curve(w0, False)]), occ().addWire([curve(w1, rev)])],
                                                    makeSolid=False, makeRuled=True) if dt[0] == 2]
        occ().synchronize()
        ends = []
        for _, c in gmsh.model.getBoundary(faces, oriented=False):
            pts = [gmsh.model.getValue(0, p, []) for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)]
            if len(pts) == 2:
                ends.append(pts)
        crossed = any((near(p, w0[1]) and near(q, w1[2])) or (near(q, w0[1]) and near(p, w1[2])) for p, q in ends)
        if not crossed:
            return faces
        occ().remove(faces, recursive=True)
    raise RuntimeError(f"ruled surface between {w0} and {w1} stays twisted")


def rotated_copies(tools: list[Dt], n: int) -> list[Dt]:
    """tools plus n - 1 copies turned by 360/n about +x."""
    out = list(tools)
    for k in range(1, n):
        copy = occ().copy(tools)
        occ().rotate(copy, 0, 0, 0, 1, 0, 0, k * 2 * math.pi / n)
        out += copy
    return out


def body_solid(body: Body, angles: list[float]) -> Dt:
    """The body (apex, ogive arc, cylinder, base) as non-periodic patches
    revolved between the given angles (rad, increasing; block edges) and
    sewn: a periodic face's parameter seam folds the structured faces that
    end on it."""
    R, total = body.R, body.total
    a = list(angles) + [angles[0] + 2 * math.pi]
    faces = []
    for k in range(len(a) - 1):
        p_apex, p_sh, p_rim = (point(q) for q in ((0, 0, 0), (body.l_ogive, R, 0), (total, R, 0)))
        p_c = point((body.xc, body.yc, 0))
        prof = [(1, occ().addCircleArc(p_apex, p_c, p_sh)), (1, occ().addLine(p_sh, p_rim))]
        occ().remove([(0, p_c)])
        occ().rotate(prof, 0, 0, 0, 1, 0, 0, a[k])
        faces += [t for d, t in occ().revolve(prof, 0, 0, 0, 1, 0, 0, a[k + 1] - a[k]) if d == 2]
    c = point((total, 0, 0))
    p = [point((total, R * math.cos(t), R * math.sin(t))) for t in a[:-1]]
    arcs = [occ().addCircleArc(p[k], c, p[(k + 1) % len(p)]) for k in range(len(p))]
    occ().remove([(0, c)])
    faces.append(occ().addPlaneSurface([occ().addCurveLoop(arcs)]))
    return (3, occ().addVolume([occ().addSurfaceLoop(faces, sewing=True)]))


def fin_solids(spec: FinSpec) -> list[Dt]:
    """The N fin solids, sewn from the scad's faces: ruled LE/TE wedges between
    the edge arc and the inner/outer arcs, cylindrical barrels, and planar
    root and tip caps."""
    z0, z1, z2, z3 = spec.z
    Ri, Re, Ro = spec.R_in, spec.R_edge, spec.R_out

    def arc(r, x):
        c = point((x, spec.Yc, spec.Zc))
        a = occ().addCircleArc(point(spec.point(r, spec.th[r], x)), c, point(spec.point(r, spec.th_tip, x)))
        occ().remove([(0, c)])
        return occ().addWire([a])

    faces = []
    for (ra, xa), (rb, xb) in (((Re, z0), (Ro, z1)), ((Re, z0), (Ri, z1)), ((Ro, z1), (Ro, z2)),
                               ((Ri, z1), (Ri, z2)), ((Ro, z2), (Re, z3)), ((Ri, z2), (Re, z3))):
        faces += [t for d, t in occ().addThruSections([arc(ra, xa), arc(rb, xb)], makeSolid=False,
                                                      makeRuled=True) if d == 2]
    for a_of in (lambda r: spec.th[r], lambda r: spec.th_tip):        # root cap, tip cap
        ring = [(Re, z0), (Ro, z1), (Ro, z2), (Re, z3), (Ri, z2), (Ri, z1)]
        faces.append(loop_face([spec.point(r, a_of(r), x) for r, x in ring])[1])
    fin = (3, occ().addVolume([occ().addSurfaceLoop(faces, sewing=True)]))
    out = rotated_copies([fin], spec.N)
    occ().synchronize()
    return out


def fragment(zones: dict[int, str], tools: list[Dt]) -> dict[int, str]:
    """Fragment zone volumes with tools; return {new volume: zone}, dropping
    tool-only pieces and leftover lower-dimensional entities."""
    objects = list(zones)
    _, children = occ().fragment([(3, v) for v in objects], tools)
    new = {}
    for i, v in enumerate(objects):
        for d, t in children[i]:
            if d == 3:
                new[t] = zones[v]
    occ().synchronize()
    stray = sorted({(3, t) for kids in children[len(objects):] for d, t in kids if d == 3 and t not in new})
    if stray:
        occ().remove(stray)
    occ().synchronize()
    for dim in (2, 1, 0):
        free = [(dim, t) for _, t in gmsh.model.getEntities(dim) if len(gmsh.model.getAdjacencies(dim, t)[0]) == 0]
        if free:
            occ().remove(free)
        occ().synchronize()
    return new


def real_edges(face: Dt) -> list[Dt]:
    """Edges of a face, without OCC's zero-length pole edges (the apex)."""
    return [e for e in gmsh.model.getBoundary([face], oriented=False) if occ().getMass(1, e[1]) > 1e-9]


def loft_solid(sections: list[list]) -> Dt:
    """Ruled solid through closed polygons of 3-D points (one per section)."""
    wires = []
    for pts in sections:
        p = [point(q) for q in pts]
        wires.append(occ().addWire([occ().addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))]))
    return next(dt for dt in occ().addThruSections(wires, makeSolid=True, makeRuled=True) if dt[0] == 3)
