"""Zones, their hexahedral blocks, and the named faces.

Groups are meshed separately and joined by flat, non-conformal seams; the
zones of one group share their faces. Seam rules (from the stitchMesh
failures of the cylindrical-seam layout):
  1. a seam's outline is straight (polygons), so both sides trace it alike
     whatever their nodes;
  2. one side of every seam is coarse;
  3. no seam touches a wall: every zone at a wall is in the group "inner".

Bare body: inner = nose | body | aft inside the polygon P (r = rZone),
far outside it. With fins: inner = the sleeve (nose | body inside the
polygon Gamma, up to xAftStart) and the fin region (inside P,
xAftStart..xSlabEnd); outer = between Gamma and P ahead of the fins;
wake = inside P behind xSlabEnd; far = outside P. See finlayout for the
fin region's blocks.
"""

from __future__ import annotations

import math

import gmsh

from . import occ as O
from .finlayout import FinLayout
from .polar import Core, at, polygon_area, ray_polygon
from .shapes import Body, FinSpec, bisect

SQ2 = math.sqrt(2.0)
RANK = {"far": 0, "outer": 1, "wake": 1, "inner": 2}       # seam master: the lower rank (coarser side)


def starting_at(rays: list[float], i: int) -> list[float]:
    """The rays rotated to start at index i (angles kept increasing)."""
    return rays[i:] + [a + 2 * math.pi for a in rays[:i]]


def inside(p, poly) -> bool:
    """Point-in-polygon (even-odd)."""
    y, z = p
    c = False
    for (y0, z0), (y1, z1) in zip(poly, poly[1:] + poly[:1]):
        if (z0 > z) != (z1 > z) and y < y0 + (z - z0) * (y1 - y0) / (z1 - z0):
            c = not c
    return c


class Layout:
    """The cross-section geometry shared by the builder and the checks."""

    def __init__(self, body: Body, P: dict, spec: FinSpec | None):
        self.body, self.P, self.spec = body, P, spec
        self.x0, self.xn, self.xs, self.x1 = P["xZoneStart"], P["xNoseEnd"], P["xAftStart"], P["xZoneEnd"]
        self.xin, self.xout, self.rff = P["xInlet"], P["xOutlet"], P["rFar"]
        if spec:
            fl = self.fins = FinLayout(spec, P)
            self.rays, self.poly_p, self.gamma = fl.rays, fl.poly_p, fl.gamma
            self.p_rays = fl.p_rays
            self.core_rays = starting_at(self.rays, 1)          # corners on the bisector rays
            self.p_core_rays = starting_at(self.p_rays, 2)      # corners on C
            self.seam = fl.fin_rays[3]                          # body revolve seam on the gap ray
            self.xe = fl.xe
            if not body.total < self.xe < self.x1:
                raise ValueError("xSlabEnd must lie between the base and xZoneEnd")
        else:
            self.fins = None
            self.rays = [k * math.pi / 4 for k in range(8)]
            self.poly_p = [at(P["rZone"], a) for a in self.rays]
            self.gamma = None
            self.p_rays = self.rays
            self.core_rays = self.p_core_rays = starting_at(self.rays, 1)
            self.seam, self.xe = 0.0, self.x1
        self.sleeve = self.gamma or self.poly_p
        # nose core: a square frustum whose sides leave the tip patch at
        # 45 - beta from the axis; half-width a where it meets the ogive
        a = P["coreTip"]
        beta = math.radians(body.slope_deg(0.0))
        k = math.tan(math.pi / 4 - beta)
        self.x_a = bisect(lambda x: body.r(x) - a, 0.0, body.l_ogive)
        self.h = lambda x: a - k * (x - self.x_a)
        self.x_45 = bisect(lambda x: body.r(x) - SQ2 * self.h(x), self.x_a, self.x_a + a / k)
        self.x_end = self.x_45 + 0.5 * (self.x_a + a / k - self.x_45)
        self.beta = math.degrees(beta)

    def region(self, x: float, y: float, z: float) -> str:
        """The group whose region contains (x, y, z) (body and fins included)."""
        if not self.x0 < x < self.x1 or not inside((y, z), self.poly_p):
            return "far"
        if not self.spec:
            return "inner"
        if x > self.xe:
            return "wake"
        if x < self.xs and not inside((y, z), self.gamma):
            return "outer"
        return "inner"


def build(body: Body, P: dict, spec: FinSpec | None = None) -> dict:
    L = Layout(body, P, spec)
    occ = gmsh.model.occ
    R, total = body.R, body.total

    body_v = O.body_solid(body, L.seam)
    fins = O.fin_solids(spec) if spec else []
    occ.synchronize()
    v_body = occ.getMass(*body_v)
    fin_info = {}
    if fins:
        whole, _ = occ.fuse(occ.copy([body_v]), occ.copy(fins))
        occ.synchronize()
        fin_info = {"exposed_volume": occ.getMass(*whole[0]) - v_body,
                    "wall_area": sum(occ.getMass(*f) for f in gmsh.model.getBoundary(whole, oriented=False))}
        occ.remove(whole, recursive=True)
    solids = [body_v] + fins

    def overlapping(x0, x1):
        out = []
        for v in solids:
            bb = occ.getBoundingBox(*v)
            if bb[0] < x1 - O.PAD and bb[3] > x0 + O.PAD:
                out.append(v)
        return out

    def minus_solids(vol, x0, x1):
        tools = overlapping(x0, x1)
        return occ.cut([vol], tools, removeTool=False)[0] if tools else [vol]

    # zones: (name, group, x0, x1, kind) with kind the cross-section for the checks
    if spec:
        zl = [("nose", "inner", L.x0, L.xn, "sleeve"), ("body", "inner", L.xn, L.xs, "sleeve"),
              ("fins", "inner", L.xs, L.xe, "p"), ("outer", "outer", L.x0, L.xs, "p-gamma"),
              ("wake", "wake", L.xe, L.x1, "p"), ("far", "far", L.xin, L.xout, "far")]
    else:
        zl = [("nose", "inner", L.x0, L.xn, "p"), ("body", "inner", L.xn, L.xs, "p"),
              ("aft", "inner", L.xs, L.x1, "p"), ("far", "far", L.xin, L.xout, "far")]
    zones: dict[int, str] = {}
    for name, group, x0, x1, kind in zl:
        if kind == "sleeve":
            pieces = minus_solids(O.x_prism(L.sleeve, x0, x1), x0, x1)
        elif kind == "p" and group == "inner":
            pieces = minus_solids(O.x_prism(L.poly_p, x0, x1), x0, x1)
        elif kind == "p":
            pieces = [O.x_prism(L.poly_p, x0, x1)]
        elif kind == "p-gamma":
            pieces, _ = occ.cut([O.x_prism(L.poly_p, x0, x1)], [O.x_prism(L.gamma, x0, x1)])
        else:
            cyl = O.x_cylinder(x0, x1 - x0, L.rff, L.p_rays)
            pieces, _ = occ.cut([cyl], [O.x_prism(L.poly_p, L.x0, L.x1)])
        zones |= {v: name for d, v in pieces if d == 3}
    occ.remove(solids, recursive=True)
    occ.synchronize()
    zone_group = {name: g for name, g, *_ in zl}

    def of_group(g):
        return {v: n for v, n in zones.items() if zone_group[n] == g}

    def others(g):
        return {v: n for v, n in zones.items() if zone_group[n] != g}

    # ── inner group, step 1: the nose core and the tip patch ────────────
    core = Core(L.core_rays, 1.0)
    h0, he = L.h(L.x0), L.h(L.x_end)

    def scaled(p, s):
        return (p[0] * s, p[1] * s)

    frustum = O.loft_solid([[(x, *scaled(c, h)) for c in core.corners()] for x, h in ((L.x0, h0), (L.x_end, he))])
    grid = [O.loop_face([(L.x0, *scaled(p, h0)), (L.x0, *scaled(q, h0)), (L.x_end, *scaled(q, he)),
                         (L.x_end, *scaled(p, he))]) for p, q in core.grid_lines()]
    nose_rays = [O.meridional(phi, [(L.x0, core.hit(phi) * h0), (L.x_end, core.hit(phi) * he), (L.x_end, 0.0),
                                    (L.xn, 0.0), (L.xn, ray_polygon(phi, L.sleeve)), (L.x0, ray_polygon(phi, L.sleeve))])
                 for phi in L.rays]
    nose = O.fragment({v: n for v, n in zones.items() if n == "nose"}, [frustum] + grid + nose_rays)
    zones = {v: n for v, n in zones.items() if n != "nose"} | nose
    funnel = tip_funnel(L, P["xTipInterface"])

    # ── inner group, step 2: everything else ────────────────────────────
    tools = list(funnel)
    tools.append(O.section_face(L.sleeve, body.l_ogive))                      # shoulder
    core_w = Core(L.core_rays, P["coreWake"])
    for i, phi in enumerate(L.rays):
        r_out = ray_polygon(phi, L.sleeve)
        if spec and L.fins.ray_kinds[i % 4] != "g":
            x_t = L.fins.x_t
            tools.append(O.meridional(phi, [(L.xn, 0.0), (x_t, 0.0), (x_t, r_out), (L.xn, r_out)]))
        elif spec:                     # the gap ray: meridional all the way, to the base core behind the base
            d_b = Core(starting_at(L.fins.root_rays, 2), P["coreWake"]).hit(phi)
            tools.append(O.meridional(phi, [(L.xn, 0.0), (total - 0.5 * R, 0.0), (total, d_b), (L.xe, d_b),
                                            (L.xe, r_out), (L.xn, r_out)]))
        else:
            d_w = core_w.hit(phi)
            tools.append(O.meridional(phi, [(L.xn, 0.0), (total - 0.5 * R, 0.0), (total, d_w), (L.x1, d_w),
                                            (L.x1, r_out), (L.xn, r_out)]))
    if spec:
        tools += fin_region_tools(L, P, core_w)
    else:
        tools += wake_core_tools(core_w, total, L.x1, R, L.rays)
        tools.append(O.disk_hole(L.poly_p, total, R, L.rays))                 # base plane outside the base
        tools.append(O.section_face(L.poly_p, P["xWakeSplit"]))
    zones = others("inner") | O.fragment(of_group("inner"), tools)

    # ── outer, wake and far groups ──────────────────────────────────────
    if spec:
        fl = L.fins
        lines = [O.x_strip(p, q, L.x0, L.xs) for p, q in fl.outer_up_lines]
        zones = others("outer") | O.fragment(of_group("outer"), O.rotated_copies(lines, spec.N))
        core_k = Core(L.p_core_rays, P["coreWake"])
        tools = [O.meridional(phi, [(L.xe, core_k.hit(phi)), (L.x1, core_k.hit(phi)), (L.x1, P["rZone"]),
                                    (L.xe, P["rZone"])]) for phi in L.p_rays]
        tools += wake_core_tools(core_k, L.xe, L.x1, R, L.p_rays)
        zones = others("wake") | O.fragment(of_group("wake"), tools)
    core_f = Core(L.p_core_rays, P["coreFar"])
    tools = [O.meridional(phi, [(L.xin, core_f.hit(phi)), (L.xout, core_f.hit(phi)), (L.xout, L.rff),
                                (L.xin, L.rff)]) for phi in L.p_rays]
    for xa, xb in ((L.xin, L.x0), (L.x1, L.xout)):
        tools.append(O.x_prism(core_f.corners(), xa, xb))
        tools += [O.x_strip(p, q, xa, xb) for p, q in core_f.grid_lines()]
        tools.append(O.x_prism(L.poly_p, xa, xb))
    for x in (L.x0, L.x1):
        tools.append(ring_face(x, L.rff, L.p_rays, L.poly_p))
    zones = others("far") | O.fragment(of_group("far"), tools)

    groups, shared = name_faces(L, zones, zone_group)
    for name, *_ in zl:
        gmsh.model.addPhysicalGroup(3, [v for v, n in zones.items() if n == name], name=name)
    for key, surfs in groups.items():
        gmsh.model.addPhysicalGroup(2, surfs, name=key)
    return {"layout": L, "zones": zones, "zone_list": zl, "zone_group": zone_group, "groups": groups,
            "shared": shared, "fins": fin_info, "body_volume": v_body, "fins_n": spec.N if spec else 0,
            "tip": {"beta": L.beta, "x_a": L.x_a, "x_45": L.x_45, "x_end": L.x_end, "h_in": L.h(L.x0)}}


def wake_core_tools(core: Core, x0: float, x1: float, R: float, rays: list[float]) -> list:
    """Behind the base: the square core, its grid, and the cylinder r = R
    (the base rim's shear layer) between the rays."""
    tools = [O.x_prism(core.corners(), x0, x1)]
    tools += [O.x_strip(p, q, x0, x1) for p, q in core.grid_lines()]
    tools.append(O.x_cylinder(x0, x1 - x0, R, rays))
    return tools


def ring_face(x: float, r: float, angles: list[float], poly: list) -> tuple:
    """The disk of radius r at x (arcs between angles) minus the polygon."""
    occ = gmsh.model.occ
    c = O.point((x, 0, 0))
    a = list(angles) + [angles[0] + 2 * math.pi]
    p = [O.point((x, r * math.cos(t), r * math.sin(t))) for t in a[:-1]]
    arcs = [occ.addCircleArc(p[k], c, p[(k + 1) % len(p)]) for k in range(len(p))]
    occ.remove([(0, c)])
    disk = (2, occ.addPlaneSurface([occ.addCurveLoop(arcs)]))
    (face,), _ = occ.cut([disk], [O.section_face(poly, x)])
    return face


def fin_region_tools(L: Layout, P: dict, core_w: Core) -> list:
    """Cutting surfaces from the sleeve twist to xSlabEnd (finlayout)."""
    fl, s, body = L.fins, L.spec, L.body
    R, total = body.R, body.total
    xs, xe, xt, xu = fl.xs, fl.xe, fl.x_t, fl.x_u
    z0, z1, z2, z3 = s.z
    Ri, Re, Ro = s.R_in, s.R_edge, s.R_out
    tools = [O.section_face(L.gamma, x) for x in (xt, xu)]                    # sleeve twist ends
    tools += [O.section_face(L.poly_p, x) for x in (fl.xf, z0, z1, z2, fl.xb)]
    tools.append(O.disk_hole(L.poly_p, total, R, fl.root_rays))               # base plane
    # behind the base: the core, its grid, r = R, and radial separators
    # from the core to where the arcs leave r = R
    core_b = Core(starting_at(fl.root_rays, 2), P["coreWake"])
    tools += wake_core_tools(core_b, total, xe, R, fl.root_rays)
    for i, phi in enumerate(fl.root_rays):
        if fl.ray_kinds[i % 4] == "g":
            continue                                   # the gap ray runs through (built with the sleeve's)
        d = core_b.hit(phi)
        tools.append(O.meridional(phi, [(total, d), (xe, d), (xe, R), (total, R)]))
    g = L.gamma
    tools += [O.x_strip(g[i], g[(i + 1) % len(g)], xs, xe) for i in range(len(g))]

    # per fin (fin 0, then turned copies)
    c = (s.Yc, s.Zc)
    t = []
    for (rho, q), phi in zip(fl.seps, fl.fin_rays):
        a0 = s.a_body(rho)
        end = at(fl.ray_end(phi), phi)
        t += O.ruled(("line", (xt, *at(R, phi)), (xt, *end)),
                     ("arc", (xu, *s.yz(rho, a0)), (xu, *q), (xu, *c)))              # sleeve twist, from the wall
        top = s.yz(rho, fl.tip)
        if rho == Re:                                                              # LE/TE bisectors
            t += O.x_arc_sweep(c, s.yz(rho, a0 - 0.03), top, xu, z0)
        else:
            t += O.x_arc_sweep(c, s.yz(rho, a0 - 0.03), top, xu, z3)
        t += O.x_arc_sweep(c, s.yz(rho, a0), top, z3, xe)
    tip_strip = O.loop_face([s.point(fl.rho_i, fl.tip, xs), s.point(fl.rho_o, fl.tip, xs),
                             s.point(fl.rho_o, fl.tip, xe), s.point(fl.rho_i, fl.tip, xe)])
    hexagon = O.loop_face([s.point(r, fl.tip, x) for r, x in ((Re, z0), (Ro, z1), (Ro, z2), (Re, z3), (Ri, z2),
                                                               (Ri, z1))])
    t += gmsh.model.occ.cut([tip_strip], [hexagon])[0]                         # tip plane, off the tip wall
    # the column above the tip plane: its lines carried straight out to P
    img = fl.chord_image
    q = s.yz(Re, fl.tip)
    t += [O.x_strip(q, img(Re), xs, z0), O.x_strip(q, img(Re), z3, xe)]
    for rho in (Ri, Ro):
        t.append(O.x_strip(s.yz(rho, fl.tip), img(rho), z1, z2))
    for (xa, ra), (xb, rb) in (((z0, Re), (z1, Ri)), ((z0, Re), (z1, Ro)), ((z2, Ro), (z3, Re)),
                               ((z2, Ri), (z3, Re))):
        t += O.ruled(("line", (xa, *s.yz(ra, fl.tip)), (xb, *s.yz(rb, fl.tip))),
                     ("line", (xa, *img(ra)), (xb, *img(rb))))
    t += [O.x_strip(p, q_, xs, xe) for p, q_ in fl.outer_lines]
    return tools + O.rotated_copies(t, s.N)


def tip_funnel(L: Layout, x_f: float) -> list:
    """Closed ruled surface from the tip-patch loop (split at every ray) to
    the sleeve polygon at x = x_f: it separates the ring around the core
    (upstream) from the ring along the ogive."""
    occ = gmsh.model.occ
    a = L.P["coreTip"]
    tol = 0.01 * a
    n = len(L.rays)
    core = Core(L.core_rays, 1.0)
    edges = []
    for _, c in gmsh.model.getEntities(1):
        bb = gmsh.model.getBoundingBox(1, c)
        if (bb[0] > L.x_a - tol and bb[3] < L.x_45 + tol
                and max(abs(v) for v in bb[1:3] + bb[4:6]) < SQ2 * a + tol and occ.getMass(1, c) > 1e-9):
            lo, hi = gmsh.model.getParametrizationBounds(1, c)
            x, y, z = gmsh.model.getValue(1, c, [0.5 * (lo[0] + hi[0])])
            if abs(math.hypot(y, z) - core.hit(math.atan2(z, y)) * L.h(x)) < tol:   # on the frustum's side
                edges.append(c)
    if len(edges) != n:
        raise RuntimeError(f"expected {n} tip-patch edges, found {len(edges)}")

    def ends(c):
        lo, hi = gmsh.model.getParametrizationBounds(1, c)
        return [gmsh.model.getValue(1, c, [u]) for u in (lo[0], hi[0])]

    def same(p, q):
        return math.dist(p, q) < 1e-3 * a

    todo = sorted(edges, key=lambda c: math.atan2(*reversed(ends(c)[0][1:])))
    chain = [todo.pop(0)]
    tip = ends(chain[0])[1]
    while todo:
        c = next(c for c in todo if any(same(tip, e) for e in ends(c)))
        todo.remove(c)
        head, tail = ends(c)
        chain.append(c if same(tip, head) else -c)
        tip = tail if same(tip, head) else head
    loop = occ.addWire(chain)

    def on_polygon(p):
        phi = math.atan2(p[2], p[1])
        return O.point((x_f, *at(ray_polygon(phi, L.sleeve), phi)))

    first = ends(abs(chain[0]))[0 if chain[0] > 0 else 1]
    ring = [on_polygon(first)]
    for c in chain[:-1]:
        head, tail = ends(abs(c))
        ring.append(on_polygon(tail if c > 0 else head))
    # OCC orients and starts the loop wire itself: try every start and both
    # directions for the outer wire and keep the untwisted surface
    best = None
    rings = [ring[k:] + ring[:k] for k in range(n)]
    rings += [r[::-1] for r in rings]
    for r in rings:
        pts = r + [r[0]]
        segs = [occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
        surf = [dt for dt in occ.addThruSections([loop, occ.addWire(segs)], makeSolid=False, makeRuled=True)
                if dt[0] == 2]
        occ.synchronize()
        twist = 0.0
        for _, sf in surf:
            (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, sf)
            for i in range(5):
                u = u0 + (u1 - u0) * i / 4
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


def on_body(s: int, body: Body) -> bool:
    """True if face s lies on the body (ogive, cylinder or base disk)."""
    pts = sample(s)
    for x, y, z in pts:
        r = math.hypot(y, z)
        on_base = abs(x - body.total) < 1e-6 and r < body.R + 1e-6
        if not on_base and abs(r - body.r(x)) > 1e-4:
            return False
    return bool(pts)


def sample(s: int) -> list:
    """A few points inside face s (a trimmed face's parameter box reaches past the face)."""
    (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, s)
    uv = [c for f in (0.21, 0.5, 0.79) for g in (0.23, 0.5, 0.77) for c in (u0 + (u1 - u0) * f, v0 + (v1 - v0) * g)]
    pts = gmsh.model.getValue(2, s, uv)
    out = []
    for i in range(0, len(pts), 3):
        p = pts[i:i + 3]
        if gmsh.model.isInside(2, s, p, parametric=False):
            out.append(tuple(p))
    return out


def name_faces(L: Layout, zones: dict, zone_group: dict) -> tuple[dict, set]:
    """Physical face groups: walls (fuselage, stabilizers), outer patches,
    shared faces between zones of a group (dropped before writing) and
    seams seam_<master>_<slave>_<i>, one pair per plane."""
    owners: dict[int, list[int]] = {}
    for v in zones:
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            owners.setdefault(s, []).append(v)
    groups: dict[str, list[int]] = {}
    shared = set()
    seams: dict[tuple, list] = {}
    for s, own in owners.items():
        names = sorted({zones[v] for v in own})
        if len(own) == 2:
            if len(names) == 2:
                key = "shared_" + "_".join(names)
                shared.add(key)
                groups.setdefault(key, []).append(s)
            continue
        group = zone_group[names[0]]
        bb = gmsh.model.getBoundingBox(2, s)
        if bb[3] < L.xin + O.PAD:
            groups.setdefault("inlet", []).append(s)
            continue
        if bb[0] > L.xout - O.PAD:
            groups.setdefault("outlet", []).append(s)
            continue
        pts = sample(s)
        if pts and all(abs(math.hypot(p[1], p[2]) - L.rff) < O.PAD for p in pts):
            groups.setdefault("farfield", []).append(s)
            continue
        other = None
        if gmsh.model.getType(2, s) == "Plane" and pts:
            p = pts[len(pts) // 2]
            n = gmsh.model.getNormal(s, gmsh.model.getParametrization(2, s, list(p)))
            sides = {L.region(*(p[i] + sg * 1e-3 * n[i] for i in range(3))) for sg in (1, -1)}
            if len(sides) == 2:
                other = next(g for g in sides if g != group)
                nn = [round(x, 4) for x in n]
                if next(x for x in nn if abs(x) > 1e-6) < 0:
                    nn = [-x for x in nn]
                plane = (tuple(nn), round(sum(nn[i] * p[i] for i in range(3)), 2))
                pair = tuple(sorted((group, other), key=lambda g: RANK[g]))
                seams.setdefault((pair, plane), []).append((s, "master" if group == pair[0] else "slave"))
        if other is None:
            groups.setdefault("fuselage" if on_body(s, L.body) else "stabilizers", []).append(s)
    by_pair: dict[tuple, list] = {}
    for pair, plane in seams:
        by_pair.setdefault(pair, []).append(plane)
    for pair, planes in by_pair.items():
        for i, plane in enumerate(sorted(planes)):
            for s, role in seams[(pair, plane)]:
                groups.setdefault(f"seam_{pair[0]}_{pair[1]}_{i}_{role}", []).append(s)
    return groups, shared


def zone_volume(L: Layout, name: str, kind: str, x0: float, x1: float, fins: dict) -> float:
    """Analytic volume of a zone."""
    body = L.body
    lo, hi = max(0.0, x0), min(body.total, x1)
    vb = body.integrate(lo, hi)[0] if hi > lo else 0.0
    if kind == "sleeve":
        return polygon_area(L.sleeve) * (x1 - x0) - vb
    if kind == "p" and name in ("fins", "nose", "body", "aft"):
        return polygon_area(L.poly_p) * (x1 - x0) - vb - (fins.get("exposed_volume", 0.0) if name == "fins" else 0.0)
    if kind == "p":
        return polygon_area(L.poly_p) * (x1 - x0)
    if kind == "p-gamma":
        return (polygon_area(L.poly_p) - polygon_area(L.gamma)) * (x1 - x0)
    return math.pi * L.rff**2 * (x1 - x0) - polygon_area(L.poly_p) * (L.x1 - L.x0)
