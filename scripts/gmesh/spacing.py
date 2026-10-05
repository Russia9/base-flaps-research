"""Cell spacing on every block edge.

Edges in one family (opposite edges of block faces, the two edges at a
prism triangle's collapsed corner) share a cell count: the most any member
needs. Each edge's distribution comes from its two end conditions, grown
at `growth` up to the edge's size cap:
  - a wall stack where the edge leaves a wall face of its own block (or a
    flat face continuing a flat wall in its plane), or
    (through the fluid) leaves a wall vertex steeply (the 20 um tip-patch
    cell near the apex);
  - the shear-layer grading where it leaves r = R behind the base (the
    rim's stack, relaxed to hWakeRadial past xWakeSplit), hLE where it
    leaves a fin's trailing-edge wake surface;
  - the far side's seam size where it leaves a seam face steeply;
  - hShoulder on axial edges at the shoulder, hLE on axial edges at a fin's
    LE or TE apex line;
  - otherwise free (the cap).
Seams are non-conformal, so no counts or nodes are matched across them.
"""

from __future__ import annotations

import math

import gmsh

from .occ import PAD, real_edges

WALLS = ("fuselage", "stabilizers")


class End:
    """Cell sizes growing away from one end: cell k is
    h * qLayer^min(k, nLayers) * q^max(0, k - nLayers)."""

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
        na, nb = min(a.size(len(A)), hmax), min(b.size(len(B)), hmax)
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
    lo, hi = min(a.n_layers, len(A)), len(sizes) - min(b.n_layers, len(B))
    free = sum(sizes[lo:hi])
    if free <= 0:
        return [x * length / sum(sizes) for x in sizes]
    scale = 1 + gap / free
    return sizes[:lo] + [x * scale for x in sizes[lo:hi]] + sizes[hi:]


def two_sided(L: float, n: int, a: End, b: End) -> tuple[list[float], bool]:
    """n cell sizes over L from end condition a to b. Returns the sizes and
    whether both were met (False: too few cells, stacks were merged)."""
    lo, hi = min(a.h, b.h, L), L
    if len(size_driven(L, a, b, lo)) < n:                 # short edge: split the largest cells
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
    sz = size_driven(L, a, b, lo)                          # one or two too many: merge in the middle
    ok = True
    while len(sz) > n:
        la, lb = min(a.n_layers, len(sz) // 2), min(b.n_layers, len(sz) // 2)
        rng = range(la, len(sz) - lb - 1)
        if not rng:
            rng, ok = range(len(sz) - 1), False
        i = min(rng, key=lambda k: sz[k] + sz[k + 1])
        sz[i:i + 2] = [sz[i] + sz[i + 1]]
    return sz, ok


class Curve:
    """Endpoints, length and an arc-length table of one model curve."""

    def __init__(self, c: int):
        (lo,), (hi,) = gmsh.model.getParametrizationBounds(1, c)
        self.u = [lo + (hi - lo) * i / 400 for i in range(401)]
        pts = gmsh.model.getValue(1, c, self.u)
        self.pts = [pts[3 * i:3 * i + 3] for i in range(401)]
        self.s = [0.0]
        for i in range(400):
            self.s.append(self.s[-1] + math.dist(self.pts[i], self.pts[i + 1]))
        self.length = self.s[-1]
        self.p, self.q = self.pts[0], self.pts[-1]
        self.mid = self.pts[200]

    def u_at(self, s: float) -> float:
        i = min(max(1, next((k for k, v in enumerate(self.s) if v >= s), 400)), 400)
        f = (s - self.s[i - 1]) / max(self.s[i] - self.s[i - 1], 1e-300)
        return self.u[i - 1] + f * (self.u[i] - self.u[i - 1])


class Topo:
    """Block edges with their end vertices (parametric start first), the
    zones they border, and the special faces at each vertex."""

    def __init__(self, info: dict):
        self.info = info
        self.group_of = {s: k for k, ss in info["groups"].items() for s in ss}
        self.blocks: dict[int, list[int]] = {}
        self.loops: dict[int, list[int]] = {}
        self.curves: dict[int, Curve] = {}
        self.ends: dict[int, tuple[int, int]] = {}
        self.edge_zones: dict[int, set] = {}
        self.xyz: dict[int, list[float]] = {}
        for v, name in info["zones"].items():
            faces = [f for _, f in gmsh.model.getBoundary([(3, v)], oriented=False)]
            self.blocks[v] = faces
            for f in faces:
                if f not in self.loops:
                    self.loops[f] = [c for _, c in real_edges((2, f))]
                for c in self.loops[f]:
                    self.edge_zones.setdefault(c, set()).add(name)
                    if c in self.curves:
                        continue
                    cv = Curve(c)
                    pts = [p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)]
                    for p in pts:
                        self.xyz.setdefault(p, gmsh.model.getValue(0, p, []))
                    start = min(pts, key=lambda p: math.dist(self.xyz[p], cv.p))
                    end = next((p for p in pts if p != start), start)
                    self.curves[c], self.ends[c] = cv, (start, end)
        self.at_vertex: dict[int, list[int]] = {}           # wall and seam faces at each vertex
        for f, loop in self.loops.items():
            key = self.group_of.get(f, "")
            if key in WALLS or key.startswith("seam_"):
                for c in loop:
                    for p in self.ends[c]:
                        if f not in self.at_vertex.setdefault(p, []):
                            self.at_vertex[p].append(f)

    def block_edges(self, v: int) -> set[int]:
        return {c for f in self.blocks[v] for c in self.loops[f]}

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

    def axial(self, c: int) -> float:
        """|dx| / length of the chord (1: along the axis)."""
        cv = self.curves[c]
        return abs(cv.q[0] - cv.p[0]) / max(math.dist(cv.p, cv.q), 1e-12)

    def tangent(self, c: int, p: int) -> list[float]:
        pts = self.curves[c].pts
        a, b = (pts[0], pts[3]) if p == self.ends[c][0] else (pts[-1], pts[-4])
        d = [b[i] - a[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in d)) or 1.0
        return [x / n for x in d]

    def steep(self, c: int, p: int, prefix: str) -> list[int]:
        """Faces whose group starts with prefix at vertex p that edge c leaves
        within 37 deg of their normal (|cos| > 0.8), excluding faces it lies
        in. A 45-degree wedge apex (fin LE/TE) is not steep."""
        out, t = [], self.tangent(c, p)
        for f in self.at_vertex.get(p, []):
            if not self.group_of.get(f, "").startswith(prefix) or c in self.loops[f]:
                continue
            n = gmsh.model.getNormal(f, gmsh.model.getParametrization(2, f, self.xyz[p]))
            if abs(sum(n[i] * t[i] for i in range(3))) > 0.8:
                out.append(f)
        return out


class Families:
    """Union-find over block edges with orientation: one family shares a
    cell count; members' parametric directions are related by a sign."""

    def __init__(self):
        self.parent: dict[int, int] = {}
        self.sign: dict[int, int] = {}
        self.conflicts: list[tuple[int, int]] = []

    def find(self, c: int) -> tuple[int, int]:
        if self.parent.setdefault(c, c) == c:
            self.sign.setdefault(c, 1)
            return c, 1
        root, s = self.find(self.parent[c])
        self.sign[c] *= s
        self.parent[c] = root
        return root, self.sign[c]

    def union(self, a: int, b: int, s: int):
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
    u, w = [a[j] - p[j] for j in range(3)], [b[j] - p[j] for j in range(3)]
    return math.acos(max(-1.0, min(1.0, sum(u[j] * w[j] for j in range(3)) / math.dist(a, p) / math.dist(b, p))))


def edge_families(topo: Topo) -> Families:
    """Opposite edges of every quad face; the two edges at a triangle's
    collapsed corner (the widest, as export.prism_corners collapses it; the
    apex for faces with OCC's pole edge)."""
    fam = Families()
    for c in topo.curves:
        fam.find(c)
    for f, loop in topo.loops.items():
        corners, edges, signs = topo.cycle(f)
        if len(edges) == 4:
            for i in (0, 1):
                fam.union(edges[i], edges[i + 2], -signs[i] * signs[i + 2])
        elif len(edges) == 3:
            if len(gmsh.model.getBoundary([(2, f)], oriented=False)) != len(loop):     # pole edge
                k = corners.index(min(corners, key=lambda p: math.hypot(*topo.xyz[p][1:]) + abs(topo.xyz[p][0])))
            else:
                pts = [topo.xyz[p] for p in corners]
                k = max(range(3), key=lambda i: corner_angle(pts[i], pts[i - 1], pts[(i + 1) % 3]))
            fam.union(edges[k], edges[k - 1], -signs[k] * signs[k - 1])
    return fam


class Rules:
    """End conditions and size caps from gmshParams and the layout."""

    def __init__(self, topo: Topo, P: dict, info: dict):
        self.topo, self.P, self.info = topo, P, info
        L = self.L = info["layout"]
        self.body, self.spec, self.fl = L.body, L.spec, L.fins
        q = P["growth"]
        self.wall = End(P["firstLayer"], q, P["nLayers"], P["layerRatio"])
        self.tip = End(P["firstLayerTip"], q)
        self.x45 = info["tip"]["x_45"]
        self.kinds = {f: self.face_kind(f) for f in topo.loops}
        self.on_wall = {c for f, k in self.kinds.items() if k == "wall" for c in topo.loops[f]}
        self.mark_wall_extensions()

    # ── faces that grade the edges leaving them ────────────────────────
    def face_kind(self, f: int) -> str | None:
        g = self.topo.group_of.get(f, "")
        if g in WALLS:
            return "wall"
        if g.startswith("seam_"):
            return "seam"
        bb = gmsh.model.getBoundingBox(2, f)
        total, R = self.body.total, self.body.R
        if bb[0] < total - PAD:
            return None
        pts = face_points(f)
        if pts and all(abs(math.hypot(p[1], p[2]) - R) < PAD for p in pts):
            return "shear"
        if self.spec and pts and all(self.on_te_surface(p) for p in pts):
            return "te"
        return None

    def mark_wall_extensions(self):
        """Flat faces that continue a flat wall in its own plane past a shared
        edge (the tip plane beside the tip wall, the base plane past the rim)
        grade like the wall for the blocks on the wall's fluid side: those get
        the wall stack on every edge leaving the face, consistent with the
        block on the wall beside them (else their first cells fan out from
        the corner and fold)."""
        topo = self.topo
        self.ext_side: dict[int, tuple] = {}
        walls_at: dict[int, list[int]] = {}
        for f, k in self.kinds.items():
            if k == "wall" and gmsh.model.getType(2, f) == "Plane":
                for c in topo.loops[f]:
                    walls_at.setdefault(c, []).append(f)

        def normal(f):
            pts = face_points(f) or [gmsh.model.getValue(2, f, [0.5, 0.5])]
            return gmsh.model.getNormal(f, gmsh.model.getParametrization(2, f, list(pts[0])))

        for f, loop in topo.loops.items():
            if self.kinds[f] or gmsh.model.getType(2, f) != "Plane":
                continue
            ws = {w for c in loop for w in walls_at.get(c, [])}
            n = None
            for w in ws:
                n = n or normal(f)
                nw = normal(w)
                if abs(sum(n[i] * nw[i] for i in range(3))) > 0.9999:
                    # the wall's fluid side: towards its block's centre
                    v = next(v for v, faces in topo.blocks.items() if w in faces)
                    cw = gmsh.model.occ.getCenterOfMass(3, v)
                    pw = face_points(w)[0]
                    side = 1 if sum((cw[i] - pw[i]) * nw[i] for i in range(3)) > 0 else -1
                    self.kinds[f] = "wall_ext"
                    self.ext_side[f] = (pw, [side * x for x in nw])
                    break

    def on_te_surface(self, p) -> bool:
        if p[0] < self.spec.z[3] - PAD:
            return False
        _, rho, a = self.spec.frame(p[1], p[2])
        return abs(rho - self.spec.R_edge) < 1e-4 and a < self.spec.th_tip + 1e-6

    def in_hood(self, p) -> bool:
        """Inside a fin's hood (between the wrap arcs, up to the tip plane)."""
        if not self.spec:
            return False
        _, rho, a = self.spec.frame(p[1], p[2])
        return self.fl.rho_i - 1e-3 <= rho <= self.fl.rho_o + 1e-3 and a <= self.fl.tip + 1e-3

    # ── size caps ───────────────────────────────────────────────────────
    def h_axial(self, x: float, zone: str, near_fin: bool) -> float:
        P, total = self.P, self.body.total
        if zone == "far":
            return P["hFarMid"] if self.L.x0 - PAD < x < self.L.x1 + PAD else P["hFar"]
        if zone == "outer":
            return P["hRingOut"]
        if zone == "fins":
            if x > total:
                return P["hRelax"]
            return P["hFinChord"] if near_fin else P["hSlabOuter"]
        if x < 0:
            return P["hUpstream"]
        if x < self.body.l_ogive:
            return P["hOgive"]
        if x < total:
            return P["hWall"]
        x_a = self.L.xe if self.spec else P["xWakeSplit"]
        t = max(0.0, min(1.0, (x - x_a) / (self.L.x1 - x_a)))
        return P["hRelax"] + t * (P["hWake"] - P["hRelax"])

    def h_cross(self, c: int, zone: str) -> float:
        """Cap on a cross-section edge: circumferential edges r * 45 deg /
        nTheta (r at least R, so the cores match their rings), radial ones
        hRingOut; in the fin region hFinSpan inside a hood, else hSlabOuter."""
        P = self.P
        cv = self.topo.curves[c]
        if zone == "far":
            return P["hFar"]
        if zone == "outer":
            return P["hRingOut"]
        if zone == "fins":
            return P["hFinSpan"] if self.in_hood(cv.mid) else P["hSlabOuter"]
        r = math.hypot(*cv.mid[1:])
        d = [cv.q[i] - cv.p[i] for i in range(3)]
        er = [0.0, cv.mid[1] / max(r, 1e-9), cv.mid[2] / max(r, 1e-9)]
        radial = abs(sum(d[i] * er[i] for i in range(3))) / max(math.sqrt(sum(x * x for x in d)), 1e-12)
        if radial > 0.7 and r > 1e-3:
            return P["hRingOut"]
        return max(r, self.body.R) * (math.pi / 4) / P["nTheta"]

    def cap(self, c: int) -> float:
        """Size cap of edge c: the smallest over its zones and, for an edge on
        a wall, hWall (wall-tangential cells next to the 3 um first cell)."""
        if c in self.on_wall:
            return min(self._cap(c), self.P["hWall"])
        return self._cap(c)

    def _cap(self, c: int) -> float:
        topo = self.topo
        cv = topo.curves[c]
        zones = topo.edge_zones[c]
        x = cv.mid[0]
        r_tip = self.spec.tip_radius() + self.fl.wrap if self.spec else 0.0
        near_fin = bool(self.spec) and self.fl.xf - PAD <= x <= self.fl.xb + PAD and math.hypot(*cv.mid[1:]) < r_tip
        hs = []
        for z in zones:
            ha, hc = self.h_axial(x, z, near_fin), self.h_cross(c, z)
            ax = topo.axial(c)
            hs.append(ha if ax > 0.85 else (hc if ax < 0.6 else min(ha, hc)))
        return min(hs)

    # ── end conditions ──────────────────────────────────────────────────
    def from_faces(self) -> dict[tuple[int, int], str]:
        """(edge, vertex) -> kind of a graded face of the edge's own block that
        the edge leaves at that vertex (one end on the face, not in it)."""
        topo, out = self.topo, {}
        order = {"wall": 0, "wall_ext": 0, "shear": 1, "te": 2, "seam": 3}
        for v, faces in topo.blocks.items():
            edges = topo.block_edges(v)
            for f in faces:
                k = self.kinds[f]
                if not k:
                    continue
                if k == "wall_ext":                      # only for blocks on the wall's fluid side
                    pw, nf = self.ext_side[f]
                    cb = gmsh.model.occ.getCenterOfMass(3, v)
                    if sum((cb[i] - pw[i]) * nf[i] for i in range(3)) <= 0:
                        continue
                corners = {p for c in topo.loops[f] for p in topo.ends[c]}
                for c in edges - set(topo.loops[f]):
                    on = [p for p in topo.ends[c] if p in corners]
                    if len(on) == 1:
                        cur = out.get((c, on[0]))
                        if cur is None or order[k] < order[cur]:
                            out[(c, on[0])] = k
        return out

    def seam_size(self, c: int) -> End | None:
        zones = self.topo.edge_zones[c]
        return End(self.P["hSeam"], self.P["growth"]) if zones == {"far"} else None

    def needs(self) -> dict[int, tuple]:
        """Per edge: (end condition at its start, at its end, size cap)."""
        topo, P, q = self.topo, self.P, self.P["growth"]
        via = self.from_faces()
        on_wall = self.on_wall
        out = {}
        for c in topo.curves:
            conds = []
            for p in topo.ends[c]:
                x = topo.xyz[p][0]
                r = math.hypot(*topo.xyz[p][1:])
                k = via.get((c, p))
                # The steep rule is for edges through the fluid only: an edge on
                # a wall that meets another wall steeply (the body at the base
                # rim) runs behind that wall's boundary layer.
                steep_wall = c not in on_wall and (topo.steep(c, p, "fuselage") or topo.steep(c, p, "stabilizers"))
                if k in ("wall", "wall_ext") or steep_wall:
                    cond = self.tip if x < self.x45 + 1e-6 else self.wall
                elif k == "shear":
                    cond = self.wall if x <= P["xWakeSplit"] + 1e-6 else End(P["hWakeRadial"], q)
                elif k == "te":
                    cond = End(P["hLE"], q)
                elif k == "seam" or topo.steep(c, p, "seam_"):
                    cond = self.seam_size(c)
                elif topo.axial(c) > 0.9 and abs(x - self.body.l_ogive) < 1e-6 and "far" not in topo.edge_zones[c]:
                    cond = End(P["hShoulder"], q)
                elif (self.spec and topo.axial(c) > 0.9 and r < self.spec.tip_radius() + 1
                      and any(abs(x - z) < 1e-6 for z in (self.spec.z[0], self.spec.z[3]))
                      and any(topo.group_of.get(f) == "stabilizers" for f in topo.at_vertex.get(p, []))):
                    cond = End(P["hLE"], q)                              # chordwise at the LE/TE apex
                else:
                    cond = None
                conds.append(cond)
            out[c] = (conds[0], conds[1], self.cap(c))
        return out


def spacing(info: dict, P: dict, coarsen: float = 1.0) -> dict:
    """Family counts and every edge's cell sizes. coarsen > 1 multiplies
    every size cap (a quick test mesh with the real wall stacks)."""
    topo = Topo(info)
    fam = edge_families(topo)
    rules = Rules(topo, P, info)
    needs = rules.needs()
    groups = fam.groups()
    counts, sizes, short = {}, {}, []

    if coarsen != 1.0:
        needs = {c: (a, b, h * coarsen) for c, (a, b, h) in needs.items()}

    def ends_of(c):
        a, b, h = needs[c]
        return (a or End(h)), (b or End(h)), h

    for root, members in groups.items():
        n = 1
        for c in members:
            a, b, h = ends_of(c)
            n = max(n, len(size_driven(topo.curves[c].length, a, b, h)))
        counts[root] = n
        for c in members:
            a, b, _ = ends_of(c)
            sizes[c], ok = two_sided(topo.curves[c].length, n, a, b)
            if not ok:
                short.append(c)
    return {"topo": topo, "fam": fam, "needs": needs, "counts": counts, "sizes": sizes, "short": short}


def block_cells(topo: Topo, fam: Families, counts: dict) -> dict[str, int]:
    """Cells per zone: per block the product of the counts on three edges
    meeting at a corner (a prism counts half its hexahedron)."""
    cells: dict[str, int] = {}
    for v, name in topo.info["zones"].items():
        corner_edges: dict[int, list[int]] = {}
        for c in topo.block_edges(v):
            for p in topo.ends[c]:
                corner_edges.setdefault(p, []).append(c)
        three = next(es for es in corner_edges.values() if len(es) == 3)
        n = math.prod(counts[fam.find(c)[0]] for c in three)
        if len(topo.blocks[v]) == 5:
            n //= 2
        cells[name] = cells.get(name, 0) + n
    return cells


def face_points(f: int) -> list:
    """A few points inside face f."""
    (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, f)
    uv = [c for a in (0.25, 0.5, 0.75) for b in (0.3, 0.6) for c in (u0 + (u1 - u0) * a, v0 + (v1 - v0) * b)]
    pts = gmsh.model.getValue(2, f, uv)
    return [tuple(pts[i:i + 3]) for i in range(0, len(pts), 3)
            if gmsh.model.isInside(2, f, pts[i:i + 3], parametric=False)]


def report(gs: dict, top: int = 25) -> list[str]:
    topo, fam, counts, needs = gs["topo"], gs["fam"], gs["counts"], gs["needs"]
    cells = block_cells(topo, fam, counts)
    lines = [(f"  {len(topo.curves)} block edges in {len(counts)} families; {len(fam.conflicts)} orientation "
             f"conflicts; {len(gs['short'])} edges too short for both end conditions")]
    for name, n in sorted(cells.items(), key=lambda kv: -kv[1]):
        lines.append(f"  zone {name:6s} {n / 1e6:7.3f} M cells")
    lines.append(f"  total        {sum(cells.values()) / 1e6:7.3f} M cells")
    lines.append("  largest families (count, members, the edge needing most: ends W/h/-, zones, from -> to):")
    groups = fam.groups()
    for root, n in sorted(counts.items(), key=lambda kv: -kv[1])[:top]:
        def need(c):
            a, b, h = needs[c]
            return len(size_driven(topo.curves[c].length, a or End(h), b or End(h), h))
        c = max(groups[root], key=need)
        cv = topo.curves[c]
        a, b, h = needs[c]
        lab = "".join("W" if e is not None and e.n_layers else ("h" if e is not None else "-") for e in (a, b))
        lines.append(f"    {n:5d} x{len(groups[root]):4d} {lab} cap {h:5.2f} {','.join(sorted(topo.edge_zones[c]))[:18]:18s} "
                     f"({cv.p[0]:.1f}, r {math.hypot(*cv.p[1:]):.1f}) -> ({cv.q[0]:.1f}, r {math.hypot(*cv.q[1:]):.1f})"
                     f" L {cv.length:.1f}")
    return lines
