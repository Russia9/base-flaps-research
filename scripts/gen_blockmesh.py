#!/usr/bin/env python3
"""Generate a fully structured blockMeshDict for the bare ogive-cylinder body.

Usage:
    python3 scripts/gen_blockmesh.py openfoam/<case>

Reads D and N from <case>/constant/caseProperties and the resolution knobs
from <case>/system/blockMeshParams, and writes <case>/system/blockMeshDict.
The written dictionary is generated: edit the params and rerun, never the
output. Stage 1 supports N = 0 only.

Geometry follows geometry/arc_stabilizers.scad (secant ogive, 10 D body)
except the nose: the tip is sharp, as in the report. The R_nose cap was a
snappy-layers concession only.

Topology, all units mm. Eight 45-degree sectors, a 2x2 core where the axis
is free (upstream of the tip and behind the base):
  Z0  inlet -> tip:  core tube whose downstream face is the ogive tip patch
                     (apex is an ordinary vertex shared by the four
                     sub-blocks), plus a ring out to the farfield
  Z1  tip -> base:   ring from the body to the farfield, one block per
                     ogive station (<= ogiveBlock long, which keeps
                     blockMesh's edge interpolation within ~0.3 um of the
                     surface), then one cylinder block
  Z2  base -> outlet: core, inner ring (core -> r = R), outer ring (R -> farfield)
The Z0/Z1 interface leaves the tip patch along the wall normal; the core
side leaves it at the bisector of the remaining angle.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from create_case import parse_scalar

NSECT = 8
DTH = 2 * math.pi / NSECT
SQ2 = math.sqrt(2.0)


# ── geometry ─────────────────────────────────────────────────────────────────

class Body:
    def __init__(self, D: float):
        self.R = D / 2
        self.rho = 8.5 * D
        self.l_ogive = 2.0 * D
        self.total = 10.0 * D
        # Secant ogive through (0, 0) and (l_ogive, R), as in the scad.
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

    def slope(self, x: float) -> float:
        """Meridian angle dr/dx as an angle, rad."""
        if x >= self.l_ogive:
            return 0.0
        return math.atan((self.xc - x) / math.sqrt(self.rho**2 - (x - self.xc) ** 2))

    def x_of_r(self, r: float) -> float:
        lo, hi = 0.0, self.l_ogive
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.r(mid) < r:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def cyl(x: float, r: float, th: float) -> tuple[float, float, float]:
    return (x, r * math.cos(th), r * math.sin(th))


def theta(j: int) -> float:
    return (j % NSECT) * DTH


def sq_r(s: float, j: int) -> float:
    """Radius of the core square's j-th point (even: side midpoint, odd: corner)."""
    return s if j % 2 == 0 else s * SQ2


def sq_pt(s: float, j: int) -> tuple[float, float]:
    r = sq_r(s, j)
    return (r * math.cos(theta(j)), r * math.sin(theta(j)))


def lerp(a, b, t):
    return tuple(ai + (bi - ai) * t for ai, bi in zip(a, b))


def dist(a, b) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def polyline_length(pts) -> float:
    return sum(dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1))


# ── cell spacing ─────────────────────────────────────────────────────────────

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
    # Absorb the remainder in the cells outside both layer stacks.
    lo = min(a.n_layers, len(A))
    hi = len(sizes) - min(b.n_layers, len(B))
    free = sum(sizes[lo:hi])
    if free <= 0:
        return [s * length / sum(sizes) for s in sizes]
    scale = 1 + gap / free
    return sizes[:lo] + [s * scale for s in sizes[lo:hi]] + sizes[hi:]


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


def grading(sizes: list[float]) -> str:
    """Compress a size list into blockMesh multiGrading runs."""
    n = len(sizes)
    length = sum(sizes)
    runs = []
    i = 0
    while i < n:
        j = i + 1
        if j < n:
            q = sizes[j] / sizes[i]
            while j + 1 < n and abs(sizes[j + 1] / sizes[j] - q) < 1e-4 * q:
                j += 1
        runs.append((i, j))  # cells i..j inclusive
        i = j + 1
    if len(runs) == 1:
        exp = sizes[-1] / sizes[0]
        return "1" if abs(exp - 1) < 1e-9 else f"{exp:.8g}"
    parts = []
    for i, j in runs:
        seg = sizes[i:j + 1]
        parts.append(f"({sum(seg) / length:.10g} {len(seg) / n:.10g} {seg[-1] / seg[0]:.8g})")
    return "(" + " ".join(parts) + ")"


# ── mesh assembly ────────────────────────────────────────────────────────────

# Hex vertex pairs in edgeGrading order.
HEX_EDGES = [(0, 1), (3, 2), (7, 6), (4, 5),
             (0, 3), (1, 2), (5, 6), (4, 7),
             (0, 4), (1, 5), (2, 6), (3, 7)]
FACES = {"x1min": (0, 4, 7, 3), "x1max": (1, 2, 6, 5),
         "x2min": (0, 1, 5, 4), "x2max": (3, 7, 6, 2),
         "x3min": (0, 3, 2, 1), "x3max": (4, 5, 6, 7)}


class Mesh:
    def __init__(self):
        self.points: list[tuple[float, float, float]] = []
        self.index: dict[tuple, int] = {}
        self.edges: dict[tuple[int, int], tuple[list[float], tuple | None]] = {}
        self.blocks: list[tuple[list[int], tuple[int, int, int], list[str]]] = []
        self.patches: dict[str, list[tuple[int, ...]]] = {}

    def v(self, key: tuple, p=None) -> int:
        if len(key) == 2:   # (name, sector): sectors wrap around
            key = (key[0], key[1] % NSECT)
        if key not in self.index:
            if p is None:
                raise KeyError(f"vertex {key} used before it was defined")
            self.index[key] = len(self.points)
            self.points.append(tuple(p))
        return self.index[key]

    def p(self, key: tuple):
        return self.points[self.v(key)]

    def edge(self, a: tuple, b: tuple, sizes: list[float], curve: tuple | None = None):
        ia, ib = self.v(a), self.v(b)
        if (ia, ib) in self.edges or (ib, ia) in self.edges:
            raise ValueError(f"edge {a}-{b} defined twice")
        self.edges[(ia, ib)] = (sizes, curve)

    def sizes(self, ia: int, ib: int) -> list[float]:
        if (ia, ib) in self.edges:
            return self.edges[(ia, ib)][0]
        if (ib, ia) in self.edges:
            return self.edges[(ib, ia)][0][::-1]
        raise KeyError(f"edge {ia}-{ib} has no spacing")

    def block(self, keys: list[tuple], patches: dict[str, str]):
        vs = [self.v(k) for k in keys]
        grads = []
        counts = []
        for d in range(3):
            ns = {len(self.sizes(vs[a], vs[b])) for a, b in HEX_EDGES[4 * d:4 * d + 4]}
            if len(ns) != 1:
                raise ValueError(f"block {keys}: direction {d} edge counts differ: {ns}")
            counts.append(ns.pop())
        for a, b in HEX_EDGES:
            grads.append(grading(self.sizes(vs[a], vs[b])))
        if signed_volume([self.points[i] for i in vs]) <= 0:
            raise ValueError(f"block {keys} is inside-out or degenerate")
        self.blocks.append((vs, tuple(counts), grads))
        for face, patch in patches.items():
            self.patches.setdefault(patch, []).append(tuple(vs[i] for i in FACES[face]))


def signed_volume(p) -> float:
    """Hex volume from straight edges (6-tet split); positive if right-handed."""
    def tet(a, b, c, d):
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        w = [d[i] - a[i] for i in range(3)]
        return (u[0] * (v[1] * w[2] - v[2] * w[1])
                - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6
    return sum(tet(p[a], p[b], p[c], p[d]) for a, b, c, d in
               [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)])


# ── params ───────────────────────────────────────────────────────────────────

PARAMS = {
    "nTheta": 13, "firstLayer": 0.003, "layerRatio": 1.2, "nLayers": 20,
    "outerRatio": 1.1, "tipRadius": 0.5, "coreInlet": 40.0, "coreWake": 20.0,
    "xInlet": -800.0, "xOutlet": 5600.0, "rFar": 600.0,
    "hInlet": 20.0, "hOgive": 1.0, "hShoulder": 0.5, "hCyl": 2.5,
    "hWake": 50.0, "hFar": 30.0, "hOutletShear": 1.0, "ogiveBlock": 5.0,
}


def read_params(case: Path) -> dict[str, float]:
    text = (case / "system" / "blockMeshParams").read_text()
    params = {}
    for name, default in PARAMS.items():
        try:
            params[name] = parse_scalar(text, name)
        except ValueError:
            params[name] = default
    params["nTheta"] = int(params["nTheta"])
    params["nLayers"] = int(params["nLayers"])
    return params


# ── topology ─────────────────────────────────────────────────────────────────

def build(body: Body, P: dict) -> tuple[Mesh, dict]:
    m = Mesh()
    nt = P["nTheta"]
    q = P["outerRatio"]
    R, x_in, x_out, r_ff = body.R, P["xInlet"], P["xOutlet"], P["rFar"]
    a, A_in, s_w = P["tipRadius"], P["coreInlet"], P["coreWake"]
    wall = End(P["firstLayer"], q, P["nLayers"], P["layerRatio"])
    J = range(NSECT)

    def ring_arc(x, r, j):
        return ("arc", cyl(x, r, theta(j) + DTH / 2))

    # ── tip patch and Z0/Z1 interface ──
    xP = {j: body.x_of_r(sq_r(a, j)) for j in J}
    beta = body.slope(xP[0])
    m.v(("apex",), (0.0, 0.0, 0.0))
    for j in J:
        m.v(("P", j), (xP[j], *sq_pt(a, j)))
    # Interface: straight along the tip normal, all farfield ends at one x.
    x_F = xP[0] - (r_ff - a) * math.tan(beta)
    for j in J:
        m.v(("F", j), cyl(x_F, r_ff, theta(j)))

    # ── axial stations on the body ──
    h_tip = math.hypot(xP[0], a) / nt
    ogive = size_driven(body.l_ogive - xP[0], End(h_tip, q), End(P["hShoulder"], q), P["hOgive"])
    nodes = [xP[0]]
    for s in ogive:
        nodes.append(nodes[-1] + s)
    nodes[-1] = body.l_ogive
    st_idx = [0]
    for i, x in enumerate(nodes[1:-1], 1):
        if x - nodes[st_idx[-1]] >= P["ogiveBlock"] and x > xP[1] + 0.5 * P["ogiveBlock"]:
            st_idx.append(i)
    st_idx.append(len(nodes) - 1)
    cyl_sizes = size_driven(body.total - body.l_ogive, End(P["hShoulder"], q), wall, P["hCyl"])
    st_x = [nodes[i] for i in st_idx] + [body.total]
    nst = len(st_x)

    def B(st, j):
        return ("P", j % NSECT) if st == 0 else ("B", st, j % NSECT)

    def FF(st, j):
        return ("F", j % NSECT) if st == 0 else ("FF", st, j % NSECT)

    for st in range(1, nst):
        for j in J:
            m.v(B(st, j), cyl(st_x[st], body.r(st_x[st]), theta(j)))
            m.v(FF(st, j), cyl(st_x[st], r_ff, theta(j)))

    # Radial count from the base-plane edge, wall-clustered on every body station.
    radial = size_driven(r_ff - R, wall, End(P["hFar"], q), P["hFar"])
    nR = len(radial)
    for st in range(nst):
        for j in J:
            L = dist(m.p(B(st, j)), m.p(FF(st, j)))
            m.edge(B(st, j), FF(st, j), fixed_n(L, nR, wall))

    # Circumferential edges on stations: body (tip curve or arcs) and farfield arcs.
    for j in J:
        pa, pb = sq_pt(a, j), sq_pt(a, j + 1)
        curve = []
        for i in range(1, 20):
            y, z = lerp(pa, pb, i / 20)
            curve.append((body.x_of_r(math.hypot(y, z)), y, z))
        m.edge(B(0, j), B(0, j + 1), [1.0] * nt, ("poly", curve))
        m.edge(FF(0, j), FF(0, j + 1), [1.0] * nt, ring_arc(x_F, r_ff, j))
        for st in range(1, nst):
            m.edge(B(st, j), B(st, j + 1), [1.0] * nt, ring_arc(st_x[st], body.r(st_x[st]), j))
            m.edge(FF(st, j), FF(st, j + 1), [1.0] * nt, ring_arc(st_x[st], r_ff, j))

    # Axial edges along the body (meridians) and farfield, Z1.
    for st in range(nst - 1):
        if st < nst - 2:
            seg = ogive[st_idx[st]:st_idx[st + 1]]
        else:
            seg = cyl_sizes
        for j in J:
            x0 = m.p(B(st, j))[0]
            x1 = st_x[st + 1]
            curve = None
            if x1 <= body.l_ogive:
                curve = ("poly", [cyl(x, body.r(x), theta(j))
                                  for x in (x0 + (x1 - x0) * i / 24 for i in range(1, 24))])
            m.edge(B(st, j), B(st + 1, j), list(seg), curve)
            m.edge(FF(st, j), FF(st + 1, j), [1.0] * len(seg))

    for st in range(nst - 1):
        for j in J:
            patches = {"x1min": "fuselage", "x1max": "farfield"}
            m.block([B(st, j), FF(st, j), FF(st, j + 1), B(st, j + 1),
                     B(st + 1, j), FF(st + 1, j), FF(st + 1, j + 1), B(st + 1, j + 1)], patches)

    # ── Z0: inlet core and ring ──
    m.v(("I_ax",), (x_in, 0.0, 0.0))
    for j in J:
        m.v(("I_sq", j), (x_in, *sq_pt(A_in, j)))
        m.v(("I_ff", j), cyl(x_in, r_ff, theta(j)))

    axis0 = size_driven(-x_in, End(P["hInlet"], q), wall, P["hInlet"])
    nX0 = len(axis0)
    m.edge(("I_ax",), ("apex",), axis0)
    t_core = math.pi * 3 / 4 + beta   # bisector of normal and upstream wall
    for j in J:
        p0 = (xP[j], sq_r(a, j))
        p3 = (x_in, sq_r(A_in, j))
        d = dist(p0, p3) / 4
        p1 = (p0[0] + d * math.cos(t_core), p0[1] + d * math.sin(t_core))
        p2 = (p3[0] + d, p3[1])
        bez = []
        for i in range(48, 0, -1):   # from inlet to tip
            t = i / 48
            u = 1 - t
            x = u**3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t**3 * p3[0]
            r = u**3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t**3 * p3[1]
            bez.append(cyl(x, r, theta(j)))
        L = polyline_length([m.p(("I_sq", j))] + bez[1:] + [m.p(("P", j))])
        m.edge(("I_sq", j), ("P", j), fixed_n(L, nX0, None, wall), ("poly", bez[1:]))
        m.edge(("I_ff", j), ("F", j), fixed_n(x_F - x_in, nX0))
        L = r_ff - sq_r(A_in, j)
        m.edge(("I_sq", j), ("I_ff", j), fixed_n(L, nR, End(A_in / nt, q)))
        m.edge(("I_sq", j), ("I_sq", j + 1), [1.0] * nt)
        m.edge(("I_ff", j), ("I_ff", j + 1), [1.0] * nt, ring_arc(x_in, r_ff, j))
    for k in range(4):
        jm = 2 * k
        m.edge(("I_ax",), ("I_sq", jm), [1.0] * nt)
        curve = [cyl(x, body.r(x), theta(jm)) for x in (xP[jm] * i / 12 for i in range(1, 12))]
        m.edge(("apex",), ("P", jm), [1.0] * nt, ("poly", curve))

    for k in range(4):
        j = 2 * k
        m.block([("I_ax",), ("I_sq", j), ("I_sq", j + 1), ("I_sq", (j + 2) % NSECT),
                 ("apex",), ("P", j), ("P", j + 1), ("P", (j + 2) % NSECT)],
                {"x3min": "inlet", "x3max": "fuselage"})
    for j in J:
        jn = (j + 1) % NSECT
        m.block([("I_sq", j), ("I_ff", j), ("I_ff", jn), ("I_sq", jn),
                 ("P", j), ("F", j), ("F", jn), ("P", jn)],
                {"x3min": "inlet", "x1max": "farfield"})

    # ── Z2: base disk and wake ──
    last = nst - 1
    m.v(("W_ax",), (body.total, 0.0, 0.0))
    m.v(("O_ax",), (x_out, 0.0, 0.0))
    for j in J:
        m.v(("W_sq", j), (body.total, *sq_pt(s_w, j)))
        m.v(("O_sq", j), (x_out, *sq_pt(s_w, j)))
        m.v(("O_r", j), cyl(x_out, R, theta(j)))
        m.v(("O_ff", j), cyl(x_out, r_ff, theta(j)))

    wake = size_driven(x_out - body.total, wall, End(P["hWake"], q), P["hWake"])
    nX2 = len(wake)
    L = x_out - body.total
    h_ff_base = m.sizes(m.v(FF(last - 1, 0)), m.v(FF(last, 0)))[-1]
    m.edge(("W_ax",), ("O_ax",), wake)
    inner = size_driven(R - s_w, End(s_w / nt, q), wall, R)
    nRw = len(inner)
    for j in J:
        m.edge(("W_sq", j), ("O_sq", j), wake)
        m.edge(B(last, j), ("O_r", j), wake)
        m.edge(FF(last, j), ("O_ff", j), fixed_n(L, nX2, End(h_ff_base)))
        m.edge(("W_sq", j), B(last, j), fixed_n(R - sq_r(s_w, j), nRw, None, wall))
        m.edge(("O_sq", j), ("O_r", j), fixed_n(R - sq_r(s_w, j), nRw))
        m.edge(("O_r", j), ("O_ff", j), fixed_n(r_ff - R, nR, End(P["hOutletShear"], q)))
        for side in ("W_sq", "O_sq"):
            m.edge((side, j), (side, j + 1), [1.0] * nt)
        m.edge(("O_r", j), ("O_r", j + 1), [1.0] * nt, ring_arc(x_out, R, j))
        m.edge(("O_ff", j), ("O_ff", j + 1), [1.0] * nt, ring_arc(x_out, r_ff, j))
    for k in range(4):
        for ax, side in (("W_ax", "W_sq"), ("O_ax", "O_sq")):
            m.edge((ax,), (side, 2 * k), [1.0] * nt)

    for k in range(4):
        j = 2 * k
        jn = (j + 2) % NSECT
        m.block([("W_ax",), ("W_sq", j), ("W_sq", j + 1), ("W_sq", jn),
                 ("O_ax",), ("O_sq", j), ("O_sq", j + 1), ("O_sq", jn)],
                {"x3min": "fuselage", "x3max": "outlet"})
    for j in J:
        jn = (j + 1) % NSECT
        m.block([("W_sq", j), B(last, j), B(last, jn), ("W_sq", jn),
                 ("O_sq", j), ("O_r", j), ("O_r", jn), ("O_sq", jn)],
                {"x3min": "fuselage", "x3max": "outlet"})
        m.block([B(last, j), FF(last, j), FF(last, jn), B(last, jn),
                 ("O_r", j), ("O_ff", j), ("O_ff", jn), ("O_r", jn)],
                {"x3max": "outlet", "x1max": "farfield"})

    info = {
        "tip half-angle (deg)": math.degrees(beta),
        "tip patch x (mm)": xP[0],
        "interface farfield x (mm)": x_F,
        "ogive stations": nst - 2,
        "nTheta x 8": nt * NSECT,
        "nR (body -> farfield)": nR,
        "nR wake core -> R": nRw,
        "nX upstream": nX0,
        "nX ogive": len(ogive),
        "nX cylinder": len(cyl_sizes),
        "nX wake": nX2,
        "max radial cell (mm)": max(radial),
        "wall stack (mm)": sum(radial[:P["nLayers"]]),
    }
    return m, info


# ── output ───────────────────────────────────────────────────────────────────

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2512                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// GENERATED by scripts/gen_blockmesh.py from constant/caseProperties and
// system/blockMeshParams. Do not edit; change the params and regenerate.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def fmt(p) -> str:
    return "(" + " ".join(f"{c:.9g}" for c in p) + ")"


def write(m: Mesh, path: Path) -> int:
    out = [HEADER, "scale   0.001;\n\nvertices\n("]
    out += [f"    {fmt(p)}  // {i}" for i, p in enumerate(m.points)]
    out.append(");\n\nblocks\n(")
    cells = 0
    for vs, n, grads in m.blocks:
        cells += n[0] * n[1] * n[2]
        out.append(f"    hex ({' '.join(map(str, vs))}) ({n[0]} {n[1]} {n[2]})")
        out.append(f"        edgeGrading ({' '.join(grads)})")
    out.append(");\n\nedges\n(")
    for (ia, ib), (_, curve) in m.edges.items():
        if curve is None:
            continue
        if curve[0] == "arc":
            out.append(f"    arc {ia} {ib} {fmt(curve[1])}")
        else:
            out.append(f"    polyLine {ia} {ib} ({' '.join(fmt(p) for p in curve[1])})")
    out.append(");\n\nboundary\n(")
    kinds = {"fuselage": "wall"}
    for name in ("inlet", "outlet", "farfield", "fuselage"):
        out.append(f"    {name}\n    {{\n        type {kinds.get(name, 'patch')};\n        faces\n        (")
        out += [f"            ({' '.join(map(str, f))})" for f in m.patches[name]]
        out.append("        );\n    }")
    out.append(");\n\nmergePatchPairs\n(\n);\n\n"
               "// ************************************************************************* //\n")
    path.write_text("\n".join(out))
    return cells


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    args = ap.parse_args(argv)
    case = args.case
    try:
        props = (case / "constant" / "caseProperties").read_text()
        D = parse_scalar(props, "D")
        N = int(parse_scalar(props, "N"))
        if N != 0:
            raise ValueError(f"N = {N}: the structured generator supports the bare body (N = 0) only")
        params = read_params(case)
        mesh, info = build(Body(D), params)
        cells = write(mesh, case / "system" / "blockMeshDict")
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for k, v in info.items():
        print(f"{k:28s}: {v:.4g}" if isinstance(v, float) else f"{k:28s}: {v}")
    print(f"{'blocks':28s}: {len(mesh.blocks)}")
    print(f"{'cells':28s}: {cells}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
