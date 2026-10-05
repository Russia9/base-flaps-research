"""Cross-section and chordwise layout of the finned body (wall-conformal design).

Per fin, in its own polar frame about the arc centre (rho, a):
  - the fin (R_in..R_out, root..tip) sits in a hood finWrap thick, between
    the inner and outer wrap arcs rho_i = R_in - wrap, rho_o = R_out + wrap
    and the tip plane a = tip; its corners there are It and Ot;
  - three separators run from the body out to the sleeve polygon Gamma:
    the two wrap arcs and the LE/TE bisector arc rho = R_edge. Gamma joins
    the wrap arcs' points at the span station a_s (V, W; the bisector
    crosses V-W at M) and runs straight across the gaps W -> next V;
  - outside Gamma: the column above the tip plane between the lines It -> C
    (-45 deg) and Ot -> E (0 deg), and W -> the next fin's C (+45 deg). P,
    the polygon bounding the inner zones, is the octagon through C and E.
Upstream, in the sleeve, the separators are meridional rays (they must
reach the nose core); between xTwist0 and xTwist1 (on the cylinder) each
twists from its ray into its arc.
Chordwise the hood runs xf = z0 - wrap .. xb = z3 + wrap; planes at xf, z0,
z1, z2, z3 (the base plane) and xb cross the whole fin region. Above the
tip plane the column carries the tip wall's layout out to P (the wall's LE
and TE triangles become prism columns, ending on P).
"""

from __future__ import annotations

import itertools
import math

from .polar import Pt, at, lerp, ray_polygon
from .shapes import FinSpec, bisect, polar, rotate

TWIST = (80.0, 30.0)                    # sleeve twist from xAftStart - 80 to xAftStart - 30


class FinLayout:
    def __init__(self, spec: FinSpec, P: dict):
        s = self.spec = spec
        d = self.wrap = P["finWrap"]
        self.rho_i, self.rho_o, self.Re = s.R_in - d, s.R_out + d, s.R_edge
        self.tip = s.th_tip
        z0, z1, z2, z3 = s.z
        self.xs, self.xe = P["xAftStart"], P["xSlabEnd"]
        self.xf, self.xb = z0 - d, z3 + d
        self.x_t, self.x_u = self.xs - TWIST[0], self.xs - TWIST[1]
        self.stations = [self.xf, z0, z1, z2, z3, self.xb]
        if not self.xs < self.xf or not self.xb < self.xe:
            raise ValueError("the hood (z0 - finWrap .. z3 + finWrap) must lie inside xAftStart..xSlabEnd")
        if self.x_t <= s.body.l_ogive:
            raise ValueError("the sleeve twist must lie on the cylinder")

        # The span station a_s (Gamma's vertices V, W on the wrap arcs) is where
        # the outer wrap comes back round to the angle of its own root: that
        # separator is then one straight line on the wall from the nose to the
        # base, and the body's revolve seam sits on it.
        root_o = polar(s.yz(self.rho_o, s.a_body(self.rho_o)))[1]
        a_peak = max((s.a_body(self.rho_o) + i * 0.002 for i in range(1500)),
                     key=lambda a: polar(s.yz(self.rho_o, a))[1])
        a_s = self.a_s = bisect(lambda a: polar(s.yz(self.rho_o, a))[1] - root_o, a_peak, self.tip)
        self.V, self.M, self.W = (s.yz(r, a_s) for r in (self.rho_i, self.Re, self.rho_o))
        self.It, self.Ot = s.yz(self.rho_i, self.tip), s.yz(self.rho_o, self.tip)
        self.seps = ((self.rho_i, self.V), (self.Re, self.M), (self.rho_o, self.W))   # (arc radius, Gamma point)
        # Sleeve rays through V and W, the bisector's halfway between them (it
        # meets Gamma near M and twists into its arc), and a fourth separator,
        # a meridional ray through the gap 45 deg past the bisector's. The nose
        # core's corners sit on the bisector rays: each side's other rays are
        # then symmetric about its middle (the gap ray), the grid lines through
        # the middles cross at the axis, the apex is a block corner, and the
        # body's revolve seam runs along the gap ray.
        self.ray_kinds = ("a", "m", "b", "g")
        roots = [polar(s.yz(r, s.a_body(r)))[1] for r, _ in self.seps]
        phi_v, phi_w = polar(self.V)[1], polar(self.W)[1]
        phi_m = 0.5 * (phi_v + phi_w)
        phi_g = phi_m + math.pi / 4
        self.fin_rays = [phi_v, phi_m, phi_w, phi_g]
        self.rays = [phi + s.turn(k) for k in range(s.N) for phi in self.fin_rays]
        # behind the base the separators reach the base core radially from
        # where their arcs leave r = R
        self.root_rays = [phi + s.turn(k) for k in range(s.N) for phi in roots + [phi_g]]
        for rays in (self.rays, self.root_rays):
            if any(b <= a for a, b in itertools.pairwise(rays)) or rays[-1] - rays[0] >= 2 * math.pi:
                raise ValueError("separators out of order")
        self.gamma = [rotate(p, s.turn(k)) for k in range(s.N) for p in (self.V, self.W)]
        self.G = at(ray_polygon(phi_g, self.gamma), phi_g)

        # P: per fin C (-45 deg), E (0) and D (30); the column above the tip
        # plane spans C-E
        r_p = P["rZone"]
        self.p_rays = [math.radians(a) + s.turn(k) for k in range(s.N) for a in (0.0, 30.0, 45.0)]
        self.poly_p = [at(r_p, a) for a in self.p_rays]
        self.C, self.E, self.D, self.C1 = (at(r_p, math.radians(a)) for a in (-45.0, 0.0, 30.0, 45.0))
        self.outer_lines = [(self.It, self.C), (self.Ot, self.E), (self.W, self.D), (self.G, self.C1)]   # fin 0
        self.outer_up_lines = [(self.V, self.E), (self.W, self.D), (self.G, self.C1)]

    def ray_end(self, phi: float) -> float:
        """Distance along a ray to Gamma."""
        return ray_polygon(phi, self.gamma)

    def chord_image(self, rho: float) -> Pt:
        """Where a tip-plane point at rho lands on P's edge C-E: the column
        above the tip plane carries its layout straight out to P."""
        return lerp(self.C, self.E, (rho - self.rho_i) / (self.rho_o - self.rho_i))

    def summary(self) -> list[str]:
        deg = math.degrees
        pts = (("V", self.V), ("M", self.M), ("W", self.W), ("It", self.It), ("Ot", self.Ot))
        g = self.gamma
        gap = min(math.hypot(*lerp(g[1], g[2], t / 200)) for t in range(201))
        return ["  fin layout points (fin 0, r mm / phi deg): "
                + ", ".join(f"{n} {polar(p)[0]:.1f}/{deg(polar(p)[1]):.1f}" for n, p in pts)
                + f"; Gamma's gap edges come within {gap - self.spec.R:.1f} mm of the wall",
                "  separators per fin (deg, at Gamma / at the wall): "
                + ", ".join(f"{k} {deg(a):.2f}/{deg(b):.2f}" for k, a, b in zip(self.ray_kinds, self.fin_rays,
                                                                               self.root_rays)),
                f"  sleeve twist x {self.x_t:g}..{self.x_u:g}; hood x {self.xf:g}..{self.xb:g}, planes "
                + ", ".join(f"{x:.2f}" for x in self.stations)]
