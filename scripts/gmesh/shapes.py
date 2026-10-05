"""The body (report B1) and the arc fins (report F1), analytically."""

from __future__ import annotations

import math


def bisect(f, lo: float, hi: float, n: int = 200) -> float:
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rotate(p: tuple[float, float], t: float) -> tuple[float, float]:
    c, s = math.cos(t), math.sin(t)
    return (p[0] * c - p[1] * s, p[0] * s + p[1] * c)


def polar(p: tuple[float, float]) -> tuple[float, float]:
    return math.hypot(*p), math.atan2(p[1], p[0])


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

    def integrate(self, x0: float, x1: float, n: int = 20000) -> tuple[float, float]:
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


class FinSpec:
    """Arc fin exactly as geometry/arc_stabilizers.scad builds it (report F1).

    Radii are absolute (Table 4: R3/R1/R2); the scad asserts R_out == R.
    The fin frame of fin 0 is polar about its arc centre (Yc, Zc): a point
    at radius rho and angle a is (Y, Z) = (Yc + rho sin a, Zc + rho cos a);
    the scad's mirror makes the bow bulge towards +Z, and fin k is fin 0
    rotated by k * 360 / N about +X.
    """

    R_in, R_edge, R_out = 36.0, 38.0, 40.0
    root_embed, wedge = 2.0, 45.0           # wedge: LE/TE included angle (deg)

    def __init__(self, body: Body, N: int, xi_deg: float, L: float):
        if abs(self.R_out - body.R) > 1e-9:
            raise ValueError("R_out must equal the body radius (scad assert)")
        self.body, self.N, self.L = body, N, L
        self.R = body.R
        self.xi = math.radians(xi_deg)
        self.root_y = body.R - self.root_embed
        self.Yc = self.root_y + self.R_edge * math.sin(self.xi / 2)
        self.Zc = -self.R_edge * math.cos(self.xi / 2)
        self.chamfer = (self.R_out - self.R_in) / 2 / math.tan(math.radians(self.wedge) / 2)
        self.z = [body.total - L, body.total - L + self.chamfer, body.total - self.chamfer, body.total]
        self.th = {r: math.atan2(self.root_y - self.Yc, self.Zc_at(r) - self.Zc)
                   for r in (self.R_in, self.R_edge, self.R_out)}
        self.th_tip = self.th[self.R_edge] + self.xi

    def Zc_at(self, r: float) -> float:
        """Z where the circle of radius r about the centre meets Y = root_y (root side)."""
        return self.Zc + math.sqrt(r * r - (self.root_y - self.Yc) ** 2)

    def turn(self, k: int) -> float:
        return k * 2 * math.pi / self.N

    def yz(self, rho: float, a: float, k: int = 0) -> tuple[float, float]:
        return rotate((self.Yc + rho * math.sin(a), self.Zc + rho * math.cos(a)), self.turn(k))

    def point(self, rho: float, a: float, x: float, k: int = 0) -> tuple[float, float, float]:
        return (x, *self.yz(rho, a, k))

    def centre(self, k: int = 0) -> tuple[float, float]:
        return rotate((self.Yc, self.Zc), self.turn(k))

    def frame(self, y: float, z: float) -> tuple[int, float, float]:
        """(fin k, rho, a) of a cross-section point in its nearest fin's frame."""
        k = round(math.atan2(z, y) / self.turn(1)) % self.N
        y, z = rotate((y, z), -self.turn(k))
        return k, math.hypot(y - self.Yc, z - self.Zc), math.atan2(y - self.Yc, z - self.Zc)

    def a_body(self, rho: float) -> float:
        """Fin-frame angle where the arc of radius rho leaves the body (r = R) on
        its way to the tip: walked down from the tip, since an arc can dip into
        the body and out again further round."""
        step = math.radians(0.25)
        a = self.th_tip
        while math.hypot(*self.yz(rho, a)) >= self.R:
            a -= step
            if a < self.th_tip - math.pi:
                raise ValueError(f"arc {rho} never enters the body")
        return bisect(lambda t: math.hypot(*self.yz(rho, t)) - self.R, a, a + step)

    def a_at_radius(self, rho: float, r: float) -> float:
        """First fin-frame angle past the body where the arc rho reaches radius r."""
        a0 = self.a_body(rho)
        a = a0
        step = math.radians(0.25)
        while math.hypot(*self.yz(rho, a)) < r:
            a += step
            if a > self.th_tip + math.pi / 2:
                raise ValueError(f"arc {rho} never reaches r = {r}")
        return bisect(lambda t: math.hypot(*self.yz(rho, t)) - r, a - step, a)

    def tip_radius(self) -> float:
        return math.hypot(*self.yz(self.R_edge, self.th_tip))
