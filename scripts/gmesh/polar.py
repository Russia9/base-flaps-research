"""Cross-section layout in the (Y, Z) plane: rays, polygons, square cores.

A zone around the axis is split by rays (meridional half-planes) into
sectors; where the axis is in the fluid, a square core takes the centre.
The core's corners sit on four of the rays (a quarter of them apart), the
other rays meet its sides, and the core is cut into an m x m grid of
blocks by lines joining the ray points on opposite sides. Under the 90
degree symmetry of the rays the grid needs no parity constraint: column j
joins side 0's point j to side 2's point m - j.
"""

from __future__ import annotations

import math

Pt = tuple[float, float]


def at(r: float, phi: float) -> Pt:
    return (r * math.cos(phi), r * math.sin(phi))


def lerp(p: Pt, q: Pt, t: float) -> Pt:
    return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))


def ray_hit(phi: float, a: Pt, b: Pt) -> float | None:
    """Distance from the origin along the ray at phi to segment a-b, or None."""
    d = (math.cos(phi), math.sin(phi))
    e = (b[0] - a[0], b[1] - a[1])
    den = d[0] * e[1] - d[1] * e[0]
    if abs(den) < 1e-15:
        return None
    t = (a[0] * e[1] - a[1] * e[0]) / den          # along the ray
    s = (a[0] * d[1] - a[1] * d[0]) / den          # along the segment
    return t if t > 0 and -1e-12 <= s <= 1 + 1e-12 else None


def ray_polygon(phi: float, poly: list[Pt]) -> float:
    """Distance along the ray at phi to a closed polygon around the origin."""
    hits = [h for i in range(len(poly)) if (h := ray_hit(phi, poly[i], poly[(i + 1) % len(poly)])) is not None]
    if not hits:
        raise ValueError(f"ray at {math.degrees(phi):.2f} deg misses the polygon")
    return min(hits)


def polygon_area(poly: list[Pt]) -> float:
    return 0.5 * abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(poly, poly[1:] + poly[:1])))


class Core:
    """Square core of inradius h with corners on the rays rays[0::m], and
    the m x m block grid joining the other rays' points on opposite sides.
    rays: n angles (rad, increasing), n = 4 m."""

    def __init__(self, rays: list[float], h: float):
        if len(rays) % 4:
            raise ValueError("the rays must come in four symmetric quarters")
        self.rays, self.h, self.m = rays, h, len(rays) // 4
        self.c = rays[0]                                  # first corner's angle

    def corners(self) -> list[Pt]:
        return [at(self.h * math.sqrt(2), self.c + k * math.pi / 2) for k in range(4)]

    def hit(self, phi: float) -> float:
        """Distance along the ray at phi to the square's boundary."""
        rel = (phi - self.c) % (math.pi / 2)              # 0..90 deg from the side's first corner
        return self.h / math.cos(rel - math.pi / 4)

    def boundary_point(self, i: int) -> Pt:
        return at(self.hit(self.rays[i]), self.rays[i])

    def side(self, k: int) -> list[Pt]:
        """Side k's ray points from corner k to corner k + 1 (both corners included)."""
        n = len(self.rays)
        return [self.boundary_point((k * self.m + j) % n) for j in range(self.m + 1)]

    def grid_lines(self) -> list[tuple[Pt, Pt]]:
        """Interior grid lines: m - 1 'columns' between sides 0 and 2, m - 1
        'rows' between sides 1 and 3."""
        out = []
        for a, b in ((0, 2), (1, 3)):
            sa, sb = self.side(a), self.side(b)
            for j in range(1, self.m):
                out.append((sa[j], sb[self.m - j]))
        return out
