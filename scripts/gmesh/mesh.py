"""The volume mesh: every block edge meshed with its spacing, then the
transfinite faces and blocks."""

from __future__ import annotations

import math
from pathlib import Path

import gmsh

from .export import set_transfinite_blocks, write_case


def volume(info: dict, gs: dict, out: Path, convert: bool) -> list[str]:
    topo, sizes = gs["topo"], gs["sizes"]
    for c, sz in sizes.items():
        gmsh.model.mesh.setTransfiniteCurve(c, len(sz) + 1)
    set_transfinite_blocks(info)
    gmsh.model.mesh.generate(1)
    for c, sz in sizes.items():                    # move every curve's nodes to its sizes
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
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.model.mesh.generate(3)
    lines = []
    worst: dict[str, tuple] = {}
    low: dict[str, int] = {}
    for v, name in info["zones"].items():
        for _, tags, _ in zip(*gmsh.model.mesh.getElements(3, v)):
            q = gmsh.model.mesh.getElementQualities(tags, "minSJ")
            low[name] = low.get(name, 0) + sum(1 for x in q if x < 0.1)
            i = min(range(len(q)), key=q.__getitem__)
            if q[i] < worst.get(name, (2.0,))[0]:
                c = gmsh.model.occ.getCenterOfMass(3, v)
                worst[name] = (q[i], c)
    for name, (q, c) in sorted(worst.items(), key=lambda kv: kv[1][0]):
        lines.append(f"  zone {name:6s}: minSJ {q:7.4f} (block at x {c[0]:.1f}, r {math.hypot(c[1], c[2]):.1f}), "
                     f"{low[name]} cells below 0.1  {'ok' if q > 0 else 'FAIL'}")
    lines.append("  " + write_case(info, out, convert))
    return lines
