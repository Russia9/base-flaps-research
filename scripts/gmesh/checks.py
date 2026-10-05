"""Checks of the built model against the analytic geometry."""

from __future__ import annotations

import math
import shutil
import struct
import subprocess
from pathlib import Path

import gmsh

from .layout import zone_volume
from .occ import real_edges
from .shapes import Body, FinSpec


def check(body: Body, P: dict, info: dict) -> list[str]:
    occ = gmsh.model.occ
    L = info["layout"]
    lines, ok = [], True

    def row(label, got, want, tol_rel=1e-5):
        nonlocal ok
        err = abs(got - want) / abs(want)
        ok &= err < tol_rel
        lines.append(f"  {label:34s} {got:16.6f} {want:16.6f}  {err:.1e} {'ok' if err < tol_rel else 'FAIL'}")

    lines.append(f"  {'quantity':34s} {'OCC':>16s} {'analytic':>16s}  rel.err")
    v_body, a_lat = body.integrate(0.0, body.total)
    fins = info["fins"]
    a_wall = sum(occ.getMass(2, s) for k in ("fuselage", "stabilizers") for s in info["groups"].get(k, []))
    row("wall area (mm^2)", a_wall, fins["wall_area"] if fins else a_lat + math.pi * body.R**2)
    v_fluid = sum(occ.getMass(3, v) for v in info["zones"])
    row("fluid volume (mm^3)", v_fluid,
        math.pi * L.rff**2 * (L.xout - L.xin) - v_body - fins.get("exposed_volume", 0.0))
    for name, _, x0, x1, kind in info["zone_list"]:
        got = sum(occ.getMass(3, v) for v, n in info["zones"].items() if n == name)
        row(f"zone {name} volume (mm^3)", got, zone_volume(L, name, kind, x0, x1, fins))

    # Block topology: every block a hexahedron (6 faces of 4 edges, 12
    # edges, 8 corners) or a prism (2 triangles, 3 quads).
    for name, *_ in info["zone_list"]:
        blocks = [v for v, n in info["zones"].items() if n == name]
        bad, prisms = [], 0
        for v in blocks:
            faces = gmsh.model.getBoundary([(3, v)], oriented=False)
            edges = {c for f in faces for _, c in real_edges(f)}
            corners = {p for c in edges for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)}
            per_face = sorted(len(real_edges(f)) for f in faces)
            hexa = len(faces) == 6 and len(edges) == 12 and len(corners) == 8 and set(per_face) == {4}
            prism = len(faces) == 5 and len(edges) == 9 and len(corners) == 6 and per_face == [3, 3, 4, 4, 4]
            prisms += prism
            if not (hexa or prism):
                c = occ.getCenterOfMass(3, v)
                bad.append(f"{v}: {len(faces)}f/{len(edges)}e/{len(corners)}c {per_face} at "
                           f"({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f})")
        ok &= not bad
        lines.append(f"  zone {name:6s}: {len(blocks):4d} blocks, {len(blocks) - len(bad) - prisms} hexahedral, "
                     f"{prisms} prisms, {len(bad)} other  {'ok' if not bad else 'FAIL'}")
        for b in bad[:8]:
            lines.append(f"      not a block: {b}")

    # Block shapes: the corner angles of every planar block face, measured
    # inside the face (topology alone passes a reflex block whose cells fold).
    worst = {}
    for v, name in info["zones"].items():
        for f in gmsh.model.getBoundary([(3, v)], oriented=False):
            if gmsh.model.getType(*f) != "Plane":
                continue
            for ang, where in face_corner_angles(f[1]):
                if ang > worst.get(name, (0.0, None))[0]:
                    worst[name] = (ang, where)
    for name, (ang, where) in sorted(worst.items(), key=lambda kv: -kv[1][0]):
        good = ang < 178.0
        ok &= good
        lines.append(f"  zone {name:6s}: largest block-face corner {ang:6.1f} deg at "
                     f"({where[0]:.1f}, {where[1]:.1f}, {where[2]:.1f})  {'ok' if good else 'FAIL'}")

    worst_r, n = 0.0, 0
    for s in info["groups"]["fuselage"]:
        if gmsh.model.getType(2, s) == "Plane":
            continue
        (u0, v0), (u1, v1) = gmsh.model.getParametrizationBounds(2, s)
        uv = [c for i in range(21) for j in range(21) for c in (u0 + (u1 - u0) * i / 20, v0 + (v1 - v0) * j / 20)]
        pts = gmsh.model.getValue(2, s, uv)
        for k in range(0, len(pts), 3):
            x, y, z = pts[k:k + 3]
            worst_r = max(worst_r, abs(math.hypot(y, z) - body.r(x)))
            n += 1
    ok &= worst_r < 1e-6
    lines.append(f"  wall radius vs profile, {n} points: max |OCC - analytic| = {worst_r:.2e} mm  "
                 f"{'ok' if worst_r < 1e-6 else 'FAIL'}")
    lines.append(f"  tip half-angle {body.slope_deg(0.0):.3f} deg, shoulder kink "
                 f"{body.slope_deg(body.l_ogive - 1e-9):.3f} deg")
    lines.append("  " + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    return lines


def face_corner_angles(f: int) -> list[tuple[float, tuple]]:
    """(angle in degrees, corner) at each corner of planar face f, measured
    inside the face from the edge tangents (over 180 = reflex)."""
    edges = [c for _, c in real_edges((2, f))]
    ends = {c: [p for _, p in gmsh.model.getBoundary([(1, c)], oriented=False)] for c in edges}
    out = []
    for p in {p for c in edges for p in ends[c]}:
        cs = [c for c in edges if p in ends[c]]
        if len(cs) != 2:
            continue
        p0 = gmsh.model.getValue(0, p, [])
        dirs = []
        for c in cs:
            (t0,), (t1,) = gmsh.model.getParametrizationBounds(1, c)
            a, b = gmsh.model.getValue(1, c, [t0]), gmsh.model.getValue(1, c, [t1])
            t = t0 + (t1 - t0) * (1e-4 if math.dist(a, p0) < math.dist(b, p0) else 1 - 1e-4)
            q = gmsh.model.getValue(1, c, [t])
            d = [q[i] - p0[i] for i in range(3)]
            nd = math.sqrt(sum(x * x for x in d))
            dirs.append([x / nd for x in d])
        ang = math.degrees(math.acos(max(-1.0, min(1.0, sum(dirs[0][i] * dirs[1][i] for i in range(3))))))
        bis = [dirs[0][i] + dirs[1][i] for i in range(3)]
        nb = math.sqrt(sum(x * x for x in bis))
        if nb > 1e-9 and not gmsh.model.isInside(2, f, [p0[i] + 0.02 * bis[i] / nb for i in range(3)]):
            ang = 360.0 - ang
        out.append((ang, tuple(p0)))
    return out


def stl_mass(path: Path) -> tuple[float, float, float]:
    """Volume, area and max radius about the x axis of an ASCII or binary STL."""
    data = path.read_bytes()
    tris = []
    if data[:5] == b"solid" and b"facet" in data[:400]:
        v = [tuple(map(float, ln.split()[1:4])) for ln in data.decode().splitlines() if ln.strip().startswith("vertex")]
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


def fin_check(spec: FinSpec, fins: list, case: Path, D: float) -> list[str]:
    """Compare the OCC fins with the scad's own stabilizers.stl."""
    occ = gmsh.model.occ
    vol = sum(occ.getMass(*f) for f in fins)
    area = sum(occ.getMass(*s) for f in fins for s in gmsh.model.getBoundary([f], oriented=False))
    lines = [(f"  fins: N = {spec.N}, xi = {math.degrees(spec.xi):.2f} deg, chord {spec.L:g} mm "
              f"(x {spec.z[0]:g}-{spec.z[3]:g}), wedge run {spec.chamfer:.3f} mm, tip at r = {spec.tip_radius():.3f} mm")]
    scad = Path(__file__).resolve().parents[2] / "geometry" / "arc_stabilizers.scad"
    stl = case / "fins_scad.stl"
    openscad = next((e for e in ("openscad-nightly", "openscad") if shutil.which(e)), None)
    if not openscad:
        return lines + [f"  OCC fin volume {vol:.4f} mm^3, area {area:.4f} mm^2 (no OpenSCAD to compare)"]
    defs = f"D={D:g}; N={spec.N}; xi={math.degrees(spec.xi):.6g}; L={spec.L:g}; EXPORT=\"stabilizers\";"
    subprocess.run([openscad, "-o", str(stl), "-D", defs, str(scad)], check=True, capture_output=True)
    sv, sa, sr = stl_mass(stl)
    r_occ = max(math.hypot(*gmsh.model.getValue(0, p, [])[1:])
                for f in fins for _, p in gmsh.model.getBoundary([f], recursive=True))
    ok = True
    for label, got, want, tol in (("fin volume (mm^3)", vol, sv, 1e-3), ("fin area (mm^2)", area, sa, 1e-3),
                                  ("max radius (mm)", r_occ, sr, 1e-4)):
        err = abs(got - want) / want
        ok &= err < tol
        lines.append(f"  {label:20s} OCC {got:14.4f}  scad STL {want:14.4f}  rel.err {err:.1e} "
                     f"{'ok' if err < tol else 'FAIL'}")
    return lines + ["  " + ("FIN CHECKS PASSED" if ok else "FIN CHECKS FAILED")]
