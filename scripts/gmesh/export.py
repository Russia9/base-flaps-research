"""Meshing helpers shared by every mesh stage: transfinite blocks (hexahedra
and prisms) and writing a mesh as MSH 2.2 plus an optional OpenFOAM preview
case."""

from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path

import gmsh

from .occ import real_edges


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
    # -keepOrientation: gmshToFoam's own "fix" flips sound 3 um wall cells on curved walls
    for cmd in (["gmshToFoam", "-keepOrientation", "mesh.msh"], ["checkMesh", "-constant", "-noZero"]):
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
        if gmsh.model.occ.getMass(1, c) > 1e-9:        # not OCC's pole edge at the apex
            gmsh.model.mesh.setTransfiniteCurve(c, n + 1)
    set_transfinite_blocks(info)
    gmsh.model.mesh.generate(3)
    return write_case(info, out)

