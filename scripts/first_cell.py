#!/usr/bin/env python3
"""First wall-cell height on the fins, mesh only (no solution needed).

Usage: python3 scripts/first_cell.py <case> [<case> ...]   (repo root, OpenFOAM env sourced)

Writes the cell centres (writeCellCentres), samples them onto the fin patch
from the wall cells (patchInternalField), and reports h1 = 2 x the wall-normal
distance from each face to its cell centre: by surface (LE bevel, flat, TE
bevel) over all fins, then mapped over the last 8 mm of fin 0 like te_map.py.
"""
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "scripts")
from wall_breakdown import read_faces, read_vectors

SAMPLE_DICT = """FoamFile { version 2.0; format ascii; class dictionary; object firstCell; }
functions
{
    firstCell
    {
        type            surfaces;
        libs            (sampling);
        surfaceFormat   foam;
        fields          (C);
        interpolationScheme cell;
        surfaces
        {
            stabilizers
            {
                type patchInternalField; patches (stabilizers);
                offsetMode normal; distance 1e-9; interpolate false;
            }
        }
    }
}
"""

X_TE = 0.800
DX = 0.0005
NX = 16          # last 8 mm of chord
R_BINS = [0.036, 0.050, 0.065, 0.080, 0.090]
PCT = (10, 50, 90, 99)


def sample(case: Path) -> Path:
    first = min((t for t in case.iterdir() if t.name.replace(".", "").isdigit()), key=lambda t: float(t.name))
    run = ["-case", str(case), "-time", first.name]
    subprocess.run(["postProcess", "-func", "writeCellCentres", *run], check=True, stdout=subprocess.DEVNULL)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "firstCell"
            dict_path.write_text(SAMPLE_DICT)
            subprocess.run(["postProcess", "-dict", str(dict_path), "-fields", "(C)", *run],
                           check=True, stdout=subprocess.DEVNULL)
    finally:
        for name in ("C", "Cx", "Cy", "Cz"):
            (first / name).unlink(missing_ok=True)
    return case / "postProcessing" / "firstCell" / first.name / "stabilizers"


def fin_faces(d: Path):
    """Yield (face centre, unit normal into the wall, h1) per face."""
    pts = read_vectors(d / "points")
    cells = read_vectors(d / "vectorField" / "C")
    for f, cc in zip(read_faces(d / "faces"), cells, strict=True):
        sx = sy = sz = 0.0
        for a, b in zip(f, f[1:] + f[:1]):
            (x1, y1, z1), (x2, y2, z2) = pts[a], pts[b]
            sx += y1 * z2 - z1 * y2
            sy += z1 * x2 - x1 * z2
            sz += x1 * y2 - y1 * x2
        area = math.sqrt(sx * sx + sy * sy + sz * sz)
        n = (sx / area, sy / area, sz / area)
        c = tuple(sum(pts[i][k] for i in f) / len(f) for k in range(3))
        yield c, n, area, 2.0 * abs(sum((c[k] - cc[k]) * n[k] for k in range(3)))


def pct(vals, p):
    vals = sorted(vals)
    return vals[min(len(vals) - 1, int(p / 100 * len(vals)))]


def main():
    for arg in sys.argv[1:]:
        case = Path(arg).resolve()
        faces = list(fin_faces(sample(case)))
        fin_le = min(c[0] for c, _, _, _ in faces)
        print(f"\n### {case.name}  ({len(faces)} fin faces)")
        print(f"first-cell height h1 [um], all fins   {'  '.join(f'p{p:<4d}' for p in PCT)}  area > 10 um")
        groups = {}
        for c, n, area, h in faces:
            g = "LE bevel" if n[0] > 0.2 else "TE bevel" if n[0] < -0.2 else "flat"
            if g == "flat" and c[0] > X_TE - 0.005:
                g = "flat, last 5 mm"
            elif g == "flat" and c[0] < fin_le + 0.010:
                g = "flat, first 10 mm"
            groups.setdefault(g, []).append((h, area))
        for g in ("LE bevel", "flat, first 10 mm", "flat", "flat, last 5 mm", "TE bevel"):
            hs = groups.get(g, [])
            if not hs:
                continue
            frac = sum(a for h, a in hs if h > 10e-6) / sum(a for _, a in hs)
            print(f"  {g:34s} " + "  ".join(f"{pct([h for h, _ in hs], p) * 1e6:6.1f}" for p in PCT)
                  + f"   {frac * 100:5.1f} %")

        # fin 0 (+Y): bow direction from the area-weighted mean z
        fin0 = [f for f in faces if abs(math.atan2(f[0][2], f[0][1])) < math.pi / 4]
        bow = math.copysign(1.0, sum(c[2] * a for c, _, a, _ in fin0))
        acc = {}
        for c, n, area, h in fin0:
            dx = X_TE - c[0]
            if dx < 0 or dx >= NX * DX:
                continue
            side = ("convex" if -n[2] * bow > 0 else "concave") + ("_TE" if abs(n[0]) > 0.2 else "_flat")
            r = math.hypot(c[1], c[2])
            rb = next((i for i in range(len(R_BINS) - 1) if r < R_BINS[i + 1]), len(R_BINS) - 2)
            a = acc.setdefault((side, int(dx / DX), rb), [0.0, 0.0])
            a[0] += h * area
            a[1] += area
        hdr = " ".join(f"r{R_BINS[i] * 1e3:.0f}-{R_BINS[i + 1] * 1e3:.0f}" for i in range(len(R_BINS) - 1))
        for side in ("convex_TE", "concave_TE", "convex_flat", "concave_flat"):
            print(f"{side:13s} h1_um  dx_mm  {hdr}")
            for ix in range(NX):
                vals = [f"{a[0] / a[1] * 1e6:9.1f}" if (a := acc.get((side, ix, rb))) else "        ."
                        for rb in range(len(R_BINS) - 1)]
                if any(v.strip() != "." for v in vals):
                    print(f"{'':20s} {(ix + 0.5) * DX * 1e3:5.2f} {' '.join(vals)}")


if __name__ == "__main__":
    main()
