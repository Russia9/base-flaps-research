#!/usr/bin/env python3
"""Map Cp and axial skin friction near the fin TE (fin 0 only, the +Y fin).

Usage: python3 scripts/te_map.py <case> [<case> ...]   (run from the repo root, after wall_breakdown.py
has written postProcessing/wallSample/<t>/stabilizers)

Cf_x > 0: shear drags the fin downstream (attached); < 0: reversed flow.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, "scripts")
from post_process import parse_freestream
from wall_breakdown import patch_faces

X_TE = 0.800
DX = 0.0005
NX = 16          # last 8 mm of chord
R_BINS = [0.036, 0.050, 0.065, 0.080, 0.090]


def main():
    for arg in sys.argv[1:]:
        case = Path(arg).resolve()
        fs = parse_freestream(case)
        p_inf, q = fs["pInf"], fs["qInf"]
        root = case / "postProcessing" / "wallSample"
        sample = max(root.iterdir(), key=lambda t: float(t.name))
        faces = [f for f in patch_faces(sample, "stabilizers") if abs(math.atan2(f[0][2], f[0][1])) < math.pi / 4]
        # bow direction of fin 0: sign of the area-weighted mean z
        bow = math.copysign(1.0, sum(c[2] * math.sqrt(sum(si * si for si in s)) for c, s, _, _ in faces))
        acc = {}
        for c, s, pv, tv in faces:
            dx = X_TE - c[0]
            if dx < 0 or dx >= NX * DX:
                continue
            area = math.sqrt(sum(si * si for si in s))
            n_lat = -s[2] / area * bow
            side = ("convex" if n_lat > 0 else "concave") + ("_TE" if abs(s[0]) > 0.2 * area else "_flat")
            r = math.hypot(c[1], c[2])
            rb = next((i for i in range(len(R_BINS) - 1) if r < R_BINS[i + 1]), len(R_BINS) - 2)
            key = (side, int(dx / DX), rb)
            a = acc.setdefault(key, [0.0, 0.0, 0.0])
            a[0] += (pv - p_inf) / q * area
            a[1] += -tv[0] / q * area
            a[2] += area
        print(f"\n### {case.name}  (t = {sample.name}, {len(faces)} fin-0 faces)")
        hdr = " ".join(f"r{R_BINS[i]*1e3:.0f}-{R_BINS[i+1]*1e3:.0f}" for i in range(len(R_BINS) - 1))
        for side in ("convex_TE", "concave_TE", "convex_flat", "concave_flat"):
            for qty, idx, scale in (("Cp", 0, 1.0), ("Cf_x*1e3", 1, 1e3)):
                print(f"{side:13s} {qty:9s} dx_mm  {hdr}")
                for ix in range(NX):
                    vals = []
                    for rb in range(len(R_BINS) - 1):
                        a = acc.get((side, ix, rb))
                        vals.append(f"{a[idx] / a[2] * scale:+9.3f}" if a and a[2] > 0 else "        .")
                    if any(v.strip() != "." for v in vals):
                        print(f"{'':23s} {(ix + 0.5) * DX * 1e3:5.2f} {' '.join(vals)}")


if __name__ == "__main__":
    main()
