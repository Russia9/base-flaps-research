import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(sys.argv[1]); out = Path(sys.argv[2])
def last100(case):
    with open(root / "results" / case / "coefficients.csv") as fh:
        rows = list(csv.DictReader(fh))[-100:]
    return sum(float(r["Mx"]) for r in rows) / len(rows)
sweep = [(1.2, "gm_arc_M1.6_s26_M1p2"), (1.6, "gm_arc_M1.6_s13_sector"), (1.9, "gm_arc_M1.6_s24_M1p9"),
         (2.2, "gm_arc_M1.6_s25_M2p2"), (2.5, "gm_arc_M1.6_s22_M2p5"), (3.0, "gm_arc_M1.6_s23_M3p0")]
with open(root / "literature" / "B1_B1F1_digitized.csv") as fh:
    lit = [r for r in csv.DictReader(fh) if r["valid"] == "1"]
def series(ds, q):
    return sorted((float(r["Mach"]), -float(r["value"])) for r in lit if r["dataset"] == ds and r["quantity"] == q)
fin, tot = series("CX5C04", "Cl_a0_fin"), series("DX5B04", "Cl_a0")
by = {}
for ds in ("CX5C04", "CX5C05", "CX5C06"):          # fin-only at phi = 0, 22.5, 45
    for m, v in series(ds, "Cl_a0_fin"): by.setdefault(round(m, 3), []).append(v)
bm = sorted(by); blo = [min(by[m]) for m in bm]; bhi = [max(by[m]) for m in bm]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3de", "axes.edgecolor": "#c8c7c0"})
fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.axhline(0, color="#8a8984", lw=0.8, ls=(0, (3, 4)))
ax.plot(*zip(*tot), "-o", color="#1baf7a", lw=1.8, ms=5, label="experiment, total roll (MDAC)")
ax.fill_between(bm, blo, bhi, color="#eb6834", alpha=0.18, lw=0, label="fin-only roll, scatter over roll angle φ = 0, 22.5, 45°")
ax.plot(*zip(*fin), "-o", color="#eb6834", lw=1.8, ms=5, label="experiment, fin-only roll (MDAC, φ = 0)")
ax.plot([m for m, _ in sweep], [last100(c) for _, c in sweep], "-o", color="#2a78d6", lw=2.2, ms=6, label="CFD (RANS, k-ω SST)")
ax.set_xlabel("Mach number"); ax.set_ylabel("rolling moment coefficient at α = 0")
ax.set_xlim(0.7, 3.1); ax.set_xticks([0.8, 1.0, 1.2, 1.6, 1.9, 2.2, 2.5, 3.0])
ax.set_title("Wrap-around fins B1F1: induced roll vs Mach (preliminary)", loc="left", fontsize=12)
ax.legend(frameon=False, loc="lower left", fontsize=9.5)
fig.tight_layout(); fig.savefig(out, dpi=200); print("wrote", out)
