"""Mx and C_A against Mach: sector CFD vs the MDAC S-256 data (our sign, Mx = -C_l)."""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(sys.argv[1]); out = Path(sys.argv[2])

def last100(case, key):
    with open(root / "results" / case / "coefficients.csv") as fh:
        rows = list(csv.DictReader(fh))[-100:]
    return sum(float(r[key]) for r in rows) / len(rows)

def base_cab(case):
    with open(root / "results" / case / "wall_breakdown.csv") as fh:
        wb = [r for r in csv.DictReader(l for l in fh if not l.startswith("#")) if r["group"] == "region" and r["key"] == "base"]
    return float(wb[0]["Cx_p"]) + float(wb[0]["Cx_v"])

sweep = [(1.2, "gm_arc_M1.6_s26_M1p2"), (1.6, "gm_arc_M1.6_s13_sector"), (1.9, "gm_arc_M1.6_s24_M1p9"),
         (2.2, "gm_arc_M1.6_s25_M2p2"), (2.5, "gm_arc_M1.6_s22_M2p5"), (3.0, "gm_arc_M1.6_s23_M3p0")]
rematched = [(1.6, "gm_arc_M1.6_s17_re"), (2.5, "gm_arc_M1.6_s28_M2p5_re")]
for m, c in [(3.0, "gm_arc_M1.6_s29_M3p0_re")]:
    if (root / "results" / c / "coefficients.csv").exists(): rematched.append((m, c))

# report, in our sign
with open(root / "literature" / "B1_B1F1_digitized.csv") as fh:
    lit = [r for r in csv.DictReader(fh) if r["valid"] == "1" and r["config"] == "B1F1"]
def series(quantity, dataset):
    pts = sorted((float(r["Mach"]), -float(r["value"])) for r in lit if r["quantity"] == quantity and r["dataset"] == dataset)
    return pts
def band(quantity, datasets):
    by = {}
    for ds in datasets:
        for m, v in series(quantity, ds): by.setdefault(round(m, 3), []).append(v)
    ms = sorted(by); return ms, [min(by[m]) for m in ms], [max(by[m]) for m in ms]
fin0 = series("Cl_a0_fin", "CX5C04"); tot0 = series("Cl_a0", "DX5B04")
fb = band("Cl_a0_fin", ("CX5C04", "CX5C05", "CX5C06")); tb = band("Cl_a0", ("DX5B04", "DX5B05", "DX5B06"))
ca = {1.6: (0.640, 0.150), 2.5: (0.459, 0.085), 3.0: (0.403, 0.060)}   # MDAC DX5B04 C_A, C_Ab

CFD, FIN, TOT = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3de", "grid.linewidth": 0.8, "axes.edgecolor": "#c8c7c0"})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8.6), sharex=True, gridspec_kw={"height_ratios": [1.45, 1]})

# -- roll --
ax1.fill_between(tb[0], tb[1], tb[2], color=TOT, alpha=0.14, lw=0, label="report total, φ scatter 0/22.5/45° (DX5B04-06)")
ax1.fill_between(fb[0], fb[1], fb[2], color=FIN, alpha=0.16, lw=0, label="report fin-only, φ scatter 0/22.5/45° (CX5C04-06)")
ax1.plot(*zip(*tot0), "o", ms=5, color=TOT, mec="white", mew=1, label="report total, φ = 0 (DX5B04)")
ax1.plot(*zip(*fin0), "o", ms=5, color=FIN, mec="white", mew=1, label="report fin-only, φ = 0 (CX5C04)")
ax1.axhline(0, color="#8a8984", lw=0.8, ls=(0, (3, 4)))
xs = [m for m, _ in sweep]; ys = [last100(c, "Mx") for _, c in sweep]
ax1.plot(xs, ys, "-o", color=CFD, lw=2, ms=6, mec="white", mew=1, label="CFD sector, p∞ fixed (Re_D ∝ Ma)")
rx = [m for m, _ in rematched]; ry = [last100(c, "Mx") for _, c in rematched]
ax1.plot(rx, ry, "D", ms=6, mfc="white", mec=CFD, mew=1.8, label="CFD sector, Re_D 2.33e6 (MDAC)")
for x, y in zip(xs, ys):
    dx, dy = (-10, 6) if x == 1.6 else (0, -14 if y < 0.004 else 9)
    ax1.annotate(f"{y:+.4f}", (x, y), textcoords="offset points", xytext=(dx, dy), ha="right" if x == 1.6 else "center", fontsize=8, color="#52514e")
for x, y in zip(rx, ry): ax1.annotate(f"{y:+.4f} (Re 2.33e6)", (x, y), textcoords="offset points", xytext=(10, -12), ha="left", fontsize=8, color="#52514e")
ax1.annotate("not settled", (1.2, ys[0]), textcoords="offset points", xytext=(8, 10), fontsize=8, color="#52514e")
ax1.set_ylabel("Mx  (= −C_ℓ of the report)"); ax1.set_ylim(-0.010, 0.048)
ax1.set_title("B1F1 at α = 0: rolling moment against Mach, 90° sector CFD vs MDAC S-256", loc="left", fontsize=11, fontweight="semibold")
ax1.legend(loc="upper right", fontsize=8, frameon=False)

# -- axial --
cas = [last100(c, "Cx") for _, c in sweep]; cabs = [base_cab(c) for _, c in sweep]
ax2.plot(xs, cas, "-o", color=CFD, lw=2, ms=6, mec="white", mew=1, label="CFD C_A")
ax2.plot(xs, cabs, "-o", color=FIN, lw=2, ms=6, mec="white", mew=1, label="CFD base C_Ab")
ax2.plot(list(ca), [v[0] for v in ca.values()], "s", ms=6, color=TOT, mec="white", mew=1, label="MDAC C_A (DX5B04)")
ax2.plot(list(ca), [v[1] for v in ca.values()], "s", ms=6, color=FIN, mfc="white", mec=FIN, mew=1.8, label="MDAC C_Ab (DX5B04)")
for x, y in zip(xs, cas): ax2.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8, color="#52514e")
ax2.set_ylim(0, 0.78); ax2.set_ylabel("axial force coefficient"); ax2.set_xlabel("Mach number")
ax2.set_xlim(0.7, 3.1); ax2.set_xticks([0.8, 1.0, 1.2, 1.6, 1.9, 2.2, 2.5, 3.0])
ax2.set_title("Axial force: total and base", loc="left", fontsize=11, fontweight="semibold")
ax2.legend(loc="upper right", fontsize=8, frameon=False, ncol=2)
fig.text(0.01, 0.004, "CFD: HiSA, k-ω SST, 90° periodic sector of the Gmsh mesh (2.84 M cells), mean of the last 100 of 2000 pseudo-time iterations;\n"
         "Ma 1.2 ran with the 1.2 m far field and had not settled. Report values digitized from the MDAC S-256 charts; our sign is Mx = −C_ℓ.", fontsize=7.5, color="#52514e")
fig.tight_layout(rect=(0, 0.035, 1, 1))
fig.savefig(out, dpi=200)
print("wrote", out, "| sweep Mx:", [f"{m}:{y:+.4f}" for m, y in zip(xs, ys)], "| re-matched:", [f"{m}:{y:+.4f}" for m, y in zip(rx, ry)])
