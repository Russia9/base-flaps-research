"""Our roll moment against Paek et al. 1999 (JSR 36:53): Fig. 5 (roll vs Mach, with the
JPL experiment and Edge's Navier-Stokes) and Fig. 12 (sectional roll along the chord).

Usage: plot_paek_comparison.py <repo root> <out.png> [ru]
Sign: Paek's (positive nose to tail), which is our Mx. Digitized data: literature/Koreans_Fig*.csv.
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

root = Path(sys.argv[1]); out = Path(sys.argv[2]); ru = "ru" in sys.argv[3:]
CHORD, BIN = 140.0, 10.0        # mm; wall_breakdown.py bins the fin roll in 10 mm of chord
SWEEP = [(1.2, "s26_M1p2"), (1.6, "s13_sector"), (1.9, "s24_M1p9"), (2.2, "s25_M2p2"), (2.5, "s22_M2p5"), (3.0, "s23_M3p0")]
FINE = [(1.6, "s34_mdac"), (2.5, "s35_M2p5_mdac"), (3.0, "s36_M3p0_mdac")]     # s33 mesh, MDAC freestream
EULER = [(1.6, "s37_euler")]
T = {"jpl": "эксперимент JPL", "edge": "Edge, Навье — Стокс (ламинарный)", "paek": "Paek и др., Эйлер",
     "cfd": "расчет, k-ω SST", "fine": "расчет, k-ω SST, подробная сетка", "eul": "расчет, Эйлер",
     "mach": r"$M_\infty$", "mx": r"$m_x$, $\alpha = 0$", "chord": "относительная хорда",
     "sec": r"$dm_x/d(x/b)$", "paekm": "Paek и др., Эйлер, M = {}", "ourm": "расчет, {}, M = {}",
     "euler": "Эйлер", "sst": "k-ω SST"} if ru else {
     "jpl": "JPL experiment", "edge": "Edge, Navier-Stokes (laminar)", "paek": "Paek et al., Euler",
     "cfd": "CFD, k-ω SST", "fine": "CFD, k-ω SST, finer mesh", "eul": "CFD, Euler",
     "mach": "Mach number", "mx": r"$C_l$, $\alpha = 0$", "chord": "relative chord",
     "sec": "sectional roll coefficient", "paekm": "Paek et al., Euler, M = {}", "ourm": "CFD, {}, M = {}",
     "euler": "Euler", "sst": "k-ω SST"}


def digitized(name):
    """WebPlotDigitizer export: series names in row 1, X,Y pairs per series below; sorted by X."""
    with open(root / "literature" / name) as fh:
        rows = list(csv.reader(fh))
    names = [n for n in rows[0] if n]
    data = {n: [] for n in names}
    for r in rows[2:]:
        for i, n in enumerate(names):
            if len(r) > 2 * i + 1 and r[2 * i] and r[2 * i + 1]:
                data[n].append((float(r[2 * i]), float(r[2 * i + 1])))
    return {n: sorted(v) for n, v in data.items()}


def last100(case):
    with open(root / "results" / f"gm_arc_M1.6_{case}" / "coefficients.csv") as fh:
        rows = list(csv.DictReader(fh))[-100:]
    return sum(float(r["Mx"]) for r in rows) / len(rows)


def chordwise(case):
    """Fin roll per unit relative chord, at the bin centres."""
    bins = {}
    with open(root / "results" / f"gm_arc_M1.6_{case}" / "wall_breakdown.csv") as fh:
        for r in csv.DictReader(line for line in fh if not line.startswith("#")):
            if r["group"] == "fin_chord":
                k = float(r["key"])
                bins[k] = bins.get(k, 0.0) + float(r["Mx_p"]) + float(r["Mx_v"])
    return [(k + BIN / 2) / CHORD for k in sorted(bins)], [bins[k] / (BIN / CHORD) for k in sorted(bins)]


BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#6f6e69"
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3de", "axes.edgecolor": "#c8c7c0"})
fig = plt.figure(figsize=(11, 6.4))
gs = fig.add_gridspec(2, 2, width_ratios=(1.25, 1), hspace=0.38, wspace=0.24)
ax, a1, a2 = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
for a in (ax, a1, a2):
    a.axhline(0, color="#8a8984", lw=0.8, ls=(0, (3, 4)))

f5 = digitized("Koreans_Fig5.csv")
ax.plot(*zip(*f5["Edge"]), "--^", color=GREY, lw=1.3, ms=5, mfc="white", label=T["edge"])
ax.plot(*zip(*f5["Present"]), "-o", color=ORANGE, lw=1.6, ms=4.5, mfc="white", label=T["paek"])
ax.plot(*zip(*f5["JPL Experiment"]), "s", color="#1b1b1a", ms=6.5, label=T["jpl"])
ax.plot([m for m, _ in SWEEP], [last100(c) for _, c in SWEEP], "-o", color=BLUE, lw=2.2, ms=6, label=T["cfd"])
ax.plot([m for m, _ in FINE], [last100(c) for _, c in FINE], "D", color=BLUE, ms=7, mfc="white", mew=1.8, label=T["fine"])
ax.plot([m for m, _ in EULER], [last100(c) for _, c in EULER], "*", color=GREEN, ms=14, mec="white", mew=0.8, label=T["eul"])
ax.set_xlabel(T["mach"]); ax.set_ylabel(T["mx"]); ax.set_xlim(0.5, 3.2)
ax.legend(frameon=False, fontsize=9.5, loc="upper left")

f12 = digitized("Koreans_Fig12.csv")
for a, mach, cases in ((a1, "1.5", (("s37_euler", GREEN, T["euler"], "1.6"), ("s34_mdac", BLUE, T["sst"], "1.6"))),
                       (a2, "2.5", (("s35_M2p5_mdac", BLUE, T["sst"], "2.5"),))):
    a.plot(*zip(*f12[f"M{mach}"]), "-", color=ORANGE, lw=1.5, label=T["paekm"].format(mach))
    for case, colour, model, m in cases:
        a.plot(*chordwise(case), "-o", color=colour, lw=1.6, ms=4, label=T["ourm"].format(model, m))
    a.set_xlim(0, 1); a.set_ylabel(T["sec"]); a.legend(frameon=False, fontsize=8.5)
a2.set_xlabel(T["chord"])
if ru:
    comma = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−"))
    for a in (ax, a1, a2):
        a.xaxis.set_major_formatter(comma); a.yaxis.set_major_formatter(comma)
    for a in (ax, a1, a2):
        for t in a.get_legend().get_texts():
            t.set_text(t.get_text().replace("1.5", "1,5").replace("1.6", "1,6").replace("2.5", "2,5"))
fig.savefig(out, dpi=200, bbox_inches="tight"); print("wrote", out)
