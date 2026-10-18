# Experiment log: Ma 1.6 validation

This log tracks the validation against the wrap-around-fin report (B1F1, Ma 1.6,
α = 0). Targets are listed in `CLAUDE.md` → "Validation target". Metrics are
given as the mean over the last 100 iterations of `results/<case>/coefficients.csv`,
plus `results/<case>/wall_breakdown.csv` for the base/forebody split.

| Metric | Target | Acceptance |
|---|---|---|
| `Mx` (= −C_ℓ) | +0.016 | ≥ +0.0128 (80 %) |
| C_Af (forebody) | 0.40 | report only |
| C_Ab (base, full base area) | 0.10 | report only |

## Policy

1. **One change per step.** Each step changes exactly one thing relative to the
   current baseline case. If a result is neutral or better, that case becomes
   the new baseline; if not, the baseline stays.
2. **Cases go through the normal pipeline.** First clone the baseline with
   `python3 scripts/create_case.py --template openfoam/<baseline> --case openfoam/<new> --N 4 --xi 79.61 --L 140 --Mach 1.6`.
   Then edit dictionaries, run `./rebuild-mesh.sh openfoam/<new>` and
   `./run-simulation.sh openfoam/<new>`. No hand-copied case directories and no
   hand-edited meshes.
3. **Everything is committed locally and pushed to `origin`**, in two commits
   per step. The CFD server never commits. It only syncs with
   `git fetch && git reset --hard origin/master`. That is safe because the
   server holds no edits to tracked files, and result files copied back to it
   are identical. It also keeps untracked cases, meshes and runs. New case
   dictionaries are created locally with `create_case.py` and reach the
   server by that sync. `results/<case>/` is copied back from the server before
   the result commit. The two commits are:
   - *Setup:* `test(sNN): <what changes>`, containing the new case's
     dictionaries (tracked by `.gitignore`: `0/`, `system/`,
     `constant/*Properties`, `case.foam`) and this log's entry with the
     hypothesis. Commit this before meshing.
   - *Result:* `test(sNN): results — <one-line outcome>`, containing
     `results/<case>/` and the filled-in log entry.
4. **Naming:** `openfoam/arc_M1.6_sNN_<slug>`, where `NN` is the step number.
5. Solves run one at a time. A ~13 M-cell case fills the server's memory.

## Reference cases (before this log)

These cases existed before the policy. Their dictionaries were committed as
found, and only their results are logged.

| Case | Change | Mcells | Mx | % | Cx | Notes |
|---|---|---|---|---|---|---|
| `arc_M1.6_buf0` | wall-layer fix, `nBufferCellsNoExtrude 0` | 11.04 | +0.00292 | 18 | 0.670 | first clean Ma 1.6 run |
| `arc_M1.6_finvol` | + per-fin boxes (tight level 5 at 4 mm, wide level 4 at 12 mm) | 12.66 | +0.00784 | 49 | 0.668 | gain at the fin TE / base |
| `arc_M1.6_finwake` | + boxes to x = 0.83/0.89, base-wake cylinder level 4 | 13.30 | +0.00852 | 53 | 0.654 | **baseline**. C_Af 0.454, C_Ab 0.202 |
| `arc_M1.6_fin6` | finwake with fins (6 7), fuselage (4 5), fin layers 10 @ 1.6 | 12.87 | ≈ 0 | 0 | 0.648 | three changes at once; unexplained, stopped at it. 968 |
| `gm_arc_M1.6` | wall-conformal zoned Gmsh mesh (no snappy), stopped at it. 1754 | 10.86 | −0.00211 | −13 | 0.630 | **Gmsh baseline**. C_Af 0.432, C_Ab 0.198; TE bevels +0.0019 vs finwake +0.0106 |

## Steps

### s00: diagnostics on `finwake` (no solve)
Hypothesis: none. This step locates the roll deficit and the high base drag
using `scripts/wall_breakdown.py`, checks the outer boundary for reflected
waves, and compares fin pressure between `fin6` and `finwake`.

**Result (2026-09-29).** Data: `results/arc_M1.6_{finwake,fin6}/wall_breakdown.csv`.
`fin6` was sampled at it. 900.

Fin roll (`Mx`) split into three pieces, each a small net of large opposing
concave and convex terms:

| | finwake | fin6 |
|---|---|---|
| LE bevels (concave + convex) | +0.1275 − 0.1135 = **+0.0140** | +0.0135 |
| flat surfaces | −0.0806 + 0.0647 = **−0.0159** | −0.0155 |
| TE bevels | −0.0409 + 0.0515 = **+0.0106** | **+0.0022** |
| fins total | +0.0086 | +0.0001 |

- **`fin6` is explained by the trailing edge.** Its LE and flat-surface roll
  are within 4 % of `finwake`. The whole collapse is in the TE bevels, which
  sit exactly in the base plane inside the base expansion. Earlier, the
  `finvol` gain also came from the TE. The TE bevel term is the one that swings
  with the mesh, and not monotonically.
- **Base:** C_Ab = 0.202 in both cases; the fin-surface mesh doesn't move it.
  Base Cp goes from −0.156 on the axis to −0.21 at the rim. The target is
  C_Ab 0.10, i.e. mean Cp −0.10.
- **Forebody:** C_Af = 0.454 (nose 0.181, cylinder friction 0.060,
  fins 0.213) against a target of 0.40.
- **Farfield (p on lines at r = 0.1 to 0.59 m, x up to 3 m):** the nose wave
  leaves through the farfield at x ≈ 0.62 with +6 % amplitude. An inward
  reflection would reach r = 0.1 only at x ≈ 1.25. That is downstream of the
  wake recompression (x ≈ 1.0–1.05), so it cannot reach the fins or the base
  recirculation. A larger domain is not expected to change Mx or C_Ab.

### s01: larger domain (supervisor's suggestion), deferred to the end
Planned: farfield radius 0.6 → 1.2 m and outlet 5.6 → 3.2 m, with near-body
cells kept identical. s00 predicts no effect (reflections reach the axis only
downstream of the wake recompression), so it runs last as a confirmation.

### s02: coarser fuselage surface upstream of x = 0.6 m (enabler)
Case `arc_M1.6_s02_fuscoarse`, baseline `arc_M1.6_finwake`. One change in
`snappyHexMeshDict`: `fuselage` surface level (5 6) → (4 5), plus an
`aftBodySleeve` cylinder (x 0.60–0.81, r ≤ 50 mm) at level 5, so the fin
region, base rim and base disk keep their old surface level.

Hypothesis: neutral. `Mx` and C_Ab stay within 2 % while about 3 M fuselage
layer cells are freed, making room for the base refinement (s03). This is
supported by s00: `fin6` coarsened the whole fuselage and its LE and
flat-surface roll did not move.
Pass: |ΔMx| ≤ 2 %, |ΔC_Ab| ≤ 2 %; then s02 becomes the baseline.

**Result (2026-09-29): FAIL, not neutral. Baseline stays `finwake`.**
10.24 M cells, 3.06 M fewer than finwake. Fuselage faces went from 446 k to
198 k, and fin layer coverage is 98 %. `checkMesh` reported **Mesh OK**, the
first time for this mesh family. The solve ran 1000 its.

| | finwake | s02 | Δ |
|---|---|---|---|
| Mx (last 100) | +0.00852 | +0.00588 | −31 % |
| fin LE / flat / TE | +0.0140 / −0.0159 / +0.0106 | +0.0128 / −0.0149 / +0.0078 | TE −26 % |
| C_Ab | 0.202 | 0.199 | −1.6 % |
| C_Af | 0.454 | 0.438 | −3.5 % |
| friction nose / cylinder | 0.0129 / 0.0602 | 0.0074 / 0.0475 | −43 % / −21 % |
| Cz, My (should be 0) | 0.0000, 0.0000 | −0.0014, +0.014 | asymmetric |

- **Not converged:** Cx is still drifting (0.614 at it. 500, 0.637 at it. 1000),
  and My peaked at 0.018 around it. 700. The arm My/Cz ≈ 10 D places the
  spurious lateral force at the base, so the base wake went asymmetric.
- **Friction is surface-level sensitive.** Level 4 on the forebody cuts skin
  friction by 20–40 %, and `fin6` shows the same drop (nose 0.0074, cylinder
  0.044). The forebody boundary layer is not resolved at level 4, and may not
  be converged at level 5 either. With a 0.58 mm layer stack, most of it sits
  in the volume cells.
- **Confound (my setup error):** `aftBodySleeve` runs to x = 0.81, so it also
  refined the first 10 mm behind the base from level 4 (`baseWake`) to level 5.
  That puts a 5 → 4 level jump inside the base recirculation. The asymmetry and
  the TE change cannot be attributed to the forebody change alone.

### s03: less dissipative reconstruction limiter
Case `arc_M1.6_s03_wmuscl`, baseline `arc_M1.6_finwake`. One change in
`fvSchemes`: `reconstruct(rho|U|T)` goes from `wVanLeer gradTVD` to
`wMUSCL gradTVD`. That is HiSA's weighted monotonized-central limiter,
`max(min(r/cd, (1+r)/2, 1/cd), 0)`, which clips less than van Leer. Same mesh
dictionaries, so the mesh must reproduce finwake's 13,304,176 cells. It runs
from scratch for 1000 its, like the baseline. It is not a restart, because
finwake's Cx is still drifting at it. 1000 and a restart would mix the limiter
with extra convergence.

Hypothesis: the roll lives in ±0.1 concave/convex terms at the 45° LE/TE
bevels (s00), where the flow turns 22.5° through a detached LE shock and a TE
expansion into the base. van Leer's clipping at those extrema smears both. A
less dissipative limiter should sharpen them and move the LE and TE pieces.
The direction is not predicted.
Pass: none. This is a sensitivity test. |ΔMx| > 10 % means the limiter
matters and the numerics join the list of error sources.

Result: *pending*

### s04: finer axial cells in the near wake behind the TE (Gmsh mesh)
Case `gm_arc_M1.6_s04_tewake`, baseline `gm_arc_M1.6` (wall-conformal zoned
Gmsh mesh, 10.86 M cells; at it. 1500 Mx −0.0021, Cx 0.630, C_Af 0.432,
C_Ab 0.198). One change in `system/gmshParams`: `hRelax` 2.0 → 0.5 mm. That
caps the axial cells from the base plane to `xSlabEnd` (x 800–830, across the
whole slab) and sets the start of the wake zone's axial ramp at x = 830
(0.5 mm there, `hWake` 4 mm at x = 1200). The bevel itself is already at
0.2–0.31 mm and does not change. Budget: 13.07 M cells (fins 7.58 M, wake
0.56 M), close to the 13.3 M that HiSA has survived on 12 ranks.

Hypothesis: the fin LE and flat-surface roll agree across finwake, fin6 and
the Gmsh mesh, but the TE bevels give +0.0106 in finwake against about
+0.002 in fin6 and Gmsh. finwake's fin boxes reach x = 830/890 at level 5/4
(0.42/0.83 mm axially); the Gmsh near wake has 2 mm cells from x = 806. The
TE bevels sit in the base expansion, so a coarse near wake may smear the
expansion and the base pressure that load them.
Pass: TE-bevel roll moves towards finwake's +0.0106 (|Δ| > 0.002), or C_Ab
moves. If neither moves, the near-wake axial resolution is ruled out.

Result (2026-10-06): **not run.** The mesh built cleanly (13,069,674
hexahedra, the same checkMesh flags as the baseline), but HiSA was killed by
the kernel's out-of-memory killer at it. 52 (rank 7 at 4.5 GB). The baseline
peaks at about 46 GB for 10.86 M cells, about 4.3 GB per million, so about
12 M cells is the limit on 12 ranks with this mesh. Replaced by s06 (0.7 mm).

### s05: near-wake axial cell 1.0 mm (ladder point between baseline and s04)
Case `gm_arc_M1.6_s05_tewake1`, baseline `gm_arc_M1.6`. One change in
`system/gmshParams`: `hRelax` 2.0 → 1.0 mm (s04 has 0.5). Budget 11.58 M
cells (fins 6.15 M, wake 0.50 M). It runs on a second machine at the same time
as s04 (README step 5b).

Hypothesis: as s04. With the baseline (2.0) and s04 (0.5), this gives three
points at a constant refinement ratio of 2. Those give the order of
convergence and an extrapolated value for the TE-bevel roll and C_Ab.
Pass: none on its own. It feeds the ladder.

Result: *pending*

### s06: near-wake axial cell 0.7 mm (s04 within the memory limit)
Case `gm_arc_M1.6_s06_tewake07`, baseline `gm_arc_M1.6`. One change in
`system/gmshParams`: `hRelax` 2.0 → 0.7 mm. Budget 12.20 M cells, an
estimated 53 GB peak on the CFD server, which is borderline. With s05 (1.0)
it gives the ladder 2.0 / 1.0 / 0.7.

Hypothesis and pass: as s04.

Result (2026-10-08): **no effect.** Stopped at it. 4391, flat since about
it. 1500. Mean of the last 100 its: Cx 0.6310, Mx −0.00217 (baseline 0.6301,
−0.00211). C_Af 0.432, C_Ab 0.199 (baseline 0.432, 0.198). Fin roll by chord
is identical to the baseline to four decimals: LE +0.0132, flats −0.0172,
TE bevels +0.0019. Base Cp −0.164 at the axis to −0.209 at r = 25 mm, as before.
The axial cell behind the TE is ruled out, both for the roll and for C_Ab.
s05 (1.0 mm) would add nothing and is not needed. The remaining difference
from finwake is the cross-section resolution of the fin region (finwake
0.94 mm at level 5, Gmsh 2.25 mm in the hood and 3 mm outside it).

### s07: TE-bevel flow on `gm_arc_M1.6` and `finwake` (no solve)
Hypothesis: the TE-bevel viscous Cx is negative on the Gmsh mesh (−0.0001
on both bevels) and positive on finwake, so the Gmsh flow separates on the
bevels. `scripts/te_map.py` maps Cp and axial skin friction on fin 0 over the
last 8 mm of chord (0.5 mm bins), in four span bands. It reads the
`wallSample` written by `wall_breakdown.py`. The bevel runs 4.83 mm, from
the shoulder at dx ≈ 4.8 mm to the TE at x = 800.

**Result (2026-10-11): confirmed. The Gmsh bevels separate on both sides;
on finwake only the concave one does.**

| | Gmsh convex | Gmsh concave | finwake convex | finwake concave |
|---|---|---|---|---|
| Cf_x on the bevel (×10³) | −0.5 … −1.4 | −0.4 … −1.5 | +0.2 … +4.0 | −0.1 … −1.4 (r > 50 mm) |
| reversed from | ~0.5 mm past the shoulder, all spans | same | root corner only (r < 50, last 2 mm) | ~1 mm past the shoulder |
| Cp on the bevel | −0.36 → −0.24 | −0.36 → −0.23 | −0.37 … −0.38, held to the TE | −0.37 → −0.19 |

- On Gmsh the expansion at the shoulder decays straight back to the base
  pressure (Cp ≈ −0.24) and is the same on both sides, so the TE bevels
  cancel (+0.0019). On finwake the convex bevel stays attached and holds
  Cp ≈ −0.37 over the outer three span bands, while the concave bevel
  separates. That asymmetry is finwake's +0.0106.
- The Gmsh shoulder also lacks finwake's skin-friction peak (Cf_x ≈ +4·10⁻³
  over the first 1 mm of the bevel and the last 1 mm of the flat), which
  suggests the expansion fan at the shoulder is under-resolved.
- So the roll gap is a separation question on the convex TE bevel, decided
  in the first ~1 mm past the shoulder. The cells normal to the bevel and
  across the shoulder are the target of the s09+ ladder. Which mesh is
  physically right is not known yet.

### s08: chordwise grading at the TE bevel's shoulder (Gmsh mesh)
Case `gm_arc_M1.6_s08_teshoulder`, baseline `gm_arc_M1.6`. One change in
`system/gmshParams`: the new key `hTEShoulder` 0 → 0.2 mm. Axial edges
leaving the fin's TE shoulder (x = z2 = 795.17) now start at 0.2 mm and
grow at 1.2. Before, the shoulder was a free end: on the bevel the cells
grew from 0.2 mm at the TE to about 1 mm at the shoulder, and on the flat
ahead of it they were at the 2.5 mm `hFinChord` cap. The bevel goes from
10 to 14 chordwise cells. That count is coupled to a sleeve circumferential
family through the tip column, so the nose and body zones grow as well.
Budget: 12.02 M cells (baseline 10.86), below s06's 12.20 M, which ran.

Hypothesis: s07 places the separation within ~0.5 mm past the shoulder
on both bevels. The Gmsh shoulder also lacks finwake's skin-friction peak,
so the expansion fan there is under-resolved. The cells normal to the
bevel are not the cause: inside the 6 mm wrap they stay well below the
`hFinSpan` cap. Resolving the corner chordwise should keep the convex
bevel attached and move the TE-bevel roll towards finwake's +0.0106.
Pass: TE-bevel roll moves by more than +0.002, or the convex bevel's
Cf_x turns positive (`scripts/te_map.py`). If neither happens, the
shoulder's chordwise resolution is ruled out.

**Result (2026-10-13): no effect. Ruled out.** 12,017,570 cells, with
the same checkMesh flags as the baseline. The run was flat from it. 1300
and stopped at it. 2768; the peak memory was 50 GB.

| | s08 | baseline |
|---|---|---|
| Mx (last 100) | −0.00210 | −0.00211 |
| Cx | 0.6297 | 0.6301 |
| C_Af / C_Ab | 0.431 / 0.199 | 0.432 / 0.198 |
| TE bevels concave + convex | −0.0389 + 0.0406 = **+0.0017** | +0.0019 |

- The grading does resolve the shoulder. Cf_x now peaks at +3…+6·10⁻³
  over the last 1 mm of the flat and the first 0.5 mm of the bevel, as on
  finwake (`te_map.txt`).
- Both bevels still reverse from about 0.75 mm past the shoulder, at every
  span, as in the baseline. The separation is not caused by the
  under-resolved expansion fan.
- Remaining differences from finwake at the TE: the cells normal to the
  bevel past the layer stack (Gmsh grows at 1.2 to ~1 mm inside the
  wrap), the wake behind the TE across the span, and the numerics
  (s03 limiter, never run). Gmsh separating symmetrically while finwake
  stays attached on the convex side could also be a snappy artefact. That
  is not ruled out.

### s09: support sting behind the base (Gmsh mesh)
Case `gm_arc_M1.6_s09_sting`, baseline `gm_arc_M1.6`. One change: a sting
of diameter 1.12 in (28.448 mm) from the base to x = 1150, i.e. 4.4 D behind
the base and past the wake recompression at x ≈ 1000–1050. The size is
AEDC's sting 4S, "S-1.12F-35.38-1.850M", read as 1.12 in at the model end
(report PDF p.15, printed p.11). The report gives no sting size for LaRC
or MDAC.

In `gmshParams`: `stingDiameter 28.448`, `xStingEnd 1150`. The sting is the
wall patch `sting`, in its own sleeve zone that replaces the axis cores
behind the base. The `forces` function object still integrates only
`fuselage` and `stabilizers`, so the sting is not in the coefficients. The
base becomes an annulus, 87.4 % of the base area. Budget: 11.08 M cells
(baseline 10.86).

Why: the base pressure was measured on the sting (PDF p.15). The report's
C_Af = C_A − C_Ab removes the base load. C_Ab and the roll cannot be
corrected for the sting, so they include it. Our sting-free C_Ab is 0.20;
the tunnel gives 0.10 (LaRC) and 0.15 (MDAC). The TE bevels, our volatile
roll term (s07, s08), sit in that base flow.

Hypothesis: the sting raises the base pressure (C_Ab falls towards
0.10–0.15). That weakens the expansion behind the TE and changes the
TE-bevel separation and loading. The direction of the roll change is not
predicted.
Pass: none; this is a sensitivity test. If |ΔMx| > 0.002, or the TE-bevel
term changes, the sting matters for the validation. C_Ab is reported two
ways: the annulus force over S, and −Cp_b over the full base area.

Result: *pending*

### Planned
- 90° periodic sector of the Gmsh mesh (4× fewer cells): deferred. Every
  zone needs one shared side surface at the seam planes x0 and xSlabEnd,
  so the wake and far zones would need the slab's cross-section layout
- extend finwake from it. 1000 (restart, +500–1000 its) to converge Cx / C_Ab
- s02 redo with the sleeve stopping at x = 0.80 (removes the base-region confound)
- downstream trim (outlet 5.6 → 2.4 m, ~−0.4 M cells), then a small base refinement (level-5 disk x 0.80–0.82)
- fin-region cell-size ladder (fin boxes one level coarser), then Richardson extrapolation
- fin layer handover (16 layers, last layer ≈ 0.47 mm)
- directional y/z refinement in a fin-following sleeve
- fin TE/LE feature edges at level 7
- s01 (larger domain) as a confirmation
