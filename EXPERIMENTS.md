# Experiment log: Ma 1.6 validation

This log tracks the validation against the wrap-around-fin report (B1F1, Ma 1.6,
α = 0). Targets are listed in `CLAUDE.md` → "Validation target". Metrics are
given as the mean over the last 100 iterations of `results/<case>/coefficients.csv`,
plus `results/<case>/wall_breakdown.csv` for the base/forebody split.

| Metric | Target | Acceptance |
|---|---|---|
| `Mx` (= −C_ℓ), fin-only roll (MDAC CX5C04) | +0.0094 | +0.006 … +0.013 (φ scatter) |
| C_Af (forebody), MDAC DX5B04 | 0.490 | report only (LaRC CT9S08: 0.400) |
| C_Ab (base, full base area), MDAC DX5B04 | 0.150 | report only (LaRC CT9S08: 0.097) |

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

**Result (2026-10-20): no effect, on the roll or on the base pressure.**
11,081,062 cells. All 40 seams stitched, including the sleeve's. checkMesh
gave no broken cells; max skewness was 6.46 (baseline 6.03). The run was
flat from it. 1000 and stopped at it. 1762, at 38.5 s/it and 44 GB.

| | s09 | baseline |
|---|---|---|
| Mx (last 100) | −0.00209 | −0.00211 |
| C_A (gauge pressure, `wall_breakdown`) | 0.610 | 0.630 |
| C_Af | 0.4325 | 0.432 |
| C_Ab, annulus force / S (area 0.873 S) | 0.177 | 0.198 |
| base Cp, mean (= C_Ab over the full base area) | −0.20 (−0.194 … −0.208) | −0.198 (−0.16 … −0.21) |
| LE / TE bevels concave + convex | +0.0133 / +0.0019 | +0.0133 / +0.0019 |

- The base pressure is unchanged. The lower C_Ab is only the area the
  sting covers. The convex TE bevel still has reversed flow at every span
  (`te_map.txt`, the same as s08's table).
- **`coefficients.csv` Cx (0.680) is wrong by +0.0706.** The `forces`
  object integrates absolute p (`pRef` 0). Without the sting's disk the
  walls are not closed, which adds p∞·πr_s²/(q∞S) = 0.0706. The viscous
  parts agree (0.0891) and Mx is unaffected (the disk is on the axis). Use
  the `wall_breakdown` C_A. The case's `forces` now has `pRef $pInf`, so
  clones of it are correct.
- A sting of d/D = 0.36 is ruled out as the cause of the C_Ab gap
  (0.20 vs 0.10–0.15) and of the roll deficit.
- `wall_breakdown.py` failed at first on stitched faces with more than 10
  vertices, which OpenFOAM writes over several lines; the parser is fixed.

### s10: bare body B1 against the report over Mach, Ma 1.6 first
Case `gm_b1_M1.6`: `gm_no_stab_M1.6`'s bare-body case (BCs, `postProcess`)
with `gm_arc_M1.6`'s `gmshParams`, i.e. the finned baseline's mesh without
the fins. Budget 6.20 M cells (the old `gm_no_stab_M1.6` is 9.6 M, with
nTheta 16 and a finer layer stack).

Why: our base pressure is C_Ab ≈ 0.20 with fins and 0.175 on the bare body,
against the report's B1F1 0.10 (LaRC) or 0.15 (MDAC), and B1 alone 0.086.
In both series the report's Ma 1.6 point is a dip: B1 alone goes 0.22 at
Ma 1.3, 0.086 at 1.6, 0.14 at 1.9. Reference: B1 alone, AEDC DT9B01,
digitized by the user in `literature/B1_MainBalance_{1,2,3}.csv`
(report PDF p.55–57).

Hypothesis: if our C_Ab and C_Af follow the tunnel at the Mach numbers
around 1.6, the 1.6 point is a tunnel artefact and the base flow needs no
numerics or turbulence-model change. If we sit about 2× high everywhere,
the base flow itself is wrong.
Plan: Ma 1.6 first, which also checks the grid against `gm_no_stab_M1.6`
(C_Ab 0.175, C_Af 0.230). Its cost decides how many more Mach numbers to
run (candidates 1.2, 1.3, 1.9, 2.36, 2.86, then transonic). Each is a
`create_case.py --template openfoam/gm_b1_M1.6 --Mach M` clone. This is the
first run with core pinning and jemalloc.

**Result, Ma 1.6 (2026-10-21):** 6,199,936 cells, 10 seams stitched, checkMesh
with one aspect-ratio flag (552 cells) and no broken cells. Flat from it. 1300
(Cx 0.4079 → 0.4091 at it. 1847), 19.9 s/it, so about 7 h to converge.

| B1 alone, Ma 1.6 | C_Af | C_Ab (= −Cp_b mean) | C_A |
|---|---|---|---|
| s10, 6.2 M | 0.228 | 0.181 | 0.409 |
| `gm_no_stab_M1.6`, 9.6 M | 0.230 | 0.175 | 0.405 |
| report B1 (AEDC DT9B01) | 0.167 | 0.086 | 0.253 |

- Grid-independent at this level (C_Af 1 %, C_Ab 3.5 %). The gap to the
  tunnel is not resolution.
- Base Cp −0.151 on the axis to −0.190 at r = 25 mm.
- Run history: hung at it. 248 on an Open MPI TCP handshake (fixed with
  `--mca btl self,sm`, c24e931, resumed with `--resume`). Then SIGKILL at
  it. 1847, almost certainly OOM: with jemalloc's decay disabled, freed
  memory never went back to the kernel (42 GB at it. 592 on a ~27 GB case).
  Reconstructed and post-processed by hand at it. 1800.

Other Mach numbers: `gm_b1_M1.3` and `gm_b1_M1.9`, the neighbours of the
report's Ma 1.6 dip, run in sequence (`create_case.py --template
openfoam/gm_b1_M1.6 --Mach M`; same mesh, same freestream p and T, so
Re ∝ U). Report B1 C_Ab / C_Af: Ma 1.3 0.220 / 0.239, Ma 1.9 0.140 / 0.096.
If ours land near those values, the 1.6 point is the tunnel's outlier.

**Result, Ma 1.3 and 1.9 (2026-10-22):** same 6.2 M mesh, 20.5 s/it. Ma 1.3
was flat from it. ~2400 (Fx +0.02 % to it. 3500) and was stopped and
post-processed at it. 3500. Ma 1.9 was stopped at it. 1100, with Fx still
creeping +0.2 N per 100 iterations (≈ 0.15 % left on C_A). Cy, Cz and the
moments are ≤ 1e-4, as the axisymmetry requires.

| B1 alone | C_Af ours / report | C_Ab ours / report | C_A ours / report |
|---|---|---|---|
| Ma 1.3 (it. 3500) | 0.248 / 0.239 | 0.222 / 0.220 | 0.470 / 0.459 (+2 %) |
| Ma 1.6 (it. 1800) | 0.228 / 0.167 | 0.181 / 0.086 | 0.409 / 0.253 (+62 %) |
| Ma 1.9 (it. 1100) | 0.213 / 0.096 | 0.163 / 0.140 | 0.376 / 0.236 (+59 %) |

Report C_A = C_Af + C_Ab (both on S_ref 12.566 in²; C_A itself is not
plotted). Our split: nose pressure 0.181 / 0.167 / 0.157 and friction
0.067 / 0.061 / 0.056 at Ma 1.3 / 1.6 / 1.9.

- Ma 1.3 matches on both base and forebody. The same mesh and numerics miss
  Ma 1.6 and 1.9 by about 60 % in C_A.
- The report's C_A drops 45 % from Ma 1.3 to 1.6, ours 13 %. Its C_Af at
  Ma 1.6 equals our nose pressure drag alone, leaving nothing for friction,
  and at Ma 1.9 (0.096) it is below our nose pressure drag (0.157). No split
  of our C_A fits the tunnel at Ma 1.6–1.9.
- The B1F1 data split the same way by facility. At Ma 1.6, LaRC CT9S08 gives
  C_Af 0.400, C_Ab 0.097, C_A 0.497, and MDAC DX5B04 gives 0.490, 0.150,
  0.640. The base dip at Ma 1.6 is only in the C/D-prefixed series (B1
  DT9B01, B1F1 CT9S08/S09); MDAC has none. Our finned C_A is 0.656
  (`finwake`), 0.630 (`gm_arc_M1.6`) and 0.610 (s09), i.e. on MDAC.
- The fin increment in C_Af agrees: LaRC 0.400 − 0.167 = 0.233, ours
  0.454 − 0.228 = 0.226.

**Conclusion: the base flow is not a numerics problem, and nothing changes in
the base modelling.** The LaRC/AEDC axial data at Ma 1.6–1.9 are
inconsistent with the same series at Ma 1.3 and with MDAC. How C_Af and C_Ab
were reduced is in the original data reports (refs 1–6, p.17), which we
don't have. The axial targets move to MDAC DX5B04 (table above), with
LaRC as a lower reference. Results: `results/gm_b1_M1.{3,6,9}/`.

### s11: first wall cell on the TE bevels, `finwake` vs Gmsh (no solve)
Hypothesis: `finwake`'s convex TE bevel stays attached because snappy drops
the layer stack at the sharp fin edges (the ~4 mm staircase band known since
the September layer fix), so its wall cell there is coarse (high y⁺), while
the Gmsh mesh keeps 3 µm and separates.
Method: `finwake` re-meshed unchanged on the new server (13,298,786 cells,
fin layers 12.1/13, 97.1 %). `scripts/first_cell.py` samples the wall cells'
centres onto the fin patch and reports h1 = 2 × the face-to-centre wall
distance. Validated on s09 (Gmsh): 2.4–4.0 µm on every fin surface, no face
above 10 µm.

**Result (2026-10-22): rejected.** `finwake`, fin 0, r = 50–80 mm, h1 in µm:

| Distance before TE | Convex bevel | Concave bevel |
|---|---|---|
| 4.75 mm (shoulder at 4.83) | 2.8–3.1 | 2.9–3.1 |
| 4.25 mm (Gmsh separates here, s07/s08) | 3.3–4.1 | 3.9–4.3 |
| 3.25 mm | 7–9 | 8 |
| 1.75 mm | 19–22 | 22 |
| 0.25 mm (TE edge) | 140–170 | 160–170 |

- Where the Gmsh flow separates, both meshes have a 3–4 µm wall cell.
- The layer collapse is real (median 33 µm, p90 147 µm over all TE
  bevels), but it builds up downstream of the separation point and is the
  same on both sides, while only the convex side stays attached on finwake.
- Still different between the meshes: chordwise cells along the bevel
  (snappy 0.2–0.4 mm; Gmsh graded at the shoulder only, growing towards
  ~1 mm), bevel-normal cells past the 0.6 mm stack (snappy 0.2–0.5 mm; Gmsh
  ratio 1.2 to the 2.25 mm cap, ~1 mm at the wrap edge), and the turbulence
  model question.

Next: s12, refine the Gmsh TE region (along and off the bevel) to snappy's
level, to test whether the Gmsh separation is grid-converged before any
turbulence-model comparison (SA, then EBRSM). In parallel, the user checks
the TE geometry in the report (bevel sides, shoulder radius, BL trip).
Results: `results/arc_M1.6_finwake/first_cell.txt`,
`results/gm_arc_M1.6_s09_sting/first_cell.txt`.

### s12: rounded fin LE, r = 0.16 mm (Gmsh)
Case `gm_arc_M1.6_s12_leround`: `gm_arc_M1.6` with `finLERadius 0.16` and
`hLERound 0.07` in `gmshParams`. Fig. 2 and Table 4 give the LE and TE
apexes a radius r = 0.008 in (0.16 mm here); the bevel shoulders are sharp
(both confirmed by the user). The apex stays at x = 660 (C_R apex to apex),
so the inscribed round moves the wedge and its shoulder 0.26 mm forward.

The TE stays sharp, as a documented deviation. With a rounded TE ending on
the base plane, the round, the base face and the rim meet tangentially at
every fin root, which hexahedra cannot fill. Shifting the fin 0.16 mm aft
avoids that, but the round then crosses the rim's shear-layer surface
r = R behind the base and breaks the wake-core blocks under each fin. The
TE needs its own re-layout of the base region; deferred until we know
whether the radius matters at all.

Geometry (`scripts/gmesh`): each LE half (the round's arc and the bevel)
is one exact quadratic NURBS edge, and the fin is its section revolved
about the arc axis. The scad's ruled wedges pair the arcs' root angles;
the revolved wedges are cones, which differ from them by < 8 um on the
exposed fin. Exposed fin volume +199 mm^3 (4 x 49.4 mm^3 in the LE slab:
the barrel is 0.26 mm longer), wall area +20 mm^2. All geometry checks
pass. The wall-area check tolerance is 5e-5 for this geometry: OCC's area
of trimmed NURBS pieces is off by up to 0.5 mm^2 (checked against direct
integration).

Mesh: 4 cells on each half of the round (first cell 0.07 mm at the apex).
The LE's chordwise count is tied, through the LE C-grid, to the radial
count at the hood front and from there to the sleeve's circumferential
count. So the finer LE also refines the body and nose rings: 12.07 M
cells (baseline 10.86; fins 5.69, body 3.59, nose 1.26). That is the only
side effect. Coarsened test mesh: near the LE the minimum scaled Jacobian
and the number of cells below 0.1 are the same as the baseline's (0.0071,
107).

Hypothesis: the round changes the LE bevel roll by a few percent (the LE
bow shock is detached either way), and the TE separation, which decides
the sign, is unaffected. Baseline `gm_arc_M1.6`: Mx −0.00204, LE +0.01319,
flat −0.01717, TE +0.00194.
Result: *pending*

### Target change (2026-10-21): compare the fin-only roll
The report's own data (user-digitized B1F1 in `literature/`) show:
- At α = 0 the 4-fin model should not depend on φ, but MDAC's total reads
  −0.0227 / −0.0287 / −0.0229 at φ = 0 / 22.5 / 45 at Ma 1.6, and LaRC reads
  −0.0163 / −0.0149 at φ = 0 / 45. The measured repeatability is ±0.003 or
  more.
- At every Mach the fin-only roll ("rolling moment due to fin", fin-balance
  sum) is 22–77 % of the total. Our smooth body carries < 2 % of the roll,
  so the CFD roll is fin roll and the total includes hardware we don't model.

The target is now the MDAC fin-only C_ℓ = −0.0094 at Ma 1.6, φ = 0
(`Mx = +0.0094`), with the φ scatter (`Mx` +0.006 … +0.013) as the pass band.
`CLAUDE.md` → Validation target lists the assumptions: S_ref 12.566 in²,
and C_NF perpendicular to the root radius line, which makes the p.18
formula equal to the exact fin roll; checked on s09.

Re-scored at Ma 1.6:

| Case | Mx (fin roll) | vs +0.0094 |
|---|---|---|
| `finwake` | +0.0086 | 92 %, **pass** |
| `finvol` | +0.0078 | 83 %, pass |
| `s02` | +0.0059 | 63 %, borderline (edge of the band) |
| `fin6` | ≈ 0 | fail |
| Gmsh `gm_arc_M1.6`, s06, s08, s09 | −0.0020 … −0.0022 | wrong sign, fail |
| Gmsh sector s13 (= full body) | −0.0022 | wrong sign, fail |
| s14 rFar 1.2 m / s15 growth 1.2 / s16 SA | −0.0007 (oscillating) / −0.0017 / −0.0004 | wrong sign, fail |

The open question is unchanged: snappy keeps the convex TE bevel attached and
lands in the band; Gmsh separates both bevels (s07, s08) and gets the wrong
sign. Which TE flow is physical decides it.

### s13–s16: one unattended queue (2026-11-03)
Prompted by Sharma & Kumar, "A Ready Reckoner of CFD for Wrap-around Fins"
(INCAS Bull. 11(2) 2019, `literature/mesh-literature.pdf`), and the
professor's suggestion to try Spalart-Allmaras. The four steps form a
ladder, each one change from the previous, run back to back by
`scripts/run_queue.sh` on the server (`endTime 2000` for all four; the Gmsh
family was read at it. 1500 so far). Template: `gm_arc_M1.6` (Mx −0.00204,
Cx 0.630, 10.86 M cells).

### s13: 90° periodic sector of the Gmsh mesh
Case `gm_arc_M1.6_s13_sector`: `caseProperties` `sectors 4`. At α = 0 the
flow is 4-fold periodic, so one fin in a 90° wedge with cyclic sides is the
same problem at a quarter of the cells. `gen_mesh.py` tags every block with
its wedge (cell zones `<zone>_w<k>`); `rebuild-mesh.sh` keeps wedge 0 (fin
0, +Y) with `topoSet`/`subsetMesh`, pairs the exposed faces as a rotational
`cyclic` (`wedgeA` z < 0, `wedgeB` z > 0), then stitches the seams inside
the wedge. `post_process.py` and `wall_breakdown.py` scale Cx, Mx and the
areas by `sectors`; Cy, Cz, My, Mz are one fin's and are no longer a check.
Layout change forced by the cut: the wedge boundary is the gap-ray plane
(52.79° from the fin) inside Γ and the strip G → C1 (to 45° on P) outside
it, and the wake and far end caps used 20 radial rays with the axis inside
a core block. They now take the sleeve's ring (16 rays, square core on the
bisector rays, rim cylinder) plus the outer zone's lines to P, leaving W
twice (one prism per fin, W–E–D). The equal-wedge-volume check in
`checks.py` confirms the classification. Full mesh 11.32 M (was 10.86;
wake 0.79 vs 0.56 M, far 0.43 vs 0.33 M), wedge ≈ 2.8 M.
Hypothesis: the sector reproduces the full body. Pass: Mx within ±0.0003
of −0.00204, Cx within ±0.003 of 0.630 (s04–s06 moved the wake cells and
changed Mx by < 0.0002, so the wake layout should not matter). This is
also the first HiSA run with a cyclic pair; a failure here fails the queue.
First launch (2026-11-02 23:15 UTC): the cut worked (2,835,124 cells, 49,775
faces per side) but createPatch could not order the cyclic faces: the edge
distributions broke ties towards the parametric start, so rotated copies
of an edge with opposite parametrisation got mirror-asymmetric nodes (2 mm
at the far-field end cap). `size_driven`/`two_sided` are now symmetric
under swapping the ends. Checked on a coarsened test mesh: every node on
a wedge boundary has its 90° image within 1e-11 mm. Interior nodes in the
fin hoods and columns still differ between fins by up to 0.4 mm (same in
every earlier full-body Gmsh run; not on the cyclic faces).
Second and third launch: the pair had to be made after the stitch (a
cyclic is re-ordered after every topology change and the first Γ seam
broke it), and then a plain `cyclic` still failed: the stitch moved 18
seam-edge points by 45 µm and split one face on one side only (5 vertices
against 4), and a cyclic needs identical faces. The pair is therefore a
`cyclicAMI` (geometric overlap, the identity up to those faces).
**Result (2026-11-03): pass.** 2,835,124 cells (a quarter of 11.32 M),
9.7 s/it on 12 ranks, 2000 its in 5.3 h, 12 GB. AMI weights 1.000 on both
sides. Mean of the last 100 its: Cx 0.6293 (p 0.5402, v 0.0892), Mx −0.00218
(p −0.00204, v −0.00014), flat from it. 1000 (±0.00002). Against
`gm_arc_M1.6` (Cx 0.630, Mx −0.00204): ΔMx −0.00014, ΔCx −0.0007, inside the
band. Region split (`wall_breakdown.csv`, ×4): nose 0.179, cylinder 0.051,
fins 0.201, base 0.198; fin roll −0.00208, body −0.00009. The sector is the
Gmsh family's baseline from here on (4× cheaper; Cy, Cz, My, Mz are one
fin's and no longer a check).

### s14: far-field radius 0.6 → 1.2 m (s01, on the sector)
Case `gm_arc_M1.6_s14_rfar`: s13 with `rFar 1200`. The review's rule of
thumb is a radial extent of several body lengths; ours is 0.75 L. s00's
wave-path argument says the nose wave leaves at x ≈ 0.62 m and a reflection
could reach the body only behind the base, and the characteristic far field
passed 6 % of the amplitude. Full mesh 12.16 M (far 1.28 M), wedge ≈ 3.0 M.
Hypothesis: no change. Pass: same band as s13.
**Result (2026-11-03): inconclusive, Mx does not settle.** 3,048,172 cells,
10.7 s/it, 2000 its in 6.0 h. Cx 0.6290 (p 0.5398, v 0.0892): unchanged.
Mx is not steady: the 100-it means wander between −0.0003 and −0.0011 from
it. 1000 on (s13: −0.00218 ± 0.00002 over the same span), and within the
last 300 its the instantaneous value sweeps from 0.0000 to −0.0018 with a
period of a few tens of iterations. The oscillation is in the fin loading
(Cx steady to 0.0003; `wall_breakdown` at it. 2000: fins −0.00073, body
−0.00009). The only input change is `rFar`, but it also doubles the far
zone's circumferential count (that family runs through the whole far zone),
so the far side of every seam with the inner groups is finer too. Neither a
pass nor a measured change: a pseudo-time oscillation of the residual
quantity. Options: resume s14 for +1000 its to see whether it decays, or
hold `hFar`'s circumferential count at s13's and move only the radius.
Decision (user, 2026-11-03): let s15 and s16 finish, then resume s14 to
it. 3000 (`endTime 3000`, `run-simulation.sh --resume`; the resumed part
lands in `postProcessing/forces/2000`, so `coefficients.csv` then covers
its. 2001–3000 and the first 2000 are in the result commit 3d6825c).
**Result at it. 3000 (2026-11-04): still oscillating, mean about −0.0007.**
100-it means over its. 2000–3000: −0.00105, −0.00027, −0.00066, −0.00110,
−0.00025; instantaneous range −0.0016 … +0.0006. Cx 0.6289. The
extension did not help: the larger domain turns s13's steady −0.00218
into an oscillation around −0.0007. The wall split at it. 3000 shows
where: the flat concave side −0.0703 (s13 −0.0713) and the flat convex
side +0.0559 (s13 +0.0548), the LE/TE bevels unchanged to 0.0001. So the
fin's flat-surface loading, not the bevels, moves with the far field, by
about 1.5 % of each side, and does not settle. Whether that is the
boundary distance or the far zone's finer seams is open (both changed).
`results/gm_arc_M1.6_s14_rfar/coefficients.csv` holds all 3000 its.

### s15: wall growth 1.4 → 1.2
Case `gm_arc_M1.6_s15_layer12`: s14 with `layerRatio 1.2`. The review's
rule for y⁺ ≈ 1 meshes is a growth ≤ 1.2; our 13-cell stack grows at 1.4
(0.59 mm thick, then 1.2). With `layerRatio = growth` the stack is one 1.2
progression from 3 µm: ~10 more cells on every wall-normal edge. Full mesh
16.04 M (fins 8.95 vs 5.51 M), wedge ≈ 4.0 M. s11 left the bevel-normal
cells past the 0.6 mm stack as one of the remaining differences to
`finwake`.
Hypothesis: a change here means the wall-normal resolution past the stack
matters for the TE-bevel separation. Pass: unchanged from s14 means the
stack is converged.
**Result (2026-11-03): small change, inside the s14 noise.** 4,019,396
cells, 14.4 s/it, 2000 its in 8.1 h. Cx 0.6272 (p 0.5379, v 0.0893; base
0.1962 vs 0.1972), Mx −0.00171 (p −0.00156, v −0.00014), 100-it means
−0.0016 … −0.0019 from it. 750, instantaneous range −0.0026 … −0.0009
(noisier than s13's ±0.00002, quieter than s14). Fin roll −0.00192. The TE
map is s13's: both bevels separate ~0.5 mm past the shoulder, convex
Cf_x −0.6 … −1.0·10⁻³. The finer wall-normal stack changes neither the
separation nor the sign; it moves Cx by −0.002 and Mx by about +0.0005
against s13, less than s14's oscillation. The stack is converged for this
question.

### s16: Spalart-Allmaras
Case `gm_arc_M1.6_s16_sa`: s15 with `RASModel SpalartAllmaras`, `0/nuTilda`
(freestream 3 ν∞ = 1.556e-4 m²/s per Spalart & Rumsey 2007, 0 at the wall,
`k`/`omega` removed), `div(phi,nuTilda)` and `grad(nuTilda)` as for k/ω,
relaxation 0.5. `run-simulation.sh` now derives the required fields from
the RAS model. SST's shear-stress limiter separates earlier under an adverse
pressure gradient than SA; the roll sign is decided by whether the convex
TE bevel separates within 1 mm of the shoulder (s07).
Hypothesis: SA keeps the convex bevel attached and Mx moves towards
`finwake`'s +0.0086. Report Mx against +0.0094 and the TE map
(`results/<case>/te_map.txt`, Cf_x sign on both bevels). A pass is not
proof that SA is right; the comparison is the result.
**Result (2026-11-04): fail. SA separates both bevels too and raises the
base drag.** 4,019,396 cells (s15's mesh), 2000 its in 7.8 h. Mx −0.00040
(p −0.00026, v −0.00014), 100-it means −0.0004 … −0.0006 from it. 1000,
instantaneous range −0.0014 … +0.0003; fin roll −0.00071. The convex TE
bevel still separates: Cf_x −0.8 … −1.4·10⁻³ over the last 2.5 mm, SA only
moves the separation ~1 mm aft (reattached-looking +0.07 … +0.17·10⁻³ at
3.75 mm before the TE where SST has −0.8). Both sides stay symmetric, so
the roll stays a small residual of the wrong sign. The axial force is
worse: Cx 0.6711 (p 0.5718, v 0.0993), base C_Ab 0.229 against SST's
0.196 and the MDAC target 0.150 (base Cp −0.247 at r = 25 mm vs −0.206),
cylinder friction 0.0573 vs 0.0511, viscous Cx +11 %. Turbulence model
is not what separates `finwake` from the Gmsh family; the finwake convex
bevel stays attached for a reason SA does not reproduce either.

### s17–s23: second sector queue, prepared 2026-11-06 (not started)
All clone `gm_arc_M1.6_s13_sector` (one change each, `endTime 2000`), to run
back to back with `scripts/run_queue.sh` in the order below (cheapest and
most informative first, the 31 M-cell full mesh last). Times assume 3.3 s/it
per million sector cells on 12 ranks. Pass bands are s13's unless stated.

### s17: Reynolds number of the MDAC tunnel
Case `gm_arc_M1.6_s17_re`: `pInf` 28536.7 → 79278.6 Pa, so Re_D 0.84e6 →
2.33e6 (MDAC 7.0e6/ft × 4 in); `omegaInf` follows (595460). Same mesh, so
y⁺ rises ×2.8 (wall average ~1.3, nose-tip maximum ~20): to be read from
`yPlus` before trusting wall-sensitive numbers; a 1 µm first cell is the
follow-up if the roll moves. September's Re test (snappy, bevel attached)
moved the roll < 2 %; on the Gmsh mesh the bevel separation decides the
sign, and the bevel sits inside the body's boundary layer, whose thickness
the Reynolds number sets. 2.84 M cells, ~5.5 h.
Hypothesis: if the separation is Re-sensitive, the roll moves by more than
0.0005; if not, Re is cleared for the whole comparison.
Result: *pending*

### s18: far-field radius alone
Case `gm_arc_M1.6_s18_rfarclean`: `rFar 1200` with `hFar 120` (scaled with
the radius, so the far zone's circumferential count stays at s13's: far
zone 0.448 M vs 0.434 M full-mesh cells, only the radial ramp grows). s14
changed both the radius and the far seams' cell size and oscillated. 2.83 M
cells, ~5.5 h.
Hypothesis: s13's steady −0.00218 returns if s14's oscillation came from
the finer far seams; it oscillates if the boundary distance itself matters.
Result: *pending*

### s22: Ma 2.5 on the sector
Case `gm_arc_M1.6_s22_M2p5`: `create_case.py --Mach 2.5` (UInf 850.66 m/s,
kInf/omegaInf recomputed). Report fin-only roll at Ma 2.5 (CX5C04)
C_ℓ −0.0056, i.e. Mx +0.0056; the report's roll magnitude has a minimum at
Ma 1.6, so a second Mach point tests whether the sign problem is a Ma 1.6
feature. Same mesh, ~5.5 h.
Hypothesis: the Gmsh family keeps the wrong sign at Ma 2.5 too (the TE
bevel separation is set by the base expansion, which is Mach-weak), or
flips, which would make Ma 1.6 the special case the report also shows.
Result: *pending*

### s23: Ma 3.0 on the sector
Case `gm_arc_M1.6_s23_M3p0`: `--Mach 3.0` (UInf 1020.79 m/s). Report
fin-only at Ma 3.0: C_ℓ −0.0192, Mx +0.0192, the largest fin-only value in
the dataset. Same mesh, ~5.5 h.
Hypothesis: as s22; a large positive target makes a sign error unmistakable.
Result: *pending*

### s19–s21: uniform refinement ladder on the sector
`gmshParams` `coarsen` (new): every size cap is multiplied by it, the 3 µm
wall stacks and the seam/LE end conditions stay, so the ladder refines
everything off the wall uniformly. Ratio 1.25 per step:
- s19 `gm_arc_M1.6_s19_k125`: coarsen 1.25, 7.91 M full → 1.98 M cells, ~4 h
- s13: coarsen 1.0, 11.32 M → 2.84 M (done)
- s20 `gm_arc_M1.6_s20_k080`: coarsen 0.8, 18.53 M → 4.63 M, ~8.5 h
- s21 `gm_arc_M1.6_s21_k063`: coarsen 0.63, 31.44 M → 7.86 M, ~14.5 h; the
  31 M full mesh goes through Gmsh, gmshToFoam and subsetMesh before the
  cut (s15's 16 M did), so this one runs last.
Hypothesis: the Gmsh family's −0.0022 is grid-converged (Richardson over
four points), or it drifts towards the finwake value with refinement,
which would make the gap a resolution question after all.
Result: *pending*

### Planned
- extend finwake from it. 1000 (restart, +500–1000 its) to converge Cx / C_Ab
- s02 redo with the sleeve stopping at x = 0.80 (removes the base-region confound)
- downstream trim (outlet 5.6 → 2.4 m, ~−0.4 M cells), then a small base refinement (level-5 disk x 0.80–0.82)
- fin-region cell-size ladder (fin boxes one level coarser), then Richardson extrapolation
- fin layer handover (16 layers, last layer ≈ 0.47 mm)
- directional y/z refinement in a fin-following sleeve
- fin TE/LE feature edges at level 7
- s01 (larger domain) as a confirmation
