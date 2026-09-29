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

### Planned
- s03: base-region refinement (`aftBaseCylinder` 3 → 4, `baseWake` 4 → 5), on the s02 budget
- s04: fin-region cell-size ladder (fin boxes one level coarser), then Richardson extrapolation
- s05: less dissipative reconstruction limiter (restart from baseline)
- s06: fin layer handover (16 layers, last layer ≈ 0.47 mm)
- s07: directional y/z refinement in a fin-following sleeve
- s08: fin TE/LE feature edges at level 7
- s09: s01 (larger domain) as a confirmation
