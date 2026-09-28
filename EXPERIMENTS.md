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
3. **Everything is committed on the CFD server**, in two commits per step:
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

Result: *pending*

### s01: larger domain (supervisor's suggestion)
*planned*: farfield radius 0.6 → 1.2 m and outlet 5.6 → 3.2 m, with near-body
cells kept identical.

### Planned
- s02: base-region refinement (`aftBaseCylinder` 3 → 4, `baseWake` 4 → 5)
- s03: fin-region cell-size ladder (fin boxes one level coarser), then Richardson extrapolation
- s04: less dissipative reconstruction limiter (restart from baseline)
- s05: coarser fuselage surface upstream of x ≈ 0.6 m (must be neutral; frees cells)
- s06: fin layer handover (16 layers, last layer ≈ 0.47 mm)
- s07: directional y/z refinement in a fin-following sleeve
- s08: fin feature edges at level 7
