# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Parametric CFD study of arc-shaped aft-base fins on a supersonic ogive-cylinder
fuselage. Pipeline is
OpenSCAD → STL → blockMesh → snappyHexMesh → HiSA → Python post-processor.
See `README.md` for the parameter space and physical setup.

Work is currently on two hand-tuned tracked cases at Ma 1.5, not an automated
sweep: `openfoam/arc` (N = 4, ξ = 79.61°, L = 140 mm) and `openfoam/arc_no_stab`
(N = 0, the clean-body baseline — `arc` with every stabilizer reference removed).
The 108-case `(N, ξ, L/D, Ma)` factorial is the eventual target but no script
drives it; `scripts/sweep.py` and `openfoam/template/` were removed.

The active work is **validating the induced roll at Ma 1.6** against the report
(see "Validation target" below), using untracked `openfoam/arc_M1.6_*` clones
on the CFD server. The tracked templates lag behind those clones:
- `openfoam/arc` still has `nBufferCellsNoExtrude -1`, which blows up HiSA at
  the base-rim corner at Ma 1.6. The working value is `0`.
- `openfoam/arc` has none of the per-fin refinement boxes (`finTight*`,
  `finWide*`, `baseWake`) that produced the best roll so far
  (`arc_M1.6_finwake`, 13.3 M cells).

Clone from the server case, not from the tracked template, when continuing that
work.

## Experiment tracking (policy)

Validation experiments follow `EXPERIMENTS.md`. The short version:
- Each step changes **one thing** relative to the current best case.
- New cases are created **only** with `scripts/create_case.py --template <baseline>`,
  then built with `rebuild-mesh.sh` and solved with `run-simulation.sh`.
- Each step gets two commits **on the CFD server**. The setup commit holds the
  case dictionaries and the log entry, and is made before meshing. The result
  commit holds `results/<case>/` and the outcome.
- `.gitignore` tracks every case's `0/`, `system/`, `constant/*Properties` and
  `case.foam`, so each experiment can be reproduced from git history. Meshes,
  `processor*/`, time directories, logs and `postProcessing/` stay untracked.
- `create_case.py` copies only those dictionaries, so cloning a solved
  multi-GB case is cheap.
- Do not start meshing or solving a step until the user has approved it.

`scripts/wall_breakdown.py <case>` splits the wall force into nose, cylinder,
fins and base (C_Af, C_Ab), roll by fin chord and side, and base Cp by
radius. It writes `results/<case>/wall_breakdown.csv` and needs the OpenFOAM
environment for its `postProcess` sampling pass.

## Validation target

Report B1F1 at Ma 1.6, α = 0, φ = 0 gives `C_ℓ = −0.016`, with
`S_ref = πd²/4` and `L_ref = d`. That is the same normalisation as
`post_process.py`, so the target is **`Mx = +0.016`** in our sign convention
(see the sign note below). Acceptance is at least 80 % of that,
`Mx ≥ +0.0128`. A second facility (MDAC S-256, fin-only roll) reads about
−0.010, so the experimental spread is itself large.

Axial force targets from the same report (circle series, Ma 1.6, α = 0):
`C_Af ≈ 0.40` (forebody, `C_A − C_Ab`) and `C_Ab ≈ 0.10` (base, over the full
base area πd²/4). `post_process.py` gives only the total `Cx`. Splitting off
the base disk (fuselage faces at x = 0.8 m with the normal along x) needs a
wall-field sample. The best case so far gives C_Af 0.454 and C_Ab 0.202, so the
base is the main axial-force discrepancy.

The roll is a small residual of two large opposing chordwise contributions
(front chord ≈ +0.040, mid chord ≈ −0.046). Small errors in either one swing
the result by tens of percent.

## Reference configuration — report body B1 + fin F1

Every dimension in `geometry/arc_stabilizers.scad` traces to a wrap-around-fin
(WAF) wind-tunnel report and is **not** free to tidy up. Scale is 1 in = 20 mm
(report D = 4.00 in → our D = 80 mm). Figure 1 gives the body, Figure 2 the fin
layout, Table 4 row F1 the fin numbers.

| Scad | mm | Report | in |
|---|---|---|---|
| `D` | 80 | body diameter | 4.00 |
| `total_len` = 10·D | 800 | overall length | 40.00 |
| `ogive_rho` = 8.5·D | 680 | nose arc radius | 34.00 |
| `ogive_len` = 2.0·D | 160 | nose length | 8.00 |
| `R_in` | 36 | R₃ inner surface | 1.800 |
| `R_edge` | 38 | R₁ mid-thickness | 1.900 |
| `R_out` | 40 | R₂ outer surface | 2.000 |
| `R_out - R_in` | 4 | t, thickness | 0.200 |
| `L` | 140 | C_R, root chord | 7.0 |
| `R - root_embed` | 38 | A, anchor radius | 1.900 |
| `delta` | 45° | δ*, LE/TE wedge | 45 |

`xi` is the one value the report does not state directly. Table 4 gives
`b/2 = 2.64 in`, and `S_f = C_R × b/2 = 18.48 in²` with taper ratio 1.00 — so
`b/2` is the span of the unwrapped rectangular panel, i.e. the developed arc
length at R₁, giving `xi = 2.64/1.900 rad = 79.61°`. Measuring the Figure 2
rear-view sketch independently gives 81.7°. Do **not** round it to 90°: four
fins at 90° tile 360° exactly, leaving zero stow clearance.

The nose is a **secant** ogive, not tangent. The report's three nose dimensions
(34.00R, centre 12.19 in aft and 31.74 in below the axis) are mutually
consistent to four decimals and put the shoulder at 7.99 ≈ 8.00 in. The scad
solves the centre from `ogive_rho` and `ogive_len` so it stays parametric in D;
the 7.07° slope kink where the ogive meets the cylinder is real geometry, not a
bug. Forcing tangency instead (the pre-2026-09 form) stretches the nose to
2.87 D.

Body B1 has **no** aft step. Report configurations B2/B3 reduce the last
7.0/4.0 in to 3.60 in diameter so stowed fins sit flush at 4.00 — we model B1,
so that step is deliberately absent.

### Geometry gotchas

- **`caseProperties` wins over the scad defaults.** `rebuild-mesh.sh` reads
  `D, N, xi, L` from `<case>/constant/caseProperties` and passes them to
  OpenSCAD with `-D`. Editing the defaults at the top of the `.scad` alone
  changes nothing for a real case — change both, or the mesh silently rebuilds
  with the old value.
- **`wing_center` is coupled to `xi`.** The fin is a symmetric bow: the arc
  centre sits at `[-R_edge·cos(xi/2), root_y + R_edge·sin(xi/2)]` so the
  mid-thickness root *and* tip both land on the φ = 0 axis, the tip at
  `root_y + 2·R_edge·sin(xi/2)`. Changing `xi` without that coupling leaves the
  tip short of the axis. At `xi = 90` both components collapse to the older
  `R_edge/sqrt(2)`.
- **`mirror([1,0,0])` in `stabilizers()` sets the roll direction**, making each
  fin sweep toward +Z. Dropping it mirrors the whole model and flips the sign of
  the induced roll.
- **The cylinder rungs in `body_profile()` are currently dead code.** They emit
  ~2800 intermediate stations to avoid a sliver quad, but OpenSCAD's Manifold
  backend (2026.09.06 here) merges the coplanar axial quads back — the exported
  fuselage STL has no vertices between `ogive_len` and `total_len`. Confirm
  before relying on them.

## Solver

The solver is **HiSA** (steady-state pseudo-transient density-based compressible
solver), not `rhoCentralFoam`. In `controlDict`, "time" is the outer pseudo-time
iteration counter, not physical seconds: `endTime 5000` is the iteration budget,
`deltaT 1`, `writeInterval 100`, and `startFrom latestTime` makes interrupted
runs resumable. The `time` column in `coefficients.csv` is therefore an
iteration index — read converged coefficients from the tail, not any single row.

## Pipeline (two stages)

```bash
./rebuild-mesh.sh openfoam/arc_no_stab           # OpenSCAD → STL → blockMesh → parallel snappyHexMesh -overwrite → decompose
./run-simulation.sh --dry-run openfoam/arc_no_stab
./run-simulation.sh openfoam/arc_no_stab         # hisa -parallel, reconstruct latest, run post_process.py
```

`rebuild-mesh.sh` preserves case dictionaries and only removes generated
mesh/run artifacts. The mesh runs `mpirun -np 12` by default; override with
`NP=<n>`.
`MAX_CELLS=11000000` is enforced after `checkMesh` by default; set `MAX_CELLS=0`
only for exploratory runs — both tracked cases use `addLayers true` with
`nSurfaceLayers 13` (first layer 3 µm, ratio 1.4, wall-resolved rather than
wall functions), so they sit close to the guard. The fin-refined server cases
exceed 11 M and are run with the guard raised. The hard limit is the server's
memory: HiSA on 12 ranks survives 13.3 M cells and is OOM-killed at 14.6 M.
`checkMesh` has never reported `Mesh OK` for this mesh family, because of
aspect-ratio and skewness flags at the fin tips.

`create_case.py` creates a *new* case by cloning an existing one
(`--template`, default `openfoam/arc`) and rewriting `constant/caseProperties`
plus the Mach-dependent freestream entries. It does not touch mesh
dictionaries, so a `--N 0` clone of `openfoam/arc` still references
`stabilizers.stl` and will fail to mesh — clone `openfoam/arc_no_stab` for
bare-body variants.

**Pass an explicit case path.** Both scripts now default to `openfoam/arc`, but
the quickstart still passes the path explicitly so the two stages provably act
on the same case.

HiSA has no `-dry-run`, so `run-simulation.sh --dry-run` instead validates the
decomposed mesh and patches with parallel `checkMesh` (the usual startup failure
mode), tee'd to `log.checkMesh.dryRun`. A solve writes `log.hisa` and
`log.reconstructPar`; `run-simulation.sh` clears these plus `postProcessing/`
before each run while preserving the mesh and `0/` fields.

Iterating on solver dicts only: edit the case's `system/*` in place — no
regeneration step is needed, since `rebuild-mesh.sh` preserves dictionaries.
Geometry changes always require `rebuild-mesh.sh`.

## Linting and validation

The Python scripts are **stdlib-only** — no `requirements.txt`, virtualenv, or
install step; run them with `python3` directly. Lint with `rtk ruff check .`
(a `.ruff_cache/` is present). There is **no unit-test framework**; validation
is the dry-run path plus sanity checks: `checkMesh -constant -noZero` must report
`Mesh OK`, `postProcessing/forces` must be non-empty, `results/<case>/coefficients.csv`
must be populated, and converged coefficients should fall in physically sane
ranges (e.g. supersonic `Cx` order 0.4–0.7, pressure-dominated).

**Roll is not zero for the finned cases.** Wrap-around fins induce roll at 0°
AoA — that is the effect under study, not a mesh error. Only the clean-body
baseline (`N = 0`) should show `Cy`, `Cz` and all three moments ≈ 0 by
axisymmetry; use that case, not the finned ones, as the force-integration sanity
check. See README "Validation Expectations".

Sign convention: the report draws `+C_ℓ` clockwise in the rear view looking
upstream, which is a right-hand rotation about **−X**, so report `+C_ℓ` = `−Mx`
in our OpenFOAM convention.

Geometry-only check, no OpenFOAM needed — renders the solid, prints the derived
constants and confirms closure:

```bash
openscad -o /tmp/chk.stl -D 'D=80;N=4;xi=79.61;L=140;EXPORT="";' \
  geometry/arc_stabilizers.scad 2>&1 | grep -E 'ECHO|Status|manifold'
```

Pass condition is `3D object (manifold)` with `Status: NoError`. The `ECHO`
lines report `ogive_len`, the ogive centre, the shoulder kink, the edge wedge
run and the fin tip's on-axis radius, for cross-checking against the report.

## Freestream constants — single source of truth

`<case>/constant/freestreamProperties` is the only place freestream
values are written. `0/U`, `0/p`, `0/T`, and `system/postProcess` all
`#include "../constant/freestreamProperties"` and reference `$UInf`, `$pInf`,
`$TInf`, `$rhoInf`, `$qInf`. The Python post-processor reads the same file.

**Important constraint**: OpenFOAM v2512's `#eval{}` cannot construct vectors
(`vector(x,y,z)` or scalar × vector both fail). `UInf` is therefore stored as
a literal `(510 0 0)` alongside the scalar `UInfMag 510;`, and the two must be
updated in lockstep when sweeping Mach. `rhoInf` and `qInf` are derived via
scalar `#eval{}`, which does work.

Per-case overrides are managed by `scripts/create_case.py`. From the requested
Mach it recomputes `UInfMag = Mach·√(γ·RGas·TInf)` and rewrites both `UInfMag`
and the literal `UInf`. It also recomputes the freestream turbulence `kInf` and
`omegaInf`: intensity `I` and eddy-viscosity ratio `muRatio` are held constant
across cases (Sutherland μ read from `thermophysicalProperties`), so k-ω
scales with Mach. `create_case.py` additionally writes `constant/caseProperties`
(`D, N, xi, L, Mach, gamma`), which `rebuild-mesh.sh` reads back for the
OpenSCAD geometry parameters.

**`L` is absolute millimetres, not `L/D`.** `arc_stabilizers.scad` takes
`D, N, xi, L` and has no `LD` or `TD` parameter; `R_in`, `R_edge`, and `R_out`
(36/38/40 mm — Table 4's R₃/R₁/R₂) are likewise absolute, and
`assert(R_out == R)` pins `D` to 80 mm unless all three change together. Earlier `caseProperties` carried dead `LD 1`
and `TD 0.02` keys that never reached the geometry — the arc case's real L/D is
140/80 = 1.75.

## Coefficient extraction

`system/postProcess` (included into `controlDict`'s `functions {}` block) runs
the `forces` function object at every time step against the wall patches
(`fuselage`, plus `stabilizers` when present) with
CofR `(0 0 0)` (nose tip) and `rho rho` — HiSA carries a live `rho` field and
real pressure in Pa, so the object integrates true forces. `force.dat`/`moment.dat`
are therefore in N / N·m (**no coefficient math in OpenFOAM**), and
`scripts/post_process.py` does the `q∞` normalization with `D = 0.08 m`,
`S = πD²/4`, and per-case `q∞` from `freestreamProperties`. It writes two files
under `results/<case-name>/`:
- `coefficients.csv` — `Cx, Cy, Cz, Mx, My, Mz` plus split pressure/viscous
  columns (`_p`, `_v` suffix), normalized by `q∞·S` / `q∞·S·D`.
- `forces.csv` — raw N / N·m from `force.dat`/`moment.dat`, same column layout.

The split exists so the same forces log can be re-normalized against different
reference quantities without rerunning the solver.

**`force.dat`/`moment.dat` format (v2512):** whitespace-separated columns
`time | total(x y z) | pressure(x y z) | viscous(x y z)` — no parentheses, and
no porous column unless porosity is enabled (it isn't). `total = pressure +
viscous`. The parser in `post_process.py` reads `total` directly and keeps the
pressure/viscous split; do not reintroduce a parenthesized or
pressure/viscous/porous-ordered parser (an earlier version assumed the old
format and failed to parse any data row).

## Case naming

Generated cases are untracked and named freely; the convention is to append the
varied parameter with decimal points replaced by `p`, e.g.
`openfoam/arc_M2p0`, `openfoam/arc_no_stab_M2p5`.

## Coordinate system

X axial (nose at origin, base at `x = 10·D = 0.8 m`), +Y is where the first
fin sits, Z lateral. Defined in `geometry/arc_stabilizers.scad` (the profile is
built in the XY plane, then `rotate([0, 90, 0])` puts the axis on X). The
OpenFOAM mesh uses the same axes; `D = 80 mm` in `blockMeshDict`
(scale = 0.001) and `D = 0.08 m` in the post-processor — keep these consistent
if `D` ever changes.

## Generated artifacts — do not hand-edit

- `<case>/constant/triSurface/fuselage.stl` and `stabilizers.stl` — produced by
  OpenSCAD via `rebuild-mesh.sh`, one call per part using the scad's
  `EXPORT="fuselage"|"stabilizers"` selector. `stabilizers.stl` is skipped when
  `N = 0`. Edit `geometry/arc_stabilizers.scad` and regenerate; never `sed` the
  STL. Use `--geometry path/to/alt.scad` to override which geometry file
  `rebuild-mesh.sh` passes to OpenSCAD.
- `<case>/processor*/` — produced by `decomposePar`. Don't edit per-processor
  files; change `system/` dicts on master and re-decompose.
- `<case>/constant/polyMesh/` — must be the final snapped mesh with the
  `fuselage` wall patch. `snappyHexMesh` runs with `-overwrite`, then the case
  is redecomposed from that final mesh.
- `run-simulation.sh` clears `log.hisa`, `log.reconstructPar`,
  `log.checkMesh.dryRun`, and `postProcessing/` before each run, while
  preserving the mesh (`processor*/constant/`) and `0/` initial fields.
- Under `openfoam/`, only case dictionaries are tracked (`0/`, `system/`,
  `constant/*Properties`, `case.foam`). Everything generated is gitignored
  (see `.gitignore` and "Experiment tracking").

## ParaView visualization

`system/postProcess` writes these viz fields at each `writeTime`:
`Ma` (Mach number), `grad(rho)` + `schlieren` (|∇ρ|), `grad(p)` + `magGradP`,
and `Cp` on the wall patches. Open `<case>/case.foam` in ParaView after
`reconstructPar -latestTime` runs.

## OpenFOAM version

Templates are written for **OpenFOAM v2512** (note the `#eval` constraint
above and the `(forces)` / `(fieldFunctionObjects)` libs syntax). Older
versions may need the explicit `.so` suffix or different `forceCoeffs` keys.
