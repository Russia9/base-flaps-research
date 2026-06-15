# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Parametric CFD study of arc-shaped aft-base fins on a supersonic ogive-cylinder
fuselage, sweeping 108 cases over `(N, ξ, L/D, Ma)`. Pipeline is
OpenSCAD → STL → blockMesh → snappyHexMesh → HiSA → Python post-processor.
See `README.md` for the full parameter space and physical setup.

The solver is **HiSA** (steady-state pseudo-transient density-based compressible
solver), not `rhoCentralFoam`. In `controlDict`, "time" is the outer pseudo-time
iteration counter, not physical seconds: `endTime 5000` is the iteration budget,
`deltaT 1`, `writeInterval 100`, and `startFrom latestTime` makes interrupted
runs resumable. The `time` column in `coefficients.csv` is therefore an
iteration index — read converged coefficients from the tail, not any single row.

## Pipeline (two stages)

```bash
python3 scripts/create_case.py --force --case openfoam/test --N 2 --xi 45 --LD 1.0 --TD 0.02 --Mach 1.5
./rebuild-mesh.sh openfoam/test           # OpenSCAD → STL → blockMesh → parallel snappyHexMesh -overwrite → decompose
./run-simulation.sh --dry-run openfoam/test
./run-simulation.sh openfoam/test         # hisa -parallel, reconstruct latest, run post_process.py
```

`create_case.py` is the step that replaces or creates a case from
`openfoam/template/`; `rebuild-mesh.sh` preserves case dictionaries and only
removes generated mesh/run artifacts. The mesh runs `mpirun -np 12` by default;
override with `NP=<n>`.
`MAX_CELLS=11000000` is enforced after `checkMesh` by default; set `MAX_CELLS=0`
only for exploratory runs.

**Pass an explicit case path.** The bare-invocation defaults differ between
scripts: `rebuild-mesh.sh` defaults to `openfoam/cases/test`, `run-simulation.sh`
to `openfoam/test`. The quickstart always passes the path explicitly so the two
stages act on the same case.

HiSA has no `-dry-run`, so `run-simulation.sh --dry-run` instead validates the
decomposed mesh and patches with parallel `checkMesh` (the usual startup failure
mode), tee'd to `log.checkMesh.dryRun`. A solve writes `log.hisa` and
`log.reconstructPar`; `run-simulation.sh` clears these plus `postProcessing/`
before each run while preserving the mesh and `0/` fields.

Iterating on solver dicts only: edit `openfoam/template/system/*`, regenerate
the case with `create_case.py --force`, then run `rebuild-mesh.sh`. Geometry
changes always require `rebuild-mesh.sh`.

## Linting and validation

The Python scripts are **stdlib-only** — no `requirements.txt`, virtualenv, or
install step; run them with `python3` directly. Lint with `rtk ruff check .`
(a `.ruff_cache/` is present). There is **no unit-test framework**; validation
is the dry-run path plus sanity checks: `checkMesh -constant -noZero` must report
`Mesh OK`, `postProcessing/forces` must be non-empty, `results/<case>/coefficients.csv`
must be populated, and converged coefficients should fall in physically sane
ranges (e.g. supersonic `Cx` order 0.4–0.7, pressure-dominated; `Cy/Cz` and
moments ≈ 0 for symmetric N = 2, 4 at 0° AoA). See README "Validation Expectations".

## Freestream constants — single source of truth

`openfoam/template/constant/freestreamProperties` is the only place freestream
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
across the sweep (Sutherland μ read from `thermophysicalProperties`), so k-ω
scales with Mach. `create_case.py` additionally writes `constant/caseProperties`
(`D, N, xi, LD, TD, Mach, gamma`), which `rebuild-mesh.sh` reads back for the
OpenSCAD geometry parameters.

## Coefficient extraction

`system/postProcess` (included into `controlDict`'s `functions {}` block) runs
the `forces` function object at every time step against the `body` patch with
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

## Sweep case naming

`scripts/sweep.py` generates 108 cases under `openfoam/cases/` (default `--root`)
with names in the form `N{N}_xi{xi}_LD{ld}_M{mach}`, where decimal points are
replaced by `p` (e.g., `N2_xi45_LD1p0_M1p5`). `sweep.py` must be run from the
repo root so its `from create_case import` resolves the sibling module.

```bash
python3 scripts/sweep.py --dry-run   # print all 108 command triplets
python3 scripts/sweep.py --force     # create all case directories (no mesh/solve)
```

## Coordinate system

X axial (nose at origin, base at `x = 10·D = 0.8 m`), +Y is where the first
fin sits, Z lateral. Defined in `geometry/model.scad`. The OpenFOAM mesh uses
the same axes; `D = 80 mm` in `blockMeshDict` (scale = 0.001) and `D = 0.08 m`
in the post-processor — keep these consistent if `D` ever changes.

## Generated artifacts — do not hand-edit

- `<case>/constant/triSurface/body.stl` — produced by OpenSCAD via
  `rebuild-mesh.sh`. Edit `geometry/model.scad` (or `geometry/arc_stabilizers.scad`)
  and regenerate; never `sed` the STL. Use `--geometry path/to/alt.scad` to
  override which geometry file `rebuild-mesh.sh` passes to OpenSCAD.
- `<case>/processor*/` — produced by `decomposePar`. Don't edit per-processor
  files; change `system/` dicts on master and re-decompose.
- `<case>/constant/polyMesh/` — must be the final snapped mesh with the `body`
  wall patch. `snappyHexMesh` runs with `-overwrite`, then the case is
  redecomposed from that final mesh.
- The default mesh keeps `addLayers false` to avoid increasing the sweep cell
  count. Boundary-layer meshes belong in a separate profile with an explicit
  y+ target and cell budget.
- `run-simulation.sh` clears `log.hisa`, `log.reconstructPar`,
  `log.checkMesh.dryRun`, and `postProcessing/` before each run, while
  preserving the mesh (`processor*/constant/`) and `0/` initial fields.
- Anything under `openfoam/` other than `openfoam/template/` is gitignored
  (see `.gitignore`).

## ParaView visualization

`system/postProcess` writes these viz fields at each `writeTime`:
`Ma` (Mach number), `grad(rho)` + `schlieren` (|∇ρ|), `grad(p)` + `magGradP`,
and `Cp` on the body patch. Open `openfoam/test/case.foam` in ParaView after
`reconstructPar -latestTime` runs.

## OpenFOAM version

Templates are written for **OpenFOAM v2512** (note the `#eval` constraint
above and the `(forces)` / `(fieldFunctionObjects)` libs syntax). Older
versions may need the explicit `.so` suffix or different `forceCoeffs` keys.
