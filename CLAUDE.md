# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Parametric CFD study of arc-shaped aft-base fins on a supersonic ogive-cylinder
fuselage, sweeping 108 cases over `(N, ξ, L/D, Ma)`. Pipeline is
OpenSCAD → STL → blockMesh → snappyHexMesh → rhoCentralFoam → Python post-processor.
See `README.md` for the full parameter space and physical setup.

## Pipeline (two stages)

```bash
python3 scripts/create_case.py --force --case openfoam/test --N 2 --xi 45 --LD 1.0 --TD 0.02 --Mach 1.5
./rebuild-mesh.sh openfoam/test           # OpenSCAD → STL → blockMesh → parallel snappyHexMesh -overwrite → decompose
./run-simulation.sh --dry-run openfoam/test
./run-simulation.sh openfoam/test         # rhoCentralFoam -parallel, reconstruct latest, run post_process.py
```

Scripts default to `openfoam/test` where practical. `create_case.py` is the
step that replaces or creates a case from `openfoam/template/`; `rebuild-mesh.sh`
preserves case dictionaries and only removes generated mesh/run artifacts. The
mesh runs `mpirun -np 6` by default to leave two cores free on an 8-core
workstation; override with `NP=<n>`. `MAX_CELLS=2000000` is enforced after
`checkMesh` by default; set `MAX_CELLS=0` only for exploratory runs.
On OpenFOAM builds where parallel `rhoCentralFoam -dry-run` fails with
`MPI_ERR_TRUNCATE`, `run-simulation.sh --dry-run` falls back to a serial dry-run
on the reconstructed master mesh and appends both attempts to the dry-run log.

Iterating on solver dicts only: edit `openfoam/template/system/*`, regenerate
the case with `create_case.py --force`, then run `rebuild-mesh.sh`. Geometry
changes always require `rebuild-mesh.sh`.

## Freestream constants — single source of truth

`openfoam/template/constant/freestreamProperties` is the only place freestream
values are written. `0/U`, `0/p`, `0/T`, and `system/postProcess` all
`#include "../constant/freestreamProperties"` and reference `$UInf`, `$pInf`,
`$TInf`, `$rhoInf`, `$qInf`. The Python post-processor reads the same file.

**Important constraint**: OpenFOAM v2512's `#eval{}` cannot construct vectors
(`vector(x,y,z)` or scalar × vector both fail). `UInf` is therefore stored as
a literal `(510 0 0)` alongside the scalar `UInfMag 510;`. The future
case-generation script must update both lines in lockstep when sweeping Mach.

Per-case overrides are managed by `scripts/create_case.py`; it updates both
`UInfMag` and `UInf` from Mach, `TInf`, `RGas`, and gamma.

## Coefficient extraction

`system/postProcess` (included into `controlDict`'s `functions {}` block) runs
the `forces` function object at every time step against the `body` patch with
CofR `(0 0 0)` (nose tip). It writes raw `force.dat`/`moment.dat` —
**no coefficient math in OpenFOAM**. `scripts/post_process.py` does the
normalization with `D = 0.08 m`, `S = πD²/4`, and per-case `q∞` derived from
`freestreamProperties`. Output: `results/<case-name>/coefficients.csv` with
`Cx, Cy, Cz, Mx, My, Mz` plus split pressure/viscous components.

This split exists so the same forces log can be re-normalized against different
reference quantities without rerunning the solver.

## Coordinate system

X axial (nose at origin, base at `x = 10·D = 0.8 m`), +Y is where the first
fin sits, Z lateral. Defined in `geometry/model.scad`. The OpenFOAM mesh uses
the same axes; `D = 80 mm` in `blockMeshDict` (scale = 0.001) and `D = 0.08 m`
in the post-processor — keep these consistent if `D` ever changes.

## Generated artifacts — do not hand-edit

- `openfoam/test/constant/triSurface/body.stl` — produced by OpenSCAD via
  `rebuild-mesh.sh`. Edit `geometry/model.scad` and regenerate; never `sed` the STL.
- `openfoam/test/processor*/` — produced by `decomposePar`. Don't edit
  per-processor files; change `system/` dicts on master and re-decompose.
- `openfoam/test/constant/polyMesh/` — must be the final snapped mesh with the
  `body` wall patch. `snappyHexMesh` runs with `-overwrite`, then the case is
  redecomposed from that final mesh.
- The default mesh keeps `addLayers false` to avoid increasing the sweep cell
  count. Boundary-layer meshes belong in a separate profile with an explicit
  y+ target and cell budget.
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
