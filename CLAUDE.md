# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Parametric CFD study of arc-shaped aft-base fins on a supersonic ogive-cylinder
fuselage. Pipeline is
OpenSCAD → STL → blockMesh → snappyHexMesh → HiSA → Python post-processor.
See `README.md` for the parameter space and physical setup.

Work is currently on two hand-tuned tracked cases at Ma 1.5, not an automated
sweep: `openfoam/arc` (N = 4, ξ = 90°, L = 140 mm) and `openfoam/arc_no_stab`
(N = 0, the clean-body baseline — `arc` with every stabilizer reference removed).
The 108-case `(N, ξ, L/D, Ma)` factorial is the eventual target but no script
drives it; `scripts/sweep.py` and `openfoam/template/` were removed.

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
`nSurfaceLayers 15`, so they sit close to the guard.

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
ranges (e.g. supersonic `Cx` order 0.4–0.7, pressure-dominated; `Cy/Cz` and
moments ≈ 0 for symmetric N = 0, 2, 4 at 0° AoA). See README "Validation Expectations".

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
(36/38/40 mm) are likewise absolute, and `assert(R_out == R)` pins `D` to 80 mm
unless all three change together. Earlier `caseProperties` carried dead `LD 1`
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
- Anything under `openfoam/` other than `openfoam/arc/` and
  `openfoam/arc_no_stab/` is gitignored (see `.gitignore`).

## ParaView visualization

`system/postProcess` writes these viz fields at each `writeTime`:
`Ma` (Mach number), `grad(rho)` + `schlieren` (|∇ρ|), `grad(p)` + `magGradP`,
and `Cp` on the wall patches. Open `<case>/case.foam` in ParaView after
`reconstructPar -latestTime` runs.

## OpenFOAM version

Templates are written for **OpenFOAM v2512** (note the `#eval` constraint
above and the `(forces)` / `(fieldFunctionObjects)` libs syntax). Older
versions may need the explicit `.so` suffix or different `forceCoeffs` keys.
