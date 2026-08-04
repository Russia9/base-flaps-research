# Repository Guidelines

## Project Structure & Module Organization

This repository is a local CFD workflow for aft-base flap studies. `geometry/arc_stabilizers.scad` defines the parametric OpenSCAD body and arc stabilizers, selecting the exported part via `EXPORT=fuselage|stabilizers`. `openfoam/arc/` (N = 4) and `openfoam/arc_no_stab/` (N = 0 clean-body baseline) are the only tracked OpenFOAM case trees; anything else under `openfoam/` is ignored. Python helpers live in `scripts/`: `create_case.py` clones a case and retargets Mach, and `post_process.py` writes coefficient CSVs under `results/<case>/`. The shell entry points are `rebuild-mesh.sh` and `run-simulation.sh`.

## Build, Test, and Development Commands

Source OpenFOAM v2512 before running mesh or solver commands. In Codex/agent sessions, prefix shell commands with `rtk`.

```bash
rtk ./rebuild-mesh.sh openfoam/arc_no_stab
rtk ./run-simulation.sh --dry-run openfoam/arc_no_stab
rtk ./run-simulation.sh openfoam/arc_no_stab
rtk python3 scripts/create_case.py --force --case openfoam/arc_M2p0 --N 4 --xi 90 --L 140 --Mach 2.0
```

Use `NP=<n>` to override the default 12 MPI ranks (matching `numberOfSubdomains` in the tracked cases). `MAX_CELLS=0` disables the mesh cell-count guard only for exploratory work.

## Coding Style & Naming Conventions

Python uses stdlib-only scripts, 4-space indentation, type hints where useful, `pathlib.Path` for paths, and clear snake_case names. Shell scripts use Bash with `set -euo pipefail`, uppercase configuration variables such as `NP` and `MAX_CELLS`, and explicit validation before destructive cleanup. Keep OpenFOAM dictionary keys and units consistent across the tracked cases; `D` and `L` appear in millimeters for meshing and `D` in meters in post-processing.

## Testing Guidelines

There is no dedicated unit-test framework yet. Validate changes by creating a representative case, rebuilding the mesh, and running `./run-simulation.sh --dry-run`. Before trusting results, require `checkMesh -constant -noZero` to report `Mesh OK`, non-empty `postProcessing/forces` logs, and a populated `results/<case>/coefficients.csv`.

## Commit & Pull Request Guidelines

Recent history uses short, imperative commits, sometimes with Conventional Commit prefixes such as `feat:`. Prefer focused subjects like `feat: add clean-body baseline case` or `mesh: tune snappy refinement`. Pull requests should describe the physical or workflow change, list validation commands run, note OpenFOAM version assumptions, and mention any generated outputs intentionally excluded from git.

## Agent-Specific Instructions

Do not hand-edit generated STL, `processor*/`, `constant/polyMesh/`, or solver time directories. Change `geometry/`, the tracked cases under `openfoam/`, or `scripts/`, then regenerate. Never overwrite unrelated dirty files.
