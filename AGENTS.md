# Repository Guidelines

## Project Structure & Module Organization

This repository is a local CFD workflow for aft-base flap studies. `geometry/model.scad` defines the parametric OpenSCAD body and fins. `openfoam/template/` is the only tracked OpenFOAM case tree; generated cases such as `openfoam/test/` and `openfoam/cases/*` are ignored. Python helpers live in `scripts/`: `create_case.py` builds one case, `sweep.py` enumerates the 108-case matrix, and `post_process.py` writes coefficient CSVs under `results/<case>/`. The shell entry points are `rebuild-mesh.sh` and `run-simulation.sh`.

## Build, Test, and Development Commands

Source OpenFOAM v2512 before running mesh or solver commands. In Codex/agent sessions, prefix shell commands with `rtk`.

```bash
rtk python3 scripts/create_case.py --force --case openfoam/test --N 2 --xi 45 --LD 1.0 --TD 0.02 --Mach 1.5
rtk ./rebuild-mesh.sh openfoam/test
rtk ./run-simulation.sh --dry-run openfoam/test
rtk ./run-simulation.sh openfoam/test
rtk python3 scripts/sweep.py --dry-run
```

Use `NP=<n>` to override the default 6 MPI ranks. `MAX_CELLS=0` disables the mesh cell-count guard only for exploratory work.

## Coding Style & Naming Conventions

Python uses stdlib-only scripts, 4-space indentation, type hints where useful, `pathlib.Path` for paths, and clear snake_case names. Shell scripts use Bash with `set -euo pipefail`, uppercase configuration variables such as `NP` and `MAX_CELLS`, and explicit validation before destructive cleanup. Keep OpenFOAM dictionary keys and units consistent with the template; `D` appears in millimeters for meshing and meters in post-processing.

## Testing Guidelines

There is no dedicated unit-test framework yet. Validate changes by creating a representative case, rebuilding the mesh, and running `./run-simulation.sh --dry-run`. Before trusting results, require `checkMesh -constant -noZero` to report `Mesh OK`, non-empty `postProcessing/forces` logs, and a populated `results/<case>/coefficients.csv`.

## Commit & Pull Request Guidelines

Recent history uses short, imperative commits, sometimes with Conventional Commit prefixes such as `feat:`. Prefer focused subjects like `feat: add sweep case generator` or `mesh: tune snappy refinement`. Pull requests should describe the physical or workflow change, list validation commands run, note OpenFOAM version assumptions, and mention any generated outputs intentionally excluded from git.

## Agent-Specific Instructions

Do not hand-edit generated STL, `processor*/`, `constant/polyMesh/`, or solver time directories. Change `geometry/`, `openfoam/template/`, or `scripts/`, then regenerate cases. Never overwrite unrelated dirty files.
