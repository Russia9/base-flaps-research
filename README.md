# Aft-Base Flap Aerodynamics Research

Parametric CFD study of arc-shaped aft-base fins as aerodynamic control elements on a supersonic cylindrical body. The goal is to characterize control authority and aerodynamic penalties across the fin parameter space.

## Geometry

### Fuselage

Cylindrical body with a tangent ogive nose. All dimensions are normalized by the base diameter **D**.

| Section | Length |
|---|---|
| Tangent ogive nose (ρ = 8.5 D) | 2.8723 D |
| Cylindrical section | 7.1277 D |
| **Total** | **10 D** |

The ogive is tangent to the cylinder at the junction (no shoulder discontinuity).

### Fins

Arc-shaped fins attached to the aft base of the fuselage, extending axially rearward. The outer surface is flush with the fuselage base (radius = D/2); the inner surface is offset inward by the fin thickness.

```
Cross-section view (perpendicular to axis):

        ← xi →
    ___________
   /           \   ← outer arc, radius D/2
  |  _________ |
  | /         \|   ← inner arc, radius D/2 - t
  |/           |
```

**Fin parameters:**

| Parameter | Symbol | Values |
|---|---|---|
| Number of fins | N | 1, 2, 3, 4 |
| Arc angle | ξ | 30°, 45°, 90° |
| Fin length / diameter | L/D | 0.5, 1.0, 1.5 |
| Thickness / diameter | t/D | 0.02 (parametric) |

**Fin placement:** The first fin is always centered on the +Y axis. Additional fins are placed at equal angular spacing (360°/N). For odd N the configuration is laterally asymmetric.

## Parameter Space

Full factorial sweep:

```
N ∈ {1, 2, 3, 4}
ξ ∈ {30°, 45°, 90°}
L/D ∈ {0.5, 1.0, 1.5}
Ma ∈ {1.5, 2.0, 2.5}

Total: 4 × 3 × 3 × 3 = 108 cases
```

## Flow Conditions

| Parameter | Value |
|---|---|
| Mach number | 1.5 / 2.0 / 2.5 |
| Angle of attack | 0° |
| Regime | Supersonic only |

## Output Coefficients

All coefficients use **D** as the reference length and **πD²/4** as the reference area. The moment reference point is the **nose tip** at the origin. Body-axis convention: +X axial (drag), +Y lateral (first fin), +Z lateral.

| Coefficient | Description |
|---|---|
| C_x | Axial force coefficient (drag) |
| C_y, C_z | Lateral force coefficients |
| M_x | Rolling moment coefficient |
| M_y, M_z | Pitching / yawing moment coefficients |

All six are written per time step by `scripts/post_process.py` to `results/<case>/coefficients.csv`, along with split pressure/viscous components. For symmetric configurations (N = 2, 4 at 0° AoA) the off-axial components vanish by symmetry; non-zero values are expected for N = 1 and N = 3.

### Extracting coefficients from a solved case

`run-simulation.sh` runs the post-processor automatically after `reconstructPar`. To (re-)generate the CSV from an already-solved case without rerunning the solver:

```bash
python3 scripts/post_process.py openfoam/test
```

The script reads:

- `<case>/constant/freestreamProperties` — `pInf`, `TInf`, `UInfMag`, `RGas`, from which it derives `rhoInf` and the dynamic pressure `q∞`.
- the newest `<case>/postProcessing/forces/<time>/force.dat` and `moment.dat` pair. OpenFOAM v2512 writes these as whitespace-separated columns `total | pressure | viscous` (each an `x y z` triple, no parentheses); HiSA integrates against the live `rho` field, so the values are true N / N·m and the script applies the `q∞` normalization.

It writes `results/<case-name>/coefficients.csv`, one row per solver time step:

| Columns | Meaning |
|---|---|
| `time` | solver iteration |
| `Cx, Cy, Cz` | total force coefficients |
| `Mx, My, Mz` | total moment coefficients |
| `Cx_p … Cz_p`, `Mx_p … Mz_p` | pressure contribution |
| `Cx_v … Cz_v`, `Mx_v … Mz_v` | viscous contribution |

It also prints a convergence summary (mean over the last 10% of samples). Use that mean only once `Cx` has plateaued — the early iterations are start-up transients.

## Toolchain

| Component | Tool |
|---|---|
| Parametric geometry | OpenSCAD |
| Surface mesh export | STL via OpenSCAD |
| Background mesh | blockMesh (OpenFOAM v2512) |
| Volume mesh | snappyHexMesh (OpenFOAM v2512) |
| CFD solver | rhoCentralFoam (OpenFOAM v2512) |
| Turbulence model | k-ω SST |
| Post-processing | Python 3 (stdlib only) |
| Visualization | ParaView (open `case.foam`) |
| Primary runtime | Linux with OpenFOAM v2512 sourced |

## Infrastructure

This repository currently implements the local OpenFOAM workflow only. Cloud
workers, job queues, object-store uploads, and Terraform provisioning are future
work unless corresponding files are added.

## Repository Structure

```
base-flaps-research/
├── geometry/
│   └── model.scad              # Parametric fuselage + fins (CLI-overridable: N, xi, LD, TD, D)
├── openfoam/
│   └── template/               # Base case (solver settings, BCs, function objects)
│       ├── 0/                  # U, p, T, k, omega, nut, alphat initial/boundary fields
│       ├── constant/
│       │   ├── caseProperties        # geometry/Mach parameters consumed by rebuild-mesh.sh
│       │   ├── freestreamProperties   # SINGLE source of truth (pInf, TInf, UInfMag, UInf, RGas)
│       │   ├── thermophysicalProperties
│       │   └── turbulenceProperties
│       └── system/
│           ├── controlDict             # functions { #include "postProcess" }
│           ├── postProcess             # forces, MachNo, schlieren, magGradP, Cp
│           ├── blockMeshDict, snappyHexMeshDict, decomposeParDict, …
├── scripts/
│   ├── create_case.py          # template → parameterized case
│   ├── sweep.py                # enumerate/create the 108-case sweep
│   └── post_process.py         # forces.dat → results/<case>/coefficients.csv (Cx..Mz)
├── rebuild-mesh.sh             # OpenSCAD → STL → blockMesh + parallel snappyHexMesh -overwrite + decompose
├── run-simulation.sh           # dry-run/solve → reconstructPar → post_process
└── results/                    # Per-case coefficient CSVs (written by post_process.py)
```

Anything under `openfoam/` other than `openfoam/template/` is gitignored. Generated cases live under `openfoam/test` or `openfoam/cases/*`.

## Quickstart

Use a Linux shell with OpenFOAM v2512 sourced, for example:

```bash
source /path/to/OpenFOAM-v2512/etc/bashrc
```

Create, mesh, validate, and solve the baseline case:

```bash
python3 scripts/create_case.py --force --case openfoam/test --N 2 --xi 45 --LD 1.0 --TD 0.02 --Mach 1.5
./rebuild-mesh.sh openfoam/test
./run-simulation.sh --dry-run openfoam/test
./run-simulation.sh openfoam/test
```

`rebuild-mesh.sh` reads geometry from `constant/caseProperties`, runs parallel `snappyHexMesh -overwrite`, reconstructs the final snapped mesh into `constant/polyMesh`, then decomposes that final mesh for the solver. The scripts default to `NP=6` so two cores remain free on an 8-core workstation; override with `NP=<n>` if needed. `MAX_CELLS=11000000` is enforced by default after `checkMesh`; set `MAX_CELLS=0` to disable the guard for exploratory runs. `run-simulation.sh` cleans prior run outputs from the case while preserving the mesh and `0/` fields. If OpenFOAM's parallel dry-run path hits the known `MPI_ERR_TRUNCATE` failure, `run-simulation.sh --dry-run` retries a serial dry-run on the reconstructed master mesh and keeps both attempts in `log.rhoCentralFoam.dryRun`.

The default template intentionally keeps `addLayers false`. This is the bounded-cell validation mesh for the sweep. Boundary-layer meshes should be introduced as a separate higher-cost profile with an explicit y+ target and cell budget.

List the full 108-case command set:

```bash
python3 scripts/sweep.py --dry-run
```

Create all case directories without meshing or solving them:

```bash
python3 scripts/sweep.py --force
```

Each generated case stores Mach in `constant/caseProperties`; `scripts/create_case.py` updates both `UInfMag` and the literal `UInf` vector in `constant/freestreamProperties`.

## Validation Expectations

Before trusting coefficients from a case, require:

- `checkMesh -constant -noZero` reports `Mesh OK`.
- `./run-simulation.sh --dry-run <case>` exits cleanly.
- `postProcessing/forces` contains non-empty force and moment logs.
- `results/<case>/coefficients.csv` is non-empty.
- Wall-function y+ is reviewed on representative cases before using wall-sensitive quantities; the default sweep mesh does not add boundary-layer cells.
