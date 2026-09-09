# Aft-Base Flap Aerodynamics Research

Parametric CFD study of arc-shaped aft-base fins as aerodynamic control elements on a supersonic cylindrical body. The goal is to characterize control authority and aerodynamic penalties across the fin parameter space.

## Geometry

### Fuselage

Cylindrical body with a secant ogive nose. All dimensions are normalized by the base diameter **D**.

| Section | Length |
|---|---|
| Secant ogive nose (ρ = 8.5 D) | 2 D |
| Cylindrical section | 8 D |
| **Total** | **10 D** |

The ogive arc has radius ρ = 8.5 D but is truncated to a 2 D nose, so it is
*secant*, not tangent: it meets the cylinder at a **7.07° slope discontinuity**.
That shoulder kink is real geometry from the reference model, not an artifact —
`geometry/arc_stabilizers.scad` solves the arc centre from ρ and the nose length
so the construction stays parametric in D.

### Reference configuration

The geometry reproduces a wrap-around-fin (WAF) wind-tunnel model — report body
**B1** with fin **F1** — at a scale of 1 in = 20 mm (report D = 4.00 in → 80 mm
here). The scad ↔ report mapping and the derivation of `xi` are documented in
`CLAUDE.md`; these are not free parameters.

### Fins

Arc-shaped fins attached to the aft base of the fuselage, extending axially rearward. Each fin is a constant-thickness circular arc. Its mid-thickness arc (`R_edge`) is anchored on the +Y axis at radius `R - root_embed` and sweeps `xi` toward +Z, so the tip returns to the same axis — a symmetric bow. The section tapers to a sharp edge at the leading and trailing arcs through a `delta` included wedge.

```
Cross-section view (perpendicular to axis):

        ← xi →
    ___________
   /           \   ← outer arc,  R_out  = 40 mm (= D/2)
  |  _________ |   ← edge arc,   R_edge = 38 mm (knife edge)
  | /         \|   ← inner arc,  R_in   = 36 mm
  |/           |
```

**Fin parameters** (`geometry/arc_stabilizers.scad`):

| Parameter | Symbol | Scad variable | Current value |
|---|---|---|---|
| Number of fins | N | `N` | 4 (`openfoam/arc`), 0 (`openfoam/arc_no_stab`) |
| Arc angle | ξ | `xi` | 79.61° |
| Fin chord | L | `L` | 140 mm (L/D = 1.75) |
| Arc radii | — | `R_in`, `R_edge`, `R_out` | 36 / 38 / 40 mm (fixed) |
| Edge wedge angle | δ | `delta` | 45° (included) |
| Root anchor depth | R − A | `root_embed` | 2 mm |

`L` is absolute millimetres, not a ratio — `R_in`, `R_edge`, and `R_out` are likewise absolute, and `assert(R_out == R)` pins `D` to 80 mm unless all three are changed together.

ξ is not stated directly in the source report — it is derived from the tabulated
fin semi-span; see `CLAUDE.md` for the derivation and why it must not be rounded
to 90°. The arc centre is tied to ξ/2, so changing `xi` moves it automatically.

**Fin placement:** The first fin is always rooted on the +Y axis. Additional fins are placed at equal angular spacing (360°/N), each sweeping the same rotational sense — which is what produces the induced roll. For odd N the configuration is laterally asymmetric. `N = 0` emits the bare fuselage.

## Parameter Space

The full 108-case factorial sweep (N × ξ × L/D × Ma) is the eventual target but is not currently driven by any script. Work is presently on two hand-tuned cases at Ma 1.5:

| Case | N | ξ | L | Purpose |
|---|---|---|---|---|
| `openfoam/arc` | 4 | 79.61° | 140 mm | Finned configuration |
| `openfoam/arc_no_stab` | 0 | — | — | Clean-body baseline |

Mach variants are produced with `scripts/create_case.py --template openfoam/arc --Mach <M>`.

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

All six are written per time step by `scripts/post_process.py` to `results/<case>/coefficients.csv`, along with split pressure/viscous components.

At 0° AoA the clean body (N = 0) has all off-axial components ≈ 0 by axisymmetry. **Finned cases do not.** Wrap-around fins induce a rolling moment at zero incidence — that induced roll is the effect this study exists to measure, so a non-zero `M_x` is the expected result, not a symmetry violation. `C_y`, `C_z`, `M_y` and `M_z` still vanish for any equally spaced N ≥ 2 by rotational symmetry; N = 1 is the asymmetric case.

Sign convention: the source report draws `+C_ℓ` clockwise in the rear view looking upstream, a right-hand rotation about **−X** — so report `+C_ℓ` corresponds to `−M_x` here.

### Extracting coefficients from a solved case

`run-simulation.sh` runs the post-processor automatically after `reconstructPar`. To (re-)generate the CSV from an already-solved case without rerunning the solver:

```bash
python3 scripts/post_process.py openfoam/arc_no_stab
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
| CFD solver | HiSA (density-based, pseudo-transient steady state) |
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
│   └── arc_stabilizers.scad    # Parametric fuselage + arc stabilizers
│                               # (CLI-overridable: D, N, xi, L, EXPORT)
├── openfoam/
│   ├── arc/                    # Finned case, N = 4 — also the source case for create_case.py
│   │   ├── 0/                  # U, p, T, k, omega, nut, alphat initial/boundary fields
│   │   ├── constant/
│   │   │   ├── caseProperties        # D, N, xi, L, Mach, gamma — consumed by rebuild-mesh.sh
│   │   │   ├── freestreamProperties   # SINGLE source of truth (pInf, TInf, UInfMag, UInf, RGas)
│   │   │   ├── thermophysicalProperties
│   │   │   └── turbulenceProperties
│   │   └── system/
│   │       ├── controlDict             # functions { #include "postProcess" }
│   │       ├── postProcess             # forces, MachNo, schlieren, magGradP, Cp
│   │       ├── blockMeshDict, snappyHexMeshDict, decomposeParDict, …
│   └── arc_no_stab/            # Clean-body baseline, N = 0 (arc minus the stabilizers)
├── scripts/
│   ├── create_case.py          # clone a case, retarget Mach + geometry parameters
│   └── post_process.py         # forces.dat → results/<case>/coefficients.csv (Cx..Mz)
├── rebuild-mesh.sh             # OpenSCAD → STL → blockMesh + parallel snappyHexMesh -overwrite + decompose
├── run-simulation.sh           # dry-run/solve → reconstructPar → post_process
└── results/                    # Per-case coefficient CSVs (written by post_process.py)
```

Anything under `openfoam/` other than `openfoam/arc/` and `openfoam/arc_no_stab/` is gitignored, so generated Mach variants stay untracked.

## Quickstart

Use a Linux shell with OpenFOAM v2512 sourced, for example:

```bash
source /path/to/OpenFOAM-v2512/etc/bashrc
```

Mesh, validate, and solve the clean-body baseline:

```bash
./rebuild-mesh.sh openfoam/arc_no_stab
./run-simulation.sh --dry-run openfoam/arc_no_stab
./run-simulation.sh openfoam/arc_no_stab
```

The finned case is the same three commands against `openfoam/arc`. To create a case at a different Mach (or a different fin count), clone an existing one:

```bash
python3 scripts/create_case.py --force --case openfoam/arc_M2p0 --N 4 --xi 79.61 --L 140 --Mach 2.0
./rebuild-mesh.sh openfoam/arc_M2p0
```

`create_case.py` rewrites `constant/caseProperties` and recomputes the Mach-dependent freestream entries — `UInfMag`, the literal `UInf` vector, `kInf`, and `omegaInf` — in `constant/freestreamProperties`. It does not touch the mesh dictionaries, so the clone inherits the source case's refinement setup. Pass `--N 0` for a clean-body case; `rebuild-mesh.sh` then skips the stabilizer STL entirely. Note that the case's `system/snappyHexMeshDict` and `system/surfaceFeatureExtractDict` must agree with `N` — a clone with `--N 0` from `openfoam/arc` still references `stabilizers.stl`, so clone from `openfoam/arc_no_stab` instead when you want a bare body.

`rebuild-mesh.sh` reads geometry from `constant/caseProperties`, exports `fuselage.stl` (plus `stabilizers.stl` when N > 0) via OpenSCAD's `EXPORT` selector, runs parallel `snappyHexMesh -overwrite`, reconstructs the final snapped mesh into `constant/polyMesh`, then decomposes that final mesh for the solver. Both scripts default to `NP=12`, matching `numberOfSubdomains` in the tracked cases; override with `NP=<n>`. `MAX_CELLS=11000000` is enforced by default after `checkMesh`; set `MAX_CELLS=0` to disable the guard for exploratory runs. `run-simulation.sh` cleans prior run outputs from the case while preserving the mesh and `0/` fields. HiSA has no `-dry-run`, so `run-simulation.sh --dry-run` validates the decomposed mesh and patches with parallel `checkMesh`, tee'd to `log.checkMesh.dryRun`.

Both tracked cases set `addLayers true` with `nSurfaceLayers 15`, so the mesh sits close to the `MAX_CELLS` guard. Review y+ before trusting wall-sensitive quantities.

## Validation Expectations

Before trusting coefficients from a case, require:

- `checkMesh -constant -noZero` reports `Mesh OK`.
- `./run-simulation.sh --dry-run <case>` exits cleanly.
- `postProcessing/forces` contains non-empty force and moment logs.
- `results/<case>/coefficients.csv` is non-empty.
- Wall-function y+ is reviewed on representative cases before using wall-sensitive quantities.
- For the clean-body baseline (`openfoam/arc_no_stab`), `Cy`, `Cz`, and all three moments should be ≈ 0 by axisymmetry at 0° AoA — a useful check on the mesh and the force integration before trusting the finned cases. Use the baseline, not a finned case, for this check: a finned case is *expected* to carry non-zero `Mx`.
- Geometry stage only (no OpenFOAM): `openscad -o /tmp/chk.stl -D 'D=80;N=4;xi=79.61;L=140;EXPORT="";' geometry/arc_stabilizers.scad` must report `3D object (manifold)` and `Status: NoError`.
