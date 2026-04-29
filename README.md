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

All coefficients use **D** as the reference length and **D²** (π D²/4) as the reference area. The moment reference point is the **nose tip**.

| Coefficient | Description |
|---|---|
| C_x | Axial force (drag) coefficient |
| C_y | Normal force (lift) coefficient |
| C_b | Base pressure coefficient |
| m_z | Pitching moment coefficient |

For symmetric configurations (N = 2, 4 at 0° AoA) C_y and m_z will be zero by symmetry. Non-zero values are expected for N = 1 and N = 3.

## Toolchain

| Component | Tool |
|---|---|
| Parametric geometry | OpenSCAD |
| Surface mesh export | STL via OpenSCAD |
| Background mesh | blockMesh (OpenFOAM) |
| Volume mesh | snappyHexMesh (OpenFOAM) |
| CFD solver | rhoCentralFoam (OpenFOAM) |
| Turbulence model | k-ω SST |
| Cloud provider | Scaleway |
| Infrastructure | Terraform |

## Infrastructure

Scaleway instances act as stateless workers. Each worker pulls one case at a time from a shared job queue, runs the full pipeline (geometry → mesh → solve → post-process), uploads results, and picks up the next case. Instances are provisioned at the start of a batch run and terminated when the queue is empty.

```
[Job queue] → [Worker 1]  →  [Results store]
            → [Worker 2]  →
            → [Worker N]  →
```

## Repository Structure (planned)

```
base-flaps-research/
├── geometry/          # OpenSCAD models
├── openfoam/
│   ├── template/      # Base case template (solver settings, BCs)
│   └── cases/         # Generated cases (gitignored)
├── scripts/
│   ├── generate_cases.py   # Expand parameter space → case directories
│   ├── run_mesh.sh          # blockMesh + snappyHexMesh for one case
│   └── post_process.py      # Extract C_x, C_y, C_b, m_z from results
├── infra/             # Terraform + worker bootstrap scripts
└── results/           # Aggregated coefficient tables (CSV)
```
