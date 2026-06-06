# oxiphysics-softbody

**Status: Stable** | Version 0.1.3 | Part of [OxiPhysics](https://github.com/cool-japan/oxiphysics)

Soft body simulation for the OxiPhysics engine — Pure Rust, no C/Fortran dependencies.

## Domain Highlights

- **Cloth & textiles**: cloth, advanced cloth, membranes, textile, wrinkling, hair/fur
- **Constraint-based**: PBD (Position-Based Dynamics), XPBD, distance/bending/volume constraints
- **FEM soft bodies**: co-rotational FEM elements, elastic waves, peridynamics
- **Fracture & tearing**: fracture dynamics, crack propagation, brittle/ductile modes
- **Ropes & cables**: rope, XPBD rope, Cosserat rod model
- **Inflatables & volumes**: inflatable bodies, volume preservation
- **Biological & surgical**: muscle simulation, tissue/surgical simulation, tendon, biomechanics, morphogenesis
- **Advanced materials**: smart materials, metamaterials, origami, active matter, food physics
- **Aerodynamics**: aerodynamics model for soft-body flight surfaces
- **Haptics & neural**: haptic feedback integration, neural deformation networks
- **Topology optimization**: structural shape optimization under soft constraints
- **Position-Based Fluids**: unified PBF within the soft-body solver

## Publicly Exported Modules (48+)

`aerodynamics`, `cloth`, `cloth_advanced`, `constraint`, `fem_soft`, `fracture`,
`fracture_dynamics`, `inflatable`, `muscle_simulation`, `numerical_softbody`,
`particle`, `pbd_system`, `rope`, `surgery_simulation`, `volume`, `xpbd`, `xpbd_rope`,
`cosserat`, `crack_propagation`, `elastic_waves`, `peridynamics`, `shape_matching`,
`hair`, `membranes`, `metamaterials`, `origami`, `active_matter`, `bio`, `biomechanics`,
`morphogenesis`, `food_physics`, `haptics`, `neural_deformation`, `pbf`,
`smart_materials`, `topology_optimization`, `tissue`, `tendon`, `textile`, `wrinkling`

## Key APIs

```rust
use oxiphysics_softbody::{SoftBodySolver, AerodynamicsModel, FemSoftBody, XpbdSolver};
use oxiphysics_softbody::{DistanceConstraint, BendingConstraint, VolumeConstraint};
use oxiphysics_softbody::{HairStrand, HairSystem, Rope, RopeSolver, XpbdRope, ShapeMatching};
use oxiphysics_softbody::{cloth::*, fracture_dynamics::*, inflatable::*, particle::*};
```

- `SoftBodySolver` trait — unified interface for all soft-body backends
- `AerodynamicsModel` — aerodynamic force model for cloth/membrane surfaces
- `FemSoftBody` / `CorotationalElement` — finite-element soft bodies
- `XpbdSolver` — Extended Position-Based Dynamics solver
- `HairStrand` / `HairSystem` — strand-based hair and fur simulation

## Statistics

| Metric | Count |
|---|---|
| Public items | 4,524 |
| Tests | 3,473 |
| Stubs | 0 |

## License

Apache-2.0 — Copyright COOLJAPAN OU (Team Kitasan)

