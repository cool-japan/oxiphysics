# oxiphysics-rigid

**Status: Alpha** | Version 0.1.2 | Part of [OxiPhysics](https://github.com/cool-japan/oxiphysics)

Rigid body dynamics for the OxiPhysics engine — Pure Rust, no C/Fortran dependencies.

## Domain Highlights

- **Core dynamics**: impulse-based rigid body integration, broadphase/narrowphase pipeline
- **Articulated systems**: multibody hierarchies, motors, ragdolls, robotic systems
- **Sleep system**: island-based sleeping/waking with activation events
- **Continuous collision detection (CCD)**: swept-volume CCD for fast objects
- **Advanced domains**: aerospace, aircraft, buoyancy, cables, crowd simulation,
  fracture/impact dynamics, granular flow, gyroscopic effects, kinematics,
  locomotion, marine, mechanisms, orbital/spacecraft, satellite, swarm, terrain,
  neural control, fluid↔rigid coupling, deformable coupling
- **Fluid interaction**: fluid coupling, fluid dynamics, fluid-structure interaction

## Publicly Exported Modules (40+)

`body`, `articulated`, `sleeping`, `ccd`, `impulse`, `motors`, `pipeline`, `world`,
`ragdoll`, `fluid_coupling`, `fluid_dynamics`, `fluid_structure`, `multibody`,
`neural_control`, `swarm_dynamics`, `robotic_systems`, `constraint_forces`,
`deformable_coupling`, `aerospace`, `aircraft`, `buoyancy`, `cables`, `crowd`,
`fracture`, `granular`, `gyroscopic`, `impact`, `kinematics`, `locomotion`,
`marine`, `mechanisms`, `orbital`, `spacecraft`, `satellite`, `swarm`, `terrain`, `thruster`

## Key APIs

```rust
use oxiphysics_rigid::{RigidBodySet, integrate_bodies, update_sleeping};
use oxiphysics_rigid::{apply_force_and_wake, BodyActivationEvent};
use oxiphysics_rigid::{body::*, collider::*, sets::*, sleeping::*, ragdoll::*, fluid_coupling::*};
```

- `RigidBodySet` — collection and management of rigid bodies
- `integrate_bodies` — semi-implicit Euler integration step
- `update_sleeping` / `apply_force_and_wake` — sleep/wake lifecycle
- `BodyActivationEvent` — event fired when a body transitions sleep state

## Statistics

| Metric | Count |
|---|---|
| Public items | 4,815 |
| Tests | 3,820 |
| Stubs | 0 |

## License

Apache-2.0 — Copyright COOLJAPAN OU (Team Kitasan)

