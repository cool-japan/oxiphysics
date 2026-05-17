# oxiphysics-articulated TODO

Last updated: 2026-05-13

## Phase 1: Foundation ✅ (shipped 2026-05-11)

- [x] Core types: `ArticulatedModel`, `Body`, `Joint`, `SpatialVec6`, `SpatialInertia`, `SpatialTransform`
- [x] Spatial algebra (`spatial.rs`, 878 LoC): velocity/force/inertia transforms, Featherstone cross-product, 6×6 inertia tensor manipulation
- [x] Joint types: revolute, prismatic (with `JointAxis` enum and `JointState { q, qd, qdd }`)
- [x] RNEA — Recursive Newton-Euler Algorithm (inverse dynamics: given q/qd/qdd, compute joint torques)
- [x] ABA — Articulated-Body Algorithm (forward dynamics: given q/qd/torques, compute qdd)
- [x] `Featherstone` driver struct: `rnea`, `aba`, `compute_mass_matrix`, `compute_gravity_torques`
- [x] Integration tests: 1-DoF pendulum, 2-DoF arm, 6-DoF robot energy conservation
- [x] Wired into umbrella `oxiphysics` crate (2026-05-13): `pub use oxiphysics_articulated as articulated`

## Phase 2: Extended Joint Types (v0.2.0)

- [x] Universal joint (2-DoF, cardan/Hooke joint) — needed for steering columns, driveshafts (planned 2026-05-13)
- [x] Spherical joint (3-DoF ball-and-socket) — needed for shoulder/hip joints in humanoid robots (planned 2026-05-13)
- [x] 6-DoF free joint — floating-base robots (SLAM / legged locomotion) (planned 2026-05-13)
- [x] Helical (screw) joint — lead-screw actuators (planned 2026-05-13)

## Phase 3: Dynamics Extensions (v0.2.0)

- [x] Centroidal momentum matrix (CMM) — legged locomotion planning (planned 2026-05-13)
- [x] Composite rigid-body algorithm (CRBA) — faster mass-matrix computation for dense problems (planned 2026-05-13)
- [x] Operational-space control (task-space Jacobian + inertia projection) (planned 2026-05-13)
- [x] Soft joint limits via penalty forces (planned 2026-05-13)
