# oxiphysics-rigid TODO

Last updated: 2026-06-06 | Version: 0.1.2

## Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Implement primary algorithms
- [x] Add integration tests
- [x] Performance benchmarks

## Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization

## Implemented Modules
- [x] Rigid body integration (impulse-based, semi-implicit Euler)
- [x] Broadphase / narrowphase collision pipeline
- [x] Articulated body / multibody hierarchy
- [x] Motors and joint constraints
- [x] Island-based sleeping system with `BodyActivationEvent`
- [x] Continuous collision detection (CCD)
- [x] Ragdoll system
- [x] Robotic systems and neural control
- [x] Fluid coupling, fluid dynamics, fluid-structure interaction
- [x] Deformable coupling
- [x] Aerospace, aircraft, spacecraft, satellite, orbital mechanics
- [x] Buoyancy and marine dynamics
- [x] Cable and rope dynamics
- [x] Crowd and swarm simulation
- [x] Fracture and impact dynamics
- [x] Granular flow
- [x] Gyroscopic effects
- [x] Kinematics and locomotion
- [x] Mechanisms and thruster systems
- [x] Terrain interaction
- [x] Constraint forces
- [x] 3,820 tests passing, 0 stubs

## v0.1.2 correctness fixes (2026-06-01)
- [x] `MotionPlanning` — Added `PlanningObstacle` (n-D C-space sphere) and `obstacles: Vec<PlanningObstacle>` field. `is_collision_free` now checks Euclidean distance to all obstacles; `is_segment_collision_free` added (10-sample discretisation) to prevent tunnelling in RRT/PRM edge expansion. Builder `with_obstacles(...)` provided. Previously `is_collision_free` always returned `true`. 6 new tests.

