# oxiphysics-constraints TODO

Last updated: 2026-06-06 | Version: 0.1.3

## Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Implement primary algorithms (PGS, TGS, island solver)
- [x] Add integration tests (2,188 tests passing)
- [x] Performance benchmarks

## Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization (warm-start, Baumgarte bias)

## Phase 4: Advanced Features
- [x] PBD (Position-Based Dynamics)
- [x] Variational constraints
- [x] Holonomic constraints
- [x] CCD constraints
- [x] Cable constraints
- [x] Haptic constraints
- [x] Trajectory optimization
- [x] Optimal control
- [x] Multi-agent coordination
- [x] Multibody dynamics
- [x] Robot control / locomotion
- [x] PID motor, servo, 6-DOF constraints
- [x] Game physics constraint sets

## Future / Post-0.1
- [x] GPU-accelerated constraint solving (`gpu_constraint_solver` — WgpuBackend dispatch path + CPU fallback, WGSL block-PGS kernel)
- [x] SIMD-optimized PGS inner loops (`simd_pgs` — SoA layout + auto-vectorized batch solver)
- [x] Real-time deformable body coupling (`deformable_coupling` — `DeformableBodyState` trait + `RigidDeformableCoupling`)
