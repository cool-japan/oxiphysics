# OxiPhysics Development Roadmap

> **Status (2026-04-06):** All 12 phases complete. 49 of 49 items done.

## Phase 1: Foundation
- [x] Project scaffold and workspace setup
- [x] Core types (Vec3, Quat, Transform, AABB)
- [x] Core traits (PhysicsBody, Collidable, Steppable)
- [x] Geometry primitives (Sphere, Box, Capsule)
- [x] Basic collision types (CollisionPair, Contact, ContactManifold)
- [x] Material system
- [x] Error handling across all crates

## Phase 2: Rigid Body Engine
- [x] GJK algorithm implementation
- [x] EPA algorithm implementation
- [x] Broad-phase collision (sweep-and-prune)
- [x] Sequential impulse constraint solver
- [x] Joint types (hinge, ball, slider, fixed)
- [x] Island-based sleeping

## Phase 3: Fluid Simulation
- [x] SPH kernel functions
- [x] SPH neighbor search
- [x] SPH pressure/viscosity solvers
- [x] LBM D3Q19/D3Q27 lattice
- [x] LBM boundary conditions

## Phase 4: Structural / FEM
- [x] Tetrahedral mesh generation
- [x] Linear elastic FEM solver
- [x] Nonlinear FEM
- [x] CalculiX-compatible I/O

## Phase 5: Molecular Dynamics
- [x] Lennard-Jones potential
- [x] Verlet integration
- [x] Neighbor lists
- [x] LAMMPS-compatible I/O

## Phase 6: Soft Body
- [x] Mass-spring systems
- [x] Position-based dynamics
- [x] Cloth simulation

## Phase 7: Vehicle Dynamics
- [x] Tire model (Pacejka)
- [x] Suspension
- [x] Drivetrain

## Phase 8: GPU & Bindings
- [x] GPU broad-phase
- [x] GPU SPH
- [x] Python bindings (PyO3)
- [x] WASM bindings

## Phase 9: Visualization & I/O
- [x] VTK export
- [x] glTF export
- [x] Real-time debug renderer

## Phase 10: Advanced Integration
- [x] Multi-physics coupling (FEM-SPH, FEM-LBM)
- [x] Adaptive time stepping across solvers
- [x] Parallel solver orchestration

## Phase 11: Performance & Optimization
- [x] SIMD-accelerated kernels (SoA batch Vec3 operations)
- [x] Cache-optimized data layouts (ParticleSoA, Morton Z-curve sorting)
- [x] Benchmarking suite with criterion

## Phase 12: Production Readiness
- [x] API stabilization and documentation (stability markers, enhanced docs)
- [x] Example gallery (9 examples: rigid, SPH, FEM, MD, LBM, vehicle, cloth, FSI, benchmark)
- [x] Performance regression tests (microbench, JSON baselines, regression detection)
- [x] Cross-platform CI/CD (GitHub Actions: check, test, clippy, fmt, bench on Linux/macOS/Windows)

---

Last Updated: 2026-04-06 — version 0.1.0
