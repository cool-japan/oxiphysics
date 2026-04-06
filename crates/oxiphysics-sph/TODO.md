# oxiphysics-sph TODO

Last updated: 2026-04-06 | Version: 0.1.0

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
- [x] Pressure solvers: IISPH, PCISPH, DFSPH, DFSPH (full), WCSPH
- [x] Kernel functions and neighbor search
- [x] Adaptive smoothing length (`adaptive_h`, `adaptive_sph`)
- [x] Boundary SPH and open boundary conditions
- [x] Free-surface tracking
- [x] Multiphase / immiscible fluid support
- [x] Surface tension (CSF model)
- [x] Turbulence and turbulence SPH models
- [x] Granular flow SPH
- [x] Viscosity models
- [x] Coupling layer (SPH↔rigid, SPH↔FEM, SPH↔DEM)
- [x] CFL-based adaptive timestepping
- [x] Simulation orchestration module
- [x] 4,367 tests passing, 0 stubs

