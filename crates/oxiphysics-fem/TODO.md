# oxiphysics-fem TODO

Last updated: 2026-04-06 / v0.1.0

## Phase 1: Foundation
- [x] Define core types and traits (CsrMatrix, TetrahedralMesh, etc.)
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Linear static analysis (LinearTetrahedron, LinearElasticMaterial, LinearStaticAnalysis)
- [x] Nonlinear solvers: Newton-Raphson, BFGS, Riks arc-length
- [x] Adaptive refinement: h-, p-, hp-refinement
- [x] Beam, shell, truss element types
- [x] Constitutive models: NeoHookean, MooneyRivlin, J2Plasticity, Perzyna, ThermoElastic, Orthotropic
- [x] Dynamic FEM: explicit time integration, modal analysis, eigenvalue analysis, wave propagation
- [x] Coupled physics: thermo-mechanical, FEM-LBM, electrochemical FEM
- [x] Extended formulations: XFEM fracture, discontinuous Galerkin, isogeometric analysis, spectral FEM
- [x] Failure & damage: cohesive zone, damage mechanics, fatigue, buckling
- [x] Multi-scale & ROM: homogenization, multiscale FEM, reduced order modeling
- [x] Meshfree methods, level-set interface tracking
- [x] Geomechanics FEM, soil FEM, probabilistic/stochastic FEM, reliability FEM
- [x] Biomechanics FEM, fluid-structure interaction (FSI)
- [x] Topology optimization FEM
- [x] Crystal plasticity FEM
- [x] Sparse linear algebra: CsrMatrix (CSR), PcgSolver (PCG)
- [x] Integration tests (4,561 tests, 0 stubs, 124 source files)
- [x] Performance benchmarks (basic)

## Phase 3: Polish
- [x] Documentation (rustdoc on all public items)
- [ ] Extended examples
- [ ] Benchmark suite expansion
- [ ] Parallel sparse solver integration
- [ ] Further SIMD / parallel assembly optimization
