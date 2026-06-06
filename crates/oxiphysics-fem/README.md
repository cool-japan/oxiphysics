# oxiphysics-fem

**Status: [Stable]** — v0.1.3 (2026-06-06)

[![Tests](https://img.shields.io/badge/tests-4620-brightgreen)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-fem)](https://docs.rs/oxiphysics-fem)

Full-featured Finite Element Method (FEM) library for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.
124 source files, 5,273 public API items, and 4,620 tests covering the breadth of modern computational mechanics.

## Features

- **Element types**: linear/quadratic tetrahedral, beam, shell, truss elements
- **Nonlinear solvers**: Newton-Raphson, BFGS, Riks arc-length continuation
- **Adaptive refinement**: h-, p-, and hp-refinement strategies
- **Constitutive models**: NeoHookean, MooneyRivlin, J2 plasticity, Perzyna viscoplasticity, ThermoElastic, Orthotropic
- **Coupled physics**: thermo-mechanical FEM, FEM-LBM coupling, electrochemical FEM
- **Advanced formulations**: XFEM fracture, discontinuous Galerkin (DG), isogeometric analysis (IGA), spectral FEM
- **Dynamics**: explicit FEM, modal analysis, eigenvalue analysis, wave propagation
- **Failure & damage**: cohesive zone model, damage mechanics, fatigue, buckling analysis
- **Multi-scale & ROM**: homogenization, multiscale FEM, reduced order modeling (ROM)
- **Meshfree & level-set**: meshfree methods, level-set interface tracking
- **Geomechanics & soil**: geomechanics FEM, soil FEM, probabilistic/stochastic FEM, reliability FEM
- **Biomechanics**: biomechanics FEM, fluid-structure interaction (FSI)
- **Topology optimization**: density-based topology optimization FEM
- **Sparse linear algebra**: `CsrMatrix` (CSR format), `PcgSolver` (preconditioned conjugate gradient)

## Key Exports

```rust
use oxiphysics_fem::{
    CsrMatrix, TetrahedralMesh, LinearTetrahedron,
    LinearElasticMaterial, LinearStaticAnalysis, PcgSolver,
};
```

## Notes

- 124 source files, 5,273 public items, 4,620 tests — 0 stubs
- Depends on `oxiphysics-core`, `oxiphysics-geometry`, `oxiphysics-materials`

## License

Apache-2.0 — Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project by COOLJAPAN OU (Team KitaSan)
