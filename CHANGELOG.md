# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-04-06

### Added
- Initial release of the OxiPhysics unified physics engine
- `oxiphysics-core`: Core types, math primitives, PDE solvers, numerical ODE, statistics, Bayesian optimization
- `oxiphysics-geometry`: B-splines, mesh geometry, computational geometry algorithms
- `oxiphysics-collision`: Broad-phase (SAP) and narrow-phase (EPA/GJK) collision detection
- `oxiphysics-rigid`: Rigid body dynamics, kinematics, mechanism simulation
- `oxiphysics-constraints`: Constraint solving, robot control
- `oxiphysics-vehicle`: Vehicle dynamics simulation
- `oxiphysics-sph`: Smoothed Particle Hydrodynamics fluid simulation
- `oxiphysics-lbm`: Lattice Boltzmann Method fluid simulation
- `oxiphysics-fem`: Finite Element Method structural analysis
- `oxiphysics-md`: Molecular dynamics simulation
- `oxiphysics-softbody`: Soft body dynamics, crack propagation, bio-mechanics
- `oxiphysics-materials`: Material models, smart materials
- `oxiphysics-gpu`: GPU acceleration support
- `oxiphysics-viz`: Visualization and rendering
- `oxiphysics-io`: I/O support for VTK, OpenFOAM, HDF5, medical imaging, and more
- `oxiphysics-python`: Python bindings via PyO3
- `oxiphysics-wasm`: WebAssembly bindings

[0.1.0]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.0
