# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2026-06-06

### Fixed
- `oxiphysics-io`: Re-exported `particle_formats` types (`DcdWriter`, `DcdReader`, `XyzWriter`,
  `XyzReader`, `ParticleFrame`, `ParticleTrajectory`, `TrajectoryStats`, `BinaryFrameReader`,
  `BinaryFrameWriter`, `GroReader`, `GroWriter`, `DcdHeader`) at the crate root; a doctest for
  `DcdWriter` was previously failing with an unresolved import.
- `oxiphysics-vehicle`: Fixed `hil.rs` module doctest — replaced non-existent `HilBridge` alias
  with `SimHilBridge` + `HilInterface` trait import; corrected `get_output` return type
  (`Option<f64>`) assertion.
- `oxiphysics-wasm`: Fixed two module-level doctests — `SharedStateBuffer` needed `let mut` to
  call `write_body_position`; `WebGlRenderFrame.positions` is private, call `.positions()` method.

### Added
- `oxiphysics-python` (KF-5 — Python bindings): `PyCsg::union/intersection/subtraction` now
  delegate to the real `oxiphysics-geometry::mesh_boolean::mesh_boolean` algorithm (proper
  winding-number/ray-cast inside-outside classification + `cleanup_mesh`/`weld_vertices`)
  instead of the previous AABB-approximation stubs.
- `oxiphysics-python` (KF-5): `PyPointCloud::poisson_reconstruct` now performs real **Implicit
  Moving Least Squares (IMLS)** surface reconstruction: PCA-estimated normals when absent,
  Gaussian-weighted tangent-plane signed-distance implicit, marching-cubes isosurface extraction
  via `oxiphysics-geometry::signed_distance_field::MarchingCubes`. Previously returned an empty
  mesh silently.
- `oxiphysics-rigid`: `MotionPlanning` now carries a real `obstacles: Vec<PlanningObstacle>` field
  (n-dimensional C-space spheres). `is_collision_free` checks Euclidean distance to every
  obstacle; `is_segment_collision_free` prevents tunnelling in RRT/PRM edge expansion (10-sample
  discretisation). Builder `with_obstacles(...)` provided. Previously the check always returned
  `true`.
- `oxiphysics-lbm`: `LbmGrid3D` now has a real `step(omega: f64)` method: BGK collision
  (equilibrium computed from cached `rho`/`ux`/`uy`/`uz`) + pull-scheme periodic streaming;
  `compute_macroscopic()` updates macroscopic fields after each step. Removed the dead empty
  `step_placeholder` (zero call sites).
- `oxiphysics-md`: New `qm_mm` module — hybrid quantum-mechanics / molecular-mechanics support.
  QM region selectable via `QmMethod` (`Pm3` semi-empirical NDDO, `SccDftb` density-functional
  tight-binding, `Hf` Hartree-Fock/STO-3G, `DftB3lyp` Kohn-Sham LDA/VWN), MM region via
  force-field point charges, with mechanical/electrostatic embedding and hydrogen link-atom
  boundary capping. Each engine runs a real SCF loop (Löwdin-orthogonalised generalised
  eigensolver) and provides numerical Hellmann-Feynman forces.

### Changed
- Workspace-wide structural refactor (splitrs): oversized modules (>2000 LoC) — e.g.
  `oxiphysics-core` `types.rs`/`pde`/`linalg`, `oxiphysics-geometry` `bspline`, `oxiphysics-rigid`
  `kinematics`/`mechanism_rigid`, `oxiphysics-lbm` `mixing_lbm`, `oxiphysics-md`
  `quantum_chemistry`, `oxiphysics-viz` `multiphysics_viz` — were split into focused
  `types`/`functions` submodules, and redundant `#[allow(...)]` Clippy attributes were removed
  across the workspace. Public APIs are unchanged.
- `oxiphysics-sph`: Refactored particle data structures across the SPH crate for clearer field
  layout and reduced duplication (no behavioural change).

### Tests
- 16 new unit/integration tests across `oxiphysics-python`, `oxiphysics-rigid`, `oxiphysics-lbm`
  covering the correctness fixes above (CSG operations, IMLS reconstruction, obstacle collision,
  segment collision, RRT path, LBM mass conservation, equilibrium fixed-point).

## [0.1.1] - 2026-05-17

### Changed
- Feature-gated Python bindings (`oxiphysics-python`) for optional PyO3 dependency
- Version bump for internal crate consistency across workspace
- Continued Pure Rust implementation (100% C/Fortran-free)

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

[0.1.2]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.2
[0.1.1]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.1
[0.1.0]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.0
