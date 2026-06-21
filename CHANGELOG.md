# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.3] - Unreleased

### Added

### Changed

### Removed
- `oxiphysics-gpu`: Deleted four orphaned, never-compiled CPU-reference kernel modules
  (`src/kernels/{md,fem,collision,lbm}_kernels.rs`, ~155 KB) that were not declared in
  `kernels/mod.rs` (zero references workspace-wide) and had been superseded by the wired,
  `splitrs`-refactored `kernels/md_force/`, `kernels/rigid/`, `kernels/sph.rs`,
  `kernels/broadphase.rs`. Dead duplicate code left over from a prior refactor.
- `oxiphysics-core` `pde`: removed a dead duplicate `src/pde/types/types.rs` — an orphaned
  pre-refactor module never declared in `pde/types/mod.rs` (which declares only
  `types_2`/`types_3`/`types_impl`) and referenced nowhere in the workspace.

### Fixed
Honesty audit (`/strict-check`): eradicated **silent fabrications** — code that compiled and
returned plausible-but-fake values with no loud marker, where callers could not tell a real
computation never happened. Each was replaced with a real implementation (or an honest error),
plus a regression test asserting the real behavior instead of the fabricated constant. The
workspace has **zero** `todo!`/`unimplemented!` stubs and the CUDA backend (cudarc) was
confirmed honest (real device queries, `NotAvailable`/`FeatureNotEnabled` errors).

- `oxiphysics-core` `neural_ode`: `AdjointMethod::backward` returned the negated input
  (`-loss_grad`) — a fake gradient that silently broke training. Now computes a real
  transpose-Jacobian–vector product of one RK4 step via central finite differences of the
  actual dynamics.
- `oxiphysics-core` `pde`: `HeatEquation3D::step_implicit_split` merely called `step_explicit`
  while advertising operator-split ADI (so it inherited the explicit CFL limit it claimed to
  escape). Now a real LOD backward-Euler ADI with an O(n) Thomas tridiagonal solver in x/y/z;
  unconditionally stable.
- `oxiphysics-core` `collision`: EPA returned a hardcoded `depth: 0.0, normal: [0,1,0]` for
  GJK-confirmed overlaps (the whole EPA loop also failed to enclose the origin for smooth
  shapes), and GJK reported `intersecting: true`/zero-distance on non-convergence. EPA was
  rewritten (origin-enclosing seed tetrahedron, outward winding, horizon-edge expansion, real
  convergence) to return true penetration depth/normal; GJK now reports its real best estimate.
- `oxiphysics-core` `numerical_linear_algebra`: `hessenberg_eigenvalues` returned the raw
  Hessenberg diagonal as "eigenvalues" with `converged: true`, feeding wrong Ritz values to the
  Arnoldi eigensolver. Now a real shifted-QR (Francis) iteration with deflation.
- `oxiphysics-materials` `construction`: `Geogrid::interaction_coefficient` was
  `0.8 * tan(φ)/tan(φ)` — identically 0.8 (NaN at φ=0). Now the real `Ci = tan(δ)/tan(φ)` from a
  stored soil–geogrid interface friction angle; varies with inputs and is finite at 0.
- `oxiphysics-materials` `battery_materials`: `apparent_activation_energy_ev` ignored both
  temperature arguments and returned the stored Ea. Now a real two-point Arrhenius fit from
  `ionic_conductivity(T1)`/`(T2)`; honest `NaN` for degenerate inputs.
- `oxiphysics-fem` `modal`: the Lanczos modal eigensolver returned raw Krylov basis vectors as
  "mode shapes" and diagonalized a mis-indexed tridiagonal matrix (off-by-one sub-diagonal)
  while reporting convergence. Replaced with a real EISPACK-`tql2` symmetric-QL solver with
  eigenvector accumulation and proper Ritz vectors `V·y` (M-normalized); fixed a start vector
  that silently dropped antisymmetric modes.
- `oxiphysics-constraints` `islands`: `build_island_dependency_graph` discarded its work and
  returned an empty `Vec` for all inputs (its tests passed vacuously). Now takes the constraints
  and builds the real cross-island edge set (deduplicated, with honest per-edge counts).
- `oxiphysics-sph` `sph_analysis`: `SphDiagnostics::compute` hardcoded `disorder_parameter: 0.0`
  despite a real helper existing. Now computes the data-driven neighbor-count disorder.
- `oxiphysics-sph` `dfsph_solver`: `PressureSolveIter::iterate` computed the kernel gradient then
  discarded it (`let _ = g;`), applying no velocity correction while claiming to. Implemented the
  missing DFSPH constant-density pressure projection; verified to drive density error below
  tolerance end-to-end.
- `oxiphysics-io` `particle_data_io`: `H5partReader::from_bytes` checked the magic then returned
  an empty reader, silently dropping everything `H5partWriter::to_bytes` wrote. Now the real
  inverse parse (bounds-checked, honest `Err` on truncation); full round-trip restored.
- `oxiphysics-io` `mesh_io`: `GltfMeshReader::parse_json` returned empty geometry while the writer
  emitted real base64 vertex/index buffers. Now decodes the base64 buffer and accessors (via the
  existing `decode_f32`/`decode_u32` plus a pure-Rust base64 decoder); writer↔reader round-trip
  recovers vertices and indices.
- `oxiphysics-wasm` `fluid_bridge`: `compute_ftle` returned a hardcoded `|sin·cos/t|` pattern
  independent of the flow. Now a real finite-time Lyapunov exponent — tracer advection (RK4) →
  flow-map Jacobian → `(1/t)·ln√λ_max(FᵀF)`; honest `Err` when no velocity field is set.
  Validated against analytic flows (uniform → 0, linear strain → strain rate).
- `oxiphysics-wasm` `sim_controls`/`web_worker`: `StepResult.perf_ms` and the worker's `step_us`
  were fabricated constants (`frame_time*0.01`, `100`) and `contact_count` was hardcoded `0`. Now
  measured with a real monotonic clock (`Instant` native / `performance.now()` on wasm32), and the
  worker preview is honestly documented as contactless ballistic integration.

Follow-up honesty sweep (the deferred latent tier — unreachable-but-fabricating GPU paths,
zero-returning "placeholder" methods, and docs/names that overclaimed). Same discipline: each
became a real implementation (or an honest error) with a regression test asserting real behavior.

- `oxiphysics-constraints` `gpu_constraint_solver`: `solve_gpu` looped a **no-op** `WgpuBackend`
  dispatch and read the *unchanged* uploaded buffers back, returning the inputs as a "solved"
  system with `used_gpu: true`. Now runs the real projected-Gauss-Seidel sweep through a shared
  `run_pgs` helper that `solve_cpu` also calls (GPU-path ≡ CPU-path, asserted to 1e-6), and
  `used_gpu` honestly reflects `backend.is_available()`. Corrected `oxiphysics-gpu` `WgpuBackend`
  docs that falsely claimed the stub executed kernels / stored shaders (it is a no-op CPU
  emulation; the real on-device path is `WgpuBackendReal`).
- `oxiphysics-fem` `boundary_element`: `DualBem::assemble` returned an all-zero matrix and RHS
  ("Simplified placeholder"). Now a real dual-BEM system — displacement BIE on boundary elements
  (Kelvin `U`/`T` kernels) plus the hypersingular traction BIE on crack elements (derived `D`/`S`
  kernels with a Hadamard finite-part self-term). Validated against the closed-form penny-crack
  opening `Δu_z = p(1−ν)a/μ` and `K_I` load-scaling. Split into `boundary_element_dual.rs` to keep
  both files <2000 lines.
- `oxiphysics-fem` `fluid_structure`: `MonolithicFsi` exposed only a placeholder `zero_matrix`.
  Added a real coupled saddle-point tangent `assemble_tangent` (continuity `B`/`−Bᵀ` blocks,
  symmetric added-mass interface coupling, PSPG `−τ·L_pp` stabilization) that validates block
  sizes and returns an honest `Error` on mismatch.
- `oxiphysics-wasm` `web_worker`: `SimCommand::ApplyImpulse` discarded the impulse entirely and
  reported `StepDone { step_us: 0 }`. The preview now tracks per-body velocity; `ApplyImpulse`
  applies a real `Δv = impulse/m` (documented unit-mass preview) that bends the ballistic path,
  `step_us` is really measured, and an out-of-range handle returns an honest `Error`.
  `RequestSnapshot` now reports the real tracked velocities instead of hardcoded zeros.
- `oxiphysics-wasm` `engine`: removed `webgpu_compute_placeholder` — a "WebGPU compute" mock that
  did no GPU work, returning body-count/time padded with zeros; its one unique value is now a
  clean `time()` query with no false GPU framing.
- `oxiphysics-md` `quantum_chemistry_md`: `ionic_forces` returned `()` while computing and
  discarding the `Vec` its doc promised (with an inverted, attractive sign); now returns the real
  repulsive forces obeying Newton's third law. `neb_force`/`climbing_image_neb` fabricated a
  `-0.1·x` "true force"; both now take caller-supplied atomic forces and perform the real NEB
  tangential projection (the climbing image inverts the tangential component).
- `oxiphysics-sph` `thermal_sph`: `von_mises_thermal` discarded its argument and returned `0`
  behind a doc stating a nonzero formula. Replaced with a real `von_mises_voigt` (genuinely 0 for
  the isotropic thermal-stress state, nonzero for anisotropic stresses) with corrected docs.
- `oxiphysics-gpu` `scheduler`: `AsyncCompute::tick` fabricated `output = vec![0u8; 4]` on a
  lifecycle transition; removed — the simulated lifecycle runs no kernel and so produces no output.
- `oxiphysics-geometry` `origami`: `gaussian_curvature_at` ignored its vertex index and returned a
  bare `f64` while documenting an `Option`. Now the real Descartes/Gauss-Bonnet angle deficit
  `2π − Σθ` over incident facets, returning `Option<f64>` (`None` for boundary/isolated vertices).
- `oxiphysics-md` `solvation`: the Born `entropy_contribution` returning `0` is physically correct
  (constant ε ⇒ ∂ΔG/∂T = 0), so the "returns 0 as a placeholder" wording was corrected and a real
  temperature-dependent `entropy_contribution_with_dielectric_slope` was added.
- Doc/name honesty where the math was already real but the documentation overclaimed:
  `oxiphysics-md` `electrostatics/coulomb` (a genuine NGP-grid FFT Ewald reciprocal estimate, not
  the "full B-spline PME" claimed), `simulation/{core_sim,plain_sim}`, `rare_event` (a
  deterministic `round(n·p)`, not a Bernoulli trial), and `ab_initio_md` Ehrenfest coupling.

Test robustness (failures pre-existing under the optimized parallel test profile, not fabrications):
- `oxiphysics-core` `parallel_orchestrator::test_timing_accumulation`: `black_box` moved inside the
  busy-loop (the optimizer can no longer fold it to a closed form) and the inter-stage work gap
  widened to 10×, so the wall-clock ordering assertion holds under a contended test run.
- `oxiphysics-lbm` `microfluidics::test_hp_mean_velocity`: the unphysical `< 1e-18` absolute
  tolerance (below a ULP at the values' magnitude) replaced with a relative `1e-12` bound — the two
  sides are the exact analytic identity `R²ΔP/(8μL)`, differing only by float rounding.

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

[0.1.3]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.3
[0.1.2]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.2
[0.1.1]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.1
[0.1.0]: https://github.com/cool-japan/oxiphysics/releases/tag/v0.1.0
