# oxiphysics-articulated TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 2,780 SLoC | ~25 tests (not yet listed in the README status table — see v0.2.0 housekeeping)

## Completed (v0.1.0 – v0.1.3)

### Foundation (shipped 2026-05-11)

- [x] Core types: `ArticulatedModel`, `Body`, `Joint`, `SpatialVec6`, `SpatialInertia`, `SpatialTransform`
- [x] Spatial algebra (`spatial.rs`, 878 LoC): velocity/force/inertia transforms, Featherstone cross-product, 6×6 inertia tensor manipulation
- [x] Joint types: revolute, prismatic (with `JointAxis` enum and `JointState { q, qd, qdd }`)
- [x] RNEA — Recursive Newton-Euler Algorithm (inverse dynamics: given q/qd/qdd, compute joint torques)
- [x] ABA — Articulated-Body Algorithm (forward dynamics: given q/qd/torques, compute qdd)
- [x] `Featherstone` driver struct: `rnea`, `aba`, `compute_mass_matrix`, `compute_gravity_torques`
- [x] Integration tests: 1-DoF pendulum, 2-DoF arm, 6-DoF robot energy conservation
- [x] Wired into umbrella `oxiphysics` crate (2026-05-13): `pub use oxiphysics_articulated as articulated`

### Extended joint types (planned 2026-05-13, shipped within v0.1.x)

- [x] Universal joint (2-DoF, cardan/Hooke joint) — needed for steering columns, driveshafts
- [x] Spherical joint (3-DoF ball-and-socket) — needed for shoulder/hip joints in humanoid robots
- [x] 6-DoF free joint — floating-base robots (SLAM / legged locomotion)
- [x] Helical (screw) joint — lead-screw actuators

### Dynamics extensions (planned 2026-05-13, shipped within v0.1.x)

- [x] Centroidal momentum matrix (CMM) — legged locomotion planning
- [x] Composite rigid-body algorithm (CRBA) — faster mass-matrix computation for dense problems
- [x] Operational-space control (task-space Jacobian + inertia projection) — `osc.rs`
- [x] Soft joint limits via penalty forces

### Baseline verification (2026-06-11)

- [x] Roadmap code audit: `rnea.rs` / `aba.rs` confirmed forward passes only — no derivative code anywhere in the crate; this gap drives the v0.2.0 theme below

## v0.2.0 — Robotics parity (Pinocchio-class derivatives, IK, import)

Everything below builds on the existing forward-pass Featherstone stack (`spatial.rs`, `rnea.rs`, `aba.rs`); no breaking changes to `ArticulatedModel` expected.

### Housekeeping

- [ ] Add oxiphysics-articulated to README status table — the crate is currently absent from the workspace "Implementation Status" table; add its row (2,780 SLoC, ~25 tests) alongside the other crates

### Analytical derivatives

- [x] Analytical derivatives of RNEA (∂τ/∂q, ∂τ/∂q̇) (shipped 2026-06-12)
  - **Goal:** matches central finite differences to 1e-6 and self-consistent adjoint identities to 1e-10 on a 7-DoF arm; ≥10x faster than finite differences.
  - **Design:** Carpentier & Mansard 2018 recursive analytic derivatives over the existing spatial algebra (`spatial.rs`); crate currently has no derivative code (verified: `rnea.rs`/`aba.rs` forward passes only).
  - **Files:** `spatial.rs`, `rnea.rs`, new derivatives module
  - **Files:** new `src/rnea_derivatives.rs` (analytic forward/backward derivative passes) + `src/spatial.rs` helpers + `src/lib.rs` registration
  - **Tests:** 1-DoF pendulum, 2-DoF arm, 7-DoF serial chain; ∂τ/∂q and ∂τ/∂q̇ vs central finite differences (h=1e-6) < 1e-6; spatial-derivative helpers unit-tested vs FD in isolation
  - **Risk:** spatial cross-product derivative operators (∂(X·v)/∂q = −S×(X·v)) are the crux; honest-split FD fallback above verified DoF if analytic parity not reached
  - **Verified 2026-06-12:** analytic forward-mode (tangent) recursion in `src/rnea_derivatives.rs` matches central FD (h=1e-6) to <1e-6 on 1-DoF pendulum, 2-DoF arm, AND 7-DoF serial chain (all 49 entries strict); crux transform-tangent helpers (`joint_transform_tangent`, `d_apply_velocity`, `d_apply_force`) unit-tested vs FD; FD fallback (`rnea_derivatives_fd_h`) retained for multi-DoF joints (universal/spherical/free). 8 new tests, 0 warnings.

- [x] ABA derivatives (∂q̈/∂q, ∂q̇, ∂τ) (shipped 2026-06-12)
  - **Goal:** same tolerances via the identity ∂q̈ = −M⁻¹·∂RNEA.
  - **Design:** reuse RNEA derivatives + Cholesky of CRBA mass matrix; cache factorization per step.
  - **Files:** `aba.rs`, CRBA module, shared derivatives module
  - **Files:** new `src/aba_derivatives.rs` (∂q̈/∂x = −M⁻¹·∂RNEA/∂x, ∂q̈/∂τ = M⁻¹) reusing crba `compute_mass_matrix_crba` + local SPD Cholesky; `src/lib.rs` registration
  - **Tests:** ∂q̈/∂q = −M⁻¹·∂τ/∂q to 1e-8; ∂q̈/∂τ = M⁻¹ matches central FD of aba wrt τ; symmetry/consistency check on M
  - **Risk:** depends on RNEA-derivative correctness; Cholesky must guard SPD (return Result, no unwrap in production)
  - **Verified 2026-06-12:** `src/aba_derivatives.rs` implements ∂q̈/∂τ = M⁻¹, ∂q̈/∂q = −M⁻¹·∂RNEA/∂q, ∂q̈/∂q̇ = −M⁻¹·∂RNEA/∂q̇ via local SPD Cholesky (returns `AbaDerivativeError`, no unwrap). ∂q̈/∂q vs −M⁻¹·∂τ/∂q to 1e-8; ∂q̈/∂τ and ∂q̈/∂q,∂q̈/∂q̇ vs central FD of aba to 1e-6 (2-link), 7-DoF ∂q̈/∂τ vs FD to 1e-5; M⁻¹ symmetry to 1e-10. 9 new tests, 0 warnings.

- [x] Sparse LTL/LTDL mass-matrix factorization (shipped 2026-06-13)
  - **Goal:** factorization+solve beats dense Cholesky ≥2x on a 36-DoF humanoid tree.
  - **Design:** Featherstone 2008 ch. 8 branching-induced sparsity (expanded-parent array), in-place LTDL; optional oxiblas (new dep) only for the dense fallback.
    - **Files:** NEW src/ltdl.rs (312 LoC); MODIFY src/lib.rs
    - **Tests:** λ correctness on branched tree; factor+solve vs dense Cholesky to 1e-10; no-fill-in sparsity; LᵀDL reconstruction to 1e-9; dim-mismatch + non-SPD error paths; solve_columns parity; perf-ratio gate ≥2x vs dense on a 36-DoF tree
    - **Risk:** 0-based index translation in λ — dense-Cholesky parity gate catches index slips
    - **Verified 2026-06-13:** `src/ltdl.rs` implements `build_lambda` (expanded-parent array, 1-based; skips 0-DOF ancestors), in-place `ltdl_factor_inplace` (H = LᵀDL, zero fill-in), three-pass `ltdl_solve`, and the `SparseMassFactorization` struct (+`factor_crba` convenience) with a no-thiserror `LtdlError` (plain enum, manual Display/Error matching `AbaDerivativeError`). Crate stays dependency-free (no criterion: perf test uses a local dense Cholesky + `std::time::Instant`). Hard gates: H·x=b parity to 1e-10 (`test_ltdl_factor_solve_vs_dense_cholesky`); dense/sparse factor+solve ratio ≈7.9x at n=36 (`test_sparse_beats_dense_cholesky_2x`). 9 new tests, 0 warnings, clippy `-D warnings --all-targets` clean. Solve pass-directions corrected vs the original sketch: Phase 1 (Lᵀz=b) descends, Phase 3 (Lx=y) ascends — deep 36-DOF branched tree exposed the off-by-direction (residual 6.6e3 → 2.1e-14); LᵀDL reconstruction ‖LᵀDL−H‖∞ ≈ 4e-16.

### Kinematics, import, and actuation

- [x] IK: damped least squares + null-space projection (shipped 2026-06-12)
  - **Goal:** 7-DoF arm reaches random reachable poses <1 mm / <0.1° in <50 iters; secondary posture task in null space verified.
  - **Design:** Levenberg-Marquardt-damped pseudo-inverse on the OSC Jacobian (exists in `osc.rs`); selectively-damped SVD near singularities; io already defines `IkSolution` types (`robotics_io/types.rs:236`).
  - **Cross-crate:** result types shared with oxiphysics-io `robotics_io`.

- [ ] URDF import
  - **Goal:** load a 7-DoF arm URDF into `ArticulatedModel`; RNEA gravity torques match the hand-built model to 1e-12.
  - **Design:** consume oxiphysics-io `UrdfRobot`/`UrdfLink`/`UrdfJoint`/`UrdfInertial` (exist at `robotics_io/types.rs:181-1068`); map to the 7 existing joint types.
  - **Risk:** parser completeness in io may need finishing (types exist, end-to-end XML parse unverified).
  - **Cross-crate:** oxiphysics-io `robotics_io` URDF parser.

- [ ] Actuator/transmission models
  - **Goal:** reflected rotor inertia (n²·I_rotor) and Stribeck joint friction change ABA results matching MuJoCo `armature`+friction semantics on a pendulum.
  - **Design:** per-joint armature term added to the articulated inertia, friction as torque-level model in the ABA loop.

## v0.3.0 — Loops, contact, and optimal control

Unlocks parallel mechanisms, humanoid contact, and trajectory optimization; loop closure is the item other crates wait on.

- [ ] Kinematic loop closure
  - **Goal:** four-bar linkage simulates with constraint violation <1e-8; validates against the analytic `FourBarLinkage` in oxiphysics-rigid (`mechanism_rigid/types.rs:723`).
  - **Design:** spanning-tree + loop-closure Lagrange multipliers with Baumgarte (or constraint embedding per Featherstone ch. 8); enables parallel robots and the vehicle-crate suspension linkages.
  - **Cross-crate:** oxiphysics-vehicle v0.2.0 "Multibody suspension from real linkage geometry" explicitly depends on this item.

- [ ] Contact-aware articulated dynamics
  - **Goal:** humanoid foot-ground contact resolves with complementarity residual <1e-8; energy decay monotone.
  - **Design:** Delassus operator via sparse LTDL solves + cone-constrained solve delegated to oxiphysics-constraints (Lemke/NCP); Featherstone+LCP per Todorov 2014 (MuJoCo-style soft option too).
  - **Cross-crate:** consumes oxiphysics-constraints direct LCP solvers (its v0.2.0) and NCP (its v0.3.0).

- [ ] DDP/iLQR on articulated dynamics
  - **Goal:** cart-pole and acrobot swing-up converge <100 iterations using analytical derivatives.
  - **Design:** iLQR/DDP already exist generically in oxiphysics-constraints `optimal_control/` (verified: `DifferentialDynamicProgramming`, iLQR, LQR/DARE) — this item wires ABA + analytic derivatives in as the `DynamicsModel`, not a new solver.
  - **Cross-crate:** oxiphysics-constraints `optimal_control/` provides the solver side.

- [ ] MJCF import + floating-base state estimation
  - **Goal:** load a MuJoCo humanoid MJCF; contact-aided floating-base estimator drift <1% over 10 s walking replay.
  - **Design:** MJCF parser in oxiphysics-io (new module alongside robotics_io), complementary filter / EKF on `FreeFloatingJoint` state.
  - **Cross-crate:** new MJCF module lands in oxiphysics-io.

## v1.0 — Production & validation

Freeze the robotics API surface only after derivatives and loop closure have proven the generalized-coordinate layout.

- [ ] Dynamics cross-validation suite
  - **Goal:** RNEA↔ABA roundtrip τ→q̈→τ to 1e-10; CRBA M equals column-wise RNEA M to 1e-12; 10 s free-swing energy drift <1e-8 (already partially covered by existing energy tests — extend to all 7 joint types and floating base).
  - **Design:** golden-model tests, randomized trees via proptest (≤64 bodies, no NaN).

- [ ] Performance gates
  - **Goal:** ABA 7-DoF <1 µs, 36-DoF humanoid <8 µs per call on reference hardware; zero heap allocation in hot loops.
  - **Design:** criterion benches + allocation-counting test harness; gates run as local scripts (no new CI workflow yamls).

- [ ] API stability for the robotics surface
  - **Goal:** `Joint` trait sealed or `#[non_exhaustive]`; semver-checked (cargo-semver-checks local script).
  - **Design:** finalize generalized-coordinate layout docs (q ordering for quaternion joints).

---

History: the pre-2026-06-11 TODO (last updated 2026-06-06) is preserved in full under "Completed"; the forward roadmap derives from the code-verified dynamics design brief of 2026-06-11.
