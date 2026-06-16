# oxiphysics-constraints TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 59,585 SLoC | 2,188 tests

## Completed (v0.1.0 – v0.1.3)

All v0.1.x history is preserved below; new work tracks under the version sections that follow.

### Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation
- [x] Implement primary algorithms (PGS, TGS, island solver)
- [x] Add integration tests (2,188 tests passing)
- [x] Performance benchmarks

### Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization (warm-start, Baumgarte bias)

### Phase 4: Advanced Features
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

### Former "Future / Post-0.1" items — shipped in v0.1.x
- [x] GPU-accelerated constraint solving (`gpu_constraint_solver` — WgpuBackend dispatch path + CPU fallback, WGSL block-PGS kernel)
- [x] SIMD-optimized PGS inner loops (`simd_pgs` — SoA layout + auto-vectorized batch solver)
- [x] Real-time deformable body coupling (`deformable_coupling` — `DeformableBodyState` trait + `RigidDeformableCoupling`)

### Roadmap foundations verified in code (2026-06-11)
- [x] Two-parameter Baumgarte softness exists, bias-only (`BaumgarteVariant::TwoParam{zeta, omega_n}`, `tgs_solver/types.rs`) — full (hz, ζ) parameterization is the open v0.2.0 item
- [x] XPBD groundwork exists (`pbd/` module; `joints/xpbdjoint_traits.rs`) — the unified v2 framework is the open v0.2.0 item
- [x] Kinematic-only `GearPair` exists (`mechanism/types.rs:237`) — the dynamic solver row is the open v0.2.0 item
- [x] Stribeck friction and backlash parameter structs exist standalone (`friction/`, `motor_constraints/`, `mechanism/`) — PGS/TGS row integration is the open v0.2.0 item
- [x] Generic optimal-control stack exists (`optimal_control/` — `DifferentialDynamicProgramming`, iLQR, LQR/DARE); oxiphysics-articulated v0.3.0 wires ABA + analytic derivatives into it as a `DynamicsModel` (no new solver needed here)
- [x] Joint force metering exists (`jointforcemeter_traits.rs`) — reused by the v1.0 joint conformance suite

## v0.2.0 — Soft unification and exact small solvers

Theme: one softness language (hz, ζ) across every solver path, XPBD as a first-class unified backend, and exact direct solvers for small systems.
Exit gate: dt-invariant softness on joints and contacts, 50-joint rigid-limit XPBD chain, and LCP residual <1e-10 on ≤64-contact scenes.

### Softness & XPBD unification

- [x] Soft-constraint parameterization (frequency/damping-ratio) across PGS/TGS (planned 2026-06-12)
  - **Files:** `tgs_solver/types.rs` (`SoftParams` + `BlockConstraint6::solve_iteration_soft`)
  - **Tests:** under-damped overshoots / critically-damped ζ=1 no overshoot; dt-invariant settle within 5% over dt∈{1/30..1/240}; rigid SoftParams reproduces Baumgarte to 1e-9
  - **Risk:** sign/scale conventions — rigid-limit equivalence asserted first
  - **Shipped:** `SoftParams::from_frequency(hz, zeta, h)` with Catto TGS-Soft coefficients; `SoftParams::rigid()` for zero-softening fallback; 15 tests added (2203 total); dt-invariance, damping-ratio step-response, and rigid-equivalence all passing; clippy clean.

- [ ] XPBD unified framework v2
  - **Goal:** contacts, all 9+ joint types, and motors expressible as compliance constraints in one substep loop; chain of 50 joints stable at compliance 0 (rigid limit).
  - **Design:** Müller et al. 2020 "Detailed Rigid Body Simulation with XPBD"; build on existing `pbd/` module and `joints/xpbdjoint_traits.rs`; λ accumulation per substep, no velocity pass needed.
  - **Files:** `pbd/`, `joints/xpbdjoint_traits.rs`
  - **Tests:** 50-joint chain at compliance 0; per-joint-type XPBD-vs-TGS parity scenes
  - **Note:** keep solver choice uniform — PGS/TGS/XPBD selectable per island through one configuration surface

### Exact small-system solvers

- [x] Direct LCP solvers (Lemke + Dantzig) + comparison harness (2026-06-13)
  - **Goal:** exact solve for ≤64-contact systems; harness reports iteration count/residual for PGS vs TGS vs LCP on canonical scenes; LCP residual <1e-10.
  - **Design:** principal pivoting (Dantzig per ODE's `dSolveLCP`) and Lemke with lexicographic anti-cycling; hand-rolled dense `f64` pivoting tableau — nalgebra-free, matching the crate's plain-`f64`-array house style.
  - **Files:** `lcp/{mod,dense,lemke,dantzig,harness}.rs`; tests in `tests/lcp.rs`
  - **Tests:** canonical harness scenes (stack, wedge, high-mass-ratio) with residual/iteration reporting
  - **Risk:** Lemke cycling on degenerate problems — lexicographic ordering must be property-tested
  - **Note:** the harness doubles as the engine for the v1.0 solver conformance matrix
  - **Shipped:** Lemke + Dantzig both shipped; lexicographic anti-cycling property-tested (degenerate + identical-q + 20 random PSD trials); all three canonical scenes solve to residual 0.0 (gate 1e-10); Lemke/Dantzig agree to 1e-10; 2236 crate tests pass, clippy clean.

### Joint-level dynamics

- [x] Mimic/gear dynamic constraint
  - **Goal:** two revolute joints coupled with ratio r maintain |θ₁ − r·θ₂| < 1e-6 rad under load.
  - **Design:** single-row Jacobian constraint (Bullet btGearConstraint); note kinematic-only `GearPair` exists in `mechanism/types.rs:237` — this adds the solver row.
  - **Files:** `mechanism/gear_constraint.rs` (new `GearConstraint` + `BodyState`), `mechanism/mod.rs` (re-export)
  - **Tests:** `rigid_gear_ratio_correctness`, `backlash_dead_zone`, `ratio_one_identity`, `no_energy_injection` — all 4 pass; 2240 crate tests pass, clippy clean.
  - **Note:** kinematic `GearPair` stays for pure-kinematic mechanism use; shipped 2026-06-15

- [ ] Joint friction + backlash in the solver
  - **Goal:** Stribeck breakaway and dead-zone backlash reproduce hysteresis loop on a driven pendulum.
  - **Design:** friction as bounded-force constraint row, backlash as unilateral pair; Stribeck/backlash parameter structs already exist standalone (`friction/`, `motor_constraints/`, `mechanism/`) — integrate into PGS/TGS rows.
  - **Files:** `friction/`, `motor_constraints/`, `mechanism/` (existing parameter structs), PGS/TGS row assembly
  - **Tests:** driven-pendulum hysteresis-loop golden curve

- [ ] SPD (stable PD) implicit drives
  - **Goal:** PD-tracked 20-DoF ragdoll stable at kp=1e4, dt=1/60 (explicit PD diverges).
  - **Design:** Tan et al. 2011 implicit PD: solve with (M + dt·Kd) on the left side; integrate with motor constraint rows.
  - **Files:** `motor_constraints/` (implicit drive rows)
  - **Tests:** 20-DoF ragdoll tracking at kp=1e4, with the explicit-PD divergence case as the control

## v0.3.0 — Exact cones and GPU scale

Theme: exact Coulomb cones beyond pyramid approximations, granular-scale convergence, and a graph-colored GPU backend.
Exit gate: NCP residual 1e-8 where PGS stalls, 10k-contact APGD convergence, and conflict-free GPU batches at 50k contacts.

- [ ] NCP / semismooth Newton contact solver
  - **Goal:** exact Coulomb cone (no pyramid) convergence to residual 1e-8 on 100-contact scenes where PGS stalls at 1e-3.
  - **Design:** Fischer-Burmeister reformulation, semismooth Newton with line search (Erleben 2013 numerical methods for contact); sparse linear solves via internal CG (optional scirs2-sparse).
  - **Tests:** PGS-stall scene suite with NCP-vs-PGS residual comparison
  - **Note:** scirs2-sparse stays optional so the default build remains dependency-light

- [ ] Cone complementarity (CCP) solver à la Chrono
  - **Goal:** 10k-contact granular box converges (residual <1e-4) in <50 iterations with APGD; matches Chrono published convergence behavior qualitatively.
  - **Design:** Anitescu-Tasora CCP, Jacobi + accelerated projected gradient descent (APGD, Mazhar 2015); rayon-parallel projection.
  - **Tests:** 10k-contact granular-box iteration budget; APGD-vs-Jacobi convergence curves
  - **Note:** rayon enters as an optional dep, shared with oxiphysics-rigid's v0.3.0 parallel island solve

- [ ] GPU graph-coloring batch solver v2
  - **Goal:** 50k contacts solved on GPU with zero write conflicts; ≥5x over SIMD PGS on desktop GPU.
  - **Design:** greedy graph coloring on the constraint graph, one WGSL dispatch per color; extends existing `gpu_constraint_solver.rs` (currently batch PGS + CPU fallback) and `simd_pgs.rs`; cudarc path stays feature-gated non-default.
  - **Files:** `gpu_constraint_solver.rs`, `simd_pgs.rs`
  - **Tests:** CPU/GPU result parity on colored batches; coloring-correctness property test (no two constraints in one color share a body)
  - **Note:** wgpu/WGSL remains the default GPU backend (Pure Rust); cudarc only via the deferred track

## v1.0 — Production & validation

Theme: documented solver conformance, joint-level guarantees, and a frozen, fuzz-hardened public API.
Exit gate: conformance matrix published, joint suite green over 60 s scenarios, 1e6 fuzz cases with zero panics.

- [ ] Solver conformance matrix
  - **Goal:** documented accuracy/speed table (PGS/TGS/XPBD/LCP/NCP/CCP) on a fixed scene suite; regressions CI-gated (local script) at ±10%.
  - **Design:** criterion benches + residual assertions; publish as rustdoc page.
  - **Note:** reuses the v0.2.0 LCP comparison harness as its measurement engine

- [ ] Joint conformance suite
  - **Goal:** all joint types hold position drift <1e-6 m / 1e-6 rad over 60 s with limits+motors active; breakable joints fire events deterministically.
  - **Design:** parameterized integration tests reusing `jointforcemeter_traits.rs`.
  - **Tests:** per-joint-type 60 s drift scenarios; breakable-joint deterministic event-replay check

- [ ] Constraint-assembly fuzzing + API freeze
  - **Goal:** proptest over degenerate anchors/axes (collinear, coincident, NaN inputs rejected via `error.rs`) — zero panics in 1e6 cases; public trait surface marked `#[non_exhaustive]` where evolving.
  - **Design:** proptest (workspace dep) + sealed-trait pass.
  - **Note:** aligns with the workspace production-readiness backlog (`#[non_exhaustive]` sweep)

## Deferred / research track
- [~] cudarc (CUDA) backend for the graph-coloring GPU batch solver (feature-gated, non-default per Pure Rust policy) — Ready when: local NVIDIA RTX hardware is available to validate cudarc kernels; the wgpu/WGSL path in v0.3.0 ships independently of this.
