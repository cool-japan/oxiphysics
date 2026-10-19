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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 182 (critical 2, high 46, medium 101, low 33); 30 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (2)
- [x] **F0381** (critical, bug, medium/small) `joints/fixedjoint_traits.rs:136` — Joint orientation correction has inverted sign and wrong units
- [x] **F0481** (critical, bug, easy/small) `six_dof_constraint.rs:610` — SixDof solve_position correction has wrong sign, amplifies error

### High (46)
- [x] **F0312** (high, bug, easy/small) `contact_constraints.rs:285` — Baumgarte bias sign inverted: penetration correction drives bodies together
- [x] **F0316** (high, fabrication, hard/medium) `contact_stability.rs:255` — force_closure_test_full is the per-axis sign heuristic, not full test
- [x] **F0318** (high, fabrication, hard/medium) `contact_stability.rs:323` — grasp_quality_epsilon is identically 1.0; Q1 is min wrench norm
- [x] **F0319** (high, bug, easy/small) `contact_stability.rs:833` — lcp_gauss_seidel_step never converges; normal force grows per sweep
- [x] **F0332** (high, bug, medium/medium) `deformable_coupling.rs:188` — Rigid–deformable coupling impulse has wrong sign and bogus bias
- [x] **F0334** (high, stub, medium/small) `deformable_coupling.rs:336` — CoupledSimulation couples everything to rigid body 0, one-way, single pass
- [x] **F0346** (high, bug, easy/small) `game_physics_constraints.rs:737` — RopeConstraint never reconciles velocity with projected positions
- [x] **F0363** (high, bug, medium/small) `holonomic_constraints.rs:713` — dae_step_index1 applies velocity-level λ as a force scaled by dt
- [x] **F0380** (high, bug, easy/small) `joints/breakablejoint_traits.rs:56` — BreakableJoint accumulates impulse forever and reapplies spring per iteration
- [x] **F0384** (high, bug, medium/medium) `joints/revolutejoint_traits.rs:45` — Joint point constraints use wrong effective mass (linear-only or inconsistent)
- [x] **F0385** (high, bug, medium/small) `joints/springjoint_traits.rs:57` — SpringJoint reapplies explicit spring impulse every PGS iteration
- [x] **F0387** (high, stub, medium/medium) `joints/types.rs:218` — RevoluteJoint::with_limits stores limits the solver never enforces
- [x] **F0389** (high, bug, medium/small) `joints/xpbdjoint_traits.rs:55` — XpbdJoint velocity pass applies 1e9 explicit spring and pollutes lambda
- [x] **F0393** (high, bug, easy/trivial) `locomotion.rs:613` — BalanceController::zmp uses CoM velocity instead of acceleration
- [x] **F0394** (high, fabrication, hard/medium) `locomotion.rs:850` — ORCA implementations are simplified fakes; overlap case pushes agents together
- [x] **F0401** (high, bug, easy/trivial) `mechanism/types.rs:461` — FourBarLinkage::follower_angle swaps Freudenstein constants K1 and K2
- [ ] **F0403** (high, bug, easy/small) `motion_planning.rs:620` — PRM build can loop forever and never collision-checks edges
- [ ] **F0404** (high, bug, easy/small) `motion_planning.rs:627` — PRM links nodes without edge collision checks; build can loop forever
- [ ] **F0407** (high, bug, easy/small) `motion_planning.rs:932` — CBiRRT returns goal-to-start path on odd iterations; edges unchecked
- [ ] **F0415** (high, fabrication, hard/medium) `multibody.rs:898` — ABA never propagates articulated inertia, omits velocity-product bias and f_ext
- [ ] **F0416** (high, fabrication, medium/small) `multibody.rs:953` — CRBA sums child composite inertia without spatial transform
- [ ] **F0420** (high, fabrication, hard/medium) `multibody.rs:1303` — hybrid_dynamics runs ABA and RNEA independently, not coupled
- [ ] **F0421** (high, fabrication, easy/small) `multibody.rs:1367` — pseudo_inverse_fat / operational_space_inertia skip inversion for m>2
- [x] **F0424** (high, fabrication, medium/medium) `optimal_control/types.rs:143` — iLQR hardcodes identity cost Hessians; DDP is just iLQR
- [x] **F0425** (high, perf, medium/medium) `optimal_control/types.rs:494` — MPC: finite-difference gradient, fixed step, unenforced state bounds
- [x] **F0430** (high, fabrication, hard/medium) `optimal_control/types.rs:970` — h_inf_upper_bound returns DC-gain Frobenius norm, not H∞ bound
- [x] **F0431** (high, fabrication, hard/medium) `optimal_control/types.rs:985` — h_inf_upper_bound returns DC Frobenius norm, a lower bound
- [x] **F0436** (high, fabrication, medium/small) `pbd/types.rs:182` — XpbdDihedralBending::project rigidly translates all vertices along −Y
- [x] **F0439** (high, fabrication, easy/small) `pbd/types.rs:552` — PbdSimulation::add_volume_constraint builds then discards the constraint
- [x] **F0441** (high, fabrication, easy/small) `pbd/types.rs:764` — XpbdPressureConstraint uses edge lengths and one normal as gradients
- [x] **F0442** (high, fabrication, medium/small) `pbd/types.rs:1428` — XpbdBendingConstraint gradient is not dθ/dp; momentum not conserved
- [x] **F0443** (high, bug, easy/small) `pbd/types.rs:1460` — PbdSubstepSolver never updates velocities or applies forces
- [x] **F0445** (high, fabrication, medium/medium) `pgs_solver.rs:307` — PgsSolver residual is a constraint count; trait hooks never honoured
- [x] **F0446** (high, bug, medium/medium) `pgs_solver.rs:354` — PGS warm start applies same-sign linear impulse to both bodies
- [x] **F0447** (high, fabrication, medium/medium) `pgs_solver.rs:869` — SorPgsSolver fakes SOR by scaling all body velocities by omega
- [x] **F0448** (high, fabrication, medium/small) `pgs_solver.rs:883` — SorPgsSolver 'SOR' multiplies absolute body velocities by omega
- [x] **F0461** (high, bug, easy/small) `robot_control/types/types_impl.rs:171` — Safety monitor overwrites Warning with Normal; NaN passes
- [x] **F0462** (high, bug, easy/small) `robot_control/types/types_impl.rs:392` — Jacobian::manipulability multiplies JJᵀ diagonal instead of determinant
- [x] **F0468** (high, fabrication, medium/small) `servo_constraints.rs:394` — Relay auto-tune ultimate gain is constant 4/π
- [x] **F0474** (high, fabrication, medium/medium) `simd_pgs.rs:369` — Batch PGS ignores lever arms and angular velocity entirely
- [x] **F0476** (high, bug, easy/small) `simd_pgs.rs:499` — SoA PGS never applies warm-start λ; builder fabricates warm value
- [x] **F0478** (high, stub, medium/medium) `six_dof_constraint.rs:426` — 6-DOF angular lock/limit uses hardcoded zero orientation error
- [x] **F0479** (high, bug, hard/medium) `six_dof_constraint.rs:492` — 6-DOF limits clamp impulses by position bounds; no angular error
- [x] **F0491** (high, fabrication, medium/small) `tgs_solver/types.rs:587` — TGS always reports converged and passes dt=0 to solve_position
- [x] **F0494** (high, bug, easy/trivial) `tgs_solver/types.rs:682` — TGS passes dt=0 to solve_position, disabling XPBD joints
- [ ] **F0503** (high, bug, easy/small) `trajectory_optimization/functions.rs:375` — sqp_step finite-difference gradient corrupted by in-place updates

### Medium (101)
- [ ] **F0304** (medium, api, easy/small) `crates/oxiphysics-constraints/Cargo.toml:19` — oxiphysics-gpu dependency pulls wgpu for an always-CPU stub
- [ ] **F0305** (medium, bug, medium/small) `cable_constraints.rs:198` — Catenary sag formula wrong and inclined-span vertex discarded
- [ ] **F0309** (medium, bug, medium/medium) `contact.rs:559` — Contact position correction reuses stale penetration every call
- [ ] **F0310** (medium, bug, easy/small) `contact.rs:783` — Manifold warm-start cache is never written back after solving
- [ ] **F0311** (medium, perf, easy/small) `contact.rs:836` — ContactCache manifold lookup is a linear scan per pair
- [ ] **F0317** (medium, fabrication, hard/medium) `contact_stability.rs:302` — Grasp Q1/epsilon/volume metrics are fabricated proxies
- [ ] **F0321** (medium, bug, easy/small) `contact_theory/types.rs:102` — MachineLearningContact zero-init makes hidden layer untrainable
- [ ] **F0323** (medium, panic, easy/small) `contact_theory/types.rs:277` — connected_components recursive find can overflow stack
- [ ] **F0324** (medium, fabrication, medium/medium) `contact_theory/types.rs:555` — LcpSolver::solve_murty silently runs PGS, not principal pivoting
- [ ] **F0325** (medium, bug, easy/small) `contact_theory/types.rs:767` — Hunt-Crossley damping force lacks stiffness factor (wrong units)
- [ ] **F0326** (medium, bug, easy/trivial) `contact_theory/types.rs:943` — Boussinesq stress_r wrong formula; dead *0.0 terms
- [ ] **F0327** (medium, bug, easy/trivial) `control_theory/functions.rs:396` — discretize_zoh B_d series multiplies extra powers of dt
- [ ] **F0328** (medium, fabrication, easy/small) `control_theory/functions.rs:913` — gain_margin_double_lag returns a hybrid, fabricated finite margin
- [ ] **F0329** (medium, bug, medium/small) `control_theory/types.rs:59` — Smith predictor uses delayed input as model output
- [ ] **F0330** (medium, bug, medium/small) `control_theory/types.rs:68` — Smith predictor subtracts delayed control input as model output
- [ ] **F0333** (medium, bug, easy/trivial) `deformable_coupling.rs:211` — apply_to_deformable passes Δv to trait method documented as impulse
- [ ] **F0335** (medium, bug, easy/small) `energy_management_constraints.rs:270` — Battery/storage update reports power that the SOC clamp did not deliver
- [ ] **F0340** (medium, bug, easy/small) `friction/types.rs:231` — Friction cache eviction removes the newest entry
- [ ] **F0341** (medium, bug, easy/small) `friction/types.rs:383` — PairwiseFrictionSolver accumulates impulse past friction budget
- [ ] **F0343** (medium, bug, easy/small) `friction/types.rs:762` — FrictionSolver::solve_contact ignores masses/dt, always applies full μN
- [ ] **F0344** (medium, bug, easy/small) `friction/types.rs:1043` — Anisotropic friction force omits mass; kinetic clamp is a box
- [ ] **F0345** (medium, bug, easy/trivial) `friction/types.rs:1078` — PatchFriction translational force wrongly scaled by 2/3
- [ ] **F0347** (medium, bug, easy/small) `game_physics_constraints.rs:831` — Wheel brake clamp compares rad/s² with rad/s
- [ ] **F0350** (medium, bug, easy/trivial) `game_physics_constraints.rs:1352` — Warm start applies scaled impulse but keeps unscaled accumulator
- [ ] **F0355** (medium, panic, easy/trivial) `gpu_constraint_solver.rs:477` — GpuConstraintSolver::solve panics on first empty frame with backend
- [ ] **F0357** (medium, panic, easy/small) `gpu_constraint_solver.rs:600` — Public solve() panics on bad body index or lambda length
- [ ] **F0359** (medium, bug, medium/medium) `gpu_constraint_solver.rs:690` — Angular impulse uses scalar inverse mass instead of inverse inertia
- [ ] **F0361** (medium, bug, easy/small) `holonomic_constraints.rs:365` — Null-space projection omits M⁻¹; only correct for unit masses
- [ ] **F0364** (medium, fabrication, easy/small) `holonomic_constraints.rs:926` — compare_accuracy returns hard-coded (C, 0.0)
- [ ] **F0366** (medium, bug, easy/small) `island.rs:214` — Island order follows HashMap iteration; sleep counters misattributed
- [ ] **F0367** (medium, bug, easy/small) `island_solver/functions.rs:46` — Island contact Baumgarte bias adds a length to a velocity (no dt)
- [ ] **F0368** (medium, panic, easy/trivial) `island_solver/types.rs:140` — Self-constraint (body_a == body_b) panics in split_at_mut indexing
- [ ] **F0370** (medium, fabrication, medium/medium) `island_solver/types.rs:218` — ParallelIslandSolver serial, O(C·B) lookups, Joint/Distance silently skipped
- [ ] **F0371** (medium, panic, easy/small) `island_solver/types.rs:249` — Island solvers panic on self-constraint; base_solver unused, O(C·B)
- [ ] **F0374** (medium, perf, easy/small) `island_solver/types.rs:648` — build_islands is O(C·N) and merges islands through static bodies
- [ ] **F0375** (medium, bug, medium/small) `island_solver/types.rs:684` — build_islands sleeps islands instantly; static bodies merge islands
- [ ] **F0377** (medium, bug, medium/small) `islands/functions.rs:88` — Wake propagation crosses static bodies and rebuilds full handle map
- [ ] **F0378** (medium, fabrication, easy/small) `islands/functions.rs:393` — split_island assigns constraints round-robin, contradicting its doc
- [ ] **F0379** (medium, perf, easy/small) `islands/types.rs:424` — build_islands unions through static bodies, merging grounded bodies
- [ ] **F0382** (medium, bug, easy/small) `joints/jointforcemeter_traits.rs:56` — JointForceMeter reports last-iteration delta, double-counts bodies
- [ ] **F0383** (medium, bug, easy/trivial) `joints/motorjoint_traits.rs:44` — MotorJoint clamps per iteration, exceeding max_torque
- [ ] **F0386** (medium, bug, easy/trivial) `joints/types.rs:159` — Lead-screw impulse applies angular part with wrong sign
- [ ] **F0388** (medium, missing-feature, medium/medium) `joints/types.rs:823` — Gear/Pulley/RackPinion joints exported but cannot act in solver
- [ ] **F0395** (medium, stub, easy/trivial) `locomotion_constraints.rs:862` — DcmState::required_vrp ignores omega and kp
- [ ] **F0396** (medium, bug, easy/trivial) `mechanism/functions.rs:357` — cam_double_harmonic reaches only half the lift
- [ ] **F0397** (medium, bug, easy/small) `mechanism/functions.rs:367` — cam_modified_trapezoid is discontinuous and ends at wrong values
- [ ] **F0398** (medium, bug, easy/trivial) `mechanism/gear_constraint.rs:205` — Backlash gear bias drives full drift instead of penetration
- [ ] **F0399** (medium, bug, easy/trivial) `mechanism/types.rs:140` — Ratchet snaps angle to nearest tooth; small advances are lost
- [ ] **F0406** (medium, panic, easy/trivial) `motion_planning.rs:799` — partial_cmp().expect panics on NaN in RRT and warm-start
- [ ] **F0408** (medium, bug, easy/trivial) `motor_constraints/types.rs:855` — Friction compensation subtracted, adding friction instead of cancelling
- [ ] **F0409** (medium, bug, medium/small) `motor_constraints/types.rs:1084` — HarmonicDrive flexspline spring anchored to ground, not input
- [ ] **F0412** (medium, fabrication, hard/medium) `multi_agent_constraints.rs:730` — ORCA half-planes wrong in non-conflict and overlap branches
- [ ] **F0413** (medium, bug, easy/trivial) `multi_agent_constraints.rs:1039` — Empty Voronoi cell returns origin, pulling agents to (0,0)
- [ ] **F0414** (medium, bug, medium/small) `multi_agent_constraints.rs:1138` — Nash best-response cycles on mixed games; panics on empty payoff
- [ ] **F0417** (medium, bug, easy/small) `multibody.rs:1058` — operational_space_inertia returns un-inverted matrix for m>=3 or singular
- [ ] **F0418** (medium, fabrication, easy/small) `multibody.rs:1064` — operational_space_inertia silently skips inversion for task dim > 2
- [ ] **F0419** (medium, bug, medium/small) `multibody.rs:1223` — Cable damping/tensions unused; cable forces never reach ABA
- [ ] **F0422** (medium, bug, medium/small) `optimal_control/functions.rs:150` — spectral_radius power iteration gives wrong answers; is_stable unreliable
- [ ] **F0423** (medium, panic, easy/small) `optimal_control/functions.rs:190` — solve_dare / LQR / Kalman panic on singular matrices
- [ ] **F0427** (medium, bug, easy/trivial) `optimal_control/types.rs:783` — MinimumTimeControl::min_time gives wrong bang-bang time
- [ ] **F0428** (medium, bug, easy/trivial) `optimal_control/types.rs:846` — Pontryagin costate integrates forward; Hamiltonian factor-2 inconsistent
- [ ] **F0429** (medium, fabrication, hard/medium) `optimal_control/types.rs:906` — DifferentialDynamicProgramming is just iLQR despite second-order claim
- [ ] **F0434** (medium, fabrication, medium/medium) `optimization_constraints/functions.rs:418` — socp_solve claims interior-point Newton but is tiny gradient descent
- [ ] **F0435** (medium, bug, medium/small) `optimization_constraints/functions.rs:719` — ADMM z-update uses 10 stale cyclic projections, not true projection
- [ ] **F0450** (medium, bug, easy/trivial) `pid_motor.rs:482` — SetpointProfile::step produces a ramp, not a step
- [ ] **F0451** (medium, panic, easy/small) `pid_motor.rs:571` — from_step_response panics on mismatch and fabricates defaults
- [ ] **F0453** (medium, bug, medium/small) `robot_control/functions.rs:620` — Singular matrix inverses silently return identity or partial result
- [ ] **F0455** (medium, bug, easy/small) `robot_control/types/types_2.rs:123` — ResolvedRateController: fixed 50-sweep Gauss-Seidel, panics on length mismatch
- [ ] **F0457** (medium, bug, easy/small) `robot_control/types/types_2.rs:771` — AdaptiveRobotControl: kp never read, sliding gain hardcoded 10
- [ ] **F0458** (medium, fabrication, medium/small) `robot_control/types/types_2.rs:897` — Gravity feedforward ignores com_offsets; inertia has magic +0.1
- [ ] **F0459** (medium, fabrication, medium/medium) `robot_control/types/types_2.rs:919` — inertia_torques uses invented diagonal inertia m·a²+0.1
- [ ] **F0460** (medium, fabrication, medium/medium) `robot_control/types/types_impl.rs:83` — JointMpc is a heuristic PD, not MPC; rollout slice panics
- [ ] **F0463** (medium, stub, medium/small) `robot_control/types/types_impl.rs:577` — Screw Cartesian trajectory multiplies rotation by 0.0
- [ ] **F0464** (medium, fabrication, hard/medium) `robot_control/types/types_impl.rs:696` — RneaSolver: zero Coriolis, invented 0.1 damping, decoupled gravity
- [ ] **F0465** (medium, fabrication, medium/medium) `robot_control/types/types_impl.rs:708` — RneaSolver is not RNEA: invented damping, no coupling
- [ ] **F0466** (medium, fabrication, medium/small) `robot_control/types/types_impl.rs:988` — Singularity measure is Frobenius norm; never detects singularity
- [ ] **F0470** (medium, bug, easy/trivial) `servo_constraints.rs:602` — home_sequence resets PID every step (derivative kick, no integral)
- [ ] **F0472** (medium, fabrication, hard/medium) `simd_pgs.rs:200` — SIMD PGS: body indices hardwired to 0, angular terms ignored
- [ ] **F0475** (medium, fabrication, medium/medium) `simd_pgs.rs:424` — SoA PGS ignores lever arms and angular velocity entirely
- [ ] **F0477** (medium, bug, easy/small) `six_dof_constraint.rs:403` — Solvers never validate dt; dt=0 (pause) yields NaN in 6-DOF
- [ ] **F0480** (medium, panic, easy/small) `six_dof_constraint.rs:547` — f64/usize clamp panics on inverted or NaN user bounds
- [ ] **F0483** (medium, stub, medium/small) `swarm_robotics.rs:319` — FormationControl ignores formation_type and is undamped
- [ ] **F0484** (medium, fabrication, easy/small) `swarm_robotics.rs:402` — Fiedler value 'approximation' is min_degree/n, not an eigenvalue
- [ ] **F0485** (medium, perf, easy/small) `swarm_robotics.rs:508` — Coverage step rescans grid per agent (O(N²·G)); centroid can panic
- [ ] **F0487** (medium, bug, easy/small) `swarm_robotics.rs:853` — Pheromone diffusion ignores cell size and explicit-Euler stability
- [ ] **F0489** (medium, fabrication, easy/small) `tgs_solver/functions.rs:196` — TGS energy-gain monitor condition can never be true (no-op)
- [ ] **F0492** (medium, fabrication, medium/medium) `tgs_solver/types.rs:594` — TGS step_warm ignores warm-start state; stats hardcoded converged
- [ ] **F0493** (medium, fabrication, medium/medium) `tgs_solver/types.rs:603` — TgsSolver::step_warm scales a cache never applied to constraints
- [ ] **F0495** (medium, bug, easy/trivial) `tgs_solver/types.rs:889` — TGS-Soft divides impulse_scale·λ term by K
- [ ] **F0500** (medium, bug, medium/small) `traits/types.rs:188` — XPBD damped Δλ lacks (1+γ) denominator and mis-scales damping
- [ ] **F0502** (medium, panic, easy/small) `trajectory_optimization/functions.rs:187` — Trajectory-optimisation APIs panic on empty or short inputs
- [ ] **F0504** (medium, panic, medium/medium) `trajectory_optimization/functions.rs:609` — cost_gradient_controls: O(N²) re-simulation plus underflow/index panics
- [ ] **F0506** (medium, fabrication, medium/small) `trajectory_optimization/types.rs:859` — fit_to_waypoints not least squares; panics on one waypoint
- [ ] **F0507** (medium, bug, easy/small) `variational_constraints/functions.rs:160` — log_so3 returns zero/garbage for rotations at or near π
- [ ] **F0509** (medium, fabrication, medium/small) `variational_constraints/types.rs:503` — LieGroupIntegrator claims implicit midpoint but is explicit Euler
- [ ] **F0510** (medium, panic, easy/trivial) `warm_start/functions.rs:946` — select_best_candidate panics on NaN impulse magnitude
- [ ] **F0511** (medium, bug, easy/trivial) `warm_start/types.rs:299` — ImpulseAging decays cumulatively decay^(k(k+1)/2) instead of decay^k
- [ ] **F0512** (medium, bug, easy/small) `warm_start/types.rs:885` — ContactVelocitySolver ignores its restitution/slop/baumgarte/friction fields
- [ ] **F0513** (medium, bug, easy/small) `warm_start/types.rs:922` — ContactVelocitySolver ignores its fields; friction clamps delta only
- [ ] **F0514** (medium, fabrication, easy/small) `warm_start/types.rs:1321` — WarmStartStrategy::Adaptive is a silent no-op
- [ ] **F0515** (medium, bug, easy/trivial) `wearable_constraints.rs:316` — Bowden cable capstan model amplifies delivered tension

### Low (33) — tracked in the external ledger
F0313, F0331, F0338, F0339, F0348, F0351, F0352, F0353, F0356, F0358, F0360, F0362, F0365, F0369, F0373, F0376, F0390, F0391, F0392, F0405, F0411, F0426, F0449, F0456, F0469, F0473, F0482, F0490, F0496, F0497, F0499, F0505, F0508

<!-- AUDIT-2026-10-12 END -->
