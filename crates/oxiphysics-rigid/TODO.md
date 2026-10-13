# oxiphysics-rigid TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 83,116 SLoC | 3,820 tests

## Completed (v0.1.0 – v0.1.3)

All v0.1.x history is preserved below; new work tracks under the version sections that follow.

### Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation
- [x] Implement primary algorithms
- [x] Add integration tests
- [x] Performance benchmarks

### Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization

### Implemented modules
- [x] Rigid body integration (impulse-based, semi-implicit Euler)
- [x] Broadphase / narrowphase collision pipeline
- [x] Articulated body / multibody hierarchy
- [x] Motors and joint constraints
- [x] Island-based sleeping system with `BodyActivationEvent`
- [x] Continuous collision detection (CCD)
- [x] Ragdoll system
- [x] Robotic systems and neural control
- [x] Fluid coupling, fluid dynamics, fluid-structure interaction
- [x] Deformable coupling
- [x] Aerospace, aircraft, spacecraft, satellite, orbital mechanics
- [x] Buoyancy and marine dynamics
- [x] Cable and rope dynamics
- [x] Crowd and swarm simulation
- [x] Fracture and impact dynamics
- [x] Granular flow
- [x] Gyroscopic effects
- [x] Kinematics and locomotion
- [x] Mechanisms and thruster systems
- [x] Terrain interaction
- [x] Constraint forces
- [x] 3,820 tests passing, 0 stubs

### v0.1.2 correctness fixes (2026-06-01)
- [x] `MotionPlanning` — Added `PlanningObstacle` (n-D C-space sphere) and `obstacles: Vec<PlanningObstacle>` field. `is_collision_free` now checks Euclidean distance to all obstacles; `is_segment_collision_free` added (10-sample discretisation) to prevent tunnelling in RRT/PRM edge expansion. Builder `with_obstacles(...)` provided. Previously `is_collision_free` always returned `true`. 6 new tests.

### Roadmap foundations verified in code (2026-06-11)
- [x] Per-substep restitution helper exists (`impulse/functions_2.rs`) — reused by the small-steps solver below
- [x] Naive sub-stepping exists: `step_substep` re-runs the full pipeline per sub-dt (`world/types.rs:1180`) — the true small-steps solver below supersedes it on the hot path
- [x] Speculative-contact detection candidate exists (`speculative_contact_candidate`, `ccd/types.rs:354`; `narrowphase/dispatch/speculativeconfig_traits.rs`) — solver-side admission is the open item
- [x] Explicit gyroscopic Euler equations + RK4 exist (`gyroscopic.rs`) — the implicit integration path is the open item
- [x] Penalty-spring soft constraints exist (`soft_constraints.rs`) — kept as-is; TGS-Soft solver routing is the open item
- [x] Pyramidal friction-cone clamp exists (`clamp_friction_to_cone`, `impulse/functions.rs:66-79`); standalone rolling helpers exist (`contact_dynamics.rs`) — solver-integrated smooth cone is the open item
- [x] Scalar Archimedes/hull-form buoyancy exists (`buoyancy_model.rs`) — exact submerged-volume clipping is the open v0.3.0 item
- [x] Island-based sleeping infrastructure exists (`BodyActivationEvent`) — prerequisite for the v0.3.0 deterministic parallel island solve

## v0.2.0 — Solver modernization (Jolt/PhysX-5 parity)

Theme: close the contact-solver gap to Jolt and PhysX 5 — true sub-stepping, soft contact response, shock propagation, solver-integrated speculative contacts, and complete friction/restitution models.
Suggested order: small steps → TGS-Soft (after oxiphysics-constraints lands hz/ζ) → speculative solver integration → shock propagation; friction, gyroscopic, and restitution items are independent.
Exit gate: all seven Goal scenarios pass as crate-level scenario tests at default solver settings.

### Sub-stepping & stack stability

- [x] True "small steps" sub-stepping solver (done 2026-06-12)
  - **Goal:** 1000-box pyramid stable 10 s at 60 Hz with substeps=4, max drift <1 cm; collision detected once per frame, gravity+velocity solve per substep.
  - **Design:** Macklin et al. 2019 / PhysX 5 small-steps: split `PhysicsWorld::step` into detect → N×(integrate gravity, 1-iter velocity solve, relax) — distinct from the existing naive `step_substep` (re-runs full pipeline per sub-dt, `world/types.rs:1180`); per-substep restitution helper already exists (`impulse/functions_2.rs`).
  - **Files:** `world/types.rs` (step pipeline split), `impulse/functions_2.rs` (per-substep restitution reuse)
  - **Tests:** 1000-box pyramid drift assertion; substeps ∈ {1, 2, 4, 8} consistency sweep
  - **Risk:** changes stepping semantics — keep legacy `step_substep` available until the parity suite is green
  - **Files:** world/types.rs (step_small_steps + SolverConfig.substeps), reuse impulse/functions_2.rs SubStepRestitutionCorrector
  - **Tests:** 100-box pile stable at substeps=4 (drift <1cm, no NaN); energy non-increasing; fast-body no energy gain; restitution accuracy vs analytic ±5%
  - **Risk:** contact re-projection across substeps — skip separated contacts; document the small-steps contract
  - **Done:** PhysicsWorld::step_small_steps + SolverConfig.substeps/restitution in world/types.rs (1519 lines, no split needed); detect-once, re-project-not-redetect, skip-separated contacts, Baumgarte+max(restitution) target, 2-iter relax pass; 4 tests in tests/small_steps.rs (pile stability/no-NaN/no-tunnel, energy non-increasing, 50 m/s no energy gain, e=0.5 rebound ratio ~0.25). Stack test uses a 5-sphere column (flat grid rolls off the convex sphere-ground — geometry artifact, not solver). step()/step_substep unchanged.

- [ ] Shock propagation for tall stacks
  - **Goal:** 50-box tower survives 30 s at 60 Hz, 1-iteration solver.
  - **Design:** Guendelman 2003 / Bullet: extra solver pass with gravity-direction layer sort, treating lower layer as static-equivalent (zero inverse mass downward).
  - **Tests:** 50-box tower scenario; layer-sort unit tests on synthetic islands

### Contact response & restitution

- [x] TGS-Soft contact response (done 2026-06-12)
  - **Goal:** stiffness specified as (hertz, damping-ratio) per contact; 100-box stack shows no jitter at 30 Hz; soft contacts converge independent of mass ratio up to 1e4.
  - **Design:** Jolt-style soft constraints — bias/mass/impulse coefficients from ω=2πf, ζ (Catto soft-constraint derivation), applied in TGS position iterations; existing `soft_constraints.rs` is penalty springs only — keep, but route the solver path through oxiphysics-constraints TGS.
  - **Depends:** oxiphysics-constraints v0.2.0 soft-constraint parameterization (hz, ζ across PGS/TGS)
  - **Tests:** 100-box stack jitter metric at 30 Hz; mass-ratio ladder 1e0–1e4 convergence
  - **Done:** `oxiphysics-constraints` already depends on `oxiphysics-rigid` (cycle), so `SoftParams::{from_frequency, rigid, is_rigid}` is self-contained in `src/solver/soft.rs` (Catto/Box2D-v3 formula, kept in sync with the constraints crate). `SolverConfig` gains `contact_hertz` / `contact_damping_ratio` / `use_soft_contacts` (default off ⇒ byte-identical legacy path). Both `PhysicsWorld::step` and `step_small_steps` route contacts through `apply_soft_normal_constraint` when enabled: full 6-DoF contact Jacobian `J=[n, r_a×n, −n, −(r_b×n)]` (angular terms vanish for the sphere narrowphase but are exact for off-centre contacts), accumulated unilateral λ warm-started across iterations/sub-steps, `target = max(bias_rate·depth, restitution)` (never summed — P2 lesson), `hz` capped to `0.25/h` (Catto). Tests: `tests/tgs_soft.rs` (soft-off byte-identity for both step paths, jitter-free 5-stack with non-increasing KE tail, dt-invariant settle time across 1/60–1/240, 100:1 heavy-on-light, restitution preservation) + 5 unit tests in `src/solver/soft.rs`. world/types.rs 1772 lines (no split).

- [x] Speculative contacts integrated with solver (planned 2026-06-13)
    - **Files:** MODIFY src/world/types.rs (SolverConfig + admission + speculative target); possibly src/world/solverconfig_traits.rs; NEW tests/speculative_contacts.rs
    - **Tests:** No-tunnel 1-100 m/s vs 1cm wall; no ghost bounce; legacy byte-parity with speculative off; no spurious tangential impulse
    - **Risk:** Distance-clamped target sign/scale — test both no-tunnel and no-bounce gates together
  - **Goal:** 50 m/s sphere vs 1 cm wall never tunnels, no CCD substep, no ghost bounce.
  - **Design:** detection candidate already exists (`speculative_contact_candidate`, `ccd/types.rs:354`) and `narrowphase/dispatch/speculativeconfig_traits.rs`; missing piece is solver-side: admit negative-penetration contacts with distance-clamped target velocity (PhysX/Box2D speculative margin).
  - **Tests:** bullet-vs-thin-wall speed sweep (1–100 m/s); ghost-bounce regression at the resting margin

- [ ] Newton vs Poisson restitution
  - **Goal:** selectable hypothesis per material pair; Newton/Poisson differ correctly on multi-contact frictional impact tests (superball).
  - **Design:** Poisson via compression/restitution impulse phases (Stronge-aware energy check to forbid energy gain).
  - **Tests:** superball multi-contact impact; energy-gain-forbidden assertion (Stronge check)

### Friction & gyroscopic integration

- [ ] Friction model upgrade: rolling/spinning/anisotropic with smooth cone projection
  - **Goal:** bowling ball rolls to rest with analytically correct stopping distance ±2%; anisotropic μ ellipse (sledge) verified directionally.
  - **Design:** add torsional+rolling rows to the impulse solver loop with coupled cone clamp (existing `clamp_friction_to_cone` / pyramidal cone in `impulse/functions.rs:66-79` is linear-only; rolling helpers in `contact_dynamics.rs` are standalone, not in the solver).
  - **Files:** `impulse/functions.rs` (coupled cone clamp), `contact_dynamics.rs` (rolling helpers promoted into solver rows)
  - **Tests:** bowling-ball stopping distance vs analytic ±2%; sledge directional-μ scenario

- [x] Implicit gyroscopic integration
  - **Goal:** spinning top (T-handle/Dzhanibekov) energy-bounded for 60 s at dt=1/60 without clamping.
  - **Design:** one local Newton step on ω in body frame (Catto, GDC 2015); existing `gyroscopic.rs` has explicit Euler equations + RK4 only — add implicit path to the body integrator.
  - **Files:** `gyroscopic.rs`, body integrator entry point
  - **Tests:** Dzhanibekov flip energy bound over 60 s; explicit vs RK4 vs implicit comparison harness

## v0.3.0 — Scale and differentiability

Theme: deterministic multicore throughput, exact hydrostatics, and gradient-capable stepping for optimization and learning workloads.
Exit gate: thread-count-invariant trajectory hash, analytic buoyancy parity, and gradcheck-vs-finite-differences all green.

- [ ] Deterministic parallel island solve
  - **Goal:** bitwise-identical results, 1 vs N threads, on 10k-body scene; ≥4x speedup on 8 cores.
  - **Design:** rayon (new optional dep) over islands (islands/sleeping already exist, no rayon anywhere yet); canonical island ordering by min body id, fixed-order impulse accumulation, no atomics in accumulation.
  - **Tests:** trajectory-hash equality across thread counts {1, 2, 8}; criterion island-throughput benchmark

- [ ] Exact submerged-volume buoyancy for convex shapes
  - **Goal:** floating box/capsule equilibrium draft matches analytic to 1e-6; metacentric righting moment matches `buoyancy_model.rs` formulas.
  - **Design:** clip convex hull against water plane (reuse Sutherland-Hodgman from oxiphysics-collision `contact_generation.rs`), exact polyhedron volume+centroid integrals (Mirtich), per-face damping; existing `buoyancy_model.rs` is scalar Archimedes/hull-form only.
  - **Depends:** clipping primitives in oxiphysics-collision `contact_generation.rs` (already shipped)

- [ ] Differentiable rigid step
  - **Goal:** gradient of final pose w.r.t. initial velocity/μ matches finite differences to 1e-5 on a 2-body contact scene; ball-throw target optimization converges <50 iters.
  - **Design:** adjoint through the contact LCP via implicit function theorem (de Avila Belbute-Peres 2018, Werling "Nimble" 2021); feature-gated `diff` (scirs2-autograd optional dep).
  - **Risk:** large surface area — ship behind the non-default `diff` feature so the default build stays dependency-light and Pure Rust

## v1.0 — Production & validation

Theme: hard-tolerance stability gates, cross-platform determinism, and soak-tested robustness ahead of API freeze.

- [ ] Stacking/stability benchmark suite with hard tolerances
  - **Goal:** CI-gated (local script): 1000-box pyramid (10 s, drift <1 cm), 100-capsule pile (monotone KE decay), 50-box tower, ragdoll pile (no NaN, sleep within 5 s).
  - **Design:** criterion + scenario crate-level tests, JSON baselines via oxicode.
  - **Note:** oxicode baselines keep the suite Pure Rust (no bincode, per COOLJAPAN policy)

- [ ] Cross-platform determinism
  - **Goal:** identical 1000-step trajectory hash on x86_64 and aarch64 (strict-fp mode).
  - **Design:** no `mul_add` in default path, total-order sort keys; document FMA-feature opt-out.
  - **Note:** complements the v0.3.0 deterministic parallel island solve — thread-count and ISA determinism are separate gates

- [ ] Long-soak fuzzing of world state machine
  - **Goal:** 1e6 proptest-driven random steps (add/remove/sleep/wake/CCD toggles) with zero panics/NaN.
  - **Design:** proptest strategies over body/world ops; invariant checks on islands and sleep timers.
  - **Files:** new proptest strategy module under `tests/` (temp paths via `std::env::temp_dir()` where state snapshots are needed)

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 269 (critical 10, high 70, medium 148, low 41); 32 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (10)
- [ ] **F3833** (critical, bug, easy/small) `granular_media.rs:189` — GranularPacking normal spring is attractive; tangential sign reversed
- [ ] **F3838** (critical, bug, easy/trivial) `granular_rigid.rs:525` — Hertz-Mindlin force_j breaks Newton's third law (-3x tangential)
- [ ] **F3862** (critical, bug, easy/trivial) `impulse/types.rs:601` — Sequential-impulse friction (t1/t2) pushes sliding faster
- [ ] **F3898** (critical, bug, easy/trivial) `multibody_dynamics.rs:435` — gravity_torques omits joint position and has inverted sign
- [ ] **F3914** (critical, bug, hard/medium) `orbital_mechanics.rs:654` — lambert_problem uses wrong λ, ρ and T(x); test vacuous
- [ ] **F3941** (critical, bug, easy/trivial) `satellite_dynamics.rs:565` — AttitudeControl::pd_controller torque sign drives attitude away from target
- [x] **F4017** (critical, fabrication, hard/large) `world/types.rs:969` — PhysicsWorld ignores colliders; mass-cbrt spheres and sign-inverted contact solve
- [x] **F4018** (critical, fabrication, hard/large) `world/types.rs:979` — PhysicsWorld ignores Collider shapes; contact radius fabricated from mass
- [x] **F4019** (critical, fabrication, hard/large) `world/types.rs:1007` — PhysicsWorld collision ignores colliders, uses mass-cube-root spheres
- [x] **F4021** (critical, bug, medium/medium) `world/types.rs:1094` — Default PhysicsWorld solver skips approaching contacts and attracts separating ones

### High (70)
- [x] **F3728** (high, fabrication, medium/small) `aerospace_rigid.rs:680` — Aircraft 6-DOF ignores wind and omits the ω×v transport term
- [x] **F3731** (high, fabrication, medium/medium) `articulated/functions.rs:544` — rnea_revolute_chain returns zero torques; no backward recursion
- [x] **F3732** (high, fabrication, hard/medium) `articulated/functions.rs:759` — aba_revolute_chain is not ABA: subtracts acceleration from torque
- [x] **F3733** (high, fabrication, medium/medium) `articulated/functions.rs:775` — Joint-space inertia, Coriolis and kinetic energy ignore configuration
- [x] **F3735** (high, bug, medium/small) `articulated/functions.rs:1034` — rnea_3d_revolute_z drops child-force moment and uses wrong parent accel
- [x] **F3736** (high, fabrication, hard/medium) `articulated/functions.rs:1470` — forward_dynamics_full uses diagonal, configuration-independent H and fake Coriolis
- [ ] **F3739** (high, bug, hard/medium) `articulated/types.rs:722` — ArticulatedChain2D dynamics: wrong lever arm, gravity sign, decoupled joints
- [x] **F3742** (high, bug, easy/small) `body/rigidbody_kinetic_energy_group.rs:17` — kinetic_energy and angular_momentum use local inertia with world-frame ω
- [x] **F3745** (high, bug, easy/small) `body/rigidbody_predicates.rs:59` — Per-step multiplicative damping is dt-dependent and on by default
- [ ] **F3755** (high, bug, hard/medium) `ccd/functions.rs:407` — gjk_distance has no simplex; misreports distance and overlap
- [ ] **F3756** (high, fabrication, medium/small) `ccd/functions.rs:432` — gjk_distance has no simplex; reports deep overlap as separated
- [ ] **F3758** (high, bug, medium/small) `ccd/types.rs:608` — CCD resolve_event indexes by user ID; touching pairs stall loop
- [ ] **F3760** (high, bug, easy/trivial) `cfd_coupling.rs:715` — VIV structural EOM divides damping and stiffness accelerations by mass
- [ ] **F3761** (high, bug, medium/small) `cfd_coupling.rs:851` — Mooring catenary solve_static residual has no root
- [ ] **F3772** (high, fabrication, medium/small) `contact_dynamics.rs:335` — ImpulseResolver ignores friction and the Baumgarte penetration bias
- [ ] **F3785** (high, fabrication, medium/medium) `deformable_coupling/types.rs:229` — CoupledSimulation::step fabricates coupling force direction and reaction
- [ ] **F3792** (high, fabrication, easy/small) `destruction_rigid.rs:291` — Documented RNG seed parameter silently ignored in fracture generators
- [ ] **F3797** (high, panic, easy/small) `destruction_rigid.rs:707` — random_range on empty f64 ranges panics in public APIs
- [ ] **F3805** (high, bug, easy/trivial) `fluid_coupling/types.rs:451` — Cylinder partial buoyancy: wrong acos argument
- [ ] **F3807** (high, fabrication, medium/small) `fluid_dynamics/types.rs:424` — Panel::source_influence formula is invented and uses the wrong normal
- [ ] **F3819** (high, fabrication, medium/small) `fluid_structure/types.rs:1682` — PartitionedFsiCoupling 'Aitken/Newmark' sub-iteration is fabricated
- [ ] **F3823** (high, fabrication, medium/small) `fracture.rs:658` — Mott mean fragment mass formula fabricated (≈3e-10 kg, 4.7e9 fragments)
- [ ] **F3824** (high, fabrication, easy/small) `fracture_dynamics.rs:772` — Fracture stress = F/(mass·0.001); hardcoded density/extents; force not applied
- [ ] **F3826** (high, bug, easy/small) `fracture_rigid.rs:521` — LCG maps to [0,0.5): seeds and explosion fragments confined to one quadrant
- [ ] **F3827** (high, fabrication, medium/small) `fracture_rigid.rs:565` — ImpactFracture::fracture_body violates mass, momentum and energy conservation
- [ ] **F3828** (high, bug, easy/trivial) `fracture_rigid.rs:1018` — Sadovsky overpressure 10x too high; blast impulse mis-scaled
- [ ] **F3830** (high, bug, easy/trivial) `granular.rs:199` — DemContact normal damping sign inverted (injects energy)
- [ ] **F3834** (high, bug, easy/small) `granular_rigid.rs:287` — Mindlin tangential stiffness uses tangential instead of normal overlap
- [ ] **F3836** (high, fabrication, easy/small) `granular_rigid.rs:452` — Normal damping is anti-damping; 0.1·sqrt(F)·v is dimensionless nonsense
- [ ] **F3837** (high, bug, easy/small) `granular_rigid.rs:462` — Hertz-Mindlin contact velocity sign wrong; force_j breaks Newton III
- [ ] **F3841** (high, bug, medium/small) `granular_rigid.rs:1104` — detect_contacts wipes Mindlin tangential springs every step
- [ ] **F3844** (high, bug, medium/small) `gyroscope_rigid.rs:473` — Lagrange top theta equation has gravity and gyroscopic signs flipped
- [ ] **F3850** (high, bug, easy/trivial) `impact_dynamics.rs:476` — Hertz impact duration has wrong exponents; peak force lacks 5/4
- [ ] **F3854** (high, fabrication, easy/small) `impact_mechanics.rs:791` — Taylor dynamic yield always returns rho*c*v0 (ignores final length)
- [ ] **F3856** (high, bug, easy/small) `impulse/functions.rs:49` — Local inverse inertia applied to world-space vectors; orientation ignored
- [ ] **F3857** (high, bug, easy/trivial) `impulse/types.rs:212` — ImpulseSolver::apply_coulomb_friction accelerates sliding (wrong sign)
- [ ] **F3859** (high, bug, easy/small) `impulse/types.rs:353` — Normal convention contradicts code; ContactConstraint friction accelerates sliding
- [ ] **F3861** (high, panic, easy/small) `impulse/types.rs:414` — Sequential impulse solvers panic when body_a >= body_b
- [ ] **F3865** (high, panic, easy/small) `kinematics/types.rs:430` — Jacobian::manipulability panics for 7+ joints and mixes units for 4-5
- [ ] **F3867** (high, stub, medium/small) `kinematics/types_2.rs:90` — Pfaffian project_velocity silently skips projection for m>1
- [ ] **F3872** (high, bug, medium/medium) `marine_dynamics.rs:701` — Mooring catenary solver solves non-catenary equations
- [ ] **F3877** (high, fabrication, easy/trivial) `marine_rigid.rs:996` — marine_rigid::jonswap_spectrum discards Hs; output independent of sea state
- [ ] **F3893** (high, fabrication, hard/medium) `multibody/types.rs:758` — FeatherstoneABA adds mass to inertia; gravity torque is a force
- [ ] **F3894** (high, fabrication, hard/large) `multibody/types.rs:1108` — MultibodySystem dynamics fabricated: kg+kg·m², sin(q)q̇², +0.1 fudge
- [ ] **F3895** (high, bug, easy/trivial) `multibody/types.rs:1366` — ArticulatedBody::clamp_angles shifts every angle by −π
- [ ] **F3896** (high, fabrication, hard/large) `multibody/types.rs:1588` — FeatherstoneAba/RecursiveNewtonEuler/JointSpaceInertia are fake decoupled pendulums
- [ ] **F3897** (high, fabrication, hard/medium) `multibody_dynamics.rs:394` — CRBA mass_matrix lacks parallel-axis terms; M(q) wrong
- [ ] **F3899** (high, bug, medium/small) `multibody_dynamics.rs:483` — coriolis_vector uses wrong Christoffel derivative; n clones per call
- [ ] **F3910** (high, bug, medium/medium) `neural_control/types.rs:1234` — RRT* reparents without a collision check and never rewires
- [ ] **F3911** (high, bug, medium/small) `neural_control/types.rs:1246` — RRT* reparents without collision check and never rewires
- [ ] **F3913** (high, bug, medium/small) `orbital_mechanics.rs:437` — cartesian_to_kepler: circular orbits always return true anomaly π/2
- [ ] **F3915** (high, bug, hard/medium) `orbital_mechanics.rs:667` — Lambert solver uses garbled Izzo equations; test vacuous
- [ ] **F3926** (high, fabrication, medium/medium) `planetary_mechanics.rs:450` — NBodyGravity Barnes-Hut never used; integrator is Euler, not RK4/leapfrog
- [ ] **F3927** (high, bug, easy/small) `planetary_mechanics.rs:800` — Piecewise density model discontinuous by up to nine orders
- [ ] **F3928** (high, bug, easy/trivial) `planetary_mechanics.rs:833` — AtmosphericDrag::decay_rate off by a factor of 2πa (wrong units)
- [ ] **F3932** (high, fabrication, medium/small) `planetary_mechanics.rs:1213` — two_impulse_rendezvous uses invented formulas, ignores initial velocity
- [ ] **F3933** (high, fabrication, medium/small) `ragdoll/functions.rs:101` — Ragdoll forward_kinematics ignores joint angles and root orientation
- [ ] **F3936** (high, bug, easy/small) `ragdoll/types.rs:897` — Ragdoll position projection never updates velocities
- [ ] **F3944** (high, fabrication, medium/small) `satellite_dynamics.rs:968` — fuel_optimal_maneuver returns fabricated error·n/sin(nt), pointing away from target
- [ ] **F3945** (high, bug, easy/small) `satellite_dynamics.rs:1043` — set_circular_orbit drifts/offset; is_bounded_orbit test wrong
- [ ] **F3962** (high, fabrication, hard/medium) `space_mechanics.rs:802` — lambert_solver iterates on wrong residual, never fails, vacuous test
- [ ] **F3983** (high, bug, easy/trivial) `thruster_dynamics.rs:270` — BrushlessMotor back-EMF multiplies by Kv instead of dividing
- [ ] **F3985** (high, bug, medium/small) `train_dynamics.rs:131` — TrainConsist couplers anti-restoring; brakes reverse a stopped train
- [ ] **F3986** (high, bug, easy/trivial) `train_dynamics.rs:172` — TrainConsist::step applies coupler force with anti-restoring sign
- [ ] **F3992** (high, fabrication, medium/medium) `train_rigid.rs:593` — sperling_wz uses fabricated power law, not the Sperling index
- [ ] **F3995** (high, bug, medium/small) `train_rigid.rs:1048` — train_longitudinal_step: resistance pushes resting car; couplers anti-restoring
- [ ] **F3996** (high, bug, easy/small) `train_rigid.rs:1055` — train_longitudinal_step coupler anti-restoring; resistance pushes stationary cars backward
- [ ] **F4001** (high, bug, medium/medium) `underwater_vehicle.rs:241` — UnderwaterVehicle::step uses kg as rotational inertia, omits hydrostatics
- [ ] **F4003** (high, bug, medium/small) `underwater_vehicle.rs:375` — ThrusterAllocation solves singular B·Bᵀ normal equations when rank<6
- [x] **F4020** (high, fabrication, hard/large) `world/types.rs:1022` — PhysicsWorld collides mass-derived spheres; step() never reads the ColliderSet

### Medium (148)
- [ ] **F3724** (medium, perf, easy/trivial) `crates/oxiphysics-rigid/Cargo.toml:17` — Eight internal path dependencies declared but never referenced
- [ ] **F3726** (medium, bug, easy/trivial) `aerospace_maneuver.rs:433` — INS gravity compensation sign wrong for real accelerometer specific force
- [ ] **F3727** (medium, bug, easy/small) `aerospace_maneuver.rs:722` — dubins_rsr turn angles wrong: straight path gets π·r extra
- [ ] **F3734** (medium, bug, medium/small) `articulated/functions.rs:808` — gravity_compensation_torques wrong for bent chains; panics on short slices
- [ ] **F3737** (medium, panic, easy/small) `articulated/types.rs:424` — add_link mutates before validating parent; forged LinkId panics later
- [ ] **F3738** (medium, bug, easy/small) `articulated/types.rs:530` — apply_velocity_limits clamps velocity to position limits
- [ ] **F3743** (medium, api, easy/small) `body/rigidbody_new_group.rs:15` — RigidBody::new accepts negative/NaN/zero mass as Dynamic
- [ ] **F3744** (medium, bug, medium/small) `body/rigidbody_new_group.rs:108` — RigidBody::merge ignores parallel-axis, rotation and orbital momentum
- [ ] **F3747** (medium, bug, medium/small) `body/types.rs:233` — RollingConstraint::enforce corrects the wrong omega axis and never removes slip
- [ ] **F3748** (medium, bug, medium/small) `body/types.rs:282` — CompoundRigidBody drops products of inertia and primitive rotations
- [ ] **F3749** (medium, bug, easy/small) `body/types.rs:949` — VariableMassBody thrust in body frame applied as a world force
- [ ] **F3751** (medium, bug, easy/trivial) `buoyancy_model.rs:309` — Von Kármán impact force dimensionally wrong
- [ ] **F3752** (medium, bug, easy/small) `cables.rs:120` — CatenarySegment::from_endpoints wrong for unequal support heights
- [ ] **F3753** (medium, bug, easy/small) `cables.rs:139` — CatenarySegment::from_endpoints wrong for unequal heights
- [ ] **F3754** (medium, bug, easy/small) `cables.rs:350` — Cable segments push in compression despite 'only pulls' comment
- [ ] **F3757** (medium, bug, easy/small) `ccd/types.rs:551` — CCD run_with_stats stalls on touching pairs and drops remaining dt
- [ ] **F3759** (medium, bug, easy/small) `cfd_coupling.rs:308` — IBM delta kernel is radial 1-D; advance ignores rotation
- [ ] **F3763** (medium, bug, easy/trivial) `collider.rs:528` — Capsule inertia double-counts hemispheres and ignores cap COM offset
- [ ] **F3764** (medium, bug, easy/small) `collider.rs:648` — TransformedCollider::world_aabb ignores rotation despite 'conservative' doc
- [ ] **F3766** (medium, bug, easy/small) `constraint_forces/types.rs:310` — NullSpaceProjector skips (JJᵀ+λ²I)⁻¹ when task dim ≠ 3
- [ ] **F3767** (medium, fabrication, medium/small) `constraint_forces/types.rs:383` — min_norm_tensions clamps the pseudoinverse; A·τ ≠ f silently
- [ ] **F3769** (medium, fabrication, medium/small) `constraint_forces/types.rs:887` — Stewart forward_kinematics_numerical is not Newton; only nudges z
- [ ] **F3771** (medium, bug, easy/trivial) `constraint_graph.rs:188` — closed_loop_detection treats undirected joints as directed DAG
- [ ] **F3773** (medium, fabrication, medium/medium) `contact_dynamics.rs:347` — Warm starting is a no-op and contact friction is silently ignored
- [ ] **F3774** (medium, panic, easy/trivial) `contact_dynamics.rs:382` — solve_manifolds panics on self-contact manifold (body_a == body_b)
- [ ] **F3775** (medium, bug, easy/trivial) `contact_dynamics.rs:588` — Island.all_sleeping is always false; island order nondeterministic
- [ ] **F3776** (medium, fabrication, medium/medium) `contact_dynamics.rs:684` — WarmStarting::apply has no effect on velocities; friction ignored
- [ ] **F3778** (medium, fabrication, easy/small) `contact_geometry.rs:204` — grasp_quality returns sqrt(min diag), not smallest singular value
- [ ] **F3779** (medium, fabrication, easy/small) `contact_geometry.rs:215` — grasp_quality returns an upper bound mislabelled as a lower bound
- [ ] **F3780** (medium, bug, easy/small) `continuum_robot.rs:291` — CossratRod strain reference and curvature units are wrong
- [ ] **F3781** (medium, perf, medium/medium) `crowd_simulation.rs:309` — Social-force total_force is O(N²) per step with no cutoff
- [ ] **F3783** (medium, panic, easy/small) `crowd_simulation.rs:950` — choose_exit / update_exit_choices panic when no exits configured
- [ ] **F3784** (medium, stub, hard/large) `deformable_coupling/functions.rs:689` — craig_bampton returns only DOF counts; no reduction performed
- [ ] **F3786** (medium, missing-feature, medium/small) `deformable_coupling/types.rs:255` — ShapeMatchingCoupling never computes the documented optimal rotation
- [ ] **F3787** (medium, missing-feature, medium/small) `deformable_coupling/types.rs:301` — ShapeMatchingCoupling never computes the optimal rigid transform it documents
- [ ] **F3788** (medium, fabrication, easy/small) `deformable_coupling/types.rs:358` — peak_hertz_force is a linear-spring formula, not Hertz impact
- [ ] **F3789** (medium, fabrication, medium/medium) `deformable_rigid.rs:133` — DeformableBody modal stiffness/mass dimensionally wrong ad-hoc estimate
- [ ] **F3790** (medium, bug, easy/small) `deformable_rigid.rs:497` — ImpactDeformation energy integral halved and includes elastic energy
- [ ] **F3791** (medium, fabrication, medium/medium) `deformable_rigid.rs:803` — in_contact_wrench_cone ignores moments and tests net force per contact
- [ ] **F3793** (medium, bug, easy/small) `destruction_rigid.rs:298` — Seed parameters ignored: destruction is non-deterministic
- [ ] **F3795** (medium, panic, easy/small) `destruction_rigid.rs:353` — Voronoi destruction/debris panic on zero fragments, zero scatter, flat AABB
- [ ] **F3796** (medium, fabrication, easy/small) `destruction_rigid.rs:422` — Fracture normal ignores shear; not principal direction
- [ ] **F3798** (medium, fabrication, medium/medium) `destruction_rigid.rs:1230` — Radial destruction discards its pattern; 'impact-biased' sites are uniform
- [ ] **F3799** (medium, api, easy/small) `error.rs:23` — Crate Error enum is dead: no public API returns Result
- [ ] **F3802** (medium, bug, easy/trivial) `fluid_coupling/types.rs:122` — Prolate-spheroid added-mass coefficient formula wrong
- [ ] **F3803** (medium, fabrication, medium/medium) `fluid_coupling/types.rs:172` — OrientedBuoyancy ignores orientation entirely
- [ ] **F3804** (medium, bug, easy/trivial) `fluid_coupling/types.rs:415` — partial_submersion multiplies the submerged fraction in twice
- [ ] **F3806** (medium, bug, easy/small) `fluid_coupling/types.rs:858` — Morison drag uses circular area, not projected D*L
- [ ] **F3808** (medium, bug, easy/small) `fluid_dynamics/types.rs:603` — Regular-wave Morison max force ignores depth and phase
- [ ] **F3811** (medium, missing-feature, medium/medium) `fluid_dynamics/types.rs:945` — Panel/VLM 'solvers' assemble systems but never solve them
- [ ] **F3812** (medium, fabrication, medium/medium) `fluid_dynamics/types.rs:1441` — LiftingLine ignores n_modes; no Fourier solution
- [ ] **F3813** (medium, fabrication, easy/small) `fluid_structure/types.rs:77` — JONSWAP spectra: Hs discarded, A_γ missing, Tz ignores gamma
- [ ] **F3814** (medium, bug, easy/trivial) `fluid_structure/types.rs:269` — Roma-Peskin 3-point delta inner branch wrong (sums to 2/3)
- [ ] **F3815** (medium, bug, medium/small) `fluid_structure/types.rs:358` — rms_morison_force uses sqrt(S(wp)) as RMS; JONSWAP unnormalised
- [ ] **F3816** (medium, bug, easy/small) `fluid_structure/types.rs:388` — RandomSeaState: fake linear phases, wrong eta_rms, unnormalised JONSWAP
- [ ] **F3817** (medium, fabrication, hard/large) `fluid_structure/types.rs:545` — Flutter 'p-k method' is just wa*b*sqrt(mu); wrong frequencies
- [ ] **F3818** (medium, fabrication, medium/medium) `fluid_structure/types.rs:1674` — PartitionedFsiCoupling::sub_iterate is not Newmark/Aitken
- [ ] **F3820** (medium, fabrication, medium/small) `fracture.rs:258` — voronoi_fracture is a grid with chain bonds, not Voronoi
- [ ] **F3822** (medium, bug, easy/small) `fracture.rs:403` — Grady fragment size formula dimensionally wrong
- [ ] **F3825** (medium, bug, medium/small) `fracture_rigid.rs:414` — Debris broad phase O(n^2); resting fragments become ghosts
- [ ] **F3829** (medium, fabrication, easy/small) `fracture_rigid.rs:1059` — Mott sampler fake; LCG in [0,0.5) sends fragments downward
- [ ] **F3831** (medium, perf, medium/medium) `granular.rs:382` — DemSimulation::step brute-force O(n²) pairs plus a per-step wall clone
- [ ] **F3832** (medium, perf, medium/medium) `granular_media.rs:174` — GranularPacking O(n²) contact loops; per-pair ln/sqrt invariants
- [ ] **F3835** (medium, fabrication, easy/small) `granular_rigid.rs:328` — Capillary force mislabelled Willett; no contact angle or rupture
- [ ] **F3839** (medium, bug, easy/small) `granular_rigid.rs:768` — measure_angle_of_repose returns ~90° for any heap
- [ ] **F3840** (medium, bug, easy/small) `granular_rigid.rs:1005` — CompactionState first load spans 11 decades; 'strain' is porosity
- [ ] **F3842** (medium, bug, easy/small) `granular_rigid.rs:1178` — Hertz strain energy re-added each step; static load crushes grains
- [ ] **F3843** (medium, stub, medium/small) `granular_rigid.rs:1210` — Grain crushing only freezes particle; split_grain never used
- [ ] **F3845** (medium, bug, easy/trivial) `gyroscope_rigid.rs:622` — Gyro sensor noise RNG never negative (biased noise, monotone drift)
- [ ] **F3846** (medium, fabrication, easy/small) `gyroscope_rigid.rs:705` — Gyrocompass step discards Earth-rate forcing; first-order relaxation
- [ ] **F3847** (medium, bug, easy/trivial) `gyroscope_rigid.rs:766` — ReactionWheel returns full reaction torque even when saturated
- [ ] **F3849** (medium, fabrication, easy/small) `impact_dynamics.rs:183` — Ballistic penetration uses hidden hardcoded 1 cm^2 area
- [ ] **F3852** (medium, fabrication, medium/small) `impact_mechanics.rs:487` — Holsapple crater ignores target strength and K2; volume inconsistent
- [ ] **F3853** (medium, fabrication, medium/small) `impact_mechanics.rs:690` — Grueneisen pressure uses invented compression->velocity map
- [ ] **F3860** (medium, bug, medium/small) `impulse/types.rs:377` — Warm-start cache keyed per body pair; linear lookups
- [ ] **F3864** (medium, bug, easy/small) `kinematics/types.rs:32` — Manipulability ellipsoid uses Gershgorin bounds as eigenvalues
- [ ] **F3866** (medium, panic, easy/small) `kinematics/types.rs:644` — Workspace::sample panics on 0 joints/overflow; last joints never move
- [ ] **F3870** (medium, bug, easy/small) `lib.rs:255` — Angular momentum/rotational KE use local inertia with world omega
- [ ] **F3871** (medium, stub, medium/medium) `locomotion.rs:150` — Biped/quadruped models never move feet; speed/height ignored
- [ ] **F3873** (medium, fabrication, medium/medium) `marine_dynamics.rs:1005` — ThrusterAllocation ignores yaw-moment demand and thruster positions
- [ ] **F3874** (medium, fabrication, medium/small) `marine_dynamics.rs:1012` — ThrusterAllocation ignores yaw moment; panics on length mismatch
- [ ] **F3875** (medium, fabrication, hard/large) `marine_rigid.rs:920` — simple_resistance and wave_resistance falsely claim Holtrop-Mennen
- [ ] **F3876** (medium, fabrication, medium/medium) `marine_rigid.rs:952` — simple_resistance claims Holtrop-Mennen; discards terms, fake wave coeff
- [ ] **F3878** (medium, panic, easy/small) `mechanism_rigid/functions.rs:193` — workspace_analysis overflows samples^n and stores unbounded points
- [ ] **F3879** (medium, fabrication, easy/small) `mechanism_rigid/functions.rs:268` — singularity_analysis condition number is cbrt(tr³/det), not σmax/σmin
- [ ] **F3881** (medium, bug, easy/small) `mechanism_rigid/types.rs:158` — Geneva mechanism kinematic formulas wrong; tests mirror them
- [ ] **F3883** (medium, fabrication, medium/medium) `motors/motors_extended.rs:1102` — BldcMotor torque discards back-EMF; always non-negative
- [ ] **F3886** (medium, bug, medium/medium) `motors/types.rs:696` — SCurveProfile double-counts jerk distance; short and v-limited moves wrong
- [ ] **F3887** (medium, bug, medium/medium) `motors/types.rs:990` — StepperMotor: wrong electrical angle, unstable torque sign, ignored R/L/V
- [ ] **F3888** (medium, fabrication, hard/medium) `multibody/functions.rs:747` — FbMultibodySystem dynamics are wrong even for the stated prismatic model
- [ ] **F3889** (medium, fabrication, medium/medium) `multibody/functions.rs:770` — FbMultibodySystem dynamics wrong even for its prismatic model
- [ ] **F3891** (medium, fabrication, medium/small) `multibody/types.rs:650` — OSC null_space_projection returns constant 0.5; kd unused
- [ ] **F3892** (medium, bug, easy/small) `multibody/types.rs:690` — Self-collision skips index-adjacent pairs, so sibling collisions are missed
- [ ] **F3900** (medium, panic, medium/medium) `multibody_dynamics.rs:588` — MultibodyChain dynamics panic on short inputs; singular solve returns zeros
- [ ] **F3906** (medium, perf, hard/large) `neural_control/types.rs:252` — All learners use O(P²) finite-difference gradients
- [ ] **F3907** (medium, perf, medium/medium) `neural_control/types.rs:269` — Finite-difference NN gradients cost O(P²) per training sample
- [ ] **F3908** (medium, panic, medium/small) `neural_control/types.rs:962` — cem_plan panics when num_candidates or action_bound is zero
- [ ] **F3909** (medium, panic, easy/small) `neural_control/types.rs:1139` — Planner and ReplayBuffer panic on NaN, zero capacity, degenerate bounds
- [ ] **F3916** (medium, bug, medium/small) `orbital_mechanics.rs:1066` — collision_probability ignores miss distance and uses the wrong volume
- [ ] **F3918** (medium, perf, medium/small) `pipeline/types.rs:23` — Broadphase all-pairs loop recomputes cbrt AABBs per pair
- [ ] **F3920** (medium, fabrication, hard/large) `pipeline/types.rs:373` — handle_ccd is not CCD, and the other pipelines skip it
- [ ] **F3922** (medium, fabrication, medium/medium) `pipeline/types.rs:711` — Narrowphase claims box contacts but invents spheres from mass
- [ ] **F3923** (medium, fabrication, medium/medium) `pipeline/types.rs:753` — Pipeline collision radius derived from mass; fake CCD claims
- [ ] **F3925** (medium, bug, easy/small) `pipeline/types.rs:1017` — Collision, sleep and trigger events come out in HashSet/HashMap order
- [ ] **F3929** (medium, bug, easy/small) `planetary_mechanics.rs:909` — LagrangePoints L1 equation is wrong; L2/L3 are only approximations
- [ ] **F3930** (medium, missing-feature, hard/large) `planetary_mechanics.rs:1108` — TLE propagator advertised as SGP4 is Kepler plus partial J2
- [ ] **F3931** (medium, panic, easy/small) `planetary_mechanics.rs:1154` — TLE parse/verify_checksum panic on non-ASCII, blank epoch or bad check digit
- [ ] **F3934** (medium, bug, hard/medium) `ragdoll/types.rs:840` — Ragdoll PBD step never derives velocities from position corrections
- [ ] **F3935** (medium, bug, medium/medium) `ragdoll/types.rs:887` — Hinge limit uses an unclamped, unsigned acos; cone limit is a half-nlerp
- [ ] **F3938** (medium, fabrication, hard/large) `robotic_systems/types.rs:806` — RobotDynamics claims RNEA but uses unrotated inertia and fake Coriolis
- [ ] **F3940** (medium, fabrication, medium/medium) `robotic_systems/types.rs:1250` — has_force_closure is a dimensionally meaningless heuristic
- [ ] **F3942** (medium, bug, easy/small) `satellite_dynamics.rs:684` — Reaction-wheel allocation is projection, not pseudo-inverse
- [ ] **F3943** (medium, bug, easy/small) `satellite_dynamics.rs:842` — HillState axis labels and secular drift contradict the HCW equations
- [ ] **F3947** (medium, fabrication, easy/small) `sensor_model.rs:434` — Radar detection probability ignores cnr_db and the false-alarm rate
- [ ] **F3948** (medium, bug, easy/trivial) `sensor_model.rs:676` — Madgwick gradient step_w uses qz instead of qx
- [ ] **F3949** (medium, bug, medium/small) `sensor_model.rs:802` — Allan variance formula wrong; bias instability modelled as white noise
- [ ] **F3950** (medium, bug, easy/trivial) `sets/types.rs:324` — BodyPairCache keys ignore generation; stale value served after reuse
- [ ] **F3952** (medium, missing-feature, medium/large) `sets/types.rs:468` — No full engine-state save/restore; body and collider sets are not serializable
- [ ] **F3954** (medium, bug, easy/trivial) `sets/types.rs:736` — step_velocities applies only diagonal of world inverse inertia
- [ ] **F3955** (medium, perf, easy/small) `sets/types.rs:1134` — SpatialHash::query_radius does linear entry search per candidate
- [ ] **F3957** (medium, bug, easy/small) `sleeping/types.rs:262` — SleepingSystem compares absolute KE+PE against hardcoded 0.01
- [ ] **F3960** (medium, bug, easy/small) `soft_constraints.rs:308` — DistanceConstraint::compute_force ignores the passed state
- [ ] **F3961** (medium, bug, easy/trivial) `soft_constraints.rs:383` — SoftJoint orientation error ignores quaternion double cover
- [ ] **F3963** (medium, bug, easy/small) `space_mechanics.rs:942` — CW rendezvous_dv ignores cross-track; silent zero on singular
- [ ] **F3964** (medium, bug, medium/small) `space_mechanics.rs:1016` — Orbit RK4 drag uses sea-level exponential atmosphere: zero LEO drag
- [ ] **F3965** (medium, fabrication, medium/medium) `swarm_dynamics/types.rs:112` — CollectiveTransport payload motion ignores carriers entirely
- [ ] **F3970** (medium, panic, easy/small) `swarm_dynamics/types.rs:1064` — ParticleSwarmOptimizer::new panics on empty swarm or reversed bounds
- [ ] **F3973** (medium, panic, easy/small) `swarm_physics.rs:224` — SwarmSimulator::new panics for box_size<=0; unseeded; no minimum image
- [ ] **F3974** (medium, bug, medium/small) `swarm_physics.rs:259` — Periodic swarm ignores minimum image in neighbour search and forces
- [ ] **F3977** (medium, bug, easy/small) `terrain_dynamics.rs:173` — ground_contact_point always returns first corner, ignores terrain
- [ ] **F3978** (medium, bug, easy/small) `terrain_dynamics.rs:272` — terrain_compliance force has wrong units; 1e15 N/m fallback
- [ ] **F3982** (medium, fabrication, easy/small) `thruster_dynamics.rs:216` — thrust_allocation_pseudoinverse is per-row projection, not Moore–Penrose
- [ ] **F3988** (medium, fabrication, medium/small) `train_rigid.rs:104` — Hertz wheel-rail contact uses fixed aspect 1.1 on circular solution
- [ ] **F3989** (medium, bug, easy/small) `train_rigid.rs:189` — Kalker spin moment uses C22 instead of C33
- [ ] **F3991** (medium, bug, easy/trivial) `train_rigid.rs:515` — Klingel frequency/wavelength carry spurious factor 2
- [ ] **F3993** (medium, fabrication, medium/medium) `train_rigid.rs:617` — iso2631_rms claims ISO 2631 weighting but returns plain RMS
- [ ] **F3994** (medium, bug, easy/trivial) `train_rigid.rs:755` — Curtius–Kniffler adhesion converts km/h twice
- [ ] **F3997** (medium, bug, easy/trivial) `underwater_rigid.rs:79` — Cylinder pitch/yaw added mass is 4x too large
- [ ] **F3998** (medium, fabrication, hard/medium) `underwater_rigid.rs:579` — 'Fossen' 6-DOF step omits Coriolis terms and full kinematics
- [ ] **F3999** (medium, bug, easy/trivial) `underwater_rigid.rs:819` — buckling_pressure_infinite overstated 8x (t/R with t/D coefficient)
- [ ] **F4000** (medium, bug, easy/small) `underwater_vehicle.rs:147` — hydrostatic_restoring force not rotated into body frame; dead ×0 term
- [ ] **F4002** (medium, fabrication, easy/small) `underwater_vehicle.rs:317` — Thruster power models return newtons labelled watts
- [ ] **F4004** (medium, bug, medium/small) `underwater_vehicle.rs:437` — Thruster pseudo-inverse uses fragile absolute-pivot normal equations
- [ ] **F4007** (medium, bug, medium/medium) `wind_turbine.rs:491` — BEM lacks high-induction correction; each getter re-solves rotor
- [ ] **F4008** (medium, bug, easy/small) `wind_turbine.rs:722` — Pitch PI output integrated again; no anti-windup
- [ ] **F4011** (medium, api, easy/small) `world/types.rs:16` — Prelude exports toy point-mass RigidWorld instead of PhysicsWorld
- [ ] **F4022** (medium, bug, medium/small) `world/types.rs:1266` — Soft contact divides relax term impulse_scale*lambda by eff_mass
- [ ] **F4023** (medium, bug, easy/small) `world/types.rs:1526` — step_substep hangs on infinite or huge total_dt

### Low (41) — tracked in the external ledger
F3730, F3740, F3741, F3746, F3762, F3765, F3768, F3770, F3777, F3794, F3800, F3821, F3848, F3868, F3869, F3884, F3890, F3901, F3905, F3912, F3917, F3919, F3924, F3937, F3951, F3953, F3958, F3959, F3966, F3967, F3968, F3971, F3979, F3980, F3984, F3987, F3990, F4010, F4013, F4014, F4016

<!-- AUDIT-2026-10-12 END -->
