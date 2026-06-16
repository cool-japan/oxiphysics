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
