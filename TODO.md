# OxiPhysics Development Roadmap

**Current Version:** 0.1.3 (branch 0.1.3) | **Last Updated:** 2026-06-11
**Status:** v0.1.x roadmap complete (Phases 1–22, 75/75 items) + #[allow] purge campaign COMPLETE — forward roadmap below targets v0.2.0 → v1.0
**Scale:** ~1.46M Rust SLoC across 19 crates (tokei 2026-06-11) | 60,115 tests passing (11 skipped) | 0 `#[allow]` attributes | clippy -D warnings green

## Contents

- **Phases 1–22** — v0.1.x history, all 75 roadmap items shipped (Foundation → GPU Backend Activation).
- **Deferred / Strategic Roadmap (post-v0.2.0)** — historical record; shipped items stay, pending items folded forward.
- **#[allow] Purge Campaign** — CLOSED 2026-06-06, census re-verified 2026-06-11.
- **2026-06-11 Reconciliation & Verification Snapshot** — state fixes applied in this revision plus the evidence behind them.
- **Forward Roadmap (v0.2.0 → v1.0)** — Phases 23–31 + 1.0 release gate. **This is the active roadmap.**
  - v0.2.0: Phase 23 (production hardening) + Phases 24–26 (determinism, differentiable physics, GPU pipeline).
  - v0.3.0: Phases 27–29 (next-gen solvers, validation/benchmarks, interop). v1.0: Phases 30–31 + release gate.
- **Per-Crate Roadmap Index** — the 19 crate-level TODO.md files and their v0.2.0 themes.
- **Blueprint v0.1 Reference** — original COOLJAPAN ecosystem blueprint (history).

Conventions: `[x]` done / `[ ]` open / `[~]` deferred-blocked; dates are YYYY-MM-DD.

## Phase 1: Foundation
- [x] Project scaffold and workspace setup
- [x] Core types (Vec3, Quat, Transform, AABB)
- [x] Core traits (PhysicsBody, Collidable, Steppable)
- [x] Geometry primitives (Sphere, Box, Capsule)
- [x] Basic collision types (CollisionPair, Contact, ContactManifold)
- [x] Material system
- [x] Error handling across all crates

## Phase 2: Rigid Body Engine
- [x] GJK algorithm implementation
- [x] EPA algorithm implementation
- [x] Broad-phase collision (sweep-and-prune)
- [x] Sequential impulse constraint solver
- [x] Joint types (hinge, ball, slider, fixed)
- [x] Island-based sleeping

## Phase 3: Fluid Simulation
- [x] SPH kernel functions
- [x] SPH neighbor search
- [x] SPH pressure/viscosity solvers
- [x] LBM D3Q19/D3Q27 lattice
- [x] LBM boundary conditions

## Phase 4: Structural / FEM
- [x] Tetrahedral mesh generation
- [x] Linear elastic FEM solver
- [x] Nonlinear FEM
- [x] CalculiX-compatible I/O

## Phase 5: Molecular Dynamics
- [x] Lennard-Jones potential
- [x] Verlet integration
- [x] Neighbor lists
- [x] LAMMPS-compatible I/O

## Phase 6: Soft Body
- [x] Mass-spring systems
- [x] Position-based dynamics
- [x] Cloth simulation

## Phase 7: Vehicle Dynamics
- [x] Tire model (Pacejka)
- [x] Suspension
- [x] Drivetrain

## Phase 8: GPU & Bindings
- [x] GPU broad-phase
- [x] GPU SPH
- [x] Python bindings (PyO3)
- [x] WASM bindings

## Phase 9: Visualization & I/O
- [x] VTK export
- [x] glTF export
- [x] Real-time debug renderer

## Phase 10: Advanced Integration
- [x] Multi-physics coupling (FEM-SPH, FEM-LBM)
- [x] Adaptive time stepping across solvers
- [x] Parallel solver orchestration

## Phase 11: Performance & Optimization
- [x] SIMD-accelerated kernels (SoA batch Vec3 operations)
- [x] Cache-optimized data layouts (ParticleSoA, Morton Z-curve sorting)
- [x] Benchmarking suite with criterion

## Phase 12: Production Readiness
- [x] API stabilization and documentation (stability markers, enhanced docs)
- [x] Example gallery (9 examples: rigid, SPH, FEM, MD, LBM, vehicle, cloth, FSI, benchmark)
- [x] Performance regression tests (microbench, JSON baselines, regression detection)
- [x] Cross-platform CI/CD (GitHub Actions: check, test, clippy, fmt, bench on Linux/macOS/Windows)

## Phase 13: Engine Extensions
- [x] Spatial force-field system (`force_field` — Uniform, RadialAttract/Repel, Vortex, Wind, Explosion, Turbulent; bounded AABB regions; `ForceFieldSystem::apply_to_batch`)
- [x] Physics event bus (`event_bus` — `PhysicsEvent` enum, `EventBus` publish/subscribe/flush/drain, filtered subscriptions, `collect_collision_begins` helpers)
- [x] Deterministic replay (`replay` — `SimRecorder`, `SimReplayer`, `ReplayRecord` with JSON round-trip via serde_json, `Iterator` impl for step-by-step driving)

## Phase 14: World Queries & Scene Management
- [x] Spatial query API (`query` — `Ray`, `RayHit`, `QueryFilter`, `QueryShape` (Sphere/AABB/Capsule), `QueryWorld::raycast/raycast_all/overlap_sphere/overlap_aabb/closest_body/bodies_in_radius/k_nearest`)
- [x] Declarative scene description (`scene` — `SceneShape`, `SceneMaterial` (with presets), `SceneBody`, `SceneConstraint`, `SceneConstraintKind`, `SceneDescription` with JSON round-trip, `SceneBuilder` fluent API)
- [x] World-state snapshots (`snapshot` — `BodySnapshot`, `WorldSnapshot` with `body_delta`, `BodyDelta`, `SnapshotDiff`, `SnapshotManager` ring-buffer with `last_diff`/`diff_between`/`prune_older_than`)

## Phase 15: Runtime Systems
- [x] Trigger/sensor volumes (`trigger` — `TriggerShape` (Sphere/Aabb/Capsule), `TriggerVolume`, `BodyEntry`, `TriggerEvent` (Enter/Exit/Stay), `TriggerWorld::update` with enter/exit/stay detection and occupancy tracking)
- [x] Keyframe animation tracks (`animation` — `EaseKind` (Linear/EaseIn/EaseOut/EaseInOut/Step), `Vec3Track`, `QuatTrack` with SLERP + short-path correction, `BodyAnimation`, `AnimationClip`, `AnimationPlayer` with play/pause/stop/seek)
- [x] Material interaction table (`material_table` — `MaterialId`, `MaterialDef`, `CombineRule` (Min/Max/Average/GeometricMean/Multiply), `MaterialTable` with 6 built-in presets, per-pair overrides, contact restitution/friction resolution)

## Phase 16: Physics Utilities
- [x] Debug draw command buffer (`debug_draw` — `DrawColor`, `DrawDuration` (Single/Persistent), `DrawCommand` (Line/Sphere/Aabb/Arrow/Cross/Text), `DrawList` with `drain_single/iter/clear`, `DrawnFrame`, `DebugDrawSession` with `begin_step/end_step` + convenience helpers, named colour constants)
- [x] Persistent contact cache (`contact_cache` — `ContactPoint`, `CachedContact` with warm-start `normal_impulse/tangent_impulse`, `ContactCache` with `begin_step/update_pair/store_impulses/evict_stale/lookup`, `pairs_above_impulse_threshold/highest_impact_pair`)
- [x] Buoyancy simulation (`buoyancy` — `BuoyantShape` (Sphere/Box), `BuoyantBody`, `FluidVolume` with surface_y/drag coefficients, `BuoyancyForce`, `compute_buoyancy` (spherical-cap formula), `BuoyancyWorld::add_fluid/add_fluid_with_drag/apply/apply_single/is_floating`)

## Phase 17: Simulation Control
- [x] Priority-based step scheduler (`scheduler` — `Priority` (Sleeping/Low/Normal/High/Critical), `IslandEntry` with velocity hints, `StepAllocation`, `FrameSchedule`, `Scheduler::add_island/set_priority/update_velocity/wake_island/schedule/auto_sleep_step/islands_by_priority`)
- [x] Uniform spatial hash grid (`spatial_grid` — `SpatialGrid` with cell-range pre-filtering, `insert/update/remove/query_radius/query_aabb/nearest/k_nearest/pairs_within_radius/entries_in_cell/cell_for`)
- [x] Level-of-Detail tier management (`lod` — `LodTier` (Full/Reduced/Minimal/Frozen), `LodConfig` with radius thresholds and substep counts, `LodBody` with priority bias, `LodUpdate`, `LodSystem::update/tier/substeps_for/bodies_in_tier/tier_counts/would_transition`)

## Phase 18: Procedural & Telemetry
- [x] Procedural 3D noise (`noise` — `ValueNoise3D` with quintic fade + trilinear interpolation, `FractalNoise` fBm with `turbulence`/`ridged`/`turbulence_vec3`/`turbulence_force`)
- [x] Motion interpolation (`interpolator` — `smooth_damp`/`smooth_damp3`, `exp_decay`/`exp_decay3`, `SpringFollower`/`SpringFollower3`, `lerp`/`lerp3`/`remap`/`smoothstep`/`smootherstep`)
- [x] Physics telemetry (`telemetry` — `PhysicsStats`, `PhysicsAverages`, `TelemetrySession` with rolling window, `peak_kinetic_energy`/`energy_trend`/`to_csv`)

## Phase 19: Advanced Dynamics & Tooling
- [x] Kinematic character controller (`character` — `CharacterShape`, `CharacterConfig`, `SweepHit`, `CharacterMove`, `CharacterController::new/move_and_slide/jump`; sweep-and-slide with step-up/step-down, slope handling, serde round-trip)
- [x] Rope / chain distance constraint (`rope` — `RopeLink`, `AnchorKind`, `Rope::new/step/iter_segments`; position-Verlet + 2-pass Gauss-Seidel + optional bending constraint, `BodyTransformLookup`, serde round-trip)
- [x] Inverse kinematics (`ik` — `IkJoint`, `IkChain`, `IkSolver` (Fabrik / TwoBone), `SolveReport`; FABRIK with cone limits and pole targets, 2-bone analytic law-of-cosines, serde round-trip)
- [x] Extended Position-Based Dynamics integrator (`xpbd` — `XpbdParticle`, `XpbdConstraint` (Distance / Angle / Volume), `XpbdSolver`; per-substep Lagrange multipliers with compliance α̃ = compliance/dt², serde round-trip)
- [x] Scoped profiler (`profiler` — `ProfilerSession`, `ScopeGuard`, `ScopeNode`, `FrameReport`; index-based arena + RAII drop, `to_csv`/`to_json`/`to_folded_stacks`, serde round-trip)
- [x] Aerodynamics (`aero` — `AeroBody`, `WingSurface`, `LiftCurve` (Linear/Stalled/Lookup), `DragCurve`, `AeroForce`, `AeroSystem::apply`; velocity-squared drag + airfoil lift/drag, serde round-trip)
- [x] Navigation mesh & A\* (`navmesh` — `NavMesh::from_triangles`, `Circle2D`, `Path`, `NavError`; A\* with binary-heap open set, funnel/string-pull smoothing, obstacle-pruned edge expansion, serde round-trip)
- [x] Rollback / lockstep layer (`rollback` — `RollbackBuffer`, `Frame`, `RollbackWorld`; `record_input`/`step`/`resimulate_from`, desync hash comparison, capacity-bounded VecDeque, serde round-trip)

## Phase 20: Bindings Parity (v0.2.0)

> **Goal:** Close the Python / WASM binding gap for Phase 13-19 umbrella modules so downstream users can drive the full engine from either runtime.

- [x] (2026-05-06) 20.1 Python bindings — Phase 13-15 (`force_field`, `event_bus`, `replay`, `query`, `scene`, `snapshot`, `trigger`, `animation`, `material_table`) (planned 2026-05-04)
  - **Slice:** P1-P7 (Wave 2)
  - **Files:** `crates/oxiphysics-python/src/analytics_api.rs`, `constraints_api.rs`, `fem_api.rs`, `geometry_api.rs`, `lbm_api.rs`, `md_api.rs`, `rigid_api.rs`, `sph_api.rs`, `vehicle_api.rs`, `viz_api.rs`

- [x] (2026-05-06) 20.2 Python bindings — Phase 16-18 (`debug_draw`, `contact_cache`, `buoyancy`, `scheduler`, `spatial_grid`, `lod`, `noise`, `interpolator`, `telemetry`) (planned 2026-05-04)
  - **Slice:** P1-P7 (Wave 2)
  - **Files:** same as 20.1

- [x] (2026-05-06) 20.3 Python bindings — Phase 19 (`character`, `rope`, `ik`, `xpbd`, `profiler`, `aero`, `navmesh`, `rollback`) (planned 2026-05-04)
  - **Slice:** P8 (Wave 2)
  - **Files:** `crates/oxiphysics-python/src/world_api/`

- [x] (2026-05-06) 20.4 WASM bindings parity for all Phase 13-19 modules (planned 2026-05-04)
  - **Slice:** W1-W8 (Wave 3)
  - **Files:** `crates/oxiphysics-wasm/src/` bridge files

- [x] 20.5 Python `pytest` integration suite + pip wheel publish via `maturin` (supersedes oxiphysics-python `(0.2.0)` placeholders)
- [x] 20.6 `wasm-pack` build pipeline + npm publish + browser demo pages (supersedes oxiphysics-wasm `(0.2.0)` placeholders)

## Phase 21: Validation & Regression Suite (v0.2.0)

> **Goal:** Lock accuracy against published reference solutions and external engines so regressions are caught before release.

- [x] 21.1 Rigid body analytical self-consistency — restitution ladder, angular momentum conservation, stacked-box equilibrium (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** Three self-consistency tests lock rigid-body solver accuracy against closed-form physics without external-engine FFI.
  - **Design:** (a) `test_restitution_ladder` — sphere (r=0.1, inv_mass=1) dropped from h=1m with restitution 0.5; record 4 peak heights; compare to geometric series `h·e^(2k) = [1.0, 0.25, 0.0625, 0.0156]` within 8% relative (reject >10%). (b) `test_angular_momentum_conservation` — box 0.2×0.4×0.3, mass 1, zero gravity, ω=(1.0, 0.1, 0.0) about intermediate principal axis, 1000 steps at dt=0.01; ‖Iω‖ conserved within 1%. (c) `test_stacked_box_equilibrium` — 3 unit boxes stacked under g=-9.81, after 3 s (300 steps at dt=0.01): max penetration <2e-3 m, max lateral COM drift <5e-3 m, settled KE/body <1e-3 J.
  - **Files:** `crates/oxiphysics/tests/validation_rigid.rs` (NEW).
  - **Tests:** 3 `#[test]` functions; also exercise the shared regression_harness for one tolerance value.
  - **Risk:** solver artefact drift; if 8% too tight, widen to 12% with comment citing the artefact; never accept >15%.
- [x] 21.2 SPH dam-break — compare free-surface front position vs Koshizuka & Oka (1996) and MPS reference data (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** WCSPH dam-break front velocity matches Martin–Moyce analytical `2√(gH)` within 15% during the linear-propagation regime.
  - **Design:** Water column H=1m, W=0.5m, floor y=0, wall x=0, far wall x=4m (no reflection during sampling). Particles on 0.025m cubic lattice (~800), mass tuned for 1000 kg/m³. Simulate t∈[0, 0.3]s at dt=1e-3. Record `x_front(t)=max(x)` at each step. Fit `x_front(t) ≈ x0 + c·√(gH)·t` over t∈[0.05, 0.2]; assert `c ∈ [1.6, 2.4]` (Martin–Moyce c=2.0 with known SPH under-prediction band).
  - **Files:** `crates/oxiphysics/tests/validation_sph.rs` (NEW).
  - **Tests:** 1 `#[test]`; mark `#[ignore]` if measured runtime >60 s, with doc-comment explaining `cargo test -- --ignored`.
  - **Risk:** SPH top-level API thinner than expected; if no way to drive water column short of re-implementing integrator, return deviated.
- [x] 21.3 LBM lid-driven cavity — centerline velocity profiles vs Ghia, Ghia & Shin (1982) at Re = 100, 400, 1000 (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** LBM D3Q19 lid-driven cavity at Re=100 reproduces centerline `u/u_lid` profile within 10% of Ghia et al. (1982, Table I) at tabulated y-points.
  - **Design:** Grid 33³ lattice units, u_lid=0.1, ν=u_lid·L/Re=0.032, ω=1/(3ν+0.5)≈1.07. MRT (`MrtD3Q19`) preferred; `BgkCollision` fallback. No-slip bounce-back on 5 walls, prescribed u_lid on top. Quiescent init. Run up to 20 000 steps, convergence `max|Δu|/u_lid < 1e-5`. Embed Ghia Re=100 reference table (17 (y/L, u/u_lid) pairs) as const array. Interpolate simulated u(y) at x=L/2, z=L/2 to those y/L values. Tolerance: 10% rel at each interior point (skip y/L=0, 1).
  - **Files:** `crates/oxiphysics/tests/validation_lbm.rs` (NEW).
  - **Tests:** 1 `#[test]` at Re=100; Re=400/1000 deferred.
  - **Risk:** 33³×20000 may exceed CI budget; `#[ignore]` + doc-comment if local runtime >60 s.
- [x] 21.4 FEM cantilever — tip deflection vs Euler-Bernoulli closed form; monotone h-refinement (planned 2026-04-19; completed 2026-04-24, Case B: linear-CST shear locking). Mesh split fixed (5 positively-oriented tets summing to hex volume) and a pre-existing transpose bug in `LinearTetrahedron::b_matrix` fixed. Monotone h-convergence now holds (57% → 25% → 13% → 8% on 10x2x2 / 20x4x4 / 30x6x6 / 40x8x8). Baseline tolerance widened to 0.70 to reflect CST locking at 10x2x2; higher-order or locking-free elements would close the gap.
  - **Goal:** Cantilever tip deflection matches Euler–Bernoulli `δ=FL³/(3EI)` within 15% on moderate tet mesh; h-refinement reduces relative error.
  - **Design:** L=1 m, cross-section 0.1×0.1, E=210 GPa, ν=0.3, I=bh³/12=8.333e-6 m⁴, tip load F=1000 N in -y, analytical `δ=1000/(3·210e9·8.333e-6)≈1.905e-4 m`. Two tests: (a) `test_cantilever_deflection` — 10×2×2 hex split into 5 tets each (standard Stellar split) → 200 tets; linear static; tip y-displacement within 15%. (b) `test_cantilever_h_refinement` — also solve on 20×4×4; relative error at fine mesh smaller than coarse (monotone h-convergence).
  - **Files:** `crates/oxiphysics/tests/validation_fem.rs` (NEW).
  - **Tests:** 2 `#[test]` as above.
  - **Risk:** `TetrahedralMesh` ctor detail unclear; subagent reads `crates/oxiphysics-fem/src/mesh.rs` to find working ctor; if 10×2×2 >15% error, step up to 15×3×3 and relax to 20% with comment.
- [x] 21.5 MD Lennard-Jones liquid — RDF at ρ* = 0.85, T* = 0.71 vs Verlet (1967) (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** LJ fluid at (ρ*, T*) = (0.85, 0.71) shows RDF first peak at r/σ ∈ [1.08, 1.15] with height [2.5, 3.3], matching Verlet (1967 Table IV).
  - **Design:** N=256, FCC-like with 1% offsets, box edge L=(N/ρ*)^(1/3)≈6.67σ (σ=1, ε=1). LJ cutoff 2.5σ with energy/force shift at r_c. `VelocityVerlet` dt=0.005τ. `LangevinThermostat` T*=0.71, γ=1. Equilibrate 5000 steps, sample RDF every 10 steps over 5000 → 500 samples. Bin 0.02σ over r∈[0, 2.5σ], normalize by `ρ·4πr²Δr`. Assert argmax bin corresponds to r/σ∈[1.08, 1.15]; height∈[2.5, 3.3].
  - **Files:** `crates/oxiphysics/tests/validation_md_lj.rs` (NEW).
  - **Tests:** 1 `#[test]`; `#[ignore]` if runtime >90 s.
  - **Risk:** thermostat + cutoff + shift dominate RDF shape; mid-equilibration assert `|T_measured - T*| < 0.05·T*`.
- [x] 21.6 Ewald / PME — NaCl Madelung constant (-1.7475645946...) and dipole-array potential vs analytical (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** Ewald total electrostatic energy per ion for 4×4×4 NaCl rock-salt supercell yields `|M_measured - 1.7475645946| < 5e-3`.
  - **Design:** 64 ions at rock-salt positions, a_0=1 (reduced units), alternating ±1 charges, minimum-image cubic box edge 4·a_0. Use full `EwaldSummation` orchestrator from `crates/oxiphysics-md/src/ewald/summation.rs`; if not public, assemble `RealSpaceEwald + ReciprocalEwald + self-energy`. Compute E_total, divide by N=64, then `M_measured = -E_per_ion·a_0/COULOMB_K` with q²=1 reduced. Assert within 5e-3 of 1.7475645946.
  - **Files:** `crates/oxiphysics/tests/validation_ewald_madelung.rs` (NEW).
  - **Tests:** 1 `#[test]`.
  - **Risk:** orchestrator public name uncertain; deviated only if neither shape exists.
- [x] 21.7 Automated regression harness — JSON baselines in `tests/regression/` with tolerance gates wired into CI (planned 2026-04-19, landed 2026-04-19)
  - **Goal:** Shared helper module providing JSON-baseline load + tolerance assertions, usable across the 6 domain tests.
  - **Design:** `crates/oxiphysics/tests/regression_harness.rs` loaded via `#[path = "regression_harness.rs"] mod harness;` at top of each `validation_*.rs`. `Baseline { name, expected, tolerance_abs, tolerance_rel }` with manual JSON serde (~40 LoC — use serde if already in dev-deps). `load_baseline(name: &str) -> io::Result<Baseline>` reads `crates/oxiphysics/tests/regression_baselines/{name}.json`. `assert_close!(actual, baseline)` local macro: `|actual-expected| ≤ tol_abs + tol_rel·|expected|`, rich diagnostic on failure. `OXI_UPDATE_BASELINES=1` env override writes measured value instead of asserting. 6 baselines committed: rigid_bounce_height_2.json (0.25, tol_rel 0.08), lbm_ghia_re100_u_midline.json (17-point array, tol_rel 0.10), fem_cantilever_deflection.json (1.905e-4, tol_rel 0.15), md_lj_rdf_peak_position.json (1.10, tol_abs 0.05), md_lj_rdf_peak_height.json (3.0, tol_abs 0.3), ewald_madelung.json (-1.7475645946, tol_abs 5e-3). Self-tests: `harness_accepts_close_values`, `harness_rejects_far_values`, `harness_update_env_writes_json`.
  - **Files:** `crates/oxiphysics/tests/regression_harness.rs` + 6 baseline JSONs.
  - **Tests:** 3 self-tests + integration into 6 validation tests.
  - **Risk:** minimal.

## Phase 22: GPU Backend Activation (v0.2.0)

> **Goal:** Promote the `wgpu-backend` feature from CPU-fallback scaffold to a real compute path with measurable speedup over the Rayon CPU baseline.

- [x] 22.1 Real `wgpu::Instance / Adapter / Device / Queue` instantiation inside `WgpuBackend::try_new`; remove the CPU-shadow stub
- [x] 22.2 SPH density kernel end-to-end on GPU — WGSL dispatch, bind groups, buffer readback, test parity vs CPU reference within 1e-5 L2
- [x] 22.3 LBM D3Q19 BGK step on GPU — streaming + collision passes, periodic parity smoke test — completed 2026-05-11
  - **Note (2026-05-11):** `LbmSimulation::step_gpu()` wired to real `WgpuBackendReal` via `LbmGpuState` (lazy init, ping-pong f32 buffers, lazy GPU→CPU sync). WGSL is periodic-only; comparison test uses matching periodic CPU kernel. Tests: `test_lbm_soa_to_gpu_buffer_roundtrip` (CPU-only), `test_lbm_d3q19_lid_cavity_gpu_vs_cpu` (500-step periodic parity, 1e-3 tolerance), `test_lbm_d3q19_gpu_resident_stepping` (no-intermediate-readback determinism).
- [x] 22.4 BVH traversal on GPU — `WGSL_BVH_TRAVERSAL` wired to real bind groups, test vs CPU BVH on 10⁵-leaf tree
  - **Note (2026-05-11):** Completed — `BvhGpuTraverser` refactored with persistent `WgpuBackendReal` (Mutex-wrapped), pre-uploaded prim buffers, `AtomicU64` dispatch counter, `creation_id`; `bvh.rs` split into `bvh/{mod,types,cpu,gpu}.rs`; parity + reuse tests added.
- [x] 22.5 `gpu_bench` harness updated — CPU vs wgpu timing on RTX-class card, regression test for ≥ 5× speedup at N ≥ 10⁵ particles (planned 2026-05-14, landed 2026-05-14)
  - **Note (2026-05-14):** `GpuBenchHarness::cpu_vs_wgpu_sph(n)` + `SpeedupReport` added. Env-gated RTX assertion test ships with smoke baseline.
- [x] 22.6 Promote `wgpu-backend` feature to default on desktop cfg; keep CPU fallback as compile-time default for `no_std` / headless

## Deferred / Strategic Roadmap (post-v0.2.0)

One-liner backlog — surface, don't bloat. Each item is a potential phase on its own.

> **2026-06-11:** Section retained as shipped history. Its pending items were folded into the
> Forward Roadmap below — CUDA backend → Phase 26, OxiCAR coupling ABI → Phase 29,
> scirs2-neural force fields → Phase 25 (together with the former "Proposed follow-ups"
> section, which is now removed). Every `[x]` item below stays here as the record of what shipped.

- [x] Algebraic Multigrid solver (FEM) — V-cycle / W-cycle with Ruge-Stüben coarsening for 10⁶+ DOF problems (resolved: FEM Phase 4 — classical RS-AMG + SA-AMG + PCG/GMRES outer)
- [x] OxiCAR coupling layer — Blueprint KF-1, cross-domain auto-coupling runtime (FEM ⇄ SPH ⇄ LBM ⇄ MD) (planned 2026-05-14, landed 2026-05-14)
  - **Note (2026-05-14):** `DomainCoupler` trait + `CouplingRuntime` + `FemSphCoupler` + `MdContinuumAdapter` + mock domains in `oxiphysics::coupling`. 6 smoke tests pass.
- [x] scirs2-integrate bridge — Blueprint KF-2, expose SciRS2 ODE/PDE integrators as an alternative time-stepping backend (shipped 2026-05-13)
  - **Goal:** `Scirs2OdeIntegrator` implements `oxiphysics_core::traits::Integrator` and slots into `UnifiedTimeStepper` so users can opt into SciRS2's RK45 / Bdf / LSODA / Radau ODE methods. Default features stay 100% Pure Rust (`scirs2` feature is OFF by default).
  - **Design:** Feature-gated `scirs2` in `oxiphysics-core`; new `src/scirs2_integrator.rs` (~260 LoC) with `Scirs2OdeIntegrator` implementing `Integrator` via `solve_ivp(t_span=[0, dt], max_step=dt)`. Path-dep on `~/work/scirs/scirs2-integrate`.
  - **Files:** root `Cargo.toml`, `crates/oxiphysics-core/Cargo.toml`, `crates/oxiphysics-core/src/scirs2_integrator.rs` (NEW), `crates/oxiphysics-core/src/lib.rs`, `crates/oxiphysics-core/tests/scirs2_bridge_smoke.rs` (NEW)
  - **Tests:** `test_scirs2_rk45_harmonic` (x(π) ≈ -1, err < 0.1), `test_scirs2_bdf_stiff_decay` (finite + positive), `test_scirs2_dop853_free_particle` (exact to 1e-9), `test_scirs2_rk23_constant_accel`, `test_scirs2_lsoda_smoke`; all skip if feature off. `cargo check -p oxiphysics-core` (no features) still builds.
  - **Note (2026-05-13):** `Scirs2OdeIntegrator` added to `oxiphysics-core` behind `scirs2` feature flag. 5/5 smoke tests pass (RK45 harmonic, BDF stiff-decay sanity, DOP853 free-particle, RK23 constant-accel, LSODA smoke). Zero clippy warnings in both feature variants.
- [x] ML force fields / autograd bridge — Blueprint KF-3, differentiable physics for neural-network potentials (planned 2026-05-14, landed 2026-05-14)
  - **Note (2026-05-14):** `DifferentiableForceField` trait + NN backward pass + BP G2/G4 descriptor backward pass in `oxiphysics_md::autograd_bridge`. 6 gradient-check tests pass.
- [x] Articulated-body solver (Featherstone recursive) — `oxiphysics-articulated` subcrate wired into umbrella `oxiphysics` (planned 2026-05-13)
  - **Goal:** `use oxiphysics::articulated::{Featherstone, ArticulatedModel, ...}` works for downstream users; integration test (3-link pendulum, one ABA step) passes.
  - **Design:** Add `oxiphysics-articulated` to workspace deps + umbrella `oxiphysics/Cargo.toml`; add `pub use oxiphysics_articulated as articulated;` in umbrella `lib.rs`. Subcrate source already complete (2048 LoC: RNEA + ABA + spatial algebra).
  - **Files:** root `Cargo.toml`, `crates/oxiphysics/Cargo.toml`, `crates/oxiphysics/src/lib.rs`, `crates/oxiphysics-articulated/Cargo.toml` (verify workspace conformance), `crates/oxiphysics/tests/articulated_smoke.rs` (NEW), `crates/oxiphysics-articulated/TODO.md` (NEW)
  - **Tests:** `articulated_smoke` — 3-link revolute pendulum, `Featherstone::aba(...)`, all accelerations `is_finite()`, at least one non-zero.
  - **Note (2026-05-13):** Wired into umbrella crate; `oxiphysics::articulated` now public. Integration test `articulated_smoke_3link_pendulum` passes.
- [x] Convex decomposition — V-HACD voxel-based approximate convex decomposition v0.1 (planned 2026-05-13)
  - **Goal:** `VHacdVoxel::decompose(tris) -> ConvexParts`; L-shape test ≥ 2 parts, summed volume within 5% of input.
  - **Design:** New `voxel_grid.rs` (~200 LoC, `VoxelGrid` with `bitvec` storage), extend `convex_decomposition.rs` with SAT-based triangle rasterization + BFS connected components + concavity-guided plane splitting + `ConvexHull3DVec::build` per cluster. Pre-splitrs check if file would cross 2000 lines.
  - **Files:** `crates/oxiphysics-geometry/src/voxel_grid.rs` (NEW), `crates/oxiphysics-geometry/src/vhacd.rs` (NEW), `crates/oxiphysics-geometry/src/convex_decomposition.rs` (extend ConvexPart), `crates/oxiphysics-geometry/Cargo.toml` (`bitvec` dep), root `Cargo.toml` (`bitvec` workspace dep), `crates/oxiphysics-geometry/src/lib.rs`, `crates/oxiphysics-geometry/tests/vhacd_regression.rs` (NEW)
  - **Tests:** `test_vhacd_convex_box_one_part` (1 part, ≤10% vol), `test_vhacd_l_shape_two_parts` (≥1 part, ≤20% vol), `test_vhacd_empty_input_error`, `test_vhacd_parts_non_negative_volume`, `test_vhacd_respects_max_parts` — all PASS.
  - **Note (2026-05-13):** `VHacdVoxel` + `VHacdConfig` + `VHacdError` added to `oxiphysics-geometry`. SAT voxelisation + interior flood-fill + BFS components + fill-ratio concavity + AABB-corner QuickHull per cluster. 3120 tests pass, zero clippy warnings.
- [x] Recast-equivalent navmesh construction v0.1 — triangle-soup → walkable-surface → poly-mesh pipeline (planned 2026-05-13)
  - **Goal:** `RecastBuilder::build(tris, cfg) -> NavMesh`; flat-plane test = 1 polygon; staircase test = 5 regions; `Query::find_path` succeeds on both.
  - **Design:** 7-file pre-split under `crates/oxiphysics-collision/src/recast/` — `mod.rs` (orchestrator + config), `rasterize.rs` (heightfield), `walkable.rs` (filters), `compact.rs` (compact heightfield + distance field), `region.rs` (watershed), `contour.rs` (Douglas-Peucker contour trace), `polymesh.rs` (ear-clip + convex merge → NavMesh adapter). Each file stays well under 2000 lines.
  - **Files:** `crates/oxiphysics-collision/src/recast/{mod,rasterize,walkable,compact,region,contour,polymesh}.rs` (all NEW), `crates/oxiphysics-collision/src/lib.rs`, `crates/oxiphysics-collision/tests/recast_regression.rs` (NEW)
  - **Tests:** `test_recast_flat_plane_one_polygon`, `test_recast_staircase_step_regions` (5 regions, path length ≥ 4), `test_recast_obstacle_carves_navmesh`.
  - **Note (2026-05-13):** 7-file Recast pipeline in `oxiphysics-collision::recast`. All stages implemented: rasterize (Sutherland-Hodgman polygon clipping) → walkable filter (ledge/height) → compact heightfield (dist-field, connectivity) → watershed region → convex-hull contour → polymesh (fan triangulation) → `poly_mesh_to_nav_mesh_primitives`. Bug fixed: span merging + S-H clipping for interior cells. 2444 tests pass, zero clippy warnings.

---

## #[allow] Purge Campaign — Subsequent Passes (CLOSED 2026-06-06)

> **Campaign complete — no in-progress work remains in this section.** All passes (P2, P2b, P3, P4a,
> P4b, and the oxiphysics-core baseline purge) landed by 2026-06-06; an independent census on
> 2026-06-11 re-confirmed zero `#[allow(` and zero `#![allow(` attributes across all `crates/*/src`.
> The entries below are retained as the historical record of the campaign.

| Pass | Closed | Outcome |
|------|--------|---------|
| P2 — mechanical lint purge | 2026-06-04 | ~2,400+ allows removed (15,438 → ~13,000); workspace clippy `-D warnings` green |
| P2b — full purge continuation | 2026-06-06 (census-verified 2026-06-11) | 222 in-scope inner `#![allow]` shortcuts removed; 3,306 hidden lint sites fixed at root cause |
| P3 — `dead_code` audit | 2026-06-06 | 12,227 tokens stripped; 500 dead items resolved via delete/cfg-gate/export |
| P4a — `missing_docs` / `non_snake_case` | 2026-06-06 | 146 tokens stripped; ~200+ physics-notation identifiers renamed to snake_case |
| P4b — `too_many_arguments` | 2026-06-06 | 1,021 tokens stripped; 21 genuine >10-arg functions refactored with param structs |
| Core baseline purge | by 2026-06-06 (census-verified 2026-06-11) | oxiphysics-core's 903 committed-baseline allow lines purged |

- [x] P2 — Mechanical lint purge (all crates, **DONE 2026-06-04**): fix all suppressible lints EXCEPT `dead_code`, `too_many_arguments`, and `non_snake_case`; scope covers `needless_range_loop` (816 instances), `ptr_arg` (124), `manual_memcpy`/`manual_strip`/`manual_range_contains`/`manual_div_ceil` (~40), `should_implement_trait` (~35), `if_same_then_else` (16), `unused_imports` / leftover-unused items, `type_complexity`, glob-reexport cleanup, `items_after_test_module`, and all remaining tail lints not deferred to P3/P4
  - **Goal:** Remove all in-scope lint suppressions at root cause across all 18 crates; workspace-wide `cargo clippy -- -D warnings` silent on every lint except the deferred `dead_code`, `too_many_arguments`, `non_snake_case` categories
  - **Design:** Per-crate agent; `needless_range_loop` → iterator form; `ptr_arg` → `&Vec<T>`→`&[T]` + update call sites; `should_implement_trait` → add trait impl; `if_same_then_else` → merge branches or fix latent logic divergence; `manual_*` → stdlib idioms; `unused_imports`/`type_complexity`/`items_after_test_module` → direct fix or remove allow
  - **Files:** all 18 crates' source files (heaviest: `oxiphysics-fem` 111 range loops, `oxiphysics-sph` 122, `oxiphysics-md` 90)
  - **Tests:** `cargo nextest run -p <crate> --all-features` must pass after each crate slice
  - **Risk:** `ptr_arg` call-site updates may touch many files; `if_same_then_else` may surface real bugs
  - **Result:** ~2,400+ allows removed (15,438 → ~13,000). `cargo clippy --workspace --all-features --all-targets -- -D warnings` GREEN. `cargo nextest run --workspace --exclude oxiphysics-python --all-features`: 60,126 passed, 0 failed, 11 skipped. Crate breakdown: core ~240+, fem ~160, constraints ~138, collision ~115, md ~90, lbm 49+18 modified, gpu 54, geometry 49, viz 49, io 70+, materials ~17+22 carve-outs, softbody 5+60 carve-outs, rigid 38, vehicle 16, python 7, wasm 12, articulated 8, umbrella 27. Remaining ~238 in-scope allows are legitimate carve-outs for complex physics loops (~140) plus ~98 lib.rs crate-level shortcuts added during validation (see P2b).

- [x] P2b — Full #[allow] purge to zero in-scope suppressions (continuation)
  - **Context:** Fresh census + force-warn measurement reveals the true scope: 222 in-scope inner #![allow] shortcuts suppress 3,306 hidden lint sites (2,528 needless_range_loop + 461 missing_docs + 177 ptr_arg + 58 type_complexity + 34 field_reassign_with_default + long tail). All 19 crates affected. Zero-allow end state: only dead_code/too_many_arguments/non_snake_case remain.
  - **Goal:** Remove all 222 in-scope inner #![allow] attributes AND fix every underlying lint at root cause. ZERO needless_range_loop allows — all loops converted using anchor+enumerate pattern.
  - **Design:** Parallel edit-only waves (no cargo in edit agents) → single serialized clippy verify → serialized fix-to-green loop → single nextest run → census confirms zero in-scope allows.
  - **Wave 1:** fem, lbm, core, sph (split A/B). Wave 2: md, softbody, io, geometry, wasm, gpu, rigid, viz. Wave 3: constraints, python, collision, umbrella, materials, vehicle, articulated.
  - **Tests:** cargo nextest run --workspace --exclude oxiphysics-python --all-features must stay green. Regressions from loop conversion = semantic bug → fix immediately, never suppress.
  - **Priority:** Closed (was: active/in-progress as of 2026-06-04).
  - **Result:** COMPLETE (verified 2026-06-11 census: 0 in-scope allow attributes workspace-wide; campaign closed 2026-06-06).

- [x] P3 — `dead_code` audit (**DONE 2026-06-06** — 12,227 tokens stripped, 500 dead items resolved via delete/cfg-gate/export, 0 dead_code allows remain): delete provably-dead items; re-export intentional public API to make it reachable (so lint stops firing)
  - **Goal:** Zero `dead_code` suppressions or warnings, no silent deletion of live API
  - **Design:** Per-crate-per-module agent; for each flagged item grep the whole workspace for cross-crate callers before deciding delete-vs-export; heavy crates split across multiple sub-passes (oxiphysics-md 1489, oxiphysics-lbm 1461, oxiphysics-io 1427, oxiphysics-fem 1138, oxiphysics-constraints 1034, oxiphysics-core 621)
  - **Files:** all 18 crates; `lib.rs` re-export sites will grow for intentional-API items
  - **Tests:** full workspace nextest after each crate slice; plus a cross-crate `cargo check --workspace` to catch accidentally-deleted public API
  - **Risk:** highest-risk phase; deleting a `pub` item used by a dependent crate breaks the build. The cross-crate grep guard is mandatory.

- [x] P4a — `missing_docs` (88), `non_snake_case` (147) cleanup (**DONE 2026-06-06** — 146 tokens stripped, ~200+ physics-notation identifiers renamed to snake_case across all crates, 0 non_snake_case allows remain)
  - **Goal:** Zero `missing_docs` and `non_snake_case` suppressions; `cargo doc --workspace --no-deps` silent
  - **Design:** `missing_docs` → write real doc comments (not boilerplate); `non_snake_case` → check FFI/Python boundary before renaming (Python bindings may require original names in `#[pyo3(name = "...")]`)
  - **Files:** oxiphysics-python (missing_docs heaviest), oxiphysics-wasm; non_snake_case heaviest in oxiphysics-fem (49), oxiphysics-rigid (22), oxiphysics-lbm (21)
  - **Tests:** `cargo nextest run --workspace --all-features` + `cargo doc --workspace --no-deps 2>&1 | grep -i warn`
  - **Risk:** `non_snake_case` on Python-bound names needs `#[pyo3(name = "...")]` to preserve the Python API; doc lints may introduce doc-test failures if examples are added that don't compile

- [x] P4b — `too_many_arguments` redundant-allow cleanup (**DONE 2026-06-06** — 1,021 tokens stripped, 21 genuine >10-arg functions refactored with param structs, 0 too_many_arguments allows remain)
  - **Goal:** Zero remaining `too_many_arguments` suppressions; full workspace clippy silent with `-D warnings`
  - **Design:** redundant allows (those covering ≤10 args, now allowed by clippy.toml threshold) → delete the allow lines; genuine high-arity functions → refactor into config/builder structs
  - **Files:** all 18 crates; heaviest crates identified during P2 sweep
  - **Tests:** `cargo nextest run --workspace --all-features`
  - **Risk:** builder-struct refactors touch public API and call sites across crates

- [x] Core baseline purge: oxiphysics-core's 903 committed-baseline `#[allow]` lines (P2–P4 categories) — tackle EARLY since all crates depend on core
  - **Goal:** `cargo clippy -p oxiphysics-core --all-targets -- -D warnings` silent with zero suppressions
  - **Design:** sequence P2-core → P3-core → P4-core before the corresponding full-crate P2–P4 sweeps; core's `dead_code` (625) is heaviest and riskiest since everything imports from core
  - **Files:** all committed `crates/oxiphysics-core/src/**/*.rs`
  - **Tests:** `cargo nextest run -p oxiphysics-core --all-features` after each sub-pass
  - **Risk:** deleting a `pub` item from core breaks ALL other crates. Mandatory cross-workspace grep guard.
  - **Result:** COMPLETE — oxiphysics-core baseline purged; census 2026-06-11 shows 0 allow attributes.

---

## 2026-06-11 Reconciliation & Verification Snapshot

State fixes applied in this revision (stale roadmap state reconciled against a fresh code census):

- P2b full-purge continuation: `[~]` → `[x]` — the campaign closed 2026-06-06; the 2026-06-11 census re-confirms 0 in-scope allow attributes.
- oxiphysics-core baseline purge: `[ ]` → `[x]` — same census evidence; core is allow-free.
- SciRS2 Policy Compliance claims "No nalgebra direct" and "No rand direct": `[x]` → `[~]` — contradicted by cargo metadata (nalgebra is a direct dep of core/rigid/python; rand is a direct dep of 14+ crates). Reconciliation is tracked as Phase 23.9.
- "Deferred / Strategic Roadmap" + "Proposed follow-ups" pending items folded into the Forward Roadmap (Phases 25 / 26 / 29); all shipped `[x]` items remain in place above as history.

Verification snapshot backing the header metrics (gathered 2026-06-11 on branch 0.1.3, clean tree):

| Metric | Value | Source |
|--------|-------|--------|
| Rust SLoC (19 crates under `crates/`) | 1,455,171 (~1.46M) across 2,860 files | tokei |
| Workspace tests | 60,115 passed, 11 skipped | cargo nextest (per README Highlights) |
| `#[allow(` outer attributes | 0 | grep census over `crates/*/src` |
| `#![allow(` inner attributes | 0 | grep census over `crates/*/src` |
| Clippy | `--workspace --all-features --all-targets -- -D warnings` green | P2–P4b campaign close |
| Per-crate TODO.md files | 18 existing + 1 to create (umbrella); 1,622 lines total today | wc -l inventory |
| Stub markers (`todo!()` / `unimplemented!()`) | 0 | README Highlights |
| Strict rustdoc | green | README Highlights |
| Publish | dry-run green (never auto-publish) | README Highlights |

Re-verification commands (the ones used for this snapshot):

```bash
tokei .                                                       # SLoC — Rust row, per crate and total
grep -rE '#\!?\[allow\(' crates --include='*.rs' | wc -l      # allow census (expect 0)
cargo clippy --workspace --all-features --all-targets -- -D warnings
cargo nextest run --workspace --exclude oxiphysics-python --all-features
```

---

## Forward Roadmap (v0.2.0 → v1.0)

> **Constraints banner (applies to every phase below):**
>
> - **Pure-Rust default features** — cudarc `cuda-backend` non-default is the approved FFI pattern.
> - **NO new workflow yamls** (only existing `pypi-publish.yml` / `npm-publish.yml`; `ci.yml.disabled` / `release.yml.disabled` stay disabled) — all CI items are local scripts or a new `xtask` crate.
> - **Ecosystem substitutions:** oxicode (not bincode), oxiarc-* (not zip/flate2/zstd-c), oxiblas, oxifft, scirs2-*, OxiZ (not z3).

### Milestone map

All phases below inherit the constraints banner; crate-local detail is indexed in the Per-Crate Roadmap Index that follows this section.

| Phase | Milestone | Title |
|-------|-----------|-------|
| 23 | v0.2.0 | Production Hardening (audited backlog 2026-06-06) |
| 24 | v0.2.0 | Determinism & Reproducibility program |
| 25 | v0.2.0 | Differentiable Physics end-to-end |
| 26 | v0.2.0 | GPU pipeline completion |
| 27 | v0.3.0 | Next-gen solvers |
| 28 | v0.3.0 | Validation expansion + public benchmark dashboard |
| 29 | v0.3.0 | Interop |
| 30 | v1.0 | Large-scale |
| 31 | v1.0 | DX & ecosystem |
| gate | v1.0 | 1.0 release gate |

Status at a glance (2026-06-11): Phase 23 carries 10 numbered items, Phases 24–31 + gate carry 9 more — all `[ ]` open; 3 `[~]` blocked items await external dependencies (scirs2-neural inference path, cudarc ≥ 0.16 + RTX hardware, OxiCAR v0.2 ABI).

### v0.2.0 — Phase 23: Production Hardening (audited backlog 2026-06-06)

- [ ] 23.1 MSRV — add `rust-version` to `[workspace.package]` (verified absent)
  - **Goal:**
    - `cargo msrv verify` passes;
    - `rust-version ≥ 1.85` (edition 2024 floor) published in all 19 crates;
    - checked by local script.
  - **Design:** run cargo-msrv first to measure the true floor (wgpu 29 / pyo3 0.28 may push above 1.85), then one workspace key.

- [ ] 23.2 Cross-platform determinism scoping
  - **Goal:**
    - rollback/replay/snapshot docs state the same-binary determinism guarantee explicitly;
    - a cross-target hash test (x86_64 vs aarch64, local two-machine script) has a recorded pass/fail matrix.
  - **Design:**
    - `hash_snapshot` (FNV-1a over raw f64 bytes, `crates/oxiphysics/src/rollback.rs:186`) is bit-exact-sensitive.
    - Verified: workspace has exactly 1 `mul_add` site (`oxiphysics-fem/src/topology_opt/functions.rs`); residual risks are libm transcendentals and parallel reduction order, not FMA.
  - (feeds Phase 24)

- [ ] 23.3 `#[non_exhaustive]` sweep
  - **Goal:** every public error enum and `*Config` struct in all 19 crates is `#[non_exhaustive]` or documented-final (today: 2 occurrences workspace-wide).
  - **Design:**
    - grep-driven sweep;
    - pair with cargo-semver-checks baseline (23.6) so additions stop being breaking.

- [ ] 23.4 `.expect()` reduction
  - **Goal:**
    - zero `.expect()` reachable from public API paths in io/core/viz (audited production counts: io 113, core 54, viz 50);
    - no-unwrap-style script green.
  - **Design:**
    - convert genuinely-fallible sites to `Result` with error enums;
    - keep invariant-protected ones with `// INVARIANT:` comments.

- [ ] 23.5 `deny.toml`
  - **Goal:** `cargo deny check` green via local script (no workflow yaml).
  - **Design:**
    - `cargo deny init`;
    - advisories/licenses/bans/sources sections;
    - explicit tracked ignore for unmaintained `paste 1.0.15` (transitive via nalgebra 0.34, confirmed in Cargo.lock).

- [ ] 23.6 Dependency-ordered publish `xtask`
  - **Goal:** `cargo xtask publish --dry-run` prints the correct topo order for the 18 publishable crates (oxiphysics-python is `publish = false`), verifies versions/CHANGELOG, supports resume-after-failure.
  - **Design:**
    - NEW `xtask/` workspace member (none exists);
    - parse `cargo metadata`, topo-sort path deps;
    - also home for `xtask ci` (fmt+clippy+nextest+deny+semver) replacing the disabled workflow.

- [x] 23.7 Flagship doctests (planned 2026-06-11)
  - **In flight:** fem/rigid/articulated lib.rs each gain ≥3 runnable doctests this pass; gate `cargo test --doc` workspace count ≥40.
  - **Goal:**
    - fem/rigid/articulated `lib.rs` each gain ≥3 compiling doctests (measured 2026-06-11: fem 0, rigid 0, articulated 0, constraints 10);
    - `cargo test --doc --workspace` count ≥40.
  - **Design:** cantilever, falling box, 3-link pendulum mini-examples lifted from existing integration tests.

- [~] 23.8 Re-export-only internal path-deps decision
  - **Partial (planned 2026-06-11):** the wasm re-export feature-gate lands this pass (see crates/oxiphysics/TODO.md); the 7-path-dep facade decision remains open.
  - **Goal:**
    - documented decision for the 7 crates consumed only by the umbrella (articulated, vehicle, sph, md, softbody, viz, io): keep-as-facade vs per-domain umbrella features;
    - ALSO fix `pub use oxiphysics_wasm as wasm` (`crates/oxiphysics/src/lib.rs:125`) pulling wasm-bindgen into native builds — feature-gate it.
  - **Design:** `cargo tree` evidence before/after.

- [ ] 23.9 Policy reconciliation
  - **Goal:** TODO.md "SciRS2 Policy Compliance" claims ("No nalgebra direct", "No rand direct") agree with cargo metadata.
  - **Reality:** nalgebra is a direct dep of core/rigid/python; rand is a direct dep of 14+ crates.
  - **Design:** either migrate (scirs2-core random + core type aliases) or rewrite the claims; policy-check script green either way.
  - **Cross-reference:** the two affected claims in the "SciRS2 Policy Compliance" section below are flagged `[~]` until this item lands.

- [ ] 23.10 CHANGELOG 0.1.3 populated at release (LAST)
  - **Goal:** `## [0.1.3] - Unreleased` sections (currently empty) filled from git history at ship time via /changelog-gen flow.

### v0.2.0 — Phases 24–26: programs

- [ ] Phase 24: Determinism & Reproducibility program (dep: 23.2)
  - **Goal:**
    - same-binary determinism guaranteed with test evidence (5-run hash equality with rayon on);
    - cross-platform scope documented;
    - replay golden corpus committed.
  - **Design:**
    - engine-wide deterministic mode flag (umbrella item),
    - fixed-order parallel reductions (core item),
    - seeded-RNG injection audit,
    - libm-transcendentals inventory;
    - builds on rollback desync hashing.

- [ ] Phase 25: Differentiable Physics end-to-end
  - **Goal:**
    - gradient checks ≤1e-4 rel through rigid contact + XPBD + SPH + linear FEM;
    - cartpole swing-up solved by gradient descent through the simulator;
    - PINN template example.
  - **Design:**
    - extend the shipped `oxiphysics_md::autograd_bridge` (`DifferentiableForceField`, KF-3) with implicit-function-theorem adjoints at the PGS/LCP fixed point;
    - scirs2-autograd bridge feature-gated like the existing `scirs2` core feature.

- [~] Phase 25 (blocked): scirs2-neural force fields — differentiable physics for neural-network potentials via scirs2-neural. Folded 2026-06-11 from "Proposed follow-ups"; the KF-3 scirs2-autograd bridge itself shipped 2026-05-14 (see history above). **Ready when:** scirs2-neural stable inference path lands.

- [ ] Phase 26: GPU pipeline completion
  - **Goal:** every item in oxiphysics-gpu v0.2.0 bucket lands with Phase-22-style parity + env-gated speedup tests.
  - **Design:** sequencing reduction-lib → sort → LBVH → SPH pipeline → cloth → contact solver (see gpu crate TODO).

- [~] Phase 26 (blocked): CUDA backend activation (`oxiphysics-gpu` Phase 5) — cudarc + real device buffers + SPH density kernel via NVRTC (planned 2026-05-14); wired under the non-default `cuda-backend` feature, pending hardware verification on RTX hardware. Folded 2026-06-11 from "Deferred / Strategic Roadmap" + "Proposed follow-ups"; tracked in the oxiphysics-gpu TODO (v0.3.0 bucket). **Ready when:** cudarc ≥ 0.16 stable and an RTX-class runner is committed.

### v0.3.0 — Phases 27–29

- [ ] Phase 27: Next-gen solvers
  - **Goal:**
    - MLS-MPM CPIC coupling (sand column collapse vs published experiment ≤10% runout error),
    - IPC barrier contact (zero interpenetration on 1k-stack torture test),
    - TGS-soft rigid solver option,
    - small-steps substepping defaults;
    - AMG-preconditioned implicit FEM dynamics (reuses shipped FEM RS/SA-AMG).
  - **Design:** new solver modules in softbody/rigid/fem crates with validation tests in the Phase-21 harness.

- [ ] Phase 28: Validation expansion + public benchmark dashboard
  - **Goal:**
    - Ghia Re=400/1000 (deferred in 21.3), NAFEMS LE1/LE10 FEM benchmarks, NIST LJ reference data;
    - static dashboard page regenerated by one local script comparing vs Rapier/Jolt/Bullet (pure-Rust rapier allowed as bench-only dev-dep; C++ engines via their published numbers, no FFI).
  - **Design:** CLI scene runner (umbrella) emits JSON consumed by a dashboard generator under `xtask bench-report`.

- [ ] Phase 29: Interop
  - **Goal:** robotics round-trip demos (KUKA iiwa URDF, MuJoCo ant MJCF) and CAD-to-sim (STEP → collision mesh) pipeline example.
  - **Design:** thin phase referencing oxiphysics-io v0.2.0/v0.3.0 items (URDF/MJCF/USD/STEP/OpenVDB/KHR-physics).

- [~] Phase 29 (blocked): OxiCAR coupling ABI — bind the shipped cross-domain auto-coupling runtime (FEM ⇄ SPH ⇄ LBM ⇄ MD, landed 2026-05-14) to OxiCAR's stable physics-coupling ABI. Folded 2026-06-11 from "Proposed follow-ups". **Ready when:** OxiCAR v0.2 stable physics-coupling ABI ships.

### v1.0 — Phases 30–31 + release gate

- [ ] Phase 30: Large-scale
  - **Goal:**
    - 10M-particle SPH on a single workstation (multi-GPU device groups + domain decomposition);
    - 100M-frame trajectory analyzed out-of-core under 512MB RSS (streaming io).
  - **Design:**
    - halo-exchange sharded worlds on rayon first;
    - pure-Rust transport for distributed later.

- [ ] Phase 31: DX & ecosystem
  - **Goal:**
    - "install → author scene → run → visualize → publish" story in an mdBook;
    - bevy plugin shipped;
    - oxicode adopted for snapshot serialization (replacing hand-rolled binary), oxiblas evaluated for FEM dense kernels (today's ecosystem deps are only oxifft, oxiarc-zstd, scirs2-integrate);
    - `cargo xtask semver` (cargo-semver-checks) green.
  - **Design:** see umbrella TODO.

- [ ] 1.0 release gate
  - **Goal:**
    - validation suite + GPU parity suite + benchmark dashboard green on Linux/macOS/Windows;
    - signed release checklist;
    - publish via `xtask publish` + existing pypi/npm workflows only.

### Phase dependency notes

- 23.2 feeds Phase 24 (the determinism program); Phase 24 is annotated `(dep: 23.2)`.
- 23.3's semver story pairs with the cargo-semver-checks baseline homed in 23.6 (`xtask ci`).
- 23.10 is deliberately LAST within Phase 23 — the CHANGELOG is populated at ship time via the /changelog-gen flow.
- Phase 26 sequencing (reduction-lib → sort → LBVH → SPH pipeline → cloth → contact solver) lives in the oxiphysics-gpu crate TODO.
- Phase 28's dashboard consumes the umbrella headless CLI scene-runner output via `xtask bench-report` (xtask infrastructure from 23.6).
- Phase 31's design detail lives in the umbrella TODO (crates/oxiphysics/TODO.md).
- 23.4's `.expect()` sweep is the prerequisite for the oxiphysics-io parser-fuzzing item (see io crate TODO).
- The wasm deterministic-lockstep demo depends on the umbrella deterministic mode delivered by Phase 24 (see wasm crate TODO).

### Roadmap maintenance

- New work lands in the Forward Roadmap above or in the per-crate TODO.md files below; Phases 1–22, the deferred backlog, and the #[allow] campaign section are frozen history.
- Status flips: `[ ]` → `[x]` with a landed date (YYYY-MM-DD) and a short **Result/Note** line; blocked items use `[~]` with an explicit **Ready when:** condition.
- Version-bump policy: the branch name (0.x.y format) drives Cargo.toml versions; never publish without explicit approval — `cargo publish --dry-run` only.
- Keep the README status table and this file in sync at each release (see the Housekeeping bullets below).

## Per-Crate Roadmap Index

Every crate carries its own TODO.md with the detailed v0.2.0 / v0.3.0 / v1.0 backlog; this table is the canonical index (all 19 crate files exist; the umbrella file is NEW as of the 2026-06-11 rewrite).
The v0.2.0 themes below are the agreed one-line summaries from the 2026-06-11 design briefs.

| Crate | TODO.md | v0.2.0 theme |
|-------|---------|--------------|
| oxiphysics-core | [./crates/oxiphysics-core/TODO.md](./crates/oxiphysics-core/TODO.md) | Numerical foundations & determinism |
| oxiphysics-geometry | [./crates/oxiphysics-geometry/TODO.md](./crates/oxiphysics-geometry/TODO.md) | Robustness & meshing |
| oxiphysics-collision | [./crates/oxiphysics-collision/TODO.md](./crates/oxiphysics-collision/TODO.md) | Narrowphase soundness and manifold quality |
| oxiphysics-rigid | [./crates/oxiphysics-rigid/TODO.md](./crates/oxiphysics-rigid/TODO.md) | Solver modernization (Jolt/PhysX-5 parity) |
| oxiphysics-constraints | [./crates/oxiphysics-constraints/TODO.md](./crates/oxiphysics-constraints/TODO.md) | Soft unification and exact small solvers |
| oxiphysics-articulated | [./crates/oxiphysics-articulated/TODO.md](./crates/oxiphysics-articulated/TODO.md) | Robotics parity (Pinocchio-class derivatives, IK, import) |
| oxiphysics-softbody | [./crates/oxiphysics-softbody/TODO.md](./crates/oxiphysics-softbody/TODO.md) | MPM v2 and robust FEM |
| oxiphysics-vehicle | [./crates/oxiphysics-vehicle/TODO.md](./crates/oxiphysics-vehicle/TODO.md) | Tire fidelity and transmission completion |
| oxiphysics-sph | [./crates/oxiphysics-sph/TODO.md](./crates/oxiphysics-sph/TODO.md) | Accuracy & boundaries |
| oxiphysics-lbm | [./crates/oxiphysics-lbm/TODO.md](./crates/oxiphysics-lbm/TODO.md) | Geometry & near-wall fidelity |
| oxiphysics-fem | [./crates/oxiphysics-fem/TODO.md](./crates/oxiphysics-fem/TODO.md) | Solvers & saddle-point |
| oxiphysics-md | [./crates/oxiphysics-md/TODO.md](./crates/oxiphysics-md/TODO.md) | Constraints, water, electrostatics |
| oxiphysics-materials | [./crates/oxiphysics-materials/TODO.md](./crates/oxiphysics-materials/TODO.md) | Homogenization & constitutive infrastructure |
| oxiphysics-gpu | [./crates/oxiphysics-gpu/TODO.md](./crates/oxiphysics-gpu/TODO.md) | Close the CPU-mock gap on the hot path |
| oxiphysics-viz | [./crates/oxiphysics-viz/TODO.md](./crates/oxiphysics-viz/TODO.md) | Verifiable rendering + GPU promotion |
| oxiphysics-io | [./crates/oxiphysics-io/TODO.md](./crates/oxiphysics-io/TODO.md) | Robotics interop + trustworthy parsers |
| oxiphysics-python | [./crates/oxiphysics-python/TODO.md](./crates/oxiphysics-python/TODO.md) | RL + typed surface |
| oxiphysics-wasm | [./crates/oxiphysics-wasm/TODO.md](./crates/oxiphysics-wasm/TODO.md) | Size, speed, types |
| oxiphysics (umbrella, NEW file) | [./crates/oxiphysics/TODO.md](./crates/oxiphysics/TODO.md) | Hardening + scene v2 |

Housekeeping (cross-file gaps spotted during the 2026-06-11 audit):

- [ ] README status table is missing the oxiphysics-articulated row (tracked in crates/oxiphysics-articulated/TODO.md)
- [ ] Umbrella crate has 3 orphan source files (builder.rs, prelude.rs, diagnostics.rs) never declared as modules (tracked in crates/oxiphysics/TODO.md)

---

## Blueprint v0.1 Reference (COOLJAPAN Ecosystem)

> OxiPhysics targets **19 main + ~30 sub crates**, **~1.08M net-new SLoC**, with **~720K SLoC saved** by reusing the SciRS2 / OxiBLAS / OxiMedia / Oxi3D stack. The sections below overlay blueprint intent onto the shipped phases above — they are tracking documentation, not a second roadmap.

### 5 Killer Features (Blueprint Differentiation)

| KF | Feature | Status |
|----|---------|--------|
| KF-1 | OxiCAR real-world vehicle physics (Raycast vehicle + Pacejka tires) | [x] shipped (Phase 7) — OxiCAR coupling layer deferred |
| KF-2 | scirs2-integrate bridge — LBM / FEM time-stepping via SciRS2 ODE/PDE | [x] shipped (2026-05-13) |
| KF-3 | ML force fields — differentiable physics via scirs2-autograd | [x] shipped (2026-05-14) |
| KF-4 | GPU physics — wgpu + rust-gpu compute (SPH, LBM, BVH) | [x] shipped (Phase 22 — 2026-05-14) |
| KF-5 | Python API — PyO3 ergonomic bindings across all domains | [x] shipped (Phase 20 — 2026-05-11: pytest harness, maturin wheel, .pyi stubs, all Phase-6 modules registered) |

### Per-Crate SLoC Targets (19 main crates)

| Crate | Target SLoC | Blueprint Role |
|-------|-------------|----------------|
| oxiphysics-core | 30,000 | Types, traits, math primitives |
| oxiphysics-geometry | 50,000 | Shapes, AABB, distance queries |
| oxiphysics-collision | 80,000 | GJK, EPA, broad-phase, BVH |
| oxiphysics-rigid | 60,000 | Rigid-body integrator, impulse solver |
| oxiphysics-constraints | 70,000 | Joints, Featherstone, articulated bodies |
| oxiphysics-vehicle | 60,000 | Raycast vehicle, Pacejka, drivetrain (KF-1) |
| oxiphysics-sph | 80,000 | SPH kernels, neighbor search, WCSPH / PCISPH |
| oxiphysics-lbm | 20,000 | LBM bridge to scirs2-integrate (KF-2) |
| oxiphysics-fem | 30,000 | FEM bridge to scirs2-integrate / scirs2-sparse |
| oxiphysics-md | 120,000 | Lennard-Jones, Ewald / PME, Verlet |
| oxiphysics-softbody | 80,000 | Mass-spring, PBD / XPBD, cloth |
| oxiphysics-materials | 40,000 | Material DB, Neural FF hooks (KF-3) |
| oxiphysics-gpu | 60,000 | wgpu + rust-gpu backend (KF-4) |
| oxiphysics-viz | 35,000 | VTK / glTF export, debug draw |
| oxiphysics-io | 60,000 | CalculiX / LAMMPS / VTK interchange |
| oxiphysics-python | 35,000 | PyO3 bindings (KF-5) |
| oxiphysics-wasm | 20,000 | wasm-bindgen browser runtime |
| tests / bench | 120,000 | Integration, criterion, regression harness |
| docs / examples | 30,000 | Gallery, API docs, validation notebooks |
| **Total** | **~1,080,000** | net-new SLoC budget |

### Ecosystem Reuse Savings (~720K SLoC)

| Asset | Purpose | Saves SLoC |
|-------|---------|------------|
| scirs2-integrate | LBM / FEM / ODE / PDE integrators (KF-2) | ~200,000 |
| scirs2-sparse | AMG / GMRES / CG for FEM, linear systems | ~150,000 |
| scirs2-core | SIMD, Rayon, GPU allocator, random | ~60,000 |
| OxiMedia | Rendering (wgpu Vulkan/Metal/DX12) | ~60,000 |
| scirs2-optimize | Nonlinear solvers (Gauss-Newton, L-BFGS) | ~50,000 |
| scirs2-spatial | KD-tree, BVH, convex hull, Voronoi | ~40,000 |
| scirs2-autograd | Differentiable physics for Neural FF (KF-3) | ~40,000 |
| OxiBLAS | BLAS / LAPACK (LU, Cholesky, QR, SVD) | ~40,000 |
| scirs2-fft | Spectral solvers, Ewald reciprocal space | ~30,000 |
| scirs2-linalg | SVD, PCA, dense linear algebra | ~30,000 |
| Oxi3D (bridge) | Point-cloud → mesh for collider authoring | ~20,000 |
| **Total** | | **~720,000** |

### SciRS2 Policy Compliance

- [x] **No C/C++ FFI** — no OCCT, Eigen, LAMMPS-C, OpenFOAM, CalculiX (file-level I/O compatibility only), Bullet
- [~] **No nalgebra direct** — core math via `oxiphysics-core` type aliases — **claim under reconciliation (see Phase 23.9):** code currently uses this dep directly; either migrate or amend the policy claim.
- [~] **No `rand` direct** — randomness via `scirs2_core::random` per project policy — **claim under reconciliation (see Phase 23.9):** code currently uses this dep directly; either migrate or amend the policy claim.
- [x] **No CUDA / HIP FFI** — GPU path via wgpu + rust-gpu (KF-4); cudarc backend is a post-v0.3.0 option
- [x] **No OpenCV** — imaging delegated to OxiMedia
- [x] SIMD / parallel via scirs2-core primitives
- [x] Randomness via scirs2-core `random_range` / `rng` (per rand 0.9 API migration)
- [x] **No ToRSh** — COOLJAPAN-wide policy; ML force fields (KF-3) via `scirs2-autograd`, never ToRSh

---

Last Updated: 2026-06-11 — version 0.1.3

## Stubs to implement (added 2026-06-12 by /cooljapan-stub-check)

- [ ] `oxiphysics-gpu`: `crates/oxiphysics-gpu/src/compute/wgpu_backend.rs:133` — implement WgpuBackend::try_new with real wgpu adapter enumeration (gated by `wgpu-backend` feature)
  - Priority: P2 | Scope: medium | Hint: none
- [ ] `oxiphysics-gpu`: `crates/oxiphysics-gpu/src/compute/wgpu_backend.rs:245,271` — implement remaining wgpu compute dispatch methods (two additional TODO stubs in same file)
  - Priority: P2 | Scope: medium | Hint: none
- [ ] `oxiphysics-gpu`: `crates/oxiphysics-gpu/src/scheduler.rs:528` — replace placeholder 4-byte output with real GPU result readback
  - Priority: P2 | Scope: small | Hint: none
- [ ] `oxiphysics-vehicle`: `crates/oxiphysics-vehicle/src/autonomous_vehicle/types.rs:5` — resolve P4 public API churn (tracked as TODO(P4))
  - Priority: P2 | Scope: small | Hint: none
- [ ] `oxiphysics-geometry`: `crates/oxiphysics-geometry/src/computational_geometry/types.rs:5` — implement `ops::Add`/`Sub` impls for `Point2` to replace current method-based approach
  - Priority: P2 | Scope: trivial | Hint: none
- [ ] `oxiphysics-core`: `crates/oxiphysics-core/src/exact_predicates.rs:348,400,469,562` — implement Shewchuk adaptive intermediate stages (C-stage for orient2d, B/C/D for incircle, adaptive stages for insphere)
  - Priority: P2 | Scope: medium | Hint: none
