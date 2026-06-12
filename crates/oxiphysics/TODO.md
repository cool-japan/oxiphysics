# oxiphysics (umbrella crate) TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 18,452 SLoC | 341 tests + 9 examples + criterion bench

Legend: `[x]` shipped · `[ ]` planned · `[~]` deferred/blocked

> Constraints: Pure-Rust default features; no new workflow yamls (only the existing `pypi-publish.yml` / `npm-publish.yml`) — all CI/dashboard items are local scripts or `xtask` subcommands; ecosystem substitutions apply (oxicode not bincode, oxiarc-* not zip/flate2/zstd-c, oxiblas, oxifft, scirs2-*, OxiZ not z3).

## Completed (v0.1.0 – v0.1.3)

- [x] Sub-crate facade — 17 `pub use` re-exports for single-dependency access (`src/lib.rs`):
  core, geometry, collision, rigid, articulated, constraints, vehicle, sph, lbm, fem, md, softbody, materials, gpu, viz, io, wasm
  (the unconditional `wasm` re-export is the v0.2.0 priority fix below)
- [x] `pipeline` — full simulation pipeline (forces → broadphase → narrowphase → solve → integrate)
  - `PhysicsConfig`-driven (gravity, solver iterations, sleeping thresholds); in-crate two-run determinism regression test
- [x] Runtime & scene layer — 9 modules (root TODO Phases 13–15):
  - `force_field` — spatial force fields (gravity wells, vortex, wind, explosion)
  - `event_bus` — physics event publish/subscribe
  - `replay` — deterministic simulation replay (record → serialise → replay)
  - `query` — raycasting, sphere sweeps, overlap tests
  - `scene` — declarative scene description, JSON round-trip + fluent `SceneBuilder`
  - `snapshot` — world-state snapshots with delta tracking and ring-buffer history
  - `trigger` — sensor volumes with enter/exit/stay events
  - `animation` — keyframe tracks (Vec3 lerp + quaternion SLERP)
  - `material_table` — runtime material interaction table with combine rules and per-pair overrides
- [x] Physics utilities layer — 9 modules (root TODO Phases 16–18):
  - `debug_draw` — renderer-agnostic debug draw command buffer
  - `contact_cache` — persistent contact pair cache with warm-start impulses
  - `buoyancy` — Archimedes buoyancy + viscous drag for fluid volumes
  - `scheduler` — priority-based physics step budget allocation per island
  - `spatial_grid` — uniform spatial hash grid for neighbourhood queries
  - `lod` — Level-of-Detail simulation tier management
  - `noise` — procedural 3D value noise and fractal Brownian motion
  - `interpolator` — smooth damp, exponential decay, spring followers, lerp utilities
  - `telemetry` — per-step stats, rolling averages, CSV export
- [x] Advanced dynamics & tooling layer — 8 modules (root TODO Phase 19):
  - `character` — kinematic capsule controller with sweep-and-slide, step-up, slope handling
  - `rope` — distance constraints with Verlet integration + Gauss-Seidel projection
  - `xpbd` — Extended Position-Based Dynamics with compliance and Lagrange multipliers
  - `ik` — FABRIK and 2-bone analytic IK with joint cone limits
  - `profiler` — hierarchical RAII scopes; `to_csv` / `to_json` / `to_folded_stacks` exports
  - `aero` — velocity-squared drag and airfoil lift/drag forces
  - `navmesh` — A* pathfinding with funnel-algorithm smoothing
  - `rollback` — snapshot rollback, deterministic lockstep input buffer, FNV-1a desync hash
- [x] Cross-domain coupling runtime (OxiCAR, Blueprint KF-1) — `coupling` module
  - `DomainCoupler` framework linking FEM/SPH/LBM/MD regions at shared interfaces; `CouplingRuntime`, `CouplingDomain`, `InterfaceSite/State/Force`, `FemSphCoupler`, `MdContinuumAdapter`
  - common types re-exported at the crate root via `#[doc(inline)]`
- [x] `perf_regression` — performance regression testing infrastructure
  - committed baselines under `tests/regression_baselines/` driven by `tests/regression_harness.rs`
- [x] Validation & integration suite — 341 tests total
  - `tests/validation_{rigid,sph,lbm,fem,fem_amg,md_lj,ewald_madelung}.rs` (Phase-21 validation harness)
  - `tests/articulated_smoke.rs`, `tests/coupling_smoke.rs`, `tests/fem_diagnostic.rs`
  - in-crate physics-law tests: momentum/KE conservation (elastic + inelastic), projectile range/height, pendulum period, restitution e²H, freefall kinematics, torque → angular velocity, center-of-mass invariance, LBM density conservation, FEM cantilever deflection, SPH boundedness, MD energy finiteness
- [x] 9 runnable examples (`cargo run --example <name>`):
  - `falling_boxes` — rigid-body stacking
  - `cantilever` — FEM beam under load
  - `dam_break` — SPH free-surface flow
  - `lbm_channel` — lattice-Boltzmann channel flow
  - `argon_md` — Lennard-Jones argon molecular dynamics
  - `cloth_simulation` — soft-body cloth
  - `vehicle_dynamics` — vehicle simulation
  - `coupled_fsi` — cross-domain fluid-structure coupling
  - `benchmark_demo` — performance demonstration
- [x] Criterion benchmark — `benches/physics_bench.rs` (`harness = false`)
- [x] Stability policy — every public API annotated with a stability level (`core::stability::HasStability`); stable APIs follow semver

## v0.2.0 — Hardening + scene v2

Headline items: deterministic mode, scene format v2, CLI scene runner, plugin registry — plus the priority wasm re-export fix.

### Priority fix
- [x] **PRIORITY:** Feature-gate the wasm re-export (planned 2026-06-11)
  - **Goal:** default native build has no wasm-bindgen in `cargo tree` (today `lib.rs:125` re-exports oxiphysics-wasm unconditionally); decision recorded for the 7 umbrella-only path-deps (root Phase 23.8: articulated, vehicle, sph, md, softbody, viz, io — keep-as-facade vs per-domain umbrella features).
  - **Design:** `cargo tree` evidence before/after.
  - **Files:** Cargo.toml (oxiphysics-wasm → optional, `[features] wasm = ["dep:oxiphysics-wasm"]`), src/lib.rs (#[cfg(feature = "wasm")] on the re-export + any intra-doc links)
  - **Tests:** `cargo tree -p oxiphysics` default shows zero wasm-bindgen; `cargo check -p oxiphysics --features wasm` still green; full crate test suite
  - **Risk:** internal `crate::wasm` references or doc links — grep first; workspace --all-features runs will still compile wasm (unchanged from today)

### Determinism (headline; couples with root Phases 23.2 / 24)
- [ ] Engine-wide deterministic mode flag
  - **Goal:** `DeterminismConfig` on builder/pipeline forces fixed iteration order, seeded RNG, ordered parallel reductions; 10k-step state hash identical across 5 runs with rayon enabled (hash via existing `rollback.rs:186` FNV path).
- [ ] Cross-target repro evidence
  - **Goal:** recorded x86_64 vs aarch64 hash matrix for rigid+xpbd scenes via local script; failures triaged to libm/ordering (workspace `mul_add` count is 1 — fem topology_opt — so FMA is not the main suspect).
  - Refs: root Phase 23.2, feeding root Phase 24.

### Scene v2 & multi-physics authoring (headline)
- [ ] Declarative multi-physics scene format v2
  - **Goal:** dam-break-on-elastic-gate (SPH+FEM) runs from pure JSON with cross-domain coupling described in data (`domains[]` + `couplings[]`); golden trajectory test.
  - **Design:** extends `scene.rs` and the shipped `coupling/` runtime (`DomainCoupler`, `CouplingRuntime`, `FemSphCoupler` all exist).
- [ ] Headless CLI scene runner
  - **Goal:** `oxiphysics-cli run scene.json --out vtu,dcd,telemetry.csv` executes all 9 example scenes headless; becomes the engine for the Phase-28 dashboard.
  - **Design:** new bin target or xtask subcommand; uses io writers + telemetry.
  - Refs: root Phase 28 (validation expansion + public benchmark dashboard).
- [ ] Co-simulation scheduler with per-domain subcycling
  - **Goal:** FEM(1x) + SPH(4x) coupled run conserves exchanged momentum to 1e-6/step with interpolated boundary exchange.
  - **Design:** marry `scheduler.rs` priority model with `CouplingRuntime` rate ratios.
- [ ] Plugin registry (custom forces/materials/solver stages as trait objects)
  - **Goal:** an external crate registers a custom force without forking (example included).
  - **Design:** registries generalizing the concrete `force_field` / `material_table` modules; object-safe `SolverStage` hooks in `pipeline.rs`.

### API hardening (root Phase 23)
- [ ] Flagship doctests for fem/rigid/constraints/articulated re-exports
  - **Goal:** each re-exported domain has ≥1 compiling doctest at the umbrella level; fem/rigid/articulated lib.rs currently have 0 examples.
  - Refs: root Phase 23.7.
- [ ] `#[non_exhaustive]` + config-struct sweep for umbrella public surface
  - **Goal:** all pub error enums/configs across the 33 modules non_exhaustive; cargo-semver-checks baseline recorded.
  - Refs: root Phases 23.3 (sweep) and 23.6 (semver baseline lives in xtask).
- [x] Wire or remove orphan source files (verified 2026-06-11) (planned 2026-06-11)
  - `src/builder.rs`, `src/prelude.rs`, `src/diagnostics.rs` exist in `src/` but have no `mod` declarations anywhere in `lib.rs` — they are not compiled into the crate.
  - **Goal:** zero unreferenced `.rs` files; if wired, the prelude becomes part of the surface the v1.0 semver freeze covers.
  - **Files:** src/lib.rs (`pub mod builder; pub mod diagnostics; pub mod prelude;`), the three orphan files (fix latent breakage if wiring reveals any)
  - **Tests:** smoke test per module (builder constructs a world; `use oxiphysics::prelude::*;` compiles; diagnostics timer round-trip); crate suite green
  - **Risk:** name collisions with existing root modules/re-exports; stale internal APIs — recon says they compile, verify

## v0.3.0 — Scale & ecosystem

- [ ] bevy plugin (feature-gated `bevy`, non-default)
  - **Goal:** bevy example with 1k cubes at 60 fps; `RigidBody`/`Collider` components synced to bevy `Transform`; pinned to one bevy minor; big ecosystem win.
- [ ] Large-world coordinates (origin shifting)
  - **Goal:** stable simulation 1e7 m from origin with jitter <1 mm at camera; shift events published through `event_bus`, snapshot/replay-safe.
- [ ] Chrome-tracing profiler export
  - **Goal:** `FrameReport::to_chrome_trace()` (today: `to_csv`/`to_json`/`to_folded_stacks` at `profiler.rs:97–118`) opens in Perfetto with nested scopes incl. GPU timestamp spans from oxiphysics-gpu.
- [ ] mdBook user guide
  - **Goal:** 12 chapters (quickstart, each domain, coupling, determinism, GPU, bindings) building clean via `mdbook build` local script; code blocks tested with `mdbook test`.
- [ ] Property-based API fuzzing
  - **Goal:** 10k random world-building op sequences (proptest — already a workspace dep) never panic/NaN across rigid/xpbd/query/snapshot.

## v1.0 — Production & validation

- [ ] Semver freeze of prelude + module surface
  - **Goal:** `cargo xtask semver` (cargo-semver-checks, local) green; 1.0 API audit doc; zero accidental breaks vs final 0.x.
- [ ] Ecosystem adoption pass (root Phase 31 — "Design: see umbrella TODO")
  - **Goal:** oxicode adopted for snapshot serialization (replacing hand-rolled binary); oxiblas evaluated for FEM dense kernels (today's ecosystem deps are only oxifft, oxiarc-zstd, scirs2-integrate); "install → author scene → run → visualize → publish" story lands in the mdBook.
- [ ] 1.0 acceptance gate
  - **Goal:** Phase-21 validation + GPU parity + golden-image + benchmark dashboard all green on Linux/macOS/Windows; release executed via publish xtask + existing two workflows only.

## Deferred / research track

- [~] Large-scale worlds (root Phase 30: 10M-particle SPH on a single workstation via multi-GPU domain decomposition; 100M-frame trajectory analyzed out-of-core under 512MB RSS) — Ready when: oxiphysics-gpu full-pipeline residency + multi-GPU device groups and oxiphysics-io streaming readers land.
- [~] Differentiable end-to-end pipeline at the umbrella level (root Phase 25: gradient checks through rigid contact + XPBD + SPH + linear FEM) — Ready when: the `oxiphysics_md::autograd_bridge` IFT-adjoint extension and the feature-gated scirs2-autograd bridge land in the domain crates.
- [~] GPU kernel spans inside umbrella flamegraphs — Ready when: oxiphysics-gpu wgpu timestamp-query scopes export into `FrameReport` (gpu crate v0.2.0); the v0.3.0 chrome-trace item then consumes them.

## Root-phase cross-reference

| Umbrella item | Root TODO phase |
|---|---|
| Feature-gate wasm re-export / path-dep decision | Phase 23.8 |
| DeterminismConfig + cross-target repro | Phases 23.2, 24 |
| Flagship doctests | Phase 23.7 |
| non_exhaustive sweep + semver baseline | Phases 23.3, 23.6 |
| CLI scene runner → benchmark dashboard | Phase 28 |
| bevy plugin, mdBook, oxicode/oxiblas adoption | Phase 31 |
| Large-scale worlds (deferred) | Phase 30 |
| Differentiable pipeline (deferred) | Phase 25 |
| 1.0 acceptance gate | root 1.0 release gate |
