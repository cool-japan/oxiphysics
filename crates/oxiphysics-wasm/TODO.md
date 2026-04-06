# oxiphysics-wasm TODO

Updated: 2026-04-06 | Version: 0.1.0

---

## Phase 1: Foundation ✅

- [x] Define error types (`error` module) and shared primitives (`types`)
- [x] Implement math helpers (`math_helpers` submodules: vectors, quaternions, transforms)
- [x] Add unit tests (826 tests passing)

## Phase 2: Engine Logic ✅

- [x] `engine` module — self-contained simulation engine (step/state/world submodules)
- [x] `physics_config` — global simulation parameters
- [x] `simulation_api` — top-level simulation entry points (`WasmPhysicsEngine`)
- [x] `sim_controls` — play/pause/reset/step controls
- [x] `body_query` — rigid body state queries
- [x] `particle_system` — particle update and management
- [x] `events` — collision/sensor event queue

## Phase 3: Bridge Modules ✅

- [x] `analytics_bridge` — performance metrics export
- [x] `constraint_bridge` — constraint/joint configuration bridge
- [x] `debug_tools` — debug info extraction (`DebugInfo`)
- [x] `fluid_bridge` — SPH/LBM fluid state bridge
- [x] `io_bridge` — scene serialization in/out
- [x] `js_api` — JavaScript-facing API surface
- [x] `material_bridge` — material parameter bridge
- [x] `renderer` — render data extraction for WebGL

## Phase 4: Documentation ✅

- [x] All 1,276 public items documented
- [x] Key type exports: `WasmPhysicsEngine`, `Vec3Wasm`, `QuatWasm`, `TransformWasm`, `RigidBodyConfig`, `SimulationConfig`, `ColliderConfig`, `ContactResult`, `RaycastResult`, `BodyState`, `DebugInfo`

## Phase 5: wasm-bindgen Packaging 🔲 (planned 0.2.0)

- [ ] Add `wasm-bindgen` dependency to `Cargo.toml`
- [ ] Annotate exported types and functions with `#[wasm_bindgen]`
- [ ] Configure `wasm-pack` build pipeline
- [ ] Generate JS/TS glue code
- [ ] Publish npm package

## Phase 6: Runtime Integration 🔲 (planned 0.2.0+)

- [ ] Web Worker / `SharedArrayBuffer` multi-threaded support
- [ ] WGPU/WebGL render integration
- [ ] Performance benchmarks vs. Rapier WASM
- [ ] Examples: browser-based rigid body demo, fluid demo
