# oxiphysics-wasm

**Status: Partial** — engine logic and web-friendly API complete; wasm-bindgen not yet linked.

WebAssembly frontend layer for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.  
Version: **0.1.0** | Updated: **2026-04-06**

---

## Architecture

This crate provides a **self-contained physics engine** with a web-oriented API surface.  
It does **not** depend on any other `oxiphysics-*` crate — all engine logic is embedded here.  
`wasm-bindgen` is **not** a current dependency; `.wasm` + JS glue generation is planned for 0.2.0.

> **Note:** WASM API and physics engine logic are complete.  
> wasm-bindgen integration is planned for **0.2.0**.

---

## Public API Surface

1,276 public items · 826 tests · 0 stubs

### Modules

| Module | Description |
|---|---|
| `analytics_bridge` | Performance metrics and analytics export |
| `body_query` | Query API for rigid body state |
| `constraint_bridge` | Constraint/joint configuration bridge |
| `debug_tools` | Debug info extraction and inspection |
| `engine` | Core simulation engine (submodules: step, state, world, …) |
| `error` | `Error` and `Result` types |
| `events` | Collision/sensor event queue |
| `fluid_bridge` | SPH/LBM fluid state bridge |
| `io_bridge` | Scene serialization in/out |
| `js_api` | JavaScript-facing API surface (future `#[wasm_bindgen]` targets) |
| `material_bridge` | Material parameter bridge |
| `math_helpers` | WASM-friendly math utilities (submodules) |
| `particle_system` | Particle update and management |
| `physics_config` | Global simulation configuration |
| `renderer` | Render data extraction for WebGL |
| `sim_controls` | Play/pause/reset/step controls |
| `simulation_api` | Top-level simulation entry points |
| `types` | Shared primitive types |

### Key exported types

`WasmPhysicsEngine`, `Error`, `Result`;  
`BodyState`, `ColliderConfig`, `ColliderShapeType`, `ContactResult`, `DebugInfo`,
`QuatWasm`, `RaycastResult`, `RigidBodyConfig`, `SimulationConfig`,
`TransformWasm`, `Vec3Wasm`

---

## Roadmap

| Milestone | Target |
|---|---|
| Self-contained engine + web API complete | ✅ 0.1.0 |
| wasm-bindgen `#[wasm_bindgen]` annotation pass | 🔲 0.2.0 |
| wasm-pack / wasm-bindgen-cli build pipeline | 🔲 0.2.0 |
| npm package publish | 🔲 0.2.0 |
| WGPU/WebGL render integration | 🔲 0.3.0 |

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
