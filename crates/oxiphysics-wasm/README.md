# oxiphysics-wasm

## Status: Stable (v0.1.3)

Full wasm-bindgen 0.2 bindings: 929 `#[wasm_bindgen]` annotations across
27 bridge files exposing the entire physics surface to JavaScript.
Build via `cargo build --target wasm32-unknown-unknown -p oxiphysics-wasm`,
then run `wasm-bindgen` / `wasm-pack` over the resulting artifact for
JS glue generation.

WebAssembly frontend layer for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.  
Version: **0.1.3** | Updated: **2026-06-06**

---

## Architecture

This crate provides a **self-contained physics engine** annotated with
`#[wasm_bindgen]` for direct JavaScript interop. It does **not** depend
on any other `oxiphysics-*` crate — all engine logic is embedded here.
The crate compiles cleanly to `wasm32-unknown-unknown` and exposes a
flat, JS-friendly API surface (primitive arguments, flat-array returns,
no Rust references across the boundary).

---

## Quick start (JavaScript)

```javascript
import init, {
    WasmPhysicsEngine,
    WasmDebugDraw,
    WasmWorld,
} from "./pkg/oxiphysics_wasm.js";

await init();

// Create an engine with Earth gravity.
const engine = new WasmPhysicsEngine(0.0, -9.81, 0.0);

// Add a 1 kg sphere at height 10 m.
const ball = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
engine.add_sphere_collider(ball, 0.5);

// Add a static ground plane at y=0.
const ground = engine.add_static_body(0.0, 0.0, 0.0);
engine.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);

// Simulate for one second (with internal substeps).
engine.step(1.0);

const [x, y, z] = engine.get_position(ball);
console.log(`Ball fell to y=${y.toFixed(3)}`);
```

---

## Public API Surface

929 `#[wasm_bindgen]` annotations · ~27k Rust SLoC · 922 host tests · 0 stubs

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
| wasm-bindgen `#[wasm_bindgen]` annotation pass | ✅ 0.1.1 |
| wasm-pack / wasm-bindgen-cli build pipeline | 🔲 0.2.0 |
| npm package publish | 🔲 0.2.0 |
| WGPU/WebGL render integration | 🔲 0.3.0 |

---

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
