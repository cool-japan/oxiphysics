# OxiPhysics Browser Demos

Four runnable WebGL2 / Canvas2D demos for the OxiPhysics WASM engine.

## Build

```sh
# From the crate root (crates/oxiphysics-wasm/)
wasm-pack build --target web --out-dir pkg-web --release
```

This produces `pkg-web/oxiphysics.js` (the ES module entry point) and
`pkg-web/oxiphysics_bg.wasm`.

## Serve

```sh
# From the repo root — serves all static files including web-demos/
python3 -m http.server -d . 8080
# Then open: http://localhost:8080/crates/oxiphysics-wasm/web-demos/
```

Or serve from the crate directory:

```sh
python3 -m http.server -d crates/oxiphysics-wasm 8080
# Open: http://localhost:8080/web-demos/
```

## Demos

| Demo | Directory | Description |
|------|-----------|-------------|
| Character Controller | `character/` | Kinematic capsule on 64×64 sine heightfield. WASD + mouse-look. |
| Rope Ragdoll | `rope/` | 32-segment Verlet rope, pinned at top, mouse drag impulse. |
| SPH Fluid Pool | `sph-pool/` | 600 SPH particles, gravity, click-to-splash. |
| NavMesh Path-Finding | `navmesh/` | 10×10 grid navmesh, A* funnel smoothing, click to set goal. |

## WASM API used

| Demo | WASM type |
|------|-----------|
| Character | `WasmCharacterControllerJs` (thin `#[wasm_bindgen]` wrapper around `WasmCharacterController`) |
| Rope | `WasmRopeJs` (thin `#[wasm_bindgen]` wrapper around `WasmRope`) |
| SPH pool | `WasmSphSim` (directly `#[wasm_bindgen]` in `engine/sph.rs`) |
| NavMesh | `WasmNavMeshJs` (thin `#[wasm_bindgen]` wrapper around `WasmNavMesh`) |

## Notes

- All demos are ESM-only (`.mjs`), no bundler required.
- No external JS dependencies (no Three.js, no Babylon.js).
- Each demo is ≤300 lines of JavaScript.
- Requires a browser with WebGL2 support (Chrome 56+, Firefox 51+, Safari 15+).
