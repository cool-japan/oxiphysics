# oxiphysics-viz TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 102,114 SLoC | 4,114 tests

> **Constraints:** Pure-Rust default features; the GPU path stays behind the feature-gated wgpu renderer. OpenXR (v1.0) is feature-gated non-default — its FFI loader mirrors the cudarc pattern; the default build remains pure Rust. Test/CI items are local scripts only — no new workflow yamls.

> **Two code-verified facts that shape this roadmap (2026-06-11):**
> 1. `GpuRenderer` core ALREADY EXISTS (`src/wgpu_renderer.rs` — Phong + particle + blit WGSL, `GpuRenderTarget`, CPU fallback, `to_hdr_framebuffer`). The open v0.2.0 work is COMPLETION (depth/MSAA/materials/camera), not creation.
> 2. The crate has NO `tests/` dir — all 4,114 tests are in-module. The golden-image harness is therefore the first v0.2.0 item: it creates the crate's first integration-test surface and locks every later rendering item.

### Performance acceptance targets at a glance (consolidated from Goals below)

| Item | Budget | Hardware / gating |
|---|---|---|
| GpuRenderer completion (v0.2.0) | 100k-triangle mesh @ 60 fps, 1080p offscreen | M-series; skip-not-fail headless |
| Instanced particles (v0.2.0) | 1M impostor spheres @ 60 fps | M-series |
| Screen-space fluid (v0.2.0) | 250k particles @ 30 fps, 1080p | M-series |
| Raymarched volume (v0.2.0) | 256³ volume @ 60 fps | M-series |
| GPU marching cubes (v0.2.0) | 128³ SDF → mesh < 10 ms | M-series; dep: gpu reduction lib |
| 3D LIC (v0.2.0) | 64³ field < 2 s | CPU |
| GPU picking (v0.2.0) | 1M-instance pick < 1 ms | ID buffer |
| OIT (v0.3.0) | 10k translucent particles, order-stable | golden-locked |
| Remote viewer (v0.3.0) | 10k-body scene @ 60 fps | localhost stream |
| OpenXR (v1.0) | stereo @ 72 fps | reference HMD; env-gated |

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation
- [x] Implement primary algorithms (rasterizer, post-processing, mesh gen)
- [x] Add integration tests (4,114 tests passing)
- [x] Performance benchmarks

### Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization

### Phase 4: Advanced Features
- [x] CPU software rasterizer (Phong, wireframe, framebuffer)
- [x] Post-processing pipeline (Bloom, DoF, SSAO, tone mapping, vignette)
- [x] Colormap & transfer functions
- [x] Streamline tracing and rendering
- [x] Stress/FEM viz (principal stress glyphs, von Mises)
- [x] Volume rendering and isosurface extraction
- [x] Particle renderer, effects, trails, instancing
- [x] Physics animation and real-time viz
- [x] Scientific plotting (chart, heatmap, statistical, uncertainty)
- [x] Medical imaging and molecular visualization
- [x] Fluid visualization
- [x] Neural rendering
- [x] Metaball and procedural texture
- [x] Topology, tensor, and phase-field visualization
- [x] Font / text / annotation rendering
- [x] Graph and network visualization
- [x] Glyph renderer, gizmos, debug overlay
- [x] Terrain renderer
- [x] VR visualization pipeline

### Post-0.1 features (shipped during v0.1.x)
- [x] GPU/wgpu backend (feature-gated) — `wgpu_renderer` (`GpuRenderer`, Phong+particle+blit WGSL shaders, `GpuRenderTarget`, CPU fallback)
- [x] WebAssembly canvas rendering path (`wasm_canvas` — `CanvasBuffer`, `CanvasRasterizer`, Porter-Duff blending)
- [x] Interactive picking / selection in rasterizer (`picking` — `Picker`, Möller–Trumbore, `SelectionSet`)
- [x] HDR framebuffer and wider color spaces (`hdr_framebuffer` — f32 RGBA, ACES/Reinhard/Filmic/AgX tone mapping)

### Verified module baseline the roadmap builds on (audited 2026-06-11)
- `src/wgpu_renderer.rs` — `GpuRenderer` + `GpuRenderTarget` + `to_hdr_framebuffer` (core exists; completion pending)
- `instancing.rs` — CPU instancing (parity reference for GPU instanced particles)
- `metaball.rs` — CPU fluid-surface path (visual reference for screen-space fluid)
- `volume_renderer.rs` + `transfer_functions` — CPU volume path (parity reference for GPU raymarch)
- `isosurface.rs` — full 256-config CPU marching cubes (parity reference for GPU MC)
- `flow_visualization.rs:459` — 2D LIC (3D LIC pending)
- `picking.rs` — Möller–Trumbore CPU picking (reference for GPU ID-buffer picking)
- `vr_visualization.rs` — CPU stereo model (basis for the OpenXR pure-Rust fallback)

## v0.2.0 — Verifiable rendering + GPU promotion

Sequencing: golden-image harness first (everything after it lands golden-locked) → GpuRenderer completion → instanced particles → screen-space fluid; GPU marching cubes waits on the oxiphysics-gpu reduction lib. Exit criteria: ≥12 golden scenes green for both CPU rasterizer and `GpuRenderer` (SSIM ≥0.99); GpuRenderer feature-complete (depth/MSAA/materials/camera); 1M instanced particles at 60 fps on M-series.

Definition of done for every rendering item from here on:
1. Golden-locked — at least one scene in the golden matrix exercises it (SSIM ≥0.99).
2. Named CPU parity/visual reference module documented in the item.
3. GPU tests skip-not-fail on headless machines; fps assertions env-gated, never default.
4. No new workflow yamls — harness invocations are local scripts.
5. No `unwrap()`; new files stay under 2,000 lines (splitrs when exceeded).

Out of scope for v0.2.0: windowing/event-loop ownership (the crate stays offscreen-first; swapchain/window integration remains the consumer's job) and any scene-graph rewrite.

### Test infrastructure first
- [ ] Headless offscreen golden-image test harness — FIRST item; crate currently has NO `tests/` dir
  - **Goal:** ≥12 golden scenes (mesh/particles/volume/plots) rendered by CPU rasterizer and `GpuRenderer`, SSIM ≥0.99 vs committed PNGs, env-var baseline refresh.
  - **Design:** offscreen `GpuRenderTarget` + `to_ldr_bytes()` already exist; SSIM comparator implemented in-crate (pure Rust) — no external image-diff tool.
  - **Files:** `tests/golden_images.rs` (NEW), committed baselines under `tests/golden/`; GPU half skip-not-fail on headless machines.
  - **Risk:** cross-platform float drift in the CPU rasterizer — SSIM threshold (not byte equality) is the mitigation; baseline refresh stays env-var-gated so drift is always an explicit decision.
  - **Candidate scene matrix (each maps to a shipped Phase-4 feature; GPU column fills in as the completion item below lands):**

    | # | Scene | Exercises | Backends |
    |---|---|---|---|
    | 1 | mesh_phong | CPU rasterizer Phong | CPU + GPU |
    | 2 | mesh_wireframe | wireframe path | CPU |
    | 3 | mesh_depth_msaa | depth + MSAA (new GPU pipeline) | GPU |
    | 4 | particles_instanced | `instancing.rs` / GPU instancing | CPU + GPU |
    | 5 | particle_trails | effects + trails | CPU |
    | 6 | volume_raymarch | `volume_renderer` + `transfer_functions` | CPU + GPU |
    | 7 | isosurface_mc | `isosurface.rs` 256-config MC | CPU |
    | 8 | streamlines | streamline tracing | CPU |
    | 9 | tensor_glyphs | stress/FEM glyphs (von Mises) | CPU |
    | 10 | chart_heatmap | scientific plotting | CPU |
    | 11 | statistical_uncertainty | statistical + uncertainty plots | CPU |
    | 12 | postfx_bloom_ssao | post-processing chain | CPU |
    | 13 | hdr_tonemap_aces | `hdr_framebuffer` tone mapping | CPU |

### GpuRenderer completion & GPU promotion
- [ ] wgpu render backend completion — note: core already exists (`src/wgpu_renderer.rs` `GpuRenderer` with Phong+particle+blit WGSL, CPU fallback, `to_hdr_framebuffer`)
  - **Remaining:** depth-tested triangle pipeline, MSAA, texture/material binding, camera uniform path.
  - **Goal:** 100k-triangle mesh at 60 fps, 1080p offscreen, on M-series.
  - **Files:** `src/wgpu_renderer.rs` (depth/MSAA/material/camera passes); **Tests:** golden scenes 1, 3, 4, 6 above in GPU mode.
  - **Risk:** wgpu major-version churn (the `pass` reserved-word episode in wgpu 29 is the precedent) — keep WGSL identifiers conservative.
- [ ] Instanced particle rendering on GPU
  - **Goal:** 1M instanced spheres (impostor quads) at 60 fps on M-series.
  - **Design:** storage-buffer instance arrays; CPU `instancing.rs` becomes the parity reference.
  - **Files:** instanced pipeline in `src/wgpu_renderer.rs`; **Tests:** golden scene 4 + instance-count stress run (fps assertion env-gated).
- [ ] Screen-space fluid rendering for SPH (dep: instanced particles)
  - **Goal:** 250k-particle fluid surface at 30 fps 1080p.
  - **Design:** depth impostor → bilateral blur → normal reconstruction + thickness; golden-image tested; `metaball.rs` CPU path as visual reference.
  - **Files:** new screen-space fluid pass over the instanced pipeline; **Tests:** new golden scene (fluid surface) once stable.
  - **Risk:** bilateral-blur parameter tuning makes goldens brittle — lock parameters inside the scene definition before committing baselines.
- [ ] GPU raymarched volume rendering
  - **Goal:** 256³ volume at 60 fps with the existing `transfer_functions` module; CPU `volume_renderer.rs` kept as parity reference.
  - **Files:** raymarch WGSL + 3D-texture upload in `src/wgpu_renderer.rs`; **Tests:** golden scene 6 in GPU mode.
  - **Risk:** transfer-function sampling differences CPU vs GPU — share one LUT-generation code path between both backends.
- [ ] GPU marching cubes (dep: oxiphysics-gpu reduction lib)
  - **Goal:** 128³ SDF → mesh <10 ms; triangle count within 1% of CPU `isosurface.rs` (full 256-config CPU MC exists).
  - **Design:** scan/compact from oxiphysics-gpu reduction lib.
  - **Files:** GPU MC pass consuming the oxiphysics-gpu `kernels_wgsl` primitives; `isosurface.rs` stays the parity reference.
  - **Tests:** triangle-count parity vs `isosurface.rs` on seeded SDFs; timing assertion env-gated.

### Fields, camera, interaction
- [ ] Vector/tensor field viz consolidation
  - **Goal:** 3D LIC (2D LIC exists at `flow_visualization.rs:459`), oriented tensor glyphs, stream ribbons/tubes exportable to glTF; 3D LIC on 64³ field <2 s CPU.
  - **Files:** extend `flow_visualization.rs` (3D LIC), tensor glyph module, glTF export hook via oxiphysics-io.
  - **Tests:** golden scene (LIC slice) + CPU timing budget (not env-gated — CPU-only item).
- [ ] Unified camera + GPU picking
  - **Goal:** orbit/fly controller; pick in a 1M-instance scene <1 ms via ID buffer (CPU `picking.rs` Möller–Trumbore exists as reference).
  - **Files:** camera controller module (NEW), ID-buffer pass in `src/wgpu_renderer.rs`; **Tests:** pick-parity vs `picking.rs` on a shared scene.

## v0.3.0 — Presentation & remote

Exit criteria: OIT and both export paths golden-locked; remote viewer streams a live scene to a browser at 60 fps on localhost; the golden matrix grows to roughly 16 scenes (adding fluid, LIC, OIT, SVG goldens from the items above and below).

- [ ] Order-independent transparency (weighted-blended OIT)
  - **Goal:** 10k overlapping translucent particles render order-stable, locked by a golden image.
  - **Files:** OIT accumulation/resolve passes in `src/wgpu_renderer.rs`; new golden scene.
- [ ] Websocket remote viewer streaming to browser
  - **Goal:** 60 fps streaming of a 10k-body scene over localhost to a browser page that reuses the wasm `webgl_bridge` renderer.
  - **Design:** feature-gated `remote-viewer`; pure-Rust WebSocket (no native TLS deps in default).
  - **Files:** `remote-viewer` server module (NEW, feature-gated); browser page lives with the wasm crate's `web-demos/`.
  - **Tests:** localhost stream smoke test (frames received and decoded), skip-not-fail without a browser; fps measured by the local script, not asserted by default.
  - **Risk:** frame-encoding cost dominating at 60 fps — start with raw LDR frames on localhost; add compression only if the budget demands it (oxiarc-* if so).
- [ ] Animation/sequence export
  - **Goal:** 300-frame turntable or sim playback to PNG sequence + glTF animation in one API call (drives docs and the benchmark dashboard — root Phase 28).
  - **Files:** sequence-export API on the renderer facade; glTF animation written through oxiphysics-io.
- [ ] SVG/PDF vector export for scientific plots
  - **Goal:** pixel-perfect SVG goldens for 6 chart types (chart/heatmap/statistical/uncertainty).
  - **Files:** SVG/PDF backends for the plotting modules.
  - **Tests:** SVG golden files alongside the PNG goldens in `tests/golden/`.

## v1.0 — Production & validation

Exit criteria: one renderer abstraction, semver-clean; OpenXR feature builds with its pure-Rust fallback tests green (HMD run env-gated).

- [ ] OpenXR VR backend, feature-gated non-default (FFI loader mirrors the cudarc pattern; default stays pure Rust)
  - **Goal:** stereo render of any GpuRenderer scene at 72 fps on a reference HMD; pure-Rust fallback = side-by-side stereo offscreen pair extending the existing `vr_visualization.rs` model.
  - **Note:** the pure-Rust stereo fallback and the feature plumbing can land hardware-free; the 72 fps acceptance run is hardware-gated (see Deferred).
  - **Files:** `openxr-backend` feature module (NEW, non-default); stereo fallback extends `vr_visualization.rs`.
- [ ] Renderer API freeze
  - **Goal:** one `Renderer` trait over CPU rasterizer / GpuRenderer / wasm canvas; cargo-semver-checks clean across 1.x.
  - **Design:** trait extraction is the last move — it freezes only after the v0.2.0/v0.3.0 GPU surface stops churning.
  - **Files:** renderer trait module unifying the rasterizer / `GpuRenderer` / `wasm_canvas` facades.

### Cross-crate dependencies (for sequencing)
- GPU marching cubes ← oxiphysics-gpu v0.2.0 reduction lib (scan/compact)
- Screen-space fluid ← instanced particles (this crate)
- Remote viewer ← oxiphysics-wasm `webgl_bridge`
- Animation export → feeds root Phase 28 benchmark dashboard and the mdBook docs
- OpenXR feature gate ← precedent: oxiphysics-gpu `cuda-backend` (cudarc pattern)

## Deferred / research track

- [~] OpenXR on-HMD acceptance run (72 fps verification) — Ready when: reference HMD hardware is available; until then the env-gated test stays skip-not-fail (same pattern as `OXIPHYSICS_RTX_BENCH`) and the side-by-side stereo offscreen pair is the testable target.
- [~] Neural-rendering GPU promotion (the shipped Phase-4 neural rendering module is CPU-only) — Ready when: GpuRenderer completion + instancing have landed and a concrete consumer (e.g. learned denoising for the raymarcher) is identified; research track until then.
- [~] In-browser WebGPU rendering path — Ready when: oxiphysics-wasm v0.3.0 "WebGPU compute in browser" lands; viz then reuses it through the remote-viewer / `webgl_bridge` integration rather than growing a second browser renderer.

---
Conventions: `[x]` done / `[ ]` planned / `[~]` deferred until the stated condition. Dates are YYYY-MM-DD. **Goal:**/**Design:** wording is the per-item contract; acceptance numbers are testable budgets (env-gated where hardware-dependent), not aspirations.
