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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 246 (critical 2, high 49, medium 138, low 57); 89 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (2)
- [ ] **F5176** (critical, fabrication, easy/small) `scientific_viz.rs:466` — bootstrap_ci circular 'resampling' always yields zero-width CI
- [ ] **F5204** (critical, fabrication, easy/small) `stress_viz/functions.rs:376` — principal_eigenvectors discards eigen-solve, always returns identity axes

### High (49)
- [ ] **F4971** (high, bug, easy/trivial) `chart/types.rs:461` — ContourPlot::marching_squares cases 5/10 and 6/9 swapped
- [ ] **F4976** (high, bug, easy/trivial) `cross_section_viz.rs:434` — extract_contour marching-squares cases 5/10 and 6/9 swapped
- [ ] **F4990** (high, panic, easy/trivial) `data_viz.rs:403` — Histogram Scott/FD bin count unbounded; constant data aborts allocation
- [ ] **F4994** (high, bug, easy/small) `flow_visualization.rs:345` — VectorField central differences divide by 4·dx; derivatives halved
- [ ] **F5000** (high, bug, easy/trivial) `flow_viz.rs:537` — middle_eigenvalue_3x3_sym reflects spectrum; lambda2_criterion values wrong
- [ ] **F5003** (high, bug, medium/small) `flow_viz.rs:885` — classify_from_gradient applies 2D trace/det rules to 3x3 Jacobian
- [ ] **F5006** (high, fabrication, easy/small) `fluid_visualization.rs:764` — marching_cubes uses partial triangle table; most cases dropped
- [ ] **F5011** (high, fabrication, medium/small) `fluid_viz/functions.rs:120` — compute_tke_spectrum bins samples by index; no Fourier transform
- [ ] **F5012** (high, fabrication, medium/small) `fluid_viz/functions.rs:210` — level_set_mean_curvature returns |∇φ| instead of curvature
- [ ] **F5013** (high, fabrication, medium/small) `fluid_viz/functions.rs:224` — level_set_mean_curvature returns |∇φ| instead of curvature
- [ ] **F5014** (high, fabrication, medium/small) `fluid_viz/types.rs:85` — FreeSurfaceExtractor emits fixed +X slivers from x-edges only
- [ ] **F5026** (high, bug, easy/trivial) `gpu_viz.rs:1014` — FrustumCull::from_perspective left/right plane normals point outward
- [ ] **F5039** (high, panic, easy/trivial) `isosurface.rs:368` — TRI_TABLE row 203 lacks -1 terminator; marching_cubes panics
- [ ] **F5049** (high, bug, easy/trivial) `medical_viz/types_core.rs:229` — DoseDistribution::dose_at_volume (D_V) quantile index inverted
- [ ] **F5052** (high, perf, easy/small) `medical_viz/types_ext.rs:142` — VolumeRenderer::cast_ray rescans entire volume every sample; march unclipped
- [ ] **F5053** (high, bug, easy/trivial) `medical_viz/types_ext.rs:164` — Volume ray caster composites far samples over near ones
- [ ] **F5054** (high, panic, easy/trivial) `medical_viz/types_ext.rs:448` — VolumeData::get/sample index past data when iz >= nz
- [ ] **F5055** (high, fabrication, easy/small) `medical_viz/types_ext.rs:1014` — IsosurfaceRenderer::marching_cube_triangles groups crossed-edge vertices in triples, no case table
- [ ] **F5062** (high, fabrication, easy/small) `metaball.rs:415` — Metaball marching cubes chunks crossed edges instead of using tri-table
- [ ] **F5066** (high, bug, easy/small) `molecular_viz/types.rs:237` — RDF ideal count uses N_A·N_B·ρ, so g(r)→1/N_B
- [x] **F5083** (high, perf, medium/medium) `network_visualization.rs:1001` — greedy_modularity recomputes full O(n²·deg) modularity per candidate move
- [x] **F5084** (high, bug, medium/small) `network_visualization.rs:1130` — spectral_clustering power-iterates L and gets the largest eigenvector, not the Fiedler vector
- [x] **F5085** (high, bug, medium/small) `network_visualization.rs:1139` — approximate_fiedler power-iterates L, converging to largest eigenvector
- [x] **F5094** (high, fabrication, easy/small) `neural_rendering.rs:764` — Gaussian3D::eval_density discards covariance, ignores rotation quaternion
- [x] **F5097** (high, bug, easy/small) `particle_effects.rs:124` — Fire/Smoke emitters emit zero particles at real frame rates
- [x] **F5099** (high, panic, easy/small) `particle_effects.rs:245` — Empty-range random_range panics in emitters and seeding_uniform
- [x] **F5104** (high, fabrication, medium/small) `particle_renderer/types.rs:488` — MetaballField::extract_isosurface emits fixed dummy triangle, not marching cubes
- [x] **F5129** (high, bug, easy/small) `picking.rs:62` — Ray::from_screen forces view z=-1, wrong direction unless near=1
- [x] **F5154** (high, bug, easy/small) `real_time_viz.rs:417` — HiZBuffer::is_occluded tests region MIN depth, culling visible objects
- [x] **F5161** (high, bug, medium/small) `scene/types.rs:589` — Shadow-map frustum omits light-view rotation; depth maps outside clip
- [ ] **F5175** (high, bug, easy/small) `scientific_viz.rs:214` — eigen3_sym flattens eigenvector rows (transposed) and leaves eigenvalues unsorted
- [ ] **F5178** (high, bug, easy/trivial) `scientific_viz.rs:1021` — classify_fixed_points halves a12 and uses one-sided a22
- [ ] **F5185** (high, bug, easy/small) `slice_renderer.rs:248` — scalar_to_colormap_slice ignores nu/nv, reshapes to sqrt square
- [ ] **F5189** (high, bug, easy/small) `streamlines/functions.rs:89` — lambda2_eigenvalues Cardano formula wrong; λ2 vortex detection broken
- [ ] **F5194** (high, bug, easy/small) `streamlines/functions.rs:339` — Flow topology det(J) omits off-diagonal term; cell_size<=0 hangs
- [ ] **F5202** (high, bug, easy/small) `stress_viz/functions.rs:173` — approximate_eigenvectors returns arbitrary perpendiculars, misordered vs magnitudes
- [ ] **F5207** (high, bug, easy/small) `structural_viz.rs:655` — eigenvector_sym3 returns [1,0,0] for every repeated eigenvalue
- [ ] **F5221** (high, fabrication, hard/large) `topology_viz.rs:654` — ReebGraph::from_triangle_mesh chains critical points by scalar value
- [ ] **F5222** (high, bug, medium/medium) `topology_viz.rs:882` — Ascending manifolds treat minima and critical-point ids as maximum vertices
- [ ] **F5223** (high, fabrication, hard/medium) `topology_viz.rs:1161` — ContourTree arcs fabricated as value-sorted chain; join/split trees discarded
- [ ] **F5224** (high, fabrication, medium/small) `topology_viz_ext/types.rs:180` — Jacobi set hard-codes cross partials to zero
- [ ] **F5230** (high, bug, easy/trivial) `uncertainty_viz.rs:198` — Bootstrap uses LCG low bits; zero-width CI for n=2^k
- [ ] **F5234** (high, panic, easy/trivial) `volume_renderer.rs:196` — TransferFunction eval panics on NaN scalar (usize underflow)
- [ ] **F5236** (high, bug, medium/small) `volume_renderer.rs:447` — IsosurfaceExtractor fan-triangulates crossings in edge-index order
- [ ] **F5237** (high, fabrication, easy/small) `volume_rendering/functions.rs:151` — extract_isosurface_simple emits fixed fake triangle per crossed cell
- [ ] **F5255** (high, bug, medium/small) `volumetric_rendering.rs:357` — iso_surface_value: world-origin atan2 sort and stale fan pivot
- [ ] **F5258** (high, bug, easy/small) `vr_visualization.rs:122` — VrCamera::view_matrix mixes R and R^T; culling wrong when rotated
- [ ] **F5266** (high, panic, easy/trivial) `wasm_canvas.rs:240` — draw_filled_circle wraps negative span to usize and panics
- [ ] **F5274** (high, bug, easy/small) `wgpu_renderer.rs:612` — GpuRenderer::draw_particles is a silent no-op (no CPU path)

### Medium (138)
- [ ] **F4954** (medium, api, easy/small) `crates/oxiphysics-viz/Cargo.toml:17` — wgpu forced on every dependent; viz README says no wgpu dependency
- [ ] **F4955** (medium, missing-feature, medium/medium) `advanced_rendering.rs:594` — RaytracingBvh claims SAH build and closest-hit; neither exists
- [ ] **F4957** (medium, fabrication, easy/small) `advanced_rendering.rs:764` — SH9Probe::evaluate silently drops all band-2 coefficients
- [ ] **F4961** (medium, bug, easy/trivial) `animation.rs:222` — Catmull-Rom interior tangents ignore time spacing
- [ ] **F4962** (medium, bug, medium/medium) `animation_system/functions.rs:283` — Pose blending lerps 4x4 TRS matrices element-wise; Additive not additive
- [ ] **F4963** (medium, perf, easy/trivial) `animation_system/types.rs:70` — Animator::update deep-clones every clip map each frame
- [ ] **F4964** (medium, fabrication, medium/medium) `animation_system/types.rs:142` — Blend tree lerps Mat4 element-wise; Additive is a cascaded lerp
- [ ] **F4965** (medium, bug, easy/small) `animation_system/types.rs:855` — State machine ignores blend_duration; automatic transitions never fire
- [ ] **F4967** (medium, panic, easy/small) `annotation_viz/functions.rs:275` — wrap_text splits at byte index; panics on multibyte text
- [ ] **F4968** (medium, bug, easy/trivial) `camera/functions.rs:55` — Catmull-Rom camera path extrapolates in the first segment
- [ ] **F4969** (medium, bug, easy/trivial) `camera/types.rs:134` — MultiTargetCamera smoothing inverted; smoothing=0 never moves
- [ ] **F4972** (medium, bug, easy/small) `colormap/functions.rs:50` — Viridis 0.25 stop wrong; data_viz/chart_3d viridis polynomials far off
- [ ] **F4974** (medium, bug, easy/small) `colormap/functions.rs:628` — CIEDE2000/Lab metrics feed gamma-encoded sRGB into linear-RGB rgb_to_lab
- [ ] **F4975** (medium, fabrication, easy/trivial) `colormap/functions.rs:858` — sample_turbo uses invented coefficients; upper half renders near-white
- [ ] **F4977** (medium, perf, easy/small) `cross_section_viz.rs:457` — Contour chaining is O(S²) and extends only forward, fragmenting lines
- [ ] **F4980** (medium, fabrication, easy/small) `cross_section_viz.rs:977` — cross_section_intersection uses bounding box; doc claims convex hull
- [ ] **F4986** (medium, pure-rust, easy/small) `data_pipeline.rs:187` — DataSeries::fft_power_spectrum is an O(n²) DFT, not an FFT
- [ ] **F4989** (medium, fabrication, easy/small) `data_viz.rs:119` — data_viz viridis polynomial is pink at mid-range
- [ ] **F4991** (medium, panic, easy/trivial) `data_viz.rs:628` — Heatmap Percentile range panics on empty data or p>100
- [ ] **F4995** (medium, fabrication, medium/small) `flow_visualization.rs:854` — FlowStatistics::lambda2_criterion returns −Q, not the λ2 eigenvalue
- [ ] **F4997** (medium, pure-rust, easy/small) `flow_visualization.rs:1159` — energy_spectrum_1d is an O(n²) DFT; oxifft unused
- [ ] **F4998** (medium, perf, easy/small) `flow_visualization.rs:1206` — energy_spectrum_2d is a brute-force O((nx·ny)²) DFT and panics on short input
- [ ] **F4999** (medium, bug, easy/small) `flow_viz.rs:407` — Clamped central differences halve boundary derivatives
- [ ] **F5001** (medium, panic, easy/trivial) `flow_viz.rs:558` — lic_2d asserts on input lengths and uses unseeded thread RNG
- [ ] **F5002** (medium, panic, easy/small) `flow_viz.rs:645` — seeding_uniform panics on degenerate bounds and is non-reproducible
- [ ] **F5004** (medium, fabrication, medium/small) `flow_viz.rs:934` — trace_separatrices seeds along the gradient's first row, not eigenvectors
- [ ] **F5009** (medium, bug, easy/small) `fluid_visualization.rs:1252` — sph_density_grid sums rho*W unnormalised; not a density field
- [ ] **F5017** (medium, bug, easy/small) `fluid_viz/types.rs:501` — LicRenderer backward pass retraces forward path; vy index unchecked
- [ ] **F5018** (medium, fabrication, easy/small) `fluid_viz/types.rs:679` — VolumeSlice::extract claims bilinear but copies node indices
- [ ] **F5019** (medium, bug, easy/small) `fluid_viz/types.rs:1001` — extract_crossings scans only x-edges, missing y/z crossings
- [ ] **F5023** (medium, bug, easy/small) `glyph_renderer.rs:376` — streamribbon twist rotates width into tangent, collapsing ribbon
- [ ] **F5027** (medium, bug, easy/small) `graph_viz.rs:202` — Out-of-range edges clamped to last node; cycle handling misdocumented
- [ ] **F5028** (medium, bug, easy/small) `graph_viz.rs:330` — Out-of-range edge endpoints silently clamped onto last node
- [ ] **F5035** (medium, panic, easy/small) `heatmap.rs:988` — HeatmapAnimation::sample_at panics on empty or mismatched animation
- [ ] **F5036** (medium, bug, easy/trivial) `instancing.rs:339` — Billboard basis is left-handed, so quads render mirrored
- [ ] **F5037** (medium, bug, easy/small) `interactive_viz.rs:429` — select_at picking ignores object rotation and scale
- [ ] **F5041** (medium, fabrication, medium/medium) `isosurface.rs:579` — dual_contouring_2d is a 1-D zero-crossing finder
- [ ] **F5044** (medium, api, easy/small) `lib.rs:61` — Root glob re-export leaks volume_rendering Camera/Vec3; Ray is shadowed
- [ ] **F5048** (medium, bug, easy/trivial) `medical_viz/types_core.rs:150` — DoseDistribution::sample_at extrapolates at negative coordinates
- [ ] **F5050** (medium, fabrication, easy/small) `medical_viz/types_core.rs:467` — SUV BSA hardcodes 170 cm patient height
- [ ] **F5057** (medium, bug, easy/small) `mesh_gen/functions.rs:79` — box_mesh ±Y faces and UV spheres wound inward
- [ ] **F5059** (medium, panic, medium/small) `mesh_render.rs:758` — Unclipped Bresenham hangs or overflows on off-screen projected vertices
- [ ] **F5060** (medium, perf, easy/small) `mesh_render.rs:986` — generate_lod is O(n²) greedy merge mislabelled edge collapse; can panic
- [ ] **F5061** (medium, bug, easy/small) `mesh_render.rs:1198` — MeshSkinData::skin_position ignores inverse bind matrices
- [ ] **F5063** (medium, perf, easy/small) `molecular_viz/functions.rs:93` — detect_bonds is O(n²) all-pairs with sqrt per pair
- [ ] **F5065** (medium, bug, easy/trivial) `molecular_viz/functions.rs:257` — marching_squares saddle case 5 emits intersecting segments
- [ ] **F5067** (medium, bug, medium/small) `molecular_viz/types.rs:1525` — TrajectoryFrame::rmsd returns 0.0 (perfect match) on length mismatch
- [ ] **F5068** (medium, bug, easy/trivial) `multiphysics_viz/functions.rs:158` — Interface normal is a rotated tangent; documented sort missing
- [ ] **F5069** (medium, perf, medium/small) `multiphysics_viz/functions.rs:345` — Field interpolation allocates and fully sorts all sources per target
- [ ] **F5070** (medium, fabrication, medium/medium) `multiphysics_viz/functions.rs:367` — 'Barycentric' interpolation is actually 3-nearest IDW
- [ ] **F5071** (medium, panic, easy/small) `multiphysics_viz/types.rs:343` — sync_grid OOM-aborts on duplicate or non-monotonic timestamps
- [ ] **F5072** (medium, panic, easy/trivial) `network_visualization.rs:218` — Recursive dfs_rec overflows stack on long paths
- [ ] **F5074** (medium, panic, easy/small) `network_visualization.rs:378` — Edges to unknown nodes panic in layout and metrics
- [ ] **F5075** (medium, fabrication, medium/medium) `network_visualization.rs:447` — barnes_hut_repulsion drops near multi-node cells; unused by apply()
- [ ] **F5076** (medium, bug, medium/small) `network_visualization.rs:589` — Sugiyama crossing minimization no-op for directed graphs, nondeterministic order
- [ ] **F5077** (medium, bug, medium/small) `network_visualization.rs:605` — Sugiyama barycenter uses successors; crossing minimisation is a no-op on DAGs
- [ ] **F5080** (medium, bug, easy/trivial) `network_visualization.rs:857` — density() halves undirected density (triangle K3 gives 0.5)
- [ ] **F5081** (medium, bug, medium/medium) `network_visualization.rs:961` — modularity mixes weighted A and m with unweighted degree; greedy is O(n^3)+
- [ ] **F5082** (medium, perf, medium/medium) `network_visualization.rs:991` — Greedy modularity is O(n³) per pass and uses unweighted degree
- [ ] **F5086** (medium, fabrication, medium/medium) `network_visualization.rs:1184` — spectral_clustering: fake bisection and power iteration finds λ_max
- [ ] **F5087** (medium, bug, easy/small) `network_visualization.rs:1300` — NetworkAnimation::frame_at replays events in insertion order
- [ ] **F5090** (medium, fabrication, medium/medium) `neural_rendering.rs:385` — NeRfVolume ignores num_fine_samples and white_background flags
- [ ] **F5091** (medium, fabrication, medium/small) `neural_rendering.rs:422` — NeRF render_ray ignores num_fine_samples and white_background
- [ ] **F5095** (medium, fabrication, medium/small) `neural_rendering.rs:1182` — NeuralDenoiser ignores kernel type and sigmas: always plain Gaussian
- [ ] **F5096** (medium, fabrication, medium/medium) `neural_rendering.rs:1373` — LatentSpaceViz::project_pca is not PCA; t-SNE params unused
- [ ] **F5098** (medium, panic, easy/small) `particle_effects.rs:191` — Smoke/explosion emitters panic on empty random ranges
- [ ] **F5102** (medium, panic, easy/trivial) `particle_renderer/types.rs:62` — ParticleTrail::add_point panics when max_points == 0
- [ ] **F5103** (medium, fabrication, medium/small) `particle_renderer/types.rs:393` — SplatRenderer::splat_to_buffer ignores position, radius and cutoff
- [ ] **F5105** (medium, bug, easy/trivial) `particle_renderer/types.rs:590` — ColorMap::Viridis has wrong RGB stops and decaying alpha
- [ ] **F5106** (medium, perf, medium/small) `particle_renderer/types.rs:678` — SphDensityField::slice_z0 is O(grid x particles) despite compact kernel
- [ ] **F5110** (medium, bug, medium/small) `path_viz.rs:366` — find_path_intersections compares midpoints only, misses true crossings
- [ ] **F5116** (medium, bug, easy/trivial) `performance_viz.rs:1384` — ascii_convergence_graph renders blank for residuals below 1
- [ ] **F5118** (medium, bug, medium/small) `phase_field_viz.rs:521` — DendriteTips::detect marks every flat-front cell as a tip
- [ ] **F5120** (medium, bug, easy/trivial) `physics_animation.rs:130` — Over-damped Spring easing returns constant 1.0 (coefficient sign error)
- [ ] **F5121** (medium, bug, easy/trivial) `physics_animation.rs:331` — AnimationBlender::crossfade always samples both clips at t=0
- [ ] **F5122** (medium, bug, easy/trivial) `physics_dashboard/types.rs:66` — Momentum drift compares magnitudes only; rotated momentum reads conserved
- [ ] **F5124** (medium, bug, easy/trivial) `physics_dashboard/types.rs:320` — Energy anomaly detection disabled for negative total energy
- [ ] **F5128** (medium, panic, easy/trivial) `physics_dashboard/types_core.rs:1267` — VelocityHistogram with zero bins panics on add
- [ ] **F5130** (medium, bug, easy/trivial) `picking.rs:170` — Picker::default() has max_distance 0, so it never reports a hit
- [ ] **F5132** (medium, panic, easy/trivial) `picking.rs:238` — pick_object panics on out-of-range mesh index
- [ ] **F5136** (medium, bug, easy/trivial) `post_processing/effects.rs:22` — GaussianBlur with sigma=0 produces an all-NaN image
- [ ] **F5140** (medium, bug, easy/trivial) `post_processing/image_filters.rs:67` — FXAA blends along the edge instead of across it (no-op)
- [ ] **F5145** (medium, bug, easy/trivial) `post_processing/processor.rs:127` — compute_lens_flare returns black image when flare disabled (default)
- [ ] **F5147** (medium, fabrication, easy/small) `postprocess.rs:715` — compute_ssao_sample ignores fragment xy and sample orientation
- [ ] **F5150** (medium, bug, easy/small) `primitives.rs:194` — Mesh primitives have inconsistent triangle winding
- [ ] **F5151** (medium, bug, easy/small) `primitives.rs:1055` — DrawCallBatch::merge_consecutive changes drawn geometry
- [ ] **F5160** (medium, panic, easy/small) `scene/types.rs:497` — Scene::deserialize trusts parent ids: index panic and stack-overflow abort
- [ ] **F5162** (medium, bug, medium/small) `scene/types.rs:613` — cull_invisible_objects uses local position and wrong occluder extents
- [ ] **F5163** (medium, missing-feature, medium/medium) `scientific_plot.rs:552` — AxisScale::Log and tick_spacing silently ignored by SVG renderer
- [ ] **F5165** (medium, panic, easy/trivial) `scientific_plots.rs:630` — correlation_dimension divides by zero when sample == 0
- [ ] **F5166** (medium, perf, easy/small) `scientific_plots.rs:735` — STFT uses naive O(N^2) DFT with per-sample trig instead of oxifft
- [ ] **F5168** (medium, perf, easy/small) `scientific_plots.rs:1082` — STFT uses naive O(N²) DFT per frame instead of oxifft rfft
- [ ] **F5174** (medium, fabrication, easy/trivial) `scientific_viz.rs:122` — ellipsoid_surface_points ignores eigenvectors; glyphs always axis-aligned
- [ ] **F5179** (medium, bug, easy/trivial) `shader/types.rs:121` — Shadow view matrix degenerate for vertical light direction
- [ ] **F5181** (medium, fabrication, medium/small) `shader/types.rs:889` — dipole_reflectance is not the dipole model; singular at r→0
- [ ] **F5182** (medium, bug, easy/trivial) `shader/types.rs:973` — PbrMaterial::ggx_distribution returns NaN for roughness 0 at n.h=1
- [ ] **F5183** (medium, stub, medium/medium) `shader/types.rs:1010` — PbrMaterial::to_material emits empty stub GLSL; Phong GLSL uncompilable
- [ ] **F5187** (medium, panic, easy/medium) `statistical_viz.rs:29` — NaN-inconsistent sort comparator; convex_hull_2d expect panics on NaN
- [ ] **F5190** (medium, bug, easy/trivial) `streamlines/functions.rs:207` — LIC backward pass integrates forward: step.abs() drops sign
- [ ] **F5193** (medium, fabrication, medium/small) `streamlines/functions.rs:291` — Enhanced LIC 'passes' only reapply a sigmoid; backward half broken
- [ ] **F5196** (medium, panic, easy/trivial) `streamlines/functions.rs:561` — arrow_glyphs_along_streamline infinite loop/OOM for interval <= 0
- [ ] **F5197** (medium, perf, medium/medium) `streamlines/functions.rs:619` — Jobard-Lefer distance test brute-forces all points of all streamlines
- [ ] **F5198** (medium, bug, medium/medium) `streamlines/functions.rs:652` — Jobard–Lefer: brute-force distance, whole-line rejection, hard cap 200
- [ ] **F5199** (medium, bug, easy/trivial) `streamlines/types.rs:342` — Stream ribbon twists at full vorticity, missing factor 1/2
- [ ] **F5200** (medium, perf, easy/small) `streamlines/types.rs:471` — Streakline::compute_steady re-integrates every particle (O(K^2)); hangs for dt=0
- [ ] **F5201** (medium, bug, easy/trivial) `stress_viz/functions.rs:108` — StressInvariant::J3 colours by I3=det(sigma), not deviatoric J3
- [ ] **F5205** (medium, bug, easy/trivial) `stress_viz/types.rs:562` — von_mises_strain off by sqrt(3/2); SQRT_2/SQRT_2 no-op
- [ ] **F5206** (medium, bug, easy/small) `structural_viz.rs:392` — FEM wireframe uses cyclic node order for solid elements
- [ ] **F5208** (medium, bug, easy/small) `structural_viz.rs:679` — Repeated eigenvalues give coincident [1,0,0] eigenvectors; scale dropped
- [ ] **F5212** (medium, bug, medium/small) `tensor_visualization/types.rs:418` — Poincaré index treats eigenvector line field as vector field
- [ ] **F5216** (medium, panic, easy/small) `terrain_viz.rs:127` — Heightfield normal/slope/aspect panic for grids narrower than 3
- [ ] **F5217** (medium, bug, easy/trivial) `terrain_viz.rs:758` — Aspect computed as upslope angle from +X, not downslope from North
- [ ] **F5219** (medium, stub, easy/small) `text_renderer.rs:120` — Bitmap font is placeholder: all letters render as hollow boxes
- [ ] **F5220** (medium, bug, hard/medium) `topology_viz.rs:151` — Bottleneck/Wasserstein use greedy, cross-dimension, asymmetric matching
- [ ] **F5225** (medium, fabrication, hard/medium) `topology_viz_ext/types.rs:1058` — bottleneck_distance_approx zip-truncates and compares persistence only
- [ ] **F5227** (medium, bug, medium/small) `transfer_functions/types.rs:323` — Pre-integration integrates over scalar delta, discontinuous and unordered
- [ ] **F5229** (medium, bug, easy/small) `uncertainty_viz.rs:107` — quantile uses rounded nearest rank; even-n median wrong
- [ ] **F5235** (medium, bug, easy/small) `volume_renderer.rs:270` — RaycastVolume::render_ray marches max_steps without AABB clipping
- [ ] **F5240** (medium, bug, easy/trivial) `volume_rendering/functions.rs:214` — iso_surface_normal returns un-negated gradient, contradicting its doc
- [ ] **F5241** (medium, bug, easy/small) `volume_rendering/functions.rs:799` — Volume ray marchers loop forever when step_size <= 0
- [ ] **F5243** (medium, bug, medium/small) `volume_rendering/types.rs:110` — AdaptiveRayMarcher lacks opacity correction; skipping marcher overshoots bricks
- [ ] **F5244** (medium, bug, easy/trivial) `volume_rendering/types.rs:228` — OccupancyGrid bricks lack 1-voxel apron; skipping drops interpolated density
- [ ] **F5248** (medium, bug, easy/small) `volume_rendering/types.rs:490` — Ray marchers loop forever for step_size <= 0
- [ ] **F5249** (medium, bug, medium/small) `volume_rendering/types.rs:545` — Empty-space skipping: fixed coarse step overshoots bricks; double lookup
- [ ] **F5252** (medium, bug, easy/small) `volume_rendering/types.rs:1081` — Trilinear samplers off by half a voxel vs cell-centred data
- [ ] **F5254** (medium, bug, easy/small) `volumetric_rendering.rs:185` — ray_march yields nothing for origins outside volume; no early-out
- [ ] **F5262** (medium, bug, easy/small) `vr_visualization.rs:484` — doppler_shift wrong listener sign and always uses source 0
- [ ] **F5263** (medium, bug, easy/trivial) `vr_visualization.rs:501` — attenuation applies 1/d² to amplitude, contradicting its doc
- [ ] **F5267** (medium, panic, easy/small) `wasm_canvas.rs:309` — draw_rect_outline: unclamped x0/y0 corrupt rows/panic; width-1 underflow
- [ ] **F5269** (medium, panic, easy/small) `wgpu_renderer.rs:362` — Mutating pub GpuRenderer.config desyncs buffers and panics
- [ ] **F5270** (medium, fabrication, hard/medium) `wgpu_renderer.rs:388` — GpuRenderer GPU path is unreachable; module docs claim GPU rendering
- [ ] **F5273** (medium, perf, easy/trivial) `wgpu_renderer.rs:584` — CPU draw_mesh clones all mesh buffers on every draw call
- [ ] **F5275** (medium, bug, medium/medium) `wgpu_renderer.rs:893` — CPU rasteriser has no near-plane clipping; w<=0 vertices mis-project
- [ ] **F5276** (medium, bug, easy/small) `wgpu_renderer.rs:919` — CPU rasteriser light sign inverted and light rotated with model
- [ ] **F5278** (medium, panic, easy/small) `wireframe/mesh_ops.rs:145` — Silhouette/feature/dihedral edge functions panic on out-of-range indices
- [ ] **F5280** (medium, fabrication, easy/small) `wireframe/mesh_ops.rs:453` — AA wireframe opacity and depth_test params silently dropped
- [ ] **F5282** (medium, bug, easy/small) `wireframe/primitives.rs:274` — from_capsule/draw_capsule ignore p0→p1 axis and omit hemispherical caps
- [ ] **F5284** (medium, perf, medium/medium) `wireframe/renderer.rs:413` — classify_hidden_lines is O(edges × triangles) brute-force ray casting

### Low (57) — tracked in the external ledger
F4953, F4958, F4960, F4970, F4978, F4979, F4981, F4982, F4985, F4987, F4988, F4993, F5005, F5010, F5038, F5042, F5043, F5045, F5046, F5073, F5079, F5093, F5108, F5112, F5115, F5119, F5126, F5131, F5135, F5137, F5144, F5146, F5153, F5155, F5159, F5167, F5169, F5170, F5171, F5177, F5180, F5184, F5191, F5195, F5211, F5228, F5245, F5246, F5250, F5253, F5259, F5260, F5261, F5268, F5271, F5277, F5285

<!-- AUDIT-2026-10-12 END -->
