# oxiphysics-viz TODO

Last updated: 2026-06-06 | Version: 0.1.2

## Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Implement primary algorithms (rasterizer, post-processing, mesh gen)
- [x] Add integration tests (4,114 tests passing)
- [x] Performance benchmarks

## Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization

## Phase 4: Advanced Features
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

## Future / Post-0.1
- [x] GPU/wgpu backend (feature-gated) — `wgpu_renderer` (`GpuRenderer`, Phong+particle+blit WGSL shaders, `GpuRenderTarget`, CPU fallback)
- [x] WebAssembly canvas rendering path (`wasm_canvas` — `CanvasBuffer`, `CanvasRasterizer`, Porter-Duff blending)
- [x] Interactive picking / selection in rasterizer (`picking` — `Picker`, Möller–Trumbore, `SelectionSet`)
- [x] HDR framebuffer and wider color spaces (`hdr_framebuffer` — f32 RGBA, ACES/Reinhard/Filmic/AgX tone mapping)
