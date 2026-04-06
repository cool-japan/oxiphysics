# oxiphysics-viz

**Status:** Alpha (CPU-only soft renderer) | **Version:** 0.1.0 | **Tests:** 4,072

Largest visualization crate in the OxiPhysics ecosystem — a pure-CPU software renderer with
scientific, engineering, and physics-specific visualization tools.
Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine. No GPU/wgpu dependencies.

## Domains

- **Rendering Core:** Software rasterizer, Phong shading, wireframe, framebuffer management
- **Post-Processing:** Bloom, depth-of-field, Gaussian blur, SSAO, tone mapping, vignette
- **Mesh Generation:** Arrow, box, plane, sphere primitives
- **Scene & Camera:** Scene graph, camera control, gizmos, debug overlay
- **Colormaps:** Scalar-to-color mapping, transfer functions
- **Streamlines:** 3D streamline tracing and rendering
- **Stress / FEM Viz:** Principal stress glyphs, von Mises colour maps, structural visualization
- **Volume Rendering:** Isosurface, slice renderer, cross-section, volumetric rendering
- **Particle Systems:** Particle renderer, particle effects, particle trails, instancing
- **Physics Animation:** Physics-driven animation, real-time viz, performance viz
- **Scientific Plotting:** Charts, heatmaps, statistical viz, uncertainty viz, phase-field viz
- **Medical / Molecular:** Medical imaging viz, molecular visualization
- **Fluid Viz:** Fluid visualization, LBM rendering
- **Advanced:** Neural rendering, metaball, procedural texture, topology viz, tensor visualization
- **Text & Fonts:** Font rendering, text renderer, annotation visualization
- **Network / Graph:** Graph viz, network viz, glyph renderer
- **VR:** VR-ready visualization pipeline

## Key Exports

```rust
use oxiphysics_viz::{
    // Colormaps
    Colormap, map_scalar,
    // Mesh generation
    // (arrow, box, plane, sphere builders)
    // Rendering
    Framebuffer, PhongShader, SoftwareRasterizer, WireframeRenderer,
    // Post-processing
    PostProcessPipeline, Bloom, DepthOfField, GaussianBlur,
    Ssao, ToneMapping, Vignette,
    // Streamlines
    trace_streamline,
    // Stress visualization
    principal_stress_glyphs, von_mises_colors,
    // Config
    VizConfig, Renderer,
};
```

## Modules (45+)

`advanced_rendering`, `animation_system`, `annotation_viz`, `camera`, `chart`, `colormap`,
`cross_section`, `data_viz`, `debug_overlay`, `fluid_viz`, `font_rendering`, `gizmos`,
`glyph_renderer`, `graph_viz`, `heatmap`, `instancing`, `isosurface`, `medical_viz`,
`mesh_gen`, `metaball`, `molecular_viz`, `network_viz`, `neural_rendering`, `particle_effects`,
`particle_renderer`, `particle_trails`, `performance_viz`, `phase_field_viz`, `physics_animation`,
`physics_dashboard`, `post_processing`, `procedural_texture`, `real_time_viz`, `scene`,
`scientific_plotting`, `shader`, `slice_renderer`, `statistical_viz`, `streamlines`,
`stress_viz`, `structural_viz`, `tensor_visualization`, `terrain_renderer`, `text_renderer`,
`topology_viz`, `transfer_functions`, `uncertainty_viz`, `volume_rendering`, `volumetric_rendering`,
`vr_visualization`, `wireframe`

## Statistics

| Metric | Count |
|--------|-------|
| Public items | 5,229 |
| Tests | 4,072 |
| Stubs | 0 |
| Modules | 45+ |

## License

MIT OR Apache-2.0 — Copyright © COOLJAPAN OU (Team Kitasan)
