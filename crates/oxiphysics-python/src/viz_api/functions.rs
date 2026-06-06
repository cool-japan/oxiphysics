//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use pyo3::prelude::*;

use super::types::{
    ColormapKind, PyCamera, PyColormap, PyDebugOverlay, PyLight, PyMaterial, PyParticleRenderer,
    PyPostProcessor, PySceneNode, PyStreamlineTracer, PyStressVisualizer, PyTransferFunction,
    PyVolumeRenderer, Streamline, StressTensor, StressVisOutput, ToneMapping, TransferPoint,
};
use super::types_3::{PyMeshRenderer, PySceneGraph};

/// Convert HSV color to RGB.
///
/// All inputs and outputs are in \[0, 1\] range.
#[pyfunction]
pub fn hsv_to_rgb(h: f64, s: f64, v: f64) -> [f64; 3] {
    if s == 0.0 {
        return [v, v, v];
    }
    let h6 = h * 6.0;
    let i = h6.floor() as i32;
    let f = h6 - h6.floor();
    let p = v * (1.0 - s);
    let q = v * (1.0 - s * f);
    let t = v * (1.0 - s * (1.0 - f));
    match i % 6 {
        0 => [v, t, p],
        1 => [q, v, p],
        2 => [p, v, t],
        3 => [p, q, v],
        4 => [t, p, v],
        _ => [v, p, q],
    }
}
/// Compute smooth normals from a triangle mesh.
///
/// `vertices` is a flat array of `[x, y, z]` values (length = 3 * num_vertices).
/// `indices` is a flat array of triangle vertex indices (length = 3 * num_triangles).
/// Returns a flat array of normals, same length as `vertices`.
#[pyfunction]
pub fn compute_normals_from_mesh(vertices: Vec<f64>, indices: Vec<u32>) -> Vec<f64> {
    let nv = vertices.len() / 3;
    let mut normals = vec![0.0f64; vertices.len()];
    let ntri = indices.len() / 3;
    for t in 0..ntri {
        let i0 = indices[3 * t] as usize;
        let i1 = indices[3 * t + 1] as usize;
        let i2 = indices[3 * t + 2] as usize;
        let p0 = [vertices[3 * i0], vertices[3 * i0 + 1], vertices[3 * i0 + 2]];
        let p1 = [vertices[3 * i1], vertices[3 * i1 + 1], vertices[3 * i1 + 2]];
        let p2 = [vertices[3 * i2], vertices[3 * i2 + 1], vertices[3 * i2 + 2]];
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        for idx in [i0, i1, i2] {
            normals[3 * idx] += n[0];
            normals[3 * idx + 1] += n[1];
            normals[3 * idx + 2] += n[2];
        }
    }
    for v in 0..nv {
        let nx = normals[3 * v];
        let ny = normals[3 * v + 1];
        let nz = normals[3 * v + 2];
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len > 1e-12 {
            normals[3 * v] /= len;
            normals[3 * v + 1] /= len;
            normals[3 * v + 2] /= len;
        }
    }
    normals
}
/// Generate a UV sphere mesh.
///
/// Returns `(vertices, normals, indices)` as flat `f64`/`u32` arrays.
#[pyfunction]
pub fn generate_sphere_mesh(
    radius: f64,
    lat_segments: u32,
    lon_segments: u32,
) -> (Vec<f64>, Vec<f64>, Vec<u32>) {
    let mut verts = Vec::new();
    let mut norms = Vec::new();
    let mut idxs = Vec::new();
    for lat in 0..=lat_segments {
        let theta = std::f64::consts::PI * lat as f64 / lat_segments as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        for lon in 0..=lon_segments {
            let phi = 2.0 * std::f64::consts::PI * lon as f64 / lon_segments as f64;
            let x = sin_t * phi.cos();
            let y = cos_t;
            let z = sin_t * phi.sin();
            verts.push(radius * x);
            verts.push(radius * y);
            verts.push(radius * z);
            norms.push(x);
            norms.push(y);
            norms.push(z);
        }
    }
    for lat in 0..lat_segments {
        for lon in 0..lon_segments {
            let row = lon_segments + 1;
            let a = lat * row + lon;
            let b = a + row;
            let c = b + 1;
            let d = a + 1;
            idxs.push(a);
            idxs.push(b);
            idxs.push(c);
            idxs.push(a);
            idxs.push(c);
            idxs.push(d);
        }
    }
    (verts, norms, idxs)
}
/// Generate a box (cuboid) mesh centered at the origin.
///
/// `half_extents` = \[hx, hy, hz\].
/// Returns `(vertices, normals, indices)`.
#[pyfunction]
pub fn generate_box_mesh(half_extents: [f64; 3]) -> (Vec<f64>, Vec<f64>, Vec<u32>) {
    let [hx, hy, hz] = half_extents;
    let face_data: &[([f64; 3], [f64; 3])] = &[
        ([hx, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([-hx, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        ([0.0, hy, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, -hy, 0.0], [0.0, -1.0, 0.0]),
        ([0.0, 0.0, hz], [0.0, 0.0, 1.0]),
        ([0.0, 0.0, -hz], [0.0, 0.0, -1.0]),
    ];
    let offsets: &[[f64; 2]] = &[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    let mut verts = Vec::new();
    let mut norms = Vec::new();
    let mut idxs = Vec::new();
    for (face_idx, (center, normal)) in face_data.iter().enumerate() {
        let base = (face_idx * 4) as u32;
        let (t1, t2) = if normal[0].abs() > 0.5 {
            ([0.0, hy, 0.0], [0.0, 0.0, hz])
        } else if normal[1].abs() > 0.5 {
            ([hx, 0.0, 0.0], [0.0, 0.0, hz])
        } else {
            ([hx, 0.0, 0.0], [0.0, hy, 0.0])
        };
        for off in offsets {
            let x = center[0] + off[0] * t1[0] + off[1] * t2[0];
            let y = center[1] + off[0] * t1[1] + off[1] * t2[1];
            let z = center[2] + off[0] * t1[2] + off[1] * t2[2];
            verts.push(x);
            verts.push(y);
            verts.push(z);
            norms.push(normal[0]);
            norms.push(normal[1]);
            norms.push(normal[2]);
        }
        idxs.push(base);
        idxs.push(base + 1);
        idxs.push(base + 2);
        idxs.push(base);
        idxs.push(base + 2);
        idxs.push(base + 3);
    }
    (verts, norms, idxs)
}
pub(super) fn viridis_sample(t: f64) -> [f64; 3] {
    let r = 0.2777 * t.powi(3) - 0.8673 * t.powi(2) + 0.4756 * t + 0.2665;
    let g = -0.0955 * t.powi(3) + 0.5925 * t.powi(2) + 0.2843 * t + 0.0038;
    let b = -1.0803 * t.powi(3) + 1.1215 * t.powi(2) - 0.6424 * t + 0.5305;
    [r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0)]
}
pub(super) fn plasma_sample(t: f64) -> [f64; 3] {
    let r = (0.9 * t + 0.05).clamp(0.0, 1.0);
    let g = (4.0 * t * (1.0 - t)).clamp(0.0, 1.0);
    let b = (1.0 - t * 1.2).clamp(0.0, 1.0);
    [r, g, b]
}
pub(super) fn magma_sample(t: f64) -> [f64; 3] {
    let r = (t * 1.1).clamp(0.0, 1.0);
    let g = (t * t * 0.9).clamp(0.0, 1.0);
    let b = if t < 0.5 { t * 1.5 } else { 1.5 - t * 1.5 };
    [r, g, b.clamp(0.0, 1.0)]
}
pub(super) fn inferno_sample(t: f64) -> [f64; 3] {
    let r = (t * 1.2).clamp(0.0, 1.0);
    let g = (t * t * 0.8).clamp(0.0, 1.0);
    let b = ((1.0 - t) * 0.4).clamp(0.0, 1.0);
    [r, g, b]
}
pub(super) fn turbo_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb((1.0 - t) * 0.667, 1.0, 1.0)
}
pub(super) fn rdbu_sample(t: f64) -> [f64; 3] {
    if t < 0.5 {
        let s = t * 2.0;
        [s * 0.7 + 0.3, s * 0.7 + 0.3, 1.0]
    } else {
        let s = (t - 0.5) * 2.0;
        [1.0, (1.0 - s) * 0.7 + 0.3, (1.0 - s) * 0.7 + 0.3]
    }
}
pub(super) fn spectral_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb(t * 0.833, 0.9, 0.9)
}
pub(super) fn coolwarm_sample(t: f64) -> [f64; 3] {
    if t < 0.5 {
        let s = (0.5 - t) * 2.0;
        [0.3 + s * 0.4, 0.3 + s * 0.4, 0.9]
    } else {
        let s = (t - 0.5) * 2.0;
        [0.9, 0.3 + (1.0 - s) * 0.4, 0.3 + (1.0 - s) * 0.4]
    }
}
pub(super) fn hot_sample(t: f64) -> [f64; 3] {
    let r = (t * 3.0).clamp(0.0, 1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    [r, g, b]
}
pub(super) fn jet_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb((1.0 - t) * 0.667, 1.0, 1.0)
}
/// Register all `viz` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
/// Wave-2 annotation pass fills in the `add_class` / `add_function` calls.
pub fn register_viz_module(parent: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = pyo3::types::PyModule::new(parent.py(), "viz")?;
    child.add_class::<ColormapKind>()?;
    child.add_class::<PyColormap>()?;
    child.add_class::<PyCamera>()?;
    child.add_class::<PyMaterial>()?;
    child.add_class::<PyMeshRenderer>()?;
    child.add_class::<PyParticleRenderer>()?;
    child.add_class::<TransferPoint>()?;
    child.add_class::<PyTransferFunction>()?;
    child.add_class::<PyVolumeRenderer>()?;
    child.add_class::<PySceneNode>()?;
    child.add_class::<PyLight>()?;
    child.add_class::<PySceneGraph>()?;
    child.add_class::<ToneMapping>()?;
    child.add_class::<PyPostProcessor>()?;
    child.add_class::<StressTensor>()?;
    child.add_class::<StressVisOutput>()?;
    child.add_class::<PyStressVisualizer>()?;
    child.add_class::<Streamline>()?;
    child.add_class::<PyStreamlineTracer>()?;
    child.add_class::<PyDebugOverlay>()?;
    child.add_function(wrap_pyfunction!(hsv_to_rgb, &child)?)?;
    child.add_function(wrap_pyfunction!(compute_normals_from_mesh, &child)?)?;
    child.add_function(wrap_pyfunction!(generate_sphere_mesh, &child)?)?;
    child.add_function(wrap_pyfunction!(generate_box_mesh, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_hsv_to_rgb_red() {
        let rgb = hsv_to_rgb(0.0, 1.0, 1.0);
        assert!((rgb[0] - 1.0).abs() < 1e-9);
        assert!(rgb[1] < 1e-9);
        assert!(rgb[2] < 1e-9);
    }
    #[test]
    fn test_hsv_to_rgb_green() {
        let rgb = hsv_to_rgb(1.0 / 3.0, 1.0, 1.0);
        assert!(rgb[1] > 0.99);
    }
    #[test]
    fn test_hsv_to_rgb_achromatic() {
        let rgb = hsv_to_rgb(0.0, 0.0, 0.7);
        assert!((rgb[0] - 0.7).abs() < 1e-9);
        assert!((rgb[1] - 0.7).abs() < 1e-9);
        assert!((rgb[2] - 0.7).abs() < 1e-9);
    }
    #[test]
    fn test_colormap_normalize() {
        let cm = PyColormap::new(ColormapKind::Viridis, 0.0, 10.0);
        assert!((cm.normalize(5.0) - 0.5).abs() < 1e-9);
        assert!((cm.normalize(-1.0) - 0.0).abs() < 1e-9);
        assert!((cm.normalize(11.0) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_colormap_viridis_alpha() {
        let cm = PyColormap::viridis();
        let rgba = cm.map_value(0.5);
        assert_eq!(rgba[3], 1.0);
    }
    #[test]
    fn test_colormap_all_kinds() {
        let kinds = [
            ColormapKind::Viridis,
            ColormapKind::Plasma,
            ColormapKind::Magma,
            ColormapKind::Inferno,
            ColormapKind::Turbo,
            ColormapKind::Greys,
            ColormapKind::RdBu,
            ColormapKind::Spectral,
            ColormapKind::Coolwarm,
            ColormapKind::Hot,
            ColormapKind::Jet,
        ];
        for kind in kinds {
            let cm = PyColormap::new(kind, 0.0, 1.0);
            let rgba = cm.map_value(0.5);
            for &c in &rgba {
                assert!((0.0..=1.0).contains(&c), "component out of range");
            }
        }
    }
    #[test]
    fn test_camera_default() {
        let cam = PyCamera::default_perspective();
        assert!(cam.distance_to_target() > 0.0);
    }
    #[test]
    fn test_camera_orbit() {
        let mut cam = PyCamera::default_perspective();
        let d0 = cam.distance_to_target();
        cam.orbit(0.1, 0.05);
        let d1 = cam.distance_to_target();
        assert!((d0 - d1).abs() < 1e-6);
    }
    #[test]
    fn test_camera_zoom() {
        let mut cam = PyCamera::default_perspective();
        let d0 = cam.distance_to_target();
        cam.zoom(0.5);
        let d1 = cam.distance_to_target();
        assert!((d1 - d0 * 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_camera_pan() {
        let mut cam = PyCamera::default_perspective();
        let pos0 = cam.position;
        cam.pan([1.0, 2.0, 3.0]);
        assert!((cam.position[0] - pos0[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_camera_view_direction() {
        let cam = PyCamera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0]);
        let d = cam.view_direction();
        assert!((d[2] + 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_compute_normals() {
        let verts = vec![0.0f64, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let idx = vec![0u32, 1, 2];
        let n = compute_normals_from_mesh(verts, idx);
        assert_eq!(n.len(), 9);
        for i in 0..3 {
            assert!((n[3 * i + 2] - 1.0).abs() < 1e-9);
        }
    }
    #[test]
    fn test_generate_sphere_mesh() {
        let (v, n, i) = generate_sphere_mesh(1.0, 8, 16);
        assert!(!v.is_empty());
        assert_eq!(v.len(), n.len());
        assert!(!i.is_empty());
    }
    #[test]
    fn test_generate_box_mesh() {
        let (v, n, i) = generate_box_mesh([1.0, 1.0, 1.0]);
        assert_eq!(v.len(), 24 * 3);
        assert_eq!(n.len(), 24 * 3);
        assert_eq!(i.len(), 36);
    }
    #[test]
    fn test_mesh_renderer_vertex_count() {
        let (v, _, i) = generate_sphere_mesh(1.0, 4, 8);
        let mat = PyMaterial::red();
        let mesh = PyMeshRenderer::new(v, i, mat);
        assert!(mesh.vertex_count() > 0);
        assert!(mesh.triangle_count() > 0);
    }
    #[test]
    fn test_mesh_renderer_apply_colormap() {
        let (v, _, i) = generate_sphere_mesh(1.0, 4, 8);
        let n = v.len() / 3;
        let mat = PyMaterial::blue();
        let mut mesh = PyMeshRenderer::new(v, i, mat);
        let scalars = vec![0.5; n];
        mesh.apply_colormap(scalars, PyColormap::viridis());
        assert!(mesh.scalars.is_some());
    }
    #[test]
    fn test_particle_renderer_count() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let pr = PyParticleRenderer::new(pos);
        assert_eq!(pr.particle_count(), 2);
    }
    #[test]
    fn test_particle_set_colors_from_scalars() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let mut pr = PyParticleRenderer::new(pos);
        pr.set_colors_from_scalars(vec![0.0, 1.0], PyColormap::viridis());
        assert_eq!(pr.colors.len(), 8);
    }
    #[test]
    fn test_transfer_function_sample() {
        let tf = PyTransferFunction::simple([0.0, 0.0, 1.0, 1.0], [1.0, 0.0, 0.0, 1.0]);
        let mid = tf.sample(0.5);
        assert!((mid[0] - 0.5).abs() < 1e-9);
        assert!((mid[2] - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_volume_renderer_sample() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let tf = PyTransferFunction::simple([0.0; 4], [1.0; 4]);
        let vr = PyVolumeRenderer::new(data, [2, 2, 2], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], tf);
        assert_eq!(vr.sample_at(0, 0, 0), 1.0);
        assert_eq!(vr.sample_at(1, 1, 1), 8.0);
        assert_eq!(vr.voxel_count(), 8);
    }
    #[test]
    fn test_scene_graph_add_node() {
        let mut sg = PySceneGraph::new();
        let id = sg.add_node("root".to_string());
        assert_eq!(id, 0);
        assert_eq!(sg.node_count(), 1);
    }
    #[test]
    fn test_scene_graph_visibility() {
        let mut sg = PySceneGraph::new();
        let id = sg.add_node("node".to_string());
        sg.set_visible(id, false);
        assert_eq!(sg.traverse_visible().len(), 0);
        sg.set_visible(id, true);
        assert_eq!(sg.traverse_visible().len(), 1);
    }
    #[test]
    fn test_post_processor_tone_mapping() {
        let pp = PyPostProcessor::new();
        let v = pp.apply_tone_mapping(1.0);
        assert!((0.0..=1.0).contains(&v));
    }
    #[test]
    fn test_post_processor_passthrough() {
        let pp = PyPostProcessor::passthrough();
        assert_eq!(pp.tone_mapping, ToneMapping::Linear);
        assert!(pp.ssao.is_none());
    }
    #[test]
    fn test_stress_tensor_von_mises() {
        let t = StressTensor::new(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let vm = t.von_mises();
        assert!((vm - 3.0_f64.sqrt()).abs() < 1e-9);
    }
    #[test]
    fn test_stress_tensor_hydrostatic() {
        let t = StressTensor::new(3.0, 3.0, 3.0, 0.0, 0.0, 0.0);
        assert!((t.hydrostatic() - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_stress_visualizer_compute() {
        let t = StressTensor::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let vis = PyStressVisualizer::new(vec![t]);
        let out = vis.compute(0).unwrap();
        assert!(out.von_mises >= 0.0);
    }
    #[test]
    fn test_streamline_tracer_trace() {
        let nx = 4u32;
        let ny = 4u32;
        let nz = 4u32;
        let n = (nx * ny * nz) as usize;
        let mut vf = vec![0.0f64; n * 3];
        for i in 0..n {
            vf[3 * i] = 1.0;
        }
        let tr = PyStreamlineTracer::new(vf, [nx, ny, nz], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let sl = tr.trace([0.1, 0.5, 0.5]);
        assert!(!sl.points.is_empty());
        assert!(sl.arc_length >= 0.0);
    }
    #[test]
    fn test_debug_overlay_draw() {
        let mut ov = PyDebugOverlay::new();
        ov.draw_sphere([0.0, 0.0, 0.0], 1.0, None);
        ov.draw_box([0.0; 3], [1.0; 3], None);
        ov.draw_arrow([0.0; 3], [1.0, 0.0, 0.0], None, 0.05);
        ov.draw_text_3d([0.0; 3], "hello".to_string(), None, 0.1);
        ov.draw_contact_point([0.0; 3], [0.0, 1.0, 0.0], 0.01, None);
        assert_eq!(ov.count(), 5);
        ov.clear();
        assert_eq!(ov.count(), 0);
    }
    #[test]
    fn test_mesh_renderer_render_to_buffer() {
        let (v, _, i) = generate_box_mesh([1.0, 1.0, 1.0]);
        let mat = PyMaterial::silver();
        let mesh = PyMeshRenderer::new(v, i, mat);
        let buf = mesh.render_to_buffer(4, 4, [0, 0, 0, 255]);
        assert_eq!(buf.len(), 4 * 4 * 4);
    }
    /// Output has exactly `width * height * 4` bytes.
    #[test]
    fn test_render_returns_correct_size() {
        let (v, _, i) = generate_box_mesh([1.0, 1.0, 1.0]);
        let mat = PyMaterial::red();
        let mesh = PyMeshRenderer::new(v, i, mat);
        let (w, h) = (128u32, 96u32);
        let buf = mesh.render_to_buffer(w, h, [0, 0, 0, 255]);
        assert_eq!(buf.len(), (w * h * 4) as usize);
    }
    /// After adding a mesh, at least one pixel must differ from the background.
    #[test]
    fn test_render_not_all_background() {
        let (v, _, i) = generate_box_mesh([1.0, 1.0, 1.0]);
        let mat = PyMaterial::red();
        let mesh = PyMeshRenderer::new(v, i, mat);
        let bg = [20u8, 20, 20, 255];
        let buf = mesh.render_to_buffer(128, 128, bg);
        let differs = buf
            .chunks(4)
            .any(|px| px[0] != bg[0] || px[1] != bg[1] || px[2] != bg[2]);
        assert!(
            differs,
            "render produced no pixels different from background"
        );
    }
    /// Two overlapping triangles in one mesh: the one with smaller z (front) must win.
    ///
    /// This test exercises `PyMeshRenderer::render_to_buffer` end-to-end.  We build
    /// a single mesh with 6 vertices: indices 0-2 form a red triangle at z=0 (front)
    /// and indices 3-5 form a blue triangle at z=2 (back), both covering the same XY
    /// region.  After rendering, every covered pixel must have red ≥ blue.
    #[test]
    fn test_render_depth_test_correct() {
        let (w, h) = (64u32, 64u32);
        let verts = vec![
            -0.5, -0.5, 0.0_f64, 0.5, -0.5, 0.0, 0.0, 0.5, 0.0, -0.5, -0.5, 2.0, 0.5, -0.5, 2.0,
            0.0, 0.5, 2.0,
        ];
        let indices = vec![0u32, 1, 2, 3, 4, 5];
        let mat = PyMaterial::new([1.0, 0.0, 0.0, 1.0]);
        let mesh_red = PyMeshRenderer::new(verts.clone(), indices.clone(), mat);
        let mat_blue = PyMaterial::new([0.0, 0.0, 1.0, 1.0]);
        let mesh_blue = PyMeshRenderer::new(verts, indices, mat_blue);
        let bg = [0u8, 0, 0, 255];
        let buf_red = mesh_red.render_to_buffer(w, h, bg);
        let buf_blue = mesh_blue.render_to_buffer(w, h, bg);
        let pixel_count = (w * h) as usize;
        let mut checked = false;
        for i in 0..pixel_count {
            let r = buf_red[4 * i];
            let b = buf_red[4 * i + 2];
            if r > 20 || b > 20 {
                assert!(
                    r >= b,
                    "red mesh: blue channel wins at idx {i}: r={r} b={b}"
                );
                checked = true;
            }
        }
        assert!(checked, "no non-background pixels found in red mesh render");
        let mut checked_blue = false;
        for i in 0..pixel_count {
            let r = buf_blue[4 * i];
            let b = buf_blue[4 * i + 2];
            if r > 20 || b > 20 {
                assert!(
                    b >= r,
                    "blue mesh: red channel wins at idx {i}: r={r} b={b}"
                );
                checked_blue = true;
            }
        }
        assert!(
            checked_blue,
            "no non-background pixels found in blue mesh render"
        );
    }
    /// Normal aligned to light → higher brightness than perpendicular normal.
    #[test]
    fn test_render_phong_shade_diffuse_increases_with_alignment() {
        use crate::rasterizer::{normalize3, phong_shade};
        let light_dir = normalize3([0.0, 0.0, 1.0]);
        let view_dir = [0.0, 0.0, 1.0_f64];
        let white = [1.0, 1.0, 1.0];
        let aligned = phong_shade([0.0, 0.0, 1.0], view_dir, light_dir, white, 0.0, 0.0, 1.0);
        let perp = phong_shade([1.0, 0.0, 0.0], view_dir, light_dir, white, 0.0, 0.0, 1.0);
        assert!(
            aligned[0] > perp[0],
            "aligned ({}) should be brighter than perpendicular ({})",
            aligned[0],
            perp[0]
        );
    }
    /// A single back-placed triangle should still render (z_buffer init to MAX).
    #[test]
    fn test_render_z_buffer_initialized() {
        let verts = vec![-0.5f64, -0.5, 100.0, 0.5, -0.5, 100.0, 0.0, 0.5, 100.0];
        let indices = vec![0u32, 1, 2];
        let mat = PyMaterial::new([0.0, 1.0, 0.0, 1.0]);
        let mut mesh = PyMeshRenderer::new(verts, indices, mat);
        mesh.recompute_normals();
        let bg = [0u8, 0, 0, 255];
        let buf = mesh.render_to_buffer(64, 64, bg);
        let has_green = buf
            .chunks(4)
            .any(|px| px[1] > 50 && px[0] < 50 && px[2] < 50);
        assert!(
            has_green,
            "back triangle not visible – z_buffer may not have been initialized to MAX"
        );
    }
}
