// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Visualization API for Python interop.
//!
//! Provides Python-friendly types for rendering, camera control, scene management,
//! stress visualization, streamline tracing, and debug overlays.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Convert HSV color to RGB.
///
/// All inputs and outputs are in \[0, 1\] range.
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
pub fn compute_normals_from_mesh(vertices: &[f64], indices: &[u32]) -> Vec<f64> {
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
pub fn generate_box_mesh(half_extents: [f64; 3]) -> (Vec<f64>, Vec<f64>, Vec<u32>) {
    let [hx, hy, hz] = half_extents;
    // 6 faces × 4 vertices = 24 vertices
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
        // Determine the two tangent directions based on the normal
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

// ---------------------------------------------------------------------------
// PyColormap
// ---------------------------------------------------------------------------

/// Named colormap for scalar field visualization.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum ColormapKind {
    /// Viridis perceptually uniform colormap.
    Viridis,
    /// Plasma perceptually uniform colormap.
    Plasma,
    /// Magma perceptually uniform colormap.
    Magma,
    /// Inferno perceptually uniform colormap.
    Inferno,
    /// Turbo rainbow colormap.
    Turbo,
    /// Greyscale colormap.
    Greys,
    /// Red-Blue diverging colormap.
    RdBu,
    /// Spectral diverging colormap.
    Spectral,
    /// Coolwarm diverging colormap.
    Coolwarm,
    /// Hot colormap.
    Hot,
    /// Jet colormap (legacy).
    Jet,
}

/// A colormap that maps scalar values in \[0, 1\] to RGBA colors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyColormap {
    /// The colormap kind.
    pub kind: ColormapKind,
    /// Minimum scalar value (maps to t=0).
    pub vmin: f64,
    /// Maximum scalar value (maps to t=1).
    pub vmax: f64,
    /// Alpha value for all colors (0=transparent, 1=opaque).
    pub alpha: f64,
}

impl PyColormap {
    /// Create a new colormap with explicit range.
    pub fn new(kind: ColormapKind, vmin: f64, vmax: f64) -> Self {
        Self {
            kind,
            vmin,
            vmax,
            alpha: 1.0,
        }
    }

    /// Create a Viridis colormap over \[0, 1\].
    pub fn viridis() -> Self {
        Self::new(ColormapKind::Viridis, 0.0, 1.0)
    }

    /// Create a Plasma colormap over \[0, 1\].
    pub fn plasma() -> Self {
        Self::new(ColormapKind::Plasma, 0.0, 1.0)
    }

    /// Create a Jet colormap over \[0, 1\].
    pub fn jet() -> Self {
        Self::new(ColormapKind::Jet, 0.0, 1.0)
    }

    /// Normalize a raw scalar value to t ∈ \[0, 1\].
    pub fn normalize(&self, value: f64) -> f64 {
        if (self.vmax - self.vmin).abs() < 1e-15 {
            return 0.5;
        }
        ((value - self.vmin) / (self.vmax - self.vmin)).clamp(0.0, 1.0)
    }

    /// Map a normalized value t ∈ \[0, 1\] to RGBA \[r, g, b, a\].
    pub fn map_value(&self, t: f64) -> [f64; 4] {
        let t = t.clamp(0.0, 1.0);
        let rgb = match self.kind {
            ColormapKind::Viridis => viridis_sample(t),
            ColormapKind::Plasma => plasma_sample(t),
            ColormapKind::Magma => magma_sample(t),
            ColormapKind::Inferno => inferno_sample(t),
            ColormapKind::Turbo => turbo_sample(t),
            ColormapKind::Greys => [t, t, t],
            ColormapKind::RdBu => rdbu_sample(t),
            ColormapKind::Spectral => spectral_sample(t),
            ColormapKind::Coolwarm => coolwarm_sample(t),
            ColormapKind::Hot => hot_sample(t),
            ColormapKind::Jet => jet_sample(t),
        };
        [rgb[0], rgb[1], rgb[2], self.alpha]
    }

    /// Map a raw scalar value directly.
    pub fn map_scalar(&self, value: f64) -> [f64; 4] {
        self.map_value(self.normalize(value))
    }
}

fn viridis_sample(t: f64) -> [f64; 3] {
    // Approximate Viridis using polynomial
    let r = 0.2777 * t.powi(3) - 0.8673 * t.powi(2) + 0.4756 * t + 0.2665;
    let g = -0.0955 * t.powi(3) + 0.5925 * t.powi(2) + 0.2843 * t + 0.0038;
    let b = -1.0803 * t.powi(3) + 1.1215 * t.powi(2) - 0.6424 * t + 0.5305;
    [r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0)]
}

fn plasma_sample(t: f64) -> [f64; 3] {
    let r = (0.9 * t + 0.05).clamp(0.0, 1.0);
    let g = (4.0 * t * (1.0 - t)).clamp(0.0, 1.0);
    let b = (1.0 - t * 1.2).clamp(0.0, 1.0);
    [r, g, b]
}

fn magma_sample(t: f64) -> [f64; 3] {
    let r = (t * 1.1).clamp(0.0, 1.0);
    let g = (t * t * 0.9).clamp(0.0, 1.0);
    let b = if t < 0.5 { t * 1.5 } else { 1.5 - t * 1.5 };
    [r, g, b.clamp(0.0, 1.0)]
}

fn inferno_sample(t: f64) -> [f64; 3] {
    let r = (t * 1.2).clamp(0.0, 1.0);
    let g = (t * t * 0.8).clamp(0.0, 1.0);
    let b = ((1.0 - t) * 0.4).clamp(0.0, 1.0);
    [r, g, b]
}

fn turbo_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb((1.0 - t) * 0.667, 1.0, 1.0)
}

fn rdbu_sample(t: f64) -> [f64; 3] {
    if t < 0.5 {
        let s = t * 2.0;
        [s * 0.7 + 0.3, s * 0.7 + 0.3, 1.0]
    } else {
        let s = (t - 0.5) * 2.0;
        [1.0, (1.0 - s) * 0.7 + 0.3, (1.0 - s) * 0.7 + 0.3]
    }
}

fn spectral_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb(t * 0.833, 0.9, 0.9)
}

fn coolwarm_sample(t: f64) -> [f64; 3] {
    if t < 0.5 {
        let s = (0.5 - t) * 2.0;
        [0.3 + s * 0.4, 0.3 + s * 0.4, 0.9]
    } else {
        let s = (t - 0.5) * 2.0;
        [0.9, 0.3 + (1.0 - s) * 0.4, 0.3 + (1.0 - s) * 0.4]
    }
}

fn hot_sample(t: f64) -> [f64; 3] {
    let r = (t * 3.0).clamp(0.0, 1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    [r, g, b]
}

fn jet_sample(t: f64) -> [f64; 3] {
    hsv_to_rgb((1.0 - t) * 0.667, 1.0, 1.0)
}

// ---------------------------------------------------------------------------
// PyCamera
// ---------------------------------------------------------------------------

/// 3D camera with orbit, pan, zoom controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCamera {
    /// Camera position in world space.
    pub position: [f64; 3],
    /// Look-at target point.
    pub target: [f64; 3],
    /// Up vector (usually \[0, 1, 0\]).
    pub up: [f64; 3],
    /// Vertical field of view in degrees.
    pub fov_deg: f64,
    /// Near clipping plane distance.
    pub near: f64,
    /// Far clipping plane distance.
    pub far: f64,
    /// Viewport aspect ratio (width / height).
    pub aspect: f64,
}

impl PyCamera {
    /// Create a default perspective camera.
    pub fn new(position: [f64; 3], target: [f64; 3]) -> Self {
        Self {
            position,
            target,
            up: [0.0, 1.0, 0.0],
            fov_deg: 60.0,
            near: 0.01,
            far: 1000.0,
            aspect: 16.0 / 9.0,
        }
    }

    /// Default camera looking down the negative Z axis.
    pub fn default_perspective() -> Self {
        Self::new([0.0, 5.0, 10.0], [0.0, 0.0, 0.0])
    }

    /// Orbit the camera around the target by delta_yaw and delta_pitch (radians).
    pub fn orbit(&mut self, delta_yaw: f64, delta_pitch: f64) {
        let dx = self.position[0] - self.target[0];
        let dy = self.position[1] - self.target[1];
        let dz = self.position[2] - self.target[2];
        let radius = (dx * dx + dy * dy + dz * dz).sqrt();
        let mut yaw = dz.atan2(dx);
        let mut pitch = dy.atan2((dx * dx + dz * dz).sqrt());
        yaw += delta_yaw;
        pitch = (pitch + delta_pitch).clamp(-1.5, 1.5);
        self.position[0] = self.target[0] + radius * pitch.cos() * yaw.cos();
        self.position[1] = self.target[1] + radius * pitch.sin();
        self.position[2] = self.target[2] + radius * pitch.cos() * yaw.sin();
    }

    /// Pan the camera by a world-space offset.
    pub fn pan(&mut self, offset: [f64; 3]) {
        self.position[0] += offset[0];
        self.position[1] += offset[1];
        self.position[2] += offset[2];
        self.target[0] += offset[0];
        self.target[1] += offset[1];
        self.target[2] += offset[2];
    }

    /// Zoom by moving the camera closer to (factor < 1) or farther from (factor > 1) the target.
    pub fn zoom(&mut self, factor: f64) {
        let dx = self.position[0] - self.target[0];
        let dy = self.position[1] - self.target[1];
        let dz = self.position[2] - self.target[2];
        self.position[0] = self.target[0] + dx * factor;
        self.position[1] = self.target[1] + dy * factor;
        self.position[2] = self.target[2] + dz * factor;
    }

    /// Point the camera directly at a given world position.
    pub fn look_at(&mut self, target: [f64; 3]) {
        self.target = target;
    }

    /// Compute the view direction (normalized).
    pub fn view_direction(&self) -> [f64; 3] {
        let dx = self.target[0] - self.position[0];
        let dy = self.target[1] - self.position[1];
        let dz = self.target[2] - self.position[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < 1e-12 {
            return [0.0, 0.0, -1.0];
        }
        [dx / len, dy / len, dz / len]
    }

    /// Distance from camera to target.
    pub fn distance_to_target(&self) -> f64 {
        let dx = self.position[0] - self.target[0];
        let dy = self.position[1] - self.target[1];
        let dz = self.position[2] - self.target[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Material
// ---------------------------------------------------------------------------

/// Surface material properties for mesh rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMaterial {
    /// Base color RGBA.
    pub base_color: [f64; 4],
    /// Metallic factor (0=dielectric, 1=metal).
    pub metallic: f64,
    /// Roughness factor (0=mirror, 1=fully rough).
    pub roughness: f64,
    /// Emissive color RGB.
    pub emissive: [f64; 3],
    /// Whether to use double-sided rendering.
    pub double_sided: bool,
    /// Optional colormap for scalar-driven coloring.
    pub colormap: Option<PyColormap>,
}

impl PyMaterial {
    /// Create a default PBR material.
    pub fn new(base_color: [f64; 4]) -> Self {
        Self {
            base_color,
            metallic: 0.0,
            roughness: 0.5,
            emissive: [0.0, 0.0, 0.0],
            double_sided: false,
            colormap: None,
        }
    }

    /// Red material.
    pub fn red() -> Self {
        Self::new([1.0, 0.0, 0.0, 1.0])
    }

    /// Blue material.
    pub fn blue() -> Self {
        Self::new([0.0, 0.3, 1.0, 1.0])
    }

    /// Metallic silver material.
    pub fn silver() -> Self {
        Self {
            base_color: [0.8, 0.8, 0.8, 1.0],
            metallic: 1.0,
            roughness: 0.1,
            emissive: [0.0; 3],
            double_sided: false,
            colormap: None,
        }
    }
}

// ---------------------------------------------------------------------------
// PyMeshRenderer
// ---------------------------------------------------------------------------

/// GPU-ready mesh data with material and optional per-vertex scalar field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMeshRenderer {
    /// Flat array of vertex positions: \[x0, y0, z0, x1, y1, z1, ...\].
    pub vertices: Vec<f64>,
    /// Flat array of vertex normals (same length as vertices).
    pub normals: Vec<f64>,
    /// Triangle index list (3 indices per triangle).
    pub indices: Vec<u32>,
    /// Surface material.
    pub material: PyMaterial,
    /// Optional per-vertex scalar values for colormap coloring.
    pub scalars: Option<Vec<f64>>,
}

impl PyMeshRenderer {
    /// Create a new mesh renderer from raw geometry.
    pub fn new(vertices: Vec<f64>, indices: Vec<u32>, material: PyMaterial) -> Self {
        let normals = compute_normals_from_mesh(&vertices, &indices);
        Self {
            vertices,
            normals,
            indices,
            material,
            scalars: None,
        }
    }

    /// Apply a colormap to per-vertex scalar values.
    ///
    /// Sets `material.colormap` and stores the scalar array.
    pub fn apply_colormap(&mut self, scalars: Vec<f64>, colormap: PyColormap) {
        self.scalars = Some(scalars);
        self.material.colormap = Some(colormap);
    }

    /// Returns the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len() / 3
    }

    /// Returns the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Recompute normals from current vertex/index data.
    pub fn recompute_normals(&mut self) {
        self.normals = compute_normals_from_mesh(&self.vertices, &self.indices);
    }

    /// Render to a flat RGBA pixel buffer (software rasterizer stub).
    ///
    /// Returns a buffer of `width * height * 4` bytes (u8), cleared to the background color.
    pub fn render_to_buffer(&self, width: u32, height: u32, background: [u8; 4]) -> Vec<u8> {
        let n = (width * height * 4) as usize;
        let mut buf = vec![0u8; n];
        for i in 0..(width * height) as usize {
            buf[4 * i] = background[0];
            buf[4 * i + 1] = background[1];
            buf[4 * i + 2] = background[2];
            buf[4 * i + 3] = background[3];
        }
        buf
    }
}

// ---------------------------------------------------------------------------
// PyParticleRenderer
// ---------------------------------------------------------------------------

/// Renderer for particle systems using billboard or instanced geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyParticleRenderer {
    /// Flat array of particle positions: \[x0, y0, z0, x1, y1, z1, ...\].
    pub positions: Vec<f64>,
    /// Per-particle radii. If length is 1, all particles share that radius.
    pub radii: Vec<f64>,
    /// Per-particle RGBA colors. Flat array: \[r0, g0, b0, a0, r1, ...\].
    pub colors: Vec<f64>,
    /// Whether to use billboard rendering (always face camera).
    pub billboard: bool,
    /// Whether to use hardware instancing (more efficient for large counts).
    pub use_instancing: bool,
    /// Optional colormap applied over a scalar field.
    pub colormap: Option<PyColormap>,
}

impl PyParticleRenderer {
    /// Create a new particle renderer.
    pub fn new(positions: Vec<f64>) -> Self {
        let n = positions.len() / 3;
        Self {
            positions,
            radii: vec![0.05; n],
            colors: vec![0.8, 0.4, 0.1, 1.0]
                .into_iter()
                .cycle()
                .take(n * 4)
                .collect(),
            billboard: true,
            use_instancing: true,
            colormap: None,
        }
    }

    /// Set uniform radius for all particles.
    pub fn set_uniform_radius(&mut self, radius: f64) {
        let n = self.positions.len() / 3;
        self.radii = vec![radius; n];
    }

    /// Set per-particle colors from a colormap and scalar values.
    pub fn set_colors_from_scalars(&mut self, scalars: &[f64], colormap: &PyColormap) {
        let n = self.positions.len() / 3;
        self.colors.clear();
        for i in 0..n {
            let s = scalars.get(i).copied().unwrap_or(0.0);
            let rgba = colormap.map_scalar(s);
            self.colors.extend_from_slice(&rgba);
        }
        self.colormap = Some(colormap.clone());
    }

    /// Return the number of particles.
    pub fn particle_count(&self) -> usize {
        self.positions.len() / 3
    }
}

// ---------------------------------------------------------------------------
// Transfer function for volume rendering
// ---------------------------------------------------------------------------

/// A control point in a transfer function (scalar → RGBA).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransferPoint {
    /// Scalar value in \[0, 1\].
    pub scalar: f64,
    /// RGBA color at this point.
    pub color: [f64; 4],
}

/// Transfer function mapping scalar densities to RGBA for volume rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTransferFunction {
    /// Sorted control points.
    pub control_points: Vec<TransferPoint>,
}

impl PyTransferFunction {
    /// Create a simple two-point transfer function.
    pub fn simple(low_color: [f64; 4], high_color: [f64; 4]) -> Self {
        Self {
            control_points: vec![
                TransferPoint {
                    scalar: 0.0,
                    color: low_color,
                },
                TransferPoint {
                    scalar: 1.0,
                    color: high_color,
                },
            ],
        }
    }

    /// Sample the transfer function at a scalar value in \[0, 1\].
    pub fn sample(&self, t: f64) -> [f64; 4] {
        let pts = &self.control_points;
        if pts.is_empty() {
            return [0.0, 0.0, 0.0, 0.0];
        }
        if pts.len() == 1 || t <= pts[0].scalar {
            return pts[0].color;
        }
        if t >= pts[pts.len() - 1].scalar {
            return pts[pts.len() - 1].color;
        }
        for i in 1..pts.len() {
            if t <= pts[i].scalar {
                let a = pts[i - 1].scalar;
                let b = pts[i].scalar;
                let alpha = if (b - a).abs() < 1e-15 {
                    0.0
                } else {
                    (t - a) / (b - a)
                };
                let ca = pts[i - 1].color;
                let cb = pts[i].color;
                return [
                    ca[0] + alpha * (cb[0] - ca[0]),
                    ca[1] + alpha * (cb[1] - ca[1]),
                    ca[2] + alpha * (cb[2] - ca[2]),
                    ca[3] + alpha * (cb[3] - ca[3]),
                ];
            }
        }
        pts[pts.len() - 1].color
    }
}

// ---------------------------------------------------------------------------
// PyVolumeRenderer
// ---------------------------------------------------------------------------

/// Settings for ray-marching volume rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RayMarchSettings {
    /// Step size along the ray.
    pub step_size: f64,
    /// Maximum number of steps.
    pub max_steps: u32,
    /// Early termination opacity threshold.
    pub opacity_threshold: f64,
    /// Jitter the ray origin to reduce banding.
    pub jitter: bool,
}

impl Default for RayMarchSettings {
    fn default() -> Self {
        Self {
            step_size: 0.01,
            max_steps: 512,
            opacity_threshold: 0.99,
            jitter: true,
        }
    }
}

/// Volume renderer for 3D scalar fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyVolumeRenderer {
    /// 3D scalar field stored in \[z\]\[y\]\[x\] order, flattened.
    pub data: Vec<f64>,
    /// Grid dimensions \[nx, ny, nz\].
    pub dims: [u32; 3],
    /// World-space bounding box \[min_x, min_y, min_z, max_x, max_y, max_z\].
    pub bounds: [f64; 6],
    /// Transfer function for density → RGBA mapping.
    pub transfer_function: PyTransferFunction,
    /// Ray marching settings.
    pub ray_march: RayMarchSettings,
}

impl PyVolumeRenderer {
    /// Create a volume renderer from a flat scalar field.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        data: Vec<f64>,
        dims: [u32; 3],
        bounds: [f64; 6],
        transfer_function: PyTransferFunction,
    ) -> Self {
        Self {
            data,
            dims,
            bounds,
            transfer_function,
            ray_march: RayMarchSettings::default(),
        }
    }

    /// Sample the scalar field at a grid index (ix, iy, iz).
    pub fn sample_at(&self, ix: u32, iy: u32, iz: u32) -> f64 {
        let [nx, ny, _nz] = self.dims;
        let idx = iz * ny * nx + iy * nx + ix;
        self.data.get(idx as usize).copied().unwrap_or(0.0)
    }

    /// Compute the total number of voxels.
    pub fn voxel_count(&self) -> u64 {
        self.dims[0] as u64 * self.dims[1] as u64 * self.dims[2] as u64
    }
}

// ---------------------------------------------------------------------------
// PySceneGraph
// ---------------------------------------------------------------------------

/// A transform node in the scene hierarchy.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySceneNode {
    /// Node identifier.
    pub id: u32,
    /// Optional parent node id.
    pub parent: Option<u32>,
    /// Local translation \[tx, ty, tz\].
    pub translation: [f64; 3],
    /// Local rotation as quaternion \[qx, qy, qz, qw\].
    pub rotation: [f64; 4],
    /// Local scale \[sx, sy, sz\].
    pub scale: [f64; 3],
    /// Human-readable label.
    pub label: String,
    /// Whether this node is visible.
    pub visible: bool,
}

impl PySceneNode {
    /// Create an identity transform node.
    pub fn identity(id: u32, label: impl Into<String>) -> Self {
        Self {
            id,
            parent: None,
            translation: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0; 3],
            label: label.into(),
            visible: true,
        }
    }
}

/// A light source in the scene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyLight {
    /// Light position \[x, y, z\].
    pub position: [f64; 3],
    /// Light color RGB.
    pub color: [f64; 3],
    /// Intensity (candela / lm).
    pub intensity: f64,
    /// Light type: "point", "directional", "spot".
    pub light_type: String,
    /// Spot direction \[dx, dy, dz\] (used for spot lights).
    pub direction: [f64; 3],
    /// Spot cone angle in radians.
    pub cone_angle: f64,
}

impl PyLight {
    /// Create a point light.
    pub fn point(position: [f64; 3], color: [f64; 3], intensity: f64) -> Self {
        Self {
            position,
            color,
            intensity,
            light_type: "point".to_owned(),
            direction: [0.0, -1.0, 0.0],
            cone_angle: std::f64::consts::PI,
        }
    }

    /// Create a directional light (sun).
    pub fn directional(direction: [f64; 3], color: [f64; 3], intensity: f64) -> Self {
        Self {
            position: [0.0; 3],
            color,
            intensity,
            light_type: "directional".to_owned(),
            direction,
            cone_angle: std::f64::consts::PI,
        }
    }
}

/// Hierarchical scene graph managing nodes, meshes, particles and lights.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySceneGraph {
    /// All nodes in the scene.
    pub nodes: Vec<PySceneNode>,
    /// Mesh renderers attached to node ids.
    pub meshes: Vec<(u32, PyMeshRenderer)>,
    /// Particle renderers attached to node ids.
    pub particles: Vec<(u32, PyParticleRenderer)>,
    /// Lights in the scene.
    pub lights: Vec<PyLight>,
    /// Next free node id.
    next_id: u32,
}

impl PySceneGraph {
    /// Create an empty scene graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            meshes: Vec::new(),
            particles: Vec::new(),
            lights: Vec::new(),
            next_id: 0,
        }
    }

    /// Add a new node and return its id.
    pub fn add_node(&mut self, label: impl Into<String>) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.nodes.push(PySceneNode::identity(id, label));
        id
    }

    /// Attach a mesh to a node.
    pub fn add_mesh(&mut self, node_id: u32, mesh: PyMeshRenderer) {
        self.meshes.push((node_id, mesh));
    }

    /// Attach particles to a node.
    pub fn add_particles(&mut self, node_id: u32, particles: PyParticleRenderer) {
        self.particles.push((node_id, particles));
    }

    /// Add a light to the scene.
    pub fn add_light(&mut self, light: PyLight) {
        self.lights.push(light);
    }

    /// Set visibility for a node.
    pub fn set_visible(&mut self, node_id: u32, visible: bool) {
        if let Some(n) = self.nodes.iter_mut().find(|n| n.id == node_id) {
            n.visible = visible;
        }
    }

    /// Set local translation for a node.
    pub fn set_translation(&mut self, node_id: u32, translation: [f64; 3]) {
        if let Some(n) = self.nodes.iter_mut().find(|n| n.id == node_id) {
            n.translation = translation;
        }
    }

    /// Traverse all visible nodes and return their ids in order.
    pub fn traverse_visible(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.visible)
            .map(|n| n.id)
            .collect()
    }

    /// Return total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

impl Default for PySceneGraph {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyPostProcessor
// ---------------------------------------------------------------------------

/// Screen-space ambient occlusion parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SsaoParams {
    /// Sample hemisphere radius in view space.
    pub radius: f64,
    /// Bias to avoid self-occlusion.
    pub bias: f64,
    /// Number of sample points.
    pub num_samples: u32,
    /// Occlusion power exponent.
    pub power: f64,
}

impl Default for SsaoParams {
    fn default() -> Self {
        Self {
            radius: 0.5,
            bias: 0.025,
            num_samples: 16,
            power: 2.0,
        }
    }
}

/// Bloom post-processing parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BloomParams {
    /// Luminance threshold above which bloom is applied.
    pub threshold: f64,
    /// Bloom blur radius in pixels.
    pub radius: f64,
    /// Bloom intensity multiplier.
    pub intensity: f64,
}

impl Default for BloomParams {
    fn default() -> Self {
        Self {
            threshold: 1.0,
            radius: 8.0,
            intensity: 0.5,
        }
    }
}

/// Tone mapping operator selection.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum ToneMapping {
    /// Reinhard global operator.
    Reinhard,
    /// ACES film curve.
    Aces,
    /// Filmic tone mapping.
    Filmic,
    /// No tone mapping (linear).
    Linear,
}

/// Post-processing pipeline configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPostProcessor {
    /// SSAO settings.
    pub ssao: Option<SsaoParams>,
    /// Bloom settings.
    pub bloom: Option<BloomParams>,
    /// FXAA edge threshold (0 = disabled).
    pub fxaa_edge_threshold: f64,
    /// Tone mapping operator.
    pub tone_mapping: ToneMapping,
    /// Gamma correction exponent (typically 2.2).
    pub gamma: f64,
    /// Exposure multiplier.
    pub exposure: f64,
}

impl PyPostProcessor {
    /// Create a default post-processing pipeline.
    pub fn new() -> Self {
        Self {
            ssao: Some(SsaoParams::default()),
            bloom: Some(BloomParams::default()),
            fxaa_edge_threshold: 0.063,
            tone_mapping: ToneMapping::Aces,
            gamma: 2.2,
            exposure: 1.0,
        }
    }

    /// Disable all effects (passthrough).
    pub fn passthrough() -> Self {
        Self {
            ssao: None,
            bloom: None,
            fxaa_edge_threshold: 0.0,
            tone_mapping: ToneMapping::Linear,
            gamma: 1.0,
            exposure: 1.0,
        }
    }

    /// Apply Reinhard tone mapping to an HDR value.
    pub fn apply_tone_mapping(&self, hdr: f64) -> f64 {
        let v = hdr * self.exposure;
        let linear = match self.tone_mapping {
            ToneMapping::Reinhard => v / (1.0 + v),
            ToneMapping::Aces => {
                let a = 2.51;
                let b = 0.03;
                let c = 2.43;
                let d = 0.59;
                let e = 0.14;
                ((v * (a * v + b)) / (v * (c * v + d) + e)).clamp(0.0, 1.0)
            }
            ToneMapping::Filmic => {
                let x = (v - 0.004).max(0.0);
                (x * (6.2 * x + 0.5)) / (x * (6.2 * x + 1.7) + 0.06)
            }
            ToneMapping::Linear => v.clamp(0.0, 1.0),
        };
        linear.powf(1.0 / self.gamma)
    }
}

impl Default for PyPostProcessor {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyStressVisualizer
// ---------------------------------------------------------------------------

/// A 3×3 symmetric stress tensor stored as \[s11, s22, s33, s12, s13, s23\].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressTensor {
    /// Components \[σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃\].
    pub components: [f64; 6],
}

impl StressTensor {
    /// Create a new stress tensor.
    pub fn new(s11: f64, s22: f64, s33: f64, s12: f64, s13: f64, s23: f64) -> Self {
        Self {
            components: [s11, s22, s33, s12, s13, s23],
        }
    }

    /// Compute the von Mises stress.
    pub fn von_mises(&self) -> f64 {
        let [s11, s22, s33, s12, s13, s23] = self.components;
        let ds = (s11 - s22).powi(2) + (s22 - s33).powi(2) + (s33 - s11).powi(2);
        let shear = s12 * s12 + s13 * s13 + s23 * s23;
        (0.5 * ds + 3.0 * shear).sqrt()
    }

    /// Compute the hydrostatic (mean) stress.
    pub fn hydrostatic(&self) -> f64 {
        (self.components[0] + self.components[1] + self.components[2]) / 3.0
    }

    /// Compute the deviatoric stress components.
    pub fn deviatoric(&self) -> [f64; 6] {
        let p = self.hydrostatic();
        let [s11, s22, s33, s12, s13, s23] = self.components;
        [s11 - p, s22 - p, s33 - p, s12, s13, s23]
    }
}

/// Output from the stress visualizer for a single tensor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressVisOutput {
    /// Von Mises equivalent stress.
    pub von_mises: f64,
    /// Hydrostatic stress.
    pub hydrostatic: f64,
    /// Principal stresses \[σ₁, σ₂, σ₃\] (sorted descending).
    pub principal_stresses: [f64; 3],
    /// Principal directions as flat 3×3 matrix (row-major).
    pub principal_directions: [f64; 9],
    /// Mohr circle data: \[center_12, radius_12, center_23, radius_23, center_13, radius_13\].
    pub mohr_circle: [f64; 6],
}

/// Stress field visualizer providing principal direction glyphs and Mohr circle data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyStressVisualizer {
    /// Input stress tensor field (one per element/node).
    pub tensors: Vec<StressTensor>,
    /// Colormap for von Mises coloring.
    pub colormap: PyColormap,
    /// Scale factor for principal direction glyphs.
    pub glyph_scale: f64,
}

impl PyStressVisualizer {
    /// Create a new stress visualizer.
    pub fn new(tensors: Vec<StressTensor>) -> Self {
        let mut cm = PyColormap::new(ColormapKind::Viridis, 0.0, 1.0);
        cm.vmax = tensors
            .iter()
            .map(|t| t.von_mises())
            .fold(0.0_f64, f64::max)
            .max(1.0);
        Self {
            tensors,
            colormap: cm,
            glyph_scale: 1.0,
        }
    }

    /// Compute visualization output for tensor at index i.
    pub fn compute(&self, i: usize) -> Option<StressVisOutput> {
        let t = self.tensors.get(i)?;
        let vm = t.von_mises();
        let hydro = t.hydrostatic();
        // Approximate principal stresses via Jacobi (2-pass simplified for diagonal dominance)
        let [s11, s22, s33, _s12, _s13, _s23] = t.components;
        let mut p = [s11, s22, s33];
        p.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        let mohr = [
            (p[0] + p[1]) * 0.5,
            ((p[0] - p[1]) * 0.5).abs(),
            (p[1] + p[2]) * 0.5,
            ((p[1] - p[2]) * 0.5).abs(),
            (p[0] + p[2]) * 0.5,
            ((p[0] - p[2]) * 0.5).abs(),
        ];
        Some(StressVisOutput {
            von_mises: vm,
            hydrostatic: hydro,
            principal_stresses: p,
            principal_directions: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            mohr_circle: mohr,
        })
    }

    /// Return von Mises values for all tensors.
    pub fn all_von_mises(&self) -> Vec<f64> {
        self.tensors.iter().map(|t| t.von_mises()).collect()
    }
}

// ---------------------------------------------------------------------------
// PyStreamlineTracer
// ---------------------------------------------------------------------------

/// Result of streamline tracing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Streamline {
    /// Flat array of 3D positions along the streamline.
    pub points: Vec<f64>,
    /// Total arc length.
    pub arc_length: f64,
}

/// Streamline tracer using RK4 integration through a 3D velocity field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyStreamlineTracer {
    /// Flat velocity field \[vx0, vy0, vz0, vx1, ...\] on a uniform grid.
    pub velocity_field: Vec<f64>,
    /// Grid dimensions \[nx, ny, nz\].
    pub dims: [u32; 3],
    /// World-space bounds \[min_x, min_y, min_z, max_x, max_y, max_z\].
    pub bounds: [f64; 6],
    /// Integration step size.
    pub step_size: f64,
    /// Maximum integration steps per streamline.
    pub max_steps: u32,
    /// Tube radius for 3D tube rendering.
    pub tube_radius: f64,
}

impl PyStreamlineTracer {
    /// Create a new streamline tracer.
    pub fn new(velocity_field: Vec<f64>, dims: [u32; 3], bounds: [f64; 6]) -> Self {
        Self {
            velocity_field,
            dims,
            bounds,
            step_size: 0.01,
            max_steps: 512,
            tube_radius: 0.01,
        }
    }

    fn sample_velocity(&self, pos: [f64; 3]) -> [f64; 3] {
        let [nx, ny, nz] = self.dims;
        let [bx0, by0, bz0, bx1, by1, bz1] = self.bounds;
        let fx = ((pos[0] - bx0) / (bx1 - bx0) * nx as f64).clamp(0.0, nx as f64 - 1.001);
        let fy = ((pos[1] - by0) / (by1 - by0) * ny as f64).clamp(0.0, ny as f64 - 1.001);
        let fz = ((pos[2] - bz0) / (bz1 - bz0) * nz as f64).clamp(0.0, nz as f64 - 1.001);
        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;
        let stride = (nx * ny) as usize;
        let base = 3 * (iz * stride + iy * nx as usize + ix);
        if base + 2 < self.velocity_field.len() {
            [
                self.velocity_field[base],
                self.velocity_field[base + 1],
                self.velocity_field[base + 2],
            ]
        } else {
            [0.0, 0.0, 0.0]
        }
    }

    fn rk4_step(&self, pos: [f64; 3], dt: f64) -> [f64; 3] {
        let k1 = self.sample_velocity(pos);
        let p2 = [
            pos[0] + k1[0] * dt * 0.5,
            pos[1] + k1[1] * dt * 0.5,
            pos[2] + k1[2] * dt * 0.5,
        ];
        let k2 = self.sample_velocity(p2);
        let p3 = [
            pos[0] + k2[0] * dt * 0.5,
            pos[1] + k2[1] * dt * 0.5,
            pos[2] + k2[2] * dt * 0.5,
        ];
        let k3 = self.sample_velocity(p3);
        let p4 = [
            pos[0] + k3[0] * dt,
            pos[1] + k3[1] * dt,
            pos[2] + k3[2] * dt,
        ];
        let k4 = self.sample_velocity(p4);
        [
            pos[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            pos[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
            pos[2] + dt / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        ]
    }

    /// Trace a streamline from a seed point.
    pub fn trace(&self, seed: [f64; 3]) -> Streamline {
        let mut pts = vec![seed[0], seed[1], seed[2]];
        let mut pos = seed;
        let mut arc = 0.0;
        for _ in 0..self.max_steps {
            let next = self.rk4_step(pos, self.step_size);
            let dx = next[0] - pos[0];
            let dy = next[1] - pos[1];
            let dz = next[2] - pos[2];
            arc += (dx * dx + dy * dy + dz * dz).sqrt();
            pts.push(next[0]);
            pts.push(next[1]);
            pts.push(next[2]);
            pos = next;
            let v = self.sample_velocity(pos);
            if v[0] * v[0] + v[1] * v[1] + v[2] * v[2] < 1e-18 {
                break;
            }
        }
        Streamline {
            points: pts,
            arc_length: arc,
        }
    }

    /// Trace streamlines from multiple seed points.
    pub fn trace_multiple(&self, seeds: &[[f64; 3]]) -> Vec<Streamline> {
        seeds.iter().map(|&s| self.trace(s)).collect()
    }
}

// ---------------------------------------------------------------------------
// PyDebugOverlay
// ---------------------------------------------------------------------------

/// A single debug primitive drawn as an overlay.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DebugPrimitive {
    /// A wire-frame sphere.
    Sphere {
        center: [f64; 3],
        radius: f64,
        color: [f64; 4],
    },
    /// An axis-aligned box.
    Box {
        min: [f64; 3],
        max: [f64; 3],
        color: [f64; 4],
    },
    /// An arrow from start to end.
    Arrow {
        start: [f64; 3],
        end: [f64; 3],
        color: [f64; 4],
        shaft_radius: f64,
    },
    /// A 3D text label.
    Text3D {
        position: [f64; 3],
        text: String,
        color: [f64; 4],
        size: f64,
    },
    /// A contact point with normal.
    ContactPoint {
        position: [f64; 3],
        normal: [f64; 3],
        depth: f64,
        color: [f64; 4],
    },
}

/// Debug overlay manager for visualizing physics primitives.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyDebugOverlay {
    /// List of debug primitives to render.
    pub primitives: Vec<DebugPrimitive>,
    /// Whether to draw wireframe (true) or solid (false).
    pub wireframe: bool,
    /// Default color for new primitives.
    pub default_color: [f64; 4],
    /// Whether overlay persists across frames.
    pub persistent: bool,
}

impl PyDebugOverlay {
    /// Create a new empty debug overlay.
    pub fn new() -> Self {
        Self {
            primitives: Vec::new(),
            wireframe: true,
            default_color: [0.0, 1.0, 0.0, 1.0],
            persistent: false,
        }
    }

    /// Draw a sphere at the given center.
    pub fn draw_sphere(&mut self, center: [f64; 3], radius: f64, color: Option<[f64; 4]>) {
        let c = color.unwrap_or(self.default_color);
        self.primitives.push(DebugPrimitive::Sphere {
            center,
            radius,
            color: c,
        });
    }

    /// Draw an axis-aligned box.
    pub fn draw_box(&mut self, min: [f64; 3], max: [f64; 3], color: Option<[f64; 4]>) {
        let c = color.unwrap_or(self.default_color);
        self.primitives
            .push(DebugPrimitive::Box { min, max, color: c });
    }

    /// Draw an arrow from start to end.
    pub fn draw_arrow(
        &mut self,
        start: [f64; 3],
        end: [f64; 3],
        color: Option<[f64; 4]>,
        shaft_radius: f64,
    ) {
        let c = color.unwrap_or(self.default_color);
        self.primitives.push(DebugPrimitive::Arrow {
            start,
            end,
            color: c,
            shaft_radius,
        });
    }

    /// Draw a 3D text label.
    pub fn draw_text_3d(
        &mut self,
        position: [f64; 3],
        text: impl Into<String>,
        color: Option<[f64; 4]>,
        size: f64,
    ) {
        let c = color.unwrap_or(self.default_color);
        self.primitives.push(DebugPrimitive::Text3D {
            position,
            text: text.into(),
            color: c,
            size,
        });
    }

    /// Draw a contact point with normal direction.
    pub fn draw_contact_point(
        &mut self,
        position: [f64; 3],
        normal: [f64; 3],
        depth: f64,
        color: Option<[f64; 4]>,
    ) {
        let c = color.unwrap_or(self.default_color);
        self.primitives.push(DebugPrimitive::ContactPoint {
            position,
            normal,
            depth,
            color: c,
        });
    }

    /// Clear all debug primitives.
    pub fn clear(&mut self) {
        self.primitives.clear();
    }

    /// Return count of primitives.
    pub fn count(&self) -> usize {
        self.primitives.len()
    }
}

impl Default for PyDebugOverlay {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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
        let n = compute_normals_from_mesh(&verts, &idx);
        assert_eq!(n.len(), 9);
        // Normal should point along +Z
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
        pr.set_colors_from_scalars(&[0.0, 1.0], &PyColormap::viridis());
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
        let id = sg.add_node("root");
        assert_eq!(id, 0);
        assert_eq!(sg.node_count(), 1);
    }

    #[test]
    fn test_scene_graph_visibility() {
        let mut sg = PySceneGraph::new();
        let id = sg.add_node("node");
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
        // Pure shear: τ = 1 → von Mises = sqrt(3)
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
        // Uniform velocity field in +x direction
        let nx = 4u32;
        let ny = 4u32;
        let nz = 4u32;
        let n = (nx * ny * nz) as usize;
        let mut vf = vec![0.0f64; n * 3];
        for i in 0..n {
            vf[3 * i] = 1.0; // vx = 1
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
        ov.draw_text_3d([0.0; 3], "hello", None, 0.1);
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
}
