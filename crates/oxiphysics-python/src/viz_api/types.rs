//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

use super::functions::{
    coolwarm_sample, hot_sample, inferno_sample, jet_sample, magma_sample, plasma_sample,
    rdbu_sample, spectral_sample, turbo_sample, viridis_sample,
};

/// A colormap that maps scalar values in \[0, 1\] to RGBA colors.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyColormap {
    /// Create a new colormap with explicit range.
    #[new]
    pub fn new(kind: ColormapKind, vmin: f64, vmax: f64) -> Self {
        Self {
            kind,
            vmin,
            vmax,
            alpha: 1.0,
        }
    }
    /// Create a Viridis colormap over \[0, 1\].
    #[staticmethod]
    pub fn viridis() -> Self {
        Self::new(ColormapKind::Viridis, 0.0, 1.0)
    }
    /// Create a Plasma colormap over \[0, 1\].
    #[staticmethod]
    pub fn plasma() -> Self {
        Self::new(ColormapKind::Plasma, 0.0, 1.0)
    }
    /// Create a Jet colormap over \[0, 1\].
    #[staticmethod]
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
/// A 3×3 symmetric stress tensor stored as \[s11, s22, s33, s12, s13, s23\].
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressTensor {
    /// Components \[σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃\].
    pub components: [f64; 6],
}
#[pymethods]
impl StressTensor {
    /// Create a new stress tensor.
    #[new]
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
/// Debug overlay manager for visualizing physics primitives.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyDebugOverlay {
    /// Create a new empty debug overlay.
    #[new]
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
        text: String,
        color: Option<[f64; 4]>,
        size: f64,
    ) {
        let c = color.unwrap_or(self.default_color);
        self.primitives.push(DebugPrimitive::Text3D {
            position,
            text,
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
// Default impl is in trait_impls.rs to avoid duplication
/// Result of streamline tracing.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Streamline {
    /// Flat array of 3D positions along the streamline.
    pub points: Vec<f64>,
    /// Total arc length.
    pub arc_length: f64,
}
/// A control point in a transfer function (scalar → RGBA).
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransferPoint {
    /// Scalar value in \[0, 1\].
    pub scalar: f64,
    /// RGBA color at this point.
    pub color: [f64; 4],
}
/// Output from the stress visualizer for a single tensor.
#[pyclass(from_py_object)]
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
/// Named colormap for scalar field visualization.
#[pyclass(eq, eq_int, from_py_object)]
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
/// Post-processing pipeline configuration.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyPostProcessor {
    /// Create a default post-processing pipeline.
    #[new]
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
    #[staticmethod]
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
// Default impl is in trait_impls.rs to avoid duplication
/// Transfer function mapping scalar densities to RGBA for volume rendering.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTransferFunction {
    /// Sorted control points.
    pub control_points: Vec<TransferPoint>,
}
#[pymethods]
impl PyTransferFunction {
    /// Create a simple two-point transfer function.
    #[staticmethod]
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
/// Stress field visualizer providing principal direction glyphs and Mohr circle data.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyStressVisualizer {
    /// Input stress tensor field (one per element/node).
    pub tensors: Vec<StressTensor>,
    /// Colormap for von Mises coloring.
    pub colormap: PyColormap,
    /// Scale factor for principal direction glyphs.
    pub glyph_scale: f64,
}
#[pymethods]
impl PyStressVisualizer {
    /// Create a new stress visualizer.
    #[new]
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
// Default impl is in trait_impls.rs to avoid duplication
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
// Default impl is in trait_impls.rs to avoid duplication
/// Surface material properties for mesh rendering.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyMaterial {
    /// Create a default PBR material.
    #[new]
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
    #[staticmethod]
    pub fn red() -> Self {
        Self::new([1.0, 0.0, 0.0, 1.0])
    }
    /// Blue material.
    #[staticmethod]
    pub fn blue() -> Self {
        Self::new([0.0, 0.3, 1.0, 1.0])
    }
    /// Metallic silver material.
    #[staticmethod]
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
/// Renderer for particle systems using billboard or instanced geometry.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyParticleRenderer {
    /// Create a new particle renderer.
    #[new]
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
    pub fn set_colors_from_scalars(&mut self, scalars: Vec<f64>, colormap: PyColormap) {
        let n = self.positions.len() / 3;
        self.colors.clear();
        for i in 0..n {
            let s = scalars.get(i).copied().unwrap_or(0.0);
            let rgba = colormap.map_scalar(s);
            self.colors.extend_from_slice(&rgba);
        }
        self.colormap = Some(colormap);
    }
    /// Return the number of particles.
    pub fn particle_count(&self) -> usize {
        self.positions.len() / 3
    }
}
/// 3D camera with orbit, pan, zoom controls.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyCamera {
    /// Create a default perspective camera.
    #[new]
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
    #[staticmethod]
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
/// Volume renderer for 3D scalar fields.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyVolumeRenderer {
    /// Create a volume renderer from a flat scalar field.
    #[new]
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
/// Tone mapping operator selection.
#[pyclass(eq, eq_int, from_py_object)]
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
/// Streamline tracer using RK4 integration through a 3D velocity field.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyStreamlineTracer {
    /// Create a new streamline tracer.
    #[new]
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
    pub fn trace_multiple(&self, seeds: Vec<[f64; 3]>) -> Vec<Streamline> {
        seeds.into_iter().map(|s| self.trace(s)).collect()
    }
}
/// A light source in the scene.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PyLight {
    /// Create a point light.
    #[staticmethod]
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
    #[staticmethod]
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
// Default impl is in trait_impls.rs to avoid duplication
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
/// A transform node in the scene hierarchy.
#[pyclass(from_py_object)]
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
#[pymethods]
impl PySceneNode {
    /// Create an identity transform node.
    #[staticmethod]
    pub fn identity(id: u32, label: String) -> Self {
        Self {
            id,
            parent: None,
            translation: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0; 3],
            label,
            visible: true,
        }
    }
}
