//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::*;
use super::functions::{Tensor3, Vec3f};

/// Classification of a point in a tensor field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyPointType {
    /// Regular point — eigenvalues are all distinct.
    Regular,
    /// Isotropic point — all three eigenvalues equal (λ₁ = λ₂ = λ₃).
    Isotropic,
    /// Planar degenerate — λ₂ = λ₃ (disk-like degeneracy).
    PlanarDegenerate,
    /// Linear degenerate — λ₁ = λ₂ (needle-like degeneracy).
    LinearDegenerate,
    /// Umbilical point — all eigenvalues equal (same as isotropic in 3-D).
    Umbilical,
    /// Saddle-like point in 2-D slice.
    Saddle,
    /// Node-like point in 2-D slice.
    Node,
    /// Center / spiral point in 2-D slice.
    Center,
}
/// A single diffusion tensor imaging voxel with full DTI metrics.
#[derive(Debug, Clone)]
pub struct DtiVoxel {
    /// Diffusion tensor in Voigt notation.
    pub tensor: Tensor3,
    /// Eigenvalues λ₁ ≥ λ₂ ≥ λ₃.
    pub eigenvalues: [f64; 3],
    /// Major eigenvector (principal diffusion direction).
    pub principal_direction: Vec3f,
    /// Fractional anisotropy ∈ \[0, 1\].
    pub fa: f64,
    /// Mean diffusivity.
    pub md: f64,
    /// Axial diffusivity.
    pub ad: f64,
    /// Radial diffusivity.
    pub rd: f64,
    /// Relative anisotropy.
    pub ra: f64,
    /// Volume ratio.
    pub vr: f64,
}
impl DtiVoxel {
    /// Compute all DTI metrics from a diffusion tensor.
    pub fn from_tensor(tensor: Tensor3) -> Self {
        let eig = eig3(&tensor);
        let fa = fractional_anisotropy(&eig.values);
        let md = mean_diffusivity(&eig.values);
        let ad = axial_diffusivity(&eig.values);
        let rd = radial_diffusivity(&eig.values);
        let ra = relative_anisotropy(&eig.values);
        let vr = volume_ratio(&eig.values);
        Self {
            tensor,
            eigenvalues: eig.values,
            principal_direction: eig.vectors[0],
            fa,
            md,
            ad,
            rd,
            ra,
            vr,
        }
    }
}
/// Integrates hyperstreamlines through a tensor field following the major eigenvector.
///
/// Uses 4th-order Runge-Kutta integration with adaptive step-size control.
#[derive(Debug, Clone)]
pub struct HyperstreamlineTracer {
    /// Integration step size in spatial units.
    pub step_size: f64,
    /// Maximum number of integration steps.
    pub max_steps: usize,
    /// Minimum FA threshold — stop integrating below this value.
    pub min_fa: f64,
    /// Minimum angle between consecutive steps (radians) — prevents sharp turns.
    pub max_curvature_angle: f64,
    /// Whether to trace in both forward and backward directions.
    pub bidirectional: bool,
}
impl HyperstreamlineTracer {
    /// Create a tracer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Trace a hyperstreamline from the given seed point.
    ///
    /// `field_fn` maps a 3-D position to the tensor at that point (sampled from
    /// a grid or analytical field).
    pub fn trace<F>(&self, seed: Vec3f, field_fn: &F) -> Hyperstreamline
    where
        F: Fn(Vec3f) -> Tensor3,
    {
        let forward = self.trace_direction(seed, field_fn, 1.0);
        if self.bidirectional {
            let backward = self.trace_direction(seed, field_fn, -1.0);
            self.merge(backward, forward)
        } else {
            forward
        }
    }
    /// Trace multiple hyperstreamlines from a list of seeds.
    pub fn trace_many<F>(&self, seeds: &[Vec3f], field_fn: &F) -> Vec<Hyperstreamline>
    where
        F: Fn(Vec3f) -> Tensor3,
    {
        seeds.iter().map(|&s| self.trace(s, field_fn)).collect()
    }
    fn trace_direction<F>(&self, seed: Vec3f, field_fn: &F, dir_sign: f64) -> Hyperstreamline
    where
        F: Fn(Vec3f) -> Tensor3,
    {
        let mut points = vec![seed];
        let mut eigenvalues_list: Vec<[f64; 3]> = Vec::new();
        let mut fa_list = Vec::new();
        let mut md_list = Vec::new();
        let mut pos = seed;
        let t0 = field_fn(seed);
        let eig0 = eig3(&t0);
        let mut prev_dir = scale3(&eig0.major_vec(), dir_sign);
        let fa0 = fractional_anisotropy(&eig0.values);
        eigenvalues_list.push(eig0.values);
        fa_list.push(fa0);
        md_list.push(mean_diffusivity(&eig0.values));
        for _ in 0..self.max_steps {
            let tensor = field_fn(pos);
            let eig = eig3(&tensor);
            let fa = fractional_anisotropy(&eig.values);
            if fa < self.min_fa {
                break;
            }
            let raw_dir = eig.major_vec();
            let aligned = if dot3(&raw_dir, &prev_dir) >= 0.0 {
                raw_dir
            } else {
                scale3(&raw_dir, -1.0)
            };
            let dir = normalize3(&aligned);
            let angle = dot3(&dir, &normalize3(&prev_dir)).clamp(-1.0, 1.0).acos();
            if angle > self.max_curvature_angle {
                break;
            }
            let k1 = dir;
            let p1 = add3(&pos, &scale3(&k1, self.step_size * 0.5));
            let t1 = field_fn(p1);
            let e1 = eig3(&t1);
            let d1_raw = e1.major_vec();
            let d1 = if dot3(&d1_raw, &k1) >= 0.0 {
                d1_raw
            } else {
                scale3(&d1_raw, -1.0)
            };
            let k2 = normalize3(&d1);
            let p2 = add3(&pos, &scale3(&k2, self.step_size * 0.5));
            let t2 = field_fn(p2);
            let e2 = eig3(&t2);
            let d2_raw = e2.major_vec();
            let d2 = if dot3(&d2_raw, &k2) >= 0.0 {
                d2_raw
            } else {
                scale3(&d2_raw, -1.0)
            };
            let k3 = normalize3(&d2);
            let p3 = add3(&pos, &scale3(&k3, self.step_size));
            let t3 = field_fn(p3);
            let e3 = eig3(&t3);
            let d3_raw = e3.major_vec();
            let d3 = if dot3(&d3_raw, &k3) >= 0.0 {
                d3_raw
            } else {
                scale3(&d3_raw, -1.0)
            };
            let k4 = normalize3(&d3);
            let rk4_dir = normalize3(&[
                (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
                (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
                (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
            ]);
            pos = add3(&pos, &scale3(&rk4_dir, self.step_size));
            prev_dir = rk4_dir;
            points.push(pos);
            eigenvalues_list.push(eig.values);
            fa_list.push(fa);
            md_list.push(mean_diffusivity(&eig.values));
        }
        Hyperstreamline {
            points,
            eigenvalues: eigenvalues_list,
            fa: fa_list,
            md: md_list,
        }
    }
    fn merge(
        &self,
        mut backward: Hyperstreamline,
        mut forward: Hyperstreamline,
    ) -> Hyperstreamline {
        backward.points.reverse();
        backward.eigenvalues.reverse();
        backward.fa.reverse();
        backward.md.reverse();
        if !forward.points.is_empty() {
            forward.points.remove(0);
            if !forward.eigenvalues.is_empty() {
                forward.eigenvalues.remove(0);
            }
            if !forward.fa.is_empty() {
                forward.fa.remove(0);
            }
            if !forward.md.is_empty() {
                forward.md.remove(0);
            }
        }
        backward.points.extend(forward.points);
        backward.eigenvalues.extend(forward.eigenvalues);
        backward.fa.extend(forward.fa);
        backward.md.extend(forward.md);
        backward
    }
}
/// Summary statistics over a DTI volume.
#[derive(Debug, Clone, Default)]
pub struct DtiStatistics {
    /// Mean fractional anisotropy.
    pub mean_fa: f64,
    /// Mean mean diffusivity.
    pub mean_md: f64,
    /// Mean axial diffusivity.
    pub mean_ad: f64,
    /// Mean radial diffusivity.
    pub mean_rd: f64,
    /// Maximum FA in the volume.
    pub max_fa: f64,
    /// Minimum FA in the volume.
    pub min_fa: f64,
    /// Number of voxels processed.
    pub n_voxels: usize,
}
/// Log-Euclidean tensor interpolation on the SPD (symmetric positive-definite) manifold.
///
/// Performs interpolation in the matrix logarithm space, which ensures the
/// result remains on the SPD manifold (no negative eigenvalues).
#[derive(Debug, Clone)]
pub struct TensorInterpolation;
impl TensorInterpolation {
    /// Compute the matrix logarithm of an SPD tensor.
    ///
    /// Uses the spectral decomposition: log(T) = V log(Λ) Vᵀ.
    pub fn log(t: &Tensor3) -> Tensor3 {
        let eig = eig3(t);
        let log_vals: [f64; 3] = std::array::from_fn(|i| eig.values[i].max(1e-15).ln());
        Self::recompose(&eig.vectors, &log_vals)
    }
    /// Compute the matrix exponential of a symmetric tensor.
    ///
    /// Uses the spectral decomposition: exp(S) = V exp(Λ) Vᵀ.
    pub fn exp(s: &Tensor3) -> Tensor3 {
        let eig = eig3(s);
        let exp_vals: [f64; 3] = std::array::from_fn(|i| eig.values[i].exp());
        Self::recompose(&eig.vectors, &exp_vals)
    }
    /// Linearly interpolate two SPD tensors in log-Euclidean space.
    ///
    /// `t = (1 - alpha) * T0 + alpha * T1` in the log domain.
    pub fn lerp(t0: &Tensor3, t1: &Tensor3, alpha: f64) -> Tensor3 {
        let s0 = Self::log(t0);
        let s1 = Self::log(t1);
        let s = voigt_lerp(&s0, &s1, alpha);
        Self::exp(&s)
    }
    /// Geodesic interpolation on the SPD manifold at parameter `alpha ∈ [0, 1]`.
    ///
    /// This is equivalent to log-Euclidean interpolation for this implementation.
    pub fn geodesic(t0: &Tensor3, t1: &Tensor3, alpha: f64) -> Tensor3 {
        Self::lerp(t0, t1, alpha)
    }
    /// Weighted Fréchet mean of a list of SPD tensors with given weights.
    ///
    /// Weights need not be normalized — they are normalized internally.
    pub fn weighted_mean(tensors: &[Tensor3], weights: &[f64]) -> Tensor3 {
        assert_eq!(tensors.len(), weights.len());
        let weight_sum: f64 = weights.iter().sum();
        let mut mean_log = [0.0_f64; 6];
        for (t, &w) in tensors.iter().zip(weights.iter()) {
            let lt = Self::log(t);
            for k in 0..6 {
                mean_log[k] += (w / weight_sum) * lt[k];
            }
        }
        Self::exp(&mean_log)
    }
    /// Compute the Riemannian distance between two SPD tensors (log-Euclidean metric).
    pub fn distance(t0: &Tensor3, t1: &Tensor3) -> f64 {
        let s0 = Self::log(t0);
        let s1 = Self::log(t1);
        let diff: [f64; 6] = std::array::from_fn(|i| s0[i] - s1[i]);
        let m = voigt_to_mat3(&diff);
        m.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    fn recompose(vecs: &[Vec3f; 3], vals: &[f64; 3]) -> Tensor3 {
        let mut m = [0.0_f64; 9];
        for k in 0..3 {
            let v = &vecs[k];
            for i in 0..3 {
                for j in 0..3 {
                    m[i * 3 + j] += vals[k] * v[i] * v[j];
                }
            }
        }
        mat3_to_voigt(&m)
    }
}
/// Analyzes tensor field topology to find degenerate lines, umbilical points,
/// and classify the topological structure of the field.
#[derive(Debug, Clone)]
pub struct TensorTopology {
    /// Tolerance for considering two eigenvalues equal (degenerate).
    pub degeneracy_tolerance: f64,
    /// Whether to search for 2-D degenerate features in XY slices.
    pub analyze_2d_slices: bool,
}
impl TensorTopology {
    /// Create a topology analyzer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Classify a single tensor point.
    pub fn classify_point(&self, tensor: &Tensor3, position: Vec3f) -> TopologyPoint {
        let eig = eig3(tensor);
        let [l1, l2, l3] = eig.values;
        let tol = self.degeneracy_tolerance;
        let kind = if (l1 - l2).abs() < tol && (l2 - l3).abs() < tol {
            TopologyPointType::Isotropic
        } else if (l2 - l3).abs() < tol {
            TopologyPointType::PlanarDegenerate
        } else if (l1 - l2).abs() < tol {
            TopologyPointType::LinearDegenerate
        } else {
            TopologyPointType::Regular
        };
        let fa = fractional_anisotropy(&eig.values);
        TopologyPoint {
            position,
            kind,
            eigenvalues: eig.values,
            anisotropy: fa,
        }
    }
    /// Classify all points in a tensor field.
    pub fn classify_field(&self, tensors: &[Tensor3], positions: &[Vec3f]) -> Vec<TopologyPoint> {
        assert_eq!(tensors.len(), positions.len());
        tensors
            .iter()
            .zip(positions.iter())
            .map(|(t, &p)| self.classify_point(t, p))
            .collect()
    }
    /// Find degenerate points in a tensor field (where two eigenvalues are nearly equal).
    pub fn find_degenerate_points<'a>(
        &self,
        points: &'a [TopologyPoint],
    ) -> Vec<&'a TopologyPoint> {
        points
            .iter()
            .filter(|p| p.kind != TopologyPointType::Regular)
            .collect()
    }
    /// Compute the degenerate line measure (how close two eigenvalues are).
    ///
    /// Returns 0 for fully isotropic, 1 for maximally anisotropic (no degeneracy).
    pub fn degeneracy_measure(&self, tensor: &Tensor3) -> f64 {
        let eig = eig3(tensor);
        let [l1, l2, l3] = eig.values;
        let denom = l1.abs() + l3.abs();
        if denom < 1e-15 {
            return 0.0;
        }
        let d12 = (l1 - l2).abs() / denom;
        let d23 = (l2 - l3).abs() / denom;
        d12.min(d23)
    }
    /// Compute the 2-D topology index (Poincaré index) for a tensor field in the XY plane.
    ///
    /// Samples the tensor field on a square contour and counts the winding number
    /// of the angle field.
    pub fn poincare_index_2d<F>(
        &self,
        center: Vec3f,
        radius: f64,
        field_fn: &F,
        n_samples: usize,
    ) -> f64
    where
        F: Fn(Vec3f) -> Tensor3,
    {
        let n = n_samples.max(8);
        let mut angles = Vec::with_capacity(n);
        for i in 0..n {
            let theta = 2.0 * PI * i as f64 / n as f64;
            let p = [
                center[0] + radius * theta.cos(),
                center[1] + radius * theta.sin(),
                center[2],
            ];
            let t = field_fn(p);
            let eig = eig3(&t);
            let v = eig.vectors[0];
            angles.push(v[1].atan2(v[0]));
        }
        let mut total = 0.0_f64;
        for i in 0..n {
            let next = (i + 1) % n;
            let mut da = angles[next] - angles[i];
            while da > PI {
                da -= 2.0 * PI;
            }
            while da < -PI {
                da += 2.0 * PI;
            }
            total += da;
        }
        total / PI
    }
    /// Classify a 2-D tensor topology point using eigenvalue signs and ratios.
    pub fn classify_2d(&self, tensor: &Tensor3) -> TopologyPointType {
        let [l1, l2, _l3] = eig3(tensor).values;
        let tol = self.degeneracy_tolerance;
        if (l1 - l2).abs() < tol {
            return TopologyPointType::Umbilical;
        }
        if l1 * l2 < 0.0 {
            TopologyPointType::Saddle
        } else if l1 * l2 > 0.0 {
            TopologyPointType::Node
        } else {
            TopologyPointType::Center
        }
    }
}
/// A topologically significant point in a tensor field.
#[derive(Debug, Clone)]
pub struct TopologyPoint {
    /// Position of this topological feature.
    pub position: Vec3f,
    /// Classification of this point.
    pub kind: TopologyPointType,
    /// Eigenvalues at this point.
    pub eigenvalues: [f64; 3],
    /// Normalized anisotropy measure ∈ \[0, 1\].
    pub anisotropy: f64,
}
/// Controls how [`TensorGlyph`]s are generated from a tensor field.
#[derive(Debug, Clone)]
pub struct TensorGlyphBuilder {
    /// Uniform scale factor applied to all semi-axes.
    pub scale: f64,
    /// Minimum semi-axis length (clamp, prevents degenerate glyphs).
    pub min_axis: f64,
    /// Maximum semi-axis length (clamp).
    pub max_axis: f64,
    /// Whether to use absolute eigenvalue magnitudes for axis lengths.
    pub abs_eigenvalues: bool,
    /// Color for positive eigenvalue direction.
    pub positive_color: [f32; 4],
    /// Color for negative eigenvalue direction.
    pub negative_color: [f32; 4],
}
impl TensorGlyphBuilder {
    /// Create a builder with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Build a single ellipsoid glyph from a tensor at the given position.
    pub fn build(&self, tensor: &Tensor3, position: Vec3f) -> TensorGlyph {
        let eig = eig3(tensor);
        let semi_axes: [f64; 3] = std::array::from_fn(|i| {
            let raw = if self.abs_eigenvalues {
                eig.values[i].abs()
            } else {
                eig.values[i].max(0.0)
            };
            (raw * self.scale).clamp(self.min_axis, self.max_axis)
        });
        let color = if eig.values[0] >= 0.0 {
            self.positive_color
        } else {
            self.negative_color
        };
        let scalar = fractional_anisotropy(&eig.values);
        TensorGlyph {
            position,
            semi_axes,
            eigenvectors: eig.vectors,
            eigenvalues: eig.values,
            color,
            scalar,
        }
    }
    /// Build glyphs for a list of (tensor, position) pairs.
    pub fn build_field(&self, tensors: &[(Tensor3, Vec3f)]) -> Vec<TensorGlyph> {
        tensors.iter().map(|(t, p)| self.build(t, *p)).collect()
    }
    /// Build a glyph for every point in a regular 3-D grid.
    pub fn build_grid(
        &self,
        field: &[Tensor3],
        nx: usize,
        ny: usize,
        nz: usize,
        origin: Vec3f,
        spacing: Vec3f,
    ) -> Vec<TensorGlyph> {
        assert_eq!(field.len(), nx * ny * nz);
        let mut glyphs = Vec::with_capacity(field.len());
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = iz * ny * nx + iy * nx + ix;
                    let pos = [
                        origin[0] + ix as f64 * spacing[0],
                        origin[1] + iy as f64 * spacing[1],
                        origin[2] + iz as f64 * spacing[2],
                    ];
                    glyphs.push(self.build(&field[idx], pos));
                }
            }
        }
        glyphs
    }
}
/// Stress tensor field analysis including principal trajectories and decomposition.
#[derive(Debug, Clone)]
pub struct StressTensorField {
    /// Scale factor for visualizing principal stress vectors.
    pub glyph_scale: f64,
    /// Color for tensile (positive) principal stresses.
    pub tension_color: [f32; 4],
    /// Color for compressive (negative) principal stresses.
    pub compression_color: [f32; 4],
    /// Whether to normalize eigenvectors to unit length (disable to show magnitude).
    pub normalize_glyphs: bool,
}
impl StressTensorField {
    /// Create a stress field analyzer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Decompose a stress tensor into principal stresses and hydrostatic/deviatoric parts.
    pub fn decompose(&self, tensor: &Tensor3) -> PrincipalStress {
        let eig = eig3(tensor);
        let hydrostatic = (eig.values[0] + eig.values[1] + eig.values[2]) / 3.0;
        let mut dev = *tensor;
        dev[0] -= hydrostatic;
        dev[1] -= hydrostatic;
        dev[2] -= hydrostatic;
        let von_mises = von_mises_stress(tensor);
        let tresca = eig.values[0] - eig.values[2];
        let mean_normal = hydrostatic;
        PrincipalStress {
            values: eig.values,
            directions: eig.vectors,
            hydrostatic,
            deviatoric: dev,
            von_mises,
            tresca,
            mean_normal,
        }
    }
    /// Generate line segments for the principal stress trajectories at a grid of points.
    ///
    /// Returns pairs of `(start, end)` points for each principal direction.
    pub fn principal_trajectories(
        &self,
        tensors: &[Tensor3],
        positions: &[Vec3f],
    ) -> Vec<(Vec3f, Vec3f, [f32; 4])> {
        assert_eq!(tensors.len(), positions.len());
        let mut lines = Vec::with_capacity(tensors.len() * 3);
        for (tensor, &pos) in tensors.iter().zip(positions.iter()) {
            let ps = self.decompose(tensor);
            for k in 0..3 {
                let dir = ps.directions[k];
                let mag = if self.normalize_glyphs {
                    1.0
                } else {
                    ps.values[k].abs()
                };
                let half = scale3(&dir, 0.5 * mag * self.glyph_scale);
                let start = sub3(&pos, &half);
                let end = add3(&pos, &half);
                let color = if ps.values[k] >= 0.0 {
                    self.tension_color
                } else {
                    self.compression_color
                };
                lines.push((start, end, color));
            }
        }
        lines
    }
    /// Compute von Mises stress for a list of tensors.
    pub fn von_mises_field(&self, tensors: &[Tensor3]) -> Vec<f64> {
        tensors.iter().map(von_mises_stress).collect()
    }
    /// Compute Tresca stress (maximum shear) for a list of tensors.
    pub fn tresca_field(&self, tensors: &[Tensor3]) -> Vec<f64> {
        tensors
            .iter()
            .map(|t| {
                let eig = eig3(t);
                eig.values[0] - eig.values[2]
            })
            .collect()
    }
    /// Compute the hydrostatic stress for a list of tensors.
    pub fn hydrostatic_field(&self, tensors: &[Tensor3]) -> Vec<f64> {
        tensors.iter().map(|t| (t[0] + t[1] + t[2]) / 3.0).collect()
    }
}
/// An ellipsoid glyph representing a symmetric tensor at a spatial location.
///
/// The three semi-axes of the ellipsoid correspond to the three eigenvectors,
/// with lengths proportional to the absolute eigenvalues (optionally clamped).
#[derive(Debug, Clone)]
pub struct TensorGlyph {
    /// Center of the glyph.
    pub position: Vec3f,
    /// Semi-axis lengths \[a₁, a₂, a₃\] (non-negative, eigenvalue magnitudes × scale).
    pub semi_axes: [f64; 3],
    /// Unit eigenvectors forming the local frame \[e₁, e₂, e₃\].
    pub eigenvectors: [Vec3f; 3],
    /// The three eigenvalues λ₁ ≥ λ₂ ≥ λ₃.
    pub eigenvalues: [f64; 3],
    /// Color of the glyph (RGBA in \[0, 1\]).
    pub color: [f32; 4],
    /// Scalar label for coloring (e.g., FA for DTI glyphs).
    pub scalar: f64,
}
/// Diffusion tensor imaging analysis and color-coded FA map generation.
#[derive(Debug, Clone)]
pub struct DiffusionTensorImaging {
    /// Whether to clamp FA values to \[0, 1\] before coloring.
    pub clamp_fa: bool,
    /// Minimum FA for fiber tracking (voxels below this are not tracked).
    pub tracking_fa_threshold: f64,
    /// Color mode for FA maps.
    pub color_mode: FaColorMode,
}
impl DiffusionTensorImaging {
    /// Create a DTI processor with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Process a volume of diffusion tensors and return per-voxel DTI metrics.
    pub fn process(&self, tensors: &[Tensor3]) -> Vec<DtiVoxel> {
        tensors.iter().map(|t| DtiVoxel::from_tensor(*t)).collect()
    }
    /// Generate a color-coded FA map as RGBA pixels.
    ///
    /// Each pixel is colored by FA × |principal direction| (directional encoding),
    /// scaled to \[0, 255\] as `u8`.
    pub fn fa_color_map(&self, voxels: &[DtiVoxel]) -> Vec<[u8; 4]> {
        voxels.iter().map(|v| self.voxel_color(v)).collect()
    }
    /// Generate a grayscale FA map as `u8` intensity values.
    pub fn fa_grayscale_map(&self, voxels: &[DtiVoxel]) -> Vec<u8> {
        voxels
            .iter()
            .map(|v| {
                let fa = if self.clamp_fa {
                    v.fa.clamp(0.0, 1.0)
                } else {
                    v.fa
                };
                (fa * 255.0) as u8
            })
            .collect()
    }
    /// Identify voxels suitable for fiber tracking (FA above threshold).
    pub fn tracking_mask(&self, voxels: &[DtiVoxel]) -> Vec<bool> {
        voxels
            .iter()
            .map(|v| v.fa >= self.tracking_fa_threshold)
            .collect()
    }
    /// Compute summary statistics over the DTI volume.
    pub fn statistics(&self, voxels: &[DtiVoxel]) -> DtiStatistics {
        if voxels.is_empty() {
            return DtiStatistics::default();
        }
        let n = voxels.len() as f64;
        let mean_fa = voxels.iter().map(|v| v.fa).sum::<f64>() / n;
        let mean_md = voxels.iter().map(|v| v.md).sum::<f64>() / n;
        let mean_ad = voxels.iter().map(|v| v.ad).sum::<f64>() / n;
        let mean_rd = voxels.iter().map(|v| v.rd).sum::<f64>() / n;
        let max_fa = voxels
            .iter()
            .map(|v| v.fa)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_fa = voxels.iter().map(|v| v.fa).fold(f64::INFINITY, f64::min);
        DtiStatistics {
            mean_fa,
            mean_md,
            mean_ad,
            mean_rd,
            max_fa,
            min_fa,
            n_voxels: voxels.len(),
        }
    }
    fn voxel_color(&self, v: &DtiVoxel) -> [u8; 4] {
        let fa = if self.clamp_fa {
            v.fa.clamp(0.0, 1.0)
        } else {
            v.fa
        };
        match self.color_mode {
            FaColorMode::Directional => {
                let pd = v.principal_direction;
                let r = (pd[0].abs() * fa * 255.0) as u8;
                let g = (pd[1].abs() * fa * 255.0) as u8;
                let b = (pd[2].abs() * fa * 255.0) as u8;
                [r, g, b, 255]
            }
            FaColorMode::Grayscale => {
                let g = (fa * 255.0) as u8;
                [g, g, g, 255]
            }
            FaColorMode::Spectral => {
                let t = fa.clamp(0.0, 1.0);
                let (r, g, b) = spectral_color(t);
                [r, g, b, 255]
            }
        }
    }
}
/// Result of eigendecomposition: eigenvalues and eigenvectors.
#[derive(Debug, Clone)]
pub struct Eigen3 {
    /// Eigenvalues in descending order: λ₁ ≥ λ₂ ≥ λ₃.
    pub values: [f64; 3],
    /// Corresponding unit eigenvectors (columns of V in M = V Λ Vᵀ).
    pub vectors: [Vec3f; 3],
}
impl Eigen3 {
    /// Return the major (largest) eigenvalue.
    pub fn major(&self) -> f64 {
        self.values[0]
    }
    /// Return the minor (smallest) eigenvalue.
    pub fn minor(&self) -> f64 {
        self.values[2]
    }
    /// Return the major eigenvector.
    pub fn major_vec(&self) -> Vec3f {
        self.vectors[0]
    }
    /// Return the minor eigenvector.
    pub fn minor_vec(&self) -> Vec3f {
        self.vectors[2]
    }
    /// Return the intermediate eigenvalue.
    pub fn mid(&self) -> f64 {
        self.values[1]
    }
    /// Return the intermediate eigenvector.
    pub fn mid_vec(&self) -> Vec3f {
        self.vectors[1]
    }
}
/// Color mode for FA maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaColorMode {
    /// Standard color-coded FA: |RGB| = FA × |principal direction|.
    Directional,
    /// Grayscale: intensity proportional to FA.
    Grayscale,
    /// Spectral colormap (blue=low, red=high).
    Spectral,
}
/// Result of principal stress decomposition at a point.
#[derive(Debug, Clone)]
pub struct PrincipalStress {
    /// Principal stress values σ₁ ≥ σ₂ ≥ σ₃.
    pub values: [f64; 3],
    /// Principal stress directions (eigenvectors).
    pub directions: [Vec3f; 3],
    /// Hydrostatic (volumetric) stress = trace/3.
    pub hydrostatic: f64,
    /// Deviatoric stress tensor (Voigt notation).
    pub deviatoric: Tensor3,
    /// Von Mises stress.
    pub von_mises: f64,
    /// Tresca stress (max shear).
    pub tresca: f64,
    /// Mean normal stress.
    pub mean_normal: f64,
}
/// A single hyperstreamline path through a tensor field.
#[derive(Debug, Clone)]
pub struct Hyperstreamline {
    /// 3-D points along the hyperstreamline.
    pub points: Vec<Vec3f>,
    /// Eigenvalues at each integration point.
    pub eigenvalues: Vec<[f64; 3]>,
    /// Fractional anisotropy at each integration point.
    pub fa: Vec<f64>,
    /// Mean diffusivity at each integration point.
    pub md: Vec<f64>,
}
