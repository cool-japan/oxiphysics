// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural rendering and implicit representations for scientific visualization.
//!
//! Implements NeRF (Neural Radiance Fields), Gaussian Splatting, Neural SDFs,
//! and related neural rendering techniques for advanced visualization pipelines.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers (no nalgebra — f64 arrays only)
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length of a 3-vector.
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Normalize a 3-vector; returns the zero vector if near-zero length.
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0; 3]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

/// Add two 3-vectors.
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Clamp a value to \[lo, hi\].
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

/// Sigmoid activation: 1 / (1 + exp(-x)).
fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// ReLU activation: max(0, x).
fn relu(x: f64) -> f64 {
    x.max(0.0)
}

/// Softplus activation: log(1 + exp(x)).
fn softplus(x: f64) -> f64 {
    (1.0 + x.exp()).ln()
}

// ---------------------------------------------------------------------------
// PositionalEncoding
// ---------------------------------------------------------------------------

/// Fourier positional encoding for neural radiance fields.
///
/// Maps a coordinate `x` to a higher-dimensional space using sinusoidal
/// frequency bands, enabling neural networks to represent high-frequency
/// details: `[sin(2^0*π*x), cos(2^0*π*x), ..., sin(2^L*π*x), cos(2^L*π*x)]`.
#[derive(Debug, Clone)]
pub struct PositionalEncoding {
    /// Number of frequency bands L.
    pub num_frequencies: usize,
    /// Whether to include the original input coordinates in the output.
    pub include_input: bool,
    /// Pre-computed frequency multipliers: `[1, 2, 4, ..., 2^(L-1)]`.
    pub frequencies: Vec<f64>,
}

impl PositionalEncoding {
    /// Create a new positional encoding with `num_frequencies` frequency bands.
    pub fn new(num_frequencies: usize, include_input: bool) -> Self {
        let frequencies = (0..num_frequencies)
            .map(|k| (2.0_f64).powi(k as i32))
            .collect();
        Self {
            num_frequencies,
            include_input,
            frequencies,
        }
    }

    /// Encode a scalar value `x` into a vector of sinusoidal features.
    ///
    /// Output dimensionality: `2 * num_frequencies` (+ 1 if `include_input`).
    pub fn encode_scalar(&self, x: f64) -> Vec<f64> {
        let mut out = Vec::new();
        if self.include_input {
            out.push(x);
        }
        for &freq in &self.frequencies {
            out.push((freq * PI * x).sin());
            out.push((freq * PI * x).cos());
        }
        out
    }

    /// Encode a 3D point `p` into concatenated per-dimension sinusoidal features.
    ///
    /// Output dimensionality: `3 * 2 * num_frequencies` (+ 3 if `include_input`).
    pub fn encode_point(&self, p: [f64; 3]) -> Vec<f64> {
        let mut out = Vec::new();
        for &coord in &p {
            out.extend(self.encode_scalar(coord));
        }
        out
    }

    /// Encode a viewing direction (unit 3-vector) into sinusoidal features.
    pub fn encode_direction(&self, d: [f64; 3]) -> Vec<f64> {
        self.encode_point(d)
    }

    /// Output dimensionality for a 3D point encoding.
    pub fn point_output_dim(&self) -> usize {
        3 * (2 * self.num_frequencies + if self.include_input { 1 } else { 0 })
    }

    /// Output dimensionality for a direction encoding.
    pub fn direction_output_dim(&self) -> usize {
        self.point_output_dim()
    }
}

// ---------------------------------------------------------------------------
// MlpLayer
// ---------------------------------------------------------------------------

/// A single fully-connected layer with optional activation.
#[derive(Debug, Clone)]
pub struct MlpLayer {
    /// Weight matrix stored row-major: `weights[i * in_dim + j]`.
    pub weights: Vec<f64>,
    /// Bias vector of length `out_dim`.
    pub biases: Vec<f64>,
    /// Input dimensionality.
    pub in_dim: usize,
    /// Output dimensionality.
    pub out_dim: usize,
}

impl MlpLayer {
    /// Create a new MLP layer initialized with small random weights.
    pub fn new(in_dim: usize, out_dim: usize) -> Self {
        let n = in_dim * out_dim;
        // Xavier initialization scale.
        let scale = (2.0 / (in_dim + out_dim) as f64).sqrt();
        let weights: Vec<f64> = (0..n)
            .map(|i| {
                // Deterministic pseudo-random seeded by index.

                ((i as f64 * 1.6180339887 + 0.5).fract() - 0.5) * 2.0 * scale
            })
            .collect();
        let biases = vec![0.0; out_dim];
        Self {
            weights,
            biases,
            in_dim,
            out_dim,
        }
    }

    /// Forward pass: apply linear transform + ReLU activation.
    pub fn forward_relu(&self, input: &[f64]) -> Vec<f64> {
        assert_eq!(input.len(), self.in_dim);
        (0..self.out_dim)
            .map(|o| {
                let sum: f64 = (0..self.in_dim)
                    .map(|i| input[i] * self.weights[o * self.in_dim + i])
                    .sum();
                relu(sum + self.biases[o])
            })
            .collect()
    }

    /// Forward pass: linear transform only (no activation).
    pub fn forward_linear(&self, input: &[f64]) -> Vec<f64> {
        assert_eq!(input.len(), self.in_dim);
        (0..self.out_dim)
            .map(|o| {
                let sum: f64 = (0..self.in_dim)
                    .map(|i| input[i] * self.weights[o * self.in_dim + i])
                    .sum();
                sum + self.biases[o]
            })
            .collect()
    }

    /// Forward pass with sigmoid output activation.
    pub fn forward_sigmoid(&self, input: &[f64]) -> Vec<f64> {
        self.forward_linear(input)
            .into_iter()
            .map(sigmoid)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// NeuralRadianceField
// ---------------------------------------------------------------------------

/// Neural Radiance Field (NeRF) model for volumetric scene representation.
///
/// Implements the core NeRF architecture: a coordinate-based MLP that maps
/// 3D position + viewing direction to RGB color and volume density.
///
/// Architecture:
/// - Coarse MLP: positional-encoded (x,y,z) → density σ + feature vector
/// - Fine MLP: feature + encoded direction → RGB color
#[derive(Debug, Clone)]
pub struct NeuralRadianceField {
    /// Positional encoding for 3D position.
    pub pos_encoding: PositionalEncoding,
    /// Directional encoding for viewing direction.
    pub dir_encoding: PositionalEncoding,
    /// Coarse network layers processing position.
    pub coarse_layers: Vec<MlpLayer>,
    /// Fine network layers incorporating view direction.
    pub fine_layers: Vec<MlpLayer>,
    /// Hidden layer width.
    pub hidden_dim: usize,
    /// Near clipping distance for ray marching.
    pub near: f64,
    /// Far clipping distance for ray marching.
    pub far: f64,
}

impl NeuralRadianceField {
    /// Construct a NeRF model with standard architecture.
    ///
    /// - `num_pos_freqs`: frequency bands for position encoding (default: 10)
    /// - `num_dir_freqs`: frequency bands for direction encoding (default: 4)
    /// - `hidden_dim`: width of hidden layers (default: 256)
    /// - `near`, `far`: ray integration bounds
    pub fn new(
        num_pos_freqs: usize,
        num_dir_freqs: usize,
        hidden_dim: usize,
        near: f64,
        far: f64,
    ) -> Self {
        let pos_encoding = PositionalEncoding::new(num_pos_freqs, true);
        let dir_encoding = PositionalEncoding::new(num_dir_freqs, true);

        let pos_dim = pos_encoding.point_output_dim();
        let dir_dim = dir_encoding.direction_output_dim();

        // Coarse network: pos → hidden → hidden → (density=1, feature=hidden)
        let coarse_layers = vec![
            MlpLayer::new(pos_dim, hidden_dim),
            MlpLayer::new(hidden_dim, hidden_dim),
            MlpLayer::new(hidden_dim, hidden_dim + 1), // +1 for density output
        ];

        // Fine network: (feature + dir) → hidden → RGB(3)
        let fine_layers = vec![
            MlpLayer::new(hidden_dim + dir_dim, hidden_dim / 2),
            MlpLayer::new(hidden_dim / 2, 3), // RGB output
        ];

        Self {
            pos_encoding,
            dir_encoding,
            coarse_layers,
            fine_layers,
            hidden_dim,
            near,
            far,
        }
    }

    /// Query the NeRF at position `pos` with view direction `dir`.
    ///
    /// Returns `(density, [r, g, b])`.
    pub fn query(&self, pos: [f64; 3], dir: [f64; 3]) -> (f64, [f64; 3]) {
        // Encode inputs.
        let pos_enc = self.pos_encoding.encode_point(pos);
        let dir_enc = self.dir_encoding.encode_direction(normalize3(dir));

        // Coarse forward pass.
        let mut x = pos_enc;
        for (i, layer) in self.coarse_layers.iter().enumerate() {
            if i < self.coarse_layers.len() - 1 {
                x = layer.forward_relu(&x);
            } else {
                x = layer.forward_linear(&x);
            }
        }

        // Split density and feature vector.
        let density = softplus(x[0]); // non-negative density via softplus
        let feature: Vec<f64> = x[1..].to_vec();

        // Concatenate feature + direction encoding.
        let mut fine_input = feature;
        fine_input.extend_from_slice(&dir_enc);

        // Fine forward pass for color.
        let mut c = fine_input;
        for (i, layer) in self.fine_layers.iter().enumerate() {
            if i < self.fine_layers.len() - 1 {
                c = layer.forward_relu(&c);
            } else {
                c = layer.forward_sigmoid(&c);
            }
        }

        let color = [
            clamp(c[0], 0.0, 1.0),
            clamp(c[1], 0.0, 1.0),
            clamp(c[2], 0.0, 1.0),
        ];

        (density, color)
    }

    /// Return the expected output dimensionality for position encoding.
    pub fn pos_enc_dim(&self) -> usize {
        self.pos_encoding.point_output_dim()
    }
}

// ---------------------------------------------------------------------------
// NeRfVolume  (ray marching)
// ---------------------------------------------------------------------------

/// Sample point along a ray.
#[derive(Debug, Clone)]
pub struct RaySample {
    /// 3D world-space position of this sample.
    pub position: [f64; 3],
    /// Distance along the ray (t parameter).
    pub t: f64,
    /// RGB radiance at this sample.
    pub color: [f64; 3],
    /// Volume density at this sample.
    pub density: f64,
}

/// Volume rendering result for a single ray.
#[derive(Debug, Clone)]
pub struct RayIntegral {
    /// Accumulated RGB color.
    pub color: [f64; 3],
    /// Expected depth (weighted by transmittance × alpha).
    pub depth: f64,
    /// Accumulated alpha (opacity).
    pub alpha: f64,
    /// Number of samples used.
    pub num_samples: usize,
}

/// Neural volume renderer using ray marching through a [`NeuralRadianceField`].
#[derive(Debug, Clone)]
pub struct NeRfVolume {
    /// Number of coarse stratified samples per ray.
    pub num_coarse_samples: usize,
    /// Number of additional importance-sampled (fine) samples per ray.
    pub num_fine_samples: usize,
    /// White background blending factor.
    pub white_background: bool,
}

impl NeRfVolume {
    /// Create a new ray marcher with given sample counts.
    pub fn new(num_coarse_samples: usize, num_fine_samples: usize) -> Self {
        Self {
            num_coarse_samples,
            num_fine_samples,
            white_background: true,
        }
    }

    /// Stratified sampling along a ray from `near` to `far`.
    ///
    /// Divides the range into `n` equal bins and picks one sample per bin.
    pub fn stratified_sample(&self, near: f64, far: f64, n: usize, jitter: bool) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let t0 = near + (far - near) * (i as f64) / (n as f64);
                let t1 = near + (far - near) * ((i + 1) as f64) / (n as f64);
                if jitter {
                    // Use deterministic pseudo-jitter.
                    let frac = ((i as f64 * std::f64::consts::E + 0.5).fract()).abs();
                    t0 + frac * (t1 - t0)
                } else {
                    (t0 + t1) * 0.5
                }
            })
            .collect()
    }

    /// Integrate radiance along a ray using the NeRF volume rendering integral.
    ///
    /// Implements: C(r) = ∫ T(t) σ(t) c(t) dt where T(t) = exp(-∫₀ᵗ σ(s) ds).
    pub fn render_ray(
        &self,
        nerf: &NeuralRadianceField,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
    ) -> RayIntegral {
        let near = nerf.near;
        let far = nerf.far;
        let n = self.num_coarse_samples;
        let ts = self.stratified_sample(near, far, n, false);

        let dir_norm = normalize3(ray_dir);
        let mut samples: Vec<RaySample> = ts
            .iter()
            .map(|&t| {
                let pos = add3(ray_origin, scale3(dir_norm, t));
                let (density, color) = nerf.query(pos, dir_norm);
                RaySample {
                    position: pos,
                    t,
                    color,
                    density,
                }
            })
            .collect();

        // Sort by t (should already be sorted).
        samples.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));

        self.volume_render_integral(&samples)
    }

    /// Numerical quadrature of the volume rendering integral over sorted samples.
    pub fn volume_render_integral(&self, samples: &[RaySample]) -> RayIntegral {
        if samples.is_empty() {
            return RayIntegral {
                color: [1.0, 1.0, 1.0],
                depth: 0.0,
                alpha: 0.0,
                num_samples: 0,
            };
        }

        let mut acc_color = [0.0f64; 3];
        let mut acc_alpha = 0.0f64;
        let mut acc_depth = 0.0f64;
        let mut transmittance = 1.0f64;

        let n = samples.len();
        for (i, s) in samples.iter().enumerate() {
            // Delta t: distance to next sample (or a small step for last).
            let delta = if i + 1 < n {
                samples[i + 1].t - s.t
            } else {
                1e-3
            };

            let alpha = 1.0 - (-s.density * delta).exp();
            let weight = transmittance * alpha;

            for (ac, sc) in acc_color.iter_mut().zip(s.color.iter()) {
                *ac += weight * sc;
            }
            acc_depth += weight * s.t;
            acc_alpha += weight;
            transmittance *= 1.0 - alpha;

            if transmittance < 1e-6 {
                break;
            }
        }

        // White background blend.
        for ac in acc_color.iter_mut() {
            *ac += transmittance; // white = 1.0
        }

        RayIntegral {
            color: [
                clamp(acc_color[0], 0.0, 1.0),
                clamp(acc_color[1], 0.0, 1.0),
                clamp(acc_color[2], 0.0, 1.0),
            ],
            depth: acc_depth,
            alpha: clamp(acc_alpha, 0.0, 1.0),
            num_samples: n,
        }
    }
}

// ---------------------------------------------------------------------------
// InstantNgp
// ---------------------------------------------------------------------------

/// Hash grid level for Instant-NGP multiresolution hash encoding.
#[derive(Debug, Clone)]
pub struct HashGridLevel {
    /// Grid resolution at this level.
    pub resolution: usize,
    /// Feature vectors stored at each hash bucket.
    pub table: Vec<[f32; 2]>,
    /// Hash table size (number of entries).
    pub table_size: usize,
}

impl HashGridLevel {
    /// Create a new hash grid level.
    pub fn new(resolution: usize, table_size: usize) -> Self {
        let table = vec![[0.0f32; 2]; table_size];
        Self {
            resolution,
            table,
            table_size,
        }
    }

    /// Spatial hash function (Teschner et al., 2003).
    pub fn hash_coords(&self, ix: usize, iy: usize, iz: usize) -> usize {
        const P1: usize = 2_654_435_761;
        const P2: usize = 805_459_861;
        const P3: usize = 3_674_653_429;
        (ix.wrapping_mul(P1) ^ iy.wrapping_mul(P2) ^ iz.wrapping_mul(P3)) % self.table_size
    }

    /// Trilinearly interpolated lookup at normalized coordinates `[0,1]^3`.
    pub fn lookup(&self, x: f64, y: f64, z: f64) -> [f32; 2] {
        let r = self.resolution as f64;
        let fx = (x * r).max(0.0).min(r - 1.0001);
        let fy = (y * r).max(0.0).min(r - 1.0001);
        let fz = (z * r).max(0.0).min(r - 1.0001);

        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;

        let wx = fx - ix as f64;
        let wy = fy - iy as f64;
        let wz = fz - iz as f64;

        // Trilinear interpolation over 8 corners.
        let mut result = [0.0f32; 2];
        for dz in 0..=1usize {
            for dy in 0..=1usize {
                for dx in 0..=1usize {
                    let h = self.hash_coords(ix + dx, iy + dy, iz + dz);
                    let w = (if dx == 0 { 1.0 - wx } else { wx })
                        * (if dy == 0 { 1.0 - wy } else { wy })
                        * (if dz == 0 { 1.0 - wz } else { wz });
                    for (res_k, tbl_k) in result.iter_mut().zip(self.table[h].iter()) {
                        *res_k += (w * *tbl_k as f64) as f32;
                    }
                }
            }
        }
        result
    }
}

/// Instant-NGP multiresolution hash encoding.
///
/// Maps a 3D position to a concatenated feature vector using multiple
/// resolution hash grids for fast scene representation.
#[derive(Debug, Clone)]
pub struct InstantNgp {
    /// Hash grid levels from coarse to fine.
    pub levels: Vec<HashGridLevel>,
    /// Number of resolution levels.
    pub num_levels: usize,
    /// Base (coarsest) resolution.
    pub base_resolution: usize,
    /// Growth factor per level.
    pub level_scale: f64,
    /// Hash table size per level.
    pub hash_table_size: usize,
    /// Final MLP for density/color prediction.
    pub mlp: Vec<MlpLayer>,
}

impl InstantNgp {
    /// Create a new Instant-NGP with specified resolution pyramid.
    pub fn new(
        num_levels: usize,
        base_resolution: usize,
        level_scale: f64,
        hash_table_size: usize,
    ) -> Self {
        let levels = (0..num_levels)
            .map(|l| {
                let res = (base_resolution as f64 * level_scale.powi(l as i32)) as usize;
                HashGridLevel::new(res, hash_table_size)
            })
            .collect();

        let feature_dim = num_levels * 2; // 2 features per level
        let mlp = vec![
            MlpLayer::new(feature_dim, 64),
            MlpLayer::new(64, 16),
            MlpLayer::new(16, 4), // density + RGB
        ];

        Self {
            levels,
            num_levels,
            base_resolution,
            level_scale,
            hash_table_size,
            mlp,
        }
    }

    /// Encode a 3D position via the multiresolution hash grid.
    pub fn encode(&self, pos: [f64; 3]) -> Vec<f32> {
        let mut features = Vec::with_capacity(self.num_levels * 2);
        for level in &self.levels {
            let f = level.lookup(pos[0], pos[1], pos[2]);
            features.push(f[0]);
            features.push(f[1]);
        }
        features
    }

    /// Query density and color at a given position.
    pub fn query(&self, pos: [f64; 3]) -> (f64, [f64; 3]) {
        let enc = self.encode(pos);
        let input: Vec<f64> = enc.iter().map(|&v| v as f64).collect();

        let mut x = input;
        for (i, layer) in self.mlp.iter().enumerate() {
            if i < self.mlp.len() - 1 {
                x = layer.forward_relu(&x);
            } else {
                x = layer.forward_linear(&x);
            }
        }

        let density = softplus(x[0]);
        let color = [sigmoid(x[1]), sigmoid(x[2]), sigmoid(x[3])];
        (density, color)
    }

    /// Total number of trainable parameters.
    pub fn param_count(&self) -> usize {
        let hash_params: usize = self.levels.iter().map(|l| l.table_size * 2).sum();
        let mlp_params: usize = self
            .mlp
            .iter()
            .map(|l| l.weights.len() + l.biases.len())
            .sum();
        hash_params + mlp_params
    }
}

// ---------------------------------------------------------------------------
// GaussianSplatting
// ---------------------------------------------------------------------------

/// A single 3D Gaussian primitive for Gaussian Splatting rendering.
#[derive(Debug, Clone)]
pub struct Gaussian3D {
    /// 3D center position.
    pub position: [f64; 3],
    /// RGB color (0–1).
    pub color: [f64; 3],
    /// Opacity (0–1).
    pub opacity: f64,
    /// Scale vector for anisotropic Gaussian: `[sx, sy, sz]`.
    pub scale: [f64; 3],
    /// Rotation as a quaternion `[w, x, y, z]`.
    pub rotation: [f64; 4],
}

impl Gaussian3D {
    /// Create a new Gaussian primitive.
    pub fn new(
        position: [f64; 3],
        color: [f64; 3],
        opacity: f64,
        scale: [f64; 3],
        rotation: [f64; 4],
    ) -> Self {
        Self {
            position,
            color,
            opacity,
            scale,
            rotation,
        }
    }

    /// Rotation matrix `R` (row-major, `R[row][col]`) described by the
    /// quaternion `[w, x, y, z]` stored in [`rotation`](Self::rotation).
    ///
    /// The quaternion is normalised first, so a scaled quaternion describes
    /// the same orientation as its unit counterpart.  A zero-length or
    /// non-finite quaternion carries no orientation and yields the identity.
    /// The columns of `R` are the principal axes of the Gaussian in world
    /// space; column `k` carries the standard deviation `scale[k]`.
    pub fn rotation_matrix(&self) -> [[f64; 3]; 3] {
        let [qw, qx, qy, qz] = self.rotation;
        let norm = (qw * qw + qx * qx + qy * qy + qz * qz).sqrt();
        if !norm.is_finite() || norm < 1e-12 {
            return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        }
        let (w, x, y, z) = (qw / norm, qx / norm, qy / norm, qz / norm);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }

    /// Compute 3×3 covariance matrix from scale and rotation.
    ///
    /// Σ = R * S * S^T * R^T where S = diag(scale) and `R` is
    /// [`rotation_matrix`](Self::rotation_matrix) (the quaternion is normalised).
    pub fn covariance(&self) -> [[f64; 3]; 3] {
        let [sx, sy, sz] = self.scale;
        // S^2 (diagonal).
        let s2 = [sx * sx, sy * sy, sz * sz];

        // Σ = R * diag(s2) * R^T.
        let mut cov = [[0.0f64; 3]; 3];
        let r = self.rotation_matrix();

        // rs = R * diag(s2).
        let mut rs = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rs[i][j] = r[i][j] * s2[j];
            }
        }
        // cov = rs * R^T.
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    cov[i][j] += rs[i][k] * r[j][k];
                }
            }
        }
        cov
    }

    /// Evaluate the (unnormalised) Gaussian density at a 3D offset from the
    /// centre, honouring both the anisotropic scale and the orientation.
    ///
    /// Returns `exp(-d²/2)` with `d² = offsetᵀ Σ⁻¹ offset` the squared
    /// Mahalanobis distance.  Since `Σ = R diag(s²) Rᵀ`, its inverse is
    /// `R diag(1/s²) Rᵀ`, so `d² = Σₖ (Rᵀ offset)ₖ² / sₖ²`.  A small epsilon is
    /// added to each `sₖ²` so that a collapsed axis stays finite.
    pub fn eval_density(&self, offset: [f64; 3]) -> f64 {
        let r = self.rotation_matrix();
        // Offset expressed in the Gaussian's principal-axis frame: Rᵀ · offset.
        let mut md2 = 0.0_f64;
        for (k, s) in self.scale.iter().enumerate() {
            let local = r[0][k] * offset[0] + r[1][k] * offset[1] + r[2][k] * offset[2];
            md2 += local * local / (s * s + 1e-10);
        }
        (-0.5 * md2).exp()
    }
}

/// 3D Gaussian Splatting scene renderer.
///
/// Renders a collection of 3D Gaussians via alpha compositing in depth-sorted order.
#[derive(Debug, Clone)]
pub struct GaussianSplatting {
    /// Collection of Gaussian primitives.
    pub gaussians: Vec<Gaussian3D>,
    /// Camera near plane distance.
    pub near: f64,
    /// Spherical harmonics degree for view-dependent color.
    pub sh_degree: usize,
}

impl GaussianSplatting {
    /// Create a new Gaussian splatting scene.
    pub fn new() -> Self {
        Self {
            gaussians: Vec::new(),
            near: 0.1,
            sh_degree: 0,
        }
    }

    /// Add a Gaussian primitive to the scene.
    pub fn add(&mut self, g: Gaussian3D) {
        self.gaussians.push(g);
    }

    /// Alpha-composite all Gaussians for a given ray.
    ///
    /// Sorts Gaussians by depth along the ray direction and composites front-to-back.
    pub fn render_ray(&self, ray_origin: [f64; 3], ray_dir: [f64; 3]) -> [f64; 3] {
        let dir = normalize3(ray_dir);

        // Sort Gaussians by depth (dot product with ray direction).
        let mut indexed: Vec<(usize, f64)> = self
            .gaussians
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let offset = sub3(g.position, ray_origin);
                (i, dot3(offset, dir))
            })
            .collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Alpha compositing (front-to-back).
        let mut acc_color = [0.0f64; 3];
        let mut acc_alpha = 0.0f64;

        for (idx, depth) in &indexed {
            if *depth < self.near {
                continue;
            }
            let g = &self.gaussians[*idx];

            // Project ray-Gaussian distance onto closest point.
            let to_g = sub3(g.position, ray_origin);
            let t = dot3(to_g, dir);
            let closest = sub3(add3(ray_origin, scale3(dir, t)), g.position);
            let density = g.eval_density(closest);
            let alpha = g.opacity * density;

            // Front-to-back blending.
            let contribution = alpha * (1.0 - acc_alpha);
            for (ac, gc) in acc_color.iter_mut().zip(g.color.iter()) {
                *ac += contribution * gc;
            }
            acc_alpha += contribution;

            if acc_alpha > 0.9999 {
                break;
            }
        }

        // Blend with white background.
        let bg = 1.0 - acc_alpha;
        [acc_color[0] + bg, acc_color[1] + bg, acc_color[2] + bg]
    }

    /// Bounding sphere of all Gaussians.
    pub fn bounding_sphere(&self) -> ([f64; 3], f64) {
        if self.gaussians.is_empty() {
            return ([0.0; 3], 0.0);
        }
        let n = self.gaussians.len() as f64;
        let cx = self.gaussians.iter().map(|g| g.position[0]).sum::<f64>() / n;
        let cy = self.gaussians.iter().map(|g| g.position[1]).sum::<f64>() / n;
        let cz = self.gaussians.iter().map(|g| g.position[2]).sum::<f64>() / n;
        let center = [cx, cy, cz];
        let radius = self
            .gaussians
            .iter()
            .map(|g| len3(sub3(g.position, center)) + len3(g.scale) * 3.0)
            .fold(0.0f64, f64::max);
        (center, radius)
    }
}

impl Default for GaussianSplatting {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// NeuralSdf
// ---------------------------------------------------------------------------

/// Signed distance function query result.
#[derive(Debug, Clone)]
pub struct SdfResult {
    /// Signed distance (negative = inside surface).
    pub distance: f64,
    /// Surface normal estimated from SDF gradient.
    pub normal: [f64; 3],
    /// Whether the point is inside the surface.
    pub is_inside: bool,
}

/// Neural implicit surface via a learned signed distance function (SDF).
///
/// Uses an MLP trained to represent the SDF of a surface: positive values
/// outside, negative inside, zero on the surface.
#[derive(Debug, Clone)]
pub struct NeuralSdf {
    /// Positional encoding.
    pub encoding: PositionalEncoding,
    /// Network layers.
    pub layers: Vec<MlpLayer>,
    /// Surface iso-level (default 0.0).
    pub iso_level: f64,
    /// Finite difference epsilon for gradient estimation.
    pub grad_eps: f64,
}

impl NeuralSdf {
    /// Create a new Neural SDF network.
    pub fn new(num_freqs: usize, hidden_dim: usize) -> Self {
        let encoding = PositionalEncoding::new(num_freqs, true);
        let input_dim = encoding.point_output_dim();

        let layers = vec![
            MlpLayer::new(input_dim, hidden_dim),
            MlpLayer::new(hidden_dim, hidden_dim),
            MlpLayer::new(hidden_dim, 1), // single SDF output
        ];

        Self {
            encoding,
            layers,
            iso_level: 0.0,
            grad_eps: 1e-3,
        }
    }

    /// Evaluate the SDF at a 3D position.
    pub fn sdf(&self, pos: [f64; 3]) -> f64 {
        let enc = self.encoding.encode_point(pos);
        let mut x = enc;
        for (i, layer) in self.layers.iter().enumerate() {
            if i < self.layers.len() - 1 {
                x = layer.forward_relu(&x);
            } else {
                x = layer.forward_linear(&x);
            }
        }
        x[0]
    }

    /// Estimate surface normal via finite differences of the SDF.
    pub fn normal(&self, pos: [f64; 3]) -> [f64; 3] {
        let eps = self.grad_eps;
        let dx =
            self.sdf([pos[0] + eps, pos[1], pos[2]]) - self.sdf([pos[0] - eps, pos[1], pos[2]]);
        let dy =
            self.sdf([pos[0], pos[1] + eps, pos[2]]) - self.sdf([pos[0], pos[1] - eps, pos[2]]);
        let dz =
            self.sdf([pos[0], pos[1], pos[2] + eps]) - self.sdf([pos[0], pos[1], pos[2] - eps]);
        normalize3([dx, dy, dz])
    }

    /// Full SDF query including normal and inside/outside flag.
    pub fn query(&self, pos: [f64; 3]) -> SdfResult {
        let distance = self.sdf(pos) - self.iso_level;
        let normal = self.normal(pos);
        SdfResult {
            distance,
            normal,
            is_inside: distance < 0.0,
        }
    }

    /// Sphere tracing from `ray_origin` along `ray_dir`.
    ///
    /// Returns the intersection point and number of steps, or `None` if no
    /// intersection within `max_steps`.
    pub fn sphere_trace(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_steps: usize,
        tolerance: f64,
    ) -> Option<([f64; 3], usize)> {
        let dir = normalize3(ray_dir);
        let mut pos = ray_origin;

        for step in 0..max_steps {
            let d = self.sdf(pos) - self.iso_level;
            if d.abs() < tolerance {
                return Some((pos, step));
            }
            if d > 100.0 {
                return None; // Escaped the scene.
            }
            pos = add3(pos, scale3(dir, d));
        }
        None
    }
}

// ---------------------------------------------------------------------------
// DeepRenderImage
// ---------------------------------------------------------------------------

/// Rendered image accumulation buffer from neural rendering.
#[derive(Debug, Clone)]
pub struct DeepRenderImage {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// RGB pixel data, row-major: `pixels[y * width + x] = [r, g, b]`.
    pub pixels: Vec<[f32; 3]>,
    /// Depth buffer, one value per pixel.
    pub depth: Vec<f32>,
    /// Alpha buffer, one value per pixel.
    pub alpha: Vec<f32>,
}

impl DeepRenderImage {
    /// Create a new blank image.
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            pixels: vec![[0.0; 3]; n],
            depth: vec![f32::INFINITY; n],
            alpha: vec![0.0; n],
        }
    }

    /// Set a pixel value.
    pub fn set_pixel(&mut self, x: usize, y: usize, color: [f32; 3], depth: f32, alpha: f32) {
        if x < self.width && y < self.height {
            let idx = y * self.width + x;
            self.pixels[idx] = color;
            self.depth[idx] = depth;
            self.alpha[idx] = alpha;
        }
    }

    /// Get a pixel color (clamped).
    pub fn get_pixel(&self, x: usize, y: usize) -> [f32; 3] {
        if x < self.width && y < self.height {
            self.pixels[y * self.width + x]
        } else {
            [0.0; 3]
        }
    }

    /// Compute mean luminance of the image.
    pub fn mean_luminance(&self) -> f32 {
        if self.pixels.is_empty() {
            return 0.0;
        }
        let sum: f32 = self
            .pixels
            .iter()
            .map(|&[r, g, b]| 0.2126 * r + 0.7152 * g + 0.0722 * b)
            .sum();
        sum / self.pixels.len() as f32
    }

    /// Apply gamma correction in-place (sRGB: γ = 2.2).
    pub fn gamma_correct(&mut self, gamma: f32) {
        let inv_gamma = 1.0 / gamma;
        for px in &mut self.pixels {
            for ch in px[0..3].iter_mut() {
                *ch = ch.max(0.0).powf(inv_gamma);
            }
        }
    }

    /// Convert to 8-bit RGBA bytes.
    pub fn to_rgba8(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.width * self.height * 4);
        for (i, &[r, g, b]) in self.pixels.iter().enumerate() {
            out.push((r.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((g.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((b.clamp(0.0, 1.0) * 255.0) as u8);
            out.push((self.alpha[i].clamp(0.0, 1.0) * 255.0) as u8);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// NeuralDenoiser
// ---------------------------------------------------------------------------

/// Denoising kernel type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DenoiserKernel {
    /// Bilateral filter-style denoiser.
    Bilateral,
    /// A-trous wavelet denoiser (Intel OIDN-like).
    ATrous,
    /// Simple Gaussian blur denoiser.
    Gaussian,
}

/// OIDN-inspired neural denoiser for rendered images.
///
/// Uses a simplified A-trous wavelet convolution cascade to remove Monte Carlo
/// noise from path-traced images while preserving feature edges.
#[derive(Debug, Clone)]
pub struct NeuralDenoiser {
    /// Denoiser kernel type.
    pub kernel: DenoiserKernel,
    /// Number of filter passes (A-trous iterations).
    pub num_passes: usize,
    /// Spatial sigma for Gaussian/bilateral falloff.
    pub sigma_spatial: f64,
    /// Color sigma for bilateral edge preservation.
    pub sigma_color: f64,
    /// Filter weights (precomputed 5×5 kernel).
    pub kernel_weights: Vec<f64>,
}

impl NeuralDenoiser {
    /// Create a new A-trous denoiser.
    pub fn new(kernel: DenoiserKernel, num_passes: usize) -> Self {
        // Precompute 5×5 Gaussian kernel weights.
        let sigma = 1.0f64;
        let kernel_weights: Vec<f64> = (0..25)
            .map(|idx| {
                let kx = (idx % 5) as f64 - 2.0;
                let ky = (idx / 5) as f64 - 2.0;
                (-(kx * kx + ky * ky) / (2.0 * sigma * sigma)).exp()
            })
            .collect();
        let sum: f64 = kernel_weights.iter().sum();
        let kernel_weights = kernel_weights.iter().map(|&w| w / sum).collect();

        Self {
            kernel,
            num_passes,
            sigma_spatial: 1.0,
            sigma_color: 0.1,
            kernel_weights,
        }
    }

    /// Apply one Gaussian blur pass to an image.
    pub fn apply_gaussian_pass(&self, image: &DeepRenderImage, step: usize) -> DeepRenderImage {
        let w = image.width;
        let h = image.height;
        let mut out = DeepRenderImage::new(w, h);
        let step = step as i64;

        for y in 0..h {
            for x in 0..w {
                let mut acc = [0.0f64; 3];
                let mut weight_sum = 0.0f64;

                for ky in -2i64..=2 {
                    for kx in -2i64..=2 {
                        let nx = x as i64 + kx * step;
                        let ny = y as i64 + ky * step;
                        if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                            continue;
                        }
                        let kidx = ((ky + 2) * 5 + (kx + 2)) as usize;
                        let w_k = self.kernel_weights[kidx];
                        let src_px = image.get_pixel(nx as usize, ny as usize);
                        for k in 0..3 {
                            acc[k] += w_k * src_px[k] as f64;
                        }
                        weight_sum += w_k;
                    }
                }

                if weight_sum > 1e-10 {
                    let color = [
                        (acc[0] / weight_sum) as f32,
                        (acc[1] / weight_sum) as f32,
                        (acc[2] / weight_sum) as f32,
                    ];
                    let depth = image.depth[y * w + x];
                    let alpha = image.alpha[y * w + x];
                    out.set_pixel(x, y, color, depth, alpha);
                }
            }
        }
        out
    }

    /// Denoise an image using A-trous multi-scale filtering.
    pub fn denoise(&self, image: &DeepRenderImage) -> DeepRenderImage {
        let mut current = image.clone();
        for pass in 0..self.num_passes {
            current = self.apply_gaussian_pass(&current, 1 << pass);
        }
        current
    }
}

// ---------------------------------------------------------------------------
// StyleTransferViz
// ---------------------------------------------------------------------------

/// Style transfer configuration for scientific visualization.
#[derive(Debug, Clone)]
pub struct StyleTransferViz {
    /// Content weight for preserving scientific data features.
    pub content_weight: f64,
    /// Style weight for artistic texture transfer.
    pub style_weight: f64,
    /// Learning rate for iterative optimization.
    pub learning_rate: f64,
    /// Number of optimization iterations.
    pub num_iterations: usize,
    /// Style image gram matrix (flattened, simplified).
    pub gram_matrix: Vec<f64>,
}

impl StyleTransferViz {
    /// Create a new style transfer configuration.
    pub fn new(content_weight: f64, style_weight: f64, num_iterations: usize) -> Self {
        Self {
            content_weight,
            style_weight,
            learning_rate: 0.01,
            num_iterations,
            gram_matrix: vec![0.0; 64 * 64], // 64-dim feature gram matrix
        }
    }

    /// Compute gram matrix of a feature map (texture descriptor).
    ///
    /// G\[i\]\[j\] = sum_k F\[i\]\[k\] * F\[j\]\[k\] (inner product of feature channels).
    pub fn compute_gram(&self, features: &[Vec<f64>]) -> Vec<f64> {
        if features.is_empty() {
            return vec![];
        }
        let n = features.len();
        let k = features[0].len();
        let mut gram = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let dot: f64 = (0..k).map(|l| features[i][l] * features[j][l]).sum();
                gram[i * n + j] = dot / (k as f64);
            }
        }
        gram
    }

    /// Style loss: Frobenius norm between two gram matrices.
    pub fn style_loss(&self, gram_a: &[f64], gram_b: &[f64]) -> f64 {
        assert_eq!(gram_a.len(), gram_b.len());
        gram_a
            .iter()
            .zip(gram_b.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
    }

    /// Content loss: mean squared difference between two feature tensors.
    pub fn content_loss(&self, feat_a: &[f64], feat_b: &[f64]) -> f64 {
        assert_eq!(feat_a.len(), feat_b.len());
        let n = feat_a.len() as f64;
        feat_a
            .iter()
            .zip(feat_b.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            / n
    }

    /// Combined loss for optimization.
    pub fn total_loss(
        &self,
        content_a: &[f64],
        content_b: &[f64],
        gram_a: &[f64],
        gram_b: &[f64],
    ) -> f64 {
        self.content_weight * self.content_loss(content_a, content_b)
            + self.style_weight * self.style_loss(gram_a, gram_b)
    }
}

// ---------------------------------------------------------------------------
// LatentSpaceViz
// ---------------------------------------------------------------------------

/// A 2D or 3D point in a latent space projection.
#[derive(Debug, Clone)]
pub struct LatentPoint {
    /// Original high-dimensional embedding vector.
    pub embedding: Vec<f64>,
    /// Projected 2D position for visualization.
    pub projected: [f64; 2],
    /// Cluster label for color coding.
    pub cluster: usize,
    /// Optional trajectory index (for dynamic systems).
    pub trajectory_id: Option<usize>,
}

/// Latent space visualization using t-SNE-like dimensionality reduction.
///
/// Projects high-dimensional physics simulation states (e.g., particle
/// configurations, neural activations) to 2D for interactive exploration.
#[derive(Debug, Clone)]
pub struct LatentSpaceViz {
    /// Collection of latent points.
    pub points: Vec<LatentPoint>,
    /// Perplexity parameter for t-SNE (controls neighborhood size).
    pub perplexity: f64,
    /// Number of optimization iterations.
    pub num_iterations: usize,
    /// Learning rate for gradient descent.
    pub learning_rate: f64,
    /// Number of clusters for k-means coloring.
    pub num_clusters: usize,
}

impl LatentSpaceViz {
    /// Create a new latent space visualization.
    pub fn new(perplexity: f64, num_iterations: usize, num_clusters: usize) -> Self {
        Self {
            points: Vec::new(),
            perplexity,
            num_iterations,
            learning_rate: 200.0,
            num_clusters,
        }
    }

    /// Add a point to the latent space.
    pub fn add_point(&mut self, embedding: Vec<f64>, cluster: usize) {
        let n = self.points.len();
        // Initialize with pseudo-random 2D coordinates.
        let angle = (n as f64 * 2.399963) % (2.0 * PI); // golden angle spiral
        let radius = (n as f64 + 1.0).sqrt() * 0.1;
        let projected = [radius * angle.cos(), radius * angle.sin()];
        self.points.push(LatentPoint {
            embedding,
            projected,
            cluster,
            trajectory_id: None,
        });
    }

    /// Compute pairwise Euclidean distance between two embeddings.
    pub fn embedding_distance(&self, a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi) * (ai - bi))
            .sum::<f64>()
            .sqrt()
    }

    /// Compute pairwise distance matrix for all points.
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.points.len();
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            0.0
                        } else {
                            self.embedding_distance(
                                &self.points[i].embedding,
                                &self.points[j].embedding,
                            )
                        }
                    })
                    .collect()
            })
            .collect()
    }

    /// Simple PCA-like projection: project onto first 2 principal axes.
    ///
    /// Uses the first two dimensions of the embedding for a baseline projection
    /// (full PCA would require eigendecomposition).
    pub fn project_pca(&mut self) {
        if self.points.is_empty() {
            return;
        }
        for point in &mut self.points {
            let x = if !point.embedding.is_empty() {
                point.embedding[0]
            } else {
                0.0
            };
            let y = if point.embedding.len() > 1 {
                point.embedding[1]
            } else {
                0.0
            };
            point.projected = [x, y];
        }
        // Center the projection.
        let n = self.points.len() as f64;
        let cx = self.points.iter().map(|p| p.projected[0]).sum::<f64>() / n;
        let cy = self.points.iter().map(|p| p.projected[1]).sum::<f64>() / n;
        for p in &mut self.points {
            p.projected[0] -= cx;
            p.projected[1] -= cy;
        }
    }

    /// Assign cluster labels using a simplified k-means (2D projected space).
    pub fn assign_clusters(&mut self) {
        if self.points.is_empty() || self.num_clusters == 0 {
            return;
        }
        let k = self.num_clusters.min(self.points.len());

        // Initialize centroids at evenly-spaced angles.
        let mut centroids: Vec<[f64; 2]> = (0..k)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / k as f64;
                [a.cos(), a.sin()]
            })
            .collect();

        // Iterate k-means.
        for _ in 0..10 {
            // Assignment step.
            for p in &mut self.points {
                let closest = centroids
                    .iter()
                    .enumerate()
                    .map(|(ci, c)| {
                        let dx = p.projected[0] - c[0];
                        let dy = p.projected[1] - c[1];
                        (ci, dx * dx + dy * dy)
                    })
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|(ci, _)| ci)
                    .unwrap_or(0);
                p.cluster = closest;
            }

            // Update step.
            for (ci, centroid) in centroids.iter_mut().enumerate() {
                let members: Vec<_> = self.points.iter().filter(|p| p.cluster == ci).collect();
                if members.is_empty() {
                    continue;
                }
                let n = members.len() as f64;
                centroid[0] = members.iter().map(|p| p.projected[0]).sum::<f64>() / n;
                centroid[1] = members.iter().map(|p| p.projected[1]).sum::<f64>() / n;
            }
        }
    }

    /// Map a cluster index to an RGB color.
    pub fn cluster_color(&self, cluster: usize) -> [f32; 3] {
        // Evenly distributed hues.
        let hue = (cluster as f64 / self.num_clusters.max(1) as f64) * 360.0;
        let (r, g, b) = hsl_to_rgb(hue, 0.7, 0.5);
        [r as f32, g as f32, b as f32]
    }
}

/// Convert HSL color to RGB.
fn hsl_to_rgb(h: f64, s: f64, l: f64) -> (f64, f64, f64) {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = l - c / 2.0;
    let (r, g, b) = if h < 60.0 {
        (c, x, 0.0)
    } else if h < 120.0 {
        (x, c, 0.0)
    } else if h < 180.0 {
        (0.0, c, x)
    } else if h < 240.0 {
        (0.0, x, c)
    } else if h < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    (r + m, g + m, b + m)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;
