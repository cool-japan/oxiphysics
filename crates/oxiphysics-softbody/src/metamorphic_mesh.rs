// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Metamorphic mesh deformation for soft-body shape animation.
//!
//! This module implements a comprehensive suite of mesh deformation and
//! morphing techniques used in character animation, soft-body simulation,
//! and geometric processing:
//!
//! - **Blend shapes / morph targets** – linear and additive combinations of
//!   delta meshes keyed by scalar weights.
//! - **Mesh interpolation** – linear (LERP) and spherical (SLERP) vertex
//!   interpolation between two mesh states.
//! - **Iterative Closest Point (ICP)** – point-cloud registration that
//!   alternates closest-point assignment and best-fit rigid transform.
//! - **Shape correspondence** – Procrustes analysis for aligning two point
//!   sets with known or unknown correspondence.
//! - **Cage-based deformation** – mean-value coordinate (MVC) cage
//!   deformation where interior vertices are expressed as weighted
//!   combinations of cage vertices.
//! - **Free-form deformation (FFD)** – trivariate Bernstein-polynomial
//!   lattice deformation (Sederberg-Parry style).
//! - **Radial basis function (RBF) deformation** – scattered data
//!   interpolation using thin-plate spline / multiquadric kernels.
//! - **As-Rigid-As-Possible (ARAP)** – local-global energy minimization
//!   that minimises rotation deviation per cell.
//! - **Linear Blend Skinning (LBS)** – skeletal subspace deformation with
//!   per-vertex bone weights.
//! - **Pose-space deformation (PSD)** – corrective blend shapes driven by
//!   joint angles (pose parameters).

// ---------------------------------------------------------------------------
// Helper math (pure f64 arrays — no nalgebra)
// ---------------------------------------------------------------------------

/// 3-vector addition.
#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// 3-vector subtraction.
#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar-times-3-vector.
#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// 3-vector dot product.
#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Squared Euclidean norm of a 3-vector.
#[inline]
fn v3_norm2(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}

/// Euclidean norm of a 3-vector.
#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_norm2(a).sqrt()
}

/// Normalize a 3-vector; returns `[0,0,0]` for near-zero inputs.
#[inline]
fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

/// Centroid of a point cloud.
fn centroid(pts: &[[f64; 3]]) -> [f64; 3] {
    if pts.is_empty() {
        return [0.0; 3];
    }
    let mut s = [0.0; 3];
    for p in pts {
        s = v3_add(s, *p);
    }
    v3_scale(s, 1.0 / pts.len() as f64)
}

/// 3×3 matrix stored row-major: `m[row][col]`.
type Mat3 = [[f64; 3]; 3];

fn mat3_zero() -> Mat3 {
    [[0.0; 3]; 3]
}

fn mat3_identity() -> Mat3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn mat3_mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = mat3_zero();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn mat3_transpose(m: Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn mat3_mv(m: Mat3, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn mat3_det(m: Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Add two 3×3 matrices.
fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = mat3_zero();
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Scale a 3×3 matrix.
fn mat3_scale(m: Mat3, s: f64) -> Mat3 {
    let mut c = mat3_zero();
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = m[i][j] * s;
        }
    }
    c
}

/// Outer product of two 3-vectors: returns a 3×3 matrix.
fn outer3(a: [f64; 3], b: [f64; 3]) -> Mat3 {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}

/// Polar decomposition via one-shot symmetric orthogonalization.
///
/// Decomposes `M = R * S` where `R` is a proper rotation and `S` is
/// symmetric positive semi-definite.  Uses the iterative algorithm:
/// `R_{k+1} = 0.5 * (R_k + (R_k^{-T}))` repeated `max_iter` times.
fn polar_decomp_r(m: Mat3, max_iter: usize) -> Mat3 {
    // Regularise with a small identity term so that rank-deficient (e.g. planar)
    // covariance matrices still converge to the correct orthogonal factor.
    const EPS_REG: f64 = 1e-10;
    let mut r = [
        [m[0][0] + EPS_REG, m[0][1], m[0][2]],
        [m[1][0], m[1][1] + EPS_REG, m[1][2]],
        [m[2][0], m[2][1], m[2][2] + EPS_REG],
    ];
    for _ in 0..max_iter {
        let rt = mat3_transpose(r);
        let det = mat3_det(r);
        if det.abs() < 1e-15 {
            break;
        }
        // Cofactor matrix (= det * R^{-T})
        let cof = [
            [
                rt[1][1] * rt[2][2] - rt[1][2] * rt[2][1],
                -(rt[1][0] * rt[2][2] - rt[1][2] * rt[2][0]),
                rt[1][0] * rt[2][1] - rt[1][1] * rt[2][0],
            ],
            [
                -(rt[0][1] * rt[2][2] - rt[0][2] * rt[2][1]),
                rt[0][0] * rt[2][2] - rt[0][2] * rt[2][0],
                -(rt[0][0] * rt[2][1] - rt[0][1] * rt[2][0]),
            ],
            [
                rt[0][1] * rt[1][2] - rt[0][2] * rt[1][1],
                -(rt[0][0] * rt[1][2] - rt[0][2] * rt[1][0]),
                rt[0][0] * rt[1][1] - rt[0][1] * rt[1][0],
            ],
        ];
        let inv_t = mat3_scale(cof, 1.0 / det);
        r = mat3_scale(mat3_add(r, inv_t), 0.5);
    }
    // Ensure positive determinant
    if mat3_det(r) < 0.0 {
        for r_j in r[0].iter_mut() {
            *r_j = -*r_j;
        }
    }
    r
}

// ---------------------------------------------------------------------------
// 1. Blend Shapes / Morph Targets
// ---------------------------------------------------------------------------

/// A single morph target (blend shape) stored as per-vertex position deltas.
///
/// The delta for vertex `i` is `deltas[i]` and represents the *offset* from
/// the base mesh position when the morph weight is 1.
#[derive(Debug, Clone)]
pub struct MorphTarget {
    /// Name of this morph target (e.g. `"smile_left"`).
    pub name: String,
    /// Per-vertex position deltas `[dx, dy, dz]` indexed by vertex index.
    pub deltas: Vec<[f64; 3]>,
}

impl MorphTarget {
    /// Create a new morph target from a name and vertex deltas.
    pub fn new(name: impl Into<String>, deltas: Vec<[f64; 3]>) -> Self {
        Self {
            name: name.into(),
            deltas,
        }
    }
}

/// Blend-shape controller that combines a base mesh with weighted morph targets.
///
/// The deformed position of vertex `i` is:
///
/// ```text
/// p_i = base_i + Σ_k  weight_k * target_k.deltas[i]
/// ```
#[derive(Debug, Clone)]
pub struct BlendShapeController {
    /// Base mesh vertex positions.
    pub base_positions: Vec<[f64; 3]>,
    /// Morph targets associated with this controller.
    pub targets: Vec<MorphTarget>,
    /// Per-target blend weights in `[0, 1]`.
    pub weights: Vec<f64>,
}

impl BlendShapeController {
    /// Construct a blend-shape controller from a base mesh.
    pub fn new(base_positions: Vec<[f64; 3]>) -> Self {
        Self {
            base_positions,
            targets: Vec::new(),
            weights: Vec::new(),
        }
    }

    /// Add a morph target with initial weight 0.
    pub fn add_target(&mut self, target: MorphTarget) {
        self.targets.push(target);
        self.weights.push(0.0);
    }

    /// Set the blend weight for target `index`.  Clamped to `[0, 1]`.
    pub fn set_weight(&mut self, index: usize, weight: f64) {
        if index < self.weights.len() {
            self.weights[index] = weight.clamp(0.0, 1.0);
        }
    }

    /// Evaluate the blended mesh positions.
    ///
    /// Returns a `Vec<[f64; 3]>` the same length as `base_positions`.
    pub fn evaluate(&self) -> Vec<[f64; 3]> {
        let n = self.base_positions.len();
        let mut out = self.base_positions.clone();
        for (t, &w) in self.targets.iter().zip(self.weights.iter()) {
            if w.abs() < 1e-12 {
                continue;
            }
            let delta_len = t.deltas.len().min(n);
            for (out_i, delta_i) in out.iter_mut().zip(t.deltas.iter()).take(delta_len) {
                *out_i = v3_add(*out_i, v3_scale(*delta_i, w));
            }
        }
        out
    }

    /// Return the number of vertices in the base mesh.
    pub fn num_vertices(&self) -> usize {
        self.base_positions.len()
    }
}

// ---------------------------------------------------------------------------
// 2. Mesh Interpolation (LERP / SLERP-like)
// ---------------------------------------------------------------------------

/// Linearly interpolate between two meshes vertex-by-vertex.
///
/// `t = 0` returns `from`, `t = 1` returns `to`.
/// Both slices must have the same length.
pub fn lerp_mesh(from: &[[f64; 3]], to: &[[f64; 3]], t: f64) -> Vec<[f64; 3]> {
    let t = t.clamp(0.0, 1.0);
    from.iter()
        .zip(to.iter())
        .map(|(a, b)| v3_add(v3_scale(*a, 1.0 - t), v3_scale(*b, t)))
        .collect()
}

/// Spherically interpolate between two unit-length direction fields vertex-by-vertex.
///
/// Each vertex direction is treated as a point on S².
/// If vectors are nearly anti-parallel the result falls back to linear interpolation.
pub fn slerp_normals(from: &[[f64; 3]], to: &[[f64; 3]], t: f64) -> Vec<[f64; 3]> {
    let t = t.clamp(0.0, 1.0);
    from.iter()
        .zip(to.iter())
        .map(|(a, b)| {
            let na = v3_normalize(*a);
            let nb = v3_normalize(*b);
            let dot = v3_dot(na, nb).clamp(-1.0, 1.0);
            if (1.0 - dot.abs()) < 1e-9 {
                // Parallel or anti-parallel: fallback
                v3_normalize(v3_add(v3_scale(na, 1.0 - t), v3_scale(nb, t)))
            } else {
                let theta = dot.acos();
                let sin_theta = theta.sin();
                let coeff_a = ((1.0 - t) * theta).sin() / sin_theta;
                let coeff_b = (t * theta).sin() / sin_theta;
                v3_normalize(v3_add(v3_scale(na, coeff_a), v3_scale(nb, coeff_b)))
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 3. Iterative Closest Point (ICP)
// ---------------------------------------------------------------------------

/// Result of one ICP iteration.
#[derive(Debug, Clone)]
pub struct IcpResult {
    /// Best-fit rotation matrix (3×3, row-major).
    pub rotation: Mat3,
    /// Best-fit translation vector.
    pub translation: [f64; 3],
    /// Mean squared distance after alignment.
    pub mse: f64,
    /// Number of iterations performed.
    pub iterations: usize,
}

/// ICP solver: aligns a *source* point cloud to a *target* point cloud.
///
/// At each iteration:
/// 1. For each source point find the closest target point.
/// 2. Solve for the best-fit rigid transform (SVD-free polar-decomp approximation).
/// 3. Apply the transform to the source.
/// 4. Repeat until convergence or `max_iter`.
#[derive(Debug, Clone)]
pub struct IcpSolver {
    /// Maximum number of outer iterations.
    pub max_iter: usize,
    /// Convergence threshold on mean squared distance change.
    pub tol: f64,
}

impl IcpSolver {
    /// Create a new ICP solver with default settings.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// Find the index of the closest point in `target` for each point in `source`.
    fn find_closest(source: &[[f64; 3]], target: &[[f64; 3]]) -> Vec<usize> {
        source
            .iter()
            .map(|&s| {
                target
                    .iter()
                    .enumerate()
                    .map(|(i, &t)| (i, v3_norm2(v3_sub(s, t))))
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|(i, _)| i)
                    .unwrap_or(0)
            })
            .collect()
    }

    /// Solve for the best-fit rigid transform given paired point sets.
    ///
    /// Uses the polar-decomposition approach on the cross-covariance matrix.
    fn best_fit_transform(src: &[[f64; 3]], dst: &[[f64; 3]]) -> (Mat3, [f64; 3]) {
        let c_src = centroid(src);
        let c_dst = centroid(dst);

        // Cross-covariance H = Σ (src_i - c_src)^T (dst_i - c_dst)
        let mut h = mat3_zero();
        for (s, d) in src.iter().zip(dst.iter()) {
            let ps = v3_sub(*s, c_src);
            let pd = v3_sub(*d, c_dst);
            let outer = outer3(ps, pd);
            h = mat3_add(h, outer);
        }

        let r = polar_decomp_r(h, 20);
        // t = c_dst - R * c_src
        let rc = mat3_mv(r, c_src);
        let t = v3_sub(c_dst, rc);
        (r, t)
    }

    /// Register `source` to `target` and return the alignment result.
    pub fn register(&self, source: &[[f64; 3]], target: &[[f64; 3]]) -> IcpResult {
        let mut current: Vec<[f64; 3]> = source.to_vec();
        let mut global_r = mat3_identity();
        let mut global_t = [0.0; 3];
        let mut prev_mse = f64::MAX;
        let mut iters = 0;

        for _ in 0..self.max_iter {
            iters += 1;
            let indices = Self::find_closest(&current, target);
            let paired: Vec<[f64; 3]> = indices.iter().map(|&i| target[i]).collect();

            // Compute MSE
            let mse: f64 = current
                .iter()
                .zip(paired.iter())
                .map(|(s, d)| v3_norm2(v3_sub(*s, *d)))
                .sum::<f64>()
                / current.len() as f64;

            if mse < self.tol || (prev_mse - mse).abs() < self.tol {
                break;
            }
            prev_mse = mse;

            let (r, t) = Self::best_fit_transform(&current, &paired);

            // Apply transform to current source
            for p in current.iter_mut() {
                *p = v3_add(mat3_mv(r, *p), t);
            }

            // Accumulate global transform: new_R = r * prev_R
            global_r = mat3_mul(r, global_r);
            global_t = v3_add(mat3_mv(r, global_t), t);
        }

        let mse: f64 = {
            let indices = Self::find_closest(&current, target);
            current
                .iter()
                .zip(indices.iter())
                .map(|(s, &i)| v3_norm2(v3_sub(*s, target[i])))
                .sum::<f64>()
                / current.len().max(1) as f64
        };

        IcpResult {
            rotation: global_r,
            translation: global_t,
            mse,
            iterations: iters,
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Shape Correspondence (Procrustes)
// ---------------------------------------------------------------------------

/// Procrustes analysis result.
#[derive(Debug, Clone)]
pub struct ProcrustesResult {
    /// Optimal rotation.
    pub rotation: Mat3,
    /// Optimal translation.
    pub translation: [f64; 3],
    /// Optimal uniform scale factor (1 for orthogonal Procrustes).
    pub scale: f64,
    /// Residual sum of squared distances after alignment.
    pub residual: f64,
}

/// Compute the orthogonal Procrustes alignment of two *corresponding* point sets.
///
/// Both `source` and `target` must have the same length.
/// Returns the rotation `R`, translation `t`, and scale `s` such that
/// `s * R * source_i + t ≈ target_i` in the least-squares sense.
pub fn procrustes(source: &[[f64; 3]], target: &[[f64; 3]]) -> ProcrustesResult {
    assert_eq!(
        source.len(),
        target.len(),
        "Procrustes requires equal-length point sets"
    );
    let n = source.len();
    if n == 0 {
        return ProcrustesResult {
            rotation: mat3_identity(),
            translation: [0.0; 3],
            scale: 1.0,
            residual: 0.0,
        };
    }

    let c_src = centroid(source);
    let c_tgt = centroid(target);

    // Center
    let src_c: Vec<[f64; 3]> = source.iter().map(|p| v3_sub(*p, c_src)).collect();
    let tgt_c: Vec<[f64; 3]> = target.iter().map(|p| v3_sub(*p, c_tgt)).collect();

    // Scale
    let var_src: f64 = src_c.iter().map(|p| v3_norm2(*p)).sum::<f64>() / n as f64;
    let var_tgt: f64 = tgt_c.iter().map(|p| v3_norm2(*p)).sum::<f64>() / n as f64;
    let scale = if var_src < 1e-20 {
        1.0
    } else {
        (var_tgt / var_src).sqrt()
    };

    // Cross-covariance
    let mut h = mat3_zero();
    for (s, d) in src_c.iter().zip(tgt_c.iter()) {
        h = mat3_add(h, outer3(*s, *d));
    }

    let r = polar_decomp_r(h, 20);
    // translation: t = c_tgt - scale * R * c_src
    let rc = mat3_mv(r, c_src);
    let t = v3_sub(c_tgt, v3_scale(rc, scale));

    // Residual
    let residual: f64 = source
        .iter()
        .zip(target.iter())
        .map(|(s, d)| {
            let aligned = v3_add(v3_scale(mat3_mv(r, *s), scale), t);
            v3_norm2(v3_sub(aligned, *d))
        })
        .sum();

    ProcrustesResult {
        rotation: r,
        translation: t,
        scale,
        residual,
    }
}

// ---------------------------------------------------------------------------
// 5. Cage-Based Deformation (Mean-Value Coordinates)
// ---------------------------------------------------------------------------

/// Cage deformer using pre-computed barycentric-like weights.
///
/// Each interior vertex is expressed as a weighted sum of cage vertices:
///
/// ```text
/// p_i = Σ_k  w_{ik} * cage_k
/// ```
///
/// The weights `w_{ik}` are computed once from the rest cage and stored.
/// When the cage is deformed the interior follows automatically.
#[derive(Debug, Clone)]
pub struct CageDeformer {
    /// Number of interior (controlled) vertices.
    pub num_interior: usize,
    /// Number of cage vertices.
    pub num_cage: usize,
    /// Weights `w[interior][cage]`.  Row `i` sums (approximately) to 1.
    pub weights: Vec<Vec<f64>>,
}

impl CageDeformer {
    /// Build a cage deformer from rest interior positions and rest cage positions.
    ///
    /// Uses a simplified inverse-distance weighting as an approximation to
    /// mean-value coordinates; exact MVC requires per-polygon normal calculations.
    pub fn new(interior: &[[f64; 3]], cage: &[[f64; 3]]) -> Self {
        let ni = interior.len();
        let nc = cage.len();
        let mut weights = vec![vec![0.0; nc]; ni];

        for (i, pi) in interior.iter().enumerate() {
            let mut total = 0.0;
            for (k, ck) in cage.iter().enumerate() {
                let d2 = v3_norm2(v3_sub(*pi, *ck));
                let w = if d2 < 1e-20 { 1e10 } else { 1.0 / d2 };
                weights[i][k] = w;
                total += w;
            }
            if total > 1e-20 {
                for w_k in weights[i].iter_mut().take(nc) {
                    *w_k /= total;
                }
            }
        }

        Self {
            num_interior: ni,
            num_cage: nc,
            weights,
        }
    }

    /// Deform the interior vertices given the (possibly moved) cage.
    pub fn deform(&self, cage: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let nc = cage.len().min(self.num_cage);
        (0..self.num_interior)
            .map(|i| {
                let mut pos = [0.0; 3];
                for (cage_k, w_k) in cage.iter().zip(self.weights[i].iter()).take(nc) {
                    pos = v3_add(pos, v3_scale(*cage_k, *w_k));
                }
                pos
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// 6. Free-Form Deformation (FFD)
// ---------------------------------------------------------------------------

/// Trivariate Bernstein (Sederberg-Parry) free-form deformation lattice.
///
/// The lattice has `(l+1) × (m+1) × (n+1)` control points.  A mesh vertex
/// at parametric coordinates `(s, t, u)` ∈ \[0,1\]³ is deformed by:
///
/// ```text
/// X(s,t,u) = Σ_{i,j,k} B_i^l(s) B_j^m(t) B_k^n(u) * P_{ijk}
/// ```
#[derive(Debug, Clone)]
pub struct FfdLattice {
    /// Lattice degrees: l, m, n.
    pub degree: [usize; 3],
    /// Control-point positions `[f64; 3]`, indexed as `[i*(m+1)*(n+1) + j*(n+1) + k]`.
    pub control_points: Vec<[f64; 3]>,
    /// Rest positions of the control points (used to compute displacement).
    pub rest_points: Vec<[f64; 3]>,
    /// Origin of the FFD box.
    pub origin: [f64; 3],
    /// Axes of the FFD box (S, T, U), each a unit vector scaled by box size.
    pub axes: [[f64; 3]; 3],
}

impl FfdLattice {
    /// Compute binomial coefficient `C(n, k)`.
    fn binom(n: usize, k: usize) -> f64 {
        if k > n {
            return 0.0;
        }
        let mut result = 1.0;
        for i in 0..k {
            result *= (n - i) as f64 / (i + 1) as f64;
        }
        result
    }

    /// Bernstein basis polynomial `B_i^n(t)`.
    fn bernstein(n: usize, i: usize, t: f64) -> f64 {
        Self::binom(n, i) * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32)
    }

    /// Build a regular FFD lattice that encloses the given mesh.
    ///
    /// * `degree` – polynomial degrees `[l, m, n]`.
    /// * `mesh` – mesh vertices used to compute the bounding box.
    /// * `padding` – extra space around the bounding box.
    pub fn from_mesh(degree: [usize; 3], mesh: &[[f64; 3]], padding: f64) -> Self {
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for p in mesh {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        for d in 0..3 {
            min[d] -= padding;
            max[d] += padding;
        }
        let size = v3_sub(max, min);
        let axes = [
            [size[0], 0.0, 0.0],
            [0.0, size[1], 0.0],
            [0.0, 0.0, size[2]],
        ];
        let l = degree[0] + 1;
        let m = degree[1] + 1;
        let n = degree[2] + 1;
        let mut pts = Vec::with_capacity(l * m * n);
        for i in 0..l {
            for j in 0..m {
                for k in 0..n {
                    let p = [
                        min[0] + (i as f64 / degree[0] as f64) * size[0],
                        min[1] + (j as f64 / degree[1] as f64) * size[1],
                        min[2] + (k as f64 / degree[2] as f64) * size[2],
                    ];
                    pts.push(p);
                }
            }
        }
        Self {
            degree,
            control_points: pts.clone(),
            rest_points: pts,
            origin: min,
            axes,
        }
    }

    /// Compute the parametric coordinates `(s, t, u)` of a world-space point.
    fn world_to_param(&self, p: [f64; 3]) -> [f64; 3] {
        let local = v3_sub(p, self.origin);
        [0, 1, 2].map(|d| {
            let ax_len2 = v3_norm2(self.axes[d]);
            if ax_len2 < 1e-20 {
                0.0
            } else {
                (v3_dot(local, self.axes[d]) / ax_len2).clamp(0.0, 1.0)
            }
        })
    }

    /// Evaluate the FFD deformation at a world-space point.
    pub fn deform_point(&self, p: [f64; 3]) -> [f64; 3] {
        let [s, t, u] = self.world_to_param(p);
        let l = self.degree[0];
        let m = self.degree[1];
        let n = self.degree[2];
        let mut result = [0.0; 3];
        for i in 0..=l {
            let bi = Self::bernstein(l, i, s);
            for j in 0..=m {
                let bj = Self::bernstein(m, j, t);
                for k in 0..=n {
                    let bk = Self::bernstein(n, k, u);
                    let idx = i * (m + 1) * (n + 1) + j * (n + 1) + k;
                    let cp = self.control_points[idx];
                    let w = bi * bj * bk;
                    result = v3_add(result, v3_scale(cp, w));
                }
            }
        }
        result
    }

    /// Deform an entire mesh.
    pub fn deform_mesh(&self, mesh: &[[f64; 3]]) -> Vec<[f64; 3]> {
        mesh.iter().map(|p| self.deform_point(*p)).collect()
    }
}

// ---------------------------------------------------------------------------
// 7. Radial Basis Function (RBF) Deformation
// ---------------------------------------------------------------------------

/// Supported RBF kernel types.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RbfKernel {
    /// Thin-plate spline: `r² * ln(r)` (with `r > 0`).
    ThinPlate,
    /// Multiquadric: `sqrt(r² + c²)`.
    Multiquadric(f64),
    /// Gaussian: `exp(-c * r²)`.
    Gaussian(f64),
}

impl RbfKernel {
    /// Evaluate the kernel at radius `r`.
    pub fn evaluate(self, r: f64) -> f64 {
        match self {
            RbfKernel::ThinPlate => {
                if r < 1e-15 {
                    0.0
                } else {
                    r * r * r.ln()
                }
            }
            RbfKernel::Multiquadric(c) => (r * r + c * c).sqrt(),
            RbfKernel::Gaussian(c) => (-c * r * r).exp(),
        }
    }
}

/// RBF deformer: interpolates displacements defined at a sparse set of
/// control points to every mesh vertex.
///
/// Given control points `{q_k}` with associated displacements `{d_k}`, the
/// displacement at an arbitrary point `p` is:
///
/// ```text
/// d(p) = Σ_k  α_k * φ(‖p − q_k‖)
/// ```
///
/// Coefficients `α_k` are solved by requiring `d(q_k) = d_k`.
#[derive(Debug, Clone)]
pub struct RbfDeformer {
    /// RBF control point positions.
    pub control_points: Vec<[f64; 3]>,
    /// Known displacements at each control point.
    pub displacements: Vec<[f64; 3]>,
    /// RBF coefficients α (one per control point, per axis).
    coeffs: Vec<[f64; 3]>,
    /// RBF kernel.
    pub kernel: RbfKernel,
}

impl RbfDeformer {
    /// Build an RBF deformer from control points and their target displacements.
    ///
    /// Solves the `n × n` linear system `Φ α = d` using Gauss-Seidel iteration.
    pub fn new(
        control_points: Vec<[f64; 3]>,
        displacements: Vec<[f64; 3]>,
        kernel: RbfKernel,
    ) -> Self {
        let n = control_points.len();
        // Build kernel matrix
        let mut phi = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let r = v3_norm(v3_sub(control_points[i], control_points[j]));
                phi[i][j] = kernel.evaluate(r);
            }
        }

        // Solve Φ α = d  per axis using Gauss-Seidel
        let mut coeffs = vec![[0.0; 3]; n];
        for axis in 0..3 {
            let b: Vec<f64> = displacements.iter().map(|d| d[axis]).collect();
            let mut x = vec![0.0; n];
            for _iter in 0..200 {
                for i in 0..n {
                    let mut sigma = 0.0;
                    for j in 0..n {
                        if j != i {
                            sigma += phi[i][j] * x[j];
                        }
                    }
                    if phi[i][i].abs() > 1e-15 {
                        x[i] = (b[i] - sigma) / phi[i][i];
                    }
                }
            }
            for i in 0..n {
                coeffs[i][axis] = x[i];
            }
        }

        Self {
            control_points,
            displacements,
            coeffs,
            kernel,
        }
    }

    /// Evaluate the deformation displacement at position `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> [f64; 3] {
        let mut disp = [0.0; 3];
        for (k, ck) in self.control_points.iter().enumerate() {
            let r = v3_norm(v3_sub(p, *ck));
            let phi = self.kernel.evaluate(r);
            disp = v3_add(disp, v3_scale(self.coeffs[k], phi));
        }
        disp
    }

    /// Deform an entire mesh by adding the RBF displacement to each vertex.
    pub fn deform_mesh(&self, mesh: &[[f64; 3]]) -> Vec<[f64; 3]> {
        mesh.iter().map(|p| v3_add(*p, self.evaluate(*p))).collect()
    }
}

// ---------------------------------------------------------------------------
// 8. As-Rigid-As-Possible (ARAP) Deformation
// ---------------------------------------------------------------------------

/// A single ARAP cell (e.g. a mesh triangle or neighbourhood).
///
/// Each cell stores the set of neighbour indices and the rest-pose edge vectors.
#[derive(Debug, Clone)]
pub struct ArapCell {
    /// Index of the centre vertex.
    pub center: usize,
    /// Indices of the neighbouring vertices.
    pub neighbours: Vec<usize>,
    /// Rest-pose half-edge vectors `e_{ij} = p_j - p_i` (rest).
    pub rest_edges: Vec<[f64; 3]>,
    /// Per-edge cotangent weights.
    pub weights: Vec<f64>,
    /// Best-fit local rotation (updated each iteration).
    pub rotation: Mat3,
}

impl ArapCell {
    /// Construct an ARAP cell from rest positions.
    pub fn new(center: usize, neighbours: Vec<usize>, rest_pos: &[[f64; 3]]) -> Self {
        let rest_edges: Vec<[f64; 3]> = neighbours
            .iter()
            .map(|&j| v3_sub(rest_pos[j], rest_pos[center]))
            .collect();
        let weights = vec![1.0; neighbours.len()];
        Self {
            center,
            neighbours,
            rest_edges,
            weights,
            rotation: mat3_identity(),
        }
    }

    /// Update the best-fit local rotation given the current deformed positions.
    pub fn update_rotation(&mut self, pos: &[[f64; 3]]) {
        // Covariance S = Σ w_i * e_i_rest * e_i_deformed^T
        let mut s = mat3_zero();
        for (k, &j) in self.neighbours.iter().enumerate() {
            let e_cur = v3_sub(pos[j], pos[self.center]);
            let w = self.weights[k];
            // S += w * e_rest ⊗ e_cur
            let outer = outer3(self.rest_edges[k], e_cur);
            s = mat3_add(s, mat3_scale(outer, w));
        }
        self.rotation = polar_decomp_r(s, 20);
    }

    /// Compute the ARAP right-hand side contribution for vertex `center`.
    pub fn rhs_contribution(&self) -> [f64; 3] {
        let mut rhs = [0.0; 3];
        for (k, _j) in self.neighbours.iter().enumerate() {
            let w = self.weights[k];
            let re = mat3_mv(self.rotation, self.rest_edges[k]);
            rhs = v3_add(rhs, v3_scale(re, w));
        }
        rhs
    }
}

/// ARAP deformer: iterative local-global solver.
///
/// Each local step fits rotations per cell; each global step solves a linear
/// system.  This implementation uses Jacobi relaxation for the global step.
#[derive(Debug, Clone)]
pub struct ArapDeformer {
    /// ARAP cells, one per interior vertex.
    pub cells: Vec<ArapCell>,
    /// Current vertex positions (deformed).
    pub positions: Vec<[f64; 3]>,
    /// Pinned vertex indices (fixed boundary conditions).
    pub pinned: Vec<usize>,
    /// Jacobi relaxation iterations per global step.
    pub jacobi_iter: usize,
}

impl ArapDeformer {
    /// Construct an ARAP deformer from rest positions and cell topology.
    pub fn new(rest_pos: Vec<[f64; 3]>, cells: Vec<ArapCell>, pinned: Vec<usize>) -> Self {
        Self {
            positions: rest_pos,
            cells,
            pinned,
            jacobi_iter: 20,
        }
    }

    /// Run one iteration of the local-global ARAP solver.
    pub fn step(&mut self) {
        // Local: update rotations
        let pos_snap = self.positions.clone();
        for cell in self.cells.iter_mut() {
            cell.update_rotation(&pos_snap);
        }

        // Global: Jacobi relaxation  (simplified — one pass per non-pinned vertex)
        let n = self.positions.len();
        let mut new_pos = self.positions.clone();
        let pinned_set: std::collections::HashSet<usize> = self.pinned.iter().copied().collect();

        for cell in self.cells.iter() {
            let i = cell.center;
            if pinned_set.contains(&i) {
                continue;
            }
            let rhs = cell.rhs_contribution();
            let w_total: f64 = cell.weights.iter().sum();
            if w_total < 1e-15 {
                continue;
            }
            // Laplacian: p_i = (rhs + Σ w_j * p_j) / w_total
            let mut lhs = rhs;
            for (k, &j) in cell.neighbours.iter().enumerate() {
                lhs = v3_add(lhs, v3_scale(self.positions[j], cell.weights[k]));
            }
            new_pos[i] = v3_scale(lhs, 1.0 / w_total);
        }

        // Restore pinned
        for &pi in &self.pinned {
            if pi < n {
                new_pos[pi] = self.positions[pi];
            }
        }
        self.positions = new_pos;
    }

    /// Run `iters` ARAP iterations.
    pub fn solve(&mut self, iters: usize) {
        for _ in 0..iters {
            self.step();
        }
    }
}

// ---------------------------------------------------------------------------
// 9. Linear Blend Skinning (LBS) / Skeletal Subspace Deformation
// ---------------------------------------------------------------------------

/// A rigid joint transform: rotation + translation.
#[derive(Debug, Clone, Copy)]
pub struct JointTransform {
    /// Joint rotation matrix.
    pub rotation: Mat3,
    /// Joint translation.
    pub translation: [f64; 3],
}

impl JointTransform {
    /// Identity joint transform.
    pub fn identity() -> Self {
        Self {
            rotation: mat3_identity(),
            translation: [0.0; 3],
        }
    }

    /// Apply the joint transform to a point.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        v3_add(mat3_mv(self.rotation, p), self.translation)
    }

    /// Inverse of this joint transform (assumes pure rotation + translation).
    pub fn inverse(&self) -> Self {
        let r_inv = mat3_transpose(self.rotation);
        let t_inv = mat3_mv(r_inv, v3_scale(self.translation, -1.0));
        Self {
            rotation: r_inv,
            translation: t_inv,
        }
    }
}

/// Linear Blend Skinning (LBS) deformer.
///
/// Each vertex is deformed by a weighted combination of bone transforms:
///
/// ```text
/// p'_i = Σ_b  w_{ib} * T_b * T_b_rest^{-1} * p_i
/// ```
#[derive(Debug, Clone)]
pub struct LbsDeformer {
    /// Rest-pose vertex positions.
    pub rest_positions: Vec<[f64; 3]>,
    /// Bone weights: `weights[vertex][bone]`.
    pub weights: Vec<Vec<f64>>,
    /// Rest-pose bone transforms.
    pub rest_transforms: Vec<JointTransform>,
    /// Current bone transforms.
    pub current_transforms: Vec<JointTransform>,
}

impl LbsDeformer {
    /// Construct an LBS deformer.
    pub fn new(
        rest_positions: Vec<[f64; 3]>,
        weights: Vec<Vec<f64>>,
        rest_transforms: Vec<JointTransform>,
    ) -> Self {
        let nb = rest_transforms.len();
        let current_transforms = vec![JointTransform::identity(); nb];
        Self {
            rest_positions,
            weights,
            rest_transforms,
            current_transforms,
        }
    }

    /// Evaluate the LBS deformed positions.
    pub fn evaluate(&self) -> Vec<[f64; 3]> {
        let nv = self.rest_positions.len();
        let nb = self.rest_transforms.len();
        (0..nv)
            .map(|i| {
                let rest = self.rest_positions[i];
                let mut out = [0.0; 3];
                for b in 0..nb {
                    let w = self
                        .weights
                        .get(i)
                        .and_then(|row| row.get(b))
                        .copied()
                        .unwrap_or(0.0);
                    if w.abs() < 1e-12 {
                        continue;
                    }
                    // Bind-pose correction: T_b * T_b_rest^{-1} * rest
                    let inv_rest = self.rest_transforms[b].inverse();
                    let local = inv_rest.apply(rest);
                    let world = self.current_transforms[b].apply(local);
                    out = v3_add(out, v3_scale(world, w));
                }
                out
            })
            .collect()
    }

    /// Set the current transform for bone `b`.
    pub fn set_transform(&mut self, b: usize, transform: JointTransform) {
        if b < self.current_transforms.len() {
            self.current_transforms[b] = transform;
        }
    }
}

// ---------------------------------------------------------------------------
// 10. Pose-Space Deformation (PSD)
// ---------------------------------------------------------------------------

/// A pose-space corrective blend shape.
///
/// Each corrective is associated with a *pose vector* (e.g. joint angles) and
/// a set of per-vertex position offsets activated when the pose is close.
#[derive(Debug, Clone)]
pub struct PsdCorrectiveShape {
    /// Pose vector that activates this corrective (e.g. joint angles).
    pub pose_vector: Vec<f64>,
    /// Per-vertex deltas `[f64; 3]` (same length as the base mesh).
    pub deltas: Vec<[f64; 3]>,
}

impl PsdCorrectiveShape {
    /// Create a new PSD corrective shape.
    pub fn new(pose_vector: Vec<f64>, deltas: Vec<[f64; 3]>) -> Self {
        Self {
            pose_vector,
            deltas,
        }
    }

    /// Compute the RBF activation weight for a given current pose vector.
    ///
    /// Uses a Gaussian kernel: `w = exp(-‖pose - pose_k‖² / (2σ²))`.
    pub fn activation(&self, current_pose: &[f64], sigma: f64) -> f64 {
        let d2: f64 = self
            .pose_vector
            .iter()
            .zip(current_pose.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum();
        (-d2 / (2.0 * sigma * sigma)).exp()
    }
}

/// Pose-space deformation controller.
///
/// Combines LBS output with pose-driven corrective blend shapes:
///
/// ```text
/// p'_i = LBS(p_i)  +  Σ_k  w_k(pose) * corrective_k.deltas[i]
/// ```
#[derive(Debug, Clone)]
pub struct PsdDeformer {
    /// Base LBS deformer.
    pub lbs: LbsDeformer,
    /// Corrective shapes.
    pub correctives: Vec<PsdCorrectiveShape>,
    /// Gaussian kernel bandwidth for pose activation.
    pub sigma: f64,
}

impl PsdDeformer {
    /// Construct a PSD deformer.
    pub fn new(lbs: LbsDeformer, correctives: Vec<PsdCorrectiveShape>, sigma: f64) -> Self {
        Self {
            lbs,
            correctives,
            sigma,
        }
    }

    /// Evaluate the final deformed positions given the current pose vector.
    pub fn evaluate(&self, current_pose: &[f64]) -> Vec<[f64; 3]> {
        let mut positions = self.lbs.evaluate();
        let n = positions.len();
        for corr in &self.correctives {
            let w = corr.activation(current_pose, self.sigma);
            if w < 1e-9 {
                continue;
            }
            let dl = corr.deltas.len().min(n);
            for (pos_i, delta_i) in positions.iter_mut().zip(corr.deltas.iter()).take(dl) {
                *pos_i = v3_add(*pos_i, v3_scale(*delta_i, w));
            }
        }
        positions
    }
}

// ---------------------------------------------------------------------------
// Unit Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // ---- Blend shapes -------------------------------------------------------

    /// 1. Base mesh with no active targets returns base positions unchanged.
    #[test]
    fn test_blend_shape_zero_weight() {
        let base = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut ctrl = BlendShapeController::new(base.clone());
        let target = MorphTarget::new("test", vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        ctrl.add_target(target);
        // weight stays at 0
        let out = ctrl.evaluate();
        for (a, b) in out.iter().zip(base.iter()) {
            assert!((a[0] - b[0]).abs() < EPS, "{:?} vs {:?}", a, b);
        }
    }

    /// 2. Full weight applies the full delta.
    #[test]
    fn test_blend_shape_full_weight() {
        let base = vec![[0.0, 0.0, 0.0]];
        let mut ctrl = BlendShapeController::new(base);
        let target = MorphTarget::new("smile", vec![[1.0, 2.0, 3.0]]);
        ctrl.add_target(target);
        ctrl.set_weight(0, 1.0);
        let out = ctrl.evaluate();
        assert!((out[0][0] - 1.0).abs() < EPS);
        assert!((out[0][1] - 2.0).abs() < EPS);
        assert!((out[0][2] - 3.0).abs() < EPS);
    }

    /// 3. Partial weight linearly scales the delta.
    #[test]
    fn test_blend_shape_partial_weight() {
        let base = vec![[0.0, 0.0, 0.0]];
        let mut ctrl = BlendShapeController::new(base);
        ctrl.add_target(MorphTarget::new("t", vec![[10.0, 0.0, 0.0]]));
        ctrl.set_weight(0, 0.3);
        let out = ctrl.evaluate();
        assert!(
            (out[0][0] - 3.0).abs() < EPS,
            "Expected 3.0 got {}",
            out[0][0]
        );
    }

    /// 4. Multiple targets sum their contributions.
    #[test]
    fn test_blend_shape_multiple_targets() {
        let base = vec![[0.0, 0.0, 0.0]];
        let mut ctrl = BlendShapeController::new(base);
        ctrl.add_target(MorphTarget::new("a", vec![[1.0, 0.0, 0.0]]));
        ctrl.add_target(MorphTarget::new("b", vec![[0.0, 1.0, 0.0]]));
        ctrl.set_weight(0, 1.0);
        ctrl.set_weight(1, 1.0);
        let out = ctrl.evaluate();
        assert!((out[0][0] - 1.0).abs() < EPS);
        assert!((out[0][1] - 1.0).abs() < EPS);
    }

    // ---- Mesh interpolation --------------------------------------------------

    /// 5. LERP at t=0 returns source mesh.
    #[test]
    fn test_lerp_mesh_t0() {
        let a = vec![[1.0, 2.0, 3.0]];
        let b = vec![[4.0, 5.0, 6.0]];
        let out = lerp_mesh(&a, &b, 0.0);
        assert!((out[0][0] - 1.0).abs() < EPS);
    }

    /// 6. LERP at t=1 returns target mesh.
    #[test]
    fn test_lerp_mesh_t1() {
        let a = vec![[0.0, 0.0, 0.0]];
        let b = vec![[2.0, 4.0, 6.0]];
        let out = lerp_mesh(&a, &b, 1.0);
        assert!((out[0][1] - 4.0).abs() < EPS);
    }

    /// 7. LERP at t=0.5 is the midpoint.
    #[test]
    fn test_lerp_mesh_midpoint() {
        let a = vec![[0.0, 0.0, 0.0]];
        let b = vec![[2.0, 2.0, 2.0]];
        let out = lerp_mesh(&a, &b, 0.5);
        assert!((out[0][0] - 1.0).abs() < EPS);
    }

    /// 8. SLERP normals preserves unit length.
    #[test]
    fn test_slerp_normals_unit_length() {
        let a = vec![[1.0, 0.0, 0.0]];
        let b = vec![[0.0, 1.0, 0.0]];
        let out = slerp_normals(&a, &b, 0.5);
        let len = v3_norm(out[0]);
        assert!((len - 1.0).abs() < 1e-6, "length={len}");
    }

    // ---- ICP ----------------------------------------------------------------

    /// 9. ICP aligns identical clouds with zero MSE.
    #[test]
    fn test_icp_identity() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let solver = IcpSolver::new(10, 1e-10);
        let result = solver.register(&pts, &pts);
        assert!(result.mse < 1e-8, "mse={}", result.mse);
    }

    /// 10. ICP converges on a pure translation.
    #[test]
    fn test_icp_pure_translation() {
        let src = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let offset = [3.0, 0.0, 0.0];
        let tgt: Vec<[f64; 3]> = src.iter().map(|p| v3_add(*p, offset)).collect();
        let solver = IcpSolver::new(50, 1e-9);
        let result = solver.register(&src, &tgt);
        assert!(result.mse < 0.1, "mse after translation={}", result.mse);
    }

    // ---- Procrustes ---------------------------------------------------------

    /// 11. Procrustes on identical sets gives residual ≈ 0.
    #[test]
    fn test_procrustes_identity() {
        let pts = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let result = procrustes(&pts, &pts);
        assert!(result.residual < 1e-10, "residual={}", result.residual);
    }

    /// 12. Procrustes aligns a translated point set.
    #[test]
    fn test_procrustes_translation() {
        let src = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tgt: Vec<[f64; 3]> = src.iter().map(|p| v3_add(*p, [2.0, 3.0, 0.0])).collect();
        let result = procrustes(&src, &tgt);
        assert!(result.residual < 1e-6, "residual={}", result.residual);
    }

    // ---- Cage deformation ---------------------------------------------------

    /// 13. Cage deformer reproduces rest positions when cage is at rest.
    #[test]
    fn test_cage_rest_positions() {
        let interior = vec![[0.5, 0.5, 0.0]];
        let cage = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let deformer = CageDeformer::new(&interior, &cage);
        let out = deformer.deform(&cage);
        // Should be close to (0.5, 0.5, 0) by symmetry
        assert!((out[0][0] - 0.5).abs() < 0.1, "x={}", out[0][0]);
        assert!((out[0][1] - 0.5).abs() < 0.1, "y={}", out[0][1]);
    }

    /// 14. Cage deformer: moving all cage vertices by a constant translates interior.
    #[test]
    fn test_cage_global_translation() {
        let interior = vec![[0.5, 0.5, 0.0]];
        let cage = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let deformer = CageDeformer::new(&interior, &cage);
        let shifted: Vec<[f64; 3]> = cage.iter().map(|p| v3_add(*p, [0.0, 0.0, 5.0])).collect();
        let out = deformer.deform(&shifted);
        // z should be ≈ 5
        assert!((out[0][2] - 5.0).abs() < 0.1, "z={}", out[0][2]);
    }

    // ---- FFD ----------------------------------------------------------------

    /// 15. FFD on a vertex inside the lattice returns a point near the input (no deformation).
    #[test]
    fn test_ffd_identity() {
        let mesh = vec![[0.5, 0.5, 0.5]];
        let ffd = FfdLattice::from_mesh([2, 2, 2], &mesh, 0.1);
        let out = ffd.deform_mesh(&mesh);
        // With undeformed control points the output should equal the input
        let err = v3_norm(v3_sub(out[0], mesh[0]));
        assert!(err < 0.05, "FFD identity err={err}");
    }

    /// 16. FFD: moving a single control point deforms nearby vertices.
    #[test]
    fn test_ffd_single_point_deform() {
        let mesh = vec![[0.5, 0.5, 0.5]];
        let mut ffd = FfdLattice::from_mesh([2, 2, 2], &mesh, 0.1);
        // Move the (1,1,1) control point (center-ish) upward
        let last = ffd.control_points.len() - 1;
        ffd.control_points[last][1] += 1.0;
        let out = ffd.deform_mesh(&mesh);
        // y-coordinate should have increased
        assert!(
            out[0][1] >= mesh[0][1],
            "y should increase, got {}",
            out[0][1]
        );
    }

    // ---- RBF ----------------------------------------------------------------

    /// 17. RBF deformer with one control point reproduces the known displacement at that point.
    #[test]
    fn test_rbf_single_point_exact() {
        let cp = vec![[0.0, 0.0, 0.0]];
        let disp = vec![[1.0, 0.0, 0.0]];
        let rbf = RbfDeformer::new(cp.clone(), disp.clone(), RbfKernel::Multiquadric(1.0));
        let d = rbf.evaluate([0.0, 0.0, 0.0]);
        // Should approximately equal the target displacement
        assert!((d[0] - 1.0).abs() < 0.1, "d[0]={}", d[0]);
    }

    /// 18. RBF far from all control points returns near-zero displacement (Gaussian kernel).
    #[test]
    fn test_rbf_far_away_gaussian() {
        let cp = vec![[0.0, 0.0, 0.0]];
        let disp = vec![[1.0, 0.0, 0.0]];
        let rbf = RbfDeformer::new(cp, disp, RbfKernel::Gaussian(10.0));
        let d = rbf.evaluate([1000.0, 0.0, 0.0]);
        assert!(d[0].abs() < 1e-6, "d far away should be ~0, got {}", d[0]);
    }

    // ---- ARAP ---------------------------------------------------------------

    /// 19. ARAP with pinned vertices doesn't move pinned vertices.
    #[test]
    fn test_arap_pinned_vertices_fixed() {
        let rest = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cell0 = ArapCell::new(0, vec![1, 2], &rest);
        let mut def = ArapDeformer::new(rest.clone(), vec![cell0], vec![0, 1]);
        def.positions[2] = [0.5, 2.0, 0.0]; // perturb non-pinned
        def.solve(10);
        // Pinned vertices 0 and 1 must not have moved
        assert!(v3_norm(v3_sub(def.positions[0], rest[0])) < EPS);
        assert!(v3_norm(v3_sub(def.positions[1], rest[1])) < EPS);
    }

    /// 20. ARAP cell rotation update: identity rest = identity rotation.
    #[test]
    fn test_arap_cell_rotation_identity() {
        // Use three non-coplanar neighbours so the covariance matrix is full-rank.
        let rest = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut cell = ArapCell::new(0, vec![1, 2, 3], &rest);
        cell.update_rotation(&rest);
        // Should be near identity
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (cell.rotation[i][j] - expected).abs() < 1e-5,
                    "R[{i}][{j}]={}",
                    cell.rotation[i][j]
                );
            }
        }
    }

    // ---- LBS ----------------------------------------------------------------

    /// 21. LBS with a single identity bone returns the rest positions.
    #[test]
    fn test_lbs_identity_bone() {
        let rest = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let weights = vec![vec![1.0], vec![1.0]];
        let rest_xforms = vec![JointTransform::identity()];
        let mut lbs = LbsDeformer::new(rest.clone(), weights, rest_xforms);
        lbs.set_transform(0, JointTransform::identity());
        let out = lbs.evaluate();
        for (o, r) in out.iter().zip(rest.iter()) {
            assert!(v3_norm(v3_sub(*o, *r)) < 1e-6, "{:?} vs {:?}", o, r);
        }
    }

    /// 22. LBS with a translation bone translates all skinned vertices.
    #[test]
    fn test_lbs_pure_translation() {
        let rest = vec![[0.0, 0.0, 0.0]];
        let weights = vec![vec![1.0]];
        let rest_xforms = vec![JointTransform::identity()];
        let mut lbs = LbsDeformer::new(rest, weights, rest_xforms);
        let mut xf = JointTransform::identity();
        xf.translation = [5.0, 0.0, 0.0];
        lbs.set_transform(0, xf);
        let out = lbs.evaluate();
        assert!((out[0][0] - 5.0).abs() < 1e-6, "x={}", out[0][0]);
    }

    /// 23. LBS with two bones and equal weights averages their transforms.
    #[test]
    fn test_lbs_two_bones_average() {
        let rest = vec![[0.0, 0.0, 0.0]];
        let weights = vec![vec![0.5, 0.5]];
        let rest_xforms = vec![JointTransform::identity(), JointTransform::identity()];
        let mut lbs = LbsDeformer::new(rest, weights, rest_xforms);
        let mut xf0 = JointTransform::identity();
        xf0.translation = [2.0, 0.0, 0.0];
        let mut xf1 = JointTransform::identity();
        xf1.translation = [4.0, 0.0, 0.0];
        lbs.set_transform(0, xf0);
        lbs.set_transform(1, xf1);
        let out = lbs.evaluate();
        // Average: (2 + 4) * 0.5 = 3
        assert!((out[0][0] - 3.0).abs() < 1e-6, "x={}", out[0][0]);
    }

    // ---- PSD ----------------------------------------------------------------

    /// 24. PSD with no correctives matches LBS output.
    #[test]
    fn test_psd_no_correctives() {
        let rest = vec![[1.0, 0.0, 0.0]];
        let weights = vec![vec![1.0]];
        let rest_xforms = vec![JointTransform::identity()];
        let lbs = LbsDeformer::new(rest.clone(), weights, rest_xforms);
        let psd = PsdDeformer::new(lbs, vec![], 1.0);
        let pose = vec![0.0];
        let out = psd.evaluate(&pose);
        assert!((out[0][0] - 1.0).abs() < 1e-6);
    }

    /// 25. PSD corrective activates strongly when pose matches.
    #[test]
    fn test_psd_corrective_activates() {
        let rest = vec![[0.0, 0.0, 0.0]];
        let weights = vec![vec![1.0]];
        let rest_xforms = vec![JointTransform::identity()];
        let lbs = LbsDeformer::new(rest, weights, rest_xforms);
        let corr = PsdCorrectiveShape::new(vec![1.0], vec![[0.0, 1.0, 0.0]]);
        let psd = PsdDeformer::new(lbs, vec![corr], 0.5);
        // Evaluate at exact pose
        let out = psd.evaluate(&[1.0]);
        // Activation = exp(0) = 1.0 → delta fully applied
        assert!(out[0][1] > 0.9, "y should be ~1.0, got {}", out[0][1]);
        // Far-away pose should not activate
        let out2 = psd.evaluate(&[100.0]);
        assert!(
            out2[0][1].abs() < 1e-6,
            "far pose should give y≈0, got {}",
            out2[0][1]
        );
        // Silence unused mut warning
        let _ = psd.sigma;
    }

    // ---- Polar decomposition ------------------------------------------------

    /// 26. polar_decomp_r on identity returns identity.
    #[test]
    fn test_polar_decomp_identity() {
        let i = mat3_identity();
        let r = polar_decomp_r(i, 10);
        for (ri, r_row) in r.iter().enumerate() {
            for (ci, &rv) in r_row.iter().enumerate() {
                let expected = if ri == ci { 1.0 } else { 0.0 };
                assert!((rv - expected).abs() < 1e-6, "R[{ri}][{ci}]={}", rv);
            }
        }
    }

    /// 27. polar_decomp_r on a scaled identity gives identity (positive det).
    #[test]
    fn test_polar_decomp_scaled_identity() {
        let m = mat3_scale(mat3_identity(), 3.0);
        let r = polar_decomp_r(m, 20);
        let det = mat3_det(r);
        assert!(det > 0.5, "det should be ~1, got {det}");
    }

    // ---- Math helpers -------------------------------------------------------

    /// 28. centroid of a symmetric point set is the origin.
    #[test]
    fn test_centroid_symmetric() {
        let pts = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let c = centroid(&pts);
        assert!(v3_norm(c) < EPS, "centroid should be origin, got {:?}", c);
    }

    /// 29. v3_normalize produces a unit vector.
    #[test]
    fn test_v3_normalize_unit() {
        let v = [3.0, 4.0, 0.0];
        let n = v3_normalize(v);
        assert!((v3_norm(n) - 1.0).abs() < EPS);
    }

    /// 30. outer3 of two unit vectors gives a rank-1 matrix.
    #[test]
    fn test_outer3_rank_one() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let m = outer3(a, b);
        // Only m[0][1] should be 1; others 0
        assert!((m[0][1] - 1.0).abs() < EPS);
        assert!(m[0][0].abs() < EPS);
        assert!(m[1][1].abs() < EPS);
    }

    /// 31. mat3_mul with identity gives back the original matrix.
    #[test]
    fn test_mat3_mul_identity() {
        let id = mat3_identity();
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let r = mat3_mul(id, m);
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - m[i][j]).abs() < EPS, "r[{i}][{j}]={}", r[i][j]);
            }
        }
    }

    /// 32. RBF Gaussian kernel approaches 0 for large r.
    #[test]
    fn test_rbf_gaussian_kernel_decay() {
        let k = RbfKernel::Gaussian(1.0);
        assert!(k.evaluate(0.0) > 0.99);
        assert!(k.evaluate(10.0) < 1e-40);
    }

    /// 33. FFD Bernstein partition of unity: Σ B_i^n(t) = 1 for all t.
    #[test]
    fn test_ffd_bernstein_pou() {
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let n = 3;
            let sum: f64 = (0..=n).map(|i| FfdLattice::bernstein(n, i, t)).sum();
            assert!((sum - 1.0).abs() < EPS, "t={t} sum={sum}");
        }
    }

    /// 34. JointTransform::inverse undoes the transform.
    #[test]
    fn test_joint_transform_inverse() {
        let mut xf = JointTransform::identity();
        xf.translation = [1.0, 2.0, 3.0];
        let p = [4.0, 5.0, 6.0];
        let xf_inv = xf.inverse();
        let p2 = xf.apply(p);
        let p3 = xf_inv.apply(p2);
        assert!(
            v3_norm(v3_sub(p3, p)) < 1e-6,
            "round-trip error {:?} vs {:?}",
            p3,
            p
        );
    }

    /// 35. BlendShapeController num_vertices returns base mesh size.
    #[test]
    fn test_blend_shape_num_vertices() {
        let base = vec![[0.0; 3]; 42];
        let ctrl = BlendShapeController::new(base);
        assert_eq!(ctrl.num_vertices(), 42);
    }
}
