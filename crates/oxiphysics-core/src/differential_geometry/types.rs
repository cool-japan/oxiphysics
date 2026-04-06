//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::functions::*;

/// A 2-form on R³.
#[derive(Debug, Clone, Copy)]
pub struct TwoForm {
    /// Component of dx∧dy.
    pub dxdy: f64,
    /// Component of dy∧dz.
    pub dydz: f64,
    /// Component of dz∧dx.
    pub dzdx: f64,
}
impl TwoForm {
    /// Create a 2-form from components.
    pub fn new(dxdy: f64, dydz: f64, dzdx: f64) -> Self {
        Self { dxdy, dydz, dzdx }
    }
    /// Zero 2-form.
    pub fn zero() -> Self {
        Self {
            dxdy: 0.0,
            dydz: 0.0,
            dzdx: 0.0,
        }
    }
    /// Evaluate on two tangent vectors.
    pub fn evaluate(&self, u: Vec3, v: Vec3) -> f64 {
        self.dxdy * (u[0] * v[1] - u[1] * v[0])
            + self.dydz * (u[1] * v[2] - u[2] * v[1])
            + self.dzdx * (u[2] * v[0] - u[0] * v[2])
    }
    /// Dual vector field (Hodge dual in R³).
    pub fn hodge_dual(&self) -> Vec3 {
        [self.dydz, self.dzdx, self.dxdy]
    }
}
/// A metric tensor g_{ij} for an N-dimensional manifold (N <= MAX_DIM).
///
/// Stored as a symmetric matrix in row-major order.
#[derive(Debug, Clone)]
pub struct MetricTensorND {
    /// Dimension of the manifold.
    pub dim: usize,
    /// Lower-index components g_{ij}, stored as \[dim\]\[dim\].
    pub g: [[f64; MAX_DIM]; MAX_DIM],
}
impl MetricTensorND {
    /// Create a new metric tensor from a dim x dim slice.
    pub fn new(dim: usize, components: &[[f64; MAX_DIM]; MAX_DIM]) -> Self {
        assert!(dim <= MAX_DIM, "Dimension exceeds MAX_DIM");
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        for i in 0..dim {
            for j in 0..dim {
                g[i][j] = components[i][j];
            }
        }
        Self { dim, g }
    }
    /// Create a flat (Euclidean/Minkowski) metric of given signature.
    ///
    /// `signature` slice: +1.0 for spacelike, -1.0 for timelike.
    pub fn from_signature(signature: &[f64]) -> Self {
        let dim = signature.len();
        assert!(dim <= MAX_DIM);
        let mut g = [[0.0; MAX_DIM]; MAX_DIM];
        for i in 0..dim {
            g[i][i] = signature[i];
        }
        Self { dim, g }
    }
    /// Create a Euclidean metric (identity) of given dimension.
    pub fn euclidean(dim: usize) -> Self {
        Self::from_signature(&vec![1.0; dim])
    }
    /// Create a Minkowski metric diag(-1,1,1,1) for dim=4.
    pub fn minkowski() -> Self {
        Self::from_signature(&[-1.0, 1.0, 1.0, 1.0])
    }
    /// Compute the inverse metric g^{ij}.
    ///
    /// Uses Gauss-Jordan elimination for general dimension.
    pub fn inverse(&self) -> [[f64; MAX_DIM]; MAX_DIM] {
        let n = self.dim;
        let mut aug = [[0.0_f64; 2 * MAX_DIM]; MAX_DIM];
        for i in 0..n {
            for j in 0..n {
                aug[i][j] = self.g[i][j];
            }
            aug[i][n + i] = 1.0;
        }
        for col in 0..n {
            let mut max_row = col;
            let mut max_val = aug[col][col].abs();
            for row in (col + 1)..n {
                if aug[row][col].abs() > max_val {
                    max_val = aug[row][col].abs();
                    max_row = row;
                }
            }
            aug.swap(col, max_row);
            let pivot = aug[col][col];
            assert!(pivot.abs() > 1e-15, "Singular metric tensor");
            let inv_pivot = 1.0 / pivot;
            for j in 0..(2 * n) {
                aug[col][j] *= inv_pivot;
            }
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = aug[row][col];
                for j in 0..(2 * n) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }
        let mut inv = [[0.0_f64; MAX_DIM]; MAX_DIM];
        for i in 0..n {
            for j in 0..n {
                inv[i][j] = aug[i][n + j];
            }
        }
        inv
    }
    /// Determinant of the metric tensor.
    pub fn determinant(&self) -> f64 {
        match self.dim {
            1 => self.g[0][0],
            2 => self.g[0][0] * self.g[1][1] - self.g[0][1] * self.g[1][0],
            3 => {
                let m = self.g;
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
            }
            4 => {
                let m = self.g;
                let mut det = 0.0;
                for j in 0..4 {
                    let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
                    let minor = self.minor_3x3(0, j);
                    det += sign * m[0][j] * minor;
                }
                det
            }
            _ => panic!("Unsupported dimension"),
        }
    }
    /// 3x3 minor obtained by deleting row `skip_r` and column `skip_c`.
    fn minor_3x3(&self, skip_r: usize, skip_c: usize) -> f64 {
        let mut sub = [[0.0_f64; 3]; 3];
        let mut si = 0;
        for i in 0..4 {
            if i == skip_r {
                continue;
            }
            let mut sj = 0;
            for j in 0..4 {
                if j == skip_c {
                    continue;
                }
                sub[si][sj] = self.g[i][j];
                sj += 1;
            }
            si += 1;
        }
        sub[0][0] * (sub[1][1] * sub[2][2] - sub[1][2] * sub[2][1])
            - sub[0][1] * (sub[1][0] * sub[2][2] - sub[1][2] * sub[2][0])
            + sub[0][2] * (sub[1][0] * sub[2][1] - sub[1][1] * sub[2][0])
    }
    /// Volume element sqrt(|det(g)|).
    pub fn volume_element(&self) -> f64 {
        self.determinant().abs().sqrt()
    }
    /// Raise an index: v^i = g^{ij} v_j.
    pub fn raise_index(&self, v_lower: &[f64; MAX_DIM]) -> [f64; MAX_DIM] {
        let ginv = self.inverse();
        let mut result = [0.0; MAX_DIM];
        for i in 0..self.dim {
            for j in 0..self.dim {
                result[i] += ginv[i][j] * v_lower[j];
            }
        }
        result
    }
    /// Lower an index: v_i = g_{ij} v^j.
    pub fn lower_index(&self, v_upper: &[f64; MAX_DIM]) -> [f64; MAX_DIM] {
        let mut result = [0.0; MAX_DIM];
        for i in 0..self.dim {
            for j in 0..self.dim {
                result[i] += self.g[i][j] * v_upper[j];
            }
        }
        result
    }
    /// Inner product of two vectors using this metric: g_{ij} u^i v^j.
    pub fn inner_product(&self, u: &[f64; MAX_DIM], v: &[f64; MAX_DIM]) -> f64 {
        let mut result = 0.0;
        for i in 0..self.dim {
            for j in 0..self.dim {
                result += self.g[i][j] * u[i] * v[j];
            }
        }
        result
    }
    /// Norm of a vector: sqrt(|g_{ij} v^i v^j|).
    pub fn norm(&self, v: &[f64; MAX_DIM]) -> f64 {
        self.inner_product(v, v).abs().sqrt()
    }
}
/// Second fundamental form coefficients (L, M, N).
#[derive(Debug, Clone, Copy)]
pub struct SecondFundamentalForm {
    /// L = -∂r/∂u · ∂n/∂u = ruu · n.
    pub l: f64,
    /// M = -∂r/∂u · ∂n/∂v = ruv · n.
    pub m: f64,
    /// N = -∂r/∂v · ∂n/∂v = rvv · n.
    pub n: f64,
}
impl SecondFundamentalForm {
    /// Compute from second derivatives and unit normal.
    pub fn from_derivatives(ruu: Vec3, ruv: Vec3, rvv: Vec3, unit_normal: Vec3) -> Self {
        Self {
            l: dot3(ruu, unit_normal),
            m: dot3(ruv, unit_normal),
            n: dot3(rvv, unit_normal),
        }
    }
    /// Gaussian curvature K = (LN - M²) / (EG - F²).
    pub fn gaussian_curvature(&self, g: &RiemannianMetric) -> f64 {
        let det_g = g.det();
        if det_g.abs() < 1e-15 {
            return 0.0;
        }
        (self.l * self.n - self.m * self.m) / det_g
    }
    /// Mean curvature H = (EN - 2FM + GL) / (2(EG - F²)).
    pub fn mean_curvature(&self, g: &RiemannianMetric) -> f64 {
        let det_g = g.det();
        if det_g.abs() < 1e-15 {
            return 0.0;
        }
        (g.e * self.n - 2.0 * g.f * self.m + g.g * self.l) / (2.0 * det_g)
    }
    /// Principal curvatures κ₁, κ₂ (eigenvalues of shape operator).
    pub fn principal_curvatures(&self, g: &RiemannianMetric) -> (f64, f64) {
        let h = self.mean_curvature(g);
        let k = self.gaussian_curvature(g);
        let disc = (h * h - k).max(0.0);
        let sqrt_disc = disc.sqrt();
        (h - sqrt_disc, h + sqrt_disc)
    }
    /// Shape index S ∈ \[-1, 1\] classification.
    pub fn shape_index(&self, g: &RiemannianMetric) -> f64 {
        let (k1, k2) = self.principal_curvatures(g);
        if (k1 - k2).abs() < 1e-15 {
            return 0.0;
        }
        (2.0 / PI) * ((k2 + k1) / (k2 - k1)).atan()
    }
}
/// State for geodesic ODE on a surface: (u, v, du/ds, dv/ds).
#[derive(Debug, Clone, Copy)]
pub struct GeodesicState {
    /// Parameter u.
    pub u: f64,
    /// Parameter v.
    pub v: f64,
    /// Tangent du/ds.
    pub du: f64,
    /// Tangent dv/ds.
    pub dv: f64,
}
impl GeodesicState {
    /// Create initial geodesic state.
    pub fn new(u: f64, v: f64, du: f64, dv: f64) -> Self {
        Self { u, v, du, dv }
    }
}
/// Riemannian metric on a parametric surface, defined by the first fundamental form.
#[derive(Debug, Clone, Copy)]
pub struct RiemannianMetric {
    /// g₁₁ = E: E = ∂r/∂u · ∂r/∂u
    pub e: f64,
    /// g₁₂ = F: F = ∂r/∂u · ∂r/∂v
    pub f: f64,
    /// g₂₂ = G: G = ∂r/∂v · ∂r/∂v
    pub g: f64,
}
impl RiemannianMetric {
    /// Create metric from first fundamental form coefficients.
    pub fn new(e: f64, f: f64, g: f64) -> Self {
        Self { e, f, g }
    }
    /// Compute from tangent vectors.
    pub fn from_tangents(ru: Vec3, rv: Vec3) -> Self {
        Self {
            e: dot3(ru, ru),
            f: dot3(ru, rv),
            g: dot3(rv, rv),
        }
    }
    /// Determinant of the metric tensor (EG - F²).
    pub fn det(&self) -> f64 {
        self.e * self.g - self.f * self.f
    }
    /// Area element √(EG - F²).
    pub fn area_element(&self) -> f64 {
        self.det().sqrt()
    }
    /// Length of a tangent vector given its components (du, dv).
    pub fn length(&self, du: f64, dv: f64) -> f64 {
        (self.e * du * du + 2.0 * self.f * du * dv + self.g * dv * dv).sqrt()
    }
    /// Angle between two tangent vectors in parameter space.
    pub fn angle(&self, du1: f64, dv1: f64, du2: f64, dv2: f64) -> f64 {
        let dot = self.e * du1 * du2 + self.f * (du1 * dv2 + dv1 * du2) + self.g * dv1 * dv2;
        let l1 = self.length(du1, dv1);
        let l2 = self.length(du2, dv2);
        if l1 < 1e-15 || l2 < 1e-15 {
            return 0.0;
        }
        (dot / (l1 * l2)).clamp(-1.0, 1.0).acos()
    }
    /// Inverse metric components g^ij.
    pub fn inverse(&self) -> (f64, f64, f64) {
        let d = self.det();
        if d.abs() < 1e-15 {
            return (0.0, 0.0, 0.0);
        }
        (self.g / d, -self.f / d, self.e / d)
    }
    /// Christoffel symbol Γ^1_11 for the metric (requires second derivatives).
    ///
    /// Parameters are second-order derivatives of the position vector.
    pub fn christoffel_111(&self, ruu: Vec3, _ruv: Vec3, _rvv: Vec3, ru: Vec3, _rv: Vec3) -> f64 {
        let (g11_inv, g12_inv, _g22_inv) = self.inverse();
        let eu = 2.0 * dot3(ruu, ru);
        let guu = dot3(ruu, _rvv);
        let e_u = eu;
        let f_u = 2.0 * dot3(ruu, _ruv);
        let e_v = 2.0 * dot3(_ruv, ru);
        let _ = guu;
        0.5 * g11_inv * e_u + g12_inv * (e_v * 0.5 - f_u)
    }
}
/// Element of the SO(3) rotation group, stored as a rotation matrix.
#[derive(Debug, Clone, Copy)]
pub struct So3 {
    /// Underlying rotation matrix (row-major).
    pub mat: Mat3,
}
impl So3 {
    /// Identity rotation.
    pub fn identity() -> Self {
        Self {
            mat: mat3_identity(),
        }
    }
    /// Construct from a rotation matrix (no orthogonality check).
    pub fn from_matrix(m: Mat3) -> Self {
        Self { mat: m }
    }
    /// Rodrigues' rotation formula: exp(\[ω\]×) for axis-angle vector ω.
    ///
    /// `omega` – axis-angle vector; its norm is the rotation angle in radians.
    pub fn exp(omega: Vec3) -> Self {
        let theta = norm3(omega);
        if theta < 1e-12 {
            return Self::identity();
        }
        let axis = scale3(omega, 1.0 / theta);
        let k = skew3(axis);
        let k2 = mat3_mul(k, k);
        let i = mat3_identity();
        let term1 = mat3_scale(k, theta.sin());
        let term2 = mat3_scale(k2, 1.0 - theta.cos());
        Self {
            mat: mat3_add(mat3_add(i, term1), term2),
        }
    }
    /// Logarithm map: returns axis-angle vector ω such that exp(ω) ≈ self.
    pub fn log(&self) -> Vec3 {
        let r = self.mat;
        let trace = mat3_trace(r);
        let cos_theta = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        if theta.abs() < 1e-12 {
            return [0.0, 0.0, 0.0];
        }
        if (theta - PI).abs() < 1e-6 {
            let rxx = r[0][0];
            let ryy = r[1][1];
            let rzz = r[2][2];
            let (axis_x, axis_y, axis_z) = if rxx >= ryy && rxx >= rzz {
                let x = ((rxx + 1.0) / 2.0).sqrt();
                (
                    x,
                    (r[0][1] + r[1][0]) / (4.0 * x),
                    (r[0][2] + r[2][0]) / (4.0 * x),
                )
            } else if ryy >= rzz {
                let y = ((ryy + 1.0) / 2.0).sqrt();
                (
                    (r[0][1] + r[1][0]) / (4.0 * y),
                    y,
                    (r[1][2] + r[2][1]) / (4.0 * y),
                )
            } else {
                let z = ((rzz + 1.0) / 2.0).sqrt();
                (
                    (r[0][2] + r[2][0]) / (4.0 * z),
                    (r[1][2] + r[2][1]) / (4.0 * z),
                    z,
                )
            };
            return scale3([axis_x, axis_y, axis_z], PI);
        }
        let factor = theta / (2.0 * theta.sin());
        let rt = mat3_transpose(r);
        let skew_part = mat3_scale(mat3_add(r, mat3_scale(rt, -1.0)), 0.5);
        scale3(vee3(skew_part), factor * 2.0)
    }
    /// Compose two rotations: self · other.
    pub fn compose(&self, other: &So3) -> So3 {
        So3 {
            mat: mat3_mul(self.mat, other.mat),
        }
    }
    /// Inverse rotation (transpose for SO(3)).
    pub fn inverse(&self) -> So3 {
        So3 {
            mat: mat3_transpose(self.mat),
        }
    }
    /// Rotate a 3-vector.
    pub fn rotate(&self, v: Vec3) -> Vec3 {
        mat3_mul_vec3(self.mat, v)
    }
    /// Adjoint representation of SO(3) — same as the rotation matrix.
    pub fn adjoint(&self) -> Mat3 {
        self.mat
    }
    /// Geodesic distance between two rotations.
    pub fn dist(&self, other: &So3) -> f64 {
        let rel = self.inverse().compose(other);
        norm3(rel.log())
    }
    /// Interpolation between two rotations (SLERP via Lie algebra).
    ///
    /// `t` = 0 → self, `t` = 1 → other.
    pub fn slerp(&self, other: &So3, t: f64) -> So3 {
        let rel = self.inverse().compose(other);
        let omega = rel.log();
        let step = So3::exp(scale3(omega, t));
        self.compose(&step)
    }
    /// Convert to unit quaternion \[w, x, y, z\].
    pub fn to_quaternion(&self) -> Quat {
        let m = self.mat;
        let trace = mat3_trace(m);
        if trace > 0.0 {
            let s = (trace + 1.0).sqrt() * 2.0;
            let w = 0.25 * s;
            let x = (m[2][1] - m[1][2]) / s;
            let y = (m[0][2] - m[2][0]) / s;
            let z = (m[1][0] - m[0][1]) / s;
            [w, x, y, z]
        } else if m[0][0] > m[1][1] && m[0][0] > m[2][2] {
            let s = (1.0 + m[0][0] - m[1][1] - m[2][2]).sqrt() * 2.0;
            let w = (m[2][1] - m[1][2]) / s;
            let x = 0.25 * s;
            let y = (m[0][1] + m[1][0]) / s;
            let z = (m[0][2] + m[2][0]) / s;
            [w, x, y, z]
        } else if m[1][1] > m[2][2] {
            let s = (1.0 + m[1][1] - m[0][0] - m[2][2]).sqrt() * 2.0;
            let w = (m[0][2] - m[2][0]) / s;
            let x = (m[0][1] + m[1][0]) / s;
            let y = 0.25 * s;
            let z = (m[1][2] + m[2][1]) / s;
            [w, x, y, z]
        } else {
            let s = (1.0 + m[2][2] - m[0][0] - m[1][1]).sqrt() * 2.0;
            let w = (m[1][0] - m[0][1]) / s;
            let x = (m[0][2] + m[2][0]) / s;
            let y = (m[1][2] + m[2][1]) / s;
            let z = 0.25 * s;
            [w, x, y, z]
        }
    }
    /// Construct from unit quaternion \[w, x, y, z\].
    pub fn from_quaternion(q: Quat) -> Self {
        let [w, x, y, z] = q;
        let mat = [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - z * w),
                2.0 * (x * z + y * w),
            ],
            [
                2.0 * (x * y + z * w),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - x * w),
            ],
            [
                2.0 * (x * z - y * w),
                2.0 * (y * z + x * w),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ];
        Self { mat }
    }
    /// Left Jacobian of SO(3) — maps Lie algebra to tangent space.
    pub fn left_jacobian(omega: Vec3) -> Mat3 {
        let theta = norm3(omega);
        if theta < 1e-10 {
            return mat3_identity();
        }
        let axis = scale3(omega, 1.0 / theta);
        let k = skew3(axis);
        let k2 = mat3_mul(k, k);
        let a = theta.sin() / theta;
        let b = (1.0 - theta.cos()) / theta;
        let i = mat3_identity();
        mat3_add(mat3_add(i, mat3_scale(k, b)), mat3_scale(k2, 1.0 - a))
    }
    /// Right Jacobian of SO(3) — transpose of left Jacobian.
    pub fn right_jacobian(omega: Vec3) -> Mat3 {
        mat3_transpose(So3::left_jacobian(omega))
    }
}
/// Levi-Civita connection (Christoffel symbols of the second kind).
///
/// Computes Gamma^sigma_{mu nu} from the metric and its numerical derivatives.
///
/// `metric_fn` returns the metric tensor components at a given coordinate point.
/// `point` is the coordinate at which to evaluate.
/// `dim` is the manifold dimension.
/// `h` is the finite-difference step size.
#[allow(clippy::too_many_arguments)]
pub struct LeviCivitaConnection {
    /// Christoffel symbols Gamma^sigma_{mu nu}, stored as \[sigma\]\[mu\]\[nu\].
    pub christoffel: [[[f64; MAX_DIM]; MAX_DIM]; MAX_DIM],
    /// Dimension of the manifold.
    pub dim: usize,
}
impl LeviCivitaConnection {
    /// Compute the Levi-Civita connection from a metric function at a point.
    ///
    /// `metric_fn(x)` returns the metric tensor g_{ij} at coordinates x.
    pub fn from_metric_fn<F>(metric_fn: &F, point: &[f64; MAX_DIM], dim: usize, h: f64) -> Self
    where
        F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
    {
        let mut dg = [[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM];
        for k in 0..dim {
            let mut p_plus = *point;
            let mut p_minus = *point;
            p_plus[k] += h;
            p_minus[k] -= h;
            let g_plus = metric_fn(&p_plus);
            let g_minus = metric_fn(&p_minus);
            for i in 0..dim {
                for j in 0..dim {
                    dg[k][i][j] = (g_plus[i][j] - g_minus[i][j]) / (2.0 * h);
                }
            }
        }
        let g_here = metric_fn(point);
        let metric = MetricTensorND { dim, g: g_here };
        let ginv = metric.inverse();
        let mut christoffel = [[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM];
        for sigma in 0..dim {
            for mu in 0..dim {
                for nu in 0..dim {
                    let mut val = 0.0;
                    for rho in 0..dim {
                        val += ginv[sigma][rho]
                            * (dg[nu][rho][mu] + dg[mu][rho][nu] - dg[rho][mu][nu]);
                    }
                    christoffel[sigma][mu][nu] = 0.5 * val;
                }
            }
        }
        Self { christoffel, dim }
    }
    /// Access Gamma^sigma_{mu nu}.
    pub fn gamma(&self, sigma: usize, mu: usize, nu: usize) -> f64 {
        self.christoffel[sigma][mu][nu]
    }
}
/// Riemann curvature tensor R^rho_{sigma mu nu}.
///
/// Computed from Christoffel symbols and their derivatives.
#[derive(Debug, Clone)]
pub struct RiemannTensor {
    /// Components R^rho_{sigma mu nu}.
    /// Indexed as \[rho\]\[sigma\]\[mu\]\[nu\].
    pub components: [[[[f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM],
    /// Dimension of the manifold.
    pub dim: usize,
}
impl RiemannTensor {
    /// Compute the Riemann tensor from a metric function at a point.
    ///
    /// Uses numerical differentiation of Christoffel symbols.
    pub fn from_metric_fn<F>(metric_fn: &F, point: &[f64; MAX_DIM], dim: usize, h: f64) -> Self
    where
        F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
    {
        let conn = LeviCivitaConnection::from_metric_fn(metric_fn, point, dim, h);
        let mut dchristoffel = [[[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM];
        for lambda in 0..dim {
            let mut p_plus = *point;
            let mut p_minus = *point;
            p_plus[lambda] += h;
            p_minus[lambda] -= h;
            let conn_plus = LeviCivitaConnection::from_metric_fn(metric_fn, &p_plus, dim, h);
            let conn_minus = LeviCivitaConnection::from_metric_fn(metric_fn, &p_minus, dim, h);
            for sigma in 0..dim {
                for mu in 0..dim {
                    for nu in 0..dim {
                        dchristoffel[lambda][sigma][mu][nu] = (conn_plus.christoffel[sigma][mu]
                            [nu]
                            - conn_minus.christoffel[sigma][mu][nu])
                            / (2.0 * h);
                    }
                }
            }
        }
        let mut components = [[[[0.0_f64; MAX_DIM]; MAX_DIM]; MAX_DIM]; MAX_DIM];
        for rho in 0..dim {
            for sigma in 0..dim {
                for mu in 0..dim {
                    for nu in 0..dim {
                        let mut val =
                            dchristoffel[mu][rho][nu][sigma] - dchristoffel[nu][rho][mu][sigma];
                        for lambda in 0..dim {
                            val += conn.christoffel[rho][mu][lambda]
                                * conn.christoffel[lambda][nu][sigma];
                            val -= conn.christoffel[rho][nu][lambda]
                                * conn.christoffel[lambda][mu][sigma];
                        }
                        components[rho][sigma][mu][nu] = val;
                    }
                }
            }
        }
        Self { components, dim }
    }
    /// Access R^rho_{sigma mu nu}.
    pub fn get(&self, rho: usize, sigma: usize, mu: usize, nu: usize) -> f64 {
        self.components[rho][sigma][mu][nu]
    }
    /// Check the first Bianchi identity: R^rho_{\[sigma mu nu\]} = 0.
    ///
    /// Returns the maximum violation.
    pub fn bianchi_identity_violation(&self) -> f64 {
        let mut max_viol = 0.0_f64;
        for rho in 0..self.dim {
            for sigma in 0..self.dim {
                for mu in 0..self.dim {
                    for nu in 0..self.dim {
                        let cyclic = self.components[rho][sigma][mu][nu]
                            + self.components[rho][mu][nu][sigma]
                            + self.components[rho][nu][sigma][mu];
                        max_viol = max_viol.max(cyclic.abs());
                    }
                }
            }
        }
        max_viol
    }
}
/// Element of se(3): a 6-vector \[ω; v\] representing angular + linear velocity.
#[derive(Debug, Clone, Copy, Default)]
pub struct Se3Algebra {
    /// Angular velocity component ω ∈ R³.
    pub omega: Vec3,
    /// Linear velocity component v ∈ R³.
    pub v: Vec3,
}
impl Se3Algebra {
    /// Create a new se(3) element.
    pub fn new(omega: Vec3, v: Vec3) -> Self {
        Self { omega, v }
    }
    /// Hat operator: convert to 4×4 matrix form.
    pub fn hat(&self) -> Mat4 {
        let k = skew3(self.omega);
        [
            [k[0][0], k[0][1], k[0][2], self.v[0]],
            [k[1][0], k[1][1], k[1][2], self.v[1]],
            [k[2][0], k[2][1], k[2][2], self.v[2]],
            [0.0, 0.0, 0.0, 0.0],
        ]
    }
    /// Norm (weighted) of the twist.
    pub fn norm(&self) -> f64 {
        (dot3(self.omega, self.omega) + dot3(self.v, self.v)).sqrt()
    }
    /// Scale the twist.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            omega: scale3(self.omega, s),
            v: scale3(self.v, s),
        }
    }
}
/// Ricci tensor R_{mu nu} obtained by contracting the Riemann tensor.
///
/// R_{mu nu} = R^lambda_{mu lambda nu}.
#[derive(Debug, Clone)]
pub struct RicciTensor {
    /// Components R_{mu nu}.
    pub components: [[f64; MAX_DIM]; MAX_DIM],
    /// Dimension of the manifold.
    pub dim: usize,
}
impl RicciTensor {
    /// Compute from a Riemann tensor by contraction.
    pub fn from_riemann(riemann: &RiemannTensor) -> Self {
        let dim = riemann.dim;
        let mut components = [[0.0_f64; MAX_DIM]; MAX_DIM];
        for mu in 0..dim {
            for nu in 0..dim {
                let mut val = 0.0;
                for lambda in 0..dim {
                    val += riemann.components[lambda][mu][lambda][nu];
                }
                components[mu][nu] = val;
            }
        }
        Self { components, dim }
    }
    /// Compute directly from a metric function at a point.
    pub fn from_metric_fn<F>(metric_fn: &F, point: &[f64; MAX_DIM], dim: usize, h: f64) -> Self
    where
        F: Fn(&[f64; MAX_DIM]) -> [[f64; MAX_DIM]; MAX_DIM],
    {
        let riemann = RiemannTensor::from_metric_fn(metric_fn, point, dim, h);
        Self::from_riemann(&riemann)
    }
    /// Access R_{mu nu}.
    pub fn get(&self, mu: usize, nu: usize) -> f64 {
        self.components[mu][nu]
    }
    /// Check symmetry: R_{mu nu} should equal R_{nu mu}.
    ///
    /// Returns maximum asymmetry.
    pub fn symmetry_violation(&self) -> f64 {
        let mut max_viol = 0.0_f64;
        for mu in 0..self.dim {
            for nu in 0..self.dim {
                max_viol = max_viol.max((self.components[mu][nu] - self.components[nu][mu]).abs());
            }
        }
        max_viol
    }
}
/// Element of the SE(3) group: a rotation plus translation.
#[derive(Debug, Clone, Copy)]
pub struct Se3 {
    /// Rotation component.
    pub rotation: So3,
    /// Translation component.
    pub translation: Vec3,
}
impl Se3 {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            rotation: So3::identity(),
            translation: [0.0, 0.0, 0.0],
        }
    }
    /// Construct from rotation and translation.
    pub fn from_rt(rotation: So3, translation: Vec3) -> Self {
        Self {
            rotation,
            translation,
        }
    }
    /// Exponential map: exp(ξ̂) where ξ = \[ω; v\] is a se(3) twist.
    pub fn exp(twist: Se3Algebra) -> Self {
        let omega = twist.omega;
        let v = twist.v;
        let theta = norm3(omega);
        let rot = So3::exp(omega);
        if theta < 1e-12 {
            return Self {
                rotation: rot,
                translation: v,
            };
        }
        let axis = scale3(omega, 1.0 / theta);
        let k = skew3(axis);
        let k2 = mat3_mul(k, k);
        let a = (1.0 - theta.cos()) / (theta * theta);
        let b = (theta - theta.sin()) / (theta * theta * theta);
        let i = mat3_identity();
        let v_mat = mat3_add(
            mat3_add(i, mat3_scale(k, a * theta)),
            mat3_scale(k2, b * theta * theta),
        );
        let t = mat3_mul_vec3(v_mat, v);
        Self {
            rotation: rot,
            translation: t,
        }
    }
    /// Logarithm map: returns se(3) twist ξ such that exp(ξ) ≈ self.
    pub fn log(&self) -> Se3Algebra {
        let omega = self.rotation.log();
        let theta = norm3(omega);
        if theta < 1e-12 {
            return Se3Algebra::new(omega, self.translation);
        }
        let k = skew3(scale3(omega, 1.0 / theta));
        let k2 = mat3_mul(k, k);
        let a = 1.0 / theta;
        let b =
            (2.0 * theta.sin() - theta * (1.0 + theta.cos())) / (2.0 * theta * theta * theta.sin());
        let i = mat3_identity();
        let v_inv = mat3_add(
            mat3_add(i, mat3_scale(k, -0.5 * theta)),
            mat3_scale(k2, b * theta * theta),
        );
        let v = mat3_mul_vec3(v_inv, self.translation);
        let _ = a;
        Se3Algebra::new(omega, v)
    }
    /// Compose two SE(3) elements: self · other.
    pub fn compose(&self, other: &Se3) -> Se3 {
        Se3 {
            rotation: self.rotation.compose(&other.rotation),
            translation: add3(self.translation, self.rotation.rotate(other.translation)),
        }
    }
    /// Inverse transform.
    pub fn inverse(&self) -> Se3 {
        let r_inv = self.rotation.inverse();
        Se3 {
            rotation: r_inv,
            translation: scale3(r_inv.rotate(self.translation), -1.0),
        }
    }
    /// Apply transform to a point.
    pub fn transform_point(&self, p: Vec3) -> Vec3 {
        add3(self.rotation.rotate(p), self.translation)
    }
    /// Apply transform to a direction vector (rotation only).
    pub fn transform_direction(&self, d: Vec3) -> Vec3 {
        self.rotation.rotate(d)
    }
    /// Adjoint representation (6×6 matrix in 3+3 block form) — returns (R, p×R; 0, R).
    /// Returned as two Mat3 blocks: (top_left, top_right).
    pub fn adjoint_blocks(&self) -> (Mat3, Mat3) {
        let r = self.rotation.mat;
        let p = skew3(self.translation);
        let pr = mat3_mul(p, r);
        (r, pr)
    }
    /// Geodesic distance to another SE(3) element.
    pub fn dist(&self, other: &Se3) -> f64 {
        let twist = self.inverse().compose(other).log();
        twist.norm()
    }
    /// SE(3) interpolation via Lie algebra.
    pub fn interp(&self, other: &Se3, t: f64) -> Se3 {
        let diff = self.inverse().compose(other);
        let twist = diff.log();
        let step = Se3::exp(twist.scale(t));
        self.compose(&step)
    }
    /// Convert to 4×4 homogeneous matrix.
    pub fn to_mat4(&self) -> Mat4 {
        let r = self.rotation.mat;
        let t = self.translation;
        [
            [r[0][0], r[0][1], r[0][2], t[0]],
            [r[1][0], r[1][1], r[1][2], t[1]],
            [r[2][0], r[2][1], r[2][2], t[2]],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}
/// A 1-form (covector field) on R³, given as coefficient functions at a point.
#[derive(Debug, Clone, Copy)]
pub struct OneForm {
    /// Components \[a₁, a₂, a₃\] such that ω = a₁dx + a₂dy + a₃dz.
    pub components: Vec3,
}
impl OneForm {
    /// Create a 1-form from components.
    pub fn new(a: Vec3) -> Self {
        Self { components: a }
    }
    /// Evaluate on a tangent vector v.
    pub fn evaluate(&self, v: Vec3) -> f64 {
        dot3(self.components, v)
    }
    /// Exterior derivative: d(α) = 0 for a constant 1-form.
    /// Returns the 2-form ∂α_j/∂x_i - ∂α_i/∂x_j (approximated at a point).
    pub fn exterior_derivative_constant(&self) -> TwoForm {
        TwoForm::zero()
    }
}
