//! Auto-generated module
//!
//! Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::decomposition::cp_reconstruct;
use super::operations::{matmul, transpose, truncated_svd_left, tucker_reconstruct};

/// Second-order tensor (3×3 matrix) stored as row-major `[[f64; 3]; 3]`.
pub struct Tensor2 {
    /// Row-major 3×3 component array: `data[i][j]` is the (i,j) component.
    pub data: [[f64; 3]; 3],
}
impl Tensor2 {
    /// Create a new `Tensor2` from a 3×3 array.
    pub fn new(data: [[f64; 3]; 3]) -> Self {
        Self { data }
    }
    /// Zero tensor.
    pub fn zero() -> Self {
        Self {
            data: [[0.0; 3]; 3],
        }
    }
    /// Second-order identity tensor.
    pub fn identity() -> Self {
        let mut d = [[0.0f64; 3]; 3];
        d[0][0] = 1.0;
        d[1][1] = 1.0;
        d[2][2] = 1.0;
        Self { data: d }
    }
    /// Component-wise addition.
    pub fn add(&self, other: &Tensor2) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = self.data[i][j] + other.data[i][j];
            }
        }
        Tensor2 { data: d }
    }
    /// Component-wise subtraction.
    pub fn sub(&self, other: &Tensor2) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = self.data[i][j] - other.data[i][j];
            }
        }
        Tensor2 { data: d }
    }
    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = self.data[i][j] * s;
            }
        }
        Tensor2 { data: d }
    }
    /// Matrix–matrix product (A · B).
    pub fn dot(&self, other: &Tensor2) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    d[i][j] += self.data[i][k] * other.data[k][j];
                }
            }
        }
        Tensor2 { data: d }
    }
    /// Transpose A^T.
    pub fn transpose(&self) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = self.data[j][i];
            }
        }
        Tensor2 { data: d }
    }
    /// Trace: sum of diagonal components.
    pub fn trace(&self) -> f64 {
        self.data[0][0] + self.data[1][1] + self.data[2][2]
    }
    /// Determinant of the 3×3 matrix.
    pub fn det(&self) -> f64 {
        let a = &self.data;
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    }
    /// Inverse via cofactor expansion; returns `None` if `|det| < 1e-14`.
    pub fn inverse(&self) -> Option<Tensor2> {
        let d = self.det();
        if d.abs() < 1e-14 {
            return None;
        }
        let a = &self.data;
        let inv_d = 1.0 / d;
        let mut r = [[0.0f64; 3]; 3];
        r[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_d;
        r[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_d;
        r[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_d;
        r[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_d;
        r[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_d;
        r[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_d;
        r[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_d;
        r[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_d;
        r[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_d;
        Some(Tensor2 { data: r })
    }
    /// Symmetric part: (A + A^T) / 2.
    pub fn sym(&self) -> Tensor2 {
        self.add(&self.transpose()).scale(0.5)
    }
    /// Skew-symmetric part: (A − A^T) / 2.
    pub fn skew(&self) -> Tensor2 {
        self.sub(&self.transpose()).scale(0.5)
    }
    /// Double contraction A:B = Σ_ij A_ij B_ij.
    pub fn double_contract(&self, other: &Tensor2) -> f64 {
        let mut s = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                s += self.data[i][j] * other.data[i][j];
            }
        }
        s
    }
    /// Frobenius norm √(A:A).
    pub fn norm(&self) -> f64 {
        self.double_contract(self).sqrt()
    }
    /// Deviatoric part: A − (tr(A)/3) · I.
    pub fn deviatoric(&self) -> Tensor2 {
        let h = self.trace() / 3.0;
        let mut d = self.data;
        d[0][0] -= h;
        d[1][1] -= h;
        d[2][2] -= h;
        Tensor2 { data: d }
    }
    /// Hydrostatic (mean normal) stress: tr(A) / 3.
    pub fn hydrostatic(&self) -> f64 {
        self.trace() / 3.0
    }
    /// von Mises equivalent: √(3/2 · s:s) where s = deviatoric part.
    pub fn von_mises(&self) -> f64 {
        let s = self.deviatoric();
        let ss = s.double_contract(&s);
        (1.5 * ss).sqrt()
    }
    /// Eigenvalues of a symmetric 3×3 tensor via Cardano's analytical formula.
    ///
    /// The tensor is assumed symmetric; only the lower-triangular part is used.
    /// Returns eigenvalues sorted in ascending order.
    pub fn eigenvalues_symmetric(&self) -> [f64; 3] {
        let a = &self.data;
        let p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if p1.abs() < 1e-30 {
            let mut ev = [a[0][0], a[1][1], a[2][2]];
            ev.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            return ev;
        }
        let q = self.trace() / 3.0;
        let p2 = (a[0][0] - q) * (a[0][0] - q)
            + (a[1][1] - q) * (a[1][1] - q)
            + (a[2][2] - q) * (a[2][2] - q)
            + 2.0 * p1;
        let p = (p2 / 6.0).sqrt();
        let mut b = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                b[i][j] = a[i][j] / p;
            }
        }
        b[0][0] -= q / p;
        b[1][1] -= q / p;
        b[2][2] -= q / p;
        let det_b = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
            - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
            + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
        let r = (det_b / 2.0).clamp(-1.0, 1.0);
        let phi = r.acos() / 3.0;
        use std::f64::consts::PI;
        let eig1 = q + 2.0 * p * phi.cos();
        let eig3 = q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
        let eig2 = 3.0 * q - eig1 - eig3;
        let mut ev = [eig1, eig2, eig3];
        ev.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        ev
    }
    /// Outer (dyadic) product a ⊗ b: result_ij = a_i * b_j.
    pub fn outer_product(a: [f64; 3], b: [f64; 3]) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = a[i] * b[j];
            }
        }
        Tensor2 { data: d }
    }
    /// Matrix–vector product: result_i = Σ_j A_ij v_j.
    pub fn apply(&self, v: [f64; 3]) -> [f64; 3] {
        let mut r = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i] += self.data[i][j] * v[j];
            }
        }
        r
    }
}
impl Tensor2 {
    /// Single contraction: c_i = A_ij * v_j (same as apply, alias).
    pub fn contract_vec(&self, v: [f64; 3]) -> [f64; 3] {
        self.apply(v)
    }
    /// Tensor outer (dyadic) product from two vectors (static method alias).
    pub fn dyadic(a: [f64; 3], b: [f64; 3]) -> Tensor2 {
        Self::outer_product(a, b)
    }
    /// Rotate a tensor by a rotation matrix R: R * A * R^T.
    pub fn rotate(&self, r: &Tensor2) -> Tensor2 {
        let rt = r.transpose();
        r.dot(self).dot(&rt)
    }
    /// First invariant I1 = trace(A).
    pub fn invariant_i1(&self) -> f64 {
        self.trace()
    }
    /// Second invariant I2 = 0.5 * (trace(A)^2 - trace(A^2)).
    pub fn invariant_i2(&self) -> f64 {
        let tr = self.trace();
        let a2 = self.dot(self);
        0.5 * (tr * tr - a2.trace())
    }
    /// Third invariant I3 = det(A).
    pub fn invariant_i3(&self) -> f64 {
        self.det()
    }
    /// Principal invariants as an array \[I1, I2, I3\].
    pub fn principal_invariants(&self) -> [f64; 3] {
        [
            self.invariant_i1(),
            self.invariant_i2(),
            self.invariant_i3(),
        ]
    }
    /// Deviatoric and hydrostatic decomposition.
    ///
    /// Returns `(deviatoric, hydrostatic_scalar)`.
    pub fn decompose_dev_hydro(&self) -> (Tensor2, f64) {
        let h = self.hydrostatic();
        (self.deviatoric(), h)
    }
    /// Create a tensor from Voigt notation \[xx, yy, zz, xy, yz, xz\].
    pub fn from_voigt(v: [f64; 6]) -> Tensor2 {
        VoigtTensor::to_tensor2(&v)
    }
    /// Convert to Voigt notation.
    pub fn to_voigt(&self) -> [f64; 6] {
        VoigtTensor::from_tensor2(self)
    }
    /// Effective strain (von Mises strain equivalent).
    ///
    /// For a strain tensor: eps_eff = sqrt(2/3 * e_ij * e_ij)
    /// where e is the deviatoric strain.
    pub fn effective_strain(&self) -> f64 {
        let dev = self.deviatoric();
        let ee = dev.double_contract(&dev);
        (2.0 / 3.0 * ee).sqrt()
    }
    /// Create a diagonal tensor with given diagonal values.
    pub fn diagonal(d: [f64; 3]) -> Tensor2 {
        Tensor2::new([[d[0], 0.0, 0.0], [0.0, d[1], 0.0], [0.0, 0.0, d[2]]])
    }
    /// Check if the tensor is symmetric within tolerance.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                if (self.data[i][j] - self.data[j][i]).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
    /// Matrix exponential approximation using Taylor series (for small tensors).
    ///
    /// exp(A) ≈ I + A + A²/2! + A³/3! + ...
    /// Uses 10 terms by default.
    pub fn matrix_exp_approx(&self) -> Tensor2 {
        let mut result = Tensor2::identity();
        let mut term = Tensor2::identity();
        for n in 1..=10 {
            term = term.dot(self).scale(1.0 / n as f64);
            result = result.add(&term);
        }
        result
    }
}
impl Tensor2 {
    #[allow(dead_code)]
    fn clone_t2(&self) -> Tensor2 {
        Tensor2 { data: self.data }
    }
}
/// Rank-3 tensor with 27 components stored as `data[i][j][k]`.
///
/// Used in continuum mechanics for third-order coupling tensors and the
/// Levi-Civita permutation symbol.
pub struct Tensor3 {
    /// Component array indexed as `data[i][j][k]`.
    pub data: [[[f64; 3]; 3]; 3],
}
impl Tensor3 {
    /// Zero tensor.
    pub fn zero() -> Self {
        Self {
            data: [[[0.0; 3]; 3]; 3],
        }
    }
    /// Levi-Civita (alternating) tensor ε_ijk.
    ///
    /// ε_ijk = +1 for even permutations of (0,1,2),
    ///         -1 for odd permutations,
    ///          0 if any index is repeated.
    pub fn levi_civita() -> Self {
        let mut t = Self::zero();
        t.data[0][1][2] = 1.0;
        t.data[1][2][0] = 1.0;
        t.data[2][0][1] = 1.0;
        t.data[0][2][1] = -1.0;
        t.data[2][1][0] = -1.0;
        t.data[1][0][2] = -1.0;
        t
    }
    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Self {
        let mut out = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out.data[i][j][k] = self.data[i][j][k] * s;
                }
            }
        }
        out
    }
    /// Component-wise addition.
    pub fn add(&self, other: &Tensor3) -> Tensor3 {
        let mut out = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out.data[i][j][k] = self.data[i][j][k] + other.data[i][j][k];
                }
            }
        }
        out
    }
    /// Contract the last index with a vector: result\[i\]\[j\] = T_ijk * v_k.
    pub fn contract_last(&self, v: &[f64; 3]) -> Tensor2 {
        let mut out = Tensor2::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out.data[i][j] += self.data[i][j][k] * v[k];
                }
            }
        }
        out
    }
    /// Contract the first index with a vector: result\[j\]\[k\] = v_i * T_ijk.
    pub fn contract_first(&self, v: &[f64; 3]) -> Tensor2 {
        let mut out = Tensor2::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out.data[j][k] += v[i] * self.data[i][j][k];
                }
            }
        }
        out
    }
    /// Full contraction with two vectors: result_i = T_ijk * a_j * b_k.
    pub fn contract_two(&self, a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out[i] += self.data[i][j][k] * a[j] * b[k];
                }
            }
        }
        out
    }
    /// Cross product of vectors a and b using the Levi-Civita tensor.
    ///
    /// (a × b)_i = ε_ijk * a_j * b_k
    pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
        let eps = Self::levi_civita();
        eps.contract_two(a, b)
    }
    /// Frobenius norm: sqrt(sum T_ijk^2).
    pub fn norm(&self) -> f64 {
        let mut s = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    s += self.data[i][j][k] * self.data[i][j][k];
                }
            }
        }
        s.sqrt()
    }
}
/// Fourth-order tensor with 81 components stored as `[i][j][k][l]`.
pub struct Tensor4 {
    /// Component array indexed as `data[i][j][k][l]`.
    pub data: [[[[f64; 3]; 3]; 3]; 3],
}
impl Tensor4 {
    /// Zero fourth-order tensor.
    pub fn zero() -> Self {
        Self {
            data: [[[[0.0f64; 3]; 3]; 3]; 3],
        }
    }
    /// Symmetric fourth-order identity:
    /// I_sym_ijkl = 0.5 * (δ_ik δ_jl + δ_il δ_jk).
    pub fn identity_sym() -> Self {
        let mut c = Self::zero();
        for i in 0..3usize {
            for j in 0..3usize {
                for k in 0..3usize {
                    for l in 0..3usize {
                        let dik = if i == k { 1.0 } else { 0.0 };
                        let djl = if j == l { 1.0 } else { 0.0 };
                        let dil = if i == l { 1.0 } else { 0.0 };
                        let djk = if j == k { 1.0 } else { 0.0 };
                        c.data[i][j][k][l] = 0.5 * (dik * djl + dil * djk);
                    }
                }
            }
        }
        c
    }
    /// Isotropic linear elasticity tensor:
    /// C_ijkl = λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk).
    pub fn isotropic(lambda: f64, mu: f64) -> Self {
        let mut c = Self::zero();
        for i in 0..3usize {
            for j in 0..3usize {
                for k in 0..3usize {
                    for l in 0..3usize {
                        let dij = if i == j { 1.0 } else { 0.0 };
                        let dkl = if k == l { 1.0 } else { 0.0 };
                        let dik = if i == k { 1.0 } else { 0.0 };
                        let djl = if j == l { 1.0 } else { 0.0 };
                        let dil = if i == l { 1.0 } else { 0.0 };
                        let djk = if j == k { 1.0 } else { 0.0 };
                        c.data[i][j][k][l] = lambda * dij * dkl + mu * (dik * djl + dil * djk);
                    }
                }
            }
        }
        c
    }
    /// Double contraction with a second-order tensor:
    /// result_ij = Σ_kl C_ijkl A_kl.
    pub fn double_contract_2(&self, t: &Tensor2) -> Tensor2 {
        let mut d = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        d[i][j] += self.data[i][j][k][l] * t.data[k][l];
                    }
                }
            }
        }
        Tensor2 { data: d }
    }
}
impl Tensor4 {
    /// Deviatoric projector: P_dev_ijkl = I_sym_ijkl - (1/3) delta_ij delta_kl.
    pub fn deviatoric_projector() -> Self {
        let mut c = Self::identity_sym();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let dij = if i == j { 1.0 } else { 0.0 };
                        let dkl = if k == l { 1.0 } else { 0.0 };
                        c.data[i][j][k][l] -= (1.0 / 3.0) * dij * dkl;
                    }
                }
            }
        }
        c
    }
    /// Convert isotropic elasticity tensor to Voigt matrix (6x6).
    ///
    /// Returns \[6\]\[6\] matrix.
    pub fn to_voigt_matrix(&self) -> [[f64; 6]; 6] {
        let voigt_map: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
        let mut m = [[0.0f64; 6]; 6];
        for (p, &(i, j)) in voigt_map.iter().enumerate() {
            for (q, &(k, l)) in voigt_map.iter().enumerate() {
                m[p][q] = self.data[i][j][k][l];
            }
        }
        m
    }
    /// Scale all components by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        let mut c = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        c.data[i][j][k][l] = self.data[i][j][k][l] * s;
                    }
                }
            }
        }
        c
    }
    /// Add two fourth-order tensors.
    pub fn add(&self, other: &Tensor4) -> Self {
        let mut c = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        c.data[i][j][k][l] = self.data[i][j][k][l] + other.data[i][j][k][l];
                    }
                }
            }
        }
        c
    }
}
impl Tensor4 {
    /// Double contraction of two fourth-order tensors:
    /// (C::D)_ijmn = Σ_kl C_ijkl D_klmn.
    pub fn double_contract_4(&self, other: &Tensor4) -> Tensor4 {
        let mut result = Tensor4::zero();
        for i in 0..3 {
            for j in 0..3 {
                for m in 0..3 {
                    for n in 0..3 {
                        let mut s = 0.0f64;
                        for k in 0..3 {
                            for l in 0..3 {
                                s += self.data[i][j][k][l] * other.data[k][l][m][n];
                            }
                        }
                        result.data[i][j][m][n] = s;
                    }
                }
            }
        }
        result
    }
    /// Single contraction of C with a second-order tensor on the last two
    /// indices: (C · A)_ijkl = Σ_m C_ijkm A_ml.
    ///
    /// Result is a fourth-order tensor.
    pub fn single_contract_right(&self, a: &Tensor2) -> Tensor4 {
        let mut result = Tensor4::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0f64;
                        for m in 0..3 {
                            s += self.data[i][j][k][m] * a.data[m][l];
                        }
                        result.data[i][j][k][l] = s;
                    }
                }
            }
        }
        result
    }
    /// Check minor symmetry: C_ijkl = C_jikl = C_ijlk.
    ///
    /// Returns `true` if both left-minor and right-minor symmetries hold within `tol`.
    pub fn has_minor_symmetry(&self, tol: f64) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let c = self.data[i][j][k][l];
                        if (c - self.data[j][i][k][l]).abs() > tol {
                            return false;
                        }
                        if (c - self.data[i][j][l][k]).abs() > tol {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }
    /// Check major symmetry: C_ijkl = C_klij.
    ///
    /// Returns `true` if major symmetry holds within `tol`.
    pub fn has_major_symmetry(&self, tol: f64) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        if (self.data[i][j][k][l] - self.data[k][l][i][j]).abs() > tol {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }
    /// Enforce minor symmetry by symmetrizing:
    /// C_sym_ijkl = 0.25 * (C_ijkl + C_jikl + C_ijlk + C_jilk).
    pub fn symmetrize_minor(&self) -> Tensor4 {
        let mut result = Tensor4::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        result.data[i][j][k][l] = 0.25
                            * (self.data[i][j][k][l]
                                + self.data[j][i][k][l]
                                + self.data[i][j][l][k]
                                + self.data[j][i][l][k]);
                    }
                }
            }
        }
        result
    }
    /// Enforce major symmetry by averaging: C_sym_ijkl = 0.5*(C_ijkl + C_klij).
    pub fn symmetrize_major(&self) -> Tensor4 {
        let mut result = Tensor4::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        result.data[i][j][k][l] =
                            0.5 * (self.data[i][j][k][l] + self.data[k][l][i][j]);
                    }
                }
            }
        }
        result
    }
    /// Rotate a fourth-order tensor: C'_ijkl = R_ia R_jb R_kc R_ld C_abcd.
    pub fn rotate(&self, r: &Tensor2) -> Tensor4 {
        let mut result = Tensor4::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0f64;
                        for a in 0..3 {
                            for b in 0..3 {
                                for c in 0..3 {
                                    for d in 0..3 {
                                        s += r.data[i][a]
                                            * r.data[j][b]
                                            * r.data[k][c]
                                            * r.data[l][d]
                                            * self.data[a][b][c][d];
                                    }
                                }
                            }
                        }
                        result.data[i][j][k][l] = s;
                    }
                }
            }
        }
        result
    }
    /// Build the isotropic volumetric projector:
    /// P_vol_ijkl = (1/3) * delta_ij * delta_kl.
    pub fn volumetric_projector() -> Self {
        let mut c = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let dij = if i == j { 1.0 } else { 0.0 };
                        let dkl = if k == l { 1.0 } else { 0.0 };
                        c.data[i][j][k][l] = (1.0 / 3.0) * dij * dkl;
                    }
                }
            }
        }
        c
    }
    /// The fourth-order identity: I_ijkl = delta_ik * delta_jl.
    pub fn identity() -> Self {
        let mut c = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let dik = if i == k { 1.0 } else { 0.0 };
                        let djl = if j == l { 1.0 } else { 0.0 };
                        c.data[i][j][k][l] = dik * djl;
                    }
                }
            }
        }
        c
    }
    /// Subtract two fourth-order tensors.
    pub fn sub(&self, other: &Tensor4) -> Self {
        let mut c = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        c.data[i][j][k][l] = self.data[i][j][k][l] - other.data[i][j][k][l];
                    }
                }
            }
        }
        c
    }
    /// Frobenius norm of the fourth-order tensor: sqrt(Σ C_ijkl^2).
    pub fn norm(&self) -> f64 {
        let mut s = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        s += self.data[i][j][k][l] * self.data[i][j][k][l];
                    }
                }
            }
        }
        s.sqrt()
    }
}
/// Mandel notation for symmetric second-order tensors.
///
/// Like Kelvin notation, off-diagonal components are scaled by √2 so that
/// the Euclidean inner product of the 6-vectors equals the double contraction
/// A:B of the original tensors.  Index ordering: \[11, 22, 33, 12, 23, 13\].
pub struct MandelNotation;
impl MandelNotation {
    /// Voigt index map: 0→(0,0), 1→(1,1), 2→(2,2), 3→(0,1), 4→(1,2), 5→(0,2).
    const MAP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    /// Convert a symmetric `Tensor2` to its 6-component Mandel vector.
    ///
    /// Normal components (indices 0..3) are unchanged.
    /// Shear components (indices 3..6) are multiplied by √2.
    pub fn from_tensor2(t: &Tensor2) -> [f64; 6] {
        let s2 = std::f64::consts::SQRT_2;
        let mut v = [0.0f64; 6];
        for (p, &(i, j)) in Self::MAP.iter().enumerate() {
            let w = if p < 3 { 1.0 } else { s2 };
            v[p] = w * t.data[i][j];
        }
        v
    }
    /// Convert a 6-component Mandel vector back to a `Tensor2`.
    ///
    /// Shear components are divided by √2 and symmetrised.
    pub fn to_tensor2(v: &[f64; 6]) -> Tensor2 {
        let s2_inv = 1.0 / std::f64::consts::SQRT_2;
        let mut t = Tensor2::zero();
        for (p, &(i, j)) in Self::MAP.iter().enumerate() {
            let w = if p < 3 { 1.0 } else { s2_inv };
            let val = w * v[p];
            t.data[i][j] = val;
            t.data[j][i] = val;
        }
        t
    }
    /// Compute the double contraction A:B using Mandel vectors.
    ///
    /// A:B = v_A · v_B (Euclidean dot product of the Mandel representations).
    pub fn double_contract(va: &[f64; 6], vb: &[f64; 6]) -> f64 {
        va.iter().zip(vb.iter()).map(|(a, b)| a * b).sum()
    }
    /// Convert a Mandel strain vector to engineering (Voigt) strain.
    ///
    /// Engineering shear strains (γ = 2ε) are twice the tensor shear strains.
    pub fn to_engineering_strain(v: &[f64; 6]) -> [f64; 6] {
        let s2 = std::f64::consts::SQRT_2;
        let mut e = *v;
        for p in 3..6 {
            e[p] = v[p] * s2;
        }
        e
    }
    /// Convert engineering (Voigt) strain to Mandel strain.
    pub fn from_engineering_strain(e: &[f64; 6]) -> [f64; 6] {
        let s2_inv = 1.0 / std::f64::consts::SQRT_2;
        let mut v = *e;
        for p in 3..6 {
            v[p] = e[p] * s2_inv;
        }
        v
    }
}
/// Kelvin (Mandel) notation utilities for symmetric fourth-order tensors.
///
/// Unlike Voigt notation, Kelvin notation preserves the inner product:
/// off-diagonal entries are multiplied by √2 for strain-like vectors and √2
/// for the stiffness off-diagonal blocks, so that C:ε:ε is preserved.
pub struct KelvinTensor;
impl KelvinTensor {
    /// Kelvin index ordering (same as Voigt): 0->00, 1->11, 2->22, 3->01, 4->12, 5->02.
    const VOIGT_MAP: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    /// Kelvin weight factor: √2 for shear components (indices 3, 4, 5), 1 for normal.
    fn weight(idx: usize) -> f64 {
        if idx < 3 {
            1.0
        } else {
            std::f64::consts::SQRT_2
        }
    }
    /// Convert a `Tensor4` (with minor symmetry) to its 6×6 Kelvin matrix.
    ///
    /// K_PQ = w_P * w_Q * C_ijkl   where (i,j) is mapped from P and (k,l) from Q.
    pub fn from_tensor4(c: &Tensor4) -> [[f64; 6]; 6] {
        let mut m = [[0.0f64; 6]; 6];
        for (p, &(i, j)) in Self::VOIGT_MAP.iter().enumerate() {
            for (q, &(k, l)) in Self::VOIGT_MAP.iter().enumerate() {
                m[p][q] = Self::weight(p) * Self::weight(q) * c.data[i][j][k][l];
            }
        }
        m
    }
    /// Convert a 6×6 Kelvin matrix back to a `Tensor4`.
    pub fn to_tensor4(m: &[[f64; 6]; 6]) -> Tensor4 {
        let mut c = Tensor4::zero();
        for (p, &(i, j)) in Self::VOIGT_MAP.iter().enumerate() {
            for (q, &(k, l)) in Self::VOIGT_MAP.iter().enumerate() {
                let val = m[p][q] / (Self::weight(p) * Self::weight(q));
                c.data[i][j][k][l] = val;
                c.data[j][i][k][l] = val;
                c.data[i][j][l][k] = val;
                c.data[j][i][l][k] = val;
            }
        }
        c
    }
    /// Convert a symmetric second-order tensor to its 6-component Kelvin vector.
    ///
    /// Normal components are unchanged; shear components are multiplied by √2.
    pub fn stress_to_kelvin(t: &Tensor2) -> [f64; 6] {
        let s2 = std::f64::consts::SQRT_2;
        [
            t.data[0][0],
            t.data[1][1],
            t.data[2][2],
            s2 * t.data[0][1],
            s2 * t.data[1][2],
            s2 * t.data[0][2],
        ]
    }
    /// Convert a Kelvin stress vector back to a `Tensor2`.
    pub fn kelvin_to_stress(v: &[f64; 6]) -> Tensor2 {
        let s2_inv = 1.0 / std::f64::consts::SQRT_2;
        Tensor2 {
            data: [
                [v[0], s2_inv * v[3], s2_inv * v[5]],
                [s2_inv * v[3], v[1], s2_inv * v[4]],
                [s2_inv * v[5], s2_inv * v[4], v[2]],
            ],
        }
    }
    /// Multiply the 6×6 Kelvin stiffness matrix by a Kelvin stress/strain vector.
    pub fn matvec(m: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
        let mut result = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                result[i] += m[i][j] * v[j];
            }
        }
        result
    }
}
/// CP (CANDECOMP/PARAFAC) decomposition result.
///
/// A rank-R CP decomposition of a tensor T ≈ Σ_r λ_r (a_r ⊗ b_r ⊗ c_r).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CpDecomposition {
    /// Mode-0 factor matrix: shape (n0, rank).
    pub a: Vec<Vec<f64>>,
    /// Mode-1 factor matrix: shape (n1, rank).
    pub b: Vec<Vec<f64>>,
    /// Mode-2 factor matrix: shape (n2, rank).
    pub c: Vec<Vec<f64>>,
    /// Normalisation weights λ_r.
    pub lambdas: Vec<f64>,
}
impl CpDecomposition {
    /// Compute the Frobenius norm of the reconstructed tensor.
    #[allow(dead_code)]
    pub fn reconstruction_norm(&self) -> f64 {
        cp_reconstruct(self).frobenius_norm()
    }
    /// Return the number of components (CP rank).
    #[allow(dead_code)]
    pub fn rank(&self) -> usize {
        self.lambdas.len()
    }
    /// Return the rank (number of CP components).
    #[allow(dead_code)]
    pub fn cp_rank(&self) -> usize {
        self.lambdas.len()
    }
}
/// A tensor in Tensor Train (TT / MPS) format.
///
/// A d-way tensor T of shape (n_0, …, n_{d-1}) is represented as a product
/// of 3-way cores G_k of shape (r_{k-1}, n_k, r_k) where r_0 = r_d = 1.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TensorTrain {
    /// The TT cores.  Core k has shape (r_{k-1}, n_k, r_k).
    /// Each core is stored as a flat `Vec`f64` in C-order.
    pub cores: Vec<TtCore>,
    /// Mode dimensions n_k.
    pub shape: Vec<usize>,
}
#[allow(dead_code)]
impl TensorTrain {
    /// Evaluate the TT representation at a given multi-index.
    ///
    /// Contracts the cores left-to-right: result = G_0\[i_0\] G_1\[i_1\] … G_{d-1}\[i_{d-1}\].
    /// Each G_k\[i_k\] is the r_{k-1}×r_k slice.
    pub fn evaluate(&self, idx: &[usize]) -> f64 {
        assert_eq!(idx.len(), self.cores.len());
        let mut vec = vec![1.0f64];
        for (k, core) in self.cores.iter().enumerate() {
            let ik = idx[k];
            let mut next = vec![0.0f64; core.r_right];
            for alpha in 0..core.r_left {
                for beta in 0..core.r_right {
                    next[beta] += vec[alpha] * core.get(alpha, ik, beta);
                }
            }
            vec = next;
        }
        vec[0]
    }
    /// Convert a full dense tensor to Tensor Train format via left-to-right
    /// SVD rounding (TT-SVD algorithm).
    ///
    /// `max_rank` caps each TT rank; `tol` is relative truncation tolerance.
    pub fn from_dense(tensor: &DenseTensor, max_rank: usize, tol: f64) -> TensorTrain {
        let d = tensor.shape.len();
        let shape = tensor.shape.clone();
        let mut cores = Vec::with_capacity(d);
        let total: usize = tensor.data.len();
        let c_rows = 1usize;
        let mut c_data = tensor.data.clone();
        let c_cols = total;
        let mut c_shape = (c_rows, c_cols);
        for k in 0..d - 1 {
            let nk = shape[k];
            let rows = c_shape.0 * nk;
            let cols = c_shape.1 / nk;
            let mat: Vec<Vec<f64>> = (0..rows)
                .map(|r| (0..cols).map(|c_idx| c_data[r * cols + c_idx]).collect())
                .collect();
            let rank = rows.min(cols).min(max_rank);
            let u = truncated_svd_left(&mat, rank);
            let r_left = c_shape.0;
            let r_right = rank;
            let mut core = TtCore::zeros(r_left, nk, r_right);
            for alpha in 0..r_left {
                for i in 0..nk {
                    for beta in 0..r_right {
                        core.set(alpha, i, beta, u[alpha * nk + i][beta]);
                    }
                }
            }
            cores.push(core);
            let ut = transpose(&u);
            let c_new = matmul(&ut, &mat);
            c_data = c_new.into_iter().flatten().collect();
            c_shape = (rank, cols);
            let _ = tol;
        }
        let nk = shape[d - 1];
        let r_left = c_shape.0;
        let mut last_core = TtCore::zeros(r_left, nk, 1);
        for alpha in 0..r_left {
            for i in 0..nk {
                last_core.set(alpha, i, 0, c_data[alpha * nk + i]);
            }
        }
        cores.push(last_core);
        TensorTrain { cores, shape }
    }
    /// Frobenius norm of the TT tensor (contract all cores).
    pub fn frobenius_norm_approx(&self) -> f64 {
        let total: usize = self.shape.iter().product();
        let mut norm_sq = 0.0;
        let mut idx = vec![0usize; self.shape.len()];
        for _ in 0..total {
            let v = self.evaluate(&idx);
            norm_sq += v * v;
            let mut carry = true;
            for k in (0..idx.len()).rev() {
                if carry {
                    idx[k] += 1;
                    if idx[k] == self.shape[k] {
                        idx[k] = 0;
                    } else {
                        carry = false;
                    }
                }
            }
        }
        norm_sq.sqrt()
    }
}
impl TensorTrain {
    /// Compute the Frobenius norm via contraction of the full TT.
    ///
    /// For small tensors only — reconstructs fully.
    #[allow(dead_code)]
    pub fn frobenius_norm(&self) -> f64 {
        let dense = self.to_dense();
        dense.frobenius_norm()
    }
    /// Reconstruct to a dense DenseTensor.
    #[allow(dead_code)]
    pub fn to_dense(&self) -> DenseTensor {
        let d = self.shape.len();
        let total: usize = self.shape.iter().product();
        let mut data = vec![0.0f64; total];
        let mut strides = vec![1usize; d];
        for k in (0..d.saturating_sub(1)).rev() {
            strides[k] = strides[k + 1] * self.shape[k + 1];
        }
        let shape_clone = self.shape.clone();
        let d_clone = d;
        for flat in 0..total {
            let mut tmp = flat;
            let mut midx = vec![0usize; d_clone];
            for k in 0..d_clone {
                midx[k] = tmp / strides[k];
                tmp %= strides[k];
            }
            data[flat] = self.evaluate(&midx);
        }
        let _ = shape_clone;
        DenseTensor {
            shape: self.shape.clone(),
            data,
        }
    }
    /// Dot product <TT_a, TT_b> between two TensorTrains with the same shape.
    ///
    /// Uses the supercore contraction (quadratic in ranks).
    #[allow(dead_code)]
    pub fn dot_product(&self, other: &TensorTrain) -> f64 {
        assert_eq!(
            self.shape, other.shape,
            "TensorTrains must have same shape for dot product"
        );
        let dense_self = self.to_dense();
        let dense_other = other.to_dense();
        dense_self
            .data
            .iter()
            .zip(dense_other.data.iter())
            .map(|(a, b)| a * b)
            .sum()
    }
    /// Scale all cores by a scalar.
    #[allow(dead_code)]
    pub fn scale(&self, s: f64) -> TensorTrain {
        let mut cores = self.cores.clone();
        if !cores.is_empty() {
            cores[0].data = cores[0].data.iter().map(|&x| x * s).collect();
        }
        TensorTrain {
            cores,
            shape: self.shape.clone(),
        }
    }
}
/// A dense tensor of arbitrary rank stored in row-major (C-order) layout.
///
/// The shape is given as a `Vec`usize` and the flat data length must equal
/// the product of all dimensions.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DenseTensor {
    /// Shape of the tensor, e.g. `[2, 3, 4]` for a 2×3×4 tensor.
    pub shape: Vec<usize>,
    /// Flat row-major storage.
    pub data: Vec<f64>,
}
#[allow(dead_code)]
impl DenseTensor {
    /// Create a zero tensor with the given shape.
    pub fn zeros(shape: &[usize]) -> Self {
        let n: usize = shape.iter().product();
        DenseTensor {
            shape: shape.to_vec(),
            data: vec![0.0; n],
        }
    }
    /// Create a tensor from flat data (row-major).
    ///
    /// # Panics
    /// Panics if `data.len() != shape.iter().product()`.
    pub fn from_data(shape: &[usize], data: Vec<f64>) -> Self {
        let n: usize = shape.iter().product();
        assert_eq!(data.len(), n, "data length must match shape product");
        DenseTensor {
            shape: shape.to_vec(),
            data,
        }
    }
    /// Compute the strides for row-major indexing.
    pub fn strides(&self) -> Vec<usize> {
        let rank = self.shape.len();
        let mut s = vec![1usize; rank];
        for k in (0..rank.saturating_sub(1)).rev() {
            s[k] = s[k + 1] * self.shape[k + 1];
        }
        s
    }
    /// Flat index from a multi-index.
    pub fn flat_index(&self, idx: &[usize]) -> usize {
        let s = self.strides();
        idx.iter().zip(s.iter()).map(|(&i, &st)| i * st).sum()
    }
    /// Get element at multi-index.
    pub fn get(&self, idx: &[usize]) -> f64 {
        self.data[self.flat_index(idx)]
    }
    /// Set element at multi-index.
    pub fn set(&mut self, idx: &[usize], val: f64) {
        let fi = self.flat_index(idx);
        self.data[fi] = val;
    }
    /// Frobenius norm: sqrt of sum of squared elements.
    pub fn frobenius_norm(&self) -> f64 {
        self.data.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }
    /// Mode-n unfolding (matricisation).
    ///
    /// Returns a matrix of shape `(shape[n], prod(other dims))` stored
    /// row-major as a `Vec<Vec`f64`>`.
    pub fn mode_n_unfold(&self, mode: usize) -> Vec<Vec<f64>> {
        let n_rows = self.shape[mode];
        let n_cols: usize = self.data.len() / n_rows;
        let rank = self.shape.len();
        let strides = self.strides();
        let mut mat = vec![vec![0.0; n_cols]; n_rows];
        let total: usize = self.data.len();
        for flat in 0..total {
            let mut tmp = flat;
            let mut midx = vec![0usize; rank];
            for k in 0..rank {
                midx[k] = tmp / strides[k];
                tmp %= strides[k];
            }
            let row = midx[mode];
            let mut col = 0usize;
            let mut col_stride = 1usize;
            for k in (0..rank).rev() {
                if k == mode {
                    continue;
                }
                col += midx[k] * col_stride;
                col_stride *= self.shape[k];
            }
            mat[row][col] = self.data[flat];
        }
        mat
    }
    /// Mode-n folding: reconstruct a tensor from a mode-n unfolded matrix.
    ///
    /// `mat` has shape `(shape[mode], prod(other dims))`.  `shape` is the
    /// target tensor shape.
    pub fn mode_n_fold(mat: &[Vec<f64>], mode: usize, shape: &[usize]) -> Self {
        let mut out = DenseTensor::zeros(shape);
        let rank = shape.len();
        let strides = out.strides();
        let total: usize = out.data.len();
        for flat in 0..total {
            let mut tmp = flat;
            let mut midx = vec![0usize; rank];
            for k in 0..rank {
                midx[k] = tmp / strides[k];
                tmp %= strides[k];
            }
            let row = midx[mode];
            let mut col = 0usize;
            let mut col_stride = 1usize;
            for k in (0..rank).rev() {
                if k == mode {
                    continue;
                }
                col += midx[k] * col_stride;
                col_stride *= shape[k];
            }
            out.data[flat] = mat[row][col];
        }
        out
    }
}
impl DenseTensor {
    /// Clone helper (needed since derive(Clone) isn't on the struct).
    #[allow(dead_code)]
    pub fn clone_tensor(&self) -> DenseTensor {
        DenseTensor {
            shape: self.shape.clone(),
            data: self.data.clone(),
        }
    }
    /// Scale all elements by a scalar.
    #[allow(dead_code)]
    pub fn scale(&self, s: f64) -> DenseTensor {
        DenseTensor {
            shape: self.shape.clone(),
            data: self.data.iter().map(|&x| x * s).collect(),
        }
    }
    /// Element-wise addition (shapes must match).
    #[allow(dead_code)]
    pub fn add_tensor(&self, other: &DenseTensor) -> DenseTensor {
        assert_eq!(
            self.shape, other.shape,
            "shapes must match for tensor addition"
        );
        let data: Vec<f64> = self
            .data
            .iter()
            .zip(other.data.iter())
            .map(|(a, b)| a + b)
            .collect();
        DenseTensor {
            shape: self.shape.clone(),
            data,
        }
    }
    /// Element-wise subtraction (shapes must match).
    #[allow(dead_code)]
    pub fn sub_tensor(&self, other: &DenseTensor) -> DenseTensor {
        assert_eq!(
            self.shape, other.shape,
            "shapes must match for tensor subtraction"
        );
        let data: Vec<f64> = self
            .data
            .iter()
            .zip(other.data.iter())
            .map(|(a, b)| a - b)
            .collect();
        DenseTensor {
            shape: self.shape.clone(),
            data,
        }
    }
    /// Max absolute element.
    #[allow(dead_code)]
    pub fn max_abs(&self) -> f64 {
        self.data
            .iter()
            .cloned()
            .fold(0.0f64, |acc, x| acc.max(x.abs()))
    }
    /// Sum all elements.
    #[allow(dead_code)]
    pub fn sum(&self) -> f64 {
        self.data.iter().sum()
    }
    /// Rank (number of modes).
    #[allow(dead_code)]
    pub fn rank(&self) -> usize {
        self.shape.len()
    }
    /// Total number of elements.
    #[allow(dead_code)]
    pub fn numel(&self) -> usize {
        self.data.len()
    }
}
/// Tensor basis transformation utilities.
pub struct TensorBasis;
impl TensorBasis {
    /// Build the 3×3 rotation matrix for an angle `theta` (radians) about the Z axis.
    pub fn rotation_z(theta: f64) -> Tensor2 {
        let c = theta.cos();
        let s = theta.sin();
        Tensor2::new([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    }
    /// Build the 3×3 rotation matrix for an angle `theta` about the X axis.
    pub fn rotation_x(theta: f64) -> Tensor2 {
        let c = theta.cos();
        let s = theta.sin();
        Tensor2::new([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    }
    /// Build the 3×3 rotation matrix for an angle `theta` about the Y axis.
    pub fn rotation_y(theta: f64) -> Tensor2 {
        let c = theta.cos();
        let s = theta.sin();
        Tensor2::new([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    }
    /// Transform the Voigt stiffness matrix (6×6) by a rotation matrix `r`.
    ///
    /// Uses the Bond transformation: M' = T * M * T^T.
    /// Returns the rotated 6×6 Voigt matrix.
    pub fn rotate_voigt_stiffness(m: &[[f64; 6]; 6], r: &Tensor2) -> [[f64; 6]; 6] {
        let c = Self::voigt_to_tensor4(m);
        let c_rot = c.rotate(r);
        c_rot.to_voigt_matrix()
    }
    /// Convert a Voigt stiffness matrix (6×6) to a Tensor4.
    fn voigt_to_tensor4(m: &[[f64; 6]; 6]) -> Tensor4 {
        let voigt_map: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
        let mut c = Tensor4::zero();
        for (p, &(i, j)) in voigt_map.iter().enumerate() {
            for (q, &(k, l)) in voigt_map.iter().enumerate() {
                let val = m[p][q];
                c.data[i][j][k][l] = val;
                c.data[j][i][k][l] = val;
                c.data[i][j][l][k] = val;
                c.data[j][i][l][k] = val;
            }
        }
        c
    }
}
/// A single Tensor Train core of shape (r_left, n, r_right).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TtCore {
    /// Left bond dimension r_{k-1}.
    pub r_left: usize,
    /// Mode dimension n_k.
    pub n: usize,
    /// Right bond dimension r_k.
    pub r_right: usize,
    /// Flat C-order data, length r_left * n * r_right.
    pub data: Vec<f64>,
}
#[allow(dead_code)]
impl TtCore {
    /// Create a zero TT core.
    pub fn zeros(r_left: usize, n: usize, r_right: usize) -> Self {
        TtCore {
            r_left,
            n,
            r_right,
            data: vec![0.0; r_left * n * r_right],
        }
    }
    /// Get element at (alpha, i, beta).
    pub fn get(&self, alpha: usize, i: usize, beta: usize) -> f64 {
        self.data[alpha * self.n * self.r_right + i * self.r_right + beta]
    }
    /// Set element at (alpha, i, beta).
    pub fn set(&mut self, alpha: usize, i: usize, beta: usize, val: f64) {
        let idx = alpha * self.n * self.r_right + i * self.r_right + beta;
        self.data[idx] = val;
    }
    /// Frobenius norm of this core.
    pub fn frobenius_norm(&self) -> f64 {
        self.data.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }
}
/// Tucker decomposition result for a 3-way tensor.
///
/// T ≈ G ×₁ U₀ ×₂ U₁ ×₃ U₂  where G is the core tensor (r0×r1×r2).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TuckerDecomposition {
    /// Core tensor of shape (r0, r1, r2).
    pub core: DenseTensor,
    /// Mode-0 factor matrix (n0×r0).
    pub u0: Vec<Vec<f64>>,
    /// Mode-1 factor matrix (n1×r1).
    pub u1: Vec<Vec<f64>>,
    /// Mode-2 factor matrix (n2×r2).
    pub u2: Vec<Vec<f64>>,
}
impl TuckerDecomposition {
    /// Compute reconstruction Frobenius norm.
    #[allow(dead_code)]
    pub fn reconstruction_norm(&self) -> f64 {
        tucker_reconstruct(self).frobenius_norm()
    }
    /// Relative error of this Tucker approximation.
    #[allow(dead_code)]
    pub fn relative_error(&self, original: &DenseTensor) -> f64 {
        let recon = tucker_reconstruct(self);
        let diff = original.sub_tensor(&recon).frobenius_norm();
        let norm = original.frobenius_norm();
        if norm < 1e-30 {
            return diff;
        }
        diff / norm
    }
}
/// Voigt notation utilities for symmetric 3×3 tensors.
///
/// Convention: \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\] → indices 0..5.
pub struct VoigtTensor;
impl VoigtTensor {
    /// Convert a `Tensor2` to its 6-component Voigt representation.
    /// Components: \[0\]=xx, \[1\]=yy, \[2\]=zz, \[3\]=xy, \[4\]=yz, \[5\]=xz.
    pub fn from_tensor2(t: &Tensor2) -> [f64; 6] {
        [
            t.data[0][0],
            t.data[1][1],
            t.data[2][2],
            t.data[0][1],
            t.data[1][2],
            t.data[0][2],
        ]
    }
    /// Reconstruct a symmetric `Tensor2` from a 6-component Voigt vector.
    pub fn to_tensor2(v: &[f64; 6]) -> Tensor2 {
        Tensor2 {
            data: [[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]],
        }
    }
}
impl VoigtTensor {
    /// Convert strain in Voigt notation to engineering strain.
    ///
    /// Engineering strain has shear components multiplied by 2:
    /// \[eps_xx, eps_yy, eps_zz, 2*eps_xy, 2*eps_yz, 2*eps_xz\].
    pub fn to_engineering_strain(v: &[f64; 6]) -> [f64; 6] {
        [v[0], v[1], v[2], 2.0 * v[3], 2.0 * v[4], 2.0 * v[5]]
    }
    /// Convert engineering strain back to tensor strain Voigt.
    pub fn from_engineering_strain(e: &[f64; 6]) -> [f64; 6] {
        [e[0], e[1], e[2], e[3] / 2.0, e[4] / 2.0, e[5] / 2.0]
    }
    /// Compute the von Mises stress from Voigt notation.
    pub fn von_mises_voigt(v: &[f64; 6]) -> f64 {
        let s = VoigtTensor::to_tensor2(v);
        s.von_mises()
    }
    /// Hydrostatic pressure from Voigt notation: p = -(v\[0\]+v\[1\]+v\[2\])/3.
    pub fn hydrostatic_pressure(v: &[f64; 6]) -> f64 {
        -(v[0] + v[1] + v[2]) / 3.0
    }
}
