//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// Dense n×n geometric stiffness matrix.
pub struct GeometricStiffness {
    pub(super) data: Vec<f64>,
    pub(super) n: usize,
}
impl GeometricStiffness {
    /// Create a new n×n zero geometric stiffness matrix.
    pub fn new(n: usize) -> Self {
        Self {
            data: vec![0.0; n * n],
            n,
        }
    }
    /// Set entry (i, j).
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }
    /// Get entry (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }
}
/// Euler column buckling: P_cr = π² E I / (K L)²
pub struct EulerBuckling {
    /// Young's modulus (Pa).
    pub e: f64,
    /// Second moment of area (m⁴).
    pub i: f64,
    /// Column length (m).
    pub l: f64,
    /// Effective length factor K (1.0 = pinned-pinned, 0.5 = fixed-fixed, 2.0 = fixed-free).
    pub k: f64,
}
impl EulerBuckling {
    /// Critical load P_cr = π² E I / (K L)²
    pub fn critical_load(&self) -> f64 {
        PI * PI * self.e * self.i / (self.k * self.l).powi(2)
    }
    /// Slenderness ratio λ = L_eff / r where L_eff = K * L.
    ///
    /// `_a` – cross-sectional area (unused here, kept for API completeness)
    /// `r`  – radius of gyration
    pub fn slenderness_ratio(&self, _a: f64, r: f64) -> f64 {
        (self.k * self.l) / r
    }
}
/// Buckling eigenvalue problem: find lambda such that (K - lambda*Kg)*phi = 0.
pub struct BucklingProblem {
    pub(super) k: StiffnessMatrix,
    pub(super) kg: GeometricStiffness,
    pub(super) n: usize,
}
impl BucklingProblem {
    /// Create a new buckling problem of size n×n.
    pub fn new(n: usize) -> Self {
        Self {
            k: StiffnessMatrix::new(n),
            kg: GeometricStiffness::new(n),
            n,
        }
    }
    /// Set elastic stiffness entry (i, j).
    pub fn set_k(&mut self, i: usize, j: usize, v: f64) {
        self.k.set(i, j, v);
    }
    /// Set geometric stiffness entry (i, j).
    pub fn set_kg(&mut self, i: usize, j: usize, v: f64) {
        self.kg.set(i, j, v);
    }
    /// Solve for the smallest `num_modes` buckling load factors via inverse power iteration.
    ///
    /// Returns load factors sorted in ascending order.
    pub fn solve_buckling_load_factors(&self, num_modes: usize) -> Vec<f64> {
        let n = self.n;
        let num_modes = num_modes.min(n);
        let mut factors = Vec::with_capacity(num_modes);
        for mode_idx in 0..num_modes {
            let shift = if mode_idx == 0 {
                0.0
            } else {
                factors[mode_idx - 1] * 1.01 + 1e-6
            };
            let lf = self.inverse_power_iteration(shift, 500, 1e-10);
            factors.push(lf);
        }
        factors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        factors
    }
    /// Perform inverse power iteration with shift `sigma` to find eigenvalue nearest sigma.
    fn inverse_power_iteration(&self, sigma: f64, max_iter: usize, tol: f64) -> f64 {
        let n = self.n;
        let mut a = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                a[i * n + j] = self.k.get(i, j) - sigma * self.kg.get(i, j);
            }
        }
        let mut v = vec![1.0f64 / (n as f64).sqrt(); n];
        let mut lambda = sigma + 1.0;
        for _ in 0..max_iter {
            let mut w = vec![0.0f64; n];
            for (i, w_i) in w.iter_mut().enumerate().take(n) {
                for (j, v_j) in v.iter().enumerate().take(n) {
                    *w_i += self.kg.get(i, j) * v_j;
                }
            }
            let z = match gaussian_solve(&a, &w, n) {
                Some(x) => x,
                None => return sigma + 1.0,
            };
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            for i in 0..n {
                let mut kv_i = 0.0f64;
                let mut kgv_i = 0.0f64;
                for (j, z_j) in z.iter().enumerate().take(n) {
                    kv_i += self.k.get(i, j) * z_j;
                    kgv_i += self.kg.get(i, j) * z_j;
                }
                num += z[i] * kv_i;
                den += z[i] * kgv_i;
            }
            let new_lambda = if den.abs() > 1e-300 {
                num / den
            } else {
                sigma + 1.0
            };
            let norm: f64 = z.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-300 {
                v = z.iter().map(|x| x / norm).collect();
            }
            if (new_lambda - lambda).abs() < tol * (1.0 + lambda.abs()) {
                return new_lambda;
            }
            lambda = new_lambda;
        }
        lambda
    }
}
/// Koiter's asymptotic post-buckling expansion for initial post-buckling.
///
/// In the vicinity of the bifurcation point, the load-displacement relation is:
///
///   λ = λ_c + a · ξ + b · ξ²
///
/// where ξ is the buckling mode amplitude, and a, b are Koiter coefficients.
///
/// For symmetric structures a = 0, and the behavior is governed by b.
pub struct KoiterExpansion {
    /// Critical load factor λ_c.
    pub lambda_c: f64,
    /// Koiter asymmetry coefficient a (0 for symmetric structures).
    pub a: f64,
    /// Koiter curvature coefficient b.
    pub b: f64,
}
impl KoiterExpansion {
    /// Create a new Koiter expansion.
    pub fn new(lambda_c: f64, a: f64, b: f64) -> Self {
        Self { lambda_c, a, b }
    }
    /// Load factor at mode amplitude ξ.
    ///
    /// λ(ξ) = λ_c + a·ξ + b·ξ²
    pub fn load_factor(&self, xi: f64) -> f64 {
        self.lambda_c + self.a * xi + self.b * xi * xi
    }
    /// Sensitivity to imperfections (stable vs unstable post-buckling).
    ///
    /// Returns `true` if the post-buckling path is stable (b > 0 for symmetric).
    pub fn is_stable(&self) -> bool {
        if self.a.abs() > 1e-14 {
            false
        } else {
            self.b > 0.0
        }
    }
    /// Maximum load factor including imperfection (Koiter's formula).
    ///
    /// For symmetric structures (a = 0):
    ///
    ///   λ_max = λ_c · (1 − c · √|δ|)    (unstable, b < 0)
    ///   λ_max = λ_c                        (stable,   b > 0)
    ///
    /// For asymmetric structures (a ≠ 0):
    ///
    ///   λ_max ≈ λ_c − (a · δ)^{2/3} / (6b)^{1/3}  (approximate)
    pub fn max_load_with_imperfection(&self, delta: f64, c: f64) -> f64 {
        if self.a.abs() < 1e-14 {
            if self.b < 0.0 {
                (self.lambda_c * (1.0 - c * delta.abs().sqrt())).max(0.0)
            } else {
                self.lambda_c
            }
        } else {
            let sign = if self.a * delta > 0.0 { 1.0 } else { -1.0 };
            let b_abs = self.b.abs().max(1e-30);
            let denom = (6.0 * b_abs).powf(1.0 / 3.0);
            (self.lambda_c - sign * (self.a.abs() * delta.abs()).powf(2.0 / 3.0) / denom).max(0.0)
        }
    }
    /// Post-buckling equilibrium amplitude ξ* at load factor λ > λ_c.
    ///
    /// Solves λ = λ_c + b·ξ² for symmetric case (a = 0).
    /// Returns None if b ≤ 0 (unstable) or λ < λ_c.
    pub fn equilibrium_amplitude(&self, lambda: f64) -> Option<f64> {
        if self.a.abs() < 1e-14 {
            if self.b <= 0.0 || lambda < self.lambda_c {
                return None;
            }
            Some(((lambda - self.lambda_c) / self.b).sqrt())
        } else {
            let delta_lambda = lambda - self.lambda_c;
            let discriminant = self.a * self.a + 4.0 * self.b * delta_lambda;
            if discriminant < 0.0 {
                return None;
            }
            Some((-self.a + discriminant.sqrt()) / (2.0 * self.b))
        }
    }
}
/// Dense n×n elastic stiffness matrix (internal helper).
pub(super) struct StiffnessMatrix {
    pub(super) data: Vec<f64>,
    pub(super) n: usize,
}
impl StiffnessMatrix {
    fn new(n: usize) -> Self {
        Self {
            data: vec![0.0; n * n],
            n,
        }
    }
    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }
    fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }
}
/// Sparse-triplet representation of the global buckling eigenvalue problem.
///
/// Stores (row, col, value) triplets for both the elastic and geometric
/// stiffness matrices.
pub struct BucklingAnalysis {
    /// Number of global degrees of freedom.
    pub n_dof: usize,
    /// Elastic stiffness triplets (row, col, value).
    pub k_global: Vec<(usize, usize, f64)>,
    /// Geometric stiffness triplets (row, col, value).
    pub kg_global: Vec<(usize, usize, f64)>,
}
impl BucklingAnalysis {
    /// Create a new analysis with `n_dof` degrees of freedom.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            k_global: Vec::new(),
            kg_global: Vec::new(),
        }
    }
    /// Assemble the 4×4 consistent geometric stiffness for the beam element
    /// defined by `element_nodes` (2 node indices) carrying axial load `N`
    /// and having length `L`.
    ///
    /// Each node contributes 2 DOFs: \[2*node, 2*node+1\].
    pub fn add_geometric_stiffness(&mut self, element_nodes: [usize; 2], n: f64, l: f64) {
        let kge = geometric_stiffness_beam(n, l);
        let dofs = [
            2 * element_nodes[0],
            2 * element_nodes[0] + 1,
            2 * element_nodes[1],
            2 * element_nodes[1] + 1,
        ];
        for a in 0..4 {
            for b in 0..4 {
                self.kg_global.push((dofs[a], dofs[b], kge[a][b]));
            }
        }
    }
    /// Rayleigh quotient estimate of the critical load factor.
    ///
    /// Assembles dense matrices from triplets and returns
    /// (v^T K v) / (v^T Kg v) for the uniform vector v = \[1, 1, …, 1\].
    pub fn critical_load_factor_estimate(&self) -> f64 {
        let n = self.n_dof;
        let mut k_dense = vec![0.0f64; n * n];
        let mut kg_dense = vec![0.0f64; n * n];
        for &(r, c, v) in &self.k_global {
            if r < n && c < n {
                k_dense[r * n + c] += v;
            }
        }
        for &(r, c, v) in &self.kg_global {
            if r < n && c < n {
                kg_dense[r * n + c] += v;
            }
        }
        let v_vec = vec![1.0f64; n];
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for i in 0..n {
            let mut kv_i = 0.0f64;
            let mut kgv_i = 0.0f64;
            for j in 0..n {
                kv_i += k_dense[i * n + j] * v_vec[j];
                kgv_i += kg_dense[i * n + j] * v_vec[j];
            }
            num += v_vec[i] * kv_i;
            den += v_vec[i] * kgv_i;
        }
        if den.abs() > 1e-300 {
            num / den
        } else {
            f64::INFINITY
        }
    }
}
/// Buckling mode: load factor and associated mode shape.
pub struct BucklingMode {
    /// The buckling load factor (eigenvalue).
    pub load_factor: f64,
    /// The normalized mode shape (eigenvector).
    pub mode_shape: Vec<f64>,
}
/// Slenderness classification:
/// * Short:  λ < λ_c → failure by yielding
/// * Medium: λ_c ≤ λ < λ_E → inelastic buckling (Johnson range)
/// * Long:   λ ≥ λ_E → elastic Euler buckling
pub enum SlendernessClass {
    /// Short column, failure by yielding.
    Short,
    /// Medium column, inelastic (Johnson) buckling.
    Medium,
    /// Long (slender) column, elastic Euler buckling.
    Long,
}
/// Buckling mode from an eigenvalue analysis: mode number, eigenvalue, and shape.
pub struct BucklingModeResult {
    /// Mode number (1-based).
    pub mode_num: usize,
    /// Eigenvalue (buckling load factor).
    pub eigenvalue: f64,
    /// Normalised mode shape vector.
    pub shape: Vec<f64>,
}
