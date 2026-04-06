//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::f64::consts::PI;

/// Dense n×n geometric stiffness matrix.
#[allow(dead_code)]
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
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct EulerBuckling {
    /// Young's modulus (Pa).
    pub E: f64,
    /// Second moment of area (m⁴).
    pub I: f64,
    /// Column length (m).
    pub L: f64,
    /// Effective length factor K (1.0 = pinned-pinned, 0.5 = fixed-fixed, 2.0 = fixed-free).
    pub K: f64,
}
impl EulerBuckling {
    /// Critical load P_cr = π² E I / (K L)²
    #[allow(dead_code)]
    pub fn critical_load(&self) -> f64 {
        PI * PI * self.E * self.I / (self.K * self.L).powi(2)
    }
    /// Slenderness ratio λ = L_eff / r where L_eff = K * L.
    ///
    /// `_a` – cross-sectional area (unused here, kept for API completeness)
    /// `r`  – radius of gyration
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn slenderness_ratio(&self, _a: f64, r: f64) -> f64 {
        (self.K * self.L) / r
    }
}
/// Buckling eigenvalue problem: find lambda such that (K - lambda*Kg)*phi = 0.
#[allow(dead_code)]
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
            for i in 0..n {
                for j in 0..n {
                    w[i] += self.kg.get(i, j) * v[j];
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
                for j in 0..n {
                    kv_i += self.k.get(i, j) * z[j];
                    kgv_i += self.kg.get(i, j) * z[j];
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(lambda_c: f64, a: f64, b: f64) -> Self {
        Self { lambda_c, a, b }
    }
    /// Load factor at mode amplitude ξ.
    ///
    /// λ(ξ) = λ_c + a·ξ + b·ξ²
    #[allow(dead_code)]
    pub fn load_factor(&self, xi: f64) -> f64 {
        self.lambda_c + self.a * xi + self.b * xi * xi
    }
    /// Sensitivity to imperfections (stable vs unstable post-buckling).
    ///
    /// Returns `true` if the post-buckling path is stable (b > 0 for symmetric).
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct BucklingAnalysis {
    /// Number of global degrees of freedom.
    pub n_dof: usize,
    /// Elastic stiffness triplets (row, col, value).
    pub K_global: Vec<(usize, usize, f64)>,
    /// Geometric stiffness triplets (row, col, value).
    pub Kg_global: Vec<(usize, usize, f64)>,
}
impl BucklingAnalysis {
    /// Create a new analysis with `n_dof` degrees of freedom.
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            K_global: Vec::new(),
            Kg_global: Vec::new(),
        }
    }
    /// Assemble the 4×4 consistent geometric stiffness for the beam element
    /// defined by `element_nodes` (2 node indices) carrying axial load `N`
    /// and having length `L`.
    ///
    /// Each node contributes 2 DOFs: \[2*node, 2*node+1\].
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn add_geometric_stiffness(&mut self, element_nodes: [usize; 2], N: f64, L: f64) {
        let kge = geometric_stiffness_beam(N, L);
        let dofs = [
            2 * element_nodes[0],
            2 * element_nodes[0] + 1,
            2 * element_nodes[1],
            2 * element_nodes[1] + 1,
        ];
        for a in 0..4 {
            for b in 0..4 {
                self.Kg_global.push((dofs[a], dofs[b], kge[a][b]));
            }
        }
    }
    /// Rayleigh quotient estimate of the critical load factor.
    ///
    /// Assembles dense matrices from triplets and returns
    /// (v^T K v) / (v^T Kg v) for the uniform vector v = \[1, 1, …, 1\].
    #[allow(dead_code)]
    pub fn critical_load_factor_estimate(&self) -> f64 {
        let n = self.n_dof;
        let mut k_dense = vec![0.0f64; n * n];
        let mut kg_dense = vec![0.0f64; n * n];
        for &(r, c, v) in &self.K_global {
            if r < n && c < n {
                k_dense[r * n + c] += v;
            }
        }
        for &(r, c, v) in &self.Kg_global {
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
#[allow(dead_code)]
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
#[allow(dead_code)]
pub enum SlendernessClass {
    /// Short column, failure by yielding.
    Short,
    /// Medium column, inelastic (Johnson) buckling.
    Medium,
    /// Long (slender) column, elastic Euler buckling.
    Long,
}
/// Buckling mode from an eigenvalue analysis: mode number, eigenvalue, and shape.
#[allow(dead_code)]
pub struct BucklingModeResult {
    /// Mode number (1-based).
    pub mode_num: usize,
    /// Eigenvalue (buckling load factor).
    pub eigenvalue: f64,
    /// Normalised mode shape vector.
    pub shape: Vec<f64>,
}
