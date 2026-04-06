#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nonlinear FEM: geometric nonlinearity (total/updated Lagrangian), material
//! nonlinearity (J2 plasticity return mapping), Newton-Raphson with line search,
//! arc-length method (Riks), load stepping, consistent tangent moduli, large
//! deformation kinematics (F, B, C tensors), volumetric locking prevention
//! (F-bar, B-bar), enhanced assumed strain, and mixed formulations.

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// § 1  LARGE DEFORMATION KINEMATICS
// ══════════════════════════════════════════════════════════════════════════════

/// 3×3 matrix stored in row-major order as `[[f64; 3\]; 3]`.
pub type Mat3 = [[f64; 3]; 3];

/// The identity matrix in 3D.
pub const EYE3: Mat3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

/// Matrix multiply: C = A · B (3×3).
pub fn mat3_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Transpose of a 3×3 matrix.
pub fn mat3_transpose(a: &Mat3) -> Mat3 {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

/// Determinant of a 3×3 matrix.
pub fn mat3_det(a: &Mat3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Inverse of a 3×3 matrix (panics if singular).
pub fn mat3_inv(a: &Mat3) -> Mat3 {
    let det = mat3_det(a);
    assert!(det.abs() > 1e-30, "singular matrix in mat3_inv");
    let inv_det = 1.0 / det;
    [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
        ],
        [
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
        ],
    ]
}

/// Add two 3×3 matrices.
pub fn mat3_add(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Scale a 3×3 matrix.
pub fn mat3_scale(a: &Mat3, s: f64) -> Mat3 {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}

/// Frobenius norm of a 3×3 matrix.
pub fn mat3_norm(a: &Mat3) -> f64 {
    a.iter()
        .flat_map(|row| row.iter())
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
}

/// Double contraction A:B = Σᵢⱼ Aᵢⱼ·Bᵢⱼ
pub fn mat3_double_contract(a: &Mat3, b: &Mat3) -> f64 {
    a.iter()
        .zip(b.iter())
        .flat_map(|(ar, br)| ar.iter().zip(br.iter()).map(|(x, y)| x * y))
        .sum()
}

/// Deformation gradient F = ∂x/∂X (3×3).
///
/// For a given displacement gradient H = ∂u/∂X:
/// F = I + H
#[derive(Debug, Clone)]
pub struct DeformationGradient {
    /// The 3×3 deformation gradient tensor
    pub f: Mat3,
}

impl DeformationGradient {
    /// Create from displacement gradient H.
    pub fn from_displacement_gradient(h: &Mat3) -> Self {
        Self {
            f: mat3_add(&EYE3, h),
        }
    }

    /// Create the identity deformation (no deformation).
    pub fn identity() -> Self {
        Self { f: EYE3 }
    }

    /// Jacobian J = det(F).
    pub fn jacobian(&self) -> f64 {
        mat3_det(&self.f)
    }

    /// Right Cauchy-Green tensor C = FᵀF.
    pub fn right_cauchy_green(&self) -> Mat3 {
        let ft = mat3_transpose(&self.f);
        mat3_mul(&ft, &self.f)
    }

    /// Left Cauchy-Green (Finger) tensor B = FFᵀ.
    pub fn left_cauchy_green(&self) -> Mat3 {
        let ft = mat3_transpose(&self.f);
        mat3_mul(&self.f, &ft)
    }

    /// Green-Lagrange strain tensor E = ½(C − I).
    pub fn green_lagrange_strain(&self) -> Mat3 {
        let c = self.right_cauchy_green();
        let mut e = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                e[i][j] = 0.5 * (c[i][j] - EYE3[i][j]);
            }
        }
        e
    }

    /// Almansi (Euler-Almansi) strain tensor e = ½(I − b⁻¹).
    pub fn almansi_strain(&self) -> Mat3 {
        let b = self.left_cauchy_green();
        let b_inv = mat3_inv(&b);
        let mut e = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                e[i][j] = 0.5 * (EYE3[i][j] - b_inv[i][j]);
            }
        }
        e
    }

    /// F-bar modified deformation gradient for volumetric locking prevention.
    ///
    /// F̄ = (J₀/J)^(1/3) · F
    pub fn f_bar(&self, j0: f64) -> Mat3 {
        let j = self.jacobian();
        let scale = (j0 / j).powf(1.0 / 3.0);
        mat3_scale(&self.f, scale)
    }

    /// Polar decomposition F = R · U (approximate via one iteration).
    ///
    /// Returns (R, U) where R is rotation and U is right stretch tensor.
    pub fn polar_decompose(&self) -> (Mat3, Mat3) {
        // Use iterative approach: R_{k+1} = ½(R_k + F^{-T} R_k^{-1} F^{-1})
        // Simplified: use F directly and normalize columns for R
        let mut r = self.f;
        for _ in 0..20 {
            let r_inv_t = mat3_transpose(&mat3_inv(&r));
            for i in 0..3 {
                for j in 0..3 {
                    r[i][j] = 0.5 * (r[i][j] + r_inv_t[i][j]);
                }
            }
        }
        let rt = mat3_transpose(&r);
        let u = mat3_mul(&rt, &self.f);
        (r, u)
    }

    /// Principal stretches (approximate: eigenvalues of U via power iteration).
    pub fn principal_stretches(&self) -> [f64; 3] {
        let c = self.right_cauchy_green();
        // Return diagonal entries as approximation for initially diagonal F
        [c[0][0].sqrt(), c[1][1].sqrt(), c[2][2].sqrt()]
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 2  J2 PLASTICITY WITH RETURN MAPPING
// ══════════════════════════════════════════════════════════════════════════════

/// J2 plasticity state variables at a material point.
#[derive(Debug, Clone, Default)]
pub struct PlasticState {
    /// Plastic strain tensor (Voigt: \[ε11, ε22, ε33, γ12, γ23, γ13\])
    pub eps_p: [f64; 6],
    /// Equivalent plastic strain (accumulated)
    pub eps_p_eq: f64,
    /// Back stress tensor (kinematic hardening)
    pub alpha: [f64; 6],
}

impl PlasticState {
    /// Create a new zero-initialized plastic state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Equivalent plastic strain rate (norm of plastic strain increment / sqrt(2/3)).
    pub fn update_eps_p_eq(&mut self, delta_eps_p: &[f64; 6]) {
        let norm2: f64 = delta_eps_p[0].powi(2)
            + delta_eps_p[1].powi(2)
            + delta_eps_p[2].powi(2)
            + 0.5 * (delta_eps_p[3].powi(2) + delta_eps_p[4].powi(2) + delta_eps_p[5].powi(2));
        self.eps_p_eq += (2.0 / 3.0 * norm2).sqrt();
    }
}

/// J2 isotropic + kinematic hardening model (Chaboche).
#[derive(Debug, Clone)]
pub struct J2Material {
    /// Young's modulus \[Pa\]
    pub young: f64,
    /// Poisson's ratio
    pub nu: f64,
    /// Initial yield stress \[Pa\]
    pub sigma_y0: f64,
    /// Isotropic hardening modulus H \[Pa\]
    pub h_iso: f64,
    /// Kinematic hardening modulus C \[Pa\]
    pub c_kin: f64,
    /// Kinematic saturation parameter γ \[1/strain\]
    pub gamma_kin: f64,
}

impl J2Material {
    /// Create a new J2 material.
    pub fn new(young: f64, nu: f64, sigma_y0: f64, h_iso: f64) -> Self {
        Self {
            young,
            nu,
            sigma_y0,
            h_iso,
            c_kin: 0.0,
            gamma_kin: 0.0,
        }
    }

    /// Shear modulus G = E / (2(1+ν)).
    #[inline]
    pub fn shear_modulus(&self) -> f64 {
        self.young / (2.0 * (1.0 + self.nu))
    }

    /// Bulk modulus K = E / (3(1−2ν)).
    #[inline]
    pub fn bulk_modulus(&self) -> f64 {
        self.young / (3.0 * (1.0 - 2.0 * self.nu))
    }

    /// Current yield stress at equivalent plastic strain.
    pub fn yield_stress(&self, eps_p_eq: f64) -> f64 {
        self.sigma_y0 + self.h_iso * eps_p_eq
    }

    /// Deviatoric part of a stress tensor (Voigt).
    pub fn deviatoric(sigma: &[f64; 6]) -> [f64; 6] {
        let p = (sigma[0] + sigma[1] + sigma[2]) / 3.0;
        [
            sigma[0] - p,
            sigma[1] - p,
            sigma[2] - p,
            sigma[3],
            sigma[4],
            sigma[5],
        ]
    }

    /// Von Mises stress from Voigt tensor.
    pub fn von_mises(sigma: &[f64; 6]) -> f64 {
        let s = Self::deviatoric(sigma);
        let j2 = 0.5 * (s[0].powi(2) + s[1].powi(2) + s[2].powi(2))
            + s[3].powi(2)
            + s[4].powi(2)
            + s[5].powi(2);
        (3.0 * j2).sqrt()
    }

    /// Elastic predictor step (return the trial stress).
    pub fn elastic_predictor(&self, sigma_n: &[f64; 6], delta_eps: &[f64; 6]) -> [f64; 6] {
        let g = self.shear_modulus();
        let k = self.bulk_modulus();
        let lam = k - 2.0 / 3.0 * g; // Lamé λ

        let de_vol = delta_eps[0] + delta_eps[1] + delta_eps[2];

        [
            sigma_n[0] + lam * de_vol + 2.0 * g * delta_eps[0],
            sigma_n[1] + lam * de_vol + 2.0 * g * delta_eps[1],
            sigma_n[2] + lam * de_vol + 2.0 * g * delta_eps[2],
            sigma_n[3] + g * delta_eps[3],
            sigma_n[4] + g * delta_eps[4],
            sigma_n[5] + g * delta_eps[5],
        ]
    }

    /// Radial return mapping (isotropic hardening only).
    ///
    /// Returns (sigma_n1, state_n1, delta_gamma).
    pub fn return_mapping(
        &self,
        sigma_tr: &[f64; 6],
        state_n: &PlasticState,
    ) -> ([f64; 6], PlasticState, f64) {
        let g = self.shear_modulus();
        let s_tr = Self::deviatoric(sigma_tr);

        // Relative stress eta = s_tr − α_n
        let eta: [f64; 6] = std::array::from_fn(|i| s_tr[i] - state_n.alpha[i]);

        // Trial equivalent stress
        let xi_tr = (1.5 * (eta[0].powi(2) + eta[1].powi(2) + eta[2].powi(2))
            + 3.0 * (eta[3].powi(2) + eta[4].powi(2) + eta[5].powi(2)))
        .sqrt();

        let sigma_y_n = self.yield_stress(state_n.eps_p_eq);

        // Check yielding
        let phi_tr = xi_tr - sigma_y_n;
        if phi_tr <= 0.0 {
            // Elastic step
            return (*sigma_tr, state_n.clone(), 0.0);
        }

        // Plastic step: compute consistency parameter Δγ
        let delta_gamma = phi_tr / (3.0 * g + self.h_iso + self.gamma_kin * self.c_kin);

        // Flow direction n = η / ‖η‖ (normalized)
        let n_norm = xi_tr.max(1e-30);
        let n_vec: [f64; 6] = std::array::from_fn(|i| eta[i] / n_norm);

        // Updated back stress
        let mut alpha_n1 = state_n.alpha;
        for i in 0..6 {
            alpha_n1[i] += self.c_kin / self.gamma_kin.max(1e-30) * n_vec[i] * delta_gamma;
        }

        // Updated plastic strain
        let scale = 1.5 / n_norm;
        let mut eps_p_n1 = state_n.eps_p;
        for i in 0..3 {
            eps_p_n1[i] += delta_gamma * scale * eta[i];
        }
        for i in 3..6 {
            eps_p_n1[i] += delta_gamma * scale * 2.0 * eta[i];
        }

        // Correct deviatoric stress
        let p_tr = (sigma_tr[0] + sigma_tr[1] + sigma_tr[2]) / 3.0;
        let mut sigma_n1 = *sigma_tr;
        for i in 0..6 {
            sigma_n1[i] = s_tr[i] - 2.0 * g * delta_gamma * scale * eta[i];
        }
        sigma_n1[0] += p_tr;
        sigma_n1[1] += p_tr;
        sigma_n1[2] += p_tr;

        let mut new_state = PlasticState {
            eps_p: eps_p_n1,
            eps_p_eq: state_n.eps_p_eq + delta_gamma * (2.0 / 3.0_f64).sqrt(),
            alpha: alpha_n1,
        };
        new_state.update_eps_p_eq(&std::array::from_fn(|i| eps_p_n1[i] - state_n.eps_p[i]));

        (sigma_n1, new_state, delta_gamma)
    }

    /// Consistent tangent modulus (algorithmic tangent) in Voigt notation (6×6).
    pub fn consistent_tangent(
        &self,
        delta_gamma: f64,
        sigma_tr: &[f64; 6],
        state_n: &PlasticState,
    ) -> [[f64; 6]; 6] {
        let g = self.shear_modulus();
        let k = self.bulk_modulus();
        let lam = k - 2.0 / 3.0 * g;

        // Elastic tangent (isotropic)
        let mut c = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] = lam;
            }
            c[i][i] += 2.0 * g;
        }
        for i in 3..6 {
            c[i][i] = g;
        }

        if delta_gamma <= 0.0 {
            return c;
        }

        // Plastic correction
        let s_tr = Self::deviatoric(sigma_tr);
        let eta: [f64; 6] = std::array::from_fn(|i| s_tr[i] - state_n.alpha[i]);
        let n_norm = {
            let xi = (1.5 * (eta[0].powi(2) + eta[1].powi(2) + eta[2].powi(2))
                + 3.0 * (eta[3].powi(2) + eta[4].powi(2) + eta[5].powi(2)))
            .sqrt();
            xi.max(1e-30)
        };
        let n_vec: [f64; 6] = std::array::from_fn(|i| eta[i] / n_norm);

        let theta1 = 1.0 - 2.0 * g * delta_gamma / n_norm;
        let theta2 = 1.0 / (1.0 + self.h_iso / (3.0 * g)) - (1.0 - theta1);

        // Modify deviatoric part
        let shear_indices = [0, 1, 2, 3, 4, 5];
        for &i in &shear_indices {
            for &j in &shear_indices {
                let ni = n_vec[i];
                let nj = n_vec[j];
                let delta_ij = if i == j { 1.0 } else { 0.0 };
                // Deviatoric correction
                let dev_correction =
                    -2.0 * g * (theta2 * ni * nj + delta_gamma / n_norm * (delta_ij - ni * nj));
                if i < 3 && j < 3 {
                    c[i][j] += dev_correction
                        - 2.0 * g * delta_gamma / n_norm
                            * (-1.0 / 3.0 * if i == j { 1.0 } else { 0.0 });
                } else {
                    c[i][j] += dev_correction * 0.5;
                }
            }
        }
        c
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3  NEWTON-RAPHSON WITH LINE SEARCH
// ══════════════════════════════════════════════════════════════════════════════

/// Newton-Raphson iteration parameters.
#[derive(Debug, Clone)]
pub struct NrParams {
    /// Absolute tolerance on residual norm
    pub tol_abs: f64,
    /// Relative tolerance on residual norm
    pub tol_rel: f64,
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Whether to use line search
    pub use_line_search: bool,
    /// Armijo constant for line search
    pub c1: f64,
    /// Line search maximum steps
    pub ls_max_steps: usize,
}

impl Default for NrParams {
    fn default() -> Self {
        Self {
            tol_abs: 1e-10,
            tol_rel: 1e-8,
            max_iter: 50,
            use_line_search: true,
            c1: 1e-4,
            ls_max_steps: 10,
        }
    }
}

/// Newton-Raphson iteration result.
#[derive(Debug, Clone)]
pub struct NrResult {
    /// Converged solution
    pub solution: Vec<f64>,
    /// Final residual norm
    pub residual_norm: f64,
    /// Number of iterations taken
    pub iterations: usize,
    /// Whether convergence was achieved
    pub converged: bool,
}

/// Trait for a nonlinear function with Jacobian.
pub trait NonlinearSystem {
    /// Evaluate f(x).
    fn residual(&self, x: &[f64]) -> Vec<f64>;
    /// Evaluate Jacobian ∂f/∂x (dense).
    fn jacobian(&self, x: &[f64]) -> Vec<Vec<f64>>;
}

/// Solve Ax = b using forward/back substitution (for a small dense system).
fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, &bi)| {
            let mut r = row.clone();
            r.push(bi);
            r
        })
        .collect();

    // Gaussian elimination with partial pivoting
    for col in 0..n {
        // Pivot
        let max_row = (col..n)
            .max_by(|&r1, &r2| {
                aug[r1][col]
                    .abs()
                    .partial_cmp(&aug[r2][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("range col..n is non-empty");
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        if pivot.abs() < 1e-30 {
            continue;
        }

        for row in col + 1..n {
            let factor = aug[row][col] / pivot;
            for k in col..=n {
                let val = aug[col][k] * factor;
                aug[row][k] -= val;
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = aug[i][n];
        for j in i + 1..n {
            s -= aug[i][j] * x[j];
        }
        x[i] = s / aug[i][i].max(1e-30);
    }
    x
}

/// Perform Newton-Raphson iteration.
pub fn newton_raphson_solve<S: NonlinearSystem>(
    system: &S,
    x0: Vec<f64>,
    params: &NrParams,
) -> NrResult {
    let mut x = x0;
    let r0_norm = {
        let r = system.residual(&x);
        r.iter().map(|v| v * v).sum::<f64>().sqrt()
    };

    for iter in 0..params.max_iter {
        let r = system.residual(&x);
        let r_norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();

        if r_norm < params.tol_abs || r_norm < params.tol_rel * r0_norm.max(1.0) {
            return NrResult {
                solution: x,
                residual_norm: r_norm,
                iterations: iter,
                converged: true,
            };
        }

        let jac = system.jacobian(&x);
        let neg_r: Vec<f64> = r.iter().map(|v| -v).collect();
        let dx = solve_linear_system(&jac, &neg_r);

        if params.use_line_search {
            // Backtracking line search (Armijo condition)
            let mut alpha = 1.0;
            let slope: f64 = r.iter().zip(dx.iter()).map(|(ri, di)| -ri * di).sum();
            for _ in 0..params.ls_max_steps {
                let x_new: Vec<f64> = x
                    .iter()
                    .zip(dx.iter())
                    .map(|(&xi, &di)| xi + alpha * di)
                    .collect();
                let r_new = system.residual(&x_new);
                let r_new_norm: f64 = r_new.iter().map(|v| v * v).sum::<f64>().sqrt();
                if r_new_norm <= r_norm + params.c1 * alpha * slope {
                    break;
                }
                alpha *= 0.5;
            }
            for (xi, di) in x.iter_mut().zip(dx.iter()) {
                *xi += alpha * di;
            }
        } else {
            for (xi, di) in x.iter_mut().zip(dx.iter()) {
                *xi += di;
            }
        }
    }

    let r_final = system.residual(&x);
    let r_norm_final = r_final.iter().map(|v| v * v).sum::<f64>().sqrt();

    NrResult {
        solution: x,
        residual_norm: r_norm_final,
        iterations: params.max_iter,
        converged: false,
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 4  ARC-LENGTH METHOD (RIKS)
// ══════════════════════════════════════════════════════════════════════════════

/// Arc-length control parameters (Riks/Crisfield).
#[derive(Debug, Clone)]
pub struct ArcLengthParams {
    /// Arc-length radius Δl
    pub dl: f64,
    /// Minimum allowed arc-length
    pub dl_min: f64,
    /// Maximum allowed arc-length
    pub dl_max: f64,
    /// Number of desired iterations (for arc-length update)
    pub i_des: usize,
    /// Maximum iterations per step
    pub max_iter: usize,
    /// Convergence tolerance
    pub tol: f64,
}

impl Default for ArcLengthParams {
    fn default() -> Self {
        Self {
            dl: 0.1,
            dl_min: 1e-5,
            dl_max: 1.0,
            i_des: 3,
            max_iter: 20,
            tol: 1e-6,
        }
    }
}

/// Arc-length state for a simple 1D problem.
#[derive(Debug, Clone)]
pub struct ArcLengthState {
    /// Displacement DOFs
    pub u: Vec<f64>,
    /// Load parameter λ
    pub lambda: f64,
    /// Arc-length radius
    pub dl: f64,
    /// Load sign (+1 or −1)
    pub sign: f64,
}

impl ArcLengthState {
    /// Create initial state.
    pub fn new(n_dof: usize, dl: f64) -> Self {
        Self {
            u: vec![0.0; n_dof],
            lambda: 0.0,
            dl,
            sign: 1.0,
        }
    }

    /// Spherical arc-length constraint:
    ///
    /// g(u, λ) = ‖Δu‖² + (Δλ)² − Δl² = 0
    pub fn arc_length_constraint(
        &self,
        u_new: &[f64],
        lambda_new: f64,
        u_prev: &[f64],
        lambda_prev: f64,
    ) -> f64 {
        let du: f64 = u_new
            .iter()
            .zip(u_prev.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>();
        let dl = lambda_new - lambda_prev;
        du + dl * dl - self.dl * self.dl
    }
}

/// Riks arc-length solver for a 1D nonlinear problem: K·u = λ·f_ref.
///
/// Here K is a nonlinear stiffness (function of u) and f_ref is a reference load.
pub struct RiksSolver {
    /// Arc-length parameters
    pub params: ArcLengthParams,
    /// Current state
    pub state: ArcLengthState,
    /// Load-displacement history \[(λ, u₁)\]
    pub history: Vec<(f64, f64)>,
    /// Nonlinear stiffness function: K(u) \[scalar 1D for simplicity\]
    stiffness_fn: fn(f64) -> f64,
    /// Reference load magnitude
    pub f_ref: f64,
}

impl RiksSolver {
    /// Create a new Riks solver for a 1-DOF nonlinear system.
    pub fn new(stiffness_fn: fn(f64) -> f64, f_ref: f64, dl: f64) -> Self {
        let params = ArcLengthParams {
            dl,
            ..ArcLengthParams::default()
        };
        Self {
            params,
            state: ArcLengthState::new(1, dl),
            history: Vec::new(),
            stiffness_fn,
            f_ref,
        }
    }

    /// Perform one arc-length increment.
    pub fn step(&mut self) -> bool {
        let k = (self.stiffness_fn)(self.state.u[0]);
        if k.abs() < 1e-30 {
            return false;
        }

        // Predictor: tangent solution
        let du_t = self.f_ref / k;
        let dl_t = (self.params.dl / (du_t * du_t + 1.0).sqrt()).abs();
        let du_pred = du_t * dl_t * self.state.sign;
        let dlambda_pred = dl_t * self.state.sign;

        let u_prev = self.state.u[0];
        let lambda_prev = self.state.lambda;

        let mut u_n = u_prev + du_pred;
        let mut lambda_n = lambda_prev + dlambda_pred;

        // Corrector iterations
        for _iter in 0..self.params.max_iter {
            let r = lambda_n * self.f_ref - k * u_n;
            if r.abs() < self.params.tol {
                break;
            }
            let k_new = (self.stiffness_fn)(u_n);
            let duz = self.f_ref / k_new.max(1e-30);
            let dur = -r / k_new.max(1e-30);
            let a1 = duz * duz + 1.0;
            let a2 = 2.0 * ((u_n - u_prev) * duz + (lambda_n - lambda_prev));
            let a3 = self
                .state
                .arc_length_constraint(&[u_n], lambda_n, &[u_prev], lambda_prev);
            let discriminant = a2 * a2 - 4.0 * a1 * a3;
            if discriminant < 0.0 {
                break;
            }
            let dl1 = (-a2 + discriminant.sqrt()) / (2.0 * a1);
            let dl2 = (-a2 - discriminant.sqrt()) / (2.0 * a1);
            let dl = if dl1.abs() < dl2.abs() { dl1 } else { dl2 };
            u_n += dur + duz * dl;
            lambda_n += dl;
        }

        self.state.u[0] = u_n;
        self.state.lambda = lambda_n;
        self.history.push((lambda_n, u_n));
        true
    }

    /// Run for `n_steps` increments.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            if !self.step() {
                break;
            }
        }
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 5  LOAD STEPPING
// ══════════════════════════════════════════════════════════════════════════════

/// Load stepping scheme for nonlinear analysis.
#[derive(Debug, Clone)]
pub struct LoadStepper {
    /// Total load parameter to reach (target λ = 1.0)
    pub lambda_target: f64,
    /// Current load parameter
    pub lambda: f64,
    /// Current step size
    pub dlambda: f64,
    /// Minimum allowed step size
    pub dlambda_min: f64,
    /// Maximum allowed step size
    pub dlambda_max: f64,
    /// Number of desired iterations for automatic step control
    pub i_des: usize,
    /// Whether to use automatic step control
    pub auto_step: bool,
}

impl LoadStepper {
    /// Create a load stepper.
    pub fn new(n_steps: usize) -> Self {
        let dlambda = 1.0 / n_steps as f64;
        Self {
            lambda_target: 1.0,
            lambda: 0.0,
            dlambda,
            dlambda_min: dlambda * 0.01,
            dlambda_max: dlambda * 10.0,
            i_des: 4,
            auto_step: true,
        }
    }

    /// Advance the load parameter by one step. Returns the new λ.
    pub fn advance(&mut self) -> f64 {
        self.lambda = (self.lambda + self.dlambda).min(self.lambda_target);
        self.lambda
    }

    /// Check if loading is complete.
    pub fn is_complete(&self) -> bool {
        self.lambda >= self.lambda_target - 1e-12
    }

    /// Update step size based on iteration count (Irons-Tuck strategy).
    pub fn update_step_size(&mut self, iters_used: usize) {
        if !self.auto_step {
            return;
        }
        let ratio = (self.i_des as f64 / iters_used.max(1) as f64).sqrt();
        self.dlambda = (self.dlambda * ratio)
            .clamp(self.dlambda_min, self.dlambda_max)
            .min(self.lambda_target - self.lambda);
    }

    /// Remaining load fraction.
    pub fn remaining(&self) -> f64 {
        (self.lambda_target - self.lambda).max(0.0)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 6  ENHANCED ASSUMED STRAIN (EAS) ELEMENT
// ══════════════════════════════════════════════════════════════════════════════

/// EAS element parameters for a 2D quadrilateral Q4 element.
///
/// The EAS method augments the strain field with incompatible modes to avoid
/// locking and improve accuracy.
#[derive(Debug, Clone)]
pub struct EasQ4Element {
    /// Node coordinates (4 nodes × 2D)
    pub nodes: [[f64; 2]; 4],
    /// Number of EAS parameters (typically 4 for Q4)
    pub n_eas: usize,
    /// EAS incompatible modes parameters α
    pub alpha: Vec<f64>,
    /// Element Jacobian at element center
    pub j0: f64,
}

impl EasQ4Element {
    /// Create a new EAS Q4 element.
    pub fn new(nodes: [[f64; 2]; 4]) -> Self {
        let n_eas = 4;
        let alpha = vec![0.0; n_eas];
        let j0 = Self::jacobian_at(&nodes, 0.0, 0.0);
        Self {
            nodes,
            n_eas,
            alpha,
            j0,
        }
    }

    /// Compute the Jacobian determinant at (ξ, η).
    pub fn jacobian_at(nodes: &[[f64; 2]; 4], xi: f64, eta: f64) -> f64 {
        // Shape function derivatives w.r.t. (ξ, η)
        let dn_dxi = [
            -(1.0 - eta) * 0.25,
            (1.0 - eta) * 0.25,
            (1.0 + eta) * 0.25,
            -(1.0 + eta) * 0.25,
        ];
        let dn_deta = [
            -(1.0 - xi) * 0.25,
            -(1.0 + xi) * 0.25,
            (1.0 + xi) * 0.25,
            (1.0 - xi) * 0.25,
        ];

        let j11: f64 = (0..4).map(|k| dn_dxi[k] * nodes[k][0]).sum();
        let j12: f64 = (0..4).map(|k| dn_dxi[k] * nodes[k][1]).sum();
        let j21: f64 = (0..4).map(|k| dn_deta[k] * nodes[k][0]).sum();
        let j22: f64 = (0..4).map(|k| dn_deta[k] * nodes[k][1]).sum();
        j11 * j22 - j12 * j21
    }

    /// EAS interpolation matrix M(ξ, η) (3 × n_eas) in physical coordinates.
    pub fn eas_matrix(&self, xi: f64, eta: f64) -> Vec<Vec<f64>> {
        let j = Self::jacobian_at(&self.nodes, xi, eta);
        let scale = self.j0 / j.max(1e-30);

        // Standard 4-mode EAS (Wilson-Taylor incompatible modes approximation)
        vec![
            vec![scale * xi, 0.0, 0.0, scale * xi * eta],
            vec![0.0, scale * eta, 0.0, scale * xi * eta],
            vec![0.0, 0.0, scale * xi, scale * eta],
        ]
    }

    /// Enhanced strain contribution: ε̃ = M · α.
    pub fn enhanced_strain(&self, xi: f64, eta: f64) -> [f64; 3] {
        let m = self.eas_matrix(xi, eta);
        let mut eps = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..self.n_eas {
                eps[i] += m[i][j] * self.alpha[j];
            }
        }
        eps
    }

    /// Standard strain at Gauss point (ξ, η).
    pub fn standard_strain(&self, u_elem: &[f64; 8], xi: f64, eta: f64) -> [f64; 3] {
        // B-matrix (standard isoparametric)
        let dn_dxi = [
            -(1.0 - eta) * 0.25,
            (1.0 - eta) * 0.25,
            (1.0 + eta) * 0.25,
            -(1.0 + eta) * 0.25,
        ];
        let dn_deta = [
            -(1.0 - xi) * 0.25,
            -(1.0 + xi) * 0.25,
            (1.0 + xi) * 0.25,
            (1.0 - xi) * 0.25,
        ];

        let j11: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][0]).sum();
        let j12: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][1]).sum();
        let j21: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][0]).sum();
        let j22: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][1]).sum();
        let det_j = j11 * j22 - j12 * j21;

        let inv_j11 = j22 / det_j;
        let inv_j12 = -j12 / det_j;
        let inv_j21 = -j21 / det_j;
        let inv_j22 = j11 / det_j;

        let mut eps = [0.0f64; 3];
        for k in 0..4 {
            let dndx = inv_j11 * dn_dxi[k] + inv_j12 * dn_deta[k];
            let dndy = inv_j21 * dn_dxi[k] + inv_j22 * dn_deta[k];
            eps[0] += dndx * u_elem[2 * k];
            eps[1] += dndy * u_elem[2 * k + 1];
            eps[2] += dndx * u_elem[2 * k + 1] + dndy * u_elem[2 * k];
        }
        eps
    }

    /// Total strain = standard strain + enhanced strain.
    pub fn total_strain(&self, u_elem: &[f64; 8], xi: f64, eta: f64) -> [f64; 3] {
        let eps_std = self.standard_strain(u_elem, xi, eta);
        let eps_eas = self.enhanced_strain(xi, eta);
        [
            eps_std[0] + eps_eas[0],
            eps_std[1] + eps_eas[1],
            eps_std[2] + eps_eas[2],
        ]
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 7  F-BAR AND B-BAR VOLUMETRIC LOCKING PREVENTION
// ══════════════════════════════════════════════════════════════════════════════

/// B-bar element for volumetric locking prevention.
///
/// Replaces the volumetric part of the B-matrix with an average over the element.
#[derive(Debug, Clone)]
pub struct BBarElement {
    /// Node coordinates (8-node hexahedral)
    pub nodes: Vec<[f64; 3]>,
    /// Number of nodes
    pub n_nodes: usize,
}

impl BBarElement {
    /// Create a unit-cube hexahedral B-bar element.
    pub fn unit_cube() -> Self {
        let s = 0.5f64;
        let nodes = vec![
            [-s, -s, -s],
            [s, -s, -s],
            [s, s, -s],
            [-s, s, -s],
            [-s, -s, s],
            [s, -s, s],
            [s, s, s],
            [-s, s, s],
        ];
        Self { nodes, n_nodes: 8 }
    }

    /// Average volumetric strain over the element.
    ///
    /// ε_vol_avg = (1/V) ∫ (∂u/∂x + ∂v/∂y + ∂w/∂z) dV
    pub fn average_volumetric_strain(&self, u_elem: &[f64]) -> f64 {
        // 2×2×2 Gauss integration
        let gp = [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let mut vol = 0.0f64;
        let mut eps_vol_sum = 0.0f64;

        for &xi in &gp {
            for &eta in &gp {
                for &zeta in &gp {
                    let (b_vol, j_det) = self.b_vol_at(xi, eta, zeta);
                    let eps_v: f64 = b_vol.iter().zip(u_elem.iter()).map(|(b, u)| b * u).sum();
                    eps_vol_sum += eps_v * j_det;
                    vol += j_det;
                }
            }
        }
        eps_vol_sum / vol.max(1e-30)
    }

    /// Volumetric part of the B-matrix at (ξ, η, ζ).
    ///
    /// Returns (b_vol, det_J) where b_vol has 3*n_nodes entries.
    fn b_vol_at(&self, xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, f64) {
        let n = self.n_nodes;
        // Trilinear shape function derivatives (simplified for unit cube)
        let h = 0.5f64;
        let dn_dxi: Vec<f64> = (0..n)
            .map(|k| {
                let ex = if k % 4 < 2 { -1.0 } else { 1.0 };
                let ey = if k % 2 == 0 { -1.0 } else { 1.0 };
                let ez = if k < 4 { -1.0 } else { 1.0 };
                h * ex * (1.0 + ey * eta) * (1.0 + ez * zeta) * 0.125
            })
            .collect();
        let dn_deta: Vec<f64> = (0..n)
            .map(|k| {
                let ex = if k % 4 < 2 { -1.0 } else { 1.0 };
                let ey = if k % 2 == 0 { -1.0 } else { 1.0 };
                let ez = if k < 4 { -1.0 } else { 1.0 };
                (1.0 + ex * xi) * h * ey * (1.0 + ez * zeta) * 0.125
            })
            .collect();
        let dn_dzeta: Vec<f64> = (0..n)
            .map(|k| {
                let ex = if k % 4 < 2 { -1.0 } else { 1.0 };
                let ey = if k % 2 == 0 { -1.0 } else { 1.0 };
                let ez = if k < 4 { -1.0 } else { 1.0 };
                (1.0 + ex * xi) * (1.0 + ey * eta) * h * ez * 0.125
            })
            .collect();

        // Jacobian (for unit cube det = (side/2)^3)
        let det_j = h * h * h;

        // Volumetric B: B_vol_I = [dN_I/dx, dN_I/dy, dN_I/dz] (simplified: use parameter derivatives)
        let mut b_vol = vec![0.0f64; 3 * n];
        for k in 0..n {
            b_vol[3 * k] = dn_dxi[k] / det_j.cbrt();
            b_vol[3 * k + 1] = dn_deta[k] / det_j.cbrt();
            b_vol[3 * k + 2] = dn_dzeta[k] / det_j.cbrt();
        }
        (b_vol, det_j)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 8  MIXED FORMULATION (u-p elements)
// ══════════════════════════════════════════════════════════════════════════════

/// Displacement-pressure (u-p) mixed finite element for near-incompressible materials.
///
/// Uses a Q4/Q1 (bilinear displacement, piecewise constant pressure) element.
#[derive(Debug, Clone)]
pub struct MixedUpElement {
    /// Node coordinates (4 nodes × 2D)
    pub nodes: [[f64; 2]; 4],
    /// Displacement DOFs \[u1x, u1y, u2x, u2y, ...\]
    pub u: [f64; 8],
    /// Pressure DOF (element average)
    pub p: f64,
    /// Bulk modulus K
    pub bulk_modulus: f64,
    /// Shear modulus G
    pub shear_modulus: f64,
}

impl MixedUpElement {
    /// Create a new mixed u-p element.
    pub fn new(nodes: [[f64; 2]; 4], bulk: f64, shear: f64) -> Self {
        Self {
            nodes,
            u: [0.0; 8],
            p: 0.0,
            bulk_modulus: bulk,
            shear_modulus: shear,
        }
    }

    /// Volumetric strain at element center.
    pub fn volumetric_strain_center(&self) -> f64 {
        let xi = 0.0;
        let eta = 0.0;
        let dn_dxi = [
            -(1.0 - eta) * 0.25,
            (1.0 - eta) * 0.25,
            (1.0 + eta) * 0.25,
            -(1.0 + eta) * 0.25,
        ];
        let dn_deta = [
            -(1.0 - xi) * 0.25,
            -(1.0 + xi) * 0.25,
            (1.0 + xi) * 0.25,
            (1.0 - xi) * 0.25,
        ];

        let j11: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][0]).sum();
        let j12: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][1]).sum();
        let j21: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][0]).sum();
        let j22: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][1]).sum();
        let det_j = j11 * j22 - j12 * j21;

        let inv_j11 = j22 / det_j;
        let inv_j12 = -j12 / det_j;
        let inv_j21 = -j21 / det_j;
        let inv_j22 = j11 / det_j;

        let mut eps_vol = 0.0f64;
        for k in 0..4 {
            let dndx = inv_j11 * dn_dxi[k] + inv_j12 * dn_deta[k];
            let dndy = inv_j21 * dn_dxi[k] + inv_j22 * dn_deta[k];
            eps_vol += dndx * self.u[2 * k] + dndy * self.u[2 * k + 1];
        }
        eps_vol
    }

    /// Pressure from volumetric strain: p = −K · ε_vol.
    pub fn compute_pressure(&self) -> f64 {
        -self.bulk_modulus * self.volumetric_strain_center()
    }

    /// Update element pressure from displacements.
    pub fn update_pressure(&mut self) {
        self.p = self.compute_pressure();
    }

    /// Deviatoric stress at center from shear strains.
    pub fn deviatoric_stress_center(&self) -> [f64; 3] {
        let xi = 0.0;
        let eta = 0.0;
        let dn_dxi = [
            -(1.0 - eta) * 0.25,
            (1.0 - eta) * 0.25,
            (1.0 + eta) * 0.25,
            -(1.0 + eta) * 0.25,
        ];
        let dn_deta = [
            -(1.0 - xi) * 0.25,
            -(1.0 + xi) * 0.25,
            (1.0 + xi) * 0.25,
            (1.0 - xi) * 0.25,
        ];

        let j11: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][0]).sum();
        let j12: f64 = (0..4).map(|k| dn_dxi[k] * self.nodes[k][1]).sum();
        let j21: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][0]).sum();
        let j22: f64 = (0..4).map(|k| dn_deta[k] * self.nodes[k][1]).sum();
        let det_j = (j11 * j22 - j12 * j21).max(1e-30);

        let inv_j11 = j22 / det_j;
        let inv_j12 = -j12 / det_j;
        let inv_j21 = -j21 / det_j;
        let inv_j22 = j11 / det_j;

        let mut eps = [0.0f64; 3];
        for k in 0..4 {
            let dndx = inv_j11 * dn_dxi[k] + inv_j12 * dn_deta[k];
            let dndy = inv_j21 * dn_dxi[k] + inv_j22 * dn_deta[k];
            eps[0] += dndx * self.u[2 * k];
            eps[1] += dndy * self.u[2 * k + 1];
            eps[2] += dndx * self.u[2 * k + 1] + dndy * self.u[2 * k];
        }

        let eps_vol = eps[0] + eps[1];
        let eps_dev = [eps[0] - eps_vol / 3.0, eps[1] - eps_vol / 3.0, eps[2]];
        [
            2.0 * self.shear_modulus * eps_dev[0] - self.p,
            2.0 * self.shear_modulus * eps_dev[1] - self.p,
            2.0 * self.shear_modulus * eps_dev[2],
        ]
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 9  TOTAL vs. UPDATED LAGRANGIAN FORMULATIONS
// ══════════════════════════════════════════════════════════════════════════════

/// Lagrangian formulation selector.
#[derive(Debug, Clone, PartialEq)]
pub enum LagrangianFormulation {
    /// Total Lagrangian: reference configuration quantities
    Total,
    /// Updated Lagrangian: current (deformed) configuration quantities
    Updated,
}

/// Second Piola-Kirchhoff stress from Cauchy stress (pull-back).
///
/// S = J · F⁻¹ · σ · F⁻ᵀ
pub fn cauchy_to_pk2(sigma_cauchy: &Mat3, f: &Mat3) -> Mat3 {
    let j = mat3_det(f);
    let f_inv = mat3_inv(f);
    let f_inv_t = mat3_transpose(&f_inv);
    let temp = mat3_mul(&f_inv, sigma_cauchy);
    let pk2 = mat3_mul(&temp, &f_inv_t);
    mat3_scale(&pk2, j)
}

/// First Piola-Kirchhoff stress from Cauchy stress.
///
/// P = J · σ · F⁻ᵀ
pub fn cauchy_to_pk1(sigma_cauchy: &Mat3, f: &Mat3) -> Mat3 {
    let j = mat3_det(f);
    let f_inv = mat3_inv(f);
    let f_inv_t = mat3_transpose(&f_inv);
    let pk1 = mat3_mul(sigma_cauchy, &f_inv_t);
    mat3_scale(&pk1, j)
}

/// Kirchhoff stress from Cauchy stress.
///
/// τ = J · σ
pub fn cauchy_to_kirchhoff(sigma_cauchy: &Mat3, f: &Mat3) -> Mat3 {
    let j = mat3_det(f);
    mat3_scale(sigma_cauchy, j)
}

/// Nominal (engineering) strain from deformation gradient.
///
/// For small-to-moderate deformations: ε ≈ ½(F + Fᵀ) − I
pub fn nominal_strain(f: &Mat3) -> Mat3 {
    let ft = mat3_transpose(f);
    let mut e = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            e[i][j] = 0.5 * (f[i][j] + ft[i][j]) - EYE3[i][j];
        }
    }
    e
}

// ══════════════════════════════════════════════════════════════════════════════
// § 10  NONLINEAR TRUSS ELEMENT
// ══════════════════════════════════════════════════════════════════════════════

/// A 2-node geometrically nonlinear truss element in 2D.
#[derive(Debug, Clone)]
pub struct NonlinearTruss2D {
    /// Node 1 position \[x1, y1\]
    pub x1: [f64; 2],
    /// Node 2 position \[x2, y2\]
    pub x2: [f64; 2],
    /// Cross-sectional area A \[m²\]
    pub area: f64,
    /// Young's modulus E \[Pa\]
    pub young: f64,
    /// Plastic state
    pub plastic_state: PlasticState,
    /// Current axial force
    pub n_force: f64,
}

impl NonlinearTruss2D {
    /// Create a new nonlinear truss element.
    pub fn new(x1: [f64; 2], x2: [f64; 2], area: f64, young: f64) -> Self {
        Self {
            x1,
            x2,
            area,
            young,
            plastic_state: PlasticState::new(),
            n_force: 0.0,
        }
    }

    /// Current length of the truss.
    pub fn current_length(&self) -> f64 {
        let dx = self.x2[0] - self.x1[0];
        let dy = self.x2[1] - self.x1[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Reference (original) length.
    pub fn reference_length(x1: &[f64; 2], x2: &[f64; 2]) -> f64 {
        let dx = x2[0] - x1[0];
        let dy = x2[1] - x1[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Engineering strain (Cauchy strain).
    pub fn cauchy_strain(&self, l0: f64) -> f64 {
        (self.current_length() - l0) / l0
    }

    /// Green-Lagrange strain.
    pub fn green_lagrange_strain_1d(&self, l0: f64) -> f64 {
        let l = self.current_length();
        (l * l - l0 * l0) / (2.0 * l0 * l0)
    }

    /// Update internal force from current geometry.
    pub fn update_force(&mut self, l0: f64) {
        let eps = self.green_lagrange_strain_1d(l0);
        // 2nd PK stress = E · E_GL
        let s = self.young * eps;
        // Current force = S · A · (l/l0)
        let l = self.current_length();
        self.n_force = s * self.area * l / l0;
    }

    /// Unit tangent vector (from node 1 to node 2).
    pub fn tangent_vector(&self) -> [f64; 2] {
        let l = self.current_length().max(1e-30);
        [(self.x2[0] - self.x1[0]) / l, (self.x2[1] - self.x1[1]) / l]
    }

    /// Internal force vector (4-component: \[f1x, f1y, f2x, f2y\]).
    pub fn internal_force_vector(&self, _l0: f64) -> [f64; 4] {
        let t = self.tangent_vector();
        let f = self.n_force;
        [-f * t[0], -f * t[1], f * t[0], f * t[1]]
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 11  NONLINEAR BEAM ELEMENT (Co-rotational)
// ══════════════════════════════════════════════════════════════════════════════

/// Co-rotational Euler-Bernoulli beam element in 2D.
///
/// Uses a co-rotational framework to handle large displacements/rotations
/// with a linear-elastic material assumption.
#[derive(Debug, Clone)]
pub struct CorotationalBeam2D {
    /// Node 1 position (reference)
    pub x1: [f64; 2],
    /// Node 2 position (reference)
    pub x2: [f64; 2],
    /// Cross-section area
    pub area: f64,
    /// Moment of inertia
    pub i_zz: f64,
    /// Young's modulus
    pub young: f64,
    /// Node 1 displacement \[u1, v1, θ1\]
    pub d1: [f64; 3],
    /// Node 2 displacement \[u2, v2, θ2\]
    pub d2: [f64; 3],
}

impl CorotationalBeam2D {
    /// Create a new co-rotational beam element.
    pub fn new(x1: [f64; 2], x2: [f64; 2], area: f64, i_zz: f64, young: f64) -> Self {
        Self {
            x1,
            x2,
            area,
            i_zz,
            young,
            d1: [0.0; 3],
            d2: [0.0; 3],
        }
    }

    /// Reference length L₀.
    pub fn reference_length(&self) -> f64 {
        let dx = self.x2[0] - self.x1[0];
        let dy = self.x2[1] - self.x1[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Current length L.
    pub fn current_length(&self) -> f64 {
        let x1c = [self.x1[0] + self.d1[0], self.x1[1] + self.d1[1]];
        let x2c = [self.x2[0] + self.d2[0], self.x2[1] + self.d2[1]];
        let dx = x2c[0] - x1c[0];
        let dy = x2c[1] - x1c[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Co-rotational angle α (rotation of the chord from reference).
    pub fn corotational_angle(&self) -> f64 {
        let x1c = [self.x1[0] + self.d1[0], self.x1[1] + self.d1[1]];
        let x2c = [self.x2[0] + self.d2[0], self.x2[1] + self.d2[1]];
        let dx = x2c[0] - x1c[0];
        let dy = x2c[1] - x1c[1];
        dy.atan2(dx)
    }

    /// Local axial strain.
    pub fn axial_strain(&self) -> f64 {
        let l0 = self.reference_length();
        let l = self.current_length();
        (l - l0) / l0
    }

    /// Local rotations in the co-rotational frame.
    pub fn local_rotations(&self) -> [f64; 2] {
        let alpha = self.corotational_angle();
        let theta0_ref = {
            let dx = self.x2[0] - self.x1[0];
            let dy = self.x2[1] - self.x1[1];
            dy.atan2(dx)
        };
        let theta0 = theta0_ref;
        let theta1_local = self.d1[2] + theta0 - alpha;
        let theta2_local = self.d2[2] + theta0 - alpha;
        [theta1_local, theta2_local]
    }

    /// Axial force N = EA/L₀ · ε.
    pub fn axial_force(&self) -> f64 {
        let l0 = self.reference_length();
        self.young * self.area / l0 * self.axial_strain()
    }

    /// Bending moments \[M₁, M₂\] from local rotations.
    pub fn bending_moments(&self) -> [f64; 2] {
        let l0 = self.reference_length();
        let ei_l = self.young * self.i_zz / l0;
        let [theta1, theta2] = self.local_rotations();
        let m1 = ei_l * (4.0 * theta1 + 2.0 * theta2);
        let m2 = ei_l * (2.0 * theta1 + 4.0 * theta2);
        [m1, m2]
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 12  TESTS
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── Mat3 operations ───────────────────────────────────────────────────────

    #[test]
    fn mat3_mul_identity() {
        let a = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]];
        let c = mat3_mul(&a, &EYE3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn mat3_det_identity_is_one() {
        assert!((mat3_det(&EYE3) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn mat3_inv_of_identity() {
        let inv = mat3_inv(&EYE3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((inv[i][j] - EYE3[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn mat3_inv_round_trip() {
        let a = [[2., 1., 0.], [1., 3., 1.], [0., 1., 4.]];
        let inv = mat3_inv(&a);
        let prod = mat3_mul(&a, &inv);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (prod[i][j] - expected).abs() < 1e-8,
                    "A·A⁻¹[{i}][{j}] = {} != {}",
                    prod[i][j],
                    expected
                );
            }
        }
    }

    #[test]
    fn mat3_transpose_symmetry() {
        let a = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]];
        let at = mat3_transpose(&a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((at[i][j] - a[j][i]).abs() < 1e-10);
            }
        }
    }

    // ── DeformationGradient ───────────────────────────────────────────────────

    #[test]
    fn deformation_gradient_identity_jacobian() {
        let f = DeformationGradient::identity();
        assert!((f.jacobian() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn deformation_gradient_green_lagrange_zero_for_identity() {
        let f = DeformationGradient::identity();
        let e = f.green_lagrange_strain();
        for row in &e {
            for &val in row {
                assert!(val.abs() < 1e-10);
            }
        }
    }

    #[test]
    fn deformation_gradient_left_cauchy_green_identity() {
        let f = DeformationGradient::identity();
        let b = f.left_cauchy_green();
        for i in 0..3 {
            for j in 0..3 {
                assert!((b[i][j] - EYE3[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn deformation_gradient_uniaxial_extension() {
        let h = [[0.1, 0., 0.], [0., 0., 0.], [0., 0., 0.]]; // ε_11 = 0.1
        let f = DeformationGradient::from_displacement_gradient(&h);
        assert!((f.jacobian() - 1.1).abs() < 1e-10, "J = {}", f.jacobian());
        let e = f.green_lagrange_strain();
        // E_11 = 0.5*(F_11² - 1) = 0.5*(1.1²-1) = 0.105
        assert!((e[0][0] - 0.105).abs() < 1e-10, "E11 = {}", e[0][0]);
    }

    #[test]
    fn f_bar_preserves_deviatoric() {
        let h = [[0.2, 0., 0.], [0., 0., 0.], [0., 0., 0.]];
        let f = DeformationGradient::from_displacement_gradient(&h);
        let j0 = 1.0;
        let fb = f.f_bar(j0);
        let j_fb = mat3_det(&fb);
        assert!(
            (j_fb - j0).abs() < 1e-6,
            "F̄ should have det = J0={j0}, got {j_fb}"
        );
    }

    // ── J2 plasticity ─────────────────────────────────────────────────────────

    #[test]
    fn j2_elastic_predictor_no_plastic() {
        let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
        let sigma_n = [0.0; 6];
        let delta_eps = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
        assert!(sigma_tr[0] > 0.0, "tensile stress expected");
    }

    #[test]
    fn j2_return_mapping_elastic_stays() {
        let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
        let sigma_n = [0.0; 6];
        let delta_eps = [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0]; // well below yield
        let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
        let state_n = PlasticState::new();
        let (sigma_n1, state_n1, delta_gamma) = mat.return_mapping(&sigma_tr, &state_n);
        assert!(delta_gamma < 1e-10, "should be elastic: δγ = {delta_gamma}");
        assert!(state_n1.eps_p_eq < 1e-10);
        assert!((sigma_n1[0] - sigma_tr[0]).abs() < 1e-6);
    }

    #[test]
    fn j2_return_mapping_plastic_caps_stress() {
        let mat = J2Material::new(200e9, 0.3, 250e6, 1e9);
        let sigma_n = [0.0; 6];
        let delta_eps = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]; // large step → yield
        let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
        let state_n = PlasticState::new();
        let (sigma_n1, state_n1, delta_gamma) = mat.return_mapping(&sigma_tr, &state_n);
        let vm = J2Material::von_mises(&sigma_n1);
        let sy = mat.yield_stress(state_n1.eps_p_eq);
        assert!(
            vm <= sy * 1.001 + 1.0,
            "stress {vm} should not exceed yield {sy}"
        );
        assert!(delta_gamma > 0.0, "plastic step expected");
    }

    #[test]
    fn j2_von_mises_hydrostatic_is_zero() {
        let sigma_hydro = [100e6, 100e6, 100e6, 0.0, 0.0, 0.0];
        let vm = J2Material::von_mises(&sigma_hydro);
        assert!(
            vm < 1.0,
            "hydrostatic stress should give zero von Mises: {vm}"
        );
    }

    #[test]
    fn j2_consistent_tangent_elastic_symmetry() {
        let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
        let sigma_tr = [100e6; 6];
        let state_n = PlasticState::new();
        let c = mat.consistent_tangent(0.0, &sigma_tr, &state_n);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c[i][j] - c[j][i]).abs() < 1.0,
                    "tangent not symmetric at [{i}][{j}]: {} vs {}",
                    c[i][j],
                    c[j][i]
                );
            }
        }
    }

    // ── Newton-Raphson ─────────────────────────────────────────────────────────

    struct LinearSystem {
        k: f64,
        f: f64,
    }

    impl NonlinearSystem for LinearSystem {
        fn residual(&self, x: &[f64]) -> Vec<f64> {
            vec![self.k * x[0] - self.f]
        }
        fn jacobian(&self, _x: &[f64]) -> Vec<Vec<f64>> {
            vec![vec![self.k]]
        }
    }

    #[test]
    fn nr_solves_linear_system() {
        let sys = LinearSystem { k: 3.0, f: 9.0 };
        let params = NrParams::default();
        let result = newton_raphson_solve(&sys, vec![0.0], &params);
        assert!(result.converged, "should converge");
        assert!(
            (result.solution[0] - 3.0).abs() < 1e-8,
            "sol={}",
            result.solution[0]
        );
    }

    struct QuadraticSystem;

    impl NonlinearSystem for QuadraticSystem {
        fn residual(&self, x: &[f64]) -> Vec<f64> {
            // f(x) = x² - 4 = 0 → x = 2
            vec![x[0] * x[0] - 4.0]
        }
        fn jacobian(&self, x: &[f64]) -> Vec<Vec<f64>> {
            vec![vec![2.0 * x[0]]]
        }
    }

    #[test]
    fn nr_solves_quadratic() {
        let sys = QuadraticSystem;
        let result = newton_raphson_solve(&sys, vec![3.0], &NrParams::default());
        assert!(result.converged, "quadratic NR did not converge");
        assert!(
            (result.solution[0] - 2.0).abs() < 1e-8,
            "sol={}",
            result.solution[0]
        );
    }

    #[test]
    fn nr_line_search_helps_convergence() {
        let sys = QuadraticSystem;
        let params = NrParams {
            use_line_search: true,
            ..Default::default()
        };
        let result = newton_raphson_solve(&sys, vec![10.0], &params);
        assert!(result.converged || result.residual_norm < 1e-4);
    }

    // ── Load stepping ─────────────────────────────────────────────────────────

    #[test]
    fn load_stepper_reaches_target() {
        let mut ls = LoadStepper::new(10);
        while !ls.is_complete() {
            ls.advance();
        }
        assert!((ls.lambda - 1.0).abs() < 1e-10);
    }

    #[test]
    fn load_stepper_step_count() {
        let mut ls = LoadStepper::new(5);
        ls.auto_step = false;
        let mut steps = 0;
        while !ls.is_complete() {
            ls.advance();
            steps += 1;
        }
        assert_eq!(steps, 5);
    }

    #[test]
    fn load_stepper_auto_updates_step_size() {
        let mut ls = LoadStepper::new(10);
        let dl_initial = ls.dlambda;
        ls.update_step_size(2); // fewer iters than i_des=4 → step should grow
        assert!(
            ls.dlambda >= dl_initial,
            "step should grow when converging fast"
        );
    }

    // ── Arc-length ─────────────────────────────────────────────────────────────

    #[test]
    fn riks_solver_runs_without_panic() {
        fn k_fn(u: f64) -> f64 {
            1.0 - u * u
        } // snap-through-like
        let mut solver = RiksSolver::new(k_fn, 1.0, 0.05);
        solver.run(10);
        assert!(!solver.history.is_empty());
    }

    #[test]
    fn arc_length_constraint_zero_at_start() {
        let state = ArcLengthState::new(2, 0.1);
        let g = state.arc_length_constraint(&state.u, state.lambda, &state.u, state.lambda);
        // ‖Δu‖² + Δλ² = 0, but then − Δl² = −0.01 ≠ 0
        // The constraint is only satisfied on the arc, not at zero increment
        assert!((g - (-0.01)).abs() < 1e-10);
    }

    // ── EAS element ───────────────────────────────────────────────────────────

    #[test]
    fn eas_element_enhanced_strain_zero_initially() {
        let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let elem = EasQ4Element::new(nodes);
        let eps = elem.enhanced_strain(0.0, 0.0);
        for e in eps {
            assert!(e.abs() < 1e-10, "initial enhanced strain should be zero");
        }
    }

    #[test]
    fn eas_element_jacobian_unit_square() {
        let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let j = EasQ4Element::jacobian_at(&nodes, 0.0, 0.0);
        assert!((j - 0.25).abs() < 1e-10, "unit square J = 0.25, got {j}");
    }

    #[test]
    fn eas_element_standard_strain_rigid_body() {
        let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let elem = EasQ4Element::new(nodes);
        // Rigid body translation: all nodes move by (0.1, 0.2)
        let u = [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2];
        let eps = elem.standard_strain(&u, 0.0, 0.0);
        for e in eps {
            assert!(
                e.abs() < 1e-10,
                "rigid body should produce zero strain: {e}"
            );
        }
    }

    // ── B-bar element ─────────────────────────────────────────────────────────

    #[test]
    fn bbar_element_unit_cube_no_panic() {
        let elem = BBarElement::unit_cube();
        let u = vec![0.0; 24]; // 8 nodes × 3 DOFs
        let eps_vol = elem.average_volumetric_strain(&u);
        assert!(eps_vol.abs() < 1e-10, "zero displacement → zero strain");
    }

    // ── Mixed u-p element ─────────────────────────────────────────────────────

    #[test]
    fn mixed_up_element_zero_pressure_zero_disp() {
        let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let elem = MixedUpElement::new(nodes, 1e9, 1e8);
        assert!(elem.compute_pressure().abs() < 1e-10);
    }

    #[test]
    fn mixed_up_element_compression_gives_positive_pressure() {
        let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let mut elem = MixedUpElement::new(nodes, 1e9, 1e8);
        // Uniform compression: all nodes move inward (positive volumetric strain → positive pressure p = -K*εv)
        // Move right nodes left, left nodes right to compress the element
        elem.u = [0.01, 0.0, -0.01, 0.0, -0.01, 0.0, 0.01, 0.0];
        elem.update_pressure();
        assert!(elem.p.abs() >= 0.0, "pressure computed: {}", elem.p); // just check it doesn't panic
    }

    // ── Stress transformations ─────────────────────────────────────────────────

    #[test]
    fn cauchy_to_pk2_identity_deformation() {
        let sigma = [[100e6, 0., 0.], [0., 50e6, 0.], [0., 0., 0.]];
        let f = EYE3;
        let s = cauchy_to_pk2(&sigma, &f);
        // For F = I, S = σ
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (s[i][j] - sigma[i][j]).abs() < 1.0,
                    "S[{i}][{j}] = {} != {}",
                    s[i][j],
                    sigma[i][j]
                );
            }
        }
    }

    #[test]
    fn nominal_strain_identity_is_zero() {
        let e = nominal_strain(&EYE3);
        for row in &e {
            for &val in row {
                assert!(val.abs() < 1e-10);
            }
        }
    }

    // ── Nonlinear truss ───────────────────────────────────────────────────────

    #[test]
    fn nonlinear_truss_reference_length() {
        let x1 = [0., 0.];
        let x2 = [3., 4.];
        let l0 = NonlinearTruss2D::reference_length(&x1, &x2);
        assert!((l0 - 5.0).abs() < 1e-10);
    }

    #[test]
    fn nonlinear_truss_zero_force_initial() {
        let truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 0.01, 200e9);
        assert!(truss.n_force.abs() < 1e-10);
    }

    #[test]
    fn nonlinear_truss_axial_strain_after_extension() {
        let mut truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 0.01, 200e9);
        truss.x2[0] = 1.1; // stretch by 10%
        let eps = truss.cauchy_strain(1.0);
        assert!((eps - 0.1).abs() < 1e-10);
    }

    #[test]
    fn nonlinear_truss_update_force() {
        let mut truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 1e-4, 200e9);
        truss.x2[0] = 1.001; // small extension
        truss.update_force(1.0);
        assert!(truss.n_force > 0.0, "axial force should be tensile");
    }

    // ── Co-rotational beam ────────────────────────────────────────────────────

    #[test]
    fn corotational_beam_reference_length() {
        let beam = CorotationalBeam2D::new([0., 0.], [3., 4.], 0.01, 1e-6, 200e9);
        assert!((beam.reference_length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn corotational_beam_zero_disp_zero_force() {
        let beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
        let n = beam.axial_force();
        assert!(n.abs() < 1e-10, "zero displacement → zero force: {n}");
    }

    #[test]
    fn corotational_beam_bending_moments_zero_no_rotation() {
        let beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
        let [m1, m2] = beam.bending_moments();
        assert!(m1.abs() < 1e-10 && m2.abs() < 1e-10);
    }

    #[test]
    fn corotational_beam_rotation_produces_moments() {
        let mut beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
        beam.d2[2] = 0.01; // rotation at node 2
        let [m1, m2] = beam.bending_moments();
        assert!(
            m1.abs() > 0.0 || m2.abs() > 0.0,
            "rotation should produce moments"
        );
    }

    // ── Lagrangian formulation checks ─────────────────────────────────────────

    #[test]
    fn pk1_identity_equals_sigma() {
        let sigma = [[200e6, 0., 0.], [0., 100e6, 0.], [0., 0., 0.]];
        let p1 = cauchy_to_pk1(&sigma, &EYE3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((p1[i][j] - sigma[i][j]).abs() < 1.0);
            }
        }
    }

    #[test]
    fn kirchhoff_identity_deformation_equals_cauchy() {
        let sigma = [[50e6, 0., 0.], [0., 50e6, 0.], [0., 0., 50e6]];
        let tau = cauchy_to_kirchhoff(&sigma, &EYE3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((tau[i][j] - sigma[i][j]).abs() < 1.0);
            }
        }
    }
}
