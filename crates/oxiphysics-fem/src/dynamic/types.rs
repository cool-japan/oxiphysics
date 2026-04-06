//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    compute_rayleigh_coefficients, modal_truncation_correction, rayleigh_damping_matrix,
};
#[allow(unused_imports)]
use super::functions_2::*;
use crate::sparse::CsrMatrix;

/// State for the Bathe two-sub-step implicit time integrator.
///
/// The Bathe method splits each time step Δt into two sub-steps (Δt/2 each)
/// using the trapezoidal rule for the first sub-step and the 3-point Euler
/// backward scheme for the second.  This gives second-order accuracy with
/// effective numerical dissipation of spurious high-frequency content.
///
/// Reference: Bathe & Baig (2005) "On a composite implicit time integration
/// procedure for nonlinear dynamics", *Computers & Structures* 83, 2513–2524.
#[allow(dead_code)]
pub struct BatheIntegrator {
    /// Number of DOF.
    pub n_dof: usize,
    /// Displacement at end of last full step.
    pub u: Vec<f64>,
    /// Velocity at end of last full step.
    pub v: Vec<f64>,
    /// Acceleration at end of last full step.
    pub a: Vec<f64>,
}
impl BatheIntegrator {
    /// Create a new Bathe integrator with zero initial conditions.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
        }
    }
    /// Advance one full step Δt using the two-sub-step Bathe scheme.
    ///
    /// Solves `M a + C v + K u = f` (no full matrix solve – diagonal mass/stiffness).
    ///
    /// For a diagonal (lumped) mass and diagonal stiffness this reduces to N
    /// independent 2×2 systems, making it suitable for demonstration.
    ///
    /// # Arguments
    /// * `m_diag`  – diagonal mass vector (length n_dof)
    /// * `k_diag`  – diagonal stiffness vector (length n_dof)
    /// * `c_diag`  – diagonal damping vector (length n_dof)
    /// * `f`       – force vector at time t + Δt
    /// * `dt`      – full step size Δt
    #[allow(clippy::too_many_arguments)]
    pub fn step_diagonal(
        &mut self,
        m_diag: &[f64],
        k_diag: &[f64],
        c_diag: &[f64],
        f: &[f64],
        dt: f64,
    ) {
        let h = dt / 2.0;
        let n = self.n_dof;
        let mut u_mid = vec![0.0_f64; n];
        let mut v_mid = vec![0.0_f64; n];
        let mut a_mid = vec![0.0_f64; n];
        for i in 0..n {
            let k_eff = 4.0 / (h * h) * m_diag[i] + 2.0 / h * c_diag[i] + k_diag[i];
            let f_eff = m_diag[i] * (4.0 / (h * h) * self.u[i] + 4.0 / h * self.v[i] + self.a[i])
                + c_diag[i] * (2.0 / h * self.u[i] + self.v[i])
                + f[i];
            if k_eff.abs() < 1e-60 {
                u_mid[i] = self.u[i];
                v_mid[i] = self.v[i];
                a_mid[i] = self.a[i];
            } else {
                u_mid[i] = f_eff / k_eff;
                v_mid[i] = 2.0 / h * (u_mid[i] - self.u[i]) - self.v[i];
                a_mid[i] = 2.0 / h * (v_mid[i] - self.v[i]) - self.a[i];
            }
        }
        let g0 = 4.0 / (h * h);
        let g1 = 2.0 / h;
        for i in 0..n {
            let k_eff = g0 * m_diag[i] + g1 * c_diag[i] + k_diag[i];
            let f_eff2 = m_diag[i] * (g0 * u_mid[i] + g1 * v_mid[i] + a_mid[i])
                + c_diag[i] * (g1 * u_mid[i] + v_mid[i])
                + f[i];
            if k_eff.abs() < 1e-60 {
            } else {
                let u_new = f_eff2 / k_eff;
                let v_new = g1 * (u_new - u_mid[i]) - v_mid[i];
                let a_new = g0 * (u_new - u_mid[i]) - g1 * v_mid[i] - a_mid[i];
                self.u[i] = u_new;
                self.v[i] = v_new;
                self.a[i] = a_new;
            }
        }
    }
}
/// Newmark linear acceleration parameters (β = 1/6, γ = 1/2).
///
/// Conditionally stable: Δt ≤ √3 / ω_max.
#[allow(dead_code)]
pub struct NewmarkLinearAcceleration;
impl NewmarkLinearAcceleration {
    /// Single SDOF step using linear acceleration assumption.
    ///
    /// Returns (u_new, v_new, a_new).
    pub fn step(
        u: f64,
        v: f64,
        a: f64,
        f_new: f64,
        omega: f64,
        zeta: f64,
        dt: f64,
    ) -> (f64, f64, f64) {
        let beta = 1.0 / 6.0_f64;
        let gamma = 0.5_f64;
        let k_eff =
            omega * omega + 2.0 * zeta * omega * gamma / (beta * dt) + 1.0 / (beta * dt * dt);
        let f_eff = f_new
            + u / (beta * dt * dt)
            + v / (beta * dt)
            + (1.0 / (2.0 * beta) - 1.0) * a
            + 2.0
                * zeta
                * omega
                * (gamma / (beta * dt) * u
                    + (gamma / beta - 1.0) * v
                    + dt * (gamma / (2.0 * beta) - 1.0) * a);
        let u_new = if k_eff.abs() < 1e-60 {
            u
        } else {
            f_eff / k_eff
        };
        let v_new = gamma / (beta * dt) * (u_new - u)
            - (gamma / beta - 1.0) * v
            - dt * (gamma / (2.0 * beta) - 1.0) * a;
        let a_new =
            (u_new - u) / (beta * dt * dt) - v / (beta * dt) - (1.0 / (2.0 * beta) - 1.0) * a;
        (u_new, v_new, a_new)
    }
}
/// State for the central difference explicit time integrator.
///
/// This is a conditionally stable explicit method requiring:
/// `dt < 2/omega_max` where `omega_max` is the highest natural frequency.
///
/// The update equations are:
/// ```text
/// a_n = M^{-1} * (F_n - K * u_n)
/// u_{n+1} = 2*u_n - u_{n-1} + dt^2 * a_n
/// v_n = (u_{n+1} - u_{n-1}) / (2*dt)
/// ```
///
/// Requires a diagonal (lumped) mass matrix for efficiency.
#[allow(dead_code)]
pub struct CentralDifferenceState {
    /// Current displacements.
    pub u: Vec<f64>,
    /// Previous displacements (u_{n-1}).
    pub u_prev: Vec<f64>,
    /// Current velocities (computed from central difference).
    pub v: Vec<f64>,
    /// Current accelerations.
    pub a: Vec<f64>,
    /// Diagonal of the lumped mass matrix (inverse).
    pub inv_m_diag: Vec<f64>,
}
#[allow(dead_code)]
impl CentralDifferenceState {
    /// Create a new central difference state with zero initial conditions.
    ///
    /// `inv_m_diag` is the diagonal of the inverse lumped mass matrix.
    pub fn new(n_dof: usize, inv_m_diag: Vec<f64>) -> Self {
        Self {
            u: vec![0.0; n_dof],
            u_prev: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
            inv_m_diag,
        }
    }
    /// Advance the state by one explicit time step.
    pub fn step(&mut self, stiffness: &CsrMatrix, force: &[f64], dt: f64, fixed_dofs: &[usize]) {
        let n = self.u.len();
        let f_int = stiffness.mul_vec(&self.u);
        for i in 0..n {
            self.a[i] = self.inv_m_diag[i] * (force[i] - f_int[i]);
        }
        let dt2 = dt * dt;
        let mut u_new = vec![0.0; n];
        for i in 0..n {
            u_new[i] = 2.0 * self.u[i] - self.u_prev[i] + dt2 * self.a[i];
        }
        let inv_2dt = 0.5 / dt;
        for i in 0..n {
            self.v[i] = (u_new[i] - self.u_prev[i]) * inv_2dt;
        }
        for &dof in fixed_dofs {
            u_new[dof] = 0.0;
            self.v[dof] = 0.0;
            self.a[dof] = 0.0;
        }
        self.u_prev = self.u.clone();
        self.u = u_new;
    }
}
/// State for the generalized-alpha time integrator.
///
/// Provides user-controlled high-frequency dissipation through the
/// spectral radius parameter `rho_inf` in `[0, 1]`:
/// - `rho_inf = 1.0`: no dissipation (equivalent to Newmark trapezoidal rule)
/// - `rho_inf = 0.0`: asymptotic annihilation of highest mode
#[allow(dead_code)]
pub struct GeneralizedAlphaState {
    /// Current displacements.
    pub u: Vec<f64>,
    /// Current velocities.
    pub v: Vec<f64>,
    /// Current accelerations.
    pub a: Vec<f64>,
    /// Alpha_m parameter.
    pub alpha_m: f64,
    /// Alpha_f parameter.
    pub alpha_f: f64,
    /// Newmark beta.
    pub beta: f64,
    /// Newmark gamma.
    pub gamma: f64,
}
#[allow(dead_code)]
impl GeneralizedAlphaState {
    /// Create a new generalized-alpha state from the spectral radius
    /// at infinity `rho_inf` in `[0, 1]`.
    pub fn new(n_dof: usize, rho_inf: f64) -> Self {
        let rho = rho_inf.clamp(0.0, 1.0);
        let alpha_m = (2.0 * rho - 1.0) / (rho + 1.0);
        let alpha_f = rho / (rho + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f).powi(2);
        Self {
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
            alpha_m,
            alpha_f,
            beta,
            gamma,
        }
    }
    /// Advance the state by one time step.
    ///
    /// The generalized-alpha method evaluates the equation of motion
    /// at interpolated time levels:
    /// ```text
    /// M*a_{n+1-alpha_m} + K*u_{n+1-alpha_f} = F_{n+1-alpha_f}
    /// ```
    pub fn step(
        &mut self,
        mass: &CsrMatrix,
        stiffness: &CsrMatrix,
        force: &[f64],
        dt: f64,
        fixed_dofs: &[usize],
    ) {
        let n = self.u.len();
        let beta = self.beta;
        let gamma = self.gamma;
        let alpha_m = self.alpha_m;
        let alpha_f = self.alpha_f;
        let c0 = 1.0 / (beta * dt * dt);
        let _c1 = gamma / (beta * dt);
        let u_pred: Vec<f64> = (0..n)
            .map(|i| self.u[i] + dt * self.v[i] + dt * dt * (0.5 - beta) * self.a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n)
            .map(|i| self.v[i] + dt * (1.0 - gamma) * self.a[i])
            .collect();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        for row in 0..stiffness.nrows {
            let start = stiffness.row_ptr[row];
            let end = stiffness.row_ptr[row + 1];
            for idx in start..end {
                let val = (1.0 - alpha_f) * stiffness.values[idx];
                triplets.push((row, stiffness.col_indices[idx], val));
            }
        }
        for row in 0..mass.nrows {
            let start = mass.row_ptr[row];
            let end = mass.row_ptr[row + 1];
            for idx in start..end {
                let val = (1.0 - alpha_m) * c0 * mass.values[idx];
                triplets.push((row, mass.col_indices[idx], val));
            }
        }
        let mut k_eff = CsrMatrix::from_triplets(n, n, &triplets);
        let k_u_n = stiffness.mul_vec(&self.u);
        let a_tilde: Vec<f64> = (0..n)
            .map(|i| self.u[i] * c0 + self.v[i] / (beta * dt) + self.a[i] * (0.5 / beta - 1.0))
            .collect();
        let m_a_tilde = mass.mul_vec(&a_tilde);
        let m_a_n = mass.mul_vec(&self.a);
        let mut f_eff: Vec<f64> = (0..n)
            .map(|i| {
                (1.0 - alpha_f) * force[i] + alpha_f * (force[i]) + (1.0 - alpha_m) * m_a_tilde[i]
                    - alpha_m * m_a_n[i]
                    - alpha_f * k_u_n[i]
            })
            .collect();
        for &dof in fixed_dofs {
            let start = k_eff.row_ptr[dof];
            let end = k_eff.row_ptr[dof + 1];
            for idx in start..end {
                k_eff.values[idx] = if k_eff.col_indices[idx] == dof {
                    1.0
                } else {
                    0.0
                };
            }
            for row in 0..n {
                if row == dof {
                    continue;
                }
                let rs = k_eff.row_ptr[row];
                let re = k_eff.row_ptr[row + 1];
                for idx in rs..re {
                    if k_eff.col_indices[idx] == dof {
                        k_eff.values[idx] = 0.0;
                    }
                }
            }
            f_eff[dof] = 0.0;
        }
        let u_new = crate::solvers::conjugate_gradient(&k_eff, &f_eff, &u_pred, 10000, 1e-10);
        let a_new: Vec<f64> = (0..n).map(|i| (u_new[i] - u_pred[i]) * c0).collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + dt * gamma * a_new[i]).collect();
        self.u = u_new;
        self.v = v_new;
        self.a = a_new;
        for &dof in fixed_dofs {
            self.u[dof] = 0.0;
            self.v[dof] = 0.0;
            self.a[dof] = 0.0;
        }
    }
}
/// State for the Wilson-θ time integration method.
///
/// The Wilson-θ method evaluates equilibrium at time t + θ Δt (θ ≥ 1.37
/// for unconditional stability) and interpolates back.  It is first-order
/// accurate in damping but preserves low-frequency content well.
#[allow(dead_code)]
pub struct WilsonThetaState {
    /// Number of DOF.
    pub n_dof: usize,
    /// Current displacements.
    pub u: Vec<f64>,
    /// Current velocities.
    pub v: Vec<f64>,
    /// Current accelerations.
    pub a: Vec<f64>,
    /// θ parameter (default 1.4).
    pub theta: f64,
}
impl WilsonThetaState {
    /// Create a new Wilson-θ integrator (default θ = 1.4).
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
            theta: 1.4,
        }
    }
    /// Single step for diagonal (lumped) systems.
    ///
    /// Solves the modified equation at t + θ Δt then interpolates.
    pub fn step_diagonal(
        &mut self,
        m_diag: &[f64],
        k_diag: &[f64],
        c_diag: &[f64],
        f_n: &[f64],
        f_np1: &[f64],
        dt: f64,
    ) {
        let th = self.theta;
        let n = self.n_dof;
        let f_th: Vec<f64> = (0..n).map(|i| f_n[i] + th * (f_np1[i] - f_n[i])).collect();
        let a0 = 6.0 / (th * th * dt * dt);
        let a1 = 3.0 / (th * dt);
        let mut a_th = vec![0.0_f64; n];
        for i in 0..n {
            let k_eff = k_diag[i] + a0 * m_diag[i] + a1 * c_diag[i];
            let f_eff = f_th[i]
                + m_diag[i] * (a0 * self.u[i] + 3.0 * a1 * self.v[i] + 2.0 * self.a[i])
                + c_diag[i]
                    * (3.0 * self.u[i] / (th * dt) + 2.0 * self.v[i] + th * dt / 2.0 * self.a[i]);
            let _ = f_eff;
            let f_eff2 = f_th[i]
                + m_diag[i] * (a0 * self.u[i] + 3.0 * a1 * self.v[i] + 2.0 * self.a[i])
                + c_diag[i] * (a1 * self.u[i] + 2.0 * self.v[i] + th * dt / 2.0 * self.a[i]);
            let u_th = if k_eff.abs() < 1e-60 {
                self.u[i]
            } else {
                f_eff2 / k_eff
            };
            a_th[i] = (u_th - self.u[i]) * a0 - 3.0 * a1 * self.v[i] - 2.0 * self.a[i];
        }
        for i in 0..n {
            let a_new = self.a[i] + (a_th[i] - self.a[i]) / th;
            let v_new = self.v[i] + dt / 2.0 * (a_new + self.a[i]);
            let u_new = self.u[i] + dt * self.v[i] + dt * dt / 6.0 * (2.0 * self.a[i] + a_new);
            self.a[i] = a_new;
            self.v[i] = v_new;
            self.u[i] = u_new;
        }
    }
}
/// State for seismic analysis including relative displacements.
///
/// In seismic analysis, the equation of motion is expressed in relative coordinates:
/// `u = u_abs - u_g`
/// where `u_g` is the ground displacement.
#[allow(dead_code)]
pub struct SeismicState {
    /// Relative displacements (m).
    pub u_rel: Vec<f64>,
    /// Relative velocities (m/s).
    pub v_rel: Vec<f64>,
    /// Relative accelerations (m/s^2).
    pub a_rel: Vec<f64>,
    /// Ground displacement history.
    pub u_ground: Vec<f64>,
}
#[allow(dead_code)]
impl SeismicState {
    /// Create a new seismic state with zero initial conditions.
    pub fn new(n_dof: usize) -> Self {
        Self {
            u_rel: vec![0.0; n_dof],
            v_rel: vec![0.0; n_dof],
            a_rel: vec![0.0; n_dof],
            u_ground: vec![0.0; n_dof],
        }
    }
    /// Compute the total (absolute) displacement.
    pub fn absolute_displacement(&self) -> Vec<f64> {
        self.u_rel
            .iter()
            .zip(self.u_ground.iter())
            .map(|(ur, ug)| ur + ug)
            .collect()
    }
}
/// Parameters for the linear acceleration Newmark variant (β = 1/6, γ = 1/2).
///
/// This method is conditionally stable with:
/// Ω_crit = sqrt(3) ≈ 1.732
/// (where Ω = ω * dt)
#[allow(dead_code)]
pub struct NewmarkLinearAccelerationParams;
#[allow(dead_code)]
impl NewmarkLinearAccelerationParams {
    /// Beta parameter for linear acceleration.
    pub const BETA: f64 = 1.0 / 6.0;
    /// Gamma parameter.
    pub const GAMMA: f64 = 0.5;
    /// Critical time step: dt_crit = sqrt(3) / omega_max.
    pub fn critical_dt(omega_max: f64) -> f64 {
        if omega_max < f64::EPSILON {
            return f64::MAX;
        }
        3.0_f64.sqrt() / omega_max
    }
}
/// State for the HHT-alpha (Hilber-Hughes-Taylor) time integrator.
///
/// The HHT-alpha method introduces numerical dissipation via a parameter
/// `alpha_hht` in `[-1/3, 0]`. Setting `alpha_hht = 0` recovers the
/// standard Newmark method.
///
/// The method modifies the equation of motion to:
/// ```text
/// M*a_{n+1} + (1+alpha)*C*v_{n+1} - alpha*C*v_n
///   + (1+alpha)*K*u_{n+1} - alpha*K*u_n = (1+alpha)*F_{n+1} - alpha*F_n
/// ```
#[allow(dead_code)]
pub struct HhtAlphaState {
    /// Current displacements.
    pub u: Vec<f64>,
    /// Current velocities.
    pub v: Vec<f64>,
    /// Current accelerations.
    pub a: Vec<f64>,
    /// HHT alpha parameter, typically in `[-1/3, 0]`.
    pub alpha_hht: f64,
    /// Newmark beta parameter (derived from alpha).
    pub beta: f64,
    /// Newmark gamma parameter (derived from alpha).
    pub gamma: f64,
}
#[allow(dead_code)]
impl HhtAlphaState {
    /// Create a new HHT-alpha state with optimal dissipation parameters.
    ///
    /// `alpha_hht` should be in `[-1/3, 0]`.
    /// The Newmark parameters are chosen for optimal accuracy:
    /// `gamma = 0.5 - alpha_hht`, `beta = (1 - alpha_hht)^2 / 4`.
    pub fn new(n_dof: usize, alpha_hht: f64) -> Self {
        let gamma = 0.5 - alpha_hht;
        let beta = (1.0 - alpha_hht).powi(2) / 4.0;
        Self {
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
            alpha_hht,
            beta,
            gamma,
        }
    }
    /// Advance the state by one time step.
    ///
    /// This is a simplified implementation using the effective stiffness
    /// approach with HHT modification. Without damping, the effective
    /// equation becomes:
    /// ```text
    /// M*a_{n+1} + (1+alpha)*K*u_{n+1} = (1+alpha)*F_{n+1} - alpha*F_n + alpha*K*u_n
    /// ```
    pub fn step(
        &mut self,
        mass: &CsrMatrix,
        stiffness: &CsrMatrix,
        force_new: &[f64],
        force_old: &[f64],
        dt: f64,
        fixed_dofs: &[usize],
    ) {
        let n = self.u.len();
        let alpha = self.alpha_hht;
        let beta = self.beta;
        let gamma = self.gamma;
        let inv_beta_dt2 = 1.0 / (beta * dt * dt);
        let u_pred: Vec<f64> = (0..n)
            .map(|i| self.u[i] + dt * self.v[i] + dt * dt * (0.5 - beta) * self.a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n)
            .map(|i| self.v[i] + dt * (1.0 - gamma) * self.a[i])
            .collect();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        for row in 0..stiffness.nrows {
            let start = stiffness.row_ptr[row];
            let end = stiffness.row_ptr[row + 1];
            for idx in start..end {
                let val = (1.0 + alpha) * stiffness.values[idx];
                triplets.push((row, stiffness.col_indices[idx], val));
            }
        }
        for row in 0..mass.nrows {
            let start = mass.row_ptr[row];
            let end = mass.row_ptr[row + 1];
            for idx in start..end {
                let val = mass.values[idx] * inv_beta_dt2;
                triplets.push((row, mass.col_indices[idx], val));
            }
        }
        let mut k_eff = CsrMatrix::from_triplets(n, n, &triplets);
        let k_u_n = stiffness.mul_vec(&self.u);
        let a_tilde: Vec<f64> = (0..n)
            .map(|i| {
                self.u[i] * inv_beta_dt2 + self.v[i] / (beta * dt) + self.a[i] * (0.5 / beta - 1.0)
            })
            .collect();
        let m_a_tilde = mass.mul_vec(&a_tilde);
        let mut f_eff: Vec<f64> = (0..n)
            .map(|i| {
                (1.0 + alpha) * force_new[i] - alpha * force_old[i]
                    + alpha * k_u_n[i]
                    + m_a_tilde[i]
            })
            .collect();
        for &dof in fixed_dofs {
            let start = k_eff.row_ptr[dof];
            let end = k_eff.row_ptr[dof + 1];
            for idx in start..end {
                k_eff.values[idx] = if k_eff.col_indices[idx] == dof {
                    1.0
                } else {
                    0.0
                };
            }
            for row in 0..n {
                if row == dof {
                    continue;
                }
                let rs = k_eff.row_ptr[row];
                let re = k_eff.row_ptr[row + 1];
                for idx in rs..re {
                    if k_eff.col_indices[idx] == dof {
                        k_eff.values[idx] = 0.0;
                    }
                }
            }
            f_eff[dof] = 0.0;
        }
        let u_new = crate::solvers::conjugate_gradient(&k_eff, &f_eff, &u_pred, 10000, 1e-10);
        let a_new: Vec<f64> = (0..n)
            .map(|i| (u_new[i] - u_pred[i]) * inv_beta_dt2)
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + dt * gamma * a_new[i]).collect();
        self.u = u_new;
        self.v = v_new;
        self.a = a_new;
        for &dof in fixed_dofs {
            self.u[dof] = 0.0;
            self.v[dof] = 0.0;
            self.a[dof] = 0.0;
        }
    }
}
/// State for the Newmark-beta implicit time integrator.
///
/// Implements the average acceleration method (beta=0.25, gamma=0.5),
/// which is unconditionally stable for linear systems without damping.
pub struct NewmarkState {
    /// Current displacements (length = n_dof).
    pub u: Vec<f64>,
    /// Current velocities (length = n_dof).
    pub v: Vec<f64>,
    /// Current accelerations (length = n_dof).
    pub a: Vec<f64>,
    /// Newmark beta parameter (0.25 for average acceleration).
    pub beta: f64,
    /// Newmark gamma parameter (0.5 for average acceleration).
    pub gamma: f64,
}
impl NewmarkState {
    /// Create a new `NewmarkState` with zero initial conditions.
    pub fn new(n_dof: usize) -> Self {
        Self {
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
            beta: 0.25,
            gamma: 0.5,
        }
    }
    /// Advance the state by one time step using the Newmark-beta method.
    pub fn step(
        &mut self,
        mass: &CsrMatrix,
        stiffness: &CsrMatrix,
        force: &[f64],
        dt: f64,
        fixed_dofs: &[usize],
    ) {
        let n = self.u.len();
        assert_eq!(mass.nrows, n);
        assert_eq!(stiffness.nrows, n);
        assert_eq!(force.len(), n);
        let beta = self.beta;
        let gamma = self.gamma;
        let inv_beta_dt2 = 1.0 / (beta * dt * dt);
        let inv_beta_dt = 1.0 / (beta * dt);
        let one_over_2beta_minus1 = 0.5 / beta - 1.0;
        let u_pred: Vec<f64> = (0..n)
            .map(|i| self.u[i] + dt * self.v[i] + dt * dt * (0.5 - beta) * self.a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n)
            .map(|i| self.v[i] + dt * (1.0 - gamma) * self.a[i])
            .collect();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        for row in 0..stiffness.nrows {
            let start = stiffness.row_ptr[row];
            let end = stiffness.row_ptr[row + 1];
            for idx in start..end {
                triplets.push((row, stiffness.col_indices[idx], stiffness.values[idx]));
            }
        }
        for row in 0..mass.nrows {
            let start = mass.row_ptr[row];
            let end = mass.row_ptr[row + 1];
            for idx in start..end {
                let val = mass.values[idx] * inv_beta_dt2;
                triplets.push((row, mass.col_indices[idx], val));
            }
        }
        let mut k_eff = CsrMatrix::from_triplets(n, n, &triplets);
        let a_tilde: Vec<f64> = (0..n)
            .map(|i| {
                self.u[i] * inv_beta_dt2
                    + self.v[i] * inv_beta_dt
                    + self.a[i] * one_over_2beta_minus1
            })
            .collect();
        let m_a_tilde = mass.mul_vec(&a_tilde);
        let mut f_eff: Vec<f64> = (0..n).map(|i| force[i] + m_a_tilde[i]).collect();
        for &dof in fixed_dofs {
            assert!(dof < n, "fixed DOF {dof} out of range");
            let start = k_eff.row_ptr[dof];
            let end = k_eff.row_ptr[dof + 1];
            for idx in start..end {
                if k_eff.col_indices[idx] == dof {
                    k_eff.values[idx] = 1.0;
                } else {
                    k_eff.values[idx] = 0.0;
                }
            }
            for row in 0..n {
                if row == dof {
                    continue;
                }
                let rs = k_eff.row_ptr[row];
                let re = k_eff.row_ptr[row + 1];
                for idx in rs..re {
                    if k_eff.col_indices[idx] == dof {
                        k_eff.values[idx] = 0.0;
                    }
                }
            }
            f_eff[dof] = 0.0;
        }
        let u_new = crate::solvers::conjugate_gradient(&k_eff, &f_eff, &u_pred, 10000, 1e-10);
        let a_new: Vec<f64> = (0..n)
            .map(|i| (u_new[i] - u_pred[i]) * inv_beta_dt2)
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + dt * gamma * a_new[i]).collect();
        self.u = u_new;
        self.v = v_new;
        self.a = a_new;
        for &dof in fixed_dofs {
            self.u[dof] = 0.0;
            self.v[dof] = 0.0;
            self.a[dof] = 0.0;
        }
    }
}
/// Input for a multi-mode response spectrum analysis.
#[allow(dead_code)]
pub struct ResponseSpectrumInput {
    /// Natural frequencies \[rad/s\].
    pub omega_n: Vec<f64>,
    /// Modal damping ratios.
    pub zeta: Vec<f64>,
    /// Mode shapes: `phi[mode][dof]`.
    pub phi: Vec<Vec<f64>>,
    /// Modal participation factors (one per mode).
    pub gamma: Vec<f64>,
    /// Spectral displacement values Sd(omega, zeta) \[m\].
    pub sd: Vec<f64>,
    /// Spectral acceleration values Sa(omega, zeta) \[m/s²\].
    pub sa: Vec<f64>,
}
/// Result of a modal superposition analysis.
#[allow(dead_code)]
pub struct ModalSuperpositionResult {
    /// Time history of displacement for each DOF: `disp[time_step][dof]`.
    pub disp: Vec<Vec<f64>>,
    /// Time history of modal coordinates: `q[time_step][mode]`.
    pub modal_coords: Vec<Vec<f64>>,
}
/// Lightweight dynamic solver that wraps state vectors for time integration.
#[allow(dead_code)]
pub struct DynamicSolver {
    /// Number of DOF.
    pub n_dof: usize,
    /// Displacement vector.
    pub u: Vec<f64>,
    /// Velocity vector.
    pub v: Vec<f64>,
    /// Acceleration vector.
    pub a: Vec<f64>,
}
#[allow(dead_code)]
impl DynamicSolver {
    /// Create a new solver with zero initial conditions.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            u: vec![0.0; n_dof],
            v: vec![0.0; n_dof],
            a: vec![0.0; n_dof],
        }
    }
    /// Perform one generalized-alpha time step.
    ///
    /// The generalized-alpha method (Chung & Hulbert, 1993) achieves high-frequency
    /// numerical dissipation while maintaining second-order accuracy.
    ///
    /// # Parameters
    /// * `rho_inf` – spectral radius at infinity ∈ \[0, 1\].
    ///   - `rho_inf = 1.0`: no dissipation (Newmark average acceleration).
    ///   - `rho_inf = 0.0`: maximum dissipation.
    #[allow(clippy::too_many_arguments)]
    pub fn generalized_alpha(
        &mut self,
        mass: &CsrMatrix,
        stiffness: &CsrMatrix,
        force: &[f64],
        dt: f64,
        rho_inf: f64,
    ) {
        let n = self.n_dof;
        let rho = rho_inf.clamp(0.0, 1.0);
        let alpha_m = (2.0 * rho - 1.0) / (rho + 1.0);
        let alpha_f = rho / (rho + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f).powi(2);
        let c0 = 1.0 / (beta * dt * dt);
        let c1 = gamma / (beta * dt);
        let c2 = 1.0 / (beta * dt);
        let c3 = 1.0 / (2.0 * beta) - 1.0;
        let u_pred: Vec<f64> = (0..n)
            .map(|i| self.u[i] + dt * self.v[i] + dt * dt * (0.5 - beta) * self.a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n)
            .map(|i| self.v[i] + dt * (1.0 - gamma) * self.a[i])
            .collect();
        let scale_m = (1.0 - alpha_m) * c0;
        let scale_k = 1.0 - alpha_f;
        let mut k_eff_diag = vec![0.0_f64; n];
        for row in 0..mass.nrows {
            for idx in mass.row_ptr[row]..mass.row_ptr[row + 1] {
                if mass.col_indices[idx] == row {
                    k_eff_diag[row] += scale_m * mass.values[idx];
                }
            }
        }
        for row in 0..stiffness.nrows {
            for idx in stiffness.row_ptr[row]..stiffness.row_ptr[row + 1] {
                if stiffness.col_indices[idx] == row {
                    k_eff_diag[row] += scale_k * stiffness.values[idx];
                }
            }
        }
        let mut f_eff = force.to_vec();
        for row in 0..mass.nrows {
            for idx in mass.row_ptr[row]..mass.row_ptr[row + 1] {
                let col = mass.col_indices[idx];
                let m_val = mass.values[idx];
                f_eff[row] += m_val * (c0 * u_pred[col] + c2 * v_pred[col] + c3 * self.a[col]);
            }
        }
        let mut u_new = vec![0.0_f64; n];
        for i in 0..n {
            u_new[i] = if k_eff_diag[i].abs() > 1e-60 {
                f_eff[i] / k_eff_diag[i]
            } else {
                u_pred[i]
            };
        }
        let a_new: Vec<f64> = (0..n)
            .map(|i| c0 * (u_new[i] - u_pred[i]) - c2 * v_pred[i] - c3 * self.a[i])
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + dt * gamma * a_new[i]).collect();
        let _ = c1;
        self.u = u_new;
        self.v = v_new;
        self.a = a_new;
    }
    /// Compute Rayleigh damping matrix C = alpha * M + beta * K using
    /// two-mode targeting.
    ///
    /// Returns the `CsrMatrix` damping matrix.
    pub fn compute_rayleigh_damping(
        &self,
        mass: &CsrMatrix,
        stiffness: &CsrMatrix,
        omega1: f64,
        zeta1: f64,
        omega2: f64,
        zeta2: f64,
    ) -> CsrMatrix {
        let (alpha, beta) = compute_rayleigh_coefficients(omega1, zeta1, omega2, zeta2);
        rayleigh_damping_matrix(mass, stiffness, alpha, beta)
    }
    /// Apply modal truncation correction (static correction for residual modes).
    ///
    /// See `modal_truncation_correction`.
    pub fn modal_truncation_correction(
        &self,
        omega_n: &[f64],
        phi: &[Vec<f64>],
        force: &[f64],
        n_modes_included: usize,
    ) -> Vec<f64> {
        modal_truncation_correction(omega_n, phi, force, n_modes_included)
    }
}
/// Monitor energy conservation over a sequence of time steps.
///
/// Returns a vector of `(time, kinetic_energy, strain_energy, total_energy)`.
#[allow(dead_code)]
pub struct EnergyMonitor {
    /// History: `(time, T, U, T+U)`.
    pub history: Vec<(f64, f64, f64, f64)>,
}
#[allow(dead_code)]
impl EnergyMonitor {
    /// Create a new energy monitor.
    pub fn new() -> Self {
        Self {
            history: Vec::new(),
        }
    }
    /// Record an energy sample.
    pub fn record(&mut self, time: f64, kinetic: f64, strain: f64) {
        self.history.push((time, kinetic, strain, kinetic + strain));
    }
    /// Compute the maximum relative energy drift: `max(|E_i - E_0|) / |E_0|`.
    pub fn max_relative_drift(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        let e0 = self.history[0].3;
        if e0.abs() < 1e-60 {
            return 0.0;
        }
        self.history
            .iter()
            .map(|&(_, _, _, e)| ((e - e0) / e0).abs())
            .fold(0.0_f64, f64::max)
    }
    /// Get the number of recorded samples.
    pub fn len(&self) -> usize {
        self.history.len()
    }
    /// Check if the monitor has no recorded samples.
    pub fn is_empty(&self) -> bool {
        self.history.is_empty()
    }
}

impl Default for EnergyMonitor {
    fn default() -> Self {
        Self::new()
    }
}
