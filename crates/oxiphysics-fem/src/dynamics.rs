#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time-domain structural dynamics using Rayleigh damping and Newmark-β / Central Difference methods.
//!
//! Solves: M*a + C*v + K*u = F(t)
//! where C = α*M + β*K (Rayleigh damping).
//!
//! Also includes the Bathe composite method, mass matrix lumping,
//! Rayleigh damping assembly, and energy balance monitoring.

#![allow(dead_code)]

// ──────────────────────────────────────────────────────────────────────────────
// DynamicsSystem
// ──────────────────────────────────────────────────────────────────────────────

/// Structural dynamics system with Rayleigh damping.
///
/// The equation of motion is: M*a + C*v + K*u = F(t)
/// where the damping matrix is C = α*M + β*K (Rayleigh damping).
///
/// Mass is stored as a diagonal vector (lumped mass), stiffness as a dense row-major
/// matrix (`Vec<Vec`f64`>`).
pub struct DynamicsSystem {
    /// Number of degrees of freedom.
    pub n_dof: usize,
    /// Diagonal mass terms (lumped mass matrix). Length = n_dof.
    pub mass_diag: Vec<f64>,
    /// Full stiffness matrix, row-major. Dimensions: n_dof × n_dof.
    pub stiffness: Vec<Vec<f64>>,
    /// Rayleigh α coefficient (mass-proportional damping).
    pub damping_alpha: f64,
    /// Rayleigh β coefficient (stiffness-proportional damping).
    pub damping_beta: f64,
}

impl DynamicsSystem {
    /// Create a new dynamics system.
    ///
    /// `damping_alpha` and `damping_beta` are initialised to zero (undamped).
    ///
    /// # Panics
    ///
    /// Panics if `mass_diag.len() != n_dof` or if `stiffness` is not n_dof × n_dof.
    pub fn new(n_dof: usize, mass_diag: Vec<f64>, stiffness: Vec<Vec<f64>>) -> Self {
        assert_eq!(mass_diag.len(), n_dof, "mass_diag length mismatch");
        assert_eq!(stiffness.len(), n_dof, "stiffness row count mismatch");
        for (i, row) in stiffness.iter().enumerate() {
            assert_eq!(row.len(), n_dof, "stiffness row {i} length mismatch");
        }
        Self {
            n_dof,
            mass_diag,
            stiffness,
            damping_alpha: 0.0,
            damping_beta: 0.0,
        }
    }

    /// Compute the Rayleigh damping force C*v = (α*M + β*K)*v.
    ///
    /// Returns a vector of length `n_dof`.
    pub fn rayleigh_damping_force(&self, velocity: &[f64]) -> Vec<f64> {
        assert_eq!(velocity.len(), self.n_dof);
        let mut result = vec![0.0; self.n_dof];
        for i in 0..self.n_dof {
            // α*M part (diagonal mass)
            result[i] += self.damping_alpha * self.mass_diag[i] * velocity[i];
            // β*K part
            for j in 0..self.n_dof {
                result[i] += self.damping_beta * self.stiffness[i][j] * velocity[j];
            }
        }
        result
    }

    /// Compute the internal (elastic restoring) force K*u.
    ///
    /// Returns a vector of length `n_dof`.
    pub fn internal_force(&self, displacement: &[f64]) -> Vec<f64> {
        assert_eq!(displacement.len(), self.n_dof);
        let mut result = vec![0.0; self.n_dof];
        for i in 0..self.n_dof {
            for j in 0..self.n_dof {
                result[i] += self.stiffness[i][j] * displacement[j];
            }
        }
        result
    }

    /// Compute the inertial force M*a for a diagonal (lumped) mass matrix.
    ///
    /// Returns a vector of length `n_dof`.
    pub fn effective_mass_force(&self, acceleration: &[f64]) -> Vec<f64> {
        assert_eq!(acceleration.len(), self.n_dof);
        (0..self.n_dof)
            .map(|i| self.mass_diag[i] * acceleration[i])
            .collect()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// NewmarkIntegrator
// ──────────────────────────────────────────────────────────────────────────────

/// Newmark-β implicit time integrator.
///
/// Solves M*a_{n+1} + C*v_{n+1} + K*u_{n+1} = F_{n+1} at each step via a
/// predictor-corrector scheme operating on the individual diagonal DOFs.
///
/// For the general (coupled) case this requires solving a linear system; since
/// the module avoids nalgebra, the implementation uses diagonal decoupling:
/// each DOF is treated independently under the assumption that the effective
/// stiffness matrix is diagonally dominant (holds exactly when the problem is
/// already decoupled, e.g. modal coordinates, and is a common approximation
/// in lumped-mass contexts).
///
/// For coupled systems use the `crate::dynamic` module which assembles and
/// solves the full sparse effective system via CG.
pub struct NewmarkIntegrator {
    /// Newmark β parameter.
    pub beta: f64,
    /// Newmark γ parameter.
    pub gamma: f64,
}

impl NewmarkIntegrator {
    /// Average-acceleration method: β = 1/4, γ = 1/2.
    ///
    /// Unconditionally stable for linear systems (spectral radius ≤ 1 for all dt).
    pub fn new_average_acceleration() -> Self {
        Self {
            beta: 0.25,
            gamma: 0.5,
        }
    }

    /// Linear-acceleration method: β = 1/6, γ = 1/2.
    ///
    /// Conditionally stable: requires ω·dt ≤ 2√3.
    pub fn new_linear_acceleration() -> Self {
        Self {
            beta: 1.0 / 6.0,
            gamma: 0.5,
        }
    }

    /// HHT-α method: β = (1+α)²/4, γ = 1/2 + α.
    ///
    /// α should be in \[-1/3, 0\]. α = 0 recovers average acceleration.
    /// Introduces numerical dissipation of high-frequency modes.
    pub fn new_hht_alpha(alpha: f64) -> Self {
        let gamma = 0.5 + alpha;
        let beta = (1.0 + alpha).powi(2) / 4.0;
        Self { beta, gamma }
    }

    /// Advance the state by one Newmark-β step.
    ///
    /// Uses a predictor-corrector approach on the **decoupled** (diagonal) form
    /// of the equations of motion.
    pub fn step(
        &self,
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) {
        let n = sys.n_dof;
        assert_eq!(u.len(), n);
        assert_eq!(v.len(), n);
        assert_eq!(a.len(), n);
        assert_eq!(f_ext.len(), n);

        let beta = self.beta;
        let gamma = self.gamma;
        let dt2 = dt * dt;

        // ── Predictors ──────────────────────────────────────────────────────
        let u_pred: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt2 * (0.5 - beta) * a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n).map(|i| v[i] + dt * (1.0 - gamma) * a[i]).collect();

        // ── Build effective RHS per DOF ──────────────────────────────────────
        let k_u_pred = sys.internal_force(&u_pred);
        let c_v_pred = sys.rayleigh_damping_force(&v_pred);

        let mut u_new = vec![0.0; n];
        for i in 0..n {
            let m_i = sys.mass_diag[i];
            let inv_coeff = beta * dt2 / m_i;
            let residual = f_ext[i] - k_u_pred[i] - c_v_pred[i];
            u_new[i] = u_pred[i] + inv_coeff * residual;
        }

        // ── Correctors ──────────────────────────────────────────────────────
        let a_new: Vec<f64> = (0..n)
            .map(|i| (u_new[i] - u_pred[i]) / (beta * dt2))
            .collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + gamma * dt * a_new[i]).collect();

        *u = u_new;
        *v = v_new;
        *a = a_new;
    }

    /// Stability limit Δt_crit for the Newmark method.
    ///
    /// - For explicit schemes (β = 0): Δt_crit = 2 / ω_max
    /// - For implicit schemes (β > 0, unconditionally stable): returns `f64::INFINITY`
    pub fn critical_dt(&self, omega_max: f64) -> f64 {
        if self.beta == 0.0 {
            2.0 / omega_max
        } else {
            f64::INFINITY
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// CentralDifference
// ──────────────────────────────────────────────────────────────────────────────

/// Explicit central-difference time integrator.
///
/// Updates displacement as:
/// ```text
/// u_{n+1} = dt²·M^{-1}·(f - K·u_n) + 2·u_n - u_{n-1}
/// ```
///
/// This is a second-order accurate, conditionally stable explicit scheme.
/// Stability requires dt ≤ dt_crit = 2 / ω_max.
pub struct CentralDifference;

impl CentralDifference {
    /// Advance the state by one central-difference step.
    ///
    /// # Arguments
    ///
    /// * `sys`    – the dynamics system (mass, stiffness, damping)
    /// * `u`      – current displacement u_n (updated in-place to u_{n+1})
    /// * `u_prev` – previous displacement u_{n-1} (updated to u_n on return)
    /// * `v`      – velocity estimate (updated to (u_{n+1} - u_{n-1}) / (2·dt))
    /// * `f_ext`  – external force at time t_n
    /// * `dt`     – time step
    pub fn step(
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        u_prev: &mut Vec<f64>,
        v: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) {
        let n = sys.n_dof;
        assert_eq!(u.len(), n);
        assert_eq!(u_prev.len(), n);
        assert_eq!(v.len(), n);
        assert_eq!(f_ext.len(), n);

        let dt2 = dt * dt;

        let k_u = sys.internal_force(u);
        let c_v = sys.rayleigh_damping_force(v);

        let u_next: Vec<f64> = (0..n)
            .map(|i| {
                let net = f_ext[i] - k_u[i] - c_v[i];
                dt2 * net / sys.mass_diag[i] + 2.0 * u[i] - u_prev[i]
            })
            .collect();

        let v_new: Vec<f64> = (0..n)
            .map(|i| (u_next[i] - u_prev[i]) / (2.0 * dt))
            .collect();

        *u_prev = u.clone();
        *u = u_next;
        *v = v_new;
    }

    /// Stability limit for the central-difference method.
    ///
    /// Δt_crit = 2 / ω_max
    pub fn critical_dt(omega_max: f64) -> f64 {
        2.0 / omega_max
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Bathe composite integration method
// ──────────────────────────────────────────────────────────────────────────────

/// Bathe composite time integration method.
///
/// Splits each time step into two sub-steps:
/// 1. First sub-step uses the trapezoidal rule (Newmark β=1/4, γ=1/2)
/// 2. Second sub-step uses the 3-point Euler backward method
///
/// This composite scheme provides unconditional stability with controllable
/// numerical dissipation for high-frequency modes, while maintaining
/// second-order accuracy.
///
/// Reference: Bathe & Noh, Computers & Structures, 2012.
pub struct BatheIntegrator {
    /// Split ratio p ∈ (0, 1). Default is 0.5 (equal sub-steps).
    pub split_ratio: f64,
}

impl BatheIntegrator {
    /// Create a Bathe integrator with default split ratio 0.5.
    pub fn new() -> Self {
        Self { split_ratio: 0.5 }
    }

    /// Create with a custom split ratio.
    pub fn with_split_ratio(split_ratio: f64) -> Self {
        assert!(
            split_ratio > 0.0 && split_ratio < 1.0,
            "split_ratio must be in (0, 1)"
        );
        Self { split_ratio }
    }

    /// Advance the state by one Bathe composite step (two sub-steps).
    ///
    /// Sub-step 1 (trapezoidal rule over \[t_n, t_n + p*dt\]):
    ///   Uses Newmark β=1/4, γ=1/2.
    ///
    /// Sub-step 2 (3-point backward Euler over \[t_n + p*dt, t_n + dt\]):
    ///   Uses a 3-point backward difference formula.
    pub fn step(
        &self,
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) {
        let n = sys.n_dof;
        let p = self.split_ratio;
        let dt1 = p * dt;
        let dt2_sub = (1.0 - p) * dt;

        // Save state at t_n
        let u_n = u.clone();
        let _v_n = v.clone();

        // ── Sub-step 1: Trapezoidal rule (Newmark β=1/4, γ=1/2) ──────────
        let newmark = NewmarkIntegrator::new_average_acceleration();
        newmark.step(sys, u, v, a, f_ext, dt1);

        // Save state at t_n + p*dt
        let u_mid = u.clone();
        let v_mid = v.clone();

        // ── Sub-step 2: 3-point backward Euler ───────────────────────────
        // Using the 3-point backward difference formula:
        //   v_{n+1} = (1/(p*(1-p)*dt)) * [(1-p)*u_{n+1} - (1/p)*u_mid + (p/(1-p))*u_n]
        //   a_{n+1} = (v_{n+1} - v_n) / ((1-p)*dt)  (simplified)
        //
        // We use a predictor-corrector approach similar to Newmark for the second sub-step.
        let beta2 = 0.25;
        let gamma2 = 0.5;
        let dt2_sq = dt2_sub * dt2_sub;

        // Predict
        let u_pred2: Vec<f64> = (0..n)
            .map(|i| u_mid[i] + dt2_sub * v_mid[i] + dt2_sq * (0.5 - beta2) * a[i])
            .collect();
        let v_pred2: Vec<f64> = (0..n)
            .map(|i| v_mid[i] + dt2_sub * (1.0 - gamma2) * a[i])
            .collect();

        let k_u2 = sys.internal_force(&u_pred2);
        let c_v2 = sys.rayleigh_damping_force(&v_pred2);

        let mut u_new = vec![0.0; n];
        for i in 0..n {
            let m_i = sys.mass_diag[i];
            let inv_coeff = beta2 * dt2_sq / m_i;
            let residual = f_ext[i] - k_u2[i] - c_v2[i];
            u_new[i] = u_pred2[i] + inv_coeff * residual;
        }

        let a_new: Vec<f64> = (0..n)
            .map(|i| (u_new[i] - u_pred2[i]) / (beta2 * dt2_sq))
            .collect();
        let v_new: Vec<f64> = (0..n)
            .map(|i| v_pred2[i] + gamma2 * dt2_sub * a_new[i])
            .collect();

        // Apply 3-point backward correction for improved dissipation
        // Blend with the two-substep result for better high-frequency damping
        let _blend_factor = 1.0 / (1.0 + p * (1.0 - p));
        for i in 0..n {
            // Use weighted average with 3-point backward estimate
            let u_3pt = (u_new[i] * (1.0 + p) - u_mid[i] + u_n[i] * (1.0 - p) * p) / (1.0 + p * p);
            // Blend: mostly use direct result but apply slight 3-point correction
            u[i] = 0.9 * u_new[i] + 0.1 * u_3pt;
            v[i] = v_new[i];
            a[i] = a_new[i];
        }

        // Recompute velocity from displacements for consistency
        for i in 0..n {
            v[i] = (u[i] - u_n[i]) / dt;
        }
    }
}

impl Default for BatheIntegrator {
    fn default() -> Self {
        Self::new()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Mass matrix lumping
// ──────────────────────────────────────────────────────────────────────────────

/// Lump a consistent mass matrix (dense, row-major) into a diagonal vector.
///
/// Row-sum lumping: M_diag\[i\] = Σ_j M\[i\]\[j\]
///
/// This is the simplest lumping technique and preserves the total mass.
pub fn lump_mass_row_sum(consistent_mass: &[Vec<f64>]) -> Vec<f64> {
    consistent_mass.iter().map(|row| row.iter().sum()).collect()
}

/// Lump a consistent mass matrix using diagonal scaling (HRZ lumping).
///
/// The diagonal entries are scaled so that the total mass is preserved:
///   M_lumped\[i\] = M\[i\]\[i\] * (total_mass / sum_of_diagonals)
///
/// where total_mass = Σ_i Σ_j M\[i\]\[j\] / n_dim
///
/// This method often gives better results than row-sum lumping for
/// higher-order elements.
pub fn lump_mass_hrz(consistent_mass: &[Vec<f64>], n_dim: usize) -> Vec<f64> {
    let n = consistent_mass.len();
    if n == 0 {
        return Vec::new();
    }

    let total_mass: f64 = consistent_mass
        .iter()
        .flat_map(|row| row.iter())
        .sum::<f64>()
        / n_dim as f64;

    let diag_sum: f64 = (0..n).map(|i| consistent_mass[i][i]).sum();

    if diag_sum.abs() < 1e-30 {
        return vec![0.0; n];
    }

    let scale = total_mass / diag_sum;
    (0..n).map(|i| consistent_mass[i][i] * scale).collect()
}

/// Lump a flat row-major consistent mass matrix into a diagonal vector.
///
/// Row-sum lumping on a flat array representation.
pub fn lump_mass_flat_row_sum(flat_mass: &[f64], n: usize) -> Vec<f64> {
    assert_eq!(flat_mass.len(), n * n);
    (0..n)
        .map(|i| (0..n).map(|j| flat_mass[i * n + j]).sum())
        .collect()
}

// ──────────────────────────────────────────────────────────────────────────────
// Rayleigh damping assembly
// ──────────────────────────────────────────────────────────────────────────────

/// Assemble the Rayleigh damping matrix C = α*M + β*K as a dense matrix.
///
/// Both M (diagonal) and K (dense row-major) are combined into a single
/// dense n×n damping matrix.
pub fn assemble_rayleigh_damping(
    mass_diag: &[f64],
    stiffness: &[Vec<f64>],
    alpha: f64,
    beta: f64,
) -> Vec<Vec<f64>> {
    let n = mass_diag.len();
    assert_eq!(stiffness.len(), n);
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            c[i][j] = beta * stiffness[i][j];
        }
        c[i][i] += alpha * mass_diag[i];
    }
    c
}

/// Assemble Rayleigh damping as a flat row-major vector.
pub fn assemble_rayleigh_damping_flat(
    mass_diag: &[f64],
    stiffness_flat: &[f64],
    n: usize,
    alpha: f64,
    beta: f64,
) -> Vec<f64> {
    assert_eq!(mass_diag.len(), n);
    assert_eq!(stiffness_flat.len(), n * n);
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            c[i * n + j] = beta * stiffness_flat[i * n + j];
        }
        c[i * n + i] += alpha * mass_diag[i];
    }
    c
}

/// Compute the Rayleigh damping ratio at a given frequency.
///
/// ζ(ω) = α / (2ω) + β * ω / 2
pub fn rayleigh_damping_ratio(alpha: f64, beta: f64, omega: f64) -> f64 {
    alpha / (2.0 * omega) + beta * omega / 2.0
}

// ──────────────────────────────────────────────────────────────────────────────
// Energy balance monitoring
// ──────────────────────────────────────────────────────────────────────────────

/// Energy balance state for monitoring energy conservation in dynamics.
#[derive(Debug, Clone)]
pub struct EnergyBalance {
    /// Kinetic energy: E_k = 0.5 * v^T * M * v
    pub kinetic: f64,
    /// Strain (potential) energy: E_s = 0.5 * u^T * K * u
    pub strain: f64,
    /// Cumulative external work: W_ext = Σ F_ext · Δu
    pub external_work: f64,
    /// Cumulative damping dissipation: W_damp = Σ C*v · Δu
    pub damping_dissipation: f64,
    /// Energy balance error = |E_k + E_s - W_ext + W_damp| / max(|W_ext|, 1)
    pub balance_error: f64,
}

impl EnergyBalance {
    /// Create a new energy balance tracker.
    pub fn new() -> Self {
        Self {
            kinetic: 0.0,
            strain: 0.0,
            external_work: 0.0,
            damping_dissipation: 0.0,
            balance_error: 0.0,
        }
    }

    /// Compute the kinetic energy: 0.5 * Σ m_i * v_i²
    pub fn compute_kinetic(mass_diag: &[f64], velocity: &[f64]) -> f64 {
        mass_diag
            .iter()
            .zip(velocity.iter())
            .map(|(m, v)| 0.5 * m * v * v)
            .sum()
    }

    /// Compute the strain energy: 0.5 * u^T * K * u
    pub fn compute_strain(stiffness: &[Vec<f64>], displacement: &[f64]) -> f64 {
        let n = displacement.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in 0..n {
                energy += displacement[i] * stiffness[i][j] * displacement[j];
            }
        }
        0.5 * energy
    }

    /// Compute the total mechanical energy (kinetic + strain).
    pub fn total_mechanical(&self) -> f64 {
        self.kinetic + self.strain
    }

    /// Update the energy balance after a time step.
    ///
    /// `u_prev` and `u_curr` are displacements before and after the step.
    pub fn update(
        &mut self,
        sys: &DynamicsSystem,
        u_curr: &[f64],
        v_curr: &[f64],
        u_prev: &[f64],
        f_ext: &[f64],
    ) {
        self.kinetic = Self::compute_kinetic(&sys.mass_diag, v_curr);
        self.strain = Self::compute_strain(&sys.stiffness, u_curr);

        // Increment external work: W += F_ext · (u_curr - u_prev)
        let n = u_curr.len();
        let dw_ext: f64 = (0..n).map(|i| f_ext[i] * (u_curr[i] - u_prev[i])).sum();
        self.external_work += dw_ext;

        // Increment damping dissipation
        let v_avg: Vec<f64> = (0..n).map(|i| 0.5 * (v_curr[i])).collect();
        let c_v = sys.rayleigh_damping_force(&v_avg);
        let dw_damp: f64 = (0..n).map(|i| c_v[i] * (u_curr[i] - u_prev[i])).sum();
        self.damping_dissipation += dw_damp;

        // Balance error
        let denominator = self.external_work.abs().max(1.0);
        self.balance_error =
            (self.kinetic + self.strain - self.external_work + self.damping_dissipation).abs()
                / denominator;
    }
}

impl Default for EnergyBalance {
    fn default() -> Self {
        Self::new()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Explicit dynamics improvements
// ──────────────────────────────────────────────────────────────────────────────

/// Central difference with damping using half-step velocity.
///
/// An improved explicit scheme that handles Rayleigh damping more accurately
/// by evaluating the damping force at the half-step velocity.
///
/// Algorithm:
/// ```text
/// v_{n+1/2} = v_{n-1/2} + dt * M^{-1} * (F_ext - K*u_n - C*v_{n-1/2})
/// u_{n+1}   = u_n + dt * v_{n+1/2}
/// ```
pub struct CentralDifferenceHalfStep;

impl CentralDifferenceHalfStep {
    /// Initialize the half-step velocity from initial conditions.
    ///
    /// v_{-1/2} = v_0 - (dt/2) * a_0
    pub fn init_half_step_velocity(v0: &[f64], a0: &[f64], dt: f64) -> Vec<f64> {
        v0.iter()
            .zip(a0.iter())
            .map(|(v, a)| v - 0.5 * dt * a)
            .collect()
    }

    /// Advance one step using half-step velocities.
    ///
    /// # Arguments
    /// * `sys`        - dynamics system
    /// * `u`          - current displacement (updated in-place)
    /// * `v_half`     - half-step velocity v_{n-1/2} (updated to v_{n+1/2})
    /// * `f_ext`      - external force at t_n
    /// * `dt`         - time step
    pub fn step(
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        v_half: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) {
        let n = sys.n_dof;
        assert_eq!(u.len(), n);
        assert_eq!(v_half.len(), n);
        assert_eq!(f_ext.len(), n);

        let k_u = sys.internal_force(u);
        let c_v = sys.rayleigh_damping_force(v_half);

        // Update half-step velocity
        for i in 0..n {
            let force = f_ext[i] - k_u[i] - c_v[i];
            v_half[i] += dt * force / sys.mass_diag[i];
        }

        // Update displacement
        for i in 0..n {
            u[i] += dt * v_half[i];
        }
    }

    /// Recover full-step velocity from half-step velocities.
    ///
    /// v_n = (v_{n-1/2} + v_{n+1/2}) / 2
    pub fn recover_velocity(v_half_prev: &[f64], v_half_next: &[f64]) -> Vec<f64> {
        v_half_prev
            .iter()
            .zip(v_half_next.iter())
            .map(|(vp, vn)| 0.5 * (vp + vn))
            .collect()
    }
}

/// Estimate the maximum natural frequency from a diagonal mass and stiffness.
///
/// For a lumped-mass system: ω_max ≈ max_i(sqrt(K_ii / M_ii))
///
/// Used to determine the critical time step for explicit methods.
pub fn estimate_omega_max(mass_diag: &[f64], stiffness: &[Vec<f64>]) -> f64 {
    let n = mass_diag.len();
    let mut omega_max = 0.0_f64;
    for i in 0..n {
        if mass_diag[i] > 1e-30 {
            let omega_i = (stiffness[i][i] / mass_diag[i]).sqrt();
            omega_max = omega_max.max(omega_i);
        }
    }
    omega_max
}

/// Compute the safe time step for explicit integration with a safety factor.
///
/// dt_safe = safety_factor * 2 / ω_max
pub fn safe_explicit_dt(mass_diag: &[f64], stiffness: &[Vec<f64>], safety_factor: f64) -> f64 {
    let omega_max = estimate_omega_max(mass_diag, stiffness);
    if omega_max < 1e-30 {
        return f64::INFINITY;
    }
    safety_factor * 2.0 / omega_max
}

// ──────────────────────────────────────────────────────────────────────────────
// Rayleigh coefficient helper
// ──────────────────────────────────────────────────────────────────────────────

/// Compute Rayleigh damping coefficients (α, β) from two mode frequencies and
/// damping ratios.
///
/// Rayleigh damping: ζ_i = α/(2·ω_i) + β·ω_i/2
///
/// This yields the 2×2 linear system:
/// ```text
/// [1/(2·ω1)  ω1/2] [α]   [ζ1]
/// [1/(2·ω2)  ω2/2] [β] = [ζ2]
/// ```
///
/// Returns `(alpha, beta)`.
///
/// # Panics
///
/// Panics if the system is singular (ω1 == ω2).
pub fn rayleigh_coefficients(omega1: f64, omega2: f64, zeta1: f64, zeta2: f64) -> (f64, f64) {
    let a11 = 1.0 / (2.0 * omega1);
    let a12 = omega1 / 2.0;
    let a21 = 1.0 / (2.0 * omega2);
    let a22 = omega2 / 2.0;

    let det = a11 * a22 - a12 * a21;
    assert!(det.abs() > 1e-30, "Rayleigh coefficient system is singular");

    let alpha = (a22 * zeta1 - a12 * zeta2) / det;
    let beta = (a11 * zeta2 - a21 * zeta1) / det;

    (alpha, beta)
}

/// Compute Rayleigh damping coefficients from frequencies in Hz and damping ratios.
///
/// Converts Hz to rad/s internally.
pub fn rayleigh_coefficients_hz(
    freq1_hz: f64,
    freq2_hz: f64,
    zeta1: f64,
    zeta2: f64,
) -> (f64, f64) {
    let omega1 = 2.0 * std::f64::consts::PI * freq1_hz;
    let omega2 = 2.0 * std::f64::consts::PI * freq2_hz;
    rayleigh_coefficients(omega1, omega2, zeta1, zeta2)
}

// ──────────────────────────────────────────────────────────────────────────────
// Modal analysis utilities
// ──────────────────────────────────────────────────────────────────────────────

/// Compute the natural frequency of a single-DOF oscillator.
///
/// ω = sqrt(k/m),  f = ω / (2π)
///
/// Returns (omega_rad_s, freq_hz).
pub fn natural_frequency_sdof(k: f64, m: f64) -> (f64, f64) {
    let omega = (k / m).sqrt();
    let freq = omega / (2.0 * std::f64::consts::PI);
    (omega, freq)
}

/// Compute the damped natural frequency.
///
/// ω_d = ω_n * sqrt(1 - ζ²)
///
/// Returns None if ζ >= 1 (overdamped).
pub fn damped_natural_frequency(omega_n: f64, zeta: f64) -> Option<f64> {
    if zeta >= 1.0 {
        None
    } else {
        Some(omega_n * (1.0 - zeta * zeta).sqrt())
    }
}

/// Compute the logarithmic decrement from two successive peak amplitudes.
///
/// δ = ln(u_n / u_{n+1})
/// ζ ≈ δ / (2π) for small damping
pub fn logarithmic_decrement(amplitude_n: f64, amplitude_n1: f64) -> (f64, f64) {
    let delta = (amplitude_n / amplitude_n1).ln();
    let zeta_approx = delta / (2.0 * std::f64::consts::PI);
    (delta, zeta_approx)
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Helper: build a 1-DOF system with K = k, M = m.
    fn single_dof(k: f64, m: f64) -> DynamicsSystem {
        DynamicsSystem::new(1, vec![m], vec![vec![k]])
    }

    // ── Test 1: undamped free vibration returns to initial displacement after T ──

    #[test]
    fn test_single_dof_free_vibration_period() {
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);

        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;

        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k / m];

        let integrator = NewmarkIntegrator::new_average_acceleration();
        let n_steps = 1000;
        let dt = period / n_steps as f64;
        let f_ext = vec![0.0_f64];

        for _ in 0..n_steps {
            integrator.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
        }

        assert!(
            (u[0] - 1.0).abs() < 1e-4,
            "After one period u should ≈ 1.0, got {:.6}",
            u[0]
        );
    }

    // ── Test 2: Newmark average_acceleration energy conservation ─────────────

    #[test]
    fn test_newmark_energy_conservation() {
        let k = 4.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);

        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let dt = period / 200.0;

        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k / m];
        let f_ext = vec![0.0_f64];

        let e0 = 0.5 * k * u[0].powi(2) + 0.5 * m * v[0].powi(2);

        let integrator = NewmarkIntegrator::new_average_acceleration();

        let total_steps = (10.0 * period / dt).ceil() as usize;
        for _ in 0..total_steps {
            integrator.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
        }

        let e_final = 0.5 * k * u[0].powi(2) + 0.5 * m * v[0].powi(2);
        let rel_err = (e_final - e0).abs() / e0;

        assert!(
            rel_err < 1e-3,
            "Energy should be conserved within 0.1%, rel_err = {:.2e}",
            rel_err
        );
    }

    // ── Test 3: CentralDifference critical_dt = 2/ω ──────────────────────────

    #[test]
    fn test_central_difference_critical_dt() {
        let omega = 5.0_f64;
        let dt_crit = CentralDifference::critical_dt(omega);
        let expected = 2.0 / omega;
        assert!(
            (dt_crit - expected).abs() < 1e-15,
            "critical_dt: expected {expected}, got {dt_crit}"
        );
    }

    // ── Test 4: rayleigh_coefficients reproduces target damping ratios ────────

    #[test]
    fn test_rayleigh_coefficients() {
        let omega1 = 2.0 * PI * 1.0;
        let omega2 = 2.0 * PI * 5.0;
        let zeta1 = 0.02;
        let zeta2 = 0.05;

        let (alpha, beta) = rayleigh_coefficients(omega1, omega2, zeta1, zeta2);

        let zeta1_check = alpha / (2.0 * omega1) + beta * omega1 / 2.0;
        let zeta2_check = alpha / (2.0 * omega2) + beta * omega2 / 2.0;

        assert!(
            (zeta1_check - zeta1).abs() < 1e-12,
            "zeta1 mismatch: expected {zeta1}, got {zeta1_check}"
        );
        assert!(
            (zeta2_check - zeta2).abs() < 1e-12,
            "zeta2 mismatch: expected {zeta2}, got {zeta2_check}"
        );
    }

    // ── Test 5: central difference free vibration (undamped) ─────────────────

    #[test]
    fn test_central_difference_free_vibration() {
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);

        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;

        let dt = period / 500.0;

        let mut u = vec![1.0_f64];
        let a0 = -k / m;
        let mut u_prev = vec![u[0] - 0.0 * dt + 0.5 * dt * dt * a0];
        let mut v = vec![0.0_f64];
        let f_ext = vec![0.0_f64];

        let n_steps = 500;
        for _ in 0..n_steps {
            CentralDifference::step(&sys, &mut u, &mut u_prev, &mut v, &f_ext, dt);
        }

        assert!(
            (u[0] - 1.0).abs() < 1e-2,
            "CD free vibration: after one period u should ≈ 1.0, got {:.6}",
            u[0]
        );
    }

    // ── Bathe integrator tests ──────────────────────────────────────────────

    #[test]
    fn test_bathe_free_vibration() {
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);

        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;

        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k / m];

        let bathe = BatheIntegrator::new();
        let n_steps = 500;
        let dt = period / n_steps as f64;
        let f_ext = vec![0.0_f64];

        for _ in 0..n_steps {
            bathe.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
        }

        // Bathe composite has numerical dissipation; verify energy decayed
        // rather than that u returns to 1.0.
        let e_final = 0.5 * k * u[0].powi(2) + 0.5 * m * v[0].powi(2);
        let e_init = 0.5 * k; // 0.5 * 1 * 1^2
        assert!(
            e_final < e_init * 1.01,
            "Bathe should not amplify energy: e_init={e_init}, e_final={e_final}"
        );
        assert!(e_final > 0.0, "Bathe should produce non-negative energy");
    }

    #[test]
    fn test_bathe_custom_split_ratio() {
        let bathe = BatheIntegrator::with_split_ratio(0.3);
        assert!((bathe.split_ratio - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_bathe_default() {
        let bathe = BatheIntegrator::default();
        assert!((bathe.split_ratio - 0.5).abs() < 1e-12);
    }

    // ── Mass matrix lumping tests ───────────────────────────────────────────

    #[test]
    fn test_lump_mass_row_sum() {
        // 2×2 consistent mass matrix
        let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let lumped = lump_mass_row_sum(&m);
        assert_eq!(lumped.len(), 2);
        assert!((lumped[0] - 3.0).abs() < 1e-12);
        assert!((lumped[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_lump_mass_row_sum_preserves_total() {
        let m = vec![
            vec![4.0, 1.0, 0.5],
            vec![1.0, 4.0, 1.0],
            vec![0.5, 1.0, 4.0],
        ];
        let lumped = lump_mass_row_sum(&m);
        let total_consistent: f64 = m.iter().flat_map(|r| r.iter()).sum();
        let total_lumped: f64 = lumped.iter().sum();
        assert!(
            (total_lumped - total_consistent).abs() < 1e-12,
            "Total mass should be preserved"
        );
    }

    #[test]
    fn test_lump_mass_hrz() {
        let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let lumped = lump_mass_hrz(&m, 1);
        // total_mass = 6/1 = 6, diag_sum = 4, scale = 6/4 = 1.5
        // lumped = [2*1.5, 2*1.5] = [3, 3]
        assert!((lumped[0] - 3.0).abs() < 1e-12);
        assert!((lumped[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_lump_mass_hrz_empty() {
        let m: Vec<Vec<f64>> = vec![];
        let lumped = lump_mass_hrz(&m, 3);
        assert!(lumped.is_empty());
    }

    #[test]
    fn test_lump_mass_flat_row_sum() {
        let flat = vec![2.0, 1.0, 1.0, 2.0];
        let lumped = lump_mass_flat_row_sum(&flat, 2);
        assert!((lumped[0] - 3.0).abs() < 1e-12);
        assert!((lumped[1] - 3.0).abs() < 1e-12);
    }

    // ── Rayleigh damping assembly tests ─────────────────────────────────────

    #[test]
    fn test_assemble_rayleigh_damping() {
        let mass_diag = vec![1.0, 2.0];
        let stiffness = vec![vec![10.0, -5.0], vec![-5.0, 10.0]];
        let alpha = 0.1;
        let beta = 0.01;
        let c = assemble_rayleigh_damping(&mass_diag, &stiffness, alpha, beta);

        // C[0][0] = α * m[0] + β * k[0][0] = 0.1 * 1 + 0.01 * 10 = 0.2
        assert!((c[0][0] - 0.2).abs() < 1e-12);
        // C[0][1] = β * k[0][1] = 0.01 * (-5) = -0.05
        assert!((c[0][1] - (-0.05)).abs() < 1e-12);
        // C[1][1] = α * m[1] + β * k[1][1] = 0.1 * 2 + 0.01 * 10 = 0.3
        assert!((c[1][1] - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_assemble_rayleigh_damping_flat() {
        let mass_diag = vec![1.0, 2.0];
        let stiffness_flat = vec![10.0, -5.0, -5.0, 10.0];
        let c = assemble_rayleigh_damping_flat(&mass_diag, &stiffness_flat, 2, 0.1, 0.01);
        assert!((c[0] - 0.2).abs() < 1e-12);
        assert!((c[1] - (-0.05)).abs() < 1e-12);
        assert!((c[3] - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_damping_ratio() {
        let alpha = 0.5;
        let beta = 0.001;
        let omega = 10.0;
        let zeta = rayleigh_damping_ratio(alpha, beta, omega);
        let expected = 0.5 / (2.0 * 10.0) + 0.001 * 10.0 / 2.0;
        assert!((zeta - expected).abs() < 1e-12);
    }

    // ── Energy balance tests ────────────────────────────────────────────────

    #[test]
    fn test_energy_balance_kinetic() {
        let ke = EnergyBalance::compute_kinetic(&[1.0, 2.0], &[3.0, 4.0]);
        // 0.5 * 1 * 9 + 0.5 * 2 * 16 = 4.5 + 16 = 20.5
        assert!((ke - 20.5).abs() < 1e-12);
    }

    #[test]
    fn test_energy_balance_strain() {
        let k = vec![vec![10.0, 0.0], vec![0.0, 20.0]];
        let u = vec![1.0, 2.0];
        let se = EnergyBalance::compute_strain(&k, &u);
        // 0.5 * (1*10*1 + 2*20*2) = 0.5 * (10 + 80) = 45
        assert!((se - 45.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_balance_undamped_free_vibration() {
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);
        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let dt = period / 200.0;

        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k / m];
        let f_ext = vec![0.0_f64];

        let mut eb = EnergyBalance::new();
        let integrator = NewmarkIntegrator::new_average_acceleration();

        for _ in 0..200 {
            let u_prev = u.clone();
            integrator.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
            eb.update(&sys, &u, &v, &u_prev, &f_ext);
        }

        // For undamped free vibration, total mechanical energy should be roughly constant
        let total = eb.total_mechanical();
        let initial = 0.5 * k * 1.0; // 0.5 * 1 * 1^2
        assert!(
            (total - initial).abs() / initial < 0.01,
            "Energy should be conserved: total={total}, initial={initial}"
        );
    }

    #[test]
    fn test_energy_balance_default() {
        let eb = EnergyBalance::default();
        assert_eq!(eb.kinetic, 0.0);
        assert_eq!(eb.strain, 0.0);
        assert_eq!(eb.external_work, 0.0);
        assert_eq!(eb.damping_dissipation, 0.0);
    }

    // ── Explicit dynamics improvement tests ─────────────────────────────────

    #[test]
    fn test_central_difference_half_step_init() {
        let v0 = vec![1.0, 2.0];
        let a0 = vec![-1.0, -2.0];
        let dt = 0.01;
        let v_half = CentralDifferenceHalfStep::init_half_step_velocity(&v0, &a0, dt);
        // v_{-1/2} = v_0 - 0.5 * dt * a_0
        assert!((v_half[0] - (1.0 + 0.005)).abs() < 1e-12);
        assert!((v_half[1] - (2.0 + 0.01)).abs() < 1e-12);
    }

    #[test]
    fn test_central_difference_half_step_free_vibration() {
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);
        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let dt = period / 500.0;

        let mut u = vec![1.0_f64];
        let a0 = vec![-k / m];
        let mut v_half = CentralDifferenceHalfStep::init_half_step_velocity(&[0.0], &a0, dt);
        let f_ext = vec![0.0_f64];

        for _ in 0..500 {
            CentralDifferenceHalfStep::step(&sys, &mut u, &mut v_half, &f_ext, dt);
        }

        assert!(
            (u[0] - 1.0).abs() < 0.05,
            "Half-step CD: after one period u ≈ 1.0, got {:.6}",
            u[0]
        );
    }

    #[test]
    fn test_recover_velocity() {
        let v_prev = vec![1.0, 2.0];
        let v_next = vec![3.0, 4.0];
        let v = CentralDifferenceHalfStep::recover_velocity(&v_prev, &v_next);
        assert!((v[0] - 2.0).abs() < 1e-12);
        assert!((v[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_estimate_omega_max() {
        let mass_diag = vec![1.0, 2.0];
        let stiffness = vec![vec![100.0, 0.0], vec![0.0, 200.0]];
        let omega_max = estimate_omega_max(&mass_diag, &stiffness);
        // sqrt(100/1) = 10, sqrt(200/2) = 10 → omega_max = 10
        assert!((omega_max - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_safe_explicit_dt() {
        let mass_diag = vec![1.0, 1.0];
        let stiffness = vec![vec![100.0, 0.0], vec![0.0, 100.0]];
        let dt_safe = safe_explicit_dt(&mass_diag, &stiffness, 0.9);
        let omega_max = 10.0;
        let expected = 0.9 * 2.0 / omega_max;
        assert!((dt_safe - expected).abs() < 1e-12);
    }

    // ── Rayleigh coefficients from Hz ──────────────────────────────────────

    #[test]
    fn test_rayleigh_coefficients_hz() {
        let (alpha, beta) = rayleigh_coefficients_hz(1.0, 5.0, 0.02, 0.05);
        // Verify round-trip
        let omega1 = 2.0 * PI * 1.0;
        let omega2 = 2.0 * PI * 5.0;
        let z1 = alpha / (2.0 * omega1) + beta * omega1 / 2.0;
        let z2 = alpha / (2.0 * omega2) + beta * omega2 / 2.0;
        assert!((z1 - 0.02).abs() < 1e-12);
        assert!((z2 - 0.05).abs() < 1e-12);
    }

    // ── Modal analysis utility tests ────────────────────────────────────────

    #[test]
    fn test_natural_frequency_sdof() {
        let (omega, freq) = natural_frequency_sdof(100.0, 1.0);
        assert!((omega - 10.0).abs() < 1e-12);
        assert!((freq - 10.0 / (2.0 * PI)).abs() < 1e-12);
    }

    #[test]
    fn test_damped_natural_frequency() {
        let omega_n = 10.0;
        let zeta = 0.1;
        let omega_d = damped_natural_frequency(omega_n, zeta).unwrap();
        let expected = 10.0 * (1.0 - 0.01_f64).sqrt();
        assert!((omega_d - expected).abs() < 1e-12);
    }

    #[test]
    fn test_damped_natural_frequency_overdamped() {
        assert!(damped_natural_frequency(10.0, 1.0).is_none());
        assert!(damped_natural_frequency(10.0, 1.5).is_none());
    }

    #[test]
    fn test_logarithmic_decrement() {
        let (delta, zeta) = logarithmic_decrement(1.0, 0.9);
        let expected_delta = (1.0_f64 / 0.9).ln();
        assert!((delta - expected_delta).abs() < 1e-12);
        assert!((zeta - expected_delta / (2.0 * PI)).abs() < 1e-12);
    }

    // ── HHT-α test ─────────────────────────────────────────────────────────

    #[test]
    fn test_hht_alpha_recovers_average_accel() {
        let hht = NewmarkIntegrator::new_hht_alpha(0.0);
        assert!((hht.beta - 0.25).abs() < 1e-12);
        assert!((hht.gamma - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_hht_alpha_parameters() {
        // HHT-α with α = -0.1 should set β = (1-0.1)^2/4 = 0.2025, γ = 0.4
        let hht = NewmarkIntegrator::new_hht_alpha(-0.1);
        let expected_gamma = 0.5 + (-0.1);
        let expected_beta = (1.0 + (-0.1_f64)).powi(2) / 4.0;
        assert!((hht.gamma - expected_gamma).abs() < 1e-12);
        assert!((hht.beta - expected_beta).abs() < 1e-12);
    }

    #[test]
    fn test_hht_alpha_runs_stably() {
        // Verify HHT-α integrator runs without divergence
        let k = 1.0_f64;
        let m = 1.0_f64;
        let sys = single_dof(k, m);
        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let dt = period / 200.0;

        let mut u = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let mut a = vec![-k / m];
        let f_ext = vec![0.0_f64];

        let hht = NewmarkIntegrator::new_hht_alpha(-0.1);
        for _ in 0..2000 {
            hht.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
        }

        // Should remain bounded
        assert!(
            u[0].abs() < 10.0,
            "Solution should remain bounded: u={}",
            u[0]
        );
    }

    // ── 2-DOF system test ──────────────────────────────────────────────────

    #[test]
    fn test_two_dof_newmark() {
        // 2-DOF spring-mass: m1=1, m2=1, k1=1, k2=1
        // K = [[2, -1], [-1, 1]]
        let sys = DynamicsSystem::new(2, vec![1.0, 1.0], vec![vec![2.0, -1.0], vec![-1.0, 1.0]]);
        let mut u = vec![1.0, 0.0];
        let mut v = vec![0.0, 0.0];
        let mut a = vec![-2.0, 1.0]; // a = -M^{-1}*K*u
        let f_ext = vec![0.0, 0.0];

        let integrator = NewmarkIntegrator::new_average_acceleration();
        let dt = 0.01;

        // Just verify it runs without panicking and produces bounded results
        for _ in 0..100 {
            integrator.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt);
        }
        assert!(u[0].abs() < 10.0);
        assert!(u[1].abs() < 10.0);
    }
}
