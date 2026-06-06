// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM integration helpers for non-Newtonian fluid models.
//!
//! Includes the `NonNewtonianLBM` helper, viscosity fields, local tau lattice,
//! thixotropic models, Papanastasiou viscoplastic, turbulence helpers,
//! collision routines, and utility functions.

use super::{CS2, LocalViscosityModel, MU_MAX, MU_MIN, NonNewtonianFluid, NonNewtonianModel};

// ---------------------------------------------------------------------------
// NonNewtonianLBM
// ---------------------------------------------------------------------------

/// Helper for computing an effective relaxation frequency from a
/// non-Newtonian rheology model in LBM simulations.
#[derive(Debug, Clone, Copy)]
pub struct NonNewtonianLBM {
    /// Base relaxation frequency (Newtonian reference).
    pub omega_base: f64,
}

impl NonNewtonianLBM {
    /// Create a new helper with the given base omega.
    pub fn new(omega_base: f64) -> Self {
        Self { omega_base }
    }

    /// Compute the effective omega from the local shear rate and a fluid model.
    ///
    /// Uses `tau = 0.5 + nu_eff / cs^2` where `nu_eff` comes from the model.
    pub fn effective_omega(&self, shear_rate: f64, model: &dyn NonNewtonianFluid) -> f64 {
        let nu_eff = model.effective_viscosity(shear_rate);
        let tau = 0.5 + nu_eff / CS2;
        1.0 / tau
    }
}

// ---------------------------------------------------------------------------
// Apparent viscosity iteration (Picard)
// ---------------------------------------------------------------------------

/// Iterate the apparent viscosity using a Picard (fixed-point) method.
///
/// Given the current shear rate and a non-Newtonian model, iterate:
///
/// ```text
/// mu_{n+1} = (1 - alpha) * mu_n + alpha * model.viscosity(gamma)
/// ```
///
/// Returns the converged viscosity. This is useful for implicit non-Newtonian
/// LBM where the shear rate depends on the viscosity (and vice versa).
///
/// # Arguments
/// * `initial_mu` -- initial guess for the viscosity
/// * `shear_rate` -- scalar shear rate
/// * `model` -- non-Newtonian fluid model
/// * `alpha` -- under-relaxation factor (0 < alpha <= 1)
/// * `max_iter` -- maximum number of iterations
/// * `tol` -- convergence tolerance
pub fn iterate_apparent_viscosity(
    initial_mu: f64,
    shear_rate: f64,
    model: &dyn NonNewtonianFluid,
    alpha: f64,
    max_iter: usize,
    tol: f64,
) -> f64 {
    let mut mu = initial_mu;
    for _ in 0..max_iter {
        let mu_new_raw = model.effective_viscosity(shear_rate);
        let mu_new = (1.0 - alpha) * mu + alpha * mu_new_raw;
        if (mu_new - mu).abs() < tol {
            return mu_new;
        }
        mu = mu_new;
    }
    mu
}

// ---------------------------------------------------------------------------
// ViscosityField
// ---------------------------------------------------------------------------

/// A field storing per-node effective viscosity values.
///
/// Useful for post-processing and visualization of non-Newtonian flows.
pub struct ViscosityField {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Per-node viscosity values.
    pub values: Vec<f64>,
}

impl ViscosityField {
    /// Create a new viscosity field with uniform initial viscosity.
    pub fn new(n_nodes: usize, initial: f64) -> Self {
        Self {
            n_nodes,
            values: vec![initial; n_nodes],
        }
    }

    /// Update from shear rates and a non-Newtonian model.
    pub fn update(&mut self, shear_rates: &[f64], model: &dyn NonNewtonianFluid) {
        assert_eq!(shear_rates.len(), self.n_nodes);
        for (k, &gamma) in shear_rates.iter().enumerate() {
            self.values[k] = model.effective_viscosity(gamma);
        }
    }

    /// Return the minimum viscosity in the field.
    pub fn min_viscosity(&self) -> f64 {
        self.values.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Return the maximum viscosity in the field.
    pub fn max_viscosity(&self) -> f64 {
        self.values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Return the mean viscosity.
    pub fn mean_viscosity(&self) -> f64 {
        if self.n_nodes == 0 {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.n_nodes as f64
    }
}

// ---------------------------------------------------------------------------
// LocalTauLattice
// ---------------------------------------------------------------------------

/// Helper for maintaining a spatially-varying tau field on a lattice.
///
/// Each node carries its own relaxation time derived from the local shear
/// rate, enabling non-Newtonian LBM simulations.
pub struct LocalTauLattice {
    /// Number of lattice nodes.
    pub n_nodes: usize,
    /// Per-node relaxation time.
    pub tau: Vec<f64>,
}

impl LocalTauLattice {
    /// Create a new `LocalTauLattice` with all tau initialised to `tau_init`.
    pub fn new(n_nodes: usize, tau_init: f64) -> Self {
        Self {
            n_nodes,
            tau: vec![tau_init; n_nodes],
        }
    }

    /// Update the tau field from per-node shear rates using the supplied fluid
    /// model.
    ///
    /// # Arguments
    /// * `shear_rates` -- scalar shear rate at each node (length = `n_nodes`)
    /// * `fluid`       -- any `NonNewtonianFluid` implementation
    pub fn update_tau_field(
        &mut self,
        shear_rates: &[f64],
        fluid: &dyn NonNewtonianFluid,
    ) -> &[f64] {
        assert_eq!(
            shear_rates.len(),
            self.n_nodes,
            "shear_rates length must equal n_nodes"
        );
        for (k, &gamma) in shear_rates.iter().enumerate() {
            self.tau[k] = fluid.local_tau(gamma);
        }
        &self.tau
    }

    /// Compute the scalar shear rate from the second invariant of the
    /// strain-rate tensor.
    ///
    /// Given the non-equilibrium part of the stress tensor `Q_neq`
    /// (a 3x3 symmetric matrix), the strain-rate magnitude is:
    ///
    /// ```text
    /// |S| = sqrt( 2 * sum_{a,b} S_{ab}^2 )
    /// ```
    pub fn shear_rate_from_strain_rate_tensor(s: [[f64; 3]; 3]) -> f64 {
        let mut sum_sq = 0.0;
        for row in &s {
            for &val in row {
                sum_sq += val * val;
            }
        }
        (2.0 * sum_sq).sqrt()
    }
}

// ---------------------------------------------------------------------------
// ThixotropicFluid (Mewis-Wagner)
// ---------------------------------------------------------------------------

/// Thixotropic fluid model (simplified Mewis-Wagner / Faltas-Volkov).
///
/// A structural parameter `lambda` (0 = fully broken, 1 = fully structured)
/// evolves over time:
///
/// ```text
/// d(lambda)/dt = a * (1 - lambda) - b * lambda * |gamma_dot|^m
/// ```
///
/// The effective viscosity is:
/// `mu(lambda) = mu_inf + (mu_0 - mu_inf) * lambda`
///
/// where `mu_0` and `mu_inf` are the viscosities of fully-structured and
/// fully-broken microstructure, respectively.
#[derive(Debug, Clone, Copy)]
pub struct ThixotropicFluid {
    /// Rest (fully-structured) viscosity.
    pub mu_0: f64,
    /// High-shear (fully-broken) viscosity.
    pub mu_inf: f64,
    /// Build-up rate coefficient `a` (1/s).
    pub a_buildup: f64,
    /// Break-down rate coefficient `b` (s^(m-1) / s).
    pub b_breakdown: f64,
    /// Breakdown exponent `m` (dimensionless).
    pub m_exp: f64,
    /// Current structural parameter (0..1).
    pub lambda: f64,
}

impl ThixotropicFluid {
    /// Create a new thixotropic fluid at full structure (lambda=1).
    pub fn new(mu_0: f64, mu_inf: f64, a: f64, b: f64, m: f64) -> Self {
        Self {
            mu_0,
            mu_inf,
            a_buildup: a,
            b_breakdown: b,
            m_exp: m,
            lambda: 1.0,
        }
    }

    /// Current effective viscosity.
    pub fn viscosity(&self) -> f64 {
        self.mu_inf + (self.mu_0 - self.mu_inf) * self.lambda
    }

    /// Advance the structural parameter by one Euler step.
    pub fn step(&mut self, shear_rate: f64, dt: f64) {
        let gamma = shear_rate.abs();
        let dl = self.a_buildup * (1.0 - self.lambda)
            - self.b_breakdown * self.lambda * gamma.powf(self.m_exp);
        self.lambda = (self.lambda + dt * dl).clamp(0.0, 1.0);
    }

    /// Equilibrium structural parameter at steady-state shear rate.
    ///
    /// At steady state: `a * (1 - lambda_eq) = b * lambda_eq * gamma^m`
    /// -> `lambda_eq = a / (a + b * gamma^m)`
    pub fn equilibrium_lambda(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let denom = self.a_buildup + self.b_breakdown * gamma.powf(self.m_exp);
        if denom < 1e-30 {
            return 1.0;
        }
        self.a_buildup / denom
    }

    /// Equilibrium viscosity at the given steady-state shear rate.
    pub fn equilibrium_viscosity(&self, shear_rate: f64) -> f64 {
        let lambda_eq = self.equilibrium_lambda(shear_rate);
        self.mu_inf + (self.mu_0 - self.mu_inf) * lambda_eq
    }

    /// Thixotropic time scale (inverse of build-up rate).
    pub fn time_scale(&self) -> f64 {
        if self.a_buildup < 1e-30 {
            f64::INFINITY
        } else {
            1.0 / self.a_buildup
        }
    }
}

// ---------------------------------------------------------------------------
// ThixotropicModel -- structural parameter with LBM integration
// ---------------------------------------------------------------------------

/// Thixotropic fluid model with explicit structural parameter `lambda` and LBM
/// integration helpers.
///
/// The structural parameter evolves as:
///
/// ```text
/// d(lambda)/dt = a * (1 - lambda) - b * lambda * |gamma|^m
/// ```
///
/// The effective viscosity depends on lambda:
///
/// ```text
/// mu_eff = mu_inf + (mu_0 - mu_inf) * lambda^k
/// ```
///
/// where `k` is the structural exponent (typically 1).
#[derive(Debug, Clone, Copy)]
pub struct ThixotropicModel {
    /// Viscosity at full structure (lambda = 1).
    pub mu_0: f64,
    /// Viscosity at broken structure (lambda = 0).
    pub mu_inf: f64,
    /// Build-up rate coefficient `a`.
    pub a: f64,
    /// Breakdown rate coefficient `b`.
    pub b: f64,
    /// Breakdown exponent `m`.
    pub m: f64,
    /// Structural exponent `k` for viscosity coupling.
    pub k_exp: f64,
    /// Current structural parameter lambda in \[0, 1\].
    pub lambda: f64,
}

impl ThixotropicModel {
    /// Create a new thixotropic model at full structure (lambda = 1).
    pub fn new(mu_0: f64, mu_inf: f64, a: f64, b: f64, m: f64) -> Self {
        Self {
            mu_0,
            mu_inf,
            a,
            b,
            m,
            k_exp: 1.0,
            lambda: 1.0,
        }
    }

    /// Create with non-unity structural exponent.
    pub fn with_structural_exponent(
        mu_0: f64,
        mu_inf: f64,
        a: f64,
        b: f64,
        m: f64,
        k_exp: f64,
    ) -> Self {
        Self {
            mu_0,
            mu_inf,
            a,
            b,
            m,
            k_exp,
            lambda: 1.0,
        }
    }

    /// Current effective viscosity: `mu = mu_inf + (mu_0 - mu_inf) * lambda^k`.
    pub fn current_viscosity(&self) -> f64 {
        self.mu_inf + (self.mu_0 - self.mu_inf) * self.lambda.powf(self.k_exp)
    }

    /// Advance the structural parameter by one explicit Euler step.
    pub fn advance(&mut self, shear_rate: f64, dt: f64) {
        let gamma = shear_rate.abs();
        let dl = self.a * (1.0 - self.lambda) - self.b * self.lambda * gamma.powf(self.m);
        self.lambda = (self.lambda + dt * dl).clamp(0.0, 1.0);
    }

    /// Steady-state structural parameter at the given shear rate.
    ///
    /// Obtained by setting d(lambda)/dt = 0:
    /// `lambda_ss = a / (a + b * |gamma|^m)`
    pub fn steady_state_lambda(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let denom = self.a + self.b * gamma.powf(self.m);
        if denom < 1e-30 {
            return 1.0;
        }
        self.a / denom
    }

    /// Steady-state viscosity at the given shear rate.
    pub fn steady_state_viscosity(&self, shear_rate: f64) -> f64 {
        let lss = self.steady_state_lambda(shear_rate);
        self.mu_inf + (self.mu_0 - self.mu_inf) * lss.powf(self.k_exp)
    }

    /// Thixotropic time scale (time to relax to steady state from rest):
    /// `t_thix = 1 / a`.
    pub fn thixotropic_time_scale(&self) -> f64 {
        if self.a < 1e-30 {
            return f64::INFINITY;
        }
        1.0 / self.a
    }

    /// LBM relaxation time from current structural parameter.
    ///
    /// `tau = 0.5 + mu_eff / cs^2`
    pub fn lbm_relaxation_time(&self) -> f64 {
        0.5 + self.current_viscosity() / CS2
    }

    /// Effective LBM relaxation frequency from current structural parameter.
    pub fn lbm_omega(&self) -> f64 {
        1.0 / self.lbm_relaxation_time()
    }

    /// Reset structural parameter to fully structured state.
    pub fn reset(&mut self) {
        self.lambda = 1.0;
    }
}

// ---------------------------------------------------------------------------
// Papanastasiou yield stress fluids (full form)
// ---------------------------------------------------------------------------

/// Extended Papanastasiou regularization: exponential smoothing of yield stress.
///
/// ```text
/// mu(gamma) = mu_inf + (tau_y / |gamma|) * (1 - exp(-m * |gamma|))
///           + K * |gamma|^(n-1)   (optional power-law plastic term)
/// ```
///
/// This combines the Papanastasiou regularization with a power-law viscoplastic
/// contribution, giving a general viscoplastic model.
#[derive(Debug, Clone, Copy)]
pub struct PapanastasiouViscoplastic {
    /// Yield stress tau_y.
    pub tau_y: f64,
    /// Infinite-shear viscosity (Newtonian plateau).
    pub mu_inf: f64,
    /// Power-law consistency K.
    pub consistency_k: f64,
    /// Power-law exponent n.
    pub n: f64,
    /// Regularization parameter m.
    pub m: f64,
}

impl PapanastasiouViscoplastic {
    /// Create with purely Papanastasiou (no power-law term, n=1, K=0).
    pub fn papanastasiou_only(tau_y: f64, mu_inf: f64, m: f64) -> Self {
        Self {
            tau_y,
            mu_inf,
            consistency_k: 0.0,
            n: 1.0,
            m,
        }
    }

    /// Create with full viscoplastic Papanastasiou + power-law.
    pub fn full(tau_y: f64, mu_inf: f64, consistency_k: f64, n: f64, m: f64) -> Self {
        Self {
            tau_y,
            mu_inf,
            consistency_k,
            n,
            m,
        }
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        if gamma < 1e-15 {
            // Limit: tau_y * m + mu_inf (from L'Hopital)
            return (self.tau_y * self.m + self.mu_inf).clamp(MU_MIN, MU_MAX);
        }
        let yield_term = self.tau_y * (1.0 - (-self.m * gamma).exp()) / gamma;
        let plastic_term = if self.consistency_k > 0.0 {
            self.consistency_k * gamma.powf(self.n - 1.0)
        } else {
            0.0
        };
        (self.mu_inf + yield_term + plastic_term).clamp(MU_MIN, MU_MAX)
    }

    /// Shear stress at the given shear rate.
    pub fn stress(&self, shear_rate: f64) -> f64 {
        self.viscosity(shear_rate) * shear_rate.abs()
    }

    /// Return true if the material is approximately yielded.
    ///
    /// The Papanastasiou stress at `gamma` should exceed `tau_y * (1 - 1/e) ~ 0.63 * tau_y`
    /// at the "yield point" `gamma = 1/m`.
    pub fn is_yielded(&self, shear_rate: f64) -> bool {
        shear_rate.abs() > 1.0 / self.m.max(1e-30)
    }
}

// ---------------------------------------------------------------------------
// Non-Newtonian turbulence helpers
// ---------------------------------------------------------------------------

/// Compute the effective viscosity for a turbulent non-Newtonian flow
/// using a mixing length approach.
///
/// In turbulent flow, the total viscosity is:
/// `mu_total = mu_eff_laminar + rho * l_m^2 * |S|`
///
/// where `l_m` is the mixing length and `S` is the strain-rate invariant.
pub fn turbulent_effective_viscosity(
    mu_laminar: f64,
    rho: f64,
    mixing_length: f64,
    strain_rate_invariant: f64,
) -> f64 {
    mu_laminar + rho * mixing_length * mixing_length * strain_rate_invariant
}

/// van Driest mixing-length model for the turbulent viscosity.
///
/// `l_m(y) = kappa * y * (1 - exp(-y_plus / A_plus))`
///
/// where `y_plus = y * u_tau / nu`, `kappa = 0.41` (von Karman constant),
/// `A_plus = 26` (van Driest constant).
pub fn van_driest_mixing_length(
    y_wall_distance: f64,
    u_tau: f64,
    nu: f64,
    kappa: f64,
    a_plus: f64,
) -> f64 {
    if nu < 1e-30 || u_tau < 1e-30 {
        return kappa * y_wall_distance;
    }
    let y_plus = y_wall_distance * u_tau / nu;
    kappa * y_wall_distance * (1.0 - (-y_plus / a_plus).exp())
}

/// Compute the turbulent kinematic viscosity using the Smagorinsky LES model.
///
/// `nu_t = (C_s * delta)^2 * |S|`
///
/// where `C_s` is the Smagorinsky constant (~0.1-0.2), `delta` is the
/// grid spacing, and `|S| = sqrt(2 * S_ij * S_ij)` is the strain-rate magnitude.
pub fn smagorinsky_turbulent_viscosity(
    cs_constant: f64,
    delta: f64,
    strain_rate_magnitude: f64,
) -> f64 {
    let l = cs_constant * delta;
    l * l * strain_rate_magnitude
}

/// Kolmogorov micro-scale for non-Newtonian turbulence.
///
/// `eta_k = (nu_eff^3 / epsilon)^0.25`
///
/// where `epsilon` is the turbulent dissipation rate.
pub fn kolmogorov_scale(nu_eff: f64, epsilon: f64) -> f64 {
    if epsilon < 1e-30 {
        return f64::INFINITY;
    }
    (nu_eff.powi(3) / epsilon).powf(0.25)
}

/// Non-Newtonian flow transition criterion.
///
/// The generalized Reynolds number for power-law fluids:
/// `Re_gen = rho * U^(2-n) * L^n / (K * 8^(n-1))`
///
/// where `U` is velocity, `L` characteristic length, `K` consistency.
pub fn generalized_reynolds_power_law(
    rho: f64,
    velocity: f64,
    length: f64,
    consistency_k: f64,
    n: f64,
) -> f64 {
    if consistency_k < 1e-30 {
        return f64::INFINITY;
    }
    rho * velocity.powf(2.0 - n) * length.powf(n) / (consistency_k * 8_f64.powf(n - 1.0))
}

/// Metzner-Reed generalized Reynolds number for pipe flow:
///
/// `Re_MR = rho * U^2 / (K' * (U/D)^n')`
pub fn metzner_reed_reynolds(
    rho: f64,
    velocity: f64,
    diameter: f64,
    k_prime: f64,
    n_prime: f64,
) -> f64 {
    if k_prime < 1e-30 || diameter < 1e-30 {
        return 0.0;
    }
    rho * velocity.powi(2) / (k_prime * (velocity / diameter).powf(n_prime))
}

// ---------------------------------------------------------------------------
// shear_rate_from_strain_tensor (free function)
// ---------------------------------------------------------------------------

/// Compute the scalar shear rate from a 3x3 symmetric strain-rate tensor S.
///
/// The second invariant is:
///
/// ```text
/// |S| = sqrt(2 * S:S) = sqrt(2 * sum_{i,j} S_{ij}^2)
/// ```
///
/// This is a free-function wrapper complementing
/// `LocalTauLattice::shear_rate_from_strain_rate_tensor`.
pub fn shear_rate_from_strain_tensor(s: [[f64; 3]; 3]) -> f64 {
    let mut sum_sq = 0.0_f64;
    for row in &s {
        for &v in row {
            sum_sq += v * v;
        }
    }
    (2.0 * sum_sq).sqrt()
}

// ---------------------------------------------------------------------------
// update_relaxation_field
// ---------------------------------------------------------------------------

/// Update the per-cell BGK relaxation frequencies from a spatially varying
/// viscosity field.
///
/// For each cell `k` the relaxation frequency is:
///
/// ```text
/// omega_k = 1 / (0.5 + nu_k / cs^2)
/// ```
///
/// where `nu_k = viscosity_field[k]`.
///
/// Returns a `Vec`f64` of per-cell omega values, length `n_cells = f.len()`.
pub fn update_relaxation_field(
    f: &[[f64; 19]],
    viscosity_field: &[f64],
    _dx: f64,
    _dt: f64,
) -> Vec<f64> {
    let n = f.len();
    assert_eq!(
        viscosity_field.len(),
        n,
        "viscosity_field length must match number of cells"
    );
    viscosity_field
        .iter()
        .map(|&nu| {
            let nu_clamped = nu.clamp(MU_MIN, MU_MAX);
            1.0 / (0.5 + nu_clamped / CS2)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// effective_relaxation_time: LBM integration helper
// ---------------------------------------------------------------------------

/// Compute the effective LBM relaxation time from a scalar shear rate and a
/// non-Newtonian fluid model.
///
/// ```text
/// tau_eff(gamma) = 0.5 + mu_eff(gamma) / cs^2
/// ```
///
/// Returns the LBM relaxation time `tau` clamped to `\[0.5 + MU_MIN/CS2, 0.5 + MU_MAX/CS2\]`.
pub fn effective_relaxation_time(shear_rate: f64, model: &dyn NonNewtonianFluid) -> f64 {
    let mu_eff = model.effective_viscosity(shear_rate).clamp(MU_MIN, MU_MAX);
    0.5 + mu_eff / CS2
}

/// Compute the effective LBM relaxation frequency from local shear rate.
///
/// `omega_eff = 1 / tau_eff`
pub fn effective_relaxation_frequency(shear_rate: f64, model: &dyn NonNewtonianFluid) -> f64 {
    1.0 / effective_relaxation_time(shear_rate, model)
}

/// Update a per-node relaxation time array from shear rates using any
/// `NonNewtonianFluid` model.
///
/// For each node `k`:
/// `tau\[k\] = 0.5 + mu_eff(shear_rates\[k\]) / cs^2`
pub fn update_tau_from_shear_rates(shear_rates: &[f64], model: &dyn NonNewtonianFluid) -> Vec<f64> {
    shear_rates
        .iter()
        .map(|&gamma| effective_relaxation_time(gamma, model))
        .collect()
}

// ---------------------------------------------------------------------------
// effective_tau helper
// ---------------------------------------------------------------------------

/// Compute the effective LBM relaxation time from a viscosity value.
///
/// In LBM the relaxation time is linked to kinematic viscosity via:
///
/// ```text
/// tau = 0.5 + nu_eff / cs^2
/// ```
///
/// # Arguments
/// * `viscosity` -- effective kinematic viscosity nu_eff (lattice units)
/// * `cs2` -- speed of sound squared (lattice units; usually 1/3)
/// * `dt` -- time step (lattice units; usually 1.0)
pub fn effective_tau(viscosity: f64, cs2: f64, dt: f64) -> f64 {
    0.5 * dt + viscosity / (cs2 * dt.max(1e-30))
}

// ---------------------------------------------------------------------------
// apply_non_newtonian_collision
// ---------------------------------------------------------------------------

/// D2Q9 weights (local copy for collision helper).
pub(crate) const W9_NN: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 velocity vectors (local copy for collision helper).
pub(crate) const C9_NN: [[i32; 2]; 9] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];

/// Apply a non-Newtonian BGK collision to a single D2Q9 node.
///
/// The local shear rate is estimated from the second invariant of the
/// non-equilibrium stress tensor (Chapman-Enskog approximation).
///
/// # Arguments
/// * `f` -- D2Q9 distribution functions for this node (modified in-place)
/// * `rho` -- macroscopic density
/// * `u` -- macroscopic velocity `\[ux, uy\]`
/// * `model` -- any `NonNewtonianModel` implementation
pub fn apply_non_newtonian_collision(
    f: &mut [f64; 9],
    rho: f64,
    u: [f64; 2],
    model: &dyn NonNewtonianModel,
) {
    // Compute equilibrium distributions
    let ux = u[0];
    let uy = u[1];
    let u_sq = ux * ux + uy * uy;
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cx = C9_NN[i][0] as f64;
        let cy = C9_NN[i][1] as f64;
        let eu = cx * ux + cy * uy;
        feq[i] =
            W9_NN[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
    }

    // Estimate local shear rate from velocity magnitude (simple approximation)
    let shear_rate = u_sq.sqrt();

    // Get effective viscosity and compute omega
    let nu_eff = model.viscosity(shear_rate);
    let tau = effective_tau(nu_eff, CS2, 1.0);
    let omega = 1.0 / tau;

    // BGK relaxation
    for i in 0..9 {
        f[i] -= omega * (f[i] - feq[i]);
    }
}

// ---------------------------------------------------------------------------
// NonNewtonianLbmStep -- combined collision step with local viscosity update
// ---------------------------------------------------------------------------

/// Perform a single BGK collision step on a 2D grid using a spatially varying
/// non-Newtonian viscosity model.
///
/// The relaxation time tau(x,y) is computed from the local shear rate, which is
/// estimated from the non-equilibrium part of the distributions.
pub fn non_newtonian_bgk_step_2d<M: LocalViscosityModel>(
    pop: &mut [[f64; 9]],
    nx: usize,
    ny: usize,
    model: &M,
) {
    use crate::lattice::{CS2, D2Q9_VELOCITIES, equilibrium_d2q9, macros_from_d2q9};
    let n = nx * ny;
    for cell in pop[..n].iter_mut() {
        let f_copy = *cell;
        let (rho, ux, uy) = macros_from_d2q9(&f_copy);
        let feq = equilibrium_d2q9(rho, ux, uy);
        // Estimate shear rate from |fneq|
        let mut pi_xy = 0.0f64;
        let mut pi_xx = 0.0f64;
        for (q, c) in D2Q9_VELOCITIES.iter().enumerate() {
            let cx = c[0] as f64;
            let cy = c[1] as f64;
            let fneq_q = f_copy[q] - feq[q];
            pi_xy += cx * cy * fneq_q;
            pi_xx += cx * cx * fneq_q;
        }
        // Approximate |S| = sqrt(pi_xy^2 + ...) / (2 rho cs^2 tau) <- iterative
        let gamma_approx = (pi_xy * pi_xy + pi_xx * pi_xx).sqrt();
        let mu = model.viscosity(gamma_approx);
        let tau = 0.5 + mu / CS2;
        let omega = 1.0 / tau;
        for (fi, &feq_i) in cell.iter_mut().zip(feq.iter()) {
            *fi += omega * (feq_i - *fi);
        }
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Compute the viscosity index (VI) for a fluid given viscosities at two temperatures.
///
/// Simplified: VI ~ (mu_40 - mu_100) / mu_100 * 100.
pub fn viscosity_index(mu_40: f64, mu_100: f64) -> f64 {
    if mu_100 < 1e-30 {
        return 0.0;
    }
    (mu_40 - mu_100) / mu_100 * 100.0
}

/// Compute the Deborah number De = lambda / t_process.
pub fn deborah_number(lambda: f64, t_process: f64) -> f64 {
    if t_process < 1e-30 {
        return f64::INFINITY;
    }
    lambda / t_process
}

/// Compute the Weissenberg number Wi = lambda * gamma.
pub fn weissenberg_number(lambda: f64, shear_rate: f64) -> f64 {
    lambda * shear_rate.abs()
}

/// Compute the Oldroyd B effective viscosity at steady shear.
///
/// eta_eff = eta_s + eta_p (for UCM/Oldroyd-B, viscosity is shear-rate-independent).
pub fn oldroyd_b_viscosity(eta_s: f64, eta_p: f64) -> f64 {
    eta_s + eta_p
}

/// Compute the Trouton ratio for an extensional flow:
/// Tr = eta_extensional / eta_shear.
///
/// For Newtonian fluids Tr = 3; for polymers Tr can be much larger.
pub fn trouton_ratio(eta_ext: f64, eta_shear: f64) -> f64 {
    if eta_shear < 1e-30 {
        return 0.0;
    }
    eta_ext / eta_shear
}
