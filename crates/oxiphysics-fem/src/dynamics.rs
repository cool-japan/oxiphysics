// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time-domain structural dynamics using Rayleigh damping and Newmark-β / HHT-α / Bathe /
//! Central Difference methods.
//!
//! Solves: M*a + C*v + K*u = F(t)
//! where C = α*M + β*K (Rayleigh damping).
//!
//! The implicit integrators ([`NewmarkIntegrator`], including HHT-α, and [`BatheIntegrator`])
//! assemble and factorise the full, generally non-diagonal effective matrix of each step, so
//! they keep the stability properties of the methods; the explicit [`CentralDifference`] schemes
//! need a lumped mass matrix and a time step below the limit of [`safe_explicit_dt`].
//!
//! Also includes mass matrix lumping, Rayleigh damping assembly, and energy balance monitoring.

use crate::dense_lu::LuFactor;
use crate::error::{Error, Result};

/// Check that `dt` is a usable (positive, finite) time step.
fn check_time_step(dt: f64) -> Result<()> {
    if dt.is_finite() && dt > 0.0 {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!(
            "time step must be positive and finite, got {dt}"
        )))
    }
}

/// Check that a vector has the expected length.
fn check_len(name: &str, got: usize, expected: usize) -> Result<()> {
    if got == expected {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!(
            "{name} has {got} entries, expected {expected}"
        )))
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// DynamicsSystem
// ──────────────────────────────────────────────────────────────────────────────

/// Structural dynamics system with Rayleigh damping.
///
/// The equation of motion is: M*a + C*v + K*u = F(t)
/// where the damping matrix is C = α*M + β*K (Rayleigh damping).
///
/// Mass is stored as a diagonal vector (lumped mass), stiffness as a dense row-major
/// matrix (`Vec<Vec<f64>>`).
#[derive(Debug, Clone)]
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
            for (j, &vel_j) in velocity.iter().enumerate().take(self.n_dof) {
                result[i] += self.damping_beta * self.stiffness[i][j] * vel_j;
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
        for (i, result_i) in result.iter_mut().enumerate().take(self.n_dof) {
            for (j, &disp_j) in displacement.iter().enumerate().take(self.n_dof) {
                *result_i += self.stiffness[i][j] * disp_j;
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

    /// Acceleration consistent with the equation of motion at the given state:
    /// `a = M⁻¹ (f − C·v − K·u)`.
    ///
    /// The Newmark-type integrators ([`NewmarkIntegrator`], [`BatheIntegrator`]) advance
    /// `(u, v, a)`, and their first step is only accurate when the supplied `a` is the
    /// acceleration implied by `u`, `v` and the load at the initial time.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when a vector does not have `n_dof` entries, when the system is
    /// inconsistent, or when a mass is not strictly positive (a massless DOF has no
    /// well-defined acceleration).
    pub fn initial_acceleration(&self, u: &[f64], v: &[f64], f_ext: &[f64]) -> Result<Vec<f64>> {
        self.validate()?;
        check_len("displacement vector", u.len(), self.n_dof)?;
        check_len("velocity vector", v.len(), self.n_dof)?;
        check_len("load vector", f_ext.len(), self.n_dof)?;
        let k_u = self.internal_force(u);
        let c_v = self.rayleigh_damping_force(v);
        self.mass_diag
            .iter()
            .enumerate()
            .map(|(i, &mass)| {
                if mass > 0.0 {
                    Ok((f_ext[i] - k_u[i] - c_v[i]) / mass)
                } else {
                    Err(Error::InvalidInput(format!(
                        "DOF {i} has mass {mass}; the initial acceleration needs a positive mass"
                    )))
                }
            })
            .collect()
    }

    /// Check that the public fields describe a consistent system: matching sizes, finite
    /// non-negative masses, and finite stiffness and damping coefficients.
    fn validate(&self) -> Result<()> {
        let n = self.n_dof;
        check_len("mass_diag", self.mass_diag.len(), n)?;
        check_len("stiffness (row count)", self.stiffness.len(), n)?;
        for (i, row) in self.stiffness.iter().enumerate() {
            check_len(&format!("stiffness row {i}"), row.len(), n)?;
            if row.iter().any(|k| !k.is_finite()) {
                return Err(Error::InvalidInput(format!(
                    "stiffness row {i} contains a non-finite entry"
                )));
            }
        }
        for (i, &mass) in self.mass_diag.iter().enumerate() {
            if !(mass.is_finite() && mass >= 0.0) {
                return Err(Error::InvalidInput(format!(
                    "mass of DOF {i} must be finite and non-negative, got {mass}"
                )));
            }
        }
        if !(self.damping_alpha.is_finite() && self.damping_beta.is_finite()) {
            return Err(Error::InvalidInput(
                "Rayleigh damping coefficients must be finite".into(),
            ));
        }
        Ok(())
    }

    /// The dense matrix `m_coeff·diag(M) + k_coeff·K`, the building block of the effective
    /// matrices of the implicit integrators.
    fn mass_stiffness_combination(&self, m_coeff: f64, k_coeff: f64) -> Vec<Vec<f64>> {
        let mut combined: Vec<Vec<f64>> = self
            .stiffness
            .iter()
            .map(|row| row.iter().map(|k| k_coeff * k).collect())
            .collect();
        for (i, row) in combined.iter_mut().enumerate() {
            row[i] += m_coeff * self.mass_diag[i];
        }
        combined
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// NewmarkIntegrator
// ──────────────────────────────────────────────────────────────────────────────

/// Newmark-β / Hilber–Hughes–Taylor (HHT-α) implicit time integrator.
///
/// The Newmark kinematics
///
/// ```text
/// u_{n+1} = u_n + Δt·v_n + Δt²·[(1/2 − β)·a_n + β·a_{n+1}]
/// v_{n+1} = v_n + Δt·[(1 − γ)·a_n + γ·a_{n+1}]
/// ```
///
/// are combined with the α-weighted equation of motion
///
/// ```text
/// M·a_{n+1} + (1+α)(C·v_{n+1} + K·u_{n+1}) − α(C·v_n + K·u_n) = (1+α)·f_{n+1} − α·f_n
/// ```
///
/// where α = [`hht_alpha`](Self::hht_alpha) ∈ [−1/3, 0] (α = 0 is the classical Newmark
/// method) and `C = α_R·M + β_R·K` is the Rayleigh damping of the [`DynamicsSystem`].
/// Substituting the kinematics for `u_{n+1}`, `v_{n+1}` leaves one linear system for the new
/// acceleration,
///
/// ```text
/// [M + (1+α)(γΔt·C + βΔt²·K)]·a_{n+1}
///     = (1+α)·f_{n+1} − α·f_n − (1+α)(K·ũ + C·ṽ) + α(K·u_n + C·v_n)
/// ```
///
/// with the predictors `ũ = u_n + Δt·v_n + Δt²(1/2 − β)·a_n` and `ṽ = v_n + Δt(1 − γ)·a_n`.
/// The full (generally non-diagonal) effective matrix is assembled and factorised with a
/// dense LU decomposition, so coupled systems are integrated exactly as the method prescribes
/// and the scheme has the stability properties of the Newmark/HHT family: with 2β ≥ γ ≥ 1/2
/// (classical Newmark) or the optimal HHT-α parameters it is unconditionally stable for linear
/// systems; otherwise see [`critical_dt`](Self::critical_dt).
///
/// The factorisation costs `O(n³)`.  For many steps with a fixed `Δt` build a
/// [`NewmarkStepper`] once with [`stepper`](Self::stepper) and reuse it; each step then costs
/// `O(n²)`.
///
/// For large sparse systems use the CSR-based `crate::dynamic` module.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NewmarkIntegrator {
    /// Newmark β parameter.
    pub beta: f64,
    /// Newmark γ parameter.
    pub gamma: f64,
    /// HHT-α parameter in `[-1/3, 0]`; `0` gives the classical Newmark method.
    pub hht_alpha: f64,
}

impl NewmarkIntegrator {
    /// Average-acceleration method: β = 1/4, γ = 1/2.
    ///
    /// Unconditionally stable for linear systems (spectral radius ≤ 1 for all dt) and
    /// non-dissipative.
    pub fn new_average_acceleration() -> Self {
        Self {
            beta: 0.25,
            gamma: 0.5,
            hht_alpha: 0.0,
        }
    }

    /// Linear-acceleration method: β = 1/6, γ = 1/2.
    ///
    /// Conditionally stable: requires ω·dt ≤ 2√3.
    pub fn new_linear_acceleration() -> Self {
        Self {
            beta: 1.0 / 6.0,
            gamma: 0.5,
            hht_alpha: 0.0,
        }
    }

    /// HHT-α method (Hilber, Hughes & Taylor, 1977): β = (1−α)²/4, γ = 1/2 − α.
    ///
    /// For α ∈ [−1/3, 0] the method is unconditionally stable, second-order accurate, and
    /// dissipates high-frequency modes (asymptotic spectral radius (1+α)/(1−α)) while leaving
    /// the low-frequency response almost untouched.  α = 0 recovers the average-acceleration
    /// method.  Note γ = 1/2 − α ≥ 1/2: a *negative* α gives positive numerical damping.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when `alpha` is not in `[-1/3, 0]` (this includes NaN).
    pub fn new_hht_alpha(alpha: f64) -> Result<Self> {
        if !(-1.0 / 3.0..=0.0).contains(&alpha) {
            return Err(Error::InvalidInput(format!(
                "HHT-α must lie in [-1/3, 0], got {alpha}"
            )));
        }
        Ok(Self {
            beta: (1.0 - alpha).powi(2) / 4.0,
            gamma: 0.5 - alpha,
            hht_alpha: alpha,
        })
    }

    /// Check that the integrator parameters define a solvable step.
    fn validate(&self) -> Result<()> {
        if !(self.beta.is_finite() && self.beta >= 0.0) {
            return Err(Error::InvalidInput(format!(
                "Newmark β must be finite and non-negative, got {}",
                self.beta
            )));
        }
        if !(self.gamma.is_finite() && self.gamma >= 0.0) {
            return Err(Error::InvalidInput(format!(
                "Newmark γ must be finite and non-negative, got {}",
                self.gamma
            )));
        }
        if !(self.hht_alpha.is_finite() && self.hht_alpha > -1.0) {
            return Err(Error::InvalidInput(format!(
                "HHT-α must be finite and greater than -1, got {}",
                self.hht_alpha
            )));
        }
        Ok(())
    }

    /// Bind the integrator to a system and a time step and factorise the effective matrix.
    ///
    /// The returned [`NewmarkStepper`] can advance the state many times at `O(n²)` per step.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidInput`] when `dt` is not positive and finite, when the integrator
    ///   parameters are unusable, or when the system is inconsistent (sizes, negative or
    ///   non-finite masses, non-finite stiffness or damping).
    /// * [`Error::SingularSystem`] when the effective matrix
    ///   `M + (1+α)(γΔt·C + βΔt²·K)` is singular, e.g. a rigid-body mode of a massless DOF.
    pub fn stepper<'a>(&self, sys: &'a DynamicsSystem, dt: f64) -> Result<NewmarkStepper<'a>> {
        self.validate()?;
        check_time_step(dt)?;
        sys.validate()?;
        let weight = 1.0 + self.hht_alpha;
        let m_coeff = 1.0 + weight * self.gamma * dt * sys.damping_alpha;
        let k_coeff = weight * (self.beta * dt * dt + self.gamma * dt * sys.damping_beta);
        let effective = sys.mass_stiffness_combination(m_coeff, k_coeff);
        let factor = LuFactor::new(&effective).ok_or(Error::SingularSystem)?;
        Ok(NewmarkStepper {
            sys,
            integrator: *self,
            dt,
            factor,
        })
    }

    /// Advance the state by one step of size `dt`, treating the load as constant over the step.
    ///
    /// `f_ext` is the external load `f_{n+1}`; it is also used for the `f_n` of the HHT
    /// weighting, which is exact for loads that are constant over the step and for the
    /// classical Newmark method (α = 0).  For time-varying loads with α ≠ 0 build a
    /// [`NewmarkStepper`] with [`stepper`](Self::stepper) and pass both loads to
    /// [`NewmarkStepper::step_with_loads`] to retain second-order accuracy.
    ///
    /// `a` must be the acceleration consistent with `u`, `v` and the load at the start of the
    /// step (see [`DynamicsSystem::initial_acceleration`]).  On error the state is unchanged.
    ///
    /// # Errors
    ///
    /// See [`stepper`](Self::stepper); additionally [`Error::InvalidInput`] when a state or load
    /// vector does not have `sys.n_dof` entries or the update becomes non-finite.
    pub fn step(
        &self,
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) -> Result<()> {
        self.stepper(sys, dt)?
            .step_with_loads(u, v, a, f_ext, f_ext)
    }

    /// Normalised amplification matrix of one step for the undamped modal equation
    /// `u'' + ω²u = 0` with `ω·Δt = omega_dt`, acting on `(u, Δt·v, Δt²·a)`.
    fn modal_amplification(&self, omega_dt: f64) -> [[f64; 3]; 3] {
        let k = omega_dt * omega_dt;
        let weight = 1.0 + self.hht_alpha;
        let denominator = 1.0 + weight * self.beta * k;
        let advance = |state: [f64; 3]| -> [f64; 3] {
            let [u, v, a] = state;
            let u_pred = u + v + (0.5 - self.beta) * a;
            let v_pred = v + (1.0 - self.gamma) * a;
            let a_new = (-weight * k * u_pred + self.hht_alpha * k * u) / denominator;
            [
                u_pred + self.beta * a_new,
                v_pred + self.gamma * a_new,
                a_new,
            ]
        };
        let columns = [
            advance([1.0, 0.0, 0.0]),
            advance([0.0, 1.0, 0.0]),
            advance([0.0, 0.0, 1.0]),
        ];
        let mut matrix = [[0.0; 3]; 3];
        for (j, column) in columns.iter().enumerate() {
            for (i, &entry) in column.iter().enumerate() {
                matrix[i][j] = entry;
            }
        }
        matrix
    }

    /// Whether every eigenvalue of the modal amplification matrix lies in the closed unit disk.
    fn is_modally_stable(&self, omega_dt: f64) -> bool {
        schur_cohn_cubic_stable(&self.modal_amplification(omega_dt))
    }

    /// Largest `ω·Δt` up to which the undamped modal recurrence is stable.
    ///
    /// Returns `f64::INFINITY` for unconditional stability and `0.0` when no step is stable.
    fn critical_omega_dt(&self) -> f64 {
        if self.hht_alpha == 0.0 {
            // Classical Newmark (Hughes, The Finite Element Method, §9.1).
            if self.gamma < 0.5 {
                return 0.0;
            }
            if 2.0 * self.beta >= self.gamma {
                return f64::INFINITY;
            }
            return 1.0 / (0.5 * self.gamma - self.beta).sqrt();
        }
        // HHT-α with its optimal parameters is unconditionally stable for α ∈ [-1/3, 0].
        let optimal_gamma = 0.5 - self.hht_alpha;
        let optimal_beta = (1.0 - self.hht_alpha).powi(2) / 4.0;
        if (-1.0 / 3.0..=0.0).contains(&self.hht_alpha)
            && (self.gamma - optimal_gamma).abs() <= 1.0e-12
            && self.beta >= optimal_beta - 1.0e-12
        {
            return f64::INFINITY;
        }
        // Any other (β, γ, α) combination: the plain Newmark criterion does not apply.
        self.scan_critical_omega_dt()
    }

    /// Largest `ω·Δt` up to which the modal amplification matrix keeps its eigenvalues in the
    /// closed unit disk, found by scanning a geometric grid for the first unstable step and
    /// bisecting the bracket.  Valid for every `(β, γ, α)`; `f64::INFINITY` when the whole
    /// grid up to `ω·Δt = 10⁸` is stable and `0.0` when already `ω·Δt = 10⁻³` is not.
    fn scan_critical_omega_dt(&self) -> f64 {
        const FIRST: f64 = 1.0e-3;
        const LAST: f64 = 1.0e8;
        const GROWTH: f64 = 1.02;
        if !self.is_modally_stable(FIRST) {
            return 0.0;
        }
        let mut stable = FIRST;
        while stable < LAST {
            let next = (stable * GROWTH).min(LAST);
            if !self.is_modally_stable(next) {
                let (mut low, mut high) = (stable, next);
                for _ in 0..60 {
                    let mid = 0.5 * (low + high);
                    if self.is_modally_stable(mid) {
                        low = mid;
                    } else {
                        high = mid;
                    }
                }
                return low;
            }
            stable = next;
        }
        f64::INFINITY
    }

    /// Stability limit Δt_crit of the scheme for a system whose highest natural frequency is
    /// `omega_max`.
    ///
    /// Obtained from the undamped modal recurrence, which is also a conservative bound for
    /// positive Rayleigh damping.  For the classical Newmark family (α = 0):
    ///
    /// - `2β ≥ γ ≥ 1/2`: unconditionally stable, returns `f64::INFINITY`;
    /// - `γ ≥ 1/2` and `β < γ/2`: conditionally stable with `ω·Δt ≤ 1/√(γ/2 − β)`
    ///   (2 for central difference, 2√3 for linear acceleration);
    /// - `γ < 1/2`: no time step is stable, returns `0.0`.
    ///
    /// The optimal HHT-α parameters of [`new_hht_alpha`](Self::new_hht_alpha) are
    /// unconditionally stable.  For any other combination with α ≠ 0 the limit is found from the
    /// spectral radius of the amplification matrix.
    ///
    /// A zero `omega_max` imposes no limit (`f64::INFINITY`); a negative or NaN `omega_max`
    /// yields NaN.
    pub fn critical_dt(&self, omega_max: f64) -> f64 {
        if omega_max.is_nan() || omega_max < 0.0 {
            return f64::NAN;
        }
        if omega_max == 0.0 {
            return f64::INFINITY;
        }
        self.critical_omega_dt() / omega_max
    }
}

/// Whether the roots of the characteristic polynomial of the 3×3 matrix `a` all lie in the
/// closed unit disk (Schur–Cohn/Jury test for a real cubic, with a small relative tolerance).
fn schur_cohn_cubic_stable(a: &[[f64; 3]; 3]) -> bool {
    let trace = a[0][0] + a[1][1] + a[2][2];
    let minors = (a[0][0] * a[1][1] - a[0][1] * a[1][0])
        + (a[0][0] * a[2][2] - a[0][2] * a[2][0])
        + (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    // p(λ) = λ³ + c2·λ² + c1·λ + c0
    let (c2, c1, c0) = (-trace, minors, -det);
    if !(c2.is_finite() && c1.is_finite() && c0.is_finite()) {
        return false;
    }
    let tol = 1.0e-9 * (1.0 + c2.abs() + c1.abs() + c0.abs());
    let p_at_one = 1.0 + c2 + c1 + c0;
    let p_at_minus_one = -1.0 + c2 - c1 + c0;
    p_at_one >= -tol
        && -p_at_minus_one >= -tol
        && c0.abs() <= 1.0 + tol
        && (c0 * c0 - 1.0).abs() >= (c0 * c2 - c1).abs() - tol
}

/// A [`NewmarkIntegrator`] bound to a [`DynamicsSystem`] and a time step, with the effective
/// matrix already factorised.
///
/// Created by [`NewmarkIntegrator::stepper`].  Each call to
/// [`step_with_loads`](Self::step_with_loads) then costs `O(n²)`.
#[derive(Debug)]
pub struct NewmarkStepper<'a> {
    /// The system being integrated.
    sys: &'a DynamicsSystem,
    /// Integration parameters.
    integrator: NewmarkIntegrator,
    /// Time step.
    dt: f64,
    /// LU factors of `M + (1+α)(γΔt·C + βΔt²·K)`.
    factor: LuFactor,
}

impl NewmarkStepper<'_> {
    /// Advance the state by one step.
    ///
    /// `f_prev` is the external load `f_n` at the start of the step and `f_next` the load
    /// `f_{n+1}` at its end; `f_prev` only matters for HHT-α (α ≠ 0).  `a` must be consistent
    /// with `u`, `v` and `f_prev` (see [`DynamicsSystem::initial_acceleration`]).  On error the
    /// state is left unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when a vector does not have `n_dof` entries or the update
    /// becomes non-finite (non-finite state or load).
    pub fn step_with_loads(
        &self,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_prev: &[f64],
        f_next: &[f64],
    ) -> Result<()> {
        let sys = self.sys;
        let n = sys.n_dof;
        check_len("displacement vector", u.len(), n)?;
        check_len("velocity vector", v.len(), n)?;
        check_len("acceleration vector", a.len(), n)?;
        check_len("previous load vector", f_prev.len(), n)?;
        check_len("load vector", f_next.len(), n)?;

        let NewmarkIntegrator {
            beta,
            gamma,
            hht_alpha: alpha,
        } = self.integrator;
        let dt = self.dt;
        let dt2 = dt * dt;
        let weight = 1.0 + alpha;

        let u_pred: Vec<f64> = (0..n)
            .map(|i| u[i] + dt * v[i] + dt2 * (0.5 - beta) * a[i])
            .collect();
        let v_pred: Vec<f64> = (0..n).map(|i| v[i] + dt * (1.0 - gamma) * a[i]).collect();
        let k_u_pred = sys.internal_force(&u_pred);
        let c_v_pred = sys.rayleigh_damping_force(&v_pred);
        let (k_u_old, c_v_old) = if alpha == 0.0 {
            (vec![0.0; n], vec![0.0; n])
        } else {
            (sys.internal_force(u), sys.rayleigh_damping_force(v))
        };

        let rhs: Vec<f64> = (0..n)
            .map(|i| {
                weight * f_next[i] - alpha * f_prev[i] - weight * (k_u_pred[i] + c_v_pred[i])
                    + alpha * (k_u_old[i] + c_v_old[i])
            })
            .collect();
        let a_new = self.factor.solve(&rhs);
        if a_new.iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput(
                "Newmark update is not finite: check the state and the load for NaN or infinity"
                    .into(),
            ));
        }

        let u_new: Vec<f64> = (0..n).map(|i| u_pred[i] + beta * dt2 * a_new[i]).collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + gamma * dt * a_new[i]).collect();
        *u = u_new;
        *v = v_new;
        *a = a_new;
        Ok(())
    }

    /// Advance the state by one step with the load held constant over the step
    /// (`f_prev = f_next = f_ext`); see [`NewmarkIntegrator::step`].
    ///
    /// # Errors
    ///
    /// See [`step_with_loads`](Self::step_with_loads).
    pub fn step(
        &self,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
    ) -> Result<()> {
        self.step_with_loads(u, v, a, f_ext, f_ext)
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

/// Bathe composite time integration method (Bathe & Noh, 2012).
///
/// Each step of size Δt is split at the ratio p ∈ (0, 1) into two implicit sub-steps:
///
/// 1. **Trapezoidal rule** over \[t_n, t_n + pΔt\] (Newmark β = 1/4, γ = 1/2) gives the
///    intermediate state `(u_p, v_p)`.
/// 2. **Three-point backward difference** over \[t_n + pΔt, t_{n+1}\] uses the states at
///    `t_n`, `t_n + pΔt` and `t_{n+1}`:
///
///    ```text
///    v_{n+1} = c1·u_n + c2·u_p + c3·u_{n+1}
///    a_{n+1} = c1·v_n + c2·v_p + c3·v_{n+1}
///    c1 = (1 − p)/(p·Δt),  c2 = −1/(p(1 − p)Δt),  c3 = (2 − p)/((1 − p)Δt)
///    ```
///
///    Together with the equation of motion at `t_{n+1}` this yields
///    `(c3²·M + c3·C + K)·u_{n+1} = f_{n+1} − M(c1·v_n + c2·v_p + c3(c1·u_n + c2·u_p))
///    − C(c1·u_n + c2·u_p)`.
///
/// The composite scheme is second-order accurate, unconditionally stable for linear systems for
/// every p ∈ (0, 1), and dissipates high-frequency content strongly (the spectral radius tends
/// to zero as ω·Δt → ∞) while preserving the low-frequency response; p = 1/2 is the standard
/// choice.  Both sub-steps solve the full (generally non-diagonal) effective system with a dense
/// LU decomposition.
///
/// For many steps with a fixed `Δt` build a [`BatheStepper`] with [`stepper`](Self::stepper).
///
/// Reference: K.-J. Bathe, G. Noh, "Insight into an implicit time integration scheme for
/// structural dynamics", Computers & Structures 98–99 (2012) 1–6.
#[derive(Debug, Clone, Copy, PartialEq)]
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
    ///
    /// # Panics
    ///
    /// Panics if `split_ratio` is not in the open interval (0, 1).
    pub fn with_split_ratio(split_ratio: f64) -> Self {
        assert!(
            split_ratio > 0.0 && split_ratio < 1.0,
            "split_ratio must be in (0, 1)"
        );
        Self { split_ratio }
    }

    /// Bind the integrator to a system and a time step and factorise both sub-step matrices.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidInput`] when `split_ratio` is not in (0, 1), when `dt` is not positive
    ///   and finite, or when the system is inconsistent (see [`NewmarkIntegrator::stepper`]).
    /// * [`Error::SingularSystem`] when one of the sub-step matrices is singular.
    pub fn stepper<'a>(&self, sys: &'a DynamicsSystem, dt: f64) -> Result<BatheStepper<'a>> {
        let p = self.split_ratio;
        if !(p.is_finite() && p > 0.0 && p < 1.0) {
            return Err(Error::InvalidInput(format!(
                "Bathe split ratio must be in (0, 1), got {p}"
            )));
        }
        check_time_step(dt)?;
        let first = NewmarkIntegrator::new_average_acceleration().stepper(sys, p * dt)?;
        let c1 = (1.0 - p) / (p * dt);
        let c2 = -1.0 / (p * (1.0 - p) * dt);
        let c3 = (2.0 - p) / ((1.0 - p) * dt);
        let m_coeff = c3 * c3 + c3 * sys.damping_alpha;
        let k_coeff = 1.0 + c3 * sys.damping_beta;
        let effective = sys.mass_stiffness_combination(m_coeff, k_coeff);
        let second = LuFactor::new(&effective).ok_or(Error::SingularSystem)?;
        Ok(BatheStepper {
            sys,
            first,
            coefficients: [c1, c2, c3],
            second,
        })
    }

    /// Advance the state by one Bathe composite step, treating the load as constant over the
    /// step.
    ///
    /// `f_ext` serves as the load at both `t_n + p·Δt` and `t_{n+1}`.  For time-varying loads use
    /// [`BatheStepper::step_with_loads`].  `a` must be the acceleration consistent with `u`,
    /// `v` and the load at the start of the step (see [`DynamicsSystem::initial_acceleration`]).
    /// On error the state is unchanged.
    ///
    /// # Errors
    ///
    /// See [`stepper`](Self::stepper); additionally [`Error::InvalidInput`] when a state or load
    /// vector does not have `sys.n_dof` entries or the update becomes non-finite.
    pub fn step(
        &self,
        sys: &DynamicsSystem,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
        dt: f64,
    ) -> Result<()> {
        self.stepper(sys, dt)?
            .step_with_loads(u, v, a, f_ext, f_ext)
    }
}

impl Default for BatheIntegrator {
    fn default() -> Self {
        Self::new()
    }
}

/// A [`BatheIntegrator`] bound to a [`DynamicsSystem`] and a time step, with both sub-step
/// matrices already factorised.
///
/// Created by [`BatheIntegrator::stepper`].
#[derive(Debug)]
pub struct BatheStepper<'a> {
    /// The system being integrated.
    sys: &'a DynamicsSystem,
    /// Trapezoidal first sub-step over `p·Δt`.
    first: NewmarkStepper<'a>,
    /// Backward-difference coefficients `[c1, c2, c3]` of the second sub-step.
    coefficients: [f64; 3],
    /// LU factors of `c3²·M + c3·C + K`.
    second: LuFactor,
}

impl BatheStepper<'_> {
    /// Advance the state by one composite step.
    ///
    /// `f_mid` is the external load at `t_n + p·Δt` (first sub-step) and `f_end` the load at
    /// `t_{n+1}` (second sub-step).  `a` must be consistent with `u`, `v` and the load at
    /// `t_n`.  On error the state is left unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when a vector does not have `n_dof` entries or the update
    /// becomes non-finite.
    pub fn step_with_loads(
        &self,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_mid: &[f64],
        f_end: &[f64],
    ) -> Result<()> {
        let sys = self.sys;
        let n = sys.n_dof;
        check_len("displacement vector", u.len(), n)?;
        check_len("velocity vector", v.len(), n)?;
        check_len("acceleration vector", a.len(), n)?;
        check_len("mid-step load vector", f_mid.len(), n)?;
        check_len("end-of-step load vector", f_end.len(), n)?;
        let [c1, c2, c3] = self.coefficients;

        // Sub-step 1: trapezoidal rule over p·Δt, on copies so that an error leaves the caller's
        // state untouched.
        let mut u_mid = u.clone();
        let mut v_mid = v.clone();
        let mut a_mid = a.clone();
        self.first
            .step_with_loads(&mut u_mid, &mut v_mid, &mut a_mid, f_mid, f_mid)?;

        // Sub-step 2: three-point backward difference over (1 − p)·Δt.
        let history_u: Vec<f64> = (0..n).map(|i| c1 * u[i] + c2 * u_mid[i]).collect();
        let history_v: Vec<f64> = (0..n).map(|i| c1 * v[i] + c2 * v_mid[i]).collect();
        let c_history_u = sys.rayleigh_damping_force(&history_u);
        let rhs: Vec<f64> = (0..n)
            .map(|i| {
                f_end[i] - sys.mass_diag[i] * (history_v[i] + c3 * history_u[i]) - c_history_u[i]
            })
            .collect();
        let u_new = self.second.solve(&rhs);
        if u_new.iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput(
                "Bathe update is not finite: check the state and the load for NaN or infinity"
                    .into(),
            ));
        }
        let v_new: Vec<f64> = (0..n).map(|i| history_u[i] + c3 * u_new[i]).collect();
        let a_new: Vec<f64> = (0..n).map(|i| history_v[i] + c3 * v_new[i]).collect();
        *u = u_new;
        *v = v_new;
        *a = a_new;
        Ok(())
    }

    /// Advance the state by one composite step with the load held constant over the step
    /// (`f_mid = f_end = f_ext`); see [`BatheIntegrator::step`].
    ///
    /// # Errors
    ///
    /// See [`step_with_loads`](Self::step_with_loads).
    pub fn step(
        &self,
        u: &mut Vec<f64>,
        v: &mut Vec<f64>,
        a: &mut Vec<f64>,
        f_ext: &[f64],
    ) -> Result<()> {
        self.step_with_loads(u, v, a, f_ext, f_ext)
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
/// `n_dim` is the number of degrees of freedom per node; DOF `i` belongs to direction
/// `i % n_dim`.  The diagonal entries are scaled *per direction* so that the mass carried in
/// every direction is preserved:
///
/// ```text
/// M_lumped[i] = M[i][i] · S_d / D_d,    d = i mod n_dim
/// S_d = Σ_{i ≡ d, j ≡ d} M[i][j]        (mass of direction d)
/// D_d = Σ_{i ≡ d} M[i][i]               (diagonal sum of direction d)
/// ```
///
/// Hence `Σ_{i ≡ d} M_lumped[i] = S_d` for every direction.  For a consistent mass matrix built
/// from a scalar element mass matrix, such as the 4-node tetrahedron with `n_dim = 3`, every
/// lumped entry is `ρV/4` and the lumped entries sum to `3ρV` (the element mass once per
/// direction).  A direction whose diagonal sum vanishes is lumped to zero, and `n_dim == 0` is
/// treated like `1`.
///
/// This method often gives better results than row-sum lumping for higher-order elements.
pub fn lump_mass_hrz(consistent_mass: &[Vec<f64>], n_dim: usize) -> Vec<f64> {
    if consistent_mass.is_empty() {
        return Vec::new();
    }
    let directions = n_dim.max(1);
    let mut direction_mass = vec![0.0_f64; directions];
    let mut diagonal_sum = vec![0.0_f64; directions];
    for (i, row) in consistent_mass.iter().enumerate() {
        let d = i % directions;
        direction_mass[d] += row
            .iter()
            .enumerate()
            .filter(|(j, _)| j % directions == d)
            .map(|(_, m)| m)
            .sum::<f64>();
        diagonal_sum[d] += row.get(i).copied().unwrap_or(0.0);
    }
    consistent_mass
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let d = i % directions;
            if diagonal_sum[d].abs() < 1e-30 {
                0.0
            } else {
                row.get(i).copied().unwrap_or(0.0) * direction_mass[d] / diagonal_sum[d]
            }
        })
        .collect()
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
    pub fn step(sys: &DynamicsSystem, u: &mut [f64], v_half: &mut [f64], f_ext: &[f64], dt: f64) {
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

/// Estimate an upper bound of the maximum natural frequency of `K·φ = ω²·M·φ` for a lumped
/// (diagonal) mass matrix.
///
/// The bound is Gershgorin's circle theorem applied to `M⁻¹K`: every eigenvalue satisfies
/// `ω² ≤ max_i (Σ_j |K_ij|) / M_ii`, so `ω_max ≤ max_i sqrt(Σ_j |K_ij| / M_ii)`.  Using only the
/// diagonal ratio `K_ii/M_ii` underestimates `ω_max` for any coupled stiffness (a uniform
/// spring chain has `K_ii = 2k` but `ω_max² → 4k/m`) and so yields a time step that is not
/// stable.  DOFs with a (numerically) zero mass are skipped.
///
/// Used to determine the critical time step for explicit methods.
pub fn estimate_omega_max(mass_diag: &[f64], stiffness: &[Vec<f64>]) -> f64 {
    mass_diag
        .iter()
        .zip(stiffness.iter())
        .filter(|(mass, _)| **mass > 1e-30)
        .fold(0.0_f64, |omega_max, (mass, row)| {
            let row_sum: f64 = row.iter().map(|k| k.abs()).sum();
            omega_max.max((row_sum / mass).sqrt())
        })
}

/// Compute the safe time step for explicit integration with a safety factor.
///
/// dt_safe = safety_factor * 2 / ω_max, with ω_max the Gershgorin bound of
/// [`estimate_omega_max`]; any `safety_factor` ≤ 1 therefore keeps the central-difference
/// scheme stable.
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
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("Newmark step");
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
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("Newmark step");
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
            bathe
                .step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("Bathe step");
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
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("Newmark step");
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
        let hht = NewmarkIntegrator::new_hht_alpha(0.0).expect("alpha in range");
        assert!((hht.beta - 0.25).abs() < 1e-12);
        assert!((hht.gamma - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_hht_alpha_parameters() {
        // Hilber–Hughes–Taylor: γ = 1/2 − α, β = (1 − α)²/4.  With α = −0.1 this is
        // γ = 0.6 ≥ 1/2 (positive numerical damping) and β = 1.1²/4 = 0.3025.
        let hht = NewmarkIntegrator::new_hht_alpha(-0.1).expect("alpha in range");
        assert!((hht.gamma - 0.6).abs() < 1e-12, "gamma = {}", hht.gamma);
        assert!((hht.beta - 0.3025).abs() < 1e-12, "beta = {}", hht.beta);
        assert!((hht.hht_alpha + 0.1).abs() < 1e-15);
    }

    #[test]
    fn test_hht_alpha_rejects_out_of_range_values() {
        for alpha in [0.1, -0.34, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    NewmarkIntegrator::new_hht_alpha(alpha),
                    Err(Error::InvalidInput(_))
                ),
                "alpha = {alpha} must be rejected"
            );
        }
        assert!(NewmarkIntegrator::new_hht_alpha(-1.0 / 3.0).is_ok());
        assert!(NewmarkIntegrator::new_hht_alpha(0.0).is_ok());
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

        let hht = NewmarkIntegrator::new_hht_alpha(-0.1).expect("alpha in range");
        for _ in 0..2000 {
            hht.step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("HHT step");
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
            integrator
                .step(&sys, &mut u, &mut v, &mut a, &f_ext, dt)
                .expect("Newmark step");
        }
        assert!(u[0].abs() < 10.0);
        assert!(u[1].abs() < 10.0);
    }
}

#[cfg(test)]
mod stability_tests;
