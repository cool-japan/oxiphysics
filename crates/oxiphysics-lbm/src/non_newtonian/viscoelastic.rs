// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Viscoelastic fluid models: Oldroyd-B, Maxwell, PTT, Giesekus,
//! Johnson-Segalman, Rolie-Poly, and viscoelastic relaxation field.

use super::{CS2, LocalViscosityModel, NonNewtonianFluid};

// ---------------------------------------------------------------------------
// Oldroyd-B viscoelastic fluid model
// ---------------------------------------------------------------------------

/// Oldroyd-B viscoelastic fluid model.
///
/// The Oldroyd-B model captures both viscous (Newtonian) and elastic
/// (Maxwell) contributions to the stress:
///
/// ```text
/// T + lambda_1 * dT/dt = 2 * eta * (D + lambda_2 * dD/dt)
/// ```
///
/// where:
/// - `T` is the extra stress tensor
/// - `D` is the rate-of-strain tensor
/// - `lambda_1` is the relaxation time
/// - `lambda_2` is the retardation time
/// - `eta` is the total viscosity
///
/// In the limit `lambda_2 = 0` it reduces to the Upper-Convected Maxwell (UCM) model.
/// In the limit `lambda_1 = 0` it reduces to a Newtonian fluid.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct OldroydB {
    /// Total (zero-shear-rate) viscosity eta (polymer + solvent).
    pub eta: f64,
    /// Relaxation time lambda_1 (s).
    pub lambda_1: f64,
    /// Retardation time lambda_2 (s). Must satisfy 0 <= lambda_2 < lambda_1.
    pub lambda_2: f64,
    /// Solvent viscosity eta_s.
    pub eta_s: f64,
}

#[allow(dead_code)]
impl OldroydB {
    /// Create a new Oldroyd-B model.
    pub fn new(eta: f64, lambda_1: f64, lambda_2: f64, eta_s: f64) -> Self {
        Self {
            eta,
            lambda_1,
            lambda_2,
            eta_s,
        }
    }

    /// Polymer (elastic) viscosity: `eta_p = eta - eta_s`.
    pub fn eta_polymer(&self) -> f64 {
        (self.eta - self.eta_s).max(0.0)
    }

    /// Weissenberg number: `Wi = lambda_1 * gamma_dot`.
    pub fn weissenberg(&self, shear_rate: f64) -> f64 {
        self.lambda_1 * shear_rate.abs()
    }

    /// Effective steady-state shear viscosity for simple shear.
    ///
    /// For steady shear the Oldroyd-B shear viscosity is constant = `eta`.
    pub fn steady_shear_viscosity(&self) -> f64 {
        self.eta
    }

    /// First normal stress coefficient `Psi_1 = 2 * eta_p * lambda_1`.
    pub fn first_normal_stress_coefficient(&self) -> f64 {
        2.0 * self.eta_polymer() * self.lambda_1
    }

    /// Second normal stress coefficient `Psi_2 = -2 * eta_p * lambda_2` (simplified).
    pub fn second_normal_stress_coefficient(&self) -> f64 {
        -2.0 * self.eta_polymer() * self.lambda_2
    }

    /// First normal stress difference `N1 = Psi_1 * gamma_dot^2`.
    pub fn first_normal_stress_difference(&self, shear_rate: f64) -> f64 {
        self.first_normal_stress_coefficient() * shear_rate * shear_rate
    }

    /// Relaxation modulus in the linear viscoelastic regime:
    /// `G(t) = eta_p / lambda_1 * exp(-t / lambda_1)`
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        if self.lambda_1 < 1e-30 {
            return 0.0;
        }
        (self.eta_polymer() / self.lambda_1) * (-t / self.lambda_1).exp()
    }

    /// Storage modulus G'(omega): `G' = eta_p * omega^2 * lambda_1^2 / (1 + (omega*lambda_1)^2)`.
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let wl = omega * self.lambda_1;
        self.eta_polymer() * omega * wl / (1.0 + wl * wl)
    }

    /// Loss modulus G''(omega): `G'' = eta_s * omega + eta_p * omega / (1 + (omega*lambda_1)^2)`.
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let wl = omega * self.lambda_1;
        self.eta_s * omega + self.eta_polymer() * omega / (1.0 + wl * wl)
    }

    /// Advance the conformation tensor one step using a simple Euler scheme.
    ///
    /// Returns `(Axx, Axy, Ayy)` after one step `dt`.
    pub fn step_conformation_tensor(
        &self,
        axx: f64,
        axy: f64,
        ayy: f64,
        shear_rate: f64,
        dt: f64,
    ) -> (f64, f64, f64) {
        if self.lambda_1 < 1e-30 {
            return (1.0, 0.0, 1.0);
        }
        let daxx = 2.0 * shear_rate * axy + (1.0 - axx) / self.lambda_1;
        let daxy = shear_rate * ayy + (0.0 - axy) / self.lambda_1;
        let dayy = (1.0 - ayy) / self.lambda_1;
        (axx + dt * daxx, axy + dt * daxy, ayy + dt * dayy)
    }
}

// ---------------------------------------------------------------------------
// Maxwell Viscoelastic Fluid
// ---------------------------------------------------------------------------

/// Upper-Convected Maxwell (UCM) viscoelastic fluid model.
///
/// The UCM model captures elastic stress relaxation:
///
/// ```text
/// T + lambda * delta_upper(T)/delta_t = 2 * eta * D
/// ```
///
/// where `lambda` is the relaxation time, `eta` is the viscosity, and `D` is the
/// rate-of-strain tensor.  For simple shear the steady-state viscosity equals `eta`
/// and the first normal stress difference is `N1 = 2 * eta * lambda * gamma^2`.
///
/// The complex viscosity at angular frequency `omega` is:
///
/// ```text
/// eta*(omega) = eta / (1 + i*omega*lambda)
/// ```
///
/// giving storage modulus `G'(omega) = eta * omega^2 * lambda^2 / (1 + omega^2 * lambda^2)`
/// and loss modulus  `G''(omega) = eta * omega / (1 + omega^2 * lambda^2)`.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct MaxwellFluid {
    /// Relaxation time lambda (s).
    pub relaxation_time: f64,
    /// Elastic modulus G (Pa) -- relates to viscosity via G = eta / lambda.
    pub elastic_modulus: f64,
}

#[allow(dead_code)]
impl MaxwellFluid {
    /// Create a new Maxwell fluid with relaxation time `lambda` and modulus `G`.
    pub fn new(relaxation_time: f64, elastic_modulus: f64) -> Self {
        Self {
            relaxation_time,
            elastic_modulus,
        }
    }

    /// Zero-shear viscosity: `eta = G * lambda`.
    pub fn viscosity(&self) -> f64 {
        self.elastic_modulus * self.relaxation_time
    }

    /// Weissenberg number: `Wi = lambda * |gamma|`.
    pub fn weissenberg(&self, shear_rate: f64) -> f64 {
        self.relaxation_time * shear_rate.abs()
    }

    /// Storage modulus at angular frequency `omega`:
    ///
    /// `G'(omega) = G * (omega*lambda)^2 / (1 + (omega*lambda)^2)`
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let wl = omega * self.relaxation_time;
        let wl2 = wl * wl;
        self.elastic_modulus * wl2 / (1.0 + wl2)
    }

    /// Loss modulus at angular frequency `omega`:
    ///
    /// `G''(omega) = G * omega*lambda / (1 + (omega*lambda)^2)`
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let wl = omega * self.relaxation_time;
        let wl2 = wl * wl;
        self.elastic_modulus * wl / (1.0 + wl2)
    }

    /// Magnitude of complex viscosity `|eta*(omega)| = eta / sqrt(1 + (omega*lambda)^2)`.
    pub fn complex_viscosity_magnitude(&self, omega: f64) -> f64 {
        let wl = omega * self.relaxation_time;
        self.viscosity() / (1.0 + wl * wl).sqrt()
    }

    /// Phase angle `delta = arctan(G'' / G')` (rad).  delta = pi/2 for purely viscous,
    /// delta = 0 for purely elastic.
    pub fn phase_angle(&self, omega: f64) -> f64 {
        let g_prime = self.storage_modulus(omega);
        let g_double_prime = self.loss_modulus(omega);
        if g_prime < 1e-30 {
            return std::f64::consts::FRAC_PI_2;
        }
        (g_double_prime / g_prime).atan()
    }

    /// First normal stress coefficient `Psi_1 = 2 * eta * lambda`.
    pub fn first_normal_stress_coefficient(&self) -> f64 {
        2.0 * self.viscosity() * self.relaxation_time
    }

    /// First normal stress difference in steady shear:
    /// `N1 = Psi_1 * gamma^2 = 2 * eta * lambda * gamma^2`.
    pub fn first_normal_stress_difference(&self, shear_rate: f64) -> f64 {
        self.first_normal_stress_coefficient() * shear_rate * shear_rate
    }

    /// Effective LBM relaxation time for the viscous part.
    ///
    /// `tau_LBM = 0.5 + eta / cs^2`
    pub fn lbm_relaxation_time(&self) -> f64 {
        0.5 + self.viscosity() / CS2
    }

    /// Relaxation modulus: `G(t) = G * exp(-t / lambda)`.
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        if self.relaxation_time < 1e-30 {
            return 0.0;
        }
        self.elastic_modulus * (-t / self.relaxation_time).exp()
    }

    /// Creep compliance: `J(t) = (1/G) * (1 - exp(-t/lambda)) + t/eta`.
    ///
    /// For the UCM model (no spring in series), the compliance is:
    /// `J(t) = t / eta`.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let eta = self.viscosity();
        if eta < 1e-30 {
            return f64::INFINITY;
        }
        t / eta
    }

    /// Advance the elastic stress component by one explicit Euler step.
    ///
    /// For simple shear the stress evolution equation is:
    /// `d(tau_xy)/dt = eta * gamma - tau_xy / lambda`
    ///
    /// Returns the new `tau_xy`.
    pub fn step_stress(&self, tau_xy: f64, shear_rate: f64, dt: f64) -> f64 {
        if self.relaxation_time < 1e-30 {
            return self.viscosity() * shear_rate;
        }
        let dtau = self.elastic_modulus * shear_rate - tau_xy / self.relaxation_time;
        tau_xy + dt * dtau
    }

    /// Steady-state shear stress: `tau_xy = eta * gamma`.
    pub fn steady_shear_stress(&self, shear_rate: f64) -> f64 {
        self.viscosity() * shear_rate
    }
}

impl NonNewtonianFluid for MaxwellFluid {
    fn effective_viscosity(&self, _shear_rate: f64) -> f64 {
        // Maxwell model: steady-state shear viscosity is constant = eta
        self.viscosity()
    }
}

impl LocalViscosityModel for MaxwellFluid {
    fn viscosity(&self, _gamma_dot: f64) -> f64 {
        MaxwellFluid::viscosity(self)
    }
}

// ---------------------------------------------------------------------------
// Phan-Thien-Tanner (PTT) viscoelastic model
// ---------------------------------------------------------------------------

/// Phan-Thien-Tanner viscoelastic fluid model.
///
/// The PTT model modifies the upper-convected Maxwell (UCM) model with a
/// linear or exponential stress coefficient function:
///
/// f(tau) * tau + lambda * tau^nabla = 2 eta_p D
///
/// where for the linear PTT: f(tau) = 1 + epsilon * lambda * tr(tau) / eta_p
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct PhanThienTanner {
    /// Polymeric viscosity eta_p.
    pub eta_p: f64,
    /// Solvent viscosity eta_s.
    pub eta_s: f64,
    /// Relaxation time lambda.
    pub lambda: f64,
    /// PTT parameter epsilon (extensibility).
    pub epsilon: f64,
    /// Slip parameter xi (0 = upper-convected, 0.5 = corotational).
    pub xi: f64,
}

impl PhanThienTanner {
    /// Create a new PTT fluid with given parameters.
    #[allow(dead_code)]
    pub fn new(eta_p: f64, eta_s: f64, lambda: f64, epsilon: f64, xi: f64) -> Self {
        Self {
            eta_p,
            eta_s,
            lambda,
            epsilon,
            xi,
        }
    }

    /// Total viscosity eta_0 = eta_p + eta_s.
    #[allow(dead_code)]
    pub fn total_viscosity(&self) -> f64 {
        self.eta_p + self.eta_s
    }

    /// Viscosity ratio beta = eta_s / (eta_p + eta_s).
    #[allow(dead_code)]
    pub fn viscosity_ratio(&self) -> f64 {
        self.eta_s / self.total_viscosity()
    }

    /// Linear PTT stress coefficient f(tr(tau)) = 1 + epsilon * lambda * tr(tau) / eta_p.
    #[allow(dead_code)]
    pub fn stress_function_linear(&self, trace_stress: f64) -> f64 {
        1.0 + self.epsilon * self.lambda * trace_stress / self.eta_p
    }

    /// Exponential PTT stress coefficient f(tr(tau)) = exp(epsilon * lambda * tr(tau) / eta_p).
    #[allow(dead_code)]
    pub fn stress_function_exp(&self, trace_stress: f64) -> f64 {
        (self.epsilon * self.lambda * trace_stress / self.eta_p).exp()
    }

    /// Weissenberg number Wi = lambda * |gamma|.
    #[allow(dead_code)]
    pub fn weissenberg_number(&self, shear_rate: f64) -> f64 {
        self.lambda * shear_rate.abs()
    }

    /// Effective shear viscosity at steady simple shear rate gamma (linear PTT).
    ///
    /// eta_eff = eta_p / f(gamma) + eta_s,  f(gamma) = 1 + 2*epsilon*(1-xi)*xi * Wi^2
    #[allow(dead_code)]
    pub fn steady_shear_viscosity(&self, shear_rate: f64) -> f64 {
        let wi = self.weissenberg_number(shear_rate);
        let f = 1.0 + 2.0 * self.epsilon * (1.0 - self.xi) * self.xi * wi * wi;
        self.eta_p / f + self.eta_s
    }

    /// Normal stress difference N1 = tau_xx - tau_yy at steady shear.
    #[allow(dead_code)]
    pub fn normal_stress_difference(&self, shear_rate: f64) -> f64 {
        let wi = self.weissenberg_number(shear_rate);
        let f = 1.0 + 2.0 * self.epsilon * (1.0 - self.xi) * self.xi * wi * wi;
        2.0 * self.eta_p * self.lambda * shear_rate * shear_rate * (1.0 - self.xi) / (f * f)
    }
}

// ---------------------------------------------------------------------------
// Giesekus viscoelastic model
// ---------------------------------------------------------------------------

/// Giesekus viscoelastic model.
///
/// Extends UCM with a quadratic stress term parameterised by alpha:
///
/// tau + lambda * tau^nabla + alpha * lambda / eta_p * tau * tau = 2 eta_p D
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct GiesekusFluid {
    /// Polymeric viscosity.
    pub eta_p: f64,
    /// Relaxation time.
    pub lambda: f64,
    /// Giesekus mobility factor alpha in \[0, 0.5\].
    pub alpha: f64,
}

impl GiesekusFluid {
    /// Create a new Giesekus fluid.
    #[allow(dead_code)]
    pub fn new(eta_p: f64, lambda: f64, alpha: f64) -> Self {
        Self {
            eta_p,
            lambda,
            alpha: alpha.clamp(0.0, 0.5),
        }
    }

    /// Weissenberg number.
    #[allow(dead_code)]
    pub fn weissenberg(&self, shear_rate: f64) -> f64 {
        self.lambda * shear_rate.abs()
    }

    /// Approximate steady shear viscosity eta(Wi).
    ///
    /// For alpha->0 recovers UCM; for finite alpha gives shear thinning.
    #[allow(dead_code)]
    pub fn steady_shear_viscosity(&self, shear_rate: f64) -> f64 {
        let wi = self.weissenberg(shear_rate);
        // Approximate closed-form for Giesekus
        let f = 1.0 + 2.0 * self.alpha * wi * wi;
        self.eta_p / f.sqrt()
    }

    /// Relaxation time for the LBM effective viscosity.
    #[allow(dead_code)]
    pub fn effective_tau(&self, shear_rate: f64, cs2: f64) -> f64 {
        let nu = self.steady_shear_viscosity(shear_rate);
        0.5 + nu / cs2
    }
}

// ---------------------------------------------------------------------------
// RheologyLookupTable -- tabulated viscosity for complex fluids
// ---------------------------------------------------------------------------

/// Pre-tabulated viscosity as a function of shear rate.
///
/// Useful for fluids whose rheology is given by experimental data.
/// Viscosity is linearly interpolated between table points.
#[allow(dead_code)]
pub struct RheologyLookupTable {
    /// Shear rate sample points (must be sorted ascending).
    pub gamma_dot: Vec<f64>,
    /// Corresponding viscosity values.
    pub viscosity: Vec<f64>,
}

impl RheologyLookupTable {
    /// Create a lookup table.
    ///
    /// Panics if `gamma_dot` and `viscosity` have different lengths or are empty.
    #[allow(dead_code)]
    pub fn new(gamma_dot: Vec<f64>, viscosity: Vec<f64>) -> Self {
        assert_eq!(gamma_dot.len(), viscosity.len(), "lengths must match");
        assert!(!gamma_dot.is_empty(), "table must not be empty");
        Self {
            gamma_dot,
            viscosity,
        }
    }

    /// Look up viscosity at the given shear rate via linear interpolation.
    #[allow(dead_code)]
    pub fn viscosity_at(&self, gamma: f64) -> f64 {
        let n = self.gamma_dot.len();
        let g = gamma.abs();
        if g <= self.gamma_dot[0] {
            return self.viscosity[0];
        }
        if g >= self.gamma_dot[n - 1] {
            return self.viscosity[n - 1];
        }
        // Binary search
        let pos = self.gamma_dot.partition_point(|&x| x < g);
        let i = pos - 1;
        let t = (g - self.gamma_dot[i]) / (self.gamma_dot[i + 1] - self.gamma_dot[i]);
        self.viscosity[i] * (1.0 - t) + self.viscosity[i + 1] * t
    }

    /// Compute LBM relaxation time tau from tabulated viscosity.
    #[allow(dead_code)]
    pub fn tau_at(&self, gamma: f64, cs2: f64) -> f64 {
        0.5 + self.viscosity_at(gamma) / cs2
    }
}

// ---------------------------------------------------------------------------
// YieldCriterion -- 2D/3D yield surface checking
// ---------------------------------------------------------------------------

/// 2D and 3D yield criteria for viscoplastic fluids.
#[allow(dead_code)]
pub struct YieldCriterion {
    /// Yield stress tau_y.
    pub tau_yield: f64,
}

impl YieldCriterion {
    /// Create a yield criterion with given yield stress.
    #[allow(dead_code)]
    pub fn new(tau_yield: f64) -> Self {
        Self { tau_yield }
    }

    /// Von Mises yield condition in 2D from stress components (sigma_xx, sigma_yy, sigma_xy).
    ///
    /// Returns `true` if the material has yielded.
    #[allow(dead_code)]
    pub fn is_yielded_2d(&self, s_xx: f64, s_yy: f64, s_xy: f64) -> bool {
        let vm = ((s_xx - s_yy) * (s_xx - s_yy) / 2.0 + s_xy * s_xy).sqrt();
        vm > self.tau_yield
    }

    /// Von Mises equivalent stress in 2D.
    #[allow(dead_code)]
    pub fn von_mises_2d(&self, s_xx: f64, s_yy: f64, s_xy: f64) -> f64 {
        (0.5 * ((s_xx - s_yy) * (s_xx - s_yy) + 2.0 * s_xy * s_xy + (s_xx * s_xx + s_yy * s_yy)))
            .sqrt()
    }

    /// Full 3x3 symmetric deviatoric stress von Mises equivalent.
    #[allow(dead_code)]
    pub fn von_mises_3d(&self, s: [[f64; 3]; 3]) -> f64 {
        let ds = [s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2]];
        let vm2 = 0.5
            * ((ds[0] - ds[1]).powi(2) + (ds[1] - ds[2]).powi(2) + (ds[2] - ds[0]).powi(2))
            + 3.0 * (ds[3] * ds[3] + ds[4] * ds[4] + ds[5] * ds[5]);
        vm2.sqrt()
    }

    /// Bingham correction factor: effective viscosity multiplier.
    ///
    /// Returns `(1 - tau_y / |tau|)` clipped to `[0, 1]`.
    #[allow(dead_code)]
    pub fn bingham_factor(&self, stress_magnitude: f64) -> f64 {
        if stress_magnitude < self.tau_yield {
            0.0
        } else {
            (1.0 - self.tau_yield / stress_magnitude).max(0.0)
        }
    }
}

// ---------------------------------------------------------------------------
// ViscoelasticRelaxation -- simple Maxwell relaxation field
// ---------------------------------------------------------------------------

/// Tracks a spatially varying viscoelastic stress field that relaxes toward zero.
///
/// Useful for LBM coupling with polymer stress: sigma += -sigma/lambda * dt.
#[allow(dead_code)]
pub struct ViscoelasticRelaxation {
    /// Number of cells.
    pub n: usize,
    /// Stress component sigma_xx.
    pub sigma_xx: Vec<f64>,
    /// Stress component sigma_yy.
    pub sigma_yy: Vec<f64>,
    /// Stress component sigma_xy.
    pub sigma_xy: Vec<f64>,
    /// Relaxation time.
    pub lambda: f64,
    /// Polymeric viscosity.
    pub eta_p: f64,
}

impl ViscoelasticRelaxation {
    /// Create a new viscoelastic relaxation field with `n` cells.
    #[allow(dead_code)]
    pub fn new(n: usize, lambda: f64, eta_p: f64) -> Self {
        Self {
            n,
            sigma_xx: vec![0.0; n],
            sigma_yy: vec![0.0; n],
            sigma_xy: vec![0.0; n],
            lambda,
            eta_p,
        }
    }

    /// Apply Maxwell relaxation: sigma_new = sigma * exp(-dt/lambda).
    #[allow(dead_code)]
    pub fn relax(&mut self, dt: f64) {
        let factor = (-dt / self.lambda).exp();
        for i in 0..self.n {
            self.sigma_xx[i] *= factor;
            self.sigma_yy[i] *= factor;
            self.sigma_xy[i] *= factor;
        }
    }

    /// Add elastic contribution from strain rate: sigma += 2 eta_p D * dt.
    #[allow(dead_code)]
    pub fn add_strain_rate(&mut self, dxx: &[f64], dyy: &[f64], dxy: &[f64], dt: f64) {
        for i in 0..self.n {
            self.sigma_xx[i] += 2.0 * self.eta_p * dxx[i] * dt;
            self.sigma_yy[i] += 2.0 * self.eta_p * dyy[i] * dt;
            self.sigma_xy[i] += 2.0 * self.eta_p * dxy[i] * dt;
        }
    }

    /// Compute the trace of the stress tensor at cell `i`.
    #[allow(dead_code)]
    pub fn trace(&self, i: usize) -> f64 {
        self.sigma_xx[i] + self.sigma_yy[i]
    }

    /// Compute |sigma| = sqrt(sigma_xx^2 + sigma_yy^2 + 2 sigma_xy^2) at cell `i`.
    #[allow(dead_code)]
    pub fn magnitude(&self, i: usize) -> f64 {
        (self.sigma_xx[i] * self.sigma_xx[i]
            + self.sigma_yy[i] * self.sigma_yy[i]
            + 2.0 * self.sigma_xy[i] * self.sigma_xy[i])
            .sqrt()
    }

    /// Total stress magnitude summed over all cells.
    #[allow(dead_code)]
    pub fn total_magnitude(&self) -> f64 {
        (0..self.n).map(|i| self.magnitude(i)).sum()
    }
}

// ---------------------------------------------------------------------------
// ShearBandingFluid -- Johnson-Segalman model for shear banding
// ---------------------------------------------------------------------------

/// Simplified Johnson-Segalman model capturing shear banding.
///
/// The model has a non-monotonic flow curve, leading to banding at intermediate
/// stresses. Parameterised by eta_1, eta_2, lambda, and slip parameter xi.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct JohnsonSegalman {
    /// First polymeric viscosity component eta_1.
    pub eta_1: f64,
    /// Solvent viscosity component eta_2.
    pub eta_2: f64,
    /// Relaxation time lambda.
    pub lambda: f64,
    /// Slip parameter xi in \[0, 1\].
    pub xi: f64,
}

impl JohnsonSegalman {
    /// Create a Johnson-Segalman fluid.
    #[allow(dead_code)]
    pub fn new(eta_1: f64, eta_2: f64, lambda: f64, xi: f64) -> Self {
        Self {
            eta_1,
            eta_2,
            lambda,
            xi: xi.clamp(0.0, 1.0),
        }
    }

    /// Total viscosity.
    #[allow(dead_code)]
    pub fn total_viscosity(&self) -> f64 {
        self.eta_1 + self.eta_2
    }

    /// Steady-state shear stress at given shear rate (approximate).
    #[allow(dead_code)]
    pub fn shear_stress(&self, gamma_dot: f64) -> f64 {
        let wi = self.lambda * gamma_dot.abs();
        let f = 1.0 + (1.0 - self.xi) * self.xi * wi * wi;
        self.eta_1 * gamma_dot / f + self.eta_2 * gamma_dot
    }

    /// Compute apparent viscosity.
    #[allow(dead_code)]
    pub fn apparent_viscosity(&self, gamma_dot: f64) -> f64 {
        if gamma_dot.abs() < 1e-30 {
            return self.total_viscosity();
        }
        self.shear_stress(gamma_dot) / gamma_dot
    }
}

// ---------------------------------------------------------------------------
// Wormlike Micelle model (Rolie-Poly)
// ---------------------------------------------------------------------------

/// Simplified Rolie-Poly model for wormlike micelles.
///
/// dA/dt = kappa * A + A * kappa^T - (1/tau_d) * (A - I) - (1/tau_R) * (tr(A)/3 - 1) * A
///
/// where A is the conformation tensor, tau_d the disengagement time and tau_R
/// the Rouse time.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct RoliePolyModel {
    /// Polymeric viscosity eta_p.
    pub eta_p: f64,
    /// Reptation (disengagement) time tau_d.
    pub tau_d: f64,
    /// Rouse time tau_R.
    pub tau_r: f64,
    /// CCR parameter beta_ccr.
    pub beta_ccr: f64,
}

impl RoliePolyModel {
    /// Create a Rolie-Poly model.
    #[allow(dead_code)]
    pub fn new(eta_p: f64, tau_d: f64, tau_r: f64, beta_ccr: f64) -> Self {
        Self {
            eta_p,
            tau_d,
            tau_r,
            beta_ccr,
        }
    }

    /// Steady-state effective viscosity approximation.
    #[allow(dead_code)]
    pub fn steady_viscosity(&self, shear_rate: f64) -> f64 {
        let wi_d = self.tau_d * shear_rate.abs();
        let factor = 1.0 + 3.0 * wi_d * wi_d;
        self.eta_p / factor
    }

    /// Zero-shear viscosity = eta_p.
    #[allow(dead_code)]
    pub fn zero_shear_viscosity(&self) -> f64 {
        self.eta_p
    }

    /// Weissenberg number for disengagement time.
    #[allow(dead_code)]
    pub fn weissenberg_d(&self, shear_rate: f64) -> f64 {
        self.tau_d * shear_rate.abs()
    }
}
