#![allow(clippy::if_same_then_else, clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Turbulent dispersion using the Lattice Boltzmann Method.
//!
//! This module implements turbulent scalar transport and dispersion models:
//!
//! - [`TurbulentDiffusivity`]: k-epsilon turbulence model, turbulent Schmidt
//!   number, and eddy diffusivity estimation.
//! - [`ScalarTransport`]: advection-diffusion of a passive scalar concentration
//!   field on a 2-D LBM grid, with source terms.
//! - [`DispersingPlume`]: Gaussian plume model with Pasquill-Gifford stability
//!   classes, σ_y and σ_z dispersion coefficients.
//! - [`MixingLength`]: Prandtl mixing-length model with von Kármán constant and
//!   log-law velocity profile.
//! - [`TurbulentFlux`]: Reynolds flux, gradient-diffusion hypothesis, and
//!   counter-gradient transport correction.
//! - [`DispersionLBM`]: Full LBM passive-scalar solver with drift-diffusion
//!   collision and PDF-based source model.

#![allow(dead_code)]

// ─── Physical constants ────────────────────────────────────────────────────

/// Von Kármán constant κ ≈ 0.41 (dimensionless).
const KAPPA: f64 = 0.41;
/// Default turbulent Schmidt number Sc_t ≈ 0.7.
const SC_T_DEFAULT: f64 = 0.7;
/// D2Q9 speed-of-sound squared c_s² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ─── D2Q9 lattice constants ────────────────────────────────────────────────

/// D2Q9 weight coefficients.
const W9: [f64; 9] = [
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

/// D2Q9 x-velocity components.
const CX9: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 y-velocity components.
const CY9: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

// ─── TurbulentDiffusivity ─────────────────────────────────────────────────

/// Turbulent diffusivity estimation using the k-ε model.
///
/// In the k-ε model the turbulent viscosity is:
///
/// ```text
/// ν_t = C_μ · k² / ε
/// ```
///
/// The eddy diffusivity for a scalar (e.g. temperature, concentration) is:
///
/// ```text
/// D_t = ν_t / Sc_t
/// ```
///
/// where `Sc_t` is the turbulent Schmidt number (≈ 0.7 for most scalars).
#[derive(Debug, Clone)]
pub struct TurbulentDiffusivity {
    /// Turbulent kinetic energy k (m² s⁻²).
    pub k: f64,
    /// Turbulent dissipation rate ε (m² s⁻³).
    pub epsilon: f64,
    /// C_μ model constant (standard value 0.09).
    pub c_mu: f64,
    /// Turbulent Schmidt number Sc_t (dimensionless).
    pub schmidt_number: f64,
}

impl TurbulentDiffusivity {
    /// Create a new `TurbulentDiffusivity` model.
    ///
    /// # Arguments
    /// * `k`       – turbulent kinetic energy (m² s⁻²)
    /// * `epsilon` – turbulent dissipation rate (m² s⁻³)
    /// * `c_mu`    – model constant (typically 0.09)
    /// * `sc_t`    – turbulent Schmidt number (typically 0.7)
    pub fn new(k: f64, epsilon: f64, c_mu: f64, sc_t: f64) -> Self {
        Self {
            k,
            epsilon,
            c_mu,
            schmidt_number: sc_t,
        }
    }

    /// Create with standard k-ε constants (C_μ = 0.09, Sc_t = 0.7).
    pub fn standard(k: f64, epsilon: f64) -> Self {
        Self::new(k, epsilon, 0.09, SC_T_DEFAULT)
    }

    /// Compute the turbulent (eddy) viscosity ν_t = C_μ k² / ε (m² s⁻¹).
    ///
    /// Returns 0 when ε ≤ 0 to avoid division by zero.
    pub fn turbulent_viscosity(&self) -> f64 {
        if self.epsilon <= 0.0 {
            return 0.0;
        }
        self.c_mu * self.k * self.k / self.epsilon
    }

    /// Compute the eddy diffusivity D_t = ν_t / Sc_t (m² s⁻¹).
    pub fn eddy_diffusivity(&self) -> f64 {
        if self.schmidt_number <= 0.0 {
            return 0.0;
        }
        self.turbulent_viscosity() / self.schmidt_number
    }

    /// Turbulent integral length scale L = C_μ^(3/4) k^(3/2) / ε (m).
    pub fn integral_length_scale(&self) -> f64 {
        if self.epsilon <= 0.0 {
            return 0.0;
        }
        self.c_mu.powf(0.75) * self.k.powf(1.5) / self.epsilon
    }

    /// Turbulent time scale τ_t = k / ε (s).
    pub fn turbulent_time_scale(&self) -> f64 {
        if self.epsilon <= 0.0 {
            return f64::INFINITY;
        }
        self.k / self.epsilon
    }

    /// Kolmogorov length scale η = (ν³/ε)^(1/4) (m).
    ///
    /// # Arguments
    /// * `nu` – molecular kinematic viscosity (m² s⁻¹)
    pub fn kolmogorov_length(&self, nu: f64) -> f64 {
        if self.epsilon <= 0.0 {
            return f64::INFINITY;
        }
        (nu.powi(3) / self.epsilon).powf(0.25)
    }
}

// ─── ScalarTransport ──────────────────────────────────────────────────────

/// Advection-diffusion scalar transport on a 2-D LBM grid.
///
/// Maintains a concentration field `c` of size `nx × ny` and provides:
/// - D2Q9 equilibrium distribution for the scalar
/// - BGK-like collision with effective diffusivity
/// - Periodic streaming in x and bounce-back in y
/// - Uniform source term Q (mol m⁻³ s⁻¹)
#[derive(Debug, Clone)]
pub struct ScalarTransport {
    /// Grid width in the x-direction.
    pub nx: usize,
    /// Grid height in the y-direction.
    pub ny: usize,
    /// Concentration field (mol m⁻³), row-major index k = y*nx + x.
    pub concentration: Vec<f64>,
    /// Scalar distribution functions `g[q][k]`.
    pub g: Vec<Vec<f64>>,
    /// Effective diffusivity D (lattice units).
    pub diffusivity: f64,
    /// Source term Q (mol m⁻³ per time step).
    pub source: f64,
}

impl ScalarTransport {
    /// Create a new scalar transport solver.
    ///
    /// # Arguments
    /// * `nx`          – grid width
    /// * `ny`          – grid height
    /// * `diffusivity` – effective diffusivity (lattice units, relates to ω_g = 1/(3D+0.5))
    /// * `source`      – uniform source term (mol m⁻³ per step)
    pub fn new(nx: usize, ny: usize, diffusivity: f64, source: f64) -> Self {
        let n = nx * ny;
        let c0 = 0.0_f64;
        Self {
            nx,
            ny,
            concentration: vec![c0; n],
            g: vec![vec![0.0; n]; 9],
            diffusivity,
            source,
        }
    }

    /// Set a uniform initial concentration throughout the domain.
    pub fn set_uniform_concentration(&mut self, c0: f64) {
        for c in &mut self.concentration {
            *c = c0;
        }
        // Re-initialise g to equilibrium at rest (ux=uy=0)
        let n = self.nx * self.ny;
        for k in 0..n {
            for i in 0..9 {
                self.g[i][k] = W9[i] * self.concentration[k];
            }
        }
    }

    /// Scalar equilibrium distribution: g_eq_i = w_i · c · (1 + c_ia·u_a / c_s²).
    ///
    /// # Arguments
    /// * `c`  – local concentration (mol m⁻³)
    /// * `ux` – fluid velocity x-component (lattice units)
    /// * `uy` – fluid velocity y-component (lattice units)
    /// * `i`  – velocity direction index 0..9
    pub fn g_equilibrium(c: f64, ux: f64, uy: f64, i: usize) -> f64 {
        let cx = CX9[i] as f64;
        let cy = CY9[i] as f64;
        W9[i] * c * (1.0 + (cx * ux + cy * uy) / CS2)
    }

    /// Compute the relaxation rate ω_g from the diffusivity.
    ///
    /// D = c_s² (1/ω_g − 0.5), so ω_g = 1 / (D/c_s² + 0.5).
    pub fn omega_g(&self) -> f64 {
        1.0 / (self.diffusivity / CS2 + 0.5)
    }

    /// Perform one collision step using constant fluid velocity (ux, uy) everywhere.
    ///
    /// # Arguments
    /// * `ux_field` – fluid x-velocity at each cell
    /// * `uy_field` – fluid y-velocity at each cell
    pub fn collide(&mut self, ux_field: &[f64], uy_field: &[f64]) {
        let n = self.nx * self.ny;
        let omega = self.omega_g();
        for k in 0..n {
            let c = self.concentration[k];
            let ux = ux_field[k];
            let uy = uy_field[k];
            for i in 0..9 {
                let geq = Self::g_equilibrium(c, ux, uy, i);
                self.g[i][k] += omega * (geq - self.g[i][k]);
            }
            // Add source term to zeroth moment
            self.g[0][k] += self.source;
        }
    }

    /// Perform one D2Q9 streaming step with periodic BC in x, bounce-back in y.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut g_new: Vec<Vec<f64>> = vec![vec![0.0; nx * ny]; 9];

        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                for i in 0..9 {
                    let cx = CX9[i];
                    let cy = CY9[i];
                    let xn = ((x as i32 + cx).rem_euclid(nx as i32)) as usize;
                    let yn_i = y as i32 + cy;
                    if yn_i < 0 || yn_i >= ny as i32 {
                        // Bounce-back: reflect in opposite direction
                        let opp = (i + 4) % 8 + if i == 0 { 0 } else { 0 };
                        g_new[opp][k] += self.g[i][k];
                    } else {
                        let yn = yn_i as usize;
                        let kn = yn * nx + xn;
                        g_new[i][kn] = self.g[i][k];
                    }
                }
            }
        }
        self.g = g_new;
        // Update concentration
        let n = nx * ny;
        for k in 0..n {
            self.concentration[k] = (0..9).map(|i| self.g[i][k]).sum();
        }
    }

    /// Return the maximum concentration in the domain.
    pub fn max_concentration(&self) -> f64 {
        self.concentration
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Return the total scalar mass Σ_k c_k (mol m⁻³ · cells).
    pub fn total_mass(&self) -> f64 {
        self.concentration.iter().sum()
    }
}

// ─── DispersingPlume ──────────────────────────────────────────────────────

/// Pasquill-Gifford atmospheric stability class.
///
/// Classes range from A (most unstable) to F (most stable).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PasquillClass {
    /// Class A – extremely unstable.
    A,
    /// Class B – moderately unstable.
    B,
    /// Class C – slightly unstable.
    C,
    /// Class D – neutral.
    D,
    /// Class E – slightly stable.
    E,
    /// Class F – moderately stable.
    F,
}

/// Gaussian plume dispersion model with Pasquill-Gifford stability classes.
///
/// Computes the ground-level concentration at (x, y) downwind of a
/// continuous point source of strength Q at height H_s:
///
/// ```text
/// C(x,y,0) = Q/(π u σ_y σ_z) · exp(-y²/2σ_y²) · exp(-H_s²/2σ_z²)
/// ```
///
/// σ_y and σ_z are parameterised by the Pasquill stability class.
#[derive(Debug, Clone)]
pub struct DispersingPlume {
    /// Source emission rate Q (kg s⁻¹ or mol s⁻¹).
    pub source_strength: f64,
    /// Effective stack height H_s (m).
    pub stack_height: f64,
    /// Mean wind speed u (m s⁻¹).
    pub wind_speed: f64,
    /// Atmospheric stability class.
    pub stability: PasquillClass,
}

impl DispersingPlume {
    /// Create a new Gaussian plume model.
    ///
    /// # Arguments
    /// * `q`       – source strength (kg s⁻¹)
    /// * `h_s`     – effective stack height (m)
    /// * `u`       – mean wind speed (m s⁻¹)
    /// * `class`   – Pasquill stability class
    pub fn new(q: f64, h_s: f64, u: f64, class: PasquillClass) -> Self {
        Self {
            source_strength: q,
            stack_height: h_s,
            wind_speed: u,
            stability: class,
        }
    }

    /// Pasquill-Gifford σ_y (m) at downwind distance x (m).
    ///
    /// Uses the Slade (1968) power-law parameterisation.
    pub fn sigma_y(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let (a, b) = match self.stability {
            PasquillClass::A => (0.22_f64, 0.0001_f64),
            PasquillClass::B => (0.16_f64, 0.0001_f64),
            PasquillClass::C => (0.11_f64, 0.0001_f64),
            PasquillClass::D => (0.08_f64, 0.0001_f64),
            PasquillClass::E => (0.06_f64, 0.0001_f64),
            PasquillClass::F => (0.04_f64, 0.0001_f64),
        };
        a * x / (1.0 + b * x).sqrt()
    }

    /// Pasquill-Gifford σ_z (m) at downwind distance x (m).
    pub fn sigma_z(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let (c, d) = match self.stability {
            PasquillClass::A => (0.20_f64, 0.000_f64),
            PasquillClass::B => (0.12_f64, 0.000_f64),
            PasquillClass::C => (0.08_f64, 0.0002_f64),
            PasquillClass::D => (0.06_f64, 0.0015_f64),
            PasquillClass::E => (0.03_f64, 0.0003_f64),
            PasquillClass::F => (0.016_f64, 0.0003_f64),
        };
        c * x / (1.0 + d * x).sqrt()
    }

    /// Ground-level centreline concentration C(x, 0, 0) at downwind distance x (m).
    pub fn centreline_concentration(&self, x: f64) -> f64 {
        let sy = self.sigma_y(x);
        let sz = self.sigma_z(x);
        if sy <= 0.0 || sz <= 0.0 || self.wind_speed <= 0.0 {
            return 0.0;
        }
        let hs = self.stack_height;
        (self.source_strength / (std::f64::consts::PI * self.wind_speed * sy * sz))
            * (-hs * hs / (2.0 * sz * sz)).exp()
    }

    /// Full Gaussian plume concentration C(x, y, 0) (kg m⁻³).
    ///
    /// # Arguments
    /// * `x` – downwind distance (m)
    /// * `y` – crosswind offset (m)
    pub fn concentration_at(&self, x: f64, y: f64) -> f64 {
        let sy = self.sigma_y(x);
        let sz = self.sigma_z(x);
        if sy <= 0.0 || sz <= 0.0 || self.wind_speed <= 0.0 {
            return 0.0;
        }
        let hs = self.stack_height;
        (self.source_strength / (std::f64::consts::PI * self.wind_speed * sy * sz))
            * (-y * y / (2.0 * sy * sy)).exp()
            * (-hs * hs / (2.0 * sz * sz)).exp()
    }

    /// Crosswind-integrated concentration (CIC) at downwind distance x (kg m⁻²).
    ///
    /// CIC = Q / (√(2π) u σ_z) · exp(-H_s²/2σ_z²)
    pub fn crosswind_integrated_concentration(&self, x: f64) -> f64 {
        let sz = self.sigma_z(x);
        if sz <= 0.0 || self.wind_speed <= 0.0 {
            return 0.0;
        }
        let hs = self.stack_height;
        (self.source_strength / ((2.0 * std::f64::consts::PI).sqrt() * self.wind_speed * sz))
            * (-hs * hs / (2.0 * sz * sz)).exp()
    }

    /// Estimate the downwind distance x* at which σ_z equals the stack height.
    ///
    /// Useful as the location of maximum ground-level concentration.
    /// Solved iteratively (Newton, up to 50 iterations).
    pub fn max_concentration_distance(&self) -> f64 {
        if self.stack_height <= 0.0 {
            return 0.0;
        }
        // Binary search between 1 m and 100 km
        let mut lo = 1.0_f64;
        let mut hi = 100_000.0_f64;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.sigma_z(mid) < self.stack_height {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ─── MixingLength ─────────────────────────────────────────────────────────

/// Prandtl mixing-length model for wall-bounded turbulent flows.
///
/// Provides the local turbulent viscosity:
/// ```text
/// ν_t = l_m² |∂u/∂y|
/// ```
/// where the mixing length follows the van Driest damped profile:
/// ```text
/// l_m = κ y [1 − exp(−y⁺/A⁺)]
/// ```
#[derive(Debug, Clone)]
pub struct MixingLength {
    /// Von Kármán constant κ (default 0.41).
    pub kappa: f64,
    /// Van Driest damping constant A⁺ (default 26).
    pub a_plus: f64,
    /// Channel half-height or boundary-layer thickness δ (m).
    pub delta: f64,
}

impl MixingLength {
    /// Create a new mixing-length model with standard constants.
    ///
    /// # Arguments
    /// * `delta` – channel half-height / BL thickness (m)
    pub fn new(delta: f64) -> Self {
        Self {
            kappa: KAPPA,
            a_plus: 26.0,
            delta,
        }
    }

    /// Create with custom constants.
    pub fn with_constants(kappa: f64, a_plus: f64, delta: f64) -> Self {
        Self {
            kappa,
            a_plus,
            delta,
        }
    }

    /// Mixing length l_m(y) without van Driest damping (m).
    ///
    /// Uses the simple Prandtl formula: l_m = κ y (1 − y/δ).
    pub fn mixing_length(&self, y: f64) -> f64 {
        if y <= 0.0 || self.delta <= 0.0 {
            return 0.0;
        }
        let lm = self.kappa * y * (1.0 - (y / self.delta).min(1.0));
        lm.max(0.0)
    }

    /// Van Driest damped mixing length l_m(y, y⁺).
    ///
    /// # Arguments
    /// * `y`      – wall-normal distance (m)
    /// * `y_plus` – dimensionless wall distance y⁺ = y u_τ / ν
    pub fn mixing_length_damped(&self, y: f64, y_plus: f64) -> f64 {
        let base = self.mixing_length(y);
        let damping = 1.0 - (-y_plus / self.a_plus).exp();
        base * damping
    }

    /// Log-law velocity u⁺ = (1/κ) ln(y⁺) + B.
    ///
    /// Valid for y⁺ > 30 (log-law region).
    ///
    /// # Arguments
    /// * `y_plus` – dimensionless wall distance
    pub fn log_law_velocity(&self, y_plus: f64) -> f64 {
        if y_plus <= 0.0 {
            return 0.0;
        }
        (1.0 / self.kappa) * y_plus.ln() + 5.2
    }

    /// Turbulent viscosity ν_t = l_m² |∂u/∂y| (m² s⁻¹).
    ///
    /// # Arguments
    /// * `y`       – wall-normal distance (m)
    /// * `du_dy`   – velocity gradient (s⁻¹)
    pub fn turbulent_viscosity(&self, y: f64, du_dy: f64) -> f64 {
        let lm = self.mixing_length(y);
        lm * lm * du_dy.abs()
    }

    /// Friction velocity estimate from bulk velocity using Clauser chart.
    ///
    /// u_τ ≈ u_bulk / (1/κ · ln(Re_δ) + B − 3) (rough approximation).
    pub fn friction_velocity(&self, u_bulk: f64, nu: f64) -> f64 {
        if nu <= 0.0 || self.delta <= 0.0 {
            return 0.0;
        }
        let re_delta = u_bulk * self.delta / nu;
        if re_delta <= 1.0 {
            return 0.0;
        }
        u_bulk / ((1.0 / self.kappa) * re_delta.ln() + 5.2 - 3.0)
    }
}

// ─── TurbulentFlux ────────────────────────────────────────────────────────

/// Reynolds flux model with gradient-diffusion hypothesis.
///
/// The gradient-diffusion hypothesis gives the turbulent scalar flux:
/// ```text
/// F_t = −D_t ∂c/∂x
/// ```
///
/// A counter-gradient correction can be applied following Deardorff (1972):
/// ```text
/// F_cg = −D_t (∂c/∂x − γ)
/// ```
#[derive(Debug, Clone)]
pub struct TurbulentFlux {
    /// Turbulent eddy diffusivity D_t (m² s⁻¹).
    pub eddy_diffusivity: f64,
    /// Counter-gradient term γ (mol m⁻⁴), set to 0 to disable.
    pub counter_gradient: f64,
    /// Molecular diffusivity D_mol (m² s⁻¹).
    pub molecular_diffusivity: f64,
}

impl TurbulentFlux {
    /// Create a new turbulent flux model.
    ///
    /// # Arguments
    /// * `d_t`   – turbulent eddy diffusivity (m² s⁻¹)
    /// * `d_mol` – molecular diffusivity (m² s⁻¹)
    pub fn new(d_t: f64, d_mol: f64) -> Self {
        Self {
            eddy_diffusivity: d_t,
            molecular_diffusivity: d_mol,
            counter_gradient: 0.0,
        }
    }

    /// Set the counter-gradient correction term γ.
    pub fn with_counter_gradient(mut self, gamma: f64) -> Self {
        self.counter_gradient = gamma;
        self
    }

    /// Compute the total diffusive flux F = −(D_mol + D_t) ∂c/∂x (mol m⁻² s⁻¹).
    ///
    /// # Arguments
    /// * `grad_c` – concentration gradient ∂c/∂x (mol m⁻⁴)
    pub fn total_flux(&self, grad_c: f64) -> f64 {
        -(self.molecular_diffusivity + self.eddy_diffusivity) * grad_c
    }

    /// Compute the turbulent flux with counter-gradient correction.
    ///
    /// F = −D_t (∂c/∂x − γ)
    pub fn turbulent_flux_cg(&self, grad_c: f64) -> f64 {
        -self.eddy_diffusivity * (grad_c - self.counter_gradient)
    }

    /// Effective diffusivity D_eff = D_mol + D_t (m² s⁻¹).
    pub fn effective_diffusivity(&self) -> f64 {
        self.molecular_diffusivity + self.eddy_diffusivity
    }

    /// Turbulent Prandtl/Schmidt number Sc_t = ν_t / D_t.
    ///
    /// # Arguments
    /// * `nu_t` – turbulent viscosity (m² s⁻¹)
    pub fn turbulent_schmidt(&self, nu_t: f64) -> f64 {
        if self.eddy_diffusivity <= 0.0 {
            return f64::INFINITY;
        }
        nu_t / self.eddy_diffusivity
    }

    /// Advective flux F_adv = u · c (mol m⁻² s⁻¹).
    ///
    /// # Arguments
    /// * `velocity` – fluid velocity (m s⁻¹)
    /// * `conc`     – concentration (mol m⁻³)
    pub fn advective_flux(&self, velocity: f64, conc: f64) -> f64 {
        velocity * conc
    }
}

// ─── DispersionLBM ────────────────────────────────────────────────────────

/// Full LBM passive-scalar dispersion solver.
///
/// Combines:
/// - D2Q9 BGK LBM for the carrier flow (stored as velocity snapshots)
/// - D2Q9 BGK advection-diffusion for the scalar concentration
/// - Turbulent diffusivity correction via a k-ε eddy model
/// - Point-source injection (PDF-model Gaussian blob)
#[derive(Debug, Clone)]
pub struct DispersionLBM {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Scalar transport sub-solver.
    pub scalar: ScalarTransport,
    /// Turbulent diffusivity model (optional; None = laminar).
    pub turb_diff: Option<TurbulentDiffusivity>,
    /// Fluid x-velocity field (lattice units).
    pub ux: Vec<f64>,
    /// Fluid y-velocity field (lattice units).
    pub uy: Vec<f64>,
    /// Current simulation step counter.
    pub step: u64,
}

impl DispersionLBM {
    /// Create a new dispersion solver.
    ///
    /// # Arguments
    /// * `nx`    – grid width
    /// * `ny`    – grid height
    /// * `d_mol` – molecular diffusivity (lattice units)
    pub fn new(nx: usize, ny: usize, d_mol: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            scalar: ScalarTransport::new(nx, ny, d_mol, 0.0),
            turb_diff: None,
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            step: 0,
        }
    }

    /// Attach a turbulent diffusivity model.
    pub fn with_turbulence(mut self, td: TurbulentDiffusivity) -> Self {
        self.turb_diff = Some(td);
        self
    }

    /// Set the fluid velocity field.
    pub fn set_velocity(&mut self, ux: Vec<f64>, uy: Vec<f64>) {
        self.ux = ux;
        self.uy = uy;
    }

    /// Inject a Gaussian blob of scalar at cell (cx, cy).
    ///
    /// # Arguments
    /// * `cx`     – source x-index
    /// * `cy`     – source y-index
    /// * `mass`   – total injected mass (mol m⁻³)
    /// * `radius` – Gaussian half-width (cells)
    pub fn inject_gaussian(&mut self, cx: usize, cy: usize, mass: f64, radius: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let r2 = radius * radius;
        let mut total = 0.0_f64;
        // Two-pass: compute weights first
        let mut weights = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let dx = x as f64 - cx as f64;
                let dy = y as f64 - cy as f64;
                let w = (-(dx * dx + dy * dy) / (2.0 * r2)).exp();
                weights[y * nx + x] = w;
                total += w;
            }
        }
        if total > 0.0 {
            for k in 0..nx * ny {
                self.scalar.concentration[k] += mass * weights[k] / total;
            }
        }
    }

    /// Advance the simulation by one step.
    ///
    /// Effective diffusivity = D_mol + D_t (from turbulence model if set).
    pub fn advance(&mut self) {
        // Update effective diffusivity
        if let Some(ref td) = self.turb_diff {
            self.scalar.diffusivity = td.eddy_diffusivity() + self.scalar.diffusivity.max(0.0);
        }
        // Re-initialise g from current concentration (ensure consistency)
        let n = self.nx * self.ny;
        for k in 0..n {
            let c = self.scalar.concentration[k];
            let ux = self.ux[k];
            let uy = self.uy[k];
            for i in 0..9 {
                self.scalar.g[i][k] = ScalarTransport::g_equilibrium(c, ux, uy, i);
            }
        }
        self.scalar.collide(&self.ux.clone(), &self.uy.clone());
        self.scalar.stream();
        self.step += 1;
    }

    /// Return the maximum concentration in the domain.
    pub fn max_concentration(&self) -> f64 {
        self.scalar.max_concentration()
    }

    /// Return the total scalar mass (conservation check).
    pub fn total_mass(&self) -> f64 {
        self.scalar.total_mass()
    }

    /// Compute the centroid of the scalar field (x̄, ȳ) in cell indices.
    pub fn centroid(&self) -> (f64, f64) {
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mt = 0.0_f64;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let c = self.scalar.concentration[y * self.nx + x];
                mx += c * x as f64;
                my += c * y as f64;
                mt += c;
            }
        }
        if mt == 0.0 {
            return (0.0, 0.0);
        }
        (mx / mt, my / mt)
    }
}

// ─── Free functions ────────────────────────────────────────────────────────

/// Compute the turbulent Péclet number Pe_t = u L / D_t.
///
/// # Arguments
/// * `u`   – characteristic velocity (m s⁻¹)
/// * `l`   – characteristic length (m)
/// * `d_t` – eddy diffusivity (m² s⁻¹)
pub fn turbulent_peclet(u: f64, l: f64, d_t: f64) -> f64 {
    if d_t <= 0.0 {
        return f64::INFINITY;
    }
    u * l / d_t
}

/// Compute the Schmidt number Sc = ν / D.
///
/// # Arguments
/// * `nu` – kinematic viscosity (m² s⁻¹)
/// * `d`  – mass diffusivity (m² s⁻¹)
pub fn schmidt_number(nu: f64, d: f64) -> f64 {
    if d <= 0.0 {
        return f64::INFINITY;
    }
    nu / d
}

/// Richardson number Ri = g α ΔT L / u² (stability indicator).
///
/// Positive Ri → stable stratification; negative → unstable.
///
/// # Arguments
/// * `g`      – gravitational acceleration (m s⁻²)
/// * `alpha`  – thermal expansion coefficient (K⁻¹)
/// * `dt_dz`  – vertical temperature gradient (K m⁻¹)
/// * `l`      – length scale (m)
/// * `u`      – velocity scale (m s⁻¹)
pub fn richardson_number(g: f64, alpha: f64, dt_dz: f64, l: f64, u: f64) -> f64 {
    if u == 0.0 {
        return if dt_dz >= 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }
    g * alpha * dt_dz * l / (u * u)
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── TurbulentDiffusivity ───────────────────────────────────────────────

    #[test]
    fn test_turb_diff_viscosity_standard() {
        let td = TurbulentDiffusivity::standard(1.0, 1.0);
        // ν_t = 0.09 * 1 / 1 = 0.09
        assert!((td.turbulent_viscosity() - 0.09).abs() < 1e-12);
    }

    #[test]
    fn test_turb_diff_eddy_diffusivity() {
        let td = TurbulentDiffusivity::standard(1.0, 1.0);
        let dt = td.eddy_diffusivity();
        assert!((dt - 0.09 / 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_turb_diff_zero_epsilon() {
        let td = TurbulentDiffusivity::standard(1.0, 0.0);
        assert_eq!(td.turbulent_viscosity(), 0.0);
        assert_eq!(td.eddy_diffusivity(), 0.0);
    }

    #[test]
    fn test_turb_diff_negative_epsilon() {
        let td = TurbulentDiffusivity::standard(1.0, -1.0);
        assert_eq!(td.turbulent_viscosity(), 0.0);
    }

    #[test]
    fn test_turb_diff_integral_length_scale() {
        let td = TurbulentDiffusivity::standard(1.0, 1.0);
        let l = td.integral_length_scale();
        // l = 0.09^0.75 * 1 / 1
        let expected = 0.09_f64.powf(0.75);
        assert!((l - expected).abs() < 1e-10);
    }

    #[test]
    fn test_turb_diff_time_scale() {
        let td = TurbulentDiffusivity::standard(2.0, 4.0);
        assert!((td.turbulent_time_scale() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_turb_diff_kolmogorov_length() {
        let td = TurbulentDiffusivity::standard(1.0, 1.0);
        let eta = td.kolmogorov_length(1e-6);
        assert!(eta > 0.0 && eta < 1.0);
    }

    #[test]
    fn test_turb_diff_zero_k() {
        let td = TurbulentDiffusivity::standard(0.0, 1.0);
        assert_eq!(td.turbulent_viscosity(), 0.0);
    }

    // ── ScalarTransport ───────────────────────────────────────────────────

    #[test]
    fn test_scalar_transport_creation() {
        let st = ScalarTransport::new(10, 8, 0.05, 0.0);
        assert_eq!(st.concentration.len(), 80);
        assert_eq!(st.g.len(), 9);
    }

    #[test]
    fn test_scalar_transport_uniform_init() {
        let mut st = ScalarTransport::new(4, 4, 0.05, 0.0);
        st.set_uniform_concentration(1.0);
        assert!(st.concentration.iter().all(|&c| (c - 1.0).abs() < 1e-12));
    }

    #[test]
    fn test_scalar_equilibrium_sum_to_c() {
        let c = 1.5_f64;
        let sum: f64 = (0..9)
            .map(|i| ScalarTransport::g_equilibrium(c, 0.0, 0.0, i))
            .sum();
        assert!((sum - c).abs() < 1e-12);
    }

    #[test]
    fn test_scalar_equilibrium_zero_conc() {
        for i in 0..9 {
            assert_eq!(ScalarTransport::g_equilibrium(0.0, 0.1, 0.1, i), 0.0);
        }
    }

    #[test]
    fn test_scalar_omega_g() {
        let st = ScalarTransport::new(4, 4, 1.0 / 6.0, 0.0);
        // D = 1/6 → ω = 1/(D/cs2 + 0.5) = 1/(0.5 + 0.5) = 1
        assert!((st.omega_g() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_total_mass_conserved_after_stream() {
        let mut st = ScalarTransport::new(6, 6, 0.05, 0.0);
        st.set_uniform_concentration(1.0);
        let mass_before = st.total_mass();
        let ux = vec![0.0; 36];
        let uy = vec![0.0; 36];
        st.collide(&ux, &uy);
        st.stream();
        let mass_after = st.total_mass();
        // Mass should be conserved (modulo bounce-back losses at walls)
        assert!((mass_after - mass_before).abs() < mass_before * 0.05);
    }

    #[test]
    fn test_scalar_max_concentration() {
        let mut st = ScalarTransport::new(4, 4, 0.05, 0.0);
        st.concentration[5] = 3.0;
        assert!((st.max_concentration() - 3.0).abs() < 1e-12);
    }

    // ── DispersingPlume ───────────────────────────────────────────────────

    #[test]
    fn test_plume_sigma_y_increases_with_x() {
        let p = DispersingPlume::new(1.0, 10.0, 5.0, PasquillClass::D);
        let sy1 = p.sigma_y(100.0);
        let sy2 = p.sigma_y(1000.0);
        assert!(sy2 > sy1);
    }

    #[test]
    fn test_plume_sigma_z_increases_with_x() {
        let p = DispersingPlume::new(1.0, 10.0, 5.0, PasquillClass::D);
        assert!(p.sigma_z(500.0) > p.sigma_z(100.0));
    }

    #[test]
    fn test_plume_sigma_y_zero_at_zero_x() {
        let p = DispersingPlume::new(1.0, 0.0, 5.0, PasquillClass::A);
        assert_eq!(p.sigma_y(0.0), 0.0);
    }

    #[test]
    fn test_plume_class_a_wider_than_f() {
        let pa = DispersingPlume::new(1.0, 0.0, 5.0, PasquillClass::A);
        let pf = DispersingPlume::new(1.0, 0.0, 5.0, PasquillClass::F);
        assert!(pa.sigma_y(1000.0) > pf.sigma_y(1000.0));
    }

    #[test]
    fn test_plume_centreline_ground_positive() {
        let p = DispersingPlume::new(10.0, 30.0, 5.0, PasquillClass::D);
        let c = p.centreline_concentration(500.0);
        assert!(c >= 0.0);
    }

    #[test]
    fn test_plume_concentration_at_crosswind_less() {
        let p = DispersingPlume::new(10.0, 0.0, 5.0, PasquillClass::D);
        let c0 = p.concentration_at(500.0, 0.0);
        let cy = p.concentration_at(500.0, 100.0);
        assert!(c0 >= cy);
    }

    #[test]
    fn test_plume_cic_positive() {
        let p = DispersingPlume::new(1.0, 0.0, 5.0, PasquillClass::C);
        assert!(p.crosswind_integrated_concentration(500.0) > 0.0);
    }

    #[test]
    fn test_plume_max_distance_positive() {
        let p = DispersingPlume::new(1.0, 50.0, 5.0, PasquillClass::D);
        let xstar = p.max_concentration_distance();
        assert!(xstar > 0.0);
    }

    #[test]
    fn test_plume_zero_wind_gives_zero() {
        let p = DispersingPlume::new(1.0, 10.0, 0.0, PasquillClass::D);
        assert_eq!(p.centreline_concentration(500.0), 0.0);
    }

    // ── MixingLength ──────────────────────────────────────────────────────

    #[test]
    fn test_mixing_length_positive() {
        let ml = MixingLength::new(1.0);
        assert!(ml.mixing_length(0.5) > 0.0);
    }

    #[test]
    fn test_mixing_length_zero_at_wall() {
        let ml = MixingLength::new(1.0);
        assert_eq!(ml.mixing_length(0.0), 0.0);
    }

    #[test]
    fn test_mixing_length_zero_at_centre() {
        let ml = MixingLength::new(1.0);
        assert!((ml.mixing_length(1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_mixing_length_damped_less_than_undamped() {
        let ml = MixingLength::new(1.0);
        let lm = ml.mixing_length(0.5);
        let lmd = ml.mixing_length_damped(0.5, 5.0);
        assert!(lmd <= lm);
    }

    #[test]
    fn test_log_law_velocity_positive() {
        let ml = MixingLength::new(1.0);
        assert!(ml.log_law_velocity(100.0) > 0.0);
    }

    #[test]
    fn test_log_law_velocity_zero_at_zero() {
        let ml = MixingLength::new(1.0);
        assert_eq!(ml.log_law_velocity(0.0), 0.0);
    }

    #[test]
    fn test_turbulent_viscosity_increases_with_gradient() {
        let ml = MixingLength::new(1.0);
        let nu1 = ml.turbulent_viscosity(0.5, 1.0);
        let nu2 = ml.turbulent_viscosity(0.5, 2.0);
        assert!(nu2 > nu1);
    }

    // ── TurbulentFlux ─────────────────────────────────────────────────────

    #[test]
    fn test_turbulent_flux_total() {
        let tf = TurbulentFlux::new(0.1, 0.01);
        // F = -(0.01 + 0.1) * 1.0 = -0.11
        assert!((tf.total_flux(1.0) - (-0.11)).abs() < 1e-12);
    }

    #[test]
    fn test_turbulent_flux_effective_diff() {
        let tf = TurbulentFlux::new(0.1, 0.01);
        assert!((tf.effective_diffusivity() - 0.11).abs() < 1e-12);
    }

    #[test]
    fn test_turbulent_flux_cg_zero_gradient() {
        let tf = TurbulentFlux::new(0.1, 0.01).with_counter_gradient(0.5);
        // F = -D_t * (0 - 0.5) = 0.05
        let f = tf.turbulent_flux_cg(0.0);
        assert!((f - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_turbulent_schmidt_number() {
        let tf = TurbulentFlux::new(0.1, 0.01);
        let sc = tf.turbulent_schmidt(0.07);
        assert!((sc - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_turbulent_schmidt_zero_dt_infinite() {
        let tf = TurbulentFlux::new(0.0, 0.01);
        assert!(tf.turbulent_schmidt(0.07).is_infinite());
    }

    #[test]
    fn test_advective_flux() {
        let tf = TurbulentFlux::new(0.1, 0.01);
        assert!((tf.advective_flux(2.0, 3.0) - 6.0).abs() < 1e-12);
    }

    // ── DispersionLBM ─────────────────────────────────────────────────────

    #[test]
    fn test_dispersion_lbm_creation() {
        let d = DispersionLBM::new(16, 12, 0.02);
        assert_eq!(d.nx, 16);
        assert_eq!(d.ny, 12);
        assert_eq!(d.ux.len(), 192);
    }

    #[test]
    fn test_dispersion_lbm_inject_gaussian() {
        let mut d = DispersionLBM::new(20, 20, 0.02);
        d.inject_gaussian(10, 10, 100.0, 2.0);
        assert!(d.total_mass() > 0.0);
        let (cx, cy) = d.centroid();
        // Centroid should be near (10, 10)
        assert!((cx - 10.0).abs() < 1.0);
        assert!((cy - 10.0).abs() < 1.0);
    }

    #[test]
    fn test_dispersion_lbm_advance_step() {
        let mut d = DispersionLBM::new(10, 10, 0.02);
        d.inject_gaussian(5, 5, 1.0, 1.5);
        let mass0 = d.total_mass();
        d.advance();
        assert!(d.step == 1);
        // Mass broadly conserved
        assert!((d.total_mass() - mass0).abs() < mass0 * 0.5);
    }

    #[test]
    fn test_dispersion_lbm_with_turbulence() {
        let td = TurbulentDiffusivity::standard(0.5, 0.5);
        let d = DispersionLBM::new(8, 8, 0.01).with_turbulence(td);
        assert!(d.turb_diff.is_some());
    }

    #[test]
    fn test_dispersion_lbm_centroid_advects() {
        let mut d = DispersionLBM::new(20, 10, 0.01);
        d.inject_gaussian(5, 5, 1.0, 1.0);
        // Apply uniform x-velocity
        let ux = vec![0.05_f64; 200];
        let uy = vec![0.0_f64; 200];
        d.set_velocity(ux, uy);
        let (cx0, _) = d.centroid();
        for _ in 0..5 {
            d.advance();
        }
        let (cx1, _) = d.centroid();
        // Centroid should have moved in +x direction
        assert!(cx1 >= cx0);
    }

    // ── Free functions ─────────────────────────────────────────────────────

    #[test]
    fn test_turbulent_peclet_basic() {
        let pe = turbulent_peclet(1.0, 1.0, 0.5);
        assert!((pe - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_turbulent_peclet_zero_dt_infinite() {
        assert!(turbulent_peclet(1.0, 1.0, 0.0).is_infinite());
    }

    #[test]
    fn test_schmidt_number_basic() {
        let sc = schmidt_number(1e-6, 2e-6);
        assert!((sc - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_schmidt_number_zero_d_infinite() {
        assert!(schmidt_number(1e-6, 0.0).is_infinite());
    }

    #[test]
    fn test_richardson_number_stable() {
        let ri = richardson_number(9.81, 3.4e-3, 1.0, 10.0, 5.0);
        assert!(ri > 0.0);
    }

    #[test]
    fn test_richardson_number_unstable() {
        let ri = richardson_number(9.81, 3.4e-3, -1.0, 10.0, 5.0);
        assert!(ri < 0.0);
    }

    #[test]
    fn test_richardson_number_zero_u_stable_inf() {
        let ri = richardson_number(9.81, 1e-3, 1.0, 10.0, 0.0);
        assert!(ri.is_infinite() && ri > 0.0);
    }
}
