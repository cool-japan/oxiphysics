// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced shape memory alloy (SMA) constitutive modeling.
//!
//! This module provides a comprehensive suite of SMA models including the
//! Brinson, Tanaka (exponential kinetics), and Liang-Rogers (cosine kinetics)
//! frameworks. It covers:
//!
//! - **Phase fraction evolution** – martensite volume fraction ξ
//! - **Transformation temperatures** – Ms, Mf, As, Af with Clausius-Clapeyron shifts
//! - **Superelasticity** – stress-induced transformation above Af
//! - **Shape memory effect** – thermal recovery under constant load
//! - **Two-way shape memory** – training-dependent spontaneous transformation
//! - **Pseudoelastic stress-strain curve** generation
//! - **Training / cyclic loading** – evolution of internal state
//! - **Thermomechanical coupling** – latent heat during phase change
//! - **Material presets** – NiTi, CuAlNi, CuZnAl

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Default Clausius-Clapeyron slope for NiTi (Pa / K).
const DEFAULT_CC_SLOPE: f64 = 6.5e6;

/// Default latent heat of transformation (J / kg).
const DEFAULT_LATENT_HEAT: f64 = 24_000.0;

/// Default density of NiTi (kg / m^3).
const DEFAULT_DENSITY: f64 = 6450.0;

/// Default specific heat capacity (J / (kg K)).
const DEFAULT_SPECIFIC_HEAT: f64 = 837.0;

/// Small tolerance for numerical comparisons.
const EPS: f64 = 1.0e-15;

// ─────────────────────────────────────────────────────────────────────────────
// TransformationTemperatures
// ─────────────────────────────────────────────────────────────────────────────

/// The four characteristic transformation temperatures of an SMA.
///
/// Convention: `mf < ms < as_ < af` for a well-posed alloy.
#[derive(Debug, Clone, Copy)]
pub struct TransformationTemperatures {
    /// Martensite start temperature (K).
    pub ms: f64,
    /// Martensite finish temperature (K).
    pub mf: f64,
    /// Austenite start temperature (K).
    pub as_: f64,
    /// Austenite finish temperature (K).
    pub af: f64,
}

impl TransformationTemperatures {
    /// Create a new set of transformation temperatures.
    pub fn new(ms: f64, mf: f64, as_: f64, af: f64) -> Self {
        Self { ms, mf, as_, af }
    }

    /// Return temperatures shifted by the Clausius-Clapeyron relation.
    ///
    /// `delta_T = sigma / C_m` where `C_m` is the slope (Pa / K).
    pub fn shifted(&self, stress: f64, cm: f64, ca: f64) -> Self {
        let dt_m = if cm.abs() > EPS { stress / cm } else { 0.0 };
        let dt_a = if ca.abs() > EPS { stress / ca } else { 0.0 };
        Self {
            ms: self.ms + dt_m,
            mf: self.mf + dt_m,
            as_: self.as_ + dt_a,
            af: self.af + dt_a,
        }
    }

    /// Temperature range for martensitic transformation.
    pub fn martensite_range(&self) -> f64 {
        (self.ms - self.mf).abs()
    }

    /// Temperature range for austenitic transformation.
    pub fn austenite_range(&self) -> f64 {
        (self.af - self.as_).abs()
    }

    /// Hysteresis width `af - ms`.
    pub fn hysteresis(&self) -> f64 {
        self.af - self.ms
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SmaPhaseState
// ─────────────────────────────────────────────────────────────────────────────

/// Internal state variables for an SMA material point.
#[derive(Debug, Clone, Copy)]
pub struct SmaPhaseState {
    /// Total martensite volume fraction ξ ∈ \[0, 1\].
    pub xi: f64,
    /// Stress-induced martensite fraction ξ_S ∈ \[0, ξ\].
    pub xi_s: f64,
    /// Temperature-induced martensite fraction ξ_T = ξ − ξ_S.
    pub xi_t: f64,
    /// Current temperature (K).
    pub temperature: f64,
    /// Current uniaxial stress (Pa).
    pub stress: f64,
    /// Current uniaxial strain (dimensionless).
    pub strain: f64,
    /// Accumulated transformation strain.
    pub transformation_strain: f64,
    /// Transformation memory of the forward (loading / detwinning) branch of
    /// [`BrinsonModel`]: the progress `s ∈ [0, 1]` through its stress window
    /// that the current cosine ramp has reached.  The branch only transforms
    /// further when its progress exceeds this value, so elastic unloading and
    /// reloading retrace the same path.  `None` until the first update (the
    /// progress at the current stress and temperature is used instead).
    pub progress_forward: Option<f64>,
    /// Transformation memory of the reverse (heating / unloading) branch; see
    /// [`SmaPhaseState::progress_forward`].
    pub progress_reverse: Option<f64>,
    /// Transformation memory of the cooling branch; see
    /// [`SmaPhaseState::progress_forward`].
    pub progress_cooling: Option<f64>,
}

impl SmaPhaseState {
    /// Create a default phase state at a given temperature.
    pub fn new(temperature: f64) -> Self {
        Self {
            xi: 0.0,
            xi_s: 0.0,
            xi_t: 0.0,
            temperature,
            stress: 0.0,
            strain: 0.0,
            transformation_strain: 0.0,
            progress_forward: None,
            progress_reverse: None,
            progress_cooling: None,
        }
    }

    /// Create a fully martensitic state.
    pub fn fully_martensite(temperature: f64) -> Self {
        Self {
            xi: 1.0,
            xi_s: 0.0,
            xi_t: 1.0,
            temperature,
            stress: 0.0,
            strain: 0.0,
            transformation_strain: 0.0,
            progress_forward: None,
            progress_reverse: None,
            progress_cooling: None,
        }
    }

    /// Create a fully austenitic state.
    pub fn fully_austenite(temperature: f64) -> Self {
        Self::new(temperature)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TanakaModel
// ─────────────────────────────────────────────────────────────────────────────

/// Tanaka exponential kinetics model for SMA phase transformation.
///
/// The martensite fraction evolves as:
///
/// * Forward (A→M): `ξ = 1 − exp(−[a_M (Ms − T) + b_M σ])`
/// * Reverse (M→A): `ξ = exp(a_A (As − T) + b_A σ)`
///
/// with `a_M = ln(100)/(Ms − Mf)`, `a_A = ln(100)/(Af − As)` (so that the
/// transformation is 99 % complete at `Mf` / `Af`) and the stress constants
/// `b_M = a_M/C_M`, `b_A = a_A/C_A`.  A tensile stress raises the transformation
/// temperatures by `σ/C` (Clausius-Clapeyron), i.e. `Ms(σ) = Ms + σ/C_M`, which
/// is exactly what the `b σ` terms express: stress promotes the forward
/// transformation and retards the reverse one.
#[derive(Debug, Clone)]
pub struct TanakaModel {
    /// Transformation temperatures.
    pub temps: TransformationTemperatures,
    /// Clausius-Clapeyron slope for martensite, `C_M` (Pa / K).
    pub cm: f64,
    /// Clausius-Clapeyron slope for austenite, `C_A` (Pa / K).
    pub ca: f64,
    /// Exponential constant for martensitic transformation.
    pub a_m: f64,
    /// Stress constant for martensitic transformation.
    pub b_m: f64,
    /// Exponential constant for austenitic transformation.
    pub a_a: f64,
    /// Stress constant for austenitic transformation.
    pub b_a: f64,
    /// Maximum transformation strain.
    pub h_max: f64,
    /// Austenite elastic modulus (Pa).
    pub e_a: f64,
    /// Martensite elastic modulus (Pa).
    pub e_m: f64,
    /// Thermoelastic tensor Θ (Pa / K).
    pub theta: f64,
}

impl TanakaModel {
    /// Create a Tanaka model from transformation temperatures and material constants.
    ///
    /// # Parameters
    /// * `temps` – transformation temperatures
    /// * `cm` – Clausius-Clapeyron slope for martensite (Pa / K)
    /// * `ca` – Clausius-Clapeyron slope for austenite (Pa / K)
    /// * `h_max` – maximum transformation strain
    /// * `e_a` – austenite modulus (Pa)
    /// * `e_m` – martensite modulus (Pa)
    /// * `theta` – thermoelastic tensor (Pa / K)
    pub fn new(
        temps: TransformationTemperatures,
        cm: f64,
        ca: f64,
        h_max: f64,
        e_a: f64,
        e_m: f64,
        theta: f64,
    ) -> Self {
        // ln(100) = -ln(0.01): the exponential has decayed to 1 % at Mf / Af.
        let ln_100 = -(0.01_f64).ln();
        let dm = temps.ms - temps.mf;
        let da = temps.af - temps.as_;
        let a_m = if dm.abs() > EPS { ln_100 / dm } else { 0.0 };
        let a_a = if da.abs() > EPS { ln_100 / da } else { 0.0 };
        let b_m = if cm.abs() > EPS { a_m / cm } else { 0.0 };
        let b_a = if ca.abs() > EPS { a_a / ca } else { 0.0 };
        Self {
            temps,
            cm,
            ca,
            a_m,
            b_m,
            a_a,
            b_a,
            h_max,
            e_a,
            e_m,
            theta,
        }
    }

    /// Forward transformation (A → M): martensite fraction
    /// `ξ = 1 − exp(−[a_M (Ms − T) + b_M σ])`, clamped to `[0, 1]`.
    pub fn forward_fraction(&self, temp: f64, stress: f64) -> f64 {
        let arg = self.a_m * (self.temps.ms - temp) + self.b_m * stress;
        (1.0 - (-arg).exp()).clamp(0.0, 1.0)
    }

    /// Reverse transformation (M → A): martensite fraction
    /// `ξ = exp(a_A (As − T) + b_A σ)`, clamped to `[0, 1]`.
    pub fn reverse_fraction(&self, temp: f64, stress: f64) -> f64 {
        let arg = self.a_a * (self.temps.as_ - temp) + self.b_a * stress;
        arg.exp().clamp(0.0, 1.0)
    }

    /// Determine the phase fraction for a given (T, σ) state.
    ///
    /// The transformation temperatures are shifted by the model's own
    /// Clausius-Clapeyron slopes `C_M`, `C_A` ([`TanakaModel::cm`],
    /// [`TanakaModel::ca`]):
    ///
    /// * `T ≤ Ms(σ)`: forward branch.  Martensite already present is never
    ///   destroyed by temperatures below `As(σ)`, so the result is
    ///   `max(ξ_prev, forward_fraction)`.
    /// * `T ≥ As(σ)`: reverse branch, `min(ξ_prev, reverse_fraction)`
    ///   (austenite is not re-transformed above `Ms(σ)`).
    /// * `Ms(σ) < T < As(σ)`: the thermal hysteresis window; nothing
    ///   transforms and the previous fraction `xi_prev` is returned.
    pub fn phase_fraction(&self, temp: f64, stress: f64, xi_prev: f64) -> f64 {
        let shifted = self.temps.shifted(stress, self.cm, self.ca);
        if temp <= shifted.ms {
            xi_prev.max(self.forward_fraction(temp, stress))
        } else if temp >= shifted.as_ {
            xi_prev.min(self.reverse_fraction(temp, stress))
        } else {
            xi_prev
        }
    }

    /// Effective elastic modulus at martensite fraction ξ.
    pub fn effective_modulus(&self, xi: f64) -> f64 {
        self.e_a + xi.clamp(0.0, 1.0) * (self.e_m - self.e_a)
    }

    /// Constitutive stress for given strain, temperature, and state.
    pub fn stress_response(&self, strain: f64, temp: f64, xi: f64, temp_ref: f64) -> f64 {
        let e = self.effective_modulus(xi);
        e * (strain - self.h_max * xi) + self.theta * (temp - temp_ref)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LiangRogersModel
// ─────────────────────────────────────────────────────────────────────────────

/// Liang-Rogers cosine kinetics model for SMA phase transformation.
///
/// The martensite fraction evolves as:
///
/// * Forward (A→M): `ξ = (1 − ξ_0)/2 · cos[a_M(T − Mf) + b_M σ] + (1 + ξ_0)/2`
/// * Reverse (M→A): `ξ = ξ_0/2 · {cos[a_A(T − As) + b_A σ] + 1}`
///
/// where `a_M = π/(Ms − Mf)`, `a_A = π/(Af − As)`.
#[derive(Debug, Clone)]
pub struct LiangRogersModel {
    /// Transformation temperatures.
    pub temps: TransformationTemperatures,
    /// Cosine coefficient for martensitic transformation.
    pub a_m: f64,
    /// Stress shift for martensitic transformation.
    pub b_m: f64,
    /// Cosine coefficient for austenitic transformation.
    pub a_a: f64,
    /// Stress shift for austenitic transformation.
    pub b_a: f64,
    /// Maximum transformation strain.
    pub h_max: f64,
    /// Austenite elastic modulus (Pa).
    pub e_a: f64,
    /// Martensite elastic modulus (Pa).
    pub e_m: f64,
    /// Thermoelastic tensor (Pa / K).
    pub theta: f64,
}

impl LiangRogersModel {
    /// Create a Liang-Rogers model.
    pub fn new(
        temps: TransformationTemperatures,
        cm: f64,
        ca: f64,
        h_max: f64,
        e_a: f64,
        e_m: f64,
        theta: f64,
    ) -> Self {
        let dm = temps.ms - temps.mf;
        let da = temps.af - temps.as_;
        let a_m = if dm.abs() > EPS { PI / dm } else { 0.0 };
        let a_a = if da.abs() > EPS { PI / da } else { 0.0 };
        let b_m = if cm.abs() > EPS { -a_m / cm } else { 0.0 };
        let b_a = if ca.abs() > EPS { -a_a / ca } else { 0.0 };
        Self {
            temps,
            a_m,
            b_m,
            a_a,
            b_a,
            h_max,
            e_a,
            e_m,
            theta,
        }
    }

    /// Forward transformation fraction (A→M) given previous ξ₀.
    pub fn forward_fraction(&self, temp: f64, stress: f64, xi_prev: f64) -> f64 {
        let arg = self.a_m * (temp - self.temps.mf) + self.b_m * stress;
        let cos_val = arg.cos().clamp(-1.0, 1.0);
        let xi = (1.0 - xi_prev) / 2.0 * cos_val + (1.0 + xi_prev) / 2.0;
        xi.clamp(0.0, 1.0)
    }

    /// Reverse transformation fraction (M→A) given previous ξ₀.
    pub fn reverse_fraction(&self, temp: f64, stress: f64, xi_prev: f64) -> f64 {
        let arg = self.a_a * (temp - self.temps.as_) + self.b_a * stress;
        let cos_val = arg.cos().clamp(-1.0, 1.0);
        let xi = xi_prev / 2.0 * (cos_val + 1.0);
        xi.clamp(0.0, 1.0)
    }

    /// Determine the phase fraction for a given (T, σ) state.
    pub fn phase_fraction(&self, temp: f64, stress: f64, xi_prev: f64) -> f64 {
        let shifted = self
            .temps
            .shifted(stress, DEFAULT_CC_SLOPE, DEFAULT_CC_SLOPE);
        if temp <= shifted.ms && temp >= shifted.mf {
            self.forward_fraction(temp, stress, xi_prev)
        } else if temp >= shifted.as_ && temp <= shifted.af {
            self.reverse_fraction(temp, stress, xi_prev)
        } else if temp < shifted.mf {
            1.0
        } else if temp > shifted.af {
            0.0
        } else {
            xi_prev
        }
    }

    /// Effective elastic modulus at martensite fraction ξ.
    pub fn effective_modulus(&self, xi: f64) -> f64 {
        self.e_a + xi.clamp(0.0, 1.0) * (self.e_m - self.e_a)
    }

    /// Constitutive stress from the Liang-Rogers framework.
    pub fn stress_response(&self, strain: f64, temp: f64, xi: f64, temp_ref: f64) -> f64 {
        let e = self.effective_modulus(xi);
        e * (strain - self.h_max * xi) + self.theta * (temp - temp_ref)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BrinsonModel
// ─────────────────────────────────────────────────────────────────────────────

/// Remaining (untransformed) fraction `(1 + cos πs)/2` of a cosine-kinetics
/// transformation that has advanced to progress `s ∈ [0, 1]` of its window.
fn cosine_remaining(s: f64) -> f64 {
    0.5 * (1.0 + (PI * s).cos())
}

/// Fraction of the material that is still untransformed after a cosine-kinetics
/// branch advances from progress `s_old` to `s_new ≥ s_old`:
/// `(1 + cos π s_new) / (1 + cos π s_old)`.
///
/// Because the remaining fraction of a cosine ramp factorises, applying this
/// ratio step by step gives exactly the closed-form cosine kinetics from the
/// start of the transformation, independent of the step size.
fn cosine_survival_ratio(s_old: f64, s_new: f64) -> f64 {
    let old = cosine_remaining(s_old);
    if old <= 1.0e-30 {
        return 0.0;
    }
    (cosine_remaining(s_new) / old).clamp(0.0, 1.0)
}

/// Linear progress of `value` through the window `[start, finish]`, clamped to
/// `[0, 1]`.  A vanishing window acts as a step located at `finish`.
fn window_progress(value: f64, start: f64, finish: f64) -> f64 {
    let width = finish - start;
    if width <= EPS {
        return if value >= finish { 1.0 } else { 0.0 };
    }
    ((value - start) / width).clamp(0.0, 1.0)
}

/// Upper limit of bisection steps used to invert the constitutive law for the
/// stress (the bracket halves every step, so this is far beyond round-off).
const BISECTION_MAX_ITERATIONS: usize = 200;

/// Progress increase below which a kinetic branch counts as not advanced
/// (round-off in the stress must not erase the transformation memory).
const PROGRESS_TOLERANCE: f64 = 1.0e-12;

/// Result of one kinetic update of the [`BrinsonModel`].
struct KineticUpdate {
    /// Stress-induced martensite fraction ξ_S.
    xi_s: f64,
    /// Temperature-induced martensite fraction ξ_T.
    xi_t: f64,
    /// Transformation memory `[forward, reverse, cooling]` after the update.
    progress: [f64; 3],
}

/// Brinson model for SMA — separates stress-induced (ξ_S) and
/// temperature-induced (ξ_T) martensite fractions.
///
/// The total martensite fraction is ξ = ξ_S + ξ_T, where each component
/// follows cosine kinetics with separate evolution rules. This allows
/// capturing both superelasticity and shape memory effect within a single
/// unified framework.
///
/// **Constitutive law:** `σ = E(ξ)·(ε − ε_L·ξ_S) + Θ·(T − T_ref)` with the
/// rule-of-mixtures modulus `E(ξ) = E_A + ξ (E_M − E_A)`.
///
/// **Kinetics (Brinson 1993).**  Three branches, each a cosine ramp in its own
/// progress variable `s ∈ [0, 1]`:
///
/// * *forward* (loading, A→M and detwinning M_twinned→M_detwinned): `s` runs
///   over the stress window `σ_s^cr + C_M (T − Ms)⁺ … σ_f^cr + C_M (T − Ms)⁺`;
///   it drives `(ξ_S, ξ_T)` towards `(1, 0)`.  Below `Ms` this is the
///   detwinning of twinned martensite between the critical stresses
///   `σ_s^cr` and `σ_f^cr`.
/// * *reverse* (heating or unloading, M→A): `s` runs over
///   `As + σ/C_A … Af + σ/C_A` in temperature; it drives `(ξ_S, ξ_T)` to `(0, 0)`,
///   i.e. complete recovery above `Af(σ)`.
/// * *cooling* (A→twinned M): `s` runs over `Ms … Mf`; the remaining austenite
///   becomes twinned martensite, `ξ_T → 1 − ξ_S`, below `Mf`.  Detwinned
///   martensite is not re-twinned by cooling (the literal formula of Brinson's
///   paper would also drive `ξ_S` to zero); the window is stress free because
///   stress-induced martensite is described by the forward branch.
///
/// A branch acts only when its progress exceeds the progress stored in the
/// state (the *transformation memory*, [`SmaPhaseState::progress_forward`] and
/// friends); the transformable material then decays by the cosine ratio of
/// the old and new progress.  The result for a monotone path therefore equals
/// Brinson's closed form from the transformation-start fractions and does not
/// depend on the step size, and elastic unloading / reloading retraces the
/// same path.  Whenever a branch transforms, the memory of all branches
/// restarts from the current progress, because a new ramp begins from the new
/// fractions (as in Brinson's original model, an inner loop that reloads after
/// a partial reverse transformation therefore resumes forward transformation
/// with a new cosine ramp starting at the window edge).  The transformation
/// criteria use the tensile part of the stress, `max(σ, 0)`.
#[derive(Debug, Clone)]
pub struct BrinsonModel {
    /// Transformation temperatures.
    pub temps: TransformationTemperatures,
    /// Clausius-Clapeyron slope for martensite (Pa / K).
    pub cm: f64,
    /// Clausius-Clapeyron slope for austenite (Pa / K).
    pub ca: f64,
    /// Maximum transformation strain ε_L.
    pub h_max: f64,
    /// Austenite elastic modulus (Pa).
    pub e_a: f64,
    /// Martensite elastic modulus (Pa).
    pub e_m: f64,
    /// Thermoelastic tensor Θ (Pa / K).
    pub theta: f64,
    /// Critical start stress for detwinning at T < Ms (Pa).
    pub sigma_crit_s: f64,
    /// Critical finish stress for detwinning at T < Ms (Pa).
    pub sigma_crit_f: f64,
    /// Reference temperature `T_ref` (K) at which the thermoelastic stress
    /// `Θ·(T − T_ref)` vanishes.  [`BrinsonModel::new`] sets it to the austenite
    /// finish temperature `Af` (stress-free austenite is the natural undeformed
    /// configuration); use [`BrinsonModel::with_reference_temperature`] to set
    /// the temperature at which a specimen is stress free, e.g. the test
    /// temperature of an isothermal experiment.
    pub t_ref: f64,
}

impl BrinsonModel {
    /// Create a new Brinson model.
    ///
    /// The thermoelastic reference temperature [`BrinsonModel::t_ref`] is
    /// initialised to the austenite finish temperature.
    pub fn new(
        temps: TransformationTemperatures,
        cm: f64,
        ca: f64,
        h_max: f64,
        e_a: f64,
        e_m: f64,
        theta: f64,
        sigma_crit_s: f64,
        sigma_crit_f: f64,
    ) -> Self {
        let t_ref = temps.af;
        Self {
            temps,
            cm,
            ca,
            h_max,
            e_a,
            e_m,
            theta,
            sigma_crit_s,
            sigma_crit_f,
            t_ref,
        }
    }

    /// Return the model with the thermoelastic reference temperature set to
    /// `t_ref` (K), the temperature at which `Θ·(T − T_ref)` vanishes.
    pub fn with_reference_temperature(mut self, t_ref: f64) -> Self {
        self.t_ref = t_ref;
        self
    }

    /// Effective elastic modulus via the rule of mixtures.
    pub fn effective_modulus(&self, xi: f64) -> f64 {
        self.e_a + xi.clamp(0.0, 1.0) * (self.e_m - self.e_a)
    }

    /// Phase diagram critical stress at temperature T (stress-induced).
    ///
    /// Returns `(sigma_ms, sigma_mf)` — the start and finish stresses for
    /// stress-induced martensitic transformation at temperature T.
    pub fn critical_stresses_at_temp(&self, temp: f64) -> (f64, f64) {
        let sigma_ms = if temp > self.temps.ms {
            self.cm * (temp - self.temps.ms) + self.sigma_crit_s
        } else {
            self.sigma_crit_s
        };
        let sigma_mf = if temp > self.temps.ms {
            self.cm * (temp - self.temps.ms) + self.sigma_crit_f
        } else {
            self.sigma_crit_f
        };
        (sigma_ms, sigma_mf)
    }

    /// Progress of the forward (loading / detwinning) transformation through
    /// its stress window at `(temp, stress)`.
    fn forward_progress(&self, temp: f64, stress: f64) -> f64 {
        let (sig_s, sig_f) = self.critical_stresses_at_temp(temp);
        window_progress(stress.max(0.0), sig_s, sig_f)
    }

    /// Progress of the reverse (heating / unloading) transformation through its
    /// stress-shifted temperature window `As(σ) … Af(σ)`.
    fn reverse_progress(&self, temp: f64, stress: f64) -> f64 {
        let shifted = self.temps.shifted(stress.max(0.0), self.cm, self.ca);
        window_progress(temp, shifted.as_, shifted.af)
    }

    /// Progress of the stress-free cooling transformation through `Ms … Mf`.
    fn cooling_progress(&self, temp: f64) -> f64 {
        window_progress(self.temps.ms - temp, 0.0, self.temps.ms - self.temps.mf)
    }

    /// Compute stress-induced martensite fraction ξ_S for the loading branch.
    ///
    /// Uses cosine kinetics between `sigma_crit_s` and `sigma_crit_f`
    /// (shifted by `C_M (T − Ms)` above `Ms`) starting from the fraction
    /// `xi_s_prev` at the start of the transformation:
    /// below the window ξ_S stays `xi_s_prev`, above it the whole material is
    /// detwinned martensite, `ξ_S = 1`.  (Brinson's ξ_S does not depend on the
    /// twinned fraction, hence `_xi_t_prev` is unused.)
    pub fn xi_s_loading(&self, temp: f64, stress: f64, xi_s_prev: f64, _xi_t_prev: f64) -> f64 {
        let remaining = cosine_remaining(self.forward_progress(temp, stress));
        // ξ_S0·r + (1 − r) = 1 − (1 − ξ_S0)·r, exact at both ends of the window.
        (xi_s_prev * remaining + (1.0 - remaining)).clamp(0.0, 1.0)
    }

    /// Compute temperature-induced martensite fraction ξ_T for cooling.
    ///
    /// Cosine kinetics between the stress-shifted `Ms(σ)` and `Mf(σ)` starting
    /// from the fraction `xi_t_prev` at the start of the transformation: above
    /// `Ms(σ)` it stays `xi_t_prev`, below `Mf(σ)` it is 1.  ([`BrinsonModel::update_state`]
    /// uses the stress-free window, i.e. `stress = 0`, as in Brinson's
    /// original model, because stress-induced martensite is described by the
    /// forward branch.)
    pub fn xi_t_cooling(&self, temp: f64, stress: f64, xi_t_prev: f64) -> f64 {
        let shifted = self.temps.shifted(stress, self.cm, self.ca);
        let progress = window_progress(shifted.ms - temp, 0.0, shifted.ms - shifted.mf);
        let remaining = cosine_remaining(progress);
        // ξ_T0·r + (1 − r) = 1 − (1 − ξ_T0)·r, exact at both ends of the window.
        (xi_t_prev * remaining + (1.0 - remaining)).clamp(0.0, 1.0)
    }

    /// Reverse transformation (M → A): both ξ_S and ξ_T decrease.
    ///
    /// Cosine kinetics between the stress-shifted `As(σ)` and `Af(σ)`, starting
    /// from the fractions at the start of the transformation: below `As(σ)` the
    /// fractions are unchanged, above `Af(σ)` the material is fully austenitic,
    /// `(0, 0)`.  Both fractions decrease proportionally.
    pub fn reverse_fractions(
        &self,
        temp: f64,
        stress: f64,
        xi_s_prev: f64,
        xi_t_prev: f64,
    ) -> (f64, f64) {
        let remaining = cosine_remaining(self.reverse_progress(temp, stress));
        (xi_s_prev * remaining, xi_t_prev * remaining)
    }

    /// Integrated constitutive relation: σ = E(ξ)(ε − ε_L ξ_S) + Θ(T − T₀).
    pub fn stress_response(
        &self,
        strain: f64,
        xi_s: f64,
        xi_t: f64,
        temp: f64,
        temp_ref: f64,
    ) -> f64 {
        let xi = (xi_s + xi_t).clamp(0.0, 1.0);
        let e = self.effective_modulus(xi);
        e * (strain - self.h_max * xi_s) + self.theta * (temp - temp_ref)
    }

    /// Progress `[forward, reverse, cooling]` of the three kinetic branches at
    /// `(stress, temp)`.
    fn branch_progress(&self, stress: f64, temp: f64) -> [f64; 3] {
        [
            self.forward_progress(temp, stress),
            self.reverse_progress(temp, stress),
            self.cooling_progress(temp),
        ]
    }

    /// Transformation memory `[forward, reverse, cooling]` of a state: the
    /// stored branch progress, or, before the first update, the progress at the
    /// state's own stress and temperature.
    fn branch_frontiers(&self, state: &SmaPhaseState) -> [f64; 3] {
        let current = self.branch_progress(state.stress, state.temperature);
        [
            state.progress_forward.unwrap_or(current[0]),
            state.progress_reverse.unwrap_or(current[1]),
            state.progress_cooling.unwrap_or(current[2]),
        ]
    }

    /// Martensite fractions and branch memory reached at stress `stress` and
    /// temperature `temp` when starting from the state `prev`.
    ///
    /// The reverse, forward and cooling branches are applied in this order,
    /// each only if its progress exceeds the progress stored in `prev` (the
    /// transformation memory).  Every application is a convex combination of
    /// the current fractions and the branch's target, so `ξ_S, ξ_T ≥ 0` and
    /// `ξ_S + ξ_T ≤ 1` are preserved.  If any branch transformed, the memory of
    /// every branch restarts from its current progress (a new ramp begins from
    /// the new fractions); otherwise the memory is unchanged.
    fn kinetics_at(&self, prev: &SmaPhaseState, stress: f64, temp: f64) -> KineticUpdate {
        let mut xi_s = prev.xi_s.clamp(0.0, 1.0);
        let mut xi_t = prev.xi_t.max(0.0).min((1.0 - xi_s).max(0.0));
        let frontier = self.branch_frontiers(prev);
        let current = self.branch_progress(stress, temp);
        let advanced = [
            current[0] > frontier[0] + PROGRESS_TOLERANCE,
            current[1] > frontier[1] + PROGRESS_TOLERANCE,
            current[2] > frontier[2] + PROGRESS_TOLERANCE,
        ];

        // Reverse transformation: heating through, or unloading into, the austenite window.
        if advanced[1] {
            let keep = cosine_survival_ratio(frontier[1], current[1]);
            xi_s *= keep;
            xi_t *= keep;
        }
        // Forward transformation: loading / detwinning towards detwinned martensite.
        if advanced[0] {
            let keep = cosine_survival_ratio(frontier[0], current[0]);
            xi_s = xi_s * keep + (1.0 - keep);
            xi_t *= keep;
        }
        // Cooling: remaining austenite becomes twinned martensite; detwinned martensite is kept.
        if advanced[2] {
            let keep = cosine_survival_ratio(frontier[2], current[2]);
            let max_twinned = (1.0 - xi_s).max(0.0);
            xi_t = xi_t * keep + max_twinned * (1.0 - keep);
        }
        let progress = if advanced.iter().any(|&a| a) {
            current
        } else {
            frontier
        };
        KineticUpdate {
            xi_s,
            xi_t,
            progress,
        }
    }

    /// Store a new material-point state.
    fn commit_state(
        &self,
        state: &mut SmaPhaseState,
        strain: f64,
        stress: f64,
        temp: f64,
        update: &KineticUpdate,
    ) {
        state.xi_s = update.xi_s;
        state.xi_t = update.xi_t;
        state.xi = (update.xi_s + update.xi_t).clamp(0.0, 1.0);
        state.strain = strain;
        state.temperature = temp;
        state.stress = stress;
        state.transformation_strain = self.h_max * update.xi_s;
        state.progress_forward = Some(update.progress[0]);
        state.progress_reverse = Some(update.progress[1]);
        state.progress_cooling = Some(update.progress[2]);
    }

    /// Strain-driven update of an SMA material point to `(strain_new, temp_new)`.
    ///
    /// The martensite fractions depend on the stress through the transformation
    /// criteria, and the stress depends on the fractions through the
    /// constitutive law `σ = E(ξ)(ε − ε_L ξ_S) + Θ (T − T_ref)`; the two are
    /// solved together: the stress is the root of
    /// `σ − [E(ξ(σ))(ε − ε_L ξ_S(σ)) + Θ (T − T_ref)]`, found by bisection on
    /// the bracket `±(E_max (|ε| + ε_L) + |Θ (T − T_ref)|)` that contains every
    /// admissible stress.  The previous state supplies the starting fractions
    /// and the previous progress of each kinetic branch (see [`BrinsonModel`]).
    ///
    /// Consequently the model detwins twinned martensite below `Ms` under load
    /// (`σ_s^cr ≤ σ ≤ σ_f^cr`, saturating at `ξ_S = 1` above `σ_f^cr`),
    /// transforms completely to austenite above `Af(σ)`, and to twinned
    /// martensite on cooling below `Mf`.  Non-finite inputs leave the state
    /// unchanged.
    pub fn update_state(&self, state: &mut SmaPhaseState, strain_new: f64, temp_new: f64) {
        if !(strain_new.is_finite() && temp_new.is_finite()) {
            return;
        }
        let previous = *state;
        let thermal = self.theta * (temp_new - self.t_ref);
        let e_max = self.e_a.abs().max(self.e_m.abs());
        let bound = e_max * (strain_new.abs() + self.h_max.abs()) + thermal.abs();
        let residual = |sigma: f64| -> f64 {
            let update = self.kinetics_at(&previous, sigma, temp_new);
            let law =
                self.stress_response(strain_new, update.xi_s, update.xi_t, temp_new, self.t_ref);
            sigma - law
        };
        let (mut lo, mut hi) = (-bound, bound);
        for _ in 0..BISECTION_MAX_ITERATIONS {
            let mid = 0.5 * (lo + hi);
            let value = residual(mid);
            if value > 0.0 {
                hi = mid;
            } else if value < 0.0 {
                lo = mid;
            } else {
                lo = mid;
                hi = mid;
            }
            if hi - lo <= f64::EPSILON * mid.abs().max(1.0) {
                break;
            }
        }
        let sigma = 0.5 * (lo + hi);
        let update = self.kinetics_at(&previous, sigma, temp_new);
        self.commit_state(state, strain_new, sigma, temp_new, &update);
    }

    /// Stress-driven update of an SMA material point to `(stress_new, temp_new)`.
    ///
    /// With the stress prescribed the martensite fractions follow directly from
    /// the kinetics and the strain from the constitutive law,
    /// `ε = (σ − Θ (T − T_ref)) / E(ξ) + ε_L ξ_S`.  This is the natural way to
    /// simulate stress-free unloading and thermal recovery
    /// (`stress_new = 0`).  Non-finite inputs, or a non-positive modulus, leave
    /// the state unchanged.
    pub fn update_state_stress_controlled(
        &self,
        state: &mut SmaPhaseState,
        stress_new: f64,
        temp_new: f64,
    ) {
        if !(stress_new.is_finite() && temp_new.is_finite()) {
            return;
        }
        let update = self.kinetics_at(state, stress_new, temp_new);
        let modulus = self.effective_modulus((update.xi_s + update.xi_t).clamp(0.0, 1.0));
        if !(modulus.is_finite() && modulus > 0.0) {
            return;
        }
        let strain = (stress_new - self.theta * (temp_new - self.t_ref)) / modulus
            + self.h_max * update.xi_s;
        self.commit_state(state, strain, stress_new, temp_new, &update);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SuperelasticCurveGenerator
// ─────────────────────────────────────────────────────────────────────────────

/// Generates a pseudoelastic (superelastic) stress-strain curve at constant
/// temperature above Af.
///
/// The curve has four branches:
/// 1. Elastic austenite loading
/// 2. Stress-induced A→M plateau (loading)
/// 3. Elastic martensite loading beyond plateau
/// 4. Unloading with reverse M→A plateau (lower stress)
#[derive(Debug, Clone)]
pub struct SuperelasticCurveGenerator {
    /// Brinson model for constitutive response.
    pub model: BrinsonModel,
    /// Operating temperature (K), must be > Af.
    pub temperature: f64,
    /// Number of points to generate per branch.
    pub resolution: usize,
}

impl SuperelasticCurveGenerator {
    /// Create a new generator.
    pub fn new(model: BrinsonModel, temperature: f64, resolution: usize) -> Self {
        Self {
            model,
            temperature,
            resolution,
        }
    }

    /// The model used for the isothermal test: the specimen is stress free at
    /// the operating temperature, so the thermoelastic reference temperature is
    /// the operating temperature and the curve starts at the origin.
    fn isothermal_model(&self) -> BrinsonModel {
        self.model
            .clone()
            .with_reference_temperature(self.temperature)
    }

    /// Generate loading branch stress-strain data.
    ///
    /// Returns a vector of `(strain, stress)` pairs.
    pub fn loading_curve(&self, max_strain: f64) -> Vec<(f64, f64)> {
        let n = self.resolution.max(2);
        let mut result = Vec::with_capacity(n);
        let model = self.isothermal_model();
        let mut state = SmaPhaseState::new(self.temperature);
        for i in 0..n {
            let eps = max_strain * (i as f64) / (n as f64 - 1.0);
            model.update_state(&mut state, eps, self.temperature);
            result.push((eps, state.stress));
        }
        result
    }

    /// Generate unloading branch stress-strain data starting from a given state.
    pub fn unloading_curve(&self, state_at_peak: &SmaPhaseState) -> Vec<(f64, f64)> {
        let n = self.resolution.max(2);
        let mut result = Vec::with_capacity(n);
        let peak_strain = state_at_peak.strain;
        let model = self.isothermal_model();
        let mut state = *state_at_peak;
        for i in 0..n {
            let eps = peak_strain * (1.0 - (i as f64) / (n as f64 - 1.0));
            // Simple unloading: update with decreasing strain
            model.update_state(&mut state, eps, self.temperature);
            result.push((eps, state.stress.max(0.0)));
        }
        result
    }

    /// Full loading-unloading loop.
    pub fn full_loop(&self, max_strain: f64) -> Vec<(f64, f64)> {
        let loading = self.loading_curve(max_strain);
        let last_state = {
            let model = self.isothermal_model();
            let mut s = SmaPhaseState::new(self.temperature);
            for &(eps, _) in &loading {
                model.update_state(&mut s, eps, self.temperature);
            }
            s
        };
        let unloading = self.unloading_curve(&last_state);
        let mut result = loading;
        result.extend(unloading);
        result
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeMemoryEffect
// ─────────────────────────────────────────────────────────────────────────────

/// Models the shape memory effect: strain recovery upon heating.
///
/// A specimen is deformed in the martensitic state (T < Mf), then heated
/// above Af to recover its original shape.
#[derive(Debug, Clone)]
pub struct ShapeMemoryEffect {
    /// Underlying Brinson model.
    pub model: BrinsonModel,
    /// Reference temperature for the deformation step.
    pub deform_temperature: f64,
    /// Recovery temperature (must be > Af).
    pub recovery_temperature: f64,
}

impl ShapeMemoryEffect {
    /// Create a new shape memory effect model.
    pub fn new(model: BrinsonModel, deform_temp: f64, recovery_temp: f64) -> Self {
        Self {
            model,
            deform_temperature: deform_temp,
            recovery_temperature: recovery_temp,
        }
    }

    /// Simulate deformation at low temperature, then thermal recovery.
    ///
    /// 1. The twinned martensite at the deformation temperature is strained to
    ///    `applied_strain` (strain-controlled); the load detwins it.
    /// 2. The specimen is unloaded to zero stress (stress-controlled); the
    ///    detwinned martensite is retained and leaves a residual strain.
    /// 3. It is heated stress-free to the recovery temperature; the reverse
    ///    transformation (completed at `Af`) recovers the residual strain.
    ///
    /// The specimen is stress free at the deformation temperature, which is
    /// therefore the thermoelastic reference temperature of the simulation.
    ///
    /// Returns `(residual_strain_before_heating, residual_strain_after_heating)`,
    /// the transformation strains `ε_L·ξ_S` of steps 2 and 3.
    pub fn simulate_recovery(&self, applied_strain: f64) -> (f64, f64) {
        let model = self
            .model
            .clone()
            .with_reference_temperature(self.deform_temperature);

        // Step 1: deform the twinned martensite (full martensite at T < Mf).
        let mut state = SmaPhaseState::fully_martensite(self.deform_temperature);
        model.update_state(&mut state, applied_strain, self.deform_temperature);

        // Step 2: unload mechanically to zero stress.
        model.update_state_stress_controlled(&mut state, 0.0, self.deform_temperature);
        let residual_before = state.transformation_strain;

        // Step 3: heat stress free to the recovery temperature.
        model.update_state_stress_controlled(&mut state, 0.0, self.recovery_temperature);
        let residual_after = state.transformation_strain;

        (residual_before, residual_after)
    }

    /// Compute the recovery ratio η = (ε_initial − ε_residual) / ε_initial.
    pub fn recovery_ratio(&self, applied_strain: f64) -> f64 {
        let (before, after) = self.simulate_recovery(applied_strain);
        if before.abs() < EPS {
            return 0.0;
        }
        ((before - after) / before).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TwoWayShapeMemory
// ─────────────────────────────────────────────────────────────────────────────

/// Two-way shape memory effect: the material spontaneously deforms during
/// both cooling and heating without external load, after training.
///
/// The trained material has internal residual stresses from dislocation
/// arrays that bias the variant selection during cooling.
#[derive(Debug, Clone)]
pub struct TwoWayShapeMemory {
    /// Brinson model for constitutive relations.
    pub model: BrinsonModel,
    /// Two-way recoverable strain (fraction of h_max, typically 0.2-0.5).
    pub two_way_fraction: f64,
    /// Number of training cycles completed.
    pub training_cycles: u32,
    /// Maximum two-way fraction after saturation.
    pub max_two_way_fraction: f64,
    /// Training saturation rate constant.
    pub saturation_rate: f64,
}

impl TwoWayShapeMemory {
    /// Create a new two-way shape memory model.
    pub fn new(model: BrinsonModel) -> Self {
        Self {
            model,
            two_way_fraction: 0.0,
            training_cycles: 0,
            max_two_way_fraction: 0.4,
            saturation_rate: 0.1,
        }
    }

    /// Apply one training cycle (loading + unloading + thermal cycle).
    pub fn train_cycle(&mut self) {
        self.training_cycles += 1;
        // Exponential saturation: f = f_max (1 − exp(−rate × N))
        let n = self.training_cycles as f64;
        self.two_way_fraction =
            self.max_two_way_fraction * (1.0 - (-self.saturation_rate * n).exp());
    }

    /// Spontaneous strain on cooling to T < Mf (no external load).
    pub fn cooling_strain(&self, temp: f64) -> f64 {
        let xi = if temp <= self.model.temps.mf {
            1.0
        } else if temp >= self.model.temps.ms {
            0.0
        } else {
            let dm = self.model.temps.ms - self.model.temps.mf;
            if dm.abs() < EPS {
                0.0
            } else {
                0.5 * (PI * (temp - self.model.temps.mf) / dm).cos() + 0.5
            }
        };
        self.two_way_fraction * self.model.h_max * xi
    }

    /// Spontaneous strain on heating to T > Af (recovery).
    pub fn heating_strain(&self, temp: f64) -> f64 {
        let xi = if temp >= self.model.temps.af {
            0.0
        } else if temp <= self.model.temps.as_ {
            1.0
        } else {
            let da = self.model.temps.af - self.model.temps.as_;
            if da.abs() < EPS {
                0.0
            } else {
                0.5 * (PI * (temp - self.model.temps.as_) / da).cos() + 0.5
            }
        };
        self.two_way_fraction * self.model.h_max * xi
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrainingEffect
// ─────────────────────────────────────────────────────────────────────────────

/// Models the evolution of SMA behaviour under cyclic mechanical loading.
///
/// Repeated loading-unloading cycles cause:
/// - Accumulation of residual (irrecoverable) strain
/// - Stabilisation of the transformation plateau stress
/// - Reduction in hysteresis width
/// - Development of two-way shape memory
#[derive(Debug, Clone)]
pub struct TrainingEffect {
    /// Number of completed cycles.
    pub cycles: u32,
    /// Accumulated residual strain per cycle (list).
    pub residual_strains: Vec<f64>,
    /// Transformation stress start per cycle (list).
    pub plateau_stresses: Vec<f64>,
    /// Base (first cycle) residual strain per cycle.
    pub base_residual: f64,
    /// Saturation residual strain.
    pub saturated_residual: f64,
    /// Rate constant for residual strain saturation.
    pub rate_residual: f64,
    /// Base plateau stress (Pa).
    pub base_plateau: f64,
    /// Plateau stress reduction per decade of cycles.
    pub plateau_drop_per_decade: f64,
}

impl TrainingEffect {
    /// Create a new training model with defaults typical for NiTi.
    pub fn new() -> Self {
        Self {
            cycles: 0,
            residual_strains: Vec::new(),
            plateau_stresses: Vec::new(),
            base_residual: 0.005,
            saturated_residual: 0.02,
            rate_residual: 0.05,
            base_plateau: 500.0e6,
            plateau_drop_per_decade: 20.0e6,
        }
    }
}
impl Default for TrainingEffect {
    fn default() -> Self {
        Self::new()
    }
}
impl TrainingEffect {
    /// Record one additional training cycle and update residual and plateau stress histories.
    pub fn add_cycle(&mut self) {
        self.cycles += 1;
        let n = self.cycles as f64;
        // Residual strain accumulates with saturation
        let eps_r =
            self.saturated_residual * (1.0 - (-self.rate_residual * n).exp()) + self.base_residual;
        self.residual_strains.push(eps_r);
        // Plateau stress drops logarithmically
        let sigma_p = self.base_plateau - self.plateau_drop_per_decade * n.ln().max(0.0);
        self.plateau_stresses.push(sigma_p.max(0.0));
    }

    /// Current residual strain after N cycles.
    pub fn current_residual(&self) -> f64 {
        self.residual_strains.last().copied().unwrap_or(0.0)
    }

    /// Current plateau stress after N cycles.
    pub fn current_plateau_stress(&self) -> f64 {
        self.plateau_stresses
            .last()
            .copied()
            .unwrap_or(self.base_plateau)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermomechanicalCoupling
// ─────────────────────────────────────────────────────────────────────────────

/// Thermomechanical coupling: latent heat generation/absorption during
/// phase transformation affects the temperature of the SMA element.
///
/// Energy balance: `ρ c_p dT/dt = −L dξ/dt + Q_ext`
///
/// where L is the latent heat (J / kg), ρ is density, c_p is specific heat.
#[derive(Debug, Clone)]
pub struct ThermomechanicalCoupling {
    /// Density (kg / m^3).
    pub density: f64,
    /// Specific heat capacity (J / (kg K)).
    pub specific_heat: f64,
    /// Latent heat of transformation (J / kg).
    pub latent_heat: f64,
    /// Convective heat transfer coefficient (W / (m^2 K)).
    pub h_conv: f64,
    /// Surface area to volume ratio (1 / m).
    pub sa_ratio: f64,
    /// Ambient temperature (K).
    pub ambient_temp: f64,
}

impl ThermomechanicalCoupling {
    /// Create a new thermomechanical coupling model with NiTi-like defaults.
    pub fn new() -> Self {
        Self {
            density: DEFAULT_DENSITY,
            specific_heat: DEFAULT_SPECIFIC_HEAT,
            latent_heat: DEFAULT_LATENT_HEAT,
            h_conv: 10.0,
            sa_ratio: 200.0,
            ambient_temp: 300.0,
        }
    }
}
impl Default for ThermomechanicalCoupling {
    fn default() -> Self {
        Self::new()
    }
}
impl ThermomechanicalCoupling {
    ///
    /// `dT = −L / c_p · dξ`
    ///
    /// Forward transformation (dξ > 0) releases heat → temperature increases.
    pub fn temperature_increment(&self, d_xi: f64) -> f64 {
        -self.latent_heat / self.specific_heat * d_xi
    }

    /// Full temperature update including convection (explicit Euler).
    pub fn update_temperature(&self, temp: f64, d_xi: f64, dt: f64) -> f64 {
        let q_latent = -self.latent_heat * self.density * d_xi / dt;
        let q_conv = self.h_conv * self.sa_ratio * (self.ambient_temp - temp);
        let dt_temp = (q_latent + q_conv) / (self.density * self.specific_heat) * dt;
        temp + dt_temp
    }

    /// Adiabatic temperature change (no convection, instantaneous).
    pub fn adiabatic_temperature_change(&self, d_xi: f64) -> f64 {
        self.temperature_increment(d_xi)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ClausiusClapeyron
// ─────────────────────────────────────────────────────────────────────────────

/// Clausius-Clapeyron relation for stress-temperature phase diagram slopes.
///
/// `dσ/dT = −ΔH / (T₀ ε_L)` in some formulations, or empirically measured
/// as a constant slope C_m (Pa / K).
#[derive(Debug, Clone, Copy)]
pub struct ClausiusClapeyron {
    /// Martensite slope (Pa / K).
    pub cm: f64,
    /// Austenite slope (Pa / K).
    pub ca: f64,
}

impl ClausiusClapeyron {
    /// Create with equal slopes.
    pub fn symmetric(slope: f64) -> Self {
        Self {
            cm: slope,
            ca: slope,
        }
    }

    /// Create with distinct martensite and austenite slopes.
    pub fn new(cm: f64, ca: f64) -> Self {
        Self { cm, ca }
    }

    /// Stress shift for a given temperature deviation from a reference.
    pub fn stress_shift_martensite(&self, delta_t: f64) -> f64 {
        self.cm * delta_t
    }

    /// Stress shift for austenite transformation.
    pub fn stress_shift_austenite(&self, delta_t: f64) -> f64 {
        self.ca * delta_t
    }

    /// Temperature shift for a given stress (martensite side).
    pub fn temp_shift_martensite(&self, stress: f64) -> f64 {
        if self.cm.abs() > EPS {
            stress / self.cm
        } else {
            0.0
        }
    }

    /// Temperature shift for austenite side.
    pub fn temp_shift_austenite(&self, stress: f64) -> f64 {
        if self.ca.abs() > EPS {
            stress / self.ca
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Material Presets
// ─────────────────────────────────────────────────────────────────────────────

/// NiTi (Nitinol) preset parameters.
#[derive(Debug, Clone)]
pub struct NiTiPreset;

impl NiTiPreset {
    /// Transformation temperatures for a typical near-equiatomic NiTi.
    pub fn temperatures() -> TransformationTemperatures {
        TransformationTemperatures::new(291.0, 271.0, 295.0, 315.0)
    }

    /// Create a Brinson model for NiTi.
    pub fn brinson() -> BrinsonModel {
        BrinsonModel::new(
            Self::temperatures(),
            6.5e6,   // cm
            6.5e6,   // ca
            0.067,   // h_max
            67.0e9,  // E_A
            26.3e9,  // E_M
            0.55e6,  // Θ
            100.0e6, // σ_crit_s
            170.0e6, // σ_crit_f
        )
    }

    /// Create a Tanaka model for NiTi.
    pub fn tanaka() -> TanakaModel {
        let temps = Self::temperatures();
        TanakaModel::new(temps, 6.5e6, 6.5e6, 0.067, 67.0e9, 26.3e9, 0.55e6)
    }

    /// Create a Liang-Rogers model for NiTi.
    pub fn liang_rogers() -> LiangRogersModel {
        let temps = Self::temperatures();
        LiangRogersModel::new(temps, 6.5e6, 6.5e6, 0.067, 67.0e9, 26.3e9, 0.55e6)
    }

    /// Clausius-Clapeyron slopes for NiTi.
    pub fn clausius_clapeyron() -> ClausiusClapeyron {
        ClausiusClapeyron::symmetric(6.5e6)
    }

    /// Thermomechanical coupling with NiTi defaults.
    pub fn thermomechanical() -> ThermomechanicalCoupling {
        ThermomechanicalCoupling {
            density: 6450.0,
            specific_heat: 837.0,
            latent_heat: 24_000.0,
            h_conv: 10.0,
            sa_ratio: 200.0,
            ambient_temp: 300.0,
        }
    }
}

/// CuAlNi preset parameters.
#[derive(Debug, Clone)]
pub struct CuAlNiPreset;

impl CuAlNiPreset {
    /// Transformation temperatures for CuAlNi.
    pub fn temperatures() -> TransformationTemperatures {
        TransformationTemperatures::new(291.0, 253.0, 295.0, 323.0)
    }

    /// Create a Brinson model for CuAlNi.
    pub fn brinson() -> BrinsonModel {
        BrinsonModel::new(
            Self::temperatures(),
            4.5e6,   // cm
            4.5e6,   // ca
            0.05,    // h_max
            85.0e9,  // E_A
            70.0e9,  // E_M
            0.4e6,   // Θ
            80.0e6,  // σ_crit_s
            150.0e6, // σ_crit_f
        )
    }

    /// Create a Tanaka model for CuAlNi.
    pub fn tanaka() -> TanakaModel {
        TanakaModel::new(
            Self::temperatures(),
            4.5e6,
            4.5e6,
            0.05,
            85.0e9,
            70.0e9,
            0.4e6,
        )
    }

    /// Clausius-Clapeyron slopes for CuAlNi.
    pub fn clausius_clapeyron() -> ClausiusClapeyron {
        ClausiusClapeyron::symmetric(4.5e6)
    }
}

/// CuZnAl preset parameters.
#[derive(Debug, Clone)]
pub struct CuZnAlPreset;

impl CuZnAlPreset {
    /// Transformation temperatures for CuZnAl.
    pub fn temperatures() -> TransformationTemperatures {
        TransformationTemperatures::new(280.0, 255.0, 295.0, 310.0)
    }

    /// Create a Brinson model for CuZnAl.
    pub fn brinson() -> BrinsonModel {
        BrinsonModel::new(
            Self::temperatures(),
            5.0e6,   // cm
            5.0e6,   // ca
            0.04,    // h_max
            72.0e9,  // E_A
            52.0e9,  // E_M
            0.35e6,  // Θ
            90.0e6,  // σ_crit_s
            160.0e6, // σ_crit_f
        )
    }

    /// Create a Tanaka model for CuZnAl.
    pub fn tanaka() -> TanakaModel {
        TanakaModel::new(
            Self::temperatures(),
            5.0e6,
            5.0e6,
            0.04,
            72.0e9,
            52.0e9,
            0.35e6,
        )
    }

    /// Clausius-Clapeyron slopes for CuZnAl.
    pub fn clausius_clapeyron() -> ClausiusClapeyron {
        ClausiusClapeyron::symmetric(5.0e6)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the dissipated energy (area inside a hysteresis loop).
///
/// Given loading and unloading stress-strain curves as vectors of `(ε, σ)`.
pub fn hysteresis_energy(loading: &[(f64, f64)], unloading: &[(f64, f64)]) -> f64 {
    let area_load = trapz_area(loading);
    let area_unload = trapz_area(unloading);
    (area_load - area_unload).abs()
}

/// Trapezoidal integration of a curve given as `(x, y)` pairs.
fn trapz_area(curve: &[(f64, f64)]) -> f64 {
    if curve.len() < 2 {
        return 0.0;
    }
    let mut area = 0.0;
    for w in curve.windows(2) {
        let (x0, y0) = w[0];
        let (x1, y1) = w[1];
        area += 0.5 * (y0 + y1) * (x1 - x0);
    }
    area.abs()
}

/// Generate a temperature sweep stress-strain path for shape memory effect
/// demonstration.
///
/// A specimen of detwinned martensite storing the residual transformation
/// strain `initial_strain` (at most `ε_L`) at `t_start` is heated stress free
/// to `t_end` in `n_steps` temperature steps (at least two); the stored strain
/// is recovered while the temperature passes the stress-free window
/// `As … Af`.
///
/// Returns a vector of `(temperature, transformation strain)` pairs showing
/// the recovery.
pub fn temperature_sweep_recovery(
    model: &BrinsonModel,
    initial_strain: f64,
    t_start: f64,
    t_end: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let n = n_steps.max(2);
    let mut result = Vec::with_capacity(n);
    let model = model.clone().with_reference_temperature(t_start);
    let xi_s0 = if model.h_max > EPS {
        (initial_strain / model.h_max).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut state = SmaPhaseState::fully_martensite(t_start);
    state.xi_s = xi_s0;
    state.xi_t = 1.0 - xi_s0;
    state.strain = model.h_max * xi_s0;
    state.transformation_strain = model.h_max * xi_s0;

    for i in 0..n {
        let t_frac = i as f64 / (n as f64 - 1.0);
        let temp = t_start + t_frac * (t_end - t_start);
        model.update_state_stress_controlled(&mut state, 0.0, temp);
        result.push((temp, state.transformation_strain));
    }
    result
}

/// Compute the secant modulus from a stress-strain point.
pub fn secant_modulus(strain: f64, stress: f64) -> f64 {
    if strain.abs() < EPS {
        0.0
    } else {
        stress / strain
    }
}

/// Compute tangent modulus from two consecutive stress-strain points.
pub fn tangent_modulus(strain0: f64, stress0: f64, strain1: f64, stress1: f64) -> f64 {
    let de = strain1 - strain0;
    if de.abs() < EPS {
        0.0
    } else {
        (stress1 - stress0) / de
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
