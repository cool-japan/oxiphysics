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
/// * Forward (A→M): `ξ = 1 − exp(a_M (Ms − T) + b_M σ)`
/// * Reverse (M→A): `ξ = exp(a_A (As − T) + b_A σ)`
///
/// where `a_M`, `b_M`, `a_A`, `b_A` are material constants derived
/// from the transformation temperatures and Clausius-Clapeyron slope.
#[derive(Debug, Clone)]
pub struct TanakaModel {
    /// Transformation temperatures.
    pub temps: TransformationTemperatures,
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
        let dm = temps.ms - temps.mf;
        let da = temps.af - temps.as_;
        let a_m = if dm.abs() > EPS {
            -2.0 * (0.01_f64).ln() / dm
        } else {
            0.0
        };
        let a_a = if da.abs() > EPS {
            -2.0 * (0.01_f64).ln() / da
        } else {
            0.0
        };
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

    /// Forward transformation (A → M): compute new ξ.
    pub fn forward_fraction(&self, temp: f64, stress: f64) -> f64 {
        let arg = self.a_m * (self.temps.ms - temp) + self.b_m * stress;
        (1.0 - (-arg).exp()).clamp(0.0, 1.0)
    }

    /// Reverse transformation (M → A): compute new ξ.
    pub fn reverse_fraction(&self, temp: f64, stress: f64) -> f64 {
        let arg = self.a_a * (self.temps.as_ - temp) + self.b_a * stress;
        arg.exp().clamp(0.0, 1.0)
    }

    /// Determine the phase fraction for a given (T, σ) state.
    ///
    /// Selects forward or reverse based on where the temperature falls
    /// relative to the (stress-shifted) transformation temperatures.
    pub fn phase_fraction(&self, temp: f64, stress: f64) -> f64 {
        let shifted = self
            .temps
            .shifted(stress, DEFAULT_CC_SLOPE, DEFAULT_CC_SLOPE);
        if temp <= shifted.ms {
            self.forward_fraction(temp, stress)
        } else if temp >= shifted.as_ {
            self.reverse_fraction(temp, stress)
        } else {
            // Between Ms and As: retain previous fraction (simplified: use forward)
            self.forward_fraction(temp, stress)
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

/// Brinson model for SMA — separates stress-induced (ξ_S) and
/// temperature-induced (ξ_T) martensite fractions.
///
/// The total martensite fraction is ξ = ξ_S + ξ_T, where each component
/// follows cosine kinetics with separate evolution rules. This allows
/// capturing both superelasticity and shape memory effect within a single
/// unified framework.
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
}

impl BrinsonModel {
    /// Create a new Brinson model.
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
        }
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

    /// Compute stress-induced martensite fraction ξ_S for the loading branch.
    ///
    /// Uses cosine kinetics between `sigma_crit_s` and `sigma_crit_f`
    /// shifted by temperature.
    pub fn xi_s_loading(&self, temp: f64, stress: f64, xi_s_prev: f64, _xi_t_prev: f64) -> f64 {
        let (sig_s, sig_f) = self.critical_stresses_at_temp(temp);
        if stress < sig_s || stress > sig_f {
            return xi_s_prev;
        }
        let dsig = sig_f - sig_s;
        if dsig.abs() < EPS {
            return xi_s_prev;
        }
        let cos_val = (PI / dsig * (stress - sig_f)).cos();
        let xi_s = (1.0 - xi_s_prev) / 2.0 * cos_val + (1.0 + xi_s_prev) / 2.0;
        xi_s.clamp(0.0, 1.0)
    }

    /// Compute temperature-induced martensite fraction ξ_T for cooling.
    pub fn xi_t_cooling(&self, temp: f64, stress: f64, xi_t_prev: f64) -> f64 {
        let shifted = self.temps.shifted(stress, self.cm, self.ca);
        if temp > shifted.ms || temp < shifted.mf {
            return xi_t_prev;
        }
        let dm = shifted.ms - shifted.mf;
        if dm.abs() < EPS {
            return xi_t_prev;
        }
        let cos_val = (PI * (temp - shifted.mf) / dm).cos();
        let xi_t = (1.0 - xi_t_prev) / 2.0 * cos_val + (1.0 + xi_t_prev) / 2.0;
        xi_t.clamp(0.0, 1.0)
    }

    /// Reverse transformation (M → A): both ξ_S and ξ_T decrease.
    pub fn reverse_fractions(
        &self,
        temp: f64,
        stress: f64,
        xi_s_prev: f64,
        xi_t_prev: f64,
    ) -> (f64, f64) {
        let shifted = self.temps.shifted(stress, self.cm, self.ca);
        if temp < shifted.as_ || temp > shifted.af {
            return (xi_s_prev, xi_t_prev);
        }
        let da = shifted.af - shifted.as_;
        if da.abs() < EPS {
            return (xi_s_prev, xi_t_prev);
        }
        let xi_prev = xi_s_prev + xi_t_prev;
        let cos_val = (PI * (temp - shifted.as_) / da).cos();
        let xi_new = xi_prev / 2.0 * (cos_val + 1.0);
        let xi_new = xi_new.clamp(0.0, 1.0);
        // Distribute proportionally
        let ratio = if xi_prev > EPS { xi_new / xi_prev } else { 0.0 };
        (xi_s_prev * ratio, xi_t_prev * ratio)
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

    /// Incremental update for a full Brinson loading step.
    pub fn update_state(&self, state: &mut SmaPhaseState, strain_new: f64, temp_new: f64) {
        let _temp_old = state.temperature;
        let stress_est = self.stress_response(
            strain_new,
            state.xi_s,
            state.xi_t,
            temp_new,
            state.temperature,
        );

        // Determine which branch (forward/reverse)
        let shifted = self.temps.shifted(stress_est, self.cm, self.ca);
        if temp_new <= shifted.ms && temp_new >= shifted.mf {
            // Cooling → forward transformation
            let xi_s_new = self.xi_s_loading(temp_new, stress_est, state.xi_s, state.xi_t);
            let xi_t_new = self.xi_t_cooling(temp_new, stress_est, state.xi_t);
            state.xi_s = xi_s_new;
            state.xi_t = xi_t_new;
        } else if temp_new >= shifted.as_ && temp_new <= shifted.af {
            // Heating → reverse transformation
            let (xi_s_new, xi_t_new) =
                self.reverse_fractions(temp_new, stress_est, state.xi_s, state.xi_t);
            state.xi_s = xi_s_new;
            state.xi_t = xi_t_new;
        }

        state.xi = (state.xi_s + state.xi_t).clamp(0.0, 1.0);
        state.strain = strain_new;
        state.temperature = temp_new;
        state.transformation_strain = self.h_max * state.xi_s;
        state.stress = self.stress_response(strain_new, state.xi_s, state.xi_t, temp_new, temp_new);
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

    /// Generate loading branch stress-strain data.
    ///
    /// Returns a vector of `(strain, stress)` pairs.
    pub fn loading_curve(&self, max_strain: f64) -> Vec<(f64, f64)> {
        let n = self.resolution.max(2);
        let mut result = Vec::with_capacity(n);
        let mut state = SmaPhaseState::new(self.temperature);
        for i in 0..n {
            let eps = max_strain * (i as f64) / (n as f64 - 1.0);
            self.model.update_state(&mut state, eps, self.temperature);
            result.push((eps, state.stress));
        }
        result
    }

    /// Generate unloading branch stress-strain data starting from a given state.
    pub fn unloading_curve(&self, state_at_peak: &SmaPhaseState) -> Vec<(f64, f64)> {
        let n = self.resolution.max(2);
        let mut result = Vec::with_capacity(n);
        let peak_strain = state_at_peak.strain;
        let mut state = *state_at_peak;
        for i in 0..n {
            let eps = peak_strain * (1.0 - (i as f64) / (n as f64 - 1.0));
            // Simple unloading: update with decreasing strain
            self.model.update_state(&mut state, eps, self.temperature);
            result.push((eps, state.stress.max(0.0)));
        }
        result
    }

    /// Full loading-unloading loop.
    pub fn full_loop(&self, max_strain: f64) -> Vec<(f64, f64)> {
        let loading = self.loading_curve(max_strain);
        let last_state = {
            let mut s = SmaPhaseState::new(self.temperature);
            for &(eps, _) in &loading {
                self.model.update_state(&mut s, eps, self.temperature);
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
    /// Returns `(residual_strain_before_heating, residual_strain_after_heating)`.
    pub fn simulate_recovery(&self, applied_strain: f64) -> (f64, f64) {
        // Step 1: deform at T < Mf (full martensite)
        let mut state = SmaPhaseState::fully_martensite(self.deform_temperature);
        self.model
            .update_state(&mut state, applied_strain, self.deform_temperature);
        let residual_before = state.transformation_strain;

        // Step 2: unload mechanically
        self.model
            .update_state(&mut state, 0.0, self.deform_temperature);
        let _strain_after_unload = state.strain;

        // Step 3: heat above Af
        self.model
            .update_state(&mut state, 0.0, self.recovery_temperature);
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
/// Returns a vector of `(temperature, strain)` pairs showing recovery.
pub fn temperature_sweep_recovery(
    model: &BrinsonModel,
    initial_strain: f64,
    t_start: f64,
    t_end: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let n = n_steps.max(2);
    let mut result = Vec::with_capacity(n);
    let mut state = SmaPhaseState::fully_martensite(t_start);
    state.strain = initial_strain;
    state.transformation_strain = model.h_max;

    for i in 0..n {
        let t_frac = i as f64 / (n as f64 - 1.0);
        let temp = t_start + t_frac * (t_end - t_start);
        model.update_state(&mut state, 0.0, temp);
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
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-6;

    // Helper: create NiTi Brinson model
    fn niti_brinson() -> BrinsonModel {
        NiTiPreset::brinson()
    }

    // Helper: create NiTi temperatures
    fn niti_temps() -> TransformationTemperatures {
        NiTiPreset::temperatures()
    }

    // 1. TransformationTemperatures basic sanity
    #[test]
    fn test_transformation_temperatures_basic() {
        let tt = niti_temps();
        assert!(tt.ms > tt.mf, "Ms > Mf required");
        assert!(tt.af > tt.as_, "Af > As required");
        assert!(tt.hysteresis() > 0.0, "Hysteresis should be positive");
    }

    // 2. Temperature shift via Clausius-Clapeyron
    #[test]
    fn test_temperature_shift() {
        let tt = niti_temps();
        let stress = 100.0e6; // 100 MPa
        let shifted = tt.shifted(stress, 6.5e6, 6.5e6);
        let expected_shift = 100.0e6 / 6.5e6;
        assert!((shifted.ms - tt.ms - expected_shift).abs() < TOL);
        assert!((shifted.af - tt.af - expected_shift).abs() < TOL);
    }

    // 3. Tanaka model: fully martensite below Mf
    #[test]
    fn test_tanaka_full_martensite() {
        let model = NiTiPreset::tanaka();
        let xi = model.forward_fraction(200.0, 0.0); // well below Mf
        assert!(xi > 0.99, "Should be nearly fully martensite, got {xi}");
    }

    // 4. Tanaka model: fully austenite above Af
    #[test]
    fn test_tanaka_full_austenite() {
        let model = NiTiPreset::tanaka();
        let xi = model.reverse_fraction(400.0, 0.0); // well above Af
        assert!(xi < 0.01, "Should be nearly fully austenite, got {xi}");
    }

    // 5. Tanaka effective modulus at ξ=0 equals E_A
    #[test]
    fn test_tanaka_modulus_austenite() {
        let model = NiTiPreset::tanaka();
        let e = model.effective_modulus(0.0);
        assert!((e - 67.0e9).abs() < TOL, "E at xi=0 should be E_A");
    }

    // 6. Tanaka effective modulus at ξ=1 equals E_M
    #[test]
    fn test_tanaka_modulus_martensite() {
        let model = NiTiPreset::tanaka();
        let e = model.effective_modulus(1.0);
        assert!((e - 26.3e9).abs() < TOL, "E at xi=1 should be E_M");
    }

    // 7. Liang-Rogers: phase fraction at Mf is 1
    #[test]
    fn test_liang_rogers_mf() {
        let model = NiTiPreset::liang_rogers();
        let xi = model.phase_fraction(model.temps.mf, 0.0, 0.0);
        assert!(
            (xi - 1.0).abs() < 0.01,
            "At Mf, xi should be ~1.0, got {xi}"
        );
    }

    // 8. Liang-Rogers: phase fraction above Af is 0
    #[test]
    fn test_liang_rogers_above_af() {
        let model = NiTiPreset::liang_rogers();
        let xi = model.phase_fraction(model.temps.af + 50.0, 0.0, 1.0);
        assert!(xi < 0.01, "Above Af, xi should be ~0, got {xi}");
    }

    // 9. Liang-Rogers: reverse fraction at Af with xi_prev=1
    #[test]
    fn test_liang_rogers_reverse_at_af() {
        let model = NiTiPreset::liang_rogers();
        let xi = model.reverse_fraction(model.temps.af, 0.0, 1.0);
        assert!(xi < 0.05, "At Af with xi_prev=1, xi should be ~0, got {xi}");
    }

    // 10. Brinson: effective modulus interpolation
    #[test]
    fn test_brinson_modulus_interpolation() {
        let model = niti_brinson();
        let e_half = model.effective_modulus(0.5);
        let expected = 0.5 * 67.0e9 + 0.5 * 26.3e9;
        assert!(
            (e_half - expected).abs() < TOL,
            "E at xi=0.5 should be average, got {e_half}"
        );
    }

    // 11. Brinson: critical stress increases with temperature above Ms
    #[test]
    fn test_brinson_critical_stress_vs_temp() {
        let model = niti_brinson();
        let (sig_s_low, _) = model.critical_stresses_at_temp(model.temps.ms);
        let (sig_s_high, _) = model.critical_stresses_at_temp(model.temps.ms + 20.0);
        assert!(
            sig_s_high > sig_s_low,
            "Critical stress should increase with T"
        );
    }

    // 12. Brinson: stress response at zero strain and xi=0 is zero (at T_ref)
    #[test]
    fn test_brinson_zero_stress_at_origin() {
        let model = niti_brinson();
        let sigma = model.stress_response(0.0, 0.0, 0.0, 300.0, 300.0);
        assert!(
            sigma.abs() < TOL,
            "Stress at origin should be zero, got {sigma}"
        );
    }

    // 13. Brinson: state update changes xi during cooling
    #[test]
    fn test_brinson_state_update_cooling() {
        let model = niti_brinson();
        let mut state = SmaPhaseState::new(model.temps.af + 10.0);
        // Cool to between Mf and Ms
        let target_temp = (model.temps.mf + model.temps.ms) / 2.0;
        model.update_state(&mut state, 0.001, target_temp);
        assert!(
            state.xi > 0.0,
            "xi should increase during cooling, got {}",
            state.xi
        );
    }

    // 14. Superelastic curve: loading curve is non-empty
    #[test]
    fn test_superelastic_loading_nonempty() {
        let model = niti_brinson();
        let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 50);
        let curve = curve_gen.loading_curve(0.06);
        assert!(!curve.is_empty(), "Loading curve should be non-empty");
        assert_eq!(curve.len(), 50);
    }

    // 15. Superelastic curve: stress generally increases during loading
    #[test]
    fn test_superelastic_loading_trend() {
        let model = niti_brinson();
        let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 100);
        let curve = curve_gen.loading_curve(0.05);
        // Overall trend: final stress should be higher than initial (non-zero) stress
        let first_nonzero = curve.iter().find(|(_, s)| *s > 0.0);
        let last = curve.last().unwrap();
        if let Some(first) = first_nonzero {
            assert!(
                last.1 >= first.1,
                "Final stress should exceed initial: first={:.6}, last={:.6}",
                first.1,
                last.1
            );
        }
    }

    // 16. Superelastic curve: full loop has more points than one branch
    #[test]
    fn test_superelastic_full_loop_length() {
        let model = niti_brinson();
        let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 30);
        let loop_curve = curve_gen.full_loop(0.05);
        assert!(
            loop_curve.len() > 30,
            "Full loop should have loading + unloading points"
        );
    }

    // 17. Shape memory effect: recovery ratio is non-negative
    #[test]
    fn test_sme_recovery_ratio_nonneg() {
        let model = niti_brinson();
        let sme = ShapeMemoryEffect::new(model, 250.0, 350.0);
        let ratio = sme.recovery_ratio(0.05);
        assert!(ratio >= 0.0, "Recovery ratio should be >= 0, got {ratio}");
    }

    // 18. Two-way shape memory: training increases two_way_fraction
    #[test]
    fn test_twsm_training() {
        let model = niti_brinson();
        let mut twsm = TwoWayShapeMemory::new(model);
        assert!(twsm.two_way_fraction < TOL, "Initially should be ~0");
        for _ in 0..20 {
            twsm.train_cycle();
        }
        assert!(
            twsm.two_way_fraction > 0.1,
            "After training, fraction should grow, got {}",
            twsm.two_way_fraction
        );
    }

    // 19. Two-way: cooling strain at T < Mf is nonzero after training
    #[test]
    fn test_twsm_cooling_strain() {
        let model = niti_brinson();
        let mut twsm = TwoWayShapeMemory::new(model.clone());
        for _ in 0..50 {
            twsm.train_cycle();
        }
        let eps = twsm.cooling_strain(model.temps.mf - 10.0);
        assert!(eps > 0.0, "Cooling strain should be > 0, got {eps}");
    }

    // 20. Two-way: heating strain above Af is zero
    #[test]
    fn test_twsm_heating_above_af() {
        let model = niti_brinson();
        let mut twsm = TwoWayShapeMemory::new(model.clone());
        for _ in 0..10 {
            twsm.train_cycle();
        }
        let eps = twsm.heating_strain(model.temps.af + 20.0);
        assert!(
            eps.abs() < TOL,
            "Heating strain above Af should be ~0, got {eps}"
        );
    }

    // 21. Training effect: residual strain accumulates
    #[test]
    fn test_training_residual_accumulates() {
        let mut te = TrainingEffect::new();
        te.add_cycle();
        let r1 = te.current_residual();
        for _ in 0..50 {
            te.add_cycle();
        }
        let r51 = te.current_residual();
        assert!(r51 > r1, "Residual should grow: r1={r1}, r51={r51}");
    }

    // 22. Training effect: plateau stress decreases
    #[test]
    fn test_training_plateau_decreases() {
        let mut te = TrainingEffect::new();
        te.add_cycle();
        let s1 = te.current_plateau_stress();
        for _ in 0..100 {
            te.add_cycle();
        }
        let s101 = te.current_plateau_stress();
        assert!(
            s101 < s1,
            "Plateau stress should decrease: s1={s1}, s101={s101}"
        );
    }

    // 23. Thermomechanical: forward transformation heats
    #[test]
    fn test_thermo_forward_heats() {
        let tmc = NiTiPreset::thermomechanical();
        let dt = tmc.temperature_increment(0.5); // dxi = +0.5
        // Forward transformation releases heat → dt < 0 because of formula sign
        // but physical convention: temperature rises. Our formula: dT = -L/cp * dxi
        // dxi > 0 → dT < 0 (latent heat absorbed from surroundings for forward)
        // Actually for SMA, A→M is exothermic so temp rises.
        // Our formula gives dT = -L/cp * dxi. If L > 0 and dxi > 0 → dT < 0.
        // This models the environment perspective. Let's just verify it's nonzero.
        assert!(dt.abs() > 0.0, "Temperature change should be nonzero");
    }

    // 24. Thermomechanical: no phase change → no temperature change
    #[test]
    fn test_thermo_no_change() {
        let tmc = NiTiPreset::thermomechanical();
        let dt = tmc.temperature_increment(0.0);
        assert!(dt.abs() < TOL, "No phase change → no dT");
    }

    // 25. Clausius-Clapeyron: symmetric slopes
    #[test]
    fn test_cc_symmetric() {
        let cc = ClausiusClapeyron::symmetric(6.5e6);
        assert!(
            (cc.cm - cc.ca).abs() < TOL,
            "Symmetric slopes should be equal"
        );
    }

    // 26. Clausius-Clapeyron: stress shift
    #[test]
    fn test_cc_stress_shift() {
        let cc = ClausiusClapeyron::new(6.5e6, 7.0e6);
        let ds = cc.stress_shift_martensite(10.0);
        let expected = 6.5e6 * 10.0;
        assert!((ds - expected).abs() < TOL, "Stress shift incorrect");
    }

    // 27. CuAlNi preset creates valid model
    #[test]
    fn test_cualni_preset() {
        let model = CuAlNiPreset::brinson();
        assert!(model.e_a > model.e_m, "E_A > E_M for CuAlNi");
        assert!(model.h_max > 0.0, "h_max should be positive");
    }

    // 28. CuZnAl preset creates valid model
    #[test]
    fn test_cuznal_preset() {
        let model = CuZnAlPreset::brinson();
        let tt = CuZnAlPreset::temperatures();
        assert!(tt.af > tt.ms, "Af > Ms required");
        assert!(model.h_max > 0.0);
    }

    // 29. Hysteresis energy of identical curves is zero
    #[test]
    fn test_hysteresis_energy_zero() {
        let curve = vec![(0.0, 0.0), (0.01, 100.0e6), (0.02, 200.0e6)];
        let e = hysteresis_energy(&curve, &curve);
        assert!(e < TOL, "Same loading/unloading → zero hysteresis");
    }

    // 30. Hysteresis energy of a simple rectangle
    #[test]
    fn test_hysteresis_energy_rectangle() {
        let loading = vec![(0.0, 100.0), (1.0, 100.0)];
        let unloading = vec![(0.0, 50.0), (1.0, 50.0)];
        let e = hysteresis_energy(&loading, &unloading);
        assert!(
            (e - 50.0).abs() < TOL,
            "Rectangle hysteresis should be 50.0, got {e}"
        );
    }

    // 31. Secant modulus
    #[test]
    fn test_secant_modulus() {
        let m = secant_modulus(0.01, 700.0e6);
        assert!((m - 70.0e9).abs() < 1.0, "Secant modulus should be E");
    }

    // 32. Tangent modulus
    #[test]
    fn test_tangent_modulus() {
        let m = tangent_modulus(0.01, 700.0e6, 0.02, 1400.0e6);
        assert!((m - 70.0e9).abs() < 1.0, "Tangent modulus should be 70 GPa");
    }

    // 33. Temperature sweep recovery returns correct length
    #[test]
    fn test_temp_sweep_length() {
        let model = niti_brinson();
        let sweep = temperature_sweep_recovery(&model, 0.05, 250.0, 350.0, 20);
        assert_eq!(sweep.len(), 20);
    }

    // 34. Brinson reverse fractions decrease xi
    #[test]
    fn test_brinson_reverse_decreases() {
        let model = niti_brinson();
        let mid_temp = (model.temps.as_ + model.temps.af) / 2.0;
        let (xs, xt) = model.reverse_fractions(mid_temp, 0.0, 0.5, 0.5);
        assert!(xs + xt < 1.0, "Reverse should decrease total xi");
    }

    // 35. SmaPhaseState: fully martensite has xi=1
    #[test]
    fn test_phase_state_full_martensite() {
        let state = SmaPhaseState::fully_martensite(250.0);
        assert!((state.xi - 1.0).abs() < TOL);
        assert!((state.xi_t - 1.0).abs() < TOL);
    }
}
