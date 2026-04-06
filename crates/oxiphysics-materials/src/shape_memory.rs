// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shape memory alloy (SMA) materials — Nitinol, Cu-Zn-Al, and related alloys.
//!
//! Provides constitutive models for shape memory alloys based on the
//! Brinson–Liang–Rogers thermomechanical framework, Clausius-Clapeyron
//! thermoelastic martensite theory, pseudoelastic loading/unloading, and
//! two-way shape memory effect.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// ShapeMemoryAlloy
// ─────────────────────────────────────────────────────────────────────────────

/// Thermomechanical model of a shape memory alloy.
///
/// Transformation temperatures are in kelvin (K); stresses in pascal (Pa);
/// elastic moduli in pascal (Pa); strains dimensionless.
///
/// The martensite fraction ξ ∈ \[0, 1\] describes the phase state:
/// * ξ = 1 → fully martensitic
/// * ξ = 0 → fully austenitic
#[derive(Debug, Clone)]
pub struct ShapeMemoryAlloy {
    /// Martensite start temperature (Ms) in K.
    pub ms: f64,
    /// Martensite finish temperature (Mf) in K.
    pub mf: f64,
    /// Austenite start temperature (As) in K.
    pub as_: f64,
    /// Austenite finish temperature (Af) in K.
    pub af: f64,
    /// Elastic modulus of the austenite phase in Pa.
    pub e_austenite: f64,
    /// Elastic modulus of the martensite phase in Pa.
    pub e_martensite: f64,
    /// Maximum recoverable (transformation) strain ε_L (dimensionless).
    pub h_max: f64,
}

impl ShapeMemoryAlloy {
    /// Create a new SMA model with the given transformation temperatures and
    /// elastic moduli.
    ///
    /// # Parameters
    /// * `ms` – martensite start temperature \[K\]
    /// * `mf` – martensite finish temperature \[K\]
    /// * `as_` – austenite start temperature \[K\]
    /// * `af` – austenite finish temperature \[K\]
    /// * `e_austenite` – austenite elastic modulus \[Pa\]
    /// * `e_martensite` – martensite elastic modulus \[Pa\]
    /// * `h_max` – maximum transformation strain
    pub fn new(
        ms: f64,
        mf: f64,
        as_: f64,
        af: f64,
        e_austenite: f64,
        e_martensite: f64,
        h_max: f64,
    ) -> Self {
        Self {
            ms,
            mf,
            as_,
            af,
            e_austenite,
            e_martensite,
            h_max,
        }
    }

    /// Create a standard Nitinol (NiTi) SMA with typical reported parameters.
    ///
    /// * Ms = 291 K, Mf = 273 K, As = 307 K, Af = 325 K
    /// * E_A = 75 GPa, E_M = 28 GPa, ε_L = 0.08
    pub fn new_nitinol() -> Self {
        Self::new(
            291.0,  // Ms
            273.0,  // Mf
            307.0,  // As
            325.0,  // Af
            75.0e9, // E_austenite
            28.0e9, // E_martensite
            0.08,   // h_max
        )
    }

    /// Alias for [`new_nitinol`](Self::new_nitinol).
    pub fn nitinol() -> Self {
        Self::new_nitinol()
    }

    /// Martensite volume fraction ξ for the given temperature and stress.
    ///
    /// Uses the Liang-Rogers cosine model:
    /// * Cooling (martensitic): ξ = 0.5 · cos(π · (T − Mf)/(Ms − Mf)) + 0.5
    /// * Heating (austenitic): ξ = 0.5 · cos(π · (T − As)/(Af − As)) + 0.5  (then 1 − value)
    ///
    /// The stress shifts the transformation temperatures via the
    /// Clausius-Clapeyron slope (~10 MPa/K for NiTi).
    pub fn phase_fraction(&self, temp: f64, stress: f64) -> f64 {
        // Clausius-Clapeyron stress correction (MPa/K → K shift)
        let cc_slope = 10.0e6; // Pa/K
        let dt = stress / cc_slope;

        let ms = self.ms + dt;
        let mf = self.mf + dt;
        let as_ = self.as_ + dt;
        let af = self.af + dt;

        if temp <= mf {
            // Fully martensitic
            1.0_f64
        } else if temp <= ms {
            // Martensitic transformation
            0.5 * (PI * (temp - mf) / (ms - mf)).cos() + 0.5
        } else if temp < as_ {
            // Stable mixed region between Ms and As
            1.0_f64
        } else if temp <= af {
            // Austenitic transformation
            1.0 - (0.5 * (PI * (temp - as_) / (af - as_)).cos() + 0.5)
        } else {
            // Fully austenitic
            0.0_f64
        }
    }

    /// Effective elastic modulus at a given martensite fraction ξ.
    ///
    /// Uses the rule of mixtures: `E(ξ) = E_A + ξ · (E_M − E_A)`.
    pub fn elastic_modulus(&self, xi: f64) -> f64 {
        self.e_austenite + xi.clamp(0.0, 1.0) * (self.e_martensite - self.e_austenite)
    }

    /// One-dimensional constitutive stress response (Pa) for a given mechanical
    /// strain and temperature.
    ///
    /// `σ = E(ξ) · ε − E(ξ) · ε_L · ξ`
    pub fn constitutive_response(&self, strain: f64, temp: f64) -> f64 {
        let xi = self.phase_fraction(temp, 0.0);
        let e = self.elastic_modulus(xi);
        e * strain - e * self.h_max * xi
    }

    /// Recoverable strain between two martensite fraction states.
    ///
    /// `Δε = h_max · (ξ_start − ξ_end)`
    pub fn recovery_strain(&self, xi_start: f64, xi_end: f64) -> f64 {
        self.h_max * (xi_start - xi_end)
    }

    /// Critical transformation stress at a given temperature.
    ///
    /// Uses the Clausius-Clapeyron relation:
    /// `σ_cr = C_M · (T − Ms)` for T > Ms (stress-induced martensite).
    /// Returns 0 for T ≤ Ms.
    pub fn critical_stress(&self, temp: f64) -> f64 {
        let cm = 10.0e6; // Pa/K, typical NiTi value
        if temp > self.ms {
            cm * (temp - self.ms)
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BraisbirdAuricchioModel
// ─────────────────────────────────────────────────────────────────────────────

/// Brinson–Auricchio 3D SMA constitutive model parameters.
///
/// Encodes the forward and reverse transformation stress slopes
/// (Ca for austenite, Cm for martensite) in Pa/K.
#[derive(Debug, Clone)]
pub struct BraisbirdAuricchioModel {
    /// Stress influence coefficient for austenite transformation in Pa/K.
    pub ca: f64,
    /// Stress influence coefficient for martensite transformation in Pa/K.
    pub cm: f64,
    /// Reference temperature above which austenite is stable \[K\].
    pub t0: f64,
    /// Plateau stress offset at the reference temperature \[Pa\].
    pub sigma0: f64,
}

impl BraisbirdAuricchioModel {
    /// Create a new Brinson-Auricchio model.
    pub fn new(ca: f64, cm: f64, t0: f64, sigma0: f64) -> Self {
        Self { ca, cm, t0, sigma0 }
    }

    /// Standard NiTi Brinson-Auricchio parameters.
    pub fn nitinol() -> Self {
        Self::new(13.0e6, 8.0e6, 291.0, 100.0e6)
    }

    /// Forward transformation (austenite → martensite) critical stress.
    ///
    /// `σ_fwd = σ₀ + Cm · (T − T₀)`
    pub fn forward_transformation_stress(&self) -> f64 {
        self.sigma0
    }

    /// Reverse transformation (martensite → austenite) critical stress.
    ///
    /// The reverse (austenite) plateau stress is lower than the forward
    /// plateau.  It is estimated as `σ₀ · Cm / Ca` — when Ca > Cm (as in
    /// typical NiTi), this gives a value smaller than `σ₀`.
    pub fn reverse_transformation_stress(&self) -> f64 {
        if self.ca.abs() < 1e-15 {
            return 0.0;
        }
        self.sigma0 * (self.cm / self.ca)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermoelasticMartensite
// ─────────────────────────────────────────────────────────────────────────────

/// Thermoelastic martensite model for Clausius-Clapeyron analysis.
///
/// Characterises the thermodynamic driving force for martensitic
/// transformation in terms of the transformation strain and entropy change.
#[derive(Debug, Clone)]
pub struct ThermoelasticMartensite {
    /// Transformation strain ε_L (dimensionless).
    pub epsilon_l: f64,
    /// Entropy change per unit volume ΔS \[J/(m³·K)\].
    pub delta_s: f64,
}

impl ThermoelasticMartensite {
    /// Create a new thermoelastic martensite model.
    pub fn new(epsilon_l: f64, delta_s: f64) -> Self {
        Self { epsilon_l, delta_s }
    }

    /// Standard NiTi thermoelastic martensite parameters.
    pub fn nitinol() -> Self {
        // ε_L = 0.08; ΔS ≈ −0.22 J/(cm³·K) ≈ −2.2e5 J/(m³·K)
        Self::new(0.08, -2.2e5)
    }

    /// Clausius-Clapeyron slope dσ/dT \[Pa/K\].
    ///
    /// From thermodynamic equilibrium:
    /// `dσ/dT = −ΔS / ε_L`
    pub fn clausius_clapeyron(&self) -> f64 {
        if self.epsilon_l.abs() < 1e-15 {
            return 0.0;
        }
        -self.delta_s / self.epsilon_l
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SMAPseudoelasticity
// ─────────────────────────────────────────────────────────────────────────────

/// Pseudoelastic (superelastic) behaviour of a shape memory alloy.
///
/// Describes the stress-strain loop with forward and reverse transformation
/// plateau stresses and the resulting hysteresis.
#[derive(Debug, Clone)]
pub struct SMAPseudoelasticity {
    /// Upper plateau (loading) stress in Pa.
    pub sigma_f: f64,
    /// Lower plateau (unloading) stress in Pa.
    pub sigma_r: f64,
    /// Maximum transformation strain ε_L.
    pub epsilon_l: f64,
    /// Elastic modulus (austenite) in Pa.
    pub e_a: f64,
}

impl SMAPseudoelasticity {
    /// Create a pseudoelastic model.
    pub fn new(sigma_f: f64, sigma_r: f64, epsilon_l: f64, e_a: f64) -> Self {
        Self {
            sigma_f,
            sigma_r,
            epsilon_l,
            e_a,
        }
    }

    /// Standard Nitinol pseudoelastic parameters at 310 K.
    pub fn nitinol_room_temp() -> Self {
        Self::new(500.0e6, 200.0e6, 0.06, 75.0e9)
    }

    /// Forward transformation (loading) plateau stress in Pa.
    pub fn loading_plateau_stress(&self) -> f64 {
        self.sigma_f
    }

    /// Reverse transformation (unloading) plateau stress in Pa.
    pub fn unloading_plateau_stress(&self) -> f64 {
        self.sigma_r
    }

    /// Mechanical energy dissipated per cycle (hysteresis area) in J/m³.
    ///
    /// Area of the stress-strain hysteresis loop:
    /// `W = (σ_f − σ_r) · ε_L`
    pub fn hysteresis_area(&self) -> f64 {
        (self.sigma_f - self.sigma_r) * self.epsilon_l
    }

    /// Stress-strain response during loading at fractional strain `eps / eps_L`.
    ///
    /// Returns stress in Pa.  Assumes linear elastic loading up to σ_f,
    /// then a flat plateau at σ_f until full transformation.
    pub fn loading_response(&self, eps: f64) -> f64 {
        let eps_onset = self.sigma_f / self.e_a;
        if eps <= eps_onset {
            self.e_a * eps
        } else {
            self.sigma_f
        }
    }

    /// Stress-strain response during unloading from full transformation.
    ///
    /// Returns stress in Pa.
    pub fn unloading_response(&self, eps: f64) -> f64 {
        // At full transformation strain eps_L the stress is sigma_f.
        // On unloading, elastic until sigma_r, then reverse plateau.
        let eps_end = self.epsilon_l;
        let eps_reverse_end = (self.sigma_f - self.sigma_r) / self.e_a;
        let eps_elastic_end = eps_end - eps_reverse_end;
        if eps >= eps_elastic_end {
            // Elastic unloading from plateau
            self.sigma_f - self.e_a * (eps_end - eps)
        } else {
            // Reverse transformation plateau
            self.sigma_r.max(0.0)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TwoWaySMA
// ─────────────────────────────────────────────────────────────────────────────

/// Two-way shape memory effect (TWSME) model.
///
/// After repeated thermomechanical training cycles, SMAs develop an intrinsic
/// martensite preferred orientation that allows them to actuate bidirectionally
/// without external stress.
#[derive(Debug, Clone)]
pub struct TwoWaySMA {
    /// Number of completed thermomechanical training cycles.
    pub training_cycles: usize,
    /// Maximum trained strain achievable (dimensionless).
    pub max_trained_strain: f64,
}

impl TwoWaySMA {
    /// Create a new two-way SMA model.
    pub fn new(training_cycles: usize, max_trained_strain: f64) -> Self {
        Self {
            training_cycles,
            max_trained_strain,
        }
    }

    /// Standard Nitinol two-way model.
    pub fn nitinol() -> Self {
        Self::new(0, 0.04) // typical max TWSME strain ~4%
    }

    /// Add training cycles and return new total.
    pub fn train(&mut self, cycles: usize) {
        self.training_cycles += cycles;
    }

    /// Recoverable trained strain after `training_cycles` cycles.
    ///
    /// The TWSME strain saturates with an empirical logarithmic law:
    /// `ε_trained = ε_max · (1 − exp(−n / 10))`
    /// where n is the number of training cycles.
    pub fn trained_strain(&self) -> f64 {
        let n = self.training_cycles as f64;
        self.max_trained_strain * (1.0 - (-n / 10.0).exp())
    }

    /// Fraction of maximum trained strain achieved.
    pub fn saturation_fraction(&self) -> f64 {
        if self.max_trained_strain.abs() < 1e-15 {
            return 0.0;
        }
        self.trained_strain() / self.max_trained_strain
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Nitinol phase diagram
// ─────────────────────────────────────────────────────────────────────────────

/// Return the four characteristic transformation temperatures for NiTi Nitinol.
///
/// The return value is `[(label_temp, stress_free_temp); 4]` for
/// `[Ms, Mf, As, Af]` in kelvin at zero stress.
///
/// Values are typical for equiatomic Nitinol near room temperature.
pub fn nitinol_phase_diagram() -> [(f64, f64); 4] {
    [
        (1.0, 291.0), // Ms  — martensite start
        (2.0, 273.0), // Mf  — martensite finish
        (3.0, 307.0), // As  — austenite start
        (4.0, 325.0), // Af  — austenite finish
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: smoothstep between two phases
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth cosine interpolation used in transformation models, ∈ \[0, 1\].
fn cosine_interpolate(x: f64) -> f64 {
    0.5 * (1.0 - (PI * x).cos())
}

/// Effective thermal conductivity of an SMA biphasic mixture \[W/(m·K)\].
///
/// Uses Voigt (rule of mixtures) averaging:
/// `k = ξ · k_M + (1 − ξ) · k_A`
pub fn thermal_conductivity(xi: f64, k_martensite: f64, k_austenite: f64) -> f64 {
    let xi = xi.clamp(0.0, 1.0);
    xi * k_martensite + (1.0 - xi) * k_austenite
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ShapeMemoryAlloy ─────────────────────────────────────────────────────

    #[test]
    fn test_nitinol_temperatures() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        assert!(sma.ms > sma.mf, "Ms must be above Mf");
        assert!(sma.af > sma.as_, "Af must be above As");
        assert!(sma.as_ > sma.ms, "As must be above Ms for NiTi");
    }

    #[test]
    fn test_phase_fraction_fully_martensitic() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        // Below Mf: fully martensitic
        let xi = sma.phase_fraction(250.0, 0.0);
        assert!(
            (xi - 1.0).abs() < 1e-10,
            "should be xi=1 below Mf, got {xi}"
        );
    }

    #[test]
    fn test_phase_fraction_fully_austenitic() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        // Above Af: fully austenitic
        let xi = sma.phase_fraction(340.0, 0.0);
        assert!(xi < 0.01, "should be xi≈0 above Af, got {xi}");
    }

    #[test]
    fn test_phase_fraction_transition_region() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        // Between Mf and Ms: partially transformed
        let xi = sma.phase_fraction(282.0, 0.0);
        assert!(
            xi > 0.0 && xi < 1.0,
            "xi should be in (0,1) in transformation region, got {xi}"
        );
    }

    #[test]
    fn test_elastic_modulus_pure_austenite() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let e = sma.elastic_modulus(0.0);
        assert!(
            (e - sma.e_austenite).abs() < 1.0,
            "pure austenite modulus wrong"
        );
    }

    #[test]
    fn test_elastic_modulus_pure_martensite() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let e = sma.elastic_modulus(1.0);
        assert!(
            (e - sma.e_martensite).abs() < 1.0,
            "pure martensite modulus wrong"
        );
    }

    #[test]
    fn test_elastic_modulus_mixed() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let e = sma.elastic_modulus(0.5);
        let expected = 0.5 * (sma.e_austenite + sma.e_martensite);
        assert!((e - expected).abs() < 1.0, "mixed modulus wrong");
    }

    #[test]
    fn test_constitutive_response_sign() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        // Positive strain at high temperature (austenitic) should give positive stress
        let sigma = sma.constitutive_response(0.01, 350.0);
        assert!(
            sigma > 0.0,
            "positive strain should give positive stress in austenite"
        );
    }

    #[test]
    fn test_recovery_strain_positive() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let eps = sma.recovery_strain(1.0, 0.0);
        assert!(
            (eps - sma.h_max).abs() < 1e-10,
            "full phase recovery should give h_max"
        );
    }

    #[test]
    fn test_recovery_strain_zero_change() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let eps = sma.recovery_strain(0.5, 0.5);
        assert!(eps.abs() < 1e-10, "no phase change → zero recovery strain");
    }

    #[test]
    fn test_critical_stress_above_ms() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        // At T > Ms the critical stress should be positive
        let sigma = sma.critical_stress(310.0);
        assert!(sigma > 0.0, "critical stress should be positive above Ms");
    }

    #[test]
    fn test_critical_stress_below_ms() {
        let sma = ShapeMemoryAlloy::new_nitinol();
        let sigma = sma.critical_stress(250.0);
        assert!(sigma == 0.0, "critical stress should be zero below Ms");
    }

    // ── BraisbirdAuricchioModel ───────────────────────────────────────────────

    #[test]
    fn test_braisbird_forward_stress() {
        let m = BraisbirdAuricchioModel::nitinol();
        let sf = m.forward_transformation_stress();
        assert!(sf > 0.0, "forward transformation stress must be positive");
    }

    #[test]
    fn test_braisbird_reverse_stress() {
        let m = BraisbirdAuricchioModel::nitinol();
        let sr = m.reverse_transformation_stress();
        let sf = m.forward_transformation_stress();
        // Reverse stress should be less than forward for typical SMAs
        assert!(sr <= sf, "reverse stress should not exceed forward stress");
    }

    #[test]
    fn test_braisbird_nitinol_ca_cm() {
        let m = BraisbirdAuricchioModel::nitinol();
        assert!(m.ca > 0.0);
        assert!(m.cm > 0.0);
    }

    // ── ThermoelasticMartensite ───────────────────────────────────────────────

    #[test]
    fn test_clausius_clapeyron_sign() {
        let tm = ThermoelasticMartensite::nitinol();
        let ds_dt = tm.clausius_clapeyron();
        // For NiTi, dσ/dT > 0 (stress increases with temperature above Ms)
        assert!(
            ds_dt > 0.0,
            "Clausius-Clapeyron slope should be positive for NiTi, got {ds_dt}"
        );
    }

    #[test]
    fn test_clausius_clapeyron_magnitude() {
        let tm = ThermoelasticMartensite::nitinol();
        let ds_dt = tm.clausius_clapeyron();
        // Typical NiTi: ~6-8 MPa/K
        assert!(
            ds_dt > 1.0e6 && ds_dt < 50.0e6,
            "CC slope should be in the MPa/K range, got {ds_dt}"
        );
    }

    #[test]
    fn test_clausius_clapeyron_zero_strain() {
        let tm = ThermoelasticMartensite::new(0.0, -2.2e5);
        assert_eq!(tm.clausius_clapeyron(), 0.0);
    }

    // ── SMAPseudoelasticity ───────────────────────────────────────────────────

    #[test]
    fn test_pseudoelastic_loading_plateau() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        assert!((pe.loading_plateau_stress() - 500.0e6).abs() < 1.0);
    }

    #[test]
    fn test_pseudoelastic_unloading_plateau() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        assert!((pe.unloading_plateau_stress() - 200.0e6).abs() < 1.0);
    }

    #[test]
    fn test_pseudoelastic_hysteresis_positive() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        let area = pe.hysteresis_area();
        assert!(area > 0.0, "hysteresis area should be positive");
    }

    #[test]
    fn test_pseudoelastic_hysteresis_magnitude() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        // (500e6 - 200e6) * 0.06 = 18e6 J/m³
        let expected = (500.0e6 - 200.0e6) * 0.06;
        assert!((pe.hysteresis_area() - expected).abs() < 1.0);
    }

    #[test]
    fn test_pseudoelastic_loading_elastic() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        // Small strain: elastic regime
        let eps = 0.001;
        let sigma = pe.loading_response(eps);
        let expected = pe.e_a * eps;
        assert!((sigma - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn test_pseudoelastic_loading_plateau_region() {
        let pe = SMAPseudoelasticity::nitinol_room_temp();
        // Large strain: plateau
        let eps = 0.05;
        let sigma = pe.loading_response(eps);
        assert!((sigma - pe.sigma_f).abs() < 1.0);
    }

    // ── TwoWaySMA ────────────────────────────────────────────────────────────

    #[test]
    fn test_twsma_zero_cycles() {
        let sma = TwoWaySMA::new(0, 0.04);
        assert!(sma.trained_strain() < 1e-10, "no training → zero strain");
    }

    #[test]
    fn test_twsma_increases_with_cycles() {
        let sma1 = TwoWaySMA::new(10, 0.04);
        let sma2 = TwoWaySMA::new(50, 0.04);
        assert!(
            sma2.trained_strain() > sma1.trained_strain(),
            "more training cycles should give more strain"
        );
    }

    #[test]
    fn test_twsma_saturates() {
        let sma = TwoWaySMA::new(1000, 0.04);
        // After many cycles, should be very close to max
        assert!(
            sma.saturation_fraction() > 0.99,
            "should saturate near 1.0 after many cycles"
        );
    }

    #[test]
    fn test_twsma_train_method() {
        let mut sma = TwoWaySMA::new(0, 0.04);
        sma.train(20);
        assert_eq!(sma.training_cycles, 20);
        let eps = sma.trained_strain();
        assert!(eps > 0.0);
    }

    // ── Nitinol phase diagram ────────────────────────────────────────────────

    #[test]
    fn test_nitinol_phase_diagram_temps() {
        let pd = nitinol_phase_diagram();
        // Check order: Mf < Ms < As < Af
        let mf = pd[1].1;
        let ms = pd[0].1;
        let as_ = pd[2].1;
        let af = pd[3].1;
        assert!(mf < ms, "Mf should be below Ms");
        assert!(ms < as_, "Ms should be below As");
        assert!(as_ < af, "As should be below Af");
    }

    #[test]
    fn test_nitinol_phase_diagram_length() {
        assert_eq!(nitinol_phase_diagram().len(), 4);
    }

    // ── Helper functions ─────────────────────────────────────────────────────

    #[test]
    fn test_cosine_interpolate_bounds() {
        assert!((cosine_interpolate(0.0) - 0.0).abs() < 1e-10);
        assert!((cosine_interpolate(1.0) - 1.0).abs() < 1e-10);
        assert!((cosine_interpolate(0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_conductivity_pure_phases() {
        let k = thermal_conductivity(0.0, 10.0, 20.0);
        assert!((k - 20.0).abs() < 1e-10, "pure austenite should give k_A");
        let k = thermal_conductivity(1.0, 10.0, 20.0);
        assert!((k - 10.0).abs() < 1e-10, "pure martensite should give k_M");
    }

    #[test]
    fn test_thermal_conductivity_mixed() {
        let k = thermal_conductivity(0.5, 10.0, 20.0);
        assert!((k - 15.0).abs() < 1e-10, "mixed should give average");
    }

    #[test]
    fn test_thermal_conductivity_clamps() {
        // xi > 1 should clamp to 1
        let k = thermal_conductivity(2.0, 10.0, 20.0);
        assert!((k - 10.0).abs() < 1e-10);
        // xi < 0 should clamp to 0
        let k2 = thermal_conductivity(-1.0, 10.0, 20.0);
        assert!((k2 - 20.0).abs() < 1e-10);
    }
}
