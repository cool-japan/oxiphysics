// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface adsorption and desorption models for molecular dynamics.
//!
//! This module provides:
//! - [`LangmuirIsotherm`]: Single-layer adsorption equilibrium.
//! - [`FreundlichIsotherm`]: Empirical multi-site adsorption isotherm.
//! - [`BETIsotherm`]: Brunauer–Emmett–Teller multi-layer isotherm.
//! - [`AdsorptionKinetics`]: Rate equations for adsorption/desorption dynamics.
//! - [`SurfaceEnergetics`]: Arrhenius-based surface reaction rates.
//! - [`StickingCoefficient`]: Coverage- and temperature-dependent sticking probability.
//! - [`henry_constant`]: Temperature-corrected Henry's law constant.
//! - [`dubinin_radushkevich`]: Dubinin–Radushkevich micropore adsorption model.
//!
//! # References
//! - Langmuir, I. (1918). *J. Am. Chem. Soc.* 40, 1361.
//! - Freundlich, H. (1907). *Z. Phys. Chem.* 57, 385.
//! - Brunauer, S., Emmett, P. H., & Teller, E. (1938). *J. Am. Chem. Soc.* 60, 309.
//! - Dubinin, M. M., & Radushkevich, L. V. (1947). *Dokl. Akad. Nauk SSSR* 55, 331.

#![allow(dead_code)]

/// Universal gas constant R (J mol⁻¹ K⁻¹).
pub const R_GAS: f64 = 8.314_462_618;
/// Boltzmann constant k_B (J K⁻¹).
pub const K_B: f64 = 1.380_649e-23;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Temperature-corrected Henry's law constant.
///
/// H(T) = k_henry · exp(−k_henry / (R · T))
///
/// where `k_henry` is treated both as the base constant and as an activation
/// energy scale (J mol⁻¹) in this simplified single-parameter form.
///
/// ```no_run
/// use oxiphysics_md::adsorption::henry_constant;
/// let h = henry_constant(1.0, 300.0);
/// assert!(h > 0.0);
/// ```
pub fn henry_constant(k_henry: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    k_henry * (-k_henry / (R_GAS * temperature)).exp()
}

/// Dubinin–Radushkevich (DR) adsorption model: amount adsorbed (mol g⁻¹).
///
/// q = q_m · exp(−K · ε²)
///
/// where ε is the Polanyi adsorption potential (J mol⁻¹) and K (mol² J⁻²) is
/// the DR constant related to the mean free energy of adsorption.
///
/// ```no_run
/// use oxiphysics_md::adsorption::dubinin_radushkevich;
/// let q = dubinin_radushkevich(1.0, 1e-6, 0.0);
/// assert!((q - 1.0).abs() < 1e-12);
/// ```
pub fn dubinin_radushkevich(q_m: f64, k_dr: f64, epsilon: f64) -> f64 {
    q_m * (-k_dr * epsilon * epsilon).exp()
}

// ---------------------------------------------------------------------------
// LangmuirIsotherm
// ---------------------------------------------------------------------------

/// Langmuir adsorption isotherm.
///
/// Models monolayer adsorption on a homogeneous surface:
///
/// ```text
/// θ = K_eq · c / (1 + K_eq · c)
/// q = n_max · θ
/// ```
///
/// where `K_eq` is the equilibrium (Langmuir) constant (L mol⁻¹ or m³ mol⁻¹),
/// `c` is the adsorbate concentration in solution, and `n_max` is the maximum
/// surface coverage (mol m⁻²).
#[derive(Debug, Clone)]
pub struct LangmuirIsotherm {
    /// Equilibrium adsorption constant K_eq (L mol⁻¹).
    pub k_eq: f64,
    /// Maximum surface coverage n_max (mol m⁻²).
    pub n_max: f64,
}

impl LangmuirIsotherm {
    /// Create a new [`LangmuirIsotherm`].
    pub fn new(k_eq: f64, n_max: f64) -> Self {
        Self { k_eq, n_max }
    }

    /// Fractional surface coverage θ ∈ \[0, 1) at a given concentration c.
    ///
    /// ```
    /// use oxiphysics_md::adsorption::LangmuirIsotherm;
    /// let iso = LangmuirIsotherm::new(1.0, 1.0);
    /// let theta = iso.coverage(1.0);
    /// assert!((theta - 0.5).abs() < 1e-12);
    /// ```
    pub fn coverage(&self, concentration: f64) -> f64 {
        if concentration <= 0.0 || self.k_eq <= 0.0 {
            return 0.0;
        }
        let kc = self.k_eq * concentration;
        kc / (1.0 + kc)
    }

    /// Amount adsorbed q = n_max · θ (mol m⁻²).
    pub fn adsorbed_amount(&self, concentration: f64) -> f64 {
        self.n_max * self.coverage(concentration)
    }

    /// Concentration c that gives a fractional coverage θ.
    ///
    /// c = θ / (K_eq · (1 − θ))
    ///
    /// Returns `f64::INFINITY` if θ ≥ 1 and 0 if θ ≤ 0.
    ///
    /// ```
    /// use oxiphysics_md::adsorption::LangmuirIsotherm;
    /// let iso = LangmuirIsotherm::new(2.0, 1.0);
    /// let c = iso.concentration_at_coverage(0.5);
    /// assert!((c - 0.5).abs() < 1e-12);
    /// ```
    pub fn concentration_at_coverage(&self, theta: f64) -> f64 {
        if theta <= 0.0 {
            return 0.0;
        }
        if theta >= 1.0 {
            return f64::INFINITY;
        }
        theta / (self.k_eq * (1.0 - theta))
    }
}

// ---------------------------------------------------------------------------
// FreundlichIsotherm
// ---------------------------------------------------------------------------

/// Freundlich adsorption isotherm (empirical, heterogeneous surfaces).
///
/// ```text
/// q = K_f · c^(1/n)
/// ```
///
/// where `K_f` is the Freundlich capacity factor (mol^(1−1/n) L^(1/n) m⁻²) and
/// `n` is the Freundlich intensity parameter (dimensionless, n > 1 for favourable).
#[derive(Debug, Clone)]
pub struct FreundlichIsotherm {
    /// Freundlich capacity factor K_f.
    pub k_f: f64,
    /// Freundlich exponent n (dimensionless).
    pub n: f64,
}

impl FreundlichIsotherm {
    /// Create a new [`FreundlichIsotherm`].
    pub fn new(k_f: f64, n: f64) -> Self {
        Self { k_f, n }
    }

    /// Surface coverage q = K_f · c^(1/n).
    ///
    /// ```no_run
    /// use oxiphysics_md::adsorption::FreundlichIsotherm;
    /// let iso = FreundlichIsotherm::new(1.0, 1.0);
    /// // n = 1 → linear: q = K_f * c
    /// assert!((iso.coverage(2.0) - 2.0).abs() < 1e-12);
    /// ```
    pub fn coverage(&self, concentration: f64) -> f64 {
        if concentration <= 0.0 || self.k_f <= 0.0 || self.n <= 0.0 {
            return 0.0;
        }
        self.k_f * concentration.powf(1.0 / self.n)
    }
}

// ---------------------------------------------------------------------------
// BETIsotherm
// ---------------------------------------------------------------------------

/// Brunauer–Emmett–Teller (BET) multi-layer adsorption isotherm.
///
/// ```text
/// V / V_m = C · (p/p0) / [(1 − p/p0)(1 − p/p0 + C · p/p0)]
/// ```
///
/// where `C` is the BET constant (dimensionless) and `V_m` is the monolayer
/// volume (same units as V).
#[derive(Debug, Clone)]
pub struct BETIsotherm {
    /// BET adsorption constant C (dimensionless).
    pub c: f64,
    /// Monolayer adsorbed volume V_m (m³ kg⁻¹ or equivalent).
    pub v_m: f64,
}

impl BETIsotherm {
    /// Create a new [`BETIsotherm`].
    pub fn new(c: f64, v_m: f64) -> Self {
        Self { c, v_m }
    }

    /// Adsorbed volume V at relative pressure p/p0.
    ///
    /// Returns 0 for unphysical inputs (p/p0 outside (0, 1)).
    ///
    /// ```
    /// use oxiphysics_md::adsorption::BETIsotherm;
    /// let bet = BETIsotherm::new(100.0, 1.0);
    /// let v = bet.adsorbed_volume(0.1, 1.0);
    /// assert!(v > 0.0);
    /// ```
    pub fn adsorbed_volume(&self, p: f64, p0: f64) -> f64 {
        if p0 <= 0.0 || p <= 0.0 {
            return 0.0;
        }
        let x = p / p0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        let denom = (1.0 - x) * (1.0 - x + self.c * x);
        if denom <= 0.0 {
            return 0.0;
        }
        self.v_m * self.c * x / denom
    }

    /// BET surface area can be extracted from the linear BET plot.
    /// Returns the slope of the linearized form: (p/p0) / \[V(1 − p/p0)\].
    pub fn bet_linear(&self, p: f64, p0: f64) -> f64 {
        if p0 <= 0.0 || p <= 0.0 {
            return 0.0;
        }
        let x = p / p0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        let v = self.adsorbed_volume(p, p0);
        if v <= 0.0 {
            return 0.0;
        }
        x / (v * (1.0 - x))
    }
}

// ---------------------------------------------------------------------------
// AdsorptionKinetics
// ---------------------------------------------------------------------------

/// First-order adsorption/desorption kinetics (Langmuir kinetics).
///
/// Rate equation:
///
/// ```text
/// dθ/dt = k_ads · c · (1 − θ) − k_des · θ
/// ```
///
/// Equilibrium coverage: θ_eq = k_ads · c / (k_ads · c + k_des)
#[derive(Debug, Clone)]
pub struct AdsorptionKinetics {
    /// Adsorption rate constant k_ads (L mol⁻¹ s⁻¹).
    pub k_ads: f64,
    /// Desorption rate constant k_des (s⁻¹).
    pub k_des: f64,
}

impl AdsorptionKinetics {
    /// Create a new [`AdsorptionKinetics`] model.
    pub fn new(k_ads: f64, k_des: f64) -> Self {
        Self { k_ads, k_des }
    }

    /// Net adsorption rate dθ/dt (s⁻¹).
    ///
    /// ```no_run
    /// use oxiphysics_md::adsorption::AdsorptionKinetics;
    /// let kin = AdsorptionKinetics::new(1.0, 0.0);
    /// let rate = kin.compute_rate(0.0, 1.0);
    /// assert!((rate - 1.0).abs() < 1e-12);
    /// ```
    pub fn compute_rate(&self, theta: f64, concentration: f64) -> f64 {
        let theta_clamped = theta.clamp(0.0, 1.0);
        self.k_ads * concentration * (1.0 - theta_clamped) - self.k_des * theta_clamped
    }

    /// Equilibrium coverage θ_eq for a given concentration.
    ///
    /// θ_eq = K_eq · c / (1 + K_eq · c),  K_eq = k_ads / k_des
    ///
    /// Returns 1 if k_des == 0 and concentration > 0, else 0.
    ///
    /// ```no_run
    /// use oxiphysics_md::adsorption::AdsorptionKinetics;
    /// let kin = AdsorptionKinetics::new(2.0, 2.0);
    /// let teq = kin.equilibrium_coverage(1.0);
    /// assert!((teq - 0.5).abs() < 1e-12);
    /// ```
    pub fn equilibrium_coverage(&self, concentration: f64) -> f64 {
        if concentration <= 0.0 {
            return 0.0;
        }
        if self.k_des <= 0.0 {
            return if concentration > 0.0 { 1.0 } else { 0.0 };
        }
        let k_eq = self.k_ads / self.k_des;
        let kc = k_eq * concentration;
        kc / (1.0 + kc)
    }

    /// Integrate the rate equation by a simple Euler step of size `dt` (s).
    pub fn step(&self, theta: f64, concentration: f64, dt: f64) -> f64 {
        (theta + self.compute_rate(theta, concentration) * dt).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// SurfaceEnergetics
// ---------------------------------------------------------------------------

/// Surface energetics and Arrhenius desorption rate.
#[derive(Debug, Clone)]
pub struct SurfaceEnergetics {
    /// Adsorption energy (J mol⁻¹, negative for exothermic binding).
    pub adsorption_energy: f64,
    /// Activation barrier for desorption (J mol⁻¹, positive).
    pub activation_barrier: f64,
    /// Pre-exponential (attempt) frequency ν₀ (s⁻¹).
    pub attempt_frequency: f64,
}

impl SurfaceEnergetics {
    /// Create a new [`SurfaceEnergetics`] model.
    pub fn new(adsorption_energy: f64, activation_barrier: f64, attempt_frequency: f64) -> Self {
        Self {
            adsorption_energy,
            activation_barrier,
            attempt_frequency,
        }
    }

    /// Arrhenius desorption rate k_des = ν₀ · exp(−E_a / (R·T)) (s⁻¹).
    ///
    /// ```no_run
    /// use oxiphysics_md::adsorption::SurfaceEnergetics;
    /// let se = SurfaceEnergetics::new(-5e4, 5e4, 1e13);
    /// let k = se.arrhenius_rate(300.0);
    /// assert!(k > 0.0);
    /// ```
    pub fn arrhenius_rate(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 || self.activation_barrier < 0.0 {
            return 0.0;
        }
        self.attempt_frequency * (-self.activation_barrier / (R_GAS * temperature)).exp()
    }

    /// Surface residence time τ = 1 / k_des (s).
    pub fn residence_time(&self, temperature: f64) -> f64 {
        let k = self.arrhenius_rate(temperature);
        if k <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / k
    }
}

// ---------------------------------------------------------------------------
// StickingCoefficient
// ---------------------------------------------------------------------------

/// Coverage- and temperature-dependent sticking coefficient.
///
/// ```text
/// S(θ, T) = s0 · (1 − θ)^n_cov · exp(−E_acc / (k_B T))
/// ```
///
/// where `s0` is the bare sticking probability, `n_cov` is the coverage
/// dependence exponent, and `E_acc` is the energy accommodation barrier (J).
#[derive(Debug, Clone)]
pub struct StickingCoefficient {
    /// Bare sticking probability S₀ ∈ \[0, 1\].
    pub s0: f64,
    /// Energy accommodation barrier E_acc (J).
    pub energy_accommodation: f64,
    /// Coverage dependence exponent (typically 1 for Langmuir, 2 for precursor).
    pub coverage_dependence: f64,
}

impl StickingCoefficient {
    /// Create a new [`StickingCoefficient`].
    pub fn new(s0: f64, energy_accommodation: f64, coverage_dependence: f64) -> Self {
        Self {
            s0,
            energy_accommodation,
            coverage_dependence,
        }
    }

    /// Sticking coefficient S(θ, T) ∈ \[0, 1\].
    ///
    /// ```
    /// use oxiphysics_md::adsorption::StickingCoefficient;
    /// let sc = StickingCoefficient::new(1.0, 0.0, 1.0);
    /// // At zero coverage and zero barrier: S = 1
    /// assert!((sc.compute(0.0, 300.0) - 1.0).abs() < 1e-12);
    /// ```
    pub fn compute(&self, theta: f64, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let theta_c = theta.clamp(0.0, 1.0);
        let coverage_factor = (1.0 - theta_c).powf(self.coverage_dependence);
        let thermal_factor = (-self.energy_accommodation / (K_B * temperature)).exp();
        (self.s0 * coverage_factor * thermal_factor).clamp(0.0, 1.0)
    }

    /// Sticking coefficient at full coverage (θ = 1) — should be 0.
    pub fn saturated(&self, temperature: f64) -> f64 {
        self.compute(1.0, temperature)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- henry_constant --------------------------------------------------

    #[test]
    fn test_henry_constant_positive() {
        let h = henry_constant(1.0, 300.0);
        assert!(h > 0.0, "h={h}");
    }

    #[test]
    fn test_henry_constant_zero_temperature() {
        assert_eq!(henry_constant(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_henry_constant_finite() {
        let h = henry_constant(1e4, 300.0);
        assert!(h.is_finite(), "h={h}");
    }

    // ---- dubinin_radushkevich -------------------------------------------

    #[test]
    fn test_dr_zero_epsilon() {
        let q = dubinin_radushkevich(1.0, 1e-6, 0.0);
        assert!((q - 1.0).abs() < 1e-12, "q={q}");
    }

    #[test]
    fn test_dr_positive() {
        let q = dubinin_radushkevich(2.0, 1e-6, 100.0);
        assert!(q > 0.0 && q <= 2.0, "q={q}");
    }

    #[test]
    fn test_dr_decays_with_epsilon() {
        let q1 = dubinin_radushkevich(1.0, 1e-6, 100.0);
        let q2 = dubinin_radushkevich(1.0, 1e-6, 200.0);
        assert!(q2 < q1, "q1={q1}, q2={q2}");
    }

    // ---- LangmuirIsotherm -----------------------------------------------

    #[test]
    fn test_langmuir_half_coverage() {
        // K=1, c=1 → θ = 0.5
        let iso = LangmuirIsotherm::new(1.0, 1.0);
        assert!((iso.coverage(1.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_langmuir_zero_concentration() {
        let iso = LangmuirIsotherm::new(1.0, 1.0);
        assert_eq!(iso.coverage(0.0), 0.0);
    }

    #[test]
    fn test_langmuir_coverage_approaches_one() {
        let iso = LangmuirIsotherm::new(1.0, 1.0);
        let theta = iso.coverage(1e10);
        assert!(theta > 0.999, "theta={theta}");
    }

    #[test]
    fn test_langmuir_adsorbed_amount() {
        let iso = LangmuirIsotherm::new(1.0, 2.0);
        let q = iso.adsorbed_amount(1.0);
        assert!((q - 1.0).abs() < 1e-12, "q={q}");
    }

    #[test]
    fn test_langmuir_concentration_at_half() {
        let iso = LangmuirIsotherm::new(2.0, 1.0);
        let c = iso.concentration_at_coverage(0.5);
        // θ = K_eq*c / (1 + K_eq*c) = 0.5 → K_eq*c = 1 → c = 0.5
        assert!((c - 0.5).abs() < 1e-12, "c={c}");
    }

    #[test]
    fn test_langmuir_concentration_zero_theta() {
        let iso = LangmuirIsotherm::new(1.0, 1.0);
        assert_eq!(iso.concentration_at_coverage(0.0), 0.0);
    }

    #[test]
    fn test_langmuir_concentration_full_coverage() {
        let iso = LangmuirIsotherm::new(1.0, 1.0);
        assert!(iso.concentration_at_coverage(1.0).is_infinite());
    }

    #[test]
    fn test_langmuir_roundtrip() {
        let iso = LangmuirIsotherm::new(3.0, 1.0);
        let c0 = 0.25;
        let theta = iso.coverage(c0);
        let c1 = iso.concentration_at_coverage(theta);
        assert!((c1 - c0).abs() < 1e-10, "c0={c0}, c1={c1}");
    }

    // ---- FreundlichIsotherm ---------------------------------------------

    #[test]
    fn test_freundlich_linear_n1() {
        let iso = FreundlichIsotherm::new(1.0, 1.0);
        assert!((iso.coverage(2.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_freundlich_zero_concentration() {
        let iso = FreundlichIsotherm::new(1.0, 2.0);
        assert_eq!(iso.coverage(0.0), 0.0);
    }

    #[test]
    fn test_freundlich_positive() {
        let iso = FreundlichIsotherm::new(1.5, 2.0);
        assert!(iso.coverage(1.0) > 0.0);
    }

    #[test]
    fn test_freundlich_increases_with_concentration() {
        let iso = FreundlichIsotherm::new(1.0, 2.0);
        assert!(iso.coverage(2.0) > iso.coverage(1.0));
    }

    #[test]
    fn test_freundlich_n2() {
        // n=2 → q = K_f * c^0.5
        let iso = FreundlichIsotherm::new(1.0, 2.0);
        let q = iso.coverage(4.0);
        assert!((q - 2.0).abs() < 1e-12, "q={q}");
    }

    // ---- BETIsotherm ----------------------------------------------------

    #[test]
    fn test_bet_positive_volume() {
        let bet = BETIsotherm::new(100.0, 1.0);
        let v = bet.adsorbed_volume(0.1, 1.0);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_bet_zero_p() {
        let bet = BETIsotherm::new(100.0, 1.0);
        assert_eq!(bet.adsorbed_volume(0.0, 1.0), 0.0);
    }

    #[test]
    fn test_bet_zero_p0() {
        let bet = BETIsotherm::new(100.0, 1.0);
        assert_eq!(bet.adsorbed_volume(0.1, 0.0), 0.0);
    }

    #[test]
    fn test_bet_increases_with_pressure() {
        let bet = BETIsotherm::new(50.0, 1.0);
        let v1 = bet.adsorbed_volume(0.1, 1.0);
        let v2 = bet.adsorbed_volume(0.5, 1.0);
        assert!(v2 > v1, "v1={v1}, v2={v2}");
    }

    #[test]
    fn test_bet_diverges_at_saturation() {
        let bet = BETIsotherm::new(100.0, 1.0);
        let v = bet.adsorbed_volume(1.0, 1.0);
        assert!(v.is_infinite(), "v={v}");
    }

    #[test]
    fn test_bet_linear_positive() {
        let bet = BETIsotherm::new(100.0, 1.0);
        let val = bet.bet_linear(0.2, 1.0);
        assert!(val > 0.0, "val={val}");
    }

    // ---- AdsorptionKinetics ---------------------------------------------

    #[test]
    fn test_kinetics_rate_zero_theta() {
        let kin = AdsorptionKinetics::new(1.0, 0.0);
        let rate = kin.compute_rate(0.0, 1.0);
        assert!((rate - 1.0).abs() < 1e-12, "rate={rate}");
    }

    #[test]
    fn test_kinetics_rate_full_coverage() {
        let kin = AdsorptionKinetics::new(1.0, 1.0);
        let rate = kin.compute_rate(1.0, 1.0);
        // k_ads * c * 0 - k_des * 1 = -1
        assert!((rate + 1.0).abs() < 1e-12, "rate={rate}");
    }

    #[test]
    fn test_kinetics_equilibrium() {
        let kin = AdsorptionKinetics::new(2.0, 2.0);
        let theta_eq = kin.equilibrium_coverage(1.0);
        assert!((theta_eq - 0.5).abs() < 1e-12, "theta_eq={theta_eq}");
    }

    #[test]
    fn test_kinetics_equilibrium_zero_concentration() {
        let kin = AdsorptionKinetics::new(2.0, 2.0);
        assert_eq!(kin.equilibrium_coverage(0.0), 0.0);
    }

    #[test]
    fn test_kinetics_equilibrium_no_desorption() {
        let kin = AdsorptionKinetics::new(1.0, 0.0);
        assert!((kin.equilibrium_coverage(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_kinetics_step_convergence() {
        // dθ/dt = k_ads*c*(1-θ) - k_des*θ  with k_ads=k_des=1, c=1
        // τ = 1/(k_ads*c + k_des) = 0.5 s → need ~10 τ = 5 s → 5000 steps at dt=1e-3
        let kin = AdsorptionKinetics::new(1.0, 1.0);
        let mut theta = 0.0;
        for _ in 0..5000 {
            theta = kin.step(theta, 1.0, 1e-3);
        }
        let theta_eq = kin.equilibrium_coverage(1.0);
        assert!(
            (theta - theta_eq).abs() < 0.01,
            "theta={theta}, eq={theta_eq}"
        );
    }

    #[test]
    fn test_kinetics_step_clamps() {
        let kin = AdsorptionKinetics::new(1e6, 0.0);
        let theta = kin.step(0.0, 1.0, 1e3);
        assert!(theta <= 1.0, "theta={theta}");
    }

    // ---- SurfaceEnergetics ----------------------------------------------

    #[test]
    fn test_arrhenius_rate_positive() {
        let se = SurfaceEnergetics::new(-5e4, 5e4, 1e13);
        let k = se.arrhenius_rate(300.0);
        assert!(k > 0.0, "k={k}");
    }

    #[test]
    fn test_arrhenius_rate_zero_temperature() {
        let se = SurfaceEnergetics::new(-5e4, 5e4, 1e13);
        assert_eq!(se.arrhenius_rate(0.0), 0.0);
    }

    #[test]
    fn test_arrhenius_rate_increases_with_temperature() {
        let se = SurfaceEnergetics::new(-5e4, 5e4, 1e13);
        let k1 = se.arrhenius_rate(300.0);
        let k2 = se.arrhenius_rate(600.0);
        assert!(k2 > k1, "k1={k1}, k2={k2}");
    }

    #[test]
    fn test_residence_time_positive() {
        let se = SurfaceEnergetics::new(-5e4, 5e4, 1e13);
        let tau = se.residence_time(300.0);
        assert!(tau > 0.0 && tau.is_finite(), "tau={tau}");
    }

    // ---- StickingCoefficient --------------------------------------------

    #[test]
    fn test_sticking_no_barrier_zero_coverage() {
        let sc = StickingCoefficient::new(1.0, 0.0, 1.0);
        assert!((sc.compute(0.0, 300.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_sticking_full_coverage_zero() {
        let sc = StickingCoefficient::new(1.0, 0.0, 1.0);
        assert!((sc.saturated(300.0)).abs() < 1e-12);
    }

    #[test]
    fn test_sticking_decreases_with_coverage() {
        let sc = StickingCoefficient::new(1.0, 0.0, 1.0);
        let s1 = sc.compute(0.2, 300.0);
        let s2 = sc.compute(0.8, 300.0);
        assert!(s2 < s1, "s1={s1}, s2={s2}");
    }

    #[test]
    fn test_sticking_zero_temperature() {
        let sc = StickingCoefficient::new(1.0, 1e-21, 1.0);
        assert_eq!(sc.compute(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_sticking_thermal_activation() {
        // With a finite barrier, higher T → higher S
        let sc = StickingCoefficient::new(1.0, 1e-20, 1.0);
        let s1 = sc.compute(0.0, 300.0);
        let s2 = sc.compute(0.0, 1000.0);
        assert!(s2 > s1, "s1={s1}, s2={s2}");
    }

    #[test]
    fn test_sticking_clipped_to_one() {
        let sc = StickingCoefficient::new(2.0, 0.0, 1.0);
        let s = sc.compute(0.0, 300.0);
        assert!(s <= 1.0, "s={s}");
    }

    // ---- integration / cross-model checks --------------------------------

    #[test]
    fn test_langmuir_kinetics_consistency() {
        // Equilibrium coverage from kinetics == Langmuir isotherm coverage
        let k_ads = 3.0;
        let k_des = 1.5;
        let c = 2.0;
        let kin = AdsorptionKinetics::new(k_ads, k_des);
        let iso = LangmuirIsotherm::new(k_ads / k_des, 1.0);
        let theta_kin = kin.equilibrium_coverage(c);
        let theta_iso = iso.coverage(c);
        assert!(
            (theta_kin - theta_iso).abs() < 1e-12,
            "kin={theta_kin}, iso={theta_iso}"
        );
    }

    #[test]
    fn test_bet_monolayer_volume_limit() {
        // For very small p/p0 BET should approach Langmuir-like behaviour
        let bet = BETIsotherm::new(1e6, 1.0);
        let v = bet.adsorbed_volume(1e-6, 1.0);
        // At very low coverage, V ≈ V_m * C * x  (≈ V_m for C >> 1, x << 1)
        assert!(v > 0.0 && v <= bet.v_m * 2.0, "v={v}");
    }
}
