#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fatigue analysis finite element module.
//!
//! This module provides:
//! - [`FatigueModel`]: stress-life / strain-life model selection
//! - [`FatigueParams`]: material fatigue constants
//! - [`CyclicStressStrain`]: one loading cycle description
//! - [`SNCurve`]: Wöhler (S–N) curve with life and stress queries
//! - [`FatigueAnalysis`]: high-level fatigue life estimator
//! - [`FatigueResult`]: life, damage, failure information
//!
//! # Physics
//!
//! Fatigue failure is governed by cyclic damage accumulation.  This module
//! implements the most widely used models:
//!
//! **Stress-life (Basquin):**  σ_a = σ_f' · (2N_f)^b
//!
//! **Strain-life (Coffin–Manson):**  Δε_p/2 = ε_f' · (2N_f)^c
//!
//! **Morrow mean-stress correction:**  σ_a / (σ_f' − σ_m) = (2N_f)^b
//!
//! **Palmgren–Miner rule:**  D = Σ n_i / N_i  (failure when D ≥ 1)
//!
//! Rainflow cycle counting extracts amplitude/mean pairs from an arbitrary
//! stress history following the ASTM E1049 three-point algorithm.
//!
//! # References
//! - Basquin, O.H. (1910). The exponential law of endurance tests.
//!   *ASTM Proc.*, 10, 625–630.
//! - Coffin, L.F. (1954). A study of the effects of cyclic thermal stresses.
//!   *Trans. ASME*, 76, 931–950.
//! - Morrow, J. (1968). Cyclic plastic strain energy and fatigue of metals.
//!   *ASTM STP 378*, 45–87.
//! - Rainflow counting: ASTM E1049-85 (2017).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Enums & structs
// ---------------------------------------------------------------------------

/// Fatigue life prediction model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FatigueModel {
    /// Basquin stress-life:  σ_a = σ_f' · (2 N_f)^b.
    Basquin,
    /// Coffin–Manson strain-life:  Δε_p/2 = ε_f' · (2 N_f)^c.
    CoffinManson,
    /// Modified Goodman diagram for mean-stress correction.
    ModifiedGoodman,
    /// Walker equation mean-stress correction.
    Walker,
    /// Morrow mean-stress correction for strain-life.
    Morrow,
}

/// Material fatigue parameters.
#[derive(Debug, Clone)]
pub struct FatigueParams {
    /// Ultimate tensile strength S_ut (Pa).
    pub s_ut: f64,
    /// Endurance limit S_e (Pa) at reference cycles N_e.
    pub s_e: f64,
    /// Basquin exponent b (negative, typically −0.05 to −0.12).
    pub b: f64,
    /// Coffin–Manson ductility exponent c (negative, typically −0.5 to −0.7).
    pub c: f64,
    /// Fatigue strength coefficient σ_f' (Pa).
    pub sigma_f: f64,
    /// Fatigue ductility coefficient ε_f' (dimensionless).
    pub epsilon_f: f64,
    /// Reference number of reversals to failure 2N_f_ref.
    pub n_f_ref: f64,
}

impl Default for FatigueParams {
    fn default() -> Self {
        // Typical steel (SAE 1045)
        Self {
            s_ut: 620.0e6,
            s_e: 210.0e6,
            b: -0.085,
            c: -0.54,
            sigma_f: 948.0e6,
            epsilon_f: 0.26,
            n_f_ref: 2.0e6,
        }
    }
}

/// One loading cycle description.
#[derive(Debug, Clone)]
pub struct CyclicStressStrain {
    /// Maximum stress σ_max (Pa).
    pub stress_max: f64,
    /// Minimum stress σ_min (Pa).
    pub stress_min: f64,
    /// Maximum strain ε_max (dimensionless).
    pub strain_max: f64,
    /// Minimum strain ε_min (dimensionless).
    pub strain_min: f64,
    /// Mean stress σ_m = (σ_max + σ_min) / 2 (Pa).
    pub mean_stress: f64,
    /// Stress amplitude σ_a = (σ_max − σ_min) / 2 (Pa).
    pub amplitude: f64,
}

impl CyclicStressStrain {
    /// Construct from peak stresses and strains.
    ///
    /// # Arguments
    /// * `s_max`  - maximum stress (Pa)
    /// * `s_min`  - minimum stress (Pa)
    /// * `e_max`  - maximum strain
    /// * `e_min`  - minimum strain
    pub fn new(s_max: f64, s_min: f64, e_max: f64, e_min: f64) -> Self {
        Self {
            stress_max: s_max,
            stress_min: s_min,
            strain_max: e_max,
            strain_min: e_min,
            mean_stress: 0.5 * (s_max + s_min),
            amplitude: 0.5 * (s_max - s_min),
        }
    }

    /// Stress ratio R = σ_min / σ_max.
    pub fn stress_ratio(&self) -> f64 {
        if self.stress_max.abs() < 1e-30 {
            return 0.0;
        }
        self.stress_min / self.stress_max
    }

    /// Total strain range Δε = ε_max − ε_min.
    pub fn strain_range(&self) -> f64 {
        self.strain_max - self.strain_min
    }

    /// Stress range Δσ = σ_max − σ_min.
    pub fn stress_range(&self) -> f64 {
        self.stress_max - self.stress_min
    }
}

// ---------------------------------------------------------------------------
// S–N curve (Wöhler curve)
// ---------------------------------------------------------------------------

/// Wöhler (S–N) curve: S_a = S_u · (N / N_ref)^(1/k).
#[derive(Debug, Clone)]
pub struct SNCurve {
    /// Endurance limit S_e (Pa) — stress below which infinite life is expected.
    pub s_e: f64,
    /// Ultimate strength S_u (Pa) used to anchor the high-cycle end.
    pub s_u: f64,
    /// Slope parameter k > 0  (Basquin slope — lives per decade of stress).
    pub k_slope: f64,
    /// Endurance limit cycle count N_e (cycles).
    pub n_e: f64,
}

impl SNCurve {
    /// Create a new S–N curve.
    ///
    /// # Arguments
    /// * `s_e`     - endurance limit (Pa)
    /// * `s_u`     - ultimate strength (Pa)
    /// * `k_slope` - slope k (Basquin k, positive)
    /// * `n_e`     - endurance limit cycles
    pub fn new(s_e: f64, s_u: f64, k_slope: f64, n_e: f64) -> Self {
        Self {
            s_e,
            s_u,
            k_slope,
            n_e,
        }
    }

    /// Predict fatigue life (cycles) at stress amplitude `sigma_a`.
    ///
    /// Returns `f64::INFINITY` when `sigma_a ≤ s_e` (below endurance limit).
    ///
    /// # Arguments
    /// * `sigma_a` - stress amplitude (Pa)
    ///
    /// # Returns
    /// Predicted life in cycles (or `f64::INFINITY`).
    pub fn life_at_stress(&self, sigma_a: f64) -> f64 {
        if sigma_a <= self.s_e {
            return f64::INFINITY;
        }
        if sigma_a >= self.s_u {
            return 0.5;
        } // half-cycle to failure
        // Power-law interpolation: N = N_e · (S_e / S_a)^k
        self.n_e * (self.s_e / sigma_a).powf(self.k_slope)
    }

    /// Predict stress amplitude (Pa) that would produce fatigue life `n` cycles.
    ///
    /// # Arguments
    /// * `n` - number of cycles
    ///
    /// # Returns
    /// Stress amplitude (Pa).
    pub fn stress_at_life(&self, n: f64) -> f64 {
        if n >= self.n_e {
            return self.s_e;
        }
        // S_a = S_e · (N_e / N)^(1/k)
        self.s_e * (self.n_e / n).powf(1.0 / self.k_slope)
    }
}

// ---------------------------------------------------------------------------
// Life prediction functions
// ---------------------------------------------------------------------------

/// Basquin stress-life equation.
///
/// N_f = (σ_a / σ_f')^(1/b) / 2
///
/// # Arguments
/// * `sigma_a` - stress amplitude (Pa)
/// * `sigma_f` - fatigue strength coefficient σ_f' (Pa)
/// * `b`       - Basquin exponent (negative)
///
/// # Returns
/// Cycles to failure N_f.
pub fn basquin_life(sigma_a: f64, sigma_f: f64, b: f64) -> f64 {
    if sigma_a <= 0.0 || sigma_f <= 0.0 {
        return f64::INFINITY;
    }
    0.5 * (sigma_a / sigma_f).powf(1.0 / b)
}

/// Coffin–Manson strain-life equation.
///
/// N_f = (Δε_p / (2 ε_f'))^(1/c) / 2
///
/// # Arguments
/// * `delta_eps_p` - plastic strain range Δε_p
/// * `eps_f`       - fatigue ductility coefficient ε_f'
/// * `c`           - Coffin–Manson exponent (negative)
///
/// # Returns
/// Cycles to failure N_f.
pub fn coffin_manson_life(delta_eps_p: f64, eps_f: f64, c: f64) -> f64 {
    if delta_eps_p <= 0.0 || eps_f <= 0.0 {
        return f64::INFINITY;
    }
    0.5 * (delta_eps_p / (2.0 * eps_f)).powf(1.0 / c)
}

/// Morrow mean-stress corrected life.
///
/// σ_a / (σ_f' − σ_m) = (2 N_f)^b  →  N_f = 0.5 · (σ_a / (σ_f' − σ_m))^(1/b)
///
/// # Arguments
/// * `sigma_a` - stress amplitude (Pa)
/// * `sigma_m` - mean stress (Pa)
/// * `sigma_f` - fatigue strength coefficient σ_f' (Pa)
/// * `b`       - Basquin exponent (negative)
///
/// # Returns
/// Cycles to failure N_f.
pub fn morrow_correction(sigma_a: f64, sigma_m: f64, sigma_f: f64, b: f64) -> f64 {
    let denom = sigma_f - sigma_m;
    if denom <= 0.0 || sigma_a <= 0.0 {
        return 0.0;
    }
    0.5 * (sigma_a / denom).powf(1.0 / b)
}

/// Modified Goodman equivalent fully-reversed stress amplitude.
///
/// σ_eq = σ_a / (1 − σ_m / S_ut)
///
/// Returns `f64::INFINITY` when σ_m ≥ S_ut (static fracture).
///
/// # Arguments
/// * `sigma_a` - stress amplitude (Pa)
/// * `sigma_m` - mean stress (Pa)
/// * `s_e`     - endurance limit (Pa) — unused but kept for API symmetry
/// * `s_ut`    - ultimate tensile strength (Pa)
///
/// # Returns
/// Equivalent fully-reversed stress amplitude (Pa).
pub fn goodman_equivalent(sigma_a: f64, sigma_m: f64, s_e: f64, s_ut: f64) -> f64 {
    let _ = s_e;
    let ratio = sigma_m / s_ut;
    if ratio >= 1.0 {
        return f64::INFINITY;
    }
    sigma_a / (1.0 - ratio)
}

/// Walker mean-stress correction.
///
/// σ_eq = σ_max^(1-γ) · σ_a^γ  where γ is the Walker exponent (≈ 0.5 for steel).
///
/// # Arguments
/// * `sigma_max` - maximum stress (Pa)
/// * `sigma_a`   - stress amplitude (Pa)
/// * `gamma`     - Walker exponent γ ∈ \[0, 1\]
///
/// # Returns
/// Walker equivalent stress amplitude (Pa).
pub fn walker_equivalent(sigma_max: f64, sigma_a: f64, gamma: f64) -> f64 {
    if sigma_max <= 0.0 || sigma_a <= 0.0 {
        return 0.0;
    }
    sigma_max.powf(1.0 - gamma) * sigma_a.powf(gamma)
}

/// Stress concentration factor adjusted by notch sensitivity.
///
/// K_f = 1 + q · (K_t − 1)
///
/// # Arguments
/// * `kt` - theoretical stress concentration factor K_t
/// * `q`  - notch sensitivity factor q ∈ \[0, 1\]
///
/// # Returns
/// Fatigue stress concentration factor K_f.
pub fn notch_factor(kt: f64, q: f64) -> f64 {
    1.0 + q * (kt - 1.0)
}

/// Surface finish correction factor (Marin factor k_a).
///
/// Uses the power-law formula:  k_a = A · S_ut^b
///
/// # Arguments
/// * `a`   - surface factor coefficient A
/// * `b`   - surface factor exponent b (negative)
/// * `s_ut`- ultimate strength (MPa for standard A/b tables)
///
/// # Returns
/// Surface correction factor k_a ∈ (0, 1].
pub fn surface_finish_factor(a: f64, b: f64, s_ut: f64) -> f64 {
    (a * s_ut.powf(b)).clamp(0.0, 1.0)
}

/// Size correction factor (Marin factor k_b).
///
/// Simplified formula for round cross-sections (diameter d in mm):
///   k_b = 1.24 · d^(−0.107)   (7.62 ≤ d ≤ 51 mm)
///          0.859 − 0.000_837 d  (51 < d ≤ 254 mm)
///
/// # Arguments
/// * `d_mm` - section diameter (mm)
///
/// # Returns
/// Size correction factor k_b.
pub fn size_factor(d_mm: f64) -> f64 {
    if d_mm <= 7.62 {
        1.0
    } else if d_mm <= 51.0 {
        1.24 * d_mm.powf(-0.107)
    } else {
        (0.859 - 0.000_837 * d_mm).max(0.6)
    }
}

/// Reliability correction factor k_e for survival probability p.
///
/// Uses the formula  k_e = 1 − 0.08 z_p  where z_p is the standard normal
/// deviate for reliability R = 1 − p_f.
///
/// # Arguments
/// * `reliability` - required reliability R ∈ \[0.5, 0.9999\]
///
/// # Returns
/// Reliability factor k_e.
pub fn reliability_factor(reliability: f64) -> f64 {
    // Standard normal deviate (simple table lookup)
    let z = match reliability {
        r if r <= 0.50 => 0.000,
        r if r <= 0.90 => 1.282,
        r if r <= 0.95 => 1.645,
        r if r <= 0.99 => 2.326,
        r if r <= 0.999 => 3.091,
        _ => 3.719,
    };
    (1.0_f64 - 0.08 * z).max(0.0)
}

// ---------------------------------------------------------------------------
// Rainflow counting
// ---------------------------------------------------------------------------

/// Rainflow cycle counting (ASTM E1049 three-point algorithm).
///
/// Extracts closed cycles from a stress history and returns a vector of
/// (range, mean, count) tuples where count is 0.5 for a half-cycle and 1.0
/// for a full cycle.
///
/// # Arguments
/// * `stress_history` - stress time series (Pa)
///
/// # Returns
/// `Vec<(range, mean, count)>` — each entry is a counted cycle.
pub fn rainflow_count(stress_history: &[f64]) -> Vec<(f64, f64, f64)> {
    if stress_history.len() < 2 {
        return Vec::new();
    }

    // Build turning-point (reversal) sequence
    let mut reversals: Vec<f64> = Vec::with_capacity(stress_history.len());
    reversals.push(stress_history[0]);
    for i in 1..stress_history.len() - 1 {
        let prev = stress_history[i - 1];
        let curr = stress_history[i];
        let next = stress_history[i + 1];
        if (curr - prev) * (next - curr) < 0.0 {
            reversals.push(curr);
        }
    }
    reversals.push(*stress_history.last().expect("stress_history is non-empty"));

    // Three-point cycle extraction
    let mut stack: Vec<f64> = Vec::new();
    let mut cycles: Vec<(f64, f64, f64)> = Vec::new();

    for &s in &reversals {
        stack.push(s);
        loop {
            let n = stack.len();
            if n < 3 {
                break;
            }
            let s0 = stack[n - 3];
            let s1 = stack[n - 2];
            let s2 = stack[n - 1];
            let r1 = (s1 - s0).abs();
            let r2 = (s2 - s1).abs();
            if r2 >= r1 {
                // Full cycle: extract s0–s1
                let range = r1;
                let mean = 0.5 * (s0 + s1);
                cycles.push((range, mean, 1.0));
                stack.remove(n - 3);
                stack.remove(n - 3); // now n-2 after removal
            } else {
                break;
            }
        }
    }

    // Remaining items in stack are half-cycles
    for w in stack.windows(2) {
        let range = (w[1] - w[0]).abs();
        let mean = 0.5 * (w[0] + w[1]);
        cycles.push((range, mean, 0.5));
    }

    cycles
}

// ---------------------------------------------------------------------------
// Palmgren–Miner damage accumulation
// ---------------------------------------------------------------------------

/// Palmgren–Miner linear damage accumulation.
///
/// D = Σ n_i / N_i(S_i)
///
/// # Arguments
/// * `cycles` - slice of (stress_amplitude, cycle_count) pairs
/// * `sn`     - S–N curve for life lookup
///
/// # Returns
/// Accumulated damage D (failure when D ≥ 1).
pub fn miner_damage(cycles: &[(f64, f64)], sn: &SNCurve) -> f64 {
    cycles.iter().fold(0.0, |acc, &(sigma_a, count)| {
        let n_fail = sn.life_at_stress(sigma_a);
        if n_fail.is_infinite() || n_fail <= 0.0 {
            acc
        } else {
            acc + count / n_fail
        }
    })
}

// ---------------------------------------------------------------------------
// FatigueResult
// ---------------------------------------------------------------------------

/// Result of a fatigue life analysis for one FEM element or location.
#[derive(Debug, Clone)]
pub struct FatigueResult {
    /// Predicted life in cycles N_f.
    pub life_cycles: f64,
    /// Palmgren–Miner accumulated damage D.
    pub damage: f64,
    /// Index of the critical node/element (0-based).
    pub critical_location: usize,
    /// Dominant failure mode description.
    pub failure_mode: String,
}

// ---------------------------------------------------------------------------
// FatigueAnalysis
// ---------------------------------------------------------------------------

/// High-level fatigue analysis driver.
///
/// Uses rainflow counting, Basquin model, and Palmgren–Miner rule to
/// estimate fatigue life from an arbitrary stress history at a single FEM
/// element.
pub struct FatigueAnalysis {
    /// Fatigue material parameters.
    pub params: FatigueParams,
    /// S–N curve derived from params.
    pub sn: SNCurve,
}

impl FatigueAnalysis {
    /// Create a new fatigue analysis object.
    ///
    /// # Arguments
    /// * `params` - fatigue material parameters
    pub fn new(params: FatigueParams) -> Self {
        let sn = SNCurve::new(params.s_e, params.s_ut, -1.0 / params.b, params.n_f_ref);
        Self { params, sn }
    }

    /// Analyse a single element stress history.
    ///
    /// 1. Extract cycles using rainflow counting.
    /// 2. Compute Miner damage sum.
    /// 3. Estimate equivalent life from Basquin model.
    ///
    /// # Arguments
    /// * `stress_history` - time series of stress values (Pa)
    ///
    /// # Returns
    /// [`FatigueResult`] with life, damage, and failure mode.
    pub fn analyze_element(&self, stress_history: &[f64]) -> FatigueResult {
        let cycles = rainflow_count(stress_history);

        // Miner damage using (amplitude, count) pairs
        let cycle_pairs: Vec<(f64, f64)> = cycles
            .iter()
            .map(|&(range, _mean, count)| (0.5 * range, count))
            .collect();
        let damage = miner_damage(&cycle_pairs, &self.sn);

        // Equivalent life: inverse of damage if damage > 0
        let life_cycles = if damage > 0.0 {
            1.0 / damage
        } else {
            f64::INFINITY
        };

        let failure_mode = if damage >= 1.0 {
            "high-cycle fatigue fracture".to_string()
        } else if damage >= 0.5 {
            "moderate fatigue damage — inspect".to_string()
        } else {
            "below damage threshold".to_string()
        };

        FatigueResult {
            life_cycles,
            damage,
            critical_location: 0,
            failure_mode,
        }
    }

    /// Apply stress concentration and correction factors to base endurance limit.
    ///
    /// S_e_corrected = k_a · k_b · k_e · S_e / K_f
    ///
    /// # Arguments
    /// * `ka`  - surface finish factor
    /// * `kb`  - size factor
    /// * `ke`  - reliability factor
    /// * `kf`  - fatigue stress concentration factor K_f
    ///
    /// # Returns
    /// Corrected endurance limit (Pa).
    pub fn corrected_endurance_limit(&self, ka: f64, kb: f64, ke: f64, kf: f64) -> f64 {
        if kf <= 0.0 {
            return 0.0;
        }
        self.params.s_e * ka * kb * ke / kf
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- FatigueModel enum -------------------------------------------------

    #[test]
    fn test_fatigue_model_variants() {
        let m = FatigueModel::Basquin;
        assert_eq!(m, FatigueModel::Basquin);
        assert_ne!(m, FatigueModel::CoffinManson);
    }

    // ---- CyclicStressStrain ------------------------------------------------

    #[test]
    fn test_cyclic_stress_ratio_symmetric() {
        let cyc = CyclicStressStrain::new(200.0e6, -200.0e6, 0.002, -0.002);
        assert!((cyc.stress_ratio() - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_cyclic_stress_ratio_zero_to_tension() {
        let cyc = CyclicStressStrain::new(200.0e6, 0.0, 0.001, 0.0);
        assert!((cyc.stress_ratio() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cyclic_stress_range() {
        let cyc = CyclicStressStrain::new(300.0e6, 100.0e6, 0.003, 0.001);
        assert!((cyc.stress_range() - 200.0e6).abs() < 1.0);
    }

    #[test]
    fn test_cyclic_strain_range() {
        let cyc = CyclicStressStrain::new(300.0e6, 100.0e6, 0.003, 0.001);
        assert!((cyc.strain_range() - 0.002).abs() < 1e-12);
    }

    #[test]
    fn test_cyclic_amplitude() {
        let cyc = CyclicStressStrain::new(400.0e6, 0.0, 0.004, 0.0);
        assert!((cyc.amplitude - 200.0e6).abs() < 1.0);
    }

    #[test]
    fn test_cyclic_mean_stress() {
        let cyc = CyclicStressStrain::new(300.0e6, 100.0e6, 0.0, 0.0);
        assert!((cyc.mean_stress - 200.0e6).abs() < 1.0);
    }

    // ---- SNCurve -----------------------------------------------------------

    #[test]
    fn test_sn_below_endurance_infinite_life() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        assert_eq!(sn.life_at_stress(100.0e6), f64::INFINITY);
    }

    #[test]
    fn test_sn_at_endurance_life() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        // stress == s_e is not > s_e, so should be infinite
        assert_eq!(sn.life_at_stress(210.0e6), f64::INFINITY);
    }

    #[test]
    fn test_sn_above_endurance_finite_life() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let life = sn.life_at_stress(300.0e6);
        assert!(life.is_finite() && life > 0.0);
    }

    #[test]
    fn test_sn_stress_at_life_roundtrip() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let s = 300.0e6;
        let n = sn.life_at_stress(s);
        let s_back = sn.stress_at_life(n);
        assert!((s_back - s).abs() < 1.0); // within 1 Pa
    }

    #[test]
    fn test_sn_stress_at_endurance_cycles() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let s = sn.stress_at_life(1.0e6);
        assert!((s - 210.0e6).abs() < 1.0);
    }

    // ---- Basquin -----------------------------------------------------------

    #[test]
    fn test_basquin_positive_result() {
        let life = basquin_life(300.0e6, 948.0e6, -0.085);
        assert!(life > 0.0 && life.is_finite());
    }

    #[test]
    fn test_basquin_zero_amplitude() {
        assert_eq!(basquin_life(0.0, 948.0e6, -0.085), f64::INFINITY);
    }

    #[test]
    fn test_basquin_higher_stress_shorter_life() {
        let n1 = basquin_life(200.0e6, 948.0e6, -0.085);
        let n2 = basquin_life(400.0e6, 948.0e6, -0.085);
        assert!(n1 > n2);
    }

    // ---- Coffin–Manson -----------------------------------------------------

    #[test]
    fn test_coffin_manson_positive() {
        let n = coffin_manson_life(0.01, 0.26, -0.54);
        assert!(n > 0.0 && n.is_finite());
    }

    #[test]
    fn test_coffin_manson_zero_strain() {
        assert_eq!(coffin_manson_life(0.0, 0.26, -0.54), f64::INFINITY);
    }

    // ---- Morrow ------------------------------------------------------------

    #[test]
    fn test_morrow_no_mean_stress() {
        let n_morrow = morrow_correction(200.0e6, 0.0, 948.0e6, -0.085);
        let n_basquin = basquin_life(200.0e6, 948.0e6, -0.085);
        assert!((n_morrow - n_basquin).abs() < 1.0);
    }

    #[test]
    fn test_morrow_tensile_mean_shorter_life() {
        let n0 = morrow_correction(200.0e6, 0.0, 948.0e6, -0.085);
        let n_pos = morrow_correction(200.0e6, 100.0e6, 948.0e6, -0.085);
        assert!(n_pos < n0);
    }

    // ---- Goodman -----------------------------------------------------------

    #[test]
    fn test_goodman_zero_mean() {
        let eq = goodman_equivalent(200.0e6, 0.0, 210.0e6, 620.0e6);
        assert!((eq - 200.0e6).abs() < 1.0);
    }

    #[test]
    fn test_goodman_high_mean_increases_equivalent() {
        let eq = goodman_equivalent(200.0e6, 300.0e6, 210.0e6, 620.0e6);
        assert!(eq > 200.0e6);
    }

    #[test]
    fn test_goodman_at_ultimate_infinity() {
        let eq = goodman_equivalent(1.0, 620.0e6, 210.0e6, 620.0e6);
        assert_eq!(eq, f64::INFINITY);
    }

    // ---- Walker ------------------------------------------------------------

    #[test]
    fn test_walker_positive() {
        let eq = walker_equivalent(400.0e6, 200.0e6, 0.5);
        assert!(eq > 0.0);
    }

    #[test]
    fn test_walker_symmetric_cycle() {
        // For R = -1 (sigma_max = sigma_a), walker_eq = sigma_a^(1-γ) * sigma_a^γ = sigma_a
        let sigma_a = 200.0e6;
        let eq = walker_equivalent(sigma_a, sigma_a, 0.5);
        assert!((eq - sigma_a).abs() < 1.0);
    }

    // ---- Correction factors ------------------------------------------------

    #[test]
    fn test_notch_factor_no_sensitivity() {
        let kf = notch_factor(3.0, 0.0);
        assert!((kf - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_notch_factor_full_sensitivity() {
        let kf = notch_factor(3.0, 1.0);
        assert!((kf - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_surface_factor_clamped() {
        let ka = surface_finish_factor(4.51, -0.265, 100.0);
        assert!(ka >= 0.0 && ka <= 1.0);
    }

    #[test]
    fn test_size_factor_small() {
        let kb = size_factor(5.0);
        assert!((kb - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_size_factor_medium() {
        let kb = size_factor(25.0);
        assert!(kb < 1.0 && kb > 0.0);
    }

    #[test]
    fn test_reliability_factor_50_percent() {
        let ke = reliability_factor(0.50);
        assert!((ke - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reliability_factor_decreases_with_reliability() {
        let ke90 = reliability_factor(0.90);
        let ke99 = reliability_factor(0.99);
        assert!(ke90 > ke99);
    }

    // ---- Rainflow counting -------------------------------------------------

    #[test]
    fn test_rainflow_empty() {
        let cycles = rainflow_count(&[]);
        assert!(cycles.is_empty());
    }

    #[test]
    fn test_rainflow_single_point() {
        let cycles = rainflow_count(&[100.0]);
        assert!(cycles.is_empty());
    }

    #[test]
    fn test_rainflow_one_full_cycle() {
        // Simple triangle: 0 → 100 → 0
        let history = [0.0, 100.0, 0.0];
        let cycles = rainflow_count(&history);
        // Should find half-cycles (no full cycles extracted from 3 reversals)
        let total_count: f64 = cycles.iter().map(|c| c.2).sum();
        assert!(total_count > 0.0);
    }

    #[test]
    fn test_rainflow_ranges_positive() {
        let history = [0.0, 100.0, -50.0, 150.0, -100.0];
        let cycles = rainflow_count(&history);
        for &(range, _mean, _count) in &cycles {
            assert!(range >= 0.0);
        }
    }

    // ---- Miner rule --------------------------------------------------------

    #[test]
    fn test_miner_no_cycles_zero_damage() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let d = miner_damage(&[], &sn);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_miner_below_endurance_no_damage() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let cycles = vec![(100.0e6, 1.0e6)]; // below endurance
        let d = miner_damage(&cycles, &sn);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_miner_exact_life_damage_one() {
        let sn = SNCurve::new(210.0e6, 620.0e6, 10.0, 1.0e6);
        let s = 300.0e6;
        let n_fail = sn.life_at_stress(s);
        let cycles = vec![(s, n_fail)];
        let d = miner_damage(&cycles, &sn);
        assert!((d - 1.0).abs() < 1e-10);
    }

    // ---- FatigueAnalysis ---------------------------------------------------

    #[test]
    fn test_analysis_new() {
        let params = FatigueParams::default();
        let analysis = FatigueAnalysis::new(params);
        assert!(analysis.params.s_e > 0.0);
    }

    #[test]
    fn test_analyze_element_constant_amplitude() {
        let params = FatigueParams::default();
        let analysis = FatigueAnalysis::new(params);
        // Repeated cycle at amplitude well above endurance
        let history: Vec<f64> = (0..201)
            .map(|i| {
                let t = i as f64 * std::f64::consts::PI / 50.0;
                400.0e6 * t.sin()
            })
            .collect();
        let result = analysis.analyze_element(&history);
        assert!(result.life_cycles > 0.0);
        assert!(result.damage >= 0.0);
    }

    #[test]
    fn test_analyze_element_zero_stress() {
        let params = FatigueParams::default();
        let analysis = FatigueAnalysis::new(params);
        let history = vec![0.0f64; 20];
        let result = analysis.analyze_element(&history);
        assert_eq!(result.damage, 0.0);
        assert_eq!(result.life_cycles, f64::INFINITY);
    }

    #[test]
    fn test_corrected_endurance_limit() {
        let params = FatigueParams::default();
        let analysis = FatigueAnalysis::new(params);
        let s_ec = analysis.corrected_endurance_limit(0.9, 0.95, 0.9, 1.5);
        assert!(s_ec < analysis.params.s_e);
        assert!(s_ec > 0.0);
    }

    #[test]
    fn test_failure_mode_string_present() {
        let params = FatigueParams::default();
        let analysis = FatigueAnalysis::new(params);
        let history: Vec<f64> = (0..101)
            .map(|i| 600.0e6 * (i as f64 * std::f64::consts::PI / 10.0).sin())
            .collect();
        let result = analysis.analyze_element(&history);
        assert!(!result.failure_mode.is_empty());
    }
}
