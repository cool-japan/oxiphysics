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

/// Reduce a sampled history to its turning points (peaks and valleys).
///
/// Non-finite samples are ignored.  Runs of equal consecutive samples are
/// collapsed to one sample first, so a flat peak or valley (plateau) is kept as
/// a single reversal instead of being lost.  The first and last samples are
/// always retained.
fn turning_points(history: &[f64]) -> Vec<f64> {
    let mut samples: Vec<f64> = Vec::with_capacity(history.len());
    for &s in history {
        if s.is_finite() && samples.last() != Some(&s) {
            samples.push(s);
        }
    }
    let (Some(&first), Some(&last)) = (samples.first(), samples.last()) else {
        return samples;
    };
    if samples.len() < 3 {
        return samples;
    }
    let mut turning: Vec<f64> = Vec::with_capacity(samples.len());
    turning.push(first);
    for w in samples.windows(3) {
        // Consecutive samples differ, so the direction of each step is a
        // plain comparison: a reversal is a change of direction.
        if (w[1] > w[0]) != (w[2] > w[1]) {
            turning.push(w[1]);
        }
    }
    turning.push(last);
    turning
}

/// Rainflow cycle counting (ASTM E1049-85 three-point algorithm).
///
/// Reduces the history to its reversals (plateaus count as a single reversal,
/// non-finite samples are ignored) and counts cycles with the stack algorithm
/// of ASTM E1049-85 §5.4.4.  Whenever the latest range `X` is at least as large
/// as the preceding range `Y`:
///
/// - if `Y` contains the starting point of the history (the first point still
///   on the stack), `Y` is counted as a **half cycle** and only that starting
///   point is discarded;
/// - otherwise `Y` is a closed cycle: it is counted as **one full cycle** and
///   both of its points are discarded.
///
/// The ranges left on the stack at the end of the history (the residue) are
/// counted as half cycles.  Every reversal-to-reversal step is therefore
/// accounted for exactly once: the counts sum to `(reversals − 1) / 2` and
/// `Σ 2·count·range` equals the total variation of the reversal sequence.
///
/// For example `[−2, 1, −3, 5, −1, 3, −4, 4, −2]` gives ranges 3 (0.5 cycle),
/// 4 (1.5), 6 (0.5), 8 (1.0) and 9 (0.5), as in the ASTM E1049 example.
///
/// # Arguments
/// * `stress_history` - stress time series (Pa)
///
/// # Returns
/// `Vec<(range, mean, count)>` — each entry is a counted cycle, with `count`
/// equal to 0.5 for a half cycle and 1.0 for a full cycle.  The vector is empty
/// when the history has fewer than two distinct finite samples.
pub fn rainflow_count(stress_history: &[f64]) -> Vec<(f64, f64, f64)> {
    let reversals = turning_points(stress_history);
    if reversals.len() < 2 {
        return Vec::new();
    }

    let mut stack: Vec<f64> = Vec::with_capacity(reversals.len());
    // Position in `stack` of the starting point S of the history.  Discarding
    // S just advances this index, so no element is ever shifted.
    let mut start = 0usize;
    let mut cycles: Vec<(f64, f64, f64)> = Vec::new();

    for &point in &reversals {
        stack.push(point);
        while stack.len() - start >= 3 {
            let n = stack.len();
            let (s0, s1, s2) = (stack[n - 3], stack[n - 2], stack[n - 1]);
            let y = (s1 - s0).abs();
            let x = (s2 - s1).abs();
            if x < y {
                break;
            }
            let mean = 0.5 * (s0 + s1);
            if n - 3 == start {
                cycles.push((y, mean, 0.5));
                start += 1;
            } else {
                cycles.push((y, mean, 1.0));
                stack.truncate(n - 3);
                stack.push(s2);
            }
        }
    }

    // The residue is counted as half cycles.
    for w in stack[start..].windows(2) {
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
        assert!((0.0..=1.0).contains(&ka));
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

    /// Total count per (integer) range, as a sorted list of `(range, count)`.
    fn counts_by_range(cycles: &[(f64, f64, f64)]) -> Vec<(i64, f64)> {
        let mut table: Vec<(i64, f64)> = Vec::new();
        for &(range, _mean, count) in cycles {
            let key = range.round() as i64;
            assert!(
                (range - key as f64).abs() < 1e-9,
                "non-integer range {range}"
            );
            match table.iter_mut().find(|(k, _)| *k == key) {
                Some(entry) => entry.1 += count,
                None => table.push((key, count)),
            }
        }
        table.sort_by_key(|&(k, _)| k);
        table
    }

    #[test]
    fn test_rainflow_astm_e1049_example() {
        // ASTM E1049-85 Fig. 6: -2, 1, -3, 5, -1, 3, -4, 4, -2.
        let history = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0];
        let cycles = rainflow_count(&history);
        assert_eq!(
            counts_by_range(&cycles),
            vec![(3, 0.5), (4, 1.5), (6, 0.5), (8, 1.0), (9, 0.5)]
        );
        let total: f64 = cycles.iter().map(|c| c.2).sum();
        assert!((total - 4.0).abs() < 1e-12, "8 reversal steps = 4 cycles");
    }

    #[test]
    fn test_rainflow_start_point_range_is_half_cycle() {
        // [0, 5, -10, 10]: the 0-5 range contains the starting point and must be
        // a half cycle, so no range is lost (the old code reported 5:1.0, 20:0.5).
        let cycles = rainflow_count(&[0.0, 5.0, -10.0, 10.0]);
        assert_eq!(
            counts_by_range(&cycles),
            vec![(5, 0.5), (15, 0.5), (20, 0.5)]
        );
    }

    #[test]
    fn test_rainflow_means_of_counted_cycles() {
        let cycles = rainflow_count(&[0.0, 5.0, -10.0, 10.0]);
        let mean_of = |range: f64| {
            cycles
                .iter()
                .find(|c| (c.0 - range).abs() < 1e-12)
                .map(|c| c.1)
        };
        assert_eq!(mean_of(5.0), Some(2.5));
        assert_eq!(mean_of(15.0), Some(-2.5));
        assert_eq!(mean_of(20.0), Some(0.0));
    }

    #[test]
    fn test_rainflow_plateau_peak_is_kept() {
        // The flat top must count as one peak: 0 -> 10 -> 0 is one full cycle
        // expressed as two half cycles of range 10 about the mean 5.  The old
        // reversal test dropped the plateau and reported nothing for [0, 10, 10, 0].
        for history in [
            vec![0.0, 10.0, 10.0, 0.0],
            vec![0.0, 10.0, 10.0, 10.0, 10.0, 0.0],
            vec![0.0, 0.0, 10.0, 10.0, 0.0, 0.0],
        ] {
            let cycles = rainflow_count(&history);
            assert_eq!(cycles.len(), 2, "history {history:?}");
            for &(range, mean, count) in &cycles {
                assert_eq!(
                    (range, mean, count),
                    (10.0, 5.0, 0.5),
                    "history {history:?}"
                );
            }
        }
    }

    #[test]
    fn test_rainflow_plateau_valley_and_interior_plateau() {
        // 3,3,0,0,0,10,10 collapses to 3,0,10: a half cycle of range 3 at the
        // start point and a residue half cycle of range 10.
        let cycles = rainflow_count(&[3.0, 3.0, 0.0, 0.0, 0.0, 10.0, 10.0]);
        assert_eq!(counts_by_range(&cycles), vec![(3, 0.5), (10, 0.5)]);
        // A monotone ramp with a plateau in the middle is one half cycle.
        let ramp = rainflow_count(&[0.0, 1.0, 2.0, 2.0, 3.0, 4.0]);
        assert_eq!(ramp, vec![(4.0, 2.0, 0.5)]);
    }

    #[test]
    fn test_rainflow_degenerate_histories() {
        assert!(rainflow_count(&[5.0, 5.0, 5.0, 5.0]).is_empty());
        assert!(rainflow_count(&[f64::NAN, 1.0, f64::INFINITY]).is_empty());
        assert_eq!(rainflow_count(&[1.0, 2.0]), vec![(1.0, 1.5, 0.5)]);
        assert_eq!(rainflow_count(&[2.0, 2.0, 1.0]), vec![(1.0, 1.5, 0.5)]);
    }

    #[test]
    fn test_rainflow_ignores_non_finite_samples() {
        let clean = rainflow_count(&[0.0, 10.0, -4.0, 6.0]);
        let dirty = rainflow_count(&[0.0, f64::NAN, 10.0, -4.0, f64::INFINITY, 6.0]);
        assert_eq!(clean, dirty);
    }

    #[test]
    fn test_rainflow_conserves_counts_and_total_variation() {
        // Deterministic pseudo-random history with deliberate plateaus.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 33) as f64) / ((1u64 << 31) as f64)
        };
        let mut history = Vec::new();
        let mut value = 0.0_f64;
        for _ in 0..400 {
            // One draw in five repeats the previous sample (a plateau).
            if next() >= 0.2 {
                value = 100.0 * (next() - 0.5);
            }
            history.push(value);
        }
        let turning = turning_points(&history);
        let cycles = rainflow_count(&history);
        let total_count: f64 = cycles.iter().map(|c| c.2).sum();
        assert!(
            (total_count - 0.5 * (turning.len() as f64 - 1.0)).abs() < 1e-9,
            "count {total_count} vs {} reversals",
            turning.len()
        );
        let variation: f64 = turning.windows(2).map(|w| (w[1] - w[0]).abs()).sum();
        let counted: f64 = cycles.iter().map(|c| 2.0 * c.2 * c.0).sum();
        assert!(
            (variation - counted).abs() < 1e-8 * variation,
            "total variation {variation} vs counted {counted}"
        );
        assert!(
            cycles
                .iter()
                .all(|c| c.0 > 0.0 && (c.2 == 0.5 || c.2 == 1.0))
        );
    }

    #[test]
    fn test_analyze_element_counts_start_point_half_cycle() {
        let analysis = FatigueAnalysis::new(FatigueParams::default());
        let history = [0.0, 500.0e6, -500.0e6, 500.0e6];
        let result = analysis.analyze_element(&history);
        // ASTM counts: ranges 500, 1000, 1000 MPa, each half a cycle.
        let pairs = [(250.0e6, 0.5), (500.0e6, 0.5), (500.0e6, 0.5)];
        let expected = miner_damage(&pairs, &analysis.sn);
        assert!(expected > 0.0);
        assert!((result.damage - expected).abs() <= 1e-12 * expected);
        // The old (start-point blind) count lost the 1000 MPa half cycle and
        // so under-predicted the damage.
        let old_pairs = [(250.0e6, 1.0), (500.0e6, 0.5)];
        assert!(result.damage > miner_damage(&old_pairs, &analysis.sn));
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
