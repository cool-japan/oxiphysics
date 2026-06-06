// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Chassis resonance and noise-vibration-harshness (NVH) analysis.
//!
//! Provides SDOF vibration models, frequency response functions, road roughness
//! PSD (ISO 8608), tire envelopment filter, dynamic vibration absorber design,
//! and insertion loss calculation.
//!
//! # Overview
//!
//! - [`ChassisMode`] — single resonant mode (frequency, damping, mode shape, description).
//! - [`ChassisModel`] — multi-DOF chassis model with stiffness and damping matrices.
//! - [`NvhSpectrum`] — amplitude spectrum with peak detection and octave-band analysis.
//! - [`VibrationAbsorber`] — dynamic vibration absorber (DVA) optimal tuning.
//! - [`undamped_natural_frequency`] — `ω_n = √(k/m)`.
//! - [`damped_natural_frequency`] — `ω_d = ω_n · √(1 − ζ²)`.
//! - [`sdof_frequency_response`] — SDOF dynamic magnification factor |H(jω)|.
//! - [`resonance_amplitude`] — peak amplitude at resonance for lightly damped SDOF.
//! - [`resonance_frequency_damped`] — `ω_r = ω_n · √(1 − 2ζ²)`.
//! - [`road_roughness_psd`] — ISO 8608 road roughness power spectral density.
//! - [`tire_envelopment_filter`] — tire contact-patch envelopment filter gain.
//! - [`insertion_loss_db`] — treatment insertion loss in dB.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// ChassisMode
// ─────────────────────────────────────────────────────────────────────────────

/// A single structural resonant mode of the chassis.
#[derive(Debug, Clone)]
pub struct ChassisMode {
    /// Natural frequency \[Hz\].
    pub frequency_hz: f64,
    /// Damping ratio ζ \[-\].
    pub damping_ratio: f64,
    /// Mode shape vector (normalised participation factors per DOF).
    pub mode_shape: Vec<f64>,
    /// Human-readable description of this mode.
    pub description: String,
}

// ─────────────────────────────────────────────────────────────────────────────
// ChassisModel
// ─────────────────────────────────────────────────────────────────────────────

/// Multi-DOF chassis model storing stiffness and damping matrices.
#[derive(Debug, Clone)]
pub struct ChassisModel {
    /// Total sprung mass \[kg\].
    pub mass: f64,
    /// Principal moments of inertia `[Ixx, Iyy, Izz]` \[kg·m²\].
    pub inertia: [f64; 3],
    /// Stiffness matrix (n_dof × n_dof), stored row-major.
    pub stiffness_matrix: Vec<Vec<f64>>,
    /// Damping matrix (n_dof × n_dof), stored row-major.
    pub damping_matrix: Vec<Vec<f64>>,
    /// Number of degrees of freedom.
    pub n_dof: usize,
}

impl ChassisModel {
    /// Construct a new chassis model with `n_dof` DOF and total mass `mass` \[kg\].
    ///
    /// All matrix entries default to zero; populate with application-specific values.
    pub fn new(n_dof: usize, mass: f64) -> Self {
        let mat = vec![vec![0.0_f64; n_dof]; n_dof];
        Self {
            mass,
            inertia: [0.0; 3],
            stiffness_matrix: mat.clone(),
            damping_matrix: mat,
            n_dof,
        }
    }

    /// Approximate undamped natural frequencies \[rad/s\] for each DOF.
    ///
    /// Assumes a diagonal stiffness matrix; off-diagonal coupling is ignored in
    /// this simplified implementation.
    pub fn natural_frequencies(&self) -> Vec<f64> {
        (0..self.n_dof)
            .map(|i| {
                let k = self.stiffness_matrix[i][i].max(0.0);
                let m = (self.mass / self.n_dof as f64).max(1e-15);
                (k / m).sqrt()
            })
            .collect()
    }

    /// Critical damping coefficient \[N·s/m\] for a mode with natural frequency
    /// `omega_n` \[rad/s\].
    ///
    /// `c_cr = 2 · m · ω_n`
    pub fn critical_damping(&self, omega_n: f64) -> f64 {
        2.0 * self.mass * omega_n
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions — vibration fundamentals
// ─────────────────────────────────────────────────────────────────────────────

/// Undamped natural frequency \[rad/s\].
///
/// `ω_n = √(k / m)`
pub fn undamped_natural_frequency(k: f64, m: f64) -> f64 {
    (k / m.max(1e-15)).sqrt()
}

/// Damped natural frequency \[rad/s\].
///
/// `ω_d = ω_n · √(1 − ζ²)`
///
/// Returns 0 for over-critically-damped systems (ζ ≥ 1).
pub fn damped_natural_frequency(omega_n: f64, zeta: f64) -> f64 {
    if zeta >= 1.0 {
        return 0.0;
    }
    omega_n * (1.0 - zeta * zeta).sqrt()
}

/// SDOF dynamic magnification factor |H(jω)| at excitation frequency `omega` \[rad/s\].
///
/// `|H| = 1 / √[(1 − r²)² + (2ζr)²]`
///
/// where `r = ω / ω_n`.
pub fn sdof_frequency_response(omega_n: f64, zeta: f64, omega: f64) -> f64 {
    if omega_n < 1e-12 {
        return 1.0;
    }
    let r = omega / omega_n;
    let denom = ((1.0 - r * r).powi(2) + (2.0 * zeta * r).powi(2)).sqrt();
    1.0 / denom.max(1e-15)
}

/// Peak amplitude at resonance for an underdamped SDOF (ζ < 1/√2 ≈ 0.707).
///
/// `A_peak = 1 / (2ζ · √(1 − ζ²))`
///
/// Returns `f64::INFINITY` for ζ = 0.
pub fn resonance_amplitude(zeta: f64) -> f64 {
    if zeta <= 0.0 {
        return f64::INFINITY;
    }
    if zeta >= 1.0 {
        return 1.0;
    }
    1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt())
}

/// Frequency of peak response (resonance frequency) for a damped SDOF \[rad/s\].
///
/// `ω_r = ω_n · √(1 − 2ζ²)`
///
/// Returns 0 if ζ ≥ 1/√2 (no resonance peak exists).
pub fn resonance_frequency_damped(omega_n: f64, zeta: f64) -> f64 {
    let threshold = (0.5_f64).sqrt(); // 1/√2
    if zeta >= threshold {
        return 0.0;
    }
    omega_n * (1.0 - 2.0 * zeta * zeta).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// NvhSpectrum
// ─────────────────────────────────────────────────────────────────────────────

/// NVH amplitude spectrum (frequencies and corresponding amplitudes).
#[derive(Debug, Clone)]
pub struct NvhSpectrum {
    /// Frequency vector \[Hz\].
    pub frequencies: Vec<f64>,
    /// Amplitude vector (linear, e.g., m/s² or Pa).
    pub amplitudes: Vec<f64>,
}

impl NvhSpectrum {
    /// Construct an NVH spectrum from parallel frequency and amplitude vectors.
    pub fn new(freqs: Vec<f64>, amps: Vec<f64>) -> Self {
        Self {
            frequencies: freqs,
            amplitudes: amps,
        }
    }

    /// Frequency \[Hz\] at which the amplitude is maximum.
    ///
    /// Returns `0.0` if the spectrum is empty.
    pub fn peak_frequency(&self) -> f64 {
        if self.amplitudes.is_empty() {
            return 0.0;
        }
        let idx = self
            .amplitudes
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        self.frequencies.get(idx).copied().unwrap_or(0.0)
    }

    /// Overall RMS level computed as `√(Σ aᵢ²)`.
    pub fn overall_level(&self) -> f64 {
        let sum_sq: f64 = self.amplitudes.iter().map(|a| a * a).sum();
        sum_sq.sqrt()
    }

    /// Octave-band RMS levels.
    ///
    /// Returns RMS amplitude in standard octave bands centred at 31.5, 63, 125,
    /// 250, 500, 1000, 2000, 4000, 8000, 16000 Hz.
    pub fn octave_band_levels(&self) -> Vec<f64> {
        let centres = [
            31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0,
        ];
        centres
            .iter()
            .map(|&fc| {
                let f_lo = fc / 2.0_f64.sqrt();
                let f_hi = fc * 2.0_f64.sqrt();
                let sum_sq: f64 = self
                    .frequencies
                    .iter()
                    .zip(self.amplitudes.iter())
                    .filter(|&(&f, _)| f >= f_lo && f < f_hi)
                    .map(|(_, a)| a * a)
                    .sum();
                sum_sq.sqrt()
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Road roughness PSD
// ─────────────────────────────────────────────────────────────────────────────

/// ISO 8608 road roughness power spectral density \[m²/(cycles/m)\].
///
/// `Gd(n) = Gd(n₀) · (n / n₀)^(−w)`
///
/// where `n₀ = 0.1 cycles/m` is the reference spatial frequency,
/// `w = 2` (road roughness exponent), and
/// `Gd(n₀)` depends on `road_class`:
///
/// | class | description | Gd(n₀) \[m²/(cycles/m)\] |
/// |-------|-------------|------------------------|
/// | 1     | very good   | 1e-6                   |
/// | 2     | good        | 4e-6                   |
/// | 3     | average     | 16e-6                  |
/// | 4     | poor        | 64e-6                  |
/// | other | very poor   | 256e-6                 |
///
/// * `frequency` — spatial frequency \[cycles/m\]
/// * `road_class` — road quality class (1–4+)
pub fn road_roughness_psd(frequency: f64, road_class: u8) -> f64 {
    let n0 = 0.1_f64; // reference spatial frequency [cycles/m]
    let gd_n0 = match road_class {
        1 => 1e-6,
        2 => 4e-6,
        3 => 16e-6,
        4 => 64e-6,
        _ => 256e-6,
    };
    gd_n0 * (frequency.max(1e-9) / n0).powi(-2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tire envelopment filter
// ─────────────────────────────────────────────────────────────────────────────

/// Tire contact-patch envelopment filter gain (low-pass) \[-\].
///
/// Models the spatial averaging of road inputs over the contact patch.
///
/// `H(f) = sinc(π · f · L_c / v)` where `L_c` is the contact patch length and
/// `v` is vehicle speed.
///
/// Returns 1.0 if speed is near zero (DC pass-through).
///
/// * `freq` — temporal frequency \[Hz\]
/// * `contact_length` — tire contact patch length \[m\]
/// * `speed` — vehicle speed \[m/s\]
pub fn tire_envelopment_filter(freq: f64, contact_length: f64, speed: f64) -> f64 {
    if speed < 1e-6 {
        return 1.0;
    }
    let x = PI * freq * contact_length / speed;
    if x.abs() < 1e-9 {
        1.0
    } else {
        (x.sin() / x).abs()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VibrationAbsorber
// ─────────────────────────────────────────────────────────────────────────────

/// Dynamic vibration absorber (DVA / tuned mass damper) parameters.
#[derive(Debug, Clone)]
pub struct VibrationAbsorber {
    /// Absorber mass \[kg\].
    pub mass: f64,
    /// Absorber stiffness \[N/m\].
    pub stiffness: f64,
    /// Absorber damping coefficient \[N·s/m\].
    pub damping: f64,
    /// Tuning ratio `f = ω_a / ω_n` \[-\].
    pub tuning_ratio: f64,
}

impl VibrationAbsorber {
    /// Optimal tuning ratio for a Den Hartog DVA.
    ///
    /// `f_opt = 1 / (1 + μ)`
    ///
    /// where `μ = m_absorber / m_primary` is the mass ratio.
    pub fn optimal_tuning_ratio(mass_ratio: f64) -> f64 {
        1.0 / (1.0 + mass_ratio)
    }

    /// Optimal damping ratio for a Den Hartog DVA.
    ///
    /// `ζ_opt = √(3μ / (8(1 + μ)³))`
    pub fn optimal_damping(mass_ratio: f64) -> f64 {
        let mu = mass_ratio.max(0.0);
        let numerator = 3.0 * mu;
        let denominator = 8.0 * (1.0 + mu).powi(3);
        (numerator / denominator).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Insertion loss
// ─────────────────────────────────────────────────────────────────────────────

/// Insertion loss \[dB\] of an acoustic treatment.
///
/// `IL = 20 · log10(L_original / L_treated)`
///
/// * `original_level` — untreated sound level (linear, e.g., Pa or m/s²)
/// * `treated_level` — treated sound level (same units)
pub fn insertion_loss_db(original_level: f64, treated_level: f64) -> f64 {
    if treated_level < 1e-15 {
        return f64::INFINITY;
    }
    20.0 * (original_level / treated_level).log10()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── undamped_natural_frequency ───────────────────────────────────────

    #[test]
    fn undamped_freq_formula() {
        let omega = undamped_natural_frequency(1000.0, 10.0);
        let expected = (100.0_f64).sqrt(); // sqrt(1000/10)
        assert!((omega - expected).abs() < 1e-9);
    }

    #[test]
    fn undamped_freq_positive() {
        assert!(undamped_natural_frequency(500.0, 5.0) > 0.0);
    }

    #[test]
    fn undamped_freq_doubles_when_k_quadruples() {
        let w1 = undamped_natural_frequency(100.0, 1.0);
        let w2 = undamped_natural_frequency(400.0, 1.0);
        assert!((w2 / w1 - 2.0).abs() < EPS);
    }

    #[test]
    fn undamped_freq_halves_when_m_quadruples() {
        let w1 = undamped_natural_frequency(100.0, 1.0);
        let w2 = undamped_natural_frequency(100.0, 4.0);
        assert!((w2 / w1 - 0.5).abs() < EPS);
    }

    // ── damped_natural_frequency ─────────────────────────────────────────

    #[test]
    fn damped_freq_less_than_undamped() {
        let omega_n = 10.0;
        let zeta = 0.2;
        let omega_d = damped_natural_frequency(omega_n, zeta);
        assert!(omega_d < omega_n);
    }

    #[test]
    fn damped_freq_equals_undamped_at_zero_damping() {
        let omega_d = damped_natural_frequency(10.0, 0.0);
        assert!((omega_d - 10.0).abs() < EPS);
    }

    #[test]
    fn damped_freq_zero_at_critical_damping() {
        let omega_d = damped_natural_frequency(10.0, 1.0);
        assert_eq!(omega_d, 0.0);
    }

    #[test]
    fn damped_freq_decreases_with_zeta() {
        let w1 = damped_natural_frequency(10.0, 0.1);
        let w2 = damped_natural_frequency(10.0, 0.5);
        assert!(w2 < w1);
    }

    // ── sdof_frequency_response ──────────────────────────────────────────

    #[test]
    fn sdof_response_one_at_dc() {
        // ω → 0 → r → 0 → |H| → 1
        let h = sdof_frequency_response(10.0, 0.1, 0.0);
        assert!((h - 1.0).abs() < EPS);
    }

    #[test]
    fn sdof_response_greater_than_one_at_resonance() {
        // At ω ≈ ω_n, underdamped → |H| > 1
        let omega_n = 10.0;
        let zeta = 0.1;
        let h = sdof_frequency_response(omega_n, zeta, omega_n);
        assert!(h > 1.0);
    }

    #[test]
    fn sdof_response_positive() {
        let h = sdof_frequency_response(10.0, 0.3, 8.0);
        assert!(h > 0.0);
    }

    // ── resonance_amplitude ──────────────────────────────────────────────

    #[test]
    fn resonance_amplitude_increases_as_zeta_decreases() {
        let a1 = resonance_amplitude(0.3);
        let a2 = resonance_amplitude(0.1);
        assert!(a2 > a1);
    }

    #[test]
    fn resonance_amplitude_greater_than_one_for_underdamped() {
        let a = resonance_amplitude(0.1);
        assert!(a > 1.0);
    }

    #[test]
    fn resonance_amplitude_infinite_at_zero_damping() {
        assert_eq!(resonance_amplitude(0.0), f64::INFINITY);
    }

    // ── resonance_frequency_damped ───────────────────────────────────────

    #[test]
    fn resonance_freq_less_than_undamped() {
        let omega_n = 10.0;
        let omega_r = resonance_frequency_damped(omega_n, 0.1);
        assert!(omega_r < omega_n);
    }

    #[test]
    fn resonance_freq_zero_at_high_damping() {
        // ζ ≥ 1/√2 → no resonance peak
        let omega_r = resonance_frequency_damped(10.0, 0.75);
        assert_eq!(omega_r, 0.0);
    }

    #[test]
    fn resonance_freq_positive_for_low_zeta() {
        assert!(resonance_frequency_damped(10.0, 0.1) > 0.0);
    }

    // ── NvhSpectrum ──────────────────────────────────────────────────────

    #[test]
    fn nvh_peak_frequency_correct() {
        let freqs = vec![10.0, 20.0, 30.0];
        let amps = vec![1.0, 5.0, 2.0];
        let spec = NvhSpectrum::new(freqs, amps);
        assert!((spec.peak_frequency() - 20.0).abs() < EPS);
    }

    #[test]
    fn nvh_overall_level_positive() {
        let spec = NvhSpectrum::new(vec![10.0, 20.0], vec![1.0, 1.0]);
        assert!(spec.overall_level() > 0.0);
    }

    #[test]
    fn nvh_overall_level_formula() {
        let spec = NvhSpectrum::new(vec![10.0, 20.0], vec![3.0, 4.0]);
        assert!((spec.overall_level() - 5.0).abs() < EPS); // sqrt(9+16)=5
    }

    #[test]
    fn nvh_octave_bands_count() {
        let spec = NvhSpectrum::new(vec![1000.0], vec![1.0]);
        assert_eq!(spec.octave_band_levels().len(), 10);
    }

    #[test]
    fn nvh_empty_peak_returns_zero() {
        let spec = NvhSpectrum::new(vec![], vec![]);
        assert_eq!(spec.peak_frequency(), 0.0);
    }

    // ── road_roughness_psd ───────────────────────────────────────────────

    #[test]
    fn road_roughness_psd_positive() {
        assert!(road_roughness_psd(1.0, 2) > 0.0);
    }

    #[test]
    fn road_roughness_psd_decreases_with_frequency() {
        let psd1 = road_roughness_psd(0.1, 2);
        let psd2 = road_roughness_psd(1.0, 2);
        assert!(psd2 < psd1);
    }

    #[test]
    fn road_roughness_psd_increases_with_class() {
        let psd_good = road_roughness_psd(1.0, 1);
        let psd_poor = road_roughness_psd(1.0, 4);
        assert!(psd_poor > psd_good);
    }

    // ── tire_envelopment_filter ──────────────────────────────────────────

    #[test]
    fn tire_filter_unity_at_dc() {
        let h = tire_envelopment_filter(0.0, 0.2, 20.0);
        assert!((h - 1.0).abs() < EPS);
    }

    #[test]
    fn tire_filter_unity_at_zero_speed() {
        let h = tire_envelopment_filter(100.0, 0.2, 0.0);
        assert!((h - 1.0).abs() < EPS);
    }

    #[test]
    fn tire_filter_positive() {
        let h = tire_envelopment_filter(50.0, 0.2, 20.0);
        assert!((0.0..=1.0 + EPS).contains(&h));
    }

    // ── VibrationAbsorber ────────────────────────────────────────────────

    #[test]
    fn optimal_tuning_ratio_mu_01() {
        // μ = 0.1 → f_opt = 1/(1+0.1) ≈ 0.9091
        let f = VibrationAbsorber::optimal_tuning_ratio(0.1);
        assert!((f - 1.0 / 1.1).abs() < 1e-9);
    }

    #[test]
    fn optimal_tuning_ratio_less_than_one() {
        let f = VibrationAbsorber::optimal_tuning_ratio(0.05);
        assert!(f < 1.0);
    }

    #[test]
    fn optimal_tuning_ratio_approaches_one_at_zero_mu() {
        let f = VibrationAbsorber::optimal_tuning_ratio(0.0);
        assert!((f - 1.0).abs() < EPS);
    }

    #[test]
    fn optimal_damping_positive() {
        assert!(VibrationAbsorber::optimal_damping(0.1) > 0.0);
    }

    #[test]
    fn optimal_damping_increases_with_mu() {
        let d1 = VibrationAbsorber::optimal_damping(0.05);
        let d2 = VibrationAbsorber::optimal_damping(0.2);
        assert!(d2 > d1);
    }

    // ── insertion_loss_db ────────────────────────────────────────────────

    #[test]
    fn insertion_loss_zero_when_same() {
        assert!((insertion_loss_db(1.0, 1.0)).abs() < EPS);
    }

    #[test]
    fn insertion_loss_positive_when_treated_less() {
        assert!(insertion_loss_db(2.0, 1.0) > 0.0);
    }

    #[test]
    fn insertion_loss_6db_for_half_amplitude() {
        // 20*log10(2/1) = 6.02 dB
        let il = insertion_loss_db(2.0, 1.0);
        assert!((il - 20.0 * 2.0_f64.log10()).abs() < 1e-9);
    }

    // ── ChassisModel ─────────────────────────────────────────────────────

    #[test]
    fn chassis_model_natural_frequencies_count() {
        let mut model = ChassisModel::new(4, 1500.0);
        model.stiffness_matrix[0][0] = 20000.0;
        model.stiffness_matrix[1][1] = 22000.0;
        model.stiffness_matrix[2][2] = 20000.0;
        model.stiffness_matrix[3][3] = 22000.0;
        let freqs = model.natural_frequencies();
        assert_eq!(freqs.len(), 4);
    }

    #[test]
    fn chassis_model_natural_frequencies_positive() {
        let mut model = ChassisModel::new(2, 1000.0);
        model.stiffness_matrix[0][0] = 10000.0;
        model.stiffness_matrix[1][1] = 15000.0;
        let freqs = model.natural_frequencies();
        assert!(freqs.iter().all(|&f| f > 0.0));
    }

    #[test]
    fn chassis_critical_damping_positive() {
        let model = ChassisModel::new(1, 1000.0);
        let cc = model.critical_damping(10.0);
        assert!(cc > 0.0);
    }
}
