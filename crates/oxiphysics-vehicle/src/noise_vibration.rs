// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Noise, vibration, and harshness (NVH) analysis for vehicle systems.
//!
//! Provides single-degree-of-freedom (SDOF) vibration models, frequency
//! response functions, A-weighted sound pressure levels, vibration isolator
//! transmissibility, road roughness PSD, and engine order analysis.
//!
//! # Overview
//!
//! - [`FrequencyResponse`] — FRF data: frequencies, magnitudes, and phases.
//! - [`Sdof`] — mass-spring-damper SDOF system with analytical FRF.
//! - [`NoiseSpl`] — octave-band SPL data with A-weighting.
//! - [`VibrationIsolator`] — transmissibility and insertion loss of an isolator mount.
//! - [`RoadRoughness`] — ISO 8608 road PSD and RMS amplitude.
//! - [`HarmonicOrder`] — engine harmonic-order frequency computation.
//! - [`octave_band_center_frequencies`] — centre frequencies for n octave bands.
//! - [`a_weighting_db`] — A-weighting correction in dB at a given frequency.
//! - [`sone_to_phon`] — loudness conversion from sone to phon.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// FrequencyResponse
// ---------------------------------------------------------------------------

/// Frequency response function (FRF) data.
///
/// Stores the magnitude and phase spectrum computed at a set of frequencies.
#[derive(Debug, Clone)]
pub struct FrequencyResponse {
    /// Frequency vector \[Hz\].
    pub frequencies: Vec<f64>,
    /// Magnitude values (linear, not dB).
    pub magnitudes: Vec<f64>,
    /// Phase values \[rad\].
    pub phases: Vec<f64>,
}

impl FrequencyResponse {
    /// Construct an FRF from parallel frequency, magnitude, and phase vectors.
    pub fn new(frequencies: Vec<f64>, magnitudes: Vec<f64>, phases: Vec<f64>) -> Self {
        Self {
            frequencies,
            magnitudes,
            phases,
        }
    }

    /// Frequency \[Hz\] at which the magnitude is largest.
    ///
    /// Returns `0.0` if the magnitude vector is empty.
    pub fn peak_frequency(&self) -> f64 {
        if self.magnitudes.is_empty() {
            return 0.0;
        }
        let idx = self
            .magnitudes
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        *self.frequencies.get(idx).unwrap_or(&0.0)
    }

    /// Peak magnitude value.
    ///
    /// Returns `0.0` for an empty spectrum.
    pub fn peak_magnitude(&self) -> f64 {
        self.magnitudes
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }
}

// ---------------------------------------------------------------------------
// Sdof
// ---------------------------------------------------------------------------

/// Single-degree-of-freedom (SDOF) mass-spring-damper system.
///
/// Equation of motion: m·ẍ + c·ẋ + k·x = F(t)
#[derive(Debug, Clone)]
pub struct Sdof {
    /// Mass m \[kg\].
    pub mass: f64,
    /// Viscous damping coefficient c \[N·s/m\].
    pub damping: f64,
    /// Spring stiffness k \[N/m\].
    pub stiffness: f64,
}

impl Sdof {
    /// Construct an SDOF system.
    pub fn new(mass: f64, damping: f64, stiffness: f64) -> Self {
        Self {
            mass,
            damping,
            stiffness,
        }
    }

    /// Undamped natural frequency ω_n = √(k/m) \[rad/s\].
    pub fn natural_frequency_rad(&self) -> f64 {
        (self.stiffness / self.mass).sqrt()
    }

    /// Undamped natural frequency f_n \[Hz\].
    pub fn natural_frequency(&self) -> f64 {
        self.natural_frequency_rad() / (2.0 * PI)
    }

    /// Critical damping coefficient c_cr = 2√(km).
    pub fn critical_damping(&self) -> f64 {
        2.0 * (self.stiffness * self.mass).sqrt()
    }

    /// Damping ratio ζ = c / c_cr.
    pub fn damping_ratio(&self) -> f64 {
        self.damping / self.critical_damping()
    }

    /// Frequency response at forcing frequency `freq` \[Hz\].
    ///
    /// Returns `(magnitude, phase_rad)` of the displacement FRF
    /// H(ω) = 1 / (k − mω² + j·cω).
    pub fn frequency_response(&self, freq: f64) -> (f64, f64) {
        let omega = 2.0 * PI * freq;
        let re = self.stiffness - self.mass * omega * omega;
        let im = self.damping * omega;
        let mag = 1.0 / (re * re + im * im).sqrt();
        let phase = (-im).atan2(re);
        (mag, phase)
    }

    /// Static deflection under unit load x_st = 1/k \[m/N\].
    pub fn static_deflection(&self) -> f64 {
        1.0 / self.stiffness
    }

    /// Compute the displacement FRF over a frequency sweep.
    pub fn build_frf(&self, f_min: f64, f_max: f64, n_points: usize) -> FrequencyResponse {
        let df = if n_points > 1 {
            (f_max - f_min) / (n_points - 1) as f64
        } else {
            0.0
        };
        let mut freqs = Vec::with_capacity(n_points);
        let mut mags = Vec::with_capacity(n_points);
        let mut phases = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let f = f_min + i as f64 * df;
            let (m, p) = self.frequency_response(f);
            freqs.push(f);
            mags.push(m);
            phases.push(p);
        }
        FrequencyResponse::new(freqs, mags, phases)
    }
}

// ---------------------------------------------------------------------------
// NoiseSpl
// ---------------------------------------------------------------------------

/// Sound pressure level (SPL) spectrum in octave frequency bands.
///
/// SPL values are in decibels \[dB re 20 µPa\].
#[derive(Debug, Clone)]
pub struct NoiseSpl {
    /// SPL in each frequency band \[dB\].
    pub spl_db: Vec<f64>,
    /// Centre frequency of each band \[Hz\].
    pub frequency_bands: Vec<f64>,
}

impl NoiseSpl {
    /// Construct from parallel SPL and band-centre frequency vectors.
    pub fn new(spl_db: Vec<f64>, frequency_bands: Vec<f64>) -> Self {
        Self {
            spl_db,
            frequency_bands,
        }
    }

    /// Overall SPL by energy summation: L_total = 10 log₁₀(Σ 10^(Lᵢ/10)) \[dB\].
    pub fn overall_spl(&self) -> f64 {
        let sum: f64 = self.spl_db.iter().map(|l| 10.0_f64.powf(l / 10.0)).sum();
        if sum <= 0.0 {
            return f64::NEG_INFINITY;
        }
        10.0 * sum.log10()
    }

    /// A-weighted overall SPL.
    ///
    /// Applies the A-weighting correction to each band before energy summation.
    pub fn a_weighted_spl(&self) -> f64 {
        let sum: f64 = self
            .spl_db
            .iter()
            .zip(self.frequency_bands.iter())
            .map(|(l, f)| {
                let la = l + a_weighting_db(*f);
                10.0_f64.powf(la / 10.0)
            })
            .sum();
        if sum <= 0.0 {
            return f64::NEG_INFINITY;
        }
        10.0 * sum.log10()
    }
}

// ---------------------------------------------------------------------------
// VibrationIsolator
// ---------------------------------------------------------------------------

/// Passive vibration isolator (spring-damper in parallel).
///
/// Models the force transmissibility between source and receiver.
#[derive(Debug, Clone)]
pub struct VibrationIsolator {
    /// Isolator stiffness k \[N/m\].
    pub stiffness: f64,
    /// Isolator viscous damping c \[N·s/m\].
    pub damping: f64,
    /// Supported mass m \[kg\] (for natural-frequency computation).
    pub mass: f64,
}

impl VibrationIsolator {
    /// Construct a vibration isolator.
    pub fn new(stiffness: f64, damping: f64, mass: f64) -> Self {
        Self {
            stiffness,
            damping,
            mass,
        }
    }

    /// Undamped natural frequency of the isolation system \[Hz\].
    pub fn natural_frequency(&self) -> f64 {
        (self.stiffness / self.mass).sqrt() / (2.0 * PI)
    }

    /// Damping ratio ζ.
    pub fn damping_ratio(&self) -> f64 {
        let cc = 2.0 * (self.stiffness * self.mass).sqrt();
        self.damping / cc
    }

    /// Force transmissibility T(f) = |H_tr(jω)| (dimensionless, linear).
    ///
    /// T = √\[(1 + (2ζr)²) / ((1 − r²)² + (2ζr)²)\] where r = f / f_n.
    pub fn transmissibility(&self, freq: f64) -> f64 {
        let fn_ = self.natural_frequency();
        if fn_ < 1e-15 {
            return 1.0;
        }
        let r = freq / fn_;
        let zeta = self.damping_ratio();
        let two_zeta_r = 2.0 * zeta * r;
        let num = (1.0 + two_zeta_r * two_zeta_r).sqrt();
        let denom = ((1.0 - r * r).powi(2) + two_zeta_r * two_zeta_r).sqrt();
        if denom < 1e-15 {
            return f64::INFINITY;
        }
        num / denom
    }

    /// Insertion loss in dB at frequency `freq`.
    ///
    /// IL = −20 log₁₀(T) (positive IL means attenuation).
    pub fn insertion_loss_db(&self, freq: f64) -> f64 {
        let t = self.transmissibility(freq);
        if t <= 0.0 {
            return f64::INFINITY;
        }
        -20.0 * t.log10()
    }
}

// ---------------------------------------------------------------------------
// RoadRoughness
// ---------------------------------------------------------------------------

/// Road roughness model based on ISO 8608 power spectral density.
///
/// The PSD is modelled as a power law: Φ(n) = Φ₀ · (n/n₀)^{−w} where
/// n is spatial frequency \[1/m\], n₀ = 0.1 cyc/m, and w ≈ 2.
#[derive(Debug, Clone)]
pub struct RoadRoughness {
    /// International Roughness Index (IRI) \[m/km\].
    pub iri: f64,
    /// ISO 8608 road roughness coefficient Φ₀ \[(m³/cycle)\] at n₀ = 0.1 cyc/m.
    pub roughness_coefficient: f64,
    /// Spectral exponent w (typically 2.0).
    pub spectral_exponent: f64,
}

impl RoadRoughness {
    /// Construct from IRI.  The roughness coefficient is estimated as
    /// `Φ₀ ≈ 5e-6 * IRI²` (empirical fit).
    pub fn from_iri(iri: f64) -> Self {
        let phi0 = 5.0e-6 * iri * iri;
        Self {
            iri,
            roughness_coefficient: phi0,
            spectral_exponent: 2.0,
        }
    }

    /// Construct with explicit parameters.
    pub fn new(iri: f64, roughness_coefficient: f64, spectral_exponent: f64) -> Self {
        Self {
            iri,
            roughness_coefficient,
            spectral_exponent,
        }
    }

    /// One-sided PSD profile Φ(n) \[m²/(cyc/m)\] at spatial frequency `n` \[cyc/m\].
    ///
    /// Φ(n) = Φ₀ · (n / n₀)^{−w}   with n₀ = 0.1 cyc/m.
    pub fn psd_profile(&self, spatial_freq: f64) -> f64 {
        if spatial_freq <= 0.0 {
            return 0.0;
        }
        let n0 = 0.1_f64;
        self.roughness_coefficient * (spatial_freq / n0).powf(-self.spectral_exponent)
    }

    /// Approximate RMS amplitude \[m\] by integrating the PSD over a standard
    /// spatial frequency range \[0.01, 10\] cyc/m using a log-space trapezoidal rule.
    pub fn rms_amplitude(&self) -> f64 {
        let n_min = 0.01_f64;
        let n_max = 10.0_f64;
        let n_pts = 500_usize;
        let log_min = n_min.ln();
        let log_max = n_max.ln();
        let d_log = (log_max - log_min) / (n_pts - 1) as f64;

        let mut integral = 0.0_f64;
        let mut prev_n = n_min;
        let mut prev_phi = self.psd_profile(n_min);
        for i in 1..n_pts {
            let n = (log_min + i as f64 * d_log).exp();
            let phi = self.psd_profile(n);
            integral += 0.5 * (prev_phi + phi) * (n - prev_n);
            prev_n = n;
            prev_phi = phi;
        }
        integral.sqrt()
    }
}

// ---------------------------------------------------------------------------
// HarmonicOrder
// ---------------------------------------------------------------------------

/// Engine harmonic order analysis.
///
/// Firing orders of internal-combustion engines produce vibration at integer
/// and half-integer multiples of the fundamental RPM frequency.
#[derive(Debug, Clone)]
pub struct HarmonicOrder {
    /// Engine rotational speed \[RPM\].
    pub engine_rpm: f64,
    /// Harmonic order (e.g., 0.5, 1.0, 1.5, 2.0 …).
    pub order: f64,
}

impl HarmonicOrder {
    /// Construct a harmonic-order descriptor.
    pub fn new(engine_rpm: f64, order: f64) -> Self {
        Self { engine_rpm, order }
    }

    /// Frequency \[Hz\] of this order at the current RPM.
    ///
    /// f = order · RPM / 60
    pub fn frequency(&self) -> f64 {
        self.order * self.engine_rpm / 60.0
    }

    /// Compute all harmonic-order frequencies up to `max_order` at `rpm`.
    ///
    /// Returns a vector of frequencies in Hz for orders 0.5, 1.0, 1.5 … up to
    /// `max_order` in steps of 0.5 (half-orders).
    pub fn compute_orders(rpm: f64, max_order: f64) -> Vec<f64> {
        let mut freqs = Vec::new();
        let mut o = 0.5_f64;
        while o <= max_order + 1e-12 {
            freqs.push(o * rpm / 60.0);
            o += 0.5;
        }
        freqs
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Centre frequencies of `n_octaves` octave bands starting at 31.5 Hz.
///
/// The standard 1-octave band series: 31.5, 63, 125, 250, 500, 1000, 2000,
/// 4000, 8000, 16000 Hz …
pub fn octave_band_center_frequencies(n_octaves: usize) -> Vec<f64> {
    let f0 = 31.5_f64;
    (0..n_octaves)
        .map(|i| f0 * 2.0_f64.powi(i as i32))
        .collect()
}

/// A-weighting correction \[dB\] at frequency `freq` \[Hz\].
///
/// Implements the IEC 61672-1 formula.  Returns very negative values below
/// ~10 Hz where the weighting is undefined in practice.
pub fn a_weighting_db(freq: f64) -> f64 {
    if freq <= 0.0 {
        return -200.0;
    }
    let f2 = freq * freq;
    let f4 = f2 * f2;

    // A-weighting: numerator and denominator poles
    let num = 12194.0_f64.powi(2) * f4;
    let d1 = f2 + 20.6_f64.powi(2);
    let d2 = (f2 + 107.7_f64.powi(2)).sqrt() * (f2 + 737.9_f64.powi(2)).sqrt();
    let d3 = f2 + 12194.0_f64.powi(2);

    let ra = num / (d1 * d2 * d3);
    // Normalise so that RA(1000 Hz) = 1.0  (0 dB at 1 kHz)
    let ra_1k = {
        let f = 1000.0_f64;
        let f2k = f * f;
        let f4k = f2k * f2k;
        let n = 12194.0_f64.powi(2) * f4k;
        let dd1 = f2k + 20.6_f64.powi(2);
        let dd2 = (f2k + 107.7_f64.powi(2)).sqrt() * (f2k + 737.9_f64.powi(2)).sqrt();
        let dd3 = f2k + 12194.0_f64.powi(2);
        n / (dd1 * dd2 * dd3)
    };

    20.0 * (ra / ra_1k).log10()
}

/// Convert loudness in sone to phon.
///
/// Uses the Stevens' power law:
/// - For sone ≥ 1:  phon = 40 + 33.22 log₂(sone)
/// - For sone < 1:  phon = 40 + 10 log₂(sone) / log₂(2) (linear region)
///
/// Returns `0.0` for non-positive sone values.
pub fn sone_to_phon(sone: f64) -> f64 {
    if sone <= 0.0 {
        return 0.0;
    }
    if sone >= 1.0 {
        40.0 + 33.22 * sone.log2()
    } else {
        // Linear (approximately: phon ≈ 40 * sone for small values)
        40.0 * sone
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // ── SDOF natural frequency ────────────────────────────────────────────

    #[test]
    fn sdof_natural_frequency_known() {
        // k=4π²m → f_n = 1 Hz
        let m = 1.0_f64;
        let k = 4.0 * PI * PI;
        let sys = Sdof::new(m, 0.0, k);
        assert!((sys.natural_frequency() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn sdof_natural_frequency_10hz() {
        let k = (2.0 * PI * 10.0).powi(2);
        let sys = Sdof::new(1.0, 0.0, k);
        assert!((sys.natural_frequency() - 10.0).abs() < 1e-6);
    }

    // ── SDOF damping ratio ────────────────────────────────────────────────

    #[test]
    fn sdof_undamped_ratio_zero() {
        let sys = Sdof::new(1.0, 0.0, 100.0);
        assert!(sys.damping_ratio().abs() < 1e-12);
    }

    #[test]
    fn sdof_critical_damping_ratio_one() {
        let sys = Sdof::new(1.0, 0.0, 100.0);
        let cc = sys.critical_damping();
        let sys2 = Sdof::new(1.0, cc, 100.0);
        assert!((sys2.damping_ratio() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sdof_damping_ratio_05() {
        let m = 1.0_f64;
        let k = 100.0_f64;
        let cc = 2.0 * (k * m).sqrt();
        let sys = Sdof::new(m, 0.5 * cc, k);
        assert!((sys.damping_ratio() - 0.5).abs() < 1e-10);
    }

    // ── SDOF frequency response ───────────────────────────────────────────

    #[test]
    fn sdof_static_response_equals_1_over_k() {
        let k = 500.0;
        let sys = Sdof::new(1.0, 10.0, k);
        let (mag, _phase) = sys.frequency_response(0.0);
        assert!((mag - 1.0 / k).abs() < 1e-9, "mag={mag}");
    }

    #[test]
    fn sdof_resonance_magnitude_exceeds_static() {
        let sys = Sdof::new(1.0, 1.0, (2.0 * PI * 10.0).powi(2));
        let fn_ = sys.natural_frequency();
        let (mag_res, _) = sys.frequency_response(fn_);
        let (mag_0, _) = sys.frequency_response(0.001); // near-static
        assert!(mag_res > mag_0, "resonance mag should exceed static");
    }

    #[test]
    fn sdof_phase_at_resonance_near_minus_90() {
        // At fn, phase ≈ -π/2 for underdamped system
        let m = 1.0_f64;
        let k = (2.0 * PI * 5.0).powi(2);
        let cc = 2.0 * (k * m).sqrt();
        let sys = Sdof::new(m, 0.05 * cc, k);
        let fn_ = sys.natural_frequency();
        let (_, phase) = sys.frequency_response(fn_);
        assert!((phase + PI / 2.0).abs() < 0.1, "phase={phase}");
    }

    #[test]
    fn sdof_build_frf_correct_length() {
        let sys = Sdof::new(1.0, 5.0, 1000.0);
        let frf = sys.build_frf(1.0, 100.0, 50);
        assert_eq!(frf.frequencies.len(), 50);
        assert_eq!(frf.magnitudes.len(), 50);
        assert_eq!(frf.phases.len(), 50);
    }

    // ── FrequencyResponse ─────────────────────────────────────────────────

    #[test]
    fn frf_peak_frequency_correct() {
        let frf = FrequencyResponse::new(
            vec![10.0, 20.0, 30.0],
            vec![0.1, 0.5, 0.2],
            vec![0.0, 0.0, 0.0],
        );
        assert!((frf.peak_frequency() - 20.0).abs() < EPS);
    }

    #[test]
    fn frf_peak_magnitude_correct() {
        let frf = FrequencyResponse::new(vec![1.0, 2.0, 3.0], vec![0.3, 0.9, 0.6], vec![0.0; 3]);
        assert!((frf.peak_magnitude() - 0.9).abs() < EPS);
    }

    #[test]
    fn frf_empty_peak_frequency_zero() {
        let frf = FrequencyResponse::new(vec![], vec![], vec![]);
        assert_eq!(frf.peak_frequency(), 0.0);
    }

    #[test]
    fn frf_peak_magnitude_empty_zero() {
        let frf = FrequencyResponse::new(vec![], vec![], vec![]);
        assert_eq!(frf.peak_magnitude(), 0.0);
    }

    // ── A-weighting ───────────────────────────────────────────────────────

    #[test]
    fn a_weighting_1khz_is_zero_db() {
        let aw = a_weighting_db(1000.0);
        assert!(aw.abs() < 0.1, "A(1kHz)={aw}");
    }

    #[test]
    fn a_weighting_low_freq_negative() {
        let aw = a_weighting_db(100.0);
        assert!(aw < -10.0, "A(100Hz)={aw} should be < -10 dB");
    }

    #[test]
    fn a_weighting_zero_freq_very_negative() {
        let aw = a_weighting_db(0.0);
        assert!(aw < -100.0);
    }

    #[test]
    fn a_weighting_4khz_near_1db() {
        // At 4 kHz A-weighting ≈ +1 dB
        let aw = a_weighting_db(4000.0);
        assert!(aw > -5.0 && aw < 10.0, "A(4kHz)={aw}");
    }

    #[test]
    fn a_weighting_negative_freq_same_as_zero() {
        let aw = a_weighting_db(-50.0);
        assert!(aw < -100.0);
    }

    // ── NoiseSpl ──────────────────────────────────────────────────────────

    #[test]
    fn spl_overall_single_band() {
        let spl = NoiseSpl::new(vec![80.0], vec![1000.0]);
        assert!((spl.overall_spl() - 80.0).abs() < 1e-6);
    }

    #[test]
    fn spl_overall_two_equal_bands() {
        // Two bands at 80 dB each → overall = 80 + 3.01 ≈ 83.01 dB
        let spl = NoiseSpl::new(vec![80.0, 80.0], vec![500.0, 1000.0]);
        let overall = spl.overall_spl();
        assert!((overall - 83.01).abs() < 0.02, "overall={overall}");
    }

    #[test]
    fn spl_a_weighted_at_1khz_band_unchanged() {
        // Single band at 1 kHz: A-weighting ≈ 0 dB, so LA ≈ L
        let spl = NoiseSpl::new(vec![75.0], vec![1000.0]);
        let la = spl.a_weighted_spl();
        assert!((la - 75.0).abs() < 0.15, "LA={la}");
    }

    #[test]
    fn spl_a_weighted_less_than_unweighted_at_low_freq() {
        let spl = NoiseSpl::new(vec![80.0], vec![100.0]);
        let la = spl.a_weighted_spl();
        let l = spl.overall_spl();
        assert!(
            la < l,
            "A-weighted should be less at low freq: la={la}, l={l}"
        );
    }

    // ── VibrationIsolator / transmissibility ──────────────────────────────

    #[test]
    fn isolator_transmissibility_static_equals_one() {
        let iso = VibrationIsolator::new(10_000.0, 50.0, 10.0);
        let t = iso.transmissibility(0.0);
        assert!((t - 1.0).abs() < 1e-6, "T(0)={t}");
    }

    #[test]
    fn isolator_transmissibility_above_resonance_less_than_one() {
        let iso = VibrationIsolator::new(10_000.0, 10.0, 10.0);
        let fn_ = iso.natural_frequency();
        let t = iso.transmissibility(fn_ * 5.0);
        assert!(t < 1.0, "T={t} should be <1 above resonance");
    }

    #[test]
    fn isolator_insertion_loss_positive_above_resonance() {
        let iso = VibrationIsolator::new(10_000.0, 10.0, 10.0);
        let fn_ = iso.natural_frequency();
        let il = iso.insertion_loss_db(fn_ * 5.0);
        assert!(il > 0.0, "IL={il} should be positive above resonance");
    }

    #[test]
    fn isolator_insertion_loss_zero_at_dc() {
        let iso = VibrationIsolator::new(10_000.0, 50.0, 10.0);
        let il = iso.insertion_loss_db(0.0);
        assert!(il.abs() < 0.01, "IL(0)={il}");
    }

    #[test]
    fn isolator_natural_frequency_matches_sdof() {
        let k = 5_000.0;
        let m = 2.0;
        let iso = VibrationIsolator::new(k, 0.0, m);
        let expected = (k / m).sqrt() / (2.0 * PI);
        assert!((iso.natural_frequency() - expected).abs() < 1e-9);
    }

    // ── RoadRoughness ─────────────────────────────────────────────────────

    #[test]
    fn road_roughness_psd_positive() {
        let road = RoadRoughness::from_iri(2.0);
        let psd = road.psd_profile(0.1);
        assert!(psd > 0.0, "PSD={psd}");
    }

    #[test]
    fn road_roughness_psd_decreases_with_frequency() {
        let road = RoadRoughness::from_iri(2.0);
        let psd_low = road.psd_profile(0.05);
        let psd_high = road.psd_profile(1.0);
        assert!(psd_low > psd_high, "PSD should decrease with spatial freq");
    }

    #[test]
    fn road_roughness_zero_spatial_freq_returns_zero() {
        let road = RoadRoughness::from_iri(2.0);
        assert_eq!(road.psd_profile(0.0), 0.0);
    }

    #[test]
    fn road_roughness_rms_positive() {
        let road = RoadRoughness::from_iri(3.0);
        let rms = road.rms_amplitude();
        assert!(rms > 0.0, "RMS={rms}");
    }

    #[test]
    fn road_roughness_rms_scales_with_iri() {
        let rms_low = RoadRoughness::from_iri(1.0).rms_amplitude();
        let rms_high = RoadRoughness::from_iri(4.0).rms_amplitude();
        assert!(rms_high > rms_low, "rougher road should have higher RMS");
    }

    // ── HarmonicOrder ─────────────────────────────────────────────────────

    #[test]
    fn harmonic_order_frequency_first_order_1000rpm() {
        // Order 1 at 1000 RPM = 1000/60 Hz ≈ 16.67 Hz
        let ho = HarmonicOrder::new(1000.0, 1.0);
        let expected = 1000.0 / 60.0;
        assert!((ho.frequency() - expected).abs() < 1e-9);
    }

    #[test]
    fn harmonic_order_frequency_half_order() {
        let ho = HarmonicOrder::new(3000.0, 0.5);
        let expected = 0.5 * 3000.0 / 60.0;
        assert!((ho.frequency() - expected).abs() < 1e-9);
    }

    #[test]
    fn harmonic_order_compute_orders_count() {
        let orders = HarmonicOrder::compute_orders(1000.0, 4.0);
        // Half-orders 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 → 8 entries
        assert_eq!(orders.len(), 8);
    }

    #[test]
    fn harmonic_order_compute_orders_values() {
        let orders = HarmonicOrder::compute_orders(600.0, 2.0);
        // 0.5 * 600/60 = 5 Hz, 1.0 → 10, 1.5 → 15, 2.0 → 20
        let expected = [5.0, 10.0, 15.0, 20.0];
        assert_eq!(orders.len(), 4);
        for (a, b) in orders.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-9, "got {a}, expected {b}");
        }
    }

    #[test]
    fn harmonic_order_compute_single_order() {
        let orders = HarmonicOrder::compute_orders(6000.0, 0.5);
        assert_eq!(orders.len(), 1);
        assert!((orders[0] - 50.0).abs() < 1e-9);
    }

    // ── octave_band_center_frequencies ────────────────────────────────────

    #[test]
    fn octave_bands_count() {
        let bands = octave_band_center_frequencies(10);
        assert_eq!(bands.len(), 10);
    }

    #[test]
    fn octave_bands_start_at_31_5() {
        let bands = octave_band_center_frequencies(5);
        assert!((bands[0] - 31.5).abs() < 1e-9);
    }

    #[test]
    fn octave_bands_double_each_step() {
        let bands = octave_band_center_frequencies(5);
        for i in 1..bands.len() {
            assert!((bands[i] / bands[i - 1] - 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn octave_bands_empty_for_zero() {
        let bands = octave_band_center_frequencies(0);
        assert!(bands.is_empty());
    }

    // ── sone_to_phon ──────────────────────────────────────────────────────

    #[test]
    fn sone_to_phon_one_sone_is_40_phon() {
        assert!((sone_to_phon(1.0) - 40.0).abs() < EPS);
    }

    #[test]
    fn sone_to_phon_zero_sone_is_zero() {
        assert_eq!(sone_to_phon(0.0), 0.0);
    }

    #[test]
    fn sone_to_phon_negative_sone_is_zero() {
        assert_eq!(sone_to_phon(-1.0), 0.0);
    }

    #[test]
    fn sone_to_phon_two_sone_greater_than_one_sone() {
        assert!(sone_to_phon(2.0) > sone_to_phon(1.0));
    }

    #[test]
    fn sone_to_phon_monotone() {
        let values: Vec<f64> = (1..=10).map(|i| sone_to_phon(i as f64)).collect();
        for w in values.windows(2) {
            assert!(
                w[1] > w[0],
                "sone_to_phon not monotone: {} >= {}",
                w[1],
                w[0]
            );
        }
    }
}
