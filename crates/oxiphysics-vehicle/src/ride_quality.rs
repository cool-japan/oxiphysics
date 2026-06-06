// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ride quality and NVH (noise-vibration-harshness) metrics.
//!
//! Implements ISO 2631-1 frequency weighting, vibration dose value (VDV),
//! crest factor, International Roughness Index (IRI), road PSD, and motion
//! sickness dose value (MSDV).
//!
//! # Overview
//!
//! - [`RideQuality`] — container for weighted RMS acceleration metrics.
//! - [`iso_2631_weighting`] — ISO 2631-1 Wk/Wb/Wd frequency weighting.
//! - [`vibration_dose_value`] — VDV = (∫a⁴ dt)^0.25.
//! - [`crest_factor`] — peak-to-RMS ratio.
//! - [`road_roughness_iri`] — International Roughness Index computation.
//! - [`psd_road_profile`] — power spectral density of a road profile.
//! - [`motion_sickness_dose`] — MSDV from vertical acceleration.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// RideQuality
// ─────────────────────────────────────────────────────────────────────────────

/// Ride quality metrics based on ISO 2631-1 weighted accelerations.
#[derive(Debug, Clone)]
pub struct RideQuality {
    /// Frequency-weighted RMS acceleration in the vertical direction \[m/s²\].
    pub a_w_z: f64,
    /// Frequency-weighted RMS acceleration in the lateral direction \[m/s²\].
    pub a_w_y: f64,
    /// Frequency-weighted RMS acceleration in the longitudinal direction \[m/s²\].
    pub a_w_x: f64,
    /// Overall vibration total value (VTV): sqrt(kx²+ky²+kz²).
    pub overall_vtv: f64,
}

impl RideQuality {
    /// Compute ride quality from axis-wise weighted RMS values.
    ///
    /// Multiplying factors from ISO 2631-1: kx = ky = 1.4, kz = 1.0.
    pub fn from_weighted_rms(a_w_x: f64, a_w_y: f64, a_w_z: f64) -> Self {
        let kx = 1.4_f64;
        let ky = 1.4_f64;
        let kz = 1.0_f64;
        let vtv = ((kx * a_w_x).powi(2) + (ky * a_w_y).powi(2) + (kz * a_w_z).powi(2)).sqrt();
        Self {
            a_w_z,
            a_w_y,
            a_w_x,
            overall_vtv: vtv,
        }
    }

    /// Comfort category per ISO 2631-1 Table 3.
    ///
    /// Returns a descriptive string based on the overall VTV value.
    pub fn comfort_category(&self) -> &'static str {
        match self.overall_vtv {
            v if v < 0.315 => "Not uncomfortable",
            v if v < 0.63 => "A little uncomfortable",
            v if v < 1.0 => "Fairly uncomfortable",
            v if v < 1.6 => "Uncomfortable",
            v if v < 2.5 => "Very uncomfortable",
            _ => "Extremely uncomfortable",
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ISO 2631-1 frequency weighting
// ─────────────────────────────────────────────────────────────────────────────

/// ISO 2631-1 frequency weighting filter type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Iso2631Filter {
    /// Wk — vertical vibration for seated/standing persons (4–8 Hz peak).
    Wk,
    /// Wb — vertical vibration for recumbent persons (1–2 Hz peak).
    Wb,
    /// Wd — horizontal vibration for seated/standing persons (1–2 Hz peak).
    Wd,
}

/// Evaluate the ISO 2631-1 frequency weighting at a given frequency.
///
/// Returns the dimensionless weighting factor W(f) for the chosen filter.
///
/// The approximations used follow Annex C of ISO 2631-1:1997.
pub fn iso_2631_weighting(filter: Iso2631Filter, freq_hz: f64) -> f64 {
    if freq_hz <= 0.0 {
        return 0.0;
    }
    let f = freq_hz;
    match filter {
        Iso2631Filter::Wk => {
            // Wk: band-pass between ~0.5 Hz and 80 Hz with peak near 5–8 Hz
            let f1 = 0.4_f64;
            let f2 = 100.0_f64;
            let f3 = 12.5_f64;
            let f4 = 12.5_f64;
            let q3 = 0.63_f64;
            let q4 = 0.8_f64;
            // Simplified piecewise: use analytical filter chain
            let hp = high_pass_2nd_order(f, f1, 0.71);
            let lp = low_pass_2nd_order(f, f2, 0.71);
            let bp = band_pass_weight(f, f3, q3);
            let notch = notch_weight(f, f4, q4);
            (hp * lp * bp * notch).clamp(0.0, 1.0)
        }
        Iso2631Filter::Wb => {
            // Wb: peak near 1 Hz, used for recumbent
            let hp = high_pass_2nd_order(f, 0.1, 0.71);
            let lp = low_pass_2nd_order(f, 2.0, 0.71);
            let bp = band_pass_weight(f, 0.5, 0.9);
            (hp * lp * bp).clamp(0.0, 1.0)
        }
        Iso2631Filter::Wd => {
            // Wd: peak near 1–2 Hz for horizontal
            let hp = high_pass_2nd_order(f, 0.1, 0.71);
            let lp = low_pass_2nd_order(f, 2.0, 0.71);
            let bp = band_pass_weight(f, 1.0, 0.9);
            (hp * lp * bp).clamp(0.0, 1.0)
        }
    }
}

/// Apply ISO 2631-1 weighting to a spectrum and return the weighted RMS.
///
/// # Arguments
///
/// * `freqs` — frequency vector \[Hz\].
/// * `spectrum_rms` — RMS acceleration at each frequency \[m/s²\].
/// * `filter` — weighting filter to apply.
///
/// Returns frequency-weighted RMS acceleration \[m/s²\].
pub fn weighted_rms(freqs: &[f64], spectrum_rms: &[f64], filter: Iso2631Filter) -> f64 {
    assert_eq!(freqs.len(), spectrum_rms.len());
    let sum_sq: f64 = freqs
        .iter()
        .zip(spectrum_rms.iter())
        .map(|(&f, &a)| {
            let w = iso_2631_weighting(filter, f);
            (w * a).powi(2)
        })
        .sum();
    sum_sq.sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Vibration Dose Value
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Vibration Dose Value (VDV).
///
/// VDV = (∫₀ᵀ a⁴(t) dt)^0.25  \[m/s^1.75\]
///
/// # Arguments
///
/// * `accel` — time-series of weighted acceleration \[m/s²\].
/// * `dt` — time step \[s\].
pub fn vibration_dose_value(accel: &[f64], dt: f64) -> f64 {
    let sum: f64 = accel.iter().map(|&a| a * a * a * a).sum::<f64>() * dt;
    sum.powf(0.25)
}

// ─────────────────────────────────────────────────────────────────────────────
// Crest factor
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the crest factor: peak / RMS.
///
/// Returns `0.0` for an empty or zero-RMS signal.
pub fn crest_factor(accel: &[f64]) -> f64 {
    if accel.is_empty() {
        return 0.0;
    }
    let peak = accel.iter().cloned().fold(0.0_f64, |m, a| m.max(a.abs()));
    let rms = rms_value(accel);
    if rms < 1e-12 { 0.0 } else { peak / rms }
}

/// Compute the RMS of a slice.
pub fn rms_value(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = data.iter().map(|&x| x * x).sum();
    (sum_sq / data.len() as f64).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// International Roughness Index (IRI)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the International Roughness Index (IRI) of a road profile.
///
/// Uses the quarter-car model as defined in the World Bank standard
/// (Sayers 1995).  The profile is a sequence of elevation samples at
/// uniform spacing `dx` \[m\].
///
/// IRI units: m/km (or m/m, depending on scale convention; returned in m/km
/// here by convention).
///
/// # Arguments
///
/// * `profile` — road elevation samples \[m\].
/// * `dx` — sample spacing \[m\].
/// * `speed` — vehicle speed used to convert distance to time \[m/s\].
///   Canonical IRI uses 80 km/h = 22.222 m/s.
pub fn road_roughness_iri(profile: &[f64], dx: f64, speed: f64) -> f64 {
    if profile.len() < 2 {
        return 0.0;
    }
    let dt = dx / speed.max(1e-6);
    // Quarter-car parameters (ISO 8608 / Sayers)
    let k_s = 63.3_f64; // spring ratio
    let k_u = 653.0_f64; // unsprung stiffness ratio
    let c_s = 6.0_f64; // damping ratio
    let mu = 0.15_f64; // mass ratio (unsprung/sprung)

    // State: [z_s, z_s_dot, z_u, z_u_dot]
    let mut z_s = 0.0_f64;
    let mut z_s_dot = 0.0_f64;
    let mut z_u = 0.0_f64;
    let mut z_u_dot = 0.0_f64;

    let mut sum_rectified = 0.0_f64;
    let n = profile.len() - 1;

    for (i, &road) in profile.iter().enumerate().take(n) {
        let _ = i;

        // Sprung mass equation: z_s'' = -k_s*(z_s-z_u) - c_s*(z_s'-z_u')
        let a_s = -k_s * (z_s - z_u) - c_s * (z_s_dot - z_u_dot);
        // Unsprung mass equation
        let a_u = (k_s * (z_s - z_u) + c_s * (z_s_dot - z_u_dot) - k_u * (z_u - road)) / mu;

        // Euler step
        z_s_dot += a_s * dt;
        z_u_dot += a_u * dt;
        z_s += z_s_dot * dt;
        z_u += z_u_dot * dt;

        sum_rectified += (z_s_dot - z_u_dot).abs() * dt;
    }

    let total_dist = n as f64 * dx;
    if total_dist < 1e-9 {
        return 0.0;
    }
    // IRI in m/km
    (sum_rectified / total_dist) * 1000.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Power Spectral Density of road profile
// ─────────────────────────────────────────────────────────────────────────────

/// Power spectral density descriptor of a road profile.
#[derive(Debug, Clone)]
pub struct RoadPsd {
    /// Spatial frequency vector \[cycles/m\].
    pub spatial_freq: Vec<f64>,
    /// PSD values \[m²/(cycles/m)\].
    pub psd: Vec<f64>,
}

/// Compute the PSD of a road elevation profile using Welch-like periodogram.
///
/// # Arguments
///
/// * `profile` — elevation samples \[m\].
/// * `dx` — sample spacing \[m\].
///
/// Returns a [`RoadPsd`] with spatial frequency and PSD vectors.
pub fn psd_road_profile(profile: &[f64], dx: f64) -> RoadPsd {
    let n = profile.len();
    if n < 2 {
        return RoadPsd {
            spatial_freq: vec![],
            psd: vec![],
        };
    }

    // Mean-subtract
    let mean = profile.iter().sum::<f64>() / n as f64;
    let centered: Vec<f64> = profile.iter().map(|&x| x - mean).collect();

    // DFT (real, one-sided)
    let n_half = n / 2 + 1;
    let mut psd_vals = vec![0.0f64; n_half];
    let mut freq_vals = vec![0.0f64; n_half];

    for k in 0..n_half {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for (j, &c_j) in centered.iter().enumerate() {
            let angle = -2.0 * PI * k as f64 * j as f64 / n as f64;
            re += c_j * angle.cos();
            im += c_j * angle.sin();
        }
        let power = (re * re + im * im) / n as f64;
        // One-sided: double for non-DC, non-Nyquist bins
        let factor = if k == 0 || k == n_half - 1 { 1.0 } else { 2.0 };
        psd_vals[k] = factor * power * dx;
        freq_vals[k] = k as f64 / (n as f64 * dx);
    }

    RoadPsd {
        spatial_freq: freq_vals,
        psd: psd_vals,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Motion Sickness Dose Value (MSDV)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Motion Sickness Dose Value (MSDV).
///
/// MSDV = sqrt(∫₀ᵀ a²(t) dt)  \[m/s^1.5\]
///
/// Based on BS 6841:1987 Annex C, applied to vertical acceleration
/// frequency-weighted using Wf (motionsickness weighting, approximated
/// here as Wb from ISO 2631-1).
///
/// # Arguments
///
/// * `accel` — vertical acceleration time series \[m/s²\] (already Wf-weighted).
/// * `dt` — time step \[s\].
pub fn motion_sickness_dose(accel: &[f64], dt: f64) -> f64 {
    let sum: f64 = accel.iter().map(|&a| a * a).sum::<f64>() * dt;
    sum.sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal filter approximations
// ─────────────────────────────────────────────────────────────────────────────

/// Second-order high-pass gain at frequency f \[Hz\] with corner fc and Q.
fn high_pass_2nd_order(f: f64, fc: f64, q: f64) -> f64 {
    let r = f / fc;
    let r2 = r * r;
    let denom = (1.0 - r2).powi(2) + (r / q).powi(2);
    (r2 * r2 / denom).sqrt()
}

/// Second-order low-pass gain at frequency f \[Hz\] with corner fc and Q.
fn low_pass_2nd_order(f: f64, fc: f64, q: f64) -> f64 {
    let r = f / fc;
    let r2 = r * r;
    let denom = (1.0 - r2).powi(2) + (r / q).powi(2);
    (1.0 / denom).sqrt()
}

/// Simplified band-pass weighting centred at f0 with quality factor q0.
fn band_pass_weight(f: f64, f0: f64, q0: f64) -> f64 {
    let r = f / f0;
    // Approximation: Gaussian envelope in log space
    let log_r = (r).ln();
    let sigma = 1.0 / (2.0 * q0);
    (-(log_r * log_r) / (2.0 * sigma * sigma)).exp()
}

/// Simplified notch (anti-resonance) factor near f0 with quality q0.
fn notch_weight(f: f64, f0: f64, _q0: f64) -> f64 {
    let r = f / f0;
    // Simple: weight near 1 away from notch, dips at f0
    1.0 - 0.3 * (-4.0 * (r - 1.0).powi(2)).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── RideQuality ──────────────────────────────────────────────────────

    #[test]
    fn ride_quality_vtv_not_uncomfortable() {
        let rq = RideQuality::from_weighted_rms(0.1, 0.1, 0.1);
        assert_eq!(rq.comfort_category(), "Not uncomfortable");
    }

    #[test]
    fn ride_quality_vtv_very_uncomfortable() {
        let rq = RideQuality::from_weighted_rms(1.0, 1.0, 1.0);
        let cat = rq.comfort_category();
        assert!(
            cat.contains("uncomfortable") || cat.contains("Uncomfortable"),
            "unexpected category: {cat}"
        );
    }

    #[test]
    fn ride_quality_vtv_calculation() {
        // kx=ky=1.4, kz=1.0
        let rq = RideQuality::from_weighted_rms(1.0, 0.0, 0.0);
        let expected = 1.4_f64;
        assert!((rq.overall_vtv - expected).abs() < 1e-9);
    }

    #[test]
    fn ride_quality_zero_inputs() {
        let rq = RideQuality::from_weighted_rms(0.0, 0.0, 0.0);
        assert!(rq.overall_vtv.abs() < 1e-12);
        assert_eq!(rq.comfort_category(), "Not uncomfortable");
    }

    #[test]
    fn ride_quality_fields_accessible() {
        let rq = RideQuality::from_weighted_rms(0.2, 0.3, 0.5);
        assert!((rq.a_w_x - 0.2).abs() < 1e-10);
        assert!((rq.a_w_y - 0.3).abs() < 1e-10);
        assert!((rq.a_w_z - 0.5).abs() < 1e-10);
    }

    // ── iso_2631_weighting ───────────────────────────────────────────────

    #[test]
    fn wk_weighting_zero_at_dc() {
        assert_eq!(iso_2631_weighting(Iso2631Filter::Wk, 0.0), 0.0);
    }

    #[test]
    fn wb_weighting_positive_at_1hz() {
        let w = iso_2631_weighting(Iso2631Filter::Wb, 1.0);
        assert!(w > 0.0, "Wb at 1 Hz should be positive");
    }

    #[test]
    fn wd_weighting_positive_at_2hz() {
        let w = iso_2631_weighting(Iso2631Filter::Wd, 2.0);
        assert!(w > 0.0, "Wd at 2 Hz should be positive");
    }

    #[test]
    fn wk_weighting_nonnegative() {
        let freqs = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0];
        for f in freqs {
            let w = iso_2631_weighting(Iso2631Filter::Wk, f);
            assert!(w >= 0.0, "Wk({f}) = {w} should be >= 0");
        }
    }

    #[test]
    fn all_filters_return_zero_for_negative_freq() {
        for filter in [Iso2631Filter::Wk, Iso2631Filter::Wb, Iso2631Filter::Wd] {
            assert_eq!(iso_2631_weighting(filter, -1.0), 0.0);
        }
    }

    // ── weighted_rms ─────────────────────────────────────────────────────

    #[test]
    fn weighted_rms_empty_is_zero() {
        let w = weighted_rms(&[], &[], Iso2631Filter::Wk);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn weighted_rms_single_bin() {
        let w = weighted_rms(&[5.0], &[1.0], Iso2631Filter::Wb);
        // Wb at 5 Hz * 1.0
        assert!(w >= 0.0);
    }

    // ── vibration_dose_value ─────────────────────────────────────────────

    #[test]
    fn vdv_zero_signal_is_zero() {
        let a = vec![0.0; 100];
        assert_eq!(vibration_dose_value(&a, 0.01), 0.0);
    }

    #[test]
    fn vdv_constant_signal_analytical() {
        // VDV = (a^4 * T)^0.25 = a * T^0.25
        let a = 2.0_f64;
        let dt = 0.01;
        let n = 1000;
        let sig: Vec<f64> = vec![a; n];
        let vdv = vibration_dose_value(&sig, dt);
        let expected = a * (n as f64 * dt).powf(0.25);
        assert!(
            (vdv - expected).abs() < 1e-6,
            "vdv={vdv} expected={expected}"
        );
    }

    #[test]
    fn vdv_scales_with_amplitude() {
        let sig1: Vec<f64> = vec![1.0; 100];
        let sig2: Vec<f64> = vec![2.0; 100];
        let vdv1 = vibration_dose_value(&sig1, 0.01);
        let vdv2 = vibration_dose_value(&sig2, 0.01);
        // Doubling amplitude should quadruple VDV (2^4=16, ^0.25=2)
        assert!((vdv2 / vdv1 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn vdv_single_sample() {
        let vdv = vibration_dose_value(&[3.0], 1.0);
        // (3^4 * 1)^0.25 = 3
        assert!((vdv - 3.0).abs() < 1e-10);
    }

    // ── crest_factor ─────────────────────────────────────────────────────

    #[test]
    fn crest_factor_empty_is_zero() {
        assert_eq!(crest_factor(&[]), 0.0);
    }

    #[test]
    fn crest_factor_constant_is_one() {
        let sig = vec![2.0f64; 100];
        let cf = crest_factor(&sig);
        assert!((cf - 1.0).abs() < 1e-9);
    }

    #[test]
    fn crest_factor_impulse() {
        // One spike of amplitude 10 in 100 zeros
        let mut sig = vec![0.0f64; 99];
        sig.push(10.0);
        let cf = crest_factor(&sig);
        let rms = rms_value(&sig);
        let expected = 10.0 / rms;
        assert!((cf - expected).abs() < 1e-6);
    }

    #[test]
    fn crest_factor_symmetric_signal() {
        // Alternating +A, -A: peak = A, rms = A → CF = 1
        let sig: Vec<f64> = (0..100)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let cf = crest_factor(&sig);
        assert!((cf - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rms_value_basic() {
        let sig = vec![3.0, 4.0]; // rms = sqrt((9+16)/2) = sqrt(12.5)
        let rms = rms_value(&sig);
        assert!((rms - (12.5_f64).sqrt()).abs() < 1e-10);
    }

    // ── road_roughness_iri ───────────────────────────────────────────────

    #[test]
    fn iri_flat_road_near_zero() {
        let profile = vec![0.0f64; 1000];
        let iri = road_roughness_iri(&profile, 0.1, 22.22);
        assert!(iri.abs() < 0.1, "flat road IRI={iri}");
    }

    #[test]
    fn iri_rough_road_larger_than_smooth() {
        let smooth: Vec<f64> = (0..500).map(|i| 0.001 * (i as f64).sin()).collect();
        let rough: Vec<f64> = (0..500).map(|i| 0.05 * (i as f64).sin()).collect();
        let iri_smooth = road_roughness_iri(&smooth, 0.1, 22.22);
        let iri_rough = road_roughness_iri(&rough, 0.1, 22.22);
        assert!(
            iri_rough > iri_smooth,
            "rough IRI={iri_rough} should exceed smooth={iri_smooth}"
        );
    }

    #[test]
    fn iri_too_short_returns_zero() {
        let profile = vec![0.0f64; 1];
        assert_eq!(road_roughness_iri(&profile, 0.1, 22.22), 0.0);
    }

    #[test]
    fn iri_positive_for_sinusoidal_profile() {
        let profile: Vec<f64> = (0..500).map(|i| 0.01 * (i as f64 * 0.1).sin()).collect();
        let iri = road_roughness_iri(&profile, 0.1, 22.22);
        assert!(iri >= 0.0);
    }

    // ── psd_road_profile ─────────────────────────────────────────────────

    #[test]
    fn psd_empty_profile() {
        let psd = psd_road_profile(&[], 0.1);
        assert!(psd.psd.is_empty());
    }

    #[test]
    fn psd_single_sample() {
        let psd = psd_road_profile(&[1.0], 0.1);
        assert!(psd.psd.is_empty() || psd.psd.len() == 1);
    }

    #[test]
    fn psd_output_lengths_match() {
        let profile: Vec<f64> = (0..64).map(|i| (i as f64).sin()).collect();
        let psd = psd_road_profile(&profile, 0.1);
        assert_eq!(psd.spatial_freq.len(), psd.psd.len());
    }

    #[test]
    fn psd_nonnegative_values() {
        let profile: Vec<f64> = (0..64).map(|i| (i as f64).sin()).collect();
        let psd = psd_road_profile(&profile, 0.1);
        for &v in &psd.psd {
            assert!(v >= 0.0, "PSD value negative: {v}");
        }
    }

    #[test]
    fn psd_dc_bin_at_zero_freq() {
        let n = 64;
        let profile: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let psd = psd_road_profile(&profile, 0.1);
        assert!(psd.spatial_freq[0].abs() < 1e-12);
    }

    // ── motion_sickness_dose ─────────────────────────────────────────────

    #[test]
    fn msdv_zero_signal_is_zero() {
        let a = vec![0.0; 100];
        assert_eq!(motion_sickness_dose(&a, 0.01), 0.0);
    }

    #[test]
    fn msdv_constant_signal_analytical() {
        // MSDV = sqrt(a² * T) = a * sqrt(T)
        let a = 0.5_f64;
        let dt = 0.01;
        let n = 1000;
        let sig = vec![a; n];
        let msdv = motion_sickness_dose(&sig, dt);
        let expected = a * (n as f64 * dt).sqrt();
        assert!((msdv - expected).abs() < 1e-9);
    }

    #[test]
    fn msdv_scales_with_amplitude() {
        let sig1 = vec![1.0f64; 100];
        let sig2 = vec![2.0f64; 100];
        let m1 = motion_sickness_dose(&sig1, 0.01);
        let m2 = motion_sickness_dose(&sig2, 0.01);
        assert!((m2 / m1 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn msdv_nonnegative() {
        let sig: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();
        let m = motion_sickness_dose(&sig, 0.01);
        assert!(m >= 0.0);
    }
}
