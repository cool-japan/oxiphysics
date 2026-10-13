// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ride quality and NVH (noise-vibration-harshness) metrics.
//!
//! Implements the ISO 2631-1 Wk/Wd and BS 6841 / ISO 2631-4 Wb frequency
//! weightings, vibration dose value (VDV), crest factor, International
//! Roughness Index (IRI), road PSD, and motion sickness dose value (MSDV).
//!
//! # Overview
//!
//! - [`RideQuality`] — container for weighted RMS acceleration metrics.
//! - [`iso_2631_weighting`] — Wk/Wb/Wd whole-body frequency weighting.
//! - [`vibration_dose_value`] — VDV = (∫a⁴ dt)^0.25.
//! - [`crest_factor`] — peak-to-RMS ratio.
//! - [`road_roughness_iri`] — International Roughness Index computation.
//! - [`psd_road_profile`] — power spectral density of a road profile.
//! - [`motion_sickness_dose`] — MSDV from vertical acceleration.

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

/// Whole-body vibration frequency weighting filter type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Iso2631Filter {
    /// Wk — vertical (z-axis) vibration of seated, standing and recumbent
    /// persons (ISO 2631-1).  Peaks at about 6.3 Hz with a gain of 1.054.
    Wk,
    /// Wb — vertical whole-body weighting of BS 6841:1987, also used by
    /// ISO 2631-4 for ride comfort in fixed-guideway (railway) transport.
    /// Roughly unit gain on a plateau between 5 and 8 Hz.
    Wb,
    /// Wd — horizontal (x- and y-axis) whole-body vibration (ISO 2631-1).
    /// Roughly unit gain between 0.8 and 1.6 Hz, falling off above 2 Hz.
    Wd,
}

/// Butterworth quality factor `1/√2` of both band-limiting filters.
const BAND_LIMITING_Q: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// Upward-step stage of a weighting (`f5`, `Q5`, `f6`, `Q6` of ISO 2631-1).
#[derive(Debug, Clone, Copy)]
struct UpwardStep {
    f5: f64,
    q5: f64,
    f6: f64,
    q6: f64,
}

/// Corner-frequency parameters of a weighting (ISO 2631-1:1997, Annex A).
#[derive(Debug, Clone, Copy)]
struct WeightingParameters {
    /// Band-limiting high-pass corner \[Hz\].
    f1: f64,
    /// Band-limiting low-pass corner \[Hz\].
    f2: f64,
    /// Zero of the acceleration-velocity transition \[Hz\].
    f3: f64,
    /// Pole frequency of the acceleration-velocity transition \[Hz\].
    f4: f64,
    /// Quality factor of the transition pole pair.
    q4: f64,
    /// Upward step; `None` when the weighting has no such stage (Wd).
    step: Option<UpwardStep>,
}

impl Iso2631Filter {
    /// Corner frequencies and quality factors of this weighting.
    fn parameters(self) -> WeightingParameters {
        match self {
            Iso2631Filter::Wk => WeightingParameters {
                f1: 0.4,
                f2: 100.0,
                f3: 12.5,
                f4: 12.5,
                q4: 0.63,
                step: Some(UpwardStep {
                    f5: 2.37,
                    q5: 0.91,
                    f6: 3.35,
                    q6: 0.91,
                }),
            },
            Iso2631Filter::Wd => WeightingParameters {
                f1: 0.4,
                f2: 100.0,
                f3: 2.0,
                f4: 2.0,
                q4: 0.63,
                step: None,
            },
            Iso2631Filter::Wb => WeightingParameters {
                f1: 0.4,
                f2: 100.0,
                f3: 16.0,
                f4: 16.0,
                q4: 0.55,
                step: Some(UpwardStep {
                    f5: 2.5,
                    q5: 0.9,
                    f6: 4.0,
                    q6: 0.95,
                }),
            },
        }
    }
}

/// Evaluate the whole-body vibration frequency weighting at a given frequency.
///
/// Returns the dimensionless weighting factor `W(f) = |H(j·2πf)|` for the
/// chosen filter, where the transfer function is the product of four stages
/// (ISO 2631-1:1997, Annex A; `s` is the Laplace variable, `ωᵢ = 2π·fᵢ`):
///
/// * band-limiting high-pass `H_h = s² / (s² + (ω1/Q1)·s + ω1²)`
/// * band-limiting low-pass `H_l = ω2² / (s² + (ω2/Q2)·s + ω2²)`, with
///   `Q1 = Q2 = 1/√2`
/// * acceleration-velocity transition
///   `H_t = (1 + s/ω3) / (1 + s/(Q4·ω4) + s²/ω4²)`
/// * upward step
///   `H_s = (ω5/ω6)² · (1 + s/(Q5·ω5) + s²/ω5²) / (1 + s/(Q6·ω6) + s²/ω6²)`
///   (absent for Wd, `H_s = 1`)
///
/// | filter | f1 | f2 | f3 = f4 | Q4 | f5 | Q5 | f6 | Q6 |
/// |--------|----|----|---------|----|----|----|----|----|
/// | Wk | 0.4 | 100 | 12.5 | 0.63 | 2.37 | 0.91 | 3.35 | 0.91 |
/// | Wd | 0.4 | 100 | 2.0 | 0.63 | – | – | – | – |
/// | Wb | 0.4 | 100 | 16 | 0.55 | 2.5 | 0.9 | 4.0 | 0.95 |
///
/// The result is *not* limited to 1: Wk reaches 1.054 at 6.3 Hz, and
/// reproduces the tabulated ISO 2631-1 values (Wk(1 Hz) = 0.482,
/// Wd(1 Hz) = 1.011).  Wb uses the parameters of BS 6841:1987 with the exact
/// step gain `(f5/f6)² = 0.3906`; BS 6841 quotes this gain rounded as
/// `K = 0.4`, which would make Wb about 2.4 % larger.
///
/// Returns `0.0` for non-positive and for non-finite (infinite or NaN)
/// frequencies.
pub fn iso_2631_weighting(filter: Iso2631Filter, freq_hz: f64) -> f64 {
    if freq_hz <= 0.0 || !freq_hz.is_finite() {
        return 0.0;
    }
    let f = freq_hz;
    let p = filter.parameters();
    let high_pass = high_pass_magnitude(f, p.f1, BAND_LIMITING_Q);
    let low_pass = low_pass_magnitude(f, p.f2, BAND_LIMITING_Q);
    let transition = transition_magnitude(f, p.f3, p.f4, p.q4);
    let step = match p.step {
        Some(st) => upward_step_magnitude(f, st),
        None => 1.0,
    };
    high_pass * low_pass * transition * step
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

/// Compute the one-sided PSD of a road elevation profile as a single
/// (unwindowed, unaveraged) periodogram.
///
/// The profile is mean-subtracted and transformed with a real-input FFT
/// ([`oxifft::rfft`]), so the cost is `O(n log n)` for any length `n`.  With
/// `X_k` the DFT of the centred profile, the one-sided density is
///
/// `S(n_k) = c_k · |X_k|² · dx / N`,   `n_k = k / (N · dx)`
///
/// where `c_k = 1` for the DC bin and (for even `N`) the Nyquist bin, and
/// `c_k = 2` for every other bin because the negative-frequency half of the
/// spectrum is folded onto the positive bins.  The normalisation satisfies
/// Parseval's theorem: `Σ S(n_k) · Δn` equals the mean square of the centred
/// profile, with `Δn = 1 / (N · dx)`.
///
/// # Arguments
///
/// * `profile` — elevation samples \[m\].
/// * `dx` — sample spacing \[m\].
///
/// Returns a [`RoadPsd`] with spatial frequency and PSD vectors (`N / 2 + 1`
/// bins), or empty vectors when fewer than two samples are supplied.
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

    // Real-input FFT: the n/2 + 1 non-redundant bins.
    let spectrum = oxifft::rfft(&centered);
    let n_f = n as f64;
    let n_is_even = n.is_multiple_of(2);

    let mut psd_vals = Vec::with_capacity(spectrum.len());
    let mut freq_vals = Vec::with_capacity(spectrum.len());
    for (k, bin) in spectrum.iter().enumerate() {
        let power = (bin.re * bin.re + bin.im * bin.im) * dx / n_f;
        // One-sided: double every bin that has a mirror image.  Only DC and,
        // for even n, the Nyquist bin are their own mirror.
        let unpaired = k == 0 || (n_is_even && k == n / 2);
        let factor = if unpaired { 1.0 } else { 2.0 };
        psd_vals.push(factor * power);
        freq_vals.push(k as f64 / (n_f * dx));
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
/// Based on BS 6841:1987, applied to vertical acceleration that has already
/// been frequency-weighted with the motion-sickness weighting Wf.  This module
/// does not implement Wf (it is not one of [`Iso2631Filter`]'s variants), so
/// the caller must supply the weighted signal.
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
// Internal transfer-function stages
// ─────────────────────────────────────────────────────────────────────────────

/// Magnitude of the normalised second-order polynomial
/// `1 + s/(Q·ωn) + (s/ωn)²` at `s = j·2πf`, i.e. `|1 − r² + j·r/Q|` with
/// `r = f / f_n`.
fn quadratic_section_magnitude(f: f64, f_n: f64, q: f64) -> f64 {
    let r = f / f_n;
    (1.0 - r * r).hypot(r / q)
}

/// `|s² / (s² + (ω1/Q)·s + ω1²)|` at frequency `f` \[Hz\] (second-order
/// high-pass with corner `fc` and quality factor `q`).
fn high_pass_magnitude(f: f64, fc: f64, q: f64) -> f64 {
    let r = f / fc;
    r * r / quadratic_section_magnitude(f, fc, q)
}

/// `|ω2² / (s² + (ω2/Q)·s + ω2²)|` at frequency `f` \[Hz\] (second-order
/// low-pass with corner `fc` and quality factor `q`).
fn low_pass_magnitude(f: f64, fc: f64, q: f64) -> f64 {
    1.0 / quadratic_section_magnitude(f, fc, q)
}

/// Acceleration-velocity transition `|(1 + s/ω3) / (1 + s/(Q4·ω4) + s²/ω4²)|`.
fn transition_magnitude(f: f64, f3: f64, f4: f64, q4: f64) -> f64 {
    let lead = 1.0_f64.hypot(f / f3);
    lead / quadratic_section_magnitude(f, f4, q4)
}

/// Upward step
/// `(ω5/ω6)² · |1 + s/(Q5·ω5) + s²/ω5²| / |1 + s/(Q6·ω6) + s²/ω6²|`.
fn upward_step_magnitude(f: f64, step: UpwardStep) -> f64 {
    let ratio = step.f5 / step.f6;
    ratio * ratio * quadratic_section_magnitude(f, step.f5, step.q5)
        / quadratic_section_magnitude(f, step.f6, step.q6)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── independent complex-arithmetic reference for the weightings ────────

    type Cx = (f64, f64);

    fn cx_add(a: Cx, b: Cx) -> Cx {
        (a.0 + b.0, a.1 + b.1)
    }
    fn cx_mul(a: Cx, b: Cx) -> Cx {
        (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
    }
    fn cx_div(a: Cx, b: Cx) -> Cx {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    }
    fn cx_scale(a: Cx, k: f64) -> Cx {
        (a.0 * k, a.1 * k)
    }
    fn cx_abs(a: Cx) -> f64 {
        a.0.hypot(a.1)
    }

    /// Parameters `(f1, f2, f3, f4, Q4, Option<(f5, Q5, f6, Q6)>)` as printed in
    /// ISO 2631-1:1997 Annex A (Wk, Wd) and BS 6841:1987 (Wb).
    type Table = (f64, f64, f64, f64, f64, Option<(f64, f64, f64, f64)>);
    const WK_TABLE: Table = (0.4, 100.0, 12.5, 12.5, 0.63, Some((2.37, 0.91, 3.35, 0.91)));
    const WD_TABLE: Table = (0.4, 100.0, 2.0, 2.0, 0.63, None);
    const WB_TABLE: Table = (0.4, 100.0, 16.0, 16.0, 0.55, Some((2.5, 0.9, 4.0, 0.95)));

    /// `|H(j·2πf)|` evaluated by complex arithmetic directly from the
    /// Laplace-domain definition of the four weighting stages.
    fn reference_weighting(f: f64, table: Table) -> f64 {
        let (f1, f2, f3, f4, q4, step) = table;
        let w = |fc: f64| 2.0 * PI * fc;
        let q = std::f64::consts::FRAC_1_SQRT_2;
        let s: Cx = (0.0, w(f));
        let s2 = cx_mul(s, s);
        let one: Cx = (1.0, 0.0);
        // H_h = s² / (s² + (ω1/Q1)·s + ω1²)
        let den_h = cx_add(cx_add(s2, cx_scale(s, w(f1) / q)), (w(f1).powi(2), 0.0));
        let h_h = cx_div(s2, den_h);
        // H_l = ω2² / (s² + (ω2/Q2)·s + ω2²)
        let den_l = cx_add(cx_add(s2, cx_scale(s, w(f2) / q)), (w(f2).powi(2), 0.0));
        let h_l = cx_div((w(f2).powi(2), 0.0), den_l);
        // H_t = (1 + s/ω3) / (1 + s/(Q4·ω4) + s²/ω4²)
        let num_t = cx_add(one, cx_scale(s, 1.0 / w(f3)));
        let den_t = cx_add(
            cx_add(one, cx_scale(s, 1.0 / (q4 * w(f4)))),
            cx_scale(s2, 1.0 / w(f4).powi(2)),
        );
        let h_t = cx_div(num_t, den_t);
        // H_s = (ω5/ω6)² (1 + s/(Q5 ω5) + s²/ω5²) / (1 + s/(Q6 ω6) + s²/ω6²)
        let h_s = match step {
            Some((f5, q5, f6, q6)) => {
                let num = cx_add(
                    cx_add(one, cx_scale(s, 1.0 / (q5 * w(f5)))),
                    cx_scale(s2, 1.0 / w(f5).powi(2)),
                );
                let den = cx_add(
                    cx_add(one, cx_scale(s, 1.0 / (q6 * w(f6)))),
                    cx_scale(s2, 1.0 / w(f6).powi(2)),
                );
                cx_scale(cx_div(num, den), (w(f5) / w(f6)).powi(2))
            }
            None => one,
        };
        cx_abs(cx_mul(cx_mul(h_h, h_l), cx_mul(h_t, h_s)))
    }

    /// Log-spaced frequency grid from `lo` to `hi` (Hz) with `n` points.
    fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| lo * (hi / lo).powf(i as f64 / (n - 1) as f64))
            .collect()
    }

    /// Frequency at which `filter` is largest on a fine grid, and that maximum.
    fn weighting_peak(filter: Iso2631Filter) -> (f64, f64) {
        let mut best = (0.0, 0.0);
        for f in log_grid(0.1, 100.0, 4000) {
            let w = iso_2631_weighting(filter, f);
            if w > best.1 {
                best = (f, w);
            }
        }
        best
    }

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

    #[test]
    fn all_filters_return_zero_for_infinite_freq() {
        for filter in [Iso2631Filter::Wk, Iso2631Filter::Wb, Iso2631Filter::Wd] {
            assert_eq!(iso_2631_weighting(filter, f64::INFINITY), 0.0);
        }
    }

    #[test]
    fn all_filters_return_zero_for_nan_freq() {
        for filter in [Iso2631Filter::Wk, Iso2631Filter::Wb, Iso2631Filter::Wd] {
            assert_eq!(iso_2631_weighting(filter, f64::NAN), 0.0);
        }
    }

    #[test]
    fn wk_matches_iso_2631_1_tabulated_values() {
        // ISO 2631-1:1997 Table 3 (one-third-octave band centre frequencies).
        let table = [
            (0.5, 0.418),
            (1.0, 0.482),
            (2.0, 0.531),
            (4.0, 0.967),
            (6.3, 1.054),
            (8.0, 1.036),
            (16.0, 0.768),
            (100.0, 0.0887),
        ];
        for (f, expected) in table {
            let w = iso_2631_weighting(Iso2631Filter::Wk, f);
            assert!(
                (w - expected).abs() <= 0.005 * expected,
                "Wk({f} Hz) = {w}, ISO table {expected}"
            );
        }
    }

    #[test]
    fn wd_matches_iso_2631_1_tabulated_values() {
        let table = [(0.5, 0.853), (1.0, 1.011), (2.0, 0.890), (4.0, 0.512)];
        for (f, expected) in table {
            let w = iso_2631_weighting(Iso2631Filter::Wd, f);
            assert!(
                (w - expected).abs() <= 0.005 * expected,
                "Wd({f} Hz) = {w}, ISO table {expected}"
            );
        }
    }

    #[test]
    fn wk_one_hertz_is_not_the_old_gaussian_value() {
        // Regression: the invented Gaussian/notch curve gave Wk(1 Hz) ≈ 0.006.
        let w = iso_2631_weighting(Iso2631Filter::Wk, 1.0);
        assert!((w - 0.482).abs() < 0.005 * 0.482, "Wk(1 Hz) = {w}");
        assert!(w > 0.4, "Wk(1 Hz) must not be near zero: {w}");
    }

    #[test]
    fn weightings_are_not_clamped_to_unity() {
        let w = iso_2631_weighting(Iso2631Filter::Wk, 6.3);
        assert!(w > 1.05, "Wk peaks at 1.054, got {w}");
        let wd = iso_2631_weighting(Iso2631Filter::Wd, 1.0);
        assert!(wd > 1.0, "Wd(1 Hz) = 1.011 exceeds unity, got {wd}");
    }

    #[test]
    fn wk_peak_is_near_six_hertz() {
        // Regression: the invented curve peaked at 12.5 Hz.
        let (f_peak, w_peak) = weighting_peak(Iso2631Filter::Wk);
        assert!(
            (5.0..=8.0).contains(&f_peak),
            "Wk peak at {f_peak} Hz, expected near 6.3 Hz"
        );
        assert!((w_peak - 1.054).abs() < 0.01, "peak gain {w_peak}");
    }

    #[test]
    fn wd_is_flat_in_the_one_to_two_hertz_band() {
        let (f_peak, w_peak) = weighting_peak(Iso2631Filter::Wd);
        assert!((0.8..=1.6).contains(&f_peak), "Wd peak at {f_peak} Hz");
        assert!((1.0..1.02).contains(&w_peak), "Wd peak gain {w_peak}");
        for f in [0.8, 1.0, 1.25, 1.6] {
            let w = iso_2631_weighting(Iso2631Filter::Wd, f);
            assert!((0.96..=1.02).contains(&w), "Wd({f} Hz) = {w}");
        }
        // Horizontal weighting falls off well before the vertical one.
        assert!(
            iso_2631_weighting(Iso2631Filter::Wd, 8.0) < 0.3,
            "Wd must be small at 8 Hz"
        );
    }

    #[test]
    fn wb_has_a_vertical_plateau_between_five_and_eight_hertz() {
        let (f_peak, w_peak) = weighting_peak(Iso2631Filter::Wb);
        assert!(
            (5.0..=8.0).contains(&f_peak),
            "Wb peak at {f_peak} Hz, expected on the 5-8 Hz plateau"
        );
        assert!(w_peak > 1.0 && w_peak < 1.05, "Wb peak gain {w_peak}");
        // Rising edge: strongly attenuated at 1 Hz, monotone up to the plateau.
        let mut previous = iso_2631_weighting(Iso2631Filter::Wb, 1.0);
        assert!(previous < 0.45, "Wb(1 Hz) = {previous}");
        for f in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] {
            let w = iso_2631_weighting(Iso2631Filter::Wb, f);
            assert!(w > previous, "Wb must rise below the plateau at {f} Hz");
            previous = w;
        }
        // Falling edge above the plateau.
        assert!(iso_2631_weighting(Iso2631Filter::Wb, 31.5) < 0.6);
    }

    #[test]
    fn weightings_attenuate_both_band_edges() {
        for filter in [Iso2631Filter::Wk, Iso2631Filter::Wb, Iso2631Filter::Wd] {
            let low = iso_2631_weighting(filter, 0.05);
            let high = iso_2631_weighting(filter, 1000.0);
            assert!(low < 0.02, "{filter:?}(0.05 Hz) = {low}");
            assert!(high < 0.02, "{filter:?}(1 kHz) = {high}");
        }
    }

    #[test]
    fn weightings_match_laplace_domain_reference_evaluation() {
        // The closed-form magnitudes must agree with a direct complex-arithmetic
        // evaluation of the transfer function over four decades.
        let cases = [
            (Iso2631Filter::Wk, WK_TABLE),
            (Iso2631Filter::Wd, WD_TABLE),
            (Iso2631Filter::Wb, WB_TABLE),
        ];
        for (filter, table) in cases {
            for f in log_grid(0.02, 500.0, 400) {
                let got = iso_2631_weighting(filter, f);
                let expected = reference_weighting(f, table);
                assert!(
                    (got - expected).abs() <= 1e-10 * expected.max(1e-3),
                    "{filter:?}({f} Hz): got {got}, reference {expected}"
                );
            }
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
        assert!((w - iso_2631_weighting(Iso2631Filter::Wb, 5.0)).abs() < 1e-12);
    }

    #[test]
    fn weighted_rms_combines_bins_in_quadrature() {
        // One unit of acceleration at 1 Hz and one at 6.3 Hz through Wk.
        let w = weighted_rms(&[1.0, 6.3], &[1.0, 1.0], Iso2631Filter::Wk);
        let expected = (0.482_f64.powi(2) + 1.054_f64.powi(2)).sqrt();
        assert!(
            (w - expected).abs() < 0.005 * expected,
            "weighted rms {w}, expected {expected}"
        );
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

    /// Deterministic multi-tone test profile with a non-trivial mean.
    fn test_profile(n: usize) -> Vec<f64> {
        (0..n)
            .map(|j| {
                let x = j as f64;
                0.013 * (0.37 * x).sin()
                    + 0.004 * (1.9 * x + 0.5).sin()
                    + 0.002 * (((j * 37) % 11) as f64 - 5.0) / 5.0
                    + 0.05
            })
            .collect()
    }

    /// Independent `O(n²)` one-sided periodogram (the algorithm that used to be
    /// the implementation), with the correct unpaired-bin rule.
    fn naive_psd(profile: &[f64], dx: f64) -> RoadPsd {
        let n = profile.len();
        let mean = profile.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = profile.iter().map(|&x| x - mean).collect();
        let mut psd = Vec::new();
        let mut freq = Vec::new();
        for k in 0..(n / 2 + 1) {
            let (mut re, mut im) = (0.0_f64, 0.0_f64);
            for (j, &c) in centered.iter().enumerate() {
                let angle = -2.0 * PI * ((k * j) % n) as f64 / n as f64;
                re += c * angle.cos();
                im += c * angle.sin();
            }
            let unpaired = k == 0 || (n.is_multiple_of(2) && k == n / 2);
            let factor = if unpaired { 1.0 } else { 2.0 };
            psd.push(factor * (re * re + im * im) / n as f64 * dx);
            freq.push(k as f64 / (n as f64 * dx));
        }
        RoadPsd {
            spatial_freq: freq,
            psd,
        }
    }

    /// Mean square of the mean-subtracted profile.
    fn centred_mean_square(profile: &[f64]) -> f64 {
        let n = profile.len() as f64;
        let mean = profile.iter().sum::<f64>() / n;
        profile.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n
    }

    /// `Σ S(n_k) · Δn` of a one-sided PSD.
    fn psd_integral(psd: &RoadPsd, n: usize, dx: f64) -> f64 {
        let df = 1.0 / (n as f64 * dx);
        psd.psd.iter().sum::<f64>() * df
    }

    #[test]
    fn psd_matches_naive_dft_for_even_length() {
        let (n, dx) = (1000, 0.1);
        let profile = test_profile(n);
        let fast = psd_road_profile(&profile, dx);
        let slow = naive_psd(&profile, dx);
        assert_eq!(fast.psd.len(), n / 2 + 1);
        assert_eq!(fast.psd.len(), slow.psd.len());
        let scale = slow.psd.iter().cloned().fold(0.0_f64, f64::max);
        for (k, (a, b)) in fast.psd.iter().zip(slow.psd.iter()).enumerate() {
            assert!(
                (a - b).abs() <= 1e-9 * scale,
                "bin {k}: fft {a} vs naive {b} (scale {scale})"
            );
            assert!((fast.spatial_freq[k] - slow.spatial_freq[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn psd_matches_naive_dft_for_odd_and_prime_lengths() {
        for n in [3_usize, 5, 63, 999, 1009] {
            let dx = 0.25;
            let profile = test_profile(n);
            let fast = psd_road_profile(&profile, dx);
            let slow = naive_psd(&profile, dx);
            assert_eq!(fast.psd.len(), n / 2 + 1, "n = {n}");
            let scale = slow.psd.iter().cloned().fold(0.0_f64, f64::max);
            for (k, (a, b)) in fast.psd.iter().zip(slow.psd.iter()).enumerate() {
                assert!(
                    (a - b).abs() <= 1e-9 * scale,
                    "n = {n}, bin {k}: fft {a} vs naive {b}"
                );
            }
        }
    }

    #[test]
    fn psd_satisfies_parseval_for_even_and_odd_lengths() {
        // The integral of the one-sided PSD equals the variance of the profile.
        for n in [2_usize, 3, 4, 5, 64, 63, 1000, 999, 1009] {
            let dx = 0.1;
            let profile = test_profile(n);
            let psd = psd_road_profile(&profile, dx);
            let integral = psd_integral(&psd, n, dx);
            let variance = centred_mean_square(&profile);
            assert!(
                (integral - variance).abs() <= 1e-9 * variance.max(1e-30),
                "n = {n}: ∫PSD = {integral}, variance = {variance}"
            );
        }
    }

    #[test]
    fn psd_pure_sinusoid_concentrates_power_in_its_bin() {
        let (n, dx, k0, amp) = (512_usize, 0.1, 20_usize, 0.03);
        let profile: Vec<f64> = (0..n)
            .map(|j| amp * (2.0 * PI * (k0 * j) as f64 / n as f64).sin() + 0.2)
            .collect();
        let psd = psd_road_profile(&profile, dx);
        let df = 1.0 / (n as f64 * dx);
        // A sinusoid of amplitude A carries a mean-square power of A²/2.
        let power_in_bin = psd.psd[k0] * df;
        assert!(
            (power_in_bin - amp * amp / 2.0).abs() < 1e-9 * amp * amp,
            "bin power {power_in_bin}"
        );
        assert!((psd.spatial_freq[k0] - k0 as f64 * df).abs() < 1e-12);
        for (k, &v) in psd.psd.iter().enumerate() {
            if k != k0 {
                assert!(v * df < 1e-12 * amp * amp, "leakage in bin {k}: {v}");
            }
        }
    }

    #[test]
    fn psd_nyquist_bin_is_not_doubled_for_even_length() {
        // x_j = (-1)^j has all its power in the Nyquist bin: mean square 1.
        let (n, dx) = (8_usize, 0.5);
        let profile: Vec<f64> = (0..n)
            .map(|j| if j % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let psd = psd_road_profile(&profile, dx);
        assert_eq!(psd.psd.len(), n / 2 + 1);
        let nyq = psd.psd[n / 2] * (1.0 / (n as f64 * dx));
        assert!((nyq - 1.0).abs() < 1e-12, "Nyquist power {nyq}");
        assert!((psd_integral(&psd, n, dx) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn psd_last_bin_of_odd_length_is_doubled() {
        // Regression: for odd N the last bin k = (N-1)/2 has a mirror image and
        // must be doubled; it used to be treated as an unpaired Nyquist bin.
        let (n, dx) = (5_usize, 1.0);
        let profile: Vec<f64> = (0..n)
            .map(|j| (2.0 * PI * (2 * j) as f64 / n as f64).cos())
            .collect();
        let psd = psd_road_profile(&profile, dx);
        assert_eq!(psd.psd.len(), 3);
        let power = psd.psd[2] * (1.0 / (n as f64 * dx));
        // cos² has mean square 1/2.
        assert!((power - 0.5).abs() < 1e-12, "power in last bin {power}");
        assert!((psd_integral(&psd, n, dx) - centred_mean_square(&profile)).abs() < 1e-12);
    }

    #[test]
    fn psd_of_constant_profile_is_zero() {
        let psd = psd_road_profile(&vec![0.37; 129], 0.2);
        assert_eq!(psd.psd.len(), 65);
        assert!(psd.psd.iter().all(|&v| v < 1e-20), "{:?}", &psd.psd[..4]);
    }

    #[test]
    fn psd_scales_linearly_with_sample_spacing() {
        let profile = test_profile(256);
        let a = psd_road_profile(&profile, 0.1);
        let b = psd_road_profile(&profile, 0.4);
        for k in 1..a.psd.len() {
            assert!((b.psd[k] - 4.0 * a.psd[k]).abs() <= 1e-12 * a.psd[k].max(1e-30));
            assert!((a.spatial_freq[k] - 4.0 * b.spatial_freq[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn psd_handles_long_profiles_in_n_log_n_time() {
        // 100 000 samples would need 5e9 sin/cos pairs with the old O(n²) DFT.
        let (n, dx) = (100_000_usize, 0.05);
        let profile = test_profile(n);
        let psd = psd_road_profile(&profile, dx);
        assert_eq!(psd.psd.len(), n / 2 + 1);
        let integral = psd_integral(&psd, n, dx);
        let variance = centred_mean_square(&profile);
        assert!(
            (integral - variance).abs() <= 1e-8 * variance,
            "∫PSD = {integral}, variance = {variance}"
        );
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
