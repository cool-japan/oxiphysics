// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rational sample-rate conversion (polyphase FIR resampling).
//!
//! Converting a signal from `fs_in` to `fs_out = fs_in · up / down` is the
//! textbook chain "insert `up - 1` zeros between the samples, low-pass filter,
//! keep every `down`-th sample".  The low-pass is a symmetric Kaiser-windowed
//! sinc of `2 · 10 · max(up, down) + 1` taps whose cutoff `π / max(up, down)` (at
//! the intermediate rate `up · fs_in`) is the Nyquist frequency of the slower of
//! the two rates, so it both interpolates (images above `fs_in / 2`) and
//! anti-aliases (components above `fs_out / 2`).  The taps are scaled to sum to
//! `up`, which restores the amplitude lost by the zero stuffing.
//!
//! Only the kept output samples are computed, and only from the non-zero input
//! samples (polyphase evaluation):
//!
//! ```text
//! y[m] = Σ_k h[k] · u[m · down + half − k],     u[j · up] = x[j],  u[other] = 0,
//! ```
//!
//! where `half = (taps − 1) / 2` removes the group delay of the symmetric filter.
//! Output sample `m` therefore sits at the input time `m · down / up`: the first
//! output is aligned with the first input and the output has
//! `ceil(len · up / down)` samples.
//!
//! Beyond both ends the input is extended by repeating its first and last
//! sample, so a constant signal stays constant up to the filter's pass-band
//! ripple (below `1e-3`, which makes the DC gain depend slightly on the output
//! phase) all the way to the edges.

use super::bessel_i0;
use crate::error::{PhysicsError, PhysicsResult};
use std::f64::consts::PI;

/// Zero crossings of the windowed sinc on each side, in units of `max(up, down)` samples.
const HALF_LENGTH_PER_FACTOR: usize = 10;
/// Shape parameter of the Kaiser window of the anti-aliasing filter
/// (roughly 50 dB of stop-band attenuation).
const KAISER_BETA: f64 = 5.0;
/// Largest denominator of the rational approximation of a sample-rate ratio that
/// is not an exact quotient of two integer rates.
const MAX_APPROXIMATION_DENOMINATOR: u64 = 1000;
/// Largest relative error accepted when a sample-rate ratio is replaced by its
/// rational approximation.
const MAX_RELATIVE_APPROXIMATION_ERROR: f64 = 1.0e-2;
/// Largest up- or down-sampling factor; it bounds the filter to `20 · 2^16 + 1`
/// taps (about 10 MB).
const MAX_FACTOR: u64 = 1 << 16;
/// Largest integer (`2^53`) that an `f64` holds exactly.
const MAX_EXACT_INTEGER_RATE: f64 = 9_007_199_254_740_992.0;

fn invalid(field: &str, reason: String) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: field.to_string(),
        reason,
    }
}

/// Greatest common divisor (`gcd(a, 0) = a`).
fn gcd(a: u64, b: u64) -> u64 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A sample rate that is an exact positive integer below `2^53`.
fn integer_rate(rate: f64) -> Option<u64> {
    (rate.fract() == 0.0 && (1.0..=MAX_EXACT_INTEGER_RATE).contains(&rate)).then_some(rate as u64)
}

/// The fraction `p / q` closest to `value > 0`, in lowest terms, with `p >= 1`
/// and `q` at most `max_denominator`.
///
/// The continued-fraction convergents of `value` are generated until the
/// denominator would exceed the bound; the result is the better of the last
/// convergent and the closest semiconvergent (the usual "limit denominator"
/// construction).  `None` if `value` is not positive and finite, or if the best
/// fraction would be `0 / 1`.
fn best_rational_approximation(value: f64, max_denominator: u64) -> Option<(u64, u64)> {
    if !(value.is_finite() && value > 0.0) || max_denominator == 0 {
        return None;
    }
    // (p0 / q0, p1 / q1): the two latest convergents of the continued fraction.
    let (mut p0, mut q0, mut p1, mut q1) = (0_u64, 1_u64, 1_u64, 0_u64);
    let mut rest = value;
    loop {
        let term = rest.floor();
        // The cast saturates for absurdly large terms; the checked arithmetic
        // below then ends the expansion.
        let a = term as u64;
        let next_q = a.checked_mul(q1).and_then(|v| v.checked_add(q0));
        let next_p = a.checked_mul(p1).and_then(|v| v.checked_add(p0));
        match (next_p, next_q) {
            (Some(p2), Some(q2)) if q2 <= max_denominator => {
                (p0, q0, p1, q1) = (p1, q1, p2, q2);
            }
            _ => break,
        }
        let fraction = rest - term;
        if fraction == 0.0 {
            // The expansion terminates: `value` is exactly p1 / q1.
            return (p1 > 0).then_some((p1, q1));
        }
        rest = fraction.recip();
    }
    // Closest semiconvergent (p0 + k p1) / (q0 + k q1) with the largest admissible k.
    let k = (max_denominator - q0) / q1;
    let semiconvergent = k
        .checked_mul(p1)
        .and_then(|v| v.checked_add(p0))
        .zip(k.checked_mul(q1).and_then(|v| v.checked_add(q0)));
    let error = |(p, q): (u64, u64)| (value - p as f64 / q as f64).abs();
    let convergent = (p1, q1);
    let best = match semiconvergent {
        Some(candidate) if error(candidate) < error(convergent) => candidate,
        _ => convergent,
    };
    (best.0 > 0).then_some(best)
}

/// Convert a bounded factor to `usize`.
fn factor_to_usize(factor: u64, field: &str) -> PhysicsResult<usize> {
    usize::try_from(factor).map_err(|_| invalid(field, format!("factor {factor} is too large")))
}

/// Up- and down-sampling factors that convert from `fs_in` to `fs_out`.
///
/// Integer rates give the exact reduced ratio (`44100 -> 48000` is `160 / 147`);
/// otherwise, or when the exact ratio exceeds [`MAX_FACTOR`], the ratio is
/// replaced by its best rational approximation with a denominator of at most
/// [`MAX_APPROXIMATION_DENOMINATOR`].
fn rate_ratio_factors(fs_in: f64, fs_out: f64) -> PhysicsResult<(usize, usize)> {
    for (field, rate) in [("fs_in", fs_in), ("fs_out", fs_out)] {
        if !(rate.is_finite() && rate > 0.0) {
            return Err(invalid(
                field,
                format!("sample rate must be positive and finite, got {rate}"),
            ));
        }
    }
    if let (Some(rate_in), Some(rate_out)) = (integer_rate(fs_in), integer_rate(fs_out)) {
        let divisor = gcd(rate_in, rate_out);
        let (up, down) = (rate_out / divisor, rate_in / divisor);
        if up.max(down) <= MAX_FACTOR {
            return Ok((
                factor_to_usize(up, "fs_out")?,
                factor_to_usize(down, "fs_in")?,
            ));
        }
    }
    let ratio = fs_out / fs_in;
    let (up, down) = best_rational_approximation(ratio, MAX_APPROXIMATION_DENOMINATOR)
        .filter(|&(p, q)| {
            (p as f64 / q as f64 - ratio).abs() <= MAX_RELATIVE_APPROXIMATION_ERROR * ratio
        })
        .ok_or_else(|| {
            invalid(
                "fs_out",
                format!(
                    "the rate ratio {ratio} (fs_out {fs_out} / fs_in {fs_in}) cannot be \
                     approximated by a fraction with a denominator of at most \
                     {MAX_APPROXIMATION_DENOMINATOR}; use resample_poly with explicit factors"
                ),
            )
        })?;
    if up > MAX_FACTOR {
        return Err(invalid(
            "fs_out",
            format!("the rate ratio {ratio} needs an upsampling factor above {MAX_FACTOR}"),
        ));
    }
    Ok((
        factor_to_usize(up, "fs_out")?,
        factor_to_usize(down, "fs_in")?,
    ))
}

/// Anti-aliasing / interpolation filter of the polyphase resampler.
///
/// A symmetric Kaiser-windowed sinc with `2 · 10 · max(up, down) + 1` taps and
/// cutoff `π / max(up, down)`, scaled to sum to `up`.
fn design_resampling_filter(up: usize, down: usize) -> Vec<f64> {
    let factor = up.max(down);
    let half = HALF_LENGTH_PER_FACTOR * factor;
    let factor_f = factor as f64;
    let half_f = half as f64;
    let window_scale = bessel_i0(KAISER_BETA).recip();
    let mut taps = vec![0.0_f64; 2 * half + 1];
    // The filter is symmetric about `half`: evaluate the right half and mirror it.
    for offset in 0..=half {
        let t = offset as f64;
        let ideal = if offset == 0 {
            1.0 / factor_f
        } else {
            (PI * t / factor_f).sin() / (PI * t)
        };
        let ratio = t / half_f;
        let window = bessel_i0(KAISER_BETA * (1.0 - ratio * ratio).max(0.0).sqrt()) * window_scale;
        let tap = ideal * window;
        taps[half + offset] = tap;
        taps[half - offset] = tap;
    }
    let gain = up as f64 / taps.iter().sum::<f64>();
    for tap in &mut taps {
        *tap *= gain;
    }
    taps
}

/// Resample `signal` by the rational factor `up / down` with a polyphase FIR filter.
///
/// The output has `ceil(len · up / down)` samples; output sample `m` is the
/// band-limited interpolation of the input at the time `m · down / up` (in input
/// samples), so the first output coincides with the first input.  Content above
/// the Nyquist frequency of the lower of the two rates is removed by a
/// Kaiser-windowed sinc anti-aliasing filter (see the module documentation for
/// the design and the boundary treatment).  `up` and `down` are reduced by their
/// greatest common divisor first, and equal factors return the input unchanged.
///
/// The cost is proportional to `20 · max(len, output length)` multiply-adds,
/// independent of the factors.
///
/// # Errors
///
/// [`PhysicsError::InvalidInput`] if `up` or `down` is zero, if the reduced
/// factors exceed `65536`, or if the output length overflows `usize`.
pub fn resample_poly(signal: &[f64], up: usize, down: usize) -> PhysicsResult<Vec<f64>> {
    if up == 0 || down == 0 {
        return Err(invalid(
            "up",
            format!("resampling factors must be positive, got up = {up}, down = {down}"),
        ));
    }
    let divisor = factor_to_usize(gcd(up as u64, down as u64), "up")?;
    let (up, down) = (up / divisor, down / divisor);
    if up.max(down) as u64 > MAX_FACTOR {
        return Err(invalid(
            "up",
            format!("resampling factors {up} / {down} exceed the limit {MAX_FACTOR}"),
        ));
    }
    if up == down {
        return Ok(signal.to_vec());
    }
    let n = signal.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let out_len = n
        .checked_mul(up)
        .map(|stuffed| stuffed.div_ceil(down))
        .ok_or_else(|| {
            invalid(
                "signal",
                format!("the output length of {n} samples upsampled by {up} overflows"),
            )
        })?;
    let taps = design_resampling_filter(up, down);
    let half = (taps.len() - 1) / 2;
    let last = n - 1;
    let mut out = Vec::with_capacity(out_len);
    for m in 0..out_len {
        // Output `m` sits at `centre / up` input samples; its filter taps are
        // those of the intermediate-rate phase `centre mod up`.
        let centre = m * down;
        let phase = centre % up;
        let lead = (half + phase) / up;
        let newest = centre / up + lead;
        let first_tap = (half + phase) % up;
        let value: f64 = taps[first_tap..]
            .iter()
            .step_by(up)
            .enumerate()
            .map(|(t, &tap)| tap * signal[newest.saturating_sub(t).min(last)])
            .sum();
        out.push(value);
    }
    Ok(out)
}

/// Resample `signal` from the sample rate `fs_in` to `fs_out`.
///
/// The rate ratio `fs_out / fs_in` is realised as an up/down factor pair for
/// [`resample_poly`]: exactly (reduced by the greatest common divisor) when both
/// rates are integers, for example `44100 -> 48000` is `160 / 147` and
/// `1000 -> 441` is `441 / 1000`; otherwise as the best rational approximation
/// with a denominator of at most 1000 (found with continued fractions), so the
/// effective output rate is `fs_in · up / down`.  See [`resample_poly`] for the
/// filter, the alignment and the output length.
///
/// # Errors
///
/// [`PhysicsError::InvalidInput`] if a rate is not positive and finite, if the
/// ratio cannot be approximated within 1 % by a fraction with a denominator of
/// at most 1000 or needs a factor above `65536`, or for the reasons listed at
/// [`resample_poly`].
pub fn try_resample(signal: &[f64], fs_in: f64, fs_out: f64) -> PhysicsResult<Vec<f64>> {
    let (up, down) = rate_ratio_factors(fs_in, fs_out)?;
    resample_poly(signal, up, down)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `sin(0.05 i) + 0.5 cos(0.11 i)`, the input of the SciPy reference values.
    fn reference_signal(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| (0.05 * i as f64).sin() + 0.5 * (0.11 * i as f64).cos())
            .collect()
    }

    fn sine(frequency: f64, fs: f64, n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| (2.0 * PI * frequency * i as f64 / fs).sin())
            .collect()
    }

    /// Brute-force closest fraction with `p >= 1` and `q <= max_denominator`.
    fn brute_force_fraction(value: f64, max_denominator: u64) -> (u64, u64) {
        let mut best = (1_u64, 1_u64);
        let mut best_error = f64::INFINITY;
        for q in 1..=max_denominator {
            let p = ((value * q as f64).round() as u64).max(1);
            let error = (value - p as f64 / q as f64).abs();
            if error < best_error {
                best_error = error;
                best = (p, q);
            }
        }
        best
    }

    // ── rate ratio ────────────────────────────────────────────────────────

    #[test]
    fn integer_rates_give_the_exact_reduced_ratio() {
        let cases = [
            (1000.0, 441.0, (441, 1000)),
            (44100.0, 48000.0, (160, 147)),
            (48000.0, 44100.0, (147, 160)),
            (8000.0, 44100.0, (441, 80)),
            (22050.0, 44100.0, (2, 1)),
            (96000.0, 32000.0, (1, 3)),
            (500.0, 500.0, (1, 1)),
            (1.0, 40000.0, (40000, 1)),
        ];
        for (fs_in, fs_out, expected) in cases {
            assert_eq!(
                rate_ratio_factors(fs_in, fs_out).unwrap(),
                expected,
                "{fs_in} -> {fs_out}"
            );
        }
    }

    #[test]
    fn non_integer_rates_use_a_continued_fraction_approximation() {
        // 441.5 / 1000 = 883 / 2000 does not fit the denominator bound.
        let (up, down) = rate_ratio_factors(1000.0, 441.5).unwrap();
        assert!(down <= 1000);
        assert_eq!((up as u64, down as u64), brute_force_fraction(0.4415, 1000));
        // Fraction(0.4415).limit_denominator(1000) in Python.
        assert_eq!((up, down), (200, 453));
        assert!((up as f64 / down as f64 - 0.4415).abs() < 2e-6);
        let (up, down) = rate_ratio_factors(8000.0, 11025.5).unwrap();
        assert_eq!(
            (up as u64, down as u64),
            brute_force_fraction(11025.5 / 8000.0, 1000)
        );
        assert_eq!((up, down), (973, 706));
        // Ratios that are exact quotients of small integers survive the float round trip.
        assert_eq!(rate_ratio_factors(0.5, 0.25).unwrap(), (1, 2));
        assert_eq!(rate_ratio_factors(2.5, 7.5).unwrap(), (3, 1));
    }

    #[test]
    fn exact_integer_ratios_beyond_the_factor_limit_fall_back_to_the_approximation() {
        // 1_000_003 and 999_983 are coprime: the exact ratio would need factors of ~1e6.
        let (up, down) = rate_ratio_factors(999_983.0, 1_000_003.0).unwrap();
        assert!(down <= 1000);
        assert_eq!((up, down), (1, 1));
    }

    #[test]
    fn continued_fraction_matches_brute_force_search() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        for step in 0..3000 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let uniform = (state >> 11) as f64 / (1u64 << 53) as f64;
            // Ratios between 0.02 and 40, log-uniform.
            let value = 0.02 * (2000.0f64).powf(uniform);
            let (p, q) = best_rational_approximation(value, 1000).unwrap();
            let (brute_p, brute_q) = brute_force_fraction(value, 1000);
            let cf_error = (value - p as f64 / q as f64).abs();
            let brute_error = (value - brute_p as f64 / brute_q as f64).abs();
            assert!(q <= 1000 && gcd(p, q) == 1, "step {step}: {p}/{q}");
            assert!(
                cf_error <= brute_error * (1.0 + 1e-12) + 1e-18,
                "step {step}, value {value}: {p}/{q} (error {cf_error}) vs {brute_p}/{brute_q} (error {brute_error})"
            );
        }
    }

    #[test]
    fn best_rational_approximation_handles_exact_and_degenerate_inputs() {
        assert_eq!(best_rational_approximation(0.5, 1000), Some((1, 2)));
        assert_eq!(best_rational_approximation(3.0, 1000), Some((3, 1)));
        assert_eq!(best_rational_approximation(0.441, 1000), Some((441, 1000)));
        assert_eq!(
            best_rational_approximation(std::f64::consts::PI, 1000),
            Some((355, 113))
        );
        assert_eq!(
            best_rational_approximation(std::f64::consts::PI, 7),
            Some((22, 7))
        );
        assert_eq!(best_rational_approximation(0.0, 1000), None);
        assert_eq!(best_rational_approximation(-1.0, 1000), None);
        assert_eq!(best_rational_approximation(f64::NAN, 1000), None);
        assert_eq!(best_rational_approximation(f64::INFINITY, 1000), None);
        assert_eq!(best_rational_approximation(1.0e-4, 1000), None);
        assert_eq!(best_rational_approximation(0.5, 0), None);
    }

    #[test]
    fn invalid_rates_are_rejected() {
        let signal = [1.0; 8];
        for (fs_in, fs_out) in [
            (0.0, 100.0),
            (100.0, 0.0),
            (-1.0, 100.0),
            (100.0, -5.0),
            (f64::NAN, 100.0),
            (100.0, f64::NAN),
            (f64::INFINITY, 100.0),
            (100.0, f64::INFINITY),
            // ratio far below 1/1000 and far above the factor limit
            (1.0e6, 1.0),
            (1.0, 1.0e7),
            (1.0e300, 1.0e-300),
        ] {
            assert!(
                matches!(
                    try_resample(&signal, fs_in, fs_out),
                    Err(PhysicsError::InvalidInput { .. })
                ),
                "{fs_in} -> {fs_out} should be rejected"
            );
            assert!(super::super::resample(&signal, fs_in, fs_out).is_empty());
        }
    }

    // ── filter design ─────────────────────────────────────────────────────

    #[test]
    fn filter_is_symmetric_with_total_gain_up() {
        for (up, down) in [(3, 2), (2, 3), (441, 1000), (160, 147), (1, 7), (5, 1)] {
            let taps = design_resampling_filter(up, down);
            assert_eq!(taps.len(), 20 * up.max(down) + 1);
            let sum: f64 = taps.iter().sum();
            assert!(
                (sum - up as f64).abs() < 1e-9 * up as f64,
                "{up}/{down}: {sum}"
            );
            for (k, tap) in taps.iter().enumerate() {
                assert_eq!(*tap, taps[taps.len() - 1 - k], "{up}/{down}: tap {k}");
            }
        }
    }

    #[test]
    fn every_polyphase_branch_has_nearly_unit_dc_gain() {
        for (up, down) in [(441, 1000), (160, 147), (3, 2), (2, 1), (10, 1)] {
            let taps = design_resampling_filter(up, down);
            for phase in 0..up {
                let branch: f64 = taps.iter().skip(phase).step_by(up).sum();
                assert!(
                    (branch - 1.0).abs() < 1e-3,
                    "{up}/{down}: branch {phase} sums to {branch}"
                );
            }
        }
    }

    // ── resampling ────────────────────────────────────────────────────────

    #[test]
    fn output_length_is_the_ceiling_of_the_scaled_length() {
        for (up, down) in [
            (3, 2),
            (2, 3),
            (441, 1000),
            (160, 147),
            (1, 4),
            (5, 1),
            (7, 7),
        ] {
            for n in [1usize, 2, 3, 10, 99, 100, 257] {
                let y = resample_poly(&vec![1.0; n], up, down).unwrap();
                let divisor = gcd(up as u64, down as u64) as usize;
                let (u, d) = (up / divisor, down / divisor);
                assert_eq!(y.len(), (n * u).div_ceil(d), "n={n}, {up}/{down}");
            }
        }
        assert!(resample_poly(&[], 3, 2).unwrap().is_empty());
    }

    #[test]
    fn equal_rates_return_the_input_unchanged() {
        let x = reference_signal(50);
        assert_eq!(resample_poly(&x, 1, 1).unwrap(), x);
        assert_eq!(resample_poly(&x, 6, 6).unwrap(), x);
        assert_eq!(try_resample(&x, 48000.0, 48000.0).unwrap(), x);
    }

    #[test]
    fn matches_scipy_resample_poly_away_from_the_edges() {
        // scipy.signal.resample_poly(x, up, down, window=('kaiser', 5.0)) with
        // x[i] = sin(0.05 i) + 0.5 cos(0.11 i); SciPy zero-pads, so only samples
        // whose filter support lies inside the signal are comparable.
        let cases = [
            (
                (3, 2),
                400,
                100,
                [
                    0.05842575223639863,
                    -0.007000742969422204,
                    -0.07292290819811924,
                    -0.13927601899025321,
                    -0.2068024747296573,
                ],
            ),
            (
                (441, 1000),
                3000,
                600,
                [
                    -0.6754040405546236,
                    -0.5115453859983623,
                    -0.35657991744519146,
                    -0.21684319728792928,
                    -0.09712580849685896,
                ],
            ),
            (
                (160, 147),
                3000,
                1000,
                [
                    1.358308977160471,
                    1.3121836154709365,
                    1.2600565577051843,
                    1.2022981244325777,
                    1.1393253466638047,
                ],
            ),
            (
                (1, 4),
                400,
                40,
                [
                    1.14827589784317,
                    1.286781185773335,
                    1.3219938270540588,
                    1.2342197756739723,
                    1.0220567490524401,
                ],
            ),
        ];
        for ((up, down), n, first, expected) in cases {
            let y = resample_poly(&reference_signal(n), up, down).unwrap();
            for (offset, want) in expected.iter().enumerate() {
                let got = y[first + offset];
                assert!(
                    (got - want).abs() < 1e-12,
                    "{up}/{down}, output {}: {got} vs {want}",
                    first + offset
                );
            }
        }
    }

    #[test]
    fn polyphase_output_equals_zero_stuffing_filtering_and_decimation() {
        // Textbook chain: stuff zeros, convolve with the same taps, drop the
        // group delay, keep every down-th sample.  Edge-replication only changes
        // the outputs whose support reaches beyond the signal.
        for (up, down) in [(3, 2), (2, 3), (5, 1), (1, 4), (4, 3)] {
            let n = 150;
            let x = reference_signal(n);
            let taps = design_resampling_filter(up, down);
            let half = (taps.len() - 1) / 2;
            let mut stuffed = vec![0.0; n * up];
            for (j, v) in x.iter().enumerate() {
                stuffed[j * up] = *v;
            }
            let mut convolved = vec![0.0; stuffed.len() + taps.len() - 1];
            for (i, s) in stuffed.iter().enumerate() {
                if *s != 0.0 {
                    for (k, h) in taps.iter().enumerate() {
                        convolved[i + k] += s * h;
                    }
                }
            }
            let y = resample_poly(&x, up, down).unwrap();
            let margin = half / up + 2;
            let mut compared = 0;
            for (m, got) in y.iter().enumerate() {
                let centre = m * down;
                if centre / up < margin || centre / up + margin >= n {
                    continue;
                }
                let want = convolved[centre + half];
                assert!(
                    (got - want).abs() < 1e-12,
                    "{up}/{down}, m={m}: {got} vs {want}"
                );
                compared += 1;
            }
            assert!(
                compared > 10,
                "{up}/{down}: only {compared} samples compared"
            );
        }
    }

    #[test]
    fn sine_keeps_amplitude_and_phase_when_resampled_from_1000_to_441_hz() {
        // Regression: the old resampler returned ~0.4 % of the amplitude.
        let n = 2000;
        let y = try_resample(&sine(50.0, 1000.0, n), 1000.0, 441.0).unwrap();
        assert_eq!(y.len(), (n * 441).div_ceil(1000));
        let margin = 40;
        let interior = margin..y.len() - margin;
        // Pointwise: y[m] = sin(2 pi 50 m / 441), which also pins the alignment.
        let mut worst = 0.0_f64;
        for m in interior.clone() {
            let want = (2.0 * PI * 50.0 * m as f64 / 441.0).sin();
            worst = worst.max((y[m] - want).abs());
        }
        assert!(worst < 2e-3, "pointwise error {worst}");
        // Amplitude from the least-squares projection onto the 50 Hz sine and cosine.
        let (mut c, mut s) = (0.0, 0.0);
        for m in interior.clone() {
            let phase = 2.0 * PI * 50.0 * m as f64 / 441.0;
            c += y[m] * phase.cos();
            s += y[m] * phase.sin();
        }
        let amplitude = 2.0 * c.hypot(s) / interior.len() as f64;
        assert!((amplitude - 1.0).abs() < 0.01, "amplitude {amplitude}");
    }

    #[test]
    fn constants_stay_constant_including_the_edges() {
        for (up, down) in [
            (441, 1000),
            (160, 147),
            (147, 160),
            (3, 2),
            (2, 3),
            (1, 2),
            (2, 1),
            (10, 1),
            (1, 10),
        ] {
            for value in [1.0, -3.5] {
                let y = resample_poly(&[value; 300], up, down).unwrap();
                let worst = y
                    .iter()
                    .map(|v| (v / value - 1.0).abs())
                    .fold(0.0, f64::max);
                assert!(
                    worst < 1e-3,
                    "{up}/{down}, value {value}: deviation {worst}"
                );
            }
        }
    }

    #[test]
    fn components_above_the_new_nyquist_frequency_are_suppressed() {
        // A 300 Hz tone sampled at 1000 Hz would alias to 141 Hz at 441 Hz.
        let n = 2000;
        for frequency in [300.0, 350.0, 450.0] {
            let y = try_resample(&sine(frequency, 1000.0, n), 1000.0, 441.0).unwrap();
            let residual = y[40..y.len() - 40]
                .iter()
                .fold(0.0_f64, |m, v| m.max(v.abs()));
            assert!(residual < 5e-3, "{frequency} Hz leaves {residual}");
        }
    }

    #[test]
    fn interpolation_keeps_the_original_samples_and_fills_in_between() {
        // 100 -> 300 Hz: every third output is an input sample, the others lie on the sine.
        let n = 400;
        let x = sine(7.0, 100.0, n);
        let y = try_resample(&x, 100.0, 300.0).unwrap();
        assert_eq!(y.len(), 3 * n);
        for (m, value) in y.iter().enumerate().take(y.len() - 150).skip(150) {
            let want = (2.0 * PI * 7.0 * m as f64 / 300.0).sin();
            assert!((value - want).abs() < 1.5e-3, "m={m}: {value} vs {want}");
        }
        for (j, sample) in x.iter().enumerate().take(n - 60).skip(60) {
            assert!((y[3 * j] - sample).abs() < 1.5e-3, "sample {j}");
        }
    }

    #[test]
    fn a_pulse_keeps_its_position_when_resampled() {
        // Gaussian pulse centred on input sample 700: after 1000 -> 441 Hz it must
        // peak at output time 700 * 0.441 = 308.7 (a misaligned filter shifts it).
        let x: Vec<f64> = (0..2000)
            .map(|i| (-0.5 * ((i as f64 - 700.0) / 20.0).powi(2)).exp())
            .collect();
        let y = resample_poly(&x, 441, 1000).unwrap();
        let k = (0..y.len()).fold(0, |best, m| if y[m] > y[best] { m } else { best });
        let (a, b, c) = (y[k - 1], y[k], y[k + 1]);
        let peak = k as f64 + 0.5 * (a - c) / (a - 2.0 * b + c);
        assert!((peak - 308.7).abs() < 0.02, "peak at {peak}");
        assert!((b - 1.0).abs() < 0.01, "peak height {b}");
    }

    #[test]
    fn decimation_by_an_integer_factor_is_aligned_and_antialiased() {
        let n = 1200;
        let x: Vec<f64> = sine(20.0, 1000.0, n)
            .iter()
            .zip(sine(480.0, 1000.0, n).iter())
            .map(|(low, high)| low + high)
            .collect();
        let y = resample_poly(&x, 1, 4).unwrap();
        assert_eq!(y.len(), 300);
        for (m, value) in y.iter().enumerate().take(280).skip(20) {
            let want = (2.0 * PI * 20.0 * (4 * m) as f64 / 1000.0).sin();
            assert!((value - want).abs() < 5e-3, "m={m}: {value} vs {want}");
        }
    }

    #[test]
    fn single_sample_signals_are_extended_by_repetition() {
        let y = resample_poly(&[2.0], 3, 2).unwrap();
        assert_eq!(y.len(), 2);
        for v in &y {
            assert!((v - 2.0).abs() < 3e-3, "{v}");
        }
        assert_eq!(resample_poly(&[2.0], 1, 4).unwrap().len(), 1);
    }

    #[test]
    fn resample_poly_rejects_unusable_factors() {
        let x = [1.0, 2.0, 3.0];
        for (up, down) in [(0, 1), (1, 0), (0, 0), (65537, 1), (1, 65537), (131_072, 3)] {
            assert!(
                matches!(
                    resample_poly(&x, up, down),
                    Err(PhysicsError::InvalidInput { .. })
                ),
                "{up}/{down} should be rejected"
            );
        }
        // Common factors are removed first, so large equal multiples are fine.
        assert_eq!(resample_poly(&x, 100_000, 100_000).unwrap(), x.to_vec());
        assert_eq!(resample_poly(&x, 131_072, 65_536).unwrap().len(), 6);
    }
}
