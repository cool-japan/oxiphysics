//! Frequency-response regression tests for the Butterworth filters.
//!
//! The digital filters are designed with the bilinear transform and frequency
//! pre-warping, so their magnitude responses follow closed forms:
//!
//! * low-pass: `|H|² = 1 / (1 + (tan(π f/fs) / tan(π fc/fs))^(2N))`
//! * band-pass: `|H|² = 1 / (1 + ((w² − w0²) / (w · bw))^(2N))` with
//!   `w = 2 fs tan(π f/fs)`.

use super::types::{BiquadCoeff, ButterworthBandpass, ButterworthLowPass};
use crate::error::PhysicsError;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

const FS: f64 = 1000.0;

/// Pre-warped analogue centre frequency `w0 = √(w_low · w_high)` of a band.
fn analogue_centre(f_low: f64, f_high: f64, fs: f64) -> f64 {
    let w_low = 2.0 * fs * (PI * f_low / fs).tan();
    let w_high = 2.0 * fs * (PI * f_high / fs).tan();
    (w_low * w_high).sqrt()
}

/// Frequency (Hz) at which the digital band-pass has unit gain.
fn bandpass_centre_hz(f_low: f64, f_high: f64, fs: f64) -> f64 {
    fs / PI * (analogue_centre(f_low, f_high, fs) / (2.0 * fs)).atan()
}

/// Closed-form magnitude of the order-`n` Butterworth band-pass.
fn bandpass_closed_form(order: usize, f_low: f64, f_high: f64, f: f64, fs: f64) -> f64 {
    let w0 = analogue_centre(f_low, f_high, fs);
    let bw = 2.0 * fs * ((PI * f_high / fs).tan() - (PI * f_low / fs).tan());
    let w = 2.0 * fs * (PI * f / fs).tan();
    let x = (w * w - w0 * w0) / (w * bw);
    1.0 / (1.0 + x.powi(2 * order as i32)).sqrt()
}

/// Closed-form magnitude of the order-`n` Butterworth low-pass.
fn lowpass_closed_form(order: usize, cutoff: f64, f: f64, fs: f64) -> f64 {
    let ratio = (PI * f / fs).tan() / (PI * cutoff / fs).tan();
    1.0 / (1.0 + ratio.powi(2 * order as i32)).sqrt()
}

/// Amplitude of the component of `y` at `freq` Hz (least-squares projection).
fn amplitude_at(y: &[f64], freq: f64, fs: f64) -> f64 {
    let (mut c, mut s) = (0.0f64, 0.0f64);
    for (k, &v) in y.iter().enumerate() {
        let phase = 2.0 * PI * freq * k as f64 / fs;
        c += v * phase.cos();
        s += v * phase.sin();
    }
    2.0 * c.hypot(s) / y.len() as f64
}

fn sine(freq: f64, n: usize, fs: f64) -> Vec<f64> {
    (0..n)
        .map(|i| (2.0 * PI * freq * i as f64 / fs).sin())
        .collect()
}

// ── ButterworthLowPass ─────────────────────────────────────────────────────

#[test]
fn lowpass_odd_orders_settle_to_unit_dc_gain() {
    // Regression: the first-order section had the wrong feedback sign, giving
    // a DC gain of alpha / (1 - alpha) = 0.325 and a pole at z = -0.51.
    for order in [1usize, 3, 5] {
        let mut lp = ButterworthLowPass::new(order, 100.0, FS);
        let out = lp.apply(&vec![1.0; 2000]);
        assert!(
            (out[1999] - 1.0).abs() < 1e-6,
            "order {order}: step response settles at {}",
            out[1999]
        );
        assert!((lp.magnitude_response(0.0, FS) - 1.0).abs() < 1e-12);
    }
}

#[test]
fn lowpass_odd_order_sine_at_cutoff_has_half_power_amplitude() {
    for order in [1usize, 3, 5] {
        let mut lp = ButterworthLowPass::new(order, 100.0, FS);
        let out = lp.apply(&sine(100.0, 3000, FS));
        let amplitude = amplitude_at(&out[1000..], 100.0, FS);
        assert!(
            (amplitude - FRAC_1_SQRT_2).abs() < 1e-3,
            "order {order}: steady-state amplitude {amplitude}"
        );
    }
}

#[test]
fn lowpass_magnitude_matches_the_butterworth_closed_form_for_every_order() {
    for order in 1..=8usize {
        let lp = ButterworthLowPass::new(order, 100.0, FS);
        for f in [0.0, 10.0, 50.0, 100.0, 150.0, 200.0, 400.0] {
            let measured = lp.magnitude_response(f, FS);
            let expected = lowpass_closed_form(order, 100.0, f, FS);
            assert!(
                (measured - expected).abs() < 1e-10,
                "order {order}, f = {f}: {measured} vs {expected}"
            );
        }
        let at_cutoff = lp.magnitude_response(100.0, FS);
        assert!((at_cutoff - FRAC_1_SQRT_2).abs() < 1e-10, "order {order}");
    }
}

#[test]
fn lowpass_matches_scipy_reference_magnitudes() {
    // scipy.signal.butter(N, 100, 'low', fs=1000, output='sos') sampled at
    // 0, 50, 100, 200 and 400 Hz.
    let reference: [(usize, [f64; 5]); 3] = [
        (
            1,
            [
                1.0,
                0.8988915633764917,
                FRAC_1_SQRT_2,
                0.4082482904638629,
                0.1049893446976228,
            ],
        ),
        (
            3,
            [
                1.0,
                0.9933588135581032,
                FRAC_1_SQRT_2,
                0.08908708063747477,
                0.001176673385731213,
            ],
        ),
        (
            5,
            [
                1.0,
                0.9996214806165118,
                FRAC_1_SQRT_2,
                0.01788568233972406,
                1.311476114628855e-05,
            ],
        ),
    ];
    for (order, values) in reference {
        let lp = ButterworthLowPass::new(order, 100.0, FS);
        for (&f, &expected) in [0.0, 50.0, 100.0, 200.0, 400.0].iter().zip(values.iter()) {
            let measured = lp.magnitude_response(f, FS);
            assert!(
                (measured - expected).abs() < 1e-9,
                "order {order}, f = {f}: {measured} vs {expected}"
            );
        }
    }
}

// ── ButterworthBandpass ────────────────────────────────────────────────────

#[test]
fn bandpass_has_unit_centre_gain_and_half_power_at_both_edges() {
    // Regression: the old ad-hoc RBJ section reached only 0.41/0.50 (order 1)
    // and 0.004/0.010 (order 4) at f_low/f_high instead of 1/sqrt(2).
    let bands = [(100.0, 200.0), (50.0, 400.0), (10.0, 450.0), (5.0, 20.0)];
    for order in [1usize, 2, 3, 4, 6] {
        for (f_low, f_high) in bands {
            let bp = ButterworthBandpass::new(order, f_low, f_high, FS);
            assert_eq!(bp.sections.len(), order);
            let f0 = bandpass_centre_hz(f_low, f_high, FS);
            let centre = bp.magnitude_response(f0, FS);
            let lower = bp.magnitude_response(f_low, FS);
            let upper = bp.magnitude_response(f_high, FS);
            assert!(
                (centre - 1.0).abs() < 1e-9,
                "order {order}, band [{f_low}, {f_high}]: centre gain {centre}"
            );
            assert!(
                (lower - FRAC_1_SQRT_2).abs() < 1e-9,
                "order {order}, band [{f_low}, {f_high}]: lower edge {lower}"
            );
            assert!(
                (upper - FRAC_1_SQRT_2).abs() < 1e-9,
                "order {order}, band [{f_low}, {f_high}]: upper edge {upper}"
            );
        }
    }
}

#[test]
fn bandpass_magnitude_matches_the_butterworth_closed_form() {
    for order in [1usize, 2, 3, 4, 5, 6] {
        for (f_low, f_high) in [(100.0, 200.0), (50.0, 400.0), (10.0, 450.0)] {
            let bp = ButterworthBandpass::new(order, f_low, f_high, FS);
            for f in [20.0, 60.0, 99.0, 143.0, 180.0, 250.0, 330.0, 480.0] {
                let measured = bp.magnitude_response(f, FS);
                let expected = bandpass_closed_form(order, f_low, f_high, f, FS);
                assert!(
                    (measured - expected).abs() < 1e-9,
                    "order {order}, band [{f_low}, {f_high}], f = {f}: {measured} vs {expected}"
                );
            }
        }
    }
}

#[test]
fn bandpass_matches_scipy_reference_magnitudes() {
    // scipy.signal.butter(N, [100, 200], 'bandpass', fs=1000, output='sos')
    // sampled at 25, 50, 100, 150, 200, 300 and 450 Hz.
    let frequencies = [25.0, 50.0, 100.0, 150.0, 200.0, 300.0, 450.0];
    let reference: [(usize, [f64; 7]); 3] = [
        (
            1,
            [
                0.1362214060281352,
                0.2886636387023943,
                FRAC_1_SQRT_2,
                0.9934441359103479,
                FRAC_1_SQRT_2,
                0.3162277660168381,
                0.06385914116738943,
            ],
        ),
        (
            2,
            [
                0.01890373851335164,
                0.09052794439783991,
                FRAC_1_SQRT_2,
                0.9999123388597768,
                FRAC_1_SQRT_2,
                0.1104315260748466,
                0.004094653680716482,
            ],
        ),
        (
            4,
            [
                0.0003574790525627846,
                0.008262744696811703,
                FRAC_1_SQRT_2,
                0.9999999846270076,
                FRAC_1_SQRT_2,
                0.01234473828166941,
                1.676646987244717e-05,
            ],
        ),
    ];
    for (order, values) in reference {
        let bp = ButterworthBandpass::new(order, 100.0, 200.0, FS);
        for (&f, &expected) in frequencies.iter().zip(values.iter()) {
            let measured = bp.magnitude_response(f, FS);
            assert!(
                (measured - expected).abs() < 1e-9 * expected.max(1e-3),
                "order {order}, f = {f}: {measured} vs {expected}"
            );
        }
    }
}

#[test]
fn bandpass_rejects_dc_and_nyquist_exactly_and_rolls_off_monotonically() {
    for order in [1usize, 2, 4] {
        let bp = ButterworthBandpass::new(order, 100.0, 200.0, FS);
        assert!(bp.magnitude_response(0.0, FS) < 1e-12);
        assert!(bp.magnitude_response(FS / 2.0, FS) < 1e-12);
        let f0 = bandpass_centre_hz(100.0, 200.0, FS);
        // The response is flat to rounding error (1 - gain ~ 1e-23 for order 4)
        // within a fraction of a hertz of the centre, hence the slack.
        let slack = 1e-12;
        let mut previous = 0.0f64;
        for f in (1..f0 as usize).map(|f| f as f64) {
            let gain = bp.magnitude_response(f, FS);
            assert!(
                gain >= previous - slack,
                "order {order}: not monotone at {f} Hz"
            );
            previous = gain;
        }
        let mut previous = bp.magnitude_response(f0, FS);
        for f in (f0 as usize + 1..499).map(|f| f as f64) {
            let gain = bp.magnitude_response(f, FS);
            assert!(
                gain <= previous + slack,
                "order {order}: not monotone at {f} Hz"
            );
            previous = gain;
        }
    }
}

#[test]
fn bandpass_sections_have_zeros_at_plus_minus_one_and_stable_poles() {
    for order in [1usize, 2, 3, 4, 7] {
        for (f_low, f_high) in [(100.0, 200.0), (10.0, 450.0), (5.0, 20.0)] {
            let bp = ButterworthBandpass::new(order, f_low, f_high, FS);
            for section in &bp.sections {
                let c: BiquadCoeff = section.coeff;
                assert_eq!(c.b1, 0.0);
                assert!((c.b0 + c.b2).abs() < 1e-15 * c.b0.abs().max(1.0));
                assert!(c.a2.abs() < 1.0, "order {order}: |a2| = {}", c.a2.abs());
                assert!(
                    c.a1.abs() < 1.0 + c.a2,
                    "order {order}: poles outside the unit circle (a1 = {}, a2 = {})",
                    c.a1,
                    c.a2
                );
            }
        }
    }
}

#[test]
fn bandpass_time_domain_amplitudes_follow_the_analytic_response() {
    for order in [1usize, 2, 4] {
        let f0 = bandpass_centre_hz(100.0, 200.0, FS);
        for f in [f0, 100.0, 200.0, 40.0, 400.0] {
            let mut bp = ButterworthBandpass::new(order, 100.0, 200.0, FS);
            let out = bp.apply(&sine(f, 6000, FS));
            let amplitude = amplitude_at(&out[2000..], f, FS);
            let expected = bp.magnitude_response(f, FS);
            assert!(
                (amplitude - expected).abs() < 2e-3,
                "order {order}, f = {f}: measured {amplitude}, analytic {expected}"
            );
        }
    }
}

#[test]
fn bandpass_reset_restores_the_zero_state() {
    let mut bp = ButterworthBandpass::new(2, 100.0, 200.0, FS);
    let input = sine(150.0, 200, FS);
    let first = bp.apply(&input);
    bp.reset();
    let second = bp.apply(&input);
    assert_eq!(first, second);
}

#[test]
fn bandpass_try_new_rejects_invalid_specifications() {
    let invalid: [(usize, f64, f64, f64); 11] = [
        (0, 100.0, 200.0, 1000.0),
        (2, 200.0, 100.0, 1000.0),
        (2, 100.0, 100.0, 1000.0),
        (2, 0.0, 100.0, 1000.0),
        (2, -5.0, 100.0, 1000.0),
        (2, 100.0, 500.0, 1000.0),
        (2, 100.0, 600.0, 1000.0),
        (2, f64::NAN, 200.0, 1000.0),
        (2, 100.0, f64::INFINITY, 1000.0),
        (2, 100.0, 200.0, 0.0),
        (2, 100.0, 200.0, f64::NAN),
    ];
    for (order, f_low, f_high, fs) in invalid {
        let result = ButterworthBandpass::try_new(order, f_low, f_high, fs);
        assert!(
            matches!(result, Err(PhysicsError::InvalidInput { .. })),
            "({order}, {f_low}, {f_high}, {fs}) should be rejected"
        );
    }
    assert!(ButterworthBandpass::try_new(3, 100.0, 200.0, 1000.0).is_ok());
}

#[test]
fn bandpass_new_falls_back_to_a_pass_through_for_invalid_specifications() {
    let mut bp = ButterworthBandpass::new(2, 300.0, 100.0, FS);
    assert!(bp.sections.is_empty());
    let input = [1.0, -2.0, 3.0, 0.5];
    assert_eq!(bp.apply(&input), input.to_vec());
    // An order of zero is treated as one section.
    assert_eq!(
        ButterworthBandpass::new(0, 100.0, 200.0, FS).sections.len(),
        1
    );
}

#[test]
fn biquad_magnitude_response_matches_dc_and_nyquist_formulas() {
    let c = BiquadCoeff::new(0.2, 0.3, 0.1, -0.4, 0.25);
    let dc = (0.2 + 0.3 + 0.1) / (1.0 - 0.4 + 0.25);
    assert!((c.magnitude_response(0.0) - dc).abs() < 1e-12);
    let nyquist: f64 = ((0.2_f64 - 0.3 + 0.1) / (1.0 + 0.4 + 0.25)).abs();
    assert!((c.magnitude_response(PI) - nyquist).abs() < 1e-12);
}
