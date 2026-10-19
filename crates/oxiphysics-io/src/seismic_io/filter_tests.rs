//! Tests for the zero-phase Butterworth bandpass of [`FrequencyFilter`].
//!
//! The responses are checked against the closed-form Butterworth magnitude,
//! against sine waves filtered in the time domain, and against reference
//! values produced by `scipy.signal.sosfiltfilt` with the same design.

use super::*;
use std::f64::consts::PI;

const FS: f64 = 100.0;

fn sine(freq: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| (2.0 * PI * freq * i as f64 / FS).sin())
        .collect()
}

fn band(order: usize) -> FrequencyFilter {
    FrequencyFilter::try_bandpass(1.0, 10.0, order, FS).unwrap()
}

/// Amplitude of the `freq` component of `window`, which starts at sample
/// `first_index` of the record (least-squares projection on sine and cosine).
fn amplitude(window: &[f64], freq: f64, first_index: usize) -> f64 {
    let (mut c, mut s) = (0.0f64, 0.0f64);
    for (k, &v) in window.iter().enumerate() {
        let phase = 2.0 * PI * freq * (first_index + k) as f64 / FS;
        c += v * phase.cos();
        s += v * phase.sin();
    }
    2.0 * c.hypot(s) / window.len() as f64
}

/// Closed-form squared magnitude `1 / (1 + x^(2N))` of the Butterworth
/// band-pass designed with pre-warped band edges, `x = (w^2 - w0^2)/(w bw)`.
fn closed_form_power(order: usize, fmin: f64, fmax: f64, freq: f64) -> f64 {
    let warp = |hz: f64| 2.0 * FS * (PI * hz / FS).tan();
    let (w_low, w_high, w) = (warp(fmin), warp(fmax), warp(freq));
    let w0_sq = w_low * w_high;
    let x = (w * w - w0_sq) / (w * (w_high - w_low));
    1.0 / (1.0 + x.powi(2 * order as i32))
}

#[test]
fn bandpass_gain_matches_the_butterworth_closed_form() {
    // Regression: the old "bandpass" was a difference of two moving averages
    // and ignored `order`.
    let freqs = [
        0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 20.0, 30.0, 40.0, 49.0,
    ];
    for order in 1..=6usize {
        let filter = band(order);
        for &f in &freqs {
            let got = filter.zero_phase_gain(f).unwrap();
            let want = closed_form_power(order, 1.0, 10.0, f);
            assert!(
                (got - want).abs() <= 1e-10 + 1e-8 * want,
                "order {order}, f = {f}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn corners_are_minus_three_db_per_pass_and_minus_six_db_zero_phase() {
    let filter = band(4);
    assert!((filter.zero_phase_gain(1.0).unwrap() - 0.5).abs() < 1e-9);
    assert!((filter.zero_phase_gain(10.0).unwrap() - 0.5).abs() < 1e-9);
    let single = filter.design().unwrap().magnitude_response(1.0, FS);
    assert!((single - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-9);
}

#[test]
fn passband_sine_has_unit_gain_and_no_phase_lag() {
    // 100 s at 100 Hz; compare away from the edges, where the record is long
    // enough for start-up transients to have died out.
    let n = 10_000;
    let input = sine(3.0, n);
    let output = band(4).apply(&input).unwrap();
    assert_eq!(output.len(), n);
    let worst = (2000..8000)
        .map(|i| (output[i] - input[i]).abs())
        .fold(0.0, f64::max);
    // A causal filter delays a 3 Hz wave by tens of degrees, which would show
    // up as errors of order 1 here.
    assert!(worst < 1e-6, "max |y - x| = {worst}");
}

#[test]
fn out_of_band_sines_are_attenuated_by_more_than_40_db() {
    let n = 10_000;
    for freq in [0.1, 40.0] {
        let input = sine(freq, n);
        let output = band(4).apply(&input).unwrap();
        let out_amp = amplitude(&output[2000..8000], freq, 2000);
        let ratio = out_amp / amplitude(&input[2000..8000], freq, 2000);
        assert!(ratio < 0.01, "f = {freq}: gain {ratio} is above -40 dB");
        assert!(ratio < 1e-6, "f = {freq}: gain {ratio}");
    }
}

#[test]
fn sine_gain_in_the_time_domain_matches_the_zero_phase_gain() {
    let n = 10_000;
    for order in [3usize, 4] {
        let filter = band(order);
        for freq in [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 25.0] {
            let output = filter.apply(&sine(freq, n)).unwrap();
            let measured = amplitude(&output[2000..8000], freq, 2000);
            let expected = filter.zero_phase_gain(freq).unwrap();
            assert!(
                (measured - expected).abs() <= 1e-6 + 1e-4 * expected,
                "order {order}, f = {freq}: measured {measured}, expected {expected}"
            );
        }
    }
}

#[test]
fn odd_order_with_a_wide_band_is_a_real_butterworth_filter() {
    // Order 3 over 1-10 Hz at 100 Hz: the real low-pass pole maps to two real
    // band-pass poles, the case that breaks naive conjugate pairing.
    let filter = band(3);
    assert_eq!(filter.design().unwrap().magnitude_response(0.0, FS), 0.0);
    let n = 10_000;
    let input = sine(3.0, n);
    let output = filter.apply(&input).unwrap();
    let worst = (2000..8000)
        .map(|i| (output[i] - input[i]).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-3, "max |y - x| = {worst}");
    let low = filter.zero_phase_gain(0.1).unwrap();
    assert!((low - closed_form_power(3, 1.0, 10.0, 0.1)).abs() < 1e-12);
    assert!(low < 1e-5);
}

#[test]
fn impulse_response_is_symmetric_about_the_impulse() {
    let n = 10_000;
    let centre = 5000;
    let mut impulse = vec![0.0; n];
    impulse[centre] = 1.0;
    let output = band(4).apply(&impulse).unwrap();
    let peak = output
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .map(|(i, _)| i);
    assert_eq!(peak, Some(centre));
    for k in 1..=400 {
        let (before, after) = (output[centre - k], output[centre + k]);
        assert!(
            (before - after).abs() <= 1e-12 * output[centre].abs(),
            "k = {k}: {before} vs {after}"
        );
    }
}

#[test]
fn dc_offset_is_removed_exactly_and_does_not_disturb_the_signal() {
    let filter = band(4);
    let flat = filter.apply(&vec![5.0; 500]).unwrap();
    assert!(flat.iter().all(|v| v.abs() < 1e-12), "{flat:?}");

    let wave = sine(3.0, 1000);
    let shifted: Vec<f64> = wave.iter().map(|v| v + 1.0e3).collect();
    let plain = filter.apply(&wave).unwrap();
    let offset = filter.apply(&shifted).unwrap();
    for (a, b) in plain.iter().zip(&offset) {
        assert!((a - b).abs() < 1e-9, "{a} vs {b}");
    }
}

#[test]
fn linear_trend_is_suppressed_away_from_the_edges() {
    let ramp: Vec<f64> = (0..10_000).map(|i| 0.01 * i as f64).collect();
    let output = band(4).apply(&ramp).unwrap();
    let worst = output[2000..8000]
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(worst < 1e-8, "max |y| = {worst}");
}

#[test]
fn short_records_are_filtered_without_panicking() {
    let filter = band(4);
    for n in [1usize, 2, 3, 10, 50, 199, 201] {
        let output = filter.apply(&sine(3.0, n)).unwrap();
        assert_eq!(output.len(), n);
        assert!(output.iter().all(|v| v.is_finite()), "n = {n}");
    }
    assert!(filter.apply(&[]).unwrap().is_empty());
}

/// `scipy.signal.sosfiltfilt(butter(N, [1, 10], 'band', fs=100, output='sos'),
/// x, padtype='odd', padlen=200)` for
/// `x[i] = 2 + sin(2 pi 3 i/fs) + 0.5 sin(2 pi 0.3 i/fs + 0.4) + 0.25 cos(2 pi 20 i/fs)`,
/// `i < 300`.
fn scipy_input() -> Vec<f64> {
    (0..300)
        .map(|i| {
            let t = i as f64;
            2.0 + (2.0 * PI * 3.0 * t / FS).sin()
                + 0.5 * (2.0 * PI * 0.3 * t / FS + 0.4).sin()
                + 0.25 * (2.0 * PI * 20.0 * t / FS).cos()
        })
        .collect()
}

fn assert_matches_reference(order: usize, reference: &[(usize, f64)]) {
    let output = band(order).apply(&scipy_input()).unwrap();
    for &(index, want) in reference {
        assert!(
            (output[index] - want).abs() < 1e-8,
            "order {order}, sample {index}: {} vs scipy {want}",
            output[index]
        );
    }
}

#[test]
fn order_4_output_matches_scipy_sosfiltfilt() {
    assert_matches_reference(
        4,
        &[
            (0, 0.001296857504944543),
            (1, 0.0880445331452252),
            (2, 0.18627025811078812),
            (5, 0.5640442530362475),
            (17, -0.15136618529516033),
            (50, 0.035302989712801214),
            (100, -0.01098738360304663),
            (150, 0.0019127625918425872),
            (250, -0.020550558755266846),
            (297, -0.4127872672169408),
            (298, -0.21109863609474108),
            (299, 0.005326432488262756),
        ],
    );
}

#[test]
fn order_3_output_matches_scipy_sosfiltfilt() {
    assert_matches_reference(
        3,
        &[
            (0, 0.0004233628322545535),
            (1, 0.08391746127996597),
            (2, 0.18220733461712568),
            (5, 0.5710836331060891),
            (17, -0.15045112735692934),
            (50, 0.03296796191744322),
            (100, -0.006283060513675387),
            (150, 0.0019642728982733892),
            (250, -0.01605902133904613),
            (297, -0.41913929213227147),
            (298, -0.21790187599415956),
            (299, 0.00039710715264185126),
        ],
    );
}

#[test]
fn zero_phase_gain_matches_scipy_frequency_response() {
    // |H|^2 of scipy.signal.butter(N, [1, 10], 'band', fs=100, output='sos').
    let order4 = [
        (0.1, 4.454579844880259e-09),
        (0.5, 0.002096860615393764),
        (1.0, 0.5000000000000233),
        (3.0, 0.9999999999800344),
        (5.0, 0.9998900913999296),
        (10.0, 0.5000000000000007),
        (20.0, 0.0008283338519257423),
        (40.0, 6.89831103440935e-09),
    ];
    let order3 = [
        (0.1, 5.452614893036211e-07),
        (0.5, 0.009718950019013222),
        (1.0, 0.5000000000000127),
        (3.0, 0.9999999905563854),
        (5.0, 0.9989276326430813),
        (10.0, 0.49999999999999944),
        (20.0, 0.004861914562170345),
        (40.0, 7.569318027124533e-07),
    ];
    for (order, table) in [(4usize, order4), (3usize, order3)] {
        let filter = band(order);
        for (freq, want) in table {
            let got = filter.zero_phase_gain(freq).unwrap();
            assert!(
                (got - want).abs() <= 1e-12 + 1e-7 * want,
                "order {order}, f = {freq}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn try_bandpass_rejects_unrealisable_parameters() {
    let bad: [(f64, f64, usize, f64); 14] = [
        (10.0, 1.0, 4, 100.0),          // reversed band
        (5.0, 5.0, 4, 100.0),           // empty band
        (0.0, 10.0, 4, 100.0),          // zero lower corner
        (-1.0, 10.0, 4, 100.0),         // negative lower corner
        (1.0, 50.0, 4, 100.0),          // upper corner at Nyquist
        (1.0, 60.0, 4, 100.0),          // upper corner above Nyquist
        (1.0, 10.0, 0, 100.0),          // no poles
        (1.0, 10.0, 65, 100.0),         // absurd order
        (1.0, 10.0, 4, 0.0),            // zero sampling rate
        (1.0, 10.0, 4, -100.0),         // negative sampling rate
        (f64::NAN, 10.0, 4, 100.0),     // NaN corner
        (1.0, f64::INFINITY, 4, 100.0), // infinite corner
        (1.0, 10.0, 4, f64::NAN),       // NaN sampling rate
        (1.0, 10.0, 4, f64::INFINITY),  // infinite sampling rate
    ];
    for (freqmin, freqmax, order, fs) in bad {
        let result = FrequencyFilter::try_bandpass(freqmin, freqmax, order, fs);
        assert!(
            matches!(result, Err(crate::Error::General(_))),
            "({freqmin}, {freqmax}, {order}, {fs}) was accepted"
        );
        let filter = FrequencyFilter::bandpass(freqmin, freqmax, order, fs);
        assert!(filter.apply(&sine(3.0, 100)).is_err());
        assert!(filter.zero_phase_gain(3.0).is_err());
    }
    assert!(FrequencyFilter::try_bandpass(1.0, 49.999, 4, 100.0).is_ok());
    assert!(FrequencyFilter::try_bandpass(0.001, 0.002, 1, 100.0).is_ok());
    assert!(FrequencyFilter::try_bandpass(1.0, 10.0, 64, 100.0).is_ok());
}

#[test]
fn unrealisable_filters_report_an_error_even_for_empty_records() {
    let filter = FrequencyFilter::bandpass(1.0, 60.0, 4, 100.0);
    assert!(filter.apply(&[]).is_err());
}
