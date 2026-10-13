// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the SURE threshold selection.
//!
//! The old risk expression was `n - 2k + Σ min(x², t²) + Σ_{|x|>t} (x² - t²)`
//! with the clamped sum cancelling against the remainder, which left
//! `n - 2k + Σ x²`: strictly decreasing in `k`, so the "optimal" threshold was
//! always the largest magnitude and soft thresholding wiped out every
//! coefficient.

use super::WaveletDenoising;

/// Deterministic uniform generator in (0, 1).
struct Uniform(u64);

impl Uniform {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    /// Standard normal sample (Box-Muller).
    fn normal(&mut self) -> f64 {
        let radius = (-2.0 * self.next().ln()).sqrt();
        radius * (2.0 * std::f64::consts::PI * self.next()).cos()
    }
}

/// Direct evaluation of `SURE(t) = n - 2 #{|z| <= t} + Σ min(z², t²)`.
fn direct_sure_risk(z: &[f64], t: f64) -> f64 {
    let below = z.iter().filter(|v| v.abs() <= t).count() as f64;
    let clipped: f64 = z.iter().map(|v| (v * v).min(t * t)).sum();
    z.len() as f64 - 2.0 * below + clipped
}

fn noise_with_spikes(n: usize, spikes: &[(usize, f64)], seed: u64) -> Vec<f64> {
    let mut rng = Uniform(seed);
    let mut signal: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    for &(index, amplitude) in spikes {
        signal[index] += amplitude;
    }
    signal
}

#[test]
fn documented_example_selects_the_noise_level_not_the_largest_coefficient() {
    // Risks of the candidates 0, 0.1, 0.2, 5 and 6 (sigma = 1):
    // 4, 2.04, 0.13, 48.05, 57.05.
    let t = WaveletDenoising::sure_threshold_sigma(&[0.1, 0.2, 5.0, 6.0], 1.0);
    assert!((t - 0.2).abs() < 1e-12, "threshold {t}");
    // Order and sign of the coefficients are irrelevant.
    let shuffled = WaveletDenoising::sure_threshold_sigma(&[6.0, -0.2, 5.0, -0.1], 1.0);
    assert!((shuffled - 0.2).abs() < 1e-12, "threshold {shuffled}");
}

#[test]
fn threshold_is_returned_in_the_units_of_the_input() {
    let t = WaveletDenoising::sure_threshold_sigma(&[0.2, 0.4, 10.0, 12.0], 2.0);
    assert!((t - 0.4).abs() < 1e-12, "threshold {t}");
}

#[test]
fn spikes_in_gaussian_noise_survive_the_selected_threshold() {
    let spikes = [
        (10, 10.0),
        (100, -10.0),
        (200, 10.0),
        (300, -10.0),
        (400, 10.0),
    ];
    let signal = noise_with_spikes(512, &spikes, 0x9e37_79b9_7f4a_7c15);
    let t = WaveletDenoising::sure_threshold(&signal);
    assert!(t > 0.5, "threshold {t} would keep a lot of noise");
    assert!(t < 5.0, "threshold {t} would also remove the spikes");
    for &(index, _) in &spikes {
        assert!(signal[index].abs() > t, "spike {index} would be removed");
    }
    // Most of the noise is below the threshold.
    let removed = signal
        .iter()
        .enumerate()
        .filter(|(i, _)| !spikes.iter().any(|&(s, _)| s == *i))
        .filter(|(_, v)| v.abs() <= t)
        .count();
    assert!(removed > 480, "only {removed} of 507 noise samples removed");
}

#[test]
fn soft_thresholding_with_the_selected_value_keeps_the_signal() {
    let spikes = [(7, 12.0), (70, -12.0)];
    let mut signal = noise_with_spikes(256, &spikes, 12345);
    let t = WaveletDenoising::sure_threshold(&signal);
    WaveletDenoising::soft_threshold(&mut signal, t);
    // The old threshold (the largest magnitude) zeroed everything.
    assert!(signal[7] > 8.0, "{}", signal[7]);
    assert!(signal[70] < -8.0, "{}", signal[70]);
}

#[test]
fn selected_threshold_minimises_the_risk_over_all_candidates() {
    let mut rng = Uniform(2024);
    for length in 1..=30 {
        let z: Vec<f64> = (0..length)
            .map(|i| rng.normal() * if i % 7 == 3 { 6.0 } else { 1.0 })
            .collect();
        let chosen = WaveletDenoising::sure_threshold_sigma(&z, 1.0);
        let mut best = (direct_sure_risk(&z, 0.0), 0.0);
        for value in &z {
            let candidate = value.abs();
            let risk = direct_sure_risk(&z, candidate);
            if risk < best.0 {
                best = (risk, candidate);
            }
        }
        assert!(
            (chosen - best.1).abs() < 1e-9,
            "length {length}: chose {chosen}, brute force found {} (risk {})",
            best.1,
            best.0
        );
    }
}

#[test]
fn sure_threshold_scales_linearly_with_the_input() {
    let spikes = [(5, 9.0), (50, -9.0)];
    let signal = noise_with_spikes(128, &spikes, 77);
    let reference = WaveletDenoising::sure_threshold(&signal);
    assert!(reference > 0.0);
    for scale in [1e-3, 4.0, 1e3] {
        let scaled: Vec<f64> = signal.iter().map(|v| v * scale).collect();
        let t = WaveletDenoising::sure_threshold(&scaled);
        assert!(
            (t - scale * reference).abs() < 1e-9 * scale * reference,
            "scale {scale}: {t} vs {}",
            scale * reference
        );
    }
}

#[test]
fn degenerate_inputs_give_a_zero_threshold_without_panicking() {
    assert_eq!(WaveletDenoising::sure_threshold(&[]), 0.0);
    assert_eq!(WaveletDenoising::sure_threshold(&[0.0; 16]), 0.0);
    assert_eq!(
        WaveletDenoising::sure_threshold_sigma(&[1.0, 2.0], 0.0),
        0.0
    );
    assert_eq!(
        WaveletDenoising::sure_threshold_sigma(&[1.0, 2.0], f64::NAN),
        0.0
    );
    assert_eq!(
        WaveletDenoising::sure_threshold_sigma(&[1.0, 2.0], -1.0),
        0.0
    );
    // Non-finite coefficients must not poison the search.
    let t = WaveletDenoising::sure_threshold_sigma(&[0.1, f64::NAN, 0.3, f64::INFINITY], 1.0);
    assert!(t.is_finite());
}

#[test]
fn threshold_search_is_not_quadratic() {
    // 200 000 coefficients would need 4e10 operations with the old per-candidate rescan.
    let signal = noise_with_spikes(200_000, &[(3, 15.0)], 99);
    let t = WaveletDenoising::sure_threshold(&signal);
    assert!(t.is_finite() && t > 0.0 && t < 15.0, "threshold {t}");
}
