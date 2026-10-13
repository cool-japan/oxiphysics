// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Perfect-reconstruction tests of the periodic orthogonal DWT.
//!
//! The analysis step used to emit `(n + L - 1) / 2` coefficients per band, so
//! for every filter longer than two taps (db2 to db6) the last rows repeated
//! the first ones and the synthesis added them twice. Only the Haar round trip
//! was tested, so the defect went unnoticed.

use super::*;

const FAMILIES: [WaveletFamily; 6] = [
    WaveletFamily::Haar,
    WaveletFamily::Db2,
    WaveletFamily::Db3,
    WaveletFamily::Db4,
    WaveletFamily::Db5,
    WaveletFamily::Db6,
];

/// Deterministic pseudo-random signal in [-1, 1].
fn test_signal(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed.max(1);
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            2.0 * ((state >> 11) as f64 / (1u64 << 53) as f64) - 1.0
        })
        .collect()
}

fn max_abs_error(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "length mismatch");
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn energy(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum()
}

#[test]
fn single_level_bands_are_critically_sampled() {
    for family in FAMILIES {
        let even = dwt_single(&test_signal(8, 1), family);
        assert_eq!(even.approx.len(), 4, "{family:?}");
        assert_eq!(even.detail.len(), 4, "{family:?}");
        let odd = dwt_single(&test_signal(7, 1), family);
        assert_eq!(odd.approx.len(), 4, "{family:?} odd length");
        assert_eq!(odd.detail.len(), 4, "{family:?} odd length");
    }
}

#[test]
fn haar_single_level_matches_the_hand_computation() {
    let level = dwt_single(&[1.0, 3.0, 5.0, 7.0], WaveletFamily::Haar);
    let s = std::f64::consts::FRAC_1_SQRT_2;
    assert!((level.approx[0] - 4.0 * s).abs() < 1e-14);
    assert!((level.approx[1] - 12.0 * s).abs() < 1e-14);
    assert!((level.detail[0] - 2.0 * s).abs() < 1e-14);
    assert!((level.detail[1] - 2.0 * s).abs() < 1e-14);
}

#[test]
fn single_level_round_trip_is_exact_for_every_family_and_length() {
    for family in FAMILIES {
        for n in 1..=40 {
            let signal = test_signal(n, 7 + n as u64);
            let level = dwt_single(&signal, family);
            let rebuilt = idwt_single(&level, family, n);
            assert!(
                max_abs_error(&signal, &rebuilt) < 1e-10,
                "{family:?}, n = {n}: error {}",
                max_abs_error(&signal, &rebuilt)
            );
        }
    }
}

#[test]
fn even_length_transform_preserves_energy() {
    for family in FAMILIES {
        let signal = test_signal(64, 11);
        let level = dwt_single(&signal, family);
        let coefficient_energy = energy(&level.approx) + energy(&level.detail);
        assert!(
            (coefficient_energy - energy(&signal)).abs() < 1e-10 * energy(&signal),
            "{family:?}: {coefficient_energy} vs {}",
            energy(&signal)
        );
    }
}

#[test]
fn multi_level_round_trip_is_exact_including_odd_intermediate_lengths() {
    for family in FAMILIES {
        for n in [16, 64, 37] {
            for levels in [1, 2, 3] {
                let signal = test_signal(n, 100 + n as u64);
                let decomposition = dwt(&signal, family, levels);
                let rebuilt = idwt(&decomposition);
                assert!(
                    max_abs_error(&signal, &rebuilt) < 1e-10,
                    "{family:?}, n = {n}, levels = {levels}: error {}",
                    max_abs_error(&signal, &rebuilt)
                );
            }
        }
    }
}

#[test]
fn odd_length_decomposition_records_the_lengths_it_needs_for_synthesis() {
    let signal = test_signal(37, 5);
    let decomposition = dwt(&signal, WaveletFamily::Db2, 3);
    assert_eq!(decomposition.lengths, vec![37, 19, 10]);
    let detail_lengths: Vec<usize> = decomposition.details.iter().map(Vec::len).collect();
    assert_eq!(detail_lengths, vec![19, 10, 5]);
    assert_eq!(decomposition.approx.len(), 5);
}

#[test]
fn signals_shorter_than_the_filter_are_left_undecomposed() {
    let signal = test_signal(3, 2);
    let decomposition = dwt(&signal, WaveletFamily::Db4, 2);
    assert!(decomposition.details.is_empty());
    assert_eq!(idwt(&decomposition), signal);
}

#[test]
fn denoise_with_a_zero_threshold_returns_the_input() {
    for family in FAMILIES {
        for mode in [ThresholdMode::Hard, ThresholdMode::Soft] {
            for n in [64, 37] {
                let signal = test_signal(n, 21);
                let out = wavelet_denoise(&signal, family, 3, mode, Some(0.0));
                assert!(
                    max_abs_error(&signal, &out) < 1e-10,
                    "{family:?} {mode:?} n = {n}"
                );
            }
        }
    }
}

#[test]
fn compress_keeping_every_coefficient_returns_the_input() {
    for family in FAMILIES {
        for n in [64, 37] {
            let signal = test_signal(n, 33);
            let out = wavelet_compress(&signal, family, 3, 1.0);
            assert!(max_abs_error(&signal, &out) < 1e-10, "{family:?} n = {n}");
        }
    }
}

#[test]
fn compress_with_a_smaller_fraction_loses_only_a_little_energy() {
    let signal: Vec<f64> = (0..64)
        .map(|i| (i as f64 * 0.2).sin() + 0.5 * (i as f64 * 0.05).cos())
        .collect();
    let out = wavelet_compress(&signal, WaveletFamily::Db4, 3, 0.5);
    let error: f64 = signal
        .iter()
        .zip(out.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    assert!(
        error < 0.05 * energy(&signal),
        "relative error {}",
        error / energy(&signal)
    );
}

#[test]
fn multiresolution_analysis_ends_with_the_original_signal() {
    for family in FAMILIES {
        for n in [32, 37] {
            let signal = test_signal(n, 44);
            let mra = multiresolution_analysis(&signal, family, 3);
            let levels_done = mra.detail_contributions.len();
            assert!(levels_done >= 1, "{family:?} n = {n}");
            let finest = mra
                .approximations
                .last()
                .expect("an analysis always has at least the input resolution");
            assert!(max_abs_error(&signal, finest) < 1e-10, "{family:?} n = {n}");
            // The signal splits, by linearity, into the expanded approximation
            // that enters the finest synthesis step and the finest detail
            // contribution.
            let entering = &mra.approximations[levels_done - 1];
            let smooth = idwt_single(
                &DwtLevel {
                    approx: entering.clone(),
                    detail: vec![0.0; entering.len()],
                },
                family,
                n,
            );
            let sum: Vec<f64> = smooth
                .iter()
                .zip(mra.detail_contributions[0].iter())
                .map(|(a, d)| a + d)
                .collect();
            assert!(max_abs_error(&signal, &sum) < 1e-10, "{family:?} n = {n}");
        }
    }
}

#[test]
fn wavelet_packet_nodes_halve_in_length_at_every_level() {
    // The packet tree reuses the analysis step; with the old redundant
    // output the nodes of an orthogonal wavelet grew instead of halving.
    for family in FAMILIES {
        let signal = test_signal(64, 55);
        let tree = wavelet_packet_decompose(&signal, family, 3);
        for (level, nodes) in tree.nodes.iter().enumerate() {
            assert_eq!(nodes.len(), 1 << level, "{family:?} level {level}");
            for node in nodes {
                assert_eq!(
                    node.coefficients.len(),
                    64 >> level,
                    "{family:?} level {level}"
                );
            }
        }
        // The orthogonal packet transform conserves energy at every level.
        for nodes in &tree.nodes {
            let total: f64 = nodes.iter().map(|n| energy(&n.coefficients)).sum();
            assert!(
                (total - energy(&signal)).abs() < 1e-10 * energy(&signal),
                "{family:?}: {total} vs {}",
                energy(&signal)
            );
        }
    }
}

#[test]
fn empty_input_is_handled_without_panicking() {
    assert!(convolve_downsample(&[], &[0.5, 0.5]).is_empty());
    assert!(convolve_downsample(&[1.0, 2.0], &[]).is_empty());
    assert!(upsample_convolve(&[1.0], &[0.5, 0.5], 0).is_empty());
    let level = dwt_single(&[], WaveletFamily::Db3);
    assert!(level.approx.is_empty() && level.detail.is_empty());
    assert!(idwt_single(&level, WaveletFamily::Db3, 0).is_empty());
    let decomposition = dwt(&[], WaveletFamily::Haar, 3);
    assert!(idwt(&decomposition).is_empty());
}
