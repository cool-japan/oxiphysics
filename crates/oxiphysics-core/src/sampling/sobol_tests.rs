//! Regression tests for the Sobol generators of the sampling module.
//!
//! All generators share the Joe-Kuo direction numbers of `crate::random`.  The
//! reference points were produced with SciPy's
//! `qmc.Sobol(d=128, scramble=False, bits=32)` and are stored as the exact
//! integer states (the coordinate value scaled by `2^32`).

use super::functions::{sobol_sequence, sobol_sequence_10d};
use super::functions_2::sobol_scrambled;
use super::types::{Sobol, SobolSequence};
use crate::random::{RandomSampler, SOBOL_MAX_DIM, sobol_directions, sobol_state_at};
use std::collections::HashSet;

/// Integer state (value scaled by `2^32`) of a sample in `[0, 1)`.
fn to_state(x: f64) -> u32 {
    (x * 4_294_967_296.0) as u32
}

/// Asserts that, in every coordinate, the first `2^m` points fall exactly once
/// into each of the `2^m` equal bins of `[0, 1)` for all `m <= max_m`.
fn assert_dyadic_stratification(points: &[Vec<f64>], dims: usize, max_m: u32, label: &str) {
    for dim in 0..dims {
        for m in 0..=max_m {
            let bins = 1usize << m;
            let mut counts = vec![0u32; bins];
            for point in &points[..bins] {
                let bin = (point[dim] * bins as f64) as usize;
                counts[bin] += 1;
            }
            let empty = counts.iter().filter(|&&c| c == 0).count();
            let crowded = counts.iter().filter(|&&c| c > 1).count();
            assert!(
                empty == 0 && crowded == 0,
                "{label}: coordinate {dim}, first {bins} points: {empty} empty and {crowded} \
                 over-full bins out of {bins}"
            );
        }
    }
}

/// First eight points of the first ten coordinates.
const FIRST_POINTS_10D: [[u32; 10]; 8] = [
    [0; 10],
    [0x8000_0000; 10],
    [
        0xc000_0000,
        0x4000_0000,
        0x4000_0000,
        0x4000_0000,
        0xc000_0000,
        0xc000_0000,
        0x4000_0000,
        0xc000_0000,
        0xc000_0000,
        0xc000_0000,
    ],
    [
        0x4000_0000,
        0xc000_0000,
        0xc000_0000,
        0xc000_0000,
        0x4000_0000,
        0x4000_0000,
        0xc000_0000,
        0x4000_0000,
        0x4000_0000,
        0x4000_0000,
    ],
    [
        0x6000_0000,
        0x6000_0000,
        0xa000_0000,
        0xe000_0000,
        0x6000_0000,
        0x2000_0000,
        0x6000_0000,
        0xe000_0000,
        0xe000_0000,
        0xa000_0000,
    ],
    [
        0xe000_0000,
        0xe000_0000,
        0x2000_0000,
        0x6000_0000,
        0xe000_0000,
        0xa000_0000,
        0xe000_0000,
        0x6000_0000,
        0x6000_0000,
        0x2000_0000,
    ],
    [
        0xa000_0000,
        0x2000_0000,
        0xe000_0000,
        0xa000_0000,
        0xa000_0000,
        0xe000_0000,
        0x2000_0000,
        0x2000_0000,
        0x2000_0000,
        0x6000_0000,
    ],
    [
        0x2000_0000,
        0xa000_0000,
        0x6000_0000,
        0x2000_0000,
        0x2000_0000,
        0x6000_0000,
        0xa000_0000,
        0xa000_0000,
        0xa000_0000,
        0xe000_0000,
    ],
];

/// Points 1000 to 1007 of coordinates 0 to 9.  These depend on every
/// initial direction integer of the table rows (`s <= 5` for the first ten
/// coordinates) and on the polynomial recurrences.
const POINTS_1000_FIRST_TEN: [[u32; 10]; 8] = [
    [
        0x3840_0000,
        0x18c0_0000,
        0x84c0_0000,
        0xad40_0000,
        0x47c0_0000,
        0xe840_0000,
        0x0bc0_0000,
        0xe640_0000,
        0x8040_0000,
        0x11c0_0000,
    ],
    [
        0xb840_0000,
        0x98c0_0000,
        0x04c0_0000,
        0x2d40_0000,
        0xc7c0_0000,
        0x6840_0000,
        0x8bc0_0000,
        0x6640_0000,
        0x0040_0000,
        0x91c0_0000,
    ],
    [
        0xf840_0000,
        0x58c0_0000,
        0xc4c0_0000,
        0xed40_0000,
        0x87c0_0000,
        0x2840_0000,
        0x4bc0_0000,
        0x2640_0000,
        0x4040_0000,
        0xd1c0_0000,
    ],
    [
        0x7840_0000,
        0xd8c0_0000,
        0x44c0_0000,
        0x6d40_0000,
        0x07c0_0000,
        0xa840_0000,
        0xcbc0_0000,
        0xa640_0000,
        0xc040_0000,
        0x51c0_0000,
    ],
    [
        0x5840_0000,
        0x78c0_0000,
        0x24c0_0000,
        0x4d40_0000,
        0x27c0_0000,
        0xc840_0000,
        0x6bc0_0000,
        0x0640_0000,
        0x6040_0000,
        0xb1c0_0000,
    ],
    [
        0xd840_0000,
        0xf8c0_0000,
        0xa4c0_0000,
        0xcd40_0000,
        0xa7c0_0000,
        0x4840_0000,
        0xebc0_0000,
        0x8640_0000,
        0xe040_0000,
        0x31c0_0000,
    ],
    [
        0x9840_0000,
        0x38c0_0000,
        0x64c0_0000,
        0x0d40_0000,
        0xe7c0_0000,
        0x0840_0000,
        0x2bc0_0000,
        0xc640_0000,
        0xa040_0000,
        0x71c0_0000,
    ],
    [
        0x1840_0000,
        0xb8c0_0000,
        0xe4c0_0000,
        0x8d40_0000,
        0x67c0_0000,
        0x8840_0000,
        0xabc0_0000,
        0x4640_0000,
        0x2040_0000,
        0xf1c0_0000,
    ],
];

/// Coordinates probed by [`POINTS_1000_HIGH_DIMS`].
const HIGH_DIMS: [usize; 5] = [20, 21, 63, 100, 127];

/// Points 1000 to 1007 of the coordinates listed in [`HIGH_DIMS`].
const POINTS_1000_HIGH_DIMS: [[u32; 5]; 8] = [
    [
        0x85c0_0000,
        0x8cc0_0000,
        0x7240_0000,
        0xc4c0_0000,
        0x8b40_0000,
    ],
    [
        0x05c0_0000,
        0x0cc0_0000,
        0xf240_0000,
        0x44c0_0000,
        0x0b40_0000,
    ],
    [
        0xc5c0_0000,
        0x4cc0_0000,
        0xb240_0000,
        0x84c0_0000,
        0xcb40_0000,
    ],
    [
        0x45c0_0000,
        0xccc0_0000,
        0x3240_0000,
        0x04c0_0000,
        0x4b40_0000,
    ],
    [
        0xa5c0_0000,
        0xacc0_0000,
        0x5240_0000,
        0x64c0_0000,
        0x2b40_0000,
    ],
    [
        0x25c0_0000,
        0x2cc0_0000,
        0xd240_0000,
        0xe4c0_0000,
        0xab40_0000,
    ],
    [
        0xe5c0_0000,
        0x6cc0_0000,
        0x9240_0000,
        0x24c0_0000,
        0x6b40_0000,
    ],
    [
        0x65c0_0000,
        0xecc0_0000,
        0x1240_0000,
        0xa4c0_0000,
        0xeb40_0000,
    ],
];

// -- Stratification: the property the fabricated direction numbers violated --

#[test]
fn sobol_sequence_10d_is_stratified_in_every_coordinate() {
    // Regression: the invented direction numbers filled only 64 of 256 bins in
    // coordinates 2 and 3 and 128 of 256 bins in coordinates 4 to 8.
    let points = sobol_sequence_10d(4096);
    assert_eq!(points.len(), 4096);
    assert!(points.iter().all(|p| p.len() == 10));
    assert_dyadic_stratification(&points, 10, 12, "sobol_sequence_10d");
}

#[test]
fn every_tabulated_dimension_is_stratified() {
    let points = sobol_sequence(4096, SOBOL_MAX_DIM);
    assert!(points.iter().all(|p| p.len() == SOBOL_MAX_DIM));
    assert_dyadic_stratification(&points, SOBOL_MAX_DIM, 12, "sobol_sequence");
    let sequence = SobolSequence::new(SOBOL_MAX_DIM).sample(4096);
    assert_dyadic_stratification(&sequence, SOBOL_MAX_DIM, 12, "SobolSequence::sample");
}

#[test]
fn digitally_shifted_sequence_keeps_the_stratification() {
    let points = sobol_scrambled(4096, SOBOL_MAX_DIM);
    assert!(points.iter().all(|p| p.len() == SOBOL_MAX_DIM));
    assert_dyadic_stratification(&points, SOBOL_MAX_DIM, 12, "sobol_scrambled");
}

#[test]
fn one_dimensional_sobol_is_stratified_and_is_the_first_coordinate() {
    let mut sobol = Sobol::new_1d();
    let values: Vec<f64> = (0..4096).map(|_| sobol.next_1d()).collect();
    let as_points: Vec<Vec<f64>> = values.iter().map(|&v| vec![v]).collect();
    assert_dyadic_stratification(&as_points, 1, 12, "Sobol::next_1d");
    let reference = sobol_sequence(4096, 1);
    for (value, point) in values.iter().zip(reference.iter()) {
        assert_eq!(*value, point[0]);
    }
    assert_eq!(
        Sobol::new_1d().direction_numbers,
        vec![
            sobol_directions(0)
                .expect("coordinate 0 is tabulated")
                .to_vec()
        ]
    );
}

#[test]
fn first_two_coordinates_form_a_zero_net() {
    // Joe-Kuo coordinates 0 and 1 are a (0, m, 2)-net in base 2: each of the
    // 2^a x 2^(m-a) elementary boxes holds exactly one of the first 2^m points.
    let points = sobol_sequence(1 << 10, 2);
    for m in 1u32..=10 {
        let n = 1usize << m;
        for a in 0..=m {
            let (cols, rows) = (1usize << a, 1usize << (m - a));
            let mut counts = vec![0u32; n];
            for point in &points[..n] {
                let col = (point[0] * cols as f64) as usize;
                let row = (point[1] * rows as f64) as usize;
                counts[row * cols + col] += 1;
            }
            assert!(
                counts.iter().all(|&c| c == 1),
                "m = {m}, a = {a}: some elementary box is empty or over-full"
            );
        }
    }
}

#[test]
fn pairs_of_coordinates_are_uncorrelated() {
    // E[x_a x_b] = 1/4 for independent uniforms; fabricated, mutually
    // dependent coordinates deviate visibly.
    let points = sobol_sequence_10d(4096);
    for a in 0..10 {
        for b in (a + 1)..10 {
            let mean: f64 = points.iter().map(|p| p[a] * p[b]).sum::<f64>() / 4096.0;
            assert!(
                (mean - 0.25).abs() < 2e-3,
                "coordinates {a} and {b}: E[x_a x_b] = {mean}"
            );
        }
    }
}

#[test]
fn quasi_monte_carlo_integral_of_a_product_integrand_is_accurate() {
    // Integral of prod_{k=2..8} 2 x_k over the unit cube is 1.
    let points = sobol_sequence_10d(4096);
    let estimate: f64 = points
        .iter()
        .map(|p| p[2..9].iter().map(|x| 2.0 * x).product::<f64>())
        .sum::<f64>()
        / 4096.0;
    assert!((estimate - 1.0).abs() < 1e-2, "estimate = {estimate}");
}

// -- Published values ---------------------------------------------------------

#[test]
fn leading_points_match_the_reference_implementation() {
    let points = sobol_sequence_10d(8);
    for (index, (point, expected)) in points.iter().zip(FIRST_POINTS_10D.iter()).enumerate() {
        let states: Vec<u32> = point.iter().map(|&x| to_state(x)).collect();
        assert_eq!(states, expected.to_vec(), "point {index}");
    }
}

#[test]
fn points_1000_to_1007_match_the_reference_implementation() {
    let points = sobol_sequence(1008, SOBOL_MAX_DIM);
    for (offset, expected) in POINTS_1000_FIRST_TEN.iter().enumerate() {
        let point = &points[1000 + offset];
        let states: Vec<u32> = point[..10].iter().map(|&x| to_state(x)).collect();
        assert_eq!(states, expected.to_vec(), "point {}", 1000 + offset);
    }
    for (offset, expected) in POINTS_1000_HIGH_DIMS.iter().enumerate() {
        let point = &points[1000 + offset];
        let states: Vec<u32> = HIGH_DIMS.iter().map(|&dim| to_state(point[dim])).collect();
        assert_eq!(states, expected.to_vec(), "point {}", 1000 + offset);
    }
}

// -- One generator for all entry points --------------------------------------

#[test]
fn sampling_sobol_is_bit_identical_to_the_random_module_generator() {
    for dims in [1usize, 2, 3, 6, 10, 21, 64, SOBOL_MAX_DIM] {
        let expected = RandomSampler::sobol(300, dims).expect("dimension is tabulated");
        assert_eq!(expected.len(), 300);
        assert_eq!(sobol_sequence(300, dims), expected, "dims = {dims}");
        assert_eq!(
            SobolSequence::new(dims).sample(300),
            expected,
            "dims = {dims}"
        );
    }
}

#[test]
fn ten_dimensional_sequence_is_the_prefix_of_higher_dimensional_sequences() {
    let ten = sobol_sequence_10d(257);
    assert_eq!(ten, sobol_sequence(257, 10));
    for (low, high) in ten.iter().zip(sobol_sequence(257, 40).iter()) {
        assert_eq!(&low[..], &high[..10]);
    }
}

#[test]
fn incremental_update_agrees_with_the_closed_form_gray_code_state() {
    let points = sobol_sequence(2048, SOBOL_MAX_DIM);
    for dim in [0usize, 1, 2, 9, 20, 63, 100, SOBOL_MAX_DIM - 1] {
        let directions = sobol_directions(dim).expect("dimension is tabulated");
        for (index, point) in points.iter().enumerate() {
            assert_eq!(
                to_state(point[dim]),
                sobol_state_at(&directions, index as u32),
                "coordinate {dim}, point {index}"
            );
        }
    }
}

#[test]
fn blocks_starting_inside_the_sequence_match_the_full_sequence() {
    let directions = SobolSequence::direction_table(12);
    let full = SobolSequence::generate_block(&directions, &[], 0, 200);
    for start in [1usize, 2, 37, 64, 127] {
        let block = SobolSequence::generate_block(&directions, &[], start, 20);
        assert_eq!(block, full[start..start + 20].to_vec(), "start = {start}");
    }
}

#[cfg(target_pointer_width = "64")]
#[test]
fn sequence_repeats_after_two_to_the_32_points() {
    let directions = SobolSequence::direction_table(8);
    let period = 1usize << 32;
    let block = SobolSequence::generate_block(&directions, &[], period - 2, 4);
    let head = SobolSequence::generate_block(&directions, &[], 0, 2);
    assert_eq!(block[2..], head[..]);
    for (offset, point) in block[..2].iter().enumerate() {
        let index = (period - 2 + offset) as u32;
        for (dim, numbers) in directions.iter().enumerate() {
            assert_eq!(to_state(point[dim]), sobol_state_at(numbers, index));
        }
    }
}

#[cfg(target_pointer_width = "64")]
#[test]
fn one_dimensional_generator_repeats_after_two_to_the_32_points() {
    // Without masking the index, the Gray code of 2^32 + k differs from that of k in bit 31.
    let mut fresh = Sobol::new_1d();
    let head: Vec<f64> = (0..16).map(|_| fresh.next_1d()).collect();
    let mut wrapped = Sobol::new_1d();
    wrapped.count = 1usize << 32;
    let tail: Vec<f64> = (0..16).map(|_| wrapped.next_1d()).collect();
    assert_eq!(head, tail);
    // The Gray code of 2^32 - 1 is 2^31, so that point is the last direction number v_32 = 1.
    let mut last = Sobol::new_1d();
    last.count = (1usize << 32) - 1;
    assert_eq!(to_state(last.next_1d()), 1);
}

// -- Dimension limits ----------------------------------------------------------

#[test]
fn dimension_count_is_clamped_to_the_tabulated_range() {
    assert_eq!(SobolSequence::new(0).n_dims, 1);
    assert_eq!(SobolSequence::new(5).n_dims, 5);
    assert_eq!(SobolSequence::new(1000).n_dims, SOBOL_MAX_DIM);
    assert!(sobol_sequence(4, 0).iter().all(|p| p.len() == 1));
    assert!(
        sobol_sequence(4, 1000)
            .iter()
            .all(|p| p.len() == SOBOL_MAX_DIM)
    );
    // A directly constructed generator cannot ask for untabulated coordinates.
    let oversized = SobolSequence { n_dims: 1000 };
    assert!(oversized.sample(3).iter().all(|p| p.len() == SOBOL_MAX_DIM));
    let empty = SobolSequence { n_dims: 0 };
    assert_eq!(empty.sample(3), vec![Vec::<f64>::new(); 3]);
    assert!(sobol_sequence(0, 5).is_empty());
}

#[test]
fn scrambled_dimension_count_is_clamped_to_the_tabulated_range() {
    assert_eq!(sobol_scrambled(3, 0), vec![Vec::<f64>::new(); 3]);
    assert!(
        sobol_scrambled(3, 1000)
            .iter()
            .all(|p| p.len() == SOBOL_MAX_DIM)
    );
    assert!(sobol_scrambled(0, 4).is_empty());
}

// -- Digital shift -------------------------------------------------------------

#[test]
fn shifted_sequence_is_the_plain_sequence_xor_a_fixed_mask_per_coordinate() {
    let dims = 20;
    let plain = sobol_sequence(256, dims);
    let shifted = sobol_scrambled(256, dims);
    let masks: Vec<u32> = (0..dims).map(|dim| to_state(shifted[0][dim])).collect();
    // Point 0 of the plain sequence is the origin, so the first shifted point is the mask.
    assert!(
        masks.iter().all(|&mask| mask > 0xFFFF),
        "masks = {masks:x?}"
    );
    let distinct: HashSet<u32> = masks.iter().copied().collect();
    assert_eq!(distinct.len(), dims, "every coordinate needs its own mask");
    for (plain_point, shifted_point) in plain.iter().zip(shifted.iter()) {
        for dim in 0..dims {
            assert_eq!(
                to_state(shifted_point[dim]),
                to_state(plain_point[dim]) ^ masks[dim]
            );
        }
    }
    assert!(shifted[0].iter().all(|&x| x > 0.0 && x < 1.0));
}

#[test]
fn shifted_points_stay_in_the_half_open_unit_interval() {
    for point in sobol_scrambled(1024, 16) {
        assert!(point.iter().all(|x| (0.0..1.0).contains(x)));
    }
}
