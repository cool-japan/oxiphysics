// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`super::distance_transform_exact`] against a
//! brute-force Euclidean distance transform.

use super::{BinaryImage, distance_transform_exact};

/// Deterministic xorshift64* generator so the random images are reproducible.
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

/// Brute-force EDT: for each pixel the minimum Euclidean distance to any
/// background pixel (exact integer squared distance, then `sqrt`).
fn brute_force_edt(img: &BinaryImage) -> Vec<f64> {
    let (w, h) = (img.width, img.height);
    let background: Vec<(i64, i64)> = (0..h)
        .flat_map(|r| (0..w).map(move |c| (c, r)))
        .filter(|&(c, r)| !img.get(c, r))
        .map(|(c, r)| (c as i64, r as i64))
        .collect();
    let mut out = vec![0.0; w * h];
    for r in 0..h {
        for c in 0..w {
            let best = background
                .iter()
                .map(|&(bc, br)| {
                    let dx = bc - c as i64;
                    let dy = br - r as i64;
                    dx * dx + dy * dy
                })
                .min();
            out[r * w + c] = match best {
                Some(d2) => (d2 as f64).sqrt(),
                None => f64::INFINITY,
            };
        }
    }
    out
}

fn random_image(rng: &mut XorShift, w: usize, h: usize, fg_per_mille: u64) -> BinaryImage {
    let mut img = BinaryImage::new(w, h, false);
    for r in 0..h {
        for c in 0..w {
            img.set(c, r, rng.next_u64() % 1000 < fg_per_mille);
        }
    }
    img
}

fn assert_matches_brute_force(img: &BinaryImage) {
    let fast = distance_transform_exact(img);
    let slow = brute_force_edt(img);
    assert_eq!(fast.len(), slow.len());
    for (idx, (a, b)) in fast.iter().zip(slow.iter()).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "pixel {idx} (col {}, row {}) of {}x{} image: meijster {a} vs brute force {b}",
            idx % img.width,
            idx / img.width,
            img.width,
            img.height
        );
    }
}

#[test]
fn exact_edt_matches_brute_force_on_random_16x16_images() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    for density in [300_u64, 600, 850, 950] {
        for _ in 0..25 {
            let img = random_image(&mut rng, 16, 16, density);
            if img.foreground_count() == img.width * img.height {
                continue;
            }
            assert_matches_brute_force(&img);
        }
    }
}

#[test]
fn exact_edt_matches_brute_force_on_non_square_images() {
    let mut rng = XorShift(0x0123_4567_89AB_CDEF);
    for &(w, h) in &[
        (7_usize, 13_usize),
        (13, 7),
        (1, 9),
        (9, 1),
        (2, 17),
        (17, 2),
    ] {
        for _ in 0..20 {
            let mut img = random_image(&mut rng, w, h, 800);
            // Guarantee at least one background pixel.
            img.set(w / 2, h / 2, false);
            assert_matches_brute_force(&img);
        }
    }
}

#[test]
fn exact_edt_single_background_pixel_in_corner() {
    let mut img = BinaryImage::new(9, 6, true);
    img.set(0, 0, false);
    let dt = distance_transform_exact(&img);
    for r in 0..6 {
        for c in 0..9 {
            let expected = (((c * c + r * r) as f64).sqrt()).to_bits();
            assert_eq!(dt[r * 9 + c].to_bits(), expected, "pixel ({c},{r})");
        }
    }
    // Far corner is at distance sqrt(8² + 5²).
    assert!((dt[5 * 9 + 8] - 89.0_f64.sqrt()).abs() < 1e-12);
}

#[test]
fn exact_edt_single_background_pixel_in_far_corner() {
    let mut img = BinaryImage::new(5, 11, true);
    img.set(4, 10, false);
    assert_matches_brute_force(&img);
    let dt = distance_transform_exact(&img);
    assert!((dt[0] - (16.0_f64 + 100.0).sqrt()).abs() < 1e-12);
}

#[test]
fn exact_edt_known_square_values() {
    // 7x7 image with a 5x5 foreground block in the middle (1-pixel border of
    // background). The centre pixel is 3 pixels from the border ring.
    let mut img = BinaryImage::new(7, 7, false);
    for r in 1..6 {
        for c in 1..6 {
            img.set(c, r, true);
        }
    }
    let dt = distance_transform_exact(&img);
    assert_eq!(dt[3 * 7 + 3], 3.0);
    assert_eq!(dt[7 + 1], 1.0);
    assert_eq!(dt[2 * 7 + 2], 2.0);
    assert_eq!(dt[0], 0.0);
    assert_matches_brute_force(&img);
}

#[test]
fn exact_edt_all_background_is_zero_and_all_foreground_is_infinite() {
    let img = BinaryImage::new(4, 3, false);
    assert!(distance_transform_exact(&img).iter().all(|&d| d == 0.0));

    let img = BinaryImage::new(4, 3, true);
    let dt = distance_transform_exact(&img);
    assert_eq!(dt.len(), 12);
    assert!(dt.iter().all(|d| d.is_infinite()));

    let img = BinaryImage::new(0, 5, true);
    assert!(distance_transform_exact(&img).is_empty());
}
