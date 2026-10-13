// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crate-private helpers for random sampling.
//!
//! `rand` 0.10 panics inside `random_range` when the range is empty or its span
//! is not finite.  The simulation code samples from user-supplied bounds that
//! are legitimately degenerate in many configurations (zero velocity scatter, a
//! flat bounding box, a joint pinned at `lower == upper`, a zero-sized spawn
//! box), so every floating-point range goes through [`uniform`] instead of
//! calling `random_range` directly.

use rand::RngExt;

/// Sample a value uniformly from the half-open interval `[lo, hi)`.
///
/// A degenerate interval has no interior to sample from: `hi <= lo`, a NaN
/// bound, or a span `hi - lo` that overflows to infinity.  In that case `lo` is
/// returned instead of panicking, which is the exact answer when `lo == hi`.
pub(crate) fn uniform<R: rand::Rng + ?Sized>(rng: &mut R, lo: f64, hi: f64) -> f64 {
    if hi > lo && (hi - lo).is_finite() {
        rng.random_range(lo..hi)
    } else {
        lo
    }
}

#[cfg(test)]
mod tests {
    use super::uniform;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn samples_stay_inside_a_proper_interval() {
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..1000 {
            let x = uniform(&mut rng, -2.0, 3.0);
            assert!((-2.0..3.0).contains(&x), "x = {x}");
        }
    }

    #[test]
    fn degenerate_intervals_return_the_lower_bound() {
        let mut rng = StdRng::seed_from_u64(7);
        assert_eq!(uniform(&mut rng, 1.5, 1.5), 1.5);
        assert_eq!(uniform(&mut rng, 4.0, -4.0), 4.0);
        assert_eq!(uniform(&mut rng, -0.0, 0.0), 0.0);
        assert_eq!(uniform(&mut rng, 2.0, f64::NAN), 2.0);
        assert!(uniform(&mut rng, f64::NAN, 1.0).is_nan());
    }

    #[test]
    fn non_finite_spans_do_not_panic() {
        let mut rng = StdRng::seed_from_u64(7);
        assert_eq!(uniform(&mut rng, 0.0, f64::INFINITY), 0.0);
        assert_eq!(uniform(&mut rng, f64::NEG_INFINITY, 0.0), f64::NEG_INFINITY);
        // The span overflows even though both bounds are finite.
        assert_eq!(uniform(&mut rng, -f64::MAX, f64::MAX), -f64::MAX);
    }

    #[test]
    fn seeded_generators_reproduce_the_same_samples() {
        let mut a = StdRng::seed_from_u64(99);
        let mut b = StdRng::seed_from_u64(99);
        for _ in 0..32 {
            assert_eq!(uniform(&mut a, 0.0, 1.0), uniform(&mut b, 0.0, 1.0));
        }
    }
}
