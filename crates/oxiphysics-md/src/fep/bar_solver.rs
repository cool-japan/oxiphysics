// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact solver for Bennett's acceptance-ratio (BAR) equation in reduced units.
//!
//! With `w = β ΔU`, `ΔU = U_B − U_A`, `n_A` samples `w_i` drawn from state A,
//! `n_B` samples `w_j` drawn from state B, the Fermi function
//! `f(y) = 1 / (1 + eʸ)` and `M = ln(n_A / n_B)`, Bennett's optimal estimator
//! `x = β ΔF` is the root of
//!
//! ```text
//!   S_A(c) = Σ_i f(w_i − c)  =  Σ_j f(c − w_j) = S_B(c),     x = c + M.
//! ```
//!
//! The shift `c = x − M` makes the equation independent of the sample sizes:
//! `S_A` increases and `S_B` decreases with `c`, so `G(c) = ln S_A − ln S_B`
//! is strictly increasing with `0 < G'(c) < 2`, has exactly one root, and
//! sign-changing brackets can be written down from the extreme samples.  The
//! solver below combines that bracket with Newton steps (bisection whenever a
//! Newton step leaves the bracket), which converges for any overlap between the
//! two states, however poor, where the classical fixed-point iteration
//! `c ← c + G(c)` only contracts by `1 − G'(c)` per sweep.

use std::f64::consts::LN_2;

use super::numerics::{log_sum_exp, softplus};

/// Relative tolerance on the shift `c` (scaled by `1 + |c|`).
const ROOT_TOLERANCE: f64 = 1.0e-13;

/// Newton/bisection steps before the best available `c` is returned.
///
/// Pure bisection from a bracket of width `W` needs `log2(W / tolerance)`
/// halvings; 400 covers brackets up to `1e100` wide.
const MAX_ITERATIONS: usize = 400;

/// Outward doublings of a bracket end that failed the sign test (a safeguard:
/// the analytic bracket already satisfies it when every sample is finite).
const MAX_BRACKET_EXPANSIONS: usize = 64;

/// One side of Bennett's equation for the reduced energy differences `w` of one
/// state at shift `c`.
///
/// With `y_i = s (w_i − c)` it returns `ln Σ_i f(y_i)` and its derivative with
/// respect to `c`, `s Σ_i f(y_i) (1 − f(y_i)) / Σ_i f(y_i)`.  `s = +1` is state
/// A (`S_A`), `s = −1` state B (`S_B`).  The sum is `−∞` (and the slope 0) when
/// every term vanishes.
fn fermi_side(w: &[f64], c: f64, sign: f64) -> (f64, f64) {
    let ln_sum = log_sum_exp(w.iter().map(|&wi| -softplus(sign * (wi - c))));
    if !ln_sum.is_finite() {
        return (ln_sum, 0.0);
    }
    // 1 − f(y) = f(−y); the normalised weights p_i = f(y_i) / S keep every
    // exponent non-positive.
    let mean_complement: f64 = w
        .iter()
        .map(|&wi| {
            let y = sign * (wi - c);
            (-softplus(y) - ln_sum).exp() * (-softplus(-y)).exp()
        })
        .sum();
    (ln_sum, sign * mean_complement)
}

/// Residual `G(c) = ln S_A(c) − ln S_B(c)` of Bennett's equation and its slope
/// `G'(c) = G'_A + G'_B ∈ (0, 2)`.
fn residual(w_a: &[f64], w_b: &[f64], c: f64) -> (f64, f64) {
    let (ln_a, slope_a) = fermi_side(w_a, c, 1.0);
    let (ln_b, slope_b) = fermi_side(w_b, c, -1.0);
    (ln_a - ln_b, slope_a - slope_b)
}

/// BAR free-energy difference `β ΔF = β (F_B − F_A)` from the reduced energy
/// differences `w = β (U_B − U_A)` of samples drawn from state A (`w_a`) and
/// from state B (`w_b`).
///
/// * both slices empty: `0`;
/// * one slice empty: the one-sided exponential average (Zwanzig), in
///   log-sum-exp form: `−ln ⟨e^{−w}⟩_A` or `+ln ⟨e^{w}⟩_B`;
/// * otherwise the unique root of Bennett's equation (see the module
///   documentation), accurate to a relative `1e-13` in the shift `c`.
///
/// `±∞` entries are allowed (an energy wall: `w_A = +∞` or `w_B = −∞` carry
/// zero Fermi weight).  The result is NaN if any entry is NaN, if no entry is
/// finite, or if one state has no sample that overlaps the other at all.
pub(super) fn bar_reduced_free_energy(w_a: &[f64], w_b: &[f64]) -> f64 {
    if w_a.iter().chain(w_b).any(|w| w.is_nan()) {
        return f64::NAN;
    }
    let n_a = w_a.len() as f64;
    let n_b = w_b.len() as f64;
    match (w_a.is_empty(), w_b.is_empty()) {
        (true, true) => return 0.0,
        (false, true) => return n_a.ln() - log_sum_exp(w_a.iter().map(|&w| -w)),
        (true, false) => return log_sum_exp(w_b.iter().copied()) - n_b.ln(),
        (false, false) => {}
    }
    let ln_ratio = n_a.ln() - n_b.ln();

    let (w_min, w_max) = w_a
        .iter()
        .chain(w_b)
        .copied()
        .filter(|w| w.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), w| {
            (lo.min(w), hi.max(w))
        });
    if !w_min.is_finite() {
        return f64::NAN;
    }

    // Analytic bracket.  For c = w_min − L every state-A Fermi weight is below
    // e^{−L} and every state-B weight above 1/2, so G(c) ≤ M + ln 2 − L; for
    // c = w_max + L the weights swap and G(c) ≥ M − ln 2 + L.
    let mut lo = w_min - (ln_ratio + LN_2).max(0.0) - 1.0;
    let mut hi = w_max + (LN_2 - ln_ratio).max(0.0) + 1.0;

    // Infinite samples weaken the bounds above, so verify the signs and widen
    // the bracket outwards if needed.
    let mut step = hi - lo;
    for attempt in 0.. {
        let g = residual(w_a, w_b, lo).0;
        if g < 0.0 {
            break;
        }
        if g.is_nan() || attempt >= MAX_BRACKET_EXPANSIONS {
            return f64::NAN;
        }
        lo -= step;
        step *= 2.0;
    }
    let mut step = hi - lo;
    for attempt in 0.. {
        let g = residual(w_a, w_b, hi).0;
        if g > 0.0 {
            break;
        }
        if g.is_nan() || attempt >= MAX_BRACKET_EXPANSIONS {
            return f64::NAN;
        }
        hi += step;
        step *= 2.0;
    }

    // Safeguarded Newton iteration on the bracket [lo, hi], G(lo) < 0 < G(hi).
    let mut c = 0.5 * (lo + hi);
    for _ in 0..MAX_ITERATIONS {
        let (g, slope) = residual(w_a, w_b, c);
        if g == 0.0 {
            return c + ln_ratio;
        }
        if g < 0.0 {
            lo = c;
        } else if g > 0.0 {
            hi = c;
        } else {
            return f64::NAN;
        }
        let newton = c - g / slope;
        let next = if slope > 0.0 && newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };
        let moved = (next - c).abs();
        c = next;
        let scale = ROOT_TOLERANCE * (1.0 + c.abs());
        if moved <= scale || hi - lo <= scale {
            break;
        }
    }
    c + ln_ratio
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `w = d` for every sample of both states.
    fn constant(d: f64, n_a: usize, n_b: usize) -> (Vec<f64>, Vec<f64>) {
        (vec![d; n_a], vec![d; n_b])
    }

    #[test]
    fn constant_shift_is_recovered_for_every_sample_size_ratio() {
        // Equal-size sets reduce to c = d, i.e. x = d; with n_A != n_B the root is
        // c = d - M and only the +M term of x = c + M gives d back.
        for (n_a, n_b) in [(1, 1), (7, 7), (30, 10), (10, 30), (1, 500), (500, 1)] {
            for d in [-4.0, 0.0, 2.0, 11.5] {
                let (w_a, w_b) = constant(d, n_a, n_b);
                let x = bar_reduced_free_energy(&w_a, &w_b);
                assert!(
                    (x - d).abs() < 1e-9 * (1.0 + d.abs()),
                    "n_A = {n_a}, n_B = {n_b}, d = {d}: x = {x}"
                );
            }
        }
    }

    #[test]
    fn residual_is_strictly_increasing_with_slope_between_zero_and_two() {
        let w_a = [0.3, 1.1, 2.4, -0.6, 3.0];
        let w_b = [-0.9, 0.4, 1.7];
        let mut previous = f64::NEG_INFINITY;
        for step in -40..=40 {
            let c = f64::from(step) * 0.25;
            let (g, slope) = residual(&w_a, &w_b, c);
            assert!(g > previous, "G is not increasing at c = {c}");
            assert!(slope > 0.0 && slope < 2.0, "G'({c}) = {slope}");
            // The analytic slope matches a central difference.
            let h = 1e-6;
            let numeric =
                (residual(&w_a, &w_b, c + h).0 - residual(&w_a, &w_b, c - h).0) / (2.0 * h);
            assert!(
                (slope - numeric).abs() < 1e-7 * (1.0 + numeric.abs()),
                "slope {slope} vs finite difference {numeric} at c = {c}"
            );
            previous = g;
        }
    }

    #[test]
    fn root_satisfies_bennetts_equation_for_irregular_samples() {
        let w_a = [0.3, 1.1, 2.4, -0.6, 3.0, 1.9, 0.8];
        let w_b = [-0.9, 0.4, 1.7, -2.2];
        let x = bar_reduced_free_energy(&w_a, &w_b);
        let c = x - (w_a.len() as f64 / w_b.len() as f64).ln();
        let sum_a: f64 = w_a.iter().map(|&w| 1.0 / (1.0 + (w - c).exp())).sum();
        let sum_b: f64 = w_b.iter().map(|&w| 1.0 / (1.0 + (c - w).exp())).sum();
        assert!(
            (sum_a - sum_b).abs() < 1e-11,
            "S_A = {sum_a}, S_B = {sum_b} at x = {x}"
        );
    }

    #[test]
    fn extreme_reduced_energies_do_not_overflow() {
        // Joule-valued energies (beta ~ 2.4e20 per J) give |w| ~ 1e20.
        let w_a = [2.4e20, -2.4e20, 4.8e20, -4.8e20];
        let w_b = [-2.4e20, 2.4e20, -4.8e20, 4.8e20];
        let x = bar_reduced_free_energy(&w_a, &w_b);
        assert!(x.is_finite() && x.abs() < 1e6, "x = {x}");

        // Widely separated states: the estimate stays finite and at the right scale.
        let far_a = vec![1.0e6; 20];
        let far_b = vec![1.0e6; 20];
        let far = bar_reduced_free_energy(&far_a, &far_b);
        assert!((far - 1.0e6).abs() < 1e-6, "x = {far}");
    }

    #[test]
    fn one_sided_inputs_fall_back_to_the_zwanzig_formulas() {
        assert_eq!(bar_reduced_free_energy(&[], &[]), 0.0);
        // A only: -ln <e^{-w}>; constant w = d gives d exactly.
        let only_a = bar_reduced_free_energy(&[2.5, 2.5, 2.5], &[]);
        assert!((only_a - 2.5).abs() < 1e-14, "{only_a}");
        // B only: +ln <e^{w}>.
        let only_b = bar_reduced_free_energy(&[], &[-1.5, -1.5]);
        assert!((only_b + 1.5).abs() < 1e-14, "{only_b}");
        // A genuine average: -ln((e^0 + e^-4) / 2).
        let mixed = bar_reduced_free_energy(&[0.0, 4.0], &[]);
        let expected = -(0.5 * (1.0 + (-4.0_f64).exp())).ln();
        assert!((mixed - expected).abs() < 1e-14, "{mixed} vs {expected}");
        // Huge energies must not overflow the average either.
        let huge = bar_reduced_free_energy(&[2.0e20, 2.0e20], &[]);
        assert!((huge - 2.0e20).abs() < 1e7, "{huge}");
    }

    #[test]
    fn invalid_inputs_give_nan() {
        assert!(bar_reduced_free_energy(&[f64::NAN, 1.0], &[0.0]).is_nan());
        assert!(bar_reduced_free_energy(&[1.0], &[f64::NAN]).is_nan());
        assert!(bar_reduced_free_energy(&[f64::NAN], &[]).is_nan());
        // No finite entry at all: no information about the shift.
        assert!(bar_reduced_free_energy(&[f64::NEG_INFINITY], &[f64::INFINITY]).is_nan());
        // State A entirely inside an energy wall of state B: no overlap.
        assert!(bar_reduced_free_energy(&[f64::INFINITY], &[1.0, 2.0]).is_nan());
    }

    #[test]
    fn energy_walls_carry_no_fermi_weight() {
        // A samples that hit a wall of state B (w = +inf) only dilute S_A; the
        // estimate equals the one obtained when those samples are counted but
        // given w = 1e30 (weight exactly zero in double precision).
        let with_wall = bar_reduced_free_energy(&[0.5, 1.0, f64::INFINITY], &[-0.5, 0.2]);
        let with_big = bar_reduced_free_energy(&[0.5, 1.0, 1.0e30], &[-0.5, 0.2]);
        assert!(with_wall.is_finite());
        assert!(
            (with_wall - with_big).abs() < 1e-9,
            "{with_wall} vs {with_big}"
        );
    }
}
