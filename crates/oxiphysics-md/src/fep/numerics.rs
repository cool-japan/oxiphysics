// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Overflow-safe building blocks shared by the BAR and MBAR estimators of
//! [`super`].
//!
//! Reduced energies `β ΔU` of molecular systems range over many orders of
//! magnitude (a Joule-valued `ΔU` gives `β ΔU ~ 10²⁰`), so every exponential
//! sum of the free-energy estimators is evaluated in logarithmic form.

/// Softplus function `ln(1 + eʸ)`.
///
/// Finite for every finite `y` (no overflow for large positive `y`, no
/// cancellation for large negative `y`) and equal to `+∞` / `0` at `y = ±∞`.
/// The Fermi function used by BAR is `f(y) = 1 / (1 + eʸ) = exp(−softplus(y))`.
pub(super) fn softplus(y: f64) -> f64 {
    y.max(0.0) + (-y.abs()).exp().ln_1p()
}

/// `ln Σ_i exp(x_i)` over the items of `values`, evaluated without overflow.
///
/// Returns `−∞` for an empty input or one that only holds `−∞`, `+∞` if an
/// item is `+∞`, and NaN if any item is NaN.  The iterator is traversed twice
/// (maximum, then the shifted sum), hence the `Clone` bound.
pub(super) fn log_sum_exp<I>(values: I) -> f64
where
    I: Iterator<Item = f64> + Clone,
{
    let mut max = f64::NEG_INFINITY;
    for value in values.clone() {
        if value.is_nan() {
            return f64::NAN;
        }
        max = max.max(value);
    }
    if !max.is_finite() {
        return max;
    }
    max + values.map(|v| (v - max).exp()).sum::<f64>().ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn softplus_matches_the_definition_and_never_overflows() {
        for y in [-30.0_f64, -2.5, -1e-3, 0.0, 1e-3, 0.7, 3.0, 25.0] {
            let direct = (1.0 + y.exp()).ln();
            assert!(
                (softplus(y) - direct).abs() < 1e-12 * (1.0 + direct.abs()),
                "softplus({y}) = {}, direct {direct}",
                softplus(y)
            );
        }
        assert!((softplus(0.0) - std::f64::consts::LN_2).abs() < 1e-15);
        // exp(1000) overflows, softplus(1000) = 1000 and softplus(-1000) = 0 to f64.
        assert_eq!(softplus(1000.0), 1000.0);
        assert!(softplus(-1000.0) >= 0.0 && softplus(-1000.0) < 1e-300);
        assert_eq!(softplus(f64::INFINITY), f64::INFINITY);
        assert_eq!(softplus(f64::NEG_INFINITY), 0.0);
        assert!(softplus(f64::NAN).is_nan());
        // Huge magnitudes appear for Joule-valued energies (beta ~ 2e20 per J).
        assert_eq!(softplus(2.4e20), 2.4e20);
    }

    #[test]
    fn log_sum_exp_is_stable_and_handles_the_degenerate_inputs() {
        let values = [0.5_f64, -1.25, 2.0];
        let direct = values.iter().map(|v| v.exp()).sum::<f64>().ln();
        assert!((log_sum_exp(values.iter().copied()) - direct).abs() < 1e-14);
        // Terms around 1000 would overflow a naive exp; the shift keeps them finite.
        let big = log_sum_exp([1000.0_f64, 1000.0].iter().copied());
        assert!(
            (big - (1000.0 + std::f64::consts::LN_2)).abs() < 1e-12,
            "{big}"
        );
        let tiny = log_sum_exp([-1000.0_f64, -1001.0].iter().copied());
        assert!(
            tiny.is_finite() && tiny < -999.0 && tiny > -1001.0,
            "{tiny}"
        );
        assert_eq!(log_sum_exp(std::iter::empty()), f64::NEG_INFINITY);
        assert_eq!(
            log_sum_exp([f64::NEG_INFINITY, f64::NEG_INFINITY].iter().copied()),
            f64::NEG_INFINITY
        );
        // A -inf term simply drops out.
        let with_zero_weight = log_sum_exp([f64::NEG_INFINITY, 1.0].iter().copied());
        assert!((with_zero_weight - 1.0).abs() < 1e-15);
        assert_eq!(
            log_sum_exp([1.0, f64::INFINITY].iter().copied()),
            f64::INFINITY
        );
        assert!(log_sum_exp([1.0, f64::NAN].iter().copied()).is_nan());
        assert!(log_sum_exp([f64::NAN, f64::NAN].iter().copied()).is_nan());
    }
}
