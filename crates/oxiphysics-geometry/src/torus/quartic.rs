// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real roots of a quartic polynomial, used by the ray/torus intersection.
//!
//! The quartic is solved with Ferrari's method and every root is polished
//! with Newton's method on the original polynomial.
//!
//! With `t = u - a/4` the monic quartic `t^4 + a t^3 + b t^2 + c t + d` becomes
//! `u^4 + p u^2 + q u + r`.  Adding `m` to the square `(u^2 + p/2)^2` turns the
//! equation into
//!
//! ```text
//! (u^2 + p/2 + m)^2 = 2m u^2 - q u + (m^2 + p m + p^2/4 - r)
//! ```
//!
//! and the right-hand side is a perfect square, `(s u - q/(2 s))^2` with
//! `s = sqrt(2m)`, exactly when `m` is a root of the resolvent cubic
//!
//! ```text
//! 8 m^3 + 8 p m^2 + (2 p^2 - 8 r) m - q^2 = 0.
//! ```
//!
//! Taking the square root of both sides gives the two real quadratics
//!
//! ```text
//! u^2 + s u + (p/2 + m - q/(2 s)) = 0   and   u^2 - s u + (p/2 + m + q/(2 s)) = 0
//! ```
//!
//! whose roots are the four roots of the quartic.  Only the largest real root
//! of the resolvent is used: every root of the resolvent gives the same four
//! quartic roots, so using all of them would count each root several times.
//! The resolvent satisfies `g(0) = -q^2 <= 0`, so its largest root is never
//! negative and `m = 0` only when `q = 0`, where the quartic is biquadratic:
//! `u^2 = -p/2 +- sqrt(p^2/4 - r)`.

/// Upper limit of Newton iterations spent on one root.
const POLISH_ITERATIONS: usize = 12;

/// Discriminants down to this fraction (of the terms that formed them) below
/// zero count as zero.  A tangent line gives a double root whose discriminant
/// is zero in exact arithmetic but may come out as a rounding residue.
const DISCRIMINANT_TOLERANCE: f64 = 1e-12;

/// A resolvent root at most this fraction of the polynomial scale is treated as
/// zero: the quartic is then biquadratic up to a perturbation that the Newton
/// polish removes.
const RESOLVENT_ZERO_TOLERANCE: f64 = 1e-11;

/// A polished root whose relative residual is above this is discarded.
const RESIDUAL_TOLERANCE: f64 = 1e-6;

/// Monic quartic `t^4 + a t^3 + b t^2 + c t + d` and its derivative at `t`.
#[inline]
fn quartic_and_slope(coeffs: [f64; 4], t: f64) -> (f64, f64) {
    let [a, b, c, d] = coeffs;
    let value = (((t + a) * t + b) * t + c) * t + d;
    let slope = ((4.0 * t + 3.0 * a) * t + 2.0 * b) * t + c;
    (value, slope)
}

/// `|t|^4 + |a||t|^3 + |b| t^2 + |c||t| + |d|`: the size of the terms of the
/// quartic at `t`, which bounds the rounding error of its evaluation.
#[inline]
fn quartic_term_size(coeffs: [f64; 4], t: f64) -> f64 {
    let [a, b, c, d] = coeffs;
    let x = t.abs();
    (((x + a.abs()) * x + b.abs()) * x + c.abs()) * x + d.abs()
}

/// Newton iteration on the monic quartic starting at `t`.
///
/// A step is only accepted when it reduces the residual (halving it up to a
/// few times otherwise), so a root with vanishing slope, as for a tangent line,
/// is never made worse.
fn polish_root(coeffs: [f64; 4], start: f64) -> f64 {
    let mut t = start;
    let (mut value, mut slope) = quartic_and_slope(coeffs, t);
    for _ in 0..POLISH_ITERATIONS {
        if value == 0.0 || slope == 0.0 {
            break;
        }
        let step = value / slope;
        if !step.is_finite() {
            break;
        }
        let mut damping = 1.0;
        let mut accepted = false;
        for _ in 0..4 {
            let candidate = t - damping * step;
            let (candidate_value, candidate_slope) = quartic_and_slope(coeffs, candidate);
            if candidate_value.abs() < value.abs() {
                t = candidate;
                value = candidate_value;
                slope = candidate_slope;
                accepted = true;
                break;
            }
            damping *= 0.5;
        }
        if !accepted || (damping * step).abs() <= 4.0 * f64::EPSILON * t.abs() {
            break;
        }
    }
    t
}

/// Largest real root of the monic cubic `g(m) = m^3 + b m^2 + c m + d` with
/// `d <= 0`; it is never negative because `g(0) = d <= 0`.
///
/// The root is bracketed in an interval where `g` increases from a value
/// `<= 0` to a value `>= 0` (using the critical points of `g`), then found with
/// a safeguarded Newton iteration that falls back to bisection.
fn largest_cubic_root(b: f64, c: f64, d: f64) -> f64 {
    let g = |m: f64| ((m + b) * m + c) * m + d;
    let slope = |m: f64| (3.0 * m + 2.0 * b) * m + c;

    // Fujiwara's bound on the size of the roots of a monic cubic.
    let bound = 2.0 * b.abs().max(c.abs().sqrt()).max((0.5 * d.abs()).cbrt());
    let mut lo = 0.0_f64;
    let mut hi = bound;

    let critical_discriminant = b * b - 3.0 * c;
    if critical_discriminant > 0.0 {
        let root = critical_discriminant.sqrt();
        let local_min = (-b + root) / 3.0;
        let local_max = (-b - root) / 3.0;
        if g(local_min) <= 0.0 {
            // g increases on [local_min, inf): the largest root is there.
            lo = lo.max(local_min);
        } else {
            // The local minimum is above the axis: the only real root is left
            // of the local maximum, where g increases.
            hi = local_max.max(lo);
        }
    }
    // Guard against a bound that rounding made slightly too small.
    for _ in 0..64 {
        if g(hi) >= 0.0 {
            break;
        }
        hi = if hi > 0.0 { hi * 2.0 } else { 1.0 };
    }

    if g(lo) >= 0.0 {
        return lo;
    }
    let mut m = 0.5 * (lo + hi);
    for _ in 0..200 {
        let value = g(m);
        if value == 0.0 {
            break;
        }
        if value < 0.0 {
            lo = m;
        } else {
            hi = m;
        }
        let derivative = slope(m);
        let newton = m - value / derivative;
        let next = if derivative > 0.0 && newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };
        let converged = (next - m).abs() <= 4.0 * f64::EPSILON * (next.abs() + 1e-3 * bound);
        m = next;
        if converged {
            break;
        }
    }
    m
}

/// Real roots `u` of the quadratic `u^2 + b u + c = 0`, with a discriminant
/// slightly below zero accepted as a double root.
fn quadratic_roots(b: f64, c: f64, out: &mut Vec<f64>) {
    let discriminant = b * b - 4.0 * c;
    if discriminant < -DISCRIMINANT_TOLERANCE * (b * b + 4.0 * c.abs()) {
        return;
    }
    push_quadratic_roots(b, c, discriminant.max(0.0), out);
}

/// Roots of `u^2 + b u + c = 0` for a given non-negative discriminant, without
/// the cancellation of the textbook formula.
fn push_quadratic_roots(b: f64, c: f64, discriminant: f64, out: &mut Vec<f64>) {
    let root = discriminant.sqrt();
    let sign = if b < 0.0 { -1.0 } else { 1.0 };
    let q = -0.5 * (b + sign * root);
    if q == 0.0 {
        // b = 0 and the discriminant is zero, hence c = 0.
        out.push(0.0);
        out.push(0.0);
    } else {
        out.push(q);
        out.push(c / q);
    }
}

/// Real roots of `c4 t^4 + c3 t^3 + c2 t^2 + c1 t + c0 = 0`, in ascending
/// order, counted with multiplicity as far as the rounding allows.
///
/// Fewer than four real roots are returned as `NaN` padding at the end, so the
/// result is `[r0, r1, ..., NaN]`.  All entries are `NaN` when `c4` is zero or
/// a coefficient is not finite.
///
/// The coefficients are rescaled so that the tolerances are relative, and the
/// roots are polished on the original polynomial.  A double root (a tangent
/// line in the torus intersection) is returned twice and is accurate to about
/// the square root of the machine precision.
pub(super) fn solve_quartic(c4: f64, c3: f64, c2: f64, c1: f64, c0: f64) -> [f64; 4] {
    let mut result = [f64::NAN; 4];
    let finite = [c4, c3, c2, c1, c0].iter().all(|c| c.is_finite());
    if !finite || c4 == 0.0 {
        return result;
    }

    // Monic form, with the variable scaled by the size of the roots.
    let (a0, b0, c_0, d0) = (c3 / c4, c2 / c4, c1 / c4, c0 / c4);
    let sigma = a0
        .abs()
        .max(b0.abs().sqrt())
        .max(c_0.abs().cbrt())
        .max(d0.abs().sqrt().sqrt());
    let sigma = if sigma > 0.0 && sigma.is_finite() {
        sigma
    } else {
        1.0
    };
    let a = a0 / sigma;
    let b = b0 / (sigma * sigma);
    let c = c_0 / (sigma * sigma * sigma);
    let d = d0 / (sigma * sigma * sigma * sigma);
    let coeffs = [a, b, c, d];

    // Depressed quartic u^4 + p u^2 + q u + r with t = u - a/4.
    let shift = a / 4.0;
    let p = b - 3.0 * a * a / 8.0;
    let q = c - a * b / 2.0 + a * a * a / 8.0;
    let r = d - a * c / 4.0 + a * a * b / 16.0 - 3.0 * a * a * a * a / 256.0;

    let mut roots: Vec<f64> = Vec::with_capacity(4);
    let scale = p.abs().max(r.abs().sqrt()).max(q.abs().powf(2.0 / 3.0));
    let m = if q == 0.0 {
        0.0
    } else {
        largest_cubic_root(p, 0.25 * p * p - r, -0.125 * q * q)
    };
    if m > RESOLVENT_ZERO_TOLERANCE * scale && m > 0.0 {
        let s = (2.0 * m).sqrt();
        for sign in [1.0_f64, -1.0] {
            // u^2 + sign s u + (p/2 + m - sign q/(2 s)) = 0
            let constant = 0.5 * p + m - sign * q / (2.0 * s);
            let discriminant = 2.0 * (sign * q / s - m - p);
            if discriminant < -DISCRIMINANT_TOLERANCE * (2.0 * m + 4.0 * constant.abs()) {
                continue;
            }
            // Stable form of the roots of u^2 + (sign s) u + constant.
            let root = discriminant.max(0.0).sqrt();
            let first = -sign * 0.5 * (s + root);
            roots.push(first);
            roots.push(if first == 0.0 { 0.0 } else { constant / first });
        }
    } else {
        // Biquadratic: u^2 = z with z^2 + p z + r = 0.
        let mut squares: Vec<f64> = Vec::with_capacity(2);
        quadratic_roots(p, r, &mut squares);
        let tolerance = DISCRIMINANT_TOLERANCE * (p.abs() + r.abs().sqrt());
        for z in squares {
            if z >= -tolerance {
                let u = z.max(0.0).sqrt();
                roots.push(u);
                roots.push(-u);
            }
        }
    }

    let mut polished: Vec<f64> = roots
        .into_iter()
        .map(|u| polish_root(coeffs, u - shift))
        .filter(|&t| {
            let (value, _) = quartic_and_slope(coeffs, t);
            t.is_finite() && value.abs() <= RESIDUAL_TOLERANCE * quartic_term_size(coeffs, t)
        })
        .collect();
    polished.sort_by(f64::total_cmp);
    for (slot, t) in result.iter_mut().zip(polished) {
        *slot = t * sigma;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Coefficients `[1, c3, c2, c1, c0]` of `prod (t - root)`.
    fn expand(roots: &[f64]) -> [f64; 5] {
        let mut coeffs = [0.0; 5];
        coeffs[0] = 1.0;
        for (k, &root) in roots.iter().enumerate() {
            for j in (1..=k + 1).rev() {
                coeffs[j] -= root * coeffs[j - 1];
            }
        }
        coeffs
    }

    fn solve(coeffs: [f64; 5]) -> Vec<f64> {
        solve_quartic(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4])
            .into_iter()
            .filter(|t| !t.is_nan())
            .collect()
    }

    fn assert_roots(actual: &[f64], expected: &[f64], tolerance: f64) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "root count: got {actual:?}, expected {expected:?}"
        );
        for (a, e) in actual.iter().zip(expected) {
            assert!(
                (a - e).abs() <= tolerance * (1.0 + e.abs()),
                "got {actual:?}, expected {expected:?}"
            );
        }
    }

    /// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
    fn lcg(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    #[test]
    fn four_distinct_roots_with_a_biquadratic_structure() {
        let roots = [1.0, 2.0, 3.0, 4.0];
        assert_roots(&solve(expand(&roots)), &roots, 1e-12);
    }

    #[test]
    fn four_distinct_roots_with_a_linear_term() {
        let roots = [-1.0, 0.5, 2.0, 3.0];
        assert_roots(&solve(expand(&roots)), &roots, 1e-12);
    }

    #[test]
    fn two_real_roots_and_a_complex_pair() {
        // (t - 1)(t + 2)(t^2 + t + 1)
        let c = [1.0, 2.0, 0.0, -1.0, -2.0];
        // expand: (t^2 + t - 2)(t^2 + t + 1) = t^4 + 2t^3 + 0t^2 - t - 2
        assert_roots(&solve(c), &[-2.0, 1.0], 1e-12);
    }

    #[test]
    fn no_real_roots() {
        // t^4 + 1
        assert!(solve([1.0, 0.0, 0.0, 0.0, 1.0]).is_empty());
        // (t^2 + 1)(t^2 + 4)
        assert!(solve([1.0, 0.0, 5.0, 0.0, 4.0]).is_empty());
    }

    #[test]
    fn biquadratic_with_two_real_roots() {
        // t^4 - 1 = (t - 1)(t + 1)(t^2 + 1)
        assert_roots(&solve([1.0, 0.0, 0.0, 0.0, -1.0]), &[-1.0, 1.0], 1e-14);
        // (t^2 - 2)(t^2 - 3)
        let s2 = 2.0_f64.sqrt();
        let s3 = 3.0_f64.sqrt();
        assert_roots(
            &solve([1.0, 0.0, -5.0, 0.0, 6.0]),
            &[-s3, -s2, s2, s3],
            1e-13,
        );
    }

    #[test]
    fn quartic_with_a_zero_constant_term() {
        // t (t - 1)(t - 2)(t + 3)
        let roots = [-3.0, 0.0, 1.0, 2.0];
        assert_roots(&solve(expand(&roots)), &roots, 1e-12);
    }

    #[test]
    fn double_root_is_returned_twice() {
        // (t - 1)^2 (t - 3)(t + 2)
        let c = expand(&[1.0, 1.0, 3.0, -2.0]);
        let roots = solve(c);
        assert_eq!(roots.len(), 4, "got {roots:?}");
        assert!((roots[0] + 2.0).abs() < 1e-9);
        assert!((roots[1] - 1.0).abs() < 1e-6 && (roots[2] - 1.0).abs() < 1e-6);
        assert!((roots[3] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn quadruple_root() {
        let roots = solve(expand(&[2.0, 2.0, 2.0, 2.0]));
        assert_eq!(roots.len(), 4, "got {roots:?}");
        assert!(
            roots.iter().all(|r| (r - 2.0).abs() < 1e-3),
            "got {roots:?}"
        );
        let zero = solve([1.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(zero.len(), 4);
        assert!(zero.iter().all(|r| r.abs() < 1e-12));
    }

    #[test]
    fn scale_invariance() {
        let base = [-2.0, 0.3, 1.0, 5.0];
        for &scale in &[1e-6, 1e-3, 1.0, 1e3, 1e6] {
            let scaled: Vec<f64> = base.iter().map(|r| r * scale).collect();
            let found: Vec<f64> = solve(expand(&scaled)).iter().map(|t| t / scale).collect();
            assert_roots(&found, &base, 1e-9);
        }
    }

    #[test]
    fn random_real_roots_are_recovered() {
        let mut next = lcg(0xdec0_ded0);
        let mut checked = 0;
        for _ in 0..500 {
            let mut roots: Vec<f64> = (0..4).map(|_| -10.0 + 20.0 * next()).collect();
            roots.sort_by(f64::total_cmp);
            // Keep the roots well separated so that the test measures the
            // solver, not the conditioning of the problem.
            if roots.windows(2).any(|w| w[1] - w[0] < 0.2) {
                continue;
            }
            let found = solve(expand(&roots));
            assert_roots(&found, &roots, 1e-8);
            checked += 1;
        }
        assert!(checked > 100, "only {checked} well separated samples");
    }

    #[test]
    fn random_two_real_roots_and_a_complex_pair() {
        let mut next = lcg(0x1234_5678);
        let mut checked = 0;
        for _ in 0..500 {
            let x1 = -5.0 + 10.0 * next();
            let x2 = -5.0 + 10.0 * next();
            if (x1 - x2).abs() < 0.3 {
                continue;
            }
            // t^2 + u t + v with a negative discriminant.
            let u = -4.0 + 8.0 * next();
            let v = u * u / 4.0 + 0.5 + 3.0 * next();
            let linear = expand(&[x1, x2]);
            // (t^2 + l1 t + l0)(t^2 + u t + v)
            let (l1, l0) = (linear[1], linear[2]);
            let c = [1.0, l1 + u, l0 + l1 * u + v, l0 * u + l1 * v, l0 * v];
            let mut expected = [x1, x2];
            expected.sort_by(f64::total_cmp);
            assert_roots(&solve(c), &expected, 1e-8);
            checked += 1;
        }
        assert!(checked > 300);
    }

    #[test]
    fn degenerate_input_returns_no_roots() {
        assert!(solve_quartic(0.0, 1.0, 2.0, 3.0, 4.0)[0].is_nan());
        assert!(solve_quartic(1.0, f64::NAN, 0.0, 0.0, 0.0)[0].is_nan());
        assert!(solve_quartic(1.0, f64::INFINITY, 0.0, 0.0, 0.0)[0].is_nan());
    }

    #[test]
    fn largest_cubic_root_matches_known_factorisations() {
        // (m - 1)(m - 2)(m - 3) = m^3 - 6 m^2 + 11 m - 6
        assert!((largest_cubic_root(-6.0, 11.0, -6.0) - 3.0).abs() < 1e-12);
        // m^3 - m = m (m - 1)(m + 1): d = 0, largest root 1.
        assert!((largest_cubic_root(0.0, -1.0, 0.0) - 1.0).abs() < 1e-12);
        // m^3 + m^2 = m^2 (m + 1): largest root 0 (a double root).
        assert!(largest_cubic_root(1.0, 0.0, 0.0).abs() < 1e-12);
        // m^3 - 8: single real root 2.
        assert!((largest_cubic_root(0.0, 0.0, -8.0) - 2.0).abs() < 1e-12);
        // m^3: triple root at 0.
        assert_eq!(largest_cubic_root(0.0, 0.0, 0.0), 0.0);
        // (m - 4)(m^2 + 1) = m^3 - 4 m^2 + m - 4: the single real root is 4.
        assert!((largest_cubic_root(-4.0, 1.0, -4.0) - 4.0).abs() < 1e-12);
    }
}
