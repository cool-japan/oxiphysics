// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Roots of real polynomials with verified accuracy.
//!
//! Degrees one and two use closed forms.  Higher degrees use the eigenvalues of
//! the companion matrix of the monic polynomial, rescaled by the geometric mean
//! of the root moduli (real Schur form, bounded number of QR sweeps), followed
//! by safeguarded Newton iterations.  Every root is then checked against its
//! scaled backward error `|p(r)| / sum_k |c_k| |r|^(n-k)`; a root that cannot be
//! verified is reported as an error instead of being returned as a fabricated
//! number.

use crate::error::PhysicsError;
use crate::numerical_linear_algebra::general_eigen::real_matrix_eigenvalues;

/// A complex number `(re, im)`.
type Root = (f64, f64);

/// Largest scaled backward error accepted for a computed root.
const MAX_BACKWARD_ERROR: f64 = 1e-8;

/// Newton iterations used to polish a root of the companion eigenproblem.
const POLISH_ITERATIONS: usize = 6;

/// Fraction of the distance to the nearest other root that one Newton step may
/// cover.  It keeps a polished root from jumping into the basin of its neighbour
/// when roots are clustered.
const POLISH_STEP_FRACTION: f64 = 0.1;

fn invalid_coefficients(reason: &str) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: "coefficients".to_string(),
        reason: reason.to_string(),
    }
}

/// Roots of the polynomial `c[0] x^n + c[1] x^(n-1) + ... + c[n]` as `(re, im)`
/// pairs sorted by real part and then imaginary part.
///
/// Leading zero coefficients are dropped and trailing zeros become roots at the
/// origin; a nonzero constant has no roots.  Complex roots come as exactly
/// conjugate pairs.
///
/// # Errors
///
/// * [`PhysicsError::InvalidInput`] if `coefficients` is empty, identically
///   zero, holds a non-finite value, or its leading coefficient is so small that
///   the monic coefficients overflow;
/// * [`PhysicsError::NumericalDivergence`] if the companion eigenvalue iteration
///   fails or a computed root has a backward error above `1e-8`.
pub(super) fn polynomial_roots(coefficients: &[f64]) -> Result<Vec<Root>, PhysicsError> {
    if coefficients.iter().any(|c| !c.is_finite()) {
        return Err(invalid_coefficients("contains a non-finite value"));
    }
    let Some(first) = coefficients.iter().position(|&c| c != 0.0) else {
        return Err(invalid_coefficients(if coefficients.is_empty() {
            "is empty"
        } else {
            "is the zero polynomial, which has no well-defined roots"
        }));
    };
    let trimmed = &coefficients[first..];
    let Some(last) = trimmed.iter().rposition(|&c| c != 0.0) else {
        return Err(invalid_coefficients("is the zero polynomial"));
    };
    let origin_roots = trimmed.len() - 1 - last;
    let core = &trimmed[..=last];
    let monic: Vec<f64> = core.iter().map(|c| c / core[0]).collect();
    if monic.iter().any(|c| !c.is_finite()) {
        return Err(invalid_coefficients(
            "has a leading coefficient too small relative to the others",
        ));
    }

    let mut roots: Vec<Root> = Vec::with_capacity(trimmed.len() - 1);
    match monic.len() - 1 {
        0 => {}
        1 => roots.push((-monic[1], 0.0)),
        2 => roots.extend(quadratic_roots(monic[1], monic[2])),
        _ => {
            let mut found = companion_roots(&monic)?;
            polish_roots(&monic, &mut found);
            roots.extend(found);
        }
    }
    for &root in &roots {
        let error = backward_error(&monic, root);
        let verified = error <= MAX_BACKWARD_ERROR;
        if !verified {
            return Err(PhysicsError::NumericalDivergence {
                message: format!(
                    "polynomial root {:e}{:+e}i could not be verified (backward error {error:e})",
                    root.0, root.1
                ),
            });
        }
    }
    roots.extend(std::iter::repeat_n((0.0, 0.0), origin_roots));
    // `x + 0.0` turns a negative zero into a positive one.
    for root in &mut roots {
        *root = (root.0 + 0.0, root.1 + 0.0);
    }
    roots.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    Ok(roots)
}

/// Roots of `x^2 + b x + c`, without cancellation.
fn quadratic_roots(b: f64, c: f64) -> [Root; 2] {
    let discriminant = b * b - 4.0 * c;
    if discriminant < 0.0 {
        let imaginary = 0.5 * (-discriminant).sqrt();
        return [(-0.5 * b, -imaginary), (-0.5 * b, imaginary)];
    }
    let root = discriminant.sqrt();
    let q = -0.5 * (b + root.copysign(b));
    if q == 0.0 {
        return [(0.0, 0.0), (0.0, 0.0)];
    }
    [(q, 0.0), (c / q, 0.0)]
}

/// Eigenvalues of the companion matrix of the monic polynomial `monic`
/// (`monic[0] == 1`, degree at least 3, `monic.last() != 0`).
///
/// With `x = s t` the scale `s = |a_n|^(1/n)` is the geometric mean of the root
/// moduli, so the scaled roots `t` are of order one whatever the time scales of
/// the polynomial.
fn companion_roots(monic: &[f64]) -> Result<Vec<Root>, PhysicsError> {
    let degree = monic.len() - 1;
    let scale = monic[degree].abs().powf(1.0 / degree as f64);
    let scale = if scale.is_finite() && scale > 0.0 {
        scale
    } else {
        1.0
    };
    let mut matrix = vec![0.0_f64; degree * degree];
    let mut power = 1.0_f64;
    for (k, coefficient) in monic.iter().enumerate().skip(1) {
        power *= scale;
        matrix[k - 1] = -coefficient / power;
    }
    for i in 1..degree {
        matrix[i * degree + (i - 1)] = 1.0;
    }
    let eigenvalues = real_matrix_eigenvalues(&matrix, degree).ok_or_else(|| {
        PhysicsError::NumericalDivergence {
            message: "companion matrix eigenvalue iteration did not converge".to_string(),
        }
    })?;
    Ok(eigenvalues
        .into_iter()
        .map(|(re, im)| (re * scale, im * scale))
        .collect())
}

/// `a / b` for complex numbers (Smith's algorithm); `None` if `b` is zero.
fn complex_div(a: Root, b: Root) -> Option<Root> {
    let (ar, ai) = a;
    let (br, bi) = b;
    if br == 0.0 && bi == 0.0 {
        return None;
    }
    if br.abs() >= bi.abs() {
        let r = bi / br;
        let den = br + bi * r;
        Some(((ar + ai * r) / den, (ai - ar * r) / den))
    } else {
        let r = br / bi;
        let den = br * r + bi;
        Some(((ar * r + ai) / den, (ai * r - ar) / den))
    }
}

/// `(p(z), p'(z))` for the polynomial with descending coefficients `coefficients`.
fn horner_with_derivative(coefficients: &[f64], z: Root) -> (Root, Root) {
    let mut value = (0.0_f64, 0.0_f64);
    let mut derivative = (0.0_f64, 0.0_f64);
    for &c in coefficients {
        derivative = (
            derivative.0 * z.0 - derivative.1 * z.1 + value.0,
            derivative.0 * z.1 + derivative.1 * z.0 + value.1,
        );
        value = (
            value.0 * z.0 - value.1 * z.1 + c,
            value.0 * z.1 + value.1 * z.0,
        );
    }
    (value, derivative)
}

/// Scaled backward error `|p(z)| / sum_k |a_k| |z|^(n-k)` of the candidate root
/// `z` of the monic polynomial `monic`: the relative size of the coefficient
/// perturbation for which `z` is an exact root.
///
/// For `|z| > 1` the quotient is evaluated on the reversed polynomial in `1/z`,
/// which avoids overflow of `z^n`.
fn backward_error(monic: &[f64], z: Root) -> f64 {
    let modulus = z.0.hypot(z.1);
    if modulus <= 1.0 {
        let (value, _) = horner_with_derivative(monic, z);
        let bound = monic.iter().fold(0.0_f64, |acc, c| acc * modulus + c.abs());
        return value.0.hypot(value.1) / bound;
    }
    let Some(w) = complex_div((1.0, 0.0), z) else {
        return f64::INFINITY;
    };
    let w_modulus = w.0.hypot(w.1);
    let mut value = (0.0_f64, 0.0_f64);
    let mut bound = 0.0_f64;
    for &c in monic.iter().rev() {
        value = (
            value.0 * w.0 - value.1 * w.1 + c,
            value.0 * w.1 + value.1 * w.0,
        );
        bound = bound * w_modulus + c.abs();
    }
    value.0.hypot(value.1) / bound
}

/// Refines each root by Newton's method on the monic polynomial.
///
/// A step is taken only if it covers less than a tenth of the distance to the
/// nearest other root and lowers the backward error, so clustered or multiple
/// roots (where Newton converges slowly) are never made worse or merged.  Real
/// roots stay exactly real and conjugate pairs stay exactly conjugate, because
/// the arithmetic is sign-symmetric.
fn polish_roots(monic: &[f64], roots: &mut [Root]) {
    let snapshot: Vec<Root> = roots.to_vec();
    for (index, root) in roots.iter_mut().enumerate() {
        let nearest = snapshot
            .iter()
            .enumerate()
            .filter(|&(other, _)| other != index)
            .map(|(_, r)| (r.0 - root.0).hypot(r.1 - root.1))
            .fold(f64::INFINITY, f64::min);
        let mut current = *root;
        let mut current_error = backward_error(monic, current);
        for _ in 0..POLISH_ITERATIONS {
            let (value, derivative) = horner_with_derivative(monic, current);
            let Some(step) = complex_div(value, derivative) else {
                break;
            };
            let size = step.0.hypot(step.1);
            let small_step = size <= POLISH_STEP_FRACTION * nearest;
            if !small_step {
                break;
            }
            let candidate = (current.0 - step.0, current.1 - step.1);
            let candidate_error = backward_error(monic, candidate);
            let improved = candidate_error < current_error;
            if !improved {
                break;
            }
            current = candidate;
            current_error = candidate_error;
        }
        *root = current;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Polynomial `prod (x - r)` over complex roots given as conjugate-closed `(re, im)` pairs.
    fn from_roots(roots: &[Root]) -> Vec<f64> {
        let mut coefficients: Vec<Root> = vec![(1.0, 0.0)];
        for &(re, im) in roots {
            let mut next = vec![(0.0, 0.0); coefficients.len() + 1];
            for (k, &(cr, ci)) in coefficients.iter().enumerate() {
                next[k].0 += cr;
                next[k].1 += ci;
                next[k + 1].0 -= cr * re - ci * im;
                next[k + 1].1 -= cr * im + ci * re;
            }
            coefficients = next;
        }
        coefficients.iter().map(|c| c.0).collect()
    }

    fn assert_roots(got: &[Root], want: &[Root], tol: f64) {
        assert_eq!(got.len(), want.len(), "roots {got:?} vs {want:?}");
        for (g, w) in got.iter().zip(want) {
            assert!(
                (g.0 - w.0).abs() <= tol && (g.1 - w.1).abs() <= tol,
                "roots {got:?} vs {want:?}"
            );
        }
    }

    #[test]
    fn closed_forms_for_degrees_one_and_two() {
        assert_roots(
            &polynomial_roots(&[1.0, 2.0]).unwrap(),
            &[(-2.0, 0.0)],
            1e-15,
        );
        assert_roots(
            &polynomial_roots(&[4.0, -2.0]).unwrap(),
            &[(0.5, 0.0)],
            1e-15,
        );
        assert_roots(
            &polynomial_roots(&[1.0, -5.0, 6.0]).unwrap(),
            &[(2.0, 0.0), (3.0, 0.0)],
            1e-14,
        );
        // A degree-2 polynomial with complex roots used to return no roots at all.
        assert_roots(
            &polynomial_roots(&[1.0, 2.0, 5.0]).unwrap(),
            &[(-1.0, -2.0), (-1.0, 2.0)],
            1e-14,
        );
        // Double root.
        assert_roots(
            &polynomial_roots(&[1.0, 2.0, 1.0]).unwrap(),
            &[(-1.0, 0.0), (-1.0, 0.0)],
            1e-15,
        );
    }

    #[test]
    fn quadratic_formula_does_not_cancel() {
        // x² + 1e8 x + 1 has the roots ≈ -1e8 and ≈ -1e-8; the naive formula
        // loses all digits of the small one.
        let roots = polynomial_roots(&[1.0, 1e8, 1.0]).unwrap();
        let small = roots[1].0;
        assert!(
            ((small + 1e-8) / 1e-8).abs() < 1e-12,
            "small root {small:e}"
        );
        assert!(((roots[0].0 + 1e8) / 1e8).abs() < 1e-15);
    }

    #[test]
    fn cubic_with_one_real_and_two_complex_roots_is_not_fabricated() {
        // s³ + 3s² + 7s + 5 = (s + 1)(s² + 2s + 5): the Newton-and-deflate code
        // returned -1 plus two made-up real numbers.
        let roots = polynomial_roots(&[1.0, 3.0, 7.0, 5.0]).unwrap();
        assert_roots(&roots, &[(-1.0, -2.0), (-1.0, 0.0), (-1.0, 2.0)], 1e-9);
    }

    #[test]
    fn roots_of_unity_are_recovered_to_machine_precision() {
        // x⁸ - 1: the eighth roots of unity.
        let mut coefficients = vec![0.0; 9];
        coefficients[0] = 1.0;
        coefficients[8] = -1.0;
        let roots = polynomial_roots(&coefficients).unwrap();
        assert_eq!(roots.len(), 8);
        let mut want: Vec<Root> = (0..8)
            .map(|k| {
                let angle = std::f64::consts::PI * k as f64 / 4.0;
                (angle.cos(), angle.sin())
            })
            .collect();
        want.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        assert_roots(&roots, &want, 1e-12);
    }

    #[test]
    fn complex_roots_come_in_exactly_conjugate_pairs() {
        let want = [
            (-0.5, 3.0),
            (-0.5, -3.0),
            (-2.0, 1.5),
            (-2.0, -1.5),
            (-7.0, 0.0),
            (-0.1, 0.0),
            (1.0, 0.5),
            (1.0, -0.5),
        ];
        let roots = polynomial_roots(&from_roots(&want)).unwrap();
        assert_eq!(roots.len(), 8);
        for &(re, im) in roots.iter().filter(|r| r.1 > 0.0) {
            assert!(
                roots.contains(&(re, -im)),
                "no exact conjugate of ({re}, {im}) in {roots:?}"
            );
        }
        let mut sorted = want.to_vec();
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        assert_roots(&roots, &sorted, 1e-9);
    }

    #[test]
    fn widely_separated_time_scales_keep_their_relative_accuracy() {
        // Poles at -1e-3, -1, -1e3, -1e6 (nine decades apart).
        let want = [(-1e6, 0.0), (-1e3, 0.0), (-1.0, 0.0), (-1e-3, 0.0)];
        let roots = polynomial_roots(&from_roots(&want)).unwrap();
        assert_eq!(roots.len(), 4);
        for (got, exact) in roots.iter().zip(&want) {
            assert!(
                ((got.0 - exact.0) / exact.0).abs() < 1e-6 && got.1 == 0.0,
                "root {got:?} vs {exact:?}"
            );
        }
    }

    #[test]
    fn triple_root_is_returned_with_the_accuracy_the_conditioning_allows() {
        // (x + 1)³: a triple root is only determined to about eps^(1/3).
        let roots = polynomial_roots(&[1.0, 3.0, 3.0, 1.0]).unwrap();
        assert_eq!(roots.len(), 3);
        for root in &roots {
            assert!(
                (root.0 + 1.0).abs() < 1e-4 && root.1.abs() < 1e-4,
                "root {root:?}"
            );
        }
    }

    #[test]
    fn leading_and_trailing_zero_coefficients() {
        // 0 x³ + 0 x² + x + 3.
        assert_roots(
            &polynomial_roots(&[0.0, 0.0, 1.0, 3.0]).unwrap(),
            &[(-3.0, 0.0)],
            1e-15,
        );
        // s² + 3s = s (s + 3): the origin is a root.
        assert_roots(
            &polynomial_roots(&[1.0, 3.0, 0.0]).unwrap(),
            &[(-3.0, 0.0), (0.0, 0.0)],
            1e-15,
        );
        // 2 s³ = 0: a triple root at the origin.
        assert_roots(
            &polynomial_roots(&[2.0, 0.0, 0.0, 0.0]).unwrap(),
            &[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)],
            0.0,
        );
        // s⁴ + 5s³ + 8s² + 6s = s (s + 3)(s² + 2s + 2): roots 0, -3, -1 ± i.
        let roots = polynomial_roots(&[1.0, 5.0, 8.0, 6.0, 0.0]).unwrap();
        assert_roots(
            &roots,
            &[(-3.0, 0.0), (-1.0, -1.0), (-1.0, 1.0), (0.0, 0.0)],
            1e-9,
        );
    }

    #[test]
    fn constants_have_no_roots() {
        assert!(polynomial_roots(&[5.0]).unwrap().is_empty());
        assert!(polynomial_roots(&[0.0, 0.0, -7.0]).unwrap().is_empty());
    }

    #[test]
    fn invalid_polynomials_are_errors_not_panics() {
        assert!(matches!(
            polynomial_roots(&[]),
            Err(PhysicsError::InvalidInput { .. })
        ));
        assert!(matches!(
            polynomial_roots(&[0.0, 0.0]),
            Err(PhysicsError::InvalidInput { .. })
        ));
        assert!(matches!(
            polynomial_roots(&[1.0, f64::NAN, 2.0]),
            Err(PhysicsError::InvalidInput { .. })
        ));
        assert!(matches!(
            polynomial_roots(&[f64::INFINITY, 1.0]),
            Err(PhysicsError::InvalidInput { .. })
        ));
        // Monic coefficients overflow: 1e-300 x² + 1e300 x + 1e300.
        assert!(matches!(
            polynomial_roots(&[1e-300, 1e300, 1e300]),
            Err(PhysicsError::InvalidInput { .. })
        ));
    }

    #[test]
    fn every_returned_root_has_a_small_backward_error() {
        // A degree-12 polynomial with pseudo-random coefficients.
        let coefficients: Vec<f64> = (0..13)
            .map(|i| 1.0 + (1.7 * i as f64 + 0.3).sin())
            .collect();
        let roots = polynomial_roots(&coefficients).unwrap();
        assert_eq!(roots.len(), 12);
        let monic: Vec<f64> = coefficients.iter().map(|c| c / coefficients[0]).collect();
        for &root in &roots {
            let error = backward_error(&monic, root);
            assert!(error < 1e-10, "root {root:?} has backward error {error:e}");
        }
    }

    #[test]
    fn backward_error_distinguishes_roots_from_non_roots() {
        let monic = [1.0, 3.0, 7.0, 5.0];
        assert!(backward_error(&monic, (-1.0, 0.0)) < 1e-15);
        assert!(backward_error(&monic, (-1.0, 2.0)) < 1e-15);
        // The "roots" the old code made up: far from zero.
        assert!(backward_error(&monic, (0.1, 0.0)) > 0.1);
        // Large roots are judged relatively: 1e6 is a root of x - 1e6 whatever its size.
        assert!(backward_error(&[1.0, -1e6], (1e6, 0.0)) < 1e-15);
        assert!(backward_error(&[1.0, -1e6], (1.1e6, 0.0)) > 1e-2);
    }
}
