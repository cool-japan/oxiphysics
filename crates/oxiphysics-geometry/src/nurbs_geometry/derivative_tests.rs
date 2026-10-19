// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the derivatives of the B-spline basis and of NURBS curves.
//!
//! `NurbsBasis::evaluate_all_derivs` omitted the values of the lower-degree basis
//! functions and the factor `p!/(p-k)!` of Piegl & Tiller's Algorithm A2.3. It was
//! right for degree 1 and for the middle of a uniform Bézier span only, so the
//! derivative, tangent, curvature and arc length of every other curve were wrong.

use super::*;

const STEP: f64 = 1.0e-6;

/// Knot vectors with simple, double and clamped knots, with their degree and
/// number of control points.
fn knot_vectors() -> Vec<(NurbsKnotVector, usize, usize)> {
    vec![
        (NurbsKnotVector::clamped_uniform(2, 1), 1, 2),
        (NurbsKnotVector::clamped_uniform(3, 2), 2, 3),
        (NurbsKnotVector::clamped_uniform(7, 3), 3, 7),
        (
            NurbsKnotVector::new(vec![
                0.0, 0.0, 0.0, 0.0, 0.2, 0.55, 0.55, 0.9, 1.0, 1.0, 1.0, 1.0,
            ]),
            3,
            8,
        ),
        (
            NurbsKnotVector::new(vec![
                0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.45, 1.0, 1.0, 1.0, 1.0, 1.0,
            ]),
            4,
            7,
        ),
        (NurbsKnotVector::uniform(8, 2), 2, 8),
    ]
}

/// Parameters of the domain that are away from every knot.
fn parameters(knots: &NurbsKnotVector, degree: usize, n_ctrl: usize) -> Vec<f64> {
    let (start, end) = (knots.domain_start(degree), knots.domain_end(degree, n_ctrl));
    (0..40)
        .map(|i| start + (end - start) * (i as f64 + 0.37) / 40.0)
        .filter(|&t| knots.knots.iter().all(|&k| (k - t).abs() > 50.0 * STEP))
        .collect()
}

#[test]
fn derivatives_of_the_basis_functions_sum_to_zero() {
    for (knots, degree, n_ctrl) in knot_vectors() {
        let d = degree.min(3);
        for t in parameters(&knots, degree, n_ctrl) {
            let ders = NurbsBasis::evaluate_all_derivs(&knots, t, degree, n_ctrl, d);
            assert_eq!(ders.len(), d + 1);
            let value: f64 = ders[0].iter().sum();
            assert!(
                (value - 1.0).abs() < 1e-12,
                "p = {degree}, t = {t}: {value}"
            );
            for (k, row) in ders.iter().enumerate().skip(1) {
                let sum: f64 = row.iter().sum();
                assert!(
                    sum.abs() < 1e-8 * (1.0 + row.iter().map(|v| v.abs()).sum::<f64>()),
                    "p = {degree}, t = {t}, order {k}: sum {sum}"
                );
            }
        }
    }
}

#[test]
fn first_derivatives_of_the_basis_match_finite_differences() {
    for (knots, degree, n_ctrl) in knot_vectors() {
        for t in parameters(&knots, degree, n_ctrl) {
            let span = knots.find_span(t, degree, n_ctrl);
            if knots.find_span(t - STEP, degree, n_ctrl) != span
                || knots.find_span(t + STEP, degree, n_ctrl) != span
            {
                continue;
            }
            let (_, dn) = NurbsBasis::evaluate_deriv(&knots, t, degree, n_ctrl);
            let above = NurbsBasis::evaluate(&knots, t + STEP, degree, n_ctrl);
            let below = NurbsBasis::evaluate(&knots, t - STEP, degree, n_ctrl);
            for j in 0..=degree {
                let expected = (above[j] - below[j]) / (2.0 * STEP);
                assert!(
                    (dn[j] - expected).abs() < 1e-5 * (1.0 + expected.abs()),
                    "p = {degree}, t = {t}, function {j}: {} vs {expected}",
                    dn[j]
                );
            }
        }
    }
}

#[test]
fn second_derivatives_of_the_basis_match_finite_differences() {
    for (knots, degree, n_ctrl) in knot_vectors() {
        if degree < 2 {
            continue;
        }
        for t in parameters(&knots, degree, n_ctrl) {
            let span = knots.find_span(t, degree, n_ctrl);
            if knots.find_span(t - STEP, degree, n_ctrl) != span
                || knots.find_span(t + STEP, degree, n_ctrl) != span
            {
                continue;
            }
            let ders = NurbsBasis::evaluate_all_derivs(&knots, t, degree, n_ctrl, 2);
            let above = NurbsBasis::evaluate_all_derivs(&knots, t + STEP, degree, n_ctrl, 1);
            let below = NurbsBasis::evaluate_all_derivs(&knots, t - STEP, degree, n_ctrl, 1);
            for j in 0..=degree {
                let expected = (above[1][j] - below[1][j]) / (2.0 * STEP);
                assert!(
                    (ders[2][j] - expected).abs() < 1e-4 * (1.0 + expected.abs()),
                    "p = {degree}, t = {t}, function {j}: {} vs {expected}",
                    ders[2][j]
                );
            }
        }
    }
}

fn cp(x: f64, y: f64, z: f64, w: f64) -> NurbsControlPoint {
    NurbsControlPoint::new(x, y, z, w)
}

fn central_difference(curve: &NurbsCurve, t: f64) -> [f64; 3] {
    let (above, below) = (curve.evaluate(t + STEP), curve.evaluate(t - STEP));
    [
        (above[0] - below[0]) / (2.0 * STEP),
        (above[1] - below[1]) / (2.0 * STEP),
        (above[2] - below[2]) / (2.0 * STEP),
    ]
}

fn assert_derivative_matches(curve: &NurbsCurve, name: &str) {
    let (start, end) = (
        curve.knots.domain_start(curve.degree),
        curve.knots.domain_end(curve.degree, curve.n_ctrl()),
    );
    for i in 0..60 {
        let t = start + (end - start) * (i as f64 + 0.41) / 60.0;
        if curve
            .knots
            .knots
            .iter()
            .any(|&k| (k - t).abs() < 50.0 * STEP)
        {
            continue;
        }
        let analytic = curve.derivative(t);
        let numeric = central_difference(curve, t);
        let error = vec3_norm(vec3_sub(analytic, numeric));
        assert!(
            error < 1e-5 * (1.0 + vec3_norm(numeric)),
            "{name}, t = {t}: {analytic:?} vs {numeric:?}"
        );
    }
}

#[test]
fn the_derivative_of_a_polynomial_curve_matches_finite_differences() {
    let quadratic = NurbsCurve::new(
        NurbsKnotVector::clamped_uniform(3, 2),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.5, 1.0),
            cp(3.0, 0.0, -1.0, 1.0),
        ],
        2,
    );
    assert_derivative_matches(&quadratic, "quadratic Bézier");
    let cubic = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.5, 1.0),
            cp(2.5, 1.5, -1.0, 1.0),
            cp(3.0, -1.0, 0.3, 1.0),
            cp(4.5, 0.5, 1.2, 1.0),
            cp(5.0, 2.0, 0.0, 1.0),
        ],
        3,
    );
    assert_derivative_matches(&cubic, "multi-span cubic");
}

#[test]
fn the_derivative_of_a_rational_curve_matches_finite_differences() {
    assert_derivative_matches(&NurbsCurve::circle(2.0), "circle");
    let rational = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.5, 2.0),
            cp(2.5, 1.5, -1.0, 0.5),
            cp(3.0, -1.0, 0.3, 3.0),
            cp(4.5, 0.5, 1.2, 1.5),
            cp(5.0, 2.0, 0.0, 1.0),
        ],
        3,
    );
    assert_derivative_matches(&rational, "rational cubic");
}

#[test]
fn the_tangent_of_the_circle_is_perpendicular_to_the_radius() {
    let circle = NurbsCurve::circle(1.5);
    for i in 1..40 {
        let t = f64::from(i) / 40.0 + 0.0123;
        let p = circle.evaluate(t);
        let tangent = circle.tangent(t);
        assert!((vec3_norm(tangent) - 1.0).abs() < 1e-9);
        assert!(
            vec3_dot(tangent, p).abs() < 1e-8,
            "t = {t}: tangent {tangent:?} at {p:?}"
        );
    }
}

#[test]
fn the_arc_length_of_the_circle_is_two_pi_r() {
    let circle = NurbsCurve::circle(2.0);
    // The speed is smooth within each knot span only, so the Gauss rule is applied
    // to sixteen pieces that do not straddle a knot.
    let mut length = 0.0;
    for piece in 0..16 {
        let a = f64::from(piece) / 16.0;
        let b = f64::from(piece + 1) / 16.0;
        length += gauss_legendre_integral(a, b, 4, |t| vec3_norm(circle.derivative(t)));
    }
    let expected = 2.0 * std::f64::consts::PI * 2.0;
    assert!(
        (length - expected).abs() < 1e-5,
        "length {length}, expected {expected}"
    );
}
