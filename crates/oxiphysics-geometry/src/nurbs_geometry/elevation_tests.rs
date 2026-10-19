// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for Bézier decomposition and degree elevation.
//!
//! The earlier implementation sampled the curve instead of extracting control
//! points (dropping the weights), blended the control points of a Bézier
//! segment with swapped coefficients and rebuilt the elevated curve on a
//! uniform knot vector, so even a single quadratic Bézier segment changed
//! its shape.

use super::*;
use crate::Error;

const HALF_SQRT_2: f64 = std::f64::consts::FRAC_1_SQRT_2;

fn cp(x: f64, y: f64, z: f64, w: f64) -> NurbsControlPoint {
    NurbsControlPoint::new(x, y, z, w)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_norm(vec3_sub(a, b))
}

/// Bézier curve through the given control points on the parameter range `[0, 1]`.
fn bezier_curve(points: &[NurbsControlPoint]) -> NurbsCurve {
    let degree = points.len() - 1;
    NurbsCurve::new(
        NurbsKnotVector::clamped_uniform(points.len(), degree),
        points.to_vec(),
        degree,
    )
}

/// Quadratic rational Bézier segment of the unit quarter circle.
fn quarter_circle() -> NurbsCurve {
    bezier_curve(&[
        cp(1.0, 0.0, 0.0, 1.0),
        cp(1.0, 1.0, 0.0, HALF_SQRT_2),
        cp(0.0, 1.0, 0.0, 1.0),
    ])
}

/// Largest distance between the two curves over `samples + 1` equally spaced
/// parameters of `[t_start, t_end]`.
fn max_deviation(a: &NurbsCurve, b: &NurbsCurve, t_start: f64, t_end: f64, samples: usize) -> f64 {
    (0..=samples)
        .map(|i| {
            let t = t_start + (t_end - t_start) * i as f64 / samples as f64;
            distance(a.evaluate(t), b.evaluate(t))
        })
        .fold(0.0, f64::max)
}

/// Parameters of a curve on `[0, 1]` that are never a knot of the test curves.
fn generic_parameters() -> impl Iterator<Item = f64> {
    (0..60).map(|i| (i as f64 + 0.5) / 60.0)
}

fn max_derivative_deviation(a: &NurbsCurve, b: &NurbsCurve) -> f64 {
    generic_parameters()
        .map(|t| {
            let da = a.derivative(t);
            let db = b.derivative(t);
            distance(da, db) / (1.0 + vec3_norm(da))
        })
        .fold(0.0, f64::max)
}

/// Cubic curve with two simple interior knots and non-uniform weights.
fn rational_cubic() -> NurbsCurve {
    NurbsCurve::new(
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
    )
}

fn knot_multiplicities(curve: &NurbsCurve) -> Vec<(f64, usize)> {
    let mut distinct: Vec<f64> = Vec::new();
    for &k in &curve.knots.knots {
        if distinct.last().is_none_or(|&last| k != last) {
            distinct.push(k);
        }
    }
    distinct
        .into_iter()
        .map(|k| (k, curve.knots.multiplicity(k)))
        .collect()
}

// ---------------------------------------------------------------------------
// Bézier segment elevation
// ---------------------------------------------------------------------------

#[test]
fn bezier_elevation_blends_with_the_correct_coefficients() {
    // Q_i = i/(p+1) P_{i-1} + (1 - i/(p+1)) P_i for the parabola (0,0), (1,1), (2,0).
    let pts = [
        cp(0.0, 0.0, 0.0, 1.0),
        cp(1.0, 1.0, 0.0, 1.0),
        cp(2.0, 0.0, 0.0, 1.0),
    ];
    let elevated = DegreeElevation::elevate_bezier(&pts);
    let expected = [
        [0.0, 0.0],
        [2.0 / 3.0, 2.0 / 3.0],
        [4.0 / 3.0, 2.0 / 3.0],
        [2.0, 0.0],
    ];
    assert_eq!(elevated.len(), 4);
    for (q, e) in elevated.iter().zip(expected) {
        assert!(
            (q.point[0] - e[0]).abs() < 1e-15 && (q.point[1] - e[1]).abs() < 1e-15,
            "control point {:?}, expected {e:?}",
            q.point
        );
        assert!((q.weight - 1.0).abs() < 1e-15);
    }
    // The shape does not change: y(0.5) of the parabola is 0.5.
    let original = bezier_curve(&pts).evaluate(0.5);
    let raised = bezier_curve(&elevated).evaluate(0.5);
    assert!(
        (original[1] - 0.5).abs() < 1e-14,
        "y(0.5) = {}",
        original[1]
    );
    assert!((raised[1] - 0.5).abs() < 1e-14, "y(0.5) = {}", raised[1]);
}

#[test]
fn rational_bezier_elevation_keeps_the_quarter_circle() {
    let circle = quarter_circle();
    let elevated = DegreeElevation::elevate_bezier(&circle.control_points);
    assert_eq!(elevated.len(), 4);
    // Weights blend linearly: w_1 = 1/3 w_0 + 2/3 w_1.
    assert!((elevated[1].weight - (1.0 / 3.0 + 2.0 / 3.0 * HALF_SQRT_2)).abs() < 1e-15);
    assert!((elevated[2].weight - (2.0 / 3.0 * HALF_SQRT_2 + 1.0 / 3.0)).abs() < 1e-15);
    let raised = bezier_curve(&elevated);
    for k in 1..=9 {
        let t = f64::from(k) / 10.0;
        let on_original = circle.evaluate(t);
        let on_raised = raised.evaluate(t);
        assert!(
            distance(on_original, on_raised) < 1e-12,
            "t = {t}: {on_original:?} vs {on_raised:?}"
        );
        let radius = (on_raised[0].powi(2) + on_raised[1].powi(2)).sqrt();
        assert!((radius - 1.0).abs() < 1e-12, "t = {t}: radius {radius}");
    }
}

#[test]
fn bezier_elevation_of_a_cubic_matches_the_original_everywhere() {
    let pts = [
        cp(0.0, 0.0, 0.0, 1.0),
        cp(1.0, 3.0, -1.0, 2.0),
        cp(2.0, -2.0, 1.0, 0.5),
        cp(4.0, 1.0, 2.0, 1.5),
    ];
    let original = bezier_curve(&pts);
    let raised = bezier_curve(&DegreeElevation::elevate_bezier(&pts));
    assert!(max_deviation(&original, &raised, 0.0, 1.0, 200) < 1e-12);
}

#[test]
fn bezier_elevation_of_trivial_inputs() {
    assert!(DegreeElevation::elevate_bezier(&[]).is_empty());

    let single = [cp(1.0, 2.0, 3.0, 2.0)];
    let doubled = DegreeElevation::elevate_bezier(&single);
    assert_eq!(doubled.len(), 2);
    for q in &doubled {
        assert_eq!(q.point, [1.0, 2.0, 3.0]);
        assert_eq!(q.weight, 2.0);
    }

    let line = [cp(0.0, 0.0, 0.0, 1.0), cp(2.0, 4.0, 6.0, 1.0)];
    let raised = DegreeElevation::elevate_bezier(&line);
    assert_eq!(raised.len(), 3);
    assert!(distance(raised[1].point, [1.0, 2.0, 3.0]) < 1e-15);
}

// ---------------------------------------------------------------------------
// Bézier decomposition
// ---------------------------------------------------------------------------

#[test]
fn decomposing_the_circle_gives_four_weighted_quarter_arcs() {
    let circle = NurbsCurve::circle(1.0);
    let segments = BezierDecomposition::try_decompose(&circle).expect("well-formed circle");
    assert_eq!(segments.len(), 4);
    for (k, segment) in segments.iter().enumerate() {
        assert!((segment.t_start - 0.25 * k as f64).abs() < 1e-15);
        assert!((segment.t_end - 0.25 * (k + 1) as f64).abs() < 1e-15);
        let weights: Vec<f64> = segment.control_points.iter().map(|c| c.weight).collect();
        assert_eq!(weights.len(), 3);
        assert!(
            (weights[0] - 1.0).abs() < 1e-14
                && (weights[1] - HALF_SQRT_2).abs() < 1e-14
                && (weights[2] - 1.0).abs() < 1e-14,
            "segment {k} has weights {weights:?}"
        );
        // The segment is the arc of the circle between its parameters.
        let arc = bezier_curve(&segment.control_points);
        for i in 0..=20 {
            let s = f64::from(i) / 20.0;
            let t = segment.t_start + s * (segment.t_end - segment.t_start);
            assert!(distance(arc.evaluate(s), circle.evaluate(t)) < 1e-12);
        }
    }
    let flat = BezierDecomposition::decompose(&circle);
    assert_eq!(flat.len(), 4);
    assert!(flat.iter().all(|segment| segment.len() == 3));
}

#[test]
fn decomposition_inserts_the_missing_knots() {
    let curve = rational_cubic();
    let segments = BezierDecomposition::try_decompose(&curve).expect("well-formed curve");
    assert_eq!(segments.len(), 3);
    let ranges = [(0.0, 0.3), (0.3, 0.6), (0.6, 1.0)];
    for (segment, (t_start, t_end)) in segments.iter().zip(ranges) {
        assert!((segment.t_start - t_start).abs() < 1e-15);
        assert!((segment.t_end - t_end).abs() < 1e-15);
        assert_eq!(segment.control_points.len(), 4);
        assert!(segment.control_points.iter().all(|c| c.weight > 0.0));
        let piece = bezier_curve(&segment.control_points);
        for i in 0..=40 {
            let s = f64::from(i) / 40.0;
            let t = t_start + s * (t_end - t_start);
            let (expected, found) = (curve.evaluate(t), piece.evaluate(s));
            assert!(
                distance(expected, found) < 1e-12,
                "t = {t}: {expected:?} vs {found:?}"
            );
        }
    }
    // Neighbouring segments of a continuous curve share their common end point.
    for pair in segments.windows(2) {
        let (left, right) = (&pair[0], &pair[1]);
        let end = left.control_points[3];
        let start = right.control_points[0];
        assert!(distance(end.point, start.point) < 1e-12);
        assert!((end.weight - start.weight).abs() < 1e-12);
    }
}

#[test]
fn decomposing_a_bezier_curve_returns_its_control_points() {
    let curve = quarter_circle();
    let segments = BezierDecomposition::decompose(&curve);
    assert_eq!(segments.len(), 1);
    for (found, original) in segments[0].iter().zip(&curve.control_points) {
        assert!(distance(found.point, original.point) < 1e-15);
        assert!((found.weight - original.weight).abs() < 1e-15);
    }
}

#[test]
fn decomposition_handles_knots_of_mixed_multiplicity() {
    // Degree 3 with a double and a simple interior knot.
    let cubic = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.5, 1.5),
            cp(2.5, 1.5, -1.0, 0.7),
            cp(3.0, -1.0, 0.3, 2.0),
            cp(4.5, 0.5, 1.2, 1.0),
            cp(5.0, 2.0, 0.0, 3.0),
            cp(6.0, 1.0, -1.0, 0.8),
        ],
        3,
    );
    // Degree 4 with a simple and a triple interior knot.
    let quartic = NurbsCurve::new(
        NurbsKnotVector::new(vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0,
        ]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 1.0, 0.0, 2.0),
            cp(2.0, -1.0, 1.0, 0.5),
            cp(3.0, 2.0, 0.5, 1.5),
            cp(4.0, 0.0, -1.0, 1.0),
            cp(5.0, 1.0, 0.0, 2.5),
            cp(6.0, -2.0, 1.0, 1.0),
            cp(7.0, 0.5, 0.2, 0.6),
            cp(8.0, 1.0, 0.0, 1.0),
        ],
        4,
    );
    for (curve, ranges) in [
        (&cubic, vec![(0.0, 0.3), (0.3, 0.7), (0.7, 1.0)]),
        (&quartic, vec![(0.0, 0.2), (0.2, 0.5), (0.5, 1.0)]),
    ] {
        let segments = BezierDecomposition::try_decompose(curve).expect("well-formed curve");
        assert_eq!(segments.len(), ranges.len());
        for (segment, (t_start, t_end)) in segments.iter().zip(ranges) {
            assert!((segment.t_start - t_start).abs() < 1e-15);
            assert!((segment.t_end - t_end).abs() < 1e-15);
            assert_eq!(segment.control_points.len(), curve.degree + 1);
            let piece = bezier_curve(&segment.control_points);
            for i in 0..=30 {
                let s = f64::from(i) / 30.0;
                let t = t_start + s * (t_end - t_start);
                assert!(
                    distance(piece.evaluate(s), curve.evaluate(t)) < 1e-11,
                    "degree {}, t = {t}",
                    curve.degree
                );
            }
        }
    }
}

#[test]
fn decomposition_keeps_both_sides_of_a_jump() {
    // The knot 0.5 has multiplicity degree + 1: the curve jumps there.
    let pts = vec![
        cp(0.0, 0.0, 0.0, 1.0),
        cp(1.0, 1.0, 0.0, 1.0),
        cp(2.0, 0.0, 0.0, 1.0),
        cp(5.0, 0.0, 0.0, 1.0),
        cp(6.0, 1.0, 0.0, 2.0),
        cp(7.0, 0.0, 0.0, 1.0),
    ];
    let curve = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]),
        pts.clone(),
        2,
    );
    let segments = BezierDecomposition::decompose(&curve);
    assert_eq!(segments.len(), 2);
    for (segment, original) in segments.iter().zip(pts.chunks(3)) {
        for (found, expected) in segment.iter().zip(original) {
            assert!(distance(found.point, expected.point) < 1e-15);
            assert!((found.weight - expected.weight).abs() < 1e-15);
        }
    }
}

#[test]
fn decomposition_of_malformed_curves_reports_an_error() {
    let mut curve = quarter_circle();
    curve.knots = NurbsKnotVector {
        knots: vec![0.0, 0.0, 1.0, 1.0],
    };
    assert!(matches!(
        BezierDecomposition::try_decompose(&curve),
        Err(Error::LengthMismatch {
            expected: 6,
            actual: 4
        })
    ));
    assert!(BezierDecomposition::decompose(&curve).is_empty());

    let mut weights = quarter_circle();
    weights.control_points[1].weight = 0.0;
    assert!(BezierDecomposition::try_decompose(&weights).is_err());

    let mut too_few = quarter_circle();
    too_few.control_points.truncate(2);
    assert!(matches!(
        BezierDecomposition::try_decompose(&too_few),
        Err(Error::TooFewPoints { .. })
    ));
}

// ---------------------------------------------------------------------------
// Curve elevation
// ---------------------------------------------------------------------------

#[test]
fn elevating_the_quarter_circle_keeps_its_shape() {
    let circle = quarter_circle();
    let elevated = DegreeElevation::try_elevate_curve(&circle).expect("well-formed curve");
    assert_eq!(elevated.degree, 3);
    assert_eq!(
        elevated.knots.knots,
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    );
    assert_eq!(elevated.control_points.len(), 4);
    assert!(max_deviation(&circle, &elevated, 0.0, 1.0, 200) < 1e-12);
    assert!(max_derivative_deviation(&circle, &elevated) < 1e-9);
}

#[test]
fn elevating_the_full_circle_raises_every_knot_multiplicity_by_one() {
    let circle = NurbsCurve::circle(2.0);
    let elevated = DegreeElevation::elevate_curve(&circle);
    assert_eq!(elevated.degree, 3);
    assert_eq!(
        knot_multiplicities(&elevated),
        vec![(0.0, 4), (0.25, 3), (0.5, 3), (0.75, 3), (1.0, 4)]
    );
    assert_eq!(elevated.control_points.len(), 13);
    assert_eq!(
        elevated.knots.len(),
        elevated.control_points.len() + elevated.degree + 1
    );
    assert!(max_deviation(&circle, &elevated, 0.0, 1.0, 200) < 1e-10);
    assert!(max_derivative_deviation(&circle, &elevated) < 1e-8);
    for i in 0..=200 {
        let p = elevated.evaluate(f64::from(i) / 200.0);
        let radius = (p[0] * p[0] + p[1] * p[1]).sqrt();
        assert!((radius - 2.0).abs() < 1e-10, "radius {radius}");
    }
}

#[test]
fn elevating_a_multi_span_rational_curve_keeps_geometry_and_continuity() {
    let curve = rational_cubic();
    let elevated = DegreeElevation::try_elevate_curve(&curve).expect("well-formed curve");
    assert_eq!(elevated.degree, 4);
    assert_eq!(
        knot_multiplicities(&elevated),
        vec![(0.0, 5), (0.3, 2), (0.6, 2), (1.0, 5)]
    );
    assert_eq!(elevated.control_points.len(), 9);
    assert!(max_deviation(&curve, &elevated, 0.0, 1.0, 200) < 1e-10);
    // The parametrisation and the smoothness at the knots are unchanged.
    assert!(max_derivative_deviation(&curve, &elevated) < 1e-8);
}

#[test]
fn elevation_keeps_the_continuity_of_a_double_knot() {
    // A double knot in a cubic is a C^1 join; in the quartic it must stay C^1,
    // so the knot has multiplicity three.
    let curve = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.0, 1.0),
            cp(2.0, -1.0, 1.0, 1.5),
            cp(3.0, 1.0, 0.5, 1.0),
            cp(4.0, 0.0, -0.5, 2.0),
            cp(5.0, 2.0, 0.0, 1.0),
        ],
        3,
    );
    let elevated = DegreeElevation::try_elevate_curve(&curve).expect("well-formed curve");
    assert_eq!(elevated.knots.multiplicity(0.5), 3);
    assert_eq!(elevated.control_points.len(), 8);
    assert!(max_deviation(&curve, &elevated, 0.0, 1.0, 200) < 1e-10);
    assert!(max_derivative_deviation(&curve, &elevated) < 1e-8);
}

#[test]
fn elevating_twice_equals_elevating_by_two() {
    let curve = rational_cubic();
    let once = DegreeElevation::try_elevate_curve(&curve).expect("well-formed curve");
    let twice = DegreeElevation::try_elevate_curve(&once).expect("well-formed curve");
    let direct = DegreeElevation::try_elevate_curve_times(&curve, 2).expect("well-formed curve");
    assert_eq!(direct.degree, 5);
    assert_eq!(twice.degree, 5);
    assert_eq!(twice.knots.knots, direct.knots.knots);
    assert_eq!(twice.control_points.len(), direct.control_points.len());
    for (a, b) in twice.control_points.iter().zip(&direct.control_points) {
        assert!(distance(a.point, b.point) < 1e-10);
        assert!((a.weight - b.weight).abs() < 1e-10);
    }
    assert!(max_deviation(&curve, &direct, 0.0, 1.0, 200) < 1e-10);
}

#[test]
fn elevating_by_zero_degrees_returns_the_same_curve() {
    let curve = rational_cubic();
    let same = DegreeElevation::try_elevate_curve_times(&curve, 0).expect("well-formed curve");
    assert_eq!(same.degree, curve.degree);
    assert_eq!(same.knots.knots, curve.knots.knots);
    assert!(max_deviation(&curve, &same, 0.0, 1.0, 100) < 1e-13);
}

#[test]
fn elevating_an_unclamped_curve_clamps_it_first() {
    // A uniform knot vector is not clamped: the curve lives on [0.25, 0.75].
    let curve = NurbsCurve::new(
        NurbsKnotVector::uniform(6, 2),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 2.0, 0.0, 1.5),
            cp(2.0, -1.0, 0.5, 1.0),
            cp(3.0, 1.0, 1.0, 2.0),
            cp(4.0, 0.0, -1.0, 1.0),
            cp(5.0, 2.0, 0.0, 0.5),
        ],
        2,
    );
    let (t_start, t_end) = (
        curve.knots.domain_start(curve.degree),
        curve.knots.domain_end(curve.degree, curve.n_ctrl()),
    );
    assert!((t_start - 0.25).abs() < 1e-15 && (t_end - 0.75).abs() < 1e-15);

    let elevated = DegreeElevation::try_elevate_curve(&curve).expect("well-formed curve");
    assert_eq!(elevated.degree, 3);
    let knots = &elevated.knots.knots;
    assert!(knots[..4].iter().all(|&k| (k - t_start).abs() < 1e-14));
    assert!(
        knots[knots.len() - 4..]
            .iter()
            .all(|&k| (k - t_end).abs() < 1e-14)
    );
    assert!(max_deviation(&curve, &elevated, t_start, t_end, 200) < 1e-10);

    // The domain [0.25, 0.75] holds the knots 0.375, 0.5 and 0.625: four spans.
    let segments = BezierDecomposition::try_decompose(&curve).expect("well-formed curve");
    assert_eq!(segments.len(), 4);
    assert!((segments[0].t_start - t_start).abs() < 1e-14);
    assert!((segments[3].t_end - t_end).abs() < 1e-14);
    for segment in &segments {
        let piece = bezier_curve(&segment.control_points);
        for i in 0..=20 {
            let s = f64::from(i) / 20.0;
            let t = segment.t_start + s * (segment.t_end - segment.t_start);
            assert!(distance(piece.evaluate(s), curve.evaluate(t)) < 1e-12);
        }
    }
}

#[test]
fn elevating_a_curve_with_a_jump_keeps_both_pieces() {
    let curve = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
        vec![
            cp(0.0, 0.0, 0.0, 1.0),
            cp(1.0, 1.0, 0.0, 2.0),
            cp(3.0, 0.0, 0.0, 1.0),
            cp(4.0, 2.0, 0.0, 1.0),
        ],
        1,
    );
    let elevated = DegreeElevation::try_elevate_curve(&curve).expect("well-formed curve");
    assert_eq!(elevated.degree, 2);
    assert_eq!(elevated.knots.multiplicity(0.5), 3);
    for t in [0.05, 0.2, 0.35, 0.45, 0.55, 0.7, 0.85, 0.95] {
        assert!(
            distance(curve.evaluate(t), elevated.evaluate(t)) < 1e-12,
            "t = {t}"
        );
    }
}

#[test]
fn elevation_reports_malformed_curves() {
    let line = NurbsCurve::new(
        NurbsKnotVector::new(vec![0.0, 1.0]),
        vec![cp(0.0, 0.0, 0.0, 1.0)],
        0,
    );
    assert!(matches!(
        DegreeElevation::try_elevate_curve(&line),
        Err(Error::Unsupported { .. })
    ));

    let mut wrong_length = quarter_circle();
    wrong_length.knots.knots.pop();
    assert!(matches!(
        DegreeElevation::try_elevate_curve(&wrong_length),
        Err(Error::LengthMismatch { .. })
    ));

    let mut decreasing = quarter_circle();
    decreasing.knots.knots = vec![0.0, 0.0, 0.5, 0.4, 1.0, 1.0];
    assert!(DegreeElevation::try_elevate_curve(&decreasing).is_err());

    let mut bad_weight = quarter_circle();
    bad_weight.control_points[1].weight = -1.0;
    assert!(DegreeElevation::try_elevate_curve(&bad_weight).is_err());

    let mut not_finite = quarter_circle();
    not_finite.control_points[0].point[0] = f64::NAN;
    assert!(DegreeElevation::try_elevate_curve(&not_finite).is_err());

    let mut no_domain = quarter_circle();
    no_domain.knots.knots = vec![0.5; 6];
    assert!(matches!(
        DegreeElevation::try_elevate_curve(&no_domain),
        Err(Error::DegenerateGeometry { .. })
    ));

    let empty = NurbsCurve {
        knots: NurbsKnotVector { knots: Vec::new() },
        control_points: Vec::new(),
        degree: 2,
    };
    assert!(DegreeElevation::try_elevate_curve(&empty).is_err());
}

#[test]
fn infallible_elevation_never_returns_a_plausible_wrong_curve() {
    let mut bad = quarter_circle();
    bad.control_points[2].weight = f64::NAN;
    let poisoned = DegreeElevation::elevate_curve(&bad);
    assert_eq!(poisoned.degree, bad.degree);
    assert_eq!(poisoned.knots.knots, bad.knots.knots);
    assert_eq!(poisoned.control_points.len(), bad.control_points.len());
    assert!(
        poisoned
            .control_points
            .iter()
            .all(|c| c.weight.is_nan() && c.point.iter().all(|x| x.is_nan()))
    );
    assert!(poisoned.evaluate(0.5)[0].is_nan());

    let good = DegreeElevation::elevate_curve(&quarter_circle());
    assert_eq!(good.degree, 3);
    assert!(good.control_points.iter().all(|c| c.weight.is_finite()));
}
