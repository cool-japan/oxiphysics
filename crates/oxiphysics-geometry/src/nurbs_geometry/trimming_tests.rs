// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for trimmed NURBS surfaces.
//!
//! `NurbsTrimming::is_inside` used to ignore its arguments and return `true`, so
//! a trimmed surface was never trimmed.

use super::*;

fn trim_cp(u: f64, v: f64) -> NurbsControlPoint {
    NurbsControlPoint::unweighted(u, v, 0.0)
}

/// A degree-1 trim curve through the given parameter points.
fn polyline_curve(points: &[(f64, f64)]) -> NurbsTrimCurve {
    let n = points.len();
    let mut knots = vec![0.0];
    knots.extend((0..n).map(|i| i as f64));
    knots.push((n - 1) as f64);
    NurbsTrimCurve::new(
        points.iter().map(|&(u, v)| trim_cp(u, v)).collect(),
        NurbsKnotVector::new(knots),
        1,
    )
}

/// Closed counter-clockwise loop around `[lo, hi]^2`.
fn square_loop(lo: f64, hi: f64) -> NurbsTrimCurve {
    polyline_curve(&[(lo, lo), (hi, lo), (hi, hi), (lo, hi), (lo, lo)])
}

/// Counter-clockwise circle as a trim curve.
fn circle_loop(center: (f64, f64), radius: f64) -> NurbsTrimCurve {
    let circle = NurbsCurve::circle(radius);
    let control_points = circle
        .control_points
        .iter()
        .map(|c| {
            NurbsControlPoint::new(c.point[0] + center.0, c.point[1] + center.1, 0.0, c.weight)
        })
        .collect();
    NurbsTrimCurve::new(control_points, circle.knots.clone(), circle.degree)
}

fn unit_surface() -> NurbsSurface {
    NurbsSurface::flat_plane(0.0, 1.0, 0.0, 1.0)
}

/// The unit square trimmed to `[0.2, 0.8]^2` with a circular hole of radius 0.15
/// at `(0.5, 0.5)`; the hole runs clockwise.
fn square_with_hole() -> NurbsTrimming {
    let mut trimmed = NurbsTrimming::new(unit_surface());
    trimmed.add_trim_curve(square_loop(0.2, 0.8));
    trimmed.add_trim_curve(circle_loop((0.5, 0.5), 0.15).reversed());
    trimmed
}

#[test]
fn the_outer_loop_and_a_hole_trim_the_surface() {
    let trimmed = square_with_hole();
    // Inside the hole.
    assert!(trimmed.evaluate(0.5, 0.5).is_none());
    assert!(trimmed.evaluate(0.5, 0.62).is_none());
    // Outside the outer loop.
    assert!(trimmed.evaluate(0.1, 0.1).is_none());
    assert!(trimmed.evaluate(0.9, 0.5).is_none());
    assert!(trimmed.evaluate(0.5, 0.1).is_none());
    // Kept: between the outer loop and the hole.
    for (u, v) in [
        (0.3, 0.3),
        (0.7, 0.7),
        (0.25, 0.75),
        (0.5, 0.67),
        (0.34, 0.5),
    ] {
        assert_eq!(
            trimmed.evaluate(u, v),
            Some(trimmed.surface.evaluate(u, v)),
            "({u}, {v}) must be kept"
        );
    }
    assert!(!trimmed.is_inside(0.5, 0.5));
    assert!(trimmed.is_inside(0.3, 0.3));
}

#[test]
fn winding_numbers_follow_the_orientation_of_the_loops() {
    let trimmed = square_with_hole();
    assert_eq!(trimmed.winding_number(0.3, 0.3), 1);
    assert_eq!(trimmed.winding_number(0.5, 0.5), 0);
    assert_eq!(trimmed.winding_number(0.1, 0.1), 0);

    let mut clockwise_outer = NurbsTrimming::new(unit_surface());
    clockwise_outer.add_trim_curve(square_loop(0.2, 0.8).reversed());
    assert_eq!(clockwise_outer.winding_number(0.5, 0.5), -1);
    assert!(clockwise_outer.is_inside(0.5, 0.5));
    assert!(!clockwise_outer.is_inside(0.1, 0.5));

    // Both loops counter-clockwise: the non-zero winding rule keeps the "hole".
    let mut same_direction = NurbsTrimming::new(unit_surface());
    same_direction.add_trim_curve(square_loop(0.2, 0.8));
    same_direction.add_trim_curve(circle_loop((0.5, 0.5), 0.15));
    assert_eq!(same_direction.winding_number(0.5, 0.5), 2);
    assert!(same_direction.is_inside(0.5, 0.5));
}

#[test]
fn the_orientation_of_the_hole_decides_whether_it_is_cut() {
    let mut trimmed = NurbsTrimming::new(unit_surface());
    trimmed.add_trim_curve(square_loop(0.2, 0.8).reversed());
    trimmed.add_trim_curve(circle_loop((0.5, 0.5), 0.15));
    // Outer loop clockwise, hole counter-clockwise: windings -1 and +1 cancel.
    assert!(!trimmed.is_inside(0.5, 0.5));
    assert!(trimmed.is_inside(0.3, 0.3));
    assert!(!trimmed.is_inside(0.1, 0.1));
}

#[test]
fn separate_curves_are_chained_into_one_loop_in_any_order() {
    let mut trimmed = NurbsTrimming::new(unit_surface());
    // The four sides of [0.2, 0.8]^2, listed out of order.
    trimmed.add_trim_curve(polyline_curve(&[(0.8, 0.2), (0.8, 0.8)]));
    trimmed.add_trim_curve(polyline_curve(&[(0.2, 0.8), (0.2, 0.2)]));
    trimmed.add_trim_curve(polyline_curve(&[(0.2, 0.2), (0.8, 0.2)]));
    trimmed.add_trim_curve(polyline_curve(&[(0.8, 0.8), (0.2, 0.8)]));
    assert_eq!(trimmed.trim_loops().len(), 1);
    assert_eq!(trimmed.winding_number(0.5, 0.5), 1);
    assert!(trimmed.is_inside(0.5, 0.5));
    assert!(trimmed.is_inside(0.25, 0.75));
    assert!(!trimmed.is_inside(0.1, 0.5));
    assert!(!trimmed.is_inside(0.5, 0.9));
}

#[test]
fn a_loop_may_mix_curve_kinds() {
    // Lower half: a rational semicircle; upper part: three straight sides.
    let semicircle = {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsTrimCurve::new(
            vec![
                NurbsControlPoint::new(0.8, 0.5, 0.0, 1.0),
                NurbsControlPoint::new(0.8, 0.8, 0.0, w),
                NurbsControlPoint::new(0.5, 0.8, 0.0, 1.0),
                NurbsControlPoint::new(0.2, 0.8, 0.0, w),
                NurbsControlPoint::new(0.2, 0.5, 0.0, 1.0),
            ],
            NurbsKnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]),
            2,
        )
    };
    let mut trimmed = NurbsTrimming::new(unit_surface());
    trimmed.add_trim_curve(polyline_curve(&[
        (0.2, 0.5),
        (0.2, 0.2),
        (0.8, 0.2),
        (0.8, 0.5),
    ]));
    trimmed.add_trim_curve(semicircle);
    assert_eq!(trimmed.trim_loops().len(), 1);
    assert!(trimmed.is_inside(0.5, 0.5));
    assert!(trimmed.is_inside(0.5, 0.75));
    assert!(!trimmed.is_inside(0.5, 0.85));
    assert!(!trimmed.is_inside(0.25, 0.8));
    assert!(!trimmed.is_inside(0.5, 0.1));
}

#[test]
fn concave_loops_are_classified_correctly() {
    // L-shaped region, counter-clockwise.
    let mut trimmed = NurbsTrimming::new(unit_surface());
    trimmed.add_trim_curve(polyline_curve(&[
        (0.1, 0.1),
        (0.9, 0.1),
        (0.9, 0.4),
        (0.4, 0.4),
        (0.4, 0.9),
        (0.1, 0.9),
        (0.1, 0.1),
    ]));
    for (u, v) in [(0.2, 0.8), (0.5, 0.2), (0.2, 0.2), (0.8, 0.3), (0.3, 0.6)] {
        assert!(trimmed.is_inside(u, v), "({u}, {v}) is inside the L");
    }
    for (u, v) in [
        (0.7, 0.7),
        (0.95, 0.5),
        (0.5, 0.5),
        (0.05, 0.5),
        (0.5, 0.95),
    ] {
        assert!(!trimmed.is_inside(u, v), "({u}, {v}) is outside the L");
    }
}

#[test]
fn an_open_curve_is_closed_by_a_straight_segment() {
    // The lower-right triangle of [0.2, 0.8]^2 is bounded by two sides and the
    // diagonal, which closes the open curve.
    let mut trimmed = NurbsTrimming::new(unit_surface());
    trimmed.add_trim_curve(polyline_curve(&[(0.2, 0.2), (0.8, 0.2), (0.8, 0.8)]));
    assert_eq!(trimmed.trim_loops().len(), 1);
    assert!(trimmed.is_inside(0.6, 0.4));
    assert!(!trimmed.is_inside(0.4, 0.6));
    assert!(!trimmed.is_inside(0.1, 0.1));
}

#[test]
fn points_on_a_loop_belong_to_the_surface() {
    let trimmed = square_with_hole();
    // On the outer loop (edge and corner) and on the hole (vertex of the circle).
    assert!(trimmed.is_inside(0.2, 0.5));
    assert!(trimmed.is_inside(0.8, 0.8));
    assert!(trimmed.is_inside(0.5, 0.35));
    assert!(trimmed.is_inside(0.35, 0.5));
}

#[test]
fn an_untrimmed_surface_keeps_every_point() {
    let trimmed = NurbsTrimming::new(unit_surface());
    assert!(trimmed.trim_loops().is_empty());
    for (u, v) in [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (0.3, 0.9)] {
        assert_eq!(trimmed.evaluate(u, v), Some(trimmed.surface.evaluate(u, v)));
    }
}

#[test]
fn malformed_trim_curves_are_ignored() {
    let malformed = [
        NurbsTrimCurve::new(Vec::new(), NurbsKnotVector { knots: Vec::new() }, 2),
        NurbsTrimCurve::new(
            vec![trim_cp(0.0, 0.0)],
            NurbsKnotVector {
                knots: vec![0.0, 1.0],
            },
            0,
        ),
        NurbsTrimCurve::new(
            vec![trim_cp(0.0, 0.0), trim_cp(1.0, 0.0)],
            NurbsKnotVector {
                knots: vec![0.0, 1.0],
            },
            1,
        ),
        NurbsTrimCurve::new(
            vec![trim_cp(0.0, 0.0), trim_cp(f64::NAN, 0.0)],
            NurbsKnotVector {
                knots: vec![0.0, 0.0, 1.0, 1.0],
            },
            1,
        ),
        NurbsTrimCurve::new(
            vec![
                NurbsControlPoint::new(0.0, 0.0, 0.0, 0.0),
                trim_cp(1.0, 0.0),
            ],
            NurbsKnotVector {
                knots: vec![0.0, 0.0, 1.0, 1.0],
            },
            1,
        ),
    ];
    for curve in &malformed {
        assert!(!curve.is_valid());
        assert!(curve.domain().is_none());
        assert!(curve.flatten(1e-6).is_empty());
    }
    assert!(malformed[0].evaluate(0.5)[0].is_nan());

    // Only malformed curves: nothing is trimmed.
    let mut only_malformed = NurbsTrimming::new(unit_surface());
    for curve in &malformed {
        only_malformed.add_trim_curve(curve.clone());
    }
    assert!(only_malformed.trim_loops().is_empty());
    assert!(only_malformed.is_inside(0.5, 0.5));

    // Next to a valid loop they change nothing.
    let mut mixed = square_with_hole();
    for curve in &malformed {
        mixed.add_trim_curve(curve.clone());
    }
    assert!(!mixed.is_inside(0.5, 0.5));
    assert!(!mixed.is_inside(0.1, 0.1));
    assert!(mixed.is_inside(0.3, 0.3));
}

#[test]
fn the_loop_snapshot_agrees_with_the_direct_test() {
    let trimmed = square_with_hole();
    let loops = trimmed.trim_loops();
    assert_eq!(loops.len(), 2);
    let n = 80;
    let mut kept = 0usize;
    for i in 0..n {
        for j in 0..n {
            let (u, v) = ((i as f64 + 0.5) / n as f64, (j as f64 + 0.5) / n as f64);
            let inside = loops.contains(u, v);
            if i % 7 == 0 {
                assert_eq!(inside, trimmed.is_inside(u, v), "({u}, {v})");
            }
            assert_eq!(
                loops.winding_number(u, v) != 0,
                inside,
                "no sample lies on a loop: ({u}, {v})"
            );
            kept += usize::from(inside);
        }
    }
    // Area of the kept region: 0.6^2 - pi 0.15^2.
    let expected = 0.36 - std::f64::consts::PI * 0.15 * 0.15;
    let measured = kept as f64 / (n * n) as f64;
    assert!(
        (measured - expected).abs() < 0.01,
        "kept fraction {measured}, expected {expected}"
    );
}

#[test]
fn flattened_trim_curves_follow_the_curve_within_the_tolerance() {
    let circle = circle_loop((0.4, 0.6), 0.3);
    let tolerance = 1e-6;
    let polyline = circle.flatten(tolerance);
    assert!(polyline.len() > 100);
    let first = polyline[0];
    let last = polyline[polyline.len() - 1];
    assert!((first[0] - last[0]).abs() < 1e-12 && (first[1] - last[1]).abs() < 1e-12);
    let radius = |p: [f64; 2]| ((p[0] - 0.4).powi(2) + (p[1] - 0.6).powi(2)).sqrt();
    for p in &polyline {
        assert!((radius(*p) - 0.3).abs() < 1e-12, "vertex off the circle");
    }
    for pair in polyline.windows(2) {
        let mid = [
            0.5 * (pair[0][0] + pair[1][0]),
            0.5 * (pair[0][1] + pair[1][1]),
        ];
        let sagitta = 0.3 - radius(mid);
        assert!(sagitta < 4.0 * tolerance, "chord deviates by {sagitta}");
    }
}

#[test]
fn a_reversed_trim_curve_runs_backwards() {
    let circle = circle_loop((0.5, 0.5), 0.2);
    let reversed = circle.reversed();
    assert_eq!(reversed.degree, circle.degree);
    assert_eq!(reversed.domain(), circle.domain());
    for k in 0..=20 {
        let t = f64::from(k) / 20.0;
        let forward = circle.evaluate(t);
        let backward = reversed.evaluate(1.0 - t);
        assert!(
            (forward[0] - backward[0]).abs() < 1e-12 && (forward[1] - backward[1]).abs() < 1e-12,
            "t = {t}"
        );
    }
    let twice = reversed.reversed();
    assert_eq!(twice.knots.knots, circle.knots.knots);
    for (a, b) in twice.control_points.iter().zip(&circle.control_points) {
        assert_eq!(a.point, b.point);
        assert_eq!(a.weight, b.weight);
    }
    // A curve over a non-trivial parameter range mirrors within that range.
    let shifted = NurbsTrimCurve::new(
        vec![trim_cp(0.0, 0.0), trim_cp(1.0, 0.0), trim_cp(1.0, 1.0)],
        NurbsKnotVector::new(vec![2.0, 2.0, 5.0, 9.0, 9.0]),
        1,
    );
    let mirrored = shifted.reversed();
    assert_eq!(mirrored.knots.knots, vec![2.0, 2.0, 6.0, 9.0, 9.0]);
    let a = shifted.evaluate(3.0);
    let b = mirrored.evaluate(2.0 + 9.0 - 3.0);
    assert!((a[0] - b[0]).abs() < 1e-12 && (a[1] - b[1]).abs() < 1e-12);
}

#[test]
fn a_valid_trim_curve_reports_its_domain() {
    let circle = circle_loop((0.5, 0.5), 0.2);
    assert!(circle.is_valid());
    assert_eq!(circle.domain(), Some((0.0, 1.0)));
    let square = square_loop(0.0, 1.0);
    assert_eq!(square.domain(), Some((0.0, 4.0)));
}
