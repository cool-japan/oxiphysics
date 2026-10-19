// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the GJK distance solver, the simplex reduction, the EPA
//! contact points and the time of impact of [`crate::gjk_extended`].
//!
//! * `GjkSolver::query` reported distance `0` for any non-enclosing
//!   tetrahedron (`closest_to_origin` returned the origin for four vertices),
//!   compared against the oldest vertex and returned its witnesses.
//! * `EpaExpander::expand` fabricated normal `+z`, depth `0` and zero contact
//!   points, and used 1/3 centroid weights instead of the projection of the
//!   origin.
//! * `GjkTimeOfImpact::compute` ignored the relative velocity and never
//!   returned a time of impact.

use super::*;
use crate::gjk_epa::simplex_closest::closest_point_to_origin;
use crate::gjk_epa::test_support::{
    Rng, axis_separation, box_axes, box_corners, cloud_support, interval, polyhedron_distance,
    polyhedron_separation, rotation_about, translated,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

// ── shapes ───────────────────────────────────────────────────────────────────

fn sphere(center: [f64; 3], radius: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |d| {
        let n = norm3(d);
        let n = if len_sq3(n) > 0.0 { n } else { [1.0, 0.0, 0.0] };
        add3(center, scale3(n, radius))
    }
}

fn cuboid(center: [f64; 3], half: [f64; 3]) -> impl Fn([f64; 3]) -> [f64; 3] {
    let pick = |x: f64| if x >= 0.0 { 1.0 } else { -1.0 };
    move |d| {
        [
            center[0] + half[0] * pick(d[0]),
            center[1] + half[1] * pick(d[1]),
            center[2] + half[2] * pick(d[2]),
        ]
    }
}

fn hull(points: Vec<[f64; 3]>) -> impl Fn([f64; 3]) -> [f64; 3] {
    move |d| cloud_support(&points, d)
}

/// Support triple `(w, sa, sb)` of the Minkowski difference of two shapes.
fn difference<A, B>(a: A, b: B) -> impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3])
where
    A: Fn([f64; 3]) -> [f64; 3],
    B: Fn([f64; 3]) -> [f64; 3],
{
    move |d| {
        let sa = a(d);
        let sb = b(neg3(d));
        (sub3(sa, sb), sa, sb)
    }
}

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
    assert!(
        len3(sub3(actual, expected)) <= tol,
        "expected {expected:?}, got {actual:?} (tolerance {tol})"
    );
}

// ── GjkSolver::query ─────────────────────────────────────────────────────────

#[test]
fn boxes_with_a_lateral_offset_report_the_face_gap_and_face_witnesses() {
    // Boxes of half extent 1 at the origin and at (3, 0.5, 0): the faces x = 1 and
    // x = 2 face each other, so the distance is exactly 1. The old solver
    // returned the length of the oldest support point.
    let center_b = [3.0, 0.5, 0.0];
    let mut solver = GjkSolver::new(64, 1e-12);
    let result = solver.query(difference(
        cuboid([0.0; 3], [1.0; 3]),
        cuboid(center_b, [1.0; 3]),
    ));
    assert!(!result.intersect);
    assert!((result.distance - 1.0).abs() < 1e-9, "{}", result.distance);
    assert!(
        (result.point_a[0] - 1.0).abs() < 1e-9,
        "{:?}",
        result.point_a
    );
    assert!(
        (result.point_b[0] - 2.0).abs() < 1e-9,
        "{:?}",
        result.point_b
    );
    // The witnesses realise the distance and lie inside their boxes.
    assert_vec_close(sub3(result.point_a, result.point_b), [-1.0, 0.0, 0.0], 1e-9);
    for (k, &center) in center_b.iter().enumerate() {
        assert!(result.point_a[k].abs() <= 1.0 + 1e-9);
        assert!((result.point_b[k] - center).abs() <= 1.0 + 1e-9);
    }
}

#[test]
fn diagonal_box_pair_reports_the_corner_to_corner_distance() {
    // Gaps of 2, 3 and 11 along x, y and z: the closest features are two corners.
    let mut solver = GjkSolver::new(64, 1e-12);
    let result = solver.query(difference(
        cuboid([0.0; 3], [1.0; 3]),
        cuboid([4.0, 5.0, 13.0], [1.0; 3]),
    ));
    assert!(!result.intersect);
    let exact = 134.0_f64.sqrt();
    assert!(
        (result.distance - exact).abs() < 1e-9,
        "{}",
        result.distance
    );
    assert_vec_close(result.point_a, [1.0, 1.0, 1.0], 1e-9);
    assert_vec_close(result.point_b, [3.0, 4.0, 12.0], 1e-9);
}

#[test]
fn separated_spheres_report_the_surface_gap() {
    let mut solver = GjkSolver::new(64, 1e-12);
    let result = solver.query(difference(
        sphere([0.0; 3], 1.0),
        sphere([0.0, 5.0, 0.0], 1.0),
    ));
    assert!(!result.intersect);
    assert!((result.distance - 3.0).abs() < 1e-9, "{}", result.distance);
    // The witnesses realise the distance exactly; on a smooth shape their
    // position is accurate to about the square root of the tolerance.
    assert!((len3(sub3(result.point_a, result.point_b)) - result.distance).abs() < 1e-9);
    assert_vec_close(result.point_a, [0.0, 1.0, 0.0], 1e-5);
    assert_vec_close(result.point_b, [0.0, 4.0, 0.0], 1e-5);
}

#[test]
fn closest_to_origin_of_a_non_enclosing_tetrahedron_is_not_the_origin() {
    // Every point of this tetrahedron has x >= 2, so the origin is outside and
    // the nearest point is the vertex (2, 0, 0). The old code answered
    // `([0, 0, 0], false)` for any four vertices.
    let mut simplex = SimplexManager::new();
    simplex.add([2.0, 0.0, 0.0], [20.0, 0.0, 0.0], [18.0, 0.0, 0.0]);
    simplex.add([3.0, 0.0, 0.0], [30.0, 0.0, 0.0], [27.0, 0.0, 0.0]);
    simplex.add([2.5, 1.0, 0.0], [25.0, 1.0, 0.0], [22.5, 0.0, 0.0]);
    simplex.add([2.5, 0.5, 1.0], [25.0, 0.5, 1.0], [22.5, 0.0, 0.0]);
    let (closest, contains) = simplex.closest_to_origin();
    assert!(!contains);
    assert_vec_close(closest, [2.0, 0.0, 0.0], 1e-12);

    // The witnesses blend the support points of the same vertex.
    let (point_a, point_b) = simplex.witness_points();
    assert_vec_close(point_a, [20.0, 0.0, 0.0], 1e-12);
    assert_vec_close(point_b, [18.0, 0.0, 0.0], 1e-12);
    assert_vec_close(sub3(point_a, point_b), closest, 1e-12);

    // Reducing keeps exactly the supporting vertex, in all three arrays.
    simplex.reduce_to_closest();
    assert_eq!(simplex.size(), 1);
    assert_eq!(simplex.support_a, vec![[20.0, 0.0, 0.0]]);
    assert_eq!(simplex.support_b, vec![[18.0, 0.0, 0.0]]);
}

#[test]
fn closest_to_origin_over_a_face_interior_keeps_the_three_face_vertices() {
    // The origin projects to (0, 0, 1), inside the base triangle at z = 1.
    let mut simplex = SimplexManager::new();
    simplex.add([-1.0, -1.0, 1.0], [-1.0, -1.0, 11.0], [0.0, 0.0, 10.0]);
    simplex.add([2.0, -1.0, 1.0], [2.0, -1.0, 11.0], [0.0, 0.0, 10.0]);
    simplex.add([-1.0, 2.0, 1.0], [-1.0, 2.0, 11.0], [0.0, 0.0, 10.0]);
    simplex.add([0.0, 0.0, 3.0], [0.0, 0.0, 13.0], [0.0, 0.0, 10.0]);
    let (closest, contains) = simplex.closest_to_origin();
    assert!(!contains);
    assert_vec_close(closest, [0.0, 0.0, 1.0], 1e-12);
    let (point_a, point_b) = simplex.witness_points();
    assert_vec_close(sub3(point_a, point_b), closest, 1e-12);
    simplex.reduce_to_closest();
    assert_eq!(simplex.size(), 3);
    assert_eq!(simplex.support_a.len(), 3);
    assert_eq!(simplex.support_b.len(), 3);
    assert!(simplex.vertices.iter().all(|v| v[2] == 1.0));
}

#[test]
fn enclosing_tetrahedron_is_kept_whole_and_contains_the_origin() {
    let mut simplex = SimplexManager::new();
    simplex.add([1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [0.0; 3]);
    simplex.add([-1.0, -1.0, 1.0], [-1.0, -1.0, 1.0], [0.0; 3]);
    simplex.add([-1.0, 1.0, -1.0], [-1.0, 1.0, -1.0], [0.0; 3]);
    simplex.add([1.0, -1.0, -1.0], [1.0, -1.0, -1.0], [0.0; 3]);
    let (closest, contains) = simplex.closest_to_origin();
    assert!(contains);
    assert_vec_close(closest, [0.0; 3], 1e-12);
    simplex.reduce_to_closest();
    assert_eq!(simplex.size(), 4);
}

#[test]
fn distance_to_a_tetrahedron_matches_the_exact_distance() {
    // Shape A is a tetrahedron away from the origin, shape B a single point at
    // the origin: the distance is the distance from the origin to the tetrahedron.
    let tetra = vec![
        [-1.0, -1.0, 1.0],
        [2.0, -1.0, 1.0],
        [-1.0, 2.0, 1.0],
        [0.0, 0.0, 3.0],
    ];
    let mut solver = GjkSolver::new(64, 1e-12);
    let result = solver.query(difference(hull(tetra.clone()), hull(vec![[0.0; 3]])));
    assert!(!result.intersect);
    assert!((result.distance - 1.0).abs() < 1e-12, "{}", result.distance);
    assert_vec_close(sub3(result.point_a, result.point_b), [0.0, 0.0, 1.0], 1e-12);
    assert_vec_close(result.point_b, [0.0; 3], 1e-12);
    assert!((polyhedron_distance(&tetra, &[[0.0; 3]]) - 1.0).abs() < 1e-12);
}

#[test]
fn overlapping_spheres_report_a_shared_point_and_a_simplex_through_the_origin() {
    let mut solver = GjkSolver::new(64, 1e-8);
    let center_b = [1.1, 0.7, 0.4];
    let a = sphere([0.0; 3], 1.0);
    let b = sphere(center_b, 1.0);
    let result = solver.query(difference(&a, &b));
    assert!(result.intersect);
    assert_eq!(result.distance, 0.0);
    // Both witnesses are the same point, inside both spheres.
    assert_vec_close(result.point_a, result.point_b, 1e-12);
    assert!(len3(result.point_a) <= 1.0 + 1e-9);
    assert!(len3(sub3(result.point_a, center_b)) <= 1.0 + 1e-9);
    // The support points of two spheres stay in one plane through the origin,
    // so GJK ends on a triangle that contains the origin.
    let (closest, contains) = result.simplex.closest_to_origin();
    assert!(contains);
    assert_vec_close(closest, [0.0; 3], 1e-12);
}

#[test]
fn overlapping_polyhedra_end_on_a_tetrahedron_that_encloses_the_origin() {
    // Generic polyhedra leave every plane through the origin, so the simplex that
    // proves the overlap is a full tetrahedron, the seed of the EPA polytope.
    let mut rng = Rng(0xC0FF_EE00_1234_5678);
    let mut enclosed = 0;
    for case in 0..200 {
        let cloud_a: Vec<[f64; 3]> = (0..8).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..8).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            rng.vec(-0.4, 0.4),
        );
        if polyhedron_separation(&cloud_a, &cloud_b).gap > -0.05 {
            continue;
        }
        let result = GjkSolver::new(64, 1e-12)
            .query(difference(hull(cloud_a.clone()), hull(cloud_b.clone())));
        assert!(result.intersect, "case {case}");
        assert_eq!(result.simplex.size(), 4, "case {case}");
        let corners = [
            result.simplex.vertices[0],
            result.simplex.vertices[1],
            result.simplex.vertices[2],
            result.simplex.vertices[3],
        ];
        assert!(
            closest_point_to_origin(&corners).encloses_origin,
            "case {case}: the final simplex does not enclose the origin"
        );
        // The support points of the simplex stay in sync with its vertices.
        for index in 0..4 {
            assert_vec_close(
                sub3(
                    result.simplex.support_a[index],
                    result.simplex.support_b[index],
                ),
                result.simplex.vertices[index],
                1e-12,
            );
        }
        enclosed += 1;
    }
    assert!(enclosed > 100, "only {enclosed} decisive overlaps");
}

#[test]
fn symmetric_overlap_through_the_origin_still_reports_the_shared_point() {
    // Two spheres on the x axis: GJK ends on the segment of the two support
    // points, which passes through the origin, so the simplex has two vertices.
    let mut solver = GjkSolver::new(64, 1e-8);
    let result = solver.query(difference(
        sphere([0.0; 3], 1.0),
        sphere([1.5, 0.0, 0.0], 1.0),
    ));
    assert!(result.intersect);
    assert_eq!(result.distance, 0.0);
    assert_vec_close(result.point_a, result.point_b, 1e-12);
    // The shared point lies in the lens of the two spheres.
    assert!(len3(result.point_a) <= 1.0 + 1e-9);
    assert!(len3(sub3(result.point_a, [1.5, 0.0, 0.0])) <= 1.0 + 1e-9);
    assert!(result.simplex.size() >= 2);
}

#[test]
fn separated_random_polyhedra_match_the_exact_distance_and_overlap_verdict() {
    let mut rng = Rng(0x5EED_0001_C0DE_F00D);
    let mut separated = 0;
    let mut overlapping = 0;
    for case in 0..300 {
        let count_a = 5 + (rng.unit() * 4.0) as usize;
        let count_b = 5 + (rng.unit() * 4.0) as usize;
        let cloud_a: Vec<[f64; 3]> = (0..count_a).map(|_| rng.vec(-1.0, 1.0)).collect();
        let offset = rng.vec(-1.6, 1.6);
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..count_b).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            offset,
        );
        let gap = polyhedron_separation(&cloud_a, &cloud_b).gap;
        if gap.abs() < 1e-6 {
            continue; // either classification is acceptable
        }
        let mut solver = GjkSolver::new(128, 1e-12);
        let result = solver.query(difference(hull(cloud_a.clone()), hull(cloud_b.clone())));
        if gap > 0.0 {
            separated += 1;
            assert!(
                !result.intersect,
                "case {case}: gap {gap} reported as overlap"
            );
            let exact = polyhedron_distance(&cloud_a, &cloud_b);
            assert!(
                (result.distance - exact).abs() < 1e-9,
                "case {case}: distance {} vs exact {exact}",
                result.distance
            );
            assert!(
                (len3(sub3(result.point_a, result.point_b)) - result.distance).abs() < 1e-9,
                "case {case}: witnesses do not realise the distance"
            );
            // The witnesses lie inside the bounding intervals of their clouds.
            for k in 0..3 {
                let mut axis = [0.0; 3];
                axis[k] = 1.0;
                let (a_lo, a_hi) = interval(&cloud_a, axis);
                let (b_lo, b_hi) = interval(&cloud_b, axis);
                assert!(result.point_a[k] >= a_lo - 1e-9 && result.point_a[k] <= a_hi + 1e-9);
                assert!(result.point_b[k] >= b_lo - 1e-9 && result.point_b[k] <= b_hi + 1e-9);
            }
        } else {
            overlapping += 1;
            assert!(
                result.intersect,
                "case {case}: overlap {gap} reported as separated"
            );
            assert_eq!(result.distance, 0.0);
        }
    }
    assert!(separated > 40, "only {separated} separated configurations");
    assert!(
        overlapping > 40,
        "only {overlapping} overlapping configurations"
    );
}

#[test]
fn cache_keeps_the_separating_axis_and_is_reset_by_an_overlap() {
    let mut solver = GjkSolver::new(64, 1e-12);
    let a = sphere([0.0; 3], 1.0);
    let far = solver.query(difference(&a, sphere([0.0, 5.0, 0.0], 1.0)));
    assert!(!far.intersect);
    assert!(solver.cache.valid);
    // The first search direction of the next query is the separating axis.
    assert_vec_close(solver.cache.last_dir, [0.0, 1.0, 0.0], 1e-6);

    // A warm-started query agrees with a cold one.
    let warm = solver.query(difference(&a, sphere([0.0, 5.5, 0.0], 1.0)));
    let cold = GjkSolver::new(64, 1e-12).query(difference(&a, sphere([0.0, 5.5, 0.0], 1.0)));
    assert!((warm.distance - 3.5).abs() < 1e-9);
    assert!((warm.distance - cold.distance).abs() < 1e-9);

    let near = solver.query(difference(&a, sphere([0.0, 1.0, 0.0], 1.0)));
    assert!(near.intersect);
    assert!(!solver.cache.valid);
}

#[test]
fn iteration_limit_reports_an_upper_bound_never_a_false_contact() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([6.0, 4.0, 0.0], [1.0; 3]);
    let exact = 20.0_f64.sqrt();
    for max_iter in [0, 1, 2, 3] {
        let result = GjkSolver::new(max_iter, 1e-12).query(difference(&a, &b));
        assert!(!result.intersect, "max_iter {max_iter}");
        assert!(result.distance > 0.0, "max_iter {max_iter}");
        assert!(
            result.distance >= exact - 1e-9,
            "max_iter {max_iter}: {} is below the true distance {exact}",
            result.distance
        );
    }
}

// ── EpaExpander contact points ───────────────────────────────────────────────

/// Alternating cube corners: a tetrahedron that encloses the origin for unit
/// cubes that overlap along x.
const CUBE_SEED: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
    [1.0, -1.0, -1.0],
];

/// Seed tetrahedron of the Minkowski difference along four directions.
///
/// The directions must be chosen so that the tetrahedron encloses the origin.
fn polytope_of<F>(support: &F, directions: [[f64; 3]; 4]) -> ExpandingPolytope
where
    F: Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]),
{
    let mut simplex = SimplexManager::new();
    for dir in directions {
        let (v, sa, sb) = support(dir);
        simplex.add(v, sa, sb);
    }
    let corners = [
        simplex.vertices[0],
        simplex.vertices[1],
        simplex.vertices[2],
        simplex.vertices[3],
    ];
    assert!(
        closest_point_to_origin(&corners).encloses_origin,
        "the seed tetrahedron must enclose the origin"
    );
    ExpandingPolytope::from_simplex(&simplex)
}

#[test]
fn epa_expander_reports_the_contact_points_of_overlapping_boxes() {
    // Unit cubes (half extent 0.5) overlapping by 0.3 along x.
    let support = difference(
        cuboid([0.0; 3], [0.5; 3]),
        cuboid([0.7, 0.0, 0.0], [0.5; 3]),
    );
    let mut polytope = polytope_of(&support, CUBE_SEED);
    let result = EpaExpander::new(64, 1e-12)
        .expand(&mut polytope, &support)
        .expect("a seeded polytope has faces");
    assert!((result.depth - 0.3).abs() < 1e-9, "depth {}", result.depth);
    // Normal from B to A: A is on the -x side of B.
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-9);
    // Deepest points: A's +x face and B's -x face, one 0.3 apart.
    assert!(
        (result.point_a[0] - 0.5).abs() < 1e-9,
        "{:?}",
        result.point_a
    );
    assert!(
        (result.point_b[0] - 0.2).abs() < 1e-9,
        "{:?}",
        result.point_b
    );
    assert_vec_close(
        sub3(result.point_a, result.point_b),
        scale3(result.normal, -result.depth),
        1e-9,
    );
    for k in 0..3 {
        assert!(
            result.point_a[k].abs() <= 0.5 + 1e-9,
            "{:?}",
            result.point_a
        );
        assert!(
            (result.point_b[k] - [0.7, 0.0, 0.0][k]).abs() <= 0.5 + 1e-9,
            "{:?}",
            result.point_b
        );
    }
}

#[test]
fn epa_expander_contact_points_of_overlapping_spheres() {
    // Unit spheres whose centres are 1.5 apart overlap by 0.5.
    let support = difference(sphere([0.0; 3], 1.0), sphere([1.5, 0.0, 0.0], 1.0));
    // The difference is a ball of radius 2 around (-1.5, 0, 0): an apex beyond
    // the origin on +x and three base vertices around the -x side enclose it.
    let mut polytope = polytope_of(
        &support,
        [
            [1.0, 0.0, 0.0],
            [-1.0, 2.0, 0.0],
            [-1.0, -1.0, 1.732],
            [-1.0, -1.0, -1.732],
        ],
    );
    let result = EpaExpander::new(256, 1e-9)
        .expand(&mut polytope, &support)
        .expect("a seeded polytope has faces");
    assert!((result.depth - 0.5).abs() < 1e-4, "depth {}", result.depth);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert!((len3(result.normal) - 1.0).abs() < 1e-12);
    // The contact points differ by exactly the penetration vector ...
    assert_vec_close(
        sub3(result.point_a, result.point_b),
        scale3(result.normal, -result.depth),
        1e-9,
    );
    // ... and lie on (just inside) the surface of their sphere: they are blends
    // of surface points, the polytope is a close approximation of the sphere.
    let on_a = len3(result.point_a);
    let on_b = len3(sub3(result.point_b, [1.5, 0.0, 0.0]));
    assert!(on_a <= 1.0 + 1e-9 && on_a > 1.0 - 1e-3, "|a| = {on_a}");
    assert!(on_b <= 1.0 + 1e-9 && on_b > 1.0 - 1e-3, "|b| = {on_b}");
}

#[test]
fn epa_expander_contacts_match_the_exact_polyhedron_contacts() {
    // Rotated cuboids: the depth, the normal and the contact points must match
    // the separating-axis solution (the deepest vertex of A along the normal).
    let mut rng = Rng(0xE9A0_C0DE_0000_0001);
    let mut verified = 0;
    for case in 0..200 {
        let rot_a = rotation_about(rng.vec(-1.0, 1.0), rng.range(0.0, 3.0));
        let rot_b = rotation_about(rng.vec(-1.0, 1.0), rng.range(0.0, 3.0));
        let half_a = [
            rng.range(0.4, 1.0),
            rng.range(0.4, 1.0),
            rng.range(0.4, 1.0),
        ];
        let half_b = [
            rng.range(0.4, 1.0),
            rng.range(0.4, 1.0),
            rng.range(0.4, 1.0),
        ];
        let cloud_a = box_corners([0.0; 3], half_a, rot_a);
        let cloud_b = box_corners(rng.vec(-0.5, 0.5), half_b, rot_b);
        let separation = polyhedron_separation(&cloud_a, &cloud_b);
        if separation.gap > -0.05 {
            continue; // disjoint or too shallow to seed EPA robustly
        }
        let support = difference(hull(cloud_a.clone()), hull(cloud_b.clone()));
        let seed: Vec<_> = (0..200)
            .map(|_| {
                let supports: Vec<_> = (0..4).map(|_| support(rng.vec(-1.0, 1.0))).collect();
                supports
            })
            .find(|supports| {
                let points = [supports[0].0, supports[1].0, supports[2].0, supports[3].0];
                closest_point_to_origin(&points).encloses_origin
            })
            .unwrap_or_default();
        if seed.len() != 4 {
            continue;
        }
        let mut simplex = SimplexManager::new();
        for (v, sa, sb) in seed {
            simplex.add(v, sa, sb);
        }
        let mut polytope = ExpandingPolytope::from_simplex(&simplex);
        if polytope.closest_face().is_none_or(|face| face.dist < 1e-3) {
            continue;
        }
        let result = EpaExpander::new(256, 1e-10)
            .expand(&mut polytope, &support)
            .expect("a seeded polytope has faces");
        assert!(
            (result.depth + separation.gap).abs() < 1e-6,
            "case {case}: depth {} vs exact {}",
            result.depth,
            -separation.gap
        );
        assert!(
            dot3(result.normal, separation.from_b_to_a) > 1.0 - 1e-6,
            "case {case}: normal {:?} vs exact {:?}",
            result.normal,
            separation.from_b_to_a
        );
        assert_vec_close(
            sub3(result.point_a, result.point_b),
            scale3(result.normal, -result.depth),
            1e-8,
        );
        // The contact point on A lies on the face of A that is furthest along
        // -normal (the side that has to travel the most), B's on the opposite side.
        let toward_b = scale3(result.normal, -1.0);
        let (_, a_hi) = interval(&cloud_a, toward_b);
        let (b_lo, _) = interval(&cloud_b, toward_b);
        assert!(
            (dot3(result.point_a, toward_b) - a_hi).abs() < 1e-6,
            "case {case}: contact on A is not on its deepest face"
        );
        assert!(
            (dot3(result.point_b, toward_b) - b_lo).abs() < 1e-6,
            "case {case}: contact on B is not on its deepest face"
        );
        verified += 1;
    }
    assert!(
        verified > 40,
        "only {verified} configurations were decisive"
    );
}

#[test]
fn gjk_simplex_seeds_the_expander_for_random_polyhedra() {
    // The public chain GjkSolver::query -> ExpandingPolytope::from_simplex ->
    // EpaExpander::expand, with no hand-picked seed, must give the penetration
    // of the separating-axis solution, contact points included.
    let mut rng = Rng(0x5EED_C4A1_0000_0042);
    let mut verified = 0;
    for case in 0..300 {
        let cloud_a: Vec<[f64; 3]> = (0..8).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..8).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            rng.vec(-0.5, 0.5),
        );
        let separation = polyhedron_separation(&cloud_a, &cloud_b);
        if separation.gap > -0.05 {
            continue;
        }
        let support = difference(hull(cloud_a.clone()), hull(cloud_b.clone()));
        let result = GjkSolver::new(64, 1e-12).query(&support);
        assert!(result.intersect, "case {case}");
        let mut polytope = ExpandingPolytope::from_simplex(&result.simplex);
        // The faces of a polytope are oriented away from the origin, which needs
        // the origin strictly inside: skip the rare simplex that has it within
        // rounding distance of a face.
        if polytope.closest_face().is_none_or(|face| face.dist < 1e-6) {
            continue;
        }
        let penetration = EpaExpander::new(256, 1e-10)
            .expand(&mut polytope, &support)
            .unwrap_or_else(|| panic!("case {case}: the GJK simplex gave no penetration"));
        assert!(
            (penetration.depth + separation.gap).abs() < 1e-6,
            "case {case}: depth {} vs exact {}",
            penetration.depth,
            -separation.gap
        );
        assert!(
            dot3(penetration.normal, separation.from_b_to_a) > 1.0 - 1e-6,
            "case {case}: normal {:?} vs exact {:?}",
            penetration.normal,
            separation.from_b_to_a
        );
        assert_vec_close(
            sub3(penetration.point_a, penetration.point_b),
            scale3(penetration.normal, -penetration.depth),
            1e-8,
        );
        let toward_b = scale3(penetration.normal, -1.0);
        let (_, a_hi) = interval(&cloud_a, toward_b);
        let (b_lo, _) = interval(&cloud_b, toward_b);
        assert!(
            (dot3(penetration.point_a, toward_b) - a_hi).abs() < 1e-6,
            "case {case}: contact on A is not on its deepest face"
        );
        assert!(
            (dot3(penetration.point_b, toward_b) - b_lo).abs() < 1e-6,
            "case {case}: contact on B is not on its deepest face"
        );
        verified += 1;
    }
    assert!(verified > 60, "only {verified} decisive overlaps");
}

#[test]
fn epa_expander_without_faces_has_no_result() {
    // Four coincident vertices span no face: nothing can be derived, and the old
    // code answered with a fabricated +z normal and zero depth.
    let mut simplex = SimplexManager::new();
    for _ in 0..4 {
        simplex.add([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]);
    }
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    assert!(polytope.faces.is_empty());
    let result = EpaExpander::new(16, 1e-9).expand(&mut polytope, |_| {
        ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3])
    });
    assert!(result.is_none());
}

#[test]
fn short_simplex_of_overlapping_spheres_still_yields_the_penetration() {
    // Two spheres on the x axis: GJK ends on a segment through the origin, which
    // spans no volume. Seeding the polytope used to panic, and a polytope without
    // faces must not mean "no penetration" either: the expander completes the
    // simplex from the support function. Unit spheres whose centres are 1.5 apart
    // overlap by 0.5.
    let support = difference(sphere([0.0; 3], 1.0), sphere([1.5, 0.0, 0.0], 1.0));
    let result = GjkSolver::new(64, 1e-8).query(&support);
    assert!(result.intersect);
    assert!(result.simplex.size() < 4);
    let mut polytope = ExpandingPolytope::from_simplex(&result.simplex);
    assert!(polytope.faces.is_empty());
    let penetration = EpaExpander::new(256, 1e-9)
        .expand(&mut polytope, &support)
        .expect("overlapping spheres penetrate");
    assert!(
        (penetration.depth - 0.5).abs() < 1e-4,
        "depth {}",
        penetration.depth
    );
    assert_vec_close(penetration.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert_vec_close(
        sub3(penetration.point_a, penetration.point_b),
        scale3(penetration.normal, -penetration.depth),
        1e-9,
    );
    let on_a = len3(penetration.point_a);
    let on_b = len3(sub3(penetration.point_b, [1.5, 0.0, 0.0]));
    assert!(on_a <= 1.0 + 1e-9 && on_a > 1.0 - 1e-3, "|a| = {on_a}");
    assert!(on_b <= 1.0 + 1e-9 && on_b > 1.0 - 1e-3, "|b| = {on_b}");
    // The polytope handed in is left as it was.
    assert!(polytope.faces.is_empty());
}

#[test]
fn query_to_polytope_to_expander_matches_the_sphere_penetration() {
    // The whole public path (GJK, polytope from its simplex, EPA) for spheres of
    // different radii whose centres differ in a generic direction. Every support
    // point of a sphere pair stays in the plane through the origin spanned by the
    // first search direction and the line of centres, so GJK ends on a flat
    // simplex there and the polytope has no faces: the expander has to complete
    // it from the support function.
    let cases = [
        ([0.9, 0.6, -0.4], 1.0, 1.0),
        ([0.3, -1.1, 0.5], 1.2, 0.8),
        ([-0.7, 0.2, 1.0], 0.9, 1.4),
        ([0.0, 0.0, 1.3], 1.0, 0.75),
        ([-1.2, -0.5, 0.25], 0.8, 0.8),
    ];
    for (center_b, radius_a, radius_b) in cases {
        let distance = len3(center_b);
        let overlap = radius_a + radius_b - distance;
        assert!(overlap > 0.1, "test geometry must overlap");
        let support = difference(sphere([0.0; 3], radius_a), sphere(center_b, radius_b));
        let result = GjkSolver::new(64, 1e-9).query(&support);
        assert!(result.intersect, "{center_b:?}");
        let mut polytope = ExpandingPolytope::from_simplex(&result.simplex);
        let penetration = EpaExpander::new(512, 1e-10)
            .expand(&mut polytope, &support)
            .expect("overlapping spheres penetrate");
        assert!(
            (penetration.depth - overlap).abs() < 1e-4,
            "{center_b:?}: depth {} vs {overlap}",
            penetration.depth
        );
        // The normal points from B to A: against the line of centres.
        let toward_b = scale3(center_b, 1.0 / distance);
        assert_vec_close(penetration.normal, neg3(toward_b), 2e-3);
        assert_vec_close(
            sub3(penetration.point_a, penetration.point_b),
            scale3(penetration.normal, -penetration.depth),
            1e-9,
        );
        // Deepest points: A towards B, B towards A, on the surfaces up to the
        // resolution of the polytope.
        let on_a = len3(penetration.point_a);
        let on_b = len3(sub3(penetration.point_b, center_b));
        assert!(
            on_a <= radius_a + 1e-9 && on_a > radius_a - 2e-3,
            "{center_b:?}: |a| = {on_a}"
        );
        assert!(
            on_b <= radius_b + 1e-9 && on_b > radius_b - 2e-3,
            "{center_b:?}: |b - c| = {on_b}"
        );
        assert!(dot3(penetration.point_a, toward_b) > 0.0, "{center_b:?}");
    }
}

#[test]
fn faceless_polytope_of_a_flat_difference_has_no_penetration() {
    // A flat square against a point inside it: the Minkowski difference is a
    // square in the plane z = 0 that contains the origin. It has no volume to
    // penetrate, so there is nothing to derive (and nothing is fabricated).
    let square = vec![
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ];
    let support = difference(hull(square), hull(vec![[0.2, 0.1, 0.0]]));
    let result = GjkSolver::new(64, 1e-9).query(&support);
    assert!(result.intersect);
    let mut polytope = ExpandingPolytope::from_simplex(&result.simplex);
    assert!(polytope.faces.is_empty());
    assert!(
        EpaExpander::new(64, 1e-9)
            .expand(&mut polytope, &support)
            .is_none()
    );
}

#[test]
fn faceless_polytope_of_separated_shapes_has_no_penetration() {
    // Seeding from vertices must not turn a separation into a penetration.
    let support = difference(sphere([0.0; 3], 1.0), sphere([5.0, 0.0, 0.0], 1.0));
    let mut simplex = SimplexManager::new();
    let (w, sa, sb) = support([1.0, 0.0, 0.0]);
    simplex.add(w, sa, sb);
    let mut polytope = ExpandingPolytope::from_simplex(&simplex);
    assert!(polytope.faces.is_empty());
    assert!(
        EpaExpander::new(64, 1e-9)
            .expand(&mut polytope, &support)
            .is_none()
    );
}

#[test]
fn epa_expander_iteration_limit_reports_the_nearest_face() {
    let support = difference(
        cuboid([0.0; 3], [0.5; 3]),
        cuboid([0.7, 0.0, 0.0], [0.5; 3]),
    );
    let mut polytope = polytope_of(&support, CUBE_SEED);
    let nearest = polytope.closest_face().expect("faces").dist;
    let result = EpaExpander::new(0, 1e-9)
        .expand(&mut polytope, &support)
        .expect("a seeded polytope has faces");
    assert!((result.depth - nearest).abs() < 1e-15);
    assert!((len3(result.normal) - 1.0).abs() < 1e-12);
}

// ── GjkTimeOfImpact ──────────────────────────────────────────────────────────

/// Sphere A at rest at the origin, sphere B starting at `start` and moving with
/// `velocity`.
fn moving_spheres(
    radius_a: f64,
    radius_b: f64,
    start: [f64; 3],
    velocity: [f64; 3],
) -> impl Fn([f64; 3], f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    move |d, t| {
        let sa = sphere([0.0; 3], radius_a)(d);
        let sb = sphere(add3(start, scale3(velocity, t)), radius_b)(neg3(d));
        (sub3(sa, sb), sa, sb)
    }
}

#[test]
fn unit_spheres_closing_at_ten_hit_at_the_analytic_time() {
    // 5 apart, radii 1 + 1: a gap of 3 closed at 10 per second takes 0.3 s.
    let toi = GjkTimeOfImpact::new(32, 1e-9);
    let result = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [-10.0, 0.0, 0.0]),
        [-10.0, 0.0, 0.0],
        1.0,
    );
    let t = result.toi.expect("the spheres collide");
    assert!((t - 0.3).abs() < 1e-6, "toi {t}");
    // Contact on A's surface facing B, normal from A toward B.
    assert_vec_close(result.contact_point, [1.0, 0.0, 0.0], 1e-4);
    assert_vec_close(result.contact_normal, [1.0, 0.0, 0.0], 1e-6);
}

#[test]
fn time_of_impact_uses_the_relative_velocity() {
    // The same pair with the velocity ignored would never collide: the support
    // closure is evaluated at the time the solver chooses.
    let toi = GjkTimeOfImpact::new(32, 1e-9);
    let slow = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [-2.0, 0.0, 0.0]),
        [-2.0, 0.0, 0.0],
        1.0,
    );
    assert!(slow.toi.is_none(), "{slow:?}");
    let slow_long = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [-2.0, 0.0, 0.0]),
        [-2.0, 0.0, 0.0],
        2.0,
    );
    let t = slow_long.toi.expect("collides within two seconds");
    assert!((t - 1.5).abs() < 1e-6, "toi {t}");
    // A zero relative velocity never closes the gap.
    let still = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [0.0; 3]),
        [0.0; 3],
        10.0,
    );
    assert!(still.toi.is_none());
    assert_eq!(still.contact_normal, [0.0; 3]);
}

#[test]
fn receding_and_passing_shapes_have_no_time_of_impact() {
    let toi = GjkTimeOfImpact::new(64, 1e-9);
    let receding = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]),
        [10.0, 0.0, 0.0],
        1.0,
    );
    assert!(receding.toi.is_none());
    // B flies past A at a lateral distance of 3 (gap 1 at the closest approach).
    let passing = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 3.0, 0.0], [-10.0, 0.0, 0.0]),
        [-10.0, 0.0, 0.0],
        1.0,
    );
    assert!(passing.toi.is_none());
}

#[test]
fn overlapping_or_touching_shapes_hit_at_time_zero() {
    let toi = GjkTimeOfImpact::new(32, 1e-9);
    let overlapping = toi.compute(
        moving_spheres(1.0, 1.0, [1.5, 0.0, 0.0], [-10.0, 0.0, 0.0]),
        [-10.0, 0.0, 0.0],
        1.0,
    );
    assert_eq!(overlapping.toi, Some(0.0));
    // No separating direction exists: B has to be pushed back along its motion.
    assert_vec_close(overlapping.contact_normal, [1.0, 0.0, 0.0], 1e-12);
    let touching = toi.compute(
        moving_spheres(1.0, 1.0, [2.0, 0.0, 0.0], [-10.0, 0.0, 0.0]),
        [-10.0, 0.0, 0.0],
        1.0,
    );
    assert_eq!(touching.toi, Some(0.0));
}

#[test]
fn invalid_time_steps_have_no_time_of_impact() {
    let toi = GjkTimeOfImpact::new(32, 1e-9);
    for dt in [f64::NAN, f64::INFINITY, -1.0] {
        let result = toi.compute(
            moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [-10.0, 0.0, 0.0]),
            [-10.0, 0.0, 0.0],
            dt,
        );
        assert!(result.toi.is_none(), "dt {dt}");
    }
    let nan_velocity = toi.compute(
        moving_spheres(1.0, 1.0, [5.0, 0.0, 0.0], [-10.0, 0.0, 0.0]),
        [f64::NAN, 0.0, 0.0],
        1.0,
    );
    assert!(nan_velocity.toi.is_none());
    // A zero time step still detects shapes that touch right away.
    let instant = toi.compute(
        moving_spheres(1.0, 1.0, [1.0, 0.0, 0.0], [-10.0, 0.0, 0.0]),
        [-10.0, 0.0, 0.0],
        0.0,
    );
    assert_eq!(instant.toi, Some(0.0));
}

/// Earliest `t >= 0` with `|p + v t| = r` (the entry time), if any.
fn sphere_entry_time(p: [f64; 3], v: [f64; 3], r: f64) -> Option<f64> {
    let a = dot3(v, v);
    let b = 2.0 * dot3(p, v);
    let c = dot3(p, p) - r * r;
    if c <= 0.0 {
        return Some(0.0);
    }
    if a < 1e-300 {
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let t = (-b - disc.sqrt()) / (2.0 * a);
    (t >= 0.0).then_some(t)
}

#[test]
fn random_sphere_sweeps_match_the_analytic_entry_time() {
    let mut rng = Rng(0x70A1_0F1E_5EED_0005);
    let toi = GjkTimeOfImpact::new(64, 1e-10);
    let (mut hits, mut misses) = (0, 0);
    for case in 0..400 {
        let radius_a = rng.range(0.3, 1.5);
        let radius_b = rng.range(0.3, 1.5);
        let start = rng.vec(-8.0, 8.0);
        // Aim roughly at the target so that a good share of the sweeps collide.
        let velocity = add3(scale3(start, -rng.range(0.4, 2.0)), rng.vec(-3.0, 3.0));
        let dt = rng.range(0.2, 1.5);
        let exact = sphere_entry_time(start, velocity, radius_a + radius_b);
        let result = toi.compute(
            moving_spheres(radius_a, radius_b, start, velocity),
            velocity,
            dt,
        );
        match exact {
            Some(t) if t <= dt - 1e-6 => {
                hits += 1;
                let found = result
                    .toi
                    .unwrap_or_else(|| panic!("case {case}: missed a hit at {t}"));
                assert!(
                    (found - t).abs() < 1e-7,
                    "case {case}: toi {found} vs exact {t}"
                );
                // The reported normal points from A toward B at that time (shapes
                // that overlap at the start have no separating direction).
                if t > 1e-9 {
                    let center_b = add3(start, scale3(velocity, found));
                    assert_vec_close(result.contact_normal, norm3(center_b), 1e-4);
                }
            }
            Some(t) if t >= dt + 1e-6 => {
                misses += 1;
                assert!(
                    result.toi.is_none(),
                    "case {case}: hit at {:?}, exact {t}",
                    result.toi
                );
            }
            None => {
                misses += 1;
                assert!(result.toi.is_none(), "case {case}: hit at {:?}", result.toi);
            }
            Some(_) => {}
        }
    }
    assert!(hits > 60 && misses > 60, "hits {hits}, misses {misses}");
}

/// Entry time of the point `p + v t` into the box `|x_k| <= half_k` (slab method).
fn box_entry_time(p: [f64; 3], v: [f64; 3], half: [f64; 3]) -> Option<f64> {
    let (mut near, mut far) = (0.0_f64, f64::INFINITY);
    for k in 0..3 {
        if v[k].abs() < 1e-300 {
            if p[k].abs() > half[k] {
                return None;
            }
            continue;
        }
        let t1 = (-half[k] - p[k]) / v[k];
        let t2 = (half[k] - p[k]) / v[k];
        near = near.max(t1.min(t2));
        far = far.min(t1.max(t2));
    }
    (near <= far).then_some(near)
}

#[test]
fn random_box_sweeps_match_the_slab_entry_time() {
    let mut rng = Rng(0xB0C5_5EED_0000_0007);
    let toi = GjkTimeOfImpact::new(64, 1e-10);
    let (mut hits, mut misses) = (0, 0);
    for case in 0..400 {
        let half_a = [
            rng.range(0.3, 1.2),
            rng.range(0.3, 1.2),
            rng.range(0.3, 1.2),
        ];
        let half_b = [
            rng.range(0.3, 1.2),
            rng.range(0.3, 1.2),
            rng.range(0.3, 1.2),
        ];
        let start = rng.vec(-6.0, 6.0);
        let velocity = add3(scale3(start, -rng.range(0.4, 2.0)), rng.vec(-3.0, 3.0));
        let dt = rng.range(0.2, 1.5);
        // The Minkowski sum of two boxes is the box of summed half extents.
        let summed = [
            half_a[0] + half_b[0],
            half_a[1] + half_b[1],
            half_a[2] + half_b[2],
        ];
        let exact = box_entry_time(start, velocity, summed);
        let support = move |d: [f64; 3], t: f64| {
            let sa = cuboid([0.0; 3], half_a)(d);
            let sb = cuboid(add3(start, scale3(velocity, t)), half_b)(neg3(d));
            (sub3(sa, sb), sa, sb)
        };
        let result = toi.compute(support, velocity, dt);
        match exact {
            Some(t) if t <= dt - 1e-6 => {
                hits += 1;
                let found = result
                    .toi
                    .unwrap_or_else(|| panic!("case {case}: missed a hit at {t}"));
                assert!(
                    (found - t).abs() < 1e-7,
                    "case {case}: toi {found} vs exact {t}"
                );
            }
            Some(t) if t >= dt + 1e-6 => {
                misses += 1;
                assert!(
                    result.toi.is_none(),
                    "case {case}: hit at {:?}, exact {t}",
                    result.toi
                );
            }
            None => {
                misses += 1;
                assert!(result.toi.is_none(), "case {case}: hit at {:?}", result.toi);
            }
            Some(_) => {}
        }
    }
    assert!(hits > 60 && misses > 60, "hits {hits}, misses {misses}");
}

#[test]
fn rotated_polyhedra_sweep_never_tunnels_and_stops_before_contact_overshoots() {
    // Moving rotated cuboids against a reference found by sampling the exact
    // separating-axis gap along the sweep and bisecting the first crossing.
    let mut rng = Rng(0xFEED_FACE_0000_0042);
    let toi = GjkTimeOfImpact::new(96, 1e-9);
    let mut checked = 0;
    for case in 0..120 {
        let rot_a = rotation_about(rng.vec(-1.0, 1.0), rng.range(0.0, 3.0));
        let rot_b = rotation_about(rng.vec(-1.0, 1.0), rng.range(0.0, 3.0));
        let axes = box_axes(rot_a, rot_b);
        let cloud_a = box_corners(
            [0.0; 3],
            [
                rng.range(0.4, 1.0),
                rng.range(0.4, 1.0),
                rng.range(0.4, 1.0),
            ],
            rot_a,
        );
        let cloud_b = box_corners(
            [0.0; 3],
            [
                rng.range(0.4, 1.0),
                rng.range(0.4, 1.0),
                rng.range(0.4, 1.0),
            ],
            rot_b,
        );
        let start = rng.vec(-5.0, 5.0);
        let velocity = add3(scale3(start, -rng.range(0.4, 2.0)), rng.vec(-2.0, 2.0));
        let dt = 1.0;
        let gap_at = |t: f64| {
            axis_separation(
                &cloud_a,
                &translated(&cloud_b, add3(start, scale3(velocity, t))),
                &axes,
            )
            .gap
        };
        if gap_at(0.0) <= 0.0 {
            continue;
        }
        let steps = 800;
        let mut first = None;
        let mut min_gap = f64::INFINITY;
        for step in 1..=steps {
            let t = dt * step as f64 / steps as f64;
            let gap = gap_at(t);
            min_gap = min_gap.min(gap);
            if gap <= 0.0 {
                first = Some(step);
                break;
            }
        }
        let support = |d: [f64; 3], t: f64| {
            let sa = cloud_support(&cloud_a, d);
            let moved = translated(&cloud_b, add3(start, scale3(velocity, t)));
            let sb = cloud_support(&moved, neg3(d));
            (sub3(sa, sb), sa, sb)
        };
        let result = toi.compute(support, velocity, dt);
        match first {
            Some(step) => {
                let (mut lo, mut hi) = (
                    dt * (step - 1) as f64 / steps as f64,
                    dt * step as f64 / steps as f64,
                );
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if gap_at(mid) > 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let found = result
                    .toi
                    .unwrap_or_else(|| panic!("case {case}: missed the contact at {hi}"));
                assert!(
                    (found - hi).abs() < 1e-6,
                    "case {case}: toi {found} vs reference {hi}"
                );
                checked += 1;
            }
            None if min_gap > 1e-3 => {
                assert!(
                    result.toi.is_none(),
                    "case {case}: spurious hit {:?}",
                    result.toi
                );
                checked += 1;
            }
            None => {}
        }
    }
    assert!(checked > 60, "only {checked} decisive sweeps");
}

#[test]
fn identity_rotation_is_the_axis_aligned_box() {
    // Guards the rotated-box helper used above.
    let corners = box_corners([1.0, 2.0, 3.0], [0.5, 1.0, 1.5], IDENTITY);
    let (lo, hi) = interval(&corners, [0.0, 1.0, 0.0]);
    assert!((lo - 1.0).abs() < 1e-12 && (hi - 3.0).abs() < 1e-12);
}
