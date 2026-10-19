// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the warm start, the support-cache variant and the EPA of
//! [`crate::gjk_enhanced`].
//!
//! * `gjk_enhanced` stopped on a boolean separating-axis test and returned the
//!   first support point's distance (3.041 instead of 3.0 for unit boxes at
//!   `(5, 1.5, 0)`); `gjk_with_support_cache` still carried that exit.
//! * A warm start seeded the simplex with the vertices of the previous frame, so
//!   a cached intersection (four vertices) grew to five on the first support
//!   point and was reported as separated.
//! * `epa_enhanced` removed only the face it expanded, leaving the other visible
//!   faces in place and the polytope non-convex, and used the 1/3 centroid
//!   weights of the face as contact points.

use super::functions::{
    ConvexShape, closest_point_on_segment, epa_enhanced, gjk_closest_points, gjk_distance_enhanced,
    gjk_enhanced, gjk_epa_pipeline, gjk_with_support_cache, johnson_tetrahedron, vadd, vdot, vlen,
    vscale, vsub,
};
use super::types::{
    ConvexHull, EnhBox, EnhCapsule, EnhEllipsoid, EnhSimplex, EnhSphere, SupportCache,
    SupportPoint, WarmStartCache,
};
use crate::gjk_epa::test_support::{
    Rng, axis_separation, box_axes, box_corners, interval, polyhedron_separation, rotation_about,
    translated,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
    assert!(
        vlen(vsub(actual, expected)) <= tol,
        "expected {expected:?}, got {actual:?} (tolerance {tol})"
    );
}

// ── gjk_enhanced distance (the separating-axis exit) ─────────────────────────

#[test]
fn offset_boxes_report_the_face_gap_everywhere() {
    // Unit boxes (half extent 1) at the origin and at (5, 1.5, 0): the faces x = 1
    // and x = 4 face each other, 3 apart. The boolean exit returned 3.041.
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([5.0, 1.5, 0.0], [1.0; 3]);
    let result = gjk_enhanced(&a, &b, None, 64);
    assert!(!result.intersecting);
    assert!(
        (result.distance() - 3.0).abs() < 1e-9,
        "{}",
        result.distance()
    );
    assert!((gjk_distance_enhanced(&a, &b) - 3.0).abs() < 1e-9);
    let (pa, pb) = gjk_closest_points(&a, &b).expect("separated boxes");
    assert!(
        (pa[0] - 1.0).abs() < 1e-9 && (pb[0] - 4.0).abs() < 1e-9,
        "{pa:?} {pb:?}"
    );
    assert_vec_close(vsub(pb, pa), [3.0, 0.0, 0.0], 1e-9);

    let contact = gjk_epa_pipeline(&a, &b, None, 64, 32, 1e-6);
    assert!(!contact.in_contact);
    assert!((contact.signed_distance - 3.0).abs() < 1e-9);
    // The normal points from B to A.
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-9);
}

#[test]
fn sphere_and_box_report_the_exact_gap() {
    let sphere = EnhSphere::new([0.0; 3], 1.0);
    let cuboid = EnhBox::new([4.0, 0.0, 0.0], [1.0; 3]);
    let result = gjk_enhanced(&sphere, &cuboid, None, 64);
    assert!(!result.intersecting);
    assert!(
        (result.distance() - 2.0).abs() < 1e-9,
        "{}",
        result.distance()
    );
    assert_vec_close(result.closest_a, [1.0, 0.0, 0.0], 1e-6);
    assert_vec_close(result.closest_b, [3.0, 0.0, 0.0], 1e-6);
}

#[test]
fn curved_shapes_report_the_distance_to_the_documented_accuracy() {
    // Spheres and capsules have no vertices, so the iteration converges only
    // geometrically. It stops once no support point can bring the simplex closer
    // than 1e-9 * (1 + distance), and the closed form has to agree to that.
    let mut rng = Rng(0xC0DE_0000_5F3E_0004);
    let mut spheres = 0;
    for case in 0..300 {
        let a = EnhSphere::new(rng.vec(-1.0, 1.0), rng.range(0.2, 1.2));
        let b = EnhSphere::new(rng.vec(-4.0, 4.0), rng.range(0.2, 1.2));
        let between = vsub(b.centre, a.centre);
        let exact = vlen(between) - a.radius - b.radius;
        if exact < 1e-3 {
            continue;
        }
        let result = gjk_enhanced(&a, &b, None, 64);
        assert!(!result.intersecting, "case {case}");
        assert!(
            (result.distance() - exact).abs() <= 1e-9 * (1.0 + exact),
            "case {case}: {} vs {exact}",
            result.distance()
        );
        // The witnesses realise the distance and sit near the surface points on
        // the line of centres (smooth shapes: accurate to about the square root
        // of the distance tolerance).
        assert!((vlen(vsub(result.closest_a, result.closest_b)) - result.distance()).abs() < 1e-12);
        let toward_b = vscale(between, 1.0 / vlen(between));
        let expected_a = vadd(a.centre, vscale(toward_b, a.radius));
        assert!(
            vlen(vsub(result.closest_a, expected_a)) < 5e-3,
            "case {case}: witness {:?} vs {expected_a:?}",
            result.closest_a
        );
        spheres += 1;
    }
    assert!(spheres > 150, "only {spheres} separated sphere pairs");

    let mut capsules = 0;
    for case in 0..300 {
        let capsule = EnhCapsule::new(rng.vec(-2.0, 2.0), rng.vec(-2.0, 2.0), rng.range(0.1, 0.6));
        let sphere = EnhSphere::new(rng.vec(-4.0, 4.0), rng.range(0.2, 0.8));
        let (nearest, _) = closest_point_on_segment(capsule.a, capsule.b, sphere.centre);
        let exact = vlen(vsub(sphere.centre, nearest)) - capsule.radius - sphere.radius;
        if exact < 1e-3 {
            continue;
        }
        let result = gjk_enhanced(&capsule, &sphere, None, 64);
        assert!(!result.intersecting, "case {case}");
        assert!(
            (result.distance() - exact).abs() <= 1e-9 * (1.0 + exact),
            "case {case}: {} vs {exact}",
            result.distance()
        );
        capsules += 1;
    }
    assert!(capsules > 100, "only {capsules} separated capsule pairs");
}

#[test]
fn ellipsoid_support_is_the_furthest_surface_point() {
    // The support point of an ellipsoid with semi-axes a_i is
    // c + A^2 d / sqrt(d' A^2 d). The old formula returned a point inside the
    // ellipsoid (x = 1 instead of x = 2 along +x for semi-axes (2, 1, 1)).
    let ellipsoid = EnhEllipsoid::new([0.5, -1.0, 2.0], [2.0, 1.0, 0.5]);
    assert_vec_close(ellipsoid.support([1.0, 0.0, 0.0]), [2.5, -1.0, 2.0], 1e-12);
    assert_vec_close(ellipsoid.support([0.0, -3.0, 0.0]), [0.5, -2.0, 2.0], 1e-12);
    assert_vec_close(ellipsoid.support([0.0, 0.0, 0.1]), [0.5, -1.0, 2.5], 1e-12);
    let mut rng = Rng(0x00E1_11D5_0001_A5E0);
    for case in 0..200 {
        let direction = rng.vec(-1.0, 1.0);
        let length = vlen(direction);
        if length < 1e-3 {
            continue;
        }
        let point = ellipsoid.support(direction);
        // On the surface: ((p - c) / a) has unit length.
        let offset = vsub(point, ellipsoid.centre);
        let level = (0..3)
            .map(|k| (offset[k] / ellipsoid.semi_axes[k]).powi(2))
            .sum::<f64>();
        assert!((level - 1.0).abs() < 1e-12, "case {case}: level {level}");
        // Its projection is the extent sqrt(d' A^2 d) of the ellipsoid along d.
        let extent = (0..3)
            .map(|k| (ellipsoid.semi_axes[k] * direction[k]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(
            (vdot(offset, direction) - extent).abs() < 1e-12,
            "case {case}: {} vs {extent}",
            vdot(offset, direction)
        );
        // No other surface point reaches further along the direction: scaling a
        // point of the unit sphere by the semi-axes gives a point of the ellipsoid.
        for _ in 0..50 {
            let sample = rng.vec(-1.0, 1.0);
            let length = vlen(sample);
            if length < 1e-3 {
                continue;
            }
            let on_surface = [
                sample[0] / length * ellipsoid.semi_axes[0],
                sample[1] / length * ellipsoid.semi_axes[1],
                sample[2] / length * ellipsoid.semi_axes[2],
            ];
            assert!(vdot(on_surface, direction) <= extent + 1e-12, "case {case}");
        }
    }
    // A direction without extent is supported by the centre.
    assert_eq!(ellipsoid.support([0.0; 3]), ellipsoid.centre);
    let flat = EnhEllipsoid::new([1.0, 1.0, 1.0], [2.0, 0.0, 1.0]);
    assert_eq!(flat.support([0.0, 1.0, 0.0]), flat.centre);
}

#[test]
fn ellipsoid_distances_follow_the_semi_axes() {
    // Semi-axes (2, 1, 1): the extent along x is 2 and along y it is 1.
    let ellipsoid = EnhEllipsoid::new([0.0; 3], [2.0, 1.0, 1.0]);
    let along_x = gjk_enhanced(&ellipsoid, &EnhSphere::new([5.0, 0.0, 0.0], 0.5), None, 64);
    assert!(!along_x.intersecting);
    assert!(
        (along_x.distance() - 2.5).abs() < 1e-9,
        "{}",
        along_x.distance()
    );
    let along_y = gjk_enhanced(&ellipsoid, &EnhSphere::new([0.0, 5.0, 0.0], 0.5), None, 64);
    assert!(
        (along_y.distance() - 3.5).abs() < 1e-9,
        "{}",
        along_y.distance()
    );
    // Overlap is decided by the real extent: a sphere reaching x = 1.8 touches
    // the ellipsoid (x extent 2) but would miss the old, smaller shape.
    let reaching = EnhSphere::new([2.3, 0.0, 0.0], 0.5);
    assert!(gjk_enhanced(&ellipsoid, &reaching, None, 64).intersecting);
}

#[test]
fn an_empty_iteration_budget_still_measures_one_support_pair() {
    // Unit boxes whose facing faces are 3 apart. With no budget the answer used
    // to be "not intersecting, distance 0" without a single support query.
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([5.0, 1.5, 0.0], [1.0; 3]);
    let none = gjk_enhanced(&a, &b, None, 0);
    assert!(!none.intersecting);
    assert_eq!(none.metrics.iterations, 1);
    assert_eq!(none.simplex.len(), 1);
    // One support pair: an upper bound of the distance that its witnesses realise.
    assert!(none.distance() >= 3.0 - 1e-9, "{}", none.distance());
    assert!((vlen(vsub(none.closest_a, none.closest_b)) - none.distance()).abs() < 1e-12);
    // A budget of one iteration is the same query, and a generous budget
    // converges on the true distance.
    let one = gjk_enhanced(&a, &b, None, 1);
    assert_eq!(one.closest_a, none.closest_a);
    assert_eq!(one.closest_b, none.closest_b);
    assert!((gjk_enhanced(&a, &b, None, 64).distance() - 3.0).abs() < 1e-9);
}

#[test]
fn support_cache_variant_returns_the_true_distance() {
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([5.0, 1.5, 0.0], [1.0; 3]);
    let mut cache = SupportCache::new();
    let cached = gjk_with_support_cache(&a, &b, &mut cache, 64);
    assert!(!cached.intersecting);
    assert!(
        (cached.distance() - 3.0).abs() < 1e-9,
        "{}",
        cached.distance()
    );
    assert!(cache.total > 0);
    assert!(cached.metrics.support_calls > 0);

    let sphere = EnhSphere::new([0.0; 3], 1.0);
    let cuboid = EnhBox::new([4.0, 0.0, 0.0], [1.0; 3]);
    let mut cache = SupportCache::new();
    let cached = gjk_with_support_cache(&sphere, &cuboid, &mut cache, 64);
    assert!((cached.distance() - 2.0).abs() < 1e-9);

    // Overlapping shapes are still detected.
    let c = EnhBox::new([1.0, 0.3, 0.2], [1.0; 3]);
    let mut cache = SupportCache::new();
    assert!(gjk_with_support_cache(&a, &c, &mut cache, 64).intersecting);
}

#[test]
fn support_cache_variant_agrees_with_the_plain_iteration() {
    let mut rng = Rng(0x5A11_CAC4_E000_0001);
    for case in 0..200 {
        let cloud_a: Vec<[f64; 3]> = (0..6).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..6).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            rng.vec(-3.0, 3.0),
        );
        let gap = polyhedron_separation(&cloud_a, &cloud_b).gap;
        if gap.abs() < 1e-6 {
            continue;
        }
        let a = ConvexHull::new(cloud_a);
        let b = ConvexHull::new(cloud_b);
        let plain = gjk_enhanced(&a, &b, None, 64);
        let mut cache = SupportCache::new();
        let cached = gjk_with_support_cache(&a, &b, &mut cache, 64);
        assert_eq!(plain.intersecting, gap < 0.0, "case {case}");
        assert_eq!(cached.intersecting, plain.intersecting, "case {case}");
        if !plain.intersecting {
            assert!(
                (cached.distance() - plain.distance()).abs() < 1e-9,
                "case {case}: {} vs {}",
                cached.distance(),
                plain.distance()
            );
        }
    }
}

/// Exact distance between two axis-aligned boxes (zero if they overlap).
fn box_gap(a: &EnhBox, b: &EnhBox) -> f64 {
    (0..3)
        .map(|k| {
            let gap = ((a.centre[k] - b.centre[k]).abs() - a.half[k] - b.half[k]).max(0.0);
            gap * gap
        })
        .sum::<f64>()
        .sqrt()
}

#[test]
fn tie_heavy_box_pairs_have_the_exact_distance_and_verdict() {
    // Coordinates on a coarse grid make support ties, aligned faces and exactly
    // touching boxes frequent: the situations in which simplices go flat.
    let mut rng = Rng(0x71E5_0000_B0C5_0002);
    let grid = |rng: &mut Rng, span: f64| ((rng.unit() * 2.0 - 1.0) * span * 4.0).round() / 4.0;
    let mut separated = 0;
    for case in 0..800 {
        let half = |rng: &mut Rng| [0.5 + 0.5 * (rng.unit() * 2.0).floor(); 3];
        let a = EnhBox::new(
            [
                grid(&mut rng, 2.0),
                grid(&mut rng, 2.0),
                grid(&mut rng, 2.0),
            ],
            half(&mut rng),
        );
        let b = EnhBox::new(
            [
                grid(&mut rng, 4.0),
                grid(&mut rng, 4.0),
                grid(&mut rng, 4.0),
            ],
            half(&mut rng),
        );
        let exact = box_gap(&a, &b);
        let result = gjk_enhanced(&a, &b, None, 64);
        if exact > 1e-9 {
            separated += 1;
            assert!(!result.intersecting, "case {case}: {a:?} vs {b:?}");
            assert!(
                (result.distance() - exact).abs() < 1e-8,
                "case {case}: {} vs exact {exact} for {a:?} vs {b:?}",
                result.distance()
            );
        } else {
            // Overlapping or exactly touching: either the verdict or a zero gap.
            assert!(
                result.intersecting || result.distance() < 1e-8,
                "case {case}: {} for {a:?} vs {b:?}",
                result.distance()
            );
        }
    }
    assert!(separated > 150, "only {separated} separated pairs");
}

#[test]
fn flat_tetrahedron_has_no_interior() {
    // Four corners of a square in the plane x = 2: the origin is 2 away from the
    // plane. The sign tests alone pass for every origin because each apex lies in
    // the plane of the opposite face.
    let (closest, inside) = johnson_tetrahedron(
        [2.0, -1.0, -1.0],
        [2.0, 1.0, -1.0],
        [2.0, 1.0, 1.0],
        [2.0, -1.0, 1.0],
    );
    assert!(!inside);
    assert_vec_close(closest, [2.0, 0.0, 0.0], 1e-12);
    // A flat tetrahedron whose hull contains the origin in its plane touches it.
    let (closest, inside) = johnson_tetrahedron(
        [1.0, 0.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [0.5, 0.5, 0.0],
    );
    assert!(!inside);
    assert!(vlen(closest) < 1e-12, "{closest:?}");
}

// ── warm start ───────────────────────────────────────────────────────────────

/// Generic polyhedra: GJK ends on a full tetrahedron around the origin.
fn overlapping_hulls(rng: &mut Rng) -> (ConvexHull, ConvexHull) {
    loop {
        let cloud_a: Vec<[f64; 3]> = (0..8).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..8).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            rng.vec(-0.4, 0.4),
        );
        if polyhedron_separation(&cloud_a, &cloud_b).gap < -0.1 {
            return (ConvexHull::new(cloud_a), ConvexHull::new(cloud_b));
        }
    }
}

fn shifted(hull: &ConvexHull, offset: [f64; 3]) -> ConvexHull {
    ConvexHull::new(translated(&hull.vertices, offset))
}

#[test]
fn cached_intersection_still_intersects_and_then_separates() {
    let mut rng = Rng(0xA0A0_0000_0000_0007);
    for case in 0..40 {
        let (a, b) = overlapping_hulls(&mut rng);
        let first = gjk_enhanced(&a, &b, None, 64);
        assert!(first.intersecting, "case {case}");
        assert_eq!(first.simplex.len(), 4, "case {case}");
        let mut cache = WarmStartCache::new();
        cache.update(&first.simplex, vnorm_or_x(rng.vec(-1.0, 1.0)));

        // Persistent contact: the warm-started query must report the overlap
        // again. The stale four-vertex simplex plus one new support point used to
        // reach `reduce` with five vertices, which returned "not intersecting".
        let again = gjk_enhanced(&a, &b, Some(&cache), 64);
        assert!(again.metrics.warm_started, "case {case}");
        assert!(again.intersecting, "case {case}: persistent contact lost");

        // After B moved far away the same cache must give a clean separation.
        let far = shifted(&b, [10.0, 0.0, 0.0]);
        let moved = gjk_enhanced(&a, &far, Some(&cache), 64);
        assert!(!moved.intersecting, "case {case}");
        let cold = gjk_enhanced(&a, &far, None, 64);
        assert!(
            (moved.distance() - cold.distance()).abs() < 1e-8,
            "case {case}: warm {} vs cold {}",
            moved.distance(),
            cold.distance()
        );
    }
}

fn vnorm_or_x(v: [f64; 3]) -> [f64; 3] {
    let l = vlen(v);
    if l > 1e-9 {
        vscale(v, 1.0 / l)
    } else {
        [1.0, 0.0, 0.0]
    }
}

#[test]
fn cached_intersecting_spheres_stay_intersecting() {
    let a = EnhSphere::new([0.0; 3], 2.0);
    let b = EnhSphere::new([1.0, 0.0, 0.0], 2.0);
    let first = gjk_enhanced(&a, &b, None, 64);
    assert!(first.intersecting);
    let mut cache = WarmStartCache::new();
    cache.update(&first.simplex, [1.0, 0.0, 0.0]);
    assert!(gjk_enhanced(&a, &b, Some(&cache), 64).intersecting);
    let far = EnhSphere::new([10.0, 0.0, 0.0], 2.0);
    let moved = gjk_enhanced(&a, &far, Some(&cache), 64);
    assert!(!moved.intersecting);
    assert!(
        (moved.distance() - 6.0).abs() < 1e-9,
        "{}",
        moved.distance()
    );
}

#[test]
fn warm_start_reduces_to_the_cold_answer_with_the_exact_distance() {
    let s1 = EnhSphere::new([0.0; 3], 1.0);
    let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
    let res1 = gjk_enhanced(&s1, &s2, None, 64);
    let mut cache = WarmStartCache::new();
    cache.update(&res1.simplex, [1.0, 0.0, 0.0]);
    let s2b = EnhSphere::new([5.1, 0.0, 0.0], 1.0);
    let warm = gjk_enhanced(&s1, &s2b, Some(&cache), 64);
    assert!(warm.metrics.warm_started);
    assert!(!warm.intersecting);
    assert!((warm.distance() - 3.1).abs() < 1e-9, "{}", warm.distance());
}

#[test]
fn unusable_cached_directions_fall_back_to_a_cold_start() {
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([5.0, 0.4, 0.0], [1.0; 3]);
    for direction in [[0.0; 3], [f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
        let mut cache = WarmStartCache::new();
        cache.valid = true;
        cache.direction = direction;
        let result = gjk_enhanced(&a, &b, Some(&cache), 64);
        assert!(!result.metrics.warm_started, "{direction:?}");
        assert!((result.distance() - 3.0).abs() < 1e-9, "{direction:?}");
    }
    // An invalid cache is ignored whatever it holds.
    let mut cache = WarmStartCache::new();
    cache.direction = [0.0, 1.0, 0.0];
    assert!(!gjk_enhanced(&a, &b, Some(&cache), 64).metrics.warm_started);
}

#[test]
fn warm_and_cold_queries_agree_on_random_configurations() {
    let mut rng = Rng(0xCA5E_0000_C01D_0003);
    let mut compared = 0;
    for case in 0..250 {
        let cloud_a: Vec<[f64; 3]> = (0..7).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = (0..7).map(|_| rng.vec(-1.0, 1.0)).collect();
        let a = ConvexHull::new(cloud_a);
        let before = shifted(&ConvexHull::new(cloud_b.clone()), rng.vec(-3.0, 3.0));
        let after_offset = rng.vec(-3.0, 3.0);
        let after = shifted(&ConvexHull::new(cloud_b), after_offset);
        let previous = gjk_enhanced(&a, &before, None, 64);
        let mut cache = WarmStartCache::new();
        // The cache of the previous frame: its simplex and a search direction.
        let toward = vsub(previous.closest_b, previous.closest_a);
        cache.update(&previous.simplex, vnorm_or_x(toward));
        let gap = polyhedron_separation(&a.vertices, &after.vertices).gap;
        if gap.abs() < 1e-6 {
            continue;
        }
        let warm = gjk_enhanced(&a, &after, Some(&cache), 64);
        let cold = gjk_enhanced(&a, &after, None, 64);
        assert_eq!(warm.intersecting, gap < 0.0, "case {case}: warm verdict");
        assert_eq!(cold.intersecting, gap < 0.0, "case {case}: cold verdict");
        if !warm.intersecting {
            assert!(
                (warm.distance() - cold.distance()).abs() < 1e-8,
                "case {case}: warm {} vs cold {}",
                warm.distance(),
                cold.distance()
            );
        }
        compared += 1;
    }
    assert!(compared > 150, "only {compared} decisive configurations");
}

fn support_point(w: [f64; 3]) -> SupportPoint {
    SupportPoint::new(w, w, [0.0; 3])
}

/// The support point of the Minkowski difference `a - b` along `d`.
fn support_pair(a: &dyn ConvexShape, b: &dyn ConvexShape, d: [f64; 3]) -> SupportPoint {
    let sa = a.support(d);
    let sb = b.support(vscale(d, -1.0));
    SupportPoint::new(vsub(sa, sb), sa, sb)
}

#[test]
fn simplex_with_more_than_four_vertices_reduces_the_newest_four() {
    // Two stale vertices far away followed by a tetrahedron around the origin.
    let mut simplex = EnhSimplex::new();
    simplex.add(support_point([9.0, 9.0, 9.0]));
    simplex.add(support_point([8.0, 9.0, 9.0]));
    simplex.add(support_point([1.0, 1.0, 1.0]));
    simplex.add(support_point([-1.0, -1.0, 1.0]));
    simplex.add(support_point([-1.0, 1.0, -1.0]));
    simplex.add(support_point([1.0, -1.0, -1.0]));
    let (_, contains) = simplex.reduce();
    assert!(contains);
    assert_eq!(simplex.verts.len(), 4);
    assert_eq!(simplex.verts[0].w, [1.0, 1.0, 1.0]);

    // Five vertices whose newest four are on one side: the direction points to
    // the origin instead of the old constant +x.
    let mut simplex = EnhSimplex::new();
    for w in [
        [-5.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [4.0, 1.0, 0.0],
        [3.5, 0.5, 1.0],
        [3.5, 0.5, -1.0],
    ] {
        simplex.add(support_point(w));
    }
    let (direction, contains) = simplex.reduce();
    assert!(!contains);
    assert!(direction[0] < 0.0, "{direction:?}");
}

#[test]
fn from_cache_exposes_the_stale_vertices_without_changing_them() {
    let mut cache = WarmStartCache::new();
    cache.update(
        &[
            support_point([1.0, 2.0, 3.0]),
            support_point([4.0, 5.0, 6.0]),
        ],
        [1.0, 0.0, 0.0],
    );
    let simplex = EnhSimplex::from_cache(&cache);
    assert_eq!(simplex.verts, cache.points);
}

// ── epa_enhanced ─────────────────────────────────────────────────────────────

/// The deepest points of both hulls lie on the extremal faces along the normal.
fn assert_deepest_points(
    point_a: [f64; 3],
    point_b: [f64; 3],
    normal_b_to_a: [f64; 3],
    cloud_a: &[[f64; 3]],
    cloud_b: &[[f64; 3]],
    tol: f64,
    label: &str,
) {
    let toward_b = vscale(normal_b_to_a, -1.0);
    let (_, a_hi) = interval(cloud_a, toward_b);
    let (b_lo, _) = interval(cloud_b, toward_b);
    assert!(
        (vdot(point_a, toward_b) - a_hi).abs() < tol,
        "{label}: point on A is not on its deepest face"
    );
    assert!(
        (vdot(point_b, toward_b) - b_lo).abs() < tol,
        "{label}: point on B is not on its deepest face"
    );
}

#[test]
fn rotated_box_overlap_of_a_quarter_matches_the_separating_axis_depth() {
    // Box B is rotated 30 degrees about z; move it along x until the
    // separating-axis depth is exactly 0.25.
    let rotation = rotation_about([0.0, 0.0, 1.0], 30.0_f64.to_radians());
    let half = [1.0; 3];
    let axes = box_axes(IDENTITY, rotation);
    let cloud_a = box_corners([0.0; 3], half, IDENTITY);
    let depth_at =
        |x: f64| -axis_separation(&cloud_a, &box_corners([x, 0.0, 0.0], half, rotation), &axes).gap;
    let (mut near, mut far) = (1.2, 3.0);
    for _ in 0..80 {
        let mid = 0.5 * (near + far);
        if depth_at(mid) > 0.25 {
            near = mid;
        } else {
            far = mid;
        }
    }
    let center_b = [0.5 * (near + far), 0.0, 0.0];
    assert!((depth_at(center_b[0]) - 0.25).abs() < 1e-12);
    let cloud_b = box_corners(center_b, half, rotation);

    let a = ConvexHull::new(cloud_a.clone());
    let b = ConvexHull::new(cloud_b.clone());
    let gjk = gjk_enhanced(&a, &b, None, 64);
    assert!(gjk.intersecting);
    let result = epa_enhanced(&a, &b, gjk.simplex.clone(), 64, 1e-10);
    assert!((result.depth - 0.25).abs() < 1e-6, "depth {}", result.depth);
    // B-to-A normal: A is on the -x side of B.
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-6);
    // The deepest point of A is on its +x face, that of B on its rotated edge.
    let point_a = vsub(
        result.contact_point,
        vscale(result.normal, 0.5 * result.depth),
    );
    let point_b = vadd(
        result.contact_point,
        vscale(result.normal, 0.5 * result.depth),
    );
    assert!((point_a[0] - 1.0).abs() < 1e-6, "{point_a:?}");
    assert!((point_b[0] - 0.75).abs() < 1e-6, "{point_b:?}");
    assert_deepest_points(
        point_a,
        point_b,
        result.normal,
        &cloud_a,
        &cloud_b,
        1e-6,
        "30 degrees",
    );

    // The combined pipeline reports the same contact.
    let contact = gjk_epa_pipeline(&a, &b, None, 64, 64, 1e-10);
    assert!(contact.in_contact);
    assert!((contact.signed_distance + 0.25).abs() < 1e-6);
    assert_vec_close(contact.normal, [-1.0, 0.0, 0.0], 1e-6);
    assert_vec_close(
        vsub(contact.point_a, contact.point_b),
        vscale(contact.normal, -0.25),
        1e-6,
    );
    assert_deepest_points(
        contact.point_a,
        contact.point_b,
        contact.normal,
        &cloud_a,
        &cloud_b,
        1e-6,
        "pipeline",
    );
}

#[test]
fn epa_enhanced_matches_the_separating_axis_solution_on_random_polyhedra() {
    let mut rng = Rng(0xE9A0_0000_BEEF_0005);
    let mut verified = 0;
    for case in 0..250 {
        let cloud_a: Vec<[f64; 3]> = (0..(5 + (rng.unit() * 4.0) as usize))
            .map(|_| rng.vec(-1.0, 1.0))
            .collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..(5 + (rng.unit() * 4.0) as usize))
                .map(|_| rng.vec(-1.0, 1.0))
                .collect::<Vec<_>>(),
            rng.vec(-0.7, 0.7),
        );
        let separation = polyhedron_separation(&cloud_a, &cloud_b);
        if separation.gap > -0.02 {
            continue;
        }
        let a = ConvexHull::new(cloud_a.clone());
        let b = ConvexHull::new(cloud_b.clone());
        let gjk = gjk_enhanced(&a, &b, None, 64);
        assert!(gjk.intersecting, "case {case}");
        let result = epa_enhanced(&a, &b, gjk.simplex.clone(), 96, 1e-10);
        assert!(
            (result.depth + separation.gap).abs() < 1e-6,
            "case {case}: depth {} vs exact {}",
            result.depth,
            -separation.gap
        );
        assert!(
            vdot(result.normal, separation.from_b_to_a) > 1.0 - 1e-6,
            "case {case}: normal {:?} vs exact {:?}",
            result.normal,
            separation.from_b_to_a
        );
        let point_a = vsub(
            result.contact_point,
            vscale(result.normal, 0.5 * result.depth),
        );
        let point_b = vadd(
            result.contact_point,
            vscale(result.normal, 0.5 * result.depth),
        );
        assert_deepest_points(
            point_a,
            point_b,
            result.normal,
            &cloud_a,
            &cloud_b,
            1e-6,
            &format!("case {case}"),
        );
        verified += 1;
    }
    assert!(verified > 100, "only {verified} decisive overlaps");
}

#[test]
fn epa_enhanced_completes_any_gjk_simplex() {
    // Whatever GJK ended with (nothing, one point, a repeated point), the EPA
    // completes it to a tetrahedron and finds the same depth.
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([1.75, 0.2, -0.1], [1.0; 3]);
    let expected = 0.25;
    let corner = support_pair(&a, &b, [1.0, 0.0, 0.0]);
    for simplex in [vec![], vec![corner], vec![corner, corner, corner]] {
        let result = epa_enhanced(&a, &b, simplex, 64, 1e-10);
        assert!(
            (result.depth - expected).abs() < 1e-9,
            "depth {}",
            result.depth
        );
        assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-9);
    }
}

#[test]
fn epa_enhanced_grows_a_simplex_that_misses_the_origin() {
    // A tetrahedron of support points that does not contain the origin, taken
    // from separated-looking directions of an overlapping pair.
    let a = EnhSphere::new([0.0; 3], 1.0);
    let b = EnhSphere::new([1.5, 0.0, 0.0], 1.0);
    let seed = vec![
        support_pair(&a, &b, [1.0, 1.0, 1.0]),
        support_pair(&a, &b, [-1.0, -1.0, 1.0]),
        support_pair(&a, &b, [-1.0, 1.0, -1.0]),
        support_pair(&a, &b, [1.0, -1.0, -1.0]),
    ];
    let result = epa_enhanced(&a, &b, seed, 128, 1e-9);
    assert!((result.depth - 0.5).abs() < 1e-6, "depth {}", result.depth);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-3);
}

#[test]
fn epa_enhanced_without_penetration_reports_a_zero_depth() {
    // Separated shapes: no penetration exists. The result is the documented
    // degenerate one, not a fabricated depth.
    let a = EnhSphere::new([0.0; 3], 1.0);
    let b = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
    let result = epa_enhanced(&a, &b, vec![], 32, 1e-9);
    assert_eq!(result.depth, 0.0);
    assert_eq!(result.iterations, 0);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-12);
    assert_vec_close(result.contact_point, [2.5, 0.0, 0.0], 1e-12);
}

#[test]
fn pipeline_normal_points_from_b_to_a_in_both_cases() {
    let a = EnhSphere::new([0.0; 3], 1.0);
    // Separated: B on +x, so the normal from B to A is -x.
    let apart = gjk_epa_pipeline(
        &a,
        &EnhSphere::new([4.0, 0.0, 0.0], 1.0),
        None,
        64,
        64,
        1e-9,
    );
    assert!(!apart.in_contact);
    assert_vec_close(apart.normal, [-1.0, 0.0, 0.0], 1e-9);
    // Penetrating: the same side, the same normal, and A's contact point is the
    // deeper one (further along +x) with a-b equal to -normal * depth.
    let deep = gjk_epa_pipeline(
        &a,
        &EnhSphere::new([1.5, 0.0, 0.0], 1.0),
        None,
        64,
        128,
        1e-9,
    );
    assert!(deep.in_contact);
    assert!(
        (deep.signed_distance + 0.5).abs() < 1e-6,
        "{}",
        deep.signed_distance
    );
    assert_vec_close(deep.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert_vec_close(
        vsub(deep.point_a, deep.point_b),
        vscale(deep.normal, -0.5),
        1e-6,
    );
    assert!(deep.point_a[0] > deep.point_b[0]);
}
