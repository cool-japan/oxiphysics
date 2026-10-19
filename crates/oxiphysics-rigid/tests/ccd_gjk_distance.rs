// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `ccd::gjk_distance` and the conservative advancement
//! built on it.
//!
//! The old `gjk_distance` had no simplex: it only re-queried the support
//! points along the previous witness vector, so overlapping shapes (e.g. the
//! AABBs `[0,2]^3` and `[1,3]^3`, or two coincident spheres) were reported as
//! separated with a positive distance, and `gjk_conservative_advancement`
//! stepped over their time of impact.

use oxiphysics_rigid::ccd::{
    AabbSupport, CcdSweep, ConvexSupport, SphereSupport, gjk_conservative_advancement, gjk_distance,
};

fn aabb(min: [f64; 3], max: [f64; 3]) -> AabbSupport {
    AabbSupport { min, max }
}

fn sphere(center: [f64; 3], radius: f64) -> SphereSupport {
    SphereSupport { center, radius }
}

/// Deterministic generator for the randomised comparisons.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn random_aabb(&mut self) -> AabbSupport {
        let mut min = [0.0; 3];
        let mut max = [0.0; 3];
        for axis in 0..3 {
            let centre = self.range(-2.0, 2.0);
            let half = self.range(0.1, 1.5);
            min[axis] = centre - half;
            max[axis] = centre + half;
        }
        aabb(min, max)
    }
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Exact distance between two AABBs (zero if they intersect).
fn aabb_aabb_distance(a: &AabbSupport, b: &AabbSupport) -> f64 {
    let mut sq = 0.0;
    for axis in 0..3 {
        let gap = (a.min[axis] - b.max[axis]).max(b.min[axis] - a.max[axis]);
        if gap > 0.0 {
            sq += gap * gap;
        }
    }
    sq.sqrt()
}

/// Exact distance between a sphere and an AABB (zero if they intersect).
fn sphere_aabb_distance(s: &SphereSupport, b: &AabbSupport) -> f64 {
    let mut sq = 0.0;
    for axis in 0..3 {
        let clamped = s.center[axis].clamp(b.min[axis], b.max[axis]);
        let d = s.center[axis] - clamped;
        sq += d * d;
    }
    (sq.sqrt() - s.radius).max(0.0)
}

// ── gjk_distance ────────────────────────────────────────────────────────────

#[test]
fn overlapping_aabbs_are_reported_as_overlapping() {
    // The old implementation returned sq_dist = 11, overlapping = false.
    let a = aabb([0.0; 3], [2.0; 3]);
    let b = aabb([1.0; 3], [3.0; 3]);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!(result.overlapping);
    assert!(result.sq_dist < 1e-18, "sq_dist = {}", result.sq_dist);
}

#[test]
fn aabb_and_a_point_three_units_away_are_two_apart() {
    let a = aabb([-1.0; 3], [1.0; 3]);
    let point = sphere([3.0, 0.5, 0.0], 0.0);
    let result = gjk_distance(&a, &point, 1e-9, 64);
    assert!(!result.overlapping);
    assert!(
        (result.sq_dist - 4.0).abs() < 1e-9,
        "sq_dist = {}",
        result.sq_dist
    );
    for (got, want) in result.closest_a.iter().zip([1.0, 0.5, 0.0]) {
        assert!(
            (got - want).abs() < 1e-9,
            "closest_a = {:?}",
            result.closest_a
        );
    }
    for (got, want) in result.closest_b.iter().zip([3.0, 0.5, 0.0]) {
        assert!(
            (got - want).abs() < 1e-9,
            "closest_b = {:?}",
            result.closest_b
        );
    }
}

#[test]
fn coincident_spheres_overlap() {
    // The old implementation reported sq_dist = 4 r^2, overlapping = false.
    let a = sphere([0.0; 3], 1.0);
    let b = sphere([0.0; 3], 1.0);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!(result.overlapping);
    assert!(result.sq_dist < 1e-18, "sq_dist = {}", result.sq_dist);
}

#[test]
fn sphere_inside_a_larger_sphere_overlaps() {
    let outer = sphere([0.0; 3], 5.0);
    let inner = sphere([1.0, 0.5, -0.5], 0.5);
    assert!(gjk_distance(&outer, &inner, 1e-9, 64).overlapping);
    assert!(gjk_distance(&inner, &outer, 1e-9, 64).overlapping);
}

#[test]
fn face_to_face_boxes_one_apart() {
    let a = aabb([0.0; 3], [1.0; 3]);
    let b = aabb([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!(!result.overlapping);
    assert!(
        (result.sq_dist - 1.0).abs() < 1e-9,
        "sq_dist = {}",
        result.sq_dist
    );
    // The witnesses lie on the facing planes, one unit apart along x.
    assert!(
        (result.closest_a[0] - 1.0).abs() < 1e-9,
        "{:?}",
        result.closest_a
    );
    assert!(
        (result.closest_b[0] - 2.0).abs() < 1e-9,
        "{:?}",
        result.closest_b
    );
    assert!((norm(sub(result.closest_a, result.closest_b)) - 1.0).abs() < 1e-9);
}

#[test]
fn touching_boxes_count_as_overlapping() {
    let a = aabb([0.0; 3], [1.0; 3]);
    let b = aabb([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!(result.overlapping);
    assert!(result.sq_dist < 1e-18);
    let identical = gjk_distance(&a, &a.clone(), 1e-9, 64);
    assert!(identical.overlapping);
}

#[test]
fn disjoint_spheres_report_the_surface_gap_and_witness_points() {
    // Centres 5 apart along (0.8, 0.6, 0); radii 1 each: gap 3.
    let a = sphere([0.0; 3], 1.0);
    let b = sphere([4.0, 3.0, 0.0], 1.0);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!(!result.overlapping);
    assert!(
        (result.sq_dist.sqrt() - 3.0).abs() < 1e-6,
        "distance = {}",
        result.sq_dist.sqrt()
    );
    // The distance converges like the tolerance, the witness points of curved
    // shapes like its square root.
    for (got, want) in result.closest_a.iter().zip([0.8, 0.6, 0.0]) {
        assert!(
            (got - want).abs() < 1e-3,
            "closest_a = {:?}",
            result.closest_a
        );
    }
    for (got, want) in result.closest_b.iter().zip([3.2, 2.4, 0.0]) {
        assert!(
            (got - want).abs() < 1e-3,
            "closest_b = {:?}",
            result.closest_b
        );
    }
}

#[test]
fn tolerance_is_the_touching_distance() {
    let a = sphere([0.0; 3], 1.0);
    let b = sphere([2.0 + 5e-5, 0.0, 0.0], 1.0);
    // The 5e-5 gap is inside a 1e-4 tolerance ...
    let loose = gjk_distance(&a, &b, 1e-4, 64);
    assert!(loose.overlapping);
    // ... and outside a 1e-6 one, where the true gap is reported.
    let tight = gjk_distance(&a, &b, 1e-6, 64);
    assert!(!tight.overlapping);
    assert!(
        (tight.sq_dist.sqrt() - 5e-5).abs() < 1e-6,
        "distance = {}",
        tight.sq_dist.sqrt()
    );
}

#[test]
fn degenerate_arguments_do_not_panic() {
    let a = aabb([0.0; 3], [1.0; 3]);
    let b = aabb([4.0; 3], [5.0; 3]);
    // No iterations: an upper bound of the distance, finite.
    let zero_iterations = gjk_distance(&a, &b, 1e-9, 0);
    assert!(zero_iterations.sq_dist.is_finite());
    // Zero, negative and NaN tolerances run to convergence (or the cap).
    for tolerance in [0.0, -1.0, f64::NAN] {
        let result = gjk_distance(&a, &b, tolerance, 64);
        assert!(!result.overlapping, "tolerance {tolerance}");
        assert!(
            (result.sq_dist - 3.0 * 9.0).abs() < 1e-9,
            "tolerance {tolerance}: sq_dist = {}",
            result.sq_dist
        );
    }
}

#[test]
fn distance_matches_the_exact_aabb_aabb_distance() {
    let mut rng = Lcg(0x5eed_0001);
    let (mut overlapping, mut separated) = (0, 0);
    for case in 0..3000 {
        let a = rng.random_aabb();
        let b = rng.random_aabb();
        let exact = aabb_aabb_distance(&a, &b);
        let result = gjk_distance(&a, &b, 1e-9, 64);
        if exact == 0.0 {
            overlapping += 1;
            assert!(result.overlapping, "case {case}: {a:?} vs {b:?}");
            assert!(result.sq_dist < 1e-12, "case {case}: {}", result.sq_dist);
        } else {
            separated += 1;
            assert!(!result.overlapping, "case {case}: {a:?} vs {b:?}");
            assert!(
                (result.sq_dist.sqrt() - exact).abs() < 1e-7,
                "case {case}: gjk {} vs exact {exact}",
                result.sq_dist.sqrt()
            );
            // The witnesses are on the boxes and realise the distance.
            let gap = norm(sub(result.closest_a, result.closest_b));
            assert!((gap - exact).abs() < 1e-6, "case {case}: witness gap {gap}");
            for axis in 0..3 {
                assert!(
                    (a.min[axis] - 1e-9..=a.max[axis] + 1e-9).contains(&result.closest_a[axis]),
                    "case {case}: witness on A outside the box"
                );
                assert!(
                    (b.min[axis] - 1e-9..=b.max[axis] + 1e-9).contains(&result.closest_b[axis]),
                    "case {case}: witness on B outside the box"
                );
            }
        }
        // Symmetry.
        let swapped = gjk_distance(&b, &a, 1e-9, 64);
        assert_eq!(swapped.overlapping, result.overlapping, "case {case}");
        assert!(
            (swapped.sq_dist.sqrt() - result.sq_dist.sqrt()).abs() < 1e-7,
            "case {case}"
        );
    }
    assert!(
        overlapping > 300 && separated > 300,
        "{overlapping}/{separated}"
    );
}

#[test]
fn distance_matches_the_exact_sphere_aabb_distance() {
    let mut rng = Lcg(0x5eed_0002);
    let (mut overlapping, mut separated) = (0, 0);
    for case in 0..3000 {
        let centre = [
            rng.range(-3.0, 3.0),
            rng.range(-3.0, 3.0),
            rng.range(-3.0, 3.0),
        ];
        let s = sphere(centre, rng.range(0.05, 1.5));
        let b = rng.random_aabb();
        let exact = sphere_aabb_distance(&s, &b);
        let result = gjk_distance(&s, &b, 1e-9, 64);
        if exact == 0.0 {
            overlapping += 1;
            assert!(result.overlapping, "case {case}: {s:?} vs {b:?}");
        } else {
            separated += 1;
            assert!(!result.overlapping, "case {case}: {s:?} vs {b:?}");
            assert!(
                (result.sq_dist.sqrt() - exact).abs() < 1e-5,
                "case {case}: gjk {} vs exact {exact}",
                result.sq_dist.sqrt()
            );
        }
    }
    assert!(
        overlapping > 300 && separated > 300,
        "{overlapping}/{separated}"
    );
}

// ── gjk_conservative_advancement ────────────────────────────────────────────

#[test]
fn conservative_advancement_between_boxes_finds_the_impact_time() {
    // A spans x in [-11, -9], B spans [9, 11]: 18 m apart, closing at 5 m/s.
    let shape = aabb([-1.0; 3], [1.0; 3]);
    let toi = gjk_conservative_advancement(
        &shape,
        [-10.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        &shape,
        [10.0, 0.0, 0.0],
        [0.0; 3],
        10.0,
        64,
        1e-6,
    )
    .expect("the boxes collide");
    assert!(toi <= 3.6 + 1e-9, "must not overstep: toi = {toi}");
    assert!(3.6 - toi < 1e-5, "toi = {toi}");
}

#[test]
fn conservative_advancement_does_not_overstep_an_oblique_impact() {
    // Sphere A (r 0.5) flies past-and-into sphere B (r 0.5), offset 0.5 off axis.
    let a = sphere([0.0; 3], 0.5);
    let b = sphere([0.0; 3], 0.5);
    let pos_a = [-8.0, 0.0, 0.0];
    let vel_a = [2.0, 0.0, 0.0];
    let pos_b = [0.0, 0.5, 0.0];
    let vel_b = [0.0; 3];
    let exact = CcdSweep::compute_collision_time_sphere(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 10.0)
        .expect("the spheres collide");
    let toi = gjk_conservative_advancement(&a, pos_a, vel_a, &b, pos_b, vel_b, 10.0, 64, 1e-6)
        .expect("the spheres collide");
    // The GJK distance of curved shapes is an upper bound within the tolerance
    // (1e-6 m of travel at 2 m/s relative speed).
    assert!(toi <= exact + 1e-6, "overstepped: {toi} vs exact {exact}");
    assert!(exact - toi < 1e-5, "toi {toi} vs exact {exact}");
}

#[test]
fn conservative_advancement_lets_a_wide_pass_go_by() {
    // Centres pass 2 m apart, radii sum to 1 m: no contact.
    let a = sphere([0.0; 3], 0.5);
    let b = sphere([0.0; 3], 0.5);
    let toi = gjk_conservative_advancement(
        &a,
        [-8.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        &b,
        [0.0, 2.0, 0.0],
        [0.0; 3],
        10.0,
        64,
        1e-6,
    );
    assert_eq!(toi, None);
}

#[test]
fn conservative_advancement_reports_shapes_that_already_overlap() {
    let a = aabb([0.0; 3], [2.0; 3]);
    let b = aabb([1.0; 3], [3.0; 3]);
    let toi = gjk_conservative_advancement(
        &a,
        [0.0; 3],
        [1.0, 0.0, 0.0],
        &b,
        [0.0; 3],
        [0.0; 3],
        1.0,
        64,
        1e-6,
    );
    assert_eq!(toi, Some(0.0));
}

#[test]
fn conservative_advancement_of_resting_separated_shapes_is_none() {
    let a = aabb([-1.0; 3], [1.0; 3]);
    let toi = gjk_conservative_advancement(
        &a,
        [-5.0, 0.0, 0.0],
        [0.0; 3],
        &a,
        [5.0, 0.0, 0.0],
        [0.0; 3],
        10.0,
        64,
        1e-6,
    );
    assert_eq!(toi, None);
}

#[test]
fn support_functions_are_consistent_with_the_witness_points() {
    // The witness point on each shape is a support point in the direction of
    // the separating axis.
    let a = aabb([-1.0; 3], [1.0; 3]);
    let b = sphere([5.0, 0.0, 0.0], 1.0);
    let result = gjk_distance(&a, &b, 1e-9, 64);
    assert!((result.sq_dist.sqrt() - 3.0).abs() < 1e-6);
    assert!((result.closest_a[0] - a.support([1.0, 0.0, 0.0])[0]).abs() < 1e-6);
    assert!((result.closest_b[0] - 4.0).abs() < 1e-6);
}
