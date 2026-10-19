// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the EPA contact points of [`gjk_epa_contact`] and [`epa`].
//!
//! The contact points used to be `sum(l_i * shape.support(w_i))`: a Minkowski
//! difference vertex `w_i = a_i - b_i` was fed to each shape's support function
//! as if it were a direction, which recovers neither `a_i` nor `b_i`. The
//! polytope expansion also removed only the face it expanded and compared the
//! shared horizon edges by direction, which left interior edges behind whenever
//! the faces were not wound consistently.

use super::test_support::{
    Rng, axis_separation, box_axes, box_corners, interval, polyhedron_separation, rotation_about,
    translated,
};
use super::*;

const TOL: f64 = 1e-9;

fn cuboid(center: [f64; 3], half: [f64; 3]) -> GjkBox {
    GjkBox {
        center,
        half_extents: half,
    }
}

fn ball(center: [f64; 3], radius: f64) -> GjkSphere {
    GjkSphere { center, radius }
}

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
    assert!(
        len3(sub3(actual, expected)) <= tol,
        "expected {expected:?}, got {actual:?} (tolerance {tol})"
    );
}

/// `contact_a - contact_b` is exactly the penetration vector.
fn assert_penetration_identity(result: &EpaResult, tol: f64) {
    assert!((len3(result.normal) - 1.0).abs() < 1e-12);
    assert_vec_close(
        sub3(result.contact_a, result.contact_b),
        scale3(result.normal, -result.depth),
        tol,
    );
}

// ── axis-aligned boxes ───────────────────────────────────────────────────────

#[test]
fn overlapping_boxes_report_depth_and_contacts_on_the_box_surfaces() {
    // Half-extent-1 boxes whose centres are 1.75 apart on x overlap by 0.25.
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([1.75, 0.0, 0.0], [1.0; 3]);
    let result = gjk_epa_contact(&a, &b).expect("the boxes overlap");
    assert!((result.depth - 0.25).abs() < TOL, "depth {}", result.depth);
    // A is on the -x side of B: the normal from B to A is -x.
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], TOL);
    // The deepest point of A inside B is on A's +x face, B's on its -x face.
    assert!(
        (result.contact_a[0] - 1.0).abs() < TOL,
        "{:?}",
        result.contact_a
    );
    assert!(
        (result.contact_b[0] - 0.75).abs() < TOL,
        "{:?}",
        result.contact_b
    );
    for k in 0..3 {
        assert!(
            result.contact_a[k].abs() <= 1.0 + TOL,
            "{:?}",
            result.contact_a
        );
        assert!(
            (result.contact_b[k] - b.center[k]).abs() <= 1.0 + TOL,
            "{:?}",
            result.contact_b
        );
    }
    assert_penetration_identity(&result, TOL);

    // Swapping the shapes mirrors the result.
    let mirrored = gjk_epa_contact(&b, &a).expect("the boxes overlap");
    assert!((mirrored.depth - 0.25).abs() < TOL);
    assert_vec_close(mirrored.normal, [1.0, 0.0, 0.0], TOL);
    assert!((mirrored.contact_a[0] - 0.75).abs() < TOL);
    assert!((mirrored.contact_b[0] - 1.0).abs() < TOL);
    assert_penetration_identity(&mirrored, TOL);
}

#[test]
fn boxes_overlapping_along_every_axis_and_sign() {
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut offset = [0.05, -0.03, 0.02];
            offset[axis] = sign * 1.75;
            let a = cuboid([0.0; 3], [1.0; 3]);
            let b = cuboid(offset, [1.0; 3]);
            let result = gjk_epa_contact(&a, &b).expect("the boxes overlap");
            assert!(
                (result.depth - 0.25).abs() < TOL,
                "axis {axis} sign {sign}: depth {}",
                result.depth
            );
            let mut expected = [0.0; 3];
            expected[axis] = -sign;
            assert_vec_close(result.normal, expected, TOL);
            // A's contact lies on the face of A that faces B.
            assert!(
                (result.contact_a[axis] - sign).abs() < TOL,
                "axis {axis} sign {sign}: {:?}",
                result.contact_a
            );
            assert!(
                (result.contact_b[axis] - (offset[axis] - sign)).abs() < TOL,
                "axis {axis} sign {sign}: {:?}",
                result.contact_b
            );
            assert_penetration_identity(&result, TOL);
        }
    }
}

#[test]
fn coincident_boxes_penetrate_by_their_full_width() {
    // The Minkowski difference is a box of half extent 2 around the origin: every
    // facet is 2 away, and the origin is the centre of a symmetric polytope.
    let a = cuboid([0.0; 3], [1.0; 3]);
    let result = gjk_epa_contact(&a, &a).expect("identical boxes overlap");
    assert!((result.depth - 2.0).abs() < TOL, "depth {}", result.depth);
    assert_penetration_identity(&result, TOL);
}

#[test]
fn separated_and_distant_shapes_have_no_contact() {
    assert!(
        gjk_epa_contact(
            &cuboid([0.0; 3], [1.0; 3]),
            &cuboid([3.0, 0.2, 0.1], [1.0; 3])
        )
        .is_none()
    );
    assert!(gjk_epa_contact(&ball([0.0; 3], 1.0), &ball([2.5, 0.0, 0.0], 1.0)).is_none());
}

// ── round shapes ─────────────────────────────────────────────────────────────

#[test]
fn overlapping_spheres_report_depth_normal_and_surface_contacts() {
    let a = ball([0.0; 3], 1.0);
    let b = ball([1.5, 0.0, 0.0], 1.0);
    let result = gjk_epa_contact(&a, &b).expect("the spheres overlap");
    assert!((result.depth - 0.5).abs() < 1e-6, "depth {}", result.depth);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-3);
    // The contact points are blends of surface points of the sphere: just inside
    // the surface, near the line of centres.
    assert_vec_close(result.contact_a, [1.0, 0.0, 0.0], 2e-3);
    assert_vec_close(result.contact_b, [0.5, 0.0, 0.0], 2e-3);
    assert!(len3(result.contact_a) <= 1.0 + TOL);
    assert!(len3(sub3(result.contact_b, b.center)) <= 1.0 + TOL);
    assert_penetration_identity(&result, 1e-9);
}

#[test]
fn spheres_overlapping_in_a_generic_direction() {
    let a = ball([0.3, -0.2, 0.1], 1.2);
    let b = ball([1.5, 0.6, 0.8], 0.9);
    let delta = sub3(b.center, a.center);
    let distance = len3(delta);
    let result = gjk_epa_contact(&a, &b).expect("the spheres overlap");
    assert!(
        (result.depth - (a.radius + b.radius - distance)).abs() < 1e-6,
        "depth {}",
        result.depth
    );
    assert_vec_close(result.normal, scale3(delta, -1.0 / distance), 1e-3);
    assert_penetration_identity(&result, 1e-9);
    // The deepest point of A is on its surface towards B.
    let expected_a = add3(a.center, scale3(delta, a.radius / distance));
    assert_vec_close(result.contact_a, expected_a, 3e-3);
}

#[test]
fn capsule_and_sphere_overlap_matches_the_closed_form() {
    let capsule = GjkCapsule {
        p1: [0.0, -1.0, 0.0],
        p2: [0.0, 1.0, 0.0],
        radius: 0.5,
    };
    let sphere_b = ball([0.8, 0.3, 0.0], 0.5);
    let result = gjk_epa_contact(&capsule, &sphere_b).expect("the shapes overlap");
    // Distance from the sphere centre to the capsule axis is 0.8.
    assert!((result.depth - 0.2).abs() < 1e-6, "depth {}", result.depth);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert_vec_close(result.contact_a, [0.5, 0.3, 0.0], 3e-3);
    assert_vec_close(result.contact_b, [0.3, 0.3, 0.0], 3e-3);
    assert_penetration_identity(&result, 1e-9);
}

// ── polyhedra against the separating-axis solution ───────────────────────────

/// The contacts of a polyhedral pair lie on the extremal faces of the two
/// hulls along the contact normal.
fn assert_contacts_on_extremal_faces(
    result: &EpaResult,
    cloud_a: &[[f64; 3]],
    cloud_b: &[[f64; 3]],
    tol: f64,
    label: &str,
) {
    let toward_b = scale3(result.normal, -1.0);
    let (_, a_hi) = interval(cloud_a, toward_b);
    let (b_lo, _) = interval(cloud_b, toward_b);
    assert!(
        (dot3(result.contact_a, toward_b) - a_hi).abs() < tol,
        "{label}: contact on A is not on its deepest face"
    );
    assert!(
        (dot3(result.contact_b, toward_b) - b_lo).abs() < tol,
        "{label}: contact on B is not on its deepest face"
    );
}

#[test]
fn rotated_box_overlap_of_a_quarter_matches_the_separating_axis_depth() {
    // Box B is rotated 30 degrees about z; move it along x until the
    // separating-axis depth is exactly 0.25.
    let rotation = rotation_about([0.0, 0.0, 1.0], 30.0_f64.to_radians());
    let half = [1.0, 1.0, 1.0];
    let axes = box_axes(IDENTITY_ROTATION, rotation);
    let cloud_a = box_corners([0.0; 3], half, IDENTITY_ROTATION);
    let depth_at = |x: f64| {
        let cloud_b = box_corners([x, 0.0, 0.0], half, rotation);
        -axis_separation(&cloud_a, &cloud_b, &axes).gap
    };
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

    let corners_b = box_corners(center_b, half, rotation);
    let hull_a = GjkConvexHull::new(cloud_a.clone());
    let hull_b = GjkConvexHull::new(corners_b.clone());
    let result = gjk_epa_contact(&hull_a, &hull_b).expect("the boxes overlap");
    assert!((result.depth - 0.25).abs() < 1e-6, "depth {}", result.depth);
    // The rotated corner of B pokes into A through A's +x face.
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-6);
    assert_penetration_identity(&result, 1e-8);
    assert_contacts_on_extremal_faces(&result, &cloud_a, &corners_b, 1e-6, "30 degrees");
}

const IDENTITY_ROTATION: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

#[test]
fn random_polyhedra_match_the_separating_axis_solution() {
    let mut rng = Rng(0xE9A0_BEEF_0000_0003);
    let mut verified = 0;
    for case in 0..250 {
        let count_a = 5 + (rng.unit() * 4.0) as usize;
        let count_b = 5 + (rng.unit() * 4.0) as usize;
        let cloud_a: Vec<[f64; 3]> = (0..count_a).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = translated(
            &(0..count_b).map(|_| rng.vec(-1.0, 1.0)).collect::<Vec<_>>(),
            rng.vec(-0.7, 0.7),
        );
        let separation = polyhedron_separation(&cloud_a, &cloud_b);
        if separation.gap > -0.02 {
            continue; // disjoint, or too shallow to be decisive
        }
        let hull_a = GjkConvexHull::new(cloud_a.clone());
        let hull_b = GjkConvexHull::new(cloud_b.clone());
        let result = gjk_epa_contact(&hull_a, &hull_b)
            .unwrap_or_else(|| panic!("case {case}: overlap {} found no contact", -separation.gap));
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
        assert_penetration_identity(&result, 1e-8);
        assert_contacts_on_extremal_faces(
            &result,
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
fn random_rotated_boxes_match_the_separating_axis_solution() {
    let mut rng = Rng(0xB0C5_0000_5EED_0011);
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
        let cloud_b = box_corners(rng.vec(-1.2, 1.2), half_b, rot_b);
        let separation = axis_separation(&cloud_a, &cloud_b, &box_axes(rot_a, rot_b));
        if separation.gap > -0.02 {
            continue;
        }
        let result = gjk_epa_contact(
            &GjkConvexHull::new(cloud_a.clone()),
            &GjkConvexHull::new(cloud_b.clone()),
        )
        .unwrap_or_else(|| panic!("case {case}: overlap {} found no contact", -separation.gap));
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
        assert_penetration_identity(&result, 1e-8);
        assert_contacts_on_extremal_faces(
            &result,
            &cloud_a,
            &cloud_b,
            1e-6,
            &format!("case {case}"),
        );
        verified += 1;
    }
    assert!(verified > 80, "only {verified} decisive overlaps");
}

// ── epa() entry point ────────────────────────────────────────────────────────

/// A tetrahedron of Minkowski support points along the four given directions.
fn seed(a: &dyn ConvexShape, b: &dyn ConvexShape, directions: [[f64; 3]; 4]) -> Simplex {
    Simplex {
        points: directions
            .iter()
            .map(|&d| minkowski_support_point(a, b, d))
            .collect(),
    }
}

#[test]
fn support_points_carry_the_shape_pair_that_built_them() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([1.75, 0.5, 0.0], [1.0; 3]);
    for direction in [[1.0, 0.3, -0.2], [-0.4, 1.0, 0.1], [0.2, -0.3, -1.0]] {
        let point = minkowski_support_point(&a, &b, direction);
        assert_eq!(point.a, a.support(direction));
        assert_eq!(point.b, b.support(neg3(direction)));
        assert_eq!(point.w, sub3(point.a, point.b));
        assert_eq!(point.w, minkowski_support(&a, &b, direction));
    }
}

#[test]
fn epa_grows_a_seed_that_does_not_contain_the_origin() {
    // The alternating-corner tetrahedron of this Minkowski difference misses the
    // origin; the expansion must grow toward it instead of trusting the seed.
    let a = ball([0.0; 3], 1.0);
    let b = ball([1.5, 0.0, 0.0], 1.0);
    let off_center = seed(
        &a,
        &b,
        [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ],
    );
    let result = epa(&a, &b, off_center).expect("the spheres overlap");
    assert!((result.depth - 0.5).abs() < 1e-6, "depth {}", result.depth);
    assert_vec_close(result.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert_penetration_identity(&result, 1e-9);
}

#[test]
fn epa_reports_no_penetration_for_a_seed_of_separated_shapes() {
    let a = ball([0.0; 3], 1.0);
    let b = ball([3.0, 0.0, 0.0], 1.0);
    let far = seed(
        &a,
        &b,
        [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ],
    );
    assert!(epa(&a, &b, far).is_none());
}

#[test]
fn epa_rejects_degenerate_seeds() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([1.5, 0.0, 0.0], [1.0; 3]);
    let point = minkowski_support_point(&a, &b, [1.0, 0.0, 0.0]);
    // Too few points.
    assert!(
        epa(
            &a,
            &b,
            Simplex {
                points: vec![point; 3]
            }
        )
        .is_none()
    );
    // Four coincident points and four coplanar points span no volume.
    assert!(
        epa(
            &a,
            &b,
            Simplex {
                points: vec![point; 4]
            }
        )
        .is_none()
    );
    let flat = Simplex {
        points: [
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
        ]
        .iter()
        .map(|&w| SupportPoint {
            w,
            a: w,
            b: [0.0; 3],
        })
        .collect(),
    };
    assert!(epa(&a, &b, flat).is_none());
}

// ── completion of a degenerate GJK simplex ───────────────────────────────────

/// Volume of the tetrahedron spanned by four points.
fn tetrahedron_volume(points: &[SupportPoint]) -> f64 {
    dot3(
        cross3(
            sub3(points[1].w, points[0].w),
            sub3(points[2].w, points[0].w),
        ),
        sub3(points[3].w, points[0].w),
    )
    .abs()
        / 6.0
}

#[test]
fn completion_builds_a_tetrahedron_from_any_degenerate_simplex() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([1.75, 0.0, 0.0], [1.0; 3]);
    let oracle = |d: [f64; 3]| minkowski_support_point(&a, &b, d);
    let known = oracle([1.0, 0.0, 0.0]);
    for simplex in [
        vec![],
        vec![known],
        vec![known, known],
        vec![known, known, known, known],
        vec![known, oracle([-1.0, 0.0, 0.0])],
    ] {
        let completed = complete_tetrahedron(&simplex, &oracle);
        assert_eq!(completed.len(), 4, "from {} points", simplex.len());
        assert!(
            tetrahedron_volume(&completed) > 1e-6,
            "from {} points: flat tetrahedron {completed:?}",
            simplex.len()
        );
        // Every completed vertex is a point of the Minkowski difference with the
        // support pair that produced it.
        for point in &completed {
            assert_vec_close(point.w, sub3(point.a, point.b), 1e-12);
        }
    }
}

#[test]
fn completion_keeps_a_valid_tetrahedron_untouched() {
    let a = ball([0.0; 3], 1.0);
    let b = ball([1.5, 0.0, 0.0], 1.0);
    let oracle = |d: [f64; 3]| minkowski_support_point(&a, &b, d);
    let known: Vec<SupportPoint> = [
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ]
    .iter()
    .map(|&d| oracle(d))
    .collect();
    assert_eq!(complete_tetrahedron(&known, &oracle), known);
}

#[test]
fn flat_shapes_have_no_volume_to_penetrate() {
    // Two zero-thickness squares in the plane z = 0 overlap in the plane, but
    // their Minkowski difference is flat: there is no 3-D penetration to find.
    let square = |center_x: f64| {
        GjkConvexHull::new(vec![
            [center_x - 1.0, -1.0, 0.0],
            [center_x + 1.0, -1.0, 0.0],
            [center_x + 1.0, 1.0, 0.0],
            [center_x - 1.0, 1.0, 0.0],
        ])
    };
    let a = square(0.0);
    let b = square(1.5);
    assert!(gjk_intersect(&a, &b));
    assert!(gjk_epa_contact(&a, &b).is_none());
    let oracle = |d: [f64; 3]| minkowski_support_point(&a, &b, d);
    assert!(complete_tetrahedron(&[], &oracle).len() < 4);
}

#[test]
fn touching_shapes_have_no_meaningful_penetration() {
    let a = ball([0.0; 3], 1.0);
    let b = ball([2.0, 0.0, 0.0], 1.0);
    if let Some(result) = gjk_epa_contact(&a, &b) {
        assert!(result.depth < 1e-6, "depth {}", result.depth);
    }
}
