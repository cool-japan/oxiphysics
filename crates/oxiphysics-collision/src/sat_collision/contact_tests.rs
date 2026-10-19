// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the SAT contact-point generator and the OBB-capsule
//! test.

use super::functions::{add, cross, dot, length, normalize, scale, sub};
use super::types::{
    Capsule, ContactFeature, Obb, ObbCapsuleCollision, ObbCollision, SatContact,
    SatContactPointGenerator,
};

const IDENTITY_AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Images of the unit vectors under the rotation by `angle` about `axis`
/// (Rodrigues' formula): the columns of the rotation matrix.
fn rotation_axes(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let k = normalize(axis);
    let (sin, cos) = angle.sin_cos();
    let rotate = |v: [f64; 3]| {
        add(
            add(scale(v, cos), scale(cross(k, v), sin)),
            scale(k, dot(k, v) * (1.0 - cos)),
        )
    };
    [
        rotate([1.0, 0.0, 0.0]),
        rotate([0.0, 1.0, 0.0]),
        rotate([0.0, 0.0, 1.0]),
    ]
}

fn obb(center: [f64; 3], half_extents: [f64; 3], axes: [[f64; 3]; 3]) -> Obb {
    Obb {
        center,
        axes,
        half_extents,
    }
}

/// Deterministic generator in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn vector(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }

    fn unit_vector(&mut self) -> [f64; 3] {
        loop {
            let v = self.vector(-1.0, 1.0);
            let len = length(v);
            if len > 0.1 && len <= 1.0 {
                return scale(v, 1.0 / len);
            }
        }
    }

    fn random_obb(&mut self, center: [f64; 3], half_lo: f64, half_hi: f64) -> Obb {
        let axes = rotation_axes(self.unit_vector(), self.range(0.0, std::f64::consts::PI));
        obb(center, self.vector(half_lo, half_hi), axes)
    }
}

fn assert_close(actual: [f64; 3], expected: [f64; 3], tol: f64, what: &str) {
    let d = length(sub(actual, expected));
    assert!(
        d <= tol,
        "{what}: {actual:?} differs from {expected:?} by {d}"
    );
}

/// Sort points lexicographically so that patches can be compared as sets.
fn sorted(mut points: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    points.sort_by(|a, b| {
        a[0].total_cmp(&b[0])
            .then(a[1].total_cmp(&b[1]))
            .then(a[2].total_cmp(&b[2]))
    });
    points
}

fn distance_to_obb(obb: &Obb, p: [f64; 3]) -> f64 {
    length(sub(p, ObbCapsuleCollision::closest_point_on_obb(obb, p)))
}

fn assert_patch(points: Vec<[f64; 3]>, expected: Vec<[f64; 3]>, what: &str) {
    assert_eq!(points.len(), expected.len(), "{what}: {points:?}");
    for (p, e) in sorted(points).into_iter().zip(sorted(expected)) {
        assert_close(p, e, 1e-9, what);
    }
}

// ── SatContactPointGenerator ────────────────────────────────────────────────

#[test]
fn face_contact_points_lie_halfway_through_the_overlap() {
    // The audit example: A = [-1, 1]^3, B centred at x = 1.5 with half-extent 1
    // overlaps A in x in [0.5, 1]. The old generator reported x = -0.25.
    let a = Obb::from_aabb([-1.0; 3], [1.0; 3]);
    let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
    let contact = ObbCollision::obb_obb_sat(&a, &b).expect("the boxes overlap");
    assert_close(contact.normal, [-1.0, 0.0, 0.0], 1e-12, "normal B -> A");
    assert!((contact.depth - 0.5).abs() < 1e-12);

    let expected = vec![
        [0.75, -1.0, -1.0],
        [0.75, -1.0, 1.0],
        [0.75, 1.0, -1.0],
        [0.75, 1.0, 1.0],
    ];
    assert_patch(
        SatContactPointGenerator::from_obb_obb(&a, &b, &contact),
        expected.clone(),
        "A as reference",
    );
    // The same pair with the roles exchanged gives the same patch.
    let swapped = ObbCollision::obb_obb_sat(&b, &a).expect("the boxes overlap");
    assert_close(swapped.normal, [1.0, 0.0, 0.0], 1e-12, "swapped normal");
    assert_patch(
        SatContactPointGenerator::from_obb_obb(&b, &a, &swapped),
        expected,
        "B first",
    );
}

#[test]
fn partially_overlapping_faces_are_clipped_to_the_shared_patch() {
    // B is shifted by +1 in y: the shared patch is y in [0, 1], z in [-1, 1].
    let a = Obb::from_aabb([-1.0; 3], [1.0; 3]);
    let b = Obb::from_aabb([0.5, 0.0, -1.0], [2.5, 2.0, 1.0]);
    let contact = ObbCollision::obb_obb_sat(&a, &b).expect("overlap");
    let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
    assert_patch(
        points,
        vec![
            [0.75, 0.0, -1.0],
            [0.75, 0.0, 1.0],
            [0.75, 1.0, -1.0],
            [0.75, 1.0, 1.0],
        ],
        "clipped patch",
    );
}

#[test]
fn small_box_resting_on_a_plate_reports_its_own_footprint() {
    let plate = obb([0.0; 3], [3.0, 0.5, 3.0], IDENTITY_AXES);
    let block = obb([0.5, 0.99, 0.2], [0.5, 0.5, 0.5], IDENTITY_AXES);
    let contact = ObbCollision::obb_obb_sat(&plate, &block).expect("overlap");
    assert_close(contact.normal, [0.0, -1.0, 0.0], 1e-12, "normal");
    // Reference: the plate top at y = 0.5; the block bottom at y = 0.49 is
    // 0.01 below it, so the points sit at y = 0.495.
    assert_patch(
        SatContactPointGenerator::from_obb_obb(&plate, &block, &contact),
        vec![
            [0.0, 0.495, -0.3],
            [0.0, 0.495, 0.7],
            [1.0, 0.495, -0.3],
            [1.0, 0.495, 0.7],
        ],
        "footprint",
    );
    // With the block as first argument the block is the reference box; the
    // plate's top face is clipped to the block's footprint: the same patch.
    let contact = ObbCollision::obb_obb_sat(&block, &plate).expect("overlap");
    assert_patch(
        SatContactPointGenerator::from_obb_obb(&block, &plate, &contact),
        vec![
            [0.0, 0.495, -0.3],
            [0.0, 0.495, 0.7],
            [1.0, 0.495, -0.3],
            [1.0, 0.495, 0.7],
        ],
        "footprint, roles exchanged",
    );
}

#[test]
fn turned_block_keeps_its_rotated_footprint_inside_the_plate() {
    // The block is turned by 30 degrees about the vertical axis: its bottom
    // face is still parallel to the plate but rotated; the footprint corners
    // are the corners of the rotated square.
    let angle = std::f64::consts::FRAC_PI_6;
    let plate = obb([0.0; 3], [3.0, 0.5, 3.0], IDENTITY_AXES);
    let block = obb(
        [0.0, 0.99, 0.0],
        [0.5, 0.5, 0.5],
        rotation_axes([0.0, 1.0, 0.0], angle),
    );
    let contact = ObbCollision::obb_obb_sat(&plate, &block).expect("overlap");
    // The cross product of two edges parallel to the faces equals the face
    // axis; the face axis must keep the contact (a flat face contact is not an
    // edge contact).
    assert!(
        matches!(contact.contact_type, ContactFeature::FaceFace { .. }),
        "{:?}",
        contact.contact_type
    );
    let points = SatContactPointGenerator::from_obb_obb(&plate, &block, &contact);
    assert_eq!(points.len(), 4, "{points:?}");
    let (sin, cos) = angle.sin_cos();
    let corner = |sx: f64, sz: f64| {
        let (lx, lz) = (0.5 * sx, 0.5 * sz);
        [lx * cos + lz * sin, 0.495, -lx * sin + lz * cos]
    };
    let expected = vec![
        corner(1.0, 1.0),
        corner(1.0, -1.0),
        corner(-1.0, 1.0),
        corner(-1.0, -1.0),
    ];
    assert_patch(points, expected, "rotated footprint");
}

#[test]
fn edge_contact_is_the_middle_of_the_closest_points_of_the_real_edges() {
    // Cube A turned 45 degrees about z (its top edge runs along z), cube B
    // turned 45 degrees about x (its bottom edge runs along x), pushed 0.02
    // into A. The edges cross at (0, ., 0).
    let half = 0.5;
    let corner = half * 2.0_f64.sqrt();
    let overlap = 0.02;
    let a = obb(
        [0.0; 3],
        [half; 3],
        rotation_axes([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
    );
    let b = obb(
        [0.0, 2.0 * corner - overlap, 0.0],
        [half; 3],
        rotation_axes([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
    );
    let contact = ObbCollision::obb_obb_sat(&a, &b).expect("the edges cross");
    assert!(
        matches!(contact.contact_type, ContactFeature::EdgeEdge { .. }),
        "{:?}",
        contact.contact_type
    );
    assert!((contact.depth - overlap).abs() < 1e-9, "{}", contact.depth);
    assert_close(contact.normal, [0.0, -1.0, 0.0], 1e-9, "normal");
    let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
    assert_eq!(points.len(), 1);
    // Top edge of A at y = corner, bottom edge of B at y = corner - overlap.
    assert_close(
        points[0],
        [0.0, corner - 0.5 * overlap, 0.0],
        1e-9,
        "edge contact",
    );
}

#[test]
fn edge_contact_uses_edges_of_full_length_not_unit_axes() {
    // Long boxes: the edge half-lengths (5 and 4) matter for where the
    // crossing point lies; with unit-length axes the closest points would
    // be clamped near the corners.
    let a = obb(
        [0.0; 3],
        [0.5, 0.5, 5.0],
        rotation_axes([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
    );
    let corner = 0.5 * 2.0_f64.sqrt();
    let b = obb(
        [3.0, 2.0 * corner - 0.04, 0.0],
        [4.0, 0.5, 0.5],
        rotation_axes([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
    );
    let contact = ObbCollision::obb_obb_sat(&a, &b).expect("the long edges cross");
    assert!(matches!(
        contact.contact_type,
        ContactFeature::EdgeEdge { .. }
    ));
    let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
    assert_eq!(points.len(), 1);
    // A's edge runs along z through x = 0 and B's along x through z = 0.
    assert_close(
        points[0],
        [0.0, corner - 0.02, 0.0],
        1e-9,
        "crossing point of the long edges",
    );
}

#[test]
fn other_face_features_use_the_same_clipping_and_a_zero_normal_is_recovered() {
    let a = Obb::from_aabb([-1.0; 3], [1.0; 3]);
    let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
    let base = ObbCollision::obb_obb_sat(&a, &b).expect("overlap");
    let expected = sorted(SatContactPointGenerator::from_obb_obb(&a, &b, &base));
    for feature in [
        ContactFeature::VertexFace { vertex: 3, face: 1 },
        ContactFeature::FaceEdge { face: 1, edge: 4 },
    ] {
        let contact = SatContact {
            normal: base.normal,
            depth: base.depth,
            contact_type: feature,
        };
        let points = sorted(SatContactPointGenerator::from_obb_obb(&a, &b, &contact));
        assert_eq!(points, expected);
    }
    // A zero normal falls back to the direction between the centres.
    let contact = SatContact {
        normal: [0.0; 3],
        depth: base.depth,
        contact_type: base.contact_type.clone(),
    };
    let points = sorted(SatContactPointGenerator::from_obb_obb(&a, &b, &contact));
    assert_eq!(points.len(), expected.len());
    for p in &points {
        assert!(p[0] >= 0.5 - 1e-9 && p[0] <= 1.0 + 1e-9, "{p:?}");
    }
}

#[test]
fn a_contact_without_a_clipped_patch_still_yields_one_point() {
    // A fabricated contact whose normal does not match the geometry: the
    // incident face lies beside the reference face, so nothing survives the
    // clipping and the generator falls back to the support-point midpoint.
    let a = Obb::from_aabb([-1.0; 3], [1.0; 3]);
    let b = Obb::from_aabb([5.0, -1.0, -1.0], [7.0, 1.0, 1.0]);
    let contact = SatContact {
        normal: [0.0, 1.0, 0.0],
        depth: 0.5,
        contact_type: ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        },
    };
    let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
    assert_eq!(points.len(), 1);
    let deepest = a.support([0.0, -1.0, 0.0]);
    assert_close(
        points[0],
        add(deepest, [0.0, 0.25, 0.0]),
        1e-12,
        "fallback point",
    );
}

/// Penetration depth of two boxes along the unit axis `n`:
/// `r_a(n) + r_b(n) - |(c_a - c_b) . n|`.
fn depth_along(a: &Obb, b: &Obb, n: [f64; 3]) -> f64 {
    let radius = |o: &Obb| -> f64 {
        (0..3)
            .map(|k| o.half_extents[k] * dot(o.axes[k], n).abs())
            .sum()
    };
    radius(a) + radius(b) - dot(sub(a.center, b.center), n).abs()
}

#[test]
fn contact_points_of_random_overlaps_stay_inside_the_overlap_zone() {
    let mut rng = Lcg(0x5A7C_0DE5_0000_0123);
    let mut face_cases = 0;
    let mut edge_cases = 0;
    for case in 0..400 {
        let a = rng.random_obb([0.0; 3], 0.3, 1.2);
        let offset = scale(rng.unit_vector(), rng.range(0.0, 1.4));
        let b = rng.random_obb(offset, 0.3, 1.2);
        let Some(contact) = ObbCollision::obb_obb_sat(&a, &b) else {
            continue;
        };
        let points = SatContactPointGenerator::from_obb_obb(&a, &b, &contact);
        assert!(!points.is_empty(), "case {case}: no contact points");
        assert!(points.len() <= 8, "case {case}: {} points", points.len());
        match contact.contact_type {
            ContactFeature::EdgeEdge { .. } => {
                edge_cases += 1;
                assert_eq!(points.len(), 1, "case {case}");
            }
            _ => face_cases += 1,
        }
        // The true penetration along the contact axis bounds how far a point
        // half-way through the overlap can be from either box.
        let depth = depth_along(&a, &b, contact.normal).max(0.0);
        for p in &points {
            assert!(
                distance_to_obb(&a, *p) <= depth + 1e-9,
                "case {case}: point {p:?} is {} from A (depth {depth})",
                distance_to_obb(&a, *p)
            );
            assert!(
                distance_to_obb(&b, *p) <= depth + 1e-9,
                "case {case}: point {p:?} is {} from B (depth {depth})",
                distance_to_obb(&b, *p)
            );
        }
    }
    assert!(face_cases > 50, "face cases: {face_cases}");
    assert!(edge_cases > 5, "edge cases: {edge_cases}");
}

// ── ObbCollision penetration depth ──────────────────────────────────────────

#[test]
fn penetration_on_axis_is_the_separating_translation_not_the_common_interval() {
    let x = [1.0, 0.0, 0.0];
    // Partially overlapping boxes: both measures agree.
    let a = Obb::from_aabb([-1.0; 3], [1.0; 3]);
    let b = Obb::from_aabb([0.5, -1.0, -1.0], [2.5, 1.0, 1.0]);
    let common = ObbCollision::test_overlap_on_axis(&a, &b, x).expect("overlap");
    let push = ObbCollision::penetration_on_axis(&a, &b, x).expect("overlap");
    assert!((common - 0.5).abs() < 1e-12 && (push - 0.5).abs() < 1e-12);
    // B lies within A along x: the common interval is B's own width 1.0, but
    // moving A until it clears B takes 3 (to the left) or 4 (to the right).
    let big = Obb::from_aabb([-3.0; 3], [3.0; 3]);
    let inner = obb([0.5, 0.0, 0.0], [0.5; 3], IDENTITY_AXES);
    let common = ObbCollision::test_overlap_on_axis(&big, &inner, x).expect("overlap");
    let push = ObbCollision::penetration_on_axis(&big, &inner, x).expect("overlap");
    assert!((common - 1.0).abs() < 1e-12, "common interval {common}");
    assert!((push - 3.0).abs() < 1e-12, "separating translation {push}");
    // The axis need not be normalised.
    let scaled = ObbCollision::penetration_on_axis(&big, &inner, [5.0, 0.0, 0.0]).expect("overlap");
    assert!((scaled - 3.0).abs() < 1e-12);
    // Separated boxes and a degenerate axis.
    let far = Obb::from_aabb([3.0, -1.0, -1.0], [5.0, 1.0, 1.0]);
    assert!(ObbCollision::penetration_on_axis(&a, &far, x).is_none());
    assert!(ObbCollision::penetration_on_axis(&a, &b, [0.0; 3]).is_none());
    // Touching boxes have zero penetration.
    let touching = Obb::from_aabb([1.0, -1.0, -1.0], [3.0, 1.0, 1.0]);
    let touch = ObbCollision::penetration_on_axis(&a, &touching, x).expect("touching");
    assert!(touch.abs() < 1e-12);
}

#[test]
fn sat_depth_of_a_box_inside_another_is_the_distance_to_get_clear() {
    // B (half 0.5) lies inside A (half 3), off-centre toward +x. Along every
    // axis B's interval is contained in A's; pushing B out takes 3 along x, the
    // shortest way, and 3.5 along y and z.
    let big = Obb::from_aabb([-3.0; 3], [3.0; 3]);
    let inner = obb([0.5, 0.0, 0.0], [0.5; 3], IDENTITY_AXES);
    let contact = ObbCollision::obb_obb_sat(&big, &inner).expect("B lies inside A");
    assert!(
        (contact.depth - 3.0).abs() < 1e-12,
        "depth {}",
        contact.depth
    );
    assert_close(contact.normal, [-1.0, 0.0, 0.0], 1e-12, "normal B -> A");
    // The contact points stay between the two surfaces.
    let points = SatContactPointGenerator::from_obb_obb(&big, &inner, &contact);
    assert!(!points.is_empty());
    for p in &points {
        assert!(distance_to_obb(&big, *p) < 1e-9, "{p:?} is inside A");
        assert!(distance_to_obb(&inner, *p) <= contact.depth, "{p:?}");
    }
}

// ── ObbCapsuleCollision ─────────────────────────────────────────────────────

/// Distance from the box to the point of the segment `p0 -> p1` at parameter `s`.
fn sampled_distance(obb: &Obb, p0: [f64; 3], p1: [f64; 3], s: f64) -> f64 {
    let p = add(p0, scale(sub(p1, p0), s));
    distance_to_obb(obb, p)
}

/// Minimum distance from the segment to the box by dense sampling followed by
/// a ternary refinement around the best sample.
fn brute_force_distance(obb: &Obb, p0: [f64; 3], p1: [f64; 3]) -> f64 {
    let samples = 4000;
    let mut best_i = 0;
    let mut best = f64::INFINITY;
    for i in 0..=samples {
        let d = sampled_distance(obb, p0, p1, f64::from(i) / f64::from(samples));
        if d < best {
            best = d;
            best_i = i;
        }
    }
    let step = 1.0 / f64::from(samples);
    let mut lo = (f64::from(best_i) * step - step).max(0.0);
    let mut hi = (f64::from(best_i) * step + step).min(1.0);
    for _ in 0..200 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if sampled_distance(obb, p0, p1, m1) < sampled_distance(obb, p0, p1, m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    best.min(sampled_distance(obb, p0, p1, 0.5 * (lo + hi)))
}

#[test]
fn audit_counterexample_the_capsule_nearest_the_centre_is_not_the_nearest_to_the_box() {
    // OBB half (10, 0.1, 0.1) at the origin; the capsule axis runs from
    // (-20, 1, 0) to (20, 0.2, 0) with radius 0.4. Its point nearest the centre
    // is 0.4998 from the box, but near x = 10 the axis is only 0.29994 away.
    let long_box = obb([0.0; 3], [10.0, 0.1, 0.1], IDENTITY_AXES);
    let capsule = Capsule::new([-20.0, 1.0, 0.0], [20.0, 0.2, 0.0], 0.4);
    let (depth, contact, normal) =
        ObbCapsuleCollision::test(&long_box, &capsule).expect("the capsule touches the box");

    let exact = brute_force_distance(&long_box, capsule.p0, capsule.p1);
    assert!((exact - 0.29994).abs() < 1e-5, "reference distance {exact}");
    assert!(
        (depth - (0.4 - exact)).abs() < 1e-8,
        "depth {depth} vs {}",
        0.4 - exact
    );
    assert!((length(normal) - 1.0).abs() < 1e-12);
    assert!(
        normal[1] > 0.99 && normal[0] > 0.0 && normal[0] < 0.05,
        "{normal:?}"
    );
    // The contact is at the box end x = 10, half-way into the penetration.
    assert!((contact[0] - 10.0).abs() < 0.05, "contact {contact:?}");
    assert!(
        contact[1] > 0.1 - depth && contact[1] < 0.1 + 1e-9,
        "{contact:?}"
    );
    // The old point (nearest to the centre) was 0.4998 away and rejected.
    let (axis_point, box_point) =
        ObbCapsuleCollision::closest_points_segment_obb(&long_box, capsule.p0, capsule.p1);
    assert!((length(sub(axis_point, box_point)) - exact).abs() < 1e-9);
}

#[test]
fn segment_box_closest_points_match_dense_sampling_for_random_inputs() {
    let mut rng = Lcg(0xFEED_FACE_CAFE_0001);
    for case in 0..200 {
        let centre = rng.vector(-0.5, 0.5);
        let b = rng.random_obb(centre, 0.2, 1.5);
        let p0 = rng.vector(-3.0, 3.0);
        let p1 = rng.vector(-3.0, 3.0);
        let (on_segment, on_box) = ObbCapsuleCollision::closest_points_segment_obb(&b, p0, p1);
        let exact = length(sub(on_segment, on_box));
        // The reported box point is the closest box point to the segment point.
        assert_close(
            on_box,
            ObbCapsuleCollision::closest_point_on_obb(&b, on_segment),
            1e-9,
            "box point",
        );
        // The segment point lies on the segment.
        let along = sub(on_segment, p0);
        let dir = sub(p1, p0);
        let t = dot(along, dir) / dot(dir, dir);
        assert!((-1e-12..=1.0 + 1e-12).contains(&t), "case {case}: t = {t}");
        assert_close(on_segment, add(p0, scale(dir, t)), 1e-9, "on segment");
        let reference = brute_force_distance(&b, p0, p1);
        assert!(
            exact <= reference + 1e-12,
            "case {case}: exact {exact} above the sampled distance {reference}"
        );
        assert!(
            (exact - reference).abs() < 1e-8,
            "case {case}: exact {exact} vs reference {reference}"
        );
    }
}

#[test]
fn segment_crossing_the_box_has_coincident_closest_points() {
    let unit = obb([0.0; 3], [1.0, 1.0, 1.0], IDENTITY_AXES);
    let (on_segment, on_box) =
        ObbCapsuleCollision::closest_points_segment_obb(&unit, [0.5, 0.2, -3.0], [0.5, 0.2, 3.0]);
    assert_close(on_segment, on_box, 1e-12, "crossing segment");
    // A segment that ends inside the box.
    let (on_segment, on_box) =
        ObbCapsuleCollision::closest_points_segment_obb(&unit, [0.0, 0.0, 0.0], [5.0, 0.0, 0.0]);
    assert_close(on_segment, on_box, 1e-12, "segment starting inside");
}

#[test]
fn capsule_lying_over_a_face_reports_the_face_normal_at_the_middle_of_the_contact() {
    let unit = obb([0.0; 3], [1.0, 1.0, 1.0], IDENTITY_AXES);
    let capsule = Capsule::new([-0.5, 1.3, 0.0], [0.5, 1.3, 0.0], 0.5);
    let (depth, contact, normal) = ObbCapsuleCollision::test(&unit, &capsule).expect("contact");
    assert!((depth - 0.2).abs() < 1e-12, "depth {depth}");
    assert_close(normal, [0.0, 1.0, 0.0], 1e-12, "normal");
    // The box surface is at y = 1 and the capsule bottom at y = 0.8.
    assert_close(contact, [0.0, 0.9, 0.0], 1e-12, "contact");
}

#[test]
fn capsule_end_cap_against_a_box_corner() {
    let unit = obb([0.0; 3], [1.0, 1.0, 1.0], IDENTITY_AXES);
    let capsule = Capsule::new([1.5, 1.5, 1.5], [3.0, 3.0, 3.0], 1.0);
    let (depth, contact, normal) = ObbCapsuleCollision::test(&unit, &capsule).expect("contact");
    let gap = 0.75_f64.sqrt();
    assert!((depth - (1.0 - gap)).abs() < 1e-12, "depth {depth}");
    let diagonal = 1.0 / 3.0_f64.sqrt();
    assert_close(normal, [diagonal; 3], 1e-12, "normal");
    // Half-way between the corner (1, 1, 1) and the capsule's deepest point.
    let expected = sub([1.0; 3], scale([diagonal; 3], 0.5 * depth));
    assert_close(contact, expected, 1e-12, "contact");
    // Out of reach: the end cap stays farther than the radius.
    let far = Capsule::new([1.7, 1.7, 1.7], [3.0, 3.0, 3.0], 1.0);
    assert!(ObbCapsuleCollision::test(&unit, &far).is_none());
}

#[test]
fn capsule_axis_through_the_box_uses_the_minimum_translation_axis() {
    let unit = obb([0.0; 3], [1.0, 1.0, 1.0], IDENTITY_AXES);
    // The axis crosses the box close to its +x face: pushing out along +x is
    // the shortest way (1 + 0.3 - 0.8).
    let capsule = Capsule::new([0.8, 0.0, -3.0], [0.8, 0.0, 3.0], 0.3);
    let (depth, contact, normal) = ObbCapsuleCollision::test(&unit, &capsule).expect("contact");
    assert!((depth - 0.5).abs() < 1e-12, "depth {depth}");
    assert_close(normal, [1.0, 0.0, 0.0], 1e-12, "normal");
    // The middle of the part of the axis inside the box.
    assert_close(contact, [0.8, 0.0, 0.0], 1e-12, "contact");
    // Through the centre both side axes tie; the depth is the same either way.
    let centred = Capsule::new([0.0, 0.0, -3.0], [0.0, 0.0, 3.0], 0.5);
    let (depth, _, normal) = ObbCapsuleCollision::test(&unit, &centred).expect("contact");
    assert!((depth - 1.5).abs() < 1e-12, "depth {depth}");
    assert!(normal[2].abs() < 1e-12 && (length(normal) - 1.0).abs() < 1e-12);
}

#[test]
fn zero_length_capsule_behaves_like_a_sphere() {
    let unit = obb([0.0; 3], [1.0, 1.0, 1.0], IDENTITY_AXES);
    let outside = Capsule::new([1.2, 0.0, 0.0], [1.2, 0.0, 0.0], 0.5);
    let (depth, contact, normal) = ObbCapsuleCollision::test(&unit, &outside).expect("contact");
    assert!((depth - 0.3).abs() < 1e-12);
    assert_close(normal, [1.0, 0.0, 0.0], 1e-12, "normal");
    assert_close(contact, [0.85, 0.0, 0.0], 1e-12, "contact");
    // Centre inside: leave through the nearest face, the whole radius included.
    let inside = Capsule::new([0.5, 0.2, 0.0], [0.5, 0.2, 0.0], 0.5);
    let (depth, _, normal) = ObbCapsuleCollision::test(&unit, &inside).expect("contact");
    assert!((depth - 1.0).abs() < 1e-12, "depth {depth}");
    assert_close(normal, [1.0, 0.0, 0.0], 1e-12, "normal");
    // A sphere that only touches the surface does not penetrate.
    let touching = Capsule::new([1.5, 0.0, 0.0], [1.5, 0.0, 0.0], 0.5);
    assert!(ObbCapsuleCollision::test(&unit, &touching).is_none());
}

/// Smallest translation along `direction` after which the capsule no longer
/// penetrates the box (the box and the capsule overlap on an interval of
/// translations containing zero, so the end of that interval is bracketed and
/// bisected).
fn translation_to_separate(
    b: &Obb,
    p0: [f64; 3],
    p1: [f64; 3],
    radius: f64,
    direction: [f64; 3],
) -> f64 {
    let penetrates = |t: f64| {
        let shift = scale(direction, t);
        let (on_segment, on_box) =
            ObbCapsuleCollision::closest_points_segment_obb(b, add(p0, shift), add(p1, shift));
        length(sub(on_segment, on_box)) < radius
    };
    if !penetrates(0.0) {
        return 0.0;
    }
    let mut lo = 0.0;
    let mut hi = 1.0;
    while penetrates(hi) {
        lo = hi;
        hi *= 2.0;
        if hi > 1.0e6 {
            return f64::INFINITY;
        }
    }
    for _ in 0..70 {
        let mid = 0.5 * (lo + hi);
        if penetrates(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

#[test]
fn reported_depth_and_normal_are_the_minimum_translation_for_random_capsules() {
    let mut rng = Lcg(0x00DD_BA11_5EED_0042);
    let mut checked = 0;
    let mut crossing = 0;
    for case in 0..160 {
        let b = rng.random_obb([0.0; 3], 0.3, 1.5);
        let p0 = rng.vector(-2.0, 2.0);
        let p1 = rng.vector(-2.0, 2.0);
        let radius = rng.range(0.1, 0.9);
        let capsule = Capsule::new(p0, p1, radius);
        let (on_segment, on_box) = ObbCapsuleCollision::closest_points_segment_obb(&b, p0, p1);
        let result = ObbCapsuleCollision::test(&b, &capsule);
        let distance = length(sub(on_segment, on_box));
        assert_eq!(
            result.is_some(),
            distance < radius,
            "case {case}: overlap decision"
        );
        let Some((depth, contact, normal)) = result else {
            continue;
        };
        checked += 1;
        if distance < 1e-9 {
            crossing += 1;
        }
        assert!(
            (length(normal) - 1.0).abs() < 1e-9,
            "case {case}: {normal:?}"
        );
        // Moving along the normal by the depth just separates the shapes ...
        let along_normal = translation_to_separate(&b, p0, p1, radius, normal);
        assert!(
            (along_normal - depth).abs() < 1e-7,
            "case {case}: separating move {along_normal} vs depth {depth}"
        );
        // ... and no other direction does better.
        let mut shortest = f64::INFINITY;
        for _ in 0..600 {
            let u = rng.unit_vector();
            shortest = shortest.min(translation_to_separate(&b, p0, p1, radius, u));
        }
        assert!(
            shortest >= depth - 1e-7,
            "case {case}: direction with translation {shortest} beats depth {depth}"
        );
        // The contact point is in the box when the axis crosses it and within
        // half the depth of the box surface otherwise.
        let allowed = if distance < 1e-9 {
            1e-9
        } else {
            0.5 * depth + 1e-9
        };
        assert!(
            distance_to_obb(&b, contact) <= allowed,
            "case {case}: contact {contact:?}"
        );
    }
    assert!(checked >= 40, "only {checked} overlapping cases");
    assert!(crossing >= 5, "only {crossing} crossing cases");
}

#[test]
fn separated_capsule_returns_none_and_rotated_box_is_handled() {
    let turned = obb(
        [0.0; 3],
        [1.0, 0.5, 2.0],
        rotation_axes([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
    );
    // A capsule axis along x at y = 3: the turned box reaches y = 0.75 * sqrt(2)
    // = 1.06 at most, so radius 1.0 does not reach it and radius 2.5 does.
    let high = Capsule::new([-2.0, 3.0, 0.0], [2.0, 3.0, 0.0], 1.0);
    assert!(ObbCapsuleCollision::test(&turned, &high).is_none());
    let reaching = Capsule::new([-2.0, 3.0, 0.0], [2.0, 3.0, 0.0], 2.5);
    let (depth, _, normal) = ObbCapsuleCollision::test(&turned, &reaching).expect("contact");
    assert!((depth - (2.5 - (3.0 - 0.75 * 2.0_f64.sqrt()))).abs() < 1e-9);
    assert_close(normal, [0.0, 1.0, 0.0], 1e-9, "normal");
}
