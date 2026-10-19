// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the contact convention of the plain-array narrow phase:
//! `normal` points from body B toward body A, the witnesses lie on the two
//! surfaces and `point_b - point_a == normal * depth`.

use super::*;

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

fn vec_approx_eq(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    approx_eq(a[0], b[0], tol) && approx_eq(a[1], b[1], tol) && approx_eq(a[2], b[2], tol)
}

// ── sphere_sphere / capsule_capsule ──────────────────────────────────────────

#[test]
fn test_sphere_sphere_normal_points_from_b_toward_a() {
    // B sits on the +x side of A: the normal from B toward A is -x, the
    // witness of A is on A's surface facing B and vice versa.
    match sphere_sphere([0.0; 3], 1.0, [1.5, 0.0, 0.0], 0.75) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(
                vec_approx_eq(c.normal, [-1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.normal
            );
            assert!(approx_eq(c.depth, 0.25, 1e-12), "depth {}", c.depth);
            assert!(
                vec_approx_eq(c.point_a, [1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.point_a
            );
            assert!(
                vec_approx_eq(c.point_b, [0.75, 0.0, 0.0], 1e-12),
                "{:?}",
                c.point_b
            );
            assert!(vec_approx_eq(
                sub(c.point_b, c.point_a),
                scale(c.normal, c.depth),
                1e-12
            ));
        }
        other => panic!("expected a contact, got {other:?}"),
    }
    // Coincident centres: +Y by convention.
    match sphere_sphere([2.0; 3], 1.0, [2.0; 3], 0.5) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(vec_approx_eq(c.normal, [0.0, 1.0, 0.0], 1e-12));
            assert!(approx_eq(c.depth, 1.5, 1e-12));
        }
        other => panic!("expected a contact, got {other:?}"),
    }
}

#[test]
fn test_capsule_capsule_normal_points_from_b_toward_a() {
    // Parallel capsules, B at +x: the normal from B toward A is -x.
    match capsule_capsule(
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        0.5,
        [0.8, -1.0, 0.0],
        [0.8, 1.0, 0.0],
        0.5,
    ) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(
                vec_approx_eq(c.normal, [-1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.normal
            );
            assert!(approx_eq(c.depth, 0.2, 1e-12), "depth {}", c.depth);
            assert!(approx_eq(c.point_a[0], 0.5, 1e-12));
            assert!(approx_eq(c.point_b[0], 0.3, 1e-12));
        }
        other => panic!("expected a contact, got {other:?}"),
    }
}

// ── sphere_box ───────────────────────────────────────────────────────────────

#[test]
fn test_sphere_box_witnesses_lie_on_the_two_surfaces() {
    // Centre outside: the sphere presses into the +x face.
    match sphere_box([1.3, 0.2, -0.1], 0.5, [0.0; 3], [1.0; 3]) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(
                vec_approx_eq(c.normal, [1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.normal
            );
            assert!(approx_eq(c.depth, 0.2, 1e-12), "depth {}", c.depth);
            assert!(
                vec_approx_eq(c.point_a, [0.8, 0.2, -0.1], 1e-12),
                "{:?}",
                c.point_a
            );
            assert!(
                vec_approx_eq(c.point_b, [1.0, 0.2, -0.1], 1e-12),
                "{:?}",
                c.point_b
            );
        }
        other => panic!("expected a contact, got {other:?}"),
    }
    // Centre inside, nearest face +x: the box witness is on that face, not at
    // the sphere centre, and the gap relation holds.
    match sphere_box([0.2, 0.0, 0.0], 0.3, [0.0; 3], [1.0; 3]) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(
                vec_approx_eq(c.normal, [1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.normal
            );
            assert!(approx_eq(c.depth, 1.1, 1e-12), "depth {}", c.depth);
            assert!(
                vec_approx_eq(c.point_a, [-0.1, 0.0, 0.0], 1e-12),
                "{:?}",
                c.point_a
            );
            assert!(
                vec_approx_eq(c.point_b, [1.0, 0.0, 0.0], 1e-12),
                "{:?}",
                c.point_b
            );
            assert!(vec_approx_eq(
                sub(c.point_b, c.point_a),
                scale(c.normal, c.depth),
                1e-12
            ));
        }
        other => panic!("expected a contact, got {other:?}"),
    }
    // Centre inside, nearest face -z of an off-centre box.
    match sphere_box([2.0, 3.0, 4.8], 0.5, [2.0, 3.0, 5.0], [1.0, 1.0, 0.4]) {
        NarrowPhaseResult::SingleContact(c) => {
            assert!(
                vec_approx_eq(c.normal, [0.0, 0.0, -1.0], 1e-12),
                "{:?}",
                c.normal
            );
            assert!(
                vec_approx_eq(c.point_b, [2.0, 3.0, 4.6], 1e-12),
                "{:?}",
                c.point_b
            );
            assert!(approx_eq(c.depth, 0.2 + 0.5, 1e-12), "depth {}", c.depth);
        }
        other => panic!("expected a contact, got {other:?}"),
    }
}

// ── box_box_sat ──────────────────────────────────────────────────────────────

/// Local axes of the frame rotated by `angle` about `axis`, one per entry:
/// `axes[i]` is the world image of local axis `i`.
fn rotated_axes(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
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

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn box_contact(result: NarrowPhaseResult) -> Contact {
    match result {
        NarrowPhaseResult::SingleContact(c) => c,
        other => panic!("expected a single contact, got {other:?}"),
    }
}

/// `point_b - point_a == normal * depth`.
fn assert_witness_gap(c: &Contact) {
    let gap = sub(c.point_b, c.point_a);
    assert!(
        vec_approx_eq(gap, scale(c.normal, c.depth), 1e-9),
        "witness gap {gap:?} vs normal * depth {:?}",
        scale(c.normal, c.depth)
    );
}

#[test]
fn test_box_box_normal_points_from_b_toward_a_and_witnesses_lie_on_the_surfaces() {
    let half = [1.0f64, 1.0, 1.0];
    let c = box_contact(box_box_sat(
        [0.0; 3],
        half,
        IDENTITY,
        [1.5, 0.0, 0.0],
        half,
        IDENTITY,
    ));
    // B sits on the +x side of A: the normal from B toward A is -x.
    assert!(
        vec_approx_eq(c.normal, [-1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        c.normal
    );
    assert!(approx_eq(c.depth, 0.5, 1e-12), "depth {}", c.depth);
    // A's +x face and B's -x face, not points near the box centres.
    assert!(
        vec_approx_eq(c.point_a, [1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        c.point_a
    );
    assert!(
        vec_approx_eq(c.point_b, [0.5, 0.0, 0.0], 1e-12),
        "{:?}",
        c.point_b
    );
    assert_witness_gap(&c);

    // Exchanging the boxes mirrors the contact.
    let m = box_contact(box_box_sat(
        [1.5, 0.0, 0.0],
        half,
        IDENTITY,
        [0.0; 3],
        half,
        IDENTITY,
    ));
    assert!(vec_approx_eq(m.normal, [1.0, 0.0, 0.0], 1e-12));
    assert!(vec_approx_eq(m.point_a, c.point_b, 1e-12));
    assert!(vec_approx_eq(m.point_b, c.point_a, 1e-12));
    assert!(approx_eq(m.depth, c.depth, 1e-12));
}

#[test]
fn test_box_box_local_axes_are_read_as_rows_of_the_array() {
    // B is turned by +30 degrees about z; `rot_b[i]` is the world direction
    // of its local axis i. Reading the array as a row-major matrix would
    // turn B by -30 degrees instead and move the witness to y = 0.367.
    let angle = std::f64::consts::FRAC_PI_6;
    let c = box_contact(box_box_sat(
        [0.0; 3],
        [1.0, 1.0, 1.0],
        IDENTITY,
        [1.6, 0.3, 0.0],
        [1.0, 0.5, 0.9],
        rotated_axes([0.0, 0.0, 1.0], angle),
    ));
    // B's corner nearest A's +x face is c - 1.0 b0 + 0.5 b1 = (0.484, 0.233).
    let depth = 1.0 + (1.0 * angle.cos() + 0.5 * angle.sin()) - 1.6;
    assert!(
        approx_eq(c.depth, depth, 1e-9),
        "depth {} vs {depth}",
        c.depth
    );
    assert!(
        vec_approx_eq(c.normal, [-1.0, 0.0, 0.0], 1e-9),
        "{:?}",
        c.normal
    );
    let corner_y = 0.3 - 0.5 + 0.5 * angle.cos();
    assert!(
        vec_approx_eq(c.point_b, [1.6 - angle.cos() - 0.25, corner_y, 0.0], 1e-9),
        "point_b {:?}",
        c.point_b
    );
    assert!(
        vec_approx_eq(c.point_a, [1.0, corner_y, 0.0], 1e-9),
        "point_a {:?}",
        c.point_a
    );
    assert_witness_gap(&c);
}

#[test]
fn test_box_box_edge_edge_depth_is_a_true_distance() {
    // Cube A turned 45 degrees about z (its top edge runs along z), cube B
    // turned 45 degrees about x (its bottom edge runs along x), pushed
    // 0.02 into A. The edge axis cross(a2, b0) is a unit vector here, but
    // a skew crossing at another angle must not scale the depth either.
    let half = [0.5f64; 3];
    let corner = 0.5 * 2.0_f64.sqrt();
    let overlap = 0.02;
    let c = box_contact(box_box_sat(
        [0.0; 3],
        half,
        rotated_axes([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
        [0.0, 2.0 * corner - overlap, 0.0],
        half,
        rotated_axes([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_4),
    ));
    assert!(approx_eq(c.depth, overlap, 1e-9), "depth {}", c.depth);
    assert!(
        vec_approx_eq(c.normal, [0.0, -1.0, 0.0], 1e-9),
        "{:?}",
        c.normal
    );
    assert!(
        vec_approx_eq(c.point_a, [0.0, corner, 0.0], 1e-9),
        "{:?}",
        c.point_a
    );
    assert!(
        vec_approx_eq(c.point_b, [0.0, corner - overlap, 0.0], 1e-9),
        "{:?}",
        c.point_b
    );
    assert_witness_gap(&c);
}

/// Minimum overlap over the 15 separating axes, computed independently of
/// the implementation: `(least face overlap, least edge overlap)` of the
/// unit axes, `None` when an axis separates the boxes.
fn oracle_overlaps(
    center_a: [f64; 3],
    half_a: [f64; 3],
    axes_a: [[f64; 3]; 3],
    center_b: [f64; 3],
    half_b: [f64; 3],
    axes_b: [[f64; 3]; 3],
) -> Option<(f64, f64)> {
    let offset = sub(center_b, center_a);
    let overlap = |l: [f64; 3]| -> f64 {
        let radius = |half: [f64; 3], axes: [[f64; 3]; 3]| -> f64 {
            (0..3).map(|k| half[k] * dot(axes[k], l).abs()).sum()
        };
        radius(half_a, axes_a) + radius(half_b, axes_b) - dot(offset, l).abs()
    };
    let mut faces = f64::INFINITY;
    for l in axes_a.iter().chain(axes_b.iter()) {
        faces = faces.min(overlap(*l));
    }
    let mut edges = f64::INFINITY;
    for a in axes_a {
        for b in axes_b {
            let c = cross(a, b);
            if len(c) > 1e-6 {
                edges = edges.min(overlap(normalize(c)));
            }
        }
    }
    (faces.min(edges) >= 0.0).then_some((faces, edges))
}

#[test]
fn test_box_box_depth_is_the_euclidean_overlap_for_random_rotations() {
    // Random rotations make most edge pairs cross at skew angles, where the
    // raw cross product of the two edge directions has a length below one
    // and an unnormalised depth is off by that factor.
    let mut state = 0xB0C5_5A71_2345_6789u64;
    let mut next = move |lo: f64, hi: f64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + (hi - lo) * ((state >> 11) as f64 / (1u64 << 53) as f64)
    };
    let mut overlapping = 0;
    let mut clear_edge_wins = 0;
    for case in 0..400 {
        let axis_a = [next(-1.0, 1.0), next(-1.0, 1.0), next(-1.0, 1.0)];
        let axis_b = [next(-1.0, 1.0), next(-1.0, 1.0), next(-1.0, 1.0)];
        if len(axis_a) < 0.1 || len(axis_b) < 0.1 {
            continue;
        }
        let axes_a = rotated_axes(axis_a, next(0.0, 3.0));
        let axes_b = rotated_axes(axis_b, next(0.0, 3.0));
        let half_a = [next(0.3, 1.2), next(0.3, 1.2), next(0.3, 1.2)];
        let half_b = [next(0.3, 1.2), next(0.3, 1.2), next(0.3, 1.2)];
        let center_b = [next(-1.5, 1.5), next(-1.5, 1.5), next(-1.5, 1.5)];
        let oracle = oracle_overlaps([0.0; 3], half_a, axes_a, center_b, half_b, axes_b);
        let result = box_box_sat([0.0; 3], half_a, axes_a, center_b, half_b, axes_b);
        let Some((faces, edges)) = oracle else {
            assert!(
                matches!(result, NarrowPhaseResult::NoContact),
                "case {case}: an axis separates the boxes"
            );
            continue;
        };
        overlapping += 1;
        let c = box_contact(result);
        let least = faces.min(edges);
        // Never deeper than the best face axis, never shallower than the
        // true minimum over all axes.
        assert!(
            c.depth <= faces + 1e-9 && c.depth >= least - 1e-9,
            "case {case}: depth {} outside [{least}, {faces}]",
            c.depth
        );
        if edges < faces - 1e-3 {
            clear_edge_wins += 1;
            // An edge axis has to beat the running best by a small bias
            // (1e-4), so a near-tie between two edge axes may keep the
            // first one.
            assert!(
                c.depth >= edges - 1e-9 && c.depth <= edges + 2e-4,
                "case {case}: edge axis wins with {edges}, got {}",
                c.depth
            );
        }
        assert!(approx_eq(len(c.normal), 1.0, 1e-9), "case {case}");
        assert!(
            dot(c.normal, sub([0.0; 3], center_b)) >= -1e-9,
            "case {case}: the normal must point from B toward A"
        );
        // The witnesses are one depth apart along the normal (up to the
        // same near-tie slack).
        let gap = dot(sub(c.point_b, c.point_a), c.normal);
        assert!(
            approx_eq(gap, c.depth, 1e-3),
            "case {case}: gap {gap} vs depth {}",
            c.depth
        );
    }
    assert!(overlapping > 100, "overlapping cases: {overlapping}");
    assert!(clear_edge_wins >= 5, "edge cases: {clear_edge_wins}");
}

#[test]
fn test_box_box_resting_box_reports_the_centre_of_the_contact_patch() {
    // Upper cube A (half 0.5) overhangs the edge of the plate B (half
    // 1, 0.5, 1): the shared patch is x in [0.1, 1.0], z in [-0.5, 0.5].
    let c = box_contact(box_box_sat(
        [0.6, 0.99, 0.0],
        [0.5, 0.5, 0.5],
        IDENTITY,
        [0.0; 3],
        [1.0, 0.5, 1.0],
        IDENTITY,
    ));
    assert!(approx_eq(c.depth, 0.01, 1e-9), "depth {}", c.depth);
    assert!(vec_approx_eq(c.normal, [0.0, 1.0, 0.0], 1e-12));
    assert!(
        vec_approx_eq(c.point_a, [0.55, 0.49, 0.0], 1e-9),
        "{:?}",
        c.point_a
    );
    assert!(
        vec_approx_eq(c.point_b, [0.55, 0.5, 0.0], 1e-9),
        "{:?}",
        c.point_b
    );
    assert_witness_gap(&c);
}

#[test]
fn test_box_box_cube_turned_45_degrees_resting_on_a_plate_picks_the_face_axis() {
    // A cube turned 45 degrees about z rests with its lowest edge (along z)
    // 0.05 below the top face of the plate. The edge-edge axis of that edge
    // and the plate's x edges coincides with the plate's face axis; the face
    // axis must win and the witnesses are surface points: the lowest edge of
    // the cube and the plate's top face directly above it.
    let corner = 0.5 * 2.0_f64.sqrt();
    let c = box_contact(box_box_sat(
        [0.0, 0.5 + corner - 0.05, 0.0],
        [0.5; 3],
        rotated_axes([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4),
        [0.0; 3],
        [2.0, 0.5, 2.0],
        IDENTITY,
    ));
    assert!(approx_eq(c.depth, 0.05, 1e-9), "depth {}", c.depth);
    assert!(
        vec_approx_eq(c.normal, [0.0, 1.0, 0.0], 1e-9),
        "{:?}",
        c.normal
    );
    assert!(
        vec_approx_eq(c.point_a, [0.0, 0.45, 0.0], 1e-9),
        "{:?}",
        c.point_a
    );
    assert!(
        vec_approx_eq(c.point_b, [0.0, 0.5, 0.0], 1e-9),
        "{:?}",
        c.point_b
    );
    assert_witness_gap(&c);
}

#[test]
fn test_box_box_touching_faces_have_zero_depth() {
    let half = [1.0f64; 3];
    let c = box_contact(box_box_sat(
        [0.0; 3],
        half,
        IDENTITY,
        [2.0, 0.0, 0.0],
        half,
        IDENTITY,
    ));
    assert!(approx_eq(c.depth, 0.0, 1e-12), "depth {}", c.depth);
    assert!(vec_approx_eq(c.normal, [-1.0, 0.0, 0.0], 1e-12));
    assert!(vec_approx_eq(c.point_a, c.point_b, 1e-12));
}

#[test]
fn test_box_box_degenerate_input_is_no_contact() {
    let half = [1.0f64; 3];
    let zero_axes = [[0.0; 3]; 3];
    assert!(matches!(
        box_box_sat([0.0; 3], half, zero_axes, [0.5, 0.0, 0.0], half, IDENTITY),
        NarrowPhaseResult::NoContact
    ));
    assert!(matches!(
        box_box_sat(
            [f64::NAN, 0.0, 0.0],
            half,
            IDENTITY,
            [0.5, 0.0, 0.0],
            half,
            IDENTITY
        ),
        NarrowPhaseResult::NoContact
    ));
}
