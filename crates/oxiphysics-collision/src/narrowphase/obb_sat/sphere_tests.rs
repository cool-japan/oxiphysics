// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`obb_sphere_test`].
//!
//! The contact normal points from the box towards the sphere centre, the
//! penetration depth is the distance the sphere has to travel along it to get
//! clear, and the contact point is the midpoint between the box-surface point
//! the sphere penetrates through and the facing point of the sphere surface.

use super::*;

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], what: &str) {
    let gap = (0..3)
        .map(|i| (actual[i] - expected[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(gap < 1e-12, "{what}: {actual:?} differs from {expected:?}");
}

#[test]
fn centre_inside_clears_the_nearest_face_plus_the_radius() {
    let obb = Obb::axis_aligned([0.0; 3], [1.0, 1.0, 1.0]);
    let result = obb_sphere_test(&obb, [0.0; 3], 0.1).expect("centre inside the box");
    assert!(
        (result.penetration_depth - 1.1).abs() < 1e-12,
        "depth {}",
        result.penetration_depth
    );
    let n = result.contact_normal;
    assert!(((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt() - 1.0).abs() < 1e-12);
    // Box surface point: 1.0 along the normal; sphere surface point: -0.1 along
    // it; the contact point is their midpoint, 0.45 along the normal.
    let along = result.contact_point[0] * n[0]
        + result.contact_point[1] * n[1]
        + result.contact_point[2] * n[2];
    assert!(
        (along - 0.45).abs() < 1e-12,
        "contact point {:?}",
        result.contact_point
    );
}

#[test]
fn off_centre_inside_uses_the_nearest_face() {
    let obb = Obb::axis_aligned([5.0, 5.0, 5.0], [1.0, 2.0, 3.0]);
    // Local (0.3, 1.8, -0.5): 0.2 below the +y face.
    let centre = [5.3, 6.8, 4.5];
    let result = obb_sphere_test(&obb, centre, 0.4).expect("centre inside the box");
    assert!((result.penetration_depth - 0.6).abs() < 1e-12);
    assert_vec_close(result.contact_normal, [0.0, 1.0, 0.0], "normal");
    // Face point (5.3, 7.0, 4.5), facing sphere point (5.3, 6.4, 4.5).
    assert_vec_close(result.contact_point, [5.3, 6.7, 4.5], "contact point");
}

#[test]
fn rotated_box_centre_inside_uses_the_nearest_rotated_face() {
    let (c, s) = (
        std::f64::consts::FRAC_PI_6.cos(),
        std::f64::consts::FRAC_PI_6.sin(),
    );
    // Rows are the world directions of the local axes.
    let axes = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
    let centre_of_box = [1.0, -2.0, 0.5];
    let obb = Obb::from_center_axes(centre_of_box, axes, [1.0, 2.0, 3.0]);
    let local = [0.3, 1.8, -0.5];
    let world = |p: [f64; 3]| {
        [
            centre_of_box[0] + axes[0][0] * p[0] + axes[1][0] * p[1] + axes[2][0] * p[2],
            centre_of_box[1] + axes[0][1] * p[0] + axes[1][1] * p[1] + axes[2][1] * p[2],
            centre_of_box[2] + axes[0][2] * p[0] + axes[1][2] * p[1] + axes[2][2] * p[2],
        ]
    };
    let result = obb_sphere_test(&obb, world(local), 0.4).expect("centre inside the box");
    assert!((result.penetration_depth - 0.6).abs() < 1e-12);
    assert_vec_close(result.contact_normal, axes[1], "normal");
    let surface = world([0.3, 2.0, -0.5]);
    let facing = world([0.3, 1.4, -0.5]);
    let midpoint = [
        0.5 * (surface[0] + facing[0]),
        0.5 * (surface[1] + facing[1]),
        0.5 * (surface[2] + facing[2]),
    ];
    assert_vec_close(result.contact_point, midpoint, "contact point");
}

#[test]
fn outside_contact_reports_the_overlap_midpoint() {
    let obb = Obb::axis_aligned([0.0; 3], [1.0, 1.0, 1.0]);
    let result = obb_sphere_test(&obb, [1.3, 0.0, 0.0], 0.5).expect("face contact");
    assert!((result.penetration_depth - 0.2).abs() < 1e-12);
    assert_vec_close(result.contact_normal, [1.0, 0.0, 0.0], "normal");
    // Box face point (1, 0, 0), sphere surface point (0.8, 0, 0).
    assert_vec_close(result.contact_point, [0.9, 0.0, 0.0], "contact point");
}
