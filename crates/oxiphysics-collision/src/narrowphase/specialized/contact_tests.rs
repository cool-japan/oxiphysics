// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the specialised pairwise contact routines.
//!
//! Every [`Contact`] follows the crate convention: the normal points from body
//! B towards body A, `point_a` lies on A, `point_b` lies on B and
//! `point_b - point_a == normal * depth`.

use super::*;
use crate::narrowphase::epa::Epa;
use crate::narrowphase::gjk::{Gjk, GjkResult};
use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Capsule, Cylinder, Sphere, Torus};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

fn assert_close(actual: Vec3, expected: Vec3, tol: f64, what: &str) {
    assert!(
        (actual - expected).norm() <= tol,
        "{what}: {actual:?} differs from {expected:?}"
    );
}

/// The witness points must be one penetration depth apart along the normal.
fn assert_witness_gap(c: &Contact, tol: f64) {
    assert_close(
        c.point_b - c.point_a,
        c.normal * c.depth,
        tol,
        "point_b - point_a",
    );
}

/// The witness points of an EPA contact are a barycentric mix of support points
/// of one face of the Minkowski difference, so they agree with the normal and
/// the depth along the normal (their tangential offset is not meaningful).
fn assert_gap_along_normal(c: &Contact, tol: f64) {
    let gap = (c.point_b - c.point_a).dot(&c.normal);
    assert!(
        (gap - c.depth).abs() <= tol,
        "gap along the normal {gap} differs from the depth {}",
        c.depth
    );
}

fn cube_vertices(half: f64) -> Vec<Vec3> {
    let mut vertices = Vec::new();
    for sx in [-half, half] {
        for sy in [-half, half] {
            for sz in [-half, half] {
                vertices.push(Vec3::new(sx, sy, sz));
            }
        }
    }
    vertices
}

/// Deterministic LCG in `[0, 1)`.
fn next_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

// ── sphere_box ──────────────────────────────────────────────────────────────

#[test]
fn sphere_box_centre_inside_clears_the_face_plus_the_radius() {
    let sphere = Sphere::new(0.1);
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let origin = Transform::from_position(Vec3::zeros());
    let c = sphere_box(&sphere, &origin, &cube, &origin).expect("centre inside the box");
    // The centre is 1.0 from every face; the whole sphere has to clear one.
    assert!((c.depth - 1.1).abs() < 1e-12, "depth {}", c.depth);
    assert!((c.normal.norm() - 1.0).abs() < 1e-12);
    assert!((c.normal.amax() - 1.0).abs() < 1e-12, "axis-aligned normal");
    // The witness on A is on the sphere surface, the witness on B on the face.
    assert!(
        (c.point_a.norm() - 0.1).abs() < 1e-12,
        "point_a {:?}",
        c.point_a
    );
    assert!(
        (c.point_b.dot(&c.normal) - 1.0).abs() < 1e-12,
        "{:?}",
        c.point_b
    );
    assert_witness_gap(&c, 1e-12);
}

#[test]
fn sphere_box_off_centre_inside_uses_the_nearest_face_of_a_rotated_box() {
    let sphere = Sphere::new(0.4);
    let cuboid = BoxShape::new(Vec3::new(1.0, 2.0, 3.0));
    let box_t = Transform::from_axis_angle(
        Vec3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_6,
    );
    // Local (0.3, 1.8, -0.5): 0.2 below the +y face, 0.7 from +x, 2.5 from -z.
    let sphere_t = Transform::from_position(box_t.transform_point(&Vec3::new(0.3, 1.8, -0.5)));
    let c = sphere_box(&sphere, &sphere_t, &cuboid, &box_t).expect("centre inside the box");
    assert!((c.depth - 0.6).abs() < 1e-12, "depth {}", c.depth);
    let expected_normal = box_t.transform_vector(&Vec3::new(0.0, 1.0, 0.0));
    assert_close(c.normal, expected_normal, 1e-12, "normal");
    let expected_b = box_t.transform_point(&Vec3::new(0.3, 2.0, -0.5));
    assert_close(c.point_b, expected_b, 1e-12, "point_b");
    assert_close(
        c.point_a,
        sphere_t.position - c.normal * 0.4,
        1e-12,
        "point_a",
    );
    assert_witness_gap(&c, 1e-12);
}

#[test]
fn sphere_box_outside_contacts_keep_both_witnesses_on_the_surfaces() {
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let box_t = Transform::from_position(Vec3::zeros());
    // Face contact.
    let sphere = Sphere::new(0.5);
    let t = Transform::from_position(Vec3::new(1.3, 0.2, -0.1));
    let c = sphere_box(&sphere, &t, &cube, &box_t).expect("face contact");
    assert!((c.depth - 0.2).abs() < 1e-12, "depth {}", c.depth);
    assert_close(c.normal, Vec3::new(1.0, 0.0, 0.0), 1e-12, "normal");
    assert_close(c.point_b, Vec3::new(1.0, 0.2, -0.1), 1e-12, "point_b");
    assert_close(c.point_a, Vec3::new(0.8, 0.2, -0.1), 1e-12, "point_a");
    assert_witness_gap(&c, 1e-12);
    // Corner contact: the normal points along the diagonal.
    let big = Sphere::new(1.0);
    let t = Transform::from_position(Vec3::new(1.5, 1.5, 1.5));
    let c = sphere_box(&big, &t, &cube, &box_t).expect("corner contact");
    let diagonal = 0.75_f64.sqrt();
    assert!(
        (c.depth - (1.0 - diagonal)).abs() < 1e-12,
        "depth {}",
        c.depth
    );
    assert_close(
        c.normal,
        Vec3::new(1.0, 1.0, 1.0) / 3.0_f64.sqrt(),
        1e-12,
        "normal",
    );
    assert_close(c.point_b, Vec3::new(1.0, 1.0, 1.0), 1e-12, "point_b");
    assert!(((c.point_a - t.position).norm() - 1.0).abs() < 1e-12);
    assert_witness_gap(&c, 1e-12);
}

// ── box_box_sat ─────────────────────────────────────────────────────────────

#[test]
fn box_box_sat_normal_points_from_b_towards_a() {
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
    let c = box_box_sat(&cube, &ta, &cube, &tb).expect("overlap");
    // B sits on the +x side of A, so the normal from B towards A is -x.
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-12, "normal");
    assert!((c.depth - 0.5).abs() < 1e-12, "depth {}", c.depth);
    assert_close(c.point_a, Vec3::new(1.0, 0.0, 0.0), 1e-12, "A's +x face");
    assert_close(c.point_b, Vec3::new(0.5, 0.0, 0.0), 1e-12, "B's -x face");
    assert_witness_gap(&c, 1e-12);
    // Swapping the arguments mirrors the contact.
    let m = box_box_sat(&cube, &tb, &cube, &ta).expect("overlap");
    assert_close(m.normal, -c.normal, 1e-12, "mirrored normal");
    assert_close(m.point_a, c.point_b, 1e-12, "mirrored point_a");
    assert_close(m.point_b, c.point_a, 1e-12, "mirrored point_b");
    assert!((m.depth - c.depth).abs() < 1e-12);
}

#[test]
fn box_box_sat_edge_edge_depth_is_a_true_distance() {
    // Two unit cubes crossing edge to edge: A rotated 45 degrees about z (its
    // top edge runs along z), B rotated 45 degrees about x (its bottom edge
    // runs along x). The edge cross product is not a unit vector, so an
    // unnormalised depth would be wrong.
    let cube = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let overlap = 0.02;
    let corner = 0.5 * 2.0_f64.sqrt();
    let ta = Transform::from_axis_angle(Vec3::zeros(), Vec3::new(0.0, 0.0, 1.0), FRAC_PI_4);
    let tb = Transform::from_axis_angle(
        Vec3::new(0.0, 2.0 * corner - overlap, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        FRAC_PI_4,
    );
    let c = box_box_sat(&cube, &ta, &cube, &tb).expect("the edges cross");
    assert!((c.depth - overlap).abs() < 1e-9, "depth {}", c.depth);
    assert_close(c.normal, Vec3::new(0.0, -1.0, 0.0), 1e-9, "normal");
    assert_close(
        c.point_a,
        Vec3::new(0.0, corner, 0.0),
        1e-9,
        "top edge of A",
    );
    assert_close(
        c.point_b,
        Vec3::new(0.0, corner - overlap, 0.0),
        1e-9,
        "bottom edge of B",
    );
    assert_witness_gap(&c, 1e-9);
    // Independent check with GJK + EPA on the same boxes.
    let simplex = match Gjk::query(&cube, &ta, &cube, &tb) {
        GjkResult::Intersecting(simplex) => simplex,
        GjkResult::Separated { .. } => panic!("the cubes overlap by {overlap}"),
    };
    let epa = Epa::penetration_depth(&cube, &ta, &cube, &tb, &simplex).expect("EPA contact");
    assert!(
        (epa.depth - c.depth).abs() < 1e-6,
        "EPA depth {}",
        epa.depth
    );
    assert!(
        epa.normal.dot(&c.normal) > 1.0 - 1e-9,
        "EPA normal {:?}",
        epa.normal
    );
}

#[test]
fn box_box_sat_resting_box_reports_the_centre_of_the_contact_patch() {
    let upper = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let lower = BoxShape::new(Vec3::new(1.0, 0.5, 1.0));
    let tb = Transform::from_position(Vec3::zeros());
    // Centred above: the patch is A's whole bottom face.
    let ta = Transform::from_position(Vec3::new(0.0, 0.99, 0.0));
    let c = box_box_sat(&upper, &ta, &lower, &tb).expect("resting overlap");
    assert!((c.depth - 0.01).abs() < 1e-9, "depth {}", c.depth);
    assert_close(c.normal, Vec3::new(0.0, 1.0, 0.0), 1e-12, "normal");
    assert_close(
        c.point_a,
        Vec3::new(0.0, 0.49, 0.0),
        1e-9,
        "A's bottom face",
    );
    assert_close(c.point_b, Vec3::new(0.0, 0.5, 0.0), 1e-9, "B's top face");
    assert_witness_gap(&c, 1e-9);
    // Overhanging the edge of B: the patch is the part of A's face over B,
    // x in [0.1, 1.0], z in [-0.5, 0.5], whose centre is at x = 0.55.
    let ta = Transform::from_position(Vec3::new(0.6, 0.99, 0.0));
    let c = box_box_sat(&upper, &ta, &lower, &tb).expect("overhanging overlap");
    assert_close(
        c.point_a,
        Vec3::new(0.55, 0.49, 0.0),
        1e-9,
        "patch centre on A",
    );
    assert_close(
        c.point_b,
        Vec3::new(0.55, 0.5, 0.0),
        1e-9,
        "patch centre on B",
    );
    assert_witness_gap(&c, 1e-9);
}

#[test]
fn box_box_sat_separated_boxes_have_no_contact() {
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(2.1, 0.0, 0.0));
    assert!(box_box_sat(&cube, &ta, &cube, &tb).is_none());
}

// ── sphere_sphere / capsule_capsule / cylinder_cylinder ─────────────────────

/// The solver pushes body A along `+normal` only while the bodies approach
/// along the normal, i.e. while `(v_a - v_b) . normal < 0`.
fn approach_speed_along(c: &Contact, v_a: Vec3, v_b: Vec3) -> f64 {
    -(v_a - v_b).dot(&c.normal)
}

#[test]
fn sphere_sphere_normal_points_from_b_towards_a() {
    let big = Sphere::new(1.0);
    let small = Sphere::new(0.5);
    let ta = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
    let tb = Transform::from_position(Vec3::new(1.2, 0.0, 0.0));
    let c = sphere_sphere(&big, &ta, &small, &tb).expect("overlap");
    // B sits on the +x side of A, so the normal from B towards A is -x.
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-12, "normal");
    assert!((c.depth - 0.3).abs() < 1e-12, "depth {}", c.depth);
    assert_close(c.point_a, Vec3::new(1.0, 0.0, 0.0), 1e-12, "A facing B");
    assert_close(c.point_b, Vec3::new(0.7, 0.0, 0.0), 1e-12, "B facing A");
    assert_witness_gap(&c, 1e-12);
    // Approaching bodies: A moves +x, B moves -x.
    let closing = approach_speed_along(&c, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
    assert!(closing > 1.99, "closing speed {closing}");

    // Exchanging the arguments mirrors the contact.
    let m = sphere_sphere(&small, &tb, &big, &ta).expect("overlap");
    assert_close(m.normal, -c.normal, 1e-12, "mirrored normal");
    assert_close(m.point_a, c.point_b, 1e-12, "mirrored point_a");
    assert_close(m.point_b, c.point_a, 1e-12, "mirrored point_b");
}

#[test]
fn sphere_sphere_oblique_and_coincident_cases() {
    let a = Sphere::new(0.8);
    let b = Sphere::new(0.6);
    let ta = Transform::from_position(Vec3::new(1.0, -2.0, 0.5));
    let tb = Transform::from_position(Vec3::new(1.6, -1.2, 1.3));
    let c = sphere_sphere(&a, &ta, &b, &tb).expect("overlap");
    let offset = ta.position - tb.position;
    assert_close(c.normal, offset / offset.norm(), 1e-12, "normal");
    assert!((c.depth - (1.4 - offset.norm())).abs() < 1e-12);
    assert!(((c.point_a - ta.position).norm() - 0.8).abs() < 1e-12);
    assert!(((c.point_b - tb.position).norm() - 0.6).abs() < 1e-12);
    assert_witness_gap(&c, 1e-12);
    // Coincident centres: no direction, +Y by convention.
    let same = Transform::from_position(Vec3::new(3.0, 3.0, 3.0));
    let c = sphere_sphere(&a, &same, &b, &same).expect("overlap");
    assert_close(c.normal, Vec3::new(0.0, 1.0, 0.0), 1e-12, "normal");
    assert!((c.depth - 1.4).abs() < 1e-12);
}

#[test]
fn capsule_capsule_normals_point_from_b_towards_a() {
    // Two parallel upright capsules side by side.
    let a = Capsule::new(0.6, 1.0);
    let b = Capsule::new(0.6, 1.0);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
    let c = capsule_capsule(&a, &ta, &b, &tb).expect("overlap");
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-12, "normal");
    assert!((c.depth - 0.2).abs() < 1e-12, "depth {}", c.depth);
    assert!((c.point_a.x - 0.6).abs() < 1e-12 && (c.point_b.x - 0.4).abs() < 1e-12);
    assert_witness_gap(&c, 1e-12);
    let closing = approach_speed_along(&c, Vec3::new(2.0, 0.0, 0.0), Vec3::zeros());
    assert!(closing > 1.99, "closing speed {closing}");

    // A crosses B: B lies on a horizontal axis 0.6 above (+z) A's vertical axis.
    let a = Capsule::new(0.5, 2.0);
    let b = Capsule::new(0.4, 2.0);
    let tb = Transform::from_axis_angle(
        Vec3::new(0.0, 0.0, 0.6),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    let c = capsule_capsule(&a, &ta, &b, &tb).expect("crossing capsules overlap");
    assert_close(c.normal, Vec3::new(0.0, 0.0, -1.0), 1e-9, "normal");
    assert!((c.depth - 0.3).abs() < 1e-9, "depth {}", c.depth);
    assert_close(c.point_a, Vec3::new(0.0, 0.0, 0.5), 1e-9, "A facing B");
    assert_close(c.point_b, Vec3::new(0.0, 0.0, 0.2), 1e-9, "B facing A");
    assert_witness_gap(&c, 1e-9);
    // Swapped roles mirror the contact.
    let m = capsule_capsule(&b, &tb, &a, &ta).expect("overlap");
    assert_close(m.normal, -c.normal, 1e-9, "mirrored normal");
    assert_close(m.point_a, c.point_b, 1e-9, "mirrored point_a");
}

#[test]
fn cylinder_cylinder_normals_point_from_b_towards_a() {
    let a = Cylinder::new(0.5, 2.0);
    let b = Cylinder::new(0.5, 2.0);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_axis_angle(
        Vec3::new(0.0, 0.0, 0.4),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    let c = cylinder_cylinder(&a, &ta, &b, &tb).expect("crossing cylinders overlap");
    assert_close(c.normal, Vec3::new(0.0, 0.0, -1.0), 1e-9, "normal");
    assert!((c.depth - 0.6).abs() < 1e-9, "depth {}", c.depth);
    assert_witness_gap(&c, 1e-9);
    assert_close(c.point_a, Vec3::new(0.0, 0.0, 0.5), 1e-9, "A facing B");

    // The precise test agrees, in both of its branches.
    let precise = cylinder_cylinder_precise(&a, &ta, &b, &tb);
    let p = precise.contact.expect("crossing cylinders overlap");
    assert_close(p.normal, c.normal, 1e-9, "precise normal");
    assert!((p.depth - c.depth).abs() < 1e-9);
    assert_close(p.point_a, c.point_a, 1e-9, "precise point_a");
    assert_close(p.point_b, c.point_b, 1e-9, "precise point_b");

    let parallel_b = Transform::from_position(Vec3::new(0.8, 0.0, 0.0));
    let parallel = cylinder_cylinder_precise(&a, &ta, &b, &parallel_b);
    assert_eq!(parallel.mode, CylinderContactMode::Parallel);
    let p = parallel.contact.expect("side by side");
    assert_close(
        p.normal,
        Vec3::new(-1.0, 0.0, 0.0),
        1e-12,
        "parallel normal",
    );
    assert!((p.depth - 0.2).abs() < 1e-12, "depth {}", p.depth);
    assert_close(p.point_a, Vec3::new(0.5, 0.0, 0.0), 1e-12, "A facing B");
    assert_close(p.point_b, Vec3::new(0.3, 0.0, 0.0), 1e-12, "B facing A");
    assert_witness_gap(&p, 1e-12);
    let plain = cylinder_cylinder(&a, &ta, &b, &parallel_b).expect("side by side");
    assert_close(plain.normal, p.normal, 1e-12, "plain parallel normal");
    // The same pair with the roles exchanged.
    let swapped = cylinder_cylinder_precise(&b, &parallel_b, &a, &ta)
        .contact
        .expect("side by side");
    assert_close(
        swapped.normal,
        Vec3::new(1.0, 0.0, 0.0),
        1e-12,
        "swapped normal",
    );
}

// ── torus_plane ─────────────────────────────────────────────────────────────

/// Lowest value of `n . x` over a dense sampling of the torus surface.
fn sampled_lowest(torus: &Torus, t: &Transform, n: &Vec3) -> f64 {
    let axis = t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let e1 = t.rotation * Vec3::new(1.0, 0.0, 0.0);
    let e2 = axis.cross(&e1);
    let mut lowest = f64::INFINITY;
    for iu in 0..720 {
        let u = std::f64::consts::TAU * f64::from(iu) / 720.0;
        let ring = e1 * u.cos() + e2 * u.sin();
        for iv in 0..360 {
            let v = std::f64::consts::TAU * f64::from(iv) / 360.0;
            let p = t.position
                + ring * (torus.major_radius + torus.minor_radius * v.cos())
                + axis * (torus.minor_radius * v.sin());
            lowest = lowest.min(n.dot(&p));
        }
    }
    lowest
}

/// Distance of `p` from the core ring of the torus (must equal the minor radius
/// for a point on the surface).
fn ring_distance(torus: &Torus, t: &Transform, p: &Vec3) -> f64 {
    let axis = t.rotation * Vec3::new(0.0, 1.0, 0.0);
    let rel = p - t.position;
    let h = rel.dot(&axis);
    let rho = (rel - axis * h).norm();
    (rho - torus.major_radius).hypot(h)
}

#[test]
fn torus_plane_upright_ring_tests_its_lowest_point() {
    // Axis along +x: the ring stands in the y-z plane, so its lowest point is
    // R + r below the centre.
    let torus = Torus::new(1.0, 0.25);
    let up = Vec3::new(0.0, 1.0, 0.0);
    let pose = |y: f64| {
        Transform::from_axis_angle(Vec3::new(0.0, y, 0.0), Vec3::new(0.0, 0.0, 1.0), -FRAC_PI_2)
    };
    assert_close(
        pose(1.2).rotation * up,
        Vec3::new(1.0, 0.0, 0.0),
        1e-12,
        "torus axis",
    );
    let c = torus_plane(&torus, &pose(1.2), &up, 0.0).expect("lowest point is at y = -0.05");
    assert!((c.depth - 0.05).abs() < 1e-9, "depth {}", c.depth);
    assert_close(c.normal, up, 1e-12, "normal");
    assert!((c.point_a.y + 0.05).abs() < 1e-9, "point_a {:?}", c.point_a);
    assert!(c.point_b.y.abs() < 1e-9, "point_b {:?}", c.point_b);
    assert!((ring_distance(&torus, &pose(1.2), &c.point_a) - 0.25).abs() < 1e-9);
    assert_witness_gap(&c, 1e-9);
    // 0.1 higher the lowest point is above the plane.
    assert!(torus_plane(&torus, &pose(1.3), &up, 0.0).is_none());
}

#[test]
fn torus_plane_tilted_torus_matches_the_analytic_depth_and_a_surface_sampling() {
    let torus = Torus::new(1.0, 0.25);
    let up = Vec3::new(0.0, 1.0, 0.0);
    // Axis tilted by 45 degrees: the ring's lowest point is R * sin(45) + r below
    // the centre.
    let tilted = Transform::from_axis_angle(
        Vec3::new(0.0, 0.8, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        -FRAC_PI_4,
    );
    let c = torus_plane(&torus, &tilted, &up, 0.0).expect("lowest point is below the plane");
    let expected = FRAC_PI_4.sin() + 0.25 - 0.8;
    assert!((c.depth - expected).abs() < 1e-9, "depth {}", c.depth);
    assert!((ring_distance(&torus, &tilted, &c.point_a) - 0.25).abs() < 1e-9);
    assert_witness_gap(&c, 1e-9);
    // Brute-force cross-check over several orientations and plane normals: the
    // depth is exactly the lowest sampled point (to within the sampling step).
    let cases = [
        (0.0, Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0)),
        (0.3, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
        (1.1, Vec3::new(1.0, 1.0, 0.0), Vec3::new(0.0, 1.0, 1.0)),
        (2.4, Vec3::new(0.3, -0.8, 0.5), Vec3::new(1.0, 2.0, -1.0)),
        (
            FRAC_PI_2,
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
        ),
    ];
    for (angle, axis, normal) in cases {
        let pose = Transform::from_axis_angle(Vec3::new(0.1, 0.4, -0.2), axis, angle);
        let n_hat = normal / normal.norm();
        let d = 0.35;
        let contact = torus_plane(&torus, &pose, &n_hat, d)
            .unwrap_or_else(|| panic!("angle {angle}: the torus reaches the plane"));
        let lowest = sampled_lowest(&torus, &pose, &n_hat);
        let analytic_depth = d - lowest;
        assert!(
            (contact.depth - analytic_depth).abs() < 2e-4,
            "angle {angle}: depth {} vs sampled {analytic_depth}",
            contact.depth
        );
        // The sampled surface never goes below the reported lowest point.
        assert!(lowest >= d - contact.depth - 1e-9, "angle {angle}");
        assert!((ring_distance(&torus, &pose, &contact.point_a) - 0.25).abs() < 1e-9);
        assert!(
            (n_hat.dot(&contact.point_b) - d).abs() < 1e-9,
            "point_b is on the plane"
        );
    }
}

#[test]
fn torus_plane_axis_parallel_to_the_normal_depends_on_the_tube_only() {
    // Plane x = 0 with the torus axis along +x: the ring is parallel to the
    // plane, so the major radius does not matter.
    let n = Vec3::new(1.0, 0.0, 0.0);
    for major in [1.0, 3.0] {
        let torus = Torus::new(major, 0.25);
        let pose = Transform::from_axis_angle(
            Vec3::new(0.2, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            -FRAC_PI_2,
        );
        let c = torus_plane(&torus, &pose, &n, 0.0).expect("tube dips into the plane");
        assert!(
            (c.depth - 0.05).abs() < 1e-9,
            "R = {major}: depth {}",
            c.depth
        );
        assert!(
            (c.point_a.x + 0.05).abs() < 1e-9,
            "R = {major}: {:?}",
            c.point_a
        );
        assert!(c.point_b.x.abs() < 1e-9);
        assert!((ring_distance(&torus, &pose, &c.point_a) - 0.25).abs() < 1e-9);
    }
}

#[test]
fn torus_plane_flat_torus_depth_is_the_tube_radius_minus_the_height() {
    let torus = Torus::new(1.0, 0.3);
    let up = Vec3::new(0.0, 1.0, 0.0);
    let flat = Transform::from_position(Vec3::new(0.0, 0.2, 0.0));
    let c = torus_plane(&torus, &flat, &up, 0.0).expect("tube dips into the plane");
    assert!((c.depth - 0.1).abs() < 1e-12, "depth {}", c.depth);
    assert!((c.point_a.y + 0.1).abs() < 1e-12);
    assert!(
        torus_plane(
            &torus,
            &Transform::from_position(Vec3::new(0.0, 0.31, 0.0)),
            &up,
            0.0
        )
        .is_none()
    );
}

// ── convex_hull_vs_convex_hull ──────────────────────────────────────────────

#[test]
fn convex_hull_contact_reports_the_penetration_axis_and_depth() {
    let cube = cube_vertices(0.5);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(0.9, 0.1, 0.0));
    let c = convex_hull_vs_convex_hull(&cube, &ta, &cube, &tb).expect("cubes overlap by 0.1");
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-6, "normal");
    assert!((c.depth - 0.1).abs() < 1e-6, "depth {}", c.depth);
    assert!((c.point_a.x - 0.5).abs() < 1e-6, "point_a {:?}", c.point_a);
    assert!((c.point_b.x - 0.4).abs() < 1e-6, "point_b {:?}", c.point_b);
    assert_gap_along_normal(&c, 1e-6);
    // The gjk_epa entry point is the same query.
    let again = convex_hull_contact_gjk_epa(&cube, &ta, &cube, &tb).expect("same query");
    assert_close(again.normal, c.normal, 1e-12, "normal");
    assert!((again.depth - c.depth).abs() < 1e-12);
}

#[test]
fn convex_hull_contact_uses_the_rotation_of_both_hulls() {
    // B is a cube turned 45 degrees about z: its corner reaches 0.5 * sqrt(2)
    // along x, so the cubes overlap by 0.5 + 0.7071 - 1.0 along x.
    let cube = cube_vertices(0.5);
    let expected = 0.5 + 0.5 * 2.0_f64.sqrt() - 1.0;
    let still = Transform::from_position(Vec3::zeros());
    let turned = Transform::from_axis_angle(
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_4,
    );
    let c = convex_hull_vs_convex_hull(&cube, &still, &cube, &turned).expect("overlap");
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-6, "normal");
    assert!((c.depth - expected).abs() < 1e-6, "depth {}", c.depth);
    assert_gap_along_normal(&c, 1e-6);
    // The same pair with the roles of the two hulls exchanged.
    let swapped = convex_hull_vs_convex_hull(&cube, &turned, &cube, &still).expect("overlap");
    assert_close(
        swapped.normal,
        Vec3::new(1.0, 0.0, 0.0),
        1e-6,
        "swapped normal",
    );
    assert!((swapped.depth - expected).abs() < 1e-6);
    // Turning the *first* hull instead of the second gives the same overlap.
    let a_turned = Transform::from_axis_angle(Vec3::zeros(), Vec3::new(0.0, 0.0, 1.0), FRAC_PI_4);
    let b_still = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
    let c = convex_hull_vs_convex_hull(&cube, &a_turned, &cube, &b_still).expect("overlap");
    assert_close(c.normal, Vec3::new(-1.0, 0.0, 0.0), 1e-6, "normal");
    assert!((c.depth - expected).abs() < 1e-6, "depth {}", c.depth);
}

#[test]
fn convex_hull_contact_matches_box_box_sat_for_random_poses() {
    let cube = cube_vertices(0.5);
    let shape = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let mut state = 0x1234_5678_9ABC_DEF0_u64;
    let random_pose = |state: &mut u64, centre: Vec3| {
        let axis = Vec3::new(
            next_unit(state) - 0.5,
            next_unit(state) - 0.5,
            next_unit(state) - 0.5,
        ) + Vec3::new(0.0, 0.0, 1e-3);
        Transform::from_axis_angle(centre, axis, next_unit(state) * std::f64::consts::PI)
    };
    for case in 0..40 {
        let ta = random_pose(&mut state, Vec3::zeros());
        // Centres closer than 1.0 always overlap (the inscribed spheres do).
        let direction = Vec3::new(
            next_unit(&mut state) - 0.5,
            next_unit(&mut state) - 0.5,
            next_unit(&mut state) - 0.5,
        );
        let offset = direction / direction.norm() * (0.3 + 0.6 * next_unit(&mut state));
        let tb = random_pose(&mut state, offset);
        let reference = box_box_sat(&shape, &ta, &shape, &tb)
            .unwrap_or_else(|| panic!("case {case}: the boxes overlap"));
        let hull = convex_hull_vs_convex_hull(&cube, &ta, &cube, &tb)
            .unwrap_or_else(|| panic!("case {case}: the hulls overlap"));
        assert!(
            (hull.depth - reference.depth).abs() < 1e-4,
            "case {case}: hull depth {} vs SAT depth {}",
            hull.depth,
            reference.depth
        );
        assert!(
            hull.normal.dot(&reference.normal) > 1.0 - 1e-6,
            "case {case}: hull normal {:?} vs SAT normal {:?}",
            hull.normal,
            reference.normal
        );
    }
}

#[test]
fn convex_hull_contact_rejects_separated_and_empty_hulls() {
    let cube = cube_vertices(0.5);
    let ta = Transform::from_position(Vec3::zeros());
    let apart = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
    assert!(convex_hull_vs_convex_hull(&cube, &ta, &cube, &apart).is_none());
    assert!(convex_hull_vs_convex_hull(&[], &ta, &cube, &ta).is_none());
    assert!(convex_hull_vs_convex_hull(&cube, &ta, &[], &ta).is_none());
}

#[test]
fn convex_hull_contact_survives_degenerate_vertex_sets() {
    let cube = cube_vertices(0.5);
    let ta = Transform::from_position(Vec3::zeros());
    // A single point and a flat quad inside the cube: no panic, and whatever is
    // reported is a sane contact.
    let point = [Vec3::new(0.1, 0.0, 0.0)];
    let quad = [
        Vec3::new(-0.2, 0.0, -0.2),
        Vec3::new(0.2, 0.0, -0.2),
        Vec3::new(0.2, 0.0, 0.2),
        Vec3::new(-0.2, 0.0, 0.2),
    ];
    for degenerate in [&point[..], &quad[..]] {
        if let Some(c) = convex_hull_vs_convex_hull(&cube, &ta, degenerate, &ta) {
            assert!(c.depth.is_finite() && c.depth >= 0.0, "depth {}", c.depth);
            assert!(
                (c.normal.norm() - 1.0).abs() < 1e-6,
                "normal {:?}",
                c.normal
            );
        }
    }
}
