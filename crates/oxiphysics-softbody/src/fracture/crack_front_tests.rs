// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the crack-front frame of [`CrackFront`].
//!
//! A front vertex carries a right-handed local frame `(propagation_dir, normal,
//! front tangent = propagation_dir x normal)`.  Propagation must advance the
//! vertex *inside* the crack plane (along `propagation_dir`, deflected by the
//! kink angle towards `normal`) and turn the frame with the kink about the
//! front tangent.  The previous implementation advanced the vertex along the
//! plane normal itself, i.e. a pure Mode-I crack left its own plane.

use super::{CrackFront, CrackFrontVertex, norm3};
use std::f64::consts::FRAC_1_SQRT_2;

const K1C: f64 = 1.0e6;
const K_BRANCH: f64 = 4.0e6;
const DA_MIN: f64 = 1.0e-4;
const DA_MAX: f64 = 1.0e-2;
const TOL: f64 = 1.0e-12;

fn front_with(vertex: CrackFrontVertex) -> CrackFront {
    let mut front = CrackFront::new(K1C, K_BRANCH, DA_MIN, DA_MAX);
    front.add_vertex(vertex);
    front
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn lin2(a: [f64; 3], ca: f64, b: [f64; 3], cb: f64) -> [f64; 3] {
    [
        ca * a[0] + cb * b[0],
        ca * a[1] + cb * b[1],
        ca * a[2] + cb * b[2],
    ]
}

fn assert_close(actual: [f64; 3], expected: [f64; 3], what: &str) {
    for d in 0..3 {
        assert!(
            (actual[d] - expected[d]).abs() < TOL,
            "{what}: component {d} is {} but {} was expected ({actual:?} vs {expected:?})",
            actual[d],
            expected[d]
        );
    }
}

fn assert_orthonormal_frame(v: &CrackFrontVertex, what: &str) {
    assert!(
        (norm3(v.normal) - 1.0).abs() < TOL,
        "{what}: |normal| = {}",
        norm3(v.normal)
    );
    assert!(
        (norm3(v.propagation_dir) - 1.0).abs() < TOL,
        "{what}: |propagation_dir| = {}",
        norm3(v.propagation_dir)
    );
    assert!(
        dot(v.normal, v.propagation_dir).abs() < TOL,
        "{what}: normal . propagation_dir = {}",
        dot(v.normal, v.propagation_dir)
    );
}

#[test]
fn constructor_builds_an_orthonormal_frame() {
    // Non-unit normal, propagation direction with a component along the normal.
    let v = CrackFrontVertex::new([0.0; 3], [0.0, 2.0, 0.0], [3.0, 1.0, 0.0], 1.0e6, 0.0);
    assert_close(v.normal, [0.0, 1.0, 0.0], "normal");
    assert_close(v.propagation_dir, [1.0, 0.0, 0.0], "propagation_dir");
    assert_orthonormal_frame(&v, "projected frame");

    // Propagation direction parallel to the normal and the zero vector both fall
    // back to some unit vector inside the crack plane.
    for degenerate in [[0.0, 3.0, 0.0], [0.0, 0.0, 0.0], [0.0, -1.0, 0.0]] {
        let v = CrackFrontVertex::new([0.0; 3], [0.0, 1.0, 0.0], degenerate, 1.0e6, 0.0);
        assert_orthonormal_frame(&v, "degenerate propagation direction");
    }

    // A zero normal keeps the historical default normal (0, 1, 0).
    let v = CrackFrontVertex::new([0.0; 3], [0.0; 3], [1.0, 0.0, 0.0], 1.0e6, 0.0);
    assert_close(v.normal, [0.0, 1.0, 0.0], "default normal");
    assert_orthonormal_frame(&v, "default normal frame");
}

#[test]
fn mode_one_front_advances_inside_the_crack_plane() {
    let n0 = [0.0, 1.0, 0.0];
    let p0 = [1.0, 0.0, 0.0];
    let origin = [0.2, 0.3, 0.4];
    let mut front = front_with(CrackFrontVertex::new(origin, n0, p0, 3.0e6, 0.0));

    let mut previous_along = 0.0;
    for step in 0..20 {
        assert_eq!(front.propagate(), 1, "step {step}: the front must advance");
        let fv = front.vertices[0];
        let offset = sub(fv.position, origin);
        assert!(
            dot(offset, n0).abs() < 1.0e-14,
            "step {step}: the front left the crack plane, offset along the normal = {}",
            dot(offset, n0)
        );
        assert!(
            dot(offset, cross(p0, n0)).abs() < 1.0e-14,
            "step {step}: the front moved along its own tangent"
        );
        let along = dot(offset, p0);
        assert!(
            along - previous_along >= DA_MIN - 1.0e-15 && along - previous_along <= DA_MAX + 1e-15,
            "step {step}: increment {} outside [da_min, da_max]",
            along - previous_along
        );
        previous_along = along;
        assert_close(fv.normal, n0, "plane normal of a straight crack");
        assert_close(
            fv.propagation_dir,
            p0,
            "propagation direction of a straight crack",
        );
    }
}

#[test]
fn mode_two_front_kinks_and_the_frame_turns_with_it() {
    let n0 = [0.0, 1.0, 0.0];
    let p0 = [1.0, 0.0, 0.0];
    let origin = [1.0, -2.0, 0.5];
    let vertex = CrackFrontVertex::new(origin, n0, p0, 3.0e6, 1.0e6);
    let theta = vertex.max_hoop_stress_angle();
    assert!(
        theta < -0.2 && theta > -1.2,
        "K_II > 0 must kink the crack to the negative side, theta = {theta}"
    );
    let mut front = front_with(vertex);
    assert_eq!(front.propagate(), 1);

    let fv = front.vertices[0];
    let step = sub(fv.position, origin);
    let da = norm3(step);
    assert!(
        (DA_MIN - TOL..=DA_MAX + TOL).contains(&da),
        "increment {da} outside [da_min, da_max]"
    );

    // The vertex advances along cos(theta) p0 + sin(theta) n0 ...
    let expected_dir = lin2(p0, theta.cos(), n0, theta.sin());
    assert_close(
        [step[0] / da, step[1] / da, step[2] / da],
        expected_dir,
        "advance direction",
    );
    assert!(dot(step, n0) < 0.0, "the kink goes towards -n for K_II > 0");

    // ... and the local frame rotated by theta about the (unchanged) front tangent.
    assert_close(
        fv.propagation_dir,
        expected_dir,
        "new propagation direction",
    );
    assert_close(
        fv.normal,
        lin2(n0, theta.cos(), p0, -theta.sin()),
        "new plane normal",
    );
    assert_orthonormal_frame(&fv, "kinked frame");
    assert_close(
        cross(fv.propagation_dir, fv.normal),
        cross(p0, n0),
        "front tangent",
    );
}

#[test]
fn successive_kinks_accumulate_the_rotation() {
    let n0 = [0.0, 1.0, 0.0];
    let p0 = [1.0, 0.0, 0.0];
    let origin = [0.0; 3];
    let vertex = CrackFrontVertex::new(origin, n0, p0, 3.0e6, 5.0e5);
    let theta = vertex.max_hoop_stress_angle();
    let mut front = front_with(vertex);

    assert_eq!(front.propagate(), 1);
    let after_first = front.vertices[0].position;
    assert_eq!(front.propagate(), 1);
    let fv = front.vertices[0];
    let second = sub(fv.position, after_first);
    let da = norm3(second);

    // The second kink acts in the already turned frame: total rotation 2 theta.
    let expected = lin2(p0, (2.0 * theta).cos(), n0, (2.0 * theta).sin());
    assert_close(
        [second[0] / da, second[1] / da, second[2] / da],
        expected,
        "second advance direction",
    );
    assert_close(
        fv.propagation_dir,
        expected,
        "propagation direction after two kinks",
    );
    assert_close(
        fv.normal,
        lin2(n0, (2.0 * theta).cos(), p0, -(2.0 * theta).sin()),
        "normal after two kinks",
    );
    assert_orthonormal_frame(&fv, "twice kinked frame");
}

#[test]
fn kink_rotates_about_the_front_tangent_in_a_tilted_plane() {
    // Orthonormal triple with a generic orientation (rows of a Pythagorean rotation).
    let n0 = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let p0 = [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0];
    let tangent = cross(p0, n0);
    let origin = [4.0, 5.0, 6.0];
    let vertex = CrackFrontVertex::new(origin, n0, p0, 3.0e6, 1.2e6);
    assert_orthonormal_frame(&vertex, "tilted frame");
    let theta = vertex.max_hoop_stress_angle();
    let mut front = front_with(vertex);
    assert_eq!(front.propagate(), 1);

    let fv = front.vertices[0];
    let step = sub(fv.position, origin);
    assert!(
        dot(step, tangent).abs() < 1.0e-14,
        "the vertex must not move along the front line, got {}",
        dot(step, tangent)
    );
    let da = norm3(step);
    assert_close(
        [step[0] / da, step[1] / da, step[2] / da],
        lin2(p0, theta.cos(), n0, theta.sin()),
        "tilted advance direction",
    );
    assert_close(
        fv.propagation_dir,
        [step[0] / da, step[1] / da, step[2] / da],
        "frame follows advance",
    );
    assert_orthonormal_frame(&fv, "tilted kinked frame");
    assert_close(
        cross(fv.propagation_dir, fv.normal),
        tangent,
        "front tangent in a tilted plane",
    );
}

#[test]
fn branches_straddle_the_propagation_direction() {
    let n0 = [0.0, 1.0, 0.0];
    let p0 = [1.0, 0.0, 0.0];
    let origin = [0.5, 0.25, -0.75];
    let front = front_with(CrackFrontVertex::new(origin, n0, p0, 6.0e6, 0.0));
    let branches = front.check_branching();
    assert_eq!(branches.len(), 2, "K_eq above k_branch gives two branches");

    for (i, b) in branches.iter().enumerate() {
        assert_eq!(b.position, origin, "branch {i} starts at the parent vertex");
        assert_orthonormal_frame(b, "branch frame");
        assert!(
            (dot(b.propagation_dir, p0) - FRAC_1_SQRT_2).abs() < TOL,
            "branch {i} must deviate 45 degrees from the propagation direction, cos = {}",
            dot(b.propagation_dir, p0)
        );
        assert_close(
            cross(b.propagation_dir, b.normal),
            cross(p0, n0),
            "branch front tangent",
        );
        assert!(!b.arrested);
    }
    assert!(
        dot(branches[0].propagation_dir, n0) > 0.0 && dot(branches[1].propagation_dir, n0) < 0.0,
        "the two branches lie on opposite sides of the crack plane"
    );
    assert!(
        dot(branches[0].propagation_dir, branches[1].propagation_dir).abs() < TOL,
        "the branches are 90 degrees apart"
    );
    // Branching never changes the parent.
    assert_eq!(front.vertices[0].propagation_dir, p0);
    assert_eq!(front.vertices[0].normal, n0);
}
