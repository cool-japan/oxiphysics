// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the moving-anchor spring of [`SteeredMd`].
//!
//! `SteeredMd::apply_force` used to rebuild the anchor from the *current*
//! position of the pulled atom, so the spring force reduced to `k v t` whatever
//! the atom did and there was no restoring force.  The anchor is now a point of
//! its own: it starts at the atom's initial position (or at an explicit origin)
//! and moves independently of the atom.

use super::*;

/// Force that the probe adds to a single pulled atom at `atom_pos` at `time`.
fn spring_force(smd: &mut SteeredMd, atom_pos: [f64; 3], time: f64) -> [f64; 3] {
    let positions = [atom_pos];
    let mut forces = [[0.0_f64; 3]];
    smd.apply_force(&positions, &mut forces, time, 1.0e-3);
    forces[0]
}

/// Component-wise comparison with an absolute tolerance.
fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() <= tolerance,
            "axis {axis}: got {actual:?}, expected {expected:?}"
        );
    }
}

#[test]
fn explicit_anchor_pulls_with_k_v_t_along_the_direction() {
    // Atom fixed at the origin, anchor starting there at t = 0: after t = 2 the
    // anchor is 2 nm away and the spring pulls with k * 2 = 20.
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0).with_anchor([0.0; 3], 0.0);
    let force = spring_force(&mut smd, [0.0; 3], 2.0);
    assert_vec_close(force, [20.0, 0.0, 0.0], 1e-12);
}

#[test]
fn atom_sitting_on_the_anchor_feels_no_force() {
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0).with_anchor([0.0; 3], 0.0);
    // At t = 2 the anchor is at x = 2.
    let force = spring_force(&mut smd, [2.0, 0.0, 0.0], 2.0);
    assert_vec_close(force, [0.0; 3], 1e-12);
    assert_eq!(smd.work_accumulated, 0.0);
}

#[test]
fn force_is_a_restoring_spring_in_the_atom_position() {
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0).with_anchor([0.0; 3], 0.0);
    let samples = [-1.0, 0.0, 1.0, 2.0, 3.0];
    let forces: Vec<f64> = samples
        .iter()
        .map(|&x| spring_force(&mut smd, [x, 0.0, 0.0], 2.0)[0])
        .collect();
    // F(x) = k (anchor - x) with the anchor at x = 2: linear with slope -k.
    for (&x, &f) in samples.iter().zip(forces.iter()) {
        assert!((f - 10.0 * (2.0 - x)).abs() < 1e-12, "F({x}) = {f}");
    }
    for pair in forces.windows(2) {
        assert!(pair[1] < pair[0], "force must fall as the atom advances");
    }
    let slope = (forces[4] - forces[0]) / (samples[4] - samples[0]);
    assert!((slope + 10.0).abs() < 1e-12, "dF/dx = {slope}, expected -k");
    // The atom overtaking the anchor is pulled back.
    assert!(forces[4] < 0.0);
}

#[test]
fn transverse_displacement_is_pulled_back_towards_the_anchor_line() {
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0).with_anchor([0.0; 3], 0.0);
    // Anchor at (1, 0, 0) after one time unit; the atom sits at (1, 0.5, -0.25).
    let force = spring_force(&mut smd, [1.0, 0.5, -0.25], 1.0);
    assert_vec_close(force, [0.0, -5.0, 2.5], 1e-12);
}

#[test]
fn lazily_pinned_anchor_starts_at_the_first_position_and_time() {
    let mut smd = SteeredMd::new(0, [0.0, 0.0, 1.0], 4.0, 0.5);
    assert_eq!(smd.anchor_origin, None);
    assert_eq!(smd.anchor_position(3.0), None);

    // First call at t = 3 with the atom at (1, 2, 3): the spring is relaxed.
    let start = [1.0, 2.0, 3.0];
    let first = spring_force(&mut smd, start, 3.0);
    assert_vec_close(first, [0.0; 3], 1e-12);
    assert_eq!(smd.anchor_origin, Some(start));
    assert_eq!(smd.start_time, 3.0);

    // Two time units later the anchor has moved 0.5 * 2 = 1 along z while the atom
    // stayed put: F = k * 1 = 4 along z.
    let second = spring_force(&mut smd, start, 5.0);
    assert_vec_close(second, [0.0, 0.0, 4.0], 1e-12);
    assert_eq!(smd.anchor_position(5.0), Some([1.0, 2.0, 4.0]));
}

#[test]
fn anchor_does_not_follow_the_atom() {
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 10.0, 1.0);
    spring_force(&mut smd, [0.0; 3], 0.0);
    // At t = 1 the anchor is at x = 1.  An atom that has advanced to 0.5 is still
    // 0.5 behind it; one at 0.9 is 0.1 behind.  A spring whose anchor was rebuilt
    // from the atom's position would give k * v * t = 10 in both cases.
    let behind_half = spring_force(&mut smd, [0.5, 0.0, 0.0], 1.0);
    let behind_tenth = spring_force(&mut smd, [0.9, 0.0, 0.0], 1.0);
    assert_vec_close(behind_half, [5.0, 0.0, 0.0], 1e-12);
    assert_vec_close(behind_tenth, [1.0, 0.0, 0.0], 1e-12);
}

#[test]
fn direction_is_normalised_for_the_anchor_motion() {
    // Direction (3, 0, 4) normalises to (0.6, 0, 0.8); the anchor travels 2 * 5 = 10.
    let mut smd = SteeredMd::new(0, [3.0, 0.0, 4.0], 2.0, 2.0).with_anchor([0.0; 3], 0.0);
    let force = spring_force(&mut smd, [0.0; 3], 5.0);
    assert_vec_close(force, [12.0, 0.0, 16.0], 1e-12);
    let travelled = smd
        .anchor_position(5.0)
        .map(|a| (a[0] * a[0] + a[2] * a[2]).sqrt());
    assert!(
        travelled.is_some_and(|d| (d - 10.0).abs() < 1e-12),
        "anchor travelled {travelled:?}"
    );
}

#[test]
fn explicit_anchor_is_not_repinned_by_the_first_call() {
    // The anchor was placed at (1, 0, 0) at t = 0; the atom is already 0.25 beyond it
    // when the first force is requested, so the spring pulls it back immediately.
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 8.0, 0.0).with_anchor([1.0, 0.0, 0.0], 0.0);
    let force = spring_force(&mut smd, [1.25, 0.0, 0.0], 0.0);
    assert_vec_close(force, [-2.0, 0.0, 0.0], 1e-12);
    assert_eq!(smd.anchor_origin, Some([1.0, 0.0, 0.0]));
}

#[test]
fn anchor_start_time_offsets_the_motion() {
    let smd = SteeredMd::new(0, [0.0, 1.0, 0.0], 1.0, 3.0).with_anchor([1.0, 1.0, 1.0], 5.0);
    assert_eq!(smd.anchor_position(5.0), Some([1.0, 1.0, 1.0]));
    assert_eq!(smd.anchor_position(7.0), Some([1.0, 7.0, 1.0]));
}

#[test]
fn zero_velocity_is_a_position_restraint_at_the_initial_position() {
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], 6.0, 0.0);
    spring_force(&mut smd, [1.0, 2.0, 3.0], 0.0);
    let force = spring_force(&mut smd, [1.5, 2.0, 2.0], 100.0);
    assert_vec_close(force, [-3.0, 0.0, 6.0], 1e-12);
    assert_eq!(smd.work_accumulated, 0.0);
}

#[test]
fn work_of_a_held_atom_is_the_rectangle_sum_of_k_v_squared_t() {
    // Atom fixed at the origin: F(t) = k v t and dW = F v dt.  Summing the n steps
    // at t_i = i dt gives k v^2 dt^2 n (n - 1) / 2, the rectangle approximation of
    // the integral k v^2 T^2 / 2.
    let (k, v, dt, n) = (10.0, 2.0, 0.01, 100_usize);
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], k, v).with_anchor([0.0; 3], 0.0);
    let positions = [[0.0_f64; 3]];
    let mut forces = [[0.0_f64; 3]];
    for step in 0..n {
        smd.apply_force(&positions, &mut forces, step as f64 * dt, dt);
    }
    let expected = k * v * v * dt * dt * (n * (n - 1)) as f64 / 2.0;
    assert!(
        (smd.work_accumulated - expected).abs() < 1e-9,
        "W = {}, expected {expected}",
        smd.work_accumulated
    );
    let continuum = 0.5 * k * v * v * (n as f64 * dt).powi(2);
    assert!(smd.work_accumulated < continuum && smd.work_accumulated > 0.98 * continuum);
}

#[test]
fn dragged_damped_atom_lags_the_anchor_by_gamma_v_over_k() {
    // m x'' = F_spring - gamma x'.  The pulled atom settles to the anchor velocity v
    // with a constant lag, so the spring force tends to gamma * v and the lag to
    // gamma * v / k.  A spring whose anchor follows the atom would instead grow as
    // k v t without bound.
    let (mass, gamma, k, v, dt) = (1.0, 2.0, 10.0, 0.5, 1.0e-3);
    let mut smd = SteeredMd::new(0, [1.0, 0.0, 0.0], k, v);
    let mut position = [0.0_f64; 3];
    let mut velocity = [0.0_f64; 3];
    let steps = 20_000;
    let mut spring = [0.0_f64; 3];
    for step in 0..steps {
        let mut forces = [[0.0_f64; 3]];
        smd.apply_force(&[position], &mut forces, step as f64 * dt, dt);
        spring = forces[0];
        for axis in 0..3 {
            velocity[axis] += dt * (spring[axis] - gamma * velocity[axis]) / mass;
            position[axis] += dt * velocity[axis];
        }
    }
    assert!(
        (spring[0] - gamma * v).abs() < 1e-3,
        "steady spring force {} should equal gamma v = {}",
        spring[0],
        gamma * v
    );
    assert!((velocity[0] - v).abs() < 1e-3, "v = {}", velocity[0]);
    let anchor = smd.anchor_position(steps as f64 * dt);
    let lag = anchor.map(|a| a[0] - position[0]);
    assert!(
        lag.is_some_and(|l| (l - gamma * v / k).abs() < 1e-3),
        "lag = {lag:?}, expected {}",
        gamma * v / k
    );
    assert!(spring[1].abs() < 1e-12 && spring[2].abs() < 1e-12);
}

#[test]
fn missing_pulled_atom_is_ignored() {
    let mut smd = SteeredMd::new(3, [1.0, 0.0, 0.0], 10.0, 1.0);
    let positions = [[0.0_f64; 3]];
    let mut forces = [[1.0_f64, 2.0, 3.0]];
    let work = smd.apply_force(&positions, &mut forces, 4.0, 0.1);
    assert_eq!(work, 0.0);
    assert_eq!(forces, [[1.0, 2.0, 3.0]]);
    assert_eq!(smd.anchor_origin, None);
    assert_eq!(smd.time_elapsed, 0.0);
    assert_eq!(smd.work_accumulated, 0.0);
}

#[test]
fn forces_of_other_atoms_are_left_untouched() {
    let mut smd = SteeredMd::new(1, [1.0, 0.0, 0.0], 10.0, 1.0).with_anchor([0.0; 3], 0.0);
    let positions = [[5.0, 5.0, 5.0], [0.0, 0.0, 0.0], [-3.0, 1.0, 2.0]];
    let mut forces = [[0.5_f64; 3]; 3];
    smd.apply_force(&positions, &mut forces, 1.0, 0.001);
    assert_eq!(forces[0], [0.5; 3]);
    assert_eq!(forces[2], [0.5; 3]);
    assert_vec_close(forces[1], [10.5, 0.5, 0.5], 1e-12);
}
