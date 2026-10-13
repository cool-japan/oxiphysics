// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`RopeConstraint`]: position-based velocity
//! reconciliation, damping and the time-step guard.

use super::*;

/// Largest speed over the un-pinned rope nodes.
fn max_unpinned_speed(rope: &RopeConstraint) -> f64 {
    rope.nodes
        .iter()
        .filter(|n| !n.pinned)
        .map(|n| v3_len(n.velocity))
        .fold(0.0, f64::max)
}

/// Taut 5 m rope with 10 segments, both ends pinned, simulated for 10 s at
/// 60 Hz.  Returns the settled rope and its initial length.
fn settle_taut_rope(solver_iterations: usize) -> (RopeConstraint, f64) {
    let mut rope = RopeConstraint::new(
        [0.0; 3],
        [5.0, 0.0, 0.0],
        10,
        0.2,
        1.0,
        0.1,
        [0.0, -9.81, 0.0],
    );
    rope.solver_iterations = solver_iterations;
    let initial_length = rope.current_length();
    let dt = 1.0 / 60.0;
    for _ in 0..600 {
        rope.step(dt);
    }
    (rope, initial_length)
}

#[test]
fn test_rope_pinned_at_both_ends_comes_to_rest() {
    // A taut rope fixed at both ends must settle: the velocity has to be
    // reconciled with the projected positions every step, otherwise
    // `g * dt` accumulates in every node forever (98 m/s after 10 s).
    let (rope, initial_length) = settle_taut_rope(4);
    let speed = max_unpinned_speed(&rope);
    assert!(speed < 0.1, "rope must settle, max node speed {speed} m/s");
    // Position-based dynamics with the default 4 Gauss-Seidel sweeps leaves
    // a finite residual stretch (measured: 1.8 % for this taut rope).
    let length = rope.current_length();
    assert!(
        length <= 1.025 * initial_length,
        "rope over-stretched: {length} vs initial {initial_length}"
    );
    assert!(length > initial_length, "gravity must stretch a taut rope");
    // The pinned ends never moved and never gained velocity.
    assert_eq!(rope.nodes[0].position, [0.0; 3]);
    assert_eq!(rope.nodes[10].position, [5.0, 0.0, 0.0]);
    assert_eq!(rope.nodes[0].velocity, [0.0; 3]);
    assert_eq!(rope.nodes[10].velocity, [0.0; 3]);
    // Gravity pulled the middle of the rope below the anchor line.
    assert!(rope.nodes[5].position[1] < 0.0);
}

#[test]
fn test_rope_stretch_shrinks_with_more_solver_iterations() {
    // The residual stretch is a convergence property of the projection:
    // with 16 sweeps the rope stays within 1 % of its rest length, and the
    // stretch decreases monotonically with the iteration count.
    let (coarse, l0) = settle_taut_rope(1);
    let (default, _) = settle_taut_rope(4);
    let (stiff, _) = settle_taut_rope(16);
    let s_coarse = coarse.current_length() / l0;
    let s_default = default.current_length() / l0;
    let s_stiff = stiff.current_length() / l0;
    assert!(s_stiff < s_default && s_default < s_coarse);
    assert!(s_stiff <= 1.01, "16 sweeps must stay within 1 %: {s_stiff}");
    assert!(max_unpinned_speed(&stiff) < 0.1);
}

#[test]
fn test_rope_velocity_matches_projected_displacement() {
    // Undamped: v = (x_after - x_before) / dt for every un-pinned node,
    // including the effect of the constraint projection.
    let mut rope = RopeConstraint::new(
        [0.0; 3],
        [3.0, 0.0, 0.0],
        6,
        0.5,
        1.0,
        0.0,
        [0.0, -9.81, 0.0],
    );
    let dt = 1.0 / 60.0;
    for _ in 0..5 {
        let before: Vec<[f64; 3]> = rope.nodes.iter().map(|n| n.position).collect();
        rope.step(dt);
        for (node, prev) in rope.nodes.iter().zip(before.iter()) {
            if node.pinned {
                continue;
            }
            for (k, ((position, previous), velocity)) in node
                .position
                .iter()
                .zip(prev.iter())
                .zip(node.velocity.iter())
                .enumerate()
            {
                let expected = (position - previous) / dt;
                assert!(
                    (velocity - expected).abs() < 1e-12,
                    "axis {k}: v={velocity} expected={expected}"
                );
            }
        }
    }
}

#[test]
fn test_rope_damping_dissipates_motion() {
    let make = |damping: f64| {
        let mut rope = RopeConstraint::new(
            [0.0; 3],
            [2.0, 0.0, 0.0],
            8,
            0.1,
            1.0,
            damping,
            [0.0, -9.81, 0.0],
        );
        // Release the far end so the rope swings like a pendulum.
        if let Some(last) = rope.nodes.last_mut() {
            last.pinned = false;
        }
        rope
    };
    let dt = 1.0 / 120.0;
    let mut free = make(0.0);
    let mut damped = make(20.0);
    let mut free_peak = 0.0_f64;
    let mut damped_peak = 0.0_f64;
    for step in 0..600 {
        free.step(dt);
        damped.step(dt);
        if step >= 300 {
            free_peak = free_peak.max(max_unpinned_speed(&free));
            damped_peak = damped_peak.max(max_unpinned_speed(&damped));
        }
    }
    assert!(
        damped_peak < 0.5 * free_peak,
        "damping must dissipate: damped {damped_peak} vs undamped {free_peak}"
    );
}

#[test]
fn test_rope_free_fall_velocity_is_gravity_times_time() {
    // Both nodes free: the rope falls rigidly and the velocity must equal
    // g * t (the reconciliation must not alter unconstrained motion).
    let mut rope = RopeConstraint::new(
        [0.0; 3],
        [1.0, 0.0, 0.0],
        1,
        1.0,
        1.0,
        0.0,
        [0.0, -9.81, 0.0],
    );
    rope.nodes[0].pinned = false;
    rope.nodes[1].pinned = false;
    let dt = 1.0 / 60.0;
    for _ in 0..60 {
        rope.step(dt);
    }
    for node in &rope.nodes {
        assert!(
            (node.velocity[1] + 9.81).abs() < 1e-9,
            "vy={}",
            node.velocity[1]
        );
        assert!(node.velocity[0].abs() < 1e-9);
    }
    assert!((rope.current_length() - 1.0).abs() < 1e-9);
}

#[test]
fn test_rope_step_ignores_non_positive_time_step() {
    let mut rope = RopeConstraint::new(
        [0.0; 3],
        [1.0, 0.0, 0.0],
        2,
        1.0,
        1.0,
        0.1,
        [0.0, -9.81, 0.0],
    );
    let before: Vec<([f64; 3], [f64; 3])> = rope
        .nodes
        .iter()
        .map(|n| (n.position, n.velocity))
        .collect();
    rope.step(0.0);
    rope.step(-0.1);
    for (node, (pos, vel)) in rope.nodes.iter().zip(before.iter()) {
        assert_eq!(node.position, *pos);
        assert_eq!(node.velocity, *vel);
    }
}
