// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the impulse module: contact-normal convention,
//! friction signs and robustness of the sequential-impulse solvers against
//! malformed constraint indices.
//!
//! All bodies use the identity orientation and isotropic inertia so the
//! expected values are independent of how local inverse inertia is mapped to
//! the world frame.

use super::functions::*;
use super::types::*;

fn body(pos: [f64; 3], vel: [f64; 3], inv_mass: f64, inv_inertia: f64) -> RigidBodyState {
    RigidBodyState {
        position: pos,
        velocity: vel,
        angular_velocity: [0.0; 3],
        orientation: [0.0, 0.0, 0.0, 1.0],
        inv_mass,
        inv_inertia_local: [inv_inertia; 3],
    }
}

fn config(iterations: usize, warm_start: bool) -> SequentialImpulseSolverConfig {
    SequentialImpulseSolverConfig {
        velocity_iterations: iterations,
        baumgarte: 0.0,
        dt: 1.0 / 60.0,
        warm_start,
    }
}

/// Kinetic energy of a body with isotropic local inertia.
fn kinetic_energy(b: &RigidBodyState) -> f64 {
    let linear = if b.inv_mass > 0.0 {
        0.5 * dot3(b.velocity, b.velocity) / b.inv_mass
    } else {
        0.0
    };
    let angular = if b.inv_inertia_local[0] > 0.0 {
        0.5 * dot3(b.angular_velocity, b.angular_velocity) / b.inv_inertia_local[0]
    } else {
        0.0
    };
    linear + angular
}

/// Dynamic point mass above a static ground; the contact normal (from the
/// dynamic body A to the ground B) points down.
fn sliding_scene(a_velocity: [f64; 3]) -> Vec<RigidBodyState> {
    vec![
        body([0.0, 1.0, 0.0], a_velocity, 1.0, 0.0),
        body([0.0, 0.0, 0.0], [0.0; 3], 0.0, 0.0),
    ]
}

/// A ground contact between `body_a` (above) and `body_b` (below).
fn ground_constraint(
    body_a: usize,
    body_b: usize,
    contact: [f64; 3],
    friction: f64,
) -> ContactConstraint {
    ContactConstraint::new(
        body_a,
        body_b,
        [0.0, -1.0, 0.0],
        contact,
        0.0,
        0.0,
        friction,
    )
}

/// Two unit point masses on the x axis approaching each other at 1 m/s each.
fn head_on_pair() -> Vec<RigidBodyState> {
    vec![
        body([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0),
        body([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0),
    ]
}

fn x_contact(body_a: usize, body_b: usize) -> ContactConstraint {
    ContactConstraint::new(
        body_a,
        body_b,
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        0.0,
        0.0,
        0.3,
    )
}

// ── Normal convention (A -> B) ──────────────────────────────────────────

#[test]
fn one_manifold_separates_through_impulse_and_positional_correction() {
    // A at x = 0 moving +x, B at x = 1.9 moving -x; spheres of radius 1
    // overlap by 0.1.  The normal points from A to B.
    let manifold = CollisionManifold {
        normal: [1.0, 0.0, 0.0],
        depth: 0.1,
        contact_point: [0.95, 0.0, 0.0],
    };
    let mut a = body([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0);
    let mut b = body([1.9, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0);

    resolve_collision(&mut a, &mut b, &manifold, 0.5, 0.0);
    let closing = dot3(sub3(a.velocity, b.velocity), manifold.normal);
    assert!(
        (closing + 1.0).abs() < 1e-12,
        "e = 0.5 must turn the closing speed 2 into a separation speed 1, got {closing}"
    );
    assert!(a.velocity[0] < 0.0 && b.velocity[0] > 0.0);

    let gap_before = b.position[0] - a.position[0];
    PositionalCorrector::new(0.8, 0.0).correct(&mut a, &mut b, manifold.normal, manifold.depth);
    let gap_after = b.position[0] - a.position[0];
    assert!(
        (gap_after - gap_before - 0.8 * 0.1).abs() < 1e-12,
        "positional correction must open the gap by beta * depth: {gap_before} -> {gap_after}"
    );
}

#[test]
fn positional_correction_moves_dynamic_body_only_along_minus_normal() {
    let mut a = body([0.0; 3], [0.0; 3], 1.0, 0.0);
    let mut ground = body([0.0, -1.0, 0.0], [0.0; 3], 0.0, 0.0);
    // A rests on the ground (the ground is body B, the normal points down).
    PositionalCorrector::new(1.0, 0.0).correct(&mut a, &mut ground, [0.0, -1.0, 0.0], 0.2);
    assert!(
        (a.position[1] - 0.2).abs() < 1e-12,
        "a.y = {}",
        a.position[1]
    );
    assert!(a.position[0].abs() < 1e-15 && a.position[2].abs() < 1e-15);
    assert_eq!(ground.position, [0.0, -1.0, 0.0]);
}

// ── Friction signs ──────────────────────────────────────────────────────

#[test]
fn sequential_impulse_friction_decelerates_a_sliding_body() {
    let mut bodies = sliding_scene([2.0, -1.0, 0.0]);
    let mut constraints = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 0.5)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(10, false));

    // The normal impulse removes the approach speed; friction removes mu * J_n.
    assert!((constraints[0].accumulated_normal_impulse - 1.0).abs() < 1e-12);
    let vx = bodies[0].velocity[0];
    assert!(
        (vx - 1.5).abs() < 1e-12,
        "sliding speed must drop from 2 to 2 - mu * J_n = 1.5, got {vx}"
    );
    assert!(bodies[0].velocity[1].abs() < 1e-12);
    assert!(bodies[0].velocity[2].abs() < 1e-12);
    assert_eq!(bodies[1].velocity, [0.0; 3], "static ground must not move");
}

#[test]
fn sequential_impulse_friction_brakes_both_tangent_directions() {
    let mut bodies = sliding_scene([2.0, -1.0, 1.5]);
    let mut constraints = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 0.5)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(10, false));
    // Each tangent direction is limited to mu * J_n = 0.5.
    let v = bodies[0].velocity;
    assert!((v[0] - 1.5).abs() < 1e-12, "{v:?}");
    assert!((v[2] - 1.0).abs() < 1e-12, "{v:?}");
}

#[test]
fn sequential_impulse_friction_can_stop_the_slide_but_never_reverse_it() {
    let mut bodies = sliding_scene([2.0, -1.0, 0.0]);
    let mut constraints = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 5.0)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(20, false));
    assert!(
        bodies[0].velocity[0].abs() < 1e-12,
        "large friction must stick, vx = {}",
        bodies[0].velocity[0]
    );
    assert!(bodies[0].velocity[2].abs() < 1e-12);
}

#[test]
fn negative_and_nan_friction_coefficients_act_frictionless_instead_of_panicking() {
    for friction in [-0.5, f64::NAN, f64::NEG_INFINITY] {
        let mut bodies = sliding_scene([2.0, -1.0, 0.5]);
        let mut constraints = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], friction)];
        solve_sequential_impulses(&mut bodies, &mut constraints, &config(10, false));
        // The normal impulse still removes the approach speed (J_n = 1) ...
        assert!(
            (constraints[0].accumulated_normal_impulse - 1.0).abs() < 1e-12,
            "friction = {friction}"
        );
        // ... and no tangential impulse is applied.
        assert_eq!(constraints[0].accumulated_friction_impulse, [0.0, 0.0]);
        let v = bodies[0].velocity;
        assert!(
            (v[0] - 2.0).abs() < 1e-12 && v[1].abs() < 1e-12 && (v[2] - 0.5).abs() < 1e-12,
            "friction = {friction}: {v:?}"
        );
    }
}

/// Sphere of radius 1 (m = 1, I = 0.4 m r^2, so the inverse inertia is 2.5)
/// whose centre is at y = 1 so that the ground contact is the origin.
fn spinning_sphere(velocity: [f64; 3], angular_velocity: [f64; 3]) -> Vec<RigidBodyState> {
    let mut sphere = body([0.0, 1.0, 0.0], velocity, 1.0, 2.5);
    sphere.angular_velocity = angular_velocity;
    vec![sphere, body([0.0; 3], [0.0; 3], 0.0, 0.0)]
}

fn contact_slip(bodies: &[RigidBodyState]) -> f64 {
    let v = relative_velocity_at_contact(&bodies[0], &bodies[1], [0.0; 3]);
    v[0].hypot(v[2])
}

#[test]
fn friction_with_rotation_never_adds_energy_and_reduces_slip() {
    let mut bodies = spinning_sphere([3.0, -0.5, 0.5], [0.4, 0.0, -1.0]);
    let energy_before = kinetic_energy(&bodies[0]);
    // Contact slip of the initial state: v + w x r = (2, -0.5, 0.1).
    let slip_before = contact_slip(&bodies);
    assert!(
        (slip_before - 4.01f64.sqrt()).abs() < 1e-12,
        "{slip_before}"
    );

    let mut constraints = vec![ground_constraint(0, 1, [0.0; 3], 0.6)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(30, false));

    let energy_after = kinetic_energy(&bodies[0]);
    assert!(
        energy_after <= energy_before + 1e-12,
        "friction injected energy: {energy_before} -> {energy_after}"
    );
    // J_n = 0.5, so each tangent impulse is capped at 0.3.  The t1 slip
    // (0.1) is removed completely; the t2 slip (2.0) drops by 0.3 * 3.5.
    let slip_after = contact_slip(&bodies);
    assert!(
        (slip_after - 0.95).abs() < 1e-9,
        "slip after = {slip_after}"
    );
}

#[test]
fn rolling_sphere_with_large_friction_reaches_zero_contact_slip() {
    let mut bodies = spinning_sphere([3.0, -0.5, 0.0], [0.0, 0.0, -1.0]);
    let mut constraints = vec![ground_constraint(0, 1, [0.0; 3], 10.0)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(10, false));
    let slip = contact_slip(&bodies);
    assert!(slip < 1e-9, "contact slip = {slip}");
    // Spin and translation are now locked together (rolling without slipping):
    // v_x = -w_z * r with r = 1.
    let v = bodies[0].velocity;
    let w = bodies[0].angular_velocity;
    assert!((v[0] + w[2]).abs() < 1e-9, "v = {v:?}, w = {w:?}");
}

#[test]
fn coulomb_friction_on_a_spinning_sphere_brakes_the_contact_slip() {
    let manifold = CollisionManifold {
        normal: [0.0, -1.0, 0.0],
        depth: 0.0,
        contact_point: [0.0; 3],
    };
    // Contact slip 2 m/s (3 - 1); tangential inverse mass 1 + 2.5 = 3.5.
    let mut bodies = spinning_sphere([3.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
    assert!((contact_slip(&bodies) - 2.0).abs() < 1e-12);
    let limited = ImpulseSolver::new(0.0, 0.5, 0.0);
    {
        let (sphere, ground) = bodies.split_at_mut(1);
        // mu * J_n = 0.5 * 0.2 = 0.1 < 2 / 3.5: the Coulomb limit is active.
        limited.apply_coulomb_friction(&mut sphere[0], &mut ground[0], &manifold, 0.2);
    }
    let v = bodies[0].velocity;
    let w = bodies[0].angular_velocity;
    assert!(
        (v[0] - 2.9).abs() < 1e-12 && v[1].abs() < 1e-12 && v[2].abs() < 1e-12,
        "{v:?}"
    );
    assert!(
        (w[2] + 1.25).abs() < 1e-12 && w[0].abs() < 1e-12 && w[1].abs() < 1e-12,
        "{w:?}"
    );
    assert!((contact_slip(&bodies) - 1.65).abs() < 1e-12);
    assert_eq!(bodies[1].velocity, [0.0; 3], "static ground must not move");

    // A friction limit above the sticking impulse 2 / 3.5 locks the contact.
    let mut bodies = spinning_sphere([3.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
    let sticky = ImpulseSolver::new(0.0, 10.0, 0.0);
    {
        let (sphere, ground) = bodies.split_at_mut(1);
        sticky.apply_coulomb_friction(&mut sphere[0], &mut ground[0], &manifold, 100.0);
    }
    let impulse = 2.0 / 3.5;
    assert!(contact_slip(&bodies) < 1e-12);
    assert!((bodies[0].velocity[0] - (3.0 - impulse)).abs() < 1e-12);
    assert!((bodies[0].angular_velocity[2] - (-1.0 - 2.5 * impulse)).abs() < 1e-12);
}

#[test]
fn warm_started_second_frame_reproduces_the_converged_friction_solution() {
    let mut solver = WarmStartSolver::new(config(10, true));

    let mut frame1 = sliding_scene([2.0, -1.0, 0.0]);
    let mut c1 = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 0.5)];
    solver.solve(&mut frame1, &mut c1);
    assert!(
        (frame1[0].velocity[0] - 1.5).abs() < 1e-12,
        "{:?}",
        frame1[0].velocity
    );
    assert_eq!(solver.cache.len(), 1);
    assert!((solver.cache[0].lambda_n - 1.0).abs() < 1e-12);
    assert!(
        (solver.cache[0].lambda_t1.hypot(solver.cache[0].lambda_t2) - 0.5).abs() < 1e-12,
        "cached friction impulse must sit on the Coulomb limit"
    );

    // Same scene again: the cached impulses are applied by the warm-start path
    // with the same sign convention as the iteration, so the result must match.
    let mut frame2 = sliding_scene([2.0, -1.0, 0.0]);
    let mut c2 = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 0.5)];
    solver.solve(&mut frame2, &mut c2);
    for k in 0..3 {
        assert!(
            (frame2[0].velocity[k] - frame1[0].velocity[k]).abs() < 1e-12,
            "component {k}: {:?} vs {:?}",
            frame2[0].velocity,
            frame1[0].velocity
        );
    }
}

// ── Body index robustness ───────────────────────────────────────────────

#[test]
fn constraint_with_body_a_after_body_b_resolves_without_panic() {
    // bodies[1] is "A" (x = 0, moving +x); bodies[0] is "B" (x = 2, moving -x),
    // so the normal +x points from A (index 1) to B (index 0).
    let mut bodies = head_on_pair();
    bodies.swap(0, 1);
    let mut constraints = vec![x_contact(1, 0)];
    solve_sequential_impulses(&mut bodies, &mut constraints, &config(10, false));
    assert!(
        bodies[0].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[0].velocity
    );
    assert!(
        bodies[1].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[1].velocity
    );

    let mut bodies = head_on_pair();
    bodies.swap(0, 1);
    let mut constraints = vec![x_contact(1, 0)];
    let mut solver = WarmStartSolver::new(config(10, true));
    solver.solve(&mut bodies, &mut constraints);
    assert!(
        bodies[0].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[0].velocity
    );
    assert!(
        bodies[1].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[1].velocity
    );
    assert_eq!(solver.cache.len(), 1);
    assert_eq!(solver.cache[0].key, (1, 0));
}

#[test]
fn self_contacts_and_out_of_range_indices_are_skipped() {
    let malformed = [(1usize, 1usize), (0, 2), (2, 0), (7, 9), (usize::MAX, 0)];
    for &(a, b) in &malformed {
        let mut bodies = head_on_pair();
        let mut constraints = vec![x_contact(a, b)];
        solve_sequential_impulses(&mut bodies, &mut constraints, &config(5, true));
        assert_eq!(bodies[0].velocity, [1.0, 0.0, 0.0], "pair ({a}, {b})");
        assert_eq!(bodies[1].velocity, [-1.0, 0.0, 0.0], "pair ({a}, {b})");
        assert_eq!(constraints[0].accumulated_normal_impulse, 0.0);

        let mut solver = WarmStartSolver::new(config(5, true));
        let mut bodies = head_on_pair();
        let mut constraints = vec![x_contact(a, b)];
        solver.solve(&mut bodies, &mut constraints);
        assert_eq!(bodies[0].velocity, [1.0, 0.0, 0.0], "pair ({a}, {b})");
        assert_eq!(bodies[1].velocity, [-1.0, 0.0, 0.0], "pair ({a}, {b})");
        assert!(
            solver.cache.is_empty(),
            "pair ({a}, {b}) must not create a cache entry"
        );
    }
}

#[test]
fn valid_constraints_are_still_solved_next_to_malformed_ones() {
    let mut bodies = head_on_pair();
    let mut constraints = vec![x_contact(0, 0), x_contact(0, 1), x_contact(3, 1)];
    let mut solver = WarmStartSolver::new(config(10, true));
    solver.solve(&mut bodies, &mut constraints);
    assert!(
        bodies[0].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[0].velocity
    );
    assert!(
        bodies[1].velocity[0].abs() < 1e-12,
        "{:?}",
        bodies[1].velocity
    );
    assert_eq!(solver.cache.len(), 1);
    assert_eq!(solver.cache[0].key, (0, 1));
}

#[test]
fn pair_helpers_reject_equal_and_out_of_range_indices() {
    let mut bodies = vec![
        body([0.0; 3], [0.0; 3], 1.0, 0.0),
        body([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0),
        body([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0),
    ];
    assert!(pair_mut(&mut bodies, 1, 1).is_none());
    assert!(pair_mut(&mut bodies, 0, 3).is_none());
    assert!(pair_mut(&mut bodies, 3, 0).is_none());
    assert!(pair_ref(&bodies, 2, 2).is_none());
    assert!(pair_ref(&bodies, 2, 9).is_none());

    let (first, second) = pair_mut(&mut bodies, 2, 0).expect("distinct in-range indices");
    assert_eq!(first.position[0], 2.0);
    assert_eq!(second.position[0], 0.0);
    let (first, second) = pair_ref(&bodies, 0, 2).expect("distinct in-range indices");
    assert_eq!(first.position[0], 0.0);
    assert_eq!(second.position[0], 2.0);
}

#[test]
fn cached_impulses_follow_a_contact_whose_bodies_swap_roles() {
    // Frame 1: the box (index 0) is body A, the ground (index 1) is body B.
    let mut solver = WarmStartSolver::new(config(10, true));
    let mut frame1 = sliding_scene([2.0, -1.0, 0.0]);
    let mut c1 = vec![ground_constraint(0, 1, [0.0, 0.5, 0.0], 0.5)];
    solver.solve(&mut frame1, &mut c1);
    assert!((frame1[0].velocity[0] - 1.5).abs() < 1e-12);

    // Frame 2: the same two bodies, but the ground is now body A and the box
    // body B, so the normal points up.  With zero iterations the result is
    // purely the warm start, which must brake the box exactly like frame 1.
    let mut warm_only = WarmStartSolver::new(config(0, true));
    warm_only.cache = solver.cache.clone();
    let mut frame2 = sliding_scene([2.0, -1.0, 0.0]);
    let mut c2 = vec![ContactConstraint::new(
        1,
        0,
        [0.0, 1.0, 0.0],
        [0.0, 0.5, 0.0],
        0.0,
        0.0,
        0.5,
    )];
    warm_only.solve(&mut frame2, &mut c2);
    let v = frame2[0].velocity;
    assert!(
        (v[0] - 1.5).abs() < 1e-12 && v[1].abs() < 1e-12 && v[2].abs() < 1e-12,
        "swapped-order warm start must reproduce the converged solution, got {v:?}"
    );

    // Full solve in the swapped order also converges to the same velocity and
    // writes the cache back without disturbing the stored orientation.
    let mut frame3 = sliding_scene([2.0, -1.0, 0.0]);
    let mut c3 = vec![ContactConstraint::new(
        1,
        0,
        [0.0, 1.0, 0.0],
        [0.0, 0.5, 0.0],
        0.0,
        0.0,
        0.5,
    )];
    solver.solve(&mut frame3, &mut c3);
    assert!((frame3[0].velocity[0] - 1.5).abs() < 1e-12);
    assert_eq!(solver.cache.len(), 1);
    assert_eq!(solver.cache[0].key, (0, 1));
    assert!((solver.cache[0].lambda_n - 1.0).abs() < 1e-12);
    assert!(
        (solver.cache[0].lambda_t2 + 0.5).abs() < 1e-12,
        "cached lambda_t2 must keep the orientation of its key: {:?}",
        solver.cache[0]
    );
}
