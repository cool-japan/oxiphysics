// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the spring-like joints: they are soft constraints
//! whose impulse does not depend on how many solver iterations run.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::{RigidBody, RigidBodySet};

use super::test_support::{linear_momentum, step_joint, total_kinetic_energy};
use super::{BreakableJoint, SpringJoint, XpbdJoint};
use crate::pgs_solver::PgsSolver;
use crate::traits::Constraint;

const DT: Real = 1.0 / 60.0;
const GRAVITY: Real = 9.81;

/// Two undamped unit-mass bodies on the x axis, `distance` apart.
fn pair_on_x_axis(distance: Real) -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(RigidBody::new(1.0));
    let mut body_b = RigidBody::new(1.0);
    body_b.transform.position = Vec3::new(distance, 0.0, 0.0);
    let handle_b = bodies.insert(body_b);
    (bodies, handle_a, handle_b)
}

/// A static anchor at the origin and a 1 kg bob hanging `drop` below it.
fn hanging_bob(drop: Real) -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(RigidBody::new_static());
    let mut bob = RigidBody::new(1.0);
    bob.transform.position = Vec3::new(0.0, -drop, 0.0);
    let handle_b = bodies.insert(bob);
    (bodies, handle_a, handle_b)
}

fn gravity() -> Vec3 {
    Vec3::new(0.0, -GRAVITY, 0.0)
}

// ---------------------------------------------------------------------------
// F0385: SpringJoint is an implicit, iteration-independent soft constraint
// ---------------------------------------------------------------------------

#[test]
fn spring_first_iteration_applies_the_implicit_euler_impulse() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
    let (stiffness, rest_length) = (100.0, 2.0);
    let mut spring = SpringJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        rest_length,
        stiffness,
        0.0,
    );
    spring.prepare(&bodies, DT);
    spring.solve_velocity(&mut bodies, DT);
    // Δv = dt·k·C / (1 + dt²·k·K) with the inverse effective mass K = 2.
    let stretch = 1.0;
    let expected = DT * stiffness * stretch / (1.0 + DT * DT * stiffness * 2.0);
    let velocity_a = bodies.get(handle_a).unwrap().velocity;
    let velocity_b = bodies.get(handle_b).unwrap().velocity;
    assert!(
        (velocity_a.x - expected).abs() < 1e-12,
        "A should gain {expected}, got {}",
        velocity_a.x
    );
    assert!(
        (velocity_b.x + expected).abs() < 1e-12,
        "B should lose {expected}, got {}",
        velocity_b.x
    );
    // The first iteration of a lone spring is exact: the second changes nothing.
    assert!(spring.residual() > 1.0);
    spring.solve_velocity(&mut bodies, DT);
    assert!(spring.residual() < 1e-12);
    assert!((bodies.get(handle_a).unwrap().velocity.x - expected).abs() < 1e-12);
}

#[test]
fn spring_impulse_does_not_depend_on_the_iteration_count() {
    let anchors = (Vec3::new(0.2, 0.1, 0.0), Vec3::new(-0.2, 0.1, 0.0));
    let build = || {
        let (mut bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
        if let Some(body) = bodies.get_mut(handle_b) {
            body.velocity = Vec3::new(0.5, -0.3, 0.0);
            body.angular_velocity = Vec3::new(0.0, 0.0, 1.0);
        }
        let spring = SpringJoint::new(handle_a, handle_b, anchors.0, anchors.1, 2.0, 80.0, 5.0);
        (bodies, handle_a, handle_b, spring)
    };

    let (mut bodies_once, ha1, hb1, mut spring_once) = build();
    spring_once.prepare(&bodies_once, DT);
    spring_once.solve_velocity(&mut bodies_once, DT);
    let (mut bodies_many, ha2, hb2, mut spring_many) = build();
    spring_many.prepare(&bodies_many, DT);
    for _ in 0..20 {
        spring_many.solve_velocity(&mut bodies_many, DT);
    }
    for (single, many) in [(ha1, ha2), (hb1, hb2)] {
        let one = bodies_once.get(single).unwrap();
        let twenty = bodies_many.get(many).unwrap();
        assert!((one.velocity - twenty.velocity).norm() < 1e-9);
        assert!((one.angular_velocity - twenty.angular_velocity).norm() < 1e-9);
    }
}

#[test]
fn spring_through_the_pgs_solver_matches_for_one_and_twenty_iterations() {
    let run = |iterations: usize| {
        let (mut bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
        if let Some(body) = bodies.get_mut(handle_b) {
            body.velocity = Vec3::new(1.0, 0.0, 0.0);
        }
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(SpringJoint::new(
            handle_a,
            handle_b,
            Vec3::zeros(),
            Vec3::zeros(),
            2.0,
            120.0,
            3.0,
        ))];
        PgsSolver::new(iterations, 0).solve(&mut constraints, &mut bodies, DT);
        (
            bodies.get(handle_a).unwrap().velocity,
            bodies.get(handle_b).unwrap().velocity,
        )
    };
    let (a1, b1) = run(1);
    let (a20, b20) = run(20);
    assert!((a1 - a20).norm() < 1e-6, "{a1:?} vs {a20:?}");
    assert!((b1 - b20).norm() < 1e-6, "{b1:?} vs {b20:?}");
    assert!(a1.x > 0.0 && b1.x < 1.0, "the spring must act at all");
}

#[test]
fn hanging_spring_settles_at_hookes_law_whatever_the_iteration_count() {
    let (stiffness, damping, rest_length) = (500.0, 20.0, 1.0);
    let expected_extension = GRAVITY / stiffness;
    let dt = 1.0 / 120.0;
    let mut finals = Vec::new();
    for iterations in [1_usize, 4, 16] {
        let (mut bodies, handle_a, handle_b) = hanging_bob(rest_length);
        let mut spring = SpringJoint::new(
            handle_a,
            handle_b,
            Vec3::zeros(),
            Vec3::zeros(),
            rest_length,
            stiffness,
            damping,
        );
        for _ in 0..2400 {
            step_joint(&mut spring, &mut bodies, &gravity(), dt, iterations, 0);
        }
        let bob = bodies.get(handle_b).unwrap();
        let extension = -bob.transform.position.y - rest_length;
        assert!(
            (extension - expected_extension).abs() < 1e-6,
            "{iterations} iterations: extension {extension}, Hooke {expected_extension}"
        );
        assert!(bob.velocity.norm() < 1e-6);
        finals.push(bob.transform.position.y);
    }
    assert!((finals[0] - finals[1]).abs() < 1e-9 && (finals[1] - finals[2]).abs() < 1e-9);
}

#[test]
fn spring_without_stiffness_and_damping_is_inert() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
    let mut spring = SpringJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        1.0,
        0.0,
        0.0,
    );
    spring.prepare(&bodies, DT);
    spring.solve_velocity(&mut bodies, DT);
    assert_eq!(bodies.get(handle_a).unwrap().velocity, Vec3::zeros());
    assert_eq!(bodies.get(handle_b).unwrap().velocity, Vec3::zeros());
    assert_eq!(spring.residual(), 0.0);
}

#[test]
fn relaxed_spring_updates_converge_to_the_unrelaxed_solution() {
    let build = || {
        let (bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
        let spring = SpringJoint::new(
            handle_a,
            handle_b,
            Vec3::zeros(),
            Vec3::zeros(),
            2.0,
            100.0,
            0.0,
        );
        (bodies, handle_a, spring)
    };
    let (mut bodies, handle_a, mut spring) = build();
    spring.prepare(&bodies, DT);
    spring.solve_velocity(&mut bodies, DT);
    let exact = bodies.get(handle_a).unwrap().velocity.x;

    let (mut bodies, handle_a, mut spring) = build();
    spring.prepare(&bodies, DT);
    spring.solve_velocity_relaxed(&mut bodies, DT, 0.5);
    spring.solve_velocity_relaxed(&mut bodies, DT, 0.5);
    // Each half-step closes half of the remaining gap: 1 - 0.5^2 = 0.75.
    let after_two = bodies.get(handle_a).unwrap().velocity.x;
    assert!(
        (after_two - 0.75 * exact).abs() < 1e-12,
        "{after_two} vs {}",
        0.75 * exact
    );
    for _ in 0..60 {
        spring.solve_velocity_relaxed(&mut bodies, DT, 0.5);
    }
    assert!((bodies.get(handle_a).unwrap().velocity.x - exact).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// F0380: BreakableJoint accounts for the impulse of one step only
// ---------------------------------------------------------------------------

/// Per-step impulse a hanging 1 kg mass demands: `m·g·dt`.
fn hanging_load(dt: Real) -> Real {
    GRAVITY * dt
}

fn hanging_breakable(threshold: Real) -> (RigidBodySet, BreakableJoint) {
    let (stiffness, damping, natural_length) = (5000.0, 10.0, 2.0);
    let drop = natural_length + GRAVITY / stiffness;
    let (bodies, handle_a, handle_b) = hanging_bob(drop);
    let joint = BreakableJoint::new(
        handle_a,
        handle_b,
        natural_length,
        stiffness,
        damping,
        threshold,
    );
    (bodies, joint)
}

#[test]
fn breakable_joint_holding_a_load_below_the_threshold_never_breaks() {
    let dt = 1.0 / 120.0;
    let load = hanging_load(dt);
    let (mut bodies, mut joint) = hanging_breakable(1.5 * load);
    let mut largest = 0.0_f64;
    for step in 0..10_000 {
        step_joint(&mut joint, &mut bodies, &gravity(), dt, 8, 0);
        assert!(!joint.is_broken(), "broke at step {step}");
        largest = largest.max(joint.current_force());
    }
    assert!(
        (largest - load).abs() < 1e-6,
        "the per-step impulse is m·g·dt = {load}, got {largest}"
    );
}

#[test]
fn breakable_joint_breaks_on_the_first_step_above_the_threshold() {
    let dt = 1.0 / 120.0;
    let load = hanging_load(dt);
    let (mut bodies, mut joint) = hanging_breakable(0.5 * load);
    step_joint(&mut joint, &mut bodies, &gravity(), dt, 8, 0);
    assert!(joint.is_broken(), "the load {load} exceeds the threshold");
    assert!(!joint.is_active(), "solvers must skip a broken joint");
}

#[test]
fn broken_breakable_joint_applies_no_further_impulse() {
    let dt = 1.0 / 120.0;
    let (mut bodies, mut joint) = hanging_breakable(0.5 * hanging_load(dt));
    step_joint(&mut joint, &mut bodies, &gravity(), dt, 8, 0);
    assert!(joint.is_broken());
    let handle_b = joint.body_handles()[1];
    let speed_after_break = bodies.get(handle_b).unwrap().velocity.y;
    for _ in 0..10 {
        step_joint(&mut joint, &mut bodies, &gravity(), dt, 8, 0);
    }
    // Free fall: every further step adds exactly -g·dt.
    let expected = speed_after_break - 10.0 * GRAVITY * dt;
    let speed = bodies.get(handle_b).unwrap().velocity.y;
    assert!(
        (speed - expected).abs() < 1e-9,
        "free fall expected {expected}, got {speed}"
    );
    assert_eq!(joint.current_force(), 0.0);
}

#[test]
fn breakable_joint_forgets_the_previous_step_and_can_be_reset() {
    let dt = 1.0 / 120.0;
    let load = hanging_load(dt);
    let (mut bodies, mut joint) = hanging_breakable(0.5 * load);
    step_joint(&mut joint, &mut bodies, &gravity(), dt, 8, 0);
    assert!(joint.is_broken());
    joint.reset();
    assert!(!joint.is_broken() && joint.is_active());
    assert_eq!(joint.current_force(), 0.0);
    // The impulse is measured per step: a fresh prepare clears it.
    joint.broken = false;
    joint.accumulated_force = 123.0;
    joint.prepare(&bodies, dt);
    assert_eq!(joint.current_force(), 0.0);
}

#[test]
fn breakable_joint_counts_compression_like_tension() {
    let (mut bodies, handle_a, handle_b) = hanging_bob(0.5);
    let mut joint = BreakableJoint::new(handle_a, handle_b, 1.0, 400.0, 0.0, 0.05);
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    assert!(
        joint.current_force() > 0.05 && joint.is_broken(),
        "a compressed spring pushes with |impulse| {}",
        joint.current_force()
    );
}

// ---------------------------------------------------------------------------
// F0389: XpbdJoint velocity pass is a compliant row that leaves lambda alone
// ---------------------------------------------------------------------------

#[test]
fn rigid_xpbd_joint_with_a_small_stretch_never_gains_kinetic_energy() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.01);
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        2.0,
        0.0,
        0.0,
    ))];
    let solver = PgsSolver::new(8, 4);
    let mut peak_energy = 0.0_f64;
    for _ in 0..100 {
        solver.solve(&mut constraints, &mut bodies, DT);
        for (_, body) in bodies.iter_mut() {
            body.integrate_velocity(DT);
        }
        peak_energy = peak_energy.max(total_kinetic_energy(&bodies));
    }
    assert!(
        peak_energy < 1e-20,
        "kinetic energy appeared: {peak_energy}"
    );
    let distance = (bodies.get(handle_b).unwrap().transform.position
        - bodies.get(handle_a).unwrap().transform.position)
        .norm();
    assert!(
        (distance - 2.0).abs() < 1e-9,
        "the position pass must restore the rest length, got {distance}"
    );
}

#[test]
fn xpbd_velocity_pass_does_not_touch_the_position_multiplier() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.01);
    if let Some(body) = bodies.get_mut(handle_b) {
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
    }
    let mut joint = XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        2.0,
        1e-4,
        0.0,
    );
    joint.prepare(&bodies, DT);
    for _ in 0..8 {
        joint.solve_velocity(&mut bodies, DT);
    }
    assert_eq!(joint.lambda, 0.0, "lambda belongs to solve_position alone");
    joint.solve_position(&mut bodies, DT);
    assert!(joint.lambda != 0.0, "the position pass accumulates lambda");
}

#[test]
fn rigid_xpbd_velocity_pass_cancels_the_relative_speed_along_the_joint() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.0);
    if let Some(body) = bodies.get_mut(handle_b) {
        body.velocity = Vec3::new(1.5, 0.7, 0.0);
    }
    let momentum = linear_momentum(&bodies);
    let mut joint = XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        2.0,
        0.0,
        0.0,
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    let velocity_a = bodies.get(handle_a).unwrap().velocity;
    let velocity_b = bodies.get(handle_b).unwrap().velocity;
    assert!((velocity_a.x - 0.75).abs() < 1e-12, "{velocity_a:?}");
    assert!((velocity_b.x - 0.75).abs() < 1e-12, "{velocity_b:?}");
    assert!(
        (velocity_b.y - 0.7).abs() < 1e-12 && velocity_a.y.abs() < 1e-12,
        "motion across the joint is not constrained"
    );
    assert!((linear_momentum(&bodies) - momentum).norm() < 1e-12);
}

#[test]
fn compliant_xpbd_velocity_pass_is_regularised_by_alpha_over_dt_squared() {
    // Relative speed u = 1 along the joint, K = 2, γ = α/dt² = 3.6:
    // the row removes 1/(K+γ) of it per unit, leaving γ/(K+γ) = 3.6/5.6.
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.0);
    if let Some(body) = bodies.get_mut(handle_b) {
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
    }
    let compliance = 1e-3;
    let mut joint = XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        2.0,
        compliance,
        0.0,
    );
    joint.prepare(&bodies, DT);
    for _ in 0..25 {
        joint.solve_velocity(&mut bodies, DT);
    }
    let relative =
        bodies.get(handle_b).unwrap().velocity.x - bodies.get(handle_a).unwrap().velocity.x;
    let gamma = compliance / (DT * DT);
    assert!(
        (relative - gamma / (2.0 + gamma)).abs() < 1e-9,
        "relative speed {relative}, expected {}",
        gamma / (2.0 + gamma)
    );
}

#[test]
fn xpbd_damping_adds_to_the_velocity_regularisation() {
    // k = 1/α = 1000 and c = 10 give γ = 1 / (dt (c + dt k)).
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.0);
    if let Some(body) = bodies.get_mut(handle_b) {
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
    }
    let mut joint = XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        2.0,
        1e-3,
        10.0,
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    let relative =
        bodies.get(handle_b).unwrap().velocity.x - bodies.get(handle_a).unwrap().velocity.x;
    let gamma = 1.0 / (DT * (10.0 + DT * 1000.0));
    assert!(
        (relative - gamma / (2.0 + gamma)).abs() < 1e-9,
        "relative speed {relative}, expected {}",
        gamma / (2.0 + gamma)
    );
}

#[test]
fn xpbd_velocity_pass_is_independent_of_the_iteration_count() {
    let run = |iterations: usize| {
        let (mut bodies, handle_a, handle_b) = pair_on_x_axis(2.0);
        if let Some(body) = bodies.get_mut(handle_b) {
            body.velocity = Vec3::new(2.0, 0.0, 0.0);
            body.angular_velocity = Vec3::new(0.0, 0.0, 3.0);
        }
        let mut joint = XpbdJoint::new(
            handle_a,
            handle_b,
            Vec3::new(0.0, 0.2, 0.0),
            Vec3::new(0.0, -0.2, 0.0),
            1.6,
            5e-4,
            1.0,
        );
        joint.prepare(&bodies, DT);
        for _ in 0..iterations {
            joint.solve_velocity(&mut bodies, DT);
        }
        (
            bodies.get(handle_a).unwrap().velocity,
            bodies.get(handle_b).unwrap().velocity,
            bodies.get(handle_b).unwrap().angular_velocity,
        )
    };
    let (a1, b1, w1) = run(1);
    let (a9, b9, w9) = run(9);
    assert!((a1 - a9).norm() < 1e-9 && (b1 - b9).norm() < 1e-9 && (w1 - w9).norm() < 1e-9);
}

#[test]
fn xpbd_joint_without_bodies_or_distance_stays_inert() {
    let (mut bodies, handle_a, handle_b) = pair_on_x_axis(0.0);
    let mut joint = XpbdJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        1.0,
        0.0,
        0.0,
    );
    joint.prepare(&bodies, DT);
    joint.solve_velocity(&mut bodies, DT);
    joint.solve_position(&mut bodies, DT);
    assert_eq!(bodies.get(handle_a).unwrap().velocity, Vec3::zeros());
    assert_eq!(
        bodies.get(handle_b).unwrap().transform.position,
        Vec3::zeros()
    );
}

// ---------------------------------------------------------------------------
// Physical checks beyond one dimension
// ---------------------------------------------------------------------------

/// Period of a mathematical pendulum of length `length` swinging through the
/// amplitude `amplitude` (radians), from the series of the elliptic integral.
fn rigid_pendulum_period(length: Real, amplitude: Real) -> Real {
    let small_angle = std::f64::consts::TAU * (length / GRAVITY).sqrt();
    let squared = amplitude * amplitude;
    small_angle * (1.0 + squared / 16.0 + 11.0 * squared * squared / 3072.0)
}

/// Time between two successive downward zero crossings of the x coordinate of a
/// 1 kg bob that swings from `amplitude` radians on a one metre rod made of the
/// spring `(stiffness, damping)`, driven with `iterations` solver iterations.
fn spring_pendulum_period(
    stiffness: Real,
    damping: Real,
    amplitude: Real,
    iterations: usize,
) -> Real {
    let length = 1.0;
    let dt = 1.0 / 2000.0;
    let mut bodies = RigidBodySet::new();
    let handle_a = bodies.insert(RigidBody::new_static());
    let mut bob = RigidBody::new(1.0);
    bob.transform.position = Vec3::new(length * amplitude.sin(), -length * amplitude.cos(), 0.0);
    let handle_b = bodies.insert(bob);
    let mut spring = SpringJoint::new(
        handle_a,
        handle_b,
        Vec3::zeros(),
        Vec3::zeros(),
        length,
        stiffness,
        damping,
    );
    let bob_x =
        |bodies: &RigidBodySet| bodies.get(handle_b).map_or(0.0, |b| b.transform.position.x);
    let mut crossings = Vec::new();
    let mut previous_x = bob_x(&bodies);
    for step in 1..=12_000_u32 {
        step_joint(&mut spring, &mut bodies, &gravity(), dt, iterations, 0);
        let x = bob_x(&bodies);
        if previous_x > 0.0 && x <= 0.0 {
            // Linear interpolation of the crossing time inside the step.
            let fraction = previous_x / (previous_x - x);
            crossings.push((Real::from(step) - 1.0 + fraction) * dt);
        }
        previous_x = x;
    }
    assert!(
        crossings.len() >= 2,
        "the pendulum must complete a period, crossings {crossings:?}"
    );
    crossings[1] - crossings[0]
}

#[test]
fn stiff_spring_swings_like_a_rigid_pendulum_rod() {
    // The implicit spring row keeps the bob on the rod while it swings, so the
    // period is that of the mathematical pendulum up to the tiny sag of the
    // rod (the effective length is l + m g / k) and the O(dt) numerical error.
    let amplitude = 0.1;
    let expected = rigid_pendulum_period(1.0, amplitude);
    for (stiffness, damping) in [(50_000.0, 100.0), (50_000.0, 0.0), (200_000.0, 400.0)] {
        let period = spring_pendulum_period(stiffness, damping, amplitude, 1);
        assert!(
            (period - expected).abs() < 1e-3 * expected,
            "k={stiffness} c={damping}: period {period}, rigid pendulum {expected}"
        );
    }
}

#[test]
fn spring_pendulum_period_does_not_depend_on_the_iteration_count() {
    let amplitude = 0.1;
    let single = spring_pendulum_period(50_000.0, 100.0, amplitude, 1);
    let many = spring_pendulum_period(50_000.0, 100.0, amplitude, 8);
    assert!(
        (single - many).abs() < 1e-9,
        "1 iteration gave {single}, 8 iterations gave {many}"
    );
}

#[test]
fn over_relaxed_spring_updates_converge_to_the_unrelaxed_solution() {
    let build = || {
        let (bodies, handle_a, handle_b) = pair_on_x_axis(3.0);
        let spring = SpringJoint::new(
            handle_a,
            handle_b,
            Vec3::zeros(),
            Vec3::zeros(),
            2.0,
            100.0,
            0.0,
        );
        (bodies, handle_a, spring)
    };
    let (mut bodies, handle_a, mut spring) = build();
    spring.prepare(&bodies, DT);
    spring.solve_velocity(&mut bodies, DT);
    let exact = bodies.get(handle_a).unwrap().velocity.x;

    // A single-row problem contracts its error by |1 - omega| per update, so
    // omega = 1.5 overshoots to 1.5x, returns to 0.75x and then settles.
    let (mut bodies, handle_a, mut spring) = build();
    spring.prepare(&bodies, DT);
    spring.solve_velocity_relaxed(&mut bodies, DT, 1.5);
    let first = bodies.get(handle_a).unwrap().velocity.x;
    assert!(
        (first - 1.5 * exact).abs() < 1e-12,
        "{first} vs {}",
        1.5 * exact
    );
    spring.solve_velocity_relaxed(&mut bodies, DT, 1.5);
    let second = bodies.get(handle_a).unwrap().velocity.x;
    assert!(
        (second - 0.75 * exact).abs() < 1e-12,
        "{second} vs {}",
        0.75 * exact
    );
    for _ in 0..60 {
        spring.solve_velocity_relaxed(&mut bodies, DT, 1.5);
    }
    assert!((bodies.get(handle_a).unwrap().velocity.x - exact).abs() < 1e-9);
}

#[test]
fn breakable_joint_through_the_pgs_solver_breaks_and_is_then_skipped() {
    let dt = 1.0 / 120.0;
    let (mut bodies, joint) = hanging_breakable(0.5 * hanging_load(dt));
    let handle_b = joint.body_handles()[1];
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(joint)];
    let solver = PgsSolver::new(8, 0);
    for (_, body) in bodies.iter_mut() {
        body.integrate_forces(dt, &gravity());
    }
    solver.solve(&mut constraints, &mut bodies, dt);
    assert!(
        !constraints[0].is_active(),
        "the load exceeds the threshold: the solver must retire the joint"
    );
    // From now on the solver neither prepares nor solves the joint: free fall.
    let speed_after_break = bodies.get(handle_b).unwrap().velocity.y;
    for (_, body) in bodies.iter_mut() {
        body.integrate_forces(dt, &gravity());
    }
    solver.solve(&mut constraints, &mut bodies, dt);
    let speed = bodies.get(handle_b).unwrap().velocity.y;
    assert!(
        (speed - (speed_after_break - GRAVITY * dt)).abs() < 1e-12,
        "free fall expected, got {speed} after {speed_after_break}"
    );
}
