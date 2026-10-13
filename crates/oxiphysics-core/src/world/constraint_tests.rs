// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the XPBD constraint pass of [`PhysicsWorld::step`].
//!
//! Before the pass existed `add_constraint` stored its argument and the step
//! never looked at it, so pendulums fell, ropes stretched and joints came
//! apart. Every test here would fail on that behaviour.

use super::{Body, BodyMode, Constraint, ConstraintKind, PhysicsEvent, PhysicsWorld};
use crate::math::{Quat, Vec3};
use crate::types::{BodyHandle, TimeStep};

const DT: f64 = 1.0 / 60.0;

/// A dynamic body without velocity damping, so energy and momentum tests are clean.
fn undamped_body(position: [f64; 3], mass: f64) -> Body {
    let mut body = Body::new_dynamic(position, mass);
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body
}

fn zero_gravity_world() -> PhysicsWorld {
    let mut world = PhysicsWorld::new();
    world.config.gravity = Vec3::zeros();
    world
}

fn constraint(
    kind: ConstraintKind,
    body_a: BodyHandle,
    body_b: Option<BodyHandle>,
    anchors: ([f64; 3], [f64; 3]),
    target: f64,
    compliance: f64,
) -> Constraint {
    Constraint {
        body_a,
        body_b,
        kind,
        anchor_a: anchors.0,
        anchor_b: anchors.1,
        target,
        compliance,
        active: true,
    }
}

fn position(world: &PhysicsWorld, handle: BodyHandle) -> Vec3 {
    world
        .body(handle)
        .map(|b| b.transform.position)
        .unwrap_or_else(Vec3::zeros)
}

/// World-space position of a local anchor of `handle`.
fn world_anchor(world: &PhysicsWorld, handle: BodyHandle, local: [f64; 3]) -> Vec3 {
    world
        .body(handle)
        .map(|b| {
            b.transform
                .transform_point(&Vec3::new(local[0], local[1], local[2]))
        })
        .unwrap_or_else(Vec3::zeros)
}

/// A static pivot at the origin and an undamped unit-mass bob at `bob_position`,
/// optionally joined by a constraint `(kind, target, compliance)` between the
/// pivot origin and the bob centre.
fn pivot_and_bob(
    bob_position: [f64; 3],
    joint: Option<(ConstraintKind, f64, f64)>,
) -> (PhysicsWorld, BodyHandle) {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0; 3]));
    let bob = world.add_body(undamped_body(bob_position, 1.0));
    if let Some((kind, target, compliance)) = joint {
        world.add_constraint(constraint(
            kind,
            pivot,
            Some(bob),
            ([0.0; 3], [0.0; 3]),
            target,
            compliance,
        ));
    }
    (world, bob)
}

#[test]
fn distance_constraint_holds_a_swinging_pendulum_to_its_length() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0, 0.0, 0.0]));
    let bob = world.add_body(undamped_body([1.0, 0.0, 0.0], 1.0));
    world.add_constraint(Constraint::distance(pivot, bob, [0.0; 3], [0.0; 3], 1.0));

    let mut worst_length_error = 0.0_f64;
    let mut lowest = 0.0_f64;
    let mut leftmost = 0.0_f64;
    for _ in 0..1000 {
        world.step(&TimeStep::new(DT));
        let p = position(&world, bob);
        worst_length_error = worst_length_error.max((p.norm() - 1.0).abs());
        lowest = lowest.min(p.y);
        leftmost = leftmost.min(p.x);
    }
    // A rigid distance constraint to an immovable pivot is solved exactly by a
    // single Gauss-Seidel update, so the error does not accumulate over time.
    assert!(
        worst_length_error < 1e-9,
        "pendulum length error {worst_length_error}"
    );
    // Released from the horizontal it swings through the bottom of the circle
    // and up the other side instead of dropping like a free body.
    assert!(lowest > -1.0 - 1e-3, "lowest point {lowest}");
    assert!(
        lowest < -0.99,
        "the bob must reach the bottom, got {lowest}"
    );
    assert!(leftmost < -0.5, "the bob must swing across, got {leftmost}");
}

#[test]
fn velocity_after_the_solve_is_tangential_to_the_pendulum_circle() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0, 0.0, 0.0]));
    let bob = world.add_body(undamped_body([0.0, -1.0, 0.0], 1.0));
    world.add_constraint(Constraint::distance(pivot, bob, [0.0; 3], [0.0; 3], 1.0));
    world.set_velocity(bob, [2.0, 0.0, 0.0], [0.0; 3]);
    for _ in 0..30 {
        world.step(&TimeStep::new(DT));
        let p = position(&world, bob);
        let v = world
            .body(bob)
            .map(|b| {
                Vec3::new(
                    b.velocity.linear[0],
                    b.velocity.linear[1],
                    b.velocity.linear[2],
                )
            })
            .unwrap_or_else(Vec3::zeros);
        let radial = v.dot(&p) / p.norm();
        let speed = v.norm();
        assert!(
            radial.abs() < 0.05 * speed.max(1.0),
            "radial velocity {radial} of a bob moving at {speed}"
        );
    }
}

#[test]
fn offset_ball_socket_swings_a_physical_pendulum_with_the_right_period() {
    // A body of mass 1 and inertia 1 whose pivot is 1 m above its centre:
    // I_pivot = 1 + 1 * 1^2 = 2 and T = 2 pi sqrt(I_pivot / (m g L)).
    let theta0 = 0.05_f64;
    let dt = 1.0 / 480.0;
    let mut world = PhysicsWorld::new();
    let mut body = undamped_body([theta0.sin(), -theta0.cos(), 0.0], 1.0);
    body.inertia = 1.0;
    body.transform.rotation = Quat::from_axis_angle(&Vec3::z_axis(), theta0);
    let handle = world.add_body(body);
    world.add_constraint(constraint(
        ConstraintKind::BallSocket,
        handle,
        None,
        ([0.0, 1.0, 0.0], [0.0; 3]),
        0.0,
        0.0,
    ));

    let angle = |world: &PhysicsWorld| {
        world
            .body(handle)
            .map(|b| b.transform.rotation.scaled_axis().z)
            .unwrap_or(0.0)
    };
    let mut crossings: Vec<f64> = Vec::new();
    let mut previous_angle = angle(&world);
    let mut worst_pivot_error = 0.0_f64;
    for step in 1..=(6.0 / dt) as usize {
        world.step(&TimeStep::new(dt));
        let a = angle(&world);
        if previous_angle.signum() != a.signum() {
            let fraction = previous_angle / (previous_angle - a);
            crossings.push((step as f64 - 1.0 + fraction) * dt);
        }
        previous_angle = a;
        worst_pivot_error =
            worst_pivot_error.max(world_anchor(&world, handle, [0.0, 1.0, 0.0]).norm());
    }
    assert!(worst_pivot_error < 1e-4, "pivot drift {worst_pivot_error}");
    assert!(crossings.len() >= 3, "crossings: {crossings:?}");
    // The angle starts at +theta0 and crosses zero every half period, so two
    // crossings in the same direction are one full period apart.
    let measured_period = crossings[2] - crossings[0];
    let expected_period = 2.0 * std::f64::consts::PI * (2.0_f64 / 9.81).sqrt();
    assert!(
        (measured_period - expected_period).abs() < 0.03 * expected_period,
        "period {measured_period}, expected {expected_period}"
    );
}

#[test]
fn two_body_ball_socket_conserves_momentum_and_keeps_the_anchors_together() {
    let mut world = zero_gravity_world();
    let mut light = undamped_body([0.0, 0.0, 0.0], 1.0);
    light.inertia = 0.5;
    light.velocity.linear = [0.0, 1.0, 0.0];
    let mut heavy = undamped_body([2.0, 0.0, 0.0], 3.0);
    heavy.inertia = 2.0;
    let a = world.add_body(light);
    let b = world.add_body(heavy);
    world.add_constraint(Constraint::ball_socket(
        a,
        b,
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
    ));

    let initial_momentum = world.total_momentum();
    let initial_energy = world.total_kinetic_energy();
    let mut worst_gap = 0.0_f64;
    for _ in 0..240 {
        world.step(&TimeStep::new(DT));
        let gap = (world_anchor(&world, a, [1.0, 0.0, 0.0])
            - world_anchor(&world, b, [-1.0, 0.0, 0.0]))
        .norm();
        worst_gap = worst_gap.max(gap);
    }
    assert!(worst_gap < 1e-4, "anchors drifted apart by {worst_gap}");
    let momentum = world.total_momentum();
    for k in 0..3 {
        assert!(
            (momentum[k] - initial_momentum[k]).abs() < 1e-9,
            "momentum component {k}: {} vs {}",
            momentum[k],
            initial_momentum[k]
        );
    }
    // Position-based solving only ever removes energy.
    assert!(world.total_kinetic_energy() <= initial_energy + 1e-9);
    // The off-centre anchor really couples translation and rotation: both
    // bodies spin up even though nobody was given an angular velocity.
    let spin_a = world.body(a).map(|body| body.velocity.angular[2]);
    let spin_b = world.body(b).map(|body| body.velocity.angular[2]);
    assert!(
        spin_a.is_some_and(|w| w.abs() > 1e-3),
        "spin of a: {spin_a:?}"
    );
    assert!(
        spin_b.is_some_and(|w| w.abs() > 1e-3),
        "spin of b: {spin_b:?}"
    );
}

#[test]
fn ball_socket_to_the_world_pins_the_anchor_point() {
    let mut world = PhysicsWorld::new();
    let handle = world.add_body(undamped_body([0.0, 3.0, 0.0], 2.0));
    world.add_constraint(constraint(
        ConstraintKind::BallSocket,
        handle,
        None,
        ([0.0; 3], [5.0, 3.0, 0.0]),
        0.0,
        0.0,
    ));
    for _ in 0..120 {
        world.step(&TimeStep::new(DT));
    }
    let p = position(&world, handle);
    assert!(
        (p - Vec3::new(5.0, 3.0, 0.0)).norm() < 1e-4,
        "body should sit on the world anchor, got {p:?}"
    );
}

#[test]
fn max_distance_is_a_rope_that_only_pulls() {
    let (mut roped, roped_bob) = pivot_and_bob(
        [0.5, 0.0, 0.0],
        Some((ConstraintKind::DistanceMax, 1.0, 0.0)),
    );
    let (mut free, free_bob) = pivot_and_bob([0.5, 0.0, 0.0], None);

    let mut slack_steps = 0;
    let mut taut_steps = 0;
    for _ in 0..240 {
        roped.step(&TimeStep::new(DT));
        free.step(&TimeStep::new(DT));
        let p_roped = position(&roped, roped_bob);
        let p_free = position(&free, free_bob);
        if p_free.norm() < 0.999 {
            // Inside the rope length the constraint must not interfere at all.
            slack_steps += 1;
            assert!(
                (p_roped - p_free).norm() < 1e-12,
                "slack rope changed the motion: {p_roped:?} vs {p_free:?}"
            );
        } else {
            taut_steps += 1;
        }
        assert!(
            p_roped.norm() < 1.0 + 1e-3,
            "rope stretched to {}",
            p_roped.norm()
        );
    }
    assert!(slack_steps > 5, "the bob should fall freely for a while");
    assert!(taut_steps > 5, "the free body must leave the rope length");
}

#[test]
fn min_distance_pushes_apart_but_never_pulls_together() {
    let mut world = zero_gravity_world();
    let anchor = world.add_body(Body::new_static([0.0; 3]));
    let too_close = world.add_body(undamped_body([0.5, 0.0, 0.0], 1.0));
    let far_away = world.add_body(undamped_body([0.0, 2.0, 0.0], 1.0));
    world.set_velocity(far_away, [0.3, 0.0, 0.0], [0.0; 3]);
    for body in [too_close, far_away] {
        world.add_constraint(constraint(
            ConstraintKind::DistanceMin,
            anchor,
            Some(body),
            ([0.0; 3], [0.0; 3]),
            1.0,
            0.0,
        ));
    }
    world.step(&TimeStep::new(DT));
    let pushed = position(&world, too_close);
    assert!(
        (pushed - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-9,
        "the close body is pushed out to the minimum distance, got {pushed:?}"
    );
    let untouched = position(&world, far_away);
    let expected = Vec3::new(0.3 * DT, 2.0, 0.0);
    assert!(
        (untouched - expected).norm() < 1e-12,
        "a body beyond the minimum distance must not be pulled in: {untouched:?}"
    );
    let v = world.body(far_away).map(|b| b.velocity.linear);
    assert_eq!(v, Some([0.3, 0.0, 0.0]));
}

#[test]
fn fixed_distance_with_coincident_anchors_still_separates_the_bodies() {
    let mut world = zero_gravity_world();
    let a = world.add_body(undamped_body([1.0, 1.0, 1.0], 1.0));
    let b = world.add_body(undamped_body([1.0, 1.0, 1.0], 1.0));
    world.add_constraint(Constraint::distance(a, b, [0.0; 3], [0.0; 3], 1.0));
    world.step(&TimeStep::new(DT));
    let separation = (position(&world, a) - position(&world, b)).norm();
    assert!(
        (separation - 1.0).abs() < 1e-9,
        "bodies should be pushed to the target distance, got {separation}"
    );
}

#[test]
fn compliance_makes_the_constraint_softer_in_proportion_to_the_load() {
    // A hanging mass m on a spring of compliance alpha settles at an
    // extension m * g * alpha, independent of the time step. (The XPBD update is
    // implicit, so the undamped oscillation dies out by itself.)
    let stretch = |compliance: f64| {
        let mut world = PhysicsWorld::new();
        // Keep the slowly settling bob awake: a sleeping body is frozen where it
        // is, before the tail of the decaying oscillation has died out.
        world.config.time_before_sleep = f64::MAX;
        let pivot = world.add_body(Body::new_static([0.0; 3]));
        let bob = world.add_body(undamped_body([0.0, -1.0, 0.0], 1.0));
        world.add_constraint(constraint(
            ConstraintKind::DistanceFixed,
            pivot,
            Some(bob),
            ([0.0; 3], [0.0; 3]),
            1.0,
            compliance,
        ));
        for _ in 0..600 {
            world.step(&TimeStep::new(DT));
        }
        position(&world, bob).norm() - 1.0
    };
    let rigid = stretch(0.0);
    let soft = stretch(1e-3);
    let softer = stretch(1e-2);
    assert!(rigid.abs() < 1e-6, "rigid stretch {rigid}");
    assert!(
        (soft - 9.81e-3).abs() < 1e-5,
        "stretch {soft}, expected m g alpha = 9.81e-3"
    );
    assert!(
        (softer - 9.81e-2).abs() < 1e-4,
        "stretch {softer}, expected m g alpha = 9.81e-2"
    );
    assert!(softer > soft && soft > rigid);
}

#[test]
fn a_hanging_pendulum_falls_asleep_because_sleep_sees_the_solved_velocity() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0; 3]));
    let bob = world.add_body(Body::new_dynamic([0.0, -1.0, 0.0], 1.0));
    world.add_constraint(Constraint::distance(pivot, bob, [0.0; 3], [0.0; 3], 1.0));
    for _ in 0..90 {
        world.step(&TimeStep::new(DT));
    }
    assert_eq!(
        world.body(bob).map(|b| b.mode),
        Some(BodyMode::Sleeping),
        "gravity is exactly cancelled by the rope, so the bob is at rest"
    );
    let slept = world
        .drain_events()
        .into_iter()
        .any(|e| matches!(e, PhysicsEvent::BodySlept { body } if body == bob));
    assert!(slept, "a BodySlept event must have been emitted");
    let p = position(&world, bob);
    assert!((p - Vec3::new(0.0, -1.0, 0.0)).norm() < 1e-6, "{p:?}");
}

#[test]
fn inactive_constraints_do_nothing() {
    let joint = Some((ConstraintKind::DistanceFixed, 1.0, 0.0));
    let (mut constrained, constrained_bob) = pivot_and_bob([1.0, 0.0, 0.0], joint);
    constrained.constraints[0].active = false;
    let (mut free, free_bob) = pivot_and_bob([1.0, 0.0, 0.0], None);
    for _ in 0..60 {
        constrained.step(&TimeStep::new(DT));
        free.step(&TimeStep::new(DT));
    }
    assert_eq!(
        position(&constrained, constrained_bob),
        position(&free, free_bob)
    );
    // Re-activating the constraint makes the same world honour it again.
    constrained.constraints[0].active = true;
    for _ in 0..60 {
        constrained.step(&TimeStep::new(DT));
    }
    let length = position(&constrained, constrained_bob).norm();
    assert!(
        (length - 1.0).abs() < 1e-3,
        "length after reactivation {length}"
    );
}

#[test]
fn hinge_and_contact_constraints_are_not_enforced_by_the_step() {
    assert!(ConstraintKind::BallSocket.is_solved_by_step());
    assert!(ConstraintKind::DistanceMin.is_solved_by_step());
    assert!(ConstraintKind::DistanceMax.is_solved_by_step());
    assert!(ConstraintKind::DistanceFixed.is_solved_by_step());
    assert!(!ConstraintKind::Hinge.is_solved_by_step());
    assert!(!ConstraintKind::Contact.is_solved_by_step());

    for kind in [ConstraintKind::Hinge, ConstraintKind::Contact] {
        let (mut constrained, constrained_bob) =
            pivot_and_bob([1.0, 0.0, 0.0], Some((kind, 1.0, 0.0)));
        let (mut free, free_bob) = pivot_and_bob([1.0, 0.0, 0.0], None);
        for _ in 0..30 {
            constrained.step(&TimeStep::new(DT));
            free.step(&TimeStep::new(DT));
        }
        assert_eq!(
            position(&constrained, constrained_bob),
            position(&free, free_bob),
            "{kind:?} must be left alone"
        );
    }
}

#[test]
fn sleeping_and_static_bodies_are_immovable_references() {
    let mut world = zero_gravity_world();
    let mut sleeper = undamped_body([0.0; 3], 1.0);
    sleeper.mode = BodyMode::Sleeping;
    let sleeper = world.add_body(sleeper);
    let mover = world.add_body(undamped_body([2.0, 0.0, 0.0], 1.0));
    let wall = world.add_body(Body::new_static([0.0, 5.0, 0.0]));
    let other = world.add_body(Body::new_static([0.0, 8.0, 0.0]));
    world.add_constraint(Constraint::distance(
        sleeper, mover, [0.0; 3], [0.0; 3], 1.0,
    ));
    world.add_constraint(Constraint::distance(wall, other, [0.0; 3], [0.0; 3], 1.0));
    world.step(&TimeStep::new(DT));
    assert_eq!(position(&world, sleeper), Vec3::zeros());
    assert!(
        (position(&world, mover) - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-9,
        "the awake body is pulled to the sleeping one"
    );
    assert_eq!(position(&world, wall), Vec3::new(0.0, 5.0, 0.0));
    assert_eq!(position(&world, other), Vec3::new(0.0, 8.0, 0.0));
}

#[test]
fn malformed_constraints_are_skipped_without_poisoning_the_world() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0; 3]));
    let bob = world.add_body(undamped_body([1.0, 0.0, 0.0], 1.0));
    let removed = world.add_body(undamped_body([3.0, 0.0, 0.0], 1.0));
    // Push constraints directly so that the stale handle survives `remove_body`.
    let stale_a = constraint(
        ConstraintKind::DistanceFixed,
        pivot,
        Some(removed),
        ([0.0; 3], [0.0; 3]),
        1.0,
        0.0,
    );
    let _ = world.remove_body(removed);
    world.constraints.push(stale_a);
    world.constraints.push(constraint(
        ConstraintKind::DistanceFixed,
        pivot,
        Some(bob),
        ([f64::NAN, 0.0, 0.0], [0.0; 3]),
        1.0,
        0.0,
    ));
    world.constraints.push(constraint(
        ConstraintKind::DistanceFixed,
        pivot,
        Some(bob),
        ([0.0; 3], [0.0; 3]),
        f64::INFINITY,
        0.0,
    ));
    world.constraints.push(constraint(
        ConstraintKind::DistanceFixed,
        pivot,
        Some(bob),
        ([0.0; 3], [0.0; 3]),
        1.0,
        f64::NAN,
    ));
    world.constraints.push(constraint(
        ConstraintKind::DistanceFixed,
        bob,
        Some(bob),
        ([0.0; 3], [0.0; 3]),
        1.0,
        0.0,
    ));
    for _ in 0..30 {
        world.step(&TimeStep::new(DT));
    }
    let p = position(&world, bob);
    assert!(
        p.iter().all(|c| c.is_finite()),
        "bad constraints must not poison the body: {p:?}"
    );
    // None of them constrained the bob, so it simply fell.
    assert!(p.y < -1.0, "the bob was free to fall, y = {}", p.y);
}

#[test]
fn removing_a_body_drops_its_constraints_before_the_next_step() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0; 3]));
    let bob = world.add_body(undamped_body([1.0, 0.0, 0.0], 1.0));
    world.add_constraint(Constraint::distance(pivot, bob, [0.0; 3], [0.0; 3], 1.0));
    assert!(world.remove_body(bob).is_some());
    assert!(world.constraints.is_empty());
    world.step(&TimeStep::new(DT));
    assert_eq!(world.body_count(), 1);
}

#[test]
fn accumulator_substeps_enforce_the_constraints_too() {
    let mut world = PhysicsWorld::new();
    let pivot = world.add_body(Body::new_static([0.0; 3]));
    let bob = world.add_body(undamped_body([1.0, 0.0, 0.0], 1.0));
    world.add_constraint(Constraint::distance(pivot, bob, [0.0; 3], [0.0; 3], 1.0));
    let mut substeps = 0;
    for _ in 0..90 {
        // Three fixed substeps of 1/60 s per call.
        substeps += world
            .step_with_accumulator(1.0 / 20.0)
            .expect("a finite, non-negative dt is valid input");
        let length = position(&world, bob).norm();
        assert!(
            (length - 1.0).abs() < 1e-9,
            "length {length} after {substeps} substeps"
        );
    }
    assert!(substeps >= 260, "only {substeps} substeps ran");
    assert!(
        position(&world, bob).y < -0.5,
        "the bob must have swung down"
    );
}

#[test]
fn every_extra_solver_iteration_tightens_a_chain_of_rigid_links() {
    // Five unit links hang from a static pivot, started horizontal. A single
    // Gauss-Seidel sweep fixes each link at the expense of its neighbours, so
    // the stretch only vanishes as `config.solver_iterations` grows.
    let worst_stretch = |iterations: u32| {
        let mut world = PhysicsWorld::new();
        world.config.solver_iterations = iterations;
        world.config.time_before_sleep = f64::MAX;
        let mut previous = world.add_body(Body::new_static([0.0; 3]));
        let mut links: Vec<(BodyHandle, BodyHandle)> = Vec::new();
        for k in 1..=5 {
            let bob = world.add_body(undamped_body([k as f64, 0.0, 0.0], 1.0));
            world.add_constraint(Constraint::distance(previous, bob, [0.0; 3], [0.0; 3], 1.0));
            links.push((previous, bob));
            previous = bob;
        }
        let mut worst = 0.0_f64;
        for _ in 0..120 {
            world.step(&TimeStep::new(DT));
            for &(upper, lower) in &links {
                let length = (position(&world, upper) - position(&world, lower)).norm();
                worst = worst.max((length - 1.0).abs());
            }
        }
        worst
    };
    let one = worst_stretch(1);
    let four = worst_stretch(4);
    let many = worst_stretch(40);
    assert!(
        one > four && four > many,
        "stretch must shrink with more sweeps: {one} / {four} / {many}"
    );
    assert!(many < 1e-3, "40 sweeps still leave a stretch of {many}");
}
