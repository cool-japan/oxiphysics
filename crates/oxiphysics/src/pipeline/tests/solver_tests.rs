// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the sequential-impulse contact solver as driven by the pipeline:
//! accumulated impulses, restitution, friction, warm starting, split-impulse
//! position correction, sleeping and determinism.

use super::*;
use oxiphysics_core::math::Quat;

/// Gravity of the default configuration.
const G: Real = 9.81;

/// Penetration the solver tolerates without correcting it.
const SLOP: Real = 0.005;

/// A unit-mass box on a wide static floor whose top face is the plane y = 0.
struct RestScene {
    scene: Scene,
    body: BodyHandle,
    box_collider: ColliderHandle,
    floor_collider: ColliderHandle,
}

/// Build the resting-box scene. `box_first` decides which of the two colliders
/// is collider A of the pair (the order flips the sign convention the solver
/// sees); `start_y` is the height of the box centre.
fn rest_scene(box_first: bool, start_y: Real, friction: Real, restitution: Real) -> RestScene {
    let mut scene = Scene::new();
    let add_box = |scene: &mut Scene| {
        scene.dynamic(
            1.0,
            Vec3::new(0.0, start_y, 0.0),
            cuboid(0.5, 0.5, 0.5),
            friction,
            restitution,
        )
    };
    let add_floor = |scene: &mut Scene| {
        scene.static_box(
            Vec3::new(0.0, -0.5, 0.0),
            Vec3::new(10.0, 0.5, 10.0),
            friction,
            restitution,
        )
    };
    let ((body, box_collider), (_, floor_collider)) = if box_first {
        let b = add_box(&mut scene);
        let f = add_floor(&mut scene);
        (b, f)
    } else {
        let f = add_floor(&mut scene);
        let b = add_box(&mut scene);
        (b, f)
    };
    RestScene {
        scene,
        body,
        box_collider,
        floor_collider,
    }
}

fn config_with_iterations(iterations: u32) -> PhysicsConfig {
    PhysicsConfig {
        solver_iterations: iterations,
        ..awake_config()
    }
}

#[test]
fn resting_box_stays_put_for_either_collider_order() {
    for box_first in [true, false] {
        let mut pipeline = PhysicsPipeline::with_config(awake_config());
        // 3 mm of penetration: inside the slop, so nothing should move.
        let mut rest = rest_scene(box_first, 0.497, 0.5, 0.0);
        let mut max_speed: Real = 0.0;
        for step in 0..120 {
            rest.scene.step(&mut pipeline);
            if step >= 5 {
                max_speed = max_speed.max(rest.scene.velocity(rest.body).norm());
            }
        }

        let body = rest.scene.bodies.get(rest.body).unwrap();
        let penetration = 0.5 - body.transform.position.y;
        assert!(
            max_speed < 1e-3,
            "box_first={box_first}: a resting box must not move, max speed {max_speed}"
        );
        assert!(
            (-1e-3..=SLOP + 1e-3).contains(&penetration),
            "box_first={box_first}: penetration {penetration}"
        );
        assert!(body.transform.position.x.abs() < 1e-6);
        assert!(body.transform.position.z.abs() < 1e-6);
        assert!(
            body.transform.rotation.angle() < 1e-4,
            "box_first={box_first}: the box rocked by {} rad",
            body.transform.rotation.angle()
        );
        assert!(body.angular_velocity.norm() < 1e-3);
        assert_eq!(pipeline.stats.last_contact_manifolds, 1);
    }
}

#[test]
fn resting_contact_carries_the_weight_of_the_box() {
    for box_first in [true, false] {
        let mut pipeline = PhysicsPipeline::with_config(awake_config());
        let mut rest = rest_scene(box_first, 0.497, 0.5, 0.0);
        for _ in 0..30 {
            rest.scene.step(&mut pipeline);
        }

        // The accumulated normal impulses over all contact points of the pair
        // support m * g * dt per step, regardless of the argument order.
        let cache = pipeline.warm_start_cache();
        let weight_impulse = G * DT;
        let total = cache.normal_impulse(rest.box_collider, rest.floor_collider);
        assert!(
            (total - weight_impulse).abs() < 0.03 * weight_impulse,
            "box_first={box_first}: supporting impulse {total}, weight impulse {weight_impulse}"
        );
        assert_eq!(
            total,
            cache.normal_impulse(rest.floor_collider, rest.box_collider)
        );

        // Every contact point is loaded, each by a fair share (four corners
        // under one box leave one degree of freedom of the load distribution
        // open, so the shares differ slightly depending on the solve order).
        let impulses = cache.contact_normal_impulses(rest.box_collider, rest.floor_collider);
        assert!(!impulses.is_empty());
        let share = total / impulses.len() as Real;
        for impulse in &impulses {
            assert!(
                *impulse > 0.0 && (impulse - share).abs() < 0.3 * share,
                "box_first={box_first}: contact impulses {impulses:?}"
            );
        }
        assert_eq!(cache.pair_count(), 1);
        assert_eq!(cache.contact_count(), impulses.len());
    }
}

#[test]
fn penetration_is_corrected_to_the_slop_without_adding_energy() {
    for box_first in [true, false] {
        let mut pipeline = PhysicsPipeline::with_config(awake_config());
        // 2 cm of penetration: 15 mm above the slop.
        let mut rest = rest_scene(box_first, 0.48, 0.5, 0.0);

        rest.scene.step(&mut pipeline);
        // Position correction works on pseudo velocities, so the first step
        // must not kick the box (a velocity bias would give ~0.18 m/s; what
        // is left here is the unconverged tail of the velocity iterations).
        let first_speed = rest.scene.velocity(rest.body).norm();
        assert!(
            first_speed < 1e-2,
            "box_first={box_first}: position correction added velocity {first_speed}"
        );
        // It removed 20 % of the excess: 0.2 * 15 mm = 3 mm.
        let y1 = rest.scene.position(rest.body).y;
        assert!(
            (y1 - 0.483).abs() < 2e-4,
            "box_first={box_first}: after one step y = {y1}, expected 0.483"
        );

        let mut max_speed: Real = 0.0;
        for _ in 1..120 {
            rest.scene.step(&mut pipeline);
            max_speed = max_speed.max(rest.scene.velocity(rest.body).norm());
        }
        let penetration = 0.5 - rest.scene.position(rest.body).y;
        assert!(
            (penetration - SLOP).abs() < 1e-3,
            "box_first={box_first}: settled at penetration {penetration}, expected {SLOP}"
        );
        assert!(
            max_speed < 1e-3,
            "box_first={box_first}: max speed {max_speed}"
        );
    }
}

#[test]
fn position_correction_does_not_depend_on_the_iteration_count() {
    // The old solver pushed bodies apart by 0.2 * depth in *every* iteration,
    // so 50 iterations overshot ten times as far as 5.
    let height_after = |iterations: u32, steps: usize| -> Real {
        let mut pipeline = PhysicsPipeline::with_config(config_with_iterations(iterations));
        let mut rest = rest_scene(true, 0.48, 0.5, 0.0);
        for _ in 0..steps {
            rest.scene.step(&mut pipeline);
        }
        rest.scene.position(rest.body).y
    };
    for steps in [1, 3, 60] {
        let few = height_after(5, steps);
        let many = height_after(50, steps);
        assert!(
            (few - many).abs() < 1e-4,
            "after {steps} steps: 5 iterations give y = {few}, 50 give y = {many}"
        );
        assert!(
            few <= 0.5 + 1e-3,
            "no overshoot: y = {few} after {steps} steps"
        );
    }
}

#[test]
fn restitution_is_applied_once_whatever_the_iteration_count() {
    // A box lands flat at 4 m/s on a floor with e = 0.5. Without gravity the
    // rebound is exactly e * 4 = 2 m/s on all four contact points.
    for iterations in [10, 30, 100] {
        let config = PhysicsConfig {
            gravity: Vec3::zeros(),
            solver_iterations: iterations,
            ..awake_config()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut rest = rest_scene(true, 0.499, 0.5, 0.5);
        if let Some(body) = rest.scene.bodies.get_mut(rest.body) {
            body.velocity = Vec3::new(0.0, -4.0, 0.0);
        }
        rest.scene.step(&mut pipeline);

        // The four corner contacts are coupled through the rotation of the
        // box, so the Gauss-Seidel sweeps leave a small residual that shrinks
        // geometrically with the iteration count: a few 1e-3 after 10 sweeps,
        // a few 1e-9 after 30.
        let residual_limit = if iterations >= 30 { 1e-6 } else { 5e-3 };
        let body = rest.scene.bodies.get(rest.body).unwrap();
        assert!(
            (body.velocity.y - 2.0).abs() < residual_limit,
            "iterations={iterations}: rebound speed {} should be 2.0",
            body.velocity.y
        );
        assert!(
            body.velocity.x.abs() < residual_limit && body.velocity.z.abs() < residual_limit,
            "iterations={iterations}: flat landing must not slide the box: {:?}",
            body.velocity
        );
        assert!(
            body.angular_velocity.norm() < residual_limit,
            "iterations={iterations}: flat landing must not spin the box: {:?}",
            body.angular_velocity
        );
    }
}

#[test]
fn sphere_rebound_height_follows_the_restitution_coefficient() {
    let dt = 1.0 / 240.0;
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    let (ball, _) = scene.dynamic(1.0, Vec3::new(0.0, 2.5, 0.0), sphere(0.5), 0.5, 0.5);
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.5,
    );

    // Drop height (sphere bottom to floor) 2.0 m: the centre rebounds to
    // 0.5 + e^2 * 2.0 = 1.0 m.
    let mut bounced = false;
    let mut apex: Real = 0.0;
    for _ in 0..600 {
        pipeline.step(dt, &mut scene.bodies, &scene.colliders);
        let vy = scene.velocity(ball).y;
        let y = scene.position(ball).y;
        if vy > 0.5 {
            bounced = true;
        }
        if bounced {
            apex = apex.max(y);
            if vy < 0.0 {
                break;
            }
        }
    }
    assert!(bounced, "the ball must bounce");
    assert!(
        (apex - 1.0).abs() < 0.05,
        "rebound apex {apex} should be about 1.0 m"
    );
}

#[test]
fn contacts_never_pull_separating_bodies_together() {
    // A ball that overlaps the floor slightly but is already moving away
    // keeps its velocity: contact impulses are one-sided (lambda >= 0).
    let config = PhysicsConfig {
        gravity: Vec3::zeros(),
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    let mut scene = Scene::new();
    let (ball, _) = scene.dynamic(1.0, Vec3::new(0.0, 0.49, 0.0), sphere(0.5), 0.5, 0.0);
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    if let Some(body) = scene.bodies.get_mut(ball) {
        body.velocity = Vec3::new(0.0, 3.0, 0.0);
    }
    scene.step(&mut pipeline);

    assert_eq!(scene.velocity(ball), Vec3::new(0.0, 3.0, 0.0));
    assert!(
        pipeline.warm_start_cache().normal_impulse_total() == 0.0,
        "a separating contact accumulates no impulse"
    );
}

#[test]
fn coulomb_friction_decelerates_a_sliding_box() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut rest = rest_scene(false, 0.4995, 0.5, 0.0);
    if let Some(body) = rest.scene.bodies.get_mut(rest.body) {
        body.velocity = Vec3::new(2.0, 0.0, 0.0);
    }

    // mu = 0.5: the deceleration is mu * g = 4.905 m/s^2.
    for _ in 0..12 {
        rest.scene.step(&mut pipeline);
    }
    let vx = rest.scene.velocity(rest.body).x;
    let expected = 2.0 - 0.5 * G * 12.0 * DT;
    assert!(
        (vx - expected).abs() < 0.03 * expected,
        "after 0.2 s the box moves at {vx} m/s, expected {expected}"
    );

    // It stops after 2 / 4.905 = 0.41 s and then stays stuck.
    for _ in 0..48 {
        rest.scene.step(&mut pipeline);
    }
    let body = rest.scene.bodies.get(rest.body).unwrap();
    assert!(body.velocity.norm() < 1e-2, "velocity {:?}", body.velocity);
    assert!(body.transform.rotation.angle() < 0.05);
    let travelled = body.transform.position.x;
    // Coulomb stopping distance v^2 / (2 mu g) = 0.407 m for continuous time;
    // the 60 Hz stepping stops it after 24 steps, 0.391 m.
    assert!(
        (travelled - 0.391).abs() < 0.01,
        "stopped after {travelled} m, expected 0.391 m"
    );
}

#[test]
fn frictionless_box_keeps_sliding_at_constant_speed() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut rest = rest_scene(true, 0.4995, 0.0, 0.0);
    if let Some(body) = rest.scene.bodies.get_mut(rest.body) {
        body.velocity = Vec3::new(2.0, 0.0, 0.0);
    }
    for _ in 0..60 {
        rest.scene.step(&mut pipeline);
    }
    let v = rest.scene.velocity(rest.body);
    assert!((v.x - 2.0).abs() < 1e-9, "vx = {}", v.x);
    assert!(v.y.abs() < 1e-3 && v.z.abs() < 1e-9);
}

/// A unit box (mass 1, half-extent 0.5) lying on a wide static plank that is
/// tilted by `angle` about +Z, so that the top face of the plank rises toward
/// +x. The box is aligned with the plank and dips 1 mm into it.
fn incline_scene(angle: Real, friction: Real) -> (Scene, BodyHandle) {
    let rotation = Quat::from_axis_angle(&Vec3::z_axis(), angle);
    let up = rotation * Vec3::y();
    let mut scene = Scene::new();

    let mut plank = RigidBody::new_static();
    plank.transform = Transform::new(Vec3::zeros(), rotation);
    scene.add(plank, cuboid(10.0, 0.5, 10.0), friction, 0.0);

    // Plank half-thickness 0.5 + box half-extent 0.5 - 1 mm of penetration.
    let mut body = RigidBody::new(1.0);
    body.transform = Transform::new(up * 0.999, rotation);
    let (box_body, _) = scene.add(body, cuboid(0.5, 0.5, 0.5), friction, 0.0);
    (scene, box_body)
}

#[test]
fn box_sticks_to_an_incline_below_the_friction_angle() {
    // tan(20 degrees) = 0.364 < mu = 0.5: static friction holds the box.
    let angle = 20.0_f64.to_radians();
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let (mut scene, body) = incline_scene(angle, 0.5);
    let start = scene.position(body);

    for _ in 0..120 {
        scene.step(&mut pipeline);
    }
    let moved = (scene.position(body) - start).norm();
    assert!(moved < 5e-3, "the box crept {moved} m down the incline");
    assert!(
        scene.velocity(body).norm() < 5e-3,
        "velocity {:?}",
        scene.velocity(body)
    );
}

#[test]
fn box_slides_down_an_incline_with_the_coulomb_acceleration() {
    // mu = 0.2 < tan(20 degrees): the box slides with a = g (sin t - mu cos t).
    let angle = 20.0_f64.to_radians();
    let friction = 0.2;
    let steps = 60_u32;
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let (mut scene, body) = incline_scene(angle, friction);
    let start = scene.position(body);
    let downhill = Vec3::new(-angle.cos(), -angle.sin(), 0.0);

    for _ in 0..steps {
        scene.step(&mut pipeline);
    }

    let acceleration = G * (angle.sin() - friction * angle.cos());
    let elapsed = DT * f64::from(steps);
    let velocity = scene.velocity(body);
    let speed = velocity.dot(&downhill);
    assert!(
        (speed - acceleration * elapsed).abs() < 0.03 * acceleration * elapsed,
        "speed {speed} m/s, expected {}",
        acceleration * elapsed
    );
    // Nothing sideways and nothing into (or out of) the plank.
    assert!(
        (velocity - downhill * speed).norm() < 0.02 * acceleration * elapsed,
        "velocity {velocity:?} leaves the incline"
    );
    // Semi-implicit Euler: x_n = a dt^2 n (n + 1) / 2.
    let n = f64::from(steps);
    let expected_distance = acceleration * DT * DT * n * (n + 1.0) / 2.0;
    let travelled = (scene.position(body) - start).dot(&downhill);
    assert!(
        (travelled - expected_distance).abs() < 0.04 * expected_distance,
        "slid {travelled} m, expected {expected_distance}"
    );
    // The box stays flat on the plank.
    let tilt = scene.bodies.get(body).map_or(Real::INFINITY, |b| {
        b.transform
            .rotation
            .angle_to(&Quat::from_axis_angle(&Vec3::z_axis(), angle))
    });
    assert!(tilt < 0.02, "the box tilted by {tilt} rad");
}

#[test]
fn warm_starting_reduces_the_residual_of_a_single_iteration() {
    let residual_after = |warm: bool| -> Real {
        let mut pipeline = PhysicsPipeline::with_config(config_with_iterations(1));
        let mut rest = rest_scene(true, 0.497, 0.5, 0.0);
        for _ in 0..8 {
            if !warm {
                pipeline.contact_cache.clear();
            }
            rest.scene.step(&mut pipeline);
        }
        rest.scene.velocity(rest.body).norm()
            + rest
                .scene
                .bodies
                .get(rest.body)
                .map_or(0.0, |b| b.angular_velocity.norm())
    };
    let cold = residual_after(false);
    let warm = residual_after(true);
    assert!(
        cold > 1e-3,
        "a single cold iteration leaves a visible residual, got {cold}"
    );
    assert!(
        warm < 0.1 * cold,
        "warm start residual {warm} should be far below the cold residual {cold}"
    );
}

#[test]
fn warm_start_cache_is_emptied_when_the_bodies_separate() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut rest = rest_scene(true, 0.497, 0.5, 0.0);
    rest.scene.step(&mut pipeline);
    assert!(
        !pipeline.warm_start_cache().is_empty(),
        "cache should be populated while the bodies are in contact"
    );

    if let Some(body) = rest.scene.bodies.get_mut(rest.body) {
        body.transform.position.y = 100.0;
    }
    rest.scene.step(&mut pipeline);
    assert!(
        pipeline.warm_start_cache().is_empty(),
        "cache should be empty after the bodies separate"
    );
}

#[test]
fn warm_start_scales_with_the_step_size() {
    // Impulses accumulated over a 1/60 s step support the weight for 1/60 s.
    // After switching to 1/120 s, the warm start must offer half of them, or
    // the box would be kicked upward by the stale impulse.
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut rest = rest_scene(true, 0.497, 0.5, 0.0);
    for _ in 0..30 {
        rest.scene.step(&mut pipeline);
    }
    let half = DT / 2.0;
    for _ in 0..30 {
        pipeline.step(half, &mut rest.scene.bodies, &rest.scene.colliders);
        assert!(
            rest.scene.velocity(rest.body).norm() < 1e-3,
            "velocity {:?} after changing the step size",
            rest.scene.velocity(rest.body)
        );
    }
    let total = pipeline
        .warm_start_cache()
        .normal_impulse(rest.box_collider, rest.floor_collider);
    assert!(
        (total - G * half).abs() < 0.03 * G * half,
        "supporting impulse {total}, expected {}",
        G * half
    );
}

#[test]
fn sleeping_box_on_the_floor_stays_asleep_and_is_not_solved() {
    // Default thresholds so that a sleeping body with zero velocity stays
    // asleep in `check_sleep`.
    let mut pipeline = PhysicsPipeline::with_config(PhysicsConfig::default());
    let mut rest = rest_scene(true, 0.5, 0.5, 0.0);
    if let Some(body) = rest.scene.bodies.get_mut(rest.body) {
        body.state = BodyState::Sleeping;
    }
    for _ in 0..30 {
        rest.scene.step(&mut pipeline);
        let body = rest.scene.bodies.get(rest.body).unwrap();
        assert_eq!(body.state, BodyState::Sleeping);
    }
    assert_eq!(rest.scene.position(rest.body), Vec3::new(0.0, 0.5, 0.0));
    assert_eq!(pipeline.stats.last_contact_manifolds, 0);
    assert!(pipeline.warm_start_cache().is_empty());
}

#[test]
fn falling_sphere_wakes_a_sleeping_box() {
    let mut pipeline = PhysicsPipeline::with_config(PhysicsConfig::default());
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    let (box_body, _) = scene.dynamic(
        1.0,
        Vec3::new(0.0, 0.5, 0.0),
        cuboid(0.5, 0.5, 0.5),
        0.5,
        0.0,
    );
    if let Some(body) = scene.bodies.get_mut(box_body) {
        body.state = BodyState::Sleeping;
    }
    // A ball dropped on top of the box from 0.5 m above its top face.
    let (ball, _) = scene.dynamic(1.0, Vec3::new(0.0, 1.75, 0.0), sphere(0.25), 0.5, 0.0);

    let mut woke_at = None;
    for step in 0..60 {
        scene.step(&mut pipeline);
        let box_state = scene.bodies.get(box_body).unwrap().state;
        if woke_at.is_none() && box_state == BodyState::Active {
            woke_at = Some(step);
        }
    }
    let woke_at = woke_at.expect("the ball must wake the box up");
    assert!(
        (15..=30).contains(&woke_at),
        "the box wakes when the ball arrives (step {woke_at})"
    );
    // Everything came to rest: the ball on the box, the box on the floor.
    assert!((scene.position(ball).y - 1.25).abs() < 0.03);
    assert!((scene.position(box_body).y - 0.5).abs() < 0.03);
}

#[test]
fn five_box_stack_holds_its_shape() {
    let mut pipeline = PhysicsPipeline::with_config(config_with_iterations(10));
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.8,
        0.0,
    );
    let mut boxes = Vec::new();
    for i in 0..5_u32 {
        // 1 mm of overlap between neighbours.
        let y = 0.5 - 0.001 + f64::from(i) * 0.999;
        let (body, _) = scene.dynamic(1.0, Vec3::new(0.0, y, 0.0), cuboid(0.5, 0.5, 0.5), 0.8, 0.0);
        boxes.push((body, y));
    }
    for _ in 0..240 {
        scene.step(&mut pipeline);
    }
    for (index, (body, start_y)) in boxes.iter().enumerate() {
        let b = scene.bodies.get(*body).unwrap();
        let p = b.transform.position;
        // Residual solver asymmetry tilts the boxes by a few milliradians, which
        // walks the upper boxes sideways by about a centimetre; a stack that
        // is falling over drifts by far more.
        assert!(
            p.x.abs() < 0.05 && p.z.abs() < 0.05,
            "box {index} slid sideways to {p:?}"
        );
        assert!(
            *start_y - p.y < 0.03 && p.y - *start_y < 0.01,
            "box {index} moved from y = {start_y} to {}",
            p.y
        );
        assert!(
            b.transform.rotation.angle() < 0.03,
            "box {index} tilted by {}",
            b.transform.rotation.angle()
        );
        assert!(b.velocity.norm() < 0.02, "box {index} is still moving");
    }
}

/// Final state of a small pile of spheres and boxes after 200 steps.
fn run_pile() -> Vec<(Vec3, Vec3)> {
    let mut pipeline = PhysicsPipeline::new();
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.2,
    );
    let mut handles = Vec::new();
    for i in 0..3_u32 {
        for j in 0..3_u32 {
            let x = f64::from(i) * 0.9 - 0.9;
            let z = f64::from(j) * 0.9 - 0.9;
            let y = 1.0 + f64::from(i + j) * 0.7;
            let shape = if (i + j) % 2 == 0 {
                sphere(0.4)
            } else {
                cuboid(0.35, 0.35, 0.35)
            };
            let (body, _) = scene.dynamic(1.0, Vec3::new(x, y, z), shape, 0.5, 0.2);
            handles.push(body);
        }
    }
    for _ in 0..200 {
        scene.step(&mut pipeline);
    }
    handles
        .iter()
        .map(|&h| (scene.position(h), scene.velocity(h)))
        .collect()
}

#[test]
fn simulation_is_deterministic() {
    let first = run_pile();
    let second = run_pile();
    assert_eq!(first.len(), 9);
    assert_eq!(first, second, "two identical runs must agree bit for bit");
    for (position, velocity) in &first {
        assert!(position.iter().all(|c| c.is_finite()));
        assert!(velocity.iter().all(|c| c.is_finite()));
        assert!(
            position.y > 0.2,
            "nothing fell through the floor: {position:?}"
        );
    }
}
