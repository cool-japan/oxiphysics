// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of continuous collision detection: time-of-impact sub-stepping,
//! resolution of several events per step, and the regular contact pass that
//! must keep running in steps that contain CCD events.

use super::*;
use crate::pipeline::ccd::swept_aabb;
use oxiphysics_core::math::Quat;
use std::f64::consts::FRAC_PI_2;

/// A configuration with CCD on and sleeping off.
fn ccd_config() -> PhysicsConfig {
    PhysicsConfig {
        ccd_enabled: true,
        ..awake_config()
    }
}

/// A thin static wall across the x axis at `x`, 3 m high and 3 m wide.
///
/// The wall is kept compact on purpose: the time-of-impact search advances in
/// steps that shrink with the size of the shapes, so against a small wall it
/// finds the impact in few iterations whatever part of the step is left.
fn wall(scene: &mut Scene, x: Real, half_thickness: Real) {
    scene.static_box(
        Vec3::new(x, 0.0, 0.0),
        Vec3::new(half_thickness, 1.5, 1.5),
        0.0,
        1.0,
    );
}

/// A gravity-free elastic bullet flying along +x.
fn bullet(scene: &mut Scene, position: Vec3, vx: Real, radius: Real) -> BodyHandle {
    let mut body = RigidBody::new(1.0);
    body.transform = Transform::from_position(position);
    body.velocity = Vec3::new(vx, 0.0, 0.0);
    body.gravity_scale = 0.0;
    scene.add(body, sphere(radius), 0.0, 1.0).0
}

/// Two walls 1 m apart (inner faces at x = -0.5 and x = 0.5) around the origin.
fn corridor(scene: &mut Scene) {
    wall(scene, -0.55, 0.05);
    wall(scene, 0.55, 0.05);
}

/// A bullet that ping-pongs between two walls 20 m from the origin, so that
/// nearly every step contains CCD events.
///
/// The inner faces of the walls are at x = 19.5 and x = 20.5; the bullet
/// (radius 0.05) flies at 150 m/s, i.e. 2.5 m per 60 Hz step, and therefore
/// stays between x = 19.55 and x = 20.45 as long as no impact is missed.
///
/// The bullet flies along the axis through the centres of both walls (1 m above
/// the floor), which makes the time-of-impact search an exactly axis-aligned
/// problem: these tests are about what the pipeline does around the events, not
/// about the accuracy of the search for oblique approaches.
fn ping_pong_bullet(scene: &mut Scene) -> BodyHandle {
    for x in [19.45, 20.55] {
        scene.static_box(Vec3::new(x, 1.0, 0.0), Vec3::new(0.05, 1.5, 1.5), 0.0, 1.0);
    }
    bullet(scene, Vec3::new(19.9, 1.0, 0.0), 150.0, 0.05)
}

#[test]
fn without_ccd_a_fast_bullet_tunnels_through_a_thin_wall() {
    let config = PhysicsConfig {
        ccd_enabled: false,
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    let mut scene = Scene::new();
    wall(&mut scene, 5.0, 0.05);
    let b = bullet(&mut scene, Vec3::new(3.0, 0.0, 0.0), 200.0, 0.1);

    scene.step(&mut pipeline);

    // 200 m/s is 3.3 m per step: the bullet jumps over the 0.1 m wall.
    assert!(scene.position(b).x > 6.0, "x = {}", scene.position(b).x);
    assert_eq!(scene.velocity(b).x, 200.0);
}

#[test]
fn with_ccd_the_same_bullet_bounces_back_at_the_time_of_impact() {
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    wall(&mut scene, 5.0, 0.05);
    let b = bullet(&mut scene, Vec3::new(3.0, 0.0, 0.0), 200.0, 0.1);

    scene.step(&mut pipeline);

    assert_eq!(pipeline.stats.last_ccd_events, 1);
    // Elastic: the speed is reversed exactly (the contact normal is x).
    assert!(
        (scene.velocity(b).x + 200.0).abs() < 1e-6,
        "vx = {}",
        scene.velocity(b).x
    );
    // The bullet reaches the wall face (x = 4.95) with its surface after
    // 1.85 m, i.e. 0.555 of the step, then travels 0.445 * 3.33 m back.
    let expected = 4.85 - 0.4444 * 200.0 / 60.0;
    assert!(
        (scene.position(b).x - expected).abs() < 0.02,
        "x = {}, expected {expected}",
        scene.position(b).x
    );
}

#[test]
fn every_time_of_impact_of_a_step_is_resolved() {
    // A bullet ping-pongs between two walls 1 m apart at 150 m/s = 2.5 m per
    // step. Starting at x = 0.4 it hits the right wall after 0.05 m, the left
    // wall after another 0.9 m and the right wall again after another 0.9 m:
    // three events in the very first step. Resolving only the earliest one
    // lets it tunnel through the left wall.
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    corridor(&mut scene);
    let b = bullet(&mut scene, Vec3::new(0.4, 0.0, 0.0), 150.0, 0.05);

    scene.step(&mut pipeline);
    assert_eq!(pipeline.stats.last_ccd_events, 3);
    assert!((scene.velocity(b).x + 150.0).abs() < 1e-6);
    assert!(
        (scene.position(b).x + 0.2).abs() < 0.01,
        "x = {}",
        scene.position(b).x
    );

    let mut total_events = 3;
    for _ in 0..120 {
        scene.step(&mut pipeline);
        total_events += pipeline.stats.last_ccd_events;
        let x = scene.position(b).x;
        assert!(
            x.abs() < 0.5,
            "the bullet escaped the corridor: x = {x} (speed {})",
            scene.velocity(b).x
        );
        assert!((scene.velocity(b).x.abs() - 150.0).abs() < 1e-4);
    }
    // 150 m/s through a 0.9 m wide corridor: about 2.8 crossings per step.
    assert!(total_events > 250, "only {total_events} events in 2 s");
}

#[test]
fn a_fast_pair_is_still_caught_in_what_is_left_of_the_step() {
    // Bullet A hits a wall 84 % into the step. That leaves 0.4 m of travel,
    // less than the displacement that makes CCD look at a pair at all, if the
    // pair were judged by the remainder of the step instead of the whole step.
    // Bullet B flies in parallel and meets a thin plate 92 % into the step,
    // i.e. in the remainder, and must not pass through it.
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(2.20, 0.0, 0.0),
        Vec3::new(0.05, 0.5, 0.5),
        0.0,
        1.0,
    );
    scene.static_box(
        Vec3::new(2.36, 3.0, 0.0),
        Vec3::new(0.01, 0.5, 0.5),
        0.0,
        1.0,
    );
    let a = bullet(&mut scene, Vec3::new(0.0, 0.0, 0.0), 150.0, 0.05);
    let b = bullet(&mut scene, Vec3::new(0.0, 3.0, 0.0), 150.0, 0.05);

    scene.step(&mut pipeline);

    assert_eq!(
        pipeline.stats.last_ccd_events, 2,
        "one impact per bullet, the second one in the remainder of the step"
    );
    // Both bounced elastically ...
    assert!((scene.velocity(a).x + 150.0).abs() < 1e-6);
    assert!(
        (scene.velocity(b).x + 150.0).abs() < 1e-6,
        "bullet B went through the plate: vx = {}",
        scene.velocity(b).x
    );
    // ... A after its impact at x = 2.10 with the 0.4 m left of the step, B
    // after its impact at x = 2.30 with the 0.2 m left of the step.
    let x_a = scene.position(a).x;
    let x_b = scene.position(b).x;
    assert!((x_a - 1.70).abs() < 0.02, "bullet A ended at x = {x_a}");
    assert!((x_b - 2.10).abs() < 0.02, "bullet B ended at x = {x_b}");
}

#[test]
fn the_event_budget_limits_the_events_per_step() {
    for (budget, expected_events) in [(0_usize, 0_usize), (1, 1), (2, 2), (8, 3)] {
        let mut pipeline = PhysicsPipeline::with_config(ccd_config());
        pipeline.max_ccd_substeps = budget;
        let mut scene = Scene::new();
        corridor(&mut scene);
        let b = bullet(&mut scene, Vec3::new(0.4, 0.0, 0.0), 150.0, 0.05);

        scene.step(&mut pipeline);
        assert_eq!(
            pipeline.stats.last_ccd_events, expected_events,
            "budget {budget}"
        );
        if budget == 0 {
            // No sub-steps at all: the bullet flies straight through the wall.
            assert!(scene.position(b).x > 1.0, "x = {}", scene.position(b).x);
        }
    }
}

#[test]
fn resting_stack_keeps_being_solved_while_a_bullet_bounces_elsewhere() {
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(50.0, 0.5, 50.0),
        0.5,
        0.0,
    );

    // Five unit spheres stacked on the floor, 1 mm of overlap between them.
    let mut stack = Vec::new();
    for i in 0..5_u32 {
        let y = 0.5 - 0.001 + f64::from(i) * 0.999;
        let (body, _) = scene.dynamic(1.0, Vec3::new(0.0, y, 0.0), sphere(0.5), 0.5, 0.0);
        stack.push((body, y));
    }

    // A bullet ping-pongs between two walls 20 m away, so nearly every step
    // contains CCD events.
    let b = ping_pong_bullet(&mut scene);

    let mut steps_with_events = 0;
    for _ in 0..120 {
        scene.step(&mut pipeline);
        if pipeline.stats.last_ccd_events > 0 {
            steps_with_events += 1;
            // The regular contact pass ran in this very step: the floor
            // contact and the four sphere-sphere contacts of the stack.
            assert!(
                pipeline.stats.last_contact_manifolds >= 5,
                "only {} manifolds in a step with CCD events",
                pipeline.stats.last_contact_manifolds
            );
        }
        let x = scene.position(b).x;
        assert!((19.5..20.5).contains(&x), "bullet escaped: x = {x}");
    }
    assert!(
        steps_with_events >= 110,
        "the bullet should produce events nearly every step ({steps_with_events})"
    );

    // The stack did not sink: its top sphere is where it started.
    for (index, (body, start_y)) in stack.iter().enumerate() {
        let p = scene.position(*body);
        assert!(
            (p.y - start_y).abs() < 0.01,
            "sphere {index} moved from y = {start_y} to {}",
            p.y
        );
        assert!(p.x.abs() < 1e-6 && p.z.abs() < 1e-6);
        assert!(scene.velocity(*body).norm() < 0.01);
    }
}

#[test]
fn ccd_events_elsewhere_do_not_disturb_a_resting_box() {
    // The same resting box with and without a bullet that produces CCD events
    // in every step. The box must not notice: the contacts are solved before
    // anything moves, so the sub-steps integrate its *solved* (zero) velocity.
    // Integrating the unsolved gravity velocity instead would sink it by a
    // couple of millimetres into the floor in every such step.
    const SLOP: Real = 0.005;
    let run = |with_bullet: bool| -> (Vec3, Vec3, usize) {
        let mut pipeline = PhysicsPipeline::with_config(ccd_config());
        let mut scene = Scene::new();
        scene.static_box(
            Vec3::new(0.0, -0.5, 0.0),
            Vec3::new(50.0, 0.5, 50.0),
            0.5,
            0.0,
        );
        // 3 mm of penetration: inside the slop, so nothing needs to move.
        let (resting_box, _) = scene.dynamic(
            1.0,
            Vec3::new(0.0, 0.497, 0.0),
            cuboid(0.5, 0.5, 0.5),
            0.5,
            0.0,
        );
        let bullet_handle = with_bullet.then(|| ping_pong_bullet(&mut scene));

        let mut steps_with_events = 0;
        for _ in 0..120 {
            scene.step(&mut pipeline);
            if pipeline.stats.last_ccd_events > 0 {
                steps_with_events += 1;
            }
            if let Some(handle) = bullet_handle {
                let x = scene.position(handle).x;
                assert!((19.5..20.5).contains(&x), "bullet escaped: x = {x}");
            }
        }
        (
            scene.position(resting_box),
            scene.velocity(resting_box),
            steps_with_events,
        )
    };

    let (quiet_position, quiet_velocity, quiet_events) = run(false);
    let (busy_position, busy_velocity, busy_events) = run(true);

    assert_eq!(quiet_events, 0, "nothing in the quiet scene is fast");
    assert!(
        busy_events >= 110,
        "the bullet should produce events nearly every step ({busy_events})"
    );
    assert!(
        (busy_position - quiet_position).norm() < 1e-9,
        "the box moved from {quiet_position:?} to {busy_position:?}"
    );
    assert!(
        (busy_velocity - quiet_velocity).norm() < 1e-9,
        "the box speed changed from {quiet_velocity:?} to {busy_velocity:?}"
    );

    // And it rests where a resting box belongs: inside the slop, at rest.
    let penetration = 0.5 - busy_position.y;
    assert!(
        (-1e-3..=SLOP + 1e-3).contains(&penetration),
        "penetration {penetration}"
    );
    assert!(busy_velocity.norm() < 1e-3);
}

#[test]
fn head_on_elastic_impact_exchanges_velocities_and_conserves_momentum() {
    let config = PhysicsConfig {
        gravity: Vec3::zeros(),
        ..ccd_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    let mut scene = Scene::new();
    let a = bullet(&mut scene, Vec3::new(-5.0, 0.0, 0.0), 100.0, 0.5);
    let b = bullet(&mut scene, Vec3::new(5.0, 0.0, 0.0), -100.0, 0.5);

    pipeline.step(0.1, &mut scene.bodies, &scene.colliders);

    assert_eq!(pipeline.stats.last_ccd_events, 1);
    let va = scene.velocity(a).x;
    let vb = scene.velocity(b).x;
    assert!((va + 100.0).abs() < 1e-6, "va = {va}");
    assert!((vb - 100.0).abs() < 1e-6, "vb = {vb}");
    assert!((va + vb).abs() < 1e-9, "momentum is conserved");
    // They bounced instead of passing through each other.
    assert!(scene.position(a).x < scene.position(b).x);
}

#[test]
fn a_bullet_passes_through_a_sensor_without_events() {
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    let mut sensor_body = RigidBody::new_static();
    sensor_body.transform = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
    let sensor_body = scene.bodies.insert(sensor_body);
    scene.colliders.insert(
        Collider::new(cuboid(0.05, 5.0, 5.0))
            .with_body(sensor_body)
            .as_sensor(),
    );
    let b = bullet(&mut scene, Vec3::new(3.0, 0.0, 0.0), 200.0, 0.1);

    scene.step(&mut pipeline);
    assert_eq!(pipeline.stats.last_ccd_events, 0);
    assert_eq!(scene.velocity(b).x, 200.0);
    assert!(scene.position(b).x > 6.0);
}

#[test]
fn a_touching_pair_that_is_already_separating_is_not_an_event() {
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    wall(&mut scene, 5.0, 0.05);
    // The bullet sits on the wall face (x = 4.85 + 0.1 = 4.95) moving away.
    let b = bullet(&mut scene, Vec3::new(4.85, 0.0, 0.0), -200.0, 0.1);

    scene.step(&mut pipeline);
    assert_eq!(pipeline.stats.last_ccd_events, 0);
    assert_eq!(scene.velocity(b).x, -200.0);
    assert!((scene.position(b).x - (4.85 - 200.0 / 60.0)).abs() < 1e-9);
}

#[test]
fn swept_bounds_include_the_rotation_over_the_step() {
    // A 2 m plank lying along x that turns upright over the step while it
    // moves 1 m to the right. Bounding only the translated, un-rotated end
    // pose would cover y in [-0.05, 0.05]; the real sweep reaches y = +-1.
    let collider = Collider::new(cuboid(1.0, 0.05, 0.05));
    let start = Transform::from_position(Vec3::zeros());
    let end = Transform::new(
        Vec3::new(1.0, 0.0, 0.0),
        Quat::from_axis_angle(&Vec3::z_axis(), FRAC_PI_2),
    );
    let swept = swept_aabb(&collider, &start, &end);

    assert!(
        (swept.min - Vec3::new(-1.0, -1.0, -0.05)).norm() < 1e-9,
        "{:?}",
        swept.min
    );
    assert!(
        (swept.max - Vec3::new(1.05, 1.0, 0.05)).norm() < 1e-9,
        "{:?}",
        swept.max
    );
}

#[test]
fn a_scripted_kinematic_body_covers_its_whole_step_even_when_ccd_splits_it() {
    // The bullet hits the wall about 55 % into the step, so only the
    // remaining 45 % of the step is left to integrate at the end. A kinematic
    // platform that a script moves 1 m during the step stood still while the
    // step was split; it must still arrive on its target, and its velocity is
    // that of a body that covered 1 m in the whole step, not in what was left
    // of it.
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    wall(&mut scene, 5.0, 0.05);
    bullet(&mut scene, Vec3::new(3.0, 0.0, 0.0), 200.0, 0.1);

    let target = Vec3::new(0.0, 20.0, 1.0);
    let mut platform_body = RigidBody::new_kinematic();
    platform_body.transform = Transform::from_position(Vec3::new(0.0, 20.0, 0.0));
    platform_body.set_kinematic_target(target, Quat::identity());
    let (platform, _) = scene.add(platform_body, cuboid(1.0, 0.1, 1.0), 0.5, 0.0);

    scene.step(&mut pipeline);

    assert_eq!(
        pipeline.stats.last_ccd_events, 1,
        "the bullet's impact splits the step"
    );
    let reached = scene.position(platform);
    assert!(
        (reached - target).norm() < 1e-12,
        "the platform stopped at {reached:?} instead of {target:?}"
    );
    let velocity = scene.velocity(platform);
    assert!(
        (velocity - Vec3::new(0.0, 0.0, 1.0 / DT)).norm() < 1e-9,
        "velocity {velocity:?} should be 1 m over the whole step ({} m/s)",
        1.0 / DT
    );

    // The same scene without CCD gives the identical platform motion.
    let mut plain = PhysicsPipeline::with_config(PhysicsConfig {
        ccd_enabled: false,
        ..awake_config()
    });
    let mut plain_scene = Scene::new();
    let mut plain_body = RigidBody::new_kinematic();
    plain_body.transform = Transform::from_position(Vec3::new(0.0, 20.0, 0.0));
    plain_body.set_kinematic_target(target, Quat::identity());
    let (plain_platform, _) = plain_scene.add(plain_body, cuboid(1.0, 0.1, 1.0), 0.5, 0.0);
    plain_scene.step(&mut plain);
    assert_eq!(plain_scene.position(plain_platform), reached);
    assert_eq!(plain_scene.velocity(plain_platform), velocity);
}

#[test]
fn ccd_is_not_confused_by_bodies_without_colliders() {
    let mut pipeline = PhysicsPipeline::with_config(ccd_config());
    let mut scene = Scene::new();
    // A body that has no collider flies past a wall: nothing to collide with.
    let mut ghost = RigidBody::new(1.0);
    ghost.velocity = Vec3::new(500.0, 0.0, 0.0);
    ghost.gravity_scale = 0.0;
    let ghost = scene.bodies.insert(ghost);
    wall(&mut scene, 5.0, 0.05);
    scene.step(&mut pipeline);
    assert_eq!(pipeline.stats.last_ccd_events, 0);
    assert!(scene.position(ghost).x > 8.0);
}
