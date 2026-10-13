// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the SPH-rigid coupling: Archimedes buoyancy applied by the
//! pipeline from the volume of each collider.

use super::*;
use std::f64::consts::PI;

const G: Real = 9.81;

/// Volume of the unit-radius-0.5 ball used below.
fn ball_volume() -> Real {
    4.0 / 3.0 * PI * 0.5_f64.powi(3)
}

/// A fluid region around the origin, large enough to hold the test bodies.
fn pool(density: Real) -> SphRigidCouplingParams {
    SphRigidCouplingParams {
        fluid_density: density,
        coupling_alpha: 1.0,
        region_half_extents: Vec3::new(5.0, 5.0, 5.0),
        region_centre: Vec3::zeros(),
    }
}

/// A pipeline without sleeping, with coupling `params` installed.
fn coupled_pipeline(params: SphRigidCouplingParams) -> PhysicsPipeline {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    pipeline.set_sph_rigid_coupling(params);
    pipeline
}

/// One ball of the given density at the origin; returns its body handle.
fn ball_of_density(scene: &mut Scene, density: Real) -> BodyHandle {
    let mass = density * ball_volume();
    scene.dynamic(mass, Vec3::zeros(), sphere(0.5), 0.5, 0.0).0
}

#[test]
fn neutrally_buoyant_body_does_not_accelerate() {
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);

    for _ in 0..30 {
        scene.step(&mut pipeline);
    }
    assert!(
        scene.velocity(ball).norm() < 1e-9,
        "weight and buoyancy cancel, got {:?}",
        scene.velocity(ball)
    );
    assert!(scene.position(ball).norm() < 1e-9);
}

#[test]
fn fluid_twice_as_dense_as_the_body_accelerates_it_upward_at_g() {
    let mut pipeline = coupled_pipeline(pool(2000.0));
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);

    scene.step(&mut pipeline);

    // Net acceleration (2 g - g) = +g for one step.
    let vy = scene.velocity(ball).y;
    assert!((vy - G * DT).abs() < 1e-9, "vy = {vy}, expected {}", G * DT);
    assert_eq!(scene.velocity(ball).x, 0.0);
    assert_eq!(scene.velocity(ball).z, 0.0);
}

#[test]
fn body_denser_than_the_fluid_still_sinks_but_more_slowly() {
    // Density ratio 4: net downward acceleration g * (1 - 1/4) = 0.75 g.
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 4000.0);

    scene.step(&mut pipeline);
    let vy = scene.velocity(ball).y;
    assert!(
        (vy + 0.75 * G * DT).abs() < 1e-9,
        "vy = {vy}, expected {}",
        -0.75 * G * DT
    );
}

#[test]
fn buoyancy_follows_the_shape_volume_not_its_bounding_sphere() {
    // A 1 m cube (volume 1 m^3) of density 1000 in water is neutral. A
    // buoyancy based on the bounding sphere of the cube (radius 0.87,
    // volume 2.7 m^3) would send it up at 1.7 g.
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();
    let (cube, _) = scene.dynamic(1000.0, Vec3::zeros(), cuboid(0.5, 0.5, 0.5), 0.5, 0.0);

    scene.step(&mut pipeline);
    assert!(
        scene.velocity(cube).norm() < 1e-9,
        "vy = {}",
        scene.velocity(cube).y
    );
}

#[test]
fn partially_submerged_body_receives_a_proportional_force() {
    // The fluid fills y in [-10, 0]; the cube's centre sits at the surface,
    // so half of it is wet and the buoyancy is half its weight.
    let params = SphRigidCouplingParams {
        fluid_density: 1000.0,
        coupling_alpha: 1.0,
        region_half_extents: Vec3::new(5.0, 5.0, 5.0),
        region_centre: Vec3::new(0.0, -5.0, 0.0),
    };
    let mut pipeline = coupled_pipeline(params);
    let mut scene = Scene::new();
    let (cube, _) = scene.dynamic(1000.0, Vec3::zeros(), cuboid(0.5, 0.5, 0.5), 0.5, 0.0);

    scene.step(&mut pipeline);
    let vy = scene.velocity(cube).y;
    assert!(
        (vy + 0.5 * G * DT).abs() < 1e-9,
        "vy = {vy}, expected {}",
        -0.5 * G * DT
    );
}

#[test]
fn body_above_the_fluid_is_not_lifted() {
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();
    let (ball, _) = scene.dynamic(1.0, Vec3::new(0.0, 20.0, 0.0), sphere(0.5), 0.5, 0.0);
    scene.step(&mut pipeline);
    assert!((scene.velocity(ball).y + G * DT).abs() < 1e-12);

    // Same for a body beside the region.
    let (aside, _) = scene.dynamic(1.0, Vec3::new(30.0, 0.0, 0.0), sphere(0.5), 0.5, 0.0);
    scene.step(&mut pipeline);
    assert!((scene.velocity(aside).y + G * DT).abs() < 1e-12);
}

#[test]
fn buoyancy_acts_against_gravity_in_any_direction() {
    // Gravity along -x: the fluid pushes toward +x, nothing else moves.
    let config = PhysicsConfig {
        gravity: Vec3::new(-G, 0.0, 0.0),
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    pipeline.set_sph_rigid_coupling(pool(2000.0));
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);

    scene.step(&mut pipeline);
    let v = scene.velocity(ball);
    assert!((v.x - G * DT).abs() < 1e-9, "v = {v:?}");
    assert!(v.y.abs() < 1e-12 && v.z.abs() < 1e-12);
}

#[test]
fn coupling_strength_scales_the_force() {
    let params = SphRigidCouplingParams {
        coupling_alpha: 0.5,
        ..pool(2000.0)
    };
    let mut pipeline = coupled_pipeline(params);
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);
    scene.step(&mut pipeline);
    // alpha = 0.5: buoyancy = g, exactly the weight.
    assert!(scene.velocity(ball).norm() < 1e-9);
}

#[test]
fn no_gravity_means_no_buoyancy() {
    let config = PhysicsConfig {
        gravity: Vec3::zeros(),
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    pipeline.set_sph_rigid_coupling(pool(5000.0));
    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);
    scene.step(&mut pipeline);
    assert_eq!(scene.velocity(ball), Vec3::zeros());
}

#[test]
fn clearing_the_coupling_restores_free_fall() {
    let mut pipeline = coupled_pipeline(pool(1000.0));
    assert!(pipeline.sph_rigid_coupling().is_some());
    pipeline.clear_sph_rigid_coupling();
    assert!(pipeline.sph_rigid_coupling().is_none());

    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 1000.0);
    scene.step(&mut pipeline);
    assert!((scene.velocity(ball).y + G * DT).abs() < 1e-12);
}

#[test]
fn static_bodies_sensors_and_colliderless_bodies_are_ignored() {
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();

    // A static box in the pool.
    let (fixed, _) = scene.static_box(Vec3::zeros(), Vec3::new(0.5, 0.5, 0.5), 0.5, 0.0);
    // A dynamic body whose only collider is a sensor: it displaces nothing.
    let mut ghost = RigidBody::new(1.0);
    ghost.transform = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
    let ghost = scene.bodies.insert(ghost);
    scene
        .colliders
        .insert(Collider::new(sphere(0.5)).with_body(ghost).as_sensor());
    // A dynamic body without any collider.
    let mut bare = RigidBody::new(1.0);
    bare.transform = Transform::from_position(Vec3::new(-2.0, 0.0, 0.0));
    let bare = scene.bodies.insert(bare);

    scene.step(&mut pipeline);
    assert_eq!(scene.velocity(fixed), Vec3::zeros());
    assert!((scene.velocity(ghost).y + G * DT).abs() < 1e-12);
    assert!((scene.velocity(bare).y + G * DT).abs() < 1e-12);
}

#[test]
fn every_collider_of_a_compound_body_displaces_fluid() {
    // Two half-size balls on one body of the same total volume as the big
    // ball: neutral in water when the body has the matching mass.
    let r = 0.5 * 0.5_f64.cbrt();
    let mut pipeline = coupled_pipeline(pool(1000.0));
    let mut scene = Scene::new();
    let body = RigidBody::new(1000.0 * ball_volume());
    let handle = scene.bodies.insert(body);
    for x in [-1.0, 1.0] {
        scene.colliders.insert(
            Collider::new(sphere(r))
                .with_body(handle)
                .with_local_transform(Transform::from_position(Vec3::new(x, 0.0, 0.0))),
        );
    }
    scene.step(&mut pipeline);
    assert!(
        scene.velocity(handle).norm() < 1e-9,
        "v = {:?}",
        scene.velocity(handle)
    );
}

#[test]
fn sleeping_body_wakes_only_when_the_fluid_can_lift_it() {
    let mut scene = Scene::new();
    let sinker = ball_of_density(&mut scene, 4000.0);
    let floater = scene
        .dynamic(
            250.0 * ball_volume(),
            Vec3::new(0.0, 0.0, 3.0),
            sphere(0.5),
            0.5,
            0.0,
        )
        .0;
    for handle in [sinker, floater] {
        if let Some(body) = scene.bodies.get_mut(handle) {
            body.state = BodyState::Sleeping;
        }
    }

    // `time_before_sleep` is huge and thresholds default so that the
    // sleeping flags are not reset by the sleep check.
    let config = PhysicsConfig {
        time_before_sleep: 1.0e9,
        ..PhysicsConfig::default()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    pipeline.set_sph_rigid_coupling(pool(1000.0));
    scene.step(&mut pipeline);

    let sinker_body = scene.bodies.get(sinker).unwrap();
    assert_eq!(sinker_body.state, BodyState::Sleeping, "buoyancy < weight");
    assert_eq!(sinker_body.velocity, Vec3::zeros());

    let floater_body = scene.bodies.get(floater).unwrap();
    assert_eq!(floater_body.state, BodyState::Active, "buoyancy > weight");
    assert!(floater_body.velocity.y > 0.0);
}

#[test]
fn coupling_stage_is_timed_and_reported() {
    let mut pipeline = coupled_pipeline(pool(1000.0));
    pipeline.enable_profiling(Some(256));
    let mut scene = Scene::new();
    ball_of_density(&mut scene, 1000.0);

    for _ in 0..20 {
        scene.step(&mut pipeline);
    }
    assert!(pipeline.stats.timings.coupling > Duration::ZERO);
    assert!(pipeline.stats.avg_stage_time("coupling") > Duration::ZERO);
    let profiles = pipeline.drain_profiles();
    assert!(
        profiles.iter().any(|p| p.stage == "coupling"),
        "the coupling stage must appear in the profile history"
    );
}

#[test]
fn builder_installs_the_coupling() {
    let params = pool(1000.0);
    let mut pipeline = PhysicsPipeline::builder()
        .sph_rigid_coupling(params)
        .finish();
    assert!(pipeline.sph_rigid_coupling().is_some());

    let mut scene = Scene::new();
    let ball = ball_of_density(&mut scene, 2000.0);
    // Sleeping is irrelevant for one step.
    scene.step(&mut pipeline);
    // Net acceleration (0.5 g - g) = -0.5 g.
    assert!((scene.velocity(ball).y + 0.5 * G * DT).abs() < 1e-9);
}

#[test]
fn sph_rigid_coupling_params_default() {
    let params = SphRigidCouplingParams::default();
    assert!((params.fluid_density - 1000.0).abs() < 1e-12);
    assert!((params.coupling_alpha - 1.0).abs() < 1e-12);
}

#[test]
fn standalone_buoyancy_helper_pushes_bodies_inside_the_region_up() {
    let coupling = SphRigidCouplingParams {
        fluid_density: 1000.0,
        coupling_alpha: 1.0,
        region_half_extents: Vec3::new(10.0, 10.0, 10.0),
        region_centre: Vec3::zeros(),
    };

    let mut bodies = RigidBodySet::new();
    let mut inside = RigidBody::new(1.0);
    inside.transform = Transform::from_position(Vec3::zeros());
    let inside = bodies.insert(inside);
    let mut outside = RigidBody::new(1.0);
    outside.transform = Transform::from_position(Vec3::new(0.0, 50.0, 0.0));
    let outside = bodies.insert(outside);

    apply_buoyancy_impulse(&mut bodies, &coupling, 0.5, G, DT);

    // F = rho g (4/3 pi r^3) for a unit mass over one step.
    let expected = 1000.0 * G * ball_volume() * DT;
    let vy = bodies.get(inside).unwrap().velocity.y;
    assert!(
        (vy - expected).abs() < 1e-9,
        "vy = {vy}, expected {expected}"
    );
    assert_eq!(bodies.get(outside).unwrap().velocity, Vec3::zeros());
}
