// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests: input validation of the engine API.
//!
//! Everything that crosses the WASM boundary is untrusted.  These tests check
//! that invalid masses, dimensions, handles, plane normals, time steps and
//! configurations are rejected with a descriptive [`Error`] instead of being
//! stored silently (a negative mass that made a body immovable, a negative
//! fixed step that integrated backwards, a collider on a missing body, ...).

#![cfg(test)]

use super::handles::decode_handle;
use super::*;
use crate::types::ColliderShapeType;

/// An engine without gravity, so that bodies stay where they are put.
fn quiet_engine() -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(0.0, 0.0, 0.0)
}

/// `true` when `result` failed with `InvalidParameter` naming `expected`.
fn is_invalid_param<T>(result: &Result<T>, expected: &str) -> bool {
    matches!(result, Err(Error::InvalidParameter { name, .. }) if name == expected)
}

// ---------------------------------------------------------------------------
// Bodies
// ---------------------------------------------------------------------------

#[test]
fn dynamic_body_rejects_non_positive_mass() {
    let mut engine = quiet_engine();
    // A zero mass used to silently create a static body, a negative mass a
    // body that was neither static nor dynamic.
    for mass in [-1.0, -1e-9, 0.0] {
        assert_eq!(
            engine.add_dynamic_body(mass, 0.0, 0.0, 0.0),
            Err(Error::InvalidMass(mass)),
            "mass {mass}"
        );
    }
    for mass in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            matches!(
                engine.add_dynamic_body(mass, 0.0, 0.0, 0.0),
                Err(Error::InvalidMass(m)) if m.is_nan() == mass.is_nan() && !m.is_finite()
            ),
            "mass {mass}"
        );
    }
    assert_eq!(engine.get_body_count(), 0, "nothing may have been stored");
    assert!(engine.bodies.is_empty());
}

#[test]
fn bodies_reject_non_finite_positions() {
    let mut engine = quiet_engine();
    assert!(is_invalid_param(
        &engine.add_dynamic_body(1.0, f64::NAN, 0.0, 0.0),
        "position"
    ));
    assert!(is_invalid_param(
        &engine.add_dynamic_body(1.0, 0.0, f64::INFINITY, 0.0),
        "position"
    ));
    assert!(is_invalid_param(
        &engine.add_static_body(0.0, 0.0, f64::NEG_INFINITY),
        "position"
    ));
    assert!(engine.bodies.is_empty());
}

#[test]
fn setters_reject_non_finite_vectors_and_leave_the_body_untouched() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(2.0, 1.0, 2.0, 3.0).expect("body");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("linear damping");
    engine
        .set_body_angular_damping(body, 0.0)
        .expect("angular damping");
    engine.set_velocity(body, 4.0, 5.0, 6.0).expect("velocity");
    engine
        .set_angular_velocity(body, 0.1, 0.2, 0.3)
        .expect("angular velocity");

    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(is_invalid_param(
            &engine.set_velocity(body, bad, 0.0, 0.0),
            "velocity"
        ));
        assert!(is_invalid_param(
            &engine.set_velocity(body, 0.0, bad, 0.0),
            "velocity"
        ));
        assert!(is_invalid_param(
            &engine.set_velocity(body, 0.0, 0.0, bad),
            "velocity"
        ));
        assert!(is_invalid_param(
            &engine.set_angular_velocity(body, bad, 0.0, 0.0),
            "angular_velocity"
        ));
        assert!(is_invalid_param(
            &engine.set_position(body, 0.0, bad, 0.0),
            "position"
        ));
        assert!(is_invalid_param(
            &engine.apply_force(body, 0.0, 0.0, bad),
            "force"
        ));
        assert!(is_invalid_param(
            &engine.apply_torque(body, bad, 0.0, 0.0),
            "torque"
        ));
        assert!(is_invalid_param(
            &engine.apply_impulse(body, bad, 0.0, 0.0),
            "impulse"
        ));
        assert!(is_invalid_param(
            &engine.apply_impulse_at_point(body, bad, 0.0, 0.0, 0.0, 0.0, 0.0),
            "impulse"
        ));
        assert!(is_invalid_param(
            &engine.apply_impulse_at_point(body, 0.0, 0.0, 0.0, 0.0, bad, 0.0),
            "offset"
        ));
    }

    // None of the rejected calls touched the body or queued a force: a step
    // moves it exactly by its (unchanged) velocity.
    assert_eq!(engine.get_position(body), [1.0, 2.0, 3.0]);
    assert_eq!(engine.get_velocity(body), [4.0, 5.0, 6.0]);
    assert_eq!(engine.get_angular_velocity(body), [0.1, 0.2, 0.3]);
    engine.step(1.0 / 60.0).expect("step");
    assert_eq!(engine.get_velocity(body), [4.0, 5.0, 6.0]);
    assert_eq!(engine.get_angular_velocity(body), [0.1, 0.2, 0.3]);
    let position = engine.get_position(body);
    assert!(position.iter().all(|c| c.is_finite()));
    assert!((position[0] - (1.0 + 4.0 / 60.0)).abs() < 1e-12);
}

#[test]
fn setters_report_a_dead_handle_before_looking_at_the_values() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.remove_body(body).expect("remove");
    assert_eq!(
        engine.set_velocity(body, f64::NAN, 0.0, 0.0),
        Err(Error::InvalidHandle(body))
    );
    assert_eq!(
        engine.set_body_linear_damping(body, f64::NAN),
        Err(Error::InvalidHandle(body))
    );
    assert_eq!(
        engine.apply_impulse_at_point(body, f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0),
        Err(Error::InvalidHandle(body))
    );
}

#[test]
fn damping_setters_reject_negative_and_non_finite_values() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.set_velocity(body, 1.0, 0.0, 0.0).expect("velocity");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("zero damping is allowed");
    for bad in [-0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(is_invalid_param(
            &engine.set_body_linear_damping(body, bad),
            "linear_damping"
        ));
        assert!(is_invalid_param(
            &engine.set_body_angular_damping(body, bad),
            "angular_damping"
        ));
    }
    // A rejected negative damping must not have been stored: it would make the
    // body accelerate on its own.
    engine.step(1.0 / 60.0).expect("step");
    assert_eq!(engine.get_velocity(body), [1.0, 0.0, 0.0]);

    engine
        .set_body_linear_damping(body, 600.0)
        .expect("a strong damping is valid");
    engine.step(1.0 / 60.0).expect("step");
    assert_eq!(
        engine.get_velocity(body),
        [0.0; 3],
        "damping above 1/dt stops the body"
    );
}

#[test]
fn rigid_body_config_is_validated() {
    let mut engine = quiet_engine();

    let mut cfg = RigidBodyConfig::new();
    cfg.mass = -2.0;
    assert_eq!(engine.add_rigid_body(&cfg), Err(Error::InvalidMass(-2.0)));

    let mut cfg = RigidBodyConfig::new();
    cfg.restitution = 1.5;
    assert_eq!(
        engine.add_rigid_body(&cfg),
        Err(Error::InvalidRestitution(1.5))
    );
    cfg.restitution = -0.1;
    assert_eq!(
        engine.add_rigid_body(&cfg),
        Err(Error::InvalidRestitution(-0.1))
    );

    let mut cfg = RigidBodyConfig::new();
    cfg.friction = -1.0;
    assert_eq!(
        engine.add_rigid_body(&cfg),
        Err(Error::InvalidFriction(-1.0))
    );

    let mut cfg = RigidBodyConfig::new();
    cfg.linear_velocity = [0.0, f64::NAN, 0.0];
    assert!(is_invalid_param(
        &engine.add_rigid_body(&cfg),
        "linear_velocity"
    ));

    let mut cfg = RigidBodyConfig::new();
    cfg.angular_velocity = [f64::INFINITY, 0.0, 0.0];
    assert!(is_invalid_param(
        &engine.add_rigid_body(&cfg),
        "angular_velocity"
    ));

    let mut cfg = RigidBodyConfig::new();
    cfg.rotation = [0.0; 4];
    assert!(is_invalid_param(&engine.add_rigid_body(&cfg), "rotation"));
    cfg.rotation = [0.0, f64::NAN, 0.0, 1.0];
    assert!(is_invalid_param(&engine.add_rigid_body(&cfg), "rotation"));

    let mut cfg = RigidBodyConfig::new();
    cfg.linear_damping = -0.5;
    assert!(is_invalid_param(
        &engine.add_rigid_body(&cfg),
        "linear_damping"
    ));
    let mut cfg = RigidBodyConfig::new();
    cfg.angular_damping = f64::NAN;
    assert!(is_invalid_param(
        &engine.add_rigid_body(&cfg),
        "angular_damping"
    ));

    let mut cfg = RigidBodyConfig::new();
    cfg.gravity_scale = f64::NAN;
    assert!(is_invalid_param(
        &engine.add_rigid_body(&cfg),
        "gravity_scale"
    ));

    assert!(
        engine.bodies.is_empty(),
        "rejected bodies must not be stored"
    );
}

#[test]
fn rigid_body_rotation_is_normalised() {
    let mut engine = quiet_engine();
    let mut cfg = RigidBodyConfig::new();
    cfg.rotation = [0.0, 0.0, 0.0, 2.0];
    let handle = engine.add_rigid_body(&cfg).expect("body");
    let q = engine.get_rotation(handle);
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-12,
        "stored rotation is a unit quaternion"
    );
    assert!((q[3] - 1.0).abs() < 1e-12);
}

#[test]
fn zero_mass_rigid_body_is_static_not_dynamic() {
    let mut engine = quiet_engine();
    let handle = engine
        .add_rigid_body(&RigidBodyConfig::static_body())
        .expect("static body");
    assert_eq!(engine.body_inv_mass(handle), 0.0);
    assert!(!engine.body_is_dynamic(handle));

    let kinematic = engine
        .add_rigid_body(&RigidBodyConfig::new().as_kinematic())
        .expect("kinematic body");
    assert!(!engine.body_is_dynamic(kinematic));

    let dynamic = engine
        .add_dynamic_body(2.0, 0.0, 0.0, 0.0)
        .expect("dynamic");
    assert!(engine.body_is_dynamic(dynamic));
    assert!((engine.body_inv_mass(dynamic) - 0.5).abs() < 1e-15);
}

// ---------------------------------------------------------------------------
// Colliders
// ---------------------------------------------------------------------------

#[test]
fn collider_requires_a_live_body() {
    let mut engine = quiet_engine();
    assert_eq!(
        engine.add_sphere_collider(0, 0.5),
        Err(Error::BodyNotFound(0))
    );
    assert_eq!(
        engine.add_collider(77, &ColliderConfig::sphere(0.5)),
        Err(Error::BodyNotFound(77))
    );

    // A removed body is rejected too, even when its slot has been reused.
    let old = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("old body");
    engine.remove_body(old).expect("remove old body");
    assert_eq!(
        engine.add_box_collider(old, 1.0, 1.0, 1.0),
        Err(Error::BodyNotFound(old))
    );
    let new = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("new body");
    assert_eq!(decode_handle(new).0, decode_handle(old).0);
    assert_eq!(
        engine.add_capsule_collider(old, 0.5, 1.0),
        Err(Error::BodyNotFound(old))
    );
    assert_eq!(engine.get_active_collider_count(), 0);
    assert!(engine.colliders.is_empty());
}

#[test]
fn sphere_radius_must_be_finite_and_positive() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    for radius in [0.0, -0.5, f64::NAN, f64::INFINITY] {
        let result = engine.add_sphere_collider(body, radius);
        assert!(
            matches!(result, Err(Error::InvalidDimension(r)) if r.is_nan() == radius.is_nan()),
            "radius {radius}: {result:?}"
        );
    }
    assert!(engine.colliders.is_empty());
}

#[test]
fn box_half_extents_must_be_finite_and_positive() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    for (hx, hy, hz) in [
        (0.0, 1.0, 1.0),
        (1.0, -2.0, 1.0),
        (1.0, 1.0, f64::NAN),
        (f64::INFINITY, 1.0, 1.0),
    ] {
        assert!(
            matches!(
                engine.add_box_collider(body, hx, hy, hz),
                Err(Error::InvalidDimension(_))
            ),
            "half extents ({hx}, {hy}, {hz})"
        );
    }
    assert!(engine.colliders.is_empty());
}

#[test]
fn capsule_and_cylinder_dimensions_must_be_finite_and_positive() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    for (radius, height) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (0.5, 0.0),
        (0.5, -1.0),
        (0.5, f64::NAN),
    ] {
        assert!(
            matches!(
                engine.add_capsule_collider(body, radius, height),
                Err(Error::InvalidDimension(_))
            ),
            "capsule ({radius}, {height})"
        );
        assert!(
            matches!(
                engine.add_cylinder_collider(body, radius, height),
                Err(Error::InvalidDimension(_))
            ),
            "cylinder ({radius}, {height})"
        );
    }
    assert!(engine.colliders.is_empty());
}

#[test]
fn plane_collider_requires_a_finite_nonzero_normal() {
    let mut engine = quiet_engine();
    let body = engine.add_static_body(0.0, 0.0, 0.0).expect("body");
    assert!(is_invalid_param(
        &engine.add_plane_collider(body, 0.0, 0.0, 0.0, 1.0),
        "plane_normal"
    ));
    assert!(is_invalid_param(
        &engine.add_plane_collider(body, f64::NAN, 1.0, 0.0, 0.0),
        "plane_normal"
    ));
    assert!(is_invalid_param(
        &engine.add_plane_collider(body, 0.0, f64::INFINITY, 0.0, 0.0),
        "plane_normal"
    ));
    assert!(is_invalid_param(
        &engine.add_plane_collider(body, 0.0, 1.0, 0.0, f64::NAN),
        "plane_offset"
    ));
    assert!(engine.colliders.is_empty());
}

#[test]
fn plane_config_normalises_the_normal_and_rescales_the_offset() {
    let cfg = ColliderConfig::plane(0.0, 2.0, 0.0, 4.0);
    assert_eq!(cfg.plane_normal, [0.0, 1.0, 0.0]);
    assert!(
        (cfg.plane_offset - 2.0).abs() < 1e-15,
        "offset 4 / |n| 2 = 2"
    );

    let cfg = ColliderConfig::plane(3.0, 0.0, 4.0, 10.0);
    assert!((cfg.plane_normal[0] - 0.6).abs() < 1e-15);
    assert!((cfg.plane_normal[2] - 0.8).abs() < 1e-15);
    assert!((cfg.plane_offset - 2.0).abs() < 1e-15);

    // A degenerate normal is kept as given so that adding it is rejected.
    let cfg = ColliderConfig::plane(0.0, 0.0, 0.0, 1.0);
    assert_eq!(cfg.plane_normal, [0.0; 3]);
}

#[test]
fn unnormalised_plane_describes_the_same_plane_as_its_unit_form() {
    let mut engine = quiet_engine();
    let ground = engine.add_static_body(0.0, 0.0, 0.0).expect("ground");
    // `2 y = 4` is the plane y = 2.
    engine
        .add_plane_collider(ground, 0.0, 2.0, 0.0, 4.0)
        .expect("plane");
    let hit = engine.raycast(0.0, 5.0, 0.0, 0.0, -1.0, 0.0, 100.0);
    assert!(hit.hit);
    assert!(
        (hit.distance - 3.0).abs() < 1e-12,
        "plane at y = 2, got {}",
        hit.distance
    );
    assert_eq!(hit.normal, [0.0, 1.0, 0.0]);
}

#[test]
fn convex_hull_colliders_are_rejected() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    let mut cfg = ColliderConfig::sphere(0.5);
    cfg.shape_type = ColliderShapeType::ConvexHull;
    assert!(matches!(
        engine.add_collider(body, &cfg),
        Err(Error::UnsupportedShape(_))
    ));
    assert!(engine.colliders.is_empty());
}

#[test]
fn collider_pose_and_density_are_validated() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");

    let mut cfg = ColliderConfig::sphere(0.5);
    cfg.local_position = [0.0, f64::NAN, 0.0];
    assert!(is_invalid_param(
        &engine.add_collider(body, &cfg),
        "local_position"
    ));

    let mut cfg = ColliderConfig::sphere(0.5);
    cfg.local_rotation = [0.0; 4];
    assert!(is_invalid_param(
        &engine.add_collider(body, &cfg),
        "local_rotation"
    ));

    for density in [0.0, -1.0, f64::NAN] {
        let mut cfg = ColliderConfig::sphere(0.5);
        cfg.density = density;
        assert!(
            is_invalid_param(&engine.add_collider(body, &cfg), "density"),
            "density {density}"
        );
    }
    assert!(engine.colliders.is_empty());
    assert!(engine.collider_free_list.is_empty());
}

#[test]
fn rejected_collider_leaves_the_body_inertia_untouched() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.add_sphere_collider(body, 1.0).expect("sphere");
    let before = engine.get_all_positions();
    assert!(engine.add_sphere_collider(body, -1.0).is_err());
    assert!(engine.add_plane_collider(body, 0.0, 0.0, 0.0, 0.0).is_err());
    assert_eq!(engine.get_active_collider_count(), 1);
    assert_eq!(engine.get_all_positions(), before);
    // The sphere's inertia is still in place: a torque still spins the body.
    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine.step(1.0 / 60.0).expect("step");
    assert!(engine.get_angular_velocity(body)[1] > 0.0);
}

// ---------------------------------------------------------------------------
// Time step and configuration
// ---------------------------------------------------------------------------

#[test]
fn set_fixed_dt_rejects_invalid_values() {
    let mut engine = quiet_engine();
    let default_dt = engine.config().fixed_dt;
    // A negative fixed step made every substep integrate backwards.
    assert_eq!(
        engine.set_fixed_dt(-0.01),
        Err(Error::InvalidTimeStep(-0.01))
    );
    assert_eq!(engine.set_fixed_dt(0.0), Err(Error::InvalidTimeStep(0.0)));
    for dt in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            matches!(engine.set_fixed_dt(dt), Err(Error::InvalidTimeStep(_))),
            "dt {dt}"
        );
    }
    assert_eq!(
        engine.config().fixed_dt,
        default_dt,
        "a rejected value is not stored"
    );

    engine.set_fixed_dt(0.01).expect("valid step");
    assert_eq!(engine.config().fixed_dt, 0.01);
}

#[test]
fn rejected_fixed_dt_does_not_reverse_the_simulation() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.set_velocity(body, 1.0, 0.0, 0.0).expect("velocity");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("linear damping");
    assert!(engine.set_fixed_dt(-0.01).is_err());
    engine.step(1.0 / 60.0).expect("step");
    assert!(
        engine.get_position(body)[0] > 0.0,
        "the body moves forward in time"
    );
    assert!(engine.time() > 0.0);
}

#[test]
fn step_rejects_invalid_time_steps_and_leaves_the_engine_untouched() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let body = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0).expect("body");
    assert_eq!(engine.step(-1.0), Err(Error::InvalidTimeStep(-1.0)));
    for dt in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            matches!(engine.step(dt), Err(Error::InvalidTimeStep(_))),
            "dt {dt}"
        );
    }
    assert_eq!(engine.time(), 0.0);
    assert_eq!(engine.interpolation_alpha(), 0.0);
    assert_eq!(engine.get_position(body), [0.0, 10.0, 0.0]);
    assert_eq!(engine.get_velocity(body), [0.0; 3]);
}

#[test]
fn from_config_accepts_valid_configurations() {
    let engine = WasmPhysicsEngine::from_config(SimulationConfig::default()).expect("default");
    assert_eq!(engine.gravity(), [0.0, -9.81, 0.0]);

    let custom = SimulationConfig {
        gravity: [0.0, -1.62, 0.0],
        fixed_dt: 0.01,
        max_substeps: 8,
        ..SimulationConfig::default()
    };
    let engine = WasmPhysicsEngine::from_config(custom).expect("custom");
    assert_eq!(engine.gravity(), [0.0, -1.62, 0.0]);
    assert_eq!(engine.config().fixed_dt, 0.01);
    assert_eq!(engine.config().max_substeps, 8);

    for preset in [
        SimulationConfig::high_accuracy(),
        SimulationConfig::performance(),
    ] {
        assert!(preset.validate().is_ok());
        assert!(WasmPhysicsEngine::from_config(preset).is_ok());
    }
}

#[test]
fn from_config_rejects_invalid_configurations() {
    let bad_dt = SimulationConfig {
        fixed_dt: -1.0,
        ..SimulationConfig::default()
    };
    assert!(matches!(
        WasmPhysicsEngine::from_config(bad_dt),
        Err(Error::InvalidTimeStep(dt)) if dt == -1.0
    ));

    for fixed_dt in [0.0, f64::NAN, f64::INFINITY] {
        let cfg = SimulationConfig {
            fixed_dt,
            ..SimulationConfig::default()
        };
        assert!(
            matches!(cfg.validate(), Err(Error::InvalidTimeStep(_))),
            "fixed_dt {fixed_dt}"
        );
    }

    let no_substeps = SimulationConfig {
        max_substeps: 0,
        ..SimulationConfig::default()
    };
    assert!(is_invalid_param(&no_substeps.validate(), "max_substeps"));

    let bad_gravity = SimulationConfig {
        gravity: [0.0, f64::NAN, 0.0],
        ..SimulationConfig::default()
    };
    assert!(is_invalid_param(&bad_gravity.validate(), "gravity"));

    let negative_slop = SimulationConfig {
        contact_slop: -0.1,
        ..SimulationConfig::default()
    };
    assert!(is_invalid_param(&negative_slop.validate(), "contact_slop"));

    let bad_sleep = SimulationConfig {
        time_before_sleep: f64::NAN,
        ..SimulationConfig::default()
    };
    assert!(is_invalid_param(&bad_sleep.validate(), "time_before_sleep"));
    assert!(WasmPhysicsEngine::from_config(bad_sleep).is_err());
}

// ---------------------------------------------------------------------------
// Atomic body + collider creation
// ---------------------------------------------------------------------------

#[test]
fn add_body_with_collider_creates_both() {
    let mut engine = quiet_engine();
    let body_cfg = RigidBodyConfig::new()
        .with_mass(2.0)
        .with_position(1.0, 2.0, 3.0);
    let (body, collider) = engine
        .add_body_with_collider(&body_cfg, &ColliderConfig::sphere(0.5))
        .expect("body with collider");
    assert_eq!(body, 0);
    assert_eq!(collider, 0);
    assert_eq!(engine.get_position(body), [1.0, 2.0, 3.0]);
    assert_eq!(engine.get_active_collider_count(), 1);
}

#[test]
fn add_body_with_collider_is_all_or_nothing() {
    let mut engine = quiet_engine();

    let body_cfg = RigidBodyConfig::new();
    assert!(matches!(
        engine.add_body_with_collider(&body_cfg, &ColliderConfig::sphere(-1.0)),
        Err(Error::InvalidDimension(_))
    ));
    assert!(is_invalid_param(
        &engine.add_body_with_collider(&body_cfg, &ColliderConfig::plane(0.0, 0.0, 0.0, 0.0)),
        "plane_normal"
    ));

    let mut bad_body = RigidBodyConfig::new();
    bad_body.mass = f64::NAN;
    assert!(matches!(
        engine.add_body_with_collider(&bad_body, &ColliderConfig::sphere(0.5)),
        Err(Error::InvalidMass(_))
    ));

    assert_eq!(engine.get_body_count(), 0, "no body may be left behind");
    assert!(engine.bodies.is_empty(), "not even a deactivated slot");
    assert!(engine.colliders.is_empty());
}
