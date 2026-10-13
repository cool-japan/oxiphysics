// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests: rotational dynamics and fixed-step time integration.
//!
//! * torques and off-centre impulses are mapped through the inertia tensor
//!   that follows from the attached colliders;
//! * the orientation follows the world-frame angular velocity;
//! * `step()` consumes real time in fixed substeps and the simulation clock
//!   advances by the time that was really integrated.

#![cfg(test)]

use std::f64::consts::{FRAC_PI_2, PI};

use super::geom::{
    Mat3, Vec3, length, mat_mul, mat_transpose, mat_vec, quat_from_axis_angle, quat_mul,
    quat_to_mat3,
};
use super::inertia::{ShapePart, compound_inertia, principal_moments, volume};
use super::shapes::Shape;
use super::*;

const DT: f64 = 1.0 / 60.0;

/// Engine without gravity whose substep length is `fixed_dt`.
fn free_engine(fixed_dt: f64, max_substeps: u32) -> WasmPhysicsEngine {
    WasmPhysicsEngine::from_config(SimulationConfig {
        gravity: [0.0; 3],
        fixed_dt,
        max_substeps,
        ..SimulationConfig::default()
    })
    .expect("valid configuration")
}

/// A body of the given `mass` with a collider and without any damping.
fn damping_free_body(engine: &mut WasmPhysicsEngine, mass: f64, collider: ColliderConfig) -> u32 {
    let body = engine
        .add_dynamic_body(mass, 0.0, 0.0, 0.0)
        .expect("dynamic body");
    engine.add_collider(body, &collider).expect("collider");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("linear damping");
    engine
        .set_body_angular_damping(body, 0.0)
        .expect("angular damping");
    body
}

/// A body without a collider and without any damping, moving with `velocity`.
fn drifting_body(engine: &mut WasmPhysicsEngine, velocity: Vec3) -> u32 {
    let body = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("dynamic body");
    engine
        .set_velocity(body, velocity[0], velocity[1], velocity[2])
        .expect("velocity");
    engine
        .set_body_linear_damping(body, 0.0)
        .expect("linear damping");
    engine
        .set_body_angular_damping(body, 0.0)
        .expect("angular damping");
    body
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: expected {expected}, got {actual} (tolerance {tolerance})"
    );
}

fn assert_vec_close(actual: Vec3, expected: Vec3, tolerance: f64, what: &str) {
    for axis in 0..3 {
        assert_close(
            actual[axis],
            expected[axis],
            tolerance,
            &format!("{what}[{axis}]"),
        );
    }
}

// ---------------------------------------------------------------------------
// Mass properties
// ---------------------------------------------------------------------------

#[test]
fn sphere_inertia_is_two_fifths_m_r_squared() {
    let moments = principal_moments(&Shape::Sphere { radius: 1.5 }, 3.0).expect("sphere inertia");
    let expected = 0.4 * 3.0 * 1.5 * 1.5;
    assert_vec_close(moments, [expected; 3], 1e-12, "sphere");
}

#[test]
fn box_inertia_follows_the_half_extents() {
    // Half-extents (1, 0.5, 0.25), mass 12: Ixx = m/3 (hy^2 + hz^2) = 1.25, ...
    let shape = Shape::Box {
        half_extents: [1.0, 0.5, 0.25],
    };
    let moments = principal_moments(&shape, 12.0).expect("box inertia");
    assert_vec_close(moments, [1.25, 4.25, 5.0], 1e-12, "box");
}

#[test]
fn cylinder_inertia_uses_the_y_axis_as_symmetry_axis() {
    let shape = Shape::Cylinder {
        radius: 0.5,
        half_height: 1.0,
    };
    let moments = principal_moments(&shape, 2.0).expect("cylinder inertia");
    let axial = 0.5 * 2.0 * 0.25;
    let transverse = 2.0 * (3.0 * 0.25 + 4.0) / 12.0;
    assert_vec_close(moments, [transverse, axial, transverse], 1e-12, "cylinder");
}

#[test]
fn planes_have_neither_volume_nor_inertia() {
    let plane = Shape::Plane {
        normal: [0.0, 1.0, 0.0],
        offset: 0.0,
    };
    assert_eq!(volume(&plane), 0.0);
    assert!(principal_moments(&plane, 1.0).is_none());
    let parts = [ShapePart {
        shape: plane,
        position: [0.0; 3],
        rotation: quat_to_mat3([0.0, 0.0, 0.0, 1.0]),
    }];
    assert!(compound_inertia(1.0, &parts).is_none());
}

/// Moments of inertia about the coordinate axes of a solid of uniform density
/// by midpoint integration over a voxel grid; `inside` tests membership.
fn numeric_moments(extent: Vec3, cells: usize, inside: impl Fn(Vec3) -> bool, mass: f64) -> Vec3 {
    let step = [
        2.0 * extent[0] / cells as f64,
        2.0 * extent[1] / cells as f64,
        2.0 * extent[2] / cells as f64,
    ];
    let cell_volume = step[0] * step[1] * step[2];
    let mut volume_sum = 0.0;
    let mut moments = [0.0; 3];
    for i in 0..cells {
        for j in 0..cells {
            for k in 0..cells {
                let p = [
                    -extent[0] + (i as f64 + 0.5) * step[0],
                    -extent[1] + (j as f64 + 0.5) * step[1],
                    -extent[2] + (k as f64 + 0.5) * step[2],
                ];
                if inside(p) {
                    volume_sum += cell_volume;
                    moments[0] += (p[1] * p[1] + p[2] * p[2]) * cell_volume;
                    moments[1] += (p[0] * p[0] + p[2] * p[2]) * cell_volume;
                    moments[2] += (p[0] * p[0] + p[1] * p[1]) * cell_volume;
                }
            }
        }
    }
    let density = mass / volume_sum;
    [
        moments[0] * density,
        moments[1] * density,
        moments[2] * density,
    ]
}

#[test]
fn capsule_inertia_matches_numerical_integration() {
    let (radius, half_height, mass) = (0.5, 1.0, 3.0);
    let shape = Shape::Capsule {
        radius,
        half_height,
    };
    let analytic = principal_moments(&shape, mass).expect("capsule inertia");
    let numeric = numeric_moments(
        [radius, half_height + radius, radius],
        120,
        |p| {
            let along = (p[1].abs() - half_height).max(0.0);
            p[0] * p[0] + p[2] * p[2] + along * along <= radius * radius
        },
        mass,
    );
    for axis in 0..3 {
        let relative = (analytic[axis] - numeric[axis]).abs() / analytic[axis];
        assert!(
            relative < 0.01,
            "capsule I[{axis}]: analytic {} vs numeric {} ({:.3} %)",
            analytic[axis],
            numeric[axis],
            100.0 * relative
        );
    }
}

#[test]
fn cylinder_inertia_matches_numerical_integration() {
    let (radius, half_height, mass) = (0.75, 0.5, 2.0);
    let analytic = principal_moments(
        &Shape::Cylinder {
            radius,
            half_height,
        },
        mass,
    )
    .expect("cylinder inertia");
    let numeric = numeric_moments(
        [radius, half_height, radius],
        120,
        |p| p[0] * p[0] + p[2] * p[2] <= radius * radius,
        mass,
    );
    for axis in 0..3 {
        let relative = (analytic[axis] - numeric[axis]).abs() / analytic[axis];
        assert!(
            relative < 0.01,
            "cylinder I[{axis}] off by {:.3} %",
            100.0 * relative
        );
    }
}

#[test]
fn off_centre_collider_adds_the_parallel_axis_term() {
    // Sphere r = 0.5 (own inertia 0.1) displaced by d = 2 along x.
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(
        &mut engine,
        1.0,
        ColliderConfig {
            local_position: [2.0, 0.0, 0.0],
            ..ColliderConfig::sphere(0.5)
        },
    );
    engine
        .apply_torque(body, 1.0, 0.0, 0.0)
        .expect("torque about x");
    engine.step(DT).expect("step");
    // Axis through the offset sphere's own centre direction: no parallel term.
    assert_close(engine.get_angular_velocity(body)[0], DT / 0.1, 1e-9, "wx");

    engine
        .set_angular_velocity(body, 0.0, 0.0, 0.0)
        .expect("reset spin");
    engine
        .apply_torque(body, 0.0, 1.0, 0.0)
        .expect("torque about y");
    engine.step(DT).expect("step");
    // Iyy = 0.1 + m d^2 = 4.1.
    assert_close(engine.get_angular_velocity(body)[1], DT / 4.1, 1e-9, "wy");
}

#[test]
fn compound_inertia_splits_the_mass_by_volume() {
    // Two identical spheres at x = +-1 share the mass of 2 kg equally.
    let mut engine = free_engine(DT, 4);
    let body = engine
        .add_dynamic_body(2.0, 0.0, 0.0, 0.0)
        .expect("dynamic body");
    for x in [-1.0, 1.0] {
        engine
            .add_collider(
                body,
                &ColliderConfig {
                    local_position: [x, 0.0, 0.0],
                    ..ColliderConfig::sphere(0.5)
                },
            )
            .expect("sphere");
    }
    engine
        .set_body_angular_damping(body, 0.0)
        .expect("angular damping");
    engine.apply_torque(body, 1.0, 0.0, 0.0).expect("torque");
    engine.step(DT).expect("step");
    // Ixx = 2 * 0.4 * 1 * 0.25 = 0.2.
    assert_close(engine.get_angular_velocity(body)[0], DT / 0.2, 1e-9, "wx");
}

// ---------------------------------------------------------------------------
// Torque (F5655)
// ---------------------------------------------------------------------------

#[test]
fn torque_spins_a_unit_sphere_by_the_inverse_inertia() {
    // tau_y = 1 on a 1 kg, r = 1 sphere (I = 0.4) for one 1/60 s step.
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine.step(DT).expect("step");

    let omega = engine.get_angular_velocity(body);
    assert_vec_close(omega, [0.0, DT / 0.4, 0.0], 1e-12, "angular velocity");

    // The orientation advanced by omega * dt about +y.
    let angle = omega[1] * DT;
    let q = engine.get_rotation(body);
    assert_vec_close(
        [q[0], q[1], q[2]],
        [0.0, (0.5 * angle).sin(), 0.0],
        1e-12,
        "rotation axis part",
    );
    assert_close(q[3], (0.5 * angle).cos(), 1e-12, "rotation w");
}

#[test]
fn torque_acts_during_every_substep_of_a_step() {
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine.step(3.0 * DT).expect("step");
    assert_close(
        engine.get_angular_velocity(body)[1],
        3.0 * DT / 0.4,
        1e-9,
        "three substeps of torque",
    );
}

#[test]
fn torque_is_cleared_after_the_step() {
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine.step(DT).expect("first step");
    let after_first = engine.get_angular_velocity(body);
    engine.step(DT).expect("second step");
    assert_vec_close(
        engine.get_angular_velocity(body),
        after_first,
        1e-15,
        "no new torque, no new spin",
    );
}

#[test]
fn torques_accumulate_before_a_step() {
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine.apply_torque(body, 0.0, 0.25, 0.0).expect("torque a");
    engine.apply_torque(body, 0.0, 0.75, 0.0).expect("torque b");
    engine.step(DT).expect("step");
    assert_close(engine.get_angular_velocity(body)[1], DT / 0.4, 1e-12, "wy");
}

#[test]
fn kinetic_energy_includes_the_rotational_part() {
    // One step of unit torque on a unit sphere (I = 0.4) leaves
    // w = DT / 0.4, i.e. the rotational energy I w^2 / 2.
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    let at_rest = engine.get_body_state(body).expect("state at rest");
    assert_eq!(at_rest.kinetic_energy, 0.0);

    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine.step(DT).expect("step");
    let omega = DT / 0.4;
    let spinning = engine.get_body_state(body).expect("spinning state");
    assert_close(
        spinning.kinetic_energy,
        0.5 * 0.4 * omega * omega,
        1e-15,
        "rotational energy",
    );

    // The translational and the rotational part add up.
    engine
        .set_velocity(body, 3.0, 0.0, 0.0)
        .expect("set velocity");
    let moving = engine.get_body_state(body).expect("moving state");
    assert_close(
        moving.kinetic_energy,
        0.5 * 9.0 + 0.5 * 0.4 * omega * omega,
        1e-12,
        "total energy",
    );
}

#[test]
fn rotational_energy_follows_the_inertia_about_the_spin_axis_in_the_world_frame() {
    // Half-extents (1, 0.5, 0.25), mass 12: I = diag(1.25, 4.25, 5.0).
    let mut engine = free_engine(DT, 4);
    let upright = damping_free_body(&mut engine, 12.0, ColliderConfig::cuboid(1.0, 0.5, 0.25));
    engine
        .set_angular_velocity(upright, 0.0, 2.0, 0.0)
        .expect("spin about y");
    assert_close(
        engine
            .get_body_state(upright)
            .expect("upright state")
            .kinetic_energy,
        0.5 * 4.25 * 4.0,
        1e-12,
        "spin about the body y axis",
    );

    // Turned by 90 degrees about z, the world y axis is the body x axis.
    let turned = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation: quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new().with_mass(12.0)
        })
        .expect("turned body");
    engine
        .add_box_collider(turned, 1.0, 0.5, 0.25)
        .expect("box collider");
    engine
        .set_angular_velocity(turned, 0.0, 2.0, 0.0)
        .expect("spin about world y");
    assert_close(
        engine
            .get_body_state(turned)
            .expect("turned state")
            .kinetic_energy,
        0.5 * 1.25 * 4.0,
        1e-9,
        "spin about the body x axis",
    );
}

#[test]
fn spinning_bodies_without_rotational_inertia_have_no_rotational_energy() {
    let mut engine = free_engine(DT, 4);
    let bare = engine
        .add_dynamic_body(2.0, 0.0, 0.0, 0.0)
        .expect("body without collider");
    engine
        .set_angular_velocity(bare, 5.0, 0.0, 0.0)
        .expect("spin");
    assert_eq!(
        engine
            .get_body_state(bare)
            .expect("bare state")
            .kinetic_energy,
        0.0
    );
}

#[test]
fn box_resists_torque_differently_about_each_axis() {
    // Half-extents (1, 0.5, 0.25), mass 12: I = diag(1.25, 4.25, 5.0).
    let mut engine = free_engine(DT, 4);
    let mut bodies = Vec::new();
    for index in 0..3 {
        let body = damping_free_body(&mut engine, 12.0, ColliderConfig::cuboid(1.0, 0.5, 0.25));
        // Keep the three boxes apart so that they do not touch each other.
        engine
            .set_position(body, 10.0 * index as f64, 0.0, 0.0)
            .expect("position");
        bodies.push(body);
    }
    engine
        .apply_torque(bodies[0], 1.0, 0.0, 0.0)
        .expect("torque x");
    engine
        .apply_torque(bodies[1], 0.0, 1.0, 0.0)
        .expect("torque y");
    engine
        .apply_torque(bodies[2], 0.0, 0.0, 1.0)
        .expect("torque z");
    engine.step(DT).expect("step");

    assert_vec_close(
        engine.get_angular_velocity(bodies[0]),
        [DT / 1.25, 0.0, 0.0],
        1e-12,
        "spin about x",
    );
    assert_vec_close(
        engine.get_angular_velocity(bodies[1]),
        [0.0, DT / 4.25, 0.0],
        1e-12,
        "spin about y",
    );
    assert_vec_close(
        engine.get_angular_velocity(bodies[2]),
        [0.0, 0.0, DT / 5.0],
        1e-12,
        "spin about z",
    );
}

#[test]
fn torque_uses_the_world_space_inertia_of_a_rotated_body() {
    // A box turned by 90 degrees about z swaps its x and y inertia in the
    // world frame: a world-x torque now meets Iyy = 4.25 of the body frame.
    let mut engine = free_engine(DT, 4);
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation: quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new().with_mass(12.0)
        })
        .expect("rotated body");
    engine
        .add_box_collider(body, 1.0, 0.5, 0.25)
        .expect("box collider");
    engine.apply_torque(body, 1.0, 0.0, 0.0).expect("torque");
    engine.step(DT).expect("step");
    assert_vec_close(
        engine.get_angular_velocity(body),
        [DT / 4.25, 0.0, 0.0],
        1e-9,
        "world-space response",
    );
}

#[test]
fn body_without_a_volume_collider_has_infinite_rotational_inertia() {
    let mut engine = free_engine(DT, 4);
    let bare = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("body without collider");
    let on_plane = engine
        .add_dynamic_body(1.0, 5.0, 0.0, 0.0)
        .expect("body with plane collider only");
    engine
        .add_plane_collider(on_plane, 0.0, 1.0, 0.0, -100.0)
        .expect("far away plane");
    engine.apply_torque(bare, 0.0, 1.0, 0.0).expect("torque");
    engine
        .apply_torque(on_plane, 0.0, 1.0, 0.0)
        .expect("torque");
    engine.step(DT).expect("step");
    assert_eq!(engine.get_angular_velocity(bare), [0.0; 3]);
    assert_eq!(engine.get_angular_velocity(on_plane), [0.0; 3]);
}

#[test]
fn static_bodies_ignore_torques_and_impulses() {
    let mut engine = free_engine(DT, 4);
    let body = engine.add_static_body(1.0, 2.0, 3.0).expect("static body");
    engine
        .add_sphere_collider(body, 1.0)
        .expect("static sphere collider");
    engine.apply_torque(body, 0.0, 1.0, 0.0).expect("torque");
    engine
        .apply_impulse_at_point(body, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
        .expect("impulse");
    engine.step(DT).expect("step");
    assert_eq!(engine.get_angular_velocity(body), [0.0; 3]);
    assert_eq!(engine.get_velocity(body), [0.0; 3]);
    assert_eq!(engine.get_position(body), [1.0, 2.0, 3.0]);
}

// ---------------------------------------------------------------------------
// Impulse at a point (F5655)
// ---------------------------------------------------------------------------

#[test]
fn off_centre_impulse_spins_a_sphere_by_the_inverse_inertia() {
    // J = (0, 0, 1) at r = (1, 0, 0): r x J = (0, -1, 0), I = 0.4.
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine
        .apply_impulse_at_point(body, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
        .expect("impulse");
    assert_vec_close(
        engine.get_velocity(body),
        [0.0, 0.0, 1.0],
        1e-12,
        "linear velocity",
    );
    assert_vec_close(
        engine.get_angular_velocity(body),
        [0.0, -2.5, 0.0],
        1e-12,
        "angular velocity",
    );
}

#[test]
fn impulse_through_the_centre_does_not_spin_the_body() {
    let mut engine = free_engine(DT, 4);
    let body = damping_free_body(&mut engine, 2.0, ColliderConfig::cuboid(1.0, 0.5, 0.25));
    engine
        .apply_impulse_at_point(body, 4.0, 0.0, 0.0, 3.0, 0.0, 0.0)
        .expect("impulse along the lever arm");
    assert_vec_close(
        engine.get_velocity(body),
        [2.0, 0.0, 0.0],
        1e-12,
        "velocity",
    );
    assert_eq!(engine.get_angular_velocity(body), [0.0; 3]);
}

#[test]
fn impulse_response_depends_on_the_axis_and_the_orientation() {
    let box_mass = 12.0;
    // J = (0, 0, 1) at r = (0, 1, 0): angular impulse (1, 0, 0).
    let mut engine = free_engine(DT, 4);
    let upright = damping_free_body(
        &mut engine,
        box_mass,
        ColliderConfig::cuboid(1.0, 0.5, 0.25),
    );
    engine
        .apply_impulse_at_point(upright, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0)
        .expect("impulse");
    assert_vec_close(
        engine.get_angular_velocity(upright),
        [1.0 / 1.25, 0.0, 0.0],
        1e-12,
        "upright box",
    );

    // The same impulse on a box turned by 90 degrees about z meets Iyy.
    let turned = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation: quat_from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2),
            ..RigidBodyConfig::new().with_mass(box_mass)
        })
        .expect("turned body");
    engine
        .add_box_collider(turned, 1.0, 0.5, 0.25)
        .expect("box collider");
    engine
        .apply_impulse_at_point(turned, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0)
        .expect("impulse");
    assert_vec_close(
        engine.get_angular_velocity(turned),
        [1.0 / 4.25, 0.0, 0.0],
        1e-9,
        "turned box",
    );
}

// ---------------------------------------------------------------------------
// Orientation and gyroscopic motion
// ---------------------------------------------------------------------------

#[test]
fn orientation_follows_the_world_frame_angular_velocity() {
    // 90 degrees about x, then spun about the *world* y axis by 90 degrees.
    let q0 = quat_from_axis_angle([1.0, 0.0, 0.0], FRAC_PI_2);
    let mut engine = free_engine(DT, 100);
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            rotation: q0,
            angular_velocity: [0.0, FRAC_PI_2, 0.0],
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new()
        })
        .expect("spinning body");
    engine.step(1.0).expect("one second");
    assert_close(engine.time(), 1.0, 1e-9, "integrated time");

    let spin = quat_from_axis_angle([0.0, 1.0, 0.0], FRAC_PI_2);
    let world_frame = quat_mul(spin, q0);
    let body_frame = quat_mul(q0, spin);
    let q = engine.get_rotation(body);
    // q and -q describe the same rotation.
    let overlap: f64 = (0..4).map(|i| q[i] * world_frame[i]).sum();
    let sign = if overlap < 0.0 { -1.0 } else { 1.0 };
    for i in 0..4 {
        assert_close(sign * q[i], world_frame[i], 1e-9, "world-frame composition");
    }
    let differs = (0..4).any(|i| (sign * q[i] - body_frame[i]).abs() > 0.1);
    assert!(
        differs,
        "the result must not equal the body-frame composition"
    );
}

#[test]
fn half_a_turn_per_second_brings_the_body_upside_down() {
    let mut engine = free_engine(DT, 100);
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            angular_velocity: [0.0, PI, 0.0],
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new()
        })
        .expect("spinning body");
    engine.step(1.0).expect("one second");
    let q = engine.get_rotation(body);
    assert_close(q[1].abs(), 1.0, 1e-9, "|qy|");
    assert_close(q[3], 0.0, 1e-9, "qw");
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    assert_close(norm, 1.0, 1e-12, "unit quaternion");
}

/// World-frame angular momentum `R I R^T w` of a body with body-frame inertia
/// `inertia`.
fn angular_momentum(rotation: [f64; 4], inertia: &Mat3, omega: Vec3) -> Vec3 {
    let r = quat_to_mat3(rotation);
    let world = mat_mul(&mat_mul(&r, inertia), &mat_transpose(&r));
    mat_vec(&world, omega)
}

#[test]
fn torque_free_asymmetric_body_conserves_angular_momentum() {
    // I = diag(1.25, 4.25, 5.0); spin mostly about the axis of largest inertia
    // with a small wobble, so the motion is a bounded precession.
    let inertia: Mat3 = [[1.25, 0.0, 0.0], [0.0, 4.25, 0.0], [0.0, 0.0, 5.0]];
    let mut engine = free_engine(1.0 / 240.0, 100);
    let body = engine
        .add_rigid_body(&RigidBodyConfig {
            angular_velocity: [0.3, 0.2, 3.0],
            linear_damping: 0.0,
            angular_damping: 0.0,
            ..RigidBodyConfig::new().with_mass(12.0)
        })
        .expect("spinning box");
    engine
        .add_box_collider(body, 1.0, 0.5, 0.25)
        .expect("box collider");

    let initial = angular_momentum(
        engine.get_rotation(body),
        &inertia,
        engine.get_angular_velocity(body),
    );
    let mut max_drift: f64 = 0.0;
    let mut direction_change: f64 = 0.0;
    let initial_omega = engine.get_angular_velocity(body);
    for _ in 0..480 {
        engine.step(1.0 / 240.0).expect("step");
        let l = angular_momentum(
            engine.get_rotation(body),
            &inertia,
            engine.get_angular_velocity(body),
        );
        max_drift = max_drift.max(length([
            l[0] - initial[0],
            l[1] - initial[1],
            l[2] - initial[2],
        ]));
        let omega = engine.get_angular_velocity(body);
        direction_change = direction_change.max(length([
            omega[0] - initial_omega[0],
            omega[1] - initial_omega[1],
            omega[2] - initial_omega[2],
        ]));
    }
    assert!(
        max_drift < 0.02 * length(initial),
        "angular momentum drifted by {max_drift} (|L| = {})",
        length(initial)
    );
    assert!(
        direction_change > 0.05,
        "an asymmetric body precesses: omega changed by only {direction_change}"
    );
}

#[test]
fn isotropic_body_has_no_gyroscopic_precession() {
    let mut engine = free_engine(DT, 100);
    let body = damping_free_body(&mut engine, 1.0, ColliderConfig::sphere(1.0));
    engine
        .set_angular_velocity(body, 1.0, 2.0, 3.0)
        .expect("spin");
    engine.step(1.0).expect("one second");
    assert_vec_close(
        engine.get_angular_velocity(body),
        [1.0, 2.0, 3.0],
        1e-9,
        "constant angular velocity",
    );
}

// ---------------------------------------------------------------------------
// Fixed-step time integration (F5659)
// ---------------------------------------------------------------------------

#[test]
fn step_of_one_second_integrates_only_the_capped_substeps() {
    // Defaults: fixed_dt = 1/60 and max_substeps = 4.
    let mut engine = free_engine(DT, 4);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(1.0).expect("step");

    assert_close(
        engine.time(),
        4.0 * DT,
        1e-12,
        "time() is the integrated time",
    );
    assert_close(
        engine.get_position(body)[0],
        4.0 * DT,
        1e-12,
        "distance travelled",
    );
    let info = engine.debug_info();
    assert_eq!(info.substeps_performed, 4);
    assert_close(info.dropped_time, 1.0 - 4.0 * DT, 1e-12, "dropped time");
    assert_eq!(
        engine.interpolation_alpha(),
        0.0,
        "the surplus is not carried over"
    );
}

#[test]
fn sub_fixed_dt_steps_accumulate_until_a_substep_fits() {
    let mut engine = free_engine(DT, 4);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);

    engine.step(0.01).expect("first partial step");
    assert_eq!(engine.time(), 0.0, "no substep has run yet");
    assert_eq!(engine.get_position(body), [0.0; 3]);
    assert_eq!(engine.debug_info().substeps_performed, 0);
    assert_close(
        engine.interpolation_alpha(),
        0.6,
        1e-12,
        "alpha after 0.01 s",
    );

    engine.step(0.01).expect("second partial step");
    assert_close(engine.time(), DT, 1e-12, "one substep has run");
    assert_close(engine.get_position(body)[0], DT, 1e-12, "position");
    assert_eq!(engine.debug_info().substeps_performed, 1);
    assert_close(
        engine.interpolation_alpha(),
        0.2,
        1e-9,
        "alpha after 0.02 s",
    );
}

#[test]
fn exact_multiples_of_the_fixed_step_are_fully_integrated() {
    let mut engine = free_engine(0.01, 10);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(0.05).expect("step");
    assert_eq!(engine.debug_info().substeps_performed, 5);
    assert_close(engine.time(), 0.05, 1e-12, "time");
    assert_close(engine.get_position(body)[0], 0.05, 1e-12, "position");
    assert_eq!(engine.debug_info().dropped_time, 0.0);
}

#[test]
fn a_whole_second_is_integrated_when_the_substep_cap_allows_it() {
    let mut engine = free_engine(DT, 100);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(1.0).expect("step");
    assert_eq!(engine.debug_info().substeps_performed, 60);
    assert_close(engine.time(), 1.0, 1e-9, "time equals the elapsed time");
    assert_close(engine.get_position(body)[0], 1.0, 1e-9, "position");
    assert_eq!(engine.debug_info().dropped_time, 0.0);
}

#[test]
fn remainders_carry_over_between_steps() {
    let mut engine = free_engine(0.1, 10);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(0.25).expect("first step");
    assert_close(engine.time(), 0.2, 1e-12, "two whole substeps");
    assert_close(
        engine.interpolation_alpha(),
        0.5,
        1e-9,
        "half a substep pending",
    );

    engine.step(0.05).expect("second step");
    assert_close(
        engine.time(),
        0.3,
        1e-12,
        "the pending half completes a third substep",
    );
    assert_close(engine.get_position(body)[0], 0.3, 1e-12, "position");
    assert!(engine.interpolation_alpha() < 1e-6);
}

#[test]
fn the_clock_does_not_drift_over_many_frames() {
    // 300 frames of 20 ms against a 1/60 s substep: 6 s of real time is
    // exactly 360 substeps.
    let mut engine = free_engine(DT, 4);
    let body = drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    let mut substeps = 0u32;
    for _ in 0..300 {
        engine.step(0.02).expect("frame");
        substeps += engine.debug_info().substeps_performed;
    }
    assert_eq!(substeps, 360);
    assert_close(engine.time(), 6.0, 1e-9, "simulation clock");
    assert_close(engine.get_position(body)[0], 6.0, 1e-9, "position");
    assert!(engine.interpolation_alpha() < 1e-6);
}

#[test]
fn max_substeps_of_one_runs_a_single_substep_per_call() {
    let mut engine = free_engine(DT, 1);
    drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(0.5).expect("step");
    assert_eq!(engine.debug_info().substeps_performed, 1);
    assert_close(engine.time(), DT, 1e-12, "time");
    assert_close(
        engine.debug_info().dropped_time,
        0.5 - DT,
        1e-12,
        "dropped time",
    );
}

#[test]
fn a_zero_time_step_is_a_no_op_that_keeps_pending_forces() {
    let mut engine = free_engine(DT, 4);
    let body = drifting_body(&mut engine, [0.0; 3]);
    engine.apply_force(body, 6.0, 0.0, 0.0).expect("force");
    engine.step(0.0).expect("zero step");
    assert_eq!(engine.time(), 0.0);
    assert_eq!(engine.get_velocity(body), [0.0; 3]);
    engine.step(DT).expect("real step");
    assert_close(
        engine.get_velocity(body)[0],
        6.0 * DT,
        1e-12,
        "the force was kept",
    );
}

#[test]
fn forces_act_during_every_substep_of_a_step() {
    let mut engine = free_engine(DT, 4);
    let body = engine
        .add_dynamic_body(2.0, 0.0, 0.0, 0.0)
        .expect("dynamic body");
    engine.set_body_linear_damping(body, 0.0).expect("damping");
    engine.apply_force(body, 4.0, 0.0, 0.0).expect("force");
    engine.step(3.0 * DT).expect("step");
    // a = F / m = 2 m/s^2 for three substeps.
    assert_close(
        engine.get_velocity(body)[0],
        2.0 * 3.0 * DT,
        1e-12,
        "velocity",
    );
}

#[test]
fn forces_are_cleared_after_every_step() {
    let mut engine = free_engine(DT, 4);
    let body = drifting_body(&mut engine, [0.0; 3]);
    engine.apply_force(body, 6.0, 0.0, 0.0).expect("force");
    engine.step(DT).expect("first step");
    let speed = engine.get_velocity(body)[0];
    assert_close(speed, 6.0 * DT, 1e-12, "first step");
    engine.step(DT).expect("second step");
    assert_close(
        engine.get_velocity(body)[0],
        speed,
        1e-15,
        "no force the second time",
    );
}

#[test]
fn set_fixed_dt_changes_the_substep_length() {
    let mut engine = WasmPhysicsEngine::new(0.0, -10.0, 0.0);
    let body = engine
        .add_dynamic_body(1.0, 0.0, 100.0, 0.0)
        .expect("dynamic body");
    engine.set_body_linear_damping(body, 0.0).expect("damping");
    engine.set_fixed_dt(0.5).expect("fixed dt");
    engine.step(0.5).expect("step");
    assert_eq!(engine.debug_info().substeps_performed, 1);
    assert_close(engine.time(), 0.5, 1e-15, "time");
    assert_close(engine.get_velocity(body)[1], -5.0, 1e-12, "v = g * dt");
}

#[test]
fn interpolation_alpha_stays_in_the_unit_interval() {
    let mut engine = free_engine(DT, 4);
    drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    for frame in 0..200 {
        engine.step(0.001 * (frame % 37) as f64).expect("step");
        let alpha = engine.interpolation_alpha();
        assert!((0.0..=1.0).contains(&alpha), "alpha {alpha} out of range");
    }
}

#[test]
fn reset_clears_the_pending_time() {
    let mut engine = free_engine(DT, 4);
    drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step(0.01).expect("partial step");
    assert!(engine.interpolation_alpha() > 0.0);
    engine.reset();
    assert_eq!(engine.interpolation_alpha(), 0.0);
    assert_eq!(engine.time(), 0.0);
    assert_eq!(engine.debug_info().substeps_performed, 0);
}

#[test]
fn contacts_are_kept_across_a_step_that_runs_no_substep() {
    let mut engine = free_engine(DT, 4);
    let a = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body a");
    let b = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0).expect("body b");
    engine.add_sphere_collider(a, 1.0).expect("collider a");
    engine.add_sphere_collider(b, 1.0).expect("collider b");
    engine.step(DT).expect("step with a substep");
    let contacts = engine.get_contact_count();
    assert!(contacts >= 1);
    engine.step(DT / 10.0).expect("step without a substep");
    assert_eq!(engine.debug_info().substeps_performed, 0);
    assert_eq!(
        engine.get_contact_count(),
        contacts,
        "the contacts of the previous step stay available"
    );
}

#[test]
fn js_step_and_alpha_wrappers_agree_with_the_rust_api() {
    let mut engine = free_engine(DT, 4);
    drifting_body(&mut engine, [1.0, 0.0, 0.0]);
    engine.step_js(0.01).expect("js step");
    assert_close(engine.interpolation_alpha_js(), 0.6, 1e-12, "alpha");
    assert_eq!(engine.substeps_performed_js(), 0);
    assert_eq!(engine.dropped_time_js(), 0.0);
    engine.step_js(1.0).expect("js step of a second");
    assert_eq!(engine.substeps_performed_js(), 4);
    assert_close(
        engine.dropped_time_js(),
        1.01 - 4.0 * DT,
        1e-9,
        "dropped time",
    );
    engine.set_fixed_dt_js(0.005).expect("js fixed dt");
    assert_eq!(engine.config().fixed_dt, 0.005);
    assert_eq!(engine.time_js(), engine.time());
}
