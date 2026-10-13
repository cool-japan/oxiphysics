// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Public-API tests for the wind and body-frame transport terms of the 6-DOF
//! aircraft model.
//!
//! * The body-frame translational equation is `dv/dt = F / m - omega x v`. With
//!   no force a velocity that is constant in the inertial frame must therefore
//!   *counter-rotate* in the body frame: it keeps its norm and turns by
//!   `-|omega| t` about the rotation axis (closed form: Rodrigues' formula).
//! * The aerodynamics see the velocity relative to the air mass, sampled at the
//!   state's own simulation time, which `integrate_aircraft` advances.

use oxiphysics_rigid::aerospace_rigid::{
    AircraftConfig, AircraftState, WindModel, body_frame_acceleration, compute_aero_forces,
    integrate_aircraft,
};
use std::f64::consts::PI;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn add_scaled(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
}

/// Rotate `v` by `angle` about the unit vector `axis` (Rodrigues' formula).
fn rotate_about(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (sin_a, cos_a) = angle.sin_cos();
    let k_cross_v = cross(axis, v);
    let k_dot_v = dot(axis, v);
    [
        v[0] * cos_a + k_cross_v[0] * sin_a + axis[0] * k_dot_v * (1.0 - cos_a),
        v[1] * cos_a + k_cross_v[1] * sin_a + axis[1] * k_dot_v * (1.0 - cos_a),
        v[2] * cos_a + k_cross_v[2] * sin_a + axis[2] * k_dot_v * (1.0 - cos_a),
    ]
}

/// Classical RK4 integration of the force-free body-frame velocity
/// `dv/dt = -omega x v` for `steps` steps of size `dt`.
fn integrate_free_velocity(
    v0: [f64; 3],
    omega: [f64; 3],
    mass: f64,
    dt: f64,
    steps: usize,
) -> [f64; 3] {
    let rate = |v: [f64; 3]| body_frame_acceleration([0.0; 3], mass, omega, v);
    let mut v = v0;
    for _ in 0..steps {
        let k1 = rate(v);
        let k2 = rate(add_scaled(v, k1, 0.5 * dt));
        let k3 = rate(add_scaled(v, k2, 0.5 * dt));
        let k4 = rate(add_scaled(v, k3, dt));
        for axis in 0..3 {
            v[axis] += dt / 6.0 * (k1[axis] + 2.0 * k2[axis] + 2.0 * k3[axis] + k4[axis]);
        }
    }
    v
}

#[test]
fn pure_pitch_rate_without_force_keeps_the_speed_and_turns_the_velocity() {
    let mass = 1500.0;
    let pitch_rate = 0.5;
    let v0 = [60.0, 0.0, 8.0];
    let (dt, steps) = (1.0e-3, 2000);
    let t = dt * steps as f64;
    let v = integrate_free_velocity(v0, [0.0, pitch_rate, 0.0], mass, dt, steps);
    assert!(
        (norm(v) - norm(v0)).abs() < 1e-9,
        "|v_b| changed from {} to {}",
        norm(v0),
        norm(v)
    );
    // A positive pitch rate raises the nose, so the inertially fixed velocity
    // ends up below the body x axis (w grows while u shrinks): seen from the
    // body it counter-rotates about +y by the angle -q t.
    let expected = rotate_about(v0, [0.0, 1.0, 0.0], -pitch_rate * t);
    for axis in 0..3 {
        assert!(
            (v[axis] - expected[axis]).abs() < 1e-9,
            "v[{axis}] = {} vs {}",
            v[axis],
            expected[axis]
        );
    }
    // The (u, w) components obey du/dt = -q w, dw/dt = q u at the start.
    let a0 = body_frame_acceleration([0.0; 3], mass, [0.0, pitch_rate, 0.0], v0);
    assert!((a0[0] + pitch_rate * v0[2]).abs() < 1e-12);
    assert!((a0[2] - pitch_rate * v0[0]).abs() < 1e-12);
    assert!(a0[1].abs() < 1e-12);
}

#[test]
fn general_body_rates_counter_rotate_the_force_free_velocity() {
    let mass = 900.0;
    for (omega, v0) in [
        ([0.1, 0.5, -0.2], [60.0, -3.0, 8.0]),
        ([-0.4, 0.05, 0.3], [35.0, 12.0, -6.0]),
        ([0.0, 0.0, 0.7], [80.0, 1.0, 0.0]),
    ] {
        let rate = norm(omega);
        let axis = [omega[0] / rate, omega[1] / rate, omega[2] / rate];
        let (dt, steps) = (1.0e-3, 2500);
        let t = dt * steps as f64;
        let v = integrate_free_velocity(v0, omega, mass, dt, steps);
        let expected = rotate_about(v0, axis, -rate * t);
        assert!(
            norm(add_scaled(v, expected, -1.0)) < 1e-9,
            "omega = {omega:?}: v = {v:?}, expected {expected:?}"
        );
        assert!((norm(v) - norm(v0)).abs() < 1e-9);
        // The component along the rotation axis is constant.
        assert!((dot(v, axis) - dot(v0, axis)).abs() < 1e-9);
    }
}

#[test]
fn force_adds_to_the_transport_term_independently_of_the_mass_scaling() {
    let omega = [0.2, -0.1, 0.4];
    let v = [50.0, 4.0, -2.0];
    let force = [3000.0, -600.0, 1200.0];
    let mass = 1200.0;
    let a = body_frame_acceleration(force, mass, omega, v);
    let w_cross_v = cross(omega, v);
    for axis in 0..3 {
        let expected = force[axis] / mass - w_cross_v[axis];
        assert!((a[axis] - expected).abs() < 1e-12);
    }
    // Doubling the mass halves the force term but leaves the transport term.
    let heavy = body_frame_acceleration(force, 2.0 * mass, omega, v);
    for axis in 0..3 {
        let expected = 0.5 * force[axis] / mass - w_cross_v[axis];
        assert!((heavy[axis] - expected).abs() < 1e-12);
    }
}

#[test]
fn gust_is_sampled_at_the_simulation_time_advanced_by_the_integrator() {
    let config = AircraftConfig::generic_ga();
    let mut wind = WindModel::calm();
    wind.gust_amplitude = 5.0;
    wind.gust_frequency = 2.0;

    let level = |time: f64| {
        let mut s = AircraftState::level_flight(1500.0, 60.0, 3000.0);
        s.time = time;
        s
    };
    // sin(2 t) = +1 at t = pi/4: a 5 m/s tailwind gust (the gust blows along
    // +x, the heading) lowers the airspeed and with it the drag.
    let drag_x = |time: f64| compute_aero_forces(&level(time), &config, &wind)[0];
    let (calm_drag, tail_drag, head_drag) = (drag_x(0.0), drag_x(PI / 4.0), drag_x(3.0 * PI / 4.0));
    assert!(
        tail_drag > calm_drag && calm_drag > head_drag,
        "Fx: tail {tail_drag}, calm {calm_drag}, head {head_drag}"
    );

    // The integrator advances the state's clock by exactly dt per step.
    let dt = 0.01;
    let mut s = level(PI / 4.0);
    for _ in 0..25 {
        integrate_aircraft(&mut s, &config, &wind, dt);
    }
    assert!(
        (s.time - (PI / 4.0 + 0.25)).abs() < 1e-12,
        "time = {}",
        s.time
    );

    // Same state, same step, different clock: the decelerations differ only
    // through the gust the integrator samples at the state's own time.
    let step_u = |time: f64| {
        let mut s = level(time);
        integrate_aircraft(&mut s, &config, &wind, dt);
        s.velocity_body[0]
    };
    assert!(step_u(PI / 4.0) > step_u(0.0));
    assert!(step_u(0.0) > step_u(3.0 * PI / 4.0));
}
