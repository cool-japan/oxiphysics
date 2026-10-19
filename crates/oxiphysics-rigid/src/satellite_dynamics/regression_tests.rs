// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the attitude PD controller and for the Hill-frame
//! formation model (bounded relative orbits, drift rate, separation accessors
//! and two-impulse transfers).

use super::*;
use crate::rng_util::uniform;
use rand::SeedableRng;
use rand::rngs::StdRng;

#[test]
fn test_pd_controller_drives_attitude_to_target() {
    // Unit inertia, kp = 1, kd = 2: a critically damped second-order loop.
    let mut attitude = AttitudeControl::new(1.0, 1.0, 1.0);
    let target = Quat::from_axis_angle(&Vec3::z_axis(), 0.3);
    for _ in 0..3000 {
        let torque = attitude.pd_controller(&target, 1.0, 2.0);
        attitude.apply_torque(torque);
        attitude.step(0.01);
    }
    let miss = attitude.quaternion.angle_to(&target);
    assert!(miss < 1e-3, "attitude error after 30 s: {miss} rad");
    assert!(attitude.angular_velocity.norm() < 1e-3);
}

#[test]
fn test_pd_controller_torque_is_expressed_in_body_axes() {
    // The body is rotated 0.4 rad about x; the target differs from it by a
    // further 0.2 rad rotation about the *body* y axis.  A world-frame
    // error axis would be rot(x, 0.4) * y = (0, cos 0.4, sin 0.4) instead.
    let mut attitude = AttitudeControl::new(100.0, 200.0, 300.0);
    attitude.quaternion = Quat::from_axis_angle(&Vec3::x_axis(), 0.4);
    let target = attitude.quaternion * Quat::from_axis_angle(&Vec3::y_axis(), 0.2);
    let torque = attitude.pd_controller(&target, 10.0, 5.0);
    assert!(
        (torque - Vec3::new(0.0, 2.0, 0.0)).norm() < 1e-12,
        "torque = {torque:?}"
    );
    // The derivative term opposes the body rate.
    attitude.angular_velocity = Vec3::new(0.1, 0.0, -0.2);
    let damped = attitude.pd_controller(&target, 10.0, 5.0);
    assert!((damped - Vec3::new(-0.5, 2.0, 1.0)).norm() < 1e-12);
}

#[test]
fn test_pd_controller_converges_from_arbitrary_attitude_with_gyroscopic_coupling() {
    let mut attitude = AttitudeControl::new(100.0, 200.0, 300.0);
    attitude.quaternion =
        Quat::from_axis_angle(&Vec3::x_axis(), 1.0) * Quat::from_axis_angle(&Vec3::y_axis(), -0.5);
    attitude.angular_velocity = Vec3::new(0.02, -0.01, 0.03);
    let target = Quat::from_scaled_axis(Vec3::new(1.0, 1.0, 0.0).normalize() * 2.5);
    for _ in 0..20_000 {
        let torque = attitude.pd_controller(&target, 50.0, 250.0);
        attitude.apply_torque(torque);
        attitude.step(0.01);
    }
    let miss = attitude.quaternion.angle_to(&target);
    assert!(miss < 1e-3, "attitude error after 200 s: {miss} rad");
    assert!(attitude.angular_velocity.norm() < 1e-3);
}

#[test]
fn test_pd_controller_takes_the_short_way_and_ignores_quaternion_sign() {
    let attitude = AttitudeControl::new(1.0, 1.0, 1.0);
    // 4 rad about +z is the same attitude as -(2 pi - 4) rad: turn the short way.
    let far = Quat::from_axis_angle(&Vec3::z_axis(), 4.0);
    let torque = attitude.pd_controller(&far, 1.0, 0.0);
    assert!(
        torque.z < 0.0,
        "must rotate in the negative sense: {torque:?}"
    );
    assert!(
        (torque.z + (2.0 * PI - 4.0)).abs() < 1e-12,
        "torque = {torque:?}"
    );
    // q and -q are one attitude.
    let near = Quat::from_axis_angle(&Vec3::y_axis(), 0.7);
    let flipped = Quat::new_unchecked(-near.into_inner());
    let a = attitude.pd_controller(&near, 3.0, 0.0);
    let b = attitude.pd_controller(&flipped, 3.0, 0.0);
    assert!((a - b).norm() < 1e-12, "{a:?} vs {b:?}");
    assert!((a - Vec3::new(0.0, 2.1, 0.0)).norm() < 1e-12);
}

#[test]
fn test_circular_orbit_stays_on_the_two_to_one_relative_ellipse() {
    let r_ref = 7e6_f64;
    let rho = 100.0;
    let phase = 0.3;
    let mut formation = SpacecraftFormation::new(
        r_ref,
        HillState::new(Vec3::zeros(), Vec3::zeros()),
        HillState::new(Vec3::zeros(), Vec3::zeros()),
    );
    formation.set_circular_orbit(rho, phase);
    assert!(formation.is_bounded_orbit());
    let n = formation.mean_motion();
    let period = 2.0 * PI / n;
    let samples_per_period = 40;
    let periods = 10;
    let mut mean_y = 0.0;
    let mut max_y = 0.0_f64;
    let mut max_x = 0.0_f64;
    for k in 0..(periods * samples_per_period) {
        let t = period * f64::from(k) / f64::from(samples_per_period);
        let state = formation.propagate_exact(t);
        let angle = n * t + phase;
        assert!(
            (state.position.x - rho * angle.sin()).abs() < 1e-6,
            "x(t) off the radial leg at sample {k}"
        );
        assert!(
            (state.position.y - 2.0 * rho * angle.cos()).abs() < 1e-6,
            "y(t) off the along-track leg at sample {k}"
        );
        assert!(state.position.z.abs() < 1e-9);
        mean_y += state.position.y;
        max_y = max_y.max(state.position.y.abs());
        max_x = max_x.max(state.position.x.abs());
    }
    mean_y /= f64::from(periods * samples_per_period);
    assert!(
        max_y <= 2.0 * rho * (1.0 + 1e-6),
        "along-track excursion {max_y}"
    );
    assert!(max_x <= rho * (1.0 + 1e-6), "radial excursion {max_x}");
    assert!(
        mean_y.abs() < 1e-6 * rho,
        "the orbit must be centred on the chief: {mean_y}"
    );
}

#[test]
fn test_circular_orbit_is_bounded_for_every_phase() {
    let r_ref = 6.9e6_f64;
    let mut formation = SpacecraftFormation::new(
        r_ref,
        HillState::new(Vec3::zeros(), Vec3::zeros()),
        HillState::new(Vec3::zeros(), Vec3::zeros()),
    );
    for k in 0..16 {
        let phase = f64::from(k) * PI / 8.0;
        formation.set_circular_orbit(250.0, phase);
        assert!(formation.is_bounded_orbit(), "phase {phase}");
        assert!(formation.secular_drift_rate().abs() < 1e-9, "phase {phase}");
    }
}

#[test]
fn test_stationary_along_track_offset_is_bounded() {
    let r_ref = 7e6_f64;
    let state = HillState::new(Vec3::new(0.0, 50.0, 0.0), Vec3::zeros());
    let formation = SpacecraftFormation::new(r_ref, state.clone(), state);
    assert!(formation.is_bounded_orbit());
    let later = formation.propagate_exact(3.7 * 2.0 * PI / formation.mean_motion());
    assert!((later.position - Vec3::new(0.0, 50.0, 0.0)).norm() < 1e-6);
}

#[test]
fn test_radial_offset_with_inertial_speed_mismatch_is_unbounded() {
    let r_ref = 7e6_f64;
    let n = (MU_EARTH / r_ref.powi(3)).sqrt();
    let rho = 100.0;
    let state = HillState::new(Vec3::new(rho, 0.0, 0.0), Vec3::new(0.0, 2.0 * n * rho, 0.0));
    let formation = SpacecraftFormation::new(r_ref, state.clone(), state);
    assert!(!formation.is_bounded_orbit());
    assert!((formation.secular_drift_rate() + 12.0 * n * rho).abs() < 1e-12);
}

#[test]
fn test_secular_drift_rate_is_the_slope_of_the_exact_solution() {
    let r_ref = 7.2e6_f64;
    let mut rng = StdRng::seed_from_u64(2718);
    for case in 0..8 {
        let position = Vec3::new(
            uniform(&mut rng, -300.0, 300.0),
            uniform(&mut rng, -300.0, 300.0),
            uniform(&mut rng, -300.0, 300.0),
        );
        let velocity = Vec3::new(
            uniform(&mut rng, -0.5, 0.5),
            uniform(&mut rng, -0.5, 0.5),
            uniform(&mut rng, -0.5, 0.5),
        );
        let state = HillState::new(position, velocity);
        let formation = SpacecraftFormation::new(r_ref, state.clone(), state);
        let period = 2.0 * PI / formation.mean_motion();
        let y1 = formation.propagate_exact(period).position.y;
        let y2 = formation.propagate_exact(2.0 * period).position.y;
        let slope = (y2 - y1) / period;
        let rate = formation.secular_drift_rate();
        assert!(
            (slope - rate).abs() <= 1e-9 * rate.abs().max(1e-6),
            "case {case}: exact slope {slope}, secular_drift_rate {rate}"
        );
    }
}

#[test]
fn test_separation_accessors_follow_the_hcw_axes() {
    let r_ref = 7e6_f64;
    let state = HillState::new(Vec3::new(3.0, 7.0, 1.0), Vec3::zeros());
    let formation = SpacecraftFormation::new(r_ref, state.clone(), state);
    assert_eq!(formation.radial_separation(), 3.0);
    assert_eq!(formation.along_track_separation(), 7.0);
    assert_eq!(formation.cross_track_separation(), 1.0);
}

#[test]
fn test_two_impulse_transfer_hits_the_target_state() {
    let r_ref = 7e6_f64;
    let mut rng = StdRng::seed_from_u64(4242);
    for case in 0..12 {
        let mut random_state = |pos_scale: f64, vel_scale: f64| {
            HillState::new(
                Vec3::new(
                    uniform(&mut rng, -pos_scale, pos_scale),
                    uniform(&mut rng, -pos_scale, pos_scale),
                    uniform(&mut rng, -pos_scale, pos_scale),
                ),
                Vec3::new(
                    uniform(&mut rng, -vel_scale, vel_scale),
                    uniform(&mut rng, -vel_scale, vel_scale),
                    uniform(&mut rng, -vel_scale, vel_scale),
                ),
            )
        };
        let deputy = random_state(1_500.0, 2.0);
        let target = random_state(1_500.0, 2.0);
        let formation = SpacecraftFormation::new(r_ref, deputy, target.clone());
        let period = 2.0 * PI / formation.mean_motion();
        let tf = 0.37 * period;
        let (dv1, dv2) = formation
            .two_impulse_transfer(tf)
            .expect("0.37 periods is not a singular transfer time");

        let mut burned = formation.clone();
        burned.apply_impulse(dv1);
        assert!((burned.delta_v_total - dv1.norm()).abs() < 1e-12);
        let arrival = burned.propagate_exact(tf);
        let miss = (arrival.position - target.position).norm();
        assert!(
            miss < 1e-6,
            "case {case}: arrival misses the target by {miss} m"
        );
        let velocity_error = (arrival.velocity + dv2 - target.velocity).norm();
        assert!(
            velocity_error < 1e-8,
            "case {case}: velocity after the second burn is off by {velocity_error} m/s"
        );
    }
}

#[test]
fn test_two_impulse_transfer_points_toward_the_target() {
    // Close up a 100 m along-track gap: the first burn must have a positive
    // along-track component (towards the target), whereas the old
    // `error * n / sin(nt)` guess, with error = deputy - target, pointed away.
    let r_ref = 7e6_f64;
    let deputy = HillState::new(Vec3::zeros(), Vec3::zeros());
    let target = HillState::new(Vec3::new(0.0, 100.0, 0.0), Vec3::zeros());
    let formation = SpacecraftFormation::new(r_ref, deputy, target.clone());
    let tf = 0.25 * 2.0 * PI / formation.mean_motion();
    let (dv1, dv2) = formation
        .two_impulse_transfer(tf)
        .expect("regular transfer");
    assert!(dv1.y > 0.0, "dv1 = {dv1:?}");
    let mut burned = formation.clone();
    burned.apply_impulse(dv1);
    let arrival = burned.propagate_exact(tf);
    assert!((arrival.position - target.position).norm() < 1e-6);
    assert!((arrival.velocity + dv2).norm() < 1e-8, "target is at rest");
}

#[test]
fn test_two_impulse_transfer_rejects_singular_times() {
    let r_ref = 7e6_f64;
    let deputy = HillState::new(Vec3::new(10.0, 20.0, 5.0), Vec3::new(0.0, 0.1, 0.0));
    let target = HillState::new(Vec3::zeros(), Vec3::zeros());
    let formation = SpacecraftFormation::new(r_ref, deputy, target);
    let period = 2.0 * PI / formation.mean_motion();
    assert!(formation.two_impulse_transfer(0.0).is_none());
    assert!(formation.two_impulse_transfer(-5.0).is_none());
    assert!(formation.two_impulse_transfer(period).is_none());
    assert!(formation.two_impulse_transfer(0.5 * period).is_none());
    assert!(formation.two_impulse_transfer(0.4 * period).is_some());
}
