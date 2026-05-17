// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the vehicle API.
//!
//! Included from `vehicle_api.rs` via `#[path = ...]` to keep the main
//! API file under the 2000-line refactor budget.

use super::vehicle_api::{
    AbsState, EngineCurve, FialaParams, PyDrivetrain, PyRacingLine, PyStabilityControl, PySteering,
    PySuspension, PyTelemetry, PyTireModel, PyTireThermal, PyVehicle, PyWheelDynamics, TcsState,
    TelemetrySample, TrackPoint, ackermann_angles, compute_slip_angle, compute_slip_ratio,
    pacejka_fx, pacejka_fy,
};

#[test]
fn test_pacejka_fy_zero_slip() {
    let fy = pacejka_fy(0.0, 10.0, 1.9, 1.0, 0.97);
    assert!(fy.abs() < 1e-9);
}

#[test]
fn test_pacejka_fy_nonzero() {
    let fy = pacejka_fy(0.1, 10.0, 1.9, 1.0, 0.97);
    assert!(fy.abs() > 0.0);
}

#[test]
fn test_pacejka_fx_peak() {
    let fx = pacejka_fx(0.2, 11.0, 1.65, 1.0, 0.0);
    assert!(fx.abs() > 0.0);
}

#[test]
fn test_compute_slip_angle() {
    let alpha = compute_slip_angle(1.0, 10.0);
    assert!(alpha.abs() < 0.2);
}

#[test]
fn test_compute_slip_ratio_no_slip() {
    let kappa = compute_slip_ratio(10.0, 10.0);
    assert!(kappa.abs() < 1e-9);
}

#[test]
fn test_compute_slip_ratio_drive_slip() {
    let kappa = compute_slip_ratio(12.0, 10.0);
    assert!(kappa > 0.0);
}

#[test]
fn test_ackermann_straight() {
    let (a, b) = ackermann_angles(0.0, 2.5, 1.5);
    assert!(a.abs() < 1e-9);
    assert!(b.abs() < 1e-9);
}

#[test]
fn test_ackermann_turn() {
    let (inner, outer) = ackermann_angles(0.3, 2.5, 1.5);
    assert!(inner.abs() > outer.abs());
}

#[test]
fn test_tire_model_forces() {
    let tm = PyTireModel::road_pacejka(0.32);
    let (fx, fy) = tm.compute_forces(0.1, 0.05, 3000.0);
    let _ = fx;
    assert!(fy.abs() > 0.0);
}

#[test]
fn test_tire_thermal_update() {
    let mut tt = PyTireThermal::new();
    let grip0 = tt.grip_scale();
    tt.update(5000.0, 10.0);
    assert!(tt.surface_temp > 25.0);
    // Check grip scale is valid
    assert!(tt.grip_scale() >= 0.0 && tt.grip_scale() <= 1.0);
    let _ = grip0;
}

#[test]
fn test_tire_thermal_optimal() {
    let mut tt = PyTireThermal::new();
    tt.surface_temp = 90.0;
    assert!((tt.grip_scale() - 1.0).abs() < 1e-9);
}

#[test]
fn test_engine_curve_flat() {
    let ec = EngineCurve::flat(300.0, 7000.0);
    let t = ec.torque_at(3500.0);
    assert!(t > 0.0);
    assert!(t <= 300.0 * 1.01);
}

#[test]
fn test_drivetrain_upshift() {
    let mut dt = PyDrivetrain::six_speed(300.0);
    dt.upshift();
    assert_eq!(dt.current_gear, 2);
}

#[test]
fn test_drivetrain_downshift() {
    let mut dt = PyDrivetrain::six_speed(300.0);
    dt.downshift();
    assert_eq!(dt.current_gear, 0);
}

#[test]
fn test_drivetrain_wheel_torque() {
    let dt = PyDrivetrain::six_speed(300.0);
    let t = dt.wheel_torque(1.0);
    assert!(t > 0.0);
}

#[test]
fn test_steering_angles() {
    let mut steer = PySteering::new(2.6, 1.5);
    steer.input = 1.0;
    let fa = steer.front_angle();
    assert!(fa > 0.0);
    assert!(fa <= steer.max_wheel_angle);
}

#[test]
fn test_stability_tcs() {
    let mut sc = PyStabilityControl::new();
    let (ts, _) = sc.update(vec![0.5, 0.5, 0.5, 0.5], 0.0, 0.0);
    assert!(ts < 1.0);
}

#[test]
fn test_stability_disabled() {
    let mut sc = PyStabilityControl::new();
    sc.enabled = false;
    let (ts, bs) = sc.update(vec![2.0, 2.0, 2.0, 2.0], 0.0, 0.0);
    assert!((ts - 1.0).abs() < 1e-9);
    assert!((bs - 1.0).abs() < 1e-9);
}

#[test]
fn test_tcs_state_update() {
    let mut tcs = TcsState::new();
    let scale = tcs.update(vec![0.0, 0.0, 0.0, 0.0]);
    assert!((scale - 1.0).abs() < 1e-9);
    assert!(!tcs.active);
}

#[test]
fn test_abs_state_update() {
    let mut abs = AbsState::new();
    let scale = abs.update(vec![0.0, 0.0, 0.0, 0.0]);
    assert!((scale - 1.0).abs() < 1e-9);
    assert!(!abs.active);
}

#[test]
fn test_suspension_force() {
    let mut s = PySuspension::double_wishbone();
    s.step(0.03, 0.01);
    let f = s.force();
    assert!(f < 0.0); // spring pushes back
}

#[test]
fn test_suspension_bump_stop() {
    let mut s = PySuspension::double_wishbone();
    s.step(s.max_compression, 0.01);
    let f = s.force();
    assert!(f < 0.0);
}

#[test]
fn test_wheel_dynamics_step() {
    let mut w = PyWheelDynamics::new(1.5, 0.32);
    w.fx = 100.0;
    w.grounded = true;
    w.step_omega(200.0, 0.0, 0.01);
    assert!(w.omega > 0.0);
}

#[test]
fn test_vehicle_sedan_step() {
    let mut v = PyVehicle::sedan();
    v.step(0.016, 0.5, 0.0, 0.0);
    // Vehicle should have moved or velocity changed
    let speed = v.forward_speed();
    let _ = speed;
    assert!(v.drivetrain.engine_rpm > 0.0);
}

#[test]
fn test_vehicle_racing_step() {
    let mut v = PyVehicle::racing_car();
    v.step(0.016, 1.0, 0.0, 0.1);
    assert!(v.gear_display >= 1);
}

#[test]
fn test_telemetry_record() {
    let mut tel = PyTelemetry::new(100);
    let mut sample = TelemetrySample::new(0.0);
    sample.speed = 50.0;
    sample.g_lon = 0.5;
    sample.g_lat = 0.3;
    sample.gear = 3;
    sample.rpm = 4000.0;
    sample.throttle = 0.8;
    sample.steer = 0.1;
    sample.tire_temps = [80.0; 4];
    tel.record(sample);
    assert_eq!(tel.sample_count(), 1);
}

#[test]
fn test_telemetry_lap_stats() {
    let mut tel = PyTelemetry::new(100);
    for i in 0..10 {
        let mut sample = TelemetrySample::new(i as f64 * 0.1);
        sample.speed = 50.0 + i as f64;
        sample.g_lon = 0.1;
        sample.g_lat = 0.2;
        sample.gear = 3;
        sample.rpm = 4000.0;
        sample.throttle = 0.8;
        sample.tire_temps = [80.0; 4];
        tel.record(sample);
    }
    let stats = tel
        .lap_stats()
        .expect("stats should be available with samples");
    assert!(stats.max_speed > 50.0);
    assert!(stats.sample_count == 10);
}

#[test]
fn test_racing_line_curvature() {
    let pts: Vec<TrackPoint> = (0..10)
        .map(|i| TrackPoint::new(i as f64, 0.0, 3.0))
        .collect();
    let rl = PyRacingLine::new(pts);
    // Straight line should have near-zero curvature
    let k = rl.curvature_at(5);
    assert!(k < 1e-9);
}

#[test]
fn test_racing_line_optimize() {
    let pts: Vec<TrackPoint> = (0..20)
        .map(|i| {
            let t = i as f64 / 20.0 * 2.0 * std::f64::consts::PI;
            TrackPoint::new(t.cos() * 10.0, t.sin() * 10.0, 2.0)
        })
        .collect();
    let mut rl = PyRacingLine::new(pts);
    let k0 = rl.total_curvature();
    rl.optimize(10, 0.1);
    let k1 = rl.total_curvature();
    // After optimization curvature should not increase dramatically
    assert!(k1 <= k0 * 1.5 + 0.1);
}

#[test]
fn test_fiala_tire_linear() {
    let fp = FialaParams::new(50000.0, 1.0);
    let fy = fp.lateral_force(0.05, 4000.0);
    let expected = -50000.0 * 0.05;
    assert!((fy - expected).abs() < 1.0);
}

#[test]
fn test_fiala_tire_saturated() {
    let fp = FialaParams::new(50000.0, 1.0);
    let fy = fp.lateral_force(1.0, 4000.0);
    assert!(fy.abs() <= 4000.0 + 1.0);
}
