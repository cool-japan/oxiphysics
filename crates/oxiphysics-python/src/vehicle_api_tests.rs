// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the vehicle API.
//!
//! Included from `vehicle_api.rs` via `#[path = ...]` to keep the main
//! API file under the 2000-line refactor budget.

use super::vehicle_api::{
    AbsState, EngineCurve, FialaParams, PacejkaCoeffs, PyDrivetrain, PyRacingLine,
    PyStabilityControl, PySteering, PySuspension, PyTelemetry, PyTireModel, PyTireThermal,
    PyVehicle, PyWheelDynamics, TcsState, TelemetrySample, TireModelKind, TrackPoint,
    ackermann_angles, compute_slip_angle, compute_slip_ratio, pacejka_fx, pacejka_fy,
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
    // Half throttle accelerates the vehicle along its heading.
    assert!(v.forward_speed() > 0.01, "speed {}", v.forward_speed());
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
    // The force has the sign of the slip angle, like the Pacejka model.
    let fp = FialaParams::new(50000.0, 1.0);
    let fy = fp.lateral_force(0.05, 4000.0);
    let expected = 50000.0 * 0.05;
    assert!((fy - expected).abs() < 1.0);
    assert!((fp.lateral_force(-0.05, 4000.0) + expected).abs() < 1.0);
}

#[test]
fn test_fiala_tire_saturated() {
    let fp = FialaParams::new(50000.0, 1.0);
    let fy = fp.lateral_force(1.0, 4000.0);
    assert!(fy.abs() <= 4000.0 + 1.0);
    assert!(fy > 0.0);
    assert!(fp.lateral_force(-1.0, 4000.0) < 0.0);
}

// ---------------------------------------------------------------------------
// Pacejka units and the tyre sign convention
// ---------------------------------------------------------------------------

/// Largest value of `f` on `[0, limit]` and where it is reached (fine grid).
fn peak_of(f: impl Fn(f64) -> f64, limit: f64) -> (f64, f64) {
    let steps = 20_000;
    let mut best = (0.0, f(0.0));
    for k in 1..=steps {
        let x = limit * k as f64 / steps as f64;
        let value = f(x);
        if value > best.1 {
            best = (x, value);
        }
    }
    best
}

/// The magic formula written out as in the textbook, for slip `x`.
fn textbook_magic_formula(x: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    d * (c * (b * x - e * (b * x - (b * x).atan())).atan()).sin()
}

#[test]
fn test_pacejka_lateral_force_peaks_at_a_realistic_slip_angle() {
    for (coeffs, expected_range) in [
        (PacejkaCoeffs::road_lateral(), 0.08..0.2),
        (PacejkaCoeffs::slick_lateral(), 0.08..0.2),
    ] {
        let (alpha, peak) = peak_of(
            |alpha| pacejka_fy(alpha, coeffs.b, coeffs.c, coeffs.d, coeffs.e),
            0.6,
        );
        assert!(
            expected_range.contains(&alpha),
            "peak at {alpha} rad for B={}",
            coeffs.b
        );
        assert!((peak - coeffs.d).abs() < 1e-6, "peak value {peak}");
        // Beyond the peak the force falls off instead of saturating early.
        let past = pacejka_fy(0.6, coeffs.b, coeffs.c, coeffs.d, coeffs.e);
        assert!(past < peak - 1e-3, "no fall-off: {past} vs {peak}");
    }
}

#[test]
fn test_pacejka_longitudinal_force_peaks_at_a_realistic_slip_ratio() {
    let c = PacejkaCoeffs::road_longitudinal();
    let (kappa, peak) = peak_of(|kappa| pacejka_fx(kappa, c.b, c.c, c.d, c.e), 1.0);
    assert!(kappa > 0.08 && kappa < 0.2, "peak at slip ratio {kappa}");
    assert!((peak - c.d).abs() < 1e-6, "peak value {peak}");
    // A 5 % slip ratio already carries most of the peak force.
    assert!(pacejka_fx(0.05, c.b, c.c, c.d, c.e) > 0.6 * peak);
}

#[test]
fn test_pacejka_is_the_textbook_formula_in_radians_and_ratio() {
    let lateral = PacejkaCoeffs::road_lateral();
    let longitudinal = PacejkaCoeffs::road_longitudinal();
    for x in [-0.7, -0.25, -0.1, -0.02, 0.0, 0.02, 0.05, 0.1, 0.25, 0.7] {
        let fy = pacejka_fy(x, lateral.b, lateral.c, lateral.d, lateral.e);
        let expected = textbook_magic_formula(x, lateral.b, lateral.c, lateral.d, lateral.e);
        assert!(
            (fy - expected).abs() < 1e-12,
            "Fy({x}) = {fy} vs {expected}"
        );
        let fx = pacejka_fx(
            x,
            longitudinal.b,
            longitudinal.c,
            longitudinal.d,
            longitudinal.e,
        );
        let expected = textbook_magic_formula(
            x,
            longitudinal.b,
            longitudinal.c,
            longitudinal.d,
            longitudinal.e,
        );
        assert!(
            (fx - expected).abs() < 1e-12,
            "Fx({x}) = {fx} vs {expected}"
        );
        // Odd in the slip.
        let mirrored = pacejka_fy(-x, lateral.b, lateral.c, lateral.d, lateral.e);
        assert!((fy + mirrored).abs() < 1e-12);
    }
}

#[test]
fn test_pacejka_with_zero_stiffness_gives_no_force() {
    assert_eq!(pacejka_fy(0.1, 0.0, 1.9, 1.0, 0.5), 0.0);
    assert_eq!(pacejka_fx(0.1, 0.0, 1.65, 1.0, 0.0), 0.0);
}

#[test]
fn test_slip_angle_sign_convention_and_reversing() {
    assert!(compute_slip_angle(-1.0, 10.0) > 0.0);
    assert!(compute_slip_angle(1.0, 10.0) < 0.0);
    assert!((compute_slip_angle(1.0, 10.0) + 0.1_f64.atan()).abs() < 1e-15);
    // Reversing does not flip the sign: the force still opposes the sliding.
    assert_eq!(
        compute_slip_angle(1.0, -10.0),
        compute_slip_angle(1.0, 10.0)
    );
    // Below 1e-4 m/s the angle is undefined.
    assert_eq!(compute_slip_angle(1.0, 5.0e-5), 0.0);
}

/// Linear tyre model of the given stiffnesses.
fn linear_tire(cornering_stiffness: f64, longitudinal_stiffness: f64) -> PyTireModel {
    PyTireModel {
        kind: TireModelKind::Linear {
            cornering_stiffness,
            longitudinal_stiffness,
        },
        radius: 0.32,
        width: 0.225,
        unloaded_radius: 0.32,
    }
}

/// Fiala tyre model.
fn fiala_tire() -> PyTireModel {
    PyTireModel {
        kind: TireModelKind::Fiala(FialaParams::new(60_000.0, 1.0)),
        radius: 0.32,
        width: 0.225,
        unloaded_radius: 0.32,
    }
}

#[test]
fn test_every_tire_model_gives_forces_with_the_sign_of_the_slips() {
    let tires = [
        PyTireModel::road_pacejka(0.32),
        PyTireModel::slick_pacejka(0.3),
        linear_tire(60_000.0, 8.0),
        fiala_tire(),
    ];
    for (index, tire) in tires.iter().enumerate() {
        let (_, fy_positive) = tire.compute_forces(0.04, 0.0, 4000.0);
        let (_, fy_negative) = tire.compute_forces(-0.04, 0.0, 4000.0);
        assert!(
            fy_positive > 0.0 && fy_negative < 0.0,
            "tire {index}: Fy = {fy_positive} / {fy_negative}"
        );
        assert!((fy_positive + fy_negative).abs() < 1e-9, "tire {index}");
    }
    // The Fiala model carries no longitudinal force, the others do.
    for (index, tire) in tires[..3].iter().enumerate() {
        let (fx_positive, _) = tire.compute_forces(0.0, 0.05, 4000.0);
        let (fx_negative, _) = tire.compute_forces(0.0, -0.05, 4000.0);
        assert!(
            fx_positive > 0.0 && fx_negative < 0.0,
            "tire {index}: Fx = {fx_positive} / {fx_negative}"
        );
    }
    let linear = tires[2].compute_forces(0.04, 0.05, 4000.0);
    assert!((linear.0 - 8.0 * 0.05 * 4000.0).abs() < 1e-9);
    assert!((linear.1 - 60_000.0 * 0.04).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// Engine, drivetrain and wheel brake
// ---------------------------------------------------------------------------

#[test]
fn test_engine_curve_interpolation_and_mismatched_lists() {
    let flat = EngineCurve::flat(300.0, 7000.0);
    assert!((flat.torque_at(-50.0) - 180.0).abs() < 1e-12);
    assert!((flat.torque_at(1050.0) - 240.0).abs() < 1e-12);
    assert!((flat.torque_at(3000.0) - 300.0).abs() < 1e-12);
    assert!((flat.torque_at(7000.0) - 240.0).abs() < 1e-12);
    assert!((flat.torque_at(9000.0) - 240.0).abs() < 1e-12);
    // Lists of different length are evaluated over their common points.
    let ragged = EngineCurve::new(vec![0.0, 1000.0, 2000.0], vec![100.0, 200.0], 7000.0, 800.0);
    assert!((ragged.torque_at(500.0) - 150.0).abs() < 1e-12);
    assert!((ragged.torque_at(5000.0) - 200.0).abs() < 1e-12);
    let no_torque = EngineCurve::new(vec![0.0, 1000.0], Vec::new(), 7000.0, 800.0);
    assert_eq!(no_torque.torque_at(500.0), 0.0);
    let empty = EngineCurve::new(Vec::new(), Vec::new(), 7000.0, 800.0);
    assert_eq!(empty.torque_at(500.0), 0.0);
}

#[test]
fn test_rev_limiter_fades_the_torque_out_at_the_redline() {
    let mut drivetrain = PyDrivetrain::six_speed(300.0);
    let redline = drivetrain.engine_curve.redline;
    let to_wheel = drivetrain.current_ratio() * drivetrain.final_drive * drivetrain.clutch;
    // Full torque up to 95 % of the redline.
    drivetrain.engine_rpm = 0.95 * redline;
    let full = drivetrain.engine_curve.torque_at(drivetrain.engine_rpm) * to_wheel;
    assert!((drivetrain.wheel_torque(1.0) - full).abs() < 1e-9 * full);
    // Halfway through the fade the torque is halved, at the redline it is gone.
    drivetrain.engine_rpm = 0.975 * redline;
    let halfway = 0.5 * drivetrain.engine_curve.torque_at(drivetrain.engine_rpm) * to_wheel;
    assert!((drivetrain.wheel_torque(1.0) - halfway).abs() < 1e-9 * halfway);
    drivetrain.engine_rpm = redline;
    assert_eq!(drivetrain.wheel_torque(1.0), 0.0);
    drivetrain.engine_rpm = 1.2 * redline;
    assert_eq!(drivetrain.wheel_torque(1.0), 0.0);
    drivetrain.engine_rpm = 3000.0;
    assert!(drivetrain.wheel_torque(1.0) > 0.0);
    assert_eq!(drivetrain.wheel_torque(0.0), 0.0);
}

/// A wheel on the ground (no air damping) turning at `omega`.
fn grounded_wheel(omega: f64) -> PyWheelDynamics {
    let mut wheel = PyWheelDynamics::new(1.5, 0.32);
    wheel.omega = omega;
    wheel.grounded = true;
    wheel
}

#[test]
fn test_wheel_brake_torque_opposes_the_rotation_and_never_reverses_it() {
    let mut forward = grounded_wheel(5.0);
    forward.step_omega(0.0, 100.0, 0.01);
    assert!((forward.omega - (5.0 - 100.0 / 1.5 * 0.01)).abs() < 1e-12);
    // A wheel turning backwards is slowed down, too.
    let mut backward = grounded_wheel(-5.0);
    backward.step_omega(0.0, 100.0, 0.01);
    assert!((backward.omega - (-5.0 + 100.0 / 1.5 * 0.01)).abs() < 1e-12);
    // Braking for a long time stops the wheel at zero.
    for _ in 0..1000 {
        forward.step_omega(0.0, 100.0, 0.01);
        backward.step_omega(0.0, 100.0, 0.01);
    }
    assert_eq!(forward.omega, 0.0);
    assert_eq!(backward.omega, 0.0);
    // The brake holds a stopped wheel against a weaker drive torque ...
    let mut held = grounded_wheel(0.0);
    held.step_omega(50.0, 100.0, 0.01);
    assert_eq!(held.omega, 0.0);
    // ... but a stronger drive torque still turns it, by the surplus only.
    held.step_omega(200.0, 100.0, 0.01);
    assert!((held.omega - 100.0 / 1.5 * 0.01).abs() < 1e-12);
    // Without brake the wheel follows the net torque as before.
    let mut free = grounded_wheel(0.0);
    free.step_omega(150.0, 0.0, 0.01);
    assert!((free.omega - 1.0).abs() < 1e-12);
    // A negative brake torque is not a drive torque.
    let mut spinning = grounded_wheel(5.0);
    spinning.step_omega(0.0, -100.0, 0.01);
    assert_eq!(spinning.omega, 5.0);
}

// ---------------------------------------------------------------------------
// Vehicle dynamics
// ---------------------------------------------------------------------------

/// Kinetic energy of the chassis (translation and yaw) and of the wheel spin.
fn vehicle_energy(vehicle: &PyVehicle) -> f64 {
    let translation =
        0.5 * vehicle.mass * (vehicle.velocity[0].powi(2) + vehicle.velocity[2].powi(2));
    let yaw = 0.5 * vehicle.inertia[2] * vehicle.angular_velocity[2].powi(2);
    let spin: f64 = vehicle
        .wheels
        .iter()
        .map(|wheel| 0.5 * wheel.inertia * wheel.omega * wheel.omega)
        .sum();
    translation + yaw + spin
}

/// A sedan rolling at `speed` m/s along its heading on free-rolling wheels.
fn rolling_sedan(speed: f64) -> PyVehicle {
    let mut vehicle = PyVehicle::sedan();
    vehicle.velocity = [speed, 0.0, 0.0];
    for wheel in &mut vehicle.wheels {
        wheel.omega = speed / wheel.radius;
    }
    vehicle
}

/// Step `vehicle` for `seconds` with constant inputs.
fn drive(vehicle: &mut PyVehicle, seconds: f64, dt: f64, throttle: f64, brake: f64, steer: f64) {
    for _ in 0..(seconds / dt).round() as usize {
        vehicle.step(dt, throttle, brake, steer);
    }
}

/// Yaw rate of a slip-free turn at the current speed and steering angle.
fn kinematic_yaw_rate(vehicle: &PyVehicle) -> f64 {
    vehicle.forward_speed() * vehicle.steering.front_angle().tan() / vehicle.steering.wheelbase
}

fn assert_vehicle_finite(vehicle: &PyVehicle) {
    let mut values: Vec<f64> = Vec::new();
    values.extend(vehicle.position);
    values.extend(vehicle.velocity);
    values.extend(vehicle.orientation);
    values.extend(vehicle.angular_velocity);
    values.push(vehicle.drivetrain.engine_rpm);
    for wheel in &vehicle.wheels {
        values.extend([
            wheel.omega,
            wheel.slip_ratio,
            wheel.slip_angle,
            wheel.fx,
            wheel.fy,
            wheel.fz,
        ]);
    }
    assert!(values.iter().all(|v| v.is_finite()), "{vehicle:?}");
}

#[test]
fn test_vehicle_throttle_accelerates_the_sedan_through_tyre_traction() {
    let mut v = PyVehicle::sedan();
    drive(&mut v, 2.0, 1.0 / 60.0, 1.0, 0.0, 0.0);
    let speed = v.forward_speed();
    assert!(speed > 5.0 && speed < 12.0, "speed after 2 s: {speed}");
    assert!(v.position[0] > 5.0, "position {}", v.position[0]);
    // The driven rear wheels carry the slip that produces the traction ...
    for i in [2, 3] {
        let wheel = &v.wheels[i];
        assert!(
            wheel.slip_ratio > 0.02 && wheel.slip_ratio < 0.2,
            "rear wheel {i} slip ratio {}",
            wheel.slip_ratio
        );
        assert!(wheel.fx > 1000.0, "rear wheel {i} Fx {}", wheel.fx);
        assert!(
            wheel.fx <= wheel.fz + 1e-9,
            "traction above the friction limit"
        );
        assert!(wheel.grounded);
        assert!((wheel.fz - v.mass * 9.81 / 4.0).abs() < 1e-9);
        // The stored slip is the slip of the ground speed under the wheel.
        let slip = compute_slip_ratio(wheel.peripheral_speed(), v.forward_speed());
        assert!(
            (wheel.slip_ratio - slip).abs() < 2e-3,
            "{} vs {slip}",
            wheel.slip_ratio
        );
    }
    // ... the free-rolling front wheels turn almost without slip, but not at
    // exactly zero.
    for i in [0, 1] {
        let slip = v.wheels[i].slip_ratio;
        assert!(
            slip != 0.0 && slip.abs() < 0.02,
            "front wheel {i} slip {slip}"
        );
    }
    // Straight-line driving of a symmetric vehicle stays straight.
    assert!(v.lateral_speed().abs() < 1e-9);
    assert!(v.orientation[2].abs() < 1e-9);
    assert!(v.angular_velocity[2].abs() < 1e-9);
    assert_vehicle_finite(&v);
}

#[test]
fn test_vehicle_launch_does_not_depend_on_the_time_step() {
    let speeds: Vec<f64> = [1.0 / 240.0, 1.0 / 60.0, 1.0 / 30.0]
        .iter()
        .map(|&dt| {
            let mut v = PyVehicle::sedan();
            drive(&mut v, 2.0, dt, 1.0, 0.0, 0.0);
            v.forward_speed()
        })
        .collect();
    for speed in &speeds {
        assert!((speed - speeds[0]).abs() < 0.02 * speeds[0], "{speeds:?}");
    }
}

#[test]
fn test_vehicle_without_throttle_does_not_accelerate() {
    let mut v = PyVehicle::sedan();
    drive(&mut v, 3.0, 1.0 / 60.0, 0.0, 0.0, 0.0);
    assert_eq!(v.forward_speed(), 0.0);
    assert_eq!(v.position, [0.0; 3]);
    assert!(v.wheels.iter().all(|wheel| wheel.omega == 0.0));
}

#[test]
fn test_vehicle_brakes_to_a_stop_without_rolling_backwards() {
    let mut v = rolling_sedan(20.0);
    let dt = 1.0 / 60.0;
    let mut slowest = f64::INFINITY;
    let mut speed_after_one_second = 0.0;
    for step in 0..360 {
        v.step(dt, 0.0, 1.0, 0.0);
        slowest = slowest.min(v.forward_speed());
        if step == 59 {
            speed_after_one_second = v.forward_speed();
        }
        assert!(
            v.wheels.iter().all(|wheel| wheel.omega > -1e-6),
            "a wheel reversed at step {step}"
        );
    }
    assert!(slowest > -1e-6, "the vehicle rolled backwards: {slowest}");
    assert!(v.forward_speed() < 0.01, "speed {}", v.forward_speed());
    assert!(v.velocity[0].abs() < 0.01);
    // Braking is effective but limited by friction (mu = 1, so at most 1 g).
    assert!(
        speed_after_one_second < 15.0 && speed_after_one_second > 20.0 - 9.81 - 0.5,
        "speed after one second of braking: {speed_after_one_second}"
    );
}

#[test]
fn test_vehicle_at_rest_stays_put_under_the_brake_and_steering() {
    let mut v = PyVehicle::sedan();
    drive(&mut v, 5.0, 1.0 / 60.0, 0.0, 1.0, 0.7);
    assert_eq!(v.position, [0.0; 3]);
    assert_eq!(v.orientation[2], 0.0);
    assert!(v.velocity.iter().all(|c| c.abs() < 1e-12));
}

#[test]
fn test_vehicle_steering_turns_it_at_the_kinematic_yaw_rate() {
    let mut v = rolling_sedan(10.0);
    drive(&mut v, 4.0, 1.0 / 60.0, 0.0, 0.0, 0.1);
    let kinematic = kinematic_yaw_rate(&v);
    let yaw_rate = v.angular_velocity[2];
    assert!(yaw_rate > 0.0 && kinematic > 0.0);
    assert!(
        (yaw_rate - kinematic).abs() < 0.05 * kinematic,
        "yaw rate {yaw_rate} vs kinematic {kinematic}"
    );
    // The heading and the path turn toward the positive lateral side.
    assert!(
        v.orientation[2] > 0.15 && v.orientation[2] < 0.35,
        "yaw {}",
        v.orientation[2]
    );
    assert!(v.position[2] > 0.5, "lateral position {}", v.position[2]);
    assert!(v.lateral_speed() > 0.0 && v.lateral_speed() < 0.2);
    // The front tyres run at the slip angle that carries the cornering force.
    assert!(v.wheels[0].slip_angle > 0.0 && v.wheels[0].slip_angle < 0.05);
    assert!(v.wheels[0].fy > 0.0);
}

#[test]
fn test_vehicle_steering_is_mirror_symmetric() {
    let mut left = rolling_sedan(15.0);
    let mut right = rolling_sedan(15.0);
    drive(&mut left, 3.0, 1.0 / 60.0, 0.3, 0.0, 0.4);
    drive(&mut right, 3.0, 1.0 / 60.0, 0.3, 0.0, -0.4);
    assert!(left.orientation[2] > 0.3, "yaw {}", left.orientation[2]);
    assert!((left.position[0] - right.position[0]).abs() < 1e-9);
    assert!((left.position[2] + right.position[2]).abs() < 1e-9);
    assert!((left.orientation[2] + right.orientation[2]).abs() < 1e-9);
    assert!((left.angular_velocity[2] + right.angular_velocity[2]).abs() < 1e-9);
    assert!((left.forward_speed() - right.forward_speed()).abs() < 1e-9);
}

#[test]
fn test_vehicle_tyres_damp_lateral_sliding() {
    let mut v = rolling_sedan(10.0);
    v.velocity[2] = 2.0;
    drive(&mut v, 1.0, 1.0 / 60.0, 0.0, 0.0, 0.0);
    assert!(
        v.lateral_speed().abs() < 0.05,
        "lateral speed {}",
        v.lateral_speed()
    );
    assert!(v.orientation[2].abs() < 0.05, "yaw {}", v.orientation[2]);
}

#[test]
fn test_vehicle_reversing_steers_the_other_way() {
    let mut v = rolling_sedan(-5.0);
    drive(&mut v, 3.0, 1.0 / 60.0, 0.0, 0.0, 0.2);
    let kinematic = kinematic_yaw_rate(&v);
    assert!(kinematic < 0.0);
    assert!(v.forward_speed() < -4.9);
    assert!(
        (v.angular_velocity[2] - kinematic).abs() < 0.05 * kinematic.abs(),
        "yaw rate {} vs kinematic {kinematic}",
        v.angular_velocity[2]
    );
}

#[test]
fn test_vehicle_with_linear_tyres_turns_the_same_way() {
    let mut v = rolling_sedan(10.0);
    let tire = linear_tire(60_000.0, 8.0);
    v.tires = [tire.clone(), tire.clone(), tire.clone(), tire];
    drive(&mut v, 3.0, 1.0 / 60.0, 0.0, 0.0, 0.1);
    let kinematic = kinematic_yaw_rate(&v);
    assert!(
        (v.angular_velocity[2] - kinematic).abs() < 0.1 * kinematic,
        "yaw rate {} vs kinematic {kinematic}",
        v.angular_velocity[2]
    );
    assert!(v.lateral_speed() > 0.0 && v.lateral_speed() < 0.3);
}

#[test]
fn test_vehicle_never_gains_energy_without_the_engine() {
    let mut v = rolling_sedan(15.0);
    let dt = 1.0 / 60.0;
    let start = vehicle_energy(&v);
    let mut previous = start;
    for step in 0..600 {
        // Coast while steering left and right, with a spell of braking.
        let steer = if step < 300 { 0.5 } else { -0.8 };
        let brake = if (200..260).contains(&step) { 0.4 } else { 0.0 };
        v.step(dt, 0.0, brake, steer);
        let energy = vehicle_energy(&v);
        assert!(
            energy <= previous + 1e-7 * start,
            "energy grew at step {step}: {previous} -> {energy}"
        );
        previous = energy;
    }
    assert!(previous < start);
}

#[test]
fn test_vehicle_step_ignores_non_positive_and_non_finite_time_steps() {
    let mut v = rolling_sedan(10.0);
    let before = format!("{v:?}");
    for dt in [0.0, -0.016, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        v.step(dt, 1.0, 0.0, 0.5);
        assert_eq!(format!("{v:?}"), before, "dt = {dt}");
    }
}

#[test]
fn test_vehicle_stays_finite_for_extreme_steps_and_inputs() {
    let mut v = PyVehicle::racing_car();
    for (k, &dt) in [1e-9, 1e-4, 0.016, 0.5, 5.0].iter().enumerate() {
        v.step(dt, 1.0, (k % 2) as f64, -1.0 + 0.5 * k as f64);
        assert_vehicle_finite(&v);
    }
    v.step(0.016, f64::NAN, f64::INFINITY, f64::NEG_INFINITY);
    assert_vehicle_finite(&v);
    // A long run with erratic controls.
    let mut v = PyVehicle::racing_car();
    for k in 0..2000 {
        let t = k as f64;
        let throttle = 0.5 + 0.5 * (0.37 * t).sin();
        let brake = (0.5 * (0.11 * t).cos()).max(0.0);
        v.step(1.0 / 30.0, throttle, brake, (0.05 * t).sin());
    }
    assert_vehicle_finite(&v);
}

#[test]
fn test_vehicle_survives_a_degenerate_configuration() {
    let mut v = PyVehicle::sedan();
    v.mass = 0.0;
    v.inertia = [0.0; 3];
    for wheel in &mut v.wheels {
        wheel.inertia = 0.0;
    }
    v.step(0.016, 1.0, 0.0, 0.5);
    assert_vehicle_finite(&v);
}

#[test]
fn test_static_wheel_loads_balance_the_weight_and_its_moments() {
    for v in [PyVehicle::sedan(), PyVehicle::racing_car()] {
        let loads = v.static_wheel_loads();
        assert_eq!(loads.len(), 4);
        let weight = v.mass * 9.81;
        assert!((loads.iter().sum::<f64>() - weight).abs() < 1e-9 * weight);
        let pitch: f64 = loads
            .iter()
            .zip(&v.wheel_positions)
            .map(|(l, p)| l * p[0])
            .sum();
        let roll: f64 = loads
            .iter()
            .zip(&v.wheel_positions)
            .map(|(l, p)| l * p[2])
            .sum();
        assert!(pitch.abs() < 1e-9 * weight && roll.abs() < 1e-9 * weight);
    }
    // The axle nearer to the centre of mass carries more of the weight ...
    let mut v = PyVehicle::sedan();
    let weight = v.mass * 9.81;
    v.wheel_positions[0][0] = 1.0;
    v.wheel_positions[1][0] = 1.0;
    v.wheel_positions[2][0] = -2.0;
    v.wheel_positions[3][0] = -2.0;
    let loads = v.static_wheel_loads();
    assert!((loads[0] + loads[1] - weight * 2.0 / 3.0).abs() < 1e-9 * weight);
    assert!((loads[2] + loads[3] - weight / 3.0).abs() < 1e-9 * weight);
    // ... and so does the side nearer to it.
    v.wheel_positions[0][2] = 0.5;
    v.wheel_positions[2][2] = 0.5;
    let loads = v.static_wheel_loads();
    assert!((loads[0] / (loads[0] + loads[1]) - 0.6).abs() < 1e-12);
    assert!((loads[2] / (loads[2] + loads[3]) - 0.6).abs() < 1e-12);
}

#[test]
fn test_vehicle_presets_have_consistent_geometry() {
    for v in [PyVehicle::sedan(), PyVehicle::racing_car()] {
        let wheelbase = v.wheel_positions[0][0] - v.wheel_positions[2][0];
        let track = v.wheel_positions[0][2] - v.wheel_positions[1][2];
        assert!(
            (wheelbase - v.steering.wheelbase).abs() < 1e-12,
            "wheelbase {wheelbase}"
        );
        assert!(
            (track - v.steering.track_width).abs() < 1e-12,
            "track {track}"
        );
        // The left wheels (FL, RL) sit on the positive lateral side.
        assert!(v.wheel_positions[0][2] > 0.0 && v.wheel_positions[1][2] < 0.0);
        assert!(v.wheel_positions[2][2] > 0.0 && v.wheel_positions[3][2] < 0.0);
    }
}
