// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`SixDofIntegrator`]: the rate-damping derivatives
//! (`Cm_q`, `Cl_q`, `Cl_p`, `Cl_r`, `Cn_p`, `Cn_r`) and the `Ixz` coupling of
//! the rigid-body rotation equations.

use super::*;

const IXX: f64 = 2_000_000.0;
const IYY: f64 = 5_000_000.0;
const IZZ: f64 = 8_000_000.0;
const IXZ: f64 = 250_000.0;
const S_REF: f64 = 120.0;
const B_REF: f64 = 28.0;
const C_MAC: f64 = 4.0;
const MASS: f64 = 70_000.0;

fn integrator(i_xz: f64) -> SixDofIntegrator {
    SixDofIntegrator {
        mass: MASS,
        inertia: [IXX, IYY, IZZ],
        i_xz,
        s_ref: S_REF,
        b_ref: B_REF,
        c_mac: C_MAC,
        long_derivs: LongitudinalDerivatives::typical_transport(),
        lat_derivs: LateralDerivatives::typical_transport(),
        polar: LiftDragPolar::typical_transport(),
    }
}

/// The same airframe with every rate-damping derivative switched off.
fn undamped(i_xz: f64) -> SixDofIntegrator {
    let mut model = integrator(i_xz);
    model.long_derivs.cm_q = 0.0;
    model.long_derivs.cl_q = 0.0;
    model.lat_derivs.cl_p = 0.0;
    model.lat_derivs.cl_r = 0.0;
    model.lat_derivs.cn_p = 0.0;
    model.lat_derivs.cn_r = 0.0;
    model
}

fn state_with_rates(speed: f64, rates: [f64; 3]) -> AircraftState6Dof {
    let mut state = AircraftState6Dof::level_flight(speed);
    state.angular_rates = rates;
    state
}

/// Solve the 3x3 system `a x = b` by Cramer's rule.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    let mut x = [0.0; 3];
    for (col, x_col) in x.iter_mut().enumerate() {
        let mut m = a;
        for (row, m_row) in m.iter_mut().enumerate() {
            m_row[col] = b[row];
        }
        *x_col = det(m) / d;
    }
    x
}

fn mat_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dynamic_pressure(speed: f64) -> f64 {
    0.5 * RHO_SL * speed * speed
}

#[test]
fn angular_acceleration_solves_euler_equation_with_full_inertia_tensor() {
    // J ω̇ = M − ω × (J ω) with J = [[Ixx, 0, −Ixz], [0, Iyy, 0], [−Ixz, 0, Izz]].
    let cases = [
        ([1.0e5, -2.0e5, 3.0e5], [0.10, -0.20, 0.30]),
        ([0.0, 0.0, 0.0], [0.5, 0.4, -0.7]),
        ([-4.0e6, 1.0e6, 2.5e6], [-0.3, 0.05, 0.15]),
        ([2.0e5, 0.0, -9.0e5], [0.0, 1.2, 0.0]),
    ];
    for i_xz in [0.0, IXZ, -IXZ] {
        let model = integrator(i_xz);
        let tensor = [[IXX, 0.0, -i_xz], [0.0, IYY, 0.0], [-i_xz, 0.0, IZZ]];
        for (moments, rates) in cases {
            let gyro = cross(rates, mat_vec(tensor, rates));
            let rhs = [
                moments[0] - gyro[0],
                moments[1] - gyro[1],
                moments[2] - gyro[2],
            ];
            let expected = solve3(tensor, rhs);
            let got = model.angular_acceleration(moments, rates);
            for axis in 0..3 {
                assert!(
                    (got[axis] - expected[axis]).abs() < 1e-9 * (1.0 + expected[axis].abs()),
                    "Ixz = {i_xz}, M = {moments:?}, w = {rates:?}, axis {axis}: {} vs {}",
                    got[axis],
                    expected[axis]
                );
            }
        }
    }
}

#[test]
fn angular_acceleration_reduces_to_diagonal_equations_without_ixz() {
    let model = integrator(0.0);
    let (moments, rates) = ([3.0e5, -1.0e5, 2.0e5], [0.2, -0.1, 0.35]);
    let [p_dot, q_dot, r_dot] = model.angular_acceleration(moments, rates);
    let [p, q, r] = rates;
    assert!((p_dot - (moments[0] + (IYY - IZZ) * q * r) / IXX).abs() < 1e-12);
    assert!((q_dot - (moments[1] + (IZZ - IXX) * p * r) / IYY).abs() < 1e-12);
    assert!((r_dot - (moments[2] + (IXX - IYY) * p * q) / IZZ).abs() < 1e-12);
}

#[test]
fn roll_and_yaw_moments_couple_through_ixz() {
    let gamma = IXX * IZZ - IXZ * IXZ;
    let model = integrator(IXZ);
    // A pure rolling moment also accelerates the yaw rate.
    let [p_dot, q_dot, r_dot] = model.angular_acceleration([1.0e5, 0.0, 0.0], [0.0; 3]);
    assert!((p_dot - IZZ * 1.0e5 / gamma).abs() < 1e-12);
    assert!((r_dot - IXZ * 1.0e5 / gamma).abs() < 1e-12);
    assert!(r_dot > 0.0 && q_dot == 0.0);
    // A pure yawing moment also accelerates the roll rate.
    let [p_dot, q_dot, r_dot] = model.angular_acceleration([0.0, 0.0, 1.0e5], [0.0; 3]);
    assert!((p_dot - IXZ * 1.0e5 / gamma).abs() < 1e-12);
    assert!((r_dot - IXX * 1.0e5 / gamma).abs() < 1e-12);
    assert!(p_dot > 0.0 && q_dot == 0.0);
    // Without the product of inertia the axes are independent.
    let [p_dot, _, r_dot] = integrator(0.0).angular_acceleration([1.0e5, 0.0, 0.0], [0.0; 3]);
    assert!((p_dot - 1.0e5 / IXX).abs() < 1e-12);
    assert_eq!(r_dot, 0.0);
}

#[test]
fn flight_model_applies_ixz_to_rudder_induced_moments() {
    // Rudder only, from zero rates, beta = 0: L = q S b Cl_dr dr, N = q S b Cn_dr dr.
    let speed = 150.0;
    let delta_r = 0.05;
    let lat = LateralDerivatives::typical_transport();
    let qsb = dynamic_pressure(speed) * S_REF * B_REF;
    let roll_moment = qsb * lat.cl_dr * delta_r;
    let yaw_moment = qsb * lat.cn_dr * delta_r;
    let state = AircraftState6Dof::level_flight(speed);

    let without = integrator(0.0).derivatives(&state, 0.0, 0.0, 0.0, delta_r, RHO_SL);
    assert!((without.angular_rates[0] - roll_moment / IXX).abs() < 1e-12);
    assert!((without.angular_rates[2] - yaw_moment / IZZ).abs() < 1e-12);

    let with = integrator(IXZ).derivatives(&state, 0.0, 0.0, 0.0, delta_r, RHO_SL);
    let gamma = IXX * IZZ - IXZ * IXZ;
    let p_expected = (IZZ * roll_moment + IXZ * yaw_moment) / gamma;
    let r_expected = (IXZ * roll_moment + IXX * yaw_moment) / gamma;
    assert!((with.angular_rates[0] - p_expected).abs() < 1e-12);
    assert!((with.angular_rates[2] - r_expected).abs() < 1e-12);
    assert!(
        (with.angular_rates[0] - without.angular_rates[0]).abs() > 1e-6,
        "Ixz must change the roll acceleration"
    );
    assert!((with.angular_rates[2] - without.angular_rates[2]).abs() > 1e-6);
}

#[test]
fn aero_coefficients_apply_rate_derivatives() {
    let model = integrator(IXZ);
    let (p_hat, q_hat, r_hat) = (0.02, -0.03, 0.04);
    let damped = model.aero_coefficients(0.05, 0.02, p_hat, q_hat, r_hat, 0.01, 0.02, 0.03);
    let static_only = model.aero_coefficients(0.05, 0.02, 0.0, 0.0, 0.0, 0.01, 0.02, 0.03);
    let long = &model.long_derivs;
    let lat = &model.lat_derivs;
    assert!((damped.cl - static_only.cl - long.cl_q * q_hat).abs() < 1e-12);
    assert!((damped.cm - static_only.cm - long.cm_q * q_hat).abs() < 1e-12);
    assert!(
        (damped.c_roll - static_only.c_roll - (lat.cl_p * p_hat + lat.cl_r * r_hat)).abs() < 1e-12
    );
    assert!((damped.cn - static_only.cn - (lat.cn_p * p_hat + lat.cn_r * r_hat)).abs() < 1e-12);
    // Side force has no rate derivative in the model.
    assert_eq!(damped.cy, static_only.cy);
}

#[test]
fn pitch_rate_produces_pitch_damping_moment_and_lift() {
    let speed = 150.0;
    let q_rate = 0.1;
    let state = state_with_rates(speed, [0.0, q_rate, 0.0]);
    let damped = integrator(IXZ).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    let reference = undamped(IXZ).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);

    // alpha = 0 and no deflections: without damping there is no pitch moment.
    assert_eq!(reference.angular_rates[1], 0.0);
    let q_hat = q_rate * C_MAC / (2.0 * speed);
    let long = LongitudinalDerivatives::typical_transport();
    let q_bar = dynamic_pressure(speed);
    let expected_q_dot = q_bar * S_REF * C_MAC * long.cm_q * q_hat / IYY;
    assert!(expected_q_dot < 0.0);
    assert!(
        (damped.angular_rates[1] - expected_q_dot).abs() < 1e-12,
        "q_dot = {}, expected {expected_q_dot}",
        damped.angular_rates[1]
    );
    // Cl_q adds lift (body z points down, so w_dot decreases).
    let expected_dw = -q_bar * S_REF * long.cl_q * q_hat / MASS;
    let dw = damped.velocity_body[2] - reference.velocity_body[2];
    assert!(
        (dw - expected_dw).abs() < 1e-12,
        "dw = {dw}, expected {expected_dw}"
    );
}

#[test]
fn roll_rate_produces_roll_damping_and_adverse_yaw() {
    let speed = 150.0;
    let p_rate = 0.1;
    let state = state_with_rates(speed, [p_rate, 0.0, 0.0]);
    let damped = integrator(0.0).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    let reference = undamped(0.0).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    assert_eq!(reference.angular_rates, [0.0, 0.0, 0.0]);

    let p_hat = p_rate * B_REF / (2.0 * speed);
    let lat = LateralDerivatives::typical_transport();
    let qsb = dynamic_pressure(speed) * S_REF * B_REF;
    let expected_p_dot = qsb * lat.cl_p * p_hat / IXX;
    let expected_r_dot = qsb * lat.cn_p * p_hat / IZZ;
    assert!(expected_p_dot < 0.0, "roll damping opposes the roll rate");
    assert!((damped.angular_rates[0] - expected_p_dot).abs() < 1e-12);
    assert!((damped.angular_rates[2] - expected_r_dot).abs() < 1e-12);
}

#[test]
fn yaw_rate_produces_yaw_damping_and_roll_due_to_yaw() {
    let speed = 150.0;
    let r_rate = 0.1;
    let state = state_with_rates(speed, [0.0, 0.0, r_rate]);
    let damped = integrator(0.0).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    let reference = undamped(0.0).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    assert_eq!(reference.angular_rates, [0.0, 0.0, 0.0]);

    let r_hat = r_rate * B_REF / (2.0 * speed);
    let lat = LateralDerivatives::typical_transport();
    let qsb = dynamic_pressure(speed) * S_REF * B_REF;
    let expected_r_dot = qsb * lat.cn_r * r_hat / IZZ;
    let expected_p_dot = qsb * lat.cl_r * r_hat / IXX;
    assert!(expected_r_dot < 0.0, "yaw damping opposes the yaw rate");
    assert!((damped.angular_rates[2] - expected_r_dot).abs() < 1e-12);
    assert!((damped.angular_rates[0] - expected_p_dot).abs() < 1e-12);
}

#[test]
fn rate_terms_are_disabled_at_very_low_airspeed() {
    // Below 1 m/s the non-dimensional rates are undefined; the model must not
    // produce non-finite values.
    let state = state_with_rates(0.5, [0.1, 0.1, 0.1]);
    let d = integrator(IXZ).derivatives(&state, 0.0, 0.0, 0.0, 0.0, RHO_SL);
    assert!(d.angular_rates.iter().all(|v| v.is_finite()));
    assert!(d.velocity_body.iter().all(|v| v.is_finite()));
}

/// Largest |Δq| between `t_from` and `t_end` of the response to an initial
/// pitch-rate impulse `q0`, measured against an unperturbed reference run so
/// that the slow drift of the (untrimmed) flight path cancels out.
fn pitch_impulse_response_peak(model: &SixDofIntegrator, q0: f64, t_from: f64, t_end: f64) -> f64 {
    let dt = 0.005;
    let mut reference = AircraftState6Dof::level_flight(150.0);
    let mut perturbed = state_with_rates(150.0, [0.0, q0, 0.0]);
    let steps = (t_end / dt).round() as usize;
    let mut peak = 0.0_f64;
    for step in 0..steps {
        reference = model.rk4_step(&reference, 0.0, 0.0, 0.0, 0.0, RHO_SL, dt);
        perturbed = model.rk4_step(&perturbed, 0.0, 0.0, 0.0, 0.0, RHO_SL, dt);
        if (step + 1) as f64 * dt >= t_from {
            let dq = perturbed.angular_rates[1] - reference.angular_rates[1];
            peak = peak.max(dq.abs());
        }
    }
    peak
}

#[test]
fn pitch_rate_impulse_decays_faster_with_damping_derivatives() {
    let q0 = 0.1;
    let damped = integrator(IXZ);
    let reference = undamped(IXZ);
    // Both airframes start with the full impulse ...
    let early_damped = pitch_impulse_response_peak(&damped, q0, 0.0, 1.0);
    let early_reference = pitch_impulse_response_peak(&reference, q0, 0.0, 1.0);
    assert!((early_damped - q0).abs() < 0.01 * q0);
    assert!((early_reference - q0).abs() < 0.01 * q0);
    // ... but with Cm_q / Cl_q the short-period oscillation dies out faster.
    let mid_damped = pitch_impulse_response_peak(&damped, q0, 4.0, 8.0);
    let mid_reference = pitch_impulse_response_peak(&reference, q0, 4.0, 8.0);
    assert!(mid_damped.is_finite() && mid_reference.is_finite());
    assert!(
        mid_damped < 0.1 * q0,
        "the impulse must have decayed: {mid_damped} vs {q0}"
    );
    assert!(
        mid_damped < 0.6 * mid_reference,
        "peak |dq| in 4..8 s with damping {mid_damped} must be well below the undamped {mid_reference}"
    );
}

#[test]
fn roll_response_to_aileron_follows_first_order_roll_subsidence() {
    // With Cl_p the roll rate approaches p_ss = -(Cl_da da / Cl_p) 2V / b with
    // the roll-subsidence time constant tau = Ixx / (q S b Cl_p b / 2V);
    // without it the roll acceleration is constant and the rate keeps growing.
    let speed = 150.0;
    let delta_a = 0.05;
    let t_end = 0.5;
    let lat = LateralDerivatives::typical_transport();
    let p_ss = -(lat.cl_da * delta_a) / lat.cl_p * 2.0 * speed / B_REF;
    let qsb = dynamic_pressure(speed) * S_REF * B_REF;
    let tau = IXX / (qsb * lat.cl_p.abs() * B_REF / (2.0 * speed));
    let expected = p_ss * (1.0 - (-t_end / tau).exp());
    let simulate = |model: &SixDofIntegrator| {
        let dt = 0.005;
        let mut state = AircraftState6Dof::level_flight(speed);
        for _ in 0..(t_end / dt).round() as usize {
            state = model.rk4_step(&state, 0.0, 0.0, delta_a, 0.0, RHO_SL, dt);
        }
        state.angular_rates[0]
    };
    let with_damping = simulate(&integrator(IXZ));
    let without = simulate(&undamped(IXZ));
    assert!(
        (with_damping - expected).abs() < 0.03 * expected,
        "roll rate after {t_end} s = {with_damping}, first-order theory {expected}"
    );
    assert!(
        without > 1.15 * expected,
        "without damping the roll rate keeps accelerating: {without} vs {expected}"
    );
}
