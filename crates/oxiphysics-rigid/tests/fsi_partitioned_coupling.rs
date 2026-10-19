// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `PartitionedFsiCoupling`.
//!
//! The coupler advances its structure with the Newmark average-acceleration
//! scheme from a committed state and relaxes the interface displacement with
//! Aitken's formula.  It used to solve `F / k_eff` without any history, which
//! cannot represent a free vibration, and it ignored Aitken altogether.

use oxiphysics_rigid::Error;
use oxiphysics_rigid::fluid_structure::PartitionedFsiCoupling;
use std::f64::consts::PI;

const MASS: f64 = 100.0;
const STIFFNESS: f64 = 1.0e4;

fn natural_frequency() -> f64 {
    (STIFFNESS / MASS).sqrt()
}

/// Newmark effective stiffness for β = 1/4, γ = 1/2.
fn effective_stiffness(m: f64, k: f64, c: f64, dt: f64) -> f64 {
    k + 0.5 / (0.25 * dt) * c + m / (0.25 * dt * dt)
}

// ── Newmark structure ───────────────────────────────────────────────────────

#[test]
fn undamped_free_vibration_follows_the_cosine_and_conserves_energy() {
    let omega = natural_frequency();
    let period = 2.0 * PI / omega;
    let dt = period / 100.0;
    let d0 = 0.02;
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, 0.0);
    coupler.displacement_tol = 1.0e-12;
    coupler.set_initial_state(d0, 0.0, 0.0);
    let energy0 = 0.5 * STIFFNESS * d0 * d0;

    for step in 1..=200 {
        let iterations = coupler
            .solve_step(dt, |_, _| 0.0)
            .expect("a force-free step converges");
        // The structure's answer does not depend on the iterate, so Aitken's
        // second factor is exactly one and the third residual vanishes.
        assert_eq!(iterations, 3, "step {step}");
        let t = f64::from(step) * dt;
        let displacement = coupler.step_displacement;
        let velocity = coupler.step_velocity;
        assert!(
            (displacement - d0 * (omega * t).cos()).abs() < 0.01 * d0,
            "step {step}: d = {displacement}"
        );
        assert!(
            (velocity + d0 * omega * (omega * t).sin()).abs() < 0.01 * d0 * omega,
            "step {step}: v = {velocity}"
        );
        let energy =
            0.5 * MASS * velocity * velocity + 0.5 * STIFFNESS * displacement * displacement;
        assert!(
            (energy - energy0).abs() < 1.0e-9 * energy0,
            "step {step}: energy {energy} drifted from {energy0}"
        );
    }
}

#[test]
fn damped_free_vibration_matches_the_analytic_decay() {
    let damping = 50.0;
    let omega = natural_frequency();
    let zeta = damping / (2.0 * (STIFFNESS * MASS).sqrt());
    let omega_d = omega * (1.0 - zeta * zeta).sqrt();
    let d0 = 0.02;
    let dt = 0.005;
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, damping);
    coupler.displacement_tol = 1.0e-12;
    coupler.set_initial_state(d0, 0.0, 0.0);
    for step in 1..=400 {
        coupler
            .solve_step(dt, |_, _| 0.0)
            .expect("a force-free step converges");
        let t = f64::from(step) * dt;
        let exact = d0
            * (-zeta * omega * t).exp()
            * ((omega_d * t).cos() + zeta * omega / omega_d * (omega_d * t).sin());
        assert!(
            (coupler.step_displacement - exact).abs() < 1.0e-3,
            "step {step}: d = {}, exact {exact}",
            coupler.step_displacement
        );
    }
}

#[test]
fn constant_load_settles_at_the_static_deflection() {
    let load = 500.0;
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, 50.0);
    coupler.displacement_tol = 1.0e-12;
    coupler.set_initial_state(0.0, 0.0, load);
    for _ in 0..7000 {
        coupler
            .solve_step(0.01, |_, _| load)
            .expect("constant loads converge");
    }
    assert!(
        (coupler.step_displacement - load / STIFFNESS).abs() < 1.0e-6,
        "d = {}",
        coupler.step_displacement
    );
    assert!(coupler.step_velocity.abs() < 1.0e-6);
}

#[test]
fn converged_state_satisfies_the_equation_of_motion() {
    let (m, k, c) = (MASS, STIFFNESS, 40.0);
    let fluid = |d: f64, v: f64| 120.0 - 800.0 * d - 30.0 * v;
    let mut coupler = PartitionedFsiCoupling::new(m, k, c);
    coupler.displacement_tol = 1.0e-13;
    coupler.set_initial_state(0.01, 0.3, fluid(0.01, 0.3));
    for step in 1..=20 {
        coupler.solve_step(0.01, fluid).expect("converges");
        let residual = m * coupler.step_acceleration
            + c * coupler.step_velocity
            + k * coupler.step_displacement
            - fluid(coupler.step_displacement, coupler.step_velocity);
        assert!(residual.abs() < 1.0e-7, "step {step}: residual {residual}");
    }
}

#[test]
fn newmark_relations_link_the_committed_states() {
    let (beta, gamma) = (0.25, 0.5);
    let dt = 0.02;
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, 30.0);
    coupler.displacement_tol = 1.0e-13;
    coupler.set_initial_state(0.004, -0.2, 90.0);
    for _ in 0..10 {
        let (d_n, v_n, a_n) = (
            coupler.step_displacement,
            coupler.step_velocity,
            coupler.step_acceleration,
        );
        coupler.solve_step(dt, |_, _| 250.0).expect("converges");
        let (d, v, a) = (
            coupler.step_displacement,
            coupler.step_velocity,
            coupler.step_acceleration,
        );
        let d_expected = d_n + dt * v_n + dt * dt * ((0.5 - beta) * a_n + beta * a);
        let v_expected = v_n + dt * ((1.0 - gamma) * a_n + gamma * a);
        assert!((d - d_expected).abs() < 1.0e-12, "{d} vs {d_expected}");
        assert!((v - v_expected).abs() < 1.0e-10, "{v} vs {v_expected}");
    }
}

#[test]
fn initial_state_fixes_the_initial_acceleration() {
    let mut coupler = PartitionedFsiCoupling::new(2.0, 50.0, 4.0);
    coupler.set_initial_state(0.1, -0.5, 7.0);
    // (7 - 4 * (-0.5) - 50 * 0.1) / 2 = 2
    assert!((coupler.step_acceleration - 2.0).abs() < 1.0e-12);
    assert!((coupler.interface_acceleration - 2.0).abs() < 1.0e-12);
    assert_eq!(coupler.step_displacement, 0.1);
    assert_eq!(coupler.step_velocity, -0.5);
    assert_eq!(coupler.interface_displacement, 0.1);
    assert_eq!(coupler.interface_velocity, -0.5);
    assert_eq!(coupler.fluid_force, 7.0);
}

// ── Aitken relaxation ───────────────────────────────────────────────────────

#[test]
fn aitken_relaxation_converges_where_a_fixed_factor_diverges() {
    let (m, k, c, dt) = (MASS, STIFFNESS, 20.0, 0.01);
    // The fluid pushes back with a stiffness far above the structure's own
    // effective stiffness: q = kf / k_eff is about 5.
    let kf = 2.0e7;
    let q = kf / effective_stiffness(m, k, c, dt);
    assert!(
        (1.0 - 0.5 * (1.0 + q)).abs() > 1.0,
        "set-up: a fixed factor of 0.5 must diverge, q = {q}"
    );
    assert!(
        (1.0 - (1.0 + q)).abs() > 1.0,
        "set-up: so must a factor of 1"
    );

    let d0 = 0.01;
    let mut coupler = PartitionedFsiCoupling::new(m, k, c);
    coupler.set_initial_state(d0, 0.1, -kf * d0);
    let mut residuals = Vec::new();
    for _ in 0..4 {
        let force = -kf * coupler.interface_displacement;
        residuals.push(coupler.sub_iterate(force, dt));
    }
    // The first factor is the initial relaxation, which makes the residual
    // grow; for a linear map Aitken's second factor is the optimal 1 / (1 + q)
    // and the third residual is zero.
    assert!(residuals[1] > residuals[0], "{residuals:?}");
    assert!(
        residuals[2] < 1.0e-9 * residuals[0],
        "Aitken must remove the residual: {residuals:?}"
    );
    assert!(residuals[3] < 1.0e-9 * residuals[0], "{residuals:?}");
}

#[test]
fn aitken_factor_is_the_optimal_one_for_a_linear_fluid() {
    let (m, k, c, dt) = (MASS, STIFFNESS, 20.0, 0.01);
    let kf = 2.0e7;
    let q = kf / effective_stiffness(m, k, c, dt);
    let mut coupler = PartitionedFsiCoupling::new(m, k, c);
    coupler.set_initial_state(0.01, 0.1, -kf * 0.01);
    assert!(coupler.previous_residual.is_none());
    coupler.sub_iterate(-kf * coupler.interface_displacement, dt);
    assert_eq!(coupler.omega, coupler.relaxation, "first factor");
    assert!(coupler.previous_residual.is_some());
    coupler.sub_iterate(-kf * coupler.interface_displacement, dt);
    assert!(
        (coupler.omega - 1.0 / (1.0 + q)).abs() < 1.0e-9,
        "omega = {}, optimal {}",
        coupler.omega,
        1.0 / (1.0 + q)
    );
}

#[test]
fn partitioned_scheme_matches_a_monolithic_added_stiffness() {
    let (m, k, c, dt) = (MASS, STIFFNESS, 20.0, 0.01);
    let kf = 2.0e7;
    let (d0, v0) = (0.01, 0.1);

    let mut partitioned = PartitionedFsiCoupling::new(m, k, c);
    partitioned.displacement_tol = 1.0e-12;
    partitioned.set_initial_state(d0, v0, -kf * d0);
    // The same structure with the fluid stiffness added to its own.
    let mut monolithic = PartitionedFsiCoupling::new(m, k + kf, c);
    monolithic.displacement_tol = 1.0e-12;
    monolithic.set_initial_state(d0, v0, 0.0);

    let mut most_iterations = 0;
    for step in 1..=100 {
        let used = partitioned
            .solve_step(dt, |d, _| -kf * d)
            .expect("Aitken converges");
        most_iterations = most_iterations.max(used);
        monolithic
            .solve_step(dt, |_, _| 0.0)
            .expect("a force-free step converges");
        let reference = monolithic.step_displacement;
        assert!(
            (partitioned.step_displacement - reference).abs()
                < 1.0e-9 * reference.abs().max(1.0e-6),
            "step {step}: {} vs {reference}",
            partitioned.step_displacement
        );
    }
    assert!(
        most_iterations <= 4,
        "needed {most_iterations} sub-iterations"
    );
}

#[test]
fn committing_a_step_restarts_the_aitken_history() {
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, 20.0);
    coupler.relaxation = 0.8;
    coupler.set_initial_state(0.0, 0.0, 0.0);
    assert_eq!(coupler.omega, 0.8);
    coupler.sub_iterate(100.0, 0.01);
    coupler.sub_iterate(100.0, 0.01);
    assert!(coupler.previous_residual.is_some());
    coupler.advance_step(0.01).expect("valid time step");
    assert!(coupler.previous_residual.is_none());
    assert_eq!(coupler.omega, 0.8);
    assert_eq!(coupler.step_displacement, coupler.interface_displacement);
    assert_eq!(coupler.step_velocity, coupler.interface_velocity);
    assert_eq!(coupler.step_acceleration, coupler.interface_acceleration);
}

// ── error handling ──────────────────────────────────────────────────────────

#[test]
fn invalid_time_steps_leave_the_state_untouched() {
    let mut coupler = PartitionedFsiCoupling::new(MASS, STIFFNESS, 20.0);
    coupler.set_initial_state(0.01, 0.2, 5.0);
    for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        let residual = coupler.sub_iterate(100.0, dt);
        assert!(residual.is_infinite(), "dt {dt}: residual {residual}");
        assert!(!coupler.is_converged(residual));
        assert_eq!(coupler.interface_displacement, 0.01);
        assert_eq!(coupler.interface_velocity, 0.2);
        assert!(coupler.previous_residual.is_none());
        assert!(matches!(coupler.advance_step(dt), Err(Error::General(_))));
        assert!(matches!(
            coupler.solve_step(dt, |_, _| 0.0),
            Err(Error::General(_))
        ));
    }
}

#[test]
fn unconverged_steps_are_reported_and_not_committed() {
    let (m, k, c, dt) = (MASS, STIFFNESS, 20.0, 0.01);
    let kf = 2.0e7;
    let mut coupler = PartitionedFsiCoupling::new(m, k, c);
    coupler.max_iterations = 1;
    coupler.set_initial_state(0.01, 0.1, -kf * 0.01);
    let before = (coupler.step_displacement, coupler.step_velocity);
    let result = coupler.solve_step(dt, |d, _| -kf * d);
    assert!(matches!(
        result,
        Err(Error::SolverNotConverged {
            max_iterations: 1,
            ..
        })
    ));
    assert_eq!((coupler.step_displacement, coupler.step_velocity), before);

    coupler.max_iterations = 0;
    assert!(matches!(
        coupler.solve_step(dt, |_, _| 0.0),
        Err(Error::SolverNotConverged { .. })
    ));
}

#[test]
fn structure_without_effective_stiffness_never_converges() {
    let mut coupler = PartitionedFsiCoupling::new(0.0, 0.0, 0.0);
    let residual = coupler.sub_iterate(10.0, 0.01);
    assert!(residual.is_infinite());
    assert_eq!(coupler.interface_displacement, 0.0);
}
