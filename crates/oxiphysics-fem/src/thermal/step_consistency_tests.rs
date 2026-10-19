//! Regression tests for the 1-D thermal time integrators.
//!
//! The explicit and implicit Euler steps discretise the same semi-discrete
//! equation `C dT/dt = Q - K T + (convection / radiation load)`, so for a tiny
//! time step they must agree.  These tests pin the convective boundary term:
//! the loss `h A (T - T_inf)` must be counted exactly once in
//! [`ThermalMesh1D::step_explicit`].

use super::types::ThermalMesh1D;

/// Uniform rod, `n` nodes, unit length / conductivity / density / heat capacity,
/// initially at `t_init`, with a convective boundary at the last node.
fn convective_rod(n: usize, t_init: f64, h: f64, t_inf: f64) -> ThermalMesh1D {
    let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
    for t in mesh.temperatures.iter_mut() {
        *t = t_init;
    }
    mesh.add_convection_bc(n - 1, h, 1.0, t_inf);
    mesh
}

#[test]
fn explicit_single_step_convective_loss_counted_once() {
    // Uniform field: conduction contributes nothing, so the end node only loses
    // heat by convection:  C_end dT = -dt * h A (T - T_inf).
    let n = 5;
    let (t0, h, t_inf) = (100.0, 10.0, 20.0);
    let mut mesh = convective_rod(n, t0, h, t_inf);
    let c_end = mesh.assemble_capacitance_vector()[n - 1];
    let dt = 1.0e-4;
    mesh.step_explicit(dt, &vec![0.0; n]);
    let expected_end = t0 - dt * h * 1.0 * (t0 - t_inf) / c_end;
    assert!(
        (mesh.temperatures[n - 1] - expected_end).abs() < 1e-12,
        "end node: got {}, expected {expected_end} (a doubled loss would give {})",
        mesh.temperatures[n - 1],
        t0 - 2.0 * dt * h * (t0 - t_inf) / c_end
    );
    for (i, &t) in mesh.temperatures.iter().enumerate().take(n - 1) {
        assert!(
            (t - t0).abs() < 1e-12,
            "node {i} must stay at {t0} (uniform field), got {t}"
        );
    }
}

#[test]
fn explicit_and_implicit_agree_for_tiny_steps_with_end_convection() {
    let n = 5;
    let dt = 1.0e-6;
    let mut explicit = convective_rod(n, 100.0, 10.0, 20.0);
    let mut implicit = convective_rod(n, 100.0, 10.0, 20.0);
    let q = vec![0.0; n];
    for _ in 0..100 {
        explicit.step_explicit(dt, &q);
        implicit.step_implicit(dt, &q);
    }
    // The convective node has cooled measurably, so the comparison is not vacuous.
    assert!(
        100.0 - implicit.temperatures[n - 1] > 0.5,
        "end node should have cooled: {}",
        implicit.temperatures[n - 1]
    );
    for i in 0..n {
        let (te, ti) = (explicit.temperatures[i], implicit.temperatures[i]);
        assert!(
            (te - ti).abs() <= 1e-6 * ti.abs(),
            "node {i}: explicit {te} vs implicit {ti}"
        );
    }
}

#[test]
fn explicit_and_implicit_agree_with_heat_source_and_dirichlet_end() {
    let n = 6;
    let dt = 1.0e-6;
    let build = || {
        let mut mesh = ThermalMesh1D::new_uniform(n, 2.0, 3.0, 2.0, 1.5);
        for (i, t) in mesh.temperatures.iter_mut().enumerate() {
            *t = 80.0 - 7.0 * i as f64;
        }
        mesh.set_temperature_bc(0, 80.0);
        mesh.add_convection_bc(n - 1, 25.0, 0.5, 15.0);
        mesh
    };
    let mut explicit = build();
    let mut implicit = build();
    let mut q = vec![0.0; n];
    q[2] = 40.0;
    q[3] = 25.0;
    for _ in 0..100 {
        explicit.step_explicit(dt, &q);
        implicit.step_implicit(dt, &q);
    }
    for i in 0..n {
        let (te, ti) = (explicit.temperatures[i], implicit.temperatures[i]);
        assert!(
            (te - ti).abs() <= 1e-6 * ti.abs(),
            "node {i}: explicit {te} vs implicit {ti}"
        );
    }
}

#[test]
fn explicit_time_marching_reaches_analytic_convective_steady_state() {
    // Rod held at 100 at x = 0 with convection (h A = 10, T_inf = 20) at x = L,
    // k A / L = 1:  conduction (100 - T_L) balances convection 10 (T_L - 20),
    // hence T_L = 300 / 11 and a linear profile in between.
    let n = 5;
    let mut mesh = convective_rod(n, 100.0, 10.0, 20.0);
    mesh.set_temperature_bc(0, 100.0);
    let q = vec![0.0; n];
    // Stability: dt < 2 C_end / K_end = 2 * 0.125 / 14.
    for _ in 0..12_000 {
        mesh.step_explicit(0.005, &q);
    }
    let t_end_exact = 300.0 / 11.0;
    for (i, &t) in mesh.temperatures.iter().enumerate() {
        let x = i as f64 / (n - 1) as f64;
        let expected = 100.0 - (100.0 - t_end_exact) * x;
        assert!(
            (t - expected).abs() < 1e-7,
            "node {i}: got {t}, expected {expected}"
        );
    }
    // The direct steady-state solve is the same discrete problem.
    let mut direct = convective_rod(n, 0.0, 10.0, 20.0);
    direct.set_temperature_bc(0, 100.0);
    direct.steady_state(&q);
    for i in 0..n {
        assert!(
            (mesh.temperatures[i] - direct.temperatures[i]).abs() < 1e-7,
            "node {i}: explicit {} vs steady_state {}",
            mesh.temperatures[i],
            direct.temperatures[i]
        );
    }
}
