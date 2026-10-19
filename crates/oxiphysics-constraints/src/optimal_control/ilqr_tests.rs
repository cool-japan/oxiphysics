//! Tests of the iLQR / DDP solver.
//!
//! The oracles are independent of the implementation: the finite-horizon LQR
//! gains for linear systems, a dense Newton method on the total cost for
//! quadratic problems with cross terms, closed-form one-step Newton updates,
//! and a nested exact value-function minimisation for the multi-stage DDP
//! second-order terms.

use super::ilqr::{
    DynamicsOrder, IlqrSettings, dynamics_derivatives, gradient_and_hessian, policy_at,
};
use super::types::{DifferentialDynamicProgramming, FiniteHorizonLqr, IlqrSolver};
use crate::dense_linalg::lu_solve;

const DT: f64 = 0.1;

fn double_integrator(x: &[f64], u: &[f64]) -> Vec<f64> {
    vec![x[0] + DT * x[1] + 0.5 * DT * DT * u[0], x[1] + DT * u[0]]
}

/// `l = 10 |x|² + 0.1 u²`, `φ = 10 |x|²` (Q = 10 I, R = 0.1, Qf = 10 I).
fn lqr_running(x: &[f64], u: &[f64]) -> f64 {
    10.0 * (x[0] * x[0] + x[1] * x[1]) + 0.1 * u[0] * u[0]
}

fn lqr_terminal(x: &[f64]) -> f64 {
    10.0 * (x[0] * x[0] + x[1] * x[1])
}

/// Total cost of the open-loop control sequence for an arbitrary problem.
fn total_cost<F, G, H>(
    x0: &[f64],
    u_seq: &[f64],
    m: usize,
    dynamics: F,
    running: G,
    terminal: H,
) -> f64
where
    F: Fn(&[f64], &[f64]) -> Vec<f64>,
    G: Fn(&[f64], &[f64]) -> f64,
    H: Fn(&[f64]) -> f64,
{
    let mut x = x0.to_vec();
    let mut cost = 0.0;
    for u in u_seq.chunks(m) {
        cost += running(&x, u);
        x = dynamics(&x, u);
    }
    cost + terminal(&x)
}

/// Dense Newton step on the total cost: `u - H⁻¹ g` with central differences.
fn newton_step(u: &[f64], cost: &dyn Fn(&[f64]) -> f64) -> Vec<f64> {
    let (grad, hess) = gradient_and_hessian(&|z: &[f64]| cost(z), u);
    let step = lu_solve(&hess, u.len(), &grad, 1).expect("regular Hessian");
    u.iter().zip(step.iter()).map(|(a, b)| a - b).collect()
}

#[test]
fn finite_difference_derivatives_match_analytic_ones() {
    // f(z) = sin(z0) z1² + exp(0.3 z0 z2) + z2⁴
    let f = |z: &[f64]| z[0].sin() * z[1] * z[1] + (0.3 * z[0] * z[2]).exp() + z[2].powi(4);
    let z = [0.7, -1.3, 0.4];
    let (grad, hess) = gradient_and_hessian(&f, &z);
    let e = (0.3 * z[0] * z[2]).exp();
    let expected_grad = [
        z[0].cos() * z[1] * z[1] + 0.3 * z[2] * e,
        2.0 * z[0].sin() * z[1],
        0.3 * z[0] * e + 4.0 * z[2].powi(3),
    ];
    let expected_hess = [
        -z[0].sin() * z[1] * z[1] + (0.3 * z[2]).powi(2) * e,
        2.0 * z[0].cos() * z[1],
        0.3 * e + 0.3 * z[2] * 0.3 * z[0] * e,
        2.0 * z[0].cos() * z[1],
        2.0 * z[0].sin(),
        0.0,
        0.3 * e + 0.3 * z[2] * 0.3 * z[0] * e,
        0.0,
        (0.3 * z[0]).powi(2) * e + 12.0 * z[2] * z[2],
    ];
    for i in 0..3 {
        assert!((grad[i] - expected_grad[i]).abs() < 1e-7, "grad[{i}]");
    }
    for i in 0..9 {
        assert!(
            (hess[i] - expected_hess[i]).abs() < 1e-6,
            "hess[{i}]: {} vs {}",
            hess[i],
            expected_hess[i]
        );
    }
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(
                hess[i * 3 + j],
                hess[j * 3 + i],
                "Hessian must be symmetric"
            );
        }
    }
}

#[test]
fn dynamics_jacobian_and_hessians_match_analytic_ones() {
    // f1 = z0 + 0.1 sin(z1) z2,   f2 = z0 z1 + z2²    (z = [x0, x1, u])
    let f = |z: &[f64]| vec![z[0] + 0.1 * z[1].sin() * z[2], z[0] * z[1] + z[2] * z[2]];
    let z = [0.5, 1.1, -0.8];
    let (jac, hess) = dynamics_derivatives(&f, &z, 2, true);
    let expected_jac = [
        1.0,
        0.1 * z[1].cos() * z[2],
        0.1 * z[1].sin(),
        z[1],
        z[0],
        2.0 * z[2],
    ];
    for (i, want) in expected_jac.iter().enumerate() {
        assert!((jac[i] - want).abs() < 1e-7, "jac[{i}]");
    }
    let h1 = [
        0.0,
        0.0,
        0.0,
        0.0,
        -0.1 * z[1].sin() * z[2],
        0.1 * z[1].cos(),
        0.0,
        0.1 * z[1].cos(),
        0.0,
    ];
    let h2 = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0];
    for i in 0..9 {
        assert!((hess[0][i] - h1[i]).abs() < 1e-6, "f1 hessian[{i}]");
        assert!((hess[1][i] - h2[i]).abs() < 1e-6, "f2 hessian[{i}]");
    }
    // The first-order variant skips the Hessians.
    let (jac1, hess1) = dynamics_derivatives(&f, &z, 2, false);
    assert_eq!(jac1, jac);
    assert!(hess1.is_empty());
}

#[test]
fn linear_quadratic_problem_reproduces_finite_horizon_lqr_gains() {
    let horizon = 12;
    let (a, b) = (vec![1.0, DT, 0.0, 1.0], vec![0.5 * DT * DT, DT]);
    let x0 = [1.0, -0.5];
    let solver = IlqrSolver::new(2, 1, horizon, 10, 1e-9, 0.0);
    let result = solver.solve(
        &x0,
        &vec![0.0; horizon],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    assert!(result.converged, "LQ problem must converge");
    assert!(
        result.iterations <= 3,
        "an LQ problem needs one step plus the stationarity check, took {}",
        result.iterations
    );
    let q = [10.0, 0.0, 0.0, 10.0];
    let lqr = FiniteHorizonLqr::new(2, 1, &a, &b, &q, &[0.1], &q, horizon);
    assert_eq!(result.gains.len(), horizon);
    for k in 0..horizon {
        for j in 0..2 {
            let (got, want) = (result.gains[k][j], lqr.gains[k][j]);
            assert!(
                (got - want).abs() < 1e-4 * want.abs().max(1.0),
                "gain K[{k}][{j}] = {got}, LQR says {want}"
            );
        }
    }
    // The open-loop optimum is the LQR closed-loop trajectory.
    let mut x = x0.to_vec();
    for k in 0..horizon {
        let u = lqr.control(k, &x);
        assert!(
            (u[0] - result.u_seq[k]).abs() < 1e-5,
            "u[{k}] = {}, LQR says {}",
            result.u_seq[k],
            u[0]
        );
        x = double_integrator(&x, &u);
    }
    // And the identity-Hessian shortcut of the previous implementation would
    // have produced different gains: Q = I/2, R = I/2 would give far smaller
    // position gains than Q = 10 I, R = 0.1.
    let weak = FiniteHorizonLqr::new(
        2,
        1,
        &a,
        &b,
        &[0.5, 0.0, 0.0, 0.5],
        &[0.5],
        &[1.0, 0.0, 0.0, 1.0],
        horizon,
    );
    assert!((weak.gains[0][0] - result.gains[0][0]).abs() > 0.5);
}

#[test]
fn cost_cross_terms_and_terminal_hessian_are_used() {
    // A positive definite quadratic cost with state-control cross terms and an
    // off-diagonal terminal Hessian.  J(u) is quadratic, so one dense Newton
    // step on the total cost is the exact minimiser.
    let horizon = 6;
    let x0 = [0.8, -0.4];
    let running = |x: &[f64], u: &[f64]| {
        x[0] * x[0] + 4.0 * x[1] * x[1] + 2.0 * x[0] * u[0] + 3.0 * u[0] * u[0] + 0.5 * x[0] * x[1]
    };
    let terminal = |x: &[f64]| 6.0 * x[0] * x[0] + 2.0 * x[1] * x[1] + x[0] * x[1];
    let cost = |u: &[f64]| total_cost(&x0, u, 1, double_integrator, running, terminal);
    let u_star = newton_step(&vec![0.0; horizon], &cost);
    let solver = IlqrSolver::new(2, 1, horizon, 20, 1e-12, 0.0);
    let result = solver.solve(
        &x0,
        &vec![0.0; horizon],
        double_integrator,
        running,
        terminal,
    );
    assert!(result.converged);
    for (k, (got, want)) in result.u_seq.iter().zip(u_star.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-5,
            "u[{k}] = {got}, dense Newton says {want}"
        );
    }
    assert!((result.cost - cost(&u_star)).abs() < 1e-8 * cost(&u_star).abs().max(1.0));
}

// ── DDP second-order terms ──────────────────────────────────────────────

fn scalar_dynamics(x: &[f64], u: &[f64]) -> Vec<f64> {
    vec![x[0] + 0.3 * u[0] + 0.4 * u[0] * u[0]]
}

#[test]
fn ddp_takes_the_exact_newton_step_where_ilqr_takes_the_gauss_newton_step() {
    // One stage: q(u) = 0.5 u² + 0.2 x0² + (f(x0, u) - 1)² with f = x0 + 0.3 u + 0.4 u².
    //   q'(u)  = u + 2 (f - 1) f_u,
    //   q''(u) = 1 + 2 f_u² + 2 (f - 1) f_uu        (DDP keeps the last term).
    let (x0, u0): (f64, f64) = (0.5, 0.2);
    let f = x0 + 0.3 * u0 + 0.4 * u0 * u0;
    let (fu, fuu) = (0.3 + 0.8 * u0, 0.8);
    let q1 = u0 + 2.0 * (f - 1.0) * fu;
    let newton = u0 - q1 / (1.0 + 2.0 * fu * fu + 2.0 * (f - 1.0) * fuu);
    let gauss_newton = u0 - q1 / (1.0 + 2.0 * fu * fu);
    assert!((newton - 0.4552094522019333).abs() < 1e-12);
    assert!((gauss_newton - 0.3335581787521079).abs() < 1e-12);

    let running = |x: &[f64], u: &[f64]| 0.5 * u[0] * u[0] + 0.2 * x[0] * x[0];
    let terminal = |x: &[f64]| (x[0] - 1.0) * (x[0] - 1.0);
    let ddp = DifferentialDynamicProgramming::new(1, 1, 1, 1, 0.0).solve(
        &[x0],
        &[u0],
        scalar_dynamics,
        running,
        terminal,
    );
    assert!(
        (ddp.u_seq[0] - newton).abs() < 1e-5,
        "DDP step {} should be the Newton step {newton}",
        ddp.u_seq[0]
    );
    let ilqr = IlqrSolver::new(1, 1, 1, 1, 0.0, 0.0).solve(
        &[x0],
        &[u0],
        scalar_dynamics,
        running,
        terminal,
    );
    assert!(
        (ilqr.u_seq[0] - gauss_newton).abs() < 1e-5,
        "iLQR step {} should be the Gauss-Newton step {gauss_newton}",
        ilqr.u_seq[0]
    );
    assert!((ddp.u_seq[0] - ilqr.u_seq[0]).abs() > 0.1);
}

/// Minimiser of `l(x1, u) + φ(f(x1, u))` over the scalar control by Newton's method.
fn stage_minimiser(
    x1: f64,
    f: &dyn Fn(f64, f64) -> f64,
    l: &dyn Fn(f64, f64) -> f64,
    phi: &dyn Fn(f64) -> f64,
) -> f64 {
    let q = |u: f64| l(x1, u) + phi(f(x1, u));
    let mut u = 0.0;
    for _ in 0..80 {
        let h = 1e-5;
        let g = (q(u + h) - q(u - h)) / (2.0 * h);
        let hess = (q(u + h) - 2.0 * q(u) + q(u - h)) / (h * h);
        let step = -g / hess;
        u += step;
        if step.abs() < 1e-15 {
            break;
        }
    }
    u
}

#[test]
fn ddp_backward_pass_matches_the_exact_nested_value_function() {
    // Two stages.  With the second control at its exact minimiser the DDP value
    // function V1(x) is the true cost-to-go to second order, so the stage-0
    // policy must equal the Newton quantities of Q0(x, u) = l + V1(f), where V1
    // is evaluated by an exact inner minimisation.
    let f = |x: f64, u: f64| x + 0.2 * x.sin() + 0.3 * u + 0.15 * u * u + 0.1 * x * u;
    let l = |x: f64, u: f64| 0.5 * u * u + 0.3 * x * x;
    let phi = |x: f64| 1.5 * (x - 0.5) * (x - 0.5);
    let (x0, u0) = (0.4, 0.1);
    let value = |x1: f64| {
        let u = stage_minimiser(x1, &f, &l, &phi);
        l(x1, u) + phi(f(x1, u))
    };
    let q0 = |z: &[f64]| l(z[0], z[1]) + value(f(z[0], z[1]));
    let (grad, hess) = gradient_and_hessian(&q0, &[x0, u0]);
    let (qu, quu, qux) = (grad[1], hess[3], hess[2]);
    let (expected_k, expected_gain) = (qu / quu, qux / quu);

    let x1 = f(x0, u0);
    let u1 = stage_minimiser(x1, &f, &l, &phi);
    let dynamics = |x: &[f64], u: &[f64]| vec![f(x[0], u[0])];
    let running = |x: &[f64], u: &[f64]| l(x[0], u[0]);
    let terminal = |x: &[f64]| phi(x[0]);
    let settings = |order| IlqrSettings {
        n: 1,
        m: 1,
        horizon: 2,
        max_iter: 1,
        tol: 0.0,
        reg: 0.0,
        order,
    };
    let ddp = policy_at(
        &settings(DynamicsOrder::Second),
        &[x0],
        &[u0, u1],
        dynamics,
        running,
        terminal,
        0.0,
    )
    .expect("positive definite");
    assert!(
        (ddp.k_ff[0][0] - expected_k).abs() < 1e-6,
        "k0 {} vs {expected_k}",
        ddp.k_ff[0][0]
    );
    assert!(
        (ddp.gains[0][0] - expected_gain).abs() < 1e-6,
        "K0 {} vs {expected_gain}",
        ddp.gains[0][0]
    );
    // The Gauss-Newton (iLQR) policy ignores the curvature of the dynamics and
    // misses the exact one by far more than the finite-difference noise.
    let gn = policy_at(
        &settings(DynamicsOrder::First),
        &[x0],
        &[u0, u1],
        dynamics,
        running,
        terminal,
        0.0,
    )
    .expect("positive definite");
    assert!((gn.k_ff[0][0] - expected_k).abs() > 1e-2);
    assert!((gn.gains[0][0] - expected_gain).abs() > 1e-2);
}

// ── Nonlinear problems ──────────────────────────────────────────────────

const PENDULUM_DT: f64 = 0.05;

fn pendulum(x: &[f64], u: &[f64]) -> Vec<f64> {
    vec![
        x[0] + PENDULUM_DT * x[1],
        x[1] + PENDULUM_DT * (-9.81 * x[0].sin() + u[0]),
    ]
}

fn pendulum_running(x: &[f64], u: &[f64]) -> f64 {
    x[0] * x[0] + 0.1 * x[1] * x[1] + 0.05 * u[0] * u[0]
}

fn pendulum_terminal(x: &[f64]) -> f64 {
    50.0 * (x[0] * x[0] + x[1] * x[1])
}

fn gradient_of_total_cost(x0: &[f64], u: &[f64]) -> Vec<f64> {
    let cost = |z: &[f64]| total_cost(x0, z, 1, pendulum, pendulum_running, pendulum_terminal);
    gradient_and_hessian(&cost, u).0
}

#[test]
fn nonlinear_solves_reach_a_stationary_point_of_the_total_cost() {
    let horizon = 12;
    let x0 = [1.2, 0.0];
    let u_init = vec![0.0; horizon];
    let initial_cost = total_cost(
        &x0,
        &u_init,
        1,
        pendulum,
        pendulum_running,
        pendulum_terminal,
    );
    let ilqr = IlqrSolver::new(2, 1, horizon, 50, 1e-10, 1e-6).solve(
        &x0,
        &u_init,
        pendulum,
        pendulum_running,
        pendulum_terminal,
    );
    let ddp = DifferentialDynamicProgramming::new(2, 1, horizon, 50, 1e-10).solve(
        &x0,
        &u_init,
        pendulum,
        pendulum_running,
        pendulum_terminal,
    );
    for (name, result) in [("iLQR", &ilqr), ("DDP", &ddp)] {
        assert!(result.converged, "{name} must converge");
        assert!(
            result.cost < 0.5 * initial_cost,
            "{name}: {} vs initial {initial_cost}",
            result.cost
        );
        let recomputed = total_cost(
            &x0,
            &result.u_seq,
            1,
            pendulum,
            pendulum_running,
            pendulum_terminal,
        );
        assert!(
            (recomputed - result.cost).abs() < 1e-9 * recomputed.abs(),
            "{name}: reported cost"
        );
        let grad = gradient_of_total_cost(&x0, &result.u_seq);
        let worst = grad.iter().fold(0.0_f64, |m, g| m.max(g.abs()));
        assert!(
            worst < 1e-4,
            "{name}: stationarity violated, max |dJ/du| = {worst}"
        );
        assert_eq!(result.x_traj.len(), (horizon + 1) * 2);
        assert_eq!(result.gains.len(), horizon);
        assert!(result.gains.iter().all(|k| k.len() == 2));
    }
    // Both methods find the same local minimum.
    assert!(
        (ilqr.cost - ddp.cost).abs() < 1e-7 * ilqr.cost,
        "{} vs {}",
        ilqr.cost,
        ddp.cost
    );
    // Second-order dynamics terms never cost extra iterations here.
    assert!(
        ddp.iterations <= ilqr.iterations,
        "DDP took {} iterations, iLQR {}",
        ddp.iterations,
        ilqr.iterations
    );
}

#[test]
fn closed_loop_gains_correct_state_deviations_along_the_solution() {
    let horizon = 12;
    let x0 = [1.2, 0.0];
    let result = DifferentialDynamicProgramming::new(2, 1, horizon, 50, 1e-10).solve(
        &x0,
        &vec![0.0; horizon],
        pendulum,
        pendulum_running,
        pendulum_terminal,
    );
    // Start the closed loop slightly off the nominal state: the terminal error
    // with feedback must be smaller than without.
    let off = [x0[0] + 0.05, x0[1] - 0.03];
    let (mut x_fb, mut x_ol) = (off.to_vec(), off.to_vec());
    for k in 0..horizon {
        let dx = [
            x_fb[0] - result.x_traj[2 * k],
            x_fb[1] - result.x_traj[2 * k + 1],
        ];
        let u_fb = result.u_seq[k] - (result.gains[k][0] * dx[0] + result.gains[k][1] * dx[1]);
        x_fb = pendulum(&x_fb, &[u_fb]);
        x_ol = pendulum(&x_ol, &[result.u_seq[k]]);
    }
    let nominal_end = [result.x_traj[2 * horizon], result.x_traj[2 * horizon + 1]];
    let err =
        |x: &[f64]| ((x[0] - nominal_end[0]).powi(2) + (x[1] - nominal_end[1]).powi(2)).sqrt();
    assert!(
        err(&x_fb) < 0.5 * err(&x_ol),
        "feedback error {} should be well below open-loop error {}",
        err(&x_fb),
        err(&x_ol)
    );
}

#[test]
fn indefinite_control_curvature_is_regularised() {
    // l = -u² + u⁴ is a double well: Quu = -2 + O(0.02) < 0 at u = 0, so the
    // plain Cholesky fails and the regularisation has to grow before the
    // iteration can leave the origin for one of the wells at |u| ≈ 0.69.
    let horizon = 8;
    let dynamics = |x: &[f64], u: &[f64]| vec![x[0] + 0.1 * u[0]];
    let running = |_x: &[f64], u: &[f64]| -u[0] * u[0] + u[0].powi(4);
    let terminal = |x: &[f64]| x[0] * x[0];
    let x0 = [0.3];
    let u_init = vec![0.0; horizon];
    let initial = total_cost(&x0, &u_init, 1, dynamics, running, terminal);
    for result in [
        IlqrSolver::new(1, 1, horizon, 60, 1e-10, 1e-6)
            .solve(&x0, &u_init, dynamics, running, terminal),
        DifferentialDynamicProgramming::new(1, 1, horizon, 60, 1e-10)
            .solve(&x0, &u_init, dynamics, running, terminal),
    ] {
        assert!(result.converged);
        assert!((initial - 0.09).abs() < 1e-12);
        assert!(
            (result.cost - -1.9321750718767874).abs() < 1e-6,
            "cost {}",
            result.cost
        );
        for u in &result.u_seq {
            assert!((u + 0.69398).abs() < 1e-4, "control {u}");
        }
    }
}

// ── Degenerate inputs ───────────────────────────────────────────────────

#[test]
fn empty_horizon_or_inputs_return_the_rollout() {
    let solver = IlqrSolver::new(2, 1, 0, 10, 1e-9, 1e-6);
    let result = solver.solve(
        &[1.0, 2.0],
        &[],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    assert_eq!(result.x_traj, vec![1.0, 2.0]);
    assert!(result.u_seq.is_empty() && result.gains.is_empty());
    assert_eq!(result.iterations, 0);
    assert!(result.converged);
    assert!((result.cost - 50.0).abs() < 1e-12);

    // No inputs: the trajectory is just the autonomous rollout.
    let free = |x: &[f64], _u: &[f64]| vec![0.5 * x[0]];
    let solver = IlqrSolver::new(1, 0, 3, 10, 1e-9, 1e-6);
    let result = solver.solve(&[8.0], &[], free, |x, _| x[0] * x[0], |x| x[0] * x[0]);
    assert_eq!(result.x_traj, vec![8.0, 4.0, 2.0, 1.0]);
    assert!(result.converged);
    assert!((result.cost - (64.0 + 16.0 + 4.0 + 1.0)).abs() < 1e-12);
}

#[test]
fn non_finite_dynamics_or_wrong_sized_output_do_not_panic() {
    let solver = IlqrSolver::new(2, 1, 4, 10, 1e-9, 1e-6);
    let nan = solver.solve(
        &[1.0, 0.0],
        &[0.0; 4],
        |_x: &[f64], _u: &[f64]| vec![f64::NAN, 0.0],
        lqr_running,
        lqr_terminal,
    );
    assert!(!nan.converged && !nan.cost.is_finite());
    assert_eq!(nan.iterations, 0);
    // A dynamics function returning too few components is padded with NaN.
    let short = solver.solve(
        &[1.0, 0.0],
        &[0.0; 4],
        |_x: &[f64], _u: &[f64]| vec![0.0],
        lqr_running,
        lqr_terminal,
    );
    assert!(!short.converged);
    // Short initial controls are zero padded, surplus ones ignored.
    let padded = solver.solve(
        &[1.0, 0.0],
        &[0.0],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    let full = solver.solve(
        &[1.0, 0.0],
        &[0.0; 4],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    assert_eq!(padded.u_seq.len(), 4);
    assert!((padded.cost - full.cost).abs() < 1e-9);
    let surplus = solver.solve(
        &[1.0, 0.0],
        &[0.0; 9],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    assert_eq!(surplus.u_seq.len(), 4);
    // max_iter = 0 returns the initial rollout untouched.
    let idle = IlqrSolver::new(2, 1, 4, 0, 1e-9, 1e-6).solve(
        &[1.0, 0.0],
        &[0.0; 4],
        double_integrator,
        lqr_running,
        lqr_terminal,
    );
    assert_eq!(idle.iterations, 0);
    assert!(!idle.converged);
    assert_eq!(idle.u_seq, vec![0.0; 4]);
}
