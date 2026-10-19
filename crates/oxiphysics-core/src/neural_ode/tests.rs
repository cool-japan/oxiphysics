//! Regression tests for the Neural ODE building blocks.

use super::*;

// ── rk4_step ─────────────────────────────────────────────────────────────

#[test]
fn test_rk4_exponential_decay() {
    // dy/dt = -y, y(0)=1 → y(t) = exp(-t)
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let y0 = vec![1.0];
    let y1 = rk4_step(&f, 0.0, &y0, 0.1);
    let exact = (-0.1_f64).exp();
    assert!(
        (y1[0] - exact).abs() < 1e-6,
        "RK4 decay: got {}, expected {}",
        y1[0],
        exact
    );
}

#[test]
fn test_rk4_harmonic_oscillator() {
    // d²x/dt² = -x → state [x, v], dz/dt = [v, -x]
    let f = |_t: f64, z: &[f64]| vec![z[1], -z[0]];
    let z0 = vec![1.0, 0.0]; // x=1, v=0 → x(t)=cos(t)
    let mut z = z0.clone();
    let dt = 0.01;
    let steps = 100; // advance to t=1.0
    for i in 0..steps {
        z = rk4_step(&f, i as f64 * dt, &z, dt);
    }
    let t = 1.0_f64;
    let exact_x = t.cos();
    assert!(
        (z[0] - exact_x).abs() < 1e-5,
        "Harmonic oscillator x: got {}",
        z[0]
    );
}

#[test]
fn test_rk4_constant_ode() {
    // dy/dt = 2, y(0)=0 → y(1)=2
    let f = |_t: f64, _y: &[f64]| vec![2.0];
    let y = rk4_step(&f, 0.0, &[0.0], 1.0);
    assert!((y[0] - 2.0).abs() < 1e-12);
}

#[test]
fn test_rk4_zero_step() {
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let y0 = vec![3.0];
    let y1 = rk4_step(&f, 0.0, &y0, 0.0);
    assert!((y1[0] - 3.0).abs() < 1e-15);
}

#[test]
fn test_rk4_linear_ode() {
    // dy/dt = t, y(0)=0 → y(2)=2
    let f = |t: f64, _y: &[f64]| vec![t];
    let mut y = vec![0.0];
    let dt = 0.01;
    for i in 0..200 {
        y = rk4_step(&f, i as f64 * dt, &y, dt);
    }
    assert!((y[0] - 2.0).abs() < 1e-6, "Linear ODE: got {}", y[0]);
}

#[test]
fn test_rk4_2d_decoupled() {
    // [dy1/dt, dy2/dt] = [-y1, -2*y2], y(0)=[1, 1]
    let f = |_t: f64, y: &[f64]| vec![-y[0], -2.0 * y[1]];
    let mut z = vec![1.0_f64, 1.0_f64];
    let dt = 0.01;
    for i in 0..50 {
        z = rk4_step(&f, i as f64 * dt, &z, dt);
    }
    let t = 0.5_f64;
    assert!((z[0] - (-t).exp()).abs() < 1e-5, "y1: {}", z[0]);
    assert!((z[1] - (-2.0 * t).exp()).abs() < 1e-5, "y2: {}", z[1]);
}

// ── dopri5_step ───────────────────────────────────────────────────────────

#[test]
fn test_dopri5_returns_three_values() {
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let (yh, yl, err) = dopri5_step(&f, 0.0, &[1.0], 0.1, 1e-3, 1e-6);
    assert_eq!(yh.len(), 1);
    assert_eq!(yl.len(), 1);
    assert!(err.is_finite());
}

#[test]
fn test_dopri5_exponential_accuracy() {
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let (yh, _yl, _err) = dopri5_step(&f, 0.0, &[1.0], 0.1, 1e-6, 1e-9);
    let exact = (-0.1_f64).exp();
    assert!(
        (yh[0] - exact).abs() < 1e-8,
        "DOPRI5 accuracy: {}",
        (yh[0] - exact).abs()
    );
}

#[test]
fn test_dopri5_zero_step_size() {
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let (yh, yl, err) = dopri5_step(&f, 0.0, &[1.0], 0.0, 1e-3, 1e-6);
    assert!((yh[0] - 1.0).abs() < 1e-12);
    assert!((yl[0] - 1.0).abs() < 1e-12);
    assert!(err < 1e-10);
}

// ── NeuralOdeFunc ─────────────────────────────────────────────────────────

#[test]
fn test_neural_ode_func_forward_shape() {
    let func = NeuralOdeFunc::new(3, 8, 42);
    let z = vec![1.0, 0.0, -1.0];
    let dz = func.forward(0.0, &z);
    assert_eq!(dz.len(), 3);
}

#[test]
fn test_neural_ode_func_forward_finite() {
    let func = NeuralOdeFunc::new(4, 16, 1234);
    let z = vec![0.5, -0.3, 1.2, -0.1];
    let dz = func.forward(1.0, &z);
    for &v in &dz {
        assert!(
            v.is_finite(),
            "NeuralOdeFunc output contains non-finite: {v}"
        );
    }
}

#[test]
fn test_neural_ode_func_deterministic() {
    let f1 = NeuralOdeFunc::new(2, 4, 99);
    let f2 = NeuralOdeFunc::new(2, 4, 99);
    let z = vec![0.1, 0.2];
    assert_eq!(f1.forward(0.0, &z), f2.forward(0.0, &z));
}

#[test]
fn test_neural_ode_func_different_seeds_differ() {
    let f1 = NeuralOdeFunc::new(2, 8, 1);
    let f2 = NeuralOdeFunc::new(2, 8, 2);
    let z = vec![1.0, 1.0];
    let d1 = f1.forward(0.0, &z);
    let d2 = f2.forward(0.0, &z);
    let diff: f64 = d1.iter().zip(d2.iter()).map(|(a, b)| (a - b).abs()).sum();
    assert!(
        diff > 1e-10,
        "Different seeds should give different outputs"
    );
}

#[test]
fn test_neural_ode_func_jvp_shape() {
    let func = NeuralOdeFunc::new(3, 6, 7);
    let z = vec![0.0, 1.0, -1.0];
    let v = vec![1.0, 0.0, 0.0];
    let jvp = func.jvp(0.5, &z, &v, 1e-5);
    assert_eq!(jvp.len(), 3);
}

// ── NeuralOdeSolver ───────────────────────────────────────────────────────

#[test]
fn test_solver_rk4_output_shape() {
    let func = NeuralOdeFunc::new(2, 4, 0);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let z0 = vec![1.0, 0.0];
    let z1 = solver.solve_rk4(&z0, 0.0, 1.0, 0.1);
    assert_eq!(z1.len(), 2);
}

#[test]
fn test_solver_rk4_zero_integration() {
    // When t0 == t1 the state should be unchanged
    let func = NeuralOdeFunc::new(2, 4, 5);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let z0 = vec![1.0, 2.0];
    let z1 = solver.solve_rk4(&z0, 0.0, 0.0, 0.1);
    // With no steps the state equals z0
    for (a, b) in z0.iter().zip(z1.iter()) {
        assert!((a - b).abs() < 1e-12);
    }
}

#[test]
fn test_solver_rk4_finite_output() {
    let func = NeuralOdeFunc::new(3, 8, 100);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let z0 = vec![0.1, -0.2, 0.3];
    let z1 = solver.solve_rk4(&z0, 0.0, 0.5, 0.05);
    for &v in &z1 {
        assert!(v.is_finite());
    }
}

#[test]
fn test_solver_trajectory_length() {
    let func = NeuralOdeFunc::new(2, 4, 3);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let z0 = vec![1.0, 0.0];
    let ts = vec![0.0, 0.25, 0.5, 0.75, 1.0];
    let traj = solver.solve_rk4_trajectory(&z0, &ts, 0.05);
    assert_eq!(traj.len(), ts.len());
}

#[test]
fn test_solver_dopri5_output_shape() {
    let func = NeuralOdeFunc::new(2, 4, 42);
    let solver = NeuralOdeSolver::new(func, 1e-4, 1e-7);
    let z0 = vec![1.0, 0.5];
    let z1 = solver.solve_dopri5(&z0, 0.0, 1.0, 0.1);
    assert_eq!(z1.len(), 2);
}

#[test]
fn test_solver_dopri5_finite_output() {
    let func = NeuralOdeFunc::new(3, 6, 77);
    let solver = NeuralOdeSolver::new(func, 1e-4, 1e-7);
    let z0 = vec![0.0, 0.5, 1.0];
    let z1 = solver.solve_dopri5(&z0, 0.0, 0.5, 0.1);
    for &v in &z1 {
        assert!(v.is_finite(), "DOPRI5 produced non-finite: {v}");
    }
}

// ── AdjointMethod ─────────────────────────────────────────────────────────

#[test]
fn test_adjoint_backward_shape() {
    let func = NeuralOdeFunc::new(4, 8, 33);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let adj = AdjointMethod::new(4);
    let z0 = vec![0.2, -0.1, 0.4, 0.0];
    let loss_grad = vec![1.0, -1.0, 0.5, -0.5];
    let grad = adj.backward(&solver, &z0, &loss_grad, 0.0, 0.1);
    assert_eq!(grad.len(), 4);
}

/// The one-step adjoint VJP must equal the gradient of the *actual* loss
/// `L(z0) = loss_grad · Φ(z0)` (Φ = one RK4 step of the real dynamics) with
/// respect to z0 — verified against an independent central finite difference
/// of L itself.  This pins down that `backward` is a real adjoint, NOT the
/// old `-loss_grad` fabrication.
#[test]
fn test_adjoint_backward_matches_finite_difference_of_loss() {
    let func = NeuralOdeFunc::new(3, 6, 101);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let adj = AdjointMethod::new(3);
    let z0 = vec![0.3, -0.2, 0.5];
    let loss_grad = vec![0.7, -1.3, 0.4];
    let t0 = 0.0;
    let dt = 0.15;

    let grad = adj.backward(&solver, &z0, &loss_grad, t0, dt);

    // Independent reference: L(z0) = loss_grad · Φ(z0); dL/dz0[j] via FD.
    let loss = |z: &[f64]| -> f64 {
        let forward = |t: f64, y: &[f64]| solver.func.forward(t, y);
        let z1 = rk4_step(&forward, t0, z, dt);
        loss_grad.iter().zip(z1.iter()).map(|(g, z)| g * z).sum()
    };
    let fd_eps = 1e-6;
    for j in 0..3 {
        let mut zp = z0.clone();
        let mut zm = z0.clone();
        zp[j] += fd_eps;
        zm[j] -= fd_eps;
        let fd = (loss(&zp) - loss(&zm)) / (2.0 * fd_eps);
        assert!(
            (grad[j] - fd).abs() < 1e-4,
            "adjoint grad[{j}]={} disagrees with FD dL/dz0={fd}",
            grad[j]
        );
    }

    // And it is emphatically NOT the old fabrication (-loss_grad).
    let negation: Vec<f64> = loss_grad.iter().map(|g| -g).collect();
    let differs = grad
        .iter()
        .zip(negation.iter())
        .any(|(a, b)| (a - b).abs() > 1e-3);
    assert!(differs, "backward must not return the negated loss_grad");
}

#[test]
fn test_adjoint_run_returns_correct_shapes() {
    let func = NeuralOdeFunc::new(2, 4, 11);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let mut adj = AdjointMethod::new(2);
    let z_final = vec![0.5, -0.5];
    let loss_grad = vec![1.0, 0.0];
    let (grad_z0, grad_params) = adj.run(&solver, &z_final, &loss_grad, 0.0, 1.0, 0.1);
    assert_eq!(grad_z0.len(), 2);
    assert!(!grad_params.is_empty());
}

#[test]
fn test_adjoint_run_finite() {
    let func = NeuralOdeFunc::new(2, 4, 22);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let mut adj = AdjointMethod::new(2);
    let z_final = vec![1.0, 1.0];
    let loss_grad = vec![0.1, -0.1];
    let (g, _) = adj.run(&solver, &z_final, &loss_grad, 0.0, 1.0, 0.1);
    for &v in &g {
        assert!(v.is_finite());
    }
}

// ── LatentOde ─────────────────────────────────────────────────────────────

#[test]
fn test_latent_ode_encode_shape() {
    let model = LatentOde::new(4, 2, 8, 55);
    let obs = vec![vec![1.0, 0.0, -1.0, 0.5], vec![0.5, 0.1, -0.5, 0.3]];
    let z = model.encode(&obs);
    assert_eq!(z.len(), 2);
}

#[test]
fn test_latent_ode_encode_empty() {
    let model = LatentOde::new(3, 2, 4, 1);
    let z = model.encode(&[]);
    assert_eq!(z.len(), 2);
    assert!(z.iter().all(|&v| v == 0.0));
}

#[test]
fn test_latent_ode_decode_single_shape() {
    let model = LatentOde::new(4, 2, 6, 88);
    let z = vec![0.5, -0.3];
    let obs = model.decode_single(&z);
    assert_eq!(obs.len(), 4);
}

#[test]
fn test_latent_ode_decode_trajectory_length() {
    let model = LatentOde::new(3, 2, 4, 33);
    let z = vec![0.1, 0.2];
    let ts = vec![0.1, 0.2, 0.5, 1.0];
    let preds = model.decode(&z, 0.0, &ts, 0.05);
    // times prepended with t0 → trajectory has len(ts)+1, then decoded: len = ts.len()+1
    assert_eq!(preds.len(), ts.len() + 1);
}

#[test]
fn test_latent_ode_encode_finite() {
    let model = LatentOde::new(3, 4, 8, 999);
    let obs: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64 * 0.1; 3]).collect();
    let z = model.encode(&obs);
    assert!(
        z.iter().all(|v| v.is_finite()),
        "Encoded latent contains non-finite"
    );
}

#[test]
fn test_latent_ode_round_trip_shape() {
    let model = LatentOde::new(2, 2, 4, 77);
    let obs = vec![vec![1.0, 0.0], vec![0.8, 0.1]];
    let z = model.encode(&obs);
    let recon = model.decode_single(&z);
    assert_eq!(recon.len(), 2);
}

// ── TimeSeriesOde ─────────────────────────────────────────────────────────

#[test]
fn test_time_series_ode_predict_shape() {
    let func = NeuralOdeFunc::new(2, 4, 13);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let times = vec![0.0, 0.5, 1.0];
    let obs = vec![vec![1.0, 0.0], vec![0.9, 0.1], vec![0.8, 0.2]];
    let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
    let pred = ts.predict(1.5);
    assert_eq!(pred.len(), 2);
}

#[test]
fn test_time_series_ode_loss_nonnegative() {
    let func = NeuralOdeFunc::new(2, 4, 14);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let times = vec![0.0, 0.5, 1.0];
    let obs = vec![vec![1.0, 0.0], vec![0.9, 0.1], vec![0.8, 0.2]];
    let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
    assert!(ts.compute_loss(0.05) >= 0.0);
}

#[test]
fn test_time_series_ode_fit_records_loss() {
    let func = NeuralOdeFunc::new(1, 4, 15);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let times = vec![0.0, 0.1, 0.2, 0.3];
    let obs: Vec<Vec<f64>> = (0..4).map(|i| vec![(-(i as f64) * 0.1).exp()]).collect();
    let mut ts = TimeSeriesOde::new(times, obs, solver, 0.001, 5);
    ts.fit();
    assert_eq!(ts.loss_history.len(), 5);
}

#[test]
fn test_time_series_ode_predict_finite() {
    let func = NeuralOdeFunc::new(2, 4, 16);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let times = vec![0.0, 0.5];
    let obs = vec![vec![1.0, 0.0], vec![0.9, -0.1]];
    let ts = TimeSeriesOde::new(times, obs, solver, 0.01, 0);
    let pred = ts.predict(0.3);
    assert!(pred.iter().all(|v| v.is_finite()));
}

#[test]
fn test_time_series_ode_empty() {
    let func = NeuralOdeFunc::new(2, 4, 17);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let ts = TimeSeriesOde::new(vec![], vec![], solver, 0.01, 0);
    let pred = ts.predict(1.0);
    assert!(pred.is_empty());
    assert_eq!(ts.compute_loss(0.1), 0.0);
}

// ── Integration accuracy tests ────────────────────────────────────────────

#[test]
fn test_rk4_logistic_growth() {
    // dy/dt = y(1-y), y(0)=0.1 → y(t) = 1/(1 + 9*exp(-t))
    let f = |_t: f64, y: &[f64]| vec![y[0] * (1.0 - y[0])];
    let mut y = vec![0.1];
    let dt = 0.01;
    let steps = 200;
    for i in 0..steps {
        y = rk4_step(&f, i as f64 * dt, &y, dt);
    }
    let t = 2.0_f64;
    let exact = 1.0 / (1.0 + 9.0 * (-t).exp());
    assert!(
        (y[0] - exact).abs() < 1e-5,
        "Logistic growth: got {}, expected {}",
        y[0],
        exact
    );
}

#[test]
fn test_rk4_accuracy_order() {
    // Compare RK4 error at h=0.1 vs h=0.05 on exponential decay
    // RK4 is 4th-order: error ~ h^4, so halving h reduces error by ~16x
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let exact = (-1.0_f64).exp();

    let y_h1 = {
        let mut y = vec![1.0];
        for i in 0..10 {
            y = rk4_step(&f, i as f64 * 0.1, &y, 0.1);
        }
        y[0]
    };
    let y_h2 = {
        let mut y = vec![1.0];
        for i in 0..20 {
            y = rk4_step(&f, i as f64 * 0.05, &y, 0.05);
        }
        y[0]
    };
    let err1 = (y_h1 - exact).abs();
    let err2 = (y_h2 - exact).abs();
    assert!(
        err2 < err1,
        "Smaller step should give smaller error: {} vs {}",
        err2,
        err1
    );
}

#[test]
fn test_rk4_system_energy_conservation() {
    // Harmonic oscillator: H = 0.5*(x^2 + v^2) = 0.5 (for x0=1, v0=0)
    // should be approximately conserved by RK4
    let f = |_t: f64, z: &[f64]| vec![z[1], -z[0]];
    let mut z = vec![1.0, 0.0];
    let dt = 0.001;
    let steps = 1000;
    for i in 0..steps {
        z = rk4_step(&f, i as f64 * dt, &z, dt);
    }
    let energy = 0.5 * (z[0].powi(2) + z[1].powi(2));
    assert!(
        (energy - 0.5).abs() < 1e-4,
        "Energy drift: {}",
        energy - 0.5
    );
}

#[test]
fn test_neural_ode_func_batch_consistency() {
    // forward(t, z) should give the same result when called twice
    let func = NeuralOdeFunc::new(3, 8, 42);
    let z = vec![0.1, -0.2, 0.3];
    let d1 = func.forward(0.5, &z);
    let d2 = func.forward(0.5, &z);
    assert_eq!(d1, d2, "forward must be deterministic");
}

#[test]
fn test_time_series_ode_fit_loss_finite() {
    let func = NeuralOdeFunc::new(1, 4, 18);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let times: Vec<f64> = (0..5).map(|i| i as f64 * 0.2).collect();
    let obs: Vec<Vec<f64>> = times.iter().map(|&t: &f64| vec![(-t).exp()]).collect();
    let mut ts = TimeSeriesOde::new(times, obs, solver, 0.001, 3);
    ts.fit();
    for &l in &ts.loss_history {
        assert!(l.is_finite(), "Loss is non-finite: {l}");
    }
}

#[test]
fn test_rk4_step_multidim() {
    // 5-dimensional decay: dy_i/dt = -i*y_i, y_i(0)=1
    let f = |_t: f64, y: &[f64]| (0..y.len()).map(|i| -(i as f64 + 1.0) * y[i]).collect();
    let y0: Vec<f64> = vec![1.0; 5];
    let mut y = y0.clone();
    let dt = 0.01;
    for k in 0..10 {
        y = rk4_step(&f, k as f64 * dt, &y, dt);
    }
    for (i, &yi) in y.iter().enumerate() {
        let exact = (-(i as f64 + 1.0) * 0.1).exp();
        assert!(
            (yi - exact).abs() < 1e-5,
            "dim {i}: got {yi}, expected {exact}"
        );
    }
}

// ── C1: DOPRI5 error-order verification ───────────────────────────────────

#[test]
fn test_dopri5_error_estimate_order() {
    // For y' = y, y(0) = 1, exact = exp(t).
    // DOPRI5 is 5th-order: halving h should reduce |y_high - exact| by ~32×.
    let f = |_t: f64, y: &[f64]| vec![y[0]];
    let rtol = 1e-12;
    let atol = 1e-12;
    let y0 = vec![1.0_f64];

    let (y_big, _, _) = dopri5_step(&f, 0.0, &y0, 0.2, rtol, atol);
    let (y_small, _, _) = dopri5_step(&f, 0.0, &y0, 0.1, rtol, atol);
    let err_big = (y_big[0] - 0.2_f64.exp()).abs();
    let err_small = (y_small[0] - 0.1_f64.exp()).abs();
    // ratio ≈ (0.2/0.1)^5 = 32; require > 10 to avoid false negatives
    let ratio = err_big / err_small.max(f64::MIN_POSITIVE);
    assert!(
        ratio > 10.0,
        "Expected ~32× error reduction when halving step; got ratio={ratio:.2}"
    );
}

#[test]
fn test_dopri5_error_norm_small_step() {
    // error_norm should be < 1 for a well-behaved problem at h=0.01
    let f = |_t: f64, y: &[f64]| vec![-y[0]];
    let (_, _, err) = dopri5_step(&f, 0.0, &[1.0], 0.01, 1e-6, 1e-8);
    assert!(err < 1.0, "error norm should be < 1 for h=0.01: {err}");
}

// ── C2: BPTT gradient parity test ─────────────────────────────────────────

#[test]
fn test_bptt_gradient_nonzero_and_finite() {
    // Verify that BPTT parameter gradients are non-zero and finite.
    let func = NeuralOdeFunc::new(2, 4, 99);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let mut adj = AdjointMethod::new(2);
    let z_final = vec![0.5, -0.3];
    let loss_grad = vec![1.0, 0.0];
    let (_, grad_params) = adj.run(&solver, &z_final, &loss_grad, 0.0, 0.5, 0.1);
    assert_eq!(grad_params.len(), solver.func.n_params());
    assert!(
        grad_params.iter().all(|v| v.is_finite()),
        "some parameter gradients are non-finite"
    );
    assert!(
        grad_params.iter().any(|v| v.abs() > 1e-15),
        "all parameter gradients are zero"
    );
}

// ── Exact derivatives, continuous adjoint and time-series fitting ─────────

/// Network whose weights are scaled up so that the tanh units are genuinely
/// nonlinear (the default `±0.1` initialisation is almost linear).
fn scaled_func(input: usize, hidden: usize, seed: u64, factor: f64) -> NeuralOdeFunc {
    let mut func = NeuralOdeFunc::new(input, hidden, seed);
    let params: Vec<f64> = func.params_flat().iter().map(|p| p * factor).collect();
    func.set_params_flat(&params);
    func
}

/// One-dimensional network with a strong explicit dependence on time
/// (time-input weights of magnitude 2 to 3).
fn time_driven_func() -> NeuralOdeFunc {
    let mut func = NeuralOdeFunc::new(1, 3, 5);
    let params = [
        // weights_in (3 x 2): [z weight, time weight] per hidden unit
        0.4, 2.5, -0.3, 3.0, 0.2, -2.0, //
        // bias_in
        0.1, -0.2, 0.3, //
        // weights_hidden (3 x 3)
        0.8, -0.4, 0.5, 0.3, 0.9, -0.6, -0.7, 0.2, 0.6, //
        // bias_hidden
        0.05, -0.1, 0.15, //
        // weights_out (1 x 3)
        1.2, -0.8, 0.9, //
        // bias_out
        0.1,
    ];
    assert_eq!(params.len(), func.n_params());
    func.set_params_flat(&params);
    func
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Central finite difference of the scalar function `f` at every coordinate of `x`.
fn central_difference(f: &dyn Fn(&[f64]) -> f64, x: &[f64], eps: f64) -> Vec<f64> {
    (0..x.len())
        .map(|j| {
            let mut plus = x.to_vec();
            let mut minus = x.to_vec();
            plus[j] += eps;
            minus[j] -= eps;
            (f(&plus) - f(&minus)) / (2.0 * eps)
        })
        .collect()
}

#[test]
fn test_backprop_matches_finite_differences() {
    let func = scaled_func(3, 5, 7, 6.0);
    let z = [0.3, -0.7, 1.1];
    let adj = [0.5, -1.2, 0.8];
    let (vjp_z, grad_theta) = func.backprop(0.4, &z, &adj);
    assert_eq!(vjp_z.len(), 3);
    assert_eq!(grad_theta.len(), func.n_params());
    // Independent references: central differences of the scalar <adj, f>.
    let scalar_in_z = |zz: &[f64]| -> f64 {
        func.forward(0.4, zz)
            .iter()
            .zip(adj.iter())
            .map(|(f, a)| f * a)
            .sum()
    };
    let reference_z = central_difference(&scalar_in_z, &z, 1e-6);
    assert!(max_abs_diff(&vjp_z, &reference_z) < 1e-8);
    let scalar_in_theta = |theta: &[f64]| -> f64 {
        let mut shifted = func.clone();
        shifted.set_params_flat(theta);
        shifted
            .forward(0.4, &z)
            .iter()
            .zip(adj.iter())
            .map(|(f, a)| f * a)
            .sum()
    };
    let reference_theta = central_difference(&scalar_in_theta, &func.params_flat(), 1e-6);
    assert!(max_abs_diff(&grad_theta, &reference_theta) < 1e-8);
    // The public finite-difference helpers agree with the exact derivatives.
    assert!(max_abs_diff(&vjp_z, &func.vjp(0.4, &z, &adj, 1e-6)) < 1e-8);
    assert!(max_abs_diff(&grad_theta, &func.param_grad_contrib(0.4, &z, &adj, 1e-6)) < 1e-8);
}

#[test]
fn test_vjp_is_the_transpose_of_the_jvp() {
    // <u, J v> = <J^T u, v>.  The old adjoint integrated J·a instead of J^T·a;
    // for a non-symmetric Jacobian the two differ.
    let func = scaled_func(3, 6, 21, 5.0);
    let z = [0.2, -0.4, 0.9];
    let u = [1.0, -0.5, 0.25];
    let v = [-0.3, 0.8, 0.6];
    let jv = func.jvp(0.3, &z, &v, 1e-7);
    let (jtu, _) = func.backprop(0.3, &z, &u);
    let lhs: f64 = u.iter().zip(jv.iter()).map(|(a, b)| a * b).sum();
    let rhs: f64 = jtu.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
    assert!(
        (lhs - rhs).abs() < 1e-6,
        "<u,Jv> = {lhs}, <J^T u,v> = {rhs}"
    );
    // J^T u differs from J u for this network, so the distinction is real.
    let ju = func.jvp(0.3, &z, &u, 1e-7);
    assert!(max_abs_diff(&ju, &jtu) > 1e-3);
}

#[test]
fn test_rk4_step_vjp_matches_finite_differences() {
    let solver = NeuralOdeSolver::new(scaled_func(3, 4, 5, 4.0), 1e-3, 1e-6);
    let (t, h) = (0.3, 0.1);
    let y = [0.2, -0.5, 0.9];
    let cotangent = [1.0, -0.5, 0.25];
    let (y_bar, theta_bar) = solver.rk4_step_vjp(t, &y, h, &cotangent);
    let step_loss = |sv: &NeuralOdeSolver, state: &[f64]| -> f64 {
        let forward = |tt: f64, yy: &[f64]| sv.func.forward(tt, yy);
        rk4_step(&forward, t, state, h)
            .iter()
            .zip(cotangent.iter())
            .map(|(a, b)| a * b)
            .sum()
    };
    let reference_y = central_difference(&|s: &[f64]| step_loss(&solver, s), &y, 1e-6);
    assert!(max_abs_diff(&y_bar, &reference_y) < 1e-8);
    let reference_theta = central_difference(
        &|theta: &[f64]| {
            let mut shifted = solver.clone();
            shifted.func.set_params_flat(theta);
            step_loss(&shifted, &y)
        },
        &solver.func.params_flat(),
        1e-6,
    );
    assert!(max_abs_diff(&theta_bar, &reference_theta) < 1e-8);
}

/// Loss `L = g · z(t1)` of the discretised forward solve, for gradient checks.
fn solve_loss(solver: &NeuralOdeSolver, z0: &[f64], g: &[f64], t1: f64, dt: f64) -> f64 {
    solver
        .solve_rk4(z0, 0.0, t1, dt)
        .iter()
        .zip(g.iter())
        .map(|(a, b)| a * b)
        .sum()
}

#[test]
fn test_adjoint_run_matches_finite_differences_of_the_solve() {
    let solver = NeuralOdeSolver::new(scaled_func(2, 6, 11, 5.0), 1e-3, 1e-6);
    let z0 = vec![0.4, -0.3];
    let g = vec![0.7, -1.1];
    let (t1, dt) = (1.0, 0.02);
    let z_final = solver.solve_rk4(&z0, 0.0, t1, dt);
    let mut adjoint = AdjointMethod::new(2);
    let (grad_z0, grad_params) = adjoint.run(&solver, &z_final, &g, 0.0, t1, dt);
    let reference_z0 =
        central_difference(&|z: &[f64]| solve_loss(&solver, z, &g, t1, dt), &z0, 1e-6);
    assert!(
        max_abs_diff(&grad_z0, &reference_z0) < 1e-6,
        "grad_z0 {grad_z0:?} vs finite differences {reference_z0:?}"
    );
    let reference_params = central_difference(
        &|theta: &[f64]| {
            let mut shifted = solver.clone();
            shifted.func.set_params_flat(theta);
            solve_loss(&shifted, &z0, &g, t1, dt)
        },
        &solver.func.params_flat(),
        1e-6,
    );
    assert_eq!(grad_params.len(), reference_params.len());
    assert!(
        max_abs_diff(&grad_params, &reference_params) < 1e-6,
        "parameter gradient deviates by {}",
        max_abs_diff(&grad_params, &reference_params)
    );
    // The gradient is not trivial.
    assert!(grad_params.iter().any(|v| v.abs() > 1e-2));
}

#[test]
fn test_adjoint_run_samples_a_time_dependent_field_at_the_right_times() {
    // Time-input weights of magnitude 2-3: evaluating the stages at forward
    // times (t1 + h/2, t1 + h ...) instead of reverse-time stage times moves
    // the parameter gradient by ~2e-3 and grad_z0 by ~4e-5 at this step size.
    let solver = NeuralOdeSolver::new(time_driven_func(), 1e-3, 1e-6);
    let (z0, g) = (vec![0.5], vec![1.0]);
    let (t1, dt) = (2.0, 0.05);
    let z_final = solver.solve_rk4(&z0, 0.0, t1, dt);
    let mut adjoint = AdjointMethod::new(1);
    let (grad_z0, grad_params) = adjoint.run(&solver, &z_final, &g, 0.0, t1, dt);
    let reference_z0 =
        central_difference(&|z: &[f64]| solve_loss(&solver, z, &g, t1, dt), &z0, 1e-6);
    let reference_params = central_difference(
        &|theta: &[f64]| {
            let mut shifted = solver.clone();
            shifted.func.set_params_flat(theta);
            solve_loss(&shifted, &z0, &g, t1, dt)
        },
        &solver.func.params_flat(),
        1e-6,
    );
    assert!(
        max_abs_diff(&grad_z0, &reference_z0) < 1e-6,
        "grad_z0 {grad_z0:?} vs {reference_z0:?}"
    );
    assert!(
        max_abs_diff(&grad_params, &reference_params) < 5e-6,
        "parameter gradient deviates by {}",
        max_abs_diff(&grad_params, &reference_params)
    );
}

#[test]
fn test_adjoint_run_linear_regime_matches_the_closed_form() {
    // f(z) = c tanh(b tanh(a z)) ~ (a b c) z for tiny a z, so
    // z(T) = z0 e^{lambda T}; with L = z(T): dL/dz0 = e^{lambda T} and
    // dL/dc = z0 a b T e^{lambda T}.  The old code produced e^{-lambda T}.
    let (a_w, b_w, c_w) = (1e-3, 1.0, 800.0);
    let lambda = a_w * b_w * c_w;
    let mut func = NeuralOdeFunc::new(1, 1, 1);
    // weights_in = [a, 0 (time)], bias_in, weights_hidden, bias_hidden, weights_out, bias_out
    func.set_params_flat(&[a_w, 0.0, 0.0, b_w, 0.0, c_w, 0.0]);
    let solver = NeuralOdeSolver::new(func, 1e-3, 1e-6);
    let (z0, big_t, dt) = (1.0, 1.5, 0.01);
    let z_final = solver.solve_rk4(&[z0], 0.0, big_t, dt);
    let mut adjoint = AdjointMethod::new(1);
    let (grad_z0, grad_params) = adjoint.run(&solver, &z_final, &[1.0], 0.0, big_t, dt);
    let growth = (lambda * big_t).exp();
    assert!(
        (grad_z0[0] / growth - 1.0).abs() < 1e-4,
        "dL/dz0 = {}, e^(lambda T) = {growth}",
        grad_z0[0]
    );
    let expected_c = z0 * a_w * b_w * big_t * growth;
    assert!(
        (grad_params[5] / expected_c - 1.0).abs() < 1e-4,
        "dL/dc = {}, expected {expected_c}",
        grad_params[5]
    );
    let expected_a = z0 * b_w * c_w * big_t * growth;
    assert!(
        (grad_params[0] / expected_a - 1.0).abs() < 1e-4,
        "dL/da = {}, expected {expected_a}",
        grad_params[0]
    );
}

#[test]
fn test_adjoint_run_reconstructs_the_initial_state_and_stores_the_result() {
    let solver = NeuralOdeSolver::new(scaled_func(2, 6, 11, 5.0), 1e-3, 1e-6);
    let z0 = vec![0.4, -0.3];
    let z_final = solver.solve_rk4(&z0, 0.0, 1.0, 0.02);
    let mut adjoint = AdjointMethod::new(2);
    let (grad_z0, grad_params) = adjoint.run(&solver, &z_final, &[0.7, -1.1], 0.0, 1.0, 0.02);
    let n_params = solver.func.n_params();
    assert_eq!(adjoint.augmented_state.len(), 2 + 2 + n_params);
    assert!(max_abs_diff(&adjoint.augmented_state[..2], &z0) < 1e-8);
    assert_eq!(&adjoint.augmented_state[2..4], grad_z0.as_slice());
    assert_eq!(&adjoint.augmented_state[4..], grad_params.as_slice());
}

#[test]
fn test_adjoint_run_handles_a_step_that_does_not_divide_the_interval() {
    // dt = 0.3 over [0, 1]: steps 0.3, 0.3, 0.3, 0.1 (the grid of solve_rk4).
    let solver = NeuralOdeSolver::new(scaled_func(2, 6, 11, 5.0), 1e-3, 1e-6);
    let z0 = vec![0.4, -0.3];
    let g = vec![0.7, -1.1];
    let (t1, dt) = (1.0, 0.3);
    let z_final = solver.solve_rk4(&z0, 0.0, t1, dt);
    let mut adjoint = AdjointMethod::new(2);
    let (grad_z0, _) = adjoint.run(&solver, &z_final, &g, 0.0, t1, dt);
    let reference = central_difference(&|z: &[f64]| solve_loss(&solver, z, &g, t1, dt), &z0, 1e-6);
    assert!(max_abs_diff(&grad_z0, &reference) < 1e-5);
}

#[test]
fn test_adjoint_run_over_an_empty_interval_is_the_identity() {
    let solver = NeuralOdeSolver::new(scaled_func(2, 4, 3, 3.0), 1e-3, 1e-6);
    let mut adjoint = AdjointMethod::new(2);
    let (grad_z0, grad_params) = adjoint.run(&solver, &[0.1, 0.2], &[0.5, -0.5], 1.0, 1.0, 0.1);
    assert_eq!(grad_z0, vec![0.5, -0.5]);
    assert_eq!(grad_params.len(), solver.func.n_params());
    assert!(grad_params.iter().all(|&v| v == 0.0));
}

#[test]
fn test_adjoint_run_with_no_state_returns_zeros() {
    let solver = NeuralOdeSolver::new(scaled_func(2, 4, 3, 3.0), 1e-3, 1e-6);
    let mut adjoint = AdjointMethod::new(0);
    let (grad_z0, grad_params) = adjoint.run(&solver, &[], &[], 0.0, 1.0, 0.1);
    assert!(grad_z0.is_empty());
    assert!(grad_params.iter().all(|&v| v == 0.0));
}

/// Observations of a two-dimensional "teacher" network on `[0, 1]`.
fn teacher_series() -> (Vec<f64>, Vec<Vec<f64>>) {
    let teacher = NeuralOdeSolver::new(scaled_func(2, 8, 99, 15.0), 1e-3, 1e-6);
    let times: Vec<f64> = (0..11).map(|i| i as f64 * 0.1).collect();
    let observations = teacher.solve_rk4_trajectory(&[0.8, 0.3], &times, 0.01);
    (times, observations)
}

#[test]
fn test_time_series_loss_gradient_matches_finite_differences() {
    let (times, observations) = teacher_series();
    let student = NeuralOdeSolver::new(scaled_func(2, 4, 21, 4.0), 1e-3, 1e-6);
    let series = TimeSeriesOde::new(times, observations, student, 0.01, 0);
    let dt = 0.05;
    let gradient = series.loss_gradient(dt);
    assert_eq!(gradient.len(), series.solver.func.n_params());
    let reference = central_difference(
        &|theta: &[f64]| {
            let mut shifted = series.clone();
            shifted.solver.func.set_params_flat(theta);
            shifted.compute_loss(dt)
        },
        &series.solver.func.params_flat(),
        1e-6,
    );
    let worst = max_abs_diff(&gradient, &reference);
    assert!(worst < 1e-6, "loss_gradient deviates by {worst}");
    // Every layer receives a gradient (the old fit only touched bias_out).
    let magnitude = gradient.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(magnitude > 1e-2, "gradient is trivially small: {magnitude}");
    let n_in = series.solver.func.weights_in.len();
    assert!(gradient[..n_in].iter().any(|v| v.abs() > 1e-6));
}

#[test]
fn test_time_series_loss_gradient_with_unsorted_extras_and_empty_inputs() {
    let student = NeuralOdeSolver::new(scaled_func(1, 3, 4, 3.0), 1e-3, 1e-6);
    let empty = TimeSeriesOde::new(vec![], vec![], student.clone(), 0.01, 0);
    let zeros = empty.loss_gradient(0.1);
    assert_eq!(zeros.len(), student.func.n_params());
    assert!(zeros.iter().all(|&v| v == 0.0));
    // More observation times than observations: the surplus times are
    // integrated but never compared with data.
    let series = TimeSeriesOde::new(
        vec![0.0, 0.2, 0.4, 0.6],
        vec![vec![0.5], vec![0.6], vec![0.7]],
        student,
        0.01,
        0,
    );
    let gradient = series.loss_gradient(0.05);
    let reference = central_difference(
        &|theta: &[f64]| {
            let mut shifted = series.clone();
            shifted.solver.func.set_params_flat(theta);
            shifted.compute_loss(0.05)
        },
        &series.solver.func.params_flat(),
        1e-6,
    );
    assert!(max_abs_diff(&gradient, &reference) < 1e-7);
}

#[test]
fn test_time_series_fit_reduces_the_loss_on_teacher_data() {
    // The old fit only shrank bias_out towards zero and never changed a
    // weight, so the loss could not follow teacher-generated data.
    let (times, observations) = teacher_series();
    let student = NeuralOdeSolver::new(NeuralOdeFunc::new(2, 8, 15), 1e-3, 1e-6);
    let before = student.func.params_flat();
    let mut series = TimeSeriesOde::new(times, observations, student, 2.0, 200);
    series.fit();
    assert_eq!(series.loss_history.len(), 200);
    let first = series.loss_history[0];
    let last = series.loss_history[199];
    assert!(
        first > 1.0,
        "the untrained student must start far away: {first}"
    );
    assert!(last < 0.01 * first, "loss only went from {first} to {last}");
    for pair in series.loss_history.windows(2) {
        assert!(
            pair[1] <= pair[0] + 1e-12,
            "loss rose from {} to {}",
            pair[0],
            pair[1]
        );
    }
    // All layers moved, not just the output bias.
    let after = series.solver.func.params_flat();
    let n_in = series.solver.func.weights_in.len();
    assert!(
        max_abs_diff(&before[..n_in], &after[..n_in]) > 1e-3,
        "weights_in unchanged"
    );
    assert!(
        max_abs_diff(
            &series.solver.func.weights_hidden,
            &NeuralOdeFunc::new(2, 8, 15).weights_hidden
        ) > 1e-3,
        "weights_hidden unchanged"
    );
    let final_dt = 1.0 / (11.0 * 10.0);
    assert!(series.compute_loss(final_dt) < 0.01 * first);
}

#[test]
fn test_time_series_fit_stops_when_the_gradient_blows_up() {
    let (times, observations) = teacher_series();
    let student = NeuralOdeSolver::new(NeuralOdeFunc::new(2, 8, 15), 1e-3, 1e-6);
    let mut series = TimeSeriesOde::new(times, observations, student, 1e9, 200);
    series.fit();
    assert!(!series.loss_history.is_empty());
    assert!(
        series.loss_history.len() < 200,
        "divergence must end the run"
    );
    assert!(
        series
            .solver
            .func
            .params_flat()
            .iter()
            .all(|p| p.is_finite()),
        "parameters must stay finite"
    );
}

#[test]
fn test_time_series_fit_with_zero_iterations_changes_nothing() {
    let (times, observations) = teacher_series();
    let student = NeuralOdeSolver::new(NeuralOdeFunc::new(2, 8, 15), 1e-3, 1e-6);
    let before = student.func.params_flat();
    let mut series = TimeSeriesOde::new(times, observations, student, 1.0, 0);
    series.fit();
    assert!(series.loss_history.is_empty());
    assert_eq!(series.solver.func.params_flat(), before);
}
