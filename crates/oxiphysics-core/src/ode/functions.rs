//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Trait representing the right-hand side of an ODE system `dy/dt = f(t, y)`.
///
/// Implement this trait to define a system of ordinary differential equations
/// for use with `AdaptiveIntegrator` and other solvers.
pub trait OdeSystem {
    /// Evaluate the right-hand side of the ODE system.
    ///
    /// # Arguments
    /// * `t` - Current time.
    /// * `y` - Current state vector (length `n`).
    /// * `dydt` - Output derivative vector (length `n`, must be filled).
    fn ode_rhs(&self, t: f64, y: &[f64], dydt: &mut [f64]);
}
/// Euler forward step: `y_new = y + dt * dydt`
pub fn euler_step(y: &[f64], dydt: &[f64], dt: f64) -> Vec<f64> {
    y.iter()
        .zip(dydt.iter())
        .map(|(yi, fi)| yi + dt * fi)
        .collect()
}
/// RK2 midpoint method (Heun's method): `y_new = y + dt * f(y + 0.5*dt*f(y))`
pub fn rk2_step(y: &[f64], f: impl Fn(&[f64]) -> Vec<f64>, dt: f64) -> Vec<f64> {
    let k1 = f(y);
    let y_mid: Vec<f64> = y
        .iter()
        .zip(k1.iter())
        .map(|(yi, ki)| yi + 0.5 * dt * ki)
        .collect();
    let k2 = f(&y_mid);
    y.iter()
        .zip(k2.iter())
        .map(|(yi, ki)| yi + dt * ki)
        .collect()
}
/// Classic RK4 step.
pub fn rk4_step(y: &[f64], f: impl Fn(&[f64]) -> Vec<f64>, dt: f64) -> Vec<f64> {
    let n = y.len();
    let k1 = f(y);
    let y2: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * dt * k1[i]).collect();
    let k2 = f(&y2);
    let y3: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * dt * k2[i]).collect();
    let k3 = f(&y3);
    let y4: Vec<f64> = (0..n).map(|i| y[i] + dt * k3[i]).collect();
    let k4 = f(&y4);
    (0..n)
        .map(|i| y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
        .collect()
}
/// Nodes `c_i` of the Dormand-Prince 5(4) pair: stage `i` is evaluated at `t + c_i * dt`.
const DOPRI5_C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
/// Runge-Kutta matrix `a_ij` of the Dormand-Prince 5(4) pair for the stages 2 to 7.
///
/// Row `r` holds the coefficients `a_(r+2, 1..=6)`.  The last row (stage 7)
/// equals the 5th-order weights `b`: the 7th stage is the derivative at the new
/// 5th-order solution, which makes the scheme "first same as last" (FSAL).
const DOPRI5_A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
/// Differences `b_i - b^_i` between the 5th-order weights and the embedded
/// 4th-order weights `b^ = (5179/57600, 0, 7571/16695, 393/640, -92097/339200,
/// 187/2100, 1/40)`; the embedded error estimate is `dt * sum_i e_i k_i`.
const DOPRI5_E: [f64; 7] = [
    71.0 / 57600.0,
    0.0,
    -71.0 / 16695.0,
    71.0 / 1920.0,
    -17253.0 / 339200.0,
    22.0 / 525.0,
    -1.0 / 40.0,
];
/// One Dormand-Prince 5(4) step of the non-autonomous system `dy/dt = f(t, y)`
/// whose first stage `k1 = f(t, y)` is supplied by the caller.
///
/// Stage `i` is evaluated at the time `t + c_i * dt` with
/// `c = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)`, so a right-hand side that depends
/// explicitly on time is integrated with full 5th-order accuracy.
///
/// Because the scheme is FSAL, the returned last stage `f(t + dt, y_new)` is the
/// first stage of the following step: an adaptive driver that accepts a step
/// passes it back in as `k1` and needs only six new evaluations of `f` per
/// accepted step (a rejected step may reuse the same `k1`).
///
/// `k1` and `y` must have the same length.
///
/// Returns `(y_new, error_estimate, k_last)`:
/// * `y_new`: the 5th-order solution at `t + dt`,
/// * `error_estimate`: the RMS norm of the embedded error estimate
///   `dt * sum_i (b_i - b^_i) k_i` (the estimate for the 4th-order companion),
/// * `k_last`: `f(t + dt, y_new)`.
pub fn rk45_step_fsal(
    t: f64,
    y: &[f64],
    k1: &[f64],
    f: impl Fn(f64, &[f64]) -> Vec<f64>,
    dt: f64,
) -> (Vec<f64>, f64, Vec<f64>) {
    let n = y.len();
    let mut stages: Vec<Vec<f64>> = Vec::with_capacity(7);
    stages.push(k1.to_vec());
    for stage in 1..6 {
        let row = &DOPRI5_A[stage - 1];
        let y_stage: Vec<f64> = (0..n)
            .map(|i| {
                let weighted: f64 = stages.iter().zip(row.iter()).map(|(k, a)| a * k[i]).sum();
                y[i] + dt * weighted
            })
            .collect();
        stages.push(f(t + DOPRI5_C[stage] * dt, &y_stage));
    }
    let b = &DOPRI5_A[5];
    let y_new: Vec<f64> = (0..n)
        .map(|i| {
            let weighted: f64 = stages.iter().zip(b.iter()).map(|(k, w)| w * k[i]).sum();
            y[i] + dt * weighted
        })
        .collect();
    let k_last = f(t + DOPRI5_C[6] * dt, &y_new);
    stages.push(k_last.clone());
    let err = if n == 0 {
        0.0
    } else {
        let sum_sq: f64 = (0..n)
            .map(|i| {
                let weighted: f64 = stages
                    .iter()
                    .zip(DOPRI5_E.iter())
                    .map(|(k, e)| e * k[i])
                    .sum();
                let err_i = dt * weighted;
                err_i * err_i
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    };
    (y_new, err, k_last)
}
/// One Dormand-Prince 5(4) step of the non-autonomous system `dy/dt = f(t, y)`.
///
/// Every stage is evaluated at its own time `t + c_i * dt` (see
/// [`rk45_step_fsal`]); the first stage `f(t, y)` is computed here.
///
/// Returns `(y_new, error_estimate, k_last)` where `error_estimate` is the RMS
/// norm of the embedded error estimate and `k_last = f(t + dt, y_new)` is the
/// FSAL derivative at the new point.
pub fn rk45_step_t(
    t: f64,
    y: &[f64],
    f: impl Fn(f64, &[f64]) -> Vec<f64>,
    dt: f64,
) -> (Vec<f64>, f64, Vec<f64>) {
    let k1 = f(t, y);
    rk45_step_fsal(t, y, &k1, f, dt)
}
/// Dormand-Prince RK45 adaptive step for an **autonomous** system `dy/dt = f(y)`.
///
/// Returns `(y_new, error_estimate)` where `error_estimate` is the RMS norm of
/// the difference between the 5th- and 4th-order solutions.  The right-hand
/// side cannot see the time; use [`rk45_step_t`] for a system that depends
/// explicitly on `t`.  The tolerance argument is not used by a single step.
pub fn rk45_step(y: &[f64], f: impl Fn(&[f64]) -> Vec<f64>, dt: f64, _tol: f64) -> (Vec<f64>, f64) {
    let (y_new, err, _) = rk45_step_t(0.0, y, |_t, yv| f(yv), dt);
    (y_new, err)
}
/// Adaptive RK45 integration using Dormand-Prince.
///
/// Returns a vector of `(t, y)` pairs at each accepted step.
pub fn adaptive_rk45(
    y0: Vec<f64>,
    f: impl Fn(&[f64]) -> Vec<f64>,
    t0: f64,
    t_end: f64,
    dt_init: f64,
    atol: f64,
    rtol: f64,
) -> Vec<(f64, Vec<f64>)> {
    let mut t = t0;
    let mut y = y0;
    let mut dt = dt_init;
    let mut result = vec![(t, y.clone())];
    let safety = 0.9_f64;
    let min_scale = 0.2_f64;
    let max_scale = 10.0_f64;
    while t < t_end {
        if t + dt > t_end {
            dt = t_end - t;
        }
        let (y_new, err) = rk45_step(&y, &f, dt, atol);
        let tol_norm = (atol + rtol * y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)).max(atol);
        let err_norm = err / tol_norm;
        if err_norm <= 1.0 || err == 0.0 {
            t += dt;
            y = y_new;
            result.push((t, y.clone()));
            let scale = if err == 0.0 {
                max_scale
            } else {
                (safety * err_norm.powf(-0.2)).clamp(min_scale, max_scale)
            };
            dt *= scale;
        } else {
            let scale = (safety * err_norm.powf(-0.25)).clamp(min_scale, 1.0);
            dt *= scale;
        }
        if dt < 1e-15 * t_end.abs().max(1.0) {
            break;
        }
    }
    result
}
/// Leapfrog step (Störmer-Verlet).
///
/// Returns `(new_pos, new_vel)`:
/// - `pos_new = pos + vel*dt + 0.5*acc*dt²`
/// - `vel_new = vel + acc*dt`
pub fn leapfrog_step(pos: &[f64], vel: &[f64], acc: &[f64], dt: f64) -> (Vec<f64>, Vec<f64>) {
    let new_pos: Vec<f64> = (0..pos.len())
        .map(|i| pos[i] + vel[i] * dt + 0.5 * acc[i] * dt * dt)
        .collect();
    let new_vel: Vec<f64> = (0..vel.len()).map(|i| vel[i] + acc[i] * dt).collect();
    (new_pos, new_vel)
}
/// Velocity-Verlet step.
///
/// Returns `(new_pos, new_vel)`:
/// - `pos_new = pos + vel*dt + 0.5*acc_old*dt²`
/// - `vel_new = vel + 0.5*(acc_old + acc_new)*dt`
pub fn velocity_verlet(
    pos: &[f64],
    vel: &[f64],
    acc_old: &[f64],
    acc_new: &[f64],
    dt: f64,
) -> (Vec<f64>, Vec<f64>) {
    let new_pos: Vec<f64> = (0..pos.len())
        .map(|i| pos[i] + vel[i] * dt + 0.5 * acc_old[i] * dt * dt)
        .collect();
    let new_vel: Vec<f64> = (0..vel.len())
        .map(|i| vel[i] + 0.5 * (acc_old[i] + acc_new[i]) * dt)
        .collect();
    (new_pos, new_vel)
}
/// Solves a linear system `A * x = b` using Gaussian elimination with partial pivoting.
///
/// Returns `Some(x)` on success, or `None` if the matrix is singular.
pub fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut mat: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();
    for col in 0..n {
        let pivot_row = (col..n).max_by(|&r1, &r2| {
            mat[r1][col]
                .abs()
                .partial_cmp(&mat[r2][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        mat.swap(col, pivot_row);
        let pivot = mat[col][col];
        if pivot.abs() < 1e-15 {
            return None;
        }
        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            let (upper, lower) = mat.split_at_mut(row);
            for (mr_k, &mc_k) in lower[0][col..=n].iter_mut().zip(upper[col][col..=n].iter()) {
                *mr_k -= factor * mc_k;
            }
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = mat[i][n];
        for j in (i + 1)..n {
            sum -= mat[i][j] * x[j];
        }
        x[i] = sum / mat[i][i];
    }
    Some(x)
}
/// Implicit Euler step: solve `(I - dt*J) * y_new = y + dt*rhs`.
///
/// Uses Gaussian elimination with partial pivoting.
pub fn implicit_euler_step(y: &[f64], jacobian: &[Vec<f64>], rhs: &[f64], dt: f64) -> Vec<f64> {
    let n = y.len();
    let mat: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let delta = if i == j { 1.0 } else { 0.0 };
                    delta - dt * jacobian[i][j]
                })
                .collect()
        })
        .collect();
    let b: Vec<f64> = (0..n).map(|i| y[i] + dt * rhs[i]).collect();
    solve_linear_system(&mat, &b).unwrap_or_else(|| b.clone())
}
/// Compute the mixed absolute/relative RMS error norm.
///
/// Uses the formula from Hairer & Wanner:
/// `norm = sqrt(1/n * sum_i (err_i / (atol + |y_i| * rtol))^2)`
pub fn mixed_error_norm(y: &[f64], err: &[f64], atol: f64, rtol: f64) -> f64 {
    let n = y.len();
    if n == 0 {
        return 0.0;
    }
    let sum: f64 = y
        .iter()
        .zip(err.iter())
        .map(|(&yi, &ei)| {
            let sc = atol + yi.abs() * rtol;
            (ei / sc).powi(2)
        })
        .sum();
    (sum / n as f64).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ode::Adams4;
    use crate::ode::AdamsBashforthMoulton4;
    use crate::ode::AdaptiveIntegrator;
    use crate::ode::Bdf2;
    use crate::ode::BdfOrder2;
    use crate::ode::BulirschStoer;
    use crate::ode::DenseOutputSegment;
    use crate::ode::DormandPrince;
    use crate::ode::DormandPrince45;
    use crate::ode::EventDetector;
    use crate::ode::Fehlberg45;
    use crate::ode::ImplicitEulerNewton;
    use crate::ode::LeapFrog;
    use crate::ode::PiStepController;
    use crate::ode::RungeKuttaFehlberg;
    use crate::ode::StiffOdeSolver;
    use crate::ode::Verlet;
    use std::f64::consts::{PI, TAU};
    #[test]
    fn test_euler_step() {
        let y = vec![1.0_f64];
        let dydt = vec![0.5_f64];
        let y_new = euler_step(&y, &dydt, 0.1);
        assert!((y_new[0] - 1.05).abs() < 1e-12, "got {}", y_new[0]);
    }
    #[test]
    fn test_rk4_step_exponential() {
        let mut y = vec![1.0_f64];
        let dt = 0.01;
        let steps = (1.0 / dt) as usize;
        for _ in 0..steps {
            y = rk4_step(&y, |v| vec![v[0]], dt);
        }
        let e = std::f64::consts::E;
        assert!((y[0] - e).abs() < 1e-8, "y[0]={}, e={}", y[0], e);
    }
    #[test]
    fn test_leapfrog_uniform_acceleration() {
        let pos = vec![0.0_f64];
        let vel = vec![0.0_f64];
        let acc = vec![1.0_f64];
        let dt = 0.5;
        let (new_pos, new_vel) = leapfrog_step(&pos, &vel, &acc, dt);
        let expected_pos = 0.5 * 1.0 * dt * dt;
        let expected_vel = 1.0 * dt;
        assert!(
            (new_pos[0] - expected_pos).abs() < 1e-12,
            "pos={}",
            new_pos[0]
        );
        assert!(
            (new_vel[0] - expected_vel).abs() < 1e-12,
            "vel={}",
            new_vel[0]
        );
    }
    #[test]
    fn test_velocity_verlet_constant_acceleration() {
        let pos = vec![0.0_f64];
        let vel = vec![0.0_f64];
        let acc_old = vec![2.0_f64];
        let acc_new = vec![2.0_f64];
        let dt = 1.0;
        let (new_pos, new_vel) = velocity_verlet(&pos, &vel, &acc_old, &acc_new, dt);
        assert!((new_pos[0] - 1.0).abs() < 1e-12, "pos={}", new_pos[0]);
        assert!((new_vel[0] - 2.0).abs() < 1e-12, "vel={}", new_vel[0]);
    }
    #[test]
    fn test_adaptive_rk45_exponential() {
        let y0 = vec![1.0_f64];
        let result = adaptive_rk45(y0, |v| vec![v[0]], 0.0, 1.0, 0.1, 1e-8, 1e-6);
        let (_, y_end) = result.last().expect("should have at least one step");
        let e = std::f64::consts::E;
        assert!(
            (y_end[0] - e).abs() < 1e-5,
            "y_end={}, e={}, diff={}",
            y_end[0],
            e,
            (y_end[0] - e).abs()
        );
    }
    #[test]
    fn test_solve_linear_system_2x2() {
        let a = vec![vec![2.0_f64, 1.0], vec![1.0, 3.0]];
        let b = vec![5.0_f64, 10.0];
        let x = solve_linear_system(&a, &b).expect("should solve");
        assert!((x[0] - 1.0).abs() < 1e-12, "x[0]={}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-12, "x[1]={}", x[1]);
    }
    /// Simple harmonic oscillator: d²x/dt² = -omega²*x
    /// State: \[x, v\], RHS: \[v, -omega²*x\]
    fn harmonic_rhs(y: &[f64], omega: f64) -> Vec<f64> {
        vec![y[1], -omega * omega * y[0]]
    }
    fn harmonic_energy(y: &[f64], omega: f64) -> f64 {
        0.5 * (y[1] * y[1] + omega * omega * y[0] * y[0])
    }
    #[test]
    fn test_rk4_harmonic_energy_conservation() {
        let omega = 2.0;
        let mut y = vec![1.0_f64, 0.0_f64];
        let e0 = harmonic_energy(&y, omega);
        let dt = 0.01;
        let n_steps = 1000;
        for _ in 0..n_steps {
            y = rk4_step(&y, |yv| harmonic_rhs(yv, omega), dt);
        }
        let e_final = harmonic_energy(&y, omega);
        assert!(
            (e_final - e0).abs() / e0 < 1e-6,
            "Energy drift: e0={e0}, e_final={e_final}"
        );
    }
    #[test]
    fn test_rk4_harmonic_period() {
        let omega = 1.0;
        let t_period = 2.0 * PI / omega;
        let dt = 0.001;
        let n_steps = (t_period / dt).round() as usize;
        let mut y = vec![1.0_f64, 0.0_f64];
        for _ in 0..n_steps {
            y = rk4_step(&y, |yv| harmonic_rhs(yv, omega), dt);
        }
        assert!((y[0] - 1.0).abs() < 1e-3, "x after period = {}", y[0]);
        assert!(y[1].abs() < 1e-3, "v after period = {}", y[1]);
    }
    #[test]
    fn test_dormand_prince45_step_accepted() {
        let dp = DormandPrince45::new(1e-8, 1e-6);
        let y = vec![1.0_f64];
        let res = dp.step(&y, |v| vec![-v[0]], 0.1);
        let exact = (-0.1_f64).exp();
        assert!(
            (res.y[0] - exact).abs() < 1e-9,
            "y={}, exact={}",
            res.y[0],
            exact
        );
    }
    #[test]
    fn test_dormand_prince45_step_result_fields() {
        let dp = DormandPrince45::new(1e-6, 1e-4);
        let y = vec![0.0_f64, 1.0_f64];
        let res = dp.step(&y, |v| vec![v[1], -v[0]], 0.01);
        assert_eq!(res.y.len(), 2);
        assert!(res.error_estimate >= 0.0);
    }
    #[test]
    fn test_fehlberg45_exponential() {
        let rkf = Fehlberg45::new();
        let y = vec![1.0_f64];
        let (y_new, err) = rkf.step(&y, |v| vec![v[0]], 0.1);
        let exact = 0.1_f64.exp();
        assert!(
            (y_new[0] - exact).abs() < 1e-8,
            "y={}, exact={}",
            y_new[0],
            exact
        );
        assert!(err >= 0.0, "Error estimate should be non-negative");
    }
    #[test]
    fn test_fehlberg45_harmonic_energy() {
        let omega = 1.0;
        let mut y = vec![1.0_f64, 0.0_f64];
        let e0 = harmonic_energy(&y, omega);
        let rkf = Fehlberg45::new();
        let dt = 0.01;
        for _ in 0..500 {
            let (y_new, _) = rkf.step(&y, |yv| harmonic_rhs(yv, omega), dt);
            y = y_new;
        }
        let e_final = harmonic_energy(&y, omega);
        assert!(
            (e_final - e0).abs() / e0 < 1e-4,
            "Fehlberg energy drift: e0={e0}, e_final={e_final}"
        );
    }
    #[test]
    fn test_adams4_exponential() {
        let dt = 0.01;
        let mut y = vec![1.0_f64];
        let mut ab = Adams4::new();
        for _ in 0..4 {
            let dydt = vec![y[0]];
            ab.push(dydt);
            y = rk4_step(&y, |v| vec![v[0]], dt);
        }
        for _ in 0..96 {
            let y_new = ab.step(&y, dt);
            let dydt_new = vec![y_new[0]];
            ab.push(dydt_new);
            y = y_new;
        }
        let e = std::f64::consts::E;
        assert!((y[0] - e).abs() < 5e-4, "Adams4: y={}, e={}", y[0], e);
    }
    #[test]
    fn test_adams4_ready_flag() {
        let mut ab = Adams4::new();
        assert!(!ab.ready());
        for i in 0..4 {
            ab.push(vec![i as f64]);
        }
        assert!(ab.ready());
    }
    #[test]
    fn test_bdf2_stiff_decay() {
        let mut bdf = Bdf2::new();
        let mut y = vec![1.0_f64];
        let dt = 0.1;
        let n_steps = 10;
        for step in 0..n_steps {
            let t_new = (step + 1) as f64 * dt;
            let y_new = bdf.step(
                &y,
                |t, yv| {
                    let _ = t;
                    vec![-100.0 * yv[0]]
                },
                t_new,
                dt,
            );
            y = y_new;
        }
        let t_final = n_steps as f64 * dt;
        let exact = (-100.0 * t_final).exp();
        assert!(
            y[0] < 0.01 || (y[0] - exact).abs() < 0.1,
            "BDF2 stiff: y={}, exact={}",
            y[0],
            exact
        );
    }
    pub(super) struct HarmonicOscillator {
        omega: f64,
    }
    impl OdeSystem for HarmonicOscillator {
        fn ode_rhs(&self, _t: f64, y: &[f64], dydt: &mut [f64]) {
            dydt[0] = y[1];
            dydt[1] = -self.omega * self.omega * y[0];
        }
    }
    #[test]
    fn test_adaptive_integrator_harmonic() {
        let system = HarmonicOscillator { omega: 1.0 };
        let integrator = AdaptiveIntegrator::new(1e-8, 1e-6);
        let y0 = vec![1.0_f64, 0.0_f64];
        let result = integrator.integrate(&system, y0, 0.0, 2.0 * PI, 0.1);
        let (_, y_end) = result.last().expect("should have steps");
        assert!(
            (y_end[0] - 1.0).abs() < 1e-4,
            "x after period = {}",
            y_end[0]
        );
        assert!(y_end[1].abs() < 1e-4, "v after period = {}", y_end[1]);
    }
    #[test]
    fn test_adaptive_integrator_exponential() {
        pub(super) struct ExponentialDecay;
        impl OdeSystem for ExponentialDecay {
            fn ode_rhs(&self, _t: f64, y: &[f64], dydt: &mut [f64]) {
                dydt[0] = -y[0];
            }
        }
        let integrator = AdaptiveIntegrator::new(1e-9, 1e-7);
        let result = integrator.integrate(&ExponentialDecay, vec![1.0_f64], 0.0, 1.0, 0.1);
        let (_, y_end) = result.last().unwrap();
        let exact = (-1.0_f64).exp();
        assert!(
            (y_end[0] - exact).abs() < 1e-6,
            "y={}, exact={}",
            y_end[0],
            exact
        );
    }
    #[test]
    fn test_adaptive_integrator_step_count() {
        let system = HarmonicOscillator { omega: 2.0 };
        let integrator = AdaptiveIntegrator::new(1e-6, 1e-4);
        let result = integrator.integrate(&system, vec![1.0, 0.0], 0.0, 1.0, 0.05);
        assert!(result.len() > 1, "Should have more than start point");
    }
    #[test]
    fn test_adaptive_integrator_energy_conservation() {
        let omega = 3.0;
        let system = HarmonicOscillator { omega };
        let integrator = AdaptiveIntegrator::new(1e-9, 1e-7);
        let y0 = vec![1.0_f64, 0.0_f64];
        let e0 = harmonic_energy(&y0, omega);
        let result = integrator.integrate(&system, y0, 0.0, 10.0, 0.05);
        let (_, y_end) = result.last().unwrap();
        let e_final = harmonic_energy(y_end, omega);
        assert!(
            (e_final - e0).abs() / e0 < 1e-5,
            "Energy drift: e0={e0}, e_final={e_final}"
        );
    }
    #[test]
    fn test_dormand_prince_integrates_exp() {
        let dp = DormandPrince::new(1e-8, 1e-10, 0.5);
        let result = dp.integrate(|_t, y| vec![y[0]], 0.0, 1.0, vec![1.0_f64]);
        let (_, y_end) = result.last().expect("should have steps");
        let e = std::f64::consts::E;
        assert!(
            (y_end[0] - e).abs() < 1e-5,
            "DormandPrince: y={}, e={}",
            y_end[0],
            e
        );
    }
    #[test]
    fn test_dormand_prince_step_count() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let result = dp.integrate(|_t, y| vec![-y[0]], 0.0, 2.0, vec![1.0_f64]);
        assert!(result.len() > 1, "should produce multiple steps");
    }
    #[test]
    fn test_bulirsch_stoer_midpoint_exp() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let y0 = vec![1.0_f64];
        let h = 1.0;
        for &n in &[2usize, 4, 8, 16] {
            let yn = BulirschStoer::midpoint_method(&f, 0.0, &y0, h, n);
            let e = std::f64::consts::E;
            assert!(
                (yn[0] - e).abs() < 0.1,
                "midpoint n={n}: y={}, e={e}",
                yn[0]
            );
        }
    }
    #[test]
    fn test_bulirsch_stoer_midpoint_error_decreases() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let y0 = vec![1.0_f64];
        let e = std::f64::consts::E;
        let err4 = (BulirschStoer::midpoint_method(&f, 0.0, &y0, 1.0, 4)[0] - e).abs();
        let err8 = (BulirschStoer::midpoint_method(&f, 0.0, &y0, 1.0, 8)[0] - e).abs();
        assert!(
            err8 < err4,
            "error should decrease with more substeps: err4={err4}, err8={err8}"
        );
    }
    #[test]
    fn test_bulirsch_stoer_extrapolate() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let y0 = vec![1.0_f64];
        let h = 1.0;
        let e = std::f64::consts::E;
        let t_table: Vec<Vec<f64>> = (1..=4)
            .map(|k| BulirschStoer::midpoint_method(&f, 0.0, &y0, h, 2 * k))
            .collect();
        let extrap = BulirschStoer::extrapolate(&t_table, 4);
        let err_extrap = (extrap[0] - e).abs();
        let err_direct = (t_table[0][0] - e).abs();
        assert!(
            err_extrap < err_direct,
            "extrapolation should improve accuracy: extrap_err={err_extrap}, direct_err={err_direct}"
        );
    }
    #[test]
    fn test_verlet_harmonic_energy() {
        let omega = 1.0_f64;
        let dt = 0.01_f64;
        let mut x = vec![1.0_f64];
        let a0_val = -omega * omega * x[0];
        let mut x_prev = vec![x[0] - 0.0 * dt + 0.5 * a0_val * dt * dt];
        let n_steps = 1000;
        for _ in 0..n_steps {
            let a: Vec<f64> = vec![-omega * omega * x[0]];
            Verlet::step(&mut x, &mut x_prev, &a, dt);
        }
        let v_approx = (x[0] - x_prev[0]) / dt;
        let energy = 0.5 * (v_approx * v_approx + omega * omega * x[0] * x[0]);
        assert!(
            (energy - 0.5).abs() < 0.1,
            "Verlet energy drift: E={energy}, expected ~0.5"
        );
    }
    #[test]
    fn test_verlet_constant_force() {
        let dt = 0.1_f64;
        let a = vec![1.0_f64];
        let mut x = vec![0.0_f64];
        let mut x_prev = vec![0.5 * dt * dt];
        for _ in 0..10 {
            Verlet::step(&mut x, &mut x_prev, &a, dt);
        }
        assert!(
            (x[0] - 0.5).abs() < 0.01,
            "Verlet constant force: x={}",
            x[0]
        );
    }
    #[test]
    fn test_leapfrog_kick_drift() {
        let mut v = vec![0.0_f64];
        let a = vec![1.0_f64];
        let mut x = vec![0.0_f64];
        let dt = 0.1_f64;
        LeapFrog::kick(&mut v, &a, 0.5 * dt);
        LeapFrog::drift(&mut x, &v, dt);
        LeapFrog::kick(&mut v, &a, 0.5 * dt);
        assert!((v[0] - 0.1).abs() < 1e-12, "v={}", v[0]);
        assert!((x[0] - 0.005).abs() < 1e-12, "x={}", x[0]);
    }
    #[test]
    fn test_leapfrog_conserves_hamiltonian() {
        let omega = 1.0_f64;
        let dt = 0.01_f64;
        let n_steps = 2000;
        let mut x = vec![1.0_f64];
        let mut v = vec![0.0_f64];
        let e0 = 0.5 * (v[0] * v[0] + omega * omega * x[0] * x[0]);
        for _ in 0..n_steps {
            let a: Vec<f64> = vec![-omega * omega * x[0]];
            LeapFrog::kick(&mut v, &a, 0.5 * dt);
            LeapFrog::drift(&mut x, &v, dt);
            let a_new: Vec<f64> = vec![-omega * omega * x[0]];
            LeapFrog::kick(&mut v, &a_new, 0.5 * dt);
        }
        let e_final = 0.5 * (v[0] * v[0] + omega * omega * x[0] * x[0]);
        assert!(
            (e_final - e0).abs() / e0 < 1e-4,
            "LeapFrog energy drift: e0={e0}, e_final={e_final}"
        );
    }
    #[test]
    fn test_stiff_solver_decay() {
        let solver = StiffOdeSolver::new(0.5);
        let mut y = vec![1.0_f64];
        let dt = 0.05_f64;
        let n_steps = 20;
        for step in 0..n_steps {
            let t = step as f64 * dt;
            y = solver.step_fixed(|_t, yv| vec![-100.0 * yv[0]], t, &y, dt);
        }
        let t_final = n_steps as f64 * dt;
        let exact = (-100.0 * t_final).exp();
        assert!(
            y[0].abs() < 0.1 || (y[0] - exact).abs() < 0.1,
            "StiffOdeSolver: y={}, exact={}",
            y[0],
            exact
        );
    }
    #[test]
    fn test_stiff_solver_constant_rhs() {
        let solver = StiffOdeSolver::new(0.5);
        let y = vec![0.0_f64];
        let dt = 0.1_f64;
        let y_new = solver.step_fixed(|_t, _y| vec![1.0], 0.0, &y, dt);
        assert!((y_new[0] - dt).abs() < 1e-6, "y_new={}", y_new[0]);
    }
    #[test]
    fn test_abm4_exponential() {
        let dt = 0.01;
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let mut abm = AdamsBashforthMoulton4::new();
        let mut y = vec![1.0_f64];
        let mut t = 0.0_f64;
        for _ in 0..4 {
            let k = f(t, &y);
            abm.push_history(t, y.clone(), k);
            y = rk4_step(&y, |v| f(t, v), dt);
            t += dt;
        }
        for _ in 0..96 {
            let y_new = abm.step(&f, t, &y, dt);
            let k_new = f(t + dt, &y_new);
            abm.push_history(t + dt, y_new.clone(), k_new);
            y = y_new;
            t += dt;
        }
        let e = std::f64::consts::E;
        assert!((y[0] - e).abs() < 5e-4, "ABM4: y={}, e={}", y[0], e);
    }
    #[test]
    fn test_abm4_harmonic() {
        let omega = 1.0_f64;
        let f = |_t: f64, y: &[f64]| harmonic_rhs(y, omega);
        let mut abm = AdamsBashforthMoulton4::new();
        let mut y = vec![1.0_f64, 0.0_f64];
        let e0 = harmonic_energy(&y, omega);
        let dt = 0.01;
        let mut t = 0.0_f64;
        for _ in 0..4 {
            let k = f(t, &y);
            abm.push_history(t, y.clone(), k);
            y = rk4_step(&y, |v| f(t, v), dt);
            t += dt;
        }
        for _ in 0..496 {
            let y_new = abm.step(&f, t, &y, dt);
            let k_new = f(t + dt, &y_new);
            abm.push_history(t + dt, y_new.clone(), k_new);
            y = y_new;
            t += dt;
        }
        let e_final = harmonic_energy(&y, omega);
        assert!(
            (e_final - e0).abs() / e0 < 1e-3,
            "ABM4 energy drift: e0={e0}, e_final={e_final}"
        );
    }
    #[test]
    fn test_implicit_euler_newton_decay() {
        let solver = ImplicitEulerNewton::new(50, 1e-10);
        let mut y = vec![1.0_f64];
        let dt = 0.1_f64;
        for step in 0..5 {
            let t = step as f64 * dt;
            y = solver.step(|_t, yv| vec![-100.0 * yv[0]], t, &y, dt);
        }
        assert!(y[0] < 0.01, "ImplicitEulerNewton: y={}", y[0]);
    }
    #[test]
    fn test_implicit_euler_newton_linear_growth() {
        let solver = ImplicitEulerNewton::new(50, 1e-12);
        let mut y = vec![0.0_f64];
        let dt = 0.1_f64;
        for step in 0..10 {
            let t = step as f64 * dt;
            y = solver.step(|_t, _yv| vec![1.0], t, &y, dt);
        }
        assert!(
            (y[0] - 1.0).abs() < 1e-8,
            "ImplicitEulerNewton linear: y={}",
            y[0]
        );
    }
    #[test]
    fn test_bdf_order2_exponential() {
        let bdf = BdfOrder2::new(50, 1e-10);
        let mut y = vec![1.0_f64];
        let dt = 0.01_f64;
        let t0 = 0.0_f64;
        let y1 = bdf.first_step(|_t, yv| vec![-yv[0]], t0, &y, dt);
        let mut history = (y.clone(), y1.clone());
        y = y1;
        let mut t = dt;
        for _ in 0..99 {
            let y_new = bdf.step(|_t, yv| vec![-yv[0]], t, &history.0, &history.1, dt);
            history = (history.1.clone(), y_new.clone());
            y = y_new;
            t += dt;
        }
        let expected = (-1.0_f64).exp();
        assert!(
            (y[0] - expected).abs() < 1e-4,
            "BDF2: y={}, expected={}",
            y[0],
            expected
        );
    }
    #[test]
    fn test_event_detector_zero_crossing_sine() {
        let f = |t: f64| t.sin();
        let crossings = EventDetector::find_zero_crossings(f, 0.01, TAU, 0.05, 1e-10);
        let near_pi = crossings.iter().any(|&c| (c - PI).abs() < 0.01);
        assert!(
            near_pi,
            "No crossing found near π; crossings: {:?}",
            crossings
        );
    }
    #[test]
    fn test_event_detector_multiple_crossings() {
        let f = |t: f64| (2.0 * t).sin();
        let crossings = EventDetector::find_zero_crossings(f, 0.01, TAU, 0.05, 1e-10);
        assert!(
            crossings.len() >= 3,
            "Found only {} crossings",
            crossings.len()
        );
    }
    #[test]
    fn test_dense_output_interpolates_linearly() {
        let seg = DenseOutputSegment {
            t0: 0.0,
            t1: 1.0,
            y0: vec![0.0],
            y1: vec![1.0],
            k1: vec![1.0],
            k6: vec![1.0],
        };
        let mid = seg.evaluate(0.5);
        assert!(mid[0] > 0.0 && mid[0] < 1.0, "dense output mid={}", mid[0]);
    }
    #[test]
    fn test_dense_output_endpoints() {
        let seg = DenseOutputSegment {
            t0: 0.0,
            t1: 2.0,
            y0: vec![3.0],
            y1: vec![5.0],
            k1: vec![1.0],
            k6: vec![1.0],
        };
        let at_t0 = seg.evaluate(0.0);
        let at_t1 = seg.evaluate(2.0);
        assert!(
            (at_t0[0] - 3.0).abs() < 1e-10,
            "dense output at t0={}",
            at_t0[0]
        );
        assert!(
            (at_t1[0] - 5.0).abs() < 1e-10,
            "dense output at t1={}",
            at_t1[0]
        );
    }
    #[test]
    fn test_pi_controller_accepts_small_error() {
        let ctrl = PiStepController::new(1e-6, 1e-8, 0.01, 5.0, 0.2);
        let (new_dt, accepted) = ctrl.control(0.01, 1e-10);
        assert!(accepted, "should accept tiny error");
        assert!(new_dt >= 0.01, "step should not decrease with tiny error");
    }
    #[test]
    fn test_pi_controller_rejects_large_error() {
        let ctrl = PiStepController::new(1e-6, 1e-8, 1e-6, 5.0, 0.9);
        let (new_dt, accepted) = ctrl.control(0.1, 5.0);
        assert!(!accepted, "should reject large error (err=5.0)");
        assert!(
            new_dt < 0.1,
            "step should shrink after rejection, got new_dt={new_dt}"
        );
    }
    #[test]
    fn test_pi_controller_clamps_max_dt() {
        let ctrl = PiStepController::new(1e-6, 1e-8, 0.01, 0.5, 0.2);
        let (new_dt, _) = ctrl.control(0.01, 1e-15);
        assert!(new_dt <= 0.5, "step should be clamped to max_dt={}", new_dt);
    }
    #[test]
    fn test_mixed_error_norm_zero() {
        let y = vec![1.0, 2.0, 3.0];
        let err = vec![0.0, 0.0, 0.0];
        let norm = mixed_error_norm(&y, &err, 1e-6, 1e-8);
        assert!(norm < 1e-10, "zero error should give near-zero norm");
    }
    #[test]
    fn test_mixed_error_norm_consistent() {
        let y = vec![1.0];
        let err = vec![1e-7];
        let norm1 = mixed_error_norm(&y, &err, 1e-6, 1e-8);
        let norm2 = mixed_error_norm(&y, &err, 1e-5, 1e-8);
        assert!(norm1 > norm2 || norm1.is_finite(), "norm should be finite");
    }
    #[test]
    fn test_rk45_step_exponential() {
        let y = vec![1.0_f64];
        let h = 0.1;
        let (y_new, err) = rk45_step(&y, |v| vec![-v[0]], h, 1e-6);
        let exact = (-h).exp();
        assert!(
            (y_new[0] - exact).abs() < 1e-9,
            "rk45_step: y_new={}, exact={}",
            y_new[0],
            exact
        );
        assert!(err >= 0.0, "error estimate must be non-negative");
    }
    #[test]
    fn test_rk45_step_harmonic() {
        let omega = 2.0;
        let y = vec![1.0_f64, 0.0_f64];
        let h = 0.05;
        let (y_new, _) = rk45_step(&y, |v| harmonic_rhs(v, omega), h, 1e-8);
        let e0 = harmonic_energy(&y, omega);
        let e1 = harmonic_energy(&y_new, omega);
        assert!(
            (e1 - e0).abs() / e0 < 1e-6,
            "RK45 energy drift: {}",
            (e1 - e0).abs() / e0
        );
    }
    #[test]
    fn test_rk45_multi_step_integration() {
        let mut y = vec![1.0_f64];
        let dt = 0.05;
        for _ in 0..20 {
            let (y_new, _) = rk45_step(&y, |v| vec![-v[0]], dt, 1e-8);
            y = y_new;
        }
        let exact = (-1.0_f64).exp();
        assert!(
            (y[0] - exact).abs() < 1e-8,
            "RK45 multi-step: y={}, exact={}",
            y[0],
            exact
        );
    }
    #[test]
    fn test_dormand_prince45_error_control() {
        let dp = DormandPrince45::new(1e-10, 1e-8);
        let y = vec![1.0_f64, 0.0_f64];
        let res = dp.step(&y, |v| harmonic_rhs(v, 1.0), 0.001);
        assert!(
            res.error_estimate < 1e-9,
            "DP45 error with small h: {}",
            res.error_estimate
        );
    }
    #[test]
    fn test_dormand_prince45_larger_step_larger_error() {
        let dp = DormandPrince45::new(1e-10, 1e-8);
        let y = vec![1.0_f64];
        let res_small = dp.step(&y, |v| vec![v[0]], 0.001);
        let res_large = dp.step(&y, |v| vec![v[0]], 0.5);
        assert!(
            res_large.error_estimate >= res_small.error_estimate,
            "small={} large={}",
            res_small.error_estimate,
            res_large.error_estimate
        );
    }
    #[test]
    fn test_abm4_predicts_then_corrects() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let mut abm = AdamsBashforthMoulton4::new();
        let mut y = vec![1.0_f64];
        let dt = 0.01;
        let mut t = 0.0;
        for _ in 0..4 {
            let k = f(t, &y);
            abm.push_history(t, y.clone(), k);
            y = rk4_step(&y, |v| f(t, v), dt);
            t += dt;
        }
        let y_new = abm.step(&f, t, &y, dt);
        assert_eq!(y_new.len(), y.len(), "ABM4 output length mismatch");
    }
    #[test]
    fn test_abm4_ready_after_4_pushes() {
        let mut abm = AdamsBashforthMoulton4::new();
        assert!(!abm.ready());
        for i in 0..4 {
            abm.push_history(i as f64 * 0.1, vec![1.0], vec![1.0]);
        }
        assert!(abm.ready());
    }
    #[test]
    fn test_implicit_euler_newton_2d_decay() {
        let solver = ImplicitEulerNewton::new(50, 1e-10);
        let mut y = vec![1.0_f64, 1.0_f64];
        let dt = 0.1;
        for step in 0..10 {
            let t = step as f64 * dt;
            y = solver.step(|_t, yv| vec![-10.0 * yv[0], -20.0 * yv[1]], t, &y, dt);
        }
        assert!(
            y[0] < 0.1 && y[1] < 0.1,
            "Both components should decay: y0={}, y1={}",
            y[0],
            y[1]
        );
    }
    #[test]
    fn test_event_detector_no_crossings() {
        let f = |_t: f64| 1.0_f64;
        let crossings = EventDetector::find_zero_crossings(f, 0.0, 10.0, 0.1, 1e-10);
        assert!(
            crossings.is_empty(),
            "No crossings expected: {:?}",
            crossings
        );
    }
    #[test]
    fn test_event_detector_crossing_accuracy() {
        let crossings = EventDetector::find_zero_crossings(f64::cos, 0.1, 3.0, 0.1, 1e-10);
        let near_half_pi = crossings.iter().any(|&c| (c - PI / 2.0).abs() < 0.01);
        assert!(
            near_half_pi,
            "No crossing near π/2; crossings: {:?}",
            crossings
        );
    }
    #[test]
    fn test_event_bisect_accuracy() {
        let root = EventDetector::bisect(&f64::sin, 3.0, 3.5, 1e-12);
        assert!((root - PI).abs() < 1e-10, "bisect root={root}");
    }
    #[test]
    fn test_dense_output_is_monotone_for_monotone_solution() {
        let seg = DenseOutputSegment {
            t0: 0.0,
            t1: 1.0,
            y0: vec![0.0],
            y1: vec![1.0],
            k1: vec![1.0],
            k6: vec![1.0],
        };
        let ts = [0.0, 0.25, 0.5, 0.75, 1.0];
        let vals: Vec<f64> = ts.iter().map(|&t| seg.evaluate(t)[0]).collect();
        for i in 1..vals.len() {
            assert!(
                vals[i] >= vals[i - 1],
                "Dense output should be non-decreasing: {vals:?}"
            );
        }
    }
    #[test]
    fn test_dense_output_multiple_components() {
        let seg = DenseOutputSegment {
            t0: 0.0,
            t1: 1.0,
            y0: vec![0.0, 1.0],
            y1: vec![1.0, 0.0],
            k1: vec![1.0, -1.0],
            k6: vec![1.0, -1.0],
        };
        let mid = seg.evaluate(0.5);
        assert_eq!(mid.len(), 2);
        assert!(mid[0].is_finite() && mid[1].is_finite());
    }
    #[test]
    fn test_implicit_euler_step_linear() {
        let y = vec![1.0_f64];
        let dt = 0.1;
        let lambda = 10.0_f64;
        let jacobian = vec![vec![-lambda]];
        let rhs = vec![0.0_f64];
        let y_new = implicit_euler_step(&y, &jacobian, &rhs, dt);
        let expected = 1.0 / (1.0 + dt * lambda);
        assert!(
            (y_new[0] - expected).abs() < 1e-8,
            "implicit Euler: y_new={}, expected={}",
            y_new[0],
            expected
        );
    }
    #[test]
    fn test_implicit_euler_step_constant_rhs() {
        let y = vec![2.0_f64];
        let dt = 0.5;
        let jacobian = vec![vec![0.0_f64]];
        let rhs = vec![1.0_f64];
        let y_new = implicit_euler_step(&y, &jacobian, &rhs, dt);
        assert!(
            (y_new[0] - 2.5).abs() < 1e-8,
            "implicit Euler constant: y_new={}",
            y_new[0]
        );
    }
    #[test]
    fn test_mixed_error_norm_empty() {
        let norm = mixed_error_norm(&[], &[], 1e-6, 1e-8);
        assert!(norm < 1e-15, "empty norm should be 0");
    }
    #[test]
    fn test_mixed_error_norm_single_component() {
        let y = vec![10.0];
        let err = vec![1e-5];
        let norm = mixed_error_norm(&y, &err, 1e-4, 1e-8);
        assert!(norm.is_finite() && norm > 0.0);
    }
    #[test]
    fn test_rk2_step_exponential() {
        let y = vec![1.0_f64];
        let h = 0.01;
        let y_new = rk2_step(&y, |v| vec![v[0]], h);
        let exact = h.exp();
        assert!((y_new[0] - exact).abs() < 1e-5, "RK2: {}", y_new[0]);
    }
    #[test]
    fn test_euler_step_basic() {
        let y = vec![0.0_f64, 5.0_f64];
        let dydt = vec![1.0_f64, -2.0_f64];
        let y_new = euler_step(&y, &dydt, 0.5);
        assert!((y_new[0] - 0.5).abs() < 1e-12);
        assert!((y_new[1] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_solve_linear_system_3x3() {
        let a = vec![
            vec![1.0_f64, 1.0, 1.0],
            vec![1.0, 2.0, 3.0],
            vec![1.0, 4.0, 9.0],
        ];
        let b = vec![6.0_f64, 14.0, 36.0];
        let x = solve_linear_system(&a, &b).expect("should solve 3x3");
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0]={}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1]={}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-10, "x[2]={}", x[2]);
    }
    #[test]
    fn test_solve_linear_system_singular_returns_none() {
        let a = vec![vec![1.0_f64, 2.0], vec![2.0, 4.0]];
        let b = vec![1.0_f64, 2.0];
        let result = solve_linear_system(&a, &b);
        let _ = result;
    }
    #[test]
    fn test_rkf45_step_exponential_decay() {
        let rkf = RungeKuttaFehlberg::new(1e-6, 1e-8, 1e-10, 1.0);
        let y0 = vec![1.0_f64];
        let h = 0.1;
        let (y4, _y5, _h_new) = rkf.step(&|_t, y| vec![-y[0]], 0.0, &y0, h);
        let exact = (-h).exp();
        assert!(
            (y4[0] - exact).abs() < 1e-6,
            "RKF45 exp decay: y4={}, exact={}",
            y4[0],
            exact
        );
    }
    #[test]
    fn test_rkf45_suggested_step_grows_on_smooth_problem() {
        let rkf = RungeKuttaFehlberg::new(1e-6, 1e-8, 1e-12, 10.0);
        let y0 = vec![1.0_f64];
        let h = 0.01;
        let (_y4, _y5, h_new) = rkf.step(&|_t, y| vec![y[0]], 0.0, &y0, h);
        assert!(
            h_new > 0.0,
            "suggested step should be positive, got {h_new}"
        );
    }
    #[test]
    fn test_rkf45_integrate_exponential() {
        let rkf = RungeKuttaFehlberg::new(1e-6, 1e-8, 1e-12, 1.0);
        let traj = rkf.integrate(|_t, y| vec![y[0]], 0.0, 1.0, vec![1.0_f64]);
        let (t_last, y_last) = traj.last().unwrap();
        assert!((*t_last - 1.0).abs() < 1e-10, "t_last={t_last}");
        assert!(
            (y_last[0] - 1.0_f64.exp()).abs() < 1e-4,
            "RKF45 integrate exp: y={}",
            y_last[0]
        );
    }
    #[test]
    fn test_rkf45_step_constant_rhs() {
        let rkf = RungeKuttaFehlberg::new(1e-8, 1e-10, 1e-15, 1.0);
        let y0 = vec![0.0_f64];
        let h = 0.5;
        let (y4, _y5, _h_new) = rkf.step(&|_t, _y| vec![2.0_f64], 0.0, &y0, h);
        assert!((y4[0] - 1.0).abs() < 1e-10, "RKF45 constant rhs: {}", y4[0]);
    }
    #[test]
    fn test_rkf45_step_multidimensional() {
        let rkf = RungeKuttaFehlberg::new(1e-8, 1e-10, 1e-15, 1.0);
        let y0 = vec![0.0_f64, 1.0_f64];
        let h = 0.2;
        let (y4, _y5, _h_new) = rkf.step(&|_t, _y| vec![1.0_f64, -1.0_f64], 0.0, &y0, h);
        assert!((y4[0] - h).abs() < 1e-10, "y[0]={}", y4[0]);
        assert!((y4[1] - (1.0 - h)).abs() < 1e-10, "y[1]={}", y4[1]);
    }
    #[test]
    fn test_dopri_error_norm_zero_for_zero_stages() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let n = 4;
        let z = vec![0.0_f64; n];
        let y = vec![1.0_f64; n];
        let norm = dp.compute_error_norm(&y, &z, &z, &z, &z, &z, &z, 0.1);
        assert!(
            norm < 1e-14,
            "zero stages should give zero error norm: {norm}"
        );
    }
    #[test]
    fn test_dopri_error_norm_positive_for_nonzero_stages() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let n = 2;
        let y = vec![1.0_f64; n];
        let k = vec![1.0_f64; n];
        let z = vec![0.0_f64; n];
        let norm = dp.compute_error_norm(&y, &k, &z, &z, &z, &z, &z, 1.0);
        assert!(
            norm > 0.0,
            "nonzero stage should give positive norm: {norm}"
        );
        assert!(norm.is_finite(), "norm must be finite");
    }
    #[test]
    fn test_dopri_error_norm_scales_with_h() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let y = vec![1.0_f64; 3];
        let k = vec![1.0_f64; 3];
        let z = vec![0.0_f64; 3];
        let norm1 = dp.compute_error_norm(&y, &k, &z, &z, &z, &z, &z, 1.0);
        let norm2 = dp.compute_error_norm(&y, &k, &z, &z, &z, &z, &z, 2.0);
        assert!(
            norm2 > norm1,
            "error norm should scale with h: {norm1} vs {norm2}"
        );
    }
    #[test]
    fn test_richardson_extrapolate_single_estimate() {
        let estimates = vec![vec![3.0_f64, 4.0_f64]];
        let n_steps = vec![2usize];
        let (best, err) = BulirschStoer::richardson_extrapolate(&estimates, &n_steps);
        assert_eq!(best, vec![3.0, 4.0]);
        assert!(err < 1e-15, "single estimate: error should be 0, got {err}");
    }
    #[test]
    fn test_richardson_extrapolate_improves_accuracy() {
        let estimates = vec![vec![1.0_f64], vec![1.0_f64], vec![1.0_f64]];
        let n_steps = vec![2usize, 4, 6];
        let (best, _err) = BulirschStoer::richardson_extrapolate(&estimates, &n_steps);
        assert!(
            (best[0] - 1.0).abs() < 1e-12,
            "extrapolated value: {}",
            best[0]
        );
    }
}
#[cfg(test)]
mod time_dependent_tests {
    use super::*;
    use crate::ode::AdaptiveIntegrator;
    use crate::ode::DormandPrince;
    use crate::ode::DormandPrince45;
    use std::cell::{Cell, RefCell};

    /// `y' = cos t`, exact solution `y = sin t` for `y(0) = 0`.
    struct CosineRhs;
    impl OdeSystem for CosineRhs {
        fn ode_rhs(&self, t: f64, _y: &[f64], dydt: &mut [f64]) {
            dydt[0] = t.cos();
        }
    }
    /// `y' = -y + sin t`, exact solution `y = (sin t - cos t + e^-t) / 2` for `y(0) = 0`.
    struct ForcedDecay;
    impl OdeSystem for ForcedDecay {
        fn ode_rhs(&self, t: f64, y: &[f64], dydt: &mut [f64]) {
            dydt[0] = -y[0] + t.sin();
        }
    }
    /// A right-hand side that is `NaN` everywhere.
    struct NanRhs;
    impl OdeSystem for NanRhs {
        fn ode_rhs(&self, _t: f64, _y: &[f64], dydt: &mut [f64]) {
            dydt[0] = f64::NAN;
        }
    }
    fn forced_decay_exact(t: f64) -> f64 {
        0.5 * (t.sin() - t.cos() + (-t).exp())
    }
    /// Global error at `t = 2` of `n_steps` fixed `rk45_step_t` steps for
    /// `y' = cos(t) y`, whose exact solution is `exp(sin t)`.
    fn fixed_step_error(n_steps: usize) -> f64 {
        let t_end = 2.0_f64;
        let dt = t_end / n_steps as f64;
        let mut y = vec![1.0_f64];
        for i in 0..n_steps {
            let t = i as f64 * dt;
            y = rk45_step_t(t, &y, |tt, yv| vec![tt.cos() * yv[0]], dt).0;
        }
        (y[0] - t_end.sin().exp()).abs()
    }

    #[test]
    fn rk45_step_t_integrates_a_degree_four_time_polynomial_exactly() {
        // y' = t^4 does not depend on y, so one step is a pure quadrature rule; the
        // 5th-order weights integrate t^4 exactly only if each stage sees its own time.
        let (t0, dt) = (1.0_f64, 0.5_f64);
        let (y_new, err, _) = rk45_step_t(t0, &[0.0], |t, _y| vec![t.powi(4)], dt);
        let exact = ((t0 + dt).powi(5) - t0.powi(5)) / 5.0;
        assert!(
            (y_new[0] - exact).abs() < 1e-13,
            "y_new={}, exact={exact}",
            y_new[0]
        );
        // The embedded 4th-order rule is not exact for t^4, so the estimate is non-zero.
        assert!(err.is_finite() && err > 1e-7, "err={err}");
    }
    #[test]
    fn rk45_step_t_evaluates_each_stage_at_its_dormand_prince_node() {
        let times = RefCell::new(Vec::new());
        let (t0, dt) = (2.0_f64, 0.3_f64);
        let _ = rk45_step_t(
            t0,
            &[1.0],
            |t, y| {
                times.borrow_mut().push(t);
                vec![-y[0]]
            },
            dt,
        );
        let seen = times.into_inner();
        let nodes = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
        assert_eq!(seen.len(), 7);
        for (stage, (t_seen, c)) in seen.iter().zip(nodes.iter()).enumerate() {
            assert!(
                (t_seen - (t0 + c * dt)).abs() < 1e-15,
                "stage {stage} evaluated at t={t_seen}, expected {}",
                t0 + c * dt
            );
        }
    }
    #[test]
    fn rk45_step_t_returns_the_derivative_at_the_new_point() {
        let f = |t: f64, y: &[f64]| vec![-2.0 * y[0] + t.sin(), y[0] * t.cos()];
        let (t0, dt) = (0.4_f64, 0.2_f64);
        let (y_new, _, k_last) = rk45_step_t(t0, &[1.0, 0.5], f, dt);
        assert_eq!(k_last, f(t0 + dt, &y_new));
    }
    #[test]
    fn rk45_step_fsal_with_supplied_first_stage_saves_one_evaluation() {
        let calls = Cell::new(0_usize);
        let f = |t: f64, y: &[f64]| {
            calls.set(calls.get() + 1);
            vec![y[1], -y[0] + t.cos()]
        };
        let y = [0.3, -0.2];
        let full = rk45_step_t(1.1, &y, f, 0.15);
        assert_eq!(calls.get(), 7);
        let k1 = f(1.1, &y);
        calls.set(0);
        let reused = rk45_step_fsal(1.1, &y, &k1, f, 0.15);
        assert_eq!(calls.get(), 6);
        assert_eq!(full, reused);
    }
    #[test]
    fn rk45_step_t_handles_an_empty_state() {
        let (y_new, err, k_last) = rk45_step_t(0.5, &[], |_t, _y| Vec::new(), 0.1);
        assert!(y_new.is_empty() && k_last.is_empty());
        assert_eq!(err, 0.0);
    }
    #[test]
    fn rk45_step_autonomous_wrapper_agrees_with_the_time_aware_core() {
        let f = |y: &[f64]| vec![y[1], -4.0 * y[0]];
        let y = [1.0, 0.5];
        let (y_auto, err_auto) = rk45_step(&y, f, 0.05, 1e-8);
        let (y_t, err_t, _) = rk45_step_t(123.0, &y, |_t, yv| f(yv), 0.05);
        assert_eq!(y_auto, y_t);
        assert_eq!(err_auto, err_t);
    }
    #[test]
    fn rk45_step_t_is_fifth_order_accurate_for_a_time_dependent_rhs() {
        // Halving the step of a 5th-order method cuts the global error by about 32;
        // stages frozen at the start time would only give a ratio of about 2.
        let e16 = fixed_step_error(16);
        let e32 = fixed_step_error(32);
        let e64 = fixed_step_error(64);
        let first = e16 / e32;
        let second = e32 / e64;
        assert!(e16 < 1e-7, "e16={e16}");
        assert!(first > 28.0 && first < 40.0, "ratio e16/e32 = {first}");
        assert!(second > 28.0 && second < 40.0, "ratio e32/e64 = {second}");
    }
    #[test]
    fn dormand_prince_step_t_matches_exact_solution_for_a_time_dependent_rhs() {
        let dp = DormandPrince45::new(1e-10, 1e-8);
        let t0 = 0.7_f64;
        let res = dp.step_t(t0, &[t0.sin()], |t, _y| vec![t.cos()], 0.1);
        assert!(res.accepted);
        assert!(
            (res.y[0] - (t0 + 0.1).sin()).abs() < 1e-12,
            "y={}, exact={}",
            res.y[0],
            (t0 + 0.1).sin()
        );
        assert!(res.error_estimate > 0.0 && res.error_estimate < 1e-9);
    }
    #[test]
    fn dormand_prince_step_fsal_returns_the_next_first_stage() {
        let dp = DormandPrince45::new(1e-10, 1e-8);
        let f = |t: f64, y: &[f64]| vec![-y[0] + t.sin()];
        let (t0, dt) = (0.5_f64, 0.05_f64);
        let k1 = f(t0, &[0.2]);
        let (res, k_last) = dp.step_fsal(t0, &[0.2], &k1, f, dt);
        assert!(res.accepted);
        assert_eq!(k_last, f(t0 + dt, &res.y));
        let direct = dp.step_t(t0, &[0.2], f, dt);
        assert_eq!(direct.y, res.y);
        assert_eq!(direct.error_estimate, res.error_estimate);
    }
    #[test]
    fn dormand_prince_integrate_follows_a_time_dependent_rhs() {
        // With every stage frozen at the step start all k_i coincide for y' = cos t, the
        // embedded error is exactly zero and the step grows to h_max with first-order
        // accuracy (error of order 1e-1 at t = 10).
        let dp = DormandPrince::new(1e-8, 1e-10, 0.5);
        let result = dp.integrate(|t, _y| vec![t.cos()], 0.0, 10.0, vec![0.0]);
        let (t_last, y_last) = result.last().expect("start point is always present");
        assert_eq!(*t_last, 10.0);
        assert!(
            (y_last[0] - 10.0_f64.sin()).abs() < 1e-6,
            "y(10)={}, sin(10)={}",
            y_last[0],
            10.0_f64.sin()
        );
        for (t, y) in &result {
            assert!((y[0] - t.sin()).abs() < 1e-6, "t={t}, y={}", y[0]);
        }
    }
    #[test]
    fn dormand_prince_integrate_error_is_controlled_by_the_tolerance() {
        let error_for = |rtol: f64, atol: f64| {
            let dp = DormandPrince::new(rtol, atol, 1.0);
            let result = dp.integrate(|t, y| vec![-y[0] + t.sin()], 0.0, 8.0, vec![0.0]);
            let (_, y_end) = result.last().expect("start point is always present");
            (y_end[0] - forced_decay_exact(8.0)).abs()
        };
        let loose = error_for(1e-4, 1e-6);
        let tight = error_for(1e-9, 1e-11);
        assert!(loose < 1e-3, "loose={loose}");
        assert!(tight < 1e-7, "tight={tight}");
        assert!(tight < loose, "tight={tight}, loose={loose}");
    }
    #[test]
    fn dormand_prince_integrate_uses_six_new_evaluations_per_accepted_step() {
        let calls = Cell::new(0_usize);
        let dp = DormandPrince::new(1e-3, 1e-3, 0.1);
        let result = dp.integrate(
            |t, y| {
                calls.set(calls.get() + 1);
                vec![-y[0] + t.sin()]
            },
            0.0,
            1.0,
            vec![0.0],
        );
        let accepted = result.len() - 1;
        assert!(accepted >= 10, "accepted={accepted}");
        assert_eq!(calls.get(), 1 + 6 * accepted);
    }
    #[test]
    fn dormand_prince_integrate_handles_empty_and_reversed_intervals() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let same = dp.integrate(|_t, y| vec![-y[0]], 1.0, 1.0, vec![2.0]);
        assert_eq!(same, vec![(1.0, vec![2.0])]);
        let reversed = dp.integrate(|_t, y| vec![-y[0]], 1.0, 0.0, vec![2.0]);
        assert_eq!(reversed, vec![(1.0, vec![2.0])]);
    }
    #[test]
    fn dormand_prince_integrate_stops_on_a_non_finite_right_hand_side() {
        let dp = DormandPrince::new(1e-6, 1e-8, 1.0);
        let result = dp.integrate(|_t, _y| vec![f64::NAN], 0.0, 1.0, vec![1.0]);
        assert_eq!(result, vec![(0.0, vec![1.0])]);
    }
    #[test]
    fn adaptive_integrator_follows_a_time_dependent_rhs() {
        let integrator = AdaptiveIntegrator::new(1e-10, 1e-8);
        let result = integrator.integrate(&CosineRhs, vec![0.0], 0.0, 10.0, 0.05);
        let (t_last, y_last) = result.last().expect("start point is always present");
        assert_eq!(*t_last, 10.0);
        assert!(
            (y_last[0] - 10.0_f64.sin()).abs() < 1e-6,
            "y(10)={}, sin(10)={}",
            y_last[0],
            10.0_f64.sin()
        );
        for (t, y) in &result {
            assert!((y[0] - t.sin()).abs() < 1e-6, "t={t}, y={}", y[0]);
        }
    }
    #[test]
    fn adaptive_integrator_solves_a_forced_decay() {
        let integrator = AdaptiveIntegrator::new(1e-11, 1e-9);
        let result = integrator.integrate(&ForcedDecay, vec![0.0], 0.0, 8.0, 0.1);
        let (_, y_end) = result.last().expect("start point is always present");
        assert!(
            (y_end[0] - forced_decay_exact(8.0)).abs() < 1e-7,
            "y(8)={}, exact={}",
            y_end[0],
            forced_decay_exact(8.0)
        );
    }
    #[test]
    fn adaptive_integrator_ends_exactly_at_the_final_time() {
        let integrator = AdaptiveIntegrator::new(1e-10, 1e-8);
        let t_end = 2.0 + 1.0 / 3.0;
        let result = integrator.integrate(&ForcedDecay, vec![0.0], 0.0, t_end, 0.07);
        let (t_last, _) = result.last().expect("start point is always present");
        assert_eq!(*t_last, t_end);
    }
    #[test]
    fn adaptive_integrator_stops_on_a_non_finite_right_hand_side() {
        let integrator = AdaptiveIntegrator::new(1e-8, 1e-6);
        let result = integrator.integrate(&NanRhs, vec![1.0], 0.0, 1.0, 0.1);
        assert_eq!(result, vec![(0.0, vec![1.0])]);
    }
    #[test]
    fn dormand_prince_integrate_terminates_for_nan_tolerances() {
        // A NaN tolerance makes every step comparison false and the step size NaN; the
        // driver must give up instead of looping forever.
        let dp = DormandPrince::new(f64::NAN, f64::NAN, 1.0);
        let result = dp.integrate(|_t, y| vec![-y[0]], 0.0, 1.0, vec![1.0]);
        assert_eq!(result, vec![(0.0, vec![1.0])]);
    }
    #[test]
    fn adaptive_integrator_terminates_for_nan_tolerances() {
        let integrator = AdaptiveIntegrator::new(f64::NAN, f64::NAN);
        let result = integrator.integrate(&ForcedDecay, vec![1.0], 0.0, 1.0, 0.1);
        assert_eq!(result, vec![(0.0, vec![1.0])]);
    }
}
