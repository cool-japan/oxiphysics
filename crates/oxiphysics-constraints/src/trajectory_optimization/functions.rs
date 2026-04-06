//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use super::types::{
    BSplineTrajectory, BoundaryConditions, CostFunction, DynamicsModel, KeepOutZone,
    ShootingResult, SqpStepResult, Trajectory, TrajectoryKnot, TrustRegionConfig,
    TrustRegionResult,
};

/// Dot product of two equal-length slices.
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Euclidean norm of a slice.
pub(super) fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}
/// Element-wise a + b.
pub(super) fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Element-wise a - b.
pub(super) fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Scalar multiply.
pub(super) fn vec_scale(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}
/// Linear interpolation between two vectors.
pub(super) fn vec_lerp(a: &[f64], b: &[f64], t: f64) -> Vec<f64> {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| x + (y - x) * t)
        .collect()
}
/// Multiply an n x n row-major matrix by a vector of length n.
pub(super) fn mat_vec(mat: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += mat[i * n + j] * x[j];
        }
    }
    y
}
/// Multiply two n x n row-major matrices.
pub(super) fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}
/// Transpose an n x n row-major matrix.
pub(super) fn mat_transpose(a: &[f64], n: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[j * n + i] = a[i * n + j];
        }
    }
    t
}
/// Identity n x n matrix (row-major).
pub(super) fn mat_eye(n: usize) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}
/// Add two n x n matrices.
pub(super) fn mat_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Subtract two n x n matrices.
pub(super) fn mat_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Scale a matrix.
pub(super) fn mat_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}
/// Solve a small n x n linear system Ax = b via Gaussian elimination with
/// partial pivoting. Returns `None` if singular.
pub(super) fn solve_linear(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }
    for col in 0..n {
        let mut pivot = col;
        let mut best = aug[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let v = aug[row * (n + 1) + col].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1e-15 {
            return None;
        }
        if pivot != col {
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, pivot * (n + 1) + j);
            }
        }
        let diag = aug[col * (n + 1) + col];
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / diag;
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            s -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = s / aug[i * (n + 1) + i];
    }
    Some(x)
}
/// Perform single-shooting forward simulation with given controls.
///
/// Integrates the dynamics from `x0` over `num_steps` steps using RK4.
pub fn direct_shooting(
    model: &DynamicsModel,
    x0: &[f64],
    controls: &[Vec<f64>],
    t0: f64,
    dt: f64,
    cost: &CostFunction,
) -> ShootingResult {
    let num_steps = controls.len();
    let mut traj = Trajectory::new(model.n_state, model.n_control);
    let mut x = x0.to_vec();
    let mut t = t0;
    let mut total_cost = 0.0;
    traj.knots.push(TrajectoryKnot {
        time: t,
        state: x.clone(),
        control: if controls.is_empty() {
            vec![0.0; model.n_control]
        } else {
            controls[0].clone()
        },
    });
    for i in 0..num_steps {
        let u = &controls[i];
        total_cost += cost.running_cost(&x, u) * dt;
        x = model.rk4_step(&x, u, t, dt);
        t += dt;
        traj.knots.push(TrajectoryKnot {
            time: t,
            state: x.clone(),
            control: if i + 1 < num_steps {
                controls[i + 1].clone()
            } else {
                vec![0.0; model.n_control]
            },
        });
    }
    total_cost += cost.terminal_cost(&x);
    ShootingResult {
        trajectory: traj,
        final_state: x,
        cost: total_cost,
        final_error: 0.0,
    }
}
/// Perform multiple shooting with `num_segments` segments.
///
/// Returns defects at segment boundaries and the trajectory.
pub fn multiple_shooting(
    model: &DynamicsModel,
    segment_states: &[Vec<f64>],
    controls: &[Vec<f64>],
    t0: f64,
    dt: f64,
    steps_per_segment: usize,
) -> (Vec<Vec<f64>>, Trajectory) {
    let num_segments = segment_states.len() - 1;
    let mut defects = Vec::with_capacity(num_segments);
    let mut traj = Trajectory::new(model.n_state, model.n_control);
    for seg in 0..num_segments {
        let mut x = segment_states[seg].clone();
        let mut t = t0 + (seg * steps_per_segment) as f64 * dt;
        for step in 0..steps_per_segment {
            let ctrl_idx = seg * steps_per_segment + step;
            let u = if ctrl_idx < controls.len() {
                &controls[ctrl_idx]
            } else {
                &controls[controls.len() - 1]
            };
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: x.clone(),
                control: u.clone(),
            });
            x = model.rk4_step(&x, u, t, dt);
            t += dt;
        }
        let defect = vec_sub(&x, &segment_states[seg + 1]);
        defects.push(defect);
    }
    let final_t = t0 + (num_segments * steps_per_segment) as f64 * dt;
    traj.knots.push(TrajectoryKnot {
        time: final_t,
        state: segment_states[num_segments].clone(),
        control: vec![0.0; model.n_control],
    });
    (defects, traj)
}
/// Defect constraints for trapezoidal transcription.
///
/// For each segment \[k, k+1\], the defect is:
///   d_k = x_{k+1} - x_k - (h/2)*(f_k + f_{k+1})
///
/// Returns a vector of defects (one per segment, each of length n_state).
pub fn trapezoidal_defects(model: &DynamicsModel, trajectory: &Trajectory) -> Vec<Vec<f64>> {
    let n = trajectory.knots.len();
    if n < 2 {
        return Vec::new();
    }
    let mut defects = Vec::with_capacity(n - 1);
    for k in 0..n - 1 {
        let xk = &trajectory.knots[k].state;
        let uk = &trajectory.knots[k].control;
        let tk = trajectory.knots[k].time;
        let xk1 = &trajectory.knots[k + 1].state;
        let uk1 = &trajectory.knots[k + 1].control;
        let tk1 = trajectory.knots[k + 1].time;
        let h = tk1 - tk;
        let fk = model.eval(xk, uk, tk);
        let fk1 = model.eval(xk1, uk1, tk1);
        let defect: Vec<f64> = (0..model.n_state)
            .map(|i| xk1[i] - xk[i] - 0.5 * h * (fk[i] + fk1[i]))
            .collect();
        defects.push(defect);
    }
    defects
}
/// Maximum defect norm across all segments.
pub fn max_defect_norm(defects: &[Vec<f64>]) -> f64 {
    defects.iter().map(|d| norm(d)).fold(0.0_f64, f64::max)
}
/// Hermite-Simpson collocation defects.
///
/// For each segment \[k, k+1\], uses a cubic Hermite interpolant with midpoint
/// collocation:
///
///   x_c = 0.5*(x_k + x_{k+1}) + (h/8)*(f_k - f_{k+1})
///   d_k = x_{k+1} - x_k - (h/6)*(f_k + 4*f_c + f_{k+1})
///
/// where f_c = f(x_c, u_c, t_c) and u_c = 0.5*(u_k + u_{k+1}).
pub fn hermite_simpson_defects(model: &DynamicsModel, trajectory: &Trajectory) -> Vec<Vec<f64>> {
    let n = trajectory.knots.len();
    if n < 2 {
        return Vec::new();
    }
    let mut defects = Vec::with_capacity(n - 1);
    for k in 0..n - 1 {
        let xk = &trajectory.knots[k].state;
        let uk = &trajectory.knots[k].control;
        let tk = trajectory.knots[k].time;
        let xk1 = &trajectory.knots[k + 1].state;
        let uk1 = &trajectory.knots[k + 1].control;
        let tk1 = trajectory.knots[k + 1].time;
        let h = tk1 - tk;
        let tc = 0.5 * (tk + tk1);
        let fk = model.eval(xk, uk, tk);
        let fk1 = model.eval(xk1, uk1, tk1);
        let xc: Vec<f64> = (0..model.n_state)
            .map(|i| 0.5 * (xk[i] + xk1[i]) + (h / 8.0) * (fk[i] - fk1[i]))
            .collect();
        let uc: Vec<f64> = (0..model.n_control)
            .map(|i| 0.5 * (uk[i] + uk1[i]))
            .collect();
        let fc = model.eval(&xc, &uc, tc);
        let defect: Vec<f64> = (0..model.n_state)
            .map(|i| xk1[i] - xk[i] - (h / 6.0) * (fk[i] + 4.0 * fc[i] + fk1[i]))
            .collect();
        defects.push(defect);
    }
    defects
}
/// Backward sweep costate estimation using the adjoint method.
///
/// Given a trajectory, cost function, and dynamics model, computes the
/// costate (lambda) at each knot by integrating backward:
///
///   d(lambda)/dt = -(df/dx)^T * lambda - dL/dx
///
/// Terminal condition: lambda(T) = dPhi/dx(x_T) = 2*Qf*x_T.
pub fn adjoint_costate(
    model: &DynamicsModel,
    trajectory: &Trajectory,
    cost: &CostFunction,
) -> Vec<Vec<f64>> {
    let n_knots = trajectory.knots.len();
    let ns = model.n_state;
    if n_knots == 0 {
        return Vec::new();
    }
    let mut costates = vec![vec![0.0; ns]; n_knots];
    let xf = &trajectory.knots[n_knots - 1].state;
    let qf_xf = mat_vec(&cost.qf_matrix, xf, ns);
    costates[n_knots - 1] = vec_scale(&qf_xf, 2.0);
    for k in (0..n_knots - 1).rev() {
        let xk = &trajectory.knots[k].state;
        let uk = &trajectory.knots[k].control;
        let tk = trajectory.knots[k].time;
        let dt = trajectory.knots[k + 1].time - tk;
        let a_mat = model.state_jacobian(xk, uk, tk, 1e-7);
        let q_x = mat_vec(&cost.q_matrix, xk, ns);
        let dl_dx = vec_scale(&q_x, 2.0);
        let at = mat_transpose(&a_mat, ns);
        let at_lam = mat_vec(&at, &costates[k + 1], ns);
        let rhs = vec_add(&at_lam, &dl_dx);
        costates[k] = vec_add(&costates[k + 1], &vec_scale(&rhs, dt));
    }
    costates
}
/// Perform a single SQP-like gradient descent step on the controls.
///
/// Uses finite-difference gradient of the augmented Lagrangian, then takes
/// a step in the negative gradient direction.
pub fn sqp_step(
    model: &DynamicsModel,
    trajectory: &Trajectory,
    cost: &CostFunction,
    boundary: &BoundaryConditions,
    step_size: f64,
    penalty: f64,
    tol: f64,
) -> SqpStepResult {
    let n_knots = trajectory.knots.len();
    let _ns = model.n_state;
    let nc = model.n_control;
    let base_cost = trajectory.evaluate_cost(cost);
    let defects = trapezoidal_defects(model, trajectory);
    let max_defect = max_defect_norm(&defects);
    let final_err = boundary.final_violation(
        &trajectory
            .knots
            .last()
            .map(|k| k.state.clone())
            .unwrap_or_default(),
    );
    let augmented_cost = base_cost + penalty * (max_defect + final_err);
    let eps = 1e-5;
    let mut new_traj = trajectory.clone();
    for k in 0..n_knots.saturating_sub(1) {
        for j in 0..nc {
            new_traj.knots[k].control[j] += eps;
            let perturbed_cost = {
                let tc = new_traj.evaluate_cost(cost);
                let defs = trapezoidal_defects(model, &new_traj);
                let md = max_defect_norm(&defs);
                let fe = boundary.final_violation(
                    &new_traj
                        .knots
                        .last()
                        .map(|kk| kk.state.clone())
                        .unwrap_or_default(),
                );
                tc + penalty * (md + fe)
            };
            let grad = (perturbed_cost - augmented_cost) / eps;
            new_traj.knots[k].control[j] = trajectory.knots[k].control[j] - step_size * grad;
        }
    }
    let controls: Vec<Vec<f64>> = new_traj
        .knots
        .iter()
        .take(n_knots.saturating_sub(1))
        .map(|k| k.control.clone())
        .collect();
    let x0 = trajectory.knots[0].state.clone();
    let t0 = trajectory.knots[0].time;
    let dt = if n_knots > 1 {
        trajectory.knots[1].time - trajectory.knots[0].time
    } else {
        0.01
    };
    let result = direct_shooting(model, &x0, &controls, t0, dt, cost);
    let new_defects = trapezoidal_defects(model, &result.trajectory);
    let new_max_defect = max_defect_norm(&new_defects);
    let new_final_err = boundary.final_violation(&result.final_state);
    let new_violation = new_max_defect + new_final_err;
    let converged = new_violation < tol && (result.cost - base_cost).abs() < tol;
    SqpStepResult {
        trajectory: result.trajectory,
        cost: result.cost,
        max_violation: new_violation,
        step_size,
        converged,
    }
}
/// Evaluate a trust region step.
///
/// Given a proposed step `delta_u` (control perturbation), evaluates the
/// actual vs predicted cost reduction and updates the trust region radius.
pub fn trust_region_evaluate(
    config: &TrustRegionConfig,
    current_radius: f64,
    current_cost: f64,
    proposed_cost: f64,
    predicted_reduction: f64,
    step_norm: f64,
) -> TrustRegionResult {
    let actual_reduction = current_cost - proposed_cost;
    let ratio = if predicted_reduction.abs() < 1e-15 {
        0.0
    } else {
        actual_reduction / predicted_reduction
    };
    let accepted = ratio > config.accept_ratio && step_norm <= current_radius * 1.01;
    let new_radius = if ratio > 0.75 && step_norm > 0.9 * current_radius {
        (current_radius * config.expand_factor).min(config.max_radius)
    } else if ratio < 0.25 {
        (current_radius * config.shrink_factor).max(config.min_radius)
    } else {
        current_radius
    };
    TrustRegionResult {
        accepted,
        new_radius,
        actual_reduction,
        predicted_reduction,
        ratio,
    }
}
/// Uniform cubic B-spline basis functions.
///
/// Returns the 4 basis function values at parameter `t` in `[0, 1)`.
pub(super) fn bspline_basis(t: f64) -> [f64; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    let b0 = (1.0 - 3.0 * t + 3.0 * t2 - t3) / 6.0;
    let b1 = (4.0 - 6.0 * t2 + 3.0 * t3) / 6.0;
    let b2 = (1.0 + 3.0 * t + 3.0 * t2 - 3.0 * t3) / 6.0;
    let b3 = t3 / 6.0;
    [b0, b1, b2, b3]
}
/// Derivative of uniform cubic B-spline basis functions at parameter `t`.
pub(super) fn bspline_basis_derivative(t: f64) -> [f64; 4] {
    let t2 = t * t;
    let db0 = (-3.0 + 6.0 * t - 3.0 * t2) / 6.0;
    let db1 = (-12.0 * t + 9.0 * t2) / 6.0;
    let db2 = (3.0 + 6.0 * t - 9.0 * t2) / 6.0;
    let db3 = 3.0 * t2 / 6.0;
    [db0, db1, db2, db3]
}
/// Smooth a trajectory using B-spline fitting.
///
/// Extracts positions from the trajectory, fits a B-spline, then re-samples.
pub fn smooth_trajectory(
    trajectory: &Trajectory,
    num_control_points: usize,
    num_output_knots: usize,
) -> Trajectory {
    let waypoints: Vec<Vec<f64>> = trajectory.knots.iter().map(|k| k.state.clone()).collect();
    let spline = BSplineTrajectory::fit_to_waypoints(&waypoints, num_control_points);
    let samples = spline.sample(num_output_knots);
    let t0 = trajectory.knots.first().map(|k| k.time).unwrap_or(0.0);
    let tf = trajectory.knots.last().map(|k| k.time).unwrap_or(1.0);
    let mut new_traj = Trajectory::new(trajectory.n_state, trajectory.n_control);
    for (i, state) in samples.into_iter().enumerate() {
        let alpha = if num_output_knots > 1 {
            i as f64 / (num_output_knots - 1) as f64
        } else {
            0.0
        };
        let t = t0 + alpha * (tf - t0);
        new_traj.knots.push(TrajectoryKnot {
            time: t,
            state,
            control: vec![0.0; trajectory.n_control],
        });
    }
    new_traj
}
/// Scale a trajectory's time grid by a factor (for free final time problems).
///
/// Returns a new trajectory with all times multiplied by `scale_factor`.
pub fn scale_time(trajectory: &Trajectory, scale_factor: f64) -> Trajectory {
    let mut traj = trajectory.clone();
    for knot in &mut traj.knots {
        knot.time *= scale_factor;
    }
    traj
}
/// Gradient of cost with respect to final time (for free final time).
///
/// Approximates dJ/dT by perturbing the time scale.
pub fn free_time_gradient(
    model: &DynamicsModel,
    trajectory: &Trajectory,
    cost: &CostFunction,
    eps: f64,
) -> f64 {
    let duration = trajectory.duration();
    if duration.abs() < 1e-15 {
        return 0.0;
    }
    let scale_plus = (duration + eps) / duration;
    let scale_minus = (duration - eps) / duration;
    let traj_plus = scale_time(trajectory, scale_plus);
    let traj_minus = scale_time(trajectory, scale_minus);
    let cost_plus = traj_plus.evaluate_cost(cost);
    let cost_minus = traj_minus.evaluate_cost(cost);
    let def_plus = max_defect_norm(&trapezoidal_defects(model, &traj_plus));
    let def_minus = max_defect_norm(&trapezoidal_defects(model, &traj_minus));
    let aug_plus = cost_plus + 100.0 * def_plus;
    let aug_minus = cost_minus + 100.0 * def_minus;
    (aug_plus - aug_minus) / (2.0 * eps)
}
/// Compute the total obstacle avoidance penalty for a trajectory.
///
/// Uses a smooth penalty: max(0, radius - distance)^2 for each knot.
pub fn obstacle_penalty(trajectory: &Trajectory, obstacles: &[KeepOutZone]) -> f64 {
    let mut total = 0.0;
    for knot in &trajectory.knots {
        for obs in obstacles {
            let v = obs.violation(&knot.state);
            if v < 0.0 {
                total += v * v;
            }
        }
    }
    total
}
/// Check if the trajectory is obstacle-free.
pub fn is_obstacle_free(trajectory: &Trajectory, obstacles: &[KeepOutZone]) -> bool {
    for knot in &trajectory.knots {
        for obs in obstacles {
            if obs.violation(&knot.state) < 0.0 {
                return false;
            }
        }
    }
    true
}
/// Compute the jerk (third derivative) cost of a trajectory.
///
/// Uses finite differences on the acceleration (which itself is finite-
/// differenced from velocity).
pub fn jerk_cost(trajectory: &Trajectory) -> f64 {
    let n_knots = trajectory.knots.len();
    if n_knots < 4 {
        return 0.0;
    }
    let ns = trajectory.n_state;
    let mut total = 0.0;
    let mut accels = Vec::with_capacity(n_knots - 2);
    for k in 1..n_knots - 1 {
        let dt1 = trajectory.knots[k].time - trajectory.knots[k - 1].time;
        let dt2 = trajectory.knots[k + 1].time - trajectory.knots[k].time;
        let dt_avg = 0.5 * (dt1 + dt2);
        if dt_avg.abs() < 1e-15 {
            accels.push(vec![0.0; ns]);
            continue;
        }
        let acc: Vec<f64> = (0..ns)
            .map(|i| {
                let v1 = (trajectory.knots[k].state[i] - trajectory.knots[k - 1].state[i]) / dt1;
                let v2 = (trajectory.knots[k + 1].state[i] - trajectory.knots[k].state[i]) / dt2;
                (v2 - v1) / dt_avg
            })
            .collect();
        accels.push(acc);
    }
    for k in 1..accels.len() {
        let dt = trajectory.knots[k + 1].time - trajectory.knots[k].time;
        if dt.abs() < 1e-15 {
            continue;
        }
        let jerk_sq: f64 = (0..ns)
            .map(|i| {
                let j = (accels[k][i] - accels[k - 1][i]) / dt;
                j * j
            })
            .sum();
        total += jerk_sq * dt;
    }
    total
}
/// Compute the gradient of total cost w.r.t. all control variables.
///
/// Controls are flattened into a single vector; the gradient is returned in
/// the same layout.
pub fn cost_gradient_controls(
    model: &DynamicsModel,
    x0: &[f64],
    controls_flat: &[f64],
    t0: f64,
    dt: f64,
    n_knots: usize,
    cost: &CostFunction,
    eps: f64,
) -> Vec<f64> {
    let nc = model.n_control;
    let num_controls = (n_knots - 1) * nc;
    let base_controls: Vec<Vec<f64>> = controls_flat.chunks(nc).map(|c| c.to_vec()).collect();
    let base_result = direct_shooting(model, x0, &base_controls, t0, dt, cost);
    let base_cost = base_result.cost;
    let mut grad = vec![0.0; num_controls];
    let mut perturbed = controls_flat.to_vec();
    for i in 0..num_controls {
        perturbed[i] += eps;
        let pert_controls: Vec<Vec<f64>> = perturbed.chunks(nc).map(|c| c.to_vec()).collect();
        let pert_result = direct_shooting(model, x0, &pert_controls, t0, dt, cost);
        grad[i] = (pert_result.cost - base_cost) / eps;
        perturbed[i] = controls_flat[i];
    }
    grad
}
