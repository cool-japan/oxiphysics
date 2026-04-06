//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Velocity and acceleration bounds.
#[derive(Debug, Clone)]
pub struct DynamicLimits {
    /// Maximum velocity per state component (absolute value).
    pub max_velocity: Vec<f64>,
    /// Maximum acceleration per state component (absolute value).
    pub max_acceleration: Vec<f64>,
    /// Maximum control per control component (absolute value).
    pub max_control: Vec<f64>,
}
impl DynamicLimits {
    /// Create uniform limits.
    pub fn uniform(n_state: usize, n_control: usize, v_max: f64, a_max: f64, u_max: f64) -> Self {
        Self {
            max_velocity: vec![v_max; n_state],
            max_acceleration: vec![a_max; n_state],
            max_control: vec![u_max; n_control],
        }
    }
    /// Check if a control vector satisfies the limits.
    pub fn control_feasible(&self, u: &[f64]) -> bool {
        for (i, &ui) in u.iter().enumerate() {
            if i < self.max_control.len() && ui.abs() > self.max_control[i] {
                return false;
            }
        }
        true
    }
    /// Clamp a control vector to feasible bounds.
    pub fn clamp_control(&self, u: &[f64]) -> Vec<f64> {
        u.iter()
            .enumerate()
            .map(|(i, &ui)| {
                if i < self.max_control.len() {
                    ui.clamp(-self.max_control[i], self.max_control[i])
                } else {
                    ui
                }
            })
            .collect()
    }
    /// Check if velocity components (extracted from state) are within bounds.
    /// Assumes velocity is in the second half of the state vector for a
    /// position-velocity system.
    pub fn velocity_feasible(&self, state: &[f64]) -> bool {
        let n_pos = state.len() / 2;
        for i in 0..n_pos.min(self.max_velocity.len()) {
            if state[n_pos + i].abs() > self.max_velocity[i] {
                return false;
            }
        }
        true
    }
}
/// Cost function configuration with weighting matrices.
#[derive(Debug, Clone)]
pub struct CostFunction {
    /// The type of cost.
    pub cost_type: CostType,
    /// State weighting matrix Q (n_state x n_state, row-major).
    pub q_matrix: Vec<f64>,
    /// Control weighting matrix R (n_control x n_control, row-major).
    pub r_matrix: Vec<f64>,
    /// Terminal cost matrix Qf (n_state x n_state, row-major).
    pub qf_matrix: Vec<f64>,
    /// State dimension.
    pub n_state: usize,
    /// Control dimension.
    pub n_control: usize,
}
impl CostFunction {
    /// Create a quadratic cost with diagonal Q and R matrices.
    pub fn quadratic(n_state: usize, n_control: usize, q_diag: &[f64], r_diag: &[f64]) -> Self {
        let mut q = vec![0.0; n_state * n_state];
        for i in 0..n_state.min(q_diag.len()) {
            q[i * n_state + i] = q_diag[i];
        }
        let mut r = vec![0.0; n_control * n_control];
        for i in 0..n_control.min(r_diag.len()) {
            r[i * n_control + i] = r_diag[i];
        }
        let qf = q.clone();
        Self {
            cost_type: CostType::Quadratic,
            q_matrix: q,
            r_matrix: r,
            qf_matrix: qf,
            n_state,
            n_control,
        }
    }
    /// Create a minimum-energy cost (Q=0, R=I).
    pub fn minimum_energy(n_state: usize, n_control: usize) -> Self {
        let q = vec![0.0; n_state * n_state];
        let mut r = vec![0.0; n_control * n_control];
        for i in 0..n_control {
            r[i * n_control + i] = 1.0;
        }
        Self {
            cost_type: CostType::MinimumEnergy,
            q_matrix: q.clone(),
            r_matrix: r,
            qf_matrix: q,
            n_state,
            n_control,
        }
    }
    /// Create a minimum-time cost.
    pub fn minimum_time(n_state: usize, n_control: usize) -> Self {
        Self {
            cost_type: CostType::MinimumTime,
            q_matrix: vec![0.0; n_state * n_state],
            r_matrix: vec![0.0; n_control * n_control],
            qf_matrix: vec![0.0; n_state * n_state],
            n_state,
            n_control,
        }
    }
    /// Evaluate running cost at a single knot: x^T Q x + u^T R u.
    pub fn running_cost(&self, x: &[f64], u: &[f64]) -> f64 {
        match self.cost_type {
            CostType::MinimumTime => 1.0,
            CostType::MinimumFuel => u.iter().map(|v| v.abs()).sum(),
            _ => {
                let qx = mat_vec(&self.q_matrix, x, self.n_state);
                let ru = mat_vec(&self.r_matrix, u, self.n_control);
                dot(x, &qx) + dot(u, &ru)
            }
        }
    }
    /// Evaluate terminal cost: x^T Qf x.
    pub fn terminal_cost(&self, x: &[f64]) -> f64 {
        let qfx = mat_vec(&self.qf_matrix, x, self.n_state);
        dot(x, &qfx)
    }
}
/// A spherical keep-out zone for obstacle avoidance.
#[derive(Debug, Clone)]
pub struct KeepOutZone {
    /// Center of the obstacle `[x, y, z]`.
    pub center: [f64; 3],
    /// Radius of the keep-out zone.
    pub radius: f64,
}
impl KeepOutZone {
    /// Create a new keep-out zone.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
    /// Check if a 3D position (taken from first 3 state components) violates
    /// this keep-out zone. Returns the violation (negative = inside).
    pub fn violation(&self, pos: &[f64]) -> f64 {
        if pos.len() < 3 {
            return 0.0;
        }
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        let dz = pos[2] - self.center[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist - self.radius
    }
    /// Gradient of the violation function w.r.t. the first 3 state components.
    pub fn violation_gradient(&self, pos: &[f64]) -> [f64; 3] {
        if pos.len() < 3 {
            return [0.0; 3];
        }
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        let dz = pos[2] - self.center[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist < 1e-15 {
            return [0.0; 3];
        }
        [dx / dist, dy / dist, dz / dist]
    }
}
/// Convergence history entry.
#[derive(Debug, Clone)]
pub struct ConvergenceEntry {
    /// Iteration number.
    pub iteration: usize,
    /// Cost at this iteration.
    pub cost: f64,
    /// Maximum constraint violation.
    pub max_violation: f64,
    /// Step size used.
    pub step_size: f64,
    /// Cost reduction from previous iteration.
    pub cost_reduction: f64,
}
/// Result of a single SQP iteration step.
#[derive(Debug, Clone)]
pub struct SqpStepResult {
    /// Updated trajectory.
    pub trajectory: Trajectory,
    /// Current cost.
    pub cost: f64,
    /// Maximum constraint violation.
    pub max_violation: f64,
    /// Step size used.
    pub step_size: f64,
    /// Converged flag.
    pub converged: bool,
}
/// Waypoint constraint: the trajectory must pass through a given state at a
/// given index (or nearest knot).
#[derive(Debug, Clone)]
pub struct WaypointConstraint {
    /// Knot index where the waypoint applies.
    pub knot_index: usize,
    /// Target state values (indexed by state component).
    pub target: Vec<f64>,
    /// Which state components are constrained (bitmask indices).
    pub mask: Vec<bool>,
    /// Tolerance for each component.
    pub tolerance: f64,
}
impl WaypointConstraint {
    /// Create a full-state waypoint constraint.
    pub fn full_state(knot_index: usize, target: Vec<f64>, tolerance: f64) -> Self {
        let mask = vec![true; target.len()];
        Self {
            knot_index,
            target,
            mask,
            tolerance,
        }
    }
    /// Create a position-only waypoint (first 3 components).
    pub fn position_only(knot_index: usize, pos: [f64; 3], n_state: usize, tolerance: f64) -> Self {
        let mut target = vec![0.0; n_state];
        target[0] = pos[0];
        target[1] = pos[1];
        target[2] = pos[2];
        let mut mask = vec![false; n_state];
        mask[0] = true;
        mask[1] = true;
        mask[2] = true;
        Self {
            knot_index,
            target,
            mask,
            tolerance,
        }
    }
    /// Evaluate waypoint violation (L2 norm of masked error).
    pub fn violation(&self, state: &[f64]) -> f64 {
        let mut sum_sq = 0.0;
        for (i, (&m, (&t, &s))) in self
            .mask
            .iter()
            .zip(self.target.iter().zip(state.iter()))
            .enumerate()
        {
            let _ = i;
            if m {
                let e = s - t;
                sum_sq += e * e;
            }
        }
        sum_sq.sqrt()
    }
    /// Check if the constraint is satisfied within tolerance.
    pub fn is_satisfied(&self, state: &[f64]) -> bool {
        self.violation(state) <= self.tolerance
    }
}
/// A complete trajectory: sequence of knots with metadata.
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Ordered knot points.
    pub knots: Vec<TrajectoryKnot>,
    /// State dimension.
    pub n_state: usize,
    /// Control dimension.
    pub n_control: usize,
}
impl Trajectory {
    /// Create a new empty trajectory.
    pub fn new(n_state: usize, n_control: usize) -> Self {
        Self {
            knots: Vec::new(),
            n_state,
            n_control,
        }
    }
    /// Total number of knot points.
    pub fn num_knots(&self) -> usize {
        self.knots.len()
    }
    /// Total trajectory time.
    pub fn duration(&self) -> f64 {
        if self.knots.len() < 2 {
            return 0.0;
        }
        self.knots
            .last()
            .expect("collection should not be empty")
            .time
            - self
                .knots
                .first()
                .expect("collection should not be empty")
                .time
    }
    /// Evaluate the trajectory cost using the given cost function.
    pub fn evaluate_cost(&self, cost: &CostFunction) -> f64 {
        if self.knots.len() < 2 {
            return 0.0;
        }
        let mut total = 0.0;
        for i in 0..self.knots.len() - 1 {
            let dt = self.knots[i + 1].time - self.knots[i].time;
            let rc = cost.running_cost(&self.knots[i].state, &self.knots[i].control);
            total += rc * dt;
        }
        total += cost.terminal_cost(
            &self
                .knots
                .last()
                .expect("collection should not be empty")
                .state,
        );
        total
    }
    /// Generate a linearly interpolated initial guess between two states.
    pub fn linear_initial_guess(
        n_state: usize,
        n_control: usize,
        x0: &[f64],
        xf: &[f64],
        t0: f64,
        tf: f64,
        num_knots: usize,
    ) -> Self {
        let mut knots = Vec::with_capacity(num_knots);
        for i in 0..num_knots {
            let alpha = if num_knots > 1 {
                i as f64 / (num_knots - 1) as f64
            } else {
                0.0
            };
            let t = t0 + alpha * (tf - t0);
            let state = vec_lerp(x0, xf, alpha);
            let control = vec![0.0; n_control];
            knots.push(TrajectoryKnot {
                time: t,
                state,
                control,
            });
        }
        Self {
            knots,
            n_state,
            n_control,
        }
    }
    /// Extract state at knot index.
    pub fn state_at(&self, idx: usize) -> &[f64] {
        &self.knots[idx].state
    }
    /// Extract control at knot index.
    pub fn control_at(&self, idx: usize) -> &[f64] {
        &self.knots[idx].control
    }
    /// Maximum control magnitude across all knots.
    pub fn max_control_magnitude(&self) -> f64 {
        self.knots
            .iter()
            .map(|k| norm(&k.control))
            .fold(0.0_f64, f64::max)
    }
    /// Maximum state component magnitude.
    pub fn max_state_magnitude(&self) -> f64 {
        self.knots
            .iter()
            .flat_map(|k| k.state.iter())
            .copied()
            .map(f64::abs)
            .fold(0.0_f64, f64::max)
    }
}
/// Trust region configuration.
#[derive(Debug, Clone)]
pub struct TrustRegionConfig {
    /// Initial trust region radius.
    pub initial_radius: f64,
    /// Maximum radius.
    pub max_radius: f64,
    /// Minimum radius (for convergence check).
    pub min_radius: f64,
    /// Ratio threshold to accept a step.
    pub accept_ratio: f64,
    /// Radius expansion factor on good steps.
    pub expand_factor: f64,
    /// Radius shrink factor on bad steps.
    pub shrink_factor: f64,
}
/// Boundary conditions (initial and/or final state).
#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    /// Initial state (if constrained).
    pub initial_state: Option<Vec<f64>>,
    /// Final state (if constrained).
    pub final_state: Option<Vec<f64>>,
    /// Whether the final time is free.
    pub free_final_time: bool,
    /// Bounds on final time (if free).
    pub time_bounds: Option<(f64, f64)>,
}
impl BoundaryConditions {
    /// Fixed initial and final state.
    pub fn fixed(x0: Vec<f64>, xf: Vec<f64>) -> Self {
        Self {
            initial_state: Some(x0),
            final_state: Some(xf),
            free_final_time: false,
            time_bounds: None,
        }
    }
    /// Only initial state is fixed.
    pub fn initial_only(x0: Vec<f64>) -> Self {
        Self {
            initial_state: Some(x0),
            final_state: None,
            free_final_time: false,
            time_bounds: None,
        }
    }
    /// Free final time with bounds.
    pub fn free_time(x0: Vec<f64>, xf: Vec<f64>, t_min: f64, t_max: f64) -> Self {
        Self {
            initial_state: Some(x0),
            final_state: Some(xf),
            free_final_time: true,
            time_bounds: Some((t_min, t_max)),
        }
    }
    /// Check initial state violation.
    pub fn initial_violation(&self, state: &[f64]) -> f64 {
        match &self.initial_state {
            Some(x0) => norm(&vec_sub(state, x0)),
            None => 0.0,
        }
    }
    /// Check final state violation.
    pub fn final_violation(&self, state: &[f64]) -> f64 {
        match &self.final_state {
            Some(xf) => norm(&vec_sub(state, xf)),
            None => 0.0,
        }
    }
}
/// Convergence monitor that tracks optimization progress.
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor {
    /// History of convergence entries.
    pub history: Vec<ConvergenceEntry>,
    /// Tolerance for cost convergence.
    pub cost_tol: f64,
    /// Tolerance for constraint satisfaction.
    pub constraint_tol: f64,
    /// Maximum allowed iterations.
    pub max_iterations: usize,
}
impl ConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(cost_tol: f64, constraint_tol: f64, max_iterations: usize) -> Self {
        Self {
            history: Vec::new(),
            cost_tol,
            constraint_tol,
            max_iterations,
        }
    }
    /// Record an iteration.
    pub fn record(&mut self, cost: f64, max_violation: f64, step_size: f64) {
        let cost_reduction = if let Some(prev) = self.history.last() {
            prev.cost - cost
        } else {
            0.0
        };
        let iteration = self.history.len();
        self.history.push(ConvergenceEntry {
            iteration,
            cost,
            max_violation,
            step_size,
            cost_reduction,
        });
    }
    /// Check if the optimization has converged.
    pub fn is_converged(&self) -> bool {
        if let Some(last) = self.history.last() {
            last.max_violation < self.constraint_tol
                && last.cost_reduction.abs() < self.cost_tol
                && self.history.len() > 1
        } else {
            false
        }
    }
    /// Check if the maximum iterations have been reached.
    pub fn max_reached(&self) -> bool {
        self.history.len() >= self.max_iterations
    }
    /// Latest cost value.
    pub fn latest_cost(&self) -> f64 {
        self.history.last().map(|e| e.cost).unwrap_or(f64::MAX)
    }
    /// Latest violation.
    pub fn latest_violation(&self) -> f64 {
        self.history
            .last()
            .map(|e| e.max_violation)
            .unwrap_or(f64::MAX)
    }
    /// Number of iterations completed.
    pub fn num_iterations(&self) -> usize {
        self.history.len()
    }
}
/// Result of a direct shooting simulation.
#[derive(Debug, Clone)]
pub struct ShootingResult {
    /// Simulated trajectory.
    pub trajectory: Trajectory,
    /// Final state.
    pub final_state: Vec<f64>,
    /// Total cost.
    pub cost: f64,
    /// Final state error (if boundary conditions are specified).
    pub final_error: f64,
}
/// A single phase in a multi-phase trajectory.
#[derive(Debug, Clone)]
pub struct TrajectoryPhase {
    /// Phase name (e.g., "boost", "coast", "descent").
    pub name: String,
    /// Dynamics model for this phase.
    pub dynamics: DynamicsModel,
    /// Trajectory within this phase.
    pub trajectory: Trajectory,
    /// Phase-specific dynamic limits.
    pub limits: Option<DynamicLimits>,
}
/// A continuous-time dynamics model: dx/dt = f(x, u, t).
///
/// Uses function pointers so the caller can supply arbitrary dynamics.
#[derive(Clone)]
pub struct DynamicsModel {
    /// State dimension.
    pub n_state: usize,
    /// Control dimension.
    pub n_control: usize,
    /// Dynamics function: f(x, u, t) -> dx/dt.
    pub f: fn(&[f64], &[f64], f64) -> Vec<f64>,
}
impl DynamicsModel {
    /// Create a new dynamics model.
    pub fn new(n_state: usize, n_control: usize, f: fn(&[f64], &[f64], f64) -> Vec<f64>) -> Self {
        Self {
            n_state,
            n_control,
            f,
        }
    }
    /// Evaluate the dynamics.
    pub fn eval(&self, x: &[f64], u: &[f64], t: f64) -> Vec<f64> {
        (self.f)(x, u, t)
    }
    /// Fourth-order Runge-Kutta integration step.
    pub fn rk4_step(&self, x: &[f64], u: &[f64], t: f64, dt: f64) -> Vec<f64> {
        let k1 = self.eval(x, u, t);
        let x2 = vec_add(x, &vec_scale(&k1, 0.5 * dt));
        let k2 = self.eval(&x2, u, t + 0.5 * dt);
        let x3 = vec_add(x, &vec_scale(&k2, 0.5 * dt));
        let k3 = self.eval(&x3, u, t + 0.5 * dt);
        let x4 = vec_add(x, &vec_scale(&k3, dt));
        let k4 = self.eval(&x4, u, t + dt);
        let combined: Vec<f64> = (0..x.len())
            .map(|i| k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            .collect();
        vec_add(x, &vec_scale(&combined, dt / 6.0))
    }
    /// Euler integration step.
    pub fn euler_step(&self, x: &[f64], u: &[f64], t: f64, dt: f64) -> Vec<f64> {
        let dx = self.eval(x, u, t);
        vec_add(x, &vec_scale(&dx, dt))
    }
    /// Numerical Jacobian df/dx via finite differences.
    pub fn state_jacobian(&self, x: &[f64], u: &[f64], t: f64, eps: f64) -> Vec<f64> {
        let n = self.n_state;
        let f0 = self.eval(x, u, t);
        let mut jac = vec![0.0; n * n];
        let mut x_pert = x.to_vec();
        for j in 0..n {
            x_pert[j] += eps;
            let fp = self.eval(&x_pert, u, t);
            x_pert[j] = x[j];
            for i in 0..n {
                jac[i * n + j] = (fp[i] - f0[i]) / eps;
            }
        }
        jac
    }
    /// Numerical Jacobian df/du via finite differences.
    pub fn control_jacobian(&self, x: &[f64], u: &[f64], t: f64, eps: f64) -> Vec<f64> {
        let n = self.n_state;
        let m = self.n_control;
        let f0 = self.eval(x, u, t);
        let mut jac = vec![0.0; n * m];
        let mut u_pert = u.to_vec();
        for j in 0..m {
            u_pert[j] += eps;
            let fp = self.eval(x, &u_pert, t);
            u_pert[j] = u[j];
            for i in 0..n {
                jac[i * m + j] = (fp[i] - f0[i]) / eps;
            }
        }
        jac
    }
}
/// Types of cost functions for trajectory optimization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CostType {
    /// Minimize total trajectory time.
    MinimumTime,
    /// Minimize total control energy (integral of u^T u).
    MinimumEnergy,
    /// Minimize jerk (integral of d^3 x / dt^3 squared).
    MinimumJerk,
    /// Minimize fuel consumption (integral of |u|).
    MinimumFuel,
    /// User-defined Lagrangian weighting via Q and R matrices.
    Quadratic,
}
/// Multi-phase trajectory (e.g., rocket launch with boost/coast/re-entry).
#[derive(Debug, Clone)]
pub struct MultiPhaseTrajectory {
    /// Ordered phases.
    pub phases: Vec<TrajectoryPhase>,
}
impl MultiPhaseTrajectory {
    /// Create a new multi-phase trajectory.
    pub fn new() -> Self {
        Self { phases: Vec::new() }
    }
    /// Add a phase.
    pub fn add_phase(&mut self, phase: TrajectoryPhase) {
        self.phases.push(phase);
    }
    /// Total number of knots across all phases.
    pub fn total_knots(&self) -> usize {
        self.phases.iter().map(|p| p.trajectory.num_knots()).sum()
    }
    /// Total duration across all phases.
    pub fn total_duration(&self) -> f64 {
        self.phases.iter().map(|p| p.trajectory.duration()).sum()
    }
    /// Concatenate all phases into a single trajectory.
    pub fn flatten(&self) -> Trajectory {
        if self.phases.is_empty() {
            return Trajectory::new(0, 0);
        }
        let ns = self.phases[0].trajectory.n_state;
        let nc = self.phases[0].trajectory.n_control;
        let mut traj = Trajectory::new(ns, nc);
        let mut time_offset = 0.0;
        for (phase_idx, phase) in self.phases.iter().enumerate() {
            for (i, knot) in phase.trajectory.knots.iter().enumerate() {
                if phase_idx > 0 && i == 0 {
                    continue;
                }
                traj.knots.push(TrajectoryKnot {
                    time: knot.time + time_offset,
                    state: knot.state.clone(),
                    control: knot.control.clone(),
                });
            }
            time_offset += phase.trajectory.duration();
        }
        traj
    }
    /// Evaluate phase linkage defects (continuity at phase boundaries).
    pub fn linkage_defects(&self) -> Vec<Vec<f64>> {
        let mut defects = Vec::new();
        for i in 0..self.phases.len().saturating_sub(1) {
            let end_state = self.phases[i]
                .trajectory
                .knots
                .last()
                .map(|k| k.state.clone())
                .unwrap_or_default();
            let start_state = self.phases[i + 1]
                .trajectory
                .knots
                .first()
                .map(|k| k.state.clone())
                .unwrap_or_default();
            defects.push(vec_sub(&end_state, &start_state));
        }
        defects
    }
    /// Maximum linkage defect norm.
    pub fn max_linkage_defect(&self) -> f64 {
        self.linkage_defects()
            .iter()
            .map(|d| norm(d))
            .fold(0.0_f64, f64::max)
    }
}
/// Trust region step result.
#[derive(Debug, Clone)]
pub struct TrustRegionResult {
    /// Whether the step was accepted.
    pub accepted: bool,
    /// New trust region radius.
    pub new_radius: f64,
    /// Actual reduction in cost.
    pub actual_reduction: f64,
    /// Predicted reduction in cost.
    pub predicted_reduction: f64,
    /// Ratio actual/predicted.
    pub ratio: f64,
}
/// B-spline trajectory parameterization.
///
/// The trajectory is defined by `n_cp` control points in n-dimensional space.
/// Evaluating the spline at parameter `s` in `[0, 1]` yields a smooth curve.
#[derive(Debug, Clone)]
pub struct BSplineTrajectory {
    /// Control points: each is a vector of dimension `dim`.
    pub control_points: Vec<Vec<f64>>,
    /// Dimension of each point.
    pub dim: usize,
}
impl BSplineTrajectory {
    /// Create a B-spline from control points.
    pub fn new(control_points: Vec<Vec<f64>>) -> Self {
        let dim = control_points.first().map(|v| v.len()).unwrap_or(0);
        Self {
            control_points,
            dim,
        }
    }
    /// Number of control points.
    pub fn num_control_points(&self) -> usize {
        self.control_points.len()
    }
    /// Number of spline segments.
    pub fn num_segments(&self) -> usize {
        if self.control_points.len() < 4 {
            0
        } else {
            self.control_points.len() - 3
        }
    }
    /// Evaluate the B-spline at parameter `s` in `[0, 1]`.
    pub fn evaluate(&self, s: f64) -> Vec<f64> {
        let n_seg = self.num_segments();
        if n_seg == 0 {
            return vec![0.0; self.dim];
        }
        let s_clamped = s.clamp(0.0, 1.0 - 1e-12);
        let seg_f = s_clamped * n_seg as f64;
        let seg_idx = (seg_f as usize).min(n_seg - 1);
        let t = seg_f - seg_idx as f64;
        let basis = bspline_basis(t);
        let mut result = vec![0.0; self.dim];
        for (b_idx, &b_val) in basis.iter().enumerate() {
            let cp = &self.control_points[seg_idx + b_idx];
            for d in 0..self.dim {
                result[d] += b_val * cp[d];
            }
        }
        result
    }
    /// Evaluate the derivative of the B-spline at parameter `s`.
    pub fn evaluate_derivative(&self, s: f64) -> Vec<f64> {
        let n_seg = self.num_segments();
        if n_seg == 0 {
            return vec![0.0; self.dim];
        }
        let s_clamped = s.clamp(0.0, 1.0 - 1e-12);
        let seg_f = s_clamped * n_seg as f64;
        let seg_idx = (seg_f as usize).min(n_seg - 1);
        let t = seg_f - seg_idx as f64;
        let dbasis = bspline_basis_derivative(t);
        let scale = n_seg as f64;
        let mut result = vec![0.0; self.dim];
        for (b_idx, &db_val) in dbasis.iter().enumerate() {
            let cp = &self.control_points[seg_idx + b_idx];
            for d in 0..self.dim {
                result[d] += db_val * cp[d] * scale;
            }
        }
        result
    }
    /// Sample the B-spline at `num_samples` uniformly spaced points.
    pub fn sample(&self, num_samples: usize) -> Vec<Vec<f64>> {
        (0..num_samples)
            .map(|i| {
                let s = if num_samples > 1 {
                    i as f64 / (num_samples - 1) as f64
                } else {
                    0.0
                };
                self.evaluate(s)
            })
            .collect()
    }
    /// Compute the arc length of the spline by numerical quadrature.
    pub fn arc_length(&self, num_samples: usize) -> f64 {
        if num_samples < 2 {
            return 0.0;
        }
        let mut length = 0.0;
        let ds = 1.0 / (num_samples - 1) as f64;
        let mut prev = self.evaluate(0.0);
        for i in 1..num_samples {
            let s = i as f64 * ds;
            let cur = self.evaluate(s);
            let seg_len: f64 = prev
                .iter()
                .zip(cur.iter())
                .map(|(a, b)| (b - a).powi(2))
                .sum::<f64>()
                .sqrt();
            length += seg_len;
            prev = cur;
        }
        length
    }
    /// Fit a B-spline to a sequence of waypoints (least-squares fit).
    ///
    /// Produces `num_cp` control points that approximate the waypoints.
    pub fn fit_to_waypoints(waypoints: &[Vec<f64>], num_cp: usize) -> Self {
        if waypoints.is_empty() || num_cp < 4 {
            return Self {
                control_points: vec![vec![0.0; 3]; num_cp.max(4)],
                dim: 3,
            };
        }
        let dim = waypoints[0].len();
        let mut cps = Vec::with_capacity(num_cp);
        for i in 0..num_cp {
            let alpha = i as f64 / (num_cp - 1) as f64;
            let idx_f = alpha * (waypoints.len() - 1) as f64;
            let idx = (idx_f as usize).min(waypoints.len() - 2);
            let t = idx_f - idx as f64;
            let pt = vec_lerp(&waypoints[idx], &waypoints[idx + 1], t);
            cps.push(pt);
        }
        Self {
            control_points: cps,
            dim,
        }
    }
}
/// A single trajectory knot point (state + control at a given time).
#[derive(Debug, Clone)]
pub struct TrajectoryKnot {
    /// Time at this knot.
    pub time: f64,
    /// State vector.
    pub state: Vec<f64>,
    /// Control vector (may be empty at the final knot).
    pub control: Vec<f64>,
}
