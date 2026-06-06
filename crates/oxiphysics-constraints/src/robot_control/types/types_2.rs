//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::functions::*;
use std::f64::consts::PI;

/// A trajectory sample: (time, positions, velocities, accelerations).
type TrajectorySample = (f64, Vec<f64>, Vec<f64>, Vec<f64>);

use super::types_impl::{Jacobian, VelocityLimits};

/// Inverse kinematics solver using Jacobian pseudoinverse and null-space optimisation.
#[derive(Debug, Clone)]
pub struct InverseKinematics {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance on position error.
    pub tolerance: f64,
    /// Step size for joint angle updates.
    pub step_size: f64,
    /// Damping coefficient for damped least squares.
    pub lambda: f64,
    /// Weight on null-space joint-limit avoidance.
    pub null_space_weight: f64,
}
impl InverseKinematics {
    /// Create a new IK solver with default parameters.
    pub fn new() -> Self {
        Self {
            max_iter: 200,
            tolerance: 1e-4,
            step_size: 0.5,
            lambda: 0.01,
            null_space_weight: 0.1,
        }
    }
    /// Solve IK for the target end-effector position `target` (3-D).
    ///
    /// Modifies `robot.joint_angles` in place.
    /// Returns `true` if converged within tolerance.
    pub fn solve_position(&self, robot: &mut SerialManipulator, target: [f64; 3]) -> bool {
        for _ in 0..self.max_iter {
            let pos = robot.end_effector_position();
            let err = [target[0] - pos[0], target[1] - pos[1], target[2] - pos[2]];
            let err_norm = norm3(err);
            if err_norm < self.tolerance {
                return true;
            }
            let jac = Jacobian::compute(robot);
            let n = robot.dof();
            let mut j3 = vec![0.0f64; 3 * n];
            for r in 0..3 {
                for c in 0..n {
                    j3[r * n + c] = jac.j[r * n + c];
                }
            }
            let jp = damped_pseudoinverse_3xn(&j3, n, self.lambda);
            let mut dq = vec![0.0f64; n];
            for i in 0..n {
                let mut s = 0.0;
                for r in 0..3 {
                    s += jp[i * 3 + r] * err[r];
                }
                dq[i] = self.step_size * s;
            }
            let ns = null_space_joint_limit_gradient(robot);
            let projected = project_null_space(&jp, &j3, n, &ns, self.null_space_weight);
            for i in 0..n {
                let (lo, hi) = robot.joint_limits[i];
                robot.joint_angles[i] = (robot.joint_angles[i] + dq[i] + projected[i])
                    .max(lo)
                    .min(hi);
            }
        }
        false
    }
}
/// Resolved-rate (Jacobian-based velocity control) at the task level.
///
/// Computes joint velocity: `q̇ = J† * ẋ_d`
/// where `J†` is the Moore-Penrose pseudoinverse of the Jacobian.
#[derive(Debug, Clone)]
pub struct ResolvedRateController {
    /// Damping factor for damped least-squares (DLS) pseudoinverse.
    pub damping: f64,
    /// Joint velocity limits.
    pub vel_limits: VelocityLimits,
}
impl ResolvedRateController {
    /// Create with DLS damping factor and per-joint velocity limits.
    pub fn new(damping: f64, vel_limits: VelocityLimits) -> Self {
        Self {
            damping,
            vel_limits,
        }
    }
    /// Compute joint velocities from desired Cartesian velocity.
    ///
    /// `jacobian` — m×n Jacobian (row-major, m task DOFs, n joint DOFs).
    /// `task_vel` — desired task-space velocity (length m).
    ///
    /// Returns joint velocity vector of length n.
    pub fn compute(&self, jacobian: &[Vec<f64>], task_vel: &[f64]) -> Vec<f64> {
        let m = jacobian.len();
        if m == 0 {
            return Vec::new();
        }
        let n = jacobian[0].len();
        let mut a = vec![vec![0.0_f64; m]; m];
        for i in 0..m {
            for j in 0..m {
                let sum: f64 = jacobian[i]
                    .iter()
                    .take(n)
                    .zip(jacobian[j].iter().take(n))
                    .map(|(a, b)| a * b)
                    .sum();
                a[i][j] = sum;
            }
            a[i][i] += self.damping * self.damping;
        }
        let y = gauss_seidel_solve(&a, task_vel, 50);
        let mut q_dot = vec![0.0; n];
        for j in 0..n {
            for i in 0..m {
                q_dot[j] += jacobian[i][j] * y[i];
            }
        }
        self.vel_limits.clamp(&q_dot)
    }
}
/// Null-space projector for redundancy resolution.
///
/// Projects secondary task velocities into the null-space of the primary task.
#[derive(Debug, Clone)]
pub struct NullSpaceProjector {
    /// Primary Jacobian (m×n row-major).
    pub jacobian: Vec<Vec<f64>>,
    /// DLS damping factor for pseudoinverse.
    pub damping: f64,
    /// Weight matrix (n×n diagonal, for weighted pseudoinverse).
    pub weights: Vec<f64>,
}
impl NullSpaceProjector {
    /// Create with Jacobian, damping, and unit weights.
    pub fn new(jacobian: Vec<Vec<f64>>, damping: f64) -> Self {
        let n = if jacobian.is_empty() {
            0
        } else {
            jacobian[0].len()
        };
        Self {
            jacobian,
            damping,
            weights: vec![1.0; n],
        }
    }
    /// Compute null-space projection matrix N = I - J† J.
    ///
    /// Returns n×n matrix (row-major flat Vec).
    pub fn null_space_matrix(&self) -> Vec<f64> {
        let n = self.weights.len();
        if n == 0 || self.jacobian.is_empty() {
            return Vec::new();
        }
        let m = self.jacobian.len();
        let mut a = vec![0.0; m * m];
        for i in 0..m {
            for j in 0..m {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.jacobian[i][k] * self.jacobian[j][k];
                }
                a[i * m + j] = s;
            }
            a[i * m + i] += self.damping * self.damping;
        }
        let a_inv = invert_dense_matrix(&a, m);
        let mut n_mat = vec![0.0; n * n];
        for i in 0..n {
            n_mat[i * n + i] = 1.0;
        }
        for i in 0..n {
            for k in 0..m {
                let mut jt_ik = 0.0;
                for kk in 0..m {
                    jt_ik += self.jacobian[kk][i] * a_inv[kk * m + k];
                }
                for l in 0..n {
                    n_mat[i * n + l] -= jt_ik * self.jacobian[k][l];
                }
            }
        }
        n_mat
    }
    /// Project a secondary velocity into the null-space.
    ///
    /// `q_dot_secondary` — secondary task velocity (length n).
    ///
    /// Returns `N * q_dot_secondary`.
    pub fn project(&self, q_dot_secondary: &[f64]) -> Vec<f64> {
        let n = q_dot_secondary.len();
        if n == 0 {
            return Vec::new();
        }
        let n_mat = self.null_space_matrix();
        if n_mat.is_empty() {
            return q_dot_secondary.to_vec();
        }
        let actual_n = (n_mat.len() as f64).sqrt() as usize;
        let mut result = vec![0.0; actual_n];
        for i in 0..actual_n {
            for j in 0..n.min(actual_n) {
                result[i] += n_mat[i * actual_n + j] * q_dot_secondary[j];
            }
        }
        result
    }
    /// Compute gradient for joint limit avoidance (used as secondary task).
    ///
    /// `q` — current joint positions.
    /// `q_min`, `q_max` — joint limits.
    ///
    /// Returns gradient ∂H/∂q where H = Σ (q_i - q_mid_i)² / (q_range_i)²
    pub fn joint_limit_gradient(q: &[f64], q_min: &[f64], q_max: &[f64]) -> Vec<f64> {
        q.iter()
            .zip(q_min.iter().zip(q_max.iter()))
            .map(|(qi, (lo, hi))| {
                let range = hi - lo;
                let mid = 0.5 * (lo + hi);
                if range > 1e-15 {
                    (qi - mid) / (range * range)
                } else {
                    0.0
                }
            })
            .collect()
    }
}
/// Task-space (Cartesian) trajectory types.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CartesianTrajectoryType {
    /// Straight line between two points.
    Linear,
    /// Circular arc: centre + radius + angle range.
    Circular,
    /// Screw motion along an axis.
    Screw,
}
/// Detects transitions between free motion and contact using F/T data.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum ContactPhase {
    /// No contact detected.
    #[default]
    FreeMotion,
    /// Contact detected but below threshold for stable grasp.
    LightContact,
    /// Stable contact with significant normal force.
    StableContact,
    /// Contact force exceeding safety threshold.
    OverforceDetected,
}
/// Result of workspace analysis for a manipulator.
#[derive(Debug, Clone)]
pub struct WorkspaceAnalysis {
    /// Approximate reachable workspace radius.
    pub max_reach: f64,
    /// Approximate minimum reach (dexterous workspace inner boundary).
    pub min_reach: f64,
    /// Number of sampled configurations in analysis.
    pub n_samples: usize,
    /// Fraction of configurations that are non-singular.
    pub non_singular_fraction: f64,
    /// Approximate workspace volume (unit³).
    pub volume_estimate: f64,
}
impl WorkspaceAnalysis {
    /// Analyze workspace from forward kinematics samples.
    ///
    /// `positions` — list of (x,y,z) end-effector positions from FK sampling.
    /// `manipulabilities` — manipulability measure for each sample.
    pub fn from_samples(positions: &[[f64; 3]], manipulabilities: &[f64]) -> Self {
        if positions.is_empty() {
            return Self {
                max_reach: 0.0,
                min_reach: 0.0,
                n_samples: 0,
                non_singular_fraction: 0.0,
                volume_estimate: 0.0,
            };
        }
        let n = positions.len();
        let reaches: Vec<f64> = positions
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
            .collect();
        let max_reach = reaches.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_reach = reaches.iter().cloned().fold(f64::INFINITY, f64::min);
        let non_singular = if manipulabilities.len() == n {
            manipulabilities.iter().filter(|&&m| m > 1e-6).count() as f64 / n as f64
        } else {
            1.0
        };
        let sphere_vol = (4.0 / 3.0) * std::f64::consts::PI * max_reach * max_reach * max_reach;
        let volume_estimate = sphere_vol * non_singular;
        Self {
            max_reach,
            min_reach,
            n_samples: n,
            non_singular_fraction: non_singular,
            volume_estimate,
        }
    }
    /// Check if a point is approximately within the workspace.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        let r = (point[0] * point[0] + point[1] * point[1] + point[2] * point[2]).sqrt();
        r >= self.min_reach && r <= self.max_reach
    }
}
/// An n-DOF serial robot manipulator defined by a chain of DH links.
///
/// Joint limits are given as `(min, max)` pairs in radians.
#[derive(Debug, Clone)]
pub struct SerialManipulator {
    /// Denavit-Hartenberg parameters for each link.
    pub links: Vec<DhParams>,
    /// Joint limits `(lo, hi)` in radians for each joint.
    pub joint_limits: Vec<(f64, f64)>,
    /// Current joint angles (radians).
    pub joint_angles: Vec<f64>,
}
impl SerialManipulator {
    /// Construct a serial manipulator from a vector of DH parameters.
    ///
    /// Joint angles are initialised to zero; limits default to `[-π, π]`.
    pub fn new(links: Vec<DhParams>) -> Self {
        let n = links.len();
        Self {
            links,
            joint_limits: vec![(-PI, PI); n],
            joint_angles: vec![0.0; n],
        }
    }
    /// Number of degrees of freedom.
    pub fn dof(&self) -> usize {
        self.links.len()
    }
    /// Set all joint angles, clamping to limits.
    pub fn set_joint_angles(&mut self, angles: &[f64]) {
        for (i, &a) in angles.iter().enumerate().take(self.dof()) {
            let (lo, hi) = self.joint_limits[i];
            self.joint_angles[i] = a.max(lo).min(hi);
        }
    }
    /// Forward kinematics: return the end-effector 4×4 transform (row-major).
    ///
    /// Each DH link's `theta` is offset by the corresponding joint angle.
    pub fn end_effector_pose(&self) -> [f64; 16] {
        let mut t = mat4_identity();
        for (i, link) in self.links.iter().enumerate() {
            let theta = link.theta + self.joint_angles[i];
            let ti = dh_matrix(link.a, link.d, link.alpha, theta);
            t = mat4_mul(t, ti);
        }
        t
    }
    /// Return the 3-D end-effector position `[x, y, z]`.
    pub fn end_effector_position(&self) -> [f64; 3] {
        let t = self.end_effector_pose();
        mat4_translation(t)
    }
}
/// Online payload mass and CoM identification for gravity compensation.
///
/// Uses recursive least-squares (RLS) to identify `[m, m*cx, m*cy, m*cz]`
/// from end-effector wrench measurements at multiple robot configurations.
#[derive(Debug, Clone)]
pub struct PayloadIdentifier {
    /// Recursive least-squares covariance matrix (4×4, flat row-major).
    pub cov: [f64; 16],
    /// Parameter vector `[m, m*cx, m*cy, m*cz]`.
    pub params: [f64; 4],
    /// RLS forgetting factor.
    pub forgetting: f64,
    /// Number of measurements processed.
    pub n_measurements: usize,
}
impl PayloadIdentifier {
    /// Create with initial covariance `cov0 * I` and zero parameters.
    pub fn new(cov0: f64, forgetting: f64) -> Self {
        let mut cov = [0.0; 16];
        for i in 0..4 {
            cov[i * 4 + i] = cov0;
        }
        Self {
            cov,
            params: [0.0; 4],
            forgetting,
            n_measurements: 0,
        }
    }
    /// Update with one measurement.
    ///
    /// `regressor` — 3×4 regression matrix for current configuration.
    /// `force_meas` — measured force \[Fx, Fy, Fz\] (gravity direction * m).
    pub fn update(&mut self, regressor: &[[f64; 4]; 3], force_meas: &[f64; 3]) {
        for axis in 0..3 {
            let phi = regressor[axis];
            let y_pred: f64 = phi.iter().zip(self.params.iter()).map(|(h, p)| h * p).sum();
            let innov = force_meas[axis] - y_pred;
            let mut p_phi = [0.0; 4];
            for (i, p_phi_i) in p_phi.iter_mut().enumerate() {
                for (j, phi_j) in phi.iter().enumerate() {
                    *p_phi_i += self.cov[i * 4 + j] * phi_j;
                }
            }
            let denom = self.forgetting
                + phi
                    .iter()
                    .zip(p_phi.iter())
                    .map(|(h, pp)| h * pp)
                    .sum::<f64>();
            let gain: Vec<f64> = p_phi.iter().map(|pp| pp / denom).collect();
            for (param_i, g_i) in self.params.iter_mut().zip(gain.iter()) {
                *param_i += g_i * innov;
            }
            let mut new_cov = [0.0; 16];
            for i in 0..4 {
                for j in 0..4 {
                    new_cov[i * 4 + j] =
                        (self.cov[i * 4 + j] - gain[i] * p_phi[j]) / self.forgetting;
                }
            }
            self.cov = new_cov;
        }
        self.n_measurements += 1;
    }
    /// Estimated payload mass (kg).
    pub fn mass(&self) -> f64 {
        self.params[0]
    }
    /// Estimated payload center of mass (m).
    pub fn com(&self) -> [f64; 3] {
        if self.params[0].abs() > 1e-6 {
            [
                self.params[1] / self.params[0],
                self.params[2] / self.params[0],
                self.params[3] / self.params[0],
            ]
        } else {
            [0.0; 3]
        }
    }
}
/// Low-pass filter for F/T sensor signals.
#[derive(Debug, Clone)]
pub struct FtSensorFilter {
    /// Filter cutoff frequency (Hz).
    pub cutoff_hz: f64,
    /// Sample rate (Hz).
    pub sample_hz: f64,
    /// Alpha coefficient for first-order IIR filter.
    pub(super) alpha: f64,
    /// Current filtered state.
    pub(super) state: [f64; 6],
    /// Offset bias (for taring).
    pub(super) bias: [f64; 6],
}
impl FtSensorFilter {
    /// Create with given cutoff and sample frequencies.
    pub fn new(cutoff_hz: f64, sample_hz: f64) -> Self {
        let rc = 1.0 / (2.0 * std::f64::consts::PI * cutoff_hz);
        let dt = 1.0 / sample_hz;
        let alpha = dt / (rc + dt);
        Self {
            cutoff_hz,
            sample_hz,
            alpha,
            state: [0.0; 6],
            bias: [0.0; 6],
        }
    }
    /// Update filter with new raw sample, returns filtered value.
    pub fn update(&mut self, raw: [f64; 6]) -> [f64; 6] {
        for (state_i, (raw_i, bias_i)) in
            self.state.iter_mut().zip(raw.iter().zip(self.bias.iter()))
        {
            let debiased = raw_i - bias_i;
            *state_i = self.alpha * debiased + (1.0 - self.alpha) * *state_i;
        }
        self.state
    }
    /// Tare sensor (set current reading as zero bias).
    pub fn tare(&mut self, raw: [f64; 6]) {
        self.bias = raw;
        self.state = [0.0; 6];
    }
    /// Get current filtered state.
    pub fn filtered(&self) -> [f64; 6] {
        self.state
    }
}
/// Multi-joint quintic trajectory.
#[derive(Debug, Clone)]
pub struct JointTrajectory {
    /// Per-joint quintic polynomials.
    pub joints: Vec<QuinticTrajectory>,
}
impl JointTrajectory {
    /// Create from start/end configurations.
    pub fn new(q0: &[f64], qf: &[f64], duration: f64) -> Self {
        let joints = q0
            .iter()
            .zip(qf.iter())
            .map(|(s, e)| QuinticTrajectory::new(*s, *e, duration))
            .collect();
        Self { joints }
    }
    /// Evaluate joint positions at time `t`.
    pub fn position(&self, t: f64) -> Vec<f64> {
        self.joints.iter().map(|j| j.position(t)).collect()
    }
    /// Evaluate joint velocities at time `t`.
    pub fn velocity(&self, t: f64) -> Vec<f64> {
        self.joints.iter().map(|j| j.velocity(t)).collect()
    }
    /// Evaluate joint accelerations at time `t`.
    pub fn acceleration(&self, t: f64) -> Vec<f64> {
        self.joints.iter().map(|j| j.acceleration(t)).collect()
    }
    /// Duration of the trajectory.
    pub fn duration(&self) -> f64 {
        self.joints.first().map(|j| j.duration).unwrap_or(0.0)
    }
    /// Sample the trajectory at uniform time intervals.
    pub fn sample(&self, n_steps: usize) -> Vec<TrajectorySample> {
        let dur = self.duration();
        if n_steps == 0 {
            return Vec::new();
        }
        (0..=n_steps)
            .map(|i| {
                let t = dur * i as f64 / n_steps as f64;
                (t, self.position(t), self.velocity(t), self.acceleration(t))
            })
            .collect()
    }
}
/// Safety categories for robot motion.
#[derive(Debug, Clone, PartialEq)]
pub enum SafetyLevel {
    /// Normal operation — all safety checks passed.
    Normal,
    /// Warning state — some limits approached.
    Warning,
    /// Emergency stop required.
    EmergencyStop,
}
/// Quintic polynomial trajectory for smooth point-to-point motion.
///
/// Satisfies boundary conditions: q(0)=q0, q(T)=qf, q'(0)=q'(T)=0, q''(0)=q''(T)=0
#[derive(Debug, Clone)]
pub struct QuinticTrajectory {
    /// Start position.
    pub q0: f64,
    /// End position.
    pub qf: f64,
    /// Motion duration (s).
    pub duration: f64,
    /// Polynomial coefficients \[a0..a5\].
    pub(super) coeffs: [f64; 6],
}
impl QuinticTrajectory {
    /// Create a quintic trajectory from `q0` to `qf` in `duration` seconds.
    ///
    /// Zero boundary velocities and accelerations.
    pub fn new(q0: f64, qf: f64, duration: f64) -> Self {
        let t = duration;
        let h = qf - q0;
        let t3 = t * t * t;
        let t4 = t3 * t;
        let t5 = t4 * t;
        let coeffs = [q0, 0.0, 0.0, 10.0 * h / t3, -15.0 * h / t4, 6.0 * h / t5];
        Self {
            q0,
            qf,
            duration,
            coeffs,
        }
    }
    /// Evaluate position at time `t` (clamped to \[0, duration\]).
    pub fn position(&self, t: f64) -> f64 {
        let t = t.max(0.0).min(self.duration);
        let c = &self.coeffs;
        c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * c[5]))))
    }
    /// Evaluate velocity at time `t`.
    pub fn velocity(&self, t: f64) -> f64 {
        let t = t.max(0.0).min(self.duration);
        let c = &self.coeffs;
        c[1] + t * (2.0 * c[2] + t * (3.0 * c[3] + t * (4.0 * c[4] + t * 5.0 * c[5])))
    }
    /// Evaluate acceleration at time `t`.
    pub fn acceleration(&self, t: f64) -> f64 {
        let t = t.max(0.0).min(self.duration);
        let c = &self.coeffs;
        2.0 * c[2] + t * (6.0 * c[3] + t * (12.0 * c[4] + t * 20.0 * c[5]))
    }
    /// Maximum velocity magnitude (approximate, evaluated at midpoint).
    pub fn peak_velocity(&self) -> f64 {
        self.velocity(self.duration / 2.0).abs()
    }
}
/// Broyden-Fletcher adaptive Jacobian update for model-free control.
///
/// Implements the rank-1 Broyden update:
///
/// `J_{k+1} = J_k + (Δy - J_k Δu) Δu^T / ||Δu||²`
#[derive(Debug, Clone)]
pub struct BroydenJacobianUpdate {
    /// Current Jacobian estimate (m×n, row-major).
    pub jacobian: Vec<Vec<f64>>,
    /// Task space dimension m.
    pub m: usize,
    /// Joint space dimension n.
    pub n: usize,
    /// Forgetting factor (0.9–1.0).
    pub forgetting: f64,
}
impl BroydenJacobianUpdate {
    /// Create with initial Jacobian estimate.
    pub fn new(jacobian: Vec<Vec<f64>>, forgetting: f64) -> Self {
        let m = jacobian.len();
        let n = if m > 0 { jacobian[0].len() } else { 0 };
        Self {
            jacobian,
            m,
            n,
            forgetting,
        }
    }
    /// Initialize with zero Jacobian.
    pub fn zeros(m: usize, n: usize) -> Self {
        Self {
            jacobian: vec![vec![0.0; n]; m],
            m,
            n,
            forgetting: 0.99,
        }
    }
    /// Perform a Broyden rank-1 update.
    ///
    /// `delta_q` — change in joint positions (Δu).
    /// `delta_x` — observed change in task positions (Δy).
    pub fn update(&mut self, delta_q: &[f64], delta_x: &[f64]) {
        let norm2 = dot_n(delta_q, delta_q);
        if norm2 < 1e-20 {
            return;
        }
        let j_delta_q: Vec<f64> = (0..self.m)
            .map(|i| {
                self.jacobian[i]
                    .iter()
                    .zip(delta_q.iter())
                    .map(|(j, dq)| j * dq)
                    .sum()
            })
            .collect();
        let residual: Vec<f64> = delta_x
            .iter()
            .zip(j_delta_q.iter())
            .map(|(dy, jdq)| dy - jdq)
            .collect();
        for (i, (jac_row, res_i)) in self
            .jacobian
            .iter_mut()
            .zip(residual.iter())
            .enumerate()
            .take(self.m)
        {
            let _ = i;
            for (jac_ij, dq_j) in jac_row.iter_mut().zip(delta_q.iter()).take(self.n) {
                *jac_ij = self.forgetting * *jac_ij + res_i * dq_j / norm2;
            }
        }
    }
    /// Get current Jacobian as a flat row-major Vec.
    pub fn flat(&self) -> Vec<f64> {
        self.jacobian
            .iter()
            .flat_map(|row| row.iter().cloned())
            .collect()
    }
}
/// Denavit-Hartenberg parameters for a single robot link.
///
/// The standard DH convention defines four parameters per joint:
/// `a` (link length), `d` (link offset), `alpha` (link twist), and `theta` (joint angle).
#[derive(Debug, Clone, Copy)]
pub struct DhParams {
    /// Link length (metres).
    pub a: f64,
    /// Link offset (metres).
    pub d: f64,
    /// Link twist (radians).
    pub alpha: f64,
    /// Joint angle (radians); for revolute joints this is the variable.
    pub theta: f64,
}
impl DhParams {
    /// Construct a new set of DH parameters.
    pub fn new(a: f64, d: f64, alpha: f64, theta: f64) -> Self {
        Self { a, d, alpha, theta }
    }
    /// Compute the 4×4 homogeneous transformation matrix for this DH link.
    ///
    /// Returns a row-major `[f64; 16]` matrix using the standard DH formula:
    /// `T = Rz(θ) · Tz(d) · Tx(a) · Rx(α)`.
    pub fn forward_kinematics(&self) -> [f64; 16] {
        dh_matrix(self.a, self.d, self.alpha, self.theta)
    }
}
/// Composite adaptive controller with online parameter estimation.
///
/// Based on the Slotine-Li regressor framework for robot manipulators.
#[derive(Debug, Clone)]
pub struct AdaptiveRobotControl {
    /// Estimated inertia parameters (one per DOF, simplified scalar).
    pub theta_hat: Vec<f64>,
    /// Adaptation gain matrix diagonal.
    pub gamma: Vec<f64>,
    /// PD control gain matrix diagonal (position).
    pub kp: Vec<f64>,
    /// PD control gain matrix diagonal (velocity).
    pub kd: Vec<f64>,
    /// Sliding variable history.
    pub sigma: Vec<f64>,
}
impl AdaptiveRobotControl {
    /// Create adaptive controller with initial parameter estimates.
    pub fn new(n: usize, gamma: f64, kp: f64, kd: f64) -> Self {
        Self {
            theta_hat: vec![1.0; n],
            gamma: vec![gamma; n],
            kp: vec![kp; n],
            kd: vec![kd; n],
            sigma: vec![0.0; n],
        }
    }
    /// Update parameter estimates and compute control torques.
    ///
    /// * `q` - current joint angles
    /// * `qdot` - current joint velocities
    /// * `q_d` - desired joint angles
    /// * `qdot_d` - desired joint velocities
    /// * `qdotdot_d` - desired joint accelerations
    /// * `dt` - time step
    ///
    /// Returns computed joint torques.
    pub fn update(
        &mut self,
        q: &[f64],
        qdot: &[f64],
        q_d: &[f64],
        qdot_d: &[f64],
        qdotdot_d: &[f64],
        dt: f64,
    ) -> Vec<f64> {
        let n = self.theta_hat.len();
        let lambda = 10.0;
        let mut torques = vec![0.0f64; n];
        for i in 0..n {
            let qi = if i < q.len() { q[i] } else { 0.0 };
            let qdoti = if i < qdot.len() { qdot[i] } else { 0.0 };
            let qdi = if i < q_d.len() { q_d[i] } else { 0.0 };
            let qdotdi = if i < qdot_d.len() { qdot_d[i] } else { 0.0 };
            let qdotdotdi = if i < qdotdot_d.len() {
                qdotdot_d[i]
            } else {
                0.0
            };
            let q_err = qdi - qi;
            let qdot_err = qdotdi - qdoti;
            let s = qdot_err + lambda * q_err;
            self.sigma[i] = s;
            let _qdot_r = qdotdi + lambda * q_err;
            let qdotdot_r = qdotdotdi + lambda * qdot_err;
            let y_i = qdotdot_r;
            self.theta_hat[i] += self.gamma[i] * y_i * s * dt;
            torques[i] = self.theta_hat[i] * y_i + self.kd[i] * s;
        }
        torques
    }
}
/// Null-space optimization criteria for redundant manipulators.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NullSpaceCriteria {
    /// Minimize joint velocity norm (minimum motion).
    MinVelocity,
    /// Maximize distance from joint limits (gradient-based).
    JointLimitAvoidance,
    /// Maximize manipulability measure.
    ManipulabilityMaximization,
    /// Avoid singularities by maximizing min singular value.
    SingularityAvoidance,
    /// Minimize kinetic energy.
    MinKineticEnergy,
}
/// Virtual spring-damper environment for stiffness rendering (haptics).
#[derive(Debug, Clone)]
pub struct StiffnessRenderer {
    /// Virtual wall position (m).
    pub wall_position: f64,
    /// Normal direction (1D: +1 or -1).
    pub wall_normal: f64,
    /// Virtual stiffness (N/m).
    pub stiffness: f64,
    /// Virtual damping (N·s/m).
    pub damping: f64,
    /// Maximum rendered force (N).
    pub max_force: f64,
}
impl StiffnessRenderer {
    /// Create a haptic wall.
    pub fn new(wall_position: f64, wall_normal: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            wall_position,
            wall_normal,
            stiffness,
            damping,
            max_force: 50.0,
        }
    }
    /// Compute rendered force at given position and velocity.
    ///
    /// Returns zero if not in contact with virtual wall.
    pub fn force(&self, pos: f64, vel: f64) -> f64 {
        let penetration = self.wall_normal * (pos - self.wall_position);
        if penetration <= 0.0 {
            return 0.0;
        }
        let f = self.wall_normal * (-self.stiffness * penetration - self.damping * vel);
        f.max(-self.max_force).min(self.max_force)
    }
    /// Check if the haptic device is inside the virtual wall.
    pub fn in_contact(&self, pos: f64) -> bool {
        self.wall_normal * (pos - self.wall_position) > 0.0
    }
}
/// Feedforward torque computation for gravity and inertia compensation.
#[derive(Debug, Clone)]
pub struct TorqueFeedforward {
    /// Mass of each link (kg).
    pub link_masses: Vec<f64>,
    /// Centre-of-mass offset for each link (along z_i).
    pub com_offsets: Vec<f64>,
    /// Gravity vector in world frame.
    pub gravity: [f64; 3],
}
impl TorqueFeedforward {
    /// Construct with given masses and gravity vector.
    pub fn new(link_masses: Vec<f64>, com_offsets: Vec<f64>, gravity: [f64; 3]) -> Self {
        Self {
            link_masses,
            com_offsets,
            gravity,
        }
    }
    /// Compute gravity compensation torques for the given robot configuration.
    ///
    /// Uses the recursive Newton-Euler method (simplified for revolute-only chains).
    pub fn gravity_torques(&self, robot: &SerialManipulator) -> Vec<f64> {
        let n = robot.dof();
        let mut torques = vec![0.0f64; n];
        let mut transforms = vec![mat4_identity(); n + 1];
        for i in 0..n {
            let theta = robot.links[i].theta + robot.joint_angles[i];
            let ti = dh_matrix(
                robot.links[i].a,
                robot.links[i].d,
                robot.links[i].alpha,
                theta,
            );
            transforms[i + 1] = mat4_mul(transforms[i], ti);
        }
        for j in 0..n {
            let z_j = [transforms[j][2], transforms[j][6], transforms[j][10]];
            let p_j = mat4_translation(transforms[j]);
            let mut tau = 0.0;
            for k in j..n {
                let m_k = if k < self.link_masses.len() {
                    self.link_masses[k]
                } else {
                    0.0
                };
                let p_k = mat4_translation(transforms[k + 1]);
                let r_jk = [p_k[0] - p_j[0], p_k[1] - p_j[1], p_k[2] - p_j[2]];
                let z_cross_r = cross3(z_j, r_jk);
                tau += m_k * dot3(self.gravity, z_cross_r);
            }
            torques[j] = -tau;
        }
        torques
    }
    /// Compute inertia compensation torques given desired joint accelerations.
    ///
    /// Simplified: uses diagonal inertia approximation `τ = I * q̈`.
    pub fn inertia_torques(&self, robot: &SerialManipulator, q_ddot: &[f64]) -> Vec<f64> {
        let n = robot.dof();
        (0..n)
            .map(|i| {
                let m = if i < self.link_masses.len() {
                    self.link_masses[i]
                } else {
                    1.0
                };
                let a = robot.links[i].a;
                let inertia = m * a * a + 0.1;
                let qdd = if i < q_ddot.len() { q_ddot[i] } else { 0.0 };
                inertia * qdd
            })
            .collect()
    }
}
/// 6-DOF Cartesian impedance control with mass-spring-damper dynamics.
///
/// Implements: F_ext = M_d * ẍ + D_d * ẋ + K_d * x_err
/// The controller computes the required joint torques via the Jacobian transpose.
#[derive(Debug, Clone)]
pub struct CartesianImpedanceController {
    /// Desired Cartesian stiffness (6×6 diagonal, \[tx,ty,tz,rx,ry,rz\]).
    pub k_d: [f64; 6],
    /// Desired Cartesian damping (6×6 diagonal).
    pub b_d: [f64; 6],
    /// Desired Cartesian inertia (6×6 diagonal).
    pub m_d: [f64; 6],
    /// Maximum wrench components.
    pub wrench_limits: [f64; 6],
}
impl CartesianImpedanceController {
    /// Create with separate stiffness, damping, and inertia diagonals.
    pub fn new(k_d: [f64; 6], b_d: [f64; 6], m_d: [f64; 6]) -> Self {
        Self {
            k_d,
            b_d,
            m_d,
            wrench_limits: [1000.0; 6],
        }
    }
    /// Critically-damped controller for given stiffness diagonal.
    pub fn critically_damped(k_d: [f64; 6], m_d: [f64; 6]) -> Self {
        let b_d = std::array::from_fn(|i| 2.0 * (m_d[i] * k_d[i]).sqrt());
        Self::new(k_d, b_d, m_d)
    }
    /// Compute Cartesian wrench from pose/velocity error.
    ///
    /// `pos_err` — 6D pose error \[dx,dy,dz,dRx,dRy,dRz\].
    /// `vel_err` — 6D velocity error.
    /// `acc_des` — 6D desired acceleration (feedforward).
    pub fn compute_wrench(
        &self,
        pos_err: &[f64; 6],
        vel_err: &[f64; 6],
        acc_des: &[f64; 6],
    ) -> [f64; 6] {
        let mut w = [0.0; 6];
        for i in 0..6 {
            w[i] = self.m_d[i] * acc_des[i] + self.b_d[i] * vel_err[i] + self.k_d[i] * pos_err[i];
            w[i] = w[i].max(-self.wrench_limits[i]).min(self.wrench_limits[i]);
        }
        w
    }
    /// Convert Cartesian wrench to joint torques via J^T * F.
    ///
    /// `jacobian` — 6×n Jacobian (row-major).
    /// `wrench` — 6D wrench.
    pub fn wrench_to_torques(&self, jacobian: &[Vec<f64>], wrench: &[f64; 6]) -> Vec<f64> {
        if jacobian.is_empty() {
            return Vec::new();
        }
        let n = jacobian[0].len();
        let mut tau = vec![0.0; n];
        for j in 0..n {
            for i in 0..6.min(jacobian.len()) {
                tau[j] += jacobian[i][j] * wrench[i];
            }
        }
        tau
    }
    /// Natural frequency for joint `i` (rad/s).
    pub fn natural_frequency(&self, i: usize) -> f64 {
        if i < 6 && self.m_d[i] > 1e-15 {
            (self.k_d[i] / self.m_d[i]).sqrt()
        } else {
            0.0
        }
    }
    /// Damping ratio for joint `i`.
    pub fn damping_ratio(&self, i: usize) -> f64 {
        let wn = self.natural_frequency(i);
        if i < 6 && wn > 1e-15 {
            self.b_d[i] / (2.0 * self.m_d[i] * wn)
        } else {
            0.0
        }
    }
}
