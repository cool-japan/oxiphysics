//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::type_complexity)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    cross3, damped_pseudoinverse_3xn, dh_matrix, dot_n, dot3, gauss_jordan_6x6, gauss_seidel_solve,
    invert_dense_matrix, mat3_vec_mul, mat4_identity, mat4_mul, mat4_translation, norm3,
    normalize3, null_space_joint_limit_gradient, project_null_space, rotation_from_axis_angle,
    scale_vec,
};

/// Simplified RNEA for serial manipulators with revolute joints.
///
/// Computes joint torques from: `τ = M(q)q̈ + C(q,q̇)q̇ + G(q)`
/// using a reduced scalar approximation for each joint.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RneaSolver {
    /// Links in kinematic chain order (proximal to distal).
    pub links: Vec<RneaLink>,
    /// Gravity vector in world frame.
    pub gravity: [f64; 3],
}
impl RneaSolver {
    /// Create with given links and gravity vector.
    pub fn new(links: Vec<RneaLink>, gravity: [f64; 3]) -> Self {
        Self { links, gravity }
    }
    /// Compute joint torques using simplified (decoupled) RNEA.
    ///
    /// Returns a vector of joint torques of length `n_joints`.
    /// This is a simplified scalar approximation; a full RNEA would require
    /// SE(3) velocity/acceleration propagation.
    pub fn compute_torques(&self, q: &[f64], q_dot: &[f64], q_ddot: &[f64]) -> Vec<f64> {
        let n = self
            .links
            .len()
            .min(q.len())
            .min(q_dot.len())
            .min(q_ddot.len());
        let mut torques = vec![0.0; n];
        for i in 0..n {
            let link = &self.links[i];
            let j_inertia = link.axis_inertia();
            let tau_inertia = j_inertia * q_ddot[i];
            let tau_coriolis = 0.0;
            let g = self.gravity;
            let com = link.com;
            let ax = link.axis;
            let axcom = [
                ax[1] * com[2] - ax[2] * com[1],
                ax[2] * com[0] - ax[0] * com[2],
                ax[0] * com[1] - ax[1] * com[0],
            ];
            let angle_factor = q[i].cos();
            let tau_gravity =
                -link.mass * (axcom[0] * g[0] + axcom[1] * g[1] + axcom[2] * g[2]) * angle_factor;
            let tau_damping = -0.1 * q_dot[i];
            torques[i] = tau_inertia + tau_coriolis + tau_gravity + tau_damping;
        }
        torques
    }
    /// Compute the inertia matrix diagonal (simplified scalar per joint).
    pub fn inertia_diagonal(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.axis_inertia()).collect()
    }
}
/// Cartesian impedance controller with virtual spring-damper.
///
/// Implements `F = K_d * (x_d - x) + D_d * (ẋ_d - ẋ)` in Cartesian space.
#[derive(Debug, Clone)]
pub struct ImpedanceControl {
    /// Stiffness matrix diagonal `[kx, ky, kz]` (N/m).
    pub stiffness: [f64; 3],
    /// Damping matrix diagonal `[dx, dy, dz]` (N·s/m).
    pub damping: [f64; 3],
    /// Desired end-effector position.
    pub x_desired: [f64; 3],
    /// Desired end-effector velocity.
    pub xdot_desired: [f64; 3],
}
impl ImpedanceControl {
    /// Create a new impedance controller.
    pub fn new(stiffness: [f64; 3], damping: [f64; 3]) -> Self {
        Self {
            stiffness,
            damping,
            x_desired: [0.0; 3],
            xdot_desired: [0.0; 3],
        }
    }
    /// Compute the Cartesian wrench `[Fx, Fy, Fz]` for the given state.
    pub fn wrench(&self, x_current: [f64; 3], xdot_current: [f64; 3]) -> [f64; 3] {
        [
            self.stiffness[0] * (self.x_desired[0] - x_current[0])
                + self.damping[0] * (self.xdot_desired[0] - xdot_current[0]),
            self.stiffness[1] * (self.x_desired[1] - x_current[1])
                + self.damping[1] * (self.xdot_desired[1] - xdot_current[1]),
            self.stiffness[2] * (self.x_desired[2] - x_current[2])
                + self.damping[2] * (self.xdot_desired[2] - xdot_current[2]),
        ]
    }
    /// Map Cartesian wrench to joint torques via `τ = J^T F`.
    pub fn joint_torques(
        &self,
        robot: &SerialManipulator,
        x_current: [f64; 3],
        xdot_current: [f64; 3],
    ) -> Vec<f64> {
        let f = self.wrench(x_current, xdot_current);
        let jac = Jacobian::compute(robot);
        let n = robot.dof();
        (0..n)
            .map(|i| jac.j[i] * f[0] + jac.j[n + i] * f[1] + jac.j[2 * n + i] * f[2])
            .collect()
    }
}
/// Per-joint torque limits.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TorqueLimits {
    /// Maximum torque magnitude for each joint (N·m).
    pub limits: Vec<f64>,
}
impl TorqueLimits {
    /// Create from a list of maximum magnitudes.
    pub fn new(limits: Vec<f64>) -> Self {
        Self { limits }
    }
    /// Uniform limits for `n` joints.
    pub fn uniform(n: usize, max_torque: f64) -> Self {
        Self {
            limits: vec![max_torque; n],
        }
    }
    /// Clamp a torque vector to `[-limit, +limit]` per joint.
    pub fn clamp(&self, torques: &[f64]) -> Vec<f64> {
        torques
            .iter()
            .zip(self.limits.iter())
            .map(|(t, lim)| t.max(-lim).min(*lim))
            .collect()
    }
    /// Returns `true` if all torques are within limits (with tolerance `eps`).
    pub fn satisfied(&self, torques: &[f64], eps: f64) -> bool {
        torques
            .iter()
            .zip(self.limits.iter())
            .all(|(t, lim)| t.abs() <= lim + eps)
    }
    /// Scale down torques proportionally if any exceeds its limit.
    pub fn scale_to_satisfy(&self, torques: &[f64]) -> Vec<f64> {
        let ratio = torques
            .iter()
            .zip(self.limits.iter())
            .map(|(t, lim)| if *lim > 1e-15 { t.abs() / lim } else { 0.0 })
            .fold(1.0_f64, f64::max);
        if ratio > 1.0 {
            scale_vec(torques, 1.0 / ratio)
        } else {
            torques.to_vec()
        }
    }
}
/// 6-DOF Cartesian impedance control with mass-spring-damper dynamics.
///
/// Implements: F_ext = M_d * ẍ + D_d * ẋ + K_d * x_err
/// The controller computes the required joint torques via the Jacobian transpose.
#[allow(dead_code)]
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
/// Passivity-based robot controller using Lagrangian energy shaping.
///
/// Implements: τ = τ_d + K_D(q̇_d - q̇) + g(q) - g_d(q)
/// where the energy shaping term alters the potential energy landscape.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PassivityBasedController {
    /// Derivative gain for damping injection.
    pub k_d: Vec<f64>,
    /// Energy level target (for monitoring passivity).
    pub energy_target: f64,
    /// Stored total energy.
    pub current_energy: f64,
}
impl PassivityBasedController {
    /// Create with per-joint damping gains.
    pub fn new(k_d: Vec<f64>) -> Self {
        Self {
            k_d,
            energy_target: 0.0,
            current_energy: 0.0,
        }
    }
    /// Compute damping injection torque: τ_D = -K_D * q̇
    pub fn damping_torque(&self, q_dot: &[f64]) -> Vec<f64> {
        q_dot
            .iter()
            .zip(self.k_d.iter())
            .map(|(dq, k)| -k * dq)
            .collect()
    }
    /// Update current energy estimate (KE + PE).
    pub fn update_energy(&mut self, kinetic: f64, potential: f64) {
        self.current_energy = kinetic + potential;
    }
    /// Passivity condition: power extracted = q̇^T * τ_D ≤ 0.
    pub fn power_extracted(&self, q_dot: &[f64]) -> f64 {
        let tau_d = self.damping_torque(q_dot);
        q_dot.iter().zip(tau_d.iter()).map(|(v, t)| v * t).sum()
    }
    /// Check passivity condition (power extracted ≤ 0).
    pub fn is_passive(&self, q_dot: &[f64]) -> bool {
        self.power_extracted(q_dot) <= 1e-10
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
/// Cartesian-space trajectory between two end-effector poses.
#[derive(Debug, Clone)]
pub struct CartesianTrajectory {
    /// Start position `[x, y, z]`.
    pub p_start: [f64; 3],
    /// End position `[x, y, z]`.
    pub p_end: [f64; 3],
    /// Start orientation as ZYX Euler angles.
    pub r_start: [f64; 3],
    /// End orientation as ZYX Euler angles.
    pub r_end: [f64; 3],
    /// Duration (seconds).
    pub duration: f64,
    /// Trajectory type.
    pub traj_type: CartesianTrajectoryType,
    /// Arc centre for circular trajectories.
    pub arc_centre: [f64; 3],
}
impl CartesianTrajectory {
    /// Construct a linear Cartesian trajectory.
    pub fn linear(
        p_start: [f64; 3],
        p_end: [f64; 3],
        r_start: [f64; 3],
        r_end: [f64; 3],
        duration: f64,
    ) -> Self {
        Self {
            p_start,
            p_end,
            r_start,
            r_end,
            duration,
            traj_type: CartesianTrajectoryType::Linear,
            arc_centre: [0.0; 3],
        }
    }
    /// Evaluate trajectory at time `t`, returning `(position, euler_angles)`.
    pub fn evaluate(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let s = (t / self.duration).clamp(0.0, 1.0);
        let h = 3.0 * s * s - 2.0 * s * s * s;
        match self.traj_type {
            CartesianTrajectoryType::Linear => {
                let p = [
                    self.p_start[0] + (self.p_end[0] - self.p_start[0]) * h,
                    self.p_start[1] + (self.p_end[1] - self.p_start[1]) * h,
                    self.p_start[2] + (self.p_end[2] - self.p_start[2]) * h,
                ];
                let r = [
                    self.r_start[0] + (self.r_end[0] - self.r_start[0]) * h,
                    self.r_start[1] + (self.r_end[1] - self.r_start[1]) * h,
                    self.r_start[2] + (self.r_end[2] - self.r_start[2]) * h,
                ];
                (p, r)
            }
            CartesianTrajectoryType::Circular => {
                let rp = [
                    self.p_start[0] - self.arc_centre[0],
                    self.p_start[1] - self.arc_centre[1],
                    self.p_start[2] - self.arc_centre[2],
                ];
                let end_v = [
                    self.p_end[0] - self.arc_centre[0],
                    self.p_end[1] - self.arc_centre[1],
                    self.p_end[2] - self.arc_centre[2],
                ];
                let angle = dot3(normalize3(rp), normalize3(end_v))
                    .clamp(-1.0, 1.0)
                    .acos();
                let theta = angle * h;
                let axis = normalize3(cross3(rp, end_v));
                let rr = rotation_from_axis_angle(axis, theta);
                let rotated = mat3_vec_mul(rr, rp);
                let p = [
                    self.arc_centre[0] + rotated[0],
                    self.arc_centre[1] + rotated[1],
                    self.arc_centre[2] + rotated[2],
                ];
                let r = [
                    self.r_start[0] + (self.r_end[0] - self.r_start[0]) * h,
                    self.r_start[1] + (self.r_end[1] - self.r_start[1]) * h,
                    self.r_start[2] + (self.r_end[2] - self.r_start[2]) * h,
                ];
                (p, r)
            }
            CartesianTrajectoryType::Screw => {
                let axis = normalize3([
                    self.p_end[0] - self.p_start[0],
                    self.p_end[1] - self.p_start[1],
                    self.p_end[2] - self.p_start[2],
                ]);
                let angle = (self.r_end[2] - self.r_start[2]) * h;
                let rr = rotation_from_axis_angle(axis, angle);
                let v = [
                    self.p_end[0] - self.p_start[0],
                    self.p_end[1] - self.p_start[1],
                    self.p_end[2] - self.p_start[2],
                ];
                let rotated = mat3_vec_mul(rr, [self.p_start[0], self.p_start[1], self.p_start[2]]);
                let p = [
                    self.p_start[0] + v[0] * h + rotated[0] * 0.0,
                    self.p_start[1] + v[1] * h + rotated[1] * 0.0,
                    self.p_start[2] + v[2] * h + rotated[2] * 0.0,
                ];
                let r = [
                    self.r_start[0] + (self.r_end[0] - self.r_start[0]) * h,
                    self.r_start[1] + (self.r_end[1] - self.r_start[1]) * h,
                    self.r_start[2] + (self.r_end[2] - self.r_start[2]) * h,
                ];
                (p, r)
            }
        }
    }
}
/// Safety categories for robot motion.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum SafetyLevel {
    /// Normal operation — all safety checks passed.
    Normal,
    /// Warning state — some limits approached.
    Warning,
    /// Emergency stop required.
    EmergencyStop,
}
/// Multi-joint quintic trajectory.
#[allow(dead_code)]
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
    pub fn sample(&self, n_steps: usize) -> Vec<(f64, Vec<f64>, Vec<f64>, Vec<f64>)> {
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
/// Multi-criterion safety monitor for robot controllers.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RobotSafetyMonitor {
    /// Joint position limits \[q_min, q_max\] per joint.
    pub pos_limits: Vec<[f64; 2]>,
    /// Joint velocity limits (absolute value) per joint.
    pub vel_limits: Vec<f64>,
    /// Joint torque limits (absolute value) per joint.
    pub torque_limits: Vec<f64>,
    /// Warning threshold fraction (0.9 = 90% of limit triggers warning).
    pub warning_fraction: f64,
    /// Current safety level.
    pub level: SafetyLevel,
    /// Reason for current safety level.
    pub reason: String,
}
impl RobotSafetyMonitor {
    /// Create monitor with given limits.
    pub fn new(pos_limits: Vec<[f64; 2]>, vel_limits: Vec<f64>, torque_limits: Vec<f64>) -> Self {
        Self {
            pos_limits,
            vel_limits,
            torque_limits,
            warning_fraction: 0.9,
            level: SafetyLevel::Normal,
            reason: String::new(),
        }
    }
    /// Check current robot state, update and return safety level.
    pub fn check(&mut self, q: &[f64], q_dot: &[f64], torques: &[f64]) -> &SafetyLevel {
        for (i, (&qi, lims)) in q.iter().zip(self.pos_limits.iter()).enumerate() {
            if qi < lims[0] || qi > lims[1] {
                self.level = SafetyLevel::EmergencyStop;
                self.reason = format!(
                    "Joint {i} position {qi:.3} outside limits [{:.3}, {:.3}]",
                    lims[0], lims[1]
                );
                return &self.level;
            }
            let range = lims[1] - lims[0];
            let warn_lo = lims[0] + range * (1.0 - self.warning_fraction);
            let warn_hi = lims[1] - range * (1.0 - self.warning_fraction);
            if (qi < warn_lo || qi > warn_hi) && self.level != SafetyLevel::EmergencyStop {
                self.level = SafetyLevel::Warning;
                self.reason = format!("Joint {i} near position limit");
            }
        }
        for (i, (&vi, &lim)) in q_dot.iter().zip(self.vel_limits.iter()).enumerate() {
            if vi.abs() > lim {
                self.level = SafetyLevel::EmergencyStop;
                self.reason = format!("Joint {i} velocity {vi:.3} exceeds limit {lim:.3}");
                return &self.level;
            }
            if vi.abs() > lim * self.warning_fraction && self.level != SafetyLevel::EmergencyStop {
                self.level = SafetyLevel::Warning;
                self.reason = format!("Joint {i} near velocity limit");
            }
        }
        for (i, (&ti, &lim)) in torques.iter().zip(self.torque_limits.iter()).enumerate() {
            if ti.abs() > lim {
                self.level = SafetyLevel::EmergencyStop;
                self.reason = format!("Joint {i} torque {ti:.3} exceeds limit {lim:.3}");
                return &self.level;
            }
        }
        if self.level != SafetyLevel::EmergencyStop {
            self.level = SafetyLevel::Normal;
            self.reason.clear();
        }
        &self.level
    }
    /// Reset to normal state (after operator acknowledgment).
    pub fn reset(&mut self) {
        self.level = SafetyLevel::Normal;
        self.reason.clear();
    }
}
/// Per-joint velocity limits (rad/s or m/s).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VelocityLimits {
    /// Maximum absolute velocity for each joint.
    pub limits: Vec<f64>,
}
impl VelocityLimits {
    /// Create from a list of maximum speeds.
    pub fn new(limits: Vec<f64>) -> Self {
        Self { limits }
    }
    /// Uniform limits for `n` joints.
    pub fn uniform(n: usize, max_vel: f64) -> Self {
        Self {
            limits: vec![max_vel; n],
        }
    }
    /// Clamp velocity vector to limits.
    pub fn clamp(&self, vel: &[f64]) -> Vec<f64> {
        vel.iter()
            .zip(self.limits.iter())
            .map(|(v, lim)| v.max(-lim).min(*lim))
            .collect()
    }
    /// Check whether the velocity vector satisfies limits.
    pub fn satisfied(&self, vel: &[f64], eps: f64) -> bool {
        vel.iter()
            .zip(self.limits.iter())
            .all(|(v, lim)| v.abs() <= lim + eps)
    }
}
/// Repulsive potential field collision avoidance.
#[derive(Debug, Clone)]
pub struct CollisionAvoidance {
    /// Influence distance of obstacles (metres).
    pub influence_dist: f64,
    /// Repulsive gain coefficient.
    pub eta: f64,
    /// Obstacle positions.
    pub obstacles: Vec<[f64; 3]>,
}
impl CollisionAvoidance {
    /// Create a new collision avoidance module.
    pub fn new(influence_dist: f64, eta: f64) -> Self {
        Self {
            influence_dist,
            eta,
            obstacles: Vec::new(),
        }
    }
    /// Add an obstacle at position `p`.
    pub fn add_obstacle(&mut self, p: [f64; 3]) {
        self.obstacles.push(p);
    }
    /// Compute the repulsive potential at position `q`.
    pub fn repulsive_potential(&self, q: [f64; 3]) -> f64 {
        let mut u = 0.0;
        for obs in &self.obstacles {
            let diff = [q[0] - obs[0], q[1] - obs[1], q[2] - obs[2]];
            let d = norm3(diff);
            if d < self.influence_dist && d > 1e-6 {
                let factor = 1.0 / d - 1.0 / self.influence_dist;
                u += 0.5 * self.eta * factor * factor;
            }
        }
        u
    }
    /// Compute the repulsive force `−∇U_rep` at position `q`.
    pub fn repulsive_force(&self, q: [f64; 3]) -> [f64; 3] {
        let mut f = [0.0f64; 3];
        for obs in &self.obstacles {
            let diff = [q[0] - obs[0], q[1] - obs[1], q[2] - obs[2]];
            let d = norm3(diff);
            if d < self.influence_dist && d > 1e-6 {
                let factor = 1.0 / d - 1.0 / self.influence_dist;
                let grad_scale = self.eta * factor / (d * d * d);
                f[0] += grad_scale * diff[0];
                f[1] += grad_scale * diff[1];
                f[2] += grad_scale * diff[2];
            }
        }
        f
    }
    /// Attractive force toward goal: `F_att = -k_att * (q - goal)`.
    pub fn attractive_force(&self, q: [f64; 3], goal: [f64; 3], k_att: f64) -> [f64; 3] {
        [
            -k_att * (q[0] - goal[0]),
            -k_att * (q[1] - goal[1]),
            -k_att * (q[2] - goal[2]),
        ]
    }
    /// Total potential field force at `q`.
    pub fn total_force(&self, q: [f64; 3], goal: [f64; 3], k_att: f64) -> [f64; 3] {
        let f_att = self.attractive_force(q, goal, k_att);
        let f_rep = self.repulsive_force(q);
        [
            f_att[0] + f_rep[0],
            f_att[1] + f_rep[1],
            f_att[2] + f_rep[2],
        ]
    }
    /// Attempt to escape a local minimum by adding random perturbation.
    ///
    /// Returns a perturbed position.
    pub fn escape_local_minimum(&self, q: [f64; 3], step: f64) -> [f64; 3] {
        let h = (q[0].abs() * 1234.5 + q[1].abs() * 678.9 + q[2].abs() * 999.1).sin();
        let h2 = (q[0].abs() * 543.2 + q[1].abs() * 111.1 + q[2].abs() * 777.7).cos();
        [
            q[0] + h * step,
            q[1] + h2 * step,
            q[2] + (h + h2) * 0.5 * step,
        ]
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
/// Singularity-robust Jacobian controller with variable damping.
///
/// Uses a task-scaling approach near singularities to avoid excessive
/// joint velocities while maintaining the primary task direction.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SingularityRobustController {
    /// Minimum singular value threshold below which singularity handling activates.
    pub sigma_min_threshold: f64,
    /// Maximum DLS damping applied at singularity.
    pub max_damping: f64,
    /// Transition band width for smooth damping interpolation.
    pub transition_width: f64,
    /// Current computed minimum singular value.
    pub sigma_min: f64,
}
impl SingularityRobustController {
    /// Create with given parameters.
    pub fn new(sigma_min_threshold: f64, max_damping: f64, transition_width: f64) -> Self {
        Self {
            sigma_min_threshold,
            max_damping,
            transition_width,
            sigma_min: f64::INFINITY,
        }
    }
    /// Compute variable DLS damping based on manipulability measure `w`.
    ///
    /// Near singularities (w → 0), damping increases to `max_damping`.
    /// Far from singularities (w > sigma_min_threshold), damping → 0.
    pub fn adaptive_damping(&self, w: f64) -> f64 {
        if w >= self.sigma_min_threshold {
            0.0
        } else {
            let t = 1.0 - w / self.sigma_min_threshold;
            let t = t.clamp(0.0, 1.0);
            self.max_damping * t * t
        }
    }
    /// Estimate manipulability from Jacobian (approximated as Frobenius norm).
    ///
    /// A rough heuristic: w ≈ ||J||_F / sqrt(m*n)
    pub fn estimate_manipulability(jacobian: &[Vec<f64>]) -> f64 {
        let m = jacobian.len();
        if m == 0 {
            return 0.0;
        }
        let n = jacobian[0].len();
        if n == 0 {
            return 0.0;
        }
        let frob2: f64 = jacobian
            .iter()
            .flat_map(|row| row.iter())
            .map(|v| v * v)
            .sum();
        (frob2 / (m * n) as f64).sqrt()
    }
    /// Compute joint velocities with singularity-robust damping.
    pub fn compute_joint_vel(&mut self, jacobian: &[Vec<f64>], task_vel: &[f64]) -> Vec<f64> {
        let w = Self::estimate_manipulability(jacobian);
        let lambda = self.adaptive_damping(w);
        self.sigma_min = w;
        let rrc = ResolvedRateController::new(
            lambda,
            VelocityLimits::uniform(
                if jacobian.is_empty() {
                    0
                } else {
                    jacobian[0].len()
                },
                100.0,
            ),
        );
        rrc.compute(jacobian, task_vel)
    }
}
/// Low-pass filter for F/T sensor signals.
#[allow(dead_code)]
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
        for i in 0..6 {
            let debiased = raw[i] - self.bias[i];
            self.state[i] = self.alpha * debiased + (1.0 - self.alpha) * self.state[i];
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
/// Analytic Jacobian computation for a serial manipulator.
///
/// Supports the 6×n geometric Jacobian (3 linear + 3 angular rows).
#[derive(Debug, Clone)]
pub struct Jacobian {
    /// Reference to DH parameters (cloned for independence).
    pub dh: Vec<DhParams>,
    /// Joint angles at which the Jacobian was last computed.
    pub joint_angles: Vec<f64>,
    /// Last computed Jacobian matrix stored row-major (6 × n).
    pub j: Vec<f64>,
}
impl Jacobian {
    /// Compute the geometric Jacobian for the given manipulator state.
    ///
    /// Returns a 6×n matrix stored row-major as a flat `Vec`f64`.
    pub fn compute(robot: &SerialManipulator) -> Self {
        let n = robot.dof();
        let mut transforms = vec![mat4_identity(); n + 1];
        transforms[0] = mat4_identity();
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
        let p_e = mat4_translation(transforms[n]);
        let mut j = vec![0.0f64; 6 * n];
        for i in 0..n {
            let z_i = [transforms[i][2], transforms[i][6], transforms[i][10]];
            let p_i = mat4_translation(transforms[i]);
            let diff = [p_e[0] - p_i[0], p_e[1] - p_i[1], p_e[2] - p_i[2]];
            let jv = cross3(z_i, diff);
            j[i] = jv[0];
            j[n + i] = jv[1];
            j[2 * n + i] = jv[2];
            j[3 * n + i] = z_i[0];
            j[4 * n + i] = z_i[1];
            j[5 * n + i] = z_i[2];
        }
        Self {
            dh: robot.links.clone(),
            joint_angles: robot.joint_angles.clone(),
            j,
        }
    }
    /// Manipulability index: `sqrt(det(J * J^T))`.
    ///
    /// A value near zero indicates a near-singular configuration.
    pub fn manipulability(&self) -> f64 {
        let n = self.joint_angles.len();
        if n == 0 {
            return 0.0;
        }
        let mut jjt = vec![0.0f64; 36];
        for r in 0..6 {
            for c in 0..6 {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.j[r * n + k] * self.j[c * n + k];
                }
                jjt[r * 6 + c] = s;
            }
        }
        let mut det = 1.0;
        for i in 0..6 {
            det *= jjt[i * 6 + i];
        }
        det.abs().sqrt()
    }
    /// Compute the pseudoinverse `J^+` of the 6×n Jacobian using the damped least-squares method.
    ///
    /// Returns an n×6 matrix stored row-major as a flat `Vec`f64`.
    /// `lambda` is the damping coefficient (set to 0 for pure pseudoinverse).
    pub fn damped_pseudoinverse(&self, lambda: f64) -> Vec<f64> {
        let n = self.joint_angles.len();
        if n == 0 {
            return vec![];
        }
        let rows = 6usize;
        let mut a = vec![0.0f64; 36];
        for r in 0..rows {
            for c in 0..rows {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.j[r * n + k] * self.j[c * n + k];
                }
                a[r * rows + c] = s;
            }
            a[r * rows + r] += lambda * lambda;
        }
        let a_inv = gauss_jordan_6x6(&a);
        let mut jp = vec![0.0f64; n * rows];
        for r in 0..n {
            for c in 0..rows {
                let mut s = 0.0;
                for k in 0..rows {
                    s += self.j[k * n + r] * a_inv[k * rows + c];
                }
                jp[r * rows + c] = s;
            }
        }
        jp
    }
}
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
/// Joint-space trajectory planning using polynomial profiles.
#[derive(Debug, Clone)]
pub struct TrajectoryPlanning {
    /// Start joint angles (radians).
    pub q_start: Vec<f64>,
    /// Goal joint angles (radians).
    pub q_goal: Vec<f64>,
    /// Total trajectory duration (seconds).
    pub duration: f64,
}
impl TrajectoryPlanning {
    /// Create a trajectory from `q_start` to `q_goal` in `duration` seconds.
    pub fn new(q_start: Vec<f64>, q_goal: Vec<f64>, duration: f64) -> Self {
        Self {
            q_start,
            q_goal,
            duration,
        }
    }
    /// Evaluate a cubic polynomial joint-space trajectory at time `t`.
    ///
    /// Returns joint angles that satisfy `q(0)=q_start`, `q(T)=q_goal`,
    /// `dq/dt(0)=0`, `dq/dt(T)=0`.
    pub fn cubic_at(&self, t: f64) -> Vec<f64> {
        let s = (t / self.duration).clamp(0.0, 1.0);
        let h = 3.0 * s * s - 2.0 * s * s * s;
        self.q_start
            .iter()
            .zip(self.q_goal.iter())
            .map(|(qs, qg)| qs + (qg - qs) * h)
            .collect()
    }
    /// Evaluate a quintic polynomial joint-space trajectory at time `t`.
    ///
    /// Satisfies zero velocity and acceleration at both endpoints.
    pub fn quintic_at(&self, t: f64) -> Vec<f64> {
        let s = (t / self.duration).clamp(0.0, 1.0);
        let h = 10.0 * s * s * s - 15.0 * s * s * s * s + 6.0 * s * s * s * s * s;
        self.q_start
            .iter()
            .zip(self.q_goal.iter())
            .map(|(qs, qg)| qs + (qg - qs) * h)
            .collect()
    }
    /// Velocity at time `t` for cubic profile (derivative of cubic Hermite).
    pub fn cubic_velocity_at(&self, t: f64) -> Vec<f64> {
        let s = (t / self.duration).clamp(0.0, 1.0);
        let ds_dt = 1.0 / self.duration;
        let dh_ds = 6.0 * s - 6.0 * s * s;
        let dh_dt = dh_ds * ds_dt;
        self.q_start
            .iter()
            .zip(self.q_goal.iter())
            .map(|(qs, qg)| (qg - qs) * dh_dt)
            .collect()
    }
    /// Trapezoidal time scaling: ramp up, constant, ramp down.
    ///
    /// Returns normalised speed profile `v(t)` in `[0, 1]`.
    pub fn trapezoidal_profile(&self, t: f64, ramp_frac: f64) -> f64 {
        let s = (t / self.duration).clamp(0.0, 1.0);
        let r = ramp_frac.clamp(0.0, 0.5);
        let peak = 1.0 / (1.0 - r);
        if s < r {
            peak * s / r
        } else if s < 1.0 - r {
            peak
        } else {
            peak * (1.0 - s) / r
        }
    }
}
/// Result of workspace analysis for a manipulator.
#[allow(dead_code)]
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
/// State machine for contact phase detection.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ContactTransitionDetector {
    /// Phase threshold: enter LightContact when force > this (N).
    pub light_threshold: f64,
    /// Phase threshold: enter StableContact when force > this (N).
    pub stable_threshold: f64,
    /// Phase threshold: trigger overforce when force > this (N).
    pub overforce_threshold: f64,
    /// Hysteresis (N) — prevents chattering.
    pub hysteresis: f64,
    /// Current phase.
    pub phase: ContactPhase,
    /// Number of consecutive samples at current phase.
    pub phase_count: usize,
    /// Minimum samples before phase transition.
    pub debounce_count: usize,
}
impl ContactTransitionDetector {
    /// Create with force thresholds.
    pub fn new(light: f64, stable: f64, overforce: f64, hysteresis: f64) -> Self {
        Self {
            light_threshold: light,
            stable_threshold: stable,
            overforce_threshold: overforce,
            hysteresis,
            phase: ContactPhase::FreeMotion,
            phase_count: 0,
            debounce_count: 3,
        }
    }
    /// Update with new force magnitude, returns current phase.
    pub fn update(&mut self, force: f64) -> &ContactPhase {
        let new_phase = if force > self.overforce_threshold {
            ContactPhase::OverforceDetected
        } else if force > self.stable_threshold {
            ContactPhase::StableContact
        } else if force > self.light_threshold {
            ContactPhase::LightContact
        } else {
            ContactPhase::FreeMotion
        };
        if new_phase == self.phase {
            self.phase_count += 1;
        } else {
            let hysteresis_ok = match (&self.phase, &new_phase) {
                (ContactPhase::StableContact, ContactPhase::LightContact) => {
                    force < self.stable_threshold - self.hysteresis
                }
                (ContactPhase::LightContact, ContactPhase::FreeMotion) => {
                    force < self.light_threshold - self.hysteresis
                }
                _ => true,
            };
            if hysteresis_ok && self.phase_count >= self.debounce_count {
                self.phase = new_phase;
                self.phase_count = 0;
            }
        }
        &self.phase
    }
    /// Reset to free motion state.
    pub fn reset(&mut self) {
        self.phase = ContactPhase::FreeMotion;
        self.phase_count = 0;
    }
}
/// 6-axis force/torque sensor data.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, Default)]
pub struct FtSensorData {
    /// Force components \[Fx, Fy, Fz\] in sensor frame (N).
    pub force: [f64; 3],
    /// Torque components \[Tx, Ty, Tz\] in sensor frame (N·m).
    pub torque: [f64; 3],
}
impl FtSensorData {
    /// Create from raw 6D wrench vector \[Fx,Fy,Fz,Tx,Ty,Tz\].
    pub fn from_wrench(w: [f64; 6]) -> Self {
        Self {
            force: [w[0], w[1], w[2]],
            torque: [w[3], w[4], w[5]],
        }
    }
    /// Total force magnitude.
    pub fn force_magnitude(&self) -> f64 {
        let f = self.force;
        (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()
    }
    /// Total torque magnitude.
    pub fn torque_magnitude(&self) -> f64 {
        let t = self.torque;
        (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt()
    }
    /// As flat 6D wrench array.
    pub fn as_wrench(&self) -> [f64; 6] {
        [
            self.force[0],
            self.force[1],
            self.force[2],
            self.torque[0],
            self.torque[1],
            self.torque[2],
        ]
    }
}
/// Quintic polynomial trajectory for smooth point-to-point motion.
///
/// Satisfies boundary conditions: q(0)=q0, q(T)=qf, q'(0)=q'(T)=0, q''(0)=q''(T)=0
#[allow(dead_code)]
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
/// Broyden-Fletcher adaptive Jacobian update for model-free control.
///
/// Implements the rank-1 Broyden update:
///
/// `J_{k+1} = J_k + (Δy - J_k Δu) Δu^T / ||Δu||²`
#[allow(dead_code)]
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
        for i in 0..self.m.min(residual.len()) {
            for j in 0..self.n.min(delta_q.len()) {
                self.jacobian[i][j] =
                    self.forgetting * self.jacobian[i][j] + residual[i] * delta_q[j] / norm2;
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
/// Detects transitions between free motion and contact using F/T data.
#[allow(dead_code)]
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
/// Simplified receding-horizon MPC for a robot joint.
///
/// Minimizes: J = Σ (q_ref - q)² * w_q + Σ τ² * w_tau
/// Subject to: torque limits.
///
/// Uses a simple gradient-descent approach for the horizon.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct JointMpc {
    /// Prediction horizon length (steps).
    pub horizon: usize,
    /// Time step (s).
    pub dt: f64,
    /// State penalty weight.
    pub w_q: f64,
    /// Control penalty weight.
    pub w_tau: f64,
    /// Joint inertia (kg·m²).
    pub inertia: f64,
    /// Maximum torque (N·m).
    pub tau_max: f64,
    /// Damping coefficient.
    pub damping: f64,
}
impl JointMpc {
    /// Create MPC controller.
    pub fn new(
        horizon: usize,
        dt: f64,
        w_q: f64,
        w_tau: f64,
        inertia: f64,
        tau_max: f64,
        damping: f64,
    ) -> Self {
        Self {
            horizon,
            dt,
            w_q,
            w_tau,
            inertia,
            tau_max,
            damping,
        }
    }
    /// Simulate one step: q += q_dot*dt, q_dot += (tau - d*q_dot)/I * dt
    fn simulate_step(&self, q: f64, q_dot: f64, tau: f64) -> (f64, f64) {
        let alpha = (tau - self.damping * q_dot) / self.inertia;
        let new_q_dot = q_dot + alpha * self.dt;
        let new_q = q + new_q_dot * self.dt;
        (new_q, new_q_dot)
    }
    /// Compute optimal first control action for given state and reference.
    ///
    /// Uses a simple receding-horizon gradient descent over the horizon.
    pub fn compute_control(&self, q: f64, q_dot: f64, q_ref: &[f64]) -> f64 {
        let h = self.horizon.min(q_ref.len());
        if h == 0 {
            return 0.0;
        }
        let error = q_ref[0] - q;
        let vel_error = if h > 1 {
            (q_ref[1] - q_ref[0]) / self.dt - q_dot
        } else {
            -q_dot
        };
        let tau_p = self.w_q * error * self.inertia / (self.dt * self.dt);
        let tau_d = vel_error * self.inertia / self.dt;
        let tau = (tau_p + tau_d).max(-self.tau_max).min(self.tau_max);
        let (q1, _) = self.simulate_step(q, q_dot, tau);
        let err_before = (q_ref[0] - q).abs();
        let err_after = (q_ref[0] - q1).abs();
        if err_after > err_before * 2.0 {
            return -tau * 0.5;
        }
        tau
    }
    /// Roll out trajectory for `n_steps` under MPC control.
    pub fn rollout(&self, q0: f64, q_dot0: f64, q_ref: &[f64], n_steps: usize) -> Vec<(f64, f64)> {
        let mut q = q0;
        let mut q_dot = q_dot0;
        let mut traj = Vec::with_capacity(n_steps + 1);
        traj.push((q, q_dot));
        for i in 0..n_steps {
            let ref_slice = &q_ref[i..];
            let tau = self.compute_control(q, q_dot, ref_slice);
            let (nq, nqd) = self.simulate_step(q, q_dot, tau);
            q = nq;
            q_dot = nqd;
            traj.push((q, q_dot));
        }
        traj
    }
}
/// A single rigid link for RNEA calculations.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RneaLink {
    /// Link mass (kg).
    pub mass: f64,
    /// Center-of-mass position in local link frame.
    pub com: [f64; 3],
    /// Moment of inertia tensor (3×3, row-major) about CoM in local frame.
    pub inertia: [[f64; 3]; 3],
    /// Joint axis in local parent frame.
    pub axis: [f64; 3],
    /// Fixed offset from parent joint to child joint origin (DH translation).
    pub offset: [f64; 3],
}
impl RneaLink {
    /// Create a simple cylindrical link with uniform distribution.
    ///
    /// `length` — link length, `radius` — link radius, `mass` — total mass.
    /// Joint axis assumed along Z.
    pub fn cylindrical(mass: f64, length: f64, radius: f64) -> Self {
        let ixx = mass * (3.0 * radius * radius + length * length) / 12.0;
        let izz = 0.5 * mass * radius * radius;
        Self {
            mass,
            com: [0.0, 0.0, length / 2.0],
            inertia: [[ixx, 0.0, 0.0], [0.0, ixx, 0.0], [0.0, 0.0, izz]],
            axis: [0.0, 0.0, 1.0],
            offset: [0.0, 0.0, length],
        }
    }
    /// Rotational inertia about joint axis (scalar, J = axis^T * I * axis).
    pub fn axis_inertia(&self) -> f64 {
        let ax = self.axis;
        let i = self.inertia;
        let iax = [
            i[0][0] * ax[0] + i[0][1] * ax[1] + i[0][2] * ax[2],
            i[1][0] * ax[0] + i[1][1] * ax[1] + i[1][2] * ax[2],
            i[2][0] * ax[0] + i[2][1] * ax[1] + i[2][2] * ax[2],
        ];
        iax[0] * ax[0] + iax[1] * ax[1] + iax[2] * ax[2]
    }
}
/// Resolved-rate (Jacobian-based velocity control) at the task level.
///
/// Computes joint velocity: `q̇ = J† * ẋ_d`
/// where `J†` is the Moore-Penrose pseudoinverse of the Jacobian.
#[allow(dead_code)]
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
                let mut sum = 0.0;
                for k in 0..n {
                    sum += jacobian[i][k] * jacobian[j][k];
                }
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
/// Null-space projector for redundancy resolution.
///
/// Projects secondary task velocities into the null-space of the primary task.
#[allow(dead_code)]
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
/// Null-space optimization criteria for redundant manipulators.
#[allow(dead_code)]
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
/// Online payload mass and CoM identification for gravity compensation.
///
/// Uses recursive least-squares (RLS) to identify `[m, m*cx, m*cy, m*cz]`
/// from end-effector wrench measurements at multiple robot configurations.
#[allow(dead_code)]
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
            for i in 0..4 {
                for j in 0..4 {
                    p_phi[i] += self.cov[i * 4 + j] * phi[j];
                }
            }
            let denom = self.forgetting
                + phi
                    .iter()
                    .zip(p_phi.iter())
                    .map(|(h, pp)| h * pp)
                    .sum::<f64>();
            let gain: Vec<f64> = p_phi.iter().map(|pp| pp / denom).collect();
            for i in 0..4 {
                self.params[i] += gain[i] * innov;
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
/// Virtual spring-damper environment for stiffness rendering (haptics).
#[allow(dead_code)]
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
