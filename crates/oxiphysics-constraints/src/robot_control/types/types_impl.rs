//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::functions::*;

use super::types_2::{
    CartesianTrajectoryType, ContactPhase, DhParams, ResolvedRateController, SafetyLevel,
    SerialManipulator,
};

/// Simplified receding-horizon MPC for a robot joint.
///
/// Minimizes: J = Σ (q_ref - q)² * w_q + Σ τ² * w_tau
/// Subject to: torque limits.
///
/// Uses a simple gradient-descent approach for the horizon.
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
/// Multi-criterion safety monitor for robot controllers.
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
/// State machine for contact phase detection.
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
/// A single rigid link for RNEA calculations.
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
/// Per-joint torque limits.
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
/// Simplified RNEA for serial manipulators with revolute joints.
///
/// Computes joint torques from: `τ = M(q)q̈ + C(q,q̇)q̇ + G(q)`
/// using a reduced scalar approximation for each joint.
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
/// 6-axis force/torque sensor data.
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
/// Passivity-based robot controller using Lagrangian energy shaping.
///
/// Implements: τ = τ_d + K_D(q̇_d - q̇) + g(q) - g_d(q)
/// where the energy shaping term alters the potential energy landscape.
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
/// Singularity-robust Jacobian controller with variable damping.
///
/// Uses a task-scaling approach near singularities to avoid excessive
/// joint velocities while maintaining the primary task direction.
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
