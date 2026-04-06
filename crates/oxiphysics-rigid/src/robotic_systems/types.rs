//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// Workspace analysis: reachable workspace sampling and singularity detection.
pub struct WorkspaceAnalysis {
    /// Number of random samples to use.
    pub n_samples: usize,
}
impl WorkspaceAnalysis {
    /// Create a new workspace analyser with a given sample count.
    pub fn new(n_samples: usize) -> Self {
        Self { n_samples }
    }
    /// Sample the reachable workspace by random joint configurations.
    ///
    /// Returns a vector of end-effector positions `[x, y, z]`.
    pub fn sample_workspace(&self, chain: &KinematicChain) -> Vec<[f64; 3]> {
        let mut rng = rand::rng();
        let n = chain.dof();
        let mut samples = Vec::with_capacity(self.n_samples);
        let mut local_chain = KinematicChain {
            params: chain.params.clone(),
            joint_types: chain.joint_types.clone(),
            q: chain.q.clone(),
            limits: chain.limits.clone(),
        };
        for _ in 0..self.n_samples {
            for (i, (lo, hi)) in chain.limits.iter().enumerate().take(n) {
                local_chain.q[i] = rng.random_range(*lo..*hi);
            }
            samples.push(local_chain.end_effector_position());
        }
        samples
    }
    /// Compute the axis-aligned bounding box of the reachable workspace.
    ///
    /// Returns `(min_corner, max_corner)`.
    pub fn bounding_box(&self, chain: &KinematicChain) -> ([f64; 3], [f64; 3]) {
        let samples = self.sample_workspace(chain);
        if samples.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = samples[0];
        let mut hi = samples[0];
        for s in &samples {
            for i in 0..3 {
                if s[i] < lo[i] {
                    lo[i] = s[i];
                }
                if s[i] > hi[i] {
                    hi[i] = s[i];
                }
            }
        }
        (lo, hi)
    }
    /// Maximum reach distance from the base origin over sampled configurations.
    pub fn max_reach(&self, chain: &KinematicChain) -> f64 {
        self.sample_workspace(chain)
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }
    /// Detect near-singular configurations by sampling and checking manipulability.
    ///
    /// Returns configurations `q` whose manipulability is below `threshold`.
    pub fn detect_singularities(&self, chain: &KinematicChain, threshold: f64) -> Vec<Vec<f64>> {
        let mut rng = rand::rng();
        let n = chain.dof();
        let mut singular_configs = Vec::new();
        let mut local_chain = KinematicChain {
            params: chain.params.clone(),
            joint_types: chain.joint_types.clone(),
            q: chain.q.clone(),
            limits: chain.limits.clone(),
        };
        for _ in 0..self.n_samples {
            for (i, (lo, hi)) in chain.limits.iter().enumerate().take(n) {
                local_chain.q[i] = rng.random_range(*lo..*hi);
            }
            if local_chain.manipulability() < threshold {
                singular_configs.push(local_chain.q.clone());
            }
        }
        for q in &mut local_chain.q {
            *q = 0.0;
        }
        if local_chain.manipulability() < threshold {
            singular_configs.push(local_chain.q.clone());
        }
        singular_configs
    }
    /// Estimate the dexterous workspace volume (fraction of samples with
    /// manipulability above `min_manipulability`).
    pub fn dexterous_workspace_fraction(
        &self,
        chain: &KinematicChain,
        min_manipulability: f64,
    ) -> f64 {
        let mut rng = rand::rng();
        let n = chain.dof();
        let mut count = 0usize;
        let mut local_chain = KinematicChain {
            params: chain.params.clone(),
            joint_types: chain.joint_types.clone(),
            q: chain.q.clone(),
            limits: chain.limits.clone(),
        };
        for _ in 0..self.n_samples {
            for (i, (lo, hi)) in chain.limits.iter().enumerate().take(n) {
                local_chain.q[i] = rng.random_range(*lo..*hi);
            }
            if local_chain.manipulability() > min_manipulability {
                count += 1;
            }
        }
        count as f64 / self.n_samples.max(1) as f64
    }
}
/// Joint type for a kinematic chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Revolute: rotates about z-axis; theta is the variable.
    Revolute,
    /// Prismatic: translates along z-axis; d is the variable.
    Prismatic,
}
/// An open kinematic chain of rigid links connected by joints.
///
/// Implements DH-convention forward kinematics and geometric Jacobian.
pub struct KinematicChain {
    /// DH parameters for each link.
    pub params: Vec<DhParam>,
    /// Joint types (revolute or prismatic).
    pub joint_types: Vec<JointType>,
    /// Current joint configuration (angles or displacements).
    pub q: Vec<f64>,
    /// Joint limits: `(q_min, q_max)` for each joint.
    pub limits: Vec<(f64, f64)>,
}
impl KinematicChain {
    /// Create a new kinematic chain with given DH parameters and joint types.
    ///
    /// `q0` sets the initial configuration; limits default to `(-π, π)` for
    /// revolute and `(-1, 1)` for prismatic joints.
    pub fn new(params: Vec<DhParam>, joint_types: Vec<JointType>, q0: Vec<f64>) -> Self {
        let n = params.len();
        let limits = joint_types
            .iter()
            .map(|jt| match jt {
                JointType::Revolute => (-PI, PI),
                JointType::Prismatic => (-1.0, 1.0),
            })
            .collect();
        let q = if q0.len() == n { q0 } else { vec![0.0; n] };
        Self {
            params,
            joint_types,
            q,
            limits,
        }
    }
    /// Return the number of degrees of freedom.
    pub fn dof(&self) -> usize {
        self.params.len()
    }
    /// Compute the 4×4 homogeneous transform of joint `i` using the current
    /// configuration `q[i]`.
    pub fn link_transform(&self, i: usize) -> [[f64; 4]; 4] {
        let p = &self.params[i];
        let qi = self.q[i];
        let (theta, d) = match self.joint_types[i] {
            JointType::Revolute => (p.theta + qi, p.d),
            JointType::Prismatic => (p.theta, p.d + qi),
        };
        dh_transform(theta, d, p.a, p.alpha)
    }
    /// Compute the forward kinematics: returns a vector of transforms T_0^i
    /// for i = 1 .. n (cumulative products from base to each joint).
    pub fn forward_kinematics(&self) -> Vec<[[f64; 4]; 4]> {
        let n = self.dof();
        let mut transforms = Vec::with_capacity(n);
        let mut t = mat4_eye();
        for i in 0..n {
            let ti = self.link_transform(i);
            t = mat4_mul(&t, &ti);
            transforms.push(t);
        }
        transforms
    }
    /// Return the end-effector pose (last cumulative transform).
    pub fn end_effector_pose(&self) -> [[f64; 4]; 4] {
        self.forward_kinematics()
            .last()
            .copied()
            .unwrap_or_else(mat4_eye)
    }
    /// End-effector position `[x, y, z]`.
    pub fn end_effector_position(&self) -> [f64; 3] {
        let t = self.end_effector_pose();
        [t[0][3], t[1][3], t[2][3]]
    }
    /// Compute the 6×n geometric Jacobian J (rows: 3 linear, 3 angular).
    ///
    /// Returns a flat row-major vector of length `6 * n`.
    pub fn jacobian(&self) -> Vec<f64> {
        let n = self.dof();
        let transforms = self.forward_kinematics();
        let pe = if let Some(t) = transforms.last() {
            [t[0][3], t[1][3], t[2][3]]
        } else {
            [0.0; 3]
        };
        let mut j = vec![0.0; 6 * n];
        let mut t_prev = mat4_eye();
        for i in 0..n {
            let z = [t_prev[0][2], t_prev[1][2], t_prev[2][2]];
            let p = [t_prev[0][3], t_prev[1][3], t_prev[2][3]];
            let pe_minus_p = vec3_sub(pe, p);
            match self.joint_types[i] {
                JointType::Revolute => {
                    let jl = vec3_cross(z, pe_minus_p);
                    j[i] = jl[0];
                    j[n + i] = jl[1];
                    j[2 * n + i] = jl[2];
                    j[3 * n + i] = z[0];
                    j[4 * n + i] = z[1];
                    j[5 * n + i] = z[2];
                }
                JointType::Prismatic => {
                    j[i] = z[0];
                    j[n + i] = z[1];
                    j[2 * n + i] = z[2];
                    j[3 * n + i] = 0.0;
                    j[4 * n + i] = 0.0;
                    j[5 * n + i] = 0.0;
                }
            }
            t_prev = transforms[i];
        }
        j
    }
    /// Extract a 6×n Jacobian as a `Vec<Vec`f64`>` for use with IK solvers.
    pub fn jacobian_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.dof();
        let j = self.jacobian();
        (0..6)
            .map(|r| (0..n).map(|c| j[r * n + c]).collect())
            .collect()
    }
    /// Clamp joint configuration to limits.
    pub fn clamp_joints(&mut self) {
        for (i, q) in self.q.iter_mut().enumerate() {
            let (lo, hi) = self.limits[i];
            *q = q.clamp(lo, hi);
        }
    }
    /// Check if a joint configuration satisfies limits.
    pub fn within_limits(&self) -> bool {
        self.q
            .iter()
            .zip(self.limits.iter())
            .all(|(&qi, &(lo, hi))| qi >= lo && qi <= hi)
    }
    /// Compute manipulability measure: sqrt(det(J_sq * J_sq^T)).
    ///
    /// Uses a square submatrix formed by selecting the `n` rows with largest
    /// norms from the 6×n Jacobian, so the result is non-zero at generic
    /// configurations regardless of robot dimension.
    ///
    /// Uses the n×n square Jacobian built from the first `n` position rows
    /// (x, y, z), giving det=0 exactly at singular configurations.
    pub fn manipulability(&self) -> f64 {
        let n = self.dof();
        if n == 0 {
            return 0.0;
        }
        let j_mat = self.jacobian_matrix();
        let k = n.min(3).min(j_mat.len());
        let j_k: Vec<Vec<f64>> = j_mat[..k].to_vec();
        let jt = mat_transpose(&j_k);
        let jjt = mat_mul_dyn(&j_k, &jt);
        det_nxn(&jjt).abs().sqrt()
    }
}
/// Impedance controller for end-effector force control.
pub struct ForceControl {
    /// Impedance parameters.
    pub params: ImpedanceParams,
    /// Desired contact force (N) in each Cartesian direction.
    pub desired_force: [f64; 3],
    /// Current position error.
    pub pos_error: [f64; 3],
    /// Current velocity error.
    pub vel_error: [f64; 3],
}
impl ForceControl {
    /// Create a force controller with given impedance parameters.
    pub fn new(params: ImpedanceParams, desired_force: [f64; 3]) -> Self {
        Self {
            params,
            desired_force,
            pos_error: [0.0; 3],
            vel_error: [0.0; 3],
        }
    }
    /// Compute the impedance control law:
    /// `M * x_ddot + B * x_dot + K * x = F_ext - F_d`.
    ///
    /// Returns the commanded acceleration vector given current errors and
    /// external force `f_ext`.
    pub fn compute_acceleration(
        &self,
        pos_err: [f64; 3],
        vel_err: [f64; 3],
        f_ext: [f64; 3],
    ) -> [f64; 3] {
        if self.params.mass < 1e-300 {
            return [0.0; 3];
        }
        let mut acc = [0.0; 3];
        for i in 0..3 {
            let force_err = f_ext[i] - self.desired_force[i];
            acc[i] =
                (force_err - self.params.damping * vel_err[i] - self.params.stiffness * pos_err[i])
                    / self.params.mass;
        }
        acc
    }
    /// Admittance control law: compute velocity from force error.
    ///
    /// `x_dot = (F_ext - F_d) / B` (simplified first-order admittance).
    pub fn admittance_velocity(&self, f_ext: [f64; 3]) -> [f64; 3] {
        let b = self.params.damping;
        if b < 1e-300 {
            return [0.0; 3];
        }
        let mut vel = [0.0; 3];
        for i in 0..3 {
            vel[i] = (f_ext[i] - self.desired_force[i]) / b;
        }
        vel
    }
    /// PI force controller: compute position correction from force tracking error.
    ///
    /// `dx = kp * (Fd - Fe) + ki * integral_error` where integral is approximated
    /// as the current error times `dt`.
    pub fn pi_force_control(
        &self,
        f_ext: [f64; 3],
        integral_error: [f64; 3],
        kp: f64,
        ki: f64,
    ) -> [f64; 3] {
        let mut dx = [0.0; 3];
        for i in 0..3 {
            let err = self.desired_force[i] - f_ext[i];
            dx[i] = kp * err + ki * integral_error[i];
        }
        dx
    }
    /// Update position and velocity errors.
    pub fn update_errors(&mut self, pos_err: [f64; 3], vel_err: [f64; 3]) {
        self.pos_error = pos_err;
        self.vel_error = vel_err;
    }
}
/// Trajectory planning with various interpolation schemes.
pub struct TrajectoryPlanning {
    /// Waypoints defining the trajectory.
    pub waypoints: Vec<Waypoint>,
}
impl TrajectoryPlanning {
    /// Create a trajectory from a list of waypoints.
    pub fn new(waypoints: Vec<Waypoint>) -> Self {
        Self { waypoints }
    }
    /// Linear interpolation between waypoints.
    ///
    /// Returns the joint configuration at time `t`.
    pub fn linear_interp(&self, t: f64) -> Option<Vec<f64>> {
        let wps = &self.waypoints;
        if wps.is_empty() {
            return None;
        }
        if t <= wps[0].time {
            return Some(wps[0].q.clone());
        }
        if t >= wps.last().expect("collection should not be empty").time {
            return Some(
                wps.last()
                    .expect("collection should not be empty")
                    .q
                    .clone(),
            );
        }
        for i in 0..wps.len() - 1 {
            let t0 = wps[i].time;
            let t1 = wps[i + 1].time;
            if t >= t0 && t <= t1 {
                let tau = (t - t0) / (t1 - t0).max(1e-300);
                let n = wps[i].q.len().min(wps[i + 1].q.len());
                let q: Vec<f64> = (0..n)
                    .map(|j| wps[i].q[j] * (1.0 - tau) + wps[i + 1].q[j] * tau)
                    .collect();
                return Some(q);
            }
        }
        None
    }
    /// Cubic polynomial trajectory between two configurations.
    ///
    /// Boundary conditions: zero velocity at start and end.
    /// Returns `(position, velocity)` at normalized time `s ∈ [0, 1]`.
    pub fn cubic_poly(q0: f64, qf: f64, s: f64) -> (f64, f64) {
        let s = s.clamp(0.0, 1.0);
        let dq = qf - q0;
        let pos = q0 + (3.0 * s * s - 2.0 * s * s * s) * dq;
        let vel = (6.0 * s - 6.0 * s * s) * dq;
        (pos, vel)
    }
    /// Quintic polynomial trajectory between two configurations.
    ///
    /// Boundary conditions: zero velocity and acceleration at start and end.
    /// Returns `(position, velocity, acceleration)` at normalized time `s ∈ [0, 1]`.
    pub fn quintic_poly(q0: f64, qf: f64, s: f64) -> (f64, f64, f64) {
        let s = s.clamp(0.0, 1.0);
        let dq = qf - q0;
        let s2 = s * s;
        let s3 = s2 * s;
        let s4 = s3 * s;
        let s5 = s4 * s;
        let h = 10.0 * s3 - 15.0 * s4 + 6.0 * s5;
        let dh = 30.0 * s2 - 60.0 * s3 + 30.0 * s4;
        let ddh = 60.0 * s - 180.0 * s2 + 120.0 * s3;
        (q0 + h * dq, dh * dq, ddh * dq)
    }
    /// Trapezoidal velocity profile for a single joint.
    ///
    /// Computes the position at time `t` given start `q0`, end `qf`, total
    /// time `T`, and acceleration time `ta` (ramp-up = ramp-down phase).
    ///
    /// Returns `None` if `ta > T/2`.
    pub fn trapezoidal_profile(q0: f64, qf: f64, t: f64, total_t: f64, ta: f64) -> Option<f64> {
        if ta > total_t / 2.0 || total_t < 1e-300 {
            return None;
        }
        let dq = qf - q0;
        let v_max = dq / (total_t - ta);
        let t = t.clamp(0.0, total_t);
        let pos = if t < ta {
            q0 + 0.5 * v_max / ta * t * t
        } else if t < total_t - ta {
            q0 + v_max * (t - ta / 2.0)
        } else {
            let dt = total_t - t;
            qf - 0.5 * v_max / ta * dt * dt
        };
        Some(pos)
    }
    /// Minimum-jerk trajectory (Hogan 1984).
    ///
    /// Returns `(position, velocity, acceleration)` at normalized time `s ∈ [0, 1]`.
    pub fn minimum_jerk(q0: f64, qf: f64, s: f64) -> (f64, f64, f64) {
        let s = s.clamp(0.0, 1.0);
        let dq = qf - q0;
        let s3 = s * s * s;
        let s4 = s3 * s;
        let s5 = s4 * s;
        let h = 10.0 * s3 - 15.0 * s4 + 6.0 * s5;
        let dh = 30.0 * s * s - 60.0 * s3 + 30.0 * s4;
        let ddh = 60.0 * s - 180.0 * s * s + 120.0 * s3;
        (q0 + h * dq, dh * dq, ddh * dq)
    }
    /// Evaluate cubic spline interpolation for a full configuration at time `t`.
    pub fn cubic_spline(&self, t: f64) -> Option<Vec<f64>> {
        let wps = &self.waypoints;
        if wps.is_empty() {
            return None;
        }
        if t <= wps[0].time {
            return Some(wps[0].q.clone());
        }
        if t >= wps.last().expect("collection should not be empty").time {
            return Some(
                wps.last()
                    .expect("collection should not be empty")
                    .q
                    .clone(),
            );
        }
        for i in 0..wps.len() - 1 {
            let t0 = wps[i].time;
            let t1 = wps[i + 1].time;
            if t >= t0 && t <= t1 {
                let dt = (t1 - t0).max(1e-300);
                let s = (t - t0) / dt;
                let n = wps[i].q.len().min(wps[i + 1].q.len());
                let q: Vec<f64> = (0..n)
                    .map(|j| Self::cubic_poly(wps[i].q[j], wps[i + 1].q[j], s).0)
                    .collect();
                return Some(q);
            }
        }
        None
    }
    /// Evaluate quintic spline interpolation for a full configuration at time `t`.
    pub fn quintic_spline(&self, t: f64) -> Option<Vec<f64>> {
        let wps = &self.waypoints;
        if wps.is_empty() {
            return None;
        }
        if t <= wps[0].time {
            return Some(wps[0].q.clone());
        }
        if t >= wps.last().expect("collection should not be empty").time {
            return Some(
                wps.last()
                    .expect("collection should not be empty")
                    .q
                    .clone(),
            );
        }
        for i in 0..wps.len() - 1 {
            let t0 = wps[i].time;
            let t1 = wps[i + 1].time;
            if t >= t0 && t <= t1 {
                let dt = (t1 - t0).max(1e-300);
                let s = (t - t0) / dt;
                let n = wps[i].q.len().min(wps[i + 1].q.len());
                let q: Vec<f64> = (0..n)
                    .map(|j| Self::quintic_poly(wps[i].q[j], wps[i + 1].q[j], s).0)
                    .collect();
                return Some(q);
            }
        }
        None
    }
    /// Compute arc length of the trajectory in joint space.
    pub fn arc_length(&self, num_samples: usize) -> f64 {
        if self.waypoints.len() < 2 || num_samples < 2 {
            return 0.0;
        }
        let t_start = self
            .waypoints
            .first()
            .expect("collection should not be empty")
            .time;
        let t_end = self
            .waypoints
            .last()
            .expect("collection should not be empty")
            .time;
        let dt = (t_end - t_start) / (num_samples - 1) as f64;
        let mut prev = self.cubic_spline(t_start);
        let mut length = 0.0;
        for step in 1..num_samples {
            let t = t_start + step as f64 * dt;
            let curr = self.cubic_spline(t);
            if let (Some(q_prev), Some(q_curr)) = (&prev, &curr) {
                let n = q_prev.len().min(q_curr.len());
                let sq_dist: f64 = (0..n).map(|j| (q_curr[j] - q_prev[j]).powi(2)).sum();
                length += sq_dist.sqrt();
            }
            prev = curr;
        }
        length
    }
}
/// Grasp stability analysis: wrench space, force closure, and stability index.
pub struct GraspStability {
    /// Contact points on the object surface.
    pub contacts: Vec<ContactPoint>,
    /// Object centre of mass.
    pub object_com: [f64; 3],
}
impl GraspStability {
    /// Create a new grasp stability object.
    pub fn new(contacts: Vec<ContactPoint>, object_com: [f64; 3]) -> Self {
        Self {
            contacts,
            object_com,
        }
    }
    /// Test for force closure.
    ///
    /// A grasp has force closure if the contact wrenches can resist any
    /// external disturbance. This implementation uses the heuristic that
    /// the vector sum of contact normals should be small relative to the
    /// individual contact magnitudes.
    pub fn force_closure(&self) -> bool {
        if self.contacts.len() < 2 {
            return false;
        }
        let mut sum = [0.0_f64; 3];
        for c in &self.contacts {
            sum = [
                sum[0] + c.normal[0],
                sum[1] + c.normal[1],
                sum[2] + c.normal[2],
            ];
        }
        let sum_norm = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
        let n = self.contacts.len() as f64;
        sum_norm < n * 0.5
    }
    /// Compute the grasp wrench space as a flat list of 6D wrench vectors.
    ///
    /// Each contact contributes a set of wrenches spanning its friction cone.
    pub fn wrench_space(&self) -> Vec<[f64; 6]> {
        let mut wrenches = Vec::new();
        for contact in &self.contacts {
            let mut ws = contact.wrench_basis(self.object_com);
            wrenches.append(&mut ws);
        }
        wrenches
    }
    /// Stability index: minimum distance from the origin to the boundary of the
    /// convex hull of contact wrenches (heuristic).
    ///
    /// A larger value indicates a more stable grasp.
    pub fn stability_index(&self) -> f64 {
        let ws = self.wrench_space();
        if ws.is_empty() {
            return 0.0;
        }
        ws.iter()
            .map(|w| w.iter().map(|x| x * x).sum::<f64>().sqrt())
            .fold(f64::INFINITY, f64::min)
    }
    /// Compute the grasp matrix G (6 × 3n) where n = number of contacts.
    ///
    /// Delegates to [`GraspPlanning::grasp_matrix`] for compatibility.
    pub fn grasp_matrix(&self) -> Vec<f64> {
        let gp = GraspPlanning::new(self.contacts.clone(), self.object_com);
        gp.grasp_matrix()
    }
    /// Check if contact forces satisfy friction cone constraints.
    ///
    /// `forces` is a flat 3n vector.
    pub fn forces_in_cones(&self, forces: &[f64]) -> bool {
        let gp = GraspPlanning::new(self.contacts.clone(), self.object_com);
        gp.forces_within_cones(forces)
    }
    /// Number of contact points.
    pub fn n_contacts(&self) -> usize {
        self.contacts.len()
    }
}
/// A sphere obstacle in workspace.
#[derive(Debug, Clone)]
pub struct SphereObstacle {
    /// Centre position.
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
}
impl SphereObstacle {
    /// Create a new sphere obstacle.
    pub fn new(centre: [f64; 3], radius: f64) -> Self {
        Self { centre, radius }
    }
    /// Check if point `p` is inside the obstacle (with optional margin).
    pub fn contains(&self, p: [f64; 3], margin: f64) -> bool {
        let dist = vec3_norm(vec3_sub(p, self.centre));
        dist < self.radius + margin
    }
    /// Signed distance from point `p` to the surface (negative = inside).
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        vec3_norm(vec3_sub(p, self.centre)) - self.radius
    }
}
/// Impedance/admittance controller parameters.
#[derive(Debug, Clone)]
pub struct ImpedanceParams {
    /// Virtual mass (kg).
    pub mass: f64,
    /// Virtual damping (N·s/m).
    pub damping: f64,
    /// Virtual stiffness (N/m).
    pub stiffness: f64,
}
impl ImpedanceParams {
    /// Create new impedance parameters.
    pub fn new(mass: f64, damping: f64, stiffness: f64) -> Self {
        Self {
            mass,
            damping,
            stiffness,
        }
    }
    /// Compute the natural frequency ω_n = sqrt(k/m).
    pub fn natural_frequency(&self) -> f64 {
        if self.mass < 1e-300 {
            return 0.0;
        }
        (self.stiffness / self.mass).sqrt()
    }
    /// Compute the damping ratio ζ = c / (2 sqrt(k m)).
    pub fn damping_ratio(&self) -> f64 {
        let denom = 2.0 * (self.stiffness * self.mass).sqrt();
        if denom < 1e-300 {
            return 0.0;
        }
        self.damping / denom
    }
    /// Check if critically damped.
    pub fn is_critically_damped(&self) -> bool {
        (self.damping_ratio() - 1.0).abs() < 0.05
    }
}
/// Rigid-body robot dynamics using the recursive Newton-Euler algorithm.
pub struct RobotDynamics {
    /// Inertial parameters for each link.
    pub inertia: Vec<LinkInertia>,
    /// Gravity vector in base frame (m/s²).
    pub gravity: [f64; 3],
}
impl RobotDynamics {
    /// Create a new dynamics object with given link inertias and gravity.
    pub fn new(inertia: Vec<LinkInertia>, gravity: [f64; 3]) -> Self {
        Self { inertia, gravity }
    }
    /// Compute the generalized gravity torques for a given configuration.
    ///
    /// Uses a simplified gravity-projection method: τ_i = Σ_j (m_j * g^T * J_v^{ij})
    /// where J_v^{ij} is the linear Jacobian column for link j with respect to joint i.
    pub fn gravity_torques(&self, chain: &KinematicChain) -> Vec<f64> {
        let n = chain.dof();
        let mut tau = vec![0.0; n];
        let transforms = chain.forward_kinematics();
        for i in 0..n {
            let mut tau_i = 0.0;
            for j in i..n {
                if j >= self.inertia.len() {
                    break;
                }
                let m_j = self.inertia[j].mass;
                let t_prev = if i == 0 {
                    mat4_eye()
                } else {
                    transforms[i - 1]
                };
                let z = [t_prev[0][2], t_prev[1][2], t_prev[2][2]];
                let t_j = transforms[j];
                let com_j = self.inertia[j].com;
                let p_com = [
                    t_j[0][3] + t_j[0][0] * com_j[0] + t_j[0][1] * com_j[1] + t_j[0][2] * com_j[2],
                    t_j[1][3] + t_j[1][0] * com_j[0] + t_j[1][1] * com_j[1] + t_j[1][2] * com_j[2],
                    t_j[2][3] + t_j[2][0] * com_j[0] + t_j[2][1] * com_j[1] + t_j[2][2] * com_j[2],
                ];
                let p_i = if i == 0 {
                    [0.0; 3]
                } else {
                    let t_im1 = transforms[i - 1];
                    [t_im1[0][3], t_im1[1][3], t_im1[2][3]]
                };
                let r = vec3_sub(p_com, p_i);
                let z_cross_r = vec3_cross(z, r);
                tau_i -= m_j * vec3_dot(self.gravity, z_cross_r);
            }
            tau[i] = tau_i;
        }
        tau
    }
    /// Compute an approximate diagonal mass matrix entry for each joint.
    ///
    /// Returns the effective inertia seen at each joint (diagonal of M(q)).
    pub fn diagonal_mass_matrix(&self, chain: &KinematicChain) -> Vec<f64> {
        let n = chain.dof();
        let transforms = chain.forward_kinematics();
        let mut m_diag = vec![0.0; n];
        for i in 0..n {
            let t_prev = if i == 0 {
                mat4_eye()
            } else {
                transforms[i - 1]
            };
            let z = [t_prev[0][2], t_prev[1][2], t_prev[2][2]];
            let p_i = if i == 0 {
                [0.0; 3]
            } else {
                let t = transforms[i - 1];
                [t[0][3], t[1][3], t[2][3]]
            };
            for j in i..n {
                if j >= self.inertia.len() {
                    break;
                }
                let t_j = transforms[j];
                let com_j = self.inertia[j].com;
                let p_com = [
                    t_j[0][3] + t_j[0][0] * com_j[0] + t_j[0][1] * com_j[1] + t_j[0][2] * com_j[2],
                    t_j[1][3] + t_j[1][0] * com_j[0] + t_j[1][1] * com_j[1] + t_j[1][2] * com_j[2],
                    t_j[2][3] + t_j[2][0] * com_j[0] + t_j[2][1] * com_j[1] + t_j[2][2] * com_j[2],
                ];
                let r = vec3_sub(p_com, p_i);
                let z_cross_r = vec3_cross(z, r);
                let r2 = vec3_dot(z_cross_r, z_cross_r);
                m_diag[i] += self.inertia[j].mass * r2;
                if j < self.inertia.len() {
                    let rot = &self.inertia[j].inertia;
                    let iz = [
                        rot[0][0] * z[0] + rot[0][1] * z[1] + rot[0][2] * z[2],
                        rot[1][0] * z[0] + rot[1][1] * z[1] + rot[1][2] * z[2],
                        rot[2][0] * z[0] + rot[2][1] * z[1] + rot[2][2] * z[2],
                    ];
                    m_diag[i] += vec3_dot(z, iz);
                }
            }
        }
        m_diag
    }
    /// Compute Coriolis/centrifugal torques (simplified velocity-product terms).
    ///
    /// Returns `C(q, qdot) * qdot` using a finite-difference approximation of the
    /// mass matrix derivative.
    pub fn coriolis_torques(&self, chain: &KinematicChain, qdot: &[f64]) -> Vec<f64> {
        let n = chain.dof();
        let h = 1e-6;
        let mut c_qdot = vec![0.0; n];
        let m0 = self.diagonal_mass_matrix(chain);
        for j in 0..n {
            let mut chain_plus = KinematicChain {
                params: chain.params.clone(),
                joint_types: chain.joint_types.clone(),
                q: chain.q.clone(),
                limits: chain.limits.clone(),
            };
            if j < chain_plus.q.len() {
                chain_plus.q[j] += h;
            }
            let m_plus = self.diagonal_mass_matrix(&chain_plus);
            for i in 0..n {
                let dm_dqj = (m_plus[i] - m0[i]) / h;
                if j < qdot.len() && i < qdot.len() {
                    c_qdot[i] += dm_dqj * qdot[j] * qdot[j] * 0.5;
                }
            }
        }
        c_qdot
    }
}
/// Collision avoidance for a kinematic chain.
pub struct CollisionAvoidance {
    /// List of workspace obstacles.
    pub obstacles: Vec<SphereObstacle>,
    /// Minimum allowed distance between link positions (self-collision).
    pub self_collision_dist: f64,
    /// Safety margin for obstacle avoidance.
    pub safety_margin: f64,
}
impl CollisionAvoidance {
    /// Create a collision avoidance module.
    pub fn new(
        obstacles: Vec<SphereObstacle>,
        self_collision_dist: f64,
        safety_margin: f64,
    ) -> Self {
        Self {
            obstacles,
            self_collision_dist,
            safety_margin,
        }
    }
    /// Check if any link position collides with workspace obstacles.
    ///
    /// Returns `true` if a collision (within safety margin) is detected.
    pub fn check_obstacle_collision(&self, chain: &KinematicChain) -> bool {
        let transforms = chain.forward_kinematics();
        for t in &transforms {
            let pos = [t[0][3], t[1][3], t[2][3]];
            for obs in &self.obstacles {
                if obs.contains(pos, self.safety_margin) {
                    return true;
                }
            }
        }
        false
    }
    /// Check for self-collision by computing pairwise distances between
    /// non-adjacent link origins.
    ///
    /// Returns `true` if any pair is closer than `self_collision_dist`.
    pub fn check_self_collision(&self, chain: &KinematicChain) -> bool {
        let transforms = chain.forward_kinematics();
        let n = transforms.len();
        for i in 0..n {
            for j in (i + 2)..n {
                let pi = [
                    transforms[i][0][3],
                    transforms[i][1][3],
                    transforms[i][2][3],
                ];
                let pj = [
                    transforms[j][0][3],
                    transforms[j][1][3],
                    transforms[j][2][3],
                ];
                let dist = vec3_norm(vec3_sub(pi, pj));
                if dist < self.self_collision_dist {
                    return true;
                }
            }
        }
        false
    }
    /// Compute a repulsive gradient for a joint configuration due to obstacles.
    ///
    /// Returns a configuration-space gradient vector pushing away from obstacles.
    pub fn obstacle_gradient(&self, chain: &KinematicChain) -> Vec<f64> {
        let n = chain.dof();
        let mut grad = vec![0.0; n];
        let transforms = chain.forward_kinematics();
        let j_mat = chain.jacobian_matrix();
        for (link_idx, t) in transforms.iter().enumerate() {
            let pos = [t[0][3], t[1][3], t[2][3]];
            for obs in &self.obstacles {
                let d = obs.signed_distance(pos);
                let d_thresh = self.safety_margin + obs.radius * 0.5;
                if d < d_thresh && d > 1e-6 {
                    let dp = vec3_sub(pos, obs.centre);
                    let dp_norm = vec3_norm(dp);
                    if dp_norm < 1e-10 {
                        continue;
                    }
                    let n_obs = vec3_scale(dp, 1.0 / dp_norm);
                    let k = 1.0;
                    let d0 = d_thresh;
                    let grad_task = vec3_scale(n_obs, k * (1.0 / d - 1.0 / d0) / (d * d));
                    for i in 0..n {
                        for dim in 0..3 {
                            if dim < j_mat.len() && i < j_mat[dim].len() {
                                let _ = link_idx;
                                grad[i] += j_mat[dim][i] * grad_task[dim];
                            }
                        }
                    }
                }
            }
        }
        grad
    }
    /// Return the minimum signed distance to all obstacles from all link positions.
    pub fn min_obstacle_distance(&self, chain: &KinematicChain) -> f64 {
        let transforms = chain.forward_kinematics();
        let mut min_dist = f64::INFINITY;
        for t in &transforms {
            let pos = [t[0][3], t[1][3], t[2][3]];
            for obs in &self.obstacles {
                let d = obs.signed_distance(pos);
                if d < min_dist {
                    min_dist = d;
                }
            }
        }
        min_dist
    }
}
/// A contact point on an object surface.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// Contact position in object frame.
    pub position: [f64; 3],
    /// Inward surface normal at contact.
    pub normal: [f64; 3],
    /// Friction coefficient at contact.
    pub mu: f64,
}
impl ContactPoint {
    /// Create a new contact point.
    pub fn new(position: [f64; 3], normal: [f64; 3], mu: f64) -> Self {
        let n = vec3_norm(normal);
        let normal = if n > 1e-10 {
            [normal[0] / n, normal[1] / n, normal[2] / n]
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            position,
            normal,
            mu,
        }
    }
    /// Compute the friction cone half-angle: atan(mu).
    pub fn cone_half_angle(&self) -> f64 {
        self.mu.atan()
    }
    /// Wrench basis for a friction contact (normal + 4 tangential directions).
    ///
    /// Returns 6 wrench vectors (force + torque) for linearized friction cone.
    pub fn wrench_basis(&self, object_com: [f64; 3]) -> Vec<[f64; 6]> {
        let r = vec3_sub(self.position, object_com);
        let t1 = tangent_vec(self.normal);
        let t2 = vec3_cross(self.normal, t1);
        let directions: Vec<[f64; 3]> = vec![
            self.normal,
            vec3_add(self.normal, vec3_scale(t1, self.mu)),
            vec3_sub(self.normal, vec3_scale(t1, self.mu)),
            vec3_add(self.normal, vec3_scale(t2, self.mu)),
            vec3_sub(self.normal, vec3_scale(t2, self.mu)),
        ];
        directions
            .iter()
            .map(|&f| {
                let torque = vec3_cross(r, f);
                [f[0], f[1], f[2], torque[0], torque[1], torque[2]]
            })
            .collect()
    }
}
/// Inverse kinematics solver using iterative Jacobian methods.
pub struct InverseKinematics {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance on position error.
    pub tol: f64,
    /// Damping coefficient for damped least squares.
    pub lambda: f64,
    /// Step size (gain) for each update.
    pub step: f64,
}
impl InverseKinematics {
    /// Create a new IK solver with default parameters.
    pub fn new() -> Self {
        Self {
            max_iter: 200,
            tol: 1e-4,
            lambda: 0.5,
            step: 0.1,
        }
    }
    /// Create an IK solver with custom parameters.
    pub fn with_params(max_iter: usize, tol: f64, lambda: f64, step: f64) -> Self {
        Self {
            max_iter,
            tol,
            lambda,
            step,
        }
    }
    /// Solve IK using the Jacobian pseudoinverse (Moore-Penrose).
    ///
    /// Solves for joint angles `q` such that the end-effector reaches
    /// position `target`. Returns `(q_solution, converged)`.
    ///
    /// Uses damped least squares `Δq = J^T (J J^T + λ²I)^{-1} e` which
    /// provides fast convergence while remaining robust near singularities.
    pub fn solve_pseudoinverse(
        &self,
        chain: &mut KinematicChain,
        target: [f64; 3],
    ) -> (Vec<f64>, bool) {
        let n = chain.dof();
        for _ in 0..self.max_iter {
            let pos = chain.end_effector_position();
            let err = vec3_sub(target, pos);
            let err_norm = vec3_norm(err);
            if err_norm < self.tol {
                return (chain.q.clone(), true);
            }
            let j_mat = chain.jacobian_matrix();
            let jp: Vec<Vec<f64>> = j_mat[..3.min(j_mat.len())].to_vec();
            let jt = mat_transpose(&jp);
            let err_vec: Vec<f64> = err.to_vec();
            let mut a = mat_mul_dyn(&jp, &jt);
            let lam2 = self.lambda * self.lambda;
            mat_add_lambda_eye(&mut a, lam2);
            if let Some(alpha) = solve_linear(&a, &err_vec) {
                let dq = mat_vec_mul(&jt, &alpha);
                for i in 0..n.min(dq.len()) {
                    chain.q[i] += self.step * dq[i];
                }
                chain.clamp_joints();
            } else {
                let dq = mat_vec_mul(&jt, &err_vec);
                for i in 0..n.min(dq.len()) {
                    chain.q[i] += self.step * 0.01 * dq[i];
                }
                chain.clamp_joints();
            }
        }
        let pos = chain.end_effector_position();
        let err = vec3_sub(target, pos);
        let converged = vec3_norm(err) < self.tol;
        (chain.q.clone(), converged)
    }
    /// Solve IK using damped least squares (Levenberg-Marquardt style).
    ///
    /// Uses `Δq = J^T (J J^T + λ²I)^{-1} e` in the task-space (m×m solve)
    /// which is numerically stable and handles near-singular Jacobians.
    pub fn solve_damped_ls(
        &self,
        chain: &mut KinematicChain,
        target: [f64; 3],
    ) -> (Vec<f64>, bool) {
        let n = chain.dof();
        for _ in 0..self.max_iter {
            let pos = chain.end_effector_position();
            let err = vec3_sub(target, pos);
            let err_norm = vec3_norm(err);
            if err_norm < self.tol {
                return (chain.q.clone(), true);
            }
            let j_mat = chain.jacobian_matrix();
            let jp: Vec<Vec<f64>> = j_mat[..3.min(j_mat.len())].to_vec();
            let jt = mat_transpose(&jp);
            let err_vec: Vec<f64> = err.to_vec();
            let mut a = mat_mul_dyn(&jp, &jt);
            let lam2 = self.lambda * self.lambda;
            mat_add_lambda_eye(&mut a, lam2);
            if let Some(alpha) = solve_linear(&a, &err_vec) {
                let dq = mat_vec_mul(&jt, &alpha);
                for i in 0..n.min(dq.len()) {
                    chain.q[i] += self.step * dq[i];
                }
                chain.clamp_joints();
            } else {
                let dq = mat_vec_mul(&jt, &err_vec);
                for i in 0..n.min(dq.len()) {
                    chain.q[i] += self.step * 0.01 * dq[i];
                }
                chain.clamp_joints();
            }
        }
        let pos = chain.end_effector_position();
        let err = vec3_sub(target, pos);
        let converged = vec3_norm(err) < self.tol;
        (chain.q.clone(), converged)
    }
    /// Null-space projection: compute a configuration update in the null space
    /// of the Jacobian so secondary objectives can be pursued without affecting
    /// the primary end-effector task.
    ///
    /// `secondary_grad` is the gradient of a secondary objective (e.g. joint
    /// centering). Returns the projected update `(I - J^+ J) * secondary_grad`.
    pub fn null_space_projection(
        &self,
        chain: &KinematicChain,
        secondary_grad: &[f64],
    ) -> Vec<f64> {
        let n = chain.dof();
        if secondary_grad.len() != n {
            return vec![0.0; n];
        }
        let j_mat = chain.jacobian_matrix();
        let jp: Vec<Vec<f64>> = j_mat[..3.min(j_mat.len())].to_vec();
        let jt = mat_transpose(&jp);
        let jjt = mat_mul_dyn(&jp, &jt);
        let jpinv = if let Some(_inv) = self.jjt_inv(&jjt) {
            mat_mul_dyn(&jt, &_inv)
        } else {
            return vec![0.0; n];
        };
        let jpj = mat_mul_dyn(&jpinv, &jp);
        let mut result = vec![0.0; n];
        for i in 0..n {
            result[i] = secondary_grad[i];
            for j in 0..n {
                if i < jpj.len() && j < jpj[i].len() {
                    result[i] -= jpj[i][j] * secondary_grad[j];
                }
            }
        }
        result
    }
    /// Compute the inverse of a square matrix (rows = cols = m).
    fn jjt_inv(&self, m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
        let n = m.len();
        if n == 0 {
            return Some(Vec::new());
        }
        let mut aug: Vec<Vec<f64>> = m
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let mut r = row.clone();
                for j in 0..n {
                    r.push(if i == j { 1.0 } else { 0.0 });
                }
                r
            })
            .collect();
        for col in 0..n {
            let pivot = (col..n).max_by(|&i, &j| {
                aug[i][col]
                    .abs()
                    .partial_cmp(&aug[j][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })?;
            aug.swap(col, pivot);
            let diag = aug[col][col];
            if diag.abs() < 1e-14 {
                return None;
            }
            for k in 0..2 * n {
                aug[col][k] /= diag;
            }
            for r in 0..n {
                if r != col {
                    let factor = aug[r][col];
                    for k in 0..2 * n {
                        let val = aug[col][k];
                        aug[r][k] -= factor * val;
                    }
                }
            }
        }
        Some(aug.iter().map(|row| row[n..].to_vec()).collect())
    }
}
/// Grasp planning and quality metrics.
pub struct GraspPlanning {
    /// Contact points of the grasp.
    pub contacts: Vec<ContactPoint>,
    /// Object centre of mass.
    pub object_com: [f64; 3],
}
impl GraspPlanning {
    /// Create a grasp planning object.
    pub fn new(contacts: Vec<ContactPoint>, object_com: [f64; 3]) -> Self {
        Self {
            contacts,
            object_com,
        }
    }
    /// Check force closure: test if the origin is in the interior of the convex
    /// hull of contact wrenches (simplified 6D test using convex combination).
    ///
    /// Returns `true` if the grasp has force closure (heuristic check).
    pub fn has_force_closure(&self) -> bool {
        if self.contacts.len() < 2 {
            return false;
        }
        let mut force_sum = [0.0f64; 3];
        for c in &self.contacts {
            force_sum = vec3_add(force_sum, c.normal);
        }
        let sum_norm = vec3_norm(force_sum);
        let avg_contact_dist: f64 = self
            .contacts
            .iter()
            .map(|c| vec3_norm(vec3_sub(c.position, self.object_com)))
            .sum::<f64>()
            / self.contacts.len() as f64;
        sum_norm < avg_contact_dist * (self.contacts.len() as f64)
    }
    /// Compute the grasp matrix G (6×3n) mapping contact forces to object wrench.
    ///
    /// Returns a flat row-major vector of size 6 × (3 * n_contacts).
    pub fn grasp_matrix(&self) -> Vec<f64> {
        let n = self.contacts.len();
        let cols = 3 * n;
        let mut g = vec![0.0f64; 6 * cols];
        for (ci, contact) in self.contacts.iter().enumerate() {
            let r = vec3_sub(contact.position, self.object_com);
            for d in 0..3 {
                let mut f = [0.0f64; 3];
                f[d] = 1.0;
                let tau = vec3_cross(r, f);
                let col = ci * 3 + d;
                g[col] = f[0];
                g[cols + col] = f[1];
                g[2 * cols + col] = f[2];
                g[3 * cols + col] = tau[0];
                g[4 * cols + col] = tau[1];
                g[5 * cols + col] = tau[2];
            }
        }
        g
    }
    /// Compute the Q1 grasp quality metric (smallest singular value of G).
    ///
    /// Uses power iteration to estimate the minimum singular value.
    /// A larger value indicates a higher quality grasp.
    pub fn quality_q1(&self) -> f64 {
        if self.contacts.is_empty() {
            return 0.0;
        }
        let n = self.contacts.len();
        let cols = 3 * n;
        let g = self.grasp_matrix();
        let mut ggt = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut s = 0.0;
                for k in 0..cols {
                    s += g[i * cols + k] * g[j * cols + k];
                }
                ggt[i][j] = s;
            }
        }
        let max_eig = self.power_iteration_max_eig(&ggt, 50);
        self.min_eigenvalue_approx(&ggt, max_eig)
    }
    /// Power iteration to estimate the maximum eigenvalue of a 6×6 matrix.
    fn power_iteration_max_eig(&self, m: &[[f64; 6]; 6], iters: usize) -> f64 {
        let mut v = [1.0f64; 6];
        let mut eig = 0.0;
        for _ in 0..iters {
            let mut mv = [0.0f64; 6];
            for i in 0..6 {
                for j in 0..6 {
                    mv[i] += m[i][j] * v[j];
                }
            }
            eig = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
            if eig < 1e-300 {
                break;
            }
            v = mv.map(|x| x / eig);
        }
        eig
    }
    /// Estimate the minimum eigenvalue using deflation.
    fn min_eigenvalue_approx(&self, m: &[[f64; 6]; 6], max_eig: f64) -> f64 {
        let mut shifted = *m;
        for i in 0..6 {
            shifted[i][i] = max_eig - m[i][i];
            for j in 0..6 {
                if i != j {
                    shifted[i][j] = -m[i][j];
                }
            }
        }
        let max_shifted = self.power_iteration_max_eig(&shifted, 50);
        (max_eig - max_shifted).max(0.0).sqrt()
    }
    /// Compute the contact force distribution minimizing total force magnitude.
    ///
    /// Solves the minimum-norm least-squares problem G * f = w_ext.
    /// Returns contact forces (3n vector) for the given external wrench `w_ext` (6 vector).
    pub fn min_norm_forces(&self, w_ext: [f64; 6]) -> Vec<f64> {
        let n = self.contacts.len();
        if n == 0 {
            return Vec::new();
        }
        let cols = 3 * n;
        let g = self.grasp_matrix();
        let g_mat: Vec<Vec<f64>> = (0..6)
            .map(|i| (0..cols).map(|j| g[i * cols + j]).collect())
            .collect();
        let gt = mat_transpose(&g_mat);
        let ggt = mat_mul_dyn(&g_mat, &gt);
        let w: Vec<f64> = w_ext.to_vec();
        if let Some(alpha) = solve_linear(&ggt, &w) {
            mat_vec_mul(&gt, &alpha)
        } else {
            vec![0.0; cols]
        }
    }
    /// Check if all contact forces are within friction cones.
    ///
    /// `forces` is a flat (3n) vector of contact forces.
    pub fn forces_within_cones(&self, forces: &[f64]) -> bool {
        let _n = self.contacts.len();
        for (i, contact) in self.contacts.iter().enumerate() {
            if i * 3 + 2 >= forces.len() {
                break;
            }
            let f = [forces[i * 3], forces[i * 3 + 1], forces[i * 3 + 2]];
            let fn_val = vec3_dot(f, contact.normal);
            if fn_val < 0.0 {
                return false;
            }
            let ft_mag = {
                let ft = vec3_sub(f, vec3_scale(contact.normal, fn_val));
                vec3_norm(ft)
            };
            if ft_mag > contact.mu * fn_val + 1e-8 {
                return false;
            }
        }
        true
    }
}
/// Geometric Jacobian matrix for velocity kinematics analysis.
///
/// Stores the 6×n Jacobian as a `Vec<Vec`f64`>` and provides
/// rank, condition number, and singular-value approximation.
pub struct JacobianMatrix {
    /// The 6×n Jacobian matrix (rows = \[vx, vy, vz, wx, wy, wz\]).
    pub data: Vec<Vec<f64>>,
    /// Number of joints (columns).
    pub n_joints: usize,
}
impl JacobianMatrix {
    /// Build a Jacobian matrix from a kinematic chain.
    pub fn from_chain(chain: &KinematicChain) -> Self {
        let n = chain.dof();
        let data = chain.jacobian_matrix();
        Self { data, n_joints: n }
    }
    /// Create from raw data (6 rows × n_joints columns).
    pub fn from_data(data: Vec<Vec<f64>>) -> Self {
        let n = data.first().map(|r| r.len()).unwrap_or(0);
        Self { data, n_joints: n }
    }
    /// Approximate numerical rank of the position sub-Jacobian (rows 0..3).
    ///
    /// Uses Gaussian elimination with a tolerance to detect near-zero pivots.
    pub fn rank(&self, tol: f64) -> usize {
        if self.data.is_empty() {
            return 0;
        }
        let rows: Vec<Vec<f64>> = self.data[..3.min(self.data.len())].to_vec();
        if rows.is_empty() {
            return 0;
        }
        let mut mat = rows;
        let n_rows = mat.len();
        let n_cols = self.n_joints;
        let mut rank = 0;
        let mut pivot_col = 0;
        for row in 0..n_rows {
            while pivot_col < n_cols {
                let max_val = mat[row..]
                    .iter()
                    .filter_map(|r| r.get(pivot_col))
                    .map(|v| v.abs())
                    .fold(0.0_f64, f64::max);
                if max_val > tol {
                    break;
                }
                pivot_col += 1;
            }
            if pivot_col >= n_cols {
                break;
            }
            let (max_row, _) = mat[row..]
                .iter()
                .enumerate()
                .map(|(i, r)| (i + row, r.get(pivot_col).copied().unwrap_or(0.0).abs()))
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or((row, 0.0));
            mat.swap(row, max_row);
            let pivot = mat[row].get(pivot_col).copied().unwrap_or(0.0);
            if pivot.abs() < tol {
                break;
            }
            for k in row + 1..n_rows {
                let factor = mat[k].get(pivot_col).copied().unwrap_or(0.0) / pivot;
                for c in 0..n_cols {
                    let val = mat[row].get(c).copied().unwrap_or(0.0);
                    if let Some(v) = mat[k].get_mut(c) {
                        *v -= factor * val;
                    }
                }
            }
            rank += 1;
            pivot_col += 1;
        }
        rank
    }
    /// Approximate condition number σ_max / σ_min of the position sub-Jacobian.
    ///
    /// Returns 0.0 if the Jacobian has zero columns.
    pub fn condition_number(&self) -> f64 {
        if self.n_joints == 0 {
            return 0.0;
        }
        let rows: Vec<Vec<f64>> = self.data[..3.min(self.data.len())].to_vec();
        let jt = mat_transpose(&rows);
        let jjt = mat_mul_dyn(&rows, &jt);
        let n = jjt.len();
        if n == 0 {
            return 0.0;
        }
        let mut v_max = vec![1.0_f64; n];
        let mut sigma_max_sq = 0.0_f64;
        for _ in 0..60 {
            let mv: Vec<f64> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            jjt.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0) * v_max[j]
                        })
                        .sum()
                })
                .collect();
            sigma_max_sq = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
            if sigma_max_sq < 1e-300 {
                break;
            }
            v_max = mv.iter().map(|x| x / sigma_max_sq).collect();
        }
        let shift = sigma_max_sq;
        let mut v_min = vec![1.0_f64; n];
        let mut sigma_min_sq = 0.0_f64;
        for _ in 0..60 {
            let mv: Vec<f64> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            let m_ij = jjt.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                            let s_ij = if i == j { shift } else { 0.0 };
                            (s_ij - m_ij) * v_min[j]
                        })
                        .sum()
                })
                .collect();
            let mv_norm = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
            if mv_norm < 1e-300 {
                break;
            }
            sigma_min_sq = mv_norm;
            v_min = mv.iter().map(|x| x / mv_norm).collect();
        }
        let sigma_min = (shift - sigma_min_sq).max(0.0).sqrt();
        let sigma_max = sigma_max_sq.sqrt();
        if sigma_min < 1e-12 {
            return 0.0;
        }
        sigma_max / sigma_min
    }
    /// Linear velocity Jacobian: the first 3 rows (Jv).
    pub fn linear_part(&self) -> Vec<Vec<f64>> {
        self.data[..3.min(self.data.len())].to_vec()
    }
    /// Angular velocity Jacobian: the last 3 rows (Jω).
    pub fn angular_part(&self) -> Vec<Vec<f64>> {
        let start = 3.min(self.data.len());
        self.data[start..].to_vec()
    }
    /// Manipulability measure: sqrt(det(Jv * Jv^T)).
    pub fn manipulability(&self) -> f64 {
        let jv = self.linear_part();
        if jv.is_empty() {
            return 0.0;
        }
        let jt = mat_transpose(&jv);
        let jjt = mat_mul_dyn(&jv, &jt);
        det_nxn(&jjt).abs().sqrt()
    }
}
/// Denavit-Hartenberg parameters for a single joint.
#[derive(Debug, Clone, Copy)]
pub struct DhParam {
    /// Link length along x-axis (metres).
    pub a: f64,
    /// Link twist about x-axis (radians).
    pub alpha: f64,
    /// Link offset along z-axis (metres).
    pub d: f64,
    /// Joint angle about z-axis (radians) — variable for revolute joints.
    pub theta: f64,
}
impl DhParam {
    /// Create a new DH parameter set.
    pub fn new(a: f64, alpha: f64, d: f64, theta: f64) -> Self {
        Self { a, alpha, d, theta }
    }
}
/// A serial-chain robot manipulator with DH parameters and utility methods.
///
/// Wraps [`KinematicChain`] and adds named presets and convenience methods.
pub struct SerialManipulator {
    /// The underlying kinematic chain.
    pub chain: KinematicChain,
    /// Human-readable name (e.g., "PUMA 560").
    pub name: String,
}
impl SerialManipulator {
    /// Create a serial manipulator from a kinematic chain.
    pub fn new(name: impl Into<String>, chain: KinematicChain) -> Self {
        Self {
            chain,
            name: name.into(),
        }
    }
    /// Preset: a 6-DOF PUMA-like robot.
    ///
    /// Uses approximate DH parameters for a classic 6-DOF revolute chain.
    pub fn puma_like() -> Self {
        let params = vec![
            DhParam::new(0.0, PI / 2.0, 0.0, 0.0),
            DhParam::new(0.432, 0.0, -0.149, 0.0),
            DhParam::new(0.0, PI / 2.0, 0.0, 0.0),
            DhParam::new(0.0, -PI / 2.0, 0.433, 0.0),
            DhParam::new(0.0, PI / 2.0, 0.0, 0.0),
            DhParam::new(0.0, 0.0, 0.0, 0.0),
        ];
        let types = vec![JointType::Revolute; 6];
        let q0 = vec![0.0; 6];
        Self::new("PUMA-like", KinematicChain::new(params, types, q0))
    }
    /// Preset: a 3-DOF spatial arm (RRR, unit-length links).
    pub fn rrr_unit() -> Self {
        let params = vec![
            DhParam::new(1.0, 0.0, 0.0, 0.0),
            DhParam::new(1.0, 0.0, 0.0, 0.0),
            DhParam::new(0.5, 0.0, 0.0, 0.0),
        ];
        let types = vec![JointType::Revolute; 3];
        Self::new("RRR-unit", KinematicChain::new(params, types, vec![0.0; 3]))
    }
    /// Forward kinematics: returns end-effector position.
    pub fn tip_position(&self) -> [f64; 3] {
        self.chain.end_effector_position()
    }
    /// Forward kinematics: returns full 4×4 pose.
    pub fn tip_pose(&self) -> [[f64; 4]; 4] {
        self.chain.end_effector_pose()
    }
    /// Degrees of freedom.
    pub fn dof(&self) -> usize {
        self.chain.dof()
    }
    /// Get current joint configuration.
    pub fn joint_config(&self) -> &[f64] {
        &self.chain.q
    }
    /// Set joint configuration (clamped to limits).
    pub fn set_config(&mut self, q: Vec<f64>) {
        if q.len() == self.chain.dof() {
            self.chain.q = q;
            self.chain.clamp_joints();
        }
    }
    /// Compute geometric Jacobian as a [`JacobianMatrix`].
    pub fn jacobian_matrix(&self) -> JacobianMatrix {
        JacobianMatrix::from_chain(&self.chain)
    }
    /// Check if current configuration is near a singularity.
    ///
    /// Returns `true` if the manipulability measure is below `threshold`.
    pub fn is_near_singularity(&self, threshold: f64) -> bool {
        self.chain.manipulability() < threshold
    }
}
/// Link inertial parameters for dynamics computation.
#[derive(Debug, Clone)]
pub struct LinkInertia {
    /// Link mass (kg).
    pub mass: f64,
    /// Centre of mass in link frame (m).
    pub com: [f64; 3],
    /// Inertia tensor about CoM in link frame (3×3, row-major).
    pub inertia: [[f64; 3]; 3],
}
impl LinkInertia {
    /// Create a uniform-density box link with given half-extents.
    pub fn box_link(mass: f64, hx: f64, hy: f64, hz: f64) -> Self {
        let ixx = mass / 12.0 * (4.0 * hy * hy + 4.0 * hz * hz);
        let iyy = mass / 12.0 * (4.0 * hx * hx + 4.0 * hz * hz);
        let izz = mass / 12.0 * (4.0 * hx * hx + 4.0 * hy * hy);
        Self {
            mass,
            com: [0.0; 3],
            inertia: [[ixx, 0.0, 0.0], [0.0, iyy, 0.0], [0.0, 0.0, izz]],
        }
    }
    /// Create a point-mass link (zero inertia tensor).
    pub fn point_mass(mass: f64, com: [f64; 3]) -> Self {
        Self {
            mass,
            com,
            inertia: [[0.0; 3]; 3],
        }
    }
}
/// A single waypoint for trajectory planning.
#[derive(Debug, Clone)]
pub struct Waypoint {
    /// Joint configuration at this waypoint.
    pub q: Vec<f64>,
    /// Time at which to reach this waypoint (s).
    pub time: f64,
}
impl Waypoint {
    /// Create a new waypoint.
    pub fn new(q: Vec<f64>, time: f64) -> Self {
        Self { q, time }
    }
}
