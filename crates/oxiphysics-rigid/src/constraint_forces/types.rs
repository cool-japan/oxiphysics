//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use std::f64::consts::PI;

/// Multiply two 4×4 homogeneous transform matrices.
fn h4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// A screw axis S = \[ω; v\] representing a joint in PoE formulation.
#[derive(Debug, Clone, Copy)]
pub struct ScrewAxis {
    /// Rotation axis ω (unit vector for revolute; zero for prismatic).
    pub omega: [f64; 3],
    /// Linear velocity part v (= -ω × q for revolute; axis direction for prismatic).
    pub linear: [f64; 3],
}
impl ScrewAxis {
    /// Create a revolute screw axis given axis ω (unit) and a point q on the axis.
    pub fn revolute(omega: [f64; 3], q: [f64; 3]) -> Self {
        let v = vec3_scale(vec3_cross(omega, q), -1.0);
        Self { omega, linear: v }
    }
    /// Create a prismatic screw axis with direction ω (unit, treated as linear).
    pub fn prismatic(direction: [f64; 3]) -> Self {
        Self {
            omega: [0.0; 3],
            linear: direction,
        }
    }
    /// Matrix exponential e^{\[S\]θ} as 4×4 SE(3) matrix.
    pub fn exp_map(&self, theta: f64) -> [[f64; 4]; 4] {
        let w = self.omega;
        let v = self.linear;
        let w_norm = vec3_norm(w);
        if w_norm < 1e-12 {
            let t = vec3_scale(v, theta);
            return [
                [1.0, 0.0, 0.0, t[0]],
                [0.0, 1.0, 0.0, t[1]],
                [0.0, 0.0, 1.0, t[2]],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }
        let r = rodrigues(w, theta);
        let k = skew(w);
        let k2 = mat3_mul(k, k);
        let g_mat = mat3_add(
            mat3_add(
                mat3_scale(mat3_identity(), theta),
                mat3_scale(k, 1.0 - theta.cos()),
            ),
            mat3_scale(k2, theta - theta.sin()),
        );
        let gv = mat3_vec(g_mat, v);
        [
            [r[0][0], r[0][1], r[0][2], gv[0]],
            [r[1][0], r[1][1], r[1][2], gv[1]],
            [r[2][0], r[2][1], r[2][2], gv[2]],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}
/// Pin wear model for revolute joint bearing (Archard's wear law).
pub struct JointWear {
    /// Archard wear coefficient K (Pa^-1 or dimensionless).
    pub wear_coefficient: f64,
    /// Contact pressure p (Pa).
    pub contact_pressure: f64,
    /// Sliding distance per cycle (m).
    pub sliding_per_cycle: f64,
    /// Hardness of softer material H (Pa).
    pub hardness: f64,
    /// Accumulated wear depth (m).
    pub wear_depth: f64,
    /// Accumulated cycles.
    pub cycles: u64,
}
impl JointWear {
    /// Create a new joint wear model.
    pub fn new(
        wear_coefficient: f64,
        contact_pressure: f64,
        sliding_per_cycle: f64,
        hardness: f64,
    ) -> Self {
        JointWear {
            wear_coefficient,
            contact_pressure,
            sliding_per_cycle,
            hardness,
            wear_depth: 0.0,
            cycles: 0,
        }
    }
    /// Archard wear volume per cycle: ΔV = K * P * s / H.
    pub fn wear_volume_per_cycle(&self, contact_area: f64) -> f64 {
        self.wear_coefficient * self.contact_pressure * contact_area * self.sliding_per_cycle
            / self.hardness
    }
    /// Advance by `n` cycles and accumulate wear depth.
    pub fn advance(&mut self, n: u64, contact_area: f64) {
        let dv = self.wear_volume_per_cycle(contact_area);
        self.wear_depth += dv * n as f64 / contact_area;
        self.cycles += n;
    }
    /// Estimated remaining life cycles before wear depth exceeds `limit_m`.
    pub fn remaining_life(&self, limit_m: f64, contact_area: f64) -> u64 {
        let dv = self.wear_volume_per_cycle(contact_area);
        if dv < 1e-30 {
            return u64::MAX;
        }
        let remaining = (limit_m - self.wear_depth) * contact_area / dv;
        remaining.max(0.0) as u64
    }
}
/// Spatial inertia (6×6) stored in compact form for a rigid body.
///
/// G = \[I_body   0; 0   m*I₃\] where I_body is the 3×3 inertia tensor.
#[derive(Debug, Clone)]
pub struct SpatialInertia {
    /// Body-frame inertia tensor rows: \[Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz, Izz\].
    pub inertia: [[f64; 3]; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Center of mass in body frame \[m\].
    pub com: [f64; 3],
}
impl SpatialInertia {
    /// Create a spatial inertia for a uniform rectangular box.
    pub fn box_inertia(mass: f64, lx: f64, ly: f64, lz: f64) -> Self {
        let ixx = mass / 12.0 * (ly * ly + lz * lz);
        let iyy = mass / 12.0 * (lx * lx + lz * lz);
        let izz = mass / 12.0 * (lx * lx + ly * ly);
        SpatialInertia {
            inertia: [[ixx, 0.0, 0.0], [0.0, iyy, 0.0], [0.0, 0.0, izz]],
            mass,
            com: [0.0; 3],
        }
    }
    /// Create a spatial inertia for a solid sphere.
    pub fn sphere_inertia(mass: f64, radius: f64) -> Self {
        let i = 2.0 / 5.0 * mass * radius * radius;
        SpatialInertia {
            inertia: [[i, 0.0, 0.0], [0.0, i, 0.0], [0.0, 0.0, i]],
            mass,
            com: [0.0; 3],
        }
    }
    /// Kinetic energy T = 0.5 * V^T * G * V.
    pub fn kinetic_energy(&self, twist: &Twist) -> f64 {
        let w = twist.omega;
        let v = twist.linear;
        let i = &self.inertia;
        let iw = [
            i[0][0] * w[0] + i[0][1] * w[1] + i[0][2] * w[2],
            i[1][0] * w[0] + i[1][1] * w[1] + i[1][2] * w[2],
            i[2][0] * w[0] + i[2][1] * w[1] + i[2][2] * w[2],
        ];
        0.5 * (vec3_dot(w, iw) + self.mass * vec3_dot(v, v))
    }
    /// Compute spatial momentum (wrench dual) G*V.
    pub fn momentum(&self, twist: &Twist) -> Wrench {
        let w = twist.omega;
        let v = twist.linear;
        let i = &self.inertia;
        let iw = [
            i[0][0] * w[0] + i[0][1] * w[1] + i[0][2] * w[2],
            i[1][0] * w[0] + i[1][1] * w[1] + i[1][2] * w[2],
            i[2][0] * w[0] + i[2][1] * w[1] + i[2][2] * w[2],
        ];
        Wrench::new(iw, vec3_scale(v, self.mass))
    }
}
/// Records joint reaction forces over time and computes FFT-based frequency content.
pub struct ReactionForceLogger {
    /// Time samples (s).
    pub times: Vec<f64>,
    /// Force samples (N) at each time step.
    pub forces: Vec<[f64; 3]>,
    /// Maximum stored samples before wrapping.
    pub max_samples: usize,
}
impl ReactionForceLogger {
    /// Create a new logger with a given capacity.
    pub fn new(max_samples: usize) -> Self {
        ReactionForceLogger {
            times: Vec::with_capacity(max_samples),
            forces: Vec::with_capacity(max_samples),
            max_samples,
        }
    }
    /// Record a new sample.
    pub fn record(&mut self, t: f64, force: [f64; 3]) {
        if self.times.len() >= self.max_samples {
            self.times.remove(0);
            self.forces.remove(0);
        }
        self.times.push(t);
        self.forces.push(force);
    }
    /// Mean force vector over all recorded samples.
    pub fn mean_force(&self) -> [f64; 3] {
        if self.forces.is_empty() {
            return [0.0; 3];
        }
        let n = self.forces.len() as f64;
        let mut sum = [0.0f64; 3];
        for f in &self.forces {
            sum[0] += f[0];
            sum[1] += f[1];
            sum[2] += f[2];
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
    /// Peak force magnitude over all recorded samples.
    pub fn peak_force_magnitude(&self) -> f64 {
        self.forces
            .iter()
            .map(|f| vec3_norm(*f))
            .fold(0.0_f64, f64::max)
    }
    /// DFT magnitude spectrum of Fx component. Returns `(frequency_hz, magnitude)` pairs.
    pub fn dft_fx(&self, sample_rate: f64) -> Vec<(f64, f64)> {
        let n = self.forces.len();
        if n == 0 {
            return Vec::new();
        }
        let fx: Vec<f64> = self.forces.iter().map(|f| f[0]).collect();
        let mut result = Vec::with_capacity(n / 2);
        for k in 0..n / 2 {
            let mut re = 0.0f64;
            let mut im = 0.0f64;
            for (j, &x) in fx.iter().enumerate() {
                let angle = -2.0 * PI * k as f64 * j as f64 / n as f64;
                re += x * angle.cos();
                im += x * angle.sin();
            }
            let mag = (re * re + im * im).sqrt() / n as f64;
            let freq = k as f64 * sample_rate / n as f64;
            result.push((freq, mag));
        }
        result
    }
}
/// Null-space projection for redundant manipulators (n > m).
///
/// Δθ = J⁺ * Δx + (I - J⁺J) * z, where z is a secondary gradient task.
pub struct NullSpaceProjector {
    /// Number of joints n.
    pub n_joints: usize,
    /// Task space dimension m (typically 3).
    pub m_task: usize,
    /// Damping for pseudoinverse.
    pub damping: f64,
}
impl NullSpaceProjector {
    /// Create a null-space projector.
    pub fn new(n_joints: usize, m_task: usize, damping: f64) -> Self {
        Self {
            n_joints,
            m_task,
            damping,
        }
    }
    /// Compute (I - J⁺J) * z.
    ///
    /// `j` is m×n row-major, `z` is the secondary task gradient of length n.
    pub fn null_space_component(&self, j: &[f64], z: &[f64]) -> Vec<f64> {
        let n = self.n_joints;
        let m = self.m_task;
        let lam2 = self.damping * self.damping;
        let mut jjt = vec![0.0f64; m * m];
        for r in 0..m {
            for c in 0..m {
                for k in 0..n {
                    jjt[r * m + c] += j[r * n + k] * j[c * n + k];
                }
                if r == c {
                    jjt[r * m + c] += lam2;
                }
            }
        }
        let mut jz = vec![0.0f64; m];
        for r in 0..m {
            for k in 0..n {
                jz[r] += j[r * n + k] * z[k];
            }
        }
        let w = if m == 3 {
            let a3 = [
                [jjt[0], jjt[1], jjt[2]],
                [jjt[3], jjt[4], jjt[5]],
                [jjt[6], jjt[7], jjt[8]],
            ];
            let sol = solve3x3(a3, [jz[0], jz[1], jz[2]]);
            vec![sol[0], sol[1], sol[2]]
        } else {
            jz.clone()
        };
        let mut jpjz = vec![0.0f64; n];
        for k in 0..n {
            for r in 0..m {
                jpjz[k] += j[r * n + k] * w[r];
            }
        }
        let mut result = vec![0.0f64; n];
        for k in 0..n {
            result[k] = z[k] - jpjz[k];
        }
        result
    }
    /// Joint-limit avoidance gradient ∂H/∂q for joint i.
    /// H = Σ ((q_i - q_mid_i) / (q_max_i - q_min_i))².
    pub fn joint_limit_gradient(q: &[f64], q_min: &[f64], q_max: &[f64]) -> Vec<f64> {
        let n = q.len();
        let mut grad = vec![0.0f64; n];
        for i in 0..n {
            let range = q_max[i] - q_min[i];
            if range.abs() < 1e-12 {
                continue;
            }
            let mid = (q_max[i] + q_min[i]) / 2.0;
            grad[i] = (q[i] - mid) / (range * range);
        }
        grad
    }
}
/// Cable-driven parallel robot with m cables and 3-DOF (translational) platform.
#[derive(Debug, Clone)]
pub struct CableDrivenRobot {
    /// Set of cables.
    pub cables: Vec<Cable>,
    /// Platform mass \[kg\].
    pub mass: f64,
    /// Current platform position \[m\].
    pub position: [f64; 3],
}
impl CableDrivenRobot {
    /// Create a cable-driven robot.
    pub fn new(cables: Vec<Cable>, mass: f64) -> Self {
        Self {
            cables,
            mass,
            position: [0.0; 3],
        }
    }
    /// Net cable force on the platform \[N\].
    pub fn net_cable_force(&self) -> [f64; 3] {
        let mut f = [0.0f64; 3];
        for c in &self.cables {
            let fc = c.force_on_platform(self.position);
            f = vec3_add(f, fc);
        }
        f
    }
    /// Structure matrix A (3 × m) where A * τ = f.
    /// Each column is the unit direction from attachment to anchor.
    pub fn structure_matrix(&self) -> Vec<f64> {
        let m = self.cables.len();
        let mut a = vec![0.0f64; 3 * m];
        for (j, c) in self.cables.iter().enumerate() {
            let dir = c.direction(self.position);
            a[j] = -dir[0];
            a[m + j] = -dir[1];
            a[2 * m + j] = -dir[2];
        }
        a
    }
    /// Minimum-norm cable tensions for a desired force wrench \[N\].
    /// Solves: A τ = f_desired with τ ≥ 0 (simplified: use pseudoinverse + offset).
    pub fn min_norm_tensions(&self, f_desired: [f64; 3]) -> Vec<f64> {
        let m = self.cables.len();
        let a = self.structure_matrix();
        let mut aat = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..m {
                    aat[r][c] += a[r * m + k] * a[c * m + k];
                }
            }
        }
        let y = solve3x3(aat, f_desired);
        let mut tau = vec![0.0f64; m];
        for k in 0..m {
            for r in 0..3 {
                tau[k] += a[r * m + k] * y[r];
            }
            tau[k] = tau[k].max(0.0);
        }
        tau
    }
}
/// Kinematic chain defined by a sequence of DH parameters.
#[derive(Debug, Clone)]
pub struct KinematicChain {
    /// Joint DH parameters (one per joint).
    pub joints: Vec<DHParams>,
    /// Number of links.
    pub n_links: usize,
    /// Current joint angles / positions.
    pub q: Vec<f64>,
    /// Link lengths (for planar chains).
    pub link_lengths: Vec<f64>,
    /// Joint lower limits.
    pub q_min: Vec<f64>,
    /// Joint upper limits.
    pub q_max: Vec<f64>,
}
impl KinematicChain {
    /// Create a kinematic chain from DH joint parameters.
    pub fn new(joints: Vec<DHParams>) -> Self {
        let n = joints.len();
        Self {
            n_links: n,
            q: vec![0.0; n],
            link_lengths: joints.iter().map(|j| j.a).collect(),
            q_min: vec![-std::f64::consts::PI; n],
            q_max: vec![std::f64::consts::PI; n],
            joints,
        }
    }
    /// Create a planar kinematic chain with `n` revolute joints, each link of
    /// length `link_length`, all rotating about the Z axis.
    pub fn planar(n: usize, link_length: f64) -> Self {
        let joints: Vec<DHParams> = (0..n)
            .map(|_| DHParams {
                theta: 0.0,
                d: 0.0,
                a: link_length,
                alpha: 0.0,
                is_revolute: true,
            })
            .collect();
        Self {
            n_links: n,
            q: vec![0.0; n],
            link_lengths: vec![link_length; n],
            q_min: vec![-std::f64::consts::PI; n],
            q_max: vec![std::f64::consts::PI; n],
            joints,
        }
    }
    /// Compute the 4×4 DH transform for joint `i` at its current angle.
    pub fn dh_transform(&self, i: usize) -> [[f64; 4]; 4] {
        let dh = &self.joints[i];
        let theta = if dh.is_revolute {
            dh.theta + self.q[i]
        } else {
            dh.theta
        };
        let d = if !dh.is_revolute {
            dh.d + self.q[i]
        } else {
            dh.d
        };
        let ct = theta.cos();
        let st = theta.sin();
        let ca = dh.alpha.cos();
        let sa = dh.alpha.sin();
        [
            [ct, -st * ca, st * sa, dh.a * ct],
            [st, ct * ca, -ct * sa, dh.a * st],
            [0.0, sa, ca, d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
    /// Compute the full 4×4 end-effector pose by chaining all joint transforms.
    pub fn end_effector_pose(&self) -> [[f64; 4]; 4] {
        let mut result = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for i in 0..self.n_links {
            let t = self.dh_transform(i);
            result = h4_mul(result, t);
        }
        result
    }
    /// End-effector position \[x, y, z\].
    pub fn end_effector_position(&self) -> [f64; 3] {
        let pose = self.end_effector_pose();
        [pose[0][3], pose[1][3], pose[2][3]]
    }
    /// Total reach (sum of link lengths).
    pub fn total_reach(&self) -> f64 {
        self.link_lengths.iter().sum()
    }
    /// Jacobian matrix as `Vec<Vec`f64`>` (3×n).
    pub fn jacobian_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_links;
        let eps = 1e-7;
        let p0 = self.end_effector_position();
        let mut jac = vec![vec![0.0; n]; 3];
        for col in 0..n {
            let mut chain_tmp = self.clone();
            chain_tmp.q[col] += eps;
            let p1 = chain_tmp.end_effector_position();
            for row in 0..3 {
                jac[row][col] = (p1[row] - p0[row]) / eps;
            }
        }
        jac
    }
    /// Flat Jacobian (3*n column-major) using stored joint angles.
    pub fn jacobian(&self) -> Vec<f64> {
        let jm = self.jacobian_matrix();
        let n = self.n_links;
        let mut flat = vec![0.0; 3 * n];
        for row in 0..3 {
            for col in 0..n {
                flat[row * n + col] = jm[row][col];
            }
        }
        flat
    }
    /// Manipulability index (sqrt of det(J J^T)).
    pub fn manipulability(&self) -> f64 {
        let jm = self.jacobian_matrix();
        let n = self.n_links;
        // J is 3×n, compute J·J^T (3×3)
        let mut jjt = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..n {
                    jjt[r][c] += jm[r][k] * jm[c][k];
                }
            }
        }
        let det = jjt[0][0] * (jjt[1][1] * jjt[2][2] - jjt[1][2] * jjt[2][1])
            - jjt[0][1] * (jjt[1][0] * jjt[2][2] - jjt[1][2] * jjt[2][0])
            + jjt[0][2] * (jjt[1][0] * jjt[2][1] - jjt[1][1] * jjt[2][0]);
        if det > 0.0 { det.sqrt() } else { 0.0 }
    }
    /// Check whether all joints are within their limits.
    pub fn within_limits(&self) -> bool {
        for i in 0..self.n_links {
            if self.q[i] < self.q_min[i] || self.q[i] > self.q_max[i] {
                return false;
            }
        }
        true
    }
    /// Clamp all joint values to their limits.
    pub fn clamp_joints(&mut self) {
        for i in 0..self.n_links {
            self.q[i] = self.q[i].clamp(self.q_min[i], self.q_max[i]);
        }
    }
    /// Forward kinematics: compute end-effector position given joint angles.
    /// Returns position \[px, py, pz\] in the base frame.
    pub fn forward_kinematics(&self, joint_angles: &[f64]) -> [f64; 3] {
        let mut pos = [0.0f64; 3];
        let mut rot = mat3_identity();
        for (i, dh) in self.joints.iter().enumerate() {
            let theta = if dh.is_revolute {
                dh.theta + joint_angles.get(i).copied().unwrap_or(0.0)
            } else {
                dh.theta
            };
            let d = if !dh.is_revolute {
                dh.d + joint_angles.get(i).copied().unwrap_or(0.0)
            } else {
                dh.d
            };
            let p_local = DHParams { theta, d, ..*dh };
            let t_local = p_local.translation();
            let r_local = p_local.rotation();
            let dp = mat3_vec(rot, t_local);
            pos = vec3_add(pos, dp);
            rot = mat3_mul(rot, r_local);
        }
        pos
    }
    /// Number of joints (degrees of freedom).
    pub fn dof(&self) -> usize {
        self.joints.len()
    }
    /// Compute the geometric Jacobian J (3×n) for the chain.
    /// Returns a flat Vec`f64` of size 3*n in column-major order.
    pub fn geometric_jacobian(&self, joint_angles: &[f64]) -> Vec<f64> {
        let n = self.dof();
        let eps = 1e-7;
        let p0 = self.forward_kinematics(joint_angles);
        let mut j = vec![0.0f64; 3 * n];
        let mut qa = joint_angles.to_vec();
        for i in 0..n {
            qa[i] += eps;
            let p1 = self.forward_kinematics(&qa);
            qa[i] -= eps;
            j[i] = (p1[0] - p0[0]) / eps;
            j[n + i] = (p1[1] - p0[1]) / eps;
            j[2 * n + i] = (p1[2] - p0[2]) / eps;
        }
        j
    }
}
/// A single cable (winch) connecting a fixed anchor to a moving point.
#[derive(Debug, Clone)]
pub struct Cable {
    /// Anchor position in world frame \[m\].
    pub anchor: [f64; 3],
    /// Attachment point on the platform \[m\] (in platform body frame).
    pub attachment: [f64; 3],
    /// Current cable length \[m\].
    pub length: f64,
    /// Minimum cable length (retracted) \[m\].
    pub min_length: f64,
    /// Maximum cable length \[m\].
    pub max_length: f64,
    /// Cable stiffness \[N/m\].
    pub stiffness: f64,
}
impl Cable {
    /// Create a cable.
    pub fn new(
        anchor: [f64; 3],
        attachment: [f64; 3],
        length: f64,
        min_length: f64,
        max_length: f64,
        stiffness: f64,
    ) -> Self {
        Self {
            anchor,
            attachment,
            length,
            min_length,
            max_length,
            stiffness,
        }
    }
    /// Compute the cable direction unit vector (anchor → platform attachment).
    pub fn direction(&self, platform_pos: [f64; 3]) -> [f64; 3] {
        let att_world = vec3_add(platform_pos, self.attachment);
        let d = vec3_sub(att_world, self.anchor);
        let n = vec3_norm(d);
        if n < 1e-12 {
            return [0.0; 3];
        }
        vec3_scale(d, 1.0 / n)
    }
    /// Current geometric length between anchor and platform attachment.
    pub fn geometric_length(&self, platform_pos: [f64; 3]) -> f64 {
        let att_world = vec3_add(platform_pos, self.attachment);
        vec3_norm(vec3_sub(att_world, self.anchor))
    }
    /// Cable tension based on length vs. current driven length \[N\].
    /// Positive tension only (cables cannot push).
    pub fn tension(&self, platform_pos: [f64; 3]) -> f64 {
        let l_geo = self.geometric_length(platform_pos);
        let ext = l_geo - self.length;
        (self.stiffness * ext).max(0.0)
    }
    /// Force vector on the platform \[N\].
    pub fn force_on_platform(&self, platform_pos: [f64; 3]) -> [f64; 3] {
        let t = self.tension(platform_pos);
        let dir = self.direction(platform_pos);
        vec3_scale(dir, -t)
    }
}
/// Static equilibrium analysis of reaction forces at joints.
pub struct JointLoadAnalysis {
    /// Positions of joints \[x, y, z\] (m).
    pub joint_positions: Vec<[f64; 3]>,
    /// Applied external forces at each node \[Fx, Fy, Fz\] (N).
    pub external_forces: Vec<[f64; 3]>,
    /// Applied external moments at each node \[Mx, My, Mz\] (N·m).
    pub external_moments: Vec<[f64; 3]>,
    /// Computed reaction forces (N) after analysis.
    pub reactions: Vec<[f64; 3]>,
}
impl JointLoadAnalysis {
    /// Create a new joint load analysis.
    pub fn new(
        joint_positions: Vec<[f64; 3]>,
        external_forces: Vec<[f64; 3]>,
        external_moments: Vec<[f64; 3]>,
    ) -> Self {
        let n = joint_positions.len();
        JointLoadAnalysis {
            joint_positions,
            external_forces,
            external_moments,
            reactions: vec![[0.0; 3]; n],
        }
    }
    /// Sum of all external forces \[Fx, Fy, Fz\].
    pub fn sum_forces(&self) -> [f64; 3] {
        let mut sum = [0.0f64; 3];
        for f in &self.external_forces {
            sum[0] += f[0];
            sum[1] += f[1];
            sum[2] += f[2];
        }
        sum
    }
    /// Sum of moments about the origin due to all forces and direct moments.
    pub fn sum_moments_about_origin(&self) -> [f64; 3] {
        let mut sum = [0.0f64; 3];
        for (i, f) in self.external_forces.iter().enumerate() {
            let r = self.joint_positions[i];
            let m = vec3_cross(r, *f);
            sum[0] += m[0] + self.external_moments[i][0];
            sum[1] += m[1] + self.external_moments[i][1];
            sum[2] += m[2] + self.external_moments[i][2];
        }
        sum
    }
    /// Check if the system is in static equilibrium (|ΣF| < tol and |ΣM| < tol).
    pub fn is_in_equilibrium(&self, tol: f64) -> bool {
        let sf = self.sum_forces();
        let sm = self.sum_moments_about_origin();
        vec3_norm(sf) < tol && vec3_norm(sm) < tol
    }
    /// Simple 2-joint reaction solver for a simply supported beam.
    /// Applied load `p` (N) at position `x` (m), span `l` (m).
    pub fn simply_supported_reactions(p: f64, x: f64, l: f64) -> ([f64; 3], [f64; 3]) {
        let r_b = p * x / l;
        let r_a = p - r_b;
        ([0.0, r_a, 0.0], [0.0, r_b, 0.0])
    }
}
/// Result of an inverse kinematics solve.
#[derive(Debug, Clone)]
pub struct IKResult {
    /// Solved joint angles \[rad or m\].
    pub joint_angles: Vec<f64>,
    /// Final position error magnitude.
    pub error: f64,
    /// Number of iterations used.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Weighted pseudoinverse for redundancy resolution.
///
/// Uses Δθ = W⁻¹Jᵀ(JW⁻¹Jᵀ)⁻¹ Δx.
pub struct WeightedPseudoinverse {
    /// Diagonal weight matrix W (n×n), stored as vector.
    pub weights: Vec<f64>,
}
impl WeightedPseudoinverse {
    /// Create a weighted pseudoinverse with given diagonal weights.
    pub fn new(weights: Vec<f64>) -> Self {
        Self { weights }
    }
    /// Compute weighted pseudoinverse step Δθ = W⁻¹ Jᵀ (J W⁻¹ Jᵀ)⁻¹ Δx.
    /// `j` is 3×n row-major. Returns Δθ of length n.
    pub fn step(&self, j: &[f64], dx: [f64; 3]) -> Vec<f64> {
        let n = self.weights.len();
        let mut winv_jt = vec![0.0f64; n * 3];
        for i in 0..n {
            let w_inv = if self.weights[i].abs() > 1e-14 {
                1.0 / self.weights[i]
            } else {
                0.0
            };
            for r in 0..3 {
                winv_jt[i * 3 + r] = w_inv * j[r * n + i];
            }
        }
        let mut jwjt = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..n {
                    jwjt[r][c] += j[r * n + k] * winv_jt[k * 3 + c];
                }
            }
        }
        let y = solve3x3(jwjt, dx);
        let mut dq = vec![0.0f64; n];
        for i in 0..n {
            for r in 0..3 {
                dq[i] += winv_jt[i * 3 + r] * y[r];
            }
        }
        dq
    }
}
/// Stewart platform (hexapod) geometry for 6-DOF parallel manipulator.
#[derive(Debug, Clone)]
pub struct StewartPlatform {
    /// Base attachment points in world frame (6 joints).
    pub base_points: Vec<[f64; 3]>,
    /// Platform attachment points in platform body frame (6 joints).
    pub platform_points: Vec<[f64; 3]>,
    /// Minimum leg length \[m\].
    pub leg_min: f64,
    /// Maximum leg length \[m\].
    pub leg_max: f64,
}
impl StewartPlatform {
    /// Create a Stewart platform with regular hexagonal geometry.
    ///
    /// `r_base`: base circle radius \[m\], `r_platform`: platform circle radius \[m\].
    pub fn regular(r_base: f64, r_platform: f64, leg_min: f64, leg_max: f64) -> Self {
        let n = 6;
        let mut base_points = Vec::with_capacity(n);
        let mut platform_points = Vec::with_capacity(n);
        for i in 0..n {
            let angle_b = 2.0 * PI * i as f64 / n as f64;
            let angle_p = 2.0 * PI * i as f64 / n as f64 + PI / n as f64;
            base_points.push([r_base * angle_b.cos(), r_base * angle_b.sin(), 0.0]);
            platform_points.push([r_platform * angle_p.cos(), r_platform * angle_p.sin(), 0.0]);
        }
        Self {
            base_points,
            platform_points,
            leg_min,
            leg_max,
        }
    }
    /// Inverse kinematics: given platform pose (position + z-rotation), compute leg lengths.
    pub fn inverse_kinematics(&self, position: [f64; 3], yaw: f64) -> Vec<f64> {
        let cos_y = yaw.cos();
        let sin_y = yaw.sin();
        self.platform_points
            .iter()
            .zip(self.base_points.iter())
            .map(|(pp, bp)| {
                let p_world = [
                    position[0] + cos_y * pp[0] - sin_y * pp[1],
                    position[1] + sin_y * pp[0] + cos_y * pp[1],
                    position[2] + pp[2],
                ];
                let leg = vec3_sub(p_world, *bp);
                vec3_norm(leg)
            })
            .collect()
    }
    /// Check if leg lengths are within joint limits.
    pub fn check_limits(&self, leg_lengths: &[f64]) -> bool {
        leg_lengths
            .iter()
            .all(|&l| l >= self.leg_min && l <= self.leg_max)
    }
    /// Forward kinematics (numerical, Newton-Raphson starting from initial guess).
    /// Returns platform position for given leg lengths.
    pub fn forward_kinematics_numerical(
        &self,
        target_lengths: &[f64],
        initial_pos: [f64; 3],
        initial_yaw: f64,
    ) -> ([f64; 3], f64, bool) {
        let mut pos = initial_pos;
        let yaw = initial_yaw;
        let tol = 1e-8;
        let max_iter = 100;
        for _iter in 0..max_iter {
            let current_lengths = self.inverse_kinematics(pos, yaw);
            let err: f64 = current_lengths
                .iter()
                .zip(target_lengths.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            if err < tol {
                return (pos, yaw, true);
            }
            let dz = 1e-6;
            let pos_dz = [pos[0], pos[1], pos[2] + dz];
            let l_dz = self.inverse_kinematics(pos_dz, yaw);
            let dzde: f64 = current_lengths
                .iter()
                .zip(target_lengths.iter())
                .zip(l_dz.iter())
                .map(|((a, b), ad)| (a - b) * (ad - a) / dz)
                .sum();
            if dzde.abs() > 1e-14 {
                pos[2] -= err * err.signum() * 0.5 / dzde.abs().max(1.0) * dz;
            }
        }
        (pos, yaw, false)
    }
    /// Workspace check: test if position is reachable.
    pub fn is_reachable(&self, position: [f64; 3], yaw: f64) -> bool {
        let lengths = self.inverse_kinematics(position, yaw);
        self.check_limits(&lengths)
    }
}
/// Impedance controller: M_d ẍ + B_d ẋ + K_d (x - x_d) = F_ext.
#[derive(Debug, Clone)]
pub struct ImpedanceController {
    /// Desired virtual mass M_d (diagonal, 6-vector) \[kg\].
    pub mass_desired: Vec<f64>,
    /// Desired virtual damping B_d (diagonal, 6-vector) \[N·s/m\].
    pub damping_desired: Vec<f64>,
    /// Desired virtual stiffness K_d (diagonal, 6-vector) \[N/m\].
    pub stiffness_desired: Vec<f64>,
    /// Desired equilibrium position x_d (6-vector).
    pub x_desired: Vec<f64>,
    /// Current velocity ẋ (6-vector).
    pub velocity: Vec<f64>,
}
impl ImpedanceController {
    /// Create an impedance controller.
    pub fn new(
        mass_desired: Vec<f64>,
        damping_desired: Vec<f64>,
        stiffness_desired: Vec<f64>,
        x_desired: Vec<f64>,
    ) -> Self {
        let n = x_desired.len();
        Self {
            mass_desired,
            damping_desired,
            stiffness_desired,
            x_desired,
            velocity: vec![0.0; n],
        }
    }
    /// Compute desired acceleration ẍ from current state and external force.
    pub fn desired_acceleration(&self, x_current: &[f64], vel: &[f64], f_ext: &[f64]) -> Vec<f64> {
        let n = self.mass_desired.len();
        let mut accel = vec![0.0f64; n];
        for i in 0..n {
            let pos_err = x_current.get(i).copied().unwrap_or(0.0) - self.x_desired[i];
            let v = vel.get(i).copied().unwrap_or(0.0);
            let f = f_ext.get(i).copied().unwrap_or(0.0);
            let m = self.mass_desired[i];
            if m.abs() > 1e-14 {
                accel[i] =
                    (f - self.damping_desired[i] * v - self.stiffness_desired[i] * pos_err) / m;
            }
        }
        accel
    }
    /// Energy stored in virtual spring: 0.5 * Σ Kᵢ (xᵢ - xdᵢ)².
    pub fn potential_energy(&self, x_current: &[f64]) -> f64 {
        self.stiffness_desired
            .iter()
            .enumerate()
            .fold(0.0, |acc, (i, &k)| {
                let e = x_current.get(i).copied().unwrap_or(0.0) - self.x_desired[i];
                acc + 0.5 * k * e * e
            })
    }
}
/// Work-energy theorem verification for a constrained rigid body system.
pub struct WorkEnergyPrinciple {
    /// Initial kinetic energy (J).
    pub ke_initial: f64,
    /// Final kinetic energy (J).
    pub ke_final: f64,
    /// Work done by external forces (J).
    pub work_external: f64,
    /// Work done by constraint forces (J) — should be zero for ideal constraints.
    pub work_constraints: f64,
    /// Energy dissipated by friction/damping (J).
    pub energy_dissipated: f64,
}
impl WorkEnergyPrinciple {
    /// Create a new work-energy principle verifier.
    pub fn new(ke_initial: f64, ke_final: f64, work_external: f64) -> Self {
        WorkEnergyPrinciple {
            ke_initial,
            ke_final,
            work_external,
            work_constraints: 0.0,
            energy_dissipated: 0.0,
        }
    }
    /// Energy balance residual: ΔKE - W_ext + W_dissipated (should be ≈ 0).
    pub fn residual(&self) -> f64 {
        (self.ke_final - self.ke_initial) - self.work_external
            + self.energy_dissipated
            + self.work_constraints
    }
    /// Check if the work-energy theorem is satisfied within tolerance.
    pub fn is_satisfied(&self, tol: f64) -> bool {
        self.residual().abs() < tol
    }
    /// Change in kinetic energy ΔKE.
    pub fn delta_ke(&self) -> f64 {
        self.ke_final - self.ke_initial
    }
}
/// A holonomic constraint C(q) = 0 defined by a Jacobian row.
#[derive(Debug, Clone)]
pub struct HolonomicConstraint {
    /// Constraint name.
    pub name: String,
    /// Constraint Jacobian row (one constraint equation, n-vector).
    pub jacobian_row: Vec<f64>,
    /// Constraint bias term (b such that J*q̈ + b = 0).
    pub bias: f64,
}
impl HolonomicConstraint {
    /// Create a holonomic constraint.
    pub fn new(name: &str, jacobian_row: Vec<f64>, bias: f64) -> Self {
        Self {
            name: name.to_string(),
            jacobian_row,
            bias,
        }
    }
    /// Constraint violation c = J*q - 0 (should be near 0).
    pub fn violation(&self, q: &[f64]) -> f64 {
        let mut v = self.bias;
        for (i, &ji) in self.jacobian_row.iter().enumerate() {
            v += ji * q.get(i).copied().unwrap_or(0.0);
        }
        v
    }
}
/// Impulsive force at a joint during collision.
pub struct ConstraintImpact {
    /// Mass of body 1 (kg).
    pub m1: f64,
    /// Mass of body 2 (kg).
    pub m2: f64,
    /// Coefficient of restitution e (0=perfectly plastic, 1=perfectly elastic).
    pub restitution: f64,
    /// Relative velocity before impact (m/s).
    pub relative_velocity: f64,
    /// Contact duration Δt (s) for peak force estimation.
    pub contact_duration: f64,
}
impl ConstraintImpact {
    /// Create a new collision impact analysis.
    pub fn new(
        m1: f64,
        m2: f64,
        restitution: f64,
        relative_velocity: f64,
        contact_duration: f64,
    ) -> Self {
        ConstraintImpact {
            m1,
            m2,
            restitution,
            relative_velocity,
            contact_duration,
        }
    }
    /// Reduced mass mr = m1 * m2 / (m1 + m2).
    pub fn reduced_mass(&self) -> f64 {
        self.m1 * self.m2 / (self.m1 + self.m2)
    }
    /// Impulse J = (1 + e) * mr * v_rel.
    pub fn impulse(&self) -> f64 {
        (1.0 + self.restitution) * self.reduced_mass() * self.relative_velocity
    }
    /// Peak contact force estimate F_peak = J / Δt (N).
    pub fn peak_force(&self) -> f64 {
        self.impulse() / self.contact_duration.max(1e-15)
    }
    /// Velocity change of body 1 after impact (m/s).
    pub fn delta_v1(&self) -> f64 {
        self.impulse() / self.m1
    }
    /// Velocity change of body 2 after impact (m/s).
    pub fn delta_v2(&self) -> f64 {
        -self.impulse() / self.m2
    }
    /// Kinetic energy lost during impact (J).
    pub fn energy_loss(&self) -> f64 {
        let j = self.impulse();
        let mr = self.reduced_mass();
        j * j / (2.0 * mr) * (1.0 - self.restitution)
    }
}
/// Check static equilibrium conditions for a rigid body system.
pub struct StaticEquilibrium {
    /// Applied forces \[Fx, Fy, Fz\] (N).
    pub forces: Vec<[f64; 3]>,
    /// Application points of forces (m).
    pub force_positions: Vec<[f64; 3]>,
    /// Applied moments \[Mx, My, Mz\] (N·m).
    pub moments: Vec<[f64; 3]>,
}
impl StaticEquilibrium {
    /// Create a new static equilibrium checker.
    pub fn new() -> Self {
        StaticEquilibrium {
            forces: Vec::new(),
            force_positions: Vec::new(),
            moments: Vec::new(),
        }
    }
    /// Add a force at a given position.
    pub fn add_force(&mut self, force: [f64; 3], position: [f64; 3]) {
        self.forces.push(force);
        self.force_positions.push(position);
    }
    /// Add a moment (couple).
    pub fn add_moment(&mut self, moment: [f64; 3]) {
        self.moments.push(moment);
    }
    /// Net force (N).
    pub fn net_force(&self) -> [f64; 3] {
        let mut sum = [0.0f64; 3];
        for f in &self.forces {
            sum[0] += f[0];
            sum[1] += f[1];
            sum[2] += f[2];
        }
        sum
    }
    /// Net moment about origin (N·m).
    pub fn net_moment(&self) -> [f64; 3] {
        let mut sum = [0.0f64; 3];
        for (i, f) in self.forces.iter().enumerate() {
            let r = self.force_positions[i];
            let m = vec3_cross(r, *f);
            sum[0] += m[0];
            sum[1] += m[1];
            sum[2] += m[2];
        }
        for m in &self.moments {
            sum[0] += m[0];
            sum[1] += m[1];
            sum[2] += m[2];
        }
        sum
    }
    /// Check if in static equilibrium (all forces and moments zero).
    pub fn check(&self, tol: f64) -> bool {
        vec3_norm(self.net_force()) < tol && vec3_norm(self.net_moment()) < tol
    }
}
/// Shear force and bending moment diagram for a rigid link.
pub struct InternalForce {
    /// Member length (m).
    pub length: f64,
    /// Distributed load intensity w (N/m), uniform.
    pub distributed_load: f64,
    /// Point loads: (position_m, force_N).
    pub point_loads: Vec<(f64, f64)>,
    /// Reaction at left end (N).
    pub reaction_left: f64,
    /// Reaction at right end (N).
    pub reaction_right: f64,
}
impl InternalForce {
    /// Create an InternalForce for a simply supported beam.
    pub fn simply_supported(
        length: f64,
        distributed_load: f64,
        point_loads: Vec<(f64, f64)>,
    ) -> Self {
        let w = distributed_load;
        let l = length;
        let wl = w * l;
        let sum_pl_moment: f64 = point_loads.iter().map(|(x, p)| p * x).sum::<f64>();
        let sum_pl: f64 = point_loads.iter().map(|(_, p)| p).sum();
        let r_b = (wl * l / 2.0 + sum_pl_moment) / l;
        let r_a = wl + sum_pl - r_b;
        InternalForce {
            length,
            distributed_load,
            point_loads,
            reaction_left: r_a,
            reaction_right: r_b,
        }
    }
    /// Shear force V(x) at position x (m).
    pub fn shear_force(&self, x: f64) -> f64 {
        let mut v = self.reaction_left - self.distributed_load * x;
        for (xi, pi) in &self.point_loads {
            if x > *xi {
                v -= pi;
            }
        }
        v
    }
    /// Bending moment M(x) at position x (m).
    pub fn bending_moment(&self, x: f64) -> f64 {
        let mut m = self.reaction_left * x - self.distributed_load * x * x / 2.0;
        for (xi, pi) in &self.point_loads {
            if x > *xi {
                m -= pi * (x - xi);
            }
        }
        m
    }
    /// Maximum bending moment and its location (m).
    pub fn max_bending_moment(&self) -> (f64, f64) {
        let steps = 1000;
        let dx = self.length / steps as f64;
        let mut max_m = 0.0f64;
        let mut max_x = 0.0f64;
        for i in 0..=steps {
            let x = i as f64 * dx;
            let m = self.bending_moment(x).abs();
            if m > max_m {
                max_m = m;
                max_x = x;
            }
        }
        (max_m, max_x)
    }
}
/// A spatial force (wrench): moment m ∈ ℝ³, force f ∈ ℝ³.
///
/// F = \[m; f\] ∈ se(3)*.
#[derive(Debug, Clone, Copy)]
pub struct Wrench {
    /// Moment part m \[N·m\].
    pub moment: [f64; 3],
    /// Force part f \[N\].
    pub force: [f64; 3],
}
impl Wrench {
    /// Create a new wrench.
    pub fn new(moment: [f64; 3], force: [f64; 3]) -> Self {
        Self { moment, force }
    }
    /// Zero wrench.
    pub fn zero() -> Self {
        Self::new([0.0; 3], [0.0; 3])
    }
    /// Add two wrenches.
    pub fn add(&self, other: &Wrench) -> Wrench {
        Wrench {
            moment: vec3_add(self.moment, other.moment),
            force: vec3_add(self.force, other.force),
        }
    }
    /// Scale wrench by scalar.
    pub fn scale(&self, s: f64) -> Wrench {
        Wrench {
            moment: vec3_scale(self.moment, s),
            force: vec3_scale(self.force, s),
        }
    }
    /// Power = F · V (wrench · twist).
    pub fn power(&self, twist: &Twist) -> f64 {
        vec3_dot(self.moment, twist.omega) + vec3_dot(self.force, twist.linear)
    }
    /// Magnitude of the force part.
    pub fn force_magnitude(&self) -> f64 {
        vec3_norm(self.force)
    }
    /// Magnitude of the moment part.
    pub fn moment_magnitude(&self) -> f64 {
        vec3_norm(self.moment)
    }
}
/// Assembles the constraint Jacobian matrix for a multi-body system.
#[derive(Debug, Clone)]
pub struct ConstraintJacobianAssembler {
    /// Number of generalized coordinates.
    pub n_dof: usize,
    /// List of constraints.
    pub constraints: Vec<HolonomicConstraint>,
}
impl ConstraintJacobianAssembler {
    /// Create a constraint Jacobian assembler.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            constraints: Vec::new(),
        }
    }
    /// Add a constraint to the system.
    pub fn add_constraint(&mut self, c: HolonomicConstraint) {
        self.constraints.push(c);
    }
    /// Get the full constraint Jacobian matrix (m × n), row-major.
    pub fn jacobian(&self) -> Vec<f64> {
        let m = self.constraints.len();
        let n = self.n_dof;
        let mut j = vec![0.0f64; m * n];
        for (row, c) in self.constraints.iter().enumerate() {
            for (col, &val) in c.jacobian_row.iter().enumerate() {
                if col < n {
                    j[row * n + col] = val;
                }
            }
        }
        j
    }
    /// Compute Lagrange multipliers λ from constraint forces.
    /// Solves (JMJ^T) λ = -J*a - bias using the 3x3 special case or iterative.
    pub fn constraint_forces(&self, q: &[f64]) -> Vec<f64> {
        self.constraints.iter().map(|c| c.violation(q)).collect()
    }
    /// Total number of constraints.
    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }
}
/// Virtual work principle for computing generalized forces from applied wrenches.
///
/// Q_i = Σ_j (F_j · (∂p_j/∂q_i) + M_j · (∂ω_j/∂q_i))
pub struct VirtualWorkSolver {
    /// Number of generalized coordinates.
    pub n_dof: usize,
    /// Applied wrenches at various body points (one per body).
    pub wrenches: Vec<Wrench>,
    /// Translational Jacobian columns ∂p/∂q (3 × n per body), stored row-major.
    pub jac_trans: Vec<Vec<f64>>,
    /// Rotational Jacobian columns ∂ω/∂q (3 × n per body), stored row-major.
    pub jac_rot: Vec<Vec<f64>>,
}
impl VirtualWorkSolver {
    /// Create a virtual work solver.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            wrenches: Vec::new(),
            jac_trans: Vec::new(),
            jac_rot: Vec::new(),
        }
    }
    /// Add a body with its applied wrench and Jacobian.
    pub fn add_body(&mut self, wrench: Wrench, jac_trans: Vec<f64>, jac_rot: Vec<f64>) {
        self.wrenches.push(wrench);
        self.jac_trans.push(jac_trans);
        self.jac_rot.push(jac_rot);
    }
    /// Compute generalized forces Q (n-vector) using virtual work.
    pub fn generalized_forces(&self) -> Vec<f64> {
        let n = self.n_dof;
        let mut q_gen = vec![0.0f64; n];
        for (b, w) in self.wrenches.iter().enumerate() {
            for i in 0..n {
                if b < self.jac_trans.len() {
                    let jt = &self.jac_trans[b];
                    if jt.len() >= 3 * n {
                        let jp = [jt[i], jt[n + i], jt[2 * n + i]];
                        q_gen[i] += vec3_dot(w.force, jp);
                    }
                }
                if b < self.jac_rot.len() {
                    let jr = &self.jac_rot[b];
                    if jr.len() >= 3 * n {
                        let jw = [jr[i], jr[n + i], jr[2 * n + i]];
                        q_gen[i] += vec3_dot(w.moment, jw);
                    }
                }
            }
        }
        q_gen
    }
    /// Check virtual work consistency: δW = Q · δq should equal F · δx + M · δω.
    pub fn virtual_work(&self, delta_q: &[f64]) -> f64 {
        let q = self.generalized_forces();
        q.iter().zip(delta_q.iter()).map(|(qi, dqi)| qi * dqi).sum()
    }
}
/// Hybrid position/force controller for constrained environments.
///
/// Partitions task space into position-controlled and force-controlled directions.
#[derive(Debug, Clone)]
pub struct HybridController {
    /// Selection matrix S (diagonal, 1 = force controlled, 0 = position controlled).
    pub selection: Vec<f64>,
    /// Force controller gains \[kp_f, ki_f, kd_f\].
    pub force_gains: [f64; 3],
    /// Position controller gains \[kp_p, ki_p, kd_p\].
    pub position_gains: [f64; 3],
    /// Desired force \[N\] (6-vector).
    pub f_desired: Vec<f64>,
    /// Desired position \[m\] (6-vector).
    pub x_desired: Vec<f64>,
    /// Force integral error accumulator.
    pub force_integral: Vec<f64>,
    /// Position integral error accumulator.
    pub pos_integral: Vec<f64>,
    /// Previous force error (for derivative term).
    pub prev_force_error: Vec<f64>,
    /// Previous position error (for derivative term).
    pub prev_pos_error: Vec<f64>,
}
impl HybridController {
    /// Create a hybrid position/force controller.
    pub fn new(
        selection: Vec<f64>,
        force_gains: [f64; 3],
        position_gains: [f64; 3],
        f_desired: Vec<f64>,
        x_desired: Vec<f64>,
    ) -> Self {
        let n = selection.len();
        Self {
            selection,
            force_gains,
            position_gains,
            f_desired,
            x_desired,
            force_integral: vec![0.0; n],
            pos_integral: vec![0.0; n],
            prev_force_error: vec![0.0; n],
            prev_pos_error: vec![0.0; n],
        }
    }
    /// Compute control output given current state.
    /// Returns task-space wrench command (6-vector).
    pub fn update(&mut self, x_current: &[f64], f_current: &[f64], dt: f64) -> Vec<f64> {
        let n = self.selection.len();
        let mut output = vec![0.0f64; n];
        for i in 0..n {
            let s = self.selection[i];
            let fe = self.f_desired[i] - f_current.get(i).copied().unwrap_or(0.0);
            self.force_integral[i] += fe * dt;
            let fd = (fe - self.prev_force_error[i]) / dt.max(1e-10);
            let f_out = self.force_gains[0] * fe
                + self.force_gains[1] * self.force_integral[i]
                + self.force_gains[2] * fd;
            self.prev_force_error[i] = fe;
            let pe = self.x_desired[i] - x_current.get(i).copied().unwrap_or(0.0);
            self.pos_integral[i] += pe * dt;
            let pd = (pe - self.prev_pos_error[i]) / dt.max(1e-10);
            let p_out = self.position_gains[0] * pe
                + self.position_gains[1] * self.pos_integral[i]
                + self.position_gains[2] * pd;
            self.prev_pos_error[i] = pe;
            output[i] = s * f_out + (1.0 - s) * p_out;
        }
        output
    }
}
/// A spatial velocity (twist): ω ∈ ℝ³ (angular), v ∈ ℝ³ (linear).
///
/// Follows the modern robotics convention: V = \[ω; v\] ∈ se(3).
#[derive(Debug, Clone, Copy)]
pub struct Twist {
    /// Angular velocity part ω \[rad/s\].
    pub omega: [f64; 3],
    /// Linear velocity part v \[m/s\] (at the reference point).
    pub linear: [f64; 3],
}
impl Twist {
    /// Create a new twist from angular and linear parts.
    pub fn new(omega: [f64; 3], linear: [f64; 3]) -> Self {
        Self { omega, linear }
    }
    /// Zero twist.
    pub fn zero() -> Self {
        Self::new([0.0; 3], [0.0; 3])
    }
    /// Add two twists.
    pub fn add(&self, other: &Twist) -> Twist {
        Twist {
            omega: vec3_add(self.omega, other.omega),
            linear: vec3_add(self.linear, other.linear),
        }
    }
    /// Scale a twist by a scalar.
    pub fn scale(&self, s: f64) -> Twist {
        Twist {
            omega: vec3_scale(self.omega, s),
            linear: vec3_scale(self.linear, s),
        }
    }
    /// Lie bracket \[V1, V2\] for two body twists.
    pub fn lie_bracket(&self, other: &Twist) -> Twist {
        let w1 = self.omega;
        let v1 = self.linear;
        let w2 = other.omega;
        let v2 = other.linear;
        Twist {
            omega: vec3_cross(w1, w2),
            linear: vec3_add(vec3_cross(w1, v2), vec3_cross(v1, w2)),
        }
    }
    /// Norm of the twist (√(ω·ω + v·v)).
    pub fn norm(&self) -> f64 {
        (vec3_dot(self.omega, self.omega) + vec3_dot(self.linear, self.linear)).sqrt()
    }
}
/// Mechanical power at a joint: P = F · v + M · ω.
pub struct MechanicalPower {
    /// Current force vector (N).
    pub force: [f64; 3],
    /// Current moment vector (N·m).
    pub moment: [f64; 3],
    /// Velocity of the joint attachment point (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity of the body (rad/s).
    pub angular_velocity: [f64; 3],
}
impl MechanicalPower {
    /// Create a new mechanical power calculation.
    pub fn new(
        force: [f64; 3],
        moment: [f64; 3],
        velocity: [f64; 3],
        angular_velocity: [f64; 3],
    ) -> Self {
        MechanicalPower {
            force,
            moment,
            velocity,
            angular_velocity,
        }
    }
    /// Translational power P_trans = F · v (W).
    pub fn translational_power(&self) -> f64 {
        vec3_dot(self.force, self.velocity)
    }
    /// Rotational power P_rot = M · ω (W).
    pub fn rotational_power(&self) -> f64 {
        vec3_dot(self.moment, self.angular_velocity)
    }
    /// Total power P = P_trans + P_rot (W).
    pub fn total_power(&self) -> f64 {
        self.translational_power() + self.rotational_power()
    }
}
/// Fatigue analysis at a constraint joint (rolling contact bearing life).
pub struct ConstraintFatigue {
    /// Radial load on bearing (N).
    pub radial_load: f64,
    /// Axial load on bearing (N).
    pub axial_load: f64,
    /// Basic dynamic load rating C (N).
    pub c_rating: f64,
    /// Rotational speed (rpm).
    pub rpm: f64,
    /// Hertz contact half-width a (m) for flat contact.
    pub hertz_half_width: f64,
    /// Contact modulus E* (GPa).
    pub contact_modulus: f64,
}
impl ConstraintFatigue {
    /// Create a constraint fatigue analysis object.
    pub fn new(radial_load: f64, c_rating: f64, rpm: f64) -> Self {
        ConstraintFatigue {
            radial_load,
            axial_load: 0.0,
            c_rating,
            rpm,
            hertz_half_width: 1e-3,
            contact_modulus: 200.0,
        }
    }
    /// Equivalent dynamic bearing load P (N) using ISO 281.
    pub fn equivalent_load(&self) -> f64 {
        self.radial_load + 0.0 * self.axial_load
    }
    /// L10 bearing life in millions of revolutions (ISO 281).
    pub fn l10_millions_revolutions(&self) -> f64 {
        let p = self.equivalent_load();
        if p < 1e-10 {
            return f64::INFINITY;
        }
        (self.c_rating / p).powi(3)
    }
    /// L10 bearing life in hours.
    pub fn l10_hours(&self) -> f64 {
        if self.rpm < 1e-10 {
            return f64::INFINITY;
        }
        self.l10_millions_revolutions() * 1e6 / (60.0 * self.rpm)
    }
    /// Hertz contact stress p_max (MPa) for a flat circular contact.
    /// Uses: p_max = (3P / (2πa²)) formula.
    pub fn hertz_contact_stress(&self, contact_radius: f64) -> f64 {
        let a = contact_radius;
        3.0 * self.radial_load / (2.0 * PI * a * a) / 1e6
    }
}
/// A single Denavit-Hartenberg joint parameter set.
#[derive(Debug, Clone, Copy)]
pub struct DHParams {
    /// Link length a \[m\].
    pub a: f64,
    /// Link twist α \[rad\].
    pub alpha: f64,
    /// Link offset d \[m\].
    pub d: f64,
    /// Joint angle θ \[rad\] (variable for revolute joints).
    pub theta: f64,
    /// True for revolute joint, false for prismatic.
    pub is_revolute: bool,
}
impl DHParams {
    /// Create a revolute joint DH parameter set.
    pub fn revolute(a: f64, alpha: f64, d: f64, theta: f64) -> Self {
        Self {
            a,
            alpha,
            d,
            theta,
            is_revolute: true,
        }
    }
    /// Create a prismatic joint DH parameter set.
    pub fn prismatic(a: f64, alpha: f64, d: f64, theta: f64) -> Self {
        Self {
            a,
            alpha,
            d,
            theta,
            is_revolute: false,
        }
    }
    /// 4×4 homogeneous transformation T_{i-1}^{i}.
    /// Returns as a flat \[16\] array in row-major order.
    pub fn transform(&self) -> [[f64; 4]; 4] {
        let ct = self.theta.cos();
        let st = self.theta.sin();
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();
        [
            [ct, -st * ca, st * sa, self.a * ct],
            [st, ct * ca, -ct * sa, self.a * st],
            [0.0, sa, ca, self.d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
    /// 3×3 rotation part of the DH transform.
    pub fn rotation(&self) -> [[f64; 3]; 3] {
        let t = self.transform();
        [
            [t[0][0], t[0][1], t[0][2]],
            [t[1][0], t[1][1], t[1][2]],
            [t[2][0], t[2][1], t[2][2]],
        ]
    }
    /// Translation vector \[px, py, pz\].
    pub fn translation(&self) -> [f64; 3] {
        let t = self.transform();
        [t[0][3], t[1][3], t[2][3]]
    }
}
/// Jacobian-transpose inverse kinematics solver.
///
/// Iterates: θ ← θ + α * Jᵀ * e  until convergence.
pub struct JacobianTransposeSolver {
    /// Step size α.
    pub alpha: f64,
    /// Convergence tolerance.
    pub tol: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}
impl JacobianTransposeSolver {
    /// Create a Jacobian-transpose IK solver.
    pub fn new(alpha: f64, tol: f64, max_iter: usize) -> Self {
        Self {
            alpha,
            tol,
            max_iter,
        }
    }
    /// Solve IK using the Jacobian-transpose method.
    pub fn solve(&self, chain: &KinematicChain, target: [f64; 3], initial: Vec<f64>) -> IKResult {
        let mut q = initial;
        let n = chain.dof();
        for iter in 0..self.max_iter {
            let pos = chain.forward_kinematics(&q);
            let e = vec3_sub(target, pos);
            let err = vec3_norm(e);
            if err < self.tol {
                return IKResult {
                    joint_angles: q,
                    error: err,
                    iterations: iter,
                    converged: true,
                };
            }
            let j = chain.geometric_jacobian(&q);
            for i in 0..n {
                let jt_e = j[i] * e[0] + j[n + i] * e[1] + j[2 * n + i] * e[2];
                q[i] += self.alpha * jt_e;
            }
        }
        let pos = chain.forward_kinematics(&q);
        let e = vec3_sub(target, pos);
        let err = vec3_norm(e);
        IKResult {
            joint_angles: q,
            error: err,
            iterations: self.max_iter,
            converged: false,
        }
    }
}
/// Pseudoinverse (Moore-Penrose) IK solver: Δθ = J⁺ * Δx.
///
/// Uses the left pseudoinverse for overdetermined (n<3) or right pseudoinverse for (n≥3).
pub struct PseudoinverseSolver {
    /// Convergence tolerance.
    pub tol: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Damping constant λ for DLS (0 = undamped pseudoinverse).
    pub damping: f64,
}
impl PseudoinverseSolver {
    /// Create a pseudoinverse IK solver.
    pub fn new(tol: f64, max_iter: usize, damping: f64) -> Self {
        Self {
            tol,
            max_iter,
            damping,
        }
    }
    /// Compute J*(Jᵀ + λ²I)⁻¹ * e for n ≥ 3 (right pseudo with DLS).
    /// J is 3×n, returns Δθ of length n.
    fn dls_step(j: &[f64], e: [f64; 3], n: usize, lambda: f64) -> Vec<f64> {
        let mut jjt = [[0.0f64; 3]; 3];
        for row in 0..3 {
            for col in 0..3 {
                for k in 0..n {
                    jjt[row][col] += j[row * n + k] * j[col * n + k];
                }
            }
        }
        for i in 0..3 {
            jjt[i][i] += lambda * lambda;
        }
        let y = solve3x3(jjt, e);
        let mut dq = vec![0.0f64; n];
        for i in 0..n {
            for row in 0..3 {
                dq[i] += j[row * n + i] * y[row];
            }
        }
        dq
    }
    /// Solve IK using damped least-squares (DLS) pseudoinverse.
    pub fn solve(&self, chain: &KinematicChain, target: [f64; 3], initial: Vec<f64>) -> IKResult {
        let mut q = initial;
        let n = chain.dof();
        for iter in 0..self.max_iter {
            let pos = chain.forward_kinematics(&q);
            let e = vec3_sub(target, pos);
            let err = vec3_norm(e);
            if err < self.tol {
                return IKResult {
                    joint_angles: q,
                    error: err,
                    iterations: iter,
                    converged: true,
                };
            }
            let j = chain.geometric_jacobian(&q);
            let dq = Self::dls_step(&j, e, n, self.damping);
            for i in 0..n {
                q[i] += dq[i];
            }
        }
        let pos = chain.forward_kinematics(&q);
        let err = vec3_norm(vec3_sub(target, pos));
        IKResult {
            joint_angles: q,
            error: err,
            iterations: self.max_iter,
            converged: false,
        }
    }
}
/// Product-of-exponentials kinematic chain.
#[derive(Debug, Clone)]
pub struct PoEChain {
    /// Screw axes in the home (zero) configuration.
    pub screws: Vec<ScrewAxis>,
    /// Home configuration: end-effector position when all θᵢ = 0.
    pub home_position: [f64; 3],
}
impl PoEChain {
    /// Create a PoE chain.
    pub fn new(screws: Vec<ScrewAxis>, home_position: [f64; 3]) -> Self {
        Self {
            screws,
            home_position,
        }
    }
    /// Multiply two 4×4 SE(3) matrices.
    pub fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
        let mut c = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }
    /// Forward kinematics using product of exponentials.
    pub fn forward_kinematics(&self, thetas: &[f64]) -> [f64; 3] {
        let mut t = [
            [1.0, 0.0, 0.0, self.home_position[0]],
            [0.0, 1.0, 0.0, self.home_position[1]],
            [0.0, 0.0, 1.0, self.home_position[2]],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for (i, screw) in self.screws.iter().enumerate() {
            let theta = thetas.get(i).copied().unwrap_or(0.0);
            let ei = screw.exp_map(theta);
            t = Self::mat4_mul(ei, t);
        }
        [t[0][3], t[1][3], t[2][3]]
    }
    /// Numerical space Jacobian (3×n).
    pub fn jacobian(&self, thetas: &[f64]) -> Vec<f64> {
        let n = self.screws.len();
        let eps = 1e-7;
        let p0 = self.forward_kinematics(thetas);
        let mut j = vec![0.0f64; 3 * n];
        let mut th = thetas.to_vec();
        for i in 0..n {
            th[i] += eps;
            let p1 = self.forward_kinematics(&th);
            th[i] -= eps;
            j[i] = (p1[0] - p0[0]) / eps;
            j[n + i] = (p1[1] - p0[1]) / eps;
            j[2 * n + i] = (p1[2] - p0[2]) / eps;
        }
        j
    }
}
/// Force and moment acting at a constraint (joint) between two bodies.
#[derive(Debug, Clone)]
pub struct ConstraintForce {
    /// Normal force component (N) along the constraint normal.
    pub normal_force: f64,
    /// Tangential force magnitude (N) perpendicular to constraint normal.
    pub tangential_force: f64,
    /// Torque at the joint (N·m).
    pub torque: f64,
    /// Force vector \[Fx, Fy, Fz\] (N).
    pub force_vec: [f64; 3],
    /// Moment vector \[Mx, My, Mz\] (N·m).
    pub moment_vec: [f64; 3],
}
impl ConstraintForce {
    /// Create a new constraint force from components.
    pub fn new(force_vec: [f64; 3], moment_vec: [f64; 3]) -> Self {
        let _mag = vec3_norm(force_vec);
        let tan = (force_vec[0].powi(2) + force_vec[1].powi(2)).sqrt();
        let torque = vec3_norm(moment_vec);
        ConstraintForce {
            normal_force: force_vec[2],
            tangential_force: tan,
            torque,
            force_vec,
            moment_vec,
        }
    }
    /// Total force magnitude (N).
    pub fn magnitude(&self) -> f64 {
        vec3_norm(self.force_vec)
    }
    /// Total moment magnitude (N·m).
    pub fn moment_magnitude(&self) -> f64 {
        vec3_norm(self.moment_vec)
    }
}
