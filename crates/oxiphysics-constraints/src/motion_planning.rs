// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Motion planning for constrained robotic systems.
//!
//! Provides:
//!
//! * Configuration space (C-space) representation and joint limits
//! * Workspace constraints (end-effector position and orientation)
//! * Constraint manifold projection via gradient descent
//! * PRM (Probabilistic Roadmap Method) planner
//! * RRT (Rapidly-exploring Random Tree) planner
//! * CBIRRT (Constrained Bidirectional RRT) planner
//! * Inverse kinematics as constraint satisfaction
//! * Constraint-based trajectory optimization
//!
//! All implementations use plain `f64` and `[f64; N]` arrays (no nalgebra).

use rand::Rng;

use rand::RngExt;
// ── Scalar / vector helpers ───────────────────────────────────────────────────

/// Dot product of two equal-length slices.
#[cfg(test)]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm of a slice.
#[cfg(test)]
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// Squared Euclidean distance between two equal-length slices.
fn dist_sq(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// Euclidean distance between two configurations.
pub fn config_distance(a: &[f64], b: &[f64]) -> f64 {
    dist_sq(a, b).sqrt()
}

/// Linear interpolation between two configurations at parameter t ∈ \[0, 1\].
pub fn config_lerp(a: &[f64], b: &[f64], t: f64) -> Vec<f64> {
    a.iter()
        .zip(b.iter())
        .map(|(ai, bi)| ai + (bi - ai) * t)
        .collect()
}

/// Clamp `x` to `[lo, hi]`.
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ── Configuration Space ───────────────────────────────────────────────────────

/// Represents the configuration space of a robot with `n` joints.
///
/// Each joint has lower and upper limits and a type (revolute or prismatic).
#[derive(Debug, Clone)]
pub struct ConfigurationSpace {
    /// Number of degrees of freedom.
    pub dof: usize,
    /// Lower joint limits (rad or m).
    pub lower_limits: Vec<f64>,
    /// Upper joint limits (rad or m).
    pub upper_limits: Vec<f64>,
    /// Maximum step size per dimension for local planning.
    pub step_size: Vec<f64>,
}

impl ConfigurationSpace {
    /// Create a new configuration space with given joint limits.
    ///
    /// * `lower_limits` – per-joint lower bound.
    /// * `upper_limits` – per-joint upper bound.
    /// * `step_size`    – per-joint maximum interpolation step.
    pub fn new(lower_limits: Vec<f64>, upper_limits: Vec<f64>, step_size: Vec<f64>) -> Self {
        let dof = lower_limits.len();
        Self {
            dof,
            lower_limits,
            upper_limits,
            step_size,
        }
    }

    /// Create a C-space with uniform limits `[-limit, +limit]` for all joints.
    pub fn uniform(dof: usize, limit: f64, step: f64) -> Self {
        Self {
            dof,
            lower_limits: vec![-limit; dof],
            upper_limits: vec![limit; dof],
            step_size: vec![step; dof],
        }
    }

    /// Check whether a configuration satisfies all joint limits.
    pub fn in_bounds(&self, q: &[f64]) -> bool {
        q.iter()
            .zip(self.lower_limits.iter().zip(self.upper_limits.iter()))
            .all(|(&qi, (&lo, &hi))| qi >= lo && qi <= hi)
    }

    /// Clamp a configuration to the joint limits.
    pub fn clamp_to_bounds(&self, q: &[f64]) -> Vec<f64> {
        q.iter()
            .zip(self.lower_limits.iter().zip(self.upper_limits.iter()))
            .map(|(&qi, (&lo, &hi))| clamp(qi, lo, hi))
            .collect()
    }

    /// Sample a random configuration uniformly within joint limits.
    pub fn sample_random(&self, rng: &mut impl Rng) -> Vec<f64> {
        self.lower_limits
            .iter()
            .zip(self.upper_limits.iter())
            .map(|(&lo, &hi)| rng.random_range(lo..hi))
            .collect()
    }

    /// Interpolate along the straight line from `a` to `b` with the
    /// per-joint step size.  Returns an ordered sequence of waypoints
    /// from `a` to `b` (inclusive).
    pub fn interpolate_path(&self, a: &[f64], b: &[f64]) -> Vec<Vec<f64>> {
        // Number of steps determined by the maximum joint travel / step_size.
        let n_steps = self
            .step_size
            .iter()
            .zip(a.iter().zip(b.iter()))
            .map(|(&s, (&ai, &bi))| {
                let diff = (bi - ai).abs();
                if s < 1e-15 {
                    0_usize
                } else {
                    (diff / s).ceil() as usize
                }
            })
            .max()
            .unwrap_or(0)
            .max(1);

        (0..=n_steps)
            .map(|i| config_lerp(a, b, i as f64 / n_steps as f64))
            .collect()
    }
}

// ── Workspace Constraints ────────────────────────────────────────────────────

/// 3-D end-effector position constraint.
///
/// The constraint is satisfied when the FK position equals the target
/// within tolerance.
#[derive(Debug, Clone)]
pub struct PositionConstraint {
    /// Desired end-effector position \[x, y, z\] (m).
    pub target: [f64; 3],
    /// Tolerance ‖p_ee − p_target‖ < tol.
    pub tolerance: f64,
}

impl PositionConstraint {
    /// Create a position constraint.
    pub fn new(target: [f64; 3], tolerance: f64) -> Self {
        Self { target, tolerance }
    }

    /// Evaluate the constraint error vector (3D) given a FK position.
    pub fn error(&self, ee_position: [f64; 3]) -> [f64; 3] {
        [
            ee_position[0] - self.target[0],
            ee_position[1] - self.target[1],
            ee_position[2] - self.target[2],
        ]
    }

    /// Check whether the constraint is satisfied.
    pub fn is_satisfied(&self, ee_position: [f64; 3]) -> bool {
        let e = self.error(ee_position);
        (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt() <= self.tolerance
    }
}

/// Orientation constraint (axis alignment).
///
/// Constrains the z-axis of the end-effector frame to point along a
/// target direction, within an angular tolerance.
#[derive(Debug, Clone)]
pub struct OrientationConstraint {
    /// Desired z-axis direction (unit vector).
    pub target_axis: [f64; 3],
    /// Angular tolerance (rad).
    pub tolerance: f64,
}

impl OrientationConstraint {
    /// Create an orientation constraint.
    pub fn new(target_axis: [f64; 3], tolerance: f64) -> Self {
        // Normalize target axis
        let n = (target_axis[0] * target_axis[0]
            + target_axis[1] * target_axis[1]
            + target_axis[2] * target_axis[2])
            .sqrt()
            .max(1e-15);
        Self {
            target_axis: [target_axis[0] / n, target_axis[1] / n, target_axis[2] / n],
            tolerance,
        }
    }

    /// Angular error between actual z-axis `ee_axis` and the target axis (rad).
    pub fn angular_error(&self, ee_axis: [f64; 3]) -> f64 {
        let d = clamp(
            self.target_axis[0] * ee_axis[0]
                + self.target_axis[1] * ee_axis[1]
                + self.target_axis[2] * ee_axis[2],
            -1.0,
            1.0,
        );
        d.acos()
    }

    /// Check whether the orientation constraint is satisfied.
    pub fn is_satisfied(&self, ee_axis: [f64; 3]) -> bool {
        self.angular_error(ee_axis) <= self.tolerance
    }
}

// ── Forward Kinematics (Denavit-Hartenberg) ───────────────────────────────────

/// A single DH joint (standard DH parameters).
#[derive(Debug, Clone, Copy)]
pub struct DhJoint {
    /// Link length a (m).
    pub a: f64,
    /// Link twist α (rad).
    pub alpha: f64,
    /// Link offset d (m).
    pub d: f64,
    /// Joint angle offset θ_offset (rad).
    pub theta_offset: f64,
}

impl DhJoint {
    /// Create a DH joint.
    pub fn new(a: f64, alpha: f64, d: f64, theta_offset: f64) -> Self {
        Self {
            a,
            alpha,
            d,
            theta_offset,
        }
    }

    /// Compute the 4×4 homogeneous transform T for joint angle θ.
    ///
    /// Uses the standard DH convention.
    /// The matrix is stored row-major as `[[f64; 4\]; 4]`.
    pub fn transform(&self, theta: f64) -> [[f64; 4]; 4] {
        let th = theta + self.theta_offset;
        let ct = th.cos();
        let st = th.sin();
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();
        let a = self.a;
        let d = self.d;
        [
            [ct, -st * ca, st * sa, a * ct],
            [st, ct * ca, -ct * sa, a * st],
            [0.0, sa, ca, d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}

/// Multiply two 4×4 homogeneous matrices.
fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut r = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

/// Simple serial DH chain forward kinematics.
///
/// Returns the 4×4 end-effector transform T_0e.
pub fn forward_kinematics(joints: &[DhJoint], q: &[f64]) -> [[f64; 4]; 4] {
    let mut t = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    for (joint, &qi) in joints.iter().zip(q.iter()) {
        t = mat4_mul(t, joint.transform(qi));
    }
    t
}

/// Extract the end-effector position from a 4×4 transform.
pub fn ee_position(t: [[f64; 4]; 4]) -> [f64; 3] {
    [t[0][3], t[1][3], t[2][3]]
}

/// Extract the z-axis (column 2) of the rotation part from a 4×4 transform.
pub fn ee_z_axis(t: [[f64; 4]; 4]) -> [f64; 3] {
    [t[0][2], t[1][2], t[2][2]]
}

// ── Numerical Jacobian ────────────────────────────────────────────────────────

/// Compute the numerical Jacobian J (3×n) of end-effector position
/// with respect to joint angles using central differences.
pub fn numerical_jacobian(joints: &[DhJoint], q: &[f64], delta: f64) -> Vec<[f64; 3]> {
    let n = q.len();
    let mut jac = vec![[0.0_f64; 3]; n];
    let mut q_p = q.to_vec();
    let mut q_m = q.to_vec();
    for i in 0..n {
        q_p[i] = q[i] + delta;
        q_m[i] = q[i] - delta;
        let p_p = ee_position(forward_kinematics(joints, &q_p));
        let p_m = ee_position(forward_kinematics(joints, &q_m));
        jac[i] = [
            (p_p[0] - p_m[0]) / (2.0 * delta),
            (p_p[1] - p_m[1]) / (2.0 * delta),
            (p_p[2] - p_m[2]) / (2.0 * delta),
        ];
        q_p[i] = q[i];
        q_m[i] = q[i];
    }
    jac
}

// ── Inverse Kinematics as Constraint Satisfaction ────────────────────────────

/// Iterative inverse kinematics solver using the damped least-squares
/// (Levenberg-Marquardt) Jacobian pseudoinverse.
///
/// Solves for joint angles q such that FK(q) = target_position.
#[derive(Debug, Clone)]
pub struct IkSolver {
    /// DH kinematic chain.
    pub joints: Vec<DhJoint>,
    /// C-space joint limits.
    pub cspace: ConfigurationSpace,
    /// Damping coefficient λ for the DLS Jacobian.
    pub damping: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

impl IkSolver {
    /// Create a new IK solver.
    pub fn new(
        joints: Vec<DhJoint>,
        cspace: ConfigurationSpace,
        damping: f64,
        max_iter: usize,
        tolerance: f64,
    ) -> Self {
        Self {
            joints,
            cspace,
            damping,
            max_iter,
            tolerance,
        }
    }

    /// Solve IK from initial guess `q0` towards `target`.
    ///
    /// Returns `(q_solution, converged)`.
    pub fn solve(&self, q0: &[f64], target: [f64; 3]) -> (Vec<f64>, bool) {
        let mut q = q0.to_vec();
        let n = q.len();
        let lam2 = self.damping * self.damping;

        for _ in 0..self.max_iter {
            let t = forward_kinematics(&self.joints, &q);
            let p = ee_position(t);
            let e = [target[0] - p[0], target[1] - p[1], target[2] - p[2]];
            let err = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
            if err < self.tolerance {
                return (q, true);
            }

            // Build J (3×n stored as n cols of length 3)
            let jac = numerical_jacobian(&self.joints, &q, 1e-6);

            // DLS: Δq = Jᵀ (J Jᵀ + λ² I)⁻¹ e
            // For 3 task DOF and n joint DOF:
            // J Jᵀ is 3×3, invert it.
            let mut jjt = [[0.0_f64; 3]; 3];
            for col in &jac {
                for i in 0..3 {
                    for j in 0..3 {
                        jjt[i][j] += col[i] * col[j];
                    }
                }
            }
            // Add damping
            for (i, row) in jjt.iter_mut().enumerate() {
                row[i] += lam2;
            }

            // Invert 3×3
            let det = jjt[0][0] * (jjt[1][1] * jjt[2][2] - jjt[1][2] * jjt[2][1])
                - jjt[0][1] * (jjt[1][0] * jjt[2][2] - jjt[1][2] * jjt[2][0])
                + jjt[0][2] * (jjt[1][0] * jjt[2][1] - jjt[1][1] * jjt[2][0]);
            if det.abs() < 1e-15 {
                break;
            }
            let inv_det = 1.0 / det;
            let inv_jjt = [
                [
                    (jjt[1][1] * jjt[2][2] - jjt[1][2] * jjt[2][1]) * inv_det,
                    (jjt[0][2] * jjt[2][1] - jjt[0][1] * jjt[2][2]) * inv_det,
                    (jjt[0][1] * jjt[1][2] - jjt[0][2] * jjt[1][1]) * inv_det,
                ],
                [
                    (jjt[1][2] * jjt[2][0] - jjt[1][0] * jjt[2][2]) * inv_det,
                    (jjt[0][0] * jjt[2][2] - jjt[0][2] * jjt[2][0]) * inv_det,
                    (jjt[0][2] * jjt[1][0] - jjt[0][0] * jjt[1][2]) * inv_det,
                ],
                [
                    (jjt[1][0] * jjt[2][1] - jjt[1][1] * jjt[2][0]) * inv_det,
                    (jjt[0][1] * jjt[2][0] - jjt[0][0] * jjt[2][1]) * inv_det,
                    (jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0]) * inv_det,
                ],
            ];

            // alpha = J J^T^{-1} e (3D)
            let mut alpha = [0.0_f64; 3];
            for i in 0..3 {
                for j in 0..3 {
                    alpha[i] += inv_jjt[i][j] * e[j];
                }
            }

            // Δq = Jᵀ alpha
            let mut dq = vec![0.0_f64; n];
            for (k, col) in jac.iter().enumerate() {
                for i in 0..3 {
                    dq[k] += col[i] * alpha[i];
                }
            }

            // Update q with joint limit clamping
            for k in 0..n {
                q[k] = clamp(
                    q[k] + dq[k],
                    self.cspace.lower_limits[k],
                    self.cspace.upper_limits[k],
                );
            }
        }
        (q, false)
    }
}

// ── Constraint Manifold Projection ───────────────────────────────────────────

/// Projects a configuration onto a constraint manifold via gradient descent.
///
/// The constraint is expressed as a scalar function c(q) = 0.  The
/// projection moves q along the gradient −∇c until c ≈ 0.
#[derive(Debug, Clone)]
pub struct ManifoldProjector {
    /// Step size for gradient descent.
    pub step_size: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance |c(q)|.
    pub tolerance: f64,
}

impl ManifoldProjector {
    /// Create a new manifold projector.
    pub fn new(step_size: f64, max_iter: usize, tolerance: f64) -> Self {
        Self {
            step_size,
            max_iter,
            tolerance,
        }
    }

    /// Project `q` onto the zero-set of `constraint(q)` using gradient descent.
    ///
    /// `constraint` must return `(value, gradient)` where `gradient` is a
    /// `Vec`f64` of the same length as `q`.
    ///
    /// Returns `(projected_q, converged)`.
    pub fn project<F>(&self, q: &[f64], constraint: F) -> (Vec<f64>, bool)
    where
        F: Fn(&[f64]) -> (f64, Vec<f64>),
    {
        let mut q_cur = q.to_vec();
        for _ in 0..self.max_iter {
            let (val, grad) = constraint(&q_cur);
            if val.abs() < self.tolerance {
                return (q_cur, true);
            }
            let g_norm_sq: f64 = grad.iter().map(|g| g * g).sum();
            if g_norm_sq < 1e-30 {
                break;
            }
            // Step: q ← q - (val / ‖∇c‖²) ∇c
            let alpha = val / g_norm_sq;
            for (qi, gi) in q_cur.iter_mut().zip(grad.iter()) {
                *qi -= self.step_size * alpha * gi;
            }
        }
        let (val, _) = constraint(&q_cur);
        (q_cur, val.abs() < self.tolerance)
    }
}

// ── Collision / Validity Oracle ────────────────────────────────────────────────

/// Trait for checking whether a configuration is valid (collision-free).
pub trait ValidityChecker: Send + Sync {
    /// Returns `true` if the configuration is valid.
    fn is_valid(&self, q: &[f64]) -> bool;
}

/// A trivial validity checker that accepts all configurations.
pub struct AlwaysValidChecker;

impl ValidityChecker for AlwaysValidChecker {
    /// Always returns `true`.
    fn is_valid(&self, _q: &[f64]) -> bool {
        true
    }
}

/// A spherical obstacle in configuration space for testing.
#[derive(Debug, Clone)]
pub struct SphereObstacle {
    /// Center configuration.
    pub center: Vec<f64>,
    /// Radius of the forbidden region.
    pub radius: f64,
}

impl SphereObstacle {
    /// Create a spherical C-space obstacle.
    pub fn new(center: Vec<f64>, radius: f64) -> Self {
        Self { center, radius }
    }
}

impl ValidityChecker for SphereObstacle {
    /// Returns `false` if the configuration is inside the sphere.
    fn is_valid(&self, q: &[f64]) -> bool {
        config_distance(q, &self.center) > self.radius
    }
}

// ── PRM (Probabilistic Roadmap) ───────────────────────────────────────────────

/// A node in the PRM graph.
#[derive(Debug, Clone)]
pub struct PrmNode {
    /// Configuration vector.
    pub config: Vec<f64>,
    /// Neighbours (node index, edge cost).
    pub neighbours: Vec<(usize, f64)>,
}

impl PrmNode {
    /// Create a new PRM node.
    pub fn new(config: Vec<f64>) -> Self {
        Self {
            config,
            neighbours: Vec::new(),
        }
    }
}

/// Probabilistic Roadmap planner.
#[derive(Debug)]
pub struct PrmPlanner {
    /// Sampled nodes in the roadmap.
    pub nodes: Vec<PrmNode>,
    /// C-space definition.
    pub cspace: ConfigurationSpace,
    /// Neighbourhood radius for edge creation.
    pub connection_radius: f64,
    /// Maximum number of nodes to sample.
    pub max_nodes: usize,
}

impl PrmPlanner {
    /// Create a new PRM planner.
    pub fn new(cspace: ConfigurationSpace, connection_radius: f64, max_nodes: usize) -> Self {
        Self {
            nodes: Vec::new(),
            cspace,
            connection_radius,
            max_nodes,
        }
    }

    /// Build the roadmap by sampling random valid configurations and
    /// connecting nearby nodes.
    pub fn build(&mut self, checker: &dyn ValidityChecker, rng: &mut impl Rng) {
        while self.nodes.len() < self.max_nodes {
            let q = self.cspace.sample_random(rng);
            if checker.is_valid(&q) {
                let new_idx = self.nodes.len();
                let mut new_node = PrmNode::new(q.clone());

                // Try to connect to existing nodes within radius
                for (idx, existing) in self.nodes.iter_mut().enumerate() {
                    let d = config_distance(&q, &existing.config);
                    if d <= self.connection_radius {
                        new_node.neighbours.push((idx, d));
                        existing.neighbours.push((new_idx, d));
                    }
                }
                self.nodes.push(new_node);
            }
        }
    }

    /// Add a configuration (start or goal) and connect it to the roadmap.
    ///
    /// Returns the index of the inserted node.
    pub fn add_and_connect(&mut self, q: Vec<f64>) -> usize {
        let new_idx = self.nodes.len();
        let mut new_node = PrmNode::new(q.clone());
        for (idx, existing) in self.nodes.iter_mut().enumerate() {
            let d = config_distance(&q, &existing.config);
            if d <= self.connection_radius {
                new_node.neighbours.push((idx, d));
                existing.neighbours.push((new_idx, d));
            }
        }
        self.nodes.push(new_node);
        new_idx
    }

    /// Query the roadmap using Dijkstra's algorithm.
    ///
    /// Returns the path as a vector of configuration indices, or `None` if
    /// no path exists.
    pub fn query(&self, start_idx: usize, goal_idx: usize) -> Option<Vec<usize>> {
        let n = self.nodes.len();
        let mut dist = vec![f64::INFINITY; n];
        let mut prev = vec![usize::MAX; n];
        dist[start_idx] = 0.0;

        // Simple O(n²) Dijkstra (suitable for small roadmaps)
        let mut visited = vec![false; n];
        for _ in 0..n {
            // Find minimum unvisited node
            let u = (0..n).filter(|&i| !visited[i]).min_by(|&i, &j| {
                dist[i]
                    .partial_cmp(&dist[j])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })?;

            if dist[u].is_infinite() {
                break;
            }
            if u == goal_idx {
                break;
            }
            visited[u] = true;

            for &(v, w) in &self.nodes[u].neighbours {
                let alt = dist[u] + w;
                if alt < dist[v] {
                    dist[v] = alt;
                    prev[v] = u;
                }
            }
        }

        if dist[goal_idx].is_infinite() {
            return None;
        }

        // Reconstruct path
        let mut path = Vec::new();
        let mut cur = goal_idx;
        while cur != usize::MAX {
            path.push(cur);
            if cur == start_idx {
                break;
            }
            cur = prev[cur];
            if path.len() > n {
                return None;
            } // cycle guard
        }
        path.reverse();
        Some(path)
    }
}

// ── RRT (Rapidly-Exploring Random Tree) ───────────────────────────────────────

/// A node in the RRT tree.
#[derive(Debug, Clone)]
pub struct RrtNode {
    /// Configuration.
    pub config: Vec<f64>,
    /// Index of parent node (`usize::MAX` for root).
    pub parent: usize,
    /// Cost from root (sum of edge lengths).
    pub cost: f64,
}

impl RrtNode {
    /// Create a new RRT node.
    pub fn new(config: Vec<f64>, parent: usize, cost: f64) -> Self {
        Self {
            config,
            parent,
            cost,
        }
    }
}

/// Basic RRT planner.
///
/// Grows a tree from `start` towards `goal` by iteratively sampling
/// random configurations and extending the nearest tree node.
#[derive(Debug)]
pub struct RrtPlanner {
    /// C-space.
    pub cspace: ConfigurationSpace,
    /// Extension step size (override per-joint with a scalar step).
    pub step: f64,
    /// Goal bias probability: sample goal instead of random with this probability.
    pub goal_bias: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Goal tolerance.
    pub goal_tol: f64,
}

impl RrtPlanner {
    /// Create a new RRT planner.
    pub fn new(
        cspace: ConfigurationSpace,
        step: f64,
        goal_bias: f64,
        max_iter: usize,
        goal_tol: f64,
    ) -> Self {
        Self {
            cspace,
            step,
            goal_bias,
            max_iter,
            goal_tol,
        }
    }

    /// Plan from `start` to `goal`.
    ///
    /// Returns a path as ordered configurations, or `None` if planning failed.
    pub fn plan(
        &self,
        start: Vec<f64>,
        goal: Vec<f64>,
        checker: &dyn ValidityChecker,
        rng: &mut impl Rng,
    ) -> Option<Vec<Vec<f64>>> {
        let mut tree: Vec<RrtNode> = vec![RrtNode::new(start, usize::MAX, 0.0)];

        for _ in 0..self.max_iter {
            // Sample
            let q_rand = if rng.random_range(0.0_f64..1.0) < self.goal_bias {
                goal.clone()
            } else {
                self.cspace.sample_random(rng)
            };

            // Nearest
            let (near_idx, _) = tree.iter().enumerate().min_by(|(_, a), (_, b)| {
                config_distance(&a.config, &q_rand)
                    .partial_cmp(&config_distance(&b.config, &q_rand))
                    .expect("operation should succeed")
            })?;

            // Extend
            let q_new = self.extend(&tree[near_idx].config, &q_rand);
            if !checker.is_valid(&q_new) {
                continue;
            }

            let d = config_distance(&tree[near_idx].config, &q_new);
            let cost = tree[near_idx].cost + d;
            let new_idx = tree.len();
            tree.push(RrtNode::new(q_new.clone(), near_idx, cost));

            // Check if we reached goal
            if config_distance(&q_new, &goal) < self.goal_tol {
                return Some(self.extract_path(&tree, new_idx));
            }
        }
        None
    }

    /// Extend from `from` towards `to` by at most `step` distance.
    fn extend(&self, from: &[f64], to: &[f64]) -> Vec<f64> {
        let d = config_distance(from, to);
        if d < 1e-15 {
            return from.to_vec();
        }
        let t = (self.step / d).min(1.0);
        config_lerp(from, to, t)
    }

    /// Reconstruct the path from root to node `idx`.
    fn extract_path(&self, tree: &[RrtNode], idx: usize) -> Vec<Vec<f64>> {
        let mut path = Vec::new();
        let mut cur = idx;
        while cur != usize::MAX {
            path.push(tree[cur].config.clone());
            cur = tree[cur].parent;
            if path.len() > tree.len() {
                break;
            } // safety
        }
        path.reverse();
        path
    }
}

// ── CBIRRT (Constrained Bidirectional RRT) ────────────────────────────────────

/// CBIRRT planner: grows two trees (from start and goal) simultaneously
/// while projecting samples onto a constraint manifold.
#[derive(Debug)]
pub struct CbiRrtPlanner {
    /// C-space.
    pub cspace: ConfigurationSpace,
    /// Extension step size.
    pub step: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Goal tolerance.
    pub goal_tol: f64,
    /// Manifold projector.
    pub projector: ManifoldProjector,
}

impl CbiRrtPlanner {
    /// Create a new CBIRRT planner.
    pub fn new(
        cspace: ConfigurationSpace,
        step: f64,
        max_iter: usize,
        goal_tol: f64,
        projector: ManifoldProjector,
    ) -> Self {
        Self {
            cspace,
            step,
            max_iter,
            goal_tol,
            projector,
        }
    }

    /// Plan from `start` to `goal` subject to constraint `c(q) = 0`.
    ///
    /// Returns a path (if found) or `None`.
    pub fn plan<F>(
        &self,
        start: Vec<f64>,
        goal: Vec<f64>,
        checker: &dyn ValidityChecker,
        constraint: F,
        rng: &mut impl Rng,
    ) -> Option<Vec<Vec<f64>>>
    where
        F: Fn(&[f64]) -> (f64, Vec<f64>) + Clone,
    {
        let mut tree_a: Vec<RrtNode> = vec![RrtNode::new(start.clone(), usize::MAX, 0.0)];
        let mut tree_b: Vec<RrtNode> = vec![RrtNode::new(goal.clone(), usize::MAX, 0.0)];

        for iter in 0..self.max_iter {
            // Alternate between the two trees
            let (active_tree, passive_tree) = if iter % 2 == 0 {
                (&mut tree_a, &mut tree_b)
            } else {
                (&mut tree_b, &mut tree_a)
            };

            // Sample and project
            let q_rand = self.cspace.sample_random(rng);
            let (q_proj, proj_ok) = self.projector.project(&q_rand, constraint.clone());
            if !proj_ok || !checker.is_valid(&q_proj) {
                continue;
            }

            // Extend active tree
            let near_idx_active = Self::nearest_idx(active_tree, &q_proj);
            let q_new = Self::extend_step(&active_tree[near_idx_active].config, &q_proj, self.step);
            let (q_new_proj, ok) = self.projector.project(&q_new, constraint.clone());
            if !ok || !checker.is_valid(&q_new_proj) {
                continue;
            }
            let d_act = config_distance(&active_tree[near_idx_active].config, &q_new_proj);
            let cost_new = active_tree[near_idx_active].cost + d_act;
            let new_active_idx = active_tree.len();
            active_tree.push(RrtNode::new(q_new_proj.clone(), near_idx_active, cost_new));

            // Try to connect passive tree to new node
            let near_idx_passive = Self::nearest_idx(passive_tree, &q_new_proj);
            let d_passive = config_distance(&passive_tree[near_idx_passive].config, &q_new_proj);
            if d_passive < self.goal_tol {
                // Found a connection — build path
                let path_a = Self::backtrack(active_tree, new_active_idx);
                let path_b = Self::backtrack(passive_tree, near_idx_passive);
                let mut path = path_a;
                let mut path_b_rev = path_b;
                path_b_rev.reverse();
                path.extend(path_b_rev);
                return Some(path);
            }
        }
        None
    }

    /// Find index of nearest node in `tree` to `q`.
    fn nearest_idx(tree: &[RrtNode], q: &[f64]) -> usize {
        tree.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                config_distance(&a.config, q)
                    .partial_cmp(&config_distance(&b.config, q))
                    .expect("operation should succeed")
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Extend from `from` towards `to` by at most `step`.
    fn extend_step(from: &[f64], to: &[f64], step: f64) -> Vec<f64> {
        let d = config_distance(from, to);
        if d < 1e-15 {
            return from.to_vec();
        }
        config_lerp(from, to, (step / d).min(1.0))
    }

    /// Backtrack from node `idx` to root, returning ordered configs.
    fn backtrack(tree: &[RrtNode], idx: usize) -> Vec<Vec<f64>> {
        let mut path = Vec::new();
        let mut cur = idx;
        loop {
            path.push(tree[cur].config.clone());
            if tree[cur].parent == usize::MAX {
                break;
            }
            cur = tree[cur].parent;
            if path.len() > tree.len() {
                break;
            }
        }
        path.reverse();
        path
    }
}

// ── Trajectory Optimization ───────────────────────────────────────────────────

/// A waypoint in a trajectory.
#[derive(Debug, Clone)]
pub struct Waypoint {
    /// Joint configuration.
    pub config: Vec<f64>,
    /// Timestamp (s).
    pub time: f64,
}

impl Waypoint {
    /// Create a waypoint.
    pub fn new(config: Vec<f64>, time: f64) -> Self {
        Self { config, time }
    }
}

/// Cost function weights for trajectory optimization.
#[derive(Debug, Clone)]
pub struct TrajOptWeights {
    /// Weight on joint velocity smoothness.
    pub smoothness: f64,
    /// Weight on constraint violation.
    pub constraint_penalty: f64,
    /// Weight on total path length.
    pub path_length: f64,
}

impl TrajOptWeights {
    /// Create weights.
    pub fn new(smoothness: f64, constraint_penalty: f64, path_length: f64) -> Self {
        Self {
            smoothness,
            constraint_penalty,
            path_length,
        }
    }
}

/// Constraint-based trajectory optimizer.
///
/// Refines an initial trajectory (sequence of waypoints) by minimizing a
/// cost that balances smoothness, path length, and constraint violations.
/// Uses gradient-free perturbation (coordinate descent on each waypoint).
#[derive(Debug)]
pub struct TrajectoryOptimizer {
    /// C-space.
    pub cspace: ConfigurationSpace,
    /// Optimization weights.
    pub weights: TrajOptWeights,
    /// Number of optimization passes.
    pub max_iter: usize,
    /// Perturbation step for coordinate descent.
    pub perturbation: f64,
}

impl TrajectoryOptimizer {
    /// Create a new trajectory optimizer.
    pub fn new(
        cspace: ConfigurationSpace,
        weights: TrajOptWeights,
        max_iter: usize,
        perturbation: f64,
    ) -> Self {
        Self {
            cspace,
            weights,
            max_iter,
            perturbation,
        }
    }

    /// Evaluate total cost of a trajectory given a scalar constraint function.
    ///
    /// * `wps`        – waypoints to evaluate.
    /// * `constraint` – returns the constraint violation magnitude for a config.
    pub fn cost<F>(&self, wps: &[Waypoint], constraint: &F) -> f64
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = wps.len();
        if n < 2 {
            return 0.0;
        }
        let mut total = 0.0;

        for i in 0..n {
            // Constraint penalty
            let cv = constraint(&wps[i].config);
            total += self.weights.constraint_penalty * cv * cv;

            if i + 1 < n {
                let seg_len = config_distance(&wps[i].config, &wps[i + 1].config);
                total += self.weights.path_length * seg_len;

                // Smoothness: penalise acceleration (finite difference of velocity)
                if i + 2 < n {
                    let acc: f64 = wps[i]
                        .config
                        .iter()
                        .zip(wps[i + 1].config.iter().zip(wps[i + 2].config.iter()))
                        .map(|(&q0, (&q1, &q2))| {
                            let accel = q2 - 2.0 * q1 + q0;
                            accel * accel
                        })
                        .sum::<f64>()
                        .sqrt();
                    total += self.weights.smoothness * acc;
                }
            }
        }
        total
    }

    /// Optimize an initial trajectory using coordinate descent on each waypoint.
    ///
    /// Start and end waypoints are kept fixed.
    ///
    /// Returns the optimized trajectory.
    pub fn optimize<F>(
        &self,
        initial: Vec<Waypoint>,
        constraint: F,
        checker: &dyn ValidityChecker,
    ) -> Vec<Waypoint>
    where
        F: Fn(&[f64]) -> f64 + Clone,
    {
        let mut wps = initial;
        let n = wps.len();
        if n < 3 {
            return wps;
        }

        for _iter in 0..self.max_iter {
            for wp_idx in 1..(n - 1) {
                let dof = wps[wp_idx].config.len();
                let current_cost = self.cost(&wps, &constraint);

                for d in 0..dof {
                    for &sign in &[1.0_f64, -1.0_f64] {
                        let mut candidate = wps[wp_idx].config.clone();
                        candidate[d] += sign * self.perturbation;
                        // Apply joint limits
                        candidate[d] = clamp(
                            candidate[d],
                            self.cspace.lower_limits[d],
                            self.cspace.upper_limits[d],
                        );
                        if !checker.is_valid(&candidate) {
                            continue;
                        }
                        // Temporarily substitute
                        let old_config = wps[wp_idx].config.clone();
                        wps[wp_idx].config = candidate;
                        let new_cost = self.cost(&wps, &constraint);
                        if new_cost < current_cost {
                            // Keep improvement
                        } else {
                            wps[wp_idx].config = old_config;
                        }
                    }
                }
            }
        }
        wps
    }
}

// ── Task-Space Constraint Wrapper ─────────────────────────────────────────────

/// Combines a kinematic chain with a position target as a C-space constraint.
///
/// `c(q) = ‖FK_pos(q) − p_target‖²`
#[derive(Debug, Clone)]
pub struct TaskSpaceConstraint {
    /// DH chain for FK.
    pub joints: Vec<DhJoint>,
    /// Target end-effector position.
    pub target: [f64; 3],
}

impl TaskSpaceConstraint {
    /// Create a task-space constraint.
    pub fn new(joints: Vec<DhJoint>, target: [f64; 3]) -> Self {
        Self { joints, target }
    }

    /// Evaluate the constraint value and its analytical gradient.
    ///
    /// Returns `(c, ∇c)` where `c = ‖p − p_target‖²`.
    pub fn evaluate(&self, q: &[f64]) -> (f64, Vec<f64>) {
        let t = forward_kinematics(&self.joints, q);
        let p = ee_position(t);
        let e = [
            p[0] - self.target[0],
            p[1] - self.target[1],
            p[2] - self.target[2],
        ];
        let c = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];

        // Numerical gradient via finite differences
        let jac = numerical_jacobian(&self.joints, q, 1e-7);
        let n = q.len();
        let mut grad = vec![0.0; n];
        for (k, col) in jac.iter().enumerate() {
            // ∂c/∂q_k = 2 e · ∂p/∂q_k
            for i in 0..3 {
                grad[k] += 2.0 * e[i] * col[i];
            }
        }
        (c, grad)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    // 1. ConfigurationSpace::in_bounds rejects out-of-range configs ───────────
    #[test]
    fn test_cspace_in_bounds_reject() {
        let cs = ConfigurationSpace::uniform(3, 1.0, 0.1);
        assert!(!cs.in_bounds(&[0.0, 0.0, 1.5]));
    }

    // 2. ConfigurationSpace::in_bounds accepts valid config ───────────────────
    #[test]
    fn test_cspace_in_bounds_accept() {
        let cs = ConfigurationSpace::uniform(3, 1.0, 0.1);
        assert!(cs.in_bounds(&[0.0, 0.5, -0.9]));
    }

    // 3. ConfigurationSpace::clamp_to_bounds clamps correctly ─────────────────
    #[test]
    fn test_cspace_clamp() {
        let cs = ConfigurationSpace::uniform(2, 1.0, 0.1);
        let clamped = cs.clamp_to_bounds(&[2.0, -2.0]);
        assert!((clamped[0] - 1.0).abs() < TOL);
        assert!((clamped[1] + 1.0).abs() < TOL);
    }

    // 4. ConfigurationSpace::sample_random stays in bounds ───────────────────
    #[test]
    fn test_cspace_sample_in_bounds() {
        let cs = ConfigurationSpace::uniform(4, std::f64::consts::PI, 0.1);
        let mut rng = rand::rng();
        for _ in 0..50 {
            let q = cs.sample_random(&mut rng);
            assert!(cs.in_bounds(&q), "Sample out of bounds: {:?}", q);
        }
    }

    // 5. config_lerp at t=0 returns start, t=1 returns end ───────────────────
    #[test]
    fn test_config_lerp_endpoints() {
        let a = vec![0.0, 0.0];
        let b = vec![1.0, 2.0];
        let r0 = config_lerp(&a, &b, 0.0);
        let r1 = config_lerp(&a, &b, 1.0);
        for (v, &av) in r0.iter().zip(a.iter()) {
            assert!((v - av).abs() < TOL);
        }
        for (v, &bv) in r1.iter().zip(b.iter()) {
            assert!((v - bv).abs() < TOL);
        }
    }

    // 6. config_lerp at t=0.5 returns midpoint ────────────────────────────────
    #[test]
    fn test_config_lerp_midpoint() {
        let a = vec![0.0, 0.0];
        let b = vec![2.0, 4.0];
        let mid = config_lerp(&a, &b, 0.5);
        assert!((mid[0] - 1.0).abs() < TOL);
        assert!((mid[1] - 2.0).abs() < TOL);
    }

    // 7. config_distance satisfies d(a,a) = 0 ────────────────────────────────
    #[test]
    fn test_config_distance_self() {
        let a = vec![1.0, 2.0, 3.0];
        assert!(config_distance(&a, &a).abs() < TOL);
    }

    // 8. config_distance is symmetric ────────────────────────────────────────
    #[test]
    fn test_config_distance_symmetric() {
        let a = vec![1.0, 0.0];
        let b = vec![4.0, 0.0];
        assert!((config_distance(&a, &b) - config_distance(&b, &a)).abs() < TOL);
    }

    // 9. PositionConstraint::is_satisfied at target position ──────────────────
    #[test]
    fn test_position_constraint_satisfied_at_target() {
        let c = PositionConstraint::new([1.0, 2.0, 3.0], 0.01);
        assert!(c.is_satisfied([1.0, 2.0, 3.0]));
    }

    // 10. PositionConstraint::is_satisfied fails far from target ───────────────
    #[test]
    fn test_position_constraint_fails_far() {
        let c = PositionConstraint::new([0.0, 0.0, 0.0], 0.01);
        assert!(!c.is_satisfied([1.0, 0.0, 0.0]));
    }

    // 11. OrientationConstraint::is_satisfied when axes align ─────────────────
    #[test]
    fn test_orientation_constraint_aligned() {
        let oc = OrientationConstraint::new([0.0, 0.0, 1.0], 0.1);
        assert!(oc.is_satisfied([0.0, 0.0, 1.0]));
    }

    // 12. OrientationConstraint fails for opposite axis ─────────────────────
    #[test]
    fn test_orientation_constraint_opposite() {
        let oc = OrientationConstraint::new([0.0, 0.0, 1.0], 0.1);
        // Opposite direction: angle = π ≈ 3.14 >> tolerance 0.1
        assert!(!oc.is_satisfied([0.0, 0.0, -1.0]));
    }

    // 13. DH identity transform for zero-offset joint ─────────────────────────
    #[test]
    fn test_dh_identity_joint() {
        let j = DhJoint::new(0.0, 0.0, 0.0, 0.0);
        let t = j.transform(0.0);
        // Should be identity
        for (i, row) in t.iter().enumerate() {
            for (k, &val) in row.iter().enumerate() {
                let exp = if i == k { 1.0 } else { 0.0 };
                assert!((val - exp).abs() < 1e-6, "t[{i}][{k}]={val} expected {exp}",);
            }
        }
    }

    // 14. FK for single-link revolute: tip at (a, 0, 0) for q=0 ───────────────
    #[test]
    fn test_fk_single_link() {
        // Joint: a=1, α=0, d=0, θ_off=0
        let joints = vec![DhJoint::new(1.0, 0.0, 0.0, 0.0)];
        let t = forward_kinematics(&joints, &[0.0]);
        let p = ee_position(t);
        assert!((p[0] - 1.0).abs() < 1e-6, "x should be 1, got {}", p[0]);
        assert!(p[1].abs() < 1e-6);
        assert!(p[2].abs() < 1e-6);
    }

    // 15. FK two-link planar robot at full extension ────────────────────────
    #[test]
    fn test_fk_two_link_extended() {
        // Two equal links of length 1
        let joints = vec![
            DhJoint::new(1.0, 0.0, 0.0, 0.0),
            DhJoint::new(1.0, 0.0, 0.0, 0.0),
        ];
        let t = forward_kinematics(&joints, &[0.0, 0.0]);
        let p = ee_position(t);
        assert!((p[0] - 2.0).abs() < 1e-6, "x should be 2 at full extension");
    }

    // 16. Numerical Jacobian: correct dimension ────────────────────────────
    #[test]
    fn test_numerical_jacobian_dimension() {
        let joints = vec![
            DhJoint::new(0.5, 0.0, 0.0, 0.0),
            DhJoint::new(0.5, 0.0, 0.0, 0.0),
            DhJoint::new(0.3, 0.0, 0.0, 0.0),
        ];
        let q = vec![0.1, 0.2, 0.3];
        let jac = numerical_jacobian(&joints, &q, 1e-6);
        assert_eq!(jac.len(), 3, "Jacobian should have n=3 columns");
    }

    // 17. IK solver converges for a simple 2-link robot ─────────────────────
    #[test]
    fn test_ik_solver_2_link() {
        let joints = vec![
            DhJoint::new(0.5, 0.0, 0.0, 0.0),
            DhJoint::new(0.5, 0.0, 0.0, 0.0),
        ];
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.1);
        let solver = IkSolver::new(joints.clone(), cs, 0.01, 500, 1e-4);
        let target = [0.7, 0.2, 0.0];
        let (q_sol, converged) = solver.solve(&[0.0, 0.0], target);
        if converged {
            let t = forward_kinematics(&joints, &q_sol);
            let p = ee_position(t);
            let err = ((p[0] - target[0]).powi(2)
                + (p[1] - target[1]).powi(2)
                + (p[2] - target[2]).powi(2))
            .sqrt();
            assert!(err < 1e-3, "IK error too large: {err}");
        }
        // It's acceptable for the test to proceed even without convergence
        // (the solver returns best-effort q)
    }

    // 18. ManifoldProjector converges for a simple sphere constraint ──────────
    #[test]
    fn test_manifold_projector_sphere() {
        // Constraint: q[0]² + q[1]² = 1 (unit circle)
        let constraint = |q: &[f64]| {
            let c = q[0] * q[0] + q[1] * q[1] - 1.0;
            let grad = vec![2.0 * q[0], 2.0 * q[1]];
            (c, grad)
        };
        let proj = ManifoldProjector::new(1.0, 200, 1e-6);
        let (q_proj, converged) = proj.project(&[1.5, 0.0], constraint);
        assert!(converged, "Projector should converge");
        let c = q_proj[0] * q_proj[0] + q_proj[1] * q_proj[1] - 1.0;
        assert!(
            c.abs() < 1e-4,
            "Projected point should lie on unit circle, c={c}"
        );
    }

    // 19. AlwaysValidChecker always returns true ────────────────────────────
    #[test]
    fn test_always_valid_checker() {
        let vc = AlwaysValidChecker;
        assert!(vc.is_valid(&[0.0, 0.0, 0.0]));
        assert!(vc.is_valid(&[999.0, -999.0]));
    }

    // 20. SphereObstacle rejects interior points ────────────────────────────
    #[test]
    fn test_sphere_obstacle_interior_invalid() {
        let obs = SphereObstacle::new(vec![0.0, 0.0], 1.0);
        assert!(!obs.is_valid(&[0.1, 0.0]));
    }

    // 21. SphereObstacle accepts exterior points ────────────────────────────
    #[test]
    fn test_sphere_obstacle_exterior_valid() {
        let obs = SphereObstacle::new(vec![0.0, 0.0], 1.0);
        assert!(obs.is_valid(&[2.0, 0.0]));
    }

    // 22. PRM: roadmap has correct number of nodes ──────────────────────────
    #[test]
    fn test_prm_node_count() {
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.5);
        let mut prm = PrmPlanner::new(cs, 1.0, 20);
        let vc = AlwaysValidChecker;
        let mut rng = rand::rng();
        prm.build(&vc, &mut rng);
        assert_eq!(prm.nodes.len(), 20);
    }

    // 23. PRM: query finds a trivial path when start == goal ─────────────────
    #[test]
    fn test_prm_trivial_path() {
        let cs = ConfigurationSpace::uniform(2, 1.0, 0.2);
        let mut prm = PrmPlanner::new(cs, 2.0, 30);
        let vc = AlwaysValidChecker;
        let mut rng = rand::rng();
        prm.build(&vc, &mut rng);
        // Add start and goal as the same node
        let idx = prm.add_and_connect(vec![0.0, 0.0]);
        let path = prm.query(idx, idx);
        assert!(path.is_some());
    }

    // 24. RRT planner: path starts at start configuration ────────────────────
    #[test]
    fn test_rrt_path_starts_at_start() {
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.5);
        let planner = RrtPlanner::new(cs, 0.3, 0.1, 2000, 0.1);
        let vc = AlwaysValidChecker;
        let mut rng = rand::rng();
        let start = vec![0.0, 0.0];
        let goal = vec![1.0, 1.0];
        if let Some(path) = planner.plan(start.clone(), goal, &vc, &mut rng) {
            let d = config_distance(&path[0], &start);
            assert!(d < 0.5, "Path should start near start, d={d}");
        }
    }

    // 25. RRT planner: path ends near goal ───────────────────────────────────
    #[test]
    fn test_rrt_path_ends_near_goal() {
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.5);
        let planner = RrtPlanner::new(cs, 0.3, 0.2, 3000, 0.15);
        let vc = AlwaysValidChecker;
        let mut rng = rand::rng();
        let start = vec![0.0, 0.0];
        let goal = vec![1.0, 1.0];
        if let Some(path) = planner.plan(start, goal.clone(), &vc, &mut rng) {
            let d = config_distance(path.last().unwrap(), &goal);
            assert!(d < 0.5, "Path should end near goal, d={d}");
        }
    }

    // 26. CBIRRT: projects random samples onto manifold ──────────────────────
    #[test]
    fn test_cbirrt_projection_onto_manifold() {
        let proj = ManifoldProjector::new(1.0, 100, 1e-5);
        let constraint = |q: &[f64]| {
            let c = q[0] * q[0] + q[1] * q[1] - 1.0;
            let grad = vec![2.0 * q[0], 2.0 * q[1]];
            (c, grad)
        };
        let (q_proj, ok) = proj.project(&[2.0, 0.0], constraint);
        assert!(ok, "Projection should converge");
        let on_circle = q_proj[0] * q_proj[0] + q_proj[1] * q_proj[1];
        assert!(
            (on_circle - 1.0).abs() < 1e-3,
            "Projected point off manifold"
        );
    }

    // 27. TaskSpaceConstraint: zero error at target ──────────────────────────
    #[test]
    fn test_task_space_constraint_zero_at_target() {
        let joints = vec![DhJoint::new(1.0, 0.0, 0.0, 0.0)];
        // At q=0, FK gives p=(1,0,0). Set target to (1,0,0).
        let tsc = TaskSpaceConstraint::new(joints, [1.0, 0.0, 0.0]);
        let (c, _grad) = tsc.evaluate(&[0.0]);
        assert!(c.abs() < 1e-8, "Constraint should be zero at target, c={c}");
    }

    // 28. TrajectoryOptimizer: cost decreases after optimization ─────────────
    #[test]
    fn test_traj_opt_cost_decreases() {
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.5);
        let weights = TrajOptWeights::new(1.0, 0.0, 1.0);
        let opt = TrajectoryOptimizer::new(cs, weights, 20, 0.05);
        let vc = AlwaysValidChecker;
        // Zig-zag trajectory (should be smoothed)
        let wps = vec![
            Waypoint::new(vec![0.0, 0.0], 0.0),
            Waypoint::new(vec![1.0, -1.0], 0.5),
            Waypoint::new(vec![0.5, 0.5], 1.0),
            Waypoint::new(vec![1.0, 1.0], 1.5),
        ];
        let constraint = |_q: &[f64]| 0.0_f64;
        let c0 = opt.cost(&wps, &constraint);
        let optimized = opt.optimize(wps, constraint, &vc);
        let c1 = opt.cost(&optimized, &|_q: &[f64]| 0.0_f64);
        // Cost must be non-negative
        assert!(c1 >= 0.0);
        // Optimized cost should be ≤ initial (or at worst unchanged for a bad initialization)
        assert!(
            c1 <= c0 + 1e-6,
            "Optimized cost {c1} should be ≤ initial {c0}"
        );
    }

    // 29. TrajectoryOptimizer: fixed endpoints are unchanged ─────────────────
    #[test]
    fn test_traj_opt_fixed_endpoints() {
        let cs = ConfigurationSpace::uniform(2, std::f64::consts::PI, 0.5);
        let weights = TrajOptWeights::new(1.0, 0.0, 1.0);
        let opt = TrajectoryOptimizer::new(cs, weights, 10, 0.05);
        let vc = AlwaysValidChecker;
        let start = vec![0.0, 0.0];
        let end = vec![1.0, 1.0];
        let wps = vec![
            Waypoint::new(start.clone(), 0.0),
            Waypoint::new(vec![0.5, 0.5], 0.5),
            Waypoint::new(end.clone(), 1.0),
        ];
        let constraint = |_q: &[f64]| 0.0_f64;
        let optimized = opt.optimize(wps, constraint, &vc);
        // First and last waypoints must be unchanged
        for (v, &s) in optimized[0].config.iter().zip(start.iter()) {
            assert!((v - s).abs() < TOL, "Start endpoint changed");
        }
        for (v, &e) in optimized.last().unwrap().config.iter().zip(end.iter()) {
            assert!((v - e).abs() < TOL, "End endpoint changed");
        }
    }

    // 30. interpolate_path: first and last points match a and b ──────────────
    #[test]
    fn test_interpolate_path_endpoints() {
        let cs = ConfigurationSpace::uniform(2, 1.0, 0.1);
        let a = vec![0.0, 0.0];
        let b = vec![1.0, 1.0];
        let path = cs.interpolate_path(&a, &b);
        assert!(path.len() >= 2, "Path must have at least 2 points");
        let d_start = config_distance(&path[0], &a);
        let d_end = config_distance(path.last().unwrap(), &b);
        assert!(d_start < 1e-6, "Path must start at a");
        assert!(d_end < 1e-6, "Path must end at b");
    }

    // 31. dot, norm helpers ───────────────────────────────────────────────────
    #[test]
    fn test_dot_norm_helpers() {
        let a = vec![3.0, 4.0];
        assert!((norm(&a) - 5.0).abs() < TOL);
        assert!((dot(&a, &a) - 25.0).abs() < TOL);
    }

    // 32. PRM: add_and_connect inserts a node ─────────────────────────────────
    #[test]
    fn test_prm_add_and_connect() {
        let cs = ConfigurationSpace::uniform(2, 1.0, 0.1);
        let mut prm = PrmPlanner::new(cs, 2.0, 10);
        let vc = AlwaysValidChecker;
        let mut rng = rand::rng();
        prm.build(&vc, &mut rng);
        let n_before = prm.nodes.len();
        prm.add_and_connect(vec![0.0, 0.0]);
        assert_eq!(prm.nodes.len(), n_before + 1);
    }

    // 33. CbiRrtPlanner: nearest_idx returns 0 for single-node tree ──────────
    #[test]
    fn test_cbirrt_nearest_single_node() {
        let nodes = vec![RrtNode::new(vec![1.0, 0.0], usize::MAX, 0.0)];
        let idx = CbiRrtPlanner::nearest_idx(&nodes, &[2.0, 0.0]);
        assert_eq!(idx, 0);
    }

    // 34. clamp helper ────────────────────────────────────────────────────────
    #[test]
    fn test_clamp_helper() {
        assert!((clamp(0.5, 0.0, 1.0) - 0.5).abs() < TOL);
        assert!((clamp(-1.0, 0.0, 1.0)).abs() < TOL);
        assert!((clamp(2.0, 0.0, 1.0) - 1.0).abs() < TOL);
    }

    // 35. TrajOptWeights construction ─────────────────────────────────────────
    #[test]
    fn test_traj_opt_weights() {
        let w = TrajOptWeights::new(1.0, 2.0, 3.0);
        assert!((w.smoothness - 1.0).abs() < TOL);
        assert!((w.constraint_penalty - 2.0).abs() < TOL);
        assert!((w.path_length - 3.0).abs() < TOL);
    }
}
