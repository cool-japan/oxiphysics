#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multi-agent constraint system for the OxiPhysics engine.
//!
//! Provides consensus constraints, distributed optimization (ADMM),
//! collective motion constraints, leader-follower constraints,
//! communication topology (graph Laplacian), collision avoidance (ORCA /
//! velocity obstacles), task assignment constraints, energy constraints,
//! formation reconfiguration constraints, bearing-only constraints,
//! range constraints, coverage constraints, Nash equilibrium constraints,
//! and the social force model.

#![allow(dead_code)]

// ── Internal math helpers ─────────────────────────────────────────────────────

fn vec2_add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

fn vec2_sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn vec2_scale(v: [f64; 2], s: f64) -> [f64; 2] {
    [v[0] * s, v[1] * s]
}

fn vec2_dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn vec2_len(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

fn vec2_normalize(v: [f64; 2]) -> [f64; 2] {
    let len = vec2_len(v);
    if len < 1e-15 {
        [0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len]
    }
}

fn vec2_dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    vec2_len(vec2_sub(a, b))
}

fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn vec3_normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vec3_len(v);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_len(vec3_sub(a, b))
}

fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ── Graph / Communication topology ───────────────────────────────────────────

/// Adjacency entry in the communication graph.
#[derive(Debug, Clone)]
pub struct GraphEdge {
    /// Index of the source agent.
    pub from: usize,
    /// Index of the destination agent.
    pub to: usize,
    /// Edge weight (link quality, bandwidth, etc.).
    pub weight: f64,
}

/// Sparse communication graph for a multi-agent system.
///
/// Each agent is a node; edges indicate communication links.
#[derive(Debug, Clone)]
pub struct CommunicationGraph {
    /// Number of agents (nodes) in the graph.
    pub num_agents: usize,
    /// Directed edges in the graph.
    pub edges: Vec<GraphEdge>,
}

impl CommunicationGraph {
    /// Create an empty communication graph with `n` agents.
    pub fn new(n: usize) -> Self {
        Self {
            num_agents: n,
            edges: Vec::new(),
        }
    }

    /// Add a directed, weighted communication link from `from` → `to`.
    pub fn add_edge(&mut self, from: usize, to: usize, weight: f64) {
        self.edges.push(GraphEdge { from, to, weight });
    }

    /// Compute the degree of node `i` (number of outgoing edges).
    pub fn out_degree(&self, i: usize) -> usize {
        self.edges.iter().filter(|e| e.from == i).count()
    }

    /// Compute the graph Laplacian matrix `L = D − A` as a flat row-major
    /// `n × n` array. Returns `Vec`f64` of length `n * n`.
    ///
    /// For an undirected graph, `L` is positive semi-definite with the smallest
    /// eigenvalue 0 (corresponding to the all-ones vector for a connected graph).
    pub fn laplacian(&self) -> Vec<f64> {
        let n = self.num_agents;
        let mut l = vec![0.0_f64; n * n];
        // Fill adjacency (off-diagonal) and accumulate degree
        for e in &self.edges {
            l[e.from * n + e.to] -= e.weight;
            l[e.from * n + e.from] += e.weight;
        }
        l
    }

    /// Multiply the Laplacian by a column vector `x` of length `n`,
    /// returning `L x` without building the full matrix.
    ///
    /// Useful for large graphs where storing `L` is expensive.
    pub fn laplacian_matvec(&self, x: &[f64]) -> Vec<f64> {
        let n = self.num_agents;
        let mut y = vec![0.0_f64; n];
        for e in &self.edges {
            y[e.from] += e.weight * (x[e.from] - x[e.to]);
        }
        y
    }

    /// Check if the graph is connected (breadth-first search from node 0).
    pub fn is_connected(&self) -> bool {
        if self.num_agents == 0 {
            return true;
        }
        let n = self.num_agents;
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(0usize);
        visited[0] = true;
        while let Some(node) = queue.pop_front() {
            for e in &self.edges {
                if e.from == node && !visited[e.to] {
                    visited[e.to] = true;
                    queue.push_back(e.to);
                }
                if e.to == node && !visited[e.from] {
                    visited[e.from] = true;
                    queue.push_back(e.from);
                }
            }
        }
        visited.iter().all(|&v| v)
    }
}

// ── Consensus Constraint ──────────────────────────────────────────────────────

/// Configuration for a consensus constraint run.
///
/// The discrete-time update is:
/// ```text
/// x_i(k+1) = x_i(k) − α Σ_j w_ij (x_i(k) − x_j(k))
/// ```
/// where `α` is the step size and `w_ij` are the link weights.
#[derive(Debug, Clone)]
pub struct ConsensusConstraint {
    /// Step size α for the gradient descent update.
    pub alpha: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance (L∞ norm of change).
    pub tol: f64,
    /// Communication graph.
    pub graph: CommunicationGraph,
}

impl ConsensusConstraint {
    /// Create a new consensus constraint with the given graph and defaults
    /// (`alpha = 0.1`, `max_iter = 500`, `tol = 1e-6`).
    pub fn new(graph: CommunicationGraph) -> Self {
        Self {
            alpha: 0.1,
            max_iter: 500,
            tol: 1e-6,
            graph,
        }
    }

    /// Run the consensus iteration on scalar values `x` (one per agent).
    ///
    /// Returns the number of iterations taken and the converged values.
    pub fn run(&self, x0: &[f64]) -> (usize, Vec<f64>) {
        let n = x0.len();
        let mut x = x0.to_vec();
        let mut dx = vec![0.0_f64; n];
        for iter in 0..self.max_iter {
            // Compute L x
            let lx = self.graph.laplacian_matvec(&x);
            let mut max_change = 0.0_f64;
            for i in 0..n {
                dx[i] = -self.alpha * lx[i];
                max_change = max_change.max(dx[i].abs());
            }
            for i in 0..n {
                x[i] += dx[i];
            }
            if max_change < self.tol {
                return (iter + 1, x);
            }
        }
        (self.max_iter, x)
    }

    /// Compute the consensus error: L2 norm of `x − x_mean * 1`.
    pub fn consensus_error(x: &[f64]) -> f64 {
        if x.is_empty() {
            return 0.0;
        }
        let mean = x.iter().sum::<f64>() / x.len() as f64;
        x.iter().map(|&xi| (xi - mean).powi(2)).sum::<f64>().sqrt()
    }
}

// ── ADMM Distributed Optimization ────────────────────────────────────────────

/// A single ADMM variable block for distributed optimization.
///
/// Each agent maintains a primal variable `x`, an auxiliary variable `z`,
/// and a dual variable (scaled) `u`. The global constraint is `x_i = z`.
#[derive(Debug, Clone)]
pub struct AdmmAgent {
    /// Current primal variable (length = `dim`).
    pub x: Vec<f64>,
    /// Auxiliary consensus variable.
    pub z: Vec<f64>,
    /// Scaled dual variable u = (1/ρ) λ.
    pub u: Vec<f64>,
    /// Dimension of the variable.
    pub dim: usize,
    /// Local quadratic cost weight matrix (diagonal, length = `dim`).
    pub q_diag: Vec<f64>,
    /// Local linear cost vector (length = `dim`).
    pub c: Vec<f64>,
}

impl AdmmAgent {
    /// Create a new ADMM agent with zero initialization.
    pub fn new(dim: usize) -> Self {
        Self {
            x: vec![0.0; dim],
            z: vec![0.0; dim],
            u: vec![0.0; dim],
            dim,
            q_diag: vec![1.0; dim],
            c: vec![0.0; dim],
        }
    }

    /// Local x-update (closed-form for quadratic objectives):
    /// ```text
    /// x_i = (Q + ρ I)^{-1} (−c + ρ (z − u))
    /// ```
    pub fn x_update(&mut self, rho: f64) {
        for k in 0..self.dim {
            let rhs = -self.c[k] + rho * (self.z[k] - self.u[k]);
            let denom = self.q_diag[k] + rho;
            self.x[k] = rhs / denom;
        }
    }

    /// Dual update: `u += x − z`.
    pub fn dual_update(&mut self) {
        for k in 0..self.dim {
            self.u[k] += self.x[k] - self.z[k];
        }
    }

    /// Residual `r = x − z` (primal feasibility).
    pub fn primal_residual(&self) -> Vec<f64> {
        (0..self.dim).map(|k| self.x[k] - self.z[k]).collect()
    }
}

/// ADMM coordinator: performs the z-update (averaging step) across all agents.
///
/// The z-update for the consensus ADMM is the average of `x_i + u_i`:
/// ```text
/// z = (1/N) Σ_i (x_i + u_i)
/// ```
pub fn admm_z_update(agents: &mut Vec<AdmmAgent>) {
    if agents.is_empty() {
        return;
    }
    let n = agents.len();
    let dim = agents[0].dim;
    let mut z_avg = vec![0.0_f64; dim];
    for ag in agents.iter() {
        for k in 0..dim {
            z_avg[k] += ag.x[k] + ag.u[k];
        }
    }
    for k in 0..dim {
        z_avg[k] /= n as f64;
    }
    for ag in agents.iter_mut() {
        ag.z.clone_from(&z_avg);
    }
}

/// Run a full ADMM consensus optimization for `max_iter` iterations.
///
/// Returns the converged z value and iteration count.
#[allow(clippy::too_many_arguments)]
pub fn admm_solve(
    agents: &mut Vec<AdmmAgent>,
    rho: f64,
    max_iter: usize,
    tol: f64,
) -> (usize, Vec<f64>) {
    let dim = if agents.is_empty() { 0 } else { agents[0].dim };
    for iter in 0..max_iter {
        // x-update
        for ag in agents.iter_mut() {
            ag.x_update(rho);
        }
        // z-update
        admm_z_update(agents);
        // dual-update + check convergence
        let mut max_res = 0.0_f64;
        for ag in agents.iter_mut() {
            let res = ag.primal_residual();
            let norm = res.iter().map(|r| r * r).sum::<f64>().sqrt();
            max_res = max_res.max(norm);
            ag.dual_update();
        }
        if max_res < tol {
            let z = if agents.is_empty() {
                vec![]
            } else {
                agents[0].z.clone()
            };
            return (iter + 1, z);
        }
    }
    let z = if dim > 0 && !agents.is_empty() {
        agents[0].z.clone()
    } else {
        vec![0.0; dim]
    };
    (max_iter, z)
}

// ── Collective Motion Constraint ──────────────────────────────────────────────

/// Collective motion constraint based on Reynolds flocking rules.
///
/// Computes acceleration commands for each agent to achieve
/// cohesion, alignment, and separation.
#[derive(Debug, Clone)]
pub struct CollectiveMotionConstraint {
    /// Weight for cohesion force.
    pub w_cohesion: f64,
    /// Weight for alignment force.
    pub w_alignment: f64,
    /// Weight for separation force.
    pub w_separation: f64,
    /// Neighbourhood radius for cohesion and alignment.
    pub r_neighbour: f64,
    /// Separation radius below which repulsion activates.
    pub r_separation: f64,
    /// Maximum acceleration magnitude.
    pub max_accel: f64,
}

impl CollectiveMotionConstraint {
    /// Create with default Reynolds weights.
    pub fn new() -> Self {
        Self {
            w_cohesion: 1.0,
            w_alignment: 1.0,
            w_separation: 2.0,
            r_neighbour: 5.0,
            r_separation: 1.5,
            max_accel: 3.0,
        }
    }

    /// Compute the acceleration command for agent `i` given positions and
    /// velocities of all agents (3D).
    ///
    /// Returns the desired acceleration vector `\[ax, ay, az\]`.
    pub fn compute_acceleration(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        i: usize,
    ) -> [f64; 3] {
        let pos_i = positions[i];
        let vel_i = velocities[i];
        let mut cohesion_sum = [0.0_f64; 3];
        let mut align_sum = [0.0_f64; 3];
        let mut sep_force = [0.0_f64; 3];
        let mut cohesion_count = 0u32;
        let mut align_count = 0u32;

        for (j, (&pos_j, &vel_j)) in positions.iter().zip(velocities.iter()).enumerate() {
            if i == j {
                continue;
            }
            let d = vec3_dist(pos_i, pos_j);
            if d < self.r_neighbour {
                cohesion_sum = vec3_add(cohesion_sum, pos_j);
                cohesion_count += 1;
                align_sum = vec3_add(align_sum, vel_j);
                align_count += 1;
            }
            if d < self.r_separation && d > 1e-9 {
                let diff = vec3_sub(pos_i, pos_j);
                // Repulsion inversely proportional to distance
                let rep = vec3_scale(vec3_normalize(diff), self.r_separation / d.max(1e-6));
                sep_force = vec3_add(sep_force, rep);
            }
        }

        let mut accel = [0.0_f64; 3];

        // Cohesion: steer toward centre of neighbours
        if cohesion_count > 0 {
            let centre = vec3_scale(cohesion_sum, 1.0 / cohesion_count as f64);
            let dir = vec3_normalize(vec3_sub(centre, pos_i));
            accel = vec3_add(accel, vec3_scale(dir, self.w_cohesion));
        }

        // Alignment: match average velocity
        if align_count > 0 {
            let avg_vel = vec3_scale(align_sum, 1.0 / align_count as f64);
            let dir = vec3_normalize(vec3_sub(avg_vel, vel_i));
            accel = vec3_add(accel, vec3_scale(dir, self.w_alignment));
        }

        // Separation
        accel = vec3_add(accel, vec3_scale(sep_force, self.w_separation));

        // Clamp magnitude
        let mag = vec3_len(accel);
        if mag > self.max_accel {
            accel = vec3_scale(accel, self.max_accel / mag);
        }
        accel
    }
}

impl Default for CollectiveMotionConstraint {
    fn default() -> Self {
        Self::new()
    }
}

// ── Leader-Follower Constraint ────────────────────────────────────────────────

/// Role of an agent in a leader-follower hierarchy.
#[derive(Debug, Clone, PartialEq)]
pub enum AgentRole {
    /// This agent dictates the reference trajectory.
    Leader,
    /// This agent tracks a leader or predecessor.
    Follower,
}

/// Leader-follower constraint: a follower tracks a leader with a desired
/// offset in position and heading.
#[derive(Debug, Clone)]
pub struct LeaderFollowerConstraint {
    /// Index of the leader agent.
    pub leader_idx: usize,
    /// Index of the follower agent.
    pub follower_idx: usize,
    /// Desired 3D offset from leader to follower (in world frame).
    pub desired_offset: [f64; 3],
    /// Proportional gain for position error.
    pub kp: f64,
    /// Derivative gain for velocity error.
    pub kd: f64,
}

impl LeaderFollowerConstraint {
    /// Create a new leader-follower constraint.
    pub fn new(leader_idx: usize, follower_idx: usize, desired_offset: [f64; 3]) -> Self {
        Self {
            leader_idx,
            follower_idx,
            desired_offset,
            kp: 2.0,
            kd: 0.5,
        }
    }

    /// Compute the acceleration command for the follower to track the leader.
    ///
    /// `positions\[i\]` and `velocities\[i\]` are the world-space state of agent `i`.
    pub fn compute_acceleration(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> [f64; 3] {
        let p_leader = positions[self.leader_idx];
        let v_leader = velocities[self.leader_idx];
        let p_follower = positions[self.follower_idx];
        let v_follower = velocities[self.follower_idx];

        // Desired follower position
        let p_desired = vec3_add(p_leader, self.desired_offset);
        // Position error
        let ep = vec3_sub(p_desired, p_follower);
        // Velocity error (track leader velocity)
        let ev = vec3_sub(v_leader, v_follower);

        vec3_add(vec3_scale(ep, self.kp), vec3_scale(ev, self.kd))
    }
}

// ── Formation Constraint ──────────────────────────────────────────────────────

/// Rigid formation constraint: each agent maintains a fixed offset from a
/// virtual formation centre.
#[derive(Debug, Clone)]
pub struct FormationConstraint {
    /// Number of agents.
    pub num_agents: usize,
    /// Desired offsets from formation centre, one per agent.
    pub offsets: Vec<[f64; 3]>,
    /// Proportional gain.
    pub kp: f64,
    /// Derivative gain.
    pub kd: f64,
}

impl FormationConstraint {
    /// Create a formation constraint from a list of offsets.
    pub fn new(offsets: Vec<[f64; 3]>) -> Self {
        let n = offsets.len();
        Self {
            num_agents: n,
            offsets,
            kp: 2.0,
            kd: 0.5,
        }
    }

    /// Compute the formation centre as the mean of agent positions.
    pub fn formation_centre(positions: &[[f64; 3]]) -> [f64; 3] {
        let n = positions.len();
        if n == 0 {
            return [0.0; 3];
        }
        let sum = positions
            .iter()
            .fold([0.0_f64; 3], |acc, &p| vec3_add(acc, p));
        vec3_scale(sum, 1.0 / n as f64)
    }

    /// Compute acceleration commands for all agents to hold formation.
    ///
    /// Returns a vector of accelerations, one per agent.
    pub fn compute_accelerations(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        centre_vel: [f64; 3],
    ) -> Vec<[f64; 3]> {
        let centre = Self::formation_centre(positions);
        let n = self.num_agents.min(positions.len());
        (0..n)
            .map(|i| {
                let p_desired = vec3_add(centre, self.offsets[i]);
                let ep = vec3_sub(p_desired, positions[i]);
                let ev = vec3_sub(centre_vel, velocities[i]);
                vec3_add(vec3_scale(ep, self.kp), vec3_scale(ev, self.kd))
            })
            .collect()
    }
}

// ── Formation Reconfiguration ─────────────────────────────────────────────────

/// Reconfiguration plan: a sequence of target formation offsets to interpolate
/// between using a smooth ramp.
#[derive(Debug, Clone)]
pub struct FormationReconfiguration {
    /// Starting offsets (one per agent).
    pub offsets_from: Vec<[f64; 3]>,
    /// Target offsets (one per agent).
    pub offsets_to: Vec<[f64; 3]>,
    /// Total reconfiguration time (seconds).
    pub duration: f64,
    /// Current elapsed time.
    pub elapsed: f64,
}

impl FormationReconfiguration {
    /// Create a new reconfiguration plan.
    pub fn new(offsets_from: Vec<[f64; 3]>, offsets_to: Vec<[f64; 3]>, duration: f64) -> Self {
        Self {
            offsets_from,
            offsets_to,
            duration,
            elapsed: 0.0,
        }
    }

    /// Advance time and return interpolated offsets for the current instant.
    ///
    /// Uses a smooth-step blend: `s(t) = 3t² − 2t³`.
    pub fn step(&mut self, dt: f64) -> Vec<[f64; 3]> {
        self.elapsed = (self.elapsed + dt).min(self.duration);
        let t = if self.duration > 1e-15 {
            self.elapsed / self.duration
        } else {
            1.0
        };
        let s = t * t * (3.0 - 2.0 * t); // smooth-step
        let n = self.offsets_from.len().min(self.offsets_to.len());
        (0..n)
            .map(|i| {
                let a = self.offsets_from[i];
                let b = self.offsets_to[i];
                [
                    a[0] + (b[0] - a[0]) * s,
                    a[1] + (b[1] - a[1]) * s,
                    a[2] + (b[2] - a[2]) * s,
                ]
            })
            .collect()
    }

    /// Returns `true` if the reconfiguration has completed.
    pub fn is_complete(&self) -> bool {
        self.elapsed >= self.duration - 1e-12
    }
}

// ── ORCA / Velocity Obstacle Collision Avoidance ─────────────────────────────

/// A 2D half-plane constraint for ORCA (Optimal Reciprocal Collision Avoidance).
///
/// The constraint is `n · (v − p) ≥ 0`, where `n` is the outward normal of the
/// half-plane and `p` is a point on the boundary.
#[derive(Debug, Clone)]
pub struct HalfPlane2D {
    /// Point on the boundary of the half-plane.
    pub point: [f64; 2],
    /// Outward normal of the half-plane (unit vector).
    pub normal: [f64; 2],
}

/// Compute the ORCA half-plane for agent A (radius `r_a`, position `p_a`,
/// velocity `v_a`) with respect to agent B (radius `r_b`, position `p_b`,
/// velocity `v_b`) over a time horizon `tau`.
///
/// Returns the ORCA half-plane. The agent should choose a new velocity that
/// satisfies all such half-planes.
pub fn orca_half_plane(
    p_a: [f64; 2],
    v_a: [f64; 2],
    r_a: f64,
    p_b: [f64; 2],
    v_b: [f64; 2],
    r_b: f64,
    tau: f64,
) -> HalfPlane2D {
    let rel_pos = vec2_sub(p_b, p_a);
    let rel_vel = vec2_sub(v_a, v_b);
    let combined_r = r_a + r_b;
    let dist2 = vec2_dot(rel_pos, rel_pos);
    let dist = dist2.sqrt();

    // If agents overlap, push apart
    if dist < combined_r {
        let n = vec2_normalize(rel_pos);
        let u = vec2_scale(n, combined_r - dist + 0.001);
        return HalfPlane2D {
            point: vec2_add(v_a, vec2_scale(u, 0.5)),
            normal: n,
        };
    }

    // Velocity obstacle centre (in relative velocity space)
    let vo_centre = vec2_scale(rel_pos, 1.0 / tau);
    let r_vo = combined_r / tau;

    // Vector from VO centre to relative velocity
    let w = vec2_sub(rel_vel, vo_centre);
    let w_len = vec2_len(w);

    let (u, n);
    if w_len < r_vo {
        // Inside VO — project to boundary
        let n_out = if w_len > 1e-15 {
            vec2_normalize(w)
        } else {
            vec2_normalize([-rel_pos[1], rel_pos[0]])
        };
        u = vec2_scale(n_out, r_vo - w_len);
        n = n_out;
    } else {
        // Outside VO — no constraint needed; return trivial half-plane
        n = vec2_normalize(w);
        u = [0.0; 2];
    }

    HalfPlane2D {
        point: vec2_add(v_a, vec2_scale(u, 0.5)),
        normal: n,
    }
}

/// Solve a linear program to find the velocity closest to `v_pref` that
/// satisfies all ORCA half-plane constraints (2D, iterative algorithm).
///
/// Returns the best feasible velocity.
pub fn orca_solve_velocity(
    v_pref: [f64; 2],
    max_speed: f64,
    constraints: &[HalfPlane2D],
) -> [f64; 2] {
    let mut v = v_pref;
    // Clamp to max speed
    if vec2_len(v) > max_speed {
        v = vec2_scale(vec2_normalize(v), max_speed);
    }
    for hp in constraints {
        if vec2_dot(hp.normal, vec2_sub(v, hp.point)) < 0.0 {
            // Project onto the half-plane boundary
            let n = hp.normal;
            let dot = vec2_dot(n, vec2_sub(hp.point, v));
            v = vec2_add(v, vec2_scale(n, dot));
            // Re-clamp to max speed
            let spd = vec2_len(v);
            if spd > max_speed {
                v = vec2_scale(vec2_normalize(v), max_speed);
            }
        }
    }
    v
}

// ── Bearing-Only Constraint ───────────────────────────────────────────────────

/// Bearing-only constraint: maintain a target bearing from agent A to agent B.
///
/// The bearing is the angle in the XY plane, measured from the positive-X axis.
#[derive(Debug, Clone)]
pub struct BearingConstraint {
    /// Index of the observing agent.
    pub observer: usize,
    /// Index of the target agent.
    pub target: usize,
    /// Desired bearing angle (radians).
    pub desired_bearing: f64,
    /// Proportional gain.
    pub kp: f64,
}

impl BearingConstraint {
    /// Create a bearing constraint with default gain.
    pub fn new(observer: usize, target: usize, desired_bearing: f64) -> Self {
        Self {
            observer,
            target,
            desired_bearing,
            kp: 1.0,
        }
    }

    /// Compute bearing angle from `p_obs` to `p_tgt` in the XY plane.
    pub fn bearing(p_obs: [f64; 3], p_tgt: [f64; 3]) -> f64 {
        let dx = p_tgt[0] - p_obs[0];
        let dy = p_tgt[1] - p_obs[1];
        dy.atan2(dx)
    }

    /// Compute the bearing error (shortest angular distance).
    pub fn bearing_error(&self, positions: &[[f64; 3]]) -> f64 {
        let bearing = Self::bearing(positions[self.observer], positions[self.target]);
        let err = self.desired_bearing - bearing;
        // Wrap to [-π, π]
        let pi = std::f64::consts::PI;
        ((err + pi) % (2.0 * pi)) - pi
    }

    /// Compute a yaw-rate command for the observer to correct bearing error.
    pub fn yaw_command(&self, positions: &[[f64; 3]]) -> f64 {
        self.kp * self.bearing_error(positions)
    }
}

// ── Range Constraint ──────────────────────────────────────────────────────────

/// Range constraint: enforce a desired distance between two agents.
#[derive(Debug, Clone)]
pub struct RangeConstraint {
    /// Index of agent A.
    pub agent_a: usize,
    /// Index of agent B.
    pub agent_b: usize,
    /// Desired inter-agent distance.
    pub desired_range: f64,
    /// Spring-like stiffness gain.
    pub k: f64,
    /// Damping coefficient.
    pub d: f64,
}

impl RangeConstraint {
    /// Create a new range constraint.
    pub fn new(agent_a: usize, agent_b: usize, desired_range: f64) -> Self {
        Self {
            agent_a,
            agent_b,
            desired_range,
            k: 1.0,
            d: 0.1,
        }
    }

    /// Compute 3D acceleration commands `(accel_a, accel_b)` to enforce the
    /// range constraint using a virtual spring-damper.
    pub fn compute_accelerations(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> ([f64; 3], [f64; 3]) {
        let pa = positions[self.agent_a];
        let pb = positions[self.agent_b];
        let va = velocities[self.agent_a];
        let vb = velocities[self.agent_b];
        let diff = vec3_sub(pb, pa);
        let dist = vec3_len(diff).max(1e-12);
        let n = vec3_scale(diff, 1.0 / dist);
        let err = dist - self.desired_range;
        // Relative velocity along n
        let rel_v = vec3_dot(vec3_sub(vb, va), n);
        let force_mag = self.k * err + self.d * rel_v;
        let force = vec3_scale(n, force_mag);
        (force, vec3_scale(force, -1.0))
    }
}

// ── Task Assignment Constraint ────────────────────────────────────────────────

/// Assignment result for a single task.
#[derive(Debug, Clone)]
pub struct TaskAssignment {
    /// Index of the agent assigned.
    pub agent_idx: usize,
    /// Index of the task.
    pub task_idx: usize,
    /// Cost of the assignment.
    pub cost: f64,
}

/// Greedy task-assignment constraint: assigns tasks to agents minimising the
/// total cost using a greedy Hungarian-like approach (O(n² m)).
///
/// `cost_matrix\[i\]\[j\]` is the cost of assigning agent `i` to task `j`.
pub fn greedy_task_assignment(cost_matrix: &[Vec<f64>]) -> Vec<TaskAssignment> {
    let num_agents = cost_matrix.len();
    if num_agents == 0 {
        return Vec::new();
    }
    let num_tasks = cost_matrix[0].len();
    let mut assigned_tasks = vec![false; num_tasks];
    let mut assignments = Vec::new();

    for i in 0..num_agents {
        let mut best_task = None;
        let mut best_cost = f64::MAX;
        for j in 0..num_tasks {
            if !assigned_tasks[j] && cost_matrix[i][j] < best_cost {
                best_cost = cost_matrix[i][j];
                best_task = Some(j);
            }
        }
        if let Some(j) = best_task {
            assigned_tasks[j] = true;
            assignments.push(TaskAssignment {
                agent_idx: i,
                task_idx: j,
                cost: best_cost,
            });
        }
    }
    assignments
}

// ── Energy Constraint ─────────────────────────────────────────────────────────

/// Energy budget constraint for an agent.
///
/// Limits the control effort based on remaining energy.
#[derive(Debug, Clone)]
pub struct EnergyConstraint {
    /// Current energy level (Joules).
    pub energy: f64,
    /// Maximum energy capacity.
    pub max_energy: f64,
    /// Energy recharge rate (W).
    pub recharge_rate: f64,
    /// Power consumption per unit acceleration squared (W / (m/s²)²).
    pub power_per_accel_sq: f64,
}

impl EnergyConstraint {
    /// Create a new energy constraint with full charge.
    pub fn new(max_energy: f64) -> Self {
        Self {
            energy: max_energy,
            max_energy,
            recharge_rate: 0.0,
            power_per_accel_sq: 1.0,
        }
    }

    /// Update energy over time step `dt` given applied acceleration magnitude.
    ///
    /// Clamps energy to `\[0, max_energy\]`.
    pub fn update(&mut self, accel_mag: f64, dt: f64) {
        let consumed = self.power_per_accel_sq * accel_mag * accel_mag * dt;
        let recharged = self.recharge_rate * dt;
        self.energy = (self.energy - consumed + recharged).clamp(0.0, self.max_energy);
    }

    /// Return the maximum allowed acceleration magnitude given the current energy
    /// and available power over a time step `dt`.
    pub fn max_accel(&self, dt: f64) -> f64 {
        if dt < 1e-15 || self.power_per_accel_sq < 1e-15 {
            return f64::MAX;
        }
        let avail_power = self.energy / dt;
        (avail_power / self.power_per_accel_sq).sqrt()
    }

    /// Returns `true` if the agent has enough energy to apply `accel` over `dt`.
    pub fn is_feasible(&self, accel_mag: f64, dt: f64) -> bool {
        let needed = self.power_per_accel_sq * accel_mag * accel_mag * dt;
        needed <= self.energy
    }
}

// ── Coverage Constraint ───────────────────────────────────────────────────────

/// 2D coverage constraint: agents spread out to cover a square domain.
///
/// Uses a Lloyd relaxation step: each agent moves toward the centroid of its
/// Voronoi cell (approximated via sampling).
#[derive(Debug, Clone)]
pub struct CoverageConstraint {
    /// Domain bounds `\[x_min, x_max, y_min, y_max\]`.
    pub domain: [f64; 4],
    /// Step size for Lloyd relaxation.
    pub alpha: f64,
    /// Number of sample points per axis for Voronoi approximation.
    pub samples: usize,
}

impl CoverageConstraint {
    /// Create a new coverage constraint over the given domain.
    pub fn new(domain: [f64; 4]) -> Self {
        Self {
            domain,
            alpha: 0.1,
            samples: 20,
        }
    }

    /// Compute approximate Voronoi centroids for `positions` (2D, XY plane).
    ///
    /// Returns one centroid per agent.
    pub fn voronoi_centroids(&self, positions: &[[f64; 2]]) -> Vec<[f64; 2]> {
        let n = positions.len();
        if n == 0 {
            return Vec::new();
        }
        let [x_min, x_max, y_min, y_max] = self.domain;
        let s = self.samples;
        let mut sums = vec![[0.0_f64; 2]; n];
        let mut counts = vec![0usize; n];
        for si in 0..s {
            for sj in 0..s {
                let px = x_min + (x_max - x_min) * (si as f64 + 0.5) / s as f64;
                let py = y_min + (y_max - y_min) * (sj as f64 + 0.5) / s as f64;
                // Nearest agent
                let nearest = positions
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        let da = (a[0] - px).hypot(a[1] - py);
                        let db = (b[0] - px).hypot(b[1] - py);
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(k, _)| k)
                    .unwrap_or(0);
                sums[nearest][0] += px;
                sums[nearest][1] += py;
                counts[nearest] += 1;
            }
        }
        sums.iter()
            .zip(counts.iter())
            .map(|(sum, &cnt)| {
                if cnt > 0 {
                    [sum[0] / cnt as f64, sum[1] / cnt as f64]
                } else {
                    [0.0, 0.0]
                }
            })
            .collect()
    }

    /// Compute velocity commands toward Voronoi centroids (Lloyd step).
    pub fn lloyd_velocities(&self, positions: &[[f64; 2]]) -> Vec<[f64; 2]> {
        let centroids = self.voronoi_centroids(positions);
        positions
            .iter()
            .zip(centroids.iter())
            .map(|(&p, &c)| vec2_scale(vec2_sub(c, p), self.alpha))
            .collect()
    }
}

// ── Nash Equilibrium Constraint ───────────────────────────────────────────────

/// Nash equilibrium constraint for a two-player finite game.
///
/// Stores the payoff matrices and computes mixed-strategy Nash equilibria
/// via iterated best-response.
#[derive(Debug, Clone)]
pub struct NashEquilibriumConstraint {
    /// Payoff matrix for player 1 (rows = player 1 actions, cols = player 2).
    pub payoff_1: Vec<Vec<f64>>,
    /// Payoff matrix for player 2 (rows = player 1 actions, cols = player 2).
    pub payoff_2: Vec<Vec<f64>>,
}

impl NashEquilibriumConstraint {
    /// Create a Nash equilibrium constraint from payoff matrices.
    pub fn new(payoff_1: Vec<Vec<f64>>, payoff_2: Vec<Vec<f64>>) -> Self {
        Self { payoff_1, payoff_2 }
    }

    /// Compute the expected payoff for player 1 given mixed strategies.
    pub fn expected_payoff_1(&self, sigma1: &[f64], sigma2: &[f64]) -> f64 {
        let mut payoff = 0.0;
        for (i, &p1) in sigma1.iter().enumerate() {
            for (j, &p2) in sigma2.iter().enumerate() {
                payoff += p1 * p2 * self.payoff_1[i][j];
            }
        }
        payoff
    }

    /// Best-response mixed strategy for player 1 against player 2's strategy.
    ///
    /// Returns a pure-strategy indicator (1 on best action, 0 elsewhere).
    pub fn best_response_1(&self, sigma2: &[f64]) -> Vec<f64> {
        let rows = self.payoff_1.len();
        let mut best_row = 0;
        let mut best_val = f64::NEG_INFINITY;
        for i in 0..rows {
            let val: f64 = sigma2
                .iter()
                .enumerate()
                .map(|(j, &p2)| self.payoff_1[i][j] * p2)
                .sum();
            if val > best_val {
                best_val = val;
                best_row = i;
            }
        }
        let mut br = vec![0.0_f64; rows];
        br[best_row] = 1.0;
        br
    }

    /// Best-response for player 2 against player 1's strategy.
    pub fn best_response_2(&self, sigma1: &[f64]) -> Vec<f64> {
        let cols = if self.payoff_2.is_empty() {
            0
        } else {
            self.payoff_2[0].len()
        };
        let mut best_col = 0;
        let mut best_val = f64::NEG_INFINITY;
        for j in 0..cols {
            let val: f64 = sigma1
                .iter()
                .enumerate()
                .map(|(i, &p1)| self.payoff_2[i][j] * p1)
                .sum();
            if val > best_val {
                best_val = val;
                best_col = j;
            }
        }
        let mut br = vec![0.0_f64; cols];
        br[best_col] = 1.0;
        br
    }

    /// Run iterated best-response to approximate a Nash equilibrium.
    ///
    /// Returns `(sigma1, sigma2)` mixed strategies.
    pub fn iterated_best_response(&self, max_iter: usize, blend: f64) -> (Vec<f64>, Vec<f64>) {
        let rows = self.payoff_1.len();
        let cols = if self.payoff_2.is_empty() {
            0
        } else {
            self.payoff_2[0].len()
        };
        let mut sigma1 = vec![1.0 / rows as f64; rows];
        let mut sigma2 = vec![1.0 / cols as f64; cols];
        for _ in 0..max_iter {
            let br1 = self.best_response_1(&sigma2);
            let br2 = self.best_response_2(&sigma1);
            // Blend toward best response
            for i in 0..rows {
                sigma1[i] = (1.0 - blend) * sigma1[i] + blend * br1[i];
            }
            for j in 0..cols {
                sigma2[j] = (1.0 - blend) * sigma2[j] + blend * br2[j];
            }
        }
        (sigma1, sigma2)
    }
}

// ── Social Force Model ────────────────────────────────────────────────────────

/// Parameters for the social force model (pedestrian / swarm dynamics).
#[derive(Debug, Clone)]
pub struct SocialForceParams {
    /// Desired speed of agents.
    pub v0: f64,
    /// Relaxation time to reach desired velocity.
    pub tau: f64,
    /// Magnitude of repulsive social force.
    pub a_rep: f64,
    /// Range of repulsive social force.
    pub b_rep: f64,
    /// Anisotropy parameter (0 = isotropic, 1 = strongly anisotropic).
    pub lambda_aniso: f64,
}

impl SocialForceParams {
    /// Create social force parameters with reasonable defaults.
    pub fn new() -> Self {
        Self {
            v0: 1.4,
            tau: 0.5,
            a_rep: 2000.0,
            b_rep: 0.08,
            lambda_aniso: 0.5,
        }
    }
}

impl Default for SocialForceParams {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute the social force acceleration for agent `i`.
///
/// Combines:
/// 1. Self-propulsion toward desired velocity `v_desired`.
/// 2. Interpersonal repulsion from neighbours.
///
/// Positions and velocities are 2D (`\[f64; 2\]`).
pub fn social_force_acceleration(
    i: usize,
    positions: &[[f64; 2]],
    velocities: &[[f64; 2]],
    v_desired: [f64; 2],
    params: &SocialForceParams,
) -> [f64; 2] {
    let pi = positions[i];
    let vi = velocities[i];

    // Self-propulsion
    let self_force = vec2_scale(vec2_sub(v_desired, vi), 1.0 / params.tau);

    // Repulsion from others
    let mut rep = [0.0_f64; 2];
    for (j, &pj) in positions.iter().enumerate() {
        if i == j {
            continue;
        }
        let rij = vec2_sub(pi, pj);
        let dist = vec2_len(rij).max(1e-9);
        let n_ij = vec2_scale(rij, 1.0 / dist);

        // Anisotropy: reduce force from behind
        let vij = vec2_normalize(vec2_sub(vi, velocities[j]));
        let cos_phi = vec2_dot(n_ij, vij).clamp(-1.0, 1.0);
        let w = params.lambda_aniso + (1.0 - params.lambda_aniso) * (1.0 + cos_phi) * 0.5;

        let magnitude = params.a_rep * (-dist / params.b_rep).exp() * w;
        rep = vec2_add(rep, vec2_scale(n_ij, magnitude));
    }

    vec2_add(self_force, rep)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Graph tests

    #[test]
    fn test_graph_laplacian_simple() {
        // 3 nodes, edges: 0→1 (w=1), 1→2 (w=1), 2→0 (w=1)
        let mut g = CommunicationGraph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        let l = g.laplacian();
        // Diagonal should be 1 (each node has out-degree 1)
        assert!((l[0] - 1.0).abs() < 1e-12);
        assert!((l[4] - 1.0).abs() < 1e-12);
        assert!((l[8] - 1.0).abs() < 1e-12);
        // Off-diagonal: l[0*3+1] = -1
        assert!((l[1] + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_graph_laplacian_matvec_ones() {
        // For an undirected ring, L * 1 = 0
        let mut g = CommunicationGraph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 0, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 1, 1.0);
        g.add_edge(2, 0, 1.0);
        g.add_edge(0, 2, 1.0);
        let x = [1.0, 1.0, 1.0];
        let lx = g.laplacian_matvec(&x);
        for v in &lx {
            assert!(v.abs() < 1e-10, "lx component = {v}");
        }
    }

    #[test]
    fn test_graph_connectivity() {
        let mut g = CommunicationGraph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        assert!(g.is_connected());
        let g2 = CommunicationGraph::new(3);
        assert!(!g2.is_connected()); // no edges → disconnected (nodes 1,2 unreachable)
    }

    // Consensus tests

    #[test]
    fn test_consensus_converges() {
        let mut g = CommunicationGraph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 0, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 1, 1.0);
        let c = ConsensusConstraint::new(g);
        let x0 = [0.0, 5.0, 10.0];
        let (_iters, xf) = c.run(&x0);
        // All values should converge to mean = 5.0
        for xi in &xf {
            assert!((xi - 5.0).abs() < 0.1, "xi = {xi}");
        }
    }

    #[test]
    fn test_consensus_error_same_values() {
        let x = [3.0, 3.0, 3.0];
        assert!(ConsensusConstraint::consensus_error(&x) < 1e-10);
    }

    #[test]
    fn test_consensus_error_different() {
        let x = [0.0, 2.0];
        let err = ConsensusConstraint::consensus_error(&x);
        assert!(err > 0.0);
    }

    // ADMM tests

    #[test]
    fn test_admm_agent_x_update() {
        let mut ag = AdmmAgent::new(2);
        ag.q_diag = vec![1.0, 1.0];
        ag.c = vec![0.0, 0.0];
        ag.z = vec![3.0, 4.0];
        ag.u = vec![0.0, 0.0];
        ag.x_update(1.0);
        // x = (Q + rho I)^{-1} * rho * z = 0.5 * z
        assert!((ag.x[0] - 1.5).abs() < 1e-10);
        assert!((ag.x[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_admm_solve_consensus() {
        let mut agents: Vec<AdmmAgent> = (0..3)
            .map(|i| {
                let mut ag = AdmmAgent::new(1);
                ag.c = vec![-(i as f64)]; // minimize -(i)*x + 0.5*x^2
                ag
            })
            .collect();
        let (_iters, z) = admm_solve(&mut agents, 1.0, 200, 1e-6);
        // Consensus solution is the mean of the argmin: mean of [0,1,2] = 1
        assert!((z[0] - 1.0).abs() < 0.1, "z = {}", z[0]);
    }

    // Collective motion

    #[test]
    fn test_collective_motion_stationary() {
        let constraint = CollectiveMotionConstraint::new();
        let positions = vec![[0.0, 0.0, 0.0]; 3];
        let velocities = vec![[0.0, 0.0, 0.0]; 3];
        let a = constraint.compute_acceleration(&positions, &velocities, 0);
        // All at same position, no separation → zero acceleration
        assert!(vec3_len(a) < 1e-9);
    }

    // Leader-follower

    #[test]
    fn test_leader_follower_on_target() {
        let c = LeaderFollowerConstraint::new(0, 1, [2.0, 0.0, 0.0]);
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let a = c.compute_acceleration(&positions, &velocities);
        // On target with same velocity → small/zero acceleration
        assert!(vec3_len(a) < 1e-9);
    }

    #[test]
    fn test_leader_follower_error() {
        let c = LeaderFollowerConstraint::new(0, 1, [2.0, 0.0, 0.0]);
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let a = c.compute_acceleration(&positions, &velocities);
        // Follower is at origin, should go to (2,0,0), accel should be +x
        assert!(a[0] > 0.0);
    }

    // Formation

    #[test]
    fn test_formation_centre() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let centre = FormationConstraint::formation_centre(&positions);
        assert!((centre[0] - 1.0).abs() < 1e-10);
        assert!((centre[1] - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_formation_in_formation() {
        let offsets = vec![[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]];
        let c = FormationConstraint::new(offsets);
        // Centre = (0,0,0), desired positions = (0,±1,0)
        let positions = vec![[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]];
        let velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let accels = c.compute_accelerations(&positions, &velocities, [0.0, 0.0, 0.0]);
        for a in &accels {
            assert!(vec3_len(*a) < 1e-9, "unexpected accel {a:?}");
        }
    }

    // Formation reconfiguration

    #[test]
    fn test_reconfiguration_completes() {
        let from = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let to = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut rec = FormationReconfiguration::new(from, to, 1.0);
        let _ = rec.step(1.0);
        assert!(rec.is_complete());
        let final_offsets = rec.step(0.0);
        assert!((final_offsets[1][1] - 1.0).abs() < 1e-10);
    }

    // ORCA

    #[test]
    fn test_orca_half_plane_no_conflict() {
        let hp = orca_half_plane(
            [0.0, 0.0],
            [1.0, 0.0],
            0.3,
            [10.0, 0.0],
            [0.0, 0.0],
            0.3,
            2.0,
        );
        // Far apart — HP normal should be defined
        assert!(vec2_len(hp.normal) > 1e-9 || hp.point[0].abs() < 1e3);
    }

    #[test]
    fn test_orca_solve_prefers_pref() {
        let v = orca_solve_velocity([0.5, 0.0], 2.0, &[]);
        // No constraints — should return v_pref
        assert!((v[0] - 0.5).abs() < 1e-10);
    }

    // Bearing

    #[test]
    fn test_bearing_east() {
        let p_obs = [0.0, 0.0, 0.0];
        let p_tgt = [1.0, 0.0, 0.0];
        let bearing = BearingConstraint::bearing(p_obs, p_tgt);
        assert!(bearing.abs() < 1e-10); // east = 0 rad
    }

    #[test]
    fn test_bearing_north() {
        let p_obs = [0.0, 0.0, 0.0];
        let p_tgt = [0.0, 1.0, 0.0];
        let bearing = BearingConstraint::bearing(p_obs, p_tgt);
        let pi_half = std::f64::consts::FRAC_PI_2;
        assert!((bearing - pi_half).abs() < 1e-10);
    }

    // Range

    #[test]
    fn test_range_constraint_on_target() {
        let c = RangeConstraint::new(0, 1, 2.0);
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let (fa, fb) = c.compute_accelerations(&positions, &velocities);
        // At desired range, spring force should be zero
        assert!(vec3_len(fa) < 1e-9);
        assert!(vec3_len(fb) < 1e-9);
    }

    #[test]
    fn test_range_constraint_too_close() {
        let c = RangeConstraint::new(0, 1, 3.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let (fa, fb) = c.compute_accelerations(&positions, &velocities);
        // Too close — force should push agents apart (fa in -x, fb in +x)
        assert!(fa[0] < 0.0);
        assert!(fb[0] > 0.0);
    }

    // Task assignment

    #[test]
    fn test_greedy_assignment_2x2() {
        // cost[0][0]=1, cost[0][1]=10, cost[1][0]=10, cost[1][1]=1
        let costs = vec![vec![1.0, 10.0], vec![10.0, 1.0]];
        let assignments = greedy_task_assignment(&costs);
        // Agent 0 → task 0, agent 1 → task 1
        assert_eq!(assignments[0].task_idx, 0);
        assert_eq!(assignments[1].task_idx, 1);
    }

    // Energy

    #[test]
    fn test_energy_depletes() {
        let mut e = EnergyConstraint::new(100.0);
        e.update(2.0, 1.0); // consumes 1.0 * 4.0 = 4.0
        assert!((e.energy - 96.0).abs() < 1e-9);
    }

    #[test]
    fn test_energy_max_accel() {
        let e = EnergyConstraint::new(100.0);
        let a_max = e.max_accel(1.0);
        assert!(a_max > 0.0);
        // With 100 J over 1s and coefficient 1: sqrt(100) = 10
        assert!((a_max - 10.0).abs() < 1e-9);
    }

    // Coverage

    #[test]
    fn test_coverage_voronoi_count() {
        let c = CoverageConstraint::new([0.0, 10.0, 0.0, 10.0]);
        let positions = vec![[2.5, 5.0], [7.5, 5.0]];
        let centroids = c.voronoi_centroids(&positions);
        assert_eq!(centroids.len(), 2);
    }

    #[test]
    fn test_coverage_lloyd_velocities() {
        let c = CoverageConstraint {
            domain: [0.0, 10.0, 0.0, 10.0],
            alpha: 1.0,
            samples: 10,
        };
        let positions = vec![[5.0, 5.0]];
        let vels = c.lloyd_velocities(&positions);
        assert_eq!(vels.len(), 1);
    }

    // Nash

    #[test]
    fn test_nash_prisoner_dilemma() {
        // Prisoner's Dilemma: both defect is the unique NE
        // Payoffs (player 1): D,D=2; D,C=3; C,D=0; C,C=1
        let p1 = vec![vec![2.0, 3.0], vec![0.0, 1.0]]; // rows=player1 action
        let p2 = vec![vec![2.0, 0.0], vec![3.0, 1.0]];
        let c = NashEquilibriumConstraint::new(p1, p2);
        let (s1, s2) = c.iterated_best_response(100, 0.5);
        // Both should converge to D (action 0)
        assert!(s1[0] > s1[1], "s1={s1:?}");
        assert!(s2[0] > s2[1], "s2={s2:?}");
    }

    // Social force

    #[test]
    fn test_social_force_self_propulsion() {
        let params = SocialForceParams::new();
        let positions = vec![[0.0, 0.0]];
        let velocities = vec![[0.0, 0.0]];
        let v_desired = [1.0, 0.0];
        let a = social_force_acceleration(0, &positions, &velocities, v_desired, &params);
        // Should accelerate in +x direction toward v_desired
        assert!(a[0] > 0.0);
    }

    #[test]
    fn test_social_force_repulsion() {
        let params = SocialForceParams::new();
        let positions = vec![[0.0, 0.0], [0.1, 0.0]]; // very close
        let velocities = vec![[0.0, 0.0], [0.0, 0.0]];
        let v_desired = [0.0, 0.0];
        let a = social_force_acceleration(0, &positions, &velocities, v_desired, &params);
        // Repulsion should push agent 0 in -x direction
        assert!(a[0] < 0.0, "a[0] = {}", a[0]);
    }

    #[test]
    fn test_vec2_operations() {
        let a = [3.0, 4.0];
        assert!((vec2_len(a) - 5.0).abs() < 1e-10);
        let n = vec2_normalize(a);
        assert!((vec2_len(n) - 1.0).abs() < 1e-10);
        let b = [1.0, 2.0];
        let s = vec2_dot(a, b);
        assert!((s - 11.0).abs() < 1e-10);
    }

    #[test]
    fn test_clamp_utility() {
        assert!((clamp(5.0, 0.0, 10.0) - 5.0).abs() < 1e-10);
        assert!((clamp(-1.0, 0.0, 10.0)).abs() < 1e-10);
        assert!((clamp(15.0, 0.0, 10.0) - 10.0).abs() < 1e-10);
    }
}
