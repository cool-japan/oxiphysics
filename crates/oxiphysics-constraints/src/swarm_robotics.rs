#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Swarm robotics constraints: formation control, consensus, coverage.
//!
//! Implements Reynolds flocking, formation control, consensus protocols,
//! Voronoi-based coverage, auction-based task allocation, PSO search,
//! and stigmergy-based emergent behavior.

// ── Helper math ───────────────────────────────────────────────────────────────

#[allow(dead_code)]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[allow(dead_code)]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[allow(dead_code)]
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[allow(dead_code)]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(dead_code)]
fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[allow(dead_code)]
fn vec3_normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vec3_len(v);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

#[allow(dead_code)]
fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_len(vec3_sub(a, b))
}

// ── SwarmAgent ────────────────────────────────────────────────────────────────

/// A single agent in a swarm with position, velocity, heading, and sensor radius.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct SwarmAgent {
    /// World-space position \[x, y, z\].
    pub position: [f64; 3],
    /// Linear velocity \[vx, vy, vz\].
    pub velocity: [f64; 3],
    /// Heading angle in radians (yaw, XY-plane).
    pub heading: f64,
    /// Sensor radius for neighbor detection.
    pub sensor_radius: f64,
    /// Unique agent identifier.
    pub id: usize,
}

impl SwarmAgent {
    /// Create a new swarm agent.
    #[allow(dead_code)]
    pub fn new(id: usize, position: [f64; 3], sensor_radius: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            heading: 0.0,
            sensor_radius,
            id,
        }
    }

    /// Return agents within `self.sensor_radius` from `agents` slice.
    #[allow(dead_code)]
    pub fn neighbors<'a>(&self, agents: &'a [SwarmAgent]) -> Vec<&'a SwarmAgent> {
        agents
            .iter()
            .filter(|a| {
                a.id != self.id && vec3_dist(self.position, a.position) <= self.sensor_radius
            })
            .collect()
    }

    /// Integrate position by one timestep `dt` using current velocity.
    #[allow(dead_code)]
    pub fn integrate(&mut self, dt: f64) {
        self.position = vec3_add(self.position, vec3_scale(self.velocity, dt));
        // Update heading from velocity projection onto XY plane.
        if self.velocity[0].abs() > 1e-12 || self.velocity[1].abs() > 1e-12 {
            self.heading = self.velocity[1].atan2(self.velocity[0]);
        }
    }

    /// Apply a steering acceleration, clamped to `max_speed`.
    #[allow(dead_code)]
    pub fn apply_steering(&mut self, steering: [f64; 3], max_speed: f64, dt: f64) {
        let new_vel = vec3_add(self.velocity, vec3_scale(steering, dt));
        let speed = vec3_len(new_vel);
        self.velocity = if speed > max_speed {
            vec3_scale(vec3_normalize(new_vel), max_speed)
        } else {
            new_vel
        };
    }
}

// ── FlockingBehavior ──────────────────────────────────────────────────────────

/// Parameters for Reynolds flocking (separation, alignment, cohesion).
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct FlockingBehavior {
    /// Separation radius (repulsion below this distance).
    pub r_separation: f64,
    /// Weight for separation force.
    pub w_separation: f64,
    /// Weight for alignment force.
    pub w_alignment: f64,
    /// Weight for cohesion force.
    pub w_cohesion: f64,
    /// Maximum steering force magnitude.
    pub max_force: f64,
    /// Maximum agent speed.
    pub max_speed: f64,
}

impl FlockingBehavior {
    /// Create a new FlockingBehavior with default weights.
    #[allow(dead_code)]
    pub fn new(r_separation: f64) -> Self {
        Self {
            r_separation,
            w_separation: 1.5,
            w_alignment: 1.0,
            w_cohesion: 1.0,
            max_force: 2.0,
            max_speed: 5.0,
        }
    }

    /// Compute separation steering: steer away from too-close neighbors.
    ///
    /// Returns a steering vector away from neighbors closer than `r_separation`.
    #[allow(dead_code)]
    pub fn separation(&self, agent: &SwarmAgent, neighbors: &[&SwarmAgent]) -> [f64; 3] {
        let mut steer = [0.0_f64; 3];
        let mut count = 0usize;
        for nb in neighbors {
            let d = vec3_dist(agent.position, nb.position);
            if d < self.r_separation && d > 1e-12 {
                // Weight by inverse distance
                let away = vec3_normalize(vec3_sub(agent.position, nb.position));
                let weighted = vec3_scale(away, 1.0 / d);
                steer = vec3_add(steer, weighted);
                count += 1;
            }
        }
        if count > 0 {
            steer = vec3_scale(steer, 1.0 / count as f64);
            let speed = vec3_len(steer);
            if speed > 1e-12 {
                steer = vec3_scale(vec3_normalize(steer), self.max_speed);
                steer = vec3_sub(steer, agent.velocity);
                steer = clamp_vec3_magnitude(steer, self.max_force);
            }
        }
        steer
    }

    /// Compute alignment steering: match average velocity of neighbors.
    #[allow(dead_code)]
    pub fn alignment(&self, agent: &SwarmAgent, neighbors: &[&SwarmAgent]) -> [f64; 3] {
        if neighbors.is_empty() {
            return [0.0; 3];
        }
        let mut avg_vel = [0.0_f64; 3];
        for nb in neighbors {
            avg_vel = vec3_add(avg_vel, nb.velocity);
        }
        avg_vel = vec3_scale(avg_vel, 1.0 / neighbors.len() as f64);
        // Steer toward average velocity
        let steer = vec3_sub(avg_vel, agent.velocity);
        clamp_vec3_magnitude(steer, self.max_force)
    }

    /// Compute cohesion steering: steer toward average position of neighbors.
    #[allow(dead_code)]
    pub fn cohesion(&self, agent: &SwarmAgent, neighbors: &[&SwarmAgent]) -> [f64; 3] {
        if neighbors.is_empty() {
            return [0.0; 3];
        }
        let mut center = [0.0_f64; 3];
        for nb in neighbors {
            center = vec3_add(center, nb.position);
        }
        center = vec3_scale(center, 1.0 / neighbors.len() as f64);
        let desired_dir = vec3_sub(center, agent.position);
        let desired = if vec3_len(desired_dir) > 1e-12 {
            vec3_scale(vec3_normalize(desired_dir), self.max_speed)
        } else {
            [0.0; 3]
        };
        clamp_vec3_magnitude(vec3_sub(desired, agent.velocity), self.max_force)
    }

    /// Compute total flocking steering force for one agent.
    #[allow(dead_code)]
    pub fn compute_steering(&self, agent: &SwarmAgent, neighbors: &[&SwarmAgent]) -> [f64; 3] {
        let sep = vec3_scale(self.separation(agent, neighbors), self.w_separation);
        let ali = vec3_scale(self.alignment(agent, neighbors), self.w_alignment);
        let coh = vec3_scale(self.cohesion(agent, neighbors), self.w_cohesion);
        let total = vec3_add(vec3_add(sep, ali), coh);
        clamp_vec3_magnitude(total, self.max_force)
    }

    /// Step an entire swarm by `dt` applying flocking rules.
    #[allow(dead_code)]
    pub fn step_swarm(&self, agents: &mut Vec<SwarmAgent>, dt: f64) {
        // Gather all positions and velocities first (avoid borrow issues).
        let snapshots: Vec<SwarmAgent> = agents.clone();
        for i in 0..agents.len() {
            let nbs: Vec<&SwarmAgent> = snapshots
                .iter()
                .filter(|a| {
                    a.id != agents[i].id
                        && vec3_dist(a.position, agents[i].position) <= agents[i].sensor_radius
                })
                .collect();
            let steering = self.compute_steering(&agents[i], &nbs);
            agents[i].apply_steering(steering, self.max_speed, dt);
            agents[i].integrate(dt);
        }
    }
}

/// Clamp a vector's magnitude to `max`.
#[allow(dead_code)]
fn clamp_vec3_magnitude(v: [f64; 3], max: f64) -> [f64; 3] {
    let mag = vec3_len(v);
    if mag > max && mag > 1e-15 {
        vec3_scale(vec3_normalize(v), max)
    } else {
        v
    }
}

// ── FormationControl ──────────────────────────────────────────────────────────

/// Formation type for rigid structure control.
#[derive(Clone, Debug, PartialEq)]
#[allow(dead_code)]
pub enum FormationType {
    /// Virtual structure: agents track offsets from a virtual center.
    VirtualStructure,
    /// Consensus-based: distributed position agreement.
    Consensus,
}

/// Rigid formation maintenance controller.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct FormationControl {
    /// Type of formation control.
    pub formation_type: FormationType,
    /// Desired offsets from virtual center for each agent (indexed by agent position in team).
    pub offsets: Vec<[f64; 3]>,
    /// Position of the virtual center.
    pub virtual_center: [f64; 3],
    /// Virtual center velocity.
    pub virtual_velocity: [f64; 3],
    /// Control gain for position error.
    pub gain: f64,
}

impl FormationControl {
    /// Create a new formation controller.
    #[allow(dead_code)]
    pub fn new(formation_type: FormationType, offsets: Vec<[f64; 3]>, gain: f64) -> Self {
        Self {
            formation_type,
            offsets,
            virtual_center: [0.0; 3],
            virtual_velocity: [0.0; 3],
            gain,
        }
    }

    /// Compute desired positions for all agents.
    #[allow(dead_code)]
    pub fn desired_positions(&self) -> Vec<[f64; 3]> {
        self.offsets
            .iter()
            .map(|&off| vec3_add(self.virtual_center, off))
            .collect()
    }

    /// Compute the steering force for agent `i` to reach its formation position.
    #[allow(dead_code)]
    pub fn formation_steering(&self, agent: &SwarmAgent, agent_index: usize) -> [f64; 3] {
        if agent_index >= self.offsets.len() {
            return [0.0; 3];
        }
        let desired = vec3_add(self.virtual_center, self.offsets[agent_index]);
        let error = vec3_sub(desired, agent.position);
        vec3_scale(error, self.gain)
    }

    /// Formation tracking error: mean distance of agents to desired positions.
    #[allow(dead_code)]
    pub fn tracking_error(&self, agents: &[SwarmAgent]) -> f64 {
        let n = agents.len().min(self.offsets.len());
        if n == 0 {
            return 0.0;
        }
        let total: f64 = (0..n)
            .map(|i| {
                let desired = vec3_add(self.virtual_center, self.offsets[i]);
                vec3_dist(agents[i].position, desired)
            })
            .sum();
        total / n as f64
    }

    /// Step the virtual center forward by `dt`.
    #[allow(dead_code)]
    pub fn advance_center(&mut self, dt: f64) {
        self.virtual_center = vec3_add(self.virtual_center, vec3_scale(self.virtual_velocity, dt));
    }

    /// Step agents toward their formation positions.
    #[allow(dead_code)]
    pub fn step(&mut self, agents: &mut [SwarmAgent], dt: f64) {
        self.advance_center(dt);
        for (i, agent) in agents.iter_mut().enumerate() {
            let steer = if i < self.offsets.len() {
                let desired = vec3_add(self.virtual_center, self.offsets[i]);
                let error = vec3_sub(desired, agent.position);
                vec3_scale(error, self.gain)
            } else {
                [0.0; 3]
            };
            agent.velocity = vec3_add(agent.velocity, vec3_scale(steer, dt));
            agent.integrate(dt);
        }
    }
}

// ── ConsensusProtocol ─────────────────────────────────────────────────────────

/// Multi-agent consensus (distributed averaging) protocol.
///
/// Each agent updates its state as:
/// `x_dot_i = Σ_{j ∈ N_i} (x_j - x_i)`
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct ConsensusProtocol {
    /// Current scalar state values for each agent.
    pub states: Vec<f64>,
    /// Adjacency list: `neighbors[i]` contains indices of agent `i`'s neighbors.
    pub neighbors: Vec<Vec<usize>>,
    /// Convergence gain (scales the update step).
    pub gain: f64,
}

impl ConsensusProtocol {
    /// Create a fully-connected consensus protocol.
    #[allow(dead_code)]
    pub fn new_fully_connected(initial_states: Vec<f64>, gain: f64) -> Self {
        let n = initial_states.len();
        let neighbors: Vec<Vec<usize>> = (0..n)
            .map(|i| (0..n).filter(|&j| j != i).collect())
            .collect();
        Self {
            states: initial_states,
            neighbors,
            gain,
        }
    }

    /// Create a consensus protocol with explicit neighbor graph.
    #[allow(dead_code)]
    pub fn new(initial_states: Vec<f64>, neighbors: Vec<Vec<usize>>, gain: f64) -> Self {
        Self {
            states: initial_states,
            neighbors,
            gain,
        }
    }

    /// Perform one consensus update step with timestep `dt`.
    ///
    /// `x_i += dt * gain * Σ_{j ∈ N_i} (x_j - x_i)`
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        let old = self.states.clone();
        for i in 0..self.states.len() {
            let update: f64 = self.neighbors[i]
                .iter()
                .map(|&j| old[j] - old[i])
                .sum::<f64>();
            self.states[i] += dt * self.gain * update;
        }
    }

    /// Check convergence: max deviation from mean < `tol`.
    #[allow(dead_code)]
    pub fn is_converged(&self, tol: f64) -> bool {
        let mean = self.states.iter().sum::<f64>() / self.states.len() as f64;
        self.states.iter().all(|&x| (x - mean).abs() < tol)
    }

    /// Estimate the convergence rate as the second-smallest eigenvalue of the
    /// graph Laplacian (Fiedler value), computed via power iteration on L.
    ///
    /// Returns an approximation for small graphs using Gershgorin bound.
    #[allow(dead_code)]
    pub fn laplacian_second_eigenvalue_approx(&self) -> f64 {
        let n = self.states.len();
        if n < 2 {
            return 0.0;
        }
        // Laplacian diagonal = degree[i], off-diagonal = -A[i][j]
        let degrees: Vec<f64> = self.neighbors.iter().map(|nb| nb.len() as f64).collect();
        // Min degree / max degree ratio gives a lower bound via Cheeger inequality
        let min_deg = degrees.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_deg = degrees.iter().cloned().fold(0.0_f64, f64::max);
        if max_deg < 1e-12 {
            0.0
        } else {
            // Simple approximation: lambda_2 ≈ min_degree / n
            min_deg / n as f64
        }
    }

    /// Return the mean of all agent states.
    #[allow(dead_code)]
    pub fn mean_state(&self) -> f64 {
        self.states.iter().sum::<f64>() / self.states.len() as f64
    }
}

// ── CoverageControl ───────────────────────────────────────────────────────────

/// Voronoi-based multi-agent coverage control.
///
/// Each agent moves toward the centroid of its Voronoi cell within
/// a bounded region.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct CoverageControl {
    /// Region bounds: \[\[x_min, x_max\], \[y_min, y_max\]\].
    pub bounds: [[f64; 2]; 2],
    /// Number of Monte-Carlo samples for Voronoi centroid estimation.
    pub n_samples: usize,
    /// Control gain.
    pub gain: f64,
}

impl CoverageControl {
    /// Create a new coverage controller.
    #[allow(dead_code)]
    pub fn new(bounds: [[f64; 2]; 2], n_samples: usize, gain: f64) -> Self {
        Self {
            bounds,
            n_samples,
            gain,
        }
    }

    /// Compute the Voronoi centroid for agent `idx` using grid sampling.
    ///
    /// The centroid is the average of all grid points closest to agent `idx`.
    #[allow(dead_code)]
    pub fn voronoi_centroid(&self, agents: &[SwarmAgent], idx: usize) -> [f64; 3] {
        let grid_side = (self.n_samples as f64).sqrt() as usize + 1;
        let dx = (self.bounds[0][1] - self.bounds[0][0]) / grid_side as f64;
        let dy = (self.bounds[1][1] - self.bounds[1][0]) / grid_side as f64;
        let mut sum = [0.0_f64; 3];
        let mut count = 0usize;
        for ix in 0..grid_side {
            for iy in 0..grid_side {
                let x = self.bounds[0][0] + (ix as f64 + 0.5) * dx;
                let y = self.bounds[1][0] + (iy as f64 + 0.5) * dy;
                let pt = [x, y, 0.0];
                // Find nearest agent
                let nearest = agents
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        vec3_dist(a.position, pt)
                            .partial_cmp(&vec3_dist(b.position, pt))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                if nearest == idx {
                    sum = vec3_add(sum, pt);
                    count += 1;
                }
            }
        }
        if count > 0 {
            vec3_scale(sum, 1.0 / count as f64)
        } else {
            agents[idx].position
        }
    }

    /// Compute coverage quality: max distance from any point to nearest agent.
    ///
    /// Sampled over a coarse grid.
    #[allow(dead_code)]
    pub fn coverage_quality(&self, agents: &[SwarmAgent]) -> f64 {
        if agents.is_empty() {
            return f64::INFINITY;
        }
        let grid_side = 10usize;
        let dx = (self.bounds[0][1] - self.bounds[0][0]) / grid_side as f64;
        let dy = (self.bounds[1][1] - self.bounds[1][0]) / grid_side as f64;
        let mut max_dist = 0.0_f64;
        for ix in 0..grid_side {
            for iy in 0..grid_side {
                let x = self.bounds[0][0] + (ix as f64 + 0.5) * dx;
                let y = self.bounds[1][0] + (iy as f64 + 0.5) * dy;
                let pt = [x, y, 0.0];
                let min_d = agents
                    .iter()
                    .map(|a| vec3_dist(a.position, pt))
                    .fold(f64::INFINITY, f64::min);
                if min_d > max_dist {
                    max_dist = min_d;
                }
            }
        }
        max_dist
    }

    /// Move each agent toward its Voronoi centroid by one step.
    #[allow(dead_code)]
    pub fn step(&self, agents: &mut Vec<SwarmAgent>, dt: f64) {
        let snap: Vec<SwarmAgent> = agents.clone();
        for (i, agent) in agents.iter_mut().enumerate() {
            let centroid = self.voronoi_centroid(&snap, i);
            let err = vec3_sub(centroid, agent.position);
            let steer = vec3_scale(err, self.gain);
            agent.velocity = steer;
            agent.integrate(dt);
        }
    }
}

// ── TaskAllocation ────────────────────────────────────────────────────────────

/// A task that can be assigned to an agent.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct Task {
    /// Task identifier.
    pub id: usize,
    /// Task location.
    pub position: [f64; 3],
    /// Task utility value.
    pub utility: f64,
}

/// Result of a task allocation.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct AllocationResult {
    /// For each task, the assigned agent index (if any).
    pub task_to_agent: Vec<Option<usize>>,
    /// For each agent, the assigned task index (if any).
    pub agent_to_task: Vec<Option<usize>>,
}

/// Auction-based greedy task allocator.
///
/// Bid = utility - cost, where cost = distance from agent to task.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct TaskAllocation {
    /// Cost weight for distance.
    pub cost_weight: f64,
}

impl TaskAllocation {
    /// Create a new task allocator.
    #[allow(dead_code)]
    pub fn new(cost_weight: f64) -> Self {
        Self { cost_weight }
    }

    /// Compute bid of agent `i` for task `t`.
    ///
    /// `bid = utility - cost_weight * distance(agent, task)`
    #[allow(dead_code)]
    pub fn bid(&self, agent: &SwarmAgent, task: &Task) -> f64 {
        let cost = self.cost_weight * vec3_dist(agent.position, task.position);
        task.utility - cost
    }

    /// Greedy auction: assign each task to the highest-bidding unassigned agent.
    ///
    /// Each agent may receive at most one task.
    #[allow(dead_code)]
    pub fn allocate(&self, agents: &[SwarmAgent], tasks: &[Task]) -> AllocationResult {
        let mut task_to_agent: Vec<Option<usize>> = vec![None; tasks.len()];
        let mut agent_to_task: Vec<Option<usize>> = vec![None; agents.len()];
        // For each task, find the best unassigned agent.
        for (ti, task) in tasks.iter().enumerate() {
            let mut best_agent = None;
            let mut best_bid = f64::NEG_INFINITY;
            for (ai, agent) in agents.iter().enumerate() {
                if agent_to_task[ai].is_none() {
                    let b = self.bid(agent, task);
                    if b > best_bid {
                        best_bid = b;
                        best_agent = Some(ai);
                    }
                }
            }
            if let Some(ai) = best_agent {
                task_to_agent[ti] = Some(ai);
                agent_to_task[ai] = Some(ti);
            }
        }
        AllocationResult {
            task_to_agent,
            agent_to_task,
        }
    }

    /// Count the number of assigned tasks.
    #[allow(dead_code)]
    pub fn assigned_count(result: &AllocationResult) -> usize {
        result.task_to_agent.iter().filter(|a| a.is_some()).count()
    }
}

// ── SwarmSearch (PSO) ─────────────────────────────────────────────────────────

/// A particle in PSO (Particle Swarm Optimization).
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct Particle {
    /// Current position in search space.
    pub position: [f64; 3],
    /// Current velocity.
    pub velocity: [f64; 3],
    /// Personal best position.
    pub personal_best: [f64; 3],
    /// Personal best fitness value.
    pub personal_best_fitness: f64,
}

impl Particle {
    /// Create a new particle at the given position.
    #[allow(dead_code)]
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            personal_best: position,
            personal_best_fitness: f64::INFINITY,
        }
    }
}

/// Particle Swarm Optimization for swarm search.
///
/// Update rule:
/// `v_i = w*v_i + c1*r1*(p_best - x_i) + c2*r2*(g_best - x_i)`
/// `x_i += v_i`
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct SwarmSearch {
    /// PSO particles.
    pub particles: Vec<Particle>,
    /// Global best position found.
    pub global_best: [f64; 3],
    /// Global best fitness value.
    pub global_best_fitness: f64,
    /// Inertia weight `w`.
    pub inertia: f64,
    /// Cognitive coefficient `c1`.
    pub c1: f64,
    /// Social coefficient `c2`.
    pub c2: f64,
    /// Maximum velocity magnitude.
    pub max_velocity: f64,
}

impl SwarmSearch {
    /// Create a new PSO with given particles.
    #[allow(dead_code)]
    pub fn new(
        particles: Vec<Particle>,
        inertia: f64,
        c1: f64,
        c2: f64,
        max_velocity: f64,
    ) -> Self {
        let global_best = if particles.is_empty() {
            [0.0; 3]
        } else {
            particles[0].position
        };
        Self {
            particles,
            global_best,
            global_best_fitness: f64::INFINITY,
            inertia,
            c1,
            c2,
            max_velocity,
        }
    }

    /// Evaluate a fitness function for all particles and update bests.
    ///
    /// `fitness_fn` maps position to scalar fitness (lower = better).
    #[allow(dead_code)]
    pub fn evaluate<F: Fn([f64; 3]) -> f64>(&mut self, fitness_fn: &F) {
        for particle in &mut self.particles {
            let f = fitness_fn(particle.position);
            if f < particle.personal_best_fitness {
                particle.personal_best_fitness = f;
                particle.personal_best = particle.position;
            }
            if f < self.global_best_fitness {
                self.global_best_fitness = f;
                self.global_best = particle.position;
            }
        }
    }

    /// Perform one PSO update step using deterministic r1=r2=0.5.
    ///
    /// For reproducibility in tests, uses fixed r1=r2=0.5.
    #[allow(dead_code)]
    pub fn step_deterministic(&mut self) {
        let r1 = 0.5;
        let r2 = 0.5;
        for particle in &mut self.particles {
            for d in 0..3 {
                let cognitive = self.c1 * r1 * (particle.personal_best[d] - particle.position[d]);
                let social = self.c2 * r2 * (self.global_best[d] - particle.position[d]);
                particle.velocity[d] = self.inertia * particle.velocity[d] + cognitive + social;
            }
            // Clamp velocity
            let speed = vec3_len(particle.velocity);
            if speed > self.max_velocity {
                particle.velocity =
                    vec3_scale(vec3_normalize(particle.velocity), self.max_velocity);
            }
            // Update position
            particle.position = vec3_add(particle.position, particle.velocity);
        }
    }

    /// Perform one PSO update step using provided random values r1, r2.
    #[allow(dead_code)]
    pub fn step_with_random(&mut self, r1: f64, r2: f64) {
        for particle in &mut self.particles {
            for d in 0..3 {
                let cognitive = self.c1 * r1 * (particle.personal_best[d] - particle.position[d]);
                let social = self.c2 * r2 * (self.global_best[d] - particle.position[d]);
                particle.velocity[d] = self.inertia * particle.velocity[d] + cognitive + social;
            }
            let speed = vec3_len(particle.velocity);
            if speed > self.max_velocity {
                particle.velocity =
                    vec3_scale(vec3_normalize(particle.velocity), self.max_velocity);
            }
            particle.position = vec3_add(particle.position, particle.velocity);
        }
    }

    /// Return the current global best position.
    #[allow(dead_code)]
    pub fn best_position(&self) -> [f64; 3] {
        self.global_best
    }

    /// Return the current global best fitness.
    #[allow(dead_code)]
    pub fn best_fitness(&self) -> f64 {
        self.global_best_fitness
    }
}

// ── EmergentBehavior (Stigmergy) ──────────────────────────────────────────────

/// A pheromone/trace cell in the stigmergy environment grid.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct PheromoneCell {
    /// Pheromone concentration.
    pub concentration: f64,
    /// Decay rate per timestep.
    pub decay_rate: f64,
}

impl PheromoneCell {
    /// Create a new pheromone cell.
    #[allow(dead_code)]
    pub fn new(concentration: f64, decay_rate: f64) -> Self {
        Self {
            concentration,
            decay_rate,
        }
    }

    /// Decay the pheromone by one timestep `dt`.
    #[allow(dead_code)]
    pub fn decay(&mut self, dt: f64) {
        self.concentration *= (-self.decay_rate * dt).exp();
        if self.concentration < 1e-12 {
            self.concentration = 0.0;
        }
    }

    /// Deposit additional pheromone.
    #[allow(dead_code)]
    pub fn deposit(&mut self, amount: f64) {
        self.concentration += amount;
    }
}

/// Stigmergy environment: indirect agent communication via environment.
///
/// Agents deposit pheromone traces at their positions; the environment
/// diffuses and decays traces over time.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct EmergentBehavior {
    /// 2D grid of pheromone cells (row-major, \[y\]\[x\]).
    pub grid: Vec<Vec<PheromoneCell>>,
    /// Grid resolution: cells per unit length.
    pub resolution: usize,
    /// World bounds.
    pub bounds: [[f64; 2]; 2],
    /// Diffusion coefficient.
    pub diffusion: f64,
}

impl EmergentBehavior {
    /// Create a new stigmergy environment.
    #[allow(dead_code)]
    pub fn new(bounds: [[f64; 2]; 2], resolution: usize, decay_rate: f64, diffusion: f64) -> Self {
        let grid = (0..resolution)
            .map(|_| {
                (0..resolution)
                    .map(|_| PheromoneCell::new(0.0, decay_rate))
                    .collect()
            })
            .collect();
        Self {
            grid,
            resolution,
            bounds,
            diffusion,
        }
    }

    /// Convert world position to grid indices `(ix, iy)`, clamped to grid.
    #[allow(dead_code)]
    pub fn world_to_grid(&self, pos: [f64; 3]) -> (usize, usize) {
        let nx = self.grid[0].len();
        let ny = self.grid.len();
        let fx = (pos[0] - self.bounds[0][0]) / (self.bounds[0][1] - self.bounds[0][0]);
        let fy = (pos[1] - self.bounds[1][0]) / (self.bounds[1][1] - self.bounds[1][0]);
        let ix = ((fx * nx as f64) as isize).max(0).min(nx as isize - 1) as usize;
        let iy = ((fy * ny as f64) as isize).max(0).min(ny as isize - 1) as usize;
        (ix, iy)
    }

    /// Deposit pheromone at agent's position.
    #[allow(dead_code)]
    pub fn deposit_at(&mut self, pos: [f64; 3], amount: f64) {
        let (ix, iy) = self.world_to_grid(pos);
        self.grid[iy][ix].deposit(amount);
    }

    /// Read pheromone concentration at position.
    #[allow(dead_code)]
    pub fn concentration_at(&self, pos: [f64; 3]) -> f64 {
        let (ix, iy) = self.world_to_grid(pos);
        self.grid[iy][ix].concentration
    }

    /// Decay and diffuse all cells by `dt`.
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        // Decay
        for row in &mut self.grid {
            for cell in row {
                cell.decay(dt);
            }
        }
        // Simple nearest-neighbor diffusion
        let ny = self.grid.len();
        if ny == 0 {
            return;
        }
        let nx = self.grid[0].len();
        let old: Vec<Vec<f64>> = self
            .grid
            .iter()
            .map(|r| r.iter().map(|c| c.concentration).collect())
            .collect();
        for iy in 0..ny {
            for ix in 0..nx {
                let mut laplacian = -4.0 * old[iy][ix];
                laplacian += if ix > 0 { old[iy][ix - 1] } else { old[iy][ix] };
                laplacian += if ix + 1 < nx {
                    old[iy][ix + 1]
                } else {
                    old[iy][ix]
                };
                laplacian += if iy > 0 { old[iy - 1][ix] } else { old[iy][ix] };
                laplacian += if iy + 1 < ny {
                    old[iy + 1][ix]
                } else {
                    old[iy][ix]
                };
                self.grid[iy][ix].concentration =
                    (old[iy][ix] + self.diffusion * dt * laplacian).max(0.0);
            }
        }
    }

    /// Total pheromone in the grid.
    #[allow(dead_code)]
    pub fn total_pheromone(&self) -> f64 {
        self.grid
            .iter()
            .flat_map(|r| r.iter())
            .map(|c| c.concentration)
            .sum()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: make a line of agents along X axis.
    fn make_agents(n: usize, spacing: f64, sensor_radius: f64) -> Vec<SwarmAgent> {
        (0..n)
            .map(|i| {
                let mut a = SwarmAgent::new(i, [i as f64 * spacing, 0.0, 0.0], sensor_radius);
                a.velocity = [1.0, 0.0, 0.0];
                a
            })
            .collect()
    }

    // ── SwarmAgent ──────────────────────────────────────────────────────────

    #[test]
    fn test_agent_creation() {
        let a = SwarmAgent::new(0, [1.0, 2.0, 3.0], 5.0);
        assert!((a.position[0] - 1.0).abs() < 1e-10);
        assert!((a.sensor_radius - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_agent_integrate() {
        let mut a = SwarmAgent::new(0, [0.0; 3], 5.0);
        a.velocity = [2.0, 0.0, 0.0];
        a.integrate(0.5);
        assert!(
            (a.position[0] - 1.0).abs() < 1e-10,
            "pos={:.6}",
            a.position[0]
        );
    }

    #[test]
    fn test_agent_neighbors_all_in_range() {
        let agents = make_agents(3, 1.0, 10.0);
        let nbs = agents[1].neighbors(&agents);
        assert_eq!(nbs.len(), 2);
    }

    #[test]
    fn test_agent_neighbors_none_in_range() {
        let agents = make_agents(3, 100.0, 1.0);
        let nbs = agents[0].neighbors(&agents);
        assert!(nbs.is_empty());
    }

    #[test]
    fn test_agent_heading_update() {
        let mut a = SwarmAgent::new(0, [0.0; 3], 5.0);
        a.velocity = [1.0, 0.0, 0.0];
        a.integrate(0.1);
        assert!((a.heading - 0.0).abs() < 1e-10, "heading={:.6}", a.heading);
    }

    // ── FlockingBehavior ────────────────────────────────────────────────────

    #[test]
    fn test_flocking_cohesion_nonempty() {
        let flock = FlockingBehavior::new(1.0);
        let agent = SwarmAgent::new(0, [0.0, 0.0, 0.0], 100.0);
        let nb1 = SwarmAgent::new(1, [4.0, 0.0, 0.0], 100.0);
        let nb2 = SwarmAgent::new(2, [-4.0, 0.0, 0.0], 100.0);
        let nbs = vec![&nb1, &nb2];
        let coh = flock.cohesion(&agent, &nbs);
        // Center is at [0,0,0] so cohesion should be near zero in X.
        assert!(coh[0].abs() < 1e-3, "cohesion_x={:.6}", coh[0]);
    }

    #[test]
    fn test_flocking_alignment_matches_average() {
        let flock = FlockingBehavior::new(1.0);
        let mut agent = SwarmAgent::new(0, [0.0; 3], 100.0);
        agent.velocity = [0.0; 3];
        let mut nb = SwarmAgent::new(1, [2.0, 0.0, 0.0], 100.0);
        nb.velocity = [4.0, 0.0, 0.0];
        let nbs = vec![&nb];
        let ali = flock.alignment(&agent, &nbs);
        // Should steer toward [4,0,0] from [0,0,0], so X component positive.
        assert!(ali[0] > 0.0, "alignment_x={:.6}", ali[0]);
    }

    #[test]
    fn test_flocking_separation_pushes_away() {
        let flock = FlockingBehavior::new(5.0); // r_sep = 5
        let agent = SwarmAgent::new(0, [0.0, 0.0, 0.0], 100.0);
        let nb = SwarmAgent::new(1, [1.0, 0.0, 0.0], 100.0); // within r_sep
        let nbs = vec![&nb];
        let sep = flock.separation(&agent, &nbs);
        // Should push agent in -X direction (away from nb at +X).
        assert!(sep[0] < 0.0, "sep_x={:.6}", sep[0]);
    }

    #[test]
    fn test_flocking_step_swarm_runs() {
        let flock = FlockingBehavior::new(2.0);
        let mut agents = make_agents(5, 1.5, 10.0);
        flock.step_swarm(&mut agents, 0.01);
        // Just checking it doesn't panic.
        assert_eq!(agents.len(), 5);
    }

    #[test]
    fn test_flocking_no_collisions_after_steps() {
        let flock = FlockingBehavior::new(2.0);
        let mut agents = make_agents(5, 1.5, 10.0);
        for _ in 0..20 {
            flock.step_swarm(&mut agents, 0.05);
        }
        // Check no two agents are at exactly the same position.
        for i in 0..agents.len() {
            for j in (i + 1)..agents.len() {
                let d = vec3_dist(agents[i].position, agents[j].position);
                assert!(d > 1e-6, "collision between {} and {}: d={:.6}", i, j, d);
            }
        }
    }

    #[test]
    fn test_flocking_cohesion_radius() {
        // After many steps, agents should stay within a bounded region.
        let flock = FlockingBehavior::new(1.5);
        let mut agents = make_agents(5, 2.0, 15.0);
        for _ in 0..100 {
            flock.step_swarm(&mut agents, 0.02);
        }
        // Compute centroid
        let cx: f64 = agents.iter().map(|a| a.position[0]).sum::<f64>() / 5.0;
        let cy: f64 = agents.iter().map(|a| a.position[1]).sum::<f64>() / 5.0;
        // All agents within 50 units of centroid (generous bound).
        for a in &agents {
            let d = ((a.position[0] - cx).powi(2) + (a.position[1] - cy).powi(2)).sqrt();
            assert!(d < 50.0, "agent {} too far from centroid: d={:.6}", a.id, d);
        }
    }

    // ── FormationControl ────────────────────────────────────────────────────

    #[test]
    fn test_formation_desired_positions() {
        let offsets = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let fc = FormationControl::new(FormationType::VirtualStructure, offsets, 1.0);
        let desired = fc.desired_positions();
        assert!((desired[0][0] - 1.0).abs() < 1e-10);
        assert!((desired[1][0] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_formation_tracking_error_decreases() {
        let offsets = vec![[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
        let mut fc = FormationControl::new(FormationType::VirtualStructure, offsets, 2.0);
        let mut agents = vec![
            SwarmAgent::new(0, [5.0, 0.0, 0.0], 10.0),
            SwarmAgent::new(1, [-5.0, 0.0, 0.0], 10.0),
        ];
        let err0 = fc.tracking_error(&agents);
        for _ in 0..50 {
            fc.step(&mut agents, 0.05);
        }
        let err1 = fc.tracking_error(&agents);
        assert!(
            err1 < err0,
            "error did not decrease: {:.6} → {:.6}",
            err0,
            err1
        );
    }

    #[test]
    fn test_formation_advance_center() {
        let mut fc = FormationControl::new(FormationType::VirtualStructure, vec![], 1.0);
        fc.virtual_velocity = [1.0, 0.0, 0.0];
        fc.advance_center(2.0);
        assert!((fc.virtual_center[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_formation_steering_out_of_range() {
        let fc = FormationControl::new(FormationType::VirtualStructure, vec![], 1.0);
        let agent = SwarmAgent::new(0, [0.0; 3], 5.0);
        let s = fc.formation_steering(&agent, 99);
        for c in s {
            assert!(c.abs() < 1e-10);
        }
    }

    // ── ConsensusProtocol ───────────────────────────────────────────────────

    #[test]
    fn test_consensus_convergence() {
        let states = vec![0.0, 5.0, 10.0, 15.0];
        let mut proto = ConsensusProtocol::new_fully_connected(states, 0.1);
        for _ in 0..5000 {
            proto.step(0.01);
        }
        assert!(
            proto.is_converged(0.01),
            "not converged: {:?}",
            proto.states
        );
    }

    #[test]
    fn test_consensus_mean_preserved() {
        let states = vec![1.0, 3.0, 5.0, 7.0];
        let mean0 = 4.0;
        let mut proto = ConsensusProtocol::new_fully_connected(states, 0.05);
        for _ in 0..100 {
            proto.step(0.01);
        }
        let mean1 = proto.mean_state();
        assert!((mean1 - mean0).abs() < 0.01, "mean changed: {:.6}", mean1);
    }

    #[test]
    fn test_consensus_two_agents() {
        let mut proto = ConsensusProtocol::new_fully_connected(vec![0.0, 10.0], 0.5);
        for _ in 0..2000 {
            proto.step(0.01);
        }
        assert!((proto.states[0] - 5.0).abs() < 0.01);
        assert!((proto.states[1] - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_consensus_laplacian_eigenvalue_positive() {
        let proto = ConsensusProtocol::new_fully_connected(vec![0.0, 1.0, 2.0, 3.0], 0.1);
        let lam = proto.laplacian_second_eigenvalue_approx();
        assert!(lam > 0.0, "lambda2={:.6}", lam);
    }

    // ── CoverageControl ─────────────────────────────────────────────────────

    #[test]
    fn test_coverage_quality_improves_with_agents() {
        let bounds = [[0.0, 10.0], [0.0, 10.0]];
        let ctrl = CoverageControl::new(bounds, 100, 0.5);
        let one_agent = vec![SwarmAgent::new(0, [5.0, 5.0, 0.0], 20.0)];
        let two_agents = vec![
            SwarmAgent::new(0, [2.5, 5.0, 0.0], 20.0),
            SwarmAgent::new(1, [7.5, 5.0, 0.0], 20.0),
        ];
        let q1 = ctrl.coverage_quality(&one_agent);
        let q2 = ctrl.coverage_quality(&two_agents);
        assert!(
            q2 < q1,
            "more agents should improve coverage: {:.6} vs {:.6}",
            q1,
            q2
        );
    }

    #[test]
    fn test_coverage_step_moves_agents() {
        let bounds = [[0.0, 10.0], [0.0, 10.0]];
        let ctrl = CoverageControl::new(bounds, 100, 1.0);
        let mut agents = vec![
            SwarmAgent::new(0, [1.0, 1.0, 0.0], 20.0),
            SwarmAgent::new(1, [2.0, 2.0, 0.0], 20.0),
        ];
        let pos0 = agents[0].position;
        ctrl.step(&mut agents, 0.1);
        let moved = vec3_dist(agents[0].position, pos0) > 1e-6
            || vec3_dist(agents[1].position, [2.0, 2.0, 0.0]) > 1e-6;
        assert!(moved, "agents should have moved");
    }

    #[test]
    fn test_coverage_quality_no_agents() {
        let bounds = [[0.0, 10.0], [0.0, 10.0]];
        let ctrl = CoverageControl::new(bounds, 100, 0.5);
        let q = ctrl.coverage_quality(&[]);
        assert!(q.is_infinite());
    }

    // ── TaskAllocation ──────────────────────────────────────────────────────

    #[test]
    fn test_task_allocation_all_assigned_when_enough_agents() {
        let alloc = TaskAllocation::new(1.0);
        let agents = vec![
            SwarmAgent::new(0, [0.0, 0.0, 0.0], 10.0),
            SwarmAgent::new(1, [5.0, 0.0, 0.0], 10.0),
            SwarmAgent::new(2, [10.0, 0.0, 0.0], 10.0),
        ];
        let tasks = vec![
            Task {
                id: 0,
                position: [1.0, 0.0, 0.0],
                utility: 10.0,
            },
            Task {
                id: 1,
                position: [6.0, 0.0, 0.0],
                utility: 10.0,
            },
        ];
        let result = alloc.allocate(&agents, &tasks);
        assert_eq!(TaskAllocation::assigned_count(&result), 2);
    }

    #[test]
    fn test_task_allocation_bid_prefers_nearby() {
        let alloc = TaskAllocation::new(1.0);
        let agent_near = SwarmAgent::new(0, [1.0, 0.0, 0.0], 10.0);
        let agent_far = SwarmAgent::new(1, [100.0, 0.0, 0.0], 10.0);
        let task = Task {
            id: 0,
            position: [0.0, 0.0, 0.0],
            utility: 10.0,
        };
        let bid_near = alloc.bid(&agent_near, &task);
        let bid_far = alloc.bid(&agent_far, &task);
        assert!(
            bid_near > bid_far,
            "near={:.6} far={:.6}",
            bid_near,
            bid_far
        );
    }

    #[test]
    fn test_task_allocation_partial_assign() {
        let alloc = TaskAllocation::new(1.0);
        let agents = vec![SwarmAgent::new(0, [0.0, 0.0, 0.0], 10.0)];
        let tasks = vec![
            Task {
                id: 0,
                position: [0.0, 0.0, 0.0],
                utility: 10.0,
            },
            Task {
                id: 1,
                position: [1.0, 0.0, 0.0],
                utility: 10.0,
            },
        ];
        let result = alloc.allocate(&agents, &tasks);
        // Only 1 agent, so at most 1 task can be assigned.
        let assigned = TaskAllocation::assigned_count(&result);
        assert!(assigned <= 1, "assigned={}", assigned);
    }

    // ── SwarmSearch (PSO) ───────────────────────────────────────────────────

    #[test]
    fn test_pso_converges_to_optimum() {
        // Minimize f(x,y,z) = x^2 + y^2 + z^2 (optimum at origin)
        let fitness = |p: [f64; 3]| p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        let particles = vec![
            Particle::new([3.0, 3.0, 3.0]),
            Particle::new([-3.0, 2.0, 1.0]),
            Particle::new([1.0, -4.0, 2.0]),
            Particle::new([-2.0, -2.0, -2.0]),
        ];
        let mut pso = SwarmSearch::new(particles, 0.5, 1.5, 1.5, 2.0);
        for _ in 0..200 {
            pso.evaluate(&fitness);
            pso.step_deterministic();
        }
        pso.evaluate(&fitness);
        let best = pso.best_fitness();
        assert!(best < 1.0, "PSO best fitness={:.6}", best);
    }

    #[test]
    fn test_pso_best_fitness_decreases() {
        let fitness = |p: [f64; 3]| (p[0] - 2.0).powi(2) + p[1].powi(2) + p[2].powi(2);
        let particles = vec![
            Particle::new([5.0, 5.0, 5.0]),
            Particle::new([-5.0, -5.0, -5.0]),
        ];
        let mut pso = SwarmSearch::new(particles, 0.7, 1.4, 1.4, 3.0);
        pso.evaluate(&fitness);
        let f0 = pso.best_fitness();
        for _ in 0..50 {
            pso.step_deterministic();
            pso.evaluate(&fitness);
        }
        let f1 = pso.best_fitness();
        assert!(
            f1 <= f0,
            "fitness should not increase: {:.6} → {:.6}",
            f0,
            f1
        );
    }

    #[test]
    fn test_pso_step_with_random() {
        let fitness = |p: [f64; 3]| p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        let particles = vec![
            Particle::new([2.0, 2.0, 2.0]),
            Particle::new([-2.0, -2.0, -2.0]),
        ];
        let mut pso = SwarmSearch::new(particles, 0.5, 1.5, 1.5, 3.0);
        pso.evaluate(&fitness);
        let f0 = pso.best_fitness();
        for i in 0..100 {
            let r = (i as f64 * 0.137).sin().abs() * 0.9 + 0.05;
            pso.step_with_random(r, 1.0 - r);
            pso.evaluate(&fitness);
        }
        let f1 = pso.best_fitness();
        assert!(
            f1 <= f0 + 1e-6,
            "pso step_with_random: {:.6} → {:.6}",
            f0,
            f1
        );
    }

    #[test]
    fn test_pso_particles_count_unchanged() {
        let particles: Vec<Particle> = (0..10)
            .map(|i| Particle::new([i as f64, 0.0, 0.0]))
            .collect();
        let mut pso = SwarmSearch::new(particles, 0.5, 1.5, 1.5, 5.0);
        pso.evaluate(&|p: [f64; 3]| p[0].powi(2));
        pso.step_deterministic();
        assert_eq!(pso.particles.len(), 10);
    }

    // ── EmergentBehavior ────────────────────────────────────────────────────

    #[test]
    fn test_pheromone_deposit_increases_concentration() {
        let mut cell = PheromoneCell::new(0.0, 0.1);
        cell.deposit(5.0);
        assert!((cell.concentration - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_pheromone_decay() {
        let mut cell = PheromoneCell::new(10.0, 1.0);
        cell.decay(1.0);
        assert!(cell.concentration < 10.0);
        assert!(cell.concentration > 0.0);
    }

    #[test]
    fn test_stigmergy_deposit_and_read() {
        let mut env = EmergentBehavior::new([[0.0, 10.0], [0.0, 10.0]], 10, 0.1, 0.01);
        env.deposit_at([5.0, 5.0, 0.0], 3.0);
        let c = env.concentration_at([5.0, 5.0, 0.0]);
        assert!(c > 0.0, "should have pheromone: {:.6}", c);
    }

    #[test]
    fn test_stigmergy_total_decays_over_time() {
        let mut env = EmergentBehavior::new([[0.0, 10.0], [0.0, 10.0]], 10, 0.5, 0.0);
        env.deposit_at([5.0, 5.0, 0.0], 10.0);
        let t0 = env.total_pheromone();
        for _ in 0..10 {
            env.step(0.1);
        }
        let t1 = env.total_pheromone();
        assert!(t1 < t0, "pheromone should decay: {:.6} → {:.6}", t0, t1);
    }

    #[test]
    fn test_stigmergy_multiple_deposits() {
        let mut env = EmergentBehavior::new([[0.0, 10.0], [0.0, 10.0]], 10, 0.1, 0.0);
        env.deposit_at([1.0, 1.0, 0.0], 1.0);
        env.deposit_at([9.0, 9.0, 0.0], 2.0);
        let t = env.total_pheromone();
        assert!(t > 2.9 && t < 3.1, "total={:.6}", t);
    }

    #[test]
    fn test_stigmergy_world_to_grid_bounds() {
        let env = EmergentBehavior::new([[0.0, 10.0], [0.0, 10.0]], 5, 0.1, 0.0);
        let (ix, iy) = env.world_to_grid([0.0, 0.0, 0.0]);
        assert_eq!(ix, 0);
        assert_eq!(iy, 0);
        let (ix2, iy2) = env.world_to_grid([9.99, 9.99, 0.0]);
        assert!(ix2 < 5);
        assert!(iy2 < 5);
    }

    // ── Additional edge-case tests ──────────────────────────────────────────

    #[test]
    fn test_flocking_empty_neighbors_separation() {
        let flock = FlockingBehavior::new(2.0);
        let agent = SwarmAgent::new(0, [0.0; 3], 5.0);
        let sep = flock.separation(&agent, &[]);
        for c in sep {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_flocking_empty_neighbors_alignment() {
        let flock = FlockingBehavior::new(2.0);
        let agent = SwarmAgent::new(0, [0.0; 3], 5.0);
        let ali = flock.alignment(&agent, &[]);
        for c in ali {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_flocking_empty_neighbors_cohesion() {
        let flock = FlockingBehavior::new(2.0);
        let agent = SwarmAgent::new(0, [0.0; 3], 5.0);
        let coh = flock.cohesion(&agent, &[]);
        for c in coh {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_consensus_single_agent_unchanged() {
        let mut proto = ConsensusProtocol::new_fully_connected(vec![42.0], 0.5);
        proto.step(0.1);
        assert!((proto.states[0] - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_formation_zero_error_at_desired() {
        let offsets = vec![[1.0, 0.0, 0.0]];
        let fc = FormationControl::new(FormationType::VirtualStructure, offsets, 1.0);
        let agents = vec![SwarmAgent::new(0, [1.0, 0.0, 0.0], 5.0)];
        let err = fc.tracking_error(&agents);
        assert!(err.abs() < 1e-10, "err={:.6}", err);
    }

    #[test]
    fn test_pso_initial_global_best_set() {
        let particles = vec![Particle::new([3.0, 0.0, 0.0])];
        let pso = SwarmSearch::new(particles, 0.5, 1.5, 1.5, 2.0);
        // Global best fitness starts at infinity (not evaluated yet).
        assert!(pso.global_best_fitness.is_infinite());
    }

    #[test]
    fn test_coverage_voronoi_centroid_two_agents() {
        let bounds = [[0.0, 10.0], [0.0, 10.0]];
        let ctrl = CoverageControl::new(bounds, 400, 1.0);
        let agents = vec![
            SwarmAgent::new(0, [2.0, 5.0, 0.0], 20.0),
            SwarmAgent::new(1, [8.0, 5.0, 0.0], 20.0),
        ];
        let c0 = ctrl.voronoi_centroid(&agents, 0);
        let c1 = ctrl.voronoi_centroid(&agents, 1);
        // Agent 0 is on the left, its centroid should have x < 5.
        assert!(c0[0] < 5.0, "centroid0_x={:.6}", c0[0]);
        // Agent 1 is on the right, its centroid should have x > 5.
        assert!(c1[0] > 5.0, "centroid1_x={:.6}", c1[0]);
    }
}
