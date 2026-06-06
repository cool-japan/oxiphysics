// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Swarm robotics simulation module.
//!
//! Provides collective behaviour primitives for multi-agent systems including
//! Reynolds boids flocking, stigmergy-based coordination, response-threshold
//! task allocation, leader-follower formation control, Levy-flight search,
//! quorum-based collective decision making, and whole-swarm analysis metrics.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector2 helper (lightweight 2-D vector, no nalgebra)
// ---------------------------------------------------------------------------

/// A simple 2-D vector used throughout the swarm module.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
}

impl Vec2 {
    /// Create a new vector.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    /// Squared length.
    pub fn len_sq(self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    /// Euclidean length.
    pub fn len(self) -> f64 {
        self.len_sq().sqrt()
    }

    /// Normalise, returning zero vector if length is near-zero.
    pub fn normalised(self) -> Self {
        let l = self.len();
        if l < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / l, self.y / l)
        }
    }

    /// Clamp magnitude to `max_len`.
    pub fn clamped(self, max_len: f64) -> Self {
        let l = self.len();
        if l > max_len && l > 1e-15 {
            let s = max_len / l;
            Self::new(self.x * s, self.y * s)
        } else {
            self
        }
    }

    /// Distance to another point.
    pub fn dist(self, other: Self) -> f64 {
        (self - other).len()
    }

    /// Dot product.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// Angle of the vector in radians (-PI..PI).
    pub fn angle(self) -> f64 {
        self.y.atan2(self.x)
    }

    /// Create a unit vector from an angle.
    pub fn from_angle(rad: f64) -> Self {
        Self::new(rad.cos(), rad.sin())
    }
}

impl std::ops::Add for Vec2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl std::ops::AddAssign for Vec2 {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl std::ops::Sub for Vec2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl std::ops::Mul<f64> for Vec2 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
}

impl std::ops::Div<f64> for Vec2 {
    type Output = Self;
    fn div(self, s: f64) -> Self {
        Self::new(self.x / s, self.y / s)
    }
}

// ---------------------------------------------------------------------------
// AgentState
// ---------------------------------------------------------------------------

/// Behavioural state of a swarm agent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AgentState {
    /// Idle / waiting.
    Idle,
    /// Searching for a target or resource.
    Searching,
    /// Transporting a resource.
    Transporting,
    /// Returning to base.
    Returning,
    /// Assigned to a task.
    Working,
}

// ---------------------------------------------------------------------------
// SwarmAgent
// ---------------------------------------------------------------------------

/// A single agent in the swarm with kinematic state and energy budget.
#[derive(Debug, Clone)]
pub struct SwarmAgent {
    /// Unique identifier.
    pub id: usize,
    /// Position in 2-D space.
    pub pos: Vec2,
    /// Velocity.
    pub vel: Vec2,
    /// Heading angle in radians.
    pub heading: f64,
    /// Current behavioural state.
    pub state: AgentState,
    /// Remaining energy (0 = depleted).
    pub energy: f64,
    /// Maximum speed.
    pub max_speed: f64,
    /// Maximum steering force magnitude.
    pub max_force: f64,
}

impl SwarmAgent {
    /// Create a new agent at the given position with full energy.
    pub fn new(id: usize, pos: Vec2, max_speed: f64, max_force: f64, energy: f64) -> Self {
        Self {
            id,
            pos,
            vel: Vec2::zero(),
            heading: 0.0,
            state: AgentState::Idle,
            energy,
            max_speed,
            max_force,
        }
    }

    /// Apply a steering force for one time-step `dt`.
    pub fn apply_steering(&mut self, force: Vec2, dt: f64) {
        let f = force.clamped(self.max_force);
        self.vel = (self.vel + f * dt).clamped(self.max_speed);
        self.pos += self.vel * dt;
        if self.vel.len() > 1e-10 {
            self.heading = self.vel.angle();
        }
        // Energy drain proportional to speed.
        self.energy = (self.energy - self.vel.len() * dt * 0.01).max(0.0);
    }

    /// Whether the agent has energy remaining.
    pub fn is_alive(&self) -> bool {
        self.energy > 0.0
    }
}

// ---------------------------------------------------------------------------
// ReynoldsBoids
// ---------------------------------------------------------------------------

/// Reynolds-style flocking model with separation, alignment, and cohesion.
#[derive(Debug, Clone)]
pub struct ReynoldsBoids {
    /// Radius for separation check.
    pub separation_radius: f64,
    /// Radius for alignment/cohesion neighbourhood.
    pub neighbour_radius: f64,
    /// Weight for separation force.
    pub w_separation: f64,
    /// Weight for alignment force.
    pub w_alignment: f64,
    /// Weight for cohesion force.
    pub w_cohesion: f64,
}

impl ReynoldsBoids {
    /// Create with sensible defaults.
    pub fn new() -> Self {
        Self {
            separation_radius: 1.5,
            neighbour_radius: 5.0,
            w_separation: 2.0,
            w_alignment: 1.0,
            w_cohesion: 1.0,
        }
    }

    /// Compute separation steering for agent `idx`.
    pub fn separation(&self, agents: &[SwarmAgent], idx: usize) -> Vec2 {
        let me = &agents[idx];
        let mut steer = Vec2::zero();
        let mut count = 0u32;
        for (i, other) in agents.iter().enumerate() {
            if i == idx {
                continue;
            }
            let d = me.pos.dist(other.pos);
            if d < self.separation_radius && d > 1e-10 {
                let diff = (me.pos - other.pos).normalised() / d;
                steer += diff;
                count += 1;
            }
        }
        if count > 0 {
            steer = steer / count as f64;
        }
        steer
    }

    /// Compute alignment steering for agent `idx`.
    pub fn alignment(&self, agents: &[SwarmAgent], idx: usize) -> Vec2 {
        let me = &agents[idx];
        let mut avg_vel = Vec2::zero();
        let mut count = 0u32;
        for (i, other) in agents.iter().enumerate() {
            if i == idx {
                continue;
            }
            if me.pos.dist(other.pos) < self.neighbour_radius {
                avg_vel += other.vel;
                count += 1;
            }
        }
        if count > 0 {
            avg_vel = avg_vel / count as f64;
            let desired = avg_vel.normalised() * me.max_speed;
            (desired - me.vel).clamped(me.max_force)
        } else {
            Vec2::zero()
        }
    }

    /// Compute cohesion steering for agent `idx`.
    pub fn cohesion(&self, agents: &[SwarmAgent], idx: usize) -> Vec2 {
        let me = &agents[idx];
        let mut center = Vec2::zero();
        let mut count = 0u32;
        for (i, other) in agents.iter().enumerate() {
            if i == idx {
                continue;
            }
            if me.pos.dist(other.pos) < self.neighbour_radius {
                center += other.pos;
                count += 1;
            }
        }
        if count > 0 {
            center = center / count as f64;
            let desired = (center - me.pos).normalised() * me.max_speed;
            (desired - me.vel).clamped(me.max_force)
        } else {
            Vec2::zero()
        }
    }

    /// Compute the total boid steering force for agent `idx`.
    pub fn compute_force(&self, agents: &[SwarmAgent], idx: usize) -> Vec2 {
        let sep = self.separation(agents, idx) * self.w_separation;
        let ali = self.alignment(agents, idx) * self.w_alignment;
        let coh = self.cohesion(agents, idx) * self.w_cohesion;
        sep + ali + coh
    }

    /// Step all agents by `dt` using the boid rules.
    pub fn step(&self, agents: &mut [SwarmAgent], dt: f64) {
        let n = agents.len();
        let mut forces = Vec::with_capacity(n);
        for i in 0..n {
            forces.push(self.compute_force(agents, i));
        }
        for (agent, force) in agents.iter_mut().zip(forces.iter()) {
            agent.apply_steering(*force, dt);
        }
    }
}

impl Default for ReynoldsBoids {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// StigmergyMap
// ---------------------------------------------------------------------------

/// Grid-based pheromone map for indirect communication (stigmergy).
#[derive(Debug, Clone)]
pub struct StigmergyMap {
    /// Number of cells along X.
    pub width: usize,
    /// Number of cells along Y.
    pub height: usize,
    /// Size of each cell in world units.
    pub cell_size: f64,
    /// Pheromone concentrations (row-major).
    pub data: Vec<f64>,
    /// Evaporation rate per second (0..1).
    pub evaporation_rate: f64,
    /// Diffusion rate per second.
    pub diffusion_rate: f64,
}

impl StigmergyMap {
    /// Create a new empty pheromone map.
    pub fn new(width: usize, height: usize, cell_size: f64) -> Self {
        Self {
            width,
            height,
            cell_size,
            data: vec![0.0; width * height],
            evaporation_rate: 0.05,
            diffusion_rate: 0.01,
        }
    }

    /// Convert world position to grid coordinates (clamped).
    pub fn world_to_grid(&self, pos: Vec2) -> (usize, usize) {
        let gx = ((pos.x / self.cell_size).floor() as isize)
            .max(0)
            .min(self.width as isize - 1) as usize;
        let gy = ((pos.y / self.cell_size).floor() as isize)
            .max(0)
            .min(self.height as isize - 1) as usize;
        (gx, gy)
    }

    /// Get pheromone value at grid cell.
    pub fn get(&self, gx: usize, gy: usize) -> f64 {
        if gx < self.width && gy < self.height {
            self.data[gy * self.width + gx]
        } else {
            0.0
        }
    }

    /// Set pheromone value at grid cell.
    pub fn set(&mut self, gx: usize, gy: usize, val: f64) {
        if gx < self.width && gy < self.height {
            self.data[gy * self.width + gx] = val;
        }
    }

    /// Deposit pheromone at a world position.
    pub fn deposit(&mut self, pos: Vec2, amount: f64) {
        let (gx, gy) = self.world_to_grid(pos);
        let cur = self.get(gx, gy);
        self.set(gx, gy, cur + amount);
    }

    /// Apply evaporation over time-step `dt`.
    pub fn evaporate(&mut self, dt: f64) {
        let factor = (1.0 - self.evaporation_rate * dt).max(0.0);
        for v in &mut self.data {
            *v *= factor;
        }
    }

    /// Apply simple diffusion (4-connected Laplacian) over time-step `dt`.
    pub fn diffuse(&mut self, dt: f64) {
        let w = self.width;
        let h = self.height;
        let old = self.data.clone();
        let rate = self.diffusion_rate * dt;
        for y in 0..h {
            for x in 0..w {
                let idx = y * w + x;
                let centre = old[idx];
                let mut lap = -4.0 * centre;
                lap += if x > 0 { old[idx - 1] } else { centre };
                lap += if x + 1 < w { old[idx + 1] } else { centre };
                lap += if y > 0 { old[idx - w] } else { centre };
                lap += if y + 1 < h { old[idx + w] } else { centre };
                self.data[idx] = (centre + rate * lap).max(0.0);
            }
        }
    }

    /// Compute the pheromone gradient direction at a world position.
    ///
    /// Returns a normalised direction vector pointing towards higher concentration.
    pub fn gradient(&self, pos: Vec2) -> Vec2 {
        let (gx, gy) = self.world_to_grid(pos);
        let c = self.get(gx, gy);
        let left = if gx > 0 { self.get(gx - 1, gy) } else { c };
        let right = if gx + 1 < self.width {
            self.get(gx + 1, gy)
        } else {
            c
        };
        let down = if gy > 0 { self.get(gx, gy - 1) } else { c };
        let up = if gy + 1 < self.height {
            self.get(gx, gy + 1)
        } else {
            c
        };
        Vec2::new(right - left, up - down).normalised()
    }

    /// Total pheromone in the map.
    pub fn total_pheromone(&self) -> f64 {
        self.data.iter().sum()
    }

    /// Step the map: evaporate then diffuse.
    pub fn step(&mut self, dt: f64) {
        self.evaporate(dt);
        self.diffuse(dt);
    }
}

// ---------------------------------------------------------------------------
// TaskAllocation
// ---------------------------------------------------------------------------

/// A task that can be assigned to agents.
#[derive(Debug, Clone)]
pub struct Task {
    /// Unique task identifier.
    pub id: usize,
    /// Priority / stimulus level.
    pub stimulus: f64,
    /// Required number of workers.
    pub required_workers: usize,
    /// Currently assigned agent ids.
    pub assigned: Vec<usize>,
    /// Whether the task is completed.
    pub completed: bool,
}

impl Task {
    /// Create a new task.
    pub fn new(id: usize, stimulus: f64, required_workers: usize) -> Self {
        Self {
            id,
            stimulus,
            required_workers,
            assigned: Vec::new(),
            completed: false,
        }
    }

    /// Whether the task has enough workers.
    pub fn is_fully_assigned(&self) -> bool {
        self.assigned.len() >= self.required_workers
    }
}

/// Response-threshold task allocation model.
///
/// Each agent has a threshold per task type. The probability of accepting
/// a task is `stimulus^2 / (stimulus^2 + threshold^2)`.
#[derive(Debug, Clone)]
pub struct TaskAllocation {
    /// Thresholds per agent per task-type: `thresholds[agent_id][task_type]`.
    pub thresholds: Vec<Vec<f64>>,
    /// Number of task types.
    pub num_task_types: usize,
    /// Threshold plasticity rate (learning speed).
    pub plasticity: f64,
}

impl TaskAllocation {
    /// Create a new allocator for `num_agents` agents and `num_task_types` task types
    /// with a uniform initial threshold.
    pub fn new(num_agents: usize, num_task_types: usize, initial_threshold: f64) -> Self {
        Self {
            thresholds: vec![vec![initial_threshold; num_task_types]; num_agents],
            num_task_types,
            plasticity: 0.1,
        }
    }

    /// Compute acceptance probability for an agent given stimulus and its threshold.
    pub fn acceptance_probability(stimulus: f64, threshold: f64) -> f64 {
        let s2 = stimulus * stimulus;
        let t2 = threshold * threshold;
        if s2 + t2 < 1e-30 { 0.0 } else { s2 / (s2 + t2) }
    }

    /// Try to assign an agent to a task. Returns true if the agent accepts.
    ///
    /// Uses a deterministic threshold comparison: accepts if probability > 0.5.
    pub fn try_assign(&mut self, agent_id: usize, task_type: usize, stimulus: f64) -> bool {
        if agent_id >= self.thresholds.len() || task_type >= self.num_task_types {
            return false;
        }
        let threshold = self.thresholds[agent_id][task_type];
        let prob = Self::acceptance_probability(stimulus, threshold);
        if prob > 0.5 {
            // Lower threshold (specialisation) when accepting.
            self.thresholds[agent_id][task_type] =
                (threshold - self.plasticity * threshold).max(0.01);
            true
        } else {
            // Raise threshold when rejecting.
            self.thresholds[agent_id][task_type] =
                threshold + self.plasticity * (1.0 - threshold).max(0.0);
            false
        }
    }

    /// Perform one round of dynamic task assignment across agents and tasks.
    ///
    /// Returns a vector of `(agent_id, task_id)` pairs for newly assigned agents.
    pub fn assign_round(
        &mut self,
        agents: &[SwarmAgent],
        tasks: &mut [Task],
    ) -> Vec<(usize, usize)> {
        let mut assignments = Vec::new();
        for task in tasks.iter_mut() {
            if task.completed || task.is_fully_assigned() {
                continue;
            }
            let task_type = task.id % self.num_task_types;
            for agent in agents {
                if task.is_fully_assigned() {
                    break;
                }
                if agent.state != AgentState::Idle || !agent.is_alive() {
                    continue;
                }
                if task.assigned.contains(&agent.id) {
                    continue;
                }
                if self.try_assign(agent.id, task_type, task.stimulus) {
                    task.assigned.push(agent.id);
                    assignments.push((agent.id, task.id));
                }
            }
        }
        assignments
    }
}

// ---------------------------------------------------------------------------
// FormationControl
// ---------------------------------------------------------------------------

/// Formation shape descriptor.
#[derive(Debug, Clone)]
pub struct FormationShape {
    /// Desired offsets relative to the virtual centre for each follower slot.
    pub offsets: Vec<Vec2>,
}

impl FormationShape {
    /// Create a line formation with `n` agents spaced by `spacing`.
    pub fn line(n: usize, spacing: f64) -> Self {
        let offsets = (0..n).map(|i| Vec2::new(i as f64 * spacing, 0.0)).collect();
        Self { offsets }
    }

    /// Create a ring formation with `n` agents at `radius`.
    pub fn ring(n: usize, radius: f64) -> Self {
        let offsets = (0..n)
            .map(|i| {
                let angle = 2.0 * PI * i as f64 / n as f64;
                Vec2::new(radius * angle.cos(), radius * angle.sin())
            })
            .collect();
        Self { offsets }
    }

    /// Create a grid formation with `cols` columns, `rows` rows, spacing `s`.
    pub fn grid(cols: usize, rows: usize, s: f64) -> Self {
        let mut offsets = Vec::with_capacity(cols * rows);
        for r in 0..rows {
            for c in 0..cols {
                offsets.push(Vec2::new(c as f64 * s, r as f64 * s));
            }
        }
        Self { offsets }
    }
}

/// Leader-follower / virtual-structure formation controller.
#[derive(Debug, Clone)]
pub struct FormationControl {
    /// The formation shape.
    pub shape: FormationShape,
    /// Position gain (proportional).
    pub kp: f64,
    /// Velocity damping gain.
    pub kd: f64,
    /// Virtual centre position (used in virtual-structure mode).
    pub centre: Vec2,
    /// Virtual centre heading.
    pub centre_heading: f64,
}

impl FormationControl {
    /// Create a new formation controller.
    pub fn new(shape: FormationShape, kp: f64, kd: f64) -> Self {
        Self {
            shape,
            kp,
            kd,
            centre: Vec2::zero(),
            centre_heading: 0.0,
        }
    }

    /// Rotate a 2-D vector by an angle.
    fn rotate(v: Vec2, angle: f64) -> Vec2 {
        let c = angle.cos();
        let s = angle.sin();
        Vec2::new(v.x * c - v.y * s, v.x * s + v.y * c)
    }

    /// Compute the desired position for follower slot `slot`.
    pub fn desired_position(&self, slot: usize) -> Vec2 {
        if slot < self.shape.offsets.len() {
            let rotated = Self::rotate(self.shape.offsets[slot], self.centre_heading);
            self.centre + rotated
        } else {
            self.centre
        }
    }

    /// Compute a PD steering force pushing `agent` towards its desired slot.
    pub fn compute_force(&self, agent: &SwarmAgent, slot: usize) -> Vec2 {
        let desired = self.desired_position(slot);
        let pos_err = desired - agent.pos;
        let vel_err = Vec2::zero() - agent.vel; // target vel = 0 relative to formation
        pos_err * self.kp + vel_err * self.kd
    }

    /// Compute formation error (RMS of position errors) for a set of agents.
    pub fn formation_error(&self, agents: &[SwarmAgent]) -> f64 {
        let n = agents.len().min(self.shape.offsets.len());
        if n == 0 {
            return 0.0;
        }
        let sum_sq: f64 = agents
            .iter()
            .take(n)
            .enumerate()
            .map(|(slot, a)| {
                let d = self.desired_position(slot);
                (a.pos - d).len_sq()
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    }

    /// Step the virtual centre forward at `speed` along its heading for `dt`.
    pub fn advance_centre(&mut self, speed: f64, dt: f64) {
        let dir = Vec2::from_angle(self.centre_heading);
        self.centre += dir * speed * dt;
    }

    /// Step all agents towards the formation using PD control.
    pub fn step(&self, agents: &mut [SwarmAgent], dt: f64) {
        let n = agents.len().min(self.shape.offsets.len());
        let forces: Vec<_> = agents
            .iter()
            .enumerate()
            .take(n)
            .map(|(slot, agent)| self.compute_force(agent, slot))
            .collect();
        for (slot, force) in forces.into_iter().enumerate() {
            agents[slot].apply_steering(force, dt);
        }
    }
}

// ---------------------------------------------------------------------------
// SwarmSearch (Levy flight + area coverage)
// ---------------------------------------------------------------------------

/// Levy-flight-inspired random search with area coverage tracking.
#[derive(Debug, Clone)]
pub struct SwarmSearch {
    /// Minimum step length for Levy flight.
    pub step_min: f64,
    /// Levy exponent (typically 1 < alpha < 3, 2 = Cauchy-like).
    pub alpha: f64,
    /// Grid for tracking coverage (cell_size = resolution).
    visited: Vec<bool>,
    /// Grid width.
    grid_w: usize,
    /// Grid height.
    grid_h: usize,
    /// Cell size in world units.
    cell_size: f64,
}

impl SwarmSearch {
    /// Create a new search tracker over a `world_w x world_h` area with `resolution` cell size.
    pub fn new(world_w: f64, world_h: f64, resolution: f64, alpha: f64) -> Self {
        let gw = (world_w / resolution).ceil() as usize;
        let gh = (world_h / resolution).ceil() as usize;
        Self {
            step_min: 0.5,
            alpha,
            visited: vec![false; gw * gh],
            grid_w: gw,
            grid_h: gh,
            cell_size: resolution,
        }
    }

    /// Generate a Levy flight step length from a deterministic seed.
    ///
    /// Uses inverse-CDF: `l = step_min * u^(-1/alpha)` where `u` in (0,1].
    pub fn levy_step_length(&self, u: f64) -> f64 {
        let u_clamped = u.clamp(1e-10, 1.0);
        self.step_min * u_clamped.powf(-1.0 / self.alpha)
    }

    /// Mark a world position as visited.
    pub fn mark_visited(&mut self, pos: Vec2) {
        let gx = ((pos.x / self.cell_size).floor() as isize)
            .max(0)
            .min(self.grid_w as isize - 1) as usize;
        let gy = ((pos.y / self.cell_size).floor() as isize)
            .max(0)
            .min(self.grid_h as isize - 1) as usize;
        self.visited[gy * self.grid_w + gx] = true;
    }

    /// Compute the coverage fraction (0..1).
    pub fn coverage_fraction(&self) -> f64 {
        let total = self.visited.len();
        if total == 0 {
            return 0.0;
        }
        let visited = self.visited.iter().filter(|&&v| v).count();
        visited as f64 / total as f64
    }

    /// Generate a search waypoint: Levy flight from `current_pos` at `angle`.
    pub fn generate_waypoint(&self, current_pos: Vec2, angle: f64, u: f64) -> Vec2 {
        let step = self.levy_step_length(u);
        let dir = Vec2::from_angle(angle);
        current_pos + dir * step
    }

    /// Mark all agent positions and return updated coverage.
    pub fn update_coverage(&mut self, agents: &[SwarmAgent]) -> f64 {
        for a in agents {
            self.mark_visited(a.pos);
        }
        self.coverage_fraction()
    }

    /// Reset coverage grid.
    pub fn reset(&mut self) {
        for v in &mut self.visited {
            *v = false;
        }
    }
}

// ---------------------------------------------------------------------------
// CollectiveDecision (quorum sensing)
// ---------------------------------------------------------------------------

/// Quorum-sensing model for collective binary decisions.
///
/// Each agent has an opinion value in \[0,1\]. When enough agents exceed
/// a quorum threshold, the swarm commits to a decision.
#[derive(Debug, Clone)]
pub struct CollectiveDecision {
    /// Quorum fraction required to commit (0..1).
    pub quorum_fraction: f64,
    /// Social influence rate per time step.
    pub influence_rate: f64,
    /// Noise level for opinion perturbation.
    pub noise_level: f64,
    /// Agent opinion values.
    pub opinions: Vec<f64>,
}

impl CollectiveDecision {
    /// Create a new decision model for `n` agents with neutral initial opinions.
    pub fn new(n: usize, quorum_fraction: f64, influence_rate: f64) -> Self {
        Self {
            quorum_fraction,
            influence_rate,
            noise_level: 0.0,
            opinions: vec![0.5; n],
        }
    }

    /// Set opinion of agent `i`.
    pub fn set_opinion(&mut self, i: usize, val: f64) {
        if i < self.opinions.len() {
            self.opinions[i] = val.clamp(0.0, 1.0);
        }
    }

    /// Run one round of social influence (mean-field coupling).
    ///
    /// Each agent's opinion moves towards the group mean.
    pub fn influence_step(&mut self, dt: f64) {
        let n = self.opinions.len();
        if n == 0 {
            return;
        }
        let mean: f64 = self.opinions.iter().sum::<f64>() / n as f64;
        let rate = self.influence_rate * dt;
        for op in &mut self.opinions {
            *op += rate * (mean - *op);
            *op = op.clamp(0.0, 1.0);
        }
    }

    /// Fraction of agents with opinion > 0.5 (voting "yes").
    pub fn yes_fraction(&self) -> f64 {
        if self.opinions.is_empty() {
            return 0.0;
        }
        let yes = self.opinions.iter().filter(|&&o| o > 0.5).count();
        yes as f64 / self.opinions.len() as f64
    }

    /// Fraction of agents with opinion <= 0.5 (voting "no").
    pub fn no_fraction(&self) -> f64 {
        1.0 - self.yes_fraction()
    }

    /// Check whether the quorum has been reached for "yes".
    pub fn quorum_reached_yes(&self) -> bool {
        self.yes_fraction() >= self.quorum_fraction
    }

    /// Check whether the quorum has been reached for "no".
    pub fn quorum_reached_no(&self) -> bool {
        self.no_fraction() >= self.quorum_fraction
    }

    /// Whether any quorum has been reached.
    pub fn is_decided(&self) -> bool {
        self.quorum_reached_yes() || self.quorum_reached_no()
    }

    /// The current decision if quorum reached, None otherwise.
    pub fn decision(&self) -> Option<bool> {
        if self.quorum_reached_yes() {
            Some(true)
        } else if self.quorum_reached_no() {
            Some(false)
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// SwarmAnalysis
// ---------------------------------------------------------------------------

/// Analysis metrics for a swarm.
pub struct SwarmAnalysis;

impl SwarmAnalysis {
    /// Compute the order parameter (polarisation) of the swarm.
    ///
    /// `Phi = |mean velocity direction|` in \[0,1\].
    /// Phi=1 means all agents move in the same direction.
    pub fn order_parameter(agents: &[SwarmAgent]) -> f64 {
        let n = agents.len();
        if n == 0 {
            return 0.0;
        }
        let mut sum = Vec2::zero();
        let mut moving = 0usize;
        for a in agents {
            let speed = a.vel.len();
            if speed > 1e-10 {
                sum += a.vel.normalised();
                moving += 1;
            }
        }
        if moving == 0 {
            return 0.0;
        }
        sum.len() / moving as f64
    }

    /// Compute the mean nearest-neighbour distance.
    pub fn mean_nn_distance(agents: &[SwarmAgent]) -> f64 {
        let n = agents.len();
        if n < 2 {
            return 0.0;
        }
        let mut total = 0.0;
        for i in 0..n {
            let mut min_d = f64::MAX;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let d = agents[i].pos.dist(agents[j].pos);
                if d < min_d {
                    min_d = d;
                }
            }
            total += min_d;
        }
        total / n as f64
    }

    /// Compute a simple clustering coefficient based on adjacency within `radius`.
    ///
    /// For each agent, count the fraction of its neighbours that are also
    /// neighbours of each other. Average over all agents.
    pub fn clustering_coefficient(agents: &[SwarmAgent], radius: f64) -> f64 {
        let n = agents.len();
        if n < 3 {
            return 0.0;
        }
        // Build adjacency
        let mut adj: Vec<Vec<bool>> = vec![vec![false; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                if agents[i].pos.dist(agents[j].pos) < radius {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }
        let mut total_cc = 0.0;
        let mut valid_count = 0usize;
        for i in 0..n {
            let neighbours: Vec<usize> = (0..n).filter(|&j| adj[i][j]).collect();
            let k = neighbours.len();
            if k < 2 {
                continue;
            }
            let mut triangles = 0usize;
            for a in 0..k {
                for b in (a + 1)..k {
                    if adj[neighbours[a]][neighbours[b]] {
                        triangles += 1;
                    }
                }
            }
            let max_edges = k * (k - 1) / 2;
            total_cc += triangles as f64 / max_edges as f64;
            valid_count += 1;
        }
        if valid_count == 0 {
            0.0
        } else {
            total_cc / valid_count as f64
        }
    }

    /// Compute the centroid of the swarm.
    pub fn centroid(agents: &[SwarmAgent]) -> Vec2 {
        if agents.is_empty() {
            return Vec2::zero();
        }
        let mut sum = Vec2::zero();
        for a in agents {
            sum += a.pos;
        }
        sum / agents.len() as f64
    }

    /// Compute the spread (standard deviation of distances from centroid).
    pub fn spread(agents: &[SwarmAgent]) -> f64 {
        let c = Self::centroid(agents);
        let n = agents.len();
        if n == 0 {
            return 0.0;
        }
        let mean_sq: f64 = agents.iter().map(|a| (a.pos - c).len_sq()).sum::<f64>() / n as f64;
        mean_sq.sqrt()
    }

    /// Compute the mean speed of the swarm.
    pub fn mean_speed(agents: &[SwarmAgent]) -> f64 {
        if agents.is_empty() {
            return 0.0;
        }
        let total: f64 = agents.iter().map(|a| a.vel.len()).sum();
        total / agents.len() as f64
    }

    /// Count agents in each state.
    pub fn state_histogram(agents: &[SwarmAgent]) -> [usize; 5] {
        let mut hist = [0usize; 5];
        for a in agents {
            match a.state {
                AgentState::Idle => hist[0] += 1,
                AgentState::Searching => hist[1] += 1,
                AgentState::Transporting => hist[2] += 1,
                AgentState::Returning => hist[3] += 1,
                AgentState::Working => hist[4] += 1,
            }
        }
        hist
    }

    /// Fraction of agents that are alive.
    pub fn alive_fraction(agents: &[SwarmAgent]) -> f64 {
        if agents.is_empty() {
            return 0.0;
        }
        let alive = agents.iter().filter(|a| a.is_alive()).count();
        alive as f64 / agents.len() as f64
    }

    /// Compute convex-hull area approximation using bounding-box area.
    ///
    /// A fast proxy for the actual convex hull area.
    pub fn bounding_box_area(agents: &[SwarmAgent]) -> f64 {
        if agents.is_empty() {
            return 0.0;
        }
        let mut min_x = f64::MAX;
        let mut max_x = f64::MIN;
        let mut min_y = f64::MAX;
        let mut max_y = f64::MIN;
        for a in agents {
            min_x = min_x.min(a.pos.x);
            max_x = max_x.max(a.pos.x);
            min_y = min_y.min(a.pos.y);
            max_y = max_y.max(a.pos.y);
        }
        (max_x - min_x) * (max_y - min_y)
    }
}

// ---------------------------------------------------------------------------
// Obstacle avoidance helper
// ---------------------------------------------------------------------------

/// Simple circular obstacle.
#[derive(Debug, Clone, Copy)]
pub struct Obstacle {
    /// Centre position.
    pub pos: Vec2,
    /// Radius.
    pub radius: f64,
}

impl Obstacle {
    /// Create a new obstacle.
    pub fn new(pos: Vec2, radius: f64) -> Self {
        Self { pos, radius }
    }
}

/// Compute an obstacle-avoidance steering force for an agent.
pub fn obstacle_avoidance_force(
    agent: &SwarmAgent,
    obstacles: &[Obstacle],
    look_ahead: f64,
    avoidance_strength: f64,
) -> Vec2 {
    let ahead = agent.pos + agent.vel.normalised() * look_ahead;
    let mut nearest: Option<(usize, f64)> = None;
    for (i, obs) in obstacles.iter().enumerate() {
        let d = ahead.dist(obs.pos) - obs.radius;
        if d < 0.0 {
            match nearest {
                None => nearest = Some((i, d)),
                Some((_, nd)) if d < nd => nearest = Some((i, d)),
                _ => {}
            }
        }
    }
    match nearest {
        Some((i, _)) => {
            let away = (ahead - obstacles[i].pos).normalised();
            away * avoidance_strength
        }
        None => Vec2::zero(),
    }
}

// ---------------------------------------------------------------------------
// Neighbourhood query utility
// ---------------------------------------------------------------------------

/// Find indices of agents within `radius` of agent `idx`.
pub fn neighbours_within(agents: &[SwarmAgent], idx: usize, radius: f64) -> Vec<usize> {
    let me = &agents[idx];
    agents
        .iter()
        .enumerate()
        .filter(|(i, a)| *i != idx && me.pos.dist(a.pos) < radius)
        .map(|(i, _)| i)
        .collect()
}

/// Compute pairwise distances between all agents.
pub fn pairwise_distances(agents: &[SwarmAgent]) -> Vec<Vec<f64>> {
    let n = agents.len();
    let mut dists = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = agents[i].pos.dist(agents[j].pos);
            dists[i][j] = d;
            dists[j][i] = d;
        }
    }
    dists
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_agents(positions: &[(f64, f64)]) -> Vec<SwarmAgent> {
        positions
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| SwarmAgent::new(i, Vec2::new(x, y), 5.0, 10.0, 100.0))
            .collect()
    }

    // ---- Vec2 ----

    #[test]
    fn vec2_len_and_normalise() {
        let v = Vec2::new(3.0, 4.0);
        assert!((v.len() - 5.0).abs() < 1e-10);
        let n = v.normalised();
        assert!((n.len() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn vec2_clamped() {
        let v = Vec2::new(10.0, 0.0);
        let c = v.clamped(3.0);
        assert!((c.len() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn vec2_from_angle() {
        let v = Vec2::from_angle(0.0);
        assert!((v.x - 1.0).abs() < 1e-10);
        assert!(v.y.abs() < 1e-10);
    }

    // ---- SwarmAgent ----

    #[test]
    fn agent_apply_steering() {
        let mut a = SwarmAgent::new(0, Vec2::zero(), 10.0, 50.0, 100.0);
        a.apply_steering(Vec2::new(10.0, 0.0), 1.0);
        assert!(a.pos.x > 0.0);
        assert!(a.energy < 100.0);
    }

    #[test]
    fn agent_energy_drains() {
        let mut a = SwarmAgent::new(0, Vec2::zero(), 10.0, 50.0, 0.001);
        a.vel = Vec2::new(5.0, 0.0);
        a.apply_steering(Vec2::zero(), 1.0);
        assert!(!a.is_alive() || a.energy < 0.001);
    }

    // ---- ReynoldsBoids ----

    #[test]
    fn boids_cohesion_pulls_together() {
        let agents = make_agents(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]);
        let boids = ReynoldsBoids {
            neighbour_radius: 10.0,
            separation_radius: 0.5,
            w_separation: 0.0,
            w_alignment: 0.0,
            w_cohesion: 1.0,
        };
        let f = boids.cohesion(&agents, 0);
        // Agent 0 at origin should be pulled towards (2,2)
        assert!(f.x > 0.0);
        assert!(f.y > 0.0);
    }

    #[test]
    fn boids_separation_pushes_apart() {
        let agents = make_agents(&[(0.0, 0.0), (0.5, 0.0)]);
        let boids = ReynoldsBoids {
            neighbour_radius: 10.0,
            separation_radius: 2.0,
            w_separation: 1.0,
            w_alignment: 0.0,
            w_cohesion: 0.0,
        };
        let f = boids.separation(&agents, 0);
        // Agent 0 should be pushed away from agent 1 (negative x)
        assert!(f.x < 0.0);
    }

    #[test]
    fn boids_alignment_matches_velocity() {
        let mut agents = make_agents(&[(0.0, 0.0), (1.0, 0.0)]);
        agents[1].vel = Vec2::new(3.0, 0.0);
        let boids = ReynoldsBoids::new();
        let f = boids.alignment(&agents, 0);
        // Agent 0 (stationary) should steer towards agent 1's velocity
        assert!(f.x > 0.0);
    }

    #[test]
    fn boids_step_moves_agents() {
        let mut agents = make_agents(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
        let boids = ReynoldsBoids::new();
        let old_pos = agents[0].pos;
        boids.step(&mut agents, 0.1);
        assert!((agents[0].pos - old_pos).len() > 0.0 || agents[0].vel.len() > 0.0);
    }

    // ---- StigmergyMap ----

    #[test]
    fn stigmergy_deposit_and_read() {
        let mut map = StigmergyMap::new(10, 10, 1.0);
        map.deposit(Vec2::new(5.5, 5.5), 1.0);
        assert!(map.get(5, 5) > 0.0);
    }

    #[test]
    fn stigmergy_evaporation_reduces() {
        let mut map = StigmergyMap::new(10, 10, 1.0);
        map.deposit(Vec2::new(5.5, 5.5), 10.0);
        let before = map.total_pheromone();
        map.evaporate(1.0);
        let after = map.total_pheromone();
        assert!(after < before);
    }

    #[test]
    fn stigmergy_full_evaporation() {
        let mut map = StigmergyMap::new(5, 5, 1.0);
        map.evaporation_rate = 1.0; // 100% per second
        map.deposit(Vec2::new(2.5, 2.5), 10.0);
        map.evaporate(1.0);
        assert!(map.total_pheromone() < 1e-10);
    }

    #[test]
    fn stigmergy_gradient_points_to_source() {
        let mut map = StigmergyMap::new(10, 10, 1.0);
        // Place pheromone at (6,5), query at (5.5,5.5) → grid (5,5)
        // Gradient checks right=(6,5) vs left=(4,5)
        map.set(6, 5, 10.0);
        let grad = map.gradient(Vec2::new(5.5, 5.5));
        // Gradient should point rightward (towards x=6)
        assert!(grad.x > 0.0, "grad.x={}", grad.x);
    }

    // ---- TaskAllocation ----

    #[test]
    fn task_acceptance_high_stimulus() {
        let prob = TaskAllocation::acceptance_probability(10.0, 1.0);
        // 100/(100+1) ≈ 0.99
        assert!(prob > 0.9);
    }

    #[test]
    fn task_acceptance_low_stimulus() {
        let prob = TaskAllocation::acceptance_probability(0.1, 10.0);
        // 0.01/(0.01+100) ≈ 0.0001
        assert!(prob < 0.01);
    }

    #[test]
    fn task_allocation_assigns() {
        let agents = make_agents(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        let mut tasks = vec![Task::new(0, 20.0, 2)]; // high stimulus, need 2 workers
        let mut alloc = TaskAllocation::new(3, 1, 1.0); // low threshold → easy to accept
        let assignments = alloc.assign_round(&agents, &mut tasks);
        assert!(!assignments.is_empty());
    }

    #[test]
    fn task_threshold_plasticity() {
        let mut alloc = TaskAllocation::new(1, 1, 5.0);
        // With stimulus=10, threshold=5: prob = 100/(100+25) = 0.8 > 0.5 → accept
        let accepted = alloc.try_assign(0, 0, 10.0);
        assert!(accepted);
        // Threshold should decrease after acceptance
        assert!(alloc.thresholds[0][0] < 5.0);
    }

    // ---- FormationControl ----

    #[test]
    fn formation_ring_positions() {
        let shape = FormationShape::ring(4, 5.0);
        let fc = FormationControl::new(shape, 1.0, 0.1);
        let p0 = fc.desired_position(0);
        assert!((p0.x - 5.0).abs() < 1e-10);
        assert!(p0.y.abs() < 1e-10);
    }

    #[test]
    fn formation_error_zero_at_target() {
        let shape = FormationShape::line(3, 2.0);
        let fc = FormationControl::new(shape, 1.0, 0.1);
        // Place agents exactly at desired positions
        let agents = make_agents(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)]);
        let err = fc.formation_error(&agents);
        assert!(err < 1e-10, "error={err}");
    }

    #[test]
    fn formation_error_nonzero_when_displaced() {
        let shape = FormationShape::line(2, 2.0);
        let fc = FormationControl::new(shape, 1.0, 0.1);
        let agents = make_agents(&[(0.0, 10.0), (2.0, 10.0)]); // y offset = 10
        let err = fc.formation_error(&agents);
        assert!(err > 5.0);
    }

    #[test]
    fn formation_pd_force_points_toward_target() {
        let shape = FormationShape::line(1, 0.0);
        let fc = FormationControl::new(shape, 1.0, 0.0);
        let agent = SwarmAgent::new(0, Vec2::new(5.0, 0.0), 5.0, 10.0, 100.0);
        let f = fc.compute_force(&agent, 0);
        // Agent is at (5,0), target is at (0,0) => force points left
        assert!(f.x < 0.0);
    }

    // ---- SwarmSearch ----

    #[test]
    fn levy_step_length_increases_with_small_u() {
        let search = SwarmSearch::new(100.0, 100.0, 1.0, 2.0);
        let l1 = search.levy_step_length(0.9);
        let l2 = search.levy_step_length(0.1);
        assert!(l2 > l1, "smaller u should give longer step");
    }

    #[test]
    fn coverage_fraction_starts_zero() {
        let search = SwarmSearch::new(100.0, 100.0, 1.0, 2.0);
        assert!(search.coverage_fraction() < 1e-10);
    }

    #[test]
    fn coverage_fraction_increases() {
        let mut search = SwarmSearch::new(10.0, 10.0, 1.0, 2.0);
        let agents = make_agents(&[(0.5, 0.5), (3.5, 3.5), (7.5, 7.5)]);
        let cov = search.update_coverage(&agents);
        assert!(cov > 0.0);
    }

    #[test]
    fn search_generate_waypoint() {
        let search = SwarmSearch::new(100.0, 100.0, 1.0, 2.0);
        let wp = search.generate_waypoint(Vec2::new(50.0, 50.0), 0.0, 0.5);
        assert!(wp.x > 50.0); // angle=0 -> positive x direction
    }

    // ---- CollectiveDecision ----

    #[test]
    fn quorum_not_reached_with_mixed_opinions() {
        let mut cd = CollectiveDecision::new(10, 0.8, 1.0);
        // Set half above 0.5, half below → neither quorum reached
        for i in 0..5 {
            cd.set_opinion(i, 0.9);
        }
        for i in 5..10 {
            cd.set_opinion(i, 0.1);
        }
        assert!(!cd.is_decided());
    }

    #[test]
    fn quorum_reached_after_strong_opinions() {
        let mut cd = CollectiveDecision::new(10, 0.7, 1.0);
        for i in 0..10 {
            cd.set_opinion(i, 0.9);
        }
        assert!(cd.quorum_reached_yes());
        assert_eq!(cd.decision(), Some(true));
    }

    #[test]
    fn influence_step_converges_opinions() {
        let mut cd = CollectiveDecision::new(4, 0.8, 5.0);
        cd.set_opinion(0, 0.0);
        cd.set_opinion(1, 0.0);
        cd.set_opinion(2, 1.0);
        cd.set_opinion(3, 1.0);
        // Mean = 0.5, all opinions should converge towards 0.5
        for _ in 0..100 {
            cd.influence_step(0.1);
        }
        for op in &cd.opinions {
            assert!((*op - 0.5).abs() < 0.1, "op={op}");
        }
    }

    #[test]
    fn quorum_no_decision() {
        let mut cd = CollectiveDecision::new(10, 0.7, 1.0);
        for i in 0..10 {
            cd.set_opinion(i, 0.1);
        }
        assert!(cd.quorum_reached_no());
        assert_eq!(cd.decision(), Some(false));
    }

    // ---- SwarmAnalysis ----

    #[test]
    fn order_parameter_all_same_direction() {
        let mut agents = make_agents(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        for a in &mut agents {
            a.vel = Vec2::new(1.0, 0.0);
        }
        let phi = SwarmAnalysis::order_parameter(&agents);
        assert!((phi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn order_parameter_opposite_directions() {
        let mut agents = make_agents(&[(0.0, 0.0), (1.0, 0.0)]);
        agents[0].vel = Vec2::new(1.0, 0.0);
        agents[1].vel = Vec2::new(-1.0, 0.0);
        let phi = SwarmAnalysis::order_parameter(&agents);
        assert!(phi < 1e-10, "phi={phi}");
    }

    #[test]
    fn mean_nn_distance_two_agents() {
        let agents = make_agents(&[(0.0, 0.0), (3.0, 4.0)]);
        let d = SwarmAnalysis::mean_nn_distance(&agents);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn clustering_coefficient_triangle() {
        // Three agents forming a complete triangle within radius
        let agents = make_agents(&[(0.0, 0.0), (1.0, 0.0), (0.5, 0.866)]);
        let cc = SwarmAnalysis::clustering_coefficient(&agents, 2.0);
        assert!((cc - 1.0).abs() < 1e-10, "cc={cc}");
    }

    #[test]
    fn centroid_symmetric() {
        let agents = make_agents(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0), (4.0, 4.0)]);
        let c = SwarmAnalysis::centroid(&agents);
        assert!((c.x - 2.0).abs() < 1e-10);
        assert!((c.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn spread_single_agent_zero() {
        let agents = make_agents(&[(5.0, 5.0)]);
        assert!(SwarmAnalysis::spread(&agents) < 1e-10);
    }

    #[test]
    fn mean_speed_stationary() {
        let agents = make_agents(&[(0.0, 0.0), (1.0, 0.0)]);
        assert!(SwarmAnalysis::mean_speed(&agents) < 1e-10);
    }

    #[test]
    fn alive_fraction_all_alive() {
        let agents = make_agents(&[(0.0, 0.0), (1.0, 0.0)]);
        assert!((SwarmAnalysis::alive_fraction(&agents) - 1.0).abs() < 1e-10);
    }

    // ---- Obstacle avoidance ----

    #[test]
    fn obstacle_avoidance_steers_away() {
        let mut agent = SwarmAgent::new(0, Vec2::new(0.0, 0.0), 5.0, 10.0, 100.0);
        agent.vel = Vec2::new(1.0, 0.0);
        // Obstacle at (3,0) radius 3 → ahead at (5,0), dist=2 < radius=3 → d < 0
        let obs = vec![Obstacle::new(Vec2::new(3.0, 0.0), 3.0)];
        let f = obstacle_avoidance_force(&agent, &obs, 5.0, 10.0);
        assert!(f.len() > 0.0, "force len={}", f.len());
    }

    // ---- Utility functions ----

    #[test]
    fn neighbours_within_radius() {
        let agents = make_agents(&[(0.0, 0.0), (1.0, 0.0), (10.0, 0.0)]);
        let nbrs = neighbours_within(&agents, 0, 2.0);
        assert_eq!(nbrs.len(), 1);
        assert_eq!(nbrs[0], 1);
    }

    #[test]
    fn pairwise_distances_symmetric() {
        let agents = make_agents(&[(0.0, 0.0), (3.0, 4.0)]);
        let dists = pairwise_distances(&agents);
        assert!((dists[0][1] - 5.0).abs() < 1e-10);
        assert!((dists[1][0] - 5.0).abs() < 1e-10);
    }
}
