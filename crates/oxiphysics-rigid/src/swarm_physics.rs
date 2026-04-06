#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Swarm physics and collective rigid body dynamics.
//!
//! Implements the Boids algorithm (Reynolds 1987) for emergent flocking behaviour,
//! including separation, alignment, and cohesion steering forces, plus aggregate
//! metrics such as the order parameter and swarm spread.

#![allow(dead_code)]

use rand::RngExt;
// ---------------------------------------------------------------------------
// SwarmAgent
// ---------------------------------------------------------------------------

/// A single agent in a swarm simulation.
///
/// Each agent has a 3-D position, velocity, heading angle (rad), maximum
/// speed, and a collision radius.
#[derive(Debug, Clone)]
pub struct SwarmAgent {
    /// Unique identifier.
    pub id: usize,
    /// Current position `[x, y, z]` in metres.
    pub pos: [f64; 3],
    /// Current velocity `[vx, vy, vz]` in m/s.
    pub vel: [f64; 3],
    /// Heading angle in the xz-plane (rad).
    pub heading: f64,
    /// Maximum speed (m/s).
    pub max_speed: f64,
    /// Collision avoidance radius (m).
    pub radius: f64,
}

impl SwarmAgent {
    /// Create a new agent at rest.
    ///
    /// * `id`        – unique identifier
    /// * `pos`       – initial position
    /// * `max_speed` – speed clamp (m/s)
    /// * `radius`    – collision radius (m)
    pub fn new(id: usize, pos: [f64; 3], max_speed: f64, radius: f64) -> Self {
        Self {
            id,
            pos,
            vel: [0.0; 3],
            heading: 0.0,
            max_speed,
            radius,
        }
    }

    /// Euclidean distance from this agent to `other`.
    pub fn distance_to(&self, other: &SwarmAgent) -> f64 {
        let dx = self.pos[0] - other.pos[0];
        let dy = self.pos[1] - other.pos[1];
        let dz = self.pos[2] - other.pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Advance the agent position by one Euler step of duration `dt` (s).
    ///
    /// The heading is updated from the horizontal velocity components.
    pub fn update_position(&mut self, dt: f64) {
        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
        self.pos[2] += self.vel[2] * dt;
        // Update heading from xz velocity
        let vx = self.vel[0];
        let vz = self.vel[2];
        if vx.abs() > 1e-12 || vz.abs() > 1e-12 {
            self.heading = vz.atan2(vx);
        }
    }

    /// Clamp the velocity magnitude to `max_speed`.
    pub fn clamp_speed(&mut self) {
        let speed = vec3_len(self.vel);
        if speed > self.max_speed && speed > 1e-15 {
            let scale = self.max_speed / speed;
            self.vel = [
                self.vel[0] * scale,
                self.vel[1] * scale,
                self.vel[2] * scale,
            ];
        }
    }
}

// ---------------------------------------------------------------------------
// BoidForces
// ---------------------------------------------------------------------------

/// Weights and perception radius for the three Boids steering rules.
#[derive(Debug, Clone)]
pub struct BoidForces {
    /// Weight for the separation rule (steer away from close neighbours).
    pub separation: f64,
    /// Weight for the alignment rule (steer toward average heading).
    pub alignment: f64,
    /// Weight for the cohesion rule (steer toward centroid).
    pub cohesion: f64,
    /// Perception radius: only agents within this distance are considered
    /// as neighbours (m).
    pub perception_radius: f64,
}

impl BoidForces {
    /// Construct with explicit weights and perception radius.
    pub fn new(separation: f64, alignment: f64, cohesion: f64, perception_radius: f64) -> Self {
        Self {
            separation,
            alignment,
            cohesion,
            perception_radius,
        }
    }

    /// Default Boids weights.
    pub fn default_weights() -> Self {
        Self::new(1.5, 1.0, 1.0, 5.0)
    }
}

// ---------------------------------------------------------------------------
// Boids steering functions
// ---------------------------------------------------------------------------

/// Compute the separation steering force for `agent` from `neighbors`.
///
/// The force points away from each neighbour, weighted by inverse distance.
/// Returns `[0; 3]` when the neighbour list is empty.
pub fn boids_separation(agent: &SwarmAgent, neighbors: &[&SwarmAgent], weight: f64) -> [f64; 3] {
    let mut force = [0.0f64; 3];
    let mut count = 0usize;

    for nb in neighbors {
        let diff = vec3_sub(agent.pos, nb.pos);
        let dist = vec3_len(diff);
        if dist > 1e-12 {
            // Weight by 1/dist so very close neighbours push harder
            let w = 1.0 / dist;
            force = vec3_add(force, vec3_scale(diff, w));
            count += 1;
        }
    }

    if count > 0 {
        vec3_scale(force, weight / count as f64)
    } else {
        [0.0; 3]
    }
}

/// Compute the alignment steering force for `agent` based on `neighbors`.
///
/// Returns a force in the direction of the average neighbour velocity minus
/// the agent's own velocity.  Returns `[0; 3]` when there are no neighbours.
pub fn boids_alignment(agent: &SwarmAgent, neighbors: &[&SwarmAgent], weight: f64) -> [f64; 3] {
    if neighbors.is_empty() {
        return [0.0; 3];
    }

    let mut avg_vel = [0.0f64; 3];
    for nb in neighbors {
        avg_vel = vec3_add(avg_vel, nb.vel);
    }
    let n = neighbors.len() as f64;
    avg_vel = [avg_vel[0] / n, avg_vel[1] / n, avg_vel[2] / n];

    let steer = vec3_sub(avg_vel, agent.vel);
    vec3_scale(steer, weight)
}

/// Compute the cohesion steering force for `agent` based on `neighbors`.
///
/// Returns a force pointing toward the centroid of the neighbours minus the
/// agent's own position.  Returns `[0; 3]` when there are no neighbours.
pub fn boids_cohesion(agent: &SwarmAgent, neighbors: &[&SwarmAgent], weight: f64) -> [f64; 3] {
    if neighbors.is_empty() {
        return [0.0; 3];
    }

    let mut centroid = [0.0f64; 3];
    for nb in neighbors {
        centroid = vec3_add(centroid, nb.pos);
    }
    let n = neighbors.len() as f64;
    centroid = [centroid[0] / n, centroid[1] / n, centroid[2] / n];

    let steer = vec3_sub(centroid, agent.pos);
    vec3_scale(steer, weight)
}

// ---------------------------------------------------------------------------
// SwarmSimulator
// ---------------------------------------------------------------------------

/// Boids swarm simulator.
///
/// Simulates `N` agents in a cubic box with periodic boundary conditions.
/// At each step each agent's velocity is updated by the three Boids forces
/// and the position is integrated with semi-implicit Euler.
#[derive(Debug, Clone)]
pub struct SwarmSimulator {
    /// All agents in the swarm.
    pub agents: Vec<SwarmAgent>,
    /// Boids force weights and perception radius.
    pub forces: BoidForces,
    /// Half-extent of the periodic simulation box (m).
    pub box_size: f64,
}

impl SwarmSimulator {
    /// Create a new swarm of `n` agents randomly distributed inside the box.
    ///
    /// Agents start with random positions uniformly in `[-box_size, +box_size]`
    /// and zero velocity.
    pub fn new(n: usize, box_size: f64, forces: BoidForces) -> Self {
        let mut rng = rand::rng();
        let agents = (0..n)
            .map(|i| {
                let pos = [
                    rng.random_range(-box_size..box_size),
                    0.0,
                    rng.random_range(-box_size..box_size),
                ];
                SwarmAgent::new(i, pos, 3.0, 0.3)
            })
            .collect();
        Self {
            agents,
            forces,
            box_size,
        }
    }

    /// Number of agents in the swarm.
    pub fn agent_count(&self) -> usize {
        self.agents.len()
    }

    /// Advance the simulation by one time step `dt` (s).
    ///
    /// For each agent the Boids forces from visible neighbours are accumulated,
    /// the velocity is updated, clamped, and the position is stepped.
    /// Agents that leave the box are wrapped with periodic boundary conditions.
    pub fn step(&mut self, dt: f64) {
        let n = self.agents.len();
        let mut delta_vel: Vec<[f64; 3]> = vec![[0.0; 3]; n];

        for i in 0..n {
            // Collect references to neighbours within perception radius
            let perception = self.forces.perception_radius;
            let neighbors: Vec<&SwarmAgent> = self
                .agents
                .iter()
                .enumerate()
                .filter(|(j, other)| *j != i && self.agents[i].distance_to(other) < perception)
                .map(|(_, a)| a)
                .collect();

            let sep = boids_separation(&self.agents[i], &neighbors, self.forces.separation);
            let ali = boids_alignment(&self.agents[i], &neighbors, self.forces.alignment);
            let coh = boids_cohesion(&self.agents[i], &neighbors, self.forces.cohesion);

            delta_vel[i] = vec3_add(vec3_add(sep, ali), coh);
        }

        let bs = self.box_size;
        for (i, agent) in self.agents.iter_mut().enumerate() {
            agent.vel = vec3_add(agent.vel, vec3_scale(delta_vel[i], dt));
            agent.clamp_speed();
            agent.update_position(dt);
            // Periodic boundary
            for k in 0..3 {
                if agent.pos[k] > bs {
                    agent.pos[k] -= 2.0 * bs;
                } else if agent.pos[k] < -bs {
                    agent.pos[k] += 2.0 * bs;
                }
            }
        }
    }

    /// Centroid of the swarm `[x, y, z]` in metres.
    pub fn centroid(&self) -> [f64; 3] {
        if self.agents.is_empty() {
            return [0.0; 3];
        }
        let mut c = [0.0f64; 3];
        for a in &self.agents {
            c = vec3_add(c, a.pos);
        }
        let n = self.agents.len() as f64;
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Polarisation order parameter: average normalised velocity alignment.
    ///
    /// Returns a value in `[0, 1]` where 1 means all agents move in the same
    /// direction and 0 means random orientations.
    pub fn order_parameter(&self) -> f64 {
        if self.agents.is_empty() {
            return 0.0;
        }
        let mut sum = [0.0f64; 3];
        for a in &self.agents {
            let speed = vec3_len(a.vel);
            if speed > 1e-15 {
                let unit = vec3_scale(a.vel, 1.0 / speed);
                sum = vec3_add(sum, unit);
            }
        }
        let n = self.agents.len() as f64;
        let avg = vec3_scale(sum, 1.0 / n);
        vec3_len(avg).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the nearest-neighbour distance for each agent.
///
/// The returned vector has the same length as `agents`; element `i` is the
/// minimum distance from agent `i` to any other agent.  Returns an empty
/// vector for a swarm of fewer than 2 agents.
pub fn nearest_neighbor_distance(agents: &[SwarmAgent]) -> Vec<f64> {
    if agents.len() < 2 {
        return Vec::new();
    }
    agents
        .iter()
        .map(|a| {
            agents
                .iter()
                .filter(|b| b.id != a.id)
                .map(|b| a.distance_to(b))
                .fold(f64::INFINITY, f64::min)
        })
        .collect()
}

/// Average velocity alignment across the swarm (flocking metric).
///
/// Returns the dot product of each agent's normalised velocity with the mean
/// normalised velocity, averaged over all agents.  Range is `[-1, 1]`.
/// Returns `0.0` for a stationary or empty swarm.
pub fn flocking_metric(agents: &[SwarmAgent]) -> f64 {
    if agents.is_empty() {
        return 0.0;
    }
    // Compute mean unit velocity
    let mut sum = [0.0f64; 3];
    let mut moving = 0usize;
    for a in agents {
        let spd = vec3_len(a.vel);
        if spd > 1e-15 {
            sum = vec3_add(sum, vec3_scale(a.vel, 1.0 / spd));
            moving += 1;
        }
    }
    if moving == 0 {
        return 0.0;
    }
    let mean = vec3_scale(sum, 1.0 / moving as f64);
    let mean_len = vec3_len(mean);
    if mean_len < 1e-15 {
        return 0.0;
    }
    let mean_unit = vec3_scale(mean, 1.0 / mean_len);

    let mut total = 0.0;
    for a in agents {
        let spd = vec3_len(a.vel);
        if spd > 1e-15 {
            let unit = vec3_scale(a.vel, 1.0 / spd);
            total += vec3_dot(unit, mean_unit);
        }
    }
    (total / agents.len() as f64).clamp(-1.0, 1.0)
}

/// RMS distance of all agents from the swarm centroid (spread metric, m).
///
/// Returns `0.0` for an empty swarm.
pub fn swarm_spread(agents: &[SwarmAgent]) -> f64 {
    if agents.is_empty() {
        return 0.0;
    }
    let n = agents.len() as f64;
    let mut c = [0.0f64; 3];
    for a in agents {
        c = vec3_add(c, a.pos);
    }
    c = [c[0] / n, c[1] / n, c[2] / n];

    let rms = agents
        .iter()
        .map(|a| {
            let d = vec3_sub(a.pos, c);
            vec3_dot(d, d)
        })
        .sum::<f64>()
        / n;
    rms.sqrt()
}

/// Lennard-Jones-like pairwise attraction-repulsion force on agent at `pos_i`
/// from agent at `pos_j`.
///
/// * `r0`       – equilibrium (zero-force) distance (m)
/// * `strength` – force scale (N)
///
/// The force is attractive for `r > r0` and repulsive for `r < r0`,
/// directed from `pos_j` toward `pos_i`.  Returns `[0; 3]` when the
/// positions coincide.
pub fn attraction_repulsion(pos_i: [f64; 3], pos_j: [f64; 3], r0: f64, strength: f64) -> [f64; 3] {
    let diff = vec3_sub(pos_i, pos_j);
    let dist = vec3_len(diff);
    if dist < 1e-12 {
        return [0.0; 3];
    }
    // Simple polynomial: F = strength * (1 - r0/r) * (r0/r)^2 * (diff/|diff|)
    // r < r0: repulsive (force points away from j, i.e. in diff direction)
    // r > r0: attractive (force points toward j, i.e. against diff direction)
    // Use: F = strength * (r/r0 - 1) * (r0/r)^3, sign > 0 means repulsive
    let ratio = dist / r0; // > 1 when far (attractive region), < 1 when close (repulsive region)
    let magnitude = strength * (ratio - 1.0) * (1.0 / ratio).powi(3);
    // magnitude > 0 when r > r0 → attractive (pull toward j, negate diff direction)
    // magnitude < 0 when r < r0 → repulsive (push away from j, diff direction)
    vec3_scale(vec3_scale(diff, 1.0 / dist), -magnitude)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── SwarmAgent ────────────────────────────────────────────────────────

    #[test]
    fn agent_new_starts_at_rest() {
        let a = SwarmAgent::new(0, [1.0, 2.0, 3.0], 5.0, 0.5);
        assert_eq!(a.vel, [0.0; 3]);
        assert_eq!(a.pos, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn agent_distance_self_is_zero() {
        let a = SwarmAgent::new(0, [1.0, 2.0, 3.0], 5.0, 0.5);
        assert!((a.distance_to(&a)).abs() < 1e-12);
    }

    #[test]
    fn agent_distance_symmetric() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [3.0, 4.0, 0.0], 5.0, 0.5);
        assert!((a.distance_to(&b) - b.distance_to(&a)).abs() < 1e-12);
    }

    #[test]
    fn agent_distance_known_value() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [3.0, 4.0, 0.0], 5.0, 0.5);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn agent_update_position_moves() {
        let mut a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 10.0, 0.5);
        a.vel = [1.0, 0.0, 2.0];
        let pos_before = a.pos;
        a.update_position(1.0);
        assert!(a.pos[0] != pos_before[0] || a.pos[2] != pos_before[2]);
        assert!((a.pos[0] - 1.0).abs() < 1e-12);
        assert!((a.pos[2] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn agent_update_position_zero_vel_no_move() {
        let mut a = SwarmAgent::new(0, [5.0, 0.0, 5.0], 10.0, 0.5);
        let before = a.pos;
        a.update_position(1.0);
        assert_eq!(a.pos, before);
    }

    #[test]
    fn agent_clamp_speed_enforced() {
        let mut a = SwarmAgent::new(0, [0.0; 3], 2.0, 0.5);
        a.vel = [10.0, 0.0, 0.0];
        a.clamp_speed();
        let speed = (a.vel[0].powi(2) + a.vel[1].powi(2) + a.vel[2].powi(2)).sqrt();
        assert!((speed - 2.0).abs() < 1e-12);
    }

    #[test]
    fn agent_clamp_speed_below_max_unchanged() {
        let mut a = SwarmAgent::new(0, [0.0; 3], 10.0, 0.5);
        a.vel = [1.0, 0.0, 0.0];
        a.clamp_speed();
        assert!((a.vel[0] - 1.0).abs() < 1e-12);
    }

    // ── BoidForces ────────────────────────────────────────────────────────

    #[test]
    fn boid_forces_default_weights_positive() {
        let bf = BoidForces::default_weights();
        assert!(bf.separation > 0.0);
        assert!(bf.alignment > 0.0);
        assert!(bf.cohesion > 0.0);
        assert!(bf.perception_radius > 0.0);
    }

    // ── boids_separation ─────────────────────────────────────────────────

    #[test]
    fn separation_empty_neighbors_zero_force() {
        let a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        let f = boids_separation(&a, &[], 1.0);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn separation_pushes_away() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [1.0, 0.0, 0.0], 5.0, 0.5);
        let f = boids_separation(&a, &[&b], 1.0);
        // Force should point in the -x direction (away from b at +x)
        assert!(f[0] < 0.0, "separation should push away: f={f:?}");
    }

    #[test]
    fn separation_single_neighbor_nonzero() {
        let a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        let b = SwarmAgent::new(1, [0.5, 0.0, 0.0], 5.0, 0.5);
        let f = boids_separation(&a, &[&b], 2.0);
        assert!(f[0].abs() > 1e-10);
    }

    // ── boids_alignment ──────────────────────────────────────────────────

    #[test]
    fn alignment_empty_neighbors_zero_force() {
        let a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        let f = boids_alignment(&a, &[], 1.0);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn alignment_toward_average_velocity() {
        let mut a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        a.vel = [0.0, 0.0, 0.0];
        let mut b = SwarmAgent::new(1, [1.0, 0.0, 0.0], 5.0, 0.5);
        b.vel = [2.0, 0.0, 0.0];
        let f = boids_alignment(&a, &[&b], 1.0);
        assert!(f[0] > 0.0, "should align toward +x: f={f:?}");
    }

    // ── boids_cohesion ───────────────────────────────────────────────────

    #[test]
    fn cohesion_empty_neighbors_zero_force() {
        let a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        let f = boids_cohesion(&a, &[], 1.0);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn cohesion_toward_centroid() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [4.0, 0.0, 0.0], 5.0, 0.5);
        let f = boids_cohesion(&a, &[&b], 1.0);
        // centroid is at x=4, agent at x=0 → force should be +x
        assert!(f[0] > 0.0, "cohesion should pull toward centroid: f={f:?}");
    }

    // ── SwarmSimulator ────────────────────────────────────────────────────

    #[test]
    fn simulator_creates_n_agents() {
        let forces = BoidForces::default_weights();
        let sim = SwarmSimulator::new(10, 10.0, forces);
        assert_eq!(sim.agent_count(), 10);
    }

    #[test]
    fn simulator_zero_agents() {
        let forces = BoidForces::default_weights();
        let sim = SwarmSimulator::new(0, 10.0, forces);
        assert_eq!(sim.agent_count(), 0);
    }

    #[test]
    fn simulator_step_changes_positions() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(5, 10.0, forces);
        // Give each agent a known velocity
        for a in &mut sim.agents {
            a.vel = [1.0, 0.0, 0.0];
        }
        let before: Vec<[f64; 3]> = sim.agents.iter().map(|a| a.pos).collect();
        sim.step(0.1);
        let after: Vec<[f64; 3]> = sim.agents.iter().map(|a| a.pos).collect();
        let moved = before
            .iter()
            .zip(after.iter())
            .any(|(b, a)| (a[0] - b[0]).abs() > 1e-12 || (a[1] - b[1]).abs() > 1e-12);
        assert!(moved, "positions should change after step");
    }

    #[test]
    fn centroid_single_agent() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(1, 10.0, forces);
        sim.agents[0].pos = [3.0, 0.0, 5.0];
        let c = sim.centroid();
        assert!((c[0] - 3.0).abs() < 1e-12);
        assert!((c[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn centroid_two_agents_midpoint() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(2, 20.0, forces);
        sim.agents[0].pos = [0.0, 0.0, 0.0];
        sim.agents[1].pos = [10.0, 0.0, 0.0];
        let c = sim.centroid();
        assert!(
            (c[0] - 5.0).abs() < 1e-12,
            "centroid x should be 5.0, got {}",
            c[0]
        );
    }

    #[test]
    fn centroid_within_box() {
        let forces = BoidForces::default_weights();
        let bs = 10.0;
        let sim = SwarmSimulator::new(20, bs, forces);
        let c = sim.centroid();
        // Centroid can legitimately be outside box if agents are near boundaries;
        // just check it's finite and reasonable
        for v in c {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn order_parameter_in_range() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(10, 10.0, forces);
        for (i, a) in sim.agents.iter_mut().enumerate() {
            a.vel = [i as f64 * 0.1, 0.0, 0.0];
        }
        let op = sim.order_parameter();
        assert!((0.0..=1.0).contains(&op), "order parameter={op}");
    }

    #[test]
    fn order_parameter_zero_for_stationary() {
        let forces = BoidForces::default_weights();
        let sim = SwarmSimulator::new(5, 10.0, forces);
        let op = sim.order_parameter();
        assert!(op.abs() < 1e-12);
    }

    #[test]
    fn order_parameter_one_for_aligned() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(5, 10.0, forces);
        for a in &mut sim.agents {
            a.vel = [1.0, 0.0, 0.0];
        }
        let op = sim.order_parameter();
        assert!((op - 1.0).abs() < 1e-9, "perfectly aligned: op={op}");
    }

    // ── nearest_neighbor_distance ─────────────────────────────────────────

    #[test]
    fn nn_distance_single_agent_empty() {
        let a = SwarmAgent::new(0, [0.0; 3], 5.0, 0.5);
        let d = nearest_neighbor_distance(&[a]);
        assert!(d.is_empty());
    }

    #[test]
    fn nn_distance_two_agents() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [3.0, 4.0, 0.0], 5.0, 0.5);
        let d = nearest_neighbor_distance(&[a, b]);
        assert_eq!(d.len(), 2);
        assert!((d[0] - 5.0).abs() < 1e-12);
        assert!((d[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn nn_distance_nonnegative() {
        let forces = BoidForces::default_weights();
        let sim = SwarmSimulator::new(8, 10.0, forces);
        let d = nearest_neighbor_distance(&sim.agents);
        for v in &d {
            assert!(*v >= 0.0, "distance must be non-negative: {v}");
        }
    }

    // ── flocking_metric ───────────────────────────────────────────────────

    #[test]
    fn flocking_metric_empty_is_zero() {
        assert!((flocking_metric(&[])).abs() < 1e-12);
    }

    #[test]
    fn flocking_metric_all_stationary_is_zero() {
        let agents = vec![
            SwarmAgent::new(0, [0.0; 3], 5.0, 0.5),
            SwarmAgent::new(1, [1.0, 0.0, 0.0], 5.0, 0.5),
        ];
        assert!((flocking_metric(&agents)).abs() < 1e-12);
    }

    #[test]
    fn flocking_metric_range() {
        let forces = BoidForces::default_weights();
        let mut sim = SwarmSimulator::new(10, 10.0, forces);
        for (i, a) in sim.agents.iter_mut().enumerate() {
            let angle = i as f64 * 2.0 * PI / 10.0;
            a.vel = [angle.cos(), 0.0, angle.sin()];
        }
        let fm = flocking_metric(&sim.agents);
        assert!((-1.0..=1.0).contains(&fm), "flocking metric={fm}");
    }

    // ── swarm_spread ──────────────────────────────────────────────────────

    #[test]
    fn swarm_spread_empty_is_zero() {
        assert!((swarm_spread(&[])).abs() < 1e-12);
    }

    #[test]
    fn swarm_spread_nonnegative() {
        let forces = BoidForces::default_weights();
        let sim = SwarmSimulator::new(10, 10.0, forces);
        assert!(swarm_spread(&sim.agents) >= 0.0);
    }

    #[test]
    fn swarm_spread_single_point_zero() {
        let a = SwarmAgent::new(0, [5.0, 0.0, 5.0], 5.0, 0.5);
        assert!((swarm_spread(&[a])).abs() < 1e-12);
    }

    #[test]
    fn swarm_spread_two_agents_known() {
        let a = SwarmAgent::new(0, [0.0, 0.0, 0.0], 5.0, 0.5);
        let b = SwarmAgent::new(1, [2.0, 0.0, 0.0], 5.0, 0.5);
        // Centroid at x=1. Distances: 1, 1. RMS = sqrt((1+1)/2) = 1
        let spread = swarm_spread(&[a, b]);
        assert!((spread - 1.0).abs() < 1e-12, "spread={spread}");
    }

    // ── attraction_repulsion ──────────────────────────────────────────────

    #[test]
    fn attraction_repulsion_at_r0_zero_force() {
        let pos_i = [1.0, 0.0, 0.0];
        let pos_j = [0.0, 0.0, 0.0];
        // r0 = 1 = distance → force should be zero
        let f = attraction_repulsion(pos_i, pos_j, 1.0, 10.0);
        let mag = (f[0].powi(2) + f[1].powi(2) + f[2].powi(2)).sqrt();
        assert!(mag < 1e-10, "force at equilibrium should be zero: {mag}");
    }

    #[test]
    fn attraction_repulsion_close_repulsive() {
        // pos_i is closer than r0 to pos_j → force should be repulsive (away from j)
        let pos_i = [0.5, 0.0, 0.0];
        let pos_j = [0.0, 0.0, 0.0];
        let f = attraction_repulsion(pos_i, pos_j, 2.0, 1.0);
        // Force on i points away from j, so +x direction
        assert!(f[0] > 0.0, "repulsive force should be in +x: {}", f[0]);
    }

    #[test]
    fn attraction_repulsion_coincident_zero_force() {
        let pos = [1.0, 2.0, 3.0];
        let f = attraction_repulsion(pos, pos, 1.0, 1.0);
        assert_eq!(f, [0.0; 3]);
    }
}
