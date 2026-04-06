//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use rand::RngExt;

/// A single PSO particle tracking its own best-known position.
#[derive(Debug, Clone)]
pub struct PsoParticle {
    /// Current position in the search space.
    pub position: Vec<f64>,
    /// Current velocity in the search space.
    pub velocity: Vec<f64>,
    /// Personal best position seen so far.
    pub personal_best: Vec<f64>,
    /// Objective value at `personal_best`.
    pub best_value: f64,
}
impl PsoParticle {
    /// Create a particle at `position` with zero velocity.
    pub fn new(position: Vec<f64>) -> Self {
        let dim = position.len();
        let best_value = f64::INFINITY;
        Self {
            personal_best: position.clone(),
            velocity: vec![0.0; dim],
            position,
            best_value,
        }
    }
    /// Update `personal_best` if `value` is better (minimisation).
    pub fn update_personal_best(&mut self, value: f64) {
        if value < self.best_value {
            self.best_value = value;
            self.personal_best = self.position.clone();
        }
    }
}
/// Collective transport system: multiple robots carry a shared payload.
///
/// Agents are arranged in a caging formation around the payload and
/// cooperatively push/pull toward a `goal`.
#[derive(Debug, Clone)]
pub struct CollectiveTransport {
    /// Carrier agents.
    pub agents: Vec<SwarmAgentFull>,
    /// Cage slot assignments.
    pub slots: Vec<CageSlot>,
    /// Payload centre position (m).
    pub payload: [f64; 3],
    /// Payload velocity (m/s).
    pub payload_velocity: [f64; 3],
    /// Payload mass (kg).
    pub payload_mass: f64,
    /// Transport goal position.
    pub goal: [f64; 3],
    /// Proportional gain for carrier slot-tracking.
    pub carrier_gain: f64,
    /// Proportional gain for payload toward goal.
    pub payload_gain: f64,
}
impl CollectiveTransport {
    /// Create a collective transport system with `n_carriers` agents arranged
    /// in a ring around the `payload` at radius `cage_radius`.
    pub fn new(
        n_carriers: usize,
        payload: [f64; 3],
        payload_mass: f64,
        goal: [f64; 3],
        cage_radius: f64,
    ) -> Self {
        use std::f64::consts::PI;
        let mut agents = Vec::with_capacity(n_carriers);
        let mut slots = Vec::with_capacity(n_carriers);
        for k in 0..n_carriers {
            let angle = 2.0 * PI * k as f64 / n_carriers as f64;
            let offset = [cage_radius * angle.cos(), 0.0, cage_radius * angle.sin()];
            let pos = vec3_add(payload, offset);
            agents.push(SwarmAgentFull::new(k as u32, pos, 2.0, 5.0));
            slots.push(CageSlot {
                offset,
                carrier_idx: k,
                load_fraction: 1.0 / n_carriers as f64,
            });
        }
        Self {
            agents,
            slots,
            payload,
            payload_velocity: [0.0; 3],
            payload_mass,
            goal,
            carrier_gain: 2.0,
            payload_gain: 1.0,
        }
    }
    /// Advance one time step: move carriers toward their cage slots, then
    /// update the payload based on aggregate carrier forces.
    pub fn step(&mut self, dt: f64) {
        let mut net_force = [0.0f64; 3];
        for slot in &self.slots {
            let desired_pos = vec3_add(self.payload, slot.offset);
            let i = slot.carrier_idx;
            if i >= self.agents.len() {
                continue;
            }
            let error = vec3_sub(desired_pos, self.agents[i].position);
            let force = vec3_scale(error, self.carrier_gain);
            self.agents[i].apply_steering(force, dt);
            let to_goal = vec3_sub(self.goal, self.payload);
            let contribution = vec3_scale(to_goal, self.payload_gain * slot.load_fraction);
            net_force = vec3_add(net_force, contribution);
        }
        let accel = vec3_scale(net_force, 1.0 / self.payload_mass.max(1e-9));
        self.payload_velocity = vec3_add(self.payload_velocity, vec3_scale(accel, dt));
        self.payload_velocity = vec3_scale(self.payload_velocity, 0.95);
        self.payload = vec3_add(self.payload, vec3_scale(self.payload_velocity, dt));
    }
    /// Distance from payload to goal (m).
    pub fn distance_to_goal(&self) -> f64 {
        vec3_dist(self.payload, self.goal)
    }
    /// Compute the load distribution variance (metric of balance).
    ///
    /// Returns 0.0 for perfectly balanced loads.
    pub fn load_variance(&self) -> f64 {
        let n = self.slots.len();
        if n == 0 {
            return 0.0;
        }
        let mean = 1.0 / n as f64;
        let var: f64 = self
            .slots
            .iter()
            .map(|s| (s.load_fraction - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        var
    }
    /// Rebalance load fractions so each carrier bears an equal share.
    pub fn rebalance_loads(&mut self) {
        let n = self.slots.len();
        if n == 0 {
            return;
        }
        let equal = 1.0 / n as f64;
        for slot in &mut self.slots {
            slot.load_fraction = equal;
        }
    }
}
/// Caging configuration for a payload: each carrier occupies a slot around it.
#[derive(Debug, Clone)]
pub struct CageSlot {
    /// Offset of this slot from the payload centre (m).
    pub offset: [f64; 3],
    /// Index of the carrier agent assigned to this slot.
    pub carrier_idx: usize,
    /// Fraction of the payload load borne by this carrier.
    pub load_fraction: f64,
}
/// A single agent (boid) in a swarm simulation.
#[derive(Debug, Clone)]
pub struct SwarmAgent {
    /// World-space position of the agent `[x, y, z]`.
    pub position: [f64; 3],
    /// Current velocity vector `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Unique agent identifier.
    pub id: u32,
    /// Maximum speed the agent may travel (m/s).
    pub max_speed: f64,
    /// Maximum steering force magnitude (N).
    pub max_force: f64,
}
impl SwarmAgent {
    /// Create a new agent at `position` with zero velocity.
    pub fn new(id: u32, position: [f64; 3], max_speed: f64, max_force: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            id,
            max_speed,
            max_force,
        }
    }
    /// Current speed (magnitude of velocity).
    pub fn speed(&self) -> f64 {
        vec3_norm(self.velocity)
    }
    /// Unit heading vector; returns zero vector if stationary.
    pub fn heading(&self) -> [f64; 3] {
        vec3_normalize(self.velocity)
    }
    /// Apply a steering force for one time-step `dt` (semi-implicit Euler).
    ///
    /// The force is clamped to `max_force` before integration, and the
    /// resulting velocity is clamped to `max_speed`.
    pub fn apply_steering(&mut self, force: [f64; 3], dt: f64) {
        let clamped = vec3_clamp_magnitude(force, self.max_force);
        self.velocity = vec3_add(self.velocity, vec3_scale(clamped, dt));
        self.velocity = vec3_clamp_magnitude(self.velocity, self.max_speed);
        self.position = vec3_add(self.position, vec3_scale(self.velocity, dt));
    }
}
/// Reynolds-style flocking simulation with obstacle avoidance.
///
/// Implements separation, alignment, cohesion, and obstacle avoidance
/// as independent steering behaviours combined by weighted summation.
#[derive(Debug, Clone)]
pub struct BoidFlocking {
    /// All agents in the flock.
    pub agents: Vec<SwarmAgentFull>,
    /// Separation neighbourhood radius (m).
    pub separation_radius: f64,
    /// Alignment neighbourhood radius (m).
    pub alignment_radius: f64,
    /// Cohesion neighbourhood radius (m).
    pub cohesion_radius: f64,
    /// Obstacle avoidance radius (m).
    pub avoidance_radius: f64,
    /// Separation weight.
    pub w_sep: f64,
    /// Alignment weight.
    pub w_ali: f64,
    /// Cohesion weight.
    pub w_coh: f64,
    /// Obstacle avoidance weight.
    pub w_avoid: f64,
    /// Static obstacle positions.
    pub obstacles: Vec<[f64; 3]>,
}
impl BoidFlocking {
    /// Create a new flock with default rule weights.
    pub fn new(agents: Vec<SwarmAgentFull>) -> Self {
        Self {
            agents,
            separation_radius: 2.0,
            alignment_radius: 5.0,
            cohesion_radius: 8.0,
            avoidance_radius: 3.0,
            w_sep: 1.5,
            w_ali: 1.0,
            w_coh: 0.8,
            w_avoid: 2.0,
            obstacles: Vec::new(),
        }
    }
    /// Add a static obstacle at `pos`.
    pub fn add_obstacle(&mut self, pos: [f64; 3]) {
        self.obstacles.push(pos);
    }
    /// Separation steering: steer away from nearby flock-mates.
    fn separation_steer(&self, idx: usize) -> [f64; 3] {
        let a = &self.agents[idx];
        let mut steer = [0.0f64; 3];
        let mut count = 0usize;
        for (j, n) in self.agents.iter().enumerate() {
            if j == idx {
                continue;
            }
            let d = vec3_dist(a.position, n.position);
            if d < self.separation_radius && d > 1e-9 {
                let diff = vec3_normalize(vec3_sub(a.position, n.position));
                steer = vec3_add(steer, vec3_scale(diff, 1.0 / d));
                count += 1;
            }
        }
        if count > 0 {
            steer = vec3_scale(steer, 1.0 / count as f64);
        }
        steer
    }
    /// Alignment steering: steer toward average heading of neighbours.
    fn alignment_steer(&self, idx: usize) -> [f64; 3] {
        let a = &self.agents[idx];
        let mut avg = [0.0f64; 3];
        let mut count = 0usize;
        for (j, n) in self.agents.iter().enumerate() {
            if j == idx {
                continue;
            }
            if vec3_dist(a.position, n.position) < self.alignment_radius {
                avg = vec3_add(avg, n.heading);
                count += 1;
            }
        }
        if count > 0 {
            avg = vec3_scale(avg, 1.0 / count as f64);
            vec3_sub(avg, a.heading)
        } else {
            [0.0; 3]
        }
    }
    /// Cohesion steering: steer toward centroid of neighbours.
    fn cohesion_steer(&self, idx: usize) -> [f64; 3] {
        let a = &self.agents[idx];
        let mut center = [0.0f64; 3];
        let mut count = 0usize;
        for (j, n) in self.agents.iter().enumerate() {
            if j == idx {
                continue;
            }
            if vec3_dist(a.position, n.position) < self.cohesion_radius {
                center = vec3_add(center, n.position);
                count += 1;
            }
        }
        if count > 0 {
            center = vec3_scale(center, 1.0 / count as f64);
            vec3_normalize(vec3_sub(center, a.position))
        } else {
            [0.0; 3]
        }
    }
    /// Obstacle avoidance steering: steer away from static obstacles.
    fn avoidance_steer(&self, idx: usize) -> [f64; 3] {
        let a = &self.agents[idx];
        let mut steer = [0.0f64; 3];
        for &obs in &self.obstacles {
            let d = vec3_dist(a.position, obs);
            if d < self.avoidance_radius && d > 1e-9 {
                let away = vec3_normalize(vec3_sub(a.position, obs));
                steer = vec3_add(
                    steer,
                    vec3_scale(away, (self.avoidance_radius - d) / self.avoidance_radius),
                );
            }
        }
        steer
    }
    /// Advance the flock by one time step `dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.agents.len();
        let forces: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let sep = self.separation_steer(i);
                let ali = self.alignment_steer(i);
                let coh = self.cohesion_steer(i);
                let avoid = self.avoidance_steer(i);
                [
                    self.w_sep * sep[0]
                        + self.w_ali * ali[0]
                        + self.w_coh * coh[0]
                        + self.w_avoid * avoid[0],
                    self.w_sep * sep[1]
                        + self.w_ali * ali[1]
                        + self.w_coh * coh[1]
                        + self.w_avoid * avoid[1],
                    self.w_sep * sep[2]
                        + self.w_ali * ali[2]
                        + self.w_coh * coh[2]
                        + self.w_avoid * avoid[2],
                ]
            })
            .collect();
        for (agent, force) in self.agents.iter_mut().zip(forces.iter()) {
            agent.apply_steering(*force, dt);
        }
    }
    /// Number of agents.
    pub fn len(&self) -> usize {
        self.agents.len()
    }
    /// Returns `true` when there are no agents.
    pub fn is_empty(&self) -> bool {
        self.agents.is_empty()
    }
}
/// Ant Colony Optimizer for the Travelling Salesman Problem (ACO-TSP).
///
/// Uses the Ant System formulation: probabilistic tour construction with
/// pheromone evaporation and deposit.
#[derive(Debug, Clone)]
pub struct AntColonyOptimizer {
    /// Number of nodes (cities).
    pub n_nodes: usize,
    /// Number of ants per iteration.
    pub n_ants: usize,
    /// Pheromone influence exponent α.
    pub alpha: f64,
    /// Heuristic influence exponent β.
    pub beta: f64,
    /// Pheromone evaporation rate ρ ∈ (0, 1).
    pub evaporation: f64,
    /// Pheromone deposit constant Q.
    pub q_deposit: f64,
    /// Pheromone matrix τ\[i\]\[j\].
    pub pheromone: Vec<Vec<f64>>,
    /// Distance matrix d\[i\]\[j\].
    pub distance: Vec<Vec<f64>>,
    /// Best tour found so far (node indices).
    pub best_tour: Vec<usize>,
    /// Length of `best_tour`.
    pub best_length: f64,
    /// Current iteration.
    pub iteration: usize,
}
impl AntColonyOptimizer {
    /// Create an ACO from a symmetric distance matrix.
    ///
    /// Initial pheromone is set to `1.0` on all edges.
    pub fn new(
        distance: Vec<Vec<f64>>,
        n_ants: usize,
        alpha: f64,
        beta: f64,
        evaporation: f64,
    ) -> Self {
        let n = distance.len();
        let pheromone = vec![vec![1.0f64; n]; n];
        let best_tour: Vec<usize> = (0..n).collect();
        let best_length = f64::INFINITY;
        Self {
            n_nodes: n,
            n_ants,
            alpha,
            beta,
            evaporation,
            q_deposit: 100.0,
            pheromone,
            distance,
            best_tour,
            best_length,
            iteration: 0,
        }
    }
    /// Compute tour length for a given tour (node permutation).
    pub fn tour_length(&self, tour: &[usize]) -> f64 {
        let n = tour.len();
        if n < 2 {
            return 0.0;
        }
        let mut total = 0.0;
        for i in 0..n - 1 {
            total += self.distance[tour[i]][tour[i + 1]];
        }
        total += self.distance[tour[n - 1]][tour[0]];
        total
    }
    /// Construct one ant's tour using probabilistic edge selection.
    fn construct_tour(&self, start: usize) -> Vec<usize> {
        let n = self.n_nodes;
        let mut visited = vec![false; n];
        let mut tour = Vec::with_capacity(n);
        let mut current = start;
        visited[current] = true;
        tour.push(current);
        let mut rng = rand::rng();
        for _ in 1..n {
            let mut weights = vec![0.0f64; n];
            let mut total = 0.0f64;
            for j in 0..n {
                if !visited[j] {
                    let tau = self.pheromone[current][j].powf(self.alpha);
                    let eta = if self.distance[current][j] > 1e-12 {
                        (1.0 / self.distance[current][j]).powf(self.beta)
                    } else {
                        1e12
                    };
                    weights[j] = tau * eta;
                    total += weights[j];
                }
            }
            let threshold: f64 = if total > 1e-12 {
                rng.random_range(0.0..total)
            } else {
                0.0
            };
            let mut cumulative = 0.0;
            let mut chosen = n;
            for j in 0..n {
                if !visited[j] {
                    cumulative += weights[j];
                    if cumulative >= threshold {
                        chosen = j;
                        break;
                    }
                }
            }
            if chosen == n {
                chosen = (0..n).find(|&j| !visited[j]).unwrap_or(0);
            }
            visited[chosen] = true;
            tour.push(chosen);
            current = chosen;
        }
        tour
    }
    /// Run one ACO iteration: construct tours, evaporate pheromones, deposit.
    pub fn step(&mut self) {
        let tours: Vec<Vec<usize>> = (0..self.n_ants)
            .map(|k| self.construct_tour(k % self.n_nodes))
            .collect();
        for i in 0..self.n_nodes {
            for j in 0..self.n_nodes {
                self.pheromone[i][j] *= 1.0 - self.evaporation;
                if self.pheromone[i][j] < 1e-6 {
                    self.pheromone[i][j] = 1e-6;
                }
            }
        }
        for tour in &tours {
            let length = self.tour_length(tour);
            let deposit = self.q_deposit / length.max(1e-12);
            let n = tour.len();
            for k in 0..n - 1 {
                self.pheromone[tour[k]][tour[k + 1]] += deposit;
                self.pheromone[tour[k + 1]][tour[k]] += deposit;
            }
            self.pheromone[tour[n - 1]][tour[0]] += deposit;
            self.pheromone[tour[0]][tour[n - 1]] += deposit;
            if length < self.best_length {
                self.best_length = length;
                self.best_tour = tour.clone();
            }
        }
        self.iteration += 1;
    }
    /// Run `n_iter` full ACO iterations.
    pub fn run(&mut self, n_iter: usize) {
        for _ in 0..n_iter {
            self.step();
        }
    }
    /// Mean pheromone level across all edges.
    pub fn mean_pheromone(&self) -> f64 {
        let n = self.n_nodes;
        if n == 0 {
            return 0.0;
        }
        let total: f64 = self
            .pheromone
            .iter()
            .flat_map(|row| row.iter().copied())
            .sum();
        total / (n * n) as f64
    }
}
/// Full Reynolds boid simulation applying separation, alignment, and cohesion.
#[derive(Debug, Clone)]
pub struct BoidSimulation {
    /// All agents in the simulation.
    pub agents: Vec<SwarmAgent>,
    /// Flocking rule parameters.
    pub rules: FlockingRules,
}
impl BoidSimulation {
    /// Create a new simulation with the given agents and rules.
    pub fn new(agents: Vec<SwarmAgent>, rules: FlockingRules) -> Self {
        Self { agents, rules }
    }
    /// Advance all agents by one time step `dt`.
    ///
    /// Forces for all agents are computed from the current positions before
    /// any agent is moved (fully explicit step).
    pub fn step(&mut self, dt: f64) {
        let snapshot: Vec<SwarmAgent> = self.agents.clone();
        let forces: Vec<[f64; 3]> = self
            .agents
            .iter()
            .map(|agent| {
                let sep = separation_force(agent, &snapshot, self.rules.separation_radius);
                let ali = alignment_force(agent, &snapshot, self.rules.alignment_radius);
                let coh = cohesion_force(agent, &snapshot, self.rules.cohesion_radius);
                let fx = self.rules.separation_weight * sep[0]
                    + self.rules.alignment_weight * ali[0]
                    + self.rules.cohesion_weight * coh[0];
                let fy = self.rules.separation_weight * sep[1]
                    + self.rules.alignment_weight * ali[1]
                    + self.rules.cohesion_weight * coh[1];
                let fz = self.rules.separation_weight * sep[2]
                    + self.rules.alignment_weight * ali[2]
                    + self.rules.cohesion_weight * coh[2];
                [fx, fy, fz]
            })
            .collect();
        for (agent, force) in self.agents.iter_mut().zip(forces.iter()) {
            agent.apply_steering(*force, dt);
        }
    }
    /// Number of agents in the simulation.
    pub fn len(&self) -> usize {
        self.agents.len()
    }
    /// Returns `true` if the simulation has no agents.
    pub fn is_empty(&self) -> bool {
        self.agents.is_empty()
    }
}
/// Swarm controller based on artificial potential fields.
///
/// Agents are attracted toward `attractive_goal` and repelled by
/// `repulsive_obstacles`.
#[derive(Debug, Clone)]
pub struct PotentialFieldSwarm {
    /// Target position that attracts all agents.
    pub attractive_goal: [f64; 3],
    /// Obstacle positions that repel agents.
    pub repulsive_obstacles: Vec<[f64; 3]>,
    /// Attractive gain (scales the pull toward the goal).
    pub attractive_gain: f64,
    /// Repulsive gain (scales the push away from obstacles).
    pub repulsive_gain: f64,
    /// Influence radius of each repulsive obstacle (m).
    pub repulsive_radius: f64,
}
impl PotentialFieldSwarm {
    /// Create a new potential field swarm controller.
    pub fn new(
        attractive_goal: [f64; 3],
        repulsive_obstacles: Vec<[f64; 3]>,
        attractive_gain: f64,
        repulsive_gain: f64,
        repulsive_radius: f64,
    ) -> Self {
        Self {
            attractive_goal,
            repulsive_obstacles,
            attractive_gain,
            repulsive_gain,
            repulsive_radius,
        }
    }
    /// Compute the net potential-field force on `agent`.
    ///
    /// Returns `[fx, fy, fz]` in N (unit-mass model).
    pub fn compute_force(&self, agent: &SwarmAgent) -> [f64; 3] {
        let diff_goal = vec3_sub(self.attractive_goal, agent.position);
        let dist_goal = vec3_norm(diff_goal);
        let attractive = if dist_goal > 1e-9 {
            vec3_scale(vec3_normalize(diff_goal), self.attractive_gain)
        } else {
            [0.0; 3]
        };
        let mut repulsive = [0.0f64; 3];
        for &obs in &self.repulsive_obstacles {
            let diff = vec3_sub(agent.position, obs);
            let d = vec3_norm(diff);
            if d < self.repulsive_radius && d > 1e-9 {
                let eta = self.repulsive_gain;
                let r0 = self.repulsive_radius;
                let magnitude = eta * (1.0 / d - 1.0 / r0) / (d * d);
                repulsive = vec3_add(repulsive, vec3_scale(vec3_normalize(diff), magnitude));
            }
        }
        vec3_add(attractive, repulsive)
    }
}
/// Consensus-based formation control: each agent steers toward its assigned
/// target position using a proportional controller.
#[derive(Debug, Clone)]
pub struct FormationControl {
    /// Desired world-space positions for each agent (indexed by agent order).
    pub target_positions: Vec<[f64; 3]>,
    /// Agents under formation control.
    pub agents: Vec<SwarmAgent>,
    /// Proportional gain for position error.
    pub gain: f64,
}
impl FormationControl {
    /// Create a new formation controller.
    ///
    /// `target_positions` and `agents` must have the same length; otherwise
    /// the shorter of the two is used.
    pub fn new(target_positions: Vec<[f64; 3]>, agents: Vec<SwarmAgent>, gain: f64) -> Self {
        Self {
            target_positions,
            agents,
            gain,
        }
    }
    /// Advance one consensus step of duration `dt`.
    ///
    /// Each agent receives a steering force proportional to its position error
    /// relative to its target formation slot.
    pub fn consensus_step(&mut self, dt: f64) {
        let n = self.agents.len().min(self.target_positions.len());
        for i in 0..n {
            let error = vec3_sub(self.target_positions[i], self.agents[i].position);
            let force = vec3_scale(error, self.gain);
            self.agents[i].apply_steering(force, dt);
        }
    }
    /// Mean position error across all agents (m).
    pub fn mean_position_error(&self) -> f64 {
        let n = self.agents.len().min(self.target_positions.len());
        if n == 0 {
            return 0.0;
        }
        let total: f64 = (0..n)
            .map(|i| vec3_dist(self.agents[i].position, self.target_positions[i]))
            .sum();
        total / n as f64
    }
    /// Number of agents in formation.
    pub fn len(&self) -> usize {
        self.agents.len()
    }
    /// Returns `true` if no agents are present.
    pub fn is_empty(&self) -> bool {
        self.agents.is_empty()
    }
}
/// Discrete state of a swarm agent's state machine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AgentState {
    /// Agent is idle and waiting for a task.
    Idle,
    /// Agent is actively moving toward a goal.
    Active,
    /// Agent is following a leader or another agent.
    Following,
    /// Agent is returning to a base or rest position.
    Returning,
    /// Agent has depleted its energy and is inactive.
    Depleted,
}
/// Extended swarm agent with energy, heading, and state machine.
#[derive(Debug, Clone)]
pub struct SwarmAgentFull {
    /// World-space position `[x, y, z]` (m).
    pub position: [f64; 3],
    /// Velocity vector `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
    /// Explicit heading direction (unit vector).  Updated from velocity each step.
    pub heading: [f64; 3],
    /// Current energy level in \[0, 1\]; depletes with motion, replenishes at rest.
    pub energy: f64,
    /// Discrete behavioural state.
    pub state: AgentState,
    /// Unique agent identifier.
    pub id: u32,
    /// Maximum speed (m/s).
    pub max_speed: f64,
    /// Maximum steering force magnitude.
    pub max_force: f64,
    /// Energy consumption rate per unit distance (1/m).
    pub energy_drain: f64,
    /// Energy recovery rate while idle (1/s).
    pub energy_recovery: f64,
}
impl SwarmAgentFull {
    /// Create a new agent at `position` with full energy and `Idle` state.
    pub fn new(id: u32, position: [f64; 3], max_speed: f64, max_force: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            heading: [1.0, 0.0, 0.0],
            energy: 1.0,
            state: AgentState::Idle,
            id,
            max_speed,
            max_force,
            energy_drain: 0.01,
            energy_recovery: 0.05,
        }
    }
    /// Current speed (scalar, m/s).
    pub fn speed(&self) -> f64 {
        vec3_norm(self.velocity)
    }
    /// Apply a steering force for `dt` seconds, update heading and energy.
    ///
    /// If energy drops to zero the agent transitions to `Depleted`.
    pub fn apply_steering(&mut self, force: [f64; 3], dt: f64) {
        if self.state == AgentState::Depleted {
            self.energy = (self.energy + self.energy_recovery * dt).min(1.0);
            if self.energy > 0.2 {
                self.state = AgentState::Idle;
            }
            return;
        }
        let clamped = vec3_clamp_magnitude(force, self.max_force);
        self.velocity = vec3_add(self.velocity, vec3_scale(clamped, dt));
        self.velocity = vec3_clamp_magnitude(self.velocity, self.max_speed);
        let dist = vec3_norm(self.velocity) * dt;
        self.position = vec3_add(self.position, vec3_scale(self.velocity, dt));
        let spd = vec3_norm(self.velocity);
        if spd > 1e-9 {
            self.heading = vec3_normalize(self.velocity);
        }
        self.energy -= self.energy_drain * dist;
        if self.energy <= 0.0 {
            self.energy = 0.0;
            self.state = AgentState::Depleted;
        }
    }
    /// Transition the agent to `state` explicitly.
    pub fn set_state(&mut self, state: AgentState) {
        self.state = state;
    }
}
/// Configuration for the three Reynolds flocking rules.
#[derive(Debug, Clone)]
pub struct FlockingRules {
    /// Neighbourhood radius for separation avoidance (m).
    pub separation_radius: f64,
    /// Neighbourhood radius for alignment matching (m).
    pub alignment_radius: f64,
    /// Neighbourhood radius for cohesion grouping (m).
    pub cohesion_radius: f64,
    /// Weight applied to the separation steering force.
    pub separation_weight: f64,
    /// Weight applied to the alignment steering force.
    pub alignment_weight: f64,
    /// Weight applied to the cohesion steering force.
    pub cohesion_weight: f64,
}
impl FlockingRules {
    /// Sensible defaults for a medium-density swarm.
    pub fn default_rules() -> Self {
        Self {
            separation_radius: 2.0,
            alignment_radius: 5.0,
            cohesion_radius: 8.0,
            separation_weight: 1.5,
            alignment_weight: 1.0,
            cohesion_weight: 0.8,
        }
    }
}
/// Swarm formation controller for V, line, and grid shapes.
///
/// Computes desired positions for each agent in the formation,
/// then applies a proportional controller toward those targets.
#[derive(Debug, Clone)]
pub struct SwarmFormationControl {
    /// Agents in the formation.
    pub agents: Vec<SwarmAgentFull>,
    /// Chosen formation shape.
    pub shape: FormationShape,
    /// Centre of the formation in world space.
    pub center: [f64; 3],
    /// Forward heading of the formation (unit vector).
    pub heading: [f64; 3],
    /// Spacing between adjacent slots (m).
    pub spacing: f64,
    /// Proportional gain for position error.
    pub gain: f64,
    /// Index of the lead agent (slot 0).
    pub leader_idx: usize,
}
impl SwarmFormationControl {
    /// Create a new formation controller.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        agents: Vec<SwarmAgentFull>,
        shape: FormationShape,
        center: [f64; 3],
        heading: [f64; 3],
        spacing: f64,
        gain: f64,
        leader_idx: usize,
    ) -> Self {
        Self {
            agents,
            shape,
            center,
            heading: vec3_normalize(heading),
            spacing,
            gain,
            leader_idx,
        }
    }
    /// Compute the desired slot positions for the current formation.
    ///
    /// Returns one target `[x, y, z]` per agent in order.
    pub fn desired_positions(&self) -> Vec<[f64; 3]> {
        let n = self.agents.len();
        if n == 0 {
            return Vec::new();
        }
        let fwd = self.heading;
        let right = vec3_normalize([-fwd[2], 0.0, fwd[0]]);
        match self.shape {
            FormationShape::Line => {
                let half = (n as f64 - 1.0) * 0.5 * self.spacing;
                (0..n)
                    .map(|i| {
                        let offset = (i as f64) * self.spacing - half;
                        vec3_add(self.center, vec3_scale(right, offset))
                    })
                    .collect()
            }
            FormationShape::V => {
                let back = vec3_scale(fwd, -1.0);
                let mut positions = vec![self.center; n];
                for i in 1..n {
                    let wing = (i as f64) * self.spacing;
                    let side = if i % 2 == 0 { 1.0f64 } else { -1.0f64 };
                    let lateral = vec3_scale(right, side * (i as f64 / 2.0).ceil() * self.spacing);
                    let backward = vec3_scale(back, wing * 0.5);
                    positions[i] = vec3_add(vec3_add(self.center, lateral), backward);
                }
                positions
            }
            FormationShape::Grid => {
                let cols = (n as f64).sqrt().ceil() as usize;
                let rows = n.div_ceil(cols);
                let half_c = (cols as f64 - 1.0) * 0.5 * self.spacing;
                let half_r = (rows as f64 - 1.0) * 0.5 * self.spacing;
                (0..n)
                    .map(|i| {
                        let row = i / cols;
                        let col = i % cols;
                        let lateral = vec3_scale(right, col as f64 * self.spacing - half_c);
                        let longi = vec3_scale(fwd, -(row as f64 * self.spacing - half_r));
                        vec3_add(vec3_add(self.center, lateral), longi)
                    })
                    .collect()
            }
        }
    }
    /// Formation error: RMS distance of agents from their desired slots (m).
    pub fn formation_error(&self) -> f64 {
        let targets = self.desired_positions();
        let n = self.agents.len().min(targets.len());
        if n == 0 {
            return 0.0;
        }
        let sum_sq: f64 = (0..n)
            .map(|i| {
                let d = vec3_dist(self.agents[i].position, targets[i]);
                d * d
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    }
    /// Step all followers toward their formation targets with the given `dt`.
    ///
    /// The leader (index `leader_idx`) is moved manually by the caller;
    /// followers apply proportional steering toward their slots.
    pub fn step(&mut self, dt: f64) {
        let targets = self.desired_positions();
        let n = self.agents.len().min(targets.len());
        for i in 0..n {
            if i == self.leader_idx {
                continue;
            }
            let error = vec3_sub(targets[i], self.agents[i].position);
            let force = vec3_scale(error, self.gain);
            self.agents[i].apply_steering(force, dt);
        }
    }
    /// Number of agents.
    pub fn len(&self) -> usize {
        self.agents.len()
    }
    /// Returns `true` if no agents are in the formation.
    pub fn is_empty(&self) -> bool {
        self.agents.is_empty()
    }
}
/// Formation shape specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormationShape {
    /// V-shaped (wedge) formation.
    V,
    /// Straight line formation.
    Line,
    /// Rectangular grid formation.
    Grid,
}
/// Pheromone trail on a directed edge between two nodes.
#[derive(Debug, Clone)]
pub struct PheromoneEdge {
    /// Source node index.
    pub from: usize,
    /// Destination node index.
    pub to: usize,
    /// Current pheromone level τ (tau).
    pub pheromone: f64,
    /// Heuristic desirability η (eta), e.g. 1/distance.
    pub heuristic: f64,
}
/// Collection of swarm-level statistical metrics.
pub struct SwarmMetrics;
impl SwarmMetrics {
    /// Order parameter Φ ∈ \[0, 1\]: mean cosine similarity of agent headings.
    ///
    /// Φ = |Σ heading_i| / N.  A value of 1 indicates perfect alignment.
    pub fn order_parameter(agents: &[SwarmAgent]) -> f64 {
        let n = agents.len();
        if n == 0 {
            return 0.0;
        }
        let sum_heading = agents
            .iter()
            .fold([0.0f64; 3], |acc, a| vec3_add(acc, a.heading()));
        vec3_norm(sum_heading) / n as f64
    }
    /// Polarization: same as `order_parameter` — magnitude of mean heading.
    ///
    /// Alias retained for domain-specific naming conventions.
    pub fn polarization(agents: &[SwarmAgent]) -> f64 {
        Self::order_parameter(agents)
    }
    /// Dispersion: root-mean-square distance from the swarm centroid (m).
    pub fn dispersion(agents: &[SwarmAgent]) -> f64 {
        let n = agents.len();
        if n == 0 {
            return 0.0;
        }
        let centroid = agents
            .iter()
            .fold([0.0f64; 3], |acc, a| vec3_add(acc, a.position));
        let centroid = vec3_scale(centroid, 1.0 / n as f64);
        let sum_sq: f64 = agents
            .iter()
            .map(|a| {
                let d = vec3_dist(a.position, centroid);
                d * d
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    }
}
/// Particle Swarm Optimizer (PSO) using the canonical velocity update rule.
///
/// Minimises an objective function over a continuous search space.
/// Uses inertia weight `w`, cognitive weight `c1`, and social weight `c2`.
#[derive(Debug, Clone)]
pub struct ParticleSwarmOptimizer {
    /// All particles in the swarm.
    pub particles: Vec<PsoParticle>,
    /// Global best position found by the entire swarm.
    pub global_best: Vec<f64>,
    /// Objective value at `global_best`.
    pub global_best_value: f64,
    /// Inertia weight (damping on previous velocity).
    pub inertia: f64,
    /// Cognitive weight (pull toward personal best).
    pub cognitive: f64,
    /// Social weight (pull toward global best).
    pub social: f64,
    /// Dimensionality of the search space.
    pub dim: usize,
    /// Current iteration count.
    pub iteration: usize,
}
impl ParticleSwarmOptimizer {
    /// Create a PSO with `n_particles` particles, each at random positions
    /// uniformly sampled in `[lo, hi]^dim`.
    pub fn new(n_particles: usize, dim: usize, lo: f64, hi: f64) -> Self {
        let mut rng = rand::rng();
        let particles: Vec<PsoParticle> = (0..n_particles)
            .map(|_| {
                let pos: Vec<f64> = (0..dim).map(|_| rng.random_range(lo..hi)).collect();
                let vel: Vec<f64> = (0..dim)
                    .map(|_| rng.random_range(-(hi - lo) * 0.1..(hi - lo) * 0.1))
                    .collect();
                let mut p = PsoParticle::new(pos);
                p.velocity = vel;
                p
            })
            .collect();
        let global_best = particles[0].position.clone();
        Self {
            particles,
            global_best,
            global_best_value: f64::INFINITY,
            inertia: 0.729,
            cognitive: 1.494,
            social: 1.494,
            dim,
            iteration: 0,
        }
    }
    /// Evaluate all particles with `objective`, update personal and global bests.
    pub fn evaluate<F: Fn(&[f64]) -> f64>(&mut self, objective: &F) {
        for p in &mut self.particles {
            let val = objective(&p.position);
            p.update_personal_best(val);
        }
        for p in &self.particles {
            if p.best_value < self.global_best_value {
                self.global_best_value = p.best_value;
                self.global_best = p.personal_best.clone();
            }
        }
    }
    /// Perform one velocity-position update step.
    pub fn step(&mut self) {
        let mut rng = rand::rng();
        let global_best = self.global_best.clone();
        for p in &mut self.particles {
            for d in 0..self.dim {
                let r1: f64 = rng.random_range(0.0..1.0);
                let r2: f64 = rng.random_range(0.0..1.0);
                let cognitive_term = self.cognitive * r1 * (p.personal_best[d] - p.position[d]);
                let social_term = self.social * r2 * (global_best[d] - p.position[d]);
                p.velocity[d] = self.inertia * p.velocity[d] + cognitive_term + social_term;
            }
            for d in 0..self.dim {
                p.position[d] += p.velocity[d];
            }
        }
        self.iteration += 1;
    }
    /// Run `n_iter` iterations, evaluating `objective` each time.
    pub fn run<F: Fn(&[f64]) -> f64>(&mut self, n_iter: usize, objective: &F) {
        for _ in 0..n_iter {
            self.evaluate(objective);
            self.step();
        }
        self.evaluate(objective);
    }
    /// Mean distance of all particles from the global best (convergence metric).
    pub fn mean_distance_to_best(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        let total: f64 = self
            .particles
            .iter()
            .map(|p| {
                let sq: f64 = p
                    .position
                    .iter()
                    .zip(self.global_best.iter())
                    .map(|(x, b)| (x - b).powi(2))
                    .sum();
                sq.sqrt()
            })
            .sum();
        total / n as f64
    }
}
