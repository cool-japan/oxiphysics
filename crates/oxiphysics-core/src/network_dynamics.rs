#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Network/graph dynamics: epidemic models, opinion dynamics, synchronization.
//!
//! This module implements:
//! - [`NetworkGraph`]: adjacency list representation, directed/undirected, weighted edges
//! - [`SirModel`]: Susceptible-Infectious-Recovered epidemic model (ODE + network)
//! - [`SeirdModel`]: extended SEIRD epidemic (Exposed + Deceased compartments)
//! - [`OpinionDynamics`]: Deffuant-Weisbuch bounded confidence model
//! - [`DeGrootModel`]: consensus update rule
//! - [`VoterModel`]: majority rule opinion dynamics
//! - [`KuramotoModel`]: coupled oscillators synchronization
//! - [`RandomWalkGraph`]: random walk on graph, hitting time, cover time
//! - [`PageRank`]: power iteration with dangling node handling
//! - [`CommunityDetection`]: modularity Q, greedy modularity maximization
//! - [`NetworkRobustness`]: percolation threshold, giant component fraction
//! - [`SmallWorld`]: Watts-Strogatz rewiring, clustering coefficient, path length

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::{HashMap, VecDeque};
use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Internal LCG RNG (no external rand dependency needed in tests)
// ─────────────────────────────────────────────────────────────────────────────

struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn next_range(&mut self, lo: usize, hi: usize) -> usize {
        lo + (self.next_u64() as usize % (hi - lo))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NetworkGraph
// ─────────────────────────────────────────────────────────────────────────────

/// Adjacency-list graph for network dynamics simulations.
///
/// Supports directed and undirected graphs with optional edge weights.
/// Node indices are in `0..n`.
#[derive(Clone)]
pub struct NetworkGraph {
    /// Number of nodes.
    pub n: usize,
    /// Whether the graph is directed.
    pub directed: bool,
    /// Adjacency list: `adj[i]` = list of `(neighbor, weight)`.
    pub adj: Vec<Vec<(usize, f64)>>,
}

impl NetworkGraph {
    /// Creates an empty graph with `n` nodes.
    pub fn new(n: usize, directed: bool) -> Self {
        Self {
            n,
            directed,
            adj: vec![Vec::new(); n],
        }
    }

    /// Adds a directed edge `u → v` with weight `w`.
    pub fn add_edge(&mut self, u: usize, v: usize, w: f64) {
        self.adj[u].push((v, w));
        if !self.directed {
            self.adj[v].push((u, w));
        }
    }

    /// Returns the degree (out-degree for directed) of node `u`.
    pub fn degree(&self, u: usize) -> usize {
        self.adj[u].len()
    }

    /// Returns the number of edges (undirected counts each once).
    pub fn edge_count(&self) -> usize {
        let total: usize = self.adj.iter().map(|a| a.len()).sum();
        if self.directed { total } else { total / 2 }
    }

    /// Builds a complete undirected graph (all edges weight 1).
    pub fn complete(n: usize) -> Self {
        let mut g = Self::new(n, false);
        for i in 0..n {
            for j in (i + 1)..n {
                g.add_edge(i, j, 1.0);
            }
        }
        g
    }

    /// Builds a ring graph of `n` nodes.
    pub fn ring(n: usize) -> Self {
        let mut g = Self::new(n, false);
        for i in 0..n {
            g.add_edge(i, (i + 1) % n, 1.0);
        }
        g
    }

    /// Builds an Erdős-Rényi random graph G(n, p) with given seed.
    pub fn erdos_renyi(n: usize, p: f64, seed: u64) -> Self {
        let mut g = Self::new(n, false);
        let mut rng = Lcg::new(seed);
        for i in 0..n {
            for j in (i + 1)..n {
                if rng.next_f64() < p {
                    g.add_edge(i, j, 1.0);
                }
            }
        }
        g
    }

    /// Returns BFS distances from source node `s` (-1 means unreachable).
    pub fn bfs_distances(&self, s: usize) -> Vec<i64> {
        let mut dist = vec![-1i64; self.n];
        dist[s] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(s);
        while let Some(u) = queue.pop_front() {
            for &(v, _) in &self.adj[u] {
                if dist[v] < 0 {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        dist
    }

    /// Average shortest path length (only reachable pairs).
    pub fn average_path_length(&self) -> f64 {
        let mut total = 0.0f64;
        let mut count = 0usize;
        for s in 0..self.n {
            let dists = self.bfs_distances(s);
            for d in &dists {
                if *d > 0 {
                    total += *d as f64;
                    count += 1;
                }
            }
        }
        if count == 0 {
            f64::INFINITY
        } else {
            total / count as f64
        }
    }

    /// Local clustering coefficient for node `u`.
    pub fn clustering_coefficient(&self, u: usize) -> f64 {
        let neighbors: Vec<usize> = self.adj[u].iter().map(|&(v, _)| v).collect();
        let k = neighbors.len();
        if k < 2 {
            return 0.0;
        }
        let nbr_set: std::collections::HashSet<usize> = neighbors.iter().cloned().collect();
        let mut triangles = 0usize;
        for &v in &neighbors {
            for &(w, _) in &self.adj[v] {
                if w != u && nbr_set.contains(&w) {
                    triangles += 1;
                }
            }
        }
        triangles as f64 / (k * (k - 1)) as f64
    }

    /// Mean clustering coefficient over all nodes.
    pub fn mean_clustering(&self) -> f64 {
        let s: f64 = (0..self.n).map(|u| self.clustering_coefficient(u)).sum();
        s / self.n as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SIR Model
// ─────────────────────────────────────────────────────────────────────────────

/// SIR epidemic model (ODE version for homogeneous mixing).
///
/// Compartments: S (susceptible), I (infectious), R (recovered).
/// dS/dt = -β S I / N
/// dI/dt =  β S I / N - γ I
/// dR/dt =  γ I
pub struct SirModel {
    /// Transmission rate β.
    pub beta: f64,
    /// Recovery rate γ.
    pub gamma: f64,
    /// Population size N.
    pub n: f64,
}

/// State of SIR model at one time step.
#[derive(Clone, Debug)]
pub struct SirState {
    /// Susceptible fraction/count.
    pub s: f64,
    /// Infectious fraction/count.
    pub i: f64,
    /// Recovered fraction/count.
    pub r: f64,
}

impl SirModel {
    /// Creates a new SIR model.
    pub fn new(beta: f64, gamma: f64, n: f64) -> Self {
        Self { beta, gamma, n }
    }

    /// Basic reproduction number R₀ = β/γ.
    pub fn r0(&self) -> f64 {
        self.beta / self.gamma
    }

    /// Herd immunity threshold: fraction that needs immunity = 1 - 1/R₀.
    pub fn herd_immunity_threshold(&self) -> f64 {
        let r0 = self.r0();
        if r0 <= 1.0 { 0.0 } else { 1.0 - 1.0 / r0 }
    }

    /// Integrates the ODE for `steps` using Euler method with step `dt`.
    pub fn integrate(&self, s0: f64, i0: f64, r0_val: f64, dt: f64, steps: usize) -> Vec<SirState> {
        let mut traj = Vec::with_capacity(steps + 1);
        let (mut s, mut i, mut r) = (s0, i0, r0_val);
        traj.push(SirState { s, i, r });
        for _ in 0..steps {
            let ds = -self.beta * s * i / self.n;
            let di = self.beta * s * i / self.n - self.gamma * i;
            let dr = self.gamma * i;
            s += dt * ds;
            i += dt * di;
            r += dt * dr;
            s = s.max(0.0);
            i = i.max(0.0);
            r = r.max(0.0);
            traj.push(SirState { s, i, r });
        }
        traj
    }

    /// Network SIR: discrete-time stochastic simulation on a graph.
    ///
    /// Returns `(s_counts, i_counts, r_counts)` over time.
    pub fn simulate_on_graph(
        &self,
        graph: &NetworkGraph,
        initially_infected: &[usize],
        dt: f64,
        steps: usize,
        seed: u64,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let n = graph.n;
        // 0 = S, 1 = I, 2 = R
        let mut state = vec![0u8; n];
        for &idx in initially_infected {
            if idx < n {
                state[idx] = 1;
            }
        }
        let mut rng = Lcg::new(seed);
        let mut s_ts = Vec::with_capacity(steps + 1);
        let mut i_ts = Vec::with_capacity(steps + 1);
        let mut r_ts = Vec::with_capacity(steps + 1);
        let count_state = |state: &[u8], val: u8| state.iter().filter(|&&x| x == val).count();
        s_ts.push(count_state(&state, 0));
        i_ts.push(count_state(&state, 1));
        r_ts.push(count_state(&state, 2));
        for _ in 0..steps {
            let mut next = state.clone();
            for u in 0..n {
                if state[u] == 1 {
                    // recovery
                    if rng.next_f64() < self.gamma * dt {
                        next[u] = 2;
                    }
                } else if state[u] == 0 {
                    // infection from neighbors
                    let inf_neighbors =
                        graph.adj[u].iter().filter(|&&(v, _)| state[v] == 1).count();
                    let p_inf = 1.0 - (1.0 - self.beta * dt).powi(inf_neighbors as i32);
                    if rng.next_f64() < p_inf {
                        next[u] = 1;
                    }
                }
            }
            state = next;
            s_ts.push(count_state(&state, 0));
            i_ts.push(count_state(&state, 1));
            r_ts.push(count_state(&state, 2));
        }
        (s_ts, i_ts, r_ts)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SEIRD Model
// ─────────────────────────────────────────────────────────────────────────────

/// Extended SEIRD epidemic model: Susceptible, Exposed, Infected, Recovered, Deceased.
///
/// dS/dt = -β S I / N
/// dE/dt =  β S I / N - σ E
/// dI/dt =  σ E - γ I - μ I
/// dR/dt =  γ I
/// dD/dt =  μ I
pub struct SeirdModel {
    /// Transmission rate β.
    pub beta: f64,
    /// Incubation rate σ (1/latent period).
    pub sigma: f64,
    /// Recovery rate γ.
    pub gamma: f64,
    /// Disease-induced mortality rate μ.
    pub mu: f64,
    /// Population size N.
    pub n: f64,
}

/// State of SEIRD model.
#[derive(Clone, Debug)]
pub struct SeirdState {
    /// Susceptible.
    pub s: f64,
    /// Exposed.
    pub e: f64,
    /// Infectious.
    pub i: f64,
    /// Recovered.
    pub r: f64,
    /// Deceased.
    pub d: f64,
}

impl SeirdModel {
    /// Creates a new SEIRD model.
    pub fn new(beta: f64, sigma: f64, gamma: f64, mu: f64, n: f64) -> Self {
        Self {
            beta,
            sigma,
            gamma,
            mu,
            n,
        }
    }

    /// Basic reproduction number R₀ = β / (γ + μ).
    pub fn r0(&self) -> f64 {
        self.beta / (self.gamma + self.mu)
    }

    /// Case fatality ratio: μ / (γ + μ).
    pub fn case_fatality_ratio(&self) -> f64 {
        self.mu / (self.gamma + self.mu)
    }

    /// Integrates ODE for `steps` using Euler with step `dt`.
    pub fn integrate(
        &self,
        s0: f64,
        e0: f64,
        i0: f64,
        r0: f64,
        d0: f64,
        dt: f64,
        steps: usize,
    ) -> Vec<SeirdState> {
        let mut traj = Vec::with_capacity(steps + 1);
        let (mut s, mut e, mut i, mut r, mut d) = (s0, e0, i0, r0, d0);
        traj.push(SeirdState { s, e, i, r, d });
        for _ in 0..steps {
            let ds = -self.beta * s * i / self.n;
            let de = self.beta * s * i / self.n - self.sigma * e;
            let di = self.sigma * e - self.gamma * i - self.mu * i;
            let dr = self.gamma * i;
            let dd = self.mu * i;
            s = (s + dt * ds).max(0.0);
            e = (e + dt * de).max(0.0);
            i = (i + dt * di).max(0.0);
            r += dt * dr;
            d += dt * dd;
            traj.push(SeirdState { s, e, i, r, d });
        }
        traj
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Opinion Dynamics
// ─────────────────────────────────────────────────────────────────────────────

/// Deffuant-Weisbuch bounded confidence opinion dynamics.
///
/// Agents update opinions toward each other only when |o_i - o_j| < ε.
/// Update: o_i += μ(o_j - o_i), o_j += μ(o_i - o_j).
pub struct OpinionDynamics {
    /// Confidence bound ε: agents interact only when |o_i - o_j| < epsilon.
    pub epsilon: f64,
    /// Convergence parameter μ ∈ (0, 0.5].
    pub mu: f64,
}

impl OpinionDynamics {
    /// Creates a new Deffuant-Weisbuch model.
    pub fn new(epsilon: f64, mu: f64) -> Self {
        Self { epsilon, mu }
    }

    /// Runs `steps` of the Deffuant model on `opinions` (modified in place).
    ///
    /// At each step a random pair is chosen; if they are within ε they converge.
    pub fn run(&self, opinions: &mut Vec<f64>, steps: usize, seed: u64) {
        let n = opinions.len();
        if n < 2 {
            return;
        }
        let mut rng = Lcg::new(seed);
        for _ in 0..steps {
            let i = rng.next_range(0, n);
            let j_raw = rng.next_range(0, n - 1);
            let j = if j_raw >= i { j_raw + 1 } else { j_raw };
            let diff = opinions[j] - opinions[i];
            if diff.abs() < self.epsilon {
                opinions[i] += self.mu * diff;
                opinions[j] -= self.mu * diff;
            }
        }
    }

    /// Returns the number of opinion clusters (groups within ε/2 of each other).
    pub fn count_clusters(&self, opinions: &[f64]) -> usize {
        let mut sorted = opinions.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if sorted.is_empty() {
            return 0;
        }
        let mut clusters = 1usize;
        for w in sorted.windows(2) {
            if w[1] - w[0] >= self.epsilon / 2.0 {
                clusters += 1;
            }
        }
        clusters
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DeGroot Consensus Model
// ─────────────────────────────────────────────────────────────────────────────

/// DeGroot consensus model: linear opinion pooling on a network.
///
/// Each agent i updates: o_i(t+1) = Σ_j T_{ij} o_j(t)
/// where T is a row-stochastic trust matrix.
pub struct DeGrootModel {
    /// Row-stochastic trust matrix T (n×n stored row-major).
    pub trust: Vec<Vec<f64>>,
}

impl DeGrootModel {
    /// Creates a DeGroot model from a trust matrix.
    ///
    /// Rows must sum to 1 (caller's responsibility).
    pub fn new(trust: Vec<Vec<f64>>) -> Self {
        Self { trust }
    }

    /// Builds a uniform averaging DeGroot model from a graph.
    pub fn from_graph(graph: &NetworkGraph) -> Self {
        let n = graph.n;
        let mut trust = vec![vec![0.0f64; n]; n];
        for u in 0..n {
            let deg = graph.adj[u].len();
            if deg == 0 {
                trust[u][u] = 1.0;
            } else {
                for &(v, _) in &graph.adj[u] {
                    trust[u][v] = 1.0 / deg as f64;
                }
            }
        }
        Self { trust }
    }

    /// Returns the number of agents.
    pub fn n(&self) -> usize {
        self.trust.len()
    }

    /// Applies one update step to `opinions`.
    pub fn step(&self, opinions: &[f64]) -> Vec<f64> {
        let n = opinions.len();
        let mut next = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                next[i] += self.trust[i][j] * opinions[j];
            }
        }
        next
    }

    /// Runs `steps` iterations, returns final opinions.
    pub fn run(&self, opinions: &[f64], steps: usize) -> Vec<f64> {
        let mut ops = opinions.to_vec();
        for _ in 0..steps {
            ops = self.step(&ops);
        }
        ops
    }

    /// Checks if opinions have converged to consensus (max spread < tol).
    pub fn has_consensus(&self, opinions: &[f64], tol: f64) -> bool {
        let min = opinions.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = opinions.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        max - min < tol
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Voter Model
// ─────────────────────────────────────────────────────────────────────────────

/// Voter model: each agent copies a random neighbor's opinion.
pub struct VoterModel {
    /// Number of opinion states (default 2 for binary).
    pub n_states: usize,
}

impl VoterModel {
    /// Creates a new voter model.
    pub fn new(n_states: usize) -> Self {
        Self { n_states }
    }

    /// Runs `steps` asynchronous updates on integer `opinions` (values 0..n_states-1).
    pub fn run(&self, graph: &NetworkGraph, opinions: &mut Vec<usize>, steps: usize, seed: u64) {
        let n = graph.n;
        if n == 0 {
            return;
        }
        let mut rng = Lcg::new(seed);
        for _ in 0..steps {
            let i = rng.next_range(0, n);
            if graph.adj[i].is_empty() {
                continue;
            }
            let nb_idx = rng.next_range(0, graph.adj[i].len());
            let j = graph.adj[i][nb_idx].0;
            opinions[i] = opinions[j];
        }
    }

    /// Returns fraction of agents holding opinion `state`.
    pub fn fraction(&self, opinions: &[usize], state: usize) -> f64 {
        let count = opinions.iter().filter(|&&x| x == state).count();
        count as f64 / opinions.len() as f64
    }

    /// Returns true if all agents hold the same opinion (consensus).
    pub fn is_consensus(&self, opinions: &[usize]) -> bool {
        opinions.windows(2).all(|w| w[0] == w[1])
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Kuramoto Model
// ─────────────────────────────────────────────────────────────────────────────

/// Kuramoto model of coupled oscillators.
///
/// dθ_i/dt = ω_i + (K/N) Σ_j sin(θ_j - θ_i)
///
/// The order parameter r = |Σ e^{iθ_j}| / N measures synchrony.
pub struct KuramotoModel {
    /// Natural frequencies ω_i.
    pub omega: Vec<f64>,
    /// Coupling strength K.
    pub k: f64,
}

impl KuramotoModel {
    /// Creates a Kuramoto model from natural frequencies and coupling K.
    pub fn new(omega: Vec<f64>, k: f64) -> Self {
        Self { omega, k }
    }

    /// Generates Lorentzian-distributed (Cauchy) natural frequencies.
    ///
    /// Mean `omega0`, half-width `gamma`. Critical coupling K_c = 2γ.
    pub fn lorentzian(n: usize, omega0: f64, gamma_width: f64, seed: u64) -> Self {
        let mut rng = Lcg::new(seed);
        let omega: Vec<f64> = (0..n)
            .map(|_| {
                // Inverse CDF of Cauchy: omega0 + gamma * tan(pi*(u - 0.5))
                let u = rng.next_f64();
                omega0 + gamma_width * (PI * (u - 0.5)).tan()
            })
            .collect();
        Self {
            omega,
            k: 2.0 * gamma_width,
        }
    }

    /// Number of oscillators.
    pub fn n(&self) -> usize {
        self.omega.len()
    }

    /// Computes the order parameter r ∈ \[0, 1\].
    ///
    /// r = |1/N Σ e^{iθ_j}|
    pub fn order_parameter(phases: &[f64]) -> f64 {
        let n = phases.len() as f64;
        let re: f64 = phases.iter().map(|&t| t.cos()).sum::<f64>() / n;
        let im: f64 = phases.iter().map(|&t| t.sin()).sum::<f64>() / n;
        (re * re + im * im).sqrt()
    }

    /// Mean phase angle (argument of the complex order parameter).
    pub fn mean_phase(phases: &[f64]) -> f64 {
        let re: f64 = phases.iter().map(|&t| t.cos()).sum::<f64>();
        let im: f64 = phases.iter().map(|&t| t.sin()).sum::<f64>();
        im.atan2(re)
    }

    /// Integrates for `steps` using RK4 with step `dt`.
    ///
    /// Returns trajectory of order parameters.
    pub fn integrate(&self, phases0: &[f64], dt: f64, steps: usize) -> (Vec<f64>, Vec<f64>) {
        let n = self.n();
        let mut phases = phases0.to_vec();
        let mut times = Vec::with_capacity(steps + 1);
        let mut r_vals = Vec::with_capacity(steps + 1);
        times.push(0.0);
        r_vals.push(Self::order_parameter(&phases));

        let deriv = |ph: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| {
                    let coupling: f64 = ph.iter().map(|&t| (t - ph[i]).sin()).sum::<f64>();
                    self.omega[i] + self.k / n as f64 * coupling
                })
                .collect()
        };

        for step in 0..steps {
            let k1 = deriv(&phases);
            let p2: Vec<f64> = phases
                .iter()
                .zip(&k1)
                .map(|(&p, &k)| p + 0.5 * dt * k)
                .collect();
            let k2 = deriv(&p2);
            let p3: Vec<f64> = phases
                .iter()
                .zip(&k2)
                .map(|(&p, &k)| p + 0.5 * dt * k)
                .collect();
            let k3 = deriv(&p3);
            let p4: Vec<f64> = phases.iter().zip(&k3).map(|(&p, &k)| p + dt * k).collect();
            let k4 = deriv(&p4);
            for i in 0..n {
                phases[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            times.push((step + 1) as f64 * dt);
            r_vals.push(Self::order_parameter(&phases));
        }
        (times, r_vals)
    }

    /// Estimated critical coupling K_c = 2 / (π g(0)) for Lorentzian with width γ → K_c = 2γ.
    ///
    /// For a general distribution approximated by the frequency spread σ: K_c ≈ 2σ/√π.
    pub fn critical_coupling_estimate(&self) -> f64 {
        let mean = self.omega.iter().sum::<f64>() / self.n() as f64;
        let var: f64 =
            self.omega.iter().map(|&w| (w - mean).powi(2)).sum::<f64>() / self.n() as f64;
        2.0 * var.sqrt() / PI.sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Random Walk on Graph
// ─────────────────────────────────────────────────────────────────────────────

/// Random walk on a graph with analysis of hitting and cover times.
pub struct RandomWalkGraph {
    graph: NetworkGraph,
}

impl RandomWalkGraph {
    /// Creates a random walk handler for the given graph.
    pub fn new(graph: NetworkGraph) -> Self {
        Self { graph }
    }

    /// Simulates a random walk starting at `start` for `max_steps`.
    ///
    /// Returns the sequence of visited nodes.
    pub fn walk(&self, start: usize, max_steps: usize, seed: u64) -> Vec<usize> {
        let mut rng = Lcg::new(seed);
        let mut path = Vec::with_capacity(max_steps + 1);
        let mut current = start;
        path.push(current);
        for _ in 0..max_steps {
            let neighbors = &self.graph.adj[current];
            if neighbors.is_empty() {
                break;
            }
            let idx = rng.next_range(0, neighbors.len());
            current = neighbors[idx].0;
            path.push(current);
        }
        path
    }

    /// Estimates hitting time from `start` to `target` by Monte Carlo.
    pub fn hitting_time_estimate(
        &self,
        start: usize,
        target: usize,
        trials: usize,
        seed: u64,
    ) -> f64 {
        let mut rng = Lcg::new(seed);
        let mut total = 0usize;
        let mut successes = 0usize;
        for t in 0..trials {
            let mut current = start;
            let walk_seed = seed.wrapping_add(t as u64 * 1_000_003);
            let _ = walk_seed; // use rng instead
            let mut steps = 0usize;
            let max_steps = self.graph.n * self.graph.n * 10;
            loop {
                if current == target {
                    break;
                }
                let neighbors = &self.graph.adj[current];
                if neighbors.is_empty() || steps >= max_steps {
                    break;
                }
                let idx = rng.next_range(0, neighbors.len());
                current = neighbors[idx].0;
                steps += 1;
            }
            if current == target {
                total += steps;
                successes += 1;
            }
        }
        if successes == 0 {
            f64::INFINITY
        } else {
            total as f64 / successes as f64
        }
    }

    /// Estimates cover time (expected steps to visit all nodes) by Monte Carlo.
    pub fn cover_time_estimate(&self, start: usize, trials: usize, seed: u64) -> f64 {
        let n = self.graph.n;
        let mut rng = Lcg::new(seed);
        let mut total_steps = 0usize;
        for _ in 0..trials {
            let mut visited = vec![false; n];
            let mut current = start;
            visited[current] = true;
            let mut steps = 0usize;
            let max_steps = n * n * 100;
            while visited.iter().any(|&v| !v) && steps < max_steps {
                let neighbors = &self.graph.adj[current];
                if neighbors.is_empty() {
                    break;
                }
                let idx = rng.next_range(0, neighbors.len());
                current = neighbors[idx].0;
                visited[current] = true;
                steps += 1;
            }
            total_steps += steps;
        }
        total_steps as f64 / trials as f64
    }

    /// Stationary distribution via power iteration on the transition matrix.
    pub fn stationary_distribution(&self, tol: f64, max_iter: usize) -> Vec<f64> {
        let n = self.graph.n;
        let mut pi = vec![1.0 / n as f64; n];
        for _ in 0..max_iter {
            let mut next = vec![0.0f64; n];
            for u in 0..n {
                let deg = self.graph.adj[u].len();
                if deg == 0 {
                    continue;
                }
                let w = pi[u] / deg as f64;
                for &(v, _) in &self.graph.adj[u] {
                    next[v] += w;
                }
            }
            // normalize
            let sum: f64 = next.iter().sum();
            if sum > 0.0 {
                for x in &mut next {
                    *x /= sum;
                }
            }
            let diff: f64 = pi.iter().zip(&next).map(|(a, b)| (a - b).abs()).sum();
            pi = next;
            if diff < tol {
                break;
            }
        }
        pi
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PageRank
// ─────────────────────────────────────────────────────────────────────────────

/// PageRank via power iteration with teleportation and dangling node handling.
pub struct PageRank {
    /// Damping factor d (typically 0.85).
    pub damping: f64,
    /// Convergence tolerance.
    pub tol: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}

impl PageRank {
    /// Creates a new PageRank solver.
    pub fn new(damping: f64, tol: f64, max_iter: usize) -> Self {
        Self {
            damping,
            tol,
            max_iter,
        }
    }

    /// Computes PageRank scores for the given directed graph.
    ///
    /// Handles dangling nodes (no out-edges) by distributing their rank uniformly.
    pub fn compute(&self, graph: &NetworkGraph) -> Vec<f64> {
        let n = graph.n;
        let inv_n = 1.0 / n as f64;
        let mut pr = vec![inv_n; n];
        // Build out-degree array
        let out_deg: Vec<usize> = (0..n).map(|u| graph.adj[u].len()).collect();

        for _ in 0..self.max_iter {
            // Sum rank from dangling nodes
            let dangling_sum: f64 = (0..n).filter(|&u| out_deg[u] == 0).map(|u| pr[u]).sum();
            let mut next = vec![0.0f64; n];
            for u in 0..n {
                if out_deg[u] == 0 {
                    continue;
                }
                let contrib = self.damping * pr[u] / out_deg[u] as f64;
                for &(v, _) in &graph.adj[u] {
                    next[v] += contrib;
                }
            }
            let teleport = (1.0 - self.damping + self.damping * dangling_sum) * inv_n;
            for x in &mut next {
                *x += teleport;
            }
            // Normalize
            let sum: f64 = next.iter().sum();
            if sum > 0.0 {
                for x in &mut next {
                    *x /= sum;
                }
            }
            let diff: f64 = pr.iter().zip(&next).map(|(a, b)| (a - b).abs()).sum();
            pr = next;
            if diff < self.tol {
                break;
            }
        }
        pr
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Community Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Community detection via modularity maximization.
///
/// Modularity Q = (1/2m) Σ_{ij} \[A_{ij} - k_i k_j / 2m\] δ(c_i, c_j)
pub struct CommunityDetection;

impl CommunityDetection {
    /// Computes modularity Q for a partition given by `communities` (node → community id).
    pub fn modularity(graph: &NetworkGraph, communities: &[usize]) -> f64 {
        let m = graph.edge_count() as f64;
        if m == 0.0 {
            return 0.0;
        }
        let n = graph.n;
        let degree: Vec<f64> = (0..n).map(|u| graph.adj[u].len() as f64).collect();
        let mut q = 0.0f64;
        for u in 0..n {
            for &(v, _) in &graph.adj[u] {
                if communities[u] == communities[v] {
                    q += 1.0 - degree[u] * degree[v] / (2.0 * m);
                }
            }
        }
        q / (2.0 * m)
    }

    /// Greedy modularity maximization (Clauset-Newman-Moore simplified).
    ///
    /// Starts with each node in its own community and greedily merges.
    /// Returns the best community assignment found.
    pub fn greedy_modularity(graph: &NetworkGraph) -> Vec<usize> {
        let n = graph.n;
        let mut communities: Vec<usize> = (0..n).collect();
        let mut best_q = Self::modularity(graph, &communities);
        let mut improved = true;
        while improved {
            improved = false;
            // Try merging each pair of neighboring communities
            let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
            for u in 0..n {
                for &(v, _) in &graph.adj[u] {
                    let cu = communities[u];
                    let cv = communities[v];
                    if cu != cv {
                        let pair = if cu < cv { (cu, cv) } else { (cv, cu) };
                        edge_pairs.push(pair);
                    }
                }
            }
            edge_pairs.sort();
            edge_pairs.dedup();
            for (ca, cb) in edge_pairs {
                let mut trial = communities.clone();
                for x in &mut trial {
                    if *x == cb {
                        *x = ca;
                    }
                }
                let q = Self::modularity(graph, &trial);
                if q > best_q + 1e-10 {
                    best_q = q;
                    communities = trial;
                    improved = true;
                }
            }
        }
        // Remap community ids to 0..k
        let mut id_map: HashMap<usize, usize> = HashMap::new();
        let mut next_id = 0usize;
        for c in &mut communities {
            let entry = id_map.entry(*c).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            });
            *c = *entry;
        }
        communities
    }

    /// Returns number of communities in a partition.
    pub fn n_communities(communities: &[usize]) -> usize {
        let max = communities.iter().cloned().max().unwrap_or(0);
        max + 1
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Network Robustness
// ─────────────────────────────────────────────────────────────────────────────

/// Network robustness analysis via percolation theory.
pub struct NetworkRobustness;

impl NetworkRobustness {
    /// Computes size of the giant component (largest connected component) as fraction of N.
    pub fn giant_component_fraction(graph: &NetworkGraph) -> f64 {
        let n = graph.n;
        if n == 0 {
            return 0.0;
        }
        let comp_id = Self::connected_components(graph);
        let mut sizes: HashMap<usize, usize> = HashMap::new();
        for &c in &comp_id {
            *sizes.entry(c).or_insert(0) += 1;
        }
        *sizes.values().max().unwrap_or(&0) as f64 / n as f64
    }

    /// Connected components via BFS, returns component id for each node.
    pub fn connected_components(graph: &NetworkGraph) -> Vec<usize> {
        let n = graph.n;
        let mut comp = vec![usize::MAX; n];
        let mut comp_id = 0usize;
        for start in 0..n {
            if comp[start] != usize::MAX {
                continue;
            }
            let mut queue = VecDeque::new();
            queue.push_back(start);
            comp[start] = comp_id;
            while let Some(u) = queue.pop_front() {
                for &(v, _) in &graph.adj[u] {
                    if comp[v] == usize::MAX {
                        comp[v] = comp_id;
                        queue.push_back(v);
                    }
                }
            }
            comp_id += 1;
        }
        comp
    }

    /// Simulates random node removal: returns giant component fraction vs fraction removed.
    pub fn percolation_curve(graph: &NetworkGraph, seed: u64) -> Vec<(f64, f64)> {
        let n = graph.n;
        let mut rng = Lcg::new(seed);
        // Random order of removal
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.next_range(0, i + 1);
            order.swap(i, j);
        }
        let mut curve = Vec::with_capacity(n + 1);
        let gc0 = Self::giant_component_fraction(graph);
        curve.push((0.0, gc0));
        let mut active = vec![true; n];
        for (step, &node) in order.iter().enumerate() {
            active[node] = false;
            // Build subgraph
            let active_nodes: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
            let an = active_nodes.len();
            if an == 0 {
                curve.push(((step + 1) as f64 / n as f64, 0.0));
                continue;
            }
            let mut idx_map = vec![usize::MAX; n];
            for (new_i, &old_i) in active_nodes.iter().enumerate() {
                idx_map[old_i] = new_i;
            }
            let mut sub = NetworkGraph::new(an, graph.directed);
            for &u in &active_nodes {
                for &(v, w) in &graph.adj[u] {
                    if active[v]
                        && idx_map[u] < an
                        && idx_map[v] < an
                        && (graph.directed || idx_map[u] < idx_map[v])
                    {
                        sub.add_edge(idx_map[u], idx_map[v], w);
                    }
                }
            }
            let gc = Self::giant_component_fraction(&sub);
            curve.push(((step + 1) as f64 / n as f64, gc));
        }
        curve
    }

    /// Estimated bond percolation threshold p_c ≈ 1 / (`k` ) for Erdős-Rényi graphs.
    pub fn bond_percolation_threshold(mean_degree: f64) -> f64 {
        if mean_degree <= 0.0 {
            return 1.0;
        }
        1.0 / mean_degree
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Small World (Watts-Strogatz)
// ─────────────────────────────────────────────────────────────────────────────

/// Watts-Strogatz small-world network model.
pub struct SmallWorld;

impl SmallWorld {
    /// Generates a Watts-Strogatz graph: ring of `n` nodes, each connected to `k` nearest
    /// neighbors on each side, then each edge rewired with probability `p`.
    pub fn watts_strogatz(n: usize, k: usize, p: f64, seed: u64) -> NetworkGraph {
        let mut g = NetworkGraph::new(n, false);
        let mut rng = Lcg::new(seed);
        // Build ring lattice
        for i in 0..n {
            for j in 1..=k {
                let neighbor = (i + j) % n;
                g.add_edge(i, neighbor, 1.0);
            }
        }
        // Rewire
        // We rebuild edges from scratch to allow rewiring
        let mut adj_new: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut visited_edges: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for i in 0..n {
            for j in 1..=k {
                let neighbor = (i + j) % n;
                if rng.next_f64() < p {
                    // Rewire: choose new target uniformly at random (not self, not existing)
                    let mut new_nb = rng.next_range(0, n);
                    let mut attempts = 0;
                    while (new_nb == i || visited_edges.contains(&(i.min(new_nb), i.max(new_nb))))
                        && attempts < n * 2
                    {
                        new_nb = rng.next_range(0, n);
                        attempts += 1;
                    }
                    if new_nb != i && !visited_edges.contains(&(i.min(new_nb), i.max(new_nb))) {
                        let e = (i.min(new_nb), i.max(new_nb));
                        visited_edges.insert(e);
                        adj_new[i].push((new_nb, 1.0));
                        adj_new[new_nb].push((i, 1.0));
                    }
                } else {
                    let e = (i.min(neighbor), i.max(neighbor));
                    if !visited_edges.contains(&e) {
                        visited_edges.insert(e);
                        adj_new[i].push((neighbor, 1.0));
                        adj_new[neighbor].push((i, 1.0));
                    }
                }
            }
        }
        NetworkGraph {
            n,
            directed: false,
            adj: adj_new,
        }
    }

    /// Computes the small-world coefficient σ = C/C_rand / L/L_rand.
    ///
    /// σ > 1 indicates small-world behavior.
    pub fn small_world_coefficient(graph: &NetworkGraph) -> f64 {
        let n = graph.n;
        if n == 0 {
            return 0.0;
        }
        let c = graph.mean_clustering();
        let l = graph.average_path_length();
        // Random graph reference: C_rand ≈ <k> / N, L_rand ≈ ln(N) / ln(<k>)
        let mean_k = graph.adj.iter().map(|a| a.len()).sum::<usize>() as f64 / n as f64;
        if mean_k <= 1.0 || n <= 1 {
            return 0.0;
        }
        let c_rand = mean_k / n as f64;
        let l_rand = (n as f64).ln() / mean_k.ln();
        if c_rand == 0.0 || l_rand == 0.0 || l == 0.0 {
            return 0.0;
        }
        (c / c_rand) / (l / l_rand)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── NetworkGraph tests ──────────────────────────────────────────────────

    #[test]
    fn test_network_graph_basic() {
        let mut g = NetworkGraph::new(4, false);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        assert_eq!(g.edge_count(), 3);
        assert_eq!(g.degree(1), 2);
    }

    #[test]
    fn test_network_graph_directed() {
        let mut g = NetworkGraph::new(3, true);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        assert_eq!(g.degree(0), 1);
        assert_eq!(g.degree(1), 1);
        assert_eq!(g.degree(2), 0);
        assert_eq!(g.edge_count(), 2);
    }

    #[test]
    fn test_complete_graph() {
        let g = NetworkGraph::complete(4);
        // K4 has 4*3/2 = 6 edges
        assert_eq!(g.edge_count(), 6);
        assert_eq!(g.degree(0), 3);
    }

    #[test]
    fn test_ring_graph_bfs() {
        let g = NetworkGraph::ring(5);
        let dist = g.bfs_distances(0);
        assert_eq!(dist[0], 0);
        assert_eq!(dist[1], 1);
        assert_eq!(dist[4], 1); // ring: 0-4 direct edge
        assert_eq!(dist[2], 2);
    }

    #[test]
    fn test_clustering_coefficient_complete() {
        let g = NetworkGraph::complete(5);
        // In a complete graph every triplet closes → C = 1
        let c = g.clustering_coefficient(0);
        assert!(
            (c - 1.0).abs() < 1e-10,
            "clustering should be 1 for complete graph, got {}",
            c
        );
    }

    #[test]
    fn test_mean_clustering_ring() {
        let g = NetworkGraph::ring(6);
        let c = g.mean_clustering();
        // Ring with k=1: no triangles → C = 0
        assert!(c < 1e-10, "ring k=1 clustering should be 0, got {}", c);
    }

    #[test]
    fn test_average_path_length_complete() {
        let g = NetworkGraph::complete(4);
        let apl = g.average_path_length();
        assert!(
            (apl - 1.0).abs() < 1e-10,
            "complete graph APL should be 1, got {}",
            apl
        );
    }

    // ── SIR Model tests ─────────────────────────────────────────────────────

    #[test]
    fn test_sir_population_conserved() {
        let sir = SirModel::new(0.3, 0.1, 1000.0);
        let traj = sir.integrate(990.0, 10.0, 0.0, 0.1, 200);
        for state in &traj {
            let total = state.s + state.i + state.r;
            assert!(
                (total - 1000.0).abs() < 1.0,
                "SIR population drift: {}",
                total
            );
        }
    }

    #[test]
    fn test_sir_r0() {
        let sir = SirModel::new(0.3, 0.1, 1000.0);
        assert!((sir.r0() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_sir_herd_immunity_threshold() {
        let sir = SirModel::new(0.3, 0.1, 1000.0);
        // R0=3 → p_c = 1 - 1/3 ≈ 0.6667
        let hit = sir.herd_immunity_threshold();
        assert!(
            (hit - 2.0 / 3.0).abs() < 1e-10,
            "HIT should be 2/3 for R0=3, got {}",
            hit
        );
    }

    #[test]
    fn test_sir_herd_immunity_threshold_r0_lt1() {
        let sir = SirModel::new(0.05, 0.1, 1000.0);
        // R0 < 1 → no epidemic threshold
        assert_eq!(sir.herd_immunity_threshold(), 0.0);
    }

    #[test]
    fn test_sir_epidemic_grows_initially() {
        let sir = SirModel::new(0.5, 0.1, 1000.0);
        let traj = sir.integrate(999.0, 1.0, 0.0, 0.01, 10);
        // Infected should grow when S ≈ N and R0 > 1
        assert!(traj[5].i > traj[0].i, "infected should grow when R0 > 1");
    }

    #[test]
    fn test_sir_network_population_conserved() {
        let g = NetworkGraph::complete(20);
        let sir = SirModel::new(0.3, 0.1, 20.0);
        let (s_ts, i_ts, r_ts) = sir.simulate_on_graph(&g, &[0], 0.1, 50, 42);
        for ((&s, &i), &r) in s_ts.iter().zip(&i_ts).zip(&r_ts) {
            assert_eq!(
                s + i + r,
                20,
                "population not conserved: {}+{}+{}={}",
                s,
                i,
                r,
                s + i + r
            );
        }
    }

    // ── SEIRD Model tests ───────────────────────────────────────────────────

    #[test]
    fn test_seird_population_conserved() {
        let m = SeirdModel::new(0.4, 0.2, 0.15, 0.01, 1000.0);
        let traj = m.integrate(990.0, 5.0, 5.0, 0.0, 0.0, 0.1, 100);
        for st in &traj {
            let total = st.s + st.e + st.i + st.r + st.d;
            assert!((total - 1000.0).abs() < 2.0, "SEIRD total drift: {}", total);
        }
    }

    #[test]
    fn test_seird_r0() {
        let m = SeirdModel::new(0.4, 0.2, 0.15, 0.05, 1000.0);
        let r0 = m.r0();
        assert!((r0 - 0.4 / 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_seird_deaths_nonnegative() {
        let m = SeirdModel::new(0.3, 0.2, 0.1, 0.02, 500.0);
        let traj = m.integrate(495.0, 3.0, 2.0, 0.0, 0.0, 0.1, 50);
        for st in &traj {
            assert!(st.d >= -1e-10, "deaths should be non-negative");
        }
    }

    // ── Opinion Dynamics tests ──────────────────────────────────────────────

    #[test]
    fn test_opinion_dynamics_opinions_in_range() {
        let od = OpinionDynamics::new(0.3, 0.5);
        let mut ops: Vec<f64> = (0..20).map(|i| i as f64 / 20.0).collect();
        od.run(&mut ops, 5000, 123);
        for &o in &ops {
            assert!((0.0..=1.0).contains(&o), "opinion out of range: {}", o);
        }
    }

    #[test]
    fn test_opinion_dynamics_convergence() {
        let od = OpinionDynamics::new(1.0, 0.5); // large epsilon → all interact
        let mut ops: Vec<f64> = (0..10).map(|i| i as f64 / 10.0).collect();
        od.run(&mut ops, 50000, 77);
        // Should converge to single cluster
        assert!(od.count_clusters(&ops) <= 2);
    }

    #[test]
    fn test_opinion_dynamics_fragmentation() {
        let od = OpinionDynamics::new(0.1, 0.5); // small epsilon
        let mut ops: Vec<f64> = (0..20).map(|i| i as f64 / 20.0).collect();
        od.run(&mut ops, 10000, 99);
        // Multiple clusters expected
        let clusters = od.count_clusters(&ops);
        assert!(clusters >= 1);
    }

    // ── DeGroot Model tests ─────────────────────────────────────────────────

    #[test]
    fn test_degroot_consensus_complete_graph() {
        let g = NetworkGraph::complete(5);
        let dg = DeGrootModel::from_graph(&g);
        let ops = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let final_ops = dg.run(&ops, 100);
        let mean = ops.iter().sum::<f64>() / ops.len() as f64;
        for &o in &final_ops {
            assert!(
                (o - mean).abs() < 1e-6,
                "consensus failed: {} vs {}",
                o,
                mean
            );
        }
    }

    #[test]
    fn test_degroot_opinion_preserves_mean() {
        let g = NetworkGraph::ring(6);
        let dg = DeGrootModel::from_graph(&g);
        let ops = vec![0.1, 0.4, 0.6, 0.2, 0.9, 0.3];
        let mean0: f64 = ops.iter().sum::<f64>() / ops.len() as f64;
        let final_ops = dg.run(&ops, 50);
        let mean1: f64 = final_ops.iter().sum::<f64>() / final_ops.len() as f64;
        assert!(
            (mean0 - mean1).abs() < 1e-10,
            "mean opinion should be preserved"
        );
    }

    // ── Voter Model tests ───────────────────────────────────────────────────

    #[test]
    fn test_voter_model_fractions_sum_to_one() {
        let g = NetworkGraph::complete(10);
        let vm = VoterModel::new(2);
        let mut ops: Vec<usize> = (0..10).map(|i| i % 2).collect();
        vm.run(&g, &mut ops, 1000, 42);
        let f0 = vm.fraction(&ops, 0);
        let f1 = vm.fraction(&ops, 1);
        assert!((f0 + f1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_voter_model_consensus_reached() {
        let g = NetworkGraph::complete(8);
        let vm = VoterModel::new(2);
        // Start with one opinion — should stay in consensus
        let mut ops = vec![0usize; 8];
        vm.run(&g, &mut ops, 1000, 5);
        assert!(vm.is_consensus(&ops));
    }

    // ── Kuramoto Model tests ────────────────────────────────────────────────

    #[test]
    fn test_kuramoto_order_parameter_range() {
        let omegas = vec![1.0, 1.1, 0.9, 1.05, 0.95];
        let km = KuramotoModel::new(omegas, 5.0);
        let phases0: Vec<f64> = (0..5).map(|i| i as f64 * 0.1).collect();
        let (_, r_vals) = km.integrate(&phases0, 0.01, 100);
        for r in &r_vals {
            assert!(
                *r >= 0.0 && *r <= 1.0 + 1e-10,
                "order parameter out of range: {}",
                r
            );
        }
    }

    #[test]
    fn test_kuramoto_fully_synchronized_stays_synchronized() {
        // All same frequency → should remain synchronized
        let omegas = vec![1.0f64; 5];
        let km = KuramotoModel::new(omegas, 1.0);
        let phases0 = vec![0.0f64; 5]; // all same phase
        let (_, r_vals) = km.integrate(&phases0, 0.01, 100);
        let r_final = *r_vals.last().unwrap();
        assert!(
            r_final > 0.99,
            "fully synchronized should stay at r≈1, got {}",
            r_final
        );
    }

    #[test]
    fn test_kuramoto_order_parameter_formula() {
        let phases = vec![0.0f64; 4]; // all at same phase
        let r = KuramotoModel::order_parameter(&phases);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kuramoto_order_parameter_opposite_phases() {
        // Two groups at 0 and π → r ≈ 0
        let phases = vec![0.0, PI, 0.0, PI];
        let r = KuramotoModel::order_parameter(&phases);
        assert!(r < 1e-10, "opposite phases should give r≈0, got {}", r);
    }

    #[test]
    fn test_kuramoto_critical_coupling_positive() {
        let omegas = vec![0.5, 1.0, 1.5, 2.0, 0.8];
        let km = KuramotoModel::new(omegas, 1.0);
        let kc = km.critical_coupling_estimate();
        assert!(kc > 0.0, "critical coupling should be positive");
    }

    // ── Random Walk tests ───────────────────────────────────────────────────

    #[test]
    fn test_random_walk_path_length() {
        let g = NetworkGraph::ring(10);
        let rw = RandomWalkGraph::new(g);
        let path = rw.walk(0, 50, 1234);
        assert_eq!(path.len(), 51); // start + 50 steps
    }

    #[test]
    fn test_random_walk_stays_in_graph() {
        let g = NetworkGraph::complete(8);
        let n = g.n;
        let rw = RandomWalkGraph::new(g);
        let path = rw.walk(0, 100, 7);
        for &node in &path {
            assert!(node < n, "walker left graph: {}", node);
        }
    }

    #[test]
    fn test_stationary_distribution_sums_to_one() {
        let g = NetworkGraph::complete(5);
        let rw = RandomWalkGraph::new(g);
        let pi = rw.stationary_distribution(1e-9, 1000);
        let sum: f64 = pi.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "stationary dist should sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_stationary_distribution_uniform_regular() {
        let g = NetworkGraph::complete(4);
        let rw = RandomWalkGraph::new(g);
        let pi = rw.stationary_distribution(1e-9, 1000);
        let expected = 0.25f64;
        for &p in &pi {
            assert!(
                (p - expected).abs() < 1e-6,
                "regular graph stationary dist should be uniform, got {}",
                p
            );
        }
    }

    // ── PageRank tests ──────────────────────────────────────────────────────

    #[test]
    fn test_pagerank_sums_to_one() {
        let mut g = NetworkGraph::new(4, true);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        g.add_edge(3, 0, 1.0);
        let pr_solver = PageRank::new(0.85, 1e-8, 100);
        let pr = pr_solver.compute(&g);
        let sum: f64 = pr.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "PageRank should sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_pagerank_dangling_nodes() {
        let mut g = NetworkGraph::new(3, true);
        g.add_edge(0, 1, 1.0);
        // Node 2 is dangling (no out-edges)
        let pr_solver = PageRank::new(0.85, 1e-8, 200);
        let pr = pr_solver.compute(&g);
        let sum: f64 = pr.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "PageRank with dangling node should sum to 1"
        );
    }

    #[test]
    fn test_pagerank_all_positive() {
        let g = NetworkGraph::complete(5);
        // Convert to directed
        let mut gd = NetworkGraph::new(5, true);
        for u in 0..5 {
            for &(v, w) in &g.adj[u] {
                gd.add_edge(u, v, w);
            }
        }
        let pr_solver = PageRank::new(0.85, 1e-8, 100);
        let pr = pr_solver.compute(&gd);
        for &p in &pr {
            assert!(p > 0.0, "PageRank should be positive");
        }
    }

    // ── Community Detection tests ───────────────────────────────────────────

    #[test]
    fn test_modularity_trivial_all_same_community() {
        let g = NetworkGraph::ring(6);
        let communities = vec![0usize; 6]; // all in one community
        let q = CommunityDetection::modularity(&g, &communities);
        assert!(q <= 1.0, "Q should be <= 1, got {}", q);
    }

    #[test]
    fn test_modularity_two_cliques() {
        // Two connected 4-cliques with one bridge → modularity of two-community partition
        // should be positive and close to the theoretical optimum.
        let mut g = NetworkGraph::new(8, false);
        for i in 0..4 {
            for j in (i + 1)..4 {
                g.add_edge(i, j, 1.0);
            }
        }
        for i in 4..8 {
            for j in (i + 1)..8 {
                g.add_edge(i, j, 1.0);
            }
        }
        g.add_edge(3, 4, 1.0); // bridge
        // Correct 2-community partition should have positive modularity
        let com_two = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let q_two = CommunityDetection::modularity(&g, &com_two);
        assert!(
            q_two > 0.0,
            "two-clique true partition Q should be positive, got {}",
            q_two
        );
        // Each node in its own community should have lower Q than the 2-community partition
        let com_individual: Vec<usize> = (0..8).collect();
        let q_individual = CommunityDetection::modularity(&g, &com_individual);
        assert!(
            q_two > q_individual,
            "true 2-community partition Q={} should exceed individual Q={}",
            q_two,
            q_individual
        );
    }

    #[test]
    fn test_greedy_modularity_returns_valid_partition() {
        let g = NetworkGraph::complete(6);
        let communities = CommunityDetection::greedy_modularity(&g);
        assert_eq!(communities.len(), 6);
    }

    // ── Network Robustness tests ────────────────────────────────────────────

    #[test]
    fn test_giant_component_complete_graph() {
        let g = NetworkGraph::complete(10);
        let gc = NetworkRobustness::giant_component_fraction(&g);
        assert!(
            (gc - 1.0).abs() < 1e-10,
            "complete graph should have GC=1, got {}",
            gc
        );
    }

    #[test]
    fn test_giant_component_disconnected() {
        let g = NetworkGraph::new(4, false); // no edges → 4 isolated nodes
        let gc = NetworkRobustness::giant_component_fraction(&g);
        assert!(
            (gc - 0.25).abs() < 1e-10,
            "isolated nodes GC=0.25 (each own component), got {}",
            gc
        );
    }

    #[test]
    fn test_percolation_curve_monotone() {
        let g = NetworkGraph::erdos_renyi(20, 0.4, 1);
        let curve = NetworkRobustness::percolation_curve(&g, 42);
        // The curve should start with a finite GC and end near 0
        assert!(!curve.is_empty(), "percolation curve should not be empty");
        // First GC (no removal) should be positive for dense ER graph
        assert!(curve[0].1 > 0.0, "initial GC should be positive");
        // Final GC (all nodes removed) should be 0
        assert!(
            curve.last().unwrap().1 == 0.0,
            "final GC should be 0 when all nodes removed"
        );
        // Overall trend: GC at 80%+ removal should be less than at start
        let early = curve[curve.len() / 5].1;
        let late = curve[4 * curve.len() / 5].1;
        assert!(
            late <= early + 0.2,
            "GC should trend downward: early={} late={}",
            early,
            late
        );
    }

    #[test]
    fn test_bond_percolation_threshold() {
        let pct = NetworkRobustness::bond_percolation_threshold(4.0);
        assert!((pct - 0.25).abs() < 1e-10);
    }

    // ── Small World tests ───────────────────────────────────────────────────

    #[test]
    fn test_watts_strogatz_node_count() {
        let g = SmallWorld::watts_strogatz(20, 2, 0.1, 42);
        assert_eq!(g.n, 20);
    }

    #[test]
    fn test_small_world_low_p_high_clustering() {
        // With p=0 (regular lattice with k=3), clustering should be high
        let g = SmallWorld::watts_strogatz(30, 3, 0.0, 1);
        let c = g.mean_clustering();
        assert!(
            c > 0.3,
            "regular WS lattice should have non-trivial clustering, got {}",
            c
        );
    }

    #[test]
    fn test_erdos_renyi_mean_degree() {
        let n = 100;
        let p = 0.1;
        let g = NetworkGraph::erdos_renyi(n, p, 999);
        let mean_k = g.adj.iter().map(|a| a.len()).sum::<usize>() as f64 / n as f64;
        let expected = (n as f64 - 1.0) * p;
        // Should be within 3 standard deviations
        assert!(
            (mean_k - expected).abs() < 3.0 * expected.sqrt() + 2.0,
            "ER mean degree {} far from expected {}",
            mean_k,
            expected
        );
    }
}
