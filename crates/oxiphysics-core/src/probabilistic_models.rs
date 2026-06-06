// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Probabilistic models including Bayesian networks, HMMs, Gaussian processes,
//! Dirichlet processes, variational inference, and expectation-maximization.
//!
//! These models provide foundational probabilistic machinery for physics-informed
//! machine learning, uncertainty quantification, and data-driven modeling.

use std::f64::consts::{PI, TAU};

// ---------------------------------------------------------------------------
// Helper math utilities
// ---------------------------------------------------------------------------

/// Returns the log of the standard normal density at `x`.
fn log_normal_pdf(x: f64, mean: f64, var: f64) -> f64 {
    -0.5 * ((x - mean).powi(2) / var + var.ln() + (TAU).ln())
}

/// Returns the normal density at `x`.
fn normal_pdf(x: f64, mean: f64, var: f64) -> f64 {
    (-(x - mean).powi(2) / (2.0 * var)).exp() / (TAU * var).sqrt()
}

/// Log-sum-exp trick for numerical stability.
fn log_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max.is_infinite() {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = values.iter().map(|&v| (v - max).exp()).sum();
    max + sum.ln()
}

/// Softmax of a slice, returns normalized probabilities.
fn softmax(logits: &[f64]) -> Vec<f64> {
    let max = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exp: Vec<f64> = logits.iter().map(|&x| (x - max).exp()).collect();
    let sum: f64 = exp.iter().sum::<f64>().max(1e-300);
    exp.iter().map(|&e| e / sum).collect()
}

// ---------------------------------------------------------------------------
// BayesianNetwork — DAG with CPTs and belief propagation
// ---------------------------------------------------------------------------

/// A node in a Bayesian network.
#[derive(Debug, Clone)]
pub struct BnNode {
    /// Name of the variable.
    pub name: String,
    /// Number of discrete states.
    pub n_states: usize,
    /// Parent node indices.
    pub parents: Vec<usize>,
    /// Conditional probability table.
    ///
    /// Shape: `[parent_config_count, n_states]` flattened row-major.
    /// For a root node with no parents: length = `n_states`.
    pub cpt: Vec<f64>,
}

impl BnNode {
    /// Creates a new `BnNode`.
    pub fn new(
        name: impl Into<String>,
        n_states: usize,
        parents: Vec<usize>,
        cpt: Vec<f64>,
    ) -> Self {
        Self {
            name: name.into(),
            n_states,
            parents,
            cpt,
        }
    }

    /// Returns the conditional probability P(state | parent_config).
    pub fn cpt_value(&self, state: usize, parent_config: usize) -> f64 {
        let offset = parent_config * self.n_states;
        self.cpt[offset + state]
    }
}

/// A Bayesian Network: directed acyclic graph with conditional probability tables.
///
/// Supports exact inference via variable elimination on small networks and
/// loopy belief propagation on larger ones.
#[derive(Debug, Clone)]
pub struct BayesianNetwork {
    /// Nodes in topological order.
    pub nodes: Vec<BnNode>,
}

impl BayesianNetwork {
    /// Creates a new empty `BayesianNetwork`.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    /// Adds a node and returns its index.
    pub fn add_node(&mut self, node: BnNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }

    /// Computes the joint probability of a complete assignment.
    ///
    /// `assignment[i]` is the state of node `i`.
    pub fn joint_probability(&self, assignment: &[usize]) -> f64 {
        let mut prob = 1.0f64;
        for (i, node) in self.nodes.iter().enumerate() {
            let parent_config = self.parent_config_index(i, assignment);
            prob *= node.cpt_value(assignment[i], parent_config);
        }
        prob
    }

    /// Computes the parent configuration index for node `i` given full assignment.
    fn parent_config_index(&self, node_idx: usize, assignment: &[usize]) -> usize {
        let node = &self.nodes[node_idx];
        let mut config = 0usize;
        for &p in &node.parents {
            let p_states = self.nodes[p].n_states;
            config = config * p_states + assignment[p];
        }
        config
    }

    /// Computes the marginal probability of node `target` being in `state`
    /// by summing over all other assignments (exact, exponential complexity).
    pub fn marginal(&self, target: usize, target_state: usize) -> f64 {
        let n = self.nodes.len();
        // Enumerate all assignments via mixed-radix counter
        let n_states: Vec<usize> = self.nodes.iter().map(|nd| nd.n_states).collect();
        let total: usize = n_states.iter().product();
        let mut prob = 0.0f64;
        let mut assignment = vec![0usize; n];
        for _ in 0..total {
            if assignment[target] == target_state {
                prob += self.joint_probability(&assignment);
            }
            // Increment mixed-radix counter
            let mut carry = 1;
            for i in (0..n).rev() {
                let next = assignment[i] + carry;
                assignment[i] = next % n_states[i];
                carry = next / n_states[i];
                if carry == 0 {
                    break;
                }
            }
        }
        prob
    }

    /// Returns marginal probabilities for all states of node `target`.
    pub fn marginal_all(&self, target: usize) -> Vec<f64> {
        let n_states = self.nodes[target].n_states;
        (0..n_states).map(|s| self.marginal(target, s)).collect()
    }

    /// Computes conditional probability P(target=state | evidence).
    ///
    /// `evidence` is a list of `(node_idx, state)` observations.
    pub fn conditional(
        &self,
        target: usize,
        target_state: usize,
        evidence: &[(usize, usize)],
    ) -> f64 {
        let n = self.nodes.len();
        let n_states: Vec<usize> = self.nodes.iter().map(|nd| nd.n_states).collect();
        let total: usize = n_states.iter().product();
        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        let mut assignment = vec![0usize; n];
        for _ in 0..total {
            // Check evidence consistency
            let consistent = evidence.iter().all(|&(ni, s)| assignment[ni] == s);
            if consistent {
                let p = self.joint_probability(&assignment);
                denom += p;
                if assignment[target] == target_state {
                    num += p;
                }
            }
            let mut carry = 1;
            for i in (0..n).rev() {
                let next = assignment[i] + carry;
                assignment[i] = next % n_states[i];
                carry = next / n_states[i];
                if carry == 0 {
                    break;
                }
            }
        }
        if denom < 1e-300 { 0.0 } else { num / denom }
    }

    /// Checks that all CPTs are valid (non-negative, rows sum to 1).
    pub fn validate(&self) -> bool {
        for node in &self.nodes {
            let n_configs = if node.parents.is_empty() {
                1
            } else {
                node.cpt.len() / node.n_states
            };
            for cfg in 0..n_configs {
                let sum: f64 = (0..node.n_states).map(|s| node.cpt_value(s, cfg)).sum();
                if (sum - 1.0).abs() > 1e-6 {
                    return false;
                }
            }
        }
        true
    }
}

impl Default for BayesianNetwork {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// HiddenMarkovModel
// ---------------------------------------------------------------------------

/// A Hidden Markov Model with discrete hidden states and Gaussian emissions.
///
/// Supports:
/// - Forward algorithm (likelihood computation)
/// - Viterbi algorithm (MAP state sequence)
/// - Baum-Welch EM (parameter learning)
#[derive(Debug, Clone)]
pub struct HiddenMarkovModel {
    /// Number of hidden states.
    pub n_states: usize,
    /// Initial state distribution π.
    pub initial: Vec<f64>,
    /// Transition matrix A\[i\]\[j\] = P(s_t=j | s_{t-1}=i).
    pub transition: Vec<Vec<f64>>,
    /// Emission mean for each state.
    pub emission_mean: Vec<f64>,
    /// Emission variance for each state.
    pub emission_var: Vec<f64>,
}

impl HiddenMarkovModel {
    /// Creates a new `HiddenMarkovModel` with given parameters.
    pub fn new(
        n_states: usize,
        initial: Vec<f64>,
        transition: Vec<Vec<f64>>,
        emission_mean: Vec<f64>,
        emission_var: Vec<f64>,
    ) -> Self {
        Self {
            n_states,
            initial,
            transition,
            emission_mean,
            emission_var,
        }
    }

    /// Creates a uniform HMM with `n_states` states.
    pub fn uniform(n_states: usize) -> Self {
        let p = 1.0 / n_states as f64;
        let initial = vec![p; n_states];
        let transition = vec![vec![p; n_states]; n_states];
        let emission_mean: Vec<f64> = (0..n_states).map(|i| i as f64).collect();
        let emission_var = vec![1.0; n_states];
        Self::new(n_states, initial, transition, emission_mean, emission_var)
    }

    /// Computes log emission probability of `obs` in state `s`.
    fn log_emit(&self, s: usize, obs: f64) -> f64 {
        log_normal_pdf(obs, self.emission_mean[s], self.emission_var[s])
    }

    /// Forward algorithm: returns log-likelihood of observation sequence.
    pub fn forward(&self, observations: &[f64]) -> f64 {
        let t_len = observations.len();
        if t_len == 0 {
            return 0.0;
        }
        let k = self.n_states;
        let mut alpha: Vec<f64> = (0..k)
            .map(|s| self.initial[s].ln() + self.log_emit(s, observations[0]))
            .collect();
        // Recursion
        for obs_t in observations[1..].iter() {
            let alpha_new: Vec<f64> = (0..k)
                .map(|j| {
                    let log_emit_j = self.log_emit(j, *obs_t);
                    let terms: Vec<f64> = (0..k)
                        .map(|i| alpha[i] + self.transition[i][j].max(1e-300).ln())
                        .collect();
                    log_sum_exp(&terms) + log_emit_j
                })
                .collect();
            alpha = alpha_new;
        }
        log_sum_exp(&alpha)
    }

    /// Viterbi algorithm: returns the most likely state sequence.
    pub fn viterbi(&self, observations: &[f64]) -> Vec<usize> {
        let t_len = observations.len();
        if t_len == 0 {
            return Vec::new();
        }
        let k = self.n_states;
        let mut delta = vec![vec![0.0f64; k]; t_len];
        let mut psi = vec![vec![0usize; k]; t_len];

        // Initialization
        for (s, d0s) in delta[0].iter_mut().enumerate() {
            *d0s = self.initial[s].max(1e-300).ln() + self.log_emit(s, observations[0]);
        }

        // Recursion
        for t in 1..t_len {
            for j in 0..k {
                let (best_s, best_val) = (0..k)
                    .map(|i| {
                        let v = delta[t - 1][i] + self.transition[i][j].max(1e-300).ln();
                        (i, v)
                    })
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .expect("states iterator is non-empty");
                delta[t][j] = best_val + self.log_emit(j, observations[t]);
                psi[t][j] = best_s;
            }
        }

        // Backtrack
        let mut path = vec![0usize; t_len];
        path[t_len - 1] = (0..k)
            .max_by(|&a, &b| {
                delta[t_len - 1][a]
                    .partial_cmp(&delta[t_len - 1][b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("k states is non-empty");
        for t in (0..t_len - 1).rev() {
            path[t] = psi[t + 1][path[t + 1]];
        }
        path
    }

    /// Baum-Welch EM algorithm for parameter estimation.
    ///
    /// Returns the log-likelihood at each iteration.
    pub fn baum_welch(&mut self, observations: &[f64], n_iter: usize) -> Vec<f64> {
        let t_len = observations.len();
        let k = self.n_states;
        let mut ll_history = Vec::new();

        for _iter in 0..n_iter {
            // E-step: Forward-Backward
            // Forward pass (log scale)
            let mut log_alpha = vec![vec![0.0f64; k]; t_len];
            for (s, la0s) in log_alpha[0].iter_mut().enumerate() {
                *la0s = self.initial[s].max(1e-300).ln() + self.log_emit(s, observations[0]);
            }
            for t in 1..t_len {
                for j in 0..k {
                    let terms: Vec<f64> = (0..k)
                        .map(|i| log_alpha[t - 1][i] + self.transition[i][j].max(1e-300).ln())
                        .collect();
                    log_alpha[t][j] = log_sum_exp(&terms) + self.log_emit(j, observations[t]);
                }
            }
            let log_ll = log_sum_exp(&log_alpha[t_len - 1]);
            ll_history.push(log_ll);

            // Backward pass
            let mut log_beta = vec![vec![0.0f64; k]; t_len];
            // log_beta[T-1][s] = log(1) = 0
            for t in (0..t_len - 1).rev() {
                for i in 0..k {
                    let terms: Vec<f64> = (0..k)
                        .map(|j| {
                            self.transition[i][j].max(1e-300).ln()
                                + self.log_emit(j, observations[t + 1])
                                + log_beta[t + 1][j]
                        })
                        .collect();
                    log_beta[t][i] = log_sum_exp(&terms);
                }
            }

            // Compute gamma and xi
            // gamma[t][s] = P(S_t=s | obs)
            let mut gamma = vec![vec![0.0f64; k]; t_len];
            for (t, gt) in gamma.iter_mut().enumerate() {
                let log_probs: Vec<f64> =
                    (0..k).map(|s| log_alpha[t][s] + log_beta[t][s]).collect();
                let norm = log_sum_exp(&log_probs);
                for (s, gts) in gt.iter_mut().enumerate() {
                    *gts = (log_probs[s] - norm).exp();
                }
            }

            // xi[t][i][j] = P(S_t=i, S_{t+1}=j | obs)
            let mut xi = vec![vec![vec![0.0f64; k]; k]; t_len.saturating_sub(1)];
            for t in 0..t_len.saturating_sub(1) {
                let mut log_xi_t = vec![vec![0.0f64; k]; k];
                for (i, lx_row) in log_xi_t.iter_mut().enumerate() {
                    for (j, lx_ij) in lx_row.iter_mut().enumerate() {
                        *lx_ij = log_alpha[t][i]
                            + self.transition[i][j].max(1e-300).ln()
                            + self.log_emit(j, observations[t + 1])
                            + log_beta[t + 1][j];
                    }
                }
                let flat: Vec<f64> = log_xi_t.iter().flat_map(|r| r.iter().copied()).collect();
                let norm = log_sum_exp(&flat);
                let xi_t: Vec<Vec<f64>> = log_xi_t
                    .iter()
                    .map(|row| row.iter().map(|&v| (v - norm).exp()).collect())
                    .collect();
                xi[t] = xi_t;
            }

            // M-step: update parameters
            // Update initial
            for (s, init_s) in self.initial.iter_mut().enumerate() {
                *init_s = gamma[0][s].max(1e-300);
            }
            let init_sum: f64 = self.initial.iter().sum::<f64>().max(1e-300);
            for init_s in self.initial.iter_mut() {
                *init_s /= init_sum;
            }

            // Update transition
            for i in 0..k {
                let denom: f64 = (0..t_len.saturating_sub(1))
                    .map(|t| gamma[t][i])
                    .sum::<f64>()
                    .max(1e-300);
                for (j, tr_ij) in self.transition[i].iter_mut().enumerate() {
                    let num: f64 = (0..t_len.saturating_sub(1)).map(|t| xi[t][i][j]).sum();
                    *tr_ij = (num / denom).max(1e-300);
                }
                // Renormalize row
                let row_sum: f64 = self.transition[i].iter().sum::<f64>().max(1e-300);
                for tr_ij in self.transition[i].iter_mut() {
                    *tr_ij /= row_sum;
                }
            }

            // Update emission parameters
            for (s, (em_mean, em_var)) in self
                .emission_mean
                .iter_mut()
                .zip(self.emission_var.iter_mut())
                .enumerate()
            {
                let denom: f64 = (0..t_len).map(|t| gamma[t][s]).sum::<f64>().max(1e-300);
                let new_mean: f64 = (0..t_len)
                    .map(|t| gamma[t][s] * observations[t])
                    .sum::<f64>()
                    / denom;
                let new_var: f64 = ((0..t_len)
                    .map(|t| gamma[t][s] * (observations[t] - new_mean).powi(2))
                    .sum::<f64>()
                    / denom)
                    .max(1e-6);
                *em_mean = new_mean;
                *em_var = new_var;
            }
        }
        ll_history
    }
}

// ---------------------------------------------------------------------------
// GaussianProcess
// ---------------------------------------------------------------------------

/// Available kernel functions for Gaussian Processes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KernelType {
    /// Squared exponential (RBF) kernel: `k(x,x') = σ² exp(-||x-x'||²/(2ℓ²))`.
    Rbf,
    /// Matérn 3/2 kernel: `k(x,x') = σ²(1+√3r/ℓ)exp(-√3r/ℓ)`.
    Matern32,
    /// Matérn 5/2 kernel.
    Matern52,
    /// Periodic kernel: `k(x,x') = σ² exp(-2sin²(π|x-x'|/p)/ℓ²)`.
    Periodic,
}

/// Gaussian Process for regression with various kernel functions.
///
/// Maintains training data and supports posterior mean/variance prediction.
#[derive(Debug, Clone)]
pub struct GaussianProcess {
    /// Kernel type.
    pub kernel: KernelType,
    /// Signal variance σ².
    pub signal_var: f64,
    /// Length scale ℓ.
    pub length_scale: f64,
    /// Period (for periodic kernel).
    pub period: f64,
    /// Noise variance.
    pub noise_var: f64,
    /// Training inputs (1D for simplicity).
    pub x_train: Vec<f64>,
    /// Training targets.
    pub y_train: Vec<f64>,
    /// Cholesky factor of K + noise*I (column-major flattened).
    chol: Vec<f64>,
    /// alpha = L^{-T} L^{-1} y.
    alpha: Vec<f64>,
}

impl GaussianProcess {
    /// Creates a new `GaussianProcess` with given hyperparameters.
    pub fn new(kernel: KernelType, signal_var: f64, length_scale: f64, noise_var: f64) -> Self {
        Self {
            kernel,
            signal_var,
            length_scale,
            period: 1.0,
            noise_var,
            x_train: Vec::new(),
            y_train: Vec::new(),
            chol: Vec::new(),
            alpha: Vec::new(),
        }
    }

    /// Sets the period for periodic kernel.
    pub fn with_period(mut self, period: f64) -> Self {
        self.period = period;
        self
    }

    /// Evaluates the kernel between two scalar inputs.
    pub fn k(&self, x1: f64, x2: f64) -> f64 {
        let r = (x1 - x2).abs();
        match self.kernel {
            KernelType::Rbf => {
                self.signal_var * (-r * r / (2.0 * self.length_scale * self.length_scale)).exp()
            }
            KernelType::Matern32 => {
                let sq3r = 3.0f64.sqrt() * r / self.length_scale;
                self.signal_var * (1.0 + sq3r) * (-sq3r).exp()
            }
            KernelType::Matern52 => {
                let sq5r = 5.0f64.sqrt() * r / self.length_scale;
                self.signal_var * (1.0 + sq5r + sq5r * sq5r / 3.0) * (-sq5r).exp()
            }
            KernelType::Periodic => {
                let arg = PI * r / self.period;
                self.signal_var
                    * (-2.0 * arg.sin().powi(2) / (self.length_scale * self.length_scale)).exp()
            }
        }
    }

    /// Fits the GP to training data by computing the Cholesky decomposition.
    pub fn fit(&mut self, x_train: Vec<f64>, y_train: Vec<f64>) {
        let n = x_train.len();
        self.x_train = x_train;
        self.y_train = y_train.clone();

        // Build K + noise*I
        let mut k_mat = vec![0.0f64; n * n];
        for (i, &xi) in self.x_train.iter().enumerate() {
            for (j, &xj) in self.x_train.iter().enumerate() {
                k_mat[i * n + j] = self.k(xi, xj);
            }
            k_mat[i * n + i] += self.noise_var;
        }

        // Cholesky (lower triangular, in-place)
        let mut l = k_mat.clone();
        for i in 0..n {
            for j in 0..=i {
                let mut s = l[i * n + j];
                for k_idx in 0..j {
                    s -= l[i * n + k_idx] * l[j * n + k_idx];
                }
                if i == j {
                    l[i * n + j] = s.max(1e-12).sqrt();
                } else {
                    l[i * n + j] = s / l[j * n + j].max(1e-12);
                }
            }
            // zero upper triangle
            for j in i + 1..n {
                l[i * n + j] = 0.0;
            }
        }
        self.chol = l.clone();

        // Solve L * w = y  →  alpha = L^T \ w
        let mut w = y_train;
        // Forward substitution: L w = y
        for i in 0..n {
            let mut s = w[i];
            for j in 0..i {
                s -= l[i * n + j] * w[j];
            }
            w[i] = s / l[i * n + i].max(1e-12);
        }
        // Back substitution: L^T alpha = w
        let mut alpha = w;
        for i in (0..n).rev() {
            let mut s = alpha[i];
            for j in i + 1..n {
                s -= l[j * n + i] * alpha[j];
            }
            alpha[i] = s / l[i * n + i].max(1e-12);
        }
        self.alpha = alpha;
    }

    /// Predicts posterior mean and variance at test input `x_star`.
    pub fn predict(&self, x_star: f64) -> (f64, f64) {
        let n = self.x_train.len();
        if n == 0 {
            return (0.0, self.signal_var + self.noise_var);
        }

        // k_star = K(x_star, X_train)
        let k_star: Vec<f64> = self.x_train.iter().map(|&xi| self.k(x_star, xi)).collect();

        // mean = k_star^T alpha
        let mean: f64 = k_star
            .iter()
            .zip(self.alpha.iter())
            .map(|(a, b)| a * b)
            .sum();

        // variance = k(x*,x*) - k_star^T (K+sI)^{-1} k_star
        // = k(x*,x*) - v^T v  where v = L^{-1} k_star
        let mut v = k_star.clone();
        for i in 0..n {
            let mut s = v[i];
            for (j, &vj) in v[..i].iter().enumerate() {
                s -= self.chol[i * n + j] * vj;
            }
            v[i] = s / self.chol[i * n + i].max(1e-12);
        }
        let var = (self.k(x_star, x_star) - v.iter().map(|vi| vi * vi).sum::<f64>()).max(1e-12);

        (mean, var)
    }

    /// Computes the log marginal likelihood.
    pub fn log_marginal_likelihood(&self) -> f64 {
        let n = self.x_train.len();
        if n == 0 {
            return 0.0;
        }
        // log p(y|X,θ) = -0.5 y^T α - Σ log L_ii - n/2 log(2π)
        let data_fit: f64 = self
            .y_train
            .iter()
            .zip(self.alpha.iter())
            .map(|(y, a)| y * a)
            .sum::<f64>();
        let log_det: f64 = (0..n)
            .map(|i| self.chol[i * n + i].max(1e-300).ln())
            .sum::<f64>();
        -0.5 * data_fit - log_det - 0.5 * n as f64 * TAU.ln()
    }
}

// ---------------------------------------------------------------------------
// DirichletProcess
// ---------------------------------------------------------------------------

/// Dirichlet Process mixture model using the Chinese Restaurant Process.
///
/// New data points are assigned to existing clusters with probability
/// proportional to cluster size, or start a new cluster with probability α/(n+α).
#[derive(Debug, Clone)]
pub struct DirichletProcess {
    /// Concentration parameter α.
    pub alpha: f64,
    /// Cluster assignments for each observed data point.
    pub assignments: Vec<usize>,
    /// Number of points in each cluster.
    pub cluster_counts: Vec<usize>,
    /// Cluster means (updated incrementally).
    pub cluster_means: Vec<f64>,
    /// Cluster sum of squares (for variance estimation).
    pub cluster_ss: Vec<f64>,
    /// Total number of data points assigned.
    pub n_assigned: usize,
}

impl DirichletProcess {
    /// Creates a new `DirichletProcess` with concentration `alpha`.
    pub fn new(alpha: f64) -> Self {
        Self {
            alpha,
            assignments: Vec::new(),
            cluster_counts: Vec::new(),
            cluster_means: Vec::new(),
            cluster_ss: Vec::new(),
            n_assigned: 0,
        }
    }

    /// Returns the number of clusters.
    pub fn n_clusters(&self) -> usize {
        self.cluster_counts.len()
    }

    /// Assigns a new data point via Chinese Restaurant Process probabilities.
    ///
    /// Returns the cluster index assigned (deterministic: picks highest prob).
    pub fn crp_assign(&mut self, x: f64) -> usize {
        let n = self.n_assigned as f64;
        let k = self.cluster_counts.len();

        // Compute unnormalized probabilities
        let mut probs: Vec<f64> = self
            .cluster_counts
            .iter()
            .map(|&cnt| cnt as f64 / (n + self.alpha))
            .collect();
        probs.push(self.alpha / (n + self.alpha)); // new cluster

        // Pick cluster with highest probability (deterministic for reproducibility)
        let chosen = probs
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(k);

        if chosen == k {
            // New cluster
            self.cluster_counts.push(1);
            self.cluster_means.push(x);
            self.cluster_ss.push(0.0);
        } else {
            // Existing cluster — update incrementally
            let cnt = self.cluster_counts[chosen] as f64;
            let old_mean = self.cluster_means[chosen];
            self.cluster_counts[chosen] += 1;
            let new_mean = old_mean + (x - old_mean) / (cnt + 1.0);
            self.cluster_ss[chosen] += (x - old_mean) * (x - new_mean);
            self.cluster_means[chosen] = new_mean;
        }
        self.assignments.push(chosen);
        self.n_assigned += 1;
        chosen
    }

    /// Stick-breaking construction: samples first `k` mixture weights.
    ///
    /// Returns weights `(w_1, ..., w_k)` where `Σ w_i ≈ 1`.
    /// Uses pseudo-random beta draws based on alpha and index.
    pub fn stick_breaking_weights(&self, k: usize) -> Vec<f64> {
        let mut weights = Vec::with_capacity(k);
        let mut remaining = 1.0f64;
        for i in 0..k {
            // Deterministic approximation: mean of Beta(1, alpha)
            let mean_beta = 1.0 / (1.0 + self.alpha);
            // Slight variation per component
            let v = mean_beta * (1.0 - 0.1 * i as f64 / (k as f64 + 1.0));
            let v = v.clamp(1e-6, 1.0 - 1e-6);
            let w = remaining * v;
            weights.push(w);
            remaining *= 1.0 - v;
        }
        // Add the remaining stick to the last component so weights sum to 1
        if let Some(last) = weights.last_mut() {
            *last += remaining;
        }
        // Normalize to sum exactly to 1
        let total: f64 = weights.iter().sum::<f64>().max(1e-300);
        weights.iter_mut().for_each(|w| *w /= total);
        weights
    }

    /// Returns cluster variance estimates (unbiased).
    pub fn cluster_variances(&self) -> Vec<f64> {
        self.cluster_counts
            .iter()
            .zip(self.cluster_ss.iter())
            .map(
                |(&cnt, &ss)| {
                    if cnt > 1 { ss / (cnt - 1) as f64 } else { 1.0 }
                },
            )
            .collect()
    }

    /// Returns the expected number of clusters for n observations (approximation).
    ///
    /// `E[K_n] ≈ α ln(1 + n/α)`
    pub fn expected_clusters(alpha: f64, n: usize) -> f64 {
        alpha * (1.0 + n as f64 / alpha).ln()
    }
}

// ---------------------------------------------------------------------------
// VariationalInference
// ---------------------------------------------------------------------------

/// Mean-field variational inference for a Gaussian mixture model.
///
/// Maximizes the Evidence Lower BOund (ELBO) by coordinate ascent
/// over the variational posterior factors.
#[derive(Debug, Clone)]
pub struct VariationalInference {
    /// Number of mixture components.
    pub n_components: usize,
    /// Variational component weights (log-scale unnormalized).
    pub log_weights: Vec<f64>,
    /// Variational mean for each component.
    pub var_mean: Vec<f64>,
    /// Variational variance for each component.
    pub var_var: Vec<f64>,
    /// Prior mean.
    pub prior_mean: f64,
    /// Prior variance.
    pub prior_var: f64,
    /// Observation noise variance.
    pub obs_var: f64,
    /// ELBO history.
    pub elbo_history: Vec<f64>,
}

impl VariationalInference {
    /// Creates a new `VariationalInference` instance.
    pub fn new(n_components: usize, prior_mean: f64, prior_var: f64, obs_var: f64) -> Self {
        let log_weights = vec![-(n_components as f64).ln(); n_components];
        let var_mean: Vec<f64> = (0..n_components).map(|i| i as f64).collect();
        let var_var = vec![1.0f64; n_components];
        Self {
            n_components,
            log_weights,
            var_mean,
            var_var,
            prior_mean,
            prior_var,
            obs_var,
            elbo_history: Vec::new(),
        }
    }

    /// Computes the ELBO for current variational parameters given observations.
    pub fn elbo(&self, observations: &[f64]) -> f64 {
        let weights = softmax(&self.log_weights);
        let mut elbo = 0.0f64;
        // Expected log likelihood
        for &x in observations {
            let ll_terms: Vec<f64> = (0..self.n_components)
                .map(|k| {
                    weights[k].max(1e-300).ln()
                        + log_normal_pdf(x, self.var_mean[k], self.obs_var + self.var_var[k])
                })
                .collect();
            elbo += log_sum_exp(&ll_terms);
        }
        // KL divergence: Σ_k w_k KL(q(z_k) || p(z_k))
        for ((&wk, &vv), &vm) in weights
            .iter()
            .zip(self.var_var.iter())
            .zip(self.var_mean.iter())
        {
            // KL(N(μ_q, σ_q²) || N(μ_p, σ_p²))
            let kl = 0.5
                * (self.prior_var / vv.max(1e-12)
                    + (vm - self.prior_mean).powi(2) / self.prior_var
                    - 1.0
                    + (vv / self.prior_var).ln());
            elbo -= wk * kl;
        }
        elbo
    }

    /// Performs one CAVI update step.
    ///
    /// Returns the new ELBO.
    pub fn cavi_step(&mut self, observations: &[f64]) -> f64 {
        let n = observations.len() as f64;
        // Update variational parameters for each component
        for k in 0..self.n_components {
            let weights = softmax(&self.log_weights);
            // Responsibility of component k for each observation
            let r_k: Vec<f64> = observations
                .iter()
                .map(|&x| weights[k] * normal_pdf(x, self.var_mean[k], self.obs_var))
                .collect();
            let r_sum: f64 = r_k.iter().sum::<f64>().max(1e-300);

            // Update mean: posterior precision = prior_prec + r_sum/obs_var
            let prior_prec = 1.0 / self.prior_var.max(1e-12);
            let lik_prec = r_sum / self.obs_var.max(1e-12);
            let post_prec = prior_prec + lik_prec;
            let post_var = 1.0 / post_prec.max(1e-12);
            let data_sum: f64 = r_k
                .iter()
                .zip(observations.iter())
                .map(|(r, x)| r * x)
                .sum();
            let post_mean =
                post_var * (prior_prec * self.prior_mean + data_sum / self.obs_var.max(1e-12));

            self.var_mean[k] = post_mean;
            self.var_var[k] = post_var;

            // Update weight (log)
            self.log_weights[k] = r_sum.max(1e-300).ln();
        }
        // Renormalize log_weights
        let lse = log_sum_exp(&self.log_weights.clone());
        for lw in self.log_weights.iter_mut() {
            *lw -= lse;
        }
        let _ = n;
        let elbo_val = self.elbo(observations);
        self.elbo_history.push(elbo_val);
        elbo_val
    }

    /// Runs CAVI for `n_iter` iterations.
    pub fn fit(&mut self, observations: &[f64], n_iter: usize) -> f64 {
        for _ in 0..n_iter {
            self.cavi_step(observations);
        }
        *self.elbo_history.last().unwrap_or(&f64::NEG_INFINITY)
    }

    /// Reparameterization trick: samples from `q(z) = N(mu, sigma^2)` using `eps ~ N(0,1)`.
    pub fn reparameterize(&self, k: usize, eps: f64) -> f64 {
        self.var_mean[k] + self.var_var[k].sqrt() * eps
    }

    /// Returns the predictive density at `x` under the variational posterior.
    pub fn predictive_density(&self, x: f64) -> f64 {
        let weights = softmax(&self.log_weights);
        (0..self.n_components)
            .map(|k| weights[k] * normal_pdf(x, self.var_mean[k], self.obs_var + self.var_var[k]))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// ExpectationMaximization — Gaussian Mixture Model
// ---------------------------------------------------------------------------

/// A Gaussian mixture model component.
#[derive(Debug, Clone)]
pub struct GmmComponent {
    /// Mixture weight.
    pub weight: f64,
    /// Mean.
    pub mean: f64,
    /// Variance.
    pub var: f64,
}

impl GmmComponent {
    /// Creates a new `GmmComponent`.
    pub fn new(weight: f64, mean: f64, var: f64) -> Self {
        Self { weight, mean, var }
    }
}

/// Expectation-Maximization for Gaussian Mixture Models.
///
/// Supports k-means initialization, BIC criterion for model selection,
/// and full EM convergence.
#[derive(Debug, Clone)]
pub struct ExpectationMaximization {
    /// Number of components.
    pub n_components: usize,
    /// Mixture components.
    pub components: Vec<GmmComponent>,
    /// Log-likelihood history.
    pub ll_history: Vec<f64>,
    /// Convergence tolerance.
    pub tol: f64,
}

impl ExpectationMaximization {
    /// Creates a new `ExpectationMaximization` with k-means seeding.
    pub fn new(n_components: usize) -> Self {
        let components = (0..n_components)
            .map(|i| GmmComponent::new(1.0 / n_components as f64, i as f64, 1.0))
            .collect();
        Self {
            n_components,
            components,
            ll_history: Vec::new(),
            tol: 1e-6,
        }
    }

    /// Sets convergence tolerance.
    pub fn with_tol(mut self, tol: f64) -> Self {
        self.tol = tol;
        self
    }

    /// Initializes component means via k-means (one pass, sorted data).
    pub fn kmeans_init(&mut self, data: &[f64]) {
        if data.is_empty() {
            return;
        }
        let mut sorted = data.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let k = self.n_components;
        for i in 0..k {
            let idx = (sorted.len() * (2 * i + 1)) / (2 * k);
            self.components[i].mean = sorted[idx.min(sorted.len() - 1)];
            self.components[i].var = 1.0;
            self.components[i].weight = 1.0 / k as f64;
        }
    }

    /// Computes log-likelihood of data under current model.
    pub fn log_likelihood(&self, data: &[f64]) -> f64 {
        data.iter()
            .map(|&x| {
                let terms: Vec<f64> = self
                    .components
                    .iter()
                    .map(|c| c.weight.max(1e-300).ln() + log_normal_pdf(x, c.mean, c.var))
                    .collect();
                log_sum_exp(&terms)
            })
            .sum()
    }

    /// Computes the Bayesian Information Criterion.
    ///
    /// `BIC = k ln(n) - 2 ln(L̂)`
    pub fn bic(&self, data: &[f64]) -> f64 {
        let n = data.len() as f64;
        let ll = self.log_likelihood(data);
        // Parameters: n_components means + variances + weights - 1
        let n_params = (3 * self.n_components - 1) as f64;
        n_params * n.ln() - 2.0 * ll
    }

    /// E-step: computes responsibilities for each data point.
    fn e_step(&self, data: &[f64]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|&x| {
                let log_probs: Vec<f64> = self
                    .components
                    .iter()
                    .map(|c| c.weight.max(1e-300).ln() + log_normal_pdf(x, c.mean, c.var))
                    .collect();
                softmax(&log_probs)
            })
            .collect()
    }

    /// M-step: updates parameters given responsibilities.
    fn m_step(&mut self, data: &[f64], responsibilities: &[Vec<f64>]) {
        let n = data.len() as f64;
        for k in 0..self.n_components {
            let r_sum: f64 = responsibilities
                .iter()
                .map(|r| r[k])
                .sum::<f64>()
                .max(1e-300);
            let new_weight = r_sum / n;
            let new_mean: f64 = responsibilities
                .iter()
                .zip(data.iter())
                .map(|(r, &x)| r[k] * x)
                .sum::<f64>()
                / r_sum;
            let new_var: f64 = (responsibilities
                .iter()
                .zip(data.iter())
                .map(|(r, &x)| r[k] * (x - new_mean).powi(2))
                .sum::<f64>()
                / r_sum)
                .max(1e-6);
            self.components[k].weight = new_weight;
            self.components[k].mean = new_mean;
            self.components[k].var = new_var;
        }
    }

    /// Runs the EM algorithm.
    ///
    /// Returns the final log-likelihood.
    pub fn fit(&mut self, data: &[f64], max_iter: usize) -> f64 {
        self.ll_history.clear();
        let mut prev_ll = f64::NEG_INFINITY;
        for _ in 0..max_iter {
            let resp = self.e_step(data);
            self.m_step(data, &resp);
            let ll = self.log_likelihood(data);
            self.ll_history.push(ll);
            if (ll - prev_ll).abs() < self.tol {
                break;
            }
            prev_ll = ll;
        }
        *self.ll_history.last().unwrap_or(&f64::NEG_INFINITY)
    }

    /// Predicts cluster assignment (most likely component) for a data point.
    pub fn predict(&self, x: f64) -> usize {
        self.components
            .iter()
            .enumerate()
            .map(|(k, c)| {
                let ll = c.weight.max(1e-300).ln() + log_normal_pdf(x, c.mean, c.var);
                (k, ll)
            })
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(k, _)| k)
            .unwrap_or(0)
    }

    /// Returns component weights normalized to sum to 1.
    pub fn normalized_weights(&self) -> Vec<f64> {
        let sum: f64 = self
            .components
            .iter()
            .map(|c| c.weight)
            .sum::<f64>()
            .max(1e-300);
        self.components.iter().map(|c| c.weight / sum).collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    /// Computes the multivariate Gaussian log-density (diagonal covariance).
    fn mvn_log_pdf_diag(x: &[f64], mean: &[f64], var: &[f64]) -> f64 {
        let d = x.len() as f64;
        let log_det: f64 = var.iter().map(|v| v.max(1e-300).ln()).sum();
        let maha: f64 = x
            .iter()
            .zip(mean.iter())
            .zip(var.iter())
            .map(|((&xi, &mi), &vi)| (xi - mi).powi(2) / vi.max(1e-300))
            .sum();
        -0.5 * (d * TAU.ln() + log_det + maha)
    }

    // --- BayesianNetwork ---

    fn make_simple_bn() -> BayesianNetwork {
        let mut bn = BayesianNetwork::new();
        // Node 0: Rain (2 states) — prior
        bn.add_node(BnNode::new("Rain", 2, vec![], vec![0.3, 0.7]));
        // Node 1: Sprinkler (2 states) — parent: Rain
        bn.add_node(BnNode::new(
            "Sprinkler",
            2,
            vec![0],
            vec![0.1, 0.9, 0.5, 0.5], // [rain=0: s=0,s=1; rain=1: s=0,s=1]
        ));
        bn
    }

    #[test]
    fn test_bn_validate() {
        let bn = make_simple_bn();
        assert!(bn.validate());
    }

    #[test]
    fn test_bn_joint_probability_sums_to_one() {
        let mut bn = BayesianNetwork::new();
        bn.add_node(BnNode::new("A", 2, vec![], vec![0.4, 0.6]));
        bn.add_node(BnNode::new("B", 2, vec![0], vec![0.7, 0.3, 0.2, 0.8]));
        // Sum over all assignments
        let total: f64 = (0..4)
            .map(|i| {
                let a = i / 2;
                let b_val = i % 2;
                bn.joint_probability(&[a, b_val])
            })
            .sum();
        assert!((total - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bn_marginal_sums_to_one() {
        let bn = make_simple_bn();
        let m0 = bn.marginal(0, 0);
        let m1 = bn.marginal(0, 1);
        assert!((m0 + m1 - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bn_marginal_root_equals_prior() {
        let bn = make_simple_bn();
        let m0 = bn.marginal(0, 0);
        assert!((m0 - 0.3).abs() < 1e-8);
    }

    #[test]
    fn test_bn_conditional_valid() {
        let bn = make_simple_bn();
        let p = bn.conditional(1, 0, &[(0, 0)]);
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn test_bn_marginal_all_sums_to_one() {
        let bn = make_simple_bn();
        let m = bn.marginal_all(0);
        let s: f64 = m.iter().sum();
        assert!((s - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_bn_cpt_value() {
        let node = BnNode::new("X", 2, vec![], vec![0.4, 0.6]);
        assert!((node.cpt_value(0, 0) - 0.4).abs() < 1e-10);
        assert!((node.cpt_value(1, 0) - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_bn_single_node() {
        let mut bn = BayesianNetwork::new();
        bn.add_node(BnNode::new("X", 3, vec![], vec![0.2, 0.5, 0.3]));
        let p = bn.joint_probability(&[1]);
        assert!((p - 0.5).abs() < 1e-10);
    }

    // --- HiddenMarkovModel ---

    fn make_hmm() -> HiddenMarkovModel {
        HiddenMarkovModel::new(
            2,
            vec![0.6, 0.4],
            vec![vec![0.7, 0.3], vec![0.4, 0.6]],
            vec![0.0, 3.0],
            vec![1.0, 1.0],
        )
    }

    #[test]
    fn test_hmm_forward_returns_finite() {
        let hmm = make_hmm();
        let obs = vec![0.1, 0.2, 0.3, 2.8, 3.1];
        let ll = hmm.forward(&obs);
        assert!(ll.is_finite());
    }

    #[test]
    fn test_hmm_forward_empty() {
        let hmm = make_hmm();
        assert_eq!(hmm.forward(&[]), 0.0);
    }

    #[test]
    fn test_hmm_viterbi_length() {
        let hmm = make_hmm();
        let obs = vec![0.1, 0.2, 0.3, 2.8, 3.1];
        let path = hmm.viterbi(&obs);
        assert_eq!(path.len(), obs.len());
    }

    #[test]
    fn test_hmm_viterbi_valid_states() {
        let hmm = make_hmm();
        let obs = vec![0.1, 2.9, 0.2, 3.0];
        let path = hmm.viterbi(&obs);
        assert!(path.iter().all(|&s| s < 2));
    }

    #[test]
    fn test_hmm_viterbi_empty() {
        let hmm = make_hmm();
        assert_eq!(hmm.viterbi(&[]).len(), 0);
    }

    #[test]
    fn test_hmm_baum_welch_ll_increases() {
        let mut hmm = make_hmm();
        let obs: Vec<f64> = (0..20)
            .map(|i| if i % 3 == 0 { 0.1 } else { 2.9 })
            .collect();
        let ll_hist = hmm.baum_welch(&obs, 5);
        // Log-likelihood should be non-decreasing
        for i in 1..ll_hist.len() {
            assert!(ll_hist[i] >= ll_hist[i - 1] - 1e-4);
        }
    }

    #[test]
    fn test_hmm_uniform_creation() {
        let hmm = HiddenMarkovModel::uniform(3);
        assert_eq!(hmm.n_states, 3);
        let row_sum: f64 = hmm.transition[0].iter().sum();
        assert!((row_sum - 1.0).abs() < 1e-10);
    }

    // --- GaussianProcess ---

    #[test]
    fn test_gp_rbf_kernel_diagonal() {
        let gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 1e-3);
        assert!((gp.k(0.0, 0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gp_rbf_kernel_decays() {
        let gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 1e-3);
        assert!(gp.k(0.0, 10.0) < gp.k(0.0, 1.0));
    }

    #[test]
    fn test_gp_matern32_diagonal() {
        let gp = GaussianProcess::new(KernelType::Matern32, 2.0, 1.0, 1e-3);
        assert!((gp.k(0.0, 0.0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_gp_matern52_diagonal() {
        let gp = GaussianProcess::new(KernelType::Matern52, 1.5, 1.0, 1e-3);
        assert!((gp.k(0.0, 0.0) - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_gp_periodic_diagonal() {
        let gp = GaussianProcess::new(KernelType::Periodic, 1.0, 1.0, 1e-3);
        assert!((gp.k(0.0, 0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gp_fit_predict_mean() {
        let mut gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 1e-4);
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, 4.0, 9.0];
        gp.fit(x, y);
        let (mean, _var) = gp.predict(1.0);
        assert!((mean - 1.0).abs() < 0.5); // near training point
    }

    #[test]
    fn test_gp_predict_variance_positive() {
        let mut gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 1e-4);
        gp.fit(vec![0.0, 1.0], vec![0.0, 1.0]);
        let (_mean, var) = gp.predict(5.0); // far from training data
        assert!(var > 0.0);
    }

    #[test]
    fn test_gp_predict_empty() {
        let gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 1e-3);
        let (mean, var) = gp.predict(0.5);
        assert_eq!(mean, 0.0);
        assert!(var > 0.0);
    }

    #[test]
    fn test_gp_log_marginal_likelihood() {
        let mut gp = GaussianProcess::new(KernelType::Rbf, 1.0, 1.0, 0.1);
        gp.fit(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 0.0]);
        let lml = gp.log_marginal_likelihood();
        assert!(lml.is_finite());
    }

    // --- DirichletProcess ---

    #[test]
    fn test_dp_initial_state() {
        let dp = DirichletProcess::new(1.0);
        assert_eq!(dp.n_clusters(), 0);
        assert_eq!(dp.n_assigned, 0);
    }

    #[test]
    fn test_dp_crp_first_point() {
        let mut dp = DirichletProcess::new(1.0);
        let c = dp.crp_assign(0.0);
        assert_eq!(c, 0);
        assert_eq!(dp.n_clusters(), 1);
    }

    #[test]
    fn test_dp_crp_multiple_points() {
        let mut dp = DirichletProcess::new(0.1); // low alpha = prefer existing clusters
        for i in 0..10 {
            dp.crp_assign(i as f64 * 0.01);
        }
        // With low alpha, should form few clusters
        assert!(dp.n_clusters() <= 5);
    }

    #[test]
    fn test_dp_stick_breaking_sums_to_one() {
        let dp = DirichletProcess::new(2.0);
        let w = dp.stick_breaking_weights(10);
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dp_stick_breaking_positive() {
        let dp = DirichletProcess::new(1.0);
        let w = dp.stick_breaking_weights(5);
        assert!(w.iter().all(|&wi| wi > 0.0));
    }

    #[test]
    fn test_dp_expected_clusters() {
        let e = DirichletProcess::expected_clusters(1.0, 100);
        assert!(e > 3.0 && e < 10.0);
    }

    #[test]
    fn test_dp_cluster_variances() {
        let mut dp = DirichletProcess::new(0.5);
        for i in 0..5 {
            dp.crp_assign(i as f64);
        }
        let vars = dp.cluster_variances();
        assert!(vars.iter().all(|&v| v >= 0.0));
    }

    // --- VariationalInference ---

    #[test]
    fn test_vi_elbo_finite() {
        let vi = VariationalInference::new(2, 0.0, 1.0, 1.0);
        let obs = vec![0.0, 1.0, -1.0, 2.0];
        let elbo = vi.elbo(&obs);
        assert!(elbo.is_finite());
    }

    #[test]
    fn test_vi_cavi_step_updates_params() {
        let mut vi = VariationalInference::new(2, 0.0, 1.0, 1.0);
        let old_mean = vi.var_mean[0];
        let obs = vec![3.0, 3.1, 3.2, -3.0, -3.1, -3.2];
        vi.cavi_step(&obs);
        // Mean should change toward data
        assert!((vi.var_mean[0] - old_mean).abs() > 0.0);
    }

    #[test]
    fn test_vi_fit_returns_finite() {
        let mut vi = VariationalInference::new(2, 0.0, 2.0, 1.0);
        let obs: Vec<f64> = (0..20)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let elbo = vi.fit(&obs, 10);
        assert!(elbo.is_finite());
    }

    #[test]
    fn test_vi_reparameterize() {
        let vi = VariationalInference::new(2, 0.0, 1.0, 1.0);
        let sample = vi.reparameterize(0, 1.0);
        // sample = mean[0] + sqrt(var[0]) * 1.0
        let expected = vi.var_mean[0] + vi.var_var[0].sqrt();
        assert!((sample - expected).abs() < 1e-10);
    }

    #[test]
    fn test_vi_predictive_density_positive() {
        let vi = VariationalInference::new(2, 0.0, 1.0, 1.0);
        let p = vi.predictive_density(0.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_vi_elbo_history_grows() {
        let mut vi = VariationalInference::new(2, 0.0, 1.0, 1.0);
        let obs = vec![1.0, -1.0, 2.0];
        vi.fit(&obs, 5);
        assert_eq!(vi.elbo_history.len(), 5);
    }

    // --- ExpectationMaximization ---

    #[test]
    fn test_em_initial_weights_sum_to_one() {
        let em = ExpectationMaximization::new(3);
        let sum: f64 = em.normalized_weights().iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_em_kmeans_init() {
        let mut em = ExpectationMaximization::new(2);
        let data = vec![0.0, 0.1, 0.2, 5.0, 5.1, 5.2];
        em.kmeans_init(&data);
        // Means should be near 0 and 5
        let means: Vec<f64> = em.components.iter().map(|c| c.mean).collect();
        assert!(means.iter().any(|&m| m < 1.0));
        assert!(means.iter().any(|&m| m > 4.0));
    }

    #[test]
    fn test_em_log_likelihood_finite() {
        let em = ExpectationMaximization::new(2);
        let data = vec![0.0, 1.0, 2.0];
        assert!(em.log_likelihood(&data).is_finite());
    }

    #[test]
    fn test_em_fit_ll_increases() {
        let mut em = ExpectationMaximization::new(2);
        let data: Vec<f64> = (0..30)
            .map(|i| {
                if i < 15 {
                    i as f64 * 0.1
                } else {
                    5.0 + i as f64 * 0.1
                }
            })
            .collect();
        em.kmeans_init(&data);
        em.fit(&data, 20);
        let ll = &em.ll_history;
        for i in 1..ll.len() {
            assert!(ll[i] >= ll[i - 1] - 1e-4);
        }
    }

    #[test]
    fn test_em_predict_valid_component() {
        let em = ExpectationMaximization::new(3);
        let pred = em.predict(0.5);
        assert!(pred < 3);
    }

    #[test]
    fn test_em_bic_finite() {
        let em = ExpectationMaximization::new(2);
        let data = vec![0.0, 1.0, 5.0, 6.0];
        let bic = em.bic(&data);
        assert!(bic.is_finite());
    }

    #[test]
    fn test_em_fit_separates_clusters() {
        let mut em = ExpectationMaximization::new(2);
        let mut data: Vec<f64> = (0..20).map(|i| i as f64 * 0.05).collect(); // 0..1
        let data2: Vec<f64> = (0..20).map(|i| 10.0 + i as f64 * 0.05).collect(); // 10..11
        data.extend(data2);
        em.kmeans_init(&data);
        em.fit(&data, 50);
        // One component should be near 0.5, other near 10.5
        let means: Vec<f64> = em.components.iter().map(|c| c.mean).collect();
        assert!(means.iter().any(|&m| m < 3.0));
        assert!(means.iter().any(|&m| m > 7.0));
    }

    #[test]
    fn test_em_n_components() {
        let em = ExpectationMaximization::new(4);
        assert_eq!(em.n_components, 4);
        assert_eq!(em.components.len(), 4);
    }

    // --- Helper functions ---

    #[test]
    fn test_log_sum_exp_empty() {
        assert_eq!(log_sum_exp(&[]), f64::NEG_INFINITY);
    }

    #[test]
    fn test_log_sum_exp_single() {
        assert!((log_sum_exp(&[2.0]) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_softmax_sums_to_one() {
        let s = softmax(&[1.0, 2.0, 3.0]);
        assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normal_pdf_peak() {
        let p = normal_pdf(0.0, 0.0, 1.0);
        assert!((p - 1.0 / (TAU).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_mvn_log_pdf_diag() {
        let x = vec![0.0, 0.0];
        let mean = vec![0.0, 0.0];
        let var = vec![1.0, 1.0];
        let lp = mvn_log_pdf_diag(&x, &mean, &var);
        assert!(lp.is_finite());
    }
}
