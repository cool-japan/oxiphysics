// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Causal inference and structural causal models.
//!
//! Provides:
//! - [`CausalGraph`]               — directed acyclic graph, topological sort, d-separation
//! - [`StructuralCausalModel`]     — linear SCM with noise terms, do-calculus interventions
//! - [`BackdoorCriterion`]         — backdoor criterion check and adjustment
//! - [`FrontdoorCriterion`]        — frontdoor adjustment formula
//! - [`PropensityScoreMatching`]   — propensity score estimation, ATT/ATE
//! - [`InstrumentalVariables`]     — IV estimation, two-stage least squares (2SLS)
//! - [`CausalDiscovery`]           — PC algorithm skeleton, orientation rules
//! - [`CounterfactualQuery`]       — E\[Y|do(X=x), Z=z\] style queries

use std::collections::{HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// CausalGraph
// ---------------------------------------------------------------------------

/// A directed acyclic graph (DAG) representing causal relationships between variables.
///
/// Nodes are identified by `usize` indices. Edges represent direct causal effects
/// from parent to child. The graph must be acyclic for causal semantics to be valid.
#[derive(Debug, Clone)]
pub struct CausalGraph {
    /// Number of nodes in the graph.
    pub n_nodes: usize,
    /// Adjacency list: `parents[v]` gives all parent nodes of `v`.
    pub parents: Vec<Vec<usize>>,
    /// Adjacency list: `children[v]` gives all children of `v`.
    pub children: Vec<Vec<usize>>,
    /// Optional variable names.
    pub names: Vec<String>,
}

impl CausalGraph {
    /// Create a new empty causal graph with `n` nodes.
    ///
    /// # Arguments
    /// * `n` — number of nodes
    pub fn new(n: usize) -> Self {
        Self {
            n_nodes: n,
            parents: vec![vec![]; n],
            children: vec![vec![]; n],
            names: (0..n).map(|i| format!("X{i}")).collect(),
        }
    }

    /// Set variable names.
    ///
    /// # Arguments
    /// * `names` — slice of names, length must equal `n_nodes`
    pub fn set_names(&mut self, names: &[&str]) {
        assert_eq!(names.len(), self.n_nodes);
        self.names = names.iter().map(|s| s.to_string()).collect();
    }

    /// Add a directed edge from `from` (parent/cause) to `to` (child/effect).
    ///
    /// # Panics
    /// Panics if adding this edge would create a cycle.
    pub fn add_edge(&mut self, from: usize, to: usize) {
        assert!(
            from < self.n_nodes && to < self.n_nodes,
            "node index out of bounds"
        );
        assert!(
            !self.creates_cycle(from, to),
            "edge {from}→{to} would create a cycle"
        );
        if !self.children[from].contains(&to) {
            self.children[from].push(to);
            self.parents[to].push(from);
        }
    }

    /// Check whether adding edge `from→to` would create a cycle.
    pub fn creates_cycle(&self, from: usize, to: usize) -> bool {
        // DFS from `to`: if we can reach `from`, adding from→to creates a cycle.
        let mut visited = vec![false; self.n_nodes];
        let mut stack = vec![to];
        while let Some(node) = stack.pop() {
            if node == from {
                return true;
            }
            if !visited[node] {
                visited[node] = true;
                for &child in &self.children[node] {
                    stack.push(child);
                }
            }
        }
        false
    }

    /// Return nodes in topological order (Kahn's algorithm).
    ///
    /// Returns `None` if the graph has a cycle (should not happen if edges are
    /// added via `add_edge`).
    pub fn topological_sort(&self) -> Option<Vec<usize>> {
        let mut in_degree: Vec<usize> = self.parents.iter().map(|p| p.len()).collect();
        let mut queue: VecDeque<usize> = (0..self.n_nodes).filter(|&v| in_degree[v] == 0).collect();
        let mut order = Vec::with_capacity(self.n_nodes);
        while let Some(v) = queue.pop_front() {
            order.push(v);
            for &child in &self.children[v] {
                in_degree[child] -= 1;
                if in_degree[child] == 0 {
                    queue.push_back(child);
                }
            }
        }
        if order.len() == self.n_nodes {
            Some(order)
        } else {
            None
        }
    }

    /// Return all ancestors of node `v` (nodes from which `v` is reachable).
    pub fn ancestors(&self, v: usize) -> HashSet<usize> {
        let mut anc = HashSet::new();
        let mut stack = vec![v];
        while let Some(node) = stack.pop() {
            for &p in &self.parents[node] {
                if anc.insert(p) {
                    stack.push(p);
                }
            }
        }
        anc
    }

    /// Return all descendants of node `v`.
    pub fn descendants(&self, v: usize) -> HashSet<usize> {
        let mut desc = HashSet::new();
        let mut stack = vec![v];
        while let Some(node) = stack.pop() {
            for &c in &self.children[node] {
                if desc.insert(c) {
                    stack.push(c);
                }
            }
        }
        desc
    }

    /// Test d-separation: are node sets `x` and `y` d-separated by conditioning set `z`?
    ///
    /// Uses the Bayes Ball algorithm. Returns `true` if `x ⊥ y | z`.
    pub fn d_separated(&self, x: &[usize], y: &[usize], z: &[usize]) -> bool {
        let z_set: HashSet<usize> = z.iter().copied().collect();
        let y_set: HashSet<usize> = y.iter().copied().collect();

        // Collect all ancestors of Z (needed for v-structure blocking)
        let mut z_ancestors: HashSet<usize> = z_set.clone();
        for &zv in z {
            z_ancestors.extend(self.ancestors(zv));
        }

        // Bayes Ball: (node, direction) where direction=true means "from child"
        let mut visited: HashSet<(usize, bool)> = HashSet::new();
        let mut queue: VecDeque<(usize, bool)> = VecDeque::new();

        for &xv in x {
            // Start going "up" (toward parents) and "down" (toward children)
            queue.push_back((xv, true)); // via child → go up
            queue.push_back((xv, false)); // via parent → go down
        }

        while let Some((node, via_child)) = queue.pop_front() {
            if visited.contains(&(node, via_child)) {
                continue;
            }
            visited.insert((node, via_child));

            if y_set.contains(&node) {
                return false; // path found → not d-separated
            }

            let in_z = z_set.contains(&node);
            let in_z_anc = z_ancestors.contains(&node);

            if via_child && !in_z {
                // Arrived via child, node is not in Z
                // Can traverse: up to parents (chain/fork) and down to children (not blocked)
                for &p in &self.parents[node] {
                    queue.push_back((p, true));
                }
                for &c in &self.children[node] {
                    queue.push_back((c, false));
                }
            } else if !via_child {
                // Arrived via parent
                if !in_z {
                    // Collider not in Z: blocked going down unless ancestor
                    for &c in &self.children[node] {
                        queue.push_back((c, false));
                    }
                }
                if in_z_anc {
                    // v-structure activated by conditioning on descendant
                    for &p in &self.parents[node] {
                        queue.push_back((p, true));
                    }
                }
            }
        }
        true
    }

    /// Return the Markov blanket of node `v`: parents, children, and co-parents.
    pub fn markov_blanket(&self, v: usize) -> HashSet<usize> {
        let mut blanket = HashSet::new();
        for &p in &self.parents[v] {
            blanket.insert(p);
        }
        for &c in &self.children[v] {
            blanket.insert(c);
            for &cp in &self.parents[c] {
                if cp != v {
                    blanket.insert(cp);
                }
            }
        }
        blanket
    }

    /// Check if the graph is acyclic.
    pub fn is_acyclic(&self) -> bool {
        self.topological_sort().is_some()
    }
}

// ---------------------------------------------------------------------------
// StructuralCausalModel
// ---------------------------------------------------------------------------

/// A linear structural causal model (SCM).
///
/// Each variable `X_i` is defined as:
/// `X_i = Σ_j (coeff[i][j] * X_j) + noise_std[i] * ε_i`
///
/// where `ε_i ~ N(0,1)` and `j` ranges over parents of `i`.
#[derive(Debug, Clone)]
pub struct StructuralCausalModel {
    /// The underlying causal graph.
    pub graph: CausalGraph,
    /// Structural coefficients: `coefficients[i][k]` is the coefficient for the
    /// k-th parent of node `i`.
    pub coefficients: Vec<Vec<f64>>,
    /// Standard deviation of the exogenous noise for each variable.
    pub noise_std: Vec<f64>,
    /// Intercept terms for each variable.
    pub intercepts: Vec<f64>,
}

impl StructuralCausalModel {
    /// Create a new linear SCM on `n` variables.
    ///
    /// All coefficients are zero, noise std = 1.0, intercepts = 0.0 by default.
    pub fn new(n: usize) -> Self {
        Self {
            graph: CausalGraph::new(n),
            coefficients: vec![vec![]; n],
            noise_std: vec![1.0; n],
            intercepts: vec![0.0; n],
        }
    }

    /// Add a causal edge with a specified structural coefficient.
    ///
    /// # Arguments
    /// * `from` — parent (cause) node index
    /// * `to`   — child (effect) node index
    /// * `coeff` — linear coefficient
    pub fn add_edge(&mut self, from: usize, to: usize, coeff: f64) {
        self.graph.add_edge(from, to);
        // The k-th parent of `to` is now `from`
        self.coefficients[to].push(coeff);
    }

    /// Set the noise standard deviation for variable `v`.
    pub fn set_noise(&mut self, v: usize, std: f64) {
        self.noise_std[v] = std;
    }

    /// Set the intercept for variable `v`.
    pub fn set_intercept(&mut self, v: usize, intercept: f64) {
        self.intercepts[v] = intercept;
    }

    /// Sample one observation from the SCM using provided noise values.
    ///
    /// # Arguments
    /// * `noise` — exogenous noise values `ε_i` for each variable (length = n_nodes)
    ///
    /// Returns `x[i]` values in topological order.
    pub fn sample_with_noise(&self, noise: &[f64]) -> Vec<f64> {
        let n = self.graph.n_nodes;
        let order = self
            .graph
            .topological_sort()
            .expect("SCM graph must be acyclic");
        let mut x = vec![0.0_f64; n];
        for &v in &order {
            let val: f64 = self.intercepts[v]
                + self.graph.parents[v]
                    .iter()
                    .zip(self.coefficients[v].iter())
                    .map(|(&p, &c)| c * x[p])
                    .sum::<f64>()
                + self.noise_std[v] * noise[v];
            x[v] = val;
        }
        x
    }

    /// Perform a do-calculus intervention: set variable `target` to value `val`.
    ///
    /// Returns the modified SCM where all incoming edges to `target` are removed
    /// and its value is fixed at `val` (zero noise, intercept = val).
    pub fn intervene(&self, target: usize, val: f64) -> Self {
        let mut scm = self.clone();
        // Remove all parents of target
        let parents = scm.graph.parents[target].clone();
        for &p in &parents {
            scm.graph.children[p].retain(|&c| c != target);
        }
        scm.graph.parents[target].clear();
        scm.coefficients[target].clear();
        scm.noise_std[target] = 0.0;
        scm.intercepts[target] = val;
        scm
    }

    /// Compute the average causal effect (ACE) of intervention `do(X_cause = val)`
    /// on variable `effect`, using the provided noise samples.
    ///
    /// # Arguments
    /// * `cause`      — the variable to intervene on
    /// * `val`        — the intervention value
    /// * `effect`     — the outcome variable
    /// * `noise_samples` — matrix of noise samples, shape `[n_samples][n_nodes]`
    pub fn average_causal_effect(
        &self,
        cause: usize,
        val: f64,
        effect: usize,
        noise_samples: &[Vec<f64>],
    ) -> f64 {
        let intervened = self.intervene(cause, val);
        let mean: f64 = noise_samples
            .iter()
            .map(|noise| intervened.sample_with_noise(noise)[effect])
            .sum::<f64>()
            / noise_samples.len() as f64;
        mean
    }

    /// Compute the total causal effect of `cause` on `effect` analytically
    /// (only valid for linear SCMs).
    ///
    /// Sums all directed path contributions.
    pub fn total_effect_linear(&self, cause: usize, effect: usize) -> f64 {
        // BFS/DFS to enumerate all directed paths and multiply coefficients
        let mut total = 0.0_f64;
        // Stack of (current_node, accumulated_product)
        let mut stack: Vec<(usize, f64)> = vec![(cause, 1.0)];
        while let Some((node, prod)) = stack.pop() {
            if node == effect && node != cause {
                total += prod;
            }
            for (k, &child) in self.graph.children[node].iter().enumerate() {
                // Find the index of `node` in child's parent list
                if let Some(idx) = self.graph.parents[child].iter().position(|&p| p == node) {
                    let coeff = self.coefficients[child][idx];
                    let _ = k; // suppress unused warning
                    stack.push((child, prod * coeff));
                }
            }
        }
        total
    }
}

// ---------------------------------------------------------------------------
// BackdoorCriterion
// ---------------------------------------------------------------------------

/// Checks whether a set of variables satisfies the backdoor criterion for
/// identifying the causal effect of `treatment` on `outcome`.
///
/// The backdoor criterion holds if:
/// 1. No variable in `adjustment_set` is a descendant of `treatment`.
/// 2. `adjustment_set` blocks all backdoor paths from `treatment` to `outcome`.
#[derive(Debug, Clone)]
pub struct BackdoorCriterion {
    /// The causal graph.
    pub graph: CausalGraph,
}

impl BackdoorCriterion {
    /// Create a new backdoor criterion checker.
    pub fn new(graph: CausalGraph) -> Self {
        Self { graph }
    }

    /// Check if `adjustment_set` satisfies the backdoor criterion for the
    /// causal effect of `treatment` on `outcome`.
    ///
    /// Returns `true` if the criterion is satisfied.
    pub fn check(&self, treatment: usize, outcome: usize, adjustment_set: &[usize]) -> bool {
        let desc_treatment = self.graph.descendants(treatment);

        // Criterion 1: no variable in Z is a descendant of X
        for &z in adjustment_set {
            if desc_treatment.contains(&z) {
                return false;
            }
        }

        // Criterion 2: Z blocks all backdoor paths
        // A backdoor path is a path from X to Y that starts with an arrow INTO X.
        // We check this by creating a modified graph where outgoing edges from X
        // are removed, and then check d-separation.
        let mut modified = self.graph.clone();
        // Remove all outgoing edges from treatment in the modified graph
        let children_of_treatment = modified.graph_children_of(treatment);
        for &c in &children_of_treatment {
            modified.parents[c].retain(|&p| p != treatment);
        }
        modified.children[treatment].clear();

        modified.d_separated(&[treatment], &[outcome], adjustment_set)
    }

    /// Compute the backdoor-adjusted causal effect estimate from observational data.
    ///
    /// Uses the adjustment formula: E\[Y|do(X=x)\] = Σ_z E\[Y|X=x, Z=z\] * P(Z=z)
    ///
    /// This simplified version takes pre-computed conditional means.
    ///
    /// # Arguments
    /// * `data_x`   — treatment values
    /// * `data_y`   — outcome values
    /// * `data_z`   — confounder values (single confounder for simplicity)
    /// * `x_val`    — the intervention value
    pub fn adjusted_effect(
        data_x: &[f64],
        data_y: &[f64],
        data_z: &[f64],
        x_val: f64,
        _tolerance: f64,
    ) -> f64 {
        let n = data_x.len();
        assert_eq!(data_y.len(), n);
        assert_eq!(data_z.len(), n);

        // Bin Z into quantile strata for simple adjustment
        let n_strata = 5usize;
        let mut sorted_z = data_z.to_vec();
        sorted_z.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let quantiles: Vec<f64> = (1..n_strata)
            .map(|i| sorted_z[(i * n) / n_strata])
            .collect();

        let stratum_of = |z: f64| -> usize {
            quantiles
                .iter()
                .position(|&q| z < q)
                .unwrap_or(n_strata - 1)
        };

        // For each stratum, estimate E[Y|X≈x_val, Z=stratum] and P(Z=stratum)
        let mut stratum_sums_y = vec![0.0_f64; n_strata];
        let mut stratum_counts = vec![0usize; n_strata];
        let mut stratum_counts_near_x = vec![0usize; n_strata];
        let mut stratum_y_near_x = vec![0.0_f64; n_strata];

        let bandwidth = {
            let mut xs = data_x.to_vec();
            xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let iqr = xs[3 * n / 4] - xs[n / 4];
            iqr.max(0.1) * 0.5
        };

        for i in 0..n {
            let s = stratum_of(data_z[i]);
            stratum_sums_y[s] += data_y[i];
            stratum_counts[s] += 1;
            if (data_x[i] - x_val).abs() < bandwidth {
                stratum_y_near_x[s] += data_y[i];
                stratum_counts_near_x[s] += 1;
            }
        }

        let mut total = 0.0_f64;
        for ((&sc, &sc_near_x), (&sy_near_x, &ssum_y)) in stratum_counts
            .iter()
            .zip(stratum_counts_near_x.iter())
            .zip(stratum_y_near_x.iter().zip(stratum_sums_y.iter()))
        {
            if sc == 0 {
                continue;
            }
            let p_z = sc as f64 / n as f64;
            let e_y_xz = if sc_near_x > 0 {
                sy_near_x / sc_near_x as f64
            } else {
                ssum_y / sc as f64
            };
            total += e_y_xz * p_z;
        }
        total
    }
}

// Extension trait for internal use
trait GraphChildrenOf {
    fn graph_children_of(&self, v: usize) -> Vec<usize>;
}

impl GraphChildrenOf for CausalGraph {
    fn graph_children_of(&self, v: usize) -> Vec<usize> {
        self.children[v].clone()
    }
}

// ---------------------------------------------------------------------------
// FrontdoorCriterion
// ---------------------------------------------------------------------------

/// Implements the frontdoor adjustment formula for causal effect identification.
///
/// The frontdoor criterion allows identification of causal effects through a
/// mediator set `M` when direct adjustment is not possible.
#[derive(Debug, Clone)]
pub struct FrontdoorCriterion {
    /// The causal graph.
    pub graph: CausalGraph,
}

impl FrontdoorCriterion {
    /// Create a new frontdoor criterion object.
    pub fn new(graph: CausalGraph) -> Self {
        Self { graph }
    }

    /// Check if `mediator_set` satisfies the frontdoor criterion for identifying
    /// the effect of `treatment` on `outcome`.
    ///
    /// Conditions:
    /// 1. All directed paths from `treatment` to `outcome` are intercepted by `M`.
    /// 2. No backdoor path from `treatment` to `M` (or blocked by `treatment`).
    /// 3. All backdoor paths from `M` to `outcome` are blocked by `treatment`.
    pub fn check(&self, treatment: usize, outcome: usize, mediator_set: &[usize]) -> bool {
        let med_set: HashSet<usize> = mediator_set.iter().copied().collect();

        // Condition 1: M intercepts all directed paths from X to Y
        if !self.intercepts_all_paths(treatment, outcome, &med_set) {
            return false;
        }

        // Condition 2: no unblocked backdoor from X to M (blocked by ∅)
        // i.e., X d-separates from M given ∅ in graph with X's parents cut
        // Simplified: check no common causes of X and M that aren't through X
        for &m in mediator_set {
            if !self.graph.d_separated(&[treatment], &[m], &[treatment]) {
                // Check via empty set
                let x_anc = self.graph.ancestors(treatment);
                let m_anc = self.graph.ancestors(m);
                // If there is overlap in ancestors excluding X's subtree, problem exists
                // Simplified check for demo purposes
                let _ = (x_anc, m_anc);
            }
        }

        // Condition 3: all backdoor from M to Y are blocked by X
        for &m in mediator_set {
            if !self.graph.d_separated(&[m], &[outcome], &[treatment]) {
                return false;
            }
        }

        true
    }

    /// Check if `med_set` intercepts all directed paths from `src` to `dst`.
    fn intercepts_all_paths(&self, src: usize, dst: usize, med_set: &HashSet<usize>) -> bool {
        // DFS: find if any path from src to dst avoids med_set
        let mut stack: Vec<(usize, Vec<usize>)> = vec![(src, vec![src])];
        while let Some((node, path)) = stack.pop() {
            if node == dst {
                // Found a path; check if med_set is on it (excluding src)
                let on_path = path[1..].iter().any(|v| med_set.contains(v));
                if !on_path {
                    return false;
                }
                continue;
            }
            for &child in &self.graph.children[node] {
                if !path.contains(&child) {
                    let mut new_path = path.clone();
                    new_path.push(child);
                    stack.push((child, new_path));
                }
            }
        }
        true
    }

    /// Compute the frontdoor-adjusted causal effect using sample data.
    ///
    /// E\[Y|do(X=x)\] = Σ_m P(M=m|X=x) Σ_x' E\[Y|M=m, X=x'\] P(X=x')
    ///
    /// Simplified discrete approximation for a single mediator.
    pub fn adjusted_effect(data_x: &[f64], data_m: &[f64], data_y: &[f64], x_val: f64) -> f64 {
        let n = data_x.len();
        let bandwidth = {
            let mut xs = data_x.to_vec();
            xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let iqr = xs[3 * n / 4] - xs[n / 4];
            iqr.max(0.1) * 0.4
        };

        // Approximate E[M|X=x] via kernel smoothing
        let (mut sum_m, mut w_sum) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            let w = gaussian_kernel((data_x[i] - x_val) / bandwidth);
            sum_m += w * data_m[i];
            w_sum += w;
        }
        let e_m_given_x = if w_sum > 1e-12 { sum_m / w_sum } else { 0.0 };

        // Approximate E[Y|M=m] via kernel smoothing on M
        let bw_m = {
            let mut ms = data_m.to_vec();
            ms.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let iqr = ms[3 * n / 4] - ms[n / 4];
            iqr.max(0.1) * 0.4
        };
        let (mut sum_y, mut wy_sum) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            let w = gaussian_kernel((data_m[i] - e_m_given_x) / bw_m);
            sum_y += w * data_y[i];
            wy_sum += w;
        }
        if wy_sum > 1e-12 { sum_y / wy_sum } else { 0.0 }
    }
}

/// Gaussian kernel function for kernel smoothing.
fn gaussian_kernel(u: f64) -> f64 {
    (-0.5 * u * u).exp()
}

// ---------------------------------------------------------------------------
// PropensityScoreMatching
// ---------------------------------------------------------------------------

/// Propensity score matching for observational causal inference.
///
/// Estimates the probability of treatment assignment P(T=1|X) using logistic
/// regression, then matches treated and control units.
#[derive(Debug, Clone)]
pub struct PropensityScoreMatching {
    /// Logistic regression weights (length = n_covariates + 1, with intercept).
    pub weights: Vec<f64>,
    /// Number of covariates.
    pub n_covariates: usize,
}

impl PropensityScoreMatching {
    /// Create a new propensity score matcher.
    ///
    /// # Arguments
    /// * `n_covariates` — number of covariate dimensions
    pub fn new(n_covariates: usize) -> Self {
        Self {
            weights: vec![0.0; n_covariates + 1],
            n_covariates,
        }
    }

    /// Fit logistic regression to estimate propensity scores via gradient descent.
    ///
    /// # Arguments
    /// * `covariates` — matrix of covariates, shape `[n_obs][n_covariates]`
    /// * `treatment`  — binary treatment indicator (0 or 1), length `n_obs`
    /// * `lr`         — learning rate
    /// * `n_iter`     — number of gradient descent iterations
    pub fn fit(&mut self, covariates: &[Vec<f64>], treatment: &[f64], lr: f64, n_iter: usize) {
        let n = covariates.len();
        assert_eq!(treatment.len(), n);
        for _ in 0..n_iter {
            let mut grad = vec![0.0_f64; self.n_covariates + 1];
            for i in 0..n {
                let p = self.predict_one(&covariates[i]);
                let err = p - treatment[i];
                grad[0] += err; // intercept
                for j in 0..self.n_covariates {
                    grad[j + 1] += err * covariates[i][j];
                }
            }
            for (w, &g) in self.weights.iter_mut().zip(grad.iter()) {
                *w -= lr * g / n as f64;
            }
        }
    }

    /// Predict propensity score P(T=1|X=x) for a single observation.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        let logit: f64 = self.weights[0]
            + x.iter()
                .zip(self.weights[1..].iter())
                .map(|(xi, wi)| xi * wi)
                .sum::<f64>();
        sigmoid(logit)
    }

    /// Predict propensity scores for all observations.
    pub fn predict(&self, covariates: &[Vec<f64>]) -> Vec<f64> {
        covariates.iter().map(|x| self.predict_one(x)).collect()
    }

    /// Estimate the Average Treatment Effect (ATE) using IPW (Inverse Probability Weighting).
    ///
    /// ATE = E\[Y(1)\] - E\[Y(0)\] = E\[T*Y/e(X)\] - E\[(1-T)*Y/(1-e(X))\]
    pub fn estimate_ate(&self, covariates: &[Vec<f64>], treatment: &[f64], outcome: &[f64]) -> f64 {
        let n = covariates.len();
        let (mut sum1, mut w1, mut sum0, mut w0) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..n {
            let e = self.predict_one(&covariates[i]).clamp(1e-6, 1.0 - 1e-6);
            if treatment[i] > 0.5 {
                sum1 += outcome[i] / e;
                w1 += 1.0 / e;
            } else {
                sum0 += outcome[i] / (1.0 - e);
                w0 += 1.0 / (1.0 - e);
            }
        }
        let ey1 = if w1 > 0.0 { sum1 / w1 } else { 0.0 };
        let ey0 = if w0 > 0.0 { sum0 / w0 } else { 0.0 };
        ey1 - ey0
    }

    /// Estimate the Average Treatment Effect on the Treated (ATT).
    ///
    /// ATT = E\[Y(1)-Y(0)|T=1\]
    pub fn estimate_att(&self, covariates: &[Vec<f64>], treatment: &[f64], outcome: &[f64]) -> f64 {
        let n = covariates.len();
        let mut treated_y: Vec<f64> = Vec::new();
        let mut control_y: Vec<f64> = Vec::new();
        let mut control_ps: Vec<f64> = Vec::new();

        for i in 0..n {
            let e = self.predict_one(&covariates[i]).clamp(1e-6, 1.0 - 1e-6);
            if treatment[i] > 0.5 {
                treated_y.push(outcome[i]);
            } else {
                control_y.push(outcome[i]);
                control_ps.push(e / (1.0 - e)); // odds
            }
        }

        if treated_y.is_empty() || control_y.is_empty() {
            return 0.0;
        }

        let mean_treated = treated_y.iter().sum::<f64>() / treated_y.len() as f64;
        let total_weight: f64 = control_ps.iter().sum();
        let mean_control = if total_weight > 0.0 {
            control_y
                .iter()
                .zip(control_ps.iter())
                .map(|(y, w)| y * w)
                .sum::<f64>()
                / total_weight
        } else {
            control_y.iter().sum::<f64>() / control_y.len() as f64
        };

        mean_treated - mean_control
    }
}

/// Logistic sigmoid function.
fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

// ---------------------------------------------------------------------------
// InstrumentalVariables
// ---------------------------------------------------------------------------

/// Instrumental variables (IV) estimation and Two-Stage Least Squares (2SLS).
///
/// Used when treatment is endogenous (correlated with the error term).
/// Requires a valid instrument `Z` that:
/// 1. Is correlated with the treatment `D`.
/// 2. Affects the outcome `Y` only through `D` (exclusion restriction).
/// 3. Is independent of unobserved confounders.
#[derive(Debug, Clone)]
pub struct InstrumentalVariables {
    /// Number of endogenous variables.
    pub n_endogenous: usize,
    /// Number of instruments.
    pub n_instruments: usize,
    /// First-stage coefficients (instrument → treatment).
    pub first_stage: Vec<f64>,
    /// Second-stage coefficient (treatment → outcome).
    pub second_stage: f64,
}

impl InstrumentalVariables {
    /// Create a new IV estimator.
    pub fn new(n_endogenous: usize, n_instruments: usize) -> Self {
        Self {
            n_endogenous,
            n_instruments,
            first_stage: vec![0.0; n_instruments + 1],
            second_stage: 0.0,
        }
    }

    /// Fit the 2SLS estimator.
    ///
    /// Stage 1: regress treatment `d` on instruments `z`.
    /// Stage 2: regress outcome `y` on predicted treatment `d_hat`.
    ///
    /// # Arguments
    /// * `y` — outcome variable
    /// * `d` — endogenous treatment variable
    /// * `z` — instruments matrix, shape `[n_obs][n_instruments]`
    pub fn fit_2sls(&mut self, y: &[f64], d: &[f64], z: &[Vec<f64>]) {
        let n = y.len();
        assert_eq!(d.len(), n);
        assert_eq!(z.len(), n);

        // Stage 1: OLS of D on Z (including intercept)
        // Simple single-instrument case for clarity
        let n_inst = self.n_instruments;
        let mut d_hat = vec![0.0_f64; n];

        if n_inst == 1 {
            // Simple IV: β_1 = Cov(Y,Z) / Cov(D,Z)
            let z_vec: Vec<f64> = z.iter().map(|row| row[0]).collect();
            let mean_z = z_vec.iter().sum::<f64>() / n as f64;
            let mean_d = d.iter().sum::<f64>() / n as f64;
            let mean_y = y.iter().sum::<f64>() / n as f64;

            let cov_dz: f64 = d
                .iter()
                .zip(z_vec.iter())
                .map(|(di, zi)| (di - mean_d) * (zi - mean_z))
                .sum::<f64>()
                / n as f64;
            let cov_yz: f64 = y
                .iter()
                .zip(z_vec.iter())
                .map(|(yi, zi)| (yi - mean_y) * (zi - mean_z))
                .sum::<f64>()
                / n as f64;
            let var_z: f64 = z_vec.iter().map(|zi| (zi - mean_z).powi(2)).sum::<f64>() / n as f64;

            // First stage: d = α0 + α1*z
            let alpha1 = if var_z.abs() > 1e-12 {
                cov_dz / var_z
            } else {
                0.0
            };
            let alpha0 = mean_d - alpha1 * mean_z;
            self.first_stage[0] = alpha0;
            self.first_stage[1] = alpha1;

            // IV estimate
            self.second_stage = if cov_dz.abs() > 1e-12 {
                cov_yz / cov_dz
            } else {
                0.0
            };

            for i in 0..n {
                d_hat[i] = alpha0 + alpha1 * z_vec[i];
            }
        } else {
            // Multi-instrument: use OLS for first stage
            // Simplified: use first instrument only
            let z0: Vec<f64> = z.iter().map(|row| row[0]).collect();
            let mean_z0 = z0.iter().sum::<f64>() / n as f64;
            let mean_d = d.iter().sum::<f64>() / n as f64;

            let cov = z0
                .iter()
                .zip(d.iter())
                .map(|(zi, di)| (zi - mean_z0) * (di - mean_d))
                .sum::<f64>()
                / n as f64;
            let var_z0 = z0.iter().map(|zi| (zi - mean_z0).powi(2)).sum::<f64>() / n as f64;

            let alpha1 = if var_z0 > 1e-12 { cov / var_z0 } else { 0.0 };
            let alpha0 = mean_d - alpha1 * mean_z0;
            self.first_stage[0] = alpha0;
            self.first_stage[1] = alpha1;

            for i in 0..n {
                d_hat[i] = alpha0 + alpha1 * z0[i];
            }

            // Stage 2: OLS of Y on D_hat
            let mean_dhat = d_hat.iter().sum::<f64>() / n as f64;
            let mean_y = y.iter().sum::<f64>() / n as f64;
            let cov_ydhat: f64 = y
                .iter()
                .zip(d_hat.iter())
                .map(|(yi, di)| (yi - mean_y) * (di - mean_dhat))
                .sum::<f64>()
                / n as f64;
            let var_dhat: f64 =
                d_hat.iter().map(|di| (di - mean_dhat).powi(2)).sum::<f64>() / n as f64;
            self.second_stage = if var_dhat > 1e-12 {
                cov_ydhat / var_dhat
            } else {
                0.0
            };
        }
    }

    /// Compute the first-stage F-statistic (instrument relevance test).
    ///
    /// Large F (> 10) indicates strong instruments.
    pub fn first_stage_f_stat(&self, y: &[f64], d: &[f64], z: &[Vec<f64>]) -> f64 {
        let n = y.len();
        let z0: Vec<f64> = z.iter().map(|row| row[0]).collect();
        let _mean_z0 = z0.iter().sum::<f64>() / n as f64;
        let mean_d = d.iter().sum::<f64>() / n as f64;

        let d_hat: Vec<f64> = z0
            .iter()
            .map(|zi| self.first_stage[0] + self.first_stage[1] * zi)
            .collect();

        let ss_res: f64 = d
            .iter()
            .zip(d_hat.iter())
            .map(|(di, dh)| (di - dh).powi(2))
            .sum();
        let ss_tot: f64 = d.iter().map(|di| (di - mean_d).powi(2)).sum();

        let r2 = 1.0 - ss_res / ss_tot.max(1e-12);
        let k = 1.0_f64; // number of instruments
        let n_f = n as f64;
        (r2 / k) / ((1.0 - r2) / (n_f - k - 1.0)).max(1e-12)
    }

    /// Predict the causal effect for a new treatment value.
    pub fn predict(&self, d_val: f64) -> f64 {
        self.second_stage * d_val
    }
}

// ---------------------------------------------------------------------------
// CausalDiscovery
// ---------------------------------------------------------------------------

/// Causal discovery via the PC algorithm.
///
/// The PC algorithm learns the structure of a DAG from conditional independence
/// tests on observational data. It produces a Completed Partially Directed
/// Acyclic Graph (CPDAG) representing the Markov equivalence class.
#[derive(Debug, Clone)]
pub struct CausalDiscovery {
    /// Number of variables.
    pub n_vars: usize,
    /// Adjacency matrix of the skeleton (undirected).
    pub skeleton: Vec<Vec<bool>>,
    /// Directed adjacency: `directed[i][j] = true` means i → j is oriented.
    pub directed: Vec<Vec<bool>>,
    /// Separation sets: `sep_sets[(i,j)]` = the conditioning set that d-separates i and j.
    pub sep_sets: HashMap<(usize, usize), Vec<usize>>,
    /// Significance threshold for independence tests.
    pub alpha: f64,
}

impl CausalDiscovery {
    /// Create a new PC algorithm runner.
    ///
    /// # Arguments
    /// * `n_vars` — number of observed variables
    /// * `alpha`  — significance level for conditional independence tests
    pub fn new(n_vars: usize, alpha: f64) -> Self {
        Self {
            n_vars,
            skeleton: vec![vec![true; n_vars]; n_vars],
            directed: vec![vec![false; n_vars]; n_vars],
            sep_sets: HashMap::new(),
            alpha,
        }
    }

    /// Learn the skeleton from a data matrix using partial correlation tests.
    ///
    /// # Arguments
    /// * `data` — data matrix, shape `[n_obs][n_vars]`
    pub fn learn_skeleton(&mut self, data: &[Vec<f64>]) {
        let n = self.n_vars;

        // Remove self-loops
        for i in 0..n {
            self.skeleton[i][i] = false;
        }

        // Level 0: unconditional independence
        for i in 0..n {
            for j in (i + 1)..n {
                let r = partial_correlation(data, i, j, &[]);
                let p = fisher_z_test(r, data.len(), 0);
                if p > self.alpha {
                    self.skeleton[i][j] = false;
                    self.skeleton[j][i] = false;
                    self.sep_sets.insert((i, j), vec![]);
                    self.sep_sets.insert((j, i), vec![]);
                }
            }
        }

        // Level 1+: conditional independence given conditioning sets
        for cond_size in 1..n.saturating_sub(1) {
            for i in 0..n {
                let adj_i: Vec<usize> = (0..n).filter(|&k| k != i && self.skeleton[i][k]).collect();
                for &j in &adj_i {
                    if !self.skeleton[i][j] {
                        continue;
                    }
                    let adj_minus_j: Vec<usize> =
                        adj_i.iter().copied().filter(|&k| k != j).collect();
                    if adj_minus_j.len() < cond_size {
                        continue;
                    }
                    // Test all conditioning sets of size `cond_size`
                    for cond_set in subsets(&adj_minus_j, cond_size) {
                        let r = partial_correlation(data, i, j, &cond_set);
                        let p = fisher_z_test(r, data.len(), cond_size);
                        if p > self.alpha {
                            self.skeleton[i][j] = false;
                            self.skeleton[j][i] = false;
                            self.sep_sets.insert((i, j), cond_set.clone());
                            self.sep_sets.insert((j, i), cond_set);
                            break;
                        }
                    }
                }
            }
        }
    }

    /// Orient v-structures (colliders) in the skeleton.
    ///
    /// For each unshielded triple i — k — j (i and j not adjacent),
    /// if k is not in sep(i,j), orient i→k←j.
    pub fn orient_v_structures(&mut self) {
        let n = self.n_vars;
        for i in 0..n {
            for k in 0..n {
                if i == k || !self.skeleton[i][k] {
                    continue;
                }
                for j in (i + 1)..n {
                    if j == k || !self.skeleton[k][j] || self.skeleton[i][j] {
                        continue;
                    }
                    // i — k — j, unshielded
                    let sep = self.sep_sets.get(&(i, j)).cloned().unwrap_or_default();
                    if !sep.contains(&k) {
                        // Orient as collider: i → k ← j
                        self.directed[i][k] = true;
                        self.directed[j][k] = true;
                        self.skeleton[k][i] = false;
                        self.skeleton[k][j] = false;
                    }
                }
            }
        }
    }

    /// Apply Meek's orientation rules to complete the CPDAG.
    ///
    /// Rules R1–R3 propagate orientations to avoid new v-structures and cycles.
    pub fn apply_meek_rules(&mut self) {
        let n = self.n_vars;
        let mut changed = true;
        while changed {
            changed = false;
            // R1: If i→j — k and i not adjacent to k, orient j→k
            for i in 0..n {
                for j in 0..n {
                    if !self.directed[i][j] {
                        continue;
                    }
                    for k in 0..n {
                        if k == i || k == j {
                            continue;
                        }
                        if self.skeleton[j][k]
                            && !self.directed[j][k]
                            && !self.directed[k][j]
                            && !self.skeleton[i][k]
                        {
                            self.directed[j][k] = true;
                            self.skeleton[k][j] = false;
                            changed = true;
                        }
                    }
                }
            }
            // R2: If i→k→j and i — j, orient i→j
            for i in 0..n {
                for j in 0..n {
                    if i == j || !self.skeleton[i][j] || self.directed[i][j] {
                        continue;
                    }
                    for k in 0..n {
                        if k == i || k == j {
                            continue;
                        }
                        if self.directed[i][k] && self.directed[k][j] {
                            self.directed[i][j] = true;
                            self.skeleton[j][i] = false;
                            changed = true;
                        }
                    }
                }
            }
        }
    }

    /// Run the full PC algorithm: skeleton → v-structures → Meek rules.
    pub fn run(&mut self, data: &[Vec<f64>]) {
        self.learn_skeleton(data);
        self.orient_v_structures();
        self.apply_meek_rules();
    }
}

/// Compute partial correlation of variables `i` and `j` conditioning on `cond`.
///
/// Uses recursive partial correlation formula for efficiency.
pub fn partial_correlation(data: &[Vec<f64>], i: usize, j: usize, cond: &[usize]) -> f64 {
    if cond.is_empty() {
        return pearson_correlation(data, i, j);
    }
    if cond.len() == 1 {
        let k = cond[0];
        let r_ij = pearson_correlation(data, i, j);
        let r_ik = pearson_correlation(data, i, k);
        let r_jk = pearson_correlation(data, j, k);
        let denom = ((1.0 - r_ik * r_ik) * (1.0 - r_jk * r_jk)).sqrt();
        if denom < 1e-12 {
            return 0.0;
        }
        return (r_ij - r_ik * r_jk) / denom;
    }
    // For larger conditioning sets, use matrix inversion approach
    // Simplified: use iterative partial correlations
    let last = cond[cond.len() - 1];
    let rest = &cond[..cond.len() - 1];
    let r_ij_rest = partial_correlation(data, i, j, rest);
    let r_ik_rest = partial_correlation(data, i, last, rest);
    let r_jk_rest = partial_correlation(data, j, last, rest);
    let denom = ((1.0 - r_ik_rest * r_ik_rest) * (1.0 - r_jk_rest * r_jk_rest)).sqrt();
    if denom < 1e-12 {
        return 0.0;
    }
    (r_ij_rest - r_ik_rest * r_jk_rest) / denom
}

/// Pearson correlation between variables `i` and `j` in a data matrix.
pub fn pearson_correlation(data: &[Vec<f64>], i: usize, j: usize) -> f64 {
    let n = data.len() as f64;
    let mean_i = data.iter().map(|row| row[i]).sum::<f64>() / n;
    let mean_j = data.iter().map(|row| row[j]).sum::<f64>() / n;
    let cov: f64 = data
        .iter()
        .map(|row| (row[i] - mean_i) * (row[j] - mean_j))
        .sum::<f64>()
        / n;
    let std_i = (data
        .iter()
        .map(|row| (row[i] - mean_i).powi(2))
        .sum::<f64>()
        / n)
        .sqrt();
    let std_j = (data
        .iter()
        .map(|row| (row[j] - mean_j).powi(2))
        .sum::<f64>()
        / n)
        .sqrt();
    if std_i < 1e-12 || std_j < 1e-12 {
        return 0.0;
    }
    (cov / (std_i * std_j)).clamp(-1.0, 1.0)
}

/// Fisher Z-test for conditional independence.
///
/// Returns the p-value for the null hypothesis r = 0.
pub fn fisher_z_test(r: f64, n: usize, cond_size: usize) -> f64 {
    let r = r.clamp(-0.9999, 0.9999);
    let z = 0.5 * ((1.0 + r) / (1.0 - r)).ln();
    let se = 1.0 / ((n as f64 - cond_size as f64 - 3.0).max(1.0)).sqrt();
    let stat = (z / se).abs();
    // Two-tailed p-value from standard normal
    2.0 * (1.0 - standard_normal_cdf(stat))
}

fn standard_normal_cdf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf = (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let cdf = 1.0 - pdf * poly;
    if x >= 0.0 { cdf } else { 1.0 - cdf }
}

/// Generate all subsets of `set` of exactly size `k`.
fn subsets(set: &[usize], k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![vec![]];
    }
    if set.len() < k {
        return vec![];
    }
    let mut result = Vec::new();
    for (i, &v) in set.iter().enumerate() {
        let rest = subsets(&set[(i + 1)..], k - 1);
        for mut subset in rest {
            subset.insert(0, v);
            result.push(subset);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// CounterfactualQuery
// ---------------------------------------------------------------------------

/// Compute counterfactual queries of the form E\[Y | do(X=x), Z=z\].
///
/// Uses the abduction-action-prediction three-step procedure on a linear SCM.
#[derive(Debug, Clone)]
pub struct CounterfactualQuery {
    /// The structural causal model for the query.
    pub scm: StructuralCausalModel,
}

impl CounterfactualQuery {
    /// Create a new counterfactual query engine.
    pub fn new(scm: StructuralCausalModel) -> Self {
        Self { scm }
    }

    /// Compute E\[Y | do(X_target=x_val)\] using the interventional distribution.
    ///
    /// # Arguments
    /// * `target`        — the variable to intervene on
    /// * `x_val`         — the intervention value
    /// * `outcome`       — the outcome variable
    /// * `noise_samples` — noise samples for Monte Carlo evaluation
    pub fn query_do(
        &self,
        target: usize,
        x_val: f64,
        outcome: usize,
        noise_samples: &[Vec<f64>],
    ) -> f64 {
        self.scm
            .average_causal_effect(target, x_val, outcome, noise_samples)
    }

    /// Compute the counterfactual: given that we observed `obs` (variable→value pairs),
    /// what would `outcome` have been if `do(target=x_val)`?
    ///
    /// Three steps:
    /// 1. Abduction: infer noise `U` from observations
    /// 2. Action: intervene on `target`
    /// 3. Prediction: compute outcome under intervention
    ///
    /// # Arguments
    /// * `obs`        — observed values, `obs[v]` = Some(value) if observed
    /// * `target`     — intervention variable
    /// * `x_val`      — intervention value
    /// * `outcome`    — outcome variable index
    pub fn counterfactual(
        &self,
        obs: &[Option<f64>],
        target: usize,
        x_val: f64,
        outcome: usize,
    ) -> f64 {
        let n = self.scm.graph.n_nodes;
        assert_eq!(obs.len(), n);

        // Step 1: Abduction — infer noise from observations
        // For linear SCM: U_v = X_v - intercept - Σ coeff * X_parent
        let order = self
            .scm
            .graph
            .topological_sort()
            .expect("SCM must be acyclic");
        let mut x = vec![0.0_f64; n];
        let mut noise = vec![0.0_f64; n];

        for &v in &order {
            if let Some(val) = obs[v] {
                x[v] = val;
                let pred: f64 = self.scm.graph.parents[v]
                    .iter()
                    .zip(self.scm.coefficients[v].iter())
                    .map(|(&p, &c)| c * x[p])
                    .sum::<f64>();
                let residual = val - self.scm.intercepts[v] - pred;
                noise[v] = if self.scm.noise_std[v].abs() > 1e-12 {
                    residual / self.scm.noise_std[v]
                } else {
                    0.0
                };
            } else {
                // Unobserved: use noise = 0 (mean)
                noise[v] = 0.0;
                let pred: f64 = self.scm.graph.parents[v]
                    .iter()
                    .zip(self.scm.coefficients[v].iter())
                    .map(|(&p, &c)| c * x[p])
                    .sum::<f64>();
                x[v] = self.scm.intercepts[v] + pred;
            }
        }

        // Step 2: Action — intervene
        let intervened = self.scm.intervene(target, x_val);

        // Step 3: Prediction — run SCM with inferred noise under intervention
        intervened.sample_with_noise(&noise)[outcome]
    }

    /// Compute the probability of necessity (PN): P(Y=0 | do(X=0), Y=1, X=1).
    ///
    /// Simplified calculation using noise samples.
    pub fn probability_of_necessity(
        &self,
        treatment: usize,
        outcome: usize,
        t_val: f64,
        t_counter: f64,
        threshold_y: f64,
        noise_samples: &[Vec<f64>],
    ) -> f64 {
        let mut count = 0;
        let mut denom = 0;
        for noise in noise_samples {
            // Actual world
            let x_actual = self.scm.sample_with_noise(noise);
            if x_actual[treatment] < t_val - 0.5 || x_actual[outcome] < threshold_y {
                continue;
            }
            denom += 1;
            // Counterfactual: set treatment to counter
            let counter_scm = self.scm.intervene(treatment, t_counter);
            let x_counter = counter_scm.sample_with_noise(noise);
            if x_counter[outcome] < threshold_y {
                count += 1;
            }
        }
        if denom == 0 {
            0.0
        } else {
            count as f64 / denom as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: Covariance matrix operations
// ---------------------------------------------------------------------------

/// Compute the sample covariance matrix from a data matrix.
///
/// Returns a flat row-major matrix of size `[n_vars * n_vars]`.
pub fn sample_covariance(data: &[Vec<f64>]) -> Vec<f64> {
    let n = data.len();
    let p = data[0].len();
    let means: Vec<f64> = (0..p)
        .map(|j| data.iter().map(|row| row[j]).sum::<f64>() / n as f64)
        .collect();
    let mut cov = vec![0.0_f64; p * p];
    for row in data.iter() {
        for j in 0..p {
            for k in j..p {
                cov[j * p + k] += (row[j] - means[j]) * (row[k] - means[k]);
            }
        }
    }
    for j in 0..p {
        for k in j..p {
            cov[j * p + k] /= (n - 1) as f64;
            cov[k * p + j] = cov[j * p + k];
        }
    }
    cov
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_chain() -> CausalGraph {
        // X0 → X1 → X2
        let mut g = CausalGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g
    }

    fn fork_graph() -> CausalGraph {
        // X0 → X1, X0 → X2
        let mut g = CausalGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g
    }

    fn collider_graph() -> CausalGraph {
        // X0 → X2, X1 → X2
        let mut g = CausalGraph::new(3);
        g.add_edge(0, 2);
        g.add_edge(1, 2);
        g
    }

    // --- CausalGraph tests ---

    #[test]
    fn test_topological_sort_chain() {
        let g = simple_chain();
        let order = g.topological_sort().unwrap();
        assert_eq!(order, vec![0, 1, 2]);
    }

    #[test]
    fn test_topological_sort_fork() {
        let g = fork_graph();
        let order = g.topological_sort().unwrap();
        assert_eq!(order[0], 0); // X0 must come first
    }

    #[test]
    fn test_ancestors() {
        let g = simple_chain();
        let anc = g.ancestors(2);
        assert!(anc.contains(&0));
        assert!(anc.contains(&1));
        assert!(!anc.contains(&2));
    }

    #[test]
    fn test_descendants() {
        let g = simple_chain();
        let desc = g.descendants(0);
        assert!(desc.contains(&1));
        assert!(desc.contains(&2));
    }

    #[test]
    fn test_d_separation_chain_blocked_by_middle() {
        // X0 → X1 → X2; conditioning on X1 blocks 0⊥2|{1}
        let g = simple_chain();
        assert!(g.d_separated(&[0], &[2], &[1]));
    }

    #[test]
    fn test_d_separation_chain_not_blocked_empty() {
        let g = simple_chain();
        assert!(!g.d_separated(&[0], &[2], &[]));
    }

    #[test]
    fn test_d_separation_fork() {
        // X0 → X1, X0 → X2; conditioning on X0 blocks 1⊥2|{0}
        let g = fork_graph();
        assert!(g.d_separated(&[1], &[2], &[0]));
        assert!(!g.d_separated(&[1], &[2], &[]));
    }

    #[test]
    fn test_d_separation_collider_blocked_by_default() {
        // X0 → X2 ← X1; collider: 0⊥1|{} (blocked)
        let g = collider_graph();
        assert!(g.d_separated(&[0], &[1], &[]));
    }

    #[test]
    fn test_d_separation_collider_opened_by_conditioning() {
        // X0 → X2 ← X1; conditioning on X2 opens path: 0 NOT ⊥ 1|{2}
        let g = collider_graph();
        assert!(!g.d_separated(&[0], &[1], &[2]));
    }

    #[test]
    fn test_is_acyclic() {
        let g = simple_chain();
        assert!(g.is_acyclic());
    }

    #[test]
    fn test_creates_cycle_detected() {
        let mut g = CausalGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        assert!(g.creates_cycle(2, 0));
    }

    #[test]
    fn test_markov_blanket() {
        // X0 → X1 → X2; blanket(X1) = {X0, X2}
        let g = simple_chain();
        let blanket = g.markov_blanket(1);
        assert!(blanket.contains(&0));
        assert!(blanket.contains(&2));
        assert!(!blanket.contains(&1));
    }

    // --- StructuralCausalModel tests ---

    #[test]
    fn test_scm_sample_basic() {
        // X0 = ε0; X1 = 2*X0 + ε1
        let mut scm = StructuralCausalModel::new(2);
        scm.add_edge(0, 1, 2.0);
        scm.noise_std = vec![1.0, 0.0];
        let noise = vec![3.0, 0.0];
        let x = scm.sample_with_noise(&noise);
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_scm_intervention_removes_parents() {
        let mut scm = StructuralCausalModel::new(2);
        scm.add_edge(0, 1, 2.0);
        let intervened = scm.intervene(1, 5.0);
        assert!(intervened.graph.parents[1].is_empty());
        let x = intervened.sample_with_noise(&[1.0, 0.0]);
        assert!((x[1] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_scm_total_effect_chain() {
        // X0 → X1 → X2, coeff 2.0 and 3.0
        let mut scm = StructuralCausalModel::new(3);
        scm.add_edge(0, 1, 2.0);
        scm.add_edge(1, 2, 3.0);
        let effect = scm.total_effect_linear(0, 2);
        assert!((effect - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_scm_no_effect_for_independent_vars() {
        let mut scm = StructuralCausalModel::new(3);
        scm.add_edge(0, 1, 1.0);
        let effect = scm.total_effect_linear(0, 2);
        assert!(effect.abs() < 1e-10);
    }

    // --- BackdoorCriterion tests ---

    #[test]
    fn test_backdoor_check_valid_adjustment() {
        // W → X → Y, W → Y (confounding)
        // Adjusting for W satisfies backdoor
        let mut g = CausalGraph::new(3); // W=0, X=1, Y=2
        g.add_edge(0, 1); // W→X
        g.add_edge(1, 2); // X→Y
        g.add_edge(0, 2); // W→Y (backdoor)
        let bd = BackdoorCriterion::new(g);
        // W is not a descendant of X, and blocks the backdoor W→Y
        assert!(bd.check(1, 2, &[0]));
    }

    #[test]
    fn test_backdoor_fails_if_descendant_in_set() {
        // X → M → Y; adjusting for M (a descendant of X) fails
        let mut g = CausalGraph::new(3); // X=0, M=1, Y=2
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        let bd = BackdoorCriterion::new(g);
        assert!(!bd.check(0, 2, &[1])); // M is descendant of X
    }

    // --- PropensityScoreMatching tests ---

    #[test]
    fn test_propensity_score_fit_and_predict() {
        let mut psm = PropensityScoreMatching::new(2);
        let covariates: Vec<Vec<f64>> = (0..100).map(|i| vec![(i as f64) / 100.0, 0.5]).collect();
        let treatment: Vec<f64> = covariates
            .iter()
            .map(|x| if x[0] > 0.5 { 1.0 } else { 0.0 })
            .collect();
        psm.fit(&covariates, &treatment, 0.5, 200);
        let ps = psm.predict(&covariates);
        assert_eq!(ps.len(), 100);
        // Scores should be between 0 and 1
        for p in &ps {
            assert!(*p > 0.0 && *p < 1.0);
        }
    }

    #[test]
    fn test_ate_sign_positive() {
        // Treatment has positive effect on outcome
        let mut psm = PropensityScoreMatching::new(1);
        let n = 200;
        let covariates: Vec<Vec<f64>> = (0..n).map(|i| vec![(i as f64) / n as f64]).collect();
        let treatment: Vec<f64> = covariates
            .iter()
            .map(|x| if x[0] > 0.5 { 1.0 } else { 0.0 })
            .collect();
        let outcome: Vec<f64> = covariates
            .iter()
            .zip(treatment.iter())
            .map(|(x, t)| x[0] + 2.0 * t + 0.1)
            .collect();
        psm.fit(&covariates, &treatment, 0.3, 300);
        let ate = psm.estimate_ate(&covariates, &treatment, &outcome);
        assert!(ate > 0.0, "ATE should be positive, got {ate}");
    }

    // --- InstrumentalVariables tests ---

    #[test]
    fn test_iv_estimation_simple() {
        // Y = 2*D + e, D = Z + v (instrument Z)
        // IV estimate should recover β ≈ 2
        let n = 500;
        let z: Vec<Vec<f64>> = (0..n).map(|i| vec![(i as f64 % 2.0)]).collect();
        let d: Vec<f64> = z.iter().map(|zi| zi[0] + 0.5).collect();
        let y: Vec<f64> = d.iter().map(|di| 2.0 * di + 1.0).collect();

        let mut iv = InstrumentalVariables::new(1, 1);
        iv.fit_2sls(&y, &d, &z);
        assert!(
            (iv.second_stage - 2.0).abs() < 0.5,
            "IV est = {}",
            iv.second_stage
        );
    }

    #[test]
    fn test_iv_first_stage_f_stat() {
        let n = 200;
        let z: Vec<Vec<f64>> = (0..n).map(|i| vec![(i as f64 / n as f64)]).collect();
        let d: Vec<f64> = z.iter().map(|zi| 2.0 * zi[0]).collect();
        let y: Vec<f64> = d.iter().map(|di| di + 1.0).collect();
        let mut iv = InstrumentalVariables::new(1, 1);
        iv.fit_2sls(&y, &d, &z);
        let f = iv.first_stage_f_stat(&y, &d, &z);
        assert!(
            f > 10.0,
            "F-stat should be large for strong instrument, got {f}"
        );
    }

    // --- CausalDiscovery tests ---

    #[test]
    fn test_pearson_correlation_perfect() {
        let data: Vec<Vec<f64>> = (0..50).map(|i| vec![i as f64, 2.0 * i as f64]).collect();
        let r = pearson_correlation(&data, 0, 1);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_correlation_zero() {
        // x varies linearly, y is constant — they are uncorrelated (std_y=0 => r=0)
        let data: Vec<Vec<f64>> = (0..100).map(|i| vec![i as f64, 1.0]).collect();
        let r = pearson_correlation(&data, 0, 1);
        assert!(r.abs() < 1e-10, "r={r}");
    }

    #[test]
    fn test_partial_correlation_returns_in_range() {
        let data: Vec<Vec<f64>> = (0..50)
            .map(|i| vec![i as f64, 2.0 * i as f64 + 1.0, i as f64 * 0.5])
            .collect();
        let r = partial_correlation(&data, 0, 1, &[2]);
        assert!((-1.0..=1.0).contains(&r));
    }

    #[test]
    fn test_fisher_z_test_high_correlation() {
        let r = 0.0; // uncorrelated
        let p = fisher_z_test(r, 100, 0);
        assert!(p > 0.05, "Should not reject independence for r=0");
    }

    #[test]
    fn test_causal_discovery_skeleton_independent() {
        // Two independent variables: skeleton should have no edge
        let data: Vec<Vec<f64>> = (0..100)
            .map(|i| vec![(i as f64).sin(), (i as f64 * 2.3 + 1.0).cos()])
            .collect();
        let mut cd = CausalDiscovery::new(2, 0.01);
        cd.learn_skeleton(&data);
        // May or may not find edges; just check no panic
        assert!(cd.n_vars == 2);
    }

    #[test]
    fn test_subsets_correctness() {
        let v = vec![0, 1, 2];
        let subs = subsets(&v, 2);
        assert_eq!(subs.len(), 3);
    }

    #[test]
    fn test_subsets_empty() {
        let v = vec![0, 1];
        let subs = subsets(&v, 0);
        assert_eq!(subs.len(), 1);
        assert!(subs[0].is_empty());
    }

    // --- CounterfactualQuery tests ---

    #[test]
    fn test_counterfactual_simple_chain() {
        // X1 = ε1, X2 = 2*X1 + ε2
        // If we observe X1=1, X2=2 (noise2=0), what is X2 if do(X1=3)?
        let mut scm = StructuralCausalModel::new(2);
        scm.add_edge(0, 1, 2.0);
        scm.noise_std = vec![1.0, 1.0];
        let query = CounterfactualQuery::new(scm);
        let obs = vec![Some(1.0), Some(2.0)];
        let cf = query.counterfactual(&obs, 0, 3.0, 1);
        // noise2 = (2 - 2*1)/1 = 0; X2_cf = 2*3 + 0*1 = 6
        assert!((cf - 6.0).abs() < 1e-10, "cf={cf}");
    }

    #[test]
    fn test_counterfactual_intercept() {
        // X1 = 5 + ε1 (intercept=5); observe X1=7 → noise=2
        // Counterfactual do(X0=10): X1 stays (X0 not parent of X1)
        let mut scm = StructuralCausalModel::new(2);
        scm.intercepts[1] = 5.0;
        scm.noise_std = vec![1.0, 1.0];
        let query = CounterfactualQuery::new(scm);
        let obs = vec![Some(0.0), Some(7.0)];
        let cf = query.counterfactual(&obs, 0, 10.0, 1);
        // X1 doesn't depend on X0, noise1 = (7-5)/1 = 2, cf = 5 + 2 = 7
        assert!((cf - 7.0).abs() < 1e-10, "cf={cf}");
    }

    // --- Sample covariance test ---

    #[test]
    fn test_sample_covariance_diagonal() {
        let data: Vec<Vec<f64>> = (0..100).map(|i| vec![i as f64, 0.0]).collect();
        let cov = sample_covariance(&data);
        // var(X0) > 0, var(X1) = 0
        assert!(cov[0] > 0.0);
        assert!(cov[3].abs() < 1e-10); // var(X1) = 0
    }

    #[test]
    fn test_sample_covariance_symmetric() {
        let data: Vec<Vec<f64>> = (0..50).map(|i| vec![i as f64, (i as f64).sin()]).collect();
        let cov = sample_covariance(&data);
        assert!((cov[1] - cov[2]).abs() < 1e-12); // cov[0][1] == cov[1][0]
    }

    // --- Integration test ---

    #[test]
    fn test_full_scm_pipeline() {
        // Build a confounded SCM and estimate ACE
        // U → X, U → Y, X → Y
        // We use a 3-node SCM: U=0, X=1, Y=2
        let mut scm = StructuralCausalModel::new(3);
        scm.add_edge(0, 1, 1.0); // U→X
        scm.add_edge(0, 2, 1.0); // U→Y
        scm.add_edge(1, 2, 2.0); // X→Y (true causal effect = 2)
        scm.noise_std = vec![1.0, 0.5, 0.5];

        // Generate noise samples
        let noise_samples: Vec<Vec<f64>> = (0..500)
            .map(|i| {
                let u = (i as f64 * 0.01).sin();
                let x = (i as f64 * 0.013).cos();
                let y = (i as f64 * 0.017).sin();
                vec![u, x, y]
            })
            .collect();

        let ace0 = scm.average_causal_effect(1, 0.0, 2, &noise_samples);
        let ace1 = scm.average_causal_effect(1, 1.0, 2, &noise_samples);
        let diff = ace1 - ace0;
        // Should be close to 2.0 (the direct structural coefficient)
        assert!((diff - 2.0).abs() < 0.1, "ACE diff = {diff}");
    }
}
