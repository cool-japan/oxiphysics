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

use nalgebra::{DMatrix, DVector, Matrix2, SVD, Vector2};

use crate::error::PhysicsError;

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

    /// Return a copy of the graph with every edge **leaving** `v` removed.
    ///
    /// Every path that starts at `v` in the result necessarily begins with an
    /// arrow *into* `v`, so these are exactly the back-door paths of `v`.  This
    /// is the graph in which the back-door and front-door criteria test
    /// d-separation.  The receiver is left untouched.
    ///
    /// # Panics
    /// Panics if `v >= n_nodes`.
    pub fn without_outgoing(&self, v: usize) -> CausalGraph {
        let mut modified = self.clone();
        let former_children = std::mem::take(&mut modified.children[v]);
        for child in former_children {
            modified.parents[child].retain(|&p| p != v);
        }
        modified
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
        // We check this by removing the outgoing edges of X and then testing
        // d-separation in the modified graph.
        self.graph
            .without_outgoing(treatment)
            .d_separated(&[treatment], &[outcome], adjustment_set)
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

// ---------------------------------------------------------------------------
// FrontdoorCriterion
// ---------------------------------------------------------------------------

/// Implements the frontdoor adjustment formula for causal effect identification.
///
/// The frontdoor criterion allows identification of causal effects through a
/// mediator set `M` when direct adjustment is not possible (for example when
/// an unobserved confounder `U` affects both `X` and `Y`).
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
    /// the effect of `treatment` on `outcome` (Pearl, *Causality*, Def. 3.3.3).
    ///
    /// Conditions:
    /// 1. All directed paths from `treatment` to `outcome` are intercepted by `M`.
    /// 2. There is no unblocked back-door path from `treatment` to `M`.
    /// 3. All back-door paths from `M` to `outcome` are blocked by `treatment`.
    ///
    /// Conditions 2 and 3 are tested per mediator by d-separation in the graph
    /// from which the outgoing edges of the path's start node were removed (see
    /// [`CausalGraph::without_outgoing`]), so only back-door paths remain.
    ///
    /// Returns `false` for out-of-range node indices and for mediator sets that
    /// contain the treatment or the outcome.
    pub fn check(&self, treatment: usize, outcome: usize, mediator_set: &[usize]) -> bool {
        let n = self.graph.n_nodes;
        if treatment >= n || outcome >= n || treatment == outcome {
            return false;
        }
        if mediator_set
            .iter()
            .any(|&m| m >= n || m == treatment || m == outcome)
        {
            return false;
        }
        let med_set: HashSet<usize> = mediator_set.iter().copied().collect();

        // Condition 1: M intercepts all directed paths from X to Y.
        if !self.intercepts_all_paths(treatment, outcome, &med_set) {
            return false;
        }

        // Condition 2: no unblocked back-door path from X to any m ∈ M.
        // Back-door paths from X are the paths of the graph without X's outgoing
        // edges; they must be blocked by the empty set.
        let without_treatment_out = self.graph.without_outgoing(treatment);
        if !mediator_set
            .iter()
            .all(|&m| without_treatment_out.d_separated(&[treatment], &[m], &[]))
        {
            return false;
        }

        // Condition 3: every back-door path from m ∈ M to Y is blocked by X.
        // The graph without m's outgoing edges keeps exactly the paths that
        // begin with an arrow into m, which X must block.
        mediator_set.iter().all(|&m| {
            self.graph
                .without_outgoing(m)
                .d_separated(&[m], &[outcome], &[treatment])
        })
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

    /// Estimate E\[Y | do(X = x_val)\] from observational data with a single
    /// mediator `M`, using the frontdoor formula.
    ///
    /// The frontdoor formula reads
    ///
    /// E\[Y | do(X=x)\] = Σ_m P(M=m | X=x) · Σ_x' E\[Y | M=m, X=x'\] · P(X=x').
    ///
    /// This estimator is semi-parametric:
    ///
    /// * the outcome regression `E[Y | M, X] = β₀ + β_m·M + β_x·X` is fitted by
    ///   ordinary least squares (centred 2×2 normal equations, Cholesky);
    ///   averaging it over the empirical distribution of `X` gives
    ///   `Σ_x' E[Y | m, x'] P(x') = β₀ + β_m·m + β_x·mean(X)`;
    /// * the mediator regression `E[M | X = x_val]` is estimated by local-linear
    ///   kernel regression (Gaussian kernel, normal-reference bandwidth), which
    ///   is exactly unbiased for a linear conditional mean regardless of the
    ///   design density.
    ///
    /// Combining the two gives `mean(Y) + β_m · (Ê[M | X = x_val] − mean(M))`.
    /// Unlike the naive regression of `Y` on `X`, this is consistent under
    /// unobserved `X`–`Y` confounding when `M` satisfies the frontdoor
    /// criterion (see [`FrontdoorCriterion::check`]) and `E[Y | M, X]` is
    /// linear.
    ///
    /// # Errors
    /// Returns [`PhysicsError::InvalidInput`] if the slices differ in length,
    /// contain fewer than four observations or non-finite values, if `x_val` is
    /// not finite, if `X` or `M` has no variation or the two are perfectly
    /// collinear (the outcome regression is singular), or if no observation
    /// carries kernel weight at `x_val` (the positivity assumption `P(x) > 0`
    /// fails).
    pub fn adjusted_effect(
        data_x: &[f64],
        data_m: &[f64],
        data_y: &[f64],
        x_val: f64,
    ) -> Result<f64, PhysicsError> {
        let n = data_x.len();
        if data_m.len() != n || data_y.len() != n {
            return Err(PhysicsError::InvalidInput {
                field: "data".to_string(),
                reason: format!(
                    "data_x, data_m and data_y must have equal lengths, got {}, {} and {}",
                    n,
                    data_m.len(),
                    data_y.len()
                ),
            });
        }
        if n < 4 {
            return Err(PhysicsError::InvalidInput {
                field: "data".to_string(),
                reason: format!("frontdoor estimation needs at least 4 observations, got {n}"),
            });
        }
        if !x_val.is_finite() {
            return Err(PhysicsError::InvalidInput {
                field: "x_val".to_string(),
                reason: "intervention value must be finite".to_string(),
            });
        }
        for (name, column) in [("data_x", data_x), ("data_m", data_m), ("data_y", data_y)] {
            if column.iter().any(|v| !v.is_finite()) {
                return Err(PhysicsError::InvalidInput {
                    field: name.to_string(),
                    reason: "observations must be finite".to_string(),
                });
            }
        }

        let nf = n as f64;
        let mean_x = data_x.iter().sum::<f64>() / nf;
        let mean_m = data_m.iter().sum::<f64>() / nf;
        let mean_y = data_y.iter().sum::<f64>() / nf;

        // Centred second moments of the regression Y ~ 1 + M + X.
        let (mut s_mm, mut s_xx, mut s_mx, mut s_my, mut s_xy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for ((&x, &m), &y) in data_x.iter().zip(data_m).zip(data_y) {
            let (dx, dm, dy) = (x - mean_x, m - mean_m, y - mean_y);
            s_mm += dm * dm;
            s_xx += dx * dx;
            s_mx += dm * dx;
            s_my += dm * dy;
            s_xy += dx * dy;
        }
        if s_xx <= 0.0 || s_mm <= 0.0 {
            return Err(PhysicsError::InvalidInput {
                field: if s_xx <= 0.0 { "data_x" } else { "data_m" }.to_string(),
                reason: "variable has no variation, the frontdoor regression is undefined"
                    .to_string(),
            });
        }
        // 1 − ρ²(M, X) must be bounded away from zero for a usable fit.
        if s_mm * s_xx - s_mx * s_mx <= 1e-12 * s_mm * s_xx {
            return Err(PhysicsError::InvalidInput {
                field: "data_m".to_string(),
                reason: "mediator is (numerically) collinear with the treatment, \
                         the outcome regression is singular"
                    .to_string(),
            });
        }
        let normal_matrix = Matrix2::new(s_mm, s_mx, s_mx, s_xx);
        let cholesky = normal_matrix
            .cholesky()
            .ok_or_else(|| PhysicsError::InvalidInput {
                field: "data_m".to_string(),
                reason: "normal equations of the outcome regression are not positive definite"
                    .to_string(),
            })?;
        let beta = cholesky.solve(&Vector2::new(s_my, s_xy));
        let beta_m = beta[0];

        let bandwidth = normal_reference_bandwidth(data_x);
        let e_m_given_x = local_linear_mean(data_x, data_m, x_val, bandwidth).ok_or_else(|| {
            PhysicsError::InvalidInput {
                field: "x_val".to_string(),
                reason: "no observations near the intervention value (positivity violated)"
                    .to_string(),
            }
        })?;

        Ok(mean_y + beta_m * (e_m_given_x - mean_m))
    }
}

/// Gaussian kernel function for kernel smoothing.
fn gaussian_kernel(u: f64) -> f64 {
    (-0.5 * u * u).exp()
}

/// Normal-reference (Silverman) bandwidth `1.06 · min(σ̂, IQR/1.34) · n^{-1/5}`.
///
/// Falls back to the sample standard deviation when the inter-quartile range
/// is zero (heavily tied data).  The input must contain at least four finite
/// values with non-zero variance.
fn normal_reference_bandwidth(x: &[f64]) -> f64 {
    let n = x.len();
    let nf = n as f64;
    let mean = x.iter().sum::<f64>() / nf;
    let variance = x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (nf - 1.0);
    let std_dev = variance.sqrt();
    let mut sorted = x.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let iqr = sorted[(3 * n) / 4] - sorted[n / 4];
    let scale = if iqr > 0.0 {
        std_dev.min(iqr / 1.34)
    } else {
        std_dev
    };
    1.06 * scale * nf.powf(-0.2)
}

/// Local-linear Gaussian-kernel estimate of `E[response | covariate = x0]`.
///
/// Fits a weighted least-squares line through the observations around `x0` and
/// evaluates it at `x0`.  Unlike the Nadaraya–Watson (local-constant)
/// estimator it has no design bias: it reproduces linear conditional means
/// exactly for any bandwidth.  When the local design is degenerate (all weight
/// on a single covariate value) it falls back to the weighted mean.  Returns
/// `None` when no observation carries kernel weight at `x0`.
fn local_linear_mean(covariate: &[f64], response: &[f64], x0: f64, bandwidth: f64) -> Option<f64> {
    let (mut s0, mut s1, mut s2, mut t0, mut t1) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for (&xi, &ri) in covariate.iter().zip(response) {
        // Work in bandwidth units for a well-conditioned 2×2 system.
        let u = (xi - x0) / bandwidth;
        let w = gaussian_kernel(u);
        s0 += w;
        s1 += w * u;
        s2 += w * u * u;
        t0 += w * ri;
        t1 += w * u * ri;
    }
    if s0.is_nan() || s0 <= 1e-12 {
        return None;
    }
    let det = s0 * s2 - s1 * s1;
    if det > 1e-12 * s0 * s2 {
        Some((s2 * t0 - s1 * t1) / det)
    } else {
        Some(t0 / s0)
    }
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

/// Relative singular-value cut-off below which the first-stage design matrix
/// `[1, Z]` (columns scaled to unit length) counts as rank deficient.
const FIRST_STAGE_RANK_TOLERANCE: f64 = 1e-12;

/// Smallest first-stage `R²` for which the instruments count as relevant.
///
/// Below it the explained variance of the treatment is rounding noise, the
/// 2SLS denominator `Cov(d̂, d)` vanishes and the effect is not identified.
const MIN_FIRST_STAGE_R_SQUARED: f64 = 1e-12;

/// First-stage rss ratio `RSS_u / RSS_r` below which the fit counts as exact.
const EXACT_FIT_RSS_RATIO: f64 = 1e-24;

/// Instrumental variables (IV) estimation and Two-Stage Least Squares (2SLS).
///
/// Used when treatment is endogenous (correlated with the error term).
/// Requires a valid instrument `Z` that:
/// 1. Is correlated with the treatment `D`.
/// 2. Affects the outcome `Y` only through `D` (exclusion restriction).
/// 3. Is independent of unobserved confounders.
///
/// The model is `Y = α + β D + ε` with a single endogenous treatment `D` and
/// `k = n_instruments ≥ 1` instruments.
#[derive(Debug, Clone)]
pub struct InstrumentalVariables {
    /// Number of endogenous variables (only a single treatment is supported).
    pub n_endogenous: usize,
    /// Number of instruments.
    pub n_instruments: usize,
    /// First-stage coefficients `[intercept, π₁, …, π_k]` of the regression of
    /// the treatment on a constant and the instruments (length `k + 1`).
    pub first_stage: Vec<f64>,
    /// Second-stage coefficient (treatment → outcome).
    pub second_stage: f64,
}

/// First-stage regression of the treatment on `[1, Z]`.
struct FirstStageFit {
    /// Coefficients `[intercept, π₁, …, π_k]`.
    coefficients: Vec<f64>,
    /// Fitted treatment `d̂ = W π`.
    fitted: Vec<f64>,
    /// Residual sum of squares of the full model, `RSS_u`.
    rss_unrestricted: f64,
    /// Residual sum of squares of the intercept-only model, `RSS_r = Σ(d − d̄)²`.
    rss_restricted: f64,
}

fn invalid_iv_input(field: &str, reason: String) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: field.to_string(),
        reason,
    }
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

    /// Validate the data and run the first-stage OLS of `d` on `[1, Z]`.
    ///
    /// The design matrix columns are scaled to unit length before the singular
    /// value decomposition so that the rank test does not depend on the units of
    /// the instruments.
    fn first_stage_fit(&self, d: &[f64], z: &[Vec<f64>]) -> Result<FirstStageFit, PhysicsError> {
        let n = d.len();
        let k = self.n_instruments;
        if self.n_endogenous != 1 {
            return Err(invalid_iv_input(
                "n_endogenous",
                format!(
                    "{} endogenous variables requested, but only a single treatment is supported",
                    self.n_endogenous
                ),
            ));
        }
        if k == 0 {
            return Err(invalid_iv_input(
                "n_instruments",
                "at least one instrument is required".to_string(),
            ));
        }
        if z.len() != n {
            return Err(invalid_iv_input(
                "z",
                format!("{} instrument rows for {n} observations", z.len()),
            ));
        }
        if let Some((row, values)) = z.iter().enumerate().find(|(_, values)| values.len() != k) {
            return Err(invalid_iv_input(
                "z",
                format!(
                    "row {row} has {} entries, expected n_instruments = {k}",
                    values.len()
                ),
            ));
        }
        if n < k + 2 {
            return Err(invalid_iv_input(
                "d",
                format!(
                    "{n} observations are too few for {k} instruments (need at least {})",
                    k + 2
                ),
            ));
        }
        if d.iter().any(|v| !v.is_finite()) {
            return Err(invalid_iv_input(
                "d",
                "contains a non-finite value".to_string(),
            ));
        }
        if z.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid_iv_input(
                "z",
                "contains a non-finite value".to_string(),
            ));
        }

        // Design matrix W = [1, Z] with unit-length columns.
        let columns = k + 1;
        let mut design = DMatrix::<f64>::zeros(n, columns);
        for (i, row) in z.iter().enumerate() {
            design[(i, 0)] = 1.0;
            for (j, &value) in row.iter().enumerate() {
                design[(i, j + 1)] = value;
            }
        }
        let column_norms: Vec<f64> = (0..columns).map(|j| design.column(j).norm()).collect();
        if column_norms.iter().any(|norm| !norm.is_finite()) {
            return Err(invalid_iv_input(
                "z",
                "values are too large in magnitude to be squared".to_string(),
            ));
        }
        if let Some(j) = column_norms.iter().position(|&norm| norm == 0.0) {
            return Err(invalid_iv_input(
                "z",
                format!("instrument {} is identically zero", j.saturating_sub(1)),
            ));
        }
        for (j, &norm) in column_norms.iter().enumerate() {
            design.column_mut(j).unscale_mut(norm);
        }
        let svd = SVD::try_new_unordered(design, true, true, f64::EPSILON, 1000 * columns)
            .ok_or_else(|| PhysicsError::NumericalDivergence {
                message: "singular value decomposition of the first-stage design did not converge"
                    .to_string(),
            })?;
        let largest = svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);
        let smallest = svd
            .singular_values
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let full_rank = smallest > FIRST_STAGE_RANK_TOLERANCE * largest;
        if !full_rank {
            return Err(invalid_iv_input(
                "z",
                "the instruments are collinear with each other or with the constant, so the \
                 first stage is not identified"
                    .to_string(),
            ));
        }
        let rhs = DVector::from_column_slice(d);
        let scaled_coefficients =
            svd.solve(&rhs, 0.0)
                .map_err(|message| PhysicsError::NumericalDivergence {
                    message: format!("first-stage least squares failed: {message}"),
                })?;
        let coefficients: Vec<f64> = scaled_coefficients
            .iter()
            .zip(&column_norms)
            .map(|(c, norm)| c / norm)
            .collect();
        if coefficients.iter().any(|c| !c.is_finite()) {
            return Err(PhysicsError::NumericalDivergence {
                message: "first-stage coefficients are not finite".to_string(),
            });
        }

        let fitted: Vec<f64> = z
            .iter()
            .map(|row| {
                coefficients[0]
                    + row
                        .iter()
                        .zip(&coefficients[1..])
                        .map(|(value, pi)| value * pi)
                        .sum::<f64>()
            })
            .collect();
        let mean_d = d.iter().sum::<f64>() / n as f64;
        let rss_unrestricted: f64 = d
            .iter()
            .zip(&fitted)
            .map(|(di, fi)| (di - fi).powi(2))
            .sum();
        let rss_restricted: f64 = d.iter().map(|di| (di - mean_d).powi(2)).sum();
        Ok(FirstStageFit {
            coefficients,
            fitted,
            rss_unrestricted,
            rss_restricted,
        })
    }

    /// Fit the 2SLS estimator.
    ///
    /// Stage 1: ordinary least squares of the treatment `d` on a constant and
    /// *all* `k` instruments, `d = π₀ + π₁ z₁ + … + π_k z_k + v`, giving the
    /// fitted treatment `d̂`; the `k + 1` coefficients are stored in
    /// [`first_stage`](Self::first_stage).
    /// Stage 2: the effect of the treatment is
    /// `β = Σ(d̂ − mean(d̂))(y − ȳ) / Σ(d̂ − mean(d̂))(d − d̄)`,
    /// stored in [`second_stage`](Self::second_stage).  With one instrument this
    /// is the Wald estimator `Cov(y, z) / Cov(d, z)`.
    ///
    /// # Arguments
    /// * `y` — outcome variable
    /// * `d` — endogenous treatment variable
    /// * `z` — instruments matrix, shape `[n_obs][n_instruments]`
    ///
    /// # Errors
    ///
    /// [`PhysicsError::InvalidInput`] if `n_endogenous != 1`, there is no
    /// instrument, the slice lengths or the row lengths of `z` disagree, fewer
    /// than `n_instruments + 2` observations are given, any value is non-finite,
    /// the treatment is constant, an instrument is identically zero, the
    /// instruments are collinear (a constant instrument is collinear with the
    /// intercept), or the instruments explain none of the variation of the
    /// treatment (the effect is then not identified).
    /// [`PhysicsError::NumericalDivergence`] if the least-squares solve fails.
    /// On error the estimator is left unchanged.
    pub fn fit_2sls(&mut self, y: &[f64], d: &[f64], z: &[Vec<f64>]) -> Result<(), PhysicsError> {
        if y.len() != d.len() {
            return Err(invalid_iv_input(
                "y",
                format!("{} outcomes for {} treatment values", y.len(), d.len()),
            ));
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(invalid_iv_input(
                "y",
                "contains a non-finite value".to_string(),
            ));
        }
        let stage = self.first_stage_fit(d, z)?;
        let n = d.len() as f64;
        if stage.rss_restricted <= 0.0 {
            return Err(invalid_iv_input(
                "d",
                "the treatment is constant, so its effect is not identified".to_string(),
            ));
        }
        let mean_hat = stage.fitted.iter().sum::<f64>() / n;
        let mean_d = d.iter().sum::<f64>() / n;
        let mean_y = y.iter().sum::<f64>() / n;
        let cov_hat_y: f64 = stage
            .fitted
            .iter()
            .zip(y)
            .map(|(hat, yi)| (hat - mean_hat) * (yi - mean_y))
            .sum();
        let cov_hat_d: f64 = stage
            .fitted
            .iter()
            .zip(d)
            .map(|(hat, di)| (hat - mean_hat) * (di - mean_d))
            .sum();
        // Cov(d̂, d) = R² Σ(d − d̄)²: the instruments must move the treatment.
        let relevant = cov_hat_d > MIN_FIRST_STAGE_R_SQUARED * stage.rss_restricted;
        if !relevant {
            return Err(invalid_iv_input(
                "z",
                "the instruments explain none of the variation of the treatment, so the effect \
                 is not identified"
                    .to_string(),
            ));
        }
        self.second_stage = cov_hat_y / cov_hat_d;
        self.first_stage = stage.coefficients;
        Ok(())
    }

    /// Compute the first-stage F-statistic (instrument relevance test).
    ///
    /// This is the joint F-test that all `k = n_instruments` instrument
    /// coefficients of the regression of `d` on `[1, Z]` vanish,
    /// `F = ((RSS_r − RSS_u) / k) / (RSS_u / (n − k − 1))`, where `RSS_r` is the
    /// residual sum of squares of the intercept-only model and `RSS_u` that of
    /// the full model.  It is computed from the data passed in, not from stored
    /// coefficients.  Large F (> 10) indicates strong instruments; an exact
    /// first-stage fit gives `f64::INFINITY`.
    ///
    /// `y` plays no role in the first stage; it is only checked to have as many
    /// entries as `d`, the sample passed to [`fit_2sls`](Self::fit_2sls).
    ///
    /// # Errors
    ///
    /// The [`PhysicsError::InvalidInput`] and [`PhysicsError::NumericalDivergence`]
    /// cases of [`fit_2sls`](Self::fit_2sls) that concern the first stage, and a
    /// constant treatment.
    pub fn first_stage_f_stat(
        &self,
        y: &[f64],
        d: &[f64],
        z: &[Vec<f64>],
    ) -> Result<f64, PhysicsError> {
        if y.len() != d.len() {
            return Err(invalid_iv_input(
                "y",
                format!("{} outcomes for {} treatment values", y.len(), d.len()),
            ));
        }
        let stage = self.first_stage_fit(d, z)?;
        if stage.rss_restricted <= 0.0 {
            return Err(invalid_iv_input(
                "d",
                "the treatment is constant, so the first-stage F-statistic is undefined"
                    .to_string(),
            ));
        }
        if stage.rss_unrestricted <= EXACT_FIT_RSS_RATIO * stage.rss_restricted {
            return Ok(f64::INFINITY);
        }
        let k = self.n_instruments as f64;
        let degrees_of_freedom = (d.len() - self.n_instruments - 1) as f64;
        let explained = (stage.rss_restricted - stage.rss_unrestricted).max(0.0);
        Ok((explained / k) / (stage.rss_unrestricted / degrees_of_freedom))
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
#[path = "causal_inference_tests.rs"]
mod tests;
