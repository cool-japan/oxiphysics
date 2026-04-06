#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Formal verification methods for physics simulation correctness.
//!
//! This module provides model checking, abstract interpretation, bisimulation,
//! propositional SAT solving, and conservation-law verifiers for physical
//! simulation states.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// IntervalDomain
// ---------------------------------------------------------------------------

/// Closed interval `[lo, hi]` used as the abstract domain for variables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntervalDomain {
    /// Lower bound of the interval.
    pub lo: f64,
    /// Upper bound of the interval.
    pub hi: f64,
}

impl IntervalDomain {
    /// Create a new interval `[lo, hi]`.
    ///
    /// Panics in debug mode if `lo > hi`.
    pub fn new(lo: f64, hi: f64) -> Self {
        debug_assert!(lo <= hi, "invalid interval: {} > {}", lo, hi);
        Self { lo, hi }
    }

    /// The "top" element — represents the universe of all reals.
    pub fn top() -> Self {
        Self {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }

    /// The "bottom" element — represents the empty set (empty interval).
    pub fn bottom() -> Self {
        Self {
            lo: f64::INFINITY,
            hi: f64::NEG_INFINITY,
        }
    }

    /// Returns `true` if the interval is non-empty.
    pub fn is_non_empty(&self) -> bool {
        self.lo <= self.hi
    }

    /// Returns `true` if `v` lies within this interval.
    pub fn contains(&self, v: f64) -> bool {
        v >= self.lo && v <= self.hi
    }

    /// Least upper bound (join): the smallest interval containing both.
    pub fn join(&self, other: &Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// Greatest lower bound (meet): the intersection.
    pub fn meet(&self, other: &Self) -> Self {
        Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        }
    }

    /// Returns `true` if `self ⊆ other`.
    pub fn is_subset_of(&self, other: &Self) -> bool {
        self.lo >= other.lo && self.hi <= other.hi
    }

    /// Interval addition.
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.lo + other.lo, self.hi + other.hi)
    }

    /// Interval subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        Self::new(self.lo - other.hi, self.hi - other.lo)
    }

    /// Interval multiplication.
    pub fn mul(&self, other: &Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        Self::new(
            products.iter().cloned().fold(f64::INFINITY, f64::min),
            products.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        )
    }

    /// Interval division (returns `top()` when divisor contains zero).
    pub fn div(&self, other: &Self) -> Self {
        if other.lo <= 0.0 && other.hi >= 0.0 {
            Self::top()
        } else {
            let quotients = [
                self.lo / other.lo,
                self.lo / other.hi,
                self.hi / other.lo,
                self.hi / other.hi,
            ];
            Self::new(
                quotients.iter().cloned().fold(f64::INFINITY, f64::min),
                quotients.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        }
    }

    /// Width of the interval.
    pub fn width(&self) -> f64 {
        (self.hi - self.lo).max(0.0)
    }

    /// Midpoint of the interval.
    pub fn midpoint(&self) -> f64 {
        (self.lo + self.hi) / 2.0
    }
}

impl std::fmt::Display for IntervalDomain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}]", self.lo, self.hi)
    }
}

// ---------------------------------------------------------------------------
// LinearTemporalLogic
// ---------------------------------------------------------------------------

/// A Linear Temporal Logic (LTL) formula represented as a string.
///
/// The current implementation uses a simplified embedded DSL rather than a
/// full LTL parser.
#[derive(Debug, Clone)]
pub struct LinearTemporalLogic {
    /// The LTL formula expressed in textual form.
    pub formula: String,
}

impl LinearTemporalLogic {
    /// Create a new LTL formula.
    pub fn new(formula: impl Into<String>) -> Self {
        Self {
            formula: formula.into(),
        }
    }

    /// Check a safety property: returns `true` if the formula is satisfied
    /// for all positions in `trace`.
    ///
    /// Safety properties have the form "always P" – here `predicate` is
    /// evaluated on every state.
    pub fn check_safety<S, F>(&self, trace: &[S], predicate: F) -> bool
    where
        F: Fn(&S) -> bool,
    {
        trace.iter().all(predicate)
    }

    /// Check a liveness property: returns `true` if `predicate` is `true` at
    /// least once in `trace`.
    ///
    /// Liveness properties have the form "eventually P".
    pub fn check_liveness<S, F>(&self, trace: &[S], predicate: F) -> bool
    where
        F: Fn(&S) -> bool,
    {
        trace.iter().any(predicate)
    }

    /// Check a "next" property: returns `true` if `predicate` holds at every
    /// state immediately following a state where `trigger` holds.
    pub fn check_next<S, F, G>(&self, trace: &[S], trigger: F, predicate: G) -> bool
    where
        F: Fn(&S) -> bool,
        G: Fn(&S) -> bool,
    {
        for i in 0..trace.len().saturating_sub(1) {
            if trigger(&trace[i]) && !predicate(&trace[i + 1]) {
                return false;
            }
        }
        true
    }

    /// Check an "until" property: `left` must hold until `right` holds.
    pub fn check_until<S, F, G>(&self, trace: &[S], left: F, right: G) -> bool
    where
        F: Fn(&S) -> bool,
        G: Fn(&S) -> bool,
    {
        let mut holding = false;
        for s in trace {
            if right(s) {
                holding = true;
                break;
            }
            if !left(s) {
                return false;
            }
        }
        holding
    }
}

// ---------------------------------------------------------------------------
// ModelChecker
// ---------------------------------------------------------------------------

/// State index type alias.
pub type StateIdx = usize;

/// A simple explicit-state model checker.
///
/// States are represented as `Vec`f64` (e.g., phase-space coordinates).
/// Transitions are labeled pairs `(from, to)`.
#[derive(Debug, Clone)]
pub struct ModelChecker {
    /// All states in the model.
    pub states: Vec<Vec<f64>>,
    /// Transitions as `(from_idx, to_idx)` pairs.
    pub transitions: Vec<(StateIdx, StateIdx)>,
}

impl ModelChecker {
    /// Create a new model checker from a set of states and transitions.
    pub fn new(states: Vec<Vec<f64>>, transitions: Vec<(StateIdx, StateIdx)>) -> Self {
        Self {
            states,
            transitions,
        }
    }

    /// Return all states reachable from `start_idx` via BFS.
    pub fn reachable_states(&self, start_idx: StateIdx) -> HashSet<StateIdx> {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(start_idx);
        visited.insert(start_idx);
        while let Some(s) = queue.pop_front() {
            for &(from, to) in &self.transitions {
                if from == s && !visited.contains(&to) {
                    visited.insert(to);
                    queue.push_back(to);
                }
            }
        }
        visited
    }

    /// Check whether all reachable states from `start_idx` satisfy `invariant`.
    pub fn satisfies_invariant<F>(&self, start_idx: StateIdx, invariant: F) -> bool
    where
        F: Fn(&[f64]) -> bool,
    {
        let reachable = self.reachable_states(start_idx);
        reachable.iter().all(|&s| invariant(&self.states[s]))
    }

    /// Find a counterexample state (as index) that violates `invariant`.
    ///
    /// Returns `None` if all reachable states satisfy the invariant.
    pub fn find_counterexample<F>(&self, start_idx: StateIdx, invariant: F) -> Option<StateIdx>
    where
        F: Fn(&[f64]) -> bool,
    {
        let reachable = self.reachable_states(start_idx);
        reachable.into_iter().find(|&s| !invariant(&self.states[s]))
    }

    /// Return the shortest path (as a sequence of state indices) from
    /// `start_idx` to `target_idx`, or `None` if unreachable.
    pub fn shortest_path(
        &self,
        start_idx: StateIdx,
        target_idx: StateIdx,
    ) -> Option<Vec<StateIdx>> {
        let mut parent: HashMap<StateIdx, StateIdx> = HashMap::new();
        let mut queue = VecDeque::new();
        queue.push_back(start_idx);
        parent.insert(start_idx, start_idx);
        while let Some(s) = queue.pop_front() {
            if s == target_idx {
                // Reconstruct path.
                let mut path = Vec::new();
                let mut cur = target_idx;
                loop {
                    path.push(cur);
                    let prev = parent[&cur];
                    if prev == cur {
                        break;
                    }
                    cur = prev;
                }
                path.reverse();
                return Some(path);
            }
            for &(from, to) in &self.transitions {
                if from == s && !parent.contains_key(&to) {
                    parent.insert(to, s);
                    queue.push_back(to);
                }
            }
        }
        None
    }

    /// Count the number of reachable states from `start_idx`.
    pub fn reachable_count(&self, start_idx: StateIdx) -> usize {
        self.reachable_states(start_idx).len()
    }
}

// ---------------------------------------------------------------------------
// AbstractInterpretation
// ---------------------------------------------------------------------------

/// Abstract interpretation engine using the interval domain.
#[derive(Debug, Clone)]
pub struct AbstractInterpretation {
    /// Abstract domain for each tracked variable.
    pub domain: HashMap<String, IntervalDomain>,
}

impl AbstractInterpretation {
    /// Create a new abstract interpreter with no variable bindings.
    pub fn new() -> Self {
        Self {
            domain: HashMap::new(),
        }
    }

    /// Bind a variable to a concrete value (singleton interval).
    pub fn bind(&mut self, var: impl Into<String>, value: f64) {
        let name = var.into();
        self.domain.insert(name, IntervalDomain::new(value, value));
    }

    /// Bind a variable to an abstract interval.
    pub fn bind_interval(&mut self, var: impl Into<String>, interval: IntervalDomain) {
        self.domain.insert(var.into(), interval);
    }

    /// Return the abstract value of a variable, or `top()` if unknown.
    pub fn get(&self, var: &str) -> IntervalDomain {
        self.domain
            .get(var)
            .copied()
            .unwrap_or_else(IntervalDomain::top)
    }

    /// Analyze a simple loop that adds `delta` to `var` each iteration for at
    /// most `max_iter` iterations, using widening to ensure termination.
    ///
    /// Returns the post-loop interval for `var`.
    pub fn analyze_loop(
        &mut self,
        var: &str,
        delta: IntervalDomain,
        max_iter: usize,
    ) -> IntervalDomain {
        let mut current = self.get(var);
        for _i in 0..max_iter {
            let next = current.add(&delta);
            let joined = current.join(&next);
            let widened = self.widening(&current, &joined);
            if widened == current {
                return current;
            }
            current = widened;
        }
        current
    }

    /// Widening operator: extends the interval to ensure convergence.
    ///
    /// If the lower bound decreases it goes to `-∞`; if the upper bound
    /// increases it goes to `+∞`.
    pub fn widening(&self, prev: &IntervalDomain, next: &IntervalDomain) -> IntervalDomain {
        let lo = if next.lo < prev.lo {
            f64::NEG_INFINITY
        } else {
            prev.lo
        };
        let hi = if next.hi > prev.hi {
            f64::INFINITY
        } else {
            prev.hi
        };
        IntervalDomain { lo, hi }
    }

    /// Narrowing operator: refines an over-approximation using additional
    /// concrete information.
    pub fn narrowing(&self, prev: &IntervalDomain, next: &IntervalDomain) -> IntervalDomain {
        let lo = if prev.lo == f64::NEG_INFINITY {
            next.lo
        } else {
            prev.lo
        };
        let hi = if prev.hi == f64::INFINITY {
            next.hi
        } else {
            prev.hi
        };
        IntervalDomain { lo, hi }
    }

    /// Check that `var` is always non-negative under the current abstraction.
    pub fn is_non_negative(&self, var: &str) -> bool {
        self.get(var).lo >= 0.0
    }

    /// Check that `var` is bounded (finite interval bounds).
    pub fn is_bounded(&self, var: &str) -> bool {
        let iv = self.get(var);
        iv.lo.is_finite() && iv.hi.is_finite()
    }
}

impl Default for AbstractInterpretation {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BisimulationChecker
// ---------------------------------------------------------------------------

/// Labeled transition system used for bisimulation checking.
///
/// Labels on transitions are `String`s.
#[derive(Debug, Clone)]
pub struct LabeledTransitionSystem {
    /// Number of states.
    pub n_states: usize,
    /// Transitions: `(from, label, to)`.
    pub transitions: Vec<(usize, String, usize)>,
    /// Optional state labels (atomic propositions).
    pub labels: Vec<Vec<String>>,
}

impl LabeledTransitionSystem {
    /// Create a new LTS with `n_states` states and no transitions.
    pub fn new(n_states: usize) -> Self {
        Self {
            n_states,
            transitions: Vec::new(),
            labels: vec![Vec::new(); n_states],
        }
    }

    /// Add a labeled transition `from --label--> to`.
    pub fn add_transition(&mut self, from: usize, label: impl Into<String>, to: usize) {
        self.transitions.push((from, label.into(), to));
    }

    /// Add an atomic proposition to a state.
    pub fn add_label(&mut self, state: usize, label: impl Into<String>) {
        if state < self.labels.len() {
            self.labels[state].push(label.into());
        }
    }
}

/// Bisimulation checker for labeled transition systems.
#[derive(Debug)]
pub struct BisimulationChecker {
    /// The LTS to analyze.
    pub lts: LabeledTransitionSystem,
}

impl BisimulationChecker {
    /// Create a new bisimulation checker for `lts`.
    pub fn new(lts: LabeledTransitionSystem) -> Self {
        Self { lts }
    }

    /// Check whether states `s` and `t` are bisimilar in the LTS.
    ///
    /// Uses a partition-refinement approximation (up to `max_rounds` rounds).
    pub fn are_bisimilar(&self, s: usize, t: usize) -> bool {
        let partition = self.compute_quotient();
        let s_class = partition.get(&s).copied().unwrap_or(s);
        let t_class = partition.get(&t).copied().unwrap_or(t);
        s_class == t_class
    }

    /// Compute the bisimulation quotient using Paige-Tarjan partition
    /// refinement (simplified version).
    ///
    /// Returns a map from state index to equivalence class index.
    pub fn compute_quotient(&self) -> HashMap<usize, usize> {
        let n = self.lts.n_states;
        // Initial partition: group by atomic propositions.
        let mut partition: Vec<usize> = (0..n).map(|_| 0).collect();
        // Assign initial classes based on label sets.
        let mut label_to_class: HashMap<Vec<String>, usize> = HashMap::new();
        let mut next_class = 0usize;
        for (i, labels) in self.lts.labels.iter().enumerate() {
            let mut sorted_labels = labels.clone();
            sorted_labels.sort();
            let class = *label_to_class.entry(sorted_labels).or_insert_with(|| {
                let c = next_class;
                next_class += 1;
                c
            });
            partition[i] = class;
        }

        // Refinement loop.
        loop {
            let old_partition = partition.clone();
            // For each state, compute a signature: (class, sorted list of (label, target_class)).
            let mut signatures: Vec<(usize, Vec<(String, usize)>)> = Vec::new();
            for i in 0..n {
                let mut sig: Vec<(String, usize)> = self
                    .lts
                    .transitions
                    .iter()
                    .filter(|(from, _, _)| *from == i)
                    .map(|(_, label, to)| (label.clone(), partition[*to]))
                    .collect();
                sig.sort();
                signatures.push((partition[i], sig));
            }
            // Re-assign classes.
            let mut sig_to_class: HashMap<(usize, Vec<(String, usize)>), usize> = HashMap::new();
            next_class = 0;
            for (i, sig) in signatures.iter().enumerate() {
                let class = *sig_to_class.entry(sig.clone()).or_insert_with(|| {
                    let c = next_class;
                    next_class += 1;
                    c
                });
                partition[i] = class;
            }
            if partition == old_partition {
                break;
            }
        }
        partition.into_iter().enumerate().collect()
    }

    /// Return the number of equivalence classes in the quotient.
    pub fn quotient_size(&self) -> usize {
        let q = self.compute_quotient();
        let classes: HashSet<usize> = q.values().copied().collect();
        classes.len()
    }
}

// ---------------------------------------------------------------------------
// Conservation law verifiers
// ---------------------------------------------------------------------------

/// Verify energy conservation over a trajectory of states.
///
/// Each state is a `\[f64; 6\]` array of the form `\[x, y, z, vx, vy, vz\]`.
/// Energy is estimated as kinetic energy `0.5 * (vx² + vy² + vz²)` (unit mass).
///
/// Returns `true` if the energy variation across all states is within `tol`.
pub fn verify_energy_conservation(states: &[[f64; 6]], tol: f64) -> bool {
    if states.len() < 2 {
        return true;
    }
    let energies: Vec<f64> = states
        .iter()
        .map(|s| 0.5 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
        .collect();
    let e0 = energies[0];
    energies.iter().all(|&e| (e - e0).abs() <= tol)
}

/// Verify linear momentum conservation over a trajectory of states.
///
/// Each state is `\[x, y, z, vx, vy, vz\]` (unit mass).
/// Returns `true` if the momentum vector changes by less than `tol` per
/// component.
pub fn verify_momentum_conservation(states: &[[f64; 6]], tol: f64) -> bool {
    if states.len() < 2 {
        return true;
    }
    let (px0, py0, pz0) = (states[0][3], states[0][4], states[0][5]);
    states.iter().all(|s| {
        (s[3] - px0).abs() <= tol && (s[4] - py0).abs() <= tol && (s[5] - pz0).abs() <= tol
    })
}

/// Verify angular momentum conservation.
///
/// Each state is `\[x, y, z, vx, vy, vz\]`.  Returns `true` if the z-component
/// of `r × v` changes by less than `tol`.
pub fn verify_angular_momentum_conservation(states: &[[f64; 6]], tol: f64) -> bool {
    if states.len() < 2 {
        return true;
    }
    // Lz = x*vy - y*vx
    let lz0 = states[0][0] * states[0][4] - states[0][1] * states[0][3];
    states
        .iter()
        .all(|s| ((s[0] * s[4] - s[1] * s[3]) - lz0).abs() <= tol)
}

// ---------------------------------------------------------------------------
// SatisfiabilityChecker (DPLL)
// ---------------------------------------------------------------------------

/// A literal in a propositional clause.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Literal {
    /// Variable index (0-based).
    pub var: usize,
    /// Polarity: `true` = positive, `false` = negative.
    pub positive: bool,
}

impl Literal {
    /// Positive literal for variable `var`.
    pub fn pos(var: usize) -> Self {
        Self {
            var,
            positive: true,
        }
    }

    /// Negative literal for variable `var`.
    pub fn neg(var: usize) -> Self {
        Self {
            var,
            positive: false,
        }
    }

    /// Return the negation of this literal.
    pub fn negate(&self) -> Self {
        Self {
            var: self.var,
            positive: !self.positive,
        }
    }
}

/// A propositional formula in Conjunctive Normal Form (CNF).
///
/// Clauses are disjunctions; the formula is the conjunction of all clauses.
#[derive(Debug, Clone)]
pub struct SatisfiabilityChecker {
    /// Number of propositional variables.
    pub n_vars: usize,
    /// CNF clauses.
    pub clauses: Vec<Vec<Literal>>,
}

impl SatisfiabilityChecker {
    /// Create an empty SAT checker with `n_vars` variables.
    pub fn new(n_vars: usize) -> Self {
        Self {
            n_vars,
            clauses: Vec::new(),
        }
    }

    /// Add a clause (disjunction of literals).
    pub fn add_clause(&mut self, clause: Vec<Literal>) {
        self.clauses.push(clause);
    }

    /// Check satisfiability using the DPLL algorithm.
    ///
    /// Returns `Some(assignment)` where `assignment\[i\]` is the truth value
    /// for variable `i`, or `None` if unsatisfiable.
    pub fn solve(&self) -> Option<Vec<bool>> {
        let mut assignment = vec![None::<bool>; self.n_vars];
        if self.dpll(&mut assignment) {
            Some(assignment.iter().map(|a| a.unwrap_or(true)).collect())
        } else {
            None
        }
    }

    fn dpll(&self, assignment: &mut Vec<Option<bool>>) -> bool {
        // Unit propagation.
        loop {
            let mut propagated = false;
            for clause in &self.clauses {
                let (unassigned, satisfied) = self.clause_status(clause, assignment);
                if satisfied {
                    continue;
                }
                if unassigned.is_empty() {
                    return false; // Conflict.
                }
                if unassigned.len() == 1 {
                    // Unit clause: force the literal.
                    let lit = unassigned[0];
                    assignment[lit.var] = Some(lit.positive);
                    propagated = true;
                }
            }
            if !propagated {
                break;
            }
        }

        // Check if all clauses are satisfied.
        if self.all_satisfied(assignment) {
            return true;
        }
        // Check for conflict.
        if self.has_conflict(assignment) {
            return false;
        }

        // Choose an unassigned variable and branch.
        if let Some(var) = (0..self.n_vars).find(|&v| assignment[v].is_none()) {
            for &value in &[true, false] {
                let mut child = assignment.clone();
                child[var] = Some(value);
                if self.dpll(&mut child) {
                    *assignment = child;
                    return true;
                }
            }
        }
        false
    }

    fn clause_status<'a>(
        &self,
        clause: &'a [Literal],
        assignment: &[Option<bool>],
    ) -> (Vec<&'a Literal>, bool) {
        let mut unassigned = Vec::new();
        for lit in clause {
            match assignment[lit.var] {
                Some(val) if val == lit.positive => return (vec![], true),
                None => unassigned.push(lit),
                _ => {}
            }
        }
        (unassigned, false)
    }

    fn all_satisfied(&self, assignment: &[Option<bool>]) -> bool {
        self.clauses.iter().all(|clause| {
            clause
                .iter()
                .any(|lit| assignment[lit.var] == Some(lit.positive))
        })
    }

    fn has_conflict(&self, assignment: &[Option<bool>]) -> bool {
        self.clauses.iter().any(|clause| {
            clause
                .iter()
                .all(|lit| assignment[lit.var].is_some_and(|v| v != lit.positive))
        })
    }

    /// Check whether the given assignment satisfies all clauses.
    pub fn check_assignment(&self, assignment: &[bool]) -> bool {
        self.clauses
            .iter()
            .all(|clause| clause.iter().any(|lit| assignment[lit.var] == lit.positive))
    }

    /// Return the number of clauses.
    pub fn clause_count(&self) -> usize {
        self.clauses.len()
    }
}

// ---------------------------------------------------------------------------
// Lyapunov stability checker
// ---------------------------------------------------------------------------

/// Check Lyapunov stability of a fixed-point under a discrete map.
///
/// `states` is a trajectory; the function checks whether `‖state‖₂` is
/// non-increasing (on average) — a sufficient condition for Lyapunov stability.
pub fn check_lyapunov_stability(states: &[Vec<f64>]) -> bool {
    if states.len() < 2 {
        return true;
    }
    let norm = |s: &Vec<f64>| s.iter().map(|x| x * x).sum::<f64>().sqrt();
    let first_norm = norm(&states[0]);
    let last_norm = norm(states.last().expect("states has at least 2 entries"));
    last_norm <= first_norm * 1.001 // 0.1% tolerance for numerical drift
}

// ---------------------------------------------------------------------------
// Invariant specification helpers
// ---------------------------------------------------------------------------

/// A named invariant with a predicate over phase-space vectors.
pub struct Invariant {
    /// Name of the invariant.
    pub name: String,
    /// Predicate: returns `true` if the state satisfies the invariant.
    pub predicate: Box<dyn Fn(&[f64]) -> bool + Send + Sync>,
}

impl Invariant {
    /// Create a new invariant.
    pub fn new(
        name: impl Into<String>,
        predicate: impl Fn(&[f64]) -> bool + Send + Sync + 'static,
    ) -> Self {
        Self {
            name: name.into(),
            predicate: Box::new(predicate),
        }
    }

    /// Evaluate the invariant on a state.
    pub fn check(&self, state: &[f64]) -> bool {
        (self.predicate)(state)
    }

    /// Verify the invariant holds for an entire trajectory.
    pub fn verify_trajectory(&self, trajectory: &[Vec<f64>]) -> bool {
        trajectory.iter().all(|s| self.check(s))
    }
}

impl std::fmt::Debug for Invariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Invariant")
            .field("name", &self.name)
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Symbolic execution engine (simple numeric version)
// ---------------------------------------------------------------------------

/// A symbolic state: each variable holds an interval.
pub type SymbolicState = HashMap<String, IntervalDomain>;

/// Symbolic execution step: apply a numeric update `x := x + delta` to the
/// symbolic state, and check whether the result remains within `bounds`.
pub fn symbolic_step(
    state: &SymbolicState,
    var: &str,
    delta: f64,
    bounds: &IntervalDomain,
) -> Option<SymbolicState> {
    let iv = state.get(var).copied().unwrap_or_else(IntervalDomain::top);
    let next_iv = IntervalDomain::new(iv.lo + delta, iv.hi + delta);
    if next_iv.meet(bounds).is_non_empty() {
        let mut next = state.clone();
        next.insert(var.to_string(), next_iv.meet(bounds));
        Some(next)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Counterexample-guided abstraction refinement (CEGAR) stub
// ---------------------------------------------------------------------------

/// A simple CEGAR loop result.
#[derive(Debug, Clone)]
pub enum CegarResult {
    /// The property was verified (no counterexample exists).
    Verified,
    /// A concrete counterexample was found (sequence of state indices).
    Counterexample(Vec<StateIdx>),
    /// The analysis exceeded the iteration budget without a definitive answer.
    Unknown,
}

/// Run a simplified CEGAR loop on `checker`.
///
/// Checks `invariant` starting from `start_idx`; returns the result.
pub fn cegar_verify<F>(
    checker: &ModelChecker,
    start_idx: StateIdx,
    invariant: F,
    _max_iterations: usize,
) -> CegarResult
where
    F: Fn(&[f64]) -> bool,
{
    match checker.find_counterexample(start_idx, &invariant) {
        None => CegarResult::Verified,
        Some(cex_state) => {
            // Build a simple path to the counterexample.
            let path = checker
                .shortest_path(start_idx, cex_state)
                .unwrap_or_else(|| vec![cex_state]);
            CegarResult::Counterexample(path)
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation trace verifier
// ---------------------------------------------------------------------------

/// Verify a collection of named properties against a simulation trace.
///
/// Returns a map from property name to `true`/`false`.
pub fn verify_trace_properties(
    trace: &[Vec<f64>],
    properties: &[(&str, Box<dyn Fn(&[Vec<f64>]) -> bool>)],
) -> HashMap<String, bool> {
    properties
        .iter()
        .map(|(name, predicate)| (name.to_string(), predicate(trace)))
        .collect()
}

// ---------------------------------------------------------------------------
// Type-state protocol checker
// ---------------------------------------------------------------------------

/// Protocol state machine for type-state verification.
#[derive(Debug, Clone)]
pub struct TypeStateProtocol {
    /// Current state of the protocol.
    pub current_state: String,
    /// Allowed transitions: `(from, event, to)`.
    pub transitions: Vec<(String, String, String)>,
}

impl TypeStateProtocol {
    /// Create a new protocol starting in `initial_state`.
    pub fn new(initial_state: impl Into<String>) -> Self {
        Self {
            current_state: initial_state.into(),
            transitions: Vec::new(),
        }
    }

    /// Register an allowed transition.
    pub fn add_transition(
        &mut self,
        from: impl Into<String>,
        event: impl Into<String>,
        to: impl Into<String>,
    ) {
        self.transitions
            .push((from.into(), event.into(), to.into()));
    }

    /// Attempt to fire `event`.  Returns `Ok(new_state)` or `Err(msg)`.
    pub fn fire(&mut self, event: &str) -> Result<String, String> {
        for (from, ev, to) in &self.transitions {
            if from == &self.current_state && ev == event {
                self.current_state = to.clone();
                return Ok(to.clone());
            }
        }
        Err(format!(
            "no transition from '{}' on event '{}'",
            self.current_state, event
        ))
    }

    /// Return `true` if `state` is reachable from the initial state.
    pub fn is_reachable(&self, initial: &str, target: &str) -> bool {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(initial.to_string());
        visited.insert(initial.to_string());
        while let Some(s) = queue.pop_front() {
            if s == target {
                return true;
            }
            for (from, _, to) in &self.transitions {
                if from == &s && !visited.contains(to) {
                    visited.insert(to.clone());
                    queue.push_back(to.clone());
                }
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Temporal property macros (as free functions)
// ---------------------------------------------------------------------------

/// Return `true` if `predicate` holds for every element of `trace`.
pub fn always<S, F: Fn(&S) -> bool>(trace: &[S], predicate: F) -> bool {
    trace.iter().all(predicate)
}

/// Return `true` if `predicate` holds for at least one element of `trace`.
pub fn eventually<S, F: Fn(&S) -> bool>(trace: &[S], predicate: F) -> bool {
    trace.iter().any(predicate)
}

/// Return `true` if `predicate` holds at every position strictly after
/// a position where `trigger` holds.
pub fn globally_after<S, F: Fn(&S) -> bool, G: Fn(&S) -> bool>(
    trace: &[S],
    trigger: F,
    predicate: G,
) -> bool {
    let mut triggered = false;
    for s in trace {
        if triggered && !predicate(s) {
            return false;
        }
        if trigger(s) {
            triggered = true;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- IntervalDomain --

    #[test]
    fn test_interval_contains() {
        let iv = IntervalDomain::new(1.0, 3.0);
        assert!(iv.contains(2.0));
        assert!(!iv.contains(4.0));
    }

    #[test]
    fn test_interval_join() {
        let a = IntervalDomain::new(1.0, 3.0);
        let b = IntervalDomain::new(2.0, 5.0);
        let j = a.join(&b);
        assert_eq!(j.lo, 1.0);
        assert_eq!(j.hi, 5.0);
    }

    #[test]
    fn test_interval_meet() {
        let a = IntervalDomain::new(1.0, 4.0);
        let b = IntervalDomain::new(2.0, 6.0);
        let m = a.meet(&b);
        assert_eq!(m.lo, 2.0);
        assert_eq!(m.hi, 4.0);
    }

    #[test]
    fn test_interval_add() {
        let a = IntervalDomain::new(1.0, 2.0);
        let b = IntervalDomain::new(3.0, 4.0);
        let c = a.add(&b);
        assert_eq!(c.lo, 4.0);
        assert_eq!(c.hi, 6.0);
    }

    #[test]
    fn test_interval_sub() {
        let a = IntervalDomain::new(3.0, 5.0);
        let b = IntervalDomain::new(1.0, 2.0);
        let c = a.sub(&b);
        assert_eq!(c.lo, 1.0);
        assert_eq!(c.hi, 4.0);
    }

    #[test]
    fn test_interval_mul() {
        let a = IntervalDomain::new(2.0, 3.0);
        let b = IntervalDomain::new(4.0, 5.0);
        let c = a.mul(&b);
        assert_eq!(c.lo, 8.0);
        assert_eq!(c.hi, 15.0);
    }

    #[test]
    fn test_interval_div_no_zero() {
        let a = IntervalDomain::new(6.0, 8.0);
        let b = IntervalDomain::new(2.0, 4.0);
        let c = a.div(&b);
        assert!(c.lo >= 1.5 - 1e-9);
        assert!(c.hi <= 4.0 + 1e-9);
    }

    #[test]
    fn test_interval_div_contains_zero() {
        let a = IntervalDomain::new(1.0, 2.0);
        let b = IntervalDomain::new(-1.0, 1.0);
        let c = a.div(&b);
        assert_eq!(c, IntervalDomain::top());
    }

    #[test]
    fn test_interval_width() {
        let iv = IntervalDomain::new(2.0, 5.0);
        assert!((iv.width() - 3.0).abs() < 1e-10);
    }

    // -- LinearTemporalLogic --

    #[test]
    fn test_ltl_check_safety() {
        let ltl = LinearTemporalLogic::new("always x >= 0");
        let trace = vec![1.0_f64, 2.0, 3.0];
        assert!(ltl.check_safety(&trace, |&x| x >= 0.0));
        let bad_trace = vec![1.0_f64, -1.0, 3.0];
        assert!(!ltl.check_safety(&bad_trace, |&x| x >= 0.0));
    }

    #[test]
    fn test_ltl_check_liveness() {
        let ltl = LinearTemporalLogic::new("eventually x > 10");
        let trace = vec![1.0_f64, 5.0, 11.0];
        assert!(ltl.check_liveness(&trace, |&x| x > 10.0));
        let bad = vec![1.0_f64, 2.0, 3.0];
        assert!(!ltl.check_liveness(&bad, |&x| x > 10.0));
    }

    #[test]
    fn test_ltl_check_next() {
        let ltl = LinearTemporalLogic::new("next");
        let trace = vec![true, true, false];
        // trigger = identity; next must also be true
        assert!(!ltl.check_next(&trace, |&x| x, |&x| x));
    }

    #[test]
    fn test_ltl_check_until() {
        let ltl = LinearTemporalLogic::new("until");
        let trace = vec![1i32, 1, 1, 2];
        assert!(ltl.check_until(&trace, |&x| x == 1, |&x| x == 2));
        let bad = vec![1i32, 0, 2];
        assert!(!ltl.check_until(&bad, |&x| x == 1, |&x| x == 2));
    }

    // -- ModelChecker --

    fn simple_model() -> ModelChecker {
        // States: 0 → 1 → 2
        ModelChecker::new(vec![vec![0.0], vec![1.0], vec![2.0]], vec![(0, 1), (1, 2)])
    }

    #[test]
    fn test_model_reachable_states() {
        let mc = simple_model();
        let r = mc.reachable_states(0);
        assert!(r.contains(&0));
        assert!(r.contains(&1));
        assert!(r.contains(&2));
    }

    #[test]
    fn test_model_satisfies_invariant() {
        let mc = simple_model();
        assert!(mc.satisfies_invariant(0, |s| s[0] >= 0.0));
        assert!(!mc.satisfies_invariant(0, |s| s[0] > 1.5));
    }

    #[test]
    fn test_model_find_counterexample() {
        let mc = simple_model();
        let cex = mc.find_counterexample(0, |s| s[0] < 2.0);
        assert_eq!(cex, Some(2));
    }

    #[test]
    fn test_model_shortest_path() {
        let mc = simple_model();
        let path = mc.shortest_path(0, 2).unwrap();
        assert_eq!(path, vec![0, 1, 2]);
    }

    // -- AbstractInterpretation --

    #[test]
    fn test_abstract_interp_bind_get() {
        let mut ai = AbstractInterpretation::new();
        ai.bind("x", 3.0);
        assert_eq!(ai.get("x"), IntervalDomain::new(3.0, 3.0));
    }

    #[test]
    fn test_abstract_interp_widening() {
        let ai = AbstractInterpretation::new();
        let prev = IntervalDomain::new(0.0, 1.0);
        let next = IntervalDomain::new(0.0, 2.0);
        let w = ai.widening(&prev, &next);
        assert_eq!(w.hi, f64::INFINITY);
    }

    #[test]
    fn test_abstract_interp_narrowing() {
        let ai = AbstractInterpretation::new();
        let prev = IntervalDomain::new(f64::NEG_INFINITY, f64::INFINITY);
        let next = IntervalDomain::new(-5.0, 5.0);
        let n = ai.narrowing(&prev, &next);
        assert_eq!(n.lo, -5.0);
        assert_eq!(n.hi, 5.0);
    }

    #[test]
    fn test_abstract_interp_is_non_negative() {
        let mut ai = AbstractInterpretation::new();
        ai.bind_interval("x", IntervalDomain::new(0.0, 10.0));
        assert!(ai.is_non_negative("x"));
        ai.bind_interval("y", IntervalDomain::new(-1.0, 10.0));
        assert!(!ai.is_non_negative("y"));
    }

    // -- BisimulationChecker --

    #[test]
    fn test_bisimulation_same_state() {
        let mut lts = LabeledTransitionSystem::new(2);
        lts.add_label(0, "a");
        lts.add_label(1, "a");
        lts.add_transition(0, "x", 0);
        lts.add_transition(1, "x", 1);
        let checker = BisimulationChecker::new(lts);
        assert!(checker.are_bisimilar(0, 1));
    }

    #[test]
    fn test_bisimulation_different_labels() {
        let mut lts = LabeledTransitionSystem::new(2);
        lts.add_label(0, "a");
        lts.add_label(1, "b");
        let checker = BisimulationChecker::new(lts);
        assert!(!checker.are_bisimilar(0, 1));
    }

    #[test]
    fn test_bisimulation_quotient_size() {
        let mut lts = LabeledTransitionSystem::new(3);
        lts.add_label(0, "a");
        lts.add_label(1, "a");
        lts.add_label(2, "b");
        let checker = BisimulationChecker::new(lts);
        assert_eq!(checker.quotient_size(), 2);
    }

    // -- Conservation verifiers --

    #[test]
    fn test_verify_energy_conservation_ok() {
        let states = vec![[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]; 5];
        assert!(verify_energy_conservation(&states, 1e-9));
    }

    #[test]
    fn test_verify_energy_conservation_fail() {
        let s1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let s2 = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        assert!(!verify_energy_conservation(&[s1, s2], 1e-9));
    }

    #[test]
    fn test_verify_momentum_conservation_ok() {
        let states = vec![[0.0, 0.0, 0.0, 1.0, 2.0, 3.0]; 4];
        assert!(verify_momentum_conservation(&states, 1e-9));
    }

    #[test]
    fn test_verify_momentum_conservation_fail() {
        let s1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let s2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(!verify_momentum_conservation(&[s1, s2], 1e-9));
    }

    // -- SatisfiabilityChecker --

    #[test]
    fn test_sat_trivial_sat() {
        // (x0) — always satisfiable.
        let mut sat = SatisfiabilityChecker::new(1);
        sat.add_clause(vec![Literal::pos(0)]);
        let result = sat.solve();
        assert!(result.is_some());
        let assignment = result.unwrap();
        assert!(sat.check_assignment(&assignment));
    }

    #[test]
    fn test_sat_trivial_unsat() {
        // (x0) ∧ (¬x0) — unsatisfiable.
        let mut sat = SatisfiabilityChecker::new(1);
        sat.add_clause(vec![Literal::pos(0)]);
        sat.add_clause(vec![Literal::neg(0)]);
        assert!(sat.solve().is_none());
    }

    #[test]
    fn test_sat_two_vars() {
        // (x0 ∨ x1) ∧ (¬x0 ∨ x1) — satisfiable with x1=true.
        let mut sat = SatisfiabilityChecker::new(2);
        sat.add_clause(vec![Literal::pos(0), Literal::pos(1)]);
        sat.add_clause(vec![Literal::neg(0), Literal::pos(1)]);
        let result = sat.solve();
        assert!(result.is_some());
        assert!(sat.check_assignment(&result.unwrap()));
    }

    #[test]
    fn test_sat_check_assignment() {
        let mut sat = SatisfiabilityChecker::new(2);
        sat.add_clause(vec![Literal::pos(0), Literal::pos(1)]);
        assert!(sat.check_assignment(&[false, true]));
        assert!(!sat.check_assignment(&[false, false]));
    }

    // -- TypeStateProtocol --

    #[test]
    fn test_type_state_fire_ok() {
        let mut proto = TypeStateProtocol::new("idle");
        proto.add_transition("idle", "start", "running");
        proto.add_transition("running", "stop", "idle");
        assert!(proto.fire("start").is_ok());
        assert_eq!(proto.current_state, "running");
    }

    #[test]
    fn test_type_state_fire_err() {
        let mut proto = TypeStateProtocol::new("idle");
        proto.add_transition("idle", "start", "running");
        assert!(proto.fire("stop").is_err());
    }

    #[test]
    fn test_type_state_is_reachable() {
        let mut proto = TypeStateProtocol::new("idle");
        proto.add_transition("idle", "start", "running");
        proto.add_transition("running", "pause", "paused");
        assert!(proto.is_reachable("idle", "paused"));
        assert!(!proto.is_reachable("idle", "done"));
    }

    // -- Temporal helpers --

    #[test]
    fn test_always() {
        assert!(always(&[1, 2, 3], |&x| x > 0));
        assert!(!always(&[1, -1, 3], |&x| x > 0));
    }

    #[test]
    fn test_eventually() {
        assert!(eventually(&[0, 0, 5], |&x| x > 4));
        assert!(!eventually(&[0, 0, 0], |&x| x > 4));
    }

    #[test]
    fn test_globally_after() {
        // After seeing 1, every element must be positive.
        let trace = vec![0i32, 1, 2, 3];
        assert!(globally_after(&trace, |&x| x == 1, |&x| x > 0));
    }

    // -- Lyapunov stability --

    #[test]
    fn test_lyapunov_stable() {
        let trace: Vec<Vec<f64>> = (0..5).map(|i| vec![1.0 / (i as f64 + 1.0)]).collect();
        assert!(check_lyapunov_stability(&trace));
    }

    // -- CEGAR --

    #[test]
    fn test_cegar_verified() {
        let mc = simple_model();
        let result = cegar_verify(&mc, 0, |s| s[0] >= 0.0, 10);
        assert!(matches!(result, CegarResult::Verified));
    }

    #[test]
    fn test_cegar_counterexample() {
        let mc = simple_model();
        let result = cegar_verify(&mc, 0, |s| s[0] < 2.0, 10);
        assert!(matches!(result, CegarResult::Counterexample(_)));
    }

    // -- Symbolic execution --

    #[test]
    fn test_symbolic_step_ok() {
        let mut state: SymbolicState = HashMap::new();
        state.insert("x".to_string(), IntervalDomain::new(0.0, 5.0));
        let bounds = IntervalDomain::new(0.0, 10.0);
        let next = symbolic_step(&state, "x", 1.0, &bounds);
        assert!(next.is_some());
    }

    #[test]
    fn test_symbolic_step_out_of_bounds() {
        let mut state: SymbolicState = HashMap::new();
        state.insert("x".to_string(), IntervalDomain::new(9.0, 10.0));
        let bounds = IntervalDomain::new(0.0, 5.0);
        let next = symbolic_step(&state, "x", 3.0, &bounds);
        assert!(next.is_none());
    }
}
