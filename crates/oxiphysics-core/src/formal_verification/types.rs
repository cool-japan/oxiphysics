//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{ConcreteState, StateIdx};
use super::ordered_float_shim::NotNan;
use std::collections::{HashMap, HashSet, VecDeque};

/// Type alias for a boxed physics-state predicate used in invariants and safety properties.
pub type StatePredicate = Box<dyn Fn(&[f64]) -> bool + Send + Sync>;

/// Configuration for the CEGAR loop.
#[derive(Debug, Clone)]
pub struct CegarConfig {
    /// Maximum number of abstraction-refinement iterations (default: 20).
    pub max_iterations: usize,
    /// BMC unrolling depth per iteration (default: 10).
    pub bmc_bound: usize,
    /// Per-iteration timeout in milliseconds (0 = no timeout).
    pub timeout_ms: u64,
    /// Predicate refinement strategy.
    pub refinement_strategy: RefinementStrategy,
}
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
/// The abstract version of a `TransitionSystem` computed from a set of predicates.
///
/// Abstract states are identified by their `AbstractState` bit-vectors.
/// We maintain a list of distinct abstract states and a transition relation.
#[derive(Debug, Clone)]
pub(super) struct AbstractSystem {
    /// All discovered abstract states.
    pub(super) states: Vec<AbstractState>,
    /// Transitions: `(abstract_from_idx, abstract_to_idx)`.
    pub(super) transitions: Vec<(usize, usize)>,
    /// Map from concrete state index → abstract state index.
    pub(super) concrete_to_abstract: HashMap<StateIdx, usize>,
}
impl AbstractSystem {
    /// Return all successor abstract-state indices from abstract state `idx`.
    pub(super) fn successors(&self, idx: usize) -> Vec<usize> {
        self.transitions
            .iter()
            .filter_map(|&(f, t)| if f == idx { Some(t) } else { None })
            .collect()
    }
}
/// Strategy for predicate refinement when a spurious counterexample is found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefinementStrategy {
    /// Syntax-guided refinement: derive new predicates from the blocking
    /// coordinate differences observed along the spurious trace.
    SyntaxGuided,
    /// Craig-interpolant-style refinement: split the spurious step using
    /// threshold midpoints between consecutive abstract states.
    CraigInterpolant,
    /// Interpolation-based refinement using the abstract state bits to derive
    /// separator predicates.
    InterpolationBased,
}
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
/// A transition system wrapper that re-uses `ModelChecker` as the concrete
/// system representation.
///
/// Provides an explicit `initial_state` index so the CEGAR loop knows where
/// to start exploration.
#[derive(Debug, Clone)]
pub struct TransitionSystem {
    /// The underlying model checker / explicit-state graph.
    pub model: ModelChecker,
    /// Index of the designated initial state.
    pub initial_idx: StateIdx,
}
impl TransitionSystem {
    /// Create a transition system from a `ModelChecker` and an initial state.
    pub fn new(model: ModelChecker, initial_idx: StateIdx) -> Self {
        Self { model, initial_idx }
    }
    /// Return all successor indices of state `idx`.
    pub fn successors(&self, idx: StateIdx) -> Vec<StateIdx> {
        self.model
            .transitions
            .iter()
            .filter_map(|&(from, to)| if from == idx { Some(to) } else { None })
            .collect()
    }
    /// Number of states in the system.
    pub fn state_count(&self) -> usize {
        self.model.states.len()
    }
}
/// A named invariant with a predicate over phase-space vectors.
pub struct Invariant {
    /// Name of the invariant.
    pub name: String,
    /// Predicate: returns `true` if the state satisfies the invariant.
    pub predicate: StatePredicate,
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
/// Result of a Bounded Model Checking run on the abstract system.
#[derive(Debug)]
pub(super) enum BmcResult {
    /// No counterexample found within the bound.
    NoCounterexample,
    /// An abstract counterexample (sequence of abstract-state indices) was found.
    AbstractCex(Vec<usize>),
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
        let mut partition: Vec<usize> = (0..n).map(|_| 0).collect();
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
        loop {
            let old_partition = partition.clone();
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
/// A predicate over phase-space vectors used to build the abstract domain.
///
/// Each predicate partitions the state space into two regions: states where the
/// predicate holds (`true`) and states where it does not (`false`). The CEGAR
/// abstraction tracks which subset of predicates holds in each abstract state.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Predicate {
    /// Human-readable name for diagnostics.
    pub name: String,
    /// Variable index this predicate tests (index into the state vector).
    pub var_idx: usize,
    /// Threshold value (NaN-safe ordering).
    pub threshold: NotNan,
    /// If `true` the predicate is `state[var_idx] >= threshold`; otherwise `<`.
    pub is_ge: bool,
}
impl Predicate {
    /// Create a `>= threshold` predicate on `state[var_idx]`.
    pub fn ge(name: impl Into<String>, var_idx: usize, threshold: f64) -> Self {
        Self {
            name: name.into(),
            var_idx,
            threshold: NotNan::new(threshold),
            is_ge: true,
        }
    }
    /// Create a `< threshold` predicate on `state[var_idx]`.
    pub fn lt(name: impl Into<String>, var_idx: usize, threshold: f64) -> Self {
        Self {
            name: name.into(),
            var_idx,
            threshold: NotNan::new(threshold),
            is_ge: false,
        }
    }
    /// Evaluate the predicate on a concrete state vector.
    pub fn eval(&self, state: &[f64]) -> bool {
        let v = state.get(self.var_idx).copied().unwrap_or(0.0);
        let t = self.threshold.into_inner();
        if self.is_ge { v >= t } else { v < t }
    }
}
/// Safety property: the invariant that must hold in every reachable state.
pub struct SafetyProperty {
    /// Human-readable label.
    pub label: String,
    /// The predicate that must hold in every reachable state.
    pub predicate: StatePredicate,
}
impl SafetyProperty {
    /// Create a new safety property.
    pub fn new(
        label: impl Into<String>,
        predicate: impl Fn(&[f64]) -> bool + Send + Sync + 'static,
    ) -> Self {
        Self {
            label: label.into(),
            predicate: Box::new(predicate),
        }
    }
    /// Evaluate the property on a concrete state.
    pub fn check(&self, state: &[f64]) -> bool {
        (self.predicate)(state)
    }
}
/// An abstract state: a bit-vector indicating which predicates hold.
///
/// `bits[i]` is `true` if the `i`-th predicate holds in this abstract state.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AbstractState {
    /// Predicate truth values — length equals the number of active predicates.
    pub bits: Vec<bool>,
    /// The concrete state index in the `ModelChecker` this abstract state was
    /// derived from (`None` if synthesized during BMC).
    pub concrete_idx: Option<StateIdx>,
}
impl AbstractState {
    /// Build an abstract state from a concrete phase-space vector and predicate list.
    pub fn from_concrete(state: &[f64], predicates: &[Predicate]) -> Self {
        Self {
            bits: predicates.iter().map(|p| p.eval(state)).collect(),
            concrete_idx: None,
        }
    }
}
/// Result of checking whether an abstract counterexample is concrete.
#[derive(Debug)]
pub(super) enum Concreteness {
    /// The abstract trace has a concrete witness — a real counterexample.
    Concrete(Vec<ConcreteState>),
    /// The abstract trace is spurious — no concrete path corresponds to it.
    Spurious,
}
/// Result of the CEGAR verification procedure.
#[derive(Debug, Clone)]
pub enum CegarResult {
    /// The property holds for all reachable states within the BMC bound.
    Verified,
    /// A concrete counterexample was found.
    Violated {
        /// The concrete execution trace leading to a bad state.
        trace: Vec<ConcreteState>,
    },
    /// The analysis terminated without a definitive answer.
    Unknown {
        /// Reason for the unknown result.
        reason: String,
    },
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
        loop {
            let mut propagated = false;
            for clause in &self.clauses {
                let (unassigned, satisfied) = self.clause_status(clause, assignment);
                if satisfied {
                    continue;
                }
                if unassigned.is_empty() {
                    return false;
                }
                if unassigned.len() == 1 {
                    let lit = unassigned[0];
                    assignment[lit.var] = Some(lit.positive);
                    propagated = true;
                }
            }
            if !propagated {
                break;
            }
        }
        if self.all_satisfied(assignment) {
            return true;
        }
        if self.has_conflict(assignment) {
            return false;
        }
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
