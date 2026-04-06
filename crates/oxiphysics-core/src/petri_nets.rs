#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Petri net modeling for concurrent discrete event systems.
//!
//! This module provides classical Petri nets, stochastic Petri nets,
//! colored Petri nets, and supporting analysis algorithms such as
//! reachability BFS and coverability trees.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// Place
// ---------------------------------------------------------------------------

/// A place in a Petri net, holding a non-negative token count.
#[derive(Debug, Clone)]
pub struct Place {
    /// Human-readable name for the place.
    pub name: String,
    /// Current number of tokens in this place.
    pub tokens: usize,
}

impl Place {
    /// Create a new place with zero tokens.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            tokens: 0,
        }
    }

    /// Create a new place with `tokens` initial tokens.
    pub fn with_tokens(name: impl Into<String>, tokens: usize) -> Self {
        Self {
            name: name.into(),
            tokens,
        }
    }
}

// ---------------------------------------------------------------------------
// Transition
// ---------------------------------------------------------------------------

/// A transition in a Petri net, optionally guarded by a predicate over the
/// token counts of the places it consumes from.
#[derive(Clone)]
pub struct Transition {
    /// Human-readable name for the transition.
    pub name: String,
    /// Optional guard function: receives current token counts for every input
    /// arc (in place-index order) and returns `true` iff the transition is
    /// allowed to fire.
    pub guard: Option<fn(&[usize]) -> bool>,
}

impl std::fmt::Debug for Transition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Transition")
            .field("name", &self.name)
            .field("guard", &self.guard.map(|_| "<fn>"))
            .finish()
    }
}

impl Transition {
    /// Create a new unguarded transition.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            guard: None,
        }
    }

    /// Create a new guarded transition.
    pub fn with_guard(name: impl Into<String>, guard: fn(&[usize]) -> bool) -> Self {
        Self {
            name: name.into(),
            guard: Some(guard),
        }
    }
}

// ---------------------------------------------------------------------------
// Arc weight type alias
// ---------------------------------------------------------------------------

/// Weight on an arc (number of tokens consumed / produced per firing).
pub type Weight = usize;

// ---------------------------------------------------------------------------
// PetriNet
// ---------------------------------------------------------------------------

/// A classical (Place/Transition) Petri net.
///
/// Arcs are stored as `(transition_idx, place_idx, weight)` triples.
#[derive(Debug, Clone)]
pub struct PetriNet {
    /// All places in this net.
    pub places: Vec<Place>,
    /// All transitions in this net.
    pub transitions: Vec<Transition>,
    /// Input arcs: `(transition, place, weight)` – tokens consumed on fire.
    pub arcs_in: Vec<(usize, usize, Weight)>,
    /// Output arcs: `(transition, place, weight)` – tokens produced on fire.
    pub arcs_out: Vec<(usize, usize, Weight)>,
}

impl PetriNet {
    /// Create an empty Petri net.
    pub fn new() -> Self {
        Self {
            places: Vec::new(),
            transitions: Vec::new(),
            arcs_in: Vec::new(),
            arcs_out: Vec::new(),
        }
    }

    /// Add a place and return its index.
    pub fn add_place(&mut self, place: Place) -> usize {
        let idx = self.places.len();
        self.places.push(place);
        idx
    }

    /// Add a transition and return its index.
    pub fn add_transition(&mut self, transition: Transition) -> usize {
        let idx = self.transitions.len();
        self.transitions.push(transition);
        idx
    }

    /// Add an input arc: transition `t` consumes `weight` tokens from place `p`.
    pub fn add_arc_in(&mut self, t: usize, p: usize, weight: Weight) {
        self.arcs_in.push((t, p, weight));
    }

    /// Add an output arc: transition `t` produces `weight` tokens into place `p`.
    pub fn add_arc_out(&mut self, t: usize, p: usize, weight: Weight) {
        self.arcs_out.push((t, p, weight));
    }

    /// Return the current marking (token counts for every place).
    pub fn marking(&self) -> Vec<usize> {
        self.places.iter().map(|p| p.tokens).collect()
    }

    /// Apply a marking vector to the net (sets token counts).
    pub fn set_marking(&mut self, marking: &[usize]) {
        for (p, &tok) in self.places.iter_mut().zip(marking.iter()) {
            p.tokens = tok;
        }
    }

    /// Return indices of all currently enabled transitions.
    ///
    /// A transition is enabled when every input arc can be satisfied and the
    /// guard (if any) returns `true`.
    pub fn enabled_transitions(&self) -> Vec<usize> {
        let mut enabled = Vec::new();
        for (t_idx, transition) in self.transitions.iter().enumerate() {
            if self.is_enabled(t_idx, transition) {
                enabled.push(t_idx);
            }
        }
        enabled
    }

    fn is_enabled(&self, t_idx: usize, transition: &Transition) -> bool {
        // Collect input token counts for this transition (for guard evaluation).
        let mut input_tokens: Vec<usize> = Vec::new();
        for &(t, p, w) in &self.arcs_in {
            if t == t_idx {
                if self.places[p].tokens < w {
                    return false;
                }
                input_tokens.push(self.places[p].tokens);
            }
        }
        if let Some(guard) = transition.guard
            && !guard(&input_tokens)
        {
            return false;
        }
        true
    }

    /// Fire transition `t_idx`, updating token counts.
    ///
    /// Returns `Ok(())` on success or `Err(String)` if the transition is not
    /// enabled.
    pub fn fire(&mut self, t_idx: usize) -> Result<(), String> {
        if t_idx >= self.transitions.len() {
            return Err(format!("transition index {} out of range", t_idx));
        }
        let transition = self.transitions[t_idx].clone();
        if !self.is_enabled(t_idx, &transition) {
            return Err(format!(
                "transition {} is not enabled",
                self.transitions[t_idx].name
            ));
        }
        // Consume tokens.
        for &(t, p, w) in &self.arcs_in {
            if t == t_idx {
                self.places[p].tokens -= w;
            }
        }
        // Produce tokens.
        for &(t, p, w) in &self.arcs_out {
            if t == t_idx {
                self.places[p].tokens += w;
            }
        }
        Ok(())
    }

    /// Returns `true` if no transition is currently enabled (deadlock).
    pub fn is_deadlocked(&self) -> bool {
        self.enabled_transitions().is_empty()
    }

    /// Sum of all tokens in the net (useful for conservation checks).
    pub fn marking_sum(&self) -> usize {
        self.places.iter().map(|p| p.tokens).sum()
    }

    /// BFS over reachable markings from the current marking.
    ///
    /// Returns the set of all reachable markings (including the initial one).
    pub fn reachability_bfs(&self) -> HashSet<Vec<usize>> {
        let mut visited: HashSet<Vec<usize>> = HashSet::new();
        let mut queue: VecDeque<Vec<usize>> = VecDeque::new();
        let initial = self.marking();
        queue.push_back(initial.clone());
        visited.insert(initial);

        // Clone the net for simulation.
        let mut sim = self.clone();

        while let Some(marking) = queue.pop_front() {
            sim.set_marking(&marking);
            for t_idx in sim.enabled_transitions() {
                let mut next_sim = sim.clone();
                let _ = next_sim.fire(t_idx);
                let next_marking = next_sim.marking();
                if !visited.contains(&next_marking) {
                    visited.insert(next_marking.clone());
                    queue.push_back(next_marking);
                }
            }
        }
        visited
    }
}

impl Default for PetriNet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CoverabilityTree
// ---------------------------------------------------------------------------

/// Omega sentinel value used in coverability analysis (represents unbounded
/// token counts).
pub const OMEGA: usize = usize::MAX;

/// A node in the coverability tree.
#[derive(Debug, Clone)]
pub struct CoverabilityNode {
    /// Marking at this node (some entries may be `OMEGA`).
    pub marking: Vec<usize>,
    /// Index of the parent node (None for root).
    pub parent: Option<usize>,
    /// Transition fired to reach this node.
    pub fired_transition: Option<usize>,
    /// Indices of child nodes.
    pub children: Vec<usize>,
}

/// Coverability tree for a Petri net.
///
/// Can determine whether the net is bounded (no `OMEGA` entries exist).
#[derive(Debug)]
pub struct CoverabilityTree {
    /// All nodes in the tree.
    pub nodes: Vec<CoverabilityNode>,
}

impl CoverabilityTree {
    /// Build the coverability tree for `net` starting from its current marking.
    pub fn build(net: &PetriNet) -> Self {
        let mut tree = CoverabilityTree { nodes: Vec::new() };
        let root_marking = net.marking();
        let root = CoverabilityNode {
            marking: root_marking,
            parent: None,
            fired_transition: None,
            children: Vec::new(),
        };
        tree.nodes.push(root);

        // Iterative BFS construction.
        let mut frontier: VecDeque<usize> = VecDeque::new();
        frontier.push_back(0);

        let mut sim = net.clone();

        while let Some(node_idx) = frontier.pop_front() {
            let current_marking = tree.nodes[node_idx].marking.clone();
            sim.set_marking(&current_marking);

            for t_idx in sim.enabled_transitions() {
                let mut next_sim = sim.clone();
                // Replace OMEGA entries with a large finite value for simulation.
                let tmp_marking: Vec<usize> = current_marking
                    .iter()
                    .map(|&x| if x == OMEGA { 1_000_000 } else { x })
                    .collect();
                next_sim.set_marking(&tmp_marking);
                if next_sim.fire(t_idx).is_err() {
                    continue;
                }
                let mut new_marking = next_sim.marking();

                // Walk ancestors to find coverability witnesses.
                let mut anc_idx = Some(node_idx);
                while let Some(ai) = anc_idx {
                    let anc_marking = &tree.nodes[ai].marking;
                    if anc_marking
                        .iter()
                        .zip(new_marking.iter())
                        .all(|(&a, &n)| a == OMEGA || a <= n)
                        && anc_marking
                            .iter()
                            .zip(new_marking.iter())
                            .any(|(&a, &n)| a != OMEGA && a < n)
                    {
                        // Strictly covered – promote to OMEGA.
                        for (m, a) in new_marking.iter_mut().zip(anc_marking.iter()) {
                            if *a != OMEGA && *a < *m {
                                *m = OMEGA;
                            }
                        }
                    }
                    anc_idx = tree.nodes[ai].parent;
                }

                // Check for duplicate markings (to terminate).
                let already_visited = tree.nodes.iter().any(|n| n.marking == new_marking);
                let child_idx = tree.nodes.len();
                tree.nodes.push(CoverabilityNode {
                    marking: new_marking,
                    parent: Some(node_idx),
                    fired_transition: Some(t_idx),
                    children: Vec::new(),
                });
                tree.nodes[node_idx].children.push(child_idx);

                if !already_visited {
                    frontier.push_back(child_idx);
                }
            }
        }

        tree
    }

    /// Returns `true` if the net is bounded (no `OMEGA` tokens in the tree).
    pub fn is_bounded(&self) -> bool {
        self.nodes
            .iter()
            .all(|n| n.marking.iter().all(|&t| t != OMEGA))
    }

    /// Returns the total number of nodes in the tree.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

// ---------------------------------------------------------------------------
// StochasticPetriNet
// ---------------------------------------------------------------------------

/// A stochastic Petri net where each transition fires with a rate drawn from
/// an exponential distribution.
#[derive(Debug, Clone)]
pub struct StochasticPetriNet {
    /// Underlying classical Petri net.
    pub net: PetriNet,
    /// Firing rates for each transition (exponential distribution parameter λ).
    pub rates: Vec<f64>,
}

impl StochasticPetriNet {
    /// Create a stochastic Petri net from an existing `PetriNet` and a rate
    /// vector (one rate per transition, in order).
    pub fn new(net: PetriNet, rates: Vec<f64>) -> Self {
        assert_eq!(
            net.transitions.len(),
            rates.len(),
            "one rate per transition required"
        );
        Self { net, rates }
    }

    /// Simulate until simulated time reaches `t_end`.
    ///
    /// Returns a vector of `(time, fired_transition_idx)` events.
    pub fn simulate_until(&mut self, t_end: f64) -> Vec<(f64, usize)> {
        let mut events = Vec::new();
        let mut t = 0.0_f64;
        use rand::RngExt as _;
        let mut rng = rand::rng();

        loop {
            let enabled = self.net.enabled_transitions();
            if enabled.is_empty() {
                break;
            }
            // Total rate of all enabled transitions.
            let total_rate: f64 = enabled.iter().map(|&i| self.rates[i]).sum();
            if total_rate == 0.0 {
                break;
            }
            // Sample time to next event (exponential distribution).
            let u: f64 = rng.random_range(0.0_f64..1.0_f64);
            let dt = -u.ln() / total_rate;
            t += dt;
            if t > t_end {
                break;
            }
            // Choose which transition fires (proportional to rates).
            let selector: f64 = rng.random_range(0.0_f64..total_rate);
            let mut cumulative = 0.0;
            let mut chosen = enabled[0];
            for &ti in &enabled {
                cumulative += self.rates[ti];
                if selector < cumulative {
                    chosen = ti;
                    break;
                }
            }
            let _ = self.net.fire(chosen);
            events.push((t, chosen));
        }
        events
    }

    /// Estimate steady-state token distribution by averaging over `n_samples`
    /// independent simulations each run to `t_end`.
    ///
    /// Returns a vector of mean token counts for each place.
    pub fn steady_state_estimate(&self, t_end: f64, n_samples: usize) -> Vec<f64> {
        let n_places = self.net.places.len();
        let mut totals = vec![0.0_f64; n_places];

        for _ in 0..n_samples {
            let mut sim = self.clone();
            let events = sim.simulate_until(t_end);
            let _ = events; // We care about the final marking.
            let marking = sim.net.marking();
            for (i, &tok) in marking.iter().enumerate() {
                totals[i] += tok as f64;
            }
        }
        totals.iter().map(|&s| s / n_samples as f64).collect()
    }
}

// ---------------------------------------------------------------------------
// ColoredPetriNet
// ---------------------------------------------------------------------------

/// A token in a colored Petri net, represented as an `f64` feature vector.
#[derive(Debug, Clone, PartialEq)]
pub struct ColoredToken {
    /// Feature vector (the "color") of the token.
    pub color: Vec<f64>,
}

impl ColoredToken {
    /// Create a token with the given color vector.
    pub fn new(color: Vec<f64>) -> Self {
        Self { color }
    }
}

/// A place in a colored Petri net, holding typed tokens.
#[derive(Debug, Clone)]
pub struct ColoredPlace {
    /// Name of the place.
    pub name: String,
    /// Tokens currently in this place.
    pub tokens: Vec<ColoredToken>,
}

impl ColoredPlace {
    /// Create an empty colored place.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            tokens: Vec::new(),
        }
    }

    /// Add a token to this place.
    pub fn add_token(&mut self, token: ColoredToken) {
        self.tokens.push(token);
    }

    /// Remove and return the first token, if any.
    pub fn take_token(&mut self) -> Option<ColoredToken> {
        if self.tokens.is_empty() {
            None
        } else {
            Some(self.tokens.remove(0))
        }
    }

    /// Number of tokens in the place.
    pub fn token_count(&self) -> usize {
        self.tokens.len()
    }
}

/// A transition in a colored Petri net.
#[derive(Debug, Clone)]
pub struct ColoredTransition {
    /// Name of the transition.
    pub name: String,
    /// Optional color guard: returns `true` iff the token is accepted.
    pub color_guard: Option<fn(&ColoredToken) -> bool>,
}

impl ColoredTransition {
    /// Create a new unguarded colored transition.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            color_guard: None,
        }
    }

    /// Create a guarded colored transition.
    pub fn with_guard(name: impl Into<String>, guard: fn(&ColoredToken) -> bool) -> Self {
        Self {
            name: name.into(),
            color_guard: Some(guard),
        }
    }
}

/// A colored Petri net where tokens carry `f64` vectors as their type.
#[derive(Debug, Clone)]
pub struct ColoredPetriNet {
    /// Places in this colored net.
    pub places: Vec<ColoredPlace>,
    /// Transitions in this colored net.
    pub transitions: Vec<ColoredTransition>,
    /// Input arcs: `(transition_idx, place_idx)`.
    pub arcs_in: Vec<(usize, usize)>,
    /// Output arcs: `(transition_idx, place_idx)`.
    pub arcs_out: Vec<(usize, usize)>,
}

impl ColoredPetriNet {
    /// Create an empty colored Petri net.
    pub fn new() -> Self {
        Self {
            places: Vec::new(),
            transitions: Vec::new(),
            arcs_in: Vec::new(),
            arcs_out: Vec::new(),
        }
    }

    /// Add a colored place and return its index.
    pub fn add_place(&mut self, place: ColoredPlace) -> usize {
        let idx = self.places.len();
        self.places.push(place);
        idx
    }

    /// Add a colored transition and return its index.
    pub fn add_transition(&mut self, transition: ColoredTransition) -> usize {
        let idx = self.transitions.len();
        self.transitions.push(transition);
        idx
    }

    /// Add an input arc.
    pub fn add_arc_in(&mut self, t: usize, p: usize) {
        self.arcs_in.push((t, p));
    }

    /// Add an output arc.
    pub fn add_arc_out(&mut self, t: usize, p: usize) {
        self.arcs_out.push((t, p));
    }

    /// Return indices of enabled transitions.
    pub fn enabled_transitions(&self) -> Vec<usize> {
        let mut enabled = Vec::new();
        'outer: for (t_idx, trans) in self.transitions.iter().enumerate() {
            for &(t, p) in &self.arcs_in {
                if t == t_idx {
                    if self.places[p].tokens.is_empty() {
                        continue 'outer;
                    }
                    if let Some(guard) = trans.color_guard
                        && !guard(&self.places[p].tokens[0])
                    {
                        continue 'outer;
                    }
                }
            }
            enabled.push(t_idx);
        }
        enabled
    }

    /// Fire transition `t_idx`, consuming one token from each input place and
    /// producing a copy into each output place.
    pub fn fire(&mut self, t_idx: usize) -> Result<(), String> {
        if t_idx >= self.transitions.len() {
            return Err(format!("transition index {} out of range", t_idx));
        }
        // Collect tokens to move.
        let mut consumed: Vec<(usize, ColoredToken)> = Vec::new();
        for &(t, p) in &self.arcs_in {
            if t == t_idx {
                let tok = self.places[p]
                    .take_token()
                    .ok_or_else(|| format!("place {} has no tokens", p))?;
                consumed.push((p, tok));
            }
        }
        // Deposit tokens into output places.
        let token_to_pass = consumed
            .first()
            .map(|(_, tok)| tok.clone())
            .unwrap_or_else(|| ColoredToken::new(vec![0.0]));
        for &(t, p) in &self.arcs_out {
            if t == t_idx {
                self.places[p].add_token(token_to_pass.clone());
            }
        }
        Ok(())
    }

    /// Total token count across all places.
    pub fn total_tokens(&self) -> usize {
        self.places.iter().map(|p| p.token_count()).sum()
    }
}

impl Default for ColoredPetriNet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PNML stub
// ---------------------------------------------------------------------------

/// Stub for parsing Petri Net Markup Language (PNML) XML files.
///
/// Returns an empty `PetriNet` – full XML parsing is outside the scope of this
/// crate (no external crate dependencies allowed).
pub fn parse_pnml(_input: &str) -> PetriNet {
    // In a full implementation this would parse the PNML XML format.
    PetriNet::new()
}

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Build a simple producer-consumer Petri net for illustration.
///
/// Returns a net with two places ("buffer" and "consumer_ready") and two
/// transitions ("produce" and "consume").
pub fn producer_consumer_net() -> PetriNet {
    let mut net = PetriNet::new();
    let buffer = net.add_place(Place::with_tokens("buffer", 0));
    let ready = net.add_place(Place::with_tokens("consumer_ready", 1));
    let produce = net.add_transition(Transition::new("produce"));
    let consume = net.add_transition(Transition::new("consume"));
    net.add_arc_in(consume, buffer, 1);
    net.add_arc_in(consume, ready, 1);
    net.add_arc_out(produce, buffer, 1);
    net.add_arc_out(consume, ready, 1);
    net
}

/// Build a simple mutual-exclusion Petri net.
///
/// Two processes compete for a single resource token.
pub fn mutex_net() -> PetriNet {
    let mut net = PetriNet::new();
    let resource = net.add_place(Place::with_tokens("resource", 1));
    let p1_idle = net.add_place(Place::with_tokens("p1_idle", 1));
    let p2_idle = net.add_place(Place::with_tokens("p2_idle", 1));
    let p1_cs = net.add_place(Place::with_tokens("p1_cs", 0));
    let p2_cs = net.add_place(Place::with_tokens("p2_cs", 0));
    let enter1 = net.add_transition(Transition::new("enter1"));
    let exit1 = net.add_transition(Transition::new("exit1"));
    let enter2 = net.add_transition(Transition::new("enter2"));
    let exit2 = net.add_transition(Transition::new("exit2"));
    net.add_arc_in(enter1, resource, 1);
    net.add_arc_in(enter1, p1_idle, 1);
    net.add_arc_out(enter1, p1_cs, 1);
    net.add_arc_in(exit1, p1_cs, 1);
    net.add_arc_out(exit1, resource, 1);
    net.add_arc_out(exit1, p1_idle, 1);
    net.add_arc_in(enter2, resource, 1);
    net.add_arc_in(enter2, p2_idle, 1);
    net.add_arc_out(enter2, p2_cs, 1);
    net.add_arc_in(exit2, p2_cs, 1);
    net.add_arc_out(exit2, resource, 1);
    net.add_arc_out(exit2, p2_idle, 1);
    net
}

// ---------------------------------------------------------------------------
// Incidence matrix helpers
// ---------------------------------------------------------------------------

/// Compute the incidence matrix of a Petri net.
///
/// Entry `[t][p]` = (weight produced by t in p) − (weight consumed by t from p).
pub fn incidence_matrix(net: &PetriNet) -> Vec<Vec<i64>> {
    let nt = net.transitions.len();
    let np = net.places.len();
    let mut matrix = vec![vec![0i64; np]; nt];
    for &(t, p, w) in &net.arcs_out {
        matrix[t][p] += w as i64;
    }
    for &(t, p, w) in &net.arcs_in {
        matrix[t][p] -= w as i64;
    }
    matrix
}

/// Compute the state equation: `m' = m + C^T * sigma` where `sigma` is a
/// firing count vector and `C` is the incidence matrix.
pub fn state_equation(net: &PetriNet, sigma: &[i64]) -> Vec<i64> {
    let matrix = incidence_matrix(net);
    let marking: Vec<i64> = net.places.iter().map(|p| p.tokens as i64).collect();
    let np = net.places.len();
    let mut result = marking;
    for (t_idx, &fires) in sigma.iter().enumerate() {
        if t_idx < matrix.len() {
            for p in 0..np {
                result[p] += matrix[t_idx][p] * fires;
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Trap and siphon detection
// ---------------------------------------------------------------------------

/// Check whether a subset of place indices forms a trap.
///
/// A trap is a set T of places such that every transition that removes a token
/// from some place in T also puts a token into some place in T.
pub fn is_trap(net: &PetriNet, place_set: &[usize]) -> bool {
    let set: HashSet<usize> = place_set.iter().copied().collect();
    for &(t, p_in, _) in &net.arcs_in {
        if set.contains(&p_in) {
            // This transition takes from the set; it must also put into the set.
            let produces_into_set = net
                .arcs_out
                .iter()
                .any(|&(tout, pout, _)| tout == t && set.contains(&pout));
            if !produces_into_set {
                return false;
            }
        }
    }
    true
}

/// Check whether a subset of place indices forms a siphon.
///
/// A siphon is a set S of places such that every transition that puts a token
/// into some place in S also removes a token from some place in S.
pub fn is_siphon(net: &PetriNet, place_set: &[usize]) -> bool {
    let set: HashSet<usize> = place_set.iter().copied().collect();
    for &(t, p_out, _) in &net.arcs_out {
        if set.contains(&p_out) {
            let takes_from_set = net
                .arcs_in
                .iter()
                .any(|&(tin, pin, _)| tin == t && set.contains(&pin));
            if !takes_from_set {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// T-invariant computation (simple linear algebra stub)
// ---------------------------------------------------------------------------

/// A T-invariant (transition invariant) is a non-negative integer vector `x`
/// such that `C * x = 0` where `C` is the incidence matrix.
///
/// This function uses a naive exhaustive search over small domains and returns
/// the first non-trivial T-invariant found (if any).
pub fn find_t_invariant(net: &PetriNet, max_count: usize) -> Option<Vec<i64>> {
    let matrix = incidence_matrix(net);
    let nt = net.transitions.len();
    let np = net.places.len();
    // Try all combinations up to max_count per transition.
    fn rec(
        matrix: &[Vec<i64>],
        nt: usize,
        np: usize,
        max_count: usize,
        idx: usize,
        current: &mut Vec<i64>,
    ) -> Option<Vec<i64>> {
        if idx == nt {
            // Check if C * current == 0.
            let mut is_trivial = true;
            for &c in current.iter() {
                if c != 0 {
                    is_trivial = false;
                    break;
                }
            }
            if is_trivial {
                return None;
            }
            for p in 0..np {
                let mut sum = 0i64;
                for t in 0..nt {
                    sum += matrix[t][p] * current[t];
                }
                if sum != 0 {
                    return None;
                }
            }
            return Some(current.clone());
        }
        for v in 0..=(max_count as i64) {
            current[idx] = v;
            if let Some(res) = rec(matrix, nt, np, max_count, idx + 1, current) {
                return Some(res);
            }
        }
        None
    }
    let mut current = vec![0i64; nt];
    rec(&matrix, nt, np, max_count, 0, &mut current)
}

// ---------------------------------------------------------------------------
// Firing sequence analysis
// ---------------------------------------------------------------------------

/// Attempt to replay a firing sequence from the current marking.
///
/// Returns `Ok(final_marking)` or `Err(step_index)` on failure.
pub fn replay_sequence(net: &mut PetriNet, sequence: &[usize]) -> Result<Vec<usize>, usize> {
    for (step, &t_idx) in sequence.iter().enumerate() {
        if net.fire(t_idx).is_err() {
            return Err(step);
        }
    }
    Ok(net.marking())
}

// ---------------------------------------------------------------------------
// Petri net metrics
// ---------------------------------------------------------------------------

/// Compute structural properties summary for a Petri net.
#[derive(Debug, Clone)]
pub struct PetriNetMetrics {
    /// Number of places.
    pub n_places: usize,
    /// Number of transitions.
    pub n_transitions: usize,
    /// Number of arcs (in + out).
    pub n_arcs: usize,
    /// Average out-degree of transitions.
    pub avg_out_degree: f64,
    /// Whether any transition has no input arcs (source transition).
    pub has_source_transition: bool,
    /// Whether any transition has no output arcs (sink transition).
    pub has_sink_transition: bool,
}

/// Compute structural metrics for the given Petri net.
pub fn compute_metrics(net: &PetriNet) -> PetriNetMetrics {
    let n_transitions = net.transitions.len();
    let total_out: usize = net
        .transitions
        .iter()
        .enumerate()
        .map(|(t, _)| net.arcs_out.iter().filter(|&&(ti, _, _)| ti == t).count())
        .sum();
    let avg_out_degree = if n_transitions == 0 {
        0.0
    } else {
        total_out as f64 / n_transitions as f64
    };
    let has_source_transition =
        (0..n_transitions).any(|t| !net.arcs_in.iter().any(|&(ti, _, _)| ti == t));
    let has_sink_transition =
        (0..n_transitions).any(|t| !net.arcs_out.iter().any(|&(ti, _, _)| ti == t));
    PetriNetMetrics {
        n_places: net.places.len(),
        n_transitions,
        n_arcs: net.arcs_in.len() + net.arcs_out.len(),
        avg_out_degree,
        has_source_transition,
        has_sink_transition,
    }
}

// ---------------------------------------------------------------------------
// Merging nets
// ---------------------------------------------------------------------------

/// Merge two Petri nets by disjoint union (all places and transitions are
/// included; no arcs are shared).
pub fn disjoint_union(a: &PetriNet, b: &PetriNet) -> PetriNet {
    let mut result = a.clone();
    let place_offset = a.places.len();
    let trans_offset = a.transitions.len();
    for p in &b.places {
        result.places.push(p.clone());
    }
    for t in &b.transitions {
        result.transitions.push(t.clone());
    }
    for &(t, p, w) in &b.arcs_in {
        result.arcs_in.push((t + trans_offset, p + place_offset, w));
    }
    for &(t, p, w) in &b.arcs_out {
        result
            .arcs_out
            .push((t + trans_offset, p + place_offset, w));
    }
    result
}

// ---------------------------------------------------------------------------
// Liveness analysis
// ---------------------------------------------------------------------------

/// Check whether transition `t_idx` is live (can always eventually fire again
/// from any reachable marking).
///
/// This is a lightweight approximation: returns `true` if `t_idx` appears in
/// the enabled set of at least one reachable marking.
pub fn is_live_approx(net: &PetriNet, t_idx: usize) -> bool {
    let reachable = net.reachability_bfs();
    let mut sim = net.clone();
    for marking in &reachable {
        sim.set_marking(marking);
        if sim.enabled_transitions().contains(&t_idx) {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Shared-resource workflow example
// ---------------------------------------------------------------------------

/// Build a workflow Petri net with a shared resource pool of size `capacity`.
pub fn workflow_with_resource(capacity: usize) -> PetriNet {
    let mut net = PetriNet::new();
    let start = net.add_place(Place::with_tokens("start", 1));
    let working = net.add_place(Place::with_tokens("working", 0));
    let done = net.add_place(Place::with_tokens("done", 0));
    let resource = net.add_place(Place::with_tokens("resource", capacity));
    let begin = net.add_transition(Transition::new("begin"));
    let finish = net.add_transition(Transition::new("finish"));
    net.add_arc_in(begin, start, 1);
    net.add_arc_in(begin, resource, 1);
    net.add_arc_out(begin, working, 1);
    net.add_arc_in(finish, working, 1);
    net.add_arc_out(finish, done, 1);
    net.add_arc_out(finish, resource, 1);
    net
}

// ---------------------------------------------------------------------------
// HashMap-based marking for large nets
// ---------------------------------------------------------------------------

/// Sparse marking: maps place names to token counts.
pub type SparseMarking = HashMap<String, usize>;

/// Convert a dense marking vector (aligned with `net.places`) into a sparse map.
pub fn dense_to_sparse(net: &PetriNet, marking: &[usize]) -> SparseMarking {
    net.places
        .iter()
        .zip(marking.iter())
        .map(|(p, &tok)| (p.name.clone(), tok))
        .collect()
}

/// Convert a sparse marking map back to a dense vector (missing places = 0).
pub fn sparse_to_dense(net: &PetriNet, sparse: &SparseMarking) -> Vec<usize> {
    net.places
        .iter()
        .map(|p| *sparse.get(&p.name).unwrap_or(&0))
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_net() -> PetriNet {
        let mut net = PetriNet::new();
        let p0 = net.add_place(Place::with_tokens("p0", 1));
        let p1 = net.add_place(Place::with_tokens("p1", 0));
        let t0 = net.add_transition(Transition::new("t0"));
        net.add_arc_in(t0, p0, 1);
        net.add_arc_out(t0, p1, 1);
        net
    }

    #[test]
    fn test_place_tokens() {
        let p = Place::with_tokens("p", 3);
        assert_eq!(p.tokens, 3);
    }

    #[test]
    fn test_transition_name() {
        let t = Transition::new("t1");
        assert_eq!(t.name, "t1");
    }

    #[test]
    fn test_add_place_transition() {
        let mut net = PetriNet::new();
        let pi = net.add_place(Place::new("p"));
        let ti = net.add_transition(Transition::new("t"));
        assert_eq!(pi, 0);
        assert_eq!(ti, 0);
    }

    #[test]
    fn test_fire_basic() {
        let mut net = simple_net();
        assert_eq!(net.marking(), vec![1, 0]);
        net.fire(0).unwrap();
        assert_eq!(net.marking(), vec![0, 1]);
    }

    #[test]
    fn test_fire_not_enabled() {
        let mut net = simple_net();
        net.fire(0).unwrap(); // consume the token
        assert!(net.fire(0).is_err());
    }

    #[test]
    fn test_enabled_transitions_empty() {
        let mut net = simple_net();
        net.fire(0).unwrap();
        assert!(net.enabled_transitions().is_empty());
    }

    #[test]
    fn test_is_deadlocked() {
        let mut net = simple_net();
        net.fire(0).unwrap();
        assert!(net.is_deadlocked());
    }

    #[test]
    fn test_marking_sum() {
        let net = simple_net();
        assert_eq!(net.marking_sum(), 1);
    }

    #[test]
    fn test_reachability_bfs() {
        let net = simple_net();
        let reachable = net.reachability_bfs();
        assert!(reachable.contains(&vec![1, 0]));
        assert!(reachable.contains(&vec![0, 1]));
        assert_eq!(reachable.len(), 2);
    }

    #[test]
    fn test_mutex_net_reachability() {
        let net = mutex_net();
        let reachable = net.reachability_bfs();
        // Mutual exclusion: p1_cs and p2_cs should never both hold a token.
        for marking in &reachable {
            // places: resource=0, p1_idle=1, p2_idle=2, p1_cs=3, p2_cs=4
            assert!(
                !(marking[3] >= 1 && marking[4] >= 1),
                "mutual exclusion violated"
            );
        }
    }

    #[test]
    fn test_producer_consumer_enabled() {
        let mut net = producer_consumer_net();
        // Initially buffer=0, consumer_ready=1; only produce is enabled.
        let enabled = net.enabled_transitions();
        assert_eq!(enabled.len(), 1); // only produce
        net.fire(enabled[0]).unwrap(); // produce a token
        assert_eq!(net.marking_sum(), 2); // buffer=1, consumer_ready=1
    }

    #[test]
    fn test_coverability_tree_bounded() {
        let net = simple_net();
        let tree = CoverabilityTree::build(&net);
        assert!(tree.is_bounded());
    }

    #[test]
    fn test_coverability_tree_unbounded() {
        // A net with a self-loop that infinitely produces tokens is unbounded.
        let mut net = PetriNet::new();
        let p0 = net.add_place(Place::with_tokens("p0", 1));
        let t0 = net.add_transition(Transition::new("t0"));
        // t0 consumes 1 from p0, produces 2.
        net.add_arc_in(t0, p0, 1);
        net.add_arc_out(t0, p0, 2);
        let tree = CoverabilityTree::build(&net);
        assert!(!tree.is_bounded());
    }

    #[test]
    fn test_stochastic_simulate() {
        let mut net = PetriNet::new();
        let p0 = net.add_place(Place::with_tokens("p0", 3));
        let p1 = net.add_place(Place::with_tokens("p1", 0));
        let t0 = net.add_transition(Transition::new("t0"));
        net.add_arc_in(t0, p0, 1);
        net.add_arc_out(t0, p1, 1);
        let mut spn = StochasticPetriNet::new(net, vec![1.0]);
        let events = spn.simulate_until(10.0);
        // Should have fired exactly 3 times (all tokens consumed).
        assert_eq!(events.len(), 3);
    }

    #[test]
    fn test_stochastic_steady_state() {
        let net = simple_net();
        let spn = StochasticPetriNet::new(net, vec![1.0]);
        let ss = spn.steady_state_estimate(100.0, 5);
        // Total tokens must be conserved in the mean.
        assert!((ss.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_colored_place_tokens() {
        let mut p = ColoredPlace::new("p");
        p.add_token(ColoredToken::new(vec![1.0, 2.0]));
        assert_eq!(p.token_count(), 1);
        let tok = p.take_token().unwrap();
        assert_eq!(tok.color, vec![1.0, 2.0]);
        assert_eq!(p.token_count(), 0);
    }

    #[test]
    fn test_colored_petri_net_fire() {
        let mut cpn = ColoredPetriNet::new();
        let p0 = cpn.add_place(ColoredPlace::new("p0"));
        let p1 = cpn.add_place(ColoredPlace::new("p1"));
        let t0 = cpn.add_transition(ColoredTransition::new("t0"));
        cpn.add_arc_in(t0, p0);
        cpn.add_arc_out(t0, p1);
        cpn.places[p0].add_token(ColoredToken::new(vec![1.0]));
        cpn.fire(t0).unwrap();
        assert_eq!(cpn.places[p0].token_count(), 0);
        assert_eq!(cpn.places[p1].token_count(), 1);
    }

    #[test]
    fn test_parse_pnml_stub() {
        let net = parse_pnml("<pnml/>");
        assert!(net.places.is_empty());
    }

    #[test]
    fn test_incidence_matrix() {
        let net = simple_net();
        let m = incidence_matrix(&net);
        // t0: -1 for p0, +1 for p1
        assert_eq!(m[0][0], -1);
        assert_eq!(m[0][1], 1);
    }

    #[test]
    fn test_state_equation() {
        let net = simple_net();
        let result = state_equation(&net, &[1]);
        assert_eq!(result, vec![0, 1]);
    }

    #[test]
    fn test_is_trap() {
        // For the simple net: {p1} is a trap (nothing removes from it).
        let simple = simple_net();
        assert!(is_trap(&simple, &[1]));
    }

    #[test]
    fn test_is_siphon() {
        let simple = simple_net();
        // {p0} is a siphon: nothing produces into p0.
        assert!(is_siphon(&simple, &[0]));
    }

    #[test]
    fn test_find_t_invariant() {
        // A loop net: t0: p0 → p1, t1: p1 → p0 should have invariant [1,1].
        let mut net = PetriNet::new();
        let p0 = net.add_place(Place::with_tokens("p0", 1));
        let p1 = net.add_place(Place::with_tokens("p1", 0));
        let t0 = net.add_transition(Transition::new("t0"));
        let t1 = net.add_transition(Transition::new("t1"));
        net.add_arc_in(t0, p0, 1);
        net.add_arc_out(t0, p1, 1);
        net.add_arc_in(t1, p1, 1);
        net.add_arc_out(t1, p0, 1);
        let inv = find_t_invariant(&net, 2);
        assert!(inv.is_some());
    }

    #[test]
    fn test_replay_sequence_ok() {
        let mut net = simple_net();
        let result = replay_sequence(&mut net, &[0]);
        assert_eq!(result, Ok(vec![0, 1]));
    }

    #[test]
    fn test_replay_sequence_fail() {
        let mut net = simple_net();
        let result = replay_sequence(&mut net, &[0, 0]);
        assert_eq!(result, Err(1));
    }

    #[test]
    fn test_compute_metrics() {
        let net = simple_net();
        let m = compute_metrics(&net);
        assert_eq!(m.n_places, 2);
        assert_eq!(m.n_transitions, 1);
        assert_eq!(m.n_arcs, 2);
    }

    #[test]
    fn test_disjoint_union() {
        let a = simple_net();
        let b = simple_net();
        let u = disjoint_union(&a, &b);
        assert_eq!(u.places.len(), 4);
        assert_eq!(u.transitions.len(), 2);
    }

    #[test]
    fn test_is_live_approx() {
        let net = simple_net();
        assert!(is_live_approx(&net, 0));
    }

    #[test]
    fn test_workflow_net() {
        let net = workflow_with_resource(2);
        assert_eq!(net.places.len(), 4);
        assert_eq!(net.transitions.len(), 2);
    }

    #[test]
    fn test_dense_sparse_roundtrip() {
        let net = simple_net();
        let marking = net.marking();
        let sparse = dense_to_sparse(&net, &marking);
        let dense = sparse_to_dense(&net, &sparse);
        assert_eq!(marking, dense);
    }

    #[test]
    fn test_guarded_transition() {
        let mut net = PetriNet::new();
        let p0 = net.add_place(Place::with_tokens("p0", 2));
        let p1 = net.add_place(Place::with_tokens("p1", 0));
        // Only fire if p0 has at least 2 tokens.
        let t0 = net.add_transition(Transition::with_guard("t0", |tokens| {
            tokens.iter().all(|&t| t >= 2)
        }));
        net.add_arc_in(t0, p0, 1);
        net.add_arc_out(t0, p1, 1);
        let enabled = net.enabled_transitions();
        assert_eq!(enabled.len(), 1);
        net.fire(0).unwrap();
        // Now p0 has 1 token, guard should fail.
        assert!(net.enabled_transitions().is_empty());
    }

    #[test]
    fn test_coverability_tree_node_count() {
        let net = simple_net();
        let tree = CoverabilityTree::build(&net);
        assert!(tree.node_count() >= 2);
    }

    #[test]
    fn test_colored_transition_guard() {
        let mut cpn = ColoredPetriNet::new();
        let p0 = cpn.add_place(ColoredPlace::new("p0"));
        let p1 = cpn.add_place(ColoredPlace::new("p1"));
        let t0 = cpn.add_transition(ColoredTransition::with_guard("t0", |tok| {
            tok.color.first().copied().unwrap_or(0.0) > 0.5
        }));
        cpn.add_arc_in(t0, p0);
        cpn.add_arc_out(t0, p1);
        // Token with low value: guard should reject.
        cpn.places[p0].add_token(ColoredToken::new(vec![0.1]));
        assert!(cpn.enabled_transitions().is_empty());
        // Token with high value: guard should accept.
        cpn.places[p0].tokens.clear();
        cpn.places[p0].add_token(ColoredToken::new(vec![0.9]));
        assert_eq!(cpn.enabled_transitions(), vec![0]);
    }
}
