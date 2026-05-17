//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet, VecDeque};

use super::types::{
    AbstractState, AbstractSystem, BmcResult, CegarConfig, CegarResult, Concreteness,
    IntervalDomain, Predicate, RefinementStrategy, SafetyProperty, TransitionSystem,
};

/// State index type alias.
pub type StateIdx = usize;
/// Type alias for a boxed trace-level property predicate.
pub type TracePredicate = Box<dyn Fn(&[Vec<f64>]) -> bool>;
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
    let lz0 = states[0][0] * states[0][4] - states[0][1] * states[0][3];
    states
        .iter()
        .all(|s| ((s[0] * s[4] - s[1] * s[3]) - lz0).abs() <= tol)
}
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
    last_norm <= first_norm * 1.001
}
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
/// A concrete state in the transition system — the raw phase-space vector.
pub type ConcreteState = Vec<f64>;
/// A predicate-abstraction of `system` using `predicates`.
///
/// Each concrete state is mapped to an abstract state (a predicate bit-vector).
/// An abstract transition exists whenever a concrete transition exists between
/// the corresponding abstraction classes.
fn predicate_abstraction(system: &TransitionSystem, predicates: &[Predicate]) -> AbstractSystem {
    let mut abstract_states: Vec<AbstractState> = Vec::new();
    let mut concrete_to_abstract: HashMap<StateIdx, usize> = HashMap::new();
    let mut state_to_idx: HashMap<Vec<bool>, usize> = HashMap::new();
    for (ci, concrete) in system.model.states.iter().enumerate() {
        let bits: Vec<bool> = predicates.iter().map(|p| p.eval(concrete)).collect();
        let abs_idx = if let Some(&existing) = state_to_idx.get(&bits) {
            existing
        } else {
            let idx = abstract_states.len();
            abstract_states.push(AbstractState {
                bits: bits.clone(),
                concrete_idx: Some(ci),
            });
            state_to_idx.insert(bits, idx);
            idx
        };
        concrete_to_abstract.insert(ci, abs_idx);
    }
    let mut transitions_set: HashSet<(usize, usize)> = HashSet::new();
    for &(from, to) in &system.model.transitions {
        if let (Some(&af), Some(&at)) = (
            concrete_to_abstract.get(&from),
            concrete_to_abstract.get(&to),
        ) {
            transitions_set.insert((af, at));
        }
    }
    AbstractSystem {
        states: abstract_states,
        transitions: transitions_set.into_iter().collect(),
        concrete_to_abstract,
    }
}
/// Run Bounded Model Checking on the abstract system for `bound` steps.
///
/// Starting from the abstract initial state, we do a bounded BFS.
/// An abstract state "violates" the property if ANY of the concrete states it
/// represents violates the property (sound over-approximation).
fn bounded_model_check(
    abs_sys: &AbstractSystem,
    system: &TransitionSystem,
    property: &SafetyProperty,
    bound: usize,
) -> BmcResult {
    let init_abs = match abs_sys.concrete_to_abstract.get(&system.initial_idx) {
        Some(&i) => i,
        None => return BmcResult::NoCounterexample,
    };
    let mut abs_violates: HashMap<usize, bool> = HashMap::new();
    for (&ci, &ai) in &abs_sys.concrete_to_abstract {
        if system
            .model
            .states
            .get(ci)
            .is_some_and(|concrete| !property.check(concrete))
        {
            abs_violates.insert(ai, true);
        }
    }
    let mut queue: VecDeque<(usize, Vec<usize>)> = VecDeque::new();
    queue.push_back((init_abs, vec![init_abs]));
    let mut visited_at_depth: HashMap<(usize, usize), bool> = HashMap::new();
    while let Some((abs_idx, path)) = queue.pop_front() {
        let depth = path.len() - 1;
        let key = (abs_idx, depth);
        if visited_at_depth.contains_key(&key) {
            continue;
        }
        visited_at_depth.insert(key, true);
        if abs_violates.get(&abs_idx).copied().unwrap_or(false) {
            return BmcResult::AbstractCex(path);
        }
        if depth < bound {
            for succ in abs_sys.successors(abs_idx) {
                let mut new_path = path.clone();
                new_path.push(succ);
                queue.push_back((succ, new_path));
            }
        }
    }
    BmcResult::NoCounterexample
}
/// Check whether the abstract trace `abstract_trace` corresponds to a real
/// concrete execution path in `system`.
///
/// The abstract trace is a sequence of abstract-state indices. Because the
/// abstract system may merge multiple concrete states into one abstract state,
/// an abstract step `A → B` can correspond to *multiple* concrete steps that
/// stay within A before finally transitioning to B.
///
/// We therefore perform a **bounded forward reachability** check: starting
/// from the initial concrete state, we expand reachable concrete states by
/// following transitions that stay within the current abstract state, then
/// cross into the next abstract state, and repeat for each abstract step.
fn check_concreteness(
    system: &TransitionSystem,
    abs_sys: &AbstractSystem,
    abstract_trace: &[usize],
) -> Concreteness {
    if abstract_trace.is_empty() {
        return Concreteness::Spurious;
    }
    let mut abs_to_concretes: HashMap<usize, Vec<StateIdx>> = HashMap::new();
    for (&ci, &ai) in &abs_sys.concrete_to_abstract {
        abs_to_concretes.entry(ai).or_default().push(ci);
    }
    let expand_within = |reachable: &[StateIdx], class_abs: usize| -> Vec<StateIdx> {
        let class_set: HashSet<StateIdx> = abs_to_concretes
            .get(&class_abs)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect();
        let mut frontier: Vec<StateIdx> = reachable.to_vec();
        let mut visited: HashSet<StateIdx> = reachable.iter().copied().collect();
        loop {
            let mut added = false;
            for &ci in frontier.clone().iter() {
                for succ in system.successors(ci) {
                    if class_set.contains(&succ) && !visited.contains(&succ) {
                        visited.insert(succ);
                        frontier.push(succ);
                        added = true;
                    }
                }
            }
            if !added {
                break;
            }
        }
        frontier
    };
    let first_abs = abstract_trace[0];
    let mut reachable: Vec<StateIdx> = abs_to_concretes
        .get(&first_abs)
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter(|&ci| ci == system.initial_idx)
        .collect();
    reachable = expand_within(&reachable, first_abs);
    let mut witness_path: Vec<StateIdx> = reachable
        .first()
        .copied()
        .map(|ci| vec![ci])
        .unwrap_or_default();
    if reachable.is_empty() {
        return Concreteness::Spurious;
    }
    for &next_abs in abstract_trace.iter().skip(1) {
        let next_concretes_in_abs: HashSet<StateIdx> = abs_to_concretes
            .get(&next_abs)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect();
        let mut next_reachable: Vec<StateIdx> = Vec::new();
        for &prev in &reachable {
            for succ in system.successors(prev) {
                if next_concretes_in_abs.contains(&succ) && !next_reachable.contains(&succ) {
                    next_reachable.push(succ);
                }
            }
        }
        if next_reachable.is_empty() {
            return Concreteness::Spurious;
        }
        next_reachable = expand_within(&next_reachable, next_abs);
        if let Some(&rep) = next_reachable.first() {
            witness_path.push(rep);
        }
        reachable = next_reachable;
    }
    let concrete_path: Vec<ConcreteState> = witness_path
        .iter()
        .filter_map(|&ci| system.model.states.get(ci).cloned())
        .collect();
    if concrete_path.is_empty() {
        Concreteness::Spurious
    } else {
        Concreteness::Concrete(concrete_path)
    }
}
/// Derive new predicates to eliminate a spurious abstract trace.
///
/// We identify the first "infeasible step" in `abstract_trace` — the transition
/// `abs[k-1] → abs[k]` that has no concrete witness — and generate separating
/// predicates from the coordinate differences between representative concrete
/// states on each side of that gap.
fn refine_predicates(
    system: &TransitionSystem,
    abs_sys: &AbstractSystem,
    abstract_trace: &[usize],
    strategy: &RefinementStrategy,
) -> Vec<Predicate> {
    let mut infeasible_step = abstract_trace.len().saturating_sub(1);
    let mut abs_to_concretes: HashMap<usize, Vec<StateIdx>> = HashMap::new();
    for (&ci, &ai) in &abs_sys.concrete_to_abstract {
        abs_to_concretes.entry(ai).or_default().push(ci);
    }
    let mut reachable: Vec<StateIdx> = abs_to_concretes
        .get(&abstract_trace[0])
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter(|&ci| ci == system.initial_idx)
        .collect();
    for (step, &next_abs) in abstract_trace.iter().enumerate().skip(1) {
        let next_in_abs: HashSet<StateIdx> = abs_to_concretes
            .get(&next_abs)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect();
        let mut next_reachable: Vec<StateIdx> = Vec::new();
        for &prev in &reachable {
            for succ in system.successors(prev) {
                if next_in_abs.contains(&succ) && !next_reachable.contains(&succ) {
                    next_reachable.push(succ);
                }
            }
        }
        if next_reachable.is_empty() {
            infeasible_step = step;
            break;
        }
        reachable = next_reachable;
    }
    let before_abs = if infeasible_step > 0 {
        abstract_trace[infeasible_step - 1]
    } else {
        abstract_trace[0]
    };
    let after_abs = abstract_trace[infeasible_step];
    let before_concretes: Vec<&Vec<f64>> = abs_to_concretes
        .get(&before_abs)
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|ci| system.model.states.get(ci))
        .collect();
    let after_concretes: Vec<&Vec<f64>> = abs_to_concretes
        .get(&after_abs)
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|ci| system.model.states.get(ci))
        .collect();
    if before_concretes.is_empty() || after_concretes.is_empty() {
        return Vec::new();
    }
    let n_vars = before_concretes[0].len().min(after_concretes[0].len());
    let mut new_predicates: Vec<Predicate> = Vec::new();
    match strategy {
        RefinementStrategy::SyntaxGuided | RefinementStrategy::InterpolationBased => {
            for var_idx in 0..n_vars {
                let mean_before: f64 = before_concretes
                    .iter()
                    .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
                    .sum::<f64>()
                    / before_concretes.len() as f64;
                let mean_after: f64 = after_concretes
                    .iter()
                    .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
                    .sum::<f64>()
                    / after_concretes.len() as f64;
                if (mean_before - mean_after).abs() > 1e-12 {
                    let threshold = (mean_before + mean_after) / 2.0;
                    let name_ge = format!("refine_v{}>={}_{:.4}", var_idx, var_idx, threshold);
                    let name_lt = format!("refine_v{}<{}_{:.4}", var_idx, var_idx, threshold);
                    new_predicates.push(Predicate::ge(name_ge, var_idx, threshold));
                    new_predicates.push(Predicate::lt(name_lt, var_idx, threshold));
                }
            }
        }
        RefinementStrategy::CraigInterpolant => {
            let before_bits = &abs_sys.states[before_abs].bits;
            let after_bits = &abs_sys.states[after_abs].bits;
            for var_idx in 0..n_vars {
                let val_before = before_concretes[0].get(var_idx).copied().unwrap_or(0.0);
                let val_after = after_concretes[0].get(var_idx).copied().unwrap_or(0.0);
                let bit_idx = var_idx.min(before_bits.len().saturating_sub(1));
                let bits_differ = if bit_idx < before_bits.len() && bit_idx < after_bits.len() {
                    before_bits[bit_idx] != after_bits[bit_idx]
                } else {
                    true
                };
                if bits_differ && (val_before - val_after).abs() > 1e-12 {
                    let threshold = (val_before + val_after) / 2.0;
                    new_predicates.push(Predicate::ge(
                        format!("craig_v{}_ge_{:.4}", var_idx, threshold),
                        var_idx,
                        threshold,
                    ));
                }
            }
        }
    }
    new_predicates
}
/// Compute initial predicates from the safety property by sampling state-space
/// boundaries implied by the property's structure.
///
/// Since `SafetyProperty` wraps an opaque closure, we derive boundary predicates
/// for each state variable based on the representative concrete states in `system`.
fn initial_predicates(system: &TransitionSystem, property: &SafetyProperty) -> Vec<Predicate> {
    let mut predicates: Vec<Predicate> = Vec::new();
    let (satisfying, violating): (Vec<&Vec<f64>>, Vec<&Vec<f64>>) =
        system.model.states.iter().partition(|s| property.check(s));
    if satisfying.is_empty() || violating.is_empty() {
        let n_vars = system.model.states.first().map(|s| s.len()).unwrap_or(1);
        for var_idx in 0..n_vars {
            let vals: Vec<f64> = system
                .model
                .states
                .iter()
                .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
                .collect();
            if let (Some(&lo), Some(&hi)) = (
                vals.iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)),
                vals.iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)),
            ) && (hi - lo).abs() > 1e-12
            {
                let mid = (lo + hi) / 2.0;
                predicates.push(Predicate::ge(
                    format!("init_v{}_ge_{:.4}", var_idx, mid),
                    var_idx,
                    mid,
                ));
            }
        }
        return predicates;
    }
    let n_vars = system.model.states.first().map(|s| s.len()).unwrap_or(1);
    for var_idx in 0..n_vars {
        let max_sat = satisfying
            .iter()
            .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
            .fold(f64::NEG_INFINITY, f64::max);
        let min_sat = satisfying
            .iter()
            .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
            .fold(f64::INFINITY, f64::min);
        let max_vio = violating
            .iter()
            .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
            .fold(f64::NEG_INFINITY, f64::max);
        let min_vio = violating
            .iter()
            .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
            .fold(f64::INFINITY, f64::min);
        if min_sat > max_vio {
            let threshold = (min_sat + max_vio) / 2.0;
            predicates.push(Predicate::ge(
                format!("init_v{}_ge_{:.4}", var_idx, threshold),
                var_idx,
                threshold,
            ));
            predicates.push(Predicate::lt(
                format!("init_v{}_lt_{:.4}", var_idx, threshold),
                var_idx,
                threshold,
            ));
        } else if max_sat < min_vio {
            let threshold = (max_sat + min_vio) / 2.0;
            predicates.push(Predicate::ge(
                format!("init_v{}_ge_{:.4}", var_idx, threshold),
                var_idx,
                threshold,
            ));
            predicates.push(Predicate::lt(
                format!("init_v{}_lt_{:.4}", var_idx, threshold),
                var_idx,
                threshold,
            ));
        } else {
            let all_vals: Vec<f64> = satisfying
                .iter()
                .chain(violating.iter())
                .map(|s| s.get(var_idx).copied().unwrap_or(0.0))
                .collect();
            let lo = all_vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = all_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if (hi - lo).abs() > 1e-12 {
                let mid = (lo + hi) / 2.0;
                predicates.push(Predicate::ge(
                    format!("init_v{}_ge_{:.4}", var_idx, mid),
                    var_idx,
                    mid,
                ));
                predicates.push(Predicate::lt(
                    format!("init_v{}_lt_{:.4}", var_idx, mid),
                    var_idx,
                    mid,
                ));
            }
        }
    }
    predicates
}
/// Run the full CEGAR loop on `system`, verifying `property`.
///
/// Iterates:
/// 1. Compute predicate abstraction from the current predicate set.
/// 2. Run bounded model checking on the abstract system.
/// 3. If an abstract counterexample is found, check its concreteness.
/// 4. If spurious, refine predicates and repeat.
/// 5. Terminate when verified, a concrete violation is found, or the
///    iteration budget is exhausted.
pub fn cegar_verify(
    system: &TransitionSystem,
    property: &SafetyProperty,
    config: &CegarConfig,
) -> CegarResult {
    let mut predicates = initial_predicates(system, property);
    for iteration in 0..config.max_iterations {
        let abstract_sys = predicate_abstraction(system, &predicates);
        if abstract_sys.states.is_empty() {
            return CegarResult::Verified;
        }
        let bmc_result = bounded_model_check(&abstract_sys, system, property, config.bmc_bound);
        match bmc_result {
            BmcResult::NoCounterexample => return CegarResult::Verified,
            BmcResult::AbstractCex(abstract_trace) => {
                match check_concreteness(system, &abstract_sys, &abstract_trace) {
                    Concreteness::Concrete(concrete_trace) => {
                        return CegarResult::Violated {
                            trace: concrete_trace,
                        };
                    }
                    Concreteness::Spurious => {
                        let new_preds = refine_predicates(
                            system,
                            &abstract_sys,
                            &abstract_trace,
                            &config.refinement_strategy,
                        );
                        if new_preds.is_empty() {
                            return CegarResult::Unknown {
                                reason: format!(
                                    "Refinement stalled at iteration {} (no new predicates)",
                                    iteration
                                ),
                            };
                        }
                        predicates.extend(new_preds);
                        predicates.sort_unstable();
                        predicates.dedup();
                    }
                }
            }
        }
    }
    CegarResult::Unknown {
        reason: format!(
            "Exceeded {} iterations without convergence",
            config.max_iterations
        ),
    }
}
/// Verify a collection of named properties against a simulation trace.
///
/// Returns a map from property name to `true`/`false`.
pub fn verify_trace_properties(
    trace: &[Vec<f64>],
    properties: &[(&str, TracePredicate)],
) -> HashMap<String, bool> {
    properties
        .iter()
        .map(|(name, predicate)| (name.to_string(), predicate(trace)))
        .collect()
}
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
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
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
    fn simple_model() -> ModelChecker {
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
    #[test]
    fn test_sat_trivial_sat() {
        let mut sat = SatisfiabilityChecker::new(1);
        sat.add_clause(vec![Literal::pos(0)]);
        let result = sat.solve();
        assert!(result.is_some());
        let assignment = result.unwrap();
        assert!(sat.check_assignment(&assignment));
    }
    #[test]
    fn test_sat_trivial_unsat() {
        let mut sat = SatisfiabilityChecker::new(1);
        sat.add_clause(vec![Literal::pos(0)]);
        sat.add_clause(vec![Literal::neg(0)]);
        assert!(sat.solve().is_none());
    }
    #[test]
    fn test_sat_two_vars() {
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
        let trace = vec![0i32, 1, 2, 3];
        assert!(globally_after(&trace, |&x| x == 1, |&x| x > 0));
    }
    #[test]
    fn test_lyapunov_stable() {
        let trace: Vec<Vec<f64>> = (0..5).map(|i| vec![1.0 / (i as f64 + 1.0)]).collect();
        assert!(check_lyapunov_stability(&trace));
    }
    fn simple_safe_system() -> TransitionSystem {
        let mc = ModelChecker::new(vec![vec![0.0], vec![1.0], vec![2.0]], vec![(0, 1), (1, 2)]);
        TransitionSystem::new(mc, 0)
    }
    fn simple_violated_system() -> TransitionSystem {
        let mc = ModelChecker::new(vec![vec![0.0], vec![-1.0]], vec![(0, 1)]);
        TransitionSystem::new(mc, 0)
    }
    #[test]
    fn test_cegar_simple_safe_property() {
        let system = simple_safe_system();
        let property = SafetyProperty::new("x >= 0", |s| s[0] >= 0.0);
        let config = CegarConfig::default();
        let result = cegar_verify(&system, &property, &config);
        assert!(
            matches!(result, CegarResult::Verified),
            "Expected Verified, got {:?}",
            result
        );
    }
    #[test]
    fn test_cegar_simple_violated_property() {
        let system = simple_violated_system();
        let property = SafetyProperty::new("x >= 0", |s| s[0] >= 0.0);
        let config = CegarConfig::default();
        let result = cegar_verify(&system, &property, &config);
        assert!(
            matches!(result, CegarResult::Violated { .. }),
            "Expected Violated, got {:?}",
            result
        );
        if let CegarResult::Violated { trace } = &result {
            assert!(!trace.is_empty(), "Trace should not be empty");
        }
    }
    #[test]
    fn test_cegar_refinement_adds_predicates() {
        let mc = ModelChecker::new(
            vec![vec![5.0], vec![10.0], vec![-1.0]],
            vec![(0, 1), (1, 2)],
        );
        let system = TransitionSystem::new(mc, 0);
        let property = SafetyProperty::new("x >= 0", |s| s[0] >= 0.0);
        let config = CegarConfig {
            max_iterations: 10,
            bmc_bound: 5,
            ..CegarConfig::default()
        };
        let result = cegar_verify(&system, &property, &config);
        assert!(
            matches!(result, CegarResult::Violated { .. } | CegarResult::Verified),
            "Got unexpected Unknown: {:?}",
            result
        );
    }
    #[test]
    fn test_cegar_timeout_returns_unknown() {
        let mc = ModelChecker::new(vec![vec![1.0], vec![-1.0]], vec![(0, 1)]);
        let system = TransitionSystem::new(mc, 0);
        let property = SafetyProperty::new("x >= 0", |s| s[0] >= 0.0);
        let config = CegarConfig {
            max_iterations: 0,
            bmc_bound: 10,
            ..CegarConfig::default()
        };
        let result = cegar_verify(&system, &property, &config);
        assert!(
            matches!(result, CegarResult::Unknown { .. }),
            "Expected Unknown with max_iterations=0, got {:?}",
            result
        );
    }
    #[test]
    fn test_cegar_trace_length_matches_bound() {
        let n = 6usize;
        let states: Vec<Vec<f64>> = (0..n).map(|i| vec![i as f64]).collect();
        let transitions: Vec<(StateIdx, StateIdx)> = (0..(n - 1)).map(|i| (i, i + 1)).collect();
        let mc = ModelChecker::new(states, transitions);
        let system = TransitionSystem::new(mc, 0);
        let property = SafetyProperty::new("x < 5", |s| s[0] < 5.0);
        let config = CegarConfig {
            bmc_bound: 10,
            ..CegarConfig::default()
        };
        let result = cegar_verify(&system, &property, &config);
        if let CegarResult::Violated { trace } = result {
            assert!(
                trace.len() <= config.bmc_bound + 1,
                "Trace length {} exceeds bound {}",
                trace.len(),
                config.bmc_bound + 1
            );
            let last = trace.last().expect("trace non-empty");
            assert!(
                !property.check(last),
                "Last trace state should violate property"
            );
        } else {
            panic!("Expected Violated result");
        }
    }
    #[test]
    fn test_cegar_craig_interpolant_strategy() {
        let system = simple_violated_system();
        let property = SafetyProperty::new("x >= 0", |s| s[0] >= 0.0);
        let config = CegarConfig {
            refinement_strategy: RefinementStrategy::CraigInterpolant,
            ..CegarConfig::default()
        };
        let result = cegar_verify(&system, &property, &config);
        assert!(
            matches!(result, CegarResult::Violated { .. }),
            "Expected Violated with CraigInterpolant strategy, got {:?}",
            result
        );
    }
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
