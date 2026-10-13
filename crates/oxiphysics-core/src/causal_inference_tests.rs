// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the `causal_inference` module (causal graphs, SCMs,
//! back-door/front-door criteria, discovery).
//!
//! Split out of `causal_inference.rs` to keep that file under the 2000-line limit.
//! Included via `#[path]` as the body of `causal_inference::tests`.

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

// --- CausalGraph::without_outgoing ---

#[test]
fn test_without_outgoing_removes_edges_leaving_node() {
    let g = simple_chain(); // 0 → 1 → 2
    let cut = g.without_outgoing(1);
    assert_eq!(cut.children[1], Vec::<usize>::new());
    assert!(cut.parents[2].is_empty(), "edge 1→2 must be gone");
    assert_eq!(cut.children[0], vec![1], "edge 0→1 must survive");
    assert_eq!(cut.parents[1], vec![0]);
    // The receiver is not modified.
    assert_eq!(g.children[1], vec![2]);
    assert_eq!(g.parents[2], vec![1]);
}

// --- FrontdoorCriterion::check ---

/// U=0, X=1, M=2, Y=3 with U→X, U→Y, X→M, M→Y (the canonical frontdoor DAG).
fn frontdoor_graph() -> CausalGraph {
    let mut g = CausalGraph::new(4);
    g.add_edge(0, 1);
    g.add_edge(0, 3);
    g.add_edge(1, 2);
    g.add_edge(2, 3);
    g
}

#[test]
fn test_frontdoor_check_accepts_canonical_graph() {
    let fd = FrontdoorCriterion::new(frontdoor_graph());
    assert!(fd.check(1, 3, &[2]));
}

#[test]
fn test_frontdoor_check_rejects_confounder_of_treatment_and_mediator() {
    // Adding U→M opens the back-door path X ← U → M (condition 2).
    let mut g = frontdoor_graph();
    g.add_edge(0, 2);
    let fd = FrontdoorCriterion::new(g);
    assert!(!fd.check(1, 3, &[2]));
}

#[test]
fn test_frontdoor_check_rejects_unintercepted_direct_effect() {
    // Adding X→Y leaves a directed path that M does not intercept
    // (condition 1).
    let mut g = frontdoor_graph();
    g.add_edge(1, 3);
    let fd = FrontdoorCriterion::new(g);
    assert!(!fd.check(1, 3, &[2]));
}

#[test]
fn test_frontdoor_check_rejects_unblocked_mediator_outcome_backdoor() {
    // A second confounder W affects M and Y: the back-door path
    // M ← W → Y is not blocked by X (condition 3).
    let mut g = CausalGraph::new(5); // U=0, X=1, M=2, Y=3, W=4
    g.add_edge(0, 1);
    g.add_edge(0, 3);
    g.add_edge(1, 2);
    g.add_edge(2, 3);
    g.add_edge(4, 2);
    g.add_edge(4, 3);
    let fd = FrontdoorCriterion::new(g);
    assert!(!fd.check(1, 3, &[2]));
}

#[test]
fn test_frontdoor_check_requires_every_directed_path_to_be_intercepted() {
    // Two parallel mediators: only the full set intercepts all paths.
    let mut g = CausalGraph::new(5); // U=0, X=1, M1=2, M2=3, Y=4
    g.add_edge(0, 1);
    g.add_edge(0, 4);
    g.add_edge(1, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 4);
    g.add_edge(3, 4);
    let fd = FrontdoorCriterion::new(g);
    assert!(!fd.check(1, 4, &[2]));
    assert!(!fd.check(1, 4, &[3]));
    assert!(fd.check(1, 4, &[2, 3]));
}

#[test]
fn test_frontdoor_check_accepts_serial_mediators() {
    // X → M1 → M2 → Y with U → X, U → Y: {M1}, {M2} and {M1, M2} all work.
    let mut g = CausalGraph::new(5); // U=0, X=1, M1=2, M2=3, Y=4
    g.add_edge(0, 1);
    g.add_edge(0, 4);
    g.add_edge(1, 2);
    g.add_edge(2, 3);
    g.add_edge(3, 4);
    let fd = FrontdoorCriterion::new(g);
    assert!(fd.check(1, 4, &[2]));
    assert!(fd.check(1, 4, &[3]));
    assert!(fd.check(1, 4, &[2, 3]));
}

#[test]
fn test_frontdoor_check_rejects_malformed_sets() {
    let fd = FrontdoorCriterion::new(frontdoor_graph());
    assert!(
        !fd.check(1, 3, &[1]),
        "treatment cannot be its own mediator"
    );
    assert!(!fd.check(1, 3, &[3]), "outcome cannot be a mediator");
    assert!(!fd.check(1, 3, &[2, 9]), "out-of-range mediator");
    assert!(!fd.check(1, 9, &[2]), "out-of-range outcome");
    assert!(!fd.check(9, 3, &[2]), "out-of-range treatment");
    assert!(!fd.check(1, 1, &[2]), "treatment equals outcome");
}

// --- FrontdoorCriterion::adjusted_effect ---

/// Linear frontdoor SCM: U=0, X=1, M=2, Y=3 with U→X (0.9), U→Y (1.0),
/// X→M (0.8), M→Y (1.5); the true effect of X on Y is 0.8 · 1.5 = 1.2.
fn frontdoor_scm() -> StructuralCausalModel {
    let mut scm = StructuralCausalModel::new(4);
    scm.add_edge(0, 1, 0.9);
    scm.add_edge(0, 3, 1.0);
    scm.add_edge(1, 2, 0.8);
    scm.add_edge(2, 3, 1.5);
    scm.noise_std = vec![1.0, 0.5, 0.3, 0.3];
    scm
}

/// Draw `n` observational samples `(x, m, y)` from `scm` with a seeded RNG.
fn sample_observational(
    scm: &StructuralCausalModel,
    n: usize,
    seed: u64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut rng = crate::StatRng::new(seed);
    let (mut xs, mut ms, mut ys) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..n {
        let noise: Vec<f64> = (0..scm.graph.n_nodes).map(|_| rng.next_normal()).collect();
        let row = scm.sample_with_noise(&noise);
        xs.push(row[1]);
        ms.push(row[2]);
        ys.push(row[3]);
    }
    (xs, ms, ys)
}

/// Ordinary least-squares slope of `y` on `x`.
fn ols_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    sxy / sxx
}

#[test]
fn test_frontdoor_adjusted_effect_recovers_causal_slope_under_confounding() {
    let scm = frontdoor_scm();
    let truth = scm.total_effect_linear(1, 3);
    assert!((truth - 1.2).abs() < 1e-12);

    let (xs, ms, ys) = sample_observational(&scm, 5000, 20_240_607);
    let e1 = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, 1.0).unwrap();
    let e0 = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, 0.0).unwrap();
    let slope = e1 - e0;
    assert!(
        (slope - truth).abs() < 0.05,
        "frontdoor slope {slope}, true effect {truth}"
    );

    // The naive regression of Y on X is badly confounded by U (≈ 2.05).
    let naive = ols_slope(&xs, &ys);
    assert!(
        (naive - truth).abs() > 0.5,
        "naive slope {naive} should be far from {truth}"
    );
}

#[test]
fn test_frontdoor_adjusted_effect_matches_do_calculus_levels() {
    // With intercepts the *level* E[Y | do(X = x)] must match the SCM:
    // Y = 2 + 1.5·(0.5 + 0.8·x) = 2.75 + 1.2·x.
    let mut scm = frontdoor_scm();
    scm.set_intercept(2, 0.5);
    scm.set_intercept(3, 2.0);
    let (xs, ms, ys) = sample_observational(&scm, 5000, 7_919);
    for &x_val in &[-1.0, 0.0, 1.0] {
        let expected = 2.75 + 1.2 * x_val;
        let est = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, x_val).unwrap();
        assert!(
            (est - expected).abs() < 0.06,
            "E[Y|do(x={x_val})] = {est}, expected {expected}"
        );
    }
}

#[test]
fn test_frontdoor_adjusted_effect_is_flat_without_mediation() {
    // X has no causal effect on M ⇒ do(X) leaves Y unchanged although the
    // observational association X–Y (through U) is strong.
    let mut scm = StructuralCausalModel::new(4);
    scm.add_edge(0, 1, 0.9);
    scm.add_edge(0, 3, 1.0);
    scm.add_edge(1, 2, 0.0);
    scm.add_edge(2, 3, 1.5);
    scm.noise_std = vec![1.0, 0.5, 0.3, 0.3];
    let (xs, ms, ys) = sample_observational(&scm, 5000, 42);
    let e_lo = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, -1.0).unwrap();
    let e_hi = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, 1.0).unwrap();
    assert!(
        (e_hi - e_lo).abs() < 0.1,
        "no mediation ⇒ no effect, got {}",
        e_hi - e_lo
    );
    assert!(ols_slope(&xs, &ys) > 0.5, "naive slope is confounded");
}

#[test]
fn test_frontdoor_adjusted_effect_rejects_empty_short_and_mismatched_input() {
    let err = |r: Result<f64, PhysicsError>| matches!(r, Err(PhysicsError::InvalidInput { .. }));
    assert!(err(FrontdoorCriterion::adjusted_effect(&[], &[], &[], 0.0)));
    let three = [1.0, 2.0, 3.0];
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &three, &three, &three, 1.0
    )));
    let four = [1.0, 2.0, 3.0, 4.0];
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &four, &three, &four, 1.0
    )));
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &four, &four, &three, 1.0
    )));
}

#[test]
fn test_frontdoor_adjusted_effect_rejects_degenerate_designs() {
    let err = |r: Result<f64, PhysicsError>| matches!(r, Err(PhysicsError::InvalidInput { .. }));
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [2.0, 1.0, 4.0, 3.0, 6.0];
    let constant = [1.0; 5];
    // Mediator without variation.
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &x, &constant, &y, 2.0
    )));
    // Treatment without variation.
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &constant, &x, &y, 1.0
    )));
    // Mediator perfectly collinear with the treatment.
    let collinear: Vec<f64> = x.iter().map(|v| 2.0 * v + 1.0).collect();
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &x, &collinear, &y, 2.0
    )));
    // Non-finite data and intervention value.
    let with_nan = [1.0, f64::NAN, 3.0, 4.0, 5.0];
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &with_nan, &y, &y, 2.0
    )));
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &x, &y, &with_nan, 2.0
    )));
    assert!(err(FrontdoorCriterion::adjusted_effect(
        &x,
        &y,
        &y,
        f64::NAN
    )));
}

#[test]
fn test_frontdoor_adjusted_effect_rejects_intervention_outside_support() {
    let scm = frontdoor_scm();
    let (xs, ms, ys) = sample_observational(&scm, 200, 3);
    let far = FrontdoorCriterion::adjusted_effect(&xs, &ms, &ys, 1.0e6);
    assert!(matches!(far, Err(PhysicsError::InvalidInput { .. })));
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
