// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the eigenvalue-based stability analysis and the
//! tangent-based bifurcation classification of the `continuation` module.
//!
//! Split out of `continuation.rs` to keep that file under the 2000-line limit.
//! Included via `#[path]` as `continuation::bifurcation_tests`.

use super::*;

fn rows(m: &[&[f64]]) -> Vec<Vec<f64>> {
    m.iter().map(|r| r.to_vec()).collect()
}

/// Fold problem F(u, λ) = u² − λ (branch λ = u², fold at the origin).
fn parabola(u: &[f64], lam: f64) -> Vec<f64> {
    vec![u[0] * u[0] - lam]
}

fn parabola_jac(u: &[f64], _lam: f64) -> Vec<Vec<f64>> {
    vec![vec![2.0 * u[0]]]
}

/// Pitchfork problem F(u, λ) = (u₀³ − λ u₀, u₁ + u₀): trivial branch u = 0.
fn pitchfork(u: &[f64], lam: f64) -> Vec<f64> {
    vec![u[0].powi(3) - lam * u[0], u[1] + u[0]]
}

fn pitchfork_jac(u: &[f64], lam: f64) -> Vec<Vec<f64>> {
    vec![vec![3.0 * u[0] * u[0] - lam, 0.0], vec![1.0, 1.0]]
}

/// The start state followed by every accepted state of a path-following run.
fn walk(
    start: ContinuationState,
    stepper: PseudoArcLengthContinuation,
    max_steps: usize,
    f: &ResidualFn,
    jac: &JacobianFn,
) -> (Vec<ContinuationState>, Vec<PathStep>) {
    let steps = PathFollowing::new(max_steps, stepper).follow(start.clone(), f, jac);
    let mut states = vec![start];
    states.extend(steps.iter().map(|step| step.state.clone()));
    (states, steps)
}

/// Indices `k` such that det J changes sign between `states[k]` and `states[k + 1]`.
fn determinant_sign_changes(states: &[ContinuationState], jac: &JacobianFn) -> Vec<usize> {
    states
        .windows(2)
        .enumerate()
        .filter(|(_, pair)| {
            let det_a = matrix_determinant(&jac(&pair[0].u, pair[0].lambda));
            let det_b = matrix_determinant(&jac(&pair[1].u, pair[1].lambda));
            det_a * det_b <= 0.0
        })
        .map(|(k, _)| k)
        .collect()
}

// ---- stability_index ------------------------------------------------------

#[test]
fn stability_index_counts_eigenvalues_not_gershgorin_discs() {
    // Spectrum {11, -4, -4}; the disc test counted no "unstable disc" at all.
    let m = rows(&[&[1.0, 5.0, 5.0], &[5.0, 1.0, 5.0], &[5.0, 5.0, 1.0]]);
    assert_eq!(stability_index(&m), 1);
}

#[test]
fn stability_index_companion_matrix_with_an_imaginary_pair() {
    // Characteristic polynomial (x - 1)(x² + 1): spectrum {1, ±i}.
    let m = rows(&[&[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0], &[1.0, -1.0, 1.0]]);
    assert_eq!(stability_index(&m), 1);
}

#[test]
fn stability_index_triangular_matrix_with_large_off_diagonal_entries() {
    // The spectrum of a triangular matrix is its diagonal {2, -3, 1}, although
    // every Gershgorin disc contains the origin.
    let m = rows(&[&[2.0, 100.0, -50.0], &[0.0, -3.0, 70.0], &[0.0, 0.0, 1.0]]);
    assert_eq!(stability_index(&m), 2);
}

#[test]
fn stability_index_complex_pair_in_the_right_half_plane() {
    // Spectrum {0.1 ± 3i, -1, -2}: the pair counts twice.
    let m = rows(&[
        &[0.1, -3.0, 0.0, 0.0],
        &[3.0, 0.1, 0.0, 0.0],
        &[0.0, 0.0, -1.0, 0.0],
        &[0.0, 0.0, 0.0, -2.0],
    ]);
    assert_eq!(stability_index(&m), 2);
}

fn matmul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

#[test]
fn stability_index_ignores_rounding_noise_on_the_imaginary_axis() {
    // Spectrum {±3i, -1}: a Hopf point.
    let hopf = rows(&[&[0.0, -3.0, 0.0], &[3.0, 0.0, 0.0], &[0.0, 0.0, -1.0]]);
    assert_eq!(stability_index(&hopf), 0);
    // The same spectrum behind a generic orthogonal similarity: the computed real
    // parts of ±3i are rounding noise of either sign and must not be counted.
    let (c1, s1) = (0.4_f64.cos(), 0.4_f64.sin());
    let (c2, s2) = (0.7_f64.cos(), 0.7_f64.sin());
    let rotation_x = [[1.0, 0.0, 0.0], [0.0, c1, -s1], [0.0, s1, c1]];
    let rotation_y = [[c2, 0.0, s2], [0.0, 1.0, 0.0], [-s2, 0.0, c2]];
    let q = matmul3(&rotation_y, &rotation_x);
    let q_t = [
        [q[0][0], q[1][0], q[2][0]],
        [q[0][1], q[1][1], q[2][1]],
        [q[0][2], q[1][2], q[2][2]],
    ];
    let block = [[0.0, -3.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
    let similar = matmul3(&matmul3(&q, &block), &q_t);
    let m: Vec<Vec<f64>> = similar.iter().map(|r| r.to_vec()).collect();
    assert_eq!(stability_index(&m), 0);
    // Nilpotent matrix with large entries: spectrum {0, 0, 0}.
    let nilpotent = rows(&[&[0.0, 1e3, 7.0], &[0.0, 0.0, 1e3], &[0.0, 0.0, 0.0]]);
    assert_eq!(stability_index(&nilpotent), 0);
}

#[test]
fn stability_index_returns_zero_without_a_spectrum() {
    assert_eq!(stability_index(&[]), 0);
    // Not square, ragged, non-finite: the old code indexed out of bounds or
    // compared NaN.
    assert_eq!(
        stability_index(&rows(&[&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]])),
        0
    );
    assert_eq!(stability_index(&rows(&[&[1.0, 2.0], &[3.0]])), 0);
    assert_eq!(stability_index(&rows(&[&[f64::NAN, 0.0], &[0.0, 1.0]])), 0);
    let infinite = rows(&[
        &[1.0, 0.0, 0.0],
        &[0.0, f64::INFINITY, 0.0],
        &[0.0, 0.0, 1.0],
    ]);
    assert_eq!(stability_index(&infinite), 0);
}

#[test]
fn stability_index_2x2_closed_form_covers_real_and_complex_spectra() {
    assert_eq!(stability_index(&rows(&[&[2.0, 7.0], &[0.0, -1.0]])), 1);
    assert_eq!(stability_index(&rows(&[&[3.0, 1.0], &[-1.0, 3.0]])), 2);
    assert_eq!(stability_index(&rows(&[&[-3.0, 1.0], &[-1.0, -3.0]])), 0);
    // Saddle with a vanishing trace: eigenvalues ±2.
    assert_eq!(stability_index(&rows(&[&[0.0, 4.0], &[1.0, 0.0]])), 1);
    // A zero matrix has no unstable direction.
    assert_eq!(stability_index(&rows(&[&[0.0, 0.0], &[0.0, 0.0]])), 0);
}

// ---- detect_bifurcation ----------------------------------------------------

#[test]
fn detect_bifurcation_sees_an_eigenvalue_crossing_in_a_3x3_jacobian() {
    // J(u) = u I + 5 (11ᵀ − I) has the spectrum {u + 10, u − 5, u − 5}: one
    // eigenvalue crosses zero between u = -11 and u = -9, but no Gershgorin disc
    // is entirely positive for u < 10.
    let jac = |u: &[f64], _lam: f64| -> Vec<Vec<f64>> {
        (0..3)
            .map(|i| (0..3).map(|j| if i == j { u[0] } else { 5.0 }).collect())
            .collect()
    };
    let far = ContinuationState::new(0.0, vec![-20.0], 0.1);
    let below = ContinuationState::new(0.0, vec![-11.0], 0.1);
    let above = ContinuationState::new(0.0, vec![-9.0], 0.1);
    assert!(detect_bifurcation(&below, &above, &jac));
    assert!(!detect_bifurcation(&far, &below, &jac));
}

// ---- StabilityAnalysis::label ---------------------------------------------

#[test]
fn label_uses_eigenvalues_for_matrices_larger_than_2x2() {
    let analysis = StabilityAnalysis::new(1e-6);
    let state = ContinuationState::new(0.0, vec![0.0], 0.1);
    // Spectrum {-1 ± 5i, -2}: asymptotically stable, but the first Gershgorin
    // disc (centre -1, radius 5) contains the origin, so the old test said Marginal.
    let stable =
        |_u: &[f64], _l: f64| rows(&[&[-1.0, 5.0, 0.0], &[-5.0, -1.0, 0.0], &[0.0, 0.0, -2.0]]);
    assert_eq!(analysis.label(&state, &stable), StabilityLabel::Stable);
    // Spectrum {11, -4, -4}: unstable, but no disc lies in the right half-plane.
    let unstable =
        |_u: &[f64], _l: f64| rows(&[&[1.0, 5.0, 5.0], &[5.0, 1.0, 5.0], &[5.0, 5.0, 1.0]]);
    assert_eq!(analysis.label(&state, &unstable), StabilityLabel::Unstable);
    // Spectrum {±5i, -2}: on the imaginary axis.
    let marginal =
        |_u: &[f64], _l: f64| rows(&[&[0.0, -5.0, 0.0], &[5.0, 0.0, 0.0], &[0.0, 0.0, -2.0]]);
    assert_eq!(analysis.label(&state, &marginal), StabilityLabel::Marginal);
}

#[test]
fn label_unstable_eigenvalue_takes_priority_over_a_marginal_one() {
    let analysis = StabilityAnalysis::new(1e-6);
    let state = ContinuationState::new(0.0, vec![0.0], 0.1);
    // Spectrum {0, 5}: the 2×2 path used to report Marginal.
    let two = |_u: &[f64], _l: f64| rows(&[&[0.0, 0.0], &[0.0, 5.0]]);
    assert_eq!(analysis.label(&state, &two), StabilityLabel::Unstable);
    let three =
        |_u: &[f64], _l: f64| rows(&[&[0.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[0.0, 0.0, -1.0]]);
    assert_eq!(analysis.label(&state, &three), StabilityLabel::Unstable);
    // Spectrum {0, -5}: marginal, not unstable.
    let neutral = |_u: &[f64], _l: f64| rows(&[&[0.0, 0.0], &[0.0, -5.0]]);
    assert_eq!(analysis.label(&state, &neutral), StabilityLabel::Marginal);
}

#[test]
fn label_marginal_tolerance_is_an_absolute_band_around_the_axis() {
    let state = ContinuationState::new(0.0, vec![0.0], 0.1);
    let eigenvalue = |_u: &[f64], _l: f64| rows(&[&[0.3]]);
    assert_eq!(
        StabilityAnalysis::new(0.5).label(&state, &eigenvalue),
        StabilityLabel::Marginal
    );
    assert_eq!(
        StabilityAnalysis::new(0.1).label(&state, &eigenvalue),
        StabilityLabel::Unstable
    );
    // A negative or NaN tolerance acts as zero instead of mislabelling everything.
    let tiny = |_u: &[f64], _l: f64| rows(&[&[1e-9]]);
    for tolerance in [-1.0, f64::NAN] {
        let analysis = StabilityAnalysis::new(tolerance);
        assert_eq!(analysis.label(&state, &tiny), StabilityLabel::Unstable);
        let zero = |_u: &[f64], _l: f64| rows(&[&[0.0]]);
        assert_eq!(analysis.label(&state, &zero), StabilityLabel::Marginal);
        let negative = |_u: &[f64], _l: f64| rows(&[&[-1e-9]]);
        assert_eq!(analysis.label(&state, &negative), StabilityLabel::Stable);
    }
}

#[test]
fn label_never_certifies_a_jacobian_without_a_spectrum() {
    let analysis = StabilityAnalysis::new(1e-6);
    let state = ContinuationState::new(0.0, vec![0.0], 0.1);
    let nan =
        |_u: &[f64], _l: f64| rows(&[&[-1.0, f64::NAN, 0.0], &[0.0, -1.0, 0.0], &[0.0, 0.0, -1.0]]);
    assert_eq!(analysis.label(&state, &nan), StabilityLabel::Unstable);
    let ragged = |_u: &[f64], _l: f64| rows(&[&[-1.0, 0.0], &[-1.0]]);
    assert_eq!(analysis.label(&state, &ragged), StabilityLabel::Unstable);
    // An empty Jacobian (no unknowns) has nothing that could be unstable.
    let empty = |_u: &[f64], _l: f64| -> Vec<Vec<f64>> { Vec::new() };
    assert_eq!(analysis.label(&state, &empty), StabilityLabel::Stable);
}

// ---- classify_bifurcation ---------------------------------------------------

#[test]
fn classify_reports_hopf_when_a_complex_pair_crosses_the_imaginary_axis() {
    // J(λ) has the spectrum {λ ± 3i, -1}: det J = -(λ² + 9) < 0 on both sides, so
    // only the stability index changes (0 → 2).
    let jac =
        |_u: &[f64], lam: f64| rows(&[&[lam, -3.0, 0.0], &[3.0, lam, 0.0], &[0.0, 0.0, -1.0]]);
    let before = ContinuationState::new(-0.1, vec![0.0], 0.1);
    let after = ContinuationState::new(0.1, vec![0.0], 0.1);
    let point = classify_bifurcation(&before, &after, &jac).expect("stability index changes");
    assert_eq!(point.bif_type, BifurcationType::Hopf);
    assert!(point.state.lambda.abs() < 1e-12);
    assert!((point.det_j + 9.0).abs() < 1e-9, "det J = {}", point.det_j);
    // No change in the index and no determinant sign change: nothing to report.
    let later = ContinuationState::new(0.3, vec![0.0], 0.1);
    assert!(classify_bifurcation(&after, &later, &jac).is_none());
}

#[test]
fn classify_separates_fold_and_branch_point_by_the_lambda_tangent() {
    let norm = 5.0_f64.sqrt();
    // Two states on either side of the fold of u² = λ with the true tangents
    // (1, 2u)/√(1 + 4u²) oriented along decreasing u: λ decreases before the
    // turning point and increases after it.
    let mut before = ContinuationState::new(0.25, vec![0.5], 0.1);
    before.tangent = vec![-1.0 / 2.0_f64.sqrt(), -1.0 / 2.0_f64.sqrt()];
    let mut after = ContinuationState::new(0.25, vec![-0.5], 0.1);
    after.tangent = vec![-1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt()];
    let fold = classify_bifurcation(&before, &after, &parabola_jac).expect("det J changes sign");
    assert_eq!(fold.bif_type, BifurcationType::Fold);
    // The mean tangent at the turning point has no λ-component and runs along -u.
    assert!(fold.state.tangent[1].abs() < 1e-12);
    assert!((fold.state.tangent[0] + 1.0).abs() < 1e-12);
    assert!(
        fold.det_j.abs() < 1e-12,
        "mid-point is the fold: {}",
        fold.det_j
    );

    // The same Jacobian signs with a tangent that keeps its λ direction is a
    // branch point instead.
    let mut straight_a = ContinuationState::new(0.25, vec![0.5], 0.1);
    straight_a.tangent = vec![-1.0 / norm, 2.0 / norm];
    let mut straight_b = ContinuationState::new(0.25, vec![-0.5], 0.1);
    straight_b.tangent = vec![-1.0 / norm, 2.0 / norm];
    let branch =
        classify_bifurcation(&straight_a, &straight_b, &parabola_jac).expect("det J changes sign");
    assert_eq!(branch.bif_type, BifurcationType::Pitchfork);

    // A vanishing λ-component marks a state sitting at the turning point.
    let mut at_turn = ContinuationState::new(0.0, vec![0.0], 0.1);
    at_turn.tangent = vec![-1.0, 0.0];
    let mut past = ContinuationState::new(0.04, vec![-0.2], 0.1);
    past.tangent = vec![-0.9, 0.4];
    let turning = classify_bifurcation(&at_turn, &past, &parabola_jac).expect("det J = 0");
    assert_eq!(turning.bif_type, BifurcationType::Fold);
}

#[test]
fn classify_reports_unknown_without_a_usable_tangent() {
    let mut empty = ContinuationState::new(0.25, vec![0.5], 0.1);
    empty.tangent.clear();
    let other = ContinuationState::new(0.25, vec![-0.5], 0.1);
    let point = classify_bifurcation(&empty, &other, &parabola_jac).expect("det J changes sign");
    assert_eq!(point.bif_type, BifurcationType::Unknown);

    let mut not_a_number = ContinuationState::new(0.25, vec![0.5], 0.1);
    not_a_number.tangent = vec![0.0, f64::NAN];
    let point =
        classify_bifurcation(&not_a_number, &other, &parabola_jac).expect("det J changes sign");
    assert_eq!(point.bif_type, BifurcationType::Unknown);
    // The midpoint falls back to the placeholder tangent.
    assert_eq!(point.state.tangent, vec![0.0, 1.0]);
}

#[test]
fn placeholder_tangents_are_never_classified_as_a_fold() {
    // ContinuationState::new stores the pure-λ placeholder: whatever the
    // geometry, no direction reversal can be read from it.
    let a = ContinuationState::new(0.25, vec![0.5], 0.1);
    let b = ContinuationState::new(0.25, vec![-0.5], 0.1);
    let point = classify_bifurcation(&a, &b, &parabola_jac).expect("det J changes sign");
    assert_eq!(point.bif_type, BifurcationType::Pitchfork);
}

#[test]
fn path_following_across_the_fold_of_a_parabola_is_classified_as_a_fold() {
    // Walk u² = λ from (1, 1) towards decreasing u, through the fold at the
    // origin (arc length ≈ 1.48 away) and on to the lower branch.
    let norm = 5.0_f64.sqrt();
    let mut start = ContinuationState::new(1.0, vec![1.0], 0.2);
    start.tangent = vec![-1.0 / norm, -2.0 / norm];
    let stepper = PseudoArcLengthContinuation::new(1e-10, 30, 4, 1e-6, 0.3);
    let (states, steps) = walk(start, stepper, 14, &parabola, &parabola_jac);
    assert!(
        states.last().map(|s| s.u[0]).unwrap_or(1.0) < -0.5,
        "the path must have crossed the fold: {:?}",
        states
            .iter()
            .map(|s| (s.u[0], s.lambda))
            .collect::<Vec<_>>()
    );

    // Every state lies on the curve and carries the exact unit tangent ∝ (1, 2u).
    for state in &states[1..] {
        let u = state.u[0];
        assert!(
            (u * u - state.lambda).abs() < 1e-9,
            "off the curve at u = {u}"
        );
        let length: f64 = state.tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((length - 1.0).abs() < 1e-9, "tangent length {length}");
        let parallel = state.tangent[1] - 2.0 * u * state.tangent[0];
        assert!(
            parallel.abs() < 1e-5,
            "tangent not parallel to (1, 2u) at u = {u}"
        );
        // Oriented along decreasing u.
        assert!(state.tangent[0] < 0.0);
    }

    let crossings = determinant_sign_changes(&states, &parabola_jac);
    assert_eq!(
        crossings.len(),
        1,
        "exactly one singular crossing: {crossings:?}"
    );
    let k = crossings[0];
    let (a, b) = (&states[k], &states[k + 1]);
    assert!(a.tangent[1] < 0.0 && b.tangent[1] > 0.0, "λ must reverse");
    let point = classify_bifurcation(a, b, &parabola_jac).expect("fold between the states");
    assert_eq!(point.bif_type, BifurcationType::Fold);
    // The step that crossed the fold is flagged by the path follower.
    assert!(steps[k].fold_detected);
    assert_eq!(steps.iter().filter(|s| s.fold_detected).count(), 1);
}

#[test]
fn path_following_through_the_trivial_pitchfork_branch_is_classified_as_a_branch_point() {
    // Along u = 0 the Jacobian [[−λ, 0], [1, 1]] has det J = −λ: it changes sign
    // at the pitchfork λ = 0 while the branch keeps running towards larger λ.
    let start = ContinuationState::new(-0.35, vec![0.0, 0.0], 0.1);
    let stepper = PseudoArcLengthContinuation::new(1e-10, 30, 4, 1e-6, 0.3);
    let (states, steps) = walk(start, stepper, 8, &pitchfork, &pitchfork_jac);
    assert!(states.last().map(|s| s.lambda).unwrap_or(-1.0) > 0.1);
    for state in &states {
        assert!(
            state.u.iter().all(|x| x.abs() < 1e-9),
            "left the trivial branch"
        );
        assert!(state.tangent[2] > 0.99, "tangent {:?}", state.tangent);
    }
    let crossings = determinant_sign_changes(&states, &pitchfork_jac);
    assert_eq!(
        crossings.len(),
        1,
        "exactly one singular crossing: {crossings:?}"
    );
    let k = crossings[0];
    let point = classify_bifurcation(&states[k], &states[k + 1], &pitchfork_jac)
        .expect("branch point between the states");
    assert_eq!(point.bif_type, BifurcationType::Pitchfork);
    assert!(point.state.lambda.abs() < 0.3);
    // The determinant flag fires (a branch point is a singular Jacobian) and so
    // does the stability index: the unstable eigenvalue −λ appears for λ < 0.
    assert!(steps[k].fold_detected);
    assert!(steps[k].bifurcation_detected);
}

// ---- branch_tangent / step ---------------------------------------------------

#[test]
fn branch_tangent_is_regular_at_a_fold_and_follows_the_reference_orientation() {
    // At the fold (0, 0) of u² = λ the branch runs along ∓u with no λ-component.
    let norm = 5.0_f64.sqrt();
    let reference = [-1.0 / norm, -2.0 / norm];
    let tangent = branch_tangent(&[0.0], 0.0, &reference, &parabola, &parabola_jac)
        .expect("the bordered system is regular at a fold");
    assert!(
        (tangent[0] + 1.0).abs() < 1e-6 && tangent[1].abs() < 1e-6,
        "{tangent:?}"
    );
    let flipped = [1.0 / norm, 2.0 / norm];
    let opposite =
        branch_tangent(&[0.0], 0.0, &flipped, &parabola, &parabola_jac).expect("regular");
    assert!(
        (opposite[0] - 1.0).abs() < 1e-6 && opposite[1].abs() < 1e-6,
        "{opposite:?}"
    );
    // Away from the fold: (1, 2u)/√(1 + 4u²).
    let off = branch_tangent(
        &[1.0],
        1.0,
        &[1.0 / norm, 2.0 / norm],
        &parabola,
        &parabola_jac,
    )
    .expect("regular");
    assert!(
        (off[0] - 1.0 / norm).abs() < 1e-6 && (off[1] - 2.0 / norm).abs() < 1e-6,
        "{off:?}"
    );
}

#[test]
fn branch_tangent_declines_at_a_branch_point_and_for_inconsistent_input() {
    // At the pitchfork the augmented Jacobian [J | F_λ] has a two-dimensional
    // kernel: the tangent is not unique.
    assert!(
        branch_tangent(
            &[0.0, 0.0],
            0.0,
            &[0.0, 0.0, 1.0],
            &pitchfork,
            &pitchfork_jac
        )
        .is_none()
    );
    // Reference of the wrong length / non-square Jacobian.
    assert!(branch_tangent(&[1.0], 1.0, &[1.0], &parabola, &parabola_jac).is_none());
    let wide = |_u: &[f64], _l: f64| rows(&[&[1.0, 2.0]]);
    assert!(branch_tangent(&[1.0], 1.0, &[0.0, 1.0], &parabola, &wide).is_none());
    // A reference that is almost orthogonal to the branch tangent (1, 2)/√5 at
    // (1, 1) means the branch turned by almost 90° within one step: refused.
    let norm = 5.0_f64.sqrt();
    let nearly_orthogonal = [-2.0 / norm, 1.002 / norm];
    assert!(branch_tangent(&[1.0], 1.0, &nearly_orthogonal, &parabola, &parabola_jac).is_none());
    // A moderate turn (about 63° between reference (1, 0) and the tangent) is fine.
    let moderate =
        branch_tangent(&[1.0], 1.0, &[1.0, 0.0], &parabola, &parabola_jac).expect("moderate turn");
    assert!((moderate[0] - 1.0 / norm).abs() < 1e-9 && (moderate[1] - 2.0 / norm).abs() < 1e-9);
}

#[test]
fn step_returns_the_branch_tangent_rather_than_the_chord() {
    // A long step along u² = λ: the chord from the start differs measurably from
    // the tangent at the new point, which must match (1, 2u)/√(1 + 4u²).
    let norm = 5.0_f64.sqrt();
    let mut state = ContinuationState::new(1.0, vec![1.0], 1.0);
    state.tangent = vec![1.0 / norm, 2.0 / norm];
    let start = state.clone();
    let stepper = PseudoArcLengthContinuation::new(1e-12, 50, 5, 1e-6, 1.0);
    let next = stepper
        .step(&mut state, &parabola, &parabola_jac)
        .expect("converges");
    let u = next.u[0];
    let length = (1.0 + 4.0 * u * u).sqrt();
    assert!(
        (next.tangent[0] - 1.0 / length).abs() < 1e-6,
        "{:?}",
        next.tangent
    );
    assert!(
        (next.tangent[1] - 2.0 * u / length).abs() < 1e-6,
        "{:?}",
        next.tangent
    );
    let chord = secant_tangent(&start, &next.u, next.lambda);
    assert!(
        (chord[1] - next.tangent[1]).abs() > 1e-3,
        "chord {chord:?} vs tangent {:?}",
        next.tangent
    );
}
