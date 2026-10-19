// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the dimer method: image/curvature conventions, the
//! Henkelman-Jónsson translation force, the modified-Newton rotation and the
//! saddle-point search.
//!
//! Before the fix the translation followed the *full* force for non-negative
//! curvature (sliding back into the minimum), the images were placed at `R₀ + L·N`
//! while the curvature formula assumed `R₀ − L·N` (so that `V = ½ k x²` gave
//! `C = −k`), and `rotation_angle`, `max_rotation_steps` and `rotation_tol` were
//! never used.

use super::*;

/// Lowest eigenvalue of [`quadratic_form`].
const LAMBDA_MIN: f64 = -1.5;

/// Symmetric 6×6 matrix `A = Q diag(λ) Qᵀ` (a two-atom quadratic form
/// `V = ½ xᵀ A x`) with `λ = (−1.5, 0.3, 1, 2, 3.5, 5)` and its lowest
/// eigenvector, `Q` being a product of Givens rotations.
fn quadratic_form() -> (Vec<Vec<f64>>, Vec<f64>) {
    let eigenvalues = [LAMBDA_MIN, 0.3, 1.0, 2.0, 3.5, 5.0];
    let dim = eigenvalues.len();
    let mut q = vec![vec![0.0_f64; dim]; dim];
    for (i, row) in q.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let rotations: [(usize, usize, f64); 8] = [
        (0, 1, 0.7),
        (1, 2, -0.4),
        (2, 3, 1.1),
        (3, 4, 0.3),
        (4, 5, -0.9),
        (0, 5, 0.5),
        (1, 4, 0.8),
        (2, 5, -0.2),
    ];
    for &(i, j, angle) in &rotations {
        let (s, c) = angle.sin_cos();
        for row in q.iter_mut() {
            let (a, b) = (row[i], row[j]);
            row[i] = c * a - s * b;
            row[j] = s * a + c * b;
        }
    }
    let a: Vec<Vec<f64>> = q
        .iter()
        .map(|row_i| {
            q.iter()
                .map(|row_j| {
                    row_i
                        .iter()
                        .zip(row_j.iter())
                        .zip(eigenvalues.iter())
                        .map(|((qi, qj), lambda)| qi * lambda * qj)
                        .sum::<f64>()
                })
                .collect()
        })
        .collect();
    let lowest: Vec<f64> = q.iter().map(|row| row[0]).collect();
    (a, lowest)
}

/// Forces `−A x` of the quadratic form for a two-atom configuration.
fn quadratic_force(a: Vec<Vec<f64>>) -> impl FnMut(&[[f64; 3]]) -> Vec<[f64; 3]> {
    move |x: &[[f64; 3]]| {
        let flat = flatten(x);
        let forces: Vec<f64> = a.iter().map(|row| -dot(row, &flat)).collect();
        forces.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
    }
}

fn two_atoms(values: [f64; 6]) -> Vec<[f64; 3]> {
    vec![
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
    ]
}

fn rayleigh_quotient(a: &[Vec<f64>], axis: &[[f64; 3]]) -> f64 {
    let n = flatten(axis);
    let an: Vec<f64> = a.iter().map(|row| dot(row, &n)).collect();
    dot(&n, &an) / dot(&n, &n)
}

#[test]
fn defaults_are_valid_and_default_matches_new() {
    let tss = TransitionStateSearch::new();
    assert!(tss.dimer_length > 0.0);
    assert!(tss.max_rotation_steps > 0);
    assert!(tss.rotation_angle > 0.0 && tss.rotation_angle < FRAC_PI_2);
    assert!(tss.rotation_tol > 0.0 && tss.translation_step > 0.0);
    let default = TransitionStateSearch::default();
    assert_eq!(default.max_rotation_steps, tss.max_rotation_steps);
    assert_eq!(default.dimer_length, tss.dimer_length);
}

#[test]
fn images_are_placed_symmetrically_about_the_midpoint() {
    let mid = vec![[1.0, 1.0, 1.0], [0.0, -2.0, 0.5]];
    let axis = vec![[1.0, 0.0, 0.0], [0.0, 0.6, 0.8]];
    let l = 0.01;
    let image1 = TransitionStateSearch::build_image1(&mid, &axis, l);
    let image2 = TransitionStateSearch::build_image2(&mid, &axis, l);
    assert!((image1[0][0] - 1.01).abs() < 1e-12);
    assert!((image2[0][0] - 0.99).abs() < 1e-12);
    for atom in 0..2 {
        for d in 0..3 {
            let average = 0.5 * (image1[atom][d] + image2[atom][d]);
            assert!((average - mid[atom][d]).abs() < 1e-14);
            assert!((image1[atom][d] - image2[atom][d] - 2.0 * l * axis[atom][d]).abs() < 1e-14);
        }
    }
}

#[test]
fn curvature_of_a_harmonic_well_is_the_force_constant() {
    // Pre-fix: image 2 sat at R0 + L N, so V = k x²/2 produced C = −k.
    let tss = TransitionStateSearch::new();
    let stiffness = 3.7;
    let mid = vec![[0.4, -0.3, 0.2]];
    let axis = vec![[0.0, 0.6, 0.8]];
    let force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter()
            .map(|r| [-stiffness * r[0], -stiffness * r[1], -stiffness * r[2]])
            .collect()
    };
    let f1 = force(&TransitionStateSearch::build_image1(
        &mid,
        &axis,
        tss.dimer_length,
    ));
    let f2 = force(&TransitionStateSearch::build_image2(
        &mid,
        &axis,
        tss.dimer_length,
    ));
    let c = tss.curvature(&f1, &f2, &axis);
    assert!((c - stiffness).abs() < 1e-9, "C = {c}, k = {stiffness}");
}

#[test]
fn curvature_is_negative_along_the_unstable_direction_of_a_saddle() {
    // V = x² − y²: curvature +2 along x, −2 along y.
    let tss = TransitionStateSearch::new();
    let force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter().map(|r| [-2.0 * r[0], 2.0 * r[1], 0.0]).collect()
    };
    let mid = vec![[0.3, 0.2, 0.0]];
    for (axis, expected) in [([1.0, 0.0, 0.0], 2.0), ([0.0, 1.0, 0.0], -2.0)] {
        let axis = vec![axis];
        let f1 = force(&TransitionStateSearch::build_image1(
            &mid,
            &axis,
            tss.dimer_length,
        ));
        let f2 = force(&TransitionStateSearch::build_image2(
            &mid,
            &axis,
            tss.dimer_length,
        ));
        let c = tss.curvature(&f1, &f2, &axis);
        assert!((c - expected).abs() < 1e-9, "C = {c}, expected {expected}");
    }
}

#[test]
fn translation_follows_minus_the_parallel_force_for_positive_curvature() {
    // F = (1, 2, 0), N = x̂: F∥ = (1, 0, 0), F⊥ = (0, 2, 0).
    // C > 0 → F† = −F∥ (the old code followed the full force F).
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.1;
    let force = vec![[1.0, 2.0, 0.0]];
    let axis = vec![[1.0, 0.0, 0.0]];
    let mut mid = vec![[0.0, 0.0, 0.0]];
    tss.translate_midpoint(&mut mid, &force, &axis, 3.0);
    assert!((mid[0][0] + 0.1).abs() < 1e-15, "x = {}", mid[0][0]);
    assert!(
        mid[0][1].abs() < 1e-15,
        "y must not move, got {}",
        mid[0][1]
    );
    // Zero curvature belongs to the same branch.
    let mut flat = vec![[0.0, 0.0, 0.0]];
    tss.translate_midpoint(&mut flat, &force, &axis, 0.0);
    assert_eq!(flat, mid);
}

#[test]
fn translation_inverts_only_the_parallel_force_for_negative_curvature() {
    // C < 0 → F† = F − 2 F∥ = (−1, 2, 0).
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.1;
    let force = vec![[1.0, 2.0, 0.0]];
    let axis = vec![[1.0, 0.0, 0.0]];
    let mut mid = vec![[0.5, 0.5, 0.5]];
    tss.translate_midpoint(&mut mid, &force, &axis, -3.0);
    assert!((mid[0][0] - 0.4).abs() < 1e-15);
    assert!((mid[0][1] - 0.7).abs() < 1e-15);
    assert!((mid[0][2] - 0.5).abs() < 1e-15);
}

#[test]
fn translation_projects_onto_the_axis_direction_only() {
    // The same result for an axis of norm 3 and for a skew axis over two atoms.
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.2;
    let force = vec![[1.0, 2.0, 0.0], [0.5, -1.0, 2.0]];
    let unit = vec![[0.0, 0.6, 0.0], [0.0, 0.0, 0.8]];
    let scaled = vec![[0.0, 1.8, 0.0], [0.0, 0.0, 2.4]];
    for curvature in [-1.0, 1.0] {
        let mut a = vec![[0.0; 3]; 2];
        let mut b = vec![[0.0; 3]; 2];
        tss.translate_midpoint(&mut a, &force, &unit, curvature);
        tss.translate_midpoint(&mut b, &force, &scaled, curvature);
        for atom in 0..2 {
            for d in 0..3 {
                assert!((a[atom][d] - b[atom][d]).abs() < 1e-14);
            }
        }
    }
    // F·N = 0.6·2 + 0.8·2 = 2.8, so for C > 0 the displacement is −0.2 · 2.8 · N.
    let mut mid = vec![[0.0; 3]; 2];
    tss.translate_midpoint(&mut mid, &force, &unit, 1.0);
    assert!((mid[0][1] + 0.2 * 2.8 * 0.6).abs() < 1e-14);
    assert!((mid[1][2] + 0.2 * 2.8 * 0.8).abs() < 1e-14);
    assert!(mid[0][0].abs() < 1e-14 && mid[1][0].abs() < 1e-14);
}

#[test]
fn rotation_aligns_the_dimer_with_the_lowest_mode_of_a_quadratic_form() {
    let (a, v_min) = quadratic_form();
    let mut tss = TransitionStateSearch::new();
    tss.max_rotation_steps = 200;
    tss.rotation_tol = 1e-9;
    let mid = two_atoms([0.3, -0.2, 0.5, 0.1, -0.4, 0.2]);
    let mut force = quadratic_force(a.clone());
    let f0 = force(&mid);
    let mut axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
    // The start axis is nearly orthogonal to the lowest mode and far above its curvature.
    let initial_curvature = rayleigh_quotient(&a, &axis);
    let start = flatten(&axis);
    let start_overlap = dot(&start, &v_min).abs() / dot(&start, &start).sqrt();
    assert!(
        initial_curvature > 0.5,
        "start curvature {initial_curvature}"
    );
    assert!(start_overlap < 0.01, "start overlap {start_overlap}");
    let rotation = tss
        .rotate(&mid, &mut axis, &f0, &mut force)
        .expect("rotation inputs are valid");
    assert!(
        rotation.converged,
        "rotation did not converge: {rotation:?}"
    );
    assert!(
        (rotation.curvature - LAMBDA_MIN).abs() < 1e-7,
        "C = {}, lowest eigenvalue {LAMBDA_MIN}",
        rotation.curvature
    );
    let n = flatten(&axis);
    assert!(
        (dot(&n, &n) - 1.0).abs() < 1e-12,
        "axis must stay normalised"
    );
    assert!(
        dot(&n, &v_min).abs() > 1.0 - 1e-9,
        "|N·v_min| = {}",
        dot(&n, &v_min).abs()
    );
    assert!(rotation.steps > 1 && rotation.force_evaluations >= 2 * rotation.steps);
}

#[test]
fn rotation_curvature_decreases_with_the_iteration_budget() {
    let (a, _) = quadratic_form();
    let mid = two_atoms([0.3, -0.2, 0.5, 0.1, -0.4, 0.2]);
    let start_axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
    let mut previous = rayleigh_quotient(&a, &start_axis);
    for budget in 1..=8 {
        let mut tss = TransitionStateSearch::new();
        tss.max_rotation_steps = budget;
        tss.rotation_tol = 1e-12;
        let mut force = quadratic_force(a.clone());
        let f0 = force(&mid);
        let mut axis = start_axis.clone();
        let rotation = tss
            .rotate(&mid, &mut axis, &f0, &mut force)
            .expect("rotation inputs are valid");
        assert_eq!(rotation.steps, budget);
        assert!(!rotation.converged);
        assert!(
            rotation.curvature < previous - 1e-9,
            "budget {budget}: C = {} did not drop below {previous}",
            rotation.curvature
        );
        // The reported curvature belongs to the rotated axis.
        let rayleigh = rayleigh_quotient(&a, &axis);
        assert!(
            (rotation.curvature - rayleigh).abs() < 1e-6,
            "{} vs {rayleigh}",
            rotation.curvature
        );
        previous = rotation.curvature;
    }
}

#[test]
fn rotation_stops_at_the_iteration_budget_and_counts_force_calls() {
    let (a, _) = quadratic_form();
    let mut tss = TransitionStateSearch::new();
    tss.max_rotation_steps = 2;
    tss.rotation_tol = 1e-12;
    let mid = two_atoms([0.3, -0.2, 0.5, 0.1, -0.4, 0.2]);
    let mut calls = 0_usize;
    let mut quadratic = quadratic_force(a);
    let f0 = quadratic(&mid);
    let mut axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
    let rotation = tss
        .rotate(&mid, &mut axis, &f0, |x: &[[f64; 3]]| {
            calls += 1;
            quadratic(x)
        })
        .expect("rotation inputs are valid");
    assert_eq!(rotation.steps, 2);
    assert!(!rotation.converged);
    // One evaluation up front, then a trial and a post-rotation evaluation per step.
    assert_eq!(rotation.force_evaluations, 1 + 2 * 2);
    assert_eq!(calls, rotation.force_evaluations);

    // A zero budget only measures the curvature along the given axis.
    tss.max_rotation_steps = 0;
    let mut untouched = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
    let norm = dot(&flatten(&untouched), &flatten(&untouched)).sqrt();
    for atom in untouched.iter_mut() {
        for c in atom.iter_mut() {
            *c /= norm;
        }
    }
    let (a, _) = quadratic_form();
    let mut force = quadratic_force(a.clone());
    let none = tss
        .rotate(&mid, &mut untouched.clone(), &f0, &mut force)
        .expect("rotation inputs are valid");
    assert_eq!(none.steps, 0);
    assert!((none.curvature - rayleigh_quotient(&a, &untouched)).abs() < 1e-6);
}

#[test]
fn loose_rotation_tolerance_stops_earlier_than_a_tight_one() {
    let (a, _) = quadratic_form();
    let mid = two_atoms([0.3, -0.2, 0.5, 0.1, -0.4, 0.2]);
    let steps_for = |tolerance: f64| {
        let mut tss = TransitionStateSearch::new();
        tss.max_rotation_steps = 500;
        tss.rotation_tol = tolerance;
        let mut force = quadratic_force(a.clone());
        let f0 = force(&mid);
        let mut axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
        let rotation = tss
            .rotate(&mid, &mut axis, &f0, &mut force)
            .expect("rotation inputs are valid");
        assert!(rotation.converged);
        rotation.steps
    };
    let loose = steps_for(0.2);
    let tight = steps_for(1e-9);
    assert!(loose < tight, "loose {loose} vs tight {tight}");
}

#[test]
fn rotation_is_immediate_for_an_isotropic_potential() {
    // Every direction has the same curvature: the rotational force vanishes.
    let stiffness = 2.3;
    let tss = TransitionStateSearch::new();
    let mut force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter()
            .map(|r| [-stiffness * r[0], -stiffness * r[1], -stiffness * r[2]])
            .collect()
    };
    let mid = vec![[0.3, -0.2, 0.5]];
    let f0 = force(&mid);
    let mut axis = vec![[1.0, 2.0, 2.0]];
    let rotation = tss
        .rotate(&mid, &mut axis, &f0, &mut force)
        .expect("rotation inputs are valid");
    assert!(rotation.converged && rotation.steps == 0);
    assert!((rotation.curvature - stiffness).abs() < 1e-9);
    // Only normalised, not rotated.
    assert!((axis[0][0] - 1.0 / 3.0).abs() < 1e-14 && (axis[0][2] - 2.0 / 3.0).abs() < 1e-14);
}

#[test]
fn rotation_validates_its_input() {
    let tss = TransitionStateSearch::new();
    let zero_force = |x: &[[f64; 3]]| vec![[0.0; 3]; x.len()];
    let mid = vec![[0.0; 3], [1.0, 0.0, 0.0]];
    let f0 = vec![[0.0; 3]; 2];

    let mut short_axis = vec![[1.0, 0.0, 0.0]];
    assert!(tss.rotate(&mid, &mut short_axis, &f0, zero_force).is_err());
    let mut zero_axis = vec![[0.0; 3]; 2];
    assert!(tss.rotate(&mid, &mut zero_axis, &f0, zero_force).is_err());
    let mut nan_axis = vec![[f64::NAN, 0.0, 0.0], [0.0; 3]];
    assert!(tss.rotate(&mid, &mut nan_axis, &f0, zero_force).is_err());
    let mut axis = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(tss.rotate(&mid, &mut axis, &f0[..1], zero_force).is_err());
    // The force routine must return one force per atom.
    assert!(
        tss.rotate(&mid, &mut axis, &f0, |_: &[[f64; 3]]| vec![[0.0; 3]])
            .is_err()
    );
    // Non-finite forces are reported rather than propagated.
    assert!(
        tss.rotate(&mid, &mut axis, &f0, |x: &[[f64; 3]]| vec![
            [f64::NAN; 3];
            x.len()
        ])
        .is_err()
    );

    let mut bad_length = TransitionStateSearch::new();
    bad_length.dimer_length = 0.0;
    assert!(bad_length.rotate(&mid, &mut axis, &f0, zero_force).is_err());
    let mut bad_angle = TransitionStateSearch::new();
    bad_angle.rotation_angle = FRAC_PI_2;
    assert!(bad_angle.rotate(&mid, &mut axis, &f0, zero_force).is_err());
    bad_angle.rotation_angle = 0.0;
    assert!(bad_angle.rotate(&mid, &mut axis, &f0, zero_force).is_err());
}

#[test]
fn trial_rotation_angle_does_not_change_the_converged_mode_of_a_quadratic_form() {
    // For a quadratic potential the fitted C(θ) is exact, so any trial angle finds the
    // same minimum; the parameter only enters the finite-difference probe.
    let (a, v_min) = quadratic_form();
    let mid = two_atoms([0.3, -0.2, 0.5, 0.1, -0.4, 0.2]);
    for &angle in &[0.005, 0.01, 0.1, 0.5, 1.2] {
        let mut tss = TransitionStateSearch::new();
        tss.rotation_angle = angle;
        tss.max_rotation_steps = 300;
        tss.rotation_tol = 1e-9;
        let mut force = quadratic_force(a.clone());
        let f0 = force(&mid);
        let mut axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
        let rotation = tss
            .rotate(&mid, &mut axis, &f0, &mut force)
            .expect("rotation inputs are valid");
        assert!(rotation.converged, "trial angle {angle}: {rotation:?}");
        assert!(
            (rotation.curvature - LAMBDA_MIN).abs() < 1e-6,
            "trial angle {angle}"
        );
        assert!(dot(&flatten(&axis), &v_min).abs() > 1.0 - 1e-8);
    }
}

#[test]
fn saddle_of_a_quadratic_form_is_found_from_a_generic_start() {
    // One negative eigenvalue → a single saddle at the origin, reached from anywhere.
    let (a, v_min) = quadratic_form();
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.3;
    let start = two_atoms([1.0, -0.6, 0.8, 0.4, -1.1, 0.7]);
    let axis = two_atoms([0.2, 0.5, -0.3, 0.8, 0.1, -0.4]);
    let result = tss
        .find_saddle(&start, &axis, 1e-7, 1000, quadratic_force(a))
        .expect("search inputs are valid");
    assert!(result.converged, "search did not converge: {result:?}");
    assert!(result.max_force <= 1e-7);
    for atom in &result.midpoint {
        for c in atom {
            assert!(c.abs() < 1e-6, "midpoint {:?}", result.midpoint);
        }
    }
    assert!(
        (result.curvature - LAMBDA_MIN).abs() < 1e-3,
        "C = {}",
        result.curvature
    );
    assert!(dot(&flatten(&result.axis), &v_min).abs() > 1.0 - 1e-6);
    assert!(result.force_evaluations > result.iterations);
}

#[test]
fn saddle_of_x_squared_minus_y_squared_is_the_origin() {
    // V = x² − y² + z² for one atom; the lowest mode is y (C = −2).
    let tss = TransitionStateSearch::new();
    let force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter()
            .map(|r| [-2.0 * r[0], 2.0 * r[1], -2.0 * r[2]])
            .collect()
    };
    let result = tss
        .find_saddle(&[[1.0, 0.5, -0.7]], &[[0.6, 0.8, 0.0]], 1e-7, 400, force)
        .expect("search inputs are valid");
    assert!(result.converged, "{result:?}");
    for c in result.midpoint[0] {
        assert!(c.abs() < 1e-6, "midpoint {:?}", result.midpoint);
    }
    assert!((result.curvature + 2.0).abs() < 1e-6);
    assert!(
        result.axis[0][1].abs() > 1.0 - 1e-9,
        "axis {:?}",
        result.axis
    );
}

#[test]
fn search_climbs_out_of_a_convex_region_along_the_lowest_mode() {
    // V = x² + z² + (y² − 1)²/4: wells at y = ±1, saddle at the origin.  Start next to the
    // y = +1 minimum where every curvature (1.43 along y, 2 along x and z) is positive.
    // With the full force the dimer would slide into the minimum; −F∥ climbs along y.
    let force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter()
            .map(|r| [-2.0 * r[0], -r[1] * (r[1] * r[1] - 1.0), -2.0 * r[2]])
            .collect()
    };
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.05;
    let start = vec![[0.3, 0.9, 0.1]];
    let axis = vec![[0.3, 0.9, 0.2]];
    let mut probe_axis = axis.clone();
    let f_start = force(&start);
    let initial = tss
        .rotate(&start, &mut probe_axis, &f_start, force)
        .expect("rotation inputs are valid");
    assert!(
        initial.curvature > 1.0,
        "start must be convex, C = {}",
        initial.curvature
    );

    let result = tss
        .find_saddle(&start, &axis, 1e-6, 2000, force)
        .expect("search inputs are valid");
    assert!(result.converged, "{result:?}");
    for c in result.midpoint[0] {
        assert!(c.abs() < 1e-5, "midpoint {:?}", result.midpoint);
    }
    // C = V'' of the quartic + O(L²) finite-difference bias: −1 + 3 L².
    assert!(
        (result.curvature + 1.0).abs() < 1e-3,
        "C = {}",
        result.curvature
    );
    assert!(result.axis[0][1].abs() > 1.0 - 1e-6);
}

#[test]
fn search_reports_missing_convergence_and_bad_input() {
    let tss = TransitionStateSearch::new();
    let force = |x: &[[f64; 3]]| -> Vec<[f64; 3]> {
        x.iter()
            .map(|r| [-2.0 * r[0], 2.0 * r[1], -2.0 * r[2]])
            .collect()
    };
    let result = tss
        .find_saddle(&[[1.0, 0.5, -0.7]], &[[0.6, 0.8, 0.0]], 1e-12, 3, force)
        .expect("search inputs are valid");
    assert!(!result.converged);
    assert_eq!(result.iterations, 3);
    assert!(result.max_force > 1e-12);

    assert!(
        tss.find_saddle(&[[0.0; 3]], &[[1.0, 0.0, 0.0], [0.0; 3]], 1e-6, 10, force)
            .is_err()
    );
    assert!(
        tss.find_saddle(
            &[[0.0; 3]],
            &[[1.0, 0.0, 0.0]],
            1e-6,
            10,
            |_: &[[f64; 3]]| { vec![[0.0; 3]; 2] }
        )
        .is_err()
    );
}

/// Deterministic linear congruential generator for reproducible random test problems.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn uniform(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.next_unit()
    }

    fn atoms(&mut self, count: usize) -> Vec<[f64; 3]> {
        (0..count)
            .map(|_| {
                [
                    self.uniform(-1.0, 1.0),
                    self.uniform(-1.0, 1.0),
                    self.uniform(-1.0, 1.0),
                ]
            })
            .collect()
    }
}

/// Random quadratic form `A = Q diag(λ) Qᵀ` in `dim` dimensions with the single negative
/// eigenvalue `λ₀ ∈ [−3, −0.5]` (eigenvector `Q e₀`) and all others in `[other_low, 6]`;
/// `Q` is a product of random Householder reflections.
fn random_quadratic_form(
    rng: &mut Lcg,
    dim: usize,
    other_low: f64,
) -> (Vec<Vec<f64>>, Vec<f64>, f64) {
    let mut eigenvalues: Vec<f64> = (0..dim).map(|_| rng.uniform(other_low, 6.0)).collect();
    let lowest = rng.uniform(-3.0, -0.5);
    eigenvalues[0] = lowest;
    let mut q: Vec<Vec<f64>> = (0..dim)
        .map(|i| (0..dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _ in 0..dim + 3 {
        let mut v: Vec<f64> = (0..dim).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let norm = dot(&v, &v).sqrt();
        for component in v.iter_mut() {
            *component /= norm;
        }
        for row in q.iter_mut() {
            let projection = dot(row, &v);
            for (x, vi) in row.iter_mut().zip(v.iter()) {
                *x -= 2.0 * projection * vi;
            }
        }
    }
    let a: Vec<Vec<f64>> = q
        .iter()
        .map(|row_i| {
            q.iter()
                .map(|row_j| {
                    row_i
                        .iter()
                        .zip(row_j.iter())
                        .zip(eigenvalues.iter())
                        .map(|((qi, qj), lambda)| qi * lambda * qj)
                        .sum::<f64>()
                })
                .collect()
        })
        .collect();
    let lowest_mode: Vec<f64> = q.iter().map(|row| row[0]).collect();
    (a, lowest_mode, lowest)
}

#[test]
fn rotation_converges_to_the_lowest_mode_of_random_quadratic_forms() {
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    let mut tss = TransitionStateSearch::new();
    tss.max_rotation_steps = 3000;
    tss.rotation_tol = 1e-9;
    for case in 0..120 {
        let atoms = 2 + case % 6;
        let (a, v_min, lambda_min) = random_quadratic_form(&mut rng, 3 * atoms, 0.0);
        let mid = rng.atoms(atoms);
        let mut axis = rng.atoms(atoms);
        let mut force = quadratic_force(a);
        let f0 = force(&mid);
        let rotation = tss
            .rotate(&mid, &mut axis, &f0, &mut force)
            .expect("rotation inputs are valid");
        let overlap = dot(&flatten(&axis), &v_min).abs();
        assert!(
            rotation.converged,
            "case {case}: no convergence after {} steps (overlap {overlap})",
            rotation.steps
        );
        assert!(
            (rotation.curvature - lambda_min).abs() < 1e-6,
            "case {case}: C = {} vs lowest eigenvalue {lambda_min}",
            rotation.curvature
        );
        assert!(overlap > 1.0 - 1e-7, "case {case}: |N·v_min| = {overlap}");
    }
}

#[test]
fn saddle_search_converges_on_random_quadratic_forms() {
    // One negative eigenvalue (the lowest mode) and positive curvature elsewhere: the unique
    // saddle is the origin, reached from random starts and random initial axes.
    let mut rng = Lcg(0xfeed_beef_cafe_1234);
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 0.25;
    for case in 0..60 {
        let atoms = 2 + case % 5;
        let (a, _, lambda_min) = random_quadratic_form(&mut rng, 3 * atoms, 0.5);
        let start = rng.atoms(atoms);
        let axis = rng.atoms(atoms);
        let result = tss
            .find_saddle(&start, &axis, 1e-7, 5000, quadratic_force(a))
            .expect("search inputs are valid");
        assert!(
            result.converged,
            "case {case}: {} iterations, max force {}, C {} (lowest {lambda_min})",
            result.iterations, result.max_force, result.curvature
        );
        let distance = result
            .midpoint
            .iter()
            .flat_map(|atom| atom.iter())
            .fold(0.0_f64, |acc, c| acc.max(c.abs()));
        assert!(
            distance < 1e-5,
            "case {case}: midpoint {distance} from the saddle"
        );
        assert!(
            (result.curvature - lambda_min).abs() < 1e-3,
            "case {case}: C = {} vs {lambda_min}",
            result.curvature
        );
    }
}

// Müller-Brown potential, the standard 2-D benchmark for saddle-point searches:
// V = Σ_i A_i exp[a_i (x − x_i)² + b_i (x − x_i)(y − y_i) + c_i (y − y_i)²].
const MB_A: [f64; 4] = [-200.0, -100.0, -170.0, 15.0];
const MB_AA: [f64; 4] = [-1.0, -1.0, -6.5, 0.7];
const MB_BB: [f64; 4] = [0.0, 0.0, 11.0, 0.6];
const MB_CC: [f64; 4] = [-10.0, -10.0, -6.5, 0.7];
const MB_X0: [f64; 4] = [1.0, 0.0, -0.5, -1.0];
const MB_Y0: [f64; 4] = [0.0, 0.5, 1.5, 1.0];
/// Stiffness of the harmonic confinement along z, large enough that the out-of-plane
/// mode never becomes the lowest-curvature mode of the dimer.
const MB_KZ: f64 = 1000.0;

fn mueller_brown_energy(x: f64, y: f64) -> f64 {
    (0..4)
        .map(|i| {
            let (dx, dy) = (x - MB_X0[i], y - MB_Y0[i]);
            MB_A[i] * (MB_AA[i] * dx * dx + MB_BB[i] * dx * dy + MB_CC[i] * dy * dy).exp()
        })
        .sum()
}

/// Forces `−∇V` of the Müller-Brown surface (plus the z confinement) for one atom.
fn mueller_brown_force(positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|r| {
            let (mut fx, mut fy) = (0.0, 0.0);
            for i in 0..4 {
                let (dx, dy) = (r[0] - MB_X0[i], r[1] - MB_Y0[i]);
                let term =
                    MB_A[i] * (MB_AA[i] * dx * dx + MB_BB[i] * dx * dy + MB_CC[i] * dy * dy).exp();
                fx -= term * (2.0 * MB_AA[i] * dx + MB_BB[i] * dy);
                fy -= term * (MB_BB[i] * dx + 2.0 * MB_CC[i] * dy);
            }
            [fx, fy, -MB_KZ * r[2]]
        })
        .collect()
}

#[test]
fn dimer_finds_both_saddle_points_of_the_mueller_brown_surface() {
    // Literature saddle points: (−0.822, 0.624) with V = −40.665 and
    // (0.212, 0.293) with V = −72.249; started next to the three minima and in between.
    let mut tss = TransitionStateSearch::new();
    tss.translation_step = 2.0e-4;
    tss.max_rotation_steps = 20;
    let axis = [[0.5, 0.5, 0.1]];
    let cases = [
        ([0.6, 0.05], (0.212_485, 0.292_988, -72.248_9)),
        ([-0.5, 1.4], (-0.822_000, 0.624_315, -40.664_9)),
        ([-0.4, 0.9], (-0.822_000, 0.624_315, -40.664_9)),
        ([0.3, 0.2], (0.212_485, 0.292_988, -72.248_9)),
    ];
    for (start, (x_ts, y_ts, v_ts)) in cases {
        let result = tss
            .find_saddle(
                &[[start[0], start[1], 0.0]],
                &axis,
                1e-3,
                100_000,
                mueller_brown_force,
            )
            .expect("search inputs are valid");
        let [x, y, z] = result.midpoint[0];
        assert!(result.converged, "start {start:?}: {result:?}");
        assert!(
            (x - x_ts).abs() < 2e-3 && (y - y_ts).abs() < 2e-3 && z.abs() < 1e-6,
            "start {start:?}: converged to ({x}, {y}, {z}), expected ({x_ts}, {y_ts})"
        );
        let energy = mueller_brown_energy(x, y);
        assert!(
            (energy - v_ts).abs() < 0.05,
            "start {start:?}: V = {energy}, expected {v_ts}"
        );
        assert!(result.curvature < -100.0, "C = {}", result.curvature);
        assert!(
            result.axis[0][2].abs() < 1e-3,
            "axis leaves the plane: {:?}",
            result.axis
        );

        // A first-order saddle: the in-plane Hessian (finite differences of the force) has
        // exactly one negative eigenvalue.
        let h = 1e-5;
        let fx = |px: f64, py: f64| mueller_brown_force(&[[px, py, 0.0]])[0];
        let (right, left) = (fx(x + h, y), fx(x - h, y));
        let (up, down) = (fx(x, y + h), fx(x, y - h));
        let (h_xx, h_xy) = (
            -(right[0] - left[0]) / (2.0 * h),
            -(right[1] - left[1]) / (2.0 * h),
        );
        let h_yy = -(up[1] - down[1]) / (2.0 * h);
        assert!(
            h_xx * h_yy - h_xy * h_xy < 0.0,
            "Hessian determinant is not negative"
        );
    }
}
