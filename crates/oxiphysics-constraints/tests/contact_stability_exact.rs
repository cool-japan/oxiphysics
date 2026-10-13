// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Public-API checks of the exact force-closure test and the Ferrari-Canny
//! grasp quality metrics of `oxiphysics_constraints::contact_stability`.

use oxiphysics_constraints::contact_stability::{
    ContactPoint, GraspQuality, Wrench, analyse_contact_stability, force_closure_exact,
    force_closure_test, force_closure_test_full, grasp_quality_epsilon, grasp_quality_q1,
    grasp_quality_q1_detailed, multi_contact_wrench_set,
};

/// The twelve wrenches `+-e_k`.
fn cross_polytope() -> Vec<Wrench> {
    let mut wrenches = Vec::new();
    for axis in 0..6 {
        for sign in [1.0, -1.0] {
            let mut w = [0.0; 6];
            w[axis] = sign;
            wrenches.push(w);
        }
    }
    wrenches
}

/// Six contacts on the coordinate axes with inward normals.
fn sphere_grasp(mu: f64) -> Vec<ContactPoint> {
    let axes = [
        ([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        ([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([0.0, 1.0, 0.0], [0.0, -1.0, 0.0]),
        ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]),
        ([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]),
    ];
    axes.iter()
        .map(|&(position, normal)| ContactPoint {
            position,
            normal,
            mu,
        })
        .collect()
}

#[test]
fn per_axis_test_is_only_a_prefilter() {
    // Three wrenches in the xy-plane that lie in an open half-plane, plus
    // +-e_2 .. +-e_5: every axis shows both signs, yet the xy-plane cannot be
    // resisted.
    let mut wrenches: Vec<Wrench> = vec![
        [1.0, 0.1, 0.0, 0.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    ];
    for axis in 2..6 {
        for sign in [1.0, -1.0] {
            let mut w = [0.0; 6];
            w[axis] = sign;
            wrenches.push(w);
        }
    }
    assert!(force_closure_test(&wrenches));
    assert!(!force_closure_exact(&wrenches));
    assert!(!force_closure_test_full(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
}

#[test]
fn cross_polytope_quality_matches_the_closed_form() {
    let wrenches = cross_polytope();
    let expected = 1.0 / 6.0_f64.sqrt();
    assert!(force_closure_exact(&wrenches));
    assert!((grasp_quality_q1(&wrenches) - expected).abs() < 1e-9);
    assert!((grasp_quality_epsilon(&wrenches) - expected).abs() < 1e-9);

    let GraspQuality {
        q1,
        weakest_direction,
        exact,
    } = grasp_quality_q1_detailed(&wrenches);
    assert!(exact);
    assert!((q1 - expected).abs() < 1e-9);
    // The support function in the weakest direction equals Q1.
    let support = wrenches
        .iter()
        .map(|w| {
            w.iter()
                .zip(weakest_direction.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>()
        })
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((support - q1).abs() < 1e-9);
}

#[test]
fn antipodal_contacts_do_not_close_and_score_zero() {
    for mu in [0.0, 0.8] {
        let contacts = [
            ContactPoint {
                position: [1.0, 0.0, 0.0],
                normal: [-1.0, 0.0, 0.0],
                mu,
            },
            ContactPoint {
                position: [-1.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
                mu,
            },
        ];
        let wrenches = multi_contact_wrench_set(&contacts, 8);
        assert!(!force_closure_exact(&wrenches), "mu = {mu}");
        assert_eq!(grasp_quality_q1(&wrenches), 0.0, "mu = {mu}");
        assert_eq!(grasp_quality_epsilon(&wrenches), 0.0, "mu = {mu}");
    }
}

#[test]
fn stability_report_reflects_the_grasp_geometry() {
    let report = analyse_contact_stability(&sphere_grasp(0.8), [0.0; 3], [0.0, 0.0, -9.81], 1.0, 8);
    assert!(report.force_closure);
    assert!(report.q1 > 0.0);
    // The old epsilon was identically 1.0.
    assert!(report.epsilon > 0.0 && report.epsilon < 0.9);

    // Moving one contact off its axis changes the quality.
    let mut moved = sphere_grasp(0.8);
    moved[4].position = [0.4, 0.2, 1.0];
    let moved_report = analyse_contact_stability(&moved, [0.0; 3], [0.0, 0.0, -9.81], 1.0, 8);
    assert!(moved_report.force_closure);
    assert!((moved_report.epsilon - report.epsilon).abs() > 1e-3);
}
