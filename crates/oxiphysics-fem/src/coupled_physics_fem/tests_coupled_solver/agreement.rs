// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cross-scheme agreement of the coupled solver for one-way coupling.
//!
//! With `C_su = 0` the one-pass staggered scheme is exact, so all three schemes must return the
//! same fields even though `C_us` is nonzero.

use super::*;

/// Largest absolute entry-wise difference between two solutions `(u, s)`.
fn max_difference(a: &CoupledResult, b: &CoupledResult) -> f64 {
    a.displacement
        .iter()
        .chain(a.scalar_field.iter())
        .zip(b.displacement.iter().chain(b.scalar_field.iter()))
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

#[test]
fn test_all_three_schemes_agree_for_one_way_coupling_with_nonzero_c() {
    let (n_u, n_s) = (4, 3);
    let criteria = ConvergenceCriteria::new(1e-13, 1e-13, 200);
    for seed in [3_u64, 17, 29] {
        let problem = Problem::random(n_u, n_s, 0.8, true, seed);
        let reference = problem.direct_solution();
        let results: Vec<(SolverScheme, CoupledResult)> = all_schemes()
            .into_iter()
            .map(|scheme| (scheme, problem.solve(scheme, criteria.clone())))
            .collect();
        for (scheme, r) in &results {
            assert!(
                r.converged,
                "{scheme:?}, seed {seed}: residuals ({:e}, {:e})",
                r.res_disp, r.res_scalar
            );
            for (i, &u) in r.displacement.iter().enumerate() {
                assert!(
                    (u - reference[i]).abs() < 1e-8,
                    "{scheme:?}, seed {seed}: u[{i}] = {u} vs direct {}",
                    reference[i]
                );
            }
            for (i, &s) in r.scalar_field.iter().enumerate() {
                assert!(
                    (s - reference[n_u + i]).abs() < 1e-8,
                    "{scheme:?}, seed {seed}: s[{i}] = {s} vs direct {}",
                    reference[n_u + i]
                );
            }
        }
        let (_, first) = &results[0];
        for (scheme, r) in results.iter().skip(1) {
            let diff = max_difference(first, r);
            assert!(
                diff < 1e-8,
                "{scheme:?} differs from the monolithic solution by {diff:e} (seed {seed})"
            );
        }
    }
}

#[test]
fn test_one_way_coupling_changes_the_displacement() {
    // Guards the test above: with C_us = 0 the displacement would be different, i.e. the coupling
    // term (and its sign) really enters the solution.
    let (n_u, n_s) = (4, 3);
    let coupled = Problem::random(n_u, n_s, 0.8, true, 3);
    let decoupled = Problem {
        k_u: coupled.k_u.clone(),
        k_s: coupled.k_s.clone(),
        c_us: constant_block(n_u, n_s, 0.0),
        c_su: constant_block(n_s, n_u, 0.0),
        f_u: coupled.f_u.clone(),
        f_s: coupled.f_s.clone(),
    };
    let criteria = ConvergenceCriteria::new(1e-13, 1e-13, 200);
    let with = coupled.solve(SolverScheme::Monolithic, criteria.clone());
    let without = decoupled.solve(SolverScheme::Monolithic, criteria);
    let shift = with
        .displacement
        .iter()
        .zip(without.displacement.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        shift > 1e-3,
        "coupling shifted the displacement by {shift:e}"
    );
    // K_u u + C_us s = f_u  =>  the shift equals -K_u^{-1} C_us s; check one row of the equation.
    for i in 0..n_u {
        let lhs = row_dot(&coupled.k_u[i], &with.displacement)
            + row_dot(&coupled.c_us[i], &with.scalar_field);
        assert!((lhs - coupled.f_u[i]).abs() < 1e-9, "row {i}: {lhs}");
    }
}
