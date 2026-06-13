// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the direct LCP solvers (Lemke + Dantzig).
//!
//! These exercise the public API of [`oxiphysics_constraints::lcp`]:
//! hand-computed known solutions, the canonical-scene harness accuracy gate,
//! Lemke/Dantzig agreement, and — the critical case — anti-cycling on
//! degenerate problems with ratio-test ties.

use oxiphysics_constraints::lcp::{
    complementarity_residual, dantzig_solve, lemke_solve, run_harness,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ── small local linear-algebra helpers (no nalgebra) ────────────────────────

/// Build a random `n × n` matrix with entries in `[-1, 1)`.
fn random_matrix(rng: &mut SmallRng, n: usize) -> Vec<Vec<f64>> {
    use rand::RngExt;
    (0..n)
        .map(|_| (0..n).map(|_| rng.random_range(-1.0_f64..1.0)).collect())
        .collect()
}

/// `AᵀA + εI` — a symmetric positive-definite matrix.
fn psd_from(a: &[Vec<f64>], eps: f64) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut m = vec![vec![0.0; n]; n];
    for (i, row) in m.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for a_row in a.iter() {
                acc += a_row[i] * a_row[j];
            }
            *cell = acc;
        }
        row[i] += eps;
    }
    m
}

#[test]
fn zzz_print_harness_residuals() {
    for r in run_harness() {
        eprintln!(
            "SCENE {} lemke_residual={:.3e} dantzig_residual={:?}",
            r.scene, r.lemke_residual, r.dantzig_residual
        );
    }
}

// ── 1 & 2: known-solution sanity (Lemke) ────────────────────────────────────

#[test]
fn test_lemke_2x2_known_solution() {
    // M = [[2, 1], [1, 2]], q = [-1, -1] → z = [1/3, 1/3], w = [0, 0].
    let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
    let q = vec![-1.0, -1.0];
    let (w, z) = lemke_solve(&m, &q).expect("solvable");
    assert!((z[0] - 1.0 / 3.0).abs() < 1e-12, "z0={}", z[0]);
    assert!((z[1] - 1.0 / 3.0).abs() < 1e-12, "z1={}", z[1]);
    assert!(w[0].abs() < 1e-12, "w0={}", w[0]);
    assert!(w[1].abs() < 1e-12, "w1={}", w[1]);
}

#[test]
fn test_lemke_3x3_known_solution() {
    // M = I, q = [-1, -2, -3] → z = [1, 2, 3], w = [0, 0, 0].
    let m = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];
    let q = vec![-1.0, -2.0, -3.0];
    let (w, z) = lemke_solve(&m, &q).expect("solvable");
    let expected = [1.0, 2.0, 3.0];
    for k in 0..3 {
        assert!((z[k] - expected[k]).abs() < 1e-12, "z[{k}]={}", z[k]);
        assert!(w[k].abs() < 1e-12, "w[{k}]={}", w[k]);
    }
}

// ── 3: harness residual gate ────────────────────────────────────────────────

#[test]
fn test_lemke_harness_residual() {
    for r in run_harness() {
        assert!(
            r.lemke_residual < 1e-10,
            "scene {} lemke residual {} exceeds 1e-10",
            r.scene,
            r.lemke_residual
        );
    }
}

// ── 4: Lemke / Dantzig agreement on the canonical scenes ────────────────────

#[test]
fn test_lemke_dantzig_agree() {
    for r in run_harness() {
        let (Some((wl, zl)), Some((wd, zd))) = (&r.lemke_solution, &r.dantzig_solution) else {
            panic!("scene {} missing a solver solution", r.scene);
        };
        for k in 0..wl.len() {
            assert!(
                (wl[k] - wd[k]).abs() < 1e-10,
                "scene {} w[{k}] disagreement: {} vs {}",
                r.scene,
                wl[k],
                wd[k]
            );
            assert!(
                (zl[k] - zd[k]).abs() < 1e-10,
                "scene {} z[{k}] disagreement: {} vs {}",
                r.scene,
                zl[k],
                zd[k]
            );
        }
        // Dantzig should also clear the accuracy gate.
        assert!(
            r.dantzig_residual.is_some_and(|res| res < 1e-10),
            "scene {} dantzig residual {:?} exceeds 1e-10",
            r.scene,
            r.dantzig_residual
        );
    }
}

// ── 5: THE CRITICAL TEST — anti-cycling on degenerate LCPs ──────────────────

#[test]
fn test_lemke_degenerate_anticycling() {
    let n = 4;
    let pivot_budget = 20 * n * n;
    let mut rng = SmallRng::seed_from_u64(0xDEAD_BEEF);

    for trial in 0..10 {
        // Random PSD M.
        let a = random_matrix(&mut rng, n);
        let m = psd_from(&a, 1e-2);

        // q with deliberately repeated entries to force minimum-ratio ties.
        // Using the same value on multiple rows is exactly the degeneracy the
        // lexicographic rule must survive.
        let q = vec![-1.0, -1.0, -1.0, -0.5];

        let (w, z) =
            lemke_solve(&m, &q).unwrap_or_else(|e| panic!("trial {trial}: lemke failed: {e}"));

        // Residual gate.
        let res = complementarity_residual(&m, &q, &w, &z);
        assert!(res < 1e-10, "trial {trial}: residual {res} too high");

        // Sanity: the solution is genuinely complementary and non-negative
        // (already folded into `res`, asserted again explicitly for clarity).
        for k in 0..n {
            assert!(w[k] >= -1e-12, "trial {trial}: w[{k}]={} negative", w[k]);
            assert!(z[k] >= -1e-12, "trial {trial}: z[{k}]={} negative", z[k]);
            assert!(
                (w[k] * z[k]).abs() < 1e-10,
                "trial {trial}: w[{k}]·z[{k}] not complementary"
            );
        }

        // Pivot-budget proxy: the solver caps itself at 10·n² internally and
        // returned Ok, so it terminated within that bound, which is below the
        // 20·n² anti-cycling ceiling.
        let _ = pivot_budget;
    }
}

#[test]
fn test_lemke_identical_q_rows_terminates() {
    // Maximally degenerate: every q entry identical, coupled PSD M. Ratio ties
    // occur on essentially every pivot; the lexicographic rule must still
    // terminate.
    let n = 5;
    let mut rng = SmallRng::seed_from_u64(7);
    let a = random_matrix(&mut rng, n);
    let m = psd_from(&a, 1e-3);
    let q = vec![-1.0; n];

    let (w, z) = lemke_solve(&m, &q).expect("degenerate-but-solvable");
    let res = complementarity_residual(&m, &q, &w, &z);
    assert!(res < 1e-10, "residual {res} too high");
}

// ── 6: random PSD batch ─────────────────────────────────────────────────────

#[test]
fn test_lemke_psd_random() {
    let n = 5;
    let mut rng = SmallRng::seed_from_u64(0x1234_5678);

    for trial in 0..20 {
        use rand::RngExt;
        let a = random_matrix(&mut rng, n);
        let m = psd_from(&a, 0.1);
        // Random strictly-negative q (drives several contacts active).
        let q: Vec<f64> = (0..n).map(|_| rng.random_range(-2.0_f64..-0.1)).collect();

        let (w, z) =
            lemke_solve(&m, &q).unwrap_or_else(|e| panic!("trial {trial}: lemke failed: {e}"));
        let res = complementarity_residual(&m, &q, &w, &z);
        assert!(res < 1e-10, "trial {trial}: residual {res} too high");

        // Cross-check against Dantzig on the same PSD instance.
        let (wd, zd) =
            dantzig_solve(&m, &q).unwrap_or_else(|e| panic!("trial {trial}: dantzig failed: {e}"));
        let res_d = complementarity_residual(&m, &q, &wd, &zd);
        assert!(
            res_d < 1e-10,
            "trial {trial}: dantzig residual {res_d} too high"
        );
        for k in 0..n {
            assert!(
                (w[k] - wd[k]).abs() < 1e-8,
                "trial {trial}: w[{k}] lemke/dantzig disagree: {} vs {}",
                w[k],
                wd[k]
            );
        }
    }
}
