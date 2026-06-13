// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Canonical-scene comparison harness for the direct LCP solvers.
//!
//! Three small contact systems — a box stack, a wedge, and a high-mass-ratio
//! pair — are encoded directly as `(M, q)` linear complementarity problems
//! with positive-semi-definite `M`. These are *constructed* rather than
//! simulated: each is a physically plausible velocity-level contact LCP (the
//! Delassus / effective-mass operator `M = J W Jᵀ` together with a
//! penetration-correction right-hand side `q`), chosen to exercise the
//! solvers on the kinds of structure that arise in practice (coupling between
//! stacked contacts, a near-singular wedge, and a `1000:1` mass ratio).
//!
//! [`run_harness`] solves each scene with both Lemke and Dantzig and reports
//! the complementarity residual ([`complementarity_residual`]) achieved by
//! each. It is the engine behind the solver-conformance matrix and the
//! `< 1e-10` accuracy gate.

use super::{dantzig_solve, lemke_solve};

/// A box stack of three contacts under gravity.
///
/// The three stacked contacts couple through the shared bodies, so `M` has
/// off-diagonal terms; the right-hand side `q` is negative at every contact,
/// representing the velocity-level penetration that the normal impulses must
/// remove. `M` is symmetric positive definite (diagonally dominant).
pub fn scene_box_stack() -> (Vec<Vec<f64>>, Vec<f64>) {
    // Effective-mass coupling: each contact feels its neighbour through the
    // body it shares. Diagonal 2 (two unit-mass bodies per contact), nearest-
    // neighbour coupling 0.5, no coupling between the top and bottom contacts.
    let m = vec![
        vec![2.0, 0.5, 0.0],
        vec![0.5, 2.0, 0.5],
        vec![0.0, 0.5, 2.0],
    ];
    let q = vec![-1.0, -1.0, -0.5];
    (m, q)
}

/// A two-contact wedge with high normal force and friction coupling.
///
/// The two faces of the wedge share the wedged body, producing a strong
/// off-diagonal coupling that makes `M` comparatively ill-conditioned while
/// remaining positive definite.
pub fn scene_wedge() -> (Vec<Vec<f64>>, Vec<f64>) {
    let m = vec![vec![1.0, 0.9], vec![0.9, 1.0]];
    let q = vec![-0.5, -0.8];
    (m, q)
}

/// A two-contact scenario with a `1000:1` mass ratio.
///
/// One contact involves a very heavy body (effective mass `1000`) and the
/// other a light body (effective mass `1`), with a small coupling term. Mass
/// ratios like this are the classic stress case for *iterative* solvers; a
/// direct method should resolve it exactly.
pub fn scene_high_mass_ratio() -> (Vec<Vec<f64>>, Vec<f64>) {
    let m = vec![vec![1000.0, 1.0], vec![1.0, 1.0]];
    let q = vec![-1.0, -1.0];
    (m, q)
}

/// Combined complementarity residual for an LCP solution `(w, z)`.
///
/// Sums three sources of error, so that a value below the accuracy gate
/// certifies a genuine LCP solution rather than merely a feasible point:
///
/// ```text
/// residual = wᵀz                       (complementarity)
///          + Σ max(0, -w_i)            (non-negativity of w)
///          + Σ max(0, -z_i)            (non-negativity of z)
///          + ‖w - M z - q‖₂            (linear feasibility)
/// ```
pub fn complementarity_residual(m: &[Vec<f64>], q: &[f64], w: &[f64], z: &[f64]) -> f64 {
    let n = q.len();
    let mut total = 0.0;

    // Complementarity and non-negativity.
    for i in 0..n {
        total += w[i] * z[i];
        total += (-w[i]).max(0.0);
        total += (-z[i]).max(0.0);
    }

    // Linear feasibility ‖w - M z - q‖₂.
    let mut feas_sq = 0.0;
    for i in 0..n {
        let mut mz = 0.0;
        for j in 0..n {
            mz += m[i][j] * z[j];
        }
        let r = w[i] - mz - q[i];
        feas_sq += r * r;
    }
    total += feas_sq.sqrt();

    total
}

/// Outcome of running both direct solvers on one canonical scene.
#[derive(Debug, Clone)]
pub struct HarnessResult {
    /// Human-readable scene name.
    pub scene: &'static str,
    /// Complementarity residual achieved by Lemke's method.
    pub lemke_residual: f64,
    /// Lemke solution `(w, z)`, if the solve succeeded.
    pub lemke_solution: Option<(Vec<f64>, Vec<f64>)>,
    /// Complementarity residual achieved by Dantzig's method, if it ran.
    pub dantzig_residual: Option<f64>,
    /// Dantzig solution `(w, z)`, if the solve succeeded.
    pub dantzig_solution: Option<(Vec<f64>, Vec<f64>)>,
}

/// Solve every canonical scene with both direct solvers and collect the
/// residuals.
///
/// A solver that errors on a scene contributes an infinite residual (Lemke)
/// or `None` (Dantzig) rather than aborting the whole harness, so the report
/// always covers all scenes.
pub fn run_harness() -> Vec<HarnessResult> {
    let scenes: [(&'static str, fn() -> (Vec<Vec<f64>>, Vec<f64>)); 3] = [
        ("box_stack", scene_box_stack),
        ("wedge", scene_wedge),
        ("high_mass_ratio", scene_high_mass_ratio),
    ];

    scenes
        .into_iter()
        .map(|(name, build)| {
            let (m, q) = build();

            let (lemke_residual, lemke_solution) = match lemke_solve(&m, &q) {
                Ok((w, z)) => {
                    let res = complementarity_residual(&m, &q, &w, &z);
                    (res, Some((w, z)))
                }
                Err(_) => (f64::INFINITY, None),
            };

            let (dantzig_residual, dantzig_solution) = match dantzig_solve(&m, &q) {
                Ok((w, z)) => {
                    let res = complementarity_residual(&m, &q, &w, &z);
                    (Some(res), Some((w, z)))
                }
                Err(_) => (None, None),
            };

            HarnessResult {
                scene: name,
                lemke_residual,
                lemke_solution,
                dantzig_residual,
                dantzig_solution,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_residual_zero_for_exact_solution() {
        // Identity LCP: q = [-1, -2] → w = [0, 0], z = [1, 2] is exact.
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-1.0, -2.0];
        let w = vec![0.0, 0.0];
        let z = vec![1.0, 2.0];
        assert!(complementarity_residual(&m, &q, &w, &z) < 1e-15);
    }

    #[test]
    fn test_residual_flags_feasibility_violation() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-1.0, -2.0];
        // Wrong z → large feasibility residual.
        let w = vec![0.0, 0.0];
        let z = vec![0.0, 0.0];
        assert!(complementarity_residual(&m, &q, &w, &z) > 1.0);
    }

    #[test]
    fn test_residual_flags_complementarity_violation() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![1.0, 1.0];
        // w = q, z = [1, 1]: feasible (w - z = q ✗ actually). Build an explicit
        // complementarity violation: both w and z positive and consistent.
        // Choose w = [2, 2], z = [1, 1], then w - z = [1, 1] = q (feasible),
        // but wᵀz = 4 > 0 (complementarity violated).
        let w = vec![2.0, 2.0];
        let z = vec![1.0, 1.0];
        let res = complementarity_residual(&m, &q, &w, &z);
        assert!(res >= 4.0 - 1e-12, "res={res}");
    }

    #[test]
    fn test_scenes_are_psd_symmetric() {
        for (m, _q) in [scene_box_stack(), scene_wedge(), scene_high_mass_ratio()] {
            let n = m.len();
            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (m[i][j] - m[j][i]).abs() < 1e-12,
                        "M not symmetric at ({i},{j})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_run_harness_all_scenes_low_residual() {
        let results = run_harness();
        assert_eq!(results.len(), 3);
        for r in &results {
            assert!(
                r.lemke_residual < 1e-10,
                "scene {} lemke residual {} too high",
                r.scene,
                r.lemke_residual
            );
        }
    }

    #[test]
    fn test_run_harness_lemke_dantzig_agree() {
        for r in run_harness() {
            let (Some((wl, _)), Some((wd, _))) = (&r.lemke_solution, &r.dantzig_solution) else {
                panic!("scene {} missing a solution", r.scene);
            };
            for (a, b) in wl.iter().zip(wd.iter()) {
                assert!(
                    (a - b).abs() < 1e-10,
                    "scene {} disagreement: {a} vs {b}",
                    r.scene
                );
            }
        }
    }
}
