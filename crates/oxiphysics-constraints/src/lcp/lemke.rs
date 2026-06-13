// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lemke's complementary pivoting algorithm for the LCP `w - M z = q`.
//!
//! The method (Lemke, 1965; presentation following Cottle, Pang & Stone,
//! *The Linear Complementarity Problem*, 1992) augments the system with a
//! covering variable `z₀` and an all-ones covering vector `d`:
//!
//! ```text
//! w - M z - d z₀ = q ,   w, z, z₀ ≥ 0 .
//! ```
//!
//! Starting from the (infeasible) basis `w = q`, it drives `z₀` in far
//! enough to make every basic variable non-negative, then performs a
//! sequence of *complementary* pivots — each one bringing in the complement
//! of the variable that just left — until `z₀` itself leaves the basis. At
//! that point `z₀ = 0`, complementarity `wᵀz = 0` holds because `w_i` and
//! `z_i` are never simultaneously basic, and `(w, z)` solves the LCP.
//!
//! Degenerate problems (ties in the minimum-ratio test) are handled by the
//! lexicographic rule in [`LcpTableau::min_ratio_lex`], which guarantees a
//! unique leaving variable and therefore termination without cycling.

use super::LcpError;
use super::dense::{EPS, LcpTableau};

/// Solve the LCP: find `w, z ≥ 0` with `w - M z = q` and `wᵀ z = 0`.
///
/// Returns `(w, z)` on success.
///
/// # Errors
///
/// * [`LcpError::DimensionMismatch`] if `m` is not `q.len() × q.len()`.
/// * [`LcpError::RayTermination`] if a secondary ray is encountered (the
///   problem is infeasible or unbounded; for a copositive-plus `M` this means
///   no solution exists).
/// * [`LcpError::MaxPivots`] if the pivot budget `10·n²` is exhausted, which
///   for a correct lexicographic implementation indicates a pathological
///   input rather than ordinary cycling.
///
/// # Examples
///
/// ```
/// use oxiphysics_constraints::lcp::lemke_solve;
///
/// // M = [[2, 1], [1, 2]], q = [-1, -1]  →  z = [1/3, 1/3], w = [0, 0].
/// let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
/// let q = vec![-1.0, -1.0];
/// let (w, z) = lemke_solve(&m, &q).expect("solvable PSD LCP");
/// assert!((z[0] - 1.0 / 3.0).abs() < 1e-12);
/// assert!((z[1] - 1.0 / 3.0).abs() < 1e-12);
/// assert!(w[0].abs() < 1e-12 && w[1].abs() < 1e-12);
/// ```
pub fn lemke_solve(m: &[Vec<f64>], q: &[f64]) -> Result<(Vec<f64>, Vec<f64>), LcpError> {
    let n = q.len();
    if m.len() != n || m.iter().any(|row| row.len() != n) {
        return Err(LcpError::DimensionMismatch { m: m.len(), n });
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    // Trivial case: q ≥ 0 ⇒ w = q, z = 0 already solves the LCP.
    if q.iter().all(|&qi| qi >= -EPS) {
        let w = q.iter().map(|&qi| qi.max(0.0)).collect();
        return Ok((w, vec![0.0; n]));
    }

    let mut tableau = LcpTableau::new(m, q);
    let z0_col = tableau.z0_col();
    let max_pivots = 10 * n * n;

    // ── Initial pivot: drive z₀ into the basis ────────────────────────────
    // The leaving row is the one with the most-negative RHS: with the all-ones
    // covering vector this is exactly argmin_i q_i / d_i = argmin_i q_i, and it
    // is the choice that lifts every basic variable to non-negativity in one
    // step. (This is Lemke's dedicated initialisation rule, distinct from the
    // ordinary ratio test, which does not apply while the basis is infeasible.)
    let leaving_row = most_negative_rhs_row(&tableau);
    // The variable leaving on this initial pivot is w at `leaving_row`
    // (basis still canonical), so its complement z enters next.
    let mut leaving_var = tableau.basic_var(leaving_row);
    tableau.pivot(leaving_row, z0_col);

    // ── Complementary pivoting loop ───────────────────────────────────────
    for _ in 0..max_pivots {
        // The variable that just left dictates the entering variable: bring
        // in its complement. If z₀ left, we are done.
        if leaving_var == z0_col {
            return Ok(tableau.extract_solution());
        }
        let entering = complement(leaving_var, n);

        // Minimum-ratio test (lexicographic) selects the leaving row.
        let Some(row) = tableau.min_ratio_lex(entering) else {
            return Err(LcpError::RayTermination);
        };
        leaving_var = tableau.basic_var(row);
        tableau.pivot(row, entering);

        // If z₀ was the variable we just displaced, the next loop iteration
        // detects termination.
    }

    Err(LcpError::MaxPivots { max_pivots })
}

/// Row index of the most-negative right-hand side (Lemke's initial-pivot
/// rule for the all-ones covering vector). Ties are broken by the lowest
/// row index, which is deterministic and harmless for the initialisation.
fn most_negative_rhs_row(tableau: &LcpTableau) -> usize {
    let rhs = tableau.rhs_col();
    let mut best_row = 0;
    let mut best_val = tableau.at(0, rhs);
    for r in 1..tableau.dim() {
        let v = tableau.at(r, rhs);
        if v < best_val {
            best_val = v;
            best_row = r;
        }
    }
    best_row
}

/// Complement of a variable column under the LCP pairing.
///
/// `w_i` (column `i`) ↔ `z_i` (column `n + i`). The covering variable `z₀`
/// has no complement and is handled by the termination check before this is
/// ever called on it.
#[inline]
fn complement(var: usize, n: usize) -> usize {
    if var < n { var + n } else { var - n }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// L2 norm of `w - M z - q`, plus complementarity and non-negativity
    /// violations, for an independent residual check.
    fn lcp_residual(m: &[Vec<f64>], q: &[f64], w: &[f64], z: &[f64]) -> f64 {
        let n = q.len();
        let mut feas_sq = 0.0;
        for i in 0..n {
            let mut mz = 0.0;
            for j in 0..n {
                mz += m[i][j] * z[j];
            }
            let r = w[i] - mz - q[i];
            feas_sq += r * r;
        }
        let mut comp = 0.0;
        for i in 0..n {
            comp += w[i] * z[i];
            comp += (-w[i]).max(0.0);
            comp += (-z[i]).max(0.0);
        }
        comp + feas_sq.sqrt()
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(0, 3), 3);
        assert_eq!(complement(1, 3), 4);
        assert_eq!(complement(3, 3), 0);
        assert_eq!(complement(5, 3), 2);
    }

    #[test]
    fn test_dimension_mismatch() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-1.0]; // wrong length
        let err = lemke_solve(&m, &q).unwrap_err();
        assert_eq!(err, LcpError::DimensionMismatch { m: 2, n: 1 });
    }

    #[test]
    fn test_trivial_nonnegative_q() {
        // q ≥ 0 ⇒ immediate w = q, z = 0.
        let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let q = vec![1.0, 3.0];
        let (w, z) = lemke_solve(&m, &q).expect("trivial");
        assert!((w[0] - 1.0).abs() < 1e-12);
        assert!((w[1] - 3.0).abs() < 1e-12);
        assert!(z[0].abs() < 1e-12 && z[1].abs() < 1e-12);
    }

    #[test]
    fn test_empty_problem() {
        let m: Vec<Vec<f64>> = Vec::new();
        let q: Vec<f64> = Vec::new();
        let (w, z) = lemke_solve(&m, &q).expect("empty ok");
        assert!(w.is_empty() && z.is_empty());
    }

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
        assert!(lcp_residual(&m, &q, &w, &z) < 1e-12);
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
        assert!((z[0] - 1.0).abs() < 1e-12);
        assert!((z[1] - 2.0).abs() < 1e-12);
        assert!((z[2] - 3.0).abs() < 1e-12);
        for wi in &w {
            assert!(wi.abs() < 1e-12);
        }
    }

    #[test]
    fn test_lemke_mixed_active_inactive() {
        // M = I, q = [-2, 3, -1]: contacts 0 and 2 active, contact 1 already
        // separated → z = [2, 0, 1], w = [0, 3, 0].
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let q = vec![-2.0, 3.0, -1.0];
        let (w, z) = lemke_solve(&m, &q).expect("solvable");
        assert!((z[0] - 2.0).abs() < 1e-12, "z0={}", z[0]);
        assert!(z[1].abs() < 1e-12, "z1={}", z[1]);
        assert!((z[2] - 1.0).abs() < 1e-12, "z2={}", z[2]);
        assert!(w[0].abs() < 1e-12);
        assert!((w[1] - 3.0).abs() < 1e-12, "w1={}", w[1]);
        assert!(w[2].abs() < 1e-12);
        assert!(lcp_residual(&m, &q, &w, &z) < 1e-12);
    }
}
