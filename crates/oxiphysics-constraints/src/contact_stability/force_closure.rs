// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact force-closure test.
//!
//! A wrench set `W = {w_1, .., w_n}` achieves **force closure** when the cone
//! of non-negative combinations `sum(lambda_i * w_i)` is the whole wrench space
//! `R^6`, i.e. when the contacts can balance an arbitrary external wrench.
//! By the positive-spanning theorem (Murray, Li & Sastry, *A Mathematical
//! Introduction to Robotic Manipulation*, 1994) this holds iff
//!
//! 1. the wrenches span `R^6` (`rank W = 6`), **and**
//! 2. a *strictly positive* combination vanishes:
//!    `sum(lambda_i * w_i) = 0` with every `lambda_i > 0`.
//!
//! Condition 1 is a rank computation; condition 2 is decided by the linear
//! program
//!
//! ```text
//! maximise   t
//! subject to sum_i lambda_i * w_i = 0,   sum_i lambda_i = 1,   lambda_i >= t
//! ```
//!
//! which has a solution with `t > 0` exactly when such a combination exists.
//! The wrenches are scaled to unit length first (positive scaling does not
//! change the positive span), so the optimum is bounded by `1 / n` and the
//! relative margin `n * t` is a meaningful scale-free quantity.

use super::linear_program::{LpStatus, maximise};
use super::wrench_set::distinct_unit;
use super::{Wrench, force_closure_test, vec6_len};

/// Pivots of the rank elimination at or below this magnitude count as zero
/// (the wrenches have unit length, so entries are bounded by one).
const RANK_TOL: f64 = 1e-9;

/// Smallest relative closure margin `n * t` that counts as strictly positive.
const CLOSURE_MARGIN_TOL: f64 = 1e-9;

/// A wrench set needs at least `dim + 1` members to span `R^dim` positively.
const MIN_WRENCHES: usize = 7;

/// Rank of a wrench set, by Gaussian elimination with complete pivoting.
fn wrench_rank(wrenches: &[Wrench]) -> usize {
    let mut rows: Vec<[f64; 6]> = wrenches.to_vec();
    let mut used_column = [false; 6];
    let mut rank = 0;
    while rank < rows.len() && rank < 6 {
        // Largest remaining entry among the unused columns.
        let mut pivot_row = rank;
        let mut pivot_col = 0;
        let mut pivot_abs = 0.0_f64;
        for (r, row) in rows.iter().enumerate().skip(rank) {
            for (c, value) in row.iter().enumerate() {
                if !used_column[c] && value.abs() > pivot_abs {
                    pivot_row = r;
                    pivot_col = c;
                    pivot_abs = value.abs();
                }
            }
        }
        if pivot_abs <= RANK_TOL {
            break;
        }
        rows.swap(rank, pivot_row);
        used_column[pivot_col] = true;
        let pivot = rows[rank];
        for row in rows.iter_mut().skip(rank + 1) {
            let factor = row[pivot_col] / pivot[pivot_col];
            for (value, p) in row.iter_mut().zip(pivot.iter()) {
                *value -= factor * p;
            }
            row[pivot_col] = 0.0;
        }
        rank += 1;
    }
    rank
}

/// Largest tolerated violation of `sum(lambda_i * w_i) = 0` and `sum(lambda_i) = 1`
/// by the LP solution before it is distrusted.
const CERTIFICATE_TOL: f64 = 1e-7;

/// Closure margin `n * min_i(lambda_i)` of the best combination found by the
/// strictly-positive-combination LP for a set of unit wrenches, or `None` if
/// the LP did not reach a verified optimum.
///
/// The returned weights are re-checked against the equations they are supposed
/// to satisfy, so a numerically broken simplex run cannot certify closure.
fn closure_margin(unit: &[Wrench]) -> Option<f64> {
    let n = unit.len();
    // Variables: x = (t, s_1, .., s_n) >= 0 with lambda_i = t + s_i, so that
    // lambda_i >= t holds by construction.
    //   sum_i (t + s_i) w_i = 0      (six equations)
    //   sum_i (t + s_i)     = 1      (one equation)
    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(7);
    for axis in 0..6 {
        let mut row = Vec::with_capacity(n + 1);
        row.push(unit.iter().map(|w| w[axis]).sum());
        row.extend(unit.iter().map(|w| w[axis]));
        rows.push(row);
    }
    let mut normalisation = vec![1.0; n + 1];
    normalisation[0] = n as f64;
    rows.push(normalisation);

    let mut rhs = vec![0.0; 7];
    rhs[6] = 1.0;
    let mut cost = vec![0.0; n + 1];
    cost[0] = 1.0;

    let solution = maximise(&cost, &rows, &rhs);
    if solution.status != LpStatus::Optimal || solution.x.len() != n + 1 {
        return None;
    }
    let t = solution.x[0];
    let mut balance = [0.0_f64; 6];
    let mut total = 0.0;
    let mut smallest = f64::INFINITY;
    for (w, s) in unit.iter().zip(solution.x.iter().skip(1)) {
        let lambda = t + s;
        total += lambda;
        smallest = smallest.min(lambda);
        for (acc, component) in balance.iter_mut().zip(w.iter()) {
            *acc += lambda * component;
        }
    }
    if (total - 1.0).abs() > CERTIFICATE_TOL || vec6_len(&balance) > CERTIFICATE_TOL {
        return None;
    }
    // At an optimum `t` cannot be raised, so the smallest weight equals it.
    if (smallest - solution.objective).abs() > CERTIFICATE_TOL {
        return None;
    }
    Some(n as f64 * smallest)
}

/// Exact force-closure test.
///
/// Returns `true` iff non-negative combinations of `wrenches` can produce
/// **every** wrench in `R^6`: the wrenches span `R^6` and some strictly
/// positive combination of them is zero (equivalently, the origin lies in the
/// interior of their convex hull).
///
/// Unlike [`force_closure_test`], which only checks that every wrench axis is
/// reached with both signs (a necessary condition that several non-closing
/// configurations pass), this test is conclusive:
///
/// * fewer than seven distinct non-zero wrenches, a rank-deficient set, a
///   non-finite entry, or an optimum margin at or below `1e-9` all give
///   `false`;
/// * the rank is computed by Gaussian elimination with complete pivoting and
///   the strict-positivity condition by a two-phase simplex method
///   (Bland's rule, so it cannot cycle on degenerate contact geometry).
///
/// Each wrench is normalised to unit length first, so the decision does not
/// depend on the magnitudes of the individual wrenches.  Force closure is also
/// independent of the torque reference point (moving the reference point is an
/// invertible linear map of wrench space); the quality metrics
/// [`grasp_quality_q1`](super::grasp_quality_q1) and
/// [`grasp_quality_epsilon`](super::grasp_quality_epsilon) are not.
pub fn force_closure_exact(wrenches: &[Wrench]) -> bool {
    let Some(unit) = distinct_unit(wrenches) else {
        return false;
    };
    if unit.len() < MIN_WRENCHES || !force_closure_test(&unit) {
        return false;
    }
    if wrench_rank(&unit) < 6 {
        return false;
    }
    closure_margin(&unit).is_some_and(|margin| margin > CLOSURE_MARGIN_TOL)
}
