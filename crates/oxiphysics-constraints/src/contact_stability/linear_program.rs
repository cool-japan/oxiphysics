// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Small dense linear-program solver for the exact force-closure test.
//!
//! Solves
//!
//! ```text
//! maximise   c · x
//! subject to A x = b,  x >= 0
//! ```
//!
//! with a two-phase tableau simplex method.  Bland's smallest-index rule is
//! used for both the entering and the leaving variable, so the method cannot
//! cycle on the (highly degenerate) problems that arise in contact analysis.
//! Those problems have a handful of rows and up to a few hundred columns, for
//! which a dense tableau is simple and robust.

/// Pivot elements at or below this magnitude are treated as zero.
const PIVOT_TOL: f64 = 1e-10;

/// Reduced costs above `-COST_TOL` count as non-negative (optimality).
const COST_TOL: f64 = 1e-10;

/// A phase-I objective below `-FEASIBILITY_TOL` means the problem is infeasible.
const FEASIBILITY_TOL: f64 = 1e-9;

/// Relative tolerance below which two ratio-test quotients count as tied.
const RATIO_TIE_TOL: f64 = 1e-12;

/// Termination status of [`maximise`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum LpStatus {
    /// An optimal solution was found.
    Optimal,
    /// No `x >= 0` satisfies `A x = b`.
    Infeasible,
    /// The objective is unbounded above on the feasible set.
    Unbounded,
    /// The iteration budget was exhausted (numerically stalled problem).
    IterationLimit,
    /// The dimensions of `c`, `A` and `b` are inconsistent.
    InvalidInput,
}

/// Result of [`maximise`].
#[derive(Debug, Clone)]
pub(super) struct LpSolution {
    /// Termination status; `objective` and `x` are meaningful only for
    /// [`LpStatus::Optimal`].
    pub(super) status: LpStatus,
    /// Optimal objective value `c · x`.
    pub(super) objective: f64,
    /// Optimal point (`n` entries).
    pub(super) x: Vec<f64>,
}

impl LpSolution {
    fn failure(status: LpStatus) -> Self {
        Self {
            status,
            objective: 0.0,
            x: Vec::new(),
        }
    }
}

/// Dense simplex tableau.
///
/// `t` has `m + 1` rows of `ncols + 1` entries: rows `0..m` are the
/// constraint rows, row `m` holds the reduced costs `z_j - c_j` (its last
/// entry is the current objective value), and the last column is the
/// right-hand side.
struct Tableau {
    t: Vec<Vec<f64>>,
    basis: Vec<usize>,
    m: usize,
    ncols: usize,
}

impl Tableau {
    /// Gauss-Jordan pivot on `(row, col)`.
    fn pivot(&mut self, row: usize, col: usize) {
        let inv = 1.0 / self.t[row][col];
        for value in self.t[row].iter_mut() {
            *value *= inv;
        }
        self.t[row][col] = 1.0;
        let pivot_row = self.t[row].clone();
        for (i, other) in self.t.iter_mut().enumerate() {
            if i == row {
                continue;
            }
            let factor = other[col];
            if factor == 0.0 {
                continue;
            }
            for (value, p) in other.iter_mut().zip(pivot_row.iter()) {
                *value -= factor * p;
            }
            other[col] = 0.0;
        }
        self.basis[row] = col;
    }

    /// Simplex iterations with Bland's rule; only columns `0..enter_limit`
    /// may enter the basis.
    fn iterate(&mut self, enter_limit: usize, max_iter: usize) -> LpStatus {
        for _ in 0..max_iter {
            let Some(col) = (0..enter_limit).find(|&j| self.t[self.m][j] < -COST_TOL) else {
                return LpStatus::Optimal;
            };
            // Minimum-ratio test; ties go to the smallest basic variable.
            let mut best: Option<(usize, f64)> = None;
            for i in 0..self.m {
                let a = self.t[i][col];
                if a <= PIVOT_TOL {
                    continue;
                }
                let ratio = self.t[i][self.ncols].max(0.0) / a;
                best = match best {
                    None => Some((i, ratio)),
                    Some((bi, br)) => {
                        let tol = RATIO_TIE_TOL * br.abs().max(1.0);
                        if ratio < br - tol || (ratio <= br + tol && self.basis[i] < self.basis[bi])
                        {
                            Some((i, ratio))
                        } else {
                            Some((bi, br))
                        }
                    }
                };
            }
            let Some((row, _)) = best else {
                return LpStatus::Unbounded;
            };
            self.pivot(row, col);
        }
        LpStatus::IterationLimit
    }
}

/// Maximise `c · x` subject to `A x = b`, `x >= 0`.
///
/// `a` is given row by row (`m` rows of `n = c.len()` entries).  Rows may be
/// linearly dependent and `b` may have any sign.  Every row is equilibrated
/// internally so the fixed tolerances are meaningful for badly scaled input.
pub(super) fn maximise(c: &[f64], a: &[Vec<f64>], b: &[f64]) -> LpSolution {
    let n = c.len();
    let m = a.len();
    if b.len() != m || a.iter().any(|row| row.len() != n) {
        return LpSolution::failure(LpStatus::InvalidInput);
    }
    if c.iter().any(|v| !v.is_finite())
        || b.iter().any(|v| !v.is_finite())
        || a.iter().any(|row| row.iter().any(|v| !v.is_finite()))
    {
        return LpSolution::failure(LpStatus::InvalidInput);
    }

    // Columns: n structural variables followed by m artificial variables.
    let ncols = n + m;
    let max_iter = 200 * (m + ncols) + 1_000;
    let mut t = vec![vec![0.0_f64; ncols + 1]; m + 1];
    for (i, (row, &rhs)) in a.iter().zip(b.iter()).enumerate() {
        let magnitude = row.iter().fold(rhs.abs(), |acc, v| acc.max(v.abs()));
        let scale = if magnitude > 0.0 {
            1.0 / magnitude
        } else {
            1.0
        };
        let sign = if rhs < 0.0 { -scale } else { scale };
        for (j, value) in row.iter().enumerate() {
            t[i][j] = sign * value;
        }
        t[i][n + i] = 1.0;
        t[i][ncols] = sign * rhs;
    }
    let mut tab = Tableau {
        t,
        basis: (n..n + m).collect(),
        m,
        ncols,
    };

    // ── Phase I: maximise -(sum of artificial variables) ─────────────────────
    for j in 0..ncols {
        let mut reduced = 0.0;
        for i in 0..m {
            reduced -= tab.t[i][j];
        }
        // Artificial columns have cost -1: z_j - c_j = -1 + 1 = 0.
        tab.t[m][j] = if j < n { reduced } else { 0.0 };
    }
    let mut phase_one_value = 0.0;
    for i in 0..m {
        phase_one_value -= tab.t[i][ncols];
    }
    tab.t[m][ncols] = phase_one_value;

    match tab.iterate(ncols, max_iter) {
        LpStatus::Optimal => {}
        LpStatus::Unbounded => return LpSolution::failure(LpStatus::InvalidInput),
        other => return LpSolution::failure(other),
    }
    if tab.t[m][ncols] < -FEASIBILITY_TOL {
        return LpSolution::failure(LpStatus::Infeasible);
    }

    // Drive artificial variables out of the basis.  A row without a usable
    // structural pivot is redundant and keeps its (zero) artificial variable.
    for i in 0..m {
        if tab.basis[i] < n {
            continue;
        }
        let mut best_col = None;
        let mut best_abs = PIVOT_TOL;
        for j in 0..n {
            let magnitude = tab.t[i][j].abs();
            if magnitude > best_abs {
                best_abs = magnitude;
                best_col = Some(j);
            }
        }
        if let Some(j) = best_col {
            tab.pivot(i, j);
        }
    }

    // ── Phase II: the real objective ──────────────────────────────────────────
    for j in 0..ncols {
        let mut reduced = if j < n { -c[j] } else { 0.0 };
        for i in 0..m {
            let cost = if tab.basis[i] < n {
                c[tab.basis[i]]
            } else {
                0.0
            };
            reduced += cost * tab.t[i][j];
        }
        tab.t[m][j] = reduced;
    }
    let mut objective_value = 0.0;
    for i in 0..m {
        if tab.basis[i] < n {
            objective_value += c[tab.basis[i]] * tab.t[i][ncols];
        }
    }
    tab.t[m][ncols] = objective_value;

    match tab.iterate(n, max_iter) {
        LpStatus::Optimal => {}
        other => return LpSolution::failure(other),
    }

    let mut x = vec![0.0_f64; n];
    for i in 0..m {
        if tab.basis[i] < n {
            x[tab.basis[i]] = tab.t[i][ncols].max(0.0);
        }
    }
    let objective = c.iter().zip(x.iter()).map(|(ci, xi)| ci * xi).sum();
    LpSolution {
        status: LpStatus::Optimal,
        objective,
        x,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(data: &[&[f64]]) -> Vec<Vec<f64>> {
        data.iter().map(|r| r.to_vec()).collect()
    }

    #[test]
    fn textbook_problem_reaches_the_known_optimum() {
        // max 3x + 2y  s.t.  x + y <= 4, x + 3y <= 6, x <= 3   (slacks s1..s3)
        let a = rows(&[
            &[1.0, 1.0, 1.0, 0.0, 0.0],
            &[1.0, 3.0, 0.0, 1.0, 0.0],
            &[1.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        let sol = maximise(&[3.0, 2.0, 0.0, 0.0, 0.0], &a, &[4.0, 6.0, 3.0]);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(
            (sol.objective - 11.0).abs() < 1e-9,
            "obj = {}",
            sol.objective
        );
        assert!((sol.x[0] - 3.0).abs() < 1e-9 && (sol.x[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn inconsistent_equalities_are_infeasible() {
        let a = rows(&[&[1.0, 1.0], &[1.0, 1.0]]);
        let sol = maximise(&[1.0, 0.0], &a, &[1.0, 2.0]);
        assert_eq!(sol.status, LpStatus::Infeasible);
    }

    #[test]
    fn negative_right_hand_side_has_no_nonnegative_solution() {
        // x + y = -1 with x, y >= 0 is infeasible.
        let sol = maximise(&[1.0, 1.0], &rows(&[&[1.0, 1.0]]), &[-1.0]);
        assert_eq!(sol.status, LpStatus::Infeasible);
    }

    #[test]
    fn negative_right_hand_side_is_normalised() {
        // -x - y = -1  <=>  x + y = 1; maximise x.
        let sol = maximise(&[1.0, 0.0], &rows(&[&[-1.0, -1.0]]), &[-1.0]);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!((sol.objective - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unbounded_objective_is_detected() {
        // max x  s.t.  x - s = 1, x, s >= 0.
        let sol = maximise(&[1.0, 0.0], &rows(&[&[1.0, -1.0]]), &[1.0]);
        assert_eq!(sol.status, LpStatus::Unbounded);
    }

    #[test]
    fn redundant_rows_are_tolerated() {
        // x + y = 1 stated twice (the second time scaled): maximise x.
        let a = rows(&[&[1.0, 1.0], &[2.0, 2.0]]);
        let sol = maximise(&[1.0, 0.0], &a, &[1.0, 2.0]);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(
            (sol.objective - 1.0).abs() < 1e-12,
            "obj = {}",
            sol.objective
        );
    }

    #[test]
    fn beale_cycling_example_terminates_at_the_optimum() {
        // Beale's classic example on which Dantzig's rule cycles; Bland's rule
        // must reach the optimum 5/4.
        let a = rows(&[
            &[1.0, 0.0, 0.0, 0.25, -8.0, -1.0, 9.0],
            &[0.0, 1.0, 0.0, 0.5, -12.0, -0.5, 3.0],
            &[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        ]);
        let c = [0.0, 0.0, 0.0, 0.75, -20.0, 0.5, -6.0];
        let sol = maximise(&c, &a, &[0.0, 0.0, 1.0]);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(
            (sol.objective - 1.25).abs() < 1e-9,
            "obj = {}",
            sol.objective
        );
    }

    #[test]
    fn badly_scaled_rows_are_equilibrated() {
        // Same feasible set at wildly different row scales.
        let a = rows(&[&[1e8, 1e8, 0.0], &[0.0, 1e-8, 1e-8]]);
        let sol = maximise(&[1.0, 0.0, 0.0], &a, &[1e8, 1e-8]);
        assert_eq!(sol.status, LpStatus::Optimal);
        assert!(
            (sol.objective - 1.0).abs() < 1e-9,
            "obj = {}",
            sol.objective
        );
    }

    #[test]
    fn empty_constraint_set_is_bounded_only_for_non_positive_costs() {
        let none: Vec<Vec<f64>> = Vec::new();
        let bounded = maximise(&[-1.0, 0.0], &none, &[]);
        assert_eq!(bounded.status, LpStatus::Optimal);
        assert!(bounded.objective.abs() < 1e-15);
        let unbounded = maximise(&[1.0, 0.0], &none, &[]);
        assert_eq!(unbounded.status, LpStatus::Unbounded);
    }

    #[test]
    fn malformed_input_is_reported() {
        let sol = maximise(&[1.0, 2.0], &rows(&[&[1.0]]), &[1.0]);
        assert_eq!(sol.status, LpStatus::InvalidInput);
        let sol = maximise(&[1.0], &rows(&[&[f64::NAN]]), &[1.0]);
        assert_eq!(sol.status, LpStatus::InvalidInput);
        let sol = maximise(&[1.0], &rows(&[&[1.0]]), &[1.0, 2.0]);
        assert_eq!(sol.status, LpStatus::InvalidInput);
    }
}
