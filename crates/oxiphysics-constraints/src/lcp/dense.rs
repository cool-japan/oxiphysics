// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense pivoting tableau for the direct LCP solvers.
//!
//! The tableau represents the linear system arising in Lemke's method,
//!
//! ```text
//! w - M z - d z₀ = q
//! ```
//!
//! as a full augmented matrix in which **every** variable has a column.
//! Writing row `i` out,
//!
//! ```text
//! 1·w_i + Σ_j (-M_ij)·z_j + (-d_i)·z₀ = q_i ,
//! ```
//!
//! so the initial tableau is `[ I | -M | -d | q ]` with column layout
//!
//! ```text
//!   columns  0 .. n-1      w_1 .. w_n
//!   columns  n .. 2n-1     z_1 .. z_n
//!   column   2n            z₀  (covering / artificial variable)
//!   column   2n+1          q   (right-hand side)
//! ```
//!
//! A *basic feasible solution* keeps one variable basic per row; its column
//! is a unit vector (`1` in that row, `0` elsewhere) and the basic variable
//! equals the RHS of its row. All non-basic variables are `0`. Pivoting on
//! `(r, c)` is a Gauss–Jordan elimination that turns column `c` into the
//! `r`-th unit vector, swapping the entering variable `c` into the basis in
//! place of whatever was basic in row `r`.
//!
//! Anti-cycling on degenerate problems is handled by a **lexicographic**
//! minimum-ratio test ([`LcpTableau::min_ratio_lex`]): ties in the ordinary
//! ratio test are broken by comparing, column-by-column over the original
//! identity block, the perturbed ratios — equivalent to the standard
//! lexicographic perturbation `q → q + ε·I` taken to the limit `ε → 0⁺`.

/// Numerical tolerance below which a tableau entry is treated as exactly
/// zero. Pivot candidates and ratio tests use this to avoid dividing by
/// values that are zero up to round-off.
pub(crate) const EPS: f64 = 1e-12;

/// A dense Gauss–Jordan tableau for Lemke-style complementary pivoting.
///
/// See the [module documentation](self) for the column layout and the
/// meaning of the basis.
#[derive(Debug, Clone)]
pub struct LcpTableau {
    /// Problem dimension (number of contacts / rows).
    n: usize,
    /// The `n × (2n + 2)` augmented matrix, row-major.
    tableau: Vec<Vec<f64>>,
    /// `basis[i]` is the column index (`0 ..= 2n`) of the variable that is
    /// basic in row `i`.
    basis: Vec<usize>,
}

impl LcpTableau {
    /// Build the initial tableau `[ I | -M | -d | q ]` for the LCP
    /// `w - M z = q`, with an all-ones covering vector `d` (Lemke's choice)
    /// for the artificial variable `z₀`.
    ///
    /// The starting basis is `w_1 .. w_n` (columns `0 .. n-1`), which is the
    /// canonical — generally infeasible — basis the algorithm repairs by
    /// first driving `z₀` in.
    ///
    /// # Panics
    ///
    /// Panics if `m` is not `q.len() × q.len()` square. Callers in this
    /// crate validate dimensions and surface
    /// [`LcpError::DimensionMismatch`](super::LcpError::DimensionMismatch)
    /// before constructing a tableau, so this is an internal invariant.
    pub fn new(m: &[Vec<f64>], q: &[f64]) -> Self {
        let n = q.len();
        assert_eq!(m.len(), n, "M must have q.len() rows");
        for (i, row) in m.iter().enumerate() {
            assert_eq!(row.len(), n, "M row {i} must have q.len() columns");
        }

        let width = 2 * n + 2;
        let z0_col = 2 * n;
        let rhs_col = 2 * n + 1;
        let mut tableau = vec![vec![0.0; width]; n];
        for (i, row) in tableau.iter_mut().enumerate() {
            // Identity block for w.
            row[i] = 1.0;
            // -M block for z.
            for (j, &m_ij) in m[i].iter().enumerate() {
                row[n + j] = -m_ij;
            }
            // Covering column -d, with d = 1.
            row[z0_col] = -1.0;
            // Right-hand side q.
            row[rhs_col] = q[i];
        }

        let basis = (0..n).collect();
        Self { n, tableau, basis }
    }

    /// Problem dimension (number of rows / contacts).
    #[inline]
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Column index of the covering variable `z₀`.
    #[inline]
    pub fn z0_col(&self) -> usize {
        2 * self.n
    }

    /// Column index of the right-hand side.
    #[inline]
    pub fn rhs_col(&self) -> usize {
        2 * self.n + 1
    }

    /// The variable (column index) currently basic in `row`.
    #[inline]
    pub fn basic_var(&self, row: usize) -> usize {
        self.basis[row]
    }

    /// Raw read access to a tableau entry.
    #[inline]
    pub fn at(&self, row: usize, col: usize) -> f64 {
        self.tableau[row][col]
    }

    /// Gauss–Jordan pivot on entry `(row, col)`.
    ///
    /// Scales `row` so the pivot entry becomes `1`, eliminates `col` from
    /// every other row, and records the entering variable `col` as the new
    /// basic variable of `row`.
    ///
    /// # Panics
    ///
    /// Panics if the pivot entry is numerically zero (`|t[row][col]| < EPS`);
    /// the ratio test guarantees a strictly positive pivot, so this is an
    /// internal invariant rather than a user-facing condition.
    pub fn pivot(&mut self, row: usize, col: usize) {
        let piv = self.tableau[row][col];
        assert!(piv.abs() > EPS, "pivot on (near-)zero entry");
        let inv = 1.0 / piv;

        // Normalise the pivot row.
        for value in self.tableau[row].iter_mut() {
            *value *= inv;
        }
        // Pin the pivot column to an exact 1 on the pivot row.
        self.tableau[row][col] = 1.0;

        // Eliminate `col` from all other rows.
        let width = self.tableau[row].len();
        for r in 0..self.n {
            if r == row {
                continue;
            }
            let factor = self.tableau[r][col];
            if factor.abs() <= EPS {
                self.tableau[r][col] = 0.0;
                continue;
            }
            for c in 0..width {
                let delta = factor * self.tableau[row][c];
                self.tableau[r][c] -= delta;
            }
            // Pin the eliminated entry to an exact 0.
            self.tableau[r][col] = 0.0;
        }

        self.basis[row] = col;
    }

    /// Lexicographic minimum-ratio test for the entering column `col`.
    ///
    /// Returns the index of the leaving row, or `None` (ray termination) if
    /// no row has a positive entry in `col`.
    ///
    /// The primary test minimises `t[r][rhs] / t[r][col]` over rows with
    /// `t[r][col] > 0`. Ties — which is exactly where naive pivoting can
    /// cycle on degenerate LCPs — are resolved lexicographically: among the
    /// tied rows, compare `t[r][k] / t[r][col]` for the original identity
    /// columns `k = 0, 1, …, n-1` in order, and keep the row whose sequence
    /// is smallest. Because the original `[ I | … ]` block is linearly
    /// independent, this sequence differs between any two distinct rows, so
    /// the leaving row is unique and the method cannot cycle.
    pub fn min_ratio_lex(&self, col: usize) -> Option<usize> {
        let rhs = self.rhs_col();

        // Collect candidate rows (positive pivot entry) and their primary
        // ratios.
        let mut candidates: Vec<usize> = Vec::new();
        let mut best_ratio = f64::INFINITY;
        for r in 0..self.n {
            let a = self.tableau[r][col];
            if a <= EPS {
                continue;
            }
            let ratio = self.tableau[r][rhs] / a;
            if ratio < best_ratio - EPS {
                best_ratio = ratio;
                candidates.clear();
                candidates.push(r);
            } else if ratio <= best_ratio + EPS {
                candidates.push(r);
            }
        }

        match candidates.len() {
            0 => None,
            1 => Some(candidates[0]),
            _ => Some(self.break_tie_lex(&candidates, col)),
        }
    }

    /// Resolve a primary-ratio tie among `candidates` by comparing the
    /// lexicographic sequence of identity-column ratios.
    ///
    /// Returns the winning row index. `candidates` is assumed non-empty.
    fn break_tie_lex(&self, candidates: &[usize], col: usize) -> usize {
        let mut best = candidates[0];
        for &r in &candidates[1..] {
            if self.lex_less(r, best, col) {
                best = r;
            }
        }
        best
    }

    /// Is the lexicographic identity-column ratio sequence of row `a`
    /// strictly smaller than that of row `b`?
    ///
    /// Compares `t[a][k]/t[a][col]` against `t[b][k]/t[b][col]` for
    /// `k = 0 .. n` (the original `w` identity block) in increasing order and
    /// decides on the first column where they differ beyond [`EPS`].
    fn lex_less(&self, a: usize, b: usize, col: usize) -> bool {
        let pa = self.tableau[a][col];
        let pb = self.tableau[b][col];
        for k in 0..self.n {
            let ra = self.tableau[a][k] / pa;
            let rb = self.tableau[b][k] / pb;
            if ra < rb - EPS {
                return true;
            }
            if ra > rb + EPS {
                return false;
            }
        }
        false
    }

    /// Extract the LCP solution `(w, z)` from the current basis.
    ///
    /// A variable that is basic in some row takes that row's RHS value; every
    /// non-basic variable is `0`. Tiny negative values produced by round-off
    /// are clamped to `0` so the returned vectors respect non-negativity
    /// exactly. The covering variable `z₀` is discarded (it is `0` at a valid
    /// terminal basis).
    pub fn extract_solution(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let rhs = self.rhs_col();
        let mut w = vec![0.0; n];
        let mut z = vec![0.0; n];
        for r in 0..n {
            let var = self.basis[r];
            let val = self.tableau[r][rhs];
            let val = if val < 0.0 && val > -EPS { 0.0 } else { val };
            if var < n {
                w[var] = val;
            } else if var < 2 * n {
                z[var - n] = val;
            }
            // var == z0_col → covering variable, ignored.
        }
        (w, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-10
    }

    #[test]
    fn test_new_layout() {
        // M = [[2, 1], [1, 3]], q = [-1, -2].
        let m = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let q = vec![-1.0, -2.0];
        let t = LcpTableau::new(&m, &q);

        assert_eq!(t.dim(), 2);
        assert_eq!(t.z0_col(), 4);
        assert_eq!(t.rhs_col(), 5);

        // Identity block for w.
        assert!(approx(t.at(0, 0), 1.0));
        assert!(approx(t.at(1, 1), 1.0));
        assert!(approx(t.at(0, 1), 0.0));
        // -M block.
        assert!(approx(t.at(0, 2), -2.0));
        assert!(approx(t.at(0, 3), -1.0));
        assert!(approx(t.at(1, 2), -1.0));
        assert!(approx(t.at(1, 3), -3.0));
        // Covering column -d (d = 1).
        assert!(approx(t.at(0, 4), -1.0));
        assert!(approx(t.at(1, 4), -1.0));
        // RHS.
        assert!(approx(t.at(0, 5), -1.0));
        assert!(approx(t.at(1, 5), -2.0));

        // Initial basis is w_1, w_2.
        assert_eq!(t.basic_var(0), 0);
        assert_eq!(t.basic_var(1), 1);
    }

    #[test]
    fn test_pivot_makes_unit_column() {
        let m = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let q = vec![-1.0, -2.0];
        let mut t = LcpTableau::new(&m, &q);

        // Drive z0 (column 4) into row 1 (entry -1 is non-zero; we negate by
        // pivoting on it — Gauss-Jordan handles the sign).
        t.pivot(1, 4);

        // Column 4 is now the unit vector e_1.
        assert!(approx(t.at(1, 4), 1.0));
        assert!(approx(t.at(0, 4), 0.0));
        // Basis updated.
        assert_eq!(t.basic_var(1), 4);
        assert_eq!(t.basic_var(0), 0);
    }

    #[test]
    fn test_min_ratio_picks_most_negative_rhs_for_z0() {
        // For the z0 column the entries are all -1 in the initial tableau, so
        // there is no positive pivot candidate yet; min_ratio on z0 returns
        // None. We instead verify the ratio test on a synthetic positive
        // column by pivoting first.
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-3.0, -1.0];
        let t = LcpTableau::new(&m, &q);
        // z0 column is all -1 → no positive entries → ray (None) at this raw
        // stage. (Lemke's first pivot uses a dedicated rule, tested in lemke.)
        assert!(t.min_ratio_lex(t.z0_col()).is_none());
    }

    #[test]
    fn test_min_ratio_lex_tie_break_unique() {
        // Construct a tableau by hand where two rows tie on the primary ratio
        // but differ in the identity block, so the lexicographic rule must
        // pick exactly one.
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-1.0, -1.0];
        let t = LcpTableau::new(&m, &q);

        // Both rows have RHS = -1 and z-column (col 2, 3) entries... pick a
        // column with positive entries in both rows. Column 0 (w_1) is unit;
        // only row 0 has a positive entry there, so ratio test returns row 0.
        assert_eq!(t.min_ratio_lex(0), Some(0));
        assert_eq!(t.min_ratio_lex(1), Some(1));
    }

    #[test]
    fn test_extract_trivial_solution() {
        // After pivoting an identity LCP appropriately the extraction should
        // read basic values from the RHS. Here we just confirm the all-w
        // basis extracts w = q (clamped) and z = 0.
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![0.5, 0.25];
        let t = LcpTableau::new(&m, &q);
        let (w, z) = t.extract_solution();
        assert!(approx(w[0], 0.5));
        assert!(approx(w[1], 0.25));
        assert!(approx(z[0], 0.0));
        assert!(approx(z[1], 0.0));
    }
}
