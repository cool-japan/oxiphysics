// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dantzig / Cottle principal pivoting for the LCP `w - M z = q`.
//!
//! Where Lemke augments the problem and pivots on an enlarged tableau, the
//! Dantzig method (Cottle & Dantzig, 1968; in the rigid-body form popularised
//! by Baraff, *Fast Contact Force Computation for Nonpenetrating Rigid
//! Bodies*, SIGGRAPH 1994, and ODE's `dSolveLCP`) works **incrementally**.
//! Contacts are introduced one at a time. Between introductions the solver
//! maintains a *principal* solution in which the index set is partitioned
//! into
//!
//! * a **clamped** set `C` with `z_i > 0` and `w_i = 0`, and
//! * a **free** set with `z_i = 0` and `w_i ≥ 0`,
//!
//! and the relation `w = M z + q` is kept exact throughout. To admit a new
//! contact `i` whose `w_i` is negative, `z_i` is increased along the
//! principal direction that holds every clamped contact at `w = 0`; the step
//! length is the largest that keeps the whole partition feasible. Whatever
//! hits its bound first (`w_i` reaching `0`, or some clamped `z_j` reaching
//! `0`) is moved across the partition, and the drive repeats until `w_i = 0`.
//!
//! For the positive-semi-definite `M` produced by contact assembly this
//! terminates in a finite number of pivots and reaches the same solution as
//! Lemke's method (verified against it in the harness tests).
//!
//! The `M_CC` subsystems are solved with a self-contained Gaussian
//! elimination with partial pivoting — again no nalgebra, in keeping with the
//! crate's plain-`f64`-array house style.

use super::LcpError;
use super::dense::EPS;

/// Solve the LCP `w - M z = q` by Dantzig/Cottle principal pivoting.
///
/// Intended for positive-semi-definite `M` (the contact case). Returns
/// `(w, z)` on success.
///
/// # Errors
///
/// * [`LcpError::DimensionMismatch`] if `m` is not `q.len() × q.len()`.
/// * [`LcpError::DantzigInfeasible`] if a contact's principal subproblem
///   cannot be made feasible (a blocked, zero-length drive — which for a PSD
///   contact LCP signals numerical breakdown rather than a real infeasibility).
///
/// # Examples
///
/// ```
/// use oxiphysics_constraints::lcp::dantzig_solve;
///
/// let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
/// let q = vec![-1.0, -1.0];
/// let (w, z) = dantzig_solve(&m, &q).expect("solvable PSD LCP");
/// assert!((z[0] - 1.0 / 3.0).abs() < 1e-10);
/// assert!((z[1] - 1.0 / 3.0).abs() < 1e-10);
/// ```
pub fn dantzig_solve(m: &[Vec<f64>], q: &[f64]) -> Result<(Vec<f64>, Vec<f64>), LcpError> {
    let n = q.len();
    if m.len() != n || m.iter().any(|row| row.len() != n) {
        return Err(LcpError::DimensionMismatch { m: m.len(), n });
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut z = vec![0.0; n];
    // w = M z + q (equivalently w - M z = q). With z = 0, w = q.
    let mut w = q.to_vec();
    // Membership flag for the clamped set C.
    let mut clamped = vec![false; n];

    for i in 0..n {
        // Contact i already feasible: separated (w_i ≥ 0, z_i = 0). Leave it
        // in the free set and move on.
        if w[i] >= -EPS {
            w[i] = w[i].max(0.0);
            continue;
        }
        drive_to_zero(m, &mut w, &mut z, &mut clamped, i)?;
    }

    // Clean up sign noise so the result respects non-negativity exactly.
    for value in w.iter_mut().chain(z.iter_mut()) {
        if *value < 0.0 && *value > -EPS {
            *value = 0.0;
        }
    }
    Ok((w, z))
}

/// Drive `z[i]` up until `w[i]` reaches zero, maintaining feasibility of the
/// clamped set throughout. On return contact `i` is clamped (`w_i = 0`,
/// `z_i ≥ 0`) and `(w, z)` is a principal solution over `C ∪ {i}`.
fn drive_to_zero(
    m: &[Vec<f64>],
    w: &mut [f64],
    z: &mut [f64],
    clamped: &mut [bool],
    i: usize,
) -> Result<(), LcpError> {
    let n = w.len();
    // A generous pivot budget; each pivot changes the clamped set, of which
    // there are finitely many, so this only trips on numerical breakdown.
    let max_pivots = 4 * n * n + 16;

    for _ in 0..max_pivots {
        // delta_z over the clamped set: how clamped forces must move so the
        // clamped rows stay at w = 0 as z_i increases by one unit.
        //   M_CC · delta_zC = -M_Ci
        let c_idx: Vec<usize> = (0..n).filter(|&j| clamped[j]).collect();
        let delta_zc = solve_principal_direction(m, &c_idx, i)?;

        // Assemble the full delta_z (driven variable + clamped responses).
        let mut delta_z = vec![0.0; n];
        delta_z[i] = 1.0;
        for (k, &j) in c_idx.iter().enumerate() {
            delta_z[j] = delta_zc[k];
        }

        // delta_w = M · delta_z. Clamped rows are 0 by construction; the free
        // rows (including i) tell us how their w moves with the step.
        let mut delta_w = vec![0.0; n];
        for (row, dw) in delta_w.iter_mut().enumerate() {
            if clamped[row] {
                continue;
            }
            let mut acc = 0.0;
            for (col, &dz) in delta_z.iter().enumerate() {
                if dz != 0.0 {
                    acc += m[row][col] * dz;
                }
            }
            *dw = acc;
        }

        // Maximum step before something hits a bound.
        let mut step = f64::INFINITY;
        let mut blocking: Option<Blocker> = None;

        // (a) Driven contact i reaches w_i = 0. delta_w[i] should be > 0
        //     (PSD ⇒ raising z_i raises w_i). The step that zeroes it is
        //     -w_i / delta_w[i].
        if delta_w[i] > EPS {
            let s = -w[i] / delta_w[i];
            if s < step {
                step = s;
                blocking = Some(Blocker::Driven);
            }
        }

        // (b) A clamped contact j has z_j driven down to 0 (delta_zC[k] < 0).
        for (k, &j) in c_idx.iter().enumerate() {
            let dz = delta_zc[k];
            if dz < -EPS {
                let s = -z[j] / dz;
                if s < step - EPS {
                    step = s;
                    blocking = Some(Blocker::Clamped(j));
                }
            }
        }

        // (c) An already-introduced free contact j (< i) would have w_j
        //     pushed below 0; it must enter the clamped set. Only contacts
        //     with index < i participate: they have all been made feasible
        //     (w_j ≥ 0) by earlier outer iterations, whereas contacts j > i
        //     are not yet part of the LCP and must not constrain this drive
        //     (their w_j is still arbitrary — typically negative — and would
        //     otherwise inject a spurious negative step). delta_w[j] < 0 with
        //     current w_j ≥ 0 gives a non-negative blocking step.
        for j in 0..i {
            if clamped[j] {
                continue;
            }
            let dw = delta_w[j];
            if dw < -EPS {
                let s = -w[j] / dw;
                if s < step - EPS {
                    step = s;
                    blocking = Some(Blocker::FreeIntoClamp(j));
                }
            }
        }

        match blocking {
            None => {
                // No bound limits the drive: the principal direction is
                // unbounded. For a PSD contact LCP this is a numerical
                // breakdown rather than a genuine ray.
                return Err(LcpError::DantzigInfeasible { i });
            }
            Some(block) => {
                if !step.is_finite() || step < 0.0 {
                    return Err(LcpError::DantzigInfeasible { i });
                }
                // Apply the step.
                apply_step(w, z, &delta_w, &delta_z, step, &c_idx, i);

                match block {
                    Blocker::Driven => {
                        // Contact i is now clamped; the drive is complete.
                        clamped[i] = true;
                        w[i] = 0.0;
                        return Ok(());
                    }
                    Blocker::Clamped(j) => {
                        // z_j hit 0: it leaves the clamped set.
                        clamped[j] = false;
                        z[j] = 0.0;
                    }
                    Blocker::FreeIntoClamp(j) => {
                        // w_j hit 0: contact j joins the clamped set.
                        clamped[j] = true;
                        w[j] = 0.0;
                    }
                }
                // Loop: recompute the direction with the updated partition.
            }
        }
    }

    Err(LcpError::DantzigInfeasible { i })
}

/// What limits the current principal-direction step.
enum Blocker {
    /// The driven contact `i` reaches `w_i = 0` — the drive succeeds.
    Driven,
    /// Clamped contact `j` reaches `z_j = 0` and leaves the clamped set.
    Clamped(usize),
    /// Free contact `j` reaches `w_j = 0` and joins the clamped set.
    FreeIntoClamp(usize),
}

/// Apply a step of length `step` along the principal direction.
fn apply_step(
    w: &mut [f64],
    z: &mut [f64],
    delta_w: &[f64],
    delta_z: &[f64],
    step: f64,
    c_idx: &[usize],
    i: usize,
) {
    z[i] += step * delta_z[i];
    for &j in c_idx {
        z[j] += step * delta_z[j];
    }
    for (row, &dw) in delta_w.iter().enumerate() {
        w[row] += step * dw;
    }
}

/// Solve `M_CC · x = -M_Ci` for the principal direction of the clamped set.
///
/// `c_idx` lists the clamped indices (in increasing order); the returned
/// vector is indexed the same way. Returns an empty vector when `C` is empty
/// (the driven contact moves on its own).
fn solve_principal_direction(
    m: &[Vec<f64>],
    c_idx: &[usize],
    i: usize,
) -> Result<Vec<f64>, LcpError> {
    let k = c_idx.len();
    if k == 0 {
        return Ok(Vec::new());
    }

    // Build the dense M_CC and right-hand side -M_Ci.
    let mut a = vec![vec![0.0; k]; k];
    let mut b = vec![0.0; k];
    for (r, &row_idx) in c_idx.iter().enumerate() {
        for (c, &col_idx) in c_idx.iter().enumerate() {
            a[r][c] = m[row_idx][col_idx];
        }
        b[r] = -m[row_idx][i];
    }

    gaussian_solve(&mut a, &mut b).map_err(|_| LcpError::DantzigInfeasible { i })
}

/// Solve `A x = b` in place by Gaussian elimination with partial pivoting.
///
/// Returns the solution `x`, or `Err(())` if `A` is singular (zero pivot).
fn gaussian_solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Result<Vec<f64>, ()> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot: largest magnitude in/below the diagonal.
        let mut pivot_row = col;
        let mut pivot_mag = a[col][col].abs();
        for r in (col + 1)..n {
            let mag = a[r][col].abs();
            if mag > pivot_mag {
                pivot_mag = mag;
                pivot_row = r;
            }
        }
        if pivot_mag < EPS {
            return Err(());
        }
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        // Eliminate below the pivot.
        let pivot = a[col][col];
        for r in (col + 1)..n {
            let factor = a[r][col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for c in col..n {
                let delta = factor * a[col][c];
                a[r][c] -= delta;
            }
            b[r] -= factor * b[col];
        }
    }

    // Back substitution.
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut acc = b[row];
        for c in (row + 1)..n {
            acc -= a[row][c] * x[c];
        }
        x[row] = acc / a[row][row];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_gaussian_solve_2x2() {
        // [[2, 1], [1, 3]] x = [3, 5] → x = [0.8, 1.4].
        let mut a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = vec![3.0, 5.0];
        let x = gaussian_solve(&mut a, &mut b).expect("nonsingular");
        assert!((x[0] - 0.8).abs() < 1e-12, "x0={}", x[0]);
        assert!((x[1] - 1.4).abs() < 1e-12, "x1={}", x[1]);
    }

    #[test]
    fn test_gaussian_solve_singular() {
        let mut a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let mut b = vec![1.0, 2.0];
        assert!(gaussian_solve(&mut a, &mut b).is_err());
    }

    #[test]
    fn test_dimension_mismatch() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let q = vec![-1.0];
        let err = dantzig_solve(&m, &q).unwrap_err();
        assert_eq!(err, LcpError::DimensionMismatch { m: 2, n: 1 });
    }

    #[test]
    fn test_dantzig_2x2_known_solution() {
        let m = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let q = vec![-1.0, -1.0];
        let (w, z) = dantzig_solve(&m, &q).expect("solvable");
        assert!((z[0] - 1.0 / 3.0).abs() < 1e-10, "z0={}", z[0]);
        assert!((z[1] - 1.0 / 3.0).abs() < 1e-10, "z1={}", z[1]);
        assert!(w[0].abs() < 1e-10 && w[1].abs() < 1e-10);
        assert!(lcp_residual(&m, &q, &w, &z) < 1e-10);
    }

    #[test]
    fn test_dantzig_3x3_identity() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let q = vec![-1.0, -2.0, -3.0];
        let (w, z) = dantzig_solve(&m, &q).expect("solvable");
        assert!((z[0] - 1.0).abs() < 1e-10);
        assert!((z[1] - 2.0).abs() < 1e-10);
        assert!((z[2] - 3.0).abs() < 1e-10);
        for wi in &w {
            assert!(wi.abs() < 1e-10);
        }
    }

    #[test]
    fn test_dantzig_mixed_active_inactive() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let q = vec![-2.0, 3.0, -1.0];
        let (w, z) = dantzig_solve(&m, &q).expect("solvable");
        assert!((z[0] - 2.0).abs() < 1e-10, "z0={}", z[0]);
        assert!(z[1].abs() < 1e-10, "z1={}", z[1]);
        assert!((z[2] - 1.0).abs() < 1e-10, "z2={}", z[2]);
        assert!((w[1] - 3.0).abs() < 1e-10, "w1={}", w[1]);
        assert!(lcp_residual(&m, &q, &w, &z) < 1e-10);
    }

    #[test]
    fn test_dantzig_coupled() {
        // A coupled PSD system where the clamped set genuinely grows.
        let m = vec![
            vec![4.0, 1.0, 1.0],
            vec![1.0, 3.0, 1.0],
            vec![1.0, 1.0, 2.0],
        ];
        let q = vec![-1.0, -1.0, -1.0];
        let (w, z) = dantzig_solve(&m, &q).expect("solvable");
        assert!(lcp_residual(&m, &q, &w, &z) < 1e-10);
        // All contacts active (q strictly negative, M strictly diag-dominant
        // PSD) ⇒ z > 0, w = 0.
        for zi in &z {
            assert!(*zi > 0.0, "expected active contact, z={zi}");
        }
        for wi in &w {
            assert!(wi.abs() < 1e-10);
        }
    }
}
