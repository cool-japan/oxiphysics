// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact small mixed-complementarity solve used by the block contact solver.
//!
//! For a contact manifold of `n ≤ 4` points sharing one normal, the normal
//! impulses `x` must satisfy the linear complementarity problem
//!
//! ```text
//! w = K x + b,   x ≥ 0,   w ≥ 0,   xᵢ wᵢ = 0
//! ```
//!
//! where `K = J M⁻¹ Jᵀ` is the manifold's normal effective-mass matrix and `b`
//! the velocity error. The problem is solved exactly by enumerating active
//! sets (at most 16 for four points), largest first. Four coplanar contact
//! points constrain only three body DOFs, so `K` is rank-deficient there; a
//! tiny Tikhonov term selects the minimum-norm (symmetric) impulse split.

use oxiphysics_core::math::Real;

/// Maximum number of points solved as one block.
pub(super) const MAX_BLOCK: usize = 4;

/// Relative Tikhonov regularisation added to the diagonal of each active-set
/// system (scaled by the largest diagonal entry of `K`).
const REGULARISATION: Real = 1e-9;

/// Feasibility tolerance on impulses and residual velocities.
const FEASIBILITY_TOL: Real = 1e-9;

/// Solve the `n × n` (`n ≤ MAX_BLOCK`) system `a x = rhs` in place by Gaussian
/// elimination with partial pivoting. Returns `None` for a singular matrix.
fn solve_dense(
    mut a: [[Real; MAX_BLOCK]; MAX_BLOCK],
    mut rhs: [Real; MAX_BLOCK],
    n: usize,
) -> Option<[Real; MAX_BLOCK]> {
    let scale = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| a[i][j].abs())
        .fold(0.0, Real::max);
    if scale <= 0.0 {
        return None;
    }
    for col in 0..n {
        let pivot = (col..n).max_by(|&p, &q| {
            a[p][col]
                .abs()
                .partial_cmp(&a[q][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if a[pivot][col].abs() <= 1e-14 * scale {
            return None;
        }
        a.swap(col, pivot);
        rhs.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            if factor == 0.0 {
                continue;
            }
            let pivot_row = a[col];
            for (dst, src) in a[row].iter_mut().zip(pivot_row.iter()).take(n).skip(col) {
                *dst -= factor * src;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut x = [0.0; MAX_BLOCK];
    for row in (0..n).rev() {
        let mut acc = rhs[row];
        for k in (row + 1)..n {
            acc -= a[row][k] * x[k];
        }
        x[row] = acc / a[row][row];
    }
    Some(x)
}

/// Solve the LCP `w = K x + b, x ≥ 0, w ≥ 0, xᵢwᵢ = 0` for `n ≤ MAX_BLOCK`
/// unknowns. Returns `None` when no active set yields a feasible solution
/// (the caller then falls back to per-point Gauss–Seidel).
pub(super) fn solve_lcp(
    k: &[[Real; MAX_BLOCK]; MAX_BLOCK],
    b: &[Real; MAX_BLOCK],
    n: usize,
) -> Option<[Real; MAX_BLOCK]> {
    if n == 0 || n > MAX_BLOCK {
        return None;
    }
    if b.iter().take(n).all(|&bi| bi >= 0.0) {
        return Some([0.0; MAX_BLOCK]);
    }
    let max_diag = (0..n).map(|i| k[i][i].abs()).fold(0.0, Real::max);
    let eps = REGULARISATION * max_diag.max(1e-30);
    // Visit active sets from largest to smallest.
    let mut masks: Vec<u32> = (1u32..(1u32 << n)).collect();
    masks.sort_by_key(|m| std::cmp::Reverse(m.count_ones()));
    for mask in masks {
        let active: Vec<usize> = (0..n).filter(|&i| mask & (1 << i) != 0).collect();
        let m = active.len();
        let mut a = [[0.0; MAX_BLOCK]; MAX_BLOCK];
        let mut rhs = [0.0; MAX_BLOCK];
        for (r, &i) in active.iter().enumerate() {
            for (c, &j) in active.iter().enumerate() {
                a[r][c] = k[i][j];
            }
            a[r][r] += eps;
            rhs[r] = -b[i];
        }
        let Some(sol) = solve_dense(a, rhs, m) else {
            continue;
        };
        if sol
            .iter()
            .take(m)
            .any(|&v| v < -FEASIBILITY_TOL || !v.is_finite())
        {
            continue;
        }
        let mut x = [0.0; MAX_BLOCK];
        for (r, &i) in active.iter().enumerate() {
            x[i] = sol[r].max(0.0);
        }
        let feasible = (0..n).filter(|i| mask & (1 << i) == 0).all(|i| {
            let w: Real = (0..n).map(|j| k[i][j] * x[j]).sum::<Real>() + b[i];
            w >= -FEASIBILITY_TOL * (1.0 + b[i].abs())
        });
        if feasible {
            return Some(x);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(
        k: &[[Real; MAX_BLOCK]; MAX_BLOCK],
        b: &[Real; MAX_BLOCK],
        x: &[Real; MAX_BLOCK],
        n: usize,
    ) -> [Real; MAX_BLOCK] {
        let mut w = [0.0; MAX_BLOCK];
        for i in 0..n {
            w[i] = (0..n).map(|j| k[i][j] * x[j]).sum::<Real>() + b[i];
        }
        w
    }

    #[test]
    fn separating_rhs_gives_zero_impulse() {
        let k = [
            [2.0, 0.5, 0.0, 0.0],
            [0.5, 2.0, 0.0, 0.0],
            [0.0; 4],
            [0.0; 4],
        ];
        let x = solve_lcp(&k, &[1.0, 0.5, 0.0, 0.0], 2).expect("feasible");
        assert_eq!(x, [0.0; 4]);
    }

    #[test]
    fn two_active_points_solve_the_linear_system_exactly() {
        let k = [
            [4.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 0.0, 0.0],
            [0.0; 4],
            [0.0; 4],
        ];
        let b = [-5.0, -4.0, 0.0, 0.0];
        let x = solve_lcp(&k, &b, 2).expect("feasible");
        let w = residual(&k, &b, &x, 2);
        assert!(w[0].abs() < 1e-6 && w[1].abs() < 1e-6, "w = {w:?}");
        assert!(x[0] > 0.0 && x[1] > 0.0);
    }

    #[test]
    fn complementarity_drops_a_pulling_point() {
        // Point 1 would need a negative (pulling) impulse; it must be released.
        let k = [
            [2.0, 1.5, 0.0, 0.0],
            [1.5, 2.0, 0.0, 0.0],
            [0.0; 4],
            [0.0; 4],
        ];
        let b = [-2.0, 1.0, 0.0, 0.0];
        let x = solve_lcp(&k, &b, 2).expect("feasible");
        assert!((x[0] - 1.0).abs() < 1e-6 && x[1] == 0.0, "x = {x:?}");
        let w = residual(&k, &b, &x, 2);
        assert!(w[1] >= 0.0, "released point must separate: w = {w:?}");
    }

    #[test]
    fn rank_deficient_coplanar_square_splits_symmetrically() {
        // Four corner contacts of a unit box (m = 1, I⁻¹ = 6) on the ground:
        // K_ij = 1 + 1.5 (sz_i sz_j + sx_i sx_j), which has rank 3.
        let corners = [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)];
        let mut k = [[0.0; MAX_BLOCK]; MAX_BLOCK];
        for (i, &(xi, zi)) in corners.iter().enumerate() {
            for (j, &(xj, zj)) in corners.iter().enumerate() {
                k[i][j] = 1.0 + 1.5 * (xi * xj + zi * zj);
            }
        }
        let b = [-1.0; MAX_BLOCK];
        let x = solve_lcp(&k, &b, 4).expect("feasible");
        for xi in x {
            assert!((xi - 0.25).abs() < 1e-6, "equal split expected, got {x:?}");
        }
    }
}
