// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Block inverse iteration with Rayleigh-Ritz extraction for dense symmetric matrices.
//!
//! This is the shared engine behind [`super::subspace_iteration`] (shift `0`) and
//! [`super::shift_invert_power`] (arbitrary shift).  A block of start vectors is repeatedly
//! multiplied by `(A - sigma I)^-1`, orthonormalised, and rotated onto the Ritz basis of the
//! projected matrix `Y^T A Y`.  The eigenvalues of `A` closest to the shift dominate the
//! iteration, and because the whole block is iterated and re-orthonormalised together
//! distinct eigenpairs are found instead of the one nearest the shift being found repeatedly.

use super::functions::{jacobi_eigen_dense, matvec_dense};
use crate::dense_lu::LuFactor;

/// Deterministic SplitMix64 generator for reproducible start vectors.
struct SplitMix64(u64);

impl SplitMix64 {
    /// Next pseudo-random number, uniformly distributed in `[-1, 1)`.
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        // Keep the 53 high-quality bits so the conversion to f64 is exact.
        let unit = (z >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * unit - 1.0
    }
}

/// Euclidean inner product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm.
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// `y += alpha * x`.
fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

/// Remove from `v` its components along the orthonormal vectors `basis`.
///
/// Modified Gram-Schmidt applied twice ("twice is enough"), which keeps the result
/// orthogonal to working precision even when `v` is almost in the span of `basis`.
fn project_out(v: &mut [f64], basis: &[Vec<f64>]) {
    for _pass in 0..2 {
        for b in basis {
            let coeff = dot(v, b);
            axpy(-coeff, b, v);
        }
    }
}

/// Replace `current` by a unit vector orthogonal to the orthonormal vectors `done`.
///
/// Pseudo-random directions are tried first; if none is clearly independent the unit
/// coordinate vector with the largest residual is used.  Because `done.len() < n` such a
/// direction always exists.
fn reseed(current: &mut [f64], done: &[Vec<f64>], rng: &mut SplitMix64) {
    let n = current.len();
    for _attempt in 0..4 {
        let mut candidate: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
        let before = norm(&candidate);
        project_out(&mut candidate, done);
        let after = norm(&candidate);
        if before > 0.0 && after > 1.0e-3 * before {
            for (c, x) in current.iter_mut().zip(candidate.iter()) {
                *c = x / after;
            }
            return;
        }
    }
    let mut best_residual = 0.0_f64;
    let mut best: Vec<f64> = vec![0.0; n];
    for k in 0..n {
        let mut unit = vec![0.0; n];
        unit[k] = 1.0;
        project_out(&mut unit, done);
        let residual = norm(&unit);
        if residual > best_residual {
            best_residual = residual;
            best = unit;
        }
    }
    if best_residual > 0.0 {
        for (c, x) in current.iter_mut().zip(best.iter()) {
            *c = x / best_residual;
        }
    }
}

/// Orthonormalise `block` in place.
///
/// A vector that is numerically dependent on its predecessors is replaced by a fresh
/// direction orthogonal to them, so a block of at most `n` vectors always keeps full rank.
fn orthonormalize_block(block: &mut [Vec<f64>], rng: &mut SplitMix64) {
    for i in 0..block.len() {
        let (done, rest) = block.split_at_mut(i);
        let current = &mut rest[0];
        let initial = norm(current);
        project_out(current, done);
        let remaining = norm(current);
        if initial > 0.0 && remaining.is_finite() && remaining > 1.0e-10 * initial {
            for x in current.iter_mut() {
                *x /= remaining;
            }
        } else {
            reseed(current, done, rng);
        }
    }
}

/// LU-factor `A - shift I`, nudging the shift when it sits on an eigenvalue.
///
/// If `A - shift I` is singular to working precision the shift is moved by a small multiple
/// of `sqrt(eps) * max(|A|, |shift|)`, trying both signs and growing the step, until the
/// factorisation succeeds.  The nudge only changes the *filter* of the inverse iteration:
/// the Rayleigh-Ritz step always uses the unshifted matrix, so eigenvalues are unaffected.
fn factor_shifted(a: &[Vec<f64>], n: usize, shift: f64, scale: f64) -> Option<LuFactor> {
    let attempt = |s: f64| -> Option<LuFactor> {
        let mut shifted: Vec<Vec<f64>> = a.to_vec();
        for (i, row) in shifted.iter_mut().enumerate().take(n) {
            row[i] -= s;
        }
        LuFactor::new(&shifted)
    };
    if let Some(factor) = attempt(shift) {
        return Some(factor);
    }
    let magnitude = scale.max(shift.abs());
    let base = if magnitude > 0.0 {
        magnitude * f64::EPSILON.sqrt()
    } else {
        f64::EPSILON.sqrt()
    };
    let mut step = base;
    for _ in 0..4 {
        for sign in [1.0, -1.0] {
            if let Some(factor) = attempt(shift + sign * step) {
                return Some(factor);
            }
        }
        step *= 10.0;
    }
    None
}

/// Give a unit vector a deterministic sign: its largest-magnitude component is positive.
fn fix_sign(v: &mut [f64]) {
    let mut largest = 0.0_f64;
    let mut sign = 1.0_f64;
    for &x in v.iter() {
        if x.abs() > largest {
            largest = x.abs();
            sign = if x < 0.0 { -1.0 } else { 1.0 };
        }
    }
    if sign < 0.0 {
        for x in v.iter_mut() {
            *x = -*x;
        }
    }
}

/// Block inverse iteration with Rayleigh-Ritz extraction.
///
/// Computes the `num_modes` eigenpairs of the symmetric part of the leading `n x n` block of
/// `a` whose eigenvalues are closest to `shift`, nearest first.  The block iterated has
/// `min(n, p + max(min(p, 8), 4))` vectors (`p = num_modes`); the guard vectors beyond `p`
/// accelerate convergence of the wanted ones.  Iteration stops when every wanted pair has
/// `||A x - theta x|| <= tol * |theta| + floor` (`floor` is a rounding-level absolute term so
/// that zero eigenvalues can converge) or after `max(max_iter, 1)` sweeps, in which case the
/// current Ritz pairs are returned.
///
/// Returns an empty vector when nothing is requested, when `a` is smaller than `n x n`, when it
/// contains non-finite entries, or when no factorisation of the shifted matrix can be formed.
pub(super) fn inverse_subspace_iteration(
    a: &[Vec<f64>],
    n: usize,
    shift: f64,
    num_modes: usize,
    max_iter: usize,
    tol: f64,
) -> Vec<(f64, Vec<f64>)> {
    let p = num_modes.min(n);
    if p == 0 || a.len() < n || a.iter().take(n).any(|row| row.len() < n) {
        return Vec::new();
    }
    // The Rayleigh-Ritz projection below is only valid for a symmetric matrix, so work with
    // the exactly symmetric part (identical to `a` for symmetric input).
    let sym: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| 0.5 * (a[i][j] + a[j][i])).collect())
        .collect();
    if !shift.is_finite() || sym.iter().flatten().any(|x| !x.is_finite()) {
        return Vec::new();
    }
    let scale = sym.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    let Some(factor) = factor_shifted(&sym, n, shift, scale) else {
        return Vec::new();
    };

    let q = (p + p.clamp(4, 8)).min(n);
    let mut rng = SplitMix64(0x5851_F42D_4C95_7F2D);
    let mut block: Vec<Vec<f64>> = (0..q)
        .map(|_| (0..n).map(|_| rng.next_f64()).collect())
        .collect();
    orthonormalize_block(&mut block, &mut rng);

    let tol = tol.max(0.0);
    let floor = 100.0 * f64::EPSILON * n as f64 * scale;
    let mut ritz: Vec<(f64, Vec<f64>)> = Vec::new();

    for _sweep in 0..max_iter.max(1) {
        let mut y: Vec<Vec<f64>> = block.iter().map(|column| factor.solve(column)).collect();
        orthonormalize_block(&mut y, &mut rng);
        let ay: Vec<Vec<f64>> = y
            .iter()
            .map(|column| matvec_dense(&sym, column, n))
            .collect();

        // Projected matrix, symmetrised and scaled to unit magnitude because the Jacobi
        // solver works with an absolute off-diagonal tolerance.
        let mut projected = vec![vec![0.0; q]; q];
        for i in 0..q {
            for j in i..q {
                let value = 0.5 * (dot(&y[i], &ay[j]) + dot(&y[j], &ay[i]));
                projected[i][j] = value;
                projected[j][i] = value;
            }
        }
        let projected_scale = projected
            .iter()
            .flatten()
            .fold(0.0_f64, |m, x| m.max(x.abs()));
        let projected_scale = if projected_scale > 0.0 {
            projected_scale
        } else {
            1.0
        };
        for row in projected.iter_mut() {
            for x in row.iter_mut() {
                *x /= projected_scale;
            }
        }
        let (theta_scaled, rotations) = jacobi_eigen_dense(&projected, q);
        let theta: Vec<f64> = theta_scaled.iter().map(|t| t * projected_scale).collect();

        let mut order: Vec<usize> = (0..q).collect();
        order.sort_by(|&i, &j| {
            (theta[i] - shift)
                .abs()
                .total_cmp(&(theta[j] - shift).abs())
        });

        let mut next_block: Vec<Vec<f64>> = Vec::with_capacity(q);
        let mut converged = true;
        ritz.clear();
        for (rank, &k) in order.iter().enumerate() {
            let mut x = vec![0.0; n];
            let mut ax = vec![0.0; n];
            for (j, &coeff) in rotations[k].iter().enumerate() {
                axpy(coeff, &y[j], &mut x);
                axpy(coeff, &ay[j], &mut ax);
            }
            if rank < p {
                let theta_k = theta[k];
                axpy(-theta_k, &x, &mut ax);
                let residual = norm(&ax);
                if residual.is_nan() || residual > tol * theta_k.abs() + floor {
                    converged = false;
                }
                ritz.push((theta_k, x.clone()));
            }
            next_block.push(x);
        }
        block = next_block;
        if converged {
            break;
        }
    }

    for (_, vector) in ritz.iter_mut() {
        let length = norm(vector);
        if length > 0.0 {
            for x in vector.iter_mut() {
                *x /= length;
            }
        }
        fix_sign(vector);
    }
    ritz
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity(n: usize) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect()
    }

    #[test]
    fn orthonormalize_block_returns_an_orthonormal_set() {
        let mut rng = SplitMix64(7);
        let mut block = vec![
            vec![1.0, 1.0, 0.0, 0.0],
            vec![1.0, 0.0, 1.0, 0.0],
            vec![0.0, 1.0, 1.0, 0.5],
        ];
        orthonormalize_block(&mut block, &mut rng);
        for i in 0..block.len() {
            for j in 0..block.len() {
                let d = dot(&block[i], &block[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((d - expected).abs() < 1e-14, "G[{i}][{j}] = {d}");
            }
        }
    }

    #[test]
    fn orthonormalize_block_replaces_dependent_and_zero_vectors() {
        let mut rng = SplitMix64(11);
        let mut block = vec![
            vec![1.0, 0.0, 0.0],
            vec![2.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];
        orthonormalize_block(&mut block, &mut rng);
        for i in 0..3 {
            for j in 0..3 {
                let d = dot(&block[i], &block[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (d - expected).abs() < 1e-12,
                    "G[{i}][{j}] = {d}; block = {block:?}"
                );
            }
        }
    }

    #[test]
    fn factor_shifted_nudges_a_shift_that_sits_on_an_eigenvalue() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 3.0]];
        // shift = 3 makes A - shift I exactly singular.
        assert!(factor_shifted(&a, 2, 3.0, 3.0).is_some());
        assert!(factor_shifted(&a, 2, 0.5, 3.0).is_some());
        // The zero matrix is singular for shift 0 too, and is still factorable after the nudge.
        let zero = vec![vec![0.0; 3]; 3];
        assert!(factor_shifted(&zero, 3, 0.0, 0.0).is_some());
    }

    #[test]
    fn identity_matrix_gives_unit_eigenvalues_and_orthonormal_vectors() {
        let a = identity(5);
        let pairs = inverse_subspace_iteration(&a, 5, 0.0, 3, 50, 1e-10);
        assert_eq!(pairs.len(), 3);
        for (i, (lambda, vector)) in pairs.iter().enumerate() {
            assert!((lambda - 1.0).abs() < 1e-12, "lambda = {lambda}");
            assert!((norm(vector) - 1.0).abs() < 1e-12);
            for (j, (_, other)) in pairs.iter().enumerate() {
                if i != j {
                    assert!(dot(vector, other).abs() < 1e-10);
                }
            }
        }
    }

    #[test]
    fn invalid_requests_return_nothing() {
        let a = identity(3);
        assert!(inverse_subspace_iteration(&a, 3, 0.0, 0, 10, 1e-8).is_empty());
        assert!(inverse_subspace_iteration(&a, 0, 0.0, 2, 10, 1e-8).is_empty());
        assert!(inverse_subspace_iteration(&a, 4, 0.0, 2, 10, 1e-8).is_empty());
        let bad = vec![vec![f64::NAN, 0.0], vec![0.0, 1.0]];
        assert!(inverse_subspace_iteration(&bad, 2, 0.0, 1, 10, 1e-8).is_empty());
    }
}
