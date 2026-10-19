// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free block eigensolver (LOBPCG) for the grid Kohn-Sham Hamiltonian.
//!
//! Locally optimal block preconditioned conjugate gradient (Knyazev, SIAM J.
//! Sci. Comput. 23, 517 (2001)): every iteration performs a Rayleigh-Ritz
//! projection of the operator onto span{X, T R, P} — current Ritz vectors,
//! preconditioned residuals and the previous search directions — with the
//! trial basis explicitly orthonormalised for stability.

use oxiphysics_core::numerical_methods::symmetric_eigen_n;

use super::error::QuantumChemistryMdError;

/// Relative norm below which a candidate basis vector is considered linearly
/// dependent on the current basis and dropped.
const DEPENDENCE_TOLERANCE: f64 = 1e-10;
/// The block H·X is recomputed from scratch every this many iterations to
/// stop round-off drift of the recurrence.
const REFRESH_INTERVAL: usize = 10;

/// Converged eigenpairs.
#[derive(Debug, Clone)]
pub(super) struct EigenPairs {
    /// Ritz values of the whole block in ascending order (the first
    /// `n_wanted` are converged eigenvalues).
    pub(super) values: Vec<f64>,
    /// Euclidean-orthonormal Ritz vectors (same order as `values`).
    pub(super) vectors: Vec<Vec<f64>>,
    /// Largest residual norm ‖Hx − λx‖ among the wanted pairs.
    pub(super) max_residual: f64,
    /// Number of LOBPCG iterations performed.
    pub(super) iterations: usize,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Append the components of `candidates` orthogonal to `basis` (two passes of
/// modified Gram-Schmidt), normalised; dependent candidates are dropped.
fn orthonormalize_into(basis: &mut Vec<Vec<f64>>, candidates: Vec<Vec<f64>>) {
    for mut v in candidates {
        let original = dot(&v, &v).sqrt();
        if !(original.is_finite() && original > 0.0) {
            continue;
        }
        for _ in 0..2 {
            for b in basis.iter() {
                let c = dot(b, &v);
                for (vi, bi) in v.iter_mut().zip(b.iter()) {
                    *vi -= c * bi;
                }
            }
        }
        let norm = dot(&v, &v).sqrt();
        if norm > DEPENDENCE_TOLERANCE * original {
            for vi in &mut v {
                *vi /= norm;
            }
            basis.push(v);
        }
    }
}

/// Deterministic pseudo-random vector (64-bit LCG) used to complete a
/// rank-deficient starting block.
fn filler_vector(len: usize, seed: u64) -> Vec<f64> {
    let mut state = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        })
        .collect()
}

/// Linear combinations Σ_a coeffs[j][a] · vectors[a] for every j.
fn combine(vectors: &[Vec<f64>], coeffs: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let len = vectors.first().map(|v| v.len()).unwrap_or(0);
    coeffs
        .iter()
        .map(|c| {
            let mut out = vec![0.0f64; len];
            for (&w, v) in c.iter().zip(vectors.iter()) {
                if w != 0.0 {
                    for (o, x) in out.iter_mut().zip(v.iter()) {
                        *o += w * x;
                    }
                }
            }
            out
        })
        .collect()
}

/// Rayleigh-Ritz on an orthonormal basis: eigen-decomposition of SᵀHS.
fn rayleigh_ritz(basis: &[Vec<f64>], h_basis: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let k = basis.len();
    let mut a = vec![vec![0.0f64; k]; k];
    for i in 0..k {
        for j in i..k {
            let v = 0.5 * (dot(&basis[i], &h_basis[j]) + dot(&basis[j], &h_basis[i]));
            a[i][j] = v;
            a[j][i] = v;
        }
    }
    symmetric_eigen_n(&a)
}

/// Lowest `n_wanted` eigenpairs of the symmetric operator `apply_h` by LOBPCG.
///
/// * `apply_h` maps a block of vectors to the block of H·v.
/// * `precondition(r, λ)` approximates (H − λ)⁻¹ r with a positive-definite operator.
/// * `initial` is the starting block; its size m ≥ `n_wanted` is the block size
///   (extra vectors accelerate convergence of the wanted ones).
///
/// Convergence: ‖H x_j − λ_j x_j‖ < `tolerance` for every wanted j, with x_j
/// normalised to unit Euclidean length.
///
/// # Errors
/// [`QuantumChemistryMdError::EigensolverNotConverged`] after `max_iter`
/// iterations, [`QuantumChemistryMdError::DimensionMismatch`] when the block is
/// smaller than `n_wanted` or the vectors are empty.
pub(super) fn lobpcg<H, P>(
    mut apply_h: H,
    mut precondition: P,
    initial: Vec<Vec<f64>>,
    n_wanted: usize,
    tolerance: f64,
    max_iter: usize,
) -> Result<EigenPairs, QuantumChemistryMdError>
where
    H: FnMut(&[Vec<f64>]) -> Vec<Vec<f64>>,
    P: FnMut(&[f64], f64) -> Vec<f64>,
{
    let m = initial.len();
    let len = initial.first().map(|v| v.len()).unwrap_or(0);
    if m < n_wanted || len < m {
        return Err(QuantumChemistryMdError::DimensionMismatch {
            what: "eigensolver block",
            expected: n_wanted.max(1),
            found: m.min(len),
        });
    }
    if n_wanted == 0 {
        return Ok(EigenPairs {
            values: Vec::new(),
            vectors: Vec::new(),
            max_residual: 0.0,
            iterations: 0,
        });
    }
    let mut x: Vec<Vec<f64>> = Vec::with_capacity(m);
    orthonormalize_into(&mut x, initial);
    let mut seed = 1u64;
    while x.len() < m {
        orthonormalize_into(&mut x, vec![filler_vector(len, seed)]);
        seed += 1;
    }
    let hx0 = apply_h(&x);
    let (values0, coeffs0) = rayleigh_ritz(&x, &hx0);
    let mut lambda: Vec<f64> = values0[..m].to_vec();
    let mut hx = combine(&hx0, &coeffs0[..m]);
    x = combine(&x, &coeffs0[..m]);
    let mut p: Vec<Vec<f64>> = Vec::new();
    let mut max_residual = f64::INFINITY;

    for iteration in 0..=max_iter {
        if iteration > 0 && iteration % REFRESH_INTERVAL == 0 {
            hx = apply_h(&x);
        }
        let residuals: Vec<Vec<f64>> = x
            .iter()
            .zip(hx.iter())
            .zip(lambda.iter())
            .map(|((xj, hxj), &l)| hxj.iter().zip(xj.iter()).map(|(h, v)| h - l * v).collect())
            .collect();
        let norms: Vec<f64> = residuals.iter().map(|r| dot(r, r).sqrt()).collect();
        max_residual = norms[..n_wanted].iter().cloned().fold(0.0, f64::max);
        if !max_residual.is_finite() {
            break;
        }
        if max_residual < tolerance {
            return Ok(EigenPairs {
                values: lambda,
                vectors: x,
                max_residual,
                iterations: iteration,
            });
        }
        if iteration == max_iter {
            break;
        }
        let w: Vec<Vec<f64>> = residuals
            .iter()
            .zip(norms.iter())
            .zip(lambda.iter())
            .filter(|&((_, &n), _)| n > 0.0)
            .map(|((r, _), &l)| precondition(r, l))
            .collect();
        let mut basis = x.clone();
        orthonormalize_into(&mut basis, w);
        orthonormalize_into(&mut basis, std::mem::take(&mut p));
        let h_new = apply_h(&basis[m..]);
        let mut h_basis = hx.clone();
        h_basis.extend(h_new);
        let (values, coeffs) = rayleigh_ritz(&basis, &h_basis);
        let selected = &coeffs[..m];
        let x_new = combine(&basis, selected);
        hx = combine(&h_basis, selected);
        let outside: Vec<Vec<f64>> = selected
            .iter()
            .map(|c| {
                let mut c = c.clone();
                for v in c.iter_mut().take(m) {
                    *v = 0.0;
                }
                c
            })
            .collect();
        p = combine(&basis, &outside);
        x = x_new;
        lambda = values[..m].to_vec();
    }
    Err(QuantumChemistryMdError::EigensolverNotConverged {
        residual: max_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::super::real_space_grid::{KineticPreconditioner, UniformGrid};
    use super::*;

    #[test]
    fn test_lobpcg_dense_symmetric_matrix() {
        // Tridiagonal 1-D Laplacian: eigenvalues 2 − 2 cos(kπ/(n+1)).
        let n = 60;
        let apply = |block: &[Vec<f64>]| -> Vec<Vec<f64>> {
            block
                .iter()
                .map(|v| {
                    (0..n)
                        .map(|i| {
                            let mut y = 2.0 * v[i];
                            if i > 0 {
                                y -= v[i - 1];
                            }
                            if i + 1 < n {
                                y -= v[i + 1];
                            }
                            y
                        })
                        .collect()
                })
                .collect()
        };
        let initial: Vec<Vec<f64>> = (0..5).map(|s| filler_vector(n, s + 10)).collect();
        let pairs = lobpcg(apply, |r, _| r.to_vec(), initial, 3, 1e-9, 2000).expect("converged");
        for (k, &value) in pairs.values.iter().take(3).enumerate() {
            let exact = 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
            assert!((value - exact).abs() < 1e-10, "{value} vs {exact}");
        }
        for i in 0..3 {
            for j in 0..3 {
                let d = dot(&pairs.vectors[i], &pairs.vectors[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((d - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_lobpcg_harmonic_oscillator_on_grid() {
        // H = −½∇² + ½r²: eigenvalues 1.5 and 2.5 (three-fold).
        let grid = UniformGrid::centered([0.0; 3], 0.3, [28, 28, 28]).expect("grid");
        let potential: Vec<f64> = grid
            .points()
            .iter()
            .map(|p| 0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]))
            .collect();
        let pre = KineticPreconditioner::new(&grid);
        let apply = |block: &[Vec<f64>]| -> Vec<Vec<f64>> {
            block
                .iter()
                .map(|v| {
                    let lap = grid.laplacian(v);
                    lap.iter()
                        .zip(v.iter())
                        .zip(potential.iter())
                        .map(|((l, x), u)| -0.5 * l + u * x)
                        .collect()
                })
                .collect()
        };
        let initial: Vec<Vec<f64>> = (0..6)
            .map(|s| filler_vector(grid.n_points(), s + 3))
            .collect();
        let pairs = lobpcg(
            apply,
            |r, l| pre.apply(r, (-l).max(0.0) + 0.5),
            initial,
            4,
            1e-6,
            500,
        )
        .expect("converged");
        assert!((pairs.values[0] - 1.5).abs() < 5e-3, "{:?}", pairs.values);
        for v in &pairs.values[1..4] {
            assert!((v - 2.5).abs() < 5e-3, "{:?}", pairs.values);
        }
        assert!(pairs.max_residual < 1e-6);
    }

    #[test]
    fn test_lobpcg_rejects_too_small_block() {
        let r = lobpcg(
            |b: &[Vec<f64>]| b.to_vec(),
            |r, _| r.to_vec(),
            vec![vec![1.0, 0.0]],
            2,
            1e-8,
            10,
        );
        assert!(matches!(
            r,
            Err(QuantumChemistryMdError::DimensionMismatch { .. })
        ));
    }
}
