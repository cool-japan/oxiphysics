// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eigen-decomposition of small dense real symmetric matrices.
//!
//! The I/O crate needs a robust eigen-solver in a few places where the matrix
//! is tiny and fixed-size: the 4×4 key matrix of Horn's quaternion method used
//! to superimpose trajectory frames, and the 3×3 normal matrix of the
//! least-squares gradient estimator used for unstructured VTU grids.  A cyclic
//! Jacobi iteration is ideal for this: it is unconditionally convergent for
//! symmetric input, delivers eigenvalues and orthonormal eigenvectors to full
//! working precision, copes with repeated and zero eigenvalues, and needs no
//! heap allocation.

/// Maximum number of full Jacobi sweeps.  The convergence of the cyclic Jacobi
/// method is quadratic, so well below ten sweeps are needed in practice for the
/// matrix sizes handled here; the cap only guards against pathological input
/// such as non-finite entries.
const MAX_SWEEPS: usize = 64;

/// Result of [`symmetric_eigen`]: `A = V · diag(values) · Vᵀ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SymmetricEigen<const N: usize> {
    /// Eigenvalues sorted in descending order.
    pub(crate) values: [f64; N],
    /// `vectors[k]` is the unit-length eigenvector that belongs to
    /// `values[k]`; the eigenvectors are mutually orthogonal.
    pub(crate) vectors: [[f64; N]; N],
}

/// Compute all eigenvalues and eigenvectors of the real symmetric matrix `a`
/// with the cyclic Jacobi method.
///
/// Only the symmetric part `(a + aᵀ) / 2` of the input is used.  Eigenvalues
/// are returned in descending order; ties keep the order produced by the
/// iteration, so a matrix that is already diagonal returns its unit vectors in
/// the original order for equal eigenvalues.
pub(crate) fn symmetric_eigen<const N: usize>(a: [[f64; N]; N]) -> SymmetricEigen<N> {
    let mut m = a;
    for (p, q) in upper_pairs::<N>() {
        let sym = 0.5 * (m[p][q] + m[q][p]);
        m[p][q] = sym;
        m[q][p] = sym;
    }

    // Accumulated rotation, columns are the eigenvectors.
    let mut v = [[0.0_f64; N]; N];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for sweep in 0..MAX_SWEEPS {
        let off_sum: f64 = upper_pairs::<N>().map(|(p, q)| m[p][q].abs()).sum();
        if off_sum == 0.0 || !off_sum.is_finite() {
            break;
        }
        // Skip tiny rotations during the first sweeps (classical threshold
        // strategy); from the fourth sweep on every non-negligible element is
        // annihilated.
        let threshold = if sweep < 3 {
            0.2 * off_sum / (N * N) as f64
        } else {
            0.0
        };

        for (p, q) in upper_pairs::<N>() {
            let apq = m[p][q];
            let guard = 100.0 * apq.abs();
            if sweep > 3
                && m[p][p].abs() + guard == m[p][p].abs()
                && m[q][q].abs() + guard == m[q][q].abs()
            {
                // Element is below the working precision of the diagonal.
                m[p][q] = 0.0;
                m[q][p] = 0.0;
                continue;
            }
            if apq.abs() <= threshold || apq == 0.0 {
                continue;
            }
            rotate(&mut m, &mut v, p, q);
        }
    }

    let mut order: [usize; N] = [0; N];
    for (k, slot) in order.iter_mut().enumerate() {
        *slot = k;
    }
    // Stable insertion sort by descending eigenvalue (N is tiny).
    for i in 1..N {
        let mut j = i;
        while j > 0 && m[order[j]][order[j]] > m[order[j - 1]][order[j - 1]] {
            order.swap(j, j - 1);
            j -= 1;
        }
    }

    let mut values = [0.0_f64; N];
    let mut vectors = [[0.0_f64; N]; N];
    for (k, &src) in order.iter().enumerate() {
        values[k] = m[src][src];
        for (r, component) in vectors[k].iter_mut().enumerate() {
            *component = v[r][src];
        }
    }
    SymmetricEigen { values, vectors }
}

/// All index pairs `(p, q)` with `p < q < N`, in row-major order.
fn upper_pairs<const N: usize>() -> impl Iterator<Item = (usize, usize)> {
    (0..N).flat_map(|p| ((p + 1)..N).map(move |q| (p, q)))
}

/// Apply one Jacobi rotation in the `(p, q)` plane that zeroes `m[p][q]`,
/// updating the matrix `m` and the accumulated eigenvector matrix `v`.
fn rotate<const N: usize>(m: &mut [[f64; N]; N], v: &mut [[f64; N]; N], p: usize, q: usize) {
    let apq = m[p][q];
    let theta = (m[q][q] - m[p][p]) / (2.0 * apq);
    let sign = if theta < 0.0 { -1.0 } else { 1.0 };
    let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
    let c = 1.0 / (t * t + 1.0).sqrt();
    let s = t * c;

    m[p][p] -= t * apq;
    m[q][q] += t * apq;
    m[p][q] = 0.0;
    m[q][p] = 0.0;
    for r in (0..N).filter(|&r| r != p && r != q) {
        let arp = m[r][p];
        let arq = m[r][q];
        m[r][p] = c * arp - s * arq;
        m[p][r] = m[r][p];
        m[r][q] = s * arp + c * arq;
        m[q][r] = m[r][q];
    }
    for row in v.iter_mut() {
        let vp = row[p];
        let vq = row[q];
        row[p] = c * vp - s * vq;
        row[q] = s * vp + c * vq;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Largest absolute entry of `V · diag(w) · Vᵀ − A`.
    fn reconstruction_error<const N: usize>(a: &[[f64; N]; N], eig: &SymmetricEigen<N>) -> f64 {
        let mut worst = 0.0_f64;
        for (i, row) in a.iter().enumerate() {
            for (j, &entry) in row.iter().enumerate() {
                let sum: f64 = eig
                    .values
                    .iter()
                    .zip(eig.vectors.iter())
                    .map(|(value, vector)| value * vector[i] * vector[j])
                    .sum();
                worst = worst.max((sum - entry).abs());
            }
        }
        worst
    }

    fn orthonormality_error<const N: usize>(eig: &SymmetricEigen<N>) -> f64 {
        let mut worst = 0.0_f64;
        for (i, vi) in eig.vectors.iter().enumerate() {
            for (j, vj) in eig.vectors.iter().enumerate() {
                let dot: f64 = vi.iter().zip(vj.iter()).map(|(x, y)| x * y).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                worst = worst.max((dot - expected).abs());
            }
        }
        worst
    }

    #[test]
    fn two_by_two_known_eigenpairs() {
        let eig = symmetric_eigen([[2.0, 1.0], [1.0, 2.0]]);
        assert!((eig.values[0] - 3.0).abs() < 1e-14);
        assert!((eig.values[1] - 1.0).abs() < 1e-14);
        let s = std::f64::consts::FRAC_1_SQRT_2;
        // Eigenvector of 3 is (1, 1)/sqrt(2) up to sign.
        assert!((eig.vectors[0][0].abs() - s).abs() < 1e-14);
        assert!((eig.vectors[0][1].abs() - s).abs() < 1e-14);
        assert!(eig.vectors[0][0] * eig.vectors[0][1] > 0.0);
        // Eigenvector of 1 is (1, -1)/sqrt(2) up to sign.
        assert!(eig.vectors[1][0] * eig.vectors[1][1] < 0.0);
    }

    #[test]
    fn diagonal_matrix_is_sorted_descending() {
        let eig = symmetric_eigen([[1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 3.0]]);
        assert_eq!(eig.values, [5.0, 3.0, 1.0]);
        assert_eq!(eig.vectors[0], [0.0, 1.0, 0.0]);
        assert_eq!(eig.vectors[1], [0.0, 0.0, 1.0]);
        assert_eq!(eig.vectors[2], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn zero_matrix_returns_identity_basis_in_original_order() {
        let eig = symmetric_eigen([[0.0; 4]; 4]);
        assert_eq!(eig.values, [0.0; 4]);
        for k in 0..4 {
            for r in 0..4 {
                let expected = if r == k { 1.0 } else { 0.0 };
                assert_eq!(eig.vectors[k][r], expected);
            }
        }
    }

    #[test]
    fn repeated_eigenvalues_still_give_an_orthonormal_basis() {
        // diag(2, 2, 7) rotated by a fixed orthogonal matrix keeps the
        // doubled eigenvalue.
        let a = [[3.0, 2.0, 2.0], [2.0, 3.0, 2.0], [2.0, 2.0, 3.0]];
        let eig = symmetric_eigen(a);
        assert!((eig.values[0] - 7.0).abs() < 1e-13);
        assert!((eig.values[1] - 1.0).abs() < 1e-13);
        assert!((eig.values[2] - 1.0).abs() < 1e-13);
        assert!(orthonormality_error(&eig) < 1e-13);
        assert!(reconstruction_error(&a, &eig) < 1e-13);
    }

    #[test]
    fn rank_one_matrix_has_a_zero_eigenvalue_pair() {
        // a = d dᵀ with d = (1, 2, 2): eigenvalues 9, 0, 0.
        let d = [1.0_f64, 2.0, 2.0];
        let mut a = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] = d[i] * d[j];
            }
        }
        let eig = symmetric_eigen(a);
        assert!((eig.values[0] - 9.0).abs() < 1e-13);
        assert!(eig.values[1].abs() < 1e-13);
        assert!(eig.values[2].abs() < 1e-13);
        // The dominant eigenvector is d / |d|.
        for (component, expected) in eig.vectors[0].iter().zip(d.iter()) {
            assert!((component.abs() - expected.abs() / 3.0).abs() < 1e-13);
        }
    }

    #[test]
    fn pseudo_random_four_by_four_reconstructs() {
        // Deterministic LCG so the test needs no external RNG.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64) / ((1_u64 << 53) as f64) * 20.0 - 10.0
        };
        for _ in 0..25 {
            let mut a = [[0.0_f64; 4]; 4];
            for (i, j) in (0..4).flat_map(|i| (i..4).map(move |j| (i, j))) {
                let x = next();
                a[i][j] = x;
                a[j][i] = x;
            }
            let eig = symmetric_eigen(a);
            assert!(reconstruction_error(&a, &eig) < 1e-12, "{a:?}");
            assert!(orthonormality_error(&eig) < 1e-13, "{a:?}");
            for pair in eig.values.windows(2) {
                assert!(pair[0] >= pair[1]);
            }
        }
    }

    #[test]
    fn badly_scaled_matrix_keeps_relative_accuracy() {
        // Entries span twelve orders of magnitude.
        let a = [[1.0e6, 1.0, 0.0], [1.0, 1.0, 1.0e-6], [0.0, 1.0e-6, 1.0e-6]];
        let eig = symmetric_eigen(a);
        assert!(reconstruction_error(&a, &eig) < 1e-9);
        assert!(orthonormality_error(&eig) < 1e-13);
    }

    #[test]
    fn asymmetric_input_uses_the_symmetric_part() {
        let eig = symmetric_eigen([[2.0, 3.0], [-1.0, 2.0]]);
        // Symmetric part is [[2, 1], [1, 2]].
        assert!((eig.values[0] - 3.0).abs() < 1e-14);
        assert!((eig.values[1] - 1.0).abs() < 1e-14);
    }
}
