// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense linear-algebra kernels shared by the modal, dynamic and buckling solvers.
//!
//! Matrices are stored row-major in flat slices of length `n * n`.  The kernels
//! are deliberately self-contained and scale-aware: every tolerance is
//! *relative* to the magnitude of the data, so matrices in SI units (stiffness
//! of order `1e9`, masses of order `1e-3`) are handled exactly like unit-sized
//! test matrices.
//!
//! * [`generalized_symmetric_eigen`] solves `A·φ = μ·B·φ` for a symmetric `A`
//!   and a symmetric positive-definite `B` by Cholesky reduction
//!   (`B = L·Lᵀ`, `C = L⁻¹·A·L⁻ᵀ`) followed by cyclic Jacobi sweeps on `C`.
//!   It backs the natural-frequency problem `K·φ = ω²·M·φ` and the linear
//!   buckling problem `Kg·φ = (1/λ)·K·φ`.
//! * [`LuFactors`] is an LU factorisation with partial pivoting whose
//!   singularity test is relative to the matrix scale.  It backs the Newmark /
//!   Wilson steppers (factorisation reuse) and the shifted inverse iteration of
//!   the mode-shape recovery.

use crate::error::{Error, Result};

/// Relative tolerance below which `|aᵢⱼ − aⱼᵢ|` is accepted as round-off.
///
/// The reference scale of a pair is `max(√|aᵢᵢ·aⱼⱼ|, |aᵢⱼ|, |aⱼᵢ|)`, which is
/// invariant under a diagonal scaling of the DOFs (translations versus
/// rotations, different unit systems).
const SYMMETRY_TOLERANCE: f64 = 1.0e-8;

/// Check that `a` is a finite `n × n` matrix.
pub(crate) fn validate_square(name: &str, a: &[f64], n: usize) -> Result<()> {
    let expected = n.checked_mul(n).ok_or_else(|| {
        Error::InvalidInput(format!("{name}: dimension {n} overflows the matrix size"))
    })?;
    if a.len() != expected {
        return Err(Error::InvalidInput(format!(
            "{name} must hold {n} x {n} = {expected} entries, got {}",
            a.len()
        )));
    }
    if let Some(position) = a.iter().position(|value| !value.is_finite()) {
        return Err(Error::InvalidInput(format!(
            "{name} contains a non-finite entry at ({}, {})",
            position / n.max(1),
            position % n.max(1)
        )));
    }
    Ok(())
}

/// Validate `a` as a symmetric matrix and return its exactly symmetric part
/// `(A + Aᵀ)/2`.
///
/// Round-off asymmetry (see [`SYMMETRY_TOLERANCE`]) is averaged away; a larger
/// asymmetry means the caller passed a non-symmetric matrix and is reported as
/// [`Error::InvalidInput`], because the Cholesky and Jacobi kernels below only
/// make sense for symmetric data.
fn symmetrized(name: &str, a: &[f64], n: usize) -> Result<Vec<f64>> {
    validate_square(name, a, n)?;
    let mut sym = a.to_vec();
    for i in 0..n {
        for j in (i + 1)..n {
            let upper = a[i * n + j];
            let lower = a[j * n + i];
            let reference = (a[i * n + i].abs() * a[j * n + j].abs())
                .sqrt()
                .max(upper.abs())
                .max(lower.abs());
            if (upper - lower).abs() > SYMMETRY_TOLERANCE * reference {
                return Err(Error::InvalidInput(format!(
                    "{name} is not symmetric: entries ({i}, {j}) = {upper:e} and ({j}, {i}) = {lower:e}"
                )));
            }
            let mean = 0.5 * (upper + lower);
            sym[i * n + j] = mean;
            sym[j * n + i] = mean;
        }
    }
    Ok(sym)
}

/// Lower Cholesky factor `L` (`A = L·Lᵀ`) of a symmetric positive-definite
/// matrix called `name` in error messages; only the lower triangle of `a` is read.
///
/// # Errors
///
/// [`Error::InvalidInput`] when a pivot is not larger than `n·ε` times the
/// corresponding diagonal entry, i.e. `a` is not (numerically) positive
/// definite — for a mass matrix this typically means a massless DOF.
fn cholesky_lower(name: &str, a: &[f64], n: usize) -> Result<Vec<f64>> {
    let mut l = vec![0.0_f64; n * n];
    let pivot_floor = n as f64 * f64::EPSILON;
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for p in 0..j {
                sum -= l[i * n + p] * l[j * n + p];
            }
            if i == j {
                let diagonal = a[i * n + i];
                // `!(x > y)` also rejects a NaN pivot.
                if !(diagonal > 0.0 && sum > pivot_floor * diagonal) {
                    return Err(Error::InvalidInput(format!(
                        "{name} is not symmetric positive definite (pivot {i} is {sum:e} for \
                         diagonal entry {diagonal:e})"
                    )));
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Ok(l)
}

/// Solve `L·X = B` in place for all columns of `B` (`n × n`, row-major) with a
/// lower-triangular `L`; whole rows are updated so the inner loops are contiguous.
fn lower_solve_rows(l: &[f64], x: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..i {
            let multiplier = l[i * n + j];
            if multiplier != 0.0 {
                for column in 0..n {
                    let update = multiplier * x[j * n + column];
                    x[i * n + column] -= update;
                }
            }
        }
        let diagonal = l[i * n + i];
        for column in 0..n {
            x[i * n + column] /= diagonal;
        }
    }
}

/// Solve `Lᵀ·X = B` in place for all columns of `B` with a lower-triangular `L`.
fn lower_transpose_solve_rows(l: &[f64], x: &mut [f64], n: usize) {
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            let multiplier = l[j * n + i];
            if multiplier != 0.0 {
                for column in 0..n {
                    let update = multiplier * x[j * n + column];
                    x[i * n + column] -= update;
                }
            }
        }
        let diagonal = l[i * n + i];
        for column in 0..n {
            x[i * n + column] /= diagonal;
        }
    }
}

/// Transpose a square row-major matrix in place.
fn transpose_in_place(x: &mut [f64], n: usize) {
    for i in 0..n {
        for j in (i + 1)..n {
            x.swap(i * n + j, j * n + i);
        }
    }
}

/// `C = L⁻¹·A·L⁻ᵀ` for a symmetric `A` and a lower-triangular `L`.
fn congruence_reduce(a: &[f64], l: &[f64], n: usize) -> Vec<f64> {
    // W = L⁻¹·A
    let mut work = a.to_vec();
    lower_solve_rows(l, &mut work, n);
    // C = W·L⁻ᵀ  ⇔  Cᵀ = L⁻¹·Wᵀ
    transpose_in_place(&mut work, n);
    lower_solve_rows(l, &mut work, n);
    transpose_in_place(&mut work, n);
    // Remove the round-off asymmetry of the two triangular solves.
    for i in 0..n {
        for j in (i + 1)..n {
            let mean = 0.5 * (work[i * n + j] + work[j * n + i]);
            work[i * n + j] = mean;
            work[j * n + i] = mean;
        }
    }
    work
}

/// Relative margin within which two component magnitudes count as tied when a
/// sign convention is applied, see [`normalise_sign`].
const SIGN_TIE_TOLERANCE: f64 = 1.0e-9;

/// Flip `shape` so that its first component of (nearly) largest magnitude is
/// positive.
///
/// Components whose magnitude is within a relative [`SIGN_TIE_TOLERANCE`] of the
/// largest one count as tied.  Modes of symmetric structures have exactly such
/// ties (`|φᵢ| = |φⱼ|` for mirrored DOFs), and picking the strictly largest
/// entry would let round-off decide the sign.
pub(crate) fn normalise_sign(shape: &mut [f64]) {
    let largest = shape.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if largest == 0.0 {
        return;
    }
    let threshold = (1.0 - SIGN_TIE_TOLERANCE) * largest;
    let reference = shape.iter().copied().find(|v| v.abs() >= threshold);
    if reference.is_some_and(|v| v < 0.0) {
        for value in shape.iter_mut() {
            *value = -*value;
        }
    }
}

/// Eigen-decomposition of a symmetric matrix.
struct SymmetricEigen {
    /// Eigenvalues in ascending order.
    values: Vec<f64>,
    /// Orthonormal eigenvectors, row-major; column `j` belongs to `values[j]`.
    vectors: Vec<f64>,
}

/// Sum of the absolute values of the strictly upper off-diagonal entries.
fn off_diagonal_abs_sum(mat: &[f64], n: usize) -> f64 {
    let mut sum = 0.0;
    for p in 0..n {
        for q in (p + 1)..n {
            sum += mat[p * n + q].abs();
        }
    }
    sum
}

/// Cyclic Jacobi eigen-solver for a symmetric matrix (consumed in place).
///
/// Rutishauser's threshold strategy: the first sweeps only rotate the larger
/// off-diagonal entries, and from the fifth sweep on an entry that is below
/// round-off relative to *both* of its diagonal neighbours is set to zero.  The
/// iteration therefore terminates with an exactly diagonal matrix, and small
/// eigenvalues keep a high relative accuracy.
///
/// # Errors
///
/// [`Error::NotConverged`] if `max_sweeps` sweeps do not diagonalise the matrix
/// (convergence is quadratic, ten sweeps are typical).
fn jacobi_eigen(mut mat: Vec<f64>, n: usize, max_sweeps: usize) -> Result<SymmetricEigen> {
    let mut vectors = vec![0.0_f64; n * n];
    for i in 0..n {
        vectors[i * n + i] = 1.0;
    }
    let mut off_sum = off_diagonal_abs_sum(&mat, n);
    let mut sweeps = 0_usize;
    while off_sum != 0.0 && sweeps < max_sweeps {
        let threshold = if sweeps < 3 {
            0.2 * off_sum / (n * n) as f64
        } else {
            0.0
        };
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = mat[p * n + q];
                let app = mat[p * n + p];
                let aqq = mat[q * n + q];
                let scaled = 100.0 * apq.abs();
                if sweeps >= 4 && app.abs() + scaled == app.abs() && aqq.abs() + scaled == aqq.abs()
                {
                    mat[p * n + q] = 0.0;
                    mat[q * n + p] = 0.0;
                    continue;
                }
                if apq.abs() <= threshold {
                    continue;
                }
                let tau = (aqq - app) / (2.0 * apq);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    1.0 / (tau - (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                mat[p * n + p] = app - t * apq;
                mat[q * n + q] = aqq + t * apq;
                mat[p * n + q] = 0.0;
                mat[q * n + p] = 0.0;
                for r in 0..n {
                    if r == p || r == q {
                        continue;
                    }
                    let arp = mat[r * n + p];
                    let arq = mat[r * n + q];
                    let new_arp = c * arp - s * arq;
                    let new_arq = s * arp + c * arq;
                    mat[r * n + p] = new_arp;
                    mat[p * n + r] = new_arp;
                    mat[r * n + q] = new_arq;
                    mat[q * n + r] = new_arq;
                }
                for r in 0..n {
                    let vrp = vectors[r * n + p];
                    let vrq = vectors[r * n + q];
                    vectors[r * n + p] = c * vrp - s * vrq;
                    vectors[r * n + q] = s * vrp + c * vrq;
                }
            }
        }
        sweeps += 1;
        off_sum = off_diagonal_abs_sum(&mat, n);
    }
    if off_sum != 0.0 {
        let scale = mat.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        return Err(Error::NotConverged {
            solver: "Jacobi symmetric eigensolver".to_string(),
            iterations: sweeps,
            residual: if scale > 0.0 {
                off_sum / scale
            } else {
                off_sum
            },
        });
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| mat[a * n + a].total_cmp(&mat[b * n + b]));
    let values: Vec<f64> = order.iter().map(|&j| mat[j * n + j]).collect();
    let mut sorted = vec![0.0_f64; n * n];
    for (new_col, &old_col) in order.iter().enumerate() {
        for row in 0..n {
            sorted[row * n + new_col] = vectors[row * n + old_col];
        }
    }
    Ok(SymmetricEigen {
        values,
        vectors: sorted,
    })
}

/// Eigenpairs of the generalized symmetric-definite problem `A·φ = μ·B·φ`.
pub(crate) struct GeneralizedEigen {
    /// Eigenvalues `μ` in ascending order.
    pub(crate) values: Vec<f64>,
    /// `B`-orthonormal eigenvectors (`φᵢᵀ·B·φⱼ = δᵢⱼ`), row-major; column `j`
    /// belongs to `values[j]`.
    pub(crate) vectors: Vec<f64>,
}

impl GeneralizedEigen {
    /// Eigenvector of the `index`-th (ascending) eigenvalue.
    pub(crate) fn eigenvector(&self, index: usize) -> Vec<f64> {
        let n = self.values.len();
        (0..n).map(|row| self.vectors[row * n + index]).collect()
    }
}

/// Solve `A·φ = μ·B·φ` for a symmetric `A` and a symmetric positive-definite `B`.
///
/// `B = L·Lᵀ` (Cholesky) reduces the problem to the standard symmetric problem
/// `C·y = μ·y` with `C = L⁻¹·A·L⁻ᵀ`, solved by cyclic Jacobi sweeps; the
/// eigenvectors are recovered as `φ = L⁻ᵀ·y`, which makes them `B`-orthonormal.
/// `max_sweeps` caps the number of Jacobi sweeps; `a_name` / `b_name` label the
/// matrices in error messages.
///
/// # Errors
///
/// * [`Error::InvalidInput`] — a wrongly sized, non-finite or non-symmetric
///   matrix, or a `B` that is not positive definite;
/// * [`Error::NotConverged`] — the Jacobi iteration needs more than `max_sweeps`.
pub(crate) fn generalized_symmetric_eigen(
    (a_name, a): (&str, &[f64]),
    (b_name, b): (&str, &[f64]),
    n: usize,
    max_sweeps: usize,
) -> Result<GeneralizedEigen> {
    let a_sym = symmetrized(a_name, a, n)?;
    let b_sym = symmetrized(b_name, b, n)?;
    if n == 0 {
        return Ok(GeneralizedEigen {
            values: Vec::new(),
            vectors: Vec::new(),
        });
    }
    let l = cholesky_lower(b_name, &b_sym, n)?;
    let reduced = congruence_reduce(&a_sym, &l, n);
    let eigen = jacobi_eigen(reduced, n, max_sweeps)?;
    // φ = L⁻ᵀ·y for every eigenvector column y.
    let mut vectors = eigen.vectors;
    lower_transpose_solve_rows(&l, &mut vectors, n);
    Ok(GeneralizedEigen {
        values: eigen.values,
        vectors,
    })
}

/// LU factorisation `P·A = L·U` with partial (row) pivoting.
pub(crate) struct LuFactors {
    n: usize,
    /// Packed factors: strict lower triangle = `L` (unit diagonal implicit),
    /// upper triangle including the diagonal = `U`.
    lu: Vec<f64>,
    /// `perm[i]` is the row of `A` that ended up in row `i` of `P·A`.
    perm: Vec<usize>,
}

impl LuFactors {
    /// Factorise the `n × n` matrix `a`.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidInput`] — wrong size or a non-finite entry;
    /// * [`Error::SingularSystem`] — a pivot is not larger than `n·ε·max|aᵢⱼ|`,
    ///   i.e. `a` is singular to working precision.  The threshold is relative to
    ///   the matrix scale, so it behaves identically for small and large entries.
    pub(crate) fn factor(a: &[f64], n: usize) -> Result<Self> {
        validate_square("matrix to factorise", a, n)?;
        let mut lu = a.to_vec();
        let scale = lu.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let pivot_tolerance = n as f64 * f64::EPSILON * scale;
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut pivot_row = k;
            let mut pivot_abs = lu[k * n + k].abs();
            for i in (k + 1)..n {
                let candidate = lu[i * n + k].abs();
                if candidate > pivot_abs {
                    pivot_abs = candidate;
                    pivot_row = i;
                }
            }
            if pivot_abs <= pivot_tolerance {
                return Err(Error::SingularSystem);
            }
            if pivot_row != k {
                for j in 0..n {
                    lu.swap(k * n + j, pivot_row * n + j);
                }
                perm.swap(k, pivot_row);
            }
            let pivot = lu[k * n + k];
            for i in (k + 1)..n {
                let multiplier = lu[i * n + k] / pivot;
                lu[i * n + k] = multiplier;
                if multiplier != 0.0 {
                    for j in (k + 1)..n {
                        let update = multiplier * lu[k * n + j];
                        lu[i * n + j] -= update;
                    }
                }
            }
        }
        Ok(Self { n, lu, perm })
    }

    /// Solve `A·x = b`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] if `b` does not have `n` entries.
    pub(crate) fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        let n = self.n;
        if b.len() != n {
            return Err(Error::InvalidInput(format!(
                "right-hand side has {} entries but the factorised matrix is {n} x {n}",
                b.len()
            )));
        }
        // Forward substitution: L·y = P·b.
        let mut x: Vec<f64> = self.perm.iter().map(|&row| b[row]).collect();
        for i in 0..n {
            let lower = &self.lu[i * n..i * n + i];
            let correction: f64 = lower.iter().zip(&x[..i]).map(|(l, y)| l * y).sum();
            x[i] -= correction;
        }
        // Back substitution: U·x = y.
        for i in (0..n).rev() {
            let upper = &self.lu[i * n + i + 1..(i + 1) * n];
            let correction: f64 = upper.iter().zip(&x[i + 1..]).map(|(u, v)| u * v).sum();
            x[i] = (x[i] - correction) / self.lu[i * n + i];
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut c = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    c[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }
        c
    }

    fn transpose(a: &[f64], n: usize) -> Vec<f64> {
        let mut t = a.to_vec();
        transpose_in_place(&mut t, n);
        t
    }

    #[test]
    fn cholesky_reproduces_the_matrix() {
        let a = [4.0, 2.0, 0.6, 2.0, 5.0, 1.0, 0.6, 1.0, 3.0];
        let l = cholesky_lower("A", &a, 3).expect("positive definite");
        let llt = matmul(&l, &transpose(&l, 3), 3);
        for (x, y) in llt.iter().zip(a.iter()) {
            assert!((x - y).abs() < 1e-14, "{x} vs {y}");
        }
        // The strict upper triangle of L stays zero.
        assert_eq!(l[1], 0.0);
        assert_eq!(l[2], 0.0);
        assert_eq!(l[5], 0.0);
    }

    #[test]
    fn cholesky_rejects_indefinite_and_semidefinite_matrices() {
        // Indefinite.
        let indefinite = [1.0, 2.0, 2.0, 1.0];
        assert!(matches!(
            cholesky_lower("A", &indefinite, 2),
            Err(Error::InvalidInput(_))
        ));
        // Singular (rank one): second pivot is exactly zero.
        let singular = [1.0, 1.0, 1.0, 1.0];
        assert!(matches!(
            cholesky_lower("A", &singular, 2),
            Err(Error::InvalidInput(_))
        ));
        // A zero diagonal entry (massless DOF).
        let massless = [1.0, 0.0, 0.0, 0.0];
        assert!(matches!(
            cholesky_lower("A", &massless, 2),
            Err(Error::InvalidInput(_))
        ));
        // NaN pivot.
        let nan = [f64::NAN, 0.0, 0.0, 1.0];
        assert!(cholesky_lower("A", &nan, 2).is_err());
    }

    #[test]
    fn symmetrized_averages_round_off_and_rejects_real_asymmetry() {
        let noisy = [2.0, 1.0, 1.0 + 1.0e-14, 3.0];
        let sym = symmetrized("m", &noisy, 2).expect("round-off asymmetry is accepted");
        assert_eq!(sym[1], sym[2]);
        assert!((sym[1] - 1.0).abs() < 1e-13);

        let skew = [2.0, 1.0, 0.5, 3.0];
        assert!(matches!(
            symmetrized("m", &skew, 2),
            Err(Error::InvalidInput(_))
        ));
        // Diagonal scaling invariance: a 1e12 : 1 graded matrix is still symmetric.
        let graded = [1.0e12, 3.0e3, 3.0e3, 1.0e-6];
        assert!(symmetrized("m", &graded, 2).is_ok());
        // Wrong size and non-finite entries.
        assert!(symmetrized("m", &[1.0, 2.0, 3.0], 2).is_err());
        assert!(symmetrized("m", &[1.0, f64::NAN, f64::NAN, 1.0], 2).is_err());
    }

    #[test]
    fn sign_convention_is_stable_for_tied_components() {
        // The largest magnitude is unique: that entry is made positive.
        let mut unique = [0.2, -0.9, 0.5];
        normalise_sign(&mut unique);
        assert_eq!(unique, [-0.2, 0.9, -0.5]);
        // A mirrored (antisymmetric) pattern: the tie is broken by position, not
        // by which of the two equal magnitudes round-off made slightly larger.
        let mut first_larger = [0.6, 0.1, -0.6 + 1.0e-13];
        let mut second_larger = [0.6 - 1.0e-13, 0.1, -0.6];
        normalise_sign(&mut first_larger);
        normalise_sign(&mut second_larger);
        assert!(first_larger[0] > 0.0 && second_larger[0] > 0.0);
        // Zero vectors are left alone.
        let mut zero = [0.0, 0.0];
        normalise_sign(&mut zero);
        assert_eq!(zero, [0.0, 0.0]);
    }

    #[test]
    fn jacobi_diagonalises_a_symmetric_matrix_and_sorts_ascending() {
        let a = vec![2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
        let eig = jacobi_eigen(a.clone(), 3, 50).expect("converges");
        let expected = [
            2.0 - std::f64::consts::SQRT_2,
            2.0,
            2.0 + std::f64::consts::SQRT_2,
        ];
        for (got, want) in eig.values.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-14, "{got} vs {want}");
        }
        // A·V = V·Λ and VᵀV = I.
        let av = matmul(&a, &eig.vectors, 3);
        for col in 0..3 {
            for row in 0..3 {
                let lhs = av[row * 3 + col];
                let rhs = eig.vectors[row * 3 + col] * eig.values[col];
                assert!((lhs - rhs).abs() < 1e-14);
            }
        }
        let vtv = matmul(&transpose(&eig.vectors, 3), &eig.vectors, 3);
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((vtv[i * 3 + j] - want).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn jacobi_handles_large_scale_matrices_with_a_relative_tolerance() {
        // Stiffness-like entries of order 1e9 must converge as quickly as unit
        // entries (an absolute 1e-14 threshold would never be reached).
        let scale = 2.0e9;
        let a: Vec<f64> = [3.0, 1.0, 0.0, 1.0, 2.0, 0.5, 0.0, 0.5, 1.0]
            .iter()
            .map(|v| v * scale)
            .collect();
        let eig = jacobi_eigen(a, 3, 30).expect("converges within 30 sweeps");
        let trace: f64 = eig.values.iter().sum();
        assert!((trace - 6.0 * scale).abs() < 1e-6 * scale);
        assert!(eig.values.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn jacobi_reports_non_convergence_instead_of_returning_a_partial_result() {
        let a = vec![2.0, 1.0, 1.0, 2.0];
        // Zero sweeps cannot diagonalise a non-diagonal matrix.
        let result = jacobi_eigen(a, 2, 0);
        assert!(matches!(result, Err(Error::NotConverged { .. })));
    }

    #[test]
    fn jacobi_handles_degenerate_and_trivial_inputs() {
        let identity = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let eig = jacobi_eigen(identity, 3, 10).expect("already diagonal");
        assert_eq!(eig.values, vec![1.0, 1.0, 1.0]);
        let one = jacobi_eigen(vec![-4.5], 1, 5).expect("1x1");
        assert_eq!(one.values, vec![-4.5]);
        assert_eq!(one.vectors, vec![1.0]);
    }

    #[test]
    fn generalized_eigen_matches_a_hand_computed_pencil() {
        // det(K - λM) = (3 - 2λ)(1 - λ) - 1 = 2λ² - 5λ + 2  ⇒  λ ∈ {1/2, 2}.
        let k = [3.0, -1.0, -1.0, 1.0];
        let m = [2.0, 0.0, 0.0, 1.0];
        let eig = generalized_symmetric_eigen(("K", &k), ("M", &m), 2, 50).expect("SPD mass");
        assert!((eig.values[0] - 0.5).abs() < 1e-14);
        assert!((eig.values[1] - 2.0).abs() < 1e-14);
        // M-orthonormality and K·φ = λ·M·φ.
        for j in 0..2 {
            let phi = eig.eigenvector(j);
            let mphi = [m[0] * phi[0] + m[1] * phi[1], m[2] * phi[0] + m[3] * phi[1]];
            let kphi = [k[0] * phi[0] + k[1] * phi[1], k[2] * phi[0] + k[3] * phi[1]];
            let norm = phi[0] * mphi[0] + phi[1] * mphi[1];
            assert!((norm - 1.0).abs() < 1e-14);
            for i in 0..2 {
                assert!((kphi[i] - eig.values[j] * mphi[i]).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn generalized_eigen_rejects_non_positive_definite_second_matrix() {
        let k = [1.0, 0.0, 0.0, 1.0];
        let m = [1.0, 2.0, 2.0, 1.0];
        assert!(matches!(
            generalized_symmetric_eigen(("K", &k), ("M", &m), 2, 50),
            Err(Error::InvalidInput(_))
        ));
        assert!(
            generalized_symmetric_eigen(("K", &[]), ("M", &[]), 0, 10)
                .expect("empty problem")
                .values
                .is_empty()
        );
    }

    #[test]
    fn lu_solves_a_system_that_needs_pivoting() {
        let a = [0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0];
        let lu = LuFactors::factor(&a, 3).expect("regular");
        let x_true = [1.0, -2.0, 0.5];
        let b: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|j| a[i * 3 + j] * x_true[j]).sum())
            .collect();
        let x = lu.solve(&b).expect("size matches");
        for (got, want) in x.iter().zip(x_true.iter()) {
            assert!((got - want).abs() < 1e-14, "{got} vs {want}");
        }
        assert!(matches!(lu.solve(&[1.0]), Err(Error::InvalidInput(_))));
    }

    #[test]
    fn lu_singularity_test_is_relative_to_the_matrix_scale() {
        // Exactly singular at two very different scales.
        assert!(matches!(
            LuFactors::factor(&[1.0, 2.0, 2.0, 4.0], 2),
            Err(Error::SingularSystem)
        ));
        assert!(matches!(
            LuFactors::factor(&[1.0e30, 2.0e30, 2.0e30, 4.0e30], 2),
            Err(Error::SingularSystem)
        ));
        // A regular matrix of tiny entries is accepted.
        assert!(LuFactors::factor(&[1.0e-30, 0.0, 0.0, 2.0e-30], 2).is_ok());
        // The zero matrix and non-finite input.
        assert!(matches!(
            LuFactors::factor(&[0.0; 4], 2),
            Err(Error::SingularSystem)
        ));
        assert!(matches!(
            LuFactors::factor(&[1.0, f64::NAN, 0.0, 1.0], 2),
            Err(Error::InvalidInput(_))
        ));
        assert!(LuFactors::factor(&[], 0).is_ok());
    }

    /// Deterministic pseudo-random numbers in `[-1, 1)` from a 64-bit linear
    /// congruential generator, so the stress tests below need no external RNG
    /// and never change between runs.
    fn lcg_stream(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64 * 2.0 - 1.0
        }
    }

    /// `G·Gᵀ + shift·I` for a pseudo-random `G`: symmetric positive definite.
    fn random_spd(n: usize, shift: f64, next: &mut impl FnMut() -> f64) -> Vec<f64> {
        let g: Vec<f64> = (0..n * n).map(|_| next()).collect();
        let mut a = matmul(&g, &transpose(&g, n), n);
        for i in 0..n {
            a[i * n + i] += shift;
        }
        a
    }

    fn apply(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
            .collect()
    }

    /// `n` eigenpairs that satisfy `K·φ = λ·M·φ` and are `M`-orthonormal form the
    /// complete eigen-decomposition of the pencil; this checks both properties
    /// on pseudo-random positive-definite pencils of several sizes and units.
    fn assert_complete_decomposition(k: &[f64], m: &[f64], n: usize) {
        let eig = generalized_symmetric_eigen(("K", k), ("M", m), n, 100).expect("SPD pencil");
        assert_eq!(eig.values.len(), n);
        assert!(
            eig.values.windows(2).all(|w| w[0] <= w[1]),
            "eigenvalues must be ascending"
        );
        let k_norm = k.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        let shapes: Vec<Vec<f64>> = (0..n).map(|j| eig.eigenvector(j)).collect();
        for (j, phi) in shapes.iter().enumerate() {
            let kphi = apply(k, phi, n);
            let mphi = apply(m, phi, n);
            let phi_max = phi.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
            let residual = kphi
                .iter()
                .zip(&mphi)
                .map(|(a, b)| (a - eig.values[j] * b).abs())
                .fold(0.0_f64, f64::max);
            assert!(
                residual <= 1e-10 * k_norm * phi_max,
                "n = {n}, mode {j}: residual {residual:e} (λ = {:e})",
                eig.values[j]
            );
            for (i, other) in shapes.iter().enumerate() {
                let inner: f64 = other.iter().zip(&mphi).map(|(a, b)| a * b).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (inner - expected).abs() < 1e-10,
                    "n = {n}: φ{i}ᵀMφ{j} = {inner:e}"
                );
            }
        }
    }

    #[test]
    fn generalized_eigen_of_random_pencils_is_a_complete_m_orthonormal_decomposition() {
        // (n, scale of K, scale of M, seed): unit sizes, SI-like units (stiffness
        // 1e9, mass 1e-3) and the reverse extreme.
        for (n, k_scale, m_scale, seed) in [
            (1, 1.0, 1.0, 7_u64),
            (5, 1.0, 1.0, 11),
            (17, 1.0, 1.0, 13),
            (24, 3.0e9, 2.0e-3, 17),
            (31, 1.0e-6, 5.0e3, 19),
        ] {
            let mut next = lcg_stream(seed);
            let k: Vec<f64> = random_spd(n, 0.5, &mut next)
                .iter()
                .map(|v| v * k_scale)
                .collect();
            let m: Vec<f64> = random_spd(n, n as f64, &mut next)
                .iter()
                .map(|v| v * m_scale)
                .collect();
            assert_complete_decomposition(&k, &m, n);
        }
    }

    /// `K = 2·M` has the single eigenvalue 2 with multiplicity `n`: the
    /// eigenvectors are an arbitrary `M`-orthonormal basis, but they must still
    /// be a complete and orthonormal set.
    #[test]
    fn generalized_eigen_handles_a_fully_degenerate_pencil() {
        let n = 12;
        let mut next = lcg_stream(23);
        let m = random_spd(n, n as f64, &mut next);
        let k: Vec<f64> = m.iter().map(|v| 2.0 * v).collect();
        assert_complete_decomposition(&k, &m, n);
        let eig = generalized_symmetric_eigen(("K", &k), ("M", &m), n, 100).expect("SPD pencil");
        for value in &eig.values {
            assert!((value - 2.0).abs() < 1e-12, "{value}");
        }
    }

    #[test]
    fn lu_solves_random_diagonally_dominant_systems_to_round_off() {
        for (n, seed) in [(3, 29_u64), (20, 31), (64, 37)] {
            let mut next = lcg_stream(seed);
            let mut a: Vec<f64> = (0..n * n).map(|_| next()).collect();
            for i in 0..n {
                a[i * n + i] += n as f64 * if i % 2 == 0 { 1.0 } else { -1.0 };
            }
            let x_true: Vec<f64> = (0..n).map(|_| next()).collect();
            let b = apply(&a, &x_true, n);
            let x = LuFactors::factor(&a, n)
                .expect("regular")
                .solve(&b)
                .expect("size matches");
            for (got, want) in x.iter().zip(&x_true) {
                assert!((got - want).abs() < 1e-12, "n = {n}: {got} vs {want}");
            }
        }
    }
}
