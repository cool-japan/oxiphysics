// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Balanced truncation of stable linear time-invariant systems.
//!
//! For a system `ẋ = A x + B u`, `y = C x` with `A` Hurwitz (all eigenvalues in
//! the open left half-plane) the **controllability** and **observability**
//! Gramians are the unique solutions of the Lyapunov equations
//!
//! ```text
//! A Wc + Wc Aᵀ + B Bᵀ = 0,        Aᵀ Wo + Wo A + Cᵀ C = 0.
//! ```
//!
//! The Hankel singular values `σ₁ ≥ σ₂ ≥ … ≥ σ_N` are the square roots of the
//! eigenvalues of `Wc Wo`; they measure how much each state contributes to the
//! input–output behaviour.  The *square-root* method of Tombs and Postlethwaite
//! (1987) / Laub et al. (1987) never forms the ill-conditioned product `Wc Wo`:
//!
//! 1. factor `Wc = Lc Lcᵀ`, `Wo = Lo Loᵀ`;
//! 2. take the SVD `Loᵀ Lc = U Σ Vᵀ`, so that `Σ` holds the Hankel singular
//!    values;
//! 3. build the biorthogonal projectors
//!    `T_r = Lc V_r Σ_r^(-1/2)` (`N × r`) and `T_l = Σ_r^(-1/2) U_rᵀ Loᵀ`
//!    (`r × N`) with `T_l T_r = I_r`;
//! 4. the reduced system is `A_r = T_l A T_r`, `B_r = T_l B`, `C_r = C T_r`.
//!
//! The reduced model is stable, is itself balanced (both of its Gramians equal
//! `Σ_r`), and satisfies the a-priori error bound
//! `‖G − G_r‖_∞ ≤ 2 (σ_{r+1} + … + σ_N)`.
//!
//! The Lyapunov equations are solved by [`solve_lyapunov`]: a vectorised
//! (Kronecker) pivoted-LU solve for small systems and the squared Smith
//! iteration on the Cayley-transformed equation for larger ones.

use super::{BalancedSystem, DenseMat, mat_mul, mat_transpose};
use crate::analysis::FemError;

/// Largest system order for which the Lyapunov equation is solved directly in
/// its vectorised form `(I⊗A + A⊗I) vec(W) = −vec(Q)` (an `N² × N²` system).
const KRONECKER_MAX_ORDER: usize = 40;
/// Maximum number of squarings of the Smith iteration (`2^64` series terms).
const SMITH_MAX_SQUARINGS: usize = 64;
/// Relative size of a negative Gramian eigenvalue that is still treated as
/// round-off (larger ones prove that `A` is not asymptotically stable).
const PSD_TOLERANCE: f64 = 1.0e-8;
/// Hankel singular values below this fraction of the largest one are treated
/// as zero: the corresponding states are not (numerically) controllable and
/// observable and are never retained.
const HSV_RANK_TOLERANCE: f64 = 1.0e-10;
/// Maximum number of Jacobi sweeps in the eigenvalue and singular value solvers.
const MAX_JACOBI_SWEEPS: usize = 100;

// ---------------------------------------------------------------------------
// Small dense helpers
// ---------------------------------------------------------------------------

/// Check that `m` is a finite `rows × cols` matrix.
fn check_matrix(name: &str, m: &[Vec<f64>], rows: usize, cols: usize) -> Result<(), FemError> {
    if m.len() != rows || m.iter().any(|row| row.len() != cols) {
        return Err(FemError::InvalidInput(format!(
            "{name} must be a {rows}x{cols} matrix"
        )));
    }
    if m.iter().flatten().any(|v| !v.is_finite()) {
        return Err(FemError::InvalidInput(format!(
            "{name} contains a non-finite entry"
        )));
    }
    Ok(())
}

/// Identity matrix of order `n`.
fn identity(n: usize) -> DenseMat {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}

/// Frobenius norm.
fn frobenius(m: &[Vec<f64>]) -> f64 {
    m.iter().flatten().map(|v| v * v).sum::<f64>().sqrt()
}

/// `(M + Mᵀ) / 2` of a square matrix.
fn symmetrized(m: &[Vec<f64>]) -> DenseMat {
    let n = m.len();
    (0..n)
        .map(|i| (0..n).map(|j| 0.5 * (m[i][j] + m[j][i])).collect())
        .collect()
}

/// LU factorisation with partial pivoting of a row-major `n × n` matrix.
struct LuFactors {
    lu: Vec<f64>,
    pivots: Vec<usize>,
    n: usize,
}

impl LuFactors {
    /// Factor `a`; fails with [`FemError::SingularSystem`] if a pivot is
    /// negligible relative to the largest entry of the matrix.
    fn new(mut a: Vec<f64>, n: usize) -> Result<Self, FemError> {
        let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let negligible = scale * n as f64 * f64::EPSILON;
        let mut pivots = vec![0usize; n];
        for col in 0..n {
            let mut pivot_row = col;
            let mut best = a[col * n + col].abs();
            for row in (col + 1)..n {
                let value = a[row * n + col].abs();
                if value > best {
                    best = value;
                    pivot_row = row;
                }
            }
            if best <= negligible {
                return Err(FemError::SingularSystem);
            }
            pivots[col] = pivot_row;
            if pivot_row != col {
                for j in 0..n {
                    a.swap(col * n + j, pivot_row * n + j);
                }
            }
            let pivot = a[col * n + col];
            for row in (col + 1)..n {
                let factor = a[row * n + col] / pivot;
                a[row * n + col] = factor;
                if factor != 0.0 {
                    for j in (col + 1)..n {
                        let update = factor * a[col * n + j];
                        a[row * n + j] -= update;
                    }
                }
            }
        }
        Ok(Self { lu: a, pivots, n })
    }

    /// Solve `A x = b` in place.
    fn solve_in_place(&self, b: &mut [f64]) {
        let n = self.n;
        for (col, &pivot_row) in self.pivots.iter().enumerate() {
            b.swap(col, pivot_row);
        }
        for i in 1..n {
            let row = &self.lu[i * n..i * n + i];
            let sum: f64 = row.iter().zip(b[..i].iter()).map(|(l, x)| l * x).sum();
            b[i] -= sum;
        }
        for i in (0..n).rev() {
            let row = &self.lu[i * n + i + 1..(i + 1) * n];
            let sum: f64 = row.iter().zip(b[i + 1..].iter()).map(|(u, x)| u * x).sum();
            b[i] = (b[i] - sum) / self.lu[i * n + i];
        }
    }

    /// Solve `A X = B` for a dense right-hand-side matrix, column by column.
    fn solve_matrix(&self, b: &[Vec<f64>]) -> DenseMat {
        let n = self.n;
        let cols = b.first().map_or(0, Vec::len);
        let mut x = vec![vec![0.0_f64; cols]; n];
        for j in 0..cols {
            let mut column: Vec<f64> = b.iter().map(|row| row[j]).collect();
            self.solve_in_place(&mut column);
            for (x_row, value) in x.iter_mut().zip(column) {
                x_row[j] = value;
            }
        }
        x
    }
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
///
/// Returns the eigenvalues in descending order and the orthogonal matrix whose
/// **columns** are the matching eigenvectors.  Convergence is judged relative
/// to the Frobenius norm of the matrix.
fn symmetric_eigen(a: &[Vec<f64>]) -> (Vec<f64>, DenseMat) {
    let n = a.len();
    let mut m: DenseMat = a.to_vec();
    let mut v = identity(n);
    for _ in 0..MAX_JACOBI_SWEEPS {
        let total: f64 = m.iter().flatten().map(|x| x * x).sum();
        let off_diagonal: f64 = m
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().enumerate().filter(move |&(j, _)| j != i))
            .map(|(_, x)| x * x)
            .sum();
        if off_diagonal <= f64::EPSILON * f64::EPSILON * total {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p][q];
                if apq == 0.0 {
                    continue;
                }
                let tau = (m[q][q] - m[p][p]) / (2.0 * apq);
                let t = if tau.abs() > 1.0e150 {
                    0.5 / tau
                } else {
                    tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                for row in m.iter_mut() {
                    let (rp, rq) = (row[p], row[q]);
                    row[p] = c * rp - s * rq;
                    row[q] = s * rp + c * rq;
                }
                let row_p = m[p].clone();
                let row_q = m[q].clone();
                for k in 0..n {
                    m[p][k] = c * row_p[k] - s * row_q[k];
                    m[q][k] = s * row_p[k] + c * row_q[k];
                }
                for row in v.iter_mut() {
                    let (vp, vq) = (row[p], row[q]);
                    row[p] = c * vp - s * vq;
                    row[q] = s * vp + c * vq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        m[j][j]
            .partial_cmp(&m[i][i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let eigenvalues: Vec<f64> = order.iter().map(|&i| m[i][i]).collect();
    let vectors: DenseMat = v
        .iter()
        .map(|row| order.iter().map(|&k| row[k]).collect())
        .collect();
    (eigenvalues, vectors)
}

/// Singular value decomposition `M = U diag(σ) Vᵀ` of a square matrix by the
/// one-sided Jacobi (Hestenes) method.
///
/// Returns `(U, σ, V)` with the singular values in descending order and the
/// singular vectors as **columns** of `U` and `V`.  Columns of `U` that belong
/// to a zero singular value are zero.  The method delivers high relative
/// accuracy for the small singular values.
fn svd_jacobi(m: &[Vec<f64>]) -> (DenseMat, Vec<f64>, DenseMat) {
    let n = m.len();
    let mut u_cols: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| m[i][j]).collect()).collect();
    let mut v_cols: Vec<Vec<f64>> = (0..n)
        .map(|j| (0..n).map(|i| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _ in 0..MAX_JACOBI_SWEEPS {
        let mut rotated = false;
        for p in 0..n {
            for q in (p + 1)..n {
                let alpha: f64 = u_cols[p].iter().map(|x| x * x).sum();
                let beta: f64 = u_cols[q].iter().map(|x| x * x).sum();
                let gamma: f64 = u_cols[p]
                    .iter()
                    .zip(u_cols[q].iter())
                    .map(|(x, y)| x * y)
                    .sum();
                if gamma == 0.0 || gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = if zeta.abs() > 1.0e150 {
                    0.5 / zeta
                } else {
                    zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for cols in [&mut u_cols, &mut v_cols] {
                    let col_p = cols[p].clone();
                    let col_q = cols[q].clone();
                    for k in 0..n {
                        cols[p][k] = c * col_p[k] - s * col_q[k];
                        cols[q][k] = s * col_p[k] + c * col_q[k];
                    }
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let sigma_unsorted: Vec<f64> = u_cols
        .iter()
        .map(|col| col.iter().map(|x| x * x).sum::<f64>().sqrt())
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        sigma_unsorted[j]
            .partial_cmp(&sigma_unsorted[i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sigma: Vec<f64> = order.iter().map(|&k| sigma_unsorted[k]).collect();
    let u: DenseMat = (0..n)
        .map(|i| {
            order
                .iter()
                .map(|&k| {
                    if sigma_unsorted[k] > 0.0 {
                        u_cols[k][i] / sigma_unsorted[k]
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let v: DenseMat = (0..n)
        .map(|i| order.iter().map(|&k| v_cols[k][i]).collect())
        .collect();
    (u, sigma, v)
}

// ---------------------------------------------------------------------------
// Lyapunov equation
// ---------------------------------------------------------------------------

/// Solve the continuous Lyapunov equation `A W + W Aᵀ + Q = 0` for `W`.
///
/// For a Hurwitz `A` and symmetric positive semi-definite `Q` the solution is
/// the symmetric positive semi-definite matrix `W = ∫₀^∞ e^{At} Q e^{Aᵀt} dt`.
/// Systems of order up to 40 are solved directly by pivoted LU on the
/// vectorised equation `(I⊗A + A⊗I) vec(W) = −vec(Q)`; larger systems use the squared Smith
/// iteration on the Cayley-transformed (Stein) equation
/// `W = A_p W A_pᵀ + Q_p` with `A_p = (A − pI)⁻¹(A + pI)`,
/// `Q_p = 2p (A − pI)⁻¹ Q (A − pI)⁻ᵀ`, which converges quadratically in the
/// number of squarings whenever `A` is Hurwitz.
///
/// # Errors
/// * [`FemError::InvalidInput`] if the shapes disagree, an entry is not finite
///   or (for the iterative path) `A` is not asymptotically stable.
/// * [`FemError::SingularSystem`] if two eigenvalues of `A` sum to zero (for
///   instance an eigenvalue at the origin or a pair `±λ`), in which case the
///   equation has no unique solution.
/// * [`FemError::SolverNotConverged`] if the Smith iteration does not converge.
pub fn solve_lyapunov(a: &[Vec<f64>], q: &[Vec<f64>]) -> Result<DenseMat, FemError> {
    let n = a.len();
    check_matrix("A", a, n, n)?;
    check_matrix("Q", q, n, n)?;
    if n == 0 {
        return Ok(Vec::new());
    }
    if n <= KRONECKER_MAX_ORDER {
        lyapunov_kronecker(a, q)
    } else {
        lyapunov_smith(a, q)
    }
}

/// Direct solution of `A W + W Aᵀ = −Q` through the vectorised equation.
fn lyapunov_kronecker(a: &[Vec<f64>], q: &[Vec<f64>]) -> Result<DenseMat, FemError> {
    let n = a.len();
    let size = n * n;
    // Unknowns are ordered row-major: index (i, j) -> i * n + j.
    let mut k = vec![0.0_f64; size * size];
    for i in 0..n {
        for j in 0..n {
            let row = i * n + j;
            for l in 0..n {
                // (A W)_ij = Σ_l A_il W_lj
                k[row * size + l * n + j] += a[i][l];
                // (W Aᵀ)_ij = Σ_l W_il A_jl
                k[row * size + i * n + l] += a[j][l];
            }
        }
    }
    let mut rhs: Vec<f64> = q.iter().flat_map(|row| row.iter().map(|v| -v)).collect();
    LuFactors::new(k, size)?.solve_in_place(&mut rhs);
    Ok(rhs.chunks(n).map(<[f64]>::to_vec).collect())
}

/// Squared Smith iteration on the Cayley-transformed Lyapunov equation.
fn lyapunov_smith(a: &[Vec<f64>], q: &[Vec<f64>]) -> Result<DenseMat, FemError> {
    let n = a.len();
    let not_stable = || {
        FemError::InvalidInput("A is not asymptotically stable (Smith iteration diverges)".into())
    };
    // Shift: the root-mean-square singular value of A.  Any positive shift
    // converges for Hurwitz A; this one keeps the spectral radius of the
    // Cayley transform well below one for typical spectra.
    let shift = frobenius(a) / (n as f64).sqrt();
    if shift <= 0.0 {
        return Err(not_stable());
    }
    let mut minus: Vec<f64> = a.iter().flat_map(|row| row.iter().copied()).collect();
    let mut plus = minus.clone();
    for i in 0..n {
        minus[i * n + i] -= shift;
        plus[i * n + i] += shift;
    }
    let lu = LuFactors::new(minus, n)?;
    let plus: DenseMat = plus.chunks(n).map(<[f64]>::to_vec).collect();
    // A_p = (A − pI)⁻¹ (A + pI)
    let mut a_p = lu.solve_matrix(&plus);
    // Q_p = 2p (A − pI)⁻¹ Q (A − pI)⁻ᵀ = 2p ((A − pI)⁻¹ ((A − pI)⁻¹ Q)ᵀ)ᵀ
    let y = lu.solve_matrix(q);
    let z = lu.solve_matrix(&mat_transpose(&y));
    let mut x: DenseMat = mat_transpose(&z)
        .into_iter()
        .map(|row| row.into_iter().map(|v| 2.0 * shift * v).collect())
        .collect();
    for _ in 0..SMITH_MAX_SQUARINGS {
        let increment = mat_mul(&mat_mul(&a_p, &x), &mat_transpose(&a_p));
        let increment_norm = frobenius(&increment);
        if !increment_norm.is_finite() {
            return Err(not_stable());
        }
        for (x_row, inc_row) in x.iter_mut().zip(increment.iter()) {
            for (x_ij, inc_ij) in x_row.iter_mut().zip(inc_row.iter()) {
                *x_ij += inc_ij;
            }
        }
        if increment_norm <= 1.0e-15 * frobenius(&x) {
            return Ok(x);
        }
        a_p = mat_mul(&a_p, &a_p);
    }
    if frobenius(&a_p) >= 1.0 {
        Err(not_stable())
    } else {
        Err(FemError::SolverNotConverged)
    }
}

/// Controllability Gramian `Wc` of `(A, B)`: the solution of
/// `A Wc + Wc Aᵀ + B Bᵀ = 0`.
///
/// # Errors
/// As for [`solve_lyapunov`], plus [`FemError::InvalidInput`] if `B` does not
/// have as many rows as `A`.
pub fn controllability_gramian(a: &[Vec<f64>], b: &[Vec<f64>]) -> Result<DenseMat, FemError> {
    let n = a.len();
    check_matrix("A", a, n, n)?;
    let inputs = b.first().map_or(0, Vec::len);
    check_matrix("B", b, n, inputs)?;
    if n > 0 && inputs == 0 {
        return Ok(vec![vec![0.0; n]; n]);
    }
    let bbt = mat_mul(b, &mat_transpose(b));
    solve_lyapunov(a, &bbt).map(|w| symmetrized(&w))
}

/// Observability Gramian `Wo` of `(A, C)`: the solution of
/// `Aᵀ Wo + Wo A + Cᵀ C = 0`.
///
/// # Errors
/// As for [`solve_lyapunov`], plus [`FemError::InvalidInput`] if `C` does not
/// have as many columns as `A` has rows.
pub fn observability_gramian(a: &[Vec<f64>], c: &[Vec<f64>]) -> Result<DenseMat, FemError> {
    let n = a.len();
    check_matrix("A", a, n, n)?;
    if n > 0 && c.is_empty() {
        return Ok(vec![vec![0.0; n]; n]);
    }
    check_matrix("C", c, c.len(), n)?;
    let ctc = mat_mul(&mat_transpose(c), c);
    solve_lyapunov(&mat_transpose(a), &ctc).map(|w| symmetrized(&w))
}

// ---------------------------------------------------------------------------
// Balancing
// ---------------------------------------------------------------------------

/// Factor `W = L Lᵀ` of a symmetric positive semi-definite Gramian through its
/// eigen-decomposition (`L = V √Λ`), which also copes with singular Gramians
/// of systems that are not fully controllable or observable.
fn psd_factor(w: &[Vec<f64>], label: &str) -> Result<DenseMat, FemError> {
    let (eigenvalues, vectors) = symmetric_eigen(&symmetrized(w));
    let largest = eigenvalues.first().copied().unwrap_or(0.0);
    let smallest = eigenvalues.last().copied().unwrap_or(0.0);
    if smallest < -PSD_TOLERANCE * largest.max(f64::MIN_POSITIVE) {
        return Err(FemError::InvalidInput(format!(
            "the {label} Gramian is not positive semi-definite, so the system is not asymptotically stable"
        )));
    }
    let roots: Vec<f64> = eigenvalues.iter().map(|&l| l.max(0.0).sqrt()).collect();
    Ok(vectors
        .iter()
        .map(|row| row.iter().zip(roots.iter()).map(|(v, r)| v * r).collect())
        .collect())
}

/// Result of the square-root balancing procedure.
#[derive(Debug, Clone)]
pub struct BalancingTransform {
    /// All `N` Hankel singular values `σ₁ ≥ … ≥ σ_N` of the full system.
    pub hankel_singular_values: Vec<f64>,
    /// Left projector `T_l` (`r × N`); `T_l T_r = I_r`.
    pub t_left: DenseMat,
    /// Right projector `T_r` (`N × r`).
    pub t_right: DenseMat,
}

impl BalancingTransform {
    /// Order `r` of the reduced model.
    pub fn reduced_dim(&self) -> usize {
        self.t_left.len()
    }

    /// A-priori bound `2 (σ_{r+1} + … + σ_N)` on the `H∞` norm of the
    /// truncation error `G − G_r`.
    pub fn error_bound(&self) -> f64 {
        2.0 * self
            .hankel_singular_values
            .iter()
            .skip(self.reduced_dim())
            .sum::<f64>()
    }

    /// Project a system onto the retained balanced states:
    /// `(T_l A T_r, T_l B, C T_r)`.
    pub fn reduce(&self, a: &[Vec<f64>], b: &[Vec<f64>], c: &[Vec<f64>]) -> BalancedSystem {
        if self.reduced_dim() == 0 {
            return (Vec::new(), Vec::new(), vec![Vec::new(); c.len()]);
        }
        let a_r = mat_mul(&self.t_left, &mat_mul(a, &self.t_right));
        let b_r = mat_mul(&self.t_left, b);
        let c_r = mat_mul(c, &self.t_right);
        (a_r, b_r, c_r)
    }
}

/// Compute the balancing transformation of the stable system `(A, B, C)` and
/// keep at most `max_modes` states.
///
/// The number of retained states is `min(max_modes, N, rank)`, where the rank
/// counts the Hankel singular values above `1e-10 σ₁`: states that are not
/// simultaneously controllable and observable are never kept.
///
/// # Errors
/// * [`FemError::InvalidInput`] if the shapes disagree or `A` is not
///   asymptotically stable (a Gramian is not positive semi-definite).
/// * [`FemError::SingularSystem`] if the Lyapunov equations have no unique
///   solution (eigenvalues of `A` on the imaginary axis).
/// * [`FemError::SolverNotConverged`] if the iterative Lyapunov solver fails.
pub fn balancing_transform(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    max_modes: usize,
) -> Result<BalancingTransform, FemError> {
    let n = a.len();
    check_matrix("A", a, n, n)?;
    let inputs = b.first().map_or(0, Vec::len);
    check_matrix("B", b, n, inputs)?;
    if n > 0 {
        check_matrix("C", c, c.len(), n)?;
    }
    if n == 0 {
        return Ok(BalancingTransform {
            hankel_singular_values: Vec::new(),
            t_left: Vec::new(),
            t_right: Vec::new(),
        });
    }
    let wc = controllability_gramian(a, b)?;
    let wo = observability_gramian(a, c)?;
    let lc = psd_factor(&wc, "controllability")?;
    let lo = psd_factor(&wo, "observability")?;
    let (u, sigma, v) = svd_jacobi(&mat_mul(&mat_transpose(&lo), &lc));
    let sigma_max = sigma.first().copied().unwrap_or(0.0);
    let rank = sigma
        .iter()
        .take_while(|&&s| sigma_max > 0.0 && s > HSV_RANK_TOLERANCE * sigma_max)
        .count();
    let r = max_modes.min(rank);
    if r == 0 {
        return Ok(BalancingTransform {
            hankel_singular_values: sigma,
            t_left: Vec::new(),
            t_right: vec![Vec::new(); n],
        });
    }
    // T_r = Lc V_r Σ_r^(-1/2),  T_l = Σ_r^(-1/2) U_rᵀ Loᵀ
    let inv_sqrt: Vec<f64> = sigma.iter().take(r).map(|s| 1.0 / s.sqrt()).collect();
    let scale_columns = |vectors: &[Vec<f64>]| -> DenseMat {
        vectors
            .iter()
            .map(|row| {
                row.iter()
                    .zip(inv_sqrt.iter())
                    .map(|(x, w)| x * w)
                    .collect()
            })
            .collect()
    };
    let t_right = mat_mul(&lc, &scale_columns(&v));
    let t_left = mat_mul(&mat_transpose(&scale_columns(&u)), &mat_transpose(&lo));
    Ok(BalancingTransform {
        hankel_singular_values: sigma,
        t_left,
        t_right,
    })
}

/// Balanced truncation of the stable linear time-invariant system `(A, B, C)`
/// to (at most) `n_reduced` states.
///
/// Returns the reduced matrices `(A_r, B_r, C_r)`.  The reduced model is
/// stable and balanced, and its transfer function deviates from the full one by
/// at most `2 (σ_{r+1} + … + σ_N)` in the `H∞` norm (see
/// [`BalancingTransform::error_bound`]).  Empty input, or `n_reduced == 0`,
/// yields empty matrices.
///
/// # Errors
/// As for [`balancing_transform`].
pub fn balanced_truncation(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    n_reduced: usize,
) -> Result<BalancedSystem, FemError> {
    if a.is_empty() || n_reduced == 0 {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    }
    let transform = balancing_transform(a, b, c, n_reduced)?;
    Ok(transform.reduce(a, b, c))
}

/// Hankel singular values `σ₁ ≥ … ≥ σ_N` of the stable system `(A, B, C)`:
/// the square roots of the eigenvalues of `Wc Wo`.
///
/// The system is not reduced; all `N` values are returned (states that are
/// not controllable and observable contribute zeros).
///
/// # Errors
/// As for [`balancing_transform`].
pub fn hankel_singular_values(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
) -> Result<Vec<f64>, FemError> {
    Ok(balancing_transform(a, b, c, 0)?.hankel_singular_values)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random value in (-1, 1).
    fn lcg_unit(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn random_matrix(rows: usize, cols: usize, seed: u64) -> DenseMat {
        let mut state = seed;
        (0..rows)
            .map(|_| (0..cols).map(|_| lcg_unit(&mut state)).collect())
            .collect()
    }

    /// Random Hurwitz matrix: Gershgorin discs all lie left of `-margin`.
    fn random_stable(n: usize, seed: u64, margin: f64) -> DenseMat {
        let mut a = random_matrix(n, n, seed);
        let row_sum = a
            .iter()
            .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max);
        for (i, row) in a.iter_mut().enumerate() {
            row[i] -= row_sum + margin;
        }
        a
    }

    fn lyapunov_residual(a: &[Vec<f64>], w: &[Vec<f64>], q: &[Vec<f64>]) -> f64 {
        let aw = mat_mul(a, w);
        let wat = mat_mul(w, &mat_transpose(a));
        let n = a.len();
        let mut sum = 0.0;
        for i in 0..n {
            for j in 0..n {
                sum += (aw[i][j] + wat[i][j] + q[i][j]).powi(2);
            }
        }
        sum.sqrt()
    }

    fn max_asymmetry(w: &[Vec<f64>]) -> f64 {
        w.iter()
            .enumerate()
            .flat_map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(move |(j, &value)| (value - w[j][i]).abs())
            })
            .fold(0.0_f64, f64::max)
    }

    fn diagonal(values: &[f64]) -> DenseMat {
        let n = values.len();
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { values[i] } else { 0.0 })
                    .collect()
            })
            .collect()
    }

    /// Frequency response `c (jωI − A)⁻¹ b` of a single-input single-output
    /// system, as `(re, im)`, through the real 2N × 2N embedding.
    fn siso_response(a: &[Vec<f64>], b: &[Vec<f64>], c: &[Vec<f64>], omega: f64) -> (f64, f64) {
        let n = a.len();
        let size = 2 * n;
        let mut m = vec![0.0_f64; size * size];
        for i in 0..n {
            for j in 0..n {
                m[i * size + j] = -a[i][j];
                m[(n + i) * size + n + j] = -a[i][j];
            }
            m[i * size + n + i] = -omega;
            m[(n + i) * size + i] = omega;
        }
        let lu = LuFactors::new(m, size).expect("(jωI - A) is regular for a stable A");
        let mut rhs = vec![0.0_f64; size];
        for i in 0..n {
            rhs[i] = b[i][0];
        }
        lu.solve_in_place(&mut rhs);
        let re: f64 = (0..n).map(|i| c[0][i] * rhs[i]).sum();
        let im: f64 = (0..n).map(|i| c[0][i] * rhs[n + i]).sum();
        (re, im)
    }

    fn response_error(
        full: &(DenseMat, DenseMat, DenseMat),
        reduced: &(DenseMat, DenseMat, DenseMat),
        omega: f64,
    ) -> f64 {
        let (fr, fi) = siso_response(&full.0, &full.1, &full.2, omega);
        let (rr, ri) = siso_response(&reduced.0, &reduced.1, &reduced.2, omega);
        ((fr - rr).powi(2) + (fi - ri).powi(2)).sqrt()
    }

    fn two_pole_system() -> (DenseMat, DenseMat, DenseMat) {
        (
            diagonal(&[-1.0, -2.0]),
            vec![vec![1.0], vec![1.0]],
            vec![vec![1.0, 1.0]],
        )
    }

    // ------------------------------------------------------------------
    // Dense helpers
    // ------------------------------------------------------------------

    #[test]
    fn lu_solves_a_pivoting_system() {
        // Requires a row exchange: the first pivot is zero.
        let a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0];
        let lu = LuFactors::new(a.clone(), 3).expect("regular matrix");
        let x_true = [1.0, -2.0, 3.0];
        let mut b: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|j| a[i * 3 + j] * x_true[j]).sum())
            .collect();
        lu.solve_in_place(&mut b);
        for (xi, ti) in b.iter().zip(x_true.iter()) {
            assert!((xi - ti).abs() < 1e-13, "x = {b:?}");
        }
        let singular = vec![1.0, 2.0, 2.0, 4.0];
        assert!(matches!(
            LuFactors::new(singular, 2),
            Err(FemError::SingularSystem)
        ));
    }

    #[test]
    fn symmetric_eigen_reconstructs_the_matrix() {
        let n = 8;
        let r = random_matrix(n, n, 21);
        let a: DenseMat = (0..n)
            .map(|i| (0..n).map(|j| r[i][j] + r[j][i]).collect())
            .collect();
        let (lambda, v) = symmetric_eigen(&a);
        for pair in lambda.windows(2) {
            assert!(pair[0] >= pair[1], "eigenvalues must be descending");
        }
        let vt = mat_transpose(&v);
        let vtv = mat_mul(&vt, &v);
        for (i, row) in vtv.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((value - want).abs() < 1e-12, "VᵀV not identity");
            }
        }
        let scaled: DenseMat = v
            .iter()
            .map(|row| row.iter().zip(lambda.iter()).map(|(x, l)| x * l).collect())
            .collect();
        let rebuilt = mat_mul(&scaled, &vt);
        for i in 0..n {
            for j in 0..n {
                assert!((rebuilt[i][j] - a[i][j]).abs() < 1e-12, "V Λ Vᵀ != A");
            }
        }
    }

    #[test]
    fn svd_jacobi_reconstructs_general_and_rank_deficient_matrices() {
        let n = 7;
        let m = random_matrix(n, n, 5);
        let (u, sigma, v) = svd_jacobi(&m);
        for pair in sigma.windows(2) {
            assert!(pair[0] >= pair[1]);
        }
        let us: DenseMat = u
            .iter()
            .map(|row| row.iter().zip(sigma.iter()).map(|(x, s)| x * s).collect())
            .collect();
        let rebuilt = mat_mul(&us, &mat_transpose(&v));
        for i in 0..n {
            for j in 0..n {
                assert!((rebuilt[i][j] - m[i][j]).abs() < 1e-12, "U Σ Vᵀ != M");
            }
        }
        // Rank two: outer products of two vector pairs.
        let x = random_matrix(5, 2, 8);
        let y = random_matrix(2, 5, 9);
        let low_rank = mat_mul(&x, &y);
        let (_, s, _) = svd_jacobi(&low_rank);
        assert!(s[1] > 1e-3);
        for tail in &s[2..] {
            assert!(*tail < 1e-13, "numerical rank must be two: {s:?}");
        }
    }

    // ------------------------------------------------------------------
    // Lyapunov equation
    // ------------------------------------------------------------------

    #[test]
    fn lyapunov_scalar_equation() {
        let w = solve_lyapunov(&[vec![-2.0]], &[vec![4.0]]).expect("stable scalar");
        assert!((w[0][0] - 1.0).abs() < 1e-14, "w = {}", w[0][0]);
    }

    #[test]
    fn lyapunov_diagonal_matches_closed_form() {
        let poles = [-1.0, -2.0, -3.0];
        let b = [1.0, 2.0, 3.0];
        let a = diagonal(&poles);
        let q: DenseMat = b
            .iter()
            .map(|bi| b.iter().map(|bj| bi * bj).collect())
            .collect();
        let w = solve_lyapunov(&a, &q).expect("stable diagonal system");
        for i in 0..3 {
            for j in 0..3 {
                let exact = b[i] * b[j] / -(poles[i] + poles[j]);
                assert!((w[i][j] - exact).abs() < 1e-13, "W[{i}][{j}] = {}", w[i][j]);
            }
        }
    }

    #[test]
    fn lyapunov_kronecker_solution_satisfies_the_equation() {
        let n = 7;
        let a = random_stable(n, 3, 0.4);
        let r = random_matrix(n, n, 4);
        let q = mat_mul(&r, &mat_transpose(&r));
        let w = solve_lyapunov(&a, &q).expect("stable system");
        assert!(lyapunov_residual(&a, &w, &q) < 1e-11 * frobenius(&q));
        assert!(max_asymmetry(&w) < 1e-12 * frobenius(&w));
    }

    #[test]
    fn lyapunov_supports_a_non_symmetric_right_hand_side() {
        let n = 5;
        let a = random_stable(n, 11, 0.5);
        let q = random_matrix(n, n, 12);
        let w = solve_lyapunov(&a, &q).expect("stable system");
        assert!(lyapunov_residual(&a, &w, &q) < 1e-11 * frobenius(&q));
        assert!(
            max_asymmetry(&w) > 1e-3,
            "solution should follow the asymmetric Q"
        );
    }

    #[test]
    fn lyapunov_smith_agrees_with_the_direct_solution() {
        let n = 9;
        let a = random_stable(n, 31, 0.3);
        let r = random_matrix(n, n, 32);
        let q = mat_mul(&r, &mat_transpose(&r));
        let direct = lyapunov_kronecker(&a, &q).expect("direct solve");
        let iterative = lyapunov_smith(&a, &q).expect("Smith iteration");
        let diff: DenseMat = direct
            .iter()
            .zip(iterative.iter())
            .map(|(d, s)| d.iter().zip(s.iter()).map(|(x, y)| x - y).collect())
            .collect();
        assert!(
            frobenius(&diff) < 1e-10 * frobenius(&direct),
            "relative difference {}",
            frobenius(&diff) / frobenius(&direct)
        );
    }

    #[test]
    fn lyapunov_large_system_uses_the_iteration_and_satisfies_the_equation() {
        let n = KRONECKER_MAX_ORDER + 8;
        let a = random_stable(n, 41, 0.3);
        let r = random_matrix(n, 3, 42);
        let q = mat_mul(&r, &mat_transpose(&r));
        let w = solve_lyapunov(&a, &q).expect("stable system");
        assert!(lyapunov_residual(&a, &w, &q) < 1e-9 * frobenius(&q));
        assert!(max_asymmetry(&w) < 1e-9 * frobenius(&w));
        let (eigenvalues, _) = symmetric_eigen(&symmetrized(&w));
        assert!(
            eigenvalues[n - 1] > -1e-10 * eigenvalues[0],
            "Gramian must be positive semi-definite: {}",
            eigenvalues[n - 1]
        );
    }

    #[test]
    fn lyapunov_reports_missing_unique_solutions() {
        let q = identity(2);
        // λ₁ + λ₂ = 0
        assert!(matches!(
            solve_lyapunov(&diagonal(&[-1.0, 1.0]), &q),
            Err(FemError::SingularSystem)
        ));
        // eigenvalue at the origin
        assert!(matches!(
            solve_lyapunov(&diagonal(&[0.0, -1.0]), &q),
            Err(FemError::SingularSystem)
        ));
        // pure imaginary pair ±i
        let rotation = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        assert!(matches!(
            solve_lyapunov(&rotation, &q),
            Err(FemError::SingularSystem)
        ));
    }

    #[test]
    fn lyapunov_iteration_rejects_unstable_matrices() {
        let n = KRONECKER_MAX_ORDER + 10;
        let growth: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let result = solve_lyapunov(&diagonal(&growth), &identity(n));
        assert!(
            matches!(result, Err(FemError::InvalidInput(_))),
            "unstable A must be rejected, got {result:?}"
        );
    }

    #[test]
    fn lyapunov_validates_shapes_and_values() {
        let a = diagonal(&[-1.0, -2.0]);
        assert!(matches!(
            solve_lyapunov(&a, &identity(3)),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            solve_lyapunov(&[vec![-1.0, 0.0]], &identity(2)),
            Err(FemError::InvalidInput(_))
        ));
        let nan = vec![vec![f64::NAN, 0.0], vec![0.0, 1.0]];
        assert!(matches!(
            solve_lyapunov(&a, &nan),
            Err(FemError::InvalidInput(_))
        ));
        assert_eq!(
            solve_lyapunov(&[], &[]).expect("empty system"),
            Vec::<Vec<f64>>::new()
        );
    }

    // ------------------------------------------------------------------
    // Gramians and Hankel singular values
    // ------------------------------------------------------------------

    #[test]
    fn gramians_of_the_two_pole_system() {
        let (a, b, c) = two_pole_system();
        let expected = [[0.5, 1.0 / 3.0], [1.0 / 3.0, 0.25]];
        let wc = controllability_gramian(&a, &b).expect("stable");
        let wo = observability_gramian(&a, &c).expect("stable");
        for i in 0..2 {
            for j in 0..2 {
                assert!((wc[i][j] - expected[i][j]).abs() < 1e-14);
                assert!((wo[i][j] - expected[i][j]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn hankel_singular_values_of_the_two_pole_system() {
        // Wc = Wo = W, hence sigma = eig(W): trace 3/4, determinant 1/72.
        let (a, b, c) = two_pole_system();
        let hsv = hankel_singular_values(&a, &b, &c).expect("stable");
        let half_trace: f64 = 0.375;
        let root = (half_trace * half_trace - 1.0 / 72.0).sqrt();
        assert_eq!(hsv.len(), 2);
        assert!(
            (hsv[0] - (half_trace + root)).abs() < 1e-13,
            "sigma1 = {}",
            hsv[0]
        );
        assert!(
            (hsv[1] - (half_trace - root)).abs() < 1e-13,
            "sigma2 = {}",
            hsv[1]
        );
        assert!((hsv[0] - 0.7310).abs() < 1e-4);
        assert!((hsv[1] - 0.0190).abs() < 1e-4);
    }

    #[test]
    fn hankel_singular_values_are_a_state_similarity_invariant() {
        let n = 5;
        let a = random_stable(n, 61, 0.5);
        let b = random_matrix(n, 2, 62);
        let c = random_matrix(3, n, 63);
        let original = hankel_singular_values(&a, &b, &c).expect("stable");
        // x = S z with S = I + 0.3 R:  (S⁻¹ A S, S⁻¹ B, C S)
        let mut s = random_matrix(n, n, 64);
        for (i, row) in s.iter_mut().enumerate() {
            for v in row.iter_mut() {
                *v *= 0.3;
            }
            row[i] += 1.0;
        }
        let flat: Vec<f64> = s.iter().flat_map(|row| row.iter().copied()).collect();
        let lu = LuFactors::new(flat, n).expect("S is regular");
        let a_t = lu.solve_matrix(&mat_mul(&a, &s));
        let b_t = lu.solve_matrix(&b);
        let c_t = mat_mul(&c, &s);
        let transformed = hankel_singular_values(&a_t, &b_t, &c_t).expect("stable");
        for (x, y) in original.iter().zip(transformed.iter()) {
            assert!((x - y).abs() < 1e-9 * original[0], "{x} vs {y}");
        }
    }

    #[test]
    fn hankel_singular_values_are_descending_and_non_negative() {
        let n = 6;
        let a = random_stable(n, 71, 0.4);
        let hsv = hankel_singular_values(&a, &random_matrix(n, 1, 72), &random_matrix(1, n, 73))
            .expect("stable");
        assert_eq!(hsv.len(), n);
        for pair in hsv.windows(2) {
            assert!(pair[0] >= pair[1]);
        }
        assert!(hsv.iter().all(|&s| s >= 0.0 && s.is_finite()));
    }

    #[test]
    fn unobservable_mode_has_a_zero_hankel_singular_value() {
        // C only sees the first pole: Wo = diag(1/2, 0) and sigma = (1/2, 0).
        let a = diagonal(&[-1.0, -2.0]);
        let b = vec![vec![1.0], vec![1.0]];
        let c = vec![vec![1.0, 0.0]];
        let hsv = hankel_singular_values(&a, &b, &c).expect("stable");
        assert!((hsv[0] - 0.5).abs() < 1e-13, "sigma1 = {}", hsv[0]);
        assert!(hsv[1].abs() < 1e-13, "sigma2 = {}", hsv[1]);
    }

    // ------------------------------------------------------------------
    // Balanced truncation
    // ------------------------------------------------------------------

    #[test]
    fn truncation_drops_a_non_minimal_state_without_error() {
        let (a, b, _) = two_pole_system();
        let c = vec![vec![1.0, 0.0]];
        let reduced = balanced_truncation(&a, &b, &c, 2).expect("stable");
        assert_eq!(reduced.0.len(), 1, "only one minimal state remains");
        let full = (a, b, c);
        for omega in [0.0, 0.3, 1.0, 5.0] {
            assert!(
                response_error(&full, &reduced, omega) < 1e-12,
                "ω = {omega}"
            );
        }
        assert!(
            (reduced.0[0][0] + 1.0).abs() < 1e-10,
            "the retained pole is -1"
        );
    }

    #[test]
    fn full_order_truncation_preserves_the_transfer_function() {
        let n = 5;
        let full = (
            random_stable(n, 81, 0.5),
            random_matrix(n, 1, 82),
            random_matrix(1, n, 83),
        );
        let reduced = balanced_truncation(&full.0, &full.1, &full.2, n).expect("stable");
        assert_eq!(reduced.0.len(), n);
        for omega in [0.0, 0.2, 1.0, 4.0, 25.0] {
            let scale = {
                let (re, im) = siso_response(&full.0, &full.1, &full.2, omega);
                re.hypot(im)
            };
            assert!(
                response_error(&full, &reduced, omega) < 1e-8 * (1.0 + scale),
                "ω = {omega}"
            );
        }
    }

    #[test]
    fn truncation_error_respects_the_hinf_bound() {
        // 1-D diffusion chain: A = tridiag(1, -2, 1), input at one end, output
        // mostly at the other.  Its Hankel singular values decay over six
        // decades, so every truncation order gives a meaningful bound.
        let n = 6;
        let mut a = diagonal(&vec![-2.0; n]);
        for i in 0..n - 1 {
            a[i][i + 1] = 1.0;
            a[i + 1][i] = 1.0;
        }
        let mut b = vec![vec![0.0]; n];
        b[0][0] = 1.0;
        let mut c = vec![vec![0.0; n]];
        c[0][0] = 0.5;
        c[0][n - 1] = 1.0;
        let full = (a, b, c);
        for r in 1..n {
            let transform = balancing_transform(&full.0, &full.1, &full.2, r).expect("stable");
            assert_eq!(transform.reduced_dim(), r);
            let reduced = transform.reduce(&full.0, &full.1, &full.2);
            let bound = transform.error_bound();
            let tail: f64 = transform.hankel_singular_values[r..].iter().sum();
            assert!((bound - 2.0 * tail).abs() < 1e-15);
            let worst = (0..300)
                .map(|k| {
                    let omega = if k == 0 {
                        0.0
                    } else {
                        10f64.powf(-2.0 + 4.0 * k as f64 / 299.0)
                    };
                    response_error(&full, &reduced, omega)
                })
                .fold(0.0_f64, f64::max);
            assert!(
                worst <= bound * (1.0 + 1e-6),
                "r = {r}: error {worst} exceeds the bound {bound}"
            );
            // The bound is not vacuous: the error is of the order of the tail.
            assert!(
                worst > 0.05 * transform.hankel_singular_values[r],
                "r = {r}"
            );
        }
    }

    #[test]
    fn reduced_model_is_balanced_with_the_leading_hankel_singular_values() {
        let n = 6;
        let (a, b, c) = (
            random_stable(n, 101, 0.4),
            random_matrix(n, 2, 102),
            random_matrix(3, n, 103),
        );
        let r = 3;
        let transform = balancing_transform(&a, &b, &c, r).expect("stable");
        let (a_r, b_r, c_r) = transform.reduce(&a, &b, &c);
        let sigma = diagonal(&transform.hankel_singular_values[..r]);
        let bbt = mat_mul(&b_r, &mat_transpose(&b_r));
        let ctc = mat_mul(&mat_transpose(&c_r), &c_r);
        let scale = transform.hankel_singular_values[0];
        // A_r Σ + Σ A_rᵀ + B_r B_rᵀ = 0 and A_rᵀ Σ + Σ A_r + C_rᵀ C_r = 0
        assert!(lyapunov_residual(&a_r, &sigma, &bbt) < 1e-9 * scale);
        assert!(lyapunov_residual(&mat_transpose(&a_r), &sigma, &ctc) < 1e-9 * scale);
    }

    #[test]
    fn projectors_are_biorthogonal_for_every_order() {
        let n = 6;
        let (a, b, c) = (
            random_stable(n, 111, 0.4),
            random_matrix(n, 1, 112),
            random_matrix(2, n, 113),
        );
        for r in 1..=n {
            let t = balancing_transform(&a, &b, &c, r).expect("stable");
            let product = mat_mul(&t.t_left, &t.t_right);
            for (i, row) in product.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!((value - want).abs() < 1e-9, "r = {r}: T_l T_r != I");
                }
            }
        }
    }

    #[test]
    fn repeated_hankel_singular_values_are_handled() {
        // Three identical decoupled modes: sigma = (1/2, 1/2, 1/2).
        let a = diagonal(&[-1.0, -1.0, -1.0]);
        let b = identity(3);
        let c = identity(3);
        let hsv = hankel_singular_values(&a, &b, &c).expect("stable");
        for s in &hsv {
            assert!((s - 0.5).abs() < 1e-13, "hsv = {hsv:?}");
        }
        let (a_r, b_r, c_r) = balanced_truncation(&a, &b, &c, 2).expect("stable");
        assert_eq!((a_r.len(), a_r[0].len()), (2, 2));
        assert_eq!((b_r.len(), b_r[0].len()), (2, 3));
        assert_eq!((c_r.len(), c_r[0].len()), (3, 2));
        for (i, row) in a_r.iter().enumerate() {
            assert!((row[i] + 1.0).abs() < 1e-10, "poles stay at -1");
        }
    }

    #[test]
    fn large_system_is_reduced_through_the_iterative_lyapunov_solver() {
        let n = KRONECKER_MAX_ORDER + 5;
        let poles: Vec<f64> = (0..n).map(|i| -(0.2 + 0.35 * i as f64)).collect();
        let a = diagonal(&poles);
        let b = random_matrix(n, 1, 121);
        let c = random_matrix(1, n, 122);
        let r = 8;
        let transform = balancing_transform(&a, &b, &c, r).expect("stable");
        let reduced = transform.reduce(&a, &b, &c);
        let full = (a, b, c);
        let bound = transform.error_bound();
        let worst = (0..200)
            .map(|k| {
                let omega = if k == 0 {
                    0.0
                } else {
                    10f64.powf(-2.0 + 3.5 * k as f64 / 199.0)
                };
                response_error(&full, &reduced, omega)
            })
            .fold(0.0_f64, f64::max);
        assert!(
            worst <= bound * (1.0 + 1e-8),
            "error {worst} vs bound {bound}"
        );
        assert!(
            bound < 1.0,
            "the leading modes should capture the system: {bound}"
        );
    }

    #[test]
    fn unstable_and_marginal_systems_are_rejected() {
        let b = vec![vec![1.0], vec![0.5]];
        let c = vec![vec![1.0, 1.0]];
        let unstable = diagonal(&[1.0, 2.0]);
        assert!(matches!(
            balanced_truncation(&unstable, &b, &c, 1),
            Err(FemError::InvalidInput(_))
        ));
        let mixed = diagonal(&[-1.0, 0.5]);
        assert!(matches!(
            hankel_singular_values(&mixed, &b, &c),
            Err(FemError::InvalidInput(_)) | Err(FemError::SingularSystem)
        ));
        let marginal = diagonal(&[0.0, -1.0]);
        assert_eq!(
            balanced_truncation(&marginal, &b, &c, 1).err(),
            Some(FemError::SingularSystem)
        );
    }

    #[test]
    fn degenerate_inputs_give_empty_or_zero_reductions() {
        let empty = balanced_truncation(&[], &[], &[], 3).expect("empty system");
        assert!(empty.0.is_empty() && empty.1.is_empty() && empty.2.is_empty());
        let (a, b, c) = two_pole_system();
        let zero_modes = balanced_truncation(&a, &b, &c, 0).expect("stable");
        assert!(zero_modes.0.is_empty() && zero_modes.1.is_empty() && zero_modes.2.is_empty());
        assert!(
            hankel_singular_values(&[], &[], &[])
                .expect("empty")
                .is_empty()
        );
        // No input at all: nothing is controllable, so nothing is retained.
        let no_input = vec![Vec::new(), Vec::new()];
        let (a_r, b_r, c_r) = balanced_truncation(&a, &no_input, &c, 2).expect("stable");
        assert!(a_r.is_empty() && b_r.is_empty());
        assert_eq!(c_r, vec![Vec::<f64>::new()]);
        let hsv = hankel_singular_values(&a, &no_input, &c).expect("stable");
        assert_eq!(hsv, vec![0.0, 0.0]);
    }

    #[test]
    fn shapes_are_validated() {
        let (a, b, c) = two_pole_system();
        assert!(matches!(
            balanced_truncation(&a, &[vec![1.0]], &c, 1),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            balanced_truncation(&a, &b, &[vec![1.0, 1.0, 1.0]], 1),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            balanced_truncation(&[vec![-1.0, 0.0]], &b, &c, 1),
            Err(FemError::InvalidInput(_))
        ));
    }
}
