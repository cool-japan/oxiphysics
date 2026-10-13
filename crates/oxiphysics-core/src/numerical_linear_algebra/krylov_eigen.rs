//! Krylov-subspace eigensolvers for sparse matrices: Arnoldi and Lanczos.
//!
//! Both solvers grow an orthonormal basis `V_m` of the Krylov space
//! `span{q_0, A q_0, ..., A^(m-1) q_0}` with the Arnoldi process, which yields
//! the factorisation
//!
//! ```text
//! A V_m = V_m H_m + beta_m q_m e_m^T
//! ```
//!
//! with an upper-Hessenberg `H_m` (tridiagonal when `A` is symmetric).  The
//! eigenpairs `(theta, y)` of `H_m` give the Ritz pairs `(theta, V_m y)` of `A`,
//! whose residual is `sqrt(||H_m y - theta y||^2 + (beta_m y[m-1])^2)`.
//!
//! Design points that keep the results honest:
//!
//! * The basis is fully reorthogonalised (Daniel-Gragg-Kaufman-Stewart
//!   "twice is enough"), so no ghost Ritz values appear.
//! * The start vector is a fixed pseudo-random vector rather than the constant
//!   vector: a constant vector is orthogonal to every antisymmetric eigenvector
//!   of a mirror-symmetric matrix such as the 1-D Laplacian and would silently
//!   hide half of the spectrum.
//! * When the Krylov space becomes invariant ("lucky breakdown") the process is
//!   restarted with a fresh vector orthogonal to the basis; the Hessenberg
//!   sub-diagonal is then exactly 0.  This is the only way a single-vector Krylov
//!   method can discover repeated eigenvalues, so once a breakdown has occurred
//!   the early exit on converged Ritz pairs is disabled and the whole step budget
//!   is used.  Before a breakdown a converged Ritz value stands for *all* copies
//!   of an eigenvalue; use LOBPCG (block method) when high multiplicities are
//!   expected.
//! * `converged` is decided from the true residual `||A x - theta x|| <= tol`
//!   (an absolute tolerance), never from a recurrence alone.

use super::functions::{axpy, dot, norm2};
use super::types::{EigenResult, EigenSolver, SparseCSR};
use crate::numerical_methods::{
    symmetric_tridiagonal_eigen, symmetric_tridiagonal_eigen_last_components,
};

/// Seed of the generator that produces start, restart and inverse-iteration vectors.
const START_VECTOR_SEED: u64 = 0x2545_F491_4F6C_DD1D;
/// DGKS threshold: a projection pass that keeps at least this fraction of the
/// norm needs no further refinement.
const REORTHOGONALISATION_ETA: f64 = std::f64::consts::FRAC_1_SQRT_2;
/// Absolute breakdown threshold in units of `eps * ||A||`.
const BREAKDOWN_FACTOR: f64 = 100.0;
/// How many random vectors are tried when restarting after a breakdown.
const RESTART_ATTEMPTS: usize = 4;
/// Imaginary parts below this (relative to `max(1, |Re|)`) count as real.
const REAL_RITZ_TOLERANCE: f64 = 1.0e-10;
/// Ritz values closer than this (relative) form a cluster whose eigenvectors
/// are kept linearly independent.
const CLUSTER_TOLERANCE: f64 = 1.0e-9;
/// Maximum number of inverse-iteration solves per Ritz vector.
const INVERSE_ITERATION_LIMIT: usize = 6;

/// Deterministic SplitMix64 generator for reproducible pseudo-random vectors.
pub(crate) struct VectorGenerator {
    state: u64,
}

impl VectorGenerator {
    /// Creates a generator from `seed`.
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Next value, uniform on `[-1, 1)`.
    pub(crate) fn next_signed(&mut self) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * unit - 1.0
    }

    /// Pseudo-random vector of length `n` scaled to unit 2-norm.
    pub(crate) fn unit_vector(&mut self, n: usize) -> Vec<f64> {
        let mut v: Vec<f64> = (0..n).map(|_| self.next_signed()).collect();
        let nrm = norm2(&v);
        if nrm > 0.0 {
            v.iter_mut().for_each(|x| *x /= nrm);
        } else if let Some(first) = v.first_mut() {
            *first = 1.0;
        }
        v
    }
}

/// Dot product with four independent accumulators.
///
/// Splitting the sum breaks the serial dependency of a left-to-right reduction so
/// the loop vectorises; the result differs from the plain sum only in the last
/// bits.  Orthogonalisation spends almost all of its time in this kernel.
pub(crate) fn dot_blocked(a: &[f64], b: &[f64]) -> f64 {
    let mut lanes = [0.0_f64; 4];
    let mut chunks_a = a.chunks_exact(4);
    let mut chunks_b = b.chunks_exact(4);
    for (ca, cb) in (&mut chunks_a).zip(&mut chunks_b) {
        for ((lane, x), y) in lanes.iter_mut().zip(ca.iter()).zip(cb.iter()) {
            *lane += x * y;
        }
    }
    let tail: f64 = chunks_a
        .remainder()
        .iter()
        .zip(chunks_b.remainder().iter())
        .map(|(x, y)| x * y)
        .sum();
    (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]) + tail
}

/// Subtracts from `w` its components along the orthonormal vectors `basis`
/// (modified Gram-Schmidt) and adds the projection coefficients to `coeffs`.
fn project_out(basis: &[Vec<f64>], w: &mut [f64], coeffs: &mut [f64]) {
    for (q, coeff) in basis.iter().zip(coeffs.iter_mut()) {
        let c = dot_blocked(q, w);
        axpy(-c, q, w);
        *coeff += c;
    }
}

/// Orthogonalises `w` against the orthonormal `basis` with DGKS refinement.
///
/// Returns the norm of `w` afterwards and whether `w` lay (numerically) inside
/// the span of `basis`, i.e. whether what is left is rounding noise.
fn orthogonalise_against(basis: &[Vec<f64>], w: &mut [f64], coeffs: &mut [f64]) -> (f64, bool) {
    let norm_in = norm2(w);
    if norm_in == 0.0 {
        return (0.0, true);
    }
    project_out(basis, w, coeffs);
    let norm_first = norm2(w);
    if norm_first >= REORTHOGONALISATION_ETA * norm_in {
        return (norm_first, false);
    }
    project_out(basis, w, coeffs);
    let norm_second = norm2(w);
    (
        norm_second,
        norm_second < REORTHOGONALISATION_ETA * norm_first,
    )
}

/// Incrementally grown, fully reorthogonalised Arnoldi factorisation.
pub(crate) struct ArnoldiFactorization {
    n: usize,
    breakdown_tol: f64,
    generator: VectorGenerator,
    /// Orthonormal basis vectors `q_0, q_1, ...`; one more than the number of
    /// completed steps unless the space is exhausted.
    vectors: Vec<Vec<f64>>,
    /// Number of invariant subspaces found so far (breakdowns).
    breakdowns: usize,
    /// `columns[j]` holds column `j` of the Hessenberg matrix: the entries
    /// `h_{0..=j, j}` followed by the sub-diagonal `h_{j+1, j}` (length `j + 2`).
    columns: Vec<Vec<f64>>,
}

impl ArnoldiFactorization {
    /// Starts a factorisation in dimension `n >= 1`; `matrix_norm` is any upper
    /// bound of `||A||` and sets the absolute breakdown threshold.
    pub(crate) fn new(n: usize, matrix_norm: f64) -> Self {
        let mut generator = VectorGenerator::new(START_VECTOR_SEED);
        let start = generator.unit_vector(n);
        Self {
            n,
            breakdown_tol: BREAKDOWN_FACTOR * f64::EPSILON * matrix_norm.max(f64::MIN_POSITIVE),
            generator,
            vectors: vec![start],
            breakdowns: 0,
            columns: Vec::new(),
        }
    }

    /// Number of completed Arnoldi steps `m`.
    pub(crate) fn steps(&self) -> usize {
        self.columns.len()
    }

    /// Whether an invariant subspace has been found (and the process restarted).
    ///
    /// Once that has happened the Ritz pairs found so far are exact eigenpairs
    /// but nothing certifies that they are the dominant ones: the part of the
    /// space the start vector did not reach has still to be explored.
    pub(crate) fn has_broken_down(&self) -> bool {
        self.breakdowns > 0
    }

    /// Norm `beta_m` of the residual vector after the last completed step (zero
    /// after a breakdown and once the basis spans the whole space).
    pub(crate) fn residual_norm(&self) -> f64 {
        self.columns
            .last()
            .and_then(|column| column.last())
            .copied()
            .unwrap_or(0.0)
    }

    /// Entry `h_{i,j}` of the Hessenberg matrix (zero outside the computed band).
    pub(crate) fn entry(&self, i: usize, j: usize) -> f64 {
        self.columns
            .get(j)
            .and_then(|column| column.get(i))
            .copied()
            .unwrap_or(0.0)
    }

    /// Leading `m x m` block of the Hessenberg matrix, row-major.
    pub(crate) fn hessenberg(&self, m: usize) -> Vec<f64> {
        let m = m.min(self.columns.len());
        let mut h = vec![0.0_f64; m * m];
        for (j, column) in self.columns.iter().take(m).enumerate() {
            for (i, &value) in column.iter().enumerate().take(m) {
                h[i * m + j] = value;
            }
        }
        h
    }

    /// Ritz vector `V_m y` for the coefficient vector `y` of length `m`.
    pub(crate) fn ritz_vector(&self, y: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0_f64; self.n];
        for (q, &coefficient) in self.vectors.iter().zip(y.iter()) {
            axpy(coefficient, q, &mut x);
        }
        x
    }

    /// Performs one Arnoldi step.  Returns `false` when no step is possible any
    /// more (the basis already spans the whole space).
    pub(crate) fn extend(&mut self, matvec: &mut dyn FnMut(&[f64]) -> Vec<f64>) -> bool {
        let j = self.columns.len();
        if j >= self.n {
            return false;
        }
        let mut w = match self.vectors.get(j) {
            Some(q) => matvec(q),
            None => return false,
        };
        if w.len() != self.n {
            return false;
        }
        let mut column = vec![0.0_f64; j + 2];
        let (norm, in_span) = orthogonalise_against(&self.vectors[..=j], &mut w, &mut column[..=j]);
        let room_for_another = self.vectors.len() < self.n;
        if room_for_another && !in_span && norm > self.breakdown_tol {
            column[j + 1] = norm;
            w.iter_mut().for_each(|x| *x /= norm);
            self.vectors.push(w);
        } else {
            // Invariant subspace (or the basis is complete): the sub-diagonal is
            // exactly zero.  Continue in the orthogonal complement if there is one.
            column[j + 1] = 0.0;
            if room_for_another {
                self.breakdowns += 1;
                if let Some(fresh) = self.restart_vector() {
                    self.vectors.push(fresh);
                }
            }
        }
        self.columns.push(column);
        true
    }

    /// A unit vector orthogonal to the current basis, if one can be found.
    fn restart_vector(&mut self) -> Option<Vec<f64>> {
        for _ in 0..RESTART_ATTEMPTS {
            let mut candidate = self.generator.unit_vector(self.n);
            let mut scratch = vec![0.0_f64; self.vectors.len()];
            let (norm, in_span) =
                orthogonalise_against(&self.vectors, &mut candidate, &mut scratch);
            if !in_span && norm > 0.0 {
                candidate.iter_mut().for_each(|x| *x /= norm);
                return Some(candidate);
            }
        }
        None
    }
}

/// Largest absolute row sum of `a`, an upper bound of its spectral radius.
pub(crate) fn infinity_norm(a: &SparseCSR) -> f64 {
    (0..a.nrows)
        .map(|i| {
            a.values[a.row_ptr[i]..a.row_ptr[i + 1]]
                .iter()
                .map(|v| v.abs())
                .sum::<f64>()
        })
        .fold(0.0, f64::max)
}

/// Indices of the `k` entries of `values` with the largest magnitude (ties go to
/// the larger signed value), returned in order of increasing value.
fn dominant_indices(values: &[f64], k: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| {
        values[b]
            .abs()
            .total_cmp(&values[a].abs())
            .then(values[b].total_cmp(&values[a]))
    });
    order.truncate(k);
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    order
}

/// A Ritz value with its coefficient vector in the Krylov basis.
struct RitzPair {
    value: f64,
    /// Coefficients `y` of the Ritz vector `V_m y` (empty for estimates only).
    y: Vec<f64>,
    /// Estimated residual norm `||A V_m y - value V_m y||`.
    estimate: f64,
}

/// Whether `pairs` holds the `wanted` pairs and all residual estimates are
/// within the absolute tolerance `tol` (a NaN estimate never is).
fn estimates_converged(pairs: &[RitzPair], wanted: usize, tol: f64) -> bool {
    pairs.len() == wanted && pairs.iter().all(|p| p.estimate <= tol)
}

/// Cheap residual estimates of the dominant Ritz values of the symmetric
/// tridiagonal projection, `O(m^2)`.
fn lanczos_estimates(fact: &ArnoldiFactorization, wanted: usize) -> Vec<RitzPair> {
    let (alpha, beta) = tridiagonal_projection(fact);
    let (theta, last) = symmetric_tridiagonal_eigen_last_components(&alpha, &beta);
    let residual = fact.residual_norm();
    dominant_indices(&theta, wanted)
        .into_iter()
        .map(|i| RitzPair {
            value: theta[i],
            y: Vec::new(),
            estimate: residual * last[i].abs(),
        })
        .collect()
}

/// Dominant Ritz pairs of the symmetric tridiagonal projection.
fn lanczos_pairs(fact: &ArnoldiFactorization, wanted: usize) -> Vec<RitzPair> {
    let (alpha, beta) = tridiagonal_projection(fact);
    let (theta, vectors) = symmetric_tridiagonal_eigen(&alpha, &beta);
    let residual = fact.residual_norm();
    dominant_indices(&theta, wanted)
        .into_iter()
        .map(|i| {
            let y = vectors[i].clone();
            let estimate = residual * y.last().copied().unwrap_or(0.0).abs();
            RitzPair {
                value: theta[i],
                y,
                estimate,
            }
        })
        .collect()
}

/// Diagonal and off-diagonal of the symmetric tridiagonal matrix `T_m`
/// extracted from the Hessenberg matrix (sub- and super-diagonal averaged).
fn tridiagonal_projection(fact: &ArnoldiFactorization) -> (Vec<f64>, Vec<f64>) {
    let m = fact.steps();
    let alpha: Vec<f64> = (0..m).map(|j| fact.entry(j, j)).collect();
    let beta: Vec<f64> = (0..m.saturating_sub(1))
        .map(|j| 0.5 * (fact.entry(j + 1, j) + fact.entry(j, j + 1)))
        .collect();
    (alpha, beta)
}

/// Dominant real Ritz pairs of the Hessenberg matrix `H_m`.
///
/// Complex-conjugate Ritz pairs have no real eigenvector and cannot be reported
/// through [`EigenResult`], so only (numerically) real Ritz values qualify.
fn arnoldi_pairs(fact: &ArnoldiFactorization, wanted: usize) -> Vec<RitzPair> {
    let m = fact.steps();
    let h = fact.hessenberg(m);
    let real_values: Vec<f64> = EigenSolver::hessenberg_eigenvalues_complex(&h, m)
        .into_iter()
        .filter(|&(re, im)| im.abs() <= REAL_RITZ_TOLERANCE * re.abs().max(1.0))
        .map(|(re, _)| re)
        .collect();
    let norm_h = (0..m)
        .map(|i| h[i * m..(i + 1) * m].iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let residual = fact.residual_norm();
    let mut pairs: Vec<RitzPair> = Vec::new();
    for (position, index) in dominant_indices(&real_values, wanted)
        .into_iter()
        .enumerate()
    {
        let theta = real_values[index];
        let cluster: Vec<&[f64]> = pairs
            .iter()
            .filter(|p| (p.value - theta).abs() <= CLUSTER_TOLERANCE * theta.abs().max(1.0))
            .map(|p| p.y.as_slice())
            .collect();
        let (y, hessenberg_residual) =
            hessenberg_eigenvector(&h, m, theta, norm_h, &cluster, position as u64);
        let last = y.last().copied().unwrap_or(0.0);
        pairs.push(RitzPair {
            value: theta,
            y,
            estimate: hessenberg_residual.hypot(residual * last.abs()),
        });
    }
    pairs
}

/// LU factorisation with partial pivoting of the shifted upper-Hessenberg matrix
/// `H - shift I` (one multiplier per column, `O(m^2)`).
struct ShiftedHessenbergLu {
    m: usize,
    /// Row-major upper-triangular factor.
    upper: Vec<f64>,
    multipliers: Vec<f64>,
    swapped: Vec<bool>,
}

impl ShiftedHessenbergLu {
    /// Factors `H - shift I`; pivots smaller than `floor` are replaced by `floor`
    /// so that an (almost) singular shift still yields a usable factorisation.
    fn factor(h: &[f64], m: usize, shift: f64, floor: f64) -> Self {
        let mut u = h.to_vec();
        for i in 0..m {
            u[i * m + i] -= shift;
        }
        let mut multipliers = vec![0.0_f64; m.saturating_sub(1)];
        let mut swapped = vec![false; m.saturating_sub(1)];
        for k in 0..m.saturating_sub(1) {
            if u[(k + 1) * m + k].abs() > u[k * m + k].abs() {
                for j in k..m {
                    u.swap(k * m + j, (k + 1) * m + j);
                }
                swapped[k] = true;
            }
            if u[k * m + k].abs() < floor {
                u[k * m + k] = floor;
            }
            let multiplier = u[(k + 1) * m + k] / u[k * m + k];
            multipliers[k] = multiplier;
            for j in k..m {
                let pivot_row_entry = u[k * m + j];
                u[(k + 1) * m + j] -= multiplier * pivot_row_entry;
            }
            u[(k + 1) * m + k] = 0.0;
        }
        if m > 0 && u[(m - 1) * m + (m - 1)].abs() < floor {
            u[(m - 1) * m + (m - 1)] = floor;
        }
        Self {
            m,
            upper: u,
            multipliers,
            swapped,
        }
    }

    /// Solves `(H - shift I) z = rhs` in place.
    fn solve(&self, rhs: &mut [f64]) {
        let m = self.m;
        for k in 0..m.saturating_sub(1) {
            if self.swapped[k] {
                rhs.swap(k, k + 1);
            }
            rhs[k + 1] -= self.multipliers[k] * rhs[k];
        }
        for i in (0..m).rev() {
            let tail: f64 = (i + 1..m).map(|j| self.upper[i * m + j] * rhs[j]).sum();
            rhs[i] = (rhs[i] - tail) / self.upper[i * m + i];
        }
    }
}

/// Unit eigenvector of the upper-Hessenberg matrix `h` (`m x m`, row-major) for
/// the real eigenvalue `theta`, by inverse iteration, together with its residual
/// `||H y - theta y||`.
///
/// `cluster` lists already computed eigenvectors of (numerically) the same
/// eigenvalue; the iterate is kept orthogonal to them so that a repeated
/// eigenvalue yields independent vectors.
fn hessenberg_eigenvector(
    h: &[f64],
    m: usize,
    theta: f64,
    norm_h: f64,
    cluster: &[&[f64]],
    salt: u64,
) -> (Vec<f64>, f64) {
    if m <= 1 {
        return (vec![1.0; m], 0.0);
    }
    let floor = f64::EPSILON * norm_h.max(f64::MIN_POSITIVE);
    let lu = ShiftedHessenbergLu::factor(h, m, theta, floor);
    let mut generator = VectorGenerator::new(START_VECTOR_SEED ^ salt.wrapping_mul(0x9E37_79B9));
    let mut iterate = generator.unit_vector(m);
    let mut best = iterate.clone();
    let mut best_residual = f64::INFINITY;
    for _ in 0..INVERSE_ITERATION_LIMIT {
        let mut z = iterate.clone();
        lu.solve(&mut z);
        for previous in cluster {
            let c = dot(previous, &z);
            axpy(-c, previous, &mut z);
        }
        let nrm = norm2(&z);
        if !nrm.is_finite() || nrm == 0.0 {
            break;
        }
        z.iter_mut().for_each(|x| *x /= nrm);
        let residual = hessenberg_residual(h, m, theta, &z);
        if residual < best_residual {
            best_residual = residual;
            best = z.clone();
        }
        if residual <= 4.0 * f64::EPSILON * norm_h {
            break;
        }
        iterate = z;
    }
    if !best_residual.is_finite() {
        best_residual = hessenberg_residual(h, m, theta, &best);
    }
    (best, best_residual)
}

/// `||H y - theta y||_2` for the row-major `m x m` matrix `h`.
fn hessenberg_residual(h: &[f64], m: usize, theta: f64, y: &[f64]) -> f64 {
    (0..m)
        .map(|i| {
            let row: f64 = h[i * m..(i + 1) * m]
                .iter()
                .zip(y.iter())
                .map(|(hij, yj)| hij * yj)
                .sum();
            let r = row - theta * y[i];
            r * r
        })
        .sum::<f64>()
        .sqrt()
}

/// Which projection the driver uses.
#[derive(Clone, Copy, PartialEq, Eq)]
enum KrylovKind {
    /// Symmetric matrix: tridiagonal projection (Lanczos).
    Symmetric,
    /// General matrix: Hessenberg projection (Arnoldi).
    General,
}

/// Turns Ritz pairs into an [`EigenResult`], measuring the true residual
/// `||A x - theta x||` of every Ritz vector `x`.
fn assemble_result(
    a: &SparseCSR,
    fact: &ArnoldiFactorization,
    pairs: &[RitzPair],
    wanted: usize,
    tol: f64,
) -> EigenResult {
    let mut eigenvalues = Vec::with_capacity(pairs.len());
    let mut eigenvectors = Vec::with_capacity(pairs.len());
    let mut converged = pairs.len() == wanted;
    for pair in pairs {
        let mut x = fact.ritz_vector(&pair.y);
        let nrm = norm2(&x);
        if nrm > 0.0 {
            x.iter_mut().for_each(|v| *v /= nrm);
        }
        let mut residual = a.matvec(&x);
        axpy(-pair.value, &x, &mut residual);
        // Written as a positive test so that a NaN residual counts as a failure.
        let within_tolerance = norm2(&residual) <= tol;
        converged &= within_tolerance;
        eigenvalues.push(pair.value);
        eigenvectors.push(x);
    }
    EigenResult {
        eigenvalues,
        eigenvectors,
        iterations: fact.steps(),
        converged,
    }
}

/// Shared driver of [`arnoldi`] and [`lanczos`].
fn krylov_eigenpairs(
    a: &SparseCSR,
    k: usize,
    max_iter: usize,
    tol: f64,
    kind: KrylovKind,
) -> EigenResult {
    let n = a.nrows;
    if a.ncols != n {
        // Eigenpairs are only defined for square matrices.
        return EigenResult {
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            iterations: 0,
            converged: false,
        };
    }
    let wanted = k.min(n);
    if wanted == 0 {
        return EigenResult {
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            iterations: 0,
            converged: true,
        };
    }
    let max_steps = n.min((2 * wanted + 1).max(max_iter));
    let mut fact = ArnoldiFactorization::new(n, infinity_norm(a));
    let mut matvec = |v: &[f64]| a.matvec(v);
    let mut next_check = wanted;
    loop {
        let advanced = fact.extend(&mut matvec);
        let m = fact.steps();
        if !advanced || m >= max_steps {
            break;
        }
        if m < next_check || fact.has_broken_down() {
            continue;
        }
        next_check = m + (m / 8).max(1);
        let candidate = match kind {
            KrylovKind::Symmetric => {
                if estimates_converged(&lanczos_estimates(&fact, wanted), wanted, tol) {
                    Some(lanczos_pairs(&fact, wanted))
                } else {
                    None
                }
            }
            KrylovKind::General => {
                let pairs = arnoldi_pairs(&fact, wanted);
                estimates_converged(&pairs, wanted, tol).then_some(pairs)
            }
        };
        if let Some(pairs) = candidate {
            let result = assemble_result(a, &fact, &pairs, wanted, tol);
            if result.converged {
                return result;
            }
        }
    }
    let pairs = match kind {
        KrylovKind::Symmetric => lanczos_pairs(&fact, wanted),
        KrylovKind::General => arnoldi_pairs(&fact, wanted),
    };
    assemble_result(a, &fact, &pairs, wanted, tol)
}

/// Arnoldi eigensolver; see [`EigenSolver::arnoldi`].
pub(super) fn arnoldi(a: &SparseCSR, k: usize, max_iter: usize, tol: f64) -> EigenResult {
    krylov_eigenpairs(a, k, max_iter, tol, KrylovKind::General)
}

/// Lanczos eigensolver; see [`EigenSolver::lanczos`].
pub(super) fn lanczos(a: &SparseCSR, k: usize, max_iter: usize, tol: f64) -> EigenResult {
    krylov_eigenpairs(a, k, max_iter, tol, KrylovKind::Symmetric)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense_matvec(rows: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
        rows.iter().map(|row| dot(row, v)).collect()
    }

    #[test]
    fn vector_generator_is_deterministic_and_normalised() {
        let mut first = VectorGenerator::new(42);
        let mut second = VectorGenerator::new(42);
        let a = first.unit_vector(25);
        let b = second.unit_vector(25);
        assert_eq!(a, b);
        assert!((norm2(&a) - 1.0).abs() < 1e-14);
        // Not a constant vector: its entries differ in sign and magnitude.
        assert!(a.iter().any(|v| *v > 0.0) && a.iter().any(|v| *v < 0.0));
        let mut third = VectorGenerator::new(43);
        assert_ne!(a, third.unit_vector(25));
        for _ in 0..1000 {
            let x = first.next_signed();
            assert!((-1.0..1.0).contains(&x));
        }
    }

    #[test]
    fn arnoldi_relation_holds_and_basis_is_orthonormal() {
        let n = 12;
        let a = SparseCSR::laplacian_1d(n);
        let mut fact = ArnoldiFactorization::new(n, infinity_norm(&a));
        let mut matvec = |v: &[f64]| a.matvec(v);
        for _ in 0..6 {
            assert!(fact.extend(&mut matvec));
        }
        assert_eq!(fact.steps(), 6);
        assert!(!fact.has_broken_down());
        let basis = &fact.vectors;
        assert_eq!(basis.len(), 7);
        for (i, qi) in basis.iter().enumerate() {
            for (j, qj) in basis.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot(qi, qj) - expected).abs() < 1e-13, "<q{i}, q{j}>");
            }
        }
        // A q_j = sum_{i <= j+1} h_{ij} q_i for every completed step.
        for j in 0..6 {
            let mut reconstructed = vec![0.0_f64; n];
            for (i, q) in basis.iter().enumerate().take(j + 2) {
                axpy(fact.entry(i, j), q, &mut reconstructed);
            }
            let direct = a.matvec(&basis[j]);
            let error = direct
                .iter()
                .zip(reconstructed.iter())
                .map(|(d, r)| (d - r).abs())
                .fold(0.0, f64::max);
            assert!(error < 1e-13, "step {j}: error {error}");
        }
        // The Laplacian is symmetric, so H is tridiagonal up to rounding.
        let h = fact.hessenberg(6);
        for i in 0..6 {
            for j in 0..6 {
                if i + 1 < j || j + 1 < i {
                    assert!(h[i * 6 + j].abs() < 1e-13, "h[{i}][{j}] = {}", h[i * 6 + j]);
                }
            }
        }
        assert!((fact.residual_norm() - fact.entry(6, 5)).abs() < 1e-15);
        assert!(fact.residual_norm() > 0.0);
    }

    #[test]
    fn breakdown_restarts_with_an_exactly_zero_subdiagonal() {
        // diag(1, 1, 1, 2, 2): the Krylov space of any vector has dimension 2.
        let rows = vec![
            vec![1.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 2.0],
        ];
        let mut fact = ArnoldiFactorization::new(5, 2.0);
        let mut matvec = |v: &[f64]| dense_matvec(&rows, v);
        while fact.extend(&mut matvec) {}
        assert_eq!(fact.steps(), 5, "the process must exhaust the whole space");
        assert!(fact.has_broken_down());
        // The sub-diagonal vanishes exactly after the first invariant subspace.
        assert_eq!(fact.entry(2, 1), 0.0);
        let basis = &fact.vectors;
        assert_eq!(basis.len(), 5);
        for (i, qi) in basis.iter().enumerate() {
            for (j, qj) in basis.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot(qi, qj) - expected).abs() < 1e-12, "<q{i}, q{j}>");
            }
        }
        // H is similar to A, so its spectrum is {1, 1, 1, 2, 2}.
        let h = fact.hessenberg(5);
        let mut ritz = EigenSolver::hessenberg_eigenvalues(&h, 5);
        ritz.sort_by(|x, y| x.total_cmp(y));
        for (got, want) in ritz.iter().zip([1.0, 1.0, 1.0, 2.0, 2.0].iter()) {
            assert!((got - want).abs() < 1e-12, "{got} vs {want}");
        }
        // No further step is possible.
        assert!(!fact.extend(&mut matvec));
    }

    #[test]
    fn orthogonalisation_flags_vectors_inside_the_span() {
        // An orthonormal basis with "awkward" floating-point entries.
        let n = 20;
        let a = SparseCSR::laplacian_1d(n);
        let mut fact = ArnoldiFactorization::new(n, infinity_norm(&a));
        let mut matvec = |v: &[f64]| a.matvec(v);
        for _ in 0..10 {
            assert!(fact.extend(&mut matvec));
        }
        let basis = &fact.vectors[..10];
        // w = sum c_i q_i lies in the span up to rounding: what is left is noise.
        let mut inside = vec![0.0_f64; n];
        for (i, q) in basis.iter().enumerate() {
            axpy(1.0 / (i as f64 + 1.0), q, &mut inside);
        }
        let mut coeffs = vec![0.0_f64; 10];
        let (norm, in_span) = orthogonalise_against(basis, &mut inside, &mut coeffs);
        assert!(in_span, "norm after projection {norm}");
        assert!(norm < 1e-14);
        for (i, c) in coeffs.iter().enumerate() {
            assert!(
                (c - 1.0 / (i as f64 + 1.0)).abs() < 1e-13,
                "coefficient {i}"
            );
        }
        // A vector with a genuine component outside the span is not flagged.
        let mut outside = vec![0.0_f64; n];
        for q in basis {
            axpy(0.5, q, &mut outside);
        }
        let fresh = fact.vectors[10].clone();
        axpy(2.0, &fresh, &mut outside);
        let mut coeffs = vec![0.0_f64; 10];
        let (norm, in_span) = orthogonalise_against(basis, &mut outside, &mut coeffs);
        assert!(!in_span);
        assert!((norm - 2.0).abs() < 1e-12, "norm {norm}");
        for q in basis {
            assert!(dot(q, &outside).abs() < 1e-13);
        }
        // The zero vector counts as inside the span.
        let mut zero = vec![0.0_f64; n];
        let mut coeffs = vec![0.0_f64; 10];
        assert_eq!(
            orthogonalise_against(basis, &mut zero, &mut coeffs),
            (0.0, true)
        );
    }

    /// Non-symmetric tridiagonal matrix: diagonal 2, sub-diagonal -1.3,
    /// super-diagonal -0.7.
    fn convection_diffusion(n: usize) -> SparseCSR {
        let (mut rows, mut cols, mut vals) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            rows.push(i);
            cols.push(i);
            vals.push(2.0);
            if i > 0 {
                rows.push(i);
                cols.push(i - 1);
                vals.push(-1.3);
            }
            if i + 1 < n {
                rows.push(i);
                cols.push(i + 1);
                vals.push(-0.7);
            }
        }
        SparseCSR::from_triplets(n, n, &rows, &cols, &vals)
    }

    fn true_residual(a: &SparseCSR, value: f64, x: &[f64]) -> f64 {
        let mut r = a.matvec(x);
        axpy(-value, x, &mut r);
        norm2(&r)
    }

    #[test]
    fn residual_estimates_equal_the_true_residuals_of_partial_krylov_spaces() {
        // The cheap estimate beta_m |y[m-1]| (Lanczos) and
        // sqrt(||H y - theta y||^2 + (beta_m y[m-1])^2) (Arnoldi) must agree with
        // ||A x - theta x|| long before the Krylov space is exhausted.
        let n = 40;
        for symmetric in [true, false] {
            let a = if symmetric {
                SparseCSR::laplacian_1d(n)
            } else {
                convection_diffusion(n)
            };
            let mut fact = ArnoldiFactorization::new(n, infinity_norm(&a));
            let mut matvec = |v: &[f64]| a.matvec(v);
            for _ in 0..9 {
                assert!(fact.extend(&mut matvec));
            }
            let pairs = if symmetric {
                lanczos_pairs(&fact, 3)
            } else {
                arnoldi_pairs(&fact, 3)
            };
            assert_eq!(pairs.len(), 3);
            for pair in &pairs {
                let mut x = fact.ritz_vector(&pair.y);
                let nrm = norm2(&x);
                x.iter_mut().for_each(|v| *v /= nrm);
                let actual = true_residual(&a, pair.value, &x);
                assert!(actual > 1e-8, "the test needs unconverged pairs: {actual}");
                assert!(
                    (pair.estimate - actual).abs() < 1e-10 * actual.max(1.0),
                    "symmetric = {symmetric}: estimate {} vs true residual {actual}",
                    pair.estimate
                );
            }
            // The cheap O(m^2) estimator of the Lanczos driver agrees as well.
            if symmetric {
                let cheap = lanczos_estimates(&fact, 3);
                for (c, p) in cheap.iter().zip(pairs.iter()) {
                    assert!((c.value - p.value).abs() < 1e-13);
                    assert!((c.estimate - p.estimate).abs() < 1e-10 * p.estimate.max(1e-3));
                }
            }
        }
    }

    #[test]
    fn dominant_indices_select_by_magnitude_and_sort_by_value() {
        let values = [3.0, -5.0, 1.0, 5.0, -2.0];
        // Magnitudes 5 and 5 tie: the larger signed value ranks first, then both
        // are reported in increasing order of value.
        assert_eq!(dominant_indices(&values, 2), vec![1, 3]);
        assert_eq!(dominant_indices(&values, 3), vec![1, 0, 3]);
        assert_eq!(dominant_indices(&values, 99).len(), 5);
        assert!(dominant_indices(&values, 0).is_empty());
    }

    #[test]
    fn hessenberg_inverse_iteration_returns_an_eigenvector() {
        // Companion matrix of (x-1)(x-2)(x-4): the eigenvector for lambda is
        // proportional to (lambda^2, lambda, 1).
        let h = [7.0, -14.0, 8.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let norm_h = 22.0;
        for lambda in [1.0_f64, 2.0, 4.0] {
            let (y, residual) = hessenberg_eigenvector(&h, 3, lambda, norm_h, &[], 0);
            assert!(residual < 1e-12, "lambda {lambda}: residual {residual}");
            assert!((norm2(&y) - 1.0).abs() < 1e-12);
            let expected = [lambda * lambda, lambda, 1.0];
            let scale = norm2(&expected);
            let sign = (y[2] / expected[2]).signum();
            for (got, want) in y.iter().zip(expected.iter()) {
                assert!((got * sign - want / scale).abs() < 1e-9, "{got} vs {want}");
            }
        }
    }

    #[test]
    fn hessenberg_inverse_iteration_separates_a_repeated_eigenvalue() {
        // diag(2, 2, 5) is Hessenberg; the eigenvalue 2 has a two-dimensional
        // eigenspace and the cluster argument must keep the second vector
        // independent of the first.
        let h = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 5.0];
        let (first, r1) = hessenberg_eigenvector(&h, 3, 2.0, 5.0, &[], 0);
        let (second, r2) = hessenberg_eigenvector(&h, 3, 2.0, 5.0, &[first.as_slice()], 1);
        assert!(r1 < 1e-12 && r2 < 1e-12);
        assert!(dot(&first, &second).abs() < 1e-10);
        assert!(first[2].abs() < 1e-10 && second[2].abs() < 1e-10);
    }
}
