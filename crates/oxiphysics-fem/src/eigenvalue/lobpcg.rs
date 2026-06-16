// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LOBPCG eigensolver with optional shift-invert (interior) spectral transform.
//!
//! Implements the Locally Optimal Block Preconditioned Conjugate Gradient method
//! (Knyazev, *Toward the Optimal Preconditioned Eigensolver*, SIAM J. Sci. Comput.
//! 2001) for the generalized symmetric eigenproblem `A x = λ B x`, where `A` and
//! `B` are symmetric and `B` is symmetric positive definite (SPD). `B = None`
//! denotes the standard eigenproblem `A x = λ x`.
//!
//! Two operating modes are provided:
//!
//! * **Standard mode** (`shift = None`): the `k` smallest eigenpairs are computed.
//!   The preconditioned residual uses either one AMG V-cycle (`T ≈ A⁻¹`) or a
//!   Jacobi (diagonal) sweep.
//! * **Interior / shift-invert mode** (`shift = Some(σ)`): the `k` eigenpairs of
//!   `(A, B)` *nearest* `σ` are computed via the genuine spectral transform
//!   `OP = (A − σ_work B)⁻¹ B`. The eigenvalues of `OP` are `μ = 1/(λ − σ_work)`,
//!   so the eigenpairs nearest `σ` are those of largest `|μ|`. `A − σ_work B` is
//!   assembled explicitly over the union sparsity pattern and is symmetric
//!   **indefinite** (CG/PCG break down), so `OP` is applied by restarted GMRES.
//!   `σ_work = σ + 1e-7·(1 + |σ|)` perturbs `σ` off any exact eigenvalue. The
//!   shift-invert operator `OP` doubles as the preconditioner for the residual,
//!   which is what drives LOBPCG to the interior of the spectrum.
//!
//! ## Numerical core
//!
//! The Rayleigh–Ritz subspace `S = [X | W | P]` (current iterate, preconditioned
//! residual, conjugate direction) is `B`-orthonormalized by a re-orthogonalizing
//! modified Gram–Schmidt in the `B`-inner product, so the projected pencil reduces
//! to a *standard* dense symmetric eigenproblem. That small problem is solved
//! completely and robustly by **cyclic Jacobi** rotations (no reliance on the
//! Lanczos-approximate generalized path).

use std::collections::BTreeMap;

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use crate::parallel_solver::{CsrMatrix, ParallelGmresSolver, ParallelPcgSolver};
use crate::solvers::amg::classical::AmgClassical;
use crate::solvers::amg::cycle::CycleKind;
use crate::solvers::amg::preconditioner::{AmgPreconditioner, Preconditioner};

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Error type returned by the LOBPCG solver (alias of the crate error).
pub type EigensolverError = crate::Error;

/// Configuration for [`lobpcg_solve`].
#[derive(Debug, Clone)]
pub struct LobpcgConfig {
    /// Number of eigenpairs to compute (block size `k`).
    pub num_eigenvalues: usize,
    /// Maximum number of LOBPCG outer iterations.
    pub max_iter: usize,
    /// Relative residual tolerance: a column is converged when
    /// `‖A x − λ B x‖₂ / (|λ| + 1) < tol`.
    pub tol: f64,
    /// Optional interior target. `None` requests the smallest eigenpairs;
    /// `Some(σ)` requests the eigenpairs nearest `σ` via shift-invert.
    pub shift: Option<f64>,
    /// In **standard** mode, selects the AMG V-cycle preconditioner (`true`) or
    /// the Jacobi/diagonal preconditioner (`false`). In **interior** mode this
    /// flag is ignored: AMG on the indefinite shifted operator is unreliable, so
    /// diagonal-preconditioned restarted GMRES is always used for the
    /// shift-invert solves.
    pub use_amg_preconditioner: bool,
}

impl Default for LobpcgConfig {
    fn default() -> Self {
        Self {
            num_eigenvalues: 1,
            max_iter: 300,
            tol: 1e-6,
            shift: None,
            use_amg_preconditioner: true,
        }
    }
}

/// Result of [`lobpcg_solve`].
///
/// The eigenpairs are returned **ascending by eigenvalue**: `eigenvalues[j]` is
/// the `j`-th eigenvalue and `eigenvectors[j]` is its corresponding eigenvector
/// (a vector of length `n`, the matrix dimension). In interior mode the
/// eigenvalues are the genuine eigenvalues of `(A, B)` nearest `σ` (recovered by
/// Rayleigh quotient), still sorted ascending.
#[derive(Debug, Clone)]
pub struct LobpcgResult {
    /// Eigenvalues, ascending.
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors; `eigenvectors[j]` (length `n`) pairs with `eigenvalues[j]`.
    pub eigenvectors: Vec<Vec<f64>>,
    /// Number of outer LOBPCG iterations performed.
    pub iterations: usize,
    /// Whether every requested column reached the tolerance.
    pub converged: bool,
}

/// Solve the symmetric (generalized) eigenproblem with LOBPCG.
///
/// * `a` — symmetric system matrix (square, `n × n`).
/// * `b` — optional symmetric positive-definite metric matrix; `None` means
///   `B = I` (standard eigenproblem).
/// * `config` — see [`LobpcgConfig`].
///
/// Returns the requested eigenpairs ascending by eigenvalue (see [`LobpcgResult`]).
///
/// # Errors
///
/// Returns [`EigensolverError`] if the dimensions are inconsistent, the matrix is
/// empty, `num_eigenvalues` is zero or exceeds `n`, or the initial block cannot be
/// made full rank.
pub fn lobpcg_solve(
    a: &CsrMatrix,
    b: Option<&CsrMatrix>,
    config: &LobpcgConfig,
) -> Result<LobpcgResult, EigensolverError> {
    let n = a.nrows;
    let k = config.num_eigenvalues;

    if n == 0 {
        return Err(EigensolverError::General(
            "LOBPCG: matrix dimension is zero".into(),
        ));
    }
    if a.ncols != n {
        return Err(EigensolverError::General(
            "LOBPCG: matrix A must be square".into(),
        ));
    }
    if k == 0 {
        return Err(EigensolverError::General(
            "LOBPCG: num_eigenvalues must be positive".into(),
        ));
    }
    if k > n {
        return Err(EigensolverError::General(
            "LOBPCG: num_eigenvalues exceeds the matrix dimension".into(),
        ));
    }
    if let Some(bm) = b {
        if bm.nrows != n || bm.ncols != n {
            return Err(EigensolverError::General(
                "LOBPCG: B must have the same dimensions as A".into(),
            ));
        }
    }

    match config.shift {
        None => solve_standard(a, b, config),
        Some(sigma) => solve_interior(a, b, config, sigma),
    }
}

/// Standard mode: the `k` smallest eigenpairs of `(A, B)`.
fn solve_standard(
    a: &CsrMatrix,
    b: Option<&CsrMatrix>,
    config: &LobpcgConfig,
) -> Result<LobpcgResult, EigensolverError> {
    let n = a.nrows;
    let spec = DriverSpec {
        k: config.num_eigenvalues,
        max_iter: config.max_iter,
        tol: config.tol,
        rr_b_weight: false,
        selection: Selection::SmallestK,
        recover_rayleigh: false,
    };

    if config.use_amg_preconditioner {
        // Build the AMG hierarchy and preconditioner ONCE and reuse every iter.
        let mut amg = AmgClassical::new();
        amg.coarse_cutoff = 8;
        let hierarchy = amg.build(a);
        let precond = AmgPreconditioner {
            hierarchy,
            cycle_kind: CycleKind::V,
            pcg: ParallelPcgSolver::new(50, 1e-8),
        };
        lobpcg_driver(
            a,
            b,
            &spec,
            |v| spmv_vec(a, v),
            |r| {
                let mut z = vec![0.0f64; n];
                precond.apply(r, &mut z);
                z
            },
        )
    } else {
        // Jacobi (diagonal) preconditioner: T r = r .* diag(A)⁻¹ (precomputed once).
        let diag_inv = a.diagonal_preconditioner();
        lobpcg_driver(
            a,
            b,
            &spec,
            |v| spmv_vec(a, v),
            |r| {
                r.iter()
                    .zip(diag_inv.iter())
                    .map(|(ri, di)| ri * di)
                    .collect()
            },
        )
    }
}

/// Interior mode: the `k` eigenpairs nearest `σ` via genuine shift-invert.
fn solve_interior(
    a: &CsrMatrix,
    b: Option<&CsrMatrix>,
    config: &LobpcgConfig,
    sigma: f64,
) -> Result<LobpcgResult, EigensolverError> {
    let n = a.nrows;
    // Perturb σ off any exact eigenvalue to keep A − σ_work B nonsingular.
    let sigma_work = sigma + 1e-7 * (1.0 + sigma.abs());
    let a_shift = assemble_a_minus_sigma_b(a, b, sigma_work);

    // A_shift is symmetric indefinite ⇒ GMRES (CG/PCG break down). AMG on an
    // indefinite operator is unreliable, so diagonal-preconditioned GMRES is
    // always used here, regardless of `config.use_amg_preconditioner`.
    let gmres = ParallelGmresSolver {
        krylov_dim: 60,
        max_restarts: 50,
        tolerance: 1e-9,
    };

    // OP v = (A − σ_work B)⁻¹ (B v). Used as BOTH the Rayleigh–Ritz operator and
    // the residual preconditioner (the shift-invert accelerator).
    let op = |v: &[f64]| -> Vec<f64> {
        let rhs = apply_b_vec(b, v);
        let mut w = vec![0.0f64; n];
        let _ = gmres.solve(&a_shift, &rhs, &mut w);
        if w.iter().all(|x| x.is_finite()) {
            w
        } else {
            vec![0.0f64; n]
        }
    };

    let spec = DriverSpec {
        k: config.num_eigenvalues,
        max_iter: config.max_iter,
        tol: config.tol,
        rr_b_weight: true,
        selection: Selection::LargestAbsK,
        recover_rayleigh: true,
    };
    lobpcg_driver(a, b, &spec, &op, &op)
}

// ─────────────────────────────────────────────────────────────────────────────
// Driver
// ─────────────────────────────────────────────────────────────────────────────

/// How Ritz pairs are selected from the projected spectrum.
#[derive(Clone, Copy)]
enum Selection {
    /// The `k` algebraically smallest values (standard mode).
    SmallestK,
    /// The `k` values of largest magnitude (shift-invert `μ`).
    LargestAbsK,
}

/// Mode-independent driver parameters (bundled to keep the arity small).
struct DriverSpec {
    k: usize,
    max_iter: usize,
    tol: f64,
    /// Weight the projected operator with `B q_i` on the left (interior mode,
    /// where `OP` is `B`-self-adjoint) instead of `q_i`.
    rr_b_weight: bool,
    selection: Selection,
    /// Recover `λ` per Ritz vector by Rayleigh quotient (interior mode) rather
    /// than taking the projected eigenvalue directly (standard mode).
    recover_rayleigh: bool,
}

/// Source block of a search-space column (used to build the `P` update).
#[derive(Clone, Copy, PartialEq, Eq)]
enum SourceTag {
    /// Current iterate block `X`.
    X,
    /// Preconditioned-residual block `W`.
    W,
    /// Conjugate-direction block `P`.
    P,
}

/// A `B`-orthonormal column together with its `B`-image and source tag.
struct OrthoCol {
    /// The column `q` (length `n`), `B`-normalized.
    q: Vec<f64>,
    /// `B q` (length `n`).
    bq: Vec<f64>,
    /// Originating block.
    tag: SourceTag,
}

/// Core LOBPCG iteration shared by all modes.
fn lobpcg_driver(
    a: &CsrMatrix,
    b: Option<&CsrMatrix>,
    spec: &DriverSpec,
    apply_op: impl Fn(&[f64]) -> Vec<f64>,
    apply_precond: impl Fn(&[f64]) -> Vec<f64>,
) -> Result<LobpcgResult, EigensolverError> {
    let n = a.nrows;
    let k = spec.k;

    // Deterministic, generic initial block: fixed-seed StdRng in [-1, 1).
    // (Not the analytic eigenvectors — a generic random start.)
    let mut rng = StdRng::seed_from_u64(0x0010_B9C6);
    let init_cols: Vec<Vec<f64>> = (0..k)
        .map(|_| (0..n).map(|_| rng.random_range(-1.0..1.0)).collect())
        .collect();
    let init_tags = vec![SourceTag::X; k];
    let init_kept = b_orthonormalize(&init_cols, &init_tags, b);
    if init_kept.len() < k {
        return Err(EigensolverError::General(
            "LOBPCG: failed to construct a full-rank initial block".into(),
        ));
    }
    let mut x_cols: Vec<Vec<f64>> = init_kept.into_iter().map(|c| c.q).collect();
    let mut lambda: Vec<f64> = x_cols.iter().map(|x| rayleigh_quotient(a, b, x)).collect();
    let mut p_cols: Vec<Vec<f64>> = Vec::new();

    let mut iterations = 0usize;
    let mut converged = false;

    loop {
        // ── 1. Residuals + active (unconverged) set. ────────────────────────
        let mut residuals: Vec<Vec<f64>> = Vec::with_capacity(k);
        let mut active: Vec<usize> = Vec::new();
        let mut all_converged = true;
        for (s, x) in x_cols.iter().enumerate() {
            let ax = spmv_vec(a, x);
            let bx = apply_b_vec(b, x);
            let r: Vec<f64> = ax
                .iter()
                .zip(bx.iter())
                .map(|(av, bv)| av - lambda[s] * bv)
                .collect();
            let rel = norm2(&r) / (lambda[s].abs() + 1.0);
            if rel >= spec.tol {
                all_converged = false;
                active.push(s);
            }
            residuals.push(r);
        }
        if all_converged {
            converged = true;
            break;
        }
        if iterations >= spec.max_iter {
            break;
        }

        // ── 2. Assemble S = [X (locked+active) | W (active only) | P]. ───────
        let mut s_cols: Vec<Vec<f64>> = Vec::new();
        let mut s_tags: Vec<SourceTag> = Vec::new();
        for x in &x_cols {
            s_cols.push(x.clone());
            s_tags.push(SourceTag::X);
        }
        for &s in &active {
            // Soft locking: only active columns spawn a W direction.
            let w = apply_precond(&residuals[s]);
            if w.iter().all(|v| v.is_finite()) {
                s_cols.push(w);
                s_tags.push(SourceTag::W);
            }
        }
        for p in &p_cols {
            s_cols.push(p.clone());
            s_tags.push(SourceTag::P);
        }

        // ── 3. B-orthonormalize the search space. ───────────────────────────
        let kept = b_orthonormalize(&s_cols, &s_tags, b);
        let m = kept.len();
        if m < k {
            // Rank collapsed below the block size — stop with the current best.
            break;
        }

        // ── 4. Rayleigh–Ritz on the (now standard) projected problem. ───────
        // op_q[j] = OP q_j (== A q_j in standard mode).
        let op_q: Vec<Vec<f64>> = kept.iter().map(|c| apply_op(&c.q)).collect();
        let mut proj = vec![0.0f64; m * m];
        for (i, ci) in kept.iter().enumerate() {
            let left = if spec.rr_b_weight { &ci.bq } else { &ci.q };
            for (j, oqj) in op_q.iter().enumerate() {
                proj[i * m + j] = dot(left, oqj);
            }
        }
        // Symmetrize to remove round-off asymmetry.
        for i in 0..m {
            for j in (i + 1)..m {
                let avg = 0.5 * (proj[i * m + j] + proj[j * m + i]);
                proj[i * m + j] = avg;
                proj[j * m + i] = avg;
            }
        }
        let (theta, evecs) = jacobi_eigensolve(&proj, m);

        // ── 5. Select k Ritz pairs. ─────────────────────────────────────────
        let sel: Vec<usize> = match spec.selection {
            Selection::SmallestK => (0..k).collect(),
            Selection::LargestAbsK => {
                let mut idx: Vec<usize> = (0..m).collect();
                idx.sort_by(|&i, &j| {
                    theta[j]
                        .abs()
                        .partial_cmp(&theta[i].abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                idx.truncate(k);
                idx
            }
        };

        // ── 6. Form X_new (= Q Y[:,sel]), P_new (W,P rows only), λ. ─────────
        let mut x_new: Vec<Vec<f64>> = Vec::with_capacity(k);
        let mut p_candidates: Vec<Vec<f64>> = Vec::with_capacity(k);
        let mut lambda_new: Vec<f64> = Vec::with_capacity(k);
        for &si in &sel {
            let y = &evecs[si];
            let mut xc = vec![0.0f64; n];
            let mut pc = vec![0.0f64; n];
            for (j, kc) in kept.iter().enumerate() {
                axpy(&mut xc, y[j], &kc.q);
                if matches!(kc.tag, SourceTag::W | SourceTag::P) {
                    axpy(&mut pc, y[j], &kc.q);
                }
            }
            let lam = if spec.recover_rayleigh {
                rayleigh_quotient(a, b, &xc)
            } else {
                theta[si]
            };
            x_new.push(xc);
            p_candidates.push(pc);
            lambda_new.push(lam);
        }

        // ── 7. B-orthonormalize P_new against X_new (drop collapsed columns). ─
        let mut po_cols: Vec<Vec<f64>> = Vec::with_capacity(2 * k);
        let mut po_tags: Vec<SourceTag> = Vec::with_capacity(2 * k);
        for xc in &x_new {
            po_cols.push(xc.clone());
            po_tags.push(SourceTag::X);
        }
        for pc in &p_candidates {
            po_cols.push(pc.clone());
            po_tags.push(SourceTag::P);
        }
        let po_kept = b_orthonormalize(&po_cols, &po_tags, b);
        p_cols = po_kept
            .into_iter()
            .filter(|c| c.tag == SourceTag::P)
            .map(|c| c.q)
            .collect();

        x_cols = x_new;
        lambda = lambda_new;
        iterations += 1;
    }

    // Sort the result ascending by eigenvalue.
    let mut order: Vec<usize> = (0..x_cols.len()).collect();
    order.sort_by(|&i, &j| {
        lambda[i]
            .partial_cmp(&lambda[j])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let eigenvalues: Vec<f64> = order.iter().map(|&i| lambda[i]).collect();
    let eigenvectors: Vec<Vec<f64>> = order.iter().map(|&i| x_cols[i].clone()).collect();

    Ok(LobpcgResult {
        eigenvalues,
        eigenvectors,
        iterations,
        converged,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear-algebra helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Sparse matvec returning a fresh vector: `y = A x`.
fn spmv_vec(a: &CsrMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0f64; a.nrows];
    a.spmv_par(x, &mut y);
    y
}

/// Apply `B` (identity when `b` is `None`): returns `B x`.
fn apply_b_vec(b: Option<&CsrMatrix>, x: &[f64]) -> Vec<f64> {
    match b {
        Some(bm) => spmv_vec(bm, x),
        None => x.to_vec(),
    }
}

/// Euclidean inner product.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean 2-norm.
fn norm2(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// `y += alpha * x`.
fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

/// Rayleigh quotient `λ = (xᵀ A x) / (xᵀ B x)` (guarded denominator).
fn rayleigh_quotient(a: &CsrMatrix, b: Option<&CsrMatrix>, x: &[f64]) -> f64 {
    let ax = spmv_vec(a, x);
    let bx = apply_b_vec(b, x);
    let num = dot(x, &ax);
    let den = dot(x, &bx);
    if den.abs() > 1e-300 {
        num / den
    } else {
        0.0
    }
}

/// `B`-orthonormalize a list of columns by re-orthogonalizing modified
/// Gram–Schmidt in the `B`-inner product.
///
/// Columns are processed in the given order. For each incoming column the running
/// `B`-image `Bs` is updated linearly as projections are subtracted, so it always
/// equals `B s` without recomputing the matvec. A single re-orthogonalization is
/// applied when more than ~50% of the norm is cancelled, and rank-deficient
/// columns (norm collapsing below `1e-10`) are dropped. The surviving columns `Q`
/// satisfy `Qᵀ B Q = I`.
fn b_orthonormalize(cols: &[Vec<f64>], tags: &[SourceTag], b: Option<&CsrMatrix>) -> Vec<OrthoCol> {
    let mut kept: Vec<OrthoCol> = Vec::with_capacity(cols.len());
    for (col, &tag) in cols.iter().zip(tags.iter()) {
        let mut s = col.clone();
        let mut bs = apply_b_vec(b, &s);
        let pre_norm = dot(&s, &bs).max(0.0).sqrt();

        // First projection pass.
        for kc in &kept {
            let proj = dot(&kc.q, &bs);
            axpy(&mut s, -proj, &kc.q);
            axpy(&mut bs, -proj, &kc.bq);
        }
        let mut new_norm = dot(&s, &bs).max(0.0).sqrt();

        // Re-orthogonalize once on heavy cancellation.
        if new_norm < 0.7 * pre_norm {
            for kc in &kept {
                let proj = dot(&kc.q, &bs);
                axpy(&mut s, -proj, &kc.q);
                axpy(&mut bs, -proj, &kc.bq);
            }
            new_norm = dot(&s, &bs).max(0.0).sqrt();
        }

        // Drop rank-deficient columns.
        if new_norm < 1e-10 * pre_norm.max(1.0) {
            continue;
        }
        let inv = 1.0 / new_norm.max(1e-300);
        for v in s.iter_mut() {
            *v *= inv;
        }
        for v in bs.iter_mut() {
            *v *= inv;
        }
        kept.push(OrthoCol { q: s, bq: bs, tag });
    }
    kept
}

/// Dense symmetric eigensolver by **cyclic Jacobi** rotations.
///
/// Takes a row-major `m × m` symmetric matrix and returns `(eigenvalues,
/// eigenvectors)` sorted ascending by eigenvalue, where `eigenvectors[s]` is the
/// `s`-th eigenvector (length `m`). Sweeps run until the off-diagonal Frobenius
/// norm drops below `1e-14 ·‖A‖_F` (guarded) or a sweep cap is reached.
fn jacobi_eigensolve(input: &[f64], m: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    if m == 0 {
        return (Vec::new(), Vec::new());
    }
    if m == 1 {
        return (vec![input[0]], vec![vec![1.0]]);
    }

    let mut a = input.to_vec();
    let mut v = vec![0.0f64; m * m];
    for i in 0..m {
        v[i * m + i] = 1.0;
    }

    let frob = a.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-300);
    let threshold = 1e-14 * frob;
    let max_sweeps = 100;

    for _sweep in 0..max_sweeps {
        // Off-diagonal Frobenius norm.
        let mut off = 0.0f64;
        for p in 0..m {
            for q in (p + 1)..m {
                off += a[p * m + q] * a[p * m + q];
            }
        }
        if off.sqrt() < threshold {
            break;
        }

        for p in 0..m {
            for q in (p + 1)..m {
                let apq = a[p * m + q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let app = a[p * m + p];
                let aqq = a[q * m + q];
                // Tangent of the rotation that zeros a[p][q]. sign(0) := +1 gives
                // the correct 45° rotation when app == aqq, with no float `==`.
                let phi = (aqq - app) / (2.0 * apq);
                let sgn = if phi >= 0.0 { 1.0 } else { -1.0 };
                let t = sgn / (phi.abs() + (phi * phi + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                let tau = s / (1.0 + c);

                a[p * m + p] = app - t * apq;
                a[q * m + q] = aqq + t * apq;
                a[p * m + q] = 0.0;
                a[q * m + p] = 0.0;

                for i in 0..m {
                    if i != p && i != q {
                        let aip = a[i * m + p];
                        let aiq = a[i * m + q];
                        let new_ip = aip - s * (aiq + tau * aip);
                        let new_iq = aiq + s * (aip - tau * aiq);
                        a[i * m + p] = new_ip;
                        a[p * m + i] = new_ip;
                        a[i * m + q] = new_iq;
                        a[q * m + i] = new_iq;
                    }
                }
                for i in 0..m {
                    let vip = v[i * m + p];
                    let viq = v[i * m + q];
                    v[i * m + p] = vip - s * (viq + tau * vip);
                    v[i * m + q] = viq + s * (vip - tau * viq);
                }
            }
        }
    }

    // Eigenvalues are the diagonal; sort ascending and permute eigenvectors.
    let mut pairs: Vec<(f64, usize)> = (0..m).map(|i| (a[i * m + i], i)).collect();
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let eigenvalues: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let eigenvectors: Vec<Vec<f64>> = pairs
        .iter()
        .map(|&(_, col)| (0..m).map(|i| v[i * m + col]).collect())
        .collect();
    (eigenvalues, eigenvectors)
}

/// Assemble `A − σ B` as an explicit CSR matrix over the union sparsity pattern.
///
/// For `b = None` this is `A` with its diagonal shifted by `−σ`, inserting a
/// diagonal entry where one is absent. Duplicate entries within a row are summed.
fn assemble_a_minus_sigma_b(a: &CsrMatrix, b: Option<&CsrMatrix>, sigma: f64) -> CsrMatrix {
    let n = a.nrows;
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices: Vec<usize> = Vec::with_capacity(a.col_indices.len() + n);
    let mut values: Vec<f64> = Vec::with_capacity(a.values.len() + n);

    for i in 0..n {
        let mut row: BTreeMap<usize, f64> = BTreeMap::new();
        for kk in a.row_offsets[i]..a.row_offsets[i + 1] {
            *row.entry(a.col_indices[kk]).or_insert(0.0) += a.values[kk];
        }
        match b {
            Some(bm) => {
                for kk in bm.row_offsets[i]..bm.row_offsets[i + 1] {
                    *row.entry(bm.col_indices[kk]).or_insert(0.0) -= sigma * bm.values[kk];
                }
            }
            None => {
                // B = I: shift the diagonal, creating it if necessary.
                *row.entry(i).or_insert(0.0) -= sigma;
            }
        }
        for (c, val) in row {
            col_indices.push(c);
            values.push(val);
        }
        row_offsets[i + 1] = col_indices.len();
    }

    CsrMatrix {
        nrows: n,
        ncols: a.ncols,
        row_offsets,
        col_indices,
        values,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests for the private numerical core
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_diagonal_matrix_sorts_ascending() {
        // diag(3, 1, 2) -> eigenvalues 1, 2, 3.
        let a = vec![3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0];
        let (ev, _vecs) = jacobi_eigensolve(&a, 3);
        assert!((ev[0] - 1.0).abs() < 1e-12, "ev0 = {}", ev[0]);
        assert!((ev[1] - 2.0).abs() < 1e-12, "ev1 = {}", ev[1]);
        assert!((ev[2] - 3.0).abs() < 1e-12, "ev2 = {}", ev[2]);
    }

    #[test]
    fn jacobi_symmetric_2x2_eigenpairs() {
        // [[2,1],[1,2]] -> eigenvalues 1 and 3.
        let a = vec![2.0, 1.0, 1.0, 2.0];
        let (ev, vecs) = jacobi_eigensolve(&a, 2);
        assert!((ev[0] - 1.0).abs() < 1e-12, "ev0 = {}", ev[0]);
        assert!((ev[1] - 3.0).abs() < 1e-12, "ev1 = {}", ev[1]);
        // Verify A v0 = ev0 v0 and unit norm.
        let v0 = &vecs[0];
        let av0 = [2.0 * v0[0] + v0[1], v0[0] + 2.0 * v0[1]];
        assert!((av0[0] - ev[0] * v0[0]).abs() < 1e-10);
        assert!((av0[1] - ev[0] * v0[1]).abs() < 1e-10);
        let nrm = (v0[0] * v0[0] + v0[1] * v0[1]).sqrt();
        assert!((nrm - 1.0).abs() < 1e-12, "norm = {nrm}");
    }

    #[test]
    fn shift_assembly_shifts_diagonal() {
        // 2x2 identity, shift by 0.5 -> diag becomes 0.5.
        let a = CsrMatrix::identity(2);
        let shifted = assemble_a_minus_sigma_b(&a, None, 0.5);
        let mut d = vec![0.0f64; 2];
        for i in 0..2 {
            for kk in shifted.row_offsets[i]..shifted.row_offsets[i + 1] {
                if shifted.col_indices[kk] == i {
                    d[i] = shifted.values[kk];
                }
            }
        }
        assert!((d[0] - 0.5).abs() < 1e-14);
        assert!((d[1] - 0.5).abs() < 1e-14);
    }
}
