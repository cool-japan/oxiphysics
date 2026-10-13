//! Krylov-subspace solvers and Ritz values for sparse symmetric systems.
//!
//! * [`minres_solve`] - MINRES (Paige & Saunders, 1975) for symmetric,
//!   possibly indefinite, systems.
//! * [`saddle_point_solve`] - Schur-complement conjugate gradients for
//!   saddle-point systems.
//! * [`ritz_values`] - eigenvalues of a Lanczos tridiagonal matrix.

use super::types::CsrMatrix;
use crate::error::PhysicsError;
use nalgebra::{DMatrix, SymmetricEigen};

/// Maximum number of times [`minres_solve`] restarts from the true residual
/// when the recurrence estimate claims convergence but the true residual
/// `b - A x` does not yet satisfy the tolerance.
const MINRES_MAX_REFINEMENTS: usize = 3;

/// Euclidean norm of a slice.
fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Euclidean inner product of two equally long slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Outcome of one MINRES cycle started from the zero vector.
struct MinresCycle {
    /// Approximate solution of `A x = rhs`.
    x: Vec<f64>,
    /// Number of Lanczos steps performed.
    iterations: usize,
}

/// Run MINRES (Paige & Saunders; SOL `minres.m`) on `A x = rhs` from `x = 0`.
///
/// The cycle stops when the recurrence estimate `|phibar|` of the residual
/// norm is at most `target` (absolute), when the Lanczos process breaks down
/// (invariant Krylov subspace, or a singular tridiagonal matrix), or after
/// `max_iter` Lanczos steps.  `rhs` must be non-zero and finite.
fn minres_cycle(a: &CsrMatrix, rhs: &[f64], target: f64, max_iter: usize) -> MinresCycle {
    let n = rhs.len();
    let beta1 = norm2(rhs);
    let mut x = vec![0.0f64; n];
    // Lanczos vectors v_k and v_{k-1}; beta couples them (beta_k, zero at the start).
    let mut v: Vec<f64> = rhs.iter().map(|&r| r / beta1).collect();
    let mut v_prev = vec![0.0f64; n];
    let mut beta = 0.0f64;
    // Running estimate of ||A|| from the Lanczos coefficients.
    let mut matrix_scale = 0.0f64;
    // Search directions w_{k-1} and w_{k-2}.
    let mut w = vec![0.0f64; n];
    let mut w_prev = vec![0.0f64; n];
    // State of the QR factorisation of the Lanczos tridiagonal matrix.
    let mut cs = -1.0f64;
    let mut sn = 0.0f64;
    let mut dbar = 0.0f64;
    let mut eps = 0.0f64;
    let mut phibar = beta1;
    let mut iterations = 0usize;
    while iterations < max_iter {
        iterations += 1;
        // Lanczos step: beta_{k+1} v_{k+1} = A v_k - alpha_k v_k - beta_k v_{k-1}.
        let av = a.matvec(&v);
        let alpha = dot(&v, &av);
        let mut v_next: Vec<f64> = av
            .iter()
            .zip(v.iter())
            .zip(v_prev.iter())
            .map(|((&avi, &vi), &vpi)| avi - alpha * vi - beta * vpi)
            .collect();
        let beta_next = norm2(&v_next);
        // Apply the previous plane rotation to the new tridiagonal column ...
        let old_eps = eps;
        let delta = cs * dbar + sn * alpha;
        let gbar = sn * dbar - cs * alpha;
        eps = sn * beta_next;
        dbar = -cs * beta_next;
        // ... and compute the next rotation.  A vanishing `gamma` means that the
        // tridiagonal matrix is numerically singular and the Lanczos process has
        // broken down at the same time (the right-hand side is outside the range
        // of a singular `A`): the previous iterate is the best one available.
        let gamma = gbar.hypot(beta_next);
        matrix_scale = matrix_scale.max(alpha.abs() + beta + beta_next);
        if !gamma.is_finite() || gamma <= 10.0 * f64::EPSILON * matrix_scale {
            break;
        }
        cs = gbar / gamma;
        sn = beta_next / gamma;
        let phi = cs * phibar;
        phibar *= sn;
        // Update the search direction and the iterate.
        let w_next: Vec<f64> = v
            .iter()
            .zip(w_prev.iter())
            .zip(w.iter())
            .map(|((&vi, &wpi), &wi)| (vi - old_eps * wpi - delta * wi) / gamma)
            .collect();
        for (xi, &wi) in x.iter_mut().zip(w_next.iter()) {
            *xi += phi * wi;
        }
        w_prev = std::mem::replace(&mut w, w_next);
        // The step above is complete even if the Lanczos process ends here, so
        // the stopping tests come after the update of `x`.
        if phibar.abs() <= target {
            break;
        }
        if beta_next <= f64::EPSILON * (alpha.abs() + beta) {
            break;
        }
        let inv_beta = 1.0 / beta_next;
        for vi in v_next.iter_mut() {
            *vi *= inv_beta;
        }
        v_prev = std::mem::replace(&mut v, v_next);
        beta = beta_next;
    }
    MinresCycle { x, iterations }
}

/// Minimum Residual (MINRES) method for symmetric (possibly indefinite) systems `A x = b`.
///
/// Implements the Paige-Saunders recurrences of SOL's `minres.m`: a Lanczos
/// process produces a tridiagonal matrix whose QR factorisation (Givens
/// rotations) yields the residual-minimising iterate in the Krylov space.
///
/// `tol` is a relative tolerance: the method is converged when
/// `||b - A x|| <= tol * ||b||` (a negative or NaN `tol` counts as zero).
/// The recurrence only supplies an *estimate* of the residual, so the true
/// residual `b - A x` is evaluated (with one extra matrix-vector product)
/// before `converged` is reported.  If the estimate and the true residual
/// disagree, the iteration is restarted from the true residual a few times,
/// which is classical iterative refinement.
///
/// Returns `(solution, converged, iterations)`, where `converged` reflects
/// the **true** residual and `iterations` counts Lanczos steps over all
/// restarts (never more than `max_iter`).
///
/// A zero right-hand side returns the zero vector at zero iterations.  For a
/// singular system whose right-hand side is outside the range of `A` the true
/// residual cannot reach the tolerance and `converged` is `false`.
///
/// # Panics
/// Panics if `b.len() != a.nrows`.
pub fn minres_solve(
    a: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, bool, usize) {
    let n = a.nrows;
    assert_eq!(b.len(), n);
    let mut x = vec![0.0f64; n];
    let b_norm = norm2(b);
    if b_norm == 0.0 {
        return (x, true, 0);
    }
    if !b_norm.is_finite() {
        return (x, false, 0);
    }
    // A negative or NaN tolerance cannot be met by any iterate: treat it as zero.
    let target = tol.max(0.0) * b_norm;
    if b_norm <= target {
        // The zero vector already satisfies the relative tolerance.
        return (x, true, 0);
    }
    let mut rhs = b.to_vec();
    let mut previous_residual = b_norm;
    let mut total_iterations = 0usize;
    for _ in 0..=MINRES_MAX_REFINEMENTS {
        let remaining = max_iter.saturating_sub(total_iterations);
        if remaining == 0 {
            break;
        }
        let cycle = minres_cycle(a, &rhs, target, remaining);
        total_iterations += cycle.iterations;
        for (xi, &di) in x.iter_mut().zip(cycle.x.iter()) {
            *xi += di;
        }
        let ax = a.matvec(&x);
        rhs = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let residual = norm2(&rhs);
        if residual <= target {
            return (x, true, total_iterations);
        }
        // Restart only while it can still help: stop on exhausted budget,
        // non-finite data, or when the last cycle brought no improvement.
        if !residual.is_finite() || residual >= previous_residual || total_iterations >= max_iter {
            break;
        }
        previous_residual = residual;
    }
    (x, false, total_iterations)
}

/// Build an [`PhysicsError::InvalidInput`] for the named argument.
fn invalid_input(field: &str, reason: String) -> PhysicsError {
    PhysicsError::InvalidInput {
        field: field.to_string(),
        reason,
    }
}

/// Conjugate gradients for a symmetric positive-definite operator given as a
/// fallible closure `v -> S v`.
///
/// Converges when the residual norm is at most `rel_tol * ||rhs||` and reports
/// [`PhysicsError::ConvergenceFailed`] when the iteration budget is exhausted,
/// the operator loses positive definiteness (`p^T S p <= 0`), or the residual
/// becomes non-finite.
pub(crate) fn cg_with_operator<F>(
    mut apply: F,
    rhs: &[f64],
    rel_tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, PhysicsError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, PhysicsError>,
{
    let mut x = vec![0.0f64; rhs.len()];
    let rhs_norm = norm2(rhs);
    if rhs_norm == 0.0 {
        return Ok(x);
    }
    let target = rel_tol * rhs_norm;
    let mut r = rhs.to_vec();
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    let mut iterations = 0usize;
    loop {
        let residual = rr.sqrt();
        if !residual.is_finite() {
            return Err(PhysicsError::ConvergenceFailed {
                iterations,
                tolerance: target,
                residual,
            });
        }
        if residual <= target {
            return Ok(x);
        }
        if iterations >= max_iter {
            return Err(PhysicsError::ConvergenceFailed {
                iterations,
                tolerance: target,
                residual,
            });
        }
        let ap = apply(&p)?;
        let pap = dot(&p, &ap);
        if pap <= 0.0 || !pap.is_finite() {
            return Err(PhysicsError::ConvergenceFailed {
                iterations,
                tolerance: target,
                residual,
            });
        }
        let alpha = rr / pap;
        for ((xi, ri), (&pi, &api)) in x.iter_mut().zip(r.iter_mut()).zip(p.iter().zip(ap.iter())) {
            *xi += alpha * pi;
            *ri -= alpha * api;
        }
        let rr_new = dot(&r, &r);
        let beta = rr_new / rr;
        for (pi, &ri) in p.iter_mut().zip(r.iter()) {
            *pi = ri + beta * *pi;
        }
        rr = rr_new;
        iterations += 1;
    }
}

/// Solve the saddle-point system
///
/// ```text
/// [ A   B^T ] [ x ]   [ f ]
/// [ B   0   ] [ p ] = [ g ]
/// ```
///
/// with a Schur-complement approach.  Eliminating `x = A^{-1} (f - B^T p)`
/// gives the pressure equation `S p = B A^{-1} f - g` with
/// `S = B A^{-1} B^T`, which is symmetric positive definite when `A` is SPD
/// and `B` has full row rank.  It is solved by matrix-free conjugate
/// gradients in which every application of `S` performs an inner CG solve
/// with `A`; afterwards `x` is recovered with one more inner solve.
///
/// # Arguments
/// * `a`        - (n x n) SPD matrix
/// * `bt`       - (n x m) matrix B^T
/// * `b`        - (m x n) matrix B
/// * `f`        - right-hand side of length n
/// * `g`        - right-hand side of length m
/// * `tol`      - relative tolerance of the pressure CG iteration; the inner
///   solves with `A` use a tighter tolerance (`tol / 1000`, but never below
///   `1e-14`)
/// * `max_iter` - iteration budget; the inner and outer solves are always
///   allowed at least `2 n + 20` and `2 m + 20` iterations respectively
///
/// Returns `(x, p)`.
///
/// # Errors
/// * [`PhysicsError::InvalidInput`] if the dimensions are inconsistent or
///   `tol` is not a positive finite number.
/// * [`PhysicsError::ConvergenceFailed`] if an inner or the outer CG
///   iteration fails to converge (for example a rank-deficient `B` with an
///   inconsistent `g`, or an indefinite `A`), or if the final residuals
///   `||f - A x - B^T p||` and `||g - B x||` exceed `10 * tol` relative to the
///   magnitude of the terms involved.
pub fn saddle_point_solve(
    a: &CsrMatrix,
    bt: &CsrMatrix,
    b: &CsrMatrix,
    f: &[f64],
    g: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, Vec<f64>), PhysicsError> {
    let n = a.nrows;
    let m = b.nrows;
    if !(tol.is_finite() && tol > 0.0) {
        return Err(invalid_input(
            "tol",
            format!("tolerance must be positive and finite, got {tol}"),
        ));
    }
    if a.ncols != n {
        return Err(invalid_input(
            "a",
            format!("matrix must be square, got {n} x {}", a.ncols),
        ));
    }
    if b.ncols != n {
        return Err(invalid_input(
            "b",
            format!("B must have {n} columns, got {}", b.ncols),
        ));
    }
    if bt.nrows != n || bt.ncols != m {
        return Err(invalid_input(
            "bt",
            format!("B^T must be {n} x {m}, got {} x {}", bt.nrows, bt.ncols),
        ));
    }
    if f.len() != n {
        return Err(invalid_input(
            "f",
            format!("expected length {n}, got {}", f.len()),
        ));
    }
    if g.len() != m {
        return Err(invalid_input(
            "g",
            format!("expected length {m}, got {}", g.len()),
        ));
    }
    let inner_tol = (tol * 1e-3).max(1e-14).min(tol);
    let inner_max_iter = max_iter.max(2 * n + 20);
    let outer_max_iter = max_iter.max(2 * m + 20);
    let solve_a = |rhs: &[f64]| -> Result<Vec<f64>, PhysicsError> {
        cg_with_operator(|v| Ok(a.matvec(v)), rhs, inner_tol, inner_max_iter)
    };
    // Pressure equation  S p = B A^{-1} f - g.
    let x_unconstrained = solve_a(f)?;
    let b_x_unconstrained = b.matvec(&x_unconstrained);
    let schur_rhs: Vec<f64> = b_x_unconstrained
        .iter()
        .zip(g.iter())
        .map(|(bx, gi)| bx - gi)
        .collect();
    let p = cg_with_operator(
        |q| {
            let bt_q = bt.matvec(q);
            let a_inv_bt_q = solve_a(&bt_q)?;
            Ok(b.matvec(&a_inv_bt_q))
        },
        &schur_rhs,
        tol,
        outer_max_iter,
    )?;
    // Velocity  x = A^{-1} (f - B^T p).
    let bt_p = bt.matvec(&p);
    let momentum_rhs: Vec<f64> = f.iter().zip(bt_p.iter()).map(|(fi, bi)| fi - bi).collect();
    let x = solve_a(&momentum_rhs)?;
    // Verify both block rows of the system.
    let ax = a.matvec(&x);
    let momentum_residual: Vec<f64> = f
        .iter()
        .zip(ax.iter())
        .zip(bt_p.iter())
        .map(|((fi, axi), bpi)| fi - axi - bpi)
        .collect();
    let bx = b.matvec(&x);
    let constraint_residual: Vec<f64> = g.iter().zip(bx.iter()).map(|(gi, bxi)| gi - bxi).collect();
    let relative = |residual: f64, scale: f64| {
        if scale > 0.0 {
            residual / scale
        } else {
            residual
        }
    };
    let momentum_error = relative(norm2(&momentum_residual), norm2(f) + norm2(&bt_p));
    let constraint_error = relative(
        norm2(&constraint_residual),
        norm2(g) + norm2(&b_x_unconstrained),
    );
    let worst = momentum_error.max(constraint_error);
    if worst > 10.0 * tol || worst.is_nan() {
        return Err(PhysicsError::ConvergenceFailed {
            iterations: outer_max_iter,
            tolerance: tol,
            residual: worst,
        });
    }
    Ok((x, p))
}

/// Compute the Ritz values of a Lanczos tridiagonal matrix.
///
/// Given `alphas` (diagonal) and `betas` (sub-diagonal; entries beyond the
/// first `alphas.len() - 1` are ignored and missing entries count as zero),
/// returns the eigenvalues of the symmetric tridiagonal matrix in ascending
/// order, obtained with the symmetric eigensolver of `nalgebra`.
///
/// If any input is non-finite, or the eigensolver does not converge, a vector
/// of `NaN` of the same length is returned.
pub fn ritz_values(alphas: &[f64], betas: &[f64]) -> Vec<f64> {
    let m = alphas.len();
    if m == 0 {
        return Vec::new();
    }
    let off = &betas[..betas.len().min(m - 1)];
    if !alphas.iter().chain(off.iter()).all(|v| v.is_finite()) {
        return vec![f64::NAN; m];
    }
    let mut t = DMatrix::<f64>::zeros(m, m);
    for (i, &alpha) in alphas.iter().enumerate() {
        t[(i, i)] = alpha;
    }
    for (i, &beta) in off.iter().enumerate() {
        t[(i + 1, i)] = beta;
        t[(i, i + 1)] = beta;
    }
    match SymmetricEigen::try_new(t, f64::EPSILON, 64 * m + 64) {
        Some(eigen) => {
            let mut values: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
            values.sort_by(|x, y| x.total_cmp(y));
            values
        }
        None => vec![f64::NAN; m],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::{add_identity_scaled, identity_csr, lanczos, laplacian_1d};
    use nalgebra::DVector;
    use std::f64::consts::PI;

    fn diagonal(values: &[f64]) -> CsrMatrix {
        let idx: Vec<usize> = (0..values.len()).collect();
        CsrMatrix::from_triplets(values.len(), values.len(), &idx, &idx, values)
    }

    /// Symmetric tridiagonal matrix with the given diagonal and off-diagonal.
    fn tridiagonal(diag: &[f64], off: &[f64]) -> CsrMatrix {
        let n = diag.len();
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for (i, &d) in diag.iter().enumerate() {
            rows.push(i);
            cols.push(i);
            vals.push(d);
        }
        for (i, &o) in off.iter().enumerate() {
            rows.push(i);
            cols.push(i + 1);
            vals.push(o);
            rows.push(i + 1);
            cols.push(i);
            vals.push(o);
        }
        CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
    }

    fn relative_residual(a: &CsrMatrix, x: &[f64], b: &[f64]) -> f64 {
        let ax = a.matvec(x);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        norm2(&r) / norm2(b)
    }

    fn dense_solve(a: &CsrMatrix, b: &[f64]) -> Vec<f64> {
        let n = a.nrows;
        let flat: Vec<f64> = a.to_dense().into_iter().flatten().collect();
        let dense = DMatrix::from_row_slice(n, n, &flat);
        let sol = dense.lu().solve(&DVector::from_row_slice(b)).unwrap();
        sol.iter().copied().collect()
    }

    // ── MINRES ───────────────────────────────────────────────────────────────

    #[test]
    fn minres_spd_diagonal_reaches_exact_solution() {
        // Regression: the old recurrence "converged" here with a relative
        // residual of 0.65.
        let a = diagonal(&[2.0, 3.0, 5.0]);
        let b = [1.0, 1.0, 1.0];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-10, 50);
        assert!(converged);
        assert!((1..=3).contains(&iterations), "iterations = {iterations}");
        for (xi, expected) in x.iter().zip([0.5, 1.0 / 3.0, 0.2]) {
            assert!((xi - expected).abs() < 1e-10, "x = {x:?}");
        }
        assert!(relative_residual(&a, &x, &b) < 1e-10);
    }

    #[test]
    fn minres_indefinite_tridiagonal_has_small_true_residual() {
        // Path-graph adjacency matrix: eigenvalues +-1.618 and +-0.618.
        let a = tridiagonal(&[0.0; 4], &[1.0; 3]);
        let b = [1.0, 2.0, 3.0, 4.0];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-10, 100);
        assert!(converged, "iterations = {iterations}");
        assert!(relative_residual(&a, &x, &b) < 1e-10);
        let reference = dense_solve(&a, &b);
        for (xi, ri) in x.iter().zip(reference.iter()) {
            assert!(
                (xi - ri).abs() < 1e-9,
                "x = {x:?}, reference = {reference:?}"
            );
        }
    }

    #[test]
    fn minres_mixed_sign_diagonal_converges() {
        let a = diagonal(&[1.0, -2.0, 3.0, -4.0, 5.0]);
        let b = [1.0, 1.0, 1.0, 1.0, 1.0];
        let (x, converged, _) = minres_solve(&a, &b, 1e-12, 50);
        assert!(converged);
        for (xi, ai) in x.iter().zip([1.0, -2.0, 3.0, -4.0, 5.0]) {
            assert!((xi - 1.0 / ai).abs() < 1e-10, "x = {x:?}");
        }
    }

    #[test]
    fn minres_shifted_laplacian_many_iterations() {
        // Indefinite: the shift lies inside the spectrum (0, 4) of the 1-D Laplacian.
        let n = 50;
        let a = add_identity_scaled(&laplacian_1d(n, 1.0), -1.7);
        let b: Vec<f64> = (0..n).map(|i| ((i * 7 + 3) % 11) as f64 - 5.0).collect();
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-10, 2000);
        assert!(converged, "iterations = {iterations}");
        assert!(
            iterations > 10,
            "should need many Lanczos steps: {iterations}"
        );
        assert!(relative_residual(&a, &x, &b) < 1e-9);
        let reference = dense_solve(&a, &b);
        for (xi, ri) in x.iter().zip(reference.iter()) {
            assert!((xi - ri).abs() < 1e-6 * (1.0 + ri.abs()), "{xi} vs {ri}");
        }
    }

    #[test]
    fn minres_spd_laplacian_converges() {
        let n = 40;
        let a = laplacian_1d(n, 1.0);
        let b = vec![1.0; n];
        let (x, converged, _) = minres_solve(&a, &b, 1e-10, 1000);
        assert!(converged);
        assert!(relative_residual(&a, &x, &b) < 1e-9);
    }

    #[test]
    fn minres_eigenvector_rhs_breaks_down_after_one_step() {
        let a = diagonal(&[2.0, 3.0, 5.0]);
        let b = [1.0, 0.0, 0.0];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-12, 10);
        assert!(converged);
        assert_eq!(iterations, 1);
        assert!((x[0] - 0.5).abs() < 1e-14);
        assert!(x[1].abs() < 1e-14 && x[2].abs() < 1e-14);
    }

    #[test]
    fn minres_identity_converges_in_one_step() {
        let a = identity_csr(4);
        let b = [1.0, 2.0, 3.0, 4.0];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-8, 100);
        assert!(converged);
        assert_eq!(iterations, 1);
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < 1e-12);
        }
    }

    #[test]
    fn minres_reports_non_convergence_honestly() {
        let n = 50;
        let a = laplacian_1d(n, 1.0);
        let b = vec![1.0; n];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-12, 3);
        assert!(!converged);
        assert_eq!(iterations, 3);
        assert!(relative_residual(&a, &x, &b) > 1e-3);
    }

    #[test]
    fn minres_inconsistent_singular_system_is_not_converged() {
        // diag(1, 0) x = (1, 1) has no solution.
        let a = diagonal(&[1.0, 0.0]);
        let b = [1.0, 1.0];
        let (x, converged, _) = minres_solve(&a, &b, 1e-10, 50);
        assert!(!converged);
        assert!(x.iter().all(|v| v.is_finite()));
        assert!(relative_residual(&a, &x, &b) > 0.5);
    }

    #[test]
    fn minres_zero_rhs_and_trivial_tolerance() {
        let a = diagonal(&[2.0, 3.0]);
        let (x, converged, iterations) = minres_solve(&a, &[0.0, 0.0], 1e-10, 10);
        assert!(converged && iterations == 0);
        assert_eq!(x, vec![0.0, 0.0]);
        // A relative tolerance of one is satisfied by x = 0.
        let (x, converged, iterations) = minres_solve(&a, &[1.0, 1.0], 1.0, 10);
        assert!(converged && iterations == 0);
        assert_eq!(x, vec![0.0, 0.0]);
    }

    #[test]
    fn minres_negative_or_nan_tolerance_never_produces_nan_iterates() {
        let a = diagonal(&[2.0, 3.0, 5.0]);
        let b = [1.0, 1.0, 1.0];
        for tol in [-1.0, f64::NAN, 0.0] {
            let (x, _, iterations) = minres_solve(&a, &b, tol, 50);
            assert!(x.iter().all(|v| v.is_finite()), "tol = {tol}: x = {x:?}");
            assert!(iterations <= 50, "tol = {tol}: {iterations} iterations");
            assert!(relative_residual(&a, &x, &b) < 1e-10, "tol = {tol}");
        }
    }

    #[test]
    fn minres_non_finite_rhs_is_not_converged() {
        let a = diagonal(&[2.0, 3.0]);
        let (_, converged, iterations) = minres_solve(&a, &[f64::NAN, 1.0], 1e-8, 10);
        assert!(!converged);
        assert_eq!(iterations, 0);
    }

    /// Symmetric matrix `Q diag(lambda) Q^T` with a known spectrum, where `Q` is
    /// the orthogonal factor of a pseudo-random matrix, stored as CSR.
    fn symmetric_with_spectrum(lambda: &[f64], seed: u64) -> CsrMatrix {
        let n = lambda.len();
        let mut state = seed.wrapping_add(1);
        let mut next_uniform = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let random = DMatrix::from_fn(n, n, |_, _| next_uniform());
        let q = random.qr().q();
        let a = &q * DMatrix::from_diagonal(&DVector::from_row_slice(lambda)) * q.transpose();
        let rows: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| 0.5 * (a[(i, j)] + a[(j, i)])).collect())
            .collect();
        CsrMatrix::from_dense(&rows)
    }

    #[test]
    fn minres_dense_symmetric_indefinite_matrices_reach_the_dense_solution() {
        let n = 24;
        // Spectrum in [-4, -0.25] and [0.25, 4]: indefinite with condition number 16.
        let lambda: Vec<f64> = (0..n)
            .map(|k| {
                let magnitude = 0.25 + 3.75 * (k / 2) as f64 / (n / 2 - 1) as f64;
                if k % 2 == 0 { magnitude } else { -magnitude }
            })
            .collect();
        for seed in [11u64, 12, 13] {
            let a = symmetric_with_spectrum(&lambda, seed);
            let b: Vec<f64> = (0..n)
                .map(|i| ((i * 5 + seed as usize) % 9) as f64 - 4.0)
                .collect();
            let (x, converged, iterations) = minres_solve(&a, &b, 1e-11, 400);
            assert!(converged, "seed {seed}: {iterations} iterations");
            assert!(iterations <= 3 * n, "seed {seed}: {iterations} iterations");
            assert!(relative_residual(&a, &x, &b) < 1e-10, "seed {seed}");
            let reference = dense_solve(&a, &b);
            for (xi, ri) in x.iter().zip(reference.iter()) {
                assert!(
                    (xi - ri).abs() < 1e-8 * (1.0 + ri.abs()),
                    "seed {seed}: {xi} vs {ri}"
                );
            }
        }
    }

    #[test]
    fn minres_residual_is_non_increasing_in_the_iteration_budget() {
        // MINRES minimises ||b - A x|| over the k-th Krylov space, which grows with k.
        let n = 24;
        let lambda: Vec<f64> = (0..n)
            .map(|k| {
                let magnitude = 0.25 + 3.75 * (k / 2) as f64 / (n / 2 - 1) as f64;
                if k % 2 == 0 { magnitude } else { -magnitude }
            })
            .collect();
        let a = symmetric_with_spectrum(&lambda, 21);
        let b: Vec<f64> = (0..n).map(|i| ((i * 3 + 1) % 7) as f64 - 3.0).collect();
        let mut previous = f64::INFINITY;
        for budget in 1..=n {
            let (x, _, iterations) = minres_solve(&a, &b, 0.0, budget);
            assert!(iterations <= budget, "budget {budget}: {iterations}");
            let residual = relative_residual(&a, &x, &b);
            // Above the round-off floor the decrease must be monotone.
            if previous > 1e-9 {
                assert!(
                    residual <= previous * (1.0 + 1e-9),
                    "budget {budget}: residual {residual} after {previous}"
                );
            }
            previous = residual;
        }
        assert!(previous < 1e-8, "full Krylov space residual = {previous}");
    }

    #[test]
    fn minres_consistent_singular_system_converges() {
        // diag(2, 0, 3) x = (1, 0, 1) is solvable: x = (1/2, free, 1/3).
        let a = diagonal(&[2.0, 0.0, 3.0]);
        let b = [1.0, 0.0, 1.0];
        let (x, converged, iterations) = minres_solve(&a, &b, 1e-12, 50);
        assert!(converged, "iterations = {iterations}");
        assert!(
            (x[0] - 0.5).abs() < 1e-12 && (x[2] - 1.0 / 3.0).abs() < 1e-12,
            "x = {x:?}"
        );
        assert!(relative_residual(&a, &x, &b) < 1e-12);
    }

    // ── saddle_point_solve ───────────────────────────────────────────────────

    fn dense_csr(rows: &[&[f64]]) -> CsrMatrix {
        let dense: Vec<Vec<f64>> = rows.iter().map(|r| r.to_vec()).collect();
        CsrMatrix::from_dense(&dense)
    }

    /// Residual norms `(||f - A x - B^T p||, ||g - B x||)`.
    fn saddle_residuals(
        a: &CsrMatrix,
        b: &CsrMatrix,
        f: &[f64],
        g: &[f64],
        x: &[f64],
        p: &[f64],
    ) -> (f64, f64) {
        let ax = a.matvec(x);
        let btp = b.transpose().matvec(p);
        let r1: Vec<f64> = f
            .iter()
            .zip(ax.iter().zip(btp.iter()))
            .map(|(fi, (axi, bpi))| fi - axi - bpi)
            .collect();
        let bx = b.matvec(x);
        let r2: Vec<f64> = g.iter().zip(bx.iter()).map(|(gi, bxi)| gi - bxi).collect();
        (norm2(&r1), norm2(&r2))
    }

    #[test]
    fn saddle_point_one_by_one_has_the_correct_signs() {
        // x = 1, p = -1 solves [1 1; 1 0] [x; p] = [0; 1].
        let a = dense_csr(&[&[1.0]]);
        let b = dense_csr(&[&[1.0]]);
        let bt = b.transpose();
        let (x, p) = saddle_point_solve(&a, &bt, &b, &[0.0], &[1.0], 1e-12, 100).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10, "x = {x:?}");
        assert!((p[0] + 1.0).abs() < 1e-10, "p = {p:?}");
    }

    #[test]
    fn saddle_point_three_by_three_with_two_constraints() {
        let a = dense_csr(&[&[4.0, 1.0, 0.0], &[1.0, 3.0, 1.0], &[0.0, 1.0, 2.0]]);
        let b = dense_csr(&[&[1.0, 0.0, 1.0], &[0.0, 1.0, -1.0]]);
        let bt = b.transpose();
        let f = [1.0, 2.0, 3.0];
        let g = [0.5, -1.0];
        let (x, p) = saddle_point_solve(&a, &bt, &b, &f, &g, 1e-12, 100).unwrap();
        let (momentum, constraint) = saddle_residuals(&a, &b, &f, &g, &x, &p);
        assert!(momentum < 1e-8, "momentum residual = {momentum}");
        assert!(constraint < 1e-8, "constraint residual = {constraint}");
        // Compare with the dense solution of the full 5 x 5 KKT system.
        let mut kkt = vec![vec![0.0f64; 5]; 5];
        for (i, row) in a.to_dense().iter().enumerate() {
            kkt[i][..3].copy_from_slice(row);
        }
        for (j, row) in b.to_dense().iter().enumerate() {
            for (i, &v) in row.iter().enumerate() {
                kkt[i][3 + j] = v;
                kkt[3 + j][i] = v;
            }
        }
        let kkt_csr = CsrMatrix::from_dense(&kkt);
        let rhs = [1.0, 2.0, 3.0, 0.5, -1.0];
        let reference = dense_solve(&kkt_csr, &rhs);
        for (xi, ri) in x.iter().chain(p.iter()).zip(reference.iter()) {
            assert!((xi - ri).abs() < 1e-8, "{xi} vs {ri}");
        }
    }

    #[test]
    fn saddle_point_larger_system_satisfies_both_block_rows() {
        let n = 24;
        let a = add_identity_scaled(&laplacian_1d(n, 1.0), 0.5);
        // Five constraints: B x = (sum over a strided subset of the unknowns).
        let m = 5;
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for r in 0..m {
            for c in (r..n).step_by(m) {
                rows.push(r);
                cols.push(c);
                vals.push(1.0 + 0.1 * c as f64);
            }
        }
        let b = CsrMatrix::from_triplets(m, n, &rows, &cols, &vals);
        let bt = b.transpose();
        let f: Vec<f64> = (0..n).map(|i| (i as f64 * 0.7).sin()).collect();
        let g: Vec<f64> = (0..m).map(|i| 0.3 * i as f64 - 0.4).collect();
        let (x, p) = saddle_point_solve(&a, &bt, &b, &f, &g, 1e-12, 200).unwrap();
        let (momentum, constraint) = saddle_residuals(&a, &b, &f, &g, &x, &p);
        assert!(momentum < 1e-8, "momentum residual = {momentum}");
        assert!(constraint < 1e-8, "constraint residual = {constraint}");
    }

    #[test]
    fn saddle_point_without_constraints_reduces_to_a_solve() {
        let a = dense_csr(&[&[2.0, 1.0], &[1.0, 3.0]]);
        let b = CsrMatrix::new(0, 2);
        let bt = CsrMatrix::new(2, 0);
        let (x, p) = saddle_point_solve(&a, &bt, &b, &[1.0, 2.0], &[], 1e-12, 50).unwrap();
        assert!(p.is_empty());
        assert!(relative_residual(&a, &x, &[1.0, 2.0]) < 1e-10);
    }

    #[test]
    fn saddle_point_rejects_inconsistent_dimensions_and_tolerance() {
        let a = dense_csr(&[&[2.0, 1.0], &[1.0, 3.0]]);
        let b = dense_csr(&[&[1.0, 1.0]]);
        let bt = b.transpose();
        let invalid = |result: Result<(Vec<f64>, Vec<f64>), PhysicsError>| {
            matches!(result, Err(PhysicsError::InvalidInput { .. }))
        };
        assert!(invalid(saddle_point_solve(
            &a,
            &bt,
            &b,
            &[1.0],
            &[1.0],
            1e-8,
            10
        )));
        assert!(invalid(saddle_point_solve(
            &a,
            &bt,
            &b,
            &[1.0, 2.0],
            &[],
            1e-8,
            10
        )));
        assert!(invalid(saddle_point_solve(
            &a,
            &b,
            &b,
            &[1.0, 2.0],
            &[1.0],
            1e-8,
            10
        )));
        for tol in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(invalid(saddle_point_solve(
                &a,
                &bt,
                &b,
                &[1.0, 2.0],
                &[1.0],
                tol,
                10
            )));
        }
    }

    #[test]
    fn saddle_point_reports_convergence_failure_for_inconsistent_constraints() {
        // Two identical constraint rows with different right-hand sides.
        let a = dense_csr(&[&[2.0, 0.0], &[0.0, 3.0]]);
        let b = dense_csr(&[&[1.0, 1.0], &[1.0, 1.0]]);
        let bt = b.transpose();
        let result = saddle_point_solve(&a, &bt, &b, &[1.0, 1.0], &[1.0, 5.0], 1e-10, 50);
        assert!(
            matches!(result, Err(PhysicsError::ConvergenceFailed { .. })),
            "result = {result:?}"
        );
    }

    #[test]
    fn saddle_point_reports_failure_for_indefinite_a() {
        let a = diagonal(&[1.0, -1.0]);
        let b = dense_csr(&[&[1.0, 1.0]]);
        let bt = b.transpose();
        let result = saddle_point_solve(&a, &bt, &b, &[1.0, 1.0], &[0.0], 1e-10, 50);
        assert!(result.is_err(), "result = {result:?}");
    }

    // ── ritz_values ──────────────────────────────────────────────────────────

    #[test]
    fn ritz_values_of_small_tridiagonal_are_exact() {
        // Regression: the old QR sweep returned [3, 3, 1].
        let ritz = ritz_values(&[2.0, 2.0, 2.0], &[1.0, 1.0]);
        let s = 2.0f64.sqrt();
        let expected = [2.0 - s, 2.0, 2.0 + s];
        assert_eq!(ritz.len(), 3);
        for (r, e) in ritz.iter().zip(expected.iter()) {
            assert!((r - e).abs() < 1e-10, "ritz = {ritz:?}");
        }
    }

    #[test]
    fn ritz_values_of_toeplitz_tridiagonal_match_the_analytic_spectrum() {
        let m = 6;
        let ritz = ritz_values(&vec![2.0; m], &vec![-1.0; m - 1]);
        for (k, &r) in ritz.iter().enumerate() {
            let exact = 2.0 - 2.0 * ((k + 1) as f64 * PI / 7.0).cos();
            assert!((r - exact).abs() < 1e-10, "k = {k}: {r} vs {exact}");
        }
    }

    #[test]
    fn ritz_values_are_sorted_ascending_for_unsorted_diagonal() {
        let ritz = ritz_values(&[5.0, -1.0, 3.0, 0.0], &[0.0, 0.0, 0.0]);
        assert_eq!(ritz, vec![-1.0, 0.0, 3.0, 5.0]);
    }

    #[test]
    fn ritz_values_of_full_lanczos_run_match_the_laplacian_spectrum() {
        let n = 10;
        let a = laplacian_1d(n, 1.0);
        // A start vector that is not orthogonal to any eigenvector.
        let b0: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let (alphas, betas, _) = lanczos(&a, &b0, n);
        let ritz = ritz_values(&alphas, &betas);
        assert_eq!(ritz.len(), n);
        for (k, &r) in ritz.iter().enumerate() {
            let exact = 2.0 - 2.0 * ((k + 1) as f64 * PI / (n as f64 + 1.0)).cos();
            assert!((r - exact).abs() < 1e-6, "k = {k}: {r} vs {exact}");
        }
    }

    #[test]
    fn ritz_values_of_partial_lanczos_runs_approach_the_spectrum_from_inside() {
        let n = 30;
        let a = laplacian_1d(n, 1.0);
        let b0: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let lambda_min = 2.0 - 2.0 * (PI / (n as f64 + 1.0)).cos();
        let lambda_max = 2.0 - 2.0 * (n as f64 * PI / (n as f64 + 1.0)).cos();
        let mut previous = (f64::INFINITY, f64::INFINITY);
        for steps in [5usize, 10, 15, 20, 25] {
            let (alphas, betas, _) = lanczos(&a, &b0, steps);
            let ritz = ritz_values(&alphas, &betas);
            assert_eq!(ritz.len(), steps);
            let smallest = ritz[0];
            let largest = ritz[ritz.len() - 1];
            // Cauchy interlacing: the Ritz values lie inside the spectrum ...
            assert!(smallest >= lambda_min - 1e-9, "{steps} steps: {smallest}");
            assert!(largest <= lambda_max + 1e-9, "{steps} steps: {largest}");
            // ... and the extreme ones improve with every enlargement of the Krylov space.
            let errors = (smallest - lambda_min, lambda_max - largest);
            assert!(
                errors.0 < previous.0,
                "{steps} steps: lower error {errors:?}"
            );
            assert!(
                errors.1 < previous.1,
                "{steps} steps: upper error {errors:?}"
            );
            previous = errors;
        }
        // Reference run (NumPy): 1.1e-5 and 5.4e-3 after 25 of the 30 steps; the
        // upper end of the Laplacian spectrum is tightly clustered and converges slowly.
        assert!(previous.0 < 1e-4, "lower end error = {}", previous.0);
        assert!(previous.1 < 1e-2, "upper end error = {}", previous.1);
    }

    #[test]
    fn ritz_values_of_isolated_extreme_eigenvalues_converge_quickly() {
        // diag(-100, 1, ..., 9, 100): both extremes are well separated from the rest.
        let mut diag: Vec<f64> = vec![-100.0];
        diag.extend((1..=9).map(f64::from));
        diag.push(100.0);
        let a = diagonal(&diag);
        let b0 = vec![1.0; diag.len()];
        let (alphas, betas, _) = lanczos(&a, &b0, 8);
        let ritz = ritz_values(&alphas, &betas);
        assert_eq!(ritz.len(), 8);
        assert!((ritz[ritz.len() - 1] - 100.0).abs() < 1e-9, "{ritz:?}");
        assert!((ritz[0] + 100.0).abs() < 1e-9, "{ritz:?}");
    }

    #[test]
    fn ritz_values_degenerate_inputs() {
        assert!(ritz_values(&[], &[]).is_empty());
        assert_eq!(ritz_values(&[3.125], &[]), vec![3.125]);
        // Missing off-diagonals decouple the matrix.
        assert_eq!(ritz_values(&[2.0, 1.0], &[]), vec![1.0, 2.0]);
        let non_finite = ritz_values(&[1.0, f64::NAN], &[0.5]);
        assert_eq!(non_finite.len(), 2);
        assert!(non_finite.iter().all(|v| v.is_nan()));
    }
}
