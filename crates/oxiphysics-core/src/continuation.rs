// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical continuation and bifurcation analysis.
//!
//! Implements pseudo-arclength continuation for tracking solution branches of
//! parameterised nonlinear systems F(u, λ) = 0, together with fold and
//! bifurcation point detection, branch switching, turning-point location,
//! stability analysis, and generic path-following with event detection.

use nalgebra::{DMatrix, DVector, SVD};

use crate::numerical_linear_algebra::general_eigen::real_matrix_eigenvalues_rows;

/// Residual function F(u, λ) → ℝⁿ.
pub type ResidualFn = dyn Fn(&[f64], f64) -> Vec<f64>;
/// Jacobian function J(u, λ) → ℝⁿˣⁿ.
pub type JacobianFn = dyn Fn(&[f64], f64) -> Vec<Vec<f64>>;

// ---------------------------------------------------------------------------
// ContinuationState
// ---------------------------------------------------------------------------

/// State along a continuation branch.
///
/// Stores the current parameter value `lambda`, the solution vector `u`, the
/// tangent direction of the branch used for the arclength predictor, and the
/// accumulated arc length `s`.
#[derive(Debug, Clone)]
pub struct ContinuationState {
    /// Continuation parameter λ.
    pub lambda: f64,
    /// Solution vector u ∈ ℝⁿ.
    pub u: Vec<f64>,
    /// Unit tangent \[du/ds, dλ/ds\] of the branch at this state.
    ///
    /// States returned by [`PseudoArcLengthContinuation::step`] carry the true
    /// branch tangent (the normalised secant from the previous state in the
    /// rare case that the bordered system defining it is singular).
    /// [`ContinuationState::new`] starts with the pure-λ direction
    /// `[0, …, 0, 1]`, which is a placeholder rather than a measured tangent.
    pub tangent: Vec<f64>,
    /// Accumulated arc length along the branch.
    pub arc_length: f64,
    /// Desired step size Δs for the next predictor.
    pub ds: f64,
}

impl ContinuationState {
    /// Create a new continuation state.
    ///
    /// - `lambda` — initial parameter value.
    /// - `u` — initial solution vector.
    /// - `ds` — initial arc-length step size.
    pub fn new(lambda: f64, u: Vec<f64>, ds: f64) -> Self {
        let n = u.len();
        // Initial tangent: move purely in the λ direction.
        let mut tangent = vec![0.0; n + 1];
        tangent[n] = 1.0;
        Self {
            lambda,
            u,
            tangent,
            arc_length: 0.0,
            ds,
        }
    }

    /// Dimension of the solution vector.
    pub fn dim(&self) -> usize {
        self.u.len()
    }

    /// Return a normalised copy of the tangent vector.
    pub fn normalised_tangent(&self) -> Vec<f64> {
        let norm: f64 = self.tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-14 {
            self.tangent.iter().map(|x| x / norm).collect()
        } else {
            self.tangent.clone()
        }
    }
}

// ---------------------------------------------------------------------------
// pseudo_arclength_step — predictor
// ---------------------------------------------------------------------------

/// Advance the continuation state by one predictor step.
///
/// Uses the current tangent direction and step size `ds` to produce a
/// predicted state.  The resulting `arc_length` is incremented by `ds`.
pub fn pseudo_arclength_step(state: &ContinuationState) -> ContinuationState {
    let n = state.u.len();
    let ds = state.ds;
    let u_pred: Vec<f64> = state
        .u
        .iter()
        .enumerate()
        .map(|(i, &ui)| ui + ds * state.tangent[i])
        .collect();
    let lambda_pred = state.lambda + ds * state.tangent[n];
    ContinuationState {
        lambda: lambda_pred,
        u: u_pred,
        tangent: state.tangent.clone(),
        arc_length: state.arc_length + ds,
        ds,
    }
}

// ---------------------------------------------------------------------------
// CorrectorResult
// ---------------------------------------------------------------------------

/// Result of a Newton corrector iteration.
#[derive(Debug, Clone)]
pub struct CorrectorResult {
    /// Whether the Newton iteration converged.
    pub converged: bool,
    /// Number of Newton iterations performed.
    pub iterations: usize,
    /// Final residual ‖F(u, λ)‖.
    pub residual: f64,
    /// Corrected solution vector.
    pub u: Vec<f64>,
    /// Corrected parameter value.
    pub lambda: f64,
}

/// Newton corrector for pseudo-arclength continuation.
///
/// Given a predicted state and functions `f` (residual) and `jac` (Jacobian
/// w.r.t. `u`), iterates Newton's method augmented by the arclength
/// constraint until convergence or `max_iter` is reached.
pub fn corrector_newton(
    predicted: &ContinuationState,
    prev: &ContinuationState,
    f: &ResidualFn,
    jac: &JacobianFn,
    tol: f64,
    max_iter: usize,
) -> CorrectorResult {
    let n = predicted.u.len();
    let mut u = predicted.u.clone();
    let mut lam = predicted.lambda;
    let tangent = &prev.tangent;

    for iter in 0..max_iter {
        let res = f(&u, lam);
        let arc_val: f64 = u
            .iter()
            .zip(prev.u.iter())
            .enumerate()
            .map(|(i, (&ui, &upi))| (ui - upi) * tangent[i])
            .sum::<f64>()
            + (lam - prev.lambda) * tangent[n]
            - prev.ds;

        let res_norm: f64 = (res.iter().map(|r| r * r).sum::<f64>() + arc_val * arc_val).sqrt();

        if res_norm < tol {
            return CorrectorResult {
                converged: true,
                iterations: iter,
                residual: res_norm,
                u,
                lambda: lam,
            };
        }

        let eps = 1e-7;
        let res_lam_plus = f(&u, lam + eps);
        let f_lam: Vec<f64> = res
            .iter()
            .zip(res_lam_plus.iter())
            .map(|(r, rp)| (rp - r) / eps)
            .collect();

        let j = jac(&u, lam);
        let m = n + 1;
        let mut mat: Vec<Vec<f64>> = Vec::with_capacity(m);
        for i in 0..n {
            let mut row = j[i].clone();
            row.push(f_lam[i]);
            row.push(-res[i]);
            mat.push(row);
        }
        {
            let mut row: Vec<f64> = tangent[..n].to_vec();
            row.push(tangent[n]);
            row.push(-arc_val);
            mat.push(row);
        }

        // Gaussian elimination with partial pivoting
        for col in 0..m {
            let mut max_row = col;
            let mut max_val = mat[col][col].abs();
            for (offset, row_vec) in mat[(col + 1)..m].iter().enumerate() {
                if row_vec[col].abs() > max_val {
                    max_val = row_vec[col].abs();
                    max_row = col + 1 + offset;
                }
            }
            mat.swap(col, max_row);
            let pivot = mat[col][col];
            if pivot.abs() < 1e-14 {
                return CorrectorResult {
                    converged: false,
                    iterations: iter,
                    residual: res_norm,
                    u,
                    lambda: lam,
                };
            }
            for row in (col + 1)..m {
                let factor = mat[row][col] / pivot;
                let col_slice: Vec<f64> = mat[col][col..=m].to_vec();
                for (cell, &cv) in mat[row][col..=m].iter_mut().zip(col_slice.iter()) {
                    *cell -= factor * cv;
                }
            }
        }
        // Back-substitution
        let mut delta = vec![0.0_f64; m];
        for i in (0..m).rev() {
            let mut s = mat[i][m];
            for jj in (i + 1)..m {
                s -= mat[i][jj] * delta[jj];
            }
            delta[i] = s / mat[i][i];
        }

        for (ui, &di) in u.iter_mut().zip(delta[..n].iter()) {
            *ui += di;
        }
        lam += delta[n];
    }

    let res = f(&u, lam);
    let res_norm: f64 = res.iter().map(|r| r * r).sum::<f64>().sqrt();
    CorrectorResult {
        converged: false,
        iterations: max_iter,
        residual: res_norm,
        u,
        lambda: lam,
    }
}

// ---------------------------------------------------------------------------
// matrix_determinant
// ---------------------------------------------------------------------------

/// Compute the determinant of a square matrix (row-major `Vec<Vec`f64`>`).
///
/// Uses LU decomposition with partial pivoting.
pub fn matrix_determinant(mat: &[Vec<f64>]) -> f64 {
    let n = mat.len();
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return mat[0][0];
    }
    if n == 2 {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    let mut sign = 1.0_f64;
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for (offset, row_vec) in a[(col + 1)..n].iter().enumerate() {
            if row_vec[col].abs() > max_val {
                max_val = row_vec[col].abs();
                max_row = col + 1 + offset;
            }
        }
        if max_row != col {
            a.swap(col, max_row);
            sign = -sign;
        }
        let pivot = a[col][col];
        if pivot.abs() < 1e-14 {
            return 0.0;
        }
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            let col_slice: Vec<f64> = a[col][col..n].to_vec();
            for (cell, &cv) in a[row][col..n].iter_mut().zip(col_slice.iter()) {
                *cell -= factor * cv;
            }
        }
    }
    let diag_product: f64 = (0..n).map(|i| a[i][i]).product();
    sign * diag_product
}

// ---------------------------------------------------------------------------
// detect_fold_point
// ---------------------------------------------------------------------------

/// Detect a singular point of the Jacobian by checking for a sign change in
/// its determinant along two consecutive states.
///
/// Returns `true` if the determinants at `state_a` and `state_b` have opposite
/// signs (or one of them is zero), indicating a fold *or a branch point*
/// between them; [`classify_bifurcation`] tells the two apart.
pub fn detect_fold_point(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &JacobianFn,
) -> bool {
    let det_a = matrix_determinant(&jac(&state_a.u, state_a.lambda));
    let det_b = matrix_determinant(&jac(&state_b.u, state_b.lambda));
    det_a * det_b <= 0.0
}

// ---------------------------------------------------------------------------
// stability_index
// ---------------------------------------------------------------------------

/// Relative tolerance (times the largest matrix entry) below which the real
/// part of an eigenvalue counts as zero when unstable eigenvalues are counted.
///
/// Computed eigenvalues on the imaginary axis (a Hopf point) carry real parts
/// of the order of the rounding error; they must not be counted as unstable.
const EIGENVALUE_SIGN_TOLERANCE: f64 = 1e-12;

/// Largest absolute entry of `mat`.
fn max_abs_entry(mat: &[Vec<f64>]) -> f64 {
    mat.iter()
        .flatten()
        .fold(0.0_f64, |acc, value| acc.max(value.abs()))
}

/// Eigenvalues `(re, im)` of `[[a, b], [c, d]]`, sorted by real part.
///
/// The eigenvalue of larger magnitude is computed first and its partner from
/// the determinant, which avoids the cancellation of `(tr ± √disc) / 2` for an
/// eigenvalue that is small next to the trace.
fn eigenvalues_2x2(a: f64, b: f64, c: f64, d: f64) -> [(f64, f64); 2] {
    let half_trace = 0.5 * (a + d);
    let half_difference = 0.5 * (a - d);
    let discriminant = half_difference * half_difference + b * c;
    if discriminant < 0.0 {
        let imaginary = (-discriminant).sqrt();
        return [(half_trace, -imaginary), (half_trace, imaginary)];
    }
    let root = discriminant.sqrt();
    let larger = if half_trace >= 0.0 {
        half_trace + root
    } else {
        half_trace - root
    };
    let smaller = if larger == 0.0 {
        0.0
    } else {
        (a * d - b * c) / larger
    };
    if smaller <= larger {
        [(smaller, 0.0), (larger, 0.0)]
    } else {
        [(larger, 0.0), (smaller, 0.0)]
    }
}

/// Spectrum `(re, im)` of the square Jacobian `mat`: closed form for `n <= 2`,
/// real Schur form beyond.
///
/// Returns `None` when `mat` is not square, holds a non-finite entry, or the
/// eigenvalue iteration does not converge.
fn jacobian_spectrum(mat: &[Vec<f64>]) -> Option<Vec<(f64, f64)>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) || mat.iter().flatten().any(|v| !v.is_finite()) {
        return None;
    }
    match n {
        0 => Some(Vec::new()),
        1 => Some(vec![(mat[0][0], 0.0)]),
        2 => Some(eigenvalues_2x2(mat[0][0], mat[0][1], mat[1][0], mat[1][1]).to_vec()),
        _ => real_matrix_eigenvalues_rows(mat),
    }
}

/// Count the eigenvalues with positive real part (unstable eigenvalues) of a
/// square real matrix.
///
/// The eigenvalues are exact: closed form for 1×1 and 2×2 matrices and the real
/// Schur form for larger ones.  An eigenvalue counts as unstable when its real
/// part exceeds `1e-12` times the largest matrix entry, so the rounding noise of
/// an eigenvalue on the imaginary axis is not counted.  Returns `0` for an empty
/// matrix and for input without a computable spectrum (not square, non-finite
/// entries, or an eigenvalue iteration that fails to converge).
pub fn stability_index(mat: &[Vec<f64>]) -> usize {
    let Some(spectrum) = jacobian_spectrum(mat) else {
        return 0;
    };
    let threshold = EIGENVALUE_SIGN_TOLERANCE * max_abs_entry(mat);
    spectrum.iter().filter(|&&(re, _)| re > threshold).count()
}

// ---------------------------------------------------------------------------
// detect_bifurcation
// ---------------------------------------------------------------------------

/// Detect a bifurcation point by checking whether the number of Jacobian
/// eigenvalues with positive real part ([`stability_index`]) differs between
/// two consecutive states.
///
/// Returns `true` if a bifurcation is detected.
pub fn detect_bifurcation(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &JacobianFn,
) -> bool {
    let ja = jac(&state_a.u, state_a.lambda);
    let jb = jac(&state_b.u, state_b.lambda);
    let idx_a = stability_index(&ja);
    let idx_b = stability_index(&jb);
    idx_a != idx_b
}

// ---------------------------------------------------------------------------
// BifurcationPoint
// ---------------------------------------------------------------------------

/// Classification of a detected bifurcation.
#[derive(Debug, Clone, PartialEq)]
pub enum BifurcationType {
    /// Fold (saddle-node / limit) point: the determinant of the Jacobian changes
    /// sign and the branch turns back in λ (the λ-component of its tangent
    /// changes sign).
    Fold,
    /// Branch point (pitchfork or transcritical): the determinant of the
    /// Jacobian changes sign while the branch keeps running in the same λ
    /// direction, so a second branch crosses it.
    Pitchfork,
    /// Hopf bifurcation: a complex-conjugate pair crosses the imaginary axis,
    /// changing the number of unstable eigenvalues without a determinant sign
    /// change.
    Hopf,
    /// Unclassified bifurcation: the determinant changes sign but the states
    /// carry no usable branch tangent to separate a fold from a branch point.
    Unknown,
}

/// A bifurcation point detected along a branch.
#[derive(Debug, Clone)]
pub struct BifurcationPoint {
    /// Interpolated state at the bifurcation.
    pub state: ContinuationState,
    /// Classification of the bifurcation type.
    pub bif_type: BifurcationType,
    /// Value of det(J) at this point.
    pub det_j: f64,
}

/// Direction of the branch in λ at `state`: the λ-component of its tangent, if
/// the tangent has the expected length `n + 1` and is finite.
fn lambda_tangent(state: &ContinuationState) -> Option<f64> {
    state
        .tangent
        .get(state.u.len())
        .copied()
        .filter(|component| component.is_finite())
}

/// Midpoint of two states with the (normalised) mean of their tangents.
///
/// For a fold the two λ-components have opposite signs and cancel, so the mean
/// tangent at the turning point lies in the `u` directions, as it should.
fn midpoint_state(state_a: &ContinuationState, state_b: &ContinuationState) -> ContinuationState {
    let u_mid: Vec<f64> = state_a
        .u
        .iter()
        .zip(&state_b.u)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let lambda_mid = 0.5 * (state_a.lambda + state_b.lambda);
    let n = u_mid.len();
    let mut mid_state = ContinuationState::new(lambda_mid, u_mid, state_a.ds);
    if state_a.tangent.len() == n + 1 && state_b.tangent.len() == n + 1 {
        let mean: Vec<f64> = state_a
            .tangent
            .iter()
            .zip(&state_b.tangent)
            .map(|(a, b)| a + b)
            .collect();
        let norm = mean.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm.is_finite() && norm > 1e-14 {
            mid_state.tangent = mean.iter().map(|x| x / norm).collect();
        }
    }
    mid_state
}

/// Classify the bifurcation between two consecutive continuation states.
///
/// A bifurcation lies between the states when the determinant of the Jacobian
/// changes sign (a real eigenvalue crosses zero) or, without a sign change,
/// when the number of unstable eigenvalues ([`stability_index`]) changes (a
/// complex-conjugate pair crosses the imaginary axis).  The type follows from
/// the branch tangents stored in the two states:
///
/// * determinant changes sign and the λ-component of the tangent changes sign
///   (or vanishes): the branch turns back in λ, a [`BifurcationType::Fold`];
/// * determinant changes sign while the tangent keeps running in the same λ
///   direction: the branch passes straight through the singular point, a
///   branch point, reported as [`BifurcationType::Pitchfork`];
/// * determinant keeps its sign but the stability index changes:
///   [`BifurcationType::Hopf`];
/// * determinant changes sign but a state has no usable tangent:
///   [`BifurcationType::Unknown`].
///
/// States returned by [`PseudoArcLengthContinuation::step`] carry the tangent of
/// the branch.  A state built with [`ContinuationState::new`] holds the placeholder
/// `[0, …, 0, 1]` (pure λ direction), so a pair of such states is never
/// classified as a fold.
///
/// Returns `None` when neither the determinant sign nor the stability index
/// differs between the states.  The returned state is the midpoint of the pair.
pub fn classify_bifurcation(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &JacobianFn,
) -> Option<BifurcationPoint> {
    let jac_a = jac(&state_a.u, state_a.lambda);
    let jac_b = jac(&state_b.u, state_b.lambda);
    let det_a = matrix_determinant(&jac_a);
    let det_b = matrix_determinant(&jac_b);
    let determinant_changes_sign = det_a * det_b <= 0.0;
    if !determinant_changes_sign && stability_index(&jac_a) == stability_index(&jac_b) {
        return None;
    }

    let mid_state = midpoint_state(state_a, state_b);
    let det_mid = matrix_determinant(&jac(&mid_state.u, mid_state.lambda));
    let bif_type = if !determinant_changes_sign {
        BifurcationType::Hopf
    } else {
        match (lambda_tangent(state_a), lambda_tangent(state_b)) {
            (Some(before), Some(after)) if before * after <= 0.0 => BifurcationType::Fold,
            (Some(_), Some(_)) => BifurcationType::Pitchfork,
            _ => BifurcationType::Unknown,
        }
    };
    Some(BifurcationPoint {
        state: mid_state,
        bif_type,
        det_j: det_mid,
    })
}

// ---------------------------------------------------------------------------
// branch_switching
// ---------------------------------------------------------------------------

/// Compute a unit vector spanning the (approximate) null space of the square
/// Jacobian `j`: the right singular vector of its smallest singular value.
///
/// The sign is fixed so that the component of largest magnitude is positive,
/// which makes the result deterministic.  Returns `None` when `j` is empty, is
/// not an `n × n` matrix, holds non-finite entries, or when the SVD does not
/// converge.
fn jacobian_null_direction(j: &[Vec<f64>], n: usize) -> Option<Vec<f64>> {
    if n == 0 || j.len() != n || j.iter().any(|row| row.len() != n) {
        return None;
    }
    if j.iter().flatten().any(|v| !v.is_finite()) {
        return None;
    }
    let matrix = DMatrix::from_fn(n, n, |r, c| j[r][c]);
    // The unordered SVD avoids nalgebra's NaN-panicking sort; the smallest
    // singular value is located explicitly below.
    let svd = SVD::try_new_unordered(matrix, false, true, f64::EPSILON * 5.0, 100 * n)?;
    let v_t = svd.v_t.as_ref()?;
    let (smallest, _) = svd
        .singular_values
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))?;
    let mut phi: Vec<f64> = (0..n).map(|c| v_t[(smallest, c)]).collect();
    let norm = phi.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !norm.is_finite() || norm < 1e-14 {
        return None;
    }
    let pivot = phi.iter().copied().fold(
        0.0_f64,
        |best, x| if x.abs() > best.abs() { x } else { best },
    );
    let sign = if pivot < 0.0 { -1.0 } else { 1.0 };
    for x in &mut phi {
        *x = sign * *x / norm;
    }
    Some(phi)
}

/// Predict a point on the secondary branch at a simple bifurcation point.
///
/// At a simple bifurcation the Jacobian `J = ∂F/∂u` is singular with a
/// one-dimensional kernel spanned by `φ`.  This function computes `φ` as the
/// right singular vector of the smallest singular value of `J(u, λ)` (so that
/// `J φ ≈ 0`, a genuine null vector rather than a column of `J`), and returns
/// the state displaced by `epsilon` along it: `u + ε φ`, with `λ` unchanged and
/// the tangent set to `[φ; 0]`.  The sign of `φ` is fixed (largest component
/// positive); pass a negative `epsilon` to switch to the opposite side.
///
/// The result is only a predictor; use [`BranchSwitching::switch`] to also
/// correct it onto the secondary branch.
///
/// Returns `None` if `epsilon` is zero or not finite, or if no null direction
/// can be computed (empty, non-square or non-finite Jacobian, SVD failure).
pub fn branch_switching(
    state: &ContinuationState,
    jac: &JacobianFn,
    epsilon: f64,
) -> Option<ContinuationState> {
    if !epsilon.is_finite() || epsilon == 0.0 {
        return None;
    }
    let n = state.u.len();
    let j = jac(&state.u, state.lambda);
    let null_dir = jacobian_null_direction(&j, n)?;
    let u_new: Vec<f64> = state
        .u
        .iter()
        .zip(&null_dir)
        .map(|(&ui, &di)| ui + epsilon * di)
        .collect();
    let mut tangent_new = null_dir;
    tangent_new.push(0.0);
    Some(ContinuationState {
        lambda: state.lambda,
        u: u_new,
        tangent: tangent_new,
        arc_length: state.arc_length,
        ds: state.ds,
    })
}

// ---------------------------------------------------------------------------
// BranchSwitching — bordered-system approach
// ---------------------------------------------------------------------------

/// Bordered-Newton branch switching at a simple bifurcation point.
///
/// The predictor displaces the bifurcation state by `epsilon` along the null
/// vector `φ` of `J = ∂F/∂u` (see [`branch_switching`]).  The corrector is
/// Newton's method on the bordered system
///
/// ```text
/// [ J    F_λ ] [Δu]     [ −F(u, λ)             ]
/// [ φᵀ    0  ] [Δλ]  =  [ −(φ·(u − u₀) − ε)    ]
/// ```
///
/// which pins the displacement along `φ` to `ε` and so converges to the
/// secondary branch rather than back to the primary one.
pub struct BranchSwitching {
    /// Perturbation magnitude ε applied in the null-space direction.
    pub epsilon: f64,
    /// Maximum Newton iterations for the corrector after switching.
    pub max_iter: usize,
    /// Newton tolerance.
    pub tol: f64,
}

impl BranchSwitching {
    /// Create a new `BranchSwitching` with given parameters.
    pub fn new(epsilon: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            epsilon,
            max_iter,
            tol,
        }
    }

    /// Switch branch at the bifurcation `state` using the bordered system.
    ///
    /// `f` is the residual `F(u, λ)` and `jac` its Jacobian with respect to
    /// `u`.  Returns the converged state on the secondary branch, whose tangent
    /// is the (normalised) secant from the bifurcation point and whose
    /// `arc_length` is advanced by the length of that secant.  Returns `None`
    /// if the null direction cannot be computed (see [`branch_switching`]) or
    /// if the bordered Newton corrector fails to converge to `tol` within
    /// `max_iter` iterations.  A negative `epsilon` selects the opposite side
    /// of the secondary branch.
    pub fn switch(
        &self,
        state: &ContinuationState,
        f: &ResidualFn,
        jac: &JacobianFn,
    ) -> Option<ContinuationState> {
        let predicted = branch_switching(state, jac, self.epsilon)?;
        // Anchor of the bordering constraint φ·(u − u₀) = ε: the bifurcation
        // point itself, with border vector [φ; 0] and step ε.
        let anchor = ContinuationState {
            lambda: state.lambda,
            u: state.u.clone(),
            tangent: predicted.tangent.clone(),
            arc_length: state.arc_length,
            ds: self.epsilon,
        };
        let corrected = corrector_newton(&predicted, &anchor, f, jac, self.tol, self.max_iter);
        if !corrected.converged {
            return None;
        }

        let mut tangent: Vec<f64> = corrected
            .u
            .iter()
            .zip(&state.u)
            .map(|(new, old)| new - old)
            .collect();
        tangent.push(corrected.lambda - state.lambda);
        let secant_length = tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
        if secant_length > 1e-14 {
            for x in &mut tangent {
                *x /= secant_length;
            }
        }
        Some(ContinuationState {
            lambda: corrected.lambda,
            u: corrected.u,
            tangent,
            arc_length: state.arc_length + secant_length,
            ds: state.ds,
        })
    }
}

// ---------------------------------------------------------------------------
// TurningPointLocator — Moore-Penrose continuation past fold points
// ---------------------------------------------------------------------------

/// Locator for turning (fold) points using bisection on the determinant.
///
/// Iteratively bisects between `state_a` (det > 0) and `state_b` (det < 0)
/// to find the fold point to within `tol` in arc length.
pub struct TurningPointLocator {
    /// Convergence tolerance on the determinant magnitude.
    pub tol: f64,
    /// Maximum bisection iterations.
    pub max_iter: usize,
}

impl TurningPointLocator {
    /// Create a new `TurningPointLocator`.
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }

    /// Locate the turning point between `state_a` and `state_b`.
    ///
    /// Returns the interpolated state closest to the fold, together with
    /// the determinant value at that state.
    pub fn locate(
        &self,
        state_a: &ContinuationState,
        state_b: &ContinuationState,
        jac: &JacobianFn,
    ) -> (ContinuationState, f64) {
        let mut alpha_lo = 0.0_f64;
        let mut alpha_hi = 1.0_f64;

        let interp = |alpha: f64| {
            let u_i: Vec<f64> = state_a
                .u
                .iter()
                .zip(&state_b.u)
                .map(|(a, b)| a + alpha * (b - a))
                .collect();
            let lam_i = state_a.lambda + alpha * (state_b.lambda - state_a.lambda);
            ContinuationState::new(lam_i, u_i, state_a.ds)
        };

        let det_a = matrix_determinant(&jac(&state_a.u, state_a.lambda));
        let mut det_lo = det_a;
        let mut mid_state = interp(0.5_f64);
        let mut det_mid = matrix_determinant(&jac(&mid_state.u, mid_state.lambda));

        for _ in 0..self.max_iter {
            if det_mid.abs() < self.tol {
                break;
            }
            let alpha_mid = (alpha_lo + alpha_hi) / 2.0;
            mid_state = interp(alpha_mid);
            det_mid = matrix_determinant(&jac(&mid_state.u, mid_state.lambda));
            if det_lo * det_mid <= 0.0 {
                alpha_hi = alpha_mid;
            } else {
                alpha_lo = alpha_mid;
                det_lo = det_mid;
            }
        }
        (mid_state, det_mid)
    }
}

// ---------------------------------------------------------------------------
// StabilityAnalysis — eigenvalue tracking
// ---------------------------------------------------------------------------

/// Assignment of stability to a solution branch point.
#[derive(Debug, Clone, PartialEq)]
pub enum StabilityLabel {
    /// All eigenvalues have negative real part (stable).
    Stable,
    /// At least one eigenvalue has positive real part (unstable).
    Unstable,
    /// Marginal stability (eigenvalue on the imaginary axis).
    Marginal,
}

/// Stability analysis along a continuation branch.
///
/// Evaluates the Jacobian eigenvalues at each `ContinuationState` in a
/// sequence and assigns a stability label.
pub struct StabilityAnalysis {
    /// Tolerance for detecting marginal (near-zero real part) eigenvalues.
    pub marginal_tol: f64,
}

impl StabilityAnalysis {
    /// Create a new stability analyser.
    pub fn new(marginal_tol: f64) -> Self {
        Self { marginal_tol }
    }

    /// Assign a stability label to a single state from the eigenvalues of its
    /// Jacobian (closed form for 1×1 and 2×2, real Schur form beyond):
    ///
    /// * [`StabilityLabel::Unstable`] if any eigenvalue has a real part above
    ///   `marginal_tol`;
    /// * otherwise [`StabilityLabel::Marginal`] if any eigenvalue has a real part
    ///   within `marginal_tol` of the imaginary axis;
    /// * otherwise [`StabilityLabel::Stable`].
    ///
    /// A Jacobian without a well-defined spectrum (not square, non-finite
    /// entries) cannot be certified stable and is labelled
    /// [`StabilityLabel::Unstable`].  A negative or NaN `marginal_tol` acts as
    /// zero.
    pub fn label(&self, state: &ContinuationState, jac: &JacobianFn) -> StabilityLabel {
        let Some(spectrum) = jacobian_spectrum(&jac(&state.u, state.lambda)) else {
            return StabilityLabel::Unstable;
        };
        let tolerance = self.marginal_tol.max(0.0);
        let mut near_axis = false;
        for &(re, _) in &spectrum {
            if re > tolerance {
                return StabilityLabel::Unstable;
            }
            if re >= -tolerance {
                near_axis = true;
            }
        }
        if near_axis {
            StabilityLabel::Marginal
        } else {
            StabilityLabel::Stable
        }
    }

    /// Label a sequence of states along a branch.
    pub fn label_branch(
        &self,
        states: &[ContinuationState],
        jac: &JacobianFn,
    ) -> Vec<StabilityLabel> {
        states.iter().map(|s| self.label(s, jac)).collect()
    }
}

// ---------------------------------------------------------------------------
// branch_tangent
// ---------------------------------------------------------------------------

/// Relative size of the central-difference step in λ used for `F_λ`.
const LAMBDA_DIFFERENCE_STEP: f64 = 1e-6;

/// Largest value of `‖τ̃‖‖t‖ = 1 / cos(angle)` between consecutive branch
/// tangents for which the bordered tangent is trusted (a turn of about 89.94°);
/// a sharper turn means the step was too long to resolve the branch direction.
const MAX_TANGENT_TURN: f64 = 1e3;

/// Unit tangent of the solution branch of `F(u, λ) = 0` at `(u, λ)`, oriented
/// along `reference`.
///
/// Solves the bordered system
///
/// ```text
/// [ J(u, λ)   F_λ(u, λ) ] [τ̃]   [0]
/// [      tᵀ             ]      = [1]
/// ```
///
/// (`t` is `reference`, `F_λ` a central difference) and normalises `τ̃`.  The
/// system stays regular at a fold, where `J` alone is singular, so the tangent
/// is exact there, unlike a secant between two states.  Returns `None` if the
/// dimensions are inconsistent, the system is (numerically) singular, as at a
/// branch point, or the branch turned too sharply since `reference`.
fn branch_tangent(
    u: &[f64],
    lambda: f64,
    reference: &[f64],
    f: &ResidualFn,
    jac: &JacobianFn,
) -> Option<Vec<f64>> {
    let n = u.len();
    let m = n + 1;
    if reference.len() != m {
        return None;
    }
    let j = jac(u, lambda);
    if j.len() != n || j.iter().any(|row| row.len() != n) {
        return None;
    }
    let h = LAMBDA_DIFFERENCE_STEP * lambda.abs().max(1.0);
    let f_plus = f(u, lambda + h);
    let f_minus = f(u, lambda - h);
    if f_plus.len() != n || f_minus.len() != n {
        return None;
    }
    let matrix = DMatrix::from_fn(m, m, |row, col| match (row < n, col < n) {
        (true, true) => j[row][col],
        (true, false) => (f_plus[row] - f_minus[row]) / (2.0 * h),
        (false, _) => reference[col],
    });
    if matrix.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let mut rhs = DVector::<f64>::zeros(m);
    rhs[n] = 1.0;
    let solution = matrix.clone().full_piv_lu().solve(&rhs)?;
    let norm = solution.norm();
    let reference_norm = reference.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !norm.is_finite() || norm * reference_norm > MAX_TANGENT_TURN {
        return None;
    }
    // A numerically singular matrix yields a solution of huge components whose
    // residual is not small relative to its size: reject it.
    let residual = (&matrix * &solution - &rhs).norm();
    let solves_the_system = residual <= 1e-8 * (1.0 + norm * matrix.norm());
    if !solves_the_system {
        return None;
    }
    Some(solution.iter().map(|x| x / norm).collect())
}

/// Normalised secant `(u − u_prev, λ − λ_prev)` between two states (left
/// unnormalised if the two states coincide).
fn secant_tangent(prev: &ContinuationState, u: &[f64], lambda: f64) -> Vec<f64> {
    let n = prev.u.len();
    let mut tangent = vec![0.0_f64; n + 1];
    for (ti, (ri, si)) in tangent[..n].iter_mut().zip(u.iter().zip(prev.u.iter())) {
        *ti = ri - si;
    }
    tangent[n] = lambda - prev.lambda;
    let norm: f64 = tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-14 {
        for x in &mut tangent {
            *x /= norm;
        }
    }
    tangent
}

// ---------------------------------------------------------------------------
// PseudoArcLengthContinuation — combined predictor-corrector stepper
// ---------------------------------------------------------------------------

/// Configuration and state for pseudo-arclength continuation.
///
/// Combines the predictor (`pseudo_arclength_step`) and corrector
/// (`corrector_newton`) into a single stepper.  Step-size adaptation
/// doubles/halves `ds` based on whether the corrector converges within
/// `max_iter_fast` Newton steps.
pub struct PseudoArcLengthContinuation {
    /// Newton convergence tolerance.
    pub tol: f64,
    /// Maximum Newton iterations per corrector call.
    pub max_iter: usize,
    /// If the corrector converges in ≤ this many steps, double `ds`.
    pub max_iter_fast: usize,
    /// Minimum allowed `ds`.
    pub ds_min: f64,
    /// Maximum allowed `ds`.
    pub ds_max: f64,
}

impl PseudoArcLengthContinuation {
    /// Create a new continuation stepper.
    pub fn new(tol: f64, max_iter: usize, max_iter_fast: usize, ds_min: f64, ds_max: f64) -> Self {
        Self {
            tol,
            max_iter,
            max_iter_fast,
            ds_min,
            ds_max,
        }
    }

    /// Advance `state` by one step, returning the accepted corrected state.
    ///
    /// The accepted state carries the unit tangent of the branch at the new
    /// point, oriented along the tangent of `state`.
    ///
    /// Returns `None` if the corrector fails to converge even after halving `ds`.
    pub fn step(
        &self,
        state: &mut ContinuationState,
        f: &ResidualFn,
        jac: &JacobianFn,
    ) -> Option<ContinuationState> {
        // Try current ds; if corrector fails, halve ds up to 5 times.
        for attempt in 0..5usize {
            let _ = attempt;
            let predicted = pseudo_arclength_step(state);
            let result = corrector_newton(&predicted, state, f, jac, self.tol, self.max_iter);
            if result.converged {
                // Tangent of the branch at the accepted point, oriented along the
                // previous tangent.  Fall back to the secant from the previous
                // state where the bordered tangent system is singular.
                let new_tangent = branch_tangent(&result.u, result.lambda, &state.tangent, f, jac)
                    .unwrap_or_else(|| secant_tangent(state, &result.u, result.lambda));

                let mut accepted = ContinuationState {
                    lambda: result.lambda,
                    u: result.u,
                    tangent: new_tangent,
                    arc_length: state.arc_length + state.ds,
                    ds: state.ds,
                };

                // Step-size adaptation
                if result.iterations <= self.max_iter_fast {
                    accepted.ds = (accepted.ds * 2.0).min(self.ds_max);
                }
                return Some(accepted);
            }
            // Halve ds and retry
            state.ds = (state.ds / 2.0).max(self.ds_min);
            if state.ds <= self.ds_min {
                break;
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// PathFollowing — generic path following with event detection
// ---------------------------------------------------------------------------

/// Result of a single path-following step.
#[derive(Debug, Clone)]
pub struct PathStep {
    /// The accepted continuation state.
    pub state: ContinuationState,
    /// Whether the Jacobian determinant changed sign between this and the
    /// previous step (a fold or a branch point; see [`classify_bifurcation`]).
    pub fold_detected: bool,
    /// Whether the number of unstable eigenvalues changed between this and the
    /// previous step.
    pub bifurcation_detected: bool,
}

/// Generic path-following algorithm with event detection.
///
/// Runs pseudo-arclength continuation for up to `max_steps` steps,
/// detecting fold and bifurcation events along the way.
pub struct PathFollowing {
    /// Maximum number of continuation steps.
    pub max_steps: usize,
    /// Continuation stepper configuration.
    pub continuation: PseudoArcLengthContinuation,
}

impl PathFollowing {
    /// Create a new path follower.
    pub fn new(max_steps: usize, continuation: PseudoArcLengthContinuation) -> Self {
        Self {
            max_steps,
            continuation,
        }
    }

    /// Follow the branch starting from `initial_state`.
    ///
    /// Returns the sequence of accepted path steps.
    pub fn follow(
        &self,
        initial_state: ContinuationState,
        f: &ResidualFn,
        jac: &JacobianFn,
    ) -> Vec<PathStep> {
        let mut steps: Vec<PathStep> = Vec::with_capacity(self.max_steps);
        let mut current = initial_state;

        for _ in 0..self.max_steps {
            let prev = current.clone();
            match self.continuation.step(&mut current, f, jac) {
                None => break,
                Some(accepted) => {
                    let fold = detect_fold_point(&prev, &accepted, jac);
                    let bif = detect_bifurcation(&prev, &accepted, jac);
                    steps.push(PathStep {
                        state: accepted.clone(),
                        fold_detected: fold,
                        bifurcation_detected: bif,
                    });
                    current = accepted;
                }
            }
        }
        steps
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Simple 1-D test problem: F(u, λ) = u² - λ = 0  →  u = ±√λ
    fn f1d(u: &[f64], lam: f64) -> Vec<f64> {
        vec![u[0] * u[0] - lam]
    }
    fn jac1d(u: &[f64], _lam: f64) -> Vec<Vec<f64>> {
        vec![vec![2.0 * u[0]]]
    }

    // 2-D pitchfork: F(u, λ) = [u0^3 - λ*u0, u1 + u0]
    fn f2d_pitch(u: &[f64], lam: f64) -> Vec<f64> {
        vec![u[0].powi(3) - lam * u[0], u[1] + u[0]]
    }
    fn jac2d_pitch(u: &[f64], lam: f64) -> Vec<Vec<f64>> {
        vec![vec![3.0 * u[0] * u[0] - lam, 0.0], vec![1.0, 1.0]]
    }

    #[test]
    fn test_continuation_state_new() {
        let s = ContinuationState::new(0.0, vec![1.0, 2.0], 0.1);
        assert_eq!(s.dim(), 2);
        assert_eq!(s.lambda, 0.0);
        assert_eq!(s.ds, 0.1);
        assert!((s.tangent[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pseudo_arclength_step_increments_arc_length() {
        let s = ContinuationState::new(1.0, vec![1.0], 0.1);
        let s2 = pseudo_arclength_step(&s);
        assert!((s2.arc_length - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_pseudo_arclength_step_lambda_moves() {
        let s = ContinuationState::new(1.0, vec![1.0], 0.5);
        let s2 = pseudo_arclength_step(&s);
        assert!((s2.lambda - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_pseudo_arclength_step_u_unchanged_when_tangent_zero() {
        let mut s = ContinuationState::new(0.0, vec![3.0, 4.0], 0.2);
        s.tangent = vec![0.0, 0.0, 1.0];
        let s2 = pseudo_arclength_step(&s);
        assert!((s2.u[0] - 3.0).abs() < 1e-12);
        assert!((s2.u[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_corrector_newton_converges_1d() {
        let prev = ContinuationState::new(1.0, vec![1.0], 0.1);
        let mut predicted = ContinuationState::new(1.0, vec![1.1], 0.1);
        predicted.tangent = prev.tangent.clone();
        let result = corrector_newton(&predicted, &prev, &f1d, &jac1d, 1e-10, 50);
        assert!(result.converged, "Newton should converge");
        assert!(result.residual < 1e-8);
    }

    #[test]
    fn test_corrector_newton_result_satisfies_equation() {
        let prev = ContinuationState::new(4.0, vec![2.0], 0.1);
        let mut predicted = ContinuationState::new(4.0, vec![2.05], 0.1);
        predicted.tangent = prev.tangent.clone();
        let result = corrector_newton(&predicted, &prev, &f1d, &jac1d, 1e-10, 50);
        if result.converged {
            let res = f1d(&result.u, result.lambda);
            assert!(res[0].abs() < 1e-6);
        }
    }

    #[test]
    fn test_matrix_determinant_2x2() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let det = matrix_determinant(&m);
        assert!((det - (1.0 * 4.0 - 2.0 * 3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_determinant_identity_3x3() {
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        assert!((matrix_determinant(&m) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_determinant_singular() {
        let m = vec![
            vec![1.0, 2.0, 3.0],
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
        ];
        assert!(matrix_determinant(&m).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_determinant_1x1() {
        let m = vec![vec![7.0]];
        assert!((matrix_determinant(&m) - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_determinant_empty() {
        let m: Vec<Vec<f64>> = vec![];
        assert!((matrix_determinant(&m) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_detect_fold_point_true() {
        let sa = ContinuationState::new(0.5, vec![0.5], 0.1);
        let sb = ContinuationState::new(-0.5, vec![-0.5], 0.1);
        let fold = detect_fold_point(&sa, &sb, &jac1d);
        assert!(
            fold,
            "fold should be detected across a sign change in det(J)"
        );
    }

    #[test]
    fn test_detect_fold_point_false() {
        let sa = ContinuationState::new(1.0, vec![1.0], 0.1);
        let sb = ContinuationState::new(2.0, vec![2.0], 0.1);
        let fold = detect_fold_point(&sa, &sb, &jac1d);
        assert!(!fold);
    }

    #[test]
    fn test_stability_index_1x1_positive() {
        let m = vec![vec![3.0]];
        assert_eq!(stability_index(&m), 1);
    }

    #[test]
    fn test_stability_index_1x1_negative() {
        let m = vec![vec![-2.0]];
        assert_eq!(stability_index(&m), 0);
    }

    #[test]
    fn test_stability_index_2x2_stable() {
        let m = vec![vec![-2.0, 0.0], vec![0.0, -3.0]];
        assert_eq!(stability_index(&m), 0);
    }

    #[test]
    fn test_stability_index_2x2_unstable() {
        let m = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        assert_eq!(stability_index(&m), 2);
    }

    #[test]
    fn test_stability_index_2x2_one_unstable() {
        let m = vec![vec![1.0, 0.0], vec![0.0, -2.0]];
        assert_eq!(stability_index(&m), 1);
    }

    #[test]
    fn test_stability_index_empty() {
        let m: Vec<Vec<f64>> = vec![];
        assert_eq!(stability_index(&m), 0);
    }

    #[test]
    fn test_detect_bifurcation_detects_change() {
        let sa = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let sb = ContinuationState::new(1.0, vec![1.0, 1.0], 0.1);
        let jac_stable = |_u: &[f64], _lam: f64| vec![vec![-1.0, 0.0], vec![0.0, -1.0]];
        assert!(!detect_bifurcation(&sa, &sb, &jac_stable));
        let jac_crossing = |u: &[f64], _lam: f64| {
            let v = u[0];
            vec![vec![v, 0.0], vec![0.0, v]]
        };
        let sa2 = ContinuationState::new(0.0, vec![-1.0, 0.0], 0.1);
        let sb2 = ContinuationState::new(1.0, vec![1.0, 0.0], 0.1);
        assert!(detect_bifurcation(&sa2, &sb2, &jac_crossing));
    }

    #[test]
    fn test_branch_switching_perturbs_solution() {
        let s = ContinuationState::new(1.0, vec![1.0, 0.0], 0.1);
        let jac2d = |_u: &[f64], _lam: f64| vec![vec![1e-15, 0.0], vec![0.0, 1.0]];
        let s_branch = branch_switching(&s, &jac2d, 0.01).expect("null direction exists");
        // The (near-)null vector of diag(1e-15, 1) is e₀.
        assert!(
            (s_branch.u[0] - 1.01).abs() < 1e-12,
            "u0 = {}",
            s_branch.u[0]
        );
        assert!(s_branch.u[1].abs() < 1e-12, "u1 = {}", s_branch.u[1]);
    }

    #[test]
    fn test_branch_switching_preserves_lambda() {
        let s = ContinuationState::new(2.5, vec![1.0], 0.05);
        let jac_id = |_u: &[f64], _lam: f64| vec![vec![1.0]];
        let s2 = branch_switching(&s, &jac_id, 0.1).expect("1x1 Jacobian has a null direction");
        assert!((s2.lambda - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_corrector_newton_diverges_on_singular() {
        let f_zero = |_u: &[f64], _lam: f64| vec![0.0];
        let jac_zero = |_u: &[f64], _lam: f64| vec![vec![0.0]];
        let prev = ContinuationState::new(0.0, vec![0.0], 0.1);
        let mut predicted = ContinuationState::new(0.0, vec![0.1], 0.1);
        predicted.tangent = prev.tangent.clone();
        let _result = corrector_newton(&predicted, &prev, &f_zero, &jac_zero, 1e-10, 5);
    }

    #[test]
    fn test_arc_length_accumulates_over_steps() {
        let s0 = ContinuationState::new(0.0, vec![0.0], 0.2);
        let s1 = pseudo_arclength_step(&s0);
        let s2 = pseudo_arclength_step(&s1);
        assert!((s2.arc_length - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_matrix_determinant_4x4_diagonal() {
        let m = vec![
            vec![2.0, 0.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 5.0, 0.0],
            vec![0.0, 0.0, 0.0, 7.0],
        ];
        let det = matrix_determinant(&m);
        assert!((det - 210.0).abs() < 1e-8);
    }

    #[test]
    fn test_continuation_state_dim() {
        let s = ContinuationState::new(0.0, vec![1.0, 2.0, 3.0], 0.1);
        assert_eq!(s.dim(), 3);
    }

    #[test]
    fn test_tangent_length_matches_n_plus_1() {
        let s = ContinuationState::new(0.0, vec![1.0, 2.0, 3.0], 0.1);
        assert_eq!(s.tangent.len(), 4);
    }

    #[test]
    fn test_classify_bifurcation_none_when_same_stability() {
        let sa = ContinuationState::new(1.0, vec![1.0], 0.1);
        let sb = ContinuationState::new(2.0, vec![2.0], 0.1);
        // Both: det(J) = 2u > 0, stability index unchanged
        let bif = classify_bifurcation(&sa, &sb, &jac1d);
        assert!(bif.is_none());
    }

    #[test]
    fn test_stability_analysis_stable_label() {
        let sa = StabilityAnalysis::new(1e-6);
        let state = ContinuationState::new(0.0, vec![0.0], 0.1);
        let jac_neg = |_u: &[f64], _lam: f64| vec![vec![-1.0]];
        assert_eq!(sa.label(&state, &jac_neg), StabilityLabel::Stable);
    }

    #[test]
    fn test_stability_analysis_unstable_label() {
        let sa = StabilityAnalysis::new(1e-6);
        let state = ContinuationState::new(0.0, vec![0.0], 0.1);
        let jac_pos = |_u: &[f64], _lam: f64| vec![vec![1.0]];
        assert_eq!(sa.label(&state, &jac_pos), StabilityLabel::Unstable);
    }

    #[test]
    fn test_stability_analysis_marginal_label() {
        let sa = StabilityAnalysis::new(1e-6);
        let state = ContinuationState::new(0.0, vec![0.0], 0.1);
        let jac_zero = |_u: &[f64], _lam: f64| vec![vec![0.0]];
        assert_eq!(sa.label(&state, &jac_zero), StabilityLabel::Marginal);
    }

    #[test]
    fn test_stability_analysis_branch_labels() {
        let sa = StabilityAnalysis::new(1e-6);
        let states = vec![
            ContinuationState::new(0.0, vec![-1.0], 0.1),
            ContinuationState::new(0.0, vec![0.0], 0.1),
            ContinuationState::new(0.0, vec![1.0], 0.1),
        ];
        let jac_sign = |u: &[f64], _lam: f64| vec![vec![u[0]]];
        let labels = sa.label_branch(&states, &jac_sign);
        assert_eq!(labels.len(), 3);
        assert_eq!(labels[0], StabilityLabel::Stable);
        assert_eq!(labels[2], StabilityLabel::Unstable);
    }

    #[test]
    fn test_turning_point_locator_basic() {
        // det(J) = 2u; sign changes between u = 0.5 and u = -0.5
        let sa = ContinuationState::new(0.25, vec![0.5], 0.1);
        let sb = ContinuationState::new(0.25, vec![-0.5], 0.1);
        let locator = TurningPointLocator::new(1e-6, 50);
        let (mid, det) = locator.locate(&sa, &sb, &jac1d);
        // midpoint u ≈ 0 → det ≈ 0
        assert!(
            det.abs() < 0.1,
            "det at turning point should be near 0, got {}",
            det
        );
        let _ = mid;
    }

    #[test]
    fn test_normalised_tangent_length() {
        let s = ContinuationState::new(0.0, vec![3.0, 4.0], 0.1);
        let nt = s.normalised_tangent();
        let norm: f64 = nt.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(
            (norm - 1.0).abs() < 1e-10,
            "normalised tangent should have unit norm"
        );
    }

    #[test]
    fn test_branch_switching_struct() {
        // 2-D pitchfork F = [u0³ − λu0, u1 + u0] at the trivial solution
        // (u = 0, λ = 0): J = [[0, 0], [1, 1]] has kernel span{[1, −1]/√2}.
        let eps = 0.01;
        let bs = BranchSwitching::new(eps, 20, 1e-10);
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let s2 = bs
            .switch(&s, &f2d_pitch, &jac2d_pitch)
            .expect("bordered Newton must converge");
        let u0 = eps / 2.0_f64.sqrt();
        assert!((s2.u[0] - u0).abs() < 1e-10, "u0 = {}", s2.u[0]);
        assert!((s2.u[1] + u0).abs() < 1e-10, "u1 = {}", s2.u[1]);
        // On the secondary branch u0² = λ.
        assert!((s2.lambda - u0 * u0).abs() < 1e-10, "λ = {}", s2.lambda);
        let res = f2d_pitch(&s2.u, s2.lambda);
        assert!(res.iter().all(|r| r.abs() < 1e-9), "residual {res:?}");
    }

    /// Pitchfork in `v` with an inert component `u`: F = (u, v³ − λv).
    fn f_pitch_uv(u: &[f64], lam: f64) -> Vec<f64> {
        vec![u[0], u[1].powi(3) - lam * u[1]]
    }
    fn jac_pitch_uv(u: &[f64], lam: f64) -> Vec<Vec<f64>> {
        vec![vec![1.0, 0.0], vec![0.0, 3.0 * u[1] * u[1] - lam]]
    }

    fn l2(x: &[f64]) -> f64 {
        x.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    #[test]
    fn test_branch_switching_null_direction_is_kernel_not_jacobian_column() {
        // Rank-1 Jacobian: range = span{[1, 2]}, kernel = span{[1, −2]}.  The
        // old implementation copied the smallest-norm *column* [1, 2]/√5 (a
        // range vector) instead of the kernel vector.
        let jac = |_u: &[f64], _lam: f64| vec![vec![2.0, 1.0], vec![4.0, 2.0]];
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let eps = 0.1;
        let b = branch_switching(&s, &jac, eps).expect("null direction exists");
        let phi: Vec<f64> = b.u.iter().map(|x| x / eps).collect();
        assert!((l2(&phi) - 1.0).abs() < 1e-12, "|φ| = {}", l2(&phi));
        let j = jac(&s.u, s.lambda);
        let j_phi: Vec<f64> = j
            .iter()
            .map(|row| row.iter().zip(&phi).map(|(a, b)| a * b).sum())
            .collect();
        assert!(l2(&j_phi) < 1e-12, "J·φ = {j_phi:?}");
        // Sign convention: the component of largest magnitude is positive.
        let expected = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt()];
        assert!((phi[0] - expected[0]).abs() < 1e-12, "φ = {phi:?}");
        assert!((phi[1] - expected[1]).abs() < 1e-12, "φ = {phi:?}");
        // The tangent is [φ; 0].
        assert_eq!(b.tangent.len(), 3);
        assert!((b.tangent[0] - phi[0]).abs() < 1e-12);
        assert!(b.tangent[2].abs() < 1e-15);
    }

    #[test]
    fn test_branch_switching_zero_column_is_not_misread_as_null_direction() {
        // J = diag(1, 0): ker J = span{e₁}.  A zero column used to trigger the
        // `null_dir[0] = 1` fallback, i.e. the range direction e₀.
        let jac = |_u: &[f64], _lam: f64| vec![vec![1.0, 0.0], vec![0.0, 0.0]];
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let b = branch_switching(&s, &jac, 0.25).expect("null direction exists");
        assert!(b.u[0].abs() < 1e-14, "u0 = {}", b.u[0]);
        assert!((b.u[1] - 0.25).abs() < 1e-14, "u1 = {}", b.u[1]);
    }

    #[test]
    fn test_branch_switching_three_dimensional_kernel() {
        // Graph-Laplacian-like Jacobian with kernel span{[1, 1, 1]/√3}.
        let jac = |_u: &[f64], _lam: f64| {
            vec![
                vec![1.0, -1.0, 0.0],
                vec![0.0, 1.0, -1.0],
                vec![-1.0, 0.0, 1.0],
            ]
        };
        let s = ContinuationState::new(0.0, vec![0.0; 3], 0.1);
        let b = branch_switching(&s, &jac, 0.3).expect("null direction exists");
        let third = 1.0 / 3.0_f64.sqrt();
        for (i, ui) in b.u.iter().enumerate() {
            assert!((ui - 0.3 * third).abs() < 1e-12, "u[{i}] = {ui}");
        }
    }

    #[test]
    fn test_branch_switching_sign_convention_and_negative_epsilon() {
        let jac = |_u: &[f64], _lam: f64| vec![vec![1.0, 0.0], vec![0.0, 0.0]];
        let s = ContinuationState::new(1.5, vec![0.0, 0.0], 0.1);
        let plus = branch_switching(&s, &jac, 0.1).expect("ok");
        let minus = branch_switching(&s, &jac, -0.1).expect("ok");
        // Largest component of φ is positive: φ = +e₁.
        assert!((plus.u[1] - 0.1).abs() < 1e-14);
        assert!((minus.u[1] + 0.1).abs() < 1e-14);
        assert!((plus.tangent[1] - 1.0).abs() < 1e-14);
        assert!((minus.tangent[1] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_branch_switching_rejects_invalid_input() {
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let good = |_u: &[f64], _lam: f64| vec![vec![1.0, 0.0], vec![0.0, 0.0]];
        assert!(branch_switching(&s, &good, 0.0).is_none(), "zero ε");
        assert!(branch_switching(&s, &good, f64::NAN).is_none(), "NaN ε");
        assert!(branch_switching(&s, &good, f64::INFINITY).is_none(), "∞ ε");

        let ragged = |_u: &[f64], _lam: f64| vec![vec![1.0, 0.0], vec![0.0]];
        assert!(branch_switching(&s, &ragged, 0.1).is_none(), "ragged J");
        let short = |_u: &[f64], _lam: f64| vec![vec![1.0, 0.0]];
        assert!(branch_switching(&s, &short, 0.1).is_none(), "1×2 J");
        let nan = |_u: &[f64], _lam: f64| vec![vec![1.0, f64::NAN], vec![0.0, 0.0]];
        assert!(branch_switching(&s, &nan, 0.1).is_none(), "NaN in J");

        let empty = ContinuationState::new(0.0, Vec::new(), 0.1);
        let empty_jac = |_u: &[f64], _lam: f64| Vec::<Vec<f64>>::new();
        assert!(branch_switching(&empty, &empty_jac, 0.1).is_none(), "n = 0");
    }

    #[test]
    fn test_branch_switching_pitchfork_reaches_secondary_branch() {
        // F(u, v, λ) = (u, v³ − λv) at the origin: J = diag(1, 0), φ = ±e_v.
        let eps = 0.05;
        let bs = BranchSwitching::new(eps, 25, 1e-12);
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let sw = bs
            .switch(&s, &f_pitch_uv, &jac_pitch_uv)
            .expect("bordered Newton must converge");
        assert!(sw.u[0].abs() < 1e-10, "u = {}", sw.u[0]);
        assert!((sw.u[1] - eps).abs() < 1e-10, "v = {}", sw.u[1]);
        // v ≠ 0 is the non-trivial branch: λ = v² = ε².
        assert!((sw.lambda - eps * eps).abs() < 1e-10, "λ = {}", sw.lambda);
        let res = f_pitch_uv(&sw.u, sw.lambda);
        assert!(l2(&res) < 1e-10, "|F| = {}", l2(&res));
        // The tangent is the unit secant from the bifurcation point.
        assert!((l2(&sw.tangent) - 1.0).abs() < 1e-12);
        assert!(sw.tangent[1] > 0.99);
        assert!(sw.arc_length > 0.0);
    }

    #[test]
    fn test_branch_switching_negative_epsilon_selects_opposite_branch() {
        let eps = 0.05;
        let bs = BranchSwitching::new(-eps, 25, 1e-12);
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let sw = bs
            .switch(&s, &f_pitch_uv, &jac_pitch_uv)
            .expect("bordered Newton must converge");
        assert!((sw.u[1] + eps).abs() < 1e-10, "v = {}", sw.u[1]);
        assert!((sw.lambda - eps * eps).abs() < 1e-10, "λ = {}", sw.lambda);
        assert!(l2(&f_pitch_uv(&sw.u, sw.lambda)) < 1e-10);
    }

    #[test]
    fn test_branch_switching_transcritical_branch() {
        // F = (u, v² − λv): trivial branch v = 0 crosses v = λ at the origin.
        let f = |u: &[f64], lam: f64| vec![u[0], u[1] * u[1] - lam * u[1]];
        let jac = |u: &[f64], lam: f64| vec![vec![1.0, 0.0], vec![0.0, 2.0 * u[1] - lam]];
        let eps = 0.1;
        let bs = BranchSwitching::new(eps, 25, 1e-12);
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let sw = bs.switch(&s, &f, &jac).expect("converges");
        assert!((sw.u[1] - eps).abs() < 1e-10, "v = {}", sw.u[1]);
        assert!((sw.lambda - eps).abs() < 1e-10, "λ = {}", sw.lambda);
    }

    #[test]
    fn test_branch_switching_honours_corrector_budget() {
        // With no Newton iterations the corrector cannot converge ⇒ None.
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let starved = BranchSwitching::new(0.05, 0, 1e-12);
        assert!(starved.switch(&s, &f_pitch_uv, &jac_pitch_uv).is_none());
        // An unreachable tolerance also fails instead of reporting success.
        let impossible = BranchSwitching::new(0.05, 3, 0.0);
        assert!(impossible.switch(&s, &f_pitch_uv, &jac_pitch_uv).is_none());
        // A zero perturbation would just return the bifurcation point.
        let zero = BranchSwitching::new(0.0, 25, 1e-12);
        assert!(zero.switch(&s, &f_pitch_uv, &jac_pitch_uv).is_none());
    }

    #[test]
    fn test_branch_switching_result_can_be_continued_along_secondary_branch() {
        let bs = BranchSwitching::new(0.05, 25, 1e-12);
        let s = ContinuationState::new(0.0, vec![0.0, 0.0], 0.1);
        let mut state = bs
            .switch(&s, &f_pitch_uv, &jac_pitch_uv)
            .expect("bordered Newton must converge");
        state.ds = 0.05;
        let stepper = PseudoArcLengthContinuation::new(1e-10, 25, 3, 1e-6, 0.2);
        let mut last_v = state.u[1];
        for _ in 0..5 {
            let next = stepper
                .step(&mut state, &f_pitch_uv, &jac_pitch_uv)
                .expect("continuation along the secondary branch must succeed");
            assert!(next.u[1] > last_v, "v must keep growing along the branch");
            // Stay on the non-trivial branch λ = v².
            assert!(
                (next.lambda - next.u[1] * next.u[1]).abs() < 1e-8,
                "λ = {}, v² = {}",
                next.lambda,
                next.u[1] * next.u[1]
            );
            last_v = next.u[1];
            state = next;
        }
    }

    #[test]
    fn test_pseudo_arc_length_continuation_step() {
        // Track u² - λ = 0 from (u=1, λ=1)
        let mut state = ContinuationState::new(1.0, vec![1.0], 0.05);
        // Set initial tangent to move along the positive u-branch
        state.tangent = vec![1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt()];
        let cont = PseudoArcLengthContinuation::new(1e-8, 30, 5, 1e-4, 0.5);
        let result = cont.step(&mut state, &f1d, &jac1d);
        assert!(result.is_some(), "continuation step should succeed");
        let accepted = result.unwrap();
        // The accepted state should still satisfy u² ≈ λ
        let residual = (accepted.u[0] * accepted.u[0] - accepted.lambda).abs();
        assert!(residual < 1e-6, "residual = {}", residual);
    }

    #[test]
    fn test_path_following_runs_multiple_steps() {
        let initial = ContinuationState::new(1.0, vec![1.0], 0.05);
        let cont = PseudoArcLengthContinuation::new(1e-8, 30, 5, 1e-4, 0.5);
        let pf = PathFollowing::new(10, cont);
        let steps = pf.follow(initial, &f1d, &jac1d);
        // Should produce at least some steps
        assert!(!steps.is_empty(), "path following should produce steps");
    }

    #[test]
    fn test_path_following_arc_length_increasing() {
        let initial = ContinuationState::new(1.0, vec![1.0], 0.05);
        let cont = PseudoArcLengthContinuation::new(1e-8, 30, 5, 1e-4, 0.5);
        let pf = PathFollowing::new(5, cont);
        let steps = pf.follow(initial, &f1d, &jac1d);
        for w in steps.windows(2) {
            assert!(
                w[1].state.arc_length >= w[0].state.arc_length,
                "arc length should be non-decreasing"
            );
        }
    }

    #[test]
    fn test_matrix_determinant_3x3_known() {
        // det([[1,2,3],[0,4,5],[1,0,6]]) = 22
        let m = vec![
            vec![1.0, 2.0, 3.0],
            vec![0.0, 4.0, 5.0],
            vec![1.0, 0.0, 6.0],
        ];
        let det = matrix_determinant(&m);
        assert!((det - 22.0).abs() < 1e-8, "det = {}", det);
    }

    #[test]
    fn test_jac2d_pitch_at_zero() {
        // At u=[0,0], lam=0: J = [[0-0, 0], [1, 1]] = [[0,0],[1,1]]
        let j = jac2d_pitch(&[0.0, 0.0], 0.0);
        assert!((j[0][0] - 0.0).abs() < 1e-12);
        assert!((j[1][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_f2d_pitch_at_trivial() {
        // At u=[0,0]: F = [0, 0]
        let res = f2d_pitch(&[0.0, 0.0], 1.0);
        assert!(res[0].abs() < 1e-12);
        assert!(res[1].abs() < 1e-12);
    }
}

#[cfg(test)]
#[path = "continuation_bifurcation_tests.rs"]
mod bifurcation_tests;
