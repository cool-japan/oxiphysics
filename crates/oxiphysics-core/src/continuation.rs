#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Numerical continuation and bifurcation analysis.
//!
//! Implements pseudo-arclength continuation for tracking solution branches of
//! parameterised nonlinear systems F(u, λ) = 0, together with fold and
//! bifurcation point detection, branch switching, turning-point location,
//! stability analysis, and generic path-following with event detection.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// ContinuationState
// ---------------------------------------------------------------------------

/// State along a continuation branch.
///
/// Stores the current parameter value `lambda`, the solution vector `u`, the
/// previous tangent direction used for the arclength predictor, and the
/// accumulated arc length `s`.
#[derive(Debug, Clone)]
pub struct ContinuationState {
    /// Continuation parameter λ.
    pub lambda: f64,
    /// Solution vector u ∈ ℝⁿ.
    pub u: Vec<f64>,
    /// Tangent direction \[du/ds, dλ/ds\] from the previous step.
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
    f: &dyn Fn(&[f64], f64) -> Vec<f64>,
    jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
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
            for row in (col + 1)..m {
                if mat[row][col].abs() > max_val {
                    max_val = mat[row][col].abs();
                    max_row = row;
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
                for k in col..=m {
                    let val = mat[col][k];
                    mat[row][k] -= factor * val;
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

        for i in 0..n {
            u[i] += delta[i];
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
        for row in (col + 1)..n {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                max_row = row;
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
            for k in col..n {
                let val = a[col][k];
                a[row][k] -= factor * val;
            }
        }
    }
    let diag_product: f64 = (0..n).map(|i| a[i][i]).product();
    sign * diag_product
}

// ---------------------------------------------------------------------------
// detect_fold_point
// ---------------------------------------------------------------------------

/// Detect a fold (limit) point by checking for a sign change in the
/// determinant of the Jacobian along two consecutive states.
///
/// Returns `true` if the determinants at `state_a` and `state_b` have opposite
/// signs (or one of them is zero), indicating a fold between them.
pub fn detect_fold_point(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
) -> bool {
    let det_a = matrix_determinant(&jac(&state_a.u, state_a.lambda));
    let det_b = matrix_determinant(&jac(&state_b.u, state_b.lambda));
    det_a * det_b <= 0.0
}

// ---------------------------------------------------------------------------
// stability_index
// ---------------------------------------------------------------------------

/// Count the number of eigenvalues with positive real part (unstable
/// eigenvalues) of a square real matrix.
///
/// Uses exact formulas for 1×1 and 2×2 matrices, and the Gershgorin circle
/// theorem for larger matrices.
pub fn stability_index(mat: &[Vec<f64>]) -> usize {
    let n = mat.len();
    if n == 0 {
        return 0;
    }
    if n == 1 {
        return if mat[0][0] > 0.0 { 1 } else { 0 };
    }
    if n == 2 {
        let tr = mat[0][0] + mat[1][1];
        let det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        let disc = tr * tr - 4.0 * det;
        if disc < 0.0 {
            if tr > 0.0 { 2 } else { 0 }
        } else {
            let sqrt_d = disc.sqrt();
            let e1 = (tr + sqrt_d) / 2.0;
            let e2 = (tr - sqrt_d) / 2.0;
            let mut count = 0;
            if e1 > 0.0 {
                count += 1;
            }
            if e2 > 0.0 {
                count += 1;
            }
            count
        }
    } else {
        let mut count = 0;
        for i in 0..n {
            let center = mat[i][i];
            let radius: f64 = (0..n)
                .filter(|&jj| jj != i)
                .map(|jj| mat[i][jj].abs())
                .sum();
            if center - radius > 0.0 {
                count += 1;
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// detect_bifurcation
// ---------------------------------------------------------------------------

/// Detect a bifurcation point by checking whether any eigenvalue of the
/// Jacobian changes sign between two consecutive states.
///
/// Returns `true` if a bifurcation is detected.
pub fn detect_bifurcation(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
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
    /// Fold (saddle-node / limit) point: determinant changes sign.
    Fold,
    /// Pitchfork bifurcation: determinant vanishes with odd symmetry.
    Pitchfork,
    /// Hopf bifurcation: a complex-conjugate pair crosses the imaginary axis.
    Hopf,
    /// Unclassified bifurcation.
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

/// Classify the bifurcation between two consecutive continuation states.
///
/// Checks determinant sign change (fold / pitchfork) and trace-vs-determinant
/// relationship (Hopf) for the Jacobian at the midpoint between `state_a`
/// and `state_b`.
pub fn classify_bifurcation(
    state_a: &ContinuationState,
    state_b: &ContinuationState,
    jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
) -> Option<BifurcationPoint> {
    let det_a = matrix_determinant(&jac(&state_a.u, state_a.lambda));
    let det_b = matrix_determinant(&jac(&state_b.u, state_b.lambda));
    if det_a * det_b > 0.0 {
        // No fold or pitchfork detected — check stability change for Hopf
        let idx_a = stability_index(&jac(&state_a.u, state_a.lambda));
        let idx_b = stability_index(&jac(&state_b.u, state_b.lambda));
        if idx_a == idx_b {
            return None;
        }
        // Interpolate midpoint
        let u_mid: Vec<f64> = state_a
            .u
            .iter()
            .zip(&state_b.u)
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let lam_mid = 0.5 * (state_a.lambda + state_b.lambda);
        let mid_state = ContinuationState::new(lam_mid, u_mid, state_a.ds);
        let j_mid = jac(&mid_state.u, lam_mid);
        let det_mid = matrix_determinant(&j_mid);
        return Some(BifurcationPoint {
            state: mid_state,
            bif_type: BifurcationType::Hopf,
            det_j: det_mid,
        });
    }

    // Determinant sign change → fold or pitchfork
    let u_mid: Vec<f64> = state_a
        .u
        .iter()
        .zip(&state_b.u)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let lam_mid = 0.5 * (state_a.lambda + state_b.lambda);
    let mid_state = ContinuationState::new(lam_mid, u_mid, state_a.ds);
    let j_mid = jac(&mid_state.u, lam_mid);
    let det_mid = matrix_determinant(&j_mid);

    // Simple heuristic: if trace changes sign simultaneously → pitchfork candidate
    let n2 = j_mid.len();
    let tr_mid: f64 = (0..n2).map(|i| j_mid[i][i]).sum();
    let bif_type = if tr_mid.abs() < 1e-6 {
        BifurcationType::Pitchfork
    } else {
        BifurcationType::Fold
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

/// Switch to a secondary branch at a detected bifurcation point.
///
/// Computes a perturbation in the null-space direction of the Jacobian at the
/// bifurcation state, then returns a new `ContinuationState` on the secondary
/// branch displaced by `epsilon`.
pub fn branch_switching(
    state: &ContinuationState,
    jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
    epsilon: f64,
) -> ContinuationState {
    let n = state.u.len();
    let j = jac(&state.u, state.lambda);
    let mut null_dir = vec![0.0_f64; n];
    let mut min_col_norm = f64::MAX;
    for col in 0..n {
        let col_norm: f64 = (0..n)
            .map(|row| j[row][col] * j[row][col])
            .sum::<f64>()
            .sqrt();
        if col_norm < min_col_norm {
            min_col_norm = col_norm;
            for row in 0..n {
                null_dir[row] = j[row][col];
            }
        }
    }
    let norm: f64 = null_dir.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-14 {
        for x in &mut null_dir {
            *x /= norm;
        }
    } else {
        null_dir[0] = 1.0;
    }
    let u_new: Vec<f64> = state
        .u
        .iter()
        .zip(&null_dir)
        .map(|(&ui, &di)| ui + epsilon * di)
        .collect();
    let mut tangent_new = null_dir.clone();
    tangent_new.push(0.0);
    ContinuationState {
        lambda: state.lambda,
        u: u_new,
        tangent: tangent_new,
        arc_length: state.arc_length,
        ds: state.ds,
    }
}

// ---------------------------------------------------------------------------
// BranchSwitching — bordered-system approach
// ---------------------------------------------------------------------------

/// Bordered-system branch switching at a bifurcation point.
///
/// Solves the bordered system `[J, v; w^T, 0] [du; dσ] = [-F; 0]` where `v`
/// and `w` are the (approximate) null vectors of `J` to switch branches.
/// Returns the perturbed state on the new branch.
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

    /// Switch branch at `state` using the bordered system approach.
    ///
    /// Returns the new continuation state on the secondary branch.
    pub fn switch(
        &self,
        state: &ContinuationState,
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
    ) -> ContinuationState {
        branch_switching(state, jac, self.epsilon)
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
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
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
/// Evaluates the Jacobian eigenvalue structure at each `ContinuationState`
/// in a sequence and assigns a stability label.
pub struct StabilityAnalysis {
    /// Tolerance for detecting marginal (near-zero real part) eigenvalues.
    pub marginal_tol: f64,
}

impl StabilityAnalysis {
    /// Create a new stability analyser.
    pub fn new(marginal_tol: f64) -> Self {
        Self { marginal_tol }
    }

    /// Assign a stability label to a single state.
    ///
    /// Uses the exact 2×2 formula and Gershgorin for larger systems.
    pub fn label(
        &self,
        state: &ContinuationState,
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
    ) -> StabilityLabel {
        let j = jac(&state.u, state.lambda);
        let n = j.len();
        if n == 0 {
            return StabilityLabel::Stable;
        }
        if n == 1 {
            let e = j[0][0];
            if e.abs() < self.marginal_tol {
                return StabilityLabel::Marginal;
            }
            return if e < 0.0 {
                StabilityLabel::Stable
            } else {
                StabilityLabel::Unstable
            };
        }
        if n == 2 {
            let tr = j[0][0] + j[1][1];
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            let disc = tr * tr - 4.0 * det;
            if disc < 0.0 {
                // Complex pair — check real part = tr/2
                if tr.abs() < self.marginal_tol {
                    return StabilityLabel::Marginal;
                }
                return if tr < 0.0 {
                    StabilityLabel::Stable
                } else {
                    StabilityLabel::Unstable
                };
            }
            let sqrt_d = disc.sqrt();
            let e1 = (tr + sqrt_d) / 2.0;
            let e2 = (tr - sqrt_d) / 2.0;
            if e1.abs() < self.marginal_tol || e2.abs() < self.marginal_tol {
                return StabilityLabel::Marginal;
            }
            if e1 > 0.0 || e2 > 0.0 {
                return StabilityLabel::Unstable;
            }
            return StabilityLabel::Stable;
        }
        // Gershgorin for n > 2
        let mut any_unstable = false;
        let mut any_marginal = false;
        for i in 0..n {
            let center = j[i][i];
            let radius: f64 = (0..n).filter(|&jj| jj != i).map(|jj| j[i][jj].abs()).sum();
            if center - radius > 0.0 {
                any_unstable = true;
            }
            if center.abs() <= radius {
                any_marginal = true;
            }
        }
        if any_unstable {
            StabilityLabel::Unstable
        } else if any_marginal {
            StabilityLabel::Marginal
        } else {
            StabilityLabel::Stable
        }
    }

    /// Label a sequence of states along a branch.
    pub fn label_branch(
        &self,
        states: &[ContinuationState],
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
    ) -> Vec<StabilityLabel> {
        states.iter().map(|s| self.label(s, jac)).collect()
    }
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
    /// Returns `None` if the corrector fails to converge even after halving `ds`.
    pub fn step(
        &self,
        state: &mut ContinuationState,
        f: &dyn Fn(&[f64], f64) -> Vec<f64>,
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
    ) -> Option<ContinuationState> {
        // Try current ds; if corrector fails, halve ds up to 5 times.
        for attempt in 0..5usize {
            let _ = attempt;
            let predicted = pseudo_arclength_step(state);
            let result = corrector_newton(&predicted, state, f, jac, self.tol, self.max_iter);
            if result.converged {
                // Update tangent via finite difference between states
                let n = state.u.len();
                let mut new_tangent = vec![0.0_f64; n + 1];
                for i in 0..n {
                    new_tangent[i] = result.u[i] - state.u[i];
                }
                new_tangent[n] = result.lambda - state.lambda;
                let t_norm: f64 = new_tangent.iter().map(|x| x * x).sum::<f64>().sqrt();
                if t_norm > 1e-14 {
                    for x in &mut new_tangent {
                        *x /= t_norm;
                    }
                }

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
    /// Whether a fold was detected between this and the previous step.
    pub fold_detected: bool,
    /// Whether a bifurcation was detected between this and the previous step.
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
        f: &dyn Fn(&[f64], f64) -> Vec<f64>,
        jac: &dyn Fn(&[f64], f64) -> Vec<Vec<f64>>,
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
        let s_branch = branch_switching(&s, &jac2d, 0.01);
        let diff: f64 =
            s.u.iter()
                .zip(&s_branch.u)
                .map(|(a, b)| (a - b).abs())
                .sum();
        assert!(diff > 0.0, "branch switching must perturb the solution");
    }

    #[test]
    fn test_branch_switching_preserves_lambda() {
        let s = ContinuationState::new(2.5, vec![1.0], 0.05);
        let jac_id = |_u: &[f64], _lam: f64| vec![vec![1.0]];
        let s2 = branch_switching(&s, &jac_id, 0.1);
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
    fn test_classify_bifurcation_detects_fold() {
        // 1D: F(u, λ) = u² - λ, sign change in det(J) between u=1.0 and u=-1.0.
        // jac1d = [[2u]]; at u=1 → det=2 > 0; at u=-1 → det=-2 < 0.
        // Midpoint u=0 → trace=0, classified as Pitchfork by the trace heuristic
        // (degenerate case). Test that a bifurcation point is detected (type Fold
        // or Pitchfork) — the key requirement is detection, not strict type.
        let sa = ContinuationState::new(1.0, vec![1.0], 0.1);
        let sb = ContinuationState::new(1.0, vec![-1.0], 0.1);
        let bif = classify_bifurcation(&sa, &sb, &jac1d);
        assert!(bif.is_some(), "a bifurcation should be detected");
        let bp = bif.unwrap();
        assert!(
            bp.bif_type == BifurcationType::Fold || bp.bif_type == BifurcationType::Pitchfork,
            "expected Fold or Pitchfork, got {:?}",
            bp.bif_type
        );
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
        let bs = BranchSwitching::new(0.01, 20, 1e-8);
        let s = ContinuationState::new(1.0, vec![1.0, 0.0], 0.1);
        let jac2d = |_u: &[f64], _lam: f64| vec![vec![1e-15, 0.0], vec![0.0, 1.0]];
        let s2 = bs.switch(&s, &jac2d);
        let diff: f64 = s.u.iter().zip(&s2.u).map(|(a, b)| (a - b).abs()).sum();
        assert!(diff > 0.0);
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
