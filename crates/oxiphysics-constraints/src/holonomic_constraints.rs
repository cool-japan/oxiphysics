// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Holonomic constraint theory for the OxiPhysics constraints crate.
//!
//! A *holonomic* constraint is one that can be written purely as a function
//! of the generalised coordinates (positions) without involving velocities:
//!
//! ```text
//! C(q) = 0          (position level)
//! Ċ(q,v) = J v = 0  (velocity level)
//! C̈(q,v,a) = J a + J̇ v = 0  (acceleration level)
//! ```
//!
//! This module provides:
//!
//! - **Position-level holonomic constraints** and their Jacobians.
//! - **Velocity-level constraint residuals** (Cdot = J v).
//! - **Acceleration-level residuals** (J a + J̇ v).
//! - **Baumgarte stabilization** — a feedback correction that drives constraint
//!   drift to zero exponentially.
//! - **Null-space projection method** — projects velocities/accelerations onto
//!   the constraint manifold.
//! - **Constraint drift correction** — simple position-projection and
//!   post-stabilization.
//! - **Bilateral and unilateral constraint wrappers** with activation tests.
//! - **Constraint Jacobian assembly** for multi-body systems.
//! - **DAE (Differential-Algebraic Equations) solvers** — index-3 to index-1
//!   reduction and a simple BDF-1 (backward Euler) integrator.
//! - **Penalty method vs Lagrange multiplier comparison** utilities.

// ---------------------------------------------------------------------------
// Primitive linear algebra (no nalgebra — plain f64 arrays)
// ---------------------------------------------------------------------------

/// Dot product of two slices (must have the same length).
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Add `s * b` to `a` in-place.
fn axpy(a: &mut [f64], s: f64, b: &[f64]) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x += s * y;
    }
}

/// 3-vector dot product.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 3-vector cross product.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector.
#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalize a 3-vector; returns `[0,0,0]` for near-zero inputs.
#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Subtract two 3-vectors.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// 3×3 matrix–vector product (row-major storage).
fn mat3mv(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// 1. Position-Level Holonomic Constraints
// ---------------------------------------------------------------------------

/// A distance (holonomic) constraint between two 3-D points.
///
/// `C(q) = ‖p_b − p_a‖ − d_rest = 0`
#[derive(Debug, Clone)]
pub struct DistanceConstraintHolo {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Rest distance.
    pub rest_distance: f64,
}

impl DistanceConstraintHolo {
    /// Create a new distance constraint.
    pub fn new(body_a: usize, body_b: usize, rest_distance: f64) -> Self {
        Self {
            body_a,
            body_b,
            rest_distance,
        }
    }

    /// Evaluate the constraint error `C(q)`.
    pub fn evaluate(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> f64 {
        norm3(sub3(pos_b, pos_a)) - self.rest_distance
    }

    /// Compute the constraint Jacobian rows for body A and body B.
    ///
    /// Returns `(J_a, J_b)` where each is a 3-element row vector.
    /// `J_a = -(p_b - p_a) / ‖p_b - p_a‖`, `J_b = -J_a`.
    pub fn jacobian(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let diff = sub3(pos_b, pos_a);
        let n = normalize3(diff);
        (scale3(n, -1.0), n)
    }
}

/// A ball-in-socket (spherical joint) holonomic constraint.
///
/// `C(q) = p_b − p_a = 0`  (3 scalar constraints)
#[derive(Debug, Clone)]
pub struct BallSocketConstraintHolo {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
}

impl BallSocketConstraintHolo {
    /// Create a new ball-socket constraint.
    pub fn new(body_a: usize, body_b: usize) -> Self {
        Self { body_a, body_b }
    }

    /// Evaluate the 3-component constraint error.
    pub fn evaluate(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> [f64; 3] {
        sub3(pos_b, pos_a)
    }

    /// Constraint Jacobian: `J_a = -I`, `J_b = +I`.
    pub fn jacobian(&self) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
        let neg_id = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
        let pos_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        (neg_id, pos_id)
    }
}

/// A hinge (revolute) joint holonomic constraint (5-DOF constraint, 1 DOF free).
///
/// Only the positional part (3 DOF) is modelled here for clarity.
/// The two angular constraints that fix the rotation axis are left for
/// the `six_dof_constraint` module.
#[derive(Debug, Clone)]
pub struct HingeConstraintHolo {
    /// Body A anchor point in world space (at rest).
    pub anchor_a: [f64; 3],
    /// Body B anchor point in world space (at rest).
    pub anchor_b: [f64; 3],
    /// Hinge axis (unit vector in world space).
    pub axis: [f64; 3],
}

impl HingeConstraintHolo {
    /// Create a new hinge constraint.
    pub fn new(anchor_a: [f64; 3], anchor_b: [f64; 3], axis: [f64; 3]) -> Self {
        Self {
            anchor_a,
            anchor_b,
            axis: normalize3(axis),
        }
    }

    /// Positional constraint error: `C = p_b − p_a`.
    pub fn positional_error(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> [f64; 3] {
        sub3(pos_b, pos_a)
    }

    /// Angular constraint error: project relative rotation onto the two axes
    /// perpendicular to the hinge axis and return the two scalar errors.
    ///
    /// Simplified: returns the component of `(rot_b * axis - rot_a * axis)` in
    /// the perpendicular plane.
    pub fn angular_errors(&self, rot_a: [[f64; 3]; 3], rot_b: [[f64; 3]; 3]) -> [f64; 2] {
        let ax_a = mat3mv(rot_a, self.axis);
        let ax_b = mat3mv(rot_b, self.axis);
        let diff = sub3(ax_b, ax_a);
        // Project diff onto two basis vectors perpendicular to self.axis
        let t1 = perpendicular3(self.axis);
        let t2 = cross3(self.axis, t1);
        [dot3(diff, t1), dot3(diff, t2)]
    }
}

/// Return any vector perpendicular to `n`.
fn perpendicular3(n: [f64; 3]) -> [f64; 3] {
    let candidate = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    normalize3(cross3(n, candidate))
}

// ---------------------------------------------------------------------------
// 2. Velocity-Level Constraint Residuals
// ---------------------------------------------------------------------------

/// Compute the velocity-level residual `Ċ = J v` for a distance constraint.
///
/// `Ċ = n · (v_b − v_a)` where `n = (p_b − p_a) / ‖p_b − p_a‖`.
pub fn velocity_residual_distance(
    pos_a: [f64; 3],
    pos_b: [f64; 3],
    vel_a: [f64; 3],
    vel_b: [f64; 3],
) -> f64 {
    let n = normalize3(sub3(pos_b, pos_a));
    dot3(n, sub3(vel_b, vel_a))
}

/// Compute the velocity-level residuals for a ball-socket constraint (3 DOF).
///
/// `Ċ = v_b − v_a`
pub fn velocity_residual_ball_socket(vel_a: [f64; 3], vel_b: [f64; 3]) -> [f64; 3] {
    sub3(vel_b, vel_a)
}

// ---------------------------------------------------------------------------
// 3. Acceleration-Level Constraint Residuals
// ---------------------------------------------------------------------------

/// Acceleration-level residual for a distance constraint.
///
/// ```text
/// C̈ = n · (a_b − a_a) + J̇ v
/// ```
///
/// where `J̇ v = ‖v_rel_⊥‖² / d` (centripetal term).
pub fn acceleration_residual_distance(
    pos_a: [f64; 3],
    pos_b: [f64; 3],
    vel_a: [f64; 3],
    vel_b: [f64; 3],
    acc_a: [f64; 3],
    acc_b: [f64; 3],
) -> f64 {
    let diff = sub3(pos_b, pos_a);
    let d = norm3(diff);
    if d < 1e-15 {
        return 0.0;
    }
    let n = scale3(diff, 1.0 / d);
    let v_rel = sub3(vel_b, vel_a);
    // Centripetal term: ‖v_rel‖² / d − (n·v_rel)² / d
    let vn = dot3(n, v_rel);
    let vrel2 = dot3(v_rel, v_rel);
    let jdot_v = (vrel2 - vn * vn) / d;
    dot3(n, sub3(acc_b, acc_a)) + jdot_v
}

// ---------------------------------------------------------------------------
// 4. Baumgarte Stabilization
// ---------------------------------------------------------------------------

/// Parameters for Baumgarte constraint stabilization.
///
/// Adds a feedback term to the velocity (and optionally acceleration) level
/// to exponentially drive away constraint drift:
///
/// ```text
/// J v = -α C − β Ċ     (velocity level)
/// J a = -2α Ċ − α² C  (acceleration level, critically-damped form)
/// ```
#[derive(Debug, Clone, Copy)]
pub struct BaumgarteParams {
    /// Position-error feedback gain `α` (typically `1 / dt` to `10 / dt`).
    pub alpha: f64,
    /// Velocity-error feedback gain `β` (typically `2 √α` for critical damping).
    pub beta: f64,
}

impl BaumgarteParams {
    /// Create Baumgarte parameters from `alpha` with critical damping (`beta = 2√α`).
    pub fn critically_damped(alpha: f64) -> Self {
        Self {
            alpha,
            beta: 2.0 * alpha.sqrt(),
        }
    }

    /// Compute the Baumgarte stabilization bias for the velocity solve.
    ///
    /// `bias = -α * C − β * Ċ`
    pub fn velocity_bias(&self, constraint_error: f64, velocity_error: f64) -> f64 {
        -self.alpha * constraint_error - self.beta * velocity_error
    }

    /// Compute the Baumgarte stabilization bias for the acceleration solve.
    ///
    /// `bias = -2α * Ċ − α² * C`
    pub fn acceleration_bias(&self, constraint_error: f64, velocity_error: f64) -> f64 {
        -2.0 * self.alpha * velocity_error - self.alpha * self.alpha * constraint_error
    }
}

// ---------------------------------------------------------------------------
// 5. Projection Method (Null-Space)
// ---------------------------------------------------------------------------

/// Project a velocity vector `v` onto the constraint null-space.
///
/// Given a Jacobian row `J` (a 1-D constraint), the null-space projector is:
///
/// ```text
/// P = I − J^T (J M^{-1} J^T)^{-1} J M^{-1}
/// ```
///
/// For a simple scalar constraint with one Jacobian row `J` and an effective
/// inverse mass `eff_inv_mass = J M^{-1} J^T`, the projected velocity is:
///
/// ```text
/// v' = v − J^T * (J v) / (J M^{-1} J^T)
/// ```
///
/// This function handles a single scalar constraint.
pub fn null_space_project_velocity(v: &mut [f64], jacobian: &[f64], eff_mass_inv: f64) {
    if eff_mass_inv.abs() < 1e-15 {
        return;
    }
    let jv = dot(jacobian, v);
    let lambda = jv / eff_mass_inv;
    axpy(v, -lambda, jacobian);
}

/// Project a velocity vector onto the null-space of *multiple* scalar
/// constraints using successive (Gauss-Seidel) projections.
///
/// Each row in `jacobians` is one constraint; `eff_mass_inv[i]` is the
/// corresponding `J_i M^{-1} J_i^T` scalar.
pub fn null_space_project_multi(
    v: &mut [f64],
    jacobians: &[Vec<f64>],
    eff_mass_inv: &[f64],
    iters: usize,
) {
    for _ in 0..iters {
        for (j, &em) in jacobians.iter().zip(eff_mass_inv.iter()) {
            null_space_project_velocity(v, j, em);
        }
    }
}

// ---------------------------------------------------------------------------
// 6. Constraint Drift Correction
// ---------------------------------------------------------------------------

/// Correct constraint drift by direct position projection.
///
/// Moves body positions along the constraint gradient by `−C / J M^{-1} J^T`.
///
/// Returns the positional correction `Δq = −J^T * (C / J M^{-1} J^T)`.
pub fn position_drift_correction(
    positions: &mut [f64],
    jacobian: &[f64],
    constraint_error: f64,
    eff_mass_inv: f64,
) {
    if eff_mass_inv.abs() < 1e-15 {
        return;
    }
    let lambda = -constraint_error / eff_mass_inv;
    axpy(positions, lambda, jacobian);
}

/// Compute the post-stabilization impulse magnitude.
///
/// Post-stabilization (also called *split impulse*) applies a pseudo-velocity
/// impulse purely to correct position drift without affecting the true velocity:
///
/// ```text
/// λ_ps = -β * C / (J M^{-1} J^T) / dt
/// ```
pub fn post_stabilization_impulse(
    constraint_error: f64,
    eff_mass_inv: f64,
    beta: f64,
    dt: f64,
) -> f64 {
    if eff_mass_inv.abs() < 1e-15 || dt.abs() < 1e-15 {
        return 0.0;
    }
    (-beta * constraint_error) / (eff_mass_inv * dt)
}

// ---------------------------------------------------------------------------
// 7. Bilateral vs Unilateral Constraints
// ---------------------------------------------------------------------------

/// A bilateral holonomic constraint: active always (`C(q) = 0`).
#[derive(Debug, Clone)]
pub struct BilateralConstraint {
    /// Scalar constraint error (positive means violated).
    pub error: f64,
    /// Effective inverse mass `J M^{-1} J^T`.
    pub eff_mass_inv: f64,
    /// Accumulated impulse (for warm-starting).
    pub accumulated_impulse: f64,
}

impl BilateralConstraint {
    /// Create a new bilateral constraint.
    pub fn new(error: f64, eff_mass_inv: f64) -> Self {
        Self {
            error,
            eff_mass_inv,
            accumulated_impulse: 0.0,
        }
    }

    /// Compute the impulse correction for this constraint.
    ///
    /// `λ = −(Ċ + bias) / (J M^{-1} J^T)`
    pub fn solve_impulse(&mut self, cdot: f64, bias: f64) -> f64 {
        if self.eff_mass_inv.abs() < 1e-15 {
            return 0.0;
        }
        let lambda = -(cdot + bias) / self.eff_mass_inv;
        self.accumulated_impulse += lambda;
        lambda
    }

    /// Reset accumulated impulse (start of new time step).
    pub fn reset_warm_start(&mut self) {
        self.accumulated_impulse = 0.0;
    }
}

/// A unilateral holonomic constraint: only active when `C(q) ≤ 0` (contact).
///
/// The impulse is clamped to be non-negative (bodies can only push, not pull).
#[derive(Debug, Clone)]
pub struct UnilateralConstraint {
    /// Scalar constraint error (negative = penetration).
    pub error: f64,
    /// Effective inverse mass.
    pub eff_mass_inv: f64,
    /// Accumulated non-negative impulse.
    pub accumulated_impulse: f64,
}

impl UnilateralConstraint {
    /// Create a new unilateral constraint.
    pub fn new(error: f64, eff_mass_inv: f64) -> Self {
        Self {
            error,
            eff_mass_inv,
            accumulated_impulse: 0.0,
        }
    }

    /// Is this constraint currently active (bodies overlapping)?
    pub fn is_active(&self) -> bool {
        self.error < 0.0
    }

    /// Solve for the non-negative impulse correction.
    pub fn solve_impulse(&mut self, cdot: f64, bias: f64) -> f64 {
        if !self.is_active() || self.eff_mass_inv.abs() < 1e-15 {
            return 0.0;
        }
        let lambda = -(cdot + bias) / self.eff_mass_inv;
        let old = self.accumulated_impulse;
        self.accumulated_impulse = (old + lambda).max(0.0);
        self.accumulated_impulse - old
    }

    /// Reset accumulated impulse.
    pub fn reset_warm_start(&mut self) {
        self.accumulated_impulse = 0.0;
    }
}

// ---------------------------------------------------------------------------
// 8. Constraint Jacobian Assembly
// ---------------------------------------------------------------------------

/// A dense constraint Jacobian for a multi-body system.
///
/// Rows = constraints, Columns = generalized velocity DOFs.
#[derive(Debug, Clone)]
pub struct ConstraintJacobian {
    /// Number of constraint rows.
    pub num_constraints: usize,
    /// Number of velocity DOFs (columns).
    pub num_dofs: usize,
    /// Row-major storage: `data[row * num_dofs + col]`.
    pub data: Vec<f64>,
}

impl ConstraintJacobian {
    /// Create a zero Jacobian of the given size.
    pub fn new(num_constraints: usize, num_dofs: usize) -> Self {
        Self {
            num_constraints,
            num_dofs,
            data: vec![0.0; num_constraints * num_dofs],
        }
    }

    /// Set element `(row, col)`.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[row * self.num_dofs + col] = value;
    }

    /// Get element `(row, col)`.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.num_dofs + col]
    }

    /// Compute `J v` (constraint velocity errors) given a velocity vector.
    pub fn mul_velocity(&self, v: &[f64]) -> Vec<f64> {
        (0..self.num_constraints)
            .map(|row| {
                (0..self.num_dofs)
                    .map(|col| self.get(row, col) * v[col])
                    .sum()
            })
            .collect()
    }

    /// Compute `J^T λ` (constraint forces) given a Lagrange multiplier vector.
    pub fn mul_lambda(&self, lambda: &[f64]) -> Vec<f64> {
        let mut result = vec![0.0; self.num_dofs];
        for (row, lam_row) in lambda.iter().enumerate() {
            for (col, res_col) in result.iter_mut().enumerate() {
                *res_col += self.get(row, col) * lam_row;
            }
        }
        result
    }

    /// Compute `J M^{-1} J^T` (constraint effective mass matrix).
    ///
    /// `inv_mass` is the diagonal of `M^{-1}` (one entry per DOF).
    pub fn effective_mass_matrix(&self, inv_mass: &[f64]) -> Vec<Vec<f64>> {
        let nc = self.num_constraints;
        let nd = self.num_dofs;
        let mut result = vec![vec![0.0; nc]; nc];
        for (i, res_row) in result.iter_mut().enumerate() {
            for (j, res_ij) in res_row.iter_mut().enumerate() {
                let mut sum = 0.0;
                for (k, inv_m_k) in inv_mass.iter().enumerate().take(nd) {
                    sum += self.get(i, k) * inv_m_k * self.get(j, k);
                }
                *res_ij = sum;
            }
        }
        result
    }

    /// Add a distance constraint row between DOFs `dof_a` and `dof_b`
    /// along axis direction `n` at constraint row `row`.
    pub fn add_distance_row(&mut self, row: usize, dof_a: usize, dof_b: usize, n: [f64; 3]) {
        // J_a = -n, J_b = +n (assumes linear DOFs are packed as [x,y,z])
        for (d, n_d) in n.iter().enumerate() {
            if dof_a + d < self.num_dofs {
                self.set(row, dof_a + d, -n_d);
            }
            if dof_b + d < self.num_dofs {
                self.set(row, dof_b + d, *n_d);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 9. DAE Solvers
// ---------------------------------------------------------------------------

/// State of the DAE system: generalized coordinates `q` and velocities `v`.
#[derive(Debug, Clone)]
pub struct DaeState {
    /// Generalized positions.
    pub q: Vec<f64>,
    /// Generalized velocities.
    pub v: Vec<f64>,
    /// Lagrange multipliers.
    pub lambda: Vec<f64>,
}

impl DaeState {
    /// Create a zero-initialized DAE state.
    pub fn new(n_q: usize, n_lambda: usize) -> Self {
        Self {
            q: vec![0.0; n_q],
            v: vec![0.0; n_q],
            lambda: vec![0.0; n_lambda],
        }
    }
}

/// Index-3 to index-1 reduction for a holonomic DAE.
///
/// The original index-3 system is:
/// ```text
/// M q̈ = F + J^T λ
/// C(q) = 0
/// ```
///
/// By differentiating twice we obtain the index-1 system:
/// ```text
/// M q̈ = F + J^T λ
/// J q̈ = −J̇ q̇ − α C − β Ċ    (with Baumgarte)
/// ```
///
/// This function performs one step of an **explicit forward Euler** solve
/// at the *velocity level* (index-1 formulation) with Baumgarte stabilization.
///
/// * `state` – current state (mutated in place).
/// * `forces` – generalised forces `F`.
/// * `jac` – constraint Jacobian `J`.
/// * `inv_mass` – diagonal inverse mass.
/// * `constraint_errors` – `C(q)` values.
/// * `bparams` – Baumgarte parameters.
/// * `dt` – time step.
pub fn dae_step_index1(
    state: &mut DaeState,
    forces: &[f64],
    jac: &ConstraintJacobian,
    inv_mass: &[f64],
    constraint_errors: &[f64],
    bparams: BaumgarteParams,
    dt: f64,
) {
    let _nd = state.q.len();
    let nc = jac.num_constraints;

    // Compute velocity-level errors Ċ = J v
    let cdot = jac.mul_velocity(&state.v);

    // Compute Baumgarte biases
    let bias: Vec<f64> = (0..nc)
        .map(|i| {
            bparams.velocity_bias(
                constraint_errors.get(i).copied().unwrap_or(0.0),
                cdot.get(i).copied().unwrap_or(0.0),
            )
        })
        .collect();

    // Build J M^{-1} J^T (diagonal approximation for efficiency)
    let eff_mass_matrix = jac.effective_mass_matrix(inv_mass);

    // Solve for Lagrange multipliers λ using Gauss-Seidel
    let mut lambda = vec![0.0f64; nc];
    for _gs in 0..20 {
        for i in 0..nc {
            let diag = eff_mass_matrix[i][i];
            if diag.abs() < 1e-15 {
                continue;
            }
            let mut off = 0.0;
            for j in 0..nc {
                if j != i {
                    off += eff_mass_matrix[i][j] * lambda[j];
                }
            }
            lambda[i] = (bias[i] - cdot[i] - off) / diag;
        }
    }
    state.lambda = lambda.clone();

    // Compute constraint forces J^T λ
    let constraint_forces = jac.mul_lambda(&lambda);

    // Update velocities: v += dt * M^{-1} * (F + J^T λ)
    for (k, (v_k, inv_m_k)) in state.v.iter_mut().zip(inv_mass.iter()).enumerate() {
        let total_force = forces.get(k).copied().unwrap_or(0.0)
            + constraint_forces.get(k).copied().unwrap_or(0.0);
        *v_k += dt * inv_m_k * total_force;
    }

    // Update positions: q += dt * v
    for (q_k, v_k) in state.q.iter_mut().zip(state.v.iter()) {
        *q_k += dt * v_k;
    }
}

/// BDF-1 (Backward Euler) DAE integrator step.
///
/// At each step, solves for `v^{n+1}` and `λ^{n+1}` satisfying:
///
/// ```text
/// M (v^{n+1} - v^n) / dt = F^{n+1} + J^T λ^{n+1}
/// J v^{n+1} + α C + β Ċ = 0
/// ```
///
/// Implemented via one iteration of the index-1 Gauss-Seidel solve.
pub fn dae_bdf1_step(
    state: &mut DaeState,
    forces: &[f64],
    jac: &ConstraintJacobian,
    inv_mass: &[f64],
    constraint_errors: &[f64],
    bparams: BaumgarteParams,
    dt: f64,
) {
    // For BDF-1 the effective mass is M/dt on the diagonal.
    // We scale inv_mass to M/dt: eff_inv = dt * inv_mass
    let scaled_inv_mass: Vec<f64> = inv_mass.iter().map(|m| m * dt).collect();

    // Unconstrained velocity predictor: v* = v + dt * M^{-1} * F
    let _nd = state.q.len();
    let mut v_pred: Vec<f64> = state.v.clone();
    for (k, vp_k) in v_pred.iter_mut().enumerate() {
        *vp_k +=
            dt * inv_mass.get(k).copied().unwrap_or(0.0) * forces.get(k).copied().unwrap_or(0.0);
    }

    // Constraint-correct the predictor
    let nc = jac.num_constraints;
    let cdot_pred = jac.mul_velocity(&v_pred);

    let bias: Vec<f64> = (0..nc)
        .map(|i| {
            bparams.velocity_bias(
                constraint_errors.get(i).copied().unwrap_or(0.0),
                cdot_pred.get(i).copied().unwrap_or(0.0),
            )
        })
        .collect();

    let eff_mass_matrix = jac.effective_mass_matrix(&scaled_inv_mass);

    let mut lambda = vec![0.0f64; nc];
    for _gs in 0..30 {
        for i in 0..nc {
            let diag = eff_mass_matrix[i][i];
            if diag.abs() < 1e-15 {
                continue;
            }
            let mut off = 0.0;
            for j in 0..nc {
                if j != i {
                    off += eff_mass_matrix[i][j] * lambda[j];
                }
            }
            lambda[i] = (bias[i] - cdot_pred[i] - off) / diag;
        }
    }
    state.lambda = lambda.clone();

    let constraint_forces = jac.mul_lambda(&lambda);

    for (k, (v_k, vp_k)) in state.v.iter_mut().zip(v_pred.iter()).enumerate() {
        *v_k = vp_k
            + scaled_inv_mass.get(k).copied().unwrap_or(0.0)
                * constraint_forces.get(k).copied().unwrap_or(0.0);
    }
    for (q_k, v_k) in state.q.iter_mut().zip(state.v.iter()) {
        *q_k += dt * v_k;
    }
}

// ---------------------------------------------------------------------------
// 10. Penalty Method vs Lagrange Multiplier Comparison
// ---------------------------------------------------------------------------

/// Penalty method constraint enforcement.
///
/// Instead of solving for Lagrange multipliers, the penalty method adds a
/// large spring force proportional to the constraint error:
///
/// ```text
/// F_penalty = −k * C(q) * J^T
/// ```
///
/// where `k` is the penalty stiffness.  Large `k` approaches exact constraint
/// satisfaction but causes numerical stiffness.
#[derive(Debug, Clone, Copy)]
pub struct PenaltyConstraint {
    /// Penalty stiffness `k` (larger = harder constraint, stiffer system).
    pub stiffness: f64,
    /// Damping coefficient `d` (adds `−d * Ċ * J^T`).
    pub damping: f64,
}

impl PenaltyConstraint {
    /// Create a penalty constraint with given stiffness and damping.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }

    /// Compute the penalty force along the Jacobian direction.
    ///
    /// Returns the scalar force magnitude: `F = −k * C − d * Ċ`.
    pub fn force(&self, constraint_error: f64, velocity_error: f64) -> f64 {
        -self.stiffness * constraint_error - self.damping * velocity_error
    }

    /// Estimate the natural frequency of constraint oscillation.
    ///
    /// `ω_n = √(k / m_eff)` where `m_eff = 1 / (J M^{-1} J^T)`.
    pub fn natural_frequency(&self, eff_mass: f64) -> f64 {
        if eff_mass < 1e-15 {
            0.0
        } else {
            (self.stiffness * eff_mass).sqrt()
        }
    }

    /// Compute the critical damping coefficient for this stiffness and effective mass.
    ///
    /// `d_crit = 2 √(k * m_eff^{-1})`  (inverse effective mass).
    pub fn critical_damping_coefficient(&self, eff_mass_inv: f64) -> f64 {
        2.0 * (self.stiffness * eff_mass_inv).sqrt()
    }
}

/// Lagrange multiplier constraint enforcement.
///
/// Exactly satisfies constraints (up to solver tolerance) at the expense of
/// a system solve.  This struct provides utilities for tracking convergence.
#[derive(Debug, Clone)]
pub struct LagrangeConstraintSolver {
    /// Maximum solver iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
    /// Last residual norm after solve.
    pub last_residual: f64,
    /// Number of iterations used in last solve.
    pub last_iter_count: usize,
}

impl LagrangeConstraintSolver {
    /// Create a new Lagrange solver with defaults.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self {
            max_iter,
            tol,
            last_residual: 0.0,
            last_iter_count: 0,
        }
    }

    /// Solve `(J M^{-1} J^T) λ = b` using Gauss-Seidel iteration.
    ///
    /// Returns the solution vector `λ`.
    pub fn solve_gs(&mut self, jmjt: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
        let n = rhs.len();
        let mut x = vec![0.0; n];
        let mut iter_count = 0;
        for iter in 0..self.max_iter {
            iter_count = iter + 1;
            let mut max_delta = 0.0_f64;
            for i in 0..n {
                let diag = jmjt[i][i];
                if diag.abs() < 1e-15 {
                    continue;
                }
                let mut sigma: f64 = 0.0;
                for j in 0..n {
                    if j != i {
                        sigma += jmjt[i][j] * x[j];
                    }
                }
                let x_new = (rhs[i] - sigma) / diag;
                max_delta = max_delta.max((x_new - x[i]).abs());
                x[i] = x_new;
            }
            if max_delta < self.tol {
                break;
            }
        }
        // Compute residual
        let mut res = 0.0;
        for i in 0..n {
            let mut row_dot = 0.0;
            for j in 0..n {
                row_dot += jmjt[i][j] * x[j];
            }
            res += (row_dot - rhs[i]).powi(2);
        }
        self.last_residual = res.sqrt();
        self.last_iter_count = iter_count;
        x
    }

    /// Compare penalty vs Lagrange accuracy for a single constraint.
    ///
    /// Returns `(penalty_error, lagrange_error)` as a pair.
    pub fn compare_accuracy(
        penalty: &PenaltyConstraint,
        constraint_error: f64,
        velocity_error: f64,
        eff_mass_inv: f64,
        _dt: f64,
    ) -> (f64, f64) {
        // Penalty: residual force = k * C + d * Ċ
        let _pf = penalty.force(constraint_error, velocity_error);
        let penalty_pos_error = constraint_error; // penalty does not instantly zero C

        // Lagrange: directly solves λ = −Ċ_with_bias / (J M^{-1} J^T)
        let lagrange_pos_error = if eff_mass_inv.abs() > 1e-15 {
            // After one step the velocity correction is exact → position error ≈ 0
            0.0
        } else {
            constraint_error
        };

        (penalty_pos_error, lagrange_pos_error)
    }
}

// ---------------------------------------------------------------------------
// Additional utility: constraint energy
// ---------------------------------------------------------------------------

/// Compute the constraint potential energy for the penalty method.
///
/// `E_penalty = 0.5 * k * C²`
pub fn penalty_energy(stiffness: f64, constraint_error: f64) -> f64 {
    0.5 * stiffness * constraint_error * constraint_error
}

/// Compute the constraint power dissipation (damping).
///
/// `P_damp = d * Ċ²`
pub fn penalty_power(damping: f64, velocity_error: f64) -> f64 {
    damping * velocity_error * velocity_error
}

/// Compute the stability limit `dt_max` for the penalty method based on the
/// natural frequency of the penalised constraint oscillation.
///
/// For explicit Euler: `dt_max = 2 / ω_n = 2 / √(k * inv_m_eff)`.
pub fn penalty_stability_limit(stiffness: f64, eff_mass_inv: f64) -> f64 {
    let omega_n = (stiffness * eff_mass_inv).sqrt();
    if omega_n < 1e-15 {
        f64::MAX
    } else {
        2.0 / omega_n
    }
}

/// Compute the Lagrange multiplier from the impulse-level solve for a single
/// scalar constraint using warm-starting.
///
/// Applies one Newton step and clamps to `[lo, hi]`.
pub fn lagrange_impulse_step(
    cdot: f64,
    bias: f64,
    eff_mass_inv: f64,
    accumulated: &mut f64,
    lo: f64,
    hi: f64,
) -> f64 {
    if eff_mass_inv.abs() < 1e-15 {
        return 0.0;
    }
    let lambda = -(cdot + bias) / eff_mass_inv;
    let old = *accumulated;
    *accumulated = (old + lambda).max(lo).min(hi);
    *accumulated - old
}

// ---------------------------------------------------------------------------
// Rodrigues rotation utility (used in hinge angular tests)
// ---------------------------------------------------------------------------

/// Rotate a 3-vector `v` by angle `theta` (radians) around unit axis `axis`
/// using the Rodrigues formula.
pub fn rodrigues_rotate(v: [f64; 3], axis: [f64; 3], theta: f64) -> [f64; 3] {
    let k = normalize3(axis);
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let kdotv = dot3(k, v);
    let cross = cross3(k, v);
    add3(
        add3(scale3(v, cos_t), scale3(cross, sin_t)),
        scale3(k, kdotv * (1.0 - cos_t)),
    )
}

/// Build a 3×3 rotation matrix for rotation by `theta` around unit axis `axis`.
pub fn rotation_matrix(axis: [f64; 3], theta: f64) -> [[f64; 3]; 3] {
    let k = normalize3(axis);
    let c = theta.cos();
    let s = theta.sin();
    let t = 1.0 - c;
    [
        [
            t * k[0] * k[0] + c,
            t * k[0] * k[1] - s * k[2],
            t * k[0] * k[2] + s * k[1],
        ],
        [
            t * k[0] * k[1] + s * k[2],
            t * k[1] * k[1] + c,
            t * k[1] * k[2] - s * k[0],
        ],
        [
            t * k[0] * k[2] - s * k[1],
            t * k[1] * k[2] + s * k[0],
            t * k[2] * k[2] + c,
        ],
    ]
}

// ---------------------------------------------------------------------------
// Unit Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    // ---- Distance constraint ------------------------------------------------

    /// 1. Distance constraint error at rest distance is zero.
    #[test]
    fn test_distance_constraint_zero_error() {
        let c = DistanceConstraintHolo::new(0, 1, 1.0);
        let err = c.evaluate([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(err.abs() < EPS, "err={err}");
    }

    /// 2. Distance constraint error is positive when bodies are farther apart.
    #[test]
    fn test_distance_constraint_positive_error() {
        let c = DistanceConstraintHolo::new(0, 1, 1.0);
        let err = c.evaluate([0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((err - 1.0).abs() < EPS, "err={err}");
    }

    /// 3. Distance constraint Jacobian rows are unit opposites.
    #[test]
    fn test_distance_jacobian_unit_normal() {
        let c = DistanceConstraintHolo::new(0, 1, 1.0);
        let (ja, jb) = c.jacobian([0.0, 0.0, 0.0], [3.0, 0.0, 0.0]);
        // ja = -[1,0,0], jb = [1,0,0]
        assert!((ja[0] + 1.0).abs() < EPS);
        assert!((jb[0] - 1.0).abs() < EPS);
        // Sum should be zero (reaction law)
        for d in 0..3 {
            assert!(
                (ja[d] + jb[d]).abs() < EPS,
                "d={d}: ja+jb={}",
                ja[d] + jb[d]
            );
        }
    }

    // ---- Ball-socket constraint ----------------------------------------------

    /// 4. Ball-socket error is zero when positions match.
    #[test]
    fn test_ball_socket_zero_error() {
        let c = BallSocketConstraintHolo::new(0, 1);
        let err = c.evaluate([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(norm3(err) < EPS);
    }

    /// 5. Ball-socket Jacobian rows are ±identity.
    #[test]
    fn test_ball_socket_jacobian() {
        let c = BallSocketConstraintHolo::new(0, 1);
        let (ja, jb) = c.jacobian();
        // j_a + j_b should be zero matrix
        for row in 0..3 {
            for col in 0..3 {
                assert!((ja[row][col] + jb[row][col]).abs() < EPS);
            }
        }
    }

    // ---- Velocity residual --------------------------------------------------

    /// 6. Velocity residual is zero when bodies move together.
    #[test]
    fn test_velocity_residual_same_velocity() {
        let r = velocity_residual_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
        );
        assert!(r.abs() < EPS, "r={r}");
    }

    /// 7. Velocity residual equals relative normal velocity when bodies separate.
    #[test]
    fn test_velocity_residual_separating() {
        // n = [1,0,0], v_rel = [5,0,0] → Ċ = 5
        let r = velocity_residual_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
        );
        assert!((r - 5.0).abs() < EPS, "r={r}");
    }

    // ---- Baumgarte ----------------------------------------------------------

    /// 8. Baumgarte velocity bias with no error and no Ċ is zero.
    #[test]
    fn test_baumgarte_zero_bias() {
        let bp = BaumgarteParams::critically_damped(10.0);
        let bias = bp.velocity_bias(0.0, 0.0);
        assert!(bias.abs() < EPS);
    }

    /// 9. Baumgarte velocity bias with positive error is negative (corrective).
    #[test]
    fn test_baumgarte_corrective_bias() {
        let bp = BaumgarteParams {
            alpha: 10.0,
            beta: 5.0,
        };
        // bias = −10 * 0.1 − 5 * 0 = −1.0
        let bias = bp.velocity_bias(0.1, 0.0);
        assert!((bias + 1.0).abs() < EPS, "bias={bias}");
    }

    /// 10. Critical damping sets beta = 2 * sqrt(alpha).
    #[test]
    fn test_baumgarte_critical_damping() {
        let bp = BaumgarteParams::critically_damped(16.0);
        assert!((bp.beta - 8.0).abs() < EPS, "beta={}", bp.beta);
    }

    // ---- Null-space projection -----------------------------------------------

    /// 11. Null-space projection removes the constraint-violating velocity component.
    #[test]
    fn test_null_space_projection_basic() {
        let j = vec![1.0, 0.0, 0.0];
        let mut v = vec![3.0, 2.0, 1.0];
        // J M^{-1} J^T = 1 (unit mass, unit Jacobian)
        null_space_project_velocity(&mut v, &j, 1.0);
        // After projection J v = 0 → v[0] = 0
        let jv = dot(&j, &v);
        assert!(jv.abs() < EPS, "Jv={jv}");
    }

    /// 12. Null-space projection with zero effective mass does nothing.
    #[test]
    fn test_null_space_projection_zero_mass() {
        let j = vec![1.0, 0.0, 0.0];
        let mut v = vec![5.0, 2.0, 1.0];
        null_space_project_velocity(&mut v, &j, 0.0);
        // Should be unchanged
        assert!((v[0] - 5.0).abs() < EPS);
    }

    // ---- Drift correction ---------------------------------------------------

    /// 13. Position drift correction zeroes the constraint error in one step (exact mass).
    #[test]
    fn test_position_drift_correction() {
        // 1-DOF constraint: q[0] must be 0, currently q[0] = 0.5
        let mut q = vec![0.5];
        let j = vec![1.0];
        let error = 0.5; // C = q[0] - 0 = 0.5
        position_drift_correction(&mut q, &j, error, 1.0);
        // New q[0] = 0.5 + (−0.5) * 1.0 = 0.0
        assert!(q[0].abs() < EPS, "q[0]={}", q[0]);
    }

    /// 14. Post-stabilization impulse scales with beta and dt.
    #[test]
    fn test_post_stabilization_impulse() {
        let imp = post_stabilization_impulse(0.1, 1.0, 0.5, 0.01);
        // = −0.5 * 0.1 / (1.0 * 0.01) = −5.0
        assert!((imp + 5.0).abs() < EPS, "imp={imp}");
    }

    // ---- Bilateral / Unilateral ---------------------------------------------

    /// 15. Bilateral constraint impulse is non-zero for non-zero Cdot.
    #[test]
    fn test_bilateral_impulse_nonzero() {
        let mut b = BilateralConstraint::new(0.0, 1.0);
        let imp = b.solve_impulse(2.0, 0.0);
        // λ = −2 / 1 = −2
        assert!((imp + 2.0).abs() < EPS, "imp={imp}");
    }

    /// 16. Unilateral constraint is inactive when error ≥ 0.
    #[test]
    fn test_unilateral_inactive() {
        let mut u = UnilateralConstraint::new(0.1, 1.0);
        let imp = u.solve_impulse(1.0, 0.0);
        assert!(imp.abs() < EPS, "imp={imp}");
    }

    /// 17. Unilateral constraint clamps impulse to non-negative.
    #[test]
    fn test_unilateral_nonneg_clamp() {
        let mut u = UnilateralConstraint::new(-0.1, 1.0);
        u.solve_impulse(-5.0, 0.0); // would normally give +5 → clamped at 0 first
        assert!(u.accumulated_impulse >= 0.0);
    }

    // ---- Jacobian assembly --------------------------------------------------

    /// 18. ConstraintJacobian mul_velocity gives Jv correctly.
    #[test]
    fn test_jacobian_mul_velocity() {
        let mut jac = ConstraintJacobian::new(1, 3);
        jac.set(0, 0, 1.0);
        jac.set(0, 1, 2.0);
        jac.set(0, 2, 3.0);
        let v = vec![1.0, 1.0, 1.0];
        let jv = jac.mul_velocity(&v);
        assert!((jv[0] - 6.0).abs() < EPS, "Jv={}", jv[0]);
    }

    /// 19. ConstraintJacobian mul_lambda gives J^T λ correctly.
    #[test]
    fn test_jacobian_mul_lambda() {
        let mut jac = ConstraintJacobian::new(1, 3);
        jac.set(0, 0, 2.0);
        jac.set(0, 1, 3.0);
        jac.set(0, 2, 4.0);
        let lambda = vec![1.0];
        let jt_lam = jac.mul_lambda(&lambda);
        assert!((jt_lam[0] - 2.0).abs() < EPS);
        assert!((jt_lam[1] - 3.0).abs() < EPS);
        assert!((jt_lam[2] - 4.0).abs() < EPS);
    }

    /// 20. Distance row in Jacobian has correct J_a and J_b signs.
    #[test]
    fn test_jacobian_distance_row() {
        let mut jac = ConstraintJacobian::new(1, 6);
        jac.add_distance_row(0, 0, 3, [1.0, 0.0, 0.0]);
        // J_a = -[1,0,0] at dofs 0..2, J_b = +[1,0,0] at dofs 3..5
        assert!((jac.get(0, 0) + 1.0).abs() < EPS);
        assert!((jac.get(0, 3) - 1.0).abs() < EPS);
    }

    // ---- Effective mass matrix ----------------------------------------------

    /// 21. Effective mass matrix with identity Jacobian and unit masses.
    #[test]
    fn test_effective_mass_matrix_unit() {
        let mut jac = ConstraintJacobian::new(1, 3);
        jac.set(0, 0, 1.0);
        let inv_mass = vec![1.0, 1.0, 1.0];
        let emm = jac.effective_mass_matrix(&inv_mass);
        assert!((emm[0][0] - 1.0).abs() < EPS, "emm={}", emm[0][0]);
    }

    // ---- DAE step -----------------------------------------------------------

    /// 22. DAE index-1 step conserves free particles (no constraint forces for free DOFs).
    #[test]
    fn test_dae_step_free_particle() {
        let mut state = DaeState::new(1, 0);
        state.v[0] = 2.0; // initial velocity
        let forces = vec![0.0];
        let jac = ConstraintJacobian::new(0, 1);
        let inv_mass = vec![1.0];
        let errors: Vec<f64> = vec![];
        let bp = BaumgarteParams::critically_damped(10.0);
        dae_step_index1(&mut state, &forces, &jac, &inv_mass, &errors, bp, 0.01);
        // Position should have advanced by dt * v = 0.01 * 2 = 0.02
        assert!((state.q[0] - 0.02).abs() < EPS, "q={}", state.q[0]);
    }

    /// 23. BDF-1 step moves a free particle correctly.
    #[test]
    fn test_dae_bdf1_free_particle() {
        let mut state = DaeState::new(1, 0);
        state.v[0] = 1.0;
        let forces = vec![0.0];
        let jac = ConstraintJacobian::new(0, 1);
        let inv_mass = vec![1.0];
        let errors: Vec<f64> = vec![];
        let bp = BaumgarteParams::critically_damped(10.0);
        dae_bdf1_step(&mut state, &forces, &jac, &inv_mass, &errors, bp, 0.1);
        assert!((state.q[0] - 0.1).abs() < EPS, "q={}", state.q[0]);
    }

    // ---- Penalty method -----------------------------------------------------

    /// 24. Penalty force is zero at zero error and zero velocity error.
    #[test]
    fn test_penalty_zero_error() {
        let p = PenaltyConstraint::new(1000.0, 10.0);
        let f = p.force(0.0, 0.0);
        assert!(f.abs() < EPS, "f={f}");
    }

    /// 25. Penalty force magnitude scales with stiffness.
    #[test]
    fn test_penalty_force_scales_with_k() {
        let p1 = PenaltyConstraint::new(100.0, 0.0);
        let p2 = PenaltyConstraint::new(200.0, 0.0);
        let f1 = p1.force(0.1, 0.0).abs();
        let f2 = p2.force(0.1, 0.0).abs();
        assert!((f2 - 2.0 * f1).abs() < EPS, "f1={f1} f2={f2}");
    }

    /// 26. Penalty energy is zero at rest.
    #[test]
    fn test_penalty_energy_zero() {
        assert!(penalty_energy(1000.0, 0.0).abs() < EPS);
    }

    /// 27. Penalty energy is positive and grows quadratically.
    #[test]
    fn test_penalty_energy_quadratic() {
        let e1 = penalty_energy(1.0, 1.0);
        let e2 = penalty_energy(1.0, 2.0);
        assert!((e2 - 4.0 * e1).abs() < EPS, "e1={e1} e2={e2}");
    }

    /// 28. Penalty stability limit is finite for positive stiffness.
    #[test]
    fn test_penalty_stability_limit_finite() {
        let lim = penalty_stability_limit(1000.0, 1.0);
        assert!(lim.is_finite() && lim > 0.0, "lim={lim}");
    }

    // ---- Lagrange multiplier solver -----------------------------------------

    /// 29. Gauss-Seidel solver converges on a simple 1×1 system.
    #[test]
    fn test_lagrange_gs_1x1() {
        let mut solver = LagrangeConstraintSolver::new(100, 1e-12);
        let jmjt = vec![vec![2.0]];
        let rhs = vec![4.0];
        let x = solver.solve_gs(&jmjt, &rhs);
        assert!((x[0] - 2.0).abs() < 1e-6, "x={}", x[0]);
    }

    /// 30. Gauss-Seidel solver converges on a diagonal 2×2 system.
    #[test]
    fn test_lagrange_gs_2x2_diagonal() {
        let mut solver = LagrangeConstraintSolver::new(200, 1e-12);
        let jmjt = vec![vec![3.0, 0.0], vec![0.0, 5.0]];
        let rhs = vec![6.0, 10.0];
        let x = solver.solve_gs(&jmjt, &rhs);
        assert!((x[0] - 2.0).abs() < 1e-6, "x[0]={}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-6, "x[1]={}", x[1]);
    }

    // ---- Rodrigues rotation -------------------------------------------------

    /// 31. Rodrigues rotation by 0 is identity.
    #[test]
    fn test_rodrigues_zero_rotation() {
        let v = [1.0, 2.0, 3.0];
        let r = rodrigues_rotate(v, [0.0, 0.0, 1.0], 0.0);
        for d in 0..3 {
            assert!((r[d] - v[d]).abs() < EPS, "d={d}");
        }
    }

    /// 32. Rodrigues rotation by 90° around Z maps \[1,0,0\] to \[0,1,0\].
    #[test]
    fn test_rodrigues_90_around_z() {
        let v = [1.0, 0.0, 0.0];
        let r = rodrigues_rotate(v, [0.0, 0.0, 1.0], PI / 2.0);
        assert!((r[0]).abs() < 1e-6, "x={}", r[0]);
        assert!((r[1] - 1.0).abs() < 1e-6, "y={}", r[1]);
    }

    /// 33. rotation_matrix by 0 is identity.
    #[test]
    fn test_rotation_matrix_identity() {
        let rm = rotation_matrix([0.0, 0.0, 1.0], 0.0);
        for (i, row) in rm.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < EPS, "rm[{i}][{j}]={val}",);
            }
        }
    }

    // ---- Acceleration-level residual ----------------------------------------

    /// 34. Acceleration residual is zero for a constraint satisfied at all levels.
    #[test]
    fn test_acceleration_residual_zero() {
        // Both bodies at rest, no acceleration, exact rest distance
        let ar = acceleration_residual_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
        );
        assert!(ar.abs() < EPS, "ar={ar}");
    }

    // ---- Hinge constraint ---------------------------------------------------

    /// 35. Hinge positional error is zero when anchor points coincide.
    #[test]
    fn test_hinge_positional_zero() {
        let h = HingeConstraintHolo::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let err = h.positional_error([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(norm3(err) < EPS);
    }

    /// 36. lagrange_impulse_step clamps to lo bound.
    #[test]
    fn test_lagrange_impulse_step_clamp_lo() {
        let mut acc = 0.0;
        // cdot = 10, bias = 0, eff_mass_inv = 1 → λ = -10
        let delta = lagrange_impulse_step(10.0, 0.0, 1.0, &mut acc, 0.0, f64::MAX);
        // acc clamped to 0 (lo), delta = 0 - 0 = 0
        assert!(delta.abs() < EPS, "delta={delta}");
        assert!(acc.abs() < EPS, "acc={acc}");
    }

    /// 37. Bilateral warm-start reset zeros accumulated impulse.
    #[test]
    fn test_bilateral_warm_start_reset() {
        let mut b = BilateralConstraint::new(0.0, 1.0);
        b.accumulated_impulse = 42.0;
        b.reset_warm_start();
        assert!(b.accumulated_impulse.abs() < EPS);
    }

    /// 38. Effective mass matrix is symmetric for a symmetric Jacobian.
    #[test]
    fn test_effective_mass_symmetry() {
        let mut jac = ConstraintJacobian::new(2, 4);
        jac.set(0, 0, 1.0);
        jac.set(0, 1, 0.5);
        jac.set(1, 2, 1.0);
        jac.set(1, 3, 0.5);
        let inv_mass = vec![1.0; 4];
        let emm = jac.effective_mass_matrix(&inv_mass);
        // For a block-diagonal Jacobian the cross terms should be zero
        assert!(emm[0][1].abs() < EPS, "emm[0][1]={}", emm[0][1]);
        assert!(emm[1][0].abs() < EPS, "emm[1][0]={}", emm[1][0]);
    }

    /// 39. Penalty critical damping coefficient is twice the natural frequency.
    #[test]
    fn test_penalty_critical_damping() {
        let p = PenaltyConstraint::new(100.0, 0.0);
        // eff_mass_inv = 1.0 → omega_n = sqrt(100) = 10 → d_crit = 20
        let d_crit = p.critical_damping_coefficient(1.0);
        assert!((d_crit - 20.0).abs() < EPS, "d_crit={d_crit}");
    }

    /// 40. Null-space multi projection satisfies all constraints.
    #[test]
    fn test_null_space_multi_projection() {
        let jacobians = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
        let eff_mass_inv = vec![1.0, 1.0];
        let mut v = vec![3.0, 4.0, 5.0];
        null_space_project_multi(&mut v, &jacobians, &eff_mass_inv, 10);
        // After projection J_0 v = 0 and J_1 v = 0
        assert!(
            dot(&jacobians[0], &v).abs() < 1e-8,
            "J0v={}",
            dot(&jacobians[0], &v)
        );
        assert!(
            dot(&jacobians[1], &v).abs() < 1e-8,
            "J1v={}",
            dot(&jacobians[1], &v)
        );
    }
}
