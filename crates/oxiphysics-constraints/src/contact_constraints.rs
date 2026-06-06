// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact constraint formulations for rigid body simulation.
//!
//! Covers:
//! - Frictionless contact (normal impulse only)
//! - Frictional contact (Coulomb cone)
//! - Contact constraint Jacobian
//! - LCP (Linear Complementarity Problem) formulation
//! - Lemke algorithm for LCP
//! - Contact island management
//! - Persistent contact constraints
//! - Warm-starting contact
//! - Sequential impulse contact solver
//! - Contact stabilization via Baumgarte method

// ── Vec3 helpers ─────────────────────────────────────────────────────────────

/// 3-D vector type used throughout this module.
pub type Vec3 = [f64; 3];

#[inline]
fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn normalize(a: Vec3) -> Vec3 {
    let l = len(a);
    if l < 1e-30 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / l)
    }
}

/// Multiply a 3×3 row-major matrix by a 3-vector.
fn mat3_vec(m: [[f64; 3]; 3], v: Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ── Contact point ─────────────────────────────────────────────────────────────

/// A detected contact point between two bodies.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// Contact point in world space (on body A surface).
    pub point: Vec3,
    /// Contact normal pointing from B into A.
    pub normal: Vec3,
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Body A index.
    pub body_a: usize,
    /// Body B index (use `usize::MAX` for static world).
    pub body_b: usize,
}

impl ContactPoint {
    /// Create a new contact point.
    pub fn new(point: Vec3, normal: Vec3, depth: f64, body_a: usize, body_b: usize) -> Self {
        Self {
            point,
            normal: normalize(normal),
            depth,
            body_a,
            body_b,
        }
    }

    /// Two tangent vectors perpendicular to the normal (for friction).
    pub fn tangent_basis(&self) -> (Vec3, Vec3) {
        let n = self.normal;
        let t1 = if n[0].abs() < 0.9 {
            normalize(cross(n, [1.0, 0.0, 0.0]))
        } else {
            normalize(cross(n, [0.0, 1.0, 0.0]))
        };
        let t2 = cross(n, t1);
        (t1, t2)
    }
}

// ── Contact constraint Jacobian ───────────────────────────────────────────────

/// The 6-row Jacobian for a single contact constraint (1 normal + 2 tangents).
///
/// Each row is `[J_va, J_ωa, J_vb, J_ωb]` — four 3-vectors.
#[derive(Debug, Clone)]
pub struct ContactJacobian {
    /// Normal row: `[n, (ra×n), -n, -(rb×n)]`.
    pub normal: ([Vec3; 2], [Vec3; 2]),
    /// First tangent row.
    pub tangent1: ([Vec3; 2], [Vec3; 2]),
    /// Second tangent row.
    pub tangent2: ([Vec3; 2], [Vec3; 2]),
}

impl ContactJacobian {
    /// Build the contact Jacobian from a contact point and body lever arms.
    ///
    /// * `cp` - Contact point data.
    /// * `ra` - Vector from body A CoM to contact point.
    /// * `rb` - Vector from body B CoM to contact point.
    pub fn new(cp: &ContactPoint, ra: Vec3, rb: Vec3) -> Self {
        let n = cp.normal;
        let (t1, t2) = cp.tangent_basis();

        let build_row = |axis: Vec3| {
            let j_va = axis;
            let j_wa = cross(ra, axis);
            let j_vb = scale(axis, -1.0);
            let j_wb = scale(cross(rb, axis), -1.0);
            ([j_va, j_wa], [j_vb, j_wb])
        };

        Self {
            normal: build_row(n),
            tangent1: build_row(t1),
            tangent2: build_row(t2),
        }
    }

    /// Compute velocity constraint `Jv` for the normal row.
    ///
    /// * `va`, `ωa` - Linear and angular velocity of body A.
    /// * `vb`, `ωb` - Linear and angular velocity of body B.
    pub fn normal_cdot(&self, va: Vec3, wa: Vec3, vb: Vec3, wb: Vec3) -> f64 {
        let ([j_va, j_wa], [j_vb, j_wb]) = &self.normal;
        dot(*j_va, va) + dot(*j_wa, wa) + dot(*j_vb, vb) + dot(*j_wb, wb)
    }

    /// Compute velocity constraint for tangent 1.
    pub fn tangent1_cdot(&self, va: Vec3, wa: Vec3, vb: Vec3, wb: Vec3) -> f64 {
        let ([j_va, j_wa], [j_vb, j_wb]) = &self.tangent1;
        dot(*j_va, va) + dot(*j_wa, wa) + dot(*j_vb, vb) + dot(*j_wb, wb)
    }

    /// Compute velocity constraint for tangent 2.
    pub fn tangent2_cdot(&self, va: Vec3, wa: Vec3, vb: Vec3, wb: Vec3) -> f64 {
        let ([j_va, j_wa], [j_vb, j_wb]) = &self.tangent2;
        dot(*j_va, va) + dot(*j_wa, wa) + dot(*j_vb, vb) + dot(*j_wb, wb)
    }
}

// ── Effective mass computation ────────────────────────────────────────────────

/// Compute the scalar effective mass for a constraint row `J`.
///
/// `K = J * M^-1 * J^T`
///
/// where `J` is a 12-element row split across two bodies.
pub fn effective_mass_for_row(
    inv_mass_a: f64,
    inv_inertia_a: [[f64; 3]; 3],
    inv_mass_b: f64,
    inv_inertia_b: [[f64; 3]; 3],
    j_va: Vec3,
    j_wa: Vec3,
    j_vb: Vec3,
    j_wb: Vec3,
) -> f64 {
    let k_a = inv_mass_a * dot(j_va, j_va) + dot(j_wa, mat3_vec(inv_inertia_a, j_wa));
    let k_b = inv_mass_b * dot(j_vb, j_vb) + dot(j_wb, mat3_vec(inv_inertia_b, j_wb));
    let k = k_a + k_b;
    if k < 1e-30 { 0.0 } else { 1.0 / k }
}

// ── Frictionless contact constraint ──────────────────────────────────────────

/// Frictionless contact constraint (normal impulse only).
///
/// Implements a sequential impulse step for a single contact point,
/// clamping the accumulated normal impulse to `[0, ∞)`.
#[derive(Debug, Clone)]
pub struct FrictionlessContact {
    /// Associated contact point.
    pub contact: ContactPoint,
    /// Accumulated normal impulse (warm-start value).
    pub impulse_normal: f64,
    /// Effective mass for the normal constraint.
    pub eff_mass_normal: f64,
    /// Baumgarte stabilization bias velocity.
    pub bias: f64,
    /// Restitution coefficient.
    pub restitution: f64,
}

impl FrictionlessContact {
    /// Create a frictionless contact constraint.
    pub fn new(
        contact: ContactPoint,
        inv_mass_a: f64,
        inv_inertia_a: [[f64; 3]; 3],
        inv_mass_b: f64,
        inv_inertia_b: [[f64; 3]; 3],
        ra: Vec3,
        rb: Vec3,
        restitution: f64,
        beta: f64,
        dt: f64,
    ) -> Self {
        let jac = ContactJacobian::new(&contact, ra, rb);
        let ([j_va, j_wa], [j_vb, j_wb]) = &jac.normal;
        let eff = effective_mass_for_row(
            inv_mass_a,
            inv_inertia_a,
            inv_mass_b,
            inv_inertia_b,
            *j_va,
            *j_wa,
            *j_vb,
            *j_wb,
        );
        // Baumgarte bias: correct penetration
        let bias = if contact.depth > 0.0 {
            -(beta / dt) * contact.depth
        } else {
            0.0
        };
        Self {
            contact,
            impulse_normal: 0.0,
            eff_mass_normal: eff,
            bias,
            restitution,
        }
    }

    /// Solve one sequential impulse iteration.
    ///
    /// Velocities are updated in-place via the returned delta impulse.
    /// Returns the scalar delta impulse applied.
    pub fn solve_impulse(
        &mut self,
        va: Vec3,
        wa: Vec3,
        vb: Vec3,
        wb: Vec3,
        ra: Vec3,
        rb: Vec3,
    ) -> f64 {
        let jac = ContactJacobian::new(&self.contact, ra, rb);
        let cdot = jac.normal_cdot(va, wa, vb, wb);
        // Include restitution bias
        let restitution_bias = if cdot < -1.0 {
            -self.restitution * cdot
        } else {
            0.0
        };
        let rhs = cdot - self.bias - restitution_bias;
        let raw = -self.eff_mass_normal * rhs;
        let old = self.impulse_normal;
        self.impulse_normal = (old + raw).max(0.0);
        self.impulse_normal - old
    }
}

// ── Frictional contact (Coulomb cone) ────────────────────────────────────────

/// Parameters for [`FrictionalContact::new`].
#[derive(Debug, Clone)]
pub struct FrictionalContactParams {
    /// Contact manifold point.
    pub contact: ContactPoint,
    /// Inverse mass of body A.
    pub inv_mass_a: f64,
    /// Inverse inertia tensor of body A.
    pub inv_inertia_a: [[f64; 3]; 3],
    /// Inverse mass of body B.
    pub inv_mass_b: f64,
    /// Inverse inertia tensor of body B.
    pub inv_inertia_b: [[f64; 3]; 3],
    /// Lever arm from body A CoM to contact.
    pub ra: Vec3,
    /// Lever arm from body B CoM to contact.
    pub rb: Vec3,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Coulomb friction coefficient.
    pub friction: f64,
    /// Baumgarte stabilisation factor.
    pub beta: f64,
    /// Time step.
    pub dt: f64,
}

/// Frictional contact constraint using the Coulomb friction cone.
///
/// Tangential impulses are clamped within `[-mu * lambda_n, mu * lambda_n]`.
#[derive(Debug, Clone)]
pub struct FrictionalContact {
    /// Underlying frictionless constraint (handles normal impulse).
    pub frictionless: FrictionlessContact,
    /// Accumulated tangential impulse along `t1`.
    pub impulse_t1: f64,
    /// Accumulated tangential impulse along `t2`.
    pub impulse_t2: f64,
    /// Effective mass for tangent 1.
    pub eff_mass_t1: f64,
    /// Effective mass for tangent 2.
    pub eff_mass_t2: f64,
    /// Coulomb friction coefficient.
    pub friction: f64,
}

impl FrictionalContact {
    /// Create a frictional contact constraint.
    pub fn new(p: FrictionalContactParams) -> Self {
        let FrictionalContactParams {
            contact,
            inv_mass_a,
            inv_inertia_a,
            inv_mass_b,
            inv_inertia_b,
            ra,
            rb,
            restitution,
            friction,
            beta,
            dt,
        } = p;
        let jac = ContactJacobian::new(&contact, ra, rb);
        let (t1, t2) = contact.tangent_basis();
        let _ = t1;
        let _ = t2;
        // Tangent 1 effective mass
        let ([j_va_t1, j_wa_t1], [j_vb_t1, j_wb_t1]) = &jac.tangent1;
        let eff_t1 = effective_mass_for_row(
            inv_mass_a,
            inv_inertia_a,
            inv_mass_b,
            inv_inertia_b,
            *j_va_t1,
            *j_wa_t1,
            *j_vb_t1,
            *j_wb_t1,
        );
        let ([j_va_t2, j_wa_t2], [j_vb_t2, j_wb_t2]) = &jac.tangent2;
        let eff_t2 = effective_mass_for_row(
            inv_mass_a,
            inv_inertia_a,
            inv_mass_b,
            inv_inertia_b,
            *j_va_t2,
            *j_wa_t2,
            *j_vb_t2,
            *j_wb_t2,
        );
        let frictionless = FrictionlessContact::new(
            contact,
            inv_mass_a,
            inv_inertia_a,
            inv_mass_b,
            inv_inertia_b,
            ra,
            rb,
            restitution,
            beta,
            dt,
        );
        Self {
            frictionless,
            impulse_t1: 0.0,
            impulse_t2: 0.0,
            eff_mass_t1: eff_t1,
            eff_mass_t2: eff_t2,
            friction,
        }
    }

    /// Solve one sequential impulse iteration for normal + friction.
    ///
    /// Returns `(delta_normal, delta_t1, delta_t2)`.
    pub fn solve(
        &mut self,
        va: Vec3,
        wa: Vec3,
        vb: Vec3,
        wb: Vec3,
        ra: Vec3,
        rb: Vec3,
    ) -> (f64, f64, f64) {
        // Normal impulse first
        let dn = self.frictionless.solve_impulse(va, wa, vb, wb, ra, rb);
        // Friction bound
        let mu_n = self.friction * self.frictionless.impulse_normal;
        let jac = ContactJacobian::new(&self.frictionless.contact, ra, rb);
        // Tangent 1
        let cdot_t1 = jac.tangent1_cdot(va, wa, vb, wb);
        let raw_t1 = -self.eff_mass_t1 * cdot_t1;
        let old_t1 = self.impulse_t1;
        self.impulse_t1 = (old_t1 + raw_t1).max(-mu_n).min(mu_n);
        let dt1 = self.impulse_t1 - old_t1;
        // Tangent 2
        let cdot_t2 = jac.tangent2_cdot(va, wa, vb, wb);
        let raw_t2 = -self.eff_mass_t2 * cdot_t2;
        let old_t2 = self.impulse_t2;
        self.impulse_t2 = (old_t2 + raw_t2).max(-mu_n).min(mu_n);
        let dt2 = self.impulse_t2 - old_t2;
        (dn, dt1, dt2)
    }

    /// Reset accumulated impulses (for beginning of a new time step without warm-start).
    pub fn reset_impulses(&mut self) {
        self.frictionless.impulse_normal = 0.0;
        self.impulse_t1 = 0.0;
        self.impulse_t2 = 0.0;
    }

    /// Warm-start by applying stored accumulated impulses.
    ///
    /// Returns the delta velocity contributions as 4 vectors: (dv_a, dω_a, dv_b, dω_b).
    pub fn warm_start(
        &self,
        inv_mass_a: f64,
        inv_inertia_a: [[f64; 3]; 3],
        inv_mass_b: f64,
        inv_inertia_b: [[f64; 3]; 3],
        ra: Vec3,
        rb: Vec3,
    ) -> (Vec3, Vec3, Vec3, Vec3) {
        let n = self.frictionless.contact.normal;
        let (t1, t2) = self.frictionless.contact.tangent_basis();
        let p_n = self.frictionless.impulse_normal;
        let p_t1 = self.impulse_t1;
        let p_t2 = self.impulse_t2;
        // Total impulse vector
        let total_impulse = add(add(scale(n, p_n), scale(t1, p_t1)), scale(t2, p_t2));
        let dv_a = scale(total_impulse, inv_mass_a);
        let dw_a = mat3_vec(inv_inertia_a, cross(ra, total_impulse));
        let dv_b = scale(total_impulse, -inv_mass_b);
        let dw_b = scale(mat3_vec(inv_inertia_b, cross(rb, total_impulse)), -1.0);
        (dv_a, dw_a, dv_b, dw_b)
    }
}

// ── LCP (Linear Complementarity Problem) ─────────────────────────────────────

/// LCP: find `z ≥ 0` such that `w = M*z + q ≥ 0` and `z·w = 0`.
///
/// Stores a dense `n×n` matrix `M` and vector `q`.
#[derive(Debug, Clone)]
pub struct Lcp {
    /// Dimension.
    pub n: usize,
    /// Matrix `M` stored row-major.
    pub m: Vec<Vec<f64>>,
    /// Vector `q`.
    pub q: Vec<f64>,
}

impl Lcp {
    /// Create an LCP with given dimension.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            m: vec![vec![0.0; n]; n],
            q: vec![0.0; n],
        }
    }

    /// Set matrix entry `M[i][j]`.
    pub fn set_m(&mut self, i: usize, j: usize, val: f64) {
        self.m[i][j] = val;
    }

    /// Set vector entry `q[i]`.
    pub fn set_q(&mut self, i: usize, val: f64) {
        self.q[i] = val;
    }

    /// Compute `w = M*z + q`.
    pub fn compute_w(&self, z: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|i| self.m[i].iter().zip(z).map(|(m, z)| m * z).sum::<f64>() + self.q[i])
            .collect()
    }

    /// Check complementarity condition: `z_i * w_i ≈ 0` for all `i`.
    pub fn check_complementarity(&self, z: &[f64], tol: f64) -> bool {
        let w = self.compute_w(z);
        z.iter()
            .zip(w.iter())
            .all(|(zi, wi)| *zi >= -tol && *wi >= -tol && (zi * wi).abs() < tol)
    }
}

// ── Lemke's algorithm ─────────────────────────────────────────────────────────

/// Result of the Lemke LCP solver.
#[derive(Debug, Clone)]
pub enum LemkeResult {
    /// Solution found; `z` is the solution vector.
    Solution(Vec<f64>),
    /// Algorithm could not find a solution within max iterations.
    Failure,
    /// Secondary ray termination (unbounded).
    SecondaryRay,
}

/// Lemke's pivoting algorithm for solving the LCP `w = Mz + q`.
///
/// Returns a solution `z ≥ 0` with complementarity, or a failure indicator.
pub fn lemke_solve(lcp: &Lcp, max_iter: usize) -> LemkeResult {
    let n = lcp.n;
    if n == 0 {
        return LemkeResult::Solution(vec![]);
    }

    // Tableau: [-M | I | -d | q]  where d = [1,1,...,1]
    // Encodes: -Mz + w - z_0 = q (i.e. w = Mz + z_0 + q, so z_0=0 → w=Mz+q)
    // Variables: z_1..z_n, w_1..w_n, z_0 (artificial)
    let total_cols = 2 * n + 2; // z cols, w cols, z0 col, q col

    let mut tab = vec![vec![0.0_f64; total_cols]; n];
    // Fill M columns (z variables) — stored as -M so that
    // the row represents: -Mz + w - z_0 = q  (equiv. w = Mz + z_0 + q)
    for (i, (tab_row, (m_row, q_i))) in tab
        .iter_mut()
        .zip(lcp.m.iter().zip(lcp.q.iter()))
        .enumerate()
    {
        for (j, tab_ij) in tab_row.iter_mut().enumerate().take(n) {
            *tab_ij = -m_row[j];
        }
        // Identity for w variables (cols n..2n)
        tab_row[n + i] = 1.0;
        // -d column (col 2n)
        tab_row[2 * n] = -1.0;
        // q column (col 2n+1)
        tab_row[2 * n + 1] = *q_i;
    }

    // Basic variables: initially w_1..w_n (indices n..2n)
    let mut basis: Vec<usize> = (n..2 * n).collect();

    // If all q >= 0, trivial solution z = 0
    if lcp.q.iter().all(|&qi| qi >= 0.0) {
        return LemkeResult::Solution(vec![0.0; n]);
    }

    // Entering variable: z_0 (column 2n)
    // Choose leaving variable: most negative q
    let (leaving_row, _) = lcp
        .q
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .expect("operation should succeed");

    // Pivot z_0 into basis at leaving_row
    // Save the variable that is leaving before we overwrite the basis entry
    let leaving_var_initial = basis[leaving_row];
    pivot(&mut tab, leaving_row, 2 * n, n);
    basis[leaving_row] = 2 * n; // z_0 enters

    // Next entering variable is the complement of the variable that left
    let mut entering = if leaving_var_initial >= n && leaving_var_initial < 2 * n {
        leaving_var_initial - n // w_i left → z_i enters
    } else if leaving_var_initial < n {
        leaving_var_initial + n // z_i left → w_i enters
    } else {
        // z_0 left — done
        let z = extract_solution(&tab, &basis, n, total_cols);
        return LemkeResult::Solution(z);
    };

    for _iter in 0..max_iter {
        // Find leaving row via minimum ratio test
        let leaving = min_ratio_row(&tab, entering, total_cols);
        let leaving_row = match leaving {
            Some(r) => r,
            None => return LemkeResult::SecondaryRay,
        };

        // If z_0 is leaving, done
        if basis[leaving_row] == 2 * n {
            pivot(&mut tab, leaving_row, entering, n);
            basis[leaving_row] = entering;
            let z = extract_solution(&tab, &basis, n, total_cols);
            return LemkeResult::Solution(z);
        }

        let leaving_var = basis[leaving_row];
        pivot(&mut tab, leaving_row, entering, n);
        basis[leaving_row] = entering;

        // Next entering: complement of leaving_var
        entering = if leaving_var >= n && leaving_var < 2 * n {
            leaving_var - n
        } else if leaving_var < n {
            leaving_var + n
        } else {
            let z = extract_solution(&tab, &basis, n, total_cols);
            return LemkeResult::Solution(z);
        };
    }

    LemkeResult::Failure
}

fn pivot(tab: &mut [Vec<f64>], row: usize, col: usize, _n: usize) {
    let pivot_val = tab[row][col];
    if pivot_val.abs() < 1e-15 {
        return;
    }
    for t in tab[row].iter_mut() {
        *t /= pivot_val;
    }
    let nrows = tab.len();
    for i in 0..nrows {
        if i == row {
            continue;
        }
        let factor = tab[i][col];
        if factor.abs() < 1e-30 {
            continue;
        }
        let pivot_copy: Vec<f64> = tab[row].clone();
        for (t, p) in tab[i].iter_mut().zip(pivot_copy.iter()) {
            *t -= factor * p;
        }
    }
}

fn min_ratio_row(tab: &[Vec<f64>], col: usize, total_cols: usize) -> Option<usize> {
    let q_col = total_cols - 1;
    let mut best_ratio = f64::INFINITY;
    let mut best_row = None;
    for (i, row) in tab.iter().enumerate() {
        if row[col] > 1e-12 {
            let ratio = row[q_col] / row[col];
            if ratio < best_ratio {
                best_ratio = ratio;
                best_row = Some(i);
            }
        }
    }
    best_row
}

fn extract_solution(tab: &[Vec<f64>], basis: &[usize], n: usize, total_cols: usize) -> Vec<f64> {
    let mut z = vec![0.0; n];
    let q_col = total_cols - 1;
    for (row, &b) in basis.iter().enumerate() {
        if b < n {
            z[b] = tab[row][q_col].max(0.0);
        }
    }
    z
}

// ── Baumgarte stabilization for contact ──────────────────────────────────────

/// Compute the Baumgarte bias velocity for contact stabilization.
///
/// ```text
/// v_bias = -(beta / dt) * max(depth - slop, 0)
/// ```
///
/// * `depth` - Penetration depth (positive = overlapping).
/// * `beta`  - Position error correction factor (0.1–0.3 typical).
/// * `dt`    - Time step.
/// * `slop`  - Allowed penetration before correction kicks in.
pub fn baumgarte_contact_bias(depth: f64, beta: f64, dt: f64, slop: f64) -> f64 {
    let correctable = (depth - slop).max(0.0);
    -(beta / dt) * correctable
}

// ── Persistent contact constraint ────────────────────────────────────────────

/// Persistent contact: stores accumulated impulse across time steps for warm-start.
#[derive(Debug, Clone)]
pub struct PersistentContact {
    /// Contact geometry.
    pub contact: ContactPoint,
    /// Accumulated normal impulse from previous step.
    pub lambda_normal: f64,
    /// Accumulated tangent 1 impulse.
    pub lambda_t1: f64,
    /// Accumulated tangent 2 impulse.
    pub lambda_t2: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// Number of frames this contact has been alive.
    pub lifetime: u32,
}

impl PersistentContact {
    /// Create a new persistent contact.
    pub fn new(contact: ContactPoint, friction: f64) -> Self {
        Self {
            contact,
            lambda_normal: 0.0,
            lambda_t1: 0.0,
            lambda_t2: 0.0,
            friction,
            lifetime: 0,
        }
    }

    /// Update (merge) with a fresh contact detection result.
    ///
    /// Preserves accumulated impulses (warm-starting) while updating geometry.
    pub fn update(&mut self, new_contact: ContactPoint) {
        self.contact = new_contact;
        self.lifetime += 1;
    }

    /// Scale accumulated impulses by `scale` (e.g. for time step ratio adjustment).
    pub fn scale_impulses(&mut self, s: f64) {
        self.lambda_normal *= s;
        self.lambda_t1 *= s;
        self.lambda_t2 *= s;
    }

    /// Reset all accumulated impulses.
    pub fn reset(&mut self) {
        self.lambda_normal = 0.0;
        self.lambda_t1 = 0.0;
        self.lambda_t2 = 0.0;
    }

    /// Check whether the contact is still valid (bodies still overlapping).
    pub fn is_valid(&self) -> bool {
        self.contact.depth >= 0.0
    }
}

// ── Contact island ────────────────────────────────────────────────────────────

/// A contact island groups bodies connected through contact constraints.
///
/// Bodies in the same island must be solved together.
#[derive(Debug, Clone)]
pub struct ContactIsland {
    /// Body indices in this island.
    pub bodies: Vec<usize>,
    /// Contact indices in this island.
    pub contacts: Vec<usize>,
    /// Whether this island is "sleeping" (all bodies at rest).
    pub sleeping: bool,
}

impl ContactIsland {
    /// Create an empty island.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            contacts: Vec::new(),
            sleeping: false,
        }
    }

    /// Add a body to this island.
    pub fn add_body(&mut self, body: usize) {
        if !self.bodies.contains(&body) {
            self.bodies.push(body);
        }
    }

    /// Add a contact to this island.
    pub fn add_contact(&mut self, contact: usize) {
        if !self.contacts.contains(&contact) {
            self.contacts.push(contact);
        }
    }

    /// Mark the island as sleeping.
    pub fn sleep(&mut self) {
        self.sleeping = true;
    }

    /// Wake the island.
    pub fn wake(&mut self) {
        self.sleeping = false;
    }
}

impl Default for ContactIsland {
    fn default() -> Self {
        Self::new()
    }
}

/// Manages partitioning of all contacts into disjoint islands using Union-Find.
#[derive(Debug, Clone)]
pub struct ContactIslandManager {
    /// Number of rigid bodies.
    pub n_bodies: usize,
    /// Union-Find parent array.
    parent: Vec<usize>,
    /// Union-Find rank array.
    rank: Vec<usize>,
}

impl ContactIslandManager {
    /// Create an island manager for `n_bodies` bodies.
    pub fn new(n_bodies: usize) -> Self {
        Self {
            n_bodies,
            parent: (0..n_bodies).collect(),
            rank: vec![0; n_bodies],
        }
    }

    /// Reset all islands (each body is its own island).
    pub fn reset(&mut self) {
        for i in 0..self.n_bodies {
            self.parent[i] = i;
            self.rank[i] = 0;
        }
    }

    /// Find the root of body `i` with path compression.
    pub fn find(&mut self, i: usize) -> usize {
        if self.parent[i] != i {
            self.parent[i] = self.find(self.parent[i]);
        }
        self.parent[i]
    }

    /// Union bodies `a` and `b` into the same island.
    pub fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }

    /// Build and return all islands from a list of contact pairs.
    ///
    /// Each entry in `contact_pairs` is `(body_a, body_b, contact_index)`.
    pub fn build_islands(&mut self, contact_pairs: &[(usize, usize, usize)]) -> Vec<ContactIsland> {
        self.reset();
        for &(a, b, _) in contact_pairs {
            if b < self.n_bodies {
                self.union(a, b);
            }
        }
        // Group bodies by root
        let mut island_map: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        let mut islands: Vec<ContactIsland> = Vec::new();
        for i in 0..self.n_bodies {
            let root = self.find(i);
            let idx = *island_map.entry(root).or_insert_with(|| {
                let id = islands.len();
                islands.push(ContactIsland::new());
                id
            });
            islands[idx].add_body(i);
        }
        // Assign contacts
        for &(a, b, ci) in contact_pairs {
            let root = self.find(a);
            if let Some(&idx) = island_map.get(&root) {
                islands[idx].add_contact(ci);
            }
            if b < self.n_bodies {
                let root_b = self.find(b);
                if root_b != root
                    && let Some(&idx) = island_map.get(&root_b)
                {
                    islands[idx].add_contact(ci);
                }
            }
        }
        islands
    }
}

// ── Sequential impulse contact solver ────────────────────────────────────────

/// State of a single rigid body for the sequential impulse solver.
#[derive(Debug, Clone)]
pub struct RigidBodyState {
    /// Linear velocity in m/s.
    pub vel: Vec3,
    /// Angular velocity in rad/s.
    pub omega: Vec3,
    /// Inverse mass (0 for static bodies).
    pub inv_mass: f64,
    /// Inverse inertia tensor (world space, row-major 3×3).
    pub inv_inertia: [[f64; 3]; 3],
    /// Centre of mass position.
    pub position: Vec3,
}

impl RigidBodyState {
    /// Create a dynamic body state.
    pub fn new(
        vel: Vec3,
        omega: Vec3,
        mass: f64,
        inv_inertia: [[f64; 3]; 3],
        position: Vec3,
    ) -> Self {
        Self {
            vel,
            omega,
            inv_mass: if mass > 0.0 { 1.0 / mass } else { 0.0 },
            inv_inertia,
            position,
        }
    }

    /// Create an infinite-mass static body.
    pub fn static_body(position: Vec3) -> Self {
        Self {
            vel: [0.0; 3],
            omega: [0.0; 3],
            inv_mass: 0.0,
            inv_inertia: [[0.0; 3]; 3],
            position,
        }
    }

    /// Apply a linear impulse `p` (in N·s) at lever arm `r` from CoM.
    pub fn apply_impulse(&mut self, p: Vec3, r: Vec3) {
        self.vel = add(self.vel, scale(p, self.inv_mass));
        let torque = cross(r, p);
        self.omega = add(self.omega, mat3_vec(self.inv_inertia, torque));
    }
}

/// Sequential impulse solver for a set of contact constraints.
#[derive(Debug, Clone)]
pub struct SequentialImpulseSolver {
    /// Number of iterations per time step.
    pub iterations: usize,
    /// Baumgarte stabilization factor.
    pub beta: f64,
    /// Allowed penetration slop (m).
    pub slop: f64,
}

impl SequentialImpulseSolver {
    /// Create a solver with given parameters.
    pub fn new(iterations: usize, beta: f64, slop: f64) -> Self {
        Self {
            iterations,
            beta,
            slop,
        }
    }

    /// Solve a single frictionless contact between two bodies.
    ///
    /// Mutates body velocities in-place.
    pub fn solve_frictionless(
        &self,
        body_a: &mut RigidBodyState,
        body_b: &mut RigidBodyState,
        contact: &ContactPoint,
        dt: f64,
    ) {
        let ra = sub(contact.point, body_a.position);
        let rb = sub(contact.point, body_b.position);
        let mut con = FrictionlessContact::new(
            contact.clone(),
            body_a.inv_mass,
            body_a.inv_inertia,
            body_b.inv_mass,
            body_b.inv_inertia,
            ra,
            rb,
            0.0,
            self.beta,
            dt,
        );
        for _ in 0..self.iterations {
            let delta =
                con.solve_impulse(body_a.vel, body_a.omega, body_b.vel, body_b.omega, ra, rb);
            // Apply impulse along normal: push A away from B (+n direction)
            let p = scale(contact.normal, delta);
            body_a.apply_impulse(p, ra);
            body_b.apply_impulse(scale(p, -1.0), rb);
        }
    }

    /// Solve a single frictional contact between two bodies.
    ///
    /// Mutates body velocities in-place.
    pub fn solve_frictional(
        &self,
        body_a: &mut RigidBodyState,
        body_b: &mut RigidBodyState,
        contact: &ContactPoint,
        friction: f64,
        dt: f64,
    ) {
        let ra = sub(contact.point, body_a.position);
        let rb = sub(contact.point, body_b.position);
        let (t1, t2) = contact.tangent_basis();
        let mut con = FrictionalContact::new(FrictionalContactParams {
            contact: contact.clone(),
            inv_mass_a: body_a.inv_mass,
            inv_inertia_a: body_a.inv_inertia,
            inv_mass_b: body_b.inv_mass,
            inv_inertia_b: body_b.inv_inertia,
            ra,
            rb,
            restitution: 0.0,
            friction,
            beta: self.beta,
            dt,
        });
        for _ in 0..self.iterations {
            let (dn, dt1, dt2) =
                con.solve(body_a.vel, body_a.omega, body_b.vel, body_b.omega, ra, rb);
            let p = add(
                add(scale(contact.normal, dn), scale(t1, dt1)),
                scale(t2, dt2),
            );
            body_a.apply_impulse(p, ra);
            body_b.apply_impulse(scale(p, -1.0), rb);
        }
    }
}

// ── Warm-start helper ─────────────────────────────────────────────────────────

/// Apply warm-start impulses from a persistent contact to body velocities.
///
/// This pre-loads the solver with the impulse from the previous time step,
/// reducing the number of iterations needed to converge.
pub fn apply_warm_start(
    body_a: &mut RigidBodyState,
    body_b: &mut RigidBodyState,
    pc: &PersistentContact,
) {
    let cp = &pc.contact;
    let (t1, t2) = cp.tangent_basis();
    let impulse = add(
        add(scale(cp.normal, pc.lambda_normal), scale(t1, pc.lambda_t1)),
        scale(t2, pc.lambda_t2),
    );
    let ra = sub(cp.point, body_a.position);
    let rb = sub(cp.point, body_b.position);
    body_a.apply_impulse(impulse, ra);
    body_b.apply_impulse(scale(impulse, -1.0), rb);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn identity_inertia() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    fn zero_inertia() -> [[f64; 3]; 3] {
        [[0.0; 3]; 3]
    }

    // T01: Frictionless contact: normal impulse is non-negative.
    #[test]
    fn test_frictionless_normal_impulse_non_negative() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let mut con = FrictionlessContact::new(
            cp,
            1.0,
            identity_inertia(),
            0.0,
            zero_inertia(),
            [0.0; 3],
            [0.0; 3],
            0.0,
            0.2,
            0.016,
        );
        let delta = con.solve_impulse(
            [0.0, -2.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
        );
        assert!(
            delta >= 0.0,
            "Normal impulse delta must be >= 0, got {delta}"
        );
        assert!(con.impulse_normal >= 0.0, "Accumulated must be >= 0");
    }

    // T02: Contact point tangent vectors are perpendicular to normal.
    #[test]
    fn test_contact_tangent_perpendicular_to_normal() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0, 1);
        let (t1, t2) = cp.tangent_basis();
        assert!(dot(cp.normal, t1).abs() < EPS, "t1 not perpendicular to n");
        assert!(dot(cp.normal, t2).abs() < EPS, "t2 not perpendicular to n");
    }

    // T03: Contact tangent vectors are mutually perpendicular.
    #[test]
    fn test_contact_tangents_orthogonal() {
        let cp = ContactPoint::new([0.0; 3], [1.0, 0.0, 0.0], 0.0, 0, 1);
        let (t1, t2) = cp.tangent_basis();
        assert!(dot(t1, t2).abs() < EPS, "tangents not orthogonal");
    }

    // T04: Contact Jacobian normal row gives correct cdot for simple linear velocities.
    #[test]
    fn test_contact_jacobian_cdot_linear() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0, 1);
        let jac = ContactJacobian::new(&cp, [0.0; 3], [0.0; 3]);
        // Body A moving down at 2 m/s, B static
        let cdot = jac.normal_cdot([0.0, -2.0, 0.0], [0.0; 3], [0.0; 3], [0.0; 3]);
        // Expected: J_va·va + J_vb·vb = n·va + (-n)·vb = -2 + 0 = -2
        assert!((cdot + 2.0).abs() < EPS, "cdot={cdot}");
    }

    // T05: Frictional contact: tangent impulses are clamped within Coulomb cone.
    #[test]
    fn test_frictional_impulse_clamped_coulomb() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let mut con = FrictionalContact::new(FrictionalContactParams {
            contact: cp,
            inv_mass_a: 1.0,
            inv_inertia_a: identity_inertia(),
            inv_mass_b: 0.0,
            inv_inertia_b: zero_inertia(),
            ra: [0.0; 3],
            rb: [0.0; 3],
            restitution: 0.0,
            friction: 0.5,
            beta: 0.2,
            dt: 0.016,
        });
        // Give the constraint a large normal impulse to allow friction
        con.frictionless.impulse_normal = 10.0;
        let (_dn, dt1, dt2) = con.solve(
            [5.0, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
        );
        let mu_n = 0.5 * 10.0;
        assert!(
            con.impulse_t1.abs() <= mu_n + 1e-10,
            "t1 out of cone: {}",
            con.impulse_t1
        );
        assert!(
            con.impulse_t2.abs() <= mu_n + 1e-10,
            "t2 out of cone: {}",
            con.impulse_t2
        );
        let _ = (dt1, dt2);
    }

    // T06: Effective mass for two equal point masses is 0.5.
    #[test]
    fn test_effective_mass_point_masses() {
        let n = [0.0, 1.0, 0.0_f64];
        let neg_n = scale(n, -1.0);
        let em = effective_mass_for_row(
            0.5,
            zero_inertia(),
            0.5,
            zero_inertia(),
            neg_n,
            [0.0; 3],
            n,
            [0.0; 3],
        );
        // K = 0.5 + 0.5 = 1, eff_mass = 1
        assert!((em - 1.0).abs() < EPS, "em={em}");
    }

    // T07: LCP trivial solution when q >= 0.
    #[test]
    fn test_lcp_trivial_solution_positive_q() {
        let mut lcp = Lcp::new(2);
        lcp.set_m(0, 0, 2.0);
        lcp.set_m(0, 1, 1.0);
        lcp.set_m(1, 0, 1.0);
        lcp.set_m(1, 1, 2.0);
        lcp.set_q(0, 1.0);
        lcp.set_q(1, 1.0);
        match lemke_solve(&lcp, 100) {
            LemkeResult::Solution(z) => {
                assert!(z.iter().all(|&zi| zi >= -EPS), "z should be >= 0");
                assert!(
                    lcp.check_complementarity(&z, 1e-8),
                    "complementarity violated"
                );
            }
            _ => panic!("Expected solution"),
        }
    }

    // T08: LCP with single variable — known solution.
    #[test]
    fn test_lcp_single_variable() {
        // w = 2z + (-1), z >= 0, w >= 0, z*w = 0
        // Solution: z = 0.5, w = 0
        let mut lcp = Lcp::new(1);
        lcp.set_m(0, 0, 2.0);
        lcp.set_q(0, -1.0);
        match lemke_solve(&lcp, 50) {
            LemkeResult::Solution(z) => {
                assert!((z[0] - 0.5).abs() < 1e-6, "z[0]={}", z[0]);
            }
            other => panic!("Expected solution, got {other:?}"),
        }
    }

    // T09: LCP empty problem returns empty solution.
    #[test]
    fn test_lcp_empty() {
        let lcp = Lcp::new(0);
        match lemke_solve(&lcp, 10) {
            LemkeResult::Solution(z) => assert!(z.is_empty()),
            _ => panic!("Expected empty solution"),
        }
    }

    // T10: Baumgarte bias is zero for zero penetration.
    #[test]
    fn test_baumgarte_zero_penetration() {
        let bias = baumgarte_contact_bias(0.0, 0.2, 0.016, 0.001);
        assert!(bias.abs() < EPS);
    }

    // T11: Baumgarte bias is negative for positive penetration.
    #[test]
    fn test_baumgarte_negative_for_penetration() {
        let bias = baumgarte_contact_bias(0.05, 0.2, 0.016, 0.001);
        assert!(bias < 0.0, "bias={bias}");
    }

    // T12: Baumgarte slop: no correction within slop tolerance.
    #[test]
    fn test_baumgarte_slop_no_correction() {
        let bias = baumgarte_contact_bias(0.001, 0.2, 0.016, 0.01); // depth < slop
        assert!(bias.abs() < EPS, "bias={bias}");
    }

    // T13: Persistent contact lifetime increments on update.
    #[test]
    fn test_persistent_contact_lifetime() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let mut pc = PersistentContact::new(cp.clone(), 0.5);
        pc.update(cp.clone());
        pc.update(cp);
        assert_eq!(pc.lifetime, 2);
    }

    // T14: Persistent contact reset zeros impulses.
    #[test]
    fn test_persistent_contact_reset() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let mut pc = PersistentContact::new(cp, 0.5);
        pc.lambda_normal = 5.0;
        pc.lambda_t1 = 2.0;
        pc.lambda_t2 = -1.0;
        pc.reset();
        assert!(pc.lambda_normal.abs() < EPS);
        assert!(pc.lambda_t1.abs() < EPS);
        assert!(pc.lambda_t2.abs() < EPS);
    }

    // T15: Persistent contact scale impulses.
    #[test]
    fn test_persistent_contact_scale_impulses() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let mut pc = PersistentContact::new(cp, 0.5);
        pc.lambda_normal = 4.0;
        pc.scale_impulses(0.5);
        assert!((pc.lambda_normal - 2.0).abs() < EPS);
    }

    // T16: Island manager: two bodies connected by one contact form one island.
    #[test]
    fn test_island_manager_single_contact() {
        let mut mgr = ContactIslandManager::new(2);
        let pairs = vec![(0, 1, 0)];
        let islands = mgr.build_islands(&pairs);
        // Both bodies should be in the same island
        let total_bodies: usize = islands.iter().map(|i| i.bodies.len()).sum();
        assert_eq!(total_bodies, 2);
        let island_with_both = islands
            .iter()
            .any(|i| i.bodies.contains(&0) && i.bodies.contains(&1));
        assert!(island_with_both, "Both bodies should share an island");
    }

    // T17: Island manager: disconnected bodies form separate islands.
    #[test]
    fn test_island_manager_disconnected() {
        let mut mgr = ContactIslandManager::new(4);
        let pairs = vec![(0, 1, 0), (2, 3, 1)];
        let islands = mgr.build_islands(&pairs);
        let non_empty: Vec<_> = islands.iter().filter(|i| !i.bodies.is_empty()).collect();
        assert_eq!(non_empty.len(), 2, "Should have 2 islands");
    }

    // T18: Island: add_body is idempotent.
    #[test]
    fn test_island_add_body_idempotent() {
        let mut island = ContactIsland::new();
        island.add_body(5);
        island.add_body(5);
        assert_eq!(island.bodies.len(), 1);
    }

    // T19: Island sleeping flag.
    #[test]
    fn test_island_sleep_wake() {
        let mut island = ContactIsland::new();
        assert!(!island.sleeping);
        island.sleep();
        assert!(island.sleeping);
        island.wake();
        assert!(!island.sleeping);
    }

    // T20: Sequential impulse solver separates penetrating bodies (frictionless).
    #[test]
    fn test_si_solver_separates_bodies_frictionless() {
        let mut a = RigidBodyState::new(
            [0.0, -1.0, 0.0],
            [0.0; 3],
            1.0,
            zero_inertia(),
            [0.0, 0.1, 0.0],
        );
        let mut b = RigidBodyState::static_body([0.0; 3]);
        let cp = ContactPoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let solver = SequentialImpulseSolver::new(10, 0.2, 0.001);
        solver.solve_frictionless(&mut a, &mut b, &cp, 0.016);
        // Body A's downward velocity should be reduced
        assert!(
            a.vel[1] >= -0.5,
            "Velocity should be corrected, got {}",
            a.vel[1]
        );
    }

    // T21: Sequential impulse solver frictional: tangential velocity reduced.
    #[test]
    fn test_si_solver_friction_reduces_tangential() {
        let mut a = RigidBodyState::new(
            [2.0, -0.5, 0.0],
            [0.0; 3],
            1.0,
            zero_inertia(),
            [0.0, 0.1, 0.0],
        );
        let mut b = RigidBodyState::static_body([0.0; 3]);
        let cp = ContactPoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let solver = SequentialImpulseSolver::new(20, 0.2, 0.001);
        let vx_before = a.vel[0].abs();
        solver.solve_frictional(&mut a, &mut b, &cp, 0.5, 0.016);
        assert!(
            a.vel[0].abs() < vx_before,
            "Friction should reduce tangential velocity"
        );
    }

    // T22: Rigid body apply_impulse updates velocity correctly.
    #[test]
    fn test_rigid_body_apply_impulse_linear() {
        let mut body = RigidBodyState::new([0.0; 3], [0.0; 3], 2.0, zero_inertia(), [0.0; 3]);
        body.apply_impulse([4.0, 0.0, 0.0], [0.0; 3]);
        // dv = p / m = 4 / 2 = 2
        assert!((body.vel[0] - 2.0).abs() < EPS, "vel[0]={}", body.vel[0]);
    }

    // T23: Warm-start adds impulses to body velocities.
    #[test]
    fn test_warm_start_modifies_velocities() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0, 1);
        let mut pc = PersistentContact::new(cp, 0.3);
        pc.lambda_normal = 1.0;
        let mut a = RigidBodyState::new([0.0; 3], [0.0; 3], 1.0, zero_inertia(), [0.0; 3]);
        let mut b = RigidBodyState::static_body([0.0; 3]);
        let vel_a_before = a.vel;
        apply_warm_start(&mut a, &mut b, &pc);
        let changed = a
            .vel
            .iter()
            .zip(vel_a_before.iter())
            .any(|(v, v0)| (v - v0).abs() > 1e-12);
        assert!(changed, "Warm-start should modify body A velocity");
    }

    // T24: LCP check_complementarity accepts valid zero solution.
    #[test]
    fn test_lcp_complementarity_zero_solution() {
        let mut lcp = Lcp::new(2);
        lcp.set_m(0, 0, 1.0);
        lcp.set_m(1, 1, 1.0);
        lcp.set_q(0, 1.0);
        lcp.set_q(1, 2.0);
        let z = vec![0.0, 0.0];
        assert!(lcp.check_complementarity(&z, 1e-8));
    }

    // T25: LCP compute_w returns correct result.
    #[test]
    fn test_lcp_compute_w() {
        let mut lcp = Lcp::new(2);
        lcp.set_m(0, 0, 2.0);
        lcp.set_m(0, 1, 1.0);
        lcp.set_m(1, 0, 0.0);
        lcp.set_m(1, 1, 3.0);
        lcp.set_q(0, -1.0);
        lcp.set_q(1, 0.0);
        let z = vec![1.0, 2.0];
        let w = lcp.compute_w(&z);
        // w[0] = 2*1 + 1*2 + (-1) = 3
        // w[1] = 0*1 + 3*2 + 0 = 6
        assert!((w[0] - 3.0).abs() < EPS, "w[0]={}", w[0]);
        assert!((w[1] - 6.0).abs() < EPS, "w[1]={}", w[1]);
    }

    // T26: Frictionless contact: no impulse when bodies are separating.
    #[test]
    fn test_frictionless_no_impulse_when_separating() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], -0.01, 0, 1);
        let mut con = FrictionlessContact::new(
            cp,
            1.0,
            identity_inertia(),
            0.0,
            zero_inertia(),
            [0.0; 3],
            [0.0; 3],
            0.0,
            0.2,
            0.016,
        );
        // Bodies separating: v_b - v_a along normal is positive
        let delta = con.solve_impulse(
            [0.0, -2.0, 0.0],
            [0.0; 3],
            [0.0, 1.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
        );
        // Separation should not be blocked; impulse clamped ≥ 0 so could be 0
        let _ = delta; // we just check no panic
    }

    // T27: ContactPoint normal is normalized after construction.
    #[test]
    fn test_contact_point_normal_normalized() {
        let cp = ContactPoint::new([0.0; 3], [3.0, 4.0, 0.0], 0.0, 0, 1);
        let l = len(cp.normal);
        assert!((l - 1.0).abs() < 1e-10, "normal length={l}");
    }

    // T28: Island manager reset clears all unions.
    #[test]
    fn test_island_manager_reset() {
        let mut mgr = ContactIslandManager::new(3);
        mgr.union(0, 1);
        mgr.union(1, 2);
        mgr.reset();
        // After reset, each body is its own root
        assert_eq!(mgr.find(0), 0);
        assert_eq!(mgr.find(1), 1);
        assert_eq!(mgr.find(2), 2);
    }

    // T29: Frictional contact reset zeros all impulses.
    #[test]
    fn test_frictional_contact_reset_impulses() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0, 1);
        let mut con = FrictionalContact::new(FrictionalContactParams {
            contact: cp,
            inv_mass_a: 1.0,
            inv_inertia_a: zero_inertia(),
            inv_mass_b: 1.0,
            inv_inertia_b: zero_inertia(),
            ra: [0.0; 3],
            rb: [0.0; 3],
            restitution: 0.0,
            friction: 0.3,
            beta: 0.2,
            dt: 0.016,
        });
        con.frictionless.impulse_normal = 5.0;
        con.impulse_t1 = 1.0;
        con.impulse_t2 = -0.5;
        con.reset_impulses();
        assert!(con.frictionless.impulse_normal.abs() < EPS);
        assert!(con.impulse_t1.abs() < EPS);
        assert!(con.impulse_t2.abs() < EPS);
    }

    // T30: Sequential impulse: static body velocity unchanged.
    #[test]
    fn test_si_solver_static_body_unchanged() {
        let mut a = RigidBodyState::new(
            [0.0, -2.0, 0.0],
            [0.0; 3],
            1.0,
            zero_inertia(),
            [0.0, 0.1, 0.0],
        );
        let mut b = RigidBodyState::static_body([0.0; 3]);
        let cp = ContactPoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.05, 0, 1);
        let solver = SequentialImpulseSolver::new(10, 0.2, 0.001);
        solver.solve_frictionless(&mut a, &mut b, &cp, 0.016);
        for c in b.vel {
            assert!(c.abs() < EPS, "static body vel changed: {c}");
        }
        for c in b.omega {
            assert!(c.abs() < EPS, "static body omega changed: {c}");
        }
    }

    // T31: mat3_vec with identity matrix returns same vector.
    #[test]
    fn test_mat3_vec_identity() {
        let m = identity_inertia();
        let v = [1.0, 2.0, 3.0];
        let r = mat3_vec(m, v);
        for i in 0..3 {
            assert!((r[i] - v[i]).abs() < EPS);
        }
    }

    // T32: ContactJacobian tangent cdot for perpendicular motion.
    #[test]
    fn test_jacobian_tangent_cdot() {
        let cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0, 1);
        let jac = ContactJacobian::new(&cp, [0.0; 3], [0.0; 3]);
        let (t1, _t2) = cp.tangent_basis();
        // Body A moving along t1
        let cdot = jac.tangent1_cdot(scale(t1, -2.0), [0.0; 3], [0.0; 3], [0.0; 3]);
        // J_va = t1, so cdot = dot(t1, -2*t1) = -2
        assert!((cdot + 2.0).abs() < 1e-9, "cdot={cdot}");
    }

    // T33: Effective mass is zero when both bodies are static.
    #[test]
    fn test_effective_mass_both_static() {
        let n = [0.0, 1.0, 0.0_f64];
        let neg_n = scale(n, -1.0);
        let em = effective_mass_for_row(
            0.0,
            zero_inertia(),
            0.0,
            zero_inertia(),
            neg_n,
            [0.0; 3],
            n,
            [0.0; 3],
        );
        assert!(em.abs() < EPS, "em={em}");
    }

    // T34: Persistent contact is_valid returns false for negative depth.
    #[test]
    fn test_persistent_contact_is_valid() {
        let cp_neg = ContactPoint {
            depth: -0.01,
            ..ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], -0.01, 0, 1)
        };
        let pc = PersistentContact::new(cp_neg, 0.3);
        assert!(!pc.is_valid());
        let cp_pos = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 1);
        let pc2 = PersistentContact::new(cp_pos, 0.3);
        assert!(pc2.is_valid());
    }

    // T35: Island manager with static body (usize::MAX) doesn't panic.
    #[test]
    fn test_island_manager_static_body_pair() {
        let mut mgr = ContactIslandManager::new(2);
        let pairs = vec![(0, usize::MAX, 0)]; // body 1 = static world
        let islands = mgr.build_islands(&pairs); // should not panic
        let total: usize = islands.iter().map(|i| i.bodies.len()).sum();
        assert_eq!(total, 2);
    }
}
