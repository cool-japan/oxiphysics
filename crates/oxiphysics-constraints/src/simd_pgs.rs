// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Batch Projected Gauss-Seidel solver on a Structure-of-Arrays (SoA) layout.
//!
//! [`BatchPgsSolver`] solves scalar constraint rows between rigid bodies with
//! the projected Gauss-Seidel iteration.  Every row has a direction `n`, lever
//! arms `r_a` / `r_b` from the centres of mass to the point where the row acts,
//! and bounds on its accumulated impulse; the solver accounts for the linear
//! *and* the angular response of both bodies.
//!
//! ## Conventions
//!
//! A positive impulse `λ` pushes body B along `+n` and body A along `-n`
//! (applied at the lever arms), so the relative velocity of the two points
//! along the row,
//!
//! ```text
//! rv = n·(v_b − v_a) + (r_b × n)·ω_b − (r_a × n)·ω_a ,
//! ```
//!
//! changes by `K·λ` with the inverse effective mass
//!
//! ```text
//! K = m_a⁻¹ + m_b⁻¹ + (r_a × n)·I_a⁻¹(r_a × n) + (r_b × n)·I_b⁻¹(r_b × n) .
//! ```
//!
//! For a contact, `n` points from A to B, so `rv < 0` means the bodies
//! approach.  A row drives `rv + bias` to zero: `Δλ = −(rv + bias) / K`, and
//! the accumulated impulse is projected onto `[lo, hi]` after every update.
//!
//! ## Data layout
//!
//! Conventional PGS stores each constraint as a struct:
//! ```text
//! [n0x,n0y,n0z,eff_mass0, n1x,n1y,n1z,eff_mass1, ...]  ← AoS
//! ```
//! This module uses SoA:
//! ```text
//! nx: [n0x, n1x, n2x, n3x, ...]
//! ny: [n0y, n1y, n2y, n3y, ...]
//! nz: [n0z, n1z, n2z, n3z, ...]
//! eff_mass: [em0, em1, em2, em3, ...]
//! lambda:   [λ0,  λ1,  λ2,  λ3, ...]
//! ```
//! The lever-arm terms `r × n` and `I⁻¹(r × n)` are precomputed per row when it
//! is pushed, so the inner loop only needs dot products with the body
//! velocities.
//!
//! ## Batched rows
//!
//! The rows are processed in groups of four.  A group whose rows touch pairwise
//! disjoint sets of moving bodies is evaluated as a batch (all four relative
//! velocities first, then the four impulse updates), which is exactly what four
//! sequential updates would compute.  A group in which two rows share a moving
//! body is processed row by row instead, so the iteration is always the plain
//! Gauss-Seidel sweep over the rows in order and gives the same result as the
//! scalar [`PgsSolver`](crate::PgsSolver).  Rows that share only static bodies
//! (zero inverse mass and inertia) do not conflict.  The velocity gathers and
//! scatters through the body index arrays stay scalar; only the per-row
//! arithmetic of independent rows is laid out for vectorisation.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_constraints::simd_pgs::{BatchPgsSolver, SoaConstraintRow, SoaConstraints};
//!
//! let mut soa = SoaConstraints::new();
//! // A unit mass (body 1) falls at 2 m/s onto a static body (body 0).
//! soa.push(SoaConstraintRow::linear([0.0, 1.0, 0.0], 0, 1, 0.0, 1.0));
//! soa.set_body_velocity(1, [0.0, -2.0, 0.0], [0.0; 3]);
//!
//! let stats = BatchPgsSolver::default().solve(&mut soa, 1.0 / 60.0);
//! assert!(stats.converged);
//! // A contact impulse of 2 N·s stops the falling body.
//! assert!((soa.lambda[0] - 2.0).abs() < 1e-12);
//! assert!(soa.body_vy[1].abs() < 1e-12);
//! ```

/// Smallest inverse effective mass for which a row is solved; below it the
/// row's effective mass is set to zero and the row leaves the bodies alone.
const MIN_INVERSE_EFFECTIVE_MASS: f64 = 1e-15;

/// Penetration depth (m) tolerated by [`SoaContactBuilder`] before position
/// correction starts, unless changed with [`SoaContactBuilder::with_slop`].
const DEFAULT_CONTACT_SLOP: f64 = 0.005;

/// A 3x3 matrix in row-major order.
type Mat3 = [[f64; 3]; 3];

/// Cross product of two 3-vectors.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3-vectors.
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Matrix-vector product `m·v`.
#[inline]
fn mat_vec(m: &Mat3, v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// Project `x` onto `[lo, hi]`.  Unlike [`f64::clamp`] this never panics on
/// inverted or NaN bounds.
#[inline]
fn project(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ── SoaConstraints ──────────────────────────────────────────────────────────

/// A single constraint row for insertion into [`SoaConstraints::push`].
///
/// The inverse effective mass of the row is computed by
/// [`push`](SoaConstraints::push) from the inverse masses, the inverse
/// inertia tensors and the lever arms, so it is always consistent with them.
#[derive(Debug, Clone, Copy)]
pub struct SoaConstraintRow {
    /// Unit direction of the row in world space `[nx, ny, nz]` (for a contact:
    /// the contact normal pointing from body A to body B).
    pub normal: [f64; 3],
    /// Lever arm from body A's centre of mass to the point where the row acts
    /// (world space).
    pub r_a: [f64; 3],
    /// Lever arm from body B's centre of mass to the point where the row acts
    /// (world space).
    pub r_b: [f64; 3],
    /// Index of body A in the body velocity arrays of [`SoaConstraints`].
    pub body_a: usize,
    /// Index of body B in the body velocity arrays of [`SoaConstraints`].
    pub body_b: usize,
    /// Inverse mass of body A (`0.0` for static bodies).
    pub inv_mass_a: f64,
    /// Inverse mass of body B.
    pub inv_mass_b: f64,
    /// World-space inverse inertia tensor of body A (all zeros for a body that
    /// cannot rotate), row-major.
    pub inv_inertia_a: [[f64; 3]; 3],
    /// World-space inverse inertia tensor of body B, row-major.
    pub inv_inertia_b: [[f64; 3]; 3],
    /// Velocity bias for Baumgarte / restitution: the row drives
    /// `rv + bias` to zero.
    pub bias: f64,
    /// Warm-start accumulated impulse (`0.0` for cold start).  It is projected
    /// onto `[lo, hi]` and applied to the bodies by
    /// [`BatchPgsSolver::solve`] before the first iteration.
    pub warm_lambda: f64,
    /// Lower clamp for the accumulated impulse.
    pub lo: f64,
    /// Upper clamp for the accumulated impulse.
    pub hi: f64,
}

impl SoaConstraintRow {
    /// A contact row that acts through the centres of mass of two bodies: zero
    /// lever arms, no rotational response, no bias, cold start, and the
    /// non-pulling bounds `[0, ∞)`.
    ///
    /// Use struct update syntax to change the remaining fields, for example
    /// `SoaConstraintRow { bias: -0.1, ..SoaConstraintRow::linear(..) }`.
    pub fn linear(
        normal: [f64; 3],
        body_a: usize,
        body_b: usize,
        inv_mass_a: f64,
        inv_mass_b: f64,
    ) -> Self {
        Self {
            normal,
            r_a: [0.0; 3],
            r_b: [0.0; 3],
            body_a,
            body_b,
            inv_mass_a,
            inv_mass_b,
            inv_inertia_a: [[0.0; 3]; 3],
            inv_inertia_b: [[0.0; 3]; 3],
            bias: 0.0,
            warm_lambda: 0.0,
            lo: 0.0,
            hi: f64::INFINITY,
        }
    }
}

/// Structure-of-Arrays storage for one frame of PGS constraint data.
///
/// All per-row arrays are parallel: index `i` refers to constraint `i` across
/// all fields; fill them with [`push`](SoaConstraints::push) so that the
/// derived arrays stay consistent.  The per-body arrays are indexed by the body
/// indices of the rows and grow automatically when a row is pushed.
#[derive(Debug, Clone, Default)]
pub struct SoaConstraints {
    // Constraint directions (unit vectors) in world space
    /// Direction X components.
    pub nx: Vec<f64>,
    /// Direction Y components.
    pub ny: Vec<f64>,
    /// Direction Z components.
    pub nz: Vec<f64>,

    // Lever-arm vectors from body centre-of-mass to the point the row acts at
    /// r_a X (body A).
    pub rax: Vec<f64>,
    /// r_a Y (body A).
    pub ray_: Vec<f64>,
    /// r_a Z (body A).
    pub raz: Vec<f64>,
    /// r_b X (body B).
    pub rbx: Vec<f64>,
    /// r_b Y (body B).
    pub rby: Vec<f64>,
    /// r_b Z (body B).
    pub rbz: Vec<f64>,

    // Angular Jacobian of each row: r × n
    /// (r_a × n) X: how a rotation of body A changes the relative velocity.
    pub rxn_ax: Vec<f64>,
    /// (r_a × n) Y.
    pub rxn_ay: Vec<f64>,
    /// (r_a × n) Z.
    pub rxn_az: Vec<f64>,
    /// (r_b × n) X: how a rotation of body B changes the relative velocity.
    pub rxn_bx: Vec<f64>,
    /// (r_b × n) Y.
    pub rxn_by: Vec<f64>,
    /// (r_b × n) Z.
    pub rxn_bz: Vec<f64>,

    // Angular response of each row: I⁻¹ (r × n)
    /// (I_a⁻¹ (r_a × n)) X: angular velocity change of body A per unit impulse.
    pub irxn_ax: Vec<f64>,
    /// (I_a⁻¹ (r_a × n)) Y.
    pub irxn_ay: Vec<f64>,
    /// (I_a⁻¹ (r_a × n)) Z.
    pub irxn_az: Vec<f64>,
    /// (I_b⁻¹ (r_b × n)) X: angular velocity change of body B per unit impulse.
    pub irxn_bx: Vec<f64>,
    /// (I_b⁻¹ (r_b × n)) Y.
    pub irxn_by: Vec<f64>,
    /// (I_b⁻¹ (r_b × n)) Z.
    pub irxn_bz: Vec<f64>,

    // Body inverse masses
    /// Inverse mass of body A.
    pub inv_mass_a: Vec<f64>,
    /// Inverse mass of body B.
    pub inv_mass_b: Vec<f64>,

    // Effective mass: `1 / K`
    /// Effective mass (computed by `push` from masses, inertias and lever arms).
    pub eff_mass: Vec<f64>,

    // Velocity bias (Baumgarte stabilization + restitution)
    /// Velocity bias for position correction and restitution.
    pub bias: Vec<f64>,

    // Accumulated impulse (warm-start / output)
    /// Accumulated impulse λ (updated in place each iteration).
    pub lambda: Vec<f64>,

    // Impulse clamp bounds
    /// Lower impulse bound (e.g. `0.0` for contact, `NEG_INFINITY` for joints).
    pub lambda_lo: Vec<f64>,
    /// Upper impulse bound (e.g. `INFINITY`).
    pub lambda_hi: Vec<f64>,

    // Per-constraint relative velocity (recorded by the solver)
    /// Relative velocity along the row, as measured when the solver last
    /// visited the row (before the row's impulse update was applied).
    pub rel_vel: Vec<f64>,

    // Body indices of every row
    /// Index into the body arrays for body A.
    pub body_a_idx: Vec<usize>,
    /// Index into the body arrays for body B.
    pub body_b_idx: Vec<usize>,

    // Body velocities, one entry per body in the scene.  Multiple rows that
    // share a body update the same entry.
    /// Linear velocity X for each body in the scene (one per body).
    pub body_vx: Vec<f64>,
    /// Linear velocity Y.
    pub body_vy: Vec<f64>,
    /// Linear velocity Z.
    pub body_vz: Vec<f64>,
    /// Angular velocity X for each body in the scene (one per body).
    pub body_wx: Vec<f64>,
    /// Angular velocity Y.
    pub body_wy: Vec<f64>,
    /// Angular velocity Z.
    pub body_wz: Vec<f64>,
}

impl SoaConstraints {
    /// Create an empty constraint set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Return the number of constraints stored.
    pub fn len(&self) -> usize {
        self.nx.len()
    }

    /// Return `true` if no constraints are stored.
    pub fn is_empty(&self) -> bool {
        self.nx.is_empty()
    }

    /// Return the number of bodies the velocity arrays cover.
    pub fn body_count(&self) -> usize {
        self.body_vx.len()
    }

    /// Append one constraint row from a [`SoaConstraintRow`].
    ///
    /// Precomputes the angular terms `r × n` and `I⁻¹ (r × n)` of both bodies
    /// and the effective mass of the row, projects the warm-start impulse onto
    /// the row's bounds, and grows the body arrays so that both body indices
    /// are valid.
    pub fn push(&mut self, row: SoaConstraintRow) {
        let SoaConstraintRow {
            normal,
            r_a,
            r_b,
            body_a,
            body_b,
            inv_mass_a,
            inv_mass_b,
            inv_inertia_a,
            inv_inertia_b,
            bias,
            warm_lambda,
            lo,
            hi,
        } = row;

        let rxn_a = cross(r_a, normal);
        let rxn_b = cross(r_b, normal);
        let irxn_a = mat_vec(&inv_inertia_a, rxn_a);
        let irxn_b = mat_vec(&inv_inertia_b, rxn_b);
        let inverse_effective_mass =
            inv_mass_a + inv_mass_b + dot(rxn_a, irxn_a) + dot(rxn_b, irxn_b);
        let eff_mass = if inverse_effective_mass > MIN_INVERSE_EFFECTIVE_MASS {
            1.0 / inverse_effective_mass
        } else {
            0.0
        };

        self.nx.push(normal[0]);
        self.ny.push(normal[1]);
        self.nz.push(normal[2]);
        self.rax.push(r_a[0]);
        self.ray_.push(r_a[1]);
        self.raz.push(r_a[2]);
        self.rbx.push(r_b[0]);
        self.rby.push(r_b[1]);
        self.rbz.push(r_b[2]);
        self.rxn_ax.push(rxn_a[0]);
        self.rxn_ay.push(rxn_a[1]);
        self.rxn_az.push(rxn_a[2]);
        self.rxn_bx.push(rxn_b[0]);
        self.rxn_by.push(rxn_b[1]);
        self.rxn_bz.push(rxn_b[2]);
        self.irxn_ax.push(irxn_a[0]);
        self.irxn_ay.push(irxn_a[1]);
        self.irxn_az.push(irxn_a[2]);
        self.irxn_bx.push(irxn_b[0]);
        self.irxn_by.push(irxn_b[1]);
        self.irxn_bz.push(irxn_b[2]);
        self.inv_mass_a.push(inv_mass_a);
        self.inv_mass_b.push(inv_mass_b);
        self.eff_mass.push(eff_mass);
        self.bias.push(bias);
        self.lambda.push(project(warm_lambda, lo, hi));
        self.lambda_lo.push(lo);
        self.lambda_hi.push(hi);
        self.rel_vel.push(0.0);
        self.body_a_idx.push(body_a);
        self.body_b_idx.push(body_b);
        self.ensure_bodies(body_a.max(body_b).saturating_add(1));
    }

    /// Clear all constraints (retains allocations).
    ///
    /// The body velocities are kept; use [`resize_bodies`](Self::resize_bodies)
    /// to reset them.
    pub fn clear(&mut self) {
        self.nx.clear();
        self.ny.clear();
        self.nz.clear();
        self.rax.clear();
        self.ray_.clear();
        self.raz.clear();
        self.rbx.clear();
        self.rby.clear();
        self.rbz.clear();
        self.rxn_ax.clear();
        self.rxn_ay.clear();
        self.rxn_az.clear();
        self.rxn_bx.clear();
        self.rxn_by.clear();
        self.rxn_bz.clear();
        self.irxn_ax.clear();
        self.irxn_ay.clear();
        self.irxn_az.clear();
        self.irxn_bx.clear();
        self.irxn_by.clear();
        self.irxn_bz.clear();
        self.inv_mass_a.clear();
        self.inv_mass_b.clear();
        self.eff_mass.clear();
        self.bias.clear();
        self.lambda.clear();
        self.lambda_lo.clear();
        self.lambda_hi.clear();
        self.rel_vel.clear();
        self.body_a_idx.clear();
        self.body_b_idx.clear();
    }

    /// Resize the body velocity arrays to `n` bodies; new entries are
    /// zero-initialised and existing ones are kept.
    pub fn resize_bodies(&mut self, n: usize) {
        self.body_vx.resize(n, 0.0);
        self.body_vy.resize(n, 0.0);
        self.body_vz.resize(n, 0.0);
        self.body_wx.resize(n, 0.0);
        self.body_wy.resize(n, 0.0);
        self.body_wz.resize(n, 0.0);
    }

    /// Make sure the body arrays cover at least `count` bodies.
    fn ensure_bodies(&mut self, count: usize) {
        if self.body_vx.len() < count {
            self.resize_bodies(count);
        }
    }

    /// Grow the body arrays so that they cover every body index the rows refer
    /// to (they can fall short after [`resize_bodies`](Self::resize_bodies)
    /// shrank them or the index arrays were edited by hand).
    fn cover_row_bodies(&mut self) {
        let required = self
            .body_a_idx
            .iter()
            .chain(self.body_b_idx.iter())
            .copied()
            .max()
            .map_or(0, |highest| highest.saturating_add(1));
        self.ensure_bodies(required);
    }

    /// Whether all per-row arrays have the same length.
    ///
    /// [`push`](Self::push) keeps them consistent; arrays that were edited by
    /// hand may not be, and [`BatchPgsSolver::solve`] refuses to solve such a
    /// set.
    pub fn is_consistent(&self) -> bool {
        let rows = self.nx.len();
        [
            self.ny.len(),
            self.nz.len(),
            self.rax.len(),
            self.ray_.len(),
            self.raz.len(),
            self.rbx.len(),
            self.rby.len(),
            self.rbz.len(),
            self.rxn_ax.len(),
            self.rxn_ay.len(),
            self.rxn_az.len(),
            self.rxn_bx.len(),
            self.rxn_by.len(),
            self.rxn_bz.len(),
            self.irxn_ax.len(),
            self.irxn_ay.len(),
            self.irxn_az.len(),
            self.irxn_bx.len(),
            self.irxn_by.len(),
            self.irxn_bz.len(),
            self.inv_mass_a.len(),
            self.inv_mass_b.len(),
            self.eff_mass.len(),
            self.bias.len(),
            self.lambda.len(),
            self.lambda_lo.len(),
            self.lambda_hi.len(),
            self.rel_vel.len(),
            self.body_a_idx.len(),
            self.body_b_idx.len(),
        ]
        .iter()
        .all(|&len| len == rows)
            && self.body_vx.len() == self.body_vy.len()
            && self.body_vx.len() == self.body_vz.len()
            && self.body_vx.len() == self.body_wx.len()
            && self.body_vx.len() == self.body_wy.len()
            && self.body_vx.len() == self.body_wz.len()
    }

    /// Set the linear and angular velocity of `body`, growing the body arrays
    /// if the index is beyond them.
    pub fn set_body_velocity(&mut self, body: usize, linear: [f64; 3], angular: [f64; 3]) {
        self.ensure_bodies(body.saturating_add(1));
        self.body_vx[body] = linear[0];
        self.body_vy[body] = linear[1];
        self.body_vz[body] = linear[2];
        self.body_wx[body] = angular[0];
        self.body_wy[body] = angular[1];
        self.body_wz[body] = angular[2];
    }

    /// Linear and angular velocity of `body`, or `None` if the index is beyond
    /// the body arrays.
    pub fn body_velocity(&self, body: usize) -> Option<([f64; 3], [f64; 3])> {
        if body >= self.body_count() {
            return None;
        }
        Some((
            [self.body_vx[body], self.body_vy[body], self.body_vz[body]],
            [self.body_wx[body], self.body_wy[body], self.body_wz[body]],
        ))
    }

    /// Whether row `i` changes the velocity of body A.
    #[inline]
    fn moves_body_a(&self, i: usize) -> bool {
        self.inv_mass_a[i] != 0.0
            || self.irxn_ax[i] != 0.0
            || self.irxn_ay[i] != 0.0
            || self.irxn_az[i] != 0.0
    }

    /// Whether row `i` changes the velocity of body B.
    #[inline]
    fn moves_body_b(&self, i: usize) -> bool {
        self.inv_mass_b[i] != 0.0
            || self.irxn_bx[i] != 0.0
            || self.irxn_by[i] != 0.0
            || self.irxn_bz[i] != 0.0
    }

    /// Whether the four rows starting at `base` can be evaluated as a batch:
    /// no body that one of them moves is read or moved by another one.
    fn lanes_are_independent(&self, base: usize) -> bool {
        for first in base..base + 4 {
            for second in first + 1..base + 4 {
                if self.rows_conflict(first, second) {
                    return false;
                }
            }
        }
        true
    }

    /// Whether rows `r` and `s` are coupled through a body that either of them
    /// moves.
    fn rows_conflict(&self, r: usize, s: usize) -> bool {
        let touches =
            |row: usize, body: usize| self.body_a_idx[row] == body || self.body_b_idx[row] == body;
        (self.moves_body_a(r) && touches(s, self.body_a_idx[r]))
            || (self.moves_body_b(r) && touches(s, self.body_b_idx[r]))
            || (self.moves_body_a(s) && touches(r, self.body_a_idx[s]))
            || (self.moves_body_b(s) && touches(r, self.body_b_idx[s]))
    }
}

// ── BatchPgsSolverStats ──────────────────────────────────────────────────────

/// Statistics returned by [`BatchPgsSolver::solve`].
#[derive(Debug, Clone, Copy)]
pub struct BatchSolverStats {
    /// Number of iterations actually performed.
    pub iterations_used: usize,
    /// Final sum of |Δλ| across all constraints in the last iteration.
    pub residual: f64,
    /// Whether the solver converged before the iteration limit.
    pub converged: bool,
}

// ── BatchPgsSolver ───────────────────────────────────────────────────────────

/// Batch PGS solver that operates on [`SoaConstraints`].
///
/// Every iteration is one Gauss-Seidel sweep over the rows in order: a row
/// measures the relative velocity of its two bodies at the lever arms (linear
/// and angular part), updates its accumulated impulse
/// `λ ← project(λ − ω · (rv + bias) / K)` and applies the change to the
/// linear and angular velocities of both bodies.  Groups of four rows that do
/// not share a moving body are evaluated as a batch, which gives exactly the
/// same result (see the [module documentation](self)).
///
/// Friction constraints are typically solved in a second pass after the normal
/// (non-penetration) pass; configure separate [`SoaConstraints`] for each.
#[derive(Debug, Clone)]
pub struct BatchPgsSolver {
    /// Maximum number of velocity-level iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on sum |Δλ|.
    pub tolerance: f64,
    /// SOR relaxation factor ω ∈ (0, 2).  1.0 = standard PGS, 1.3 typical.
    pub omega: f64,
}

impl Default for BatchPgsSolver {
    fn default() -> Self {
        Self {
            max_iterations: 10,
            tolerance: 1e-6,
            omega: 1.0,
        }
    }
}

impl BatchPgsSolver {
    /// Create a new solver with explicit parameters.
    pub fn new(max_iterations: usize, tolerance: f64, omega: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
            omega,
        }
    }

    /// Run the PGS velocity iteration on `constraints`.
    ///
    /// The accumulated impulses in [`SoaConstraints::lambda`] are the warm
    /// start of this solve: before the first iteration every `lambda[i]` is
    /// applied to the velocities of the two bodies of row `i` (`+λ` along the
    /// row on body B, `-λ` on body A, linear and angular part), so the body
    /// velocities passed in must be the velocities *without* those impulses.
    /// The iteration then changes `lambda` and the velocities consistently, and
    /// the final `lambda` can be used for warm-starting the next frame.
    ///
    /// The body arrays are grown to cover every body index of the rows.  A set
    /// whose per-row arrays differ in length (see
    /// [`SoaConstraints::is_consistent`]) is not solved: nothing is touched and
    /// the returned statistics say that no iteration ran and nothing converged.
    ///
    /// Returns solver statistics.
    pub fn solve(&self, c: &mut SoaConstraints, _dt: f64) -> BatchSolverStats {
        if !c.is_consistent() {
            return BatchSolverStats {
                iterations_used: 0,
                residual: f64::MAX,
                converged: false,
            };
        }
        let n = c.len();
        if n == 0 {
            return BatchSolverStats {
                iterations_used: 0,
                residual: 0.0,
                converged: true,
            };
        }

        c.cover_row_bodies();
        Self::warm_start(c);
        let independent: Vec<bool> = (0..n / 4)
            .map(|chunk| c.lanes_are_independent(chunk * 4))
            .collect();

        let mut residual = f64::MAX;
        let mut iters = 0usize;

        for _iter in 0..self.max_iterations {
            residual = self.velocity_iteration(c, &independent);
            iters += 1;
            if residual < self.tolerance {
                break;
            }
        }

        BatchSolverStats {
            iterations_used: iters,
            residual,
            converged: residual < self.tolerance,
        }
    }

    /// Apply the accumulated impulse of every row to the body velocities.
    fn warm_start(c: &mut SoaConstraints) {
        for i in 0..c.len() {
            let accumulated = c.lambda[i];
            if accumulated != 0.0 {
                c.apply_row_impulse(i, accumulated);
            }
        }
    }

    /// One velocity-level PGS iteration over all constraints.
    ///
    /// `independent[k]` tells whether the four rows `4k..4k+4` can be
    /// evaluated as a batch.
    fn velocity_iteration(&self, c: &mut SoaConstraints, independent: &[bool]) -> f64 {
        let n = c.len();
        let chunks = n / 4;
        let mut total_delta = 0.0f64;

        // ── Groups of four rows ──────────────────────────────────────────────
        for chunk in 0..chunks {
            let base = chunk * 4;
            if independent.get(chunk).copied().unwrap_or(false) {
                for delta in self.solve_lanes(c, base) {
                    total_delta += delta.abs();
                }
            } else {
                for i in base..base + 4 {
                    total_delta += self.solve_row(c, i);
                }
            }
        }

        // ── Scalar tail for remaining constraints ────────────────────────────
        for i in (chunks * 4)..n {
            total_delta += self.solve_row(c, i);
        }

        total_delta
    }

    /// Evaluate four rows that share no moving body as a batch: all four
    /// relative velocities first, then the four impulse updates.
    ///
    /// Returns the applied impulse change of every row.
    #[inline]
    fn solve_lanes(&self, c: &mut SoaConstraints, base: usize) -> [f64; 4] {
        let rv: [f64; 4] = std::array::from_fn(|lane| c.relative_velocity(base + lane));
        let mut applied = [0.0f64; 4];
        for (lane, (slot, rel_vel)) in applied.iter_mut().zip(rv).enumerate() {
            *slot = self.update_lambda(c, base + lane, rel_vel);
        }
        for (lane, delta) in applied.into_iter().enumerate() {
            c.apply_row_impulse(base + lane, delta);
        }
        applied
    }

    /// Solve one row on its own; returns the size of its impulse change.
    #[inline]
    fn solve_row(&self, c: &mut SoaConstraints, i: usize) -> f64 {
        let rel_vel = c.relative_velocity(i);
        let applied = self.update_lambda(c, i, rel_vel);
        c.apply_row_impulse(i, applied);
        applied.abs()
    }

    /// Update the accumulated impulse of row `i` for the measured relative
    /// velocity and return the change that was applied after clamping.
    #[inline]
    fn update_lambda(&self, c: &mut SoaConstraints, i: usize, rel_vel: f64) -> f64 {
        // delta_lambda = -eff_mass * (rv + bias) * omega
        let old = c.lambda[i];
        let delta = -c.eff_mass[i] * (rel_vel + c.bias[i]) * self.omega;
        let new = project(old + delta, c.lambda_lo[i], c.lambda_hi[i]);
        c.lambda[i] = new;
        c.rel_vel[i] = rel_vel;
        new - old
    }
}

impl SoaConstraints {
    /// Relative velocity of the two bodies of row `i` along the row, including
    /// the rotation of both bodies about their lever arms:
    /// `n·(v_b − v_a) + (r_b × n)·ω_b − (r_a × n)·ω_a`.
    #[inline]
    fn relative_velocity(&self, i: usize) -> f64 {
        let a = self.body_a_idx[i];
        let b = self.body_b_idx[i];
        let linear = self.nx[i] * (self.body_vx[b] - self.body_vx[a])
            + self.ny[i] * (self.body_vy[b] - self.body_vy[a])
            + self.nz[i] * (self.body_vz[b] - self.body_vz[a]);
        let angular_b = self.rxn_bx[i] * self.body_wx[b]
            + self.rxn_by[i] * self.body_wy[b]
            + self.rxn_bz[i] * self.body_wz[b];
        let angular_a = self.rxn_ax[i] * self.body_wx[a]
            + self.rxn_ay[i] * self.body_wy[a]
            + self.rxn_az[i] * self.body_wz[a];
        linear + angular_b - angular_a
    }

    /// Apply the impulse change `delta_lambda` of row `i` to the linear and
    /// angular velocities of its two bodies: `+Δλ` on body B, `-Δλ` on body A.
    #[inline]
    fn apply_row_impulse(&mut self, i: usize, delta_lambda: f64) {
        let a = self.body_a_idx[i];
        let b = self.body_b_idx[i];
        let nx = self.nx[i];
        let ny = self.ny[i];
        let nz = self.nz[i];
        let ima = self.inv_mass_a[i];
        let imb = self.inv_mass_b[i];
        // vb += inv_mass_b * delta_lambda * n ; wb += delta_lambda * I_b^-1 (r_b × n)
        self.body_vx[b] += imb * delta_lambda * nx;
        self.body_vy[b] += imb * delta_lambda * ny;
        self.body_vz[b] += imb * delta_lambda * nz;
        self.body_wx[b] += delta_lambda * self.irxn_bx[i];
        self.body_wy[b] += delta_lambda * self.irxn_by[i];
        self.body_wz[b] += delta_lambda * self.irxn_bz[i];
        // va -= inv_mass_a * delta_lambda * n ; wa -= delta_lambda * I_a^-1 (r_a × n)
        self.body_vx[a] -= ima * delta_lambda * nx;
        self.body_vy[a] -= ima * delta_lambda * ny;
        self.body_vz[a] -= ima * delta_lambda * nz;
        self.body_wx[a] -= delta_lambda * self.irxn_ax[i];
        self.body_wy[a] -= delta_lambda * self.irxn_ay[i];
        self.body_wz[a] -= delta_lambda * self.irxn_az[i];
    }
}

// ── Convenience builder ──────────────────────────────────────────────────────

/// One contact of a rigid-body scene, in the form
/// [`SoaContactBuilder::add_contact`] consumes.
#[derive(Debug, Clone, Copy)]
pub struct SoaContact {
    /// Index of body A (the body the normal points away from).
    pub body_a: usize,
    /// Index of body B (the body the normal points towards).
    pub body_b: usize,
    /// World-space contact normal (A→B direction).
    pub normal: [f64; 3],
    /// Signed penetration depth (positive = overlapping).
    pub penetration: f64,
    /// Lever arm from body A's centre of mass to the contact point.
    pub r_a: [f64; 3],
    /// Lever arm from body B's centre of mass to the contact point.
    pub r_b: [f64; 3],
    /// Inverse mass of body A (`0.0` for static bodies).
    pub inv_mass_a: f64,
    /// Inverse mass of body B.
    pub inv_mass_b: f64,
    /// World-space inverse inertia tensor of body A (zeros for a static body).
    pub inv_inertia_a: [[f64; 3]; 3],
    /// World-space inverse inertia tensor of body B.
    pub inv_inertia_b: [[f64; 3]; 3],
    /// Coefficient of restitution ∈ \[0,1\].
    pub restitution: f64,
    /// Relative velocity of the contact points along the normal before the
    /// solve, `n·(v_b + ω_b × r_b − v_a − ω_a × r_a)`; negative when the bodies
    /// approach.  Only an approaching contact bounces.
    pub approach_vel: f64,
    /// Impulse accumulated for this contact in the previous frame (`0.0` for a
    /// cold start).
    pub warm_lambda: f64,
}

impl SoaContact {
    /// A contact through the centres of mass of two bodies without rotational
    /// response, with no restitution and a cold start.  Use struct update
    /// syntax to set the remaining fields.
    pub fn centered(
        body_a: usize,
        body_b: usize,
        normal: [f64; 3],
        penetration: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            normal,
            penetration,
            r_a: [0.0; 3],
            r_b: [0.0; 3],
            inv_mass_a,
            inv_mass_b,
            inv_inertia_a: [[0.0; 3]; 3],
            inv_inertia_b: [[0.0; 3]; 3],
            restitution: 0.0,
            approach_vel: 0.0,
            warm_lambda: 0.0,
        }
    }
}

/// Builder that converts contacts into [`SoaConstraints`] ready for
/// [`BatchPgsSolver`].
///
/// Each contact becomes one non-pulling row (`λ ≥ 0`) whose bias combines the
/// Baumgarte position correction with restitution:
///
/// ```text
/// bias = −(β/dt)·max(penetration − slop, 0) + e·min(approach_vel, 0)
/// ```
///
/// The row drives `rv + bias` to zero, so a penetrating contact is pushed apart
/// at `(β/dt)·(penetration − slop)` and an approaching contact leaves with
/// `e·|approach_vel|` on top of that.
#[derive(Debug)]
pub struct SoaContactBuilder {
    constraints: SoaConstraints,
    slop: f64,
    restitution_threshold: f64,
}

impl Default for SoaContactBuilder {
    fn default() -> Self {
        Self {
            constraints: SoaConstraints::new(),
            slop: DEFAULT_CONTACT_SLOP,
            restitution_threshold: 0.0,
        }
    }
}

impl SoaContactBuilder {
    /// Create a new empty builder (slop 5 mm, restitution at every approach
    /// speed).
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the penetration depth (m) that is tolerated before the Baumgarte
    /// position correction starts.
    pub fn with_slop(mut self, slop: f64) -> Self {
        self.slop = slop.max(0.0);
        self
    }

    /// Set the approach speed (m/s) below which restitution is ignored, which
    /// keeps resting contacts from jittering.
    pub fn with_restitution_threshold(mut self, threshold: f64) -> Self {
        self.restitution_threshold = threshold.max(0.0);
        self
    }

    /// Add a contact constraint.
    ///
    /// * `contact`   — bodies, normal (A→B), lever arms, masses, restitution
    ///   and the approach velocity (see [`SoaContact`])
    /// * `dt`        — simulation time step
    /// * `baumgarte` — stabilization factor (typical 0.2)
    ///
    /// The row starts from the contact's `warm_lambda` (zero for a cold start);
    /// restitution enters through the bias only, never through the warm start.
    pub fn add_contact(&mut self, contact: SoaContact, dt: f64, baumgarte: f64) {
        let position_bias = if dt > 0.0 {
            -(baumgarte / dt) * (contact.penetration - self.slop).max(0.0)
        } else {
            0.0
        };
        let restitution_bias = if contact.approach_vel < -self.restitution_threshold {
            contact.restitution * contact.approach_vel.min(0.0)
        } else {
            0.0
        };
        self.constraints.push(SoaConstraintRow {
            normal: contact.normal,
            r_a: contact.r_a,
            r_b: contact.r_b,
            body_a: contact.body_a,
            body_b: contact.body_b,
            inv_mass_a: contact.inv_mass_a,
            inv_mass_b: contact.inv_mass_b,
            inv_inertia_a: contact.inv_inertia_a,
            inv_inertia_b: contact.inv_inertia_b,
            bias: position_bias + restitution_bias,
            warm_lambda: contact.warm_lambda,
            lo: 0.0, // contact: no pulling
            hi: f64::INFINITY,
        });
    }

    /// Consume the builder, returning the constraint set.
    pub fn build(self) -> SoaConstraints {
        self.constraints
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
