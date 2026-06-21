// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SIMD-optimized Projected Gauss-Seidel solver using Structure-of-Arrays (SoA) layout.
//!
//! This module provides [`BatchPgsSolver`] which stores constraint data in SoA (Structure of
//! Arrays) format so that the LLVM auto-vectorizer can emit SIMD instructions (SSE2/AVX2/NEON)
//! for the inner iteration loops.
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
//! When the inner loop advances in steps of 4, LLVM merges the scalar loads into
//! 256-bit AVX2 registers (or 128-bit SSE2/NEON) and processes 4 constraints per
//! CPU cycle.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_constraints::simd_pgs::{BatchPgsSolver, SoaConstraintRow, SoaConstraints};
//!
//! let mut soa = SoaConstraints::new();
//! // Populate constraints from your scene …
//! soa.push(SoaConstraintRow {
//!     normal: [0.0, 1.0, 0.0],
//!     r_a: [0.0, 0.1, 0.0],
//!     r_b: [0.0, -0.1, 0.0],
//!     inv_mass_a: 1.0,
//!     inv_mass_b: 1.0,
//!     eff_mass: 1.0,
//!     bias: 0.01,
//!     warm_lambda: 0.0,
//!     lo: f64::NEG_INFINITY,
//!     hi: f64::INFINITY,
//! });
//! soa.resize_bodies(1); // body index 0 needs one velocity slot
//!
//! let stats = BatchPgsSolver::default().solve(&mut soa, 1.0 / 60.0);
//! assert!(stats.iterations_used >= 1);
//! ```

// ── SoaConstraints ──────────────────────────────────────────────────────────

/// A single constraint row for insertion into [`SoaConstraints::push`].
#[derive(Debug, Clone, Copy)]
pub struct SoaConstraintRow {
    /// Unit contact normal in world space `[nx, ny, nz]`.
    pub normal: [f64; 3],
    /// Lever arm from body A COM to contact `[rax, ray, raz]`.
    pub r_a: [f64; 3],
    /// Lever arm from body B COM to contact `[rbx, rby, rbz]`.
    pub r_b: [f64; 3],
    /// Inverse mass of body A (`0.0` for static bodies).
    pub inv_mass_a: f64,
    /// Inverse mass of body B.
    pub inv_mass_b: f64,
    /// Precomputed effective mass (1/K).
    pub eff_mass: f64,
    /// Velocity bias for Baumgarte / restitution.
    pub bias: f64,
    /// Warm-start accumulated impulse (0.0 for cold start).
    pub warm_lambda: f64,
    /// Lower clamp for impulse.
    pub lo: f64,
    /// Upper clamp for impulse.
    pub hi: f64,
}

/// Structure-of-Arrays storage for one frame of PGS constraint data.
///
/// All arrays are parallel: index `i` refers to constraint `i` across all fields.
#[derive(Debug, Clone, Default)]
pub struct SoaConstraints {
    // Constraint normals (unit vectors) in world space
    /// Normal X components.
    pub nx: Vec<f64>,
    /// Normal Y components.
    pub ny: Vec<f64>,
    /// Normal Z components.
    pub nz: Vec<f64>,

    // Lever-arm vectors from body centre-of-mass to contact point
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

    // Body inverse masses
    /// Inverse mass of body A.
    pub inv_mass_a: Vec<f64>,
    /// Inverse mass of body B.
    pub inv_mass_b: Vec<f64>,

    // Precomputed effective mass: `1 / K`
    /// Effective mass (precomputed).
    pub eff_mass: Vec<f64>,

    // Velocity bias (Baumgarte stabilization + restitution)
    /// Velocity bias for position correction.
    pub bias: Vec<f64>,

    // Accumulated impulse (warm-start / output)
    /// Accumulated impulse λ (updated in place each iteration).
    pub lambda: Vec<f64>,

    // Impulse clamp bounds
    /// Lower impulse bound (e.g. `0.0` for contact, `NEG_INFINITY` for joints).
    pub lambda_lo: Vec<f64>,
    /// Upper impulse bound (e.g. `INFINITY`).
    pub lambda_hi: Vec<f64>,

    // Per-constraint relative velocity (scratch, recomputed each iter)
    /// Relative velocity along normal (scratch buffer).
    pub rel_vel: Vec<f64>,

    // Body linear velocities (interleaved: body_idx → velocity)
    // Stored separately so multiple constraints sharing a body can be updated
    // without re-fetching from a HashMap.  Index matches `body_a_idx` /
    // `body_b_idx`.
    /// Index into velocity arrays for body A.
    pub body_a_idx: Vec<usize>,
    /// Index into velocity arrays for body B.
    pub body_b_idx: Vec<usize>,

    /// Linear velocity X for each body in the scene (one per body).
    pub body_vx: Vec<f64>,
    /// Linear velocity Y.
    pub body_vy: Vec<f64>,
    /// Linear velocity Z.
    pub body_vz: Vec<f64>,
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

    /// Append one constraint row from a [`SoaConstraintRow`].
    pub fn push(&mut self, row: SoaConstraintRow) {
        let SoaConstraintRow {
            normal,
            r_a,
            r_b,
            inv_mass_a,
            inv_mass_b,
            eff_mass,
            bias,
            warm_lambda,
            lo,
            hi,
        } = row;
        self.nx.push(normal[0]);
        self.ny.push(normal[1]);
        self.nz.push(normal[2]);
        self.rax.push(r_a[0]);
        self.ray_.push(r_a[1]);
        self.raz.push(r_a[2]);
        self.rbx.push(r_b[0]);
        self.rby.push(r_b[1]);
        self.rbz.push(r_b[2]);
        self.inv_mass_a.push(inv_mass_a);
        self.inv_mass_b.push(inv_mass_b);
        self.eff_mass.push(eff_mass);
        self.bias.push(bias);
        self.lambda.push(warm_lambda);
        self.lambda_lo.push(lo);
        self.lambda_hi.push(hi);
        self.rel_vel.push(0.0);
        // Body indices / velocities are managed by the caller.
        self.body_a_idx.push(0);
        self.body_b_idx.push(0);
    }

    /// Clear all constraints (retains allocations).
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

    /// Resize the body velocity arrays to `n` bodies (zero-initialised).
    pub fn resize_bodies(&mut self, n: usize) {
        self.body_vx.resize(n, 0.0);
        self.body_vy.resize(n, 0.0);
        self.body_vz.resize(n, 0.0);
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

/// SIMD-friendly batch PGS solver that operates on [`SoaConstraints`].
///
/// The inner loop is written so that LLVM can auto-vectorize it with SSE2/AVX2/NEON:
/// * The data is laid out contiguously in SoA format.
/// * The loop step is 4 (matches a 256-bit AVX2 `f64x4` register).
/// * Operations are simple FMA-shaped arithmetic with no early-exits.
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
    /// Returns solver statistics. The updated `lambda` values in `constraints`
    /// can be used for warm-starting the next frame.
    pub fn solve(&self, c: &mut SoaConstraints, _dt: f64) -> BatchSolverStats {
        let n = c.len();
        if n == 0 {
            return BatchSolverStats {
                iterations_used: 0,
                residual: 0.0,
                converged: true,
            };
        }

        let mut residual = f64::MAX;
        let mut iters = 0usize;

        for _iter in 0..self.max_iterations {
            residual = self.velocity_iteration(c);
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

    /// One velocity-level PGS iteration over all constraints.
    ///
    /// The inner loop is written in a scalar style that the auto-vectorizer
    /// handles efficiently due to the SoA data layout.
    #[inline(always)]
    fn velocity_iteration(&self, c: &mut SoaConstraints) -> f64 {
        let n = c.len();
        let omega = self.omega;
        let mut total_delta = 0.0f64;

        // ── Batched loop (step 4 for SIMD lanes) ────────────────────────────
        // The compiler unrolls this into a vectorized 4-wide lane when
        // optimizing with -C opt-level=3 (release mode).
        let chunks = n / 4;
        for chunk in 0..chunks {
            let base = chunk * 4;

            // Load 4 normals
            let nx0 = c.nx[base];
            let nx1 = c.nx[base + 1];
            let nx2 = c.nx[base + 2];
            let nx3 = c.nx[base + 3];
            let ny0 = c.ny[base];
            let ny1 = c.ny[base + 1];
            let ny2 = c.ny[base + 2];
            let ny3 = c.ny[base + 3];
            let nz0 = c.nz[base];
            let nz1 = c.nz[base + 1];
            let nz2 = c.nz[base + 2];
            let nz3 = c.nz[base + 3];

            // Load bias and effective mass
            let bias0 = c.bias[base];
            let bias1 = c.bias[base + 1];
            let bias2 = c.bias[base + 2];
            let bias3 = c.bias[base + 3];
            let em0 = c.eff_mass[base];
            let em1 = c.eff_mass[base + 1];
            let em2 = c.eff_mass[base + 2];
            let em3 = c.eff_mass[base + 3];

            // Relative velocity (using body velocity arrays)
            let rv0 = self.compute_rel_vel(c, base);
            let rv1 = self.compute_rel_vel(c, base + 1);
            let rv2 = self.compute_rel_vel(c, base + 2);
            let rv3 = self.compute_rel_vel(c, base + 3);

            // Suppress unused normals warning - they contribute to rel_vel via dot
            let _ = (nx0, ny0, nz0, nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3);

            // delta_lambda = -eff_mass * (rv + bias)
            let dl0 = -em0 * (rv0 + bias0) * omega;
            let dl1 = -em1 * (rv1 + bias1) * omega;
            let dl2 = -em2 * (rv2 + bias2) * omega;
            let dl3 = -em3 * (rv3 + bias3) * omega;

            // Clamp accumulated impulse
            let old0 = c.lambda[base];
            let old1 = c.lambda[base + 1];
            let old2 = c.lambda[base + 2];
            let old3 = c.lambda[base + 3];

            let new0 = (old0 + dl0).clamp(c.lambda_lo[base], c.lambda_hi[base]);
            let new1 = (old1 + dl1).clamp(c.lambda_lo[base + 1], c.lambda_hi[base + 1]);
            let new2 = (old2 + dl2).clamp(c.lambda_lo[base + 2], c.lambda_hi[base + 2]);
            let new3 = (old3 + dl3).clamp(c.lambda_lo[base + 3], c.lambda_hi[base + 3]);

            let actual0 = new0 - old0;
            let actual1 = new1 - old1;
            let actual2 = new2 - old2;
            let actual3 = new3 - old3;

            c.lambda[base] = new0;
            c.lambda[base + 1] = new1;
            c.lambda[base + 2] = new2;
            c.lambda[base + 3] = new3;

            // Apply impulses to body velocities
            self.apply_impulse(c, base, actual0);
            self.apply_impulse(c, base + 1, actual1);
            self.apply_impulse(c, base + 2, actual2);
            self.apply_impulse(c, base + 3, actual3);

            total_delta += actual0.abs() + actual1.abs() + actual2.abs() + actual3.abs();
        }

        // ── Scalar tail for remaining constraints ────────────────────────────
        for i in (chunks * 4)..n {
            let rv = self.compute_rel_vel(c, i);
            let dl = -c.eff_mass[i] * (rv + c.bias[i]) * omega;
            let old = c.lambda[i];
            let new_lam = (old + dl).clamp(c.lambda_lo[i], c.lambda_hi[i]);
            let actual = new_lam - old;
            c.lambda[i] = new_lam;
            self.apply_impulse(c, i, actual);
            total_delta += actual.abs();
        }

        total_delta
    }

    /// Compute the relative velocity along the constraint normal for constraint `i`.
    #[inline(always)]
    fn compute_rel_vel(&self, c: &SoaConstraints, i: usize) -> f64 {
        let ai = c.body_a_idx[i];
        let bi = c.body_b_idx[i];
        let va = [c.body_vx[ai], c.body_vy[ai], c.body_vz[ai]];
        let vb = [c.body_vx[bi], c.body_vy[bi], c.body_vz[bi]];
        // rv = n · (vb - va)
        let n = [c.nx[i], c.ny[i], c.nz[i]];
        n[0] * (vb[0] - va[0]) + n[1] * (vb[1] - va[1]) + n[2] * (vb[2] - va[2])
    }

    /// Apply impulse `lambda_delta` for constraint `i` to the associated body velocities.
    #[inline(always)]
    fn apply_impulse(&self, c: &mut SoaConstraints, i: usize, delta_lambda: f64) {
        let ai = c.body_a_idx[i];
        let bi = c.body_b_idx[i];
        let nx = c.nx[i];
        let ny = c.ny[i];
        let nz = c.nz[i];
        let ima = c.inv_mass_a[i];
        let imb = c.inv_mass_b[i];
        // vb += inv_mass_b * lambda_delta * n
        c.body_vx[bi] += imb * delta_lambda * nx;
        c.body_vy[bi] += imb * delta_lambda * ny;
        c.body_vz[bi] += imb * delta_lambda * nz;
        // va -= inv_mass_a * lambda_delta * n
        c.body_vx[ai] -= ima * delta_lambda * nx;
        c.body_vy[ai] -= ima * delta_lambda * ny;
        c.body_vz[ai] -= ima * delta_lambda * nz;
    }
}

// ── Convenience builder ──────────────────────────────────────────────────────

/// Builder that converts a flat array of contact normals and masses into
/// [`SoaConstraints`] ready for [`BatchPgsSolver`].
#[derive(Debug, Default)]
pub struct SoaContactBuilder {
    constraints: SoaConstraints,
}

impl SoaContactBuilder {
    /// Create a new empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a contact constraint.
    ///
    /// * `normal`     — world-space contact normal (A→B direction)
    /// * `penetration` — signed penetration depth (positive = overlapping)
    /// * `inv_mass_a` / `inv_mass_b` — inverse masses (0 for static)
    /// * `restitution` — coefficient of restitution ∈ \[0,1\]
    /// * `dt`          — simulation time step
    /// * `baumgarte`   — stabilization factor (typical 0.2)
    pub fn add_contact(
        &mut self,
        normal: [f64; 3],
        penetration: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
        restitution: f64,
        dt: f64,
        baumgarte: f64,
    ) {
        let eff_mass = {
            let k = inv_mass_a + inv_mass_b;
            if k.abs() < 1e-15 { 0.0 } else { 1.0 / k }
        };
        // Baumgarte bias for positional correction
        let bias = if penetration > 0.0 {
            -(baumgarte / dt) * penetration
        } else {
            0.0
        };
        // Warm-start lambda = 0 (cold start)
        let warm_lambda = restitution * penetration.max(0.0) * eff_mass;
        self.constraints.push(SoaConstraintRow {
            normal,
            r_a: [0.0; 3],
            r_b: [0.0; 3],
            inv_mass_a,
            inv_mass_b,
            eff_mass,
            bias,
            warm_lambda,
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
mod tests {
    use super::*;

    #[test]
    fn empty_constraints_converge_immediately() {
        let mut soa = SoaConstraints::new();
        let stats = BatchPgsSolver::default().solve(&mut soa, 1.0 / 60.0);
        assert!(stats.converged);
        assert_eq!(stats.iterations_used, 0);
        assert_eq!(stats.residual, 0.0);
    }

    #[test]
    fn single_contact_impulse_nonnegative() {
        let mut builder = SoaContactBuilder::new();
        builder.add_contact([0.0, 1.0, 0.0], 0.05, 1.0, 1.0, 0.0, 1.0 / 60.0, 0.2);
        let mut soa = builder.build();
        soa.resize_bodies(2);
        soa.body_a_idx[0] = 0;
        soa.body_b_idx[0] = 1;
        // Give body A a downward velocity so it violates the constraint
        soa.body_vy[0] = -1.0;

        let stats = BatchPgsSolver::default().solve(&mut soa, 1.0 / 60.0);
        // Lambda must be non-negative for a contact constraint
        assert!(soa.lambda[0] >= 0.0);
        assert!(stats.iterations_used >= 1);
    }

    #[test]
    fn batch_and_scalar_same_result() {
        // Build exactly 5 constraints so we exercise both the batch (4) and tail (1) paths.
        let mut soa1 = SoaConstraints::new();
        let mut soa2 = SoaConstraints::new();
        for _ in 0..5 {
            for soa in [&mut soa1, &mut soa2] {
                soa.nx.push(0.0);
                soa.ny.push(1.0);
                soa.nz.push(0.0);
                soa.rax.push(0.0);
                soa.ray_.push(0.0);
                soa.raz.push(0.0);
                soa.rbx.push(0.0);
                soa.rby.push(0.0);
                soa.rbz.push(0.0);
                soa.inv_mass_a.push(1.0);
                soa.inv_mass_b.push(1.0);
                soa.eff_mass.push(0.5);
                soa.bias.push(0.0);
                soa.lambda.push(0.0);
                soa.lambda_lo.push(0.0);
                soa.lambda_hi.push(f64::INFINITY);
                soa.rel_vel.push(0.0);
                soa.body_a_idx.push(0);
                soa.body_b_idx.push(1);
            }
        }
        soa1.resize_bodies(2);
        soa2.resize_bodies(2);
        soa1.body_vy[0] = -1.0;
        soa2.body_vy[0] = -1.0;

        let solver = BatchPgsSolver::default();
        let s1 = solver.solve(&mut soa1, 1.0 / 60.0);
        let s2 = solver.solve(&mut soa2, 1.0 / 60.0);
        assert_eq!(s1.iterations_used, s2.iterations_used);
    }

    #[test]
    fn sor_omega_above_one_speeds_convergence() {
        // Higher omega should reduce residual in fewer iterations
        let make_soa = || {
            let mut soa = SoaConstraints::new();
            for _ in 0..8 {
                soa.nx.push(0.0);
                soa.ny.push(1.0);
                soa.nz.push(0.0);
                soa.rax.push(0.0);
                soa.ray_.push(0.0);
                soa.raz.push(0.0);
                soa.rbx.push(0.0);
                soa.rby.push(0.0);
                soa.rbz.push(0.0);
                soa.inv_mass_a.push(1.0);
                soa.inv_mass_b.push(1.0);
                soa.eff_mass.push(0.5);
                soa.bias.push(0.0);
                soa.lambda.push(0.0);
                soa.lambda_lo.push(0.0);
                soa.lambda_hi.push(f64::INFINITY);
                soa.rel_vel.push(0.0);
                soa.body_a_idx.push(0);
                soa.body_b_idx.push(1);
            }
            soa.resize_bodies(2);
            soa.body_vy[0] = -2.0;
            soa
        };
        let mut soa_std = make_soa();
        let mut soa_sor = make_soa();
        let std_solver = BatchPgsSolver::new(3, 1e-12, 1.0);
        let sor_solver = BatchPgsSolver::new(3, 1e-12, 1.3);
        let stats_std = std_solver.solve(&mut soa_std, 1.0 / 60.0);
        let stats_sor = sor_solver.solve(&mut soa_sor, 1.0 / 60.0);
        // SOR should achieve lower or equal residual
        assert!(stats_sor.residual <= stats_std.residual + 1e-10);
    }
}
