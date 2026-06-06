// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended Position-Based Dynamics (XPBD) integrator.
//!
//! Implements the algorithm from Macklin et al. "XPBD: Position-Based Simulation
//! of Compliant Constrained Dynamics" (2016).  XPBD extends PBD by introducing
//! per-constraint compliance (α) and accumulated Lagrange multipliers (λ),
//! allowing physically meaningful stiffness control without iteration-count
//! dependence.
//!
//! ## Types
//!
//! - `XpbdParticle` — mass point with current and previous position.
//! - `XpbdConstraint` — distance, angle, or volume constraint with compliance.
//! - `XpbdSolver` — hosts particles and constraints; drives the substep loop.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::xpbd::{XpbdConstraint, XpbdSolver};
//!
//! let mut solver = XpbdSolver::new();
//! let a = solver.add_particle([0.0, 0.0, 0.0], 0.0); // pinned
//! let b = solver.add_particle([1.0, 0.0, 0.0], 1.0); // dynamic
//! solver.add_constraint(XpbdConstraint::Distance {
//!     a,
//!     b,
//!     rest_length: 1.0,
//!     compliance: 0.0,
//! });
//! solver.step(1.0 / 60.0);
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Vector math helpers (no external dep)
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A single XPBD particle (mass point).
///
/// Velocity is never stored explicitly; it is derived each substep from the
/// difference `(pos − prev_pos) / h`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XpbdParticle {
    /// Current world-space position.
    pub pos: [f64; 3],
    /// Position at the start of the current substep (used for velocity estimate).
    pub prev_pos: [f64; 3],
    /// Inverse mass.  `0.0` means the particle is static (pinned).
    pub inv_mass: f64,
}

impl XpbdParticle {
    /// Construct a new particle.
    pub fn new(pos: [f64; 3], inv_mass: f64) -> Self {
        Self {
            pos,
            prev_pos: pos,
            inv_mass,
        }
    }
}

/// Constraint between particles in an [`XpbdSolver`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum XpbdConstraint {
    /// Two-particle distance constraint.
    Distance {
        /// Index of first particle.
        a: usize,
        /// Index of second particle.
        b: usize,
        /// Rest (target) distance.
        rest_length: f64,
        /// Compliance α (inverse stiffness).  `0.0` = rigid, `f64::INFINITY` = no-op.
        compliance: f64,
    },
    /// Three-particle angular spring at apex `b`.
    Angle {
        /// First arm endpoint.
        a: usize,
        /// Apex particle.
        b: usize,
        /// Second arm endpoint.
        c: usize,
        /// Rest angle in radians.
        rest_angle: f64,
        /// Compliance α.
        compliance: f64,
    },
    /// Volume conservation for a tetrahedron defined by four particle indices.
    Volume {
        /// Four corner particle indices `[i, j, k, l]`.
        tet: [usize; 4],
        /// Target signed volume (can be used to penalise inversion).
        rest_volume: f64,
        /// Compliance α.
        compliance: f64,
    },
}

/// XPBD solver — owns particles and constraints, drives the substep loop.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XpbdSolver {
    /// All mass points.
    pub particles: Vec<XpbdParticle>,
    /// All constraints.
    pub constraints: Vec<XpbdConstraint>,
    /// Number of substeps taken per [`XpbdSolver::step`] call.
    pub substeps: usize,
    /// Gravity acceleration applied each substep.
    pub gravity: [f64; 3],
}

impl Default for XpbdSolver {
    fn default() -> Self {
        Self {
            particles: Vec::new(),
            constraints: Vec::new(),
            substeps: 10,
            gravity: [0.0, -9.81, 0.0],
        }
    }
}

// ---------------------------------------------------------------------------
// Solver implementation
// ---------------------------------------------------------------------------

impl XpbdSolver {
    /// Create an empty solver with default gravity `(0, −9.81, 0)` and 10 substeps.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a particle and return its index.
    pub fn add_particle(&mut self, pos: [f64; 3], inv_mass: f64) -> usize {
        let idx = self.particles.len();
        self.particles.push(XpbdParticle::new(pos, inv_mass));
        idx
    }

    /// Add a constraint and return its index.
    pub fn add_constraint(&mut self, c: XpbdConstraint) -> usize {
        let idx = self.constraints.len();
        self.constraints.push(c);
        idx
    }

    /// Pin a particle (set `inv_mass` to `0.0` so it never moves).
    pub fn pin_particle(&mut self, idx: usize) {
        if let Some(p) = self.particles.get_mut(idx) {
            p.inv_mass = 0.0;
        }
    }

    /// Estimate the velocity of particle `idx` as `(pos − prev_pos) / dt`.
    ///
    /// Pass the *outer* `dt` (not the substep `h`) to get a meaningful
    /// world-time velocity.  Returns `[0.0; 3]` if `idx` is out of bounds or
    /// `dt ≤ 0`.
    pub fn particle_velocity(&self, idx: usize, dt: f64) -> [f64; 3] {
        if dt <= 0.0 {
            return [0.0; 3];
        }
        match self.particles.get(idx) {
            Some(p) => scale(sub(p.pos, p.prev_pos), 1.0 / dt),
            None => [0.0; 3],
        }
    }

    /// Advance the simulation by one outer time step, using `self.substeps`
    /// internal substeps.
    pub fn step(&mut self, dt: f64) {
        if self.substeps == 0 || dt <= 0.0 {
            return;
        }
        let h = dt / self.substeps as f64;
        let h2 = h * h;
        let grav = self.gravity;

        // Per-constraint Lagrange multiplier accumulator (reset each substep).
        let mut lambdas = vec![0.0_f64; self.constraints.len()];

        for _ in 0..self.substeps {
            // ---- 1. Predict positions (semi-implicit Euler) ----
            for p in &mut self.particles {
                if p.inv_mass == 0.0 {
                    continue;
                }
                // velocity from previous substep
                let v = sub(p.pos, p.prev_pos);
                // v = (pos - prev) / h — but we multiply by h below so divide cancels
                p.prev_pos = p.pos;
                // new_pos = pos + v + h² * gravity
                // where v here is already (pos - prev_pos) i.e. h * velocity
                p.pos = add(add(p.pos, v), scale(grav, h2));
            }

            // Reset lambdas each substep (per XPBD paper).
            lambdas.fill(0.0);

            // ---- 2. Project constraints ----
            for (ci, constraint) in self.constraints.iter().enumerate() {
                match constraint {
                    XpbdConstraint::Distance {
                        a,
                        b,
                        rest_length,
                        compliance,
                    } => {
                        let (ai, bi) = (*a, *b);
                        let (rl, alpha) = (*rest_length, *compliance);

                        // Bounds check — skip silently if invalid indices.
                        if ai >= self.particles.len() || bi >= self.particles.len() {
                            continue;
                        }

                        let pa = self.particles[ai].pos;
                        let pb = self.particles[bi].pos;
                        let ima = self.particles[ai].inv_mass;
                        let imb = self.particles[bi].inv_mass;

                        let delta = sub(pb, pa);
                        let d = len(delta);
                        if d < 1e-12 {
                            continue;
                        }

                        let c_val = d - rl;
                        let alpha_tilde = alpha / h2;

                        // grad_a = -delta/d, grad_b = +delta/d
                        // |grad|² = 1 for both
                        let w_sum = ima + imb + alpha_tilde;
                        if w_sum < 1e-12 {
                            continue;
                        }

                        let delta_lambda = (-c_val - alpha_tilde * lambdas[ci]) / w_sum;
                        lambdas[ci] += delta_lambda;

                        let dir = scale(delta, 1.0 / d);
                        if ima > 0.0 {
                            // correction = inv_mass * delta_lambda * grad_a = inv_mass * (-delta_lambda) * dir
                            self.particles[ai].pos =
                                add(self.particles[ai].pos, scale(dir, -ima * delta_lambda));
                        }
                        if imb > 0.0 {
                            self.particles[bi].pos =
                                add(self.particles[bi].pos, scale(dir, imb * delta_lambda));
                        }
                    }

                    XpbdConstraint::Angle {
                        a,
                        b,
                        c,
                        rest_angle,
                        compliance,
                    } => {
                        let (ai, bi, ci_idx) = (*a, *b, *c);
                        let (rest, alpha) = (*rest_angle, *compliance);

                        if ai >= self.particles.len()
                            || bi >= self.particles.len()
                            || ci_idx >= self.particles.len()
                        {
                            continue;
                        }

                        let pa = self.particles[ai].pos;
                        let pb = self.particles[bi].pos;
                        let pc = self.particles[ci_idx].pos;
                        let ima = self.particles[ai].inv_mass;
                        let imb = self.particles[bi].inv_mass;
                        let imc = self.particles[ci_idx].inv_mass;

                        // Arm vectors from apex b
                        let ba = sub(pa, pb);
                        let bc = sub(pc, pb);

                        let len_ba = len(ba);
                        let len_bc = len(bc);
                        if len_ba < 1e-12 || len_bc < 1e-12 {
                            continue;
                        }

                        // cos θ, clamped to [-1, 1] to avoid NaN from acos
                        let cos_theta = (dot(ba, bc) / (len_ba * len_bc)).clamp(-1.0, 1.0);
                        let theta = cos_theta.acos();
                        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-12);

                        let c_val = theta - rest;
                        let alpha_tilde = alpha / h2;

                        // Gradient of θ w.r.t. each particle position.
                        // d(θ)/d(a) = -grad_cos_a / sin_θ
                        // grad_cos_a = (bc/|bc| - ba*cosθ/|ba|) / |ba|
                        let grad_cos_a = scale(
                            sub(scale(bc, 1.0 / len_bc), scale(ba, cos_theta / len_ba)),
                            1.0 / len_ba,
                        );
                        let grad_cos_c = scale(
                            sub(scale(ba, 1.0 / len_ba), scale(bc, cos_theta / len_bc)),
                            1.0 / len_bc,
                        );

                        let grad_theta_a = scale(grad_cos_a, -1.0 / sin_theta);
                        let grad_theta_c = scale(grad_cos_c, -1.0 / sin_theta);
                        // Conservation: grad_b = -(grad_a + grad_c)
                        let grad_theta_b = scale(add(grad_theta_a, grad_theta_c), -1.0);

                        let w_sum = ima * dot(grad_theta_a, grad_theta_a)
                            + imb * dot(grad_theta_b, grad_theta_b)
                            + imc * dot(grad_theta_c, grad_theta_c)
                            + alpha_tilde;

                        if w_sum < 1e-12 {
                            continue;
                        }

                        let delta_lambda = (-c_val - alpha_tilde * lambdas[ci]) / w_sum;
                        lambdas[ci] += delta_lambda;

                        if ima > 0.0 {
                            self.particles[ai].pos = add(
                                self.particles[ai].pos,
                                scale(grad_theta_a, ima * delta_lambda),
                            );
                        }
                        if imb > 0.0 {
                            self.particles[bi].pos = add(
                                self.particles[bi].pos,
                                scale(grad_theta_b, imb * delta_lambda),
                            );
                        }
                        if imc > 0.0 {
                            self.particles[ci_idx].pos = add(
                                self.particles[ci_idx].pos,
                                scale(grad_theta_c, imc * delta_lambda),
                            );
                        }
                    }

                    XpbdConstraint::Volume {
                        tet,
                        rest_volume,
                        compliance,
                    } => {
                        let [ii, ji, ki, li] = [tet[0], tet[1], tet[2], tet[3]];
                        let (rv, alpha) = (*rest_volume, *compliance);

                        if ii >= self.particles.len()
                            || ji >= self.particles.len()
                            || ki >= self.particles.len()
                            || li >= self.particles.len()
                        {
                            continue;
                        }

                        let pi = self.particles[ii].pos;
                        let pj = self.particles[ji].pos;
                        let pk = self.particles[ki].pos;
                        let pl = self.particles[li].pos;

                        let imi = self.particles[ii].inv_mass;
                        let imj = self.particles[ji].inv_mass;
                        let imk = self.particles[ki].inv_mass;
                        let iml = self.particles[li].inv_mass;

                        // Signed volume V = (1/6) * (j-i) · ((k-i) × (l-i))
                        let ji_vec = sub(pj, pi);
                        let ki_vec = sub(pk, pi);
                        let li_vec = sub(pl, pi);

                        let vol = dot(ji_vec, cross(ki_vec, li_vec)) / 6.0;
                        let c_val = vol - rv;
                        let alpha_tilde = alpha / h2;

                        // Gradients of V w.r.t. each vertex (verified signs):
                        // dV/d(j) = (1/6) * (k-i) × (l-i)
                        // dV/d(k) = (1/6) * (l-i) × (j-i)
                        // dV/d(l) = (1/6) * (j-i) × (k-i)
                        // dV/d(i) = -(dV/d(j) + dV/d(k) + dV/d(l))
                        let grad_j = scale(cross(ki_vec, li_vec), 1.0 / 6.0);
                        let grad_k = scale(cross(li_vec, ji_vec), 1.0 / 6.0);
                        let grad_l = scale(cross(ji_vec, ki_vec), 1.0 / 6.0);
                        let grad_i = scale(add(add(grad_j, grad_k), grad_l), -1.0);

                        let w_sum = imi * dot(grad_i, grad_i)
                            + imj * dot(grad_j, grad_j)
                            + imk * dot(grad_k, grad_k)
                            + iml * dot(grad_l, grad_l)
                            + alpha_tilde;

                        if w_sum < 1e-12 {
                            continue;
                        }

                        let delta_lambda = (-c_val - alpha_tilde * lambdas[ci]) / w_sum;
                        lambdas[ci] += delta_lambda;

                        if imi > 0.0 {
                            self.particles[ii].pos =
                                add(self.particles[ii].pos, scale(grad_i, imi * delta_lambda));
                        }
                        if imj > 0.0 {
                            self.particles[ji].pos =
                                add(self.particles[ji].pos, scale(grad_j, imj * delta_lambda));
                        }
                        if imk > 0.0 {
                            self.particles[ki].pos =
                                add(self.particles[ki].pos, scale(grad_k, imk * delta_lambda));
                        }
                        if iml > 0.0 {
                            self.particles[li].pos =
                                add(self.particles[li].pos, scale(grad_l, iml * delta_lambda));
                        }
                    }
                }
            }
            // Velocity is implicit: v = (pos - prev_pos) / h — no explicit update needed.
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: distance between two particles.
    fn particle_dist(s: &XpbdSolver, a: usize, b: usize) -> f64 {
        len(sub(s.particles[a].pos, s.particles[b].pos))
    }

    // -----------------------------------------------------------------------
    // Test 1: Stiff distance constraint at rest — distance stays 1.0.
    // -----------------------------------------------------------------------
    #[test]
    fn test_distance_at_rest() {
        let mut solver = XpbdSolver::new();
        solver.gravity = [0.0, 0.0, 0.0];
        let a = solver.add_particle([0.0, 0.0, 0.0], 0.0); // pinned
        let b = solver.add_particle([1.0, 0.0, 0.0], 1.0);
        solver.add_constraint(XpbdConstraint::Distance {
            a,
            b,
            rest_length: 1.0,
            compliance: 0.0,
        });

        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            solver.step(dt);
        }

        let d = particle_dist(&solver, a, b);
        assert!(
            (d - 1.0).abs() < 1e-4,
            "Distance should remain 1.0, got {d}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: Compliance=0 rigid — distance 2→1 after one step.
    // -----------------------------------------------------------------------
    #[test]
    fn test_rigid_compliance_corrects_distance() {
        let mut solver = XpbdSolver::new();
        solver.gravity = [0.0, 0.0, 0.0];
        let a = solver.add_particle([0.0, 0.0, 0.0], 1.0);
        let b = solver.add_particle([2.0, 0.0, 0.0], 1.0);
        solver.add_constraint(XpbdConstraint::Distance {
            a,
            b,
            rest_length: 1.0,
            compliance: 0.0,
        });

        solver.step(1.0 / 60.0);

        let d = particle_dist(&solver, a, b);
        assert!(
            (d - 1.0).abs() < 1e-4,
            "Rigid constraint should correct distance to 1.0, got {d}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: Very high compliance (nearly ∞) — positions nearly unchanged.
    // -----------------------------------------------------------------------
    #[test]
    fn test_high_compliance_no_correction() {
        let mut solver = XpbdSolver::new();
        solver.gravity = [0.0, 0.0, 0.0]; // no gravity so positions are stable
        let a = solver.add_particle([0.0, 0.0, 0.0], 1.0);
        let b = solver.add_particle([2.0, 0.0, 0.0], 1.0);
        solver.add_constraint(XpbdConstraint::Distance {
            a,
            b,
            rest_length: 1.0,
            compliance: 1e30, // effectively infinite
        });

        let pos_a_before = solver.particles[a].pos;
        let pos_b_before = solver.particles[b].pos;

        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            solver.step(dt);
        }

        let da = len(sub(solver.particles[a].pos, pos_a_before));
        let db = len(sub(solver.particles[b].pos, pos_b_before));
        assert!(
            da < 0.01,
            "High compliance: particle a should barely move, moved {da}"
        );
        assert!(
            db < 0.01,
            "High compliance: particle b should barely move, moved {db}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 4: Pinned pendulum — distance stays ≈ 1.0 and Y never explodes.
    // -----------------------------------------------------------------------
    #[test]
    fn test_pinned_pendulum_bounded() {
        let mut solver = XpbdSolver::new();
        solver.gravity = [0.0, -9.81, 0.0];
        solver.substeps = 20;

        let pin = solver.add_particle([0.0, 0.0, 0.0], 0.0); // pinned at origin
        let bob = solver.add_particle([1.0, 0.0, 0.0], 1.0); // free

        solver.add_constraint(XpbdConstraint::Distance {
            a: pin,
            b: bob,
            rest_length: 1.0,
            compliance: 0.0,
        });

        let dt = 1.0 / 60.0;
        for _ in 0..120 {
            solver.step(dt);
        }

        let d = particle_dist(&solver, pin, bob);
        assert!(
            (d - 1.0).abs() < 0.1,
            "Pendulum distance should be ≈ 1.0, got {d}"
        );

        let bob_y = solver.particles[bob].pos[1];
        assert!(
            bob_y <= 0.1,
            "Bob Y should not exceed 0.1 (pendulum should swing down), got {bob_y}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 5: Volume tet recovery — inverted tet returns positive volume.
    // -----------------------------------------------------------------------
    #[test]
    fn test_volume_tet_recovery() {
        let mut solver = XpbdSolver::new();
        solver.gravity = [0.0, 0.0, 0.0];
        solver.substeps = 20;

        // Standard unit tet
        let i = solver.add_particle([0.0, 0.0, 0.0], 1.0);
        let j = solver.add_particle([1.0, 0.0, 0.0], 1.0);
        let k = solver.add_particle([0.0, 1.0, 0.0], 1.0);
        let l = solver.add_particle([0.0, 0.0, 1.0], 1.0);

        // Invert: flip j and k  → volume goes negative
        solver.particles[j].pos = [0.0, 1.0, 0.0];
        solver.particles[j].prev_pos = [0.0, 1.0, 0.0];
        solver.particles[k].pos = [1.0, 0.0, 0.0];
        solver.particles[k].prev_pos = [1.0, 0.0, 0.0];

        // Target: positive 1/6 volume
        let rest_vol = 1.0_f64 / 6.0;
        solver.add_constraint(XpbdConstraint::Volume {
            tet: [i, j, k, l],
            rest_volume: rest_vol,
            compliance: 1e-6,
        });

        let dt = 1.0 / 60.0;
        for _ in 0..20 {
            solver.step(dt);
        }

        let pi = solver.particles[i].pos;
        let pj = solver.particles[j].pos;
        let pk = solver.particles[k].pos;
        let pl = solver.particles[l].pos;
        let ji_vec = sub(pj, pi);
        let ki_vec = sub(pk, pi);
        let li_vec = sub(pl, pi);
        let vol = dot(ji_vec, cross(ki_vec, li_vec)) / 6.0;

        assert!(
            vol > 0.0,
            "Tet volume should be positive after recovery, got {vol}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 6: Serde round-trip — particle and constraint counts preserved.
    // -----------------------------------------------------------------------
    #[test]
    fn test_serde_round_trip() {
        let mut solver = XpbdSolver::new();
        let a = solver.add_particle([0.0, 0.0, 0.0], 0.0);
        let b = solver.add_particle([1.0, 0.0, 0.0], 1.0);
        let c = solver.add_particle([0.5, 1.0, 0.0], 1.0);
        solver.add_constraint(XpbdConstraint::Distance {
            a,
            b,
            rest_length: 1.0,
            compliance: 0.0,
        });
        solver.add_constraint(XpbdConstraint::Angle {
            a,
            b,
            c,
            rest_angle: std::f64::consts::FRAC_PI_2,
            compliance: 1e-4,
        });
        solver.add_constraint(XpbdConstraint::Distance {
            a: b,
            b: c,
            rest_length: 1.118,
            compliance: 0.0,
        });

        let json = serde_json::to_string(&solver).expect("serialise");
        let restored: XpbdSolver = serde_json::from_str(&json).expect("deserialise");

        assert_eq!(
            restored.particles.len(),
            solver.particles.len(),
            "Particle count mismatch after serde"
        );
        assert_eq!(
            restored.constraints.len(),
            solver.constraints.len(),
            "Constraint count mismatch after serde"
        );
    }
}
