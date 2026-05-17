// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soft joint limits using linear penalty + one-sided damping forces.
//!
//! Joint limits are enforced by injecting penalty torques/forces into the
//! generalised-force vector `τ` before calling the dynamics algorithms.
//! The penalty law is:
//!
//! ```text
//! penalty(q, q_dot) =
//!   if q < lo:  +k·(lo − q) + c·max(0, −q_dot)
//!   if q > hi:  −k·(q − hi) + c·min(0, −q_dot)
//!   otherwise:  0
//! ```
//!
//! Damping is one-sided: it only opposes motion *into* the limit.

use crate::model::ArticulatedModel;

// ─── JointLimit ──────────────────────────────────────────────────────────────

/// Soft limit for a single degree of freedom.
#[derive(Debug, Clone, Copy)]
pub struct JointLimit {
    /// Lower bound [rad or m].
    pub lo: f64,
    /// Upper bound [rad or m].
    pub hi: f64,
    /// Stiffness coefficient [N·m/rad or N/m].
    pub stiffness: f64,
    /// Damping coefficient [N·m·s/rad or N·s/m].
    pub damping: f64,
}

impl JointLimit {
    /// Compute the penalty generalised force for this DOF.
    ///
    /// Returns `0.0` when `q ∈ [lo, hi]`.
    pub fn penalty(&self, q: f64, q_dot: f64) -> f64 {
        if q < self.lo {
            let penetration = self.lo - q;
            let damp = if q_dot < 0.0 {
                -self.damping * q_dot
            } else {
                0.0
            };
            self.stiffness * penetration + damp
        } else if q > self.hi {
            let penetration = q - self.hi;
            let damp = if q_dot > 0.0 {
                -self.damping * q_dot
            } else {
                0.0
            };
            -self.stiffness * penetration + damp
        } else {
            0.0
        }
    }
}

// ─── JointLimitSet ────────────────────────────────────────────────────────────

/// Per-body soft limits indexed by joint DOF.
///
/// Created via [`JointLimitSet::for_model`] to allocate the correct storage for
/// the given model. Individual limits are set with [`JointLimitSet::set`] and
/// applied to a generalised-force vector with [`JointLimitSet::apply`].
pub struct JointLimitSet {
    /// `limits[body_index][dof_index]` — `None` means no limit for that DOF.
    limits: Vec<Vec<Option<JointLimit>>>,
}

impl JointLimitSet {
    /// Create an empty limit set for the given model.
    ///
    /// All limits are initially `None` (no limit enforced).
    pub fn for_model(model: &ArticulatedModel) -> Self {
        let limits = (0..model.num_bodies())
            .map(|i| vec![None; model.dof_count(i)])
            .collect();
        Self { limits }
    }

    /// Set the soft limit for a specific body and DOF index.
    ///
    /// `dof_idx` is zero-based within body `body`'s joint.
    ///
    /// # Panics
    ///
    /// Panics if `body` or `dof_idx` is out of range.
    pub fn set(&mut self, body: usize, dof_idx: usize, limit: JointLimit) {
        self.limits[body][dof_idx] = Some(limit);
    }

    /// Add penalty forces into `tau` in-place.
    ///
    /// Iterates over all bodies and DOFs, computing the penalty torque/force
    /// for each active limit and adding it to `tau`.
    ///
    /// `q`, `qd`, and `tau` must all have length `model.total_dof()`.
    pub fn apply(&self, model: &ArticulatedModel, q: &[f64], qd: &[f64], tau: &mut [f64]) {
        for body in 0..model.num_bodies() {
            let start = model.dof_start(body);
            for dof in 0..model.dof_count(body) {
                if let Some(ref lim) = self.limits[body][dof] {
                    let idx = start + dof;
                    tau[idx] += lim.penalty(q[idx], qd[idx]);
                }
            }
        }
    }
}

// ─── Convenience: ABA with limits ────────────────────────────────────────────

/// Run ABA forward dynamics after applying soft-limit penalty forces.
///
/// Equivalent to: `limits.apply(model, q, qd, &mut tau_copy); aba(model, q, qd, &tau_copy)`.
///
/// The caller's `tau` slice is NOT modified; a temporary copy is used.
pub fn aba_with_limits(
    model: &ArticulatedModel,
    q: &[f64],
    qd: &[f64],
    tau: &[f64],
    limits: &JointLimitSet,
) -> Vec<f64> {
    let mut tau_with_limits = tau.to_vec();
    limits.apply(model, q, qd, &mut tau_with_limits);
    crate::aba::aba(model, q, qd, &tau_with_limits)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_limit_inside_range_zero_penalty() {
        let lim = JointLimit {
            lo: -1.0,
            hi: 1.0,
            stiffness: 100.0,
            damping: 2.0,
        };
        assert_eq!(lim.penalty(0.0, 0.0), 0.0);
        assert_eq!(lim.penalty(0.5, 1.0), 0.0);
        assert_eq!(lim.penalty(-0.9, -5.0), 0.0);
    }

    #[test]
    fn test_limit_above_high_returns_negative_penalty() {
        let lim = JointLimit {
            lo: -1.0,
            hi: 1.0,
            stiffness: 100.0,
            damping: 2.0,
        };
        let p = lim.penalty(1.1, 0.0);
        // Expected: -100 * 0.1 = -10.0 (no damping since qd=0)
        assert!((p - (-10.0)).abs() < 1e-12, "got {p}");
    }

    #[test]
    fn test_limit_below_low_returns_positive_penalty() {
        let lim = JointLimit {
            lo: -1.0,
            hi: 1.0,
            stiffness: 100.0,
            damping: 2.0,
        };
        let p = lim.penalty(-1.2, 0.0);
        // Expected: +100 * 0.2 = +20.0
        assert!((p - 20.0).abs() < 1e-12, "got {p}");
    }

    #[test]
    fn test_limit_damping_only_when_moving_into_limit() {
        let lim = JointLimit {
            lo: -1.0,
            hi: 1.0,
            stiffness: 100.0,
            damping: 2.0,
        };
        // q > hi, qd > 0 (moving deeper into limit): penalty + damping
        let p_into = lim.penalty(1.1, 0.5);
        // Expected: -10.0 (stiffness) + (-2 * 0.5) = -11.0
        assert!((p_into - (-11.0)).abs() < 1e-12, "got {p_into}");
        // q > hi, qd < 0 (moving back out): only stiffness, no damping
        let p_out = lim.penalty(1.1, -0.5);
        assert!((p_out - (-10.0)).abs() < 1e-12, "got {p_out}");
    }
}
