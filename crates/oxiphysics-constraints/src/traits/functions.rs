//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

use super::types::{ConstraintKind, ConstraintPriority, ErrorNorm, SolverHints};

/// Trait for physics constraints that restrict relative motion between bodies.
///
/// Constraints are solved in three phases:
/// 1. **Prepare**: precompute Jacobians, effective masses, and bias terms.
/// 2. **Solve velocity**: apply impulses to correct velocity errors.
/// 3. **Solve position**: apply pseudo-impulses to correct positional drift.
///
/// # Thread Safety
///
/// The trait is *not* `Send + Sync` by default; solvers that parallelize must
/// ensure proper synchronization externally.
pub trait Constraint {
    /// Prepare cached data for this constraint (Jacobians, effective mass, bias).
    ///
    /// Called once per sub-step, before the velocity-solve iterations.
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64);
    /// Apply velocity-level impulses to satisfy the constraint.
    ///
    /// May be called multiple times per sub-step (once per PGS iteration).
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64);
    /// Apply position-level corrections (Baumgarte stabilization or split impulse).
    ///
    /// Called at most once per sub-step after all velocity iterations.
    fn solve_position(&mut self, bodies: &mut RigidBodySet, dt: f64);
    /// Return the body handles involved in this constraint.
    fn body_handles(&self) -> Vec<BodyHandle>;
    /// Return the number of constrained degrees of freedom.
    ///
    /// Default implementation returns `1`.
    fn dof_count(&self) -> usize {
        1
    }
    /// Returns `true` if the constraint is currently active (not sleeping/disabled).
    ///
    /// Inactive constraints are skipped by the solver.
    fn is_active(&self) -> bool {
        true
    }
    /// Called by the solver after all velocity iterations have converged.
    ///
    /// Can be used to update internal state (e.g. accumulated impulses for
    /// warm-starting the next frame).
    fn post_solve(&mut self) {}
    /// Returns the current residual (velocity-error magnitude) after the last
    /// [`solve_velocity`](Constraint::solve_velocity) call.
    ///
    /// Used by the solver to check global convergence.  Default returns 0.0.
    fn residual(&self) -> f64 {
        0.0
    }
}
/// Extension trait for constraints that drive a DOF toward a target.
///
/// Motor constraints combine a regular constraint with a controller that
/// continuously updates the impulse target (e.g. a revolute joint with a
/// servo-motor set to a target angle or angular velocity).
pub trait MotorConstraintTrait: Constraint {
    /// Set the target position (angle, distance, …) for the motor.
    fn set_target_position(&mut self, pos: f64);
    /// Set the target velocity for the motor.
    fn set_target_velocity(&mut self, vel: f64);
    /// Set the maximum force / torque the motor may exert.
    fn set_max_force(&mut self, force: f64);
    /// Returns the current motor force / torque being applied.
    fn current_force(&self) -> f64;
    /// Returns `true` if the motor has reached (or is within tolerance of)
    /// its target position.
    fn at_target(&self, tolerance: f64) -> bool;
}
/// Extension trait for constraints that enforce hard joint limits.
///
/// A limit constraint becomes active when the joint coordinate hits its lower
/// or upper bound, and applies a one-sided impulse to prevent violation.
pub trait LimitConstraintTrait: Constraint {
    /// Set the lower bound for the joint coordinate.
    fn set_lower_limit(&mut self, limit: f64);
    /// Set the upper bound for the joint coordinate.
    fn set_upper_limit(&mut self, limit: f64);
    /// Returns the current joint coordinate value.
    fn joint_coordinate(&self) -> f64;
    /// Returns `true` if the lower limit is currently active.
    fn lower_limit_active(&self) -> bool;
    /// Returns `true` if the upper limit is currently active.
    fn upper_limit_active(&self) -> bool;
    /// Returns the impulse applied by the lower limit in the last solve step.
    fn lower_limit_impulse(&self) -> f64;
    /// Returns the impulse applied by the upper limit in the last solve step.
    fn upper_limit_impulse(&self) -> f64;
}
/// Extension trait for compliant (soft) constraints following the XPBD formulation.
///
/// Soft constraints have a compliance α (inverse stiffness) and a damping
/// coefficient β.  The XPBD Lagrange multiplier update is:
///
/// ```text
/// Δλ = -(C + α̃ λ) / (∇C^T M⁻¹ ∇C + α̃)
/// α̃  = α / (dt²)
/// ```
pub trait SoftConstraint: Constraint {
    /// Returns the compliance α \[m/N or rad/(N·m)\].
    ///
    /// α → 0  ⇒ rigid constraint.
    /// α → ∞  ⇒ spring-like (infinitely compliant).
    fn compliance(&self) -> f64;
    /// Set the compliance α.
    fn set_compliance(&mut self, alpha: f64);
    /// Returns the positional damping coefficient β \[s/m or s·rad/N·m\].
    fn damping(&self) -> f64;
    /// Set the damping coefficient β.
    fn set_damping(&mut self, beta: f64);
    /// Returns the accumulated XPBD Lagrange multiplier Δλ for this frame.
    fn accumulated_lambda(&self) -> f64;
    /// Reset the accumulated Lagrange multiplier (called at the start of each frame).
    fn reset_lambda(&mut self);
    /// Compute the scaled compliance ᾱ = α / dt².
    fn scaled_compliance(&self, dt: f64) -> f64 {
        if dt.abs() < 1e-15 {
            0.0
        } else {
            self.compliance() / (dt * dt)
        }
    }
    /// Compute the XPBD Lagrange multiplier update Δλ.
    ///
    /// * `c`          — constraint residual C(x).
    /// * `jt_m_inv_j` — scalar J M⁻¹ Jᵀ (generalized inverse mass).
    /// * `lambda`     — current accumulated Lagrange multiplier.
    /// * `dt`         — time step.
    fn xpbd_delta_lambda(&self, c: f64, jt_m_inv_j: f64, lambda: f64, dt: f64) -> f64 {
        let alpha_tilde = self.scaled_compliance(dt);
        let denom = jt_m_inv_j + alpha_tilde;
        if denom.abs() < 1e-15 {
            0.0
        } else {
            -(c + alpha_tilde * lambda) / denom
        }
    }
}
/// Observer trait for solver introspection and profiling.
///
/// Implement this to collect per-constraint statistics during the solve.
pub trait ConstraintDiagnostics {
    /// Called at the start of `prepare`.
    fn on_prepare_begin(&self) {}
    /// Called at the end of `prepare`.
    fn on_prepare_end(&self) {}
    /// Called after each velocity-solve iteration.
    ///
    /// * `iteration` — 0-based iteration index.
    /// * `impulse`   — magnitude of the impulse applied in this iteration.
    fn on_velocity_iteration(&self, _iteration: usize, _impulse: f64) {}
    /// Called after position correction.
    ///
    /// * `correction` — magnitude of the position correction applied.
    fn on_position_correction(&self, _correction: f64) {}
    /// Called after `post_solve`.
    fn on_post_solve(&self) {}
}
/// Object-safe clone helper for `Box<dyn Constraint>`.
///
/// Implement this on concrete types to enable boxed-constraint cloning.
pub trait ConstraintClone {
    /// Return a heap-allocated copy of this constraint.
    fn clone_box(&self) -> Box<dyn Constraint>;
}
/// Extension trait providing priority information for a constraint.
pub trait PrioritizedConstraint: Constraint {
    /// Return the solve priority for this constraint.
    fn priority(&self) -> ConstraintPriority {
        ConstraintPriority::Normal
    }
}
/// Extension trait for constraints that support warm-starting.
///
/// Warm-starting seeds the accumulated impulse with a fraction of the previous
/// frame's impulse, dramatically improving solver convergence.
pub trait WarmStartable: Constraint {
    /// Apply warm-start impulses at the beginning of the solve.
    ///
    /// * `bodies` — mutable rigid-body set.
    /// * `warm_start_factor` — fraction of the cached impulse to apply (0–1).
    fn warm_start(&mut self, bodies: &mut RigidBodySet, warm_start_factor: f64);
    /// Save the current accumulated impulse for use in the next frame.
    fn save_impulse(&mut self);
    /// Return the cached impulse from the previous frame.
    fn cached_impulse(&self) -> f64;
}
/// Extension trait for contact constraints that expose collision geometry.
pub trait ContactPatch: Constraint {
    /// World-space contact normal (pointing from body B toward body A).
    fn contact_normal(&self) -> [f64; 3];
    /// Signed penetration depth (positive = overlapping).
    fn penetration_depth(&self) -> f64;
    /// World-space position of the contact point.
    fn contact_point(&self) -> [f64; 3];
    /// Coefficient of restitution (0 = inelastic, 1 = elastic).
    fn restitution(&self) -> f64;
    /// Returns the normal impulse applied in the last solve step.
    fn normal_impulse(&self) -> f64;
    /// Returns the tangential (friction) impulse magnitude in the last solve step.
    fn friction_impulse(&self) -> f64;
}
/// Extension trait that allows a constraint to self-report its kind.
pub trait TypedConstraint: Constraint {
    /// Return the [`ConstraintKind`] for this constraint.
    fn kind(&self) -> ConstraintKind;
}
/// Extension trait for constraints that can break when force/torque exceeds a
/// threshold.
///
/// Once broken, the constraint becomes permanently inactive and is skipped by
/// the solver.  Typical use: destructible joints, fracture simulation.
pub trait BreakableConstraint: Constraint {
    /// Maximum force magnitude before the constraint breaks \[N\].
    fn break_force(&self) -> f64;
    /// Set the break force threshold.
    fn set_break_force(&mut self, force: f64);
    /// Maximum torque magnitude before the constraint breaks \[N·m\].
    fn break_torque(&self) -> f64;
    /// Set the break torque threshold.
    fn set_break_torque(&mut self, torque: f64);
    /// Returns `true` if the constraint has permanently broken.
    fn is_broken(&self) -> bool;
    /// Force the constraint into the broken state (e.g. scripted destruction).
    fn break_now(&mut self);
    /// Check current impulse against thresholds and break if exceeded.
    ///
    /// * `force_impulse`  — linear impulse magnitude applied this step \[N·s\].
    /// * `torque_impulse` — angular impulse magnitude applied this step \[N·m·s\].
    /// * `dt`             — time step \[s\].
    ///
    /// Returns `true` if the constraint just broke.
    fn check_and_break(&mut self, force_impulse: f64, torque_impulse: f64, dt: f64) -> bool {
        if dt <= 0.0 {
            return false;
        }
        let force = force_impulse / dt;
        let torque = torque_impulse / dt;
        if force > self.break_force() || torque > self.break_torque() {
            self.break_now();
            true
        } else {
            false
        }
    }
}
/// Extension trait for constraints that can provide solver hints.
pub trait HintedConstraint: Constraint {
    /// Return solver hints for this constraint.
    fn solver_hints(&self) -> SolverHints {
        SolverHints::default()
    }
}
/// Compute the error norm for a slice of residuals.
pub fn compute_error_norm(residuals: &[f64], norm: ErrorNorm) -> f64 {
    if residuals.is_empty() {
        return 0.0;
    }
    match norm {
        ErrorNorm::L1 => residuals.iter().map(|r| r.abs()).sum(),
        ErrorNorm::L2 => {
            let sum_sq: f64 = residuals.iter().map(|r| r * r).sum();
            sum_sq.sqrt()
        }
        ErrorNorm::LInfinity => residuals.iter().map(|r| r.abs()).fold(0.0_f64, f64::max),
        ErrorNorm::WeightedL2 => {
            let sum_sq: f64 = residuals.iter().map(|r| r * r).sum();
            sum_sq.sqrt()
        }
    }
}
/// Weighted L2 norm (each component multiplied by the corresponding weight).
///
/// Panics in debug mode if `residuals.len() != weights.len()`.
pub fn compute_weighted_l2_norm(residuals: &[f64], weights: &[f64]) -> f64 {
    debug_assert_eq!(
        residuals.len(),
        weights.len(),
        "residuals and weights must have the same length"
    );
    let sum_sq: f64 = residuals
        .iter()
        .zip(weights.iter())
        .map(|(r, w)| (r * w) * (r * w))
        .sum();
    sum_sq.sqrt()
}
/// Normalize a residual vector to have unit L2 norm.
/// Returns a zero vector if the input is already zero.
pub fn normalize_residuals(residuals: &[f64]) -> Vec<f64> {
    let l2 = compute_error_norm(residuals, ErrorNorm::L2);
    if l2 < 1e-15 {
        return vec![0.0; residuals.len()];
    }
    residuals.iter().map(|r| r / l2).collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::types::ComplianceRange;
    use crate::traits::types::ConstraintEvent;
    use crate::traits::types::ConstraintGroup;
    use crate::traits::types::ConstraintMetrics;
    use crate::traits::types::ConstraintState;
    use crate::traits::types::ImpulseCache;
    #[test]
    fn test_constraint_priority_ordering() {
        assert!(ConstraintPriority::Low < ConstraintPriority::Normal);
        assert!(ConstraintPriority::Normal < ConstraintPriority::High);
        assert!(ConstraintPriority::High < ConstraintPriority::Critical);
    }
    #[test]
    fn test_constraint_priority_default() {
        assert_eq!(ConstraintPriority::default(), ConstraintPriority::Normal);
    }
    #[test]
    fn test_constraint_kind_display_fixed() {
        assert_eq!(ConstraintKind::Fixed.to_string(), "Fixed");
    }
    #[test]
    fn test_constraint_kind_display_custom() {
        assert_eq!(ConstraintKind::Custom(7).to_string(), "Custom(7)");
    }
    #[test]
    fn test_constraint_kind_eq() {
        assert_eq!(ConstraintKind::Ball, ConstraintKind::Ball);
        assert_ne!(ConstraintKind::Ball, ConstraintKind::Revolute);
        assert_eq!(ConstraintKind::Custom(1), ConstraintKind::Custom(1));
        assert_ne!(ConstraintKind::Custom(1), ConstraintKind::Custom(2));
    }
    pub(super) struct MockSoftConstraint {
        compliance: f64,
        damping: f64,
        lambda: f64,
    }
    impl MockSoftConstraint {
        pub(super) fn new(compliance: f64, damping: f64) -> Self {
            Self {
                compliance,
                damping,
                lambda: 0.0,
            }
        }
    }
    impl Constraint for MockSoftConstraint {
        fn prepare(&mut self, _bodies: &RigidBodySet, _dt: f64) {}
        fn solve_velocity(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
        fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
        fn body_handles(&self) -> Vec<BodyHandle> {
            vec![]
        }
    }
    impl SoftConstraint for MockSoftConstraint {
        fn compliance(&self) -> f64 {
            self.compliance
        }
        fn set_compliance(&mut self, alpha: f64) {
            self.compliance = alpha;
        }
        fn damping(&self) -> f64 {
            self.damping
        }
        fn set_damping(&mut self, beta: f64) {
            self.damping = beta;
        }
        fn accumulated_lambda(&self) -> f64 {
            self.lambda
        }
        fn reset_lambda(&mut self) {
            self.lambda = 0.0;
        }
    }
    #[test]
    fn test_soft_constraint_scaled_compliance() {
        let sc = MockSoftConstraint::new(1e-4, 0.01);
        let dt = 0.01;
        let alpha_tilde = sc.scaled_compliance(dt);
        assert!(
            (alpha_tilde - 1.0).abs() < 1e-10,
            "alpha_tilde={alpha_tilde}"
        );
    }
    #[test]
    fn test_soft_constraint_scaled_compliance_zero_dt() {
        let sc = MockSoftConstraint::new(1e-4, 0.01);
        let alpha_tilde = sc.scaled_compliance(0.0);
        assert!(alpha_tilde.abs() < 1e-10);
    }
    #[test]
    fn test_xpbd_delta_lambda_rigid_limit() {
        let sc = MockSoftConstraint::new(0.0, 0.0);
        let delta_lambda = sc.xpbd_delta_lambda(0.1, 2.0, 0.0, 0.01);
        assert!(
            (delta_lambda + 0.05).abs() < 1e-10,
            "delta_lambda={delta_lambda}"
        );
    }
    #[test]
    fn test_xpbd_delta_lambda_soft() {
        let sc = MockSoftConstraint::new(1e-4, 0.0);
        let dt = 0.01;
        let alpha_tilde = sc.scaled_compliance(dt);
        let c = 0.1;
        let jt_m_inv_j = 2.0;
        let lambda = 0.0;
        let delta_lambda = sc.xpbd_delta_lambda(c, jt_m_inv_j, lambda, dt);
        let expected = -(c + alpha_tilde * lambda) / (jt_m_inv_j + alpha_tilde);
        assert!(
            (delta_lambda - expected).abs() < 1e-10,
            "delta_lambda={delta_lambda}"
        );
    }
    #[test]
    fn test_xpbd_delta_lambda_degenerate_mass() {
        let sc = MockSoftConstraint::new(0.0, 0.0);
        let delta_lambda = sc.xpbd_delta_lambda(1.0, 0.0, 0.0, 0.01);
        assert!(delta_lambda.abs() < 1e-10, "delta_lambda={delta_lambda}");
    }
    #[test]
    fn test_constraint_default_dof_count() {
        let sc = MockSoftConstraint::new(1e-4, 0.0);
        assert_eq!(sc.dof_count(), 1);
    }
    #[test]
    fn test_constraint_default_is_active() {
        let sc = MockSoftConstraint::new(1e-4, 0.0);
        assert!(sc.is_active());
    }
    #[test]
    fn test_constraint_default_residual() {
        let sc = MockSoftConstraint::new(1e-4, 0.0);
        assert_eq!(sc.residual(), 0.0);
    }
    #[test]
    fn test_constraint_state_default() {
        let s = ConstraintState::default();
        assert_eq!(s, ConstraintState::Uninitialised);
    }
    #[test]
    fn test_constraint_state_display_active() {
        assert_eq!(ConstraintState::Active.to_string(), "Active");
    }
    #[test]
    fn test_constraint_state_display_disabled() {
        assert_eq!(ConstraintState::Disabled.to_string(), "Disabled");
    }
    #[test]
    fn test_constraint_state_display_removed() {
        assert_eq!(ConstraintState::Removed.to_string(), "Removed");
    }
    #[test]
    fn test_constraint_state_display_sleeping() {
        assert_eq!(ConstraintState::Sleeping.to_string(), "Sleeping");
    }
    #[test]
    fn test_constraint_state_display_uninitialised() {
        assert_eq!(ConstraintState::Uninitialised.to_string(), "Uninitialised");
    }
    #[test]
    fn test_constraint_state_eq() {
        assert_eq!(ConstraintState::Active, ConstraintState::Active);
        assert_ne!(ConstraintState::Active, ConstraintState::Sleeping);
    }
    #[test]
    fn test_constraint_metrics_default() {
        let m = ConstraintMetrics::default();
        assert_eq!(m.velocity_iterations, 0);
        assert_eq!(m.velocity_residual, 0.0);
        assert!(m.converged);
    }
    #[test]
    fn test_constraint_metrics_new() {
        let m = ConstraintMetrics::new();
        assert_eq!(m.velocity_residual, 0.0);
    }
    #[test]
    fn test_constraint_metrics_within_tolerance() {
        let m = ConstraintMetrics {
            velocity_residual: 1e-7,
            ..Default::default()
        };
        assert!(m.is_within_tolerance(1e-6));
        assert!(!m.is_within_tolerance(1e-8));
    }
    #[test]
    fn test_constraint_metrics_merge_accumulates_iterations() {
        let mut a = ConstraintMetrics {
            velocity_iterations: 5,
            velocity_residual: 0.01,
            position_correction: 0.001,
            accumulated_impulse: 1.0,
            converged: true,
        };
        let b = ConstraintMetrics {
            velocity_iterations: 3,
            velocity_residual: 0.02,
            position_correction: 0.005,
            accumulated_impulse: 2.0,
            converged: true,
        };
        a.merge(&b);
        assert_eq!(a.velocity_iterations, 8);
        assert!((a.velocity_residual - 0.02).abs() < 1e-15);
        assert!((a.accumulated_impulse - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_constraint_metrics_merge_converged_flag() {
        let mut a = ConstraintMetrics {
            converged: true,
            ..ConstraintMetrics::default()
        };
        let b = ConstraintMetrics {
            converged: false,
            ..ConstraintMetrics::default()
        };
        a.merge(&b);
        assert!(!a.converged);
    }
    #[test]
    fn test_constraint_metrics_merge_position_correction_max() {
        let mut a = ConstraintMetrics {
            position_correction: 0.1,
            ..ConstraintMetrics::default()
        };
        let b = ConstraintMetrics {
            position_correction: 0.5,
            ..ConstraintMetrics::default()
        };
        a.merge(&b);
        assert!((a.position_correction - 0.5).abs() < 1e-15);
    }
    #[test]
    fn test_solver_hints_default() {
        let h = SolverHints::default();
        assert!(h.max_velocity_iterations.is_none());
        assert!(h.needs_position_solve);
        assert!(!h.is_unilateral);
        assert!(h.priority_override.is_none());
    }
    #[test]
    fn test_solver_hints_custom() {
        let h = SolverHints {
            max_velocity_iterations: Some(10),
            needs_position_solve: false,
            is_unilateral: true,
            priority_override: Some(ConstraintPriority::Critical),
        };
        assert_eq!(h.max_velocity_iterations, Some(10));
        assert!(!h.needs_position_solve);
        assert!(h.is_unilateral);
        assert_eq!(h.priority_override, Some(ConstraintPriority::Critical));
    }
    #[test]
    fn test_constraint_group_new_empty() {
        let g = ConstraintGroup::new("structural", ConstraintPriority::High);
        assert_eq!(g.name, "structural");
        assert_eq!(g.priority, ConstraintPriority::High);
        assert!(g.is_empty());
        assert_eq!(g.len(), 0);
    }
    #[test]
    fn test_constraint_group_push() {
        let mut g = ConstraintGroup::new("contacts", ConstraintPriority::Normal);
        g.push(0);
        g.push(5);
        g.push(12);
        assert_eq!(g.len(), 3);
        assert!(!g.is_empty());
        assert_eq!(g.indices, vec![0, 5, 12]);
    }
    #[test]
    fn test_constraint_event_broken_display() {
        let e = ConstraintEvent::Broken {
            constraint_index: 3,
            force: 1000.0,
        };
        let s = e.to_string();
        assert!(s.contains("Broken"), "s={s}");
        assert!(s.contains("3"), "s={s}");
    }
    #[test]
    fn test_constraint_event_contact_begin_display() {
        let e = ConstraintEvent::ContactBegin {
            constraint_index: 7,
        };
        let s = e.to_string();
        assert!(s.contains("ContactBegin"), "s={s}");
        assert!(s.contains("7"), "s={s}");
    }
    #[test]
    fn test_constraint_event_contact_end_display() {
        let e = ConstraintEvent::ContactEnd {
            constraint_index: 2,
        };
        let s = e.to_string();
        assert!(s.contains("ContactEnd"), "s={s}");
    }
    #[test]
    fn test_constraint_event_limit_hit_display() {
        let e = ConstraintEvent::LimitHit {
            constraint_index: 1,
            is_lower: true,
        };
        let s = e.to_string();
        assert!(s.contains("LimitHit"), "s={s}");
        assert!(s.contains("true"), "s={s}");
    }
    #[test]
    fn test_constraint_event_convergence_warning_display() {
        let e = ConstraintEvent::ConvergenceWarning {
            constraint_index: 9,
            residual: 1e-3,
        };
        let s = e.to_string();
        assert!(s.contains("ConvergenceWarning"), "s={s}");
        assert!(s.contains("9"), "s={s}");
    }
    #[test]
    fn test_constraint_event_eq() {
        let a = ConstraintEvent::ContactBegin {
            constraint_index: 5,
        };
        let b = ConstraintEvent::ContactBegin {
            constraint_index: 5,
        };
        let c = ConstraintEvent::ContactEnd {
            constraint_index: 5,
        };
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
    #[test]
    fn test_impulse_cache_new_empty() {
        let cache = ImpulseCache::new(16);
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
    }
    #[test]
    fn test_impulse_cache_store_and_get() {
        let mut cache = ImpulseCache::new(8);
        cache.store(0, 1.5);
        cache.store(3, -0.7);
        assert_eq!(cache.get(0), Some(1.5));
        assert_eq!(cache.get(3), Some(-0.7));
        assert_eq!(cache.get(99), None);
    }
    #[test]
    fn test_impulse_cache_update_existing() {
        let mut cache = ImpulseCache::new(8);
        cache.store(1, 1.0);
        cache.store(1, 2.0);
        assert_eq!(cache.get(1), Some(2.0));
        assert_eq!(cache.len(), 1);
    }
    #[test]
    fn test_impulse_cache_eviction() {
        let mut cache = ImpulseCache::new(3);
        cache.store(0, 1.0);
        cache.store(1, 2.0);
        cache.store(2, 3.0);
        cache.store(3, 4.0);
        assert_eq!(cache.len(), 3);
        assert_eq!(cache.get(0), None);
        assert_eq!(cache.get(3), Some(4.0));
    }
    #[test]
    fn test_impulse_cache_clear() {
        let mut cache = ImpulseCache::new(8);
        cache.store(0, 1.0);
        cache.store(1, 2.0);
        cache.clear();
        assert!(cache.is_empty());
        assert_eq!(cache.get(0), None);
    }
    #[test]
    fn test_compliance_range_clamp_within() {
        let r = ComplianceRange::new(1e-6, 1e-2);
        assert!((r.clamp(1e-4) - 1e-4).abs() < 1e-20);
    }
    #[test]
    fn test_compliance_range_clamp_below_min() {
        let r = ComplianceRange::new(1e-6, 1e-2);
        assert!((r.clamp(1e-10) - 1e-6).abs() < 1e-20);
    }
    #[test]
    fn test_compliance_range_clamp_above_max() {
        let r = ComplianceRange::new(1e-6, 1e-2);
        assert!((r.clamp(1.0) - 1e-2).abs() < 1e-20);
    }
    #[test]
    fn test_compliance_range_lerp_zero() {
        let r = ComplianceRange::new(1e-6, 1e-2);
        assert!((r.lerp(0.0) - 1e-6).abs() < 1e-20);
    }
    #[test]
    fn test_compliance_range_lerp_one() {
        let r = ComplianceRange::new(1e-6, 1e-2);
        assert!((r.lerp(1.0) - 1e-2).abs() < 1e-20);
    }
    #[test]
    fn test_compliance_range_lerp_half() {
        let _r = ComplianceRange::new(0.0_f64.max(f64::EPSILON), 1.0);
        let r2 = ComplianceRange::new(1e-3, 1e-1);
        let mid = r2.lerp(0.5);
        assert!((mid - (1e-3 + 1e-1) * 0.5).abs() < 1e-15, "mid={mid}");
    }
    #[test]
    fn test_compliance_range_midpoint() {
        let r = ComplianceRange::new(2.0, 4.0);
        assert!((r.midpoint() - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_compliance_range_lerp_clamps_t() {
        let r = ComplianceRange::new(1.0, 10.0);
        assert!((r.lerp(2.0) - 10.0).abs() < 1e-15);
        assert!((r.lerp(-1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_constraint_kind_display_all() {
        let kinds = [
            (ConstraintKind::Fixed, "Fixed"),
            (ConstraintKind::Ball, "Ball"),
            (ConstraintKind::Revolute, "Revolute"),
            (ConstraintKind::Prismatic, "Prismatic"),
            (ConstraintKind::Spring, "Spring"),
            (ConstraintKind::Contact, "Contact"),
            (ConstraintKind::SixDof, "SixDof"),
            (ConstraintKind::Gear, "Gear"),
            (ConstraintKind::Pulley, "Pulley"),
            (ConstraintKind::RackPinion, "RackPinion"),
            (ConstraintKind::Motor, "Motor"),
            (ConstraintKind::Pbd, "Pbd"),
        ];
        for (kind, expected) in &kinds {
            assert_eq!(kind.to_string(), *expected, "kind={kind}");
        }
    }
    #[test]
    fn test_constraint_kind_custom_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(ConstraintKind::Custom(1));
        set.insert(ConstraintKind::Custom(2));
        set.insert(ConstraintKind::Custom(1));
        assert_eq!(set.len(), 2);
    }
}
