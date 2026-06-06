//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_impulse_new {
    use super::super::functions::*;
    use crate::impulse::AngularImpulseState;
    use crate::impulse::CollisionManifold;
    use crate::impulse::ContactConstraint;
    use crate::impulse::ImpulseCache;
    use crate::impulse::ImpulseSolver;
    use crate::impulse::PositionalCorrector;
    use crate::impulse::RigidBodyState;
    use crate::impulse::SequentialImpulseSolverConfig;
    use crate::impulse::SubStepRestitutionCorrector;
    use crate::impulse::WarmStartSolver;
    fn sphere(pos: [f64; 3], vel: [f64; 3], inv_mass: f64, inv_inertia: f64) -> RigidBodyState {
        RigidBodyState {
            position: pos,
            velocity: vel,
            angular_velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            inv_mass,
            inv_inertia_local: [inv_inertia; 3],
        }
    }
    #[test]
    fn test_positional_corrector_separates_overlapping_bodies() {
        let mut a = sphere([0.0; 3], [0.0; 3], 1.0, 0.0);
        let mut b = sphere([0.1, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let corrector = PositionalCorrector::new(0.5, 0.0);
        let initial_dist = b.position[0] - a.position[0];
        corrector.correct(&mut a, &mut b, [-1.0, 0.0, 0.0], 0.05);
        let final_dist = b.position[0] - a.position[0];
        assert!(
            final_dist > initial_dist,
            "bodies should be further apart after correction"
        );
    }
    #[test]
    fn test_positional_corrector_slop_prevents_correction() {
        let mut a = sphere([0.0; 3], [0.0; 3], 1.0, 0.0);
        let mut b = sphere([0.001, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let corrector = PositionalCorrector::new(0.5, 0.01);
        let before = a.position[0];
        corrector.correct(&mut a, &mut b, [1.0, 0.0, 0.0], 0.001);
        assert!(
            (a.position[0] - before).abs() < 1e-15,
            "no correction should occur below slop"
        );
    }
    #[test]
    fn test_positional_corrector_static_body_unaffected() {
        let mut a = sphere([0.0; 3], [0.0; 3], 1.0, 0.0);
        let mut b = sphere([0.05, 0.0, 0.0], [0.0; 3], 0.0, 0.0);
        let corrector = PositionalCorrector::new(0.5, 0.0);
        let b_before = b.position[0];
        corrector.correct(&mut a, &mut b, [1.0, 0.0, 0.0], 0.05);
        assert!(
            (b.position[0] - b_before).abs() < 1e-15,
            "static body should not move"
        );
    }
    #[test]
    fn test_angular_impulse_state_torque_integration() {
        let body = sphere([0.0; 3], [0.0; 3], 1.0, 1.0);
        let mut state = AngularImpulseState::new(body);
        state.apply_torque([0.0, 0.0, 10.0]);
        state.integrate_torque(0.1);
        assert!(
            (state.body.angular_velocity[2] - 1.0).abs() < 1e-12,
            "ω_z = {}",
            state.body.angular_velocity[2]
        );
    }
    #[test]
    fn test_angular_impulse_state_torque_reset_after_integrate() {
        let body = sphere([0.0; 3], [0.0; 3], 1.0, 1.0);
        let mut state = AngularImpulseState::new(body);
        state.apply_torque([5.0, 0.0, 0.0]);
        state.integrate_torque(0.1);
        assert_eq!(state.torque, [0.0; 3], "torque should be reset");
    }
    #[test]
    fn test_angular_impulse_state_direct_impulse() {
        let body = sphere([0.0; 3], [0.0; 3], 1.0, 1.0);
        let mut state = AngularImpulseState::new(body);
        state.apply_angular_impulse([0.0, 2.0, 0.0]);
        assert!(
            (state.body.angular_velocity[1] - 2.0).abs() < 1e-12,
            "ω_y = {}",
            state.body.angular_velocity[1]
        );
    }
    #[test]
    fn test_angular_kinetic_energy() {
        let mut body = sphere([0.0; 3], [0.0; 3], 1.0, 1.0);
        body.angular_velocity = [0.0, 0.0, 2.0];
        let state = AngularImpulseState::new(body);
        let ke = state.angular_kinetic_energy();
        assert!((ke - 2.0).abs() < 1e-10, "angular KE = {}", ke);
    }
    #[test]
    fn test_substep_restitution_per_substep_e() {
        let corrector = SubStepRestitutionCorrector::new(0.5, 0.01, 1);
        assert!(
            (corrector.per_substep_restitution() - 0.5).abs() < 1e-10,
            "e_sub = {}",
            corrector.per_substep_restitution()
        );
    }
    #[test]
    fn test_substep_restitution_multi_substep_less_than_full() {
        let corrector = SubStepRestitutionCorrector::new(0.8, 0.01, 4);
        let e_sub = corrector.per_substep_restitution();
        assert!(
            e_sub < 0.8,
            "per-substep e={} should be < full e=0.8",
            e_sub
        );
        assert!(e_sub > 0.0, "per-substep e should be positive");
    }
    #[test]
    fn test_substep_restitution_below_threshold_returns_zero() {
        let corrector = SubStepRestitutionCorrector::new(0.9, 1.0, 1);
        let a = sphere([0.0; 3], [0.5, 0.0, 0.0], 1.0, 0.0);
        let b = sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let j = corrector.compute_substep_impulse(&a, &b, &manifold);
        assert!(
            j.abs() < 1e-15,
            "below threshold, impulse should be 0 but got {}",
            j
        );
    }
    #[test]
    fn test_substep_restitution_above_threshold_returns_nonzero() {
        let corrector = SubStepRestitutionCorrector::new(0.5, 0.1, 1);
        let a = sphere([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.0);
        let b = sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let j = corrector.compute_substep_impulse(&a, &b, &manifold);
        assert!(
            j > 0.0,
            "fast approach should produce non-zero impulse, got {}",
            j
        );
    }
    #[test]
    fn test_impulse_cache_scale() {
        let mut entry = ImpulseCache::new(0, 1);
        entry.lambda_n = 4.0;
        entry.lambda_t1 = 2.0;
        entry.scale(0.5);
        assert!((entry.lambda_n - 2.0).abs() < 1e-15);
        assert!((entry.lambda_t1 - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_warm_start_solver_basic_convergence() {
        let mut bodies = vec![
            sphere([0.0; 3], [3.0, 0.0, 0.0], 1.0, 0.0),
            sphere([2.0, 0.0, 0.0], [-3.0, 0.0, 0.0], 1.0, 0.0),
        ];
        let mut constraints = vec![ContactConstraint::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.01,
            0.5,
            0.0,
        )];
        let config = SequentialImpulseSolverConfig {
            velocity_iterations: 15,
            baumgarte: 0.2,
            dt: 1.0 / 60.0,
            warm_start: true,
        };
        let mut solver = WarmStartSolver::new(config);
        solver.solve(&mut bodies, &mut constraints);
        let v_rel_n = bodies[0].velocity[0] - bodies[1].velocity[0];
        assert!(v_rel_n <= 0.1 + 1e-6, "not resolved: v_rel_n = {}", v_rel_n);
    }
    #[test]
    fn test_warm_start_solver_cache_populated() {
        let mut bodies = vec![
            sphere([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0),
            sphere([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0),
        ];
        let mut constraints = vec![ContactConstraint::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.0,
            0.3,
            0.0,
        )];
        let config = SequentialImpulseSolverConfig::default();
        let mut solver = WarmStartSolver::new(config);
        solver.solve(&mut bodies, &mut constraints);
        assert_eq!(solver.cache.len(), 1, "cache should have one entry");
        assert!(
            solver.cache[0].lambda_n >= 0.0,
            "normal impulse must be non-negative"
        );
    }
    #[test]
    fn test_warm_start_solver_cache_prune() {
        let config = SequentialImpulseSolverConfig::default();
        let mut solver = WarmStartSolver::new(config);
        solver.cache.push(ImpulseCache::new(0, 1));
        solver.cache.push(ImpulseCache::new(2, 3));
        solver.prune_cache(&[(0, 1)]);
        assert_eq!(solver.cache.len(), 1);
        assert_eq!(solver.cache[0].key, (0, 1));
    }
    #[test]
    fn test_restitution_threshold_below_threshold_inelastic() {
        let mut a = sphere([0.0; 3], [0.1, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let e_used = restitution_with_threshold(&mut a, &mut b, &manifold, 0.9, 0.5);
        assert_eq!(e_used, 0.0, "slow collision should use e=0");
    }
    #[test]
    fn test_restitution_threshold_above_threshold_bouncy() {
        let mut a = sphere([0.0; 3], [5.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let e_used = restitution_with_threshold(&mut a, &mut b, &manifold, 0.8, 0.5);
        assert!(
            (e_used - 0.8).abs() < 1e-12,
            "fast collision should use e=0.8"
        );
    }
    #[test]
    fn test_restitution_threshold_separating_bodies_no_impulse() {
        let mut a = sphere([0.0; 3], [-2.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let vx_before = a.velocity[0];
        restitution_with_threshold(&mut a, &mut b, &manifold, 0.8, 0.01);
        assert!(
            (a.velocity[0] - vx_before).abs() < 1e-15,
            "separating bodies should not receive impulse"
        );
    }
    #[test]
    fn test_cor_above_threshold_returns_restitution() {
        let solver = ImpulseSolver::new(0.8, 0.3, 0.5);
        assert!((solver.compute_restitution_coefficient(1.0) - 0.8).abs() < 1e-12);
    }
    #[test]
    fn test_cor_below_threshold_returns_zero() {
        let solver = ImpulseSolver::new(0.8, 0.3, 0.5);
        assert_eq!(solver.compute_restitution_coefficient(0.1), 0.0);
    }
    #[test]
    fn test_cor_exactly_at_threshold_returns_restitution() {
        let solver = ImpulseSolver::new(0.6, 0.2, 1.0);
        assert!((solver.compute_restitution_coefficient(1.0) - 0.6).abs() < 1e-12);
    }
    #[test]
    fn test_coulomb_friction_changes_tangential_velocity() {
        let solver = ImpulseSolver::new(0.0, 0.5, 0.0);
        let mut a = sphere([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere([0.0, 1.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [0.0, 1.0, 0.0],
            depth: 0.0,
            contact_point: [0.0, 0.5, 0.0],
        };
        let vx_before = a.velocity[0];
        let vy_b_before = b.velocity[1];
        solver.apply_coulomb_friction(&mut a, &mut b, &manifold, 10.0);
        let vx_changed = (a.velocity[0] - vx_before).abs() > 1e-10;
        let vy_b_changed = (b.velocity[1] - vy_b_before).abs() > 1e-10;
        let any_changed = vx_changed || vy_b_changed || (b.velocity[0].abs() > 1e-10);
        assert!(
            any_changed,
            "friction should change some velocity: a_vx={}, b_vx={}",
            a.velocity[0], b.velocity[0]
        );
    }
    #[test]
    fn test_coulomb_friction_zero_normal_impulse_no_change() {
        let solver = ImpulseSolver::new(0.0, 0.5, 0.0);
        let mut a = sphere([0.0; 3], [5.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere([0.0, 1.0, 0.0], [0.0; 3], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [0.0, 1.0, 0.0],
            depth: 0.0,
            contact_point: [0.0, 0.5, 0.0],
        };
        let vx_before = a.velocity[0];
        solver.apply_coulomb_friction(&mut a, &mut b, &manifold, 0.0);
        assert!(
            (a.velocity[0] - vx_before).abs() < 1e-12,
            "zero normal impulse should prevent friction"
        );
    }
    #[test]
    fn test_chain_constraint_first_link_receives_impulse() {
        let solver = ImpulseSolver::new(0.0, 0.0, 0.0);
        let mut bodies = vec![
            sphere([0.0; 3], [0.0; 3], 1.0, 0.0),
            sphere([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0),
        ];
        let mags = solver.solve_chain_constraint(&mut bodies, [1.0, 0.0, 0.0], 10.0);
        assert!(
            mags[0] > 0.0,
            "first link should receive impulse: {}",
            mags[0]
        );
        assert!(
            bodies[0].velocity[0] > 0.0,
            "first body vx should be positive"
        );
    }
    #[test]
    fn test_chain_constraint_returns_per_link_magnitudes() {
        let solver = ImpulseSolver::new(0.0, 0.0, 0.0);
        let mut bodies = vec![
            sphere([0.0; 3], [0.0; 3], 1.0, 0.0),
            sphere([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0),
            sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0),
        ];
        let mags = solver.solve_chain_constraint(&mut bodies, [1.0, 0.0, 0.0], 8.0);
        assert_eq!(mags.len(), 3, "should return one magnitude per link");
        assert!(
            mags[1] <= mags[0] + 1e-10,
            "second link should receive <= first: {} vs {}",
            mags[1],
            mags[0]
        );
    }
    #[test]
    fn test_chain_constraint_zero_impulse_no_velocity_change() {
        let solver = ImpulseSolver::new(0.0, 0.0, 0.0);
        let mut bodies = vec![sphere([0.0; 3], [0.0; 3], 1.0, 0.0)];
        solver.solve_chain_constraint(&mut bodies, [1.0, 0.0, 0.0], 0.0);
        assert_eq!(bodies[0].velocity[0], 0.0);
    }
}
