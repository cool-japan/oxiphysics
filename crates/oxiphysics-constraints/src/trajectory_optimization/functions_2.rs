//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::trajectory_optimization::BSplineTrajectory;
    use crate::trajectory_optimization::BoundaryConditions;
    use crate::trajectory_optimization::ConvergenceMonitor;
    use crate::trajectory_optimization::CostFunction;
    use crate::trajectory_optimization::DynamicLimits;
    use crate::trajectory_optimization::DynamicsModel;
    use crate::trajectory_optimization::KeepOutZone;
    use crate::trajectory_optimization::MultiPhaseTrajectory;
    use crate::trajectory_optimization::Trajectory;
    use crate::trajectory_optimization::TrajectoryKnot;
    use crate::trajectory_optimization::TrajectoryPhase;
    use crate::trajectory_optimization::TrustRegionConfig;
    use crate::trajectory_optimization::WaypointConstraint;
    fn double_integrator(x: &[f64], u: &[f64], _t: f64) -> Vec<f64> {
        vec![x[1], u[0]]
    }
    fn point_mass_2d(x: &[f64], u: &[f64], _t: f64) -> Vec<f64> {
        vec![x[2], x[3], u[0], u[1]]
    }
    fn point_mass_3d(x: &[f64], u: &[f64], _t: f64) -> Vec<f64> {
        vec![x[3], x[4], x[5], u[0], u[1], u[2]]
    }
    fn make_double_integrator() -> DynamicsModel {
        DynamicsModel::new(2, 1, double_integrator)
    }
    fn make_point_mass_2d() -> DynamicsModel {
        DynamicsModel::new(4, 2, point_mass_2d)
    }
    fn make_point_mass_3d() -> DynamicsModel {
        DynamicsModel::new(6, 3, point_mass_3d)
    }
    #[test]
    fn test_dot_product() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot(&a, &b) - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_norm() {
        let v = [3.0, 4.0];
        assert!((norm(&v) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec_operations() {
        let a = vec![1.0, 2.0];
        let b = vec![3.0, 4.0];
        let s = vec_add(&a, &b);
        assert!((s[0] - 4.0).abs() < 1e-10);
        assert!((s[1] - 6.0).abs() < 1e-10);
        let d = vec_sub(&b, &a);
        assert!((d[0] - 2.0).abs() < 1e-10);
        let sc = vec_scale(&a, 3.0);
        assert!((sc[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_solve_linear_2x2() {
        let a = [2.0, 1.0, 1.0, 3.0];
        let b = [5.0, 7.0];
        let x = solve_linear(&a, &b, 2).unwrap();
        assert!((x[0] - 1.6).abs() < 1e-10);
        assert!((x[1] - 1.8).abs() < 1e-10);
    }
    #[test]
    fn test_cost_function_quadratic() {
        let cost = CostFunction::quadratic(2, 1, &[1.0, 1.0], &[0.1]);
        let x = [1.0, 2.0];
        let u = [3.0];
        let rc = cost.running_cost(&x, &u);
        assert!((rc - 5.9).abs() < 1e-10, "rc = {rc}");
    }
    #[test]
    fn test_cost_minimum_energy() {
        let cost = CostFunction::minimum_energy(2, 1);
        let x = [100.0, 200.0];
        let u = [5.0];
        let rc = cost.running_cost(&x, &u);
        assert!((rc - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_cost_minimum_time() {
        let cost = CostFunction::minimum_time(2, 1);
        let x = [100.0, 200.0];
        let u = [5.0];
        assert!((cost.running_cost(&x, &u) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_trajectory_linear_guess() {
        let traj = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[10.0, 0.0], 0.0, 1.0, 11);
        assert_eq!(traj.num_knots(), 11);
        assert!((traj.duration() - 1.0).abs() < 1e-10);
        assert!((traj.knots[5].state[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_keep_out_zone_outside() {
        let zone = KeepOutZone::new([5.0, 0.0, 0.0], 1.0);
        let pos = [0.0, 0.0, 0.0];
        let v = zone.violation(&pos);
        assert!(v > 0.0, "Should be outside: v={v}");
    }
    #[test]
    fn test_keep_out_zone_inside() {
        let zone = KeepOutZone::new([0.0, 0.0, 0.0], 5.0);
        let pos = [1.0, 0.0, 0.0];
        let v = zone.violation(&pos);
        assert!(v < 0.0, "Should be inside: v={v}");
    }
    #[test]
    fn test_keep_out_zone_gradient() {
        let zone = KeepOutZone::new([0.0, 0.0, 0.0], 1.0);
        let pos = [3.0, 4.0, 0.0];
        let g = zone.violation_gradient(&pos);
        assert!((g[0] - 0.6).abs() < 1e-10);
        assert!((g[1] - 0.8).abs() < 1e-10);
    }
    #[test]
    fn test_waypoint_constraint() {
        let wp = WaypointConstraint::full_state(5, vec![10.0, 0.0], 0.01);
        let state = [10.005, 0.003];
        assert!(wp.is_satisfied(&state));
        let state_bad = [11.0, 0.0];
        assert!(!wp.is_satisfied(&state_bad));
    }
    #[test]
    fn test_waypoint_position_only() {
        let wp = WaypointConstraint::position_only(0, [1.0, 2.0, 3.0], 6, 0.1);
        let state = [1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!(wp.is_satisfied(&state));
    }
    #[test]
    fn test_dynamic_limits_feasible() {
        let limits = DynamicLimits::uniform(4, 2, 10.0, 5.0, 3.0);
        assert!(limits.control_feasible(&[2.0, -2.0]));
        assert!(!limits.control_feasible(&[4.0, 0.0]));
    }
    #[test]
    fn test_dynamic_limits_clamp() {
        let limits = DynamicLimits::uniform(2, 2, 10.0, 5.0, 3.0);
        let clamped = limits.clamp_control(&[5.0, -1.0]);
        assert!((clamped[0] - 3.0).abs() < 1e-10);
        assert!((clamped[1] + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_boundary_conditions() {
        let bc = BoundaryConditions::fixed(vec![0.0, 0.0], vec![10.0, 0.0]);
        assert!((bc.initial_violation(&[0.0, 0.0])).abs() < 1e-10);
        assert!(bc.final_violation(&[10.0, 0.0]).abs() < 1e-10);
        assert!(bc.final_violation(&[5.0, 0.0]) > 1.0);
    }
    #[test]
    fn test_dynamics_rk4_double_integrator() {
        let model = make_double_integrator();
        let x = model.rk4_step(&[0.0, 1.0], &[0.0], 0.0, 0.1);
        assert!((x[0] - 0.1).abs() < 1e-10, "pos = {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-10, "vel = {}", x[1]);
    }
    #[test]
    fn test_dynamics_rk4_with_accel() {
        let model = make_double_integrator();
        let x = model.rk4_step(&[0.0, 0.0], &[1.0], 0.0, 1.0);
        assert!((x[0] - 0.5).abs() < 1e-6, "pos = {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-6, "vel = {}", x[1]);
    }
    #[test]
    fn test_dynamics_jacobian() {
        let model = make_double_integrator();
        let jac = model.state_jacobian(&[1.0, 2.0], &[0.0], 0.0, 1e-7);
        assert!(jac[0].abs() < 1e-4);
        assert!((jac[1] - 1.0).abs() < 1e-4);
        assert!(jac[2].abs() < 1e-4);
        assert!(jac[3].abs() < 1e-4);
    }
    #[test]
    fn test_direct_shooting() {
        let model = make_double_integrator();
        let cost = CostFunction::minimum_energy(2, 1);
        let controls: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0]).collect();
        let result = direct_shooting(&model, &[0.0, 0.0], &controls, 0.0, 0.1, &cost);
        assert_eq!(result.trajectory.num_knots(), 11);
        assert!((result.final_state[1] - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_multiple_shooting() {
        let model = make_double_integrator();
        let seg_states = vec![vec![0.0, 0.0], vec![0.5, 1.0], vec![1.5, 1.0]];
        let controls: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0]).collect();
        let (defects, traj) = multiple_shooting(&model, &seg_states, &controls, 0.0, 0.1, 5);
        assert_eq!(defects.len(), 2);
        assert!(traj.num_knots() > 0);
    }
    #[test]
    fn test_trapezoidal_defects_exact() {
        let model = make_double_integrator();
        let mut traj = Trajectory::new(2, 1);
        for i in 0..5 {
            let t = i as f64 * 0.1;
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: vec![t, 1.0],
                control: vec![0.0],
            });
        }
        let defects = trapezoidal_defects(&model, &traj);
        for (k, d) in defects.iter().enumerate() {
            assert!(norm(d) < 1e-10, "defect[{k}] = {:.6e}", norm(d));
        }
    }
    #[test]
    fn test_hermite_simpson_defects_exact() {
        let model = make_double_integrator();
        let mut traj = Trajectory::new(2, 1);
        for i in 0..5 {
            let t = i as f64 * 0.1;
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: vec![t, 1.0],
                control: vec![0.0],
            });
        }
        let defects = hermite_simpson_defects(&model, &traj);
        for (k, d) in defects.iter().enumerate() {
            assert!(norm(d) < 1e-10, "HS defect[{k}] = {:.6e}", norm(d));
        }
    }
    #[test]
    fn test_adjoint_costate() {
        let model = make_double_integrator();
        let cost = CostFunction::quadratic(2, 1, &[1.0, 1.0], &[0.1]);
        let traj = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[1.0, 0.0], 0.0, 1.0, 11);
        let costates = adjoint_costate(&model, &traj, &cost);
        assert_eq!(costates.len(), 11);
        assert!(
            (costates[10][0] - 2.0).abs() < 1e-6,
            "lambda_T = {:?}",
            costates[10]
        );
    }
    #[test]
    fn test_trust_region_accept() {
        let config = TrustRegionConfig::default();
        let result = trust_region_evaluate(&config, 1.0, 10.0, 8.0, 2.5, 0.8);
        assert!(result.accepted);
        assert!((result.actual_reduction - 2.0).abs() < 1e-10);
        assert!((result.ratio - 0.8).abs() < 1e-10);
    }
    #[test]
    fn test_trust_region_reject() {
        let config = TrustRegionConfig::default();
        let result = trust_region_evaluate(&config, 1.0, 10.0, 10.5, 2.5, 0.8);
        assert!(!result.accepted);
        assert!(result.new_radius < 1.0);
    }
    #[test]
    fn test_convergence_monitor() {
        let mut monitor = ConvergenceMonitor::new(1e-6, 1e-6, 100);
        monitor.record(10.0, 1.0, 0.1);
        monitor.record(5.0, 0.1, 0.05);
        assert!(!monitor.is_converged());
        monitor.record(5.0 - 1e-8, 1e-8, 0.01);
        assert!(monitor.is_converged());
        assert_eq!(monitor.num_iterations(), 3);
    }
    #[test]
    fn test_bspline_basis_partition_of_unity() {
        for ti in 0..10 {
            let t = ti as f64 / 10.0;
            let b = bspline_basis(t);
            let sum: f64 = b.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "t={t}: sum={sum}");
        }
    }
    #[test]
    fn test_bspline_evaluate() {
        let cps = (0..8).map(|i| vec![i as f64, 0.0]).collect::<Vec<_>>();
        let spline = BSplineTrajectory::new(cps);
        assert_eq!(spline.num_segments(), 5);
        let mid = spline.evaluate(0.5);
        assert!((mid[0] - 3.5).abs() < 1.0, "mid[0] = {}", mid[0]);
    }
    #[test]
    fn test_bspline_derivative() {
        let cps = (0..8).map(|i| vec![i as f64, 0.0]).collect::<Vec<_>>();
        let spline = BSplineTrajectory::new(cps);
        let deriv = spline.evaluate_derivative(0.5);
        assert!(deriv[0].abs() > 0.1, "deriv = {:?}", deriv);
    }
    #[test]
    fn test_bspline_arc_length() {
        let cps = (0..8).map(|i| vec![i as f64, 0.0]).collect::<Vec<_>>();
        let spline = BSplineTrajectory::new(cps);
        let length = spline.arc_length(100);
        assert!(length > 1.0, "length = {length}");
    }
    #[test]
    fn test_bspline_fit_to_waypoints() {
        let waypoints: Vec<Vec<f64>> = (0..10)
            .map(|i| vec![i as f64, (i as f64 * 0.3).sin()])
            .collect();
        let spline = BSplineTrajectory::fit_to_waypoints(&waypoints, 8);
        assert_eq!(spline.num_control_points(), 8);
        let start = spline.evaluate(0.0);
        assert!((start[0] - 0.0).abs() < 2.0);
    }
    #[test]
    fn test_smooth_trajectory() {
        let mut traj = Trajectory::new(2, 1);
        for i in 0..20 {
            let t = i as f64 * 0.1;
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: vec![t, (t * 2.0).sin()],
                control: vec![0.0],
            });
        }
        let smoothed = smooth_trajectory(&traj, 10, 30);
        assert_eq!(smoothed.num_knots(), 30);
    }
    #[test]
    fn test_multi_phase_trajectory() {
        let model = make_double_integrator();
        let mut mp = MultiPhaseTrajectory::new();
        let traj1 = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[5.0, 10.0], 0.0, 1.0, 6);
        mp.add_phase(TrajectoryPhase {
            name: "boost".to_string(),
            dynamics: model.clone(),
            trajectory: traj1,
            limits: None,
        });
        let traj2 =
            Trajectory::linear_initial_guess(2, 1, &[5.0, 10.0], &[15.0, 10.0], 1.0, 2.0, 6);
        mp.add_phase(TrajectoryPhase {
            name: "coast".to_string(),
            dynamics: model.clone(),
            trajectory: traj2,
            limits: None,
        });
        assert_eq!(mp.phases.len(), 2);
        assert!(mp.total_knots() > 0);
        assert!((mp.total_duration() - 2.0).abs() < 1e-10);
        assert!(mp.max_linkage_defect() < 1e-10);
        let flat = mp.flatten();
        assert!(flat.num_knots() > 0);
    }
    #[test]
    fn test_scale_time() {
        let traj = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[1.0, 0.0], 0.0, 1.0, 5);
        let scaled = scale_time(&traj, 2.0);
        assert!((scaled.duration() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_obstacle_penalty() {
        let mut traj = Trajectory::new(3, 1);
        for i in 0..5 {
            let x = -2.0 + i as f64;
            traj.knots.push(TrajectoryKnot {
                time: i as f64,
                state: vec![x, 0.0, 0.0],
                control: vec![0.0],
            });
        }
        let obs = vec![KeepOutZone::new([0.0, 0.0, 0.0], 0.5)];
        let penalty = obstacle_penalty(&traj, &obs);
        assert!(penalty > 0.0, "penalty = {penalty}");
        assert!(!is_obstacle_free(&traj, &obs));
    }
    #[test]
    fn test_obstacle_free() {
        let mut traj = Trajectory::new(3, 1);
        for i in 0..5 {
            traj.knots.push(TrajectoryKnot {
                time: i as f64,
                state: vec![10.0 + i as f64, 0.0, 0.0],
                control: vec![0.0],
            });
        }
        let obs = vec![KeepOutZone::new([0.0, 0.0, 0.0], 1.0)];
        assert!(is_obstacle_free(&traj, &obs));
    }
    #[test]
    fn test_jerk_cost_constant_velocity() {
        let mut traj = Trajectory::new(2, 1);
        for i in 0..10 {
            let t = i as f64 * 0.1;
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: vec![t, 1.0],
                control: vec![0.0],
            });
        }
        let jc = jerk_cost(&traj);
        assert!(jc < 1e-8, "jerk_cost = {jc}");
    }
    #[test]
    fn test_jerk_cost_nonzero() {
        let mut traj = Trajectory::new(1, 1);
        for i in 0..10 {
            let t = i as f64 * 0.1;
            traj.knots.push(TrajectoryKnot {
                time: t,
                state: vec![t * t * t],
                control: vec![0.0],
            });
        }
        let jc = jerk_cost(&traj);
        assert!(jc > 0.0, "jerk_cost = {jc}");
    }
    #[test]
    fn test_cost_gradient_direction() {
        let model = make_double_integrator();
        let cost = CostFunction::minimum_energy(2, 1);
        let controls_flat = vec![1.0; 5];
        let grad = cost_gradient_controls(
            &model,
            &[0.0, 0.0],
            &controls_flat,
            0.0,
            0.1,
            6,
            &cost,
            1e-5,
        );
        for (i, &g) in grad.iter().enumerate() {
            assert!(g > 0.0, "grad[{i}] = {g} should be positive");
        }
    }
    #[test]
    fn test_sqp_step() {
        let model = make_double_integrator();
        let cost = CostFunction::quadratic(2, 1, &[1.0, 1.0], &[0.1]);
        let bc = BoundaryConditions::initial_only(vec![0.0, 0.0]);
        let traj = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[1.0, 0.0], 0.0, 1.0, 11);
        let result = sqp_step(&model, &traj, &cost, &bc, 0.001, 10.0, 1e-6);
        assert!(result.trajectory.num_knots() > 0);
        assert!(result.cost >= 0.0);
    }
    #[test]
    fn test_free_time_gradient() {
        let model = make_double_integrator();
        let cost = CostFunction::minimum_time(2, 1);
        let traj = Trajectory::linear_initial_guess(2, 1, &[0.0, 0.0], &[1.0, 0.0], 0.0, 1.0, 11);
        let grad = free_time_gradient(&model, &traj, &cost, 1e-4);
        assert!(grad > 0.0, "free_time_grad = {grad}");
    }
    #[test]
    fn test_point_mass_2d_shooting() {
        let model = make_point_mass_2d();
        let cost = CostFunction::minimum_energy(4, 2);
        let controls: Vec<Vec<f64>> = (0..20).map(|_| vec![0.5, -0.1]).collect();
        let result = direct_shooting(&model, &[0.0, 0.0, 0.0, 0.0], &controls, 0.0, 0.1, &cost);
        assert_eq!(result.trajectory.num_knots(), 21);
        assert!((result.final_state[2] - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_point_mass_3d_shooting() {
        let model = make_point_mass_3d();
        let cost = CostFunction::minimum_energy(6, 3);
        let controls: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 0.0, -9.81]).collect();
        let result = direct_shooting(
            &model,
            &[0.0, 0.0, 100.0, 10.0, 0.0, 0.0],
            &controls,
            0.0,
            0.1,
            &cost,
        );
        assert_eq!(result.trajectory.num_knots(), 11);
    }
    #[test]
    fn test_trajectory_max_control() {
        let mut traj = Trajectory::new(2, 1);
        traj.knots.push(TrajectoryKnot {
            time: 0.0,
            state: vec![0.0, 0.0],
            control: vec![3.0],
        });
        traj.knots.push(TrajectoryKnot {
            time: 0.1,
            state: vec![0.1, 1.0],
            control: vec![-5.0],
        });
        assert!((traj.max_control_magnitude() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat_operations() {
        let a = mat_eye(2);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let c = mat_mul(&a, &b, 2);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[3] - 4.0).abs() < 1e-10);
        let t = mat_transpose(&b, 2);
        assert!((t[1] - 3.0).abs() < 1e-10);
        assert!((t[2] - 2.0).abs() < 1e-10);
    }
}
