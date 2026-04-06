//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::LbmSimulation;

#[cfg(test)]
mod tests {

    use crate::boundary::channel_walls;
    use crate::lattice::CS2;
    use crate::simulation::types::*;
    #[test]
    fn test_lbm_simulation_new_d2q9() {
        let cfg = LbmConfig::d2q9(10, 8, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        assert_eq!(sim.step_count, 0);
        assert_eq!(sim.get_velocity_field().len(), 10 * 8);
        assert_eq!(sim.get_density_field().len(), 10 * 8);
    }
    #[test]
    fn test_lbm_simulation_step_count() {
        let cfg = LbmConfig::d2q9(10, 8, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(50);
        assert_eq!(sim.step_count, 50);
    }
    #[test]
    fn test_lbm_simulation_density_conservation_2d() {
        let cfg = LbmConfig::d2q9(10, 10, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let rho_before: f64 = sim.get_density_field().iter().sum();
        sim.step_n(20);
        let rho_after: f64 = sim.get_density_field().iter().sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-8,
            "{rho_before} vs {rho_after}"
        );
    }
    #[test]
    fn test_lbm_simulation_pressure_field() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let p = sim.get_pressure_field();
        let rho = sim.get_density_field();
        for (p_k, rho_k) in p.iter().zip(rho.iter()) {
            assert!((p_k - rho_k * CS2).abs() < 1e-14);
        }
    }
    #[test]
    fn test_lbm_simulation_d3q19_runs() {
        let cfg = LbmConfig::d3q19(5, 5, 5, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(10);
        assert_eq!(sim.step_count, 10);
        assert_eq!(sim.get_velocity_field().len(), 125);
    }
    #[test]
    fn test_lbm_simulation_d3q27_runs() {
        let cfg = LbmConfig::d3q27(4, 4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(10);
        assert_eq!(sim.step_count, 10);
        assert_eq!(sim.get_density_field().len(), 64);
    }
    #[test]
    fn test_lbm_simulation_body_force() {
        let mut cfg = LbmConfig::d2q9(10, 10, 1.0 / 6.0);
        cfg.body_force = Some([1e-4, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(50);
        let vel = sim.get_velocity_field();
        let mean_ux: f64 = vel.iter().map(|v| v[0]).sum::<f64>() / vel.len() as f64;
        assert!(mean_ux > 0.0, "mean_ux = {mean_ux}");
    }
    #[test]
    fn test_lbm_simulation_smagorinsky() {
        let mut cfg = LbmConfig::d2q9(8, 8, 1.0 / 6.0);
        cfg.smagorinsky_cs = Some(0.1);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(20);
        assert_eq!(sim.step_count, 20);
    }
    #[test]
    fn test_poiseuille_flow_error() {
        let nu = 1.0 / 6.0;
        let mut cfg = LbmConfig::d2q9(5, 12, nu);
        cfg.body_force = Some([1e-4, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.set_boundaries(channel_walls(5, 12));
        sim.step_n(5000);
        let err = sim.poiseuille_flow_error();
        assert!(err.is_finite(), "err = {err}");
        assert!(err < 0.5, "err = {err}");
    }
    #[test]
    fn test_lbm_config_omega() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        assert!((cfg.omega() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_step_controller_initial_omega() {
        let sc = StepController::new(0.01, 0.5);
        let omega = sc.omega();
        assert!(omega > 0.0 && omega < 2.0, "omega = {omega}");
    }
    #[test]
    fn test_step_controller_increases_viscosity_high_ma() {
        let mut sc = StepController::new(0.01, 0.5);
        sc.nu_current = 0.01;
        sc.adjust_viscosity(0.5);
        assert!(sc.nu_current > 0.01, "Viscosity should increase");
    }
    #[test]
    fn test_step_controller_decreases_viscosity_low_ma() {
        let mut sc = StepController::new(0.01, 0.5);
        sc.nu_current = 0.1;
        sc.adjust_viscosity(0.01);
        assert!(sc.nu_current < 0.1, "Viscosity should decrease");
    }
    #[test]
    fn test_step_controller_clamps_to_bounds() {
        let mut sc = StepController::new(0.01, 0.5);
        sc.nu_current = 0.49;
        sc.adjust_viscosity(1.0);
        assert!(sc.nu_current <= 0.5, "nu = {}", sc.nu_current);
        sc.nu_current = 0.011;
        sc.adjust_viscosity(0.001);
        assert!(sc.nu_current >= 0.01, "nu = {}", sc.nu_current);
    }
    #[test]
    fn test_convergence_monitor_initial_not_converged() {
        let mut mon = ConvergenceMonitor::new(10, 1e-6, 100);
        let ux = vec![1.0; 10];
        let uy = vec![0.0; 10];
        let converged = mon.check(&ux, &uy);
        assert!(!converged, "Should not converge on first check");
    }
    #[test]
    fn test_convergence_monitor_converges_same_field() {
        let mut mon = ConvergenceMonitor::new(10, 1e-6, 100);
        let ux = vec![1.0; 10];
        let uy = vec![0.5; 10];
        mon.check(&ux, &uy);
        let converged = mon.check(&ux, &uy);
        assert!(converged, "Same field twice should converge");
    }
    #[test]
    fn test_convergence_monitor_latest_residual() {
        let mut mon = ConvergenceMonitor::new(5, 1e-6, 1);
        mon.check(&[1.0; 5], &[0.0; 5]);
        let r = mon.latest_residual();
        assert!(r > 0.0, "residual = {r}");
    }
    #[test]
    fn test_convergence_monitor_is_decreasing() {
        let mut mon = ConvergenceMonitor::new(5, 1e-6, 1);
        let uy = vec![0.0; 5];
        mon.check(&[0.0; 5], &uy);
        mon.check(&[10.0; 5], &uy);
        mon.check(&[10.5; 5], &uy);
        mon.check(&[10.55; 5], &uy);
        assert!(mon.is_decreasing());
    }
    #[test]
    fn test_snapshot_preserves_step_count() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(42);
        let snap = SimulationSnapshot::from_simulation(&sim);
        assert_eq!(snap.step_count, 42);
    }
    #[test]
    fn test_snapshot_n_cells() {
        let cfg = LbmConfig::d2q9(8, 6, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        assert_eq!(snap.n_cells(), 48);
    }
    #[test]
    fn test_snapshot_density_values() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        for &r in &snap.density {
            assert!((r - 1.0).abs() < 1e-10, "density = {r}");
        }
    }
    #[test]
    fn test_stats_record() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        sim.step_n(10);
        stats.record(&sim);
        assert_eq!(stats.mean_velocity_history.len(), 1);
        assert_eq!(stats.mean_density_history.len(), 1);
        assert_eq!(stats.mach_history.len(), 1);
    }
    #[test]
    fn test_stats_initial_zero_velocity() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        stats.record(&sim);
        assert!(
            stats.latest_mean_velocity() < 1e-10,
            "mean vel = {}",
            stats.latest_mean_velocity()
        );
        assert!(stats.latest_mach() < 1e-10, "Ma = {}", stats.latest_mach());
    }
    #[test]
    fn test_stats_is_stable_at_rest() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        stats.record(&sim);
        assert!(stats.is_stable(), "Rest state should be stable");
    }
    #[test]
    fn test_stats_forced_flow() {
        let mut cfg = LbmConfig::d2q9(8, 8, 1.0 / 6.0);
        cfg.body_force = Some([1e-5, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.step_n(100);
        let mut stats = SimulationStatistics::new();
        stats.record(&sim);
        assert!(
            stats.latest_mean_velocity() > 0.0,
            "Should have nonzero velocity"
        );
        assert!(stats.latest_mach() > 0.0, "Should have nonzero Ma");
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::boundary::channel_walls;
    use crate::simulation::types::*;
    #[test]
    fn test_mrt_params_d2q9_standard() {
        let p = MrtParamsD2Q9::standard(1.0);
        assert_eq!(p.s[0], 0.0);
        assert_eq!(p.s[3], 0.0);
        assert!((p.shear_omega() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mrt_params_d2q9_max_stable() {
        let p = MrtParamsD2Q9::max_stable();
        assert_eq!(p.s[0], 0.0);
        assert!((p.s[1] - 1.8).abs() < 1e-12);
    }
    #[test]
    fn test_mrt_params_d3q19_conserved_zero() {
        let p = MrtParamsD3Q19::uniform(1.0, 1.5);
        assert_eq!(p.s[0], 0.0);
        assert_eq!(p.s[3], 0.0);
        assert_eq!(p.s[5], 0.0);
    }
    #[test]
    fn test_guo_force_prefactor() {
        let gp = GuoForceParams::new([1e-4, 0.0, 0.0], 1.0);
        assert!((gp.prefactor() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_guo_force_term_zero_velocity() {
        let gp = GuoForceParams::new([1e-4, 0.0, 0.0], 1.0);
        let f1 = gp.d2q9_force_term(1, 0.0, 0.0);
        assert!(f1.is_finite(), "force term should be finite: {f1}");
    }
    #[test]
    fn test_guo_force_term_rest_sum() {
        let gp = GuoForceParams::new([1e-4, 0.0, 0.0], 1.0);
        let sum: f64 = (0..9).map(|i| gp.d2q9_force_term(i, 0.0, 0.0)).sum();
        assert!(sum.abs() < 1e-15, "force sum at rest = {sum}");
    }
    #[test]
    fn test_run_loop_runs_to_max() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut rl = RunLoop::new(200, 1e-12, 10);
        let (steps, converged) = rl.run(&mut sim);
        assert!(steps <= 200, "steps = {steps}");
        let _ = converged;
    }
    #[test]
    fn test_run_loop_records_residuals() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut rl = RunLoop::new(50, 1e-3, 5);
        rl.run(&mut sim);
        assert!(!rl.residuals.is_empty(), "Should have recorded residuals");
    }
    #[test]
    fn test_run_loop_body_force_converges() {
        let mut cfg = LbmConfig::d2q9(5, 10, 1.0 / 6.0);
        cfg.body_force = Some([1e-5, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.set_boundaries(channel_walls(5, 10));
        let mut rl = RunLoop::new(10000, 1e-5, 100);
        let (steps, _) = rl.run(&mut sim);
        assert!(steps <= 10000);
    }
    #[test]
    fn test_run_loop_last_residual() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut rl = RunLoop::new(100, 1e-6, 10);
        rl.run(&mut sim);
        let r = rl.last_residual();
        assert!(r.is_finite(), "last residual = {r}");
    }
    #[test]
    fn test_checkpoint_store_saves() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut store = CheckpointStore::new(3);
        store.save(&sim);
        sim.step_n(10);
        store.save(&sim);
        assert_eq!(store.len(), 2);
    }
    #[test]
    fn test_checkpoint_store_evicts_oldest() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut store = CheckpointStore::new(2);
        store.save(&sim);
        sim.step_n(5);
        store.save(&sim);
        sim.step_n(5);
        store.save(&sim);
        assert_eq!(store.len(), 2);
        assert_eq!(store.latest().unwrap().step_count, 10);
    }
    #[test]
    fn test_checkpoint_store_clear() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let mut store = CheckpointStore::new(5);
        store.save(&sim);
        store.save(&sim);
        store.clear();
        assert!(store.is_empty());
    }
    #[test]
    fn test_checkpoint_store_empty_latest() {
        let store = CheckpointStore::new(5);
        assert!(store.latest().is_none());
    }
    #[test]
    fn test_convergence_monitor_min_max() {
        let mut mon = ConvergenceMonitor::new(5, 1e-6, 1);
        let uy = vec![0.0; 5];
        mon.check(&[1.0; 5], &uy);
        mon.check(&[1.1; 5], &uy);
        mon.check(&[1.1; 5], &uy);
        let min = mon.min_residual();
        let max = mon.max_residual();
        assert!(min <= max, "min={min}, max={max}");
        assert!(min >= 0.0);
    }
    #[test]
    fn test_convergence_monitor_n_records() {
        let mut mon = ConvergenceMonitor::new(4, 1e-6, 1);
        let ux = vec![1.0; 4];
        let uy = vec![0.0; 4];
        for _ in 0..5 {
            mon.check(&ux, &uy);
        }
        assert_eq!(mon.n_records(), 5);
    }
    #[test]
    fn test_convergence_monitor_stagnating() {
        let mut mon = ConvergenceMonitor::new(4, 1e-6, 1);
        let ux = vec![1.0; 4];
        let uy = vec![0.0; 4];
        mon.check(&ux, &uy);
        for _ in 0..6 {
            mon.check(&ux, &uy);
        }
        assert!(mon.is_stagnating(), "Identical fields should be stagnating");
    }
    #[test]
    fn test_stats_min_max_mach() {
        let mut cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        cfg.body_force = Some([1e-5, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        stats.record(&sim);
        sim.step_n(100);
        stats.record(&sim);
        assert!(stats.min_mach() <= stats.max_mach());
    }
    #[test]
    fn test_stats_mean_velocity_overall() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        stats.record(&sim);
        stats.record(&sim);
        let mv = stats.mean_velocity_overall();
        assert!(mv >= 0.0 && mv.is_finite());
    }
    #[test]
    fn test_stats_n_records() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let mut stats = SimulationStatistics::new();
        assert_eq!(stats.n_records(), 0);
        stats.record(&sim);
        stats.record(&sim);
        assert_eq!(stats.n_records(), 2);
    }
    #[test]
    fn test_lbm_config_viscosity_from_omega() {
        let nu = LbmConfig::viscosity_from_omega(1.0);
        assert!((nu - 1.0 / 6.0).abs() < 1e-12, "nu = {nu}");
    }
    #[test]
    fn test_lbm_config_reynolds_number() {
        let cfg = LbmConfig::d2q9(64, 32, 1.0 / 6.0);
        let re = cfg.reynolds_number(0.1, 32.0);
        assert!((re - 19.2).abs() < 1e-10, "Re = {re}");
    }
    #[test]
    fn test_lbm_config_mach_number() {
        let cfg = LbmConfig::d2q9(10, 10, 1.0 / 6.0);
        let ma = cfg.mach_number();
        assert!(ma > 0.0 && ma < 0.3, "Ma = {ma}");
        assert!(cfg.is_incompressible());
    }
    #[test]
    fn test_lbm_config_d3q19_omega() {
        let cfg = LbmConfig::d3q19(8, 8, 8, 1.0 / 6.0);
        let omega = cfg.omega();
        assert!((omega - 1.0).abs() < 1e-12, "omega = {omega}");
    }
}
/// A callback type invoked every `output_interval` steps during a run.
pub type StepCallback = Box<dyn Fn(usize, &LbmSimulation)>;
#[cfg(test)]
mod tests_simulation_extra {

    use crate::boundary::channel_walls;
    use crate::simulation::types::*;
    #[test]
    fn test_adaptive_dt_new_default_dt() {
        let ctrl = AdaptiveDtController::new(0.1, 50);
        assert!((ctrl.dt_current - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_adaptive_dt_zero_velocity_unchanged() {
        let mut ctrl = AdaptiveDtController::new(0.1, 50);
        ctrl.update(0.0);
        assert!(
            (ctrl.dt_current - 1.0).abs() < 1e-14,
            "dt should not change at zero velocity"
        );
    }
    #[test]
    fn test_adaptive_dt_high_ma_reduces_dt() {
        let mut ctrl = AdaptiveDtController::new(0.1, 50);
        ctrl.dt_current = 1.0;
        ctrl.update(0.5);
        assert!(
            ctrl.dt_current < 1.0,
            "dt should decrease when Ma is too high: {}",
            ctrl.dt_current
        );
    }
    #[test]
    fn test_adaptive_dt_n_adjustments() {
        let mut ctrl = AdaptiveDtController::new(0.1, 50);
        ctrl.update(0.1);
        ctrl.update(0.1);
        assert_eq!(ctrl.n_adjustments(), 2);
    }
    #[test]
    fn test_adaptive_dt_is_stable_low_velocity() {
        let ctrl = AdaptiveDtController::new(0.1, 50);
        assert!(ctrl.is_stable(0.01), "Low velocity should be stable");
    }
    #[test]
    fn test_adaptive_dt_is_stable_high_velocity() {
        let ctrl = AdaptiveDtController::new(0.1, 50);
        assert!(!ctrl.is_stable(1.0), "High velocity should be unstable");
    }
    #[test]
    fn test_adaptive_dt_latest_dt_after_update() {
        let mut ctrl = AdaptiveDtController::new(0.1, 50);
        ctrl.update(0.2);
        let latest = ctrl.latest_dt();
        assert!(latest.is_finite() && latest > 0.0, "latest dt = {latest}");
    }
    #[test]
    fn test_ensemble_runner_basic() {
        let cfg = LbmConfig::d2q9(6, 6, 1.0 / 6.0);
        let mut runner = EnsembleRunner::new(cfg, 3, 20);
        runner.run_all();
        assert_eq!(runner.results.len(), 3);
    }
    #[test]
    fn test_ensemble_runner_member_indices() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut runner = EnsembleRunner::new(cfg, 4, 10);
        runner.run_all();
        for (i, r) in runner.results.iter().enumerate() {
            assert_eq!(r.member_idx, i, "member index mismatch at {i}");
        }
    }
    #[test]
    fn test_ensemble_runner_mean_velocity_finite() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut runner = EnsembleRunner::new(cfg, 2, 10);
        runner.run_all();
        let mv = runner.ensemble_mean_velocity();
        assert!(mv.is_finite(), "mean velocity = {mv}");
    }
    #[test]
    fn test_ensemble_runner_mach_finite() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut runner = EnsembleRunner::new(cfg, 2, 10);
        runner.run_all();
        let ma = runner.ensemble_mean_mach();
        assert!(ma.is_finite(), "mean Ma = {ma}");
    }
    #[test]
    fn test_ensemble_runner_velocity_variance_nonneg() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut runner = EnsembleRunner::new(cfg, 3, 20);
        runner.run_all();
        let var = runner.velocity_variance();
        assert!(var >= 0.0, "variance = {var}");
    }
    #[test]
    fn test_full_step_loop_runs() {
        let cfg = LbmConfig::d2q9(5, 5, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut fsl = FullStepLoop::new(100, 1e-6, 10, 0);
        let (steps, _) = fsl.run(&mut sim);
        assert!(steps <= 100, "steps = {steps}");
    }
    #[test]
    fn test_full_step_loop_residuals_recorded() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut fsl = FullStepLoop::new(50, 1e-6, 5, 0);
        fsl.run(&mut sim);
        assert!(!fsl.residuals.is_empty(), "No residuals recorded");
    }
    #[test]
    fn test_full_step_loop_checkpoint_saved() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut fsl = FullStepLoop::new(30, 1e-8, 5, 10);
        fsl.run(&mut sim);
        assert!(!fsl.checkpoint_store.is_empty(), "No checkpoints saved");
    }
    #[test]
    fn test_full_step_loop_adaptive_dt_enabled() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut fsl = FullStepLoop::new(50, 1e-6, 10, 0);
        fsl.enable_adaptive_dt(0.1);
        let (steps, _) = fsl.run(&mut sim);
        assert!(steps <= 50);
        assert!(fsl.adaptive_dt.is_some());
    }
    #[test]
    fn test_full_step_loop_last_residual() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let mut sim = LbmSimulation::new(cfg);
        let mut fsl = FullStepLoop::new(50, 1e-6, 5, 0);
        fsl.run(&mut sim);
        let r = fsl.last_residual();
        assert!(r.is_finite(), "last residual = {r}");
    }
    #[test]
    fn test_full_step_loop_forced_flow_converges() {
        let mut cfg = LbmConfig::d2q9(5, 10, 1.0 / 6.0);
        cfg.body_force = Some([1e-5, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.set_boundaries(channel_walls(5, 10));
        let mut fsl = FullStepLoop::new(15000, 1e-5, 200, 0);
        let (steps, converged) = fsl.run(&mut sim);
        assert!(steps <= 15000);
        let _ = converged;
    }
    #[test]
    fn test_flow_analysis_l2_norm_zero_field() {
        let vel = vec![[0.0, 0.0, 0.0]; 16];
        let norm = FlowFieldAnalysis::velocity_l2_norm(&vel);
        assert!(norm.abs() < 1e-14, "norm = {norm}");
    }
    #[test]
    fn test_flow_analysis_max_velocity_uniform() {
        let vel = vec![[0.1, 0.0, 0.0]; 10];
        let mx = FlowFieldAnalysis::max_velocity(&vel);
        assert!((mx - 0.1).abs() < 1e-14, "max = {mx}");
    }
    #[test]
    fn test_flow_analysis_mean_velocity() {
        let vel = vec![[3.0, 4.0, 0.0]; 5];
        let mean = FlowFieldAnalysis::mean_velocity(&vel);
        assert!((mean - 5.0).abs() < 1e-12, "mean = {mean}");
    }
    #[test]
    fn test_flow_analysis_relative_l2_error_identical() {
        let vel = vec![[0.1, 0.0, 0.0]; 10];
        let err = FlowFieldAnalysis::relative_l2_error(&vel, &vel);
        assert!(err.abs() < 1e-14, "Error between identical fields = {err}");
    }
    #[test]
    fn test_flow_analysis_vorticity_uniform_flow() {
        let nx = 5;
        let ny = 5;
        let vel: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; nx * ny];
        let vort = FlowFieldAnalysis::vorticity_z(&vel, nx, ny);
        for &v in vort.iter().skip(nx).take(nx * (ny - 2)) {
            assert!(
                v.abs() < 1e-13,
                "Vorticity should be 0 for uniform flow: {v}"
            );
        }
    }
    #[test]
    fn test_flow_analysis_divergence_uniform_flow() {
        let nx = 5;
        let ny = 5;
        let vel: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; nx * ny];
        let div = FlowFieldAnalysis::divergence_2d(&vel, nx, ny);
        for &d in div.iter().skip(nx).take(nx * (ny - 2)) {
            assert!(
                d.abs() < 1e-13,
                "Divergence should be 0 for uniform flow: {d}"
            );
        }
    }
    #[test]
    fn test_flow_analysis_kinetic_energy_positive() {
        let vel = vec![[0.1, 0.0, 0.0]; 8];
        let rho = vec![1.0_f64; 8];
        let ke = FlowFieldAnalysis::total_kinetic_energy(&vel, &rho);
        assert!(ke > 0.0, "KE should be positive: {ke}");
    }
    #[test]
    fn test_flow_analysis_pressure_variance_uniform() {
        let p = vec![1.0 / 3.0; 16];
        let var = FlowFieldAnalysis::pressure_variance(&p);
        assert!(var.abs() < 1e-14, "Uniform pressure variance = {var}");
    }
    #[test]
    fn test_snapshot_velocity_change_self() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        let change = snap.velocity_change(&snap);
        assert!(change.abs() < 1e-14, "Change from self = {change}");
    }
    #[test]
    fn test_snapshot_is_density_valid_at_rest() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        assert!(snap.is_density_valid(), "Density should be valid at rest");
    }
    #[test]
    fn test_snapshot_is_velocity_bounded_at_rest() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        assert!(
            snap.is_velocity_bounded(0.01),
            "Velocity should be bounded at rest"
        );
    }
    #[test]
    fn test_snapshot_mean_density_at_rest() {
        let cfg = LbmConfig::d2q9(4, 4, 1.0 / 6.0);
        let sim = LbmSimulation::new(cfg);
        let snap = SimulationSnapshot::from_simulation(&sim);
        assert!(
            (snap.mean_density() - 1.0).abs() < 1e-10,
            "mean_density = {}",
            snap.mean_density()
        );
    }
}
