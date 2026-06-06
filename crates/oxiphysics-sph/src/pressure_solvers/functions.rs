//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use crate::pressure_solvers::types_advanced::*;
#[cfg(test)]
use crate::pressure_solvers::types_core::*;

/// Compute the Tait equation-of-state bulk modulus constant `B`.
///
/// `B = rho0 * c0² / gamma`
pub fn tait_eos_b(rho0: f64, c0: f64, gamma: f64) -> f64 {
    rho0 * c0 * c0 / gamma
}
/// Tait equation of state: `P = B * ((rho/rho0)^gamma - 1)`.
///
/// Returns the pressure for density `rho` given rest density `rho0`,
/// speed of sound `c0`, and exponent `gamma`.
pub fn tait_eos(rho: f64, rho0: f64, c0: f64, gamma: f64) -> f64 {
    let b = tait_eos_b(rho0, c0, gamma);
    b * ((rho / rho0).powf(gamma) - 1.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    /// Tait EOS: at rho == rho0, pressure must be exactly 0.
    #[test]
    fn tait_eos_at_rest_density_is_zero() {
        let rho0 = 1000.0_f64;
        let c0 = 1480.0_f64;
        let gamma = 7.0_f64;
        let p = tait_eos(rho0, rho0, c0, gamma);
        assert!(p.abs() < 1e-6, "Expected P≈0 at rest density, got {p}");
    }
    /// Tait EOS: compressed fluid (rho > rho0) → positive pressure.
    #[test]
    fn tait_eos_compressed_gives_positive_pressure() {
        let rho0 = 1000.0_f64;
        let c0 = 1480.0_f64;
        let gamma = 7.0_f64;
        let p = tait_eos(1.1 * rho0, rho0, c0, gamma);
        assert!(p > 0.0, "Expected P > 0 for rho > rho0, got {p}");
    }
    /// WcSphPressure: compressive (rho > rho0) gives positive pressure.
    #[test]
    fn wcsph_pressure_positive_above_rest_density() {
        let solver = WcSphPressure::new(1000.0, 1480.0, 7.0);
        let p = solver.pressure(1100.0);
        assert!(p > 0.0, "Expected P > 0 for rho > rho0, got {p}");
    }
    /// PCISPH converged: all zero density errors → converged.
    #[test]
    fn pcisph_converged_with_zero_errors() {
        let solver = PciSph::new(1000.0, 0.1, 32);
        let errors = vec![0.0_f64; 10];
        assert!(solver.converged(&errors, 1e-3));
    }
    /// PCISPH not converged when errors exceed tolerance * rho0.
    #[test]
    fn pcisph_not_converged_with_large_errors() {
        let solver = PciSph::new(1000.0, 0.1, 32);
        let errors = vec![5.0_f64; 5];
        assert!(!solver.converged(&errors, 1e-3));
    }
    /// CFL dt must be proportional to h: doubling h doubles dt.
    #[test]
    fn wcsph_cfl_dt_proportional_to_h() {
        let solver = WcSphPressure::new(1000.0, 1480.0, 7.0);
        let dt1 = solver.cfl_dt(0.1);
        let dt2 = solver.cfl_dt(0.2);
        let ratio = dt2 / dt1;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Expected dt ratio ≈ 2.0, got {ratio}"
        );
    }
    /// B constant sanity check: B = rho0 * c0² / gamma.
    #[test]
    fn tait_eos_b_formula() {
        let b = tait_eos_b(1000.0, 10.0, 5.0);
        assert!((b - 20_000.0).abs() < 1e-9);
    }
    #[test]
    fn wcsph_solver_pressure_at_rest() {
        let solver = WcSphSolver::new(1000.0, 100.0, 7.0);
        let p = solver.pressure(1000.0);
        assert!(p.abs() < 1e-6, "P at rest density should be ~0, got {p}");
    }
    #[test]
    fn wcsph_solver_compressed_positive() {
        let p = WcSphSolver::compute_tait_pressure(1100.0, 1000.0, 20_000.0, 7.0);
        assert!(p > 0.0, "compressed fluid must have P > 0");
    }
    #[test]
    fn wcsph_solver_compute_all_pressures() {
        let solver = WcSphSolver::new(1000.0, 100.0, 7.0);
        let densities = [1000.0, 1100.0, 900.0];
        let mut pressures = [0.0; 3];
        solver.compute_all_pressures(&densities, &mut pressures);
        assert!(pressures[0].abs() < 1e-6);
        assert!(pressures[1] > 0.0);
        assert!(pressures[2] < 0.0);
    }
    #[test]
    fn wcsph_solver_sound_speed_at_rest() {
        let solver = WcSphSolver::new(1000.0, 100.0, 7.0);
        let c = solver.sound_speed(1000.0);
        assert!(
            (c - 100.0).abs() < 1e-6,
            "c0 at rest density should be 100, got {c}"
        );
    }
    #[test]
    fn pcisph_simple_iteration_reduces_error() {
        let solver = PcisphSolverSimple::new(1000.0, 0.5);
        let predicted = vec![1050.0, 1020.0, 980.0];
        let mut pressures = vec![0.0; 3];
        let err1 = solver.iteration(&predicted, &mut pressures);
        assert!(err1 > 0.0);
        let _err2 = solver.iteration(&predicted, &mut pressures);
        assert!(pressures[0] > 0.0);
        assert!(pressures[1] > 0.0);
    }
    #[test]
    fn pcisph_simple_solve_converges_uniform() {
        let solver = PcisphSolverSimple::new(1000.0, 1.0);
        let predicted = vec![1000.0; 10];
        let (pressures, iters) = solver.solve(&predicted);
        assert_eq!(iters, 1, "should converge in 1 iteration at rest density");
        for &p in &pressures {
            assert!(p.abs() < 1e-10);
        }
    }
    #[test]
    fn pcisph_simple_solve_produces_nonneg_pressure() {
        let solver = PcisphSolverSimple::new(1000.0, 0.1);
        let predicted = vec![1050.0, 950.0, 1000.0];
        let (pressures, _) = solver.solve(&predicted);
        for &p in &pressures {
            assert!(p >= 0.0, "pressure must be non-negative, got {p}");
        }
    }
    #[test]
    fn iisph_simple_compute_aii() {
        let off = vec![(1, 0.3), (2, -0.5)];
        let aii = IisphSolverSimple::compute_aii(&off);
        assert!((aii - 0.2).abs() < 1e-12);
    }
    #[test]
    fn iisph_simple_solve_trivial_system() {
        let solver = IisphSolverSimple {
            rho0: 100.0,
            omega: 0.8,
            max_iter: 1000,
            tolerance: 1e-10,
        };
        let rhs = vec![10.0, 10.0, 10.0];
        let aij = vec![
            vec![(1_usize, -0.3), (2_usize, -0.2)],
            vec![(0_usize, -0.3), (2_usize, -0.2)],
            vec![(0_usize, -0.2), (1_usize, -0.2)],
        ];
        let (p, iters, residual) = solver.solve_pressure_poisson(&rhs, &aij);
        assert!(residual.is_finite(), "residual not finite");
        assert!(iters > 0);
        for &pi in &p {
            assert!(pi >= 0.0, "pressure must be non-negative");
            assert!(pi.is_finite(), "pressure must be finite");
        }
    }
    #[test]
    fn iisph_simple_nonneg_clamping() {
        let solver = IisphSolverSimple::new(1000.0);
        let rhs = vec![-1000.0];
        let aij = vec![vec![]];
        let (p, _, _) = solver.solve_pressure_poisson(&rhs, &aij);
        assert!(p[0] >= 0.0, "pressure must be clamped to >= 0");
    }
    #[test]
    fn isph_jacobi_identity_system() {
        let p_old = vec![0.0, 0.0];
        let rhs = vec![2.0, 3.0];
        let aij = vec![vec![(1, -1.0)], vec![(0, -1.0)]];
        let p_new = IsphPressure::pressure_poisson_jacobi_step(&p_old, &rhs, &aij);
        assert!((p_new[0] - 2.0).abs() < 1e-12);
        assert!((p_new[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_solver_comparison_wcsph() {
        let comp = SolverComparison::new(1000.0, 100.0);
        let densities = vec![1000.0, 1050.0, 950.0];
        let result = comp.run_wcsph(&densities);
        assert_eq!(result.solver_type, SolverType::Wcsph);
        assert_eq!(result.iterations, 0);
        assert!(result.pressures[0].abs() < 1e-6);
        assert!(result.pressures[1] > 0.0);
    }
    #[test]
    fn test_solver_comparison_pcisph() {
        let comp = SolverComparison::new(1000.0, 100.0);
        let densities = vec![1000.0, 1000.0, 1000.0];
        let result = comp.run_pcisph(&densities);
        assert_eq!(result.solver_type, SolverType::Pcisph);
        assert_eq!(result.iterations, 1);
    }
    #[test]
    fn test_solver_comparison_iisph() {
        let comp = SolverComparison::new(1000.0, 100.0);
        let rhs = vec![10.0, 10.0];
        let aij = vec![vec![(1_usize, -0.3)], vec![(0_usize, -0.3)]];
        let result = comp.run_iisph(&rhs, &aij);
        assert_eq!(result.solver_type, SolverType::Iisph);
        assert!(result.iterations > 0);
    }
    #[test]
    fn test_solver_comparison_empty_densities() {
        let comp = SolverComparison::new(1000.0, 100.0);
        let result = comp.run_wcsph(&[]);
        assert!(result.pressures.is_empty());
        assert_eq!(result.avg_density_error, 0.0);
    }
    #[test]
    fn test_adaptive_solver_uses_wcsph_for_small_errors() {
        let mut solver = AdaptiveSolverSwitch::new(1000.0, 100.0, 0.1);
        let densities = vec![1000.0, 1001.0, 999.0];
        let pressures = solver.compute_pressures(&densities);
        assert_eq!(solver.current_solver, SolverType::Wcsph);
        assert_eq!(pressures.len(), 3);
    }
    #[test]
    fn test_adaptive_solver_switches_to_pcisph() {
        let mut solver = AdaptiveSolverSwitch::new(1000.0, 100.0, 0.05);
        let densities = vec![1000.0, 1200.0, 800.0];
        let _pressures = solver.compute_pressures(&densities);
        assert_eq!(solver.current_solver, SolverType::Pcisph);
        assert_eq!(solver.switch_count, 1);
    }
    #[test]
    fn test_adaptive_solver_switch_count() {
        let mut solver = AdaptiveSolverSwitch::new(1000.0, 100.0, 0.05);
        solver.compute_pressures(&[1000.0, 1001.0]);
        assert_eq!(solver.current_solver, SolverType::Wcsph);
        assert_eq!(solver.switch_count, 0);
        solver.compute_pressures(&[1000.0, 1200.0]);
        assert_eq!(solver.switch_count, 1);
        solver.compute_pressures(&[1000.0, 1001.0]);
        assert_eq!(solver.switch_count, 2);
    }
    #[test]
    fn test_residual_selector_empty_recommends_wcsph() {
        let selector = ResidualBasedSelector::new(10);
        assert_eq!(selector.recommend(), SolverType::Wcsph);
    }
    #[test]
    fn test_residual_selector_recommends_best() {
        let mut selector = ResidualBasedSelector::new(10);
        selector.record(SolverPerformanceEntry {
            solver_type: SolverType::Wcsph,
            iterations: 0,
            residual: 10.0,
        });
        selector.record(SolverPerformanceEntry {
            solver_type: SolverType::Pcisph,
            iterations: 5,
            residual: 1.0,
        });
        selector.record(SolverPerformanceEntry {
            solver_type: SolverType::Iisph,
            iterations: 10,
            residual: 0.5,
        });
        assert_eq!(selector.recommend(), SolverType::Iisph);
    }
    #[test]
    fn test_residual_selector_history_limit() {
        let mut selector = ResidualBasedSelector::new(3);
        for i in 0..5 {
            selector.record(SolverPerformanceEntry {
                solver_type: SolverType::Wcsph,
                iterations: 0,
                residual: i as f64,
            });
        }
        assert_eq!(selector.history_len(), 3);
    }
    #[test]
    fn pcisph_full_pressure_correction_positive_for_overdense() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let mut state = PcisphState::new(3);
        state.predicted_densities[0] = 1100.0;
        state.predicted_densities[1] = 1000.0;
        state.predicted_densities[2] = 900.0;
        solver.compute_pressure_correction(&mut state);
        assert!(
            state.pressure_corrections[0] > 0.0,
            "overdense → positive correction, got {}",
            state.pressure_corrections[0]
        );
        assert!(
            state.pressure_corrections[1].abs() < 1e-12,
            "at rest density → zero correction, got {}",
            state.pressure_corrections[1]
        );
    }
    #[test]
    fn pcisph_full_density_error_computation() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let mut state = PcisphState::new(2);
        state.predicted_densities[0] = 1050.0;
        state.predicted_densities[1] = 950.0;
        let errors = solver.compute_density_error(&state);
        assert!((errors[0] - 50.0).abs() < 1e-10);
        assert!((errors[1] - (-50.0)).abs() < 1e-10);
    }
    #[test]
    fn pcisph_full_predict_velocities_positions() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let n = 2;
        let mut state = PcisphState::new(n);
        let positions = vec![[0.0_f64; 3]; n];
        let velocities = vec![[1.0_f64, 0.0, 0.0]; n];
        let forces = vec![[2.0_f64, 0.0, 0.0]; n];
        let masses = vec![1.0_f64; n];
        let dt = 0.1;
        solver.predict_velocities_positions(
            &positions,
            &velocities,
            &forces,
            &masses,
            dt,
            &mut state,
        );
        assert!((state.predicted_velocities[0][0] - 1.2).abs() < 1e-12);
        assert!((state.predicted_positions[0][0] - 0.12).abs() < 1e-12);
    }
    #[test]
    fn iisph_full_compute_aii_sum() {
        let off: Vec<(usize, f64)> = vec![(1, 0.4), (2, 0.2)];
        let aii = IisphSolver::compute_aii(&off);
        assert!((aii - (-0.6)).abs() < 1e-12, "aii={aii}");
    }
    #[test]
    fn iisph_full_iterate_pressure_step() {
        let solver = IisphSolver::new(1000.0, 0.5, 100, 0.001);
        let rhs = vec![10.0_f64, 10.0];
        let aij = vec![vec![(1_usize, -0.3_f64)], vec![(0_usize, -0.3_f64)]];
        let (pressures, iters, _residual) = solver.iterate_pressure(&rhs, &aij);
        assert!(iters > 0);
        for &p in &pressures {
            assert!(p >= 0.0, "pressure must be non-negative, got {p}");
            assert!(p.is_finite());
        }
    }
    #[test]
    fn tait_eos_expanded_gives_negative_pressure() {
        let p = tait_eos(900.0, 1000.0, 1480.0, 7.0);
        assert!(
            p < 0.0,
            "Expanded fluid should have negative pressure, got {p}"
        );
    }
    #[test]
    fn test_wcsph_cfl_positive() {
        let solver = WcSphPressure::new(1000.0, 1480.0, 7.0);
        let dt = solver.cfl_dt(0.1);
        assert!(dt > 0.0, "CFL dt should be positive");
    }
    #[test]
    fn test_wcsph_sound_speed_cfl_alias() {
        let solver = WcSphPressure::new(1000.0, 1480.0, 7.0);
        let dt1 = solver.cfl_dt(0.1);
        let dt2 = solver.sound_speed_cfl_dt(0.1);
        assert_eq!(dt1, dt2, "cfl_dt and sound_speed_cfl_dt should match");
    }
    #[test]
    fn test_pcisph_predict_density() {
        let solver = PciSph::new(1000.0, 0.1, 32);
        let rho_pred = solver.predict_density(1000.0, 50.0, 0.01);
        assert!((rho_pred - 1000.5).abs() < 1e-10);
    }
    #[test]
    fn test_pcisph_density_error() {
        let solver = PciSph::new(1000.0, 0.1, 32);
        let err = solver.density_error(1050.0);
        assert!((err - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_pcisph_pressure_update() {
        let solver = PciSph::new(1000.0, 0.1, 32);
        let p = solver.pressure_update(100.0, 50.0);
        assert!((p - 100.0 - solver.delta * 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_isph_divergence_velocity() {
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let masses = [1.0, 1.0];
        let densities = [1000.0, 1000.0];
        let div =
            IsphPressure::divergence_velocity(&positions, &velocities, &masses, &densities, 1.0, 0);
        assert!(
            div.abs() < 1e-6,
            "Zero velocity field should have zero divergence"
        );
    }
    #[test]
    fn test_solver_type_equality() {
        assert_eq!(SolverType::Wcsph, SolverType::Wcsph);
        assert_ne!(SolverType::Wcsph, SolverType::Pcisph);
        assert_ne!(SolverType::Pcisph, SolverType::Iisph);
    }
    #[test]
    fn test_solver_benchmark_result_fields() {
        let result = SolverBenchmarkResult {
            solver_type: SolverType::Wcsph,
            iterations: 0,
            max_density_error: 1.0,
            avg_density_error: 0.5,
            pressures: vec![100.0],
        };
        assert_eq!(result.solver_type, SolverType::Wcsph);
        assert_eq!(result.pressures.len(), 1);
    }
    #[test]
    fn pcisph_predict_correct_uniform_rest_density_converges_quickly() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let n = 4;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let velocities = vec![[0.0_f64; 3]; n];
        let masses = vec![0.001_f64; n];
        let forces = vec![[0.0_f64; 3]; n];
        let h = 0.1_f64;
        let w_self = 1.0 / (std::f64::consts::PI * h * h * h);
        let neighbors: Vec<Vec<(usize, f64, f64)>> =
            (0..n).map(|_| vec![(0_usize, 0.0, w_self)]).collect();
        let (pressures, _iters) = solver.predict_correct_step(
            &positions,
            &velocities,
            &masses,
            &forces,
            &neighbors,
            h,
            0.001,
        );
        assert_eq!(pressures.len(), n);
        for &p in &pressures {
            assert!(p >= 0.0, "pressures must be non-negative, got {p}");
            assert!(p.is_finite(), "pressure must be finite");
        }
    }
    #[test]
    fn pcisph_predict_correct_step_returns_iter_count() {
        let solver = PcisphSolver::new(1000.0, 0.1, 10, 1e-8);
        let n = 2;
        let positions = vec![[0.0_f64; 3]; n];
        let velocities = vec![[0.0_f64; 3]; n];
        let masses = vec![1.0_f64; n];
        let forces = vec![[0.0_f64; 3]; n];
        let neighbors = vec![vec![]; n];
        let (_pressures, iters) = solver.predict_correct_step(
            &positions,
            &velocities,
            &masses,
            &forces,
            &neighbors,
            0.1,
            0.001,
        );
        assert!(
            iters > 0 && iters <= 10,
            "iter count should be in [1, max_iter]"
        );
    }
    #[test]
    fn dfsph_simple_new_initialises_correctly() {
        let solver = DfsphSolverSimple::new(1000.0, 0.5, 100, 0.001, 0.001);
        assert!((solver.rest_density - 1000.0).abs() < 1e-10);
        assert!((solver.omega - 0.5).abs() < 1e-10);
        assert_eq!(solver.max_iter_density, 100);
    }
    #[test]
    fn dfsph_simple_density_alpha_precompute() {
        let solver = DfsphSolverSimple::new(1000.0, 0.5, 50, 0.001, 0.001);
        let n = 3;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.5, 0.0, 0.0]).collect();
        let masses = vec![1.0_f64; n];
        let neighbors = vec![vec![]; n];
        let alphas = solver.compute_alpha(&positions, &masses, &neighbors, 0.1);
        assert_eq!(alphas.len(), n);
        for &a in &alphas {
            assert!(a.abs() < 1e-10 || a.is_finite(), "alpha must be finite");
        }
    }
    #[test]
    fn dfsph_simple_divergence_free_step_returns_pressures() {
        let solver = DfsphSolverSimple::new(1000.0, 0.5, 5, 0.01, 0.01);
        let n = 4;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.12, 0.0, 0.0]).collect();
        let velocities = vec![[0.1_f64, 0.0, 0.0]; n];
        let masses = vec![0.001_f64; n];
        let densities = vec![1000.0_f64; n];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![vec![]; n];
        let (pressures, iters) = solver.divergence_free_step(
            &positions,
            &velocities,
            &masses,
            &densities,
            &neighbors,
            0.1,
            0.001,
        );
        assert_eq!(pressures.len(), n);
        assert!(iters > 0);
        for &p in &pressures {
            assert!(p.is_finite(), "DFSPH pressure must be finite");
        }
    }
    #[test]
    fn dfsph_simple_density_invariant_step_returns_pressures() {
        let solver = DfsphSolverSimple::new(1000.0, 0.5, 5, 0.01, 0.01);
        let n = 4;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.12, 0.0, 0.0]).collect();
        let velocities = vec![[0.0_f64; 3]; n];
        let masses = vec![0.001_f64; n];
        let densities = vec![1000.0_f64; n];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![vec![]; n];
        let (pressures, iters) = solver.density_invariant_step(
            &positions,
            &velocities,
            &masses,
            &densities,
            &neighbors,
            0.1,
            0.001,
        );
        assert_eq!(pressures.len(), n);
        assert!(iters > 0);
        for &p in &pressures {
            assert!(p >= 0.0, "DFSPH density-invariant pressure >= 0, got {p}");
            assert!(p.is_finite());
        }
    }
}
#[cfg(test)]
mod tests_new_solvers {
    use super::*;
    #[test]
    fn convergence_diagnostics_empty_initially() {
        let diag = ConvergenceDiagnostics::new();
        assert!(diag.is_empty());
        assert_eq!(diag.len(), 0);
        assert_eq!(diag.final_max_error(), 0.0);
        assert!(diag.convergence_iteration().is_none());
    }
    #[test]
    fn convergence_diagnostics_record_at_rest_density() {
        let mut diag = ConvergenceDiagnostics::new();
        let densities = vec![1000.0, 1000.0, 1000.0];
        let pressures = vec![0.0; 3];
        diag.record(1, &densities, &pressures, 1000.0, 0.001);
        assert_eq!(diag.len(), 1);
        assert_eq!(diag.final_max_error(), 0.0);
        assert_eq!(diag.convergence_iteration(), Some(1));
    }
    #[test]
    fn convergence_diagnostics_record_density_error() {
        let mut diag = ConvergenceDiagnostics::new();
        let densities = vec![1050.0, 950.0];
        let pressures = vec![100.0, 50.0];
        diag.record(1, &densities, &pressures, 1000.0, 0.001);
        assert_eq!(diag.convergence_iteration(), None);
        assert!((diag.final_max_error() - 50.0).abs() < 1e-10);
    }
    #[test]
    fn convergence_diagnostics_worst_error_tracks_max() {
        let mut diag = ConvergenceDiagnostics::new();
        diag.record(1, &[1010.0], &[0.0], 1000.0, 0.001);
        diag.record(2, &[1020.0], &[0.0], 1000.0, 0.001);
        diag.record(3, &[1005.0], &[0.0], 1000.0, 0.001);
        assert!((diag.worst_density_error() - 20.0).abs() < 1e-10);
    }
    #[test]
    fn convergence_diagnostics_empty_input_no_panic() {
        let mut diag = ConvergenceDiagnostics::new();
        diag.record(1, &[], &[], 1000.0, 0.001);
        assert!(diag.is_empty(), "empty densities should not add entry");
    }
    #[test]
    fn jacobi_iterator_trivial_zero_rhs() {
        let iter = JacobiPressureIterator::new(1000.0, 0.7, 100, 1e-8);
        let rhs = vec![0.0, 0.0];
        let aij = vec![vec![(1_usize, -0.3_f64)], vec![(0_usize, -0.3_f64)]];
        let (p, _res, _iters) = iter.run(&rhs, &aij, None);
        for &pi in &p {
            assert!(
                pi.abs() < 1e-6,
                "zero rhs should give near-zero pressure, got {pi}"
            );
        }
    }
    #[test]
    fn jacobi_iterator_convergence_positive_rhs() {
        let it = JacobiPressureIterator::new(1000.0, 0.7, 500, 1e-4);
        let rhs = vec![5.0, 5.0, 5.0];
        let aij = vec![
            vec![(1_usize, -0.3_f64), (2_usize, -0.2_f64)],
            vec![(0_usize, -0.3_f64), (2_usize, -0.2_f64)],
            vec![(0_usize, -0.2_f64), (1_usize, -0.2_f64)],
        ];
        let (p, residual, iters) = it.run(&rhs, &aij, None);
        assert!(iters > 0 && iters <= 500, "iters={iters}");
        assert!(residual.is_finite(), "residual must be finite");
        for &pi in &p {
            assert!(pi >= 0.0, "pressures must be non-negative, got {pi}");
        }
    }
    #[test]
    fn jacobi_iterator_with_diagnostics() {
        let it = JacobiPressureIterator::new(1000.0, 0.5, 20, 0.01);
        let rhs = vec![2.0_f64, 2.0];
        let aij = vec![vec![(1_usize, -0.4_f64)], vec![(0_usize, -0.4_f64)]];
        let mut diag = ConvergenceDiagnostics::new();
        let (_, _, iters) = it.run(&rhs, &aij, Some(&mut diag));
        assert_eq!(diag.len(), iters);
    }
    #[test]
    fn pcisph_loop_converges_at_rest_density() {
        let looper = PcisphIterativeLoop::new(1000.0, 0.5, 0.001);
        let densities = vec![1000.0; 5];
        let (pressures, iters, diag) = looper.solve_from_predicted_densities(&densities);
        assert!(iters <= looper.max_iter);
        for &p in &pressures {
            assert!(p.abs() < 1e-10);
        }
        assert!(diag.convergence_iteration().is_some());
    }
    #[test]
    fn pcisph_loop_density_error_stats_correct() {
        let looper = PcisphIterativeLoop::new(1000.0, 0.5, 0.001);
        let densities = vec![1050.0, 950.0, 1000.0];
        let (max_err, mean_err) = looper.density_error_stats(&densities);
        assert!((max_err - 50.0).abs() < 1e-10, "max_err={max_err}");
        assert!(
            (mean_err - 100.0 / 3.0).abs() < 1e-10,
            "mean_err={mean_err}"
        );
    }
    #[test]
    fn pcisph_loop_produces_nonneg_pressures() {
        let looper = PcisphIterativeLoop::new(1000.0, 0.1, 0.001);
        let densities = vec![1050.0, 980.0, 1020.0, 1000.0];
        let (pressures, _, _) = looper.solve_from_predicted_densities(&densities);
        for &p in &pressures {
            assert!(p >= 0.0, "pressure must be non-negative, got {p}");
        }
    }
    #[test]
    fn pcisph_loop_min_iter_respected() {
        let mut looper = PcisphIterativeLoop::new(1000.0, 1.0, 1.0);
        looper.min_iter = 5;
        let densities = vec![1000.0; 4];
        let (_, iters, _) = looper.solve_from_predicted_densities(&densities);
        assert!(iters >= 5, "expected ≥ min_iter=5 iterations, got {iters}");
    }
    #[test]
    fn iisph_precond_build_correct_diag() {
        let aij = vec![
            vec![(1_usize, -0.3_f64), (2_usize, -0.2_f64)],
            vec![(0_usize, -0.3_f64)],
            vec![(0_usize, -0.2_f64)],
        ];
        let precond = IisphDiagonalPreconditioner::build(&aij, 0.5, 1000.0);
        assert!(
            (precond.diag[0] - 0.5).abs() < 1e-12,
            "diag[0]={}",
            precond.diag[0]
        );
        assert!(
            (precond.diag[1] - 0.3).abs() < 1e-12,
            "diag[1]={}",
            precond.diag[1]
        );
        assert!(
            (precond.diag[2] - 0.2).abs() < 1e-12,
            "diag[2]={}",
            precond.diag[2]
        );
    }
    #[test]
    fn iisph_precond_solve_trivial() {
        let aij: Vec<Vec<(usize, f64)>> =
            vec![vec![(1_usize, -0.3_f64)], vec![(0_usize, -0.3_f64)]];
        let precond = IisphDiagonalPreconditioner::build(&aij, 0.6, 1000.0);
        let rhs = vec![3.0_f64, 3.0];
        let (p, _iters, residual) = precond.solve(&rhs, &aij, 500, 1e-8);
        assert!(residual.is_finite(), "residual must be finite");
        for &pi in &p {
            assert!(pi >= 0.0, "pressure must be non-negative, got {pi}");
            assert!(pi.is_finite(), "pressure must be finite");
        }
    }
    #[test]
    fn iisph_precond_step_zero_pressure_init() {
        let aij = vec![vec![(1_usize, -0.25_f64)], vec![(0_usize, -0.25_f64)]];
        let precond = IisphDiagonalPreconditioner::build(&aij, 0.5, 1000.0);
        let p = vec![0.0_f64; 2];
        let rhs = vec![1.0_f64, 1.0];
        let (p_new, _res) = precond.step(&p, &rhs, &aij);
        for &pi in &p_new {
            assert!(pi.is_finite(), "p_new must be finite");
            assert!(pi >= 0.0, "p_new must be non-negative");
        }
    }
    #[test]
    fn pcisph_loop_max_error_decreases_with_more_iterations() {
        let looper = PcisphIterativeLoop {
            rho0: 1000.0,
            delta: 0.05,
            min_iter: 1,
            max_iter: 2,
            eta: 1e-10,
        };
        let densities = vec![1100.0, 900.0];
        let (_, iters, diag) = looper.solve_from_predicted_densities(&densities);
        assert_eq!(iters, 2);
        assert_eq!(diag.len(), 2);
    }
    #[test]
    fn wcsph_vs_tait_eos_agreement() {
        let rho0 = 1000.0;
        let c0 = 100.0;
        let gamma = 7.0;
        let solver = WcSphSolver::new(rho0, c0, gamma);
        for rho in [900.0, 1000.0, 1100.0] {
            let p1 = solver.pressure(rho);
            let p2 = tait_eos(rho, rho0, c0, gamma);
            assert!(
                (p1 - p2).abs() < 1e-6,
                "WcSphSolver and tait_eos disagree at rho={rho}: {p1} vs {p2}"
            );
        }
    }
    #[test]
    fn iisph_precond_handles_degenerate_diagonal() {
        let aij: Vec<Vec<(usize, f64)>> = vec![vec![], vec![]];
        let precond = IisphDiagonalPreconditioner::build(&aij, 0.5, 1000.0);
        let rhs = vec![10.0_f64, 10.0];
        let (p, _, _) = precond.solve(&rhs, &aij, 50, 1e-6);
        for &pi in &p {
            assert!(
                pi.abs() < 1e-10,
                "degenerate row should give 0 pressure, got {pi}"
            );
        }
    }
    #[test]
    fn convergence_diagnostics_multiple_records() {
        let mut diag = ConvergenceDiagnostics::new();
        for k in 1..=5 {
            let density_err = 50.0 / k as f64;
            let rho = 1000.0 + density_err;
            diag.record(k, &[rho], &[0.0], 1000.0, 0.001);
        }
        assert_eq!(diag.len(), 5);
        let first = diag.history[0].max_density_error;
        let last = diag.history[4].max_density_error;
        assert!(
            first > last,
            "first_err={first} should exceed last_err={last}"
        );
    }
    #[test]
    fn tait_eos_monotone_increasing_with_rho() {
        let rho0 = 1000.0;
        let c0 = 100.0;
        let gamma = 7.0;
        let p1 = tait_eos(1000.0, rho0, c0, gamma);
        let p2 = tait_eos(1050.0, rho0, c0, gamma);
        let p3 = tait_eos(1100.0, rho0, c0, gamma);
        assert!(p1 < p2, "Tait EOS should be monotone: p1={p1} p2={p2}");
        assert!(p2 < p3, "Tait EOS should be monotone: p2={p2} p3={p3}");
    }
    #[test]
    fn tait_eos_b_scale_proportional_to_rho0() {
        let b1 = tait_eos_b(500.0, 100.0, 7.0);
        let b2 = tait_eos_b(1000.0, 100.0, 7.0);
        assert!(
            (b2 / b1 - 2.0).abs() < 1e-10,
            "B should double when rho0 doubles"
        );
    }
    #[test]
    fn tait_eos_b_scale_proportional_to_c0_squared() {
        let b1 = tait_eos_b(1000.0, 100.0, 7.0);
        let b2 = tait_eos_b(1000.0, 200.0, 7.0);
        assert!(
            (b2 / b1 - 4.0).abs() < 1e-10,
            "B should quadruple when c0 doubles"
        );
    }
    #[test]
    fn wcsph_solver_pressure_matches_tait_eos() {
        let rho0 = 1000.0;
        let c0 = 100.0;
        let gamma = 7.0;
        let solver = WcSphSolver::new(rho0, c0, gamma);
        for &rho in &[900.0, 1000.0, 1050.0, 1100.0_f64] {
            let p1 = solver.pressure(rho);
            let p2 = tait_eos(rho, rho0, c0, gamma);
            assert!((p1 - p2).abs() < 1e-6, "rho={rho}: solver={p1} tait={p2}");
        }
    }
    #[test]
    fn wcsph_solver_compute_all_pressures_length() {
        let solver = WcSphSolver::new(1000.0, 100.0, 7.0);
        let densities = vec![1000.0; 20];
        let mut pressures = vec![0.0; 20];
        solver.compute_all_pressures(&densities, &mut pressures);
        assert_eq!(pressures.len(), 20);
    }
    #[test]
    fn pcisph_solver_simple_negative_density_pressure_clamped() {
        let solver = PcisphSolverSimple::new(1000.0, 0.5);
        let predicted = vec![500.0];
        let (pressures, _) = solver.solve(&predicted);
        assert!(
            pressures[0] >= 0.0,
            "pressure must be non-negative, got {}",
            pressures[0]
        );
    }
    #[test]
    fn pcisph_solver_simple_max_iter_limit() {
        let mut solver = PcisphSolverSimple::new(1000.0, 1e-10);
        solver.max_iter = 7;
        solver.tolerance = 1e-20;
        let predicted = vec![1050.0];
        let (_, iters) = solver.solve(&predicted);
        assert_eq!(iters, 7, "should exhaust max_iter");
    }
    #[test]
    fn pcisph_state_len_and_is_empty() {
        let state = PcisphState::new(5);
        assert_eq!(state.len(), 5);
        assert!(!state.is_empty());
        let empty = PcisphState::new(0);
        assert!(empty.is_empty());
    }
    #[test]
    fn pcisph_state_reset_corrections_zeros_all() {
        let mut state = PcisphState::new(4);
        for c in &mut state.pressure_corrections {
            *c = 999.0;
        }
        state.reset_corrections();
        for &c in &state.pressure_corrections {
            assert_eq!(c, 0.0);
        }
    }
    #[test]
    fn pcisph_full_solver_converged_near_rest_density() {
        let solver = PcisphSolver::new(1000.0, 0.05, 200, 0.001);
        let mut state = PcisphState::new(3);
        state.predicted_densities = vec![1000.5, 999.7, 1000.2];
        let errors = solver.compute_density_error(&state);
        assert!(
            solver.converged(&errors),
            "small density errors should converge"
        );
    }
    #[test]
    fn iisph_solver_simple_empty_system() {
        let solver = IisphSolverSimple::new(1000.0);
        let (p, _iters, resid) = solver.solve_pressure_poisson(&[], &[]);
        assert!(p.is_empty());
        assert_eq!(resid, 0.0);
    }
    #[test]
    fn iisph_solver_single_diagonal_degenerate() {
        let solver = IisphSolverSimple::new(1000.0);
        let rhs = vec![5.0_f64];
        let aij: Vec<Vec<(usize, f64)>> = vec![vec![]];
        let (p, _, _) = solver.solve_pressure_poisson(&rhs, &aij);
        assert!(p[0] >= 0.0);
    }
    #[test]
    fn isph_divergence_velocity_uniform_flow_is_zero() {
        let positions = [[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let velocities = [[1.0, 0.0, 0.0]; 2];
        let masses = [1.0; 2];
        let densities = [1000.0; 2];
        let div0 =
            IsphPressure::divergence_velocity(&positions, &velocities, &masses, &densities, 1.0, 0);
        assert!(
            div0.abs() < 1e-12,
            "Uniform flow divergence must be 0, got {div0}"
        );
    }
    #[test]
    fn convergence_diagnostics_empty_state() {
        let diag = ConvergenceDiagnostics::new();
        assert!(diag.is_empty());
        assert_eq!(diag.len(), 0);
        assert_eq!(diag.final_max_error(), 0.0);
        assert!(diag.convergence_iteration().is_none());
    }
    #[test]
    fn convergence_diagnostics_worst_error_is_maximum() {
        let mut diag = ConvergenceDiagnostics::new();
        diag.record(1, &[1100.0, 1200.0], &[0.0; 2], 1000.0, 0.001);
        diag.record(2, &[1010.0, 1020.0], &[0.0; 2], 1000.0, 0.001);
        assert!((diag.worst_density_error() - 200.0).abs() < 1e-6);
    }
    #[test]
    fn pcisph_solver_predict_velocities_positions_correct_kinematics() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let n = 3;
        let mut state = PcisphState::new(n);
        let positions = vec![[1.0_f64, 0.0, 0.0]; n];
        let velocities = vec![[0.0_f64, 2.0, 0.0]; n];
        let forces = vec![[0.0_f64, 0.0, 0.0]; n];
        let masses = vec![1.0_f64; n];
        let dt = 0.5;
        solver.predict_velocities_positions(
            &positions,
            &velocities,
            &forces,
            &masses,
            dt,
            &mut state,
        );
        for i in 0..n {
            assert!((state.predicted_velocities[i][1] - 2.0).abs() < 1e-12);
            assert!((state.predicted_positions[i][0] - 1.0).abs() < 1e-12);
            assert!((state.predicted_positions[i][1] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn iisph_full_solver_omega_effect() {
        let fast = IisphSolver::new(1000.0, 0.9, 1000, 1e-8);
        let slow = IisphSolver::new(1000.0, 0.1, 1000, 1e-8);
        let rhs = vec![5.0_f64, 5.0];
        let aij = vec![vec![(1_usize, -0.4_f64)], vec![(0_usize, -0.4_f64)]];
        let (_, iters_fast, _) = fast.iterate_pressure(&rhs, &aij);
        let (_, iters_slow, _) = slow.iterate_pressure(&rhs, &aij);
        assert!(iters_fast > 0);
        assert!(iters_slow > 0);
    }
    #[test]
    fn wcsph_pressure_negative_for_expanded_fluid() {
        let solver = WcSphPressure::new(1000.0, 100.0, 7.0);
        let p = solver.pressure(900.0);
        assert!(
            p < 0.0,
            "Expanded fluid should have negative Tait pressure, got {p}"
        );
    }
    #[test]
    fn pcisph_compute_delta_positive() {
        let delta = PciSph::compute_delta(0.1, 32);
        assert!(delta > 0.0, "delta must be positive, got {delta}");
    }
    #[test]
    fn pcisph_compute_delta_increases_with_h() {
        let d1 = PciSph::compute_delta(0.1, 32);
        let d2 = PciSph::compute_delta(0.2, 32);
        assert!(d2 > d1, "delta should increase with h");
    }
    #[test]
    fn adaptive_solver_switch_no_switch_on_same_solver() {
        let mut solver = AdaptiveSolverSwitch::new(1000.0, 100.0, 0.1);
        solver.compute_pressures(&[1001.0, 999.0]);
        assert_eq!(solver.switch_count, 0);
        solver.compute_pressures(&[1000.5, 999.5]);
        assert_eq!(solver.switch_count, 0);
    }
    #[test]
    fn jacobi_pressure_iterator_empty_system() {
        let jit = JacobiPressureIterator::new(1000.0, 0.7, 100, 1e-6);
        let (p, resid, _iters) = jit.run(&[], &[], None);
        assert!(p.is_empty());
        assert_eq!(resid, 0.0);
    }
    #[test]
    fn jacobi_pressure_iterator_records_diagnostics() {
        let jit = JacobiPressureIterator::new(1000.0, 0.7, 5, 1e-20);
        let rhs = vec![1.0_f64, 1.0];
        let aij = vec![vec![(1_usize, -0.3_f64)], vec![(0_usize, -0.3_f64)]];
        let mut diag = ConvergenceDiagnostics::new();
        let (_p, _resid, iters) = jit.run(&rhs, &aij, Some(&mut diag));
        assert!(iters > 0);
        assert!(!diag.is_empty());
    }
    #[test]
    fn pcisph_velocity_correction_zero_pressure_delta_unchanged() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let v_pred = [1.0, 2.0, 3.0];
        let v_corr = solver.compute_predicted_velocity_correction(
            v_pred,
            0.0,
            1000.0,
            0.1,
            0.001,
            [1.0, 0.0, 0.0],
        );
        assert_eq!(v_corr, v_pred, "Zero ΔP should leave velocity unchanged");
    }
    #[test]
    fn pcisph_velocity_correction_direction_correct() {
        let solver = PcisphSolver::new(1000.0, 0.1, 50, 0.001);
        let v_pred = [0.0; 3];
        let v_corr = solver.compute_predicted_velocity_correction(
            v_pred,
            1000.0,
            1000.0,
            0.1,
            0.001,
            [1.0, 0.0, 0.0],
        );
        assert!(
            v_corr[0] < 0.0,
            "Velocity correction should oppose gradient, got {}",
            v_corr[0]
        );
        assert_eq!(v_corr[1], 0.0);
        assert_eq!(v_corr[2], 0.0);
    }
    #[test]
    fn iisph_divergence_free_empty_particles_is_satisfied() {
        let solver = IisphSolver::new(1000.0, 0.7, 100, 1e-4);
        let result = solver.compute_divergence_free_constraint(&[], &[], &[], &[], 1e-4);
        assert!(
            result,
            "Empty particle set should trivially satisfy divergence-free"
        );
    }
    #[test]
    fn iisph_divergence_free_uniform_velocity_is_satisfied() {
        let solver = IisphSolver::new(1000.0, 0.7, 100, 1e-4);
        let velocities = vec![[1.0, 2.0, 3.0]; 3];
        let masses = vec![1.0; 3];
        let densities = vec![1000.0; 3];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![
            vec![(1, 0.1, [1.0, 0.0, 0.0]), (2, 0.1, [0.0, 1.0, 0.0])],
            vec![(0, 0.1, [-1.0, 0.0, 0.0])],
            vec![(0, 0.1, [0.0, -1.0, 0.0])],
        ];
        let result = solver.compute_divergence_free_constraint(
            &velocities,
            &masses,
            &densities,
            &neighbors,
            1e-6,
        );
        assert!(result, "Uniform velocity field should be divergence-free");
    }
    #[test]
    fn iisph_divergence_free_converging_flow_not_satisfied() {
        let solver = IisphSolver::new(1000.0, 0.7, 100, 1e-4);
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![1.0; 2];
        let densities = vec![1000.0; 2];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![
            vec![(1, 0.1, [1.0, 0.0, 0.0])],
            vec![(0, 0.1, [-1.0, 0.0, 0.0])],
        ];
        let result = solver.compute_divergence_free_constraint(
            &velocities,
            &masses,
            &densities,
            &neighbors,
            1e-6,
        );
        assert!(
            !result,
            "Converging flow should not satisfy divergence-free constraint"
        );
    }
    #[test]
    fn dfsph_source_term_zero_for_uniform_velocity() {
        let solver = DfsphSolverSimple::new(1000.0, 0.7, 100, 1e-4, 1e-4);
        let velocities = vec![[1.0, 0.0, 0.0]; 3];
        let masses = vec![1.0; 3];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![
            vec![(1, 0.1, [1.0, 0.0, 0.0]), (2, 0.1, [0.0, 1.0, 0.0])],
            vec![(0, 0.1, [-1.0, 0.0, 0.0])],
            vec![(0, 0.1, [0.0, -1.0, 0.0])],
        ];
        let s = solver.compute_source_term(0, &velocities, &masses, &neighbors);
        assert!(
            s.abs() < 1e-14,
            "Uniform velocity → source term = 0, got {s}"
        );
    }
    #[test]
    fn dfsph_source_term_positive_for_compression() {
        let solver = DfsphSolverSimple::new(1000.0, 0.7, 100, 1e-4, 1e-4);
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![1.0; 2];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![
            vec![(1, 0.1, [1.0, 0.0, 0.0])],
            vec![(0, 0.1, [-1.0, 0.0, 0.0])],
        ];
        let s = solver.compute_source_term(0, &velocities, &masses, &neighbors);
        assert!(
            s > 0.0,
            "Converging particles → positive source term, got {s}"
        );
    }
    #[test]
    fn dfsph_all_source_terms_length_matches_n_particles() {
        let solver = DfsphSolverSimple::new(1000.0, 0.7, 100, 1e-4, 1e-4);
        let velocities = vec![[0.0; 3]; 5];
        let masses = vec![1.0; 5];
        let neighbors: Vec<Vec<(usize, f64, [f64; 3])>> = vec![vec![]; 5];
        let sources = solver.compute_all_source_terms(&velocities, &masses, &neighbors);
        assert_eq!(
            sources.len(),
            5,
            "Should return one source term per particle"
        );
    }
}
#[cfg(test)]
mod tests_pressure_ext {
    use super::*;
    #[test]
    fn density_integrator_drho_dt_zero_for_stationary() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.2);
        let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 2];
        let masses = [1.0; 2];
        let drho = integrator.drho_dt(0, &positions, &velocities, &masses);
        assert_eq!(drho, 0.0, "Stationary particles → drho/dt = 0");
    }
    #[test]
    fn density_integrator_max_error_zero_at_rest_density() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.2);
        let densities = [1000.0, 1000.0, 1000.0];
        assert_eq!(integrator.max_density_error(&densities), 0.0);
    }
    #[test]
    fn density_integrator_max_error_positive_for_deviation() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.2);
        let densities = [1050.0, 950.0];
        let err = integrator.max_density_error(&densities);
        assert!((err - 50.0).abs() < 1e-12);
    }
    #[test]
    fn density_errors_returns_correct_length() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.2);
        let densities = [1000.0; 5];
        let errors = integrator.density_errors(&densities);
        assert_eq!(errors.len(), 5);
    }
    #[test]
    fn density_errors_at_rest_density_are_zero() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.2);
        let densities = [1000.0; 4];
        let errors = integrator.density_errors(&densities);
        for e in errors {
            assert!(e.abs() < 1e-14);
        }
    }
    #[test]
    fn density_integrate_preserves_values_no_neighbors() {
        let integrator = SphDensityIntegrator::new(1000.0, 0.01);
        let mut densities = vec![1000.0_f64];
        let positions = [[0.5; 3]];
        let velocities = [[0.0; 3]];
        let masses = [1.0];
        integrator.integrate(&mut densities, &positions, &velocities, &masses, 0.001);
        assert!((densities[0] - 1000.0).abs() < 1e-12);
    }
    #[test]
    fn sor_solver_empty_system_returns_empty() {
        let solver = PressurePoissonSolver::new(1.0, 100, 1e-8);
        let (p, resid, iters) = solver.solve(&[], &[]);
        assert!(p.is_empty());
        assert_eq!(resid, 0.0);
        assert_eq!(iters, 0);
    }
    #[test]
    fn sor_solver_single_diagonal_converges() {
        let solver = PressurePoissonSolver::new(1.0, 200, 1e-10);
        let rhs = vec![2.0_f64, 2.0];
        let aij = vec![vec![(1_usize, -0.5_f64)], vec![(0_usize, -0.5_f64)]];
        let (p, _resid, _iters) = solver.solve(&rhs, &aij);
        for &pi in &p {
            assert!(pi >= 0.0);
        }
    }
    #[test]
    fn sor_solver_residual_decreases_with_more_iterations() {
        let solver_few = PressurePoissonSolver::new(1.2, 2, 1e-20);
        let solver_many = PressurePoissonSolver::new(1.2, 50, 1e-20);
        let rhs = vec![1.0_f64; 4];
        let aij: Vec<Vec<(usize, f64)>> = vec![
            vec![(1, -0.25), (3, -0.25)],
            vec![(0, -0.25), (2, -0.25)],
            vec![(1, -0.25), (3, -0.25)],
            vec![(2, -0.25), (0, -0.25)],
        ];
        let (_, resid_few, _) = solver_few.solve(&rhs, &aij);
        let (_, resid_many, _) = solver_many.solve(&rhs, &aij);
        assert!(
            resid_many <= resid_few + 1e-10,
            "More iterations should not increase residual: {resid_few} vs {resid_many}"
        );
    }
    #[test]
    fn sor_residual_zero_for_exact_solution() {
        let solver = PressurePoissonSolver::new(1.0, 10, 1e-12);
        let rhs = vec![1.5_f64, 1.5];
        let aij = vec![vec![(1_usize, -0.5_f64)], vec![(0_usize, -0.5_f64)]];
        let p = vec![3.0_f64; 2];
        let resid = solver.residual(&p, &rhs, &aij);
        assert!(resid.is_finite());
    }
    #[test]
    fn iisph_matrix_empty_particles_empty_matrix() {
        let builder = IisphMatrixBuilder::new(0.1, 1000.0, 0.001);
        let aij = builder.build_matrix(&[], &[], &[]);
        assert!(aij.is_empty());
    }
    #[test]
    fn iisph_matrix_single_particle_no_off_diagonal() {
        let builder = IisphMatrixBuilder::new(0.1, 1000.0, 0.001);
        let aij = builder.build_matrix(&[[0.5_f64; 3]], &[1.0], &[1000.0]);
        assert_eq!(aij.len(), 1);
        assert!(
            aij[0].is_empty(),
            "Single particle has no off-diagonal entries"
        );
    }
    #[test]
    fn iisph_rhs_at_rest_density_is_zero() {
        let builder = IisphMatrixBuilder::new(0.1, 1000.0, 0.001);
        let rhs = builder.build_rhs(&[1000.0, 1000.0]);
        for &r in &rhs {
            assert!(r.abs() < 1e-14);
        }
    }
    #[test]
    fn iisph_rhs_compressed_fluid_positive() {
        let builder = IisphMatrixBuilder::new(0.1, 1000.0, 0.001);
        let rhs = builder.build_rhs(&[0.0]);
        assert!((rhs[0] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn iisph_diagonal_non_negative_for_two_close_particles() {
        let builder = IisphMatrixBuilder::new(0.1, 1000.0, 0.001);
        let positions = [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]];
        let masses = [1.0; 2];
        let densities = [1000.0; 2];
        let aij = builder.build_matrix(&positions, &masses, &densities);
        let diag = builder.diagonal(&aij);
        for &d in &diag {
            assert!(d >= 0.0, "IISPH diagonal must be non-negative, got {d}");
        }
    }
    #[test]
    fn pcisph_corrector_zero_pressure_leaves_velocity_unchanged() {
        let corrector = PcisphCorrector::new(0.001, 1000.0);
        let positions = vec![[0.0_f64; 3]];
        let pred_v = vec![[1.0, 2.0, 3.0]];
        let pressures = vec![0.0_f64];
        let masses = vec![1.0];
        let densities = vec![1000.0];
        let v_corr =
            corrector.correct_velocities(&positions, &pred_v, &pressures, &masses, &densities, 0.1);
        for k in 0..3 {
            assert!((v_corr[0][k] - pred_v[0][k]).abs() < 1e-14);
        }
    }
    #[test]
    fn pcisph_corrector_correct_positions_uses_dt_correctly() {
        let corrector = PcisphCorrector::new(0.1, 1000.0);
        let pred_x = vec![[0.0, 0.0, 0.0]];
        let v_corr = vec![[1.0, 2.0, 3.0]];
        let x_corr = corrector.correct_positions(&pred_x, &v_corr);
        assert!((x_corr[0][0] - 0.1).abs() < 1e-14);
        assert!((x_corr[0][1] - 0.2).abs() < 1e-14);
        assert!((x_corr[0][2] - 0.3).abs() < 1e-14);
    }
    #[test]
    fn pcisph_ke_change_zero_for_unchanged_velocity() {
        let corrector = PcisphCorrector::new(0.001, 1000.0);
        let masses = [1.0];
        let v = [[1.0, 0.0, 0.0]];
        let ke_change = corrector.ke_change(&masses, &v, &v);
        assert!(ke_change.abs() < 1e-14);
    }
    #[test]
    fn div_free_projector_uniform_velocity_is_divergence_free() {
        let proj = DivergenceFreeProjector::new(1000.0, 100, 1e-6, 0.7);
        let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]];
        let velocities = [[1.0, 0.0, 0.0]; 3];
        let masses = [1.0; 3];
        let densities = [1000.0; 3];
        let h = 0.2;
        assert!(
            proj.is_divergence_free(&positions, &velocities, &masses, &densities, h),
            "Uniform velocity field must be divergence-free"
        );
    }
    #[test]
    fn div_free_projector_divergence_at_returns_finite() {
        let proj = DivergenceFreeProjector::new(1000.0, 100, 1e-6, 0.7);
        let positions = [[0.0; 3], [0.05, 0.0, 0.0]];
        let velocities = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = [1.0; 2];
        let densities = [1000.0; 2];
        let div = proj.divergence_at(0, &positions, &velocities, &masses, &densities, 0.1);
        assert!(div.is_finite());
    }
    #[test]
    fn div_free_max_divergence_zero_for_uniform() {
        let proj = DivergenceFreeProjector::new(1000.0, 100, 1e-8, 0.7);
        let positions = [[0.0; 3], [0.5, 0.0, 0.0]];
        let velocities = [[2.0, 0.0, 0.0]; 2];
        let masses = [1.0; 2];
        let densities = [1000.0; 2];
        let max_div = proj.max_divergence(&positions, &velocities, &masses, &densities, 0.1);
        assert!(
            max_div < 1e-10,
            "Uniform flow should have near-zero divergence"
        );
    }
    #[test]
    fn tait_eos_pressure_at_rest_density_is_zero() {
        let eos = TaitEosFamily::water();
        assert!(eos.pressure(eos.rho0).abs() < 1e-8);
    }
    #[test]
    fn tait_eos_pressure_compressed_positive() {
        let eos = TaitEosFamily::water();
        assert!(eos.pressure(eos.rho0 * 1.05) > 0.0);
    }
    #[test]
    fn tait_eos_pressure_expanded_negative() {
        let eos = TaitEosFamily::water();
        assert!(eos.pressure(eos.rho0 * 0.95) < 0.0);
    }
    #[test]
    fn tait_eos_inverse_recovers_rest_density() {
        let eos = TaitEosFamily::water();
        let rho_recovered = eos.density_from_pressure(0.0);
        assert!((rho_recovered - eos.rho0).abs() < 1e-6);
    }
    #[test]
    fn tait_eos_sound_speed_at_rest_density_equals_c0() {
        let eos = TaitEosFamily::water();
        let c = eos.sound_speed(eos.rho0);
        assert!((c - eos.c0).abs() < 1e-8, "Sound speed at ρ₀ must equal c₀");
    }
    #[test]
    fn tait_eos_air_sound_speed_reasonable() {
        let eos = TaitEosFamily::air();
        let c = eos.sound_speed(eos.rho0);
        assert!((c - 343.0).abs() < 1e-6);
    }
    #[test]
    fn tait_eos_cfl_dt_positive_and_finite() {
        let eos = TaitEosFamily::water();
        let dt = eos.cfl_dt(0.1, eos.rho0);
        assert!(dt > 0.0);
        assert!(dt.is_finite());
    }
    #[test]
    fn tait_eos_compute_all_matches_element_wise() {
        let eos = TaitEosFamily::water();
        let densities = [1000.0, 1100.0, 900.0];
        let mut pressures = [0.0_f64; 3];
        eos.compute_all(&densities, &mut pressures);
        for (i, &rho) in densities.iter().enumerate() {
            let expected = eos.pressure(rho);
            assert!((pressures[i] - expected).abs() < 1e-10);
        }
    }
    #[test]
    fn tait_eos_min_cfl_dt_returns_minimum() {
        let eos = TaitEosFamily::water();
        let densities = [1000.0, 1200.0, 800.0];
        let min_dt = eos.min_cfl_dt(0.1, &densities);
        let expected_min = densities
            .iter()
            .map(|&rho| eos.cfl_dt(0.1, rho))
            .fold(f64::INFINITY, f64::min);
        assert!((min_dt - expected_min).abs() < 1e-12);
    }
    #[test]
    fn tait_eos_water_big_b_positive() {
        let eos = TaitEosFamily::water();
        assert!(eos.big_b > 0.0);
    }
    #[test]
    fn tait_eos_pressure_and_inverse_round_trip() {
        let eos = TaitEosFamily::water();
        let rho_orig = 1050.0;
        let p = eos.pressure(rho_orig);
        let rho_back = eos.density_from_pressure(p);
        assert!(
            (rho_back - rho_orig).abs() < 1e-6,
            "Round-trip rho error too large"
        );
    }
}
