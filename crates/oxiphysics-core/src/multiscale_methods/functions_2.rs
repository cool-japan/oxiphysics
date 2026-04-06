//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiscale_methods::AdaptiveCriterion;
    use crate::multiscale_methods::Atom;
    use crate::multiscale_methods::CgBead;
    use crate::multiscale_methods::CoarseGrid;
    use crate::multiscale_methods::ConcurrentDomainDecomposition;
    use crate::multiscale_methods::EquationFreeIntegrator;
    use crate::multiscale_methods::HmmConfig;
    use crate::multiscale_methods::MicroDataTable;
    use crate::multiscale_methods::MultilevelErrorEstimator;
    use crate::multiscale_methods::PhaseFieldParams;
    use crate::multiscale_methods::ResolutionLevel;
    use crate::multiscale_methods::RichardsonEstimator;
    use crate::multiscale_methods::SubDomain;
    use crate::multiscale_methods::UnitCell;
    #[test]
    fn test_unit_cell_uniform_voigt() {
        let cell = UnitCell::uniform(4, 4, 4, 0.25, 2.5);
        let v = cell.voigt_average();
        assert!(
            (v - 2.5).abs() < 1e-10,
            "uniform Voigt should equal property"
        );
    }
    #[test]
    fn test_unit_cell_uniform_reuss() {
        let cell = UnitCell::uniform(4, 4, 4, 0.25, 3.0);
        let r = cell.reuss_average();
        assert!((r - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_unit_cell_geometric_average() {
        let cell = UnitCell::uniform(2, 2, 2, 0.5, std::f64::consts::E);
        let g = cell.geometric_average();
        assert!((g - std::f64::consts::E).abs() < 1e-8);
    }
    #[test]
    fn test_voigt_geq_reuss() {
        let cell = UnitCell::checkerboard(4, 4, 4, 0.25, 10.0, 1.0);
        let v = cell.voigt_average();
        let r = cell.reuss_average();
        assert!(v >= r - 1e-10, "Voigt >= Reuss");
    }
    #[test]
    fn test_checkerboard_cell_idx() {
        let cell = UnitCell::checkerboard(3, 3, 1, 1.0, 5.0, 2.0);
        assert_eq!(cell.idx(0, 0, 0), 0);
        assert_eq!(cell.idx(1, 0, 0), 1);
    }
    #[test]
    fn test_homogenize_hs_bounds_ordered() {
        let cell = UnitCell::checkerboard(4, 4, 4, 0.25, 10.0, 1.0);
        let result = homogenize_two_phase(&cell, 10.0, 1.0);
        assert!(
            result.hs_lower <= result.hs_upper + 1e-10,
            "HS lower <= HS upper"
        );
    }
    #[test]
    fn test_homogenize_bounds_contain_effective() {
        let cell = UnitCell::checkerboard(4, 4, 4, 0.25, 10.0, 1.0);
        let result = homogenize_two_phase(&cell, 10.0, 1.0);
        assert!(result.effective_property >= result.reuss_bound - 1e-10);
        assert!(result.effective_property <= result.voigt_bound + 1e-10);
    }
    #[test]
    fn test_homogenize_uniform_trivial() {
        let cell = UnitCell::uniform(4, 4, 4, 0.25, 5.0);
        let result = homogenize_two_phase(&cell, 5.0, 5.0);
        assert!((result.voigt_bound - 5.0).abs() < 1e-8);
    }
    #[test]
    fn test_cell_problem_1d_mean_zero() {
        let cell = UnitCell::uniform(8, 1, 1, 0.125, 1.0);
        let chi = cell_problem_1d(&cell, 2000, 1e-12);
        let mean: f64 = chi.iter().sum::<f64>() / chi.len() as f64;
        assert!(mean.abs() < 0.5, "cell problem solution should be bounded");
    }
    #[test]
    fn test_efi_step_scalar_decay() {
        let efi = EquationFreeIntegrator::new(0.01, 5, 5, 0.1);
        let coarse = vec![1.0f64];
        let result = efi.step(
            &coarse,
            |c| c.to_vec(),
            |fine| fine.iter().map(|&x| x * 0.99).collect(),
            |fine| fine.to_vec(),
        );
        assert!(result[0] < 1.0, "decaying system should decrease");
        assert!(result[0] > 0.0, "should remain positive");
    }
    #[test]
    fn test_efi_step_preserves_length() {
        let efi = EquationFreeIntegrator::new(0.01, 3, 3, 0.05);
        let coarse = vec![1.0, 2.0, 3.0];
        let result = efi.step(
            &coarse,
            |c| c.to_vec(),
            |f| f.iter().map(|&x| x * 0.999).collect(),
            |f| f.to_vec(),
        );
        assert_eq!(result.len(), 3);
    }
    #[test]
    fn test_hmm_flux_positive_for_positive_strain() {
        let config = HmmConfig {
            n_macro_points: 4,
            micro_steps: 10,
            dt_micro: 0.01,
            dt_macro: 0.1,
            micro_domain_size: 1.0,
        };
        let flux = hmm_estimate_flux(&config, 1.0, |_x| 1.0);
        assert!(flux.is_finite());
    }
    #[test]
    fn test_hmm_flux_zero_for_zero_strain() {
        let config = HmmConfig {
            n_macro_points: 4,
            micro_steps: 10,
            dt_micro: 0.01,
            dt_macro: 0.1,
            micro_domain_size: 1.0,
        };
        let flux = hmm_estimate_flux(&config, 0.0, |_x| 1.0);
        assert!(flux.abs() < 1e-6, "zero strain => near-zero flux");
    }
    #[test]
    fn test_subdomain_contains() {
        let dom = SubDomain::new(0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], false);
        assert!(dom.contains([0.5, 0.5, 0.5]));
        assert!(!dom.contains([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_subdomain_volume() {
        let dom = SubDomain::new(0, [0.0, 0.0, 0.0], [2.0, 3.0, 4.0], false);
        assert!((dom.volume() - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_subdomain_overlap_detection() {
        let d0 = SubDomain::new(0, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0], false);
        let d1 = SubDomain::new(1, [1.0, 1.0, 1.0], [3.0, 3.0, 3.0], true);
        let d2 = SubDomain::new(2, [5.0, 5.0, 5.0], [6.0, 6.0, 6.0], false);
        assert!(d0.overlaps(&d1), "d0 and d1 should overlap");
        assert!(!d0.overlaps(&d2), "d0 and d2 should not overlap");
    }
    #[test]
    fn test_domain_decomp_interpolation() {
        let mut dd = ConcurrentDomainDecomposition::new();
        dd.add_domain(SubDomain::new(0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], false));
        dd.add_domain(SubDomain::new(1, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0], false));
        let field = [10.0, 20.0];
        let val = dd.interpolate_field([0.5, 0.5, 0.5], &field);
        assert!(val.is_finite());
        assert!(val > 0.0);
    }
    #[test]
    fn test_coarse_grid_average_values() {
        let fine = vec![1.0, 3.0, 5.0, 7.0];
        let cg = CoarseGrid::from_fine(fine, 2);
        assert!((cg.coarse_values[0] - 2.0).abs() < 1e-10);
        assert!((cg.coarse_values[1] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_downscale_length() {
        let fine = vec![1.0f64; 8];
        let cg = CoarseGrid::from_fine(fine.clone(), 2);
        let d = cg.downscale();
        assert_eq!(d.len(), fine.len());
    }
    #[test]
    fn test_downscale_linear_length() {
        let fine = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let cg = CoarseGrid::from_fine(fine.clone(), 2);
        let d = cg.downscale_linear();
        assert_eq!(d.len(), fine.len());
    }
    #[test]
    fn test_lj_potential_at_min() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let u = lennard_jones_potential(r_min, epsilon, sigma);
        assert!(
            (u - (-epsilon)).abs() < 1e-8,
            "LJ minimum should be -epsilon"
        );
    }
    #[test]
    fn test_lj_force_zero_at_min() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let f = lennard_jones_force(r_min, epsilon, sigma);
        assert!(f.abs() < 1e-7, "LJ force at minimum should be ~0");
    }
    #[test]
    fn test_atom_kinetic_energy() {
        let atom = Atom::new([0.0; 3], [1.0, 0.0, 0.0], 2.0);
        assert!((atom.kinetic_energy() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_atomistic_temperature_positive() {
        let atoms = vec![
            Atom::new([0.0; 3], [1.0, 0.0, 0.0], 1.0),
            Atom::new([1.0; 3], [0.0, 1.0, 0.0], 1.0),
        ];
        let t = atomistic_temperature(&atoms, 1.0);
        assert!(t > 0.0);
    }
    #[test]
    fn test_virial_stress_tensor_shape() {
        let atoms = vec![
            Atom::new([0.0, 0.0, 0.0], [0.0; 3], 1.0),
            Atom::new([1.5, 0.0, 0.0], [0.0; 3], 1.0),
        ];
        let pairs = vec![(0, 1)];
        let sigma = virial_stress_tensor(&atoms, &pairs, 1.0, 1.0, 1.0);
        assert_eq!(sigma.len(), 9);
    }
    #[test]
    fn test_richardson_extrapolation() {
        let re = RichardsonEstimator::new(2.0);
        let u_c = 1.04;
        let u_f = 1.01;
        let extrap = re.extrapolate(u_c, u_f);
        assert!((extrap - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_richardson_error_estimate() {
        let re = RichardsonEstimator::new(2.0);
        let err = re.error_estimate(1.04, 1.01);
        assert!((err - 0.01).abs() < 1e-10);
    }
    #[test]
    fn test_multilevel_estimator_convergence() {
        let mut est = MultilevelErrorEstimator::new(3, 1e-3);
        est.record(0, 1e-5);
        est.record(1, 1e-5);
        est.record(2, 1e-5);
        assert!(est.all_converged());
    }
    #[test]
    fn test_multilevel_not_converged() {
        let mut est = MultilevelErrorEstimator::new(2, 1e-6);
        est.record(0, 1.0);
        assert!(!est.all_converged());
    }
    #[test]
    fn test_macro_diffusion_boundary_values() {
        let n = 10;
        let k_eff = vec![1.0f64; n];
        let source = vec![0.0f64; n];
        let result = solve_macro_diffusion_1d(n, 1.0, &k_eff, &source, 0.0, 1.0);
        assert!((result[0] - 0.0).abs() < 1e-10);
        assert!((result[n + 1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_macro_diffusion_linear_solution() {
        let n = 10;
        let k_eff = vec![1.0f64; n];
        let source = vec![0.0f64; n];
        let result = solve_macro_diffusion_1d(n, 1.0, &k_eff, &source, 0.0, 1.0);
        for i in 1..=n {
            let x = i as f64 / (n + 1) as f64;
            assert!(
                (result[i] - x).abs() < 1e-6,
                "linear solution expected at node {}",
                i
            );
        }
    }
    #[test]
    fn test_restriction_l2_sum_preserved() {
        let fine = vec![1.0f64, 3.0, 5.0, 7.0];
        let coarse = restriction_l2(&fine, 2);
        assert_eq!(coarse.len(), 2);
        assert!((coarse[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_prolongation_linear_length() {
        let coarse = vec![0.0f64, 1.0];
        let fine = prolongation_linear(&coarse, 4);
        assert_eq!(fine.len(), 8);
    }
    #[test]
    fn test_prolongation_preserves_values_at_centres() {
        let coarse = vec![0.0f64, 2.0];
        let fine = prolongation_linear(&coarse, 2);
        assert!((fine[0] - 0.0).abs() < 1e-10);
        assert!((fine[fine.len() - 1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_allen_cahn_step_bounds() {
        let params = PhaseFieldParams::new(0.01, 1.0, 1.0);
        let mut phi = vec![0.2, 0.5, 0.8, 0.5, 0.2, 0.5, 0.8, 0.5];
        allen_cahn_step(&mut phi, &params, 0.1, 0.001);
        for &p in &phi {
            assert!((0.0..=1.0).contains(&p), "phi out of [0,1]: {}", p);
        }
    }
    #[test]
    fn test_phase_field_params_interface_width() {
        let params = PhaseFieldParams::new(0.04, 1.0, 2.0);
        let w = params.interface_width();
        assert!(w > 0.0, "interface width should be positive");
    }
    #[test]
    fn test_scale_separation_pass() {
        assert!(check_scale_separation(0.001, 1.0, 100.0));
    }
    #[test]
    fn test_scale_separation_fail() {
        assert!(!check_scale_separation(0.1, 1.0, 100.0));
    }
    #[test]
    fn test_table_interpolation_boundary() {
        let table = MicroDataTable::new(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 4.0]);
        assert!((table.query(0.0) - 0.0).abs() < 1e-10);
        assert!((table.query(2.0) - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_table_interpolation_midpoint() {
        let table = MicroDataTable::new(vec![0.0, 2.0], vec![0.0, 4.0]);
        let mid = table.query(1.0);
        assert!(
            (mid - 2.0).abs() < 1e-10,
            "linear interpolation at midpoint"
        );
    }
    #[test]
    fn test_table_build() {
        let table = MicroDataTable::build(0.0, 1.0, 5, |x| x * x);
        assert_eq!(table.inputs.len(), 5);
        assert!((table.outputs[0] - 0.0).abs() < 1e-10);
        assert!((table.outputs[4] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_adaptive_refine_trigger() {
        let crit = AdaptiveCriterion::new(0.01, 0.0001, 0.5);
        let level = crit.decide(0.1, 0.0, ResolutionLevel::Coarse);
        assert_eq!(level, ResolutionLevel::Meso);
    }
    #[test]
    fn test_adaptive_coarsen_trigger() {
        let crit = AdaptiveCriterion::new(0.01, 0.001, 0.5);
        let level = crit.decide(1e-5, 1e-5, ResolutionLevel::Fine);
        assert_eq!(level, ResolutionLevel::Meso);
    }
    #[test]
    fn test_adaptive_no_change() {
        let crit = AdaptiveCriterion::new(0.1, 0.001, 1.0);
        let level = crit.decide(0.05, 0.0, ResolutionLevel::Meso);
        assert_eq!(level, ResolutionLevel::Meso);
    }
    #[test]
    fn test_cg_bond_force_magnitude() {
        let pa = [0.0, 0.0, 0.0];
        let pb = [2.0, 0.0, 0.0];
        let f = cg_bond_force(pa, pb, 10.0, 1.5);
        let mag = norm3(f);
        assert!(
            (mag - 10.0 * 0.5).abs() < 1e-8,
            "bond force magnitude should be k*(r-r0)"
        );
    }
    #[test]
    fn test_cg_bead_kinetic_energy_zero() {
        let bead = CgBead::new([0.0; 3], 1.0);
        let ke = 0.5 * bead.mass * dot3(bead.vel, bead.vel);
        assert!((ke - 0.0).abs() < 1e-14);
    }
    #[test]
    fn test_multigrid_vcycle_convergence() {
        let n = 17usize;
        let h = 1.0 / (n - 1) as f64;
        let f = vec![0.0f64; n];
        let mut u = vec![0.0f64; n];
        u[0] = 0.0;
        u[n - 1] = 1.0;
        multigrid_vcycle(&mut u, &f, h, 3, 3);
        assert!(u[n / 2].abs() < 2.0, "multigrid solution should be bounded");
    }
    #[test]
    fn test_arlequin_blend_extremes() {
        assert!((arlequin_blend_energy(10.0, 2.0, 1.0) - 10.0).abs() < 1e-10);
        assert!((arlequin_blend_energy(10.0, 2.0, 0.0) - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_arlequin_blend_midpoint() {
        let blend = arlequin_blend_energy(4.0, 2.0, 0.5);
        assert!((blend - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_hill_mandel_zero_for_consistent() {
        let s = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let e = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let ms = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let me = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let err = hill_mandel_error(&s, &e, ms, me);
        assert!(err.is_finite());
    }
    #[test]
    fn test_msd_zero_for_same_positions() {
        let atoms = vec![Atom::new([1.0, 2.0, 3.0], [0.0; 3], 1.0)];
        let refs = vec![[1.0, 2.0, 3.0]];
        let msd = mean_squared_displacement(&atoms, &refs);
        assert!(msd.abs() < 1e-14);
    }
    #[test]
    fn test_msd_unit_displacement() {
        let atoms = vec![Atom::new([1.0, 0.0, 0.0], [0.0; 3], 1.0)];
        let refs = vec![[0.0, 0.0, 0.0]];
        let msd = mean_squared_displacement(&atoms, &refs);
        assert!((msd - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cauchy_born_stress_zero_for_identity() {
        let sigma = 1.0;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let s = cauchy_born_stress_1d(1.0, 1.0, sigma, r_min);
        assert!(
            s.abs() < 1e-5,
            "stress at equilibrium should be ~0, got {}",
            s
        );
    }
    #[test]
    fn test_cauchy_born_modulus_positive() {
        let modulus = cauchy_born_modulus(&[10.0, 10.0], &[1.0, 1.0], 2.0);
        assert!(modulus > 0.0);
    }
}
