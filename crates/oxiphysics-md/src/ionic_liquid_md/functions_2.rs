//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::ionic_liquid_md::IlSimParams;
    use crate::ionic_liquid_md::IonModel;
    use crate::ionic_liquid_md::IonType;
    use crate::ionic_liquid_md::IonicLiquidSimulation;
    use crate::ionic_liquid_md::PairType;
    use crate::ionic_liquid_md::RdfResult;
    use crate::ionic_liquid_md::WolfParams;
    use std::f64::consts::PI;
    fn make_two_ion_state() -> IonState {
        let cat = IonModel::emim_cation();
        let an = IonModel::bf4_anion();
        let mut state = IonState::new();
        state.add_ion(&cat, [5.0, 5.0, 5.0], [0.1, 0.0, 0.0]);
        state.add_ion(&an, [8.0, 5.0, 5.0], [-0.1, 0.0, 0.0]);
        state
    }
    fn make_small_system() -> IonState {
        let cat = IonModel::emim_cation();
        let an = IonModel::bf4_anion();
        init_cubic_lattice(&cat, &an, 4, 20.0)
    }
    #[test]
    fn test_ion_model_creation() {
        let cat = IonModel::emim_cation();
        assert_eq!(cat.ion_type, IonType::Cation);
        assert!(cat.charge > 0.0);
        assert!(cat.mass > 0.0);
        let an = IonModel::bf4_anion();
        assert_eq!(an.ion_type, IonType::Anion);
        assert!(an.charge < 0.0);
    }
    #[test]
    fn test_ion_state_basic() {
        let state = make_two_ion_state();
        assert_eq!(state.n, 2);
        assert_eq!(state.n_cations(), 1);
        assert_eq!(state.n_anions(), 1);
    }
    #[test]
    fn test_total_charge_neutral() {
        let state = make_two_ion_state();
        let total_q = state.total_charge();
        assert!((total_q).abs() < 1e-10, "System should be charge-neutral");
    }
    #[test]
    fn test_kinetic_energy_positive() {
        let state = make_two_ion_state();
        assert!(state.kinetic_energy() > 0.0);
    }
    #[test]
    fn test_temperature_positive() {
        let state = make_two_ion_state();
        assert!(state.temperature() > 0.0);
    }
    #[test]
    fn test_center_of_mass() {
        let state = make_two_ion_state();
        let com = state.center_of_mass();
        assert!(com[0] > 5.0 && com[0] < 8.0);
    }
    #[test]
    fn test_bmh_energy_repulsive_at_short_range() {
        let bmh = BornMayerHugginsParams::nacl_like();
        let e_short = bmh.energy(1.5);
        let e_long = bmh.energy(5.0);
        assert!(
            e_short > e_long,
            "BMH should be more repulsive at short range"
        );
    }
    #[test]
    fn test_bmh_force_positive_short_range() {
        let bmh = BornMayerHugginsParams::nacl_like();
        let f = bmh.force(1.5);
        assert!(f > 0.0, "BMH force should be repulsive at short range");
    }
    #[test]
    fn test_buckingham_energy() {
        let buck = BuckinghamParams::generic_ionic();
        let e = buck.energy(3.0);
        assert!(e.is_finite());
    }
    #[test]
    fn test_buckingham_repulsive_short() {
        let buck = BuckinghamParams::generic_ionic();
        let e_short = buck.energy(1.0);
        let e_mid = buck.energy(3.0);
        assert!(
            e_short > e_mid,
            "e(1.0)={:.4e} should exceed e(3.0)={:.4e}",
            e_short,
            e_mid
        );
    }
    #[test]
    fn test_wolf_opposite_charges_attractive() {
        let wolf = WolfParams::default_il();
        let e = wolf.energy(1.0, -1.0, 4.0);
        assert!(e < 0.0, "Opposite charges should attract: e={}", e);
    }
    #[test]
    fn test_wolf_same_charges_repulsive() {
        let wolf = WolfParams::default_il();
        let e = wolf.energy(1.0, 1.0, 4.0);
        assert!(e > 0.0, "Same charges should repel: e={}", e);
    }
    #[test]
    fn test_wolf_beyond_cutoff_zero() {
        let wolf = WolfParams::default_il();
        let e = wolf.energy(1.0, -1.0, 15.0);
        assert!((e).abs() < 1e-20, "Beyond cutoff should be zero");
    }
    #[test]
    fn test_erfc_at_zero() {
        let val = erfc_approx(0.0);
        assert!((val - 1.0).abs() < 1e-5, "erfc(0) ~ 1.0, got {}", val);
    }
    #[test]
    fn test_erfc_large_argument() {
        let val = erfc_approx(5.0);
        assert!(val < 1e-6, "erfc(5) should be very small, got {}", val);
    }
    #[test]
    fn test_mixing_rules() {
        assert!((mix_sigma(4.0, 6.0) - 5.0).abs() < 1e-10);
        assert!((mix_epsilon(4.0, 9.0) - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_lj_energy_minimum() {
        let sigma = 3.5;
        let eps = 1.0;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let e_min = lj_energy(r_min, sigma, eps);
        assert!(
            (e_min - (-eps)).abs() < 1e-10,
            "LJ minimum should be -eps, got {}",
            e_min
        );
    }
    #[test]
    fn test_force_computation_runs() {
        let mut state = make_two_ion_state();
        let params = IlSimParams::default_params();
        compute_forces(&mut state, &params);
        let f_norm = norm3(state.forces[0]);
        assert!(f_norm > 0.0, "Forces should be nonzero");
    }
    #[test]
    fn test_newton_third_law() {
        let mut state = make_two_ion_state();
        let params = IlSimParams::default_params();
        compute_forces(&mut state, &params);
        for d in 0..3 {
            let diff = (state.forces[0][d] + state.forces[1][d]).abs();
            assert!(
                diff < 1e-10,
                "Newton's 3rd law violated on axis {}: diff={}",
                d,
                diff
            );
        }
    }
    #[test]
    fn test_potential_energy_finite() {
        let state = make_two_ion_state();
        let params = IlSimParams::default_params();
        let pe = compute_potential_energy(&state, &params);
        assert!(pe.is_finite());
    }
    #[test]
    fn test_velocity_verlet_conserves_particles() {
        let mut state = make_two_ion_state();
        let params = IlSimParams::default_params();
        compute_forces(&mut state, &params);
        let n_before = state.n;
        velocity_verlet_step(&mut state, &params);
        assert_eq!(state.n, n_before);
    }
    #[test]
    fn test_rdf_cation_anion() {
        let state = make_small_system();
        let params = IlSimParams {
            box_len: 20.0,
            ..IlSimParams::default_params()
        };
        let rdf = compute_rdf(&state, &params, PairType::CationAnion, 50, 10.0);
        assert_eq!(rdf.n_bins, 50);
        assert_eq!(rdf.r.len(), 50);
        assert_eq!(rdf.g_r.len(), 50);
    }
    #[test]
    fn test_rdf_all_pairs() {
        let state = make_small_system();
        let params = IlSimParams {
            box_len: 20.0,
            ..IlSimParams::default_params()
        };
        let rdf = compute_rdf(&state, &params, PairType::All, 50, 15.0);
        let sum: f64 = rdf.g_r.iter().sum();
        assert!(sum > 0.0, "RDF should have nonzero values, sum={:.6e}", sum);
    }
    #[test]
    fn test_coordination_number() {
        let rdf = RdfResult {
            r: vec![0.5, 1.5, 2.5, 3.5, 4.5],
            g_r: vec![0.0, 2.0, 1.5, 1.0, 1.0],
            n_bins: 5,
            dr: 1.0,
        };
        let cn = coordination_number(&rdf, 0.01);
        assert_eq!(cn.len(), 5);
        for i in 1..cn.len() {
            assert!(cn[i] >= cn[i - 1]);
        }
    }
    #[test]
    fn test_voronoi_coordination_basic() {
        let state = make_small_system();
        let counts = voronoi_coordination(&state, 20.0, IonType::Cation, IonType::Anion, 8.0);
        assert!(!counts.is_empty());
    }
    #[test]
    fn test_msd_zero_displacement() {
        let pos = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let types = vec![IonType::Cation, IonType::Anion];
        let msd = compute_msd(&pos, &pos, &types, None);
        assert!((msd).abs() < 1e-10, "MSD should be zero for same positions");
    }
    #[test]
    fn test_msd_known_displacement() {
        let pos0 = vec![[0.0, 0.0, 0.0]];
        let pos1 = vec![[3.0, 4.0, 0.0]];
        let types = vec![IonType::Cation];
        let msd = compute_msd(&pos0, &pos1, &types, None);
        assert!((msd - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_diffusion_from_msd() {
        let d = diffusion_from_msd(6.0, 1.0);
        assert!((d - 1.0).abs() < 1e-10, "D should be 1.0 for MSD=6 at t=1");
    }
    #[test]
    fn test_electric_current() {
        let charges = vec![1.0, -1.0];
        let vels = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let j = electric_current(&charges, &vels);
        assert!((j[0]).abs() < 1e-10);
    }
    #[test]
    fn test_current_autocorrelation_symmetric() {
        let traj = vec![[1.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let cac = current_autocorrelation(&traj, 2);
        assert_eq!(cac.len(), 2);
        assert!(cac[0] > 0.0);
    }
    #[test]
    fn test_nernst_einstein_positive() {
        let sigma = nernst_einstein_conductivity(10, 10, 1.0, -1.0, 1e-9, 1e-9, 1e-27, 300.0);
        assert!(sigma > 0.0);
    }
    #[test]
    fn test_haven_ratio_unity() {
        let hr = haven_ratio(5.0, 5.0);
        assert!((hr - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_walden_plot_generation() {
        let conds = vec![0.1, 0.5, 1.0];
        let viscs = vec![10.0, 5.0, 2.0];
        let temps = vec![300.0, 350.0, 400.0];
        let points = walden_plot(&conds, &viscs, &temps, 200.0, 1.2);
        assert_eq!(points.len(), 3);
    }
    #[test]
    fn test_walden_classification() {
        assert_eq!(walden_classification(1.0, 1.0), "good ionic liquid");
        assert_eq!(walden_classification(0.0, 1.0), "non-ionic / associated");
    }
    #[test]
    fn test_cage_radius_small_system() {
        let state = make_small_system();
        let r = cage_radius(&state, 0, 20.0, 3);
        assert!(r > 0.0, "Cage radius should be positive");
    }
    #[test]
    fn test_vft_viscosity() {
        let eta = vft_viscosity(0.01, 500.0, 200.0, 400.0);
        assert!(eta > 0.0 && eta.is_finite());
    }
    #[test]
    fn test_fragility_index() {
        let m = fragility_index(500.0, 300.0, 200.0);
        assert!(m > 0.0);
    }
    #[test]
    fn test_fragility_classification() {
        assert_eq!(fragility_class(20.0), "strong");
        assert_eq!(fragility_class(100.0), "fragile");
    }
    #[test]
    fn test_glass_transition_detection() {
        let temps: Vec<f64> = (0..20).map(|i| 200.0 + i as f64 * 20.0).collect();
        let densities: Vec<f64> = temps
            .iter()
            .map(|&t| {
                if t < 350.0 {
                    1.3 - 0.0002 * t
                } else {
                    1.3 - 0.0002 * 350.0 - 0.0008 * (t - 350.0)
                }
            })
            .collect();
        let tg = detect_glass_transition(&temps, &densities);
        assert!(tg.is_some(), "Should detect glass transition");
        let tg_val = tg.unwrap();
        assert!(
            (tg_val - 350.0).abs() < 50.0,
            "T_g should be near 350, got {:.1}",
            tg_val
        );
    }
    #[test]
    fn test_non_gaussian_parameter_gaussian() {
        let displacements: Vec<[f64; 3]> = (0..1000)
            .map(|i| {
                let x = (i as f64 / 1000.0 - 0.5) * 2.0;
                [x, 0.0, 0.0]
            })
            .collect();
        let alpha2 = non_gaussian_parameter(&displacements);
        assert!(alpha2.is_finite());
    }
    #[test]
    fn test_init_cubic_lattice() {
        let cat = IonModel::emim_cation();
        let an = IonModel::bf4_anion();
        let state = init_cubic_lattice(&cat, &an, 4, 20.0);
        assert_eq!(state.n, 8);
        assert_eq!(state.n_cations(), 4);
        assert_eq!(state.n_anions(), 4);
    }
    #[test]
    fn test_assign_velocities() {
        let mut state = make_small_system();
        assign_velocities(&mut state, 300.0);
        let mut v_com = [0.0; 3];
        let mut total_mass = 0.0;
        for i in 0..state.n {
            v_com = add3(v_com, scale3(state.velocities[i], state.masses[i]));
            total_mass += state.masses[i];
        }
        v_com = scale3(v_com, 1.0 / total_mass);
        let v_com_mag = norm3(v_com);
        assert!(
            v_com_mag < 1e-10,
            "COM velocity should be zero, got {:.6e}",
            v_com_mag
        );
    }
    #[test]
    fn test_compute_density() {
        let state = make_small_system();
        let rho = compute_density(&state, 20.0);
        assert!(
            rho > 0.0 && rho < 100.0,
            "Density should be reasonable, got {:.4}",
            rho
        );
    }
    #[test]
    fn test_charge_density_profile() {
        let state = make_small_system();
        let (pos, prof) = charge_density_profile(&state, 20.0, 2, 10);
        assert_eq!(pos.len(), 10);
        assert_eq!(prof.len(), 10);
    }
    #[test]
    fn test_number_density_profile() {
        let state = make_small_system();
        let (pos, prof) = number_density_profile(&state, 20.0, 0, 10, IonType::Cation);
        assert_eq!(pos.len(), 10);
        let total: f64 = prof.iter().sum::<f64>() * (20.0 * 20.0 * 2.0);
        assert!(total > 0.0);
    }
    #[test]
    fn test_simulation_driver() {
        let state = make_small_system();
        let params = IlSimParams {
            box_len: 20.0,
            ..IlSimParams::default_params()
        };
        let mut sim = IonicLiquidSimulation::new(state, params);
        sim.run(10, 5);
        assert!(sim.step == 10);
        assert!(sim.time > 0.0);
        assert!(!sim.trajectory.is_empty());
    }
    #[test]
    fn test_linear_fit() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let (a, b) = linear_fit(&x, &y);
        assert!((a - 2.0).abs() < 1e-10, "slope should be 2, got {:.6}", a);
        assert!(b.abs() < 1e-10, "intercept should be 0, got {:.6}", b);
    }
    #[test]
    fn test_stress_tensor_element() {
        let state = make_two_ion_state();
        let vol = 30.0_f64.powi(3);
        let sxy = stress_tensor_element(&state, vol, 0, 1);
        assert!(sxy.is_finite());
    }
    #[test]
    fn test_min_image() {
        let dr = [16.0, -14.0, 5.0];
        let wrapped = min_image(dr, 20.0);
        for d in 0..3 {
            assert!(
                wrapped[d].abs() <= 10.0,
                "Wrapped component should be within half box"
            );
        }
    }
    #[test]
    fn test_walden_product() {
        let wp = walden_product(0.1, 5.0);
        assert!((wp - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_structure_factor() {
        let positions = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let types = vec![IonType::Cation, IonType::Cation, IonType::Cation];
        let sq = structure_factor(&positions, &types, IonType::Cation, 2.0 * PI);
        assert!(sq.is_finite() && sq >= 0.0);
    }
    #[test]
    fn test_pf6_and_lithium_models() {
        let pf6 = IonModel::pf6_anion();
        assert_eq!(pf6.ion_type, IonType::Anion);
        assert!(pf6.mass > 100.0);
        let li = IonModel::lithium_cation();
        assert_eq!(li.ion_type, IonType::Cation);
        assert!(li.mass < 10.0);
    }
    #[test]
    fn test_ion_pair_lifetime_trivial() {
        let traj = vec![
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        ];
        let types = vec![IonType::Cation, IonType::Anion];
        let c_t = ion_pair_lifetime(&traj, &types, 20.0, 5.0, 2);
        assert_eq!(c_t.len(), 2);
        assert!((c_t[0] - 1.0).abs() < 1e-10, "C(0) should be 1.0");
    }
}
