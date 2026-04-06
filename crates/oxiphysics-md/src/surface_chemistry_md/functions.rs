//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Boltzmann constant (J/K).
pub(super) const KB: f64 = 1.380_649e-23;
/// Avogadro's number (mol⁻¹).
pub(super) const NA: f64 = 6.022_140_76e23;
#[cfg(test)]
mod tests {

    use crate::surface_chemistry_md::AdsorptionModel;
    use crate::surface_chemistry_md::CatalysisMd;
    use crate::surface_chemistry_md::InterfaceMd;
    use crate::surface_chemistry_md::InterfaceMolecule;
    use crate::surface_chemistry_md::IsothermType;
    use crate::surface_chemistry_md::PrecursorState;
    use crate::surface_chemistry_md::ReactionStep;
    use crate::surface_chemistry_md::SurfaceAtom;
    use crate::surface_chemistry_md::SurfaceDiffusion;
    use crate::surface_chemistry_md::SurfaceReconstruction;
    use crate::surface_chemistry_md::TriboFilm;
    use crate::surface_chemistry_md::Tribochemistry;
    use crate::surface_chemistry_md::WearParticle;
    fn default_model() -> AdsorptionModel {
        AdsorptionModel::new(1e-4, 1.0, 0.5, 2.0, 8000.0, 1e-6, 100.0, 1e5, 1.0, 300.0)
    }
    #[test]
    fn test_langmuir_zero_pressure() {
        let m = default_model();
        assert_eq!(m.langmuir(0.0), 0.0);
    }
    #[test]
    fn test_langmuir_half_coverage() {
        let m = AdsorptionModel::new(2.0, 1.0, 0.5, 2.0, 8000.0, 1e-6, 100.0, 1e5, 1.0, 300.0);
        let p = 1.0 / 2.0;
        let theta = m.langmuir(p);
        assert!((theta - 0.5).abs() < 1e-10, "theta={theta}");
    }
    #[test]
    fn test_langmuir_high_pressure() {
        let m = default_model();
        let theta = m.langmuir(1e12);
        assert!(theta > 0.99, "theta={theta}");
    }
    #[test]
    fn test_langmuir_monotone() {
        let m = default_model();
        let t1 = m.langmuir(1e3);
        let t2 = m.langmuir(1e6);
        assert!(t2 > t1);
    }
    #[test]
    fn test_freundlich_zero_pressure() {
        let m = default_model();
        assert_eq!(m.freundlich(0.0), 0.0);
    }
    #[test]
    fn test_freundlich_positive() {
        let m = default_model();
        assert!(m.freundlich(1e4) > 0.0);
    }
    #[test]
    fn test_freundlich_n1_linear() {
        let m = AdsorptionModel::new(1e-4, 1.0, 2.0, 1.0, 8000.0, 1e-6, 100.0, 1e5, 1.0, 300.0);
        let q = m.freundlich(3.0);
        assert!((q - 6.0).abs() < 1e-10, "q={q}");
    }
    #[test]
    fn test_freundlich_sublinear_n2() {
        let m = default_model();
        let q1 = m.freundlich(1.0);
        let q4 = m.freundlich(4.0);
        assert!((q4 / q1 - 2.0).abs() < 1e-6, "ratio={}", q4 / q1);
    }
    #[test]
    fn test_temkin_zero_pressure() {
        let m = default_model();
        assert_eq!(m.temkin(0.0), 0.0);
    }
    #[test]
    fn test_temkin_positive_pressure() {
        let m = AdsorptionModel::new(1e-4, 1.0, 0.5, 2.0, 8000.0, 1e-3, 100.0, 1e5, 1.0, 300.0);
        let q = m.temkin(1e4);
        assert!(q > 0.0, "q={q}");
    }
    #[test]
    fn test_temkin_increases_with_pressure() {
        let m = AdsorptionModel::new(1e-4, 1.0, 0.5, 2.0, 8000.0, 1e-3, 100.0, 1e5, 1.0, 300.0);
        let q1 = m.temkin(1e3);
        let q2 = m.temkin(1e5);
        assert!(q2 > q1, "q1={q1} q2={q2}");
    }
    #[test]
    fn test_bet_zero_pressure() {
        let m = default_model();
        assert_eq!(m.bet(0.0), 0.0);
    }
    #[test]
    fn test_bet_at_saturation() {
        let m = default_model();
        assert_eq!(m.bet(1e5), 0.0);
    }
    #[test]
    fn test_bet_below_saturation() {
        let m = default_model();
        let v = m.bet(0.5e5);
        assert!(v > 0.0, "v={v}");
    }
    #[test]
    fn test_bet_monotone() {
        let m = default_model();
        let v1 = m.bet(0.1e5);
        let v2 = m.bet(0.8e5);
        assert!(v2 > v1, "v1={v1} v2={v2}");
    }
    #[test]
    fn test_coverage_dispatch() {
        let m = default_model();
        let p = 5e4;
        assert_eq!(m.coverage(p, IsothermType::Langmuir), m.langmuir(p));
        assert_eq!(m.coverage(p, IsothermType::Freundlich), m.freundlich(p));
        assert_eq!(m.coverage(p, IsothermType::Bet), m.bet(p));
    }
    #[test]
    fn test_henry_constant_positive() {
        let m = default_model();
        assert!(m.henry_constant() > 0.0);
    }
    fn make_sd() -> SurfaceDiffusion {
        SurfaceDiffusion::new(
            vec![0.0, 1.0, 2.0],
            vec![1e-21, 5e-21, 1e-20],
            300.0,
            1e13,
            3e-10,
        )
    }
    #[test]
    fn test_sd_hop_rate_positive() {
        let sd = make_sd();
        assert!(sd.hop_rate(1e-21) > 0.0);
    }
    #[test]
    fn test_sd_hop_rate_decreases_with_barrier() {
        let sd = make_sd();
        let r1 = sd.hop_rate(1e-21);
        let r2 = sd.hop_rate(1e-20);
        assert!(r1 > r2);
    }
    #[test]
    fn test_sd_msd_zero_at_start() {
        let sd = make_sd();
        let init = vec![0.0, 1.0, 2.0];
        assert_eq!(sd.mean_square_displacement(&init), 0.0);
    }
    #[test]
    fn test_sd_msd_after_displacement() {
        let mut sd = make_sd();
        sd.positions = vec![1.0, 2.0, 3.0];
        let init = vec![0.0, 0.0, 0.0];
        let msd = sd.mean_square_displacement(&init);
        assert!((msd - 14.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_sd_mc_step_runs() {
        let mut sd = make_sd();
        let h = sd.monte_carlo_step(1e-15, 42);
        assert!(h <= 3);
    }
    #[test]
    fn test_sd_mc_step_empty() {
        let mut sd = SurfaceDiffusion::new(vec![], vec![], 300.0, 1e13, 3e-10);
        assert_eq!(sd.monte_carlo_step(1e-12, 1), 0);
    }
    #[test]
    fn test_sd_precursor_add() {
        let mut sd = make_sd();
        sd.add_precursor(1.5, -5e-21);
        assert_eq!(sd.precursors.len(), 1);
    }
    #[test]
    fn test_sd_chemisorption_fraction_empty() {
        let sd = make_sd();
        assert_eq!(sd.chemisorption_fraction(), 0.0);
    }
    #[test]
    fn test_sd_precursor_step_runs() {
        let mut sd = make_sd();
        sd.add_precursor(1.0, -5e-21);
        sd.add_precursor(2.0, -3e-21);
        let n = sd.precursor_step(1e-21, 1e-12, 99);
        assert!(n <= 2);
    }
    #[test]
    fn test_sd_diffusion_coeff_update() {
        let mut sd = make_sd();
        let init = vec![0.0, 1.0, 2.0];
        sd.positions = vec![1.0, 2.0, 3.0];
        let d = sd.update_diffusion_coefficient(&init, 1e-9);
        assert!(d >= 0.0);
    }
    fn make_cat() -> CatalysisMd {
        let mut cat = CatalysisMd::new(500.0, -1.5e-19, 1e-6, 100);
        cat.add_step(ReactionStep::new("ads", 0.5e-19, 1e13, -0.5e-19));
        cat.add_step(ReactionStep::new("rxn", 1.0e-19, 1e13, -1.0e-19));
        cat.add_step(ReactionStep::new("des", 0.3e-19, 1e13, 0.5e-19));
        cat
    }
    #[test]
    fn test_cat_compute_tof_positive() {
        let mut cat = make_cat();
        let tof = cat.compute_tof();
        assert!(tof > 0.0, "tof={tof}");
    }
    #[test]
    fn test_cat_sabatier_peak() {
        let cat = make_cat();
        let act = cat.sabatier_activity(cat.adsorption_energy, 1e-19);
        assert!((act - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cat_sabatier_off_peak() {
        let cat = make_cat();
        let act = cat.sabatier_activity(0.0, 1e-19);
        assert!(act < 1.0);
    }
    #[test]
    fn test_cat_reaction_step_rate() {
        let step = ReactionStep::new("test", 1e-20, 1e13, -1e-20);
        let rate = step.rate(300.0);
        assert!(rate > 0.0 && rate < 1e13);
    }
    #[test]
    fn test_cat_simulate_reactions() {
        let mut cat = make_cat();
        cat.simulate_reactions(100, 1e-6, 42);
        assert!(cat.total_turnovers > 0);
    }
    #[test]
    fn test_cat_measured_tof() {
        let mut cat = make_cat();
        cat.total_turnovers = 1000;
        let tof = cat.measured_tof(1.0);
        assert!((tof - 10.0).abs() < 1e-8, "tof={tof}");
    }
    #[test]
    fn test_cat_apparent_e_act() {
        let ea = CatalysisMd::apparent_activation_energy(1.0, 10.0, 300.0, 400.0);
        assert!(ea.is_finite());
    }
    #[test]
    fn test_cat_add_adsorbate() {
        let mut cat = make_cat();
        cat.add_adsorbate([0.0, 0.0, 2.0], [0.1, 0.0, -0.1], 18.0);
        assert_eq!(cat.adsorbate_positions.len(), 1);
    }
    #[test]
    fn test_cat_langevin_step() {
        let mut cat = make_cat();
        cat.add_adsorbate([0.0, 0.0, 2.0], [0.0, 0.0, 0.0], 18.0);
        cat.langevin_step(0.001, 1.0, 0.1, 12345);
        let z = cat.adsorbate_positions[0][2];
        assert!(z.is_finite());
    }
    fn make_recon() -> SurfaceReconstruction {
        let atoms: Vec<SurfaceAtom> = (0..10)
            .map(|i| SurfaceAtom::new([i as f64, 0.0, 2.0], -3.5))
            .collect();
        SurfaceReconstruction::new(atoms, 0.5, 10.0, 300.0)
    }
    #[test]
    fn test_recon_no_displacement_initially() {
        let recon = make_recon();
        for atom in &recon.atoms {
            assert_eq!(atom.displacement(), 0.0);
        }
    }
    #[test]
    fn test_recon_detect_none() {
        let mut recon = make_recon();
        let n = recon.detect_reconstruction();
        assert_eq!(n, 0);
    }
    #[test]
    fn test_recon_detect_after_displacement() {
        let mut recon = make_recon();
        recon.atoms[0].pos[0] += 1.0;
        let n = recon.detect_reconstruction();
        assert!(n >= 1);
    }
    #[test]
    fn test_recon_fraction_range() {
        let mut recon = make_recon();
        recon.atoms[0].pos[0] += 1.0;
        recon.detect_reconstruction();
        let frac = recon.reconstruction_fraction();
        assert!((0.0..=1.0).contains(&frac));
    }
    #[test]
    fn test_recon_mean_binding_energy() {
        let mut recon = make_recon();
        let mbe = recon.update_mean_binding_energy();
        assert!((mbe - (-3.5)).abs() < 1e-10);
    }
    #[test]
    fn test_recon_mc_relax_runs() {
        let mut recon = make_recon();
        let accepted = recon.mc_relax_step(0.1, 42);
        assert!(accepted <= recon.atoms.len());
    }
    #[test]
    fn test_recon_potential_energy_zero() {
        let recon = make_recon();
        assert_eq!(recon.potential_energy(), 0.0);
    }
    #[test]
    fn test_recon_relax() {
        let mut recon = make_recon();
        recon.relax(0.1);
        for atom in &recon.atoms {
            assert!(atom.pos[2].is_finite());
        }
    }
    fn make_interface() -> InterfaceMd {
        InterfaceMd::new(300.0, 0.072, 0.035, 0.072, [100.0, 100.0, 50.0], 25.0)
    }
    #[test]
    fn test_interface_contact_angle() {
        let imd = make_interface();
        let angle = imd.contact_angle_degrees();
        assert!((0.0..=180.0).contains(&angle), "angle={angle}");
    }
    #[test]
    fn test_interface_wetting_state_hydrophilic() {
        let imd = InterfaceMd::new(300.0, 0.1, 0.028, 0.072, [100.0, 100.0, 50.0], 25.0);
        let state = imd.wetting_state();
        assert_eq!(state, "superhydrophilic");
    }
    #[test]
    fn test_interface_wetting_state_hydrophobic() {
        let imd = InterfaceMd::new(300.0, 0.02, 0.065, 0.072, [100.0, 100.0, 50.0], 25.0);
        let state = imd.wetting_state();
        assert!(state == "hydrophobic" || state == "superhydrophobic");
    }
    #[test]
    fn test_interface_density_profile_empty() {
        let imd = make_interface();
        let profile = imd.density_profile(5.0);
        assert!(!profile.is_empty());
    }
    #[test]
    fn test_interface_density_profile_with_molecules() {
        let mut imd = make_interface();
        for i in 0..20 {
            imd.add_molecule(InterfaceMolecule::new(
                [5.0 * i as f64, 0.0, 10.0 + i as f64],
                [0.0; 3],
                18.0,
                0.0,
                true,
            ));
        }
        let profile = imd.density_profile(5.0);
        let total: f64 = profile.iter().map(|(_, d)| *d).sum();
        assert!(total > 0.0);
    }
    #[test]
    fn test_interface_tension_finite() {
        let mut imd = make_interface();
        let t = imd.compute_interfacial_tension(1e5, 1e5, 1.1e5);
        assert!(t.is_finite());
    }
    #[test]
    fn test_interface_integrate_step() {
        let mut imd = make_interface();
        imd.add_molecule(InterfaceMolecule::new(
            [0.0, 0.0, 30.0],
            [0.0, 0.0, 0.1],
            18.0,
            0.0,
            true,
        ));
        imd.integrate_step(0.001, 0.5);
        let z = imd.molecules[0].pos[2];
        assert!(z.is_finite());
    }
    #[test]
    fn test_interface_fraction() {
        let mut imd = make_interface();
        imd.add_molecule(InterfaceMolecule::new(
            [0.0, 0.0, 25.0],
            [0.0; 3],
            18.0,
            0.0,
            true,
        ));
        imd.add_molecule(InterfaceMolecule::new(
            [0.0, 0.0, 48.0],
            [0.0; 3],
            18.0,
            0.0,
            true,
        ));
        let frac = imd.interface_fraction(5.0);
        assert!((0.0..=1.0).contains(&frac));
    }
    fn make_tribofilm() -> TriboFilm {
        TriboFilm::new(1e-9, 1e-12, 1e-13, 1e8, 1e9)
    }
    fn make_tribo() -> Tribochemistry {
        Tribochemistry::new(
            10.0,
            0.1,
            0.3,
            300.0,
            2e-19,
            1e-29,
            1e13,
            1e9,
            make_tribofilm(),
        )
    }
    #[test]
    fn test_tribo_shear_stress_positive() {
        let tribo = make_tribo();
        let tau = tribo.shear_stress(1e-6);
        assert!(tau > 0.0, "tau={tau}");
    }
    #[test]
    fn test_tribo_shear_stress_zero_area() {
        let tribo = make_tribo();
        assert_eq!(tribo.shear_stress(0.0), 0.0);
    }
    #[test]
    fn test_tribo_bond_breaking_rate_positive() {
        let tribo = make_tribo();
        let rate = tribo.bond_breaking_rate(1e8);
        assert!(rate > 0.0, "rate={rate}");
    }
    #[test]
    fn test_tribo_bond_breaking_increases_with_stress() {
        let tribo = make_tribo();
        let r1 = tribo.bond_breaking_rate(1e7);
        let r2 = tribo.bond_breaking_rate(1e9);
        assert!(r2 > r1, "r1={r1} r2={r2}");
    }
    #[test]
    fn test_tribo_archard_wear_rate() {
        let tribo = make_tribo();
        let rate = tribo.archard_wear_rate(1e-4);
        assert!(rate > 0.0, "rate={rate}");
    }
    #[test]
    fn test_tribo_friction_force() {
        let tribo = make_tribo();
        let ff = tribo.friction_force();
        assert!((ff - 3.0).abs() < 1e-10, "ff={ff}");
    }
    #[test]
    fn test_tribo_step_advances_distance() {
        let mut tribo = make_tribo();
        let d_before = tribo.sliding_distance;
        tribo.step(0.01, 1e-6, 1e-4, 0, 42);
        assert!(tribo.sliding_distance > d_before);
    }
    #[test]
    fn test_tribo_step_tribofilm_evolves() {
        let mut tribo = make_tribo();
        let h_before = tribo.tribofilm.thickness;
        tribo.step(0.01, 1e-6, 1e-4, 0, 1234);
        assert!(tribo.tribofilm.thickness >= 0.0);
        let _ = h_before;
    }
    #[test]
    fn test_tribo_specific_wear_zero_initially() {
        let tribo = make_tribo();
        assert_eq!(tribo.specific_wear_rate(), 0.0);
    }
    #[test]
    fn test_tribofilm_steady_state() {
        let tf = make_tribofilm();
        let h_ss = tf.steady_state_thickness();
        assert!(h_ss > 0.0, "h_ss={h_ss}");
    }
    #[test]
    fn test_tribofilm_evolve() {
        let mut tf = make_tribofilm();
        let h0 = tf.thickness;
        tf.evolve(1.0);
        assert!(tf.thickness.is_finite());
        assert!(tf.thickness >= 0.0);
        let _ = h0;
    }
    #[test]
    fn test_tribofilm_evolve_no_negative() {
        let mut tf = TriboFilm::new(0.0, 1e-15, 1e-10, 1e6, 1e12);
        tf.evolve(1000.0);
        assert!(tf.thickness >= 0.0);
    }
    #[test]
    fn test_tribo_wear_particle_new() {
        let p = WearParticle::new(1e-18, 1e9, 5);
        assert_eq!(p.generation_step, 5);
        assert!(p.volume > 0.0);
    }
    #[test]
    fn test_tribo_n_wear_particles() {
        let tribo = make_tribo();
        assert_eq!(tribo.n_wear_particles(), 0);
    }
    #[test]
    fn test_tribo_mean_particle_volume_empty() {
        let tribo = make_tribo();
        assert_eq!(tribo.mean_particle_volume(), 0.0);
    }
    #[test]
    fn test_tribo_multiple_steps() {
        let mut tribo = make_tribo();
        for i in 0..10 {
            tribo.step(0.001, 1e-6, 1e-4, i, (i as u64) * 0xcafe + 1);
        }
        assert!(tribo.sliding_distance > 0.0);
        assert!(tribo.bonds_broken <= 20);
    }
    #[test]
    fn test_precursor_state_new() {
        let p = PrecursorState::new(1.5, -5e-21);
        assert!(!p.chemisorbed);
        assert_eq!(p.position, 1.5);
    }
    #[test]
    fn test_surface_atom_displacement() {
        let mut a = SurfaceAtom::new([0.0, 0.0, 0.0], -3.5);
        a.pos[0] = 1.0;
        a.pos[1] = 0.0;
        a.pos[2] = 0.0;
        assert!((a.displacement() - 1.0).abs() < 1e-10);
    }
}
