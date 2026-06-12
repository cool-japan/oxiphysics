//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::functions::*;
    use super::super::types::*;
    use std::f64::consts::PI;
    /// TIP3P charges must sum to zero: q_O + 2*q_H = 0.
    #[test]
    fn tip3p_charge_neutrality() {
        let p = WaterParams::tip3p();
        let total = p.total_charge();
        assert!(total.abs() < 1e-10, "TIP3P charge sum {total} != 0");
    }
    /// SPC charges must sum to zero.
    #[test]
    fn spc_charge_neutrality() {
        let p = WaterParams::spc();
        let total = p.total_charge();
        assert!(total.abs() < 1e-10, "SPC charge sum {total} != 0");
    }
    /// SPC/E charges must sum to zero.
    #[test]
    fn spce_charge_neutrality() {
        let p = WaterParams::spce();
        let total = p.total_charge();
        assert!(total.abs() < 1e-10, "SPC/E charge sum {total} != 0");
    }
    /// TIP4P charges must sum to zero: q_O + 2*q_H + q_M = 0.
    #[test]
    fn tip4p_charge_neutrality() {
        let p = WaterParams::tip4p();
        let total = p.total_charge();
        assert!(total.abs() < 1e-10, "TIP4P charge sum {total} != 0");
    }
    /// Hydrogen atoms must be placed at the correct bond length from oxygen.
    #[test]
    fn water_molecule_bond_length() {
        let mol = WaterMolecule::new([0.0, 0.0, 0.0]);
        let r1 = dist(&mol.oxygen, &mol.hydrogen1);
        let r2 = dist(&mol.oxygen, &mol.hydrogen2);
        let expected = WaterMolecule::bond_length();
        assert!(
            (r1 - expected).abs() < 1e-10,
            "O-H1 distance {r1} != {expected}"
        );
        assert!(
            (r2 - expected).abs() < 1e-10,
            "O-H2 distance {r2} != {expected}"
        );
    }
    /// SPC geometry has different bond length (1.0 Å) and angle (109.47°).
    #[test]
    fn spc_geometry_bond_length() {
        let mol = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let geom = WaterGeometry::spc();
        let r1 = dist(&mol.oxygen, &mol.hydrogen1);
        let r2 = dist(&mol.oxygen, &mol.hydrogen2);
        assert!(
            (r1 - geom.r_oh).abs() < 1e-10,
            "SPC O-H1 = {r1}, expected {}",
            geom.r_oh
        );
        assert!(
            (r2 - geom.r_oh).abs() < 1e-10,
            "SPC O-H2 = {r2}, expected {}",
            geom.r_oh
        );
    }
    /// SPC H-H distance should match the tetrahedral angle.
    #[test]
    fn spc_hh_distance() {
        let mol = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let geom = WaterGeometry::spc();
        let r_hh = dist(&mol.hydrogen1, &mol.hydrogen2);
        let expected = geom.r_hh();
        assert!(
            (r_hh - expected).abs() < 1e-10,
            "SPC H-H = {r_hh}, expected {expected}"
        );
    }
    /// TIP4P M-site should be present and at the correct distance from O.
    #[test]
    fn tip4p_m_site_position() {
        let mol = WaterMolecule::tip4p([0.0, 0.0, 0.0]);
        assert!(mol.m_site.is_some(), "TIP4P should have M-site");
        let m = mol.m_site.unwrap();
        let r_om = dist(&mol.oxygen, &m);
        let expected = WaterGeometry::tip4p().r_om;
        assert!(
            (r_om - expected).abs() < 1e-10,
            "O-M distance {r_om} != {expected}"
        );
    }
    /// TIP3P dipole moment must be non-zero.
    #[test]
    fn dipole_moment_nonzero() {
        let mol = WaterMolecule::new([0.0, 0.0, 0.0]);
        let params = WaterParams::tip3p();
        let p = mol.dipole_moment(&params);
        let mag = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!(mag > 1e-6, "Dipole magnitude {mag} is unexpectedly zero");
    }
    /// TIP4P dipole should also be non-zero.
    #[test]
    fn tip4p_dipole_moment_nonzero() {
        let mol = WaterMolecule::tip4p([0.0, 0.0, 0.0]);
        let params = WaterParams::tip4p();
        let d = mol.dipole_magnitude_debye(&params);
        assert!(d > 1.0, "TIP4P dipole {d} D should be > 1 D");
    }
    /// Two molecules separated by 10 Å should have small total pair energy.
    #[test]
    fn pair_energy_far_apart_is_small() {
        let coulomb_k = 1389.35;
        let mol_a = WaterMolecule::new([0.0, 0.0, 0.0]);
        let mol_b = WaterMolecule::new([0.0, 0.0, 10.0]);
        let params = WaterParams::tip3p();
        let e = water_total_pair_energy(&mol_a, &mol_b, &params, coulomb_k);
        assert!(e.is_finite(), "Pair energy is not finite: {e}");
        assert!(e.abs() < 10.0, "Pair energy at 10 Å too large: {e} kJ/mol");
    }
    /// LJ energy at r = sigma equals 0; minimum is at r = 2^(1/6)*sigma with energy -epsilon.
    #[test]
    fn lj_energy_at_minimum() {
        let params = WaterParams::tip3p();
        let r_min = 2.0_f64.powf(1.0 / 6.0) * params.sigma_o;
        let e = water_lj_energy(r_min, &params);
        assert!(
            (e - (-params.epsilon_o)).abs() < 1e-10,
            "LJ energy at minimum {e} != -epsilon = {}",
            -params.epsilon_o
        );
    }
    /// LJ energy at r = sigma should be exactly zero.
    #[test]
    fn lj_energy_at_sigma_is_zero() {
        let params = WaterParams::spc();
        let e = water_lj_energy(params.sigma_o, &params);
        assert!(e.abs() < 1e-10, "LJ energy at sigma should be 0, got {e}");
    }
    /// LJ force should be zero at r = 2^(1/6)*sigma (potential minimum).
    #[test]
    fn lj_force_at_minimum_is_zero() {
        let params = WaterParams::spce();
        let r_min = 2.0_f64.powf(1.0 / 6.0) * params.sigma_o;
        let f = water_lj_force(r_min, &params);
        assert!(f.abs() < 1e-10, "LJ force at minimum should be 0, got {f}");
    }
    /// Constraint distances should match geometry.
    #[test]
    fn constraint_distances_match_geometry() {
        for model in [
            WaterModelType::Tip3p,
            WaterModelType::Spc,
            WaterModelType::Spce,
            WaterModelType::Tip4p,
        ] {
            let params = WaterParams::from_model_type(model);
            let (r_oh, r_hh) = params.constraint_distances();
            assert!((r_oh - params.geometry.r_oh).abs() < 1e-12);
            assert!((r_hh - params.geometry.r_hh()).abs() < 1e-12);
        }
    }
    /// SETTLE should restore bond constraints after perturbation.
    #[test]
    fn settle_restores_constraints() {
        let mut mol = WaterMolecule::new([0.0, 0.0, 0.0]);
        mol.hydrogen1[0] += 0.05;
        mol.hydrogen2[1] -= 0.03;
        rigid_water_settle(&mut mol, 0.001)
            .expect("settle should succeed for ideal-reference compat wrapper");
        let geom = WaterGeometry::tip3p();
        assert!(
            mol.check_constraints(&geom, 0.01),
            "SETTLE should restore constraints"
        );
    }
    /// SETTLE with SPC geometry.
    #[test]
    fn settle_spc_geometry() {
        let geom = WaterGeometry::spc();
        let mut mol = WaterMolecule::spc([1.0, 2.0, 3.0]);
        mol.hydrogen1[0] += 0.04;
        mol.hydrogen2[1] -= 0.02;
        rigid_water_project_geom(&mut mol, &geom);
        assert!(mol.check_constraints(&geom, 0.01));
    }
    /// Water box creation should produce correct number of molecules.
    #[test]
    fn water_box_creation() {
        let geom = WaterGeometry::tip3p();
        let mols = create_water_box(8, 3.0, &geom);
        assert_eq!(mols.len(), 8, "Should create 8 molecules");
        for mol in &mols {
            assert!(mol.check_constraints(&geom, 1e-8));
        }
    }
    /// TIP4P pair energy should be finite for separated molecules.
    #[test]
    fn tip4p_pair_energy_finite() {
        let mol_a = WaterMolecule::tip4p([0.0, 0.0, 0.0]);
        let mol_b = WaterMolecule::tip4p([0.0, 0.0, 5.0]);
        let params = WaterParams::tip4p();
        let e = water_tip4p_pair_energy(&mol_a, &mol_b, &params, 1389.35);
        assert!(e.is_finite(), "TIP4P pair energy should be finite, got {e}");
    }
    /// SPC/E self-polarisation energy should be positive.
    #[test]
    fn spce_self_polarisation_positive() {
        let e = spce_self_polarisation_energy();
        assert!(e > 0.0, "SPC/E self-polarisation should be positive");
    }
    /// Center of mass should be between O and the H atoms.
    #[test]
    fn center_of_mass_position() {
        let mol = WaterMolecule::new([0.0, 0.0, 0.0]);
        let params = WaterParams::tip3p();
        let com = mol.center_of_mass(&params);
        assert!(com[1] > 0.0, "COM y should be positive (H atoms are in +y)");
        assert!(com[0].abs() < 1e-10, "COM x should be ~0 by symmetry");
    }
    /// from_model_type should return consistent parameters.
    #[test]
    fn from_model_type_consistency() {
        let tip3p = WaterParams::from_model_type(WaterModelType::Tip3p);
        assert_eq!(tip3p.name, "TIP3P");
        let spc = WaterParams::from_model_type(WaterModelType::Spc);
        assert_eq!(spc.name, "SPC");
        let spce = WaterParams::from_model_type(WaterModelType::Spce);
        assert_eq!(spce.name, "SPC/E");
        let tip4p = WaterParams::from_model_type(WaterModelType::Tip4p);
        assert_eq!(tip4p.name, "TIP4P");
    }
    /// M-site update should place the M-site along the bisector.
    #[test]
    fn m_site_update_on_bisector() {
        let mut mol = WaterMolecule::tip4p([0.0, 0.0, 0.0]);
        let geom = WaterGeometry::tip4p();
        mol.update_m_site(&geom);
        let m = mol.m_site.unwrap();
        assert!(m[1] > 0.0, "M-site y should be positive");
        assert!(m[0].abs() < 1e-10, "M-site x should be ~0 by symmetry");
    }
    #[test]
    fn tip5p_charge_neutrality() {
        let p = Tip5pParams::new();
        let total = p.total_charge();
        assert!(total.abs() < 1e-10, "TIP5P charge sum {total} != 0");
    }
    #[test]
    fn tip5p_molecule_h_bond_lengths() {
        let mol = Tip5pMolecule::new([0.0, 0.0, 0.0]);
        assert!(
            mol.check_bond_lengths(1e-10),
            "TIP5P O-H bond lengths incorrect"
        );
    }
    #[test]
    fn tip5p_molecule_lp_distance() {
        let mol = Tip5pMolecule::new([0.0, 0.0, 0.0]);
        let p = Tip5pParams::new();
        let r_lp = mol.r_olp1();
        assert!(
            (r_lp - p.r_olp).abs() < 1e-10,
            "TIP5P O-LP distance = {r_lp}, expected {}",
            p.r_olp
        );
    }
    #[test]
    fn tip5p_has_five_sites() {
        let mol = Tip5pMolecule::new([1.0, 2.0, 3.0]);
        let all_pos = [mol.oxygen, mol.hydrogen1, mol.hydrogen2, mol.lp1, mol.lp2];
        for p in &all_pos {
            assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
        }
    }
    #[test]
    fn tip5p_symmetry_lp_sites() {
        let mol = Tip5pMolecule::new([0.0, 0.0, 0.0]);
        assert!(
            (mol.lp1[0] + mol.lp2[0]).abs() < 1e-10,
            "LP sites should be symmetric"
        );
        assert!(
            (mol.lp1[1] - mol.lp2[1]).abs() < 1e-10,
            "LP sites should have same y"
        );
    }
    #[test]
    fn fspc_bond_energy_zero_at_equilibrium() {
        let p = FlexibleSpcParams::new();
        let e = p.bond_energy(p.r_oh);
        assert!(
            e.abs() < 1e-12,
            "Bond energy at equilibrium should be 0, got {e}"
        );
    }
    #[test]
    fn fspc_bond_energy_positive_for_stretch() {
        let p = FlexibleSpcParams::new();
        let e = p.bond_energy(p.r_oh + 0.1);
        assert!(e > 0.0, "Bond energy for stretched bond should be positive");
    }
    #[test]
    fn fspc_angle_energy_zero_at_equilibrium() {
        let p = FlexibleSpcParams::new();
        let theta0 = p.angle_hoh_deg * PI / 180.0;
        let e = p.angle_energy(theta0);
        assert!(
            e.abs() < 1e-12,
            "Angle energy at equilibrium should be 0, got {e}"
        );
    }
    #[test]
    fn fspc_intramolecular_energy_at_equilibrium_near_zero() {
        let p = FlexibleSpcParams::new();
        let mol = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let e = p.intramolecular_energy(&mol);
        assert!(
            e < 1.0,
            "Intramolecular energy at equilibrium should be small, got {e}"
        );
    }
    #[test]
    fn fspc_intramolecular_energy_increases_on_perturbation() {
        let p = FlexibleSpcParams::new();
        let mol_eq = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let e_eq = p.intramolecular_energy(&mol_eq);
        let mut mol_perturbed = mol_eq.clone();
        mol_perturbed.hydrogen1[0] += 0.3;
        let e_pert = p.intramolecular_energy(&mol_perturbed);
        assert!(
            e_pert > e_eq,
            "Perturbation should increase intramolecular energy"
        );
    }
    #[test]
    fn compare_water_models_returns_four() {
        let summaries = compare_water_models();
        assert_eq!(summaries.len(), 4, "Should return 4 model summaries");
    }
    #[test]
    fn compare_water_models_dipoles_nonzero() {
        let summaries = compare_water_models();
        for s in &summaries {
            assert!(
                s.dipole_debye > 0.5,
                "Dipole of {} should be > 0.5 D, got {}",
                s.name,
                s.dipole_debye
            );
        }
    }
    #[test]
    fn compare_water_models_tip4p_has_virtual_site() {
        let summaries = compare_water_models();
        let tip4p = summaries.iter().find(|s| s.name == "TIP4P").unwrap();
        assert!(tip4p.has_virtual_site, "TIP4P should have a virtual site");
    }
    #[test]
    fn rdf_length_matches_n_bins() {
        let geom = WaterGeometry::tip3p();
        let oxygens: Vec<[f64; 3]> = create_water_box(8, 3.0, &geom)
            .iter()
            .map(|m| m.oxygen)
            .collect();
        let (r, g) = oxygen_rdf(&oxygens, 9.0, 4.0, 20);
        assert_eq!(r.len(), 20);
        assert_eq!(g.len(), 20);
    }
    #[test]
    fn rdf_values_non_negative() {
        let geom = WaterGeometry::spc();
        let oxygens: Vec<[f64; 3]> = create_water_box(27, 3.1, &geom)
            .iter()
            .map(|m| m.oxygen)
            .collect();
        let box_l = 3.1 * 3.0;
        let (_, g) = oxygen_rdf(&oxygens, box_l, 5.0, 20);
        for &v in &g {
            assert!(v >= 0.0, "RDF value should be non-negative, got {v}");
        }
    }
    #[test]
    fn hydrogen_bond_count_positive() {
        let geom = WaterGeometry::tip3p();
        let mols = create_water_box(8, 3.0, &geom);
        let box_l = 3.0 * 2.5;
        let n_hb = count_hydrogen_bonds_simple(&mols, box_l, 3.5);
        assert!(n_hb >= 0.0);
    }
    #[test]
    fn tetrahedral_order_parameter_range() {
        let geom = WaterGeometry::tip3p();
        let oxygens: Vec<[f64; 3]> = create_water_box(27, 3.0, &geom)
            .iter()
            .map(|m| m.oxygen)
            .collect();
        let box_l = 3.0 * 3.0;
        let q_vals = tetrahedral_order_parameter(&oxygens, box_l);
        assert_eq!(q_vals.len(), oxygens.len());
        let mean_q: f64 = q_vals.iter().sum::<f64>() / q_vals.len() as f64;
        assert!(mean_q.is_finite(), "Mean q should be finite");
    }
    #[test]
    fn tetrahedral_order_few_molecules() {
        let geom = WaterGeometry::tip3p();
        let oxygens: Vec<[f64; 3]> = create_water_box(4, 3.0, &geom)
            .iter()
            .map(|m| m.oxygen)
            .collect();
        let q_vals = tetrahedral_order_parameter(&oxygens, 9.0);
        assert_eq!(q_vals.len(), 4);
    }
    #[test]
    fn dimer_binding_energy_tip3p_negative_at_equilibrium() {
        let coulomb_k = 1389.35_f64;
        let params = WaterParams::tip3p();
        let mol_a = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &WaterGeometry::tip3p());
        let mol_b = WaterMolecule::with_geometry([0.0, 0.0, 2.8], &WaterGeometry::tip3p());
        let e = water_dimer_binding_energy(&mol_a, &mol_b, &params, coulomb_k);
        assert!(e.is_finite(), "Dimer energy should be finite, got {e}");
    }
    #[test]
    fn dimer_binding_energy_spce_separation() {
        let e_close = spce_dimer_energy_at_separation(2.5);
        let e_far = spce_dimer_energy_at_separation(8.0);
        assert!(
            e_close.abs() >= e_far.abs(),
            "Close dimer energy |{}| should be >= far |{}|",
            e_close.abs(),
            e_far.abs()
        );
    }
    #[test]
    fn dimer_binding_energy_tip4p_finite() {
        let coulomb_k = 1389.35_f64;
        let params = WaterParams::tip4p();
        let mol_a = WaterMolecule::tip4p([0.0, 0.0, 0.0]);
        let mol_b = WaterMolecule::tip4p([0.0, 0.0, 3.0]);
        let e = water_dimer_binding_energy(&mol_a, &mol_b, &params, coulomb_k);
        assert!(
            e.is_finite(),
            "TIP4P dimer energy should be finite, got {e}"
        );
    }
    #[test]
    fn clausius_mossotti_dielectric_greater_than_one() {
        let mu_d = 2.35_f64;
        let n_mol = 512.0_f64;
        let box_l_a = 24.8_f64;
        let vol = box_l_a * box_l_a * box_l_a;
        let eps = clausius_mossotti_dielectric(mu_d, n_mol, vol, 300.0);
        assert!(eps > 1.0, "Dielectric constant should be > 1, got {eps}");
    }
    #[test]
    fn clausius_mossotti_dielectric_returns_one_on_bad_params() {
        let eps = clausius_mossotti_dielectric(2.35, 512.0, 0.0, 300.0);
        assert_eq!(eps, 1.0, "Zero volume should return 1.0");
    }
    #[test]
    fn kirkwood_g_factor_parallel_dipoles_greater_than_one() {
        let dipoles: Vec<[f64; 3]> = (0..10).map(|_| [0.0_f64, 0.0, 1.85]).collect();
        let g = kirkwood_g_factor(&dipoles);
        assert!(g > 1.0, "Parallel dipoles should give g_K > 1, got {g}");
    }
    #[test]
    fn kirkwood_g_factor_empty() {
        let g = kirkwood_g_factor(&[]);
        assert_eq!(g, 0.0);
    }
    #[test]
    fn oh_rdf_length_matches_n_bins() {
        let geom = WaterGeometry::tip3p();
        let mols = create_water_box(8, 3.0, &geom);
        let oxygens: Vec<[f64; 3]> = mols.iter().map(|m| m.oxygen).collect();
        let hydrogens: Vec<[f64; 3]> = mols
            .iter()
            .flat_map(|m| [m.hydrogen1, m.hydrogen2])
            .collect();
        let (r, g) = oh_rdf(&oxygens, &hydrogens, 9.0, 4.5, 20);
        assert_eq!(r.len(), 20, "r_centers should have 20 elements");
        assert_eq!(g.len(), 20, "g_OH should have 20 elements");
    }
    #[test]
    fn oh_rdf_values_non_negative() {
        let geom = WaterGeometry::spce();
        let mols = create_water_box(27, 3.1, &geom);
        let oxygens: Vec<[f64; 3]> = mols.iter().map(|m| m.oxygen).collect();
        let hydrogens: Vec<[f64; 3]> = mols
            .iter()
            .flat_map(|m| [m.hydrogen1, m.hydrogen2])
            .collect();
        let box_l = 3.1 * 3.0;
        let (_, g) = oh_rdf(&oxygens, &hydrogens, box_l, 5.0, 25);
        for &v in &g {
            assert!(v >= 0.0, "g_OH must be non-negative, got {v}");
        }
    }
    #[test]
    fn hh_rdf_length_correct() {
        let geom = WaterGeometry::tip3p();
        let mols = create_water_box(8, 3.0, &geom);
        let hydrogens: Vec<[f64; 3]> = mols
            .iter()
            .flat_map(|m| [m.hydrogen1, m.hydrogen2])
            .collect();
        let (r, g) = hh_rdf(&hydrogens, 9.0, 4.0, 15);
        assert_eq!(r.len(), 15);
        assert_eq!(g.len(), 15);
    }
    #[test]
    fn hh_rdf_values_non_negative() {
        let geom = WaterGeometry::spc();
        let mols = create_water_box(8, 3.0, &geom);
        let hydrogens: Vec<[f64; 3]> = mols
            .iter()
            .flat_map(|m| [m.hydrogen1, m.hydrogen2])
            .collect();
        let (_, g) = hh_rdf(&hydrogens, 9.0, 4.0, 15);
        for &v in &g {
            assert!(v >= 0.0, "g_HH must be >= 0, got {v}");
        }
    }
    #[test]
    fn ice_ih_density_close_to_physical() {
        let ice = IceIhCell::new();
        let rho = ice.density_gcc();
        assert!(
            rho > 0.85 && rho < 1.0,
            "Ice Ih density should be ~0.917 g/cm³, got {rho}"
        );
    }
    #[test]
    fn ice_ih_cell_volume_positive() {
        let ice = IceIhCell::new();
        let vol = ice.cell_volume();
        assert!(vol > 0.0, "Unit cell volume should be positive, got {vol}");
    }
    #[test]
    fn ice_ih_number_density_positive() {
        let ice = IceIhCell::new();
        let rho_n = ice.number_density();
        assert!(
            rho_n > 0.0,
            "Number density should be positive, got {rho_n}"
        );
    }
    #[test]
    fn water_model_estimated_density_spce_reasonable() {
        let params = WaterParams::spce();
        let rho = water_model_estimated_density(&params, 0.64);
        assert!(
            rho > 0.3 && rho < 3.0,
            "Estimated density should be plausible, got {rho} g/cm³"
        );
    }
    #[test]
    fn spce_polarisation_correction_positive() {
        let e = spce_polarisation_correction(2.35);
        assert!(
            e > 0.0,
            "Polarisation correction should be positive, got {e}"
        );
    }
    #[test]
    fn spce_polarisation_correction_zero_at_gas_phase_dipole() {
        let e = spce_polarisation_correction(1.85);
        assert!(
            e.abs() < 1e-12,
            "Correction should be 0 at gas-phase dipole, got {e}"
        );
    }
    #[test]
    fn spce_polarisation_correction_known_value() {
        let e = spce_polarisation_correction(2.35);
        assert!(
            (e - 5.22).abs() < 2.0,
            "SPC/E correction should be ~5.22 kJ/mol, got {e}"
        );
    }
    #[test]
    fn tip4p_force_redistribution_conserves_total() {
        let tip4p = Tip4p::default_params();
        let f_m = [1.0_f64, -2.0, 3.0];
        let err = tip4p.force_conservation_check(f_m);
        assert!(err < 1e-12, "Force conservation violated: max_err={err}");
    }
    #[test]
    fn tip4p_force_redistribution_zero_force() {
        let tip4p = Tip4p::default_params();
        let f_m = [0.0_f64, 0.0, 0.0];
        let (df_o, df_h1, df_h2) = tip4p.compute_virtual_site_force_redistribution(f_m);
        assert_eq!(df_o, [0.0; 3]);
        assert_eq!(df_h1, [0.0; 3]);
        assert_eq!(df_h2, [0.0; 3]);
    }
    #[test]
    fn tip4p_m_site_position_on_bisector() {
        let r_o = [0.0, 0.0, 0.0];
        let r_h1 = [0.756, 0.586, 0.0];
        let r_h2 = [-0.756, 0.586, 0.0];
        let tip4p = Tip4p::default_params();
        let m = tip4p.m_site_position(r_o, r_h1, r_h2);
        assert!(
            m[0].abs() < 1e-12,
            "M-site x should be 0 for symmetric geometry, got {}",
            m[0]
        );
        assert!(m[1] > 0.0, "M-site y should be positive, got {}", m[1]);
    }
    #[test]
    fn tip4p_h1_h2_get_equal_force_share() {
        let tip4p = Tip4p::default_params();
        let f_m = [5.0_f64, 3.0, -1.0];
        let (_, df_h1, df_h2) = tip4p.compute_virtual_site_force_redistribution(f_m);
        for k in 0..3 {
            assert!(
                (df_h1[k] - df_h2[k]).abs() < 1e-14,
                "H1 and H2 should get equal force shares"
            );
        }
    }
    #[test]
    fn spc_angle_energy_zero_at_equilibrium() {
        let spc = Spc::new();
        let r_o = [0.0, 0.0, 0.0];
        let theta_eq = spc.theta_0;
        let r_oh = 1.0_f64;
        let half = theta_eq / 2.0;
        let r_h1 = [r_oh * half.sin(), r_oh * half.cos(), 0.0];
        let r_h2 = [-r_oh * half.sin(), r_oh * half.cos(), 0.0];
        let e = spc.compute_angle_energy(r_o, r_h1, r_h2);
        assert!(
            e.abs() < 1e-8,
            "SPC angle energy at equilibrium should be ~0, got {e}"
        );
    }
    #[test]
    fn spc_angle_energy_positive_when_distorted() {
        let spc = Spc::new();
        let r_o = [0.0, 0.0, 0.0];
        let r_h1 = [1.0, 0.0, 0.0];
        let r_h2 = [0.0, 1.0, 0.0];
        let e = spc.compute_angle_energy(r_o, r_h1, r_h2);
        assert!(
            e > 0.0,
            "Distorted SPC angle should have positive energy, got {e}"
        );
    }
    #[test]
    fn spc_compute_angle_at_90_degrees() {
        let r_o = [0.0, 0.0, 0.0];
        let r_h1 = [1.0, 0.0, 0.0];
        let r_h2 = [0.0, 1.0, 0.0];
        let theta = Spc::compute_angle(r_o, r_h1, r_h2);
        let expected = std::f64::consts::FRAC_PI_2;
        assert!(
            (theta - expected).abs() < 1e-10,
            "90° angle: got {theta}, expected {expected}"
        );
    }
    #[test]
    fn water_cluster_single_molecule_zero_hbonds() {
        let mol = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let cluster = WaterCluster::new(vec![mol]);
        let n = cluster.compute_hydrogen_bond_count();
        assert_eq!(n, 0.0, "Single water molecule has no H-bonds");
    }
    #[test]
    fn water_cluster_dimer_has_hbond() {
        let mol1 = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let mol2 = WaterMolecule::spc([0.8165 * 2.8, 0.5774 * 2.8, 0.0]);
        let cluster = WaterCluster::new(vec![mol1, mol2]);
        let n = cluster.compute_hydrogen_bond_count();
        assert!(
            n >= 1.0,
            "Water dimer at 2.8 Å should have at least 1 H-bond, got {n}"
        );
    }
    #[test]
    fn water_cluster_distant_molecules_no_hbond() {
        let mol1 = WaterMolecule::spc([0.0, 0.0, 0.0]);
        let mol2 = WaterMolecule::spc([0.0, 10.0, 0.0]);
        let cluster = WaterCluster::new(vec![mol1, mol2]);
        let n = cluster.compute_hydrogen_bond_count();
        assert_eq!(n, 0.0, "Distant molecules should have 0 H-bonds, got {n}");
    }

    /// Deterministic LCG for reproducible perturbation tests (no `rand` dependency).
    struct Lcg {
        state: u64,
    }
    impl Lcg {
        fn new(seed: u64) -> Self {
            Self { state: seed }
        }
        fn next_f64(&mut self) -> f64 {
            self.state = self
                .state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.state >> 11) as f64) / ((1u64 << 53) as f64)
        }
        fn signed(&mut self) -> f64 {
            2.0 * self.next_f64() - 1.0
        }
    }

    /// Mass-weighted COM of a molecule's three real atoms.
    fn com_of(mol: &WaterMolecule, m_o: f64, m_h: f64) -> [f64; 3] {
        let m_t = m_o + 2.0 * m_h;
        [
            (m_o * mol.oxygen[0] + m_h * mol.hydrogen1[0] + m_h * mol.hydrogen2[0]) / m_t,
            (m_o * mol.oxygen[1] + m_h * mol.hydrogen1[1] + m_h * mol.hydrogen2[1]) / m_t,
            (m_o * mol.oxygen[2] + m_h * mol.hydrogen1[2] + m_h * mol.hydrogen2[2]) / m_t,
        ]
    }

    /// Build an ideal reference at `o` and a ≤5% randomly perturbed copy.
    fn make_reference_and_perturbed(
        o: [f64; 3],
        geom: &WaterGeometry,
        rng: &mut Lcg,
    ) -> (WaterMolecule, WaterMolecule) {
        let reference = WaterMolecule::with_geometry(o, geom);
        let mut unconstrained = reference.clone();
        let amp = 0.05 * geom.r_oh;
        for atom in [
            &mut unconstrained.oxygen,
            &mut unconstrained.hydrogen1,
            &mut unconstrained.hydrogen2,
        ] {
            for c in atom.iter_mut() {
                *c += rng.signed() * amp;
            }
        }
        (reference, unconstrained)
    }

    /// Constraint residuals are < 1e-12 for TIP3P and SPC/E after analytic SETTLE.
    #[test]
    fn settle_positions_constraint_rms_tip3p_and_spce() {
        let cases = [
            (WaterGeometry::tip3p(), WaterParams::tip3p()),
            (WaterGeometry::spce(), WaterParams::spce()),
        ];
        let n = 10_000;
        for (geom, params) in cases.iter() {
            let mut rng = Lcg::new(0x1234_5678_9abc_def0);
            let mut max_res = 0.0_f64;
            for i in 0..n {
                let o = [
                    (i as f64) * 0.001,
                    -(i as f64) * 0.0007,
                    (i as f64) * 0.0003,
                ];
                let (reference, mut unconstrained) =
                    make_reference_and_perturbed(o, geom, &mut rng);
                settle_positions(
                    &reference,
                    &mut unconstrained,
                    geom,
                    params.mass_o,
                    params.mass_h,
                )
                .expect("settle_positions should succeed for small perturbations");
                let r1 = dist(&unconstrained.oxygen, &unconstrained.hydrogen1);
                let r2 = dist(&unconstrained.oxygen, &unconstrained.hydrogen2);
                let rhh = dist(&unconstrained.hydrogen1, &unconstrained.hydrogen2);
                max_res = max_res.max((r1 - geom.r_oh).abs());
                max_res = max_res.max((r2 - geom.r_oh).abs());
                max_res = max_res.max((rhh - geom.r_hh()).abs());
            }
            assert!(
                max_res < 1e-12,
                "max constraint residual {max_res} exceeds 1e-12 for {}",
                params.name
            );
        }
    }

    /// SETTLE preserves the unconstrained centre of mass exactly.
    #[test]
    fn settle_positions_preserves_com() {
        let geom = WaterGeometry::tip3p();
        let params = WaterParams::tip3p();
        let mut rng = Lcg::new(0xdead_beef_cafe_babe);
        let mut max_dev = 0.0_f64;
        for i in 0..2_000 {
            let o = [(i as f64) * 0.002, 1.0, -(i as f64) * 0.001];
            let (reference, mut unconstrained) = make_reference_and_perturbed(o, &geom, &mut rng);
            let com_before = com_of(&unconstrained, params.mass_o, params.mass_h);
            settle_positions(
                &reference,
                &mut unconstrained,
                &geom,
                params.mass_o,
                params.mass_h,
            )
            .expect("settle_positions should succeed");
            let com_after = com_of(&unconstrained, params.mass_o, params.mass_h);
            for (before, after) in com_before.iter().zip(com_after.iter()) {
                max_dev = max_dev.max((before - after).abs());
            }
        }
        assert!(max_dev < 1e-14, "COM drift {max_dev} exceeds 1e-14");
    }

    /// Analytic SETTLE is the least-displacement rigid solution: its mass-weighted
    /// squared displacement from the unconstrained positions is <= that of the
    /// iterative SHAKE-style projector (which is not COM-preserving and not
    /// optimal). This validates that the analytic placement is the optimal rigid
    /// fit, not merely *a* rigid fit.
    #[test]
    fn settle_positions_agrees_with_shake_oracle() {
        let geom = WaterGeometry::tip3p();
        let params = WaterParams::tip3p();
        let m_o = params.mass_o;
        let m_h = params.mass_h;
        let mut rng = Lcg::new(0x0f0f_0f0f_1234_9999);
        for i in 0..1_000 {
            let o = [(i as f64) * 0.003, 0.5, (i as f64) * 0.0011];
            let (reference, unconstrained) = make_reference_and_perturbed(o, &geom, &mut rng);

            // Analytic SETTLE.
            let mut analytic = unconstrained.clone();
            settle_positions(&reference, &mut analytic, &geom, m_o, m_h)
                .expect("settle_positions should succeed");

            // Iterative SHAKE-style projection oracle.
            let mut shake = unconstrained.clone();
            rigid_water_project_geom(&mut shake, &geom);

            let msd = |c: &WaterMolecule| -> f64 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += m_o * (c.oxygen[k] - unconstrained.oxygen[k]).powi(2);
                    s += m_h * (c.hydrogen1[k] - unconstrained.hydrogen1[k]).powi(2);
                    s += m_h * (c.hydrogen2[k] - unconstrained.hydrogen2[k]).powi(2);
                }
                s
            };
            let analytic_msd = msd(&analytic);
            let shake_msd = msd(&shake);
            assert!(
                analytic_msd <= shake_msd + 1e-9,
                "analytic MSD {analytic_msd} should be <= SHAKE MSD {shake_msd} (+1e-9)"
            );
        }
    }

    /// After the velocity step every bond's projected relative velocity is ~0.
    #[test]
    fn settle_velocities_zero_bond_velocity() {
        let geom = WaterGeometry::tip3p();
        let params = WaterParams::tip3p();
        let mut rng = Lcg::new(0xabcd_1234_5678_9999);
        for _ in 0..1_000 {
            let mol = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
            let mut vel = WaterVelocities::new(
                [rng.signed(), rng.signed(), rng.signed()],
                [rng.signed(), rng.signed(), rng.signed()],
                [rng.signed(), rng.signed(), rng.signed()],
            );
            settle_velocities(&mol, &mut vel, params.mass_o, params.mass_h)
                .expect("settle_velocities should succeed");
            let bonds = [
                (mol.oxygen, mol.hydrogen1, vel.v_oxygen, vel.v_hydrogen1),
                (mol.oxygen, mol.hydrogen2, vel.v_oxygen, vel.v_hydrogen2),
                (
                    mol.hydrogen1,
                    mol.hydrogen2,
                    vel.v_hydrogen1,
                    vel.v_hydrogen2,
                ),
            ];
            for (ri, rj, vi, vj) in bonds {
                let dr = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
                let dv = [vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2]];
                let proj = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
                assert!(proj.abs() < 1e-12, "bond velocity {proj} not removed");
            }
        }
    }

    /// The velocity step conserves the centre-of-mass velocity exactly.
    #[test]
    fn settle_velocities_preserves_com_velocity() {
        let geom = WaterGeometry::tip3p();
        let params = WaterParams::tip3p();
        let m_o = params.mass_o;
        let m_h = params.mass_h;
        let m_t = m_o + 2.0 * m_h;
        let mut rng = Lcg::new(0x5555_aaaa_3333_cccc);
        let mut max_dev = 0.0_f64;
        for _ in 0..1_000 {
            let mol = WaterMolecule::with_geometry([1.0, -2.0, 0.5], &geom);
            let mut vel = WaterVelocities::new(
                [rng.signed(), rng.signed(), rng.signed()],
                [rng.signed(), rng.signed(), rng.signed()],
                [rng.signed(), rng.signed(), rng.signed()],
            );
            let com_v_before = [
                (m_o * vel.v_oxygen[0] + m_h * vel.v_hydrogen1[0] + m_h * vel.v_hydrogen2[0]) / m_t,
                (m_o * vel.v_oxygen[1] + m_h * vel.v_hydrogen1[1] + m_h * vel.v_hydrogen2[1]) / m_t,
                (m_o * vel.v_oxygen[2] + m_h * vel.v_hydrogen1[2] + m_h * vel.v_hydrogen2[2]) / m_t,
            ];
            settle_velocities(&mol, &mut vel, m_o, m_h).expect("settle_velocities should succeed");
            let com_v_after = [
                (m_o * vel.v_oxygen[0] + m_h * vel.v_hydrogen1[0] + m_h * vel.v_hydrogen2[0]) / m_t,
                (m_o * vel.v_oxygen[1] + m_h * vel.v_hydrogen1[1] + m_h * vel.v_hydrogen2[1]) / m_t,
                (m_o * vel.v_oxygen[2] + m_h * vel.v_hydrogen1[2] + m_h * vel.v_hydrogen2[2]) / m_t,
            ];
            for (before, after) in com_v_before.iter().zip(com_v_after.iter()) {
                max_dev = max_dev.max((before - after).abs());
            }
        }
        assert!(
            max_dev < 1e-14,
            "COM velocity drift {max_dev} exceeds 1e-14"
        );
    }

    /// Degenerate / non-finite reference geometries must return an error.
    #[test]
    fn settle_degenerate_reference_returns_err() {
        let geom = WaterGeometry::tip3p();
        let params = WaterParams::tip3p();
        // Collinear reference.
        let collinear = WaterMolecule {
            oxygen: [0.0, 0.0, 0.0],
            hydrogen1: [1.0, 0.0, 0.0],
            hydrogen2: [2.0, 0.0, 0.0],
            m_site: None,
        };
        let mut unconstrained = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
        assert!(
            settle_positions(
                &collinear,
                &mut unconstrained,
                &geom,
                params.mass_o,
                params.mass_h
            )
            .is_err()
        );
        // Zero-area (all identical).
        let coincident = WaterMolecule {
            oxygen: [1.0, 1.0, 1.0],
            hydrogen1: [1.0, 1.0, 1.0],
            hydrogen2: [1.0, 1.0, 1.0],
            m_site: None,
        };
        let mut u2 = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
        assert!(
            settle_positions(&coincident, &mut u2, &geom, params.mass_o, params.mass_h).is_err()
        );
        // Non-finite input.
        let mut nan_u = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
        nan_u.hydrogen1[0] = f64::NAN;
        let reference = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
        assert!(
            settle_positions(&reference, &mut nan_u, &geom, params.mass_o, params.mass_h).is_err()
        );
        // Non-positive mass.
        let mut u3 = WaterMolecule::with_geometry([0.0, 0.0, 0.0], &geom);
        assert!(settle_positions(&reference, &mut u3, &geom, 0.0, params.mass_h).is_err());
    }

    /// Random valid perturbations never produce NaN/inf output coordinates.
    #[test]
    fn settle_positions_no_nan_on_random() {
        let geom = WaterGeometry::spce();
        let params = WaterParams::spce();
        let mut rng = Lcg::new(0x7777_8888_9999_aaaa);
        for i in 0..5_000 {
            let o = [(i as f64) * 0.004, -(i as f64) * 0.002, 1.0];
            let (reference, mut unconstrained) = make_reference_and_perturbed(o, &geom, &mut rng);
            settle_positions(
                &reference,
                &mut unconstrained,
                &geom,
                params.mass_o,
                params.mass_h,
            )
            .expect("settle_positions should succeed");
            for atom in [
                unconstrained.oxygen,
                unconstrained.hydrogen1,
                unconstrained.hydrogen2,
            ] {
                for c in atom {
                    assert!(c.is_finite(), "non-finite output coordinate {c}");
                }
            }
        }
    }
}
