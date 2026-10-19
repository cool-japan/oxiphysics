//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::nano::BioMaterial;

    use crate::nano::CompositeLaminate;

    use crate::nano::FiberMatrix;

    use crate::nano::GrapheneSheet;
    use crate::nano::InterfacialZone;
    use crate::nano::MolecularCrystal;

    use crate::nano::NanoparticleDispersion;
    use crate::nano::NanotubeProperties;

    use crate::nano::Ply;
    use crate::nano::PolymerChain;

    use crate::nano::SizeEffect;

    use crate::nano::ThermoelectricMaterial;

    #[test]
    fn test_cnt_armchair_diameter() {
        let cnt = NanotubeProperties::armchair(10);
        assert!(cnt.diameter_nm > 0.0);
        assert!(cnt.diameter_nm < 5.0);
    }
    #[test]
    fn test_cnt_zigzag_diameter() {
        let cnt = NanotubeProperties::zigzag(10);
        assert!(cnt.diameter_nm > 0.0);
        assert_eq!(cnt.m, 0);
    }
    #[test]
    fn test_cnt_chiral() {
        let cnt = NanotubeProperties::chiral(8, 4);
        assert_eq!(cnt.n, 8);
        assert_eq!(cnt.m, 4);
    }
    #[test]
    fn test_cnt_is_metallic_armchair() {
        let cnt = NanotubeProperties::armchair(5);
        assert!(cnt.is_metallic());
    }
    #[test]
    fn test_cnt_axial_stiffness_positive() {
        let cnt = NanotubeProperties::armchair(10);
        assert!(cnt.axial_stiffness_n() > 0.0);
    }
    #[test]
    fn test_graphene_monolayer_modulus() {
        let g = GrapheneSheet::monolayer();
        assert!((g.youngs_modulus_tpa - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_graphene_effective_modulus_with_vacancies() {
        let mut g = GrapheneSheet::monolayer();
        g.vacancy_fraction = 0.1;
        let e = g.effective_youngs_modulus();
        assert!(e < g.youngs_modulus_tpa);
        assert!(e > 0.0);
    }
    #[test]
    fn test_graphene_fracture_strain() {
        let g = GrapheneSheet::monolayer();
        assert!((g.fracture_strain() - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_graphene_bilayer_lower_conductivity() {
        let mono = GrapheneSheet::monolayer();
        let bi = GrapheneSheet::bilayer();
        assert!(bi.effective_thermal_conductivity() < mono.effective_thermal_conductivity());
    }
    #[test]
    fn test_ply_nu21() {
        let ply = Ply::new(0.125e-3, 0.0, 140e9, 10e9, 5e9, 0.3);
        let nu21 = ply.nu21();
        assert!(nu21 > 0.0 && nu21 < 0.3);
    }
    #[test]
    fn test_ply_reduced_stiffness() {
        let ply = Ply::new(0.125e-3, 0.0, 140e9, 10e9, 5e9, 0.3);
        let [q11, q12, q22, q66] = ply.reduced_stiffness();
        assert!(q11 > 0.0);
        assert!(q12 > 0.0);
        assert!(q22 > 0.0);
        assert!(q66 > 0.0);
        assert!(q11 > q22);
    }
    #[test]
    fn test_laminate_total_thickness() {
        let plies = vec![
            Ply::new(1e-4, 0.0, 140e9, 10e9, 5e9, 0.3),
            Ply::new(1e-4, 90.0, 140e9, 10e9, 5e9, 0.3),
        ];
        let lam = CompositeLaminate::new(plies);
        assert!((lam.total_thickness() - 2e-4).abs() < 1e-15);
    }
    #[test]
    fn test_laminate_a11_positive() {
        let plies = vec![Ply::new(1e-4, 0.0, 140e9, 10e9, 5e9, 0.3)];
        let lam = CompositeLaminate::new(plies);
        assert!(lam.a11() > 0.0);
    }
    #[test]
    fn test_laminate_d11_single_ply_is_q11_t3_over_12() {
        let t = 1.5e-4;
        let ply = Ply::new(t, 0.0, 140e9, 10e9, 5e9, 0.3);
        let q11 = ply.transformed_q11();
        let lam = CompositeLaminate::new(vec![ply]);
        let expected = q11 * t.powi(3) / 12.0;
        assert!((lam.d11() - expected).abs() <= 1.0e-12 * expected);
        assert!(lam.b11().abs() <= 1.0e-12 * q11 * t * t);
    }
    #[test]
    fn test_laminate_d11_two_identical_plies_equals_single_layer() {
        let t = 1.0e-4;
        let ply = Ply::new(t, 0.0, 140e9, 10e9, 5e9, 0.3);
        let q11 = ply.transformed_q11();
        let lam = CompositeLaminate::new(vec![ply.clone(), ply]);
        // Two identical plies behave as one homogeneous layer of thickness 2t.
        let expected = q11 * (2.0 * t).powi(3) / 12.0;
        let d11 = lam.d11();
        assert!(
            (d11 - expected).abs() <= 1.0e-12 * expected,
            "D11 = {d11}, expected {expected}"
        );
        // Summing every ply about its own centroid (no parallel-axis term)
        // would give 2 Q t^3 / 12, a factor 4 too small.
        let centroid_only = 2.0 * q11 * t.powi(3) / 12.0;
        assert!((d11 / centroid_only - 4.0).abs() < 1.0e-12);
    }
    #[test]
    fn test_laminate_d11_of_n_identical_plies_grows_as_n_cubed() {
        let t = 1.25e-4;
        let ply = Ply::new(t, 0.0, 140e9, 10e9, 5e9, 0.3);
        let d11_single = CompositeLaminate::new(vec![ply.clone()]).d11();
        for n in 2_usize..=8 {
            let lam = CompositeLaminate::new(vec![ply.clone(); n]);
            let expected = d11_single * (n as f64).powi(3);
            assert!(
                (lam.d11() - expected).abs() <= 1.0e-12 * expected,
                "n = {n}: D11 = {}, expected {expected}",
                lam.d11()
            );
        }
    }
    #[test]
    fn test_laminate_b11_vanishes_for_symmetric_layup() {
        let t = 1.0e-4;
        let ply_0 = Ply::new(t, 0.0, 140e9, 10e9, 5e9, 0.3);
        let ply_90 = Ply::new(t, 90.0, 140e9, 10e9, 5e9, 0.3);
        let q11 = ply_0.transformed_q11();
        let q22 = ply_90.transformed_q11();
        let lam = CompositeLaminate::new(vec![
            ply_0.clone(),
            ply_90.clone(),
            ply_90.clone(),
            ply_0.clone(),
        ]);
        assert!(
            lam.b11().abs() <= 1.0e-12 * q11 * t * t,
            "B11 = {}",
            lam.b11()
        );
        // D11 = (14 Q11 + 2 Q22) t^3 / 3 for [0/90/90/0] with equal ply thickness.
        let expected = (14.0 * q11 + 2.0 * q22) * t.powi(3) / 3.0;
        assert!(
            (lam.d11() - expected).abs() <= 1.0e-12 * expected,
            "D11 = {}, expected {expected}",
            lam.d11()
        );
    }
    #[test]
    fn test_laminate_cross_ply_d11_and_b11_closed_form() {
        let t = 1.0e-4;
        let ply_0 = Ply::new(t, 0.0, 140e9, 10e9, 5e9, 0.3);
        let ply_90 = Ply::new(t, 90.0, 140e9, 10e9, 5e9, 0.3);
        let q11 = ply_0.transformed_q11();
        let q22 = ply_90.transformed_q11();
        assert!(q11 > q22);
        let lam = CompositeLaminate::new(vec![ply_0.clone(), ply_90.clone()]);
        // [0/90]: mid-plane at the interface, z in {-t, 0, t}.
        let d11_expected = (q11 + q22) * t.powi(3) / 3.0;
        let b11_expected = 0.5 * t * t * (q22 - q11);
        assert!((lam.d11() - d11_expected).abs() <= 1.0e-12 * d11_expected);
        assert!((lam.b11() - b11_expected).abs() <= 1.0e-12 * b11_expected.abs());
        assert!(
            lam.b11() < 0.0,
            "the stiffer 0 ply sits below the mid-plane"
        );
        // Reversing the stack mirrors z: D11 is unchanged, B11 changes sign.
        let reversed = CompositeLaminate::new(vec![ply_90, ply_0]);
        assert!((reversed.d11() - lam.d11()).abs() <= 1.0e-12 * d11_expected);
        assert!((reversed.b11() + lam.b11()).abs() <= 1.0e-12 * b11_expected.abs());
    }
    #[test]
    fn test_laminate_d11_matches_parallel_axis_theorem_for_mixed_stack() {
        let plies = vec![
            Ply::new(1.0e-4, 0.0, 140e9, 10e9, 5e9, 0.3),
            Ply::new(2.0e-4, 45.0, 140e9, 10e9, 5e9, 0.3),
            Ply::new(1.5e-4, 90.0, 140e9, 10e9, 5e9, 0.3),
            Ply::new(0.5e-4, -30.0, 140e9, 10e9, 5e9, 0.3),
        ];
        let lam = CompositeLaminate::new(plies.clone());
        // Independent evaluation: D11 = sum Qbar11 (t^3/12 + t z_c^2), B11 = sum Qbar11 t z_c
        // with z_c the distance of each ply centroid from the laminate mid-plane.
        let total: f64 = plies.iter().map(|p| p.thickness).sum();
        let mut z_bottom = -0.5 * total;
        let mut d11 = 0.0;
        let mut b11 = 0.0;
        for p in &plies {
            let z_c = z_bottom + 0.5 * p.thickness;
            let q = p.transformed_q11();
            d11 += q * (p.thickness.powi(3) / 12.0 + p.thickness * z_c * z_c);
            b11 += q * p.thickness * z_c;
            z_bottom += p.thickness;
        }
        assert!(
            (lam.d11() - d11).abs() <= 1.0e-10 * d11,
            "{} vs {d11}",
            lam.d11()
        );
        let scale = lam.a11() * total;
        assert!(
            (lam.b11() - b11).abs() <= 1.0e-10 * scale,
            "{} vs {b11}",
            lam.b11()
        );
        assert!(lam.b11().abs() > 1.0e-6 * scale, "stack is unsymmetric");
    }
    #[test]
    fn test_fiber_matrix_e1() {
        let fm = FiberMatrix::new(0.6, 230e9, 3.5e9, 0.2, 0.35);
        let e1 = fm.e1();
        assert!((e1 - (0.6 * 230e9 + 0.4 * 3.5e9)).abs() < 1e3);
    }
    #[test]
    fn test_fiber_matrix_e2_halpin_tsai() {
        let fm = FiberMatrix::new(0.6, 230e9, 3.5e9, 0.2, 0.35);
        let e2 = fm.e2_halpin_tsai();
        assert!(e2 > fm.em);
        assert!(e2 < fm.ef);
    }
    #[test]
    fn test_nanoparticle_effective_modulus() {
        let nano = NanoparticleDispersion::new(0.1, 100e9, 3e9, 0.2, 0.35);
        let e = nano.effective_youngs_modulus();
        assert!(e > nano.em);
        assert!(e < nano.ep);
    }
    #[test]
    fn test_interfacial_zone_debonding() {
        let mut iz = InterfacialZone::new(100.0, 50e6, 0.1, 1e-9);
        assert!(!iz.debonded);
        iz.check_debond(60e6);
        assert!(iz.debonded);
    }
    #[test]
    fn test_interfacial_zone_vdw() {
        let iz = InterfacialZone::new(100.0, 50e6, 0.1, 1e-9);
        let f = iz.vdw_force_per_area(0.3);
        assert!(f > 0.0);
    }
    #[test]
    fn test_size_effect_hall_petch() {
        let se = SizeEffect::new(0.5e6, 50e6, 10.0, 1e-6, 100e6);
        let sy = se.hall_petch_strength(1e-6);
        assert!(sy > se.friction_stress);
    }
    #[test]
    fn test_size_effect_weibull_strength_scaling() {
        let se = SizeEffect::new(0.5e6, 50e6, 10.0, 1e-6, 100e6);
        let s_large = se.weibull_strength(10e-6);
        let s_small = se.weibull_strength(1e-7);
        assert!(s_small > s_large);
    }
    #[test]
    fn test_molecular_crystal_modulus() {
        let mc = MolecularCrystal::new([5.0, 5.0, 5.0], 10.0, 4);
        let e = mc.youngs_modulus();
        assert!(e > 0.0);
    }
    #[test]
    fn test_polymer_chain_contour_length() {
        let pc = PolymerChain::new(100, 0.38e-9, 0.5e-9, 300.0);
        assert!((pc.contour_length() - 100.0 * 0.38e-9).abs() < 1e-20);
    }
    #[test]
    fn test_polymer_chain_entropic_force() {
        let pc = PolymerChain::new(100, 0.38e-9, 0.5e-9, 300.0);
        let f = pc.entropic_force(0.1 * pc.contour_length());
        assert!(f > 0.0);
    }
    #[test]
    fn test_polymer_flory_huggins() {
        let pc = PolymerChain::new(100, 0.38e-9, 0.5e-9, 300.0);
        let fe = pc.flory_huggins_free_energy(0.5, 0.5);
        assert!(fe.is_finite());
    }
    #[test]
    fn test_biomaterial_cortical_bone_modulus() {
        let bone = BioMaterial::cortical_bone();
        let e = bone.effective_modulus();
        assert!(e > 1e9);
    }
    #[test]
    fn test_biomaterial_mooney_rivlin() {
        let tissue = BioMaterial::soft_tissue();
        let w = tissue.mooney_rivlin_energy(3.5, 3.2);
        assert!(w > 0.0);
    }
    #[test]
    fn test_thermoelectric_zt() {
        let te = ThermoelectricMaterial::bismuth_telluride();
        let zt = te.zt();
        assert!((zt - 0.8).abs() < 0.5, "ZT={}", zt);
    }
    #[test]
    fn test_thermoelectric_peltier() {
        let te = ThermoelectricMaterial::bismuth_telluride();
        let pi = te.peltier_coefficient();
        assert!((pi - 200e-6 * 300.0).abs() < 1e-6);
    }
    #[test]
    fn test_thermoelectric_efficiency() {
        let te = ThermoelectricMaterial::bismuth_telluride();
        let eta = te.device_efficiency(50.0);
        assert!(eta > 0.0 && eta < 1.0);
    }
}
#[cfg(test)]
mod nano_extended_tests {

    use crate::nano::CntBuckling;

    use crate::nano::DislocationMechanics;

    use crate::nano::GrainBoundaryMechanics;
    use crate::nano::GrapheneElasticStiffness;

    use crate::nano::NanoGrapheneThermal;
    use crate::nano::NanoindentationOliverPharr;

    use crate::nano::PhononTransport;

    use crate::nano::QuantumDotConfinement;

    use crate::nano::SurfaceToVolumeEffects;

    use crate::nano::ThinFilmMechanics;
    use crate::nano::VirialStress;
    #[test]
    fn cnt_buckling_moi_positive() {
        let cnt = CntBuckling::swcnt(1.4, 10.0);
        assert!(cnt.moment_of_inertia() > 0.0);
    }
    #[test]
    fn cnt_euler_buckling_positive() {
        let cnt = CntBuckling::swcnt(1.4, 100.0);
        assert!(cnt.euler_buckling_load() > 0.0);
    }
    #[test]
    fn cnt_shell_buckling_stress_positive() {
        let cnt = CntBuckling::swcnt(1.4, 50.0);
        assert!(cnt.shell_buckling_stress() > 0.0);
    }
    #[test]
    fn cnt_slenderness_longer_is_higher() {
        let cnt_short = CntBuckling::swcnt(1.4, 20.0);
        let cnt_long = CntBuckling::swcnt(1.4, 100.0);
        assert!(cnt_long.slenderness_ratio() > cnt_short.slenderness_ratio());
    }
    #[test]
    fn cnt_critical_strain_positive() {
        let cnt = CntBuckling::swcnt(1.4, 50.0);
        assert!(cnt.critical_strain() > 0.0);
    }
    #[test]
    fn graphene_youngs_modulus_2d_positive() {
        let g = GrapheneElasticStiffness::monolayer();
        assert!(g.youngs_modulus_2d() > 0.0);
    }
    #[test]
    fn graphene_poisson_ratio_between_0_05() {
        let g = GrapheneElasticStiffness::monolayer();
        let nu = g.poisson_ratio();
        assert!(nu > 0.0 && nu < 0.5);
    }
    #[test]
    fn graphene_3d_modulus_near_1tpa() {
        let g = GrapheneElasticStiffness::monolayer();
        let e3d = g.youngs_modulus_3d_gpa();
        assert!(e3d > 500.0 && e3d < 2000.0, "E3D={}", e3d);
    }
    #[test]
    fn graphene_rgo_lower_modulus() {
        let g_perfect = GrapheneElasticStiffness::monolayer();
        let g_rgo = GrapheneElasticStiffness::reduced_graphene_oxide();
        assert!(g_rgo.youngs_modulus_2d() < g_perfect.youngs_modulus_2d());
    }
    #[test]
    fn qdot_confinement_energy_positive() {
        let qd = QuantumDotConfinement::cdse(2.0);
        let e = qd.confinement_energy_ev();
        assert!(e > 0.0);
    }
    #[test]
    fn qdot_bandgap_larger_than_bulk() {
        let qd = QuantumDotConfinement::cdse(2.0);
        assert!(qd.effective_bandgap_ev() > qd.bulk_bandgap_ev);
    }
    #[test]
    fn qdot_emission_wavelength_visible() {
        let qd = QuantumDotConfinement::cdse(2.5);
        let lambda = qd.emission_wavelength_nm();
        assert!(lambda > 400.0 && lambda < 1500.0, "λ={}", lambda);
    }
    #[test]
    fn qdot_smaller_radius_larger_bandgap() {
        let qd_small = QuantumDotConfinement::cdse(1.5);
        let qd_large = QuantumDotConfinement::cdse(4.0);
        assert!(qd_small.effective_bandgap_ev() > qd_large.effective_bandgap_ev());
    }
    #[test]
    fn qdot_bohr_radius_positive() {
        let qd = QuantumDotConfinement::inp(3.0);
        assert!(qd.bohr_radius_nm() > 0.0);
    }
    #[test]
    fn sv_hall_petch_normal_regime() {
        let m = SurfaceToVolumeEffects::nc_copper();
        let sy = m.yield_strength_mpa(100.0);
        assert!(sy > m.sigma_0);
    }
    #[test]
    fn sv_inverse_hp_below_breakdown() {
        let m = SurfaceToVolumeEffects::nc_copper();
        let sy_above = m.yield_strength_mpa(m.d_breakdown_nm + 5.0);
        let sy_below = m.yield_strength_mpa(m.d_breakdown_nm / 2.0);
        assert!(sy_below < sy_above);
    }
    #[test]
    fn sv_capillary_pressure_positive() {
        let m = SurfaceToVolumeEffects::nc_iron();
        let p = m.capillary_pressure_mpa(5.0);
        assert!(p > 0.0);
    }
    #[test]
    fn sv_surface_atom_fraction_increases_with_smaller_particle() {
        let m = SurfaceToVolumeEffects::nc_copper();
        let f_large = m.surface_atom_fraction(10.0);
        let f_small = m.surface_atom_fraction(2.0);
        assert!(f_small > f_large);
    }
    #[test]
    fn nanoindent_contact_area_positive() {
        let ind = NanoindentationOliverPharr::berkovich_diamond();
        let a = ind.contact_area_nm2(200.0, 1.0e-3, 1.0e-2);
        assert!(a > 0.0);
    }
    #[test]
    fn nanoindent_hardness_positive() {
        let ind = NanoindentationOliverPharr::berkovich_diamond();
        let a = ind.contact_area_nm2(200.0, 1.0e-3, 1.0e-2);
        let h = ind.hardness_gpa(1.0e-3, a);
        assert!(h > 0.0);
    }
    #[test]
    fn nanoindent_sample_modulus_positive() {
        let ind = NanoindentationOliverPharr::berkovich_diamond();
        let e_sample = ind.sample_youngs_modulus(100.0, 0.3);
        assert!(e_sample > 0.0);
    }
    #[test]
    fn nanoindent_pile_up_correction_gt_1() {
        let ind = NanoindentationOliverPharr::cube_corner();
        let f_pu = ind.pile_up_correction(150.0, 5.0);
        assert!(f_pu > 1.0);
    }
    #[test]
    fn disloc_peierls_stress_positive() {
        let dm = DislocationMechanics::fcc_copper();
        let tau = dm.peierls_nabarro_stress();
        assert!(tau > 0.0);
    }
    #[test]
    fn disloc_orowan_positive() {
        let dm = DislocationMechanics::fcc_copper();
        let d_sigma = dm.orowan_stress_mpa(100.0);
        assert!(d_sigma > 0.0);
    }
    #[test]
    fn disloc_taylor_hardening_positive() {
        let dm = DislocationMechanics::bcc_iron();
        let d_sigma = dm.taylor_hardening_mpa(1.0e14);
        assert!(d_sigma > 0.0);
    }
    #[test]
    fn disloc_taylor_hardening_copper_is_about_99_mpa() {
        // Annealed-to-moderately-worked Cu, rho = 1e14 m^-2: M alpha G b sqrt(rho)
        // = 3.06 * 0.3 * 42 GPa * 0.256 nm * 1e7 m^-1 = 98.7 MPa.
        let dm = DislocationMechanics::fcc_copper();
        let d_sigma = dm.taylor_hardening_mpa(1.0e14);
        assert!(
            (d_sigma - 98.7).abs() < 0.5,
            "Taylor hardening = {d_sigma} MPa"
        );
        // Independent evaluation entirely in SI units (Pa, m).
        let sigma_pa = 3.06 * 0.3 * 42.0e9 * 0.256e-9 * 1.0e14_f64.sqrt();
        assert!((d_sigma - sigma_pa * 1.0e-6).abs() <= 1.0e-12 * d_sigma);
    }
    #[test]
    fn disloc_taylor_hardening_iron_matches_si_evaluation() {
        let dm = DislocationMechanics::bcc_iron();
        let rho: f64 = 4.0e13;
        let sigma_pa = 3.06 * 0.3 * 82.0e9 * 0.248e-9 * rho.sqrt();
        let d_sigma = dm.taylor_hardening_mpa(rho);
        assert!((d_sigma - sigma_pa * 1.0e-6).abs() <= 1.0e-12 * d_sigma);
        // Tens to a few hundred MPa for typical dislocation densities.
        assert!(d_sigma > 50.0 && d_sigma < 500.0, "{d_sigma} MPa");
    }
    #[test]
    fn disloc_taylor_hardening_scales_with_sqrt_of_density() {
        let dm = DislocationMechanics::fcc_copper();
        let base = dm.taylor_hardening_mpa(1.0e13);
        assert!((dm.taylor_hardening_mpa(4.0e13) / base - 2.0).abs() < 1.0e-12);
        assert!((dm.taylor_hardening_mpa(1.0e15) / base - 10.0).abs() < 1.0e-12);
        assert_eq!(dm.taylor_hardening_mpa(0.0), 0.0);
    }
    #[test]
    fn disloc_line_energy_positive() {
        let dm = DislocationMechanics::bcc_iron();
        assert!(dm.line_energy() > 0.0);
    }
    #[test]
    fn gb_read_shockley_zero_at_zero_angle() {
        let gb = GrainBoundaryMechanics::aluminum_gb();
        let e = gb.read_shockley_energy(0.0);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn gb_read_shockley_positive_small_angle() {
        let gb = GrainBoundaryMechanics::aluminum_gb();
        let e = gb.read_shockley_energy(5.0);
        assert!(e > 0.0);
    }
    #[test]
    fn gb_diffusivity_increases_with_temp() {
        let gb = GrainBoundaryMechanics::nickel_gb();
        let d_low = gb.gb_diffusivity(800.0);
        let d_high = gb.gb_diffusivity(1200.0);
        assert!(d_high > d_low);
    }
    #[test]
    fn gb_coble_creep_positive() {
        let gb = GrainBoundaryMechanics::aluminum_gb();
        let eps_dot = gb.coble_creep_rate(10.0e6, 1.0e-6, 773.15, 1.66e-29);
        assert!(eps_dot > 0.0);
    }
    #[test]
    fn phonon_kn_large_for_small_length() {
        let si = PhononTransport::silicon();
        let kn = si.knudsen_number(10.0);
        assert!(kn > 1.0);
    }
    #[test]
    fn phonon_effective_conductivity_reduced_at_small_scale() {
        let si = PhononTransport::silicon();
        let k_large = si.effective_conductivity(10_000.0);
        let k_small = si.effective_conductivity(10.0);
        assert!(k_small < k_large);
    }
    #[test]
    fn phonon_ballistic_heat_flux_positive() {
        let ge = PhononTransport::germanium();
        let q = ge.ballistic_heat_flux(10.0);
        assert!(q > 0.0);
    }
    #[test]
    fn phonon_conductivity_drops_above_debye() {
        let si = PhononTransport::silicon();
        let k_debye = si.conductivity_at_temp(si.debye_temp);
        let k_high = si.conductivity_at_temp(si.debye_temp * 2.0);
        assert!(k_high < k_debye);
    }
    #[test]
    fn thin_film_biaxial_modulus_positive() {
        let tf = ThinFilmMechanics::tin_on_silicon(200.0);
        assert!(tf.biaxial_modulus() > 0.0);
    }
    #[test]
    fn thin_film_thermal_stress_nonzero() {
        let tf = ThinFilmMechanics::cu_on_sio2(500.0);
        let sigma = tf.thermal_mismatch_stress_mpa(100.0);
        assert!(sigma != 0.0);
    }
    #[test]
    fn thin_film_stoney_curvature_finite() {
        let tf = ThinFilmMechanics::tin_on_silicon(300.0);
        let kappa = tf.stoney_curvature(500.0);
        assert!(kappa.is_finite());
    }
    #[test]
    fn thin_film_channel_crack_erg_positive() {
        let tf = ThinFilmMechanics::tin_on_silicon(200.0);
        let g = tf.channel_crack_erg(1000.0);
        assert!(g > 0.0);
    }
    #[test]
    fn thin_film_delamination_energy_positive() {
        let tf = ThinFilmMechanics::cu_on_sio2(500.0);
        let g = tf.delamination_energy(500.0);
        assert!(g > 0.0);
    }
    #[test]
    fn virial_kinetic_stress_positive() {
        let vs = VirialStress::silicon_rve(10.0);
        let p = vs.kinetic_stress_pa();
        assert!(p > 0.0);
    }
    #[test]
    fn virial_silicon_rve_atom_count_is_cubic_in_side() {
        // 8 atoms per conventional cell of edge 0.543 nm.
        assert_eq!(VirialStress::silicon_rve(0.543).n_atoms, 8);
        assert_eq!(VirialStress::silicon_rve(5.43).n_atoms, 8000);
        assert_eq!(VirialStress::silicon_rve(10.86).n_atoms, 64_000);
    }
    #[test]
    fn virial_kinetic_stress_equals_n_kb_t_over_v() {
        let vs = VirialStress::silicon_rve(5.43);
        let kb = 1.380649e-23;
        let volume = (5.43e-9_f64).powi(3);
        let expected = 8000.0 * kb * 300.0 / volume;
        let p = vs.kinetic_stress_pa();
        assert!(
            (p - expected).abs() / expected < 1e-12,
            "P_kin = {p}, N kB T / V = {expected}"
        );
        // 8 kB T / a^3 with a = 0.543 nm is about 2.07e8 Pa.
        assert!((p - 2.07e8).abs() / 2.07e8 < 0.01, "P_kin = {p}");
    }
    #[test]
    fn virial_kinetic_stress_is_independent_of_atom_mass_and_rve_size() {
        let mut light = VirialStress::silicon_rve(5.43);
        let p_ref = light.kinetic_stress_pa();
        light.atom_mass *= 7.0;
        assert!((light.kinetic_stress_pa() - p_ref).abs() / p_ref < 1e-15);
        // Atom number density is the lattice value, so the pressure is intensive.
        let big = VirialStress::silicon_rve(10.86);
        assert!((big.kinetic_stress_pa() - p_ref).abs() / p_ref < 1e-9);
    }
    #[test]
    fn hardy_stress_with_zero_velocities_equals_potential_virial_stress() {
        let vs = VirialStress::silicon_rve(5.0);
        let pairs = vec![
            (-1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0),
            (2.0e-10, 1.0e-10, 0.0, -0.25e-9, 0.1e-9, 0.0),
            (-5.0e-10, 0.0, 3.0e-10, 0.2e-9, 0.0, -0.1e-9),
        ];
        let velocities = vec![[0.0; 3]; 4];
        let hardy = vs.hardy_stress_xx(&velocities, &pairs);
        let virial = vs.potential_stress_xx(&pairs);
        assert!(virial != 0.0);
        assert!((hardy - virial).abs() <= 1e-15 * virial.abs());
    }
    #[test]
    fn hardy_stress_kinetic_term_uses_thermal_velocities() {
        let vs = VirialStress::silicon_rve(5.0);
        // Two atoms with opposite x-velocities: zero mean, so both are thermal.
        let velocities = vec![[100.0, 0.0, 0.0], [-100.0, 0.0, 0.0]];
        let sigma = vs.hardy_stress_xx(&velocities, &[]);
        let expected = -2.0 * vs.atom_mass * 1.0e4 / vs.volume_m3;
        assert!(
            (sigma - expected).abs() / expected.abs() < 1e-12,
            "sigma_xx = {sigma}, expected {expected}"
        );
        assert!(sigma < 0.0, "kinetic stress is compressive (negative)");
    }
    #[test]
    fn hardy_stress_is_invariant_under_uniform_translation() {
        let vs = VirialStress::silicon_rve(5.0);
        let pairs = vec![(-1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0)];
        let base = vec![[120.0, 5.0, -3.0], [-40.0, 1.0, 2.0], [-80.0, -6.0, 1.0]];
        let shifted: Vec<[f64; 3]> = base.iter().map(|v| [v[0] + 250.0, v[1], v[2]]).collect();
        let s0 = vs.hardy_stress_xx(&base, &pairs);
        let s1 = vs.hardy_stress_xx(&shifted, &pairs);
        assert!((s0 - s1).abs() / s0.abs() < 1e-9, "{s0} vs {s1}");
        // A rigidly translating body without pair forces carries no stress at all.
        let rigid = vec![[300.0, 0.0, 0.0]; 5];
        assert_eq!(vs.hardy_stress_xx(&rigid, &[]), 0.0);
    }
    #[test]
    fn hardy_stress_is_sum_of_kinetic_and_configurational_parts() {
        let vs = VirialStress::silicon_rve(5.0);
        let pairs = vec![
            (-1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0),
            (4.0e-10, 0.0, 0.0, -0.2e-9, 0.0, 0.0),
        ];
        let velocities = vec![[60.0, 0.0, 0.0], [-20.0, 0.0, 0.0], [-40.0, 0.0, 0.0]];
        let total = vs.hardy_stress_xx(&velocities, &pairs);
        let kinetic_only = vs.hardy_stress_xx(&velocities, &[]);
        let virial_only = vs.potential_stress_xx(&pairs);
        assert!((total - (kinetic_only + virial_only)).abs() <= 1e-12 * total.abs());
        // Attractive pair (f . r < 0) is tensile, repulsive pair compressive.
        let attractive = vs.potential_stress_xx(&[(-1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0)]);
        let repulsive = vs.potential_stress_xx(&[(1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0)]);
        assert!(attractive > 0.0 && repulsive < 0.0);
    }
    #[test]
    fn hardy_kinetic_stress_matches_equipartition_pressure() {
        // 512-atom RVE (4 x 4 x 4 conventional cells) at 300 K whose x-velocities are
        // +-v_th with m v_th^2 = k_B T: the Hardy kinetic stress is exactly -N k_B T / V.
        let vs = VirialStress::silicon_rve(4.0 * 0.543);
        assert_eq!(vs.n_atoms, 512);
        let v_th = (1.380649e-23 * vs.temperature / vs.atom_mass).sqrt();
        let velocities: Vec<[f64; 3]> = (0..vs.n_atoms)
            .map(|i| [if i % 2 == 0 { v_th } else { -v_th }, 0.0, 0.0])
            .collect();
        let sigma = vs.hardy_stress_xx(&velocities, &[]);
        let expected = -vs.kinetic_stress_pa();
        assert!(
            (sigma - expected).abs() / expected.abs() < 1e-12,
            "sigma_xx = {sigma}, -P_kin = {expected}"
        );
    }
    #[test]
    fn virial_cauchy_born_identity_gives_zero_strain() {
        let f_identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let e = VirialStress::cauchy_born_strain(&f_identity);
        for &ei in &e {
            assert!(ei.abs() < 1.0e-14, "strain component {} != 0", ei);
        }
    }
    #[test]
    fn virial_potential_stress_finite() {
        let vs = VirialStress::silicon_rve(5.0);
        let pairs = vec![(1.0e-9, 0.0, 0.0, 0.3e-9, 0.0, 0.0)];
        let sigma = vs.potential_stress_xx(&pairs);
        assert!(sigma.is_finite());
    }
    #[test]
    fn graphene_thermal_suspended_higher_k() {
        let susp = NanoGrapheneThermal::suspended();
        let sio2 = NanoGrapheneThermal::on_sio2();
        let k_susp = susp.effective_conductivity(300.0);
        let k_sio2 = sio2.effective_conductivity(300.0);
        assert!(k_susp > k_sio2);
    }
    #[test]
    fn graphene_thermal_eff_mfp_less_than_intrinsic() {
        let g = NanoGrapheneThermal::on_sio2();
        let mfp_eff = g.effective_mfp_nm(300.0);
        assert!(mfp_eff < 300.0);
    }
    #[test]
    fn graphene_thermal_substrate_resistance_finite() {
        let g = NanoGrapheneThermal::on_sio2();
        let r = g.substrate_thermal_resistance(1.0);
        assert!(r.is_finite() && r > 0.0);
    }
}
