//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metal_alloy_md::AlloySystem;
    use crate::metal_alloy_md::DiffusionTracker;
    use crate::metal_alloy_md::EamParams;
    use crate::metal_alloy_md::MetalAtom;
    use crate::metal_alloy_md::MixingRule;
    use crate::metal_alloy_md::types::CrystalStructure;
    use crate::metal_alloy_md::types::SimBox;
    #[test]
    fn test_eam_copper_pair_at_equilibrium() {
        let cu = EamParams::copper();
        let phi = cu.pair_potential(cu.re);
        assert!(phi.abs() < 5.0, "pair potential at re = {:.6}", phi);
    }
    #[test]
    fn test_eam_pair_zero_beyond_cutoff() {
        let cu = EamParams::copper();
        assert_eq!(cu.pair_potential(cu.cutoff + 1.0), 0.0);
        assert_eq!(cu.electron_density(cu.cutoff + 1.0), 0.0);
    }
    #[test]
    fn test_eam_embedding_negative() {
        let cu = EamParams::copper();
        let f = cu.embedding_energy(1.0);
        assert!(f < 0.0, "embedding energy should be negative, got {:.6}", f);
    }
    #[test]
    fn test_eam_embedding_deriv() {
        let cu = EamParams::copper();
        let rho = 2.0;
        let dr = 1e-6;
        let numerical =
            (cu.embedding_energy(rho + dr) - cu.embedding_energy(rho - dr)) / (2.0 * dr);
        let analytical = cu.embedding_deriv(rho);
        assert!(
            (numerical - analytical).abs() < 1e-4,
            "embedding deriv: numerical={:.6} analytical={:.6}",
            numerical,
            analytical
        );
    }
    #[test]
    fn test_mixing_rules_symmetric() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let r = 2.6;
        for rule in [
            MixingRule::Geometric,
            MixingRule::Arithmetic,
            MixingRule::Johnson,
        ] {
            let v1 = cross_pair_potential(&cu, &ni, r, rule);
            let v2 = cross_pair_potential(&ni, &cu, r, rule);
            assert!(
                (v1 - v2).abs() < 1e-10,
                "mixing {:?} not symmetric: {:.6} vs {:.6}",
                rule,
                v1,
                v2
            );
        }
    }
    #[test]
    fn test_mixing_arithmetic_average() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let r = 2.7;
        let avg = 0.5 * (cu.pair_potential(r) + ni.pair_potential(r));
        let mixed = cross_pair_potential(&cu, &ni, r, MixingRule::Arithmetic);
        assert!((avg - mixed).abs() < 1e-12);
    }
    #[test]
    fn test_sim_box_min_image() {
        let b = SimBox::cubic(10.0);
        let dr = b.min_image([9.0, -8.0, 3.0]);
        assert!((dr[0] - (-1.0)).abs() < 1e-10);
        assert!((dr[1] - 2.0).abs() < 1e-10);
        assert!((dr[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_sim_box_wrap() {
        let b = SimBox::cubic(10.0);
        let p = b.wrap([12.0, -3.0, 5.0]);
        assert!(p[0] >= 0.0 && p[0] < 10.0);
        assert!(p[1] >= 0.0 && p[1] < 10.0);
        assert!(p[2] >= 0.0 && p[2] < 10.0);
    }
    #[test]
    fn test_fcc_generation_count() {
        let positions = generate_crystal_positions(3.6, 2, CrystalStructure::FCC);
        assert_eq!(positions.len(), 32);
    }
    #[test]
    fn test_bcc_generation_count() {
        let positions = generate_crystal_positions(2.87, 2, CrystalStructure::BCC);
        assert_eq!(positions.len(), 16);
    }
    #[test]
    fn test_hcp_generation_count() {
        let positions = generate_crystal_positions(3.2, 2, CrystalStructure::HCP);
        assert_eq!(positions.len(), 16);
    }
    #[test]
    fn test_alloy_system_basic() {
        let cu = EamParams::copper();
        let a = cu.lattice_a;
        let n = 2;
        let box_len = n as f64 * a;
        let mut sys = AlloySystem::new(
            vec![cu],
            SimBox::cubic(box_len),
            MixingRule::Geometric,
            0.001,
        );
        let positions = generate_crystal_positions(a, n, CrystalStructure::FCC);
        for pos in &positions {
            sys.add_atom(*pos, 0);
        }
        assert_eq!(sys.num_atoms(), 32);
        sys.compute_forces();
        let pe = sys.potential_energy();
        assert!(pe < 0.0, "PE should be negative, got {:.6}", pe);
    }
    #[test]
    fn test_binary_alloy_creation() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let a = 3.57;
        let box_len = 2.0 * a;
        let mut sys = AlloySystem::new(
            vec![cu, ni],
            SimBox::cubic(box_len),
            MixingRule::Johnson,
            0.001,
        );
        let positions = generate_crystal_positions(a, 2, CrystalStructure::FCC);
        for (idx, pos) in positions.iter().enumerate() {
            sys.add_atom(*pos, idx % 2);
        }
        sys.compute_forces();
        let pe = sys.potential_energy();
        assert!(pe.is_finite(), "PE should be finite");
    }
    #[test]
    fn test_warren_cowley_random() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let a = 3.57;
        let box_len = 3.0 * a;
        let mut sys = AlloySystem::new(
            vec![cu, ni],
            SimBox::cubic(box_len),
            MixingRule::Geometric,
            0.001,
        );
        let positions = generate_crystal_positions(a, 3, CrystalStructure::FCC);
        for (idx, pos) in positions.iter().enumerate() {
            sys.add_atom(*pos, idx % 2);
        }
        let sro = warren_cowley_sro(&sys, 0, 1);
        assert!(sro.is_finite());
        assert!(sro.abs() < 2.0, "SRO out of range: {:.6}", sro);
    }
    #[test]
    fn test_diffusion_tracker_initial_msd_zero() {
        let cu = EamParams::copper();
        let mut sys = AlloySystem::new(vec![cu], SimBox::cubic(7.23), MixingRule::Geometric, 0.001);
        sys.add_atom([1.0, 2.0, 3.0], 0);
        sys.add_atom([4.0, 5.0, 6.0], 0);
        let tracker = DiffusionTracker::new(&sys, None);
        let msd = tracker.msd(&sys);
        assert!(
            msd.abs() < 1e-12,
            "initial MSD should be ~0, got {:.6}",
            msd
        );
    }
    #[test]
    fn test_diffusion_tracker_after_displacement() {
        let cu = EamParams::copper();
        let mut sys = AlloySystem::new(vec![cu], SimBox::cubic(20.0), MixingRule::Geometric, 0.001);
        sys.add_atom([5.0, 5.0, 5.0], 0);
        let mut tracker = DiffusionTracker::new(&sys, None);
        sys.atoms[0].position = [6.0, 5.0, 5.0];
        tracker.update(&sys);
        let msd = tracker.msd(&sys);
        assert!(
            (msd - 1.0).abs() < 1e-10,
            "MSD should be 1.0, got {:.6}",
            msd
        );
    }
    #[test]
    fn test_mole_fractions() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let mut sys = AlloySystem::new(
            vec![cu, ni],
            SimBox::cubic(10.0),
            MixingRule::Geometric,
            0.001,
        );
        for i in 0..10 {
            sys.add_atom([i as f64, 0.0, 0.0], 0);
        }
        for i in 0..5 {
            sys.add_atom([i as f64, 1.0, 0.0], 1);
        }
        let fracs = mole_fractions(&sys);
        assert!((fracs[0] - 10.0 / 15.0).abs() < 1e-10);
        assert!((fracs[1] - 5.0 / 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_peierls_stress_positive() {
        let sigma = peierls_stress(50.0, 0.33, 2.1, 2.56);
        assert!(
            sigma > 0.0,
            "Peierls stress should be positive, got {:.6}",
            sigma
        );
    }
    #[test]
    fn test_peierls_stress_from_eam_positive() {
        let cu = EamParams::copper();
        let sigma = peierls_stress_from_eam(&cu);
        assert!(sigma > 0.0, "Peierls stress should be positive");
        assert!(sigma.is_finite());
    }
    #[test]
    fn test_dislocation_core_structure() {
        let cu = EamParams::copper();
        let core = analyze_dislocation_core(&cu, 3);
        assert!(core.burgers_vector[0] > 0.0);
        assert!(core.core_width > 0.0);
        assert!(!core.displacements.is_empty());
    }
    #[test]
    fn test_grain_boundary_energy_positive() {
        let cu = EamParams::copper();
        let e_gb = grain_boundary_energy(&cu, 15.0, 3);
        assert!(e_gb >= 0.0, "GB energy should be >= 0, got {:.6}", e_gb);
        assert!(e_gb.is_finite());
    }
    #[test]
    fn test_stacking_fault_curve_monotonic_start() {
        let cu = EamParams::copper();
        let curve = stacking_fault_energy_curve(&cu, 2, 5);
        assert_eq!(curve.len(), 5);
        assert!(
            curve[0].1.abs() < 0.1,
            "SF energy at zero displacement should be ~0, got {:.6}",
            curve[0].1
        );
    }
    #[test]
    fn test_partial_rdf_basic() {
        let cu = EamParams::copper();
        let a = cu.lattice_a;
        let mut sys = AlloySystem::new(
            vec![cu],
            SimBox::cubic(2.0 * a),
            MixingRule::Geometric,
            0.001,
        );
        let positions = generate_crystal_positions(a, 2, CrystalStructure::FCC);
        for pos in &positions {
            sys.add_atom(*pos, 0);
        }
        let (centers, g) = partial_rdf(&sys, 0, 0, 5.0, 50);
        assert_eq!(centers.len(), 50);
        assert_eq!(g.len(), 50);
        let nn_dist = a / 2.0_f64.sqrt();
        let peak_bin = (nn_dist / 0.1) as usize;
        if peak_bin < 50 {
            assert!(g[peak_bin] > 0.0, "RDF should have peak near nn distance");
        }
    }
    #[test]
    fn test_elastic_constants_positive() {
        let cu = EamParams::copper();
        let (c11, _c12, c44) = elastic_constants_eam(&cu);
        assert!(c11.is_finite());
        assert!(c44.is_finite());
        assert!(c11 != 0.0, "C11 should be nonzero");
    }
    #[test]
    fn test_heat_of_mixing_finite() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let dh = heat_of_mixing(&cu, &ni, 0.5, 2, MixingRule::Johnson);
        assert!(dh.is_finite(), "heat of mixing should be finite");
    }
    #[test]
    fn test_thermal_expansion_positive() {
        let cu = EamParams::copper();
        let alpha = thermal_expansion_coefficient(&cu, 300.0, 600.0, 2);
        assert!(
            alpha > 0.0 && alpha.is_finite(),
            "thermal expansion should be positive finite, got {:.6}",
            alpha
        );
    }
    #[test]
    fn test_coordination_numbers() {
        let cu = EamParams::copper();
        let a = cu.lattice_a;
        let mut sys = AlloySystem::new(
            vec![cu],
            SimBox::cubic(3.0 * a),
            MixingRule::Geometric,
            0.001,
        );
        let positions = generate_crystal_positions(a, 3, CrystalStructure::FCC);
        for pos in &positions {
            sys.add_atom(*pos, 0);
        }
        let nn_cut = a / 2.0_f64.sqrt() + 0.3;
        let cn = coordination_numbers(&sys, nn_cut);
        assert!(
            cn[0][0] > 8.0,
            "FCC coordination should be ~12, got {:.1}",
            cn[0][0]
        );
    }
    #[test]
    fn test_metal_atom_kinetic_energy() {
        let mut atom = MetalAtom::new([0.0; 3], 0);
        atom.velocity = [1.0, 0.0, 0.0];
        let ke = atom.kinetic_energy(63.546);
        assert!(ke > 0.0);
    }
    #[test]
    fn test_vacancy_migration_energy_nonnegative() {
        let cu = EamParams::copper();
        let e = vacancy_migration_energy(
            &[cu],
            3.615,
            CrystalStructure::FCC,
            0,
            1,
            5,
            MixingRule::Geometric,
        );
        assert!(
            e >= 0.0,
            "vacancy migration energy should be >= 0, got {:.6}",
            e
        );
        assert!(e.is_finite());
    }
    #[test]
    fn test_vector_helpers() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let s = vadd(a, b);
        assert!((s[0] - 5.0).abs() < 1e-12);
        let d = vsub(b, a);
        assert!((d[0] - 3.0).abs() < 1e-12);
        let dot = vdot(a, b);
        assert!((dot - 32.0).abs() < 1e-12);
        let cross = vcross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((cross[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_common_neighbor_analysis_runs() {
        let cu = EamParams::copper();
        let a = cu.lattice_a;
        let mut sys = AlloySystem::new(
            vec![cu],
            SimBox::cubic(2.0 * a),
            MixingRule::Geometric,
            0.001,
        );
        let positions = generate_crystal_positions(a, 2, CrystalStructure::FCC);
        for pos in &positions {
            sys.add_atom(*pos, 0);
        }
        let labels = common_neighbor_analysis(&sys, a / 2.0_f64.sqrt() + 0.3);
        assert_eq!(labels.len(), sys.num_atoms());
    }
}
