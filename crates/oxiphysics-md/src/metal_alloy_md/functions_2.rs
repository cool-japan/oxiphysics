//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::functions::*;
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
    fn test_heat_of_mixing_finite() {
        let cu = EamParams::copper();
        let ni = EamParams::nickel();
        let dh = heat_of_mixing(&cu, &ni, 0.5, 2, MixingRule::Johnson);
        assert!(dh.is_finite(), "heat of mixing should be finite");
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

    // ------------------------------------------------------------------
    // Forces, units and energy conservation
    // ------------------------------------------------------------------

    /// A slightly distorted four-atom tetrahedron (all pair distances well
    /// inside the 5.3-5.5 Å cutoffs and well beyond the 0.5 Å inner cut) whose
    /// atoms are wrapped into a cubic box of edge `box_len`, with its corner at
    /// `origin`.
    fn tetrahedron_cluster(
        species: Vec<EamParams>,
        assignment: [usize; 4],
        rule: MixingRule,
        box_len: f64,
        origin: f64,
    ) -> AlloySystem {
        let offsets = [
            [0.0, 0.0, 0.0],
            [2.62, 0.05, -0.03],
            [1.28, 2.27, 0.04],
            [1.35, 0.78, 2.15],
        ];
        let mut sys = AlloySystem::new(species, SimBox::cubic(box_len), rule, 0.001);
        for (offset, &s) in offsets.iter().zip(assignment.iter()) {
            let p = [origin + offset[0], origin + offset[1], origin + offset[2]];
            let wrapped = sys.sim_box.wrap(p);
            sys.add_atom(wrapped, s);
        }
        sys
    }

    /// Largest deviation between the analytic forces and a central finite
    /// difference of the potential energy (step 1e-5 Å).
    fn max_force_error_vs_finite_difference(sys: &AlloySystem) -> f64 {
        let h = 1.0e-5;
        let mut reference = sys.clone();
        reference.compute_forces();
        let mut worst = 0.0_f64;
        for i in 0..sys.atoms.len() {
            for d in 0..3 {
                let mut plus = sys.clone();
                plus.atoms[i].position[d] += h;
                plus.compute_forces();
                let e_plus = plus.potential_energy();
                let mut minus = sys.clone();
                minus.atoms[i].position[d] -= h;
                minus.compute_forces();
                let e_minus = minus.potential_energy();
                let fd = -(e_plus - e_minus) / (2.0 * h);
                worst = worst.max((reference.atoms[i].force[d] - fd).abs());
            }
        }
        worst
    }

    /// Force on the second atom of a dimer along the (+x) bond direction.
    fn dimer_axial_force(species: &[EamParams], pair: [usize; 2], r: f64) -> f64 {
        let mut sys = AlloySystem::new(
            species.to_vec(),
            SimBox::cubic(30.0),
            MixingRule::Geometric,
            0.001,
        );
        sys.add_atom([10.0, 10.0, 10.0], pair[0]);
        sys.add_atom([10.0 + r, 10.0, 10.0], pair[1]);
        sys.compute_forces();
        sys.atoms[1].force[0]
    }

    /// Separation at which the axial force of a copper dimer changes sign.
    fn copper_dimer_equilibrium() -> f64 {
        let species = [EamParams::copper()];
        let (mut lo, mut hi) = (0.8, 2.0);
        assert!(dimer_axial_force(&species, [0, 0], lo) > 0.0);
        assert!(dimer_axial_force(&species, [0, 0], hi) < 0.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if dimer_axial_force(&species, [0, 0], mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    #[test]
    fn test_eam_forces_match_energy_gradient_single_species() {
        let sys = tetrahedron_cluster(
            vec![EamParams::copper()],
            [0, 0, 0, 0],
            MixingRule::Geometric,
            30.0,
            10.0,
        );
        let err = max_force_error_vs_finite_difference(&sys);
        assert!(err < 1.0e-5, "max |F - (-dE/dx)| = {err:e}");
    }

    #[test]
    fn test_eam_dimer_force_at_two_angstrom_matches_central_difference() {
        // Two Cu atoms 2.0 Å apart: the embedding term dominates the pair term, so a
        // sign error in either of them shows up directly in the force.
        let dimer = |x1: f64| {
            let mut sys = AlloySystem::new(
                vec![EamParams::copper()],
                SimBox::cubic(30.0),
                MixingRule::Geometric,
                0.001,
            );
            sys.add_atom([10.0, 10.0, 10.0], 0);
            sys.add_atom([x1, 10.0, 10.0], 0);
            sys
        };
        let h = 1.0e-6;
        let mut sys = dimer(12.0);
        sys.compute_forces();
        let analytic = sys.atoms[1].force[0];
        let numeric =
            -(dimer(12.0 + h).potential_energy() - dimer(12.0 - h).potential_energy()) / (2.0 * h);
        assert!(
            analytic.abs() > 1.0,
            "force {analytic} eV/Å must be non-trivial"
        );
        assert!(
            ((analytic - numeric) / numeric).abs() < 1.0e-6,
            "analytic {analytic} vs central difference {numeric}"
        );
        // Newton's third law, and no transverse force.
        assert!((sys.atoms[0].force[0] + analytic).abs() < 1.0e-12);
        assert!(sys.atoms[1].force[1].abs() < 1.0e-12 && sys.atoms[1].force[2].abs() < 1.0e-12);
    }

    #[test]
    fn test_eam_forces_match_energy_gradient_across_periodic_boundary() {
        // Corner 11.0 Å inside a 12 Å box: the cluster straddles the boundary.
        let sys = tetrahedron_cluster(
            vec![EamParams::copper()],
            [0, 0, 0, 0],
            MixingRule::Geometric,
            12.0,
            11.0,
        );
        let crossing = sys
            .atoms
            .iter()
            .any(|a| a.position.iter().any(|&c| c < 2.0));
        let remaining = sys
            .atoms
            .iter()
            .any(|a| a.position.iter().any(|&c| c > 10.0));
        assert!(crossing && remaining, "cluster must straddle the boundary");
        let err = max_force_error_vs_finite_difference(&sys);
        assert!(err < 1.0e-5, "max |F - (-dE/dx)| = {err:e}");
    }

    #[test]
    fn test_eam_forces_match_energy_gradient_binary_alloy_all_mixing_rules() {
        for rule in [
            MixingRule::Geometric,
            MixingRule::Arithmetic,
            MixingRule::Johnson,
        ] {
            let sys = tetrahedron_cluster(
                vec![EamParams::copper(), EamParams::nickel()],
                [0, 1, 1, 0],
                rule,
                30.0,
                10.0,
            );
            let err = max_force_error_vs_finite_difference(&sys);
            assert!(err < 1.0e-5, "{rule:?}: max |F - (-dE/dx)| = {err:e}");
        }
    }

    #[test]
    fn test_eam_net_force_vanishes() {
        let mut sys = tetrahedron_cluster(
            vec![EamParams::copper(), EamParams::nickel()],
            [0, 1, 0, 1],
            MixingRule::Johnson,
            30.0,
            10.0,
        );
        sys.compute_forces();
        for d in 0..3 {
            let total: f64 = sys.atoms.iter().map(|a| a.force[d]).sum();
            assert!(total.abs() < 1.0e-10, "net force component {d} = {total:e}");
        }
    }

    #[test]
    fn test_eam_repulsive_pair_pushes_atoms_apart() {
        // Without the embedding term the pair potential is purely repulsive
        // below its zero crossing at ~1.99 Å for copper.
        let mut cu = EamParams::copper();
        cu.f0 = 0.0;
        let r = 1.8;
        assert!(cu.pair_force(r) < 0.0, "phi must decrease with r here");
        let mut sys = AlloySystem::new(
            vec![cu.clone()],
            SimBox::cubic(30.0),
            MixingRule::Geometric,
            0.001,
        );
        sys.add_atom([10.0, 10.0, 10.0], 0);
        sys.add_atom([10.0 + r, 10.0, 10.0], 0);
        sys.compute_forces();
        let expected = cu.pair_force(r);
        // atom 0 sits at smaller x, so a repulsive pair pushes it towards -x.
        assert!(sys.atoms[0].force[0] < 0.0, "{:?}", sys.atoms[0].force);
        assert!(sys.atoms[1].force[0] > 0.0, "{:?}", sys.atoms[1].force);
        assert!((sys.atoms[0].force[0] - expected).abs() < 1.0e-12);
        assert!((sys.atoms[1].force[0] + expected).abs() < 1.0e-12);
        assert!(sys.atoms[0].force[1].abs() < 1.0e-12 && sys.atoms[0].force[2].abs() < 1.0e-12);
    }

    #[test]
    fn test_eam_dimer_force_changes_sign_around_equilibrium() {
        let req = copper_dimer_equilibrium();
        let species = [EamParams::copper()];
        assert!(dimer_axial_force(&species, [0, 0], req - 0.05) > 0.0);
        assert!(dimer_axial_force(&species, [0, 0], req + 0.05) < 0.0);
        assert!(req > 0.9 && req < 1.2, "equilibrium separation {req}");
    }

    #[test]
    fn test_potential_energy_does_not_read_stale_density_cache() {
        let mut sys = tetrahedron_cluster(
            vec![EamParams::copper(), EamParams::nickel()],
            [0, 1, 0, 1],
            MixingRule::Johnson,
            30.0,
            10.0,
        );
        sys.compute_forces();
        let e_fresh = sys.potential_energy();
        // A zeroed cache (as in a system that never computed forces) ...
        for atom in &mut sys.atoms {
            atom.rho_bar = 0.0;
        }
        let e_zero_cache = sys.potential_energy();
        assert!((e_fresh - e_zero_cache).abs() < 1.0e-12);
        // ... and a cache that is stale after the atoms moved.
        sys.compute_forces();
        sys.atoms[2].position[0] += 0.07;
        let e_stale = sys.potential_energy();
        let mut reference = sys.clone();
        reference.compute_forces();
        let e_reference = reference.potential_energy();
        assert!(
            (e_stale - e_reference).abs() < 1.0e-12,
            "{e_stale} vs {e_reference}"
        );
        assert!((e_stale - e_fresh).abs() > 1.0e-6, "the move must matter");
    }

    #[test]
    fn test_kinetic_energy_unit_conversion_gives_equipartition() {
        // v_x^2 = kT/m in amu·Å²/ps² units gives KE = kT/2 exactly.
        let kb_ev = 8.617333262e-5;
        let t = 300.0;
        let mass = 63.546;
        let v = (kb_ev * t * EV_PER_ANGSTROM_AMU_TO_ANGSTROM_PER_PS2 / mass).sqrt();
        let mut atom = MetalAtom::new([0.0; 3], 0);
        atom.velocity = [v, 0.0, 0.0];
        let ke = atom.kinetic_energy(mass);
        assert!(
            (ke - 0.5 * kb_ev * t).abs() < 1.0e-9 * ke,
            "KE = {ke}, expected {}",
            0.5 * kb_ev * t
        );
        // Thermal speed of copper at 300 K is a few hundred m/s = a few Å/ps.
        assert!(v > 1.5 && v < 3.0, "v = {v} Å/ps");
    }

    #[test]
    fn test_step_applies_species_mass_and_unit_conversion() {
        let mut sys = AlloySystem::new(
            vec![EamParams::copper(), EamParams::nickel()],
            SimBox::cubic(30.0),
            MixingRule::Geometric,
            0.001,
        );
        sys.add_atom([10.0, 10.0, 10.0], 0);
        sys.add_atom([11.25, 10.0, 10.0], 1);
        sys.compute_forces();
        let f_start = [sys.atoms[0].force, sys.atoms[1].force];
        assert!(f_start[0][0].abs() > 1.0e-3, "need a non-trivial force");
        sys.step();
        let masses = [sys.species[0].mass, sys.species[1].mass];
        for i in 0..2 {
            let f_mean = [
                0.5 * (f_start[i][0] + sys.atoms[i].force[0]),
                0.5 * (f_start[i][1] + sys.atoms[i].force[1]),
                0.5 * (f_start[i][2] + sys.atoms[i].force[2]),
            ];
            for (d, f_component) in f_mean.iter().enumerate() {
                let expected =
                    f_component * sys.dt * EV_PER_ANGSTROM_AMU_TO_ANGSTROM_PER_PS2 / masses[i];
                let got = sys.atoms[i].velocity[d];
                assert!(
                    (got - expected).abs() <= 1.0e-12 * expected.abs().max(1.0e-3),
                    "atom {i} component {d}: dv = {got}, expected {expected}"
                );
            }
        }
        // Newton's third law in momentum form: m_0 v_0 + m_1 v_1 = 0, which
        // only holds when each species' own mass is used.
        for d in 0..3 {
            let p = masses[0] * sys.atoms[0].velocity[d] + masses[1] * sys.atoms[1].velocity[d];
            assert!(p.abs() < 1.0e-10, "momentum component {d} = {p:e}");
        }
    }

    #[test]
    fn test_step_treats_massless_species_as_immobile() {
        let mut ghost = EamParams::copper();
        ghost.mass = 0.0;
        let mut sys = AlloySystem::new(
            vec![ghost, EamParams::copper()],
            SimBox::cubic(30.0),
            MixingRule::Geometric,
            0.001,
        );
        sys.add_atom([10.0, 10.0, 10.0], 0);
        sys.add_atom([11.2, 10.0, 10.0], 1);
        sys.compute_forces();
        sys.step();
        assert_eq!(sys.atoms[0].velocity, [0.0; 3]);
        assert_eq!(sys.atoms[0].position, [10.0, 10.0, 10.0]);
        assert!(sys.atoms[1].velocity[0].is_finite());
    }

    /// Total energy drift of a bound dimer over 2000 velocity-Verlet steps.
    fn dimer_energy_drift(species: Vec<EamParams>, pair: [usize; 2], r0: f64) -> (f64, f64) {
        let mut sys = AlloySystem::new(species, SimBox::cubic(30.0), MixingRule::Geometric, 0.001);
        sys.add_atom([10.0, 10.0, 10.0], pair[0]);
        sys.add_atom([10.0 + r0, 10.0, 10.0], pair[1]);
        sys.compute_forces();
        let e0 = sys.kinetic_energy() + sys.potential_energy();
        let mut worst = 0.0_f64;
        let (mut r_min, mut r_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for _ in 0..2000 {
            sys.step();
            let e = sys.kinetic_energy() + sys.potential_energy();
            worst = worst.max((e - e0).abs());
            let r = sys.atoms[1].position[0] - sys.atoms[0].position[0];
            r_min = r_min.min(r);
            r_max = r_max.max(r);
        }
        (worst, r_max - r_min)
    }

    #[test]
    fn test_copper_dimer_conserves_total_energy() {
        let req = copper_dimer_equilibrium();
        let (drift, amplitude) = dimer_energy_drift(vec![EamParams::copper()], [0, 0], req + 0.02);
        assert!(
            amplitude > 0.03,
            "the dimer must actually oscillate, range = {amplitude}"
        );
        assert!(drift < 1.0e-5, "max |E_total - E_0| = {drift:e} eV");
    }

    #[test]
    fn test_heteronuclear_dimer_conserves_total_energy() {
        let species = vec![EamParams::copper(), EamParams::nickel()];
        let mut lo = 0.8;
        let mut hi = 2.0;
        assert!(dimer_axial_force(&species, [0, 1], lo) > 0.0);
        assert!(dimer_axial_force(&species, [0, 1], hi) < 0.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if dimer_axial_force(&species, [0, 1], mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let req = 0.5 * (lo + hi);
        let (drift, amplitude) = dimer_energy_drift(species, [0, 1], req + 0.02);
        assert!(amplitude > 0.03, "range = {amplitude}");
        assert!(drift < 1.0e-5, "max |E_total - E_0| = {drift:e} eV");
    }
}
