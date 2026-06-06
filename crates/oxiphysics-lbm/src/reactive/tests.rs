// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for reactive transport models.

use super::*;

#[cfg(test)]
mod reactive_tests {
    use super::*;

    /// First-order reaction: after one step, C should equal C0 * exp(-k*dt).
    #[test]
    fn test_first_order_reaction() {
        let mut c = ConcentrationField::new(4, 4, 1e-3, 1.0);
        let k: f64 = 0.5;
        let dt: f64 = 0.1;
        first_order_reaction(&mut c, k, dt);
        let expected = (-k * dt).exp();
        for &v in &c.data {
            assert!(
                (v - expected).abs() < 1e-14,
                "first_order_reaction: got {v}, expected {expected}"
            );
        }
    }

    /// Bimolecular reaction: check explicit-Euler update.
    #[test]
    fn test_bimolecular_reaction() {
        let a0 = 2.0;
        let b0 = 3.0;
        let rate = 0.1;
        let dt = 0.05;
        let mut a = ConcentrationField::new(2, 2, 1e-3, a0);
        let mut b = ConcentrationField::new(2, 2, 1e-3, b0);
        bimolecular_reaction(&mut a, &mut b, rate, dt);
        let delta = rate * a0 * b0 * dt;
        let expected_a = a0 - delta;
        let expected_b = b0 - delta;
        for &v in &a.data {
            assert!(
                (v - expected_a).abs() < 1e-14,
                "bimolecular A: got {v}, expected {expected_a}"
            );
        }
        for &v in &b.data {
            assert!(
                (v - expected_b).abs() < 1e-14,
                "bimolecular B: got {v}, expected {expected_b}"
            );
        }
    }

    /// Diffusion conserves total mass (periodic BC).
    #[test]
    fn test_diffusion_conserves_mass() {
        let nx = 8;
        let ny = 8;
        let mut c = ConcentrationField::new(nx, ny, 0.01, 0.0);
        // Place a spike in the middle.
        c.set(4, 4, 1.0);
        let mass_before = c.total_mass();

        for _ in 0..20 {
            c.diffusion_step(0.1, 1.0);
        }

        let mass_after = c.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Diffusion did not conserve mass: before={mass_before}, after={mass_after}"
        );
    }

    /// Passive scalar equilibrium distributions sum to the concentration c.
    #[test]
    fn test_passive_scalar_equilibrium_sums_to_c() {
        let c = 1.5;
        let u = [0.05, -0.02];
        let geq = LbmPassiveScalar::equilibrium(c, u);
        let sum: f64 = geq.iter().sum();
        assert!(
            (sum - c).abs() < 1e-14,
            "Passive scalar equilibrium sum = {sum}, expected {c}"
        );
    }

    // -------------------------------------------------------------------
    // Species / ArrheniusRate / BimolecularRate tests
    // -------------------------------------------------------------------

    #[test]
    fn test_species_creation() {
        let s = Species::new("O2".to_string(), 2.1e-5, 32.0);
        assert_eq!(s.name, "O2");
        assert!((s.diffusivity - 2.1e-5).abs() < 1e-20);
        assert!((s.molar_mass - 32.0).abs() < 1e-14);
    }

    #[test]
    fn test_arrhenius_rate() {
        let a = ArrheniusRate::new(1e10, 80_000.0);
        let r300 = a.rate(300.0);
        let r600 = a.rate(600.0);
        assert!(r300 > 0.0);
        assert!(r600 > r300, "rate should increase with temperature");
    }

    #[test]
    fn test_bimolecular_rate() {
        let rxn = BimolecularRate::new(0.5);
        let r = rxn.rate(2.0, 3.0);
        assert!((r - 3.0).abs() < 1e-14, "rate = {r}, expected 3.0");
    }

    #[test]
    fn test_reactive_lattice_diffusion_conserves() {
        let nx = 8;
        let ny = 8;
        let mut rl = ReactiveLattice::new(nx, ny);
        rl.add_species(Species::new("A".into(), 0.01, 28.0));
        // Place a spike.
        rl.concentrations[0][4 * nx + 4] = 5.0;
        let mass_before = rl.total_concentration(0);
        for _ in 0..10 {
            rl.advance_diffusion(0.1, 1.0);
        }
        let mass_after = rl.total_concentration(0);
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Diffusion did not conserve: before={mass_before}, after={mass_after}"
        );
    }

    #[test]
    fn test_reactive_lattice_apply_reactions() {
        let nx = 4;
        let ny = 4;
        let mut rl = ReactiveLattice::new(nx, ny);
        rl.add_species(Species::new("A".into(), 0.01, 28.0));
        rl.add_species(Species::new("B".into(), 0.01, 32.0));
        // Uniform concentrations.
        for v in rl.concentrations[0].iter_mut() {
            *v = 2.0;
        }
        for v in rl.concentrations[1].iter_mut() {
            *v = 3.0;
        }
        // stoich: A consumed, B consumed, rate = 0.1
        rl.reactions.push(Reaction {
            stoichiometry: vec![-1.0, -1.0],
            rate: 0.1,
        });
        rl.apply_reactions(0.05);
        // Check concentrations decreased.
        for &v in &rl.concentrations[0] {
            assert!(v < 2.0, "A should have decreased");
        }
    }

    #[test]
    fn test_compute_reaction_source() {
        let conc = vec![2.0, 3.0];
        let stoich = vec![-1.0, -2.0];
        let rate = 0.5;
        let src = compute_reaction_source(&conc, &stoich, rate);
        // source_i = stoich_i * rate * prod(conc) = stoich_i * 0.5 * 6.0
        assert!((src[0] - (-3.0)).abs() < 1e-14);
        assert!((src[1] - (-6.0)).abs() < 1e-14);
    }

    // -------------------------------------------------------------------
    // Bimolecular with product tests
    // -------------------------------------------------------------------

    #[test]
    fn test_bimolecular_with_product_conservation() {
        let mut a = ConcentrationField::new(4, 4, 1e-3, 2.0);
        let mut b = ConcentrationField::new(4, 4, 1e-3, 3.0);
        let mut c = ConcentrationField::new(4, 4, 1e-3, 0.0);
        let total_before: f64 = a.total_mass() + b.total_mass() + c.total_mass();
        bimolecular_reaction_with_product(&mut a, &mut b, &mut c, 0.1, 0.05);
        let total_after: f64 = a.total_mass() + b.total_mass() + c.total_mass();
        // A+B consumed = C produced: per cell, A+C is conserved
        let ac_before_cell: f64 = 2.0 + 0.0; // initial A + C per cell
        let ac_after_cell = a.data[0] + c.data[0];
        assert!(
            (ac_before_cell - ac_after_cell).abs() < 1e-12,
            "A+C should be conserved per cell"
        );
        // Verify product formed
        assert!(c.data[0] > 0.0, "Product C should have formed");
        // Ignore total_before/total_after to avoid unused warning
        let _ = (total_before, total_after);
    }

    #[test]
    fn test_bimolecular_product_amount() {
        let mut a = ConcentrationField::new(2, 2, 1e-3, 1.0);
        let mut b = ConcentrationField::new(2, 2, 1e-3, 1.0);
        let mut c = ConcentrationField::new(2, 2, 1e-3, 0.0);
        bimolecular_reaction_with_product(&mut a, &mut b, &mut c, 0.5, 0.1);
        let expected_delta = 0.5 * 1.0 * 1.0 * 0.1;
        assert!(
            (c.data[0] - expected_delta).abs() < 1e-14,
            "Product = {}, expected {}",
            c.data[0],
            expected_delta
        );
    }

    // -------------------------------------------------------------------
    // Reversible reaction tests
    // -------------------------------------------------------------------

    #[test]
    fn test_reversible_reaction_equilibrium() {
        // At equilibrium: kf*A = kr*B => A/B = kr/kf
        let kf = 0.2;
        let kr = 0.1;
        let mut a = ConcentrationField::new(2, 2, 1e-3, 1.0);
        let mut b = ConcentrationField::new(2, 2, 1e-3, 0.0);
        // Run many steps to approach equilibrium
        for _ in 0..100000 {
            reversible_reaction(&mut a, &mut b, kf, kr, 0.001);
        }
        // At eq: kf*A = kr*B => A/B = kr/kf = 0.5
        // Total = A + B = 1.0 (conserved)
        let ratio = a.data[0] / b.data[0];
        assert!(
            (ratio - 0.5).abs() < 0.1,
            "Equilibrium ratio should be ~0.5, got {ratio}"
        );
    }

    #[test]
    fn test_reversible_reaction_conserves_total() {
        let mut a = ConcentrationField::new(4, 4, 1e-3, 2.0);
        let mut b = ConcentrationField::new(4, 4, 1e-3, 1.0);
        let total_before = a.total_mass() + b.total_mass();
        for _ in 0..100 {
            reversible_reaction(&mut a, &mut b, 0.1, 0.05, 0.01);
        }
        let total_after = a.total_mass() + b.total_mass();
        assert!(
            (total_before - total_after).abs() < 1e-10,
            "Reversible reaction should conserve A+B"
        );
    }

    // -------------------------------------------------------------------
    // MultiStepReaction tests
    // -------------------------------------------------------------------

    #[test]
    fn test_multi_step_reaction_chain() {
        // A -> B -> C with rates k1 and k2
        let msr = MultiStepReaction::new(vec![0, 1, 2], vec![0.1, 0.05]);
        assert_eq!(msr.num_steps(), 2);

        let n = 4;
        let mut concentrations = vec![vec![1.0; n], vec![0.0; n], vec![0.0; n]];
        msr.apply(&mut concentrations, 0.1);

        // A should decrease, B and C should increase
        assert!(concentrations[0][0] < 1.0, "A should decrease");
        assert!(concentrations[1][0] > 0.0, "B should increase");
        // C might be 0 if B was 0 at start (second step: B->C with B=0 initially)
        // After first step: B has some value, then second step converts some B->C
        // Actually in this apply, both steps run sequentially in same call
        // Step 0: A -> B: B gets delta = 0.1 * 1.0 * 0.1 = 0.01
        // Step 1: B -> C: C gets delta = 0.05 * 0.01 * 0.1 = 0.00005
        assert!(concentrations[2][0] >= 0.0, "C should be non-negative");
    }

    #[test]
    fn test_multi_step_conservation() {
        let msr = MultiStepReaction::new(vec![0, 1, 2], vec![0.1, 0.05]);
        let n = 4;
        let mut concentrations = vec![vec![1.0; n], vec![0.5; n], vec![0.0; n]];
        let total_before: f64 = concentrations.iter().flatten().sum();
        msr.apply(&mut concentrations, 0.1);
        let total_after: f64 = concentrations.iter().flatten().sum();
        assert!(
            (total_before - total_after).abs() < 1e-10,
            "Multi-step reaction should conserve total"
        );
    }

    // -------------------------------------------------------------------
    // CatalyticSurface tests
    // -------------------------------------------------------------------

    #[test]
    fn test_catalytic_surface_effective_rate() {
        let nx = 4;
        let ny = 4;
        let mut cat = CatalyticSurface::new(nx, ny, 10.0);
        cat.set_catalytic_row(nx, 0);
        assert_eq!(cat.num_catalytic_cells(), nx);

        // Catalytic cell
        assert!(
            (cat.effective_rate(0, 1.0) - 10.0).abs() < 1e-14,
            "Catalytic rate should be 10x"
        );
        // Non-catalytic cell
        assert!(
            (cat.effective_rate(nx, 1.0) - 1.0).abs() < 1e-14,
            "Non-catalytic rate should be 1x"
        );
    }

    #[test]
    fn test_catalytic_surface_coverage_update() {
        let nx = 2;
        let ny = 2;
        let mut cat = CatalyticSurface::new(nx, ny, 5.0).with_adsorption(1.0, 0.1);
        cat.mask[0] = true;
        let conc = vec![1.0, 0.0, 0.0, 0.0];
        cat.update_coverage(&conc, 0.1);
        // Cell 0 is catalytic with c=1.0:
        // d_theta = 1.0 * 1.0 * (1-0) - 0.1 * 0 = 1.0
        // theta = 0 + 1.0 * 0.1 = 0.1
        assert!(
            (cat.coverage[0] - 0.1).abs() < 1e-12,
            "Coverage should be 0.1, got {}",
            cat.coverage[0]
        );
        // Cell 1 is not catalytic, coverage stays 0
        assert!(
            cat.coverage[1].abs() < 1e-14,
            "Non-catalytic coverage should be 0"
        );
    }

    #[test]
    fn test_catalytic_coverage_saturates() {
        let mut cat = CatalyticSurface::new(1, 1, 2.0).with_adsorption(100.0, 0.0);
        cat.mask[0] = true;
        let conc = vec![1.0];
        // Large dt to overshoot
        cat.update_coverage(&conc, 1.0);
        assert!(
            (cat.coverage[0] - 1.0).abs() < 1e-14,
            "Coverage should be clamped to 1.0"
        );
    }

    // -------------------------------------------------------------------
    // SpeciesBoundaryCondition tests
    // -------------------------------------------------------------------

    #[test]
    fn test_fixed_concentration_bc() {
        let mut field = ConcentrationField::new(8, 8, 0.01, 0.0);
        let bc = SpeciesBoundaryCondition::new(
            BoundaryEdge::Left,
            SpeciesBcType::FixedConcentration(1.0),
        );
        bc.apply(&mut field);
        for j in 0..8 {
            assert!(
                (field.get(0, j) - 1.0).abs() < 1e-14,
                "Left BC should set concentration to 1.0"
            );
        }
    }

    #[test]
    fn test_zero_flux_bc_right() {
        let mut field = ConcentrationField::new(8, 8, 0.01, 0.0);
        // Set interior to some value
        for j in 0..8 {
            field.set(6, j, 2.5);
        }
        let bc = SpeciesBoundaryCondition::new(BoundaryEdge::Right, SpeciesBcType::ZeroFlux);
        bc.apply(&mut field);
        for j in 0..8 {
            assert!(
                (field.get(7, j) - 2.5).abs() < 1e-14,
                "Right zero-flux BC should copy from interior"
            );
        }
    }

    #[test]
    fn test_convective_outflow_bc() {
        let mut field = ConcentrationField::new(8, 8, 0.01, 0.0);
        for i in 0..8 {
            field.set(i, 6, 3.0);
        }
        let bc = SpeciesBoundaryCondition::new(BoundaryEdge::Top, SpeciesBcType::ConvectiveOutflow);
        bc.apply(&mut field);
        for i in 0..8 {
            assert!(
                (field.get(i, 7) - 3.0).abs() < 1e-14,
                "Top outflow BC should extrapolate"
            );
        }
    }

    // -------------------------------------------------------------------
    // ReactionFrontTracker tests
    // -------------------------------------------------------------------

    #[test]
    fn test_reaction_front_tracker_finds_front() {
        let nx = 16;
        let ny = 4;
        let mut field = ConcentrationField::new(nx, ny, 0.01, 0.0);
        // Create a step profile: left half = 1.0, right half = 0.0
        for j in 0..ny {
            for i in 0..nx / 2 {
                field.set(i, j, 1.0);
            }
        }
        let mut tracker = ReactionFrontTracker::new(0.5);
        tracker.track(&field, 0.0);
        assert_eq!(tracker.num_records(), 1);
        // Front should be near x = 7.5 (between index 7 and 8)
        let (_, x) = tracker.history[0];
        assert!(
            (x - 7.5).abs() < 1.0,
            "Front position should be near 7.5, got {x}"
        );
    }

    #[test]
    fn test_reaction_front_speed() {
        let mut tracker = ReactionFrontTracker::new(0.5);
        tracker.history.push((0.0, 5.0));
        tracker.history.push((1.0, 8.0));
        let speed = tracker.front_speed();
        assert!(
            (speed - 3.0).abs() < 1e-14,
            "Front speed should be 3.0, got {speed}"
        );
    }

    #[test]
    fn test_reaction_front_speed_insufficient_data() {
        let tracker = ReactionFrontTracker::new(0.5);
        assert!(
            tracker.front_speed().abs() < 1e-14,
            "Speed with no data should be 0"
        );
    }

    // -------------------------------------------------------------------
    // HeatRelease tests
    // -------------------------------------------------------------------

    #[test]
    fn test_heat_release_exothermic() {
        let nx = 4;
        let ny = 4;
        let mut hr = HeatRelease::new(nx, ny, 10.0, 300.0);
        let conc = vec![1.0; nx * ny];
        let t_before = hr.mean_temperature();
        hr.apply_first_order(&conc, 0.1, 0.01);
        let t_after = hr.mean_temperature();
        assert!(
            t_after > t_before,
            "Exothermic reaction should raise temperature"
        );
    }

    #[test]
    fn test_heat_release_temperature_increase() {
        let nx = 2;
        let ny = 2;
        let mut hr = HeatRelease::new(nx, ny, 5.0, 300.0);
        let conc = vec![2.0; nx * ny];
        hr.apply_first_order(&conc, 0.1, 0.1);
        // dT = 5.0 * 0.1 * 2.0 * 0.1 = 0.1
        let expected = 300.1;
        assert!(
            (hr.temperature[0] - expected).abs() < 1e-12,
            "T = {}, expected {}",
            hr.temperature[0],
            expected
        );
    }

    #[test]
    fn test_heat_release_diffusion_conserves_energy() {
        let nx = 8;
        let ny = 8;
        let mut hr = HeatRelease::new(nx, ny, 1.0, 300.0).with_thermal_diffusivity(0.01);
        // Hot spot
        hr.temperature[4 * nx + 4] = 400.0;
        let total_before: f64 = hr.temperature.iter().sum();
        for _ in 0..10 {
            hr.diffuse_temperature(nx, ny, 0.1, 1.0);
        }
        let total_after: f64 = hr.temperature.iter().sum();
        assert!(
            (total_before - total_after).abs() < 1e-8,
            "Temperature diffusion should conserve total energy"
        );
    }

    #[test]
    fn test_heat_release_max_temperature() {
        let mut hr = HeatRelease::new(4, 4, 1.0, 300.0);
        hr.temperature[5] = 500.0;
        assert!(
            (hr.max_temperature() - 500.0).abs() < 1e-14,
            "Max temperature should be 500"
        );
    }

    // -------------------------------------------------------------------
    // Arrhenius ratio test
    // -------------------------------------------------------------------

    #[test]
    fn test_arrhenius_rate_ratio() {
        let a = ArrheniusRate::new(1e10, 80_000.0);
        let ratio = a.rate_ratio(300.0, 600.0);
        // ratio = exp((Ea/R) * (1/300 - 1/600)) = exp(Ea/(R*600))
        // Should be > 1 since higher temp gives higher rate
        assert!(
            ratio > 1.0,
            "Rate ratio (600K/300K) should be > 1, got {ratio}"
        );
        // Also verify it matches rate(600)/rate(300)
        let direct = a.rate(600.0) / a.rate(300.0);
        assert!(
            (ratio - direct).abs() / direct < 1e-6,
            "Rate ratio should match direct computation"
        );
    }

    // -------------------------------------------------------------------
    // ConcentrationField helper tests
    // -------------------------------------------------------------------

    #[test]
    fn test_concentration_field_max_min() {
        let mut c = ConcentrationField::new(4, 4, 0.01, 1.0);
        c.set(0, 0, 5.0);
        c.set(1, 1, -0.5);
        assert!((c.max_concentration() - 5.0).abs() < 1e-14);
        assert!((c.min_concentration() - (-0.5)).abs() < 1e-14);
    }

    #[test]
    fn test_concentration_field_clamp_min() {
        let mut c = ConcentrationField::new(4, 4, 0.01, 1.0);
        c.set(0, 0, -0.5);
        c.clamp_min(0.0);
        assert!(c.min_concentration() >= 0.0, "All values should be >= 0");
        assert!(
            (c.get(1, 1) - 1.0).abs() < 1e-14,
            "Non-negative values unchanged"
        );
    }

    // -------------------------------------------------------------------
    // LbmPassiveScalar additional tests
    // -------------------------------------------------------------------

    #[test]
    fn test_passive_scalar_initialize_concentration() {
        let nx = 4;
        let ny = 4;
        let mut ps = LbmPassiveScalar::new(nx, ny, 0.1);
        let conc: Vec<f64> = (0..nx * ny).map(|i| i as f64 * 0.1).collect();
        ps.initialize_concentration(&conc);
        for (k, &expected_c) in conc.iter().enumerate() {
            let c = ps.concentration(k);
            assert!(
                (c - expected_c).abs() < 1e-12,
                "Initialized concentration mismatch at {k}"
            );
        }
    }

    #[test]
    fn test_passive_scalar_total_conserved() {
        let nx = 8;
        let ny = 8;
        let mut ps = LbmPassiveScalar::new(nx, ny, 0.05);
        let mut conc = vec![0.0; nx * ny];
        conc[4 * nx + 4] = 1.0;
        ps.initialize_concentration(&conc);
        let total_before = ps.total_concentration();

        let zero_vel = vec![[0.0f64; 2]; nx * ny];
        for _ in 0..5 {
            ps.collide(&zero_vel);
            ps.stream();
        }
        let total_after = ps.total_concentration();
        assert!(
            (total_before - total_after).abs() < 1e-10,
            "LBM passive scalar should conserve total concentration"
        );
    }

    #[test]
    fn test_passive_scalar_add_source() {
        let nx = 4;
        let ny = 4;
        let mut ps = LbmPassiveScalar::new(nx, ny, 0.1);
        let conc = vec![1.0; nx * ny];
        ps.initialize_concentration(&conc);
        let before = ps.total_concentration();
        let sources = vec![0.5; nx * ny];
        ps.add_source(&sources);
        let after = ps.total_concentration();
        assert!(
            after > before,
            "Total concentration should increase after adding source"
        );
    }

    // -------------------------------------------------------------------
    // Damkoehler and Peclet number tests
    // -------------------------------------------------------------------

    #[test]
    fn test_damkoehler_number() {
        let da = damkoehler_number(0.1, 10.0, 0.01);
        // Da = 0.1 * 100 / 0.01 = 1000
        assert!((da - 1000.0).abs() < 1e-10, "Da should be 1000, got {da}");
    }

    #[test]
    fn test_peclet_number() {
        let pe = peclet_number(1.0, 10.0, 0.1);
        assert!((pe - 100.0).abs() < 1e-10, "Pe should be 100, got {pe}");
    }

    // -------------------------------------------------------------------
    // ReactiveLattice additional tests
    // -------------------------------------------------------------------

    #[test]
    fn test_reactive_lattice_set_uniform() {
        let mut rl = ReactiveLattice::new(4, 4);
        rl.add_species(Species::new("X".into(), 0.01, 10.0));
        rl.set_uniform(0, 5.0);
        assert!(
            (rl.total_concentration(0) - 5.0 * 16.0).abs() < 1e-10,
            "Total should be 5*16=80"
        );
    }

    #[test]
    fn test_reactive_lattice_num_species() {
        let mut rl = ReactiveLattice::new(4, 4);
        assert_eq!(rl.num_species(), 0);
        rl.add_species(Species::new("A".into(), 0.01, 28.0));
        rl.add_species(Species::new("B".into(), 0.02, 32.0));
        assert_eq!(rl.num_species(), 2);
    }
}

mod reactive_additions_tests {
    use super::*;

    // Arrhenius rate: high T gives higher rate than low T.
    #[test]
    fn test_arrhenius_high_vs_low_temperature() {
        let a = 1e10_f64;
        let ea = 50_000.0_f64; // J/mol
        let r = R_GAS;
        let k_high = arrhenius_rate(a, ea, 1000.0, r);
        let k_low = arrhenius_rate(a, ea, 300.0, r);
        assert!(
            k_high > k_low,
            "High-T rate {k_high} should exceed low-T rate {k_low}"
        );
    }

    // ReactionRate::Arrhenius evaluate matches free function.
    #[test]
    fn test_reaction_rate_arrhenius_evaluate() {
        let rr = ReactionRate::Arrhenius {
            a: 1e8,
            ea: 40_000.0,
            t: 500.0,
        };
        let expected = arrhenius_rate(1e8, 40_000.0, 500.0, R_GAS);
        assert!(
            (rr.evaluate() - expected).abs() < 1e-10 * expected.abs(),
            "ReactionRate evaluate mismatch"
        );
    }

    // Zero reaction: react_species with rate 0 should not change concentrations.
    #[test]
    fn test_zero_reaction_no_change() {
        let mut rlbm = ReactiveLbm::new(4, 4, 1, 1.0);
        let n = 16;
        for k in 0..n {
            rlbm.species[0][k] = 1.0;
        }
        rlbm.add_reaction(0, ReactionRate::Custom(0.0), -1.0);
        rlbm.react_species(0.1);
        for k in 0..n {
            assert!(
                (rlbm.species[0][k] - 1.0).abs() < 1e-14,
                "Concentration changed with zero rate at {k}"
            );
        }
    }

    // Species conservation: two coupled species (A -> B) should conserve total.
    #[test]
    fn test_species_conservation_coupled_reaction() {
        let mut rlbm = ReactiveLbm::new(4, 4, 2, 1.0);
        let n = 16;
        for k in 0..n {
            rlbm.species[0][k] = 2.0; // species A
            rlbm.species[1][k] = 0.0; // species B
        }
        let rate = 0.5;
        rlbm.add_reaction(0, ReactionRate::Custom(rate), -1.0);
        rlbm.add_reaction(1, ReactionRate::Custom(rate), 1.0);
        let dt = 0.1;
        let total_before: f64 = rlbm.species[0].iter().chain(rlbm.species[1].iter()).sum();
        rlbm.react_species(dt);
        let total_after: f64 = rlbm.species[0].iter().chain(rlbm.species[1].iter()).sum();
        assert!(
            (total_before - total_after).abs() < 1e-12,
            "Total species not conserved: before={total_before}, after={total_after}"
        );
    }

    // compute_mixture_density: single species with uniform c=1 → rho_mix=base_rho everywhere.
    #[test]
    fn test_compute_mixture_density_uniform() {
        let mut rlbm = ReactiveLbm::new(3, 3, 1, 1.5);
        for k in 0..9 {
            rlbm.species[0][k] = 1.0;
        }
        let rho = rlbm.compute_mixture_density();
        for (k, &r) in rho.iter().enumerate() {
            assert!((r - 1.5).abs() < 1e-14, "rho_mix[{k}] = {r}, expected 1.5");
        }
    }

    // SpeciesTransport: total mass check.
    #[test]
    fn test_species_transport_total_mass() {
        let diff = vec![0.01, 0.02];
        let mut st = SpeciesTransport::new(2, 4, diff);
        for k in 0..4 {
            st.set(0, k, 1.0);
            st.set(1, k, 2.0);
        }
        assert!((st.total_mass(0) - 4.0).abs() < 1e-14);
        assert!((st.total_mass(1) - 8.0).abs() < 1e-14);
    }
}

mod species_combustion_tests {
    use super::*;
    use crate::reactive::combustion::combustion_indices::*;

    // -----------------------------------------------------------------------
    // SpeciesGrid: equilibrium sums to concentration
    // -----------------------------------------------------------------------
    #[test]
    fn test_species_grid_equilibrium_sum() {
        let c = 2.5_f64;
        let geq = SpeciesGrid::equilibrium(c, 0.05, -0.02);
        let sum: f64 = geq.iter().sum();
        assert!(
            (sum - c).abs() < 1e-13,
            "SpeciesGrid equilibrium sum = {sum}, expected {c}"
        );
    }

    // -----------------------------------------------------------------------
    // SpeciesGrid: total concentration conserved under collide+stream
    // -----------------------------------------------------------------------
    #[test]
    fn test_species_grid_conservation() {
        let nx = 8;
        let ny = 8;
        let mut sg = SpeciesGrid::new(nx, ny, vec![0.05, 0.1]);
        let n = nx * ny;
        // Initialize species 0 to a spike.
        sg.set_concentration(0, 4 * nx + 4, 1.0);
        let total_before = sg.total_concentration(0);
        let zero_vel = vec![[0.0_f64; 2]; n];
        for _ in 0..5 {
            sg.collide(&zero_vel);
            sg.stream();
        }
        let total_after = sg.total_concentration(0);
        assert!(
            (total_before - total_after).abs() < 1e-10,
            "SpeciesGrid: total concentration not conserved: before={total_before}, after={total_after}"
        );
    }

    // -----------------------------------------------------------------------
    // SpeciesGrid: apply_source increases total concentration
    // -----------------------------------------------------------------------
    #[test]
    fn test_species_grid_apply_source() {
        let nx = 4;
        let ny = 4;
        let mut sg = SpeciesGrid::new(nx, ny, vec![0.01]);
        let n = nx * ny;
        let before = sg.total_concentration(0);
        let sources = vec![0.1_f64; n];
        sg.apply_source(0, &sources);
        let after = sg.total_concentration(0);
        assert!(
            after > before,
            "apply_source should increase total concentration"
        );
    }

    // -----------------------------------------------------------------------
    // Arrhenius: k(T) = A * exp(-Ea / (R*T)) at known T
    // -----------------------------------------------------------------------
    #[test]
    fn test_arrhenius_known_value() {
        // k(T) = 1 * exp(-8314 / (8.314 * 1000)) = exp(-1) ≈ 0.3679
        let a = 1.0_f64;
        let ea = 8314.0_f64; // Ea = R * 1000 K
        let t = 1000.0_f64;
        let k = arrhenius_rate(a, ea, t, R_GAS);
        let expected = (-1.0_f64).exp();
        assert!(
            (k - expected).abs() < 1e-10,
            "Arrhenius k at 1000K = {k}, expected {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // Arrhenius: rate increases with temperature
    // -----------------------------------------------------------------------
    #[test]
    fn test_arrhenius_rate_increases_with_temperature() {
        let rate = ArrheniusRate::new(1e12, 100_000.0);
        let r_low = rate.rate(500.0);
        let r_high = rate.rate(1500.0);
        assert!(
            r_high > r_low,
            "Arrhenius: rate at 1500K ({r_high}) should exceed rate at 500K ({r_low})"
        );
    }

    // -----------------------------------------------------------------------
    // CombustionReactor: CH4 consumed, CO2/H2O produced
    // -----------------------------------------------------------------------
    #[test]
    fn test_combustion_reactor_step_consumes_ch4() {
        use combustion_indices::*;
        let mut reactor = CombustionReactor::new(4, 1500.0); // high T to ignite
        for k in 0..4 {
            reactor.set(CH4, k, 1.0);
            reactor.set(O2, k, 2.0);
        }
        reactor.step(1e-4);
        for k in 0..4 {
            assert!(
                reactor.get(CH4, k) < 1.0,
                "CH4 should be consumed at cell {k}"
            );
            assert!(
                reactor.get(CO2, k) > 0.0,
                "CO2 should be produced at cell {k}"
            );
            assert!(
                reactor.get(H2O, k) > 0.0,
                "H2O should be produced at cell {k}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // CombustionReactor: stoichiometric conservation (approx 1 CH4 + 2 O2 → 1 CO2 + 2 H2O)
    // -----------------------------------------------------------------------
    #[test]
    fn test_combustion_stoich_ratio() {
        let mut reactor = CombustionReactor::new(1, 2000.0);
        reactor.set(CH4, 0, 1.0);
        reactor.set(O2, 0, 2.0);
        reactor.step(1e-6);
        let d_ch4: f64 = 1.0 - reactor.get(CH4, 0);
        let d_o2: f64 = 2.0 - reactor.get(O2, 0);
        let d_co2 = reactor.get(CO2, 0);
        let d_h2o = reactor.get(H2O, 0);
        // All deltas should be positive and O2 consumed ≈ 2× CH4 consumed.
        assert!(d_ch4 > 0.0, "CH4 delta = {d_ch4}");
        assert!(
            (d_o2 / d_ch4 - 2.0).abs() < 1e-8,
            "O2/CH4 consumption ratio should be 2: got {}",
            d_o2 / d_ch4
        );
        assert!(
            (d_co2 / d_ch4 - 1.0).abs() < 1e-8,
            "CO2/CH4 production ratio should be 1: got {}",
            d_co2 / d_ch4
        );
        assert!(
            (d_h2o / d_ch4 - 2.0).abs() < 1e-8,
            "H2O/CH4 production ratio should be 2: got {}",
            d_h2o / d_ch4
        );
    }

    // -----------------------------------------------------------------------
    // CombustionReactor: exothermic reaction raises temperature
    // -----------------------------------------------------------------------
    #[test]
    fn test_combustion_exothermic_heating() {
        let t0 = 2000.0;
        let mut reactor = CombustionReactor::new(1, t0);
        reactor.set(CH4, 0, 1.0);
        reactor.set(O2, 0, 2.0);
        reactor.step(1e-4);
        // ΔH < 0 → exothermic → temperature should rise.
        assert!(
            reactor.temperature[0] > t0,
            "Temperature should increase due to exothermic combustion"
        );
    }

    // -----------------------------------------------------------------------
    // CombustionReactor: mass conservation check
    // (total concentration of reactants + products is conserved by stoichiometry)
    // CH4 consumed contributes: -1 CH4 - 2 O2 + 1 CO2 + 2 H2O → net = 0 atoms?
    // Actually C is conserved: C in CH4 = C in CO2.
    // -----------------------------------------------------------------------
    #[test]
    fn test_combustion_carbon_conservation() {
        let mut reactor = CombustionReactor::new(1, 2000.0);
        reactor.set(CH4, 0, 1.0);
        reactor.set(O2, 0, 4.0); // excess O2
        for _ in 0..100 {
            reactor.step(1e-5);
        }
        // C in = CH4_initial = 1.0; C out = CH4_remaining + CO2
        let c_in = 1.0_f64;
        let c_out = reactor.get(CH4, 0) + reactor.get(CO2, 0);
        assert!(
            (c_in - c_out).abs() < 1e-8,
            "Carbon not conserved: in={c_in}, out={c_out}"
        );
    }

    // -----------------------------------------------------------------------
    // species_source_terms: shape and values
    // -----------------------------------------------------------------------
    #[test]
    fn test_species_source_terms() {
        let stoich = vec![-1.0, -2.0, 1.0, 2.0];
        let rates = vec![0.5, 1.0, 0.25];
        let src = species_source_terms(&stoich, &rates);
        assert_eq!(src.len(), 4);
        assert_eq!(src[0].len(), 3);
        // src[0][1] = -1.0 * 1.0 = -1.0
        assert!(
            (src[0][1] - (-1.0)).abs() < 1e-14,
            "src[CH4][1] = {}",
            src[0][1]
        );
        // src[3][0] = 2.0 * 0.5 = 1.0
        assert!(
            (src[3][0] - 1.0).abs() < 1e-14,
            "src[H2O][0] = {}",
            src[3][0]
        );
    }

    // -----------------------------------------------------------------------
    // apply_heat_source: temperature increases by correct amount
    // -----------------------------------------------------------------------
    #[test]
    fn test_apply_heat_source() {
        let mut temp: Vec<f64> = vec![300.0, 300.0, 300.0];
        let rates: Vec<f64> = vec![1.0, 2.0, 0.5];
        let delta_h = 10.0;
        let dt = 0.1;
        apply_heat_source(&mut temp, &rates, delta_h, dt);
        // dT[0] = 10 * 1.0 * 0.1 = 1.0
        assert!((temp[0] - 301.0).abs() < 1e-12, "temp[0] = {}", temp[0]);
        // dT[1] = 10 * 2.0 * 0.1 = 2.0
        assert!((temp[1] - 302.0).abs() < 1e-12, "temp[1] = {}", temp[1]);
    }

    // -----------------------------------------------------------------------
    // FlameFrontTracker: records a position
    // -----------------------------------------------------------------------
    #[test]
    fn test_flame_front_tracker_detects_gradient() {
        let nx = 16_usize;
        let ny = 4_usize;
        // Create a step temperature profile: left=300K, right=1200K with a steep gradient at i=7-8.
        let mut temperature = vec![300.0_f64; nx * ny];
        for j in 0..ny {
            for i in (nx / 2)..nx {
                temperature[j * nx + i] = 1200.0;
            }
        }
        let mut tracker = FlameFrontTracker::new(100.0);
        tracker.track(&temperature, nx, ny, 0.0);
        assert_eq!(tracker.num_records(), 1);
        let (_, x) = tracker.history[0];
        assert!(!x.is_nan(), "Flame front position should not be NaN");
        assert!(
            x > 5.0 && x < 11.0,
            "Flame front should be near x=7-8, got x={x}"
        );
    }

    // -----------------------------------------------------------------------
    // FlameFrontTracker: flame speed from two records
    // -----------------------------------------------------------------------
    #[test]
    fn test_flame_front_speed() {
        let mut tracker = FlameFrontTracker::new(50.0);
        tracker.history.push((0.0, 4.0));
        tracker.history.push((2.0, 10.0));
        let speed = tracker.flame_speed();
        assert!(
            (speed - 3.0).abs() < 1e-12,
            "Flame speed should be 3.0, got {speed}"
        );
    }

    // -----------------------------------------------------------------------
    // FlameFrontTracker: no records → zero speed
    // -----------------------------------------------------------------------
    #[test]
    fn test_flame_front_speed_no_records() {
        let tracker = FlameFrontTracker::new(50.0);
        assert_eq!(tracker.flame_speed(), 0.0, "No records → zero speed");
    }
}
