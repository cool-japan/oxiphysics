// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended tests for reactive transport models.

use super::*;

mod reactive_extended_tests {
    use super::*;

    // ─── ZeldovichMechanism tests ───

    #[test]
    fn test_zeldovich_rates_increase_with_temperature() {
        let mech = ZeldovichMechanism::default_air();
        let k1_low = mech.k1f(1500.0);
        let k1_high = mech.k1f(2500.0);
        assert!(
            k1_high > k1_low,
            "k1f should increase with T: {k1_high} vs {k1_low}"
        );
    }

    #[test]
    fn test_zeldovich_no_production_positive() {
        let mech = ZeldovichMechanism::default_air();
        // At high temperature with plenty of O atoms
        let rate = mech.no_production_rate(2500.0, 1e-3, 0.8, 0.2, 0.0);
        assert!(rate > 0.0, "NO production rate should be positive: {rate}");
    }

    #[test]
    fn test_zeldovich_equilibrium_no_positive() {
        let mech = ZeldovichMechanism::default_air();
        let c_no_eq = mech.equilibrium_no(2500.0, 0.8, 0.2);
        assert!(
            c_no_eq >= 0.0,
            "Equilibrium NO should be non-negative: {c_no_eq}"
        );
    }

    #[test]
    fn test_zeldovich_equilibrium_no_increases_with_temperature() {
        let mech = ZeldovichMechanism::default_air();
        let c_low = mech.equilibrium_no(1500.0, 0.8, 0.2);
        let c_high = mech.equilibrium_no(2500.0, 0.8, 0.2);
        // Zeldovich: equilibrium NO increases with T
        assert!(
            c_high > c_low,
            "Equilibrium NO should increase with T: {c_high} vs {c_low}"
        );
    }

    // ─── ChainBranching tests ───

    #[test]
    fn test_chain_branching_factor_negative_at_low_t() {
        // High termination rate, high activation energy for branching → sub-critical at low T
        let cb = ChainBranching::new(1e10, 150_000.0, 1e12, 150_000.0, 1e10);
        let phi = cb.net_branching_factor(400.0, 0.1);
        // At low T with very high k_term, phi < 0
        assert!(phi < 0.0, "Low-T chain should be sub-critical: {phi}");
    }

    #[test]
    fn test_chain_branching_factor_positive_at_high_t() {
        let cb = ChainBranching::new(1e10, 50_000.0, 1e15, 50_000.0, 1.0);
        let phi = cb.net_branching_factor(2000.0, 1.0);
        // At high T with large k_b, phi > 0
        assert!(phi > 0.0, "High-T chain should be super-critical: {phi}");
    }

    #[test]
    fn test_chain_branching_explosion_temperature() {
        // k_branch large with high Ea → phi negative at 300K, positive at high T
        let cb = ChainBranching::new(1e10, 100_000.0, 1e15, 150_000.0, 1000.0);
        let t_exp = cb.explosion_temperature(0.5, 300.0, 3000.0);
        assert!(
            t_exp > 300.0 && t_exp < 3000.0,
            "Explosion T out of range: {t_exp}"
        );
    }

    #[test]
    fn test_chain_branching_quasi_steady_oh_some_below_limit() {
        // Below explosion limit: QSS [OH] should be finite
        let cb = ChainBranching::new(1e6, 80_000.0, 1e8, 40_000.0, 1e5);
        let oh = cb.quasi_steady_oh(400.0, 0.1, 0.1);
        assert!(
            oh.is_some(),
            "Should have QSS solution below explosion limit"
        );
        let c_oh = oh.unwrap();
        assert!(c_oh >= 0.0, "QSS [OH] should be non-negative: {c_oh}");
    }

    // ─── IgnitionDelayModel tests ───

    #[test]
    fn test_ignition_delay_decreases_with_temperature() {
        let model = IgnitionDelayModel::new(1e-12, -0.1, -0.5, 125_000.0);
        let tau_low = model.delay_time(800.0, 1.0, 1.0);
        let tau_high = model.delay_time(1200.0, 1.0, 1.0);
        assert!(
            tau_high < tau_low,
            "Ignition delay should decrease with T: {tau_high} vs {tau_low}"
        );
    }

    #[test]
    fn test_ignition_delay_effective_activation_energy() {
        let model = IgnitionDelayModel::new(1e-12, -0.1, -0.5, 125_000.0);
        let t1 = 1000.0;
        let t2 = 1200.0;
        let tau1 = model.delay_time(t1, 1.0, 1.0);
        let tau2 = model.delay_time(t2, 1.0, 1.0);
        let ea_eff = model.effective_activation_energy(t1, t2, tau1, tau2);
        // Effective Ea should be close to 125 kJ/mol
        assert!(
            (ea_eff - 125_000.0).abs() / 125_000.0 < 0.05,
            "Effective Ea = {ea_eff}, expected ~125000"
        );
    }

    #[test]
    fn test_ignition_delay_damkohler_number_positive() {
        let model = IgnitionDelayModel::new(1e-12, -0.1, -0.5, 125_000.0);
        let da = model.damkohler_number(0.01, 1000.0, 1.0, 1.0);
        assert!(da > 0.0, "Da_ign should be positive: {da}");
    }

    // ─── LaminarFlameSpeed tests ───

    #[test]
    fn test_laminar_flame_speed_increases_with_temperature() {
        let model = LaminarFlameSpeed::methane_air();
        let s_low = model.speed(298.0, 101325.0, 0.0);
        let s_high = model.speed(500.0, 101325.0, 0.0);
        assert!(
            s_high > s_low,
            "Flame speed should increase with unburnt T: {s_high} vs {s_low}"
        );
    }

    #[test]
    fn test_laminar_flame_speed_decreases_with_diluent() {
        let model = LaminarFlameSpeed::methane_air();
        let s_clean = model.speed(400.0, 101325.0, 0.0);
        let s_diluted = model.speed(400.0, 101325.0, 0.1);
        assert!(
            s_diluted < s_clean,
            "Diluent should reduce flame speed: {s_diluted} vs {s_clean}"
        );
    }

    #[test]
    fn test_laminar_flame_speed_decreases_with_pressure() {
        let model = LaminarFlameSpeed::methane_air();
        let s_atm = model.speed(400.0, 101325.0, 0.0);
        let s_high = model.speed(400.0, 5e5, 0.0);
        // beta < 0 → higher pressure → lower speed
        assert!(s_high < s_atm, "Higher pressure should reduce flame speed");
    }

    #[test]
    fn test_laminar_flame_speed_thickness_positive() {
        let model = LaminarFlameSpeed::methane_air();
        let delta = model.flame_thickness(1.5e-5, 400.0, 101325.0);
        assert!(delta > 0.0, "Flame thickness should be positive: {delta}");
    }

    // ─── ElementaryReactionLbm tests ───

    #[test]
    fn test_elementary_reaction_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut rxn = ElementaryReactionLbm::new(nx, ny, 0.1, 0.05, 0.05, 0.02);
        let n = nx * ny;
        for k in 0..n {
            rxn.c_a[k] = 1.0;
            rxn.c_b[k] = 1.0;
        }
        // A + C = const (per cell)
        let a_plus_c_before = rxn.a_plus_c_total();
        for _ in 0..10 {
            rxn.step(0.01, 1.0);
        }
        let a_plus_c_after = rxn.a_plus_c_total();
        assert!(
            (a_plus_c_before - a_plus_c_after).abs() < 1e-8,
            "A + C should be conserved: before={a_plus_c_before}, after={a_plus_c_after}"
        );
    }

    #[test]
    fn test_elementary_reaction_product_formed() {
        let nx = 4;
        let ny = 4;
        let mut rxn = ElementaryReactionLbm::new(nx, ny, 0.5, 0.01, 0.01, 0.005);
        let n = nx * ny;
        for k in 0..n {
            rxn.c_a[k] = 1.0;
            rxn.c_b[k] = 1.0;
        }
        rxn.step(0.1, 1.0);
        let c_total = rxn.total_c();
        assert!(c_total > 0.0, "Product C should form: {c_total}");
    }

    #[test]
    fn test_elementary_reaction_reactants_decrease() {
        let nx = 4;
        let ny = 4;
        let mut rxn = ElementaryReactionLbm::new(nx, ny, 0.5, 0.01, 0.01, 0.005);
        let n = nx * ny;
        for k in 0..n {
            rxn.c_a[k] = 2.0;
            rxn.c_b[k] = 2.0;
        }
        let a_before = rxn.total_a();
        rxn.step(0.1, 1.0);
        let a_after = rxn.total_a();
        assert!(
            a_after < a_before,
            "Reactant A should decrease: {a_after} vs {a_before}"
        );
    }

    // ─── SemenovExplosion tests ───

    #[test]
    fn test_semenov_critical_delta() {
        let delta_cr = SemenovExplosion::critical_delta();
        assert!(
            (delta_cr - 1.0 / std::f64::consts::E).abs() < 1e-12,
            "delta_cr = {delta_cr}"
        );
    }

    #[test]
    fn test_semenov_induction_time_below_critical() {
        let delta = 0.2; // < 1/e
        let t_ind = SemenovExplosion::induction_time(delta);
        assert!(
            t_ind.is_some(),
            "Should have induction time below critical delta"
        );
        assert!(t_ind.unwrap() > 0.0, "Induction time should be positive");
    }

    #[test]
    fn test_semenov_no_explosion_above_critical() {
        let delta = 0.5; // > 1/e ~ 0.368
        let t_ind = SemenovExplosion::induction_time(delta);
        assert!(t_ind.is_none(), "No explosion above critical delta");
    }

    #[test]
    fn test_semenov_adiabatic_rise_positive() {
        let delta_t =
            SemenovExplosion::adiabatic_temperature_rise(500_000.0, 1.0, 1.2, 1000.0, 300.0);
        assert!(
            delta_t > 0.0,
            "Adiabatic rise should be positive: {delta_t}"
        );
    }

    #[test]
    fn test_frank_kamenetskii_positive() {
        let fk = SemenovExplosion::frank_kamenetskii(500_000.0, 0.01, 0.01, 80_000.0, 0.1, 600.0);
        assert!(
            fk > 0.0,
            "Frank-Kamenetskii number should be positive: {fk}"
        );
    }

    // ─── MultiSpeciesDiffusion tests ───

    #[test]
    fn test_multi_species_diffusion_conservation() {
        let nx = 6;
        let ny = 6;
        let mut msd = MultiSpeciesDiffusion::new(nx, ny, vec![0.05, 0.03]);
        let n = nx * ny;
        msd.set_concentration(0, 3 * nx + 3, 1.0);
        msd.set_concentration(1, 2 * nx + 2, 0.5);
        let total0_before = msd.total_concentration(0);
        let total1_before = msd.total_concentration(1);
        let zero_vel = vec![[0.0_f64; 2]; n];
        for _ in 0..5 {
            msd.step(&zero_vel);
        }
        let total0_after = msd.total_concentration(0);
        let total1_after = msd.total_concentration(1);
        assert!(
            (total0_before - total0_after).abs() < 1e-9,
            "Species 0 total should be conserved: {total0_before} vs {total0_after}"
        );
        assert!(
            (total1_before - total1_after).abs() < 1e-9,
            "Species 1 total should be conserved: {total1_before} vs {total1_after}"
        );
    }

    #[test]
    fn test_multi_species_diffusion_source_increases_total() {
        let nx = 4;
        let ny = 4;
        let n = nx * ny;
        let mut msd = MultiSpeciesDiffusion::new(nx, ny, vec![0.05]);
        msd.set_sources(0, vec![0.1; n]);
        let before = msd.total_concentration(0);
        let zero_vel = vec![[0.0_f64; 2]; n];
        msd.step(&zero_vel);
        let after = msd.total_concentration(0);
        assert!(
            after > before,
            "Source term should increase total: before={before}, after={after}"
        );
    }

    // ─── ZeroDReactor tests ───

    #[test]
    fn test_zero_d_reactor_temperature_increases() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        let mut reactor = ZeroDReactor::new(1500.0, 1.0, 1.0, 200_000.0, 1000.0, arr);
        let t0 = reactor.temperature;
        reactor.step(1e-6);
        assert!(
            reactor.temperature > t0,
            "Temperature should increase in exothermic reaction"
        );
    }

    #[test]
    fn test_zero_d_reactor_fuel_decreases() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        let mut reactor = ZeroDReactor::new(1500.0, 2.0, 2.0, 200_000.0, 1000.0, arr);
        let c0 = reactor.c_fuel;
        reactor.step(1e-6);
        assert!(
            reactor.c_fuel <= c0,
            "Fuel should not increase: {} vs {}",
            reactor.c_fuel,
            c0
        );
    }

    #[test]
    fn test_zero_d_reactor_ignition_delay_found() {
        let arr = ArrheniusRate::new(1e14, 80_000.0);
        let mut reactor = ZeroDReactor::new(1200.0, 1.0, 1.0, 500_000.0, 500.0, arr);
        let t_delay = reactor.ignition_delay(1500.0, 1e-7, 0.01);
        assert!(t_delay.is_some(), "Ignition should occur within max_time");
        let tau = t_delay.unwrap();
        assert!(tau > 0.0, "Ignition delay should be positive: {tau}");
    }

    #[test]
    fn test_zero_d_reactor_no_ignition_if_cold() {
        let arr = ArrheniusRate::new(1e6, 200_000.0);
        let mut reactor = ZeroDReactor::new(300.0, 1.0, 1.0, 10_000.0, 5000.0, arr);
        // Very cold start: ignition target 2000K unreachable
        let t_delay = reactor.ignition_delay(2000.0, 1e-4, 1.0);
        assert!(
            t_delay.is_none(),
            "No ignition expected at T=300K with high Ea"
        );
    }
}

mod reactive_extended_tests_2 {
    use super::*;

    // --- ConcentrationField basics ---

    #[test]
    fn test_concentration_field_new_uniform() {
        let cf = ConcentrationField::new(5, 4, 0.01, 1.0);
        for j in 0..4 {
            for i in 0..5 {
                assert!((cf.get(i, j) - 1.0).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_concentration_field_set_get() {
        let mut cf = ConcentrationField::new(5, 5, 0.01, 0.0);
        cf.set(2, 3, 2.5);
        assert!((cf.get(2, 3) - 2.5).abs() < 1e-14);
    }

    #[test]
    fn test_concentration_field_total_mass_uniform() {
        let cf = ConcentrationField::new(4, 4, 0.01, 1.0);
        assert!((cf.total_mass() - 16.0).abs() < 1e-13);
    }

    #[test]
    fn test_concentration_field_diffusion_step_conserves_mass_periodic() {
        let mut cf = ConcentrationField::new(8, 8, 0.05, 1.0);
        cf.set(4, 4, 2.0);
        let mass_before = cf.total_mass();
        cf.diffusion_step(0.1, 1.0);
        let mass_after = cf.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Diffusion mass: {mass_before} vs {mass_after}"
        );
    }

    #[test]
    fn test_concentration_field_max_min() {
        let mut cf = ConcentrationField::new(4, 4, 0.01, 0.5);
        cf.set(1, 1, 2.0);
        cf.set(3, 3, 0.1);
        assert!((cf.max_concentration() - 2.0).abs() < 1e-14);
        assert!((cf.min_concentration() - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_concentration_field_clamp_min() {
        let mut cf = ConcentrationField::new(4, 4, 0.01, 0.5);
        cf.set(2, 2, -0.5);
        cf.clamp_min(0.0);
        assert!(
            cf.get(2, 2) >= 0.0,
            "Clamp should remove negative concentration"
        );
    }

    // --- first_order_reaction ---

    #[test]
    fn test_first_order_reaction_decays_exponentially() {
        let mut cf = ConcentrationField::new(4, 4, 0.01, 1.0);
        let k: f64 = 0.1;
        let dt: f64 = 1.0;
        first_order_reaction(&mut cf, k, dt);
        let expected = (-k * dt).exp();
        assert!(
            (cf.get(0, 0) - expected).abs() < 1e-13,
            "First order decay: {} vs {}",
            cf.get(0, 0),
            expected
        );
    }

    #[test]
    fn test_first_order_reaction_zero_rate_no_change() {
        let mut cf = ConcentrationField::new(4, 4, 0.01, 1.5);
        first_order_reaction(&mut cf, 0.0, 1.0);
        assert!((cf.get(2, 2) - 1.5).abs() < 1e-14);
    }

    // --- bimolecular_reaction ---

    #[test]
    fn test_bimolecular_reaction_reduces_concentrations() {
        let mut ca = ConcentrationField::new(4, 4, 0.01, 1.0);
        let mut cb = ConcentrationField::new(4, 4, 0.01, 1.0);
        bimolecular_reaction(&mut ca, &mut cb, 0.1, 0.01);
        assert!(ca.get(0, 0) < 1.0, "A should decrease");
        assert!(cb.get(0, 0) < 1.0, "B should decrease");
    }

    #[test]
    fn test_bimolecular_reaction_zero_rate_no_change() {
        let mut ca = ConcentrationField::new(4, 4, 0.01, 2.0);
        let mut cb = ConcentrationField::new(4, 4, 0.01, 3.0);
        bimolecular_reaction(&mut ca, &mut cb, 0.0, 0.01);
        assert!((ca.get(1, 1) - 2.0).abs() < 1e-14);
        assert!((cb.get(1, 1) - 3.0).abs() < 1e-14);
    }

    // --- ArrheniusRate ---

    #[test]
    fn test_arrhenius_rate_increases_with_temperature() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        let r1 = arr.rate(800.0);
        let r2 = arr.rate(1200.0);
        assert!(
            r2 > r1,
            "Arrhenius rate should increase with T: r1={r1}, r2={r2}"
        );
    }

    #[test]
    fn test_arrhenius_rate_zero_at_zero_temperature() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        // Rate at T→0 should be essentially zero
        let r = arr.rate(1.0);
        assert!(r < 1e-100, "Arrhenius rate at T=1 K should be ~0: {r}");
    }

    #[test]
    fn test_arrhenius_rate_ratio_positive() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        // rate_ratio(t1, t2) = exp(Ea/R * (1/t1 - 1/t2))
        // To get rate(high_T)/rate(low_T), pick t1=high, t2=low → 1/t1 < 1/t2 → result < 1
        // To get > 1: t1=low, t2=high (rate at lower T relative to higher T is less)
        // Actually we just verify the ratio is finite and positive
        let ratio = arr.rate_ratio(1000.0, 500.0);
        assert!(ratio > 0.0, "Rate ratio should be positive: {ratio}");
    }

    // --- BimolecularRate ---

    #[test]
    fn test_bimolecular_rate_proportional_to_concentrations() {
        let br = BimolecularRate::new(0.5);
        let r = br.rate(2.0, 3.0);
        assert!(
            (r - 3.0).abs() < 1e-14,
            "BimolecularRate: 0.5*2*3=3, got {r}"
        );
    }

    #[test]
    fn test_bimolecular_rate_zero_if_either_zero() {
        let br = BimolecularRate::new(10.0);
        assert!(br.rate(0.0, 5.0).abs() < 1e-15);
        assert!(br.rate(5.0, 0.0).abs() < 1e-15);
    }

    // --- MultiStepReaction ---

    #[test]
    fn test_multi_step_reaction_decreases_first_species() {
        let msr = MultiStepReaction::new(vec![0, 1, 2], vec![0.1, 0.05]);
        let mut concentrations = vec![vec![1.0; 1], vec![0.5; 1], vec![0.0; 1]];
        let c0_before = concentrations[0][0];
        msr.apply(&mut concentrations, 0.01);
        assert!(
            concentrations[0][0] < c0_before,
            "First species should decrease"
        );
    }

    #[test]
    fn test_multi_step_reaction_num_steps() {
        let msr = MultiStepReaction::new(vec![0, 1, 2], vec![0.1, 0.05]);
        assert_eq!(msr.num_steps(), 2);
    }

    // --- Species struct ---

    #[test]
    fn test_species_new() {
        let sp = Species::new("H2".to_string(), 0.01, 2.016);
        assert_eq!(sp.name, "H2");
        assert!((sp.diffusivity - 0.01).abs() < 1e-15);
        assert!((sp.molar_mass - 2.016).abs() < 1e-12);
    }

    // --- Damkoehler and Peclet numbers ---

    #[test]
    fn test_damkoehler_number_positive() {
        let da = damkoehler_number(0.1, 0.1, 0.01);
        assert!(da > 0.0, "Damkoehler should be positive: {da}");
    }

    #[test]
    fn test_peclet_number_positive() {
        let pe = peclet_number(0.5, 0.1, 0.01);
        assert!(pe > 0.0, "Peclet should be positive: {pe}");
    }

    #[test]
    fn test_peclet_number_value() {
        // Pe = u*L/D = 0.5*0.1/0.01 = 5
        let pe = peclet_number(0.5, 0.1, 0.01);
        assert!((pe - 5.0).abs() < 1e-13, "Peclet = {pe}");
    }

    // --- LaminarFlameSpeed ---

    #[test]
    fn test_laminar_flame_speed_positive() {
        let lfs = LaminarFlameSpeed::methane_air();
        let s = lfs.speed(600.0, 101325.0, 0.0);
        assert!(s > 0.0, "Flame speed should be positive: {s}");
    }

    #[test]
    fn test_laminar_flame_speed_zero_diluent_equals_base_conditions() {
        let lfs = LaminarFlameSpeed::methane_air();
        let s = lfs.speed(298.0, 101325.0, 0.0);
        assert!((s - 0.37).abs() < 1e-10, "At reference conditions: s={s}");
    }

    #[test]
    fn test_flame_thickness_decreases_with_higher_speed() {
        let lfs = LaminarFlameSpeed::methane_air();
        let alpha = 2e-5;
        let thick1 = lfs.flame_thickness(alpha, 600.0, 101325.0);
        let thick2 = lfs.flame_thickness(alpha, 1200.0, 101325.0);
        assert!(
            thick2 < thick1,
            "Thickness should decrease at higher T: {thick1} vs {thick2}"
        );
    }

    #[test]
    fn test_speed_temperature_sensitivity_positive() {
        let lfs = LaminarFlameSpeed::methane_air();
        let sens = lfs.speed_temperature_sensitivity(600.0, 101325.0);
        assert!(
            sens > 0.0,
            "Temperature sensitivity should be positive: {sens}"
        );
    }

    // --- SemenovExplosion ---

    #[test]
    fn test_semenov_critical_delta_value() {
        let delta_cr = SemenovExplosion::critical_delta();
        let expected = 1.0 / std::f64::consts::E;
        assert!((delta_cr - expected).abs() < 1e-14, "delta_cr={delta_cr}");
    }

    #[test]
    fn test_semenov_induction_time_below_critical() {
        let delta = 0.1; // well below critical
        let t_ind = SemenovExplosion::induction_time(delta);
        assert!(
            t_ind.is_some(),
            "Should have induction time below critical delta"
        );
        assert!(t_ind.unwrap() > 0.0, "Induction time should be positive");
    }

    #[test]
    fn test_semenov_adiabatic_rise_linearly_proportional_to_heat() {
        let rise1 = SemenovExplosion::adiabatic_temperature_rise(1.0, 1.0, 1.0, 1.0, 1.0);
        let rise2 = SemenovExplosion::adiabatic_temperature_rise(2.0, 1.0, 1.0, 1.0, 1.0);
        assert!(
            (rise2 - 2.0 * rise1).abs() < 1e-13,
            "Linear in Q: {rise1} vs {rise2}"
        );
    }

    // --- MultiSpeciesDiffusion ---

    #[test]
    fn test_multi_species_diffusion_zero_initial_concentration() {
        let msd = MultiSpeciesDiffusion::new(4, 4, vec![0.05]);
        let total = msd.total_concentration(0);
        assert!(total.abs() < 1e-14, "Initial total = {total}");
    }

    #[test]
    fn test_multi_species_diffusion_set_concentration_scalar() {
        let mut msd = MultiSpeciesDiffusion::new(4, 4, vec![0.05]);
        msd.set_concentration(0, 5, 2.0);
        assert!(
            (msd.concentration(0, 5) - 2.0).abs() < 1e-13,
            "Concentration at cell 5: {}",
            msd.concentration(0, 5)
        );
    }

    #[test]
    fn test_multi_species_n_species() {
        let msd = MultiSpeciesDiffusion::new(4, 4, vec![0.05, 0.03, 0.01]);
        assert_eq!(msd.n_species, 3);
    }

    #[test]
    fn test_multi_species_step_preserves_zero_velocity_total() {
        let nx = 4;
        let ny = 4;
        let n = nx * ny;
        let mut msd = MultiSpeciesDiffusion::new(nx, ny, vec![0.05]);
        msd.set_concentration(0, 5, 1.0);
        let before = msd.total_concentration(0);
        let zero_vel = vec![[0.0_f64; 2]; n];
        for _ in 0..3 {
            msd.step(&zero_vel);
        }
        let after = msd.total_concentration(0);
        assert!(
            (before - after).abs() < 1e-10,
            "MultiSpeciesDiffusion conservation: {before} vs {after}"
        );
    }

    // --- ZeroDReactor ---

    #[test]
    fn test_zero_d_reactor_step_positive_temperature() {
        let arr = ArrheniusRate::new(1e10, 60_000.0);
        let mut reactor = ZeroDReactor::new(1000.0, 1.0, 1.0, 100_000.0, 1000.0, arr);
        reactor.step(1e-5);
        assert!(reactor.temperature > 0.0);
    }

    #[test]
    fn test_zero_d_reactor_fuel_non_negative() {
        let arr = ArrheniusRate::new(1e12, 80_000.0);
        let mut reactor = ZeroDReactor::new(1500.0, 1.0, 1.0, 200_000.0, 1000.0, arr);
        for _ in 0..100 {
            reactor.step(1e-7);
        }
        assert!(
            reactor.c_fuel >= 0.0,
            "Fuel cannot go negative: {}",
            reactor.c_fuel
        );
    }

    // --- HeatRelease ---

    #[test]
    fn test_heat_release_mean_temperature_initial() {
        let hr = HeatRelease::new(4, 4, 200_000.0, 300.0);
        assert!((hr.mean_temperature() - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_heat_release_max_temperature_initial_equals_t0() {
        let hr = HeatRelease::new(4, 4, 200_000.0, 500.0);
        assert!((hr.max_temperature() - 500.0).abs() < 1e-12);
    }

    #[test]
    fn test_heat_release_apply_first_order_increases_temperature() {
        let mut hr = HeatRelease::new(4, 4, 200_000.0, 300.0);
        let conc = vec![1.0; 16];
        hr.apply_first_order(&conc, 0.1, 1e-3);
        assert!(
            hr.mean_temperature() > 300.0,
            "Heat release should raise temperature"
        );
    }

    // --- LbmPassiveScalar ---

    #[test]
    fn test_lbm_passive_scalar_total_zero_initial() {
        let ps = LbmPassiveScalar::new(4, 4, 0.01);
        let total = ps.total_concentration();
        assert!(total.abs() < 1e-14, "Initial passive scalar = {total}");
    }

    #[test]
    fn test_lbm_passive_scalar_initialize_uniform() {
        let mut ps = LbmPassiveScalar::new(4, 4, 0.01);
        let conc = vec![1.0; 16];
        ps.initialize_concentration(&conc);
        assert!((ps.total_concentration() - 16.0).abs() < 1e-12);
    }

    #[test]
    fn test_lbm_passive_scalar_add_source_increases_total() {
        let mut ps = LbmPassiveScalar::new(4, 4, 0.01);
        let sources = vec![0.1; 16];
        ps.add_source(&sources);
        // After adding source, total should increase
        assert!(
            ps.total_concentration() > 0.0,
            "Source should increase total scalar"
        );
    }

    #[test]
    fn test_lbm_passive_scalar_equilibrium_sums_to_c() {
        let c = 2.0;
        let u = [0.04, -0.02];
        let geq = LbmPassiveScalar::equilibrium(c, u);
        let sum: f64 = geq.iter().sum();
        assert!((sum - c).abs() < 1e-13, "Passive scalar eq sum = {sum}");
    }
}

mod extended_reactive_tests {
    use super::*;

    // ---- ArrheniusKinetics ----

    #[test]
    fn test_arrhenius_rate_zero_at_zero_temp() {
        let kin = ArrheniusKinetics::new(1e10, 5000.0, 8.314, 0.0);
        assert_eq!(kin.rate(0.0), 0.0);
    }

    #[test]
    fn test_arrhenius_rate_positive_at_finite_temp() {
        let kin = ArrheniusKinetics::new(1e10, 5000.0, 8.314, 0.0);
        assert!(kin.rate(1000.0) > 0.0);
    }

    #[test]
    fn test_arrhenius_rate_increases_with_temp() {
        let kin = ArrheniusKinetics::new(1e10, 5000.0, 8.314, 0.0);
        assert!(kin.rate(1200.0) > kin.rate(1000.0));
    }

    #[test]
    fn test_arrhenius_d_rate_d_temp_positive() {
        let kin = ArrheniusKinetics::new(1e10, 5000.0, 8.314, 0.0);
        assert!(kin.d_rate_d_temp(1000.0) > 0.0);
    }

    #[test]
    fn test_arrhenius_activation_temperature() {
        let kin = ArrheniusKinetics::new(1.0, 8314.0, 8.314, 0.0);
        let ta = kin.activation_temperature();
        assert!((ta - 1000.0).abs() < 1e-6, "Ta = {ta}");
    }

    #[test]
    fn test_arrhenius_beta_increases_rate_at_high_t() {
        let kin0 = ArrheniusKinetics::new(1.0, 1000.0, 8.314, 0.0);
        let kin1 = ArrheniusKinetics::new(1.0, 1000.0, 8.314, 1.0);
        assert!(kin1.rate(1000.0) > kin0.rate(1000.0));
    }

    // ---- MultiComponentField ----

    #[test]
    fn test_multi_component_field_total_conc_initial() {
        let mcf = MultiComponentField::new(4, 4, 2, &[1.0, 0.5], 300.0, &[1e-5, 2e-5], 1e-5);
        assert!((mcf.total_concentration(0) - 16.0).abs() < 1e-12);
        assert!((mcf.total_concentration(1) - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_multi_component_field_mean_temp_initial() {
        let mcf = MultiComponentField::new(4, 4, 1, &[1.0], 500.0, &[1e-5], 1e-5);
        assert!((mcf.mean_temperature() - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_component_diffusion_conserves_total() {
        let mut mcf = MultiComponentField::new(4, 4, 1, &[1.0], 300.0, &[0.1], 0.1);
        let total_before = mcf.total_concentration(0);
        mcf.diffuse_species(0.01, 1.0);
        let total_after = mcf.total_concentration(0);
        assert!((total_before - total_after).abs() < 1e-10);
    }

    #[test]
    fn test_multi_component_temperature_diffusion_conserves() {
        let mut mcf = MultiComponentField::new(4, 4, 1, &[1.0], 300.0, &[0.1], 0.1);
        let sum_before: f64 = mcf.temp.iter().sum();
        mcf.diffuse_temperature(0.01, 1.0);
        let sum_after: f64 = mcf.temp.iter().sum();
        assert!((sum_before - sum_after).abs() < 1e-9);
    }

    #[test]
    fn test_multi_component_bimolecular_decreases_both() {
        let mut mcf = MultiComponentField::new(4, 4, 2, &[1.0, 1.0], 300.0, &[1e-5, 1e-5], 1e-5);
        mcf.apply_bimolecular_reaction(0, 1, 0.5, 0.1);
        assert!(mcf.total_concentration(0) < 16.0);
        assert!(mcf.total_concentration(1) < 16.0);
    }

    #[test]
    fn test_multi_component_heat_release_increases_temp() {
        let mut mcf = MultiComponentField::new(4, 4, 2, &[1.0, 1.0], 300.0, &[1e-5, 1e-5], 1e-5);
        let t_before = mcf.mean_temperature();
        mcf.apply_heat_release_bimolecular(0, 1, 0.5, 100.0, 0.1);
        let t_after = mcf.mean_temperature();
        assert!(
            t_after > t_before,
            "heat release should warm: {t_before} → {t_after}"
        );
    }

    #[test]
    fn test_arrhenius_reaction_decreases_concentration() {
        let mut mcf = MultiComponentField::new(2, 2, 1, &[1.0], 1000.0, &[1e-5], 1e-5);
        let c_before = mcf.total_concentration(0);
        mcf.apply_arrhenius_reaction(0, 1e6, 5000.0, 8.314, 0.01);
        let c_after = mcf.total_concentration(0);
        assert!(c_after < c_before, "Arrhenius should consume species");
    }

    #[test]
    fn test_hotspot_index_at_maximum() {
        let mut mcf = MultiComponentField::new(4, 4, 1, &[1.0], 300.0, &[1e-5], 1e-5);
        mcf.temp[7] = 1000.0;
        assert_eq!(mcf.hotspot_index(), 7);
    }

    // ---- SpeciesLbm2D ----

    #[test]
    fn test_species_lbm2d_total_zero_initial() {
        let sp = SpeciesLbm2D::new(4, 4, 1.5);
        assert!(sp.total_concentration().abs() < 1e-14);
    }

    #[test]
    fn test_species_lbm2d_initialize_uniform() {
        let mut sp = SpeciesLbm2D::new(4, 4, 1.5);
        sp.initialize_uniform(1.0, 0.0, 0.0);
        assert!((sp.total_concentration() - 16.0).abs() < 1e-12);
    }

    #[test]
    fn test_species_lbm2d_equilibrium_sums_to_c() {
        let geq = SpeciesLbm2D::equilibrium(2.0, 0.05, -0.02);
        let sum: f64 = geq.iter().sum();
        assert!((sum - 2.0).abs() < 1e-13, "geq sum = {sum}");
    }

    #[test]
    fn test_species_lbm2d_collide_conserves_total() {
        let nx = 4usize;
        let ny = 4usize;
        let ux = vec![0.05; nx * ny];
        let uy = vec![0.0; nx * ny];
        let mut sp = SpeciesLbm2D::new(nx, ny, 1.5);
        sp.initialize_uniform(1.0, 0.05, 0.0);
        let before = sp.total_concentration();
        sp.collide(&ux, &uy);
        let after = sp.total_concentration();
        assert!((before - after).abs() < 1e-12);
    }

    #[test]
    fn test_species_lbm2d_stream_conserves_total() {
        let nx = 4usize;
        let ny = 4usize;
        let mut sp = SpeciesLbm2D::new(nx, ny, 1.5);
        sp.initialize_uniform(1.0, 0.0, 0.0);
        let before = sp.total_concentration();
        sp.stream();
        let after = sp.total_concentration();
        assert!((before - after).abs() < 1e-12);
    }

    #[test]
    fn test_species_lbm2d_add_source_increases_total() {
        let nx = 4usize;
        let ny = 4usize;
        let mut sp = SpeciesLbm2D::new(nx, ny, 1.5);
        let source = vec![0.1; nx * ny];
        sp.add_source(&source, 1.0);
        assert!(sp.total_concentration() > 0.0);
    }

    // ---- TwoSpeciesReaction ----

    #[test]
    fn test_two_species_initial_total_zero() {
        let kin = ArrheniusKinetics::new(1.0, 0.0, 8.314, 0.0);
        let rxn = TwoSpeciesReaction::new(4, 4, 1.5, kin, 300.0);
        assert!(rxn.total_a().abs() < 1e-14);
        assert!(rxn.total_b().abs() < 1e-14);
    }

    #[test]
    fn test_two_species_react_decreases_total() {
        let kin = ArrheniusKinetics::new(1e6, 5000.0, 8.314, 0.0);
        let mut rxn = TwoSpeciesReaction::new(4, 4, 1.5, kin, 1000.0);
        rxn.species_a.initialize_uniform(1.0, 0.0, 0.0);
        rxn.species_b.initialize_uniform(1.0, 0.0, 0.0);
        let ta = rxn.total_a();
        let tb = rxn.total_b();
        rxn.react(0.01);
        assert!(rxn.total_a() < ta, "A should decrease");
        assert!(rxn.total_b() < tb, "B should decrease");
    }

    #[test]
    fn test_two_species_step_runs() {
        // Use tiny rate constant and small dt so reaction is weak and totals stay finite
        let kin = ArrheniusKinetics::new(0.01, 1000.0, 8.314, 0.0);
        let mut rxn = TwoSpeciesReaction::new(4, 4, 1.5, kin, 500.0);
        rxn.species_a.initialize_uniform(0.5, 0.0, 0.0);
        rxn.species_b.initialize_uniform(0.5, 0.0, 0.0);
        let nx = 4usize;
        let ny = 4usize;
        let ux = vec![0.02; nx * ny];
        let uy = vec![0.0; nx * ny];
        rxn.step(&ux, &uy, 0.001);
        // Verify step ran without panic and totals are finite
        assert!(rxn.total_a().is_finite());
        assert!(rxn.total_b().is_finite());
    }

    // ---- GrayScott ----

    #[test]
    fn test_gray_scott_initial_u_equals_area() {
        let gs = GrayScott::new(8, 8, 0.2, 0.1, 0.04, 0.06);
        assert!((gs.total_u() - 64.0).abs() < 1e-12);
    }

    #[test]
    fn test_gray_scott_initial_v_zero() {
        let gs = GrayScott::new(8, 8, 0.2, 0.1, 0.04, 0.06);
        assert!(gs.total_v().abs() < 1e-14);
    }

    #[test]
    fn test_gray_scott_seed_increases_v() {
        let mut gs = GrayScott::new(8, 8, 0.2, 0.1, 0.04, 0.06);
        gs.seed_centre(1, 0.25);
        assert!(gs.total_v() > 0.0, "seed should add V");
    }

    #[test]
    fn test_gray_scott_step_runs() {
        let mut gs = GrayScott::new(8, 8, 0.2, 0.1, 0.04, 0.06);
        gs.seed_centre(1, 0.25);
        for _ in 0..10 {
            gs.step(0.1);
        }
        // Check fields stay in [0,1]
        for &u in &gs.u {
            assert!((0.0..=1.0).contains(&u), "u out of range: {u}");
        }
        for &v in &gs.v {
            assert!((0.0..=1.0).contains(&v), "v out of range: {v}");
        }
    }

    #[test]
    fn test_gray_scott_mean_u_positive() {
        let gs = GrayScott::new(4, 4, 0.2, 0.1, 0.04, 0.06);
        assert!(gs.mean_u() > 0.0);
    }

    // ---- ReactiveFlowSolver ----

    #[test]
    fn test_reactive_flow_solver_initial_species_zero() {
        let solver = ReactiveFlowSolver::new(4, 4, 1.5, 2, 1.5, vec![0.1, 0.2]);
        assert!(solver.total_species(0).abs() < 1e-14);
        assert!(solver.total_species(1).abs() < 1e-14);
    }

    #[test]
    fn test_reactive_flow_solver_step_runs() {
        let mut solver = ReactiveFlowSolver::new(4, 4, 1.5, 1, 1.5, vec![0.0]);
        solver.species[0].initialize_uniform(1.0, 0.0, 0.0);
        solver.step(0.01);
        assert!(solver.total_species(0) > 0.0);
    }

    #[test]
    fn test_reactive_flow_solver_reaction_reduces_species() {
        let mut solver = ReactiveFlowSolver::new(4, 4, 1.5, 1, 1.5, vec![1.0]);
        solver.species[0].initialize_uniform(1.0, 0.0, 0.0);
        let before = solver.total_species(0);
        solver.step(0.1);
        let after = solver.total_species(0);
        assert!(
            after < before,
            "reaction should reduce species: {before} → {after}"
        );
    }

    // ---- Dimensionless numbers ----

    #[test]
    fn test_zeldovich_number_positive() {
        let ze = zeldovich_number(5000.0 * 8.314, 2000.0, 300.0, 8.314);
        assert!(ze > 0.0, "Ze = {ze}");
    }

    #[test]
    fn test_lewis_number_ratio() {
        let le = lewis_number(2e-5, 1e-5);
        assert!((le - 2.0).abs() < 1e-12, "Le = {le}");
    }

    #[test]
    fn test_schmidt_number_positive() {
        let sc = schmidt_number(1e-6, 1e-9);
        assert!(sc > 0.0 && sc.is_finite());
    }

    #[test]
    fn test_progress_variable_zero_to_one() {
        let chi = progress_variable(0.5, 0.0, 1.0);
        assert!((chi - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_progress_variable_zero_denom() {
        let chi = progress_variable(0.5, 0.5, 0.5);
        assert_eq!(chi, 0.0);
    }

    #[test]
    fn test_equivalence_ratio_stoichiometric() {
        let phi = equivalence_ratio(1.0, 1.0, 1.0);
        assert!((phi - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_equivalence_ratio_rich() {
        let phi = equivalence_ratio(2.0, 1.0, 1.0);
        assert!((phi - 2.0).abs() < 1e-14);
    }
}
