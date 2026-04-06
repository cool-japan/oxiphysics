//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::NortonCreep;
    use crate::anisotropic::PuckFailureCriteria;
    use crate::damage::BaileyNortonCreep;
    use crate::damage::CdmDamage;
    use crate::damage::CrackBandModel;
    use crate::damage::DamageEvolutionLaw;
    use crate::damage::DamageInducedAnisotropy;
    use crate::damage::GradientEnhancedDamage;
    use crate::damage::GursonDamage;
    use crate::damage::HashinFailureCriteria;
    use crate::damage::IsotropicDamage;
    use crate::damage::LemaitreDuctileDamage;
    use crate::damage::MazarsDamage;
    use crate::damage::MinerFatigueDamage;
    use crate::damage::NonLocalRegularization;
    use crate::damage::ProgressiveDamage;
    use crate::damage::ProgressiveFailureAnalysis;
    use crate::damage::RepairEffectiveness;
    use crate::damage::TsaiWuCriterion;
    #[test]
    fn test_damage_undamaged() {
        let dmg = IsotropicDamage::new(1.0e4, 1.0);
        let e = 200.0e9;
        assert!(
            (dmg.effective_stiffness(e) - e).abs() < 1.0e-6,
            "Undamaged stiffness should equal E"
        );
    }
    #[test]
    fn test_damage_initiation() {
        let y0 = 1.0e4_f64;
        let e = 200.0e9_f64;
        let mut dmg = IsotropicDamage::new(y0, 1.0);
        let sigma = 1.0e6;
        dmg.update(sigma, e);
        assert_eq!(dmg.d, 0.0, "Damage should be zero below threshold");
    }
    #[test]
    fn test_damage_evolution() {
        let y0 = 1.0e4_f64;
        let e = 200.0e9_f64;
        let mut dmg = IsotropicDamage::new(y0, 1.0);
        let sigma = 1.0e8;
        dmg.update(sigma, e);
        assert!(
            dmg.d > 0.0,
            "Damage should increase above threshold, got {}",
            dmg.d
        );
    }
    #[test]
    fn test_damage_fully_damaged() {
        let e = 200.0e9_f64;
        let mut dmg = IsotropicDamage::new(1.0e4, 1.0);
        dmg.d = 1.0 - 1.0e-10;
        let e_eff = dmg.effective_stiffness(e);
        assert!(
            e_eff / e < 1.0e-8,
            "Fully damaged stiffness should approach 0, got {e_eff}"
        );
    }
    #[test]
    fn test_mazars_tension() {
        let eps = [1.0e-3_f64, 0.5e-3, 0.0];
        let eps_tilde = MazarsDamage::equivalent_strain(&eps);
        let expected = ((1.0e-3_f64).powi(2) + (0.5e-3_f64).powi(2)).sqrt();
        assert!(
            (eps_tilde - expected).abs() < 1.0e-15,
            "Equivalent strain mismatch: got {eps_tilde}, expected {expected}"
        );
        assert!(eps_tilde > 0.0, "Equivalent strain should be positive");
    }
    #[test]
    fn test_mazars_compression() {
        let eps = [-1.0e-3_f64, -0.5e-3, -0.2e-3];
        let eps_tilde = MazarsDamage::equivalent_strain(&eps);
        assert!(
            eps_tilde.abs() < 1.0e-15,
            "Equivalent strain under pure compression should be 0, got {eps_tilde}"
        );
    }
    #[test]
    fn test_norton_creep_rate() {
        let a = 1.0e-10_f64;
        let n = 3.0_f64;
        let q = 50_000.0_f64;
        let norton = NortonCreep::new(a, n, q);
        let sigma = 1.0e8_f64;
        let temp = 1000.0_f64;
        let rate = norton.creep_rate(sigma, temp);
        let expected = a * sigma.powf(n) * (-q / (8.314 * temp)).exp();
        assert!(
            (rate - expected).abs() / expected < 1.0e-10,
            "Norton creep rate mismatch: got {rate}, expected {expected}"
        );
    }
    #[test]
    fn test_norton_creep_temperature() {
        let norton = NortonCreep::new(1.0e-10, 3.0, 50_000.0);
        let sigma = 1.0e8;
        let rate_low = norton.creep_rate(sigma, 800.0);
        let rate_high = norton.creep_rate(sigma, 1200.0);
        assert!(
            rate_high > rate_low,
            "Higher temperature should give higher creep rate: {rate_high} vs {rate_low}"
        );
    }
    #[test]
    fn test_bailey_norton_strain() {
        let b = 2.0e-12_f64;
        let m = 2.5_f64;
        let p = 0.4_f64;
        let bn = BaileyNortonCreep::new(b, m, p);
        let sigma = 5.0e7_f64;
        let time = 100.0_f64;
        let strain = bn.creep_strain(sigma, time);
        let expected = b * sigma.powf(m) * time.powf(p);
        assert!(
            (strain - expected).abs() / expected < 1.0e-10,
            "Bailey-Norton strain mismatch: got {strain}, expected {expected}"
        );
    }
    #[test]
    fn test_bailey_norton_rate_decreases() {
        let bn = BaileyNortonCreep::new(2.0e-12, 2.5, 0.4);
        let sigma = 5.0e7;
        let rate_early = bn.creep_rate(sigma, 10.0);
        let rate_late = bn.creep_rate(sigma, 1000.0);
        assert!(
            rate_late < rate_early,
            "Creep rate should decrease with time (p<1): early={rate_early}, late={rate_late}"
        );
    }
    #[test]
    fn test_cdm_linear_softening() {
        let mut cdm = CdmDamage::new(0.001, 0.01, 0.0, DamageEvolutionLaw::Linear);
        cdm.update(0.0005);
        assert_eq!(cdm.d, 0.0);
        cdm.update(0.005);
        assert!(cdm.d > 0.0 && cdm.d < 1.0, "D = {}", cdm.d);
        cdm.update(0.01);
        assert!(
            (cdm.d - 1.0).abs() < 1e-10,
            "D should be 1.0 at failure strain, got {}",
            cdm.d
        );
    }
    #[test]
    fn test_cdm_exponential_softening() {
        let mut cdm = CdmDamage::new(0.001, 0.01, 100.0, DamageEvolutionLaw::Exponential);
        cdm.update(0.005);
        assert!(cdm.d > 0.0);
        cdm.update(0.1);
        assert!(
            cdm.d > 0.99,
            "D should be near 1 for large strains, got {}",
            cdm.d
        );
    }
    #[test]
    fn test_cdm_power_law_softening() {
        let mut cdm = CdmDamage::new(0.001, 0.01, 2.0, DamageEvolutionLaw::PowerLaw);
        cdm.update(0.0055);
        let d1 = cdm.d;
        assert!(
            (d1 - 0.25).abs() < 0.01,
            "D at midpoint should be ~0.25, got {d1}"
        );
    }
    #[test]
    fn test_cdm_irreversibility() {
        let mut cdm = CdmDamage::new(0.001, 0.01, 0.0, DamageEvolutionLaw::Linear);
        cdm.update(0.005);
        let d_load = cdm.d;
        cdm.update(0.002);
        assert_eq!(cdm.d, d_load, "Damage should not decrease on unloading");
    }
    #[test]
    fn test_lemaitre_no_damage_below_threshold() {
        let mut lem = LemaitreDuctileDamage::new(0.1, 2.0e6, 1.0, 0.5, 0.3);
        lem.update(300.0e6, 0.333, 200.0e9, 0.05);
        assert_eq!(
            lem.d, 0.0,
            "No damage expected below threshold plastic strain"
        );
    }
    #[test]
    fn test_lemaitre_damage_above_threshold() {
        let mut lem = LemaitreDuctileDamage::new(0.1, 2.0e6, 1.0, 0.5, 0.3);
        lem.eps_p = 0.15;
        lem.update(300.0e6, 0.333, 200.0e9, 0.05);
        assert!(
            lem.d > 0.0,
            "Damage should grow above threshold, got {}",
            lem.d
        );
    }
    #[test]
    fn test_lemaitre_critical_damage() {
        let mut lem = LemaitreDuctileDamage::new(0.0, 1.0e3, 1.0, 0.5, 0.3);
        for _ in 0..1000 {
            lem.update(500.0e6, 0.5, 200.0e9, 0.01);
        }
        assert!(lem.d <= lem.d_critical + 1e-10, "D should not exceed D_c");
    }
    #[test]
    fn test_lemaitre_triaxiality_function() {
        let lem = LemaitreDuctileDamage::new(0.1, 2.0e6, 1.0, 0.5, 0.3);
        let rv = lem.triaxiality_function(0.333);
        let expected = (2.0 / 3.0) * 1.3 + 3.0 * 0.4 * 0.333 * 0.333;
        assert!(
            (rv - expected).abs() < 1e-6,
            "R_v mismatch: got {rv}, expected {expected}"
        );
    }
    #[test]
    fn test_gurson_initial_state() {
        let gtn = GursonDamage::new(0.001, 0.01, 0.25, 0.04, 0.3, 0.1);
        assert!((gtn.f - 0.001).abs() < 1e-15);
        assert!((gtn.q1 - 1.5).abs() < 1e-15);
        assert!((gtn.q2 - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_gurson_yield_function_undamaged() {
        let mut gtn = GursonDamage::new(0.0, 0.01, 0.25, 0.04, 0.3, 0.1);
        gtn.f = 0.0;
        let phi = gtn.yield_function(300.0e6, 0.0, 300.0e6);
        assert!(
            phi.abs() < 1e-6,
            "Yield function should be 0 at yield, got {phi}"
        );
    }
    #[test]
    fn test_gurson_void_growth() {
        let gtn = GursonDamage::new(0.001, 0.01, 0.25, 0.04, 0.3, 0.1);
        let growth = gtn.growth_rate(0.01);
        let expected = (1.0 - 0.001) * 0.01;
        assert!((growth - expected).abs() < 1e-10);
    }
    #[test]
    fn test_gurson_nucleation() {
        let gtn = GursonDamage::new(0.001, 0.01, 0.25, 0.04, 0.3, 0.1);
        let rate_at_peak = gtn.nucleation_rate(0.3, 1.0);
        let rate_away = gtn.nucleation_rate(0.6, 1.0);
        assert!(rate_at_peak > rate_away, "Nucleation should peak at ε_N");
    }
    #[test]
    fn test_gurson_effective_void_fraction() {
        let mut gtn = GursonDamage::new(0.001, 0.01, 0.25, 0.04, 0.3, 0.1);
        assert!((gtn.effective_void_fraction() - 0.001).abs() < 1e-15);
        gtn.f = 0.15;
        let f_star = gtn.effective_void_fraction();
        assert!(f_star > gtn.fc, "f* should exceed fc after coalescence");
    }
    #[test]
    fn test_gurson_update() {
        let mut gtn = GursonDamage::new(0.001, 0.01, 0.25, 0.04, 0.3, 0.1);
        let f_init = gtn.f;
        gtn.update(0.3, 0.001, 0.005, 1.0);
        assert!(gtn.f > f_init, "Void fraction should grow");
    }
    #[test]
    fn test_crack_band_initiation_strain() {
        let cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.01);
        let eps0 = cb.initiation_strain();
        assert!((eps0 - 1.0e-4).abs() < 1e-10);
    }
    #[test]
    fn test_crack_band_failure_strain() {
        let cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.01);
        let eps_f = cb.failure_strain();
        let expected = 2.0 * 100.0 / (3.0e6 * 0.01);
        assert!((eps_f - expected).abs() < 1e-10);
    }
    #[test]
    fn test_crack_band_snapback() {
        let cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.001);
        assert!(!cb.has_snapback(), "Small element should not have snapback");
        let cb_large = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 100.0);
        assert!(
            cb_large.has_snapback(),
            "Large element should have snapback"
        );
    }
    #[test]
    fn test_crack_band_max_element_size() {
        let cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.01);
        let h_max = cb.max_element_size();
        let expected = 2.0 * 30.0e9 * 100.0 / (3.0e6 * 3.0e6);
        assert!((h_max - expected).abs() < 1e-6);
    }
    #[test]
    fn test_crack_band_damage_evolution() {
        let mut cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.001);
        cb.update(0.5e-4);
        assert_eq!(cb.d, 0.0);
        cb.update(3.0e-3);
        assert!(cb.d > 0.0, "D should be > 0, got {}", cb.d);
        cb.update(1.0);
        assert!(
            (cb.d - 1.0).abs() < 1e-10,
            "D should be 1.0 at very large strain"
        );
    }
    #[test]
    fn test_crack_band_energy_dissipation() {
        let mut cb = CrackBandModel::new(100.0, 3.0e6, 30.0e9, 0.001);
        cb.update(1.0);
        let g_density = cb.dissipated_energy_density();
        let g_total = g_density * cb.h;
        assert!(
            (g_total - cb.gf).abs() / cb.gf < 0.05,
            "Dissipated energy {g_total} should ≈ G_f {}",
            cb.gf
        );
    }
    #[test]
    fn test_nonlocal_gaussian_weight() {
        let nl = NonLocalRegularization::new(1.0);
        assert!(
            (nl.weight(0.0) - 1.0).abs() < 1e-15,
            "Weight at r=0 should be 1"
        );
        let w1 = nl.weight(1.0);
        let expected = (-0.5_f64).exp();
        assert!((w1 - expected).abs() < 1e-10);
    }
    #[test]
    fn test_nonlocal_bell_weight() {
        let nl = NonLocalRegularization::new(2.0);
        assert!((nl.weight_bell(0.0) - 1.0).abs() < 1e-15);
        assert_eq!(nl.weight_bell(2.0), 0.0, "Bell weight should be 0 at r=lc");
        assert_eq!(
            nl.weight_bell(3.0),
            0.0,
            "Bell weight should be 0 beyond lc"
        );
    }
    #[test]
    fn test_nonlocal_average_uniform_field() {
        let nl = NonLocalRegularization::new(1.0);
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let values = [5.0, 5.0, 5.0];
        let avg = nl.nonlocal_average(&[0.5, 0.0, 0.0], &positions, &values);
        assert!(
            (avg - 5.0).abs() < 1e-10,
            "Non-local avg of uniform field should be the field value"
        );
    }
    #[test]
    fn test_nonlocal_average_localized() {
        let nl = NonLocalRegularization::new(0.1);
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let values = [1.0, 0.0, 0.0];
        let avg = nl.nonlocal_average(&[0.0, 0.0, 0.0], &positions, &values);
        assert!(
            avg > 0.5,
            "Average at damaged point should be high, got {avg}"
        );
    }
    #[test]
    fn test_nonlocal_average_all() {
        let nl = NonLocalRegularization::new(1.0);
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let values = [2.0, 4.0];
        let avgs = nl.nonlocal_average_all(&positions, &values);
        assert_eq!(avgs.len(), 2);
        assert!(avgs[0] >= 2.0 - 1e-10 && avgs[0] <= 4.0 + 1e-10);
        assert!(avgs[1] >= 2.0 - 1e-10 && avgs[1] <= 4.0 + 1e-10);
    }
    #[test]
    fn test_gradient_enhanced_uniform_field() {
        let ged = GradientEnhancedDamage::new(0.1);
        let local = [0.001, 0.001, 0.001, 0.001, 0.001];
        let nonlocal = ged.solve_1d(&local, 0.01);
        assert_eq!(nonlocal.len(), 5);
        for (i, &val) in nonlocal.iter().enumerate() {
            assert!(
                (val - 0.001).abs() < 1e-8,
                "Non-local[{i}] = {val} should be 0.001"
            );
        }
    }
    #[test]
    fn test_gradient_enhanced_localized_peak() {
        let ged = GradientEnhancedDamage::new(0.5);
        let local = [0.0, 0.0, 1.0, 0.0, 0.0];
        let nonlocal = ged.solve_1d(&local, 0.1);
        assert!(
            nonlocal[2] < 1.0,
            "Peak should be reduced by regularization, got {}",
            nonlocal[2]
        );
        assert!(
            nonlocal[1] > 0.0,
            "Neighbors should pick up some strain, got {}",
            nonlocal[1]
        );
    }
    #[test]
    fn test_gradient_enhanced_single_element() {
        let ged = GradientEnhancedDamage::new(0.1);
        let local = [0.005];
        let nonlocal = ged.solve_1d(&local, 0.01);
        assert_eq!(nonlocal.len(), 1);
        assert!((nonlocal[0] - 0.005).abs() < 1e-15);
    }
    #[test]
    fn test_miner_fatigue_life() {
        let miner = MinerFatigueDamage::new(1.0e20, 5.0);
        let n = miner.fatigue_life(100.0e6);
        let expected = 1.0e20 / (100.0e6_f64).powf(5.0);
        assert!((n - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_miner_accumulation() {
        let mut miner = MinerFatigueDamage::new(1.0e20, 5.0);
        let n_f = miner.fatigue_life(100.0e6);
        miner.accumulate(100.0e6, n_f * 0.5);
        assert!(
            (miner.d - 0.5).abs() < 1e-6,
            "D should be 0.5, got {}",
            miner.d
        );
    }
    #[test]
    fn test_miner_failure() {
        let mut miner = MinerFatigueDamage::new(1.0e20, 5.0);
        let n_f = miner.fatigue_life(200.0e6);
        miner.accumulate(200.0e6, n_f * 1.5);
        assert!(
            miner.is_failed(),
            "Should be failed after exceeding fatigue life"
        );
        assert!(miner.d <= 1.0, "D should be capped at 1.0");
    }
    #[test]
    fn test_miner_multi_level() {
        let mut miner = MinerFatigueDamage::new(1.0e20, 5.0);
        let n1 = miner.fatigue_life(100.0e6);
        let n2 = miner.fatigue_life(200.0e6);
        miner.accumulate(100.0e6, n1 * 0.3);
        miner.accumulate(200.0e6, n2 * 0.2);
        assert!(
            (miner.d - 0.5).abs() < 1e-6,
            "D should be 0.5, got {}",
            miner.d
        );
    }
    #[test]
    fn test_progressive_damage_stiffness_undamaged() {
        let pd = ProgressiveDamage::new(200.0e9, 0.3, 0.0, 0.0);
        let c = pd.effective_stiffness_tensor();
        assert!(c[0][0] > 0.0, "C11 should be positive, got {}", c[0][0]);
    }
    #[test]
    fn test_progressive_damage_stiffness_reduction() {
        let pd_intact = ProgressiveDamage::new(200.0e9, 0.3, 0.0, 0.0);
        let pd_damaged = ProgressiveDamage::new(200.0e9, 0.3, 0.5, 0.5);
        let c0 = pd_intact.effective_stiffness_tensor();
        let cd = pd_damaged.effective_stiffness_tensor();
        assert!(cd[0][0] < c0[0][0], "damaged stiffness should be smaller");
    }
    #[test]
    fn test_progressive_damage_update_irreversible() {
        let mut pd = ProgressiveDamage::new(200.0e9, 0.3, 0.0, 0.0);
        pd.update_damage(0.3, 0.2);
        let d11 = pd.d11;
        let d22 = pd.d22;
        pd.update_damage(0.1, 0.05);
        assert!((pd.d11 - d11).abs() < 1e-15, "D11 should not decrease");
        assert!((pd.d22 - d22).abs() < 1e-15, "D22 should not decrease");
    }
    #[test]
    fn test_progressive_damage_residual_strength() {
        let pd = ProgressiveDamage::new(200.0e9, 0.3, 0.0, 0.0);
        let rs = pd.residual_strength(500.0e6, 0.0);
        assert!(
            (rs - 500.0e6).abs() < 1.0,
            "undamaged residual strength = original"
        );
        let pd_half = ProgressiveDamage::new(200.0e9, 0.3, 0.5, 0.0);
        let rs_half = pd_half.residual_strength(500.0e6, 0.0);
        assert!(rs_half < rs, "damaged residual strength should be lower");
    }
    #[test]
    fn test_puck_fiber_failure_tension() {
        let puck = PuckFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 200.0e6, 70.0e6);
        let ff = puck.fiber_failure_tension(1500.0e6, 130.0e9, 0.0, 10.0e9);
        assert!(
            (0.99..=1.01).contains(&ff),
            "FF at tensile strength should be ~1, got {ff}"
        );
    }
    #[test]
    fn test_puck_fiber_failure_compression() {
        let puck = PuckFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 200.0e6, 70.0e6);
        let ffc = puck.fiber_failure_compression(-1500.0e6);
        assert!(
            (0.99..=1.01).contains(&ffc),
            "FFC at compressive strength should be ~1, got {ffc}"
        );
    }
    #[test]
    fn test_puck_inter_fiber_failure_tension() {
        let puck = PuckFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 200.0e6, 70.0e6);
        let iff = puck.inter_fiber_failure_mode_a(50.0e6, 0.0, 0.25, 0.2);
        assert!(iff > 0.0, "IFF should be positive, got {iff}");
    }
    #[test]
    fn test_puck_no_failure_below_threshold() {
        let puck = PuckFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 200.0e6, 70.0e6);
        let ff = puck.fiber_failure_tension(100.0e6, 130.0e9, 0.0, 10.0e9);
        assert!(ff < 1.0, "no fiber failure below strength");
    }
    #[test]
    fn test_hashin_fiber_tension_at_strength() {
        let hashin = HashinFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        let f = hashin.fiber_tension(1500.0e6, 0.0, 70.0e6);
        assert!(
            (0.99..=1.01).contains(&f),
            "Fiber tension criterion at X_T should be ~1, got {f}"
        );
    }
    #[test]
    fn test_hashin_fiber_compression_at_strength() {
        let hashin = HashinFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        let f = hashin.fiber_compression(-1200.0e6);
        assert!(
            (0.99..=1.01).contains(&f),
            "Fiber compression criterion at X_C should be ~1, got {f}"
        );
    }
    #[test]
    fn test_hashin_matrix_tension_at_strength() {
        let hashin = HashinFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        let f = hashin.matrix_tension(50.0e6, 0.0);
        assert!(
            (0.99..=1.01).contains(&f),
            "Matrix tension at Y_T should be ~1, got {f}"
        );
    }
    #[test]
    fn test_hashin_matrix_compression() {
        let hashin = HashinFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        let f = hashin.matrix_compression(-250.0e6, 0.0);
        assert!(
            (0.99..=1.01).contains(&f),
            "Matrix compression at Y_C should be ~1, got {f}"
        );
    }
    #[test]
    fn test_hashin_no_failure_at_zero_stress() {
        let hashin = HashinFailureCriteria::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        assert!(hashin.fiber_tension(0.0, 0.0, 0.0) < 1e-15);
        assert!(hashin.fiber_compression(0.0) < 1e-15);
        assert!(hashin.matrix_tension(0.0, 0.0) < 1e-15);
        assert!(hashin.matrix_compression(0.0, 0.0) < 1e-15);
    }
    #[test]
    fn test_tsai_wu_at_x_t() {
        let tw = TsaiWuCriterion::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6, -0.5);
        let f = tw.failure_index(1500.0e6, 0.0, 0.0);
        assert!(f > 0.8, "Tsai-Wu index at X_T should be near 1, got {f}");
    }
    #[test]
    fn test_tsai_wu_zero_stress() {
        let tw = TsaiWuCriterion::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6, -0.5);
        let f = tw.failure_index(0.0, 0.0, 0.0);
        assert!(f.abs() < 1e-15, "zero stress → zero failure index, got {f}");
    }
    #[test]
    fn test_tsai_wu_compression_direction() {
        let tw = TsaiWuCriterion::new(1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6, -0.5);
        let ft = tw.failure_index(500.0e6, 0.0, 0.0);
        let fc = tw.failure_index(-500.0e6, 0.0, 0.0);
        assert_ne!(
            ft, fc,
            "tension/compression should give different Tsai-Wu index"
        );
    }
    #[test]
    fn test_pfa_ply_count() {
        let pfa = ProgressiveFailureAnalysis::new(4, 1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        assert_eq!(pfa.n_plies, 4);
        assert_eq!(pfa.ply_damage.len(), 4);
    }
    #[test]
    fn test_pfa_initial_no_failure() {
        let pfa = ProgressiveFailureAnalysis::new(4, 1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        assert!(
            !pfa.is_fully_failed(),
            "fresh laminate should not be fully failed"
        );
        for &d in &pfa.ply_damage {
            assert!(d < 1e-15, "initial damage should be zero");
        }
    }
    #[test]
    fn test_pfa_degradation_increases_damage() {
        let mut pfa =
            ProgressiveFailureAnalysis::new(3, 1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        pfa.degrade_ply(0, 0.5);
        assert!(
            (pfa.ply_damage[0] - 0.5).abs() < 1e-15,
            "ply damage should be 0.5"
        );
    }
    #[test]
    fn test_pfa_full_laminate_failure() {
        let mut pfa =
            ProgressiveFailureAnalysis::new(2, 1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        pfa.degrade_ply(0, 1.0);
        pfa.degrade_ply(1, 1.0);
        assert!(pfa.is_fully_failed(), "all plies failed → laminate failed");
    }
    #[test]
    fn test_pfa_mean_damage() {
        let mut pfa =
            ProgressiveFailureAnalysis::new(4, 1500.0e6, 1200.0e6, 50.0e6, 250.0e6, 70.0e6);
        pfa.degrade_ply(0, 0.2);
        pfa.degrade_ply(1, 0.4);
        pfa.degrade_ply(2, 0.6);
        pfa.degrade_ply(3, 0.8);
        let mean = pfa.mean_damage();
        assert!((mean - 0.5).abs() < 1e-14, "mean damage = 0.5, got {mean}");
    }
    #[test]
    fn test_anisotropic_damage_initial_isotropic() {
        let dia = DamageInducedAnisotropy::new(200.0e9, 0.3);
        assert_eq!(dia.d12, 0.0);
        assert_eq!(dia.d13, 0.0);
        assert_eq!(dia.d23, 0.0);
    }
    #[test]
    fn test_anisotropic_damage_shear_reduction() {
        let mut dia = DamageInducedAnisotropy::new(200.0e9, 0.3);
        dia.update_shear_damage(0.5, 0.3, 0.2);
        let g_eff = dia.effective_shear_modulus_12();
        let g0 = 200.0e9 / (2.0 * 1.3);
        assert!(g_eff < g0, "shear modulus should decrease with damage");
    }
    #[test]
    fn test_repair_effectiveness_full_repair() {
        let rep = RepairEffectiveness::new(1.0, 0.0);
        assert!((rep.effectiveness() - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_repair_effectiveness_no_repair() {
        let rep = RepairEffectiveness::new(0.0, 0.8);
        assert!(
            rep.effectiveness() < 0.5,
            "no repair should have low effectiveness"
        );
    }
    #[test]
    fn test_repair_effectiveness_partial() {
        let rep = RepairEffectiveness::new(0.6, 0.4);
        let eff = rep.effectiveness();
        assert!(
            eff > 0.0 && eff < 1.0,
            "partial repair effectiveness in (0,1), got {eff}"
        );
    }
    #[test]
    fn test_repair_strength_recovery() {
        let rep = RepairEffectiveness::new(0.8, 0.2);
        let recovered = rep.strength_recovery(500.0e6);
        assert!(
            recovered > 0.0 && recovered <= 500.0e6,
            "recovered = {recovered}"
        );
    }
    #[test]
    fn test_triaxiality_correction_uniaxial() {
        let lem = LemaitreDuctileDamage::new(0.01, 1.0e6, 2.0, 0.3, 0.3);
        let tc = lem.compute_triaxiality_correction(1.0 / 3.0);
        assert!(
            tc > 0.0,
            "triaxiality correction should be positive, got {tc}"
        );
    }
    #[test]
    fn test_triaxiality_correction_increases_with_triaxiality() {
        let lem = LemaitreDuctileDamage::new(0.01, 1.0e6, 2.0, 0.3, 0.3);
        let tc_low = lem.compute_triaxiality_correction(0.33);
        let tc_high = lem.compute_triaxiality_correction(1.0);
        assert!(
            tc_high > tc_low,
            "higher triaxiality → larger correction: {tc_low} vs {tc_high}"
        );
    }
    #[test]
    fn test_triaxiality_correction_zero_triaxiality() {
        let lem = LemaitreDuctileDamage::new(0.01, 1.0e6, 2.0, 0.3, 0.3);
        let tc = lem.compute_triaxiality_correction(0.0);
        let expected = (2.0 / 3.0) * (1.0 + 0.3);
        assert!(
            (tc - expected).abs() < 1e-12,
            "tc={tc}, expected={expected}"
        );
    }
    #[test]
    fn test_gurson_void_growth_rate_positive() {
        let gtn = GursonDamage::new(0.01, 0.15, 0.25, 0.04, 0.3, 0.1);
        let rate = gtn.compute_void_growth_rate(200.0e6, 500.0e6, 300.0e6, 0.01);
        assert!(
            rate >= 0.0,
            "void growth rate should be non-negative, got {rate}"
        );
    }
    #[test]
    fn test_gurson_void_growth_rate_increases_with_triaxiality() {
        let gtn = GursonDamage::new(0.01, 0.15, 0.25, 0.04, 0.3, 0.1);
        let rate_low = gtn.compute_void_growth_rate(100.0e6, 500.0e6, 100.0e6, 0.01);
        let rate_high = gtn.compute_void_growth_rate(300.0e6, 500.0e6, 100.0e6, 0.01);
        assert!(
            rate_high >= rate_low,
            "higher triaxiality → faster void growth: {rate_low} vs {rate_high}"
        );
    }
    #[test]
    fn test_gurson_void_growth_rate_zero_plastic_strain() {
        let gtn = GursonDamage::new(0.01, 0.15, 0.25, 0.04, 0.3, 0.1);
        let rate = gtn.compute_void_growth_rate(200.0e6, 500.0e6, 300.0e6, 0.0);
        assert!(
            rate.abs() < 1e-30,
            "zero plastic increment → zero growth, got {rate}"
        );
    }
    #[test]
    fn test_localization_indicator_undamaged() {
        let cdm = CdmDamage::new(0.001, 0.01, 5.0, DamageEvolutionLaw::Exponential);
        let indicator = cdm.compute_localization_indicator(200.0e9, 0.3);
        assert!(
            indicator >= 0.0,
            "localization indicator should be non-negative, got {indicator}"
        );
    }
    #[test]
    fn test_localization_indicator_increases_with_damage() {
        let mut cdm_lo = CdmDamage::new(0.001, 0.01, 5.0, DamageEvolutionLaw::Exponential);
        let mut cdm_hi = CdmDamage::new(0.001, 0.01, 5.0, DamageEvolutionLaw::Exponential);
        cdm_lo.update(0.003);
        cdm_hi.update(0.008);
        let ind_lo = cdm_lo.compute_localization_indicator(200.0e9, 0.3);
        let ind_hi = cdm_hi.compute_localization_indicator(200.0e9, 0.3);
        assert!(
            ind_hi >= ind_lo,
            "more damage → higher localization indicator: {ind_lo} vs {ind_hi}"
        );
    }
    #[test]
    fn test_localization_indicator_fully_damaged() {
        let mut cdm = CdmDamage::new(0.0, 0.001, 50.0, DamageEvolutionLaw::Exponential);
        cdm.update(1.0);
        let indicator = cdm.compute_localization_indicator(200.0e9, 0.3);
        assert!(
            indicator > 1.0,
            "near-full damage → high localization indicator: {indicator}"
        );
    }
}
