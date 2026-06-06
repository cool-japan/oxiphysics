//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Von-Mises equivalent stress from a Voigt stress vector
/// \[sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz\].
pub(super) fn von_mises(s: &[f64; 6]) -> f64 {
    let sxx = s[0];
    let syy = s[1];
    let szz = s[2];
    let txy = s[3];
    let tyz = s[4];
    let txz = s[5];
    let val = 0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2))
        + 3.0 * (txy * txy + tyz * tyz + txz * txz);
    val.sqrt()
}
/// Hydrostatic (mean) stress sigma_H = (sigma_xx + sigma_yy + sigma_zz) / 3.
pub(super) fn hydrostatic(s: &[f64; 6]) -> f64 {
    (s[0] + s[1] + s[2]) / 3.0
}
/// Von-Mises equivalent stress (module-local version to avoid conflicts).
pub(super) fn von_mises_local(s: &[f64; 6]) -> f64 {
    let sxx = s[0];
    let syy = s[1];
    let szz = s[2];
    let txy = s[3];
    let tyz = s[4];
    let txz = s[5];
    let val = 0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2))
        + 3.0 * (txy * txy + tyz * tyz + txz * txz);
    val.sqrt()
}
/// Deviatoric stress tensor.
pub(super) fn deviatoric_stress(s: &[f64; 6]) -> [f64; 6] {
    let mean = (s[0] + s[1] + s[2]) / 3.0;
    [s[0] - mean, s[1] - mean, s[2] - mean, s[3], s[4], s[5]]
}
#[cfg(test)]
mod tests {

    use crate::damage::*;
    fn iso_mat() -> IsotropicDamage {
        IsotropicDamage {
            threshold: 1e-4,
            a: 5000.0,
        }
    }
    #[test]
    fn test_isotropic_damage_below_threshold() {
        let mat = iso_mat();
        let strain = [5e-5, 0.0, 0.0, 0.0, 0.0, 0.0];
        let kappa = IsotropicDamage::update_kappa(0.0, &strain);
        let d = mat.damage_variable(kappa);
        assert_eq!(d, 0.0, "expected D = 0 below threshold, got {d}");
    }
    #[test]
    fn test_isotropic_damage_above_threshold() {
        let mat = iso_mat();
        let strain = [5e-4, 0.0, 0.0, 0.0, 0.0, 0.0];
        let kappa = IsotropicDamage::update_kappa(0.0, &strain);
        let d = mat.damage_variable(kappa);
        assert!(d > 0.0, "expected D > 0 above threshold, got {d}");
        assert!(d < 1.0, "damage should not be fully saturated yet, got {d}");
    }
    #[test]
    fn test_mazars_tension_positive() {
        let maz = MazarsDamage {
            ac: 1.4,
            bc: 1900.0,
            at: 0.8,
            bt: 20_000.0,
            kappa0: 1e-4,
        };
        let kappa = 5e-4;
        let dt = maz.damage_tension(kappa);
        assert!(dt > 0.0, "expected D_t > 0 for kappa > kappa_0, got {dt}");
    }
    #[test]
    fn test_lemaitre_thermodynamic_force_positive() {
        let stress = [200e6_f64, 0.0, 0.0, 0.0, 0.0, 0.0];
        let y = LemaitreDamage::thermodynamic_force(&stress, 0.0, 200e9, 0.3);
        assert!(y > 0.0, "thermodynamic force should be positive, got {y}");
    }
    #[test]
    fn test_effective_stress_larger_than_nominal() {
        let stress = [100e6_f64, 0.0, 0.0, 0.0, 0.0, 0.0];
        let d = 0.3;
        let eff = LemaitreDamage::effective_stress(&stress, d);
        assert!(
            eff[0] > stress[0],
            "effective stress ({}) should exceed nominal ({}) when D = {}",
            eff[0],
            stress[0],
            d
        );
    }
    #[test]
    fn test_element_max_damage_increases() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(4, mat);
        assert_eq!(elem.max_damage(), 0.0);
        let high_strain = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let strains = vec![high_strain; 4];
        elem.update_all(&strains);
        assert!(
            elem.max_damage() > 0.0,
            "max_damage should be > 0 after high strain, got {}",
            elem.max_damage()
        );
    }
    #[test]
    fn test_failure_mode_intact_for_zero_damage() {
        let mode = FailureMode::classify(0.0, 0.5);
        assert_eq!(
            mode,
            FailureMode::Intact,
            "expected Intact for D = 0, got {mode:?}"
        );
    }
    #[test]
    fn test_damage_rate_limiter() {
        let limiter = DamageRateLimiter::new(0.01);
        assert!((limiter.limit(0.005) - 0.005).abs() < 1e-10);
        assert!((limiter.limit(0.05) - 0.01).abs() < 1e-10);
        assert!((limiter.limit(-0.01) - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_damage_rate_limiter_apply() {
        let limiter = DamageRateLimiter::new(0.02);
        let mut state = DamageState {
            d: 0.1,
            ..Default::default()
        };
        limiter.apply(&mut state, 0.5);
        assert!((state.d - 0.12).abs() < 1e-10);
        assert!((state.damage_rate - 0.02).abs() < 1e-10);
    }
    #[test]
    fn test_element_deletion() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(4, mat);
        for state in &mut elem.damage_states {
            state.d = 0.99;
        }
        assert!(!elem.deleted);
        elem.check_deletion(0.95);
        assert!(elem.deleted);
    }
    #[test]
    fn test_element_partial_deletion() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(4, mat);
        elem.damage_states[0].d = 0.99;
        elem.damage_states[1].d = 0.99;
        elem.damage_states[2].d = 0.1;
        elem.damage_states[3].d = 0.1;
        elem.delete_failed_points(0.95);
        assert_eq!(elem.deleted_point_count(), 2);
        assert_eq!(elem.active_point_count(), 2);
        assert!(!elem.deleted);
    }
    #[test]
    fn test_element_all_points_deleted() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(2, mat);
        for state in &mut elem.damage_states {
            state.d = 0.99;
        }
        elem.delete_failed_points(0.95);
        assert!(elem.deleted);
        assert_eq!(elem.deleted_point_count(), 2);
    }
    #[test]
    fn test_damage_visualization_colors() {
        let color_0 = DamageVisualization::damage_to_color(0.0);
        assert!((color_0[0] - 0.0).abs() < 1e-10);
        assert!((color_0[2] - 1.0).abs() < 1e-10);
        let color_1 = DamageVisualization::damage_to_color(1.0);
        assert!((color_1[0] - 1.0).abs() < 1e-10);
        assert!((color_1[1] - 0.0).abs() < 1e-10);
        assert!((color_1[2] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_damage_visualization_from_states() {
        let states = vec![
            DamageState {
                d: 0.0,
                ..DamageState::default()
            },
            DamageState {
                d: 0.5,
                ..DamageState::default()
            },
            DamageState {
                d: 1.0,
                ..DamageState::default()
            },
        ];
        let vis = DamageVisualization::from_states(&states);
        assert_eq!(vis.damage_values.len(), 3);
        assert_eq!(vis.colors.len(), 3);
        assert!((vis.min_damage() - 0.0).abs() < 1e-10);
        assert!((vis.max_damage() - 1.0).abs() < 1e-10);
        assert!((vis.avg_damage() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_nonlocal_damage_weight() {
        let nl = NonLocalDamage::new(1.0);
        assert!((nl.weight(0.0) - 1.0).abs() < 1e-10);
        assert!(nl.weight(10.0) < 1e-10);
    }
    #[test]
    fn test_nonlocal_strain_averaging() {
        let nl = NonLocalDamage::new(1.0);
        let local_strains = vec![0.1, 0.2, 0.3];
        let distances = vec![0.0, 0.5, 1.0];
        let nl_strain = nl.nonlocal_strain(&local_strains, &distances);
        assert!(nl_strain > 0.1);
        assert!(nl_strain < 0.3);
    }
    #[test]
    fn test_coupled_strain() {
        let nl = NonLocalDamage::new(1.0);
        let c1 = nl.coupled_strain(0.5, 0.3, 1.0);
        assert!((c1 - 0.5).abs() < 1e-10);
        let c0 = nl.coupled_strain(0.5, 0.3, 0.0);
        assert!((c0 - 0.3).abs() < 1e-10);
        let c05 = nl.coupled_strain(0.5, 0.3, 0.5);
        assert!((c05 - 0.4).abs() < 1e-10);
    }
    #[test]
    fn test_lemaitre_rate_limiting() {
        let lem = LemaitreDamage {
            s: 2.0,
            r: 1.0,
            d_c: 0.9,
            yield_stress: 250e6,
            hardening: 1e9,
        };
        let limiter = DamageRateLimiter::new(0.001);
        let mut state = DamageState::default();
        let stress = [300e6_f64, 0.0, 0.0, 0.0, 0.0, 0.0];
        lem.update_damage_limited(&mut state, &stress, 0.01, 200e9, 0.3, &limiter);
        assert!(state.d <= 0.001 + 1e-10);
    }
    #[test]
    fn test_element_update_skips_deleted() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(4, mat);
        elem.damage_states[0].deleted = true;
        let strains = vec![[1e-3, 0.0, 0.0, 0.0, 0.0, 0.0]; 4];
        elem.update_all(&strains);
        assert!((elem.damage_states[0].d - 0.0).abs() < 1e-10);
        assert!(elem.damage_states[1].d > 0.0);
    }
    #[test]
    fn test_element_avg_damage() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(4, mat);
        elem.damage_states[0].d = 0.2;
        elem.damage_states[1].d = 0.4;
        elem.damage_states[2].d = 0.6;
        elem.damage_states[3].d = 0.8;
        let avg = elem.avg_damage();
        assert!((avg - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_element_stiffness_reduction_factors() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(3, mat);
        elem.damage_states[0].d = 0.0;
        elem.damage_states[1].d = 0.5;
        elem.damage_states[2].d = 1.0;
        elem.damage_states[2].deleted = true;
        let factors = elem.stiffness_reduction_factors();
        assert!((factors[0] - 1.0).abs() < 1e-10);
        assert!((factors[1] - 0.5).abs() < 1e-10);
        assert!((factors[2] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_element_deletion_manager() {
        let mat = iso_mat();
        let mut elems = vec![
            DamageMechanicsElement::with_id(2, mat.clone(), 0),
            DamageMechanicsElement::with_id(2, mat.clone(), 1),
            DamageMechanicsElement::with_id(2, mat, 2),
        ];
        for s in &mut elems[0].damage_states {
            s.d = 0.99;
        }
        elems[1].damage_states[0].d = 0.99;
        elems[1].damage_states[1].d = 0.1;
        let mut mgr = ElementDeletionManager::new(0.95);
        mgr.process(&mut elems);
        assert_eq!(mgr.deleted_count, 1);
        assert_eq!(mgr.deleted_ids, vec![0]);
        assert!(elems[0].deleted);
        assert!(!elems[1].deleted);
        assert!(!elems[2].deleted);
        assert!((mgr.deletion_ratio(3) - 1.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_failure_mode_helpers() {
        let intact = FailureMode::classify(0.0, 0.5);
        assert!(!intact.is_damaged());
        assert!(!intact.is_fully_damaged());
        assert!((intact.damage_value() - 0.0).abs() < 1e-10);
        let tension = FailureMode::classify(0.5, 0.7);
        assert!(tension.is_damaged());
        assert!(!tension.is_fully_damaged());
        assert!((tension.damage_value() - 0.5).abs() < 1e-10);
        let fully = FailureMode::classify(0.99, 0.5);
        assert!(fully.is_damaged());
        assert!(fully.is_fully_damaged());
        assert!((fully.damage_value() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_isotropic_stiffness_reduction() {
        let mat = iso_mat();
        assert!((mat.stiffness_reduction_factor(5e-5) - 1.0).abs() < 1e-10);
        let factor = mat.stiffness_reduction_factor(5e-4);
        assert!(factor < 1.0);
        assert!(factor > 0.0);
    }
    #[test]
    fn test_mazars_update_rate_limited() {
        let maz = MazarsDamage {
            ac: 1.4,
            bc: 1900.0,
            at: 0.8,
            bt: 20_000.0,
            kappa0: 1e-4,
        };
        let limiter = DamageRateLimiter::new(0.001);
        let mut state = DamageState::default();
        let strain = [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0];
        maz.update_limited(&mut state, &strain, 0.5, &limiter);
        assert!(state.d <= 0.001 + 1e-10);
    }
    #[test]
    fn test_lemaitre_damaged_stiffness() {
        let mut c = [[0.0_f64; 6]; 6];
        for (i, row) in c.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        let cd = LemaitreDamage::damaged_stiffness(0.5, &c);
        for (i, row) in cd.iter().enumerate() {
            assert!((row[i] - 0.5).abs() < 1e-10);
        }
    }
    #[test]
    fn test_element_update_limited() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(2, mat);
        let limiter = DamageRateLimiter::new(0.001);
        let strains = vec![[1e-2, 0.0, 0.0, 0.0, 0.0, 0.0]; 2];
        elem.update_all_limited(&strains, &limiter);
        for state in &elem.damage_states {
            assert!(state.d <= 0.001 + 1e-10);
        }
    }
    #[test]
    fn test_visualization_empty() {
        let vis = DamageVisualization::from_states(&[]);
        assert_eq!(vis.damage_values.len(), 0);
        assert!((vis.avg_damage() - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_damage_state_default_not_deleted() {
        let state = DamageState::default();
        assert!(!state.deleted);
        assert!((state.damage_rate - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_element_visualization() {
        let mat = iso_mat();
        let mut elem = DamageMechanicsElement::new(3, mat);
        elem.damage_states[0].d = 0.0;
        elem.damage_states[1].d = 0.5;
        elem.damage_states[2].d = 1.0;
        let vis = elem.visualization();
        assert_eq!(vis.colors.len(), 3);
        assert!((vis.colors[0][2] - 1.0).abs() < 1e-10);
        assert!((vis.colors[2][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_scalar_isotropic_damage_monotone() {
        let mut model = ScalarIsotropicDamage {
            d: 0.0,
            kappa: 0.0,
            kappa_0: 1e-4,
        };
        let mut prev_d = 0.0f64;
        for i in 1..=5 {
            let eps = vec![i as f64 * 5e-4, 0.0, 0.0];
            let eps_eq = ScalarIsotropicDamage::equivalent_strain(&eps);
            let d = model.update(eps_eq);
            assert!(
                d >= prev_d,
                "D should be non-decreasing, got {d} < {prev_d}"
            );
            prev_d = d;
        }
        assert!(prev_d > 0.0, "D should be positive after loading");
    }
    #[test]
    fn test_scalar_isotropic_damage_zero_below_threshold() {
        let mut model = ScalarIsotropicDamage {
            d: 0.0,
            kappa: 0.0,
            kappa_0: 1e-3,
        };
        let eps = vec![1e-5, 0.0, 0.0];
        let eps_eq = ScalarIsotropicDamage::equivalent_strain(&eps);
        let d = model.update(eps_eq);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_gurson_yield_positive_outside() {
        let model = GursonModel {
            f: 0.01,
            q1: 1.5,
            q2: 1.0,
            q3: 2.25,
        };
        let phi = model.phi(400e6, 100e6, 200e6);
        assert!(
            phi > 0.0,
            "yield function should be positive outside surface, got {phi}"
        );
    }
    #[test]
    fn test_gurson_yield_zero_on_surface() {
        let model = GursonModel {
            f: 0.0,
            q1: 1.5,
            q2: 1.0,
            q3: 2.25,
        };
        let phi = model.phi(200e6, 0.0, 200e6);
        assert!(
            phi.abs() < 1e-6,
            "yield function on surface should be ~0, got {phi}"
        );
    }
    #[test]
    fn test_void_growth_rate_positive() {
        let model = GursonModel {
            f: 0.01,
            q1: 1.5,
            q2: 1.0,
            q3: 2.25,
        };
        let rate = model.void_growth_rate(0.01);
        assert!(
            rate > 0.0,
            "void growth rate must be positive for positive eps_p_vol"
        );
    }
    #[test]
    fn test_lemaitre_cdm_rate_positive() {
        let model = LemaitreCDM {
            d: 0.0,
            s: 1.0,
            s0: 1.0,
        };
        let rate = model.evolution_rate(0.01, 300e6, 100e6);
        assert!(
            rate >= 0.0,
            "damage rate should be non-negative, got {rate}"
        );
    }
    #[test]
    fn test_crack_band_softening_slope_negative() {
        let model = CrackBandModel {
            gf: 100.0,
            fc: 3e6,
            h: 0.02,
        };
        let slope = model.softening_slope();
        assert!(slope < 0.0, "softening slope must be negative, got {slope}");
    }
    #[test]
    fn test_crack_band_softening_formula() {
        let model = CrackBandModel {
            gf: 100.0,
            fc: 3e6,
            h: 0.02,
        };
        let expected = -(3e6_f64 * 3e6 * 0.02) / (2.0 * 100.0);
        assert!((model.softening_slope() - expected).abs() < 1.0);
    }
    #[test]
    fn test_fatigue_damage_rate_zero_below_endurance() {
        let model = FatigueDamageModel::new(400e6, 2.0, 200e6, 0.0, 0.1);
        let rate = model.damage_rate_per_cycle(100e6, 0.0, 0.1);
        assert_eq!(rate, 0.0, "No damage below endurance limit");
    }
    #[test]
    fn test_fatigue_damage_rate_positive_above_endurance() {
        let model = FatigueDamageModel::new(400e6, 2.0, 200e6, 0.0, 0.1);
        let rate = model.damage_rate_per_cycle(300e6, 0.0, 0.1);
        assert!(
            rate > 0.0,
            "Damage rate should be positive above endurance, got {rate}"
        );
    }
    #[test]
    fn test_fatigue_damage_cycles_to_failure_positive() {
        let model = FatigueDamageModel::new(400e6, 2.0, 200e6, 0.0, 0.1);
        let nf = model.cycles_to_failure(300e6, 0.0);
        assert!(
            nf > 0.0 && nf.is_finite(),
            "Cycles to failure should be positive finite, got {nf}"
        );
    }
    #[test]
    fn test_fatigue_higher_stress_fewer_cycles() {
        let model = FatigueDamageModel::new(400e6, 2.0, 200e6, 0.0, 0.1);
        let n1 = model.cycles_to_failure(250e6, 0.0);
        let n2 = model.cycles_to_failure(350e6, 0.0);
        assert!(n2 < n1, "Higher stress should give fewer cycles to failure");
    }
    #[test]
    fn test_fatigue_remaining_life_less_than_total() {
        let model = FatigueDamageModel::new(400e6, 2.0, 200e6, 0.0, 0.1);
        let nf_total = model.cycles_to_failure(300e6, 0.0);
        let nf_remain = model.remaining_life(300e6, 0.0, 0.3, 100);
        assert!(
            nf_remain < nf_total + 1.0,
            "Remaining life from D=0.3 < total life"
        );
    }
    #[test]
    fn test_coupled_damage_plasticity_yield_stress_increases() {
        let model = CoupledDamagePlasticity::new(CoupledDamagePlasticityParams {
            e: 200e9,
            nu: 0.3,
            sigma_y0: 250e6,
            h: 1e9,
            q: 100e6,
            b_hard: 5.0,
            c_kin: 1e9,
            gamma_kin: 0.5,
            s_dmg: 2.0,
            s_dmg_exp: 1.0,
            d_c: 0.8,
        });
        let s0 = model.yield_stress(0.0);
        let s1 = model.yield_stress(0.01);
        assert!(s1 > s0, "Yield stress should increase with plastic strain");
    }
    #[test]
    fn test_coupled_damage_plasticity_elastic_step() {
        let model = CoupledDamagePlasticity::new(CoupledDamagePlasticityParams {
            e: 200e9,
            nu: 0.3,
            sigma_y0: 250e6,
            h: 1e9,
            q: 100e6,
            b_hard: 5.0,
            c_kin: 1e9,
            gamma_kin: 0.5,
            s_dmg: 2.0,
            s_dmg_exp: 1.0,
            d_c: 0.8,
        });
        let state = DamagePlasticityState::default();
        let sigma_trial = [100e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (new_state, converged) = model.return_mapping(&state, &sigma_trial);
        assert!(converged, "Should converge");
        assert!(
            (new_state.d - 0.0).abs() < 1e-12,
            "No damage in elastic regime"
        );
        assert!(
            (new_state.p_bar - 0.0).abs() < 1e-12,
            "No plastic strain in elastic regime"
        );
    }
    #[test]
    fn test_coupled_damage_plasticity_plastic_step() {
        let model = CoupledDamagePlasticity::new(CoupledDamagePlasticityParams {
            e: 200e9,
            nu: 0.3,
            sigma_y0: 250e6,
            h: 1e9,
            q: 100e6,
            b_hard: 5.0,
            c_kin: 1e9,
            gamma_kin: 0.5,
            s_dmg: 2.0,
            s_dmg_exp: 1.0,
            d_c: 0.8,
        });
        let state = DamagePlasticityState::default();
        let sigma_trial = [600e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (new_state, converged) = model.return_mapping(&state, &sigma_trial);
        assert!(converged, "Should converge");
        assert!(new_state.p_bar > 0.0, "Should have plastic strain");
    }
    #[test]
    fn test_coupled_damage_plasticity_shear_modulus() {
        let model = CoupledDamagePlasticity::new(CoupledDamagePlasticityParams {
            e: 200e9,
            nu: 0.3,
            sigma_y0: 250e6,
            h: 1e9,
            q: 100e6,
            b_hard: 5.0,
            c_kin: 1e9,
            gamma_kin: 0.5,
            s_dmg: 2.0,
            s_dmg_exp: 1.0,
            d_c: 0.8,
        });
        let g = model.shear_modulus();
        let expected = 200e9 / (2.0 * 1.3);
        assert!((g - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_nonlocal_continuum_gradient_parameter() {
        let model = NonlocalContinuumDamage::new(0.02, 1e-4, 0.99, 200.0);
        let c = model.gradient_parameter();
        assert!((c - 0.02 * 0.02 / 2.0).abs() < 1e-20);
    }
    #[test]
    fn test_nonlocal_continuum_no_damage_below_threshold() {
        let model = NonlocalContinuumDamage::new(0.02, 1e-4, 0.99, 0.99);
        let d = model.damage_loading_function(5e-5, 300.0);
        assert_eq!(d, 0.0, "No damage below threshold");
    }
    #[test]
    fn test_nonlocal_continuum_damage_above_threshold() {
        let model = NonlocalContinuumDamage::new(0.02, 1e-4, 0.99, 0.5);
        let d = model.damage_loading_function(1e-3, 10.0);
        assert!(
            d > 0.0,
            "Damage should be positive above threshold, got {d}"
        );
        assert!(d < 1.0, "Damage should be less than 1.0, got {d}");
    }
    #[test]
    fn test_nonlocal_continuum_smoothing() {
        let model = NonlocalContinuumDamage::new(1.0, 1e-4, 0.99, 200.0);
        let positions = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let local_strains = vec![0.0, 0.0, 1.0, 0.0, 0.0];
        let nl_strains = model.compute_nonlocal_strains(&positions, &local_strains);
        assert!(
            nl_strains[2] <= local_strains[2] + 1e-10,
            "Nonlocal strain should not exceed local"
        );
        assert!(
            nl_strains[1] > 0.0,
            "Neighbor should have nonzero nonlocal strain"
        );
    }
    #[test]
    fn test_lemaitre_chaboche_energy_release_rate_positive() {
        let model = LemaitreChabocheDamage {
            s_strength: 2.0,
            s: 1.0,
            p_d: 0.0,
            d_c: 0.9,
            e: 200e9,
            nu: 0.3,
        };
        let sigma = [300e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let y = model.damage_energy_release_rate(&sigma, 0.0);
        assert!(
            y > 0.0,
            "Damage energy release rate should be positive, got {y}"
        );
    }
    #[test]
    fn test_lemaitre_chaboche_no_damage_below_threshold() {
        let model = LemaitreChabocheDamage {
            s_strength: 2.0,
            s: 1.0,
            p_d: 0.01,
            d_c: 0.9,
            e: 200e9,
            nu: 0.3,
        };
        let sigma = [300e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let dd = model.damage_increment(&sigma, 0.0, 0.001, 0.005);
        assert_eq!(dd, 0.0, "No damage below plastic threshold");
    }
    #[test]
    fn test_lemaitre_chaboche_damage_accumulates() {
        let model = LemaitreChabocheDamage {
            s_strength: 2.0,
            s: 1.0,
            p_d: 0.0,
            d_c: 0.9,
            e: 200e9,
            nu: 0.3,
        };
        let mut state = DamageState::default();
        let sigma = [500e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        for _ in 0..10 {
            model.update(&mut state, &sigma, 0.001);
        }
        assert!(
            state.d > 0.0,
            "Damage should accumulate after plastic loading"
        );
        assert!(
            state.accumulated_plastic_strain > 0.0,
            "Plastic strain should accumulate"
        );
    }
    #[test]
    fn test_lemaitre_chaboche_saturation_at_dc() {
        let model = LemaitreChabocheDamage {
            s_strength: 0.1,
            s: 1.0,
            p_d: 0.0,
            d_c: 0.8,
            e: 200e9,
            nu: 0.3,
        };
        let mut state = DamageState::default();
        let sigma = [1e9, 0.0, 0.0, 0.0, 0.0, 0.0];
        for _ in 0..200 {
            model.update(&mut state, &sigma, 0.01);
        }
        assert!(state.d <= model.d_c + 1e-10, "Damage should not exceed D_c");
    }
    #[test]
    fn test_lemaitre_chaboche_effective_stress_scales() {
        let model = LemaitreChabocheDamage {
            s_strength: 2.0,
            s: 1.0,
            p_d: 0.0,
            d_c: 0.9,
            e: 200e9,
            nu: 0.3,
        };
        let sigma = [100e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let d = 0.5;
        let eff = model.effective_stress(&sigma, d);
        assert!(
            (eff[0] - 200e6).abs() < 1.0,
            "Effective stress should be doubled when D=0.5"
        );
    }
    #[test]
    fn test_thermal_damage_at_reference_temp() {
        let model = ThermalDamage::new(300.0, 1500.0, 200e9, 2.0);
        let factor = model.stiffness_factor(300.0);
        assert!(
            (factor - 1.0).abs() < 1e-12,
            "No degradation at reference temperature"
        );
    }
    #[test]
    fn test_thermal_damage_at_melt_temp() {
        let model = ThermalDamage::new(300.0, 1500.0, 200e9, 2.0);
        let factor = model.stiffness_factor(1500.0);
        assert!(
            (factor - 0.0).abs() < 1e-12,
            "Full degradation at melting temperature"
        );
    }
    #[test]
    fn test_thermal_damage_monotone() {
        let model = ThermalDamage::new(300.0, 1500.0, 200e9, 2.0);
        let f1 = model.stiffness_factor(600.0);
        let f2 = model.stiffness_factor(900.0);
        assert!(f2 < f1, "Stiffness should decrease with temperature");
    }
    #[test]
    fn test_combined_thermal_mechanical_damage() {
        let model = ThermalDamage::new(300.0, 1500.0, 200e9, 2.0);
        let d_mech = 0.2;
        let temp = 900.0;
        let d_total = model.combined_damage(d_mech, temp);
        assert!(d_total > d_mech, "Combined damage > mechanical alone");
        assert!(d_total < 1.0, "Combined damage < 1");
    }
    #[test]
    fn test_localization_condition_no_softening() {
        let model = DamageBandLocalization::new(0.1, 0.0, 200e9, 0.3);
        assert!(!model.localization_condition());
    }
    #[test]
    fn test_localization_condition_softening() {
        let model = DamageBandLocalization::new(0.5, 300e9, 200e9, 0.3);
        assert!(model.localization_condition());
    }
    #[test]
    fn test_process_zone_width_positive() {
        let model = DamageBandLocalization::new(0.1, 1e9, 200e9, 0.3);
        let w = model.process_zone_width(1e-6);
        assert!(w > 0.0, "Process zone width should be positive, got {w}");
    }
    #[test]
    fn test_critical_wave_vector_none_when_no_localization() {
        let model = DamageBandLocalization::new(0.0, 0.0, 200e9, 0.3);
        assert!(model.critical_wave_vector().is_none());
    }
    #[test]
    fn test_damage_homogenization_zero_crack_density() {
        let model = DamageHomogenization::new(200e9, 0.3, 0.0);
        assert!(
            (model.effective_youngs_modulus() - 200e9).abs() < 1.0,
            "Zero crack density → E_eff = E_m"
        );
    }
    #[test]
    fn test_damage_homogenization_reduces_modulus() {
        let model = DamageHomogenization::new(200e9, 0.3, 0.5);
        assert!(
            model.effective_youngs_modulus() < 200e9,
            "Cracks should reduce Young's modulus"
        );
    }
    #[test]
    fn test_damage_from_crack_density_zero_at_no_cracks() {
        let model = DamageHomogenization::new(200e9, 0.3, 0.0);
        assert!(
            (model.damage_from_crack_density() - 0.0).abs() < 1e-12,
            "Zero crack density → D = 0"
        );
    }
    #[test]
    fn test_damage_from_crack_density_increases_with_density() {
        let m1 = DamageHomogenization::new(200e9, 0.3, 0.1);
        let m2 = DamageHomogenization::new(200e9, 0.3, 0.3);
        assert!(
            m2.damage_from_crack_density() > m1.damage_from_crack_density(),
            "Higher crack density → more damage"
        );
    }
    #[test]
    fn test_exponential_law_below_threshold_zero() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let d = law.compute_damage(5e-5);
        assert_eq!(d, 0.0, "Damage below threshold should be 0");
    }
    #[test]
    fn test_exponential_law_at_threshold_zero() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let d = law.compute_damage(1e-4);
        assert_eq!(d, 0.0, "Damage at threshold should be 0");
    }
    #[test]
    fn test_exponential_law_above_threshold_positive() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let d = law.compute_damage(5e-4);
        assert!(
            d > 0.0 && d < 1.0,
            "Damage above threshold should be in (0, 1): {}",
            d
        );
    }
    #[test]
    fn test_exponential_law_monotone_increasing() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let d1 = law.compute_damage(2e-4);
        let d2 = law.compute_damage(5e-4);
        let d3 = law.compute_damage(1e-3);
        assert!(
            d1 <= d2 && d2 <= d3,
            "Damage law must be monotone: d1={}, d2={}, d3={}",
            d1,
            d2,
            d3
        );
    }
    #[test]
    fn test_linear_law_at_ultimate_strain() {
        let kappa0 = 1e-4;
        let kappa_u = 1e-3;
        let law = DamageEvolutionLaw::linear(kappa0, kappa_u);
        let d = law.compute_damage(kappa_u + 1e-10);
        assert_eq!(d, 1.0, "Damage at ultimate strain should be 1");
    }
    #[test]
    fn test_power_law_zero_below_threshold() {
        let law = DamageEvolutionLaw::power_law(1e-4, 2.0);
        let d = law.compute_damage(1e-5);
        assert_eq!(d, 0.0, "Power-law damage below threshold should be 0");
    }
    #[test]
    fn test_damage_rate_positive_above_threshold() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let rate = law.compute_damage_rate(5e-4);
        assert!(
            rate > 0.0,
            "Damage rate above threshold should be positive: {}",
            rate
        );
    }
    #[test]
    fn test_damage_rate_zero_below_threshold() {
        let law = DamageEvolutionLaw::exponential(1e-4, 0.99, 500.0);
        let rate = law.compute_damage_rate(5e-5);
        assert_eq!(rate, 0.0, "Damage rate below threshold should be zero");
    }
    #[test]
    fn test_creep_strain_rate_positive() {
        let model = DamageCreepCoupling::new(1e-28, 3.5, 1.5, 1e-35, 3.0, 2.0);
        let rate = model.creep_strain_rate(100e6, 0.0);
        assert!(
            rate > 0.0 && rate.is_finite(),
            "Creep strain rate should be positive: {}",
            rate
        );
    }
    #[test]
    fn test_time_to_rupture_finite_for_nonzero_stress() {
        let model = DamageCreepCoupling::new(1e-28, 3.5, 1.5, 1e-35, 3.0, 2.0);
        let t_r = model.time_to_rupture(100e6);
        assert!(
            t_r > 0.0 && t_r.is_finite(),
            "Rupture time should be finite positive: {}",
            t_r
        );
    }
    #[test]
    fn test_explicit_step_damage_increases() {
        let model = DamageCreepCoupling::new(1e-28, 3.5, 1.5, 1e-20, 3.0, 2.0);
        let sigma = 200e6;
        let (_, d_new) = model.explicit_step(sigma, 0.0, 0.0, 1000.0);
        assert!(d_new > 0.0, "Damage should increase after explicit step");
    }
}
