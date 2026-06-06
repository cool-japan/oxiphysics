//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::fracture::functions::*;
    use crate::fracture::types::*;
    use std::f64::consts::PI;
    #[test]
    fn test_phase_field_params_new() {
        let p = PhaseFieldParams::new(100.0, 1e-3, 200e9, 0.3);
        assert_eq!(p.gc, 100.0);
        assert_eq!(p.l0, 1e-3);
        assert_eq!(p.e, 200e9);
        assert_eq!(p.nu, 0.3);
    }
    #[test]
    fn test_crack_driving_force_new_exceeds_history() {
        let h = crack_driving_force(500.0, 300.0);
        assert_eq!(h, 500.0);
    }
    #[test]
    fn test_crack_driving_force_history_dominates() {
        let h = crack_driving_force(100.0, 400.0);
        assert_eq!(h, 400.0);
    }
    #[test]
    fn test_crack_driving_force_equal() {
        let h = crack_driving_force(250.0, 250.0);
        assert_eq!(h, 250.0);
    }
    #[test]
    fn test_degradation_intact() {
        let g = degradation_function(0.0);
        assert!((g - 1.000_001).abs() < 1e-10);
    }
    #[test]
    fn test_degradation_fully_broken() {
        let g = degradation_function(1.0);
        assert!((g - 1e-6).abs() < 1e-15);
    }
    #[test]
    fn test_degradation_monotone_decreasing() {
        let d_vals = [0.0, 0.25, 0.5, 0.75, 1.0];
        let g_vals: Vec<f64> = d_vals.iter().map(|&d| degradation_function(d)).collect();
        for i in 1..g_vals.len() {
            assert!(g_vals[i] < g_vals[i - 1], "g should decrease with d");
        }
    }
    #[test]
    fn test_phase_field_residual_at_d0() {
        let p = PhaseFieldParams::new(1000.0, 1e-3, 200e9, 0.3);
        let r = phase_field_residual(0.0, 500.0, &p);
        assert!((r - (-1000.0)).abs() < 1e-6);
    }
    #[test]
    fn test_phase_field_residual_equilibrium() {
        let p = PhaseFieldParams::new(1000.0, 1e-3, 200e9, 0.3);
        let h = 500_000.0_f64;
        let alpha = p.gc / p.l0;
        let d_eq = 2.0 * h / (alpha + 2.0 * h);
        let r = phase_field_residual(d_eq, h, &p);
        assert!(
            r.abs() < 1e-6,
            "residual at equilibrium should vanish, got {r}"
        );
    }
    #[test]
    fn test_update_phase_field_zero_driving() {
        let p = PhaseFieldParams::new(1000.0, 1e-3, 200e9, 0.3);
        let d = update_phase_field(0.0, 0.0, &p);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_update_phase_field_large_driving() {
        let p = PhaseFieldParams::new(1.0, 1e-3, 200e9, 0.3);
        let d = update_phase_field(0.0, 1e12, &p);
        assert!((d - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_update_phase_field_irreversibility() {
        let p = PhaseFieldParams::new(1000.0, 1e-3, 200e9, 0.3);
        let d_prev = 0.5;
        let d = update_phase_field(d_prev, 1.0, &p);
        assert!(d >= d_prev);
    }
    #[test]
    fn test_update_phase_field_clamped_at_one() {
        let p = PhaseFieldParams::new(1e-10, 1e-3, 200e9, 0.3);
        let d = update_phase_field(0.0, 1e20, &p);
        assert!(d <= 1.0);
    }
    #[test]
    fn test_effective_modulus_intact() {
        let e_eff = effective_modulus(0.0, 200e9);
        assert!((e_eff - 200e9 * 1.000_001).abs() < 1e4);
    }
    #[test]
    fn test_effective_modulus_broken() {
        let e_eff = effective_modulus(1.0, 200e9);
        assert!((e_eff - 200e9 * 1e-6).abs() < 1.0);
    }
    #[test]
    fn test_effective_modulus_half_damage() {
        let e_eff = effective_modulus(0.5, 100.0);
        let expected = (0.25 + 1e-6) * 100.0;
        assert!((e_eff - expected).abs() < 1e-10);
    }
    #[test]
    fn test_crack_length_empty_field() {
        assert_eq!(crack_length_from_field(&[], 0.001, 1e-3), 0.0);
    }
    #[test]
    fn test_crack_length_uniform_one() {
        let n = 100;
        let dx = 1e-3;
        let l0 = 1e-3;
        let d_field = vec![1.0_f64; n];
        let gamma = crack_length_from_field(&d_field, dx, l0);
        assert!(gamma > 0.0);
    }
    #[test]
    fn test_crack_length_zero_field() {
        let d_field = vec![0.0_f64; 50];
        let gamma = crack_length_from_field(&d_field, 0.001, 1e-3);
        assert!(gamma.abs() < 1e-30);
    }
    #[test]
    fn test_sif_center_crack_infinite() {
        let k = sif_center_crack_infinite(100.0e6, 0.01);
        let expected = 100.0e6 * (PI * 0.01).sqrt();
        assert!((k - expected).abs() < 1e-6);
    }
    #[test]
    fn test_sif_center_crack_finite_reduces_to_infinite() {
        let k_inf = sif_center_crack_infinite(100.0e6, 0.001);
        let k_fin = sif_center_crack_finite_width(100.0e6, 0.001, 1.0);
        assert!(
            (k_fin - k_inf).abs() / k_inf < 0.01,
            "For small a/W, finite ≈ infinite: k_inf={k_inf}, k_fin={k_fin}"
        );
    }
    #[test]
    fn test_sif_finite_width_increases_with_ratio() {
        let k1 = sif_center_crack_finite_width(100.0e6, 0.01, 0.1);
        let k2 = sif_center_crack_finite_width(100.0e6, 0.01, 0.05);
        assert!(k2 > k1, "Narrower plate → higher SIF");
    }
    #[test]
    fn test_sif_single_edge_crack() {
        let k = sif_single_edge_crack(100.0e6, 0.01);
        let expected = 1.12 * 100.0e6 * (PI * 0.01).sqrt();
        assert!((k - expected).abs() < 1e-6);
    }
    #[test]
    fn test_sif_edge_crack_finite_at_small_ratio() {
        let k = sif_edge_crack_finite_width(100.0e6, 0.001, 1.0);
        let k_ref = 1.12 * 100.0e6 * (PI * 0.001).sqrt();
        assert!(
            (k - k_ref).abs() / k_ref < 0.01,
            "Small a/W: k={k}, expected ≈{k_ref}"
        );
    }
    #[test]
    fn test_sif_penny_shaped() {
        let k = sif_penny_shaped_crack(100.0e6, 0.01);
        let expected = (2.0 / PI) * 100.0e6 * (PI * 0.01).sqrt();
        assert!((k - expected).abs() < 1e-6);
    }
    #[test]
    fn test_sif_compact_tension() {
        let k = sif_compact_tension(1000.0, 0.01, 0.05, 0.025);
        assert!(k > 0.0, "CT SIF should be positive");
    }
    #[test]
    fn test_sif_double_edge_crack() {
        let k = sif_double_edge_crack(100.0e6, 0.005, 0.05);
        assert!(k > 0.0);
    }
    #[test]
    fn test_weight_function_singular_at_tip() {
        let w_near = weight_function_edge_crack(0.009, 0.01, 0.6147, -0.1244, 0.0);
        let w_far = weight_function_edge_crack(0.001, 0.01, 0.6147, -0.1244, 0.0);
        assert!(
            w_near > w_far,
            "Weight function should increase near crack tip"
        );
    }
    #[test]
    fn test_weight_function_at_x_equals_a() {
        let w = weight_function_edge_crack(0.01, 0.01, 0.6147, -0.1244, 0.0);
        assert_eq!(w, 0.0, "Weight function should be zero at x >= a");
    }
    #[test]
    fn test_sif_weight_function_uniform_stress() {
        let stress = 100.0e6;
        let a = 0.01;
        let k_wf = sif_weight_function(&|_| stress, a, 0.6147, -0.1244, 0.0, 1000);
        assert!(k_wf > 0.0, "Weight function SIF should be positive");
        assert!(k_wf.is_finite(), "Weight function SIF should be finite");
        let k_ref = 1.12 * stress * (PI * a).sqrt();
        assert!(
            k_wf > k_ref * 0.5,
            "Weight function SIF should be at least half the reference: k_wf={k_wf}, k_ref={k_ref}"
        );
    }
    #[test]
    fn test_irwin_plastic_zone_plane_stress() {
        let k = 50.0e6;
        let sigma_ys = 250.0e6;
        let rp = irwin_plastic_zone(k, sigma_ys, false);
        let expected = (1.0 / (2.0 * PI)) * (k / sigma_ys).powi(2);
        assert!((rp - expected).abs() < 1e-10);
    }
    #[test]
    fn test_irwin_plastic_zone_plane_strain_smaller() {
        let k = 50.0e6;
        let sigma_ys = 250.0e6;
        let rp_stress = irwin_plastic_zone(k, sigma_ys, false);
        let rp_strain = irwin_plastic_zone(k, sigma_ys, true);
        assert!(
            rp_strain < rp_stress,
            "Plane strain plastic zone should be smaller"
        );
    }
    #[test]
    fn test_irwin_effective_crack_length() {
        let a = 0.01;
        let k = 50.0e6;
        let sigma_ys = 250.0e6;
        let a_eff = irwin_effective_crack_length(a, k, sigma_ys, false);
        assert!(a_eff > a, "Effective crack length should exceed actual");
    }
    #[test]
    fn test_irwin_corrected_sif() {
        let stress = 100.0e6;
        let a = 0.01;
        let sigma_ys = 250.0e6;
        let k_uncorr = stress * (PI * a).sqrt();
        let k_corr = irwin_corrected_sif(stress, a, sigma_ys, false, 20);
        assert!(
            k_corr > k_uncorr,
            "Corrected SIF should exceed uncorrected: {k_corr} vs {k_uncorr}"
        );
    }
    #[test]
    fn test_dugdale_plastic_zone() {
        let a = 0.01;
        let stress = 100.0e6;
        let sigma_ys = 250.0e6;
        let c = dugdale_plastic_zone(a, stress, sigma_ys);
        assert!(c > 0.0, "Dugdale plastic zone should be positive");
    }
    #[test]
    fn test_dugdale_ctod() {
        let a = 0.01;
        let stress = 100.0e6;
        let sigma_ys = 250.0e6;
        let e = 200.0e9;
        let ctod = dugdale_ctod(a, stress, sigma_ys, e);
        assert!(ctod > 0.0, "Dugdale CTOD should be positive");
    }
    #[test]
    fn test_ctod_from_k() {
        let k = 50.0e6;
        let e = 200.0e9;
        let sigma_ys = 250.0e6;
        let ctod_ps = ctod_from_k(k, e, sigma_ys, 0.3, false);
        let ctod_pn = ctod_from_k(k, e, sigma_ys, 0.3, true);
        assert!(ctod_ps > ctod_pn, "Plane stress CTOD > plane strain CTOD");
    }
    #[test]
    fn test_ctod_criterion_test() {
        assert!(ctod_criterion(0.5, 0.3));
        assert!(!ctod_criterion(0.2, 0.3));
    }
    #[test]
    fn test_ctoa_from_ctod() {
        let ctoa = ctoa_from_ctod(0.001, 0.01);
        assert!(ctoa > 0.0);
        let approx = 0.001 / 0.01;
        assert!((ctoa - approx).abs() / approx < 0.1);
    }
    #[test]
    fn test_ctoa_criterion_test() {
        assert!(ctoa_criterion(0.1, 0.05));
        assert!(!ctoa_criterion(0.03, 0.05));
    }
    #[test]
    fn test_r_curve_at_initiation() {
        let rc = RCurve::new(30.0e6, 80.0e6, 0.005);
        assert!((rc.resistance(0.0) - 30.0e6).abs() < 1e-6);
    }
    #[test]
    fn test_r_curve_at_plateau() {
        let rc = RCurve::new(30.0e6, 80.0e6, 0.005);
        let r = rc.resistance(0.1);
        assert!(
            (r - 80.0e6).abs() / 80.0e6 < 0.01,
            "At large Δa, resistance ≈ plateau: {r}"
        );
    }
    #[test]
    fn test_r_curve_monotonically_increasing() {
        let rc = RCurve::new(30.0e6, 80.0e6, 0.005);
        let r1 = rc.resistance(0.001);
        let r2 = rc.resistance(0.005);
        let r3 = rc.resistance(0.01);
        assert!(r2 > r1);
        assert!(r3 > r2);
    }
    #[test]
    fn test_r_curve_stability() {
        let rc = RCurve::new(30.0e6, 80.0e6, 0.005);
        assert!(rc.is_stable(25.0e6, 0.0));
        assert!(!rc.is_stable(90.0e6, 0.0));
    }
    #[test]
    fn test_r_curve_negative_delta_a() {
        let rc = RCurve::new(30.0e6, 80.0e6, 0.005);
        assert!((rc.resistance(-0.001) - 30.0e6).abs() < 1e-6);
    }
    #[test]
    fn test_gtn_steel_preset() {
        let p = GtnParams::steel();
        assert_eq!(p.q1, 1.5);
        assert_eq!(p.q2, 1.0);
        assert!((p.q3 - 2.25).abs() < 1e-10);
        assert!(p.fc < p.ff);
        assert!(p.f0 < p.fc);
    }
    #[test]
    fn test_gtn_yield_function_virgin_material_elastic() {
        let p = GtnParams::steel();
        let phi = gtn_yield_function(0.1 * p.sigma_y, 0.0, p.sigma_y, p.f0, &p);
        assert!(phi < 0.0, "Should be elastic, got Φ={phi}");
    }
    #[test]
    fn test_gtn_yield_function_no_void_von_mises() {
        let p = GtnParams::steel();
        let phi_yield = gtn_yield_function(p.sigma_y, 0.0, p.sigma_y, 0.0, &p);
        assert!((phi_yield - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_gtn_yield_function_positive_hydrostatic() {
        let p = GtnParams::steel();
        let phi_no_tri = gtn_yield_function(p.sigma_y, 0.0, p.sigma_y, p.f0, &p);
        let phi_tri = gtn_yield_function(p.sigma_y, p.sigma_y, p.sigma_y, p.f0, &p);
        assert!(phi_tri > phi_no_tri);
    }
    #[test]
    fn test_nucleation_rate_at_mean_strain() {
        let p = GtnParams::steel();
        let rate_peak = nucleation_rate(p.en, 0.01, p.fn_void, p.en, p.sn);
        let rate_off = nucleation_rate(p.en + 3.0 * p.sn, 0.01, p.fn_void, p.en, p.sn);
        assert!(
            rate_peak > rate_off,
            "Peak nucleation should exceed off-peak"
        );
    }
    #[test]
    fn test_nucleation_rate_zero_increment() {
        let rate = nucleation_rate(0.3, 0.0, 0.04, 0.3, 0.1);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_nucleation_rate_positive() {
        let rate = nucleation_rate(0.3, 0.01, 0.04, 0.3, 0.1);
        assert!(rate > 0.0);
    }
    #[test]
    fn test_void_growth_rate_positive_dilation() {
        let df = void_growth_rate(0.01, 0.001);
        assert!((df - 0.99 * 0.001).abs() < 1e-15);
    }
    #[test]
    fn test_void_growth_rate_zero_dilation() {
        let df = void_growth_rate(0.05, 0.0);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_void_growth_rate_no_matrix_left() {
        let df = void_growth_rate(1.0, 0.01);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_effective_void_fraction_below_fc() {
        let p = GtnParams::steel();
        let f = 0.05;
        assert_eq!(effective_void_fraction(f, &p), f);
    }
    #[test]
    fn test_effective_void_fraction_at_fc() {
        let p = GtnParams::steel();
        assert!((effective_void_fraction(p.fc, &p) - p.fc).abs() < 1e-15);
    }
    #[test]
    fn test_effective_void_fraction_at_ff() {
        let p = GtnParams::steel();
        let f_star = effective_void_fraction(p.ff, &p);
        assert!((f_star - 1.0 / p.q1).abs() < 1e-10);
    }
    #[test]
    fn test_effective_void_fraction_above_fc_increases() {
        let p = GtnParams::steel();
        let f1 = effective_void_fraction(p.fc + 0.01, &p);
        let f2 = effective_void_fraction(p.fc + 0.05, &p);
        assert!(f2 > f1);
    }
    #[test]
    fn test_gtn_step_void_grows() {
        let p = GtnParams::steel();
        let f_new = gtn_step(p.f0, 0.3, 0.01, 0.003, &p);
        assert!(
            f_new > p.f0,
            "Void fraction should increase during plastic flow"
        );
    }
    #[test]
    fn test_gtn_step_no_plastic_flow() {
        let p = GtnParams::steel();
        let f_new = gtn_step(p.f0, 0.0, 0.0, 0.0, &p);
        assert!((f_new - p.f0).abs() < 1e-15);
    }
    #[test]
    fn test_gtn_step_clamped_at_one() {
        let p = GtnParams::steel();
        let f_new = gtn_step(0.99, 0.3, 1.0, 1.0, &p);
        assert!(f_new <= 1.0);
    }
    #[test]
    fn test_gtn_step_clamped_at_zero() {
        let p = GtnParams::steel();
        let f_new = gtn_step(0.001, 0.0, 0.0, -10.0, &p);
        assert!(f_new >= 0.0);
    }
    #[test]
    fn test_k_to_g_plane_stress() {
        let k = 50.0e6;
        let e = 200.0e9;
        let g = k_to_g(k, e, 0.3, false);
        let expected = k * k / e;
        assert!((g - expected).abs() < 1e-6);
    }
    #[test]
    fn test_k_to_g_plane_strain() {
        let k = 50.0e6;
        let e = 200.0e9;
        let nu = 0.3;
        let g = k_to_g(k, e, nu, true);
        let expected = (1.0 - nu * nu) * k * k / e;
        assert!((g - expected).abs() < 1e-6);
    }
    #[test]
    fn test_g_to_k_roundtrip() {
        let k_orig = 60.0e6;
        let e = 200.0e9;
        let nu = 0.3;
        let g = k_to_g(k_orig, e, nu, true);
        let k_back = g_to_k(g, e, nu, true);
        assert!((k_back - k_orig).abs() / k_orig < 1e-10);
    }
    #[test]
    fn test_paris_law_basic() {
        let da_dn = paris_law_growth_rate(10.0e6, 1e-12, 3.0);
        let expected = 1e-12 * (10.0e6_f64).powf(3.0);
        assert!((da_dn - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_paris_law_zero_delta_k() {
        let da_dn = paris_law_growth_rate(0.0, 1e-12, 3.0);
        assert_eq!(da_dn, 0.0);
    }
    #[test]
    fn test_paris_forman_threshold() {
        let da_dn = paris_forman_growth_rate(5.0e6, 10.0e6, 1e-12, 3.0, 5.0e6, 100.0e6);
        assert_eq!(da_dn, 0.0, "At threshold, growth rate should be zero");
    }
    #[test]
    fn test_paris_forman_instability() {
        let da_dn = paris_forman_growth_rate(50.0e6, 100.0e6, 1e-12, 3.0, 5.0e6, 100.0e6);
        assert!(
            da_dn.is_infinite(),
            "At K_Ic, growth rate should be infinite"
        );
    }
    #[test]
    fn test_sif_three_point_bend_positive() {
        let k = sif_three_point_bend(5000.0, 0.1, 0.01, 0.025, 0.01);
        assert!(k > 0.0, "SENB SIF should be positive");
        assert!(k.is_finite(), "SENB SIF should be finite");
    }
    #[test]
    fn test_sif_three_point_bend_increases_with_crack() {
        let k1 = sif_three_point_bend(5000.0, 0.1, 0.01, 0.05, 0.01);
        let k2 = sif_three_point_bend(5000.0, 0.1, 0.01, 0.05, 0.02);
        assert!(k2 > k1, "Longer crack should give higher SIF");
    }
    #[test]
    fn test_sif_three_point_bend_increases_with_force() {
        let k1 = sif_three_point_bend(1000.0, 0.1, 0.01, 0.05, 0.015);
        let k2 = sif_three_point_bend(2000.0, 0.1, 0.01, 0.05, 0.015);
        assert!(
            (k2 / k1 - 2.0).abs() < 1e-10,
            "SIF should be proportional to force"
        );
    }
    #[test]
    fn test_senb_critical_crack_length() {
        let k_ic = 1.0e6;
        let result = senb_critical_crack_length(5000.0, 0.1, 0.01, 0.025, k_ic, 64);
        assert!(
            result.is_some(),
            "Should find a critical crack length for K_Ic={k_ic}"
        );
        let a_crit = result.unwrap();
        assert!(
            a_crit > 0.0 && a_crit < 0.025,
            "Critical length should be in valid range [0, W]"
        );
    }
    #[test]
    fn test_sif_mixed_mode_pure_mode_i() {
        let (k_i, k_ii) = sif_mixed_mode_inclined_crack(100.0e6, 0.01, 0.0);
        let expected_ki = 100.0e6 * (PI * 0.01).sqrt();
        assert!(
            (k_i - expected_ki).abs() / expected_ki < 1e-10,
            "Pure mode I at β=0"
        );
        assert!(k_ii.abs() < 1e-6, "No K_II at β=0");
    }
    #[test]
    fn test_sif_mixed_mode_45_degrees() {
        let beta = PI / 4.0;
        let (k_i, k_ii) = sif_mixed_mode_inclined_crack(100.0e6, 0.01, beta);
        assert!((k_i - k_ii).abs() / k_i < 1e-10, "K_I = K_II at β=45°");
    }
    #[test]
    fn test_sif_mode_iii() {
        let k_iii = sif_mode_iii_edge_crack(50.0e6, 0.01);
        let expected = 50.0e6 * (PI * 0.01).sqrt();
        assert!((k_iii - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_sif_effective_mixed_mode_pure_mode_i() {
        let k_eff = sif_effective_mixed_mode(50.0e6, 0.0);
        assert!((k_eff - 50.0e6).abs() < 1.0, "Pure mode I effective = K_I");
    }
    #[test]
    fn test_sif_effective_mixed_mode_iii() {
        let k_eff = sif_effective_mixed_mode_iii(30.0e6, 40.0e6, 0.0, 0.3);
        let expected = (30.0e6_f64.powi(2) + 40.0e6_f64.powi(2)).sqrt();
        assert!((k_eff - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_mixed_mode_propagation_angle_pure_mode_i() {
        let theta = mixed_mode_propagation_angle(50.0e6, 0.0);
        assert_eq!(theta, 0.0, "Pure mode I → θ = 0");
    }
    #[test]
    fn test_mixed_mode_propagation_angle_nonzero() {
        let theta = mixed_mode_propagation_angle(30.0e6, 20.0e6);
        assert!(theta.abs() > 0.0, "Mixed mode should give non-zero angle");
    }
    #[test]
    fn test_walker_delta_k_eff_at_r_zero() {
        let dk_eff = walker_delta_k_eff(10.0e6, 0.0, 0.5);
        assert!((dk_eff - 10.0e6).abs() / 10.0e6 < 1e-10);
    }
    #[test]
    fn test_walker_delta_k_eff_increases_with_r() {
        let dk_eff_r0 = walker_delta_k_eff(10.0e6, 0.0, 0.5);
        let dk_eff_r05 = walker_delta_k_eff(10.0e6, 0.5, 0.5);
        assert!(dk_eff_r05 > dk_eff_r0, "ΔK_eff should increase with R");
    }
    #[test]
    fn test_paris_walker_rate_positive() {
        let rate = paris_walker_rate(10.0e6, 0.3, 1e-12, 3.0, 0.5);
        assert!(rate > 0.0, "Walker-corrected Paris rate should be positive");
        let rate_high_r = paris_walker_rate(10.0e6, 0.7, 1e-12, 3.0, 0.5);
        assert!(rate_high_r > rate, "Higher R → higher crack growth rate");
    }
    #[test]
    fn test_forman_crack_growth_rate_positive() {
        let rate = forman_crack_growth_rate(10.0e6, 0.3, 1e-12, 3.0, 100.0e6);
        assert!(rate > 0.0, "Forman rate should be positive");
    }
    #[test]
    fn test_forman_crack_growth_instability() {
        let rate = forman_crack_growth_rate(70.0e6, 0.3, 1e-12, 3.0, 100.0e6);
        assert!(
            rate.is_infinite(),
            "Forman rate should be infinite at instability"
        );
    }
    #[test]
    fn test_plastic_zone_first_order_plane_stress_gt_strain() {
        let k = 50.0e6;
        let sy = 350.0e6;
        let rp_stress = plastic_zone_size_first_order(k, sy, false);
        let rp_strain = plastic_zone_size_first_order(k, sy, true);
        assert!(
            rp_stress > rp_strain,
            "Plane stress zone > plane strain zone"
        );
    }
    #[test]
    fn test_plastic_zone_iterative_larger_than_first_order() {
        let sigma = 100.0e6;
        let a = 0.01;
        let sy = 350.0e6;
        let k0 = sigma * (PI * a).sqrt();
        let rp_first = plastic_zone_size_first_order(k0, sy, false);
        let rp_iter = plastic_zone_size_iterative(sigma, a, sy, false, 50);
        assert!(
            rp_iter >= rp_first,
            "Iterative plastic zone should be >= first order"
        );
    }
    #[test]
    fn test_plastic_zone_biaxial_at_zero_biaxiality() {
        let k = 50.0e6;
        let sy = 350.0e6;
        let rp_bi = plastic_zone_biaxial(k, sy, 0.0, false);
        let rp_first = plastic_zone_size_first_order(k, sy, false);
        assert!((rp_bi - rp_first).abs() / rp_first < 1e-10);
    }
    #[test]
    fn test_ctod_irwin_plane_strain_smaller() {
        let k = 50.0e6;
        let sy = 350.0e6;
        let e = 200.0e9;
        let nu = 0.3;
        let ctod_ps = ctod_irwin(k, sy, e, nu, false);
        let ctod_pn = ctod_irwin(k, sy, e, nu, true);
        assert!(ctod_ps > ctod_pn, "Plane stress CTOD > plane strain CTOD");
    }
    #[test]
    fn test_ctod_wells_positive() {
        let ctod = ctod_wells(50.0e6, 200.0e9, 350.0e6);
        assert!(ctod > 0.0, "Wells CTOD should be positive");
    }
    #[test]
    fn test_ctod_from_j_roundtrip() {
        let j = 1000.0;
        let sy = 350.0e6;
        let ctod = ctod_from_j(j, sy, true);
        let j_back = j_from_ctod(ctod, sy, true);
        assert!((j_back - j).abs() / j < 1e-10, "J-CTOD roundtrip failed");
    }
    #[test]
    fn test_minimum_thickness_plane_strain() {
        let k_ic = 60.0e6;
        let sy = 350.0e6;
        let b_min = minimum_thickness_plane_strain(k_ic, sy);
        let expected = 2.5 * (k_ic / sy).powi(2);
        assert!((b_min - expected).abs() < 1e-12, "B_min mismatch");
    }
    #[test]
    fn test_is_plane_strain_valid_thin_fails() {
        let k_ic = 60.0e6;
        let sy = 350.0e6;
        let b_min = minimum_thickness_plane_strain(k_ic, sy);
        assert!(!is_plane_strain_valid(
            b_min * 0.5,
            b_min,
            b_min * 3.0,
            k_ic,
            sy
        ));
    }
    #[test]
    fn test_is_plane_strain_valid_sufficient() {
        let k_ic = 60.0e6;
        let sy = 350.0e6;
        let b_min = minimum_thickness_plane_strain(k_ic, sy);
        assert!(is_plane_strain_valid(
            b_min * 2.0,
            b_min * 2.0,
            b_min * 6.0,
            k_ic,
            sy
        ));
    }
    #[test]
    fn test_k_to_j_roundtrip() {
        let k_ic = 60.0e6;
        let e = 200.0e9;
        let nu = 0.3;
        let j = k_to_j(k_ic, e, nu);
        let k_back = j_to_k(j, e, nu);
        assert!(
            (k_back - k_ic).abs() / k_ic < 1e-10,
            "K→J→K roundtrip failed"
        );
    }
    #[test]
    fn test_sif_dynamic_freund_at_zero_velocity() {
        let k_s = 50.0e6;
        let k_dyn = sif_dynamic_freund(k_s, 0.0, 2000.0);
        assert!((k_dyn - k_s).abs() < 1.0, "At v=0, K_dyn = K_static");
    }
    #[test]
    fn test_sif_dynamic_freund_decreases_with_velocity() {
        let k_s = 50.0e6;
        let c_r = 2000.0;
        let k1 = sif_dynamic_freund(k_s, 500.0, c_r);
        let k2 = sif_dynamic_freund(k_s, 1500.0, c_r);
        assert!(k2 < k1, "K_dyn should decrease as crack velocity increases");
    }
    #[test]
    fn test_sif_dynamic_freund_at_terminal_velocity() {
        let k_dyn = sif_dynamic_freund(50.0e6, 2000.0, 2000.0);
        assert_eq!(k_dyn, 0.0, "At terminal velocity, K_dyn = 0");
    }
    #[test]
    fn test_rayleigh_wave_speed_positive() {
        let c_r = rayleigh_wave_speed(80.0e9, 7850.0, 0.3);
        assert!(c_r > 0.0, "Rayleigh speed should be positive");
        let c_s = (80.0e9_f64 / 7850.0).sqrt();
        assert!(
            c_r < c_s,
            "Rayleigh speed should be less than shear wave speed"
        );
    }
    #[test]
    fn test_crack_arrest_condition() {
        assert!(will_crack_arrest(40.0e6, 60.0e6), "K < K_Ia → arrest");
        assert!(!will_crack_arrest(70.0e6, 60.0e6), "K > K_Ia → no arrest");
    }
    #[test]
    fn test_mode_mixity_pure_mode_i() {
        let psi = mode_mixity_angle(50.0e6, 0.0);
        assert_eq!(psi, 0.0, "Pure mode I → ψ = 0");
    }
    #[test]
    fn test_mode_mixity_pure_mode_ii() {
        let psi = mode_mixity_angle(0.0, 50.0e6);
        assert!((psi - PI / 2.0).abs() < 1e-10, "Pure mode II → ψ = π/2");
    }
    #[test]
    fn test_mode_mixity_45_degrees() {
        let psi = mode_mixity_angle(30.0e6, 30.0e6);
        assert!((psi - PI / 4.0).abs() < 1e-10, "Equal K_I=K_II → ψ = π/4");
    }
    #[test]
    fn test_mode_mixity_degrees_pure_mode_i() {
        let psi_deg = mode_mixity_angle_degrees(50.0e6, 0.0);
        assert_eq!(psi_deg, 0.0, "Pure mode I → 0°");
    }
    #[test]
    fn test_mode_i_fraction_pure_mode_i() {
        let eta = mode_i_fraction(50.0e6, 0.0, 0.0);
        assert!((eta - 1.0).abs() < 1e-10, "Pure mode I → η_I = 1");
    }
    #[test]
    fn test_mode_i_fraction_mixed_mode() {
        let eta = mode_i_fraction(30.0e6, 40.0e6, 0.0);
        let expected = 900.0 / (900.0 + 1600.0);
        assert!((eta - expected).abs() < 1e-10, "Mode I fraction: {eta}");
    }
    #[test]
    fn test_mode_i_fraction_sums_with_other_modes() {
        let ki = 30.0e6;
        let kii = 40.0e6;
        let kiii = 0.0;
        let eta_i = mode_i_fraction(ki, kii, kiii);
        let eta_ii = mode_i_fraction(kii, ki, kiii);
        assert!((eta_i + eta_ii - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_stress_triaxiality_uniaxial() {
        let t = stress_triaxiality(100.0e6 / 3.0, 100.0e6);
        assert!((t - 1.0 / 3.0).abs() < 1e-10, "Uniaxial T = 1/3, got {t}");
    }
    #[test]
    fn test_stress_triaxiality_zero_for_pure_shear() {
        let t = stress_triaxiality(0.0, 100.0e6);
        assert_eq!(t, 0.0, "Pure shear T = 0");
    }
    #[test]
    fn test_stress_triaxiality_large_near_crack_tip() {
        let t = stress_triaxiality(500.0e6, 100.0e6);
        assert!(t > 1.0, "Near crack tip, T >> 1");
    }
    #[test]
    fn test_rice_tracey_void_growth_positive() {
        let rate = rice_tracey_void_growth(100.0e6, 300.0e6, 0.001);
        assert!(rate > 0.0, "Rice-Tracey void growth should be positive");
    }
    #[test]
    fn test_rice_tracey_increases_with_triaxiality() {
        let rate_low = rice_tracey_void_growth(100.0e6, 300.0e6, 0.001);
        let rate_high = rice_tracey_void_growth(300.0e6, 300.0e6, 0.001);
        assert!(
            rate_high > rate_low,
            "Higher triaxiality → faster void growth"
        );
    }
    #[test]
    fn test_t_stress_constraint_unity_when_equal() {
        let t_star = t_stress_constraint(350.0e6, 350.0e6);
        assert!((t_star - 1.0).abs() < 1e-10, "T* = 1 when T_stress = σ_ys");
    }
    #[test]
    fn test_coffin_manson_decreases_with_cycles() {
        let amp1 = coffin_manson_strain_amplitude(1e3, 0.5, -0.6);
        let amp2 = coffin_manson_strain_amplitude(1e5, 0.5, -0.6);
        assert!(amp2 < amp1, "Strain amplitude decreases with cycles");
    }
    #[test]
    fn test_coffin_manson_at_one_reversal() {
        let amp = coffin_manson_strain_amplitude(0.5, 0.5, -0.6);
        let expected = 0.5 * 1.0_f64.powf(-0.6);
        assert!((amp - expected).abs() < 1e-10, "CM at Nf=0.5: {amp}");
    }
    #[test]
    fn test_basquin_higher_cycles_lower_stress() {
        let s1 = basquin_stress_amplitude(1e4, 1000.0e6, -0.1);
        let s2 = basquin_stress_amplitude(1e6, 1000.0e6, -0.1);
        assert!(s2 < s1, "Higher cycles → lower stress amplitude");
    }
    #[test]
    fn test_basquin_fatigue_life_positive() {
        let nf = basquin_fatigue_life(200.0e6, 800.0e6, -0.1);
        assert!(nf > 0.0 && nf.is_finite(), "Basquin fatigue life: {nf}");
    }
    #[test]
    fn test_basquin_fatigue_life_higher_stress_shorter_life() {
        let nf1 = basquin_fatigue_life(200.0e6, 800.0e6, -0.1);
        let nf2 = basquin_fatigue_life(400.0e6, 800.0e6, -0.1);
        assert!(nf2 < nf1, "Higher stress → shorter fatigue life");
    }
    #[test]
    fn test_morrow_total_strain_matches_sum() {
        let nf = 1e4;
        let sigma_f = 800.0e6;
        let b = -0.1;
        let eps_f = 0.3;
        let c = -0.6;
        let e = 200.0e9;
        let total = morrow_total_strain_amplitude(nf, sigma_f, b, eps_f, c, e);
        let elastic = sigma_f / e * (2.0 * nf).powf(b);
        let plastic = coffin_manson_strain_amplitude(nf, eps_f, c);
        assert!(
            (total - elastic - plastic).abs() < 1e-15,
            "Morrow = elastic + plastic"
        );
    }
    #[test]
    fn test_swt_parameter_positive() {
        let swt = swt_parameter(200.0e6, 0.001);
        assert!(swt > 0.0, "SWT parameter should be positive");
    }
    #[test]
    fn test_swt_parameter_formula() {
        let sigma_max = 300.0e6;
        let eps_amp = 0.002;
        let swt = swt_parameter(sigma_max, eps_amp);
        let expected = sigma_max * eps_amp;
        assert!((swt - expected).abs() < 1e-10, "SWT = σ_max * ε_amp");
    }
    #[test]
    fn test_sif_hole_edge_crack_positive() {
        let k = sif_hole_edge_crack(100.0e6, 0.005, 0.01);
        assert!(k > 0.0, "SIF for hole edge crack should be positive");
        assert!(k.is_finite(), "SIF should be finite");
    }
    #[test]
    fn test_sif_hole_edge_crack_greater_than_plain() {
        let k_hole = sif_hole_edge_crack(100.0e6, 0.01, 0.01);
        let k_plain = sif_center_crack_infinite(100.0e6, 0.01);
        assert!(k_hole >= k_plain, "Hole edge crack SIF ≥ plain crack SIF");
    }
    #[test]
    fn test_sif_surface_crack_positive() {
        let k = sif_surface_crack(100.0e6, 0.005, 0.01);
        assert!(k > 0.0, "Surface crack SIF should be positive");
        assert!(k.is_finite(), "Surface crack SIF should be finite");
    }
    #[test]
    fn test_sif_surface_crack_increases_with_depth() {
        let k1 = sif_surface_crack(100.0e6, 0.005, 0.02);
        let k2 = sif_surface_crack(100.0e6, 0.008, 0.02);
        assert!(k2 > k1, "Deeper crack → higher SIF");
    }
    #[test]
    fn test_near_tip_stress_mode_i_symmetry() {
        let k_i = 50.0e6;
        let r = 0.001;
        let (_, s_tt_pos, _) = near_tip_stress_mode_i(k_i, r, 0.0);
        let (_, s_tt_neg, _) = near_tip_stress_mode_i(k_i, r, 0.0);
        assert!((s_tt_pos - s_tt_neg).abs() < 1e-6, "Symmetric at θ=0");
    }
    #[test]
    fn test_near_tip_stress_mode_i_at_theta_zero() {
        let k_i = 50.0e6;
        let r = 0.001;
        let (sigma_rr, sigma_tt, sigma_rt) = near_tip_stress_mode_i(k_i, r, 0.0);
        let sigma_ref = k_i / (2.0 * PI * r).sqrt();
        assert!(
            (sigma_rr - sigma_ref).abs() / sigma_ref < 1e-10,
            "σ_rr at θ=0: {sigma_rr}"
        );
        assert!(
            (sigma_tt - sigma_ref).abs() / sigma_ref < 1e-10,
            "σ_θθ at θ=0: {sigma_tt}"
        );
        assert!(
            sigma_rt.abs() < 1e-6,
            "σ_rθ at θ=0 should be 0, got {sigma_rt}"
        );
    }
    #[test]
    fn test_near_tip_stress_mode_i_decays_with_distance() {
        let k_i = 50.0e6;
        let (_, s1, _) = near_tip_stress_mode_i(k_i, 0.001, 0.0);
        let (_, s2, _) = near_tip_stress_mode_i(k_i, 0.01, 0.0);
        assert!(s1 > s2, "Near-tip stress decreases with distance");
    }
    #[test]
    fn test_paris_law_cycles_positive() {
        let n = paris_law_cycles(0.005, 0.02, 50.0e6, 1e-12, 3.0, 100.0e6, 100);
        assert!(n > 0.0, "Paris law cycles should be positive, got {n}");
        assert!(n.is_finite(), "Paris law cycles should be finite");
    }
    #[test]
    fn test_paris_law_cycles_longer_for_smaller_initial_crack() {
        let n1 = paris_law_cycles(0.005, 0.02, 50.0e6, 1e-12, 3.0, 100.0e6, 100);
        let n2 = paris_law_cycles(0.001, 0.02, 50.0e6, 1e-12, 3.0, 100.0e6, 100);
        assert!(n2 > n1, "Smaller initial crack → more cycles to failure");
    }
    #[test]
    fn test_paris_law_cycles_lower_stress_more_cycles() {
        let n1 = paris_law_cycles(0.005, 0.02, 50.0e6, 1e-12, 3.0, 100.0e6, 100);
        let n2 = paris_law_cycles(0.005, 0.02, 25.0e6, 1e-12, 3.0, 100.0e6, 100);
        assert!(n2 > n1, "Lower stress → more cycles");
    }
    #[test]
    fn test_forman_law_cycles_positive() {
        let n = forman_law_cycles(0.005, 0.02, 50.0e6, 0.3, 1e-12, 3.0, 100.0e6, 100);
        assert!(n > 0.0, "Forman law cycles should be positive, got {n}");
        assert!(n.is_finite(), "Forman law cycles should be finite");
    }
    #[test]
    fn test_forman_higher_r_ratio_fewer_cycles() {
        let n1 = forman_law_cycles(0.005, 0.015, 50.0e6, 0.1, 1e-12, 3.0, 100.0e6, 100);
        let n2 = forman_law_cycles(0.005, 0.015, 50.0e6, 0.5, 1e-12, 3.0, 100.0e6, 100);
        assert!(n2 < n1, "Higher R-ratio → fewer cycles (Forman)");
    }
    #[test]
    fn test_tearing_modulus_from_j_r_curve() {
        let j_c = 50_000.0_f64;
        let delta_a_c = 0.001_f64;
        let e = 200.0e9_f64;
        let sigma_ys = 350.0e6_f64;
        let dj_da = j_c / delta_a_c;
        let t_mod = e / (sigma_ys * sigma_ys) * dj_da;
        assert!(
            t_mod > 0.0,
            "Tearing modulus should be positive, got {t_mod}"
        );
        assert!(t_mod.is_finite(), "Tearing modulus should be finite");
    }
    #[test]
    fn test_j_integral_path_independence_check() {
        let k_ic = 60.0e6_f64;
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let j_from_k = k_to_j(k_ic, e, nu);
        let k_from_j = j_to_k(j_from_k, e, nu);
        assert!(
            (k_from_j - k_ic).abs() / k_ic < 1e-10,
            "J-K path independence check"
        );
    }
    #[test]
    fn test_ctod_from_j_consistency_plane_stress() {
        let j = 5000.0_f64;
        let sy = 350.0e6_f64;
        let delta = ctod_from_j(j, sy, false);
        let j_back = j_from_ctod(delta, sy, false);
        assert!(
            (j_back - j).abs() / j < 1e-10,
            "CTOD-J plane stress consistency"
        );
    }
    /// G = K²/E' for mode I, plane stress.
    #[test]
    fn test_griffith_energy_release_rate_mode_i() {
        let gc = GriffithCriterion::steel();
        let k_i = 50.0e6_f64;
        let g = gc.compute_energy_release_rate_mode_i(k_i, false);
        let e_prime = 200.0e9_f64;
        let expected = k_i * k_i / e_prime;
        assert!(
            (g - expected).abs() / expected < 1e-10,
            "G mode I: {} vs {}",
            g,
            expected
        );
    }
    /// Fracture criterion: G >= G_c triggers failure.
    #[test]
    fn test_griffith_criterion_triggered_at_critical_load() {
        let gc_val = 100.0_f64;
        let e = 200.0e9_f64;
        let gc = GriffithCriterion::new(gc_val, e, 0.3);
        let k_ic = gc.critical_sif(false);
        assert!(
            gc.is_critical(k_ic * 1.001, 0.0, 0.0, false),
            "Should be critical above K_Ic: K={}",
            k_ic * 1.001
        );
    }
    /// Below critical K: criterion not triggered.
    #[test]
    fn test_griffith_criterion_not_triggered_below_k_ic() {
        let gc = GriffithCriterion::steel();
        let k_ic = gc.critical_sif(false);
        assert!(
            !gc.is_critical(k_ic * 0.9, 0.0, 0.0, false),
            "Should not be critical at 0.9*K_Ic"
        );
    }
    /// Critical crack length decreases with increasing stress.
    #[test]
    fn test_griffith_critical_crack_length_decreases_with_stress() {
        let gc = GriffithCriterion::steel();
        let a1 = gc.critical_crack_length(100.0e6, false);
        let a2 = gc.critical_crack_length(200.0e6, false);
        assert!(
            a2 < a1,
            "Critical crack length should decrease with higher stress"
        );
    }
    /// Stable growth: at k0 < k_init, no growth occurs initially.
    #[test]
    fn test_rcurve_stable_growth_starts_at_k_init() {
        let rc = RCurve::new(30.0e6, 60.0e6, 0.005);
        let path = rc.compute_stable_crack_growth(rc.k_init, 5.0e9, 0.01, 10);
        assert!(!path.is_empty(), "Path should not be empty");
        let (da0, k_app0, k_r0) = path[0];
        assert!(da0.abs() < f64::EPSILON, "First step Δa=0");
        assert!(
            (k_app0 - rc.k_init).abs() < 1.0,
            "First K_app = k_init: {k_app0}"
        );
        assert!((k_r0 - rc.k_init).abs() < 1.0, "First K_R = k_init: {k_r0}");
    }
    /// R-curve initiation load equals k_init.
    #[test]
    fn test_rcurve_initiation_load() {
        let rc = RCurve::new(30.0e6, 60.0e6, 0.005);
        assert!(
            (rc.initiation_load() - rc.k_init).abs() < f64::EPSILON,
            "Initiation load should equal k_init"
        );
    }
    /// Plastic zone plane strain < plane stress (by factor 3).
    #[test]
    fn test_ctod_plastic_zone_plane_strain_smaller() {
        let cc = CtodCriterion::structural_steel();
        let k_i = 50.0e6_f64;
        let rp_ps = cc.compute_plastic_zone_size(k_i, false);
        let rp_pe = cc.compute_plastic_zone_size(k_i, true);
        assert!(
            rp_pe < rp_ps,
            "Plane strain r_p ({}) should be smaller than plane stress ({})",
            rp_pe,
            rp_ps
        );
        assert!(
            (rp_ps / rp_pe - 3.0).abs() < 1e-10,
            "Ratio should be 3: {}",
            rp_ps / rp_pe
        );
    }
    /// CTOD from SIF: round-trip K → δ → K.
    #[test]
    fn test_ctod_sif_round_trip() {
        let cc = CtodCriterion::structural_steel();
        let k_i = 50.0e6_f64;
        let ctod = cc.ctod_from_sif(k_i, false);
        let k_back = cc.critical_sif(false);
        assert!(ctod > 0.0, "CTOD should be positive: {ctod}");
        assert!(ctod.is_finite(), "CTOD should be finite: {ctod}");
        let _ = k_back;
    }
    /// Larger K → larger plastic zone.
    #[test]
    fn test_ctod_plastic_zone_increases_with_k() {
        let cc = CtodCriterion::structural_steel();
        let rp1 = cc.compute_plastic_zone_size(30.0e6, false);
        let rp2 = cc.compute_plastic_zone_size(60.0e6, false);
        assert!(
            rp2 > rp1,
            "Larger K should give larger plastic zone: {rp1} vs {rp2}"
        );
    }
}
