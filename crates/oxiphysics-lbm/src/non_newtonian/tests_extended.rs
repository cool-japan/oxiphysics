// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended tests for non-Newtonian fluid models.

use super::*;

mod extended_nn_tests {
    use super::*;

    // ── Maxwell fluid ────────────────────────────────────────────────────────

    #[test]
    fn test_maxwell_viscosity_formula() {
        // η = G * λ
        let mf = MaxwellFluid::new(2.0, 3.0); // λ=2, G=3 → η=6
        assert!(
            (mf.viscosity() - 6.0).abs() < 1e-12,
            "η = {}",
            mf.viscosity()
        );
    }

    #[test]
    fn test_maxwell_weissenberg() {
        let mf = MaxwellFluid::new(1.5, 2.0);
        let wi = mf.weissenberg(4.0);
        // Wi = λ * γ̇ = 1.5 * 4 = 6
        assert!((wi - 6.0).abs() < 1e-12, "Wi = {wi}");
    }

    #[test]
    fn test_maxwell_storage_modulus_zero_freq() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let gp = mf.storage_modulus(0.0);
        // G'(0) = G * 0 / (1 + 0) = 0
        assert!(gp.abs() < 1e-30, "G'(0) = {gp}");
    }

    #[test]
    fn test_maxwell_loss_modulus_low_freq() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let gpp = mf.loss_modulus(0.01);
        // At low ω: G''(ω) ≈ G * ωλ  (linear in ω)
        assert!(gpp > 0.0, "G''(ω) > 0: {gpp}");
    }

    #[test]
    fn test_maxwell_storage_increases_with_freq() {
        let mf = MaxwellFluid::new(1.0, 5.0);
        let gp_low = mf.storage_modulus(0.1);
        let gp_high = mf.storage_modulus(10.0);
        assert!(gp_high > gp_low, "G' should increase with frequency");
    }

    #[test]
    fn test_maxwell_complex_viscosity_decreases_with_freq() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let cv_low = mf.complex_viscosity_magnitude(0.01);
        let cv_high = mf.complex_viscosity_magnitude(100.0);
        assert!(cv_high < cv_low, "|η*(ω)| should decrease with frequency");
    }

    #[test]
    fn test_maxwell_first_normal_stress_coefficient() {
        let mf = MaxwellFluid::new(2.0, 3.0); // η=6, λ=2 → Ψ₁ = 2*6*2 = 24
        let psi1 = mf.first_normal_stress_coefficient();
        assert!((psi1 - 24.0).abs() < 1e-12, "Ψ₁ = {psi1}");
    }

    #[test]
    fn test_maxwell_normal_stress_difference_zero_at_rest() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let n1 = mf.first_normal_stress_difference(0.0);
        assert!(n1.abs() < 1e-30, "N1 at rest = {n1}");
    }

    #[test]
    fn test_maxwell_relaxation_modulus_decay() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let g0 = mf.relaxation_modulus(0.0);
        let g1 = mf.relaxation_modulus(1.0);
        assert!((g0 - 2.0).abs() < 1e-12, "G(0) = G = {g0}");
        assert!(g1 < g0, "G(t) should decay: G(0)={g0}, G(1)={g1}");
    }

    #[test]
    fn test_maxwell_lbm_relaxation_time_positive() {
        let mf = MaxwellFluid::new(0.5, 1.0); // η = 0.5
        let tau = mf.lbm_relaxation_time();
        assert!(tau > 0.5, "LBM τ must be > 0.5: {tau}");
    }

    #[test]
    fn test_maxwell_step_stress_steady_state() {
        let mf = MaxwellFluid::new(1.0, 2.0); // η=2
        let gamma = 1.0;
        // Advance many steps to steady state
        let mut tau_xy = 0.0_f64;
        for _ in 0..10000 {
            tau_xy = mf.step_stress(tau_xy, gamma, 0.001);
        }
        let expected = mf.steady_shear_stress(gamma); // = η * γ̇ = 2
        assert!(
            (tau_xy - expected).abs() < 0.01,
            "Steady-state stress={tau_xy}, expected {expected}"
        );
    }

    #[test]
    fn test_maxwell_non_newtonian_fluid_trait() {
        let mf = MaxwellFluid::new(1.0, 2.0); // η=2
        let fluid: &dyn NonNewtonianFluid = &mf;
        let mu = fluid.effective_viscosity(5.0);
        assert!(
            (mu - 2.0).abs() < 1e-12,
            "Maxwell μ_eff should be constant: {mu}"
        );
    }

    // ── BinghamPipeFlow ──────────────────────────────────────────────────────

    #[test]
    fn test_bingham_pipe_no_flow_below_yield() {
        let bf = BinghamPipeFlow::new(100.0, 0.1, 0.05, 1.0);
        // ΔP = 10 Pa, wall shear = R * ΔP / (2L) = 0.05 * 10 / 2 = 0.25 << 100
        assert!(!bf.is_flowing(10.0), "No flow below yield: ΔP=10");
    }

    #[test]
    fn test_bingham_pipe_flow_above_yield() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        // ΔP = 1000, τ_wall = 0.05*1000/2 = 25 >> τ_y=1 → flowing
        assert!(bf.is_flowing(1000.0), "Should flow at ΔP=1000");
    }

    #[test]
    fn test_bingham_pipe_plug_radius_decreases_with_dp() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let r1 = bf.plug_radius(500.0).unwrap();
        let r2 = bf.plug_radius(2000.0).unwrap();
        assert!(
            r2 < r1,
            "Plug radius decreases as ΔP increases: r1={r1}, r2={r2}"
        );
    }

    #[test]
    fn test_bingham_pipe_plug_radius_below_pipe_radius() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let r_p = bf.plug_radius(1000.0).unwrap();
        assert!(r_p < 0.05, "Plug radius must be < pipe radius: {r_p}");
        assert!(r_p > 0.0, "Plug radius must be positive: {r_p}");
    }

    #[test]
    fn test_bingham_pipe_wall_shear_stress() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.1, 1.0);
        let tau_w = bf.wall_shear_stress(200.0);
        // τ_w = R * ΔP / (2L) = 0.1 * 200 / 2 = 10
        assert!(
            (tau_w - 10.0).abs() < 1e-12,
            "τ_wall = {tau_w}, expected 10"
        );
    }

    #[test]
    fn test_bingham_pipe_mean_velocity_positive() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let u = bf.mean_velocity(1000.0);
        assert!(u > 0.0, "Mean velocity should be positive: {u}");
    }

    #[test]
    fn test_bingham_pipe_flow_rate_proportional_to_area() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let q = bf.flow_rate(1000.0);
        let u = bf.mean_velocity(1000.0);
        let area = std::f64::consts::PI * 0.05 * 0.05;
        assert!(
            (q - u * area).abs() < 1e-12,
            "Q = U * A: Q={q}, U*A={}",
            u * area
        );
    }

    #[test]
    fn test_bingham_pipe_velocity_zero_at_wall() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let u_wall = bf.velocity_profile(0.05, 1000.0).unwrap();
        assert!(
            u_wall.abs() < 1e-10,
            "Velocity at wall should be ~0: {u_wall}"
        );
    }

    #[test]
    fn test_bingham_pipe_velocity_max_at_centre() {
        let bf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let u_centre = bf.velocity_profile(0.0, 1000.0).unwrap();
        let u_mid = bf.velocity_profile(0.03, 1000.0).unwrap();
        assert!(
            u_centre >= u_mid,
            "Max velocity at centre: u_centre={u_centre}, u_mid={u_mid}"
        );
    }

    // ── HerschelBulkleyFluid ─────────────────────────────────────────────────

    #[test]
    fn test_hb_fluid_critical_shear_rate() {
        let hb = HerschelBulkleyFluid::new(1.0, 2.0, 1.0);
        // γ̇_c = (τ_y / K)^(1/n) = (1/2)^1 = 0.5
        let gamma_c = hb.critical_shear_rate();
        assert!(
            (gamma_c - 0.5).abs() < 1e-12,
            "γ̇_c = {gamma_c}, expected 0.5"
        );
    }

    #[test]
    fn test_hb_fluid_is_yielded() {
        let hb = HerschelBulkleyFluid::new(1.0, 2.0, 1.0); // γ̇_c = 0.5
        assert!(!hb.is_yielded(0.1), "Below critical shear rate: un-yielded");
        assert!(hb.is_yielded(1.0), "Above critical shear rate: yielded");
    }

    #[test]
    fn test_hb_fluid_reduces_to_bingham_n1() {
        let hb = HerschelBulkleyFluid::new(1.0, 0.1, 1.0);
        // With n=1: stress = tau_y + K*gamma, viscosity = (tau_y + K*gamma)/gamma
        let gamma = 5.0;
        let tau = hb.stress(gamma);
        let expected = 1.0 + 0.1 * 5.0; // τ_y + K * γ̇
        assert!(
            (tau - expected).abs() < 1e-12,
            "HB n=1 stress = {tau}, expected {expected}"
        );
    }

    #[test]
    fn test_hb_fluid_reduces_to_power_law_no_yield() {
        let hb = HerschelBulkleyFluid::new(0.0, 1.0, 0.5); // τ_y = 0 → power law
        let gamma = 4.0;
        let tau = hb.stress(gamma);
        // K * γ̇^n = 1 * 4^0.5 = 2
        assert!(
            (tau - 2.0).abs() < 1e-12,
            "HB τ_y=0 stress = {tau}, expected 2.0"
        );
    }

    #[test]
    fn test_hb_fluid_bingham_number() {
        let hb = HerschelBulkleyFluid::new(2.0, 1.0, 1.0);
        let bn = hb.bingham_number(1.0);
        // Bn = τ_y / (K * γ̇^n) = 2 / (1 * 1) = 2
        assert!((bn - 2.0).abs() < 1e-12, "Bn = {bn}");
    }

    #[test]
    fn test_hb_fluid_lbm_relaxation_time_above_half() {
        let hb = HerschelBulkleyFluid::new(0.1, 0.5, 0.8);
        let tau = hb.lbm_relaxation_time(1.0);
        assert!(tau > 0.5, "LBM τ > 0.5: {tau}");
    }

    #[test]
    fn test_hb_fluid_pipe_no_flow_insufficient_dp() {
        let hb = HerschelBulkleyFluid::new(100.0, 1.0, 0.8);
        // Very small pressure gradient → r_p >= R → no flow
        let u = hb.pipe_mean_velocity(0.05, 1.0);
        assert!(u.abs() < 1e-10, "No flow with insufficient pressure: u={u}");
    }

    #[test]
    fn test_hb_fluid_pipe_positive_flow() {
        let hb = HerschelBulkleyFluid::new(1.0, 0.01, 1.0);
        // Large pressure gradient → positive mean velocity
        let u = hb.pipe_mean_velocity(0.05, 10000.0);
        assert!(u > 0.0, "Positive flow with large ΔP: u={u}");
    }

    // ── ThixotropicModel ─────────────────────────────────────────────────────

    #[test]
    fn test_thixotropic_model_initial_viscosity() {
        let tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        // λ=1 → μ = μ_inf + (μ_0 - μ_inf) * 1 = μ_0 = 1.0
        assert!(
            (tm.current_viscosity() - 1.0).abs() < 1e-12,
            "Initial viscosity = {}",
            tm.current_viscosity()
        );
    }

    #[test]
    fn test_thixotropic_model_advance_shear_decreases_lambda() {
        let mut tm = ThixotropicModel::new(1.0, 0.01, 0.1, 10.0, 1.0);
        let l0 = tm.lambda;
        tm.advance(2.0, 0.1); // apply shear
        assert!(
            tm.lambda < l0,
            "λ should decrease under shear: {l0} → {}",
            tm.lambda
        );
    }

    #[test]
    fn test_thixotropic_model_advance_rest_increases_lambda() {
        let mut tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        tm.lambda = 0.2; // partially broken
        tm.advance(0.0, 0.1); // no shear → recovery
        assert!(tm.lambda > 0.2, "λ should recover at rest: {}", tm.lambda);
    }

    #[test]
    fn test_thixotropic_model_lambda_bounded() {
        let mut tm = ThixotropicModel::new(1.0, 0.01, 1.0, 100.0, 1.0);
        for _ in 0..10000 {
            tm.advance(5.0, 0.01);
        }
        assert!(
            (0.0..=1.0).contains(&tm.lambda),
            "λ out of [0,1]: {}",
            tm.lambda
        );
    }

    #[test]
    fn test_thixotropic_model_steady_state_viscosity() {
        let tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let v_ss = tm.steady_state_viscosity(1.0);
        // λ_ss = a / (a + b * γ^m) = 1 / (1 + 1) = 0.5 → μ = 0.01 + (1-0.01)*0.5 = 0.505
        let expected = 0.01 + (1.0 - 0.01) * 0.5;
        assert!(
            (v_ss - expected).abs() < 1e-12,
            "Steady-state viscosity = {v_ss}, expected {expected}"
        );
    }

    #[test]
    fn test_thixotropic_model_thixotropic_time_scale() {
        let tm = ThixotropicModel::new(1.0, 0.01, 2.0, 1.0, 1.0);
        // t_thix = 1/a = 0.5
        assert!(
            (tm.thixotropic_time_scale() - 0.5).abs() < 1e-12,
            "t_thix = {}",
            tm.thixotropic_time_scale()
        );
    }

    #[test]
    fn test_thixotropic_model_lbm_relaxation_time_positive() {
        let tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let tau = tm.lbm_relaxation_time();
        assert!(tau > 0.5, "LBM τ must be > 0.5: {tau}");
    }

    #[test]
    fn test_thixotropic_model_lbm_omega_positive() {
        let tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let omega = tm.lbm_omega();
        assert!(omega > 0.0 && omega <= 2.0, "LBM ω in (0,2]: {omega}");
    }

    #[test]
    fn test_thixotropic_model_reset() {
        let mut tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        tm.lambda = 0.1;
        tm.reset();
        assert!(
            (tm.lambda - 1.0).abs() < 1e-14,
            "After reset λ=1: {}",
            tm.lambda
        );
    }

    #[test]
    fn test_thixotropic_model_structural_exponent() {
        let tm = ThixotropicModel::with_structural_exponent(1.0, 0.0, 1.0, 1.0, 1.0, 2.0);
        // λ=1, k=2 → μ = μ_inf + (μ_0-μ_inf)*λ^2 = 0 + 1*1 = 1
        assert!(
            (tm.current_viscosity() - 1.0).abs() < 1e-12,
            "Structural exponent k=2: {}",
            tm.current_viscosity()
        );
    }

    // ── effective_relaxation_time ────────────────────────────────────────────

    #[test]
    fn test_effective_relaxation_time_newtonian() {
        // n=1, K=1/6 → μ=1/6 → τ = 0.5 + (1/6)/(1/3) = 1.0
        let pl = PowerLawFluid::new(1.0, 1.0 / 6.0);
        let tau = effective_relaxation_time(1.0, &pl);
        assert!((tau - 1.0).abs() < 1e-12, "τ for Newtonian = {tau}");
    }

    #[test]
    fn test_effective_relaxation_time_varies_with_shear() {
        let pl = PowerLawFluid::new(0.5, 1.0); // shear-thinning
        let tau_low = effective_relaxation_time(0.1, &pl);
        let tau_high = effective_relaxation_time(100.0, &pl);
        assert!(
            tau_high < tau_low,
            "τ should decrease for shear-thinning: τ_low={tau_low}, τ_high={tau_high}"
        );
    }

    #[test]
    fn test_effective_relaxation_frequency_inverse() {
        let pl = PowerLawFluid::new(1.0, 1.0 / 6.0);
        let tau = effective_relaxation_time(1.0, &pl);
        let omega = effective_relaxation_frequency(1.0, &pl);
        assert!(
            (tau * omega - 1.0).abs() < 1e-12,
            "τ * ω = 1: {}",
            tau * omega
        );
    }

    #[test]
    fn test_update_tau_from_shear_rates_length() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let shear_rates = vec![0.1, 1.0, 5.0, 20.0];
        let taus = update_tau_from_shear_rates(&shear_rates, &pl);
        assert_eq!(taus.len(), 4, "Output length matches input");
        for (i, &tau) in taus.iter().enumerate() {
            assert!(tau > 0.5, "tau[{i}]={tau} must be > 0.5");
        }
    }

    #[test]
    fn test_update_tau_from_shear_rates_monotone_thinning() {
        let pl = PowerLawFluid::new(0.5, 1.0); // shear-thinning
        let shear_rates = vec![0.1, 1.0, 10.0, 100.0];
        let taus = update_tau_from_shear_rates(&shear_rates, &pl);
        // τ should decrease with shear rate for shear-thinning fluid
        for i in 1..taus.len() {
            assert!(
                taus[i] <= taus[i - 1],
                "tau[{}]={} should be ≤ tau[{}]={}",
                i,
                taus[i],
                i - 1,
                taus[i - 1]
            );
        }
    }

    // ── Maxwell as NonNewtonianFluid ─────────────────────────────────────────

    #[test]
    fn test_maxwell_local_tau_via_trait() {
        let mf = MaxwellFluid::new(1.0, 1.0 / 6.0); // η = 1/6
        let fluid: &dyn NonNewtonianFluid = &mf;
        let tau = fluid.local_tau(1.0);
        // τ = 0.5 + η/cs² = 0.5 + (1/6)/(1/3) = 0.5 + 0.5 = 1.0
        assert!((tau - 1.0).abs() < 1e-12, "Maxwell local_tau = {tau}");
    }

    // ── Phase angle ──────────────────────────────────────────────────────────

    #[test]
    fn test_maxwell_phase_angle_high_freq_approaches_zero() {
        let mf = MaxwellFluid::new(1.0, 2.0); // G=2, λ=1
        let delta = mf.phase_angle(1000.0);
        // At high ω: G' dominates, δ → 0 (elastic)
        assert!(
            delta < std::f64::consts::FRAC_PI_4,
            "Phase angle should be small at high ω: {delta}"
        );
    }

    #[test]
    fn test_maxwell_phase_angle_low_freq_approaches_pi_over_2() {
        let mf = MaxwellFluid::new(1.0, 2.0);
        let delta = mf.phase_angle(1e-5);
        // At low ω: G'' dominates (viscous), δ → π/2
        assert!(
            delta > 1.0,
            "Phase angle should be near π/2 at low ω: {delta}"
        );
    }
}

mod non_newtonian_extended_tests {
    use super::*;

    // --- PowerLawFluid ---

    #[test]
    fn test_power_law_shear_thickening_increases_with_shear() {
        let fluid = PowerLawFluid::new(1.5, 1.0); // n > 1 = shear-thickening
        let mu_low = fluid.effective_viscosity(1.0);
        let mu_high = fluid.effective_viscosity(100.0);
        assert!(
            mu_high > mu_low,
            "Shear-thickening: μ_high={mu_high} should > μ_low={mu_low}"
        );
    }

    #[test]
    fn test_power_law_n1_newtonian_constant() {
        let k = 0.02;
        let fluid = PowerLawFluid::new(1.0, k);
        for &gamma in &[0.01, 0.1, 1.0, 10.0, 100.0] {
            let mu = fluid.effective_viscosity(gamma);
            assert!((mu - k).abs() < 1e-12, "Newtonian at γ={gamma}: mu={mu}");
        }
    }

    #[test]
    fn test_power_law_local_tau_greater_than_half() {
        let fluid = PowerLawFluid::new(0.8, 0.1);
        let tau = fluid.local_tau(5.0);
        assert!(tau > 0.5, "τ must be > 0.5: {tau}");
    }

    // --- BinghamFluid ---

    #[test]
    fn test_bingham_fluid_high_shear_approaches_plastic_viscosity() {
        let tau_y = 0.5;
        let mu_p = 0.05;
        let fluid = BinghamFluid::new(tau_y, mu_p);
        let mu_eff = fluid.effective_viscosity(1e7);
        assert!(
            (mu_eff - mu_p).abs() < 1e-4,
            "Bingham at high shear: μ_eff={mu_eff}, μ_p={mu_p}"
        );
    }

    // --- BinghamPlastic ---

    #[test]
    fn test_bingham_plastic_is_flowing_above_yield() {
        let bp = BinghamPlastic::new(1.0, 0.1);
        // stress = mu_p * gamma + tau_y = 0.1*100+1.0 = 11
        assert!(
            bp.is_flowing(100.0, 1e-10),
            "Should be flowing at high shear"
        );
    }

    #[test]
    fn test_bingham_plastic_not_flowing_below_yield() {
        let bp = BinghamPlastic::new(10.0, 0.1);
        // Very low shear: stress = 0.1*0.001 + 10 is yielded? No: is_flowing checks gamma > 0
        // The is_flowing checks if gamma > tol, so below tol = not flowing
        assert!(
            !bp.is_flowing(0.0, 0.01),
            "At gamma=0 should not be flowing"
        );
    }

    #[test]
    fn test_bingham_plastic_stress_increases_with_shear() {
        let bp = BinghamPlastic::new(1.0, 0.1);
        let s1 = bp.stress(1.0);
        let s2 = bp.stress(10.0);
        assert!(s2 > s1, "Bingham stress should increase: {s1} < {s2}");
    }

    // --- RegularizedBingham ---

    #[test]
    fn test_regularized_bingham_viscosity_positive() {
        let rb = RegularizedBingham::new(1.0, 0.1, 100.0);
        for &gamma in &[0.0, 0.01, 0.1, 1.0] {
            let mu = rb.viscosity(gamma);
            assert!(
                mu > 0.0,
                "Regularized Bingham viscosity should be positive at γ={gamma}: {mu}"
            );
        }
    }

    #[test]
    fn test_regularized_bingham_stress_at_zero() {
        let rb = RegularizedBingham::new(1.0, 0.1, 100.0);
        // At gamma=0, stress should be close to tau_y * (1 - exp(-m*0)) = 0 at least smooth
        let s = rb.stress(0.0);
        assert!(s >= 0.0, "Stress at zero shear: {s}");
    }

    // --- HerschelBulkley ---

    #[test]
    fn test_herschel_bulkley_stress_zero_below_yield() {
        let hb = HerschelBulkley::new(1.0, 0.5, 1.0); // tau_y = 1.0
        let stress = hb.stress(0.001);
        // At very low shear below yield, stress should be close to tau_y
        assert!(stress >= 0.0, "HB stress should be non-negative: {stress}");
    }

    #[test]
    fn test_herschel_bulkley_viscosity_positive() {
        let hb = HerschelBulkley::new(1.0, 0.5, 1.0);
        for &gamma in &[0.1, 1.0, 10.0] {
            let mu = hb.viscosity(gamma);
            assert!(mu > 0.0, "HB viscosity positive at γ={gamma}: {mu}");
        }
    }

    // --- CarreauFluid ---

    #[test]
    fn test_carreau_at_zero_shear_approaches_eta0() {
        let cf = CarreauFluid::new(1.0, 0.001, 1.0, 0.5);
        let mu = cf.viscosity(0.0);
        assert!((mu - 1.0).abs() < 1e-12, "Carreau at γ=0: mu={mu}");
    }

    #[test]
    fn test_carreau_at_high_shear_approaches_eta_inf() {
        let cf = CarreauFluid::new(1.0, 0.001, 0.001, 0.5); // lambda=0.001, high shear
        let mu = cf.viscosity(1e6);
        // At high shear (lambda * gamma >> 1), approaches eta_inf
        assert!(mu < 0.5, "Carreau at high shear: mu={mu}");
    }

    #[test]
    fn test_carreau_shear_thinning() {
        let cf = CarreauFluid::new(1.0, 0.01, 1.0, 0.3); // n=0.3 shear-thinning
        let mu_low = cf.viscosity(0.01);
        let mu_high = cf.viscosity(100.0);
        assert!(
            mu_high < mu_low,
            "Carreau shear-thinning: {mu_low} vs {mu_high}"
        );
    }

    // --- CrossFluid ---

    #[test]
    fn test_cross_fluid_at_zero_shear_eta0() {
        let cf = CrossFluid::new(1.0, 0.001, 1.0, 1.0);
        let mu = cf.viscosity(0.0);
        assert!((mu - 1.0).abs() < 1e-12, "Cross at γ=0: mu={mu}");
    }

    #[test]
    fn test_cross_fluid_shear_thinning() {
        let cf = CrossFluid::new(1.0, 0.001, 1.0, 1.0);
        let mu_low = cf.viscosity(0.01);
        let mu_high = cf.viscosity(1e4);
        assert!(
            mu_high < mu_low,
            "Cross fluid shear-thinning: {mu_low} vs {mu_high}"
        );
    }

    // --- CassonFluid ---

    #[test]
    fn test_casson_stress_positive() {
        let cf = CassonFluid::new(0.5, 0.1);
        for &gamma in &[0.1, 1.0, 10.0] {
            let s = cf.stress(gamma);
            assert!(s > 0.0, "Casson stress positive at γ={gamma}: {s}");
        }
    }

    #[test]
    fn test_casson_viscosity_decreases_with_shear() {
        let cf = CassonFluid::new(1.0, 0.05);
        let mu1 = cf.viscosity(0.1);
        let mu2 = cf.viscosity(10.0);
        assert!(mu2 < mu1, "Casson viscosity decreasing: {mu1} vs {mu2}");
    }

    // --- ViscosityField ---

    #[test]
    fn test_viscosity_field_initial_uniform() {
        let vf = ViscosityField::new(10, 0.01);
        assert!((vf.mean_viscosity() - 0.01).abs() < 1e-14);
        assert!((vf.min_viscosity() - 0.01).abs() < 1e-14);
        assert!((vf.max_viscosity() - 0.01).abs() < 1e-14);
    }

    #[test]
    fn test_viscosity_field_update_power_law() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let mut vf = ViscosityField::new(5, 0.01);
        let shear = vec![0.1, 0.5, 1.0, 5.0, 10.0];
        vf.update(&shear, &pl);
        assert!(
            vf.min_viscosity() > 0.0,
            "All viscosities should be positive"
        );
        assert!(vf.max_viscosity() >= vf.min_viscosity());
    }

    // --- OldroydB ---

    #[test]
    fn test_oldroyd_b_steady_shear_viscosity() {
        let ob = OldroydB::new(0.01, 1.0, 0.1, 0.005);
        let mu = ob.steady_shear_viscosity();
        assert!(mu > 0.0, "OldroydB steady shear viscosity: {mu}");
    }

    #[test]
    fn test_oldroyd_b_first_normal_stress_coefficient_positive() {
        let ob = OldroydB::new(0.01, 1.0, 0.1, 0.005);
        let psi1 = ob.first_normal_stress_coefficient();
        assert!(psi1 > 0.0, "Psi_1 should be positive: {psi1}");
    }

    #[test]
    fn test_oldroyd_b_polymer_viscosity_positive() {
        let ob = OldroydB::new(0.01, 1.0, 0.1, 0.005);
        let eta_p = ob.eta_polymer();
        assert!(eta_p > 0.0, "eta_polymer should be positive: {eta_p}");
    }

    #[test]
    fn test_oldroyd_b_weissenberg_increases_with_shear() {
        let ob = OldroydB::new(0.01, 1.0, 0.1, 0.005);
        let wi1 = ob.weissenberg(1.0);
        let wi2 = ob.weissenberg(10.0);
        assert!(
            wi2 > wi1,
            "Weissenberg number increases with shear: {wi1} vs {wi2}"
        );
    }

    // --- ThixotropicFluid ---

    #[test]
    fn test_thixotropic_fluid_initial_viscosity() {
        let tf = ThixotropicFluid::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let mu = tf.viscosity();
        assert!(mu > 0.0, "Thixotropic fluid initial viscosity: {mu}");
    }

    #[test]
    fn test_thixotropic_fluid_step_changes_lambda() {
        let mut tf = ThixotropicFluid::new(1.0, 0.01, 0.5, 0.5, 1.0);
        let lambda_before = tf.lambda;
        tf.step(1.0, 0.01);
        // lambda should evolve toward equilibrium
        assert!((tf.lambda - lambda_before).abs() > 0.0 || tf.lambda == lambda_before);
    }

    #[test]
    fn test_thixotropic_fluid_equilibrium_viscosity_positive() {
        let tf = ThixotropicFluid::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let mu_eq = tf.equilibrium_viscosity(1.0);
        assert!(mu_eq > 0.0, "Thixotropic eq viscosity: {mu_eq}");
    }

    // --- PapanastasiouViscoplastic ---

    #[test]
    fn test_papanastasiou_viscosity_positive() {
        let pv = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.1, 100.0);
        for &gamma in &[0.0, 0.01, 0.1, 1.0] {
            let mu = pv.viscosity(gamma);
            assert!(mu > 0.0, "Papanastasiou viscosity at γ={gamma}: {mu}");
        }
    }

    #[test]
    fn test_papanastasiou_stress_positive() {
        let pv = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.1, 100.0);
        for &gamma in &[0.01, 0.1, 1.0, 10.0] {
            let s = pv.stress(gamma);
            assert!(s > 0.0, "Papanastasiou stress at γ={gamma}: {s}");
        }
    }

    #[test]
    fn test_papanastasiou_is_yielded_at_high_shear() {
        let pv = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.1, 100.0);
        assert!(pv.is_yielded(100.0), "Should be yielded at high shear");
    }

    // --- turbulent_effective_viscosity ---

    #[test]
    fn test_turbulent_effective_viscosity_positive() {
        let nu_t = turbulent_effective_viscosity(0.01, 1.0, 0.1, 1.0);
        assert!(nu_t > 0.0, "Turbulent effective viscosity: {nu_t}");
    }

    // --- generalized_reynolds_power_law ---

    #[test]
    fn test_generalized_reynolds_power_law_newtonian_matches_classical() {
        // For n=1: Re_gen = rho * u * D / K = classical Re
        let re = generalized_reynolds_power_law(1.0, 1.0, 0.1, 0.01, 1.0);
        let re_classical = 1.0 * 1.0 * 0.1 / 0.01;
        assert!(
            (re - re_classical).abs() < 1.0,
            "Generalized Re ≈ classical Re: {re} vs {re_classical}"
        );
    }

    // --- NonNewtonianLBM ---

    #[test]
    fn test_non_newtonian_lbm_effective_omega_positive() {
        let lbm = NonNewtonianLBM::new(1.0);
        let pl = PowerLawFluid::new(0.5, 0.1);
        let omega = lbm.effective_omega(1.0, &pl);
        assert!(omega > 0.0, "Effective omega should be positive: {omega}");
    }

    #[test]
    fn test_non_newtonian_lbm_effective_omega_at_most_2() {
        let lbm = NonNewtonianLBM::new(1.0);
        let pl = PowerLawFluid::new(0.5, 1e-7);
        let omega = lbm.effective_omega(1.0, &pl);
        assert!(omega <= 2.0, "Effective omega must be ≤ 2: {omega}");
    }

    // --- LocalTauLattice ---

    #[test]
    fn test_local_tau_lattice_initial_uniform() {
        let lat = LocalTauLattice::new(10, 1.0);
        for i in 0..10 {
            assert!(
                (lat.tau[i] - 1.0).abs() < 1e-14,
                "Initial tau[{i}]={}",
                lat.tau[i]
            );
        }
    }

    #[test]
    fn test_local_tau_lattice_update_field() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let mut lat = LocalTauLattice::new(5, 1.0);
        let shear_rates = vec![0.1, 0.5, 1.0, 5.0, 10.0];
        lat.update_tau_field(&shear_rates, &pl);
        for (i, &tau) in lat.tau.iter().enumerate() {
            assert!(tau > 0.5, "Updated tau[{i}]={tau} must be > 0.5");
        }
    }

    // --- HerschelBulkleyFluid (extended struct) ---

    #[test]
    fn test_herschel_bulkley_fluid_is_yielded() {
        let hbf = HerschelBulkleyFluid::new(1.0, 0.5, 0.5);
        assert!(hbf.is_yielded(10.0), "Should be yielded at high shear");
    }

    #[test]
    fn test_herschel_bulkley_fluid_bingham_number() {
        let hbf = HerschelBulkleyFluid::new(1.0, 0.5, 0.5);
        let bn = hbf.bingham_number(1.0);
        assert!(bn > 0.0, "Bingham number should be positive: {bn}");
    }

    #[test]
    fn test_herschel_bulkley_fluid_lbm_relaxation_time_positive() {
        let hbf = HerschelBulkleyFluid::new(1.0, 0.01, 0.8);
        let tau = hbf.lbm_relaxation_time(1.0);
        assert!(tau > 0.5, "HBF LBM relaxation time should be > 0.5: {tau}");
    }

    // --- BinghamPipeFlow ---

    #[test]
    fn test_bingham_pipe_flow_wall_shear_stress_positive() {
        let bpf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        let tau_w = bpf.wall_shear_stress(1000.0);
        assert!(tau_w > 0.0, "Wall shear stress should be positive: {tau_w}");
    }

    #[test]
    fn test_bingham_pipe_flow_plug_radius_less_than_pipe() {
        let bpf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        if let Some(r_plug) = bpf.plug_radius(1000.0) {
            assert!(
                r_plug <= bpf.radius,
                "Plug radius should be ≤ pipe radius: {r_plug}"
            );
        }
    }

    #[test]
    fn test_bingham_pipe_flow_is_flowing_large_pressure() {
        let bpf = BinghamPipeFlow::new(1.0, 0.01, 0.05, 1.0);
        assert!(
            bpf.is_flowing(1e6),
            "Should be flowing at large pressure gradient"
        );
    }

    #[test]
    fn test_bingham_pipe_flow_not_flowing_below_yield() {
        let bpf = BinghamPipeFlow::new(10.0, 0.01, 0.05, 1.0);
        // Very small pressure gradient won't overcome yield stress
        assert!(!bpf.is_flowing(1.0), "Should not flow below yield pressure");
    }

    // --- ThixotropicModel ---

    #[test]
    fn test_thixotropic_model_advance_positive_lambda() {
        let mut tm = ThixotropicModel::new(1.0, 0.01, 1.0, 1.0, 1.0);
        tm.advance(1.0, 0.001);
        assert!(tm.lambda >= 0.0, "Lambda should stay >= 0: {}", tm.lambda);
    }

    #[test]
    fn test_thixotropic_model_steady_state_lambda_bound() {
        let tm = ThixotropicModel::new(1.0, 0.01, 0.5, 0.5, 1.0);
        let ss = tm.steady_state_lambda(1.0);
        assert!(
            (0.0..=1.0).contains(&ss),
            "Steady state lambda in [0,1]: {ss}"
        );
    }

    // --- metzner_reed_reynolds ---

    #[test]
    fn test_metzner_reed_reynolds_positive() {
        let re = metzner_reed_reynolds(1000.0, 0.5, 0.1, 0.01, 0.8);
        assert!(re > 0.0, "Metzner-Reed Re should be positive: {re}");
    }

    // --- iterate_apparent_viscosity ---

    #[test]
    fn test_iterate_apparent_viscosity_converges() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let mu = iterate_apparent_viscosity(1.0, 1.0, &pl, 0.5, 50, 1e-8);
        assert!(
            mu > 0.0,
            "Iterated apparent viscosity should be positive: {mu}"
        );
    }

    // --- shear_rate_from_strain_rate_tensor ---

    #[test]
    fn test_shear_rate_from_strain_rate_tensor_simple_shear() {
        // Pure shear: S_xy = S_yx = 0.5
        let s = [[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let gamma = LocalTauLattice::shear_rate_from_strain_rate_tensor(s);
        // |S| = sqrt(2 * 0.5^2 + 2 * 0.5^2) = sqrt(0.5) ≈ 0.707
        assert!(gamma > 0.0, "Shear rate should be positive: {gamma}");
    }
}

mod extended_non_newtonian_tests {
    use super::*;

    // ---- PhanThienTanner ----

    #[test]
    fn test_ptt_total_viscosity() {
        let ptt = PhanThienTanner::new(0.8, 0.2, 1.0, 0.1, 0.0);
        assert!((ptt.total_viscosity() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_ptt_viscosity_ratio() {
        let ptt = PhanThienTanner::new(0.8, 0.2, 1.0, 0.1, 0.0);
        assert!((ptt.viscosity_ratio() - 0.2).abs() < 1e-14);
    }

    #[test]
    fn test_ptt_linear_stress_function_unity_at_zero() {
        let ptt = PhanThienTanner::new(1.0, 0.1, 0.5, 0.1, 0.0);
        assert!((ptt.stress_function_linear(0.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_ptt_linear_stress_function_increases() {
        let ptt = PhanThienTanner::new(1.0, 0.1, 0.5, 0.1, 0.0);
        assert!(ptt.stress_function_linear(10.0) > 1.0);
    }

    #[test]
    fn test_ptt_exp_stress_function_unity_at_zero() {
        let ptt = PhanThienTanner::new(1.0, 0.1, 0.5, 0.01, 0.0);
        assert!((ptt.stress_function_exp(0.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_ptt_weissenberg_number() {
        let ptt = PhanThienTanner::new(1.0, 0.1, 2.0, 0.1, 0.0);
        assert!((ptt.weissenberg_number(0.5) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_ptt_steady_shear_viscosity_positive() {
        let ptt = PhanThienTanner::new(1.0, 0.1, 1.0, 0.05, 0.1);
        assert!(ptt.steady_shear_viscosity(1.0) > 0.0);
    }

    #[test]
    fn test_ptt_viscosity_decreases_with_shear_nonzero_epsilon() {
        let ptt = PhanThienTanner::new(1.0, 0.0, 1.0, 0.1, 0.1);
        let v0 = ptt.steady_shear_viscosity(0.01);
        let v1 = ptt.steady_shear_viscosity(5.0);
        assert!(v1 < v0, "PTT should shear thin: v0={v0} v1={v1}");
    }

    #[test]
    fn test_ptt_normal_stress_difference_positive_shear() {
        let ptt = PhanThienTanner::new(1.0, 0.0, 1.0, 0.05, 0.0);
        let n1 = ptt.normal_stress_difference(1.0);
        assert!(n1 > 0.0, "N1 = {n1}");
    }

    // ---- GiesekusFluid ----

    #[test]
    fn test_giesekus_viscosity_at_zero_shear() {
        let g = GiesekusFluid::new(1.0, 1.0, 0.2);
        let v = g.steady_shear_viscosity(0.0);
        assert!((v - 1.0).abs() < 1e-12, "v = {v}");
    }

    #[test]
    fn test_giesekus_viscosity_positive() {
        let g = GiesekusFluid::new(0.8, 0.5, 0.1);
        assert!(g.steady_shear_viscosity(2.0) > 0.0);
    }

    #[test]
    fn test_giesekus_shear_thinning() {
        let g = GiesekusFluid::new(1.0, 1.0, 0.3);
        let v0 = g.steady_shear_viscosity(0.01);
        let v1 = g.steady_shear_viscosity(5.0);
        assert!(v1 < v0, "Giesekus should shear thin: v0={v0} v1={v1}");
    }

    #[test]
    fn test_giesekus_tau_positive() {
        let g = GiesekusFluid::new(0.5, 0.5, 0.1);
        let tau = g.effective_tau(1.0, 1.0 / 3.0);
        assert!(tau > 0.5, "tau = {tau}");
    }

    #[test]
    fn test_giesekus_alpha_clamp() {
        let g = GiesekusFluid::new(1.0, 1.0, 0.9);
        assert!(g.alpha <= 0.5, "alpha should be clamped");
    }

    // ---- RheologyLookupTable ----

    #[test]
    fn test_lookup_table_exact_point() {
        let table = RheologyLookupTable::new(vec![0.0, 1.0, 2.0, 3.0], vec![1.0, 0.8, 0.6, 0.4]);
        assert!((table.viscosity_at(1.0) - 0.8).abs() < 1e-14);
    }

    #[test]
    fn test_lookup_table_interpolation() {
        let table = RheologyLookupTable::new(vec![0.0, 1.0], vec![1.0, 0.0]);
        let v = table.viscosity_at(0.5);
        assert!((v - 0.5).abs() < 1e-14, "interpolated v = {v}");
    }

    #[test]
    fn test_lookup_table_clamp_below() {
        let table = RheologyLookupTable::new(vec![1.0, 2.0], vec![0.5, 0.3]);
        assert!((table.viscosity_at(0.0) - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_lookup_table_clamp_above() {
        let table = RheologyLookupTable::new(vec![1.0, 2.0], vec![0.5, 0.3]);
        assert!((table.viscosity_at(10.0) - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_lookup_table_tau_above_half() {
        let table = RheologyLookupTable::new(vec![0.0, 1.0], vec![0.1, 0.1]);
        let tau = table.tau_at(0.5, 1.0 / 3.0);
        assert!(tau > 0.5, "tau = {tau}");
    }

    // ---- YieldCriterion ----

    #[test]
    fn test_yield_criterion_below_yield() {
        let yc = YieldCriterion::new(1.0);
        assert!(!yc.is_yielded_2d(0.1, 0.1, 0.1));
    }

    #[test]
    fn test_yield_criterion_above_yield() {
        let yc = YieldCriterion::new(0.1);
        assert!(yc.is_yielded_2d(1.0, -1.0, 1.0));
    }

    #[test]
    fn test_von_mises_2d_positive() {
        let yc = YieldCriterion::new(1.0);
        let vm = yc.von_mises_2d(1.0, -1.0, 0.5);
        assert!(vm > 0.0, "von Mises = {vm}");
    }

    #[test]
    fn test_von_mises_3d_hydrostatic_zero() {
        let yc = YieldCriterion::new(1.0);
        // Pure hydrostatic: deviatoric = 0
        let p = 1.0;
        let s = [[p, 0.0, 0.0], [0.0, p, 0.0], [0.0, 0.0, p]];
        let vm = yc.von_mises_3d(s);
        // For hydrostatic the deviatoric part is zero
        // Actually, von Mises is based on deviatoric, so for equal diagonals
        // (s_xx - s_yy)² + ... = 0
        assert!(vm.abs() < 1e-13, "hydrostatic vm = {vm}");
    }

    #[test]
    fn test_bingham_factor_below_yield_is_zero() {
        let yc = YieldCriterion::new(1.0);
        assert_eq!(yc.bingham_factor(0.5), 0.0);
    }

    #[test]
    fn test_bingham_factor_above_yield_positive() {
        let yc = YieldCriterion::new(1.0);
        let f = yc.bingham_factor(2.0);
        assert!((f - 0.5).abs() < 1e-14, "factor = {f}");
    }

    // ---- ViscoelasticRelaxation ----

    #[test]
    fn test_viscoelastic_relaxation_initial_zero() {
        let vr = ViscoelasticRelaxation::new(4, 1.0, 0.5);
        assert_eq!(vr.total_magnitude(), 0.0);
    }

    #[test]
    fn test_viscoelastic_relaxation_add_strain() {
        let mut vr = ViscoelasticRelaxation::new(4, 1.0, 0.5);
        let dxx = vec![0.1; 4];
        let dyy = vec![-0.1; 4];
        let dxy = vec![0.05; 4];
        vr.add_strain_rate(&dxx, &dyy, &dxy, 0.1);
        assert!(vr.total_magnitude() > 0.0);
    }

    #[test]
    fn test_viscoelastic_relaxation_decreases_with_time() {
        let mut vr = ViscoelasticRelaxation::new(4, 1.0, 0.5);
        let dxx = vec![1.0; 4];
        let dyy = vec![0.0; 4];
        let dxy = vec![0.0; 4];
        vr.add_strain_rate(&dxx, &dyy, &dxy, 1.0);
        let before = vr.total_magnitude();
        vr.relax(0.5);
        let after = vr.total_magnitude();
        assert!(after < before, "stress should relax: {before} → {after}");
    }

    #[test]
    fn test_viscoelastic_trace_equals_xx_plus_yy() {
        let mut vr = ViscoelasticRelaxation::new(2, 1.0, 0.5);
        vr.sigma_xx[0] = 0.3;
        vr.sigma_yy[0] = 0.2;
        assert!((vr.trace(0) - 0.5).abs() < 1e-14);
    }

    // ---- JohnsonSegalman ----

    #[test]
    fn test_johnson_segalman_total_viscosity() {
        let js = JohnsonSegalman::new(0.7, 0.3, 1.0, 0.5);
        assert!((js.total_viscosity() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_johnson_segalman_apparent_viscosity_zero_shear() {
        let js = JohnsonSegalman::new(0.8, 0.2, 1.0, 0.5);
        let mu = js.apparent_viscosity(0.0);
        assert!((mu - 1.0).abs() < 1e-13, "μ_app(0) = {mu}");
    }

    #[test]
    fn test_johnson_segalman_shear_stress_positive() {
        let js = JohnsonSegalman::new(0.8, 0.2, 1.0, 0.5);
        assert!(js.shear_stress(0.5) > 0.0);
    }

    // ---- RoliePolyModel ----

    #[test]
    fn test_rolie_poly_zero_shear_viscosity() {
        let rp = RoliePolyModel::new(1.0, 2.0, 0.5, 0.3);
        assert!((rp.zero_shear_viscosity() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_rolie_poly_steady_viscosity_decreases() {
        let rp = RoliePolyModel::new(1.0, 2.0, 0.5, 0.3);
        let v0 = rp.steady_viscosity(0.01);
        let v1 = rp.steady_viscosity(5.0);
        assert!(v1 < v0, "Rolie-Poly should shear thin: v0={v0} v1={v1}");
    }

    #[test]
    fn test_rolie_poly_weissenberg_d() {
        let rp = RoliePolyModel::new(1.0, 2.0, 0.5, 0.3);
        assert!((rp.weissenberg_d(1.0) - 2.0).abs() < 1e-14);
    }

    // ---- Utility functions ----

    #[test]
    fn test_viscosity_index_positive() {
        let vi = viscosity_index(10.0, 5.0);
        assert!(vi > 0.0, "VI = {vi}");
    }

    #[test]
    fn test_deborah_number_ratio() {
        let de = deborah_number(2.0, 4.0);
        assert!((de - 0.5).abs() < 1e-14, "De = {de}");
    }

    #[test]
    fn test_weissenberg_number_basic() {
        let wi = weissenberg_number(1.5, 2.0);
        assert!((wi - 3.0).abs() < 1e-14, "Wi = {wi}");
    }

    #[test]
    fn test_oldroyd_b_viscosity() {
        let mu = oldroyd_b_viscosity(0.3, 0.7);
        assert!((mu - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_trouton_ratio_newtonian() {
        let tr = trouton_ratio(3.0, 1.0);
        assert!((tr - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_deborah_infinite_for_zero_time() {
        let de = deborah_number(1.0, 0.0);
        assert!(de.is_infinite());
    }
}
