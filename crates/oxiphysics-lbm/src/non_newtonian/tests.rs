// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for non-Newtonian fluid models.

use super::*;

#[allow(clippy::module_inception)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Power-law: n=1, K=μ → constant Newtonian viscosity
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_water_n1() {
        let mu = 1e-3;
        let fluid = PowerLawFluid::new(1.0, mu);
        for &gamma in &[0.01, 0.1, 1.0, 10.0, 100.0] {
            let mu_eff = fluid.effective_viscosity(gamma);
            assert!(
                (mu_eff - mu).abs() < 1e-12,
                "n=1 fluid: expected μ={mu}, got {mu_eff} at γ̇={gamma}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Power-law shear-thinning: n=0.5 → higher shear → lower viscosity
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_shear_thinning() {
        let fluid = PowerLawFluid::new(0.5, 1.0);
        let mu_low = fluid.effective_viscosity(1.0);
        let mu_high = fluid.effective_viscosity(100.0);
        assert!(
            mu_high < mu_low,
            "Shear-thinning: μ at high shear ({mu_high}) should be < μ at low shear ({mu_low})"
        );
    }

    // -----------------------------------------------------------------------
    // Bingham: shear below yield → very large (rigid) effective viscosity
    // -----------------------------------------------------------------------
    #[test]
    fn test_bingham_below_yield() {
        let tau_y = 1.0;
        let mu_p = 0.01;
        let fluid = BinghamFluid::new(tau_y, mu_p);
        // γ̇ = 0.001, so μ_p * γ̇ = 1e-5 << τ_y  → rigid
        let mu_eff = fluid.effective_viscosity(0.001);
        assert!(
            mu_eff > 1.0,
            "Below yield: effective viscosity should be very large, got {mu_eff}"
        );
    }

    // -----------------------------------------------------------------------
    // Bingham: shear >> yield → approaches plastic viscosity
    // -----------------------------------------------------------------------
    #[test]
    fn test_bingham_above_yield() {
        let tau_y = 1.0;
        let mu_p = 0.1;
        let fluid = BinghamFluid::new(tau_y, mu_p);
        // At very high shear rate the τ_y / γ̇ term → 0.
        let gamma = 1e6;
        let mu_eff = fluid.effective_viscosity(gamma);
        // μ_eff = μ_p + τ_y / γ̇ ≈ μ_p for large γ̇
        assert!(
            (mu_eff - mu_p).abs() < 1e-4,
            "Above yield: μ_eff={mu_eff} should approach μ_p={mu_p}"
        );
    }

    // -----------------------------------------------------------------------
    // LocalTauLattice: shear_rate_from_strain_rate_tensor
    // -----------------------------------------------------------------------
    #[test]
    fn test_shear_rate_from_strain_rate_tensor_zero() {
        let s = [[0.0; 3]; 3];
        let gamma = LocalTauLattice::shear_rate_from_strain_rate_tensor(s);
        assert!(
            gamma.abs() < 1e-15,
            "Zero strain tensor should give zero shear rate, got {gamma}"
        );
    }

    // -----------------------------------------------------------------------
    // LocalTauLattice: update_tau_field applies fluid model correctly
    // -----------------------------------------------------------------------
    #[test]
    fn test_local_tau_lattice_update() {
        let n_nodes = 5;
        let mut lat = LocalTauLattice::new(n_nodes, 1.0);
        let fluid = PowerLawFluid::new(1.0, 1.0 / 6.0); // ν = 1/6 → τ = 1
        let shear_rates = vec![0.1; n_nodes];
        lat.update_tau_field(&shear_rates, &fluid);
        let expected_tau = fluid.local_tau(0.1);
        for (k, &tau_k) in lat.tau.iter().enumerate() {
            assert!(
                (tau_k - expected_tau).abs() < 1e-12,
                "Node {k}: tau={tau_k}, expected {expected_tau}"
            );
        }
    }

    // -------------------------------------------------------------------
    // PowerLawFluid: viscosity / stress / is_shear_thinning
    // -------------------------------------------------------------------

    #[test]
    fn test_power_law_stress() {
        let fluid = PowerLawFluid::new(0.5, 2.0);
        let s = fluid.stress(4.0);
        // stress = K * |gamma|^n = 2.0 * 4.0^0.5 = 2 * 2 = 4
        assert!((s - 4.0).abs() < 1e-12, "stress = {s}, expected 4.0");
    }

    #[test]
    fn test_power_law_is_shear_thinning() {
        let thin = PowerLawFluid::new(0.5, 1.0);
        let thick = PowerLawFluid::new(1.5, 1.0);
        let newton = PowerLawFluid::new(1.0, 1.0);
        assert!(thin.is_shear_thinning());
        assert!(!thick.is_shear_thinning());
        assert!(!newton.is_shear_thinning());
    }

    // -------------------------------------------------------------------
    // BinghamPlastic: stress / viscosity / is_flowing
    // -------------------------------------------------------------------

    #[test]
    fn test_bingham_stress() {
        let b = BinghamPlastic::new(1.0, 0.1);
        let s = b.stress(10.0);
        // stress = tau_y + mu_p * gamma = 1.0 + 0.1 * 10 = 2.0
        assert!((s - 2.0).abs() < 1e-12, "stress = {s}, expected 2.0");
    }

    #[test]
    fn test_bingham_is_flowing() {
        let b = BinghamPlastic::new(1.0, 0.1);
        assert!(!b.is_flowing(0.001, 1e-10));
        assert!(b.is_flowing(100.0, 1e-10));
    }

    // -------------------------------------------------------------------
    // HerschelBulkley
    // -------------------------------------------------------------------

    #[test]
    fn test_herschel_bulkley_reduces_to_power_law() {
        // With tau_y = 0, HB reduces to power law.
        let hb = HerschelBulkley::new(2.0, 0.5, 0.0);
        let pl = PowerLawFluid::new(0.5, 2.0);
        let gamma = 4.0;
        let s_hb = hb.stress(gamma);
        let s_pl = pl.stress(gamma);
        assert!(
            (s_hb - s_pl).abs() < 1e-12,
            "HB stress = {s_hb}, PL stress = {s_pl}"
        );
    }

    #[test]
    fn test_herschel_bulkley_viscosity() {
        let hb = HerschelBulkley::new(1.0, 1.0, 2.0);
        let v = hb.viscosity(5.0);
        // stress = tau_y + K * gamma^n = 2.0 + 1.0 * 5.0 = 7.0
        // viscosity = stress / gamma = 7.0 / 5.0 = 1.4
        assert!((v - 1.4).abs() < 1e-12, "viscosity = {v}, expected 1.4");
    }

    // -------------------------------------------------------------------
    // CarreauFluid
    // -------------------------------------------------------------------

    #[test]
    fn test_carreau_zero_shear_rate() {
        let c = CarreauFluid::new(1.0, 0.001, 2.0, 0.5);
        let v = c.viscosity(0.0);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "viscosity at 0 = {v}, expected eta_0 = 1.0"
        );
    }

    #[test]
    fn test_carreau_high_shear_rate() {
        let c = CarreauFluid::new(1.0, 0.001, 2.0, 0.5);
        let v = c.viscosity(1e10);
        // At very high shear rate, should approach eta_inf.
        assert!(
            (v - 0.001).abs() < 0.01,
            "viscosity at high shear = {v}, should approach eta_inf = 0.001"
        );
    }

    // -------------------------------------------------------------------
    // NonNewtonianLBM
    // -------------------------------------------------------------------

    #[test]
    fn test_non_newtonian_lbm_effective_omega() {
        let lbm = NonNewtonianLBM::new(1.0);
        let pl = PowerLawFluid::new(1.0, 1.0 / 6.0);
        let omega = lbm.effective_omega(1.0, &pl);
        // nu = 1/6, tau = 0.5 + nu/cs2 = 0.5 + (1/6)/(1/3) = 0.5 + 0.5 = 1.0
        // omega = 1/tau = 1.0
        assert!((omega - 1.0).abs() < 1e-12, "omega = {omega}, expected 1.0");
    }

    // -------------------------------------------------------------------
    // Cross model
    // -------------------------------------------------------------------

    #[test]
    fn test_cross_zero_shear() {
        let c = CrossFluid::new(1.0, 0.01, 2.0, 1.0);
        let v = c.viscosity(0.0);
        // At gamma=0: denominator = 1 + 0 = 1 → mu = eta_inf + (eta_0 - eta_inf) = eta_0
        assert!(
            (v - 1.0).abs() < 1e-12,
            "Cross at gamma=0: {v}, expected 1.0"
        );
    }

    #[test]
    fn test_cross_high_shear() {
        let c = CrossFluid::new(1.0, 0.01, 2.0, 1.0);
        let v = c.viscosity(1e10);
        // At very high gamma: denominator → ∞ → mu → eta_inf
        assert!(
            (v - 0.01).abs() < 0.01,
            "Cross at high shear: {v}, expected ~0.01"
        );
    }

    #[test]
    fn test_cross_shear_thinning() {
        let c = CrossFluid::new(1.0, 0.01, 2.0, 1.0);
        let v_low = c.viscosity(0.1);
        let v_high = c.viscosity(100.0);
        assert!(v_high < v_low, "Cross model should be shear-thinning");
    }

    // -------------------------------------------------------------------
    // Casson model
    // -------------------------------------------------------------------

    #[test]
    fn test_casson_stress_at_zero() {
        let c = CassonFluid::new(1.0, 0.01);
        let s = c.stress(0.0);
        // sqrt(stress) = sqrt(tau_y) + sqrt(0) = 1 → stress = 1
        assert!((s - 1.0).abs() < 1e-12, "Casson stress at 0: {s}");
    }

    #[test]
    fn test_casson_viscosity_decreases() {
        let c = CassonFluid::new(1.0, 0.01);
        let v_low = c.viscosity(1.0);
        let v_high = c.viscosity(100.0);
        assert!(
            v_high < v_low,
            "Casson: viscosity should decrease with shear rate"
        );
    }

    #[test]
    fn test_casson_stress_increases() {
        let c = CassonFluid::new(1.0, 0.01);
        let s1 = c.stress(1.0);
        let s10 = c.stress(10.0);
        assert!(s10 > s1, "Casson stress should increase with shear rate");
    }

    // -------------------------------------------------------------------
    // Regularized Bingham
    // -------------------------------------------------------------------

    #[test]
    fn test_regularized_bingham_smooth() {
        let rb = RegularizedBingham::new(1.0, 0.1, 1000.0);
        // Check viscosity is finite at gamma=0
        let v = rb.viscosity(0.0);
        assert!(
            v.is_finite(),
            "Reg. Bingham viscosity at 0 should be finite"
        );
        assert!(v > 0.0, "Reg. Bingham viscosity at 0 should be positive");
    }

    #[test]
    fn test_regularized_bingham_approaches_bingham() {
        let rb = RegularizedBingham::new(1.0, 0.1, 10000.0);
        let bf = BinghamFluid::new(1.0, 0.1);
        // At high shear rate, regularized should approach standard
        let gamma = 1000.0;
        let v_reg = rb.viscosity(gamma);
        let v_std = bf.effective_viscosity(gamma);
        assert!(
            (v_reg - v_std).abs() < 0.01,
            "Reg. Bingham should approach standard at high shear: reg={v_reg}, std={v_std}"
        );
    }

    #[test]
    fn test_regularized_bingham_stress() {
        let rb = RegularizedBingham::new(1.0, 0.1, 1000.0);
        let s = rb.stress(10.0);
        assert!(s > 0.0, "Stress should be positive");
    }

    // -------------------------------------------------------------------
    // Apparent viscosity iteration
    // -------------------------------------------------------------------

    #[test]
    fn test_iterate_apparent_viscosity_converges() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let gamma = 10.0;
        let result = iterate_apparent_viscosity(
            1.0, // initial guess
            gamma, &pl, 1.0, // full relaxation
            100, 1e-10,
        );
        let expected = pl.effective_viscosity(gamma);
        assert!(
            (result - expected).abs() < 1e-8,
            "Iteration should converge: result={result}, expected={expected}"
        );
    }

    #[test]
    fn test_iterate_with_underrelaxation() {
        let pl = PowerLawFluid::new(0.5, 1.0);
        let gamma = 10.0;
        let result = iterate_apparent_viscosity(
            0.5, // different initial guess
            gamma, &pl, 0.5, // under-relaxation
            100, 1e-10,
        );
        let expected = pl.effective_viscosity(gamma);
        assert!(
            (result - expected).abs() < 1e-6,
            "Under-relaxed iteration should converge: result={result}, expected={expected}"
        );
    }

    // -------------------------------------------------------------------
    // ViscosityField
    // -------------------------------------------------------------------

    #[test]
    fn test_viscosity_field_uniform() {
        let vf = ViscosityField::new(10, 0.5);
        assert_eq!(vf.values.len(), 10);
        assert!((vf.min_viscosity() - 0.5).abs() < 1e-14);
        assert!((vf.max_viscosity() - 0.5).abs() < 1e-14);
        assert!((vf.mean_viscosity() - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_viscosity_field_update() {
        let mut vf = ViscosityField::new(3, 1.0);
        let pl = PowerLawFluid::new(0.5, 1.0);
        let shear_rates = vec![1.0, 10.0, 100.0];
        vf.update(&shear_rates, &pl);
        // Check that viscosities differ (shear-thinning)
        assert!(vf.values[0] > vf.values[1]);
        assert!(vf.values[1] > vf.values[2]);
    }

    #[test]
    fn test_viscosity_field_min_max() {
        let mut vf = ViscosityField::new(3, 1.0);
        vf.values = vec![0.1, 0.5, 0.3];
        assert!((vf.min_viscosity() - 0.1).abs() < 1e-14);
        assert!((vf.max_viscosity() - 0.5).abs() < 1e-14);
    }

    // ── Oldroyd-B tests ────────────────────────────────────────────────────

    #[test]
    fn test_oldroyd_b_polymer_viscosity() {
        let ob = OldroydB::new(0.1, 1.0, 0.1, 0.01);
        // eta_p = eta - eta_s = 0.1 - 0.01 = 0.09
        assert!(
            (ob.eta_polymer() - 0.09).abs() < 1e-12,
            "eta_p = {}",
            ob.eta_polymer()
        );
    }

    #[test]
    fn test_oldroyd_b_weissenberg_number() {
        let ob = OldroydB::new(0.1, 2.0, 0.1, 0.01);
        let wi = ob.weissenberg(5.0);
        // Wi = lambda_1 * gamma = 2 * 5 = 10
        assert!((wi - 10.0).abs() < 1e-12, "Wi = {wi}");
    }

    #[test]
    fn test_oldroyd_b_steady_shear_viscosity_constant() {
        let ob = OldroydB::new(0.15, 1.0, 0.1, 0.05);
        // Oldroyd-B has constant steady-shear viscosity = eta
        assert!((ob.steady_shear_viscosity() - 0.15).abs() < 1e-12);
    }

    #[test]
    fn test_oldroyd_b_first_normal_stress_coefficient() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        // Psi_1 = 2 * eta_p * lambda_1 = 2 * 0.09 * 1.0 = 0.18
        let psi1 = ob.first_normal_stress_coefficient();
        assert!((psi1 - 0.18).abs() < 1e-12, "Psi_1 = {psi1}");
    }

    #[test]
    fn test_oldroyd_b_first_normal_stress_difference() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        let n1 = ob.first_normal_stress_difference(2.0);
        let expected = ob.first_normal_stress_coefficient() * 4.0;
        assert!((n1 - expected).abs() < 1e-12, "N1 = {n1}");
    }

    #[test]
    fn test_oldroyd_b_relaxation_modulus_at_zero() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        let g0 = ob.relaxation_modulus(0.0);
        let expected = ob.eta_polymer() / ob.lambda_1;
        assert!((g0 - expected).abs() < 1e-12, "G(0) = {g0}");
    }

    #[test]
    fn test_oldroyd_b_relaxation_modulus_decays() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        let g0 = ob.relaxation_modulus(0.0);
        let g1 = ob.relaxation_modulus(1.0);
        assert!(g1 < g0, "G(t) should decay: G(0)={g0}, G(1)={g1}");
    }

    #[test]
    fn test_oldroyd_b_storage_loss_moduli() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        let gp = ob.storage_modulus(1.0);
        let gpp = ob.loss_modulus(1.0);
        assert!(gp > 0.0, "G' > 0: {gp}");
        assert!(gpp > 0.0, "G'' > 0: {gpp}");
    }

    #[test]
    fn test_oldroyd_b_conformation_tensor_equilibrium() {
        let ob = OldroydB::new(0.1, 1.0, 0.0, 0.01);
        // Advance many steps with zero shear rate → should approach identity tensor
        let (mut axx, mut axy, mut ayy): (f64, f64, f64) = (2.0, 0.5, 1.5);
        for _ in 0..20000 {
            let (a, b, c) = ob.step_conformation_tensor(axx, axy, ayy, 0.0, 1e-3);
            axx = a;
            axy = b;
            ayy = c;
        }
        assert!((axx - 1.0).abs() < 1e-3, "Axx at eq: {axx}");
        assert!((axy - 0.0).abs() < 1e-3, "Axy at eq: {axy}");
        assert!((ayy - 1.0).abs() < 1e-3, "Ayy at eq: {ayy}");
    }

    // ── Thixotropic model tests ────────────────────────────────────────────

    #[test]
    fn test_thixotropic_initial_viscosity() {
        let f = ThixotropicFluid::new(1.0, 0.01, 1.0, 1.0, 1.0);
        // lambda = 1 → mu = mu_inf + (mu_0 - mu_inf) * 1 = mu_0
        assert!(
            (f.viscosity() - 1.0).abs() < 1e-12,
            "Initial viscosity = {}",
            f.viscosity()
        );
    }

    #[test]
    fn test_thixotropic_equilibrium_lambda() {
        let f = ThixotropicFluid::new(1.0, 0.01, 1.0, 1.0, 1.0);
        let leq = f.equilibrium_lambda(1.0);
        // lambda_eq = a / (a + b * gamma^m) = 1 / (1 + 1) = 0.5
        assert!((leq - 0.5).abs() < 1e-12, "lambda_eq = {leq}");
    }

    #[test]
    fn test_thixotropic_step_decreases_lambda_under_shear() {
        let mut f = ThixotropicFluid::new(1.0, 0.01, 0.1, 10.0, 1.0);
        let l_before = f.lambda;
        f.step(1.0, 0.1); // high shear → breakdown
        assert!(f.lambda < l_before, "Lambda should decrease under shear");
    }

    #[test]
    fn test_thixotropic_step_recovers_at_rest() {
        let mut f = ThixotropicFluid::new(1.0, 0.01, 1.0, 1.0, 1.0);
        f.lambda = 0.1; // start broken
        f.step(0.0, 0.1); // no shear → recovery
        assert!(f.lambda > 0.1, "Lambda should recover at rest");
    }

    #[test]
    fn test_thixotropic_time_scale() {
        let f = ThixotropicFluid::new(1.0, 0.01, 2.0, 1.0, 1.0);
        assert!(
            (f.time_scale() - 0.5).abs() < 1e-12,
            "Time scale = {}",
            f.time_scale()
        );
    }

    #[test]
    fn test_thixotropic_lambda_stays_bounded() {
        let mut f = ThixotropicFluid::new(1.0, 0.01, 1.0, 100.0, 1.0);
        for _ in 0..1000 {
            f.step(10.0, 0.01);
        }
        assert!(
            f.lambda >= 0.0 && f.lambda <= 1.0,
            "Lambda out of bounds: {}",
            f.lambda
        );
    }

    // ── Papanastasiou viscoplastic tests ──────────────────────────────────

    #[test]
    fn test_papanastasiou_viscosity_finite_at_zero() {
        let p = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.01, 1000.0);
        let v = p.viscosity(0.0);
        assert!(
            v.is_finite() && v > 0.0,
            "Papanastasiou viscosity at 0: {v}"
        );
    }

    #[test]
    fn test_papanastasiou_approaches_newtonian_high_shear() {
        let p = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.1, 1000.0);
        let v_high = p.viscosity(1e8);
        // At very high gamma: yield_term → tau_y / gamma → 0, mu → mu_inf
        assert!(
            (v_high - 0.1).abs() < 0.01,
            "High shear viscosity → mu_inf: {v_high}"
        );
    }

    #[test]
    fn test_papanastasiou_stress_positive() {
        let p = PapanastasiouViscoplastic::full(1.0, 0.01, 0.5, 0.8, 100.0);
        let s = p.stress(5.0);
        assert!(s > 0.0, "Papanastasiou stress should be positive: {s}");
    }

    #[test]
    fn test_papanastasiou_is_yielded() {
        let p = PapanastasiouViscoplastic::papanastasiou_only(1.0, 0.01, 100.0);
        // At gamma = 0, 1/m = 0.01 → not yielded
        assert!(!p.is_yielded(0.0), "gamma=0 should be un-yielded");
        // At gamma = 1.0 >> 0.01 = 1/m → yielded
        assert!(p.is_yielded(1.0), "gamma=1 should be yielded");
    }

    // ── Non-Newtonian turbulence tests ────────────────────────────────────

    #[test]
    fn test_turbulent_effective_viscosity() {
        let mu_total = turbulent_effective_viscosity(0.001, 1.0, 0.1, 1.0);
        // mu_lam + rho * l^2 * S = 0.001 + 1 * 0.01 * 1 = 0.011
        assert!((mu_total - 0.011).abs() < 1e-12, "mu_total = {mu_total}");
    }

    #[test]
    fn test_turbulent_effective_viscosity_laminar_at_zero_length() {
        let mu_total = turbulent_effective_viscosity(0.001, 1.0, 0.0, 1.0);
        assert!(
            (mu_total - 0.001).abs() < 1e-12,
            "No mixing length → laminar"
        );
    }

    #[test]
    fn test_van_driest_mixing_length_near_wall() {
        // Very close to wall: y_plus is small → l_m → 0
        let l_m = van_driest_mixing_length(1e-6, 0.1, 1e-6, 0.41, 26.0);
        let l_log = van_driest_mixing_length(0.01, 0.1, 1e-6, 0.41, 26.0);
        assert!(l_m < l_log, "Mixing length should increase away from wall");
    }

    #[test]
    fn test_smagorinsky_turbulent_viscosity() {
        let nu_t = smagorinsky_turbulent_viscosity(0.1, 1.0, 2.0);
        // nu_t = (0.1 * 1)^2 * 2 = 0.01 * 2 = 0.02
        assert!((nu_t - 0.02).abs() < 1e-12, "nu_t = {nu_t}");
    }

    #[test]
    fn test_kolmogorov_scale_positive() {
        let eta_k = kolmogorov_scale(1e-6, 1e-4);
        assert!(
            eta_k > 0.0 && eta_k.is_finite(),
            "Kolmogorov scale = {eta_k}"
        );
    }

    #[test]
    fn test_kolmogorov_scale_zero_epsilon() {
        let eta_k = kolmogorov_scale(1e-6, 0.0);
        assert_eq!(eta_k, f64::INFINITY, "Zero dissipation → infinite scale");
    }

    #[test]
    fn test_generalized_reynolds_power_law() {
        let re = generalized_reynolds_power_law(1000.0, 1.0, 0.01, 0.001, 1.0);
        // n=1: Re = rho * U * L / K = 1000 * 1 * 0.01 / 0.001 / 8^0 = 10000
        assert!((re - 10000.0).abs() < 1e-6, "Re = {re}");
    }

    #[test]
    fn test_metzner_reed_reynolds_positive() {
        let re = metzner_reed_reynolds(1000.0, 1.0, 0.01, 0.001, 1.0);
        assert!(re > 0.0 && re.is_finite(), "Metzner-Reed Re = {re}");
    }
}

mod nn_additions_tests {
    use super::*;

    #[test]
    fn test_local_viscosity_model_power_law_newtonian() {
        // n=1 should give Newtonian: viscosity = K for any shear rate > 0.
        let mu = 0.01;
        let fluid = PowerLawFluid::new(1.0, mu);
        let gamma = 5.0;
        let nu = fluid.viscosity(gamma);
        assert!(
            (nu - mu).abs() < 1e-12,
            "Power-law n=1 should give nu=K={mu}, got {nu}"
        );
    }

    #[test]
    fn test_local_viscosity_model_bingham_below_yield() {
        // Below yield stress the Bingham viscosity should be RIGID_MU (very large).
        let tau_y = 1.0;
        let mu_p = 0.01;
        let fluid = BinghamFluid::new(tau_y, mu_p);
        let nu = fluid.viscosity(0.0); // zero shear rate → unyielded
        assert!(
            nu > 1.0,
            "Bingham below yield should return large viscosity, got {nu}"
        );
    }

    #[test]
    fn test_local_viscosity_model_carreau_limits() {
        // At very high shear rates Carreau should approach mu_inf.
        let mu_0 = 1.0;
        let mu_inf = 0.001;
        let lambda = 1.0;
        let n = 0.5;
        let fluid = CarreauFluid::new(mu_0, mu_inf, lambda, n);
        let nu_high = fluid.viscosity(1e6);
        assert!(
            nu_high < mu_0,
            "Carreau viscosity at high shear should be < mu_0, got {nu_high}"
        );
    }

    #[test]
    fn test_shear_rate_from_strain_tensor_zero() {
        let s = [[0.0; 3]; 3];
        let gamma = shear_rate_from_strain_tensor(s);
        assert!(
            gamma.abs() < 1e-14,
            "Zero tensor → zero shear rate, got {gamma}"
        );
    }

    #[test]
    fn test_shear_rate_from_strain_tensor_identity() {
        // S = I (identity), sum_sq = 3, result = sqrt(6).
        let s = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let gamma = shear_rate_from_strain_tensor(s);
        let expected = (2.0_f64 * 3.0).sqrt();
        assert!(
            (gamma - expected).abs() < 1e-12,
            "Identity tensor: expected sqrt(6)={expected}, got {gamma}"
        );
    }

    #[test]
    fn test_update_relaxation_field_uniform() {
        let n = 4;
        let f = vec![[0.0_f64; 19]; n];
        let nu = 1.0 / 6.0; // gives tau = 0.5 + 3*nu = 0.5 + 0.5 = 1.0, omega = 1.0
        let visc = vec![nu; n];
        let omegas = update_relaxation_field(&f, &visc, 1.0, 1.0);
        for (k, &omega) in omegas.iter().enumerate() {
            assert!(
                (omega - 1.0).abs() < 1e-12,
                "omega[{k}] = {omega}, expected 1.0"
            );
        }
    }

    #[test]
    fn test_update_relaxation_field_positive() {
        let n = 8;
        let f = vec![[0.0_f64; 19]; n];
        let visc: Vec<f64> = (0..n).map(|k| 0.01 * (k as f64 + 1.0)).collect();
        let omegas = update_relaxation_field(&f, &visc, 1.0, 1.0);
        for &omega in &omegas {
            assert!(omega > 0.0 && omega <= 2.0, "omega={omega} out of (0,2]");
        }
    }
}

mod carreau_yasuda_dispatch_tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Carreau-Yasuda: zero-shear limit → η₀
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_yasuda_zero_shear() {
        let cy = CarreauYasudaFluid::new(1.0, 0.01, 2.0, 2.0, 0.5);
        let v = cy.viscosity(0.0);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "Carreau-Yasuda at γ̇=0 should return η₀=1.0, got {v}"
        );
    }

    // -----------------------------------------------------------------------
    // Carreau-Yasuda: high shear approaches η_∞
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_yasuda_high_shear() {
        let cy = CarreauYasudaFluid::new(1.0, 0.001, 1.0, 2.0, 0.3);
        let v = cy.viscosity(1e8);
        assert!(
            (v - 0.001).abs() < 0.005,
            "Carreau-Yasuda at high shear should approach η_∞=0.001, got {v}"
        );
    }

    // -----------------------------------------------------------------------
    // Carreau-Yasuda: a=2 matches standard Carreau model
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_yasuda_a2_matches_carreau() {
        let eta_0 = 1.0;
        let eta_inf = 0.01;
        let lambda = 2.0;
        let n = 0.5;
        let cy = CarreauYasudaFluid::new(eta_0, eta_inf, lambda, 2.0, n);
        let cr = CarreauFluid::new(eta_0, eta_inf, lambda, n);
        for &gamma in &[0.0, 0.1, 1.0, 10.0, 100.0] {
            let v_cy = cy.viscosity(gamma);
            let v_cr = cr.viscosity(gamma);
            assert!(
                (v_cy - v_cr).abs() < 1e-12,
                "Carreau-Yasuda(a=2) should match Carreau at γ̇={gamma}: cy={v_cy}, cr={v_cr}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Carreau-Yasuda: shear-thinning (viscosity decreases with shear rate)
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_yasuda_shear_thinning() {
        let cy = CarreauYasudaFluid::new(1.0, 0.01, 1.0, 1.5, 0.4);
        let v_low = cy.viscosity(0.01);
        let v_high = cy.viscosity(100.0);
        assert!(
            v_high < v_low,
            "Carreau-Yasuda should be shear-thinning: v_low={v_low}, v_high={v_high}"
        );
    }

    // -----------------------------------------------------------------------
    // Carreau-Yasuda implements NonNewtonianFluid trait
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_yasuda_non_newtonian_trait() {
        let cy = CarreauYasudaFluid::new(1.0, 0.01, 2.0, 2.0, 0.5);
        let fluid: &dyn NonNewtonianFluid = &cy;
        let nu = fluid.effective_viscosity(1.0);
        let nu_direct = cy.viscosity(1.0);
        assert!(
            (nu - nu_direct).abs() < 1e-14,
            "Trait dispatch should match direct call"
        );
    }

    // -----------------------------------------------------------------------
    // effective_viscosity dispatch: Newtonian power-law (n=1)
    // -----------------------------------------------------------------------
    #[test]
    fn test_dispatch_power_law_newtonian() {
        let mu = 0.05_f64;
        let model = RheologyModel::PowerLaw(PowerLawFluid::new(1.0, mu));
        for &gamma in &[0.1, 1.0, 10.0] {
            let nu = effective_viscosity(model, gamma);
            assert!(
                (nu - mu).abs() < 1e-12,
                "Dispatch power-law n=1: nu={nu}, expected {mu}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // effective_viscosity dispatch: Bingham below yield is RIGID_MU
    // -----------------------------------------------------------------------
    #[test]
    fn test_dispatch_bingham_below_yield() {
        let model = RheologyModel::Bingham(BinghamFluid::new(10.0, 0.01));
        let nu = effective_viscosity(model, 0.0);
        assert!(
            nu > 1.0,
            "Bingham below yield via dispatch should return large viscosity, got {nu}"
        );
    }

    // -----------------------------------------------------------------------
    // effective_viscosity dispatch: Carreau-Yasuda zero-shear limit
    // -----------------------------------------------------------------------
    #[test]
    fn test_dispatch_carreau_yasuda_zero_shear() {
        let cy = CarreauYasudaFluid::new(1.0, 0.01, 2.0, 2.0, 0.5);
        let model = RheologyModel::CarreauYasuda(cy);
        let nu = effective_viscosity(model, 0.0);
        assert!((nu - 1.0).abs() < 1e-12, "Dispatch CY at γ̇=0: {nu}");
    }

    // -----------------------------------------------------------------------
    // relaxation_time: Newtonian reference
    // -----------------------------------------------------------------------
    #[test]
    fn test_relaxation_time_newtonian() {
        // nu = 1/6, tau = 0.5 + (1/6) / (1/3) = 0.5 + 0.5 = 1.0
        let model = RheologyModel::PowerLaw(PowerLawFluid::new(1.0, 1.0 / 6.0));
        let tau = relaxation_time(model, 1.0);
        assert!(
            (tau - 1.0).abs() < 1e-12,
            "Relaxation time = {tau}, expected 1.0"
        );
    }

    // -----------------------------------------------------------------------
    // relaxation_frequency: omega = 1/tau
    // -----------------------------------------------------------------------
    #[test]
    fn test_relaxation_frequency_newtonian() {
        let model = RheologyModel::PowerLaw(PowerLawFluid::new(1.0, 1.0 / 6.0));
        let omega = relaxation_frequency(model, 1.0);
        assert!((omega - 1.0).abs() < 1e-12, "Omega = {omega}, expected 1.0");
    }

    // -----------------------------------------------------------------------
    // is_yielded_2d: pure shear exceeding yield stress
    // -----------------------------------------------------------------------
    #[test]
    fn test_is_yielded_2d_above_yield() {
        // s_xy = 5.0, tau_y = 1.0; sigma_vm = sqrt(0 + 0 + 2 * 25) = sqrt(50) ≈ 7.07 > 1
        let yielded = is_yielded_2d(0.0, 0.0, 5.0, 1.0);
        assert!(yielded, "Pure shear above yield stress should be yielded");
    }

    // -----------------------------------------------------------------------
    // is_yielded_2d: below yield stress → un-yielded
    // -----------------------------------------------------------------------
    #[test]
    fn test_is_yielded_2d_below_yield() {
        // sigma_vm = sqrt(2 * 0.01) ≈ 0.14 < 1.0
        let yielded = is_yielded_2d(0.0, 0.0, 0.1, 1.0);
        assert!(!yielded, "Small shear should be below yield stress");
    }

    // -----------------------------------------------------------------------
    // von_mises_stress_2d: known value
    // -----------------------------------------------------------------------
    #[test]
    fn test_von_mises_stress_2d() {
        // s_xx=3, s_yy=4, s_xy=0: vm = sqrt(9 + 16 + 0) = 5
        let vm = von_mises_stress_2d(3.0, 4.0, 0.0);
        assert!((vm - 5.0).abs() < 1e-12, "von Mises = {vm}, expected 5.0");
    }

    // -----------------------------------------------------------------------
    // bingham_yield_correction: above yield → stress unchanged
    // -----------------------------------------------------------------------
    #[test]
    fn test_bingham_yield_correction_above() {
        let corrected = bingham_yield_correction(5.0, 2.0);
        assert!(
            (corrected - 5.0).abs() < 1e-14,
            "Above yield: stress unchanged"
        );
    }

    // -----------------------------------------------------------------------
    // bingham_yield_correction: below yield → zero
    // -----------------------------------------------------------------------
    #[test]
    fn test_bingham_yield_correction_below() {
        let corrected = bingham_yield_correction(1.0, 2.0);
        assert!(
            (corrected - 0.0).abs() < 1e-14,
            "Below yield: stress zeroed"
        );
    }
}

mod new_nn_tests {
    use super::*;

    // -----------------------------------------------------------------------
    // PowerLaw: Newtonian case n=1
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_newtonian_case() {
        let pl = PowerLaw::new(0.01, 1.0);
        for &gamma in &[0.1, 1.0, 10.0, 100.0] {
            let mu = pl.effective_viscosity(gamma);
            assert!(
                (mu - 0.01).abs() < 1e-12,
                "PowerLaw n=1: mu={mu} at gamma={gamma}, expected 0.01"
            );
        }
    }

    // -----------------------------------------------------------------------
    // PowerLaw: shear-thinning (n < 1)
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_shear_thinning_behavior() {
        let pl = PowerLaw::new(1.0, 0.5);
        assert!(pl.is_shear_thinning(), "n=0.5 should be shear-thinning");
        let mu_low = pl.effective_viscosity(0.1);
        let mu_high = pl.effective_viscosity(10.0);
        assert!(
            mu_high < mu_low,
            "Shear-thinning: mu at high shear ({mu_high}) should be < mu at low shear ({mu_low})"
        );
    }

    // -----------------------------------------------------------------------
    // PowerLaw: trait dispatch via NonNewtonianModel
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_non_newtonian_model_trait() {
        let pl = PowerLaw::new(0.1, 0.8);
        let model: &dyn NonNewtonianModel = &pl;
        let mu = model.viscosity(5.0);
        let expected = pl.effective_viscosity(5.0);
        assert!(
            (mu - expected).abs() < 1e-14,
            "Trait dispatch should match direct call: {mu} vs {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // PowerLaw: local_tau is physically valid
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_local_tau_positive() {
        let pl = PowerLaw::new(1.0 / 6.0, 1.0);
        let tau = pl.local_tau(1.0);
        // nu = 1/6, tau = 0.5 + (1/6) / (1/3) = 0.5 + 0.5 = 1.0
        assert!((tau - 1.0).abs() < 1e-12, "tau should be 1.0, got {tau}");
    }

    // -----------------------------------------------------------------------
    // Carreau: zero-shear viscosity = mu_0
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_zero_shear_is_mu0() {
        let c = Carreau::new(1.5, 0.01, 2.0, 0.6);
        let mu = c.effective_viscosity(0.0);
        assert!(
            (mu - 1.5).abs() < 1e-12,
            "Carreau at gamma=0 should return mu_0=1.5, got {mu}"
        );
    }

    // -----------------------------------------------------------------------
    // Carreau: shear thinning
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_shear_thinning() {
        let c = Carreau::new(1.0, 0.001, 1.0, 0.5);
        let mu_low = c.effective_viscosity(0.01);
        let mu_high = c.effective_viscosity(100.0);
        assert!(
            mu_high < mu_low,
            "Carreau shear-thinning: mu_low={mu_low}, mu_high={mu_high}"
        );
    }

    // -----------------------------------------------------------------------
    // Carreau: NonNewtonianModel trait
    // -----------------------------------------------------------------------
    #[test]
    fn test_carreau_non_newtonian_model_trait() {
        let c = Carreau::new(1.0, 0.001, 2.0, 0.5);
        let model: &dyn NonNewtonianModel = &c;
        let mu = model.viscosity(1.0);
        let expected = c.effective_viscosity(1.0);
        assert!(
            (mu - expected).abs() < 1e-14,
            "Carreau trait dispatch: {mu} vs {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // effective_tau: standard lattice units (nu=1/6, cs2=1/3, dt=1) → tau=1
    // -----------------------------------------------------------------------
    #[test]
    fn test_effective_tau_standard_lattice_units() {
        let nu = 1.0_f64 / 6.0;
        let cs2 = 1.0_f64 / 3.0;
        let dt = 1.0_f64;
        let tau = effective_tau(nu, cs2, dt);
        assert!(
            (tau - 1.0).abs() < 1e-12,
            "effective_tau = {tau}, expected 1.0"
        );
    }

    // -----------------------------------------------------------------------
    // effective_tau: tau > 0.5 always (for positive viscosity)
    // -----------------------------------------------------------------------
    #[test]
    fn test_effective_tau_always_above_half() {
        let cs2 = 1.0_f64 / 3.0;
        for &nu in &[0.001, 0.01, 0.1, 1.0] {
            let tau = effective_tau(nu, cs2, 1.0);
            assert!(
                tau > 0.5,
                "effective_tau must be > 0.5, got {tau} for nu={nu}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // apply_non_newtonian_collision: conserves mass
    // -----------------------------------------------------------------------
    #[test]
    fn test_apply_non_newtonian_collision_conserves_mass() {
        let rho = 1.2_f64;
        let u = [0.05_f64, 0.0];
        // Use a plain Carreau model
        let model = Carreau::new(1.0, 0.01, 1.0, 0.5);
        // Initialize f to feq + small perturbation
        let ux = u[0];
        let uy = u[1];
        let u_sq = ux * ux + uy * uy;
        let mut f = [0.0_f64; 9];
        for i in 0..9 {
            let cx = C9_NN[i][0] as f64;
            let cy = C9_NN[i][1] as f64;
            let eu = cx * ux + cy * uy;
            f[i] = W9_NN[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        f[0] += 0.01;
        f[1] -= 0.005;
        f[2] -= 0.005;
        let mass_before: f64 = f.iter().sum();
        apply_non_newtonian_collision(&mut f, rho, u, &model);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "Non-Newtonian collision must conserve mass: before={mass_before}, after={mass_after}"
        );
    }

    // -----------------------------------------------------------------------
    // apply_non_newtonian_collision: with Bingham model
    // -----------------------------------------------------------------------
    #[test]
    fn test_apply_non_newtonian_collision_bingham() {
        let rho = 1.0_f64;
        let u = [0.01_f64, 0.0];
        let model = BinghamPlastic::new(0.1, 0.01);
        let ux = u[0];
        let uy = u[1];
        let u_sq = ux * ux + uy * uy;
        let mut f = [0.0_f64; 9];
        for i in 0..9 {
            let cx = C9_NN[i][0] as f64;
            let cy = C9_NN[i][1] as f64;
            let eu = cx * ux + cy * uy;
            f[i] = W9_NN[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        let mass_before: f64 = f.iter().sum();
        apply_non_newtonian_collision(&mut f, rho, u, &model);
        let mass_after: f64 = f.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "Bingham collision must conserve mass: {mass_before} vs {mass_after}"
        );
    }

    // -----------------------------------------------------------------------
    // BinghamPlastic implements NonNewtonianModel
    // -----------------------------------------------------------------------
    #[test]
    fn test_bingham_plastic_non_newtonian_model() {
        let b = BinghamPlastic::new(1.0, 0.1);
        let model: &dyn NonNewtonianModel = &b;
        let mu = model.viscosity(100.0);
        assert!(mu > 0.0 && mu.is_finite(), "BinghamPlastic viscosity={mu}");
    }

    // -----------------------------------------------------------------------
    // HerschelBulkley implements NonNewtonianModel
    // -----------------------------------------------------------------------
    #[test]
    fn test_herschel_bulkley_non_newtonian_model() {
        let hb = HerschelBulkley::new(1.0, 0.5, 1.0);
        let model: &dyn NonNewtonianModel = &hb;
        let mu = model.viscosity(5.0);
        assert!(mu > 0.0 && mu.is_finite(), "HerschelBulkley viscosity={mu}");
    }

    // -----------------------------------------------------------------------
    // PowerLaw: mass normalization of equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_power_law_equilibrium_via_collision_mass() {
        // Start at equilibrium, collision should not change f.
        let rho = 1.0_f64;
        let u = [0.0_f64, 0.0];
        let pl = PowerLaw::new(1.0 / 6.0, 1.0);
        let mut f = [0.0_f64; 9];
        for i in 0..9 {
            f[i] = W9_NN[i] * rho;
        }
        let f_before = f;
        apply_non_newtonian_collision(&mut f, rho, u, &pl);
        // Should be unchanged (feq = f at equilibrium, no matter what omega is)
        let diff: f64 = f
            .iter()
            .zip(f_before.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(
            diff < 1e-13,
            "Equilibrium unchanged by collision: diff={diff}"
        );
    }
}
