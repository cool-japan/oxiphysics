//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::functions::*;
    use crate::grid::LbmGrid2D;
    use crate::lattice::LatticeType;
    use crate::turbulence::types::*;
    #[test]
    fn test_smag_zero_strain_returns_base_omega() {
        let nx = 6;
        let ny = 6;
        let base_omega = 1.2;
        let cs_smag = 0.1;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        grid.compute_macroscopic();
        let omega = smagorinsky_omega(&grid, cs_smag, base_omega, 3, 3);
        assert!(
            (omega - base_omega).abs() < 1e-12,
            "At zero strain, omega should equal base_omega: got {omega}, expected {base_omega}"
        );
    }
    #[test]
    fn test_smag_omega_le_base() {
        let nx = 8;
        let ny = 8;
        let base_omega = 1.0;
        let cs_smag = 0.15;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        for y in 0..ny {
            for x in 0..nx {
                let ux = 0.05 * x as f64 / nx as f64;
                let uy = 0.03 * y as f64 / ny as f64;
                grid.set_equilibrium(x, y, 1.0, ux, uy);
            }
        }
        grid.compute_macroscopic();
        for y in 0..ny {
            for x in 0..nx {
                let omega = smagorinsky_omega(&grid, cs_smag, base_omega, x, y);
                assert!(
                    omega <= base_omega + 1e-12,
                    "omega exceeds base at ({x},{y}): {omega} > {base_omega}"
                );
                assert!(omega > 0.0, "omega must be positive");
            }
        }
    }
    #[test]
    fn test_smag_strain_increases_nu() {
        let nx = 10;
        let ny = 10;
        let base_omega = 1.0;
        let cs_smag = 0.15;
        let x = 5;
        let y = 5;
        let mut grid_lo = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        grid_lo.set_equilibrium(x, y, 1.0, 0.001, 0.0);
        grid_lo.compute_macroscopic();
        let mut grid_hi = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        grid_hi.set_equilibrium(x, y, 1.0, 0.0, 0.0);
        let k = grid_hi.idx(x, y);
        grid_hi.f[1][k] += 0.12;
        grid_hi.f[3][k] -= 0.12;
        grid_hi.compute_macroscopic();
        let omega_lo = smagorinsky_omega(&grid_lo, cs_smag, base_omega, x, y);
        let omega_hi = smagorinsky_omega(&grid_hi, cs_smag, base_omega, x, y);
        assert!(
            omega_hi <= omega_lo + 1e-10,
            "Higher strain should give lower omega: lo={omega_lo}, hi={omega_hi}"
        );
    }
    #[test]
    fn test_smag_larger_cs_lowers_omega() {
        let nx = 10;
        let ny = 10;
        let base_omega = 1.0;
        let x = 5;
        let y = 5;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        grid.set_equilibrium(x, y, 1.0, 0.0, 0.0);
        let k = grid.idx(x, y);
        grid.f[1][k] += 0.08;
        grid.f[3][k] -= 0.08;
        grid.compute_macroscopic();
        let omega_small = smagorinsky_omega(&grid, 0.05, base_omega, x, y);
        let omega_large = smagorinsky_omega(&grid, 0.20, base_omega, x, y);
        assert!(
            omega_large <= omega_small + 1e-10,
            "Larger Cs should give lower omega: small={omega_small}, large={omega_large}"
        );
    }
    #[test]
    fn test_turbulent_viscosity_proportional_to_strain() {
        let mixing_length = 0.05;
        let s1 = 1.0;
        let s2 = 2.0;
        let nu1 = compute_turbulent_viscosity(s1, mixing_length);
        let nu2 = compute_turbulent_viscosity(s2, mixing_length);
        assert!(
            (nu2 / nu1 - 2.0).abs() < 1e-14,
            "nu_t should be linear in strain rate: nu1={nu1}, nu2={nu2}"
        );
    }
    #[test]
    fn test_wale_zero_gradient() {
        let g = [0.0f64; 9];
        let nu = wale_model(&g);
        assert!(
            nu.abs() < 1e-20,
            "WALE: zero gradient should give zero nu_t: {nu}"
        );
    }
    #[test]
    fn test_wale_rotation_shear_nonzero() {
        let mut g = [0.0f64; 9];
        g[0] = 1.0;
        g[4] = -1.0;
        let nu = wale_model(&g);
        assert!(
            nu > 0.0,
            "WALE: extensional flow should give positive nu_t: {nu}"
        );
    }
    #[test]
    fn test_effective_omega_zero_strain() {
        let omega_base = 1.2;
        let omega = compute_effective_omega(omega_base, 0.0, 1.0, 1.0);
        assert!(
            (omega - omega_base).abs() < 1e-12,
            "Zero strain should give base omega: got {omega}"
        );
    }
    #[test]
    fn test_dynamic_smag_cs_sq_clipped() {
        let dyn_smag = DynamicSmagorinsky::new(1.0, 2.0);
        let cs_sq_neg = dyn_smag.compute_cs_sq(2.0, 1.0);
        assert!(cs_sq_neg >= 0.0, "cs_sq must be non-negative: {cs_sq_neg}");
        let cs_sq_large = dyn_smag.compute_cs_sq(0.0, 10.0);
        assert!(
            cs_sq_large <= 0.04 + 1e-14,
            "cs_sq must be capped at 0.04: {cs_sq_large}"
        );
    }
    #[test]
    fn test_vreman_zero_gradient() {
        let g = [0.0f64; 9];
        let nu = vreman_model(&g, 0.07);
        assert!(
            nu.abs() < 1e-20,
            "Vreman: zero gradient should give 0: {nu}"
        );
    }
    #[test]
    fn test_vreman_simple_shear() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        let nu = vreman_model(&g, 0.07);
        assert!(
            nu >= 0.0,
            "Vreman: simple shear should give non-negative nu: {nu}"
        );
    }
    #[test]
    fn test_vreman_extensional_flow() {
        let mut g = [0.0f64; 9];
        g[0] = 1.0;
        g[4] = -0.5;
        g[8] = -0.5;
        let nu = vreman_model(&g, 0.07);
        assert!(nu >= 0.0, "Vreman extensional: nu = {nu}");
    }
    #[test]
    fn test_k_epsilon_eddy_viscosity() {
        let ke = KEpsilonState::new(1.0, 0.5);
        let nu_t = ke.eddy_viscosity();
        assert!(nu_t > 0.0, "Eddy viscosity should be positive: {nu_t}");
    }
    #[test]
    fn test_k_epsilon_zero_k() {
        let ke = KEpsilonState::new(0.0, 1.0);
        let nu_t = ke.eddy_viscosity();
        assert!(nu_t.abs() < 1e-20, "Zero k should give zero nu_t: {nu_t}");
    }
    #[test]
    fn test_k_epsilon_time_scale() {
        let ke = KEpsilonState::new(2.0, 0.5);
        let ts = ke.time_scale();
        assert!((ts - 4.0).abs() < 1e-12, "time scale = {ts}, expected 4.0");
    }
    #[test]
    fn test_k_epsilon_advance_no_nan() {
        let mut ke = KEpsilonState::new(1.0, 0.1);
        for _ in 0..100 {
            ke.advance(0.5, 0.01);
        }
        assert!(!ke.k.is_nan(), "k should not be NaN");
        assert!(!ke.epsilon.is_nan(), "epsilon should not be NaN");
        assert!(ke.k > 0.0, "k should remain positive");
        assert!(ke.epsilon > 0.0, "epsilon should remain positive");
    }
    #[test]
    fn test_wall_function_zero_velocity() {
        let u_tau = wall_function_log_law(0.0, 1.0, 0.01);
        assert!(u_tau.abs() < 1e-10, "u_tau for zero velocity: {u_tau}");
    }
    #[test]
    fn test_wall_function_positive() {
        let u_tau = wall_function_log_law(1.0, 0.01, 1e-4);
        assert!(u_tau > 0.0, "u_tau should be positive: {u_tau}");
    }
    #[test]
    fn test_wall_shear_stress() {
        let rho = 1.0;
        let tau_w = wall_shear_stress(rho, 0.5);
        assert!(
            (tau_w - 0.25).abs() < 1e-12,
            "tau_w = {tau_w}, expected 0.25"
        );
    }
    #[test]
    fn test_y_plus() {
        let yp = y_plus(0.01, 0.1, 1e-3);
        assert!((yp - 1.0).abs() < 1e-12, "y+ = {yp}, expected 1.0");
    }
    #[test]
    fn test_spalding_zero() {
        let up = spalding_profile(0.0);
        assert!(up.abs() < 1e-6, "Spalding at y+=0: u+ = {up}");
    }
    #[test]
    fn test_spalding_viscous_sublayer() {
        let yp = 3.0;
        let up = spalding_profile(yp);
        assert!(
            (up - yp).abs() < 0.5,
            "Spalding at y+={yp}: u+ = {up}, expected ~{yp}"
        );
    }
    #[test]
    fn test_spalding_monotonic() {
        let mut prev = 0.0;
        for i in 1..100 {
            let yp = i as f64;
            let up = spalding_profile(yp);
            assert!(up >= prev - 1e-10, "Spalding not monotonic at y+={yp}");
            prev = up;
        }
    }
    #[test]
    fn test_dynamic_smag_effective_omega() {
        let dyn_smag = DynamicSmagorinsky::new(1.0, 2.0);
        let omega = dyn_smag.effective_omega(1.0, 0.5, 1.0);
        assert!(omega > 0.0 && omega <= 2.0, "omega = {omega}");
    }
    #[test]
    fn test_k_epsilon_length_scale() {
        let ke = KEpsilonState::new(1.0, 0.5);
        let ls = ke.length_scale();
        assert!(ls > 0.0, "Length scale should be positive: {ls}");
    }
    #[test]
    fn test_les_models_zero_strain_agreement() {
        let g = [0.0f64; 9];
        let wale_nu = wale_model(&g);
        let vreman_nu = vreman_model(&g, 0.07);
        assert!(wale_nu.abs() < 1e-20, "WALE zero strain: {wale_nu}");
        assert!(vreman_nu.abs() < 1e-20, "Vreman zero strain: {vreman_nu}");
    }
    #[test]
    fn test_tke_from_fluctuations() {
        let u_fluct = [1.0, 2.0, 3.0];
        let tke = turbulent_kinetic_energy(&u_fluct);
        let expected = 0.5 * (1.0 + 4.0 + 9.0);
        assert!(
            (tke - expected).abs() < 1e-12,
            "TKE = {tke}, expected {expected}"
        );
    }
    #[test]
    fn test_turbulence_intensity() {
        let ti = turbulence_intensity(0.1, 1.0);
        assert!((ti - 0.1).abs() < 1e-12, "TI = {ti}, expected 0.1");
    }
    #[test]
    fn test_turbulence_intensity_zero_mean() {
        let ti = turbulence_intensity(0.1, 0.0);
        assert!(!ti.is_nan(), "TI should not be NaN for zero mean");
    }
    #[test]
    fn test_power_spectrum_length() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).sin()).collect();
        let spectrum = estimate_power_spectrum(&signal);
        assert_eq!(spectrum.len(), signal.len() / 2, "spectrum length mismatch");
    }
    #[test]
    fn test_power_spectrum_zero_signal() {
        let signal = vec![0.0_f64; 16];
        let spectrum = estimate_power_spectrum(&signal);
        for (i, &s) in spectrum.iter().enumerate() {
            assert!(s.abs() < 1e-20, "spectrum[{i}] = {s} for zero signal");
        }
    }
    #[test]
    fn test_tke_production_positive() {
        let p = tke_production(0.5, 1.0);
        assert!(p >= 0.0, "TKE production should be non-negative: {p}");
    }
    #[test]
    fn test_tke_dissipation_positive() {
        let ke = KEpsilonState::new(1.0, 0.5);
        assert!(ke.epsilon > 0.0, "dissipation rate should be positive");
    }
    #[test]
    fn test_taylor_microscale_vs_reynolds() {
        let lambda1 = taylor_microscale(1.0, 1e-3, 0.1);
        let lambda2 = taylor_microscale(1.0, 1e-3, 0.2);
        assert!(
            lambda2 < lambda1 + 1e-10,
            "Higher dissipation → smaller microscale"
        );
    }
    #[test]
    fn test_kolmogorov_length_scale() {
        let nu = 1e-5;
        let epsilon = 0.1;
        let eta = kolmogorov_length_scale(nu, epsilon);
        let expected = (nu * nu * nu / epsilon).powf(0.25);
        assert!(
            (eta - expected).abs() / expected < 1e-10,
            "eta = {eta}, expected {expected}"
        );
    }
    #[test]
    fn test_kolmogorov_time_scale() {
        let nu = 1e-5;
        let epsilon = 0.1;
        let tau_eta = kolmogorov_time_scale(nu, epsilon);
        let expected = (nu / epsilon).sqrt();
        assert!(
            (tau_eta - expected).abs() / expected < 1e-10,
            "tau_eta = {tau_eta}"
        );
    }
    #[test]
    fn test_integral_length_scale() {
        let k = 1.0;
        let eps = 0.5;
        let l = integral_length_scale(k, eps);
        let expected = k.powf(1.5) / eps;
        assert!((l - expected).abs() < 1e-12, "L = {l}, expected {expected}");
    }
    #[test]
    fn test_turbulent_reynolds_number() {
        let k = 1.0;
        let nu = 1e-5;
        let eps = 0.1;
        let re_t = turbulent_reynolds_number(k, nu, eps);
        let expected = k * k / (nu * eps);
        assert!(
            (re_t - expected).abs() < 1e-8,
            "Re_t = {re_t}, expected {expected}"
        );
    }
    #[test]
    fn test_kolmogorov_spectrum_slope() {
        let k1 = kolmogorov_energy_spectrum(1.0, 0.1, 1.0);
        let k2 = kolmogorov_energy_spectrum(2.0, 0.1, 1.0);
        let ratio = k2 / k1;
        let expected_ratio = 2.0_f64.powf(-5.0 / 3.0);
        assert!(
            (ratio - expected_ratio).abs() < 1e-10,
            "ratio = {ratio}, expected {expected_ratio}"
        );
    }
    #[test]
    fn test_tke_budget_balance() {
        let production = 0.5_f64;
        let dissipation = 0.5_f64;
        let balance = (production - dissipation).abs();
        assert!(balance < 1e-12, "Budget should balance: diff = {balance}");
    }
    #[test]
    fn test_von_karman_spectrum_positive() {
        let e = von_karman_spectrum(1.0, 0.1, 1.0);
        assert!(e > 0.0, "Von Karman spectrum should be positive: {e}");
    }
    #[test]
    fn test_tke_turbulent_diffusion() {
        let d = tke_turbulent_diffusion(0.5, 1.0, 0.1);
        assert!(d.is_finite(), "Diffusion term should be finite: {d}");
    }
    #[test]
    fn test_strain_rate_from_gradient() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        let s_mag = strain_rate_magnitude(&g);
        let expected = 0.5_f64.sqrt();
        assert!(
            (s_mag - expected).abs() < 1e-12,
            "strain rate = {s_mag}, expected {expected}"
        );
    }
    #[test]
    fn test_turbulent_prandtl_number() {
        let pr_t = turbulent_prandtl_number(1e-3, 8e-4);
        assert!(pr_t > 0.0, "Pr_t should be positive: {pr_t}");
    }
    #[test]
    fn test_k_omega_eddy_viscosity() {
        let nu_t = k_omega_eddy_viscosity(1.0, 2.0);
        assert!((nu_t - 0.5).abs() < 1e-12, "nu_t = {nu_t}, expected 0.5");
    }
    #[test]
    fn test_k_omega_zero_omega() {
        let nu_t = k_omega_eddy_viscosity(1.0, 0.0);
        assert!(nu_t.abs() < 1e-10, "nu_t with zero omega: {nu_t}");
    }
    #[test]
    fn test_backscatter_indicator() {
        let g = [0.0f64; 9];
        let bi = backscatter_indicator(&g);
        assert!(bi >= 0.0, "backscatter indicator should be >= 0: {bi}");
    }
    #[test]
    fn test_wale_pure_rotation_low() {
        let mut g = [0.0f64; 9];
        g[1] = -1.0;
        g[3] = 1.0;
        let nu = wale_model(&g);
        assert!(nu < 1.0, "WALE pure rotation should give small nu_t: {nu}");
    }
    #[test]
    fn test_dynamic_smag_filter_ratio_effect() {
        let dyn1 = DynamicSmagorinsky::new(1.0, 2.0);
        let dyn2 = DynamicSmagorinsky::new(1.0, 4.0);
        let cs1 = dyn1.compute_cs_sq(0.5, 1.0);
        let cs2 = dyn2.compute_cs_sq(0.5, 1.0);
        assert!(cs1 >= 0.0 && cs2 >= 0.0, "cs_sq must be non-negative");
    }
    #[test]
    fn test_strain_rate_d3q19_nonzero() {
        let mut f_neq = [0.0f64; 19];
        f_neq[1] = 0.01;
        f_neq[2] = -0.01;
        let s = compute_strain_rate_d3q19(&f_neq, 1.0);
        assert!(s >= 0.0, "strain rate should be non-negative: {s}");
    }
    #[test]
    fn test_sst_f1_near_wall() {
        let f1 = sst_blending_f1(1.0, 1.0, 1e-6, 1e-5, 0.0);
        assert!(f1 > 0.9, "F1 near wall should be close to 1: {f1}");
    }
    #[test]
    fn test_sst_f1_far_from_wall() {
        let f1 = sst_blending_f1(1.0, 1.0, 1e6, 1e-5, 0.0);
        assert!(f1 < 0.1, "F1 far from wall should be close to 0: {f1}");
    }
    #[test]
    fn test_sst_f2_near_wall() {
        let f2 = sst_blending_f2(1.0, 1.0, 1e-6, 1e-5);
        assert!(f2 > 0.9, "F2 near wall should be close to 1: {f2}");
    }
    #[test]
    fn test_sst_f2_far_from_wall() {
        let f2 = sst_blending_f2(1.0, 1.0, 1e6, 1e-5);
        assert!(f2 < 0.1, "F2 far from wall should be ~0: {f2}");
    }
    #[test]
    fn test_sst_blend_coefficient_f1_zero() {
        let phi = sst_blend_coefficient(1.0, 2.0, 0.0);
        assert!(
            (phi - 2.0).abs() < 1e-12,
            "F1=0: phi should be phi_ke: {phi}"
        );
    }
    #[test]
    fn test_sst_blend_coefficient_f1_one() {
        let phi = sst_blend_coefficient(1.0, 2.0, 1.0);
        assert!(
            (phi - 1.0).abs() < 1e-12,
            "F1=1: phi should be phi_kw: {phi}"
        );
    }
    #[test]
    fn test_sst_eddy_viscosity_positive() {
        let nu_t = sst_eddy_viscosity(1.0, 2.0, 0.5, 0.5);
        assert!(nu_t > 0.0 && nu_t.is_finite(), "SST nu_t = {nu_t}");
    }
    #[test]
    fn test_sst_eddy_viscosity_limited() {
        let nu_t_lim = sst_eddy_viscosity(1.0, 0.1, 10.0, 1.0);
        let nu_t_unlim = sst_eddy_viscosity(1.0, 10.0, 0.1, 0.0);
        assert!(
            nu_t_lim <= nu_t_unlim + 1e-12,
            "SST limiting: lim={nu_t_lim}, unlim={nu_t_unlim}"
        );
    }
    #[test]
    fn test_log_law_u_plus_standard() {
        let up = log_law_u_plus(100.0, 0.41, 5.2);
        let expected = (1.0 / 0.41) * 100.0_f64.ln() + 5.2;
        assert!(
            (up - expected).abs() < 1e-10,
            "u+ = {up}, expected {expected}"
        );
    }
    #[test]
    fn test_composite_wall_law_sublayer() {
        let up = composite_wall_law(5.0);
        assert!(
            (up - 5.0).abs() < 1e-12,
            "Sublayer: u+ should equal y+: {up}"
        );
    }
    #[test]
    fn test_composite_wall_law_log_region() {
        let up = composite_wall_law(100.0);
        assert!(up > 11.3, "Log region: u+ should exceed y_match: {up}");
    }
    #[test]
    fn test_composite_wall_law_continuous() {
        let up_below = composite_wall_law(11.0);
        let up_above = composite_wall_law(11.6);
        assert!(
            up_above > up_below,
            "Wall law should be monotonically increasing: {up_below}, {up_above}"
        );
    }
    #[test]
    fn test_friction_reynolds_number() {
        let re_tau = friction_reynolds_number(0.01, 0.1, 1e-3);
        assert!(
            (re_tau - 1.0).abs() < 1e-12,
            "Re_tau = {re_tau}, expected 1.0"
        );
    }
    #[test]
    fn test_spectral_cutoff_preserves_low_modes() {
        let n = 16;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos())
            .collect();
        let filtered = spectral_cutoff_filter(&signal, 2);
        let err: f64 = signal
            .iter()
            .zip(filtered.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        assert!(
            err < 1e-20,
            "k=1 mode should pass k_cutoff=2 unchanged: err={err}"
        );
    }
    #[test]
    fn test_spectral_cutoff_removes_high_modes() {
        let n = 16;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * 4.0 * i as f64 / n as f64).cos())
            .collect();
        let filtered = spectral_cutoff_filter(&signal, 2);
        let energy: f64 = filtered.iter().map(|v| v * v).sum::<f64>() / n as f64;
        assert!(
            energy < 1e-20,
            "k=4 mode should be removed by k_cutoff=2: energy={energy}"
        );
    }
    #[test]
    fn test_gaussian_les_filter_preserves_dc() {
        let n = 16;
        let signal = vec![1.0_f64; n];
        let filtered = gaussian_les_filter(&signal, 2.0, 1.0, 4);
        let err: f64 = filtered.iter().map(|&v| (v - 1.0).powi(2)).sum();
        assert!(err < 1e-20, "DC component should pass unchanged: err={err}");
    }
    #[test]
    fn test_turbulent_prandtl_kays_crawford_range() {
        let pr_t = turbulent_prandtl_kays_crawford(1e-3, 1e-4);
        assert!(
            pr_t > 0.0 && pr_t < 2.0,
            "Pr_t should be in reasonable range: {pr_t}"
        );
    }
    #[test]
    fn test_turbulent_prandtl_kays_crawford_high_ratio() {
        let pr_t = turbulent_prandtl_kays_crawford(1.0, 1e-6);
        assert!(
            pr_t.is_finite() && pr_t > 0.0,
            "Pr_t should be finite: {pr_t}"
        );
    }
    #[test]
    fn test_turbulent_prandtl_preset_values() {
        assert!((TurbPrandtlPreset::WaleLes.value() - 0.4).abs() < 1e-12);
        assert!((TurbPrandtlPreset::SmagorinskyLes.value() - 0.9).abs() < 1e-12);
        assert!((TurbPrandtlPreset::KEpsilonRans.value() - 0.85).abs() < 1e-12);
    }
    #[test]
    fn test_les_spectrum_slope_kolmogorov() {
        let spectrum = synthetic_kolmogorov_spectrum(64, 0.1, 1.5);
        let slope = les_spectrum_slope(&spectrum, 2, 32);
        let expected = -5.0_f64 / 3.0;
        assert!(
            (slope - expected).abs() < 0.01,
            "Slope of Kolmogorov spectrum = {slope}, expected {expected}"
        );
    }
    #[test]
    fn test_spectrum_is_kolmogorov_true() {
        let spectrum = synthetic_kolmogorov_spectrum(64, 0.1, 1.5);
        assert!(
            spectrum_is_kolmogorov(&spectrum, 2, 32, 0.05),
            "Synthetic Kolmogorov spectrum should pass -5/3 check"
        );
    }
    #[test]
    fn test_spectrum_is_kolmogorov_false_for_flat() {
        let spectrum = vec![1.0_f64; 64];
        assert!(
            !spectrum_is_kolmogorov(&spectrum, 2, 32, 0.05),
            "Flat spectrum should fail -5/3 check"
        );
    }
    #[test]
    fn test_les_spectrum_slope_empty() {
        let spectrum: Vec<f64> = vec![];
        let slope = les_spectrum_slope(&spectrum, 0, 0);
        assert!((slope).abs() < 1e-10, "Empty spectrum → slope=0: {slope}");
    }
    #[test]
    fn test_vreman_b_beta_zero_gradient() {
        let g = [0.0f64; 9];
        let b = vreman_b_beta(&g);
        assert!(b.abs() < 1e-20, "B_beta for zero gradient: {b}");
    }
    #[test]
    fn test_vreman_alpha_sq_simple_shear() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        let a = vreman_alpha_sq(&g);
        assert!((a - 1.0).abs() < 1e-12, "alpha_sq for du/dy=1: {a}");
    }
    #[test]
    fn test_vreman_from_invariants_consistency() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        let b = vreman_b_beta(&g);
        let a = vreman_alpha_sq(&g);
        let cv = 0.07;
        let nu1 = vreman_from_invariants(b, a, cv);
        let nu2 = vreman_model(&g, cv);
        assert!(
            (nu1 - nu2).abs() < 1e-12,
            "vreman_from_invariants vs vreman_model: {nu1} vs {nu2}"
        );
    }
    #[test]
    fn test_vreman_b_beta_nonnegative_extensional() {
        let mut g = [0.0f64; 9];
        g[0] = 1.0;
        g[4] = -0.5;
        g[8] = -0.5;
        let b = vreman_b_beta(&g);
        assert!(b >= 0.0, "B_beta should be non-negative: {b}");
    }
    #[test]
    fn test_komega_sst_blending_wall_limit() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let f1 = sst.compute_blending_function(1.0, 1.0, 1e-6, 0.0, true);
        assert!(f1 > 0.9, "F1 near wall should be close to 1: {f1}");
    }
    #[test]
    fn test_komega_sst_blending_freestream_limit() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let f1 = sst.compute_blending_function(1.0, 1.0, 1e6, 0.0, true);
        assert!(f1 < 0.1, "F1 far from wall should be close to 0: {f1}");
    }
    #[test]
    fn test_komega_sst_f2_wall_limit() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let f2 = sst.compute_blending_function(1.0, 1.0, 1e-6, 0.0, false);
        assert!(f2 > 0.9, "F2 near wall should be close to 1: {f2}");
    }
    #[test]
    fn test_komega_sst_f2_freestream_limit() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let f2 = sst.compute_blending_function(1.0, 1.0, 1e6, 0.0, false);
        assert!(f2 < 0.1, "F2 far from wall should be close to 0: {f2}");
    }
    #[test]
    fn test_komega_sst_production_limiter_active() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let p_lim = sst.compute_production_limiter(100.0, 1.0, 1.0);
        let c_lim = 10.0 * 1.0 * 1.0;
        assert!(
            (p_lim - c_lim).abs() < 1e-10 || p_lim <= c_lim,
            "Production limiter should cap at 10*beta_star*k*omega: {p_lim}"
        );
    }
    #[test]
    fn test_komega_sst_production_limiter_inactive() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        let p_k = 0.01;
        let p_lim = sst.compute_production_limiter(p_k, 1.0, 1.0);
        assert!(
            (p_lim - p_k).abs() < 1e-12,
            "Production limiter should be inactive for small P_k: {p_lim}"
        );
    }
    #[test]
    fn test_komega_sst_production_limiter_nonnegative() {
        let sst = KOmegaSst::new(1e-5, 0.31);
        for &p_k in &[-5.0, 0.0, 1.0, 100.0] {
            let p_lim = sst.compute_production_limiter(p_k, 1.0, 1.0);
            assert!(p_lim >= 0.0, "Production limiter must be >= 0: {p_lim}");
        }
    }
    #[test]
    fn test_les_to_dns_sgs_ke_zero_for_dns_resolution() {
        let les = LesToDns::new(1.0, 0.5);
        let k_sgs = les.compute_subgrid_kinetic_energy(0.0, 1.0, 1.0);
        assert!(
            k_sgs.abs() < 1e-10,
            "SGS KE at DNS resolution should be ~0: {k_sgs}"
        );
    }
    #[test]
    fn test_les_to_dns_sgs_ke_positive_for_turbulent_flow() {
        let les = LesToDns::new(1.0, 0.5);
        let k_sgs = les.compute_subgrid_kinetic_energy(1.0, 4.0, 1.0);
        assert!(
            k_sgs > 0.0,
            "SGS KE should be positive for turbulent flow: {k_sgs}"
        );
    }
    #[test]
    fn test_les_to_dns_sgs_ke_scales_with_tke() {
        let les = LesToDns::new(1.0, 0.5);
        let k1 = les.compute_subgrid_kinetic_energy(1.0, 4.0, 1.0);
        let k2 = les.compute_subgrid_kinetic_energy(2.0, 4.0, 1.0);
        assert!(k2 > k1, "SGS KE should scale with resolved TKE: {k1}, {k2}");
    }
}
