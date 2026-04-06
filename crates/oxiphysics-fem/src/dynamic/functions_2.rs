//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_dynamic_extended {
    use super::*;

    #[test]
    fn test_bathe_integrator_new_zeroed() {
        let b = BatheIntegrator::new(4);
        assert_eq!(b.n_dof, 4);
        assert!(b.u.iter().all(|&x| x == 0.0));
        assert!(b.v.iter().all(|&x| x == 0.0));
        assert!(b.a.iter().all(|&x| x == 0.0));
    }
    #[test]
    fn test_bathe_integrator_step_finite() {
        let n = 2;
        let mut b = BatheIntegrator::new(n);
        let m = vec![1.0_f64, 1.0];
        let k = vec![100.0_f64, 100.0];
        let c = vec![0.1_f64, 0.1];
        let f = vec![1.0_f64, 0.0];
        b.step_diagonal(&m, &k, &c, &f, 0.01);
        for &u in &b.u {
            assert!(u.is_finite(), "u must be finite: {u}");
        }
    }
    #[test]
    fn test_bathe_integrator_step_no_force_decays() {
        let n = 1;
        let mut b = BatheIntegrator::new(n);
        b.u[0] = 1.0;
        b.v[0] = 0.0;
        let m = vec![1.0_f64];
        let k = vec![400.0_f64];
        let c = vec![4.0_f64];
        let f = vec![0.0_f64];
        for _ in 0..200 {
            b.step_diagonal(&m, &k, &c, &f, 0.002);
        }
        assert!(b.u[0].abs() < 1.0, "displacement should decay: {}", b.u[0]);
    }
    #[test]
    fn test_bathe_integrator_velocity_finite() {
        let n = 3;
        let mut b = BatheIntegrator::new(n);
        b.v[1] = 0.5;
        let m = vec![2.0_f64; 3];
        let k = vec![500.0_f64; 3];
        let c = vec![1.0_f64; 3];
        let f = vec![0.0_f64; 3];
        b.step_diagonal(&m, &k, &c, &f, 0.005);
        for &v in &b.v {
            assert!(v.is_finite(), "v must be finite: {v}");
        }
    }
    #[test]
    fn test_modal_superposition_zero_force_finite() {
        let omega_n = vec![10.0_f64, 30.0];
        let zeta = vec![0.05_f64, 0.05];
        let phi = vec![vec![1.0_f64, 0.0], vec![0.0, 1.0]];
        let f_time: Vec<Vec<f64>> = (0..50).map(|_| vec![0.0, 0.0]).collect();
        let u0 = vec![0.01_f64, 0.0];
        let v0 = vec![0.0_f64, 0.0];
        let res = modal_superposition(&omega_n, &zeta, &phi, &f_time, 0.005, &u0, &v0);
        assert_eq!(res.disp.len(), 50);
        for row in &res.disp {
            for &d in row {
                assert!(d.is_finite(), "displacement must be finite: {d}");
            }
        }
    }
    #[test]
    fn test_modal_superposition_initial_decay() {
        let omega_n = vec![20.0_f64];
        let zeta = vec![0.1_f64];
        let phi = vec![vec![1.0_f64]];
        let f_time: Vec<Vec<f64>> = (0..200).map(|_| vec![0.0]).collect();
        let u0 = vec![1.0_f64];
        let v0 = vec![0.0_f64];
        let res = modal_superposition(&omega_n, &zeta, &phi, &f_time, 0.002, &u0, &v0);
        let u_final = res.disp.last().unwrap()[0].abs();
        assert!(
            u_final < 1.0,
            "should decay from initial condition: {u_final}"
        );
    }
    #[test]
    fn test_modal_superposition_modal_coords_length() {
        let omega_n = vec![5.0_f64, 15.0, 25.0];
        let zeta = vec![0.02_f64; 3];
        let phi = vec![
            vec![1.0_f64, 0.5, 0.0],
            vec![0.5, -1.0, 0.5],
            vec![0.0, 0.5, -1.0],
        ];
        let f_time: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 0.0, 0.0]).collect();
        let u0 = vec![0.0_f64; 3];
        let v0 = vec![0.0_f64; 3];
        let res = modal_superposition(&omega_n, &zeta, &phi, &f_time, 0.01, &u0, &v0);
        assert_eq!(res.modal_coords[0].len(), 3);
    }
    #[test]
    fn test_modal_superposition_empty_steps() {
        let omega_n = vec![10.0_f64];
        let zeta = vec![0.05_f64];
        let phi = vec![vec![1.0_f64]];
        let f_time: Vec<Vec<f64>> = vec![];
        let u0 = vec![0.0_f64];
        let v0 = vec![0.0_f64];
        let res = modal_superposition(&omega_n, &zeta, &phi, &f_time, 0.01, &u0, &v0);
        assert_eq!(res.disp.len(), 0);
    }
    #[test]
    fn test_response_spectrum_analysis_positive() {
        let input = ResponseSpectrumInput {
            omega_n: vec![10.0_f64, 30.0],
            zeta: vec![0.05_f64, 0.05],
            phi: vec![vec![1.0_f64, 0.0], vec![0.0, 1.0]],
            gamma: vec![1.0_f64, 0.5],
            sd: vec![0.01_f64, 0.005],
            sa: vec![1.0_f64, 4.5],
        };
        let (peak_u, peak_a) = response_spectrum_analysis(&input);
        for &u in &peak_u {
            assert!(
                u >= 0.0 && u.is_finite(),
                "peak_u must be non-negative finite: {u}"
            );
        }
        for &a in &peak_a {
            assert!(
                a >= 0.0 && a.is_finite(),
                "peak_a must be non-negative finite: {a}"
            );
        }
    }
    #[test]
    fn test_response_spectrum_analysis_zero_sd() {
        let input = ResponseSpectrumInput {
            omega_n: vec![10.0_f64],
            zeta: vec![0.05_f64],
            phi: vec![vec![1.0_f64, 0.0]],
            gamma: vec![1.0_f64],
            sd: vec![0.0_f64],
            sa: vec![0.0_f64],
        };
        let (peak_u, _peak_a) = response_spectrum_analysis(&input);
        assert_eq!(peak_u[0], 0.0);
    }
    #[test]
    fn test_compute_sd_positive() {
        let n = 200;
        let dt = 0.005_f64;
        let omega_exc = 5.0_f64;
        let ag: Vec<f64> = (0..n).map(|i| (omega_exc * i as f64 * dt).sin()).collect();
        let sd = compute_sd(&ag, dt, 10.0, 0.05);
        assert!(
            sd >= 0.0 && sd.is_finite(),
            "Sd must be non-negative finite: {sd}"
        );
    }
    #[test]
    fn test_compute_sa_ge_sd_times_omega2() {
        let n = 100;
        let dt = 0.01_f64;
        let ag: Vec<f64> = (0..n).map(|i| (3.0 * i as f64 * dt).sin()).collect();
        let omega = 8.0_f64;
        let zeta = 0.05_f64;
        let sd = compute_sd(&ag, dt, omega, zeta);
        let sa = compute_sa(&ag, dt, omega, zeta);
        assert!(
            (sa - omega * omega * sd).abs() < 1e-10,
            "Sa = ω² Sd should hold exactly"
        );
    }
    #[test]
    fn test_mass_proportional_alpha_positive() {
        let alpha = mass_proportional_damping_alpha(10.0, 0.05);
        assert!(alpha > 0.0 && alpha.is_finite(), "alpha = {alpha}");
    }
    #[test]
    fn test_stiffness_proportional_beta_positive() {
        let beta = stiffness_proportional_damping_beta(10.0, 0.05);
        assert!(beta > 0.0 && beta.is_finite(), "beta = {beta}");
    }
    #[test]
    fn test_rayleigh_modal_damping_recovery() {
        let omega = 20.0_f64;
        let alpha = mass_proportional_damping_alpha(omega, 0.05);
        let beta_val = stiffness_proportional_damping_beta(omega, 0.05);
        let zeta = rayleigh_modal_damping(omega, alpha, beta_val);
        assert!((zeta - 0.1).abs() < 1e-10, "combined zeta = {zeta}");
    }
    #[test]
    fn test_mass_proportional_force_scales() {
        let m = vec![2.0_f64, 3.0];
        let v = vec![1.0_f64, 1.0];
        let alpha = 0.5_f64;
        let d = mass_proportional_damping_force(&m, &v, alpha);
        assert!((d[0] - 1.0).abs() < 1e-10, "d[0] = {}", d[0]);
        assert!((d[1] - 1.5).abs() < 1e-10, "d[1] = {}", d[1]);
    }
    #[test]
    fn test_stiffness_proportional_force_scales() {
        let k = vec![100.0_f64, 200.0];
        let v = vec![0.5_f64, 0.5];
        let beta = 0.01_f64;
        let d = stiffness_proportional_damping_force(&k, &v, beta);
        assert!((d[0] - 0.5).abs() < 1e-10, "d[0] = {}", d[0]);
        assert!((d[1] - 1.0).abs() < 1e-10, "d[1] = {}", d[1]);
    }
    #[test]
    fn test_cqc_rho_self_correlation_is_one() {
        let rho = cqc_rho(10.0, 10.0, 0.05, 0.05);
        assert!((rho - 1.0).abs() < 1e-8, "rho_ii should be 1: {rho}");
    }
    #[test]
    fn test_cqc_rho_bounded() {
        let rho = cqc_rho(10.0, 30.0, 0.05, 0.05);
        assert!(rho >= 0.0 && rho <= 1.0, "rho must be in [0,1]: {rho}");
    }
    #[test]
    fn test_cqc_combination_individual_zeta_single_mode() {
        let peaks = vec![5.0_f64];
        let omega = vec![10.0_f64];
        let zeta = vec![0.05_f64];
        let result = cqc_combination_individual_zeta(&peaks, &omega, &zeta);
        assert!(
            (result - 5.0).abs() < 1e-8,
            "single mode CQC = SRSS = peak: {result}"
        );
    }
    #[test]
    fn test_cqc_combination_individual_zeta_finite() {
        let peaks = vec![3.0_f64, 4.0];
        let omega = vec![10.0_f64, 30.0];
        let zeta = vec![0.05_f64, 0.05];
        let result = cqc_combination_individual_zeta(&peaks, &omega, &zeta);
        assert!(
            result.is_finite() && result >= 0.0,
            "CQC must be finite positive: {result}"
        );
    }
    #[test]
    fn test_newmark_linear_accel_step_finite() {
        let (u, v, a) = NewmarkLinearAcceleration::step(0.0, 0.0, 0.0, 1.0, 10.0, 0.05, 0.01);
        assert!(u.is_finite(), "u finite: {u}");
        assert!(v.is_finite(), "v finite: {v}");
        assert!(a.is_finite(), "a finite: {a}");
    }
    #[test]
    fn test_newmark_linear_accel_step_from_rest() {
        let (u, _v, _a) = NewmarkLinearAcceleration::step(0.0, 0.0, 0.0, 1.0, 10.0, 0.05, 0.01);
        assert!(u > 0.0, "displacement should be positive: {u}");
    }
    #[test]
    fn test_free_vibration_decay_starts_at_one() {
        let sig = free_vibration_decay(10.0, 0.05, 0.001, 100);
        assert!((sig[0] - 1.0).abs() < 1e-10, "u(0) = 1: {}", sig[0]);
    }
    #[test]
    fn test_free_vibration_decay_decays() {
        let sig = free_vibration_decay(10.0, 0.2, 0.01, 500);
        let last = sig.last().unwrap().abs();
        assert!(last < 1.0, "signal should decay: last = {last}");
    }
    #[test]
    fn test_damping_from_successive_peaks_approx() {
        let zeta_target = 0.1_f64;
        let omega_d = (1.0 - zeta_target * zeta_target).sqrt();
        let period_d = 2.0 * std::f64::consts::PI / omega_d;
        let u1 = 1.0_f64;
        let u2 = (-zeta_target * std::f64::consts::PI * 2.0 / omega_d).exp();
        let zeta_est = damping_from_successive_peaks(u1, u2);
        let _ = period_d;
        assert!(
            (zeta_est - zeta_target).abs() < 0.02,
            "zeta est = {zeta_est} vs target {zeta_target}"
        );
    }
    #[test]
    fn test_irf_starts_near_zero() {
        let h = impulse_response_function(1.0, 10.0, 0.05, 0.001, 100);
        assert!(h[0].abs() < 1e-10, "h(0) should be ~0: {}", h[0]);
    }
    #[test]
    fn test_irf_all_finite() {
        let h = impulse_response_function(2.0, 5.0, 0.1, 0.01, 200);
        for (i, &v) in h.iter().enumerate() {
            assert!(v.is_finite(), "h[{i}] must be finite: {v}");
        }
    }
    #[test]
    fn test_wilson_theta_new_zeroed() {
        let w = WilsonThetaState::new(3);
        assert_eq!(w.n_dof, 3);
        assert!((w.theta - 1.4).abs() < 1e-10);
    }
    #[test]
    fn test_wilson_theta_step_finite() {
        let mut w = WilsonThetaState::new(2);
        let m = vec![1.0_f64, 1.0];
        let k = vec![100.0_f64, 100.0];
        let c = vec![1.0_f64, 1.0];
        let f_n = vec![0.0_f64, 0.0];
        let f_np1 = vec![1.0_f64, 0.0];
        w.step_diagonal(&m, &k, &c, &f_n, &f_np1, 0.01);
        for &u in &w.u {
            assert!(u.is_finite(), "u must be finite: {u}");
        }
    }
    #[test]
    fn test_wilson_theta_step_no_force_decays() {
        let mut w = WilsonThetaState::new(1);
        w.u[0] = 1.0;
        let m = vec![1.0_f64];
        let k = vec![100.0_f64];
        let c = vec![2.0_f64];
        let f = vec![0.0_f64];
        for _ in 0..100 {
            w.step_diagonal(&m, &k, &c, &f, &f, 0.01);
        }
        assert!(w.u[0].abs() < 1.0, "Wilson-θ should dissipate: {}", w.u[0]);
    }
    #[test]
    fn test_numerical_energy_ratio_conservative() {
        let u = vec![1.0_f64, 0.5];
        let v = vec![0.2_f64, 0.3];
        let m = vec![1.0_f64, 1.0];
        let k = vec![100.0_f64, 100.0];
        let ratio = numerical_energy_ratio(&u, &v, &u, &v, &m, &k);
        assert!((ratio - 1.0).abs() < 1e-10, "same state: ratio = {ratio}");
    }
    #[test]
    fn test_numerical_energy_ratio_dissipated() {
        let u_old = vec![1.0_f64];
        let v_old = vec![1.0_f64];
        let u_new = vec![0.5_f64];
        let v_new = vec![0.5_f64];
        let m = vec![1.0_f64];
        let k = vec![1.0_f64];
        let ratio = numerical_energy_ratio(&u_new, &v_new, &u_old, &v_old, &m, &k);
        assert!(ratio < 1.0, "dissipated: ratio = {ratio}");
    }
}
