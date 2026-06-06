//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod viscoelastic_extended_tests {

    use crate::viscoelastic_fem::*;
    #[test]
    fn test_maxwell_relaxation_time_correct() {
        let m = MaxwellModel::new(1e9, 1e6);
        let tau = m.relaxation_time();
        assert!((tau - 1e6 / 1e9).abs() < 1e-20, "tau = eta / E = 1e-3 s");
    }
    #[test]
    fn test_maxwell_stress_relaxation_at_t0() {
        let m = MaxwellModel::new(1e9, 1e7);
        let sigma = m.stress_relaxation(0.001, 0.0);
        assert!(
            (sigma - 1e9 * 0.001).abs() < 1e6,
            "sigma(0) = E * eps0 = 1e6 Pa"
        );
    }
    #[test]
    fn test_maxwell_stress_relaxation_decays() {
        let m = MaxwellModel::new(1e9, 1e7);
        let s0 = m.stress_relaxation(0.001, 0.0);
        let s1 = m.stress_relaxation(0.001, 1.0);
        let s_inf = m.stress_relaxation(0.001, 1e12);
        assert!(s1 < s0, "stress must decay: s(1s)={s1} < s(0)={s0}");
        assert!(
            s_inf < 1e-3 * s0,
            "stress must relax to near zero at long time: {s_inf}"
        );
    }
    #[test]
    fn test_maxwell_relaxation_modulus_decays() {
        let m = MaxwellModel::new(1e9, 1e7);
        let e0 = m.relaxation_modulus(0.0);
        let e1 = m.relaxation_modulus(1.0);
        assert!(e1 < e0, "relaxation modulus must decay");
    }
    #[test]
    fn test_maxwell_creep_compliance_increases() {
        let m = MaxwellModel::new(1e9, 1e7);
        let j0 = m.creep_compliance(0.0);
        let j1 = m.creep_compliance(1.0);
        let j10 = m.creep_compliance(10.0);
        assert!(
            j1 > j0 && j10 > j1,
            "Maxwell creep compliance must increase"
        );
    }
    #[test]
    fn test_maxwell_storage_modulus_zero_at_low_freq() {
        let m = MaxwellModel::new(1e9, 1e7);
        let e_prime = m.storage_modulus(1e-10);
        assert!(e_prime < 1e3, "Maxwell E'(omega→0) → 0: {e_prime}");
    }
    #[test]
    fn test_maxwell_storage_modulus_approaches_e_at_high_freq() {
        let m = MaxwellModel::new(1e9, 1e7);
        let e_prime = m.storage_modulus(1e15);
        assert!(
            (e_prime / m.modulus - 1.0).abs() < 1e-3,
            "Maxwell E'(omega→inf) → E: {e_prime}"
        );
    }
    #[test]
    fn test_maxwell_loss_tangent_decreases_at_high_freq() {
        let m = MaxwellModel::new(1e9, 1e7);
        let td1 = m.loss_tangent(1.0);
        let td_high = m.loss_tangent(1e10);
        assert!(
            td_high < td1,
            "loss tangent should decrease at high frequency"
        );
    }
    #[test]
    fn test_kv_retardation_time_correct() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let tau = kv.retardation_time();
        assert!((tau - 1e6 / 1e9).abs() < 1e-20, "tau_KV = eta / E");
    }
    #[test]
    fn test_kv_creep_compliance_zero_at_t0() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let j0 = kv.creep_compliance(0.0);
        assert!(
            j0.abs() < 1e-40,
            "KV J(0) = 0 (spring is undisplaced): {j0}"
        );
    }
    #[test]
    fn test_kv_creep_compliance_increases_monotonically() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let times = [0.0, 0.001, 0.01, 0.1, 1.0, 10.0];
        let jv: Vec<f64> = times.iter().map(|&t| kv.creep_compliance(t)).collect();
        for w in jv.windows(2) {
            assert!(w[1] >= w[0] - 1e-30, "KV compliance must be non-decreasing");
        }
    }
    #[test]
    fn test_kv_creep_compliance_approaches_1_over_e() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let j_inf = kv.creep_compliance(1e15);
        let expected = 1.0 / kv.modulus;
        assert!(
            (j_inf / expected - 1.0).abs() < 1e-4,
            "KV J(∞) → 1/E: {j_inf} vs {expected}"
        );
    }
    #[test]
    fn test_kv_storage_modulus_is_constant() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let e_prime_1 = kv.storage_modulus(1.0);
        let e_prime_2 = kv.storage_modulus(100.0);
        assert!(
            (e_prime_1 - e_prime_2).abs() < 1e-6,
            "KV E' is constant = E"
        );
        assert!((e_prime_1 - kv.modulus).abs() < 1e-6);
    }
    #[test]
    fn test_kv_loss_modulus_increases_with_freq() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let e_pp_1 = kv.loss_modulus(1.0);
        let e_pp_2 = kv.loss_modulus(10.0);
        assert!(
            e_pp_2 > e_pp_1,
            "KV E'' = eta*omega should increase with freq"
        );
    }
    #[test]
    fn test_kv_loss_tangent_increases_with_freq() {
        let kv = KelvinVoigtModel::new(1e9, 1e6);
        let td1 = kv.loss_tangent(1.0);
        let td2 = kv.loss_tangent(10.0);
        assert!(td2 > td1, "KV tan_delta should increase with frequency");
    }
    #[test]
    fn test_sls_equilibrium_modulus_correct() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let e_inf = sls.equilibrium_modulus();
        let expected = 2e9 * 3e9 / (2e9 + 3e9);
        assert!(
            (e_inf - expected).abs() / expected < 1e-10,
            "SLS E_inf: {e_inf} vs {expected}"
        );
    }
    #[test]
    fn test_sls_instantaneous_modulus_is_e1() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        assert!((sls.instantaneous_modulus() - 2e9).abs() < 1.0);
    }
    #[test]
    fn test_sls_relaxation_modulus_at_t0() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let e0 = sls.relaxation_modulus(0.0);
        assert!(
            (e0 - sls.instantaneous_modulus()).abs() / sls.instantaneous_modulus() < 1e-10,
            "SLS E(0) = E_0: {e0}"
        );
    }
    #[test]
    fn test_sls_relaxation_modulus_decays_to_e_inf() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let e_late = sls.relaxation_modulus(1e15);
        let e_inf = sls.equilibrium_modulus();
        assert!(
            (e_late - e_inf).abs() / e_inf < 1e-4,
            "SLS E(∞) → E_inf: {e_late} vs {e_inf}"
        );
    }
    #[test]
    fn test_sls_relaxation_modulus_monotone_decreasing() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let times = [0.0, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0];
        let vals: Vec<f64> = times.iter().map(|&t| sls.relaxation_modulus(t)).collect();
        for w in vals.windows(2) {
            assert!(w[1] <= w[0] + 1e-3, "SLS relaxation must be non-increasing");
        }
    }
    #[test]
    fn test_sls_creep_compliance_at_t0() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let j0 = sls.creep_compliance(0.0);
        let expected = 1.0 / sls.instantaneous_modulus();
        assert!(
            (j0 - expected).abs() / expected < 1e-10,
            "SLS J(0) = 1/E_0: {j0}"
        );
    }
    #[test]
    fn test_sls_storage_modulus_bounds() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let e_inf = sls.equilibrium_modulus();
        let e0 = sls.instantaneous_modulus();
        let e_prime_mid = sls.storage_modulus(1.0 / sls.relaxation_time());
        assert!(e_prime_mid >= e_inf, "SLS E' >= E_inf");
        assert!(e_prime_mid <= e0, "SLS E' <= E_0");
    }
    #[test]
    fn test_sls_loss_tangent_peak_near_1_over_tau() {
        let sls = StandardLinearSolid::new(2e9, 3e9, 1e7);
        let tau_r = sls.relaxation_time();
        let omega_peak = 1.0 / tau_r;
        let td_peak = sls.loss_tangent(omega_peak);
        let td_lower = sls.loss_tangent(omega_peak * 0.01);
        let td_higher = sls.loss_tangent(omega_peak * 100.0);
        assert!(
            td_peak > td_lower && td_peak > td_higher,
            "SLS loss tangent peak near 1/tau: peak={td_peak}, lower={td_lower}, higher={td_higher}"
        );
    }
    #[test]
    fn test_shear_prony_instantaneous_modulus_is_sum() {
        let gp = ShearProny::rubber();
        let g0 = gp.instantaneous_modulus();
        let expected = gp.g_inf + gp.terms.iter().map(|(g, _)| g).sum::<f64>();
        assert!(
            (g0 - expected).abs() / expected < 1e-10,
            "G_0 = G_inf + sum(G_i)"
        );
    }
    #[test]
    fn test_shear_prony_relaxation_at_t0() {
        let gp = ShearProny::rubber();
        let g0 = gp.relaxation_modulus(0.0);
        assert!((g0 - gp.instantaneous_modulus()).abs() / gp.instantaneous_modulus() < 1e-10);
    }
    #[test]
    fn test_shear_prony_relaxation_decays_to_g_inf() {
        let gp = ShearProny::rubber();
        let g_late = gp.relaxation_modulus(1e15);
        assert!(
            (g_late - gp.g_inf).abs() / gp.g_inf < 1e-3,
            "G(∞) → G_inf: {g_late} vs {}",
            gp.g_inf
        );
    }
    #[test]
    fn test_shear_prony_loss_tangent_nonnegative() {
        let gp = ShearProny::rubber();
        for omega in [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0] {
            let td = gp.loss_tangent(omega);
            assert!(
                td >= 0.0,
                "loss tangent must be non-negative at omega={omega}"
            );
        }
    }
    #[test]
    fn test_shear_prony_peak_loss_frequency_positive() {
        let gp = ShearProny::rubber();
        let f_peak = gp.peak_loss_frequency();
        assert!(
            f_peak > 0.0 && f_peak.is_finite(),
            "peak loss frequency = {f_peak}"
        );
    }
    #[test]
    fn test_hereditary_integral_ramp_stress_positive() {
        let mat = ViscoelasticMaterial::polymer();
        let dt = 0.1;
        let n = 20;
        let strain_history: Vec<f64> = (0..n).map(|i| i as f64 * 1e-4).collect();
        let sigma = hereditary_stress_integral(&mat, &strain_history, dt);
        assert!(
            sigma > 0.0,
            "ramp strain → positive stress via hereditary integral: {sigma}"
        );
    }
    #[test]
    fn test_hereditary_integral_zero_strain_zero_stress() {
        let mat = ViscoelasticMaterial::polymer();
        let dt = 0.1;
        let n = 10;
        let strain_history = vec![0.0; n];
        let sigma = hereditary_stress_integral(&mat, &strain_history, dt);
        assert!(sigma.abs() < 1e-30, "zero strain → zero stress: {sigma}");
    }
    #[test]
    fn test_frequency_sweep_correct_length() {
        let mat = ViscoelasticMaterial::polymer();
        let result = frequency_sweep(&mat, 0.01, 100.0, 20);
        assert_eq!(result.len(), 20);
    }
    #[test]
    fn test_frequency_sweep_first_last_omega() {
        let mat = ViscoelasticMaterial::polymer();
        let result = frequency_sweep(&mat, 0.01, 100.0, 10);
        assert!(
            (result[0].omega - 0.01).abs() / 0.01 < 1e-10,
            "first omega should match omega_min"
        );
        assert!(
            (result[9].omega - 100.0).abs() / 100.0 < 1e-10,
            "last omega should match omega_max"
        );
    }
    #[test]
    fn test_frequency_sweep_storage_in_bounds() {
        let mat = ViscoelasticMaterial::polymer();
        let result = frequency_sweep(&mat, 0.001, 1000.0, 30);
        let e0 = mat.instantaneous_modulus();
        let e_inf = mat.e_inf;
        for pt in &result {
            assert!(
                (e_inf - 1.0..=e0 + 1.0).contains(&pt.storage),
                "storage modulus {:.3e} out of bounds [{e_inf:.3e}, {e0:.3e}]",
                pt.storage
            );
        }
    }
    #[test]
    fn test_frequency_sweep_tan_delta_nonnegative() {
        let mat = ViscoelasticMaterial::polymer();
        let result = frequency_sweep(&mat, 0.01, 100.0, 15);
        for pt in &result {
            assert!(
                pt.tan_delta >= 0.0,
                "tan_delta must be non-negative at omega={}",
                pt.omega
            );
        }
    }
    #[test]
    fn test_frequency_sweep_magnitude_geq_storage() {
        let mat = ViscoelasticMaterial::polymer();
        let result = frequency_sweep(&mat, 0.1, 10.0, 10);
        for pt in &result {
            assert!(
                pt.magnitude >= pt.storage - 1.0,
                "|E*| >= E' at omega={}",
                pt.omega
            );
        }
    }
}
