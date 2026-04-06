//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{MaxwellModel, StandardLinearSolid};

/// Simulate a creep test on a [`MaxwellModel`] under constant stress.
///
/// Returns a vector of `(time, strain)` pairs from t=0 to `t_end` with `n_steps` steps.
/// For the Maxwell model under constant stress sigma, strain grows as J(t) = 1/E + t/eta.
pub fn creep_test(
    model: &MaxwellModel,
    stress: f64,
    t_end: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let dt = t_end / n_steps as f64;
    (0..=n_steps)
        .map(|i| {
            let t = i as f64 * dt;
            let strain = stress * model.creep_compliance(t);
            (t, strain)
        })
        .collect()
}
/// Simulate a relaxation test on a [`MaxwellModel`] under constant strain.
///
/// Returns a vector of `(time, stress)` pairs from t=0 to `t_end` with `n_steps` steps.
pub fn relaxation_test(
    model: &MaxwellModel,
    strain: f64,
    t_end: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let dt = t_end / n_steps as f64;
    (0..=n_steps)
        .map(|i| {
            let t = i as f64 * dt;
            let stress = strain * model.relaxation_modulus(t);
            (t, stress)
        })
        .collect()
}
/// Compute dynamic moduli (storage, loss) over a frequency range.
///
/// Returns vectors of `(omega, E', E'')` tuples.
pub fn dynamic_moduli_sweep(
    model: &MaxwellModel,
    omega_min: f64,
    omega_max: f64,
    n_points: usize,
) -> Vec<(f64, f64, f64)> {
    let log_min = omega_min.ln();
    let log_max = omega_max.ln();
    (0..n_points)
        .map(|i| {
            let frac = if n_points > 1 {
                i as f64 / (n_points - 1) as f64
            } else {
                0.0
            };
            let omega = (log_min + frac * (log_max - log_min)).exp();
            let ep = model.storage_modulus(omega);
            let epp = model.loss_modulus(omega);
            (omega, ep, epp)
        })
        .collect()
}
/// Lanczos approximation of Γ(z) for z > 0.
pub(super) fn gamma_approx(z: f64) -> f64 {
    pub(super) const G: f64 = 7.0;
    pub(super) const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_9,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_72,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let z = z - 1.0;
    let mut x = C[0];
    for (i, &c) in C[1..].iter().enumerate() {
        x += c / (z + (i + 1) as f64);
    }
    let t = z + G + 0.5;
    (2.0 * std::f64::consts::PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
}
/// Estimate the coefficient of restitution (COR) for a viscoelastic sphere
/// impacting a rigid wall at velocity `v0`.
///
/// Uses the SLS model approximation of Hunter (1960):
///
///   COR ≈ 1 - α·(v0/v_ref)^(1/5)   (first-order correction)
///
/// where α depends on the ratio of relaxation modulus to glassy modulus.
/// This is a simplified engineering estimate.
pub fn sls_coefficient_of_restitution(sls: &StandardLinearSolid, v0: f64, v_ref: f64) -> f64 {
    let alpha_ratio = sls.e2 / (sls.e1 + sls.e2);
    let correction = alpha_ratio * (v0 / v_ref).powf(0.2);
    (1.0 - correction).max(0.0)
}
/// Estimate energy dissipation fraction for one creep cycle.
///
/// Δ = π tan(δ) / (1 + π tan(δ)) (approximate).
pub fn energy_dissipation_fraction(tan_delta: f64) -> f64 {
    let x = std::f64::consts::PI * tan_delta;
    x / (1.0 + x)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::KelvinVoigt;
    use crate::Maxwell;
    use crate::PowerLawCreep;
    use crate::viscoelastic::ArrheniusShift;
    use crate::viscoelastic::Burgers;
    use crate::viscoelastic::FractionalSls;
    use crate::viscoelastic::GeneralizedMaxwell;
    use crate::viscoelastic::KelvinVoigtChain;
    use crate::viscoelastic::KelvinVoigtModel;
    use crate::viscoelastic::Kww;
    use crate::viscoelastic::LogNormalRetardation;
    use crate::viscoelastic::MasterCurveBuilder;
    use crate::viscoelastic::NonlinearKelvinVoigt;
    use crate::viscoelastic::PronySeries;
    use crate::viscoelastic::SpringPot;
    use crate::viscoelastic::WlfShift;
    #[test]
    fn test_kv_elastic_limit() {
        let kv = KelvinVoigt::new(1.0e6, 0.0);
        let strain = 0.01;
        let strain_rate = 100.0;
        let sigma = kv.stress(strain, strain_rate);
        assert!(
            (sigma - kv.young_modulus * strain).abs() < 1.0e-10,
            "With eta=0, stress should equal E*epsilon, got {sigma}"
        );
    }
    #[test]
    fn test_kv_creep_initial_zero() {
        let kv = KelvinVoigt::new(2.0e6, 1000.0);
        let eps = kv.creep_strain(1.0e5, 0.0);
        assert!(
            eps.abs() < 1.0e-10,
            "Creep strain at t=0 should be 0, got {eps}"
        );
    }
    #[test]
    fn test_kv_creep_at_infinity() {
        let e = 2.0e6_f64;
        let kv = KelvinVoigt::new(e, 1000.0);
        let sigma0 = 1.0e5;
        let t_large = 100.0 * kv.relaxation_time();
        let eps = kv.creep_strain(sigma0, t_large);
        let expected = sigma0 / e;
        assert!(
            (eps - expected).abs() < 1.0e-10,
            "Creep strain at t->inf should be sigma0/E = {expected}, got {eps}"
        );
    }
    #[test]
    fn test_maxwell_stress_zero_at_infinity() {
        let m = Maxwell::new(1.0e6, 500.0);
        let t_large = 100.0 * m.relaxation_time();
        let sigma = m.stress_relaxation(0.01, t_large);
        assert!(
            sigma.abs() < 1.0e-10,
            "Maxwell stress at t->inf should be ~0, got {sigma}"
        );
    }
    #[test]
    fn test_maxwell_initial_stress() {
        let e = 3.0e6_f64;
        let m = Maxwell::new(e, 200.0);
        let epsilon0 = 0.005;
        let sigma = m.stress_relaxation(epsilon0, 0.0);
        let expected = e * epsilon0;
        assert!(
            (sigma - expected).abs() < 1.0e-10,
            "Maxwell stress at t=0 should be E*epsilon0 = {expected}, got {sigma}"
        );
    }
    #[test]
    fn test_sls_modulus_bounds() {
        let sls = StandardLinearSolid::new(1.0e6, 2.0e6, 500.0);
        let e_inf = sls.long_time_modulus();
        let e0 = sls.short_time_modulus();
        assert!(e_inf > 0.0, "Long-time modulus must be positive");
        assert!(e0 > 0.0, "Short-time modulus must be positive");
        assert!(
            e_inf < e0,
            "Long-time modulus {e_inf} should be less than short-time modulus {e0}"
        );
    }
    #[test]
    fn test_kv_relaxation_time() {
        let eta = 800.0_f64;
        let e = 4.0e5_f64;
        let kv = KelvinVoigt::new(e, eta);
        let tau = kv.relaxation_time();
        let expected = eta / e;
        assert!(
            (tau - expected).abs() < 1.0e-10,
            "Relaxation time should be eta/E = {expected}, got {tau}"
        );
    }
    #[test]
    fn test_maxwell_model_relaxation_modulus_at_zero() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let e_t0 = m.relaxation_modulus(0.0);
        assert!(
            (e_t0 - m.elastic_modulus).abs() < 1.0e-10,
            "E(0) should equal E = {}, got {e_t0}",
            m.elastic_modulus
        );
    }
    #[test]
    fn test_maxwell_model_relaxation_modulus_at_infinity() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let t_large = 1000.0 * m.relaxation_time();
        let e_t_inf = m.relaxation_modulus(t_large);
        assert!(
            e_t_inf.abs() < 1.0e-10,
            "E(t->inf) should be ~0, got {e_t_inf}"
        );
    }
    #[test]
    fn test_kv_model_creep_compliance_at_zero() {
        let kv = KelvinVoigtModel::new(2.0e6, 1000.0);
        let j0 = kv.creep_compliance(0.0);
        assert!(j0.abs() < 1.0e-15, "J(0) should be 0, got {j0}");
    }
    #[test]
    fn test_kv_model_creep_compliance_at_infinity() {
        let e = 2.0e6_f64;
        let kv = KelvinVoigtModel::new(e, 1000.0);
        let t_large = 1000.0 * kv.viscosity / kv.elastic_modulus;
        let j_inf = kv.creep_compliance(t_large);
        let expected = 1.0 / e;
        assert!(
            (j_inf - expected).abs() < 1.0e-15,
            "J(inf) should be 1/E = {expected}, got {j_inf}"
        );
    }
    #[test]
    fn test_generalized_maxwell_no_elements() {
        let e_inf = 5.0e5_f64;
        let gm = GeneralizedMaxwell::new(e_inf);
        for &t in &[0.0, 1.0, 100.0, 1000.0] {
            let et = gm.relaxation_modulus(t);
            assert!(
                (et - e_inf).abs() < 1.0e-10,
                "E({t}) should equal E_inf = {e_inf}, got {et}"
            );
        }
    }
    #[test]
    fn test_maxwell_stress_update_zero_strain_rate_decays() {
        let m = MaxwellModel::new(1.0e6, 200.0);
        let tau = m.relaxation_time();
        let sigma_old = 1.0e4;
        let dt = tau;
        let sigma_new = m.stress_update(0.0, 0.0, sigma_old, dt);
        let expected = sigma_old * (-1.0_f64).exp();
        assert!(
            (sigma_new - expected).abs() < 1.0e-6,
            "Stress should decay to sigma*exp(-1) = {expected}, got {sigma_new}"
        );
    }
    #[test]
    fn test_creep_test_length() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let result = creep_test(&m, 1.0e4, 10.0, 100);
        assert_eq!(
            result.len(),
            101,
            "creep_test should return n_steps+1 points"
        );
    }
    #[test]
    fn test_relaxation_test_length() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let result = relaxation_test(&m, 0.01, 10.0, 100);
        assert_eq!(
            result.len(),
            101,
            "relaxation_test should return n_steps+1 points"
        );
    }
    #[test]
    fn test_relaxation_test_initial_stress() {
        let e = 1.0e6_f64;
        let m = MaxwellModel::new(e, 500.0);
        let strain = 0.01;
        let result = relaxation_test(&m, strain, 10.0, 100);
        let (t0, sigma0) = result[0];
        assert!((t0 - 0.0).abs() < 1.0e-15);
        assert!(
            (sigma0 - strain * e).abs() < 1.0e-6,
            "Initial stress should be E*strain = {}, got {sigma0}",
            strain * e
        );
    }
    #[test]
    fn test_maxwell_storage_modulus() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let tau = m.relaxation_time();
        let omega = 1.0 / tau;
        let ep = m.storage_modulus(omega);
        assert!(
            (ep - m.elastic_modulus * 0.5).abs() < 1.0,
            "E'(omega=1/tau) should be E/2 = {}, got {ep}",
            m.elastic_modulus * 0.5
        );
    }
    #[test]
    fn test_maxwell_loss_modulus_peak() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let tau = m.relaxation_time();
        let omega = 1.0 / tau;
        let epp = m.loss_modulus(omega);
        assert!(
            (epp - m.elastic_modulus * 0.5).abs() < 1.0,
            "E''(omega=1/tau) should be E/2 = {}, got {epp}",
            m.elastic_modulus * 0.5
        );
    }
    #[test]
    fn test_maxwell_loss_tangent() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let tau = m.relaxation_time();
        let omega = 1.0 / tau;
        let tan_d = m.loss_tangent(omega);
        assert!((tan_d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_maxwell_complex_modulus() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let tau = m.relaxation_time();
        let omega = 1.0 / tau;
        let mag = m.complex_modulus_magnitude(omega);
        let expected = m.elastic_modulus / 2.0_f64.sqrt();
        assert!(
            (mag - expected).abs() < 1.0,
            "|E*| should be E/sqrt(2) = {expected}, got {mag}"
        );
    }
    #[test]
    fn test_kv_model_storage_modulus() {
        let kv = KelvinVoigtModel::new(1.0e6, 500.0);
        assert!((kv.storage_modulus(100.0) - 1.0e6).abs() < 1e-10);
    }
    #[test]
    fn test_kv_model_loss_modulus() {
        let kv = KelvinVoigtModel::new(1.0e6, 500.0);
        let omega = 2.0;
        let epp = kv.loss_modulus(omega);
        assert!((epp - 500.0 * 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_kv_chain_creep() {
        let mut chain = KelvinVoigtChain::new();
        chain.add_element(1.0e6, 1000.0);
        chain.add_element(2.0e6, 2000.0);
        assert_eq!(chain.n_elements(), 2);
        let j0 = chain.creep_compliance(0.0);
        assert!(j0.abs() < 1e-15);
        let t_large = 1e6;
        let j_inf = chain.creep_compliance(t_large);
        let expected = 1.0 / 1.0e6 + 1.0 / 2.0e6;
        assert!((j_inf - expected).abs() < 1e-10);
    }
    #[test]
    fn test_kv_chain_retardation_times() {
        let mut chain = KelvinVoigtChain::new();
        chain.add_element(1.0e6, 1000.0);
        chain.add_element(2.0e6, 2000.0);
        let taus = chain.retardation_times();
        assert_eq!(taus.len(), 2);
        assert!((taus[0] - 1000.0 / 1.0e6).abs() < 1e-15);
        assert!((taus[1] - 2000.0 / 2.0e6).abs() < 1e-15);
    }
    #[test]
    fn test_generalized_maxwell_storage_loss() {
        let mut gm = GeneralizedMaxwell::new(1.0e5);
        gm.add_element(MaxwellModel::new(1.0e6, 500.0), 1.0);
        let ep_zero = gm.storage_modulus(1e-10);
        assert!((ep_zero - 1.0e5).abs() < 1.0);
        let epp_zero = gm.loss_modulus(1e-10);
        assert!(epp_zero.abs() < 1.0);
    }
    #[test]
    fn test_generalized_maxwell_glassy_modulus() {
        let mut gm = GeneralizedMaxwell::new(1.0e5);
        gm.add_element(MaxwellModel::new(1.0e6, 500.0), 1.0);
        gm.add_element(MaxwellModel::new(2.0e6, 1000.0), 0.5);
        let e0 = gm.glassy_modulus();
        assert!((e0 - 2.1e6).abs() < 1.0);
    }
    #[test]
    fn test_generalized_maxwell_loss_tangent() {
        let mut gm = GeneralizedMaxwell::new(1.0e5);
        gm.add_element(MaxwellModel::new(1.0e6, 500.0), 1.0);
        let tan_d = gm.loss_tangent(1.0);
        assert!(tan_d >= 0.0);
    }
    #[test]
    fn test_prony_series() {
        let mut ps = PronySeries::new(1.0e5);
        ps.add_term(5.0e5, 1.0);
        ps.add_term(3.0e5, 10.0);
        let e0 = ps.relaxation_modulus(0.0);
        assert!((e0 - 9.0e5).abs() < 1.0);
        let e_inf = ps.relaxation_modulus(1e10);
        assert!((e_inf - 1.0e5).abs() < 1.0);
    }
    #[test]
    fn test_prony_series_storage_loss() {
        let mut ps = PronySeries::new(1.0e5);
        ps.add_term(5.0e5, 1.0);
        let ep = ps.storage_modulus(1.0);
        let epp = ps.loss_modulus(1.0);
        assert!(ep > 0.0);
        assert!(epp > 0.0);
    }
    #[test]
    fn test_prony_to_generalized_maxwell() {
        let mut ps = PronySeries::new(1.0e5);
        ps.add_term(5.0e5, 1.0);
        let gm = ps.to_generalized_maxwell();
        assert_eq!(gm.n_elements(), 1);
        assert!((gm.spring_inf - 1.0e5).abs() < 1e-10);
        let e_ps = ps.relaxation_modulus(0.5);
        let e_gm = gm.relaxation_modulus(0.5);
        assert!((e_ps - e_gm).abs() < 1.0);
    }
    #[test]
    fn test_wlf_shift_at_reference() {
        let wlf = WlfShift::new(17.44, 51.6, 25.0);
        let a = wlf.shift_factor(25.0);
        assert!((a - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_wlf_shift_above_reference() {
        let wlf = WlfShift::new(17.44, 51.6, 25.0);
        let a = wlf.shift_factor(50.0);
        assert!(a < 1.0);
        assert!(a > 0.0);
    }
    #[test]
    fn test_wlf_shift_below_reference() {
        let wlf = WlfShift::new(17.44, 51.6, 25.0);
        let a = wlf.shift_factor(10.0);
        assert!(a > 1.0);
    }
    #[test]
    fn test_wlf_shifted_frequency() {
        let wlf = WlfShift::new(17.44, 51.6, 25.0);
        let omega = 10.0;
        let omega_shifted = wlf.shifted_frequency(omega, 25.0);
        assert!((omega_shifted - omega).abs() < 1e-10);
    }
    #[test]
    fn test_arrhenius_shift_at_reference() {
        let arr = ArrheniusShift::new(50000.0, 300.0);
        let a = arr.shift_factor(300.0);
        assert!((a - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_arrhenius_shift_above_reference() {
        let arr = ArrheniusShift::new(50000.0, 300.0);
        let a = arr.shift_factor(350.0);
        assert!(a < 1.0);
        assert!(a > 0.0);
    }
    #[test]
    fn test_arrhenius_shift_below_reference() {
        let arr = ArrheniusShift::new(50000.0, 300.0);
        let a = arr.shift_factor(250.0);
        assert!(a > 1.0);
    }
    #[test]
    fn test_dynamic_moduli_sweep() {
        let m = MaxwellModel::new(1.0e6, 500.0);
        let sweep = dynamic_moduli_sweep(&m, 0.01, 100.0, 50);
        assert_eq!(sweep.len(), 50);
        assert!(sweep.last().unwrap().1 > sweep.first().unwrap().1);
    }
    #[test]
    fn test_sls_relaxation_modulus() {
        let sls = StandardLinearSolid::new(1.0e6, 2.0e6, 500.0);
        let e0 = sls.relaxation_modulus(0.0);
        assert!((e0 - 3.0e6).abs() < 1.0);
        let e_inf = sls.relaxation_modulus(1e10);
        assert!((e_inf - 1.0e6).abs() < 1.0);
    }
    #[test]
    fn test_sls_creep_compliance() {
        let sls = StandardLinearSolid::new(1.0e6, 2.0e6, 500.0);
        let j0 = sls.creep_compliance(0.0);
        let expected_j0 = 1.0 / 3.0e6;
        assert!((j0 - expected_j0).abs() < 1e-15);
        let j_inf = sls.creep_compliance(1e10);
        let expected_j_inf = 1.0 / 1.0e6;
        assert!((j_inf - expected_j_inf).abs() < 1e-15);
    }
    #[test]
    fn test_generalized_maxwell_n_elements() {
        let mut gm = GeneralizedMaxwell::new(1.0e5);
        assert_eq!(gm.n_elements(), 0);
        gm.add_element(MaxwellModel::new(1.0e6, 500.0), 1.0);
        assert_eq!(gm.n_elements(), 1);
    }
    #[test]
    fn test_springpot_elastic_limit_alpha_zero() {
        let sp = SpringPot::new(1.0e6, 0.0);
        let ep = sp.storage_modulus(10.0);
        let epp = sp.loss_modulus(10.0);
        assert!((ep - 1.0e6).abs() < 1.0, "E'={ep}");
        assert!(epp.abs() < 1.0e-6, "E''={epp}");
    }
    #[test]
    fn test_springpot_viscous_limit_alpha_one() {
        let sp = SpringPot::new(1.0e3, 1.0);
        let epp = sp.loss_modulus(100.0);
        let ep = sp.storage_modulus(100.0);
        assert!(ep.abs() < 1.0e-6, "E'={ep}");
        assert!(epp > 1.0, "E''={epp}");
    }
    #[test]
    fn test_springpot_loss_tangent_frequency_independent() {
        let sp = SpringPot::new(1.0e5, 0.5);
        let td1 = sp.storage_modulus(1.0) / sp.loss_modulus(1.0);
        let td2 = sp.storage_modulus(100.0) / sp.loss_modulus(100.0);
        assert!(
            (td1 - td2).abs() < 1e-10,
            "tan δ should be frequency-independent"
        );
    }
    #[test]
    fn test_springpot_complex_modulus_magnitude() {
        let sp = SpringPot::new(1.0e5, 0.4);
        let omega = 5.0;
        let mag = sp.complex_modulus_magnitude(omega);
        let expected = 1.0e5 * omega.powf(0.4);
        assert!((mag - expected).abs() < 1.0, "magnitude={mag}");
    }
    #[test]
    fn test_springpot_phase_angle() {
        let alpha = 0.6;
        let sp = SpringPot::new(1.0e5, alpha);
        let phase = sp.phase_angle();
        let expected = alpha * std::f64::consts::FRAC_PI_2;
        assert!((phase - expected).abs() < 1e-12);
    }
    #[test]
    fn test_fractional_sls_glassy_modulus() {
        let sp = SpringPot::new(1.0e4, 0.5);
        let fsls = FractionalSls::new(1.0e6, 2.0e6, sp);
        assert!((fsls.glassy_modulus() - 3.0e6).abs() < 1.0);
    }
    #[test]
    fn test_fractional_sls_loss_modulus_positive() {
        let sp = SpringPot::new(1.0e4, 0.5);
        let fsls = FractionalSls::new(1.0e5, 5.0e5, sp);
        let epp = fsls.loss_modulus(10.0);
        assert!(epp >= 0.0, "E''={epp}");
    }
    #[test]
    fn test_fractional_sls_storage_modulus_positive() {
        let sp = SpringPot::new(1.0e4, 0.5);
        let fsls = FractionalSls::new(1.0e5, 5.0e5, sp);
        let ep = fsls.storage_modulus(10.0);
        assert!(ep > 0.0, "E'={ep}");
    }
    #[test]
    fn test_fractional_sls_low_freq_approaches_e_inf() {
        let sp = SpringPot::new(1.0e4, 0.5);
        let e_inf = 1.0e5_f64;
        let fsls = FractionalSls::new(e_inf, 5.0e5, sp);
        let ep_lo = fsls.storage_modulus(1e-6);
        assert!((ep_lo - e_inf).abs() / e_inf < 0.01, "E'(ω→0)={ep_lo}");
    }
    #[test]
    fn test_power_law_creep_initial() {
        let m = PowerLawCreep::new(1.0e-11, 1.0e-14, 0.3);
        assert!((m.creep_compliance(0.0) - 1.0e-11).abs() < 1e-25);
    }
    #[test]
    fn test_power_law_creep_increases_with_time() {
        let m = PowerLawCreep::new(1.0e-11, 1.0e-14, 0.3);
        assert!(m.creep_compliance(1000.0) > m.creep_compliance(1.0));
    }
    #[test]
    fn test_power_law_creep_rate_positive() {
        let m = PowerLawCreep::new(1.0e-11, 1.0e-14, 0.5);
        let rate = m.creep_rate(10.0);
        assert!(rate > 0.0, "creep rate={rate}");
    }
    #[test]
    fn test_power_law_creep_strain_proportional_to_stress() {
        let m = PowerLawCreep::new(1.0e-11, 1.0e-14, 0.4);
        let sigma = 1.0e6;
        let eps = m.creep_strain(sigma, 100.0);
        let eps2 = m.creep_strain(2.0 * sigma, 100.0);
        assert!(
            (eps2 - 2.0 * eps).abs() < 1e-20,
            "not linear: {eps} vs {eps2}"
        );
    }
    #[test]
    fn test_nonlinear_kv_m1_reduces_to_linear() {
        let nlkv = NonlinearKelvinVoigt::new(1.0e6, 500.0, 1.0);
        let kv = KelvinVoigt::new(1.0e6, 500.0);
        let strain = 0.01;
        let rate = 0.1;
        let sigma_nl = nlkv.stress(strain, rate);
        let sigma_lin = kv.stress(strain, rate);
        assert!(
            (sigma_nl - sigma_lin).abs() < 1.0,
            "nl={sigma_nl}, lin={sigma_lin}"
        );
    }
    #[test]
    fn test_nonlinear_kv_stress_increases_with_rate() {
        let nlkv = NonlinearKelvinVoigt::new(0.0, 500.0, 2.0);
        let s1 = nlkv.stress(0.0, 0.1);
        let s2 = nlkv.stress(0.0, 1.0);
        assert!(s2 > s1, "s1={s1}, s2={s2}");
    }
    #[test]
    fn test_nonlinear_kv_negative_rate_sign() {
        let nlkv = NonlinearKelvinVoigt::new(0.0, 100.0, 1.5);
        let sigma_pos = nlkv.stress(0.0, 1.0);
        let sigma_neg = nlkv.stress(0.0, -1.0);
        assert!(
            (sigma_pos + sigma_neg).abs() < 1e-10,
            "should be antisymmetric"
        );
    }
    #[test]
    fn test_burgers_creep_initial() {
        let b = Burgers::new(1.0e6, 1.0e9, 5.0e5, 2.0e8);
        let j0 = b.creep_compliance(0.0);
        let expected = 1.0 / 1.0e6;
        assert!(
            (j0 - expected).abs() < 1e-15,
            "J(0)={j0}, expected={expected}"
        );
    }
    #[test]
    fn test_burgers_creep_increases_with_time() {
        let b = Burgers::new(1.0e6, 1.0e9, 5.0e5, 2.0e8);
        assert!(b.creep_compliance(1.0) > b.creep_compliance(0.0));
        assert!(b.creep_compliance(10.0) > b.creep_compliance(1.0));
    }
    #[test]
    fn test_burgers_steady_state_creep_rate() {
        let eta_m = 1.0e9_f64;
        let sigma = 1.0e6;
        let b = Burgers::new(1.0e6, eta_m, 5.0e5, 2.0e8);
        let rate = b.steady_state_creep_rate(sigma);
        assert!((rate - sigma / eta_m).abs() < 1e-20, "rate={rate}");
    }
    #[test]
    fn test_burgers_retardation_time() {
        let b = Burgers::new(1.0e6, 1.0e9, 4.0e5, 8.0e7);
        let tau = b.retardation_time();
        assert!((tau - 8.0e7 / 4.0e5).abs() < 1e-10, "tau={tau}");
    }
    #[test]
    fn test_kww_initial_modulus() {
        let kww = Kww::new(3.0e9, 1.0e8, 100.0, 0.6);
        let e0 = kww.relaxation_modulus(0.0);
        assert!((e0 - 3.0e9).abs() < 1.0, "E(0)={e0}");
    }
    #[test]
    fn test_kww_long_time_modulus() {
        let e_inf = 1.0e8_f64;
        let kww = Kww::new(3.0e9, e_inf, 100.0, 0.6);
        let e_long = kww.relaxation_modulus(1e15);
        assert!((e_long - e_inf).abs() < 1.0, "E(∞)={e_long}");
    }
    #[test]
    fn test_kww_modulus_decreases_monotonically() {
        let kww = Kww::new(3.0e9, 1.0e8, 100.0, 0.6);
        let e1 = kww.relaxation_modulus(10.0);
        let e2 = kww.relaxation_modulus(1000.0);
        assert!(e1 > e2, "E(10)={e1} should be > E(1000)={e2}");
    }
    #[test]
    fn test_kww_beta_one_matches_maxwell() {
        let e0 = 1.0e9_f64;
        let e_inf = 0.0;
        let tau = 50.0_f64;
        let kww = Kww::new(e0, e_inf, tau, 1.0);
        let m = MaxwellModel::new(e0, e0 * tau);
        let t = 25.0;
        let e_kww = kww.relaxation_modulus(t);
        let e_maxwell = m.relaxation_modulus(t);
        assert!(
            (e_kww - e_maxwell).abs() / e0 < 1e-10,
            "KWW={e_kww}, Maxwell={e_maxwell}"
        );
    }
    #[test]
    fn test_kww_half_relaxation_time() {
        let kww = Kww::new(3.0e9, 0.0, 100.0, 0.5);
        let t_half = kww.half_relaxation_time();
        let phi_half = kww.phi(t_half);
        assert!((phi_half - 0.5).abs() < 1e-10, "φ(t₁/₂)={phi_half}");
    }
    #[test]
    fn test_log_normal_retardation_spectrum_peak() {
        let lnr = LogNormalRetardation::new(1e-11, 1e-9, 10.0, 1.0);
        let at_peak = lnr.spectrum(10.0);
        let off_peak = lnr.spectrum(100.0);
        assert!(
            at_peak > off_peak,
            "spectrum peak={at_peak}, off-peak={off_peak}"
        );
    }
    #[test]
    fn test_log_normal_compliance_increases_with_time() {
        let lnr = LogNormalRetardation::new(1e-11, 1e-9, 10.0, 1.0);
        let j1 = lnr.creep_compliance(1.0);
        let j2 = lnr.creep_compliance(100.0);
        assert!(j2 > j1, "J(100)={j2} should be > J(1)={j1}");
    }
    #[test]
    fn test_log_normal_loss_compliance_positive() {
        let lnr = LogNormalRetardation::new(1e-11, 1e-9, 10.0, 1.0);
        let jpp = lnr.loss_compliance(1.0 / 10.0);
        assert!(jpp > 0.0, "J''={jpp}");
    }
    #[test]
    fn test_master_curve_builder_add_and_interpolate() {
        let wlf = WlfShift::new(17.44, 51.6, 25.0);
        let mut builder = MasterCurveBuilder::new(25.0);
        builder.add_point_wlf(&wlf, 25.0, 1.0, 3.0e9);
        builder.add_point_wlf(&wlf, 25.0, 10.0, 2.0e9);
        builder.add_point_wlf(&wlf, 25.0, 100.0, 1.0e9);
        builder.sort();
        assert_eq!(builder.n_points(), 3);
        let e = builder.interpolate(1.0);
        assert!(e.is_some(), "interpolation failed");
        let e = e.unwrap();
        assert!((e - 2.0e9).abs() < 1.0, "interpolated E={e}");
    }
    #[test]
    fn test_master_curve_builder_arrhenius() {
        let arr = ArrheniusShift::new(50_000.0, 300.0);
        let mut builder = MasterCurveBuilder::new(300.0);
        builder.add_point_arrhenius(&arr, 300.0, 1.0, 1.0e9);
        builder.add_point_arrhenius(&arr, 300.0, 10.0, 8.0e8);
        builder.sort();
        assert_eq!(builder.n_points(), 2);
    }
    #[test]
    fn test_sls_cor_less_than_one() {
        let sls = StandardLinearSolid::new(1.0e6, 2.0e6, 500.0);
        let cor = sls_coefficient_of_restitution(&sls, 1.0, 1.0);
        assert!((0.0..=1.0).contains(&cor), "COR={cor}");
    }
    #[test]
    fn test_energy_dissipation_fraction_positive() {
        let f = energy_dissipation_fraction(0.3);
        assert!(f > 0.0 && f < 1.0, "dissipation fraction={f}");
    }
    #[test]
    fn test_gamma_approx_half() {
        let g = gamma_approx(0.5);
        assert!((g - std::f64::consts::PI.sqrt()).abs() < 1e-6, "Γ(0.5)={g}");
    }
    #[test]
    fn test_gamma_approx_one() {
        let g = gamma_approx(1.0);
        assert!((g - 1.0).abs() < 1e-10, "Γ(1)={g}");
    }
    #[test]
    fn test_gamma_approx_two() {
        let g = gamma_approx(2.0);
        assert!((g - 1.0).abs() < 1e-10, "Γ(2)={g}");
    }
    #[test]
    fn test_gamma_approx_three() {
        let g = gamma_approx(3.0);
        assert!((g - 2.0).abs() < 1e-8, "Γ(3)={g}");
    }
    #[test]
    fn test_prony_compute_relaxation_modulus_t0() {
        let mut ps = PronySeries::new(1.0e6);
        ps.add_term(2.0e6, 1.0);
        ps.add_term(0.5e6, 0.1);
        let g0 = ps.compute_relaxation_modulus(0.0);
        let expected = 3.5e6;
        assert!(
            (g0 - expected).abs() < 1e-6,
            "G(0)={g0}, expected={expected}"
        );
    }
    #[test]
    fn test_prony_compute_relaxation_modulus_t_inf() {
        let mut ps = PronySeries::new(0.5e6);
        ps.add_term(1.0e6, 0.01);
        let g_large = ps.compute_relaxation_modulus(1000.0);
        assert!((g_large - 0.5e6).abs() < 1.0, "G(∞)≈G_inf: got {g_large}");
    }
    #[test]
    fn test_prony_relaxation_modulus_decreases_with_time() {
        let mut ps = PronySeries::new(1.0e6);
        ps.add_term(3.0e6, 1.0);
        let g1 = ps.compute_relaxation_modulus(0.5);
        let g2 = ps.compute_relaxation_modulus(2.0);
        assert!(g1 > g2, "G should decrease with time: {g1} vs {g2}");
    }
    #[test]
    fn test_maxwell_dynamic_modulus_high_freq() {
        let m = MaxwellModel::new(1.0e6, 1000.0);
        let (gp, gpp) = m.compute_dynamic_modulus(1.0e6);
        assert!(gp > gpp, "at high freq G' > G'': {gp} vs {gpp}");
        assert!(
            (gp - 1.0e6).abs() / 1.0e6 < 0.01,
            "G' ≈ E at high freq: {gp}"
        );
    }
    #[test]
    fn test_maxwell_dynamic_modulus_low_freq() {
        let m = MaxwellModel::new(1.0e6, 1.0);
        let (gp, gpp) = m.compute_dynamic_modulus(1.0e-3);
        assert!(gp < 1.0, "G' small at low freq: {gp}");
        assert!(gpp < 1.0, "G'' small at low freq: {gpp}");
    }
    #[test]
    fn test_maxwell_dynamic_modulus_peak_loss() {
        let e = 2.0e6_f64;
        let tau = 0.5_f64;
        let m = MaxwellModel::new(e, e * tau);
        let omega_peak = 1.0 / tau;
        let (_, gpp) = m.compute_dynamic_modulus(omega_peak);
        let expected = e / 2.0;
        assert!(
            (gpp - expected).abs() / expected < 1e-10,
            "G''_max={gpp}, expected={expected}"
        );
    }
    #[test]
    fn test_tan_delta_elastic_limit() {
        let gm = GeneralizedMaxwell::new(1.0e6);
        let td = gm.compute_tan_delta(1.0);
        assert!(td.abs() < 1e-15, "elastic-only tan δ = 0, got {td}");
    }
    #[test]
    fn test_tan_delta_positive() {
        let mut gm = GeneralizedMaxwell::new(0.5e6);
        gm.add_element(MaxwellModel::new(2.0e6, 2.0e6 * 0.1), 1.0);
        let td = gm.compute_tan_delta(10.0);
        assert!(td > 0.0, "tan δ should be positive: {td}");
    }
    #[test]
    fn test_kelvin_voigt_creep_compliance_t0() {
        let kv = KelvinVoigt::new(1.0e6, 1000.0);
        let j0 = kv.compute_creep_compliance(0.0);
        assert!(j0.abs() < 1e-30, "J(0) = 0 for KV, got {j0}");
    }
    #[test]
    fn test_kelvin_voigt_creep_compliance_saturates() {
        let e = 2.0e6_f64;
        let kv = KelvinVoigt::new(e, e * 0.001);
        let j_large = kv.compute_creep_compliance(100.0);
        let expected = 1.0 / e;
        assert!(
            (j_large - expected).abs() / expected < 1e-6,
            "J(∞) ≈ 1/E: got {j_large}, expected {expected}"
        );
    }
}
