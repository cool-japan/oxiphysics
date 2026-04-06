//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    CobleCreep, LarsonMillerParameter, MonkmanGrant, NabarroHerringCreep, NortonCreep,
    OrrSherbyDornParameter, RuptureLifePrediction,
};

/// Universal gas constant R (J/mol/K).
pub(super) const GAS_CONSTANT: f64 = 8.314;
/// Determine which diffusion creep mechanism dominates.
///
/// Returns "NH" for Nabarro-Herring or "Coble" depending on which
/// gives a higher strain rate.
pub fn dominant_diffusion_creep(
    nh: &NabarroHerringCreep,
    coble: &CobleCreep,
    stress: f64,
    temperature: f64,
) -> &'static str {
    let rate_nh = nh.strain_rate(stress, temperature);
    let rate_coble = coble.strain_rate(stress, temperature);
    if rate_nh > rate_coble { "NH" } else { "Coble" }
}
/// Creep-relaxation coupling: stress relaxation under creep conditions.
///
/// In a fixed-strain scenario, the stress relaxes as creep accommodates
/// the imposed strain. The rate equation is:
///
/// d(sigma)/dt = -E * eps_dot_creep(sigma, T)
pub fn creep_stress_relaxation(
    norton: &NortonCreep,
    sigma_0: f64,
    temperature: f64,
    young_modulus: f64,
    dt: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let mut sigma = sigma_0;
    let mut history = Vec::with_capacity(n_steps + 1);
    history.push((0.0, sigma));
    for i in 1..=n_steps {
        let eps_dot = norton.creep_strain_rate(sigma, temperature);
        sigma -= young_modulus * eps_dot * dt;
        if sigma < 0.0 {
            sigma = 0.0;
        }
        history.push((i as f64 * dt, sigma));
    }
    history
}
/// Simulate creep deformation under a variable temperature history.
///
/// Given a stress (constant) and a vector of (time, temperature) pairs,
/// returns the accumulated creep strain history.
pub fn variable_temperature_creep(
    norton: &NortonCreep,
    stress: f64,
    temp_history: &[(f64, f64)],
) -> Vec<f64> {
    if temp_history.is_empty() {
        return vec![];
    }
    let mut accumulated = 0.0_f64;
    let mut strain_history = Vec::with_capacity(temp_history.len());
    for i in 0..temp_history.len() {
        if i == 0 {
            strain_history.push(accumulated);
            continue;
        }
        let dt = temp_history[i].0 - temp_history[i - 1].0;
        let temp = temp_history[i].1;
        accumulated += norton.creep_strain_rate(stress, temp) * dt.max(0.0);
        strain_history.push(accumulated);
    }
    strain_history
}
/// Predict rupture life from multiple methods for comparison.
#[allow(clippy::too_many_arguments)]
pub fn predict_rupture_life(
    temperature: f64,
    lm: &LarsonMillerParameter,
    lm_param_for_stress: f64,
    osd: &OrrSherbyDornParameter,
    osd_param_for_stress: f64,
    mg: &MonkmanGrant,
    steady_state_rate: f64,
) -> RuptureLifePrediction {
    RuptureLifePrediction {
        larson_miller: lm.rupture_time(temperature, lm_param_for_stress),
        orr_sherby_dorn: osd.rupture_time(temperature, osd_param_for_stress),
        monkman_grant: mg.rupture_life(steady_state_rate),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AndradeCreep;
    use crate::ChabocheKinematicHardening;
    use crate::CobleCreep;
    use crate::CoupledCreepDamage;
    use crate::CreepCurve;
    use crate::CreepDamage;
    use crate::CreepFatigueInteraction;
    use crate::DeformationMechanismMap;
    use crate::LarsonMillerStress;
    use crate::MansonHaferdParameter;
    use crate::MultiaxialCreep;
    use crate::NabarroHerringCreep;
    use crate::NortonBaileyTimeHardening;
    use crate::NortonCreep;
    use crate::PowerLawCreep;

    use crate::RuptureEnvelope;

    use crate::SherbyDornParameter;
    use crate::StressRelaxation;

    use crate::ThermalCreepIntegrator;
    use crate::ThermallyActivatedCreep;

    use crate::ViscoplasticModel;
    use crate::ZenerHollomon;
    #[test]
    fn test_norton_higher_stress_higher_rate() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let r1 = norton.creep_strain_rate(100e6, 800.0);
        let r2 = norton.creep_strain_rate(200e6, 800.0);
        assert!(r2 > r1, "Higher stress should give higher creep rate");
    }
    #[test]
    fn test_norton_higher_temperature_higher_rate() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let r_low = norton.creep_strain_rate(100e6, 700.0);
        let r_high = norton.creep_strain_rate(100e6, 1000.0);
        assert!(
            r_high > r_low,
            "Higher temperature should give higher creep rate"
        );
    }
    #[test]
    fn test_norton_effective_viscosity() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let eta = norton.effective_viscosity(100e6, 800.0);
        assert!(eta > 0.0, "viscosity should be positive");
        assert!(eta.is_finite(), "viscosity should be finite");
    }
    #[test]
    fn test_power_law_creep() {
        let plc = PowerLawCreep::new(1e-6, 100e6, 3.0);
        let r1 = plc.strain_rate(100e6);
        let r2 = plc.strain_rate(200e6);
        assert!((r2 / r1 - 8.0).abs() < 1e-6, "ratio = {}", r2 / r1);
    }
    #[test]
    fn test_power_law_strain_at_time() {
        let plc = PowerLawCreep::new(1e-6, 100e6, 3.0);
        let strain = plc.strain_at_time(100e6, 3600.0);
        let expected = 1e-6 * 3600.0;
        assert!((strain - expected).abs() < 1e-15);
    }
    #[test]
    fn test_nabarro_herring_rate_increases_with_temperature() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let r1 = nh.strain_rate(10e6, 800.0);
        let r2 = nh.strain_rate(10e6, 1200.0);
        assert!(r2 > r1, "NH rate should increase with temperature");
    }
    #[test]
    fn test_nabarro_herring_linear_in_stress() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let r1 = nh.strain_rate(10e6, 1000.0);
        let r2 = nh.strain_rate(20e6, 1000.0);
        assert!(
            (r2 / r1 - 2.0).abs() < 1e-6,
            "NH should be linear in stress"
        );
    }
    #[test]
    fn test_coble_creep_grain_size_dependence() {
        let c1 = CobleCreep::new(1e-5, 150_000.0, 5e-10, 1.2e-29, 10e-6);
        let c2 = CobleCreep::new(1e-5, 150_000.0, 5e-10, 1.2e-29, 5e-6);
        let r1 = c1.strain_rate(10e6, 1000.0);
        let r2 = c2.strain_rate(10e6, 1000.0);
        assert!(
            (r2 / r1 - 8.0).abs() < 0.1,
            "Coble d^-3: ratio = {}",
            r2 / r1
        );
    }
    #[test]
    fn test_dominant_diffusion_creep() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let coble = CobleCreep::new(1e-5, 150_000.0, 5e-10, 1.2e-29, 10e-6);
        let dom = dominant_diffusion_creep(&nh, &coble, 10e6, 1000.0);
        assert!(dom == "NH" || dom == "Coble");
    }
    #[test]
    fn test_creep_curve_monotonically_increasing() {
        let curve = CreepCurve::new(1e-4, 1e-6, 0.05);
        let times = [0.0, 100.0, 1000.0, 10_000.0, 100_000.0];
        let strains: Vec<f64> = times.iter().map(|&t| curve.strain_at_time(t)).collect();
        for w in strains.windows(2) {
            assert!(w[1] >= w[0]);
        }
    }
    #[test]
    fn test_creep_curve_strain_rate() {
        let curve = CreepCurve::new(1e-4, 1e-6, 0.05);
        let rate = curve.strain_rate_at_time(10000.0);
        assert!(rate > 0.0);
        assert!(rate.is_finite());
        let rate_large = curve.strain_rate_at_time(1e20);
        assert!((rate_large - 1e-6).abs() < 1e-10);
    }
    #[test]
    fn test_viscoplastic_zero_overstress_zero_rate() {
        let model = ViscoplasticModel::new(250e6, 1e10, 1e9);
        let rate = model.plastic_strain_rate(100e6, 0.0);
        assert!(rate.abs() < f64::EPSILON);
    }
    #[test]
    fn test_chaboche_back_stress_saturates() {
        let model = ChabocheKinematicHardening::new(50e3, 500.0);
        let mut x = 0.0_f64;
        let dep = 1e-4_f64;
        let dt = 1.0_f64;
        for _ in 0..10_000 {
            x = model.update_back_stress(x, dep, dt);
        }
        let x_sat = model.saturation_stress();
        assert!(
            (x - x_sat).abs() / x_sat < 0.01,
            "Back stress should saturate near C/gamma = {x_sat}, got {x}"
        );
    }
    #[test]
    fn test_thermal_creep_integrator_strain_increases() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let integrator = ThermalCreepIntegrator::new(norton);
        let history: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, 100e6)).collect();
        let strains = integrator.integrate(&history, 800.0, 3600.0);
        for w in strains.windows(2) {
            assert!(w[1] > w[0]);
        }
    }
    #[test]
    fn test_creep_damage_reaches_one() {
        let damage = CreepDamage::new(1e-6, 3.0, 2.0);
        let history = damage.integrate(100.0, 0.1, 1_000_000);
        let final_omega = *history.last().unwrap();
        assert!(
            final_omega >= 0.99,
            "Damage should reach ~1.0, got {final_omega}"
        );
    }
    #[test]
    fn test_creep_damage_effective_stress() {
        let damage = CreepDamage::new(1e-6, 3.0, 2.0);
        let sigma_eff = damage.effective_stress(100.0, 0.5);
        assert!((sigma_eff - 200.0).abs() < 1e-10, "sigma_eff = {sigma_eff}");
    }
    #[test]
    fn test_larson_miller_roundtrip() {
        let lm = LarsonMillerParameter::new(20.0);
        let t1 = 1000.0;
        let temp = 800.0;
        let p = lm.lm_parameter(temp, t1);
        let t_back = lm.rupture_time(temp, p);
        assert!(
            (t_back - t1).abs() / t1 < 1e-6,
            "t = {t_back}, expected {t1}"
        );
    }
    #[test]
    fn test_larson_miller_equivalent_time() {
        let lm = LarsonMillerParameter::new(20.0);
        let t_equiv = lm.equivalent_time(1000.0, 800.0, 900.0);
        assert!(t_equiv < 1000.0, "higher temp should give shorter life");
    }
    #[test]
    fn test_stress_relaxation_exponential() {
        let sr = StressRelaxation::new(200e9, 1e14);
        let tau = sr.relaxation_time();
        let sigma_0 = 100e6;
        let sigma_tau = sr.stress_at_time(sigma_0, tau);
        let expected = sigma_0 / std::f64::consts::E;
        assert!(
            (sigma_tau - expected).abs() < 1.0,
            "sigma(tau) = {sigma_tau}, expected {expected}"
        );
    }
    #[test]
    fn test_stress_relaxation_time_to_fraction() {
        let sr = StressRelaxation::new(200e9, 1e14);
        let tau = sr.relaxation_time();
        let t_half = sr.time_to_fraction(0.5);
        let expected = -tau * 0.5_f64.ln();
        assert!((t_half - expected).abs() < 1e-6 * expected);
    }
    #[test]
    fn test_stress_relaxation_history() {
        let sr = StressRelaxation::new(200e9, 1e14);
        let history = sr.stress_history(100e6, 1e6, 10);
        assert_eq!(history.len(), 11);
        for w in history.windows(2) {
            assert!(w[1].1 <= w[0].1, "stress should decrease");
        }
    }
    #[test]
    fn test_creep_stress_relaxation() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let history = creep_stress_relaxation(&norton, 100e6, 800.0, 200e9, 3600.0, 10);
        assert_eq!(history.len(), 11);
        for w in history.windows(2) {
            assert!(w[1].1 <= w[0].1, "stress should relax");
        }
    }
    #[test]
    fn test_creep_fatigue_zero_damage_safe() {
        let cfi = CreepFatigueInteraction::new(0.3, 1.0, 1.0);
        assert!(cfi.is_safe(0.0, 0.0), "Zero damage should be safe");
    }
    #[test]
    fn test_creep_fatigue_full_damage_unsafe() {
        let cfi = CreepFatigueInteraction::new(0.3, 1.0, 1.0);
        assert!(
            !cfi.is_safe(1.0, 1.0),
            "Full combined damage should be unsafe"
        );
    }
    #[test]
    fn test_creep_fatigue_interaction_damage_at_zero() {
        let cfi = CreepFatigueInteraction::new(0.3, 1.0, 1.0);
        let d = cfi.interaction_damage(0.0, 0.0);
        assert!(d.abs() < 1e-12, "Interaction damage at zero should be 0");
    }
    #[test]
    fn test_creep_fatigue_remaining_creep() {
        let cfi = CreepFatigueInteraction::new(0.5, 1.0, 1.0);
        let rem = cfi.remaining_creep(0.0);
        assert!(
            (rem - 0.5).abs() < 1e-10,
            "Remaining creep should be intercept, got {rem}"
        );
    }
    #[test]
    fn test_creep_fatigue_remaining_decreases_with_fatigue() {
        let cfi = CreepFatigueInteraction::new(0.5, 1.0, 1.0);
        let rem0 = cfi.remaining_creep(0.0);
        let rem1 = cfi.remaining_creep(0.5);
        assert!(
            rem1 < rem0,
            "More fatigue damage → less remaining creep allowance"
        );
    }
    #[test]
    fn test_coupled_creep_damage_increases_over_time() {
        let m = CoupledCreepDamage::new(1e-20, 5.0, 150_000.0, 1e-6, 3.0, 2.0);
        let (_eps, dam) = m.integrate(100e6, 800.0, 3600.0, 100);
        assert!(
            *dam.last().unwrap() > 0.0,
            "Damage should increase from zero"
        );
    }
    #[test]
    fn test_coupled_creep_damage_strain_positive() {
        let m = CoupledCreepDamage::new(1e-20, 5.0, 150_000.0, 1e-6, 3.0, 2.0);
        let (eps, _) = m.integrate(100e6, 800.0, 3600.0, 10);
        assert!(
            *eps.last().unwrap() > 0.0,
            "Creep strain should be positive"
        );
    }
    #[test]
    fn test_coupled_creep_damage_rate_higher_with_damage() {
        let m = CoupledCreepDamage::new(1e-20, 5.0, 150_000.0, 1e-6, 3.0, 2.0);
        let rate0 = m.creep_rate(100e6, 800.0, 0.0);
        let rate1 = m.creep_rate(100e6, 800.0, 0.5);
        assert!(rate1 > rate0, "Creep rate should be higher with damage");
    }
    #[test]
    fn test_multiaxial_von_mises_uniaxial() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let mc = MultiaxialCreep::new(norton);
        let s_eq = mc.von_mises_stress(100e6, 0.0, 0.0);
        assert!(
            (s_eq - 100e6).abs() < 1e-3,
            "Von Mises for uniaxial = σ1, got {s_eq}"
        );
    }
    #[test]
    fn test_multiaxial_hydrostatic_zero_von_mises() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let mc = MultiaxialCreep::new(norton);
        let s_eq = mc.von_mises_stress(100e6, 100e6, 100e6);
        assert!(
            s_eq.abs() < 1.0,
            "Von Mises for hydrostatic should be 0, got {s_eq}"
        );
    }
    #[test]
    fn test_multiaxial_volume_preserving() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let mc = MultiaxialCreep::new(norton);
        let vol_rate = mc.volumetric_creep_rate(200e6, 50e6, 0.0, 800.0);
        assert!(
            vol_rate.abs() < 1e-3,
            "Creep should be volume-preserving (div ≈ 0), got {vol_rate}"
        );
    }
    #[test]
    fn test_multiaxial_uniaxial_principal_rates() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let mc = MultiaxialCreep::new(norton);
        let rates = mc.principal_creep_rates(100e6, 0.0, 0.0, 800.0);
        assert!(rates[0] > 0.0, "Axial rate should be positive");
        assert!(
            rates[1] < 0.0,
            "Transverse rate should be negative (incompressible)"
        );
    }
    #[test]
    fn test_manson_haferd_roundtrip() {
        let mh = MansonHaferdParameter::new(3.0, 200.0);
        let t = 1000.0_f64;
        let temp = 800.0_f64;
        let p = mh.mh_parameter(t, temp);
        let t_back = mh.rupture_time(temp, p);
        assert!(
            (t_back - t).abs() / t < 1e-6,
            "MH roundtrip failed: {t_back} vs {t}"
        );
    }
    #[test]
    fn test_manson_haferd_higher_temp_shorter_life() {
        let mh = MansonHaferdParameter::new(3.0, 200.0);
        let p = mh.mh_parameter(100.0, 800.0);
        let t_high = mh.rupture_time(900.0, p);
        let t_low = mh.rupture_time(800.0, p);
        assert!(
            t_high < t_low,
            "Higher temperature should give shorter life: t_high={t_high}, t_low={t_low}"
        );
    }
    #[test]
    fn test_osd_roundtrip() {
        let osd = OrrSherbyDornParameter::new(250_000.0);
        let t = 500.0_f64;
        let temp = 850.0_f64;
        let p = osd.osd_parameter(t, temp);
        let t_back = osd.rupture_time(temp, p);
        assert!(
            (t_back - t).abs() / t < 1e-4,
            "OSD roundtrip: {t_back} vs {t}"
        );
    }
    #[test]
    fn test_osd_parameter_changes_with_temperature() {
        let osd = OrrSherbyDornParameter::new(250_000.0);
        let p1 = osd.osd_parameter(1000.0, 700.0);
        let p2 = osd.osd_parameter(1000.0, 900.0);
        assert!(p1 != p2, "OSD parameter should change with temperature");
    }
    #[test]
    fn test_variable_temp_creep_increases() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let history: Vec<(f64, f64)> = (0..10).map(|i| (i as f64 * 3600.0, 800.0)).collect();
        let strains = variable_temperature_creep(&norton, 100e6, &history);
        assert_eq!(strains.len(), 10);
        for w in strains.windows(2) {
            assert!(w[1] >= w[0], "Creep strain should not decrease");
        }
    }
    #[test]
    fn test_variable_temp_creep_empty_history() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let strains = variable_temperature_creep(&norton, 100e6, &[]);
        assert!(strains.is_empty());
    }
    #[test]
    fn test_variable_temp_higher_temp_more_creep() {
        let norton = NortonCreep::new(1e-20, 5.0, 150_000.0);
        let hist_low: Vec<(f64, f64)> = (0..5).map(|i| (i as f64 * 3600.0, 700.0)).collect();
        let hist_high: Vec<(f64, f64)> = (0..5).map(|i| (i as f64 * 3600.0, 900.0)).collect();
        let s_low = variable_temperature_creep(&norton, 100e6, &hist_low);
        let s_high = variable_temperature_creep(&norton, 100e6, &hist_high);
        assert!(
            s_high.last() > s_low.last(),
            "Higher temperature should produce more creep"
        );
    }
    #[test]
    fn test_andrade_strain_at_zero() {
        let m = AndradeCreep::new(1e-4, 1e-7);
        assert_eq!(m.strain(0.0), 0.0, "Andrade strain at t=0 should be 0");
    }
    #[test]
    fn test_andrade_strain_monotonically_increasing() {
        let m = AndradeCreep::new(1e-4, 1e-7);
        let s1 = m.strain(1000.0);
        let s2 = m.strain(10_000.0);
        let s3 = m.strain(100_000.0);
        assert!(
            s2 > s1 && s3 > s2,
            "Andrade strain should increase monotonically"
        );
    }
    #[test]
    fn test_andrade_strain_rate_decreasing_early() {
        let m = AndradeCreep::new(1e-4, 1e-10);
        let r1 = m.strain_rate(100.0);
        let r2 = m.strain_rate(10_000.0);
        assert!(r2 < r1, "Primary creep rate should decrease with time");
    }
    #[test]
    fn test_andrade_strain_rate_approaches_steady_state() {
        let k = 1e-6_f64;
        let m = AndradeCreep::new(1e-8, k);
        let rate = m.strain_rate(1e20);
        assert!(
            (rate - k).abs() / k < 0.01,
            "At large t, Andrade rate → k = {k}, got {rate}"
        );
    }
    #[test]
    fn test_andrade_transition_time() {
        let beta = 1e-3;
        let k = 1e-5;
        let m = AndradeCreep::new(beta, k);
        let t_trans = m.transition_time();
        let expected = (beta / k).powf(1.5);
        assert!((t_trans - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_sherby_dorn_constant_temp() {
        let sd = SherbyDornParameter::new(150_000.0);
        let p = sd.parameter_constant_temp(1000.0, 800.0);
        assert!(p > 0.0, "Sherby-Dorn parameter should be positive");
    }
    #[test]
    fn test_sherby_dorn_higher_temp_larger_param() {
        let sd = SherbyDornParameter::new(150_000.0);
        let p_low = sd.parameter_constant_temp(1000.0, 700.0);
        let p_high = sd.parameter_constant_temp(1000.0, 1000.0);
        assert!(
            p_high > p_low,
            "Higher temperature → larger Sherby-Dorn parameter"
        );
    }
    #[test]
    fn test_sherby_dorn_equivalent_time_roundtrip() {
        let sd = SherbyDornParameter::new(150_000.0);
        let t_orig = 500.0_f64;
        let temp = 800.0;
        let p = sd.parameter_constant_temp(t_orig, temp);
        let t_back = sd.equivalent_time(p, temp);
        assert!(
            (t_back - t_orig).abs() / t_orig < 1e-6,
            "Sherby-Dorn roundtrip: {t_back} vs {t_orig}"
        );
    }
    #[test]
    fn test_sherby_dorn_variable_temp() {
        let sd = SherbyDornParameter::new(150_000.0);
        let segments = vec![(100.0_f64, 1000.0_f64), (1000.0, 700.0)];
        let p = sd.parameter_variable_temp(&segments);
        assert!(p > 0.0, "Variable temp P_SD should be positive");
    }
    #[test]
    fn test_monkman_grant_rupture_life() {
        let mg = MonkmanGrant::new(0.05, 1.0);
        let rate = 1e-5;
        let t_f = mg.rupture_life(rate);
        let expected = 0.05 / rate;
        assert!((t_f - expected).abs() / expected < 1e-10, "t_f = {t_f}");
    }
    #[test]
    fn test_monkman_grant_higher_rate_shorter_life() {
        let mg = MonkmanGrant::new(0.05, 1.0);
        let t_f1 = mg.rupture_life(1e-5);
        let t_f2 = mg.rupture_life(1e-4);
        assert!(t_f2 < t_f1, "Higher creep rate → shorter life");
    }
    #[test]
    fn test_monkman_grant_roundtrip() {
        let mg = MonkmanGrant::new(0.05, 1.2);
        let rate_orig = 1e-5;
        let t_f = mg.rupture_life(rate_orig);
        let rate_back = mg.steady_state_rate_from_life(t_f);
        assert!(
            (rate_back - rate_orig).abs() / rate_orig < 1e-9,
            "MG roundtrip: {rate_back} vs {rate_orig}"
        );
    }
    #[test]
    fn test_monkman_grant_consistency() {
        let mg = MonkmanGrant::new(0.05, 1.2);
        assert!(
            mg.is_consistent(1e-5, 1e-9),
            "MG model should be internally consistent"
        );
    }
    #[test]
    fn test_lm_stress_parameter_decreases_with_stress() {
        let lms = LarsonMillerStress::new(30_000.0, 5000.0, 100.0e6, 20.0);
        let p_low = lms.lm_parameter(50.0e6);
        let p_high = lms.lm_parameter(200.0e6);
        assert!(p_high < p_low, "Higher stress should reduce LM parameter");
    }
    #[test]
    fn test_lm_stress_rupture_time_positive() {
        let lms = LarsonMillerStress::new(30_000.0, 5000.0, 100.0e6, 20.0);
        let t_r = lms.rupture_time(100.0e6, 800.0);
        assert!(t_r > 0.0, "Rupture time should be positive");
        assert!(t_r.is_finite(), "Rupture time should be finite");
    }
    #[test]
    fn test_lm_stress_max_stress_decreases_with_life() {
        let lms = LarsonMillerStress::new(30_000.0, 5000.0, 100.0e6, 20.0);
        let sigma1 = lms.max_stress(800.0, 100.0);
        let sigma2 = lms.max_stress(800.0, 10_000.0);
        assert!(sigma2 < sigma1, "Longer life → lower allowable stress");
    }
    #[test]
    fn test_zener_hollomon_z_positive() {
        let zh = ZenerHollomon::new(150_000.0);
        let z = zh.z_parameter(1.0, 1000.0);
        assert!(z > 0.0, "Z parameter should be positive");
    }
    #[test]
    fn test_zener_hollomon_higher_rate_higher_z() {
        let zh = ZenerHollomon::new(150_000.0);
        let z1 = zh.z_parameter(1.0, 1000.0);
        let z2 = zh.z_parameter(10.0, 1000.0);
        assert!(
            (z2 / z1 - 10.0).abs() < 1e-6,
            "Z should scale with strain rate"
        );
    }
    #[test]
    fn test_zener_hollomon_equivalent_rate_roundtrip() {
        let zh = ZenerHollomon::new(150_000.0);
        let eps_dot_orig = 2.0;
        let temp = 900.0;
        let z = zh.z_parameter(eps_dot_orig, temp);
        let eps_dot_back = zh.equivalent_strain_rate(z, temp);
        assert!(
            (eps_dot_back - eps_dot_orig).abs() / eps_dot_orig < 1e-9,
            "ZH roundtrip: {eps_dot_back} vs {eps_dot_orig}"
        );
    }
    #[test]
    fn test_zener_hollomon_temperature_for_z() {
        let zh = ZenerHollomon::new(150_000.0);
        let eps_dot = 1.0;
        let temp_target = 900.0;
        let z = zh.z_parameter(eps_dot, temp_target);
        let temp_back = zh.temperature_for_z(z, eps_dot);
        assert!(
            (temp_back - temp_target).abs() / temp_target < 1e-6,
            "ZH temperature roundtrip: {temp_back} vs {temp_target}"
        );
    }
    #[test]
    fn test_norton_bailey_strain_at_zero_time() {
        let nb = NortonBaileyTimeHardening::new(1e-15, 5.0, 0.5);
        assert_eq!(nb.strain(100e6, 0.0), 0.0, "Strain at t=0 should be zero");
    }
    #[test]
    fn test_norton_bailey_strain_increases_with_time() {
        let nb = NortonBaileyTimeHardening::new(1e-15, 5.0, 0.5);
        let s1 = nb.strain(100e6, 100.0);
        let s2 = nb.strain(100e6, 10_000.0);
        assert!(s2 > s1, "Strain should increase with time");
    }
    #[test]
    fn test_norton_bailey_strain_increases_with_stress() {
        let nb = NortonBaileyTimeHardening::new(1e-15, 5.0, 0.5);
        let s1 = nb.strain(100e6, 1000.0);
        let s2 = nb.strain(200e6, 1000.0);
        assert!((s2 / s1 - 32.0).abs() < 1e-6, "Stress scaling: {}", s2 / s1);
    }
    #[test]
    fn test_norton_bailey_equivalent_time() {
        let nb = NortonBaileyTimeHardening::new(1e-15, 5.0, 0.5);
        let sigma = 100e6;
        let t_orig = 1000.0;
        let eps = nb.strain(sigma, t_orig);
        let t_eq = nb.equivalent_time(eps, sigma);
        assert!(
            (t_eq - t_orig).abs() / t_orig < 1e-9,
            "NB equivalent time roundtrip: {t_eq} vs {t_orig}"
        );
    }
    #[test]
    fn test_predict_rupture_life_all_positive() {
        let lm = LarsonMillerParameter::new(20.0);
        let osd = OrrSherbyDornParameter::new(250_000.0);
        let mg = MonkmanGrant::new(0.05, 1.0);
        let temp = 800.0;
        let t_lm_ref = 1000.0;
        let lm_param = lm.lm_parameter(temp, t_lm_ref);
        let t_osd_ref = 1000.0;
        let osd_param = osd.osd_parameter(t_osd_ref, temp);
        let result = predict_rupture_life(temp, &lm, lm_param, &osd, osd_param, &mg, 1e-5);
        assert!(result.larson_miller > 0.0, "LM life should be positive");
        assert!(result.orr_sherby_dorn > 0.0, "OSD life should be positive");
        assert!(result.monkman_grant > 0.0, "MG life should be positive");
    }
    #[test]
    fn test_thermally_activated_rate_positive() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let rate = tac.strain_rate(100.0e6, 700.0);
        assert!(rate > 0.0, "Thermally activated rate should be positive");
        assert!(rate.is_finite(), "Rate should be finite");
    }
    #[test]
    fn test_thermally_activated_rate_increases_with_temperature() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let r1 = tac.strain_rate(100.0e6, 500.0);
        let r2 = tac.strain_rate(100.0e6, 800.0);
        assert!(r2 > r1, "Higher T → higher thermally activated rate");
    }
    #[test]
    fn test_thermally_activated_rate_increases_with_stress() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let r1 = tac.strain_rate(50.0e6, 700.0);
        let r2 = tac.strain_rate(100.0e6, 700.0);
        assert!(r2 > r1, "Higher stress → higher thermally activated rate");
    }
    #[test]
    fn test_thermally_activated_activation_stress_roundtrip() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let sigma = 100.0e6;
        let temp = 700.0;
        let rate = tac.strain_rate(sigma, temp);
        let sigma_back = tac.activation_stress(rate, temp);
        assert!(
            (sigma_back - sigma).abs() / sigma < 1e-6,
            "Activation stress roundtrip: {sigma_back} vs {sigma}"
        );
    }
    #[test]
    fn test_thermally_activated_apparent_activation_energy_decreases_with_stress() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let dh1 = tac.apparent_activation_energy(50.0e6);
        let dh2 = tac.apparent_activation_energy(100.0e6);
        assert!(
            dh2 < dh1,
            "Apparent activation energy decreases with stress"
        );
    }
    #[test]
    fn test_thermally_activated_decreases_sharply_at_low_temperature() {
        let tac = ThermallyActivatedCreep::new(1e12, 2.0e-19, 5.0e-29);
        let r_high = tac.strain_rate(100.0e6, 900.0);
        let r_low = tac.strain_rate(100.0e6, 300.0);
        assert!(
            r_low < r_high * 1e-5,
            "Rate at 300K should be << rate at 900K"
        );
    }
    #[test]
    fn test_dmm_high_stress_gives_dislocation_glide() {
        let dmm = DeformationMechanismMap::nickel();
        let sigma = dmm.shear_modulus * (dmm.peierls_stress_norm + 0.01);
        let mech = dmm.dominant_mechanism(sigma, 800.0);
        assert_eq!(mech, "dislocation_glide", "High stress → dislocation glide");
    }
    #[test]
    fn test_dmm_low_temperature_gives_elastic() {
        let dmm = DeformationMechanismMap::nickel();
        let mech = dmm.dominant_mechanism(1.0e6, 100.0);
        assert_eq!(mech, "elastic", "Low T → elastic");
    }
    #[test]
    fn test_dmm_high_normalised_stress_power_law() {
        let dmm = DeformationMechanismMap::nickel();
        let sigma = dmm.shear_modulus * (dmm.power_law_diff_boundary + 1e-3);
        let mech = dmm.dominant_mechanism(sigma, 800.0);
        assert_eq!(mech, "power_law_creep", "Moderate σ → power law creep");
    }
    #[test]
    fn test_dmm_high_temperature_superplastic() {
        let dmm = DeformationMechanismMap::nickel();
        let sigma_low = dmm.shear_modulus * (dmm.power_law_diff_boundary * 0.1);
        let mech = dmm.dominant_mechanism(sigma_low, 1500.0);
        assert_eq!(mech, "superplastic", "High T + low σ → superplastic");
    }
    #[test]
    fn test_dmm_combined_strain_rate_positive() {
        let dmm = DeformationMechanismMap::nickel();
        let rate = dmm.combined_strain_rate(50.0e6, 900.0, 1e6, 200_000.0, 1e4, 150_000.0);
        assert!(rate > 0.0, "Combined strain rate should be positive");
        assert!(rate.is_finite(), "Rate should be finite");
    }
    #[test]
    fn test_dmm_combined_rate_increases_with_stress() {
        let dmm = DeformationMechanismMap::nickel();
        let r1 = dmm.combined_strain_rate(50.0e6, 900.0, 1e6, 200_000.0, 1e4, 150_000.0);
        let r2 = dmm.combined_strain_rate(100.0e6, 900.0, 1e6, 200_000.0, 1e4, 150_000.0);
        assert!(r2 > r1, "Higher stress → higher combined strain rate");
    }
    #[test]
    fn test_rupture_envelope_life_positive() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let tr = envelope.rupture_life(100.0e6, 800.0);
        assert!(tr > 0.0, "Rupture life should be positive, got {tr}");
        assert!(tr.is_finite(), "Rupture life should be finite");
    }
    #[test]
    fn test_rupture_envelope_higher_stress_shorter_life() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let tr1 = envelope.rupture_life(100.0e6, 800.0);
        let tr2 = envelope.rupture_life(200.0e6, 800.0);
        assert!(tr2 < tr1, "Higher stress → shorter rupture life");
    }
    #[test]
    fn test_rupture_envelope_higher_temperature_shorter_life() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let tr1 = envelope.rupture_life(100.0e6, 700.0);
        let tr2 = envelope.rupture_life(100.0e6, 900.0);
        assert!(tr2 < tr1, "Higher temperature → shorter rupture life");
    }
    #[test]
    fn test_rupture_envelope_allowable_stress_decreases_with_temperature() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let s1 = envelope.allowable_stress(700.0, 1000.0 * 3600.0);
        let s2 = envelope.allowable_stress(900.0, 1000.0 * 3600.0);
        assert!(s2 < s1, "Higher T → lower allowable stress");
    }
    #[test]
    fn test_rupture_envelope_allowable_stress_decreases_with_life() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let s_short = envelope.allowable_stress(800.0, 100.0 * 3600.0);
        let s_long = envelope.allowable_stress(800.0, 10_000.0 * 3600.0);
        assert!(
            s_long < s_short,
            "Longer design life → lower allowable stress"
        );
    }
    #[test]
    fn test_rupture_envelope_lm_param_decreases_with_stress() {
        let envelope = RuptureEnvelope::new(20.0, [30_000.0, -2000.0, 0.0]);
        let p1 = envelope.lm_parameter_from_stress(100.0e6);
        let p2 = envelope.lm_parameter_from_stress(200.0e6);
        assert!(
            p2 < p1,
            "Higher stress → lower LM parameter (for negative slope)"
        );
    }
}
#[cfg(test)]
mod tests_creep_expansion {

    use crate::CobleCreep;

    use crate::NabarroHerringCreep;

    use crate::NortonCreep;

    use crate::PrimaryCreepModel;

    use crate::SecondaryCreepModel;

    use crate::TertiaryCreepModel;

    use crate::ThreeStageCreep;

    #[test]
    fn test_nh_effective_diffusivity_positive() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let d = nh.effective_diffusivity(1000.0);
        assert!(d > 0.0, "effective diffusivity > 0: {d}");
    }
    #[test]
    fn test_nh_diffusivity_increases_with_temperature() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let d1 = nh.effective_diffusivity(800.0);
        let d2 = nh.effective_diffusivity(1000.0);
        assert!(d2 > d1, "diffusivity increases with temperature");
    }
    #[test]
    fn test_nh_compliance_positive() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let c = nh.compliance(1000.0);
        assert!(c > 0.0, "compliance > 0: {c}");
    }
    #[test]
    fn test_nh_activation_energy_getter() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        assert!(
            (nh.activation_energy() - 200_000.0).abs() < 1.0,
            "activation energy: {}",
            nh.activation_energy()
        );
    }
    #[test]
    fn test_coble_gb_diffusivity_positive() {
        let c = CobleCreep::new(1e-5, 150_000.0, 5e-10, 1.2e-29, 10e-6);
        let d = c.gb_diffusivity(900.0);
        assert!(d > 0.0, "GB diffusivity > 0: {d}");
    }
    #[test]
    fn test_coble_fraction_in_unit_interval() {
        let nh = NabarroHerringCreep::new(1e-5, 200_000.0, 1.2e-29, 10e-6);
        let coble = CobleCreep::new(1e-5, 150_000.0, 5e-10, 1.2e-29, 10e-6);
        let frac = coble.coble_fraction(&nh, 10e6, 1000.0);
        assert!(
            (0.0..=1.0).contains(&frac),
            "Coble fraction in [0,1]: {frac}"
        );
    }
    #[test]
    fn test_primary_creep_strain_positive() {
        let pm = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let eps = pm.strain(100e6, 1000.0);
        assert!(eps > 0.0, "primary strain > 0: {eps}");
    }
    #[test]
    fn test_primary_creep_strain_increases_with_time() {
        let pm = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let e1 = pm.strain(100e6, 100.0);
        let e2 = pm.strain(100e6, 1000.0);
        assert!(e2 > e1, "primary strain increases with time");
    }
    #[test]
    fn test_primary_creep_rate_decreasing_with_time() {
        let pm = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let r1 = pm.strain_rate(100e6, 100.0);
        let r2 = pm.strain_rate(100e6, 1000.0);
        assert!(
            r2 < r1,
            "primary creep rate decreases with time (hardening)"
        );
    }
    #[test]
    fn test_primary_creep_transition_time_positive() {
        let pm = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let norton = NortonCreep::new(1e-24, 4.0, 200_000.0);
        let t_tr = pm.transition_time(&norton, 100e6, 1000.0);
        assert!(
            t_tr > 0.0 && t_tr.is_finite(),
            "transition time > 0: {t_tr}"
        );
    }
    #[test]
    fn test_secondary_creep_rate_positive() {
        let sec = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let rate = sec.strain_rate(100e6, 1000.0);
        assert!(rate > 0.0, "secondary creep rate > 0: {rate}");
    }
    #[test]
    fn test_secondary_creep_accumulated_strain() {
        let sec = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let rate = sec.strain_rate(100e6, 1000.0);
        let eps = sec.accumulated_strain(100e6, 1000.0, 1000.0);
        assert!(
            (eps - rate * 1000.0).abs() / (rate * 1000.0) < 1e-10,
            "acc. strain = rate*dt"
        );
    }
    #[test]
    fn test_secondary_creep_viscosity_positive() {
        let sec = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let eta = sec.viscosity(100e6, 1000.0);
        assert!(eta > 0.0, "viscosity > 0: {eta}");
    }
    #[test]
    fn test_tertiary_initial_no_rupture() {
        let tert = TertiaryCreepModel::new(1e-24, 4.0, 1e-15, 4.0, 4.0);
        assert!(!tert.is_ruptured(), "initial damage = 0 → not ruptured");
    }
    #[test]
    fn test_tertiary_strain_rate_increases_with_damage() {
        let mut tert = TertiaryCreepModel::new(1e-24, 4.0, 1e-8, 4.0, 4.0);
        let rate0 = tert.strain_rate(100e6);
        tert.step(100e6, 1e3);
        let rate1 = tert.strain_rate(100e6);
        assert!(rate1 >= rate0, "strain rate should increase with damage");
    }
    #[test]
    fn test_tertiary_rupture_time_positive() {
        let tert = TertiaryCreepModel::new(1e-24, 4.0, 1e-8, 4.0, 4.0);
        let tr = tert.rupture_time(100e6);
        assert!(tr > 0.0 && tr.is_finite(), "rupture time > 0: {tr}");
    }
    #[test]
    fn test_tertiary_damage_rate_positive() {
        let tert = TertiaryCreepModel::new(1e-24, 4.0, 1e-8, 4.0, 4.0);
        let dw = tert.damage_rate(100e6);
        assert!(dw > 0.0, "damage rate > 0: {dw}");
    }
    #[test]
    fn test_tertiary_reset_clears_damage() {
        let mut tert = TertiaryCreepModel::new(1e-24, 4.0, 1e-8, 4.0, 4.0);
        tert.step(100e6, 1e4);
        tert.reset();
        assert_eq!(tert.omega, 0.0, "reset should clear omega");
    }
    #[test]
    fn test_three_stage_initial_stage_primary() {
        let primary = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let secondary = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let tertiary = TertiaryCreepModel::new(1e-24, 4.0, 1e-15, 4.0, 4.0);
        let model = ThreeStageCreep::new(primary, secondary, tertiary, 0.01);
        assert_eq!(
            model.current_stage(),
            "primary",
            "initial stage should be primary"
        );
    }
    #[test]
    fn test_three_stage_step_increments_strain() {
        let primary = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let secondary = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let tertiary = TertiaryCreepModel::new(1e-24, 4.0, 1e-15, 4.0, 4.0);
        let mut model = ThreeStageCreep::new(primary, secondary, tertiary, 0.01);
        let d_eps = model.step(100e6, 1000.0, 100.0);
        assert!(d_eps >= 0.0, "incremental strain >= 0: {d_eps}");
        assert!(model.strain > 0.0, "total strain should increase");
    }
    #[test]
    fn test_three_stage_not_ruptured_initially() {
        let primary = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let secondary = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let tertiary = TertiaryCreepModel::new(1e-24, 4.0, 1e-15, 4.0, 4.0);
        let model = ThreeStageCreep::new(primary, secondary, tertiary, 0.01);
        assert!(!model.is_ruptured(), "should not be ruptured initially");
    }
    #[test]
    fn test_three_stage_time_advances() {
        let primary = PrimaryCreepModel::new(1e-12, 4.0, 0.33);
        let secondary = SecondaryCreepModel::new(1e-24, 4.0, 200_000.0);
        let tertiary = TertiaryCreepModel::new(1e-24, 4.0, 1e-15, 4.0, 4.0);
        let mut model = ThreeStageCreep::new(primary, secondary, tertiary, 0.01);
        model.step(100e6, 1000.0, 50.0);
        model.step(100e6, 1000.0, 50.0);
        assert!(
            (model.time - 100.0).abs() < 1e-10,
            "time should be 100: {}",
            model.time
        );
    }
}
