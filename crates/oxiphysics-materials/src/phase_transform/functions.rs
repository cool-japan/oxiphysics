//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::E;

// Types are re-exported via mod.rs; only import in test module where needed.

/// Gas constant R (J / mol K).
pub(super) const R_GAS: f64 = 8.314_462_618;
/// Avrami (JMAK) isothermal transformation fraction.
///
/// X = 1 - exp(-k * t^n)
///
/// # Arguments
/// * `time` - Time (s).
/// * `n`    - Avrami exponent.
/// * `k`    - Rate constant (s^-n).
pub fn avrami_transformation(time: f64, n: f64, k: f64) -> f64 {
    if time <= 0.0 {
        return 0.0;
    }
    1.0 - E.powf(-k * time.powf(n))
}
/// Estimate transformation time on a CCT (continuous cooling transformation)
/// diagram via linear interpolation in log-time space.
///
/// # Arguments
/// * `t`       - Temperature of interest (K).
/// * `t_start` - Start temperature of the CCT nose (K).
/// * `t_end`   - End temperature of the CCT nose (K).
/// * `time_0`  - Log-interpolation time at T_start (s).
/// * `time_f`  - Log-interpolation time at T_end (s).
/// * `t_ref`   - Reference temperature for interpolation (K).
pub fn cct_diagram_transformation_time(
    t: f64,
    t_start: f64,
    t_end: f64,
    time_0: f64,
    time_f: f64,
    _t_ref: f64,
) -> f64 {
    if (t_end - t_start).abs() < 1e-15 {
        return time_0;
    }
    let frac = (t - t_start) / (t_end - t_start);
    let frac = frac.clamp(0.0, 1.0);
    let log_t0 = time_0.max(1e-30).ln();
    let log_tf = time_f.max(1e-30).ln();
    (log_t0 + frac * (log_tf - log_t0)).exp()
}
/// Landau free energy: F(φ) = a/2 * φ^2 - b/4 * φ^4.
///
/// This is a quartic double-well potential when b > 0 and a < 0.
///
/// # Arguments
/// * `phi` - Order parameter φ.
/// * `a`   - Coefficient of the φ^2 term (negative for broken symmetry).
/// * `b`   - Coefficient of the φ^4 term (positive for stability).
pub fn landau_free_energy(phi: f64, a: f64, b: f64) -> f64 {
    a / 2.0 * phi * phi - b / 4.0 * phi.powi(4)
}
/// Multi-component lever rule for ternary alloy A-B-C.
///
/// Given overall composition (x_b, x_c), and end-member compositions of
/// two co-existing phases α and β, returns the fraction of the β phase.
pub fn ternary_lever_rule(
    x_b_overall: f64,
    x_c_overall: f64,
    x_b_alpha: f64,
    x_c_alpha: f64,
    x_b_beta: f64,
    x_c_beta: f64,
) -> f64 {
    let db = x_b_beta - x_b_alpha;
    let dc = x_c_beta - x_c_alpha;
    let len2 = db * db + dc * dc;
    if len2 < 1e-60 {
        return 0.0;
    }
    let f = ((x_b_overall - x_b_alpha) * db + (x_c_overall - x_c_alpha) * dc) / len2;
    f.clamp(0.0, 1.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AusteniteMartensiteKinetics;
    use crate::BainiteTransformation;
    use crate::BrinsonModel;
    use crate::CahnHilliardField;
    use crate::CalphadModel;
    use crate::CctDiagram;
    use crate::ClausiusClapeyron;

    use crate::JmakExtended;
    use crate::JmatProModel;
    use crate::JohnsonMehlAvramiKolmogorov;
    use crate::KoistinenMarburger;
    use crate::LandauFreeEnergy;
    use crate::LatentHeatModel;
    use crate::LeverRule;
    use crate::MartensiteTransformation;
    use crate::PearliteTransformation;
    use crate::PhaseDependentProperties;
    use crate::PhaseFieldModel2D;
    use crate::PhaseFieldOrder;
    use crate::PhaseFractionEvolution;
    use crate::PhaseTransformType;

    use crate::SpinocalDecomposition;

    use crate::TransformationPlasticity;
    use crate::TripEffect;

    use crate::TttDiagram;

    use crate::VolumeFractionTracker;
    #[test]
    fn test_landau_order_param_zero_above_tc() {
        let landau = LandauFreeEnergy::new(1.0, 2.0, 500.0);
        let eta = landau.equilibrium_order_parameter(600.0);
        assert!(eta.abs() < 1.0e-12, "Expected eta=0 above Tc, got {eta}");
    }
    #[test]
    fn test_landau_order_param_positive_below_tc() {
        let landau = LandauFreeEnergy::new(1.0, 2.0, 500.0);
        let eta = landau.equilibrium_order_parameter(400.0);
        assert!(eta > 0.0, "Expected eta>0 below Tc, got {eta}");
        let expected = (50.0_f64).sqrt();
        assert!(
            (eta - expected).abs() < 1.0e-10,
            "Expected eta={expected}, got {eta}"
        );
    }
    #[test]
    fn test_martensite_fraction_increases_below_ms() {
        let mt = MartensiteTransformation::new(500.0, 300.0, 550.0, 700.0, 0.011, 0.02);
        let f1 = mt.martensite_fraction_koistinen_marburger(490.0);
        let f2 = mt.martensite_fraction_koistinen_marburger(400.0);
        assert!(f2 > f1, "Fraction should increase as T decreases below Ms");
    }
    #[test]
    fn test_martensite_fraction_at_mf() {
        let mt = MartensiteTransformation::new(500.0, 200.0, 550.0, 700.0, 0.05, 0.02);
        let f = mt.martensite_fraction_koistinen_marburger(200.0);
        assert!(f > 0.99, "Fraction at Mf should be close to 1, got {f}");
    }
    #[test]
    fn test_jmak_fraction_at_zero() {
        let jmak = JohnsonMehlAvramiKolmogorov::new(0.01, 2.0);
        assert_eq!(jmak.fraction_transformed(0.0), 0.0);
    }
    #[test]
    fn test_jmak_fraction_increases_with_time() {
        let jmak = JohnsonMehlAvramiKolmogorov::new(0.01, 2.0);
        let f1 = jmak.fraction_transformed(1.0);
        let f2 = jmak.fraction_transformed(5.0);
        let f3 = jmak.fraction_transformed(10.0);
        assert!(f1 < f2, "JMAK fraction must increase with time");
        assert!(f2 < f3, "JMAK fraction must increase with time");
    }
    #[test]
    fn test_calphad_mixing_enthalpy_pure_components() {
        let model = CalphadModel::new(-5000.0, -8000.0, 2000.0, 500.0);
        let h0 = model.mixing_enthalpy(0.0);
        let h1 = model.mixing_enthalpy(1.0);
        assert!(h0.abs() < 1.0e-12, "H_mix at x_B=0 should be 0, got {h0}");
        assert!(h1.abs() < 1.0e-12, "H_mix at x_B=1 should be 0, got {h1}");
    }
    #[test]
    fn test_phase_field_solid_fraction_increases_after_step() {
        let nx = 32;
        let ny = 32;
        let mut pf = PhaseFieldModel2D::new(nx, ny, 2.0, 1.0, 100.0);
        for t in pf.temperature.iter_mut() {
            *t = 1000.0;
        }
        pf.initialize_seed(16.0, 16.0, 4.0);
        let f0 = pf.solid_fraction();
        for _ in 0..10 {
            pf.step(1.0e-4, 1.0);
        }
        let f1 = pf.solid_fraction();
        assert!(
            f1 >= f0,
            "Solid fraction should not decrease after solidification steps: f0={f0}, f1={f1}"
        );
    }
    #[test]
    fn test_km_fraction_zero_at_ms() {
        let km = KoistinenMarburger::new(500.0, 0.011);
        assert!((km.fraction(500.0)).abs() < 1e-12);
    }
    #[test]
    fn test_km_fraction_increases_below_ms() {
        let km = KoistinenMarburger::new(500.0, 0.011);
        let f1 = km.fraction(490.0);
        let f2 = km.fraction(450.0);
        assert!(f2 > f1, "fraction increases as T decreases");
    }
    #[test]
    fn test_km_temp_for_fraction_roundtrip() {
        let km = KoistinenMarburger::new(500.0, 0.011);
        let t = km.temperature_for_fraction(0.5);
        let f = km.fraction(t);
        assert!((f - 0.5).abs() < 1e-9, "roundtrip: expected 0.5 got {f}");
    }
    #[test]
    fn test_austenite_martensite_kinetics_cooling() {
        let mut k = AusteniteMartensiteKinetics::new(500.0, 0.011, 550.0, 700.0, 600.0);
        k.step(450.0);
        assert!(k.f_martensite > 0.0, "martensite should form on cooling");
        assert!(k.f_austenite() < 1.0);
    }
    #[test]
    fn test_austenite_martensite_kinetics_heating_dissolves() {
        let mut k = AusteniteMartensiteKinetics::new(500.0, 0.011, 550.0, 700.0, 600.0);
        k.step(400.0);
        let fm_before = k.f_martensite;
        k.step(700.0);
        assert!(
            k.f_martensite < fm_before,
            "martensite should dissolve on heating"
        );
        assert!(
            k.f_martensite.abs() < 1e-9,
            "should be fully austenitic above Af"
        );
    }
    #[test]
    fn test_trip_strain_increment_positive() {
        let mut trip = TripEffect::new(1e-5);
        let de = trip.trip_strain_increment(200e6, 0.3, 0.1);
        assert!(de > 0.0);
        assert!((trip.accumulated_strain - de).abs() < 1e-15);
    }
    #[test]
    fn test_trip_no_increment_for_negative_df() {
        let mut trip = TripEffect::new(1e-5);
        let de = trip.trip_strain_increment(200e6, 0.3, -0.1);
        assert_eq!(de, 0.0);
    }
    #[test]
    fn test_trip_reset() {
        let mut trip = TripEffect::new(1e-5);
        trip.trip_strain_increment(200e6, 0.3, 0.1);
        trip.reset();
        assert_eq!(trip.accumulated_strain, 0.0);
    }
    #[test]
    fn test_brinson_effective_modulus_austenite() {
        let b = BrinsonModel::new(70e9, 30e9, 300.0, 250.0, 340.0, 380.0, 0.05);
        assert!((b.effective_modulus() - 70e9).abs() < 1.0);
    }
    #[test]
    fn test_brinson_update_temperature_below_mf() {
        let mut b = BrinsonModel::new(70e9, 30e9, 300.0, 200.0, 340.0, 380.0, 0.05);
        b.update_temperature(150.0);
        assert!(
            (b.xi_t - 1.0).abs() < 1e-9,
            "should be fully martensitic below Mf"
        );
    }
    #[test]
    fn test_brinson_update_temperature_above_af() {
        let mut b = BrinsonModel::new(70e9, 30e9, 300.0, 200.0, 340.0, 380.0, 0.05);
        b.update_temperature(150.0);
        b.update_temperature(400.0);
        assert!(b.xi_t.abs() < 1e-9, "should be fully austenitic above Af");
    }
    #[test]
    fn test_transformation_plasticity_increment() {
        let mut tp = TransformationPlasticity::new(0.01, 300e6);
        let de = tp.increment(100e6, 0.1);
        assert!(de > 0.0);
        assert!((tp.total_strain() - de).abs() < 1e-15);
    }
    #[test]
    fn test_transformation_plasticity_no_dissolution_increment() {
        let mut tp = TransformationPlasticity::new(0.01, 300e6);
        let de = tp.increment(100e6, -0.01);
        assert_eq!(de, 0.0);
    }
    #[test]
    fn test_phase_fraction_evolution_record_and_query() {
        let mut evo = PhaseFractionEvolution::new(3);
        evo.record(0.0, &[1.0, 0.0, 0.0]);
        evo.record(1.0, &[0.5, 0.3, 0.2]);
        assert!((evo.current(0) - 0.5).abs() < 1e-12);
        assert!((evo.current(1) - 0.3).abs() < 1e-12);
        assert!((evo.current(2) - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_phase_fraction_evolution_rate() {
        let mut evo = PhaseFractionEvolution::new(2);
        evo.record(0.0, &[1.0, 0.0]);
        evo.record(2.0, &[0.6, 0.4]);
        let rate_0 = evo.rate(0);
        let rate_1 = evo.rate(1);
        assert!((rate_0 + 0.2).abs() < 1e-9, "phase 0 rate = {rate_0}");
        assert!((rate_1 - 0.2).abs() < 1e-9, "phase 1 rate = {rate_1}");
    }
    #[test]
    fn test_phase_fraction_total() {
        let mut evo = PhaseFractionEvolution::new(3);
        evo.record(0.0, &[0.5, 0.3, 0.2]);
        assert!((evo.total_fraction() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_koistinen_marburger_at_ms_gives_zero() {
        let xm = JmatProModel::martensite_fraction_koistinen_marburger(500.0, 500.0);
        assert!(xm.abs() < 1e-12, "X_m at T=Ms should be 0, got {xm}");
    }
    #[test]
    fn test_koistinen_marburger_increases_below_ms() {
        let xm1 = JmatProModel::martensite_fraction_koistinen_marburger(490.0, 500.0);
        let xm2 = JmatProModel::martensite_fraction_koistinen_marburger(450.0, 500.0);
        assert!(
            xm2 > xm1,
            "Martensite fraction should increase as T decreases below Ms"
        );
    }
    #[test]
    fn test_avrami_at_t_zero_gives_zero() {
        let x = avrami_transformation(0.0, 2.0, 0.01);
        assert_eq!(x, 0.0, "Avrami at t=0 should be 0");
    }
    #[test]
    fn test_avrami_increases_with_time() {
        let x1 = avrami_transformation(1.0, 2.0, 0.01);
        let x2 = avrami_transformation(5.0, 2.0, 0.01);
        let x3 = avrami_transformation(10.0, 2.0, 0.01);
        assert!(
            x1 < x2 && x2 < x3,
            "Avrami fraction should increase with time"
        );
    }
    #[test]
    fn test_allen_cahn_reduces_free_energy() {
        let mut pf = PhaseFieldOrder::new(1);
        pf.phi[0] = 0.5;
        let phi_initial = pf.phi[0];
        let df_dphi = [1.0];
        pf.update_allen_cahn(0.1, 1.0, &df_dphi);
        assert!(
            pf.phi[0] < phi_initial,
            "Allen-Cahn should decrease φ when dF/dφ > 0"
        );
    }
    #[test]
    fn test_cahn_hilliard_conserves_mass() {
        let nx = 32;
        let c0 = 0.5;
        let mut ch = CahnHilliardField::new(nx, c0);
        let total_before = ch.total_concentration();
        ch.update(0.01, 1.0, 1.0);
        let total_after = ch.total_concentration();
        assert!(
            (total_after - total_before).abs() < 1e-10,
            "Cahn-Hilliard should conserve total concentration: {total_before} → {total_after}"
        );
    }
    #[test]
    fn test_landau_free_energy_at_zero_phi() {
        let f = landau_free_energy(0.0, -1.0, 1.0);
        assert_eq!(f, 0.0, "F(φ=0) should be 0");
    }
    #[test]
    fn test_landau_free_energy_double_well() {
        let a = -1.0_f64;
        let b = 1.0_f64;
        let phi_test = 0.5_f64;
        let f_pos = landau_free_energy(phi_test, a, b);
        let f_neg = landau_free_energy(-phi_test, a, b);
        assert!(
            (f_pos - f_neg).abs() < 1e-12,
            "Landau free energy should be even: {f_pos} vs {f_neg}"
        );
    }
    #[test]
    fn test_phase_transform_type_variants() {
        let t1 = PhaseTransformType::Martensitic;
        let t2 = PhaseTransformType::Bainitic;
        assert_ne!(t1, t2);
        assert_eq!(
            PhaseTransformType::Austenitic,
            PhaseTransformType::Austenitic
        );
    }
    #[test]
    fn test_jmat_pro_model_range() {
        let model = JmatProModel::new(500.0, 300.0, 600.0);
        assert!(
            model.in_martensite_range(400.0),
            "400K should be in martensite range [300, 500]"
        );
        assert!(
            !model.in_martensite_range(200.0),
            "200K should be below martensite range"
        );
        assert!(
            model.above_bainite_start(650.0),
            "650K should be above bainite start 600K"
        );
    }
    #[test]
    fn test_bainite_fraction_zero_at_time_zero() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 80_000.0, 750.0, 500.0);
        assert_eq!(bt.fraction_at(0.0, 650.0), 0.0);
    }
    #[test]
    fn test_bainite_no_fraction_above_bs() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 80_000.0, 750.0, 500.0);
        assert_eq!(bt.fraction_at(1000.0, 800.0), 0.0, "no bainite above B_s");
    }
    #[test]
    fn test_bainite_fraction_increases_with_time() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 80_000.0, 750.0, 500.0);
        let f1 = bt.fraction_at(10.0, 650.0);
        let f2 = bt.fraction_at(100.0, 650.0);
        assert!(
            f2 > f1,
            "bainite fraction must increase with time: {f1} → {f2}"
        );
    }
    #[test]
    fn test_bainite_rate_constant_increases_with_temperature() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 50_000.0, 800.0, 400.0);
        let k_low = bt.rate_constant(500.0);
        let k_high = bt.rate_constant(700.0);
        assert!(k_high > k_low, "k should increase with temperature");
    }
    #[test]
    fn test_bainite_nose_temperature() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 80_000.0, 750.0, 550.0);
        assert!(
            (bt.nose_temperature() - 650.0).abs() < 1e-10,
            "nose at midpoint of Bs–Bf"
        );
    }
    #[test]
    fn test_bainite_transformation_rate_positive() {
        let bt = BainiteTransformation::new(2.0, 1e-3, 80_000.0, 750.0, 500.0);
        let rate = bt.transformation_rate(50.0, 650.0);
        assert!(
            rate > 0.0,
            "transformation rate should be positive, got {rate}"
        );
    }
    #[test]
    fn test_pearlite_fraction_zero_above_ae1() {
        let pt = PearliteTransformation::new(2.5, 1e-4, 90_000.0, 700.0, 1000.0);
        assert_eq!(
            pt.fraction_at(1000.0, 1100.0),
            0.0,
            "no pearlite above A_e1"
        );
    }
    #[test]
    fn test_pearlite_fraction_zero_below_ps() {
        let pt = PearliteTransformation::new(2.5, 1e-4, 90_000.0, 700.0, 1000.0);
        assert_eq!(pt.fraction_at(1000.0, 600.0), 0.0, "no pearlite below P_s");
    }
    #[test]
    fn test_pearlite_fraction_increases_with_time() {
        let pt = PearliteTransformation::new(2.5, 1e-4, 50_000.0, 600.0, 1000.0);
        let f1 = pt.fraction_at(10.0, 800.0);
        let f2 = pt.fraction_at(100.0, 800.0);
        assert!(f2 > f1, "pearlite fraction must increase with time");
    }
    #[test]
    fn test_pearlite_max_fraction_capped_at_one() {
        let pt = PearliteTransformation::new(2.5, 1.0, 1.0, 400.0, 1200.0);
        let f = pt.max_fraction_at(1e10, 800.0);
        assert!(f <= 1.0 + 1e-12, "fraction must not exceed 1.0, got {f}");
    }
    #[test]
    fn test_pearlite_incubation_time_positive() {
        let pt = PearliteTransformation::new(2.5, 1e-4, 50_000.0, 600.0, 1000.0);
        let t_inc = pt.incubation_time(800.0);
        assert!(
            t_inc > 0.0,
            "incubation time should be positive, got {t_inc}"
        );
    }
    #[test]
    fn test_cct_martensite_at_high_rate() {
        let cct = CctDiagram::new(1.0, 5.0, 50.0, 500.0, 820.0, 1.0);
        let micro = cct.predict_microstructure(100.0);
        assert_eq!(micro, "martensite", "fast cooling → martensite");
    }
    #[test]
    fn test_cct_pearlite_at_low_rate() {
        let cct = CctDiagram::new(2.0, 5.0, 50.0, 500.0, 820.0, 1.0);
        let micro = cct.predict_microstructure(0.5);
        assert_eq!(micro, "pearlite", "slow cooling → pearlite");
    }
    #[test]
    fn test_cct_bainite_intermediate_rate() {
        let cct = CctDiagram::new(1.0, 5.0, 50.0, 500.0, 820.0, 1.0);
        let micro = cct.predict_microstructure(8.0);
        assert_eq!(micro, "bainite", "intermediate rate → bainite");
    }
    #[test]
    fn test_cct_martensite_fraction_at_critical_rate() {
        let cct = CctDiagram::new(1.0, 5.0, 50.0, 500.0, 820.0, 1.0);
        let fm = cct.martensite_fraction(50.0);
        assert!(
            (fm - 1.0).abs() < 1e-12,
            "at critical rate, full martensite"
        );
    }
    #[test]
    fn test_cct_martensite_fraction_zero_at_slow_rate() {
        let cct = CctDiagram::new(1.0, 5.0, 50.0, 500.0, 820.0, 1.0);
        let fm = cct.martensite_fraction(0.1);
        assert!(fm < 1e-12, "slow cooling → no martensite");
    }
    #[test]
    fn test_ttt_start_time_infinite_above_ae1() {
        let ttt = TttDiagram::new(820.0, 1.0, 1000.0, 500.0, 750.0);
        let t = ttt.start_time(1100.0);
        assert!(t.is_infinite(), "no transformation above A_e1");
    }
    #[test]
    fn test_ttt_start_time_infinite_at_ms() {
        let ttt = TttDiagram::new(820.0, 1.0, 1000.0, 500.0, 750.0);
        let t = ttt.start_time(500.0);
        assert!(t.is_infinite(), "no austenite transformation at Ms");
    }
    #[test]
    fn test_ttt_bainite_range() {
        let ttt = TttDiagram::new(820.0, 1.0, 1000.0, 500.0, 750.0);
        let (low, high) = ttt.bainite_range();
        assert_eq!(low, 500.0);
        assert_eq!(high, 750.0);
    }
    #[test]
    fn test_ttt_pearlite_range() {
        let ttt = TttDiagram::new(820.0, 1.0, 1000.0, 500.0, 750.0);
        let (low, high) = ttt.pearlite_range();
        assert_eq!(low, 820.0);
        assert_eq!(high, 1000.0);
    }
    #[test]
    fn test_ttt_transformation_avoided_fast_cooling() {
        let ttt = TttDiagram::new(820.0, 10.0, 1000.0, 500.0, 750.0);
        assert!(
            ttt.transformation_avoided(1000.0),
            "fast quench avoids transformation"
        );
    }
    #[test]
    fn test_volume_fraction_tracker_initial_state() {
        let tracker = VolumeFractionTracker::new(vec![
            "austenite".into(),
            "martensite".into(),
            "bainite".into(),
        ]);
        assert!(
            (tracker.fraction(0) - 1.0).abs() < 1e-12,
            "start fully austenitic"
        );
        assert!(tracker.fraction(1).abs() < 1e-12);
        assert!(tracker.fraction(2).abs() < 1e-12);
    }
    #[test]
    fn test_volume_fraction_tracker_update() {
        let mut tracker = VolumeFractionTracker::new(vec!["austenite".into(), "martensite".into()]);
        tracker.update(vec![0.3, 0.7], 10.0, 450.0);
        assert!((tracker.fraction(0) - 0.3).abs() < 1e-12);
        assert!((tracker.fraction(1) - 0.7).abs() < 1e-12);
    }
    #[test]
    fn test_volume_fraction_tracker_total() {
        let mut tracker = VolumeFractionTracker::new(vec!["a".into(), "m".into(), "b".into()]);
        tracker.update(vec![0.5, 0.3, 0.2], 5.0, 600.0);
        assert!((tracker.total() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_volume_fraction_dominant_phase() {
        let mut tracker = VolumeFractionTracker::new(vec![
            "austenite".into(),
            "martensite".into(),
            "bainite".into(),
        ]);
        tracker.update(vec![0.1, 0.8, 0.1], 20.0, 420.0);
        assert_eq!(tracker.dominant_phase(), 1, "martensite should dominate");
    }
    #[test]
    fn test_volume_fraction_history_recorded() {
        let mut tracker = VolumeFractionTracker::new(vec!["a".into(), "m".into()]);
        tracker.update(vec![0.8, 0.2], 1.0, 700.0);
        tracker.update(vec![0.5, 0.5], 2.0, 600.0);
        assert_eq!(tracker.history.len(), 2);
    }
    #[test]
    fn test_phase_props_fully_martensitic() {
        let props = PhaseDependentProperties::medium_carbon_steel();
        let fractions = [1.0_f64, 0.0, 0.0, 0.0];
        let e = props.effective_youngs_modulus(&fractions);
        assert!(
            (e - 205e9).abs() < 1.0,
            "fully martensitic E = 205 GPa, got {e}"
        );
        let sy = props.effective_yield_strength(&fractions);
        assert!(
            (sy - 1500e6).abs() < 1.0,
            "martensitic yield = 1500 MPa, got {sy}"
        );
    }
    #[test]
    fn test_phase_props_rule_of_mixtures() {
        let props = PhaseDependentProperties::medium_carbon_steel();
        let fractions = [0.5_f64, 0.5, 0.0, 0.0];
        let e = props.effective_youngs_modulus(&fractions);
        let expected = 0.5 * 205e9 + 0.5 * 200e9;
        assert!((e - expected).abs() < 1.0, "rule of mixtures E, got {e}");
    }
    #[test]
    fn test_phase_props_hardness_increases_with_martensite() {
        let props = PhaseDependentProperties::medium_carbon_steel();
        let f_martensite = [1.0_f64, 0.0, 0.0, 0.0];
        let f_ferrite = [0.0_f64, 0.0, 1.0, 0.0];
        let hv_m = props.effective_hardness(&f_martensite);
        let hv_f = props.effective_hardness(&f_ferrite);
        assert!(
            hv_m > hv_f,
            "martensite harder than ferrite: {hv_m} vs {hv_f}"
        );
    }
    #[test]
    fn test_phase_props_uts_from_hardness() {
        let props = PhaseDependentProperties::medium_carbon_steel();
        let fractions = [1.0_f64, 0.0, 0.0, 0.0];
        let uts = props.uts_from_hardness(&fractions);
        let expected = 3.45e6 * 700.0;
        assert!((uts - expected).abs() < 1.0, "UTS from hardness, got {uts}");
    }
    #[test]
    fn test_clausius_clapeyron_water_ice_negative_slope() {
        let cc = ClausiusClapeyron::water_ice();
        let slope = cc.slope(cc.t0);
        assert!(
            slope < 0.0,
            "Water/ice dP/dT should be negative, got {slope}"
        );
    }
    #[test]
    fn test_clausius_clapeyron_melting_temperature_increases_with_pressure_for_normal_substance() {
        let cc = ClausiusClapeyron::new(6000.0, 300.0, 101325.0, 1.0e-6);
        let t_low = cc.melting_temperature(101325.0);
        let t_high = cc.melting_temperature(1.0e6);
        assert!(t_high > t_low, "For ΔV>0, higher P → higher T_m");
    }
    #[test]
    fn test_clausius_clapeyron_equilibrium_pressure_at_t0_equals_p0() {
        let cc = ClausiusClapeyron::water_ice();
        let p = cc.equilibrium_pressure(cc.t0);
        assert!(
            (p - cc.p0).abs() < 1e-6,
            "At T₀, equilibrium pressure = P₀, got {p}"
        );
    }
    #[test]
    fn test_clausius_clapeyron_latent_heat_at_reference() {
        let cc = ClausiusClapeyron::water_ice();
        let l = cc.latent_heat_at(cc.t0, 0.0);
        assert!((l - cc.latent_heat).abs() < 1e-10, "L at T₀ = L, got {l}");
    }
    #[test]
    fn test_clausius_clapeyron_latent_heat_increases_with_delta_cp_positive() {
        let cc = ClausiusClapeyron::new(6000.0, 300.0, 101325.0, 1.0e-6);
        let l0 = cc.latent_heat_at(300.0, 10.0);
        let l1 = cc.latent_heat_at(350.0, 10.0);
        assert!(l1 > l0, "L should increase with T when ΔCp > 0");
    }
    #[test]
    fn test_clausius_clapeyron_slope_zero_at_zero_delta_v() {
        let cc = ClausiusClapeyron::new(6000.0, 300.0, 101325.0, 0.0);
        let slope = cc.slope(300.0);
        assert_eq!(slope, 0.0, "Zero ΔV → zero slope");
    }
    #[test]
    fn test_clausius_clapeyron_roundtrip_pressure_temperature() {
        let cc = ClausiusClapeyron::new(6000.0, 300.0, 101325.0, 1.0e-6);
        let t_test = 310.0_f64;
        let p = cc.equilibrium_pressure(t_test);
        let t_back = cc.melting_temperature(p);
        assert!(
            (t_back - t_test).abs() < 1e-6,
            "P→T roundtrip: {t_back} vs {t_test}"
        );
    }
    #[test]
    fn test_lever_rule_at_alpha_composition() {
        let lr = LeverRule::new(0.1, 0.4);
        let f_beta = lr.beta_fraction(0.1);
        assert!(f_beta.abs() < 1e-10, "At x_α, f_β = 0, got {f_beta}");
    }
    #[test]
    fn test_lever_rule_at_beta_composition() {
        let lr = LeverRule::new(0.1, 0.4);
        let f_beta = lr.beta_fraction(0.4);
        assert!(
            (f_beta - 1.0).abs() < 1e-10,
            "At x_β, f_β = 1, got {f_beta}"
        );
    }
    #[test]
    fn test_lever_rule_at_midpoint() {
        let lr = LeverRule::new(0.1, 0.3);
        let f_beta = lr.beta_fraction(0.2);
        assert!(
            (f_beta - 0.5).abs() < 1e-10,
            "At midpoint, f_β = 0.5, got {f_beta}"
        );
    }
    #[test]
    fn test_lever_rule_fractions_sum_to_one() {
        let lr = LeverRule::new(0.05, 0.25);
        let x = 0.15;
        let fa = lr.alpha_fraction(x);
        let fb = lr.beta_fraction(x);
        assert!(
            (fa + fb - 1.0).abs() < 1e-10,
            "Phase fractions must sum to 1"
        );
    }
    #[test]
    fn test_lever_rule_is_two_phase() {
        let lr = LeverRule::new(0.1, 0.4);
        assert!(lr.is_two_phase(0.2), "x=0.2 is in two-phase region");
        assert!(!lr.is_two_phase(0.05), "x=0.05 is outside two-phase region");
        assert!(!lr.is_two_phase(0.5), "x=0.5 is outside two-phase region");
    }
    #[test]
    fn test_lever_rule_tie_line_length() {
        let lr = LeverRule::new(0.1, 0.4);
        assert!((lr.tie_line_length() - 0.3).abs() < 1e-10, "Tie line = 0.3");
    }
    #[test]
    fn test_lever_rule_clamped_outside_range() {
        let lr = LeverRule::new(0.1, 0.4);
        assert_eq!(lr.beta_fraction(0.0), 0.0, "Below α: f_β = 0");
        assert_eq!(lr.beta_fraction(0.5), 1.0, "Above β: f_β = 1");
    }
    #[test]
    fn test_ternary_lever_rule_at_alpha() {
        let f = ternary_lever_rule(0.1, 0.05, 0.1, 0.05, 0.3, 0.1);
        assert!(f.abs() < 1e-10, "At α, f_β = 0");
    }
    #[test]
    fn test_ternary_lever_rule_at_beta() {
        let f = ternary_lever_rule(0.3, 0.1, 0.1, 0.05, 0.3, 0.1);
        assert!((f - 1.0).abs() < 1e-10, "At β, f_β = 1");
    }
    #[test]
    fn test_spinodal_f_double_prime_positive_at_low_omega() {
        let spd = SpinocalDecomposition::new(1000.0, 1e-10, 1e-5);
        let f2 = spd.f_double_prime(0.5, 1000.0);
        assert!(f2 > 0.0, "f'' > 0 for stable mixture, got {f2}");
    }
    #[test]
    fn test_spinodal_f_double_prime_negative_at_high_omega() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        let f2 = spd.f_double_prime(0.5, 300.0);
        assert!(f2 < 0.0, "f'' < 0 in spinodal region, got {f2}");
    }
    #[test]
    fn test_spinodal_is_spinodal_at_center_high_omega() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        assert!(
            spd.is_spinodal(0.5, 300.0),
            "c=0.5 should be in spinodal at high Ω"
        );
    }
    #[test]
    fn test_spinodal_boundary_symmetric() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        if let Some((c_low, c_high)) = spd.spinodal_boundary(300.0) {
            assert!(
                (c_low + c_high - 1.0).abs() < 1e-10,
                "Spinodal boundary symmetric about 0.5"
            );
            assert!(c_low < 0.5 && c_high > 0.5, "Boundary straddles 0.5");
        } else {
            panic!("Expected spinodal boundary at T=300");
        }
    }
    #[test]
    fn test_spinodal_no_boundary_above_critical_temperature() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        let tc = spd.critical_temperature();
        let result = spd.spinodal_boundary(tc * 1.5);
        assert!(result.is_none(), "No spinodal above T_c");
    }
    #[test]
    fn test_spinodal_critical_temperature_formula() {
        let omega = 40000.0_f64;
        let spd = SpinocalDecomposition::new(omega, 1e-10, 1e-5);
        let tc = spd.critical_temperature();
        let r = 8.314_462_618;
        let expected = omega / (2.0 * r);
        assert!(
            (tc - expected).abs() / expected < 1e-10,
            "T_c formula: {tc}"
        );
    }
    #[test]
    fn test_spinodal_fastest_wavelength_exists_in_spinodal() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        let lambda = spd.fastest_growing_wavelength(0.5, 300.0);
        assert!(
            lambda.is_some(),
            "Should have fastest wavelength in spinodal region"
        );
        assert!(lambda.unwrap() > 0.0, "Wavelength must be positive");
    }
    #[test]
    fn test_spinodal_amplification_positive_in_spinodal() {
        let spd = SpinocalDecomposition::new(50000.0, 1e-10, 1e-5);
        if let Some(lambda) = spd.fastest_growing_wavelength(0.5, 300.0) {
            let q = 2.0 * std::f64::consts::PI / lambda;
            let r = spd.amplification_factor(q, 0.5, 300.0, 1.0e-20);
            assert!(
                r > 0.0,
                "Amplification factor should be positive at fastest q, got {r}"
            );
        }
    }
    #[test]
    fn test_latent_heat_volumetric_steel() {
        let lh = LatentHeatModel::steel_solidification();
        let v = lh.volumetric_latent_heat();
        let expected = 271_000.0 * 7800.0;
        assert!(
            (v - expected).abs() < 1.0,
            "Volumetric latent heat steel: {v}"
        );
    }
    #[test]
    fn test_latent_heat_enthalpy_release_full_transformation() {
        let lh = LatentHeatModel::aluminum_solidification();
        let h = lh.enthalpy_release(1.0);
        assert!(
            (h - lh.volumetric_latent_heat()).abs() < 1.0,
            "Full transformation releases all latent heat"
        );
    }
    #[test]
    fn test_latent_heat_enthalpy_release_half_transformation() {
        let lh = LatentHeatModel::steel_solidification();
        let h_half = lh.enthalpy_release(0.5);
        let h_full = lh.enthalpy_release(1.0);
        assert!(
            (h_half - 0.5 * h_full).abs() < 1.0,
            "Half transform → half latent heat"
        );
    }
    #[test]
    fn test_latent_heat_adiabatic_rise_positive() {
        let lh = LatentHeatModel::steel_solidification();
        let dt = lh.adiabatic_temperature_rise(1.0, 500.0);
        assert!(dt > 0.0, "Adiabatic temperature rise should be positive");
    }
    #[test]
    fn test_latent_heat_adiabatic_rise_formula() {
        let lh = LatentHeatModel::new(300_000.0, 8000.0);
        let cp = 500.0;
        let dt = lh.adiabatic_temperature_rise(1.0, cp);
        let expected = 300_000.0 / cp;
        assert!(
            (dt - expected).abs() < 1e-6,
            "ΔT = L/Cp: {dt} vs {expected}"
        );
    }
    #[test]
    fn test_latent_heat_fraction_released() {
        let lh = LatentHeatModel::steel_solidification();
        let df = lh.fraction_released(0.2, 0.8);
        assert!(
            (df - 0.6).abs() < 1e-10,
            "Fraction released = 0.6, got {df}"
        );
    }
    #[test]
    fn test_latent_heat_zero_cp_returns_zero() {
        let lh = LatentHeatModel::steel_solidification();
        let dt = lh.adiabatic_temperature_rise(1.0, 0.0);
        assert_eq!(dt, 0.0, "Zero Cp → zero ΔT");
    }
    #[test]
    fn test_jmak_extended_fraction_zero_at_time_zero() {
        let jmak = JmakExtended::new(1e-10, 2.0, 80_000.0, 1.0);
        assert_eq!(jmak.fraction(0.0, 900.0), 0.0, "f(0) = 0");
    }
    #[test]
    fn test_jmak_extended_fraction_approaches_one() {
        let jmak = JmakExtended::new(1e-5, 2.0, 1.0, 1.0);
        let f = jmak.fraction(1e10, 1000.0);
        assert!((f - 1.0).abs() < 1e-6, "f → 1 at large t, got {f}");
    }
    #[test]
    fn test_jmak_extended_fraction_monotonically_increasing() {
        let jmak = JmakExtended::new(1e-10, 2.0, 50_000.0, 1.0);
        let f1 = jmak.fraction(10.0, 900.0);
        let f2 = jmak.fraction(100.0, 900.0);
        let f3 = jmak.fraction(1000.0, 900.0);
        assert!(f2 > f1 && f3 > f2, "JMAK fraction must increase");
    }
    #[test]
    fn test_jmak_extended_higher_temp_faster_transform() {
        let jmak = JmakExtended::new(1e-10, 2.0, 80_000.0, 1.0);
        let f_low = jmak.fraction(100.0, 700.0);
        let f_high = jmak.fraction(100.0, 1200.0);
        assert!(f_high > f_low, "Higher T → faster transformation");
    }
    #[test]
    fn test_jmak_extended_time_to_fraction_roundtrip() {
        let jmak = JmakExtended::new(1e-10, 2.0, 50_000.0, 1.0);
        let target = 0.5;
        let t = jmak.time_to_fraction(target, 1000.0);
        let f_back = jmak.fraction(t, 1000.0);
        assert!(
            (f_back - target).abs() < 1e-6,
            "JMAK roundtrip: {f_back} vs {target}"
        );
    }
    #[test]
    fn test_jmak_extended_ttt_start_before_finish() {
        let jmak = JmakExtended::new(1e-10, 2.0, 50_000.0, 1.0);
        let t1_pct = jmak.time_to_fraction(0.01, 1000.0);
        let t99_pct = jmak.time_to_fraction(0.99, 1000.0);
        assert!(
            t1_pct.is_finite() && t99_pct.is_finite(),
            "Times must be finite"
        );
        assert!(
            t99_pct > t1_pct,
            "99% transform takes longer than 1%: t1={t1_pct}, t99={t99_pct}"
        );
    }
    #[test]
    fn test_jmak_extended_transformation_rate_positive() {
        let jmak = JmakExtended::new(1e-10, 2.0, 50_000.0, 1.0);
        let rate = jmak.transformation_rate(100.0, 1000.0);
        assert!(
            rate > 0.0,
            "Transformation rate must be positive, got {rate}"
        );
    }
    #[test]
    fn test_jmak_extended_impingement_reduces_fraction() {
        let jmak_full = JmakExtended::new(1e-10, 2.0, 50_000.0, 1.0);
        let jmak_partial = JmakExtended::new(1e-10, 2.0, 50_000.0, 0.5);
        let f_full = jmak_full.fraction(1000.0, 1000.0);
        let f_partial = jmak_partial.fraction(1000.0, 1000.0);
        assert!(
            f_partial < f_full,
            "Lower impingement → lower fraction at same time"
        );
    }
}
#[cfg(test)]
mod tests_new {

    use crate::EutecticTransformation;

    use crate::ScheilSolidification;

    use crate::StressAidedMartensite;

    use crate::TttCcurve;

    use crate::TttExtended;

    #[test]
    fn test_scheil_liquid_composition_at_zero_fraction() {
        let sch = ScheilSolidification::new(4.0, 0.17);
        let c_l = sch.liquid_composition(0.0);
        assert!((c_l - 4.0).abs() < 1e-6, "C_L(f_s=0) = C_0, got {c_l}");
    }
    #[test]
    fn test_scheil_liquid_composition_increases_with_fs() {
        let sch = ScheilSolidification::new(4.0, 0.17);
        let c_l0 = sch.liquid_composition(0.0);
        let c_l5 = sch.liquid_composition(0.5);
        assert!(
            c_l5 > c_l0,
            "C_L increases with f_s for k < 1: {c_l0} → {c_l5}"
        );
    }
    #[test]
    fn test_scheil_solid_front_composition_proportional_to_liquid() {
        let sch = ScheilSolidification::al_cu();
        let f_s = 0.3;
        let c_l = sch.liquid_composition(f_s);
        let c_s = sch.solid_composition_at_front(f_s);
        assert!(
            (c_s - sch.k_partition * c_l).abs() < 1e-10,
            "C_S* = k * C_L"
        );
    }
    #[test]
    fn test_scheil_mean_solid_composition_less_than_c0() {
        let sch = ScheilSolidification::new(4.0, 0.17);
        let c_s_mean = sch.mean_solid_composition(0.5);
        assert!(
            c_s_mean < sch.c0,
            "Mean solid composition < C_0 for k < 1: {c_s_mean}"
        );
    }
    #[test]
    fn test_scheil_eutectic_solid_fraction_increases_toward_one() {
        let sch = ScheilSolidification::al_cu();
        if let Some(f_s) = sch.eutectic_solid_fraction(33.0) {
            assert!(f_s > 0.0 && f_s < 1.0, "Eutectic f_s in (0,1): {f_s}");
        }
    }
    #[test]
    fn test_scheil_solidification_range() {
        let sch = ScheilSolidification::new(4.0, 0.17);
        let (t_l, t_e) = sch.solidification_range(933.0, -3.4, 33.0);
        assert!(t_l > t_e, "T_liquidus > T_eutectic: {t_l} > {t_e}");
    }
    #[test]
    fn test_scheil_fe_c_liquid_composition() {
        let sch = ScheilSolidification::fe_c();
        let c_l = sch.liquid_composition(0.3);
        assert!(c_l > sch.c0, "Fe-C: C_L increases with solidification");
    }
    #[test]
    fn test_ttt_extended_martensite_fraction_at_ms() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let fm = ttt.martensite_fraction(500.0);
        assert!(fm.abs() < 1e-10, "No martensite at Ms exactly, got {fm}");
    }
    #[test]
    fn test_ttt_extended_martensite_increases_below_ms() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let fm1 = ttt.martensite_fraction(490.0);
        let fm2 = ttt.martensite_fraction(400.0);
        assert!(fm2 > fm1, "More martensite at lower T: {fm1} vs {fm2}");
    }
    #[test]
    fn test_ttt_extended_pearlite_zero_above_ae1() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let fp = ttt.pearlite_fraction(1e6, 1100.0);
        assert!(fp.abs() < 1e-12, "No pearlite above A_e1");
    }
    #[test]
    fn test_ttt_extended_bainite_fraction_increases_with_time() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let fb1 = ttt.bainite_fraction(1.0, 650.0);
        let fb2 = ttt.bainite_fraction(100.0, 650.0);
        assert!(
            fb2 >= fb1,
            "Bainite fraction non-decreasing with time: {fb1} vs {fb2}"
        );
    }
    #[test]
    fn test_ttt_extended_dominant_product_austenite_above_ae1() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let product = ttt.dominant_product(1.0, 1100.0);
        assert_eq!(product, "austenite", "Above A_e1 → austenite");
    }
    #[test]
    fn test_ttt_extended_dominant_product_martensite_below_ms() {
        let ttt = TttExtended::new(820.0, 1.0, 650.0, 0.5, 500.0, 1000.0, 3.0, 2.5);
        let product = ttt.dominant_product(1.0, 400.0);
        assert_eq!(product, "martensite", "Below Ms → martensite");
    }
    #[test]
    fn test_eutectic_alpha_plus_beta_fractions_sum_to_one() {
        let eut = EutecticTransformation::al_si();
        let fa = eut.alpha_fraction();
        let fb = eut.beta_fraction();
        assert!(
            (fa + fb - 1.0).abs() < 1e-10,
            "Phase fractions sum to 1: {fa} + {fb}"
        );
    }
    #[test]
    fn test_eutectic_al_si_alpha_fraction_in_range() {
        let eut = EutecticTransformation::al_si();
        let fa = eut.alpha_fraction();
        assert!(
            fa > 0.0 && fa < 1.0,
            "Al-Si α fraction should be in (0,1): {fa}"
        );
    }
    #[test]
    fn test_eutectic_lamellar_spacing_decreases_with_velocity() {
        let eut = EutecticTransformation::al_si();
        let lam_slow = eut.lamellar_spacing(1e-6);
        let lam_fast = eut.lamellar_spacing(1e-4);
        assert!(
            lam_fast < lam_slow,
            "Faster growth → finer lamellar spacing"
        );
    }
    #[test]
    fn test_eutectic_growth_velocity_decreases_with_spacing() {
        let eut = EutecticTransformation::al_si();
        let v_fine = eut.growth_velocity(1e-7);
        let v_coarse = eut.growth_velocity(1e-6);
        assert!(v_fine > v_coarse, "Finer spacing → higher velocity");
    }
    #[test]
    fn test_eutectic_spacing_velocity_product_constant() {
        let eut = EutecticTransformation::al_si();
        let v = 1e-5;
        let lam = eut.lamellar_spacing(v);
        let product = lam * lam * v;
        assert!(
            (product - eut.k_jh).abs() / eut.k_jh < 1e-6,
            "λ²·v = K_JH: {product} vs {}",
            eut.k_jh
        );
    }
    #[test]
    fn test_eutectic_is_eutectic_composition() {
        let eut = EutecticTransformation::al_si();
        assert!(
            eut.is_eutectic_composition(eut.c_eutectic, 0.01),
            "Exact eutectic composition matches"
        );
        assert!(
            !eut.is_eutectic_composition(5.0, 0.01),
            "Non-eutectic composition does not match"
        );
    }
    #[test]
    fn test_ttt_ccurve_minimum_time_at_nose() {
        let ccurve = TttCcurve::new(820.0, 1.0, 0.001, 3.0);
        let tau_nose = ccurve.incubation_time(820.0);
        let tau_off = ccurve.incubation_time(700.0);
        assert!(
            tau_off > tau_nose,
            "Incubation time minimum at nose, τ_nose={tau_nose}, τ_off={tau_off}"
        );
    }
    #[test]
    fn test_ttt_ccurve_fraction_monotonically_increasing() {
        let ccurve = TttCcurve::new(820.0, 1.0, 0.001, 3.0);
        let f1 = ccurve.fraction(1.0, 820.0);
        let f2 = ccurve.fraction(10.0, 820.0);
        let f3 = ccurve.fraction(100.0, 820.0);
        assert!(f2 > f1 && f3 > f2, "Fraction increases with time");
    }
    #[test]
    fn test_ttt_ccurve_time_to_fraction_roundtrip() {
        let ccurve = TttCcurve::new(820.0, 1.0, 0.001, 3.0);
        let target = 0.5;
        let t = ccurve.time_to_fraction(target, 820.0);
        let f_back = ccurve.fraction(t, 820.0);
        assert!(
            (f_back - target).abs() < 1e-6,
            "TTT C-curve roundtrip: {f_back} vs {target}"
        );
    }
    #[test]
    fn test_ttt_ccurve_nose_temperature_correct() {
        let ccurve = TttCcurve::new(820.0, 1.0, 0.001, 3.0);
        assert!((ccurve.nose_temperature() - 820.0).abs() < 1e-10);
    }
    #[test]
    fn test_stress_aided_ms_shifts_with_stress() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        let ms_no_stress = sam.ms_temperature_under_stress(0.0);
        let ms_stressed = sam.ms_temperature_under_stress(100.0e6);
        assert_eq!(ms_no_stress, 500.0, "No-stress Ms = 500 K");
        assert!(ms_stressed != ms_no_stress, "Ms shifts with stress");
    }
    #[test]
    fn test_stress_aided_shear_band_fraction_zero_at_zero_strain() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        assert_eq!(sam.shear_band_fraction(0.0), 0.0);
    }
    #[test]
    fn test_stress_aided_shear_band_fraction_increases_with_strain() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        let f1 = sam.shear_band_fraction(0.1);
        let f2 = sam.shear_band_fraction(0.5);
        assert!(f2 > f1, "Shear band fraction increases with strain");
    }
    #[test]
    fn test_stress_aided_martensite_from_strain_increases() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        let fm1 = sam.martensite_fraction_from_strain(0.1);
        let fm2 = sam.martensite_fraction_from_strain(0.5);
        assert!(
            fm2 > fm1,
            "Martensite fraction increases with strain: {fm1} vs {fm2}"
        );
    }
    #[test]
    fn test_stress_aided_transformation_mode_thermal_below_ms() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        let mode = sam.transformation_mode(400.0, 0.0);
        assert_eq!(mode, "thermal", "Below Ms₀: thermal transformation");
    }
    #[test]
    fn test_stress_aided_transformation_mode_none_above_both() {
        let sam = StressAidedMartensite::new(500.0, 1.0e-7, 4.0, 4.5, 4.0);
        let mode = sam.transformation_mode(600.0, 0.0);
        assert_eq!(mode, "none", "Above Ms(σ): no transformation");
    }
}
