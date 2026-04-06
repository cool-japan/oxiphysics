//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::smart_materials::BimorphActuator;

    use crate::smart_materials::DielectricElastomer;

    use crate::smart_materials::ElectrorheologicalFluid;

    use crate::smart_materials::HydrogelActuator;

    use crate::smart_materials::MagnetorheologicalFluid;

    use crate::smart_materials::PiezoActuator;

    use crate::smart_materials::PolymerActuator;

    use crate::smart_materials::ShapeMemoryAlloy;

    use crate::smart_materials::SmaController;
    use crate::smart_materials::SmaPhase;

    use crate::smart_materials::SmartStructure;

    #[test]
    fn sma_nitinol_initial_phase() {
        let sma = ShapeMemoryAlloy::nitinol();
        assert_eq!(sma.phase(), SmaPhase::Martensite);
    }
    #[test]
    fn sma_phase_transition_heating() {
        let mut sma = ShapeMemoryAlloy::nitinol();
        sma.update_phase(400.0, 0.0);
        assert_eq!(sma.phase(), SmaPhase::Austenite);
    }
    #[test]
    fn sma_elastic_modulus_changes_with_phase() {
        let mut sma = ShapeMemoryAlloy::nitinol();
        let e_m = sma.elastic_modulus();
        sma.update_phase(400.0, 0.0);
        let e_a = sma.elastic_modulus();
        assert!(e_m != e_a, "modulus should change with phase");
    }
    #[test]
    fn sma_martensite_fraction_zero_at_high_temp() {
        let mut sma = ShapeMemoryAlloy::nitinol();
        sma.update_phase(400.0, 0.0);
        assert!(sma.xi < 0.01, "xi={}", sma.xi);
    }
    #[test]
    fn sma_recovery_stress_finite() {
        let sma = ShapeMemoryAlloy::nitinol();
        let sigma = sma.recovery_stress(0.04);
        assert!(sigma.is_finite(), "sigma={sigma}");
    }
    #[test]
    fn mr_fluid_yield_stress_zero_field() {
        let mr = MagnetorheologicalFluid::new(1.0, 100.0, 200.0, 1.5, 0.3);
        let tau_y = mr.yield_stress(0.0);
        assert!((tau_y - 100.0).abs() < 1e-10, "tau_y={tau_y}");
    }
    #[test]
    fn mr_fluid_yield_stress_increases_with_field() {
        let mr = MagnetorheologicalFluid::new(1.0, 100.0, 200.0, 1.5, 0.3);
        let t1 = mr.yield_stress(1000.0);
        let t2 = mr.yield_stress(5000.0);
        assert!(t2 > t1, "yield stress should increase with field");
    }
    #[test]
    fn mr_fluid_shear_stress_finite() {
        let mr = MagnetorheologicalFluid::new(1.0, 100.0, 200.0, 1.5, 0.3);
        let tau = mr.shear_stress(100.0, 1000.0);
        assert!(tau.is_finite() && tau > 0.0, "tau={tau}");
    }
    #[test]
    fn mr_fluid_mason_number_high_without_field() {
        let mr = MagnetorheologicalFluid::new(1.0, 100.0, 200.0, 1.5, 0.3);
        let mn = mr.mason_number(100.0, 0.0);
        assert!(mn.is_infinite(), "should be infinite without field");
    }
    #[test]
    fn er_fluid_yield_stress_increases_with_field() {
        let er = ElectrorheologicalFluid::new(0.1, 1e-5, 1e6, 100.0, 2.5);
        let t1 = er.yield_stress(1e4);
        let t2 = er.yield_stress(1e5);
        assert!(t2 > t1, "t2={t2} > t1={t1}");
    }
    #[test]
    fn er_fluid_saturation_limits_stress() {
        let er = ElectrorheologicalFluid::new(0.1, 1e-5, 1e5, 100.0, 2.5);
        let t1 = er.yield_stress(1e5);
        let t2 = er.yield_stress(1e8);
        assert!(
            (t1 - t2).abs() < 1e-6,
            "stress should saturate: t1={t1}, t2={t2}"
        );
    }
    #[test]
    fn er_fluid_beta_parameter_range() {
        let er = ElectrorheologicalFluid::new(0.1, 1e-5, 1e6, 100.0, 2.5);
        let beta = er.beta_parameter();
        assert!(beta > -1.0 && beta < 1.0, "beta={beta}");
    }
    #[test]
    fn ipmc_steady_state_deflection_proportional_to_voltage() {
        let ipmc = PolymerActuator::new(0.01, 0.1, 0.05, 2e-4, 0.005);
        let d1 = ipmc.steady_state_deflection(1.0);
        let d2 = ipmc.steady_state_deflection(2.0);
        assert!((d2 / d1 - 2.0).abs() < 1e-10, "ratio={}", d2 / d1);
    }
    #[test]
    fn ipmc_transient_approaches_steady_state() {
        let ipmc = PolymerActuator::new(0.01, 0.1, 0.05, 2e-4, 0.005);
        let d_ss = ipmc.steady_state_deflection(1.0);
        let d_t = ipmc.transient_deflection(1.0, 10.0 * 0.1);
        assert!((d_t - d_ss).abs() / d_ss < 0.01, "d_t={d_t}, d_ss={d_ss}");
    }
    #[test]
    fn ipmc_blocking_force_positive() {
        let ipmc = PolymerActuator::new(0.01, 0.1, 0.05, 2e-4, 0.005);
        assert!(ipmc.blocking_force(1.0) > 0.0);
    }
    #[test]
    fn piezo_free_stroke_linear() {
        let piezo = PiezoActuator::new(1e-8, 1.0, 1e6, 1000.0, 1e-8);
        let s1 = piezo.free_stroke(100.0);
        let s2 = piezo.free_stroke(200.0);
        assert!((s2 / s1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn piezo_blocking_force_positive() {
        let piezo = PiezoActuator::new(1e-8, 1.0, 1e6, 1000.0, 1e-8);
        assert!(piezo.blocking_force(100.0) > 0.0);
    }
    #[test]
    fn piezo_output_force_at_free_displacement() {
        let piezo = PiezoActuator::new(1e-8, 1.0, 1e6, 1000.0, 1e-8);
        let free = piezo.free_stroke(100.0);
        let f = piezo.output_force(100.0, free);
        assert!(f.abs() < 1e-8, "force at free displacement ≈ 0: {f}");
    }
    #[test]
    fn dea_actuation_strain_positive() {
        let dea = DielectricElastomer::new(3.0, 1e-4, 1e-4, 1e5, 1e8);
        let strain = dea.actuation_strain(1000.0);
        assert!(strain > 0.0, "strain={strain}");
    }
    #[test]
    fn dea_area_strain_double_thickness_strain() {
        let dea = DielectricElastomer::new(3.0, 1e-4, 1e-4, 1e5, 1e8);
        let es = dea.actuation_strain(500.0);
        let ea = dea.area_strain(500.0);
        assert!(
            (ea - 2.0 * es).abs() < 1e-10,
            "area_strain={ea}, 2*es={}",
            2.0 * es
        );
    }
    #[test]
    fn dea_pressure_below_breakdown() {
        let dea = DielectricElastomer::new(3.0, 1e-4, 1e-4, 1e5, 1e8);
        let p1 = dea.actuation_pressure(1e6);
        let p2 = dea.actuation_pressure(1e12);
        assert!(
            (p1 - p2).abs() < 1.0,
            "pressure should saturate at breakdown"
        );
    }
    #[test]
    fn hydrogel_swelling_ratio_minimum_one() {
        let hg = HydrogelActuator::new(0.5, 1e3, 10.0, 0.01, 7.0, 0.1, 1e4);
        let q = hg.swelling_ratio(300.0, 7.0);
        assert!(q >= 1.0, "q={q}");
    }
    #[test]
    fn hydrogel_force_positive() {
        let hg = HydrogelActuator::new(0.5, 1e3, 10.0, 0.01, 7.0, 0.1, 1e4);
        let f = hg.force(310.0, 9.0, 1e-4);
        assert!(f.is_finite(), "f={f}");
    }
    #[test]
    fn bimorph_deflection_positive_for_positive_dt() {
        let bimorph = BimorphActuator::new(0.1, 1e-4, 1e-4, 200e9, 70e9, 12e-6, 23e-6);
        let d = bimorph.tip_deflection(100.0);
        assert!(d.is_finite(), "d={d}");
    }
    #[test]
    fn bimorph_deflection_zero_for_zero_dt() {
        let bimorph = BimorphActuator::new(0.1, 1e-4, 1e-4, 200e9, 70e9, 12e-6, 23e-6);
        let d = bimorph.tip_deflection(0.0);
        assert!(d.abs() < 1e-20, "d={d}");
    }
    #[test]
    fn bimorph_neutral_axis_between_layers() {
        let bimorph = BimorphActuator::new(0.1, 1e-4, 1e-4, 200e9, 70e9, 12e-6, 23e-6);
        let na = bimorph.neutral_axis();
        let total = bimorph.t1 + bimorph.t2;
        assert!(na > 0.0 && na < total, "na={na}, total={total}");
    }
    #[test]
    fn sma_controller_zero_power_at_target_martensite() {
        let sma = ShapeMemoryAlloy::nitinol();
        let mut ctrl = SmaController::new(sma, 400.0, 250.0, 1.0, 100.0);
        ctrl.set_target(1.0);
        let power = ctrl.compute_power();
        assert!(power >= 0.0, "power={power}");
    }
    #[test]
    fn sma_controller_max_power_for_austenite_target() {
        let sma = ShapeMemoryAlloy::nitinol();
        let mut ctrl = SmaController::new(sma, 400.0, 250.0, 1.0, 100.0);
        ctrl.set_target(0.0);
        let power = ctrl.compute_power();
        assert!(
            power > 0.0,
            "should heat to transition to austenite: power={power}"
        );
    }
    #[test]
    fn sma_controller_step_changes_temperature() {
        let sma = ShapeMemoryAlloy::nitinol();
        let mut ctrl = SmaController::new(sma, 400.0, 250.0, 1.0, 100.0);
        ctrl.set_target(0.0);
        let t_before = ctrl.temperature;
        ctrl.step(1.0, 1.0, 250.0);
        let t_after = ctrl.temperature;
        assert!(t_after != t_before, "temperature should change");
    }
    #[test]
    fn smart_structure_displacement_zero_at_rest() {
        let ss = SmartStructure::new(
            2,
            vec![10.0, 30.0],
            vec![0.01, 0.01],
            vec![0.25, 0.75],
            vec![0.25, 0.75],
        );
        assert!(ss.displacement_at(0.5).abs() < 1e-12);
    }
    #[test]
    fn smart_structure_step_changes_state() {
        let mut ss =
            SmartStructure::new(2, vec![10.0, 30.0], vec![0.01, 0.01], vec![0.5], vec![0.5]);
        ss.step(&[1.0], 0.001);
        let d = ss.displacement_at(0.5);
        assert!(d.is_finite(), "d={d}");
    }
    #[test]
    fn smart_structure_set_control_gains() {
        let mut ss = SmartStructure::new(
            3,
            vec![5.0, 15.0, 30.0],
            vec![0.02, 0.02, 0.02],
            vec![0.5],
            vec![0.5],
        );
        ss.set_control_gains(vec![2.0, 0.5, 0.0]);
        assert!((ss.control_gains[0] - 2.0).abs() < 1e-10);
        assert!((ss.control_gains[2] - 0.0).abs() < 1e-10);
    }
}
/// Boltzmann constant used in SMA thermal computations.
#[allow(dead_code)]
pub(super) const BOLTZMANN: f64 = 1.380649e-23;
#[cfg(test)]
mod smart_new_tests {

    use crate::smart_materials::ElectroactivePoly;

    use crate::smart_materials::FerroelectricHysteresis;

    use crate::smart_materials::HydrogelSwelling;
    use crate::smart_materials::MagneticShape;

    use crate::smart_materials::MagnetostrictiveMaterial;

    use crate::smart_materials::PiezoelectricCoupling;

    use crate::smart_materials::ShapeMemoryAlloy;
    use crate::smart_materials::SmaActuator;

    use crate::smart_materials::SmartComposite;

    use crate::smart_materials::ThermochromicResponse;
    use std::f64::consts::PI;
    #[test]
    fn sma_actuator_initial_xi_matches_sma() {
        let sma = ShapeMemoryAlloy::nitinol();
        let xi = sma.xi;
        let act = SmaActuator::new(sma, 1e-6, 0.1);
        assert!((act.total_xi() - xi).abs() < 1e-10, "xi mismatch");
    }
    #[test]
    fn sma_actuator_stroke_positive() {
        let sma = ShapeMemoryAlloy::nitinol();
        let act = SmaActuator::new(sma, 1e-6, 0.1);
        let s = act.stroke(0.0);
        assert!(s > 0.0, "stroke={s}");
    }
    #[test]
    fn sma_actuator_update_changes_state() {
        let sma = ShapeMemoryAlloy::nitinol();
        let mut act = SmaActuator::new(sma, 1e-6, 0.1);
        act.update(400.0, 0.0);
        assert!(
            act.total_xi() < 0.5,
            "should be mostly austenite after heating"
        );
    }
    #[test]
    fn sma_actuator_force_finite() {
        let sma = ShapeMemoryAlloy::nitinol();
        let mut act = SmaActuator::new(sma, 1e-6, 0.1);
        act.update(300.0, 0.02);
        let f = act.force();
        assert!(f.is_finite(), "force={f}");
    }
    #[test]
    fn sma_actuator_stroke_zero_same_xi() {
        let sma = ShapeMemoryAlloy::nitinol();
        let act = SmaActuator::new(sma, 1e-6, 0.1);
        let s = act.stroke(act.total_xi());
        assert!(s.abs() < 1e-12, "stroke to same xi should be 0, got {s}");
    }
    #[test]
    fn piezo_coupling_pzt5a_d33_positive() {
        let p = PiezoelectricCoupling::pzt5a();
        assert!(p.d33 > 0.0, "d33={}", p.d33);
    }
    #[test]
    fn piezo_coupling_d31_negative() {
        let p = PiezoelectricCoupling::pzt5a();
        assert!(p.d31 < 0.0, "d31 should be negative for PZT, got {}", p.d31);
    }
    #[test]
    fn piezo_coupling_strain_33_proportional_to_voltage() {
        let p = PiezoelectricCoupling::pzt5a();
        let s1 = p.strain_from_voltage_33(100.0, 1e-3);
        let s2 = p.strain_from_voltage_33(200.0, 1e-3);
        assert!((s2 / s1 - 2.0).abs() < 1e-10, "s2/s1={}", s2 / s1);
    }
    #[test]
    fn piezo_coupling_charge_from_stress_positive() {
        let p = PiezoelectricCoupling::pzt5a();
        let q = p.charge_from_stress_33(1e6);
        assert!(q > 0.0, "q={q}");
    }
    #[test]
    fn piezo_coupling_voltage_from_stress_positive() {
        let p = PiezoelectricCoupling::pzt5a();
        let v = p.voltage_from_stress(1e6, 1e-3);
        assert!(v > 0.0, "v={v}");
    }
    #[test]
    fn piezo_coupling_harvested_energy_positive() {
        let p = PiezoelectricCoupling::pzt5a();
        let e = p.harvested_energy_density(1e6);
        assert!(e > 0.0, "e={e}");
    }
    #[test]
    fn magnetostrictive_zero_field_zero_strain() {
        let m = MagnetostrictiveMaterial::terfenol_d();
        let s = m.strain_magnitude(0.0);
        assert!(s.abs() < 1e-15, "s={s}");
    }
    #[test]
    fn magnetostrictive_strain_positive() {
        let m = MagnetostrictiveMaterial::terfenol_d();
        let s = m.strain_magnitude(50e3);
        assert!(s > 0.0, "s={s}");
    }
    #[test]
    fn magnetostrictive_strain_saturates() {
        let m = MagnetostrictiveMaterial::terfenol_d();
        let s1 = m.strain_magnitude(1e6);
        let s2 = m.strain_magnitude(1e8);
        assert!(
            (s1 - m.lambda_sat).abs() < 0.2 * m.lambda_sat,
            "s1={s1}, lambda_sat={}",
            m.lambda_sat
        );
        assert!(
            (s2 - m.lambda_sat).abs() < 0.2 * m.lambda_sat,
            "s2={s2}, lambda_sat={}",
            m.lambda_sat
        );
    }
    #[test]
    fn magnetostrictive_blocked_stress_positive() {
        let m = MagnetostrictiveMaterial::terfenol_d();
        let sigma = m.blocked_stress(50e3);
        assert!(sigma > 0.0, "sigma={sigma}");
    }
    #[test]
    fn eap_ipmc_low_voltage() {
        let eap = ElectroactivePoly::ipmc();
        let s = eap.steady_state_strain(0.0);
        assert!(
            s < eap.max_strain * 0.6,
            "strain at 0V should be below half-max: {s}"
        );
    }
    #[test]
    fn eap_strain_increases_with_voltage() {
        let eap = ElectroactivePoly::ipmc();
        let s1 = eap.steady_state_strain(0.5);
        let s2 = eap.steady_state_strain(3.0);
        assert!(s2 > s1, "s2={s2} > s1={s1}");
    }
    #[test]
    fn eap_transient_zero_at_t0() {
        let eap = ElectroactivePoly::ipmc();
        let s = eap.transient_strain(2.0, 0.0);
        assert!(s.abs() < 1e-10, "s={s}");
    }
    #[test]
    fn eap_blocking_stress_positive() {
        let eap = ElectroactivePoly::dielectric_elastomer();
        let stress = eap.blocking_stress(2000.0);
        assert!(stress > 0.0, "stress={stress}");
    }
    #[test]
    fn eap_update_changes_state() {
        let mut eap = ElectroactivePoly::ipmc();
        eap.update(2.0, 0.5);
        assert!(eap.state > 0.0, "state should have increased");
    }
    #[test]
    fn hydrogel_swelling_ratio_above_one() {
        let hg = HydrogelSwelling::pnipam();
        let q = hg.swelling_ratio();
        assert!(q > 1.0, "q={q}");
    }
    #[test]
    fn hydrogel_linear_strain_positive() {
        let hg = HydrogelSwelling::pnipam();
        let eps = hg.linear_strain();
        assert!(eps > 0.0, "eps={eps}");
    }
    #[test]
    fn hydrogel_mixing_free_energy_finite() {
        let hg = HydrogelSwelling::pnipam();
        let f = hg.mixing_free_energy(0.1);
        assert!(f.is_finite(), "f={f}");
    }
    #[test]
    fn hydrogel_osmotic_pressure_positive_at_low_phi() {
        let hg = HydrogelSwelling::new(0.5, 0.05, 100.0, 18e-6, 298.0);
        let pi = hg.osmotic_pressure(0.05);
        assert!(pi.is_finite(), "pi={pi}");
    }
    #[test]
    fn thermochromic_color_at_low_temp() {
        let tc = ThermochromicResponse::new(340.0, 5.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e4);
        let [r, g, _b] = tc.reflectance(250.0);
        assert!(r > 0.9, "should be mostly red at low temp: r={r}");
        assert!(g < 0.1, "g={g}");
    }
    #[test]
    fn thermochromic_color_at_high_temp() {
        let tc = ThermochromicResponse::new(340.0, 5.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e4);
        let [r, g, _b] = tc.reflectance(450.0);
        assert!(g > 0.9, "should be mostly green at high temp: g={g}");
        assert!(r < 0.1, "r={r}");
    }
    #[test]
    fn thermochromic_transition_fraction_half_at_t_transition() {
        let tc = ThermochromicResponse::new(340.0, 5.0, [0.0; 3], [1.0; 3], 1e4);
        let s = tc.transition_fraction(340.0);
        assert!((s - 0.5).abs() < 1e-10, "s={s}");
    }
    #[test]
    fn thermochromic_effective_heat_capacity_increased_at_transition() {
        let tc = ThermochromicResponse::new(340.0, 5.0, [0.0; 3], [1.0; 3], 1e4);
        let c_base = 500.0;
        let c_far = tc.effective_heat_capacity(250.0, c_base);
        let c_at = tc.effective_heat_capacity(340.0, c_base);
        assert!(c_at > c_far, "c_at={c_at} > c_far={c_far}");
    }
    #[test]
    fn magnetic_shape_initial_strain_zero() {
        let m = MagneticShape::ni2mnga();
        assert!(m.strain().abs() < 1e-12);
    }
    #[test]
    fn magnetic_shape_above_critical_full_switch() {
        let mut m = MagneticShape::ni2mnga();
        m.update(m.h_critical * 2.0);
        assert!((m.eta - 1.0).abs() < 1e-10, "eta={}", m.eta);
    }
    #[test]
    fn magnetic_shape_strain_increases_with_field() {
        let mut m = MagneticShape::ni2mnga();
        m.update(m.h_critical * 0.5);
        let s1 = m.strain();
        m.update(m.h_critical * 0.9);
        let s2 = m.strain();
        assert!(s2 >= s1, "s2={s2}, s1={s1}");
    }
    #[test]
    fn magnetic_shape_blocking_stress_positive() {
        let mut m = MagneticShape::ni2mnga();
        m.update(m.h_critical * 2.0);
        let sigma = m.blocking_stress();
        assert!(sigma > 0.0, "sigma={sigma}");
    }
    #[test]
    fn ferroelectric_initial_polarization_negative() {
        let fe = FerroelectricHysteresis::barium_titanate();
        assert!(fe.polarization < 0.0, "p={}", fe.polarization);
    }
    #[test]
    fn ferroelectric_positive_field_switches_positive() {
        let mut fe = FerroelectricHysteresis::barium_titanate();
        fe.update(fe.e_coercive * 5.0);
        assert!(fe.polarization > 0.0, "p={}", fe.polarization);
    }
    #[test]
    fn ferroelectric_hysteresis_loop_area_positive() {
        let mut fe = FerroelectricHysteresis::barium_titanate();
        let e_max = fe.e_coercive * 3.0;
        let mut work = 0.0_f64;
        let n = 100;
        let e_values: Vec<f64> = (0..=n)
            .map(|i| e_max * (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let mut p_prev = fe.polarization;
        for &e in e_values.iter().skip(1) {
            fe.update(e);
            let de = e - (e_values[0]);
            work += (fe.polarization + p_prev) / 2.0 * de;
            p_prev = fe.polarization;
        }
        assert!(work.is_finite(), "work={work}");
    }
    #[test]
    fn ferroelectric_susceptibility_finite() {
        let fe = FerroelectricHysteresis::barium_titanate();
        let chi = fe.susceptibility(0.0, 1e3);
        assert!(chi.is_finite(), "chi={chi}");
    }
    #[test]
    fn smart_composite_longitudinal_modulus_between_constituents() {
        let sma = ShapeMemoryAlloy::nitinol();
        let e_sma = sma.elastic_modulus();
        let e_mat = 70e9_f64;
        let sc = SmartComposite::new(0.3, sma, e_mat, 2700.0, 6450.0, 1e-3);
        let e_l = sc.longitudinal_modulus();
        assert!(
            e_l >= e_mat.min(e_sma) && e_l <= e_mat.max(e_sma),
            "e_l={e_l} should be between {e_mat} and {e_sma}"
        );
    }
    #[test]
    fn smart_composite_transverse_modulus_below_longitudinal() {
        let sma = ShapeMemoryAlloy::nitinol();
        let sc = SmartComposite::new(0.5, sma, 70e9, 2700.0, 6450.0, 1e-3);
        let e_l = sc.longitudinal_modulus();
        let e_t = sc.transverse_modulus();
        assert!(e_t <= e_l + 1.0, "e_t={e_t}, e_l={e_l}");
    }
    #[test]
    fn smart_composite_density_between_constituents() {
        let sma = ShapeMemoryAlloy::nitinol();
        let rho_m = 2700.0_f64;
        let rho_sma = 6450.0_f64;
        let sc = SmartComposite::new(0.4, sma, 70e9, rho_m, rho_sma, 1e-3);
        let rho = sc.density();
        assert!(rho >= rho_m && rho <= rho_sma, "rho={rho}");
    }
    #[test]
    fn smart_composite_recovery_stress_finite() {
        let sma = ShapeMemoryAlloy::nitinol();
        let sc = SmartComposite::new(0.3, sma, 70e9, 2700.0, 6450.0, 1e-3);
        let sigma = sc.composite_recovery_stress(0.02);
        assert!(sigma.is_finite(), "sigma={sigma}");
    }
    #[test]
    fn smart_composite_curvature_finite() {
        let sma = ShapeMemoryAlloy::nitinol();
        let sc = SmartComposite::new(0.3, sma, 70e9, 2700.0, 6450.0, 1e-3);
        let kappa = sc.actuation_curvature(350.0);
        assert!(kappa.is_finite(), "kappa={kappa}");
    }
}
#[cfg(test)]
mod smart_ext_tests {

    use crate::smart_materials::BrinsonModel;

    use crate::smart_materials::Electrostrictive;

    use crate::smart_materials::HydrogelFloryRehner;

    use crate::smart_materials::Magnetorheological;

    use crate::smart_materials::PiezoelectricSmart;

    use crate::smart_materials::SelfHealingMaterial;

    #[test]
    fn brinson_initial_fully_martensite() {
        let b = BrinsonModel::nitinol();
        assert!((b.total_xi() - 1.0).abs() < 1e-10, "xi={}", b.total_xi());
    }
    #[test]
    fn brinson_heated_above_af_fully_austenite() {
        let mut b = BrinsonModel::nitinol();
        b.update(400.0, 0.0);
        assert!(b.total_xi() < 0.01, "xi={}", b.total_xi());
    }
    #[test]
    fn brinson_elastic_modulus_between_phases() {
        let mut b = BrinsonModel::nitinol();
        let e_martensite = b.elastic_modulus();
        b.update(400.0, 0.0);
        let e_austenite = b.elastic_modulus();
        assert!(
            e_austenite > e_martensite,
            "E_A={e_austenite} should be > E_M={e_martensite}"
        );
    }
    #[test]
    fn brinson_recovery_stress_finite() {
        let b = BrinsonModel::nitinol();
        let sigma = b.recovery_stress(0.04);
        assert!(sigma.is_finite(), "sigma={sigma}");
    }
    #[test]
    fn brinson_stress_from_strain_finite() {
        let b = BrinsonModel::nitinol();
        let sigma = b.stress_from_strain(0.02, 298.0);
        assert!(sigma.is_finite(), "sigma={sigma}");
    }
    #[test]
    fn brinson_sme_stroke_positive() {
        let b = BrinsonModel::nitinol();
        let stroke = b.sme_stroke(1.0, 0.0, 0.1);
        assert!(stroke > 0.0, "stroke={stroke}");
    }
    #[test]
    fn brinson_xi_t_increases_on_cooling() {
        let mut b = BrinsonModel::nitinol();
        b.update(400.0, 0.0);
        let xi_before = b.total_xi();
        b.update(280.0, 0.0);
        let xi_after = b.total_xi();
        assert!(xi_after >= xi_before, "xi should increase on cooling");
    }
    #[test]
    fn mr_new_yield_stress_increases_with_field() {
        let mr = Magnetorheological::carbonyl_iron();
        let t1 = mr.yield_stress(1e4);
        let t2 = mr.yield_stress(1e5);
        assert!(t2 > t1, "yield stress should increase with field");
    }
    #[test]
    fn mr_new_shear_stress_finite() {
        let mr = Magnetorheological::carbonyl_iron();
        let tau = mr.shear_stress(100.0, 1e4);
        assert!(tau.is_finite() && tau > 0.0, "tau={tau}");
    }
    #[test]
    fn mr_new_mason_number_infinite_zero_field() {
        let mr = Magnetorheological::carbonyl_iron();
        let mn = mr.mason_number(100.0, 0.0);
        assert!(mn.is_infinite(), "mn should be infinite at H=0");
    }
    #[test]
    fn mr_new_apparent_viscosity_increases_with_field() {
        let mr = Magnetorheological::carbonyl_iron();
        let v1 = mr.apparent_viscosity(100.0, 0.0);
        let v2 = mr.apparent_viscosity(100.0, 1e4);
        assert!(v2 > v1, "apparent viscosity should increase with field");
    }
    #[test]
    fn mr_new_chain_fraction_zero_at_zero_field() {
        let mr = Magnetorheological::carbonyl_iron();
        let cf = mr.chain_fraction(0.0);
        assert!(cf.abs() < 1e-6, "cf={cf}");
    }
    #[test]
    fn mr_new_chain_fraction_increases_with_field() {
        let mr = Magnetorheological::carbonyl_iron();
        let cf1 = mr.chain_fraction(1e3);
        let cf2 = mr.chain_fraction(1e5);
        assert!(cf2 >= cf1, "chain fraction should not decrease with field");
    }
    #[test]
    fn mr_effective_viscosity_off_ge_eta0() {
        let mr = Magnetorheological::carbonyl_iron();
        let eta_eff = mr.effective_viscosity_off();
        assert!(eta_eff >= mr.eta0, "effective viscosity >= base viscosity");
    }
    #[test]
    fn electrostrictive_polarization_increases_with_field() {
        let es = Electrostrictive::pmn_pt();
        let p1 = es.polarization(1e5);
        let p2 = es.polarization(1e6);
        assert!(p2 > p1, "p2={p2}, p1={p1}");
    }
    #[test]
    fn electrostrictive_maxwell_stress_positive() {
        let es = Electrostrictive::pmn_pt();
        let p = es.maxwell_stress(1e6);
        assert!(p > 0.0, "p={p}");
    }
    #[test]
    fn electrostrictive_strain_33_nonnegative() {
        let es = Electrostrictive::pmn_pt();
        let eps = es.strain_33(1e6);
        assert!(eps >= 0.0, "eps={eps}");
    }
    #[test]
    fn electrostrictive_strain_31_nonpositive() {
        let es = Electrostrictive::pmn_pt();
        let eps = es.strain_31(1e6);
        assert!(eps <= 0.0, "eps={eps}");
    }
    #[test]
    fn electrostrictive_thickness_change_finite() {
        let es = Electrostrictive::pmn_pt();
        let dt = es.thickness_change(100.0);
        assert!(dt.is_finite(), "dt={dt}");
    }
    #[test]
    fn electrostrictive_blocking_stress_positive() {
        let es = Electrostrictive::pmn_pt();
        let sigma = es.blocking_stress(100.0);
        assert!(sigma >= 0.0, "sigma={sigma}");
    }
    #[test]
    fn electrostrictive_coupling_in_unit_interval() {
        let es = Electrostrictive::pmn_pt();
        let k = es.coupling_coefficient(1e6);
        assert!((0.0..=1.0).contains(&k), "k={k}");
    }
    #[test]
    fn piezosmart_k33_in_unit_interval() {
        let p = PiezoelectricSmart::pzt5h();
        let k = p.k33();
        assert!(k > 0.0 && k <= 1.0, "k33={k}");
    }
    #[test]
    fn piezosmart_kp_in_unit_interval() {
        let p = PiezoelectricSmart::pzt4();
        let kp = p.kp();
        assert!((0.0..=1.0).contains(&kp), "kp={kp}");
    }
    #[test]
    fn piezosmart_resonance_frequency_positive() {
        let p = PiezoelectricSmart::pzt5h();
        let fr = p.resonance_frequency();
        assert!(fr > 0.0, "fr={fr}");
    }
    #[test]
    fn piezosmart_antiresonance_above_resonance() {
        let p = PiezoelectricSmart::pzt5h();
        let fr = p.resonance_frequency();
        let fa = p.antiresonance_frequency();
        assert!(fa > fr, "fa={fa}, fr={fr}");
    }
    #[test]
    fn piezosmart_actuation_stroke_linear() {
        let p = PiezoelectricSmart::pzt5h();
        let s1 = p.actuation_stroke(100.0);
        let s2 = p.actuation_stroke(200.0);
        assert!((s2 / s1 - 2.0).abs() < 1e-10, "ratio={}", s2 / s1);
    }
    #[test]
    fn piezosmart_blocking_force_positive() {
        let p = PiezoelectricSmart::pzt5h();
        let f = p.blocking_force(100.0);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn piezosmart_charge_from_force_nonzero() {
        let p = PiezoelectricSmart::pzt5h();
        let q = p.charge_from_force(10.0);
        assert!(q > 0.0, "q={q}");
    }
    #[test]
    fn piezosmart_voltage_from_force_finite() {
        let p = PiezoelectricSmart::pzt5h();
        let v = p.voltage_from_force(1.0);
        assert!(v.is_finite(), "v={v}");
    }
    #[test]
    fn piezosmart_mechanical_q_positive() {
        let p = PiezoelectricSmart::pzt5h();
        let q = p.mechanical_q();
        assert!(q > 0.0, "q={q}");
    }
    #[test]
    fn hydrogel_fr_mixing_fe_finite() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let f = hg.mixing_free_energy(0.1);
        assert!(f.is_finite(), "f={f}");
    }
    #[test]
    fn hydrogel_fr_elastic_fe_positive() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let f = hg.elastic_free_energy(0.05);
        let f0 = hg.elastic_free_energy(hg.phi0);
        assert!(f0.abs() < 1e-6, "f_el at phi0 should be ~0, got {f0}");
        let _ = f;
    }
    #[test]
    fn hydrogel_fr_swelling_ratio_above_one() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let q = hg.equilibrium_swelling();
        assert!(q > 1.0, "q={q}");
    }
    #[test]
    fn hydrogel_fr_linear_strain_positive() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let eps = hg.linear_strain();
        assert!(eps > 0.0, "eps={eps}");
    }
    #[test]
    fn hydrogel_fr_osmotic_pressure_finite() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let pi = hg.osmotic_pressure(0.1);
        assert!(pi.is_finite(), "pi={pi}");
    }
    #[test]
    fn hydrogel_fr_chemical_potential_finite() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let mu = hg.chemical_potential(0.1);
        assert!(mu.is_finite(), "mu={mu}");
    }
    #[test]
    fn hydrogel_fr_total_fe_finite() {
        let hg = HydrogelFloryRehner::polyacrylamide();
        let f = hg.total_free_energy(0.1);
        assert!(f.is_finite(), "f={f}");
    }
    #[test]
    fn self_heal_initial_strength_is_virgin() {
        let m = SelfHealingMaterial::microencapsulated_epoxy();
        assert!((m.current_strength() - m.virgin_strength).abs() < 1e-6);
    }
    #[test]
    fn self_heal_damage_reduces_strength() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.5);
        let s = m.current_strength();
        assert!(s < m.virgin_strength, "s={s}");
    }
    #[test]
    fn self_heal_healing_increases_strength() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.3);
        let s_before = m.current_strength();
        for _ in 0..3600 {
            m.heal(1.0, 333.0);
        }
        let s_after = m.current_strength();
        assert!(
            s_after >= s_before,
            "s_after={s_after}, s_before={s_before}"
        );
    }
    #[test]
    fn self_heal_healing_efficiency_bounded() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.5);
        for _ in 0..100000 {
            m.heal(1.0, 333.0);
        }
        assert!(m.healing_efficiency <= m.max_efficiency + 1e-10);
    }
    #[test]
    fn self_heal_stiffness_reduction_range() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.4);
        let r = m.stiffness_reduction();
        assert!((0.0..=1.0).contains(&r), "r={r}");
    }
    #[test]
    fn self_heal_toughness_recovery_le_one() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.5);
        for _ in 0..10000 {
            m.heal(1.0, 333.0);
        }
        let k = m.toughness_recovery();
        assert!((0.0..=1.0 + 1e-10).contains(&k), "k={k}");
    }
    #[test]
    fn self_heal_healing_rate_increases_with_temperature() {
        let m = SelfHealingMaterial::microencapsulated_epoxy();
        let k1 = m.healing_rate(298.0);
        let k2 = m.healing_rate(333.0);
        assert!(k2 > k1, "k2={k2}, k1={k1}");
    }
    #[test]
    fn self_heal_no_healing_below_threshold() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.damage = 0.01;
        let eta_before = m.healing_efficiency;
        m.heal(1000.0, 333.0);
        assert!((m.healing_efficiency - eta_before).abs() < 1e-10);
    }
    #[test]
    fn self_heal_progress_monotone() {
        let mut m = SelfHealingMaterial::microencapsulated_epoxy();
        m.apply_damage(0.5);
        let mut prev = m.healing_progress();
        for _ in 0..100 {
            m.heal(10.0, 333.0);
            let prog = m.healing_progress();
            assert!(prog >= prev - 1e-12, "progress decreased: {prog} < {prev}");
            prev = prog;
        }
    }
}
