//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod piezo_extended_tests2 {
    use super::*;
    use crate::piezo::*;
    #[test]
    fn test_pvdf_capacitance_positive() {
        let pvdf = PvdfFilm::typical();
        let c = pvdf.capacitance();
        assert!(c > 0.0, "PVDF capacitance should be positive: {c}");
        assert!(c.is_finite());
    }
    #[test]
    fn test_pvdf_voltage_from_stress_d31_nonzero() {
        let pvdf = PvdfFilm::typical();
        let v = pvdf.voltage_from_stress_d31(1e6);
        assert!(
            v.abs() > 0.0,
            "voltage from stress (d31) should be non-zero"
        );
        assert!(v.is_finite());
    }
    #[test]
    fn test_pvdf_voltage_from_stress_d33_opposite_sign() {
        let pvdf = PvdfFilm::typical();
        let v31 = pvdf.voltage_from_stress_d31(1e6);
        let v33 = pvdf.voltage_from_stress_d33(1e6);
        assert!(
            v31 * v33 < 0.0,
            "d31 and d33 voltages should have opposite signs for PVDF"
        );
    }
    #[test]
    fn test_pvdf_charge_from_force_d31_nonzero() {
        let pvdf = PvdfFilm::typical();
        let q = pvdf.charge_from_force_d31(1.0);
        assert!(q.abs() > 0.0, "charge from force should be non-zero");
    }
    #[test]
    fn test_pvdf_actuation_strain_d31_sign() {
        let pvdf = PvdfFilm::typical();
        let eps = pvdf.actuation_strain_d31(100.0);
        assert!(
            eps > 0.0,
            "d31 actuation strain should be positive for positive V"
        );
    }
    #[test]
    fn test_pvdf_actuation_strain_d33_sign() {
        let pvdf = PvdfFilm::typical();
        let eps = pvdf.actuation_strain_d33(100.0);
        assert!(
            eps < 0.0,
            "d33 actuation strain should be negative for PVDF"
        );
    }
    #[test]
    fn test_pvdf_k31_coupling_factor_range() {
        let pvdf = PvdfFilm::typical();
        let k31 = pvdf.k31_coupling_factor();
        assert!(
            k31 > 0.0 && k31 < 0.5,
            "PVDF k31 = {k31}, expected ~ 0.1-0.15"
        );
    }
    #[test]
    fn test_pvdf_bending_stiffness_positive() {
        let pvdf = PvdfFilm::typical();
        let bs = pvdf.bending_stiffness_per_unit_width();
        assert!(bs > 0.0, "bending stiffness should be positive");
        assert!(bs.is_finite());
    }
    #[test]
    fn test_harvester_max_efficiency_range() {
        let h = PiezoHarvester::new(0.1, 100.0, 2.0 * std::f64::consts::PI * 100.0, 10e3);
        let eta = h.max_efficiency();
        assert!(
            eta > 0.0 && eta < 1.0,
            "efficiency should be in (0,1): {eta}"
        );
    }
    #[test]
    fn test_harvester_efficiency_increases_with_coupling() {
        let h1 = PiezoHarvester::new(0.01, 100.0, 100.0, 1e4);
        let h2 = PiezoHarvester::new(0.5, 100.0, 100.0, 1e4);
        assert!(
            h2.max_efficiency() > h1.max_efficiency(),
            "higher k^2 → higher efficiency"
        );
    }
    #[test]
    fn test_harvester_max_power_at_resonance_positive() {
        let h = PiezoHarvester::new(0.1, 100.0, 2.0 * std::f64::consts::PI * 100.0, 10e3);
        let p = h.max_power_at_resonance(0.01, 9.81);
        assert!(p > 0.0, "power should be positive: {p}");
    }
    #[test]
    fn test_harvester_figure_of_merit() {
        let h = PiezoHarvester::new(0.2, 50.0, 100.0, 1e4);
        let fom = h.figure_of_merit();
        assert!((fom - 0.2 * 50.0).abs() < 1e-10, "FOM = k^2 * Q_m");
    }
    #[test]
    fn test_harvester_bandwidth_positive() {
        let h = PiezoHarvester::new(0.1, 100.0, 2.0 * std::f64::consts::PI * 100.0, 10e3);
        let bw = h.bandwidth_hz();
        assert!(
            bw > 0.0 && bw.is_finite(),
            "bandwidth should be finite positive: {bw}"
        );
    }
    #[test]
    fn test_piezo_beam_resonance_positive() {
        let f = piezo_beam_resonance_frequency(200e9, 1e-8, 7800.0, 1e-4, 0.1);
        assert!(
            f > 0.0 && f.is_finite(),
            "resonance frequency should be positive finite: {f}"
        );
    }
    #[test]
    fn test_piezo_beam_resonance_increases_with_stiffness() {
        let f1 = piezo_beam_resonance_frequency(200e9, 1e-8, 7800.0, 1e-4, 0.1);
        let f2 = piezo_beam_resonance_frequency(400e9, 1e-8, 7800.0, 1e-4, 0.1);
        assert!(
            f2 > f1,
            "stiffer beam should have higher resonance frequency"
        );
    }
    #[test]
    fn test_piezo_beam_resonance_decreases_with_length() {
        let f1 = piezo_beam_resonance_frequency(200e9, 1e-8, 7800.0, 1e-4, 0.1);
        let f2 = piezo_beam_resonance_frequency(200e9, 1e-8, 7800.0, 1e-4, 0.2);
        assert!(f2 < f1, "longer beam should have lower resonance frequency");
    }
    #[test]
    fn test_piezo_beam_mode_n_increasing() {
        let params = (200e9, 1e-8, 7800.0, 1e-4, 0.1);
        let f1 = piezo_beam_resonance_mode_n(params.0, params.1, params.2, params.3, params.4, 1);
        let f2 = piezo_beam_resonance_mode_n(params.0, params.1, params.2, params.3, params.4, 2);
        let f3 = piezo_beam_resonance_mode_n(params.0, params.1, params.2, params.3, params.4, 3);
        assert!(
            f2 > f1 && f3 > f2,
            "mode frequencies should increase: f1={f1}, f2={f2}, f3={f3}"
        );
    }
    #[test]
    fn test_piezo_actuator_force_d33_positive_voltage() {
        let f = piezo_actuator_force_d33(100.0, 374e-12, 111e9, 1e-4, 1e-3);
        assert!(
            f > 0.0,
            "d33 actuator force should be positive for positive voltage and d33"
        );
    }
    #[test]
    fn test_piezo_actuator_force_d31_negative_for_negative_d31() {
        let f = piezo_actuator_force_d31(100.0, -171e-12, 121e9, 1e-4, 1e-3);
        assert!(f < 0.0, "d31 force should be negative for negative d31");
    }
    #[test]
    fn test_piezo_sensor_voltage_d33_nonzero() {
        let v = piezo_sensor_voltage_d33(1e-6, 374e-12, 111e9, 8.854e-12 * 1470.0, 1e-3);
        assert!(v.abs() > 0.0, "sensor voltage d33 should be non-zero");
    }
    #[test]
    fn test_piezo_sensor_voltage_d31_nonzero() {
        let v = piezo_sensor_voltage_d31(1e-6, -171e-12, 121e9, 8.854e-12 * 1700.0, 1e-3);
        assert!(v.abs() > 0.0, "sensor voltage d31 should be non-zero");
    }
    #[test]
    fn test_piezo_sensor_voltage_scales_linearly() {
        let strain = 1e-6;
        let v1 = piezo_sensor_voltage_d33(strain, 374e-12, 111e9, 8.854e-12 * 1470.0, 1e-3);
        let v2 = piezo_sensor_voltage_d33(2.0 * strain, 374e-12, 111e9, 8.854e-12 * 1470.0, 1e-3);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-10,
            "sensor output should scale linearly with strain"
        );
    }
    #[test]
    fn test_pvdf_permittivity_nonzero() {
        let pvdf = PvdfFilm::typical();
        let eps = pvdf.permittivity();
        assert!(eps > 0.0 && eps.is_finite(), "PVDF permittivity = {eps}");
    }
    #[test]
    fn test_harvester_zero_mass_zero_power() {
        let h = PiezoHarvester::new(0.1, 100.0, 100.0, 1e4);
        let p = h.max_power_at_resonance(0.0, 9.81);
        assert!(p.abs() < 1e-30, "zero mass should give zero power");
    }
    #[test]
    fn test_electromechanical_coupling_pzt5a_positive() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let k33 = elem.compute_electromechanical_coupling();
        assert!(
            k33 > 0.0 && k33 < 1.0,
            "k33 for PZT-5A should be in (0, 1): {}",
            k33
        );
    }
    #[test]
    fn test_electromechanical_coupling_zero_permittivity() {
        let mut mat = PiezoMaterial::pzt5a();
        mat.epsilon_matrix[2][2] = 0.0;
        let elem = PiezoElement::new(mat, 1e-9);
        let k33 = elem.compute_electromechanical_coupling();
        assert_eq!(k33, 0.0, "Zero permittivity should give k33 = 0");
    }
    #[test]
    fn test_resonance_frequency_pzt5a_audible_range() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let rho = 7750.0;
        let l = 1e-3;
        let f = elem.compute_resonance_frequency(l, rho);
        assert!(
            f > 1e4 && f < 1e8,
            "Resonance frequency for 1mm PZT-5A: {} Hz (expect ~100 kHz)",
            f
        );
    }
    #[test]
    fn test_resonance_frequency_scales_inverse_length() {
        let mat1 = PiezoMaterial::pzt5a();
        let mat2 = PiezoMaterial::pzt5a();
        let elem1 = PiezoElement::new(mat1, 1e-9);
        let elem2 = PiezoElement::new(mat2, 1e-9);
        let rho = 7750.0;
        let f1 = elem1.compute_resonance_frequency(1e-3, rho);
        let f2 = elem2.compute_resonance_frequency(2e-3, rho);
        assert!(
            (f1 / f2 - 2.0).abs() < 1e-10,
            "Frequency should halve when length doubles: f1={}, f2={}",
            f1,
            f2
        );
    }
    #[test]
    fn test_resonance_frequency_zero_length_returns_zero() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let f = elem.compute_resonance_frequency(0.0, 7750.0);
        assert_eq!(f, 0.0, "Zero length should return zero frequency");
    }
    #[test]
    fn test_effective_piezo_charge_pzt5a_nonzero() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let d33 = elem.compute_effective_piezo_charge();
        assert!(
            d33.abs() > 0.0 && d33.is_finite(),
            "d33 should be nonzero finite: {}",
            d33
        );
    }
    #[test]
    fn test_effective_piezo_charge_zero_stiffness() {
        let mut mat = PiezoMaterial::pzt5a();
        mat.c_matrix[2][2] = 0.0;
        let elem = PiezoElement::new(mat, 1e-9);
        let d33 = elem.compute_effective_piezo_charge();
        assert_eq!(d33, 0.0, "Zero stiffness should give d33 = 0");
    }
    #[test]
    fn test_blocked_force_proportional_to_voltage() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let f1 = elem.compute_blocked_force(10.0, 1e-4, 1e-3);
        let f2 = elem.compute_blocked_force(20.0, 1e-4, 1e-3);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-10,
            "Blocked force should double when voltage doubles: f1={}, f2={}",
            f1,
            f2
        );
    }
    #[test]
    fn test_free_displacement_proportional_to_voltage() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let d1 = elem.compute_free_displacement(100.0);
        let d2 = elem.compute_free_displacement(200.0);
        assert!(
            (d2 / d1 - 2.0).abs() < 1e-10,
            "Free displacement should double when voltage doubles: d1={}, d2={}",
            d1,
            d2
        );
    }
    #[test]
    fn test_coupling_consistency_with_charge_coefficient() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let k33 = elem.compute_electromechanical_coupling();
        let d33 = elem.compute_effective_piezo_charge();
        let c33 = elem.material.c_matrix[2][2];
        let eps33 = elem.material.epsilon_matrix[2][2];
        let k33_from_d33 = (d33 * d33 * c33 / eps33).sqrt();
        assert!(
            (k33 - k33_from_d33).abs() / k33 < 1e-9,
            "k33 from e33 vs d33: {} vs {}",
            k33,
            k33_from_d33
        );
    }
}
