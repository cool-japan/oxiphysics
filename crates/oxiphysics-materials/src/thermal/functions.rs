//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Gas constant R (J/(mol*K)).
pub(super) const R_GAS: f64 = 8.314;
/// Fourier number: Fo = alpha * t / L^2.
///
/// Dimensionless time for heat conduction problems.
pub fn fourier_number(diffusivity: f64, time: f64, length: f64) -> f64 {
    diffusivity * time / (length * length)
}
/// Semi-infinite solid surface temperature response to a step change in surface temperature.
///
/// T(x,t) = T_s + (T_i - T_s) * erf(x / (2 * sqrt(alpha * t)))
pub fn semi_infinite_temperature(
    x: f64,
    t: f64,
    diffusivity: f64,
    t_initial: f64,
    t_surface: f64,
) -> f64 {
    if t < 1e-30 {
        return t_initial;
    }
    let eta = x / (2.0 * (diffusivity * t).sqrt());
    t_surface + (t_initial - t_surface) * erf_approx(eta)
}
/// Simple erf approximation for thermal calculations.
pub(super) fn erf_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = 1.0 - poly * (-x * x).exp();
    if x >= 0.0 { result } else { -result }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermal::AblationModel;
    use crate::thermal::DebyeModel;
    use crate::thermal::EinsteinModel;
    use crate::thermal::EnthalpyMethod;
    use crate::thermal::HeatAffectedZone;
    use crate::thermal::HeatConduction1D;
    use crate::thermal::JohnsonCookModel;
    use crate::thermal::NewtonianCooling;
    use crate::thermal::PhaseChangeModel;
    use crate::thermal::TemperatureDependentConductivity;
    use crate::thermal::TemperatureDependentSpecificHeat;
    use crate::thermal::ThermalConductivityTensor;
    use crate::thermal::ThermalExpansion;
    use crate::thermal::ThermalFatigue;
    use crate::thermal::ThermalInterfaceResistance;
    use crate::thermal::ThermalMaterial;
    use crate::thermal::ThermalShockParam;
    use crate::thermal::ThermalShockResistance;
    use crate::thermal::ThermalStressCoupling;
    #[test]
    fn test_steel_diffusivity() {
        let mat = ThermalMaterial::steel();
        let alpha = mat.diffusivity();
        assert!(
            (alpha - 1.2739e-5).abs() < 1e-8,
            "steel diffusivity = {alpha}"
        );
    }
    #[test]
    fn test_aluminum_effusivity() {
        let mat = ThermalMaterial::aluminum();
        let e = mat.effusivity();
        let expected = (237.0_f64 * 2700.0 * 900.0).sqrt();
        assert!((e - expected).abs() < 1.0, "aluminium effusivity = {e}");
    }
    #[test]
    fn test_thermal_strain() {
        let mat = ThermalMaterial::steel();
        let strain = mat.thermal_strain(100.0);
        assert!((strain - 1.2e-3).abs() < 1e-12);
    }
    #[test]
    fn test_volumetric_heat_capacity() {
        let mat = ThermalMaterial::steel();
        let rhoc = mat.volumetric_heat_capacity();
        let expected = 7850.0 * 500.0;
        assert!((rhoc - expected).abs() < 1.0);
    }
    #[test]
    fn test_penetration_depth() {
        let mat = ThermalMaterial::steel();
        let d = mat.penetration_depth(1.0);
        assert!(d > 0.0);
        assert!(d.is_finite());
        let d2 = mat.penetration_depth(4.0);
        assert!((d2 / d - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_temp_dependent_conductivity_interpolation() {
        let tdc = TemperatureDependentConductivity::new(
            vec![300.0, 400.0, 500.0],
            vec![50.0, 45.0, 40.0],
        );
        let k = tdc.conductivity_at(350.0);
        assert!((k - 47.5).abs() < 1e-10);
        assert!((tdc.conductivity_at(200.0) - 50.0).abs() < 1e-10);
        assert!((tdc.conductivity_at(600.0) - 40.0).abs() < 1e-10);
    }
    #[test]
    fn test_temp_dependent_specific_heat() {
        let tdcp = TemperatureDependentSpecificHeat::new(vec![300.0, 500.0], vec![500.0, 600.0]);
        let cp = tdcp.cp_at(400.0);
        assert!((cp - 550.0).abs() < 1e-10);
    }
    #[test]
    fn test_debye_high_temperature_limit() {
        let debye = DebyeModel::new(300.0, 1.0);
        let cv = debye.molar_cv(10000.0);
        let dp = debye.dulong_petit();
        assert!(
            (cv - dp).abs() / dp < 0.01,
            "high-T Cv = {cv}, expected ~{dp}"
        );
    }
    #[test]
    fn test_debye_low_temperature() {
        let debye = DebyeModel::new(300.0, 1.0);
        let cv_low = debye.molar_cv(10.0);
        let cv_high = debye.molar_cv(1000.0);
        assert!(cv_low < cv_high, "low-T Cv should be smaller than high-T");
    }
    #[test]
    fn test_debye_zero_temperature() {
        let debye = DebyeModel::new(300.0, 1.0);
        let cv = debye.molar_cv(0.0);
        assert!(cv.abs() < 1e-12, "Cv at T=0 should be 0");
    }
    #[test]
    fn test_einstein_high_temperature_limit() {
        let einstein = EinsteinModel::new(200.0, 1.0);
        let cv = einstein.molar_cv(10000.0);
        let dp = 3.0 * 1.0 * R_GAS;
        assert!(
            (cv - dp).abs() / dp < 0.01,
            "high-T Cv = {cv}, expected ~{dp}"
        );
    }
    #[test]
    fn test_einstein_low_temperature() {
        let einstein = EinsteinModel::new(200.0, 1.0);
        let cv_10 = einstein.molar_cv(10.0);
        let cv_1000 = einstein.molar_cv(1000.0);
        assert!(cv_10 < cv_1000);
    }
    #[test]
    fn test_thermal_shock_r_parameter() {
        let tsr = ThermalShockResistance::new(400e6, 410e9, 0.14, 4e-6, 120.0);
        let r = tsr.r_parameter();
        assert!(r > 0.0, "R should be positive, got {r}");
        assert!(r > 100.0 && r < 1000.0, "R = {r} K, expected ~200");
    }
    #[test]
    fn test_thermal_shock_will_fail() {
        let tsr = ThermalShockResistance::new(400e6, 410e9, 0.14, 4e-6, 120.0);
        let r = tsr.r_parameter();
        assert!(!tsr.will_fail(r * 0.5));
        assert!(tsr.will_fail(r * 1.5));
    }
    #[test]
    fn test_thermal_shock_r_prime() {
        let tsr = ThermalShockResistance::new(400e6, 410e9, 0.14, 4e-6, 120.0);
        let rp = tsr.r_prime();
        assert!(rp > 0.0);
        assert!((rp - tsr.r_parameter() * 120.0).abs() < 1.0);
    }
    #[test]
    fn test_thermal_interface_heat_flux() {
        let tir = ThermalInterfaceResistance::metal_metal();
        let q = tir.heat_flux(10.0);
        assert!((q - 1e5).abs() < 1.0);
    }
    #[test]
    fn test_thermal_interface_roundtrip() {
        let tir = ThermalInterfaceResistance::new(5e-4);
        let delta_t = 25.0;
        let q = tir.heat_flux(delta_t);
        let dt_back = tir.temperature_drop(q);
        assert!((dt_back - delta_t).abs() < 1e-10);
    }
    #[test]
    fn test_heat_conduction_1d_steady() {
        let mat = ThermalMaterial::steel();
        let mut solver = HeatConduction1D::new(21, 1.0, 300.0, mat);
        solver.set_temperature_bc(300.0, 500.0);
        let dt = solver.critical_dt() * 0.4;
        for _ in 0..100_000 {
            solver.step_explicit(dt);
            solver.set_temperature_bc(300.0, 500.0);
        }
        let ss = solver.steady_state_temperature(300.0, 500.0);
        for (i, (&t_actual, &t_expected)) in solver.temperature.iter().zip(ss.iter()).enumerate() {
            assert!(
                (t_actual - t_expected).abs() < 0.5,
                "node {i}: T={t_actual}, expected ~{t_expected}"
            );
        }
    }
    #[test]
    fn test_heat_conduction_critical_dt() {
        let mat = ThermalMaterial::copper();
        let solver = HeatConduction1D::new(11, 0.1, 300.0, mat);
        let dt = solver.critical_dt();
        assert!(dt > 0.0);
    }
    #[test]
    fn test_heat_conduction_with_source() {
        let mat = ThermalMaterial::steel();
        let mut solver = HeatConduction1D::new(11, 1.0, 300.0, mat);
        solver.set_temperature_bc(300.0, 300.0);
        let dt = solver.critical_dt() * 0.4;
        for _ in 0..100 {
            solver.step_explicit_with_source(dt, 1e6);
            solver.set_temperature_bc(300.0, 300.0);
        }
        let mid = solver.n_nodes / 2;
        assert!(solver.temperature[mid] > 300.0, "interior should heat up");
    }
    #[test]
    fn test_newtonian_cooling_time_constant() {
        let nc = NewtonianCooling::new(1.0, 500.0, 10.0, 0.1, 400.0);
        let tau = nc.time_constant();
        assert!((tau - 500.0).abs() < 1e-10);
    }
    #[test]
    fn test_newtonian_cooling_temperature() {
        let nc = NewtonianCooling::new(1.0, 500.0, 10.0, 0.1, 500.0);
        let tau = nc.time_constant();
        let t = nc.temperature_at_time(300.0, tau);
        let expected = 300.0 + (500.0 - 300.0) / std::f64::consts::E;
        assert!((t - expected).abs() < 1e-6);
    }
    #[test]
    fn test_newtonian_cooling_time_to_target() {
        let nc = NewtonianCooling::new(1.0, 500.0, 10.0, 0.1, 500.0);
        let tau = nc.time_constant();
        let t = nc.time_to_temperature(300.0, 400.0).unwrap();
        let expected = -tau * 0.5_f64.ln();
        assert!((t - expected).abs() < 1e-6);
    }
    #[test]
    fn test_thermal_stress_formula() {
        let mat = ThermalMaterial::steel();
        let sigma = mat.thermal_stress(100.0, 200e9, 0.3);
        let expected = -200e9 * 12e-6 * 100.0 / (1.0 - 0.3);
        assert!((sigma - expected).abs() < 1.0);
        assert!(sigma < 0.0);
    }
    #[test]
    fn test_fourier_number() {
        let fo = fourier_number(1e-5, 100.0, 0.1);
        assert!((fo - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_semi_infinite_temperature_at_surface() {
        let t = semi_infinite_temperature(0.0, 1.0, 1e-5, 300.0, 500.0);
        assert!((t - 500.0).abs() < 1e-6, "T at surface should be T_s = 500");
    }
    #[test]
    fn test_semi_infinite_temperature_far_away() {
        let t = semi_infinite_temperature(1.0, 0.001, 1e-5, 300.0, 500.0);
        assert!(
            (t - 300.0).abs() < 1.0,
            "T far from surface should be near T_i = 300, got {t}"
        );
    }
    #[test]
    fn test_biot_number() {
        let nc = NewtonianCooling::new(1.0, 500.0, 10.0, 0.1, 400.0);
        let bi = nc.biot_number(0.01, 50.0);
        assert!((bi - 0.002).abs() < 1e-10);
    }
    #[test]
    fn test_sic_higher_shock_resistance_than_alumina() {
        let sic = ThermalShockResistance::new(400e6, 410e9, 0.14, 4e-6, 120.0);
        let al2o3 = ThermalShockResistance::new(300e6, 380e9, 0.22, 8e-6, 30.0);
        assert!(
            sic.r_prime() > al2o3.r_prime(),
            "SiC should have higher R' than alumina"
        );
    }
    #[test]
    fn test_isotropic_tensor_heat_flux() {
        let kt = ThermalConductivityTensor::isotropic(50.0);
        let grad_t = [100.0, 0.0, 0.0];
        let q = kt.heat_flux(grad_t);
        assert!((q[0] - (-5000.0)).abs() < 1e-8, "q_x = {}", q[0]);
        assert!(q[1].abs() < 1e-12, "q_y should be 0");
        assert!(q[2].abs() < 1e-12, "q_z should be 0");
    }
    #[test]
    fn test_orthotropic_tensor_diagonal() {
        let kt = ThermalConductivityTensor::orthotropic(10.0, 20.0, 30.0);
        let grad_t = [1.0, 1.0, 1.0];
        let q = kt.heat_flux(grad_t);
        assert!((q[0] - (-10.0)).abs() < 1e-10);
        assert!((q[1] - (-20.0)).abs() < 1e-10);
        assert!((q[2] - (-30.0)).abs() < 1e-10);
    }
    #[test]
    fn test_thermal_tensor_symmetry() {
        let kt = ThermalConductivityTensor::isotropic(50.0);
        assert!(
            kt.is_symmetric(1e-10),
            "isotropic tensor should be symmetric"
        );
    }
    #[test]
    fn test_thermal_tensor_effective_conductivity_isotropic() {
        let k = 50.0;
        let kt = ThermalConductivityTensor::isotropic(k);
        let k_eff = kt.effective_conductivity();
        assert!((k_eff - k).abs() < 1e-8, "isotropic k_eff = k, got {k_eff}");
    }
    #[test]
    fn test_phase_change_liquid_fraction_solid() {
        let pc = PhaseChangeModel::new(500.0, 550.0, 300e3, 500.0, 550.0, 2700.0);
        assert!((pc.liquid_fraction(400.0) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_phase_change_liquid_fraction_liquid() {
        let pc = PhaseChangeModel::new(500.0, 550.0, 300e3, 500.0, 550.0, 2700.0);
        assert!((pc.liquid_fraction(600.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_phase_change_liquid_fraction_mushy() {
        let pc = PhaseChangeModel::new(500.0, 600.0, 300e3, 500.0, 550.0, 2700.0);
        let fl = pc.liquid_fraction(550.0);
        assert!((fl - 0.5).abs() < 1e-10, "at midpoint, fl = 0.5, got {fl}");
    }
    #[test]
    fn test_phase_change_apparent_cp_enhanced_in_mushy() {
        let pc = PhaseChangeModel::new(500.0, 600.0, 300e3, 500.0, 550.0, 2700.0);
        let cp_solid = pc.apparent_specific_heat(400.0);
        let cp_mushy = pc.apparent_specific_heat(550.0);
        assert!(
            cp_mushy > cp_solid,
            "apparent cp should be enhanced in mushy zone"
        );
    }
    #[test]
    fn test_phase_change_enthalpy_increases_with_temp() {
        let pc = PhaseChangeModel::new(500.0, 600.0, 300e3, 500.0, 550.0, 2700.0);
        let h1 = pc.enthalpy(400.0);
        let h2 = pc.enthalpy(700.0);
        assert!(h2 > h1, "enthalpy should increase with temperature");
    }
    #[test]
    fn test_phase_change_enthalpy_jump_across_melting() {
        let pc = PhaseChangeModel::new(500.0, 500.0, 300e3, 500.0, 550.0, 2700.0);
        let h_before = pc.enthalpy(499.9);
        let h_after = pc.enthalpy(500.1);
        let jump = h_after - h_before;
        assert!(
            jump > 0.0,
            "enthalpy should jump across melting point, got {jump}"
        );
    }
    #[test]
    fn test_ablation_no_recession_below_threshold() {
        let ab = AblationModel::new(3000.0, 5e6, 1600.0, 1.0, 30000.0);
        let r = ab.recession_rate(2000.0);
        assert!(r.abs() < 1e-100, "no ablation below threshold, got {r}");
    }
    #[test]
    fn test_ablation_positive_above_threshold() {
        let ab = AblationModel::new(1000.0, 5e6, 1600.0, 1e3, 5000.0);
        let r = ab.recession_rate(2000.0);
        assert!(r > 0.0, "ablation above threshold, got {r}");
    }
    #[test]
    fn test_ablation_mass_loss_rate_positive() {
        let ab = AblationModel::new(1000.0, 5e6, 1600.0, 1e3, 5000.0);
        let mdot = ab.mass_loss_rate(2000.0);
        assert!(mdot > 0.0, "mass loss rate should be positive, got {mdot}");
    }
    #[test]
    fn test_ablation_net_flux_reduced() {
        let ab = AblationModel::new(1000.0, 5e6, 1600.0, 1e3, 5000.0);
        let q_in = 1e7;
        let q_net = ab.net_heat_flux(q_in, 2000.0);
        assert!(q_net < q_in, "ablative cooling should reduce net heat flux");
    }
    #[test]
    fn test_thermal_stress_zero_at_ref_temp() {
        let tsc = ThermalStressCoupling::new(5, 300.0, 300.0, 12e-6, 200e9, 0.3);
        let sigma = tsc.thermal_stress(0);
        assert!(
            sigma.abs() < 1.0,
            "zero ΔT → zero thermal stress, got {sigma}"
        );
    }
    #[test]
    fn test_thermal_stress_compressive_on_heating() {
        let tsc = ThermalStressCoupling::new(5, 400.0, 300.0, 12e-6, 200e9, 0.3);
        let sigma = tsc.thermal_stress(0);
        assert!(
            sigma < 0.0,
            "constrained heating → compressive stress, got {sigma}"
        );
    }
    #[test]
    fn test_thermal_stress_max_positive() {
        let tsc = ThermalStressCoupling::new(3, 400.0, 300.0, 12e-6, 200e9, 0.3);
        let max_s = tsc.max_thermal_stress();
        assert!(max_s > 0.0, "max stress should be positive");
    }
    #[test]
    fn test_thermal_stress_has_fracture() {
        let tsc = ThermalStressCoupling::new(3, 1000.0, 300.0, 12e-6, 200e9, 0.3);
        assert!(
            tsc.has_fracture(100e6),
            "high ΔT should cause fracture for low fracture stress"
        );
    }
    #[test]
    fn test_thermal_strain_formula() {
        let tsc = ThermalStressCoupling::new(1, 400.0, 300.0, 12e-6, 200e9, 0.3);
        let eps = tsc.thermal_strain(0);
        let expected = 12e-6 * 100.0;
        assert!(
            (eps - expected).abs() < 1e-15,
            "thermal strain = α·ΔT = {expected}, got {eps}"
        );
    }
    #[test]
    fn test_jc_isothermal_reference_state() {
        let jc = JohnsonCookModel::aisi_4340();
        let sig = jc.flow_stress(0.0, jc.eps_dot0, jc.t_room);
        assert!(
            (sig - jc.a).abs() / jc.a < 1e-10,
            "flow stress at reference = A, got {sig}"
        );
    }
    #[test]
    fn test_jc_strain_hardening() {
        let jc = JohnsonCookModel::aisi_4340();
        let sig0 = jc.flow_stress(0.0, jc.eps_dot0, jc.t_room);
        let sig1 = jc.flow_stress(0.5, jc.eps_dot0, jc.t_room);
        assert!(
            sig1 > sig0,
            "flow stress should increase with plastic strain"
        );
    }
    #[test]
    fn test_jc_strain_rate_sensitivity() {
        let jc = JohnsonCookModel::aisi_4340();
        let ratio = jc.strain_rate_sensitivity(1000.0 * jc.eps_dot0, 0.1, jc.t_room);
        assert!(
            ratio > 1.0,
            "higher strain rate → higher flow stress, ratio = {ratio}"
        );
    }
    #[test]
    fn test_jc_thermal_softening() {
        let jc = JohnsonCookModel::aisi_4340();
        let sig_cold = jc.flow_stress(0.1, jc.eps_dot0, jc.t_room);
        let sig_hot = jc.flow_stress(0.1, jc.eps_dot0, 0.5 * (jc.t_room + jc.t_melt));
        assert!(sig_hot < sig_cold, "higher temperature → softer material");
    }
    #[test]
    fn test_jc_at_melt_temperature_zero() {
        let jc = JohnsonCookModel::aisi_4340();
        let sig = jc.flow_stress(0.1, jc.eps_dot0, jc.t_melt);
        assert!(
            sig.abs() < 1e-6,
            "flow stress at T_melt should be ~0, got {sig}"
        );
    }
    #[test]
    fn test_jc_homologous_temperature_bounds() {
        let jc = JohnsonCookModel::aisi_4340();
        assert!((jc.homologous_temperature(jc.t_room)).abs() < 1e-12);
        assert!((jc.homologous_temperature(jc.t_melt) - 1.0).abs() < 1e-12);
        assert!(
            (jc.homologous_temperature(200.0)).abs() < 1e-12,
            "below T_room → T* = 0"
        );
        assert!(
            (jc.homologous_temperature(2000.0) - 1.0).abs() < 1e-12,
            "above T_melt → T* = 1"
        );
    }
    #[test]
    fn test_enthalpy_liquid_fraction_solid() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 900.0);
        assert!((em.liquid_fraction(800.0) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_enthalpy_liquid_fraction_liquid() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 900.0);
        assert!((em.liquid_fraction(950.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_enthalpy_liquid_fraction_midpoint() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 950.0);
        let fl = em.liquid_fraction(900.0);
        assert!(
            (fl - 0.5).abs() < 1e-10,
            "midpoint liquid fraction = 0.5, got {fl}"
        );
    }
    #[test]
    fn test_enthalpy_roundtrip() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 900.0);
        let t_test = 1100.0_f64;
        let h = em.enthalpy_at(t_test);
        let t_recovered = em.temperature_from_enthalpy(h);
        assert!(
            (t_recovered - t_test).abs() < 0.01,
            "enthalpy roundtrip: expected {t_test}, got {t_recovered}"
        );
    }
    #[test]
    fn test_enthalpy_effective_cp_enhanced_in_mushy() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 950.0);
        let cp_solid = em.effective_cp(800.0);
        let cp_mushy = em.effective_cp(900.0);
        assert!(cp_mushy > cp_solid, "effective cp enhanced in mushy zone");
    }
    #[test]
    fn test_enthalpy_increases_with_temperature() {
        let em = EnthalpyMethod::new(2700.0, 900.0, 396e3, 850.0, 900.0);
        let h1 = em.enthalpy_at(700.0);
        let h2 = em.enthalpy_at(1000.0);
        assert!(h2 > h1, "enthalpy increases with temperature");
    }
    #[test]
    fn test_thermal_shock_r_first_positive() {
        let tsp = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 120.0, 3e6);
        let r = tsp.r_first();
        assert!(r > 0.0, "R_first should be positive, got {r}");
    }
    #[test]
    fn test_thermal_shock_r_second_greater_than_r_first() {
        let tsp = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 120.0, 3e6);
        let r1 = tsp.r_first();
        let r2 = tsp.r_second();
        assert!(r2 > r1, "R' = k·R should be greater than R for k > 1");
    }
    #[test]
    fn test_thermal_shock_r_third_positive() {
        let tsp = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 120.0, 3e6);
        let r3 = tsp.r_third();
        assert!(r3 > 0.0, "R'' should be positive, got {r3}");
    }
    #[test]
    fn test_thermal_shock_hasselman_positive() {
        let tsp = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 120.0, 3e6);
        let rh = tsp.r_hasselman();
        assert!(rh > 0.0, "Hasselman R'''' should be positive, got {rh}");
    }
    #[test]
    fn test_thermal_shock_higher_k_gives_higher_r_second() {
        let tsp1 = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 50.0, 3e6);
        let tsp2 = ThermalShockParam::new(400e6, 410e9, 0.14, 4e-6, 200.0, 3e6);
        assert!(tsp2.r_second() > tsp1.r_second(), "higher k → higher R'");
    }
    #[test]
    fn test_thermal_fatigue_life_infinite_at_zero_dt() {
        let tf = ThermalFatigue::new(0.5, 0.5, 12e-6, 200e9);
        let n = tf.fatigue_life(0.0);
        assert!(n.is_infinite(), "zero ΔT → infinite life");
    }
    #[test]
    fn test_thermal_fatigue_life_decreases_with_dt() {
        let tf = ThermalFatigue::new(0.5, 0.5, 12e-6, 200e9);
        let n1 = tf.fatigue_life(200.0);
        let n2 = tf.fatigue_life(500.0);
        assert!(n1 > n2, "larger ΔT → shorter fatigue life: {n1} vs {n2}");
    }
    #[test]
    fn test_thermal_fatigue_accumulation() {
        let mut tf = ThermalFatigue::new(0.5, 0.5, 12e-6, 200e9);
        let n_f = tf.fatigue_life(300.0);
        tf.accumulate(300.0, n_f * 0.5);
        assert!(
            (tf.damage - 0.5).abs() < 1e-9,
            "half life → D=0.5, got {}",
            tf.damage
        );
    }
    #[test]
    fn test_thermal_fatigue_failure_after_full_cycle() {
        let mut tf = ThermalFatigue::new(0.5, 0.5, 12e-6, 200e9);
        let n_f = tf.fatigue_life(300.0);
        tf.accumulate(300.0, n_f * 1.5);
        assert!(
            tf.is_failed(),
            "should be failed after exceeding fatigue life"
        );
    }
    #[test]
    fn test_thermal_fatigue_critical_delta_t() {
        let tf = ThermalFatigue::new(0.5, 0.5, 12e-6, 200e9);
        let n = 1000.0;
        let dt_crit = tf.critical_delta_t(n);
        let n_at_dt = tf.fatigue_life(dt_crit);
        assert!(
            (n_at_dt - n).abs() / n < 0.01,
            "critical ΔT roundtrip: N = {n_at_dt}, expected {n}"
        );
    }
    #[test]
    fn test_haz_diffusivity() {
        let haz = HeatAffectedZone::new(1e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let alpha = haz.diffusivity();
        assert!((alpha - 50.0 / (7850.0 * 500.0)).abs() < 1e-12);
    }
    #[test]
    fn test_haz_halfwidth_positive() {
        let haz = HeatAffectedZone::new(1e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let w = haz.haz_halfwidth(1000.0);
        assert!(w > 0.0, "HAZ half-width should be positive, got {w}");
    }
    #[test]
    fn test_haz_halfwidth_increases_with_heat_input() {
        let haz1 = HeatAffectedZone::new(0.5e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let haz2 = HeatAffectedZone::new(2e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let w1 = haz1.haz_halfwidth(1000.0);
        let w2 = haz2.haz_halfwidth(1000.0);
        assert!(w2 > w1, "higher heat input → wider HAZ: {w1} vs {w2}");
    }
    #[test]
    fn test_haz_cooling_rate_positive() {
        let haz = HeatAffectedZone::new(1e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let cr = haz.cooling_rate_centreline(1000.0);
        assert!(cr > 0.0, "cooling rate should be positive, got {cr}");
    }
    #[test]
    fn test_haz_cooling_rate_increases_with_lower_heat_input() {
        let haz_low = HeatAffectedZone::new(0.1e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let haz_high = HeatAffectedZone::new(5e6, 7850.0, 500.0, 50.0, 0.01, 300.0, 1800.0);
        let cr_low = haz_low.cooling_rate_centreline(1000.0);
        let cr_high = haz_high.cooling_rate_centreline(1000.0);
        assert!(
            cr_low > cr_high,
            "lower heat input → faster cooling: {cr_low} vs {cr_high}"
        );
    }
    #[test]
    fn test_haz_martensite_detection() {
        let haz = HeatAffectedZone::new(100.0, 7850.0, 500.0, 50.0, 0.001, 300.0, 1800.0);
        let cr = haz.cooling_rate_centreline(800.0);
        assert!(
            haz.forms_martensite(1.0, 800.0),
            "fast cooling (cr={cr}) should form martensite at 1 K/s threshold"
        );
    }
    #[test]
    fn test_thermal_shock_resistance_steel() {
        let mat = ThermalMaterial::steel();
        let r = mat.compute_thermal_shock_resistance(250.0e6, 200.0e9, 0.3);
        let expected = 250.0e6 * 0.7 / (200.0e9 * 12.0e-6);
        assert!(
            (r - expected).abs() < 1e-6,
            "R_steel={r}, expected={expected}"
        );
    }
    #[test]
    fn test_thermal_shock_resistance_increases_with_strength() {
        let mat = ThermalMaterial::aluminum();
        let r1 = mat.compute_thermal_shock_resistance(100.0e6, 70.0e9, 0.33);
        let r2 = mat.compute_thermal_shock_resistance(300.0e6, 70.0e9, 0.33);
        assert!(r2 > r1, "higher strength → better TSR: {r1} vs {r2}");
    }
    #[test]
    fn test_thermal_shock_resistance_positive() {
        let mat = ThermalMaterial::silicon_carbide();
        let r = mat.compute_thermal_shock_resistance(500.0e6, 420.0e9, 0.17);
        assert!(r > 0.0, "TSR should be positive for SiC, got {r}");
    }
    #[test]
    fn test_wiedemann_franz_copper() {
        let mat = ThermalMaterial::copper();
        let kappa_el = mat.compute_wiedemann_franz(300.0, 1.7e-8);
        assert!(kappa_el > 300.0, "Cu κ_el should be ~430: {kappa_el}");
        assert!(kappa_el < 600.0, "Cu κ_el should be ~430: {kappa_el}");
    }
    #[test]
    fn test_wiedemann_franz_linear_in_temperature() {
        let mat = ThermalMaterial::steel();
        let rho_e = 1.0e-7_f64;
        let kappa1 = mat.compute_wiedemann_franz(200.0, rho_e);
        let kappa2 = mat.compute_wiedemann_franz(400.0, rho_e);
        assert!(
            (kappa2 / kappa1 - 2.0).abs() < 1e-12,
            "κ_el linear in T: ratio={}",
            kappa2 / kappa1
        );
    }
    #[test]
    fn test_wiedemann_franz_higher_resistivity_lower_conductivity() {
        let mat = ThermalMaterial::steel();
        let k1 = mat.compute_wiedemann_franz(300.0, 1.0e-7);
        let k2 = mat.compute_wiedemann_franz(300.0, 5.0e-7);
        assert!(k1 > k2, "lower ρ_e → higher κ_el: {k1} vs {k2}");
    }
    #[test]
    fn test_volumetric_strain_isotropic_identity() {
        let alpha = 12.0e-6_f64;
        let te = ThermalExpansion::isotropic(alpha, 293.0);
        let delta_t = 100.0_f64;
        let vol_strain = te.compute_volumetric_strain(293.0 + delta_t);
        let expected = 3.0 * alpha * delta_t;
        assert!(
            (vol_strain - expected).abs() < 1e-20,
            "ΔV/V={vol_strain}, expected={expected}"
        );
    }
    #[test]
    fn test_volumetric_strain_zero_at_reference() {
        let te = ThermalExpansion::isotropic(15.0e-6, 300.0);
        let vs = te.compute_volumetric_strain(300.0);
        assert!(vs.abs() < 1e-30, "no strain at reference temperature: {vs}");
    }
    #[test]
    fn test_volumetric_strain_orthotropic() {
        let te = ThermalExpansion::orthotropic(10.0e-6, 20.0e-6, 30.0e-6, 293.0);
        let vs = te.compute_volumetric_strain(293.0 + 100.0);
        let expected = (10.0e-6 + 20.0e-6 + 30.0e-6) * 100.0;
        assert!(
            (vs - expected).abs() < 1e-20,
            "ortho ΔV/V={vs}, expected={expected}"
        );
    }
    #[test]
    fn test_volumetric_strain_negative_for_cooling() {
        let te = ThermalExpansion::isotropic(12.0e-6, 400.0);
        let vs = te.compute_volumetric_strain(300.0);
        assert!(
            vs < 0.0,
            "cooling should give negative volumetric strain: {vs}"
        );
    }
}
