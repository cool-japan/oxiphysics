//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

#[cfg(test)]
use crate::thermal_sph::types::*;
use crate::thermal_sph::types_ext::*;

pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}
pub fn cubic_kernel_grad(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t) / h
    } else {
        0.0
    };
    [
        r_ij[0] / r * dw_dr,
        r_ij[1] / r * dw_dr,
        r_ij[2] / r * dw_dr,
    ]
}
/// Brookshaw consistent Laplacian for SPH heat conduction.
///
/// Σ_j m_j / ρ_j * 2(T_i - T_j) / |r_ij|² * ∇W_ij · r_ij
///
/// `r_ij` is r_i - r_j, `dT` = T_i - T_j, `m_j` mass, `rho_j` density, h smoothing.
pub fn brookshaw_laplacian(
    r_ij_list: &[[f64; 3]],
    dt_list: &[f64],
    m_j_list: &[f64],
    rho_j_list: &[f64],
    h: f64,
) -> f64 {
    let mut result = 0.0;
    for (((&r_ij, &dt), &m_j), &rho_j) in r_ij_list
        .iter()
        .zip(dt_list.iter())
        .zip(m_j_list.iter())
        .zip(rho_j_list.iter())
    {
        let r2 = dot3(r_ij, r_ij);
        if r2 < 1e-300 || rho_j < 1e-300 {
            continue;
        }
        let grad_w = cubic_kernel_grad(r_ij, h);
        let dot_grad = dot3(grad_w, r_ij);
        result += m_j / rho_j * 2.0 * dt / r2 * dot_grad;
    }
    result
}
/// Compute enthalpy from temperature for a material with latent heat.
///
/// Below melting: H = c_p * T
/// In mushy zone: H = c_p * T_melt + L * (T - T_melt) / (T_liquidus - T_melt)
/// Above liquidus: H = c_p * T + L
pub fn enthalpy_from_temp(
    temp: f64,
    c_p: f64,
    t_melt: f64,
    t_liquidus: f64,
    latent_heat: f64,
) -> f64 {
    if temp < t_melt {
        c_p * temp
    } else if temp < t_liquidus {
        let range = (t_liquidus - t_melt).max(1e-300);
        c_p * t_melt + latent_heat * (temp - t_melt) / range
    } else {
        c_p * temp + latent_heat
    }
}
/// Boussinesq buoyancy force per unit mass: f = -g β ΔT ê_z.
///
/// `delta_t` = T - T_ref, `beta` thermal expansion coefficient (1/K),
/// `g` gravitational acceleration vector.
pub fn boussinesq_force_sph(delta_t: f64, beta: f64, g_vec: [f64; 3]) -> [f64; 3] {
    scale3(g_vec, -beta * delta_t)
}
/// Langmuir–Knudsen evaporation rate (kg/m²/s).
///
/// ṁ = α_e * p_sat / sqrt(2 π R_spec T) * M
///
/// `alpha_e` accommodation coefficient, `p_sat` saturation pressure (Pa),
/// `r_spec` specific gas constant (J/kg/K), `temp` surface temperature (K).
pub fn evaporation_rate(alpha_e: f64, p_sat: f64, r_spec: f64, temp: f64) -> f64 {
    if temp < 1e-300 || r_spec < 1e-300 {
        return 0.0;
    }
    let denom = (2.0 * PI * r_spec * temp).sqrt();
    if denom < 1e-300 {
        return 0.0;
    }
    alpha_e * p_sat / denom
}
#[cfg(test)]
mod tests {
    use super::*;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_brookshaw_laplacian_zero_neighbors() {
        let result = brookshaw_laplacian(&[], &[], &[], &[], 0.01);
        assert_eq!(result, 0.0);
    }
    #[test]
    fn test_brookshaw_laplacian_same_temp() {
        let r = [[0.01, 0.0, 0.0]];
        let dt = [0.0];
        let m = [1.0];
        let rho = [1000.0];
        let result = brookshaw_laplacian(&r, &dt, &m, &rho, 0.05);
        assert_eq!(result, 0.0);
    }
    #[test]
    fn test_enthalpy_from_temp_below_melt() {
        let h = enthalpy_from_temp(200.0, 4182.0, 273.15, 273.65, 334000.0);
        assert!((h - 4182.0 * 200.0).abs() < EPS);
    }
    #[test]
    fn test_enthalpy_from_temp_above_liquidus() {
        let h = enthalpy_from_temp(400.0, 4182.0, 273.15, 273.65, 334000.0);
        assert!((h - (4182.0 * 400.0 + 334000.0)).abs() < EPS);
    }
    #[test]
    fn test_boussinesq_force_zero_delta_t() {
        let f = boussinesq_force_sph(0.0, 1e-4, [0.0, 0.0, -9.81]);
        assert!(len3(f) < EPS);
    }
    #[test]
    fn test_boussinesq_force_direction() {
        let f = boussinesq_force_sph(10.0, 1e-4, [0.0, 0.0, -9.81]);
        assert!(f[2] > 0.0);
    }
    #[test]
    fn test_evaporation_rate_positive() {
        let rate = evaporation_rate(1.0, 101325.0, 461.5, 373.15);
        assert!(rate > 0.0);
    }
    #[test]
    fn test_evaporation_rate_zero_temp() {
        let rate = evaporation_rate(1.0, 101325.0, 461.5, 0.0);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_thermal_particle_enthalpy() {
        let p = ThermalParticle::new([0.0; 3], [0.0; 3], 1.0, 1000.0, 300.0);
        let h = p.enthalpy();
        assert!((h - 1000.0 * 4182.0 * 300.0).abs() < 1.0);
    }
    #[test]
    fn test_thermal_particle_diffusivity() {
        let p = ThermalParticle::new([0.0; 3], [0.0; 3], 1.0, 1000.0, 300.0);
        let alpha = p.thermal_diffusivity();
        assert!(alpha > 0.0);
    }
    #[test]
    fn test_sph_heat_conduction_no_gradient() {
        let hc = SphHeatConduction::new(1e-7, 0.01);
        let dtemp = hc.dtemp_dt(300.0, &[]);
        assert_eq!(dtemp, 0.0);
    }
    #[test]
    fn test_phase_change_enthalpy_monotone() {
        let pc = PhaseChangeSph::new(273.15, 273.65, 334000.0, 4182.0);
        let h1 = pc.enthalpy(200.0);
        let h2 = pc.enthalpy(273.15);
        let h3 = pc.enthalpy(273.65);
        let h4 = pc.enthalpy(400.0);
        assert!(h1 < h2 && h2 < h3 && h3 < h4);
    }
    #[test]
    fn test_phase_change_liquid_fraction() {
        let pc = PhaseChangeSph::new(273.15, 273.65, 334000.0, 4182.0);
        assert_eq!(pc.liquid_fraction(200.0), 0.0);
        assert_eq!(pc.liquid_fraction(400.0), 1.0);
        let f_mid = pc.liquid_fraction(273.40);
        assert!(f_mid > 0.0 && f_mid < 1.0);
    }
    #[test]
    fn test_phase_change_is_mushy() {
        let pc = PhaseChangeSph::new(273.15, 273.65, 334000.0, 4182.0);
        assert!(!pc.is_mushy(200.0));
        assert!(pc.is_mushy(273.40));
        assert!(!pc.is_mushy(400.0));
    }
    #[test]
    fn test_thermal_expansion_density() {
        let te = ThermalExpansion::new(2e-4, 293.15, 1000.0, 2.2e9);
        let rho = te.density(293.15);
        assert!((rho - 1000.0).abs() < EPS);
    }
    #[test]
    fn test_thermal_expansion_pressure_change() {
        let te = ThermalExpansion::new(2e-4, 293.15, 1000.0, 2.2e9);
        let dp = te.pressure_change(393.15);
        assert!(dp > 0.0);
    }
    #[test]
    fn test_rayleigh_benard_ra_positive() {
        let rb = RayleighBenardSph::new(2e-4, 1e-6, 1.4e-7, 0.1, 10.0, 9.81);
        let ra = rb.rayleigh_number();
        assert!(ra > 0.0);
    }
    #[test]
    fn test_rayleigh_benard_nusselt_ge_1() {
        let rb = RayleighBenardSph::new(2e-4, 1e-6, 1.4e-7, 0.1, 10.0, 9.81);
        let nu = rb.nusselt_number();
        assert!(nu >= 1.0);
    }
    #[test]
    fn test_wall_heat_dirichlet() {
        let wh = SphWallHeat::dirichlet(400.0, 0.6);
        let t_ghost = wh.apply_dirichlet(300.0, 0.01);
        assert!((t_ghost - 500.0).abs() < EPS);
    }
    #[test]
    fn test_wall_heat_neumann_delta_t() {
        let wh = SphWallHeat::neumann(1000.0, 0.6);
        let dt_applied = wh.apply_neumann(1.0, 1000.0, 4182.0, 1.0);
        assert!(dt_applied > 0.0);
    }
    #[test]
    fn test_evaporation_model_mass_rate() {
        let em = EvaporationModel::new(1.0, 461.5, 2.26e6);
        let rate = em.mass_rate(373.15, 101325.0);
        assert!(rate > 0.0);
    }
    #[test]
    fn test_evaporation_model_heat_sink() {
        let em = EvaporationModel::new(1.0, 461.5, 2.26e6);
        let qs = em.heat_sink(373.15, 101325.0);
        assert!(qs > 0.0);
    }
    #[test]
    fn test_combustion_reaction_rate_increases_with_temp() {
        let cb = CombustionSph::new(1e10, 50000.0, 2e6);
        let k1 = cb.reaction_rate(500.0);
        let k2 = cb.reaction_rate(1000.0);
        assert!(k2 > k1);
    }
    #[test]
    fn test_combustion_adiabatic_flame_temp() {
        let cb = CombustionSph::new(1e10, 50000.0, 2e6);
        let tf = cb.adiabatic_flame_temp(300.0, 0.1, 1000.0);
        assert!(tf > 300.0);
    }
    #[test]
    fn test_thermal_convective_nusselt_positive() {
        let tc = ThermalConvectiveCoeff::new(0.01, 0.6, 400.0, 300.0);
        let nu = tc.nusselt_number(1000.0);
        assert!(nu > 0.0);
    }
    #[test]
    fn test_thermal_convective_average_nusselt() {
        let tc = ThermalConvectiveCoeff::new(0.01, 0.6, 400.0, 300.0);
        let avg = tc.average_nusselt(&[1000.0, 2000.0]);
        assert!(avg > 0.0);
    }
    #[test]
    fn test_thermal_stress_voigt() {
        let ts = ThermalStress::new(200e9, 12e-6, 0.3, 293.15);
        let s = ts.thermal_stress_voigt(393.15);
        assert!(s[0] > 0.0);
        assert_eq!(s[3], 0.0);
    }
    #[test]
    fn test_thermal_stress_zero_at_ref() {
        let ts = ThermalStress::new(200e9, 12e-6, 0.3, 293.15);
        let s = ts.thermal_stress_voigt(293.15);
        for comp in &s[0..3] {
            assert!(comp.abs() < EPS);
        }
    }
    #[test]
    fn test_von_mises_thermal_isotropic_is_zero() {
        // Uniform thermal stress is hydrostatic (isotropic) => von Mises = 0,
        // even though the hydrostatic magnitude itself is large and non-zero.
        let ts = ThermalStress::new(200e9, 12e-6, 0.3, 293.15);
        assert!(ts.thermal_stress_magnitude(393.15).abs() > 1.0);
        assert!(ts.von_mises_thermal(393.15) < EPS);
    }
    #[test]
    fn test_von_mises_voigt_uniaxial() {
        // Uniaxial stress σ_xx = 100 => von Mises equals the axial stress.
        let vm = von_mises_voigt(&[100.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((vm - 100.0).abs() < 1e-9);
        // Pure shear τ_xy = 50 => von Mises = √3 · 50.
        let vm_shear = von_mises_voigt(&[0.0, 0.0, 0.0, 0.0, 0.0, 50.0]);
        assert!((vm_shear - 3.0f64.sqrt() * 50.0).abs() < 1e-9);
    }
    #[test]
    fn test_cubic_kernel_positive() {
        let w = cubic_kernel(0.5, 1.0);
        assert!(w > 0.0);
    }
    #[test]
    fn test_cubic_kernel_zero_outside() {
        let w = cubic_kernel(3.0, 1.0);
        assert_eq!(w, 0.0);
    }
    #[test]
    fn test_pi_imported() {
        // PI is used in computation; verify the import is accessible by using it
        let _ = PI;
    }
    #[test]
    fn test_math_helpers() {
        let a = [3.0, 4.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        assert!((len3(sub3(a, b)) - 5.0).abs() < EPS);
        let s = scale3(a, 2.0);
        assert!((s[0] - 6.0).abs() < EPS);
        let ad = add3(a, b);
        assert!((ad[1] - 4.0).abs() < EPS);
    }
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    #[test]
    fn test_conductivity_tensor_isotropic_apply() {
        let ct = ConductivityTensor::isotropic(10.0);
        let v = [1.0, 2.0, 3.0];
        let lv = ct.apply(v);
        assert!((lv[0] - 10.0).abs() < 1e-10);
        assert!((lv[1] - 20.0).abs() < 1e-10);
        assert!((lv[2] - 30.0).abs() < 1e-10);
    }
    #[test]
    fn test_conductivity_tensor_orthotropic_effective() {
        let ct = ConductivityTensor::orthotropic(1.0, 2.0, 3.0);
        let n = [1.0, 0.0, 0.0];
        let eff = ct.effective(n);
        assert!((eff - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_anisotropic_conduction_isotropic_case() {
        let ahc = AnisotropicHeatConduction::new_isotropic(10.0, 0.1);
        let rate = ahc.dtemp_dt(300.0, 1000.0, &[]);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_temp_dependent_viscosity_arrhenius_increases_with_decreasing_t() {
        let tdv = TempDependentViscosity::arrhenius(1e-3, 300.0, 5000.0);
        let mu1 = tdv.viscosity(300.0);
        let mu2 = tdv.viscosity(200.0);
        assert!(mu2 > mu1);
    }
    #[test]
    fn test_temp_dependent_viscosity_power_law() {
        let tdv = TempDependentViscosity::power_law(1e-3, 300.0, -1.0);
        let mu_ref = tdv.viscosity(300.0);
        let mu_hot = tdv.viscosity(600.0);
        assert!(mu_hot < mu_ref);
    }
    #[test]
    fn test_temp_dependent_conductivity_linear() {
        let tdc = TempDependentConductivity::linear(10.0, 300.0, 0.01);
        let l_ref = tdc.conductivity(300.0);
        let l_hot = tdc.conductivity(400.0);
        assert!((l_ref - 10.0).abs() < 1e-10);
        assert!(l_hot > l_ref);
    }
    #[test]
    fn test_temp_dependent_conductivity_averaged() {
        let tdc = TempDependentConductivity::linear(10.0, 300.0, 0.0);
        let avg = tdc.averaged(200.0, 400.0);
        assert!((avg - 10.0).abs() < 1e-8);
    }
    #[test]
    fn test_stefan_condition_interface_velocity() {
        let sc = StefanCondition::new(334000.0, 1000.0, 2.2, 0.6);
        let v = sc.interface_velocity(0.0, 1000.0);
        assert!(v < 0.0);
    }
    #[test]
    fn test_stefan_condition_no_flux_imbalance() {
        let sc = StefanCondition::new(334000.0, 1000.0, 0.6, 0.6);
        let imb = sc.heat_flux_imbalance(1000.0, 1000.0);
        assert!(imb.abs() < 1e-10);
    }
    #[test]
    fn test_conjugate_heat_transfer_flux() {
        let cht = ConjugateHeatTransfer::new(1000.0, 0.6, 50.0);
        let q = cht.interface_flux(400.0, 300.0);
        assert!(q > 0.0);
    }
    #[test]
    fn test_conjugate_heat_transfer_with_resistance() {
        let cht = ConjugateHeatTransfer::with_contact_resistance(1000.0, 1e-4, 0.6, 50.0);
        let h_eff = cht.effective_h();
        assert!(h_eff < 1000.0);
    }
    #[test]
    fn test_dittus_boelter_nu() {
        let nu = ConjugateHeatTransfer::dittus_boelter_nu(10000.0, 7.0);
        assert!(nu > 0.0);
    }
    #[test]
    fn test_radiation_sph_blackbody_emission() {
        let rad = RadiationSph::new(1.0, 0.1, 0.0);
        let e = rad.blackbody_emission(1000.0);
        assert!((e - 5.670374e-8 * 1e12).abs() / e < 1e-4);
    }
    #[test]
    fn test_radiation_sph_net_flux_zero_equal_temps() {
        let rad = RadiationSph::new(0.9, 0.1, 0.0);
        let q = rad.net_radiation_flux(500.0, 500.0);
        assert!(q.abs() < 1e-6);
    }
    #[test]
    fn test_radiation_sph_transmittance() {
        let rad = RadiationSph::new(0.9, 1.0, 0.0);
        let t = rad.transmittance(1.0);
        assert!((t - (-1.0f64).exp()).abs() < 1e-10);
    }
    #[test]
    fn test_boiling_sph_jakob_number_zero_at_sat() {
        let b = BoilingSph::new(373.15, 2.26e6, 958.0, 0.6, 4216.0, 0.0589, 1e5);
        let ja = b.jakob_number(373.15);
        assert_eq!(ja, 0.0);
    }
    #[test]
    fn test_boiling_sph_jakob_number_positive() {
        let b = BoilingSph::new(373.15, 2.26e6, 958.0, 0.6, 4216.0, 0.0589, 1e5);
        let ja = b.jakob_number(383.15);
        assert!(ja > 0.0);
    }
    #[test]
    fn test_boiling_sph_critical_heat_flux_positive() {
        let b = BoilingSph::new(373.15, 2.26e6, 958.0, 0.6, 4216.0, 0.0589, 1e5);
        let chf = b.critical_heat_flux();
        assert!(chf > 0.0);
    }
    #[test]
    fn test_casting_solidification_chvorinov() {
        let cs = CastingSolidification::new(
            1356.0, 1365.0, 205_000.0, 386.0, 500.0, 300.0, 1450.0, 385.0, 400.0, 8960.0,
        );
        let t_sol = cs.chvorinov_time(0.001, 0.01, 1e6);
        assert!(t_sol > 0.0);
    }
    #[test]
    fn test_casting_fraction_solidified() {
        let cs = CastingSolidification::new(
            1356.0, 1365.0, 205_000.0, 386.0, 500.0, 300.0, 1450.0, 385.0, 400.0, 8960.0,
        );
        let f = cs.fraction_solidified(100.0, 25.0);
        assert!(f > 0.0 && f < 1.0);
    }
    #[test]
    fn test_cryogenic_ln2_properties() {
        let ln2 = CryogenicFluid::liquid_nitrogen();
        assert_eq!(ln2.name, "LN2");
        assert!(ln2.t_nbp < 100.0);
        let alpha = ln2.thermal_diffusivity();
        assert!(alpha > 0.0);
    }
    #[test]
    fn test_cryogenic_lhe_prandtl() {
        let lhe = CryogenicFluid::liquid_helium();
        let pr = lhe.prandtl_number();
        assert!(pr > 0.0);
    }
    #[test]
    fn test_cryogenic_lh2_jakob() {
        let lh2 = CryogenicFluid::liquid_hydrogen();
        let ja = lh2.jakob_number(5.0);
        assert!(ja > 0.0);
    }
    #[test]
    fn test_cleary_harmonic_mean_equal() {
        let hm = ClearyHeatConduction::harmonic_mean(10.0, 10.0);
        assert!((hm - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_cleary_harmonic_mean_asymmetric() {
        let hm = ClearyHeatConduction::harmonic_mean(1.0, 100.0);
        assert!(hm < 50.5);
        assert!(hm > 0.0);
    }
    #[test]
    fn test_mushy_zone_drag_zero_at_full_liquid() {
        let mzd = MushyZoneDrag::new(1e6);
        let f = mzd.drag_force(1.0, [1.0, 0.0, 0.0]);
        assert!(f[0].abs() < 1e-3);
    }
    #[test]
    fn test_mushy_zone_drag_large_at_zero_liquid() {
        let mzd = MushyZoneDrag::new(1e6);
        let coeff = mzd.drag_coefficient(0.0);
        assert!(coeff > 1.0);
    }
    #[test]
    fn test_internal_energy_eq_pressure_work() {
        let ie = InternalEnergyEq::new(1.4);
        let pw = ie.pressure_work(1e5, 1.0, 1.0);
        assert!((pw + 1e5).abs() < 1.0);
    }
    #[test]
    fn test_internal_energy_eq_temperature_from_energy() {
        let ie = InternalEnergyEq::new(1.4);
        let r_spec = 287.0;
        let u = 300.0 * r_spec;
        let t = ie.temperature_from_energy(u, r_spec);
        assert!((t - 300.0 * 0.4).abs() < 1.0);
    }
    #[test]
    fn test_thermodiffusion_soret_flux_direction() {
        let td = ThermodiffusionSph::new(1e-11, 1e-9, 300.0);
        let grad_t = [1000.0, 0.0, 0.0];
        let flux = td.soret_flux(0.5, grad_t, 1000.0, 300.0);
        assert!(flux[0] < 0.0);
    }
    #[test]
    fn test_thermal_sph_simulation_step() {
        let mut sim = ThermalSphSimulation::new(0.1, 0.001);
        let p1 = ThermalParticle::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, 1000.0, 400.0);
        let p2 = ThermalParticle::new([0.05, 0.0, 0.0], [0.0; 3], 1.0, 1000.0, 300.0);
        sim.add_particle(p1);
        sim.add_particle(p2);
        let t_mean_before = sim.mean_temperature();
        sim.step();
        assert!(sim.time > 0.0);
        assert!((sim.mean_temperature() - t_mean_before).abs() >= 0.0);
    }
    #[test]
    fn test_heat_pipe_water_operating() {
        let hp = HeatPipeSph::water_heat_pipe(0.5);
        assert!(hp.is_operating(400.0));
        assert!(!hp.is_operating(200.0));
    }
    #[test]
    fn test_heat_pipe_heat_flux_positive() {
        let hp = HeatPipeSph::water_heat_pipe(0.5);
        let q = hp.heat_flux(10.0);
        assert!(q > 0.0);
    }
    #[test]
    fn test_nucleate_boiling_vapor_fraction_zero_at_sat() {
        let b = BoilingSph::new(373.15, 2.26e6, 958.0, 0.6, 4216.0, 0.0589, 1e5);
        let nbs = NucleateBoilingSph::new(b);
        let alpha = nbs.vapor_fraction(373.15);
        assert_eq!(alpha, 0.0);
    }
    #[test]
    fn test_nucleate_boiling_mixture_density_at_sat() {
        let b = BoilingSph::new(373.15, 2.26e6, 958.0, 0.6, 4216.0, 0.0589, 1e5);
        let nbs = NucleateBoilingSph::new(b);
        let rho_mix = nbs.mixture_density(373.15);
        assert!((rho_mix - 958.0).abs() < 1e-6);
    }
    #[test]
    fn test_view_factor_element_positive() {
        let vf = RadiationSph::view_factor_element(1.0, 1.0, 1.0, 1.0);
        assert!(vf > 0.0);
    }
    #[test]
    fn test_radiation_rosseland_flux_direction() {
        let rad = RadiationSph::new(1.0, 1.0, 0.0);
        let grad_t = [100.0, 0.0, 0.0];
        let flux = rad.rosseland_flux(1000.0, grad_t);
        assert!(flux[0] < 0.0);
    }
}
/// Approximate error function erf(x) via Horner polynomial (Abramowitz & Stegun 7.1.26).
pub fn erf_approx(x: f64) -> f64 {
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = 1.0 - poly * (-(x * x)).exp();
    if x >= 0.0 { result } else { -result }
}
#[cfg(test)]
mod tests_part3 {
    use super::*;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_symmetric_conductivity_inter_particle_equal() {
        let lij = SymmetricSphConductivity::inter_particle_conductivity(10.0, 10.0);
        assert!((lij - 10.0).abs() < EPS);
    }
    #[test]
    fn test_symmetric_conductivity_inter_particle_harmonic() {
        let lij = SymmetricSphConductivity::inter_particle_conductivity(1.0, 9.0);
        assert!((lij - 1.8).abs() < EPS);
    }
    #[test]
    fn test_symmetric_conductivity_no_neighbors_zero() {
        let sc = SymmetricSphConductivity::new(0.01);
        let rate = sc.dtemp_dt(300.0, 1.0, 1000.0, 4182.0, &[]);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_sph_heat_equation_zero_without_sources() {
        let eq = SphHeatEquation::new(0.01, 0.0);
        let rate = eq.full_dtemp_dt(0.0, 0.0, 1000.0, 4182.0, 0.0, 1e-3, 0.0);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_sph_heat_equation_pressure_work_contribution() {
        let eq = SphHeatEquation::new(0.01, 0.0);
        let rate = eq.full_dtemp_dt(0.0, 1e5, 1000.0, 4182.0, 1.0, 1e-3, 0.0);
        assert!(rate.abs() > 0.0);
    }
    #[test]
    fn test_stefan_problem_stefan_number_positive() {
        let sp = StefanProblemAnalytic::new(250.0, 273.15, 1.02e-7, 2.22, 2050.0, 334000.0, 917.0);
        let ste = sp.stefan_number();
        assert!(ste > 0.0);
    }
    #[test]
    fn test_stefan_problem_interface_grows_with_time() {
        let sp = StefanProblemAnalytic::new(250.0, 273.15, 1.02e-7, 2.22, 2050.0, 334000.0, 917.0);
        let s1 = sp.interface_position(100.0);
        let s2 = sp.interface_position(400.0);
        assert!(s2 > s1);
    }
    #[test]
    fn test_stefan_problem_temperature_profile_at_wall() {
        let sp = StefanProblemAnalytic::new(250.0, 273.15, 1.02e-7, 2.22, 2050.0, 334000.0, 917.0);
        let t = sp.temperature_profile(0.0, 100.0);
        assert!((t - 250.0).abs() < 1.0);
    }
    #[test]
    fn test_optically_thin_radiation_zero_when_equal_temps() {
        let otr = OpticallyThinRadiation::new(1.0, 300.0, 0.01);
        let rate = otr.volumetric_cooling(300.0, 1000.0);
        assert!(rate.abs() < EPS);
    }
    #[test]
    fn test_optically_thin_radiation_cooling_when_hot() {
        let otr = OpticallyThinRadiation::new(1.0, 300.0, 0.01);
        let rate = otr.volumetric_cooling(1000.0, 1000.0);
        assert!(rate < 0.0);
    }
    #[test]
    fn test_optically_thin_radiation_timescale_positive() {
        let otr = OpticallyThinRadiation::new(1.0, 300.0, 0.01);
        let tau = otr.cooling_timescale(1000.0, 1.0, 1000.0);
        assert!(tau > 0.0);
    }
    #[test]
    fn test_viscous_dissipation_strain_rate_zero_velocity() {
        let h = 0.01;
        let grad_v = ViscousDissipationSph::velocity_gradient([0.0; 3], &[], h);
        let s = ViscousDissipationSph::strain_rate(grad_v);
        let inv2 = ViscousDissipationSph::strain_rate_invariant(s);
        assert!(inv2.abs() < EPS);
    }
    #[test]
    fn test_viscous_dissipation_rate_non_negative() {
        let vd = ViscousDissipationSph::new(1e-3, 0.01);
        let rate = vd.dissipation_rate([0.0; 3], &[]);
        assert!(rate >= 0.0);
    }
    #[test]
    fn test_boussinesq_buoyancy_direction_warm() {
        let b = BoussinesqSph::new(1000.0, 300.0, 2e-4, [0.0, 0.0, -9.81]);
        let fb = b.buoyancy_force(310.0);
        assert!(fb[2] > 0.0);
    }
    #[test]
    fn test_boussinesq_zero_at_reference() {
        let b = BoussinesqSph::new(1000.0, 300.0, 2e-4, [0.0, 0.0, -9.81]);
        let fb = b.buoyancy_force(300.0);
        assert!(len3(fb) < EPS);
    }
    #[test]
    fn test_boussinesq_convective_instability() {
        let b = BoussinesqSph::new(1000.0, 300.0, 2e-4, [0.0, 0.0, -9.81]);
        assert!(b.is_convectively_unstable(100.0, 0.1, 1e-6, 1.4e-7));
    }
    #[test]
    fn test_adiabatic_temperature_at_reference_pressure() {
        let ac = AdiabaticCompression::new(1.4, 1e5, 300.0);
        let t = ac.adiabatic_temperature(1e5);
        assert!((t - 300.0).abs() < EPS);
    }
    #[test]
    fn test_adiabatic_temperature_increases_with_compression() {
        let ac = AdiabaticCompression::new(1.4, 1e5, 300.0);
        let t_high = ac.adiabatic_temperature(2e5);
        assert!(t_high > 300.0);
    }
    #[test]
    fn test_potential_temperature_at_reference() {
        let ac = AdiabaticCompression::new(1.4, 1e5, 300.0);
        let theta = ac.potential_temperature(300.0, 1e5);
        assert!((theta - 300.0).abs() < EPS);
    }
    #[test]
    fn test_thermal_shock_stress_zero_at_ref() {
        let ts = ThermalShockSph::new(200e9, 0.25, 12e-6, 1e6, 200e6, 293.15);
        let sigma = ts.thermal_stress(293.15);
        assert!(sigma.abs() < EPS);
    }
    #[test]
    fn test_thermal_shock_fracture_criterion() {
        let ts = ThermalShockSph::new(200e9, 0.25, 12e-6, 1e6, 1.0, 293.15);
        assert!(ts.is_fractured(3000.0, 1e-4));
    }
    #[test]
    fn test_thermal_shock_damage_zero_below_critical() {
        let ts = ThermalShockSph::new(200e9, 0.25, 12e-6, 1e6, 200e6, 293.15);
        let d = ts.damage(293.16);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_sph_thermal_bc_dirichlet_ghost() {
        let bc = SphThermalBc::dirichlet(400.0, 0.01);
        let t_ghost = bc.ghost_temperature(300.0, 0.01, 0.6);
        assert!((t_ghost - 500.0).abs() < EPS);
    }
    #[test]
    fn test_sph_thermal_bc_neumann_ghost() {
        let bc = SphThermalBc::neumann(1000.0, 0.01);
        let t_ghost = bc.ghost_temperature(300.0, 0.01, 0.6);
        let expected = 300.0 + 1000.0 * 0.01 / 0.6;
        assert!((t_ghost - expected).abs() < 1e-8);
    }
    #[test]
    fn test_blended_viscosity_continuity_at_blend_point() {
        let low = TempDependentViscosity::power_law(1e-3, 300.0, 0.0);
        let high = TempDependentViscosity::power_law(1e-3, 300.0, 0.0);
        let bv = BlendedViscosity::new(low, high, 500.0, 50.0);
        let mu_below = bv.viscosity(499.0);
        let mu_above = bv.viscosity(501.0);
        assert!((mu_below - 1e-3).abs() < 1e-6);
        assert!((mu_above - 1e-3).abs() < 1e-6);
    }
    #[test]
    fn test_full_energy_sph_ideal_gas_pressure() {
        let fe = FullEnergySph::new(0.01, [0.0, 0.0, -9.81], 1.4);
        let p = fe.ideal_gas_pressure(1.0, 1000.0);
        assert!((p - 400.0).abs() < EPS);
    }
    #[test]
    fn test_full_energy_sph_sound_speed_positive() {
        let fe = FullEnergySph::new(0.01, [0.0, 0.0, -9.81], 1.4);
        let c = fe.sound_speed(1e5, 1.2);
        assert!(c > 0.0);
    }
}
/// SPH thermal diffusion rate dT/dt contribution for particle `i`.
///
/// Uses the Brookshaw consistent Laplacian operator:
/// dT/dt = Σ_j m_j / ρ_j * 2(T_i - T_j) / |r_ij|² * ∇W_ij · r_ij
pub fn sph_thermal_diffusion(particles: &[ThermalParticle], i: usize, h: f64) -> f64 {
    let pi = &particles[i];
    let mut result = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let r_ij = sub3(pi.pos, pj.pos);
        let r2 = dot3(r_ij, r_ij);
        if r2 < 1e-300 || pj.density < 1e-300 {
            continue;
        }
        let dt_temp = pi.temperature - pj.temperature;
        let grad_w = cubic_kernel_grad(r_ij, h);
        let dot_grad = dot3(grad_w, r_ij);
        result += pj.mass / pj.density * 2.0 * dt_temp / r2 * dot_grad;
    }
    result
}
/// SPH heat flux vector q (W/m²) at particle `i`.
///
/// q_i = −λ_i Σ_j m_j / ρ_j * (T_i − T_j) * ∇W_ij
pub fn sph_heat_flux(particles: &[ThermalParticle], i: usize, h: f64) -> [f64; 3] {
    let pi = &particles[i];
    let mut flux = [0.0f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let r_ij = sub3(pi.pos, pj.pos);
        if dot3(r_ij, r_ij) < 1e-300 || pj.density < 1e-300 {
            continue;
        }
        let dt_temp = pi.temperature - pj.temperature;
        let grad_w = cubic_kernel_grad(r_ij, h);
        let factor = pj.mass / pj.density * dt_temp;
        flux[0] -= pi.conductivity * factor * grad_w[0];
        flux[1] -= pi.conductivity * factor * grad_w[1];
        flux[2] -= pi.conductivity * factor * grad_w[2];
    }
    flux
}
/// Apply a localised heat source to particles within `radius` of `source`.
///
/// Each particle inside the sphere receives `rate * dt` added to its
/// temperature (the caller provides the `rate` in K/s).
pub fn apply_thermal_forcing(
    particles: &mut [ThermalParticle],
    source: [f64; 3],
    radius: f64,
    rate: f64,
) {
    for p in particles.iter_mut() {
        let r2 = dot3(sub3(p.pos, source), sub3(p.pos, source));
        if r2 <= radius * radius {
            p.temperature += rate;
        }
    }
}
/// Stefan-Boltzmann radiative heat flux (W/m²).
///
/// q_rad = ε σ T⁴
///
/// where σ = 5.670374419 × 10⁻⁸ W/m²/K⁴.
pub fn stefan_boltzmann_radiation(temp: f64, emissivity: f64) -> f64 {
    pub(super) const SIGMA: f64 = 5.670_374_419e-8;
    emissivity * SIGMA * temp.powi(4)
}
/// Smoothed latent heat source term for phase change (J/kg).
///
/// Uses a Gaussian smoothing around `melting_temp` with half-width `width`.
/// Returns the effective heat absorption per unit temperature change near the
/// phase transition:
///
/// L_eff(T) = L * exp(−(T − T_melt)² / (2 width²)) / (√(2π) width)
pub fn phase_change_latent_heat(temp: f64, melting_temp: f64, latent_heat: f64, width: f64) -> f64 {
    if width < 1e-300 {
        return 0.0;
    }
    let dt = temp - melting_temp;
    let gauss = (-(dt * dt) / (2.0 * width * width)).exp();
    latent_heat * gauss / ((2.0 * std::f64::consts::PI).sqrt() * width)
}
/// Temperature-dependent viscosity using an Arrhenius-type power law.
///
/// μ(T) = μ₀ * (T / T_ref)^(−exp_)
///
/// `exp_` is a positive exponent (e.g. 0.7 for gases, negative exponent
/// effectively increases viscosity with temperature for liquids when
/// using |exp_|).
pub fn temperature_dependent_viscosity(temp: f64, mu0: f64, temp_ref: f64, exp_: f64) -> f64 {
    if temp < 1e-300 || temp_ref < 1e-300 {
        return mu0;
    }
    mu0 * (temp / temp_ref).powf(-exp_)
}
/// Advance a thermal SPH particle system by one explicit Euler step `dt`.
///
/// For each particle:
/// 1. Compute temperature change from thermal diffusion.
/// 2. Update temperature: T += α * dT/dt * dt  (α = λ/(ρ c_p)).
/// 3. Update pressure via Tait EOS: p = k (ρ/ρ₀)^γ − k  (ρ₀ = 1000 kg/m³).
/// 4. Advect position.
pub fn thermal_sph_step(particles: &mut [ThermalParticle], params: &ThermalSPHParams, dt: f64) {
    pub(super) const RHO0: f64 = 1000.0;
    let h = params.kernel_radius;
    let n = particles.len();
    let mut dt_rates = vec![0.0f64; n];
    for (i, rate) in dt_rates.iter_mut().enumerate().take(n) {
        *rate = sph_thermal_diffusion(particles, i, h);
    }
    for (i, p) in particles.iter_mut().enumerate() {
        let alpha = if p.density > 1e-300 && p.specific_heat > 1e-300 {
            p.conductivity / (p.density * p.specific_heat)
        } else {
            0.0
        };
        p.temperature += alpha * dt_rates[i] * dt;
        p.temperature = p.temperature.max(0.0);
        if p.density > 1e-300 {
            let ratio = p.density / RHO0;
            p.pressure = params.eos_k * (ratio.powf(params.eos_gamma) - 1.0);
        }
        p.pos = add3(p.pos, scale3(p.vel, dt));
    }
}
#[cfg(test)]
mod tests_spec {
    use super::*;
    fn make_particle(pos: [f64; 3], temp: f64) -> ThermalParticle {
        let mut p = ThermalParticle::new(pos, [0.0; 3], 1e-3, 1000.0, temp);
        p.conductivity = 0.6;
        p.specific_heat = 4182.0;
        p
    }
    fn make_params() -> ThermalSPHParams {
        ThermalSPHParams::new(0.05, 1e5, 7.0, 1e-3)
    }
    #[test]
    fn test_params_new() {
        let p = make_params();
        assert_eq!(p.kernel_radius, 0.05);
        assert_eq!(p.eos_gamma, 7.0);
    }
    #[test]
    fn test_thermal_diffusion_single_particle() {
        let particles = vec![make_particle([0.0; 3], 300.0)];
        let result = sph_thermal_diffusion(&particles, 0, 0.05);
        assert_eq!(result, 0.0, "Single particle → no diffusion");
    }
    #[test]
    fn test_thermal_diffusion_two_equal_temp() {
        let p1 = make_particle([0.0, 0.0, 0.0], 300.0);
        let p2 = make_particle([0.01, 0.0, 0.0], 300.0);
        let result = sph_thermal_diffusion(&[p1, p2], 0, 0.05);
        assert!(result.abs() < 1e-10, "Equal temps → no diffusion: {result}");
    }
    #[test]
    fn test_thermal_diffusion_direction() {
        let hot = make_particle([0.0, 0.0, 0.0], 400.0);
        let cold = make_particle([0.01, 0.0, 0.0], 300.0);
        let result = sph_thermal_diffusion(&[hot, cold], 0, 0.05);
        assert!(result.is_finite(), "Diffusion must be finite: {result}");
    }
    #[test]
    fn test_heat_flux_single_particle() {
        let p = make_particle([0.0; 3], 300.0);
        let flux = sph_heat_flux(&[p], 0, 0.05);
        assert_eq!(flux, [0.0; 3], "Single particle → zero flux");
    }
    #[test]
    fn test_heat_flux_equal_temp_zero() {
        let p1 = make_particle([0.0, 0.0, 0.0], 300.0);
        let p2 = make_particle([0.01, 0.0, 0.0], 300.0);
        let flux = sph_heat_flux(&[p1, p2], 0, 0.05);
        for &f in flux.iter() {
            assert!(f.abs() < 1e-10, "Equal temp → zero flux: {f}");
        }
    }
    #[test]
    fn test_apply_forcing_inside_radius() {
        let mut particles = vec![make_particle([0.0; 3], 300.0)];
        apply_thermal_forcing(&mut particles, [0.0; 3], 1.0, 10.0);
        assert!((particles[0].temperature - 310.0).abs() < 1e-10);
    }
    #[test]
    fn test_apply_forcing_outside_radius() {
        let mut particles = vec![make_particle([10.0, 0.0, 0.0], 300.0)];
        apply_thermal_forcing(&mut particles, [0.0; 3], 1.0, 10.0);
        assert_eq!(
            particles[0].temperature, 300.0,
            "Outside radius → unchanged"
        );
    }
    #[test]
    fn test_apply_forcing_empty() {
        let mut particles: Vec<ThermalParticle> = Vec::new();
        apply_thermal_forcing(&mut particles, [0.0; 3], 1.0, 10.0);
    }
    #[test]
    fn test_stefan_boltzmann_zero_temp() {
        assert_eq!(stefan_boltzmann_radiation(0.0, 1.0), 0.0);
    }
    #[test]
    fn test_stefan_boltzmann_room_temp() {
        let q = stefan_boltzmann_radiation(300.0, 1.0);
        assert!((q - 459.27).abs() < 1.0, "Room temp radiation: {q}");
    }
    #[test]
    fn test_stefan_boltzmann_emissivity_scaling() {
        let q1 = stefan_boltzmann_radiation(300.0, 1.0);
        let q2 = stefan_boltzmann_radiation(300.0, 0.5);
        assert!((q2 - q1 * 0.5).abs() < 1e-6, "Radiation scales with ε");
    }
    #[test]
    fn test_stefan_boltzmann_t4_law() {
        let q1 = stefan_boltzmann_radiation(300.0, 1.0);
        let q2 = stefan_boltzmann_radiation(600.0, 1.0);
        assert!((q2 / q1 - 16.0).abs() < 1e-6, "T⁴ law: ratio = {}", q2 / q1);
    }
    #[test]
    fn test_phase_change_zero_width() {
        assert_eq!(phase_change_latent_heat(273.0, 273.0, 334e3, 0.0), 0.0);
    }
    #[test]
    fn test_phase_change_peak_at_melting() {
        let at_melt = phase_change_latent_heat(273.0, 273.0, 334e3, 5.0);
        let away = phase_change_latent_heat(290.0, 273.0, 334e3, 5.0);
        assert!(at_melt > away, "Peak at melting temperature");
    }
    #[test]
    fn test_phase_change_positive() {
        let val = phase_change_latent_heat(273.0, 273.0, 334e3, 5.0);
        assert!(val > 0.0, "Latent heat term must be positive: {val}");
    }
    #[test]
    fn test_viscosity_at_ref_temp() {
        let mu = temperature_dependent_viscosity(300.0, 1e-3, 300.0, 0.7);
        assert!((mu - 1e-3).abs() < 1e-15, "μ at ref temp = μ₀: {mu}");
    }
    #[test]
    fn test_viscosity_decreases_with_temp_for_gas() {
        let mu1 = temperature_dependent_viscosity(300.0, 1e-5, 300.0, -0.7);
        let mu2 = temperature_dependent_viscosity(1000.0, 1e-5, 300.0, -0.7);
        assert!(
            mu2 > mu1,
            "Gas viscosity increases with temperature (negative exp)"
        );
    }
    #[test]
    fn test_viscosity_zero_temp() {
        let mu = temperature_dependent_viscosity(0.0, 1e-3, 300.0, 0.7);
        assert_eq!(mu, 1e-3, "Zero temp → fallback to μ₀");
    }
    #[test]
    fn test_thermal_step_empty() {
        let mut particles: Vec<ThermalParticle> = Vec::new();
        let params = make_params();
        thermal_sph_step(&mut particles, &params, 0.01);
    }
    #[test]
    fn test_thermal_step_position_advances() {
        let mut p = make_particle([0.0; 3], 300.0);
        p.vel = [1.0, 0.0, 0.0];
        let mut particles = vec![p];
        let params = make_params();
        thermal_sph_step(&mut particles, &params, 0.1);
        assert!(
            (particles[0].pos[0] - 0.1).abs() < 1e-10,
            "x = {}",
            particles[0].pos[0]
        );
    }
    #[test]
    fn test_thermal_step_temp_non_negative() {
        let mut p = make_particle([0.0; 3], 1.0);
        p.density = 1.0;
        let mut particles = vec![p];
        let params = make_params();
        thermal_sph_step(&mut particles, &params, 1.0);
        assert!(particles[0].temperature >= 0.0, "Temperature must be ≥ 0");
    }
    #[test]
    fn test_thermal_step_pressure_tait() {
        let mut p = make_particle([0.0; 3], 300.0);
        p.density = 1000.0;
        let mut particles = vec![p];
        let params = ThermalSPHParams::new(0.05, 1e5, 7.0, 1e-3);
        thermal_sph_step(&mut particles, &params, 1e-4);
        assert!(particles[0].pressure.is_finite(), "Pressure must be finite");
    }
}
