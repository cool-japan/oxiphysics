//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Gravitational acceleration (m/s²).
pub(crate) const G: f64 = 9.81;
/// Kinematic viscosity of water at 20 °C (m²/s).
pub(crate) const NU_WATER: f64 = 1.0e-6;
/// Density of water (kg/m³).
pub(crate) const RHO_WATER: f64 = 1000.0;
/// Density of quartz sediment (kg/m³).
pub(crate) const RHO_SEDIMENT: f64 = 2650.0;
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Stokes settling velocity for a sphere of diameter `d` (m) in a fluid.
///
/// v_s = (ρ_s − ρ_f) g d² / (18 μ)
///
/// Parameters: `d` grain diameter (m), `rho_s` sediment density (kg/m³),
/// `rho_f` fluid density (kg/m³), `mu` dynamic viscosity (Pa·s).
pub fn settling_velocity_stokes(d: f64, rho_s: f64, rho_f: f64, mu: f64) -> f64 {
    if mu < 1e-300 || d < 0.0 {
        return 0.0;
    }
    (rho_s - rho_f) * G * d * d / (18.0 * mu)
}
/// Shields stress parameter θ = τ_b / ((ρ_s − ρ_f) g d).
///
/// `tau_b` bed shear stress (Pa), `d` grain diameter (m).
pub fn shields_stress(tau_b: f64, d: f64, rho_s: f64, rho_f: f64) -> f64 {
    let denom = (rho_s - rho_f) * G * d;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    tau_b / denom
}
/// Rouse sediment concentration profile C(z).
///
/// C(z) = C_ref * ((a / z) * ((H - z) / (H - a)))^Ro
///
/// Parameters: `z` elevation above bed (m), `h` total water depth (m),
/// `a` reference height (m), `c_ref` reference concentration at `a`,
/// `rouse_number` Ro = w_s / (κ u*).
pub fn rouse_concentration(z: f64, h: f64, a: f64, c_ref: f64, rouse_number: f64) -> f64 {
    if z <= 0.0 || h <= 0.0 || a <= 0.0 || z >= h {
        return 0.0;
    }
    let ratio = (a / z) * ((h - z) / (h - a));
    if ratio < 0.0 {
        return 0.0;
    }
    c_ref * ratio.powf(rouse_number)
}
/// Meyer–Peter–Müller bed-load transport rate q_b (m²/s).
///
/// q_b = 8 * ((θ − θ_cr)^{3/2}) * sqrt((s−1)*g*d³)
///
/// Returns 0 when θ ≤ θ_cr.
pub fn mpm_bedload(theta: f64, theta_cr: f64, d: f64, rho_s: f64, rho_f: f64) -> f64 {
    if theta <= theta_cr {
        return 0.0;
    }
    let s = rho_s / rho_f;
    let term = ((theta - theta_cr).max(0.0)).powf(1.5);
    8.0 * term * ((s - 1.0) * G * d * d * d).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sediment::types::*;
    use crate::sediment::types_advanced::*;
    use std::f64::consts::PI;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_settling_velocity_stokes_positive() {
        let mu = RHO_WATER * NU_WATER;
        let vs = settling_velocity_stokes(1e-3, RHO_SEDIMENT, RHO_WATER, mu);
        assert!(vs > 0.0);
    }
    #[test]
    fn test_settling_velocity_stokes_zero_diameter() {
        let mu = RHO_WATER * NU_WATER;
        let vs = settling_velocity_stokes(0.0, RHO_SEDIMENT, RHO_WATER, mu);
        assert_eq!(vs, 0.0);
    }
    #[test]
    fn test_settling_velocity_stokes_proportional_d2() {
        let mu = RHO_WATER * NU_WATER;
        let v1 = settling_velocity_stokes(1e-3, RHO_SEDIMENT, RHO_WATER, mu);
        let v2 = settling_velocity_stokes(2e-3, RHO_SEDIMENT, RHO_WATER, mu);
        assert!((v2 / v1 - 4.0).abs() < 1e-8);
    }
    #[test]
    fn test_shields_stress_zero_tau() {
        let theta = shields_stress(0.0, 1e-3, RHO_SEDIMENT, RHO_WATER);
        assert_eq!(theta, 0.0);
    }
    #[test]
    fn test_shields_stress_positive() {
        let theta = shields_stress(1.0, 1e-3, RHO_SEDIMENT, RHO_WATER);
        assert!(theta > 0.0);
    }
    #[test]
    fn test_rouse_concentration_zero_at_bed() {
        let c = rouse_concentration(0.0, 1.0, 0.1, 0.01, 1.2);
        assert_eq!(c, 0.0);
    }
    #[test]
    fn test_rouse_concentration_decreases_upward() {
        let c1 = rouse_concentration(0.1, 1.0, 0.1, 0.01, 1.2);
        let c2 = rouse_concentration(0.5, 1.0, 0.1, 0.01, 1.2);
        assert!(c1 >= c2);
    }
    #[test]
    fn test_rouse_concentration_at_ref_height() {
        let c = rouse_concentration(0.1, 1.0, 0.1, 0.05, 1.2);
        assert!((c - 0.05).abs() < 1e-10);
    }
    #[test]
    fn test_mpm_bedload_zero_below_critical() {
        let q = mpm_bedload(0.03, 0.047, 1e-3, RHO_SEDIMENT, RHO_WATER);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_mpm_bedload_positive_above_critical() {
        let q = mpm_bedload(0.1, 0.047, 1e-3, RHO_SEDIMENT, RHO_WATER);
        assert!(q > 0.0);
    }
    #[test]
    fn test_mpm_bedload_monotone() {
        let q1 = mpm_bedload(0.1, 0.047, 1e-3, RHO_SEDIMENT, RHO_WATER);
        let q2 = mpm_bedload(0.2, 0.047, 1e-3, RHO_SEDIMENT, RHO_WATER);
        assert!(q2 > q1);
    }
    #[test]
    fn test_sediment_particle_settling_velocity_positive() {
        let p = SedimentParticle::new([0.0; 3], [0.0; 3], RHO_WATER, 1e-3);
        assert!(p.settling_velocity > 0.0);
    }
    #[test]
    fn test_sediment_particle_kinetic_energy() {
        let p = SedimentParticle::new([0.0; 3], [1.0, 0.0, 0.0], 1000.0, 1e-3);
        assert!((p.kinetic_energy() - 500.0).abs() < 1.0);
    }
    #[test]
    fn test_rouse_profile_at_ref_height() {
        let rp = RouseProfile::new(2.0, 0.2, 0.01, 0.001, 0.05);
        let c = rp.concentration_at(0.2);
        assert!((c - 0.01).abs() < EPS);
    }
    #[test]
    fn test_rouse_profile_depth_averaged_nonnegative() {
        let rp = RouseProfile::new(2.0, 0.2, 0.01, 0.001, 0.05);
        let avg = rp.depth_averaged_concentration(20);
        assert!(avg >= 0.0);
    }
    #[test]
    fn test_shields_is_mobile() {
        let sp = ShieldsParam::new(1e-3, RHO_SEDIMENT, RHO_WATER, 0.047);
        let critical_tau = sp.theta_cr * (RHO_SEDIMENT - RHO_WATER) * G * 1e-3;
        assert!(!sp.is_mobile(critical_tau * 0.5));
        assert!(sp.is_mobile(critical_tau * 2.0));
    }
    #[test]
    fn test_shields_excess_zero() {
        let sp = ShieldsParam::new(1e-3, RHO_SEDIMENT, RHO_WATER, 0.047);
        let excess = sp.excess_shields(0.0);
        assert_eq!(excess, 0.0);
    }
    #[test]
    fn test_bedload_transport_zero_below_critical() {
        let bl = BedLoadTransport::new(1e-3, RHO_SEDIMENT, RHO_WATER, 0.047);
        let critical_tau = 0.047 * (RHO_SEDIMENT - RHO_WATER) * G * 1e-3;
        let rate = bl.transport_rate(critical_tau * 0.5);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_bedload_transport_vector_direction() {
        let bl = BedLoadTransport::new(1e-3, RHO_SEDIMENT, RHO_WATER, 0.047);
        let tau_b_vec = [2.0, 0.0, 0.0];
        let q_vec = bl.transport_vector(tau_b_vec);
        assert!(q_vec[0] >= 0.0);
        assert_eq!(q_vec[1], 0.0);
    }
    #[test]
    fn test_suspended_load_erosion_zero_below_critical() {
        let sl = SuspendedLoad::new(0.7, 1e-4, 0.001, 1e-3);
        let e = sl.erosion_flux(0.01, 0.05, 1.0);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn test_suspended_load_deposition_positive() {
        let sl = SuspendedLoad::new(0.7, 1e-4, 0.001, 1e-3);
        let d = sl.deposition_flux(0.1);
        assert!((d - 0.0001).abs() < 1e-12);
    }
    #[test]
    fn test_erosion_model_rate_zero() {
        let em = ErosionModel::new(0.1, 1e-3, 1.0);
        assert_eq!(em.erosion_rate(0.05), 0.0);
    }
    #[test]
    fn test_erosion_model_rate_positive() {
        let em = ErosionModel::new(0.1, 1e-3, 1.0);
        assert!(em.erosion_rate(0.5) > 0.0);
    }
    #[test]
    fn test_settling_velocity_stokes_method() {
        let sv = SettlingVelocity::new(1e-3, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let vs = sv.stokes();
        assert!(vs > 0.0);
    }
    #[test]
    fn test_settling_velocity_rubey_positive() {
        let sv = SettlingVelocity::new(1e-3, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let vr = sv.rubey();
        assert!(vr > 0.0);
    }
    #[test]
    fn test_settling_velocity_non_spherical() {
        let sv = SettlingVelocity::new(1e-3, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let vs = sv.stokes();
        let vns = sv.non_spherical(0.64);
        assert!(vns < vs);
    }
    #[test]
    fn test_turbulent_diffusion_diffusivity() {
        let td = TurbulentDiffusion::new(0.7, 1e-4);
        let d = td.diffusivity();
        assert!((d - 1e-4 / 0.7).abs() < 1e-16);
    }
    #[test]
    fn test_turbulent_diffusion_flux_opposes_gradient() {
        let td = TurbulentDiffusion::new(0.7, 1e-4);
        let grad_c = [1.0, 0.0, 0.0];
        let flux = td.diffusion_flux(grad_c);
        assert!(flux[0] < 0.0);
    }
    #[test]
    fn test_bed_evolution_update() {
        let bed = vec![0.0, 0.0, 0.0];
        let dx = vec![1.0, 1.0, 1.0];
        let mut be = BedEvolution::new(bed, dx, 0.4, 30.0);
        let q_b = vec![0.0, 0.1, 0.0];
        be.update(&q_b, 1.0);
        assert!(be.bed_elevation[1] < 0.0);
    }
    #[test]
    fn test_bed_evolution_avalanche() {
        let bed = vec![0.0, 10.0, 0.0];
        let dx = vec![1.0, 1.0, 1.0];
        let mut be = BedEvolution::new(bed, dx, 0.4, 30.0);
        be.avalanche();
        let slope = (be.bed_elevation[1] - be.bed_elevation[0]).abs();
        assert!(slope < 10.0);
    }
    #[test]
    fn test_coastal_morphology_longshore_transport() {
        let cm = CoastalMorphology::new(0.39, 0.05, vec![0.0; 5], 100.0);
        let q = cm.longshore_transport(1.0);
        assert!(q > 0.0);
    }
    #[test]
    fn test_coastal_morphology_longshore_zero_wave() {
        let cm = CoastalMorphology::new(0.39, 0.05, vec![0.0; 5], 100.0);
        let q = cm.longshore_transport(0.0);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_pi_imported() {
        // PI is used in computation; verify the import is accessible by using it
        let _ = PI;
    }
    #[test]
    fn test_math_helpers() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let diff = sub3(b, a);
        let l = len3(diff);
        let expected = 3.0 * 3.0_f64.sqrt();
        assert!(
            (l - expected).abs() < 1e-10,
            "l={}, expected={}",
            l,
            expected
        );
        let s = add3(a, b);
        assert!((s[0] - 5.0).abs() < EPS);
        let sc = scale3(a, 3.0);
        assert!((sc[2] - 9.0).abs() < EPS);
    }
}
pub(super) fn dot3_sed(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cubic_kernel_grad_sed(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = dot3_sed(r_ij, r_ij).sqrt();
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
/// Approximate inverse error function (Abramowitz & Stegun).
pub(super) fn erfinv_approx(x: f64) -> f64 {
    let a = 0.147_f64;
    let ln_term = (1.0 - x * x).ln();
    let term1 = 2.0 / (std::f64::consts::PI * a) + ln_term / 2.0;
    let sqrt_val = (term1 * term1 - ln_term / a).max(0.0).sqrt();
    (sqrt_val - term1).sqrt()
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::sediment::types::*;
    use crate::sediment::types_advanced::*;
    #[test]
    fn test_engelund_hansen_zero_velocity() {
        let eh = EngelundHansen::new(1e-3, RHO_SEDIMENT, RHO_WATER, 40.0);
        let q = eh.total_transport(0.0);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_engelund_hansen_positive() {
        let eh = EngelundHansen::new(1e-3, RHO_SEDIMENT, RHO_WATER, 40.0);
        let q = eh.total_transport(1.0);
        assert!(q > 0.0);
    }
    #[test]
    fn test_engelund_hansen_shields() {
        let eh = EngelundHansen::new(1e-3, RHO_SEDIMENT, RHO_WATER, 40.0);
        let theta = eh.shields_from_velocity(1.0);
        assert!(theta > 0.0);
    }
    #[test]
    fn test_flocculation_particles_per_floc() {
        let fl = FlocculationModel::new(1e-6);
        let n = fl.particles_per_floc(10e-6);
        assert!((n - 100.0).abs() < 1.0);
    }
    #[test]
    fn test_flocculation_floc_density_approaches_primary() {
        let fl = FlocculationModel::new(1e-6);
        let rho = fl.floc_density(1e-6, 2650.0, 1000.0);
        assert!((rho - 2650.0).abs() < 1.0);
    }
    #[test]
    fn test_flocculation_equilibrium_diameter_decreases_with_shear() {
        let mut fl = FlocculationModel::new(1e-6);
        fl.aggregation_rate = 1e-6;
        let d1 = fl.equilibrium_diameter(0.1);
        let d2 = fl.equilibrium_diameter(10.0);
        assert!(
            d2 < d1,
            "higher shear should give smaller flocs: d1={d1}, d2={d2}"
        );
    }
    #[test]
    fn test_kynch_flux_zero_at_zero_phi() {
        let ks = KynchSettling::new(0.001, 0.65, 4.65);
        let f = ks.flux(0.0);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn test_kynch_hindered_velocity_decreases() {
        let ks = KynchSettling::new(0.001, 0.65, 4.65);
        let v1 = ks.hindered_velocity(0.1);
        let v2 = ks.hindered_velocity(0.5);
        assert!(v1 > v2);
    }
    #[test]
    fn test_kynch_critical_concentration_positive() {
        let ks = KynchSettling::new(0.001, 0.65, 4.65);
        let phi_c = ks.critical_concentration();
        assert!(phi_c > 0.0 && phi_c < 0.65);
    }
    #[test]
    fn test_turbidity_current_reduced_gravity_positive() {
        let tc =
            TurbidityCurrentSph::new(2.0, 0.01, 1.0, 0.01, 0.004, 0.001, 0.001, 1000.0, 2650.0);
        let gp = tc.reduced_gravity();
        assert!(gp > 0.0);
    }
    #[test]
    fn test_turbidity_current_froude_positive() {
        let tc =
            TurbidityCurrentSph::new(2.0, 0.01, 1.0, 0.01, 0.004, 0.001, 0.001, 1000.0, 2650.0);
        let fr = tc.froude_number();
        assert!(fr > 0.0);
    }
    #[test]
    fn test_turbidity_current_advance() {
        let mut tc =
            TurbidityCurrentSph::new(2.0, 0.01, 1.0, 0.001, 0.004, 0.001, 0.001, 1000.0, 2650.0);
        let u_init = tc.velocity;
        tc.advance(0.0, 1.0);
        let _changed = tc.velocity != u_init;
    }
    #[test]
    fn test_scour_model_csu_pier_positive() {
        let sc = ScourModel::new_pier(1.0, 3.0);
        let ys = sc.csu_pier_scour(0.3);
        assert!(ys > 0.0);
    }
    #[test]
    fn test_scour_model_froehlich_abutment() {
        let sc = ScourModel::new_pier(1.0, 3.0);
        let ys = sc.froehlich_abutment_scour(10.0, 0.3);
        assert!(ys > 0.0);
    }
    #[test]
    fn test_dune_migration_velocity_positive() {
        let dm = DuneMigration::new(0.3, 5.0, 1e-4, 0.4, 0.1);
        let c = dm.migration_velocity();
        assert!(c > 0.0);
    }
    #[test]
    fn test_dune_van_rijn_height_zero_at_critical() {
        let dm = DuneMigration::new(0.3, 5.0, 1e-4, 0.4, 0.1);
        let h = dm.van_rijn_height(2.0, 2e-4, 0.047, 0.047);
        assert_eq!(h, 0.0);
    }
    #[test]
    fn test_dune_form_drag_positive() {
        let dm = DuneMigration::new(0.3, 5.0, 1e-4, 0.4, 0.1);
        let cf = dm.form_drag_coeff(2.0);
        assert!(cf > 0.0);
    }
    #[test]
    fn test_shoreline_evolution_advance_diffuses() {
        let pos = vec![0.0, 10.0, 0.0, 0.0, 0.0];
        let mut se = ShorelineEvolution::new(pos, 100.0, 100.0, 5.0, 0.39);
        let dt = se.max_dt() * 0.4;
        se.advance(dt);
        assert!(se.position[1] < 10.0);
    }
    #[test]
    fn test_shoreline_max_dt_positive() {
        let se = ShorelineEvolution::new(vec![0.0; 5], 100.0, 100.0, 5.0, 0.39);
        let dt = se.max_dt();
        assert!(dt > 0.0);
    }
    #[test]
    fn test_exner_equation_update() {
        let eta = vec![0.0, 0.0, 0.0, 0.0];
        let mut ex = ExnerEquation::new(eta, 1.0, 0.4);
        let q_b = vec![0.0, 1e-4, 1e-4, 0.0];
        ex.update(&q_b, 100.0);
        assert!(ex.eta[1] < 0.0);
    }
    #[test]
    fn test_exner_max_scour_positive() {
        let eta_0 = vec![0.0, 0.0, 0.0, 0.0];
        let eta = vec![0.0, -0.1, -0.2, 0.0];
        let ex = ExnerEquation {
            eta,
            dx: 1.0,
            porosity: 0.4,
        };
        let scour = ex.max_scour(&eta_0);
        assert!((scour - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_sph_sediment_advection_diffusion_rate_zero_uniform() {
        let sas = SphSedimentAdvection::new(0.1, 0.001, 1e-5);
        let rate = sas.diffusion_rate(0.05, &[([0.05, 0.0, 0.0], 0.05, 1.0, 1000.0)]);
        assert!(rate.abs() < 1e-10);
    }
    #[test]
    fn test_morphodynamic_feedback_effective_flux() {
        let mf = MorphodynamicFeedback::new(10.0, 0.05, 0.02);
        let q_eff = mf.effective_flux(1e-4);
        assert!((q_eff - 1e-3).abs() < 1e-15);
    }
    #[test]
    fn test_morphodynamic_feedback_manning_n() {
        let mf = MorphodynamicFeedback::new(10.0, 0.001, 0.02);
        let n = mf.manning_n();
        assert!(n > 0.0 && n < 0.1);
    }
    #[test]
    fn test_van_rijn_bedload_positive() {
        let vr = VanRijnBedLoad::new(2e-4, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let q = vr.transport_rate(0.05, 0.01);
        assert!(q > 0.0);
    }
    #[test]
    fn test_van_rijn_bedload_zero_below_critical() {
        let vr = VanRijnBedLoad::new(2e-4, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let q = vr.transport_rate(0.005, 0.01);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_van_rijn_dimensionless_grain_size() {
        let vr = VanRijnBedLoad::new(2e-4, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let d_star = vr.dimensionless_grain_size();
        assert!(d_star > 1.0);
    }
    #[test]
    fn test_sph_bedload_flux_vector_zero_below_critical() {
        let sbf = SphBedLoadFlux::new(1e-3, 0.047, RHO_SEDIMENT, RHO_WATER);
        let f = sbf.flux_vector([0.1, 0.0, 0.0]);
        assert_eq!(len3(f), 0.0);
    }
    #[test]
    fn test_suspended_sediment_particle_rouse_number() {
        let p = SuspendedSedimentParticle::new([0.0; 3], [0.0; 3], 1000.0, 1e-3, 0.1);
        let ro = p.rouse_number(0.05);
        assert!(ro > 0.0);
    }
    #[test]
    fn test_suspended_sediment_transport_mode_bed_load() {
        let p = SuspendedSedimentParticle::new([0.0; 3], [0.0; 3], 1000.0, 1e-3, 0.1);
        let mode = p.transport_mode(0.001);
        assert_eq!(mode, SedimentTransportMode::BedLoad);
    }
    #[test]
    fn test_grain_size_distribution_d84_greater_d16() {
        let gsd = GrainSizeDistribution::new(1e-3, 2.0);
        let d16 = gsd.d16();
        let d84 = gsd.d84();
        assert!(d84 > d16);
    }
    #[test]
    fn test_grain_size_distribution_sorting_coeff() {
        let gsd = GrainSizeDistribution::new(1e-3, 2.0);
        let psi = gsd.sorting_coeff();
        assert!(psi > 1.0);
    }
    #[test]
    fn test_time_evolving_scour_starts_zero() {
        let tes = TimeEvolvingScour::new(3.0, 3600.0);
        assert_eq!(tes.depth, 0.0);
    }
    #[test]
    fn test_time_evolving_scour_approaches_equilibrium() {
        let tes = TimeEvolvingScour::new(3.0, 3600.0);
        let d = tes.scour_at_time(10.0 * 3600.0);
        assert!(d > 2.9);
    }
    #[test]
    fn test_time_evolving_scour_advance() {
        let mut tes = TimeEvolvingScour::new(3.0, 3600.0);
        tes.advance(3600.0);
        assert!(tes.depth > 0.0);
    }
    #[test]
    fn test_flocculation_settling_velocity_positive() {
        let fl = FlocculationModel::new(1e-6);
        let mu = RHO_WATER * NU_WATER;
        let vs = fl.floc_settling_velocity(10e-6, 2650.0, 1000.0, mu);
        assert!(vs > 0.0);
    }
    #[test]
    fn test_turbidity_current_bed_shear_positive() {
        let tc =
            TurbidityCurrentSph::new(2.0, 0.01, 1.0, 0.001, 0.004, 0.001, 0.001, 1000.0, 2650.0);
        let tau = tc.bed_shear_stress();
        assert!(tau > 0.0);
    }
    #[test]
    fn test_scour_model_melville_pier_positive() {
        let sc = ScourModel::new_pier(1.0, 3.0);
        let ys = sc.melville_pier_scour(0.3, 2e-4);
        assert!(ys > 0.0);
    }
    #[test]
    fn test_shoreline_longshore_drift_positive() {
        let se = ShorelineEvolution::new(vec![0.0; 5], 100.0, 100.0, 5.0, 0.39);
        let q = se.longshore_drift(1.0, 0.1);
        assert!(q >= 0.0);
    }
}
#[cfg(test)]
mod tests_part3 {
    use super::*;
    use crate::sediment::types::*;
    use crate::sediment::types_advanced::*;
    #[test]
    fn test_mixture_density_pure_fluid() {
        let md = MixtureDensity::new(2650.0, 1000.0);
        let rho = md.density(0.0);
        assert!((rho - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_mixture_density_pure_sediment() {
        let md = MixtureDensity::new(2650.0, 1000.0);
        let rho = md.density(1.0);
        assert!((rho - 2650.0).abs() < 1e-10);
    }
    #[test]
    fn test_mixture_density_multi_fraction() {
        let md = MixtureDensity::new(2650.0, 1000.0);
        let fractions = vec![(0.05, 2650.0), (0.05, 2500.0)];
        let rho = md.multi_fraction_density(&fractions);
        assert!(rho > 1000.0 && rho < 2650.0);
    }
    #[test]
    fn test_mixture_density_submerged_specific_gravity() {
        let md = MixtureDensity::new(2650.0, 1000.0);
        let s_minus1 = md.submerged_specific_gravity();
        assert!((s_minus1 - 1.65).abs() < 1e-10);
    }
    #[test]
    fn test_multi_regime_settling_ferguson_church_positive() {
        let ms = MultiRegimeSettling::new(1e-3, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let ws = ms.ferguson_church();
        assert!(ws > 0.0, "ws={}", ws);
    }
    #[test]
    fn test_multi_regime_settling_stokes_positive() {
        let ms = MultiRegimeSettling::new(1e-4, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let ws = ms.stokes();
        assert!(ws > 0.0);
    }
    #[test]
    fn test_multi_regime_settling_newton_positive() {
        let ms = MultiRegimeSettling::new(1e-2, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let ws = ms.newton_regime();
        assert!(ws > 0.0);
    }
    #[test]
    fn test_multi_regime_settling_small_grain_uses_stokes() {
        let ms = MultiRegimeSettling::new(1e-5, RHO_SEDIMENT, RHO_WATER, NU_WATER);
        let re = ms.particle_reynolds_number();
        assert!(re < 0.5, "re={}", re);
    }
    #[test]
    fn test_armoring_model_surface_mean_diameter() {
        let am = ArmoringModel::new(vec![1e-3, 2e-3, 4e-3], 0.05);
        let dm = am.surface_mean_diameter();
        assert!(dm > 0.0);
    }
    #[test]
    fn test_armoring_model_hiding_factor_unity_at_dm() {
        let mut am = ArmoringModel::new(vec![2e-3], 0.05);
        am.surface_fractions = vec![1.0];
        let xi = am.hiding_factor(0);
        assert!((xi - 1.0).abs() < 1e-6, "xi={}", xi);
    }
    #[test]
    fn test_armoring_model_fractional_transport_nonnegative() {
        let am = ArmoringModel::new(vec![1e-3, 2e-3], 0.05);
        let q = am.fractional_transport(0, 5.0, 0.047);
        assert!(q >= 0.0);
    }
    #[test]
    fn test_armoring_model_update_surface_normalised() {
        let mut am = ArmoringModel::new(vec![1e-3, 2e-3, 4e-3], 0.05);
        am.update_surface(-0.01);
        let total: f64 = am.surface_fractions.iter().sum();
        assert!((total - 1.0).abs() < 1e-10, "total={}", total);
    }
    #[test]
    fn test_cohesive_sediment_deposition_zero_above_tau_d() {
        let cs = CohesiveSediment::new(0.1, 0.3, 1e-3, 0.001);
        let d = cs.deposition_flux(0.15, 0.01);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_cohesive_sediment_erosion_zero_below_tau_e() {
        let cs = CohesiveSediment::new(0.1, 0.3, 1e-3, 0.001);
        let e = cs.erosion_flux(0.1);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn test_cohesive_sediment_erosion_positive_above_tau_e() {
        let cs = CohesiveSediment::new(0.1, 0.3, 1e-3, 0.001);
        let e = cs.erosion_flux(0.5);
        assert!(e > 0.0);
    }
    #[test]
    fn test_cohesive_sediment_update_bed_increases() {
        let mut cs = CohesiveSediment::new(0.1, 0.3, 1e-3, 0.001);
        cs.update_bed(0.05, 0.01, 10.0);
        assert!(cs.bed_mass > 0.0);
    }
    #[test]
    fn test_sediment_concentration_1d_advances() {
        let mut sc = SedimentConcentration1D::new(5, 1.0, 0.5, 1e-4, 0.001, 2.0);
        sc.concentration[2] = 0.1;
        let dt = sc.max_dt() * 0.4;
        let zero = vec![0.0; 5];
        sc.advance(dt, &zero, &zero);
        for &c in &sc.concentration {
            assert!(c >= 0.0);
        }
    }
    #[test]
    fn test_sediment_concentration_1d_total_mass() {
        let sc = SedimentConcentration1D::new(5, 1.0, 0.5, 1e-4, 0.001, 2.0);
        let m = sc.total_mass();
        assert_eq!(m, 0.0);
    }
    #[test]
    fn test_multi_fraction_exner_update() {
        let eta = vec![0.0, 0.0, 0.0, 0.0];
        let mut mfe = MultiFractionExner::new(eta, 1.0, 0.4, 2);
        let q1 = vec![0.0, 1e-4, 1e-4, 0.0];
        let q2 = vec![0.0, 5e-5, 5e-5, 0.0];
        mfe.update(&[q1, q2], 100.0);
        assert!(mfe.eta[1] < 0.0, "eta[1]={}", mfe.eta[1]);
    }
    #[test]
    fn test_multi_fraction_exner_volume_change_zero_initially() {
        let eta = vec![0.0; 5];
        let mfe = MultiFractionExner::new(eta.clone(), 1.0, 0.4, 2);
        let vc = mfe.volume_change(&eta);
        assert_eq!(vc, 0.0);
    }
    #[test]
    fn test_vertical_diffusion_mean_concentration_zero_init() {
        let vd = VerticalDiffusion::new(10, 0.01, 1e-5, 0.001);
        assert_eq!(vd.mean_concentration(), 0.0);
    }
    #[test]
    fn test_vertical_diffusion_max_dt_positive() {
        let vd = VerticalDiffusion::new(10, 0.01, 1e-5, 0.001);
        let dt = vd.max_dt();
        assert!(dt > 0.0);
    }
    #[test]
    fn test_morphodynamic_timestep_max_dt_positive() {
        let mt = MorphodynamicTimestep::new(0.1, 1.0, 5.0);
        let eta = vec![1.0; 5];
        let q_b = vec![1e-4; 5];
        let dt = mt.max_dt(&eta, &q_b, 0.4);
        assert!(dt > 0.0);
    }
    #[test]
    fn test_shields_slope_horizontal_equals_base() {
        let ss = ShieldsSlope::new(0.047, 32.0);
        let theta_cr = ss.critical_shields_downslope(0.0);
        assert!((theta_cr - 0.047).abs() < 1e-10, "theta_cr={}", theta_cr);
    }
    #[test]
    fn test_shields_slope_downslope_less_than_base() {
        let ss = ShieldsSlope::new(0.047, 32.0);
        let theta_cr_slope = ss.critical_shields_downslope(0.1);
        assert!(theta_cr_slope < 0.047, "theta_cr={}", theta_cr_slope);
    }
}
/// Settling velocity using Stokes' law (Richardson-Zaki form for small grains).
///
/// v_s = (ρ_s − ρ_f) g d² / (18 μ),  μ = ρ_f * ν
///
/// - `d`: grain diameter (m)
/// - `rho_s`: sediment density (kg/m³)
/// - `rho_f`: fluid density (kg/m³)
/// - `nu`: fluid kinematic viscosity (m²/s)
/// - `g`: gravitational acceleration (m/s²)
pub fn settling_velocity(d: f64, rho_s: f64, rho_f: f64, nu: f64, g: f64) -> f64 {
    let mu = rho_f * nu;
    if mu < 1e-300 || d < 0.0 {
        return 0.0;
    }
    (rho_s - rho_f) * g * d * d / (18.0 * mu)
}
/// Shields parameter θ = τ_b / ((ρ_s − ρ_f) g d).
///
/// - `tau_b`: bed shear stress (Pa)
/// - `rho_s`: sediment density (kg/m³)
/// - `rho_f`: fluid density (kg/m³)
/// - `g`: gravitational acceleration (m/s²)
/// - `d`: grain diameter (m)
pub fn shields_parameter(tau_b: f64, rho_s: f64, rho_f: f64, g: f64, d: f64) -> f64 {
    let denom = (rho_s - rho_f) * g * d;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    tau_b / denom
}
/// Soulsby critical Shields parameter as a function of dimensionless grain size d*.
///
/// θ_cr = 0.3 / (1 + 1.2 d*) + 0.055 * (1 − exp(−0.02 d*))
pub fn critical_shields(d_star: f64) -> f64 {
    if d_star <= 0.0 {
        return 0.3;
    }
    0.3 / (1.0 + 1.2 * d_star) + 0.055 * (1.0 - (-0.02 * d_star).exp())
}
/// Dimensionless grain size d* = d * ((ρ_s/ρ_f − 1) g / ν²)^{1/3}.
///
/// - `d`: grain diameter (m)
/// - `rho_s`: sediment density (kg/m³)
/// - `rho_f`: fluid density (kg/m³)
/// - `nu`: kinematic viscosity (m²/s)
/// - `g`: gravitational acceleration (m/s²)
pub fn dimensionless_grain(d: f64, rho_s: f64, rho_f: f64, nu: f64, g: f64) -> f64 {
    if rho_f.abs() < 1e-300 || nu.abs() < 1e-300 {
        return 0.0;
    }
    let s = rho_s / rho_f;
    d * ((s - 1.0) * g / (nu * nu)).powf(1.0 / 3.0)
}
/// Meyer-Peter Müller bed-load transport rate (m²/s).
///
/// q_b = 8 * (θ − θ_cr)^{1.5} * sqrt((s−1) g d³),  returns 0 when θ ≤ θ_cr.
///
/// - `theta`: Shields parameter
/// - `theta_cr`: critical Shields parameter
/// - `d`: grain diameter (m)
/// - `rho_s`: sediment density (kg/m³)
/// - `rho_f`: fluid density (kg/m³)
/// - `g`: gravitational acceleration (m/s²)
pub fn bedload_transport_rate(
    theta: f64,
    theta_cr: f64,
    d: f64,
    rho_s: f64,
    rho_f: f64,
    g: f64,
) -> f64 {
    if theta <= theta_cr {
        return 0.0;
    }
    let s = rho_s / rho_f;
    let excess = (theta - theta_cr).powf(1.5);
    8.0 * excess * ((s - 1.0) * g * d * d * d).sqrt()
}
/// Rouse suspended sediment concentration profile.
///
/// C(z) = C₀ * ((a/z) * ((H−z)/(H−a)))^Z,  where a = c_ref height, Z = Rouse number.
///
/// - `c0`: reference concentration at reference height
/// - `z`: elevation above bed (m)
/// - `h`: total water depth (m)
/// - `ws`: settling velocity (m/s)
/// - `kappa`: von Kármán constant (≈ 0.41)
/// - `u_star`: shear velocity (m/s)
pub fn suspended_sediment_concentration(
    c0: f64,
    z: f64,
    h: f64,
    ws: f64,
    kappa: f64,
    u_star: f64,
) -> f64 {
    let a = h * 0.05;
    let ro = rouse_number(ws, kappa, u_star);
    if z <= 0.0 || z >= h || a >= h {
        return 0.0;
    }
    let ratio = (a / z.max(1e-300)) * ((h - z) / (h - a).max(1e-300));
    if ratio <= 0.0 {
        return 0.0;
    }
    c0 * ratio.powf(ro)
}
/// Rouse number Z = w_s / (κ u*).
///
/// - `ws`: settling velocity (m/s)
/// - `kappa`: von Kármán constant
/// - `u_star`: bed shear velocity (m/s)
pub fn rouse_number(ws: f64, kappa: f64, u_star: f64) -> f64 {
    let denom = kappa * u_star;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    ws / denom
}
/// Empirical bedform (dune) height prediction.
///
/// Δ/d ≈ 0.11 * (θ/θ_cr − 1)^{0.3}  for θ > θ_cr.  Returns 0 otherwise.
///
/// - `shields`: current Shields parameter θ
/// - `d`: grain diameter (m)
pub fn bedform_height(shields: f64, d: f64) -> f64 {
    let theta_cr = 0.047;
    if shields <= theta_cr {
        return 0.0;
    }
    let excess = shields / theta_cr - 1.0;
    0.11 * d * excess.powf(0.3)
}
#[cfg(test)]
mod sediment_new_tests {
    use super::*;
    use crate::sediment::types::*;

    #[test]
    fn test_settling_velocity_positive() {
        let vs = settling_velocity(1e-3, 2650.0, 1000.0, 1e-6, 9.81);
        assert!(vs > 0.0);
    }
    #[test]
    fn test_settling_velocity_zero_diameter() {
        let vs = settling_velocity(0.0, 2650.0, 1000.0, 1e-6, 9.81);
        assert_eq!(vs, 0.0);
    }
    #[test]
    fn test_settling_velocity_proportional_d2() {
        let v1 = settling_velocity(1e-3, 2650.0, 1000.0, 1e-6, 9.81);
        let v2 = settling_velocity(2e-3, 2650.0, 1000.0, 1e-6, 9.81);
        assert!((v2 - 4.0 * v1).abs() < 1e-10 * v1);
    }
    #[test]
    fn test_settling_velocity_zero_nu() {
        let vs = settling_velocity(1e-3, 2650.0, 1000.0, 0.0, 9.81);
        assert_eq!(vs, 0.0);
    }
    #[test]
    fn test_shields_parameter_formula() {
        let theta = shields_parameter(1.0, 2650.0, 1000.0, 9.81, 1e-3);
        let expected = 1.0 / ((2650.0 - 1000.0) * 9.81 * 1e-3);
        assert!((theta - expected).abs() < 1e-12);
    }
    #[test]
    fn test_shields_parameter_zero_stress() {
        let theta = shields_parameter(0.0, 2650.0, 1000.0, 9.81, 1e-3);
        assert_eq!(theta, 0.0);
    }
    #[test]
    fn test_shields_parameter_zero_grain() {
        let theta = shields_parameter(1.0, 2650.0, 1000.0, 9.81, 0.0);
        assert_eq!(theta, 0.0);
    }
    #[test]
    fn test_critical_shields_large_dstar() {
        let tc_small = critical_shields(1.0);
        let tc_large = critical_shields(1000.0);
        assert!(tc_large < tc_small);
    }
    #[test]
    fn test_critical_shields_zero() {
        let tc = critical_shields(0.0);
        assert!(tc > 0.0);
    }
    #[test]
    fn test_critical_shields_positive() {
        let tc = critical_shields(5.0);
        assert!(tc > 0.0);
    }
    #[test]
    fn test_dimensionless_grain_positive() {
        let dstar = dimensionless_grain(1e-3, 2650.0, 1000.0, 1e-6, 9.81);
        assert!(dstar > 0.0);
    }
    #[test]
    fn test_dimensionless_grain_zero_nu() {
        let dstar = dimensionless_grain(1e-3, 2650.0, 1000.0, 0.0, 9.81);
        assert_eq!(dstar, 0.0);
    }
    #[test]
    fn test_bedload_transport_rate_zero_at_critical() {
        let q = bedload_transport_rate(0.047, 0.047, 1e-3, 2650.0, 1000.0, 9.81);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_bedload_transport_rate_zero_below_critical() {
        let q = bedload_transport_rate(0.02, 0.047, 1e-3, 2650.0, 1000.0, 9.81);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_bedload_transport_rate_positive_above_critical() {
        let q = bedload_transport_rate(0.1, 0.047, 1e-3, 2650.0, 1000.0, 9.81);
        assert!(q > 0.0);
    }
    #[test]
    fn test_bedload_transport_rate_increases_with_shields() {
        let q1 = bedload_transport_rate(0.1, 0.047, 1e-3, 2650.0, 1000.0, 9.81);
        let q2 = bedload_transport_rate(0.2, 0.047, 1e-3, 2650.0, 1000.0, 9.81);
        assert!(q2 > q1);
    }
    #[test]
    fn test_rouse_number_formula() {
        let ro = rouse_number(0.02, 0.41, 0.1);
        assert!((ro - 0.02 / (0.41 * 0.1)).abs() < 1e-12);
    }
    #[test]
    fn test_rouse_number_zero_shear() {
        let ro = rouse_number(0.02, 0.41, 0.0);
        assert_eq!(ro, 0.0);
    }
    #[test]
    fn test_suspended_sediment_concentration_decreases_with_height() {
        let c1 = suspended_sediment_concentration(1.0, 0.1, 1.0, 0.01, 0.41, 0.05);
        let c2 = suspended_sediment_concentration(1.0, 0.4, 1.0, 0.01, 0.41, 0.05);
        assert!(c1 >= c2 || c1 == 0.0);
    }
    #[test]
    fn test_suspended_sediment_concentration_at_surface_zero() {
        let c = suspended_sediment_concentration(1.0, 1.0, 1.0, 0.01, 0.41, 0.05);
        assert_eq!(c, 0.0);
    }
    #[test]
    fn test_sediment_bed_new() {
        let bed = SedimentBed::new(5, 0.1);
        assert_eq!(bed.z_bed.len(), 5);
        assert_eq!(bed.x.len(), 5);
    }
    #[test]
    fn test_sediment_bed_erode_reduces_depth() {
        let mut bed = SedimentBed::new(3, 0.1);
        bed.z_bed = vec![1.0, 1.0, 1.0];
        let rates = vec![0.5, 0.5, 0.5];
        bed.erode(&rates, 1.0);
        for z in &bed.z_bed {
            assert!((*z - 0.5).abs() < 1e-12);
        }
    }
    #[test]
    fn test_sediment_bed_deposit_increases_depth() {
        let mut bed = SedimentBed::new(3, 0.1);
        let rates = vec![0.2, 0.2, 0.2];
        bed.deposit(&rates, 1.0);
        for z in &bed.z_bed {
            assert!((*z - 0.2).abs() < 1e-12);
        }
    }
    #[test]
    fn test_sediment_bed_total_volume() {
        let mut bed = SedimentBed::new(4, 0.5);
        bed.z_bed = vec![1.0; 4];
        assert!((bed.total_volume() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_sediment_bed_depth_at() {
        let mut bed = SedimentBed::new(3, 0.1);
        bed.z_bed[1] = 2.5;
        assert!((bed.depth_at(1) - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_sediment_bed_depth_at_out_of_bounds() {
        let bed = SedimentBed::new(3, 0.1);
        assert_eq!(bed.depth_at(100), 0.0);
    }
    #[test]
    fn test_bedform_height_zero_at_critical() {
        let h = bedform_height(0.047, 1e-3);
        assert_eq!(h, 0.0);
    }
    #[test]
    fn test_bedform_height_positive_above_critical() {
        let h = bedform_height(0.1, 1e-3);
        assert!(h > 0.0);
    }
    #[test]
    fn test_bedform_height_scales_with_d() {
        let h1 = bedform_height(0.1, 1e-3);
        let h2 = bedform_height(0.1, 2e-3);
        assert!((h2 - 2.0 * h1).abs() < 1e-12);
    }
}
