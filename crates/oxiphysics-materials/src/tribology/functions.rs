//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Complementary error function approximation.
pub(super) fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}
/// Coulomb friction force F = μ·N.
///
/// Returns the friction force magnitude \[N\] given a normal force `normal` \[N\]
/// and a coefficient of friction `mu`.
pub fn coulomb_friction_force(normal: f64, mu: f64) -> f64 {
    mu * normal
}
/// Stribeck friction coefficient as a function of sliding velocity.
///
/// Models the transition from static (`mu_static`) to kinetic (`mu_kinetic`)
/// friction using an exponential decay:
/// μ(v) = mu_kinetic + (mu_static − mu_kinetic) · exp(−|v| / v_stribeck)
pub fn stribeck_friction(velocity: f64, mu_static: f64, mu_kinetic: f64, v_stribeck: f64) -> f64 {
    let v_abs = velocity.abs();
    let vs = v_stribeck.max(1e-30);
    mu_kinetic + (mu_static - mu_kinetic) * (-v_abs / vs).exp()
}
/// Viscous (hydrodynamic) friction force using Newton's law of viscosity.
///
/// F = η · A · v / h, where `eta` is dynamic viscosity \[Pa·s\],
/// `area` is the shear area \[m²\], `velocity` is the sliding speed \[m/s\],
/// and `gap` is the fluid film thickness \[m\].
pub fn viscous_friction(velocity: f64, eta: f64, area: f64, gap: f64) -> f64 {
    if gap < 1e-30 {
        return 0.0;
    }
    eta * area * velocity / gap
}
/// Adhesive wear rate \[m³/m\] based on real contact area.
///
/// Q/s = k_adhesive · A_real / H, where `real_contact_area` \[m²\],
/// `hardness` \[Pa\], and `k_adhesive` is a dimensionless adhesive wear
/// coefficient.
pub fn adhesive_wear_rate(real_contact_area: f64, hardness: f64, k_adhesive: f64) -> f64 {
    if hardness < 1e-30 {
        return 0.0;
    }
    k_adhesive * real_contact_area / hardness
}
/// Abrasive wear factor combining particle size and hardness ratio effects.
///
/// Returns a dimensionless wear intensity factor = particle_size² · hardness_ratio.
/// Larger, harder particles produce more abrasion.
pub fn abrasive_wear_factor(particle_size: f64, hardness_ratio: f64) -> f64 {
    particle_size * particle_size * hardness_ratio
}
/// Fraction of load carried by asperity (boundary) contact as a function of
/// the lambda ratio Λ.
///
/// Returns a value in \[0, 1\] using the empirical Greenwood model:
/// f = exp(−0.6·Λ²)  clamped to \[0, 1\].
pub fn asperity_contact_fraction(lambda: f64) -> f64 {
    if lambda < 0.0 {
        return 1.0;
    }
    (-0.6 * lambda * lambda).exp()
}
/// Flash temperature rise at a sliding contact using Jaeger's formula.
///
/// ΔT = q / (k_thermal) · sqrt(a / (π · v))
///
/// where `heat_flux` q \[W/m²\], `k_thermal` \[W/(m·K)\], `v` is the sliding
/// speed \[m/s\], and `a` is the contact radius \[m\].
pub fn flash_temperature(heat_flux: f64, k_thermal: f64, v: f64, a: f64) -> f64 {
    if k_thermal < 1e-30 || v < 1e-30 {
        return 0.0;
    }
    heat_flux / k_thermal * (a / (PI * v)).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::tribology::AbrasiveWear;
    use crate::tribology::ArchardWear;
    use crate::tribology::AsperityModel;
    use crate::tribology::BearingLife;
    use crate::tribology::BowdenTaborFriction;
    use crate::tribology::BrakePad;
    use crate::tribology::CoatingMaterial;
    use crate::tribology::CylinderContact;
    use crate::tribology::DerjaguinMullerToporov;
    use crate::tribology::EhlFilm;
    use crate::tribology::ElastohydrodynamicFilm;
    use crate::tribology::FrettingFatigue;
    use crate::tribology::FrictionEvolution;
    use crate::tribology::FrictionOscillator;
    use crate::tribology::GearTooth;
    use crate::tribology::HamrockDowson;
    use crate::tribology::HertzContact;
    use crate::tribology::HertzContactParams;
    use crate::tribology::HertzElliptical;
    use crate::tribology::HolmWear;
    use crate::tribology::JohnsonKendallRoberts;
    use crate::tribology::LubricantViscosity;
    use crate::tribology::LubricationRegime;
    use crate::tribology::LubricationRegimeKind;
    use crate::tribology::ScuffingModel;
    use crate::tribology::StribeckCurve;
    use crate::tribology::SurfaceRoughness;
    use crate::tribology::ThermalWear;
    use crate::tribology::TribologicalState;
    use crate::tribology::WearModel;
    #[test]
    fn hertz_sphere_sphere_contact_radius_positive() {
        let hz = HertzContact::sphere_sphere(200e9, 0.3, 0.01, 200e9, 0.3, 0.01);
        let res = hz.compute(100.0);
        assert!(
            res.contact_radius > 0.0,
            "contact_radius={}",
            res.contact_radius
        );
    }
    #[test]
    fn hertz_sphere_flat_pressure_positive() {
        let hz = HertzContact::sphere_flat(200e9, 0.3, 0.01, 200e9, 0.3);
        let res = hz.compute(50.0);
        assert!(res.max_pressure > 0.0);
    }
    #[test]
    fn hertz_contact_radius_increases_with_load() {
        let hz = HertzContact::sphere_sphere(200e9, 0.3, 0.01, 200e9, 0.3, 0.01);
        let r1 = hz.compute(100.0).contact_radius;
        let r2 = hz.compute(1000.0).contact_radius;
        assert!(r2 > r1, "r2={r2}, r1={r1}");
    }
    #[test]
    fn hertz_force_for_radius_roundtrip() {
        let hz = HertzContact::sphere_sphere(200e9, 0.3, 0.01, 200e9, 0.3, 0.01);
        let f = 100.0;
        let a = hz.compute(f).contact_radius;
        let f2 = hz.force_for_radius(a);
        assert!((f2 - f).abs() / f < 0.01, "f={f}, f2={f2}");
    }
    #[test]
    fn hertz_cylinder_cylinder_creates_model() {
        let hz = HertzContact::cylinder_cylinder(200e9, 0.3, 0.01, 200e9, 0.3, 0.01);
        assert!(hz.effective_modulus > 0.0);
        assert!(hz.effective_radius > 0.0);
    }
    #[test]
    fn archard_wear_volume_linear_in_load() {
        let archard = ArchardWear::new(1e-4, 1e9);
        let v1 = archard.volume_loss(100.0, 1.0);
        let v2 = archard.volume_loss(200.0, 1.0);
        assert!((v2 / v1 - 2.0).abs() < 1e-10, "ratio={}", v2 / v1);
    }
    #[test]
    fn archard_wear_depth_finite() {
        let archard = ArchardWear::new(1e-4, 1e9);
        let d = archard.depth_loss(100.0, 1.0, 1e-4);
        assert!(d.is_finite() && d > 0.0, "d={d}");
    }
    #[test]
    fn archard_zero_area_returns_zero() {
        let archard = ArchardWear::new(1e-4, 1e9);
        let d = archard.depth_loss(100.0, 1.0, 0.0);
        assert!((d).abs() < 1e-30);
    }
    #[test]
    fn jkr_pull_off_force_negative() {
        let jkr = JohnsonKendallRoberts::new(1e9, 1e-6, 0.1);
        assert!(jkr.pull_off_force() < 0.0);
    }
    #[test]
    fn jkr_contact_radius_at_zero_force() {
        let jkr = JohnsonKendallRoberts::new(1e9, 1e-6, 0.1);
        let a = jkr.contact_radius(0.0);
        assert!(a >= 0.0, "a={a}");
    }
    #[test]
    fn jkr_contact_radius_increases_with_force() {
        let jkr = JohnsonKendallRoberts::new(1e9, 1e-6, 0.1);
        let a1 = jkr.contact_radius(0.0);
        let a2 = jkr.contact_radius(1e-4);
        assert!(a2 >= a1, "a2={a2}, a1={a1}");
    }
    #[test]
    fn dmt_pull_off_force_negative() {
        let dmt = DerjaguinMullerToporov::new(1e9, 1e-6, 0.1);
        assert!(dmt.pull_off_force() < 0.0);
    }
    #[test]
    fn dmt_pull_off_larger_than_jkr() {
        let jkr = JohnsonKendallRoberts::new(1e9, 1e-6, 0.1);
        let dmt = DerjaguinMullerToporov::new(1e9, 1e-6, 0.1);
        assert!(dmt.pull_off_force().abs() > jkr.pull_off_force().abs());
    }
    #[test]
    fn dmt_tabor_parameter_finite() {
        let dmt = DerjaguinMullerToporov::new(1e9, 1e-6, 0.1);
        let mu = dmt.tabor_parameter(0.3e-9);
        assert!(mu.is_finite() && mu > 0.0, "tabor={mu}");
    }
    #[test]
    fn asperity_bearing_area_at_zero_separation() {
        let model = AsperityModel::new(1e12, 1e-6, 1e-6, 1e11);
        let ba = model.bearing_area(0.0);
        assert!((ba - 0.5).abs() < 0.01, "ba={ba}");
    }
    #[test]
    fn asperity_bearing_area_decreases_with_separation() {
        let model = AsperityModel::new(1e12, 1e-6, 1e-6, 1e11);
        let ba1 = model.bearing_area(0.0);
        let ba2 = model.bearing_area(3.0e-6);
        assert!(ba2 < ba1, "ba2={ba2} should be less than ba1={ba1}");
    }
    #[test]
    fn asperity_plasticity_index_computed() {
        let model = AsperityModel::new(1e12, 1e-6, 1e-6, 1e11);
        let pi = model.plasticity_index(2e9);
        assert!(pi.is_finite() && pi > 0.0, "pi={pi}");
    }
    #[test]
    fn ehl_film_positive_thickness() {
        let ehl = EhlFilm::new(0.01, 2e-8, 200e9, 0.01);
        let h = ehl.central_film_thickness_line(1.0, 1000.0);
        assert!(h > 0.0, "h={h}");
    }
    #[test]
    fn ehl_film_increases_with_speed() {
        let ehl = EhlFilm::new(0.01, 2e-8, 200e9, 0.01);
        let h1 = ehl.central_film_thickness_line(1.0, 1000.0);
        let h2 = ehl.central_film_thickness_line(10.0, 1000.0);
        assert!(h2 > h1, "h2={h2} should be larger than h1={h1}");
    }
    #[test]
    fn ehl_lambda_ratio_large_roughness() {
        let ehl = EhlFilm::new(0.01, 2e-8, 200e9, 0.01);
        let lambda = ehl.lambda_ratio(1.0, 1000.0, 1e-3);
        assert!(lambda < 10.0, "lambda={lambda}");
    }
    #[test]
    fn lubrication_regime_boundary() {
        let regime = LubricationRegime::new(0.15, 0.005, 1e-8);
        assert_eq!(regime.classify(0.5), LubricationRegimeKind::Boundary);
    }
    #[test]
    fn lubrication_regime_full_film() {
        let regime = LubricationRegime::new(0.15, 0.005, 1e-8);
        assert_eq!(regime.classify(15.0), LubricationRegimeKind::FullFilm);
    }
    #[test]
    fn lubrication_friction_bounded() {
        let regime = LubricationRegime::new(0.15, 0.005, 1e-8);
        let mu = regime.friction_coefficient(0.01, 1.0, 1e6);
        assert!(
            mu >= regime.mu_fullfilm && mu <= regime.mu_boundary,
            "mu={mu}"
        );
    }
    #[test]
    fn coating_effective_hardness_surface() {
        let coating = CoatingMaterial::new(200e9, 400e9, 1e-6, 20e9, 10.0, 5000.0);
        let h = coating.effective_hardness(0.05e-6);
        assert!((h - 20e9).abs() < 1.0, "h={h}");
    }
    #[test]
    fn coating_delamination_check() {
        let coating = CoatingMaterial::new(200e9, 400e9, 1e-6, 20e9, 10.0, 5000.0);
        assert!(!coating.will_delaminate(5.0));
        assert!(coating.will_delaminate(15.0));
    }
    #[test]
    fn coating_modulus_ratio() {
        let coating = CoatingMaterial::new(200e9, 400e9, 1e-6, 20e9, 10.0, 5000.0);
        assert!((coating.modulus_ratio() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn friction_evolution_starts_at_initial() {
        let fe = FrictionEvolution::new(0.5, 0.1, 100.0);
        assert!((fe.friction_at(0.0) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn friction_evolution_converges_to_steady() {
        let fe = FrictionEvolution::new(0.5, 0.1, 100.0);
        let mu = fe.friction_at(1000.0);
        assert!((mu - 0.1).abs() < 0.01, "mu={mu}");
    }
    #[test]
    fn friction_evolution_step() {
        let mut fe = FrictionEvolution::new(0.5, 0.1, 100.0);
        fe.step(50.0);
        assert!((fe.time - 50.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_wear_flash_temp_above_ambient() {
        let tw = ThermalWear::new(0.3, 50.0, 50.0, 300.0, 80e3, 1e-6);
        let t = tw.flash_temperature(100.0, 1.0, 1e-3);
        assert!(t > 300.0, "t={t}");
    }
    #[test]
    fn thermal_wear_oxidative_rate_positive() {
        let tw = ThermalWear::new(0.3, 50.0, 50.0, 300.0, 80e3, 1e-6);
        let rate = tw.oxidative_wear_rate(600.0);
        assert!(rate > 0.0, "rate={rate}");
    }
    #[test]
    fn thermal_wear_volume_finite() {
        let tw = ThermalWear::new(0.3, 50.0, 50.0, 300.0, 80e3, 1e-6);
        let v = tw.total_wear_volume(100.0, 1.0, 1e-3, 1.0);
        assert!(v.is_finite() && v >= 0.0, "v={v}");
    }
    #[test]
    fn hertz_elliptical_semi_axis_a_positive() {
        let he = HertzElliptical::new(200e9, 0.3, 200e9, 0.3, 0.02, 0.01, 0.02, 0.01);
        let a = he.semi_axis_a(100.0);
        assert!(a > 0.0, "a={a}");
    }
    #[test]
    fn hertz_elliptical_semi_axis_b_le_a() {
        let he = HertzElliptical::new(200e9, 0.3, 200e9, 0.3, 0.02, 0.01, 0.02, 0.01);
        let a = he.semi_axis_a(100.0);
        let b = he.semi_axis_b(100.0);
        assert!(b <= a + 1e-15, "b={b}, a={a}");
    }
    #[test]
    fn hertz_elliptical_max_pressure_positive() {
        let he = HertzElliptical::new(200e9, 0.3, 200e9, 0.3, 0.02, 0.01, 0.02, 0.01);
        let p0 = he.max_pressure(100.0);
        assert!(p0 > 0.0, "p0={p0}");
    }
    #[test]
    fn hertz_elliptical_contact_area_positive() {
        let he = HertzElliptical::new(200e9, 0.3, 200e9, 0.3, 0.02, 0.01, 0.02, 0.01);
        let area = he.contact_area(200.0);
        assert!(area > 0.0, "area={area}");
    }
    #[test]
    fn holm_volume_linear_in_load() {
        let h = HolmWear::new(1e-15);
        let v1 = h.volume_loss(100.0, 1.0);
        let v2 = h.volume_loss(200.0, 1.0);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn holm_wear_depth_finite() {
        let h = HolmWear::new(1e-15);
        let d = h.wear_depth(1000.0, 1000.0, 1e-4);
        assert!(d.is_finite() && d > 0.0, "d={d}");
    }
    #[test]
    fn holm_zero_area_returns_zero() {
        let h = HolmWear::new(1e-15);
        let d = h.wear_depth(100.0, 1.0, 0.0);
        assert!(d.abs() < 1e-30);
    }
    #[test]
    fn abrasive_two_body_gt_three_body() {
        let a = AbrasiveWear::new(1e-4, 1e9, 10e-6, 0.5);
        let v2 = a.two_body_volume(100.0, 1.0);
        let v3 = a.three_body_volume(100.0, 1.0);
        assert!(v2 > v3, "v2={v2}, v3={v3}");
    }
    #[test]
    fn abrasive_specific_wear_rate() {
        let a = AbrasiveWear::new(1e-4, 1e9, 10e-6, 0.5);
        let swr = a.specific_wear_rate();
        assert!((swr - 1e-13).abs() < 1e-20, "swr={swr}");
    }
    #[test]
    fn abrasive_scratch_depth_positive() {
        let a = AbrasiveWear::new(1e-4, 1e9, 10e-6, 0.5);
        let d = a.scratch_depth(1e-3);
        assert!(d.is_finite() && d > 0.0, "d={d}");
    }
    #[test]
    fn roughness_rq_ra_ratio_close_to_1_25() {
        let r = SurfaceRoughness::ground_steel();
        let ratio = r.rq_ra_ratio();
        assert!(ratio > 1.0 && ratio < 2.0, "ratio={ratio}");
    }
    #[test]
    fn roughness_combined_rq_gt_individual() {
        let r1 = SurfaceRoughness::ground_steel();
        let r2 = SurfaceRoughness::lapped();
        let comb = r1.combined_rq(&r2);
        assert!(comb > r1.rq && comb > r2.rq, "comb={comb}");
    }
    #[test]
    fn roughness_bearing_area_at_zero() {
        let r = SurfaceRoughness::ground_steel();
        let ba = r.bearing_area(0.0);
        assert!((0.0..=1.0).contains(&ba), "ba={ba}");
    }
    #[test]
    fn roughness_texture_index_positive() {
        let r = SurfaceRoughness::ground_steel();
        let sti = r.surface_texture_index();
        assert!(sti > 0.0, "sti={sti}");
    }
    #[test]
    fn bearing_l10_mrv_ball_positive() {
        let b = BearingLife::ball_bearing(20_000.0, 5_000.0);
        let l10 = b.l10_mrv();
        assert!(l10 > 0.0, "l10={l10}");
    }
    #[test]
    fn bearing_l10_hours_decreases_with_rpm() {
        let b = BearingLife::ball_bearing(20_000.0, 5_000.0);
        let h1 = b.l10_hours(1000.0);
        let h2 = b.l10_hours(3000.0);
        assert!(h1 > h2, "h1={h1}, h2={h2}");
    }
    #[test]
    fn bearing_reliability_at_zero_is_one() {
        let b = BearingLife::ball_bearing(20_000.0, 5_000.0);
        let r = b.reliability_at(0.0, 1000.0);
        assert!((r - 1.0).abs() < 1e-10, "r={r}");
    }
    #[test]
    fn bearing_l10m_geq_l10_with_good_lubrication() {
        let mut b = BearingLife::ball_bearing(20_000.0, 5_000.0);
        b.a_iso = 2.0;
        b.a_skf = 1.5;
        let l10h = b.l10_hours(1000.0);
        let l10mh = b.l10m_hours(1000.0, 1.0);
        assert!(l10mh >= l10h, "l10mh={l10mh}, l10h={l10h}");
    }
    #[test]
    fn bearing_equivalent_load() {
        let p = BearingLife::equivalent_load(1000.0, 500.0, 0.56, 1.7);
        assert!(p > 1000.0 * 0.56, "p={p}");
    }
    #[test]
    fn gear_pitch_line_velocity_positive() {
        let g = GearTooth::spur_pair(3.0, 20, 40, 20.0, 30.0, 500.0, 1500.0);
        let v = g.pitch_line_velocity();
        assert!(v > 0.0, "v={v}");
    }
    #[test]
    fn gear_contact_ratio_greater_than_1() {
        let g = GearTooth::spur_pair(3.0, 20, 40, 20.0, 30.0, 500.0, 1500.0);
        let cr = g.contact_ratio();
        assert!(cr > 1.0, "cr={cr}");
    }
    #[test]
    fn gear_hertz_stress_positive() {
        let g = GearTooth::spur_pair(3.0, 20, 40, 20.0, 30.0, 500.0, 1500.0);
        let sh = g.hertz_contact_stress();
        assert!(sh > 0.0, "sh={sh}");
    }
    #[test]
    fn gear_ratio_correct() {
        let g = GearTooth::spur_pair(3.0, 20, 40, 20.0, 30.0, 500.0, 1500.0);
        assert!((g.gear_ratio() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn gear_friction_power_loss_positive() {
        let g = GearTooth::spur_pair(3.0, 20, 40, 20.0, 30.0, 500.0, 1500.0);
        let p = g.friction_power_loss(0.05);
        assert!(p >= 0.0, "p={p}");
    }
    #[test]
    fn brake_torque_positive() {
        let b = BrakePad::semi_metallic();
        assert!(b.braking_torque() > 0.0);
    }
    #[test]
    fn brake_temp_rise_increases_with_speed() {
        let b = BrakePad::semi_metallic();
        let t1 = b.temperature_rise_steady(5.0);
        let t2 = b.temperature_rise_steady(20.0);
        assert!(t2 > t1, "t2={t2}, t1={t1}");
    }
    #[test]
    fn brake_wear_depth_positive() {
        let b = BrakePad::semi_metallic();
        let d = b.wear_depth(1_000_000.0);
        assert!(d > 0.0, "d={d}");
    }
    #[test]
    fn brake_remaining_life_positive() {
        let b = BrakePad::semi_metallic();
        let life = b.remaining_life(b.tp, 0.003);
        assert!(life > 0.0, "life={life}");
    }
    #[test]
    fn fretting_max_tangential_force() {
        let f = FrettingFatigue::new(100.0, 0.3, 5e-6, 1e6, 300e6);
        assert!((f.max_tangential_force() - 30.0).abs() < 0.01);
    }
    #[test]
    fn fretting_gross_slip_detected() {
        let f = FrettingFatigue::new(100.0, 0.3, 5e-6, 1e6, 300e6);
        assert!(!f.is_gross_slip(), "should be partial slip");
        let f2 = FrettingFatigue::new(100.0, 0.3, 50e-6, 1e6, 300e6);
        assert!(f2.is_gross_slip(), "should be gross slip");
    }
    #[test]
    fn fretting_energy_per_cycle_positive() {
        let f = FrettingFatigue::new(100.0, 0.3, 50e-6, 1e6, 300e6);
        assert!(f.energy_per_cycle() > 0.0);
    }
    #[test]
    fn fretting_cycles_to_crack_positive() {
        let f = FrettingFatigue::new(100.0, 0.3, 5e-6, 1e6, 300e6);
        let n = f.cycles_to_crack(100e6, 1e-3);
        assert!(n > 0.0 && n.is_finite(), "n={n}");
    }
    #[test]
    fn hamrock_central_film_positive() {
        let hd = HamrockDowson::new(0.01, 2e-8, 200e9, 0.01, 0.01);
        let h = hd.central_film_thickness(1.0, 500.0);
        assert!(h > 0.0, "h={h}");
    }
    #[test]
    fn hamrock_min_film_positive() {
        let hd = HamrockDowson::new(0.01, 2e-8, 200e9, 0.01, 0.01);
        let h = hd.minimum_film_thickness(1.0, 500.0);
        assert!(h > 0.0, "h={h}");
    }
    #[test]
    fn hamrock_central_geq_min_film() {
        let hd = HamrockDowson::new(0.01, 2e-8, 200e9, 0.01, 0.01);
        let hc = hd.central_film_thickness(1.0, 500.0);
        let hm = hd.minimum_film_thickness(1.0, 500.0);
        assert!(hc >= hm * 0.5, "hc={hc}, hm={hm}");
    }
    #[test]
    fn hamrock_lambda_increases_with_speed() {
        let hd = HamrockDowson::new(0.01, 2e-8, 200e9, 0.01, 0.01);
        let l1 = hd.lambda_ratio(1.0, 500.0, 1e-7);
        let l2 = hd.lambda_ratio(5.0, 500.0, 1e-7);
        assert!(l2 > l1, "l2={l2}, l1={l1}");
    }
    #[test]
    fn stribeck_mu_bounded() {
        let s = StribeckCurve::new(0.15, 0.005, 1e-8);
        let mu = s.mu_at_hersey(1e-10);
        assert!((0.0..=0.15 + 0.001).contains(&mu), "mu={mu}");
    }
    #[test]
    fn stribeck_classify_boundary_low_lambda() {
        let s = StribeckCurve::new(0.15, 0.005, 1e-8);
        assert_eq!(s.classify(0.5), LubricationRegimeKind::Boundary);
    }
    #[test]
    fn stribeck_viscous_shear_zero_film() {
        let s = StribeckCurve::new(0.15, 0.005, 1e-8);
        assert_eq!(s.viscous_shear_stress(0.01, 1.0, 0.0), 0.0);
    }
    #[test]
    fn bowden_tabor_mu_adhesion() {
        let b = BowdenTaborFriction::new(2e8, 1e9, 0.05);
        assert!((b.mu_adhesion() - 0.2).abs() < 1e-10);
    }
    #[test]
    fn bowden_tabor_total_mu_gt_adhesion() {
        let b = BowdenTaborFriction::new(2e8, 1e9, 0.05);
        assert!(b.mu_total() > b.mu_adhesion());
    }
    #[test]
    fn bowden_tabor_friction_force_positive() {
        let b = BowdenTaborFriction::new(2e8, 1e9, 0.05);
        let ff = b.friction_force(100.0);
        assert!(ff > 0.0, "ff={ff}");
    }
    #[test]
    fn lubricant_viscosity_decreases_with_temperature() {
        let v = LubricantViscosity::sae10w40();
        let eta40 = v.eta_reynolds(313.15);
        let eta100 = v.eta_reynolds(373.15);
        assert!(eta40 > eta100, "eta40={eta40}, eta100={eta100}");
    }
    #[test]
    fn lubricant_barus_increases_with_pressure() {
        let v = LubricantViscosity::sae10w40();
        let eta0 = v.eta_barus(313.15, 2e-8, 0.0);
        let eta_high = v.eta_barus(313.15, 2e-8, 1e9);
        assert!(eta_high > eta0, "eta_high={eta_high}, eta0={eta0}");
    }
    #[test]
    fn lubricant_vi_positive() {
        let v = LubricantViscosity::sae10w40();
        let vi = v.viscosity_index();
        assert!(vi > 0.0, "vi={vi}");
    }
    #[test]
    fn cylinder_half_width_positive() {
        let c = CylinderContact::new(200e9, 0.3, 0.02, 200e9, 0.3, 0.02);
        let a = c.half_width(1000.0);
        assert!(a > 0.0, "a={a}");
    }
    #[test]
    fn cylinder_max_pressure_positive() {
        let c = CylinderContact::new(200e9, 0.3, 0.02, 200e9, 0.3, 0.02);
        assert!(c.max_pressure(1000.0) > 0.0);
    }
    #[test]
    fn cylinder_pressure_zero_outside_contact() {
        let c = CylinderContact::new(200e9, 0.3, 0.02, 200e9, 0.3, 0.02);
        let a = c.half_width(1000.0);
        let p = c.pressure_at(1000.0, a * 2.0);
        assert!((p).abs() < 1e-10, "p={p}");
    }
    #[test]
    fn scuffing_safety_factor_above_one_low_heat() {
        let s = ScuffingModel::steel_steel();
        let sf = s.safety_factor(1e4, 1e-3, 0.1);
        assert!(sf > 1.0, "sf={sf}");
    }
    #[test]
    fn scuffing_contact_temp_above_bulk() {
        let s = ScuffingModel::steel_steel();
        let tc = s.contact_temperature(1e6, 1e-3, 1.0);
        assert!(tc > s.t_bulk, "tc={tc}");
    }
    #[test]
    fn tribological_state_initial_lambda_zero() {
        let state = TribologicalState::new(1.0, 100.0);
        assert_eq!(state.lambda, 0.0);
    }
    #[test]
    fn tribological_state_update_film() {
        let mut state = TribologicalState::new(1.0, 100.0);
        state.update_film(1e-6, 0.5e-6);
        assert!((state.lambda - 2.0).abs() < 1e-10);
        assert_eq!(state.regime, LubricationRegimeKind::Mixed);
    }
    #[test]
    fn tribological_state_accumulate_wear() {
        let mut state = TribologicalState::new(1.0, 100.0);
        state.accumulate_wear(1e-4, 1e9, 1.0);
        assert!(state.wear_volume > 0.0, "wear_volume={}", state.wear_volume);
    }
    #[test]
    fn tribological_state_power_dissipation() {
        let state = TribologicalState::new(2.0, 100.0);
        assert!((state.power_dissipation() - 0.15 * 100.0 * 2.0).abs() < 1e-10);
    }
    #[test]
    fn coulomb_friction_force_basic() {
        let f = coulomb_friction_force(100.0, 0.3);
        assert!((f - 30.0).abs() < 1e-10, "f={f}");
    }
    #[test]
    fn coulomb_friction_force_zero_normal() {
        assert_eq!(coulomb_friction_force(0.0, 0.5), 0.0);
    }
    #[test]
    fn stribeck_friction_at_zero_velocity_gives_static() {
        let f = stribeck_friction(0.0, 0.4, 0.2, 0.1);
        assert!((f - 0.4).abs() < 1e-6, "f={f}");
    }
    #[test]
    fn stribeck_friction_high_velocity_approaches_kinetic() {
        let f = stribeck_friction(100.0, 0.4, 0.2, 0.1);
        assert!((f - 0.2).abs() < 0.01, "f={f}");
    }
    #[test]
    fn viscous_friction_basic() {
        let f = viscous_friction(1.0, 0.1, 1.0, 0.001);
        assert!((f - 100.0).abs() < 1e-6, "f={f}");
    }
    #[test]
    fn viscous_friction_zero_velocity() {
        assert_eq!(viscous_friction(0.0, 0.1, 1.0, 0.001), 0.0);
    }
    #[test]
    fn hertz_contact_params_reduced_modulus() {
        let h = HertzContactParams::new(0.01, f64::INFINITY, 200e9, 0.3, 200e9, 0.3);
        let e_star = h.reduced_modulus();
        assert!(e_star > 0.0, "e_star={e_star}");
    }
    #[test]
    fn hertz_contact_params_contact_radius_positive() {
        let h = HertzContactParams::new(0.01, f64::INFINITY, 200e9, 0.3, 200e9, 0.3);
        let a = h.contact_radius(100.0);
        assert!(a > 0.0, "a={a}");
    }
    #[test]
    fn hertz_contact_params_max_pressure_positive() {
        let h = HertzContactParams::new(0.01, f64::INFINITY, 200e9, 0.3, 200e9, 0.3);
        let p = h.max_pressure(100.0);
        assert!(p > 0.0, "p={p}");
    }
    #[test]
    fn hertz_contact_params_stiffness_positive() {
        let h = HertzContactParams::new(0.01, f64::INFINITY, 200e9, 0.3, 200e9, 0.3);
        let k = h.contact_stiffness(100.0);
        assert!(k > 0.0, "k={k}");
    }
    #[test]
    fn hertz_contact_params_max_shear_stress_positive() {
        let h = HertzContactParams::new(0.01, f64::INFINITY, 200e9, 0.3, 200e9, 0.3);
        let tau = h.max_shear_stress(100.0);
        assert!(tau > 0.0, "tau={tau}");
    }
    #[test]
    fn wear_model_volume_proportional_to_distance() {
        let wm = WearModel::new(1e-4, 1e9);
        let v1 = wm.wear_volume(100.0, 1.0);
        let v2 = wm.wear_volume(100.0, 2.0);
        assert!((v2 / v1 - 2.0).abs() < 1e-10, "v1={v1}, v2={v2}");
    }
    #[test]
    fn wear_model_rate_proportional_to_force() {
        let wm = WearModel::new(1e-4, 1e9);
        let r1 = wm.wear_rate(100.0, 1.0);
        let r2 = wm.wear_rate(200.0, 1.0);
        assert!((r2 / r1 - 2.0).abs() < 1e-10, "r1={r1}, r2={r2}");
    }
    #[test]
    fn adhesive_wear_rate_positive() {
        let r = adhesive_wear_rate(1e-6, 1e9, 0.01);
        assert!(r > 0.0, "r={r}");
    }
    #[test]
    fn abrasive_wear_factor_positive() {
        let f = abrasive_wear_factor(1e-5, 0.5);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn ehl_film_central_thickness_positive() {
        let ehl = ElastohydrodynamicFilm::new(0.04, 2e-8, 1.0, 0.01, 110e9);
        let h = ehl.central_film_thickness(1000.0);
        assert!(h > 0.0, "h={h}");
    }
    #[test]
    fn ehl_film_minimum_less_than_central() {
        let ehl = ElastohydrodynamicFilm::new(0.04, 2e-8, 1.0, 0.01, 110e9);
        let hc = ehl.central_film_thickness(1000.0);
        let hm = ehl.minimum_film_thickness(1000.0);
        assert!(hm < hc, "hm={hm}, hc={hc}");
    }
    #[test]
    fn ehl_lambda_ratio_positive() {
        let ehl = ElastohydrodynamicFilm::new(0.04, 2e-8, 1.0, 0.01, 110e9);
        let lambda = ehl.lambda_ratio(1e-7, 1000.0);
        assert!(lambda > 0.0, "lambda={lambda}");
    }
    #[test]
    fn asperity_contact_fraction_boundary_values() {
        assert!((asperity_contact_fraction(0.0) - 1.0).abs() < 1e-10);
        assert!((asperity_contact_fraction(f64::INFINITY)).abs() < 1e-10);
    }
    #[test]
    fn asperity_contact_fraction_monotone() {
        let f1 = asperity_contact_fraction(1.0);
        let f2 = asperity_contact_fraction(3.0);
        assert!(f2 < f1, "f1={f1}, f2={f2}");
    }
    #[test]
    fn friction_oscillator_initial_sticking() {
        let fo = FrictionOscillator::new(1.0, 100.0, 0.4, 0.3);
        assert!(fo.is_sticking(0.0, 0.0));
    }
    #[test]
    fn friction_oscillator_step_returns_finite() {
        let mut fo = FrictionOscillator::new(1.0, 100.0, 0.4, 0.3);
        let (x, v) = fo.step(0.0, 0.0, 0.1, 1e-4);
        assert!(x.is_finite() && v.is_finite(), "x={x}, v={v}");
    }
    #[test]
    fn flash_temperature_standalone_positive() {
        let t = flash_temperature(1e6, 50.0, 1.0, 1e-4);
        assert!(t > 0.0, "t={t}");
    }
    #[test]
    fn flash_temperature_scales_with_flux() {
        let t1 = flash_temperature(1e6, 50.0, 1.0, 1e-4);
        let t2 = flash_temperature(2e6, 50.0, 1.0, 1e-4);
        assert!(t2 > t1, "t1={t1}, t2={t2}");
    }
}
