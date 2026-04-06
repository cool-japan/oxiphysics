//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

// Types used by functions below

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-300 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / len)
    }
}
/// Added mass of a sphere: m_a = 0.5 * rho * (4/3 * pi * r^3).
pub fn added_mass_sphere(rho: f64, r: f64) -> f64 {
    let volume = (4.0 / 3.0) * PI * r * r * r;
    0.5 * rho * volume
}
/// Added mass of a cylinder: m_a = rho * pi * r^2 * L.
pub fn added_mass_cylinder(rho: f64, r: f64, l: f64) -> f64 {
    rho * PI * r * r * l
}
/// Drag coefficient for a sphere using the Schiller-Naumann correlation.
///
/// Cd = 24/Re * (1 + 0.15 * Re^0.687)  for Re <= 1000, else 0.44.
pub fn drag_coefficient_sphere(re: f64) -> f64 {
    if re < 1e-10 {
        return f64::INFINITY;
    }
    if re <= 1000.0 {
        (24.0 / re) * (1.0 + 0.15 * re.powf(0.687))
    } else {
        0.44
    }
}
/// Drag coefficient for a cylinder in cross-flow (empirical fit).
pub fn drag_coefficient_cylinder(re: f64) -> f64 {
    if re < 1e-10 {
        return f64::INFINITY;
    }
    if re < 1.0 {
        8.0 * PI / (re * (2.002 - (re / 4.0).max(0.001).ln()))
    } else {
        1.0 + 10.0 / re.powf(2.0 / 3.0)
    }
}
/// Slamming force: F = cs * rho * v^2 * A.
///
/// `cs` is the slamming coefficient (dimensionless).
pub fn slamming_force(rho: f64, v: f64, a: f64, cs: f64) -> f64 {
    cs * rho * v * v * a
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AddedMass;

    use crate::BuoyancyForce;
    use crate::FloatingBody;

    use crate::FluidProperties;

    use crate::HydrodynamicDrag;
    use crate::MorisonElement;

    use crate::SphereBuoyancy;
    use crate::SubmergedBodyCoupling;

    use crate::VivParams;
    use crate::VortexInducedVibration;
    use crate::WaveForce;
    use crate::fluid_coupling::functions::PI;
    pub(super) const WATER: FluidProperties = FluidProperties {
        density: 1000.0,
        dynamic_viscosity: 1e-3,
    };
    #[test]
    fn test_archimedes() {
        let g = [0.0, 0.0, -9.81];
        let f = BuoyancyForce::archimedes(1000.0, 1.0, g);
        assert!((f[2] - 9810.0).abs() < 1e-9);
        assert_eq!(f[0], 0.0);
        assert_eq!(f[1], 0.0);
    }
    #[test]
    fn test_stokes_drag() {
        let vel = [1.0, 0.0, 0.0];
        let f = HydrodynamicDrag::stokes_drag(0.01, vel, &WATER);
        assert!(f[0] < 0.0);
        let expected = -6.0 * PI * 1e-3 * 0.01;
        assert!((f[0] - expected).abs() < 1e-12);
    }
    #[test]
    fn test_drag_coefficient_sphere_re1() {
        let cd = HydrodynamicDrag::drag_coefficient_sphere(1.0);
        let expected = 24.0 + 6.0 / 2.0 + 0.4;
        assert!((cd - expected).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_added_mass() {
        let r = 1.0;
        let rho = 1000.0;
        let m_a = AddedMass::sphere_added_mass(r, rho);
        let expected = 0.5 * rho * (4.0 / 3.0) * PI * r * r * r;
        assert!((m_a - expected).abs() < 1e-9);
    }
    #[test]
    fn test_floating_body_stable() {
        let body = FloatingBody {
            mass: 100.0,
            volume: 0.1,
            center_of_mass: [0.0, 0.0, 0.0],
            center_of_buoyancy: [0.0, 0.0, -0.05],
            metacentric_height: 0.5,
            position: [0.0; 3],
            orientation_angle: 0.0,
            linear_vel: [0.0; 3],
            angular_vel: 0.0,
        };
        assert!(body.is_stable());
        let unstable = FloatingBody {
            metacentric_height: -0.1,
            ..body
        };
        assert!(!unstable.is_stable());
    }
    #[test]
    fn test_net_force_neutrally_buoyant() {
        let coupling = SubmergedBodyCoupling {
            body_density: 1000.0,
            body_volume: 1.0,
            fluid: FluidProperties {
                density: 1000.0,
                dynamic_viscosity: 1e-3,
            },
            cd: 0.47,
        };
        let g = [0.0, 0.0, -9.81];
        let f = coupling.net_force([0.0; 3], [0.0; 3], g);
        assert!(length(f) < 1e-9, "net force should be zero, got {:?}", f);
    }
    #[test]
    fn test_buoyancy_torque_zero_when_cob_eq_com() {
        let com = [1.0, 2.0, 3.0];
        let cob = [1.0, 2.0, 3.0];
        let force = [0.0, 0.0, 9810.0];
        let tau = BuoyancyForce::buoyancy_torque(cob, com, force);
        assert!(length(tau) < 1e-15);
    }
    #[test]
    fn test_sphere_partial_buoyancy_fully_submerged() {
        let r = 0.5;
        let g = [0.0, 0.0, -9.81];
        let f = BuoyancyForce::sphere_partial_buoyancy(r, 2.0 * r, 1000.0, g);
        let expected_vol = (4.0 / 3.0) * PI * r * r * r;
        let expected_f = 1000.0 * expected_vol * 9.81;
        assert!(
            (f[2] - expected_f).abs() < 0.01,
            "fully submerged buoyancy: got {}, expected {}",
            f[2],
            expected_f
        );
    }
    #[test]
    fn test_sphere_partial_buoyancy_half_submerged() {
        let r = 0.5;
        let g = [0.0, 0.0, -9.81];
        let f_full = BuoyancyForce::sphere_partial_buoyancy(r, 2.0 * r, 1000.0, g);
        let f_half = BuoyancyForce::sphere_partial_buoyancy(r, r, 1000.0, g);
        assert!(
            (f_half[2] - f_full[2] * 0.5).abs() < f_full[2] * 0.01,
            "half submerged: got {}, expected ~{}",
            f_half[2],
            f_full[2] * 0.5
        );
    }
    #[test]
    fn test_oseen_drag_greater_than_stokes() {
        let vel = [0.1, 0.0, 0.0];
        let f_stokes = HydrodynamicDrag::stokes_drag(0.01, vel, &WATER);
        let f_oseen = HydrodynamicDrag::oseen_drag(0.01, vel, &WATER);
        assert!(
            f_oseen[0].abs() >= f_stokes[0].abs(),
            "Oseen should give >= Stokes drag"
        );
    }
    #[test]
    fn test_drag_coefficient_cylinder() {
        let cd = HydrodynamicDrag::drag_coefficient_cylinder(1000.0);
        assert!(cd > 0.9 && cd < 2.5, "cylinder Cd at Re=1000: got {cd}");
    }
    #[test]
    fn test_drag_coefficient_flat_plate() {
        let cd = HydrodynamicDrag::drag_coefficient_flat_plate();
        assert!((cd - 1.98).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_drag_opposes_velocity() {
        let vel = [2.0, 0.0, 0.0];
        let f = HydrodynamicDrag::sphere_drag(0.1, vel, &WATER);
        assert!(f[0] < 0.0, "drag should oppose velocity");
    }
    #[test]
    fn test_cylinder_drag_opposes_velocity() {
        let vel = [2.0, 0.0, 0.0];
        let f = HydrodynamicDrag::cylinder_drag(0.05, 1.0, vel, &WATER);
        assert!(f[0] < 0.0, "drag should oppose velocity");
    }
    #[test]
    fn test_arbitrary_drag() {
        let vel = [1.0, 0.0, 0.0];
        let f = HydrodynamicDrag::arbitrary_drag(vel, 0.1, 1.0, &WATER);
        assert!((f[0] - (-50.0)).abs() < 1e-6);
    }
    #[test]
    fn test_cylinder_added_mass() {
        let m_a = AddedMass::cylinder_added_mass(0.5, 2.0, 1000.0);
        let expected = 1000.0 * PI * 0.25 * 2.0;
        assert!((m_a - expected).abs() < 1e-9);
    }
    #[test]
    fn test_ellipsoid_added_mass_sphere_limit() {
        let r = 0.5;
        let m_sphere = AddedMass::sphere_added_mass(r, 1000.0);
        let m_ellipsoid = AddedMass::ellipsoid_added_mass_a(r, r, 1000.0);
        assert!(
            (m_sphere - m_ellipsoid).abs() < m_sphere * 0.01,
            "sphere limit: sphere={m_sphere}, ellipsoid={m_ellipsoid}"
        );
    }
    #[test]
    fn test_added_mass_force() {
        let f = AddedMass::added_mass_force(100.0, [1.0, 0.0, 0.0]);
        assert!((f[0] - (-100.0)).abs() < 1e-10);
    }
    #[test]
    fn test_floating_body_roll_period() {
        let body = FloatingBody {
            mass: 1000.0,
            volume: 1.0,
            center_of_mass: [0.0, 0.0, 0.0],
            center_of_buoyancy: [0.0, 0.0, -0.1],
            metacentric_height: 0.5,
            position: [0.0; 3],
            orientation_angle: 0.0,
            linear_vel: [0.0; 3],
            angular_vel: 0.0,
        };
        let period = body.roll_period(500.0);
        let expected = 2.0 * PI * (500.0_f64 / (1000.0 * 9.81 * 0.5)).sqrt();
        assert!(
            (period - expected).abs() < 1e-6,
            "period={period}, expected={expected}"
        );
    }
    #[test]
    fn test_floating_body_roll_dynamics() {
        let mut body = FloatingBody {
            mass: 1000.0,
            volume: 1.0,
            center_of_mass: [0.0, 0.0, 0.0],
            center_of_buoyancy: [0.0, 0.0, -0.1],
            metacentric_height: 0.5,
            position: [0.0; 3],
            orientation_angle: 0.1,
            linear_vel: [0.0; 3],
            angular_vel: 0.0,
        };
        for _ in 0..10000 {
            body.step_roll(500.0, 100.0, 0.001);
        }
        assert!(
            body.orientation_angle.abs() < 0.1,
            "should decay toward zero: got {}",
            body.orientation_angle
        );
    }
    #[test]
    fn test_wave_particle_acceleration() {
        let acc = WaveForce::wave_particle_acceleration(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        assert!(acc[0].abs() < 1e-10);
        assert!((acc[2] - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_wave_velocity_with_depth_surface() {
        let v_simple = WaveForce::wave_particle_velocity(1.0, 1.0, 0.5, 0.0, 0.0, 0.0);
        let v_depth = WaveForce::wave_velocity_with_depth(1.0, 1.0, 0.5, 0.0, 0.0, 0.0, 100.0);
        assert!(
            (v_simple[0] - v_depth[0]).abs() < 0.1,
            "surface velocity should be similar: simple={}, depth={}",
            v_simple[0],
            v_depth[0]
        );
    }
    #[test]
    fn test_deep_water_dispersion() {
        let omega = WaveForce::deep_water_omega(0.1, 9.81);
        assert!((omega - (9.81 * 0.1_f64).sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_deep_water_wavelength() {
        let l = WaveForce::deep_water_wavelength(10.0, 9.81);
        let expected = 9.81 * 100.0 / (2.0 * PI);
        assert!((l - expected).abs() < 1e-6);
    }
    #[test]
    fn test_morison_max_components() {
        let (fd, fi) = WaveForce::morison_max_components(0.5, 10.0, 1025.0, 2.0, 0.6, 1.0, 2.0);
        assert!(fd > 0.0);
        assert!(fi > 0.0);
        let u_max = 2.0 * 0.6;
        let a_max = 2.0 * 0.6 * 0.6;
        let area = 0.25 * PI * 0.25;
        let vol = area * 10.0;
        let fd_expected = 0.5 * 1025.0 * 1.0 * area * u_max * u_max;
        let fi_expected = 1025.0 * vol * 2.0 * a_max;
        assert!((fd - fd_expected).abs() < 1e-6);
        assert!((fi - fi_expected).abs() < 1e-6);
    }
    #[test]
    fn test_viv_shedding_frequency() {
        let viv = VortexInducedVibration::new(0.5, 1.0, 0.01, 100.0);
        let fs = viv.shedding_frequency(5.0);
        assert!((fs - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_viv_reduced_velocity() {
        let viv = VortexInducedVibration::new(0.5, 1.0, 0.01, 100.0);
        let vr = viv.reduced_velocity(2.5);
        assert!((vr - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_viv_lock_in() {
        let viv = VortexInducedVibration::new(0.5, 1.0, 0.01, 100.0);
        assert!(viv.is_lock_in(2.5), "Vr=5 should be lock-in");
        assert!(!viv.is_lock_in(10.0), "Vr=20 should not be lock-in");
        assert!(!viv.is_lock_in(0.5), "Vr=1 should not be lock-in");
    }
    #[test]
    fn test_viv_peak_amplitude() {
        let viv = VortexInducedVibration::new(0.5, 1.0, 0.01, 100.0);
        let a_d = viv.peak_amplitude_ratio(1000.0);
        assert!(a_d > 0.0 && a_d < 2.0, "A/D ratio: {a_d}");
    }
    #[test]
    fn test_viv_lift_force() {
        let viv = VortexInducedVibration::new(0.5, 1.0, 0.01, 100.0);
        let fl = viv.lift_force_per_length(5.0, 1000.0, 0.5, 0.0);
        assert!(fl.abs() < 1e-10, "lift force at t=0 should be zero");
        let fl_max = viv.lift_force_per_length(5.0, 1000.0, 0.5, 0.125);
        assert!(
            fl_max.abs() > 0.0,
            "lift force should be nonzero at quarter cycle"
        );
    }
    #[test]
    fn test_will_float() {
        let floating = SubmergedBodyCoupling {
            body_density: 500.0,
            body_volume: 1.0,
            fluid: FluidProperties {
                density: 1000.0,
                dynamic_viscosity: 1e-3,
            },
            cd: 0.47,
        };
        assert!(floating.will_float());
        let sinking = SubmergedBodyCoupling {
            body_density: 2000.0,
            ..floating
        };
        assert!(!sinking.will_float());
    }
    #[test]
    fn test_equilibrium_submerged_fraction() {
        let coupling = SubmergedBodyCoupling {
            body_density: 500.0,
            body_volume: 1.0,
            fluid: FluidProperties {
                density: 1000.0,
                dynamic_viscosity: 1e-3,
            },
            cd: 0.47,
        };
        let frac = coupling.equilibrium_submerged_fraction();
        assert!(
            (frac - 0.5).abs() < 1e-10,
            "should float at 50%: got {frac}"
        );
    }
    #[test]
    fn test_net_force_with_added_mass() {
        let coupling = SubmergedBodyCoupling {
            body_density: 2000.0,
            body_volume: 1.0,
            fluid: FluidProperties {
                density: 1000.0,
                dynamic_viscosity: 1e-3,
            },
            cd: 0.47,
        };
        let g = [0.0, 0.0, -9.81];
        let f_no_am = coupling.net_force([0.0; 3], [0.0; 3], g);
        let f_with_am = coupling.net_force_with_added_mass([0.0; 3], [0.0; 3], [0.0, 0.0, -1.0], g);
        assert!(
            f_with_am[2] > f_no_am[2],
            "added mass opposing downward accel should give upward force contribution"
        );
    }
    #[test]
    fn test_sphere_buoyancy_scales_with_fraction() {
        let full = SphereBuoyancy::buoyancy_force(1000.0, 9.81, 1.0, 0.5);
        let half = SphereBuoyancy::buoyancy_force(1000.0, 9.81, 0.5, 0.5);
        assert!(
            (full - 2.0 * half).abs() < 1e-9,
            "buoyancy should scale linearly with fraction"
        );
        assert!(full > 0.0);
    }
    #[test]
    fn test_morison_element_inertia_at_zero_velocity() {
        let elem = MorisonElement {
            diameter: 1.0,
            length: 5.0,
            cm: 2.0,
            cd: 1.0,
        };
        let f = elem.morison_force([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1025.0);
        let vol = 0.25 * PI * 5.0;
        let expected = 1025.0 * vol * 2.0;
        assert!(
            (f[0] - expected).abs() < 1e-6,
            "inertia force: got {}, expected {}",
            f[0],
            expected
        );
        assert!(f[1].abs() < 1e-10);
        assert!(f[2].abs() < 1e-10);
    }
    #[test]
    fn test_vortex_shedding_frequency() {
        let viv = VivParams {
            st: 0.2,
            rho: 1025.0,
            d: 0.5,
        };
        let f_vs = viv.vortex_shedding_frequency(5.0);
        assert!((f_vs - 2.0).abs() < 1e-10, "f_vs: {}", f_vs);
    }
    #[test]
    fn test_viv_lock_in_range() {
        let viv = VivParams {
            st: 0.2,
            rho: 1025.0,
            d: 0.5,
        };
        assert!(viv.lock_in_range(2.0, 5.0), "should be in lock-in range");
        assert!(
            !viv.lock_in_range(2.0, 0.1),
            "should not be in lock-in range at low U"
        );
    }
    #[test]
    fn test_added_mass_sphere_formula() {
        let r = 1.0;
        let rho = 1000.0;
        let m_a = added_mass_sphere(rho, r);
        let expected = 0.5 * rho * (4.0 / 3.0) * PI;
        assert!((m_a - expected).abs() < 1e-9);
    }
    #[test]
    fn test_added_mass_cylinder_formula() {
        let r = 0.5;
        let l = 3.0;
        let rho = 1000.0;
        let m_a = added_mass_cylinder(rho, r, l);
        let expected = rho * PI * r * r * l;
        assert!((m_a - expected).abs() < 1e-9);
    }
    #[test]
    fn test_drag_coefficient_sphere_schiller_naumann() {
        let cd = drag_coefficient_sphere(1.0);
        assert!((cd - 24.0 * 1.15).abs() < 1e-9, "Cd at Re=1: {cd}");
        let cd_high = drag_coefficient_sphere(10000.0);
        assert!((cd_high - 0.44).abs() < 1e-10);
    }
    #[test]
    fn test_drag_coefficient_cylinder_moderate_re() {
        let cd = drag_coefficient_cylinder(1000.0);
        assert!(cd > 0.5 && cd < 3.0, "cylinder Cd at Re=1000: {cd}");
    }
    #[test]
    fn test_slamming_force() {
        let f = slamming_force(1025.0, 3.0, 2.0, 5.15);
        let expected = 5.15 * 1025.0 * 9.0 * 2.0;
        assert!((f - expected).abs() < 1e-6);
    }
    #[test]
    fn test_fluid_properties_water() {
        let w = FluidProperties::water();
        assert!(w.density > 990.0 && w.density < 1010.0);
    }
    #[test]
    fn test_fluid_properties_air() {
        let a = FluidProperties::air();
        assert!(a.density > 1.0 && a.density < 1.5);
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::AddedMassTensor;

    use crate::FlowRegime;

    use crate::FroudeKrylovForce;
    use crate::FroudeNumber;

    use crate::ReynoldsUtils;
    use crate::SloshingModel;

    use crate::SubmergedVolumeFraction;

    use crate::fluid_coupling::functions::PI;
    #[test]
    fn test_added_mass_tensor_sphere_isotropic() {
        let r = 0.5;
        let rho = 1000.0;
        let tensor = AddedMassTensor::sphere(r, rho);
        assert!((tensor.diagonal[0] - tensor.diagonal[1]).abs() < 1e-12);
        assert!((tensor.diagonal[1] - tensor.diagonal[2]).abs() < 1e-12);
        let expected = 0.5 * 1000.0 * (4.0 / 3.0) * PI * 0.125;
        assert!((tensor.diagonal[0] - expected).abs() < 1e-9);
    }
    #[test]
    fn test_added_mass_tensor_cylinder_transverse() {
        let tensor = AddedMassTensor::cylinder(0.5, 2.0, 1000.0);
        assert!((tensor.diagonal[0] - tensor.diagonal[1]).abs() < 1e-12);
        assert!(tensor.diagonal[0] > 0.0);
        assert_eq!(tensor.diagonal[2], 0.0);
        let expected = 1000.0 * PI * 0.25 * 2.0;
        assert!((tensor.diagonal[0] - expected).abs() < 1e-9);
    }
    #[test]
    fn test_added_mass_tensor_flat_plate_normal() {
        let tensor = AddedMassTensor::flat_plate(1.0, 2.0, 1000.0);
        assert_eq!(tensor.diagonal[0], 0.0);
        assert_eq!(tensor.diagonal[1], 0.0);
        assert!(tensor.diagonal[2] > 0.0);
    }
    #[test]
    fn test_added_mass_tensor_force_opposes_accel() {
        let tensor = AddedMassTensor::sphere(0.5, 1000.0);
        let accel = [1.0, 0.0, 0.0];
        let f = tensor.force(accel);
        assert!(f[0] < 0.0, "force should be negative: {}", f[0]);
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn test_added_mass_tensor_effective_mass() {
        let tensor = AddedMassTensor::sphere(0.5, 1000.0);
        let m_eff = tensor.effective_mass(10.0, 0);
        assert!(m_eff > 10.0, "effective mass must exceed body mass");
    }
    #[test]
    fn test_froude_krylov_force_basic() {
        let fk = FroudeKrylovForce::new(1.0, 1.0);
        let a_fluid = [2.0, 0.0, 0.0];
        let f = fk.froude_krylov(1000.0, a_fluid);
        assert!((f[0] - 2000.0).abs() < 1e-9);
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn test_froude_krylov_total_excitation() {
        let fk = FroudeKrylovForce::new(1.0, 1.0);
        let a_fluid = [1.0, 0.0, 0.0];
        let f_total = fk.total_excitation(1000.0, a_fluid);
        let f_fk = fk.froude_krylov(1000.0, a_fluid);
        assert!((f_total[0] - 2.0 * f_fk[0]).abs() < 1e-9);
    }
    #[test]
    fn test_froude_krylov_diffraction() {
        let fk = FroudeKrylovForce::new(0.5, 1.0);
        let a_fluid = [3.0, 0.0, 0.0];
        let f_diff = fk.diffraction(1000.0, a_fluid);
        assert!((f_diff[0] - 1500.0).abs() < 1e-9);
    }
    #[test]
    fn test_radiation_force_opposes_velocity() {
        let f = FroudeKrylovForce::radiation_force(50.0, 10.0, [1.0, 0.0, 0.0], [0.0; 3]);
        assert!((f[0] - (-10.0)).abs() < 1e-9);
    }
    #[test]
    fn test_radiation_force_opposes_acceleration() {
        let f = FroudeKrylovForce::radiation_force(100.0, 0.0, [0.0; 3], [2.0, 0.0, 0.0]);
        assert!((f[0] - (-200.0)).abs() < 1e-9);
    }
    #[test]
    fn test_froude_number_subcritical() {
        let fr = FroudeNumber::froude(1.0, 10.0, 9.81);
        assert!(fr < 1.0);
        assert_eq!(FroudeNumber::classify(fr), FlowRegime::SubCritical);
    }
    #[test]
    fn test_froude_number_supercritical() {
        let fr = FroudeNumber::froude(20.0, 1.0, 9.81);
        assert!(fr > 1.0);
        assert_eq!(FroudeNumber::classify(fr), FlowRegime::SuperCritical);
    }
    #[test]
    fn test_froude_number_critical() {
        assert_eq!(FroudeNumber::classify(1.0), FlowRegime::Critical);
    }
    #[test]
    fn test_froude_conjugate_depth_ratio() {
        let ratio = FroudeNumber::conjugate_depth_ratio(2.0);
        let expected = 0.5 * ((1.0 + 32.0_f64).sqrt() - 1.0);
        assert!((ratio - expected).abs() < 1e-10);
        assert!(ratio > 1.0, "downstream depth should be greater");
    }
    #[test]
    fn test_froude_critical_celerity() {
        let c = FroudeNumber::critical_celerity(1.0, 9.81);
        assert!((c - 9.81_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_froude_hydraulic_jump_detection() {
        assert!(FroudeNumber::has_hydraulic_jump(0.5, 2.0));
        assert!(FroudeNumber::has_hydraulic_jump(2.0, 0.5));
        assert!(!FroudeNumber::has_hydraulic_jump(0.5, 0.8));
        assert!(!FroudeNumber::has_hydraulic_jump(1.5, 2.5));
    }
    #[test]
    fn test_sloshing_natural_period_positive() {
        let model = SloshingModel::new(1.0, 0.5, 100.0, 9.81);
        let t = model.natural_period();
        assert!(t > 0.0 && t.is_finite());
    }
    #[test]
    fn test_sloshing_effective_mass_bounded() {
        let model = SloshingModel::new(1.0, 0.5, 100.0, 9.81);
        let m_eff = model.effective_mass();
        assert!((0.0..=100.0).contains(&m_eff));
    }
    #[test]
    fn test_sloshing_step_changes_angle() {
        let mut model = SloshingModel::new(1.0, 0.5, 100.0, 9.81);
        model.angle = 0.1;
        let angle_before = model.angle;
        model.step(0.0, 0.01);
        assert!(
            (model.angle - angle_before).abs() > 0.0,
            "angle should change during free oscillation"
        );
    }
    #[test]
    fn test_sloshing_force_zero_at_rest() {
        let model = SloshingModel::new(1.0, 0.5, 100.0, 9.81);
        assert_eq!(model.sloshing_force(), 0.0);
    }
    #[test]
    fn test_sloshing_effective_length_deep() {
        let model = SloshingModel::new(1.0, 1e-5, 10.0, 9.81);
        let l_eff = model.effective_pendulum_length();
        assert!(l_eff.is_finite());
    }
    #[test]
    fn test_sphere_fraction_fully_submerged() {
        let frac = SubmergedVolumeFraction::sphere(0.5, -10.0);
        assert!((frac - 1.0).abs() < 1e-6, "frac={frac}");
    }
    #[test]
    fn test_sphere_fraction_fully_above() {
        let frac = SubmergedVolumeFraction::sphere(0.5, 10.0);
        assert!((frac).abs() < 1e-6, "frac={frac}");
    }
    #[test]
    fn test_sphere_fraction_half_submerged() {
        let frac = SubmergedVolumeFraction::sphere(1.0, 0.0);
        assert!((frac - 0.5).abs() < 0.01, "frac={frac}");
    }
    #[test]
    fn test_box_fraction_fully_below() {
        let frac = SubmergedVolumeFraction::box_fraction(1.0, -2.0);
        assert!((frac - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_box_fraction_fully_above() {
        let frac = SubmergedVolumeFraction::box_fraction(1.0, 1.0);
        assert_eq!(frac, 0.0);
    }
    #[test]
    fn test_box_fraction_partial() {
        let frac = SubmergedVolumeFraction::box_fraction(1.0, -0.5);
        assert!((frac - 0.5).abs() < 1e-10, "frac={frac}");
    }
    #[test]
    fn test_submerged_buoyancy_force() {
        let f = SubmergedVolumeFraction::buoyancy_force(1000.0, 9.81, 0.5, 2.0);
        assert!((f - 9810.0).abs() < 1e-6);
    }
    #[test]
    fn test_reynolds_number_computation() {
        let re = ReynoldsUtils::reynolds(1000.0, 1.0, 0.1, 1e-3);
        assert!((re - 100_000.0).abs() < 1e-6);
    }
    #[test]
    fn test_strouhal_number() {
        let st = ReynoldsUtils::strouhal(2.0, 0.5, 5.0);
        assert!((st - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_keulegan_carpenter_number() {
        let kc = ReynoldsUtils::keulegan_carpenter(2.0, 5.0, 1.0);
        assert!((kc - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_pipe_regime_classification() {
        assert_eq!(ReynoldsUtils::pipe_regime(1000.0), "laminar");
        assert_eq!(ReynoldsUtils::pipe_regime(3000.0), "transitional");
        assert_eq!(ReynoldsUtils::pipe_regime(5000.0), "turbulent");
    }
    #[test]
    fn test_weber_number() {
        let we = ReynoldsUtils::weber(1000.0, 1.0, 0.01, 0.072);
        let expected = 1000.0 * 1.0 * 0.01 / 0.072;
        assert!((we - expected).abs() < 1e-6);
    }
    #[test]
    fn test_euler_number() {
        let eu = ReynoldsUtils::euler(1000.0, 1000.0, 1.0);
        assert!((eu - 2.0).abs() < 1e-10);
    }
}
#[cfg(test)]
mod tests_new_fluid {

    use crate::FluidCoupling;

    use crate::HullResistance;

    use crate::OrientedBuoyancy;
    use crate::PropellerThrust;

    #[test]
    fn test_propeller_thrust_positive_at_zero_advance() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let t = p.thrust(1025.0, 10.0, 0.0);
        assert!(t > 0.0, "thrust at zero advance: {t}");
    }
    #[test]
    fn test_propeller_kt_zero_advance_equals_kt0() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let kt = p.kt(0.0);
        assert!((kt - 0.5).abs() < 1e-12, "kt at J=0: {kt}");
    }
    #[test]
    fn test_propeller_thrust_decreases_with_advance() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let t0 = p.thrust(1025.0, 10.0, 0.0);
        let t1 = p.thrust(1025.0, 10.0, 1.0);
        assert!(t1 < t0, "thrust should decrease with advance velocity");
    }
    #[test]
    fn test_propeller_torque_positive() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let q = p.torque(1025.0, 10.0, 0.0);
        assert!(q > 0.0, "torque: {q}");
    }
    #[test]
    fn test_propeller_advance_ratio_formula() {
        let p = PropellerThrust::new(1.0, 0.5, 0.3, 0.05, 0.03);
        let j = p.advance_ratio(2.0, 4.0);
        assert!((j - 0.5).abs() < 1e-12, "J={j}");
    }
    #[test]
    fn test_propeller_efficiency_zero_at_zero_advance() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let eta = p.open_water_efficiency(0.0, 10.0);
        assert!(eta.abs() < 1e-12, "efficiency at J=0 should be 0: {eta}");
    }
    #[test]
    fn test_propeller_efficiency_positive_at_moderate_advance() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let eta = p.open_water_efficiency(0.5, 10.0);
        assert!(eta >= 0.0, "efficiency should be non-negative: {eta}");
    }
    #[test]
    fn test_propeller_no_thrust_at_zero_rps() {
        let p = PropellerThrust::new(0.3, 0.5, 0.3, 0.05, 0.03);
        let t = p.thrust(1025.0, 0.0, 0.0);
        assert_eq!(t, 0.0, "no thrust at zero rps");
    }
    #[test]
    fn test_hull_cf_ittc57_turbulent() {
        let cf = HullResistance::cf_ittc57(1e7);
        assert!((cf - 0.003).abs() < 1e-10, "cf={cf}");
    }
    #[test]
    fn test_hull_cf_ittc57_zero_re() {
        let cf = HullResistance::cf_ittc57(0.0);
        assert_eq!(cf, 0.0);
    }
    #[test]
    fn test_hull_froude_number_formula() {
        let fn_val = HullResistance::froude_number(2.0, 100.0, 9.81);
        let expected = 2.0 / (9.81 * 100.0_f64).sqrt();
        assert!((fn_val - expected).abs() < 1e-10, "fn={fn_val}");
    }
    #[test]
    fn test_hull_frictional_resistance_positive() {
        let hull = HullResistance::new(500.0, 80.0, 50.0, 0.65, 1.2);
        let rf = hull.frictional_resistance(1025.0, 5.0, 1e8);
        assert!(rf > 0.0, "frictional resistance: {rf}");
    }
    #[test]
    fn test_hull_wave_resistance_zero_at_zero_speed() {
        let hull = HullResistance::new(500.0, 80.0, 50.0, 0.65, 1.2);
        let rw = hull.wave_resistance(100_000.0, 0.0, 9.81);
        assert_eq!(rw, 0.0);
    }
    #[test]
    fn test_hull_total_resistance_greater_than_frictional() {
        let hull = HullResistance::new(500.0, 80.0, 50.0, 0.65, 1.2);
        let rf = hull.frictional_resistance(1025.0, 5.0, 1e8);
        let rt = hull.total_resistance(1025.0, 5.0, 1e8, 1_000_000.0, 9.81);
        assert!(rt >= rf, "total resistance should be >= frictional");
    }
    #[test]
    fn test_oriented_buoyancy_full_volume() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [1.0, 1.0, 1.0]);
        assert!((ob.full_volume() - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_oriented_buoyancy_fully_submerged_force() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [1.0, 1.0, 1.0]);
        let f = ob.buoyancy_force(2.0, [0.0, 1.0, 0.0]);
        let expected = 1000.0 * 9.81 * 4.0 * 2.0;
        assert!((f[1] - expected).abs() < 1.0, "f_y={}", f[1]);
    }
    #[test]
    fn test_oriented_buoyancy_zero_depth_zero_force() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [1.0, 1.0, 1.0]);
        let f = ob.buoyancy_force(0.0, [0.0, 1.0, 0.0]);
        assert_eq!(f[1], 0.0);
    }
    #[test]
    fn test_oriented_buoyancy_partial_submersion() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [1.0, 1.0, 1.0]);
        let f_half = ob.buoyancy_force(1.0, [0.0, 1.0, 0.0]);
        let f_full = ob.buoyancy_force(2.0, [0.0, 1.0, 0.0]);
        assert!(
            (f_half[1] - f_full[1] / 2.0).abs() < 1e-9,
            "half force: {}",
            f_half[1]
        );
    }
    #[test]
    fn test_oriented_buoyancy_metacentric_height_positive_for_stable() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [5.0, 1.0, 2.0]);
        let gm = ob.metacentric_height(1.0, 0.3);
        assert!(gm > 0.0, "GM should be positive for wide hull: {gm}");
    }
    #[test]
    fn test_oriented_buoyancy_force_direction() {
        let ob = OrientedBuoyancy::new(1000.0, 9.81, [1.0, 1.0, 1.0]);
        let f = ob.buoyancy_force(1.0, [1.0, 0.0, 0.0]);
        assert!(f[0] > 0.0, "force should be in x direction");
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn test_morrison_force_pure_drag_no_accel() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_morrison_force(
            0.5,
            10.0,
            2.0,
            1.0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        );
        assert!(f[0] < 0.0, "drag should oppose body motion: f[0]={}", f[0]);
    }
    #[test]
    fn test_morrison_force_pure_inertia_no_relative_vel() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let accel = [1.0, 0.0, 0.0];
        let f =
            fc.compute_morrison_force(0.5, 10.0, 2.0, 1.0, accel, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(
            f[0] > 0.0,
            "inertia force should be in fluid accel direction: {}",
            f[0]
        );
    }
    #[test]
    fn test_morrison_force_zero_velocity_zero_drag() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_morrison_force(
            1.0,
            5.0,
            2.0,
            1.0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn test_added_mass_force_opposes_relative_acceleration() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_added_mass_force(0.5, 0.5, [10.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(
            f[0] < 0.0,
            "added mass force should oppose body accel: {}",
            f[0]
        );
    }
    #[test]
    fn test_added_mass_force_zero_when_accelerations_equal() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let a = [3.0, -2.0, 1.0];
        let f = fc.compute_added_mass_force(0.5, 0.5, a, a);
        assert!(
            f[0].abs() < 1e-12 && f[1].abs() < 1e-12 && f[2].abs() < 1e-12,
            "equal accels → zero force: {:?}",
            f
        );
    }
    #[test]
    fn test_added_mass_magnitude_scales_with_radius() {
        let fc = FluidCoupling::new(1000.0, 9.81);
        let body_accel = [1.0, 0.0, 0.0];
        let f1 = fc.compute_added_mass_force(1.0, 0.5, body_accel, [0.0; 3]);
        let f2 = fc.compute_added_mass_force(2.0, 0.5, body_accel, [0.0; 3]);
        assert!(
            (f2[0].abs() / f1[0].abs() - 8.0).abs() < 1e-6,
            "force ratio should be 8: {}",
            f2[0].abs() / f1[0].abs()
        );
    }
    #[test]
    fn test_viv_zero_outside_lock_in() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_vortex_induced_vibration(0.5, 1.0, 0.001, 0.5, 1.0, [0.0, 1.0, 0.0]);
        assert_eq!(f, [0.0; 3], "no lock-in at very slow flow");
    }
    #[test]
    fn test_viv_nonzero_in_lock_in() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_vortex_induced_vibration(0.5, 1.0, 2.5, 0.5, 0.25, [0.0, 1.0, 0.0]);
        assert!(f[1] > 0.0, "VIV lock-in force should be nonzero: {}", f[1]);
    }
    #[test]
    fn test_viv_zero_for_zero_diameter() {
        let fc = FluidCoupling::new(1025.0, 9.81);
        let f = fc.compute_vortex_induced_vibration(0.0, 1.0, 5.0, 0.5, 1.0, [0.0, 1.0, 0.0]);
        assert_eq!(f, [0.0; 3], "zero diameter → zero force");
    }
}
