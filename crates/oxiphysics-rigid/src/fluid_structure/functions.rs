//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

// Types imported as needed by functions below

/// Euclidean norm of a 3D vector.
pub fn vec3_norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Dot product of two 3D vectors.
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::constraint_forces::vec3_norm;
    use crate::fluid_structure::AddedMass;

    use crate::fluid_structure::AeroelasticCoupling;

    use crate::fluid_structure::BuoyancyForce;

    use crate::fluid_structure::HydrodynamicDrag;

    use crate::fluid_structure::MorisonForce;

    use crate::fluid_structure::SlammingForce;

    use crate::fluid_structure::ThermalFluidForce;

    use crate::fluid_structure::VortexSheddingForce;
    use crate::fluid_structure::WaveLoading;
    use crate::fluid_structure::WindPressure;

    use std::f64::consts::PI;
    #[test]
    fn buoyancy_magnitude_archimedes() {
        let b = BuoyancyForce::new(1000.0, 1.0, [0.0, 0.0, 0.0]);
        assert!((b.buoyancy_magnitude() - 9810.0).abs() < 0.1);
    }
    #[test]
    fn buoyancy_vector_upward() {
        let b = BuoyancyForce::new(1000.0, 0.1, [0.0, 0.0, 0.0]);
        let v = b.buoyancy_vector();
        assert!(v[1] > 0.0);
        assert!(v[0].abs() < 1e-10);
    }
    #[test]
    fn buoyancy_displaced_mass() {
        let b = BuoyancyForce::new(1025.0, 2.0, [0.0, 0.0, 0.0]);
        assert!((b.displaced_mass() - 2050.0).abs() < 0.01);
    }
    #[test]
    fn buoyancy_net_force_positive_for_buoyant_body() {
        let b = BuoyancyForce::new(1000.0, 1.0, [0.0, 0.0, 0.0]);
        assert!(b.net_vertical_force(500.0) > 0.0);
    }
    #[test]
    fn buoyancy_net_force_negative_for_sinking() {
        let b = BuoyancyForce::new(1000.0, 0.1, [0.0, 0.0, 0.0]);
        assert!(b.net_vertical_force(200.0) < 0.0);
    }
    #[test]
    fn drag_force_increases_with_velocity() {
        let d = HydrodynamicDrag::new(1025.0, 1.0, 0.4);
        let f1 = d.drag_force(1.0);
        let f2 = d.drag_force(2.0);
        assert!(f2 > f1);
    }
    #[test]
    fn drag_force_zero_at_zero_velocity() {
        let d = HydrodynamicDrag::new(1025.0, 1.0, 0.4);
        assert!(d.drag_force(0.0) < 1e-10);
    }
    #[test]
    fn drag_vector_opposes_motion() {
        let d = HydrodynamicDrag::new(1025.0, 1.0, 0.4);
        let fv = d.drag_vector([1.0, 0.0, 0.0]);
        assert!(fv[0] < 0.0);
    }
    #[test]
    fn cd_from_re_stokes_flow() {
        let cd = HydrodynamicDrag::cd_from_reynolds(0.1);
        assert!(cd > 100.0);
    }
    #[test]
    fn cd_from_re_turbulent() {
        let cd = HydrodynamicDrag::cd_from_reynolds(1e5);
        assert!((cd - 0.44).abs() < 0.01);
    }
    #[test]
    fn added_mass_sphere() {
        let am = AddedMass::sphere(1000.0, 1.0);
        let expected = (2.0 / 3.0) * PI * 1000.0;
        assert!((am.translational_added_mass() - expected).abs() < 1.0);
    }
    #[test]
    fn added_mass_sphere_coefficient() {
        let am = AddedMass::sphere(1000.0, 0.5);
        assert!((am.added_mass_coefficient() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn added_mass_cylinder_positive() {
        let am = AddedMass::cylinder(1025.0, 0.5, 10.0);
        assert!(am.translational_added_mass() > 0.0);
    }
    #[test]
    fn morison_inertia_force_scales_with_acceleration() {
        let m = MorisonForce::new(1025.0, 1.0, 5.0);
        let f1 = m.inertia_force(1.0);
        let f2 = m.inertia_force(2.0);
        assert!((f2 / f1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn morison_drag_quadratic() {
        let m = MorisonForce::new(1025.0, 1.0, 5.0);
        let f1 = m.drag_force_component(1.0);
        let f2 = m.drag_force_component(2.0);
        assert!((f2 / f1 - 4.0).abs() < 1e-10);
    }
    #[test]
    fn morison_total_force_positive() {
        let m = MorisonForce::new(1025.0, 1.0, 5.0);
        let f = m.total_force(1.0, 1.0);
        assert!(f > 0.0);
    }
    #[test]
    fn wave_peak_frequency() {
        let w = WaveLoading::new(3.0, 10.0, 50.0);
        assert!((w.peak_frequency() - 2.0 * PI / 10.0).abs() < 1e-10);
    }
    #[test]
    fn wave_number_deep_water() {
        let w = WaveLoading::new(3.0, 10.0, 1000.0);
        let k = w.wave_number();
        let omega = w.peak_frequency();
        let k_expected = omega * omega / 9.81;
        assert!((k / k_expected - 1.0).abs() < 0.01);
    }
    #[test]
    fn wave_jonswap_spectrum_positive() {
        let w = WaveLoading::new(3.0, 10.0, 50.0);
        let s = w.jonswap_spectrum(w.peak_frequency());
        assert!(s > 0.0);
    }
    #[test]
    fn vortex_shedding_frequency() {
        let v = VortexSheddingForce::new(2.0, 0.1, 5.0, 1025.0);
        assert!((v.shedding_frequency() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn vortex_peak_transverse_force_positive() {
        let v = VortexSheddingForce::new(2.0, 0.1, 5.0, 1025.0);
        assert!(v.peak_transverse_force() > 0.0);
    }
    #[test]
    fn vortex_lock_in_detection() {
        let v = VortexSheddingForce::new(2.0, 0.1, 5.0, 1025.0);
        let fs = v.shedding_frequency();
        assert!(v.is_locked_in(fs));
        assert!(!v.is_locked_in(fs * 2.0));
    }
    #[test]
    fn slamming_pressure_positive() {
        let s = SlammingForce::new(1025.0, 5.0, 20.0, 5.0, 2.0);
        assert!(s.slamming_pressure() > 0.0);
    }
    #[test]
    fn slamming_force_positive() {
        let s = SlammingForce::new(1025.0, 5.0, 20.0, 5.0, 2.0);
        assert!(s.peak_slamming_force() > 0.0);
    }
    #[test]
    fn slamming_higher_velocity_higher_force() {
        let s1 = SlammingForce::new(1025.0, 2.0, 20.0, 5.0, 1.0);
        let s2 = SlammingForce::new(1025.0, 4.0, 20.0, 5.0, 1.0);
        assert!(s2.peak_slamming_force() > s1.peak_slamming_force());
    }
    #[test]
    fn flutter_speed_positive() {
        let ae = AeroelasticCoupling::new(1.225, 10.0, 1.5, 1e6, 5e5, 20.0, 5.0);
        assert!(ae.flutter_speed() > 0.0);
    }
    #[test]
    fn divergence_speed_positive() {
        let ae = AeroelasticCoupling::new(1.225, 10.0, 1.5, 1e6, 5e5, 20.0, 5.0);
        assert!(ae.divergence_speed(2.0 * PI) > 0.0);
    }
    #[test]
    fn gust_response_greater_than_one() {
        let ae = AeroelasticCoupling::new(1.225, 10.0, 1.5, 1e6, 5e5, 20.0, 5.0);
        let n = ae.gust_response(100.0, 5.73, 5.0, 2000.0);
        assert!(n > 1.0);
    }
    #[test]
    fn wind_speed_increases_with_height() {
        let w = WindPressure::new(30.0);
        let v10 = w.wind_speed_at(10.0);
        let v30 = w.wind_speed_at(30.0);
        assert!(v30 > v10);
    }
    #[test]
    fn wind_pressure_positive() {
        let w = WindPressure::new(30.0);
        let p = w.design_pressure(10.0);
        assert!(p > 0.0);
    }
    #[test]
    fn wind_force_scales_with_area() {
        let w = WindPressure::new(30.0);
        let f1 = w.wind_force(10.0, 1.0);
        let f2 = w.wind_force(10.0, 2.0);
        assert!((f2 / f1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_buoyancy_positive_hot_source() {
        let t = ThermalFluidForce::new(1.2, 293.15, 373.15);
        assert!(t.buoyancy_per_volume() > 0.0);
    }
    #[test]
    fn hot_density_less_than_ambient() {
        let t = ThermalFluidForce::new(1.2, 293.15, 373.15);
        assert!(t.hot_density() < t.ambient_density);
    }
    #[test]
    fn archimedes_buoyancy_positive() {
        let t = ThermalFluidForce::new(1.2, 293.15, 373.15);
        assert!(t.archimedes_buoyancy() > 0.0);
    }
    #[test]
    fn grashof_number_positive() {
        let t = ThermalFluidForce::new(1.2, 293.15, 373.15);
        let gr = t.grashof(1.0, 1.5e-5);
        assert!(gr > 0.0);
    }
    #[test]
    fn plume_rise_positive() {
        let t = ThermalFluidForce::new(1.2, 293.15, 373.15);
        let rise = t.plume_rise(1000.0, 3.0, 100.0);
        assert!(rise > 0.0);
    }
    #[test]
    fn vec3_norm_unit() {
        assert!((vec3_norm([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn vec3_dot_parallel() {
        assert!((vec3_dot([2.0, 0.0, 0.0], [3.0, 0.0, 0.0]) - 6.0).abs() < 1e-10);
    }
}
/// Cross product of two 3D vectors.
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Scale a 3D vector by a scalar.
pub fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Add two 3D vectors.
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Normalise a 3D vector (returns zero vector if near zero).
pub fn vec3_normalise(v: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(v);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Compute hydrostatic pressure at depth `z` below surface (Pa).
pub fn hydrostatic_pressure(fluid_density: f64, gravity: f64, depth: f64) -> f64 {
    fluid_density * gravity * depth.max(0.0)
}
/// Froude number: Fr = U / sqrt(g * L).
pub fn froude_number(velocity: f64, gravity: f64, length: f64) -> f64 {
    velocity / (gravity * length).max(1e-30).sqrt()
}
/// Keulegan-Carpenter number: KC = U_max * T / D.
pub fn keulegan_carpenter(u_max: f64, period: f64, diameter: f64) -> f64 {
    u_max * period / diameter.max(1e-30)
}
/// Cauchy number: Ca = ρ * U² * L³ / EI.
pub fn cauchy_number(fluid_density: f64, velocity: f64, length: f64, ei: f64) -> f64 {
    fluid_density * velocity.powi(2) * length.powi(3) / ei.max(1e-30)
}
/// Strouhal number: St = f * D / U.
pub fn strouhal_number(frequency: f64, diameter: f64, velocity: f64) -> f64 {
    frequency * diameter / velocity.max(1e-30)
}
/// Lock-off reduced velocity range for a cylinder (approximate).
///
/// Returns `(Ur_min, Ur_max)` where lock-in occurs.
pub fn viv_lockin_range(strouhal: f64) -> (f64, f64) {
    let ur_center = 1.0 / strouhal;
    (0.7 * ur_center, 1.3 * ur_center)
}
#[cfg(test)]
mod tests_ext {
    use super::*;
    use crate::constraint_forces::vec3_cross;
    use crate::constraint_forces::vec3_norm;

    use crate::fluid_structure::AddedMassMatrix;

    use crate::fluid_structure::AleDomain;
    use crate::fluid_structure::AleMeshNode;

    use crate::fluid_structure::FlutterAnalysis;
    use crate::fluid_structure::GallopingAnalysis;

    use crate::fluid_structure::HydroelasticBeam;
    use crate::fluid_structure::IbMarker;
    use crate::fluid_structure::ImmersedBoundaryBody;
    use crate::fluid_structure::MonolithicFsiCoupling;

    use crate::fluid_structure::PartitionedFsiCoupling;
    use crate::fluid_structure::RadiationDampingMatrix;
    use crate::fluid_structure::RandomSeaState;

    use crate::fluid_structure::TallBuildingWindResponse;
    use crate::fluid_structure::TankSloshing;

    use crate::fluid_structure::TuneableLiquidColumnDamper;
    use crate::fluid_structure::VivLockIn;

    use crate::fluid_structure::cauchy_number;
    use crate::fluid_structure::hydrostatic_pressure;
    use crate::fluid_structure::keulegan_carpenter;
    use crate::fluid_structure::strouhal_number;
    use crate::fluid_structure::vec3_normalise;
    use crate::fluid_structure::viv_lockin_range;
    use std::f64::consts::PI;
    #[test]
    fn fsi_partitioned_natural_frequency() {
        let c = PartitionedFsiCoupling::new(100.0, 1e4, 50.0);
        assert!(c.natural_frequency() > 0.0);
    }
    #[test]
    fn fsi_partitioned_damping_ratio() {
        let c = PartitionedFsiCoupling::new(100.0, 1e4, 10.0);
        let zeta = c.damping_ratio();
        assert!(zeta > 0.0 && zeta < 1.0);
    }
    #[test]
    fn fsi_partitioned_sub_iterate() {
        let mut c = PartitionedFsiCoupling::new(100.0, 1e4, 50.0);
        let res = c.sub_iterate(500.0, 0.01);
        assert!(res.is_finite());
    }
    #[test]
    fn fsi_monolithic_total_mass() {
        let m = MonolithicFsiCoupling::new(100.0, 1e4, 50.0, 20.0, 5.0, 0.0);
        assert!((m.total_mass() - 120.0).abs() < 1e-9);
    }
    #[test]
    fn fsi_monolithic_steady_state() {
        let mut m = MonolithicFsiCoupling::new(100.0, 1e4, 50.0, 20.0, 5.0, 0.0);
        m.forcing_amplitude = 100.0;
        m.forcing_frequency = 5.0;
        let a = m.steady_state_amplitude();
        assert!(a.is_finite() && a > 0.0);
    }
    #[test]
    fn ib_marker_force() {
        let mut mk = IbMarker::new([0.0, 0.0, 0.0], 1000.0);
        mk.target_position = [0.01, 0.0, 0.0];
        mk.compute_force();
        assert!((mk.ib_force[0] - 10.0).abs() < 1e-9);
    }
    #[test]
    fn ib_body_total_force() {
        let mut body = ImmersedBoundaryBody::new(0.01, 1025.0);
        body.add_marker([0.0; 3], 500.0);
        body.add_marker([0.1, 0.0, 0.0], 500.0);
        body.markers[0].target_position = [0.005, 0.0, 0.0];
        body.markers[1].target_position = [0.095, 0.0, 0.0];
        body.compute_forces();
        let f = body.total_force();
        assert!(f[0].abs() < 1.0);
    }
    #[test]
    fn ib_delta_function_origin() {
        let d = ImmersedBoundaryBody::delta_function(0.0);
        assert!(d > 0.0);
    }
    #[test]
    fn ib_delta_function_far() {
        let d = ImmersedBoundaryBody::delta_function(2.0);
        assert!(d.abs() < 1e-10);
    }
    #[test]
    fn ale_node_advance() {
        let mut node = AleMeshNode::new([0.0, 0.0, 0.0]);
        node.mesh_velocity = [1.0, 0.0, 0.0];
        node.advance_mesh(0.1);
        assert!((node.current_position[0] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn ale_convective_velocity() {
        let mut node = AleMeshNode::new([0.0; 3]);
        node.material_velocity = [3.0, 0.0, 0.0];
        node.mesh_velocity = [1.0, 0.0, 0.0];
        let cv = node.convective_velocity();
        assert!((cv[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn ale_domain_max_velocity() {
        let mut dom = AleDomain::new(1);
        dom.add_node([0.0; 3]);
        dom.add_node([1.0, 0.0, 0.0]);
        dom.nodes[1].mesh_velocity = [2.0, 0.0, 0.0];
        assert!((dom.max_mesh_velocity() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn added_mass_matrix_sphere_trace() {
        let am = AddedMassMatrix::sphere(1000.0, 1.0);
        let expected = 3.0 * (2.0 / 3.0) * PI * 1000.0;
        assert!((am.trace() - expected).abs() < 1.0);
    }
    #[test]
    fn added_mass_matrix_ellipsoid_diagonal() {
        let am = AddedMassMatrix::ellipsoid(1000.0, 2.0, 1.0, 0.5);
        let diag = am.translational_diagonal();
        for d in diag {
            assert!(d >= 0.0);
        }
    }
    #[test]
    fn radiation_damping_power_positive() {
        let mut rd = RadiationDampingMatrix::new(1.0);
        rd.set_diagonal([100.0, 200.0, 150.0, 10.0, 20.0, 15.0]);
        let vel = [1.0, 0.5, 0.3, 0.1, 0.05, 0.02];
        let p = rd.dissipated_power(vel);
        assert!(p >= 0.0);
    }
    #[test]
    fn random_sea_jonswap_positive() {
        let sea = RandomSeaState::new(3.0, 10.0, 50.0);
        let s = sea.jonswap(sea.omega_p());
        assert!(s > 0.0);
    }
    #[test]
    fn random_sea_hs_from_spectrum() {
        let sea = RandomSeaState::new(3.0, 10.0, 100.0);
        let hs_calc = sea.hs_from_spectrum(0.1, 5.0, 500);
        assert!(
            hs_calc > 1.5 && hs_calc < 6.0,
            "hs_calc={hs_calc} outside expected range [1.5, 6.0]"
        );
    }
    #[test]
    fn random_sea_max_wave() {
        let sea = RandomSeaState::new(3.0, 10.0, 50.0);
        let hmax = sea.most_probable_max_wave(1000.0);
        assert!(hmax > sea.hs);
    }
    #[test]
    fn random_sea_synthesise_count() {
        let sea = RandomSeaState::new(2.0, 8.0, 30.0);
        let ts = sea.synthesise(0.0, 10.0, 0.5, 0.1, 3.0);
        assert_eq!(ts.len(), 20);
    }
    #[test]
    fn random_sea_rms_morison() {
        let sea = RandomSeaState::new(2.0, 10.0, 50.0);
        let f = sea.rms_morison_force(1.0, 2.0, 1.0);
        assert!(f > 0.0);
    }
    #[test]
    fn galloping_critical_speed_positive() {
        let g = GallopingAnalysis::new(1.2, 10.0, 0.1, 5.0, 0.01, 2.0 * PI, -2.5);
        let uc = g.critical_wind_speed();
        assert!(uc > 0.0);
    }
    #[test]
    fn galloping_instability_check() {
        let g = GallopingAnalysis::new(1.2, 20.0, 0.1, 5.0, 0.001, 2.0 * PI, -10.0);
        assert!(g.is_galloping_unstable());
    }
    #[test]
    fn galloping_stable_positive_coeff() {
        let g = GallopingAnalysis::new(1.2, 10.0, 0.1, 5.0, 0.01, 2.0 * PI, 1.0);
        assert!(!g.is_galloping_unstable());
    }
    #[test]
    fn flutter_speed_positive() {
        let fa = FlutterAnalysis::new(1e6, 5e5, 20.0, 5.0, 0.75, 10.0, 1.225, 0.1);
        assert!(fa.flutter_speed() > 0.0);
    }
    #[test]
    fn flutter_frequency_ratio() {
        let fa = FlutterAnalysis::new(1e5, 5e5, 20.0, 5.0, 0.75, 10.0, 1.225, 0.0);
        let fr = fa.frequency_ratio();
        assert!(fr.is_finite() && fr > 0.0);
    }
    #[test]
    fn flutter_reduced_speed_gt_one() {
        let fa = FlutterAnalysis::new(1e6, 5e5, 20.0, 5.0, 0.75, 10.0, 1.225, 0.0);
        let vfstar = fa.reduced_flutter_speed();
        assert!(vfstar > 1.0);
    }
    #[test]
    fn viv_scruton_number() {
        let viv = VivLockIn::new(0.1, 100.0, 0.01, 2.0, 4.0);
        assert!((viv.scruton_number() - 1.0).abs() < 1e-9);
    }
    #[test]
    fn viv_max_amplitude_ratio() {
        let viv = VivLockIn::new(0.1, 100.0, 0.01, 2.0, 4.0);
        assert!(viv.max_amplitude_ratio() > 0.0);
    }
    #[test]
    fn viv_lock_in_at_ur5() {
        let viv = VivLockIn::new(0.1, 100.0, 0.01, 2.0, 4.0);
        assert!(viv.is_locked_in());
    }
    #[test]
    fn viv_not_locked_far() {
        let viv = VivLockIn::new(0.1, 100.0, 0.01, 20.0, 4.0);
        assert!(!viv.is_locked_in());
    }
    #[test]
    fn sloshing_natural_freq_positive() {
        let s = TankSloshing::new(2.0, 0.5, 1.0, 1000.0);
        assert!(s.natural_frequency() > 0.0);
    }
    #[test]
    fn sloshing_pendulum_length() {
        let s = TankSloshing::new(2.0, 0.5, 1.0, 1000.0);
        let omega = s.natural_frequency();
        let l_eq = s.equivalent_pendulum_length();
        assert!((l_eq - s.gravity / omega.powi(2)).abs() < 1e-9);
    }
    #[test]
    fn sloshing_effective_mass_lt_total() {
        let s = TankSloshing::new(2.0, 0.3, 1.0, 1000.0);
        assert!(s.effective_mass() < s.total_liquid_mass());
    }
    #[test]
    fn tlcd_natural_frequency() {
        let t = TuneableLiquidColumnDamper::new(2.0, 1.0, 0.01, 1000.0, 2.0);
        let omega = t.natural_frequency();
        let expected = (2.0_f64 * 9.81 / 2.0).sqrt();
        assert!((omega - expected).abs() < 1e-6);
    }
    #[test]
    fn hydrobeam_natural_freq() {
        let b = HydroelasticBeam::new(10.0, 1e6, 100.0, 1025.0, 0.5);
        assert!(b.natural_frequency_mode(1) > 0.0);
    }
    #[test]
    fn hydrobeam_tip_deflection_increases() {
        let mut b1 = HydroelasticBeam::new(10.0, 1e6, 100.0, 1025.0, 0.5);
        let mut b2 = HydroelasticBeam::new(10.0, 1e6, 100.0, 1025.0, 0.5);
        b1.fluid_velocity = 1.0;
        b2.fluid_velocity = 2.0;
        assert!(b2.static_tip_deflection() > b1.static_tip_deflection());
    }
    #[test]
    fn hydrobeam_added_mass_coeff() {
        let b = HydroelasticBeam::new(10.0, 1e6, 100.0, 1025.0, 0.5);
        assert!((b.added_mass_coefficient() - 1.0).abs() < 1e-6);
    }
    #[test]
    fn building_mean_along_wind() {
        let b = TallBuildingWindResponse::new(200.0, 40.0, 40.0, 0.2, 0.02, 40.0, 1e7);
        assert!(b.mean_along_wind_force() > 0.0);
    }
    #[test]
    fn building_gust_factor() {
        let b = TallBuildingWindResponse::new(200.0, 40.0, 40.0, 0.2, 0.02, 40.0, 1e7);
        let g = b.gust_factor(0.15);
        assert!(g > 1.0);
    }
    #[test]
    fn building_across_wind_displacement() {
        let b = TallBuildingWindResponse::new(200.0, 40.0, 40.0, 0.2, 0.02, 40.0, 1e7);
        let d = b.peak_across_wind_displacement(0.1);
        assert!(d.is_finite() && d > 0.0);
    }
    #[test]
    fn building_shedding_frequency() {
        let b = TallBuildingWindResponse::new(200.0, 40.0, 40.0, 0.2, 0.02, 40.0, 1e7);
        assert!(b.shedding_frequency() > 0.0);
    }
    #[test]
    fn building_resonance_check() {
        let b = TallBuildingWindResponse::new(200.0, 40.0, 40.0, 0.2, 0.02, 40.0, 1e7);
        let _ = b.is_across_wind_resonant();
    }
    #[test]
    fn helper_vec3_cross() {
        let c = vec3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-10);
        assert!(c[0].abs() < 1e-10 && c[1].abs() < 1e-10);
    }
    #[test]
    fn helper_vec3_normalise() {
        let v = vec3_normalise([3.0, 4.0, 0.0]);
        assert!((vec3_norm(v) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn helper_froude_number() {
        let fr = froude_number(5.0, 9.81, 2.0);
        assert!(fr > 0.0);
    }
    #[test]
    fn helper_keulegan_carpenter() {
        let kc = keulegan_carpenter(2.0, 10.0, 0.5);
        assert!((kc - 40.0).abs() < 1e-9);
    }
    #[test]
    fn helper_strouhal_number() {
        let st = strouhal_number(4.0, 0.1, 2.0);
        assert!((st - 0.2).abs() < 1e-10);
    }
    #[test]
    fn helper_viv_lockin_range() {
        let (lo, hi) = viv_lockin_range(0.2);
        assert!(lo < hi);
    }
    #[test]
    fn helper_hydrostatic_pressure() {
        let p1 = hydrostatic_pressure(1025.0, 9.81, 10.0);
        let p2 = hydrostatic_pressure(1025.0, 9.81, 20.0);
        assert!(p2 > p1);
    }
    #[test]
    fn helper_cauchy_number() {
        let ca = cauchy_number(1025.0, 2.0, 5.0, 1e6);
        assert!(ca.is_finite() && ca > 0.0);
    }
}
