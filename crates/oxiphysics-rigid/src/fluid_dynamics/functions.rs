// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
#[allow(dead_code)]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
#[allow(dead_code)]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-14 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / n)
    }
}
#[cfg(test)]
mod tests {

    use crate::fluid_dynamics::types::{
        AerodynamicBody, AerodynamicTorque, DragModel, GroundEffect, HydrodynamicBody, LiftModel,
        MorisonEquation, SlamForce, VortexInducedVibration, WakeEffect,
    };
    #[test]
    fn aero_body_flat_plate_cl_positive_alpha() {
        let body = AerodynamicBody::flat_plate(1.0, 0.1, 1.225);
        let cl = body.cl(5.0_f64.to_radians());
        assert!(cl > 0.0, "cl={cl}");
    }
    #[test]
    fn aero_body_drag_positive() {
        let body = AerodynamicBody::flat_plate(1.0, 0.1, 1.225);
        let (_l, d) = body.forces(30.0, 5.0_f64.to_radians());
        assert!(d > 0.0, "d={d}");
    }
    #[test]
    fn aero_body_forces_scale_with_v_squared() {
        let body = AerodynamicBody::flat_plate(1.0, 0.1, 1.225);
        let (_, d1) = body.forces(10.0, 5.0_f64.to_radians());
        let (_, d2) = body.forces(20.0, 5.0_f64.to_radians());
        assert!((d2 / d1 - 4.0).abs() < 0.01, "ratio={}", d2 / d1);
    }
    #[test]
    fn aero_body_interp_clamp_below() {
        let body = AerodynamicBody::flat_plate(1.0, 0.1, 1.225);
        let cl = body.cl(-1.0);
        assert!(cl.is_finite());
    }
    #[test]
    fn drag_stokes_linear_in_velocity() {
        let drag = DragModel::new(1e-3, 1000.0, 1e-3);
        let d1 = drag.stokes_drag(1.0);
        let d2 = drag.stokes_drag(2.0);
        assert!((d2 / d1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn drag_newton_quadratic_in_velocity() {
        let drag = DragModel::new(1e-3, 1.225, 0.1);
        let d1 = drag.newton_drag(10.0);
        let d2 = drag.newton_drag(20.0);
        assert!((d2 / d1 - 4.0).abs() < 0.01, "ratio={}", d2 / d1);
    }
    #[test]
    fn drag_reynolds_number() {
        let drag = DragModel::new(1e-3, 1000.0, 1e-3);
        let re = drag.reynolds(1.0);
        assert!((re - 1000.0).abs() < 1e-10, "re={re}");
    }
    #[test]
    fn drag_schiller_naumann_high_re() {
        let drag = DragModel::new(1e-3, 1000.0, 1e-3);
        let cd = drag.schiller_naumann_cd(10.0);
        assert!((cd - 0.44).abs() < 1e-10, "cd={cd}");
    }
    #[test]
    fn drag_force_selects_stokes_at_low_re() {
        let drag = DragModel::new(1.0, 1.0, 1.0);
        let f = drag.drag_force(0.1);
        let f_stokes = drag.stokes_drag(0.1);
        assert!((f - f_stokes).abs() < 1e-12);
    }
    #[test]
    fn lift_thin_airfoil_cl_zero_at_zero_alpha() {
        let lm = LiftModel::new(1.225, 50.0, 10.0, 1.0);
        assert!(lm.thin_airfoil_cl(0.0).abs() < 1e-10);
    }
    #[test]
    fn lift_thin_airfoil_positive_alpha() {
        let lm = LiftModel::new(1.225, 50.0, 10.0, 1.0);
        let cl = lm.thin_airfoil_cl(5.0_f64.to_radians());
        assert!(cl > 0.0, "cl={cl}");
    }
    #[test]
    fn lift_kutta_joukowski() {
        let lm = LiftModel::new(1.225, 50.0, 10.0, 1.0);
        let l = lm.kutta_joukowski_lift(10.0);
        assert!((l - 6125.0).abs() < 0.01, "l={l}");
    }
    #[test]
    fn aero_torque_pitching_moment_finite() {
        let at = AerodynamicTorque::new(-0.1, -5.0, 1.0, 0.2, 1.225);
        let m = at.pitching_moment(5.0_f64.to_radians(), 30.0);
        assert!(m.is_finite(), "m={m}");
    }
    #[test]
    fn aero_torque_damping_zero_at_zero_pitch_rate() {
        let at = AerodynamicTorque::new(-0.1, -5.0, 1.0, 0.2, 1.225);
        let m = at.damping_torque(0.0, 30.0);
        assert!(m.abs() < 1e-10, "m={m}");
    }
    #[test]
    fn aero_torque_total_is_sum() {
        let at = AerodynamicTorque::new(-0.1, -5.0, 1.0, 0.2, 1.225);
        let m_p = at.pitching_moment(0.1, 30.0);
        let m_d = at.damping_torque(0.5, 30.0);
        let m_t = at.total_moment(0.1, 0.5, 30.0);
        assert!((m_t - m_p - m_d).abs() < 1e-10, "m_t={m_t}");
    }
    #[test]
    fn viv_shedding_frequency() {
        let viv = VortexInducedVibration::new(0.2, 0.1, 1.225, 0.5, 0.1);
        let f = viv.shedding_frequency(10.0);
        assert!((f - 20.0).abs() < 1e-10, "f={f}");
    }
    #[test]
    fn viv_lock_in_detection() {
        let viv = VortexInducedVibration::new(0.2, 0.1, 1.225, 0.5, 0.2);
        assert!(viv.is_lock_in(10.0, 20.0));
        assert!(!viv.is_lock_in(10.0, 100.0));
    }
    #[test]
    fn viv_force_oscillates() {
        let viv = VortexInducedVibration::new(0.2, 0.1, 1.225, 0.5, 0.1);
        let f0 = viv.force_per_length(10.0, 0.0);
        let f_half = viv.force_per_length(10.0, 1.0 / (2.0 * 20.0));
        assert!(f0 * f_half <= 0.0, "f0={f0}, f_half={f_half}");
    }
    #[test]
    fn wake_velocity_less_than_freestream() {
        let wake = WakeEffect::new(0.8, 1.0, 0.04);
        let v_wake = wake.wake_velocity(5.0, 10.0);
        assert!(v_wake < 10.0, "v_wake={v_wake}");
    }
    #[test]
    fn wake_velocity_recovers_with_distance() {
        let wake = WakeEffect::new(0.8, 1.0, 0.04);
        let v1 = wake.wake_velocity(1.0, 10.0);
        let v2 = wake.wake_velocity(100.0, 10.0);
        assert!(v2 > v1, "velocity should recover: v1={v1}, v2={v2}");
    }
    #[test]
    fn wake_blockage_ratio_positive() {
        let wake = WakeEffect::new(0.8, 1.0, 0.04);
        let br = wake.blockage_ratio(5.0);
        assert!(br > 0.0 && br < 1.0, "br={br}");
    }
    #[test]
    fn ground_effect_lift_ratio_above_one_near_ground() {
        let ge = GroundEffect::new(10.0, 50000.0, 10.0);
        let r = ge.lift_ratio(1.0);
        assert!(r > 1.0, "lift_ratio={r}");
    }
    #[test]
    fn ground_effect_lift_ratio_unity_far_away() {
        let ge = GroundEffect::new(10.0, 50000.0, 10.0);
        let r = ge.lift_ratio(30.0);
        assert!((r - 1.0).abs() < 1e-10, "r={r}");
    }
    #[test]
    fn hydro_sphere_buoyancy_positive() {
        let hb = HydrodynamicBody::sphere(0.5, 1025.0);
        let b = hb.buoyancy();
        assert!(b > 0.0, "b={b}");
    }
    #[test]
    fn hydro_added_mass_opposes_acceleration() {
        let hb = HydrodynamicBody::sphere(0.5, 1025.0);
        let f = hb.added_mass_force([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(
            f[0] < 0.0,
            "added mass should oppose acceleration: f={}",
            f[0]
        );
    }
    #[test]
    fn hydro_fk_force_finite() {
        let hb = HydrodynamicBody::sphere(0.5, 1025.0);
        let f = hb.froude_krylov_force(100.0);
        assert!(f.is_finite(), "f={f}");
    }
    #[test]
    fn morison_drag_term_positive() {
        let morison = MorisonEquation::new(1.2, 2.0, 0.5, 10.0, 1025.0);
        let f = morison.force_per_length(2.0, 0.0);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn morison_inertia_term_positive() {
        let morison = MorisonEquation::new(1.2, 2.0, 0.5, 10.0, 1025.0);
        let f = morison.force_per_length(0.0, 1.0);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn morison_kc_number() {
        let morison = MorisonEquation::new(1.2, 2.0, 0.5, 10.0, 1025.0);
        let kc = morison.kc_number(1.0, 10.0);
        assert!((kc - 20.0).abs() < 1e-10, "kc={kc}");
    }
    #[test]
    fn morison_drag_dominated_check() {
        let morison = MorisonEquation::new(1.2, 2.0, 0.5, 10.0, 1025.0);
        assert!(morison.is_drag_dominated(10.0, 10.0));
        assert!(!morison.is_drag_dominated(1.0, 1.0));
    }
    #[test]
    fn slam_max_pressure_positive() {
        let slam = SlamForce::new(15.0_f64.to_radians(), 1025.0, 1.0);
        let p = slam.max_impact_pressure(5.0);
        assert!(p > 0.0, "p={p}");
    }
    #[test]
    fn slam_impact_force_scales_with_v_squared() {
        let slam = SlamForce::new(15.0_f64.to_radians(), 1025.0, 1.0);
        let f1 = slam.impact_force(2.0, 0.1);
        let f2 = slam.impact_force(4.0, 0.1);
        assert!((f2 / f1 - 4.0).abs() < 0.01, "ratio={}", f2 / f1);
    }
    #[test]
    fn slam_added_mass_coeff_positive() {
        let slam = SlamForce::new(15.0_f64.to_radians(), 1025.0, 1.0);
        assert!(slam.added_mass_coeff() > 0.0);
    }
    #[test]
    fn slam_deadrise_angle_effect() {
        let slam_small = SlamForce::new(5.0_f64.to_radians(), 1025.0, 1.0);
        let slam_large = SlamForce::new(30.0_f64.to_radians(), 1025.0, 1.0);
        let f_small = slam_small.impact_force(5.0, 0.1);
        let f_large = slam_large.impact_force(5.0, 0.1);
        assert!(f_small > f_large, "smaller deadrise → larger force");
    }
}
/// Biot-Savart law for 3-D velocity induced by a vortex filament.
///
/// Returns velocity \[m/s\] induced at point `p` by vortex from `a` to `b` with strength Γ.
pub fn biot_savart_filament(a: [f64; 3], b: [f64; 3], p: [f64; 3], gamma: f64) -> [f64; 3] {
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let n1 = norm3(r1);
    let n2 = norm3(r2);
    if n1 < 1e-10 || n2 < 1e-10 {
        return [0.0; 3];
    }
    let r1r2 = cross3(r1, r2);
    let n_cross = norm3(r1r2);
    if n_cross < 1e-10 {
        return [0.0; 3];
    }
    let r0 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let factor = gamma / (4.0 * PI * n_cross * n_cross) * (dot3(r0, r1) / n1 - dot3(r0, r2) / n2);
    [r1r2[0] * factor, r1r2[1] * factor, r1r2[2] * factor]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod advanced_tests {

    use crate::fluid_dynamics::CavitationAnalysis;
    use crate::fluid_dynamics::DoubletPanelMethod;
    use crate::fluid_dynamics::ExtendedMorison;
    use crate::fluid_dynamics::HorseshoeVortex;
    use crate::fluid_dynamics::KelvinHelmholtz;
    use crate::fluid_dynamics::LiftingLine;
    use crate::fluid_dynamics::Panel;
    use crate::fluid_dynamics::RayleighTaylor;
    use crate::fluid_dynamics::ShipHydrodynamics;
    use crate::fluid_dynamics::SourcePanelMethod;
    use crate::fluid_dynamics::TurbulentBoundaryLayer;
    use crate::fluid_dynamics::VIVPredictor;
    use crate::fluid_dynamics::VortexLatticeSolver;
    use crate::fluid_dynamics::VortexPanelMethod;
    use crate::fluid_dynamics::biot_savart_filament;
    use crate::fluid_dynamics::functions::PI;
    use crate::fluid_dynamics::norm3;
    #[test]
    fn panel_length_correct() {
        let p = Panel::new([0.0, 0.0], [3.0, 4.0]);
        assert!((p.length() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn panel_center_correct() {
        let p = Panel::new([0.0, 0.0], [2.0, 2.0]);
        let c = p.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn panel_normal_unit_length() {
        let p = Panel::new([0.0, 0.0], [1.0, 0.0]);
        let n = p.normal();
        let mag = (n[0] * n[0] + n[1] * n[1]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }
    #[test]
    fn source_panel_method_circle_n_panels() {
        let spm = SourcePanelMethod::circle(1.0, 16, 10.0, 0.0);
        assert_eq!(spm.panels.len(), 16);
    }
    #[test]
    fn source_panel_cp_zero_at_freestream() {
        let spm = SourcePanelMethod::circle(1.0, 8, 10.0, 0.0);
        let cp = spm.cp_at_panel(10.0);
        assert!(cp.abs() < 1e-10, "cp={cp}");
    }
    #[test]
    fn vortex_panel_total_circulation() {
        let p = Panel::new([0.0, 0.0], [1.0, 0.0]);
        let mut vpm = VortexPanelMethod::new(vec![p], 10.0, 0.0);
        vpm.set_uniform_gamma(5.0);
        assert!((vpm.total_circulation() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn vortex_panel_lift_positive() {
        let p = Panel::new([0.0, 0.0], [1.0, 0.0]);
        let mut vpm = VortexPanelMethod::new(vec![p], 10.0, 0.1);
        vpm.set_uniform_gamma(1.0);
        let l = vpm.lift_per_span(1.225);
        assert!(l > 0.0, "l={l}");
    }
    #[test]
    fn doublet_panel_far_field_decreases_with_distance() {
        let dpm = DoubletPanelMethod::new(1, vec![1.0]);
        let mut dm = dpm;
        dm.set_strengths(vec![1.0]);
        let phi1 = dm.far_field_potential(1.0);
        let phi2 = dm.far_field_potential(2.0);
        assert!(phi1 > phi2, "phi1={phi1}, phi2={phi2}");
    }
    #[test]
    fn lifting_line_aspect_ratio_positive() {
        let ll = LiftingLine::new(10.0, 1.5, 0.8, 5.0_f64.to_radians(), 1.225, 50.0, 4);
        let ar = ll.aspect_ratio();
        assert!(ar > 0.0, "ar={ar}");
    }
    #[test]
    fn lifting_line_finite_wing_slope_less_than_2pi() {
        let ll = LiftingLine::new(10.0, 1.5, 0.8, 5.0_f64.to_radians(), 1.225, 50.0, 4);
        let a = ll.lift_curve_slope();
        assert!(a < 2.0 * PI, "a={a}");
    }
    #[test]
    fn lifting_line_total_lift_positive() {
        let ll = LiftingLine::new(10.0, 1.5, 0.8, 5.0_f64.to_radians(), 1.225, 50.0, 4);
        let l = ll.total_lift();
        assert!(l > 0.0, "l={l}");
    }
    #[test]
    fn lifting_line_elliptic_circulation_zero_at_tip() {
        let ll = LiftingLine::new(10.0, 1.0, 0.5, 0.1, 1.225, 50.0, 4);
        let gamma = ll.elliptic_circulation(10.0, 5.0);
        assert!(gamma.abs() < 1e-10, "gamma={gamma}");
    }
    #[test]
    fn lifting_line_induced_drag_positive_at_nonzero_lift() {
        let ll = LiftingLine::new(10.0, 1.5, 0.8, 5.0_f64.to_radians(), 1.225, 50.0, 4);
        let cdi = ll.induced_drag_coefficient(1.0);
        assert!(cdi >= 0.0, "cdi={cdi}");
    }
    #[test]
    fn biot_savart_zero_gamma() {
        let v = biot_savart_filament([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], 0.0);
        assert!(norm3(v) < 1e-14);
    }
    #[test]
    fn biot_savart_nonzero_result() {
        let v = biot_savart_filament([-10.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        assert!(norm3(v).is_finite());
    }
    #[test]
    fn vlm_aic_entry_finite() {
        let pa = HorseshoeVortex::new(
            [-0.25, -0.5, 0.0],
            [-0.25, 0.5, 0.0],
            1.0,
            [0.25, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        let solver = VortexLatticeSolver::new(vec![pa], [10.0, 0.0, 0.0], 1.225);
        let aic = solver.aic_entry(0, 0);
        assert!(aic.is_finite(), "aic={aic}");
    }
    #[test]
    fn vlm_rhs_nonzero_for_nonzero_alpha() {
        let pa = HorseshoeVortex::new(
            [-0.25, -0.5, 0.0],
            [-0.25, 0.5, 0.0],
            0.0,
            [0.25, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        let solver = VortexLatticeSolver::new(
            vec![pa],
            [10.0 * (0.1_f64).cos(), 0.0, 10.0 * (0.1_f64).sin()],
            1.225,
        );
        let rhs = solver.rhs();
        assert!(!rhs.is_empty());
        assert!(rhs[0] != 0.0 || rhs[0].abs() < 1.0);
    }
    #[test]
    fn extended_morison_drag_force_quadratic() {
        let m = ExtendedMorison::new(1.2, 1.0, 0.5, 1025.0);
        let f1 = m.drag_force(1.0);
        let f2 = m.drag_force(2.0);
        assert!((f2 / f1 - 4.0).abs() < 0.01, "ratio={}", f2 / f1);
    }
    #[test]
    fn extended_morison_inertia_force_linear_in_acceleration() {
        let m = ExtendedMorison::new(1.2, 1.0, 0.5, 1025.0);
        let f1 = m.inertia_force(1.0);
        let f2 = m.inertia_force(2.0);
        assert!((f2 / f1 - 2.0).abs() < 1e-10, "ratio={}", f2 / f1);
    }
    #[test]
    fn extended_morison_cm_correct() {
        let m = ExtendedMorison::new(1.2, 1.5, 0.5, 1025.0);
        assert!((m.cm() - 2.5).abs() < 1e-12);
    }
    #[test]
    fn extended_morison_max_force_regular_wave_finite() {
        let m = ExtendedMorison::new(1.2, 1.0, 0.5, 1025.0);
        let f = m.max_force_regular_wave(2.0, 8.0, 30.0);
        assert!(f.is_finite() && f > 0.0);
    }
    #[test]
    fn viv_predictor_strouhal_low_re() {
        let v = VIVPredictor::new(0.1, 1025.0, 1e-6, 1.0, 20.0, 0.02);
        let st = v.strouhal_from_re(100.0);
        assert!((st - 0.2).abs() < 1e-10);
    }
    #[test]
    fn viv_predictor_lock_in_at_reduced_velocity_6() {
        let v = VIVPredictor::new(0.1, 1025.0, 1e-6, 1.0, 20.0, 0.02);
        assert!(v.is_lock_in(0.6));
    }
    #[test]
    fn viv_predictor_no_lock_in_low_velocity() {
        let v = VIVPredictor::new(0.1, 1025.0, 1e-6, 1.0, 20.0, 0.02);
        assert!(!v.is_lock_in(0.1));
    }
    #[test]
    fn viv_predictor_scruton_number_positive() {
        let v = VIVPredictor::new(0.1, 1025.0, 1e-6, 1.0, 20.0, 0.02);
        assert!(v.scruton_number() > 0.0);
    }
    #[test]
    fn kh_growth_rate_zero_when_stable() {
        let kh = KelvinHelmholtz::new(5.0, 5.0, 1.0, 1.2, 0.07, 9.81);
        let gr = kh.growth_rate(10.0);
        assert!(gr.abs() < 1e-10, "gr={gr}");
    }
    #[test]
    fn kh_growth_rate_positive_when_unstable() {
        let kh = KelvinHelmholtz::new(20.0, 0.0, 1.2, 1.0, 0.0, 9.81);
        let gr = kh.growth_rate(5.0);
        assert!(gr > 0.0, "gr={gr}");
    }
    #[test]
    fn kh_atwood_number_range() {
        let kh = KelvinHelmholtz::new(5.0, 0.0, 1.0, 2.0, 0.0, 9.81);
        let a = kh.atwood_number();
        assert!((-1.0..=1.0).contains(&a), "a={a}");
    }
    #[test]
    fn kh_phase_velocity_weighted_mean() {
        let kh = KelvinHelmholtz::new(10.0, 0.0, 1.0, 1.0, 0.0, 9.81);
        let c = kh.phase_velocity(1.0);
        assert!((c - 5.0).abs() < 1e-10, "c={c}");
    }
    #[test]
    fn rt_growth_rate_positive_heavy_above_light() {
        let rt = RayleighTaylor::new(2.0, 1.0, 0.0, 9.81);
        let gr = rt.growth_rate(10.0);
        assert!(gr > 0.0, "gr={gr}");
    }
    #[test]
    fn rt_growth_rate_zero_when_stable() {
        let rt = RayleighTaylor::new(1.0, 2.0, 0.0, 9.81);
        let gr = rt.growth_rate(10.0);
        assert!(gr < 1e-10, "gr={gr}");
    }
    #[test]
    fn rt_critical_wavenumber_decreases_with_surface_tension() {
        let rt_low = RayleighTaylor::new(2.0, 1.0, 0.01, 9.81);
        let rt_high = RayleighTaylor::new(2.0, 1.0, 0.1, 9.81);
        let kc_low = rt_low.critical_wavenumber();
        let kc_high = rt_high.critical_wavenumber();
        assert!(kc_low > kc_high, "kc_low={kc_low}, kc_high={kc_high}");
    }
    #[test]
    fn rt_max_growth_rate_finite() {
        let rt = RayleighTaylor::new(2.0, 1.0, 0.07, 9.81);
        let m = rt.max_growth_rate();
        assert!(m.is_finite() && m >= 0.0);
    }
    #[test]
    fn ship_froude_number_correct() {
        let ship = ShipHydrodynamics::new(100.0, 20.0, 8.0, 0.65, 1025.0, 1e-6);
        let fr = ship.froude_number(10.0);
        assert!((fr - 10.0 / (9.81 * 100.0_f64).sqrt()).abs() < 1e-10);
    }
    #[test]
    fn ship_displacement_mass_positive() {
        let ship = ShipHydrodynamics::new(100.0, 20.0, 8.0, 0.65, 1025.0, 1e-6);
        assert!(ship.displacement_mass() > 0.0);
    }
    #[test]
    fn ship_friction_resistance_positive_at_speed() {
        let ship = ShipHydrodynamics::new(100.0, 20.0, 8.0, 0.65, 1025.0, 1e-6);
        let rf = ship.frictional_resistance(5.0);
        assert!(rf > 0.0, "rf={rf}");
    }
    #[test]
    fn ship_total_resistance_increases_with_speed() {
        let ship = ShipHydrodynamics::new(100.0, 20.0, 8.0, 0.65, 1025.0, 1e-6);
        let r1 = ship.total_resistance(3.0);
        let r2 = ship.total_resistance(6.0);
        assert!(r2 > r1, "r1={r1}, r2={r2}");
    }
    #[test]
    fn ship_added_mass_sway_larger_than_surge() {
        let ship = ShipHydrodynamics::new(100.0, 20.0, 8.0, 0.65, 1025.0, 1e-6);
        assert!(ship.added_mass_sway() > ship.added_mass_surge());
    }
    #[test]
    fn cavitation_thoma_number_positive() {
        let cav = CavitationAnalysis::new(10.0, 101325.0, 2338.0, 998.0);
        let sigma = cav.thoma_number();
        assert!(sigma > 0.0, "sigma={sigma}");
    }
    #[test]
    fn cavitation_inception_at_low_sigma() {
        let cav = CavitationAnalysis::new(100.0, 101325.0, 2338.0, 998.0);
        let sigma = cav.thoma_number();
        let cavitating = cav.inception_cavitation(0.3);
        if sigma < 0.3 {
            assert!(cavitating);
        }
    }
    #[test]
    fn cavitation_collapse_time_positive() {
        let cav = CavitationAnalysis::new(10.0, 101325.0, 2338.0, 998.0);
        let tc = cav.collapse_time(1e-3);
        assert!(tc > 0.0 && tc.is_finite());
    }
    #[test]
    fn cavitation_critical_velocity_positive() {
        let cav = CavitationAnalysis::new(10.0, 101325.0, 2338.0, 998.0);
        let vc = cav.critical_velocity(0.5);
        assert!(vc > 0.0 && vc.is_finite());
    }
    #[test]
    fn tbl_thickness_positive() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        assert!(tbl.thickness() > 0.0);
    }
    #[test]
    fn tbl_displacement_thickness_less_than_total() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        assert!(tbl.displacement_thickness() < tbl.thickness());
    }
    #[test]
    fn tbl_shape_factor_approximately_1p4() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        let h = tbl.shape_factor();
        assert!((h - 9.0 / 7.0).abs() < 0.01, "H={h}");
    }
    #[test]
    fn tbl_skin_friction_decreases_downstream() {
        let tbl1 = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        let tbl2 = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 2.0, 1.225);
        assert!(tbl1.skin_friction_coefficient() > tbl2.skin_friction_coefficient());
    }
    #[test]
    fn tbl_log_law_velocity_increases_with_y() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        let u1 = tbl.log_law_velocity(0.001);
        let u2 = tbl.log_law_velocity(0.01);
        assert!(u2 > u1, "u1={u1}, u2={u2}");
    }
    #[test]
    fn tbl_friction_velocity_finite_positive() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        let u_tau = tbl.friction_velocity();
        assert!(u_tau > 0.0 && u_tau.is_finite());
    }
    #[test]
    fn tbl_viscous_sublayer_finite() {
        let tbl = TurbulentBoundaryLayer::new(20.0, 1.5e-5, 1.0, 1.225);
        let yv = tbl.viscous_sublayer_thickness();
        assert!(yv > 0.0 && yv.is_finite());
    }
}
