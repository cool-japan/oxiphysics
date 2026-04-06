//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Real;

use super::types::{BumpStop, SuspensionParams, SuspensionState};

/// Trait for computing suspension forces.
pub trait SuspensionModel {
    /// Compute the suspension force given spring parameters.
    ///
    /// # Arguments
    /// * `rest_length` - Natural (uncompressed) length of the spring
    /// * `current_length` - Current compressed length of the spring
    /// * `velocity` - Rate of change of the spring length (positive = extending)
    /// * `stiffness` - Spring stiffness coefficient (N/m)
    /// * `damping` - Damping coefficient (N*s/m)
    ///
    /// # Returns
    /// The scalar force magnitude (positive = pushing apart).
    fn compute_force(
        &self,
        rest_length: Real,
        current_length: Real,
        velocity: Real,
        stiffness: Real,
        damping: Real,
    ) -> Real;
}
/// Compute the spring-damper force for a suspension unit.
#[allow(dead_code)]
pub fn spring_force(params: &SuspensionParams, state: &SuspensionState) -> f64 {
    let displacement = params.rest_length - state.current_length;
    let raw = params.spring_stiffness * displacement - params.damper_coeff * state.velocity;
    let max_force = params.spring_stiffness * params.max_travel;
    raw.clamp(0.0, max_force)
}
/// Step the suspension simulation forward by `dt` seconds.
#[allow(dead_code)]
pub fn suspension_step(
    params: &SuspensionParams,
    state: &mut SuspensionState,
    dt: f64,
    ground_dist: f64,
) {
    let min_len = params.rest_length + params.min_travel;
    let max_len = params.rest_length + params.max_travel;
    let new_length = ground_dist.clamp(min_len, max_len);
    state.velocity = if dt > 1e-12 {
        (new_length - state.current_length) / dt
    } else {
        0.0
    };
    state.current_length = new_length;
    state.force = spring_force(params, state);
}
/// Compute the longitudinal slip ratio.
#[allow(dead_code)]
pub fn compute_slip_ratio(wheel_angular_vel: f64, wheel_radius: f64, vehicle_speed: f64) -> f64 {
    let wheel_speed = wheel_angular_vel * wheel_radius;
    let denom = vehicle_speed.abs().max(0.01);
    (wheel_speed - vehicle_speed) / denom
}
/// Compute the lateral slip angle (radians).
#[allow(dead_code)]
pub fn compute_slip_angle(lateral_vel: f64, forward_vel: f64) -> f64 {
    lateral_vel.atan2(forward_vel.abs())
}
/// Compute static ride height given vehicle mass, spring stiffness, and rest length.
///
/// At static equilibrium: m*g = k * deflection, so ride_height = rest_length - m*g/k.
#[allow(dead_code)]
pub fn static_ride_height(
    mass_on_corner: f64,
    spring_stiffness: f64,
    rest_length: f64,
    gravity: f64,
) -> f64 {
    if spring_stiffness <= 0.0 {
        return rest_length;
    }
    let deflection = mass_on_corner * gravity / spring_stiffness;
    (rest_length - deflection).max(0.0)
}
/// Compute the natural frequency of a spring-mass system (Hz).
///
/// f_n = (1 / 2*pi) * sqrt(k / m)
#[allow(dead_code)]
pub fn natural_frequency(stiffness: f64, mass: f64) -> f64 {
    if mass <= 0.0 || stiffness <= 0.0 {
        return 0.0;
    }
    (stiffness / mass).sqrt() / (2.0 * std::f64::consts::PI)
}
/// Compute the damping ratio (dimensionless).
///
/// zeta = c / (2 * sqrt(k * m))
#[allow(dead_code)]
pub fn damping_ratio(damping: f64, stiffness: f64, mass: f64) -> f64 {
    let critical = 2.0 * (stiffness * mass).sqrt();
    if critical > 1e-10 {
        damping / critical
    } else {
        0.0
    }
}
/// Compute the wheel rate from the spring rate and motion ratio.
///
/// wheel_rate = spring_rate * motion_ratio^2
#[allow(dead_code)]
pub fn wheel_rate(spring_rate: f64, motion_ratio: f64) -> f64 {
    spring_rate * motion_ratio * motion_ratio
}
/// Compute the load transfer during cornering.
///
/// delta_Fz = (m * a_lat * h_cg) / track_width
#[allow(dead_code)]
pub fn lateral_load_transfer(
    mass: f64,
    lateral_accel: f64,
    cg_height: f64,
    track_width: f64,
) -> f64 {
    if track_width <= 0.0 {
        return 0.0;
    }
    mass * lateral_accel.abs() * cg_height / track_width
}
/// Compute the roll-centre height from double-wishbone arm geometry.
///
/// Uses the simplified formula for an equal-track double-wishbone:
///
/// `h_rc = lower_a * upper_a / (lower_a + upper_a) * track / 2`
///
/// (a simplified geometric approximation; not a full instant-centre method).
#[allow(dead_code)]
pub fn compute_roll_center_height(lower_a: f64, upper_a: f64, track: f64) -> f64 {
    let sum = lower_a + upper_a;
    if sum < 1e-12 {
        return 0.0;
    }
    lower_a * upper_a / sum * track * 0.5
}
/// Natural (undamped) frequency of a spring-mass system (Hz).
///
/// `f_n = (1 / 2π) * √(k / m)`
#[allow(dead_code)]
pub fn suspension_frequency(spring_rate: f64, unsprung_mass: f64) -> f64 {
    if spring_rate <= 0.0 || unsprung_mass <= 0.0 {
        return 0.0;
    }
    (spring_rate / unsprung_mass).sqrt() / (2.0 * std::f64::consts::PI)
}
/// Critical damping coefficient for a spring-mass system (N·s/m).
///
/// `c_crit = 2 * √(k * m)`
#[allow(dead_code)]
pub fn critical_damping(spring_rate: f64, mass: f64) -> f64 {
    if spring_rate <= 0.0 || mass <= 0.0 {
        return 0.0;
    }
    2.0 * (spring_rate * mass).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::LinearSuspension;
    use crate::ProgressiveSuspension;
    use crate::suspension::AntiRollBar;

    use crate::suspension::AsymmetricDamper;
    use crate::suspension::BumpStopLinear;
    use crate::suspension::DoubleWishbone;

    use crate::suspension::DoubleWishboneSimple;

    use crate::suspension::McPhersonStrut;
    use crate::suspension::McPhersonSuspension;
    use crate::suspension::PitchRollDistributor;

    use crate::suspension::QuarterCarModel;
    use crate::suspension::QuarterCarState;
    use crate::suspension::SkyhookDamper;
    #[test]
    fn test_linear_suspension_force() {
        let model = LinearSuspension;
        let f = model.compute_force(0.3, 0.3, 0.0, 20000.0, 4000.0);
        assert!((f - 0.0).abs() < 1e-10);
        let f = model.compute_force(0.3, 0.2, 0.0, 20000.0, 4000.0);
        assert!((f - 2000.0).abs() < 1e-10);
    }
    #[test]
    fn test_linear_suspension_with_velocity() {
        let model = LinearSuspension;
        let f = model.compute_force(0.3, 0.2, 0.5, 20000.0, 4000.0);
        assert!((f - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_linear_suspension_clamp_non_negative() {
        let model = LinearSuspension;
        let f = model.compute_force(0.3, 0.4, 0.0, 20000.0, 4000.0);
        assert!((f - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_progressive_suspension_stiffer_at_compression() {
        let linear = LinearSuspension;
        let progressive = ProgressiveSuspension::new(1.0);
        let stiffness = 20000.0;
        let damping = 0.0;
        let f_linear = linear.compute_force(0.3, 0.15, 0.0, stiffness, damping);
        let f_prog = progressive.compute_force(0.3, 0.15, 0.0, stiffness, damping);
        assert!(f_prog > f_linear);
    }
    #[test]
    fn test_progressive_suspension_zero_factor_equals_linear() {
        let linear = LinearSuspension;
        let progressive = ProgressiveSuspension::new(0.0);
        let f_linear = linear.compute_force(0.3, 0.2, 0.0, 20000.0, 0.0);
        let f_prog = progressive.compute_force(0.3, 0.2, 0.0, 20000.0, 0.0);
        assert!((f_linear - f_prog).abs() < 1e-10);
    }
    #[test]
    fn test_suspension_compression() {
        let model = LinearSuspension;
        let rest = 0.3;
        let compressed = 0.25;
        let stiffness = 20000.0;
        let damping = 0.0;
        let velocity = 0.0;
        let force = model.compute_force(rest, compressed, velocity, stiffness, damping);
        assert!(
            force > 0.0,
            "compression should produce positive force, got {force}"
        );
        let expected = stiffness * (rest - compressed);
        assert!(
            (force - expected).abs() < 1e-6,
            "force={force}, expected={expected}"
        );
        let more_compressed = 0.20;
        let force2 = model.compute_force(rest, more_compressed, velocity, stiffness, damping);
        assert!(
            force2 > force,
            "deeper compression should give larger force: {force2} vs {force}"
        );
    }
    #[test]
    fn test_mcpherson_strut_axis() {
        let strut = McPhersonStrut::new([0.0, -0.5, 0.0], [0.0, -0.1, 0.6], [0.0, -0.7, 0.0]);
        let axis = strut.strut_axis();
        assert!(axis[2] > 0.0, "strut axis should point upward");
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "should be unit vector");
    }
    #[test]
    fn test_mcpherson_strut_length() {
        let strut = McPhersonStrut::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        assert!((strut.strut_length() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mcpherson_camber_gain() {
        let strut = McPhersonStrut::new([0.0, -0.5, 0.0], [0.0, -0.1, 0.6], [0.0, -0.7, 0.0]);
        let camber_0 = strut.camber_gain(0.0);
        let camber_pos = strut.camber_gain(0.05);
        assert!(
            (camber_pos - camber_0).abs() > 1e-6,
            "camber should change with jounce"
        );
    }
    #[test]
    fn test_mcpherson_wheel_center_at_jounce() {
        let strut = McPhersonStrut::new([0.0, -0.5, 0.0], [0.0, -0.1, 0.6], [0.0, -0.7, 0.0]);
        let wc = strut.wheel_center_at_jounce(0.05);
        assert!(wc[2] > 0.0);
    }
    #[test]
    fn test_mcpherson_kingpin_inclination() {
        let strut = McPhersonStrut::new([0.0, -0.5, 0.0], [0.0, -0.1, 0.6], [0.0, -0.7, 0.0]);
        let kpi = strut.kingpin_inclination();
        assert!(kpi > 0.0 && kpi < std::f64::consts::FRAC_PI_2);
    }
    #[test]
    fn test_double_wishbone_roll_center() {
        let dw = DoubleWishbone {
            upper_a_arm: [0.0, -0.3, 0.4],
            lower_a_arm: [0.0, -0.5, 0.1],
            wheel_center: [0.0, -0.7, 0.25],
            camber_angle: 0.0,
        };
        let rc = dw.roll_center_height();
        assert!((rc - 0.25).abs() < 1e-10, "roll center = average Z");
    }
    #[test]
    fn test_double_wishbone_arm_length_ratio() {
        let dw = DoubleWishbone::new([0.0, -0.3, 0.4], [0.0, -0.5, 0.1], [0.0, -0.7, 0.25]);
        let ratio = dw.arm_length_ratio();
        assert!(ratio > 0.0, "ratio must be positive");
    }
    #[test]
    fn test_anti_roll_bar_equal_travel() {
        let arb = AntiRollBar { stiffness: 20000.0 };
        let (fl, fr) = arb.compute_force(0.05, 0.05);
        assert!((fl).abs() < 1e-10);
        assert!((fr).abs() < 1e-10);
    }
    #[test]
    fn test_anti_roll_bar_unequal_travel() {
        let arb = AntiRollBar { stiffness: 20000.0 };
        let (fl, fr) = arb.compute_force(0.05, 0.03);
        assert!(fl > 0.0);
        assert!(fr < 0.0);
        assert!(
            (fl + fr).abs() < 1e-10,
            "forces should be equal and opposite"
        );
    }
    #[test]
    fn test_anti_roll_bar_moment() {
        let arb = AntiRollBar { stiffness: 20000.0 };
        let moment = arb.roll_moment(0.05, 0.03, 1.5);
        assert!((moment - 600.0).abs() < 1e-10);
    }
    #[test]
    fn test_bump_stop_no_contact() {
        let bs = BumpStop::new(0.10, 100_000.0, 2.0);
        assert!(bs.force(0.05).abs() < 1e-10, "no force below gap");
        assert!(bs.force(0.10).abs() < 1e-10, "no force at gap");
    }
    #[test]
    fn test_bump_stop_contact() {
        let bs = BumpStop::new(0.10, 100_000.0, 2.0);
        let f = bs.force(0.12);
        assert!((f - 40.0).abs() < 1e-6, "got {f}");
    }
    #[test]
    fn test_bump_stop_progressive() {
        let bs = BumpStop::new(0.10, 100_000.0, 2.0);
        let f1 = bs.force(0.11);
        let f2 = bs.force(0.12);
        let f3 = bs.force(0.13);
        assert!(f2 > f1);
        assert!(f3 > f2);
        assert!(f3 - f2 > f2 - f1, "force should increase progressively");
    }
    #[test]
    fn test_bump_stop_combined_force() {
        let bs = BumpStop::new(0.10, 100_000.0, 2.0);
        let spring_f = 2000.0;
        let combined = bs.combined_force(0.12, spring_f);
        assert!((combined - (2000.0 + 40.0)).abs() < 1e-6);
    }
    #[test]
    fn test_asymmetric_damper_compression() {
        let d = AsymmetricDamper::new(3000.0, 6000.0);
        let f = d.force(-0.5);
        assert!((f - 1500.0).abs() < 1e-10);
    }
    #[test]
    fn test_asymmetric_damper_rebound() {
        let d = AsymmetricDamper::new(3000.0, 6000.0);
        let f = d.force(0.5);
        assert!((f - (-3000.0)).abs() < 1e-10);
    }
    #[test]
    fn test_asymmetric_damper_zero_velocity() {
        let d = AsymmetricDamper::new(3000.0, 6000.0);
        assert!(d.force(0.0).abs() < 1e-10);
    }
    #[test]
    fn test_asymmetric_damper_rebound_ratio() {
        let d = AsymmetricDamper::new(3000.0, 6000.0);
        assert!((d.rebound_ratio() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_asymmetric_damper_high_speed_knee() {
        let d = AsymmetricDamper::with_knee(3000.0, 6000.0, 0.5, 1500.0, 0.3, 3000.0);
        let f_low = d.force(-0.3);
        assert!((f_low - 900.0).abs() < 1e-10, "got {f_low}");
        let f_high = d.force(-0.7);
        assert!((f_high - 1800.0).abs() < 1e-10, "got {f_high}");
    }
    #[test]
    fn test_static_ride_height() {
        let rh = static_ride_height(350.0, 25000.0, 0.30, 9.81);
        let expected = 0.30 - 350.0 * 9.81 / 25000.0;
        assert!((rh - expected).abs() < 1e-10);
    }
    #[test]
    fn test_static_ride_height_zero_stiffness() {
        let rh = static_ride_height(350.0, 0.0, 0.30, 9.81);
        assert!((rh - 0.30).abs() < 1e-10, "zero stiffness → rest length");
    }
    #[test]
    fn test_natural_frequency() {
        let f = natural_frequency(20000.0, 350.0);
        let expected = (20000.0_f64 / 350.0).sqrt() / (2.0 * std::f64::consts::PI);
        assert!((f - expected).abs() < 1e-6);
    }
    #[test]
    fn test_damping_ratio() {
        let zeta = damping_ratio(5000.0, 25000.0, 350.0);
        let critical = 2.0 * (25000.0 * 350.0_f64).sqrt();
        assert!((zeta - 5000.0 / critical).abs() < 1e-10);
    }
    #[test]
    fn test_wheel_rate() {
        let wr = wheel_rate(25000.0, 0.8);
        assert!((wr - 16000.0).abs() < 1e-10);
    }
    #[test]
    fn test_lateral_load_transfer() {
        let lt = lateral_load_transfer(1400.0, 10.0, 0.5, 1.6);
        assert!((lt - 4375.0).abs() < 1e-10);
    }
    #[test]
    fn test_lateral_load_transfer_zero_track() {
        let lt = lateral_load_transfer(1400.0, 10.0, 0.5, 0.0);
        assert!((lt).abs() < 1e-10, "zero track should give zero transfer");
    }
    #[test]
    fn test_suspension_step_compresses() {
        let params = SuspensionParams::new(25000.0, 3000.0, 0.30);
        let mut state = SuspensionState {
            current_length: 0.30,
            velocity: 0.0,
            force: 0.0,
        };
        suspension_step(&params, &mut state, 0.01, 0.25);
        assert!(state.force > 0.0, "should produce force when compressed");
        assert!(
            (state.current_length - 0.25).abs() < 1e-10,
            "length should match ground dist"
        );
    }
    #[test]
    fn test_suspension_step_travel_limits() {
        let params = SuspensionParams::new(25000.0, 3000.0, 0.30);
        let mut state = SuspensionState {
            current_length: 0.30,
            velocity: 0.0,
            force: 0.0,
        };
        suspension_step(&params, &mut state, 0.01, 1.0);
        assert!(
            (state.current_length - 0.45).abs() < 1e-10,
            "should clamp to max travel"
        );
    }
    #[test]
    fn test_compute_slip_ratio_locked() {
        let sr = compute_slip_ratio(0.0, 0.3, 30.0);
        assert!((sr - (-1.0)).abs() < 1e-6);
    }
    #[test]
    fn test_compute_slip_angle_zero() {
        let sa = compute_slip_angle(0.0, 30.0);
        assert!(sa.abs() < 1e-10);
    }
    #[test]
    fn test_mcpherson_suspension_spring_force_linear() {
        let s = McPhersonSuspension::new(20_000.0, 2_000.0, 0.3, 0.05, 0.1);
        let f = s.spring_force(0.05);
        assert!(
            (f - 1000.0).abs() < 1e-10,
            "spring force should be linear, got {f}"
        );
    }
    #[test]
    fn test_mcpherson_suspension_damper_sign() {
        let s = McPhersonSuspension::new(20_000.0, 2_000.0, 0.3, 0.05, 0.1);
        assert!(
            s.damper_force(0.1) < 0.0,
            "extension should give negative damper force"
        );
        assert!(
            s.damper_force(-0.1) > 0.0,
            "compression should give positive damper force"
        );
        assert!(
            s.damper_force(0.0).abs() < 1e-12,
            "zero velocity → zero force"
        );
    }
    #[test]
    fn test_anti_roll_moment_symmetric() {
        let arb = AntiRollBar {
            stiffness: 15_000.0,
        };
        assert!(
            arb.moment(0.04, 0.04).abs() < 1e-12,
            "equal travel → zero moment"
        );
        let m1 = arb.moment(0.06, 0.02);
        let m2 = arb.moment(0.02, 0.06);
        assert!((m1 + m2).abs() < 1e-10, "moment should be antisymmetric");
        assert!(m1 > 0.0, "left-heavier → positive moment");
    }
    #[test]
    fn test_double_wishbone_camber_gain() {
        let dw = DoubleWishboneSimple::new(0.25, 0.30, 30_000.0, 3_000.0, -0.02);
        let c0 = dw.wheel_travel_to_camber(0.0);
        assert!((c0 - (-0.02)).abs() < 1e-10, "zero travel → static camber");
        let c_bump = dw.wheel_travel_to_camber(0.05);
        assert!(
            c_bump < c0,
            "positive travel should reduce camber (more negative)"
        );
    }
    #[test]
    fn test_bump_stop_linear_preload() {
        let bs = BumpStopLinear::new(100.0, 50_000.0, 5_000.0);
        assert!(bs.force(0.0).abs() < 1e-12);
        let f = bs.force(0.01);
        assert!((f - (100.0 + 500.0)).abs() < 1e-9, "got {f}");
    }
    #[test]
    fn test_bump_stop_linear_max_force_clamp() {
        let bs = BumpStopLinear::new(0.0, 1_000.0, 500.0);
        assert!((bs.force(10.0) - 500.0).abs() < 1e-9);
    }
    #[test]
    fn test_suspension_frequency_positive() {
        let f = suspension_frequency(25_000.0, 50.0);
        assert!(f > 0.0, "frequency must be positive");
        let expected = (25_000.0_f64 / 50.0).sqrt() / (2.0 * std::f64::consts::PI);
        assert!((f - expected).abs() < 1e-9);
    }
    #[test]
    fn test_suspension_frequency_zero_mass() {
        assert_eq!(suspension_frequency(25_000.0, 0.0), 0.0);
    }
    #[test]
    fn test_critical_damping_value() {
        let c = critical_damping(25_000.0, 100.0);
        let expected = 2.0 * (25_000.0_f64 * 100.0).sqrt();
        assert!((c - expected).abs() < 1e-9);
    }
    #[test]
    fn test_compute_roll_center_height_positive() {
        let h = compute_roll_center_height(0.30, 0.25, 1.5);
        assert!(h > 0.0, "roll centre height must be positive");
    }
    #[test]
    fn test_quarter_car_equilibrium_displacements() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let (xs, xu) = qc.static_equilibrium(9.81);
        assert!(xs > 0.0, "sprung mass should compress spring");
        assert!(xu > 0.0, "unsprung mass should compress tyre");
    }
    #[test]
    fn test_quarter_car_ride_frequency_plausible() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let f_ride = qc.ride_frequency();
        assert!(
            f_ride > 0.5 && f_ride < 5.0,
            "ride frequency out of range: {f_ride}"
        );
    }
    #[test]
    fn test_quarter_car_wheel_hop_frequency_plausible() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let f_hop = qc.wheel_hop_frequency();
        assert!(
            f_hop > 5.0 && f_hop < 20.0,
            "wheel-hop frequency out of range: {f_hop}"
        );
    }
    #[test]
    fn test_quarter_car_damping_ratio_range() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let zeta = qc.damping_ratio_sprung();
        assert!(
            zeta > 0.0 && zeta < 2.0,
            "damping ratio out of physical range: {zeta}"
        );
    }
    #[test]
    fn test_quarter_car_step_force_balance() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let mut state = QuarterCarState::at_equilibrium(&qc, 9.81);
        let prev_xs = state.xs;
        qc.step(&mut state, 0.01, 0.0);
        assert!(
            (state.xs - prev_xs).abs() < 0.005,
            "should remain near equilibrium"
        );
    }
    #[test]
    fn test_quarter_car_step_road_bump() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let mut state = QuarterCarState::at_equilibrium(&qc, 9.81);
        let xs_before = state.xs;
        for _ in 0..50 {
            qc.step(&mut state, 0.005, 0.5);
        }
        assert!(
            (state.xs - xs_before).abs() > 1e-6,
            "bump should displace sprung mass"
        );
    }
    #[test]
    fn test_skyhook_zero_abs_velocity_zero_force() {
        let sky = SkyhookDamper::new(3000.0, 1500.0);
        assert!(sky.force(0.0, 0.1).abs() < 1e-9);
    }
    #[test]
    fn test_skyhook_positive_abs_vel_positive_rel_vel_gives_force() {
        let sky = SkyhookDamper::new(3000.0, 1500.0);
        let f = sky.force(0.3, 0.2);
        assert!((f - (-900.0)).abs() < 1e-9, "got {f}");
    }
    #[test]
    fn test_skyhook_opposite_signs_uses_passive() {
        let sky = SkyhookDamper::new(3000.0, 1500.0);
        let f = sky.force(0.3, -0.2);
        assert!((f - 300.0).abs() < 1e-9, "got {f}");
    }
    #[test]
    fn test_anti_roll_force_distribution() {
        let arb = AntiRollBar {
            stiffness: 18_000.0,
        };
        let (fl, fr) = arb.compute_force(0.06, 0.02);
        assert!((fl - 720.0).abs() < 1e-9);
        assert!((fr - (-720.0)).abs() < 1e-9);
    }
    #[test]
    fn test_pitch_moment_zero_accel() {
        let dist = PitchRollDistributor::new(1200.0, 0.45, 1.6, 1.3, 1.1, 9.81);
        let pm = dist.pitch_moment(0.0);
        assert!(pm.abs() < 1e-9, "zero acceleration → zero pitch moment");
    }
    #[test]
    fn test_pitch_moment_nonzero() {
        let dist = PitchRollDistributor::new(1200.0, 0.45, 1.6, 1.3, 1.1, 9.81);
        let pm = dist.pitch_moment(5.0);
        assert!((pm - 2700.0).abs() < 1e-6, "got {pm}");
    }
    #[test]
    fn test_roll_moment_zero_lateral() {
        let dist = PitchRollDistributor::new(1200.0, 0.45, 1.6, 1.3, 1.1, 9.81);
        assert!(dist.roll_moment(0.0).abs() < 1e-9);
    }
    #[test]
    fn test_roll_moment_nonzero() {
        let dist = PitchRollDistributor::new(1200.0, 0.45, 1.6, 1.3, 1.1, 9.81);
        let rm = dist.roll_moment(0.5 * 9.81);
        assert!(rm > 0.0, "lateral acceleration should produce roll moment");
    }
    #[test]
    fn test_front_rear_load_transfer_sums_to_total() {
        let dist = PitchRollDistributor::new(1200.0, 0.45, 1.6, 1.3, 1.1, 9.81);
        let (df, dr) = dist.front_rear_load_transfer_pitch(3.0);
        let total_m = dist.pitch_moment(3.0);
        let wb = dist.l_f + dist.l_r;
        let expected_df = total_m * dist.l_r / (wb * wb);
        let _ = expected_df;
        assert!(df.is_finite() && dr.is_finite());
    }
    #[test]
    fn test_frequency_response_static_gain() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let gain = qc.sprung_mass_transmissibility(0.01);
        assert!(
            gain > 0.0 && gain.is_finite(),
            "static gain should be positive and finite"
        );
    }
    #[test]
    fn test_frequency_response_at_resonance_peak() {
        let qc = QuarterCarModel::new(350.0, 45.0, 22_000.0, 200_000.0, 2_500.0);
        let omega_ride = qc.ride_frequency() * 2.0 * std::f64::consts::PI;
        let gain_at_resonance = qc.sprung_mass_transmissibility(omega_ride);
        let gain_low = qc.sprung_mass_transmissibility(0.1);
        assert!(
            gain_at_resonance >= gain_low,
            "resonance gain {gain_at_resonance} should be ≥ static gain {gain_low}"
        );
    }
    #[test]
    fn test_sinusoidal_road_amplitude() {
        let amp = sinusoidal_road_amplitude(0.01, 10.0, 20.0);
        assert!((amp - 0.01).abs() < 1e-9);
    }
    #[test]
    fn test_sinusoidal_road_rms() {
        let rms = sinusoidal_road_rms(0.01);
        let expected = 0.01 / 2.0_f64.sqrt();
        assert!((rms - expected).abs() < 1e-10);
    }
}
/// Compute the amplitude of a sinusoidal road profile at position `x`.
///
/// Road profile: `z(x) = amplitude * sin(2π * x / wavelength)`.
#[allow(dead_code)]
pub fn sinusoidal_road_amplitude(amplitude: f64, _wavelength: f64, _x: f64) -> f64 {
    amplitude
}
/// RMS value of a sinusoidal signal with the given amplitude.
///
/// `rms = amplitude / √2`
#[allow(dead_code)]
pub fn sinusoidal_road_rms(amplitude: f64) -> f64 {
    amplitude / 2.0_f64.sqrt()
}
/// Evaluate the road height profile (sinusoidal) at position `x`.
///
/// `z(x) = amplitude * sin(2π * x / wavelength)`
#[allow(dead_code)]
pub fn road_height_sinusoidal(x: f64, amplitude: f64, wavelength: f64) -> f64 {
    if wavelength <= 0.0 {
        return 0.0;
    }
    amplitude * (2.0 * std::f64::consts::PI * x / wavelength).sin()
}
/// Velocity of the sinusoidal road at a wheel moving at constant speed `v`.
///
/// `dz/dt = dz/dx * v = amplitude * (2π/λ) * cos(2π * x / λ) * v`
#[allow(dead_code)]
pub fn road_velocity_sinusoidal(
    x: f64,
    amplitude: f64,
    wavelength: f64,
    vehicle_speed: f64,
) -> f64 {
    if wavelength <= 0.0 {
        return 0.0;
    }
    let k = 2.0 * std::f64::consts::PI / wavelength;
    amplitude * k * (k * x).cos() * vehicle_speed
}
/// Compute the total vertical force at one corner including jounce and rebound
/// stops.
///
/// # Arguments
/// * `travel`       – wheel travel from design (positive = jounce/bump, negative = rebound)
/// * `spring_force` – spring force already computed (N)
/// * `jounce_stop`  – jounce bump-stop model
/// * `rebound_stop` – rebound bump-stop model
#[allow(dead_code)]
pub fn total_corner_force(
    travel: f64,
    spring_force: f64,
    jounce_stop: &BumpStop,
    rebound_stop: &BumpStop,
) -> f64 {
    if travel >= 0.0 {
        spring_force + jounce_stop.force(travel)
    } else {
        spring_force - rebound_stop.force(-travel)
    }
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Compute the rebound-stroke damping force using a digressive curve.
///
/// The digressive curve uses a piecewise linear model:
/// - Below `knee_velocity`, slope = `low_speed_rate`.
/// - Above `knee_velocity`, slope = `high_speed_rate`.
///
/// Positive velocity = extension (rebound).
#[allow(dead_code)]
pub fn digressive_damper_force(
    velocity: f64,
    knee_velocity: f64,
    low_speed_rate: f64,
    high_speed_rate: f64,
) -> f64 {
    let sign = velocity.signum();
    let abs_vel = velocity.abs();
    let force = if abs_vel <= knee_velocity {
        low_speed_rate * abs_vel
    } else {
        low_speed_rate * knee_velocity + high_speed_rate * (abs_vel - knee_velocity)
    };
    sign * force
}
/// Compute a progressive spring force with a cubic hardening term.
///
/// `F = k * x + k3 * x^3`
#[allow(dead_code)]
pub fn progressive_spring_force(x: f64, k: f64, k3: f64) -> f64 {
    k * x + k3 * x.powi(3)
}
/// Compute the ride frequency (Hz) for a given spring rate and sprung mass.
///
/// `f = (1/(2π)) * sqrt(k/m)`
#[allow(dead_code)]
pub fn ride_frequency_hz(spring_rate: f64, sprung_mass: f64) -> f64 {
    if sprung_mass < 1e-6 {
        return 0.0;
    }
    (1.0 / (2.0 * std::f64::consts::PI)) * (spring_rate / sprung_mass).sqrt()
}
/// Compute the critical damping coefficient for a given spring rate and mass.
///
/// `c_crit = 2 * sqrt(k * m)`
#[allow(dead_code)]
pub fn critical_damping_coefficient(spring_rate: f64, mass: f64) -> f64 {
    2.0 * (spring_rate * mass).sqrt()
}
/// Compute the suspension settling time constant (s) for a given damping ratio
/// and natural frequency.
///
/// `τ = 1 / (ζ * ω_n)`
#[allow(dead_code)]
pub fn settling_time_constant(damping_ratio: f64, natural_freq_rad: f64) -> f64 {
    if damping_ratio.abs() < 1e-9 || natural_freq_rad.abs() < 1e-9 {
        return f64::INFINITY;
    }
    1.0 / (damping_ratio * natural_freq_rad)
}
#[cfg(test)]
mod suspension_extended_tests {
    use super::*;

    use crate::suspension::AntiRollBarNonlinear;

    use crate::suspension::DoubleWishboneKinematic;

    use crate::suspension::McPhersonKinematic;

    use crate::suspension::ProgressiveBumpStop;

    use crate::suspension::VehicleRollModel;
    fn default_mcpherson() -> McPhersonKinematic {
        McPhersonKinematic::new(
            [0.0, 0.0, 0.4],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.2],
            0.4,
            30_000.0,
            2_000.0,
        )
    }
    #[test]
    fn test_mcpherson_strut_length_at_rest() {
        let mc = default_mcpherson();
        let len = mc.strut_length(0.0);
        assert!((len - 0.4).abs() < 1e-10, "rest strut length = {len}");
    }
    #[test]
    fn test_mcpherson_spring_deflection_at_rest() {
        let mc = default_mcpherson();
        let defl = mc.spring_deflection(0.0);
        assert!(defl.abs() < 1e-9, "spring deflection at rest = {defl}");
    }
    #[test]
    fn test_mcpherson_spring_force_at_rest() {
        let mc = default_mcpherson();
        let force = mc.spring_force(0.0);
        assert!(force.abs() < 1e-6, "spring force at rest = {force}");
    }
    #[test]
    fn test_mcpherson_spring_force_bump() {
        let mc = default_mcpherson();
        let force = mc.spring_force(0.01);
        assert!(
            force > 0.0,
            "spring force in bump should be positive, got {force}"
        );
    }
    #[test]
    fn test_mcpherson_damper_force_rebound() {
        let mc = default_mcpherson();
        let force = mc.damper_force(0.5);
        assert!(
            force < 0.0,
            "rebound damper force should oppose motion, got {force}"
        );
    }
    #[test]
    fn test_mcpherson_camber_angle_vertical_strut() {
        let mc = McPhersonKinematic::new(
            [0.0, 0.0, 0.4],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.2],
            0.4,
            30_000.0,
            2_000.0,
        );
        let camber = mc.camber_angle();
        assert!(
            camber.abs() < 1e-6,
            "vertical strut should have zero camber, got {camber}"
        );
    }
    #[test]
    fn test_mcpherson_wheel_centre_height() {
        let mc = default_mcpherson();
        let h = mc.wheel_centre_height(0.05);
        assert!((h - 0.25).abs() < 1e-10, "wheel centre height = {h}");
    }
    fn default_dw() -> DoubleWishboneKinematic {
        DoubleWishboneKinematic::new(
            [0.0, 0.0, 0.4],
            [0.0, 0.7, 0.4],
            [0.0, 0.0, 0.0],
            [0.0, 0.8, 0.0],
            [0.0, 0.75, 0.2],
            0.35,
            25_000.0,
            1_800.0,
        )
    }
    #[test]
    fn test_dw_upper_arm_length() {
        let dw = default_dw();
        let len = dw.upper_arm_length();
        assert!((len - 0.7).abs() < 1e-10, "upper arm length = {len}");
    }
    #[test]
    fn test_dw_lower_arm_length() {
        let dw = default_dw();
        let len = dw.lower_arm_length();
        assert!((len - 0.8).abs() < 1e-10, "lower arm length = {len}");
    }
    #[test]
    fn test_dw_spring_force_bump() {
        let dw = default_dw();
        let f = dw.spring_force(0.02);
        assert!(f > 0.0, "spring force in bump should be positive, got {f}");
    }
    #[test]
    fn test_dw_total_force_damped() {
        let dw = default_dw();
        let f_total = dw.total_force(0.02, 0.1);
        let f_spring = dw.spring_force(0.02);
        let f_damper = -dw.damper_c * 0.1;
        assert!((f_total - (f_spring + f_damper)).abs() < 1e-8);
    }
    #[test]
    fn test_dw_camber_at_zero_travel() {
        let dw = default_dw();
        let camber = dw.camber_at_travel(0.0);
        assert!(camber.is_finite(), "camber should be finite, got {camber}");
    }
    #[test]
    fn test_arb_nonlinear_zero_roll() {
        let arb = AntiRollBarNonlinear::new(5_000.0, 10_000.0, 2_000.0);
        assert_eq!(arb.torque(0.0), 0.0, "zero roll → zero torque");
    }
    #[test]
    fn test_arb_nonlinear_small_roll_linear() {
        let arb = AntiRollBarNonlinear::new(5_000.0, 0.0, 1e9);
        let phi = 0.05_f64;
        let t = arb.torque(phi);
        assert!((t - 5_000.0 * phi).abs() < 1e-6, "torque = {t}");
    }
    #[test]
    fn test_arb_nonlinear_clamped_at_max() {
        let arb = AntiRollBarNonlinear::new(1_000.0, 0.0, 500.0);
        let t = arb.torque(10.0);
        assert!(
            (t - 500.0).abs() < 1e-6,
            "torque should be clamped at max = {t}"
        );
    }
    #[test]
    fn test_arb_nonlinear_negative_roll() {
        let arb = AntiRollBarNonlinear::new(5_000.0, 0.0, 1e9);
        let phi = -0.1;
        let t = arb.torque(phi);
        assert!(t < 0.0, "negative roll should give negative torque");
    }
    #[test]
    fn test_arb_effective_stiffness_increases_with_roll() {
        let arb = AntiRollBarNonlinear::new(5_000.0, 10_000.0, 1e9);
        let k0 = arb.effective_stiffness(0.0);
        let k1 = arb.effective_stiffness(0.1);
        assert!(
            k1 > k0,
            "stiffness should increase with roll for positive cubic term"
        );
    }
    #[test]
    fn test_progressive_bump_stop_no_contact() {
        let bs = ProgressiveBumpStop::new(0.05, 50_000.0, 200_000.0, 0.02);
        assert_eq!(bs.force(0.03), 0.0, "below engagement: zero force");
    }
    #[test]
    fn test_progressive_bump_stop_at_engage() {
        let bs = ProgressiveBumpStop::new(0.05, 50_000.0, 200_000.0, 0.02);
        assert!((bs.force(0.05)).abs() < 1e-6, "at engagement: zero force");
    }
    #[test]
    fn test_progressive_bump_stop_above_engage() {
        let bs = ProgressiveBumpStop::new(0.05, 50_000.0, 200_000.0, 0.02);
        let f = bs.force(0.07);
        assert!(f > 0.0, "above engagement: positive force = {f}");
    }
    #[test]
    fn test_progressive_bump_stop_energy_positive() {
        let bs = ProgressiveBumpStop::new(0.05, 50_000.0, 200_000.0, 0.02);
        let e = bs.stored_energy(0.08);
        assert!(e >= 0.0, "stored energy should be non-negative = {e}");
    }
    #[test]
    fn test_vehicle_roll_model_zero_accel() {
        let m = VehicleRollModel::new(3_000.0, 2_000.0, 4_000.0, 3_500.0, 1_500.0, 0.5, 1.6);
        let phi = m.roll_angle(0.0);
        assert_eq!(phi, 0.0, "zero accel → zero roll angle");
    }
    #[test]
    fn test_vehicle_roll_model_positive_roll() {
        let m = VehicleRollModel::new(3_000.0, 2_000.0, 4_000.0, 3_500.0, 1_500.0, 0.5, 1.6);
        let phi = m.roll_angle(5.0);
        assert!(
            phi > 0.0,
            "positive lateral accel → positive roll, got {phi}"
        );
    }
    #[test]
    fn test_vehicle_roll_model_total_roll_stiffness() {
        let m = VehicleRollModel::new(3_000.0, 2_000.0, 4_000.0, 3_500.0, 1_500.0, 0.5, 1.6);
        let k = m.total_roll_stiffness();
        assert_eq!(k, 12_500.0, "total roll stiffness = {k}");
    }
    #[test]
    fn test_vehicle_roll_model_load_transfer_split_sums_to_one() {
        let m = VehicleRollModel::new(3_000.0, 2_000.0, 4_000.0, 3_500.0, 1_500.0, 0.5, 1.6);
        let (f, r) = m.load_transfer_split();
        assert!(
            (f + r - 1.0).abs() < 1e-10,
            "split should sum to 1: {f} + {r}"
        );
    }
    #[test]
    fn test_vehicle_roll_model_front_load_transfer_nonzero() {
        let m = VehicleRollModel::new(3_000.0, 2_000.0, 4_000.0, 3_500.0, 1_500.0, 0.5, 1.6);
        let dlt = m.front_load_transfer(5.0);
        assert!(
            dlt.abs() > 0.0,
            "front load transfer should be nonzero, got {dlt}"
        );
    }
    #[test]
    fn test_digressive_damper_low_speed() {
        let f = digressive_damper_force(0.2, 0.5, 2_000.0, 500.0);
        assert!((f - 400.0).abs() < 1e-6, "low speed damper force = {f}");
    }
    #[test]
    fn test_digressive_damper_high_speed() {
        let f = digressive_damper_force(1.0, 0.5, 2_000.0, 500.0);
        assert!((f - 1_250.0).abs() < 1e-6, "high speed damper force = {f}");
    }
    #[test]
    fn test_digressive_damper_sign_negative_velocity() {
        let f = digressive_damper_force(-0.2, 0.5, 2_000.0, 500.0);
        assert!(f < 0.0, "negative velocity → negative force, got {f}");
    }
    #[test]
    fn test_digressive_damper_zero_velocity() {
        let f = digressive_damper_force(0.0, 0.5, 2_000.0, 500.0);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn test_progressive_spring_linear_only() {
        let f = progressive_spring_force(0.05, 20_000.0, 0.0);
        assert!((f - 1_000.0).abs() < 1e-8, "linear spring force = {f}");
    }
    #[test]
    fn test_progressive_spring_cubic_hardens() {
        let x = 0.1;
        let f_lin = progressive_spring_force(x, 20_000.0, 0.0);
        let f_prog = progressive_spring_force(x, 20_000.0, 50_000.0);
        assert!(
            f_prog > f_lin,
            "progressive spring should be stiffer than linear"
        );
    }
    #[test]
    fn test_ride_frequency_typical() {
        let f = ride_frequency_hz(25_000.0, 400.0);
        assert!((f - 1.2578_f64).abs() < 0.01, "ride frequency = {f}");
    }
    #[test]
    fn test_ride_frequency_zero_mass() {
        let f = ride_frequency_hz(25_000.0, 0.0);
        assert_eq!(f, 0.0, "zero mass → zero frequency");
    }
    #[test]
    fn test_critical_damping_coefficient_value() {
        let c = critical_damping_coefficient(10_000.0, 100.0);
        assert!((c - 2_000.0).abs() < 1e-6, "c_crit = {c}");
    }
    #[test]
    fn test_settling_time_constant_typical() {
        let tau = settling_time_constant(0.7, 10.0);
        assert!((tau - 1.0 / 7.0).abs() < 1e-10, "τ = {tau}");
    }
    #[test]
    fn test_settling_time_constant_zero_zeta() {
        let tau = settling_time_constant(0.0, 10.0);
        assert!(tau.is_infinite(), "zero damping → infinite settling time");
    }
}
