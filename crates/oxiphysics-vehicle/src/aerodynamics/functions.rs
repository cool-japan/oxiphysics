//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Air density at sea level in kg/m³.
pub(super) const AIR_DENSITY_SEA_LEVEL: f64 = 1.225;
/// Dynamic viscosity of air at 20 °C in Pa·s.
#[cfg(test)]
pub(super) const AIR_VISCOSITY: f64 = 1.81e-5;
/// Linear interpolation helper for 1-D tables.
pub(super) fn interpolate_table(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if xs.is_empty() || ys.is_empty() {
        return 0.0;
    }
    if x <= xs[0] {
        return ys[0];
    }
    let last = xs.len() - 1;
    if x >= xs[last] {
        return ys[last];
    }
    for i in 0..last {
        if x >= xs[i] && x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] + t * (ys[i + 1] - ys[i]);
        }
    }
    ys[last]
}
/// Compute the aerodynamic drag force vector.
///
/// The drag force acts **opposite** to the velocity direction.
///
/// `F_drag = -0.5 · ρ · |v|² · CdA · v̂`
///
/// Returns `[f64; 3]` in N.
pub fn compute_drag(velocity: [f64; 3], air_density: f64, cda: f64) -> [f64; 3] {
    let speed_sq =
        velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    let speed = speed_sq.sqrt();
    if speed < 1e-12 {
        return [0.0; 3];
    }
    let magnitude = 0.5 * air_density * speed_sq * cda;
    [
        -magnitude * velocity[0] / speed,
        -magnitude * velocity[1] / speed,
        -magnitude * velocity[2] / speed,
    ]
}
/// Compute the aerodynamic lift force vector.
///
/// Lift acts along the **+z** axis in body frame.  A negative `cla` (downforce
/// package) produces a force in the −z direction.
///
/// `F_lift = 0.5 · ρ · v_xy² · ClA · ẑ`
///
/// `cop_offset` shifts the magnitude by a fraction (e.g. ground-effect factor).
///
/// Returns `[f64; 3]` in N.
pub fn compute_lift(velocity: [f64; 3], air_density: f64, cla: f64, cop_offset: f64) -> [f64; 3] {
    let vx = velocity[0];
    let vy = velocity[1];
    let speed_sq = vx * vx + vy * vy;
    let magnitude = 0.5 * air_density * speed_sq * cla * (1.0 + cop_offset);
    [0.0, 0.0, magnitude]
}
/// Compute the aerodynamic side-force vector.
///
/// Side force acts along the **+y** axis in body frame.
///
/// `F_side = 0.5 · ρ · v_xz² · CsA · ŷ`
///
/// Returns `[f64; 3]` in N.
pub fn compute_side_force(velocity: [f64; 3], air_density: f64, csa: f64) -> [f64; 3] {
    let vx = velocity[0];
    let vz = velocity[2];
    let speed_sq = vx * vx + vz * vz;
    let magnitude = 0.5 * air_density * speed_sq * csa;
    [0.0, magnitude, 0.0]
}
/// Returns `(index, fraction)` for a value in a sorted slice.
pub(super) fn bracket_index(xs: &[f64], x: f64) -> (usize, f64) {
    let n = xs.len();
    if n == 0 {
        return (0, 0.0);
    }
    if x <= xs[0] {
        return (0, 0.0);
    }
    let last = n - 1;
    if x >= xs[last] {
        return (last.saturating_sub(1), 1.0);
    }
    for i in 0..last {
        if x >= xs[i] && x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return (i, t);
        }
    }
    (last.saturating_sub(1), 1.0)
}
#[inline]
pub(super) fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AeroBalance;
    use crate::AeroBalanceController;
    use crate::AeroBody;

    use crate::AeroForces;
    use crate::AeroMap;

    use crate::AeroWing;
    use crate::AerodynamicBody;
    use crate::CfdSurrogateModel;
    use crate::CompressibilityCorrection;
    use crate::DiffuserModel;
    use crate::DownforceMap;
    use crate::DownforcePackage;
    use crate::DragBreakdown;
    use crate::GroundEffect;
    use crate::GroundEffectAdvanced;
    use crate::SpeedAeroMap;
    use crate::WakeInteraction;
    use crate::WindModel;
    use crate::WingProfile;
    use crate::WingStallModel;
    use std::f64::consts::PI;
    fn test_body() -> AerodynamicBody {
        AerodynamicBody {
            frontal_area: 2.0,
            drag_coefficient: 0.33,
            lift_coefficient: -0.3,
            side_force_coefficient: 0.0,
            length: 4.5,
            wheelbase: 2.7,
        }
    }
    #[test]
    fn test_drag_force_quadratic_with_velocity() {
        let body = test_body();
        let rho = AIR_DENSITY_SEA_LEVEL;
        let f1 = body.drag_force(20.0, rho);
        let f2 = body.drag_force(40.0, rho);
        let ratio = f2 / f1;
        assert!((ratio - 4.0).abs() < 1e-9, "ratio={ratio}");
    }
    #[test]
    fn test_lift_coefficient_symmetric_zero_aoa() {
        let wing = AeroWing {
            span: 1.5,
            chord: 0.3,
            angle_of_attack: 0.0,
            profile: WingProfile::Symmetric {
                max_thickness: 0.12,
            },
        };
        let cl = wing.lift_coefficient();
        assert!(
            cl.abs() < 1e-12,
            "symmetric wing at 0° AOA must have Cl=0, got {cl}"
        );
    }
    #[test]
    fn test_wing_efficiency_ratio_positive_for_nonzero_lift() {
        let wing = AeroWing {
            span: 1.5,
            chord: 0.3,
            angle_of_attack: 5.0,
            profile: WingProfile::Symmetric {
                max_thickness: 0.12,
            },
        };
        let eta = wing.efficiency_ratio();
        assert!(
            eta > 0.0,
            "efficiency ratio should be positive at AOA=5°, got {eta}"
        );
    }
    #[test]
    fn test_ground_effect_enhancement_factor_greater_than_one() {
        let ge = GroundEffect { ride_height: 0.05 };
        let f = ge.enhancement_factor();
        assert!(f > 1.0, "enhancement factor should exceed 1, got {f}");
    }
    #[test]
    fn test_aero_balance_front_fraction_in_range() {
        let body = test_body();
        let balance = AeroBalance::compute(&body, 50.0, AIR_DENSITY_SEA_LEVEL, 1.35, 2.7);
        let ff = balance.front_fraction();
        assert!(
            (0.0..=1.0).contains(&ff),
            "front_fraction must be in [0,1], got {ff}"
        );
    }
    #[test]
    fn test_drag_breakdown_total_equals_sum_of_components() {
        let body = test_body();
        let velocity = 30.0;
        let re = body.reynolds_number(velocity, AIR_VISCOSITY / AIR_DENSITY_SEA_LEVEL);
        let bd = DragBreakdown::compute(&body, velocity, re);
        let sum = bd.pressure_drag + bd.friction_drag + bd.induced_drag + bd.interference_drag;
        let diff = (bd.total() - sum).abs();
        assert!(diff < 1e-10, "total mismatch: {} vs {}", bd.total(), sum);
    }
    #[test]
    fn test_speed_aero_map_boundary_returns_exact_value() {
        let map = SpeedAeroMap {
            velocities: vec![10.0, 50.0, 100.0],
            cds: vec![0.30, 0.32, 0.35],
            cls: vec![-0.2, -0.4, -0.6],
        };
        assert!((map.interpolate_cd(10.0) - 0.30).abs() < 1e-12);
        assert!((map.interpolate_cd(100.0) - 0.35).abs() < 1e-12);
        assert!((map.interpolate_cl(10.0) - (-0.2)).abs() < 1e-12);
        assert!((map.interpolate_cl(100.0) - (-0.6)).abs() < 1e-12);
    }
    #[test]
    fn test_aero_body_forces_drag_opposes_velocity() {
        let body = AeroBody {
            cda: 0.6,
            cla: -0.5,
            csa: 0.0,
            cop: [0.0; 3],
        };
        let vel = [30.0_f64, 0.0, 0.0];
        let forces = body.forces(vel, AIR_DENSITY_SEA_LEVEL);
        assert!(
            forces.drag[0] < 0.0,
            "drag Fx should be negative, got {}",
            forces.drag[0]
        );
    }
    #[test]
    fn test_aero_forces_total_sums_components() {
        let f = AeroForces {
            drag: [1.0, 2.0, 3.0],
            lift: [0.1, 0.2, 0.3],
            side: [-0.5, 0.0, 0.5],
        };
        let t = f.total();
        assert!((t[0] - 0.6).abs() < 1e-12);
        assert!((t[1] - 2.2).abs() < 1e-12);
        assert!((t[2] - 3.8).abs() < 1e-12);
    }
    #[test]
    fn test_compute_drag_zero_velocity_returns_zero() {
        let f = compute_drag([0.0, 0.0, 0.0], AIR_DENSITY_SEA_LEVEL, 0.6);
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn test_compute_drag_magnitude_quadratic() {
        let cda = 0.6;
        let rho = AIR_DENSITY_SEA_LEVEL;
        let f1 = compute_drag([20.0, 0.0, 0.0], rho, cda);
        let f2 = compute_drag([40.0, 0.0, 0.0], rho, cda);
        let ratio = f2[0] / f1[0];
        assert!((ratio - 4.0).abs() < 1e-9, "expected ratio 4, got {ratio}");
    }
    #[test]
    fn test_compute_lift_downforce_negative_z() {
        let f = compute_lift([30.0, 0.0, 0.0], AIR_DENSITY_SEA_LEVEL, -1.0, 0.0);
        assert!(f[2] < 0.0, "downforce should be negative Fz, got {}", f[2]);
    }
    #[test]
    fn test_compute_lift_cop_offset_scales_magnitude() {
        let base = compute_lift([20.0, 0.0, 0.0], AIR_DENSITY_SEA_LEVEL, -0.8, 0.0);
        let boosted = compute_lift([20.0, 0.0, 0.0], AIR_DENSITY_SEA_LEVEL, -0.8, 0.5);
        let ratio = boosted[2] / base[2];
        assert!(
            (ratio - 1.5).abs() < 1e-9,
            "expected ratio 1.5, got {ratio}"
        );
    }
    #[test]
    fn test_downforce_package_min_level() {
        let pkg = DownforcePackage {
            front_level: 1,
            rear_level: 1,
            max_level: 7,
            cl_min: -0.5,
            cl_max: -2.0,
        };
        assert!((pkg.front_cla() - (-0.5)).abs() < 1e-9);
    }
    #[test]
    fn test_downforce_package_max_level() {
        let pkg = DownforcePackage {
            front_level: 7,
            rear_level: 7,
            max_level: 7,
            cl_min: -0.5,
            cl_max: -2.0,
        };
        assert!((pkg.rear_cla() - (-2.0)).abs() < 1e-9);
    }
    #[test]
    fn test_downforce_package_total_cla() {
        let pkg = DownforcePackage {
            front_level: 1,
            rear_level: 7,
            max_level: 7,
            cl_min: -0.5,
            cl_max: -2.0,
        };
        let total = pkg.total_cla();
        assert!(
            total < pkg.front_cla(),
            "total should be more negative than front alone"
        );
    }
    #[test]
    fn test_wind_model_no_gust_returns_mean() {
        let wind = WindModel {
            mean_velocity: [5.0, 0.0, 0.0],
            gust_amplitude: 0.0,
            gust_wavelength: 100.0,
        };
        let v = wind.velocity_at([0.0, 0.0, 0.0]);
        assert!((v[0] - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_wind_model_relative_wind_subtracted() {
        let wind = WindModel {
            mean_velocity: [10.0, 0.0, 0.0],
            gust_amplitude: 0.0,
            gust_wavelength: 100.0,
        };
        let rel = wind.relative_wind([0.0; 3], [6.0, 0.0, 0.0]);
        assert!(
            (rel[0] - 4.0).abs() < 1e-9,
            "relative wind should be 4 m/s, got {}",
            rel[0]
        );
    }
    #[test]
    fn test_aero_map_constant_returns_same_everywhere() {
        let map = AeroMap::from_constant(0.35, -0.8, 0.0);
        let (cd, cl, cs) = map.interpolate(0.0, 0.0);
        assert!((cd - 0.35).abs() < 1e-12);
        assert!((cl - (-0.8)).abs() < 1e-12);
        assert!(cs.abs() < 1e-12);
    }
    #[test]
    fn test_aero_map_interpolation_midpoint() {
        let map = AeroMap {
            alphas: vec![0.0, 10.0],
            betas: vec![0.0],
            cd: vec![0.3, 0.4],
            cl: vec![-0.5, -0.8],
            cs: vec![0.0, 0.0],
        };
        let (cd, cl, _) = map.interpolate(5.0, 0.0);
        assert!(
            (cd - 0.35).abs() < 1e-9,
            "expected Cd=0.35 at mid, got {cd}"
        );
        assert!(
            (cl - (-0.65)).abs() < 1e-9,
            "expected Cl=-0.65 at mid, got {cl}"
        );
    }
    #[test]
    fn test_ground_effect_advanced_venturi_suction_increases_at_lower_height() {
        let ge_low = GroundEffectAdvanced {
            ride_height: 0.03,
            reference_height: 0.12,
            max_amplification: 4.0,
            diffuser_angle: 15.0,
        };
        let ge_high = GroundEffectAdvanced {
            ride_height: 0.10,
            reference_height: 0.12,
            max_amplification: 4.0,
            diffuser_angle: 15.0,
        };
        assert!(
            ge_low.venturi_suction_factor() > ge_high.venturi_suction_factor(),
            "lower ride height should give more suction"
        );
    }
    #[test]
    fn test_ground_effect_advanced_amplified_cl_negative() {
        let ge = GroundEffectAdvanced::new_race_car();
        let base_cl = -2.0;
        let amp = ge.amplified_cl(base_cl);
        assert!(
            amp < base_cl,
            "amplified Cl should be more negative: amp={amp}"
        );
    }
    #[test]
    fn test_ground_effect_advanced_stall_at_extreme_height() {
        let ge = GroundEffectAdvanced {
            ride_height: 0.01,
            ..GroundEffectAdvanced::new_race_car()
        };
        assert!(ge.is_stalled(), "very low ride height should be stalled");
    }
    #[test]
    fn test_wake_deficit_zero_at_large_distance() {
        let wake = WakeInteraction {
            following_distance: 1000.0,
            lateral_offset: 0.0,
            decay_length: 3.0,
            leading_frontal_area: 2.0,
        };
        let deficit = wake.velocity_deficit(50.0);
        assert!(
            deficit < 0.01,
            "far behind should have near-zero deficit: {deficit}"
        );
    }
    #[test]
    fn test_wake_deficit_larger_when_closer() {
        let mk_wake = |d: f64| WakeInteraction {
            following_distance: d,
            lateral_offset: 0.0,
            decay_length: 3.0,
            leading_frontal_area: 2.0,
        };
        let close = mk_wake(2.0).velocity_deficit(50.0);
        let far = mk_wake(20.0).velocity_deficit(50.0);
        assert!(close > far, "closer following should have larger deficit");
    }
    #[test]
    fn test_wake_downforce_loss_fraction_in_range() {
        let wake = WakeInteraction {
            following_distance: 5.0,
            lateral_offset: 0.5,
            decay_length: 2.5,
            leading_frontal_area: 1.8,
        };
        let loss = wake.downforce_loss_fraction(40.0);
        assert!(
            (0.0..=1.0).contains(&loss),
            "loss fraction must be in [0,1]: {loss}"
        );
    }
    #[test]
    fn test_downforce_map_constant_returns_constant() {
        let map = DownforceMap::from_constant(-2.5);
        let v = map.interpolate(50.0, 0.05);
        assert!(
            (v - (-2.5)).abs() < 1e-9,
            "constant map should return -2.5, got {v}"
        );
    }
    #[test]
    fn test_downforce_map_parametric_increases_at_lower_height() {
        let map = DownforceMap::parametric(
            -2.0,
            50.0,
            0.1,
            0.5,
            vec![30.0, 50.0, 80.0],
            vec![0.03, 0.06, 0.10],
        );
        let high_h = map.interpolate(50.0, 0.10);
        let low_h = map.interpolate(50.0, 0.03);
        assert!(
            low_h < high_h,
            "lower ride height should give more negative ClA: low={low_h}, high={high_h}"
        );
    }
    #[test]
    fn test_balance_controller_increases_front_wing_when_understeer() {
        let ctrl = AeroBalanceController {
            target_front_fraction: 0.40,
            gain: 2.0,
            max_delta_level: 1.0,
        };
        let balance = AeroBalance {
            front_downforce: 500.0,
            rear_downforce: 1500.0,
        };
        let delta = ctrl.front_wing_correction(&balance);
        assert!(delta > 0.0, "should increase front wing, got delta={delta}");
    }
    #[test]
    fn test_balance_controller_apply_does_not_exceed_max_level() {
        let ctrl = AeroBalanceController {
            target_front_fraction: 0.40,
            gain: 100.0,
            max_delta_level: 3.0,
        };
        let balance = AeroBalance {
            front_downforce: 100.0,
            rear_downforce: 1900.0,
        };
        let mut pkg = DownforcePackage {
            front_level: 7,
            rear_level: 7,
            max_level: 7,
            cl_min: -0.5,
            cl_max: -2.0,
        };
        ctrl.apply_correction(&mut pkg, &balance);
        assert!(pkg.front_level <= 7, "front_level must not exceed max");
    }
    #[test]
    fn test_wing_stall_linear_below_stall() {
        let stall = WingStallModel {
            lift_slope: 2.0 * std::f64::consts::PI,
            alpha_stall: 15.0,
            cl_max: 1.5,
            cl_post_stall_fraction: 0.6,
            leading_edge_stall: true,
        };
        let cl_5 = stall.lift_coefficient(5.0);
        let cl_10 = stall.lift_coefficient(10.0);
        let ratio = cl_10 / cl_5;
        assert!(
            (ratio - 2.0).abs() < 0.01,
            "below stall ratio should be ~2, got {ratio}"
        );
    }
    #[test]
    fn test_wing_stall_abrupt_drops_above_stall() {
        let stall = WingStallModel {
            lift_slope: 2.0 * PI,
            alpha_stall: 15.0,
            cl_max: 1.5,
            cl_post_stall_fraction: 0.5,
            leading_edge_stall: true,
        };
        let cl_pre = stall.lift_coefficient(14.9);
        let cl_post = stall.lift_coefficient(20.0);
        assert!(
            cl_post < cl_pre,
            "abrupt stall should drop CL: pre={cl_pre}, post={cl_post}"
        );
    }
    #[test]
    fn test_wing_stall_trailing_edge_gradual() {
        let stall = WingStallModel {
            lift_slope: 2.0 * PI,
            alpha_stall: 15.0,
            cl_max: 1.5,
            cl_post_stall_fraction: 0.7,
            leading_edge_stall: false,
        };
        let cl_15 = stall.lift_coefficient(15.0);
        let cl_20 = stall.lift_coefficient(20.0);
        let cl_25 = stall.lift_coefficient(25.0);
        assert!(
            cl_20 <= cl_15 + 1e-3,
            "should not increase after stall onset"
        );
        assert!(cl_25 <= cl_20 + 1e-3, "gradual decay should continue");
    }
    #[test]
    fn test_cfd_surrogate_exact_at_centre() {
        let centers = vec![vec![0.0_f64, 0.0], vec![10.0, 0.0], vec![0.0, 10.0]];
        let values = vec![1.0, 2.0, 3.0];
        let model = CfdSurrogateModel::new(centers, values, 1.0);
        let pred = model.predict(&[0.0, 0.0]);
        assert!(
            (pred - 1.0).abs() < 0.2,
            "prediction near first centre: {pred}"
        );
    }
    #[test]
    fn test_cfd_surrogate_interpolation_between_centres() {
        let centers = vec![vec![0.0_f64], vec![1.0]];
        let values = vec![0.0, 10.0];
        let model = CfdSurrogateModel::new(centers, values, 1.0);
        let pred = model.predict(&[0.5]);
        assert!(
            (0.0..=10.0).contains(&pred),
            "midpoint prediction out of range: {pred}"
        );
    }
    #[test]
    fn test_pg_factor_unity_at_low_mach() {
        let cc = CompressibilityCorrection::sea_level();
        let f = cc.pg_factor(10.0);
        assert!((f - 1.0).abs() < 0.01, "PG factor near 1 at low speed: {f}");
    }
    #[test]
    fn test_pg_factor_increases_with_mach() {
        let cc = CompressibilityCorrection::sea_level();
        let f_low = cc.pg_factor(50.0);
        let f_high = cc.pg_factor(250.0);
        assert!(
            f_high > f_low,
            "PG factor should grow with Mach: low={f_low}, high={f_high}"
        );
    }
    #[test]
    fn test_corrected_cl_exceeds_base() {
        let cc = CompressibilityCorrection::sea_level();
        let cl0 = 1.0;
        let cl_c = cc.corrected_cl(cl0, 200.0);
        assert!(cl_c > cl0, "compressibility correction should increase Cl");
    }
    #[test]
    fn test_diffuser_area_ratio_greater_than_one() {
        let diff = DiffuserModel {
            inlet_height: 0.05,
            exit_height: 0.20,
            length: 0.80,
            width: 1.50,
            efficiency: 0.70,
        };
        assert!(diff.area_ratio() > 1.0, "exit larger than inlet → AR > 1");
    }
    #[test]
    fn test_diffuser_downforce_increases_with_speed() {
        let diff = DiffuserModel {
            inlet_height: 0.05,
            exit_height: 0.20,
            length: 0.80,
            width: 1.50,
            efficiency: 0.70,
        };
        let f1 = diff.downforce(30.0, AIR_DENSITY_SEA_LEVEL);
        let f2 = diff.downforce(60.0, AIR_DENSITY_SEA_LEVEL);
        assert!(f2 > f1, "downforce should increase with speed");
        assert!(
            (f2 / f1 - 4.0).abs() < 0.1,
            "quadratic scaling: ratio={}",
            f2 / f1
        );
    }
    #[test]
    fn test_diffuser_expansion_angle_nonzero() {
        let diff = DiffuserModel {
            inlet_height: 0.05,
            exit_height: 0.20,
            length: 0.80,
            width: 1.50,
            efficiency: 0.70,
        };
        let angle = diff.expansion_angle();
        assert!(angle > 0.0, "expansion angle must be positive: {angle}");
    }
}
/// Helper: find lower index and fractional position in a sorted table.
///
/// Returns `(index, fraction)` where `table[index] <= value < table[index+1]`.
pub(super) fn interp_index(table: &[f64], value: f64) -> (usize, f64) {
    let n = table.len();
    if n == 0 {
        return (0, 0.0);
    }
    if value <= table[0] {
        return (0, 0.0);
    }
    if value >= table[n - 1] {
        return (n - 1, 0.0);
    }
    for i in 0..n - 1 {
        if value >= table[i] && value < table[i + 1] {
            let dv = table[i + 1] - table[i];
            let frac = if dv.abs() > 1e-15 {
                (value - table[i]) / dv
            } else {
                0.0
            };
            return (i, frac);
        }
    }
    (n - 2, 1.0)
}
#[cfg(test)]
mod tests_aero_new {

    use crate::AeroBody;
    use crate::AeroDatabase;

    use crate::AeroSurface;

    #[test]
    fn test_drag_polar_zero_alpha_equals_cd0() {
        let body = AeroBody {
            cda: 0.6,
            cla: -1.2,
            csa: 0.0,
            cop: [0.0; 3],
        };
        let cd0 = 0.025;
        let (_, cd) = body.compute_drag_polar(0.0, 6.0, cd0);
        assert!(
            (cd - cd0).abs() < 1e-9,
            "Cd at alpha=0 should equal Cd0, got {cd}"
        );
    }
    #[test]
    fn test_drag_polar_increases_with_alpha() {
        let body = AeroBody {
            cda: 0.6,
            cla: -1.2,
            csa: 0.0,
            cop: [0.0; 3],
        };
        let cd0 = 0.02;
        let (_, cd1) = body.compute_drag_polar(5.0, 6.0, cd0);
        let (_, cd2) = body.compute_drag_polar(10.0, 6.0, cd0);
        assert!(
            cd2 > cd1,
            "Cd should increase with |alpha|: cd1={cd1}, cd2={cd2}"
        );
    }
    #[test]
    fn test_drag_polar_symmetric_alpha() {
        let body = AeroBody {
            cda: 0.6,
            cla: 0.0,
            csa: 0.0,
            cop: [0.0; 3],
        };
        let (_, cd_pos) = body.compute_drag_polar(8.0, 5.0, 0.02);
        let (_, cd_neg) = body.compute_drag_polar(-8.0, 5.0, 0.02);
        assert!(
            (cd_pos - cd_neg).abs() < 1e-9,
            "Cd should be symmetric: {cd_pos} vs {cd_neg}"
        );
    }
    #[test]
    fn test_induced_drag_zero_at_zero_cl() {
        let cdi = AeroBody::compute_induced_drag(0.0, 8.0, 0.9);
        assert_eq!(cdi, 0.0, "zero lift → zero induced drag");
    }
    #[test]
    fn test_induced_drag_positive_for_nonzero_cl() {
        let cdi = AeroBody::compute_induced_drag(1.0, 6.0, 0.85);
        assert!(
            cdi > 0.0,
            "nonzero Cl should produce positive CDi, got {cdi}"
        );
    }
    #[test]
    fn test_induced_drag_decreases_with_higher_ar() {
        let cdi_low = AeroBody::compute_induced_drag(1.0, 4.0, 0.9);
        let cdi_high = AeroBody::compute_induced_drag(1.0, 10.0, 0.9);
        assert!(
            cdi_high < cdi_low,
            "higher AR reduces induced drag: low={cdi_low}, high={cdi_high}"
        );
    }
    #[test]
    fn test_vlm_returns_correct_panel_count() {
        let surf = AeroSurface::new(2.0, 0.3, 8, 5.0, 0.85);
        let dist = surf.vortex_lattice_method(50.0);
        assert_eq!(dist.len(), 8, "should return 8 panel entries");
    }
    #[test]
    fn test_vlm_elliptic_peak_at_centre() {
        let surf = AeroSurface::new(2.0, 0.3, 10, 8.0, 0.85);
        let dist = surf.vortex_lattice_method(60.0);
        let cl_centre = (dist[4].1 + dist[5].1) * 0.5;
        let cl_tip = dist[0].1;
        assert!(
            cl_centre > cl_tip,
            "elliptic: centre Cl > tip Cl: centre={cl_centre}, tip={cl_tip}"
        );
    }
    #[test]
    fn test_vlm_total_cl_positive_for_positive_alpha() {
        let surf = AeroSurface::new(1.8, 0.25, 10, 6.0, 0.85);
        let cl_total = surf.total_lift_coefficient(50.0);
        assert!(
            cl_total > 0.0,
            "positive alpha → positive total Cl: {cl_total}"
        );
    }
    #[test]
    fn test_aero_db_exact_lookup() {
        let alphas = vec![-5.0, 0.0, 5.0, 10.0];
        let cds = vec![0.04, 0.02, 0.04, 0.10];
        let db = AeroDatabase::from_alpha_cd(alphas, cds);
        let cd = db.interpolate_cd(0.0, 1.0e6);
        assert!((cd - 0.02).abs() < 1e-6, "exact lookup at alpha=0: {cd}");
    }
    #[test]
    fn test_aero_db_interpolation_midpoint() {
        let alphas = vec![0.0, 10.0];
        let cds = vec![0.02, 0.10];
        let db = AeroDatabase::from_alpha_cd(alphas, cds);
        let cd = db.interpolate_cd(5.0, 1.0e6);
        let expected = 0.06;
        assert!(
            (cd - expected).abs() < 1e-6,
            "midpoint interpolation: expected {expected}, got {cd}"
        );
    }
    #[test]
    fn test_aero_db_clamp_below_range() {
        let alphas = vec![0.0, 5.0, 10.0];
        let cds = vec![0.02, 0.04, 0.08];
        let db = AeroDatabase::from_alpha_cd(alphas, cds);
        let cd_clamp = db.interpolate_cd(-10.0, 1.0e6);
        let cd_zero = db.interpolate_cd(0.0, 1.0e6);
        assert!(
            (cd_clamp - cd_zero).abs() < 1e-9,
            "below range should clamp to first entry"
        );
    }
}
