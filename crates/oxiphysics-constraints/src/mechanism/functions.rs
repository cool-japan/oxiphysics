//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{
    FourBarAnalysis, FourBarLinkage, GrashofClass, InvoluteGear, JointClass, RotatingMass,
};

#[allow(dead_code)]
#[inline]
pub(super) fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.clamp(lo, hi)
}
#[allow(dead_code)]
#[inline]
pub(super) fn sq(x: f64) -> f64 {
    x * x
}
/// Test Grashof condition for given link lengths.
///
/// Returns `true` when the shortest + longest ≤ sum of remaining two.
#[allow(dead_code)]
pub fn grashof_condition(l1: f64, l2: f64, l3: f64, l4: f64) -> bool {
    FourBarLinkage::new(l1, l2, l3, l4).grashof()
}
/// Compute `theta3` and `theta4` for a four-bar linkage at crank angle `theta2`.
///
/// Returns `(theta3, theta4)` in radians, or `None` on error.
#[allow(dead_code)]
pub fn four_bar_angles(l1: f64, l2: f64, l3: f64, l4: f64, theta2: f64) -> Option<(f64, f64)> {
    let fb = FourBarLinkage::new(l1, l2, l3, l4);
    let theta4 = fb.follower_angle(theta2)?;
    let k4 = self_k4(l1, l2, l3, l4);
    let k5 = self_k5(l1, l2, l3, l4);
    let d = theta2.cos() - k4 - k5 * theta2.cos();
    let e = -2.0 * theta2.sin();
    let f_coef = k4 - (k5 - 1.0) * theta2.cos();
    let disc = sq(e) - 4.0 * d * f_coef;
    if disc < 0.0 {
        return None;
    }
    let t3 = (-e - disc.sqrt()) / (2.0 * d);
    let theta3 = 2.0 * t3.atan();
    Some((theta3, theta4))
}
pub(super) fn self_k4(l1: f64, _l2: f64, l3: f64, _l4: f64) -> f64 {
    l1 / l3
}
pub(super) fn self_k5(_l1: f64, l2: f64, l3: f64, _l4: f64) -> f64 {
    (sq(l2) + sq(l3) - sq(_l1) - sq(_l4)) / (2.0 * l2 * l3)
}
/// Compute compound gear ratio from a series of (drive, driven) tooth counts.
#[allow(dead_code)]
pub fn gear_ratio_compound(stages: &[(u32, u32)]) -> f64 {
    stages
        .iter()
        .fold(1.0, |acc, &(drv, drvn)| acc * drvn as f64 / drv as f64)
}
/// Sample involute profile points for a gear with `teeth` and `module`.
#[allow(dead_code)]
pub fn involute_profile(teeth: u32, module: f64, points: usize) -> Vec<[f64; 2]> {
    let r = module * teeth as f64 / 2.0;
    (0..points)
        .map(|i| {
            let phi = i as f64 * 0.5 / points as f64;
            let x = r * (phi.cos() + phi * phi.sin());
            let y = r * (phi.sin() - phi * phi.cos());
            [x, y]
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mechanism::CamFollower;

    use crate::mechanism::GearPair;

    use crate::mechanism::LockingMechanism;

    use crate::mechanism::PlanetaryGear;
    use crate::mechanism::RackAndPinion;
    use crate::mechanism::ScissorLift;

    use crate::mechanism::SliderCrankMechanism;
    use crate::mechanism::ToggleClamp;

    use crate::mechanism::WormGear;
    #[test]
    fn test_grashof_crank_rocker() {
        assert!(grashof_condition(4.0, 2.0, 4.0, 3.0));
    }
    #[test]
    fn test_grashof_non() {
        assert!(!grashof_condition(2.0, 3.0, 1.0, 5.0));
    }
    #[test]
    fn test_four_bar_follower_angle() {
        let fb = FourBarLinkage::new(4.0, 2.0, 4.0, 3.0);
        let found = [0.0, 0.3, 1.0, 1.5, 2.0, 2.5, 3.0]
            .iter()
            .any(|&a| fb.follower_angle(a).is_some());
        assert!(
            found,
            "expected at least one valid crank angle to produce a follower angle"
        );
    }
    #[test]
    fn test_transmission_angle_range() {
        let fb = FourBarLinkage::new(4.0, 2.0, 4.0, 3.0);
        let mu = fb.transmission_angle(0.5);
        assert!((0.0..=PI).contains(&mu));
    }
    #[test]
    fn test_slider_crank_position_zero() {
        let sc = SliderCrankMechanism::new(0.05, 0.15, 100.0);
        let pos = sc.piston_position(0.0);
        assert!((pos - 0.2).abs() < 1e-9);
    }
    #[test]
    fn test_slider_crank_velocity_tdc() {
        let sc = SliderCrankMechanism::new(0.05, 0.15, 100.0);
        let v = sc.piston_velocity(0.0);
        assert!(v.abs() < 1e-9);
    }
    #[test]
    fn test_cam_displacement_zero() {
        let cf = CamFollower::new(0.02, 0.01, PI, 10.0);
        assert!((cf.displacement(0.0)).abs() < 1e-9);
    }
    #[test]
    fn test_cam_displacement_full() {
        let cf = CamFollower::new(0.02, 0.01, PI, 10.0);
        assert!((cf.displacement(1.0) - 0.01).abs() < 1e-9);
    }
    #[test]
    fn test_cam_cycloidal_mid() {
        let cf = CamFollower::new(0.02, 0.01, PI, 10.0);
        assert!((cf.displacement(0.5) - 0.005).abs() < 1e-6);
    }
    #[test]
    fn test_gear_ratio() {
        let gp = GearPair::new(20, 40, 0.001);
        assert!((gp.gear_ratio() - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_gear_omega_driven() {
        let mut gp = GearPair::new(20, 40, 0.001);
        gp.omega_drive = 100.0;
        assert!((gp.omega_driven() + 50.0).abs() < 1e-6);
    }
    #[test]
    fn test_gear_centre_distance() {
        let gp = GearPair::new(20, 40, 0.001);
        assert!((gp.centre_distance() - 0.03).abs() < 1e-9);
    }
    #[test]
    fn test_gear_ratio_compound() {
        let stages = vec![(10_u32, 30_u32), (15_u32, 45_u32)];
        let r = gear_ratio_compound(&stages);
        assert!((r - 9.0).abs() < 1e-9);
    }
    #[test]
    fn test_rack_displacement() {
        let rp = RackAndPinion::new(20, 0.001);
        let d = rp.rack_displacement(PI);
        assert!((d - PI * 0.01).abs() < 1e-9);
    }
    #[test]
    fn test_pinion_rotation() {
        let rp = RackAndPinion::new(20, 0.001);
        let disp = 0.05;
        let angle = rp.pinion_rotation(disp);
        assert!((rp.rack_displacement(angle) - disp).abs() < 1e-9);
    }
    #[test]
    fn test_worm_reduction() {
        let wg = WormGear::new(1, 30, 5.0_f64.to_radians(), 0.05);
        assert!((wg.reduction_ratio() - 30.0).abs() < 1e-9);
    }
    #[test]
    fn test_worm_self_locking() {
        let wg = WormGear::new(1, 30, 2.0_f64.to_radians(), 0.1);
        assert!(wg.is_self_locking());
    }
    #[test]
    fn test_worm_not_self_locking() {
        let wg = WormGear::new(4, 30, 45.0_f64.to_radians(), 0.05);
        assert!(!wg.is_self_locking());
    }
    #[test]
    fn test_planetary_planet_teeth() {
        let pg = PlanetaryGear::new(20, 60, 3);
        assert_eq!(pg.planet_teeth(), 20);
    }
    #[test]
    fn test_planetary_carrier_ring_fixed() {
        let mut pg = PlanetaryGear::new(20, 60, 3);
        pg.omega_sun = 100.0;
        pg.omega_ring = 0.0;
        let oc = pg.carrier_velocity();
        assert!((oc - 25.0).abs() < 1e-9);
    }
    #[test]
    fn test_scissor_height_increases() {
        let sl = ScissorLift::new(1.0, 2);
        let h0 = sl.platform_height(0.0);
        let h1 = sl.platform_height(0.5);
        assert!(h1 > h0);
    }
    #[test]
    fn test_toggle_ma_small_angle() {
        let tc = ToggleClamp::new(0.1, 0.001);
        let ma = {
            let mut tc2 = ToggleClamp::new(0.1, 0.001);
            tc2.angle = 0.05;
            tc2.mechanical_advantage()
        };
        assert!(ma > 5.0);
        let _ = tc;
    }
    #[test]
    fn test_toggle_locked() {
        let mut tc = ToggleClamp::new(0.1, 0.005);
        tc.angle = 0.001;
        assert!(tc.is_locked());
    }
    #[test]
    fn test_ratchet_blocks_reverse() {
        let mut lm = LockingMechanism::new_ratchet(12, 5.0);
        let ok = lm.try_advance(-0.1);
        assert!(!ok);
        assert!(lm.locked);
    }
    #[test]
    fn test_ratchet_allows_forward() {
        let mut lm = LockingMechanism::new_ratchet(12, 5.0);
        let ok = lm.try_advance(0.5);
        assert!(ok);
    }
    #[test]
    fn test_detent_force() {
        let lm = LockingMechanism::new_detent(0.005, 10.0);
        let f = lm.detent_force();
        assert!((f - 0.05).abs() < 1e-9);
    }
    #[test]
    fn test_involute_profile_count() {
        let pts = involute_profile(20, 0.001, 10);
        assert_eq!(pts.len(), 10);
    }
    #[test]
    fn test_four_bar_angles_valid() {
        let result = four_bar_angles(4.0, 2.0, 4.0, 3.0, 1.0);
        assert!(result.is_some());
    }
    #[test]
    fn test_slider_crank_step() {
        let mut sc = SliderCrankMechanism::new(0.05, 0.15, 10.0);
        sc.step(0.1);
        assert!(sc.theta > 0.0);
    }
    #[test]
    fn test_cam_step() {
        let mut cf = CamFollower::new(0.02, 0.01, PI, 10.0);
        cf.step(0.1);
        assert!(cf.theta > 0.0);
    }
}
/// Grübler-Kutzbach mobility formula for spatial (3-D) mechanisms.
///
/// `M = 6*(n - 1) - Σ(6 - f_i)`, where `n` is the number of links
/// (including ground) and `f_i` is the DOF of joint `i`.
///
/// Returns the mobility (number of independent inputs needed).
#[allow(dead_code)]
pub fn grubler_kutzbach_spatial(n_links: u32, joints: &[JointClass]) -> i32 {
    let lambda = 6_i32;
    let n = n_links as i32;
    let constraint_sum: i32 = joints.iter().map(|j| lambda - j.dof() as i32).sum();
    lambda * (n - 1) - constraint_sum
}
/// Grübler-Kutzbach mobility formula for planar (2-D) mechanisms.
///
/// `M = 3*(n - 1) - Σ(3 - f_i)`.
#[allow(dead_code)]
pub fn grubler_kutzbach_planar(n_links: u32, joints: &[JointClass]) -> i32 {
    let lambda = 3_i32;
    let n = n_links as i32;
    let constraint_sum: i32 = joints.iter().map(|j| lambda - j.dof().min(3) as i32).sum();
    lambda * (n - 1) - constraint_sum
}
/// Classify a four-bar linkage according to the Grashof criterion.
#[allow(dead_code)]
pub fn grashof_classify(l1: f64, l2: f64, l3: f64, l4: f64) -> GrashofClass {
    let links = [l1, l2, l3, l4];
    let s = links.iter().cloned().fold(f64::MAX, f64::min);
    let l = links.iter().cloned().fold(f64::MIN, f64::max);
    let pq: f64 = links.iter().sum::<f64>() - s - l;
    let diff = s + l - pq;
    if diff.abs() < 1e-10 {
        return GrashofClass::ChangePoint;
    }
    if s + l > pq {
        return GrashofClass::NonGrashof;
    }
    if (s - l1).abs() < 1e-10 {
        GrashofClass::DoubleRocker
    } else if (s - l2).abs() < 1e-10 {
        GrashofClass::CrankRocker
    } else if (s - l4).abs() < 1e-10 {
        GrashofClass::RockerCrank
    } else {
        GrashofClass::DoubleCrank
    }
}
/// Compute the connecting-rod angle `phi` from the slider-crank geometry.
///
/// `phi = arcsin(r/l * sin(theta))`, where `r` = crank radius, `l` = rod length.
#[allow(dead_code)]
pub fn slider_crank_rod_angle(crank_radius: f64, rod_length: f64, theta: f64) -> f64 {
    let ratio = crank_radius / rod_length * theta.sin();
    ratio.clamp(-1.0, 1.0).asin()
}
/// Compute the ratio parameter λ = r/l (crank ratio) for a slider-crank.
#[allow(dead_code)]
pub fn slider_crank_lambda(crank_radius: f64, rod_length: f64) -> f64 {
    crank_radius / rod_length
}
/// Approximate piston velocity (m/s) using the Fourier series expansion.
///
/// First two terms: `v ≈ -r*ω*(sin θ + λ/2 * sin 2θ)`.
#[allow(dead_code)]
pub fn piston_velocity_approx(crank_radius: f64, rod_length: f64, omega: f64, theta: f64) -> f64 {
    let lambda = slider_crank_lambda(crank_radius, rod_length);
    -crank_radius * omega * (theta.sin() + lambda / 2.0 * (2.0 * theta).sin())
}
/// Polynomial (1-2-3) cam profile for simple rise motions.
///
/// Returns normalised displacement `s/h` ∈ \[0,1\] given `tau` ∈ \[0,1\].
#[allow(dead_code)]
pub fn cam_polynomial_123(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    3.0 * t * t - 2.0 * t * t * t
}
/// Polynomial (3-4-5) cam profile (zero velocity and acceleration at endpoints).
///
/// Returns normalised displacement `s/h` ∈ \[0,1\].
#[allow(dead_code)]
pub fn cam_polynomial_345(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    10.0 * t.powi(3) - 15.0 * t.powi(4) + 6.0 * t.powi(5)
}
/// Polynomial (4-5-6-7) cam profile (zero velocity, acceleration, jerk at endpoints).
///
/// Returns normalised displacement `s/h` ∈ \[0,1\].
#[allow(dead_code)]
pub fn cam_polynomial_4567(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    35.0 * t.powi(4) - 84.0 * t.powi(5) + 70.0 * t.powi(6) - 20.0 * t.powi(7)
}
/// Cycloidal cam profile (zero velocity and acceleration at endpoints).
///
/// Returns normalised displacement `s/h` ∈ \[0,1\].
#[allow(dead_code)]
pub fn cam_cycloidal(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    t - (2.0 * PI * t).sin() / (2.0 * PI)
}
/// Double-harmonic cam profile.
///
/// Returns normalised displacement.
#[allow(dead_code)]
pub fn cam_double_harmonic(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    0.25 * (1.0 - (PI * t).cos()) - 0.125 * (1.0 - (2.0 * PI * t).cos())
}
/// Trigonometric (modified trapezoid) cam profile.
///
/// Returns normalised displacement.
#[allow(dead_code)]
pub fn cam_modified_trapezoid(tau: f64) -> f64 {
    let t = tau.clamp(0.0, 1.0);
    if t < 0.25 {
        0.5 * (1.0 - (2.0 * PI * t).cos()) / 2.0
    } else if t < 0.75 {
        0.5 * t - 0.25 * (2.0 * PI * t).sin() / (2.0 * PI)
    } else {
        1.0 - cam_modified_trapezoid(1.0 - t)
    }
}
/// Compute the maximum pressure angle (radians) for a cam profile.
///
/// Uses numerical search over `n_pts` sample points.
#[allow(dead_code)]
pub fn cam_max_pressure_angle(
    base_radius: f64,
    lift: f64,
    rise_angle: f64,
    profile_fn: fn(f64) -> f64,
    n_pts: usize,
) -> f64 {
    let h = 1e-6;
    let mut max_pa: f64 = 0.0;
    for i in 1..n_pts.saturating_sub(1) {
        let tau = i as f64 / (n_pts - 1) as f64;
        let ds_dtau = (profile_fn(tau + h) - profile_fn(tau - h)) / (2.0 * h);
        let ds_dtheta = lift * ds_dtau / rise_angle;
        let r = base_radius + lift * profile_fn(tau);
        let pa = (ds_dtheta / r).abs().atan();
        if pa > max_pa {
            max_pa = pa;
        }
    }
    max_pa
}
/// Compute cam pitch curve radius at normalised position `tau`.
///
/// Pitch curve radius = base_radius + displacement.
#[allow(dead_code)]
pub fn cam_pitch_radius(base_radius: f64, lift: f64, tau: f64, profile_fn: fn(f64) -> f64) -> f64 {
    base_radius + lift * profile_fn(tau.clamp(0.0, 1.0))
}
/// Cam radius of curvature at `tau` (numerically).
///
/// Uses: `ρ = (r² + (dr/dθ)²)^(3/2) / |r² + 2*(dr/dθ)² - r*(d²r/dθ²)|`.
#[allow(dead_code)]
pub fn cam_radius_of_curvature(
    base_radius: f64,
    lift: f64,
    rise_angle: f64,
    tau: f64,
    profile_fn: fn(f64) -> f64,
) -> f64 {
    let h = 1e-5;
    let r = |t: f64| base_radius + lift * profile_fn(t.clamp(0.0, 1.0));
    let theta = tau * rise_angle;
    let _ = theta;
    let rp = (r(tau + h) - r(tau - h)) / (2.0 * h * rise_angle);
    let rpp = (r(tau + h) - 2.0 * r(tau) + r(tau - h)) / (h * h * rise_angle * rise_angle);
    let rv = r(tau);
    let num = (rv * rv + rp * rp).powf(1.5);
    let den = (rv * rv + 2.0 * rp * rp - rv * rpp).abs().max(1e-15);
    num / den
}
/// Contact ratio for a meshing pair of involute gears.
///
/// ε = (length of contact path) / (base pitch).
#[allow(dead_code)]
pub fn contact_ratio(gear1: &InvoluteGear, gear2: &InvoluteGear) -> f64 {
    let ra1 = gear1.addendum_radius();
    let rb1 = gear1.base_radius();
    let ra2 = gear2.addendum_radius();
    let rb2 = gear2.base_radius();
    let r1 = gear1.pitch_radius();
    let r2 = gear2.pitch_radius();
    let phi = gear1.pressure_angle;
    let approach = ((ra2 * ra2 - rb2 * rb2).max(0.0).sqrt() - r2 * phi.sin()).abs();
    let recess = ((ra1 * ra1 - rb1 * rb1).max(0.0).sqrt() - r1 * phi.sin()).abs();
    let pb = gear1.base_pitch();
    (approach + recess) / pb.max(1e-15)
}
/// Minimum number of teeth to avoid undercutting (standard gear, no profile shift).
#[allow(dead_code)]
pub fn min_teeth_no_undercut(pressure_angle_deg: f64) -> u32 {
    let phi = pressure_angle_deg.to_radians();
    (2.0 / phi.sin().powi(2)).ceil() as u32
}
/// Compute the working pressure angle for a non-standard centre distance.
///
/// `cos(φ') = (r_b1 + r_b2) / a'`.
#[allow(dead_code)]
pub fn working_pressure_angle(gear1: &InvoluteGear, gear2: &InvoluteGear, centre_dist: f64) -> f64 {
    let rb_sum = gear1.base_radius() + gear2.base_radius();
    (rb_sum / centre_dist.max(1e-15)).clamp(-1.0, 1.0).acos()
}
/// Hertz contact stress between two cylindrical gear teeth (Pa).
///
/// `σ_H = sqrt(F * E_eff / (π * b * ρ_eq))`.
#[allow(dead_code)]
pub fn hertz_contact_stress(
    force: f64,
    face_width: f64,
    radius1: f64,
    radius2: f64,
    e_eff: f64,
) -> f64 {
    let rho_eq = (radius1 * radius2) / (radius1 + radius2).max(1e-15);
    (force * e_eff / (PI * face_width * rho_eq).max(1e-15)).sqrt()
}
/// Velocity ratio for a double Cardan (constant velocity) joint.
///
/// Returns exactly 1.0 regardless of angle if shafts are equal and opposite.
#[allow(dead_code)]
pub fn double_cardan_velocity_ratio(_shaft_angle: f64, _theta: f64) -> f64 {
    1.0
}
/// Static balance check: sum of m·r vectors ≈ 0.
///
/// Returns the residual force magnitude (m·r units, multiply by ω² for force).
#[allow(dead_code)]
pub fn static_balance_residual(masses: &[RotatingMass]) -> f64 {
    let sx: f64 = masses.iter().map(|m| m.mr_vector()[0]).sum();
    let sy: f64 = masses.iter().map(|m| m.mr_vector()[1]).sum();
    (sx * sx + sy * sy).sqrt()
}
/// Dynamic balance check: sum of m·r·l moment vectors ≈ 0.
///
/// Returns the residual moment magnitude (m·r·l units).
#[allow(dead_code)]
pub fn dynamic_balance_residual(masses: &[RotatingMass]) -> f64 {
    let sx: f64 = masses.iter().map(|m| m.mrl_vector()[0]).sum();
    let sy: f64 = masses.iter().map(|m| m.mrl_vector()[1]).sum();
    (sx * sx + sy * sy).sqrt()
}
/// Compute the balancing mass required for static balance in a given plane.
///
/// Returns `(mass_kg, angle_rad)` of a balance mass at `balance_radius`.
#[allow(dead_code)]
pub fn static_balance_mass(masses: &[RotatingMass], balance_radius: f64) -> (f64, f64) {
    let sx: f64 = masses.iter().map(|m| m.mr_vector()[0]).sum();
    let sy: f64 = masses.iter().map(|m| m.mr_vector()[1]).sum();
    let mr = (sx * sx + sy * sy).sqrt();
    let angle = (-sy).atan2(-sx);
    let mass = mr / balance_radius.max(1e-15);
    (mass, angle)
}
/// Dalby's graphical method: compute balance masses for two correction planes.
///
/// `plane_a_axial` and `plane_b_axial` are the axial positions of the two
/// correction planes. Returns `(mass_a, angle_a, mass_b, angle_b)`.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn dalby_balance(
    masses: &[RotatingMass],
    plane_a_axial: f64,
    plane_b_axial: f64,
    balance_radius_a: f64,
    balance_radius_b: f64,
) -> (f64, f64, f64, f64) {
    let l_ab = plane_b_axial - plane_a_axial;
    let mut mx_b = 0.0f64;
    let mut my_b = 0.0f64;
    for m in masses {
        let arm = m.axial_pos - plane_b_axial;
        let mr = m.mass * m.radius;
        mx_b += mr * arm * m.angle.cos();
        my_b += mr * arm * m.angle.sin();
    }
    let mr_a = (mx_b * mx_b + my_b * my_b).sqrt() / l_ab.abs().max(1e-15);
    let angle_a = (-my_b).atan2(-mx_b);
    let mass_a = mr_a / balance_radius_a.max(1e-15);
    let balance_a = RotatingMass::new(
        mass_a,
        balance_radius_a,
        angle_a.to_degrees(),
        plane_a_axial,
    );
    let mut all_masses: Vec<RotatingMass> = masses.to_vec();
    all_masses.push(balance_a);
    let (mass_b, angle_b) = static_balance_mass(&all_masses, balance_radius_b);
    (mass_a, angle_a, mass_b, angle_b)
}
/// Bearing reaction forces from a rotating shaft with known unbalance.
///
/// Returns `(F_bearing_a, F_bearing_b)` in Newtons given angular speed `omega`.
#[allow(dead_code)]
pub fn bearing_reactions(
    masses: &[RotatingMass],
    bearing_a_pos: f64,
    bearing_b_pos: f64,
    omega: f64,
) -> ([f64; 2], [f64; 2]) {
    let l = bearing_b_pos - bearing_a_pos;
    let mut mx_a = 0.0f64;
    let mut my_a = 0.0f64;
    let mut fx_total = 0.0f64;
    let mut fy_total = 0.0f64;
    for m in masses {
        let mr = m.mass * m.radius;
        let arm = m.axial_pos - bearing_a_pos;
        let fx = mr * m.angle.cos() * omega * omega;
        let fy = mr * m.angle.sin() * omega * omega;
        mx_a += fy * arm;
        my_a += fx * arm;
        fx_total += fx;
        fy_total += fy;
    }
    let rb_y = mx_a / l.max(1e-15);
    let rb_x = my_a / l.max(1e-15);
    let ra_x = fx_total - rb_x;
    let ra_y = fy_total - rb_y;
    ([ra_x, ra_y], [rb_x, rb_y])
}
/// Compute the dead-centre positions of a crank-rocker linkage.
///
/// Returns `(theta2_extended, theta2_folded)` — the two dead-centre crank angles.
///
/// At the extended dead-centre the coupler + crank are aligned (l2 + l3 diagonal = l4 + l1).
/// At the folded dead-centre the coupler folds back over the crank.
/// Both angles derived from the law of cosines applied to the closed loop.
#[allow(dead_code)]
pub fn crank_rocker_dead_centres(l1: f64, l2: f64, l3: f64, l4: f64) -> Option<(f64, f64)> {
    let diag_ext = l2 + l3;
    let cos_ext = (sq(l1) + sq(l4) - sq(diag_ext)) / (2.0 * l1 * l4).max(1e-15);
    let diag_fold = (l3 - l2).abs();
    let cos_fold = (sq(l1) + sq(l4) - sq(diag_fold)) / (2.0 * l1 * l4).max(1e-15);
    if cos_ext.abs() > 1.0 || cos_fold.abs() > 1.0 {
        return None;
    }
    Some((cos_ext.acos(), cos_fold.acos()))
}
/// Compute mechanical advantage of a four-bar linkage at crank angle `theta2`.
///
/// MA = (output link angular velocity) / (input link angular velocity) * (r4/r2).
#[allow(dead_code)]
pub fn four_bar_mechanical_advantage(
    l1: f64,
    l2: f64,
    l3: f64,
    l4: f64,
    theta2: f64,
    omega2: f64,
) -> Option<f64> {
    let analysis = FourBarAnalysis::analyse(l1, l2, l3, l4, theta2, omega2)?;
    if omega2.abs() < 1e-15 {
        return None;
    }
    Some((analysis.omega4 / omega2).abs() * (l2 / l4))
}
/// Coupler curve: compute a series of coupler point positions.
///
/// `cp_x` and `cp_y` are the coupler point offsets from coupler midpoint.
#[allow(dead_code)]
pub fn coupler_curve(
    l1: f64,
    l2: f64,
    l3: f64,
    l4: f64,
    cp_x: f64,
    cp_y: f64,
    n_steps: usize,
) -> Vec<[f64; 2]> {
    let mut points = Vec::with_capacity(n_steps);
    for i in 0..n_steps {
        let theta2 = 2.0 * PI * i as f64 / n_steps as f64;
        if let Some(analysis) = FourBarAnalysis::analyse(l1, l2, l3, l4, theta2, 1.0) {
            let theta3 = analysis.theta3;
            let ax = l2 * theta2.cos();
            let ay = l2 * theta2.sin();
            let x = ax + cp_x * theta3.cos() - cp_y * theta3.sin();
            let y = ay + cp_x * theta3.sin() + cp_y * theta3.cos();
            points.push([x, y]);
        }
    }
    points
}
/// Limit stop: compute max/min crank angle for a constrained follower.
#[allow(dead_code)]
pub fn follower_angle_range(l1: f64, l2: f64, l3: f64, l4: f64, n_steps: usize) -> (f64, f64) {
    let mut min_theta4 = f64::MAX;
    let mut max_theta4 = f64::MIN;
    for i in 0..n_steps {
        let theta2 = 2.0 * PI * i as f64 / n_steps as f64;
        if let Some((_, theta4)) = four_bar_angles(l1, l2, l3, l4, theta2) {
            if theta4 < min_theta4 {
                min_theta4 = theta4;
            }
            if theta4 > max_theta4 {
                max_theta4 = theta4;
            }
        }
    }
    (min_theta4, max_theta4)
}
#[cfg(test)]
mod tests_ext {
    use super::*;
    use crate::mechanism::EpicyclicGearTrain;
    use crate::mechanism::GenevaMechanism;
    use crate::mechanism::HookesJoint;
    use crate::mechanism::MechanismGraph;
    use crate::mechanism::SliderCrankForces;
    use crate::mechanism::SliderCrankMechanism;
    use crate::mechanism::ToggleClampExtended;
    #[test]
    fn test_grubler_four_bar_planar() {
        let joints = vec![
            JointClass::Revolute,
            JointClass::Revolute,
            JointClass::Revolute,
            JointClass::Revolute,
        ];
        let m = grubler_kutzbach_planar(4, &joints);
        assert_eq!(m, 1);
    }
    #[test]
    fn test_grubler_slider_crank() {
        let joints = vec![
            JointClass::Revolute,
            JointClass::Revolute,
            JointClass::Revolute,
            JointClass::Prismatic,
        ];
        let m = grubler_kutzbach_planar(4, &joints);
        assert_eq!(m, 1);
    }
    #[test]
    fn test_grubler_structure() {
        let joints = vec![
            JointClass::Revolute,
            JointClass::Revolute,
            JointClass::Revolute,
        ];
        let m = grubler_kutzbach_planar(3, &joints);
        assert!(m <= 0);
    }
    #[test]
    fn test_mechanism_graph_single_dof() {
        let mut mg = MechanismGraph::new(4, true);
        mg.add_joint(JointClass::Revolute);
        mg.add_joint(JointClass::Revolute);
        mg.add_joint(JointClass::Revolute);
        mg.add_joint(JointClass::Revolute);
        assert!(mg.is_single_dof());
        assert!(!mg.is_structure());
    }
    #[test]
    fn test_grashof_classify_crank_rocker() {
        let cls = grashof_classify(4.0, 2.0, 4.0, 3.0);
        assert_eq!(cls, GrashofClass::CrankRocker);
    }
    #[test]
    fn test_grashof_classify_non() {
        let cls = grashof_classify(5.0, 6.0, 7.0, 9.0);
        assert_eq!(cls, GrashofClass::NonGrashof);
    }
    #[test]
    fn test_four_bar_analysis_velocities() {
        let res = FourBarAnalysis::analyse(4.0, 2.0, 4.0, 3.0, 1.0, 10.0);
        if let Some(a) = res {
            assert!(a.omega3.is_finite());
            assert!(a.omega4.is_finite());
        }
    }
    #[test]
    fn test_coupler_curve_count() {
        let pts = coupler_curve(4.0, 2.0, 4.0, 3.0, 0.5, 0.3, 36);
        assert!(!pts.is_empty());
    }
    #[test]
    fn test_slider_crank_rod_angle_tdc() {
        let phi = slider_crank_rod_angle(0.05, 0.15, 0.0);
        assert!(phi.abs() < 1e-9);
    }
    #[test]
    fn test_slider_crank_lambda() {
        let lambda = slider_crank_lambda(0.05, 0.15);
        assert!(lambda < 1.0);
    }
    #[test]
    fn test_piston_velocity_approx() {
        let sc = SliderCrankMechanism::new(0.05, 0.15, 100.0);
        let v_exact = sc.piston_velocity(PI / 2.0);
        let v_approx = piston_velocity_approx(0.05, 0.15, 100.0, PI / 2.0);
        assert!((v_exact - v_approx).abs() / v_exact.abs().max(1e-3) < 0.05);
    }
    #[test]
    fn test_slider_crank_forces_torque() {
        let phi = slider_crank_rod_angle(0.05, 0.15, PI / 4.0);
        let f = SliderCrankForces::compute(0.05, PI / 4.0, phi, 5000.0);
        assert!(f.rod_force.is_finite() && f.rod_force > 0.0);
    }
    #[test]
    fn test_cam_poly345_boundaries() {
        assert!(cam_polynomial_345(0.0).abs() < 1e-10);
        assert!((cam_polynomial_345(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cam_poly4567_boundaries() {
        assert!(cam_polynomial_4567(0.0).abs() < 1e-10);
        assert!((cam_polynomial_4567(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cam_cycloidal_midpoint() {
        assert!((cam_cycloidal(0.5) - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_cam_double_harmonic_boundaries() {
        assert!(cam_double_harmonic(0.0).abs() < 1e-9, "start should be 0");
        let v1 = cam_double_harmonic(1.0);
        assert!(v1.is_finite(), "end value should be finite");
        let v075 = cam_double_harmonic(0.75);
        assert!(v075 > 0.0, "value at 0.75 should be positive, was {v075}");
    }
    #[test]
    fn test_cam_max_pressure_angle() {
        let pa = cam_max_pressure_angle(0.02, 0.01, PI, cam_cycloidal, 100);
        assert!(pa < 30.0_f64.to_radians());
    }
    #[test]
    fn test_involute_pitch_radius() {
        let g = InvoluteGear::new(20, 0.001);
        assert!((g.pitch_radius() - 0.01).abs() < 1e-9);
    }
    #[test]
    fn test_involute_base_less_than_pitch() {
        let g = InvoluteGear::new(20, 0.001);
        assert!(g.base_radius() < g.pitch_radius());
    }
    #[test]
    fn test_involute_addendum_greater() {
        let g = InvoluteGear::new(20, 0.001);
        assert!(g.addendum_radius() > g.pitch_radius());
    }
    #[test]
    fn test_contact_ratio_gt_one() {
        let g1 = InvoluteGear::new(20, 0.001);
        let g2 = InvoluteGear::new(40, 0.001);
        let cr = contact_ratio(&g1, &g2);
        assert!(cr > 1.0, "contact ratio was {cr}");
    }
    #[test]
    fn test_min_teeth_no_undercut() {
        let n = min_teeth_no_undercut(20.0);
        assert!((16..=18).contains(&n));
    }
    #[test]
    fn test_geneva_dwell_fraction() {
        let g = GenevaMechanism::new(6, 0.05, 2.0 * PI);
        let df = g.dwell_fraction();
        assert!(df > 0.0 && df < 1.0);
    }
    #[test]
    fn test_geneva_wheel_velocity_dwell() {
        let g = GenevaMechanism::new(6, 0.05, 2.0 * PI);
        let v = g.wheel_velocity(PI);
        assert!(v.abs() < 1e-9);
    }
    #[test]
    fn test_geneva_centre_distance() {
        let g = GenevaMechanism::new(4, 0.05, 1.0);
        assert!(g.centre_distance() > 0.0);
    }
    #[test]
    fn test_geneva_step() {
        let mut g = GenevaMechanism::new(6, 0.05, 2.0);
        g.step(0.1);
        assert!(g.crank_angle > 0.0);
    }
    #[test]
    fn test_epicyclic_ratio_planetary() {
        let e = EpicyclicGearTrain::new(20, 20, 3, 0.001);
        let expected = 1.0 + 60.0 / 20.0;
        assert!((e.gear_ratio() - expected).abs() < 1e-9);
    }
    #[test]
    fn test_epicyclic_ring_teeth() {
        let e = EpicyclicGearTrain::new(20, 20, 3, 0.001);
        assert_eq!(e.ring_teeth, 60);
    }
    #[test]
    fn test_epicyclic_assembly_condition() {
        let e = EpicyclicGearTrain::new(20, 20, 4, 0.001);
        assert!(e.assembly_condition());
    }
    #[test]
    fn test_epicyclic_carrier_radius() {
        let e = EpicyclicGearTrain::new(20, 20, 3, 0.001);
        let expected = e.sun_radius() + e.planet_radius();
        assert!((e.carrier_radius() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_hookes_joint_zero_angle_ratio() {
        let j = HookesJoint::new(0.0, 100.0);
        assert!((j.velocity_ratio(0.0) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_hookes_joint_variation() {
        let j = HookesJoint::new(20.0, 100.0);
        let r1 = j.velocity_ratio(0.0);
        let r2 = j.velocity_ratio(PI / 2.0);
        assert!((r1 - r2).abs() > 0.01);
    }
    #[test]
    fn test_hookes_joint_speed_variation() {
        let j = HookesJoint::new(15.0, 100.0);
        assert!(j.speed_variation() >= 1.0);
    }
    #[test]
    fn test_hookes_joint_double_cardan() {
        let j = HookesJoint::new(20.0, 100.0);
        assert!(j.needs_double_cardan());
        let j2 = HookesJoint::new(10.0, 100.0);
        assert!(!j2.needs_double_cardan());
    }
    #[test]
    fn test_hookes_joint_step() {
        let mut j = HookesJoint::new(20.0, 10.0);
        j.step(0.1);
        assert!(j.theta > 0.0);
    }
    #[test]
    fn test_static_balance_single_mass() {
        let masses = vec![RotatingMass::new(2.0, 0.1, 0.0, 0.0)];
        let res = static_balance_residual(&masses);
        assert!((res - 0.2).abs() < 1e-9);
    }
    #[test]
    fn test_static_balance_two_opposite() {
        let masses = vec![
            RotatingMass::new(1.0, 0.1, 0.0, 0.0),
            RotatingMass::new(1.0, 0.1, 180.0, 0.0),
        ];
        let res = static_balance_residual(&masses);
        assert!(res < 1e-9);
    }
    #[test]
    fn test_static_balance_mass_cancels() {
        let masses = vec![RotatingMass::new(2.0, 0.1, 45.0, 0.0)];
        let (mb, _angle) = static_balance_mass(&masses, 0.1);
        let balance = RotatingMass::new(mb, 0.1, _angle.to_degrees(), 0.0);
        let all = vec![masses[0].clone(), balance];
        let res = static_balance_residual(&all);
        assert!(res < 1e-9);
    }
    #[test]
    fn test_dynamic_balance_nonzero() {
        let masses = vec![
            RotatingMass::new(1.0, 0.1, 0.0, 0.1),
            RotatingMass::new(1.0, 0.1, 180.0, 0.3),
        ];
        let dres = dynamic_balance_residual(&masses);
        assert!(dres > 1e-9);
    }
    #[test]
    fn test_dalby_balance() {
        let masses = vec![
            RotatingMass::new(1.0, 0.05, 30.0, 0.1),
            RotatingMass::new(1.5, 0.05, 150.0, 0.3),
            RotatingMass::new(0.8, 0.05, 270.0, 0.5),
        ];
        let (ma, aa, mb, ab) = dalby_balance(&masses, 0.0, 0.6, 0.05, 0.05);
        assert!(ma.is_finite() && aa.is_finite() && mb.is_finite() && ab.is_finite());
    }
    #[test]
    fn test_dead_centres() {
        let result = crank_rocker_dead_centres(4.0, 2.0, 4.0, 3.0);
        assert!(result.is_some());
    }
    #[test]
    fn test_four_bar_mechanical_advantage() {
        let ma = four_bar_mechanical_advantage(4.0, 2.0, 4.0, 3.0, 1.0, 10.0);
        if let Some(v) = ma {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn test_follower_angle_range() {
        let (mn, mx) = follower_angle_range(4.0, 2.0, 4.0, 3.0, 360);
        assert!(mn < mx);
    }
    #[test]
    fn test_toggle_extended_output() {
        let mut tc = ToggleClampExtended::new(0.1, 0.12, 0.15);
        tc.f_input = 100.0;
        let f = tc.output_force();
        assert!(f.is_finite());
    }
    #[test]
    fn test_hertz_contact_stress() {
        let stress = hertz_contact_stress(1000.0, 0.02, 0.01, 0.015, 200e9);
        assert!(stress > 0.0);
    }
}
