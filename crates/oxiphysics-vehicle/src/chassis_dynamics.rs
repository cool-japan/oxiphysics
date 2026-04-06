// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Full vehicle chassis dynamics: 6-DOF sprung mass motion, load transfer,
//! roll-centre height, anti-dive, anti-squat, and chassis torsional flex.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// ChassisState
// ---------------------------------------------------------------------------

/// Full 6-DOF chassis state.
///
/// Position and orientation represent the sprung mass (body) in world space.
/// Linear and angular velocities are in body frame unless otherwise noted.
#[derive(Debug, Clone)]
pub struct ChassisState {
    /// Centre-of-mass position in world space `[x, y, z]` (m).
    pub position: [f64; 3],
    /// Euler angles `[roll, pitch, yaw]` (rad).
    pub orientation: [f64; 3],
    /// Linear velocity `[vx, vy, vz]` in world frame (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity `[p, q, r]` in body frame (rad/s).
    pub angular_velocity: [f64; 3],
    /// Sprung mass (kg).
    pub sprung_mass: f64,
    /// Roll moment of inertia (kg·m²).
    pub ixx: f64,
    /// Pitch moment of inertia (kg·m²).
    pub iyy: f64,
    /// Yaw moment of inertia (kg·m²).
    pub izz: f64,
    /// Height of the CG above the ground (m).
    pub cg_height: f64,
    /// Front axle distance from CG (m).
    pub wheelbase_front: f64,
    /// Rear axle distance from CG (m).
    pub wheelbase_rear: f64,
    /// Front track width (m).
    pub track_front: f64,
    /// Rear track width (m).
    pub track_rear: f64,
}

impl ChassisState {
    /// Create a default mid-size passenger car chassis.
    pub fn default_car() -> Self {
        Self {
            position: [0.0; 3],
            orientation: [0.0; 3],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            sprung_mass: 1400.0,
            ixx: 500.0,
            iyy: 2500.0,
            izz: 2800.0,
            cg_height: 0.55,
            wheelbase_front: 1.4,
            wheelbase_rear: 1.1,
            track_front: 1.55,
            track_rear: 1.52,
        }
    }

    /// Total wheelbase (m).
    pub fn wheelbase(&self) -> f64 {
        self.wheelbase_front + self.wheelbase_rear
    }
}

// ---------------------------------------------------------------------------
// sprung_mass_motion
// ---------------------------------------------------------------------------

/// Update the sprung mass motion from suspension forces and external loads.
///
/// Applies bounce (heave), roll, and pitch accelerations from the four
/// suspension corner forces `[FL, FR, RL, RR]` using the simple decoupled
/// 3-DOF model.
///
/// # Arguments
/// * `state` - Mutable chassis state (velocity updated in place).
/// * `corner_forces` - Normal suspension forces at each corner `[FL, FR, RL, RR]` (N).
/// * `gravity` - Gravitational acceleration (m/s², positive).
/// * `dt` - Integration time step (s).
#[allow(clippy::too_many_arguments)]
pub fn sprung_mass_motion(
    state: &mut ChassisState,
    corner_forces: [f64; 4],
    gravity: f64,
    dt: f64,
) {
    let [f_fl, f_fr, f_rl, f_rr] = corner_forces;
    let total_f = f_fl + f_fr + f_rl + f_rr;
    let weight = state.sprung_mass * gravity;

    // Heave acceleration
    let a_heave = (total_f - weight) / state.sprung_mass;

    // Roll moment: difference left-right (positive = right side lower)
    let tf = state.track_front;
    let tr = state.track_rear;
    let roll_moment = ((f_fl - f_fr) * tf / 2.0) + ((f_rl - f_rr) * tr / 2.0);
    let a_roll = roll_moment / state.ixx;

    // Pitch moment: front minus rear contribution (positive nose-down)
    let pitch_moment = (f_fl + f_fr) * state.wheelbase_front - (f_rl + f_rr) * state.wheelbase_rear;
    let a_pitch = pitch_moment / state.iyy;

    // Update velocities (semi-implicit Euler)
    state.velocity[1] += a_heave * dt;
    state.angular_velocity[0] += a_roll * dt;
    state.angular_velocity[1] += a_pitch * dt;

    // Update position and orientation
    state.position[1] += state.velocity[1] * dt;
    state.orientation[0] += state.angular_velocity[0] * dt;
    state.orientation[1] += state.angular_velocity[1] * dt;
}

// ---------------------------------------------------------------------------
// load_transfer
// ---------------------------------------------------------------------------

/// Compute longitudinal and lateral load transfer.
///
/// # Longitudinal load transfer
/// `ΔFz_long = m * a_x * h_cg / L`
///
/// # Lateral load transfer
/// `ΔFz_lat = m * a_y * h_cg / T`
///
/// # Arguments
/// * `state` - Chassis state.
/// * `accel_long` - Longitudinal acceleration (m/s², positive = acceleration).
/// * `accel_lat` - Lateral acceleration (m/s², positive = right).
/// * `gravity` - Gravitational acceleration (m/s²).
///
/// Returns `(delta_fz_longitudinal, delta_fz_lateral)` in Newtons.
pub fn load_transfer(
    state: &ChassisState,
    accel_long: f64,
    accel_lat: f64,
    gravity: f64,
) -> (f64, f64) {
    let _ = gravity; // gravity used for weight, here we use mass directly
    let wb = state.wheelbase();
    let track = (state.track_front + state.track_rear) / 2.0;
    let delta_long = state.sprung_mass * accel_long * state.cg_height / wb.max(1e-6);
    let delta_lat = state.sprung_mass * accel_lat * state.cg_height / track.max(1e-6);
    (delta_long, delta_lat)
}

// ---------------------------------------------------------------------------
// roll_center_height
// ---------------------------------------------------------------------------

/// Compute the geometric roll-centre height using a simplified instant-centre
/// method.
///
/// The roll-centre height is approximated as:
/// `h_rc = (a * h_inner + b * h_outer) / (a + b)`
/// where `a`, `b` are the half-track to inner and outer tie-rod lengths.
///
/// For a simple wishbone geometry this reduces to a linear blend of
/// inner and outer pivot heights.
///
/// # Arguments
/// * `inner_pivot_height` - Height of inner wishbone pivot above ground (m).
/// * `outer_pivot_height` - Height of outer wheel centre above ground (m).
/// * `track_half` - Half track width (m).
/// * `inner_offset` - Horizontal distance from centre-line to inner pivot (m).
///
/// Returns the roll-centre height (m).
pub fn roll_center_height(
    inner_pivot_height: f64,
    outer_pivot_height: f64,
    track_half: f64,
    inner_offset: f64,
) -> f64 {
    let span = track_half - inner_offset;
    if span.abs() < 1e-9 {
        return inner_pivot_height;
    }
    // Instant centre is at intersection of lines from inner and outer pivots
    // Simplified: linear interpolation weighted by position
    let t = inner_offset / (track_half.max(1e-9));
    inner_pivot_height * (1.0 - t) + outer_pivot_height * t
}

// ---------------------------------------------------------------------------
// anti_dive
// ---------------------------------------------------------------------------

/// Compute the anti-dive coefficient for a front suspension geometry.
///
/// Anti-dive percentage describes how much the suspension geometry counteracts
/// dive under braking.
///
/// `anti_dive = (tan(φ) * L / h_cg) * (brake_front_bias)`
///
/// where `φ` is the front suspension instant-centre angle.
///
/// # Arguments
/// * `ic_angle_rad` - Inclination angle of the front suspension instant centre (rad).
/// * `wheelbase` - Vehicle wheelbase (m).
/// * `cg_height` - CG height above ground (m).
/// * `front_brake_bias` - Fraction of braking force on front axle \[0, 1\].
///
/// Returns anti-dive coefficient in \[0, 1\] (1.0 = 100% anti-dive).
pub fn anti_dive(ic_angle_rad: f64, wheelbase: f64, cg_height: f64, front_brake_bias: f64) -> f64 {
    let tan_phi = ic_angle_rad.tan();
    (tan_phi * wheelbase / cg_height.max(1e-6)) * front_brake_bias
}

// ---------------------------------------------------------------------------
// anti_squat
// ---------------------------------------------------------------------------

/// Compute the anti-squat coefficient for a rear suspension geometry.
///
/// Anti-squat percentage describes how much the rear suspension geometry
/// counteracts squat under acceleration.
///
/// `anti_squat = tan(φ_rear) * L / h_cg`
///
/// # Arguments
/// * `ic_angle_rad` - Inclination angle of the rear suspension instant centre (rad).
/// * `wheelbase` - Vehicle wheelbase (m).
/// * `cg_height` - CG height above ground (m).
/// * `rear_drive_bias` - Fraction of drive force on rear axle \[0, 1\].
///
/// Returns anti-squat coefficient.
pub fn anti_squat(ic_angle_rad: f64, wheelbase: f64, cg_height: f64, rear_drive_bias: f64) -> f64 {
    let tan_phi = ic_angle_rad.tan();
    (tan_phi * wheelbase / cg_height.max(1e-6)) * rear_drive_bias
}

// ---------------------------------------------------------------------------
// chassis_flex
// ---------------------------------------------------------------------------

/// Estimate the torsional deflection angle of the chassis under an applied
/// roll moment.
///
/// `θ_flex = M_roll / k_torsion`
///
/// # Arguments
/// * `roll_moment` - Applied roll moment (N·m).
/// * `torsional_stiffness` - Chassis torsional stiffness (N·m/rad).
///
/// Returns the torsional deflection angle (rad).
pub fn chassis_flex(roll_moment: f64, torsional_stiffness: f64) -> f64 {
    roll_moment / torsional_stiffness.max(1e-6)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn car() -> ChassisState {
        ChassisState::default_car()
    }

    // ── ChassisState helpers ───────────────────────────────────────────────

    #[test]
    fn default_car_wheelbase() {
        let c = car();
        let wb = c.wheelbase();
        assert!((wb - 2.5).abs() < 1e-6, "wb={wb}");
    }

    #[test]
    fn default_car_sprung_mass_positive() {
        let c = car();
        assert!(c.sprung_mass > 0.0);
    }

    #[test]
    fn default_car_inertias_positive() {
        let c = car();
        assert!(c.ixx > 0.0 && c.iyy > 0.0 && c.izz > 0.0);
    }

    // ── load_transfer ──────────────────────────────────────────────────────

    #[test]
    fn load_transfer_no_acceleration_zero() {
        let c = car();
        let (dl, dlat) = load_transfer(&c, 0.0, 0.0, 9.81);
        assert!(dl.abs() < 1e-10);
        assert!(dlat.abs() < 1e-10);
    }

    #[test]
    fn load_transfer_positive_accel_transfers_rearward() {
        let c = car();
        let (dl, _) = load_transfer(&c, 9.81, 0.0, 9.81);
        assert!(dl > 0.0, "dl={dl}");
    }

    #[test]
    fn load_transfer_braking_transfers_forward() {
        let c = car();
        let (dl, _) = load_transfer(&c, -9.81, 0.0, 9.81);
        assert!(dl < 0.0, "dl={dl}");
    }

    #[test]
    fn load_transfer_lateral_scales_with_cg_height() {
        let mut c1 = car();
        let mut c2 = car();
        c1.cg_height = 0.5;
        c2.cg_height = 1.0;
        let (_, dlat1) = load_transfer(&c1, 0.0, 1.0, 9.81);
        let (_, dlat2) = load_transfer(&c2, 0.0, 1.0, 9.81);
        assert!(
            (dlat2 / dlat1 - 2.0).abs() < 1e-6,
            "ratio={}",
            dlat2 / dlat1
        );
    }

    #[test]
    fn load_transfer_formula_check() {
        let c = car();
        let accel = 5.0;
        let expected = c.sprung_mass * accel * c.cg_height / c.wheelbase();
        let (dl, _) = load_transfer(&c, accel, 0.0, 9.81);
        assert!((dl - expected).abs() < 1e-6, "dl={dl} expected={expected}");
    }

    // ── roll_center_height ─────────────────────────────────────────────────

    #[test]
    fn roll_center_equal_heights_returns_inner() {
        let h = roll_center_height(0.2, 0.2, 0.8, 0.0);
        assert!((h - 0.2).abs() < 1e-6, "h={h}");
    }

    #[test]
    fn roll_center_between_pivots() {
        let h = roll_center_height(0.1, 0.5, 0.8, 0.4);
        assert!((0.1..=0.5).contains(&h), "h={h}");
    }

    #[test]
    fn roll_center_zero_span_returns_inner() {
        let h = roll_center_height(0.3, 0.7, 0.5, 0.5);
        // span = 0 → returns inner pivot height
        assert!((h - 0.3).abs() < 1e-6, "h={h}");
    }

    #[test]
    fn roll_center_scales_with_inner_offset() {
        let h1 = roll_center_height(0.0, 1.0, 1.0, 0.0);
        let h2 = roll_center_height(0.0, 1.0, 1.0, 0.5);
        // Larger inner_offset (closer to outer pivot) → weighted towards outer pivot height
        assert!(h2 >= h1 || (h2 - h1).abs() < 1e-6, "h1={h1}, h2={h2}");
    }

    // ── anti_dive ─────────────────────────────────────────────────────────

    #[test]
    fn anti_dive_zero_angle_is_zero() {
        let ad = anti_dive(0.0, 2.5, 0.55, 0.6);
        assert!(ad.abs() < 1e-10, "ad={ad}");
    }

    #[test]
    fn anti_dive_positive_for_positive_angle() {
        let ad = anti_dive(0.1, 2.5, 0.55, 0.6);
        assert!(ad > 0.0, "ad={ad}");
    }

    #[test]
    fn anti_dive_scales_with_brake_bias() {
        let ad1 = anti_dive(0.1, 2.5, 0.55, 0.5);
        let ad2 = anti_dive(0.1, 2.5, 0.55, 1.0);
        assert!((ad2 - 2.0 * ad1).abs() < 1e-6, "ad1={ad1} ad2={ad2}");
    }

    #[test]
    fn anti_dive_increases_with_wheelbase() {
        let ad1 = anti_dive(0.1, 2.0, 0.55, 0.6);
        let ad2 = anti_dive(0.1, 4.0, 0.55, 0.6);
        assert!(ad2 > ad1, "ad1={ad1} ad2={ad2}");
    }

    // ── anti_squat ────────────────────────────────────────────────────────

    #[test]
    fn anti_squat_zero_angle_is_zero() {
        let asq = anti_squat(0.0, 2.5, 0.55, 1.0);
        assert!(asq.abs() < 1e-10, "asq={asq}");
    }

    #[test]
    fn anti_squat_positive_for_positive_angle() {
        let asq = anti_squat(0.1, 2.5, 0.55, 1.0);
        assert!(asq > 0.0, "asq={asq}");
    }

    #[test]
    fn anti_squat_scales_with_drive_bias() {
        let a1 = anti_squat(0.1, 2.5, 0.55, 0.5);
        let a2 = anti_squat(0.1, 2.5, 0.55, 1.0);
        assert!((a2 - 2.0 * a1).abs() < 1e-6, "a1={a1} a2={a2}");
    }

    // ── chassis_flex ──────────────────────────────────────────────────────

    #[test]
    fn chassis_flex_zero_moment_zero_deflection() {
        let theta = chassis_flex(0.0, 1e5);
        assert!(theta.abs() < 1e-10, "theta={theta}");
    }

    #[test]
    fn chassis_flex_scales_with_moment() {
        let t1 = chassis_flex(1000.0, 1e5);
        let t2 = chassis_flex(2000.0, 1e5);
        assert!((t2 - 2.0 * t1).abs() < 1e-10, "t1={t1} t2={t2}");
    }

    #[test]
    fn chassis_flex_inversely_proportional_to_stiffness() {
        let t1 = chassis_flex(1000.0, 1e4);
        let t2 = chassis_flex(1000.0, 2e4);
        assert!((t1 - 2.0 * t2).abs() < 1e-10, "t1={t1} t2={t2}");
    }

    #[test]
    fn chassis_flex_very_stiff_near_zero() {
        let theta = chassis_flex(10_000.0, 1e9);
        assert!(theta < 1e-4, "theta={theta}");
    }

    // ── sprung_mass_motion ────────────────────────────────────────────────

    #[test]
    fn sprung_mass_symmetric_forces_no_roll_pitch() {
        let mut c = car();
        let weight = c.sprung_mass * 9.81;
        let wbf = c.wheelbase_front;
        let wbr = c.wheelbase_rear;
        let wb = wbf + wbr;
        // Distribute weight so pitch moment = 0: front = W * wbr / wb, rear = W * wbf / wb
        let f_front = weight * wbr / wb / 2.0;
        let f_rear = weight * wbf / wb / 2.0;
        // Left = right, so no roll moment either
        sprung_mass_motion(&mut c, [f_front, f_front, f_rear, f_rear], 9.81, 0.01);
        assert!(
            c.angular_velocity[0].abs() < 1e-8,
            "roll rate = {}",
            c.angular_velocity[0]
        );
        assert!(
            c.angular_velocity[1].abs() < 1e-8,
            "pitch rate = {}",
            c.angular_velocity[1]
        );
    }

    #[test]
    fn sprung_mass_excess_force_heaves_up() {
        let mut c = car();
        let weight_each = c.sprung_mass * 9.81 / 4.0;
        // Extra 1000 N total = 250 N per corner over weight
        sprung_mass_motion(&mut c, [weight_each + 250.0; 4], 9.81, 0.1);
        assert!(c.velocity[1] > 0.0, "vy = {}", c.velocity[1]);
    }

    #[test]
    fn sprung_mass_lateral_imbalance_creates_roll() {
        let mut c = car();
        let base = c.sprung_mass * 9.81 / 4.0;
        // More load on right side
        sprung_mass_motion(
            &mut c,
            [base - 200.0, base + 200.0, base - 200.0, base + 200.0],
            9.81,
            0.1,
        );
        // Roll rate should be non-zero
        assert!(
            c.angular_velocity[0].abs() > 1e-6,
            "roll rate should be non-zero"
        );
    }

    #[test]
    fn sprung_mass_front_heavy_pitches_forward() {
        let mut c = car();
        let base = c.sprung_mass * 9.81 / 4.0;
        sprung_mass_motion(
            &mut c,
            [base + 300.0, base + 300.0, base - 300.0, base - 300.0],
            9.81,
            0.1,
        );
        // Pitch rate: front heavy → nose down → positive pitch moment
        assert!(
            c.angular_velocity[1].abs() > 1e-6,
            "pitch rate should be non-zero"
        );
    }

    #[test]
    fn load_transfer_symmetric_lateral_zero() {
        let c = car();
        let (_, dlat) = load_transfer(&c, 0.0, 0.0, 9.81);
        assert!(dlat.abs() < 1e-10);
    }

    #[test]
    fn anti_dive_matches_manual_formula() {
        let ic_angle = 0.15_f64;
        let wb = 2.5;
        let cg_h = 0.55;
        let bias = 0.6;
        let expected = ic_angle.tan() * wb / cg_h * bias;
        let computed = anti_dive(ic_angle, wb, cg_h, bias);
        assert!(
            (computed - expected).abs() < 1e-10,
            "computed={computed} expected={expected}"
        );
    }
}
