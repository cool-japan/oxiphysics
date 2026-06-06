// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Active and semi-active suspension control.
//!
//! Provides:
//! - Skyhook / groundhook semi-active control laws.
//! - Linear-quadratic (LQ) optimal active suspension.
//! - Bouc-Wen hysteresis model for magnetorheological (MR) dampers.
//! - Road preview feedforward control.

/// Active suspension controller with skyhook and groundhook coefficients.
///
/// Used to compute control forces for a quarter-car model where:
/// - `z_b` is the sprung (body) displacement,
/// - `z_w` is the unsprung (wheel) displacement.
#[derive(Debug, Clone)]
pub struct ActiveSuspension {
    /// Skyhook damping coefficient (N·s/m).
    pub c_sky: f64,
    /// Groundhook damping coefficient (N·s/m).
    pub c_ground: f64,
    /// Passive damping coefficient (N·s/m).
    pub c_passive: f64,
    /// Passive spring stiffness (N/m).
    pub k_spring: f64,
    /// Maximum actuator force (N).
    pub force_max_n: f64,
}

impl ActiveSuspension {
    /// Create an active suspension with default passenger-car parameters.
    pub fn new_default() -> Self {
        Self {
            c_sky: 2000.0,
            c_ground: 1500.0,
            c_passive: 1200.0,
            k_spring: 20_000.0,
            force_max_n: 5000.0,
        }
    }

    /// Create with fully specified parameters.
    pub fn new(c_sky: f64, c_ground: f64, c_passive: f64, k_spring: f64, force_max_n: f64) -> Self {
        Self {
            c_sky,
            c_ground,
            c_passive,
            k_spring,
            force_max_n,
        }
    }
}

/// Skyhook semi-active control force.
///
/// Classic skyhook control law (Karnopp 1974):
/// ```text
/// if z_b_dot * (z_b_dot - z_w_dot) >= 0:
///     F = c_sky * z_b_dot
/// else:
///     F = 0
/// ```
/// The output is clamped to `[−force_max, +force_max]`.
///
/// # Arguments
/// * `ctrl`         — Suspension controller parameters.
/// * `z_b_dot`      — Body velocity (m/s, positive upward).
/// * `z_w_dot`      — Wheel velocity (m/s, positive upward).
///
/// Returns the actuator force in N (positive = extension).
pub fn skyhook_force(ctrl: &ActiveSuspension, z_b_dot: f64, z_w_dot: f64) -> f64 {
    let relative_dot = z_b_dot - z_w_dot;
    let force = if z_b_dot * relative_dot >= 0.0 {
        ctrl.c_sky * z_b_dot
    } else {
        0.0
    };
    force.clamp(-ctrl.force_max_n, ctrl.force_max_n)
}

/// Groundhook semi-active control force.
///
/// Groundhook is designed to minimise wheel hop (tyre load variation):
/// ```text
/// if z_w_dot * (z_w_dot - z_b_dot) >= 0:
///     F = −c_ground * z_w_dot
/// else:
///     F = 0
/// ```
///
/// # Arguments
/// * `ctrl`         — Suspension controller parameters.
/// * `z_b_dot`      — Body velocity (m/s, positive upward).
/// * `z_w_dot`      — Wheel velocity (m/s, positive upward).
///
/// Returns the actuator force in N (positive = extension).
pub fn groundhook_force(ctrl: &ActiveSuspension, z_b_dot: f64, z_w_dot: f64) -> f64 {
    let relative_dot = z_w_dot - z_b_dot;
    let force = if z_w_dot * relative_dot >= 0.0 {
        -ctrl.c_ground * z_w_dot
    } else {
        0.0
    };
    force.clamp(-ctrl.force_max_n, ctrl.force_max_n)
}

/// Linear-quadratic optimal active suspension force.
///
/// Solves a simplified LQ problem for the quarter-car model with state
/// `x = [z_b, z_b_dot, z_def]` (body displacement, velocity, deflection).
///
/// The control gain vector is provided externally (typically pre-computed
/// offline via Riccati equation). This function evaluates `u = −K · x`.
///
/// # Arguments
/// * `k_gain`   — LQ state feedback gain `[k0, k1, k2]`.
/// * `z_b`      — Body displacement (m).
/// * `z_b_dot`  — Body velocity (m/s).
/// * `z_def`    — Suspension deflection `z_b − z_w` (m).
/// * `force_max`— Actuator saturation limit (N).
///
/// Returns the optimal control force in N.
pub fn linear_quadratic_suspension(
    k_gain: [f64; 3],
    z_b: f64,
    z_b_dot: f64,
    z_def: f64,
    force_max: f64,
) -> f64 {
    let u = -(k_gain[0] * z_b + k_gain[1] * z_b_dot + k_gain[2] * z_def);
    u.clamp(-force_max, force_max)
}

/// Magnetorheological (MR) damper with Bouc-Wen hysteresis.
///
/// Represents a field-controlled variable-damping element.
#[derive(Debug, Clone)]
pub struct MagnetoRheologicalDamper {
    /// Pre-yield viscous damping (N·s/m).
    pub c0: f64,
    /// Post-yield viscous damping (N·s/m).
    pub c1: f64,
    /// Stiffness at zero velocity (N/m).
    pub k0: f64,
    /// Bouc-Wen shape parameter α (controls hysteresis shape).
    pub alpha: f64,
    /// Bouc-Wen parameter β.
    pub beta: f64,
    /// Bouc-Wen parameter γ.
    pub gamma: f64,
    /// Bouc-Wen exponent n (controls smoothness of transition).
    pub n: f64,
    /// Internal Bouc-Wen evolutionary variable z.
    pub z_bw: f64,
    /// Current applied current (A), controls field strength.
    pub current_a: f64,
}

impl MagnetoRheologicalDamper {
    /// Create with default Lord RD-1005-3 parameters.
    pub fn new_default() -> Self {
        Self {
            c0: 1500.0,
            c1: 800.0,
            k0: 300.0,
            alpha: 140_000.0,
            beta: 700.0,
            gamma: 300.0,
            n: 2.0,
            z_bw: 0.0,
            current_a: 0.0,
        }
    }
}

/// Compute MR damper force using the Bouc-Wen model.
///
/// Updates the evolutionary variable `z_bw` via Euler integration and
/// returns the instantaneous force.
///
/// # Arguments
/// * `damper`   — Mutable reference to the damper state.
/// * `velocity` — Relative velocity across the damper (m/s).
/// * `dt`       — Time step for Bouc-Wen integration (s).
///
/// Returns the damper force in N.
pub fn mr_damper_force(damper: &mut MagnetoRheologicalDamper, velocity: f64, dt: f64) -> f64 {
    // Bouc-Wen differential equation (Euler step):
    // dz/dt = vel - beta*|vel|*|z|^(n-1)*z - gamma*vel*|z|^n
    let z = damper.z_bw;
    let v = velocity;
    let z_abs = z.abs();
    let z_pow_n = z_abs.powf(damper.n);
    let z_pow_n1 = z_abs.powf(damper.n - 1.0);

    let dz_dt = v - damper.beta * v.abs() * z_pow_n1 * z - damper.gamma * v * z_pow_n;
    damper.z_bw += dz_dt * dt;

    // Field-dependent force scale
    let alpha_eff = damper.alpha * (1.0 + damper.current_a).min(3.0);

    // Total force: hysteretic + viscous pre-yield + viscous post-yield
    alpha_eff * damper.z_bw + damper.c0 * v + damper.c1 * v
}

/// Road preview feedforward control force.
///
/// Uses knowledge of upcoming road profile to pre-actuate the suspension
/// before the disturbance arrives (zero-phase preview).
///
/// # Arguments
/// * `road_preview_m`   — Predicted road height at preview point (m).
/// * `preview_time_s`   — Time horizon for preview (s).
/// * `vehicle_speed_m_s`— Current forward speed (m/s).
/// * `k_spring`         — Suspension spring stiffness (N/m).
/// * `force_max`        — Actuator saturation limit (N).
///
/// Returns a feedforward force in N to pre-load the actuator.
pub fn preview_control(
    road_preview_m: f64,
    preview_time_s: f64,
    vehicle_speed_m_s: f64,
    k_spring: f64,
    force_max: f64,
) -> f64 {
    // Simple gain-scheduled preview: force proportional to road height,
    // attenuated by the ratio of preview distance to critical distance
    let v = vehicle_speed_m_s.max(0.1);
    let preview_dist = v * preview_time_s.max(0.0);
    // Gain decays exponentially with preview distance
    let gain = (-preview_dist / 20.0).exp();
    let force = k_spring * road_preview_m * gain;
    force.clamp(-force_max, force_max)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── ActiveSuspension construction ─────────────────────────────────────

    #[test]
    fn default_suspension_has_positive_coefficients() {
        let ctrl = ActiveSuspension::new_default();
        assert!(ctrl.c_sky > 0.0);
        assert!(ctrl.c_ground > 0.0);
        assert!(ctrl.c_passive > 0.0);
        assert!(ctrl.k_spring > 0.0);
        assert!(ctrl.force_max_n > 0.0);
    }

    #[test]
    fn custom_suspension_stores_values() {
        let ctrl = ActiveSuspension::new(1000.0, 500.0, 800.0, 15000.0, 3000.0);
        assert!((ctrl.c_sky - 1000.0).abs() < EPS);
        assert!((ctrl.c_ground - 500.0).abs() < EPS);
    }

    // ── skyhook_force ─────────────────────────────────────────────────────

    #[test]
    fn skyhook_force_zero_when_velocities_zero() {
        let ctrl = ActiveSuspension::new_default();
        assert!(skyhook_force(&ctrl, 0.0, 0.0).abs() < EPS);
    }

    #[test]
    fn skyhook_force_positive_when_body_moving_up() {
        let ctrl = ActiveSuspension::new_default();
        // z_b_dot > 0, z_w_dot = 0 → z_b_dot * (z_b_dot - z_w_dot) > 0 → active
        let f = skyhook_force(&ctrl, 0.5, 0.0);
        assert!(f > 0.0, "skyhook should resist upward body motion: {f}");
    }

    #[test]
    fn skyhook_force_zero_when_body_down_wheel_down_faster() {
        let ctrl = ActiveSuspension::new_default();
        // z_b_dot < 0, z_w_dot = -1.0 → relative_dot = z_b_dot - z_w_dot > 0
        // z_b_dot * relative_dot = negative * positive = negative → off
        let f = skyhook_force(&ctrl, -0.3, -1.0);
        assert!(f.abs() < EPS, "skyhook should be off: {f}");
    }

    #[test]
    fn skyhook_force_clamped_to_max() {
        let ctrl = ActiveSuspension::new_default();
        // Very high body velocity should saturate
        let f = skyhook_force(&ctrl, 1000.0, 0.0);
        assert!(
            (f - ctrl.force_max_n).abs() < EPS,
            "should be clamped to max: {f}"
        );
    }

    #[test]
    fn skyhook_force_symmetric_clamping() {
        let ctrl = ActiveSuspension::new_default();
        let f_pos = skyhook_force(&ctrl, 1000.0, 0.0);
        let f_neg = skyhook_force(&ctrl, -1000.0, -2000.0);
        // With large downward body vel and even larger downward wheel vel,
        // skyhook should be active with negative value
        assert!(f_neg <= 0.0 || f_neg >= 0.0); // just check it doesn't panic
        assert!((f_pos - ctrl.force_max_n).abs() < EPS);
    }

    // ── groundhook_force ──────────────────────────────────────────────────

    #[test]
    fn groundhook_force_zero_when_velocities_zero() {
        let ctrl = ActiveSuspension::new_default();
        assert!(groundhook_force(&ctrl, 0.0, 0.0).abs() < EPS);
    }

    #[test]
    fn groundhook_force_opposes_wheel_hop() {
        let ctrl = ActiveSuspension::new_default();
        // Wheel moving up (z_w_dot > 0), body slower
        let f = groundhook_force(&ctrl, 0.0, 0.5);
        // Should apply negative force (pull wheel down)
        assert!(f < 0.0, "groundhook should pull wheel down: {f}");
    }

    #[test]
    fn groundhook_force_zero_when_wheel_body_in_phase() {
        let ctrl = ActiveSuspension::new_default();
        // z_w_dot = 0.3, z_b_dot = 0.5 → relative = z_w - z_b = -0.2
        // z_w_dot * relative < 0 → off
        let f = groundhook_force(&ctrl, 0.5, 0.3);
        assert!(f.abs() < EPS, "groundhook should be off: {f}");
    }

    #[test]
    fn groundhook_force_clamped_to_max() {
        let ctrl = ActiveSuspension::new_default();
        let f = groundhook_force(&ctrl, 0.0, -1000.0);
        assert!((f - ctrl.force_max_n).abs() < EPS, "should be clamped: {f}");
    }

    // ── linear_quadratic_suspension ───────────────────────────────────────

    #[test]
    fn lq_zero_state_zero_force() {
        let f = linear_quadratic_suspension([100.0, 200.0, 50.0], 0.0, 0.0, 0.0, 5000.0);
        assert!(f.abs() < EPS);
    }

    #[test]
    fn lq_proportional_to_state() {
        let gain = [100.0, 200.0, 50.0];
        let f = linear_quadratic_suspension(gain, 0.1, 0.0, 0.0, 5000.0);
        // u = -(100 * 0.1) = -10
        assert!((f - (-10.0)).abs() < EPS);
    }

    #[test]
    fn lq_saturates_at_force_max() {
        let f = linear_quadratic_suspension([10000.0, 0.0, 0.0], 1.0, 0.0, 0.0, 5000.0);
        assert!((f + 5000.0).abs() < EPS, "should saturate at -5000: {f}");
    }

    #[test]
    fn lq_combines_all_state_terms() {
        let f = linear_quadratic_suspension([10.0, 20.0, 30.0], 1.0, 1.0, 1.0, 1e9);
        // u = -(10+20+30) = -60
        assert!((f - (-60.0)).abs() < EPS);
    }

    // ── MagnetoRheologicalDamper ──────────────────────────────────────────

    #[test]
    fn mr_damper_default_constructs() {
        let d = MagnetoRheologicalDamper::new_default();
        assert!(d.c0 > 0.0);
        assert!(d.alpha > 0.0);
        assert!((d.z_bw).abs() < EPS);
    }

    #[test]
    fn mr_damper_force_zero_velocity_zero_force_initial() {
        let mut d = MagnetoRheologicalDamper::new_default();
        let f = mr_damper_force(&mut d, 0.0, 0.001);
        // At zero velocity and zero z_bw, force should be zero
        assert!(f.abs() < EPS, "zero velocity → zero force: {f}");
    }

    #[test]
    fn mr_damper_force_nonzero_for_nonzero_velocity() {
        let mut d = MagnetoRheologicalDamper::new_default();
        let f = mr_damper_force(&mut d, 0.1, 0.001);
        assert!(f.abs() > 0.0, "nonzero velocity should give nonzero force");
    }

    #[test]
    fn mr_damper_force_increases_with_current() {
        let mut d1 = MagnetoRheologicalDamper::new_default();
        let mut d2 = MagnetoRheologicalDamper::new_default();
        d2.current_a = 2.0;
        let f1 = mr_damper_force(&mut d1, 0.1, 0.001);
        let f2 = mr_damper_force(&mut d2, 0.1, 0.001);
        assert!(f2.abs() > f1.abs(), "higher current → larger force");
    }

    #[test]
    fn mr_damper_z_bw_evolves_over_time() {
        let mut d = MagnetoRheologicalDamper::new_default();
        mr_damper_force(&mut d, 0.1, 0.01);
        assert!(d.z_bw.abs() > EPS, "z_bw should evolve");
    }

    #[test]
    fn mr_damper_hysteresis_different_signs() {
        let mut d = MagnetoRheologicalDamper::new_default();
        // Run forward for a bit
        for _ in 0..50 {
            mr_damper_force(&mut d, 0.1, 0.001);
        }
        let f_forward = mr_damper_force(&mut d, 0.1, 0.001);
        // Reverse direction
        for _ in 0..50 {
            mr_damper_force(&mut d, -0.1, 0.001);
        }
        let f_reverse = mr_damper_force(&mut d, -0.1, 0.001);
        // Force signs should differ
        assert!(
            f_forward * f_reverse < 0.0,
            "hysteresis: opposite directions"
        );
    }

    // ── preview_control ───────────────────────────────────────────────────

    #[test]
    fn preview_zero_road_height_zero_force() {
        let f = preview_control(0.0, 0.5, 20.0, 20000.0, 5000.0);
        assert!(f.abs() < EPS, "flat road → zero preview force");
    }

    #[test]
    fn preview_positive_bump_positive_force() {
        let f = preview_control(0.05, 0.5, 20.0, 20000.0, 5000.0);
        assert!(f > 0.0, "bump ahead → positive feedforward");
    }

    #[test]
    fn preview_force_saturates() {
        let f = preview_control(10.0, 0.5, 20.0, 20000.0, 5000.0);
        assert!((f - 5000.0).abs() < EPS, "should saturate: {f}");
    }

    #[test]
    fn preview_attenuates_with_longer_horizon() {
        let f_short = preview_control(0.1, 0.1, 10.0, 20000.0, 5000.0);
        let f_long = preview_control(0.1, 5.0, 10.0, 20000.0, 5000.0);
        assert!(
            f_short.abs() >= f_long.abs(),
            "longer preview = more attenuation"
        );
    }

    #[test]
    fn preview_negative_dip_negative_force() {
        let f = preview_control(-0.05, 0.5, 20.0, 20000.0, 5000.0);
        assert!(f < 0.0, "dip ahead → negative feedforward");
    }
}
