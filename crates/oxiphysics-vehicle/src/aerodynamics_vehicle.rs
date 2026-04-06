// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle aerodynamics.
//!
//! Covers drag, lift, downforce, DRS, ground effect, wind triangle,
//! crosswind force, and cooling airflow.

// ---------------------------------------------------------------------------
// AeroBody
// ---------------------------------------------------------------------------

/// Aerodynamic properties of a vehicle body.
#[derive(Debug, Clone)]
pub struct AeroBody {
    /// Drag coefficient (dimensionless).
    pub drag_coeff: f64,
    /// Lift coefficient (dimensionless, positive = upward lift).
    pub lift_coeff: f64,
    /// Frontal reference area (m²).
    pub frontal_area: f64,
    /// Downforce coefficient (dimensionless, positive = downward force).
    pub downforce_coeff: f64,
    /// Air density (kg/m³).  Use ≈ 1.225 at sea level.
    pub air_density: f64,
}

impl AeroBody {
    /// Create an `AeroBody` with typical road-car defaults.
    pub fn default_road_car() -> Self {
        Self {
            drag_coeff: 0.30,
            lift_coeff: 0.10,
            frontal_area: 2.2,
            downforce_coeff: 0.05,
            air_density: 1.225,
        }
    }

    /// Create an `AeroBody` with typical open-wheel racing defaults.
    pub fn default_race_car() -> Self {
        Self {
            drag_coeff: 0.95,
            lift_coeff: -2.5, // net lift is negative (downforce)
            frontal_area: 1.5,
            downforce_coeff: 3.0,
            air_density: 1.225,
        }
    }
}

/// Compute the aerodynamic drag force (N).
///
/// `F_drag = 0.5 · ρ · Cd · A · v²`
pub fn aerodynamic_drag(body: &AeroBody, speed: f64) -> f64 {
    0.5 * body.air_density * body.drag_coeff * body.frontal_area * speed * speed
}

/// Compute the aerodynamic lift force (N).
///
/// Positive = upward lift; negative = downforce.
/// `F_lift = 0.5 · ρ · Cl · A · v²`
pub fn aerodynamic_lift(body: &AeroBody, speed: f64) -> f64 {
    0.5 * body.air_density * body.lift_coeff * body.frontal_area * speed * speed
}

/// Compute the aerodynamic downforce (N).
///
/// Uses `downforce_coeff`; positive = downward force.
/// `F_df = 0.5 · ρ · Cdf · A · v²`
pub fn aerodynamic_downforce(body: &AeroBody, speed: f64) -> f64 {
    0.5 * body.air_density * body.downforce_coeff * body.frontal_area * speed * speed
}

/// Compute the aerodynamic yaw moment (N·m).
///
/// `M_yaw = 0.5 · ρ · Cy · A · v² · sin(β)`
/// where `β` is the body side-slip angle.
pub fn yaw_moment(body: &AeroBody, side_slip_angle: f64, speed: f64, yaw_coeff: f64) -> f64 {
    0.5 * body.air_density * yaw_coeff * body.frontal_area * speed * speed * side_slip_angle.sin()
}

// ---------------------------------------------------------------------------
// DRS (Drag Reduction System)
// ---------------------------------------------------------------------------

/// Drag Reduction System configuration and state.
#[derive(Debug, Clone)]
pub struct DrsSystem {
    /// Baseline drag force at reference speed (N).  Used as a scale reference.
    pub base_drag: f64,
    /// Fractional drag reduction when DRS is active (0–1).
    pub drs_drag_reduction: f64,
    /// Whether DRS is currently deployed.
    pub drs_active: bool,
    /// Minimum speed (m/s) for DRS activation to be permitted.
    pub speed_threshold: f64,
}

impl DrsSystem {
    /// Create a default DRS system with 20 % drag reduction above 70 m/s.
    pub fn new_default() -> Self {
        Self {
            base_drag: 500.0,
            drs_drag_reduction: 0.20,
            drs_active: false,
            speed_threshold: 70.0,
        }
    }
}

/// Compute the effective drag force (N) with or without DRS.
///
/// If `drs.drs_active` and `speed >= drs.speed_threshold`, the drag is reduced
/// by `drs_drag_reduction * base_drag`.
pub fn drs_drag(drs: &DrsSystem, speed: f64) -> f64 {
    if drs.drs_active && speed >= drs.speed_threshold {
        drs.base_drag * (1.0 - drs.drs_drag_reduction)
    } else {
        drs.base_drag
    }
}

// ---------------------------------------------------------------------------
// Ground effect
// ---------------------------------------------------------------------------

/// Compute the downforce multiplier from ground effect.
///
/// As the vehicle approaches the ground, downforce increases roughly as
/// `F = base_downforce / height²` (inverse-square law, clamped for stability).
///
/// * `height_above_ground` — clearance in metres (clamped to ≥ 0.01 m).
/// * `base_downforce` — reference downforce at height = 1 m (N).
pub fn ground_effect(height_above_ground: f64, base_downforce: f64) -> f64 {
    let h = height_above_ground.max(0.01);
    base_downforce / (h * h)
}

// ---------------------------------------------------------------------------
// Wind triangle
// ---------------------------------------------------------------------------

/// Compute the apparent wind experienced by a moving vehicle.
///
/// Given the vehicle speed, its heading (radians, 0 = north), the true wind
/// speed, and the true wind direction (radians, meteorological convention:
/// direction the wind is coming *from*), returns `(apparent_speed, apparent_angle)`.
///
/// The apparent angle is measured relative to the vehicle's nose (0 = headwind,
/// π = tailwind, π/2 = side wind from port).
pub fn wind_triangle(
    vehicle_speed: f64,
    vehicle_heading: f64,
    wind_speed: f64,
    wind_dir: f64,
) -> (f64, f64) {
    // Vehicle velocity components (east, north convention)
    let vx = vehicle_speed * vehicle_heading.sin();
    let vy = vehicle_speed * vehicle_heading.cos();

    // True wind velocity: wind_dir is direction wind comes FROM
    let wx = -wind_speed * wind_dir.sin();
    let wy = -wind_speed * wind_dir.cos();

    // Apparent wind = true_wind - vehicle_vel
    let ax = wx - vx;
    let ay = wy - vy;

    let apparent_speed = (ax * ax + ay * ay).sqrt();
    // Angle relative to vehicle nose
    let apparent_angle = ax.atan2(ay) - vehicle_heading;
    (apparent_speed, apparent_angle)
}

// ---------------------------------------------------------------------------
// Crosswind force
// ---------------------------------------------------------------------------

/// Compute the crosswind side force on the vehicle (N).
///
/// `F_side = 0.5 · ρ · A · v² · sin²(angle)`
///
/// * `apparent_wind_speed` — apparent wind speed (m/s).
/// * `apparent_angle` — angle between wind and vehicle nose (rad).
/// * `area` — side projected area (m²).
/// * `density` — air density (kg/m³).
pub fn crosswind_force(
    apparent_wind_speed: f64,
    apparent_angle: f64,
    area: f64,
    density: f64,
) -> f64 {
    let sin_a = apparent_angle.sin();
    0.5 * density * area * apparent_wind_speed * apparent_wind_speed * sin_a * sin_a
}

// ---------------------------------------------------------------------------
// Cooling airflow
// ---------------------------------------------------------------------------

/// Estimate the cooling airflow mass-flow rate (kg/s).
///
/// `ṁ = coefficient · duct_area · speed · ρ_air`
///
/// Uses standard air density (1.225 kg/m³) internally.
pub fn cooling_airflow(speed: f64, duct_area: f64, coefficient: f64) -> f64 {
    const RHO_AIR: f64 = 1.225;
    coefficient * duct_area * speed * RHO_AIR
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    fn road_body() -> AeroBody {
        AeroBody {
            drag_coeff: 0.3,
            lift_coeff: 0.1,
            frontal_area: 2.0,
            downforce_coeff: 0.05,
            air_density: 1.225,
        }
    }

    // ── aerodynamic_drag ──────────────────────────────────────────────────

    #[test]
    fn drag_zero_speed_zero() {
        let body = road_body();
        assert!(aerodynamic_drag(&body, 0.0).abs() < EPS);
    }

    #[test]
    fn drag_proportional_to_speed_squared() {
        let body = road_body();
        let d1 = aerodynamic_drag(&body, 10.0);
        let d2 = aerodynamic_drag(&body, 20.0);
        assert!((d2 - 4.0 * d1).abs() < EPS, "drag ~ v²: d1={d1}, d2={d2}");
    }

    #[test]
    fn drag_positive_at_positive_speed() {
        let body = road_body();
        assert!(aerodynamic_drag(&body, 30.0) > 0.0);
    }

    #[test]
    fn drag_formula_check() {
        let body = AeroBody {
            drag_coeff: 0.3,
            lift_coeff: 0.0,
            frontal_area: 2.0,
            downforce_coeff: 0.0,
            air_density: 1.0,
        };
        // F = 0.5 * 1.0 * 0.3 * 2.0 * 100 = 30
        assert!((aerodynamic_drag(&body, 10.0) - 30.0).abs() < EPS);
    }

    // ── aerodynamic_lift ──────────────────────────────────────────────────

    #[test]
    fn lift_zero_speed_zero() {
        let body = road_body();
        assert!(aerodynamic_lift(&body, 0.0).abs() < EPS);
    }

    #[test]
    fn lift_positive_for_positive_lift_coeff() {
        let body = road_body(); // Cl = 0.1
        assert!(aerodynamic_lift(&body, 30.0) > 0.0);
    }

    #[test]
    fn lift_negative_for_negative_lift_coeff() {
        let body = AeroBody {
            drag_coeff: 0.3,
            lift_coeff: -2.0,
            frontal_area: 2.0,
            downforce_coeff: 0.0,
            air_density: 1.225,
        };
        assert!(aerodynamic_lift(&body, 30.0) < 0.0);
    }

    #[test]
    fn lift_proportional_to_speed_squared() {
        let body = road_body();
        let l1 = aerodynamic_lift(&body, 10.0);
        let l2 = aerodynamic_lift(&body, 20.0);
        assert!((l2 - 4.0 * l1).abs() < EPS);
    }

    // ── aerodynamic_downforce ─────────────────────────────────────────────

    #[test]
    fn downforce_zero_speed_zero() {
        let body = road_body();
        assert!(aerodynamic_downforce(&body, 0.0).abs() < EPS);
    }

    #[test]
    fn downforce_positive_for_positive_coeff() {
        let body = road_body();
        assert!(aerodynamic_downforce(&body, 50.0) > 0.0);
    }

    #[test]
    fn downforce_scales_with_density() {
        let mut b1 = road_body();
        let mut b2 = road_body();
        b2.air_density = b1.air_density * 2.0;
        b1.downforce_coeff = 0.5;
        b2.downforce_coeff = 0.5;
        let df1 = aerodynamic_downforce(&b1, 30.0);
        let df2 = aerodynamic_downforce(&b2, 30.0);
        assert!((df2 - 2.0 * df1).abs() < EPS);
    }

    // ── yaw_moment ────────────────────────────────────────────────────────

    #[test]
    fn yaw_moment_zero_slip_zero() {
        let body = road_body();
        let m = yaw_moment(&body, 0.0, 30.0, 0.15);
        assert!(m.abs() < EPS);
    }

    #[test]
    fn yaw_moment_nonzero_at_nonzero_slip() {
        let body = road_body();
        let m = yaw_moment(&body, 0.1, 30.0, 0.15);
        assert!(m.abs() > 0.0);
    }

    #[test]
    fn yaw_moment_symmetric_at_pi() {
        let body = road_body();
        let m = yaw_moment(&body, PI, 30.0, 0.15);
        assert!(m.abs() < EPS, "sin(π) ≈ 0 → zero moment");
    }

    // ── DrsSystem ─────────────────────────────────────────────────────────

    #[test]
    fn drs_inactive_returns_base_drag() {
        let drs = DrsSystem {
            base_drag: 500.0,
            drs_drag_reduction: 0.20,
            drs_active: false,
            speed_threshold: 70.0,
        };
        assert!((drs_drag(&drs, 80.0) - 500.0).abs() < EPS);
    }

    #[test]
    fn drs_active_above_threshold_reduces_drag() {
        let drs = DrsSystem {
            base_drag: 500.0,
            drs_drag_reduction: 0.20,
            drs_active: true,
            speed_threshold: 70.0,
        };
        let d = drs_drag(&drs, 80.0);
        assert!((d - 400.0).abs() < EPS, "expected 400, got {d}");
    }

    #[test]
    fn drs_active_below_threshold_no_reduction() {
        let drs = DrsSystem {
            base_drag: 500.0,
            drs_drag_reduction: 0.20,
            drs_active: true,
            speed_threshold: 70.0,
        };
        let d = drs_drag(&drs, 60.0);
        assert!((d - 500.0).abs() < EPS);
    }

    #[test]
    fn drs_default_new() {
        let drs = DrsSystem::new_default();
        assert!(!drs.drs_active);
        assert!((drs.speed_threshold - 70.0).abs() < EPS);
    }

    // ── ground_effect ─────────────────────────────────────────────────────

    #[test]
    fn ground_effect_unit_height_equals_base() {
        let df = ground_effect(1.0, 1000.0);
        assert!((df - 1000.0).abs() < EPS, "at h=1m, df=base: {df}");
    }

    #[test]
    fn ground_effect_lower_height_higher_downforce() {
        let df_high = ground_effect(0.5, 1000.0);
        let df_low = ground_effect(0.1, 1000.0);
        assert!(df_low > df_high);
    }

    #[test]
    fn ground_effect_clamps_near_zero_height() {
        // Should not panic or return infinity for very small height
        let df = ground_effect(0.001, 1000.0);
        assert!(df.is_finite());
    }

    #[test]
    fn ground_effect_inverse_square_at_half_height() {
        // At h=0.5: 1/0.25 * base = 4 * base
        let df = ground_effect(0.5, 1000.0);
        assert!((df - 4000.0).abs() < EPS);
    }

    // ── wind_triangle ─────────────────────────────────────────────────────

    #[test]
    fn wind_triangle_no_wind_apparent_equals_vehicle() {
        let (speed, _angle) = wind_triangle(20.0, 0.0, 0.0, 0.0);
        assert!(
            (speed - 20.0).abs() < EPS,
            "no wind: apparent = vehicle speed"
        );
    }

    #[test]
    fn wind_triangle_pure_headwind_adds_speeds() {
        // Vehicle north at 20 m/s, wind from north at 10 m/s
        let (speed, _angle) = wind_triangle(20.0, 0.0, 10.0, 0.0);
        assert!((speed - 30.0).abs() < EPS, "headwind adds: got {speed}");
    }

    #[test]
    fn wind_triangle_tailwind_reduces_apparent() {
        // Vehicle north at 20 m/s, wind from south (180°) at 10 m/s → tailwind
        let (speed, _angle) = wind_triangle(20.0, 0.0, 10.0, PI);
        assert!((speed - 10.0).abs() < EPS, "tailwind reduces: got {speed}");
    }

    #[test]
    fn wind_triangle_returns_finite_values() {
        let (speed, angle) = wind_triangle(30.0, 0.5, 5.0, 1.2);
        assert!(speed.is_finite() && angle.is_finite());
    }

    // ── crosswind_force ───────────────────────────────────────────────────

    #[test]
    fn crosswind_zero_angle_zero_force() {
        let f = crosswind_force(20.0, 0.0, 3.0, 1.225);
        assert!(f.abs() < EPS);
    }

    #[test]
    fn crosswind_max_at_ninety_degrees() {
        let f_90 = crosswind_force(20.0, PI / 2.0, 3.0, 1.225);
        let f_45 = crosswind_force(20.0, PI / 4.0, 3.0, 1.225);
        assert!(f_90 > f_45, "max crosswind force at 90°");
    }

    #[test]
    fn crosswind_positive_for_positive_speed() {
        let f = crosswind_force(20.0, PI / 4.0, 3.0, 1.225);
        assert!(f > 0.0);
    }

    #[test]
    fn crosswind_proportional_to_area() {
        let f1 = crosswind_force(20.0, PI / 4.0, 1.0, 1.225);
        let f2 = crosswind_force(20.0, PI / 4.0, 3.0, 1.225);
        assert!((f2 - 3.0 * f1).abs() < EPS);
    }

    // ── cooling_airflow ───────────────────────────────────────────────────

    #[test]
    fn cooling_airflow_zero_speed_zero() {
        assert!(cooling_airflow(0.0, 0.1, 0.8).abs() < EPS);
    }

    #[test]
    fn cooling_airflow_positive_at_speed() {
        assert!(cooling_airflow(30.0, 0.05, 0.7) > 0.0);
    }

    #[test]
    fn cooling_airflow_scales_linearly_with_speed() {
        let m1 = cooling_airflow(10.0, 0.05, 0.7);
        let m2 = cooling_airflow(20.0, 0.05, 0.7);
        assert!((m2 - 2.0 * m1).abs() < EPS);
    }

    #[test]
    fn cooling_airflow_scales_with_area() {
        let m1 = cooling_airflow(30.0, 0.05, 0.7);
        let m2 = cooling_airflow(30.0, 0.10, 0.7);
        assert!((m2 - 2.0 * m1).abs() < EPS);
    }
}
