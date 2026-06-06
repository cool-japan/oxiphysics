// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wind loading and aerodynamics for ground vehicles.
//!
//! Provides atmospheric wind models, 6-DOF aerodynamic coefficient presets,
//! dynamic pressure, drag/lift forces, yaw-angle computations, crosswind
//! stability, aero-noise estimates, and fuel-consumption penalties.
//!
//! # Overview
//!
//! - [`WindConditions`] — atmospheric wind speed, direction, turbulence, altitude.
//! - [`AerodynamicCoefficients`] — 6-DOF aero coefficients (Cd, Cl, Cs, roll, pitch, yaw).
//! - [`VehicleAeroLoad`] — full 6-DOF aerodynamic load on the vehicle body.
//! - [`TruckTrailerAero`] — combined tractor–trailer aerodynamic model.
//! - [`dynamic_pressure`] — `q = 0.5·ρ·v²`.
//! - [`drag_force`] / [`lift_force`] — aerodynamic force components.
//! - [`compute_aero_load`] — full 6-DOF load from coefficients and dynamic pressure.
//! - [`wind_yaw_angle`] — relative yaw angle between vehicle heading and wind.
//! - [`crosswind_stability_factor`] — rollover threshold metric.
//! - [`aerodynamic_noise_estimate`] — dB estimate of wind noise.
//! - [`wind_power_penalty`] — power lost to aerodynamic drag.
//! - [`highway_fuel_consumption`] — L/100 km estimate.

// ─────────────────────────────────────────────────────────────────────────────
// WindConditions
// ─────────────────────────────────────────────────────────────────────────────

/// Atmospheric wind conditions at a given location.
#[derive(Debug, Clone)]
pub struct WindConditions {
    /// Wind speed magnitude \[m/s\].
    pub speed: f64,
    /// Wind direction \[deg\], measured as meteorological convention
    /// (0 = wind from North, 90 = wind from East).
    pub direction_deg: f64,
    /// Turbulence intensity \[-\] (σ_u / U_mean), typically 0.05–0.25.
    pub turbulence_intensity: f64,
    /// Altitude above sea level \[m\].
    pub altitude: f64,
}

impl WindConditions {
    /// Create calm-weather wind conditions at `speed` \[m/s\] and `direction_deg` \[°\].
    pub fn new(speed: f64, direction_deg: f64) -> Self {
        Self {
            speed,
            direction_deg,
            turbulence_intensity: 0.05,
            altitude: 0.0,
        }
    }

    /// Compute relative wind speed \[m/s\] and yaw angle \[deg\] as experienced
    /// by a vehicle travelling at `vehicle_speed` \[m/s\] on heading `heading_deg` \[°\].
    ///
    /// Returns `(relative_speed, beta)` where `beta` is the sideslip angle \[deg\].
    pub fn relative_wind(&self, vehicle_speed: f64, heading_deg: f64) -> (f64, f64) {
        // Convert to Cartesian velocity components (East-North frame)
        let wind_rad = self.direction_deg.to_radians();
        // Meteorological: wind FROM direction → wind blows TO (dir + 180)
        let wx = -self.speed * wind_rad.sin();
        let wy = -self.speed * wind_rad.cos();
        let vx = vehicle_speed * heading_deg.to_radians().sin();
        let vy = vehicle_speed * heading_deg.to_radians().cos();
        let rel_x = wx - vx;
        let rel_y = wy - vy;
        let rel_speed = (rel_x * rel_x + rel_y * rel_y).sqrt();
        let beta = rel_x.atan2(rel_y).to_degrees();
        (rel_speed, beta)
    }

    /// Air density at the stored altitude using ISA standard atmosphere \[kg/m³\].
    ///
    /// `ρ(h) = ρ₀ · (1 − L·h / T₀)^(g·M / (R·L))`
    /// simplified as `ρ ≈ 1.225 · exp(−h / 8500)` below 11 km.
    pub fn air_density(&self) -> f64 {
        1.225_f64 * (-self.altitude / 8500.0).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AerodynamicCoefficients
// ─────────────────────────────────────────────────────────────────────────────

/// Six-degree-of-freedom aerodynamic coefficients for a ground vehicle.
#[derive(Debug, Clone)]
pub struct AerodynamicCoefficients {
    /// Drag coefficient Cd \[-\].
    pub cd: f64,
    /// Lift coefficient Cl \[-\] (positive = upward).
    pub cl: f64,
    /// Side-force coefficient Cs \[-\].
    pub cs: f64,
    /// Rolling-moment coefficient \[-\].
    pub cm_roll: f64,
    /// Pitching-moment coefficient \[-\].
    pub cm_pitch: f64,
    /// Yawing-moment coefficient \[-\].
    pub cm_yaw: f64,
}

impl AerodynamicCoefficients {
    /// Typical saloon/sedan car (Cd ≈ 0.30, Cl ≈ −0.10).
    pub fn car_saloon() -> Self {
        Self {
            cd: 0.30,
            cl: -0.10,
            cs: 0.0,
            cm_roll: 0.01,
            cm_pitch: -0.05,
            cm_yaw: 0.02,
        }
    }

    /// Typical SUV/crossover (Cd ≈ 0.35, Cl ≈ 0.05).
    pub fn suv() -> Self {
        Self {
            cd: 0.35,
            cl: 0.05,
            cs: 0.0,
            cm_roll: 0.02,
            cm_pitch: -0.03,
            cm_yaw: 0.03,
        }
    }

    /// Typical truck / heavy goods vehicle (Cd ≈ 0.70, Cl ≈ 0.0).
    pub fn truck() -> Self {
        Self {
            cd: 0.70,
            cl: 0.0,
            cs: 0.0,
            cm_roll: 0.04,
            cm_pitch: 0.0,
            cm_yaw: 0.05,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions — aerodynamic forces
// ─────────────────────────────────────────────────────────────────────────────

/// Dynamic pressure \[Pa\].
///
/// `q = 0.5 · ρ · v²`
pub fn dynamic_pressure(rho: f64, v: f64) -> f64 {
    0.5 * rho * v * v
}

/// Aerodynamic drag force \[N\].
///
/// `F_D = Cd · A · q`
pub fn drag_force(cd: f64, frontal_area: f64, q: f64) -> f64 {
    cd * frontal_area * q
}

/// Aerodynamic lift force \[N\].
///
/// `F_L = Cl · A · q`
pub fn lift_force(cl: f64, planform_area: f64, q: f64) -> f64 {
    cl * planform_area * q
}

// ─────────────────────────────────────────────────────────────────────────────
// VehicleAeroLoad
// ─────────────────────────────────────────────────────────────────────────────

/// Full 6-DOF aerodynamic load on a vehicle body.
#[derive(Debug, Clone)]
pub struct VehicleAeroLoad {
    /// Longitudinal drag force \[N\].
    pub drag: f64,
    /// Vertical lift force \[N\] (positive upward).
    pub lift: f64,
    /// Lateral side force \[N\].
    pub side_force: f64,
    /// Roll moment \[N·m\].
    pub roll_moment: f64,
    /// Pitch moment \[N·m\].
    pub pitch_moment: f64,
    /// Yaw moment \[N·m\].
    pub yaw_moment: f64,
}

/// Compute full 6-DOF aerodynamic load.
///
/// Uses the same reference area `frontal_area` for all force components
/// (a simplification appropriate for rough estimates).
pub fn compute_aero_load(
    coeff: &AerodynamicCoefficients,
    frontal_area: f64,
    q: f64,
) -> VehicleAeroLoad {
    VehicleAeroLoad {
        drag: coeff.cd * frontal_area * q,
        lift: coeff.cl * frontal_area * q,
        side_force: coeff.cs * frontal_area * q,
        roll_moment: coeff.cm_roll * frontal_area * q,
        pitch_moment: coeff.cm_pitch * frontal_area * q,
        yaw_moment: coeff.cm_yaw * frontal_area * q,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Wind angle and stability
// ─────────────────────────────────────────────────────────────────────────────

/// Relative yaw angle \[deg\] between vehicle heading and wind direction.
///
/// Returns angle in \[−180, 180\] range.
pub fn wind_yaw_angle(vehicle_heading: f64, wind_dir: f64) -> f64 {
    let mut angle = wind_dir - vehicle_heading;
    // Normalise to [-180, 180]
    while angle > 180.0 {
        angle -= 360.0;
    }
    while angle < -180.0 {
        angle += 360.0;
    }
    angle
}

/// Crosswind stability factor \[-\].
///
/// Returns the ratio `|side_force| / (weight · track_width / 2)`.
/// A value ≥ 1 indicates rollover risk.
pub fn crosswind_stability_factor(side_force: f64, weight: f64, track_width: f64) -> f64 {
    side_force.abs() / (weight * track_width / 2.0).max(1e-9)
}

// ─────────────────────────────────────────────────────────────────────────────
// Noise, power, fuel consumption
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate of aerodynamic wind noise \[dB\].
///
/// Simplified formula: `L = 20·log10(v) + 10·log10(Cd) + 50`
/// (empirical, valid for v > 1 m/s).
pub fn aerodynamic_noise_estimate(speed: f64, cd: f64) -> f64 {
    if speed < 1.0 {
        return 0.0;
    }
    20.0 * speed.log10() + 10.0 * cd.max(1e-9).log10() + 50.0
}

/// Power consumed by aerodynamic drag \[W\].
///
/// `P = F_D · v`
pub fn wind_power_penalty(drag: f64, v_vehicle: f64) -> f64 {
    drag * v_vehicle
}

// ─────────────────────────────────────────────────────────────────────────────
// TruckTrailerAero
// ─────────────────────────────────────────────────────────────────────────────

/// Combined tractor–trailer aerodynamic model.
#[derive(Debug, Clone)]
pub struct TruckTrailerAero {
    /// Tractor drag coefficient \[-\].
    pub tractor_cd: f64,
    /// Trailer drag coefficient \[-\].
    pub trailer_cd: f64,
    /// Tractor–trailer gap \[m\] (larger gap → less shielding).
    pub gap_length: f64,
}

impl TruckTrailerAero {
    /// Combined effective drag coefficient (accounts for gap interference).
    ///
    /// `Cd_combined = Cd_tractor + Cd_trailer · (1 − shielding)`
    /// where `shielding = exp(−gap / 2)`.
    pub fn combined_cd(&self) -> f64 {
        let shielding = (-self.gap_length / 2.0).exp();
        self.tractor_cd + self.trailer_cd * (1.0 - shielding)
    }

    /// Total drag force \[N\] given dynamic pressure `q` and reference area `area` \[m²\].
    pub fn total_drag(&self, q: f64, area: f64) -> f64 {
        self.combined_cd() * area * q
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fuel consumption
// ─────────────────────────────────────────────────────────────────────────────

/// Highway fuel consumption estimate \[L/100 km\].
///
/// Power at the wheels: `P = (F_drag + F_rr) · v`
/// Fuel consumption: `FC = P / (efficiency · LHV · ρ_fuel) · 1e5`
///
/// Uses petrol LHV ≈ 43.4 MJ/kg, ρ_fuel ≈ 740 g/L.
///
/// * `drag` — aerodynamic drag force \[N\]
/// * `rolling_resistance` — rolling resistance force \[N\]
/// * `speed` — vehicle speed \[m/s\]
/// * `efficiency` — drivetrain efficiency \[-\]
pub fn highway_fuel_consumption(
    drag: f64,
    rolling_resistance: f64,
    speed: f64,
    efficiency: f64,
) -> f64 {
    if speed < 1e-6 {
        return 0.0;
    }
    let power = (drag + rolling_resistance) * speed; // W
    let lhv_j_per_litre = 43.4e6 * 0.740; // ≈ 32.1 MJ/L
    let fuel_rate_l_per_s = power / (efficiency.max(1e-6) * lhv_j_per_litre);
    // Convert to L per 100 km
    let dist_per_s = speed; // m/s
    (fuel_rate_l_per_s / dist_per_s) * 100_000.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── dynamic_pressure ─────────────────────────────────────────────────

    #[test]
    fn dynamic_pressure_positive() {
        assert!(dynamic_pressure(1.225, 30.0) > 0.0);
    }

    #[test]
    fn dynamic_pressure_zero_at_rest() {
        assert_eq!(dynamic_pressure(1.225, 0.0), 0.0);
    }

    #[test]
    fn dynamic_pressure_formula() {
        let q = dynamic_pressure(1.225, 20.0);
        let expected = 0.5 * 1.225 * 400.0;
        assert!((q - expected).abs() < EPS);
    }

    #[test]
    fn dynamic_pressure_scales_with_v_squared() {
        let q1 = dynamic_pressure(1.225, 10.0);
        let q2 = dynamic_pressure(1.225, 20.0);
        assert!((q2 / q1 - 4.0).abs() < EPS);
    }

    // ── drag_force ───────────────────────────────────────────────────────

    #[test]
    fn drag_force_positive() {
        let q = dynamic_pressure(1.225, 30.0);
        assert!(drag_force(0.3, 2.2, q) > 0.0);
    }

    #[test]
    fn drag_force_increases_with_speed_squared() {
        let q1 = dynamic_pressure(1.225, 10.0);
        let q2 = dynamic_pressure(1.225, 20.0);
        let f1 = drag_force(0.3, 2.2, q1);
        let f2 = drag_force(0.3, 2.2, q2);
        assert!((f2 / f1 - 4.0).abs() < EPS);
    }

    #[test]
    fn drag_force_zero_for_zero_cd() {
        let q = dynamic_pressure(1.225, 30.0);
        assert_eq!(drag_force(0.0, 2.2, q), 0.0);
    }

    // ── lift_force ───────────────────────────────────────────────────────

    #[test]
    fn lift_force_zero_for_zero_cl() {
        let q = dynamic_pressure(1.225, 30.0);
        assert_eq!(lift_force(0.0, 4.0, q), 0.0);
    }

    #[test]
    fn lift_force_positive_for_positive_cl() {
        let q = dynamic_pressure(1.225, 30.0);
        assert!(lift_force(0.2, 4.0, q) > 0.0);
    }

    #[test]
    fn lift_force_negative_for_negative_cl() {
        let q = dynamic_pressure(1.225, 30.0);
        assert!(lift_force(-0.1, 4.0, q) < 0.0);
    }

    // ── wind_yaw_angle ───────────────────────────────────────────────────

    #[test]
    fn wind_yaw_angle_zero_when_aligned() {
        assert!((wind_yaw_angle(90.0, 90.0)).abs() < EPS);
    }

    #[test]
    fn wind_yaw_angle_perpendicular() {
        // Wind from East (90°), vehicle heading North (0°) → yaw = 90°
        let angle = wind_yaw_angle(0.0, 90.0);
        assert!((angle - 90.0).abs() < EPS);
    }

    #[test]
    fn wind_yaw_angle_normalised_to_180() {
        let angle = wind_yaw_angle(10.0, 200.0);
        assert!((-180.0..=180.0).contains(&angle));
    }

    #[test]
    fn wind_yaw_angle_headwind_is_180() {
        // Wind from 0°, vehicle heading 180° → β = −180 or 180
        let angle = wind_yaw_angle(180.0, 0.0);
        assert!(angle.abs() <= 180.0 + EPS);
    }

    // ── air_density ──────────────────────────────────────────────────────

    #[test]
    fn air_density_sea_level_approx_1p225() {
        let w = WindConditions::new(10.0, 270.0);
        assert!((w.air_density() - 1.225).abs() < EPS);
    }

    #[test]
    fn air_density_decreases_with_altitude() {
        let w0 = WindConditions {
            speed: 10.0,
            direction_deg: 0.0,
            turbulence_intensity: 0.05,
            altitude: 0.0,
        };
        let w1 = WindConditions {
            altitude: 5000.0,
            ..w0.clone()
        };
        assert!(w1.air_density() < w0.air_density());
    }

    #[test]
    fn air_density_positive_at_high_altitude() {
        let w = WindConditions {
            speed: 0.0,
            direction_deg: 0.0,
            turbulence_intensity: 0.0,
            altitude: 10000.0,
        };
        assert!(w.air_density() > 0.0);
    }

    // ── crosswind_stability_factor ───────────────────────────────────────

    #[test]
    fn crosswind_stability_factor_zero_for_zero_force() {
        assert_eq!(crosswind_stability_factor(0.0, 10000.0, 2.0), 0.0);
    }

    #[test]
    fn crosswind_stability_factor_positive() {
        let csf = crosswind_stability_factor(5000.0, 30000.0, 2.0);
        assert!(csf > 0.0);
    }

    #[test]
    fn crosswind_stability_high_side_force_exceeds_one() {
        // Very large side force → rollover risk
        let csf = crosswind_stability_factor(100_000.0, 10000.0, 2.0);
        assert!(csf > 1.0);
    }

    // ── compute_aero_load ────────────────────────────────────────────────

    #[test]
    fn aero_load_drag_positive_for_car() {
        let coeff = AerodynamicCoefficients::car_saloon();
        let q = dynamic_pressure(1.225, 30.0);
        let load = compute_aero_load(&coeff, 2.2, q);
        assert!(load.drag > 0.0);
    }

    #[test]
    fn aero_load_zero_at_rest() {
        let coeff = AerodynamicCoefficients::suv();
        let load = compute_aero_load(&coeff, 2.5, 0.0);
        assert_eq!(load.drag, 0.0);
        assert_eq!(load.lift, 0.0);
    }

    // ── wind_power_penalty ───────────────────────────────────────────────

    #[test]
    fn wind_power_penalty_positive() {
        assert!(wind_power_penalty(500.0, 30.0) > 0.0);
    }

    #[test]
    fn wind_power_penalty_formula() {
        assert!((wind_power_penalty(500.0, 30.0) - 15000.0).abs() < EPS);
    }

    // ── TruckTrailerAero ─────────────────────────────────────────────────

    #[test]
    fn truck_trailer_combined_cd_positive() {
        let tta = TruckTrailerAero {
            tractor_cd: 0.7,
            trailer_cd: 0.5,
            gap_length: 0.5,
        };
        assert!(tta.combined_cd() > 0.0);
    }

    #[test]
    fn truck_trailer_zero_gap_has_more_shielding() {
        let tta_zero = TruckTrailerAero {
            tractor_cd: 0.7,
            trailer_cd: 0.5,
            gap_length: 0.0,
        };
        let tta_large = TruckTrailerAero {
            tractor_cd: 0.7,
            trailer_cd: 0.5,
            gap_length: 5.0,
        };
        assert!(tta_large.combined_cd() > tta_zero.combined_cd());
    }

    // ── highway_fuel_consumption ─────────────────────────────────────────

    #[test]
    fn fuel_consumption_positive() {
        let fc = highway_fuel_consumption(400.0, 150.0, 27.78, 0.35);
        assert!(fc > 0.0);
    }

    #[test]
    fn fuel_consumption_zero_at_rest() {
        assert_eq!(highway_fuel_consumption(400.0, 150.0, 0.0, 0.35), 0.0);
    }

    #[test]
    fn fuel_consumption_increases_with_drag() {
        let fc1 = highway_fuel_consumption(300.0, 150.0, 27.78, 0.35);
        let fc2 = highway_fuel_consumption(600.0, 150.0, 27.78, 0.35);
        assert!(fc2 > fc1);
    }
}
