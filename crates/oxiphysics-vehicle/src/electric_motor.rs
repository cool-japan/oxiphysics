// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electric motor modelling: torque-speed curves, efficiency maps,
//! regenerative braking, thermal model, battery SoC, and regen limits.

// ---------------------------------------------------------------------------
// Efficiency map helpers
// ---------------------------------------------------------------------------

/// Bilinear lookup in a 2-D efficiency map.
///
/// `torque_points` and `speed_points` define the grid axes.  `map` is stored
/// in row-major order `map[speed_idx * n_torque + torque_idx]`.
fn bilinear_lookup(
    map: &[f64],
    torque_points: &[f64],
    speed_points: &[f64],
    torque: f64,
    speed: f64,
) -> f64 {
    let nt = torque_points.len();
    let ns = speed_points.len();
    if nt == 0 || ns == 0 {
        return 0.0;
    }

    // Clamp and find lower-bound index for torque
    let torque_c = torque.clamp(torque_points[0], torque_points[nt - 1]);
    let ti = torque_points
        .windows(2)
        .position(|w| torque_c >= w[0] && torque_c <= w[1])
        .unwrap_or(nt - 2)
        .min(nt - 2);

    // Clamp and find lower-bound index for speed
    let speed_c = speed.clamp(speed_points[0], speed_points[ns - 1]);
    let si = speed_points
        .windows(2)
        .position(|w| speed_c >= w[0] && speed_c <= w[1])
        .unwrap_or(ns - 2)
        .min(ns - 2);

    let ft = (torque_c - torque_points[ti]) / (torque_points[ti + 1] - torque_points[ti] + 1e-15);
    let fs = (speed_c - speed_points[si]) / (speed_points[si + 1] - speed_points[si] + 1e-15);

    let h00 = map[si * nt + ti];
    let h10 = map[si * nt + ti + 1];
    let h01 = map[(si + 1) * nt + ti];
    let h11 = map[(si + 1) * nt + ti + 1];

    h00 * (1.0 - ft) * (1.0 - fs) + h10 * ft * (1.0 - fs) + h01 * (1.0 - ft) * fs + h11 * ft * fs
}

// ---------------------------------------------------------------------------
// ElectricMotor
// ---------------------------------------------------------------------------

/// Electric motor specification and runtime state.
#[derive(Debug, Clone)]
pub struct ElectricMotor {
    /// Peak torque at standstill (N·m).
    pub peak_torque: f64,
    /// Rated (base) speed at which constant-torque region ends (RPM).
    pub rated_rpm: f64,
    /// Peak mechanical power (W).
    pub peak_power: f64,
    /// Maximum allowable continuous temperature (°C).
    pub max_temp_celsius: f64,
    /// Winding thermal resistance (°C/W).
    pub thermal_resistance: f64,
    /// Winding thermal capacitance (J/°C).
    pub thermal_capacitance: f64,
    /// Current winding temperature (°C).
    pub temperature: f64,
    /// Ambient temperature (°C).
    pub ambient_temp: f64,
    /// Torque grid points for the efficiency map (N·m).
    pub eff_torque_pts: Vec<f64>,
    /// Speed grid points for the efficiency map (RPM).
    pub eff_speed_pts: Vec<f64>,
    /// Efficiency map values \[0, 1\] in row-major order (speed × torque).
    pub eff_map: Vec<f64>,
    /// Battery state of charge \[0, 1\].
    pub soc: f64,
    /// Battery capacity (A·h).
    pub battery_capacity_ah: f64,
    /// Nominal battery voltage (V).
    pub nominal_voltage: f64,
}

impl ElectricMotor {
    /// Create a simple electric motor with a constant efficiency approximation.
    ///
    /// # Arguments
    /// * `peak_torque` - Peak torque (N·m).
    /// * `rated_rpm` - Base speed (RPM).
    /// * `peak_power` - Peak power (W).
    /// * `nominal_efficiency` - Flat efficiency for the simple map \[0, 1\].
    /// * `battery_capacity_ah` - Battery capacity (A·h).
    /// * `nominal_voltage` - Nominal battery voltage (V).
    pub fn new(
        peak_torque: f64,
        rated_rpm: f64,
        peak_power: f64,
        nominal_efficiency: f64,
        battery_capacity_ah: f64,
        nominal_voltage: f64,
    ) -> Self {
        // Build a 2×2 efficiency map for simplicity
        let eff_torque_pts = vec![0.0, peak_torque];
        let eff_speed_pts = vec![0.0, rated_rpm * 2.0];
        let eff = nominal_efficiency.clamp(0.0, 1.0);
        let eff_map = vec![eff; 4];
        Self {
            peak_torque,
            rated_rpm,
            peak_power,
            max_temp_celsius: 180.0,
            thermal_resistance: 0.1,
            thermal_capacitance: 500.0,
            temperature: 25.0,
            ambient_temp: 25.0,
            eff_torque_pts,
            eff_speed_pts,
            eff_map,
            soc: 1.0,
            battery_capacity_ah,
            nominal_voltage,
        }
    }
}

// ---------------------------------------------------------------------------
// torque_speed_curve
// ---------------------------------------------------------------------------

/// Compute the available torque at a given speed using the standard
/// constant-torque / constant-power characteristic.
///
/// * Below `rated_rpm`: output torque = `peak_torque`.
/// * Above `rated_rpm`: output torque = `peak_power / omega` (constant power),
///   clamped to `[0, peak_torque]`.
///
/// # Arguments
/// * `motor` - Motor parameters.
/// * `rpm` - Current shaft speed (RPM).
///
/// Returns available torque (N·m).
pub fn torque_speed_curve(motor: &ElectricMotor, rpm: f64) -> f64 {
    let rpm = rpm.max(0.0);
    if rpm <= motor.rated_rpm {
        motor.peak_torque
    } else {
        let omega = rpm * std::f64::consts::PI / 30.0; // rad/s
        (motor.peak_power / omega.max(1e-6)).min(motor.peak_torque)
    }
}

// ---------------------------------------------------------------------------
// motor_efficiency
// ---------------------------------------------------------------------------

/// Look up the motor efficiency for a given operating point.
///
/// # Arguments
/// * `motor` - Motor with efficiency map.
/// * `torque` - Shaft torque (N·m), clamped to map range.
/// * `rpm` - Shaft speed (RPM), clamped to map range.
///
/// Returns efficiency in `[0, 1]`.
pub fn motor_efficiency(motor: &ElectricMotor, torque: f64, rpm: f64) -> f64 {
    bilinear_lookup(
        &motor.eff_map,
        &motor.eff_torque_pts,
        &motor.eff_speed_pts,
        torque.abs(),
        rpm.abs(),
    )
    .clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// regenerative_braking
// ---------------------------------------------------------------------------

/// Compute the regenerative braking torque available.
///
/// Regen is limited by:
/// 1. The maximum torque the motor can absorb (same curve as motoring).
/// 2. The battery's ability to accept charge (simplified as full acceptance
///    when SoC < `max_soc_for_regen`).
/// 3. Motor temperature — no regen above `max_temp_celsius`.
///
/// # Arguments
/// * `motor` - Motor state.
/// * `rpm` - Shaft speed (RPM).
/// * `max_soc_for_regen` - Maximum SoC above which regen is disabled \[0, 1\].
///
/// Returns the regen torque magnitude (N·m, always ≥ 0).
pub fn regenerative_braking(motor: &ElectricMotor, rpm: f64, max_soc_for_regen: f64) -> f64 {
    if motor.soc >= max_soc_for_regen {
        return 0.0;
    }
    if motor.temperature >= motor.max_temp_celsius {
        return 0.0;
    }
    torque_speed_curve(motor, rpm)
}

// ---------------------------------------------------------------------------
// thermal_model
// ---------------------------------------------------------------------------

/// Update the motor winding temperature using a first-order thermal model.
///
/// `dT/dt = (P_loss - (T - T_ambient) / R_th) / C_th`
///
/// Uses explicit Euler with time step `dt`.
///
/// # Arguments
/// * `motor` - Mutable motor state (temperature is updated in place).
/// * `torque` - Shaft torque (N·m).
/// * `rpm` - Shaft speed (RPM).
/// * `dt` - Time step (s).
pub fn thermal_model(motor: &mut ElectricMotor, torque: f64, rpm: f64, dt: f64) {
    let omega = rpm * std::f64::consts::PI / 30.0;
    let p_mech = torque.abs() * omega;
    let eff = motor_efficiency(motor, torque, rpm).max(1e-6);
    // Electrical input power
    let p_elec = p_mech / eff;
    let p_loss = (p_elec - p_mech).max(0.0);
    let cooling = (motor.temperature - motor.ambient_temp) / motor.thermal_resistance;
    let d_temp = (p_loss - cooling) / motor.thermal_capacitance;
    motor.temperature += d_temp * dt;
}

// ---------------------------------------------------------------------------
// battery_soc
// ---------------------------------------------------------------------------

/// Update battery state of charge from a current integral.
///
/// `SoC(t + dt) = SoC(t) - current * dt / (capacity_ah * 3600)`
///
/// A positive `current_a` discharges the battery; negative charges it.
///
/// # Arguments
/// * `motor` - Mutable motor state (SoC updated in place).
/// * `current_a` - Battery current (A); positive = discharge.
/// * `dt` - Time step (s).
pub fn battery_soc(motor: &mut ElectricMotor, current_a: f64, dt: f64) {
    let capacity_as = motor.battery_capacity_ah * 3600.0;
    motor.soc -= current_a * dt / capacity_as;
    motor.soc = motor.soc.clamp(0.0, 1.0);
}

// ---------------------------------------------------------------------------
// max_regen_torque
// ---------------------------------------------------------------------------

/// Compute the maximum regenerative braking torque limited by SoC and
/// temperature.
///
/// Linearly reduces available regen torque as SoC approaches
/// `max_soc_for_regen` and as temperature approaches `max_temp_celsius`.
///
/// # Arguments
/// * `motor` - Motor state.
/// * `rpm` - Shaft speed (RPM).
/// * `max_soc_for_regen` - SoC limit above which regen is zero \[0, 1\].
///
/// Returns the maximum regen torque (N·m, ≥ 0).
pub fn max_regen_torque(motor: &ElectricMotor, rpm: f64, max_soc_for_regen: f64) -> f64 {
    let base = torque_speed_curve(motor, rpm);
    // SoC derating
    let soc_factor = if motor.soc >= max_soc_for_regen {
        0.0
    } else {
        1.0 - motor.soc / max_soc_for_regen.max(1e-6)
    };
    // Temperature derating: linear from (max_temp - 20) to max_temp
    let temp_onset = (motor.max_temp_celsius - 20.0).max(0.0);
    let temp_factor = if motor.temperature >= motor.max_temp_celsius {
        0.0
    } else if motor.temperature >= temp_onset {
        1.0 - (motor.temperature - temp_onset) / 20.0_f64.max(1e-6)
    } else {
        1.0
    };
    base * soc_factor * temp_factor
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn default_motor() -> ElectricMotor {
        ElectricMotor::new(200.0, 3000.0, 60_000.0, 0.92, 60.0, 400.0)
    }

    // ── torque_speed_curve ─────────────────────────────────────────────────

    #[test]
    fn torque_at_zero_rpm_is_peak() {
        let m = default_motor();
        let t = torque_speed_curve(&m, 0.0);
        assert!((t - m.peak_torque).abs() < 1e-6, "t={t}");
    }

    #[test]
    fn torque_below_rated_rpm_is_peak() {
        let m = default_motor();
        let t = torque_speed_curve(&m, m.rated_rpm * 0.5);
        assert!((t - m.peak_torque).abs() < 1e-6, "t={t}");
    }

    #[test]
    fn torque_at_rated_rpm_is_peak() {
        let m = default_motor();
        let t = torque_speed_curve(&m, m.rated_rpm);
        assert!((t - m.peak_torque).abs() < 1e-6, "t={t}");
    }

    #[test]
    fn torque_above_rated_rpm_decreases() {
        let m = default_motor();
        let t_high = torque_speed_curve(&m, m.rated_rpm * 2.0);
        assert!(t_high < m.peak_torque, "t_high={t_high}");
    }

    #[test]
    fn torque_product_constant_power_region() {
        let m = default_motor();
        let rpm = m.rated_rpm * 3.0;
        let omega = rpm * std::f64::consts::PI / 30.0;
        let t = torque_speed_curve(&m, rpm);
        let power = t * omega;
        assert!(
            (power - m.peak_power).abs() / m.peak_power < 0.01,
            "power={power}"
        );
    }

    #[test]
    fn torque_negative_rpm_clamped() {
        let m = default_motor();
        let t = torque_speed_curve(&m, -100.0);
        assert!((t - m.peak_torque).abs() < 1e-6);
    }

    // ── motor_efficiency ───────────────────────────────────────────────────

    #[test]
    fn efficiency_in_unit_range() {
        let m = default_motor();
        let eff = motor_efficiency(&m, 100.0, 1500.0);
        assert!((0.0..=1.0).contains(&eff), "eff={eff}");
    }

    #[test]
    fn efficiency_matches_nominal() {
        let m = default_motor();
        let eff = motor_efficiency(&m, 100.0, 1500.0);
        assert!((eff - 0.92).abs() < 1e-6, "eff={eff}");
    }

    #[test]
    fn efficiency_zero_torque_in_map() {
        let m = default_motor();
        let eff = motor_efficiency(&m, 0.0, 0.0);
        assert!((0.0..=1.0).contains(&eff));
    }

    #[test]
    fn efficiency_symmetric_torque_sign() {
        let m = default_motor();
        let e1 = motor_efficiency(&m, 50.0, 1000.0);
        let e2 = motor_efficiency(&m, -50.0, 1000.0);
        assert!((e1 - e2).abs() < 1e-10);
    }

    // ── regenerative_braking ───────────────────────────────────────────────

    #[test]
    fn regen_zero_when_soc_at_limit() {
        let mut m = default_motor();
        m.soc = 0.9;
        let t = regenerative_braking(&m, 1500.0, 0.9);
        assert!(t < 1e-10, "t={t}");
    }

    #[test]
    fn regen_positive_when_soc_below_limit() {
        let mut m = default_motor();
        m.soc = 0.5; // set SoC below the regen limit
        let t = regenerative_braking(&m, 1500.0, 0.9);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn regen_zero_above_max_temperature() {
        let mut m = default_motor();
        m.temperature = m.max_temp_celsius + 1.0;
        let t = regenerative_braking(&m, 1500.0, 0.9);
        assert!(t < 1e-10, "t={t}");
    }

    #[test]
    fn regen_torque_follows_speed_curve() {
        let mut m = default_motor();
        m.soc = 0.5; // set SoC below the regen limit so regen is active
        let t_regen = regenerative_braking(&m, 1500.0, 0.9);
        let t_curve = torque_speed_curve(&m, 1500.0);
        assert!(
            (t_regen - t_curve).abs() < 1e-6,
            "t_regen={t_regen} t_curve={t_curve}"
        );
    }

    // ── thermal_model ──────────────────────────────────────────────────────

    #[test]
    fn thermal_model_increases_temp_under_load() {
        let mut m = default_motor();
        let t_initial = m.temperature;
        thermal_model(&mut m, 150.0, 2000.0, 10.0);
        assert!(m.temperature > t_initial, "temp should rise under load");
    }

    #[test]
    fn thermal_model_no_load_cools_to_ambient() {
        let mut m = default_motor();
        m.temperature = 100.0;
        // Run many steps with zero torque
        for _ in 0..1000 {
            thermal_model(&mut m, 0.0, 0.0, 1.0);
        }
        assert!(
            m.temperature < 50.0,
            "temp should cool toward ambient, got {}",
            m.temperature
        );
    }

    #[test]
    fn thermal_model_temp_stays_finite() {
        let mut m = default_motor();
        for _ in 0..100 {
            thermal_model(&mut m, 200.0, 5000.0, 0.1);
        }
        assert!(m.temperature.is_finite());
    }

    // ── battery_soc ────────────────────────────────────────────────────────

    #[test]
    fn soc_decreases_on_discharge() {
        let mut m = default_motor();
        let soc_before = m.soc;
        battery_soc(&mut m, 50.0, 1.0);
        assert!(m.soc < soc_before, "SoC should decrease");
    }

    #[test]
    fn soc_increases_on_charge() {
        let mut m = default_motor();
        m.soc = 0.5;
        let soc_before = m.soc;
        battery_soc(&mut m, -50.0, 1.0);
        assert!(
            m.soc > soc_before,
            "SoC should increase on negative current"
        );
    }

    #[test]
    fn soc_clamped_at_zero() {
        let mut m = default_motor();
        m.soc = 0.0;
        battery_soc(&mut m, 1000.0, 3600.0);
        assert!(m.soc >= 0.0, "SoC must not go negative");
    }

    #[test]
    fn soc_clamped_at_one() {
        let mut m = default_motor();
        m.soc = 1.0;
        battery_soc(&mut m, -1000.0, 3600.0);
        assert!(m.soc <= 1.0, "SoC must not exceed 1");
    }

    #[test]
    fn soc_full_cycle_returns_approximately_original() {
        let mut m = default_motor();
        m.soc = 0.5;
        let cap = m.battery_capacity_ah * 3600.0;
        let current = 10.0;
        let dt = cap / current / 10.0; // 1/10 of capacity in dt
        battery_soc(&mut m, current, dt); // discharge
        battery_soc(&mut m, -current, dt); // charge back
        assert!((m.soc - 0.5).abs() < 1e-10, "SoC={}", m.soc);
    }

    // ── max_regen_torque ───────────────────────────────────────────────────

    #[test]
    fn max_regen_zero_at_max_soc() {
        let mut m = default_motor();
        m.soc = 0.95;
        let t = max_regen_torque(&m, 2000.0, 0.95);
        assert!(t < 1e-10, "t={t}");
    }

    #[test]
    fn max_regen_positive_low_soc() {
        let mut m = default_motor();
        m.soc = 0.5; // set SoC well below the regen limit
        let t = max_regen_torque(&m, 2000.0, 0.95);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn max_regen_derated_near_max_soc() {
        let mut m = default_motor();
        m.soc = 0.0;
        let t_low = max_regen_torque(&m, 1500.0, 0.95);
        m.soc = 0.8;
        let t_high = max_regen_torque(&m, 1500.0, 0.95);
        assert!(t_high < t_low, "t_high={t_high} t_low={t_low}");
    }

    #[test]
    fn max_regen_zero_over_temp() {
        let mut m = default_motor();
        m.temperature = m.max_temp_celsius + 5.0;
        let t = max_regen_torque(&m, 2000.0, 0.95);
        assert!(t < 1e-10, "t={t}");
    }

    #[test]
    fn max_regen_lte_torque_curve() {
        let m = default_motor();
        let t_regen = max_regen_torque(&m, 1500.0, 0.95);
        let t_curve = torque_speed_curve(&m, 1500.0);
        assert!(t_regen <= t_curve + 1e-10);
    }
}
