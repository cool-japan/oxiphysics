// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fuel system and combustion modeling.
//!
//! Provides physics-based models for:
//! - BSFC (brake-specific fuel consumption) based fuel consumption.
//! - Injection timing advance angle.
//! - Fuel spray atomization (Sauter mean diameter).
//! - Lambda-based combustion efficiency.
//! - Thermal NOx (Zeldovich mechanism) emissions.
//! - Simplified pendulum model for fuel tank sloshing.

/// Fuel system parameters and state.
///
/// Models a single-tank, single-injector fuel system.
#[derive(Debug, Clone)]
pub struct FuelSystem {
    /// Total tank volume in litres.
    pub tank_volume_l: f64,
    /// Current fuel level as a fraction of tank volume (0.0–1.0).
    pub fuel_level: f64,
    /// Fuel density in kg/L (gasoline ≈ 0.745, diesel ≈ 0.832).
    pub fuel_density_kg_l: f64,
    /// Maximum injection rate in mg/stroke.
    pub max_injection_mg: f64,
    /// Brake-specific fuel consumption at rated power (g/kWh).
    pub bsfc_g_per_kwh: f64,
}

impl FuelSystem {
    /// Create a fuel system with typical gasoline parameters.
    ///
    /// # Arguments
    /// * `tank_volume_l` — Tank volume in litres.
    /// * `fuel_level`    — Initial fill fraction (0.0–1.0).
    pub fn new_gasoline(tank_volume_l: f64, fuel_level: f64) -> Self {
        Self {
            tank_volume_l,
            fuel_level: fuel_level.clamp(0.0, 1.0),
            fuel_density_kg_l: 0.745,
            max_injection_mg: 60.0,
            bsfc_g_per_kwh: 250.0,
        }
    }

    /// Create a fuel system with typical diesel parameters.
    ///
    /// # Arguments
    /// * `tank_volume_l` — Tank volume in litres.
    /// * `fuel_level`    — Initial fill fraction (0.0–1.0).
    pub fn new_diesel(tank_volume_l: f64, fuel_level: f64) -> Self {
        Self {
            tank_volume_l,
            fuel_level: fuel_level.clamp(0.0, 1.0),
            fuel_density_kg_l: 0.832,
            max_injection_mg: 80.0,
            bsfc_g_per_kwh: 210.0,
        }
    }

    /// Current fuel mass in the tank (kg).
    pub fn fuel_mass_kg(&self) -> f64 {
        self.tank_volume_l * self.fuel_level * self.fuel_density_kg_l
    }

    /// Consume `mass_kg` of fuel, updating `fuel_level`.
    ///
    /// Returns the amount of fuel actually consumed (limited by availability).
    pub fn consume_fuel(&mut self, mass_kg: f64) -> f64 {
        let available = self.fuel_mass_kg();
        let consumed = mass_kg.min(available).max(0.0);
        let consumed_l = consumed / self.fuel_density_kg_l;
        self.fuel_level -= consumed_l / self.tank_volume_l;
        self.fuel_level = self.fuel_level.max(0.0);
        consumed
    }
}

/// Compute instantaneous fuel consumption rate using BSFC.
///
/// # Arguments
/// * `power_kw`    — Engine output power in kW.
/// * `bsfc_g_kwh`  — Brake-specific fuel consumption in g/(kW·h).
///
/// Returns mass flow rate in g/s.
pub fn fuel_consumption(power_kw: f64, bsfc_g_kwh: f64) -> f64 {
    // mdot [g/s] = P [kW] * BSFC [g/kWh] / 3600
    (power_kw * bsfc_g_kwh) / 3600.0
}

/// Compute optimal injection advance angle in crank degrees.
///
/// Uses a simplified model based on engine speed and load:
/// `advance = base + speed_coeff * rpm + load_coeff * (1 - load)`.
///
/// # Arguments
/// * `rpm`          — Engine speed in rev/min.
/// * `load`         — Normalised engine load (0.0 = idle, 1.0 = full).
/// * `base_deg`     — Base advance angle in crank degrees BTDC.
///
/// Returns advance angle in crank degrees BTDC.
pub fn injection_timing(rpm: f64, load: f64, base_deg: f64) -> f64 {
    let speed_coeff = 0.004; // deg per rpm
    let load_coeff = 5.0; // deg at zero load
    base_deg + speed_coeff * rpm + load_coeff * (1.0 - load.clamp(0.0, 1.0))
}

/// Sauter mean diameter (SMD) of a fuel spray in micrometres.
///
/// Uses an empirical correlation (Lefebvre 1989):
/// `SMD = C * (μ_f / (ρ_a * U_rel))^a * (σ / (ρ_a * U_rel^2 * d))^b`
///
/// In simplified non-dimensional form here:
/// `SMD ≈ 3.08 * (μ_f^0.385) * (σ^0.737) / (ρ_a^0.737 * U_rel^1.669 * d^0.06)`
///
/// # Arguments
/// * `u_rel_m_s`   — Relative velocity of fuel jet vs air (m/s).
/// * `d_nozzle_mm` — Nozzle orifice diameter (mm).
/// * `sigma_n_m`   — Surface tension of fuel (N/m); gasoline ≈ 0.022.
/// * `mu_f_pas`    — Dynamic viscosity of fuel (Pa·s); gasoline ≈ 3e-4.
/// * `rho_air`     — Air density (kg/m³); ambient ≈ 1.2.
///
/// Returns SMD in micrometres.
pub fn fuel_atomization(
    u_rel_m_s: f64,
    d_nozzle_mm: f64,
    sigma_n_m: f64,
    mu_f_pas: f64,
    rho_air: f64,
) -> f64 {
    let d = d_nozzle_mm * 1e-3; // m
    let u = u_rel_m_s.max(1e-6);
    // Empirical Lefebvre correlation (simplified)
    3.08e6 * mu_f_pas.powf(0.385) * sigma_n_m.powf(0.737)
        / (rho_air.powf(0.737) * u.powf(1.669) * d.powf(0.06))
}

/// Lambda-based combustion efficiency.
///
/// Lambda (λ) is the actual air-fuel ratio divided by the stoichiometric AFR.
/// - λ < 1: rich mixture — efficiency < 1 due to incomplete combustion.
/// - λ = 1: stoichiometric — peak efficiency (≈ 0.98).
/// - λ > 1: lean mixture — efficiency drops at very lean conditions.
///
/// # Arguments
/// * `lambda` — Air-excess ratio (λ = AFR / AFR_stoich).
///
/// Returns combustion efficiency in \[0, 1\].
pub fn combustion_efficiency(lambda: f64) -> f64 {
    let lambda = lambda.max(0.01);
    if lambda < 1.0 {
        // Rich: linearly drops below stoich
        (0.98 * lambda).clamp(0.0, 1.0)
    } else if lambda <= 1.3 {
        // Lean but within normal range
        let t = (lambda - 1.0) / 0.3;
        0.98 - 0.15 * t
    } else {
        // Very lean — rapid drop
        let t = (lambda - 1.3).min(1.0);
        (0.83 - 0.4 * t).max(0.0)
    }
}

/// Thermal NOx formation rate (Zeldovich mechanism).
///
/// Simplified Zeldovich expression:
/// `d[NO]/dt ≈ A * exp(-E_a / (R * T)) * [O2]^0.5`
///
/// The function returns a relative NOx index (dimensionless, arbitrary units)
/// suitable for comparing operating conditions.
///
/// # Arguments
/// * `temp_k`       — Peak flame temperature in kelvin.
/// * `o2_fraction`  — Mole fraction of O₂ in the charge (0–0.23).
///
/// Returns relative NOx formation index.
pub fn emissions_nox(temp_k: f64, o2_fraction: f64) -> f64 {
    // Activation temperature for N2 + O → NO + N (≈ 38 000 K)
    let e_a_over_r = 38_000.0_f64;
    let pre_exp = 6e16_f64; // s⁻¹ mol⁻⁰·⁵ m¹·⁵ (scaled)
    let t = temp_k.max(300.0);
    let o2 = o2_fraction.clamp(0.0, 0.23);
    pre_exp * (-e_a_over_r / t).exp() * o2.sqrt()
}

/// Tank sloshing pendulum model.
///
/// Models fuel sloshing as a simple pendulum with effective length `l_eff`
/// and damping ratio `zeta`. Returns the pendulum angle in radians given:
///
/// # Arguments
/// * `lateral_accel_m_s2` — Lateral acceleration of the vehicle (m/s²).
/// * `fill_fraction`      — Fuel fill fraction (0–1); full tanks slosh less.
/// * `omega_n_rad_s`      — Natural frequency of the slosh mode (rad/s).
/// * `zeta`               — Damping ratio (0 = undamped, 1 = critically damped).
/// * `time_s`             — Time since lateral acceleration applied (s).
///
/// Returns the approximate sloshing displacement angle (rad) at time `t`.
pub fn tank_sloshing(
    lateral_accel_m_s2: f64,
    fill_fraction: f64,
    omega_n_rad_s: f64,
    zeta: f64,
    time_s: f64,
) -> f64 {
    let g = 9.81_f64;
    // Static equilibrium angle under lateral acceleration
    let theta_static = (lateral_accel_m_s2 / g).atan();
    // Fill effect: full tank has smaller effective free surface
    let fill_factor = (1.0 - fill_fraction.clamp(0.0, 1.0)).sqrt().max(0.1);
    let omega_n = omega_n_rad_s.max(0.1) * fill_factor;
    let zeta = zeta.clamp(0.0, 2.0);
    let t = time_s.max(0.0);

    if zeta < 1.0 {
        // Underdamped oscillation
        let omega_d = omega_n * (1.0 - zeta * zeta).sqrt();
        theta_static * (1.0 - (-zeta * omega_n * t).exp() * (omega_d * t).cos())
    } else {
        // Overdamped — exponential approach
        theta_static * (1.0 - (-omega_n * t).exp())
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── FuelSystem construction ───────────────────────────────────────────

    #[test]
    fn gasoline_tank_creates_correct_density() {
        let fs = FuelSystem::new_gasoline(50.0, 1.0);
        assert!((fs.fuel_density_kg_l - 0.745).abs() < EPS);
    }

    #[test]
    fn diesel_tank_creates_correct_density() {
        let fs = FuelSystem::new_diesel(60.0, 1.0);
        assert!((fs.fuel_density_kg_l - 0.832).abs() < EPS);
    }

    #[test]
    fn fuel_level_clamped_above_one() {
        let fs = FuelSystem::new_gasoline(50.0, 1.5);
        assert!((fs.fuel_level - 1.0).abs() < EPS);
    }

    #[test]
    fn fuel_level_clamped_below_zero() {
        let fs = FuelSystem::new_gasoline(50.0, -0.1);
        assert!((fs.fuel_level).abs() < EPS);
    }

    #[test]
    fn fuel_mass_full_tank() {
        let fs = FuelSystem::new_gasoline(50.0, 1.0);
        let expected = 50.0 * 0.745;
        assert!((fs.fuel_mass_kg() - expected).abs() < 1e-6);
    }

    #[test]
    fn fuel_mass_half_tank() {
        let fs = FuelSystem::new_gasoline(50.0, 0.5);
        let expected = 25.0 * 0.745;
        assert!((fs.fuel_mass_kg() - expected).abs() < 1e-6);
    }

    // ── consume_fuel ──────────────────────────────────────────────────────

    #[test]
    fn consume_fuel_reduces_level() {
        let mut fs = FuelSystem::new_gasoline(50.0, 1.0);
        let initial_mass = fs.fuel_mass_kg();
        let consumed = fs.consume_fuel(5.0);
        assert!((consumed - 5.0).abs() < 1e-6);
        assert!((fs.fuel_mass_kg() - (initial_mass - 5.0)).abs() < 1e-4);
    }

    #[test]
    fn consume_fuel_cannot_exceed_available() {
        let mut fs = FuelSystem::new_gasoline(10.0, 0.1);
        let available = fs.fuel_mass_kg();
        let consumed = fs.consume_fuel(1000.0);
        assert!((consumed - available).abs() < 1e-4);
        assert!(fs.fuel_level < EPS);
    }

    #[test]
    fn consume_zero_fuel_no_change() {
        let mut fs = FuelSystem::new_gasoline(50.0, 0.8);
        let level_before = fs.fuel_level;
        fs.consume_fuel(0.0);
        assert!((fs.fuel_level - level_before).abs() < EPS);
    }

    // ── fuel_consumption ──────────────────────────────────────────────────

    #[test]
    fn fuel_consumption_zero_power() {
        assert!(fuel_consumption(0.0, 250.0).abs() < EPS);
    }

    #[test]
    fn fuel_consumption_known_value() {
        // 100 kW * 250 g/kWh / 3600 ≈ 6.944 g/s
        let rate = fuel_consumption(100.0, 250.0);
        assert!((rate - 100.0 * 250.0 / 3600.0).abs() < 1e-6);
    }

    #[test]
    fn fuel_consumption_scales_linearly_with_power() {
        let r1 = fuel_consumption(50.0, 300.0);
        let r2 = fuel_consumption(100.0, 300.0);
        assert!((r2 - 2.0 * r1).abs() < 1e-10);
    }

    // ── injection_timing ──────────────────────────────────────────────────

    #[test]
    fn injection_timing_full_load_minimum_advance() {
        // At full load, load_coeff term = 0
        let adv = injection_timing(3000.0, 1.0, 10.0);
        assert!(adv > 10.0); // base + speed contribution
    }

    #[test]
    fn injection_timing_low_load_more_advance() {
        let adv_full = injection_timing(3000.0, 1.0, 10.0);
        let adv_low = injection_timing(3000.0, 0.0, 10.0);
        assert!(adv_low > adv_full, "low load should advance more");
    }

    #[test]
    fn injection_timing_higher_rpm_more_advance() {
        let adv_low_rpm = injection_timing(1000.0, 0.5, 5.0);
        let adv_high_rpm = injection_timing(5000.0, 0.5, 5.0);
        assert!(adv_high_rpm > adv_low_rpm);
    }

    // ── fuel_atomization ──────────────────────────────────────────────────

    #[test]
    fn atomization_positive_smd() {
        let smd = fuel_atomization(50.0, 0.2, 0.022, 3e-4, 1.2);
        assert!(smd > 0.0);
    }

    #[test]
    fn atomization_higher_velocity_smaller_droplets() {
        let smd_low = fuel_atomization(10.0, 0.2, 0.022, 3e-4, 1.2);
        let smd_high = fuel_atomization(100.0, 0.2, 0.022, 3e-4, 1.2);
        assert!(smd_high < smd_low, "faster injection → smaller SMD");
    }

    // ── combustion_efficiency ─────────────────────────────────────────────

    #[test]
    fn combustion_efficiency_at_stoich_near_peak() {
        let eff = combustion_efficiency(1.0);
        assert!((eff - 0.98).abs() < 1e-6);
    }

    #[test]
    fn combustion_efficiency_rich_below_stoich() {
        let eff_rich = combustion_efficiency(0.8);
        let eff_stoich = combustion_efficiency(1.0);
        assert!(eff_rich < eff_stoich);
    }

    #[test]
    fn combustion_efficiency_very_lean_drops() {
        let eff_lean = combustion_efficiency(1.5);
        let eff_stoich = combustion_efficiency(1.0);
        assert!(eff_lean < eff_stoich);
    }

    #[test]
    fn combustion_efficiency_clipped_to_zero() {
        let eff = combustion_efficiency(0.001);
        assert!(eff >= 0.0);
    }

    #[test]
    fn combustion_efficiency_bounded_to_one() {
        for &lambda in &[0.5, 1.0, 1.5, 2.0] {
            let eff = combustion_efficiency(lambda);
            assert!(eff <= 1.0, "efficiency={eff} for lambda={lambda}");
            assert!(eff >= 0.0);
        }
    }

    // ── emissions_nox ─────────────────────────────────────────────────────

    #[test]
    fn nox_zero_at_room_temperature() {
        // At 300 K the activation temperature >> T → exp ≈ 0
        let nox = emissions_nox(300.0, 0.21);
        assert!(nox < 1e-10, "NOx at room temp should be negligible: {nox}");
    }

    #[test]
    fn nox_increases_with_temperature() {
        let nox_low = emissions_nox(2000.0, 0.21);
        let nox_high = emissions_nox(2500.0, 0.21);
        assert!(nox_high > nox_low, "NOx should increase with temperature");
    }

    #[test]
    fn nox_increases_with_o2() {
        let nox_low = emissions_nox(2200.0, 0.05);
        let nox_high = emissions_nox(2200.0, 0.21);
        assert!(nox_high > nox_low, "more O2 → more NOx");
    }

    #[test]
    fn nox_zero_o2_gives_zero() {
        let nox = emissions_nox(2500.0, 0.0);
        assert!(nox < EPS, "zero O2 → zero NOx");
    }

    // ── tank_sloshing ─────────────────────────────────────────────────────

    #[test]
    fn sloshing_zero_at_time_zero() {
        let theta = tank_sloshing(5.0, 0.5, 3.0, 0.3, 0.0);
        assert!(theta.abs() < 1e-6, "angle at t=0 should be ~0");
    }

    #[test]
    fn sloshing_approaches_static_equilibrium() {
        let lat = 5.0_f64;
        let theta_inf = tank_sloshing(lat, 0.5, 3.0, 1.5, 100.0);
        let static_eq = (lat / 9.81).atan();
        assert!(
            (theta_inf - static_eq).abs() < 0.05,
            "should approach static eq"
        );
    }

    #[test]
    fn sloshing_larger_for_less_fuel() {
        let theta_full = tank_sloshing(5.0, 1.0, 3.0, 0.3, 0.5);
        let theta_empty = tank_sloshing(5.0, 0.1, 3.0, 0.3, 0.5);
        assert!(
            theta_empty.abs() > theta_full.abs(),
            "less fuel → more sloshing"
        );
    }

    #[test]
    fn sloshing_zero_lateral_accel() {
        let theta = tank_sloshing(0.0, 0.5, 3.0, 0.3, 1.0);
        assert!(theta.abs() < 1e-9, "no lateral accel → no slosh");
    }
}
