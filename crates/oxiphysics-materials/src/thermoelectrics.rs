// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermoelectric material models.
//!
//! Provides models for the Seebeck effect, Peltier effect, thermoelectric
//! generators (TEG), Peltier coolers, and figure-of-merit (ZT) interpolation.

// ---------------------------------------------------------------------------
// ThermoelectricProps
// ---------------------------------------------------------------------------

/// Bulk thermoelectric material properties.
pub struct ThermoelectricProps {
    /// Seebeck coefficient S \[V/K\].
    pub seebeck_coeff: f64,
    /// Electrical conductivity σ \[S/m\].
    pub electrical_conductivity: f64,
    /// Thermal conductivity κ \[W/(m·K)\].
    pub thermal_conductivity: f64,
    /// Absolute temperature T \[K\] at which the properties are evaluated.
    pub temperature: f64,
}

impl ThermoelectricProps {
    /// Create a new `ThermoelectricProps`.
    pub fn new(
        seebeck_coeff: f64,
        electrical_conductivity: f64,
        thermal_conductivity: f64,
        temperature: f64,
    ) -> Self {
        Self {
            seebeck_coeff,
            electrical_conductivity,
            thermal_conductivity,
            temperature,
        }
    }

    /// Thermoelectric power factor PF = S² · σ \[W/(m·K²)\].
    pub fn power_factor(&self) -> f64 {
        self.seebeck_coeff * self.seebeck_coeff * self.electrical_conductivity
    }

    /// Dimensionless figure of merit ZT = S² · σ · T / κ.
    pub fn figure_of_merit_zt(&self) -> f64 {
        self.power_factor() * self.temperature / self.thermal_conductivity
    }
}

// ---------------------------------------------------------------------------
// PeltierEffect
// ---------------------------------------------------------------------------

/// Peltier effect model at a junction.
pub struct PeltierEffect {
    /// Seebeck coefficient of the junction material \[V/K\].
    pub seebeck: f64,
    /// Junction absolute temperature T \[K\].
    pub temperature: f64,
}

impl PeltierEffect {
    /// Create a new `PeltierEffect`.
    pub fn new(seebeck: f64, temperature: f64) -> Self {
        Self {
            seebeck,
            temperature,
        }
    }

    /// Peltier coefficient Π = S · T \[V\].
    pub fn peltier_coefficient(&self) -> f64 {
        self.seebeck * self.temperature
    }

    /// Heat pumped per unit time Q = Π · I = S · T · I \[W\].
    pub fn heat_pumped(&self, current: f64) -> f64 {
        self.peltier_coefficient() * current
    }
}

// ---------------------------------------------------------------------------
// SeebeckEffect
// ---------------------------------------------------------------------------

/// Seebeck effect model for an open-circuit EMF.
pub struct SeebeckEffect {
    /// Seebeck coefficient S \[V/K\].
    pub seebeck_coeff: f64,
    /// Temperature difference ΔT = T_hot − T_cold \[K\].
    pub delta_t: f64,
}

impl SeebeckEffect {
    /// Create a new `SeebeckEffect`.
    pub fn new(seebeck_coeff: f64, delta_t: f64) -> Self {
        Self {
            seebeck_coeff,
            delta_t,
        }
    }

    /// Open-circuit EMF: V = S · ΔT \[V\].
    pub fn emf(&self) -> f64 {
        self.seebeck_coeff * self.delta_t
    }
}

// ---------------------------------------------------------------------------
// ThermoelectricGenerator
// ---------------------------------------------------------------------------

/// Thermoelectric generator (TEG) model.
///
/// Converts a temperature gradient into electrical power using the Seebeck
/// effect. The matched-load maximum power and efficiency are computed via
/// the standard TEG relations.
pub struct ThermoelectricGenerator {
    /// Hot-side temperature T_h \[K\].
    pub hot_temp: f64,
    /// Cold-side temperature T_c \[K\].
    pub cold_temp: f64,
    /// Seebeck coefficient S \[V/K\].
    pub seebeck: f64,
    /// Internal electrical resistance R \[Ω\].
    pub resistance: f64,
    /// Thermal conductance K \[W/K\].
    pub thermal_conductance: f64,
}

impl ThermoelectricGenerator {
    /// Create a new `ThermoelectricGenerator`.
    pub fn new(
        hot_temp: f64,
        cold_temp: f64,
        seebeck: f64,
        resistance: f64,
        thermal_conductance: f64,
    ) -> Self {
        Self {
            hot_temp,
            cold_temp,
            seebeck,
            resistance,
            thermal_conductance,
        }
    }

    /// Temperature difference ΔT = T_h − T_c \[K\].
    pub fn delta_t(&self) -> f64 {
        self.hot_temp - self.cold_temp
    }

    /// Maximum electrical power at matched load: P_max = S²·ΔT² / (4·R) \[W\].
    pub fn max_power(&self) -> f64 {
        let dt = self.delta_t();
        self.seebeck * self.seebeck * dt * dt / (4.0 * self.resistance)
    }

    /// Carnot efficiency η_C = ΔT / T_h.
    pub fn carnot_efficiency(&self) -> f64 {
        self.delta_t() / self.hot_temp
    }

    /// TEG efficiency at matched load using the standard formula:
    ///
    /// η = η_C · (√(1 + ZT_m) − 1) / (√(1 + ZT_m) + T_c / T_h)
    ///
    /// where ZT_m = S² / (R · K) · T_mean, T_mean = (T_h + T_c) / 2.
    pub fn efficiency(&self) -> f64 {
        let t_mean = (self.hot_temp + self.cold_temp) / 2.0;
        let zt_m =
            self.seebeck * self.seebeck / (self.resistance * self.thermal_conductance) * t_mean;
        let sqrt_term = (1.0 + zt_m).sqrt();
        let eta_c = self.carnot_efficiency();
        eta_c * (sqrt_term - 1.0) / (sqrt_term + self.cold_temp / self.hot_temp)
    }
}

// ---------------------------------------------------------------------------
// PeltierCooler
// ---------------------------------------------------------------------------

/// Peltier cooler (reverse-Seebeck) module.
///
/// Operates in reverse mode: electrical current drives heat from the cold
/// side to the hot side.
pub struct PeltierCooler {
    /// Seebeck coefficient S \[V/K\].
    pub seebeck: f64,
    /// Internal resistance R \[Ω\].
    pub resistance: f64,
    /// Thermal conductance K \[W/K\].
    pub thermal_conductance: f64,
    /// Hot-side temperature T_h \[K\].
    pub hot_temp: f64,
    /// Cold-side temperature T_c \[K\].
    pub cold_temp: f64,
}

impl PeltierCooler {
    /// Create a new `PeltierCooler`.
    pub fn new(
        seebeck: f64,
        resistance: f64,
        thermal_conductance: f64,
        hot_temp: f64,
        cold_temp: f64,
    ) -> Self {
        Self {
            seebeck,
            resistance,
            thermal_conductance,
            hot_temp,
            cold_temp,
        }
    }

    /// Cooling power Q_c = S·T_c·I − ½·R·I² − K·ΔT \[W\].
    ///
    /// Returns the net heat removed from the cold side at current `i` \[A\].
    pub fn cooling_power(&self, i: f64) -> f64 {
        self.seebeck * self.cold_temp * i
            - 0.5 * self.resistance * i * i
            - self.thermal_conductance * (self.hot_temp - self.cold_temp)
    }

    /// Coefficient of performance COP = Q_c / P_elec.
    ///
    /// `P_elec = S·ΔT·I + R·I²` (electrical power consumed).
    pub fn coefficient_of_performance(&self, i: f64) -> f64 {
        let q_c = self.cooling_power(i);
        let dt = self.hot_temp - self.cold_temp;
        let p_elec = self.seebeck * dt * i + self.resistance * i * i;
        if p_elec.abs() < 1e-30 {
            0.0
        } else {
            q_c / p_elec
        }
    }

    /// Maximum temperature difference achievable (Q_c = 0 solved for ΔT_max):
    ///
    /// ΔT_max = (½ · ZT · T_c²) where Z = S² / (R · K).
    pub fn max_temperature_difference(&self) -> f64 {
        let z = self.seebeck * self.seebeck / (self.resistance * self.thermal_conductance);
        0.5 * z * self.cold_temp * self.cold_temp
    }
}

// ---------------------------------------------------------------------------
// MaterialZT
// ---------------------------------------------------------------------------

/// Temperature-dependent ZT data for a material, supporting linear interpolation.
pub struct MaterialZT {
    /// Evaluation temperature \[K\] (used as a query default).
    pub temperature: f64,
    /// Tabulated (T \[K\], ZT) data pairs, sorted ascending by T.
    pub zt_values: Vec<(f64, f64)>,
}

impl MaterialZT {
    /// Create a new `MaterialZT` with the given (T, ZT) table.
    ///
    /// The table need not be pre-sorted; it will be sorted internally.
    pub fn new(temperature: f64, mut zt_values: Vec<(f64, f64)>) -> Self {
        zt_values.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self {
            temperature,
            zt_values,
        }
    }

    /// Linearly interpolate ZT at temperature `t` \[K\].
    ///
    /// Clamps to the boundary values outside the tabulated range.
    pub fn interpolate(&self, t: f64) -> f64 {
        let data = &self.zt_values;
        if data.is_empty() {
            return 0.0;
        }
        if t <= data.first().expect("collection should not be empty").0 {
            return data.first().expect("collection should not be empty").1;
        }
        if t >= data.last().expect("collection should not be empty").0 {
            return data.last().expect("collection should not be empty").1;
        }
        // Binary search for the interval.
        let pos = data.partition_point(|&(ti, _)| ti <= t);
        let (t0, z0) = data[pos - 1];
        let (t1, z1) = data[pos];
        let frac = (t - t0) / (t1 - t0);
        z0 + frac * (z1 - z0)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Thomson coefficient τ = T · dS/dT \[V/K\].
///
/// # Arguments
/// * `seebeck` – Seebeck coefficient S at temperature T \[V/K\].
/// * `t` – Absolute temperature T \[K\].
/// * `d_seebeck_dt` – Derivative dS/dT at temperature T \[V/K²\].
pub fn thomson_coefficient(seebeck: f64, t: f64, d_seebeck_dt: f64) -> f64 {
    let _ = seebeck; // S itself is not needed for τ = T·dS/dT
    t * d_seebeck_dt
}

/// Thermoelectric reduced efficiency η_r for a device with Carnot efficiency
/// η_C and figure of merit ZT.
///
/// η_r = η_C · (√(1 + ZT) − 1) / (√(1 + ZT) + 1 − η_C)
pub fn reduced_efficiency(eta_c: f64, zt: f64) -> f64 {
    let m = (1.0 + zt).sqrt();
    eta_c * (m - 1.0) / (m + 1.0 - eta_c)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- ThermoelectricProps ---

    #[test]
    fn test_power_factor_basic() {
        let tp = ThermoelectricProps::new(200e-6, 1e5, 1.5, 300.0);
        let pf = tp.power_factor();
        let expected = (200e-6_f64).powi(2) * 1e5;
        assert!((pf - expected).abs() < 1e-20, "pf={pf} expected={expected}");
    }

    #[test]
    fn test_figure_of_merit_zt() {
        let tp = ThermoelectricProps::new(200e-6, 1e5, 1.0, 300.0);
        let zt = tp.figure_of_merit_zt();
        let expected = tp.power_factor() * 300.0 / 1.0;
        assert!((zt - expected).abs() < 1e-20);
    }

    #[test]
    fn test_zt_doubles_with_doubled_seebeck_squared() {
        let tp1 = ThermoelectricProps::new(100e-6, 1e5, 1.0, 300.0);
        let tp2 = ThermoelectricProps::new(200e-6, 1e5, 1.0, 300.0);
        assert!((tp2.figure_of_merit_zt() - 4.0 * tp1.figure_of_merit_zt()).abs() < 1e-20);
    }

    #[test]
    fn test_zt_high_quality_material() {
        // Bi₂Te₃-class: S ≈ 200 µV/K, σ ≈ 1e5 S/m, κ ≈ 1.5 W/(m·K), T = 300 K → ZT ≈ 0.8
        let tp = ThermoelectricProps::new(200e-6, 1e5, 1.5, 300.0);
        let zt = tp.figure_of_merit_zt();
        assert!(zt > 0.5 && zt < 2.0, "ZT={zt} should be in 0.5–2.0 range");
    }

    #[test]
    fn test_power_factor_zero_seebeck() {
        let tp = ThermoelectricProps::new(0.0, 1e5, 1.5, 300.0);
        assert!(tp.power_factor().abs() < EPS);
        assert!(tp.figure_of_merit_zt().abs() < EPS);
    }

    // --- PeltierEffect ---

    #[test]
    fn test_peltier_coefficient() {
        let pe = PeltierEffect::new(200e-6, 300.0);
        assert!((pe.peltier_coefficient() - 200e-6 * 300.0).abs() < EPS);
    }

    #[test]
    fn test_heat_pumped_proportional_to_current() {
        let pe = PeltierEffect::new(200e-6, 300.0);
        let q1 = pe.heat_pumped(1.0);
        let q2 = pe.heat_pumped(2.0);
        assert!((q2 - 2.0 * q1).abs() < EPS);
    }

    #[test]
    fn test_peltier_zero_temperature() {
        let pe = PeltierEffect::new(200e-6, 0.0);
        assert!(pe.peltier_coefficient().abs() < EPS);
        assert!(pe.heat_pumped(5.0).abs() < EPS);
    }

    #[test]
    fn test_peltier_negative_current() {
        let pe = PeltierEffect::new(200e-6, 300.0);
        assert!(pe.heat_pumped(-1.0) < 0.0);
    }

    // --- SeebeckEffect ---

    #[test]
    fn test_seebeck_emf_basic() {
        let se = SeebeckEffect::new(200e-6, 100.0);
        assert!((se.emf() - 0.02).abs() < EPS);
    }

    #[test]
    fn test_seebeck_zero_delta_t() {
        let se = SeebeckEffect::new(200e-6, 0.0);
        assert!(se.emf().abs() < EPS);
    }

    #[test]
    fn test_seebeck_negative_gradient() {
        let se = SeebeckEffect::new(200e-6, -50.0);
        assert!(se.emf() < 0.0);
    }

    #[test]
    fn test_seebeck_proportional_to_delta_t() {
        let se1 = SeebeckEffect::new(100e-6, 50.0);
        let se2 = SeebeckEffect::new(100e-6, 100.0);
        assert!((se2.emf() - 2.0 * se1.emf()).abs() < EPS);
    }

    // --- ThermoelectricGenerator ---

    #[test]
    fn test_teg_max_power_positive() {
        let teg = ThermoelectricGenerator::new(600.0, 300.0, 200e-6, 1.0, 0.01);
        assert!(teg.max_power() > 0.0);
    }

    #[test]
    fn test_teg_max_power_formula() {
        let s = 200e-6_f64;
        let r = 1.0_f64;
        let dt = 300.0_f64;
        let teg = ThermoelectricGenerator::new(600.0, 300.0, s, r, 0.01);
        let expected = s * s * dt * dt / (4.0 * r);
        assert!((teg.max_power() - expected).abs() < 1e-20);
    }

    #[test]
    fn test_teg_carnot_efficiency() {
        let teg = ThermoelectricGenerator::new(600.0, 300.0, 200e-6, 1.0, 0.01);
        let eta_c = teg.carnot_efficiency();
        assert!((eta_c - 0.5).abs() < EPS);
    }

    #[test]
    fn test_teg_efficiency_less_than_carnot() {
        let teg = ThermoelectricGenerator::new(600.0, 300.0, 200e-6, 1.0, 0.01);
        assert!(teg.efficiency() < teg.carnot_efficiency());
    }

    #[test]
    fn test_teg_efficiency_positive() {
        let teg = ThermoelectricGenerator::new(600.0, 300.0, 200e-6, 1.0, 0.01);
        assert!(teg.efficiency() > 0.0);
    }

    #[test]
    fn test_teg_zero_delta_t_zero_power() {
        let teg = ThermoelectricGenerator::new(300.0, 300.0, 200e-6, 1.0, 0.01);
        assert!(teg.max_power().abs() < EPS);
    }

    #[test]
    fn test_teg_delta_t() {
        let teg = ThermoelectricGenerator::new(700.0, 300.0, 200e-6, 1.0, 0.01);
        assert!((teg.delta_t() - 400.0).abs() < EPS);
    }

    // --- PeltierCooler ---

    #[test]
    fn test_cooler_cooling_power_positive_at_low_current() {
        // At very low current the conduction term should dominate — check sign behaviour.
        let pc = PeltierCooler::new(200e-6, 0.1, 0.01, 320.0, 280.0);
        // At I = 0: Q_c = -K·ΔT < 0 (no pumping)
        let q0 = pc.cooling_power(0.0);
        assert!(q0 < 0.0, "No current → heat leaks in, Q_c={q0}");
    }

    #[test]
    fn test_cooler_cooling_power_large_current_positive() {
        // High current: Peltier pumping dominates.
        let pc = PeltierCooler::new(200e-6, 0.001, 1e-4, 320.0, 280.0);
        let q_c = pc.cooling_power(10.0);
        assert!(
            q_c > 0.0,
            "High current should give positive cooling, got {q_c}"
        );
    }

    #[test]
    fn test_cooler_max_temperature_difference_positive() {
        let pc = PeltierCooler::new(200e-6, 0.1, 0.01, 320.0, 280.0);
        let dt_max = pc.max_temperature_difference();
        assert!(dt_max > 0.0, "ΔT_max should be positive, got {dt_max}");
    }

    #[test]
    fn test_cooler_cop_finite() {
        let pc = PeltierCooler::new(200e-6, 0.001, 1e-4, 320.0, 280.0);
        let cop = pc.coefficient_of_performance(5.0);
        assert!(cop.is_finite(), "COP should be finite");
    }

    #[test]
    fn test_cooler_max_dt_scales_with_z() {
        // Doubling S doubles Z → doubles ΔT_max (for fixed R, K).
        let pc1 = PeltierCooler::new(100e-6, 1.0, 1.0, 320.0, 300.0);
        let pc2 = PeltierCooler::new(200e-6, 1.0, 1.0, 320.0, 300.0);
        let ratio = pc2.max_temperature_difference() / pc1.max_temperature_difference();
        assert!((ratio - 4.0).abs() < EPS, "ratio={ratio}");
    }

    // --- MaterialZT ---

    #[test]
    fn test_materialzt_interpolate_midpoint() {
        let mzt = MaterialZT::new(350.0, vec![(300.0, 0.8), (400.0, 1.2)]);
        let zt = mzt.interpolate(350.0);
        assert!((zt - 1.0).abs() < 1e-10, "zt={zt}");
    }

    #[test]
    fn test_materialzt_clamp_below() {
        let mzt = MaterialZT::new(300.0, vec![(300.0, 0.8), (400.0, 1.2)]);
        assert!((mzt.interpolate(200.0) - 0.8).abs() < EPS);
    }

    #[test]
    fn test_materialzt_clamp_above() {
        let mzt = MaterialZT::new(300.0, vec![(300.0, 0.8), (400.0, 1.2)]);
        assert!((mzt.interpolate(500.0) - 1.2).abs() < EPS);
    }

    #[test]
    fn test_materialzt_exact_point() {
        let mzt = MaterialZT::new(300.0, vec![(300.0, 0.8), (350.0, 1.0), (400.0, 1.2)]);
        assert!((mzt.interpolate(350.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_materialzt_empty() {
        let mzt = MaterialZT::new(300.0, vec![]);
        assert!(mzt.interpolate(300.0).abs() < EPS);
    }

    #[test]
    fn test_materialzt_unsorted_input() {
        // Unsorted input should be sorted internally.
        let mzt = MaterialZT::new(300.0, vec![(400.0, 1.2), (300.0, 0.8)]);
        let zt = mzt.interpolate(350.0);
        assert!((zt - 1.0).abs() < 1e-10, "zt={zt}");
    }

    // --- Free functions ---

    #[test]
    fn test_thomson_coefficient_basic() {
        // τ = T · dS/dT
        let tau = thomson_coefficient(200e-6, 300.0, 0.5e-6);
        assert!((tau - 300.0 * 0.5e-6).abs() < EPS);
    }

    #[test]
    fn test_thomson_coefficient_zero_derivative() {
        assert!(thomson_coefficient(200e-6, 300.0, 0.0).abs() < EPS);
    }

    #[test]
    fn test_thomson_coefficient_zero_temperature() {
        assert!(thomson_coefficient(200e-6, 0.0, 0.5e-6).abs() < EPS);
    }

    #[test]
    fn test_reduced_efficiency_less_than_carnot() {
        // η_r < η_C for any finite ZT.
        let eta_c = 0.5;
        let eta_r = reduced_efficiency(eta_c, 1.0);
        assert!(eta_r < eta_c, "eta_r={eta_r}, eta_c={eta_c}");
    }

    #[test]
    fn test_reduced_efficiency_increases_with_zt() {
        let eta_c = 0.5;
        let e1 = reduced_efficiency(eta_c, 1.0);
        let e2 = reduced_efficiency(eta_c, 4.0);
        assert!(e2 > e1, "Higher ZT should give higher efficiency");
    }

    #[test]
    fn test_reduced_efficiency_zero_carnot() {
        assert!(reduced_efficiency(0.0, 2.0).abs() < EPS);
    }

    #[test]
    fn test_reduced_efficiency_formula_spot_check() {
        // eta_c=0.5, zt=3: sqrt(4)=2, result = 0.5*(2-1)/(2+1-0.5) = 0.5/2.5 = 0.2
        let result = reduced_efficiency(0.5, 3.0);
        assert!((result - 0.2).abs() < EPS, "result={result}");
    }
}
