// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tire wear and degradation models: Archard wear, temperature layers,
//! compound types, pressure, grip degradation, and lap-level simulation.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// CompoundType
// ---------------------------------------------------------------------------

/// Racing or road tire compound classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundType {
    /// Soft compound — maximum grip, low durability.
    Soft,
    /// Medium compound — balanced grip and durability.
    Medium,
    /// Hard compound — low grip, high durability.
    Hard,
    /// Intermediate wet compound.
    Intermediate,
    /// Full wet compound.
    Wet,
}

impl CompoundType {
    /// Baseline grip coefficient at optimal temperature.
    pub fn grip_coefficient(self) -> f64 {
        match self {
            CompoundType::Soft => 1.60,
            CompoundType::Medium => 1.45,
            CompoundType::Hard => 1.30,
            CompoundType::Intermediate => 1.20,
            CompoundType::Wet => 1.10,
        }
    }

    /// Durability factor (higher = longer life). Relative scale.
    pub fn durability_factor(self) -> f64 {
        match self {
            CompoundType::Soft => 0.50,
            CompoundType::Medium => 1.00,
            CompoundType::Hard => 2.00,
            CompoundType::Intermediate => 1.50,
            CompoundType::Wet => 1.80,
        }
    }

    /// Optimal temperature range (K): (min, peak, max).
    pub fn optimal_temp_range(self) -> (f64, f64, f64) {
        match self {
            CompoundType::Soft => (353.15, 368.15, 383.15),
            CompoundType::Medium => (343.15, 363.15, 383.15),
            CompoundType::Hard => (333.15, 358.15, 383.15),
            CompoundType::Intermediate => (303.15, 323.15, 343.15),
            CompoundType::Wet => (283.15, 298.15, 313.15),
        }
    }

    /// Returns `true` for dry-weather compounds.
    pub fn is_dry(self) -> bool {
        matches!(
            self,
            CompoundType::Soft | CompoundType::Medium | CompoundType::Hard
        )
    }
}

// ---------------------------------------------------------------------------
// TireTemperature
// ---------------------------------------------------------------------------

/// Three-layer tire temperature model: surface, bulk, and carcass.
#[derive(Debug, Clone)]
pub struct TireTemperature {
    /// Surface temperature (K).
    pub surface: f64,
    /// Bulk (midtread) temperature (K).
    pub bulk: f64,
    /// Carcass temperature (K).
    pub carcass: f64,
    /// Ambient temperature (K).
    pub ambient: f64,
    /// Surface thermal capacitance (J/K).
    pub cap_surface: f64,
    /// Bulk thermal capacitance (J/K).
    pub cap_bulk: f64,
    /// Carcass thermal capacitance (J/K).
    pub cap_carcass: f64,
    /// Conductance surface→bulk (W/K).
    pub cond_sb: f64,
    /// Conductance bulk→carcass (W/K).
    pub cond_bc: f64,
    /// Convective conductance surface→ambient (W/K).
    pub conv_surface: f64,
    /// Convective conductance carcass→ambient (W/K).
    pub conv_carcass: f64,
}

impl TireTemperature {
    /// Initialise all layers at ambient temperature.
    pub fn new_at_ambient(ambient: f64) -> Self {
        Self {
            surface: ambient,
            bulk: ambient,
            carcass: ambient,
            ambient,
            cap_surface: 1_500.0,
            cap_bulk: 3_000.0,
            cap_carcass: 4_000.0,
            cond_sb: 8.0,
            cond_bc: 5.0,
            conv_surface: 25.0,
            conv_carcass: 10.0,
        }
    }

    /// Default for a medium-compound racing tire.
    pub fn default_medium() -> Self {
        Self::new_at_ambient(293.15)
    }

    /// Integrate temperatures over `dt` (s) with `heat_input` (W) at the surface.
    pub fn step(&mut self, dt: f64, heat_input: f64) {
        // Surface: heat input + conduction from bulk − convection − conduction to bulk
        let q_sb = self.cond_sb * (self.surface - self.bulk);
        let q_s_conv = self.conv_surface * (self.surface - self.ambient);
        let d_surface = (heat_input - q_sb - q_s_conv) / self.cap_surface;

        // Bulk: receives from surface conduction, loses to carcass
        let q_bc = self.cond_bc * (self.bulk - self.carcass);
        let d_bulk = (q_sb - q_bc) / self.cap_bulk;

        // Carcass: receives from bulk, loses to ambient
        let q_c_conv = self.conv_carcass * (self.carcass - self.ambient);
        let d_carcass = (q_bc - q_c_conv) / self.cap_carcass;

        self.surface += d_surface * dt;
        self.bulk += d_bulk * dt;
        self.carcass += d_carcass * dt;
    }

    /// Mean tire temperature (K) — simple average of three layers.
    pub fn mean_temp(&self) -> f64 {
        (self.surface + self.bulk + self.carcass) / 3.0
    }

    /// Temperature gradient surface-to-carcass (K/m) — approximate.
    ///
    /// `tread_depth_m` — tread depth used as distance proxy (m).
    pub fn gradient(&self, tread_depth_m: f64) -> f64 {
        let thickness = tread_depth_m.max(1e-3);
        (self.surface - self.carcass) / thickness
    }
}

// ---------------------------------------------------------------------------
// TreadDepth
// ---------------------------------------------------------------------------

/// Tread depth state and wear-rate estimation.
#[derive(Debug, Clone)]
pub struct TreadDepth {
    /// New-tire tread depth (m).
    pub initial_depth: f64,
    /// Current tread depth (m).
    pub current_depth: f64,
    /// Wear rate (m/km).
    pub wear_rate_per_km: f64,
}

impl TreadDepth {
    /// Construct with standard racing slick starting depth of 3 mm.
    pub fn new_slick() -> Self {
        Self {
            initial_depth: 0.003,
            current_depth: 0.003,
            wear_rate_per_km: 0.000_040,
        }
    }

    /// Construct with road tire starting depth of 8 mm.
    pub fn new_road() -> Self {
        Self {
            initial_depth: 0.008,
            current_depth: 0.008,
            wear_rate_per_km: 0.000_010,
        }
    }

    /// Wear the tire by `distance_km` (km), modifying `current_depth`.
    pub fn wear(&mut self, distance_km: f64) {
        let worn = self.wear_rate_per_km * distance_km;
        self.current_depth = (self.current_depth - worn).max(0.0);
    }

    /// Remaining tread as a fraction of initial depth (0–1).
    pub fn remaining_fraction(&self) -> f64 {
        self.current_depth / self.initial_depth.max(1e-6)
    }

    /// Estimated remaining life in km.
    pub fn remaining_km(&self) -> f64 {
        if self.wear_rate_per_km < 1e-12 {
            f64::INFINITY
        } else {
            self.current_depth / self.wear_rate_per_km
        }
    }

    /// Returns `true` when tread depth is below minimum legal limit (1.6 mm).
    pub fn is_below_legal_limit(&self) -> bool {
        self.current_depth < 0.0016
    }
}

// ---------------------------------------------------------------------------
// TireWearModel (Archard + slip energy)
// ---------------------------------------------------------------------------

/// Combined Archard wear + slip-energy wear model.
#[derive(Debug, Clone)]
pub struct TireWearModel {
    /// Archard wear coefficient K (dimensionless / hardness).
    pub archard_k: f64,
    /// Hardness of the rubber compound (Pa).
    pub hardness: f64,
    /// Slip energy wear coefficient (m/J).
    pub slip_energy_coeff: f64,
    /// Temperature sensitivity exponent.
    pub temp_exponent: f64,
    /// Reference temperature for wear rate (K).
    pub reference_temp: f64,
}

impl TireWearModel {
    /// Default wear model for a medium-compound racing tire.
    pub fn default_medium() -> Self {
        Self {
            archard_k: 1.0e-7,
            hardness: 2.5e6,
            slip_energy_coeff: 5.0e-9,
            temp_exponent: 1.5,
            reference_temp: 363.15,
        }
    }

    /// Archard wear volume (m³) for a given normal load (N) and sliding distance (m).
    pub fn archard_volume(&self, load: f64, sliding_distance: f64) -> f64 {
        self.archard_k * load * sliding_distance / self.hardness
    }

    /// Slip energy wear rate (m/s of tread depth) given slip power (W) and contact area (m²).
    pub fn slip_energy_wear_rate(&self, slip_power: f64, contact_area: f64) -> f64 {
        let area = contact_area.max(1e-6);
        self.slip_energy_coeff * slip_power / area
    }

    /// Temperature factor for wear rate — exponential above reference.
    pub fn temperature_factor(&self, temp_k: f64) -> f64 {
        let dt = (temp_k - self.reference_temp).max(0.0);
        (1.0 + dt / self.reference_temp).powf(self.temp_exponent)
    }

    /// Total wear rate (m/s of tread depth) combining Archard and slip-energy methods.
    pub fn total_wear_rate(
        &self,
        load: f64,
        slip_speed: f64,
        contact_area: f64,
        temp_k: f64,
    ) -> f64 {
        let sliding = slip_speed;
        let area = contact_area.max(1e-6);
        let v_archard = self.archard_k * load * sliding / self.hardness / area;
        let slip_power = load * slip_speed;
        let v_slip = self.slip_energy_wear_rate(slip_power, area);
        let tf = self.temperature_factor(temp_k);
        (v_archard + v_slip) * tf
    }
}

// ---------------------------------------------------------------------------
// GripDegradation
// ---------------------------------------------------------------------------

/// Tire grip degradation due to tread depth reduction and temperature.
#[derive(Debug, Clone)]
pub struct GripDegradation {
    /// Compound type.
    pub compound: CompoundType,
    /// Grip loss per unit tread fraction lost (0–1 scale).
    pub grip_loss_per_fraction: f64,
    /// Temperature sensitivity (grip reduction per K away from optimal).
    pub temp_sensitivity: f64,
}

impl GripDegradation {
    /// Construct for the given compound type with default sensitivity.
    pub fn new(compound: CompoundType) -> Self {
        Self {
            compound,
            grip_loss_per_fraction: 0.15,
            temp_sensitivity: 0.002,
        }
    }

    /// Grip factor (0–1+) given tread fraction remaining and current temperature (K).
    pub fn grip_factor(&self, tread_fraction: f64, temp_k: f64) -> f64 {
        let base = self.compound.grip_coefficient();
        let tread_factor = 1.0 - self.grip_loss_per_fraction * (1.0 - tread_fraction);

        let (_, t_opt, _) = self.compound.optimal_temp_range();
        let temp_delta = (temp_k - t_opt).abs();
        let temp_factor = (1.0 - self.temp_sensitivity * temp_delta).max(0.3);

        base * tread_factor * temp_factor
    }

    /// Returns `true` if temperature is within the optimal window.
    pub fn is_in_optimal_window(&self, temp_k: f64) -> bool {
        let (t_min, _, t_max) = self.compound.optimal_temp_range();
        temp_k >= t_min && temp_k <= t_max
    }
}

// ---------------------------------------------------------------------------
// TireAge
// ---------------------------------------------------------------------------

/// Tire aging model: ozone, UV, and storage effects.
#[derive(Debug, Clone)]
pub struct TireAge {
    /// Age in months from manufacture.
    pub age_months: f64,
    /// Exposure to UV index (0–11+ scale, hours-weighted).
    pub uv_exposure: f64,
    /// Storage temperature (K).
    pub storage_temp: f64,
    /// Ozone concentration (pphm — parts per hundred million).
    pub ozone_pphm: f64,
}

impl TireAge {
    /// Fresh tire from manufacturer.
    pub fn new() -> Self {
        Self {
            age_months: 0.0,
            uv_exposure: 0.0,
            storage_temp: 293.15,
            ozone_pphm: 25.0,
        }
    }

    /// Hardness increase index (0 = new, >1 = significantly degraded).
    pub fn hardness_index(&self) -> f64 {
        let age_factor = self.age_months / 60.0; // 5-year reference
        let uv_factor = self.uv_exposure / 1000.0;
        let ozone_factor = self.ozone_pphm / 100.0;
        let temp_factor = ((self.storage_temp - 293.15) / 20.0).max(0.0);
        age_factor + uv_factor + ozone_factor * 0.5 + temp_factor * 0.3
    }

    /// Estimated grip penalty (0 = no degradation, 1 = fully degraded).
    pub fn grip_penalty(&self) -> f64 {
        (self.hardness_index() * 0.10).min(1.0)
    }

    /// Returns `true` if the tire has exceeded its recommended shelf life (5 years).
    pub fn is_expired(&self) -> bool {
        self.age_months > 60.0
    }
}

impl Default for TireAge {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// WearPattern
// ---------------------------------------------------------------------------

/// Classification of tire wear patterns for diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WearPattern {
    /// Uniform wear — ideal condition.
    Uniform,
    /// Centre wear — overinflation.
    CentreWear,
    /// Edge wear (both shoulders) — underinflation.
    EdgeWear,
    /// One-sided wear — excess camber.
    OneSideWear,
    /// Cupping/scalloping — worn shock absorbers.
    Cupping,
    /// Flat spot — locked-wheel braking.
    FlatSpot,
}

impl WearPattern {
    /// Diagnosis string for this pattern.
    pub fn diagnosis(self) -> &'static str {
        match self {
            WearPattern::Uniform => "Good condition — uniform wear",
            WearPattern::CentreWear => "Overinflation — reduce pressure",
            WearPattern::EdgeWear => "Underinflation — increase pressure",
            WearPattern::OneSideWear => "Excess camber — check alignment",
            WearPattern::Cupping => "Worn dampers — replace shocks",
            WearPattern::FlatSpot => "Locked-wheel brake event — inspect",
        }
    }

    /// Estimate the pattern from centre/edge tread ratios and lateral imbalance.
    ///
    /// `centre_depth` — tread depth at crown centre (m).
    /// `edge_depth` — tread depth at shoulder (m).
    /// `inner_depth` — tread depth at inner shoulder (m).
    /// `outer_depth` — tread depth at outer shoulder (m).
    pub fn diagnose(
        centre_depth: f64,
        edge_depth: f64,
        inner_depth: f64,
        outer_depth: f64,
    ) -> Self {
        let lateral_diff = (inner_depth - outer_depth).abs();
        let centre_edge_diff = centre_depth - edge_depth;

        if lateral_diff > 0.001 {
            return WearPattern::OneSideWear;
        }
        // centre_depth < edge_depth → centre is more worn → CentreWear
        if centre_edge_diff < -0.001 {
            return WearPattern::CentreWear;
        }
        // centre_depth > edge_depth → edge is more worn → EdgeWear
        if centre_edge_diff > 0.001 {
            return WearPattern::EdgeWear;
        }
        WearPattern::Uniform
    }
}

// ---------------------------------------------------------------------------
// TirePressureModel
// ---------------------------------------------------------------------------

/// Tire pressure model — cold pressure, thermal rise, blowout threshold.
#[derive(Debug, Clone)]
pub struct TirePressureModel {
    /// Cold pressure at ambient temperature (Pa).
    pub cold_pressure: f64,
    /// Current pressure (Pa).
    pub current_pressure: f64,
    /// Ambient temperature (K).
    pub ambient_temp: f64,
    /// Blowout pressure threshold (Pa).
    pub blowout_threshold: f64,
    /// Pressure sensitivity to temperature (Pa/K) — Gay-Lussac approximation.
    pub pressure_temp_coeff: f64,
}

impl TirePressureModel {
    /// Default F1-style slick pressure setup.
    pub fn default_f1_slick() -> Self {
        let cold_psi = 21.0; // psi
        let cold_pa = cold_psi * 6_894.76;
        Self {
            cold_pressure: cold_pa,
            current_pressure: cold_pa,
            ambient_temp: 293.15,
            blowout_threshold: 350_000.0,
            pressure_temp_coeff: cold_pa / 293.15,
        }
    }

    /// Default road car pressure (2.2 bar).
    pub fn default_road() -> Self {
        Self {
            cold_pressure: 220_000.0,
            current_pressure: 220_000.0,
            ambient_temp: 293.15,
            blowout_threshold: 550_000.0,
            pressure_temp_coeff: 220_000.0 / 293.15,
        }
    }

    /// Update pressure based on current tire temperature (K).
    pub fn update_pressure(&mut self, temp_k: f64) {
        // Gay-Lussac's law: P/T = const
        self.current_pressure = self.pressure_temp_coeff * temp_k;
    }

    /// Pressure in bar.
    pub fn pressure_bar(&self) -> f64 {
        self.current_pressure / 100_000.0
    }

    /// Returns `true` if pressure exceeds blowout threshold.
    pub fn is_blowout_risk(&self) -> bool {
        self.current_pressure >= self.blowout_threshold
    }

    /// Pressure rise above cold pressure (Pa).
    pub fn pressure_rise(&self) -> f64 {
        self.current_pressure - self.cold_pressure
    }
}

// ---------------------------------------------------------------------------
// TireForceHistory
// ---------------------------------------------------------------------------

/// Rolling contact patch analysis and stress history.
#[derive(Debug, Clone)]
pub struct TireForceHistory {
    /// Contact patch half-length (m).
    pub half_length: f64,
    /// Contact patch half-width (m).
    pub half_width: f64,
    /// Normal load (N).
    pub normal_load: f64,
    /// Crown radius (m) — changes with wear.
    pub crown_radius: f64,
    /// Accumulated vertical force integral over distance (N·m).
    pub load_integral: f64,
    /// Accumulated lateral force integral (N·m).
    pub lateral_integral: f64,
}

impl TireForceHistory {
    /// Construct for a typical F1 slick tire.
    pub fn new_f1_slick(load: f64) -> Self {
        Self {
            half_length: 0.10,
            half_width: 0.155,
            normal_load: load,
            crown_radius: 0.330,
            load_integral: 0.0,
            lateral_integral: 0.0,
        }
    }

    /// Contact patch area (m²).
    pub fn contact_area(&self) -> f64 {
        PI * self.half_length * self.half_width
    }

    /// Average contact pressure (Pa).
    pub fn mean_contact_pressure(&self) -> f64 {
        self.normal_load / self.contact_area().max(1e-6)
    }

    /// Accumulate force history over a distance step `ds` (m).
    pub fn accumulate(&mut self, lateral_force: f64, ds: f64) {
        self.load_integral += self.normal_load * ds;
        self.lateral_integral += lateral_force.abs() * ds;
    }

    /// Update crown radius as tread wears — tread removal reduces crown radius slightly.
    pub fn update_crown_radius(&mut self, wear_depth: f64) {
        self.crown_radius = (self.crown_radius - wear_depth * 0.10).max(0.200);
    }

    /// Reset accumulated force integrals.
    pub fn reset_history(&mut self) {
        self.load_integral = 0.0;
        self.lateral_integral = 0.0;
    }
}

// ---------------------------------------------------------------------------
// TireWearIncrement
// ---------------------------------------------------------------------------

/// Output of a single wear simulation step.
#[derive(Debug, Clone)]
pub struct TireWearIncrement {
    /// Tread depth reduction in this step (m).
    pub depth_reduction: f64,
    /// Heat generated (J).
    pub heat_generated: f64,
    /// Grip factor after wear.
    pub grip_factor: f64,
    /// Estimated remaining life (km).
    pub remaining_km: f64,
}

// ---------------------------------------------------------------------------
// WearSimulation
// ---------------------------------------------------------------------------

/// Lap-level wear simulation combining all tire sub-models.
#[derive(Debug, Clone)]
pub struct WearSimulation {
    /// Tread depth state.
    pub tread: TreadDepth,
    /// Thermal state.
    pub temperature: TireTemperature,
    /// Wear model.
    pub wear_model: TireWearModel,
    /// Grip degradation model.
    pub grip_degradation: GripDegradation,
    /// Pressure model.
    pub pressure: TirePressureModel,
    /// Force history.
    pub force_history: TireForceHistory,
}

impl WearSimulation {
    /// Construct a default medium-compound F1-style simulation.
    pub fn new_medium_compound() -> Self {
        Self {
            tread: TreadDepth::new_slick(),
            temperature: TireTemperature::default_medium(),
            wear_model: TireWearModel::default_medium(),
            grip_degradation: GripDegradation::new(CompoundType::Medium),
            pressure: TirePressureModel::default_f1_slick(),
            force_history: TireForceHistory::new_f1_slick(4_000.0),
        }
    }

    /// Simulate one track segment step.
    ///
    /// # Arguments
    /// * `load` — vertical tire load (N)
    /// * `slip_speed` — lateral slip speed (m/s)
    /// * `lateral_force` — lateral force (N)
    /// * `distance_m` — segment distance (m)
    /// * `dt` — time step (s)
    pub fn step(
        &mut self,
        load: f64,
        slip_speed: f64,
        lateral_force: f64,
        distance_m: f64,
        dt: f64,
    ) -> TireWearIncrement {
        let area = self.force_history.contact_area();

        // Heat generation from slip
        let heat_input = load * slip_speed.abs();
        self.temperature.step(dt, heat_input);

        // Wear rate
        let wear_rate =
            self.wear_model
                .total_wear_rate(load, slip_speed.abs(), area, self.temperature.surface);
        let depth_reduction = wear_rate * dt;
        self.tread.current_depth = (self.tread.current_depth - depth_reduction).max(0.0);

        // Update pressure from temperature
        self.pressure.update_pressure(self.temperature.bulk);

        // Accumulate force history
        self.force_history.accumulate(lateral_force, distance_m);

        // Update crown radius
        self.force_history.update_crown_radius(depth_reduction);

        // Compute grip
        let tread_frac = self.tread.remaining_fraction();
        let grip = self
            .grip_degradation
            .grip_factor(tread_frac, self.temperature.surface);

        TireWearIncrement {
            depth_reduction,
            heat_generated: heat_input * dt,
            grip_factor: grip,
            remaining_km: self.tread.remaining_km(),
        }
    }

    /// Total distance driven in this session (km) inferred from tread loss.
    pub fn distance_driven_km(&self) -> f64 {
        let worn = self.tread.initial_depth - self.tread.current_depth;
        worn / self.wear_model.slip_energy_coeff.max(1e-12) * 0.001
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CompoundType ---

    #[test]
    fn compound_soft_highest_grip() {
        assert!(CompoundType::Soft.grip_coefficient() > CompoundType::Medium.grip_coefficient());
        assert!(CompoundType::Medium.grip_coefficient() > CompoundType::Hard.grip_coefficient());
    }

    #[test]
    fn compound_hard_highest_durability() {
        assert!(CompoundType::Hard.durability_factor() > CompoundType::Medium.durability_factor());
    }

    #[test]
    fn compound_is_dry() {
        assert!(CompoundType::Soft.is_dry());
        assert!(!CompoundType::Wet.is_dry());
    }

    #[test]
    fn compound_optimal_temp_range_ordered() {
        let (min, peak, max) = CompoundType::Medium.optimal_temp_range();
        assert!(min < peak && peak < max);
    }

    // --- TireTemperature ---

    #[test]
    fn tire_temp_starts_at_ambient() {
        let t = TireTemperature::new_at_ambient(293.15);
        assert!((t.surface - 293.15).abs() < 1e-10);
        assert!((t.bulk - 293.15).abs() < 1e-10);
        assert!((t.carcass - 293.15).abs() < 1e-10);
    }

    #[test]
    fn tire_temp_heats_with_slip() {
        let mut t = TireTemperature::new_at_ambient(293.15);
        for _ in 0..100 {
            t.step(0.01, 5000.0); // large heat input
        }
        assert!(t.surface > 293.15, "surface should heat up: {}", t.surface);
    }

    #[test]
    fn tire_temp_mean() {
        let mut t = TireTemperature::new_at_ambient(300.0);
        t.surface = 320.0;
        t.bulk = 310.0;
        t.carcass = 300.0;
        let mean = t.mean_temp();
        assert!((mean - 310.0).abs() < 1e-10);
    }

    #[test]
    fn tire_temp_gradient() {
        let mut t = TireTemperature::new_at_ambient(293.15);
        t.surface = 373.15;
        t.carcass = 313.15;
        let grad = t.gradient(0.003);
        assert!(grad > 0.0, "gradient should be positive: {grad}");
    }

    // --- TreadDepth ---

    #[test]
    fn tread_depth_wears_over_distance() {
        let mut td = TreadDepth::new_slick();
        let initial = td.current_depth;
        td.wear(50.0); // 50 km
        assert!(td.current_depth < initial);
    }

    #[test]
    fn tread_depth_remaining_fraction() {
        let mut td = TreadDepth::new_slick();
        td.wear(10.0);
        let frac = td.remaining_fraction();
        assert!(frac < 1.0 && frac > 0.0, "fraction should be < 1: {frac}");
    }

    #[test]
    fn tread_depth_remaining_km() {
        let td = TreadDepth::new_slick();
        let km = td.remaining_km();
        assert!(km > 0.0 && km < 1_000.0, "remaining km out of range: {km}");
    }

    #[test]
    fn tread_depth_legal_limit() {
        let mut td = TreadDepth::new_road();
        td.current_depth = 0.001; // below 1.6 mm
        assert!(td.is_below_legal_limit());
    }

    #[test]
    fn tread_depth_not_below_zero() {
        let mut td = TreadDepth::new_slick();
        td.wear(1_000_000.0); // extreme wear
        assert_eq!(td.current_depth, 0.0);
    }

    // --- TireWearModel ---

    #[test]
    fn archard_volume_increases_with_load() {
        let model = TireWearModel::default_medium();
        let v1 = model.archard_volume(1000.0, 1.0);
        let v2 = model.archard_volume(2000.0, 1.0);
        assert!(v2 > v1, "higher load → more wear: {v1} vs {v2}");
    }

    #[test]
    fn wear_model_temperature_factor_above_reference() {
        let model = TireWearModel::default_medium();
        let tf = model.temperature_factor(model.reference_temp + 20.0);
        assert!(tf > 1.0, "wear should increase above reference temp: {tf}");
    }

    #[test]
    fn wear_model_temperature_factor_at_reference() {
        let model = TireWearModel::default_medium();
        let tf = model.temperature_factor(model.reference_temp);
        assert!(
            (tf - 1.0).abs() < 1e-10,
            "factor should be 1 at reference: {tf}"
        );
    }

    #[test]
    fn total_wear_rate_positive() {
        let model = TireWearModel::default_medium();
        let rate = model.total_wear_rate(4000.0, 2.0, 0.02, 370.0);
        assert!(rate > 0.0, "wear rate should be positive: {rate}");
    }

    // --- GripDegradation ---

    #[test]
    fn grip_at_full_tread_near_nominal() {
        let g = GripDegradation::new(CompoundType::Medium);
        let grip = g.grip_factor(1.0, 363.15);
        let nominal = CompoundType::Medium.grip_coefficient();
        assert!(
            (grip - nominal).abs() < 0.01,
            "grip at full tread should match nominal: {grip}"
        );
    }

    #[test]
    fn grip_reduces_with_tread_loss() {
        let g = GripDegradation::new(CompoundType::Medium);
        let g1 = g.grip_factor(1.0, 363.15);
        let g2 = g.grip_factor(0.5, 363.15);
        assert!(g2 < g1, "grip should reduce with tread loss: {g1} vs {g2}");
    }

    #[test]
    fn grip_in_optimal_window() {
        let g = GripDegradation::new(CompoundType::Medium);
        let (min, peak, _) = CompoundType::Medium.optimal_temp_range();
        assert!(g.is_in_optimal_window(peak));
        assert!(!g.is_in_optimal_window(min - 20.0));
    }

    // --- TireAge ---

    #[test]
    fn new_tire_no_hardness() {
        let age = TireAge::new();
        assert!(
            age.hardness_index() < 0.5,
            "new tire should have low hardness index"
        );
    }

    #[test]
    fn old_tire_higher_hardness() {
        let mut age = TireAge::new();
        age.age_months = 72.0; // 6 years
        assert!(age.hardness_index() > 1.0);
    }

    #[test]
    fn tire_age_expiry() {
        let mut age = TireAge::new();
        assert!(!age.is_expired());
        age.age_months = 61.0;
        assert!(age.is_expired());
    }

    // --- WearPattern ---

    #[test]
    fn wear_pattern_uniform() {
        let p = WearPattern::diagnose(0.003, 0.003, 0.003, 0.003);
        assert_eq!(p, WearPattern::Uniform);
    }

    #[test]
    fn wear_pattern_centre_wear() {
        let p = WearPattern::diagnose(0.002, 0.004, 0.004, 0.004);
        assert_eq!(p, WearPattern::CentreWear);
    }

    #[test]
    fn wear_pattern_edge_wear() {
        let p = WearPattern::diagnose(0.004, 0.002, 0.002, 0.002);
        assert_eq!(p, WearPattern::EdgeWear);
    }

    #[test]
    fn wear_pattern_one_side() {
        let p = WearPattern::diagnose(0.003, 0.003, 0.001, 0.003);
        assert_eq!(p, WearPattern::OneSideWear);
    }

    #[test]
    fn wear_pattern_diagnosis_non_empty() {
        assert!(!WearPattern::Uniform.diagnosis().is_empty());
        assert!(!WearPattern::Cupping.diagnosis().is_empty());
    }

    // --- TirePressureModel ---

    #[test]
    fn pressure_rises_with_temperature() {
        let mut pm = TirePressureModel::default_road();
        let p0 = pm.current_pressure;
        pm.update_pressure(pm.ambient_temp + 30.0);
        assert!(
            pm.current_pressure > p0,
            "pressure should rise: {p0} vs {}",
            pm.current_pressure
        );
    }

    #[test]
    fn pressure_bar_in_range() {
        let pm = TirePressureModel::default_road();
        let bar = pm.pressure_bar();
        assert!(bar > 1.5 && bar < 4.0, "pressure out of range: {bar} bar");
    }

    #[test]
    fn blowout_risk_at_high_pressure() {
        let mut pm = TirePressureModel::default_road();
        pm.current_pressure = pm.blowout_threshold + 1.0;
        assert!(pm.is_blowout_risk());
    }

    // --- TireForceHistory ---

    #[test]
    fn contact_area_positive() {
        let fh = TireForceHistory::new_f1_slick(4000.0);
        assert!(fh.contact_area() > 0.0);
    }

    #[test]
    fn mean_contact_pressure_reasonable() {
        let fh = TireForceHistory::new_f1_slick(4000.0);
        let p = fh.mean_contact_pressure();
        // Typical range 100–400 kPa
        assert!(
            p > 50_000.0 && p < 1_000_000.0,
            "pressure out of range: {p}"
        );
    }

    #[test]
    fn force_history_accumulates() {
        let mut fh = TireForceHistory::new_f1_slick(4000.0);
        fh.accumulate(1000.0, 100.0);
        assert!(fh.load_integral > 0.0);
        assert!(fh.lateral_integral > 0.0);
    }

    #[test]
    fn force_history_reset() {
        let mut fh = TireForceHistory::new_f1_slick(4000.0);
        fh.accumulate(1000.0, 100.0);
        fh.reset_history();
        assert_eq!(fh.load_integral, 0.0);
        assert_eq!(fh.lateral_integral, 0.0);
    }

    // --- WearSimulation ---

    #[test]
    fn wear_simulation_step_reduces_tread() {
        let mut sim = WearSimulation::new_medium_compound();
        let initial = sim.tread.current_depth;
        for _ in 0..100 {
            sim.step(4000.0, 2.0, 800.0, 10.0, 0.1);
        }
        assert!(sim.tread.current_depth < initial, "tread should reduce");
    }

    #[test]
    fn wear_simulation_step_grip_factor_positive() {
        let mut sim = WearSimulation::new_medium_compound();
        let inc = sim.step(4000.0, 2.0, 800.0, 10.0, 0.1);
        assert!(
            inc.grip_factor > 0.0,
            "grip factor should be positive: {}",
            inc.grip_factor
        );
    }

    #[test]
    fn wear_simulation_step_heat_generated_positive() {
        let mut sim = WearSimulation::new_medium_compound();
        let inc = sim.step(4000.0, 2.0, 800.0, 10.0, 0.1);
        assert!(inc.heat_generated > 0.0);
    }

    #[test]
    fn wear_simulation_pressure_updates() {
        let mut sim = WearSimulation::new_medium_compound();
        // Heat up tire to change pressure
        for _ in 0..500 {
            sim.step(4000.0, 2.0, 800.0, 10.0, 0.1);
        }
        assert!(sim.pressure.current_pressure != sim.pressure.cold_pressure);
    }
}
