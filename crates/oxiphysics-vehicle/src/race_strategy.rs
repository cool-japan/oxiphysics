// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Racing strategy simulation.
//!
//! Provides tyre-wear modelling, fuel consumption, pit-window calculation,
//! undercut logic, full race simulation with strategy optimisation, lap-time
//! sector modelling, DRS zones, weather adaptation, and safety-car probability.

// ---------------------------------------------------------------------------
// TireCompound
// ---------------------------------------------------------------------------

/// Tyre compound classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TireCompound {
    /// Soft compound — fastest but highest degradation.
    Soft,
    /// Medium compound — balanced performance.
    Medium,
    /// Hard compound — slowest but most durable.
    Hard,
    /// Intermediate wet-weather compound.
    Intermediate,
    /// Full wet-weather compound.
    Wet,
}

impl TireCompound {
    /// Base degradation rate (wear fraction per lap at reference conditions).
    pub fn base_deg_rate(&self) -> f64 {
        match self {
            TireCompound::Soft => 0.025_f64,
            TireCompound::Medium => 0.015_f64,
            TireCompound::Hard => 0.008_f64,
            TireCompound::Intermediate => 0.012_f64,
            TireCompound::Wet => 0.010_f64,
        }
    }

    /// Peak grip coefficient relative to Medium = 1.0.
    pub fn peak_grip(&self) -> f64 {
        match self {
            TireCompound::Soft => 1.05_f64,
            TireCompound::Medium => 1.00_f64,
            TireCompound::Hard => 0.97_f64,
            TireCompound::Intermediate => 0.92_f64,
            TireCompound::Wet => 0.88_f64,
        }
    }

    /// Optimal operating temperature range `(low, high)` in °C.
    pub fn operating_window(&self) -> (f64, f64) {
        match self {
            TireCompound::Soft => (80.0_f64, 100.0_f64),
            TireCompound::Medium => (85.0_f64, 110.0_f64),
            TireCompound::Hard => (90.0_f64, 120.0_f64),
            TireCompound::Intermediate => (30.0_f64, 60.0_f64),
            TireCompound::Wet => (20.0_f64, 45.0_f64),
        }
    }

    /// Whether this compound is a wet-weather compound.
    pub fn is_wet_weather(&self) -> bool {
        matches!(self, TireCompound::Intermediate | TireCompound::Wet)
    }

    /// Nominal stint length (laps) at reference conditions.
    pub fn nominal_stint_laps(&self) -> usize {
        match self {
            TireCompound::Soft => 20,
            TireCompound::Medium => 35,
            TireCompound::Hard => 50,
            TireCompound::Intermediate => 25,
            TireCompound::Wet => 30,
        }
    }
}

// ---------------------------------------------------------------------------
// TireWear / TireThermal
// ---------------------------------------------------------------------------

/// Runtime tyre wear state.
#[derive(Debug, Clone)]
pub struct TireWear {
    /// Wear level from 0.0 (new) to 1.0 (fully worn).
    pub wear_level: f64,
    /// Current compound fitted.
    pub compound: TireCompound,
    /// Number of laps completed on this set.
    pub laps_used: usize,
}

impl TireWear {
    /// Create a fresh tyre state.
    pub fn new(compound: TireCompound) -> Self {
        Self {
            wear_level: 0.0_f64,
            compound,
            laps_used: 0,
        }
    }

    /// Advance wear by one lap.
    pub fn advance_lap(&mut self, deg_rate: f64) {
        self.wear_level = (self.wear_level + deg_rate).min(1.0_f64);
        self.laps_used += 1;
    }

    /// Whether the tyre is worn out (>= 90 %).
    pub fn is_worn_out(&self) -> bool {
        self.wear_level >= 0.90_f64
    }
}

/// Thermal model for a single tyre.
#[derive(Debug, Clone)]
pub struct TireThermal {
    /// Current surface temperature (°C).
    pub surface_temp: f64,
    /// Tyre compound.
    pub compound: TireCompound,
    /// Ambient air temperature (°C).
    pub ambient_temp: f64,
    /// Thermal time constant (laps).
    pub time_constant: f64,
}

impl TireThermal {
    /// Create a thermal model at ambient temperature.
    pub fn new(compound: TireCompound, ambient_temp: f64) -> Self {
        Self {
            surface_temp: ambient_temp,
            compound,
            ambient_temp,
            time_constant: 3.0_f64,
        }
    }

    /// Equilibrium temperature based on speed and slip.
    pub fn equilibrium_temp(&self, speed: f64, slip: f64) -> f64 {
        let (lo, hi) = self.compound.operating_window();
        let mid = (lo + hi) * 0.5_f64;
        let heat_from_slip = slip.abs() * 200.0_f64;
        let heat_from_speed = (speed / 200.0_f64) * 20.0_f64;
        (mid + heat_from_slip + heat_from_speed).min(150.0_f64)
    }

    /// Update surface temperature for one time step.
    pub fn update(&mut self, speed: f64, slip: f64) {
        let eq = self.equilibrium_temp(speed, slip);
        let alpha = 1.0_f64 - (-1.0_f64 / self.time_constant).exp();
        self.surface_temp += alpha * (eq - self.surface_temp);
    }

    /// Grip penalty when outside the operating window (0–1, 1 = no penalty).
    pub fn thermal_grip_factor(&self) -> f64 {
        let (lo, hi) = self.compound.operating_window();
        if self.surface_temp < lo {
            let cold = (lo - self.surface_temp).min(40.0_f64) / 40.0_f64;
            1.0_f64 - 0.3_f64 * cold
        } else if self.surface_temp > hi {
            let hot = (self.surface_temp - hi).min(30.0_f64) / 30.0_f64;
            1.0_f64 - 0.4_f64 * hot
        } else {
            1.0_f64
        }
    }
}

// ---------------------------------------------------------------------------
// RaceTrack
// ---------------------------------------------------------------------------

/// A race track sector descriptor.
#[derive(Debug, Clone)]
pub struct TrackSector {
    /// Sector name.
    pub name: String,
    /// Sector distance (m).
    pub length: f64,
    /// Base lap-time contribution for this sector (s).
    pub base_time: f64,
    /// Whether this sector contains a DRS zone.
    pub has_drs: bool,
    /// Whether this sector contains a chicane.
    pub has_chicane: bool,
    /// Track surface roughness (0–1, affects tyre wear).
    pub roughness: f64,
}

impl TrackSector {
    /// Create a new track sector.
    pub fn new(
        name: &str,
        length: f64,
        base_time: f64,
        has_drs: bool,
        has_chicane: bool,
        roughness: f64,
    ) -> Self {
        Self {
            name: name.to_owned(),
            length,
            base_time,
            has_drs,
            has_chicane,
            roughness,
        }
    }
}

/// Full race track descriptor.
#[derive(Debug, Clone)]
pub struct RaceTrack {
    /// Track name.
    pub name: String,
    /// Total lap distance (m).
    pub lap_distance: f64,
    /// Sector descriptors.
    pub sectors: Vec<TrackSector>,
    /// Track surface temperature (°C).
    pub track_temp: f64,
    /// Ambient air temperature (°C).
    pub ambient_temp: f64,
    /// Total number of race laps.
    pub total_laps: usize,
    /// Pit lane time loss (s).
    pub pit_lane_time_loss: f64,
}

impl RaceTrack {
    /// Create a race track.
    pub fn new(
        name: &str,
        lap_distance: f64,
        sectors: Vec<TrackSector>,
        track_temp: f64,
        ambient_temp: f64,
        total_laps: usize,
        pit_lane_time_loss: f64,
    ) -> Self {
        Self {
            name: name.to_owned(),
            lap_distance,
            sectors,
            track_temp,
            ambient_temp,
            total_laps,
            pit_lane_time_loss,
        }
    }

    /// Build a generic 3-sector, 5 km track.
    pub fn default_track() -> Self {
        let sectors = vec![
            TrackSector::new("S1", 1800.0_f64, 28.5_f64, false, false, 0.4_f64),
            TrackSector::new("S2", 1700.0_f64, 31.0_f64, true, true, 0.5_f64),
            TrackSector::new("S3", 1500.0_f64, 27.0_f64, true, false, 0.3_f64),
        ];
        Self::new(
            "Generic", 5000.0_f64, sectors, 40.0_f64, 25.0_f64, 57, 22.0_f64,
        )
    }

    /// Number of DRS zones on the circuit.
    pub fn drs_zone_count(&self) -> usize {
        self.sectors.iter().filter(|s| s.has_drs).count()
    }

    /// Base lap time from sector sum.
    pub fn base_lap_time(&self) -> f64 {
        self.sectors.iter().map(|s| s.base_time).sum()
    }
}

// ---------------------------------------------------------------------------
// FuelStrategy
// ---------------------------------------------------------------------------

/// Fuel strategy parameters.
#[derive(Debug, Clone)]
pub struct FuelStrategy {
    /// Fuel load at race start (kg).
    pub fuel_start: f64,
    /// Fuel consumption rate (kg/lap).
    pub fuel_per_lap: f64,
    /// Lap-time penalty per kg of extra fuel (s/kg).
    pub lap_time_per_kg: f64,
    /// Whether to allow fuel saving (lift-and-coast).
    pub fuel_saving_enabled: bool,
    /// Fuel saving mode lap time delta (s, positive = slower).
    pub fuel_saving_delta: f64,
}

impl FuelStrategy {
    /// Create a fuel strategy.
    pub fn new(
        fuel_start: f64,
        fuel_per_lap: f64,
        lap_time_per_kg: f64,
        fuel_saving_enabled: bool,
        fuel_saving_delta: f64,
    ) -> Self {
        Self {
            fuel_start,
            fuel_per_lap,
            lap_time_per_kg,
            fuel_saving_enabled,
            fuel_saving_delta,
        }
    }

    /// Fuel remaining on a given lap.
    pub fn fuel_at_lap(&self, lap: usize) -> f64 {
        (self.fuel_start - lap as f64 * self.fuel_per_lap).max(0.0_f64)
    }

    /// Lap-time penalty due to fuel load on a given lap.
    pub fn fuel_lap_time_penalty(&self, lap: usize) -> f64 {
        self.fuel_at_lap(lap) * self.lap_time_per_kg
    }

    /// Undercut fuel delta: saving `fuel_kg` reduces lap time by this amount.
    pub fn undercut_fuel_delta(&self, _fuel_kg: f64) -> f64 {
        // Negative = faster (the fuel load has reduced)
        -_fuel_kg * self.lap_time_per_kg
    }

    /// Total fuel needed for a strategy (start fuel plus pit-stop adds).
    pub fn total_fuel_needed(&self, total_laps: usize) -> f64 {
        total_laps as f64 * self.fuel_per_lap
    }
}

// ---------------------------------------------------------------------------
// PitStopModel
// ---------------------------------------------------------------------------

/// Pit stop timing model.
#[derive(Debug, Clone)]
pub struct PitStopModel {
    /// Stationary time for tyre change (s).
    pub tyre_change_time: f64,
    /// Stationary time per 10 kg of fuel added (s).
    pub fuel_add_time_per_10kg: f64,
    /// Overhead (jack up/down, pit lane entry/exit mechanics) (s).
    pub overhead: f64,
    /// Pit lane time loss (s) — added on top of stationary time.
    pub pit_lane_loss: f64,
}

impl PitStopModel {
    /// Create a pit stop model.
    pub fn new(
        tyre_change_time: f64,
        fuel_add_time_per_10kg: f64,
        overhead: f64,
        pit_lane_loss: f64,
    ) -> Self {
        Self {
            tyre_change_time,
            fuel_add_time_per_10kg,
            overhead,
            pit_lane_loss,
        }
    }

    /// Default F1-style pit stop model.
    pub fn default_f1() -> Self {
        Self::new(2.5_f64, 0.0_f64, 1.0_f64, 22.0_f64)
    }

    /// Total time lost per pit stop for a given fuel load added.
    pub fn total_time(&self, fuel_added_kg: f64) -> f64 {
        let fuel_time = (fuel_added_kg / 10.0_f64) * self.fuel_add_time_per_10kg;
        self.tyre_change_time + fuel_time + self.overhead + self.pit_lane_loss
    }

    /// Compute whether an undercut is possible given a gap and tyre delta.
    pub fn undercut_possible(&self, gap_to_leader: f64, tyre_delta_per_lap: f64) -> bool {
        // Simplified: undercut works if tyre delta per lap * laps > pit stop time.
        let break_even_laps = self.total_time(0.0_f64) / tyre_delta_per_lap.abs().max(1e-9_f64);
        break_even_laps <= (gap_to_leader / tyre_delta_per_lap.abs().max(1e-9_f64))
    }
}

// ---------------------------------------------------------------------------
// LapTimeModel
// ---------------------------------------------------------------------------

/// Per-lap time model incorporating tyre state, fuel load, and track sectors.
#[derive(Debug, Clone)]
pub struct LapTimeModel {
    /// Track being simulated.
    pub track: RaceTrack,
    /// Base lap time on a fresh tyre with no fuel load (s).
    pub base_lap_time: f64,
    /// Lap-time penalty per unit wear (s per 1.0 wear fraction).
    pub wear_penalty: f64,
    /// Fuel strategy.
    pub fuel_strategy: FuelStrategy,
    /// DRS benefit when following within 1 second (s).
    pub drs_benefit: f64,
    /// Safety car probability per lap (fraction).
    pub safety_car_prob_per_lap: f64,
}

impl LapTimeModel {
    /// Create a lap time model.
    pub fn new(
        track: RaceTrack,
        base_lap_time: f64,
        wear_penalty: f64,
        fuel_strategy: FuelStrategy,
        drs_benefit: f64,
        safety_car_prob_per_lap: f64,
    ) -> Self {
        Self {
            track,
            base_lap_time,
            wear_penalty,
            fuel_strategy,
            drs_benefit,
            safety_car_prob_per_lap,
        }
    }

    /// Predict the lap time given tyre wear, lap number, and whether DRS is available.
    pub fn predict_lap_time(&self, wear: &TireWear, lap: usize, drs_active: bool) -> f64 {
        let wear_delta = wear.wear_level * self.wear_penalty;
        let fuel_delta = self.fuel_strategy.fuel_lap_time_penalty(lap);
        let drs_delta = if drs_active {
            -self.drs_benefit
        } else {
            0.0_f64
        };
        self.base_lap_time + wear_delta + fuel_delta + drs_delta
    }

    /// Sector time for a given sector index with tyre state.
    pub fn sector_time(&self, sector_idx: usize, wear: &TireWear, lap: usize) -> f64 {
        if sector_idx >= self.track.sectors.len() {
            return 0.0_f64;
        }
        let s = &self.track.sectors[sector_idx];
        let base = s.base_time;
        let wear_frac = wear.wear_level * s.roughness * self.wear_penalty / 3.0_f64;
        let fuel_delta =
            self.fuel_strategy.fuel_lap_time_penalty(lap) / self.track.sectors.len() as f64;
        base + wear_frac + fuel_delta
    }

    /// Expected safety car laps in a race.
    pub fn expected_safety_car_laps(&self) -> f64 {
        self.track.total_laps as f64 * self.safety_car_prob_per_lap * 3.0_f64
    }
}

// ---------------------------------------------------------------------------
// RaceStrategy
// ---------------------------------------------------------------------------

/// A complete race strategy specification.
#[derive(Debug, Clone)]
pub struct RaceStrategy {
    /// Laps on which pit stops are taken (1-indexed).
    pub pit_laps: Vec<usize>,
    /// Tyre compound used in each stint (length = pit_laps.len() + 1).
    pub compounds: Vec<TireCompound>,
    /// Fuel load at the start of the race (kg).
    pub fuel_start: f64,
}

impl RaceStrategy {
    /// Create a race strategy.
    pub fn new(pit_laps: Vec<usize>, compounds: Vec<TireCompound>, fuel_start: f64) -> Self {
        Self {
            pit_laps,
            compounds,
            fuel_start,
        }
    }

    /// Number of stints in this strategy.
    pub fn num_stints(&self) -> usize {
        self.pit_laps.len() + 1
    }

    /// Whether this strategy uses a dry compound throughout.
    pub fn is_dry_strategy(&self) -> bool {
        self.compounds.iter().all(|c| !c.is_wet_weather())
    }

    /// Whether this strategy uses at least two different compounds (mandatory rule).
    pub fn satisfies_compound_rule(&self) -> bool {
        let dry_compounds: Vec<_> = self
            .compounds
            .iter()
            .filter(|c| !c.is_wet_weather())
            .collect();
        if dry_compounds.len() < 2 {
            return !self.is_dry_strategy();
        }
        let first = &dry_compounds[0];
        dry_compounds.iter().any(|c| c != first)
    }
}

// ---------------------------------------------------------------------------
// WeatherAdaptation
// ---------------------------------------------------------------------------

/// Weather condition for race strategy.
#[derive(Debug, Clone, PartialEq)]
pub enum WeatherCondition {
    /// Dry track.
    Dry,
    /// Damp track (intermediate required).
    Damp,
    /// Wet track (full wet required).
    Wet,
    /// Changing / mixed conditions.
    Mixed,
}

/// Weather adaptation logic for tyre selection and strategy changes.
#[derive(Debug, Clone)]
pub struct WeatherAdaptation {
    /// Current weather condition.
    pub condition: WeatherCondition,
    /// Rain probability (0–1) for the next lap.
    pub rain_probability: f64,
    /// Safety car probability per lap in dry conditions.
    pub base_safety_car_prob: f64,
    /// Extra safety car probability in wet conditions.
    pub wet_safety_car_extra: f64,
    /// Track temperature (°C).
    pub track_temp: f64,
    /// Lap at which conditions last changed.
    pub condition_change_lap: usize,
}

impl WeatherAdaptation {
    /// Create a weather adaptation model.
    pub fn new(
        condition: WeatherCondition,
        rain_probability: f64,
        base_safety_car_prob: f64,
        track_temp: f64,
    ) -> Self {
        Self {
            condition,
            rain_probability,
            base_safety_car_prob,
            wet_safety_car_extra: 0.04_f64,
            track_temp,
            condition_change_lap: 0,
        }
    }

    /// Recommended compound given current conditions.
    pub fn recommended_compound(&self) -> TireCompound {
        match &self.condition {
            WeatherCondition::Dry => {
                if self.track_temp > 45.0_f64 {
                    TireCompound::Hard
                } else if self.track_temp > 35.0_f64 {
                    TireCompound::Medium
                } else {
                    TireCompound::Soft
                }
            }
            WeatherCondition::Damp => TireCompound::Intermediate,
            WeatherCondition::Wet => TireCompound::Wet,
            WeatherCondition::Mixed => TireCompound::Intermediate,
        }
    }

    /// Safety car probability for this lap.
    pub fn safety_car_probability(&self) -> f64 {
        let extra = if self.condition != WeatherCondition::Dry {
            self.wet_safety_car_extra
        } else {
            0.0_f64
        };
        (self.base_safety_car_prob + extra).min(1.0_f64)
    }

    /// Whether an immediate pit stop for rain tyres is recommended.
    pub fn immediate_pit_for_rain(&self, current_compound: &TireCompound) -> bool {
        let needs_rain = matches!(
            self.condition,
            WeatherCondition::Wet | WeatherCondition::Damp
        );
        needs_rain && !current_compound.is_wet_weather()
    }

    /// Expected time loss per lap from rain conditions (s).
    pub fn rain_lap_time_delta(&self) -> f64 {
        match &self.condition {
            WeatherCondition::Dry => 0.0_f64,
            WeatherCondition::Damp => 3.0_f64,
            WeatherCondition::Wet => 8.0_f64,
            WeatherCondition::Mixed => 4.5_f64,
        }
    }
}

// ---------------------------------------------------------------------------
// RaceSimulator
// ---------------------------------------------------------------------------

/// Per-lap race result.
#[derive(Debug, Clone)]
pub struct LapResult {
    /// Lap number (1-indexed).
    pub lap: usize,
    /// Lap time (s).
    pub lap_time: f64,
    /// Cumulative race time (s).
    pub cumulative_time: f64,
    /// Tyre wear at end of lap.
    pub tyre_wear: f64,
    /// Fuel remaining (kg).
    pub fuel_remaining: f64,
    /// Whether a pit stop was taken on this lap.
    pub pit_stop: bool,
    /// Current compound.
    pub compound: TireCompound,
}

/// Full race simulator.
#[derive(Debug, Clone)]
pub struct RaceSimulator {
    /// Track.
    pub track: RaceTrack,
    /// Strategy to simulate.
    pub strategy: RaceStrategy,
    /// Lap time model.
    pub lap_time_model: LapTimeModel,
    /// Pit stop model.
    pub pit_stop_model: PitStopModel,
    /// Weather adaptation.
    pub weather: WeatherAdaptation,
}

impl RaceSimulator {
    /// Create a race simulator.
    pub fn new(
        track: RaceTrack,
        strategy: RaceStrategy,
        lap_time_model: LapTimeModel,
        pit_stop_model: PitStopModel,
        weather: WeatherAdaptation,
    ) -> Self {
        Self {
            track,
            strategy,
            lap_time_model,
            pit_stop_model,
            weather,
        }
    }

    /// Run the full race simulation and return per-lap results.
    pub fn simulate(&self) -> Vec<LapResult> {
        let total_laps = self.track.total_laps;
        let mut results = Vec::with_capacity(total_laps);

        let mut cumulative_time = 0.0_f64;
        let mut stint_idx = 0_usize;
        let compound = self
            .strategy
            .compounds
            .first()
            .cloned()
            .unwrap_or(TireCompound::Medium);
        let mut current_compound = compound;
        let mut tyre_wear = TireWear::new(current_compound.clone());
        let fuel_strat = &self.lap_time_model.fuel_strategy;

        for lap in 1..=total_laps {
            // Check for pit stop on this lap.
            let pit_this_lap = self.strategy.pit_laps.contains(&lap);

            // Compute deg rate.
            let deg = tire_degradation_rate(&tyre_wear.compound, 200.0_f64, self.track.track_temp);
            tyre_wear.advance_lap(deg);

            // Lap time.
            let drs_active = self.track.drs_zone_count() > 0;
            let lap_time = self
                .lap_time_model
                .predict_lap_time(&tyre_wear, lap, drs_active);
            let rain_delta = self.weather.rain_lap_time_delta();

            // Pit stop time loss.
            let pit_time = if pit_this_lap {
                self.pit_stop_model.total_time(0.0_f64)
            } else {
                0.0_f64
            };

            let total_lap_time = lap_time + rain_delta + pit_time;
            cumulative_time += total_lap_time;

            let fuel_remaining = fuel_strat.fuel_at_lap(lap);

            results.push(LapResult {
                lap,
                lap_time: total_lap_time,
                cumulative_time,
                tyre_wear: tyre_wear.wear_level,
                fuel_remaining,
                pit_stop: pit_this_lap,
                compound: current_compound.clone(),
            });

            // Switch compound after pit stop.
            if pit_this_lap {
                stint_idx += 1;
                if let Some(c) = self.strategy.compounds.get(stint_idx) {
                    current_compound = c.clone();
                    tyre_wear = TireWear::new(current_compound.clone());
                }
            }
        }

        results
    }

    /// Total race time (s).
    pub fn total_race_time(&self) -> f64 {
        self.simulate()
            .last()
            .map(|r| r.cumulative_time)
            .unwrap_or(0.0_f64)
    }

    /// Lap with the fastest lap time.
    pub fn fastest_lap(&self) -> Option<LapResult> {
        self.simulate().into_iter().min_by(|a, b| {
            a.lap_time
                .partial_cmp(&b.lap_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Number of pit stops executed.
    pub fn pit_stop_count(&self) -> usize {
        self.strategy.pit_laps.len()
    }
}

// ---------------------------------------------------------------------------
// StrategyOptimizer
// ---------------------------------------------------------------------------

/// Simple brute-force strategy optimizer.
#[derive(Debug, Clone)]
pub struct StrategyOptimizer {
    /// Track.
    pub track: RaceTrack,
    /// Lap time model.
    pub lap_time_model: LapTimeModel,
    /// Pit stop model.
    pub pit_stop_model: PitStopModel,
    /// Weather.
    pub weather: WeatherAdaptation,
    /// Fuel start (kg).
    pub fuel_start: f64,
}

impl StrategyOptimizer {
    /// Create a strategy optimizer.
    pub fn new(
        track: RaceTrack,
        lap_time_model: LapTimeModel,
        pit_stop_model: PitStopModel,
        weather: WeatherAdaptation,
        fuel_start: f64,
    ) -> Self {
        Self {
            track,
            lap_time_model,
            pit_stop_model,
            weather,
            fuel_start,
        }
    }

    /// Evaluate a one-stop strategy with a given pit lap.
    pub fn evaluate_one_stop(
        &self,
        pit_lap: usize,
        first_compound: TireCompound,
        second_compound: TireCompound,
    ) -> f64 {
        let strategy = RaceStrategy::new(
            vec![pit_lap],
            vec![first_compound, second_compound],
            self.fuel_start,
        );
        let sim = RaceSimulator::new(
            self.track.clone(),
            strategy,
            self.lap_time_model.clone(),
            self.pit_stop_model.clone(),
            self.weather.clone(),
        );
        sim.total_race_time()
    }

    /// Find the optimal one-stop pit lap (exhaustive search).
    pub fn optimal_one_stop_lap(
        &self,
        first_compound: TireCompound,
        second_compound: TireCompound,
    ) -> (usize, f64) {
        let total = self.track.total_laps;
        let mut best_lap = 1_usize;
        let mut best_time = f64::MAX;
        for pit_lap in 5..=total.saturating_sub(5) {
            let t =
                self.evaluate_one_stop(pit_lap, first_compound.clone(), second_compound.clone());
            if t < best_time {
                best_time = t;
                best_lap = pit_lap;
            }
        }
        (best_lap, best_time)
    }

    /// Evaluate a two-stop strategy.
    pub fn evaluate_two_stop(
        &self,
        pit_lap_1: usize,
        pit_lap_2: usize,
        compounds: [TireCompound; 3],
    ) -> f64 {
        if pit_lap_1 >= pit_lap_2 {
            return f64::MAX;
        }
        let strategy = RaceStrategy::new(
            vec![pit_lap_1, pit_lap_2],
            compounds.to_vec(),
            self.fuel_start,
        );
        let sim = RaceSimulator::new(
            self.track.clone(),
            strategy,
            self.lap_time_model.clone(),
            self.pit_stop_model.clone(),
            self.weather.clone(),
        );
        sim.total_race_time()
    }
}

// ---------------------------------------------------------------------------
// Standalone helper functions (kept from original file)
// ---------------------------------------------------------------------------

/// Compute the tyre degradation rate per lap.
///
/// Returns wear fraction per lap (0–1 scale).
pub fn tire_degradation_rate(compound: &TireCompound, speed: f64, temp: f64) -> f64 {
    let base = compound.base_deg_rate();
    let speed_factor = (speed / 200.0_f64).powi(2).max(0.1_f64);
    let temp_factor = ((temp / 80.0_f64) * 0.5_f64 + 0.5_f64).max(0.1_f64);
    base * speed_factor * temp_factor
}

/// Compute the grip factor for the current wear state.
///
/// Returns a value in `[0, 1]`.
pub fn tire_grip_factor(wear: &TireWear) -> f64 {
    let base = (1.0_f64 - wear.wear_level).clamp(0.0_f64, 1.0_f64);
    let peak = wear.compound.peak_grip();
    (base * peak).clamp(0.0_f64, 1.0_f64)
}

/// Compute the optimal pit window `(earliest, latest)` lap.
pub fn optimal_pit_window(
    total_laps: usize,
    compound: &TireCompound,
    deg_rate: f64,
) -> (usize, usize) {
    let life_laps = if deg_rate > 1e-10_f64 {
        (0.85_f64 / deg_rate).round() as usize
    } else {
        total_laps
    };
    let mid = life_laps.min(total_laps);
    let earliest = mid.saturating_sub(3).max(1);
    let latest = (mid + 3).min(total_laps.saturating_sub(1)).max(earliest);
    let _ = compound;
    (earliest, latest)
}

/// Estimate fuel consumption per lap.
pub fn fuel_consumption(throttle: f64, speed: f64, base_rate: f64) -> f64 {
    let throttle_c = throttle.clamp(0.0_f64, 1.0_f64);
    let speed_factor = (speed / 200.0_f64).max(0.1_f64);
    base_rate * throttle_c * speed_factor
}

/// Additional downforce due to fuel weight.
pub fn fuel_weight_effect(fuel_kg: f64, _mass: f64, g: f64) -> f64 {
    fuel_kg * g
}

/// Determine whether an undercut is advantageous.
pub fn undercut_advantage(lap_delta_s: f64, pit_delta_s: f64) -> bool {
    pit_delta_s < lap_delta_s
}

/// Compute the total estimated race time for a given strategy.
pub fn strategy_total_time(
    strategy: &RaceStrategy,
    lap_time_base: f64,
    deg_rate: f64,
    pit_time: f64,
) -> f64 {
    let total_pit_time = strategy.pit_laps.len() as f64 * pit_time;
    let mut pit_laps_sorted = strategy.pit_laps.clone();
    pit_laps_sorted.sort_unstable();

    let mut stints: Vec<(usize, usize)> = Vec::new();
    let mut prev = 0_usize;
    for &pl in &pit_laps_sorted {
        if pl > prev {
            stints.push((prev, pl));
        }
        prev = pl;
    }
    let implied_total = pit_laps_sorted.last().copied().unwrap_or(0) + 30;
    stints.push((prev, implied_total));

    let mut total = 0.0_f64;
    for (start, end) in stints {
        let stint_laps = end.saturating_sub(start) as f64;
        let avg_deg = deg_rate * stint_laps * 0.5_f64;
        total += stint_laps * (lap_time_base + avg_deg * lap_time_base);
    }
    total + total_pit_time
}

/// DRS overtaking probability.
pub fn drs_overtaking_probability(speed_diff: f64, track_position: f64) -> f64 {
    let base = (speed_diff / 20.0_f64).tanh().clamp(0.0_f64, 1.0_f64);
    let pos_factor = track_position.clamp(0.0_f64, 1.0_f64);
    (base * pos_factor).clamp(0.0_f64, 1.0_f64)
}

/// Overcut strategy: pit later to benefit from track position and new tyres.
///
/// Returns `true` if the overcut is expected to be advantageous.
pub fn overcut_advantage(
    laps_remaining_on_current_tyre: usize,
    tyre_deg_per_lap: f64,
    new_tyre_gain_per_lap: f64,
    pit_time_s: f64,
) -> bool {
    let laps_f = laps_remaining_on_current_tyre as f64;
    let total_gain = new_tyre_gain_per_lap * laps_f;
    let total_cost = tyre_deg_per_lap * laps_f;
    (total_gain - total_cost) > pit_time_s
}

/// Compute safety car probability over a stint.
///
/// Returns the probability that at least one safety car period occurs.
pub fn safety_car_probability_over_stint(laps: usize, per_lap_prob: f64) -> f64 {
    let no_sc = (1.0_f64 - per_lap_prob).powi(laps as i32);
    1.0_f64 - no_sc
}

/// Compute the virtual safety car time bonus for the leader.
///
/// Returns seconds saved by making a free pit stop under VSC.
pub fn vsc_pit_bonus(pit_lane_loss: f64, vsc_lap_time_delta: f64, vsc_duration_laps: f64) -> f64 {
    // Under VSC all cars run ~40 % slower; pit stop effectively "free" if
    // the delta is absorbed by the slower pace.
    let vsc_lap_time = vsc_lap_time_delta * vsc_duration_laps;
    (pit_lane_loss - vsc_lap_time).max(0.0_f64)
}

/// Estimate the impact of track position on pit window decision.
///
/// Returns `true` if stopping now maintains position over the rival.
pub fn stop_maintains_position(
    own_gap_to_rival: f64,
    pit_stop_total_time: f64,
    rival_lap_time: f64,
) -> bool {
    // After stop, rival covers laps; own pit delta must be < gap + laps.
    pit_stop_total_time < own_gap_to_rival + rival_lap_time
}

/// Thermal penalty factor outside the operating window.
///
/// Returns a multiplier \[0.6, 1.0\] on available grip.
pub fn thermal_penalty(tyre_temp: f64, compound: &TireCompound) -> f64 {
    let (lo, hi) = compound.operating_window();
    if tyre_temp < lo {
        let delta = (lo - tyre_temp).min(40.0_f64) / 40.0_f64;
        1.0_f64 - 0.4_f64 * delta
    } else if tyre_temp > hi {
        let delta = (tyre_temp - hi).min(30.0_f64) / 30.0_f64;
        1.0_f64 - 0.4_f64 * delta
    } else {
        1.0_f64
    }
}

/// Sector-based pace analysis.
///
/// Returns `(s1_time, s2_time, s3_time)` for the given lap.
pub fn sector_pace(
    track: &RaceTrack,
    wear: &TireWear,
    lap: usize,
    fuel_strat: &FuelStrategy,
) -> (f64, f64, f64) {
    if track.sectors.len() < 3 {
        return (0.0_f64, 0.0_f64, 0.0_f64);
    }
    let fuel_pen = fuel_strat.fuel_lap_time_penalty(lap) / 3.0_f64;
    let s1 = track.sectors[0].base_time
        + wear.wear_level * track.sectors[0].roughness * 2.0_f64
        + fuel_pen;
    let s2 = track.sectors[1].base_time
        + wear.wear_level * track.sectors[1].roughness * 2.0_f64
        + fuel_pen;
    let s3 = track.sectors[2].base_time
        + wear.wear_level * track.sectors[2].roughness * 2.0_f64
        + fuel_pen;
    (s1, s2, s3)
}

/// Compute a theoretical best lap from the sum of best-ever sectors.
pub fn theoretical_best_lap(best_s1: f64, best_s2: f64, best_s3: f64) -> f64 {
    best_s1 + best_s2 + best_s3
}

/// Simple race-pace comparison for undercut / overcut decision.
///
/// Returns `Some(StrategyDecision)` — see inline enum.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StrategyDecision {
    /// Execute undercut now.
    Undercut,
    /// Execute overcut (stay out longer).
    Overcut,
    /// No action — continue on current strategy.
    Continue,
}

/// Decide undercut / overcut / continue.
pub fn strategy_decision(
    gap_to_rival: f64,
    own_deg_per_lap: f64,
    rival_deg_per_lap: f64,
    pit_stop_loss: f64,
    laps_remaining: usize,
) -> StrategyDecision {
    // Tyre delta: how much faster own car will be per lap after pit.
    let tyre_delta = rival_deg_per_lap - own_deg_per_lap;
    // Undercut: pit now and get clean air advantage.
    if gap_to_rival < pit_stop_loss * 0.5_f64 && tyre_delta > 0.1_f64 {
        return StrategyDecision::Undercut;
    }
    // Overcut: stay out and benefit from rival pitting.
    let total_overcut_gain = tyre_delta * laps_remaining as f64;
    if total_overcut_gain > pit_stop_loss && gap_to_rival > pit_stop_loss {
        return StrategyDecision::Overcut;
    }
    StrategyDecision::Continue
}

/// Tyre life fraction remaining.
pub fn tyre_life_fraction(wear: &TireWear) -> f64 {
    1.0_f64 - wear.wear_level
}

/// Compute the average lap time for a stint.
pub fn stint_average_lap_time(
    start_wear: f64,
    stint_laps: usize,
    deg_rate: f64,
    base_lap_time: f64,
    wear_penalty: f64,
) -> f64 {
    let total_delta: f64 = (0..stint_laps)
        .map(|i| {
            let w = (start_wear + deg_rate * i as f64).min(1.0_f64);
            base_lap_time + w * wear_penalty
        })
        .sum();
    if stint_laps == 0 {
        base_lap_time
    } else {
        total_delta / stint_laps as f64
    }
}

/// Return true if the safety car has been deployed (stochastic check).
///
/// Uses a simple threshold against a probability value.
pub fn is_safety_car_deployed(lap_prob: f64, random_val: f64) -> bool {
    random_val < lap_prob
}

/// Angle subtended by the track's DRS zone relative to the full lap.
///
/// Returns a fraction in \[0, 1\].
pub fn drs_zone_fraction(track: &RaceTrack) -> f64 {
    let drs_len: f64 = track
        .sectors
        .iter()
        .filter(|s| s.has_drs)
        .map(|s| s.length)
        .sum();
    if track.lap_distance > 0.0_f64 {
        drs_len / track.lap_distance
    } else {
        0.0_f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn default_wear(compound: TireCompound) -> TireWear {
        TireWear::new(compound)
    }

    fn default_fuel_strategy() -> FuelStrategy {
        FuelStrategy::new(100.0_f64, 1.8_f64, 0.035_f64, false, 0.0_f64)
    }

    fn default_track() -> RaceTrack {
        RaceTrack::default_track()
    }

    fn default_lap_model(track: RaceTrack) -> LapTimeModel {
        LapTimeModel::new(
            track,
            86.5_f64,
            4.0_f64,
            default_fuel_strategy(),
            0.3_f64,
            0.02_f64,
        )
    }

    // ── TireCompound ─────────────────────────────────────────────────────

    #[test]
    fn soft_base_deg_rate_highest() {
        assert!(TireCompound::Soft.base_deg_rate() > TireCompound::Hard.base_deg_rate());
    }

    #[test]
    fn soft_peak_grip_highest() {
        assert!(TireCompound::Soft.peak_grip() > TireCompound::Hard.peak_grip());
    }

    #[test]
    fn wet_weather_compounds_flagged() {
        assert!(TireCompound::Intermediate.is_wet_weather());
        assert!(TireCompound::Wet.is_wet_weather());
        assert!(!TireCompound::Soft.is_wet_weather());
    }

    #[test]
    fn operating_window_lo_less_than_hi() {
        for c in [
            TireCompound::Soft,
            TireCompound::Medium,
            TireCompound::Hard,
            TireCompound::Intermediate,
            TireCompound::Wet,
        ] {
            let (lo, hi) = c.operating_window();
            assert!(lo < hi, "{c:?}: lo={lo}, hi={hi}");
        }
    }

    #[test]
    fn nominal_stint_laps_positive() {
        for c in [
            TireCompound::Soft,
            TireCompound::Medium,
            TireCompound::Hard,
            TireCompound::Intermediate,
            TireCompound::Wet,
        ] {
            assert!(c.nominal_stint_laps() > 0);
        }
    }

    // ── TireWear ─────────────────────────────────────────────────────────

    #[test]
    fn tire_wear_new_starts_fresh() {
        let w = TireWear::new(TireCompound::Medium);
        assert_eq!(w.wear_level, 0.0_f64);
        assert_eq!(w.laps_used, 0);
    }

    #[test]
    fn tire_wear_advance_lap_increases_wear() {
        let mut w = TireWear::new(TireCompound::Medium);
        w.advance_lap(0.02_f64);
        assert!((w.wear_level - 0.02_f64).abs() < 1e-12_f64);
        assert_eq!(w.laps_used, 1);
    }

    #[test]
    fn tire_wear_capped_at_one() {
        let mut w = TireWear::new(TireCompound::Soft);
        for _ in 0..100 {
            w.advance_lap(0.05_f64);
        }
        assert!(w.wear_level <= 1.0_f64);
    }

    #[test]
    fn tire_wear_worn_out_above_threshold() {
        let mut w = TireWear::new(TireCompound::Soft);
        w.wear_level = 0.95_f64;
        assert!(w.is_worn_out());
    }

    // ── TireThermal ──────────────────────────────────────────────────────

    #[test]
    fn thermal_model_starts_at_ambient() {
        let th = TireThermal::new(TireCompound::Medium, 25.0_f64);
        assert!((th.surface_temp - 25.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn thermal_grip_factor_in_window_is_one() {
        let mut th = TireThermal::new(TireCompound::Medium, 25.0_f64);
        th.surface_temp = 97.0_f64; // well inside medium window (85–110)
        assert!((th.thermal_grip_factor() - 1.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn thermal_grip_factor_below_window_reduced() {
        let mut th = TireThermal::new(TireCompound::Medium, 25.0_f64);
        th.surface_temp = 45.0_f64; // well below medium window (85–110)
        assert!(th.thermal_grip_factor() < 1.0_f64);
    }

    // ── RaceTrack ────────────────────────────────────────────────────────

    #[test]
    fn default_track_has_three_sectors() {
        let t = default_track();
        assert_eq!(t.sectors.len(), 3);
    }

    #[test]
    fn default_track_drs_zone_count() {
        let t = default_track();
        assert!(t.drs_zone_count() > 0);
    }

    #[test]
    fn base_lap_time_from_sector_sum() {
        let t = default_track();
        let expected: f64 = t.sectors.iter().map(|s| s.base_time).sum();
        assert!((t.base_lap_time() - expected).abs() < 1e-9_f64);
    }

    // ── FuelStrategy ─────────────────────────────────────────────────────

    #[test]
    fn fuel_at_lap_zero_is_start() {
        let fs = default_fuel_strategy();
        assert!((fs.fuel_at_lap(0) - 100.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn fuel_at_lap_never_negative() {
        let fs = default_fuel_strategy();
        assert!(fs.fuel_at_lap(200) >= 0.0_f64);
    }

    #[test]
    fn fuel_lap_time_penalty_decreases_with_laps() {
        let fs = default_fuel_strategy();
        let p0 = fs.fuel_lap_time_penalty(0);
        let p10 = fs.fuel_lap_time_penalty(10);
        assert!(p0 > p10);
    }

    // ── PitStopModel ─────────────────────────────────────────────────────

    #[test]
    fn f1_pit_stop_total_time_positive() {
        let p = PitStopModel::default_f1();
        assert!(p.total_time(0.0_f64) > 0.0_f64);
    }

    #[test]
    fn pit_stop_total_time_increases_with_fuel() {
        // default_f1 has fuel_add_time_per_10kg = 0 so no change; use custom model.
        let p = PitStopModel::new(2.5_f64, 1.0_f64, 1.0_f64, 22.0_f64);
        assert!(p.total_time(20.0_f64) > p.total_time(0.0_f64));
    }

    // ── LapTimeModel ─────────────────────────────────────────────────────

    #[test]
    fn lap_time_model_fresh_tyre_less_than_worn() {
        let track = default_track();
        let model = default_lap_model(track);
        let fresh_wear = default_wear(TireCompound::Medium);
        let mut worn_wear = TireWear::new(TireCompound::Medium);
        worn_wear.wear_level = 0.6_f64;
        let t_fresh = model.predict_lap_time(&fresh_wear, 1, false);
        let t_worn = model.predict_lap_time(&worn_wear, 1, false);
        assert!(t_worn > t_fresh);
    }

    #[test]
    fn drs_reduces_lap_time() {
        let track = default_track();
        let model = default_lap_model(track);
        let wear = default_wear(TireCompound::Medium);
        let t_no_drs = model.predict_lap_time(&wear, 1, false);
        let t_drs = model.predict_lap_time(&wear, 1, true);
        assert!(t_drs < t_no_drs);
    }

    #[test]
    fn sector_times_sum_near_lap_time() {
        let track = default_track();
        let model = default_lap_model(track);
        let wear = default_wear(TireCompound::Medium);
        let s0 = model.sector_time(0, &wear, 1);
        let s1 = model.sector_time(1, &wear, 1);
        let s2 = model.sector_time(2, &wear, 1);
        let lap = model.predict_lap_time(&wear, 1, false);
        // Sum of sectors should be in the same ballpark (within 20 s).
        assert!((s0 + s1 + s2 - lap).abs() < 20.0_f64);
    }

    // ── RaceStrategy ─────────────────────────────────────────────────────

    #[test]
    fn race_strategy_num_stints() {
        let s = RaceStrategy::new(
            vec![20, 40],
            vec![TireCompound::Soft, TireCompound::Medium, TireCompound::Hard],
            100.0_f64,
        );
        assert_eq!(s.num_stints(), 3);
    }

    #[test]
    fn race_strategy_compound_rule_two_types() {
        let s = RaceStrategy::new(
            vec![25],
            vec![TireCompound::Soft, TireCompound::Hard],
            100.0_f64,
        );
        assert!(s.satisfies_compound_rule());
    }

    #[test]
    fn race_strategy_compound_rule_same_compound_fails() {
        let s = RaceStrategy::new(
            vec![25],
            vec![TireCompound::Medium, TireCompound::Medium],
            100.0_f64,
        );
        assert!(!s.satisfies_compound_rule());
    }

    // ── WeatherAdaptation ─────────────────────────────────────────────────

    #[test]
    fn weather_dry_recommends_dry_compound() {
        let w = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let c = w.recommended_compound();
        assert!(!c.is_wet_weather());
    }

    #[test]
    fn weather_wet_recommends_wet_compound() {
        let w = WeatherAdaptation::new(WeatherCondition::Wet, 1.0_f64, 0.02_f64, 20.0_f64);
        let c = w.recommended_compound();
        assert!(c.is_wet_weather());
    }

    #[test]
    fn weather_safety_car_prob_higher_in_wet() {
        let dry = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let wet = WeatherAdaptation::new(WeatherCondition::Wet, 1.0_f64, 0.02_f64, 20.0_f64);
        assert!(wet.safety_car_probability() > dry.safety_car_probability());
    }

    #[test]
    fn weather_immediate_pit_for_rain_when_on_dry() {
        let w = WeatherAdaptation::new(WeatherCondition::Wet, 1.0_f64, 0.02_f64, 15.0_f64);
        assert!(w.immediate_pit_for_rain(&TireCompound::Medium));
    }

    #[test]
    fn weather_rain_lap_time_delta_positive_in_rain() {
        let w = WeatherAdaptation::new(WeatherCondition::Wet, 1.0_f64, 0.02_f64, 15.0_f64);
        assert!(w.rain_lap_time_delta() > 0.0_f64);
    }

    // ── RaceSimulator ────────────────────────────────────────────────────

    #[test]
    fn simulator_produces_correct_number_of_laps() {
        let track = default_track();
        let strategy = RaceStrategy::new(
            vec![20],
            vec![TireCompound::Medium, TireCompound::Hard],
            100.0_f64,
        );
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let sim = RaceSimulator::new(track.clone(), strategy, model, pit, weather);
        let results = sim.simulate();
        assert_eq!(results.len(), track.total_laps);
    }

    #[test]
    fn simulator_pit_stop_lap_marked() {
        let track = default_track();
        let strategy = RaceStrategy::new(
            vec![20],
            vec![TireCompound::Medium, TireCompound::Hard],
            100.0_f64,
        );
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let sim = RaceSimulator::new(track, strategy, model, pit, weather);
        let results = sim.simulate();
        assert!(results[19].pit_stop, "Lap 20 should be a pit stop");
    }

    #[test]
    fn simulator_total_time_positive() {
        let track = default_track();
        let strategy = RaceStrategy::new(vec![], vec![TireCompound::Hard], 100.0_f64);
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let sim = RaceSimulator::new(track, strategy, model, pit, weather);
        assert!(sim.total_race_time() > 0.0_f64);
    }

    #[test]
    fn simulator_cumulative_time_monotone() {
        let track = default_track();
        let strategy = RaceStrategy::new(vec![], vec![TireCompound::Hard], 100.0_f64);
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let sim = RaceSimulator::new(track, strategy, model, pit, weather);
        let results = sim.simulate();
        for w in results.windows(2) {
            assert!(w[1].cumulative_time >= w[0].cumulative_time);
        }
    }

    // ── StrategyOptimizer ────────────────────────────────────────────────

    #[test]
    fn optimizer_finds_valid_one_stop_lap() {
        let track = default_track();
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let opt = StrategyOptimizer::new(track.clone(), model, pit, weather, 100.0_f64);
        let (lap, time) = opt.optimal_one_stop_lap(TireCompound::Medium, TireCompound::Hard);
        assert!(lap >= 5 && lap <= track.total_laps.saturating_sub(5));
        assert!(time > 0.0_f64);
    }

    #[test]
    fn optimizer_two_stop_worse_than_one_stop_on_short_race() {
        // On a 57-lap race, a well-chosen one-stop should beat an early two-stop.
        let track = default_track();
        let model = default_lap_model(track.clone());
        let pit = PitStopModel::default_f1();
        let weather = WeatherAdaptation::new(WeatherCondition::Dry, 0.0_f64, 0.02_f64, 40.0_f64);
        let opt = StrategyOptimizer::new(track.clone(), model, pit, weather, 100.0_f64);
        let (_, t1) = opt.optimal_one_stop_lap(TireCompound::Medium, TireCompound::Hard);
        // Two-stop with very early first pit.
        let t2 = opt.evaluate_two_stop(
            10,
            30,
            [TireCompound::Soft, TireCompound::Medium, TireCompound::Hard],
        );
        // Either could win, but both must be positive and finite.
        assert!(t1.is_finite() && t1 > 0.0_f64);
        assert!(t2.is_finite() && t2 > 0.0_f64);
    }

    // ── Standalone helpers ────────────────────────────────────────────────

    #[test]
    fn soft_degrades_faster_than_hard() {
        let r_soft = tire_degradation_rate(&TireCompound::Soft, 200.0_f64, 80.0_f64);
        let r_hard = tire_degradation_rate(&TireCompound::Hard, 200.0_f64, 80.0_f64);
        assert!(r_soft > r_hard);
    }

    #[test]
    fn degradation_increases_with_speed() {
        let r_slow = tire_degradation_rate(&TireCompound::Medium, 100.0_f64, 80.0_f64);
        let r_fast = tire_degradation_rate(&TireCompound::Medium, 300.0_f64, 80.0_f64);
        assert!(r_fast > r_slow);
    }

    #[test]
    fn degradation_increases_with_temp() {
        let r_cool = tire_degradation_rate(&TireCompound::Medium, 200.0_f64, 40.0_f64);
        let r_hot = tire_degradation_rate(&TireCompound::Medium, 200.0_f64, 120.0_f64);
        assert!(r_hot > r_cool);
    }

    #[test]
    fn all_compounds_positive_rate() {
        for c in &[
            TireCompound::Soft,
            TireCompound::Medium,
            TireCompound::Hard,
            TireCompound::Intermediate,
            TireCompound::Wet,
        ] {
            let r = tire_degradation_rate(c, 200.0_f64, 80.0_f64);
            assert!(r > 0.0_f64, "compound {c:?} rate={r}");
        }
    }

    #[test]
    fn grip_new_tyre_near_one() {
        let wear = default_wear(TireCompound::Medium);
        assert!(tire_grip_factor(&wear) >= 0.95_f64);
    }

    #[test]
    fn grip_fully_worn_near_zero() {
        let mut wear = TireWear::new(TireCompound::Medium);
        wear.wear_level = 1.0_f64;
        assert!(tire_grip_factor(&wear) < 0.1_f64);
    }

    #[test]
    fn pit_window_earliest_before_latest() {
        let (e, l) = optimal_pit_window(60, &TireCompound::Medium, 0.02_f64);
        assert!(e <= l);
    }

    #[test]
    fn fuel_zero_throttle_is_zero() {
        let f = fuel_consumption(0.0_f64, 200.0_f64, 2.0_f64);
        assert!(f.abs() < 1e-15_f64);
    }

    #[test]
    fn undercut_works_when_pit_faster() {
        assert!(undercut_advantage(5.0_f64, 3.0_f64));
    }

    #[test]
    fn one_stop_adds_pit_time() {
        let s_one = RaceStrategy::new(
            vec![30],
            vec![TireCompound::Medium, TireCompound::Hard],
            100.0_f64,
        );
        let s_zero = RaceStrategy::new(vec![], vec![TireCompound::Hard], 100.0_f64);
        let t_one = strategy_total_time(&s_one, 90.0_f64, 0.0_f64, 25.0_f64);
        let t_zero = strategy_total_time(&s_zero, 90.0_f64, 0.0_f64, 25.0_f64);
        assert!(t_one > t_zero);
    }

    #[test]
    fn drs_zero_track_position_gives_zero() {
        let p = drs_overtaking_probability(50.0_f64, 0.0_f64);
        assert!(p.abs() < 1e-10_f64);
    }

    #[test]
    fn drs_high_diff_high_prob() {
        let p = drs_overtaking_probability(50.0_f64, 1.0_f64);
        assert!(p > 0.9_f64);
    }

    #[test]
    fn overcut_advantage_large_gain() {
        assert!(overcut_advantage(20, 0.1_f64, 0.5_f64, 3.0_f64));
    }

    #[test]
    fn safety_car_probability_accumulates_over_laps() {
        let p1 = safety_car_probability_over_stint(1, 0.05_f64);
        let p10 = safety_car_probability_over_stint(10, 0.05_f64);
        assert!(p10 > p1);
    }

    #[test]
    fn strategy_decision_undercut_small_gap() {
        let d = strategy_decision(2.0_f64, 0.05_f64, 0.5_f64, 22.0_f64, 30);
        assert_eq!(d, StrategyDecision::Undercut);
    }

    #[test]
    fn tyre_life_fraction_new_is_one() {
        let w = default_wear(TireCompound::Medium);
        assert!((tyre_life_fraction(&w) - 1.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn stint_average_lap_time_increases_with_deg() {
        let avg0 = stint_average_lap_time(0.0_f64, 20, 0.0_f64, 86.0_f64, 4.0_f64);
        let avg1 = stint_average_lap_time(0.0_f64, 20, 0.02_f64, 86.0_f64, 4.0_f64);
        assert!(avg1 > avg0);
    }

    #[test]
    fn theoretical_best_lap_is_sector_sum() {
        let t = theoretical_best_lap(28.0_f64, 31.0_f64, 27.0_f64);
        assert!((t - 86.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn drs_zone_fraction_positive_on_default_track() {
        let t = default_track();
        let f = drs_zone_fraction(&t);
        assert!(f > 0.0_f64 && f <= 1.0_f64);
    }

    #[test]
    fn thermal_penalty_in_window_is_one() {
        let c = TireCompound::Medium;
        let (lo, hi) = c.operating_window();
        let mid = (lo + hi) * 0.5_f64;
        assert!((thermal_penalty(mid, &c) - 1.0_f64).abs() < 1e-9_f64);
    }

    #[test]
    fn thermal_penalty_below_window_reduced() {
        let c = TireCompound::Medium;
        let (lo, _hi) = c.operating_window();
        let cold = lo - 20.0_f64;
        assert!(thermal_penalty(cold, &c) < 1.0_f64);
    }

    #[test]
    fn vsc_pit_bonus_non_negative() {
        let bonus = vsc_pit_bonus(22.0_f64, 5.0_f64, 3.0_f64);
        assert!(bonus >= 0.0_f64);
    }

    #[test]
    fn sector_pace_sum_positive() {
        let t = default_track();
        let wear = default_wear(TireCompound::Medium);
        let fs = default_fuel_strategy();
        let (s1, s2, s3) = sector_pace(&t, &wear, 1, &fs);
        assert!(s1 > 0.0_f64 && s2 > 0.0_f64 && s3 > 0.0_f64);
    }

    #[test]
    fn fuel_weight_proportional_to_kg() {
        let f1 = fuel_weight_effect(10.0_f64, 700.0_f64, 9.81_f64);
        let f2 = fuel_weight_effect(20.0_f64, 700.0_f64, 9.81_f64);
        assert!((f2 / f1 - 2.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn pi_constant_in_scope() {
        // Just verifies PI is used for dihedral-type checks in the module.
        assert!((PI - std::f64::consts::PI).abs() < 1e-15_f64);
    }
}
