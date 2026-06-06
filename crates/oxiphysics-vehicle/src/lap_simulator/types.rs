//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
/// Complete race track layout.
#[derive(Debug, Clone)]
pub struct TrackLayout {
    /// Ordered list of track segments.
    pub segments: Vec<TrackSegment>,
    /// Total track length (m) — computed from segments.
    pub total_length: f64,
    /// Sector boundary indices into `segments`.
    pub sector_boundaries: Vec<usize>,
    /// Indices of DRS activation zones.
    pub drs_zones: Vec<usize>,
}
impl TrackLayout {
    /// Construct a track layout from segments.
    pub fn new(segments: Vec<TrackSegment>) -> Self {
        let total_length: f64 = segments.iter().map(|s| s.length()).sum();
        let drs_zones: Vec<usize> = segments
            .iter()
            .enumerate()
            .filter(|(_, s)| s.is_drs_zone())
            .map(|(i, _)| i)
            .collect();
        Self {
            sector_boundaries: vec![segments.len() / 3, 2 * segments.len() / 3],
            segments,
            total_length,
            drs_zones,
        }
    }
    /// Create a simplified 3-sector layout representative of a modern F1 circuit.
    pub fn f1_example() -> Self {
        let segments = vec![
            TrackSegment::Straight {
                length: 600.0,
                drs_available: true,
            },
            TrackSegment::Corner {
                radius: 50.0,
                angle: std::f64::consts::FRAC_PI_2,
                direction: 1.0,
            },
            TrackSegment::Straight {
                length: 150.0,
                drs_available: false,
            },
            TrackSegment::Chicane {
                radii: vec![30.0, 30.0],
                angles: vec![std::f64::consts::FRAC_PI_4; 2],
            },
            TrackSegment::Straight {
                length: 400.0,
                drs_available: true,
            },
            TrackSegment::Corner {
                radius: 120.0,
                angle: std::f64::consts::PI,
                direction: -1.0,
            },
            TrackSegment::Elevation {
                length: 200.0,
                gradient: 0.05,
            },
            TrackSegment::Corner {
                radius: 80.0,
                angle: std::f64::consts::FRAC_PI_2,
                direction: 1.0,
            },
            TrackSegment::Straight {
                length: 300.0,
                drs_available: true,
            },
        ];
        Self::new(segments)
    }
    /// Total number of corners on the track.
    pub fn corner_count(&self) -> usize {
        self.segments
            .iter()
            .filter(|s| {
                matches!(
                    s,
                    TrackSegment::Corner { .. } | TrackSegment::Chicane { .. }
                )
            })
            .count()
    }
    /// Number of DRS zones.
    pub fn drs_zone_count(&self) -> usize {
        self.drs_zones.len()
    }
}
/// Vehicle performance map: acceleration limits vs speed.
#[derive(Debug, Clone)]
pub struct VehiclePerformanceMap {
    /// Speed breakpoints (m/s).
    pub speed_points: Vec<f64>,
    /// Maximum lateral acceleration at each speed (m/s²).
    pub max_lat_g: Vec<f64>,
    /// Maximum longitudinal braking deceleration (m/s²) — positive value.
    pub max_brake_g: Vec<f64>,
    /// Maximum longitudinal traction acceleration (m/s²).
    pub max_traction_g: Vec<f64>,
    /// Aerodynamic drag force (N) at each speed.
    pub drag_n: Vec<f64>,
    /// Aerodynamic lift (negative = downforce) at each speed (N).
    pub lift_n: Vec<f64>,
    /// Vehicle mass (kg).
    pub mass: f64,
    /// Maximum speed (m/s) — terminal velocity.
    pub v_max: f64,
}
impl VehiclePerformanceMap {
    /// Default F1 car performance map (approximate 2024 car).
    pub fn default_f1() -> Self {
        let speeds: Vec<f64> = (0..=10).map(|i| i as f64 * 10.0).collect();
        let lat_g: Vec<f64> = speeds
            .iter()
            .map(|&v| {
                let down = 0.5 * 1.2 * 3.5 * v * v;
                (4.0 * 9.81 + down / 800.0 * 9.81).min(55.0)
            })
            .collect();
        let brake_g: Vec<f64> = speeds
            .iter()
            .map(|&v| (5.5 * 9.81 + 0.5 * 1.2 * 0.8 * v * v / 800.0).min(65.0))
            .collect();
        let traction_g: Vec<f64> = speeds
            .iter()
            .map(|&v| {
                if v < 20.0 {
                    1.8 * 9.81
                } else {
                    (1.4 * 9.81 - (v - 20.0) * 0.02).max(0.5 * 9.81)
                }
            })
            .collect();
        let drag: Vec<f64> = speeds.iter().map(|&v| 0.5 * 1.2 * 1.0 * v * v).collect();
        let lift: Vec<f64> = speeds.iter().map(|&v| -0.5 * 1.2 * 3.5 * v * v).collect();
        Self {
            speed_points: speeds,
            max_lat_g: lat_g,
            max_brake_g: brake_g,
            max_traction_g: traction_g,
            drag_n: drag,
            lift_n: lift,
            mass: 800.0,
            v_max: 95.0,
        }
    }
    /// Interpolate a value from a map given speed.
    fn interpolate(speed: f64, speeds: &[f64], values: &[f64]) -> f64 {
        if speeds.is_empty() || values.is_empty() {
            return 0.0;
        }
        let s = speed
            .max(speeds[0])
            .min(*speeds.last().expect("collection should not be empty"));
        for i in 1..speeds.len() {
            if s <= speeds[i] {
                let t = (s - speeds[i - 1]) / (speeds[i] - speeds[i - 1]).max(1e-12);
                return values[i - 1] + t * (values[i] - values[i - 1]);
            }
        }
        *values.last().expect("collection should not be empty")
    }
    /// Maximum lateral acceleration (m/s²) at the given speed.
    pub fn lat_g_at(&self, speed: f64) -> f64 {
        Self::interpolate(speed, &self.speed_points, &self.max_lat_g)
    }
    /// Maximum braking deceleration (m/s²) at the given speed.
    pub fn brake_g_at(&self, speed: f64) -> f64 {
        Self::interpolate(speed, &self.speed_points, &self.max_brake_g)
    }
    /// Maximum traction acceleration (m/s²) at the given speed.
    pub fn traction_g_at(&self, speed: f64) -> f64 {
        Self::interpolate(speed, &self.speed_points, &self.max_traction_g)
    }
    /// Drag force (N) at the given speed.
    pub fn drag_at(&self, speed: f64) -> f64 {
        Self::interpolate(speed, &self.speed_points, &self.drag_n)
    }
}
/// Sector time record for a single lap.
#[derive(Debug, Clone)]
pub struct SectorTiming {
    /// Sector index (0-based).
    pub sector: usize,
    /// Sector time this lap (s).
    pub time: f64,
    /// Personal best sector time (s).
    pub personal_best: f64,
    /// Theoretical best sector time (absolute minimum from all laps, s).
    pub theoretical_best: f64,
}
impl SectorTiming {
    /// Create a new sector timing record.
    pub fn new(sector: usize, time: f64) -> Self {
        Self {
            sector,
            time,
            personal_best: time,
            theoretical_best: time,
        }
    }
    /// Update with a new sector time. Updates personal best if improved.
    pub fn update(&mut self, new_time: f64) {
        self.time = new_time;
        if new_time < self.personal_best {
            self.personal_best = new_time;
        }
        if new_time < self.theoretical_best {
            self.theoretical_best = new_time;
        }
    }
    /// Delta from personal best (positive = slower than PB).
    pub fn delta_to_pb(&self) -> f64 {
        self.time - self.personal_best
    }
    /// Delta from theoretical best (positive = slower).
    pub fn delta_to_theoretical(&self) -> f64 {
        self.time - self.theoretical_best
    }
    /// Returns `true` if current time equals the personal best.
    pub fn is_personal_best(&self) -> bool {
        (self.time - self.personal_best).abs() < 1e-9
    }
}
/// Time delta between two laps at each sampled position.
#[derive(Debug, Clone)]
pub struct LapDelta {
    /// Track positions of sample points (m).
    pub positions: Vec<f64>,
    /// Cumulative time delta at each position (comparison − reference, s).
    /// Negative = comparison car is ahead.
    pub delta: Vec<f64>,
    /// Final lap time delta (s).
    pub final_delta: f64,
}
impl LapDelta {
    /// Compute a lap delta from two sequences of (position, time) pairs.
    ///
    /// Both sequences must be sorted by position and of equal length.
    pub fn from_segments(positions: Vec<f64>, ref_times: &[f64], cmp_times: &[f64]) -> Self {
        assert_eq!(ref_times.len(), cmp_times.len());
        let delta: Vec<f64> = ref_times
            .iter()
            .zip(cmp_times.iter())
            .map(|(r, c)| c - r)
            .collect();
        let final_delta = delta.last().copied().unwrap_or(0.0);
        Self {
            positions,
            delta,
            final_delta,
        }
    }
    /// Position at which the comparison car first gains time (delta goes negative).
    pub fn first_gain_position(&self) -> Option<f64> {
        for (i, &d) in self.delta.iter().enumerate() {
            if d < 0.0 {
                return self.positions.get(i).copied();
            }
        }
        None
    }
    /// Maximum gain achieved by comparison (most negative delta, s).
    pub fn max_gain(&self) -> f64 {
        self.delta.iter().cloned().fold(0.0_f64, f64::min)
    }
    /// Maximum loss suffered by comparison (most positive delta, s).
    pub fn max_loss(&self) -> f64 {
        self.delta.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Number of sample points.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Returns `true` if no sample points are stored.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}
/// Find the fastest lap from a collection of lap results.
#[derive(Debug, Clone)]
pub struct FastestLapFinder {
    /// Lap times collected (s).
    pub lap_times: Vec<f64>,
}
impl FastestLapFinder {
    /// Create an empty finder.
    pub fn new() -> Self {
        Self {
            lap_times: Vec::new(),
        }
    }
    /// Add a lap time (s).
    pub fn add_lap(&mut self, time: f64) {
        self.lap_times.push(time);
    }
    /// Return `(lap_index, lap_time)` of the fastest lap (0-based index).
    pub fn fastest(&self) -> Option<(usize, f64)> {
        self.lap_times
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, &t)| (i, t))
    }
    /// Mean lap time (s).
    pub fn mean_lap_time(&self) -> f64 {
        if self.lap_times.is_empty() {
            return 0.0;
        }
        self.lap_times.iter().sum::<f64>() / self.lap_times.len() as f64
    }
    /// Standard deviation of lap times (s).
    pub fn lap_time_std_dev(&self) -> f64 {
        if self.lap_times.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_lap_time();
        let variance = self
            .lap_times
            .iter()
            .map(|&t| (t - mean).powi(2))
            .sum::<f64>()
            / self.lap_times.len() as f64;
        variance.sqrt()
    }
    /// Slowest lap (0-based index, lap_time) or `None` if empty.
    pub fn slowest(&self) -> Option<(usize, f64)> {
        self.lap_times
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, &t)| (i, t))
    }
    /// Number of laps recorded.
    pub fn lap_count(&self) -> usize {
        self.lap_times.len()
    }
}
/// Full race pit stop strategy.
#[derive(Debug, Clone)]
pub struct PitStopStrategy {
    /// Ordered list of stints.
    pub stints: Vec<Stint>,
    /// Time lost per pit stop (s).
    pub pit_loss_time: f64,
    /// Total race laps.
    pub total_laps: usize,
}
impl PitStopStrategy {
    /// Create a 1-stop strategy: soft → hard.
    pub fn one_stop_soft_hard(total_laps: usize, pit_loss: f64, pit_lap: usize) -> Self {
        let pit_lap = pit_lap.clamp(1, total_laps - 1);
        Self {
            stints: vec![
                Stint::new(TireCompound::Soft, pit_lap, 1),
                Stint::new(TireCompound::Hard, total_laps - pit_lap, pit_lap + 1),
            ],
            pit_loss_time: pit_loss,
            total_laps,
        }
    }
    /// Create a 2-stop strategy: soft → medium → soft.
    pub fn two_stop_soft_medium_soft(
        total_laps: usize,
        pit_loss: f64,
        pit1: usize,
        pit2: usize,
    ) -> Self {
        let p1 = pit1.clamp(1, total_laps - 2);
        let p2 = pit2.clamp(p1 + 1, total_laps - 1);
        Self {
            stints: vec![
                Stint::new(TireCompound::Soft, p1, 1),
                Stint::new(TireCompound::Medium, p2 - p1, p1 + 1),
                Stint::new(TireCompound::Soft, total_laps - p2, p2 + 1),
            ],
            pit_loss_time: pit_loss,
            total_laps,
        }
    }
    /// Estimate total race time for a given base lap time (s).
    pub fn estimated_race_time(&self, base_laptime: f64) -> f64 {
        let mut total = 0.0;
        for stint in &self.stints {
            let pace_adj = stint.pace_advantage_per_lap();
            for lap_in_stint in 0..stint.laps {
                let deg = lap_in_stint as f64 * stint.compound.degradation_rate() * 2.0;
                total += base_laptime + pace_adj + deg;
            }
        }
        let pit_count = self.stints.len().saturating_sub(1);
        total += self.pit_loss_time * pit_count as f64;
        total
    }
    /// Find the optimal single pit stop lap to minimise race time.
    ///
    /// Sweeps pit lap from `min_lap` to `max_lap` inclusive.
    /// Returns `(best_pit_lap, best_race_time)`.
    pub fn optimise_one_stop(
        total_laps: usize,
        base_laptime: f64,
        pit_loss: f64,
        min_lap: usize,
        max_lap: usize,
    ) -> (usize, f64) {
        let mut best_lap = min_lap;
        let mut best_time = f64::INFINITY;
        for pit_lap in min_lap..=max_lap.min(total_laps - 1) {
            let strategy = PitStopStrategy::one_stop_soft_hard(total_laps, pit_loss, pit_lap);
            let t = strategy.estimated_race_time(base_laptime);
            if t < best_time {
                best_time = t;
                best_lap = pit_lap;
            }
        }
        (best_lap, best_time)
    }
    /// Number of pit stops planned.
    pub fn pit_count(&self) -> usize {
        self.stints.len().saturating_sub(1)
    }
    /// Compound used on lap `lap_number` (1-indexed).
    pub fn compound_at_lap(&self, lap_number: usize) -> Option<TireCompound> {
        for stint in &self.stints {
            if lap_number >= stint.start_lap && lap_number <= stint.end_lap() {
                return Some(stint.compound);
            }
        }
        None
    }
}
/// G-G diagram: combined acceleration envelope.
#[derive(Debug, Clone)]
pub struct GGDiagram {
    /// Lateral acceleration axis samples (m/s²) — symmetric about zero.
    pub lat_samples: Vec<f64>,
    /// Longitudinal acceleration axis samples (m/s²).
    pub lon_samples: Vec<f64>,
    /// Feasible region boundary: `envelope[i]` = max |a_lon| at `lat_samples[i]`.
    pub envelope: Vec<f64>,
    /// Speed at which this diagram was computed (m/s).
    pub speed: f64,
}
impl GGDiagram {
    /// Build a G-G diagram at the given speed for a vehicle performance map.
    ///
    /// The friction circle is approximated with an ellipse:
    /// `(a_lat/a_lat_max)² + (a_lon/a_lon_max)² ≤ 1`.
    pub fn build(perf: &VehiclePerformanceMap, speed: f64, n: usize) -> Self {
        let a_lat_max = perf.lat_g_at(speed);
        let a_lon_accel = perf.traction_g_at(speed);
        let a_lon_brake = perf.brake_g_at(speed);
        let mut lat_samples = Vec::with_capacity(2 * n + 1);
        let mut envelope = Vec::with_capacity(2 * n + 1);
        let mut lon_samples = Vec::with_capacity(2 * n + 1);
        for i in 0..=(2 * n) {
            let a_lat = -a_lat_max + (2.0 * a_lat_max * i as f64) / (2 * n) as f64;
            let lat_frac_sq = (a_lat / a_lat_max.max(1e-12)).powi(2);
            let remaining = (1.0 - lat_frac_sq).max(0.0).sqrt();
            let max_accel = a_lon_accel * remaining;
            let max_brake = -a_lon_brake * remaining;
            lat_samples.push(a_lat);
            lon_samples.push(max_accel);
            envelope.push(max_accel.abs().max(max_brake.abs()));
        }
        Self {
            lat_samples,
            lon_samples,
            envelope,
            speed,
        }
    }
    /// Check whether a `(a_lat, a_lon)` combination is within the envelope.
    pub fn is_feasible(&self, a_lat: f64, a_lon: f64) -> bool {
        if self.lat_samples.is_empty() {
            return false;
        }
        let mut best_idx = 0;
        let mut best_dist = f64::INFINITY;
        for (i, &al) in self.lat_samples.iter().enumerate() {
            let d = (al - a_lat).abs();
            if d < best_dist {
                best_dist = d;
                best_idx = i;
            }
        }
        let max_lon = self.envelope[best_idx];
        a_lon.abs() <= max_lon
    }
    /// Maximum combined acceleration magnitude (m/s²) at index `i`.
    pub fn combined_magnitude(&self, i: usize) -> f64 {
        if i >= self.lat_samples.len() {
            return 0.0;
        }
        (self.lat_samples[i].powi(2) + self.lon_samples[i].powi(2)).sqrt()
    }
    /// Number of points in the G-G diagram.
    pub fn len(&self) -> usize {
        self.lat_samples.len()
    }
    /// Returns `true` if the G-G diagram has no samples.
    pub fn is_empty(&self) -> bool {
        self.lat_samples.is_empty()
    }
}
/// Tire compound type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TireCompound {
    /// Soft compound: fast but high degradation.
    Soft,
    /// Medium compound: balanced.
    Medium,
    /// Hard compound: durable but slower.
    Hard,
    /// Intermediate wet compound.
    Intermediate,
    /// Full wet compound.
    Wet,
}
impl TireCompound {
    /// Peak grip multiplier relative to medium.
    pub fn peak_grip(&self) -> f64 {
        match self {
            TireCompound::Soft => 1.06,
            TireCompound::Medium => 1.0,
            TireCompound::Hard => 0.97,
            TireCompound::Intermediate => 0.85,
            TireCompound::Wet => 0.75,
        }
    }
    /// Degradation rate per lap (fraction of grip lost per lap at full load).
    pub fn degradation_rate(&self) -> f64 {
        match self {
            TireCompound::Soft => 0.035,
            TireCompound::Medium => 0.018,
            TireCompound::Hard => 0.010,
            TireCompound::Intermediate => 0.020,
            TireCompound::Wet => 0.012,
        }
    }
    /// Optimal operating temperature window (°C): (min, max).
    pub fn temperature_window(&self) -> (f64, f64) {
        match self {
            TireCompound::Soft => (85.0, 105.0),
            TireCompound::Medium => (90.0, 120.0),
            TireCompound::Hard => (95.0, 135.0),
            TireCompound::Intermediate => (50.0, 90.0),
            TireCompound::Wet => (30.0, 70.0),
        }
    }
    /// Maximum recommended stint length (laps).
    pub fn max_stint_laps(&self) -> usize {
        match self {
            TireCompound::Soft => 18,
            TireCompound::Medium => 32,
            TireCompound::Hard => 50,
            TireCompound::Intermediate => 20,
            TireCompound::Wet => 30,
        }
    }
}
/// Minimum-time manoeuvring geometry and limits.
#[derive(Debug, Clone)]
pub struct OptimalLap {
    /// Vehicle performance map.
    pub perf: VehiclePerformanceMap,
    /// Driver model.
    pub driver: DriverModel,
}
impl OptimalLap {
    /// Construct with default F1 performance and experienced driver.
    pub fn new_f1() -> Self {
        Self {
            perf: VehiclePerformanceMap::default_f1(),
            driver: DriverModel::experienced(),
        }
    }
    /// Maximum geometric corner speed (m/s) for a given radius (m).
    pub fn v_corner(&self, radius: f64) -> f64 {
        let v_uncapped =
            (self.perf.lat_g_at(50.0) * self.driver.effective_grip_fraction() * radius).sqrt();
        v_uncapped.min(self.perf.v_max)
    }
    /// Braking distance from `v_entry` (m/s) to `v_corner` (m/s).
    pub fn braking_distance(&self, v_entry: f64, v_exit: f64) -> f64 {
        let a_brake = self.perf.brake_g_at(v_entry);
        if a_brake < 1e-3 || v_entry <= v_exit {
            return 0.0;
        }
        (v_entry * v_entry - v_exit * v_exit) / (2.0 * a_brake)
    }
    /// Acceleration distance from `v_exit` to `v_max` limited by segment length.
    pub fn accel_distance(&self, v_exit: f64, v_max: f64, available_length: f64) -> f64 {
        let a_trac = self.perf.traction_g_at(v_exit);
        if a_trac < 1e-3 {
            return available_length;
        }
        let d = (v_max * v_max - v_exit * v_exit) / (2.0 * a_trac);
        d.min(available_length)
    }
    /// Time to travel `distance` (m) at constant speed `v` (m/s).
    pub fn time_at_speed(&self, distance: f64, speed: f64) -> f64 {
        if speed < 1e-3 {
            return f64::INFINITY;
        }
        distance / speed
    }
    /// Estimate segment time (s) for a straight of given length and entry speed.
    pub fn straight_time(&self, length: f64, v_entry: f64, v_limit: f64) -> f64 {
        let a = self.perf.traction_g_at(v_entry);
        let v_lim = v_limit.min(self.perf.v_max);
        let d_accel = if v_entry < v_lim && a > 1e-3 {
            ((v_lim * v_lim - v_entry * v_entry) / (2.0 * a)).min(length)
        } else {
            0.0
        };
        let t_accel = if a > 1e-3 && v_entry < v_lim {
            (v_lim - v_entry) / a
        } else {
            0.0
        };
        let d_const = (length - d_accel).max(0.0);
        let v_avg = (v_entry + v_lim) / 2.0;
        let t_const = if v_avg > 1e-3 { d_const / v_avg } else { 0.0 };
        t_accel + t_const
    }
    /// Estimate corner time (s) at geometric speed.
    pub fn corner_time(&self, radius: f64, angle: f64) -> f64 {
        let v = self.v_corner(radius);
        let arc_length = radius * angle.abs();
        self.time_at_speed(arc_length, v)
    }
}
/// Comparison between two lap results.
#[derive(Debug, Clone)]
pub struct LapCompare {
    /// Reference lap.
    pub reference: LapResult,
    /// Comparison lap.
    pub comparison: LapResult,
}
impl LapCompare {
    /// Construct a comparison.
    pub fn new(reference: LapResult, comparison: LapResult) -> Self {
        Self {
            reference,
            comparison,
        }
    }
    /// Overall lap time delta (comparison − reference, s). Negative = comparison faster.
    pub fn lap_time_delta(&self) -> f64 {
        self.comparison.lap_time - self.reference.lap_time
    }
    /// Sector time deltas (comparison − reference, s) for each sector.
    pub fn sector_deltas(&self) -> Vec<f64> {
        self.reference
            .sector_times
            .iter()
            .zip(self.comparison.sector_times.iter())
            .map(|(r, c)| c - r)
            .collect()
    }
    /// Speed delta (m/s): comparison max speed − reference max speed.
    pub fn max_speed_delta(&self) -> f64 {
        self.comparison.max_speed - self.reference.max_speed
    }
    /// Tire degradation delta.
    pub fn tire_delta(&self) -> f64 {
        self.comparison.tire_degradation - self.reference.tire_degradation
    }
    /// Returns `true` if the comparison lap is faster.
    pub fn is_comparison_faster(&self) -> bool {
        self.lap_time_delta() < 0.0
    }
}
/// A single segment of a race track.
#[derive(Debug, Clone)]
pub enum TrackSegment {
    /// Straight section.
    Straight {
        /// Length of the straight (m).
        length: f64,
        /// DRS activation zone flag.
        drs_available: bool,
    },
    /// Constant-radius corner.
    Corner {
        /// Corner radius (m).
        radius: f64,
        /// Subtended angle (rad).
        angle: f64,
        /// Direction: positive = left-hand, negative = right-hand.
        direction: f64,
    },
    /// Chicane — two or more alternating corners.
    Chicane {
        /// Radii of each mini-corner (m).
        radii: Vec<f64>,
        /// Angles of each mini-corner (rad).
        angles: Vec<f64>,
    },
    /// Elevation change section.
    Elevation {
        /// Length along the road (m).
        length: f64,
        /// Gradient (rise/run, positive = uphill).
        gradient: f64,
    },
}
impl TrackSegment {
    /// Length of this segment along the road centreline (m).
    pub fn length(&self) -> f64 {
        match self {
            TrackSegment::Straight { length, .. } => *length,
            TrackSegment::Corner { radius, angle, .. } => radius * angle.abs(),
            TrackSegment::Chicane { radii, angles } => radii
                .iter()
                .zip(angles.iter())
                .map(|(r, a)| r * a.abs())
                .sum(),
            TrackSegment::Elevation { length, .. } => *length,
        }
    }
    /// Minimum geometric corner radius (m). Returns `None` for non-corner segments.
    pub fn min_radius(&self) -> Option<f64> {
        match self {
            TrackSegment::Corner { radius, .. } => Some(*radius),
            TrackSegment::Chicane { radii, .. } => radii.iter().cloned().reduce(f64::min),
            _ => None,
        }
    }
    /// Returns `true` if this segment is a DRS zone.
    pub fn is_drs_zone(&self) -> bool {
        matches!(
            self,
            TrackSegment::Straight {
                drs_available: true,
                ..
            }
        )
    }
}
/// KERS/ERS energy recovery and deployment model.
#[derive(Debug, Clone)]
pub struct EnergyRecovery {
    /// Maximum energy storage capacity (J).
    pub capacity_j: f64,
    /// Current state of charge (J).
    pub soc_j: f64,
    /// Maximum deployment power (W).
    pub max_deploy_power: f64,
    /// Maximum harvest power (W).
    pub max_harvest_power: f64,
    /// Total energy budget per lap (J).
    pub lap_energy_budget: f64,
    /// Energy deployed this lap (J).
    pub deployed_this_lap: f64,
    /// Energy harvested this lap (J).
    pub harvested_this_lap: f64,
}
impl EnergyRecovery {
    /// Default MGU-K model for a modern F1 car.
    pub fn default_f1_mguk() -> Self {
        Self {
            capacity_j: 4_000_000.0,
            soc_j: 4_000_000.0,
            max_deploy_power: 120_000.0,
            max_harvest_power: 120_000.0,
            lap_energy_budget: 4_000_000.0,
            deployed_this_lap: 0.0,
            harvested_this_lap: 0.0,
        }
    }
    /// Deploy `power` (W) for `dt` (s). Returns actual energy deployed (J).
    pub fn deploy(&mut self, power: f64, dt: f64) -> f64 {
        let p = power.min(self.max_deploy_power);
        let e = (p * dt)
            .min(self.soc_j)
            .min(self.lap_energy_budget - self.deployed_this_lap);
        let e = e.max(0.0);
        self.soc_j -= e;
        self.deployed_this_lap += e;
        e
    }
    /// Harvest `power` (W) for `dt` (s). Returns actual energy harvested (J).
    pub fn harvest(&mut self, power: f64, dt: f64) -> f64 {
        let p = power.min(self.max_harvest_power);
        let e = (p * dt).min(self.capacity_j - self.soc_j);
        let e = e.max(0.0);
        self.soc_j += e;
        self.harvested_this_lap += e;
        e
    }
    /// State of charge as a fraction (0–1).
    pub fn soc_fraction(&self) -> f64 {
        self.soc_j / self.capacity_j.max(1.0)
    }
    /// Reset lap counters at the start of a new lap.
    pub fn reset_lap(&mut self) {
        self.deployed_this_lap = 0.0;
        self.harvested_this_lap = 0.0;
    }
    /// Returns `true` if the energy budget for this lap has been exhausted.
    pub fn budget_exhausted(&self) -> bool {
        self.deployed_this_lap >= self.lap_energy_budget
    }
}
/// Configuration for a full lap simulation run.
#[derive(Debug, Clone)]
pub struct LapSimConfig {
    /// Vehicle performance map.
    pub perf: VehiclePerformanceMap,
    /// Driver model.
    pub driver: DriverModel,
    /// Fuel strategy.
    pub fuel: FuelEffect,
    /// Energy recovery model.
    pub ers: EnergyRecovery,
    /// Minimum corner speed clamp (m/s).
    pub min_speed: f64,
}
impl LapSimConfig {
    /// Default F1 lap simulation configuration.
    pub fn default_f1() -> Self {
        Self {
            perf: VehiclePerformanceMap::default_f1(),
            driver: DriverModel::experienced(),
            fuel: FuelEffect::default_f1_one_stop(),
            ers: EnergyRecovery::default_f1_mguk(),
            min_speed: 10.0,
        }
    }
}
/// A single stint description.
#[derive(Debug, Clone)]
pub struct Stint {
    /// Tire compound used.
    pub compound: TireCompound,
    /// Number of laps in this stint.
    pub laps: usize,
    /// Lap number this stint starts on (1-indexed).
    pub start_lap: usize,
}
impl Stint {
    /// Create a new stint.
    pub fn new(compound: TireCompound, laps: usize, start_lap: usize) -> Self {
        Self {
            compound,
            laps,
            start_lap,
        }
    }
    /// End lap of this stint (inclusive).
    pub fn end_lap(&self) -> usize {
        self.start_lap + self.laps - 1
    }
    /// Estimated stint pace advantage relative to medium (s/lap).
    pub fn pace_advantage_per_lap(&self) -> f64 {
        match self.compound {
            TireCompound::Soft => -0.4,
            TireCompound::Medium => 0.0,
            TireCompound::Hard => 0.25,
            TireCompound::Intermediate => 5.0,
            TireCompound::Wet => 8.0,
        }
    }
}
/// Map of ERS deploy/harvest zones over the track.
#[derive(Debug, Clone)]
pub struct ErsZoneMap {
    /// Zone classification per segment.
    pub zones: Vec<ErsZoneType>,
    /// Energy available for deployment this lap (J).
    pub available_j: f64,
    /// Energy harvested this lap (J).
    pub harvested_j: f64,
}
impl ErsZoneMap {
    /// Build an ERS zone map from a track layout.
    ///
    /// Straights are deploy zones; pre-corner zones are harvest zones.
    pub fn from_track(track: &TrackLayout) -> Self {
        let n = track.segments.len();
        let mut zones = vec![ErsZoneType::Neutral; n];
        for (i, seg) in track.segments.iter().enumerate() {
            match seg {
                TrackSegment::Straight { .. } => {
                    zones[i] = ErsZoneType::Deploy;
                }
                TrackSegment::Corner { .. } | TrackSegment::Chicane { .. } => {
                    if i > 0 {
                        zones[i - 1] = ErsZoneType::Harvest;
                    }
                    zones[i] = ErsZoneType::Neutral;
                }
                TrackSegment::Elevation { gradient, .. } => {
                    zones[i] = if *gradient < 0.0 {
                        ErsZoneType::Harvest
                    } else {
                        ErsZoneType::Deploy
                    };
                }
            }
        }
        Self {
            zones,
            available_j: 4_000_000.0,
            harvested_j: 0.0,
        }
    }
    /// Number of deploy zones.
    pub fn deploy_zone_count(&self) -> usize {
        self.zones
            .iter()
            .filter(|&&z| z == ErsZoneType::Deploy)
            .count()
    }
    /// Number of harvest zones.
    pub fn harvest_zone_count(&self) -> usize {
        self.zones
            .iter()
            .filter(|&&z| z == ErsZoneType::Harvest)
            .count()
    }
    /// Zone type at segment index `i`.
    pub fn zone_at(&self, i: usize) -> ErsZoneType {
        self.zones.get(i).copied().unwrap_or(ErsZoneType::Neutral)
    }
    /// Reset lap energy counters.
    pub fn reset(&mut self) {
        self.available_j = 4_000_000.0;
        self.harvested_j = 0.0;
    }
}
/// Result of a single lap simulation.
#[derive(Debug, Clone)]
pub struct LapResult {
    /// Total lap time (s).
    pub lap_time: f64,
    /// Sector times (s) — one entry per sector.
    pub sector_times: Vec<f64>,
    /// Maximum speed achieved (m/s).
    pub max_speed: f64,
    /// Average speed (m/s).
    pub avg_speed: f64,
    /// Tire degradation index (0 = new, 1 = fully worn).
    pub tire_degradation: f64,
    /// Fuel consumed this lap (kg).
    pub fuel_consumption: f64,
    /// ERS energy deployed this lap (J).
    pub ers_deployed: f64,
}
impl LapResult {
    /// Lap time in minutes:seconds.milliseconds string.
    pub fn lap_time_str(&self) -> String {
        let minutes = (self.lap_time / 60.0) as u32;
        let seconds = self.lap_time - minutes as f64 * 60.0;
        format!("{minutes}:{seconds:06.3}")
    }
    /// Return sector time for sector index (0-based).
    pub fn sector_time(&self, sector: usize) -> Option<f64> {
        self.sector_times.get(sector).copied()
    }
}
/// Input conditions for compound selection.
#[derive(Debug, Clone)]
pub struct RaceConditions {
    /// Track surface temperature (°C).
    pub track_temp_c: f64,
    /// Remaining race laps.
    pub remaining_laps: usize,
    /// Whether the track is wet.
    pub wet: bool,
    /// Current stint lap count (laps already done on current set).
    pub laps_on_current: usize,
}
/// Recommendation from the compound optimiser.
#[derive(Debug, Clone)]
pub struct CompoundRecommendation {
    /// Recommended tire compound.
    pub compound: TireCompound,
    /// Estimated time advantage over medium (s total over remaining stint).
    pub time_advantage_s: f64,
    /// Recommended stint length (laps).
    pub recommended_stint: usize,
}
/// Compound selection optimiser.
pub struct CompoundOptimiser;
impl CompoundOptimiser {
    /// Select the best tire compound given race conditions.
    pub fn select(conditions: &RaceConditions) -> CompoundRecommendation {
        if conditions.wet {
            let compound = if conditions.track_temp_c < 15.0 {
                TireCompound::Wet
            } else {
                TireCompound::Intermediate
            };
            return CompoundRecommendation {
                compound,
                time_advantage_s: 0.0,
                recommended_stint: compound.max_stint_laps(),
            };
        }
        let candidates = [TireCompound::Soft, TireCompound::Medium, TireCompound::Hard];
        let mut best_compound = TireCompound::Medium;
        let mut best_advantage = f64::NEG_INFINITY;
        for &compound in &candidates {
            let stint_len = compound.max_stint_laps().min(conditions.remaining_laps);
            let pace_adv = -Stint::new(compound, 1, 1).pace_advantage_per_lap();
            let mut total_adv = 0.0;
            for lap in 0..stint_len {
                let deg_penalty = lap as f64 * compound.degradation_rate() * 1.5;
                total_adv += pace_adv - deg_penalty;
            }
            if total_adv > best_advantage {
                best_advantage = total_adv;
                best_compound = compound;
            }
        }
        CompoundRecommendation {
            compound: best_compound,
            time_advantage_s: best_advantage,
            recommended_stint: best_compound
                .max_stint_laps()
                .min(conditions.remaining_laps),
        }
    }
}
/// Result of a multi-lap race simulation.
#[derive(Debug, Clone)]
pub struct RaceSimResult {
    /// Lap times for each lap (s).
    pub lap_times: Vec<f64>,
    /// Tire degradation at end of each lap.
    pub tire_degradation: Vec<f64>,
    /// Fuel mass at start of each lap (kg).
    pub fuel_masses: Vec<f64>,
    /// ERS deployed per lap (J).
    pub ers_deployed: Vec<f64>,
    /// Total race time (s).
    pub total_race_time: f64,
    /// Lap number of fastest lap (1-indexed).
    pub fastest_lap_number: usize,
    /// Fastest lap time (s).
    pub fastest_lap_time: f64,
}
impl RaceSimResult {
    /// Average lap time (s).
    pub fn average_lap_time(&self) -> f64 {
        if self.lap_times.is_empty() {
            return 0.0;
        }
        self.lap_times.iter().sum::<f64>() / self.lap_times.len() as f64
    }
    /// Lap time spread (max − min, s).
    pub fn lap_time_spread(&self) -> f64 {
        if self.lap_times.is_empty() {
            return 0.0;
        }
        let max = self
            .lap_times
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min = self.lap_times.iter().cloned().fold(f64::INFINITY, f64::min);
        max - min
    }
}
/// Lap-by-lap tire degradation state.
#[derive(Debug, Clone)]
pub struct TireDegradationState {
    /// Tire compound in use.
    pub compound: TireCompound,
    /// Current degradation level (0 = new, 1 = cliff).
    pub degradation: f64,
    /// Laps on this set of tires.
    pub laps_on_tire: usize,
    /// Track temperature (°C) — affects degradation rate.
    pub track_temp: f64,
}
impl TireDegradationState {
    /// Create a fresh tire set.
    pub fn new(compound: TireCompound, track_temp: f64) -> Self {
        Self {
            compound,
            degradation: 0.0,
            laps_on_tire: 0,
            track_temp,
        }
    }
    /// Effective grip fraction (0–1) after degradation.
    pub fn grip_fraction(&self) -> f64 {
        let base = self.compound.peak_grip();
        let degrade_factor = (1.0 - self.degradation * 0.8).max(0.0);
        (base * degrade_factor).min(1.15)
    }
    /// Laptime delta (s) from tire degradation relative to fresh tire.
    ///
    /// A fully degraded soft at cliff contributes ~1.5 s/lap.
    pub fn laptime_delta_s(&self) -> f64 {
        let grip_loss = 1.0 - self.grip_fraction() / self.compound.peak_grip();
        grip_loss * 3.0
    }
    /// Advance degradation by one lap. Returns degradation added.
    pub fn advance_lap(&mut self) -> f64 {
        let temp_factor = {
            let (t_min, t_max) = self.compound.temperature_window();
            let t_mid = (t_min + t_max) * 0.5;
            let deviation = (self.track_temp - t_mid).abs() / (t_max - t_min).max(1.0);
            1.0 + deviation * 0.5
        };
        let rate = self.compound.degradation_rate() * temp_factor;
        let added = rate.min(1.0 - self.degradation);
        self.degradation = (self.degradation + added).min(1.0);
        self.laps_on_tire += 1;
        added
    }
    /// Returns `true` if tires have reached the cliff (degradation > 0.8).
    pub fn at_cliff(&self) -> bool {
        self.degradation > 0.8
    }
}
/// Fuel mass and strategy model.
#[derive(Debug, Clone)]
pub struct FuelEffect {
    /// Starting fuel mass (kg).
    pub start_fuel: f64,
    /// Fuel consumption rate (kg/lap).
    pub fuel_per_lap: f64,
    /// Laptime sensitivity to fuel mass (s/kg).
    pub laptime_per_kg: f64,
    /// Number of pit stops planned.
    pub pit_stop_count: usize,
    /// Time lost per pit stop (s).
    pub pit_stop_time: f64,
}
impl FuelEffect {
    /// Default F1 race fuel strategy (50-lap race, 1-stop).
    pub fn default_f1_one_stop() -> Self {
        Self {
            start_fuel: 110.0,
            fuel_per_lap: 2.2,
            laptime_per_kg: 0.035,
            pit_stop_count: 1,
            pit_stop_time: 22.0,
        }
    }
    /// Default 2-stop strategy.
    pub fn default_f1_two_stop() -> Self {
        Self {
            start_fuel: 105.0,
            fuel_per_lap: 2.2,
            laptime_per_kg: 0.035,
            pit_stop_count: 2,
            pit_stop_time: 22.0,
        }
    }
    /// Fuel mass at the start of lap `lap_number` (1-indexed).
    pub fn fuel_at_lap(&self, lap_number: usize) -> f64 {
        let consumed = self.fuel_per_lap * (lap_number.saturating_sub(1)) as f64;
        (self.start_fuel - consumed).max(0.0)
    }
    /// Laptime penalty from fuel mass relative to a zero-fuel baseline (s).
    pub fn laptime_penalty(&self, fuel_kg: f64) -> f64 {
        fuel_kg * self.laptime_per_kg
    }
    /// Estimated total race time for `total_laps` laps (s), including pit stops.
    pub fn total_race_time(&self, base_laptime: f64, total_laps: usize) -> f64 {
        let mut total = 0.0;
        for lap in 1..=total_laps {
            let fuel = self.fuel_at_lap(lap);
            total += base_laptime + self.laptime_penalty(fuel);
        }
        total += self.pit_stop_time * self.pit_stop_count as f64;
        total
    }
    /// Undercut potential: time saved by pitting one lap earlier than opponent.
    pub fn undercut_gain(&self, current_lap: usize, opponent_lap: usize) -> f64 {
        let my_fuel_now = self.fuel_at_lap(current_lap);
        let opp_fuel_now = self.fuel_at_lap(opponent_lap);
        self.laptime_penalty(opp_fuel_now) - self.laptime_penalty(my_fuel_now)
    }
}
/// Classifies each track segment as a harvest or deploy zone for ERS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErsZoneType {
    /// Deploy ERS power on this segment.
    Deploy,
    /// Harvest energy by regenerative braking on this segment.
    Harvest,
    /// Neither (mid-corner, elevation).
    Neutral,
}
/// Detected apex information for a corner.
#[derive(Debug, Clone)]
pub struct CornerApex {
    /// Segment index in the track layout.
    pub segment_index: usize,
    /// Apex speed (m/s).
    pub apex_speed: f64,
    /// Distance from segment start to apex along centreline (m).
    pub apex_distance: f64,
    /// Corner radius at apex (m).
    pub radius: f64,
    /// Braking zone start distance before apex (m).
    pub braking_zone_start: f64,
    /// Acceleration zone end distance after apex (m).
    pub accel_zone_end: f64,
}
impl CornerApex {
    /// Detect all corner apexes in a track layout using the optimal lap model.
    pub fn detect_all(track: &TrackLayout, opt: &OptimalLap) -> Vec<Self> {
        let mut apexes = Vec::new();
        let mut cumulative_dist = 0.0_f64;
        for (idx, seg) in track.segments.iter().enumerate() {
            match seg {
                TrackSegment::Corner { radius, angle, .. } => {
                    let apex_speed = opt.v_corner(*radius);
                    let arc = radius * angle.abs();
                    let apex_dist = arc * 0.5;
                    let braking_zone = opt.braking_distance(opt.perf.v_max, apex_speed);
                    let accel_zone = opt.accel_distance(apex_speed, opt.perf.v_max, 300.0);
                    apexes.push(CornerApex {
                        segment_index: idx,
                        apex_speed,
                        apex_distance: cumulative_dist + apex_dist,
                        radius: *radius,
                        braking_zone_start: (cumulative_dist - braking_zone).max(0.0),
                        accel_zone_end: cumulative_dist + arc + accel_zone,
                    });
                }
                TrackSegment::Chicane { radii, angles } => {
                    for (i, (r, a)) in radii.iter().zip(angles.iter()).enumerate() {
                        let apex_speed = opt.v_corner(*r);
                        let arc = r * a.abs();
                        let apex_dist = arc * 0.5;
                        let braking_zone = opt.braking_distance(opt.perf.v_max, apex_speed);
                        let accel_zone = opt.accel_distance(apex_speed, opt.perf.v_max, 200.0);
                        apexes.push(CornerApex {
                            segment_index: idx,
                            apex_speed,
                            apex_distance: cumulative_dist
                                + radii[..i]
                                    .iter()
                                    .zip(angles[..i].iter())
                                    .map(|(rr, aa)| rr * aa.abs())
                                    .sum::<f64>()
                                + apex_dist,
                            radius: *r,
                            braking_zone_start: (cumulative_dist - braking_zone).max(0.0),
                            accel_zone_end: cumulative_dist + arc + accel_zone,
                        });
                    }
                }
                _ => {}
            }
            cumulative_dist += seg.length();
        }
        apexes
    }
    /// Returns `true` if a given track position (m) is within the braking zone.
    pub fn in_braking_zone(&self, track_pos: f64) -> bool {
        track_pos >= self.braking_zone_start && track_pos < self.apex_distance
    }
    /// Returns `true` if a given track position (m) is within the acceleration zone.
    pub fn in_accel_zone(&self, track_pos: f64) -> bool {
        track_pos > self.apex_distance && track_pos <= self.accel_zone_end
    }
}
/// Racing line optimised by minimum curvature (sum of 1/r²).
#[derive(Debug, Clone)]
pub struct RacingLine {
    /// Track width available (m).
    pub track_width: f64,
    /// Number of optimisation iterations.
    pub iterations: usize,
    /// Curvature cost weights per segment (1/m²).
    pub curvature_costs: Vec<f64>,
    /// Effective radii after optimisation (m).
    pub effective_radii: Vec<f64>,
}
impl RacingLine {
    /// Initialise a racing line from a track layout with given track width.
    pub fn new(track: &TrackLayout, track_width: f64) -> Self {
        let n = track.segments.len();
        let mut curvature_costs = vec![0.0_f64; n];
        let mut effective_radii = vec![f64::INFINITY; n];
        for (i, seg) in track.segments.iter().enumerate() {
            match seg {
                TrackSegment::Corner { radius, .. } => {
                    let r_opt = radius + track_width * 0.5;
                    effective_radii[i] = r_opt;
                    curvature_costs[i] = 1.0 / (r_opt * r_opt);
                }
                TrackSegment::Chicane { radii, .. } => {
                    let r_min = radii.iter().cloned().fold(f64::INFINITY, f64::min);
                    let r_opt = r_min + track_width * 0.3;
                    effective_radii[i] = r_opt;
                    curvature_costs[i] = 1.0 / (r_opt * r_opt);
                }
                _ => {
                    curvature_costs[i] = 0.0;
                }
            }
        }
        Self {
            track_width,
            iterations: 0,
            curvature_costs,
            effective_radii,
        }
    }
    /// Total curvature cost (sum of 1/r² over corners).
    pub fn total_curvature_cost(&self) -> f64 {
        self.curvature_costs.iter().sum()
    }
    /// Perform one round of gradient descent to reduce curvature.
    /// Returns the change in total curvature cost.
    pub fn optimise_step(&mut self, step_size: f64) -> f64 {
        let before = self.total_curvature_cost();
        for i in 0..self.effective_radii.len() {
            if self.effective_radii[i].is_finite() {
                let r = self.effective_radii[i];
                let gradient = -2.0 / (r * r * r);
                let new_r = (r - step_size * gradient).max(1.0);
                self.effective_radii[i] = new_r;
                self.curvature_costs[i] = 1.0 / (new_r * new_r);
            }
        }
        self.iterations += 1;
        let after = self.total_curvature_cost();
        before - after
    }
    /// Run `n` optimisation steps. Returns final total curvature cost.
    pub fn optimise(&mut self, n: usize, step_size: f64) -> f64 {
        for _ in 0..n {
            self.optimise_step(step_size);
        }
        self.total_curvature_cost()
    }
    /// Effective corner speed (m/s) at segment `i` after optimisation.
    pub fn effective_corner_speed(&self, seg_idx: usize, opt: &OptimalLap) -> f64 {
        let r = self
            .effective_radii
            .get(seg_idx)
            .copied()
            .unwrap_or(f64::INFINITY);
        if r.is_finite() {
            opt.v_corner(r)
        } else {
            opt.perf.v_max
        }
    }
}
/// Driver model parameters that modify the theoretical optimal lap.
#[derive(Debug, Clone)]
pub struct DriverModel {
    /// Braking point offset (m) — positive = later braking.
    pub braking_point_offset: f64,
    /// Throttle application commitment (0–1, 1 = instant full throttle).
    pub throttle_application: f64,
    /// Cornering commitment (0–1, 1 = uses full lateral limit).
    pub cornering_commitment: f64,
    /// Random error sigma (s per corner — lap time variability).
    pub error_sigma: f64,
}
impl DriverModel {
    /// Perfect theoretical driver.
    pub fn perfect() -> Self {
        Self {
            braking_point_offset: 0.0,
            throttle_application: 1.0,
            cornering_commitment: 1.0,
            error_sigma: 0.0,
        }
    }
    /// Typical experienced racing driver.
    pub fn experienced() -> Self {
        Self {
            braking_point_offset: -5.0,
            throttle_application: 0.92,
            cornering_commitment: 0.96,
            error_sigma: 0.05,
        }
    }
    /// Grip multiplier applied to the vehicle performance map.
    pub fn effective_grip_fraction(&self) -> f64 {
        self.cornering_commitment
    }
}
