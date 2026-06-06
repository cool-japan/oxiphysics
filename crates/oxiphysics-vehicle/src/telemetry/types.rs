//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Records a bounded ring of [`TelemetryFrame`] samples and provides
/// summary statistics and CSV export.
pub struct TelemetryRecorder {
    pub(super) frames: Vec<TelemetryFrame>,
    pub(super) max_frames: usize,
}
impl TelemetryRecorder {
    /// Create a new recorder that keeps at most `max_frames` samples.
    pub fn new(max_frames: usize) -> Self {
        Self {
            frames: Vec::with_capacity(max_frames.min(4096)),
            max_frames,
        }
    }
    /// Append a frame, discarding the oldest if the buffer is full.
    pub fn record(&mut self, frame: TelemetryFrame) {
        if self.frames.len() >= self.max_frames {
            self.frames.remove(0);
        }
        self.frames.push(frame);
    }
    /// Peak absolute lateral G across all recorded frames.
    pub fn max_lateral_g(&self) -> f64 {
        self.frames
            .iter()
            .map(|f| f.lateral_g.abs())
            .fold(0.0_f64, f64::max)
    }
    /// Peak absolute longitudinal G across all recorded frames.
    pub fn max_longitudinal_g(&self) -> f64 {
        self.frames
            .iter()
            .map(|f| f.longitudinal_g.abs())
            .fold(0.0_f64, f64::max)
    }
    /// Mean speed (magnitude of velocity vector) over all frames (m/s).
    pub fn average_speed(&self) -> f64 {
        if self.frames.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.frames.iter().map(|f| vec3_len(f.velocity)).sum();
        sum / self.frames.len() as f64
    }
    /// Total distance travelled, estimated by integrating speed over time
    /// using the trapezoidal rule.
    pub fn distance_traveled(&self) -> f64 {
        if self.frames.len() < 2 {
            return 0.0;
        }
        let mut dist = 0.0_f64;
        for w in self.frames.windows(2) {
            let dt = w[1].time - w[0].time;
            let v0 = vec3_len(w[0].velocity);
            let v1 = vec3_len(w[1].velocity);
            dist += 0.5 * (v0 + v1) * dt;
        }
        dist
    }
    /// Export all frames as a CSV string (header + one row per frame).
    pub fn export_csv(&self) -> String {
        let mut out = String::from(
            "time,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,\
             acc_x,acc_y,acc_z,steering_angle,throttle,brake,gear,rpm,\
             wheel_fl,wheel_fr,wheel_rl,wheel_rr,lateral_g,longitudinal_g\n",
        );
        for f in &self.frames {
            out.push_str(&format!(
                "{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},\
                 {:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{:.2},\
                 {:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                f.time,
                f.position[0],
                f.position[1],
                f.position[2],
                f.velocity[0],
                f.velocity[1],
                f.velocity[2],
                f.acceleration[0],
                f.acceleration[1],
                f.acceleration[2],
                f.steering_angle,
                f.throttle,
                f.brake,
                f.gear,
                f.rpm,
                f.wheel_speeds[0],
                f.wheel_speeds[1],
                f.wheel_speeds[2],
                f.wheel_speeds[3],
                f.lateral_g,
                f.longitudinal_g,
            ));
        }
        out
    }
}
/// Records per-tyre slip data over time for post-session analysis.
///
/// Stores both slip angles (lateral) and longitudinal slip ratios for all
/// four corners in a rolling window.
pub struct TireSlipRecorder {
    /// Timestamps (s).
    pub times: Vec<f64>,
    /// Lateral slip angles `[FL, FR, RL, RR]` (rad) at each time step.
    pub slip_angles: Vec<[f64; 4]>,
    /// Longitudinal slip ratios `[FL, FR, RL, RR]` (dimensionless, −1..1).
    pub long_slip: Vec<[f64; 4]>,
    pub(super) max_samples: usize,
}
impl TireSlipRecorder {
    /// Create a new recorder with a rolling window.
    pub fn new(max_samples: usize) -> Self {
        Self {
            times: Vec::new(),
            slip_angles: Vec::new(),
            long_slip: Vec::new(),
            max_samples: max_samples.max(2),
        }
    }
    /// Record one time step.
    ///
    /// * `t`            – simulation time (s)
    /// * `slip_ang`     – lateral slip angles `[FL, FR, RL, RR]` (rad)
    /// * `long_slip_r`  – longitudinal slip ratios `[FL, FR, RL, RR]`
    pub fn record(&mut self, t: f64, slip_ang: [f64; 4], long_slip_r: [f64; 4]) {
        if self.times.len() >= self.max_samples {
            self.times.remove(0);
            self.slip_angles.remove(0);
            self.long_slip.remove(0);
        }
        self.times.push(t);
        self.slip_angles.push(slip_ang);
        self.long_slip.push(long_slip_r);
    }
    /// Peak lateral slip angle magnitude across all tyres and all frames.
    pub fn peak_slip_angle(&self) -> f64 {
        self.slip_angles
            .iter()
            .flat_map(|row| row.iter())
            .map(|a| a.abs())
            .fold(0.0_f64, f64::max)
    }
    /// Peak longitudinal slip magnitude across all tyres and all frames.
    pub fn peak_long_slip(&self) -> f64 {
        self.long_slip
            .iter()
            .flat_map(|row| row.iter())
            .map(|s| s.abs())
            .fold(0.0_f64, f64::max)
    }
    /// Mean slip angle magnitude for a single tyre index (0..4) over the window.
    pub fn mean_slip_angle(&self, wheel: usize) -> f64 {
        if self.slip_angles.is_empty() {
            return 0.0;
        }
        let wheel = wheel.min(3);
        let sum: f64 = self.slip_angles.iter().map(|row| row[wheel].abs()).sum();
        sum / self.slip_angles.len() as f64
    }
    /// Fraction of time that tyre `wheel` is operating above `threshold` slip angle (rad).
    pub fn time_above_slip_threshold(&self, wheel: usize, threshold: f64) -> f64 {
        let n = self.times.len();
        if n < 2 {
            return 0.0;
        }
        let wheel = wheel.min(3);
        let mut total = 0.0_f64;
        for i in 1..n {
            if self.slip_angles[i][wheel].abs() >= threshold {
                total += self.times[i] - self.times[i - 1];
            }
        }
        total
    }
}
/// Analyzes driver input quality: smoothness, overlap, and peak usage.
///
/// Records throttle and brake inputs over a session and computes metrics
/// that are used in driver coaching systems.
pub struct DriverInputAnalyzer {
    pub(super) throttle_history: Vec<f64>,
    pub(super) brake_history: Vec<f64>,
    pub(super) steering_history: Vec<f64>,
    pub(super) times: Vec<f64>,
    pub(super) max_samples: usize,
}
impl DriverInputAnalyzer {
    /// Create a new input analyzer with a rolling window.
    pub fn new(max_samples: usize) -> Self {
        Self {
            throttle_history: Vec::new(),
            brake_history: Vec::new(),
            steering_history: Vec::new(),
            times: Vec::new(),
            max_samples: max_samples.max(2),
        }
    }
    /// Record one step of driver inputs.
    pub fn record(&mut self, t: f64, throttle: f64, brake: f64, steering: f64) {
        if self.times.len() >= self.max_samples {
            self.throttle_history.remove(0);
            self.brake_history.remove(0);
            self.steering_history.remove(0);
            self.times.remove(0);
        }
        self.throttle_history.push(throttle.clamp(0.0, 1.0));
        self.brake_history.push(brake.clamp(0.0, 1.0));
        self.steering_history.push(steering);
        self.times.push(t);
    }
    /// Throttle-brake overlap duration (s): time when both are non-zero.
    pub fn overlap_duration(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let mut total = 0.0_f64;
        for i in 1..self.times.len() {
            if self.throttle_history[i] > 0.01 && self.brake_history[i] > 0.01 {
                total += self.times[i] - self.times[i - 1];
            }
        }
        total
    }
    /// Throttle smoothness: mean absolute rate of change (per second).
    /// Lower is smoother.
    pub fn throttle_smoothness(&self) -> f64 {
        self.mean_rate_of_change(&self.throttle_history)
    }
    /// Steering smoothness: mean absolute rate of change (rad/s).
    pub fn steering_smoothness(&self) -> f64 {
        self.mean_rate_of_change(&self.steering_history)
    }
    /// Peak throttle usage (0–1) in the window.
    pub fn peak_throttle(&self) -> f64 {
        self.throttle_history
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }
    /// Peak brake usage (0–1) in the window.
    pub fn peak_brake(&self) -> f64 {
        self.brake_history.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Fraction of time at full throttle (> 0.95).
    pub fn full_throttle_fraction(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let total_t = self.times.last().expect("collection should not be empty")
            - self.times.first().expect("collection should not be empty");
        if total_t < 1e-12 {
            return 0.0;
        }
        let mut ft = 0.0_f64;
        for i in 1..self.times.len() {
            if self.throttle_history[i] > 0.95 {
                ft += self.times[i] - self.times[i - 1];
            }
        }
        ft / total_t
    }
    fn mean_rate_of_change(&self, values: &[f64]) -> f64 {
        let n = values.len();
        if n < 2 || self.times.len() < 2 {
            return 0.0;
        }
        let mut sum = 0.0_f64;
        for i in 1..n {
            let dt = self.times[i] - self.times[i - 1];
            if dt > 1e-12 {
                sum += (values[i] - values[i - 1]).abs() / dt;
            }
        }
        sum / (n - 1) as f64
    }
}
/// Detects pitlane entry and exit events from speed data and a speed-limit
/// threshold.
///
/// The pitlane is assumed to be entered when speed drops below `pit_speed_limit`
/// and exited when it rises above it again.
pub struct PitlaneDetector {
    /// Pitlane speed limit in m/s.
    pub pit_speed_limit: f64,
    /// Whether the vehicle is currently in the pitlane.
    pub(super) in_pit: bool,
    /// Time of the last pitlane entry (s).
    pub(super) entry_time: Option<f64>,
    /// Completed pit stops as (entry_time, exit_time).
    pub pit_stops: Vec<(f64, f64)>,
}
impl PitlaneDetector {
    /// Create a new pitlane detector.
    ///
    /// `pit_speed_limit` – typical: 60 km/h = 16.67 m/s for Formula cars
    pub fn new(pit_speed_limit: f64) -> Self {
        Self {
            pit_speed_limit,
            in_pit: false,
            entry_time: None,
            pit_stops: Vec::new(),
        }
    }
    /// Feed a telemetry sample.
    ///
    /// Returns `Some(PitEvent::Entry)` or `Some(PitEvent::Exit)` if an event
    /// occurred this step, `None` otherwise.
    pub fn update(&mut self, t: f64, speed: f64) -> Option<PitEvent> {
        if !self.in_pit && speed < self.pit_speed_limit {
            self.in_pit = true;
            self.entry_time = Some(t);
            return Some(PitEvent::Entry);
        }
        if self.in_pit && speed >= self.pit_speed_limit {
            self.in_pit = false;
            if let Some(entry) = self.entry_time.take() {
                self.pit_stops.push((entry, t));
            }
            return Some(PitEvent::Exit);
        }
        None
    }
    /// Number of completed pit stops.
    pub fn num_pit_stops(&self) -> usize {
        self.pit_stops.len()
    }
    /// Mean pit stop duration in seconds.
    pub fn mean_stop_duration(&self) -> f64 {
        if self.pit_stops.is_empty() {
            return 0.0;
        }
        let total: f64 = self.pit_stops.iter().map(|(e, x)| x - e).sum();
        total / self.pit_stops.len() as f64
    }
    /// Total time spent in the pitlane.
    pub fn total_pit_time(&self) -> f64 {
        self.pit_stops.iter().map(|(e, x)| x - e).sum()
    }
}
/// Pitlane event type emitted by [`PitlaneDetector`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PitEvent {
    /// Vehicle entered the pitlane (speed dropped below limit).
    Entry,
    /// Vehicle exited the pitlane (speed rose above limit).
    Exit,
}
/// Summary statistics computed from a full telemetry session.
#[derive(Debug, Clone)]
pub struct TelemetryStatistics {
    /// Total distance covered (m).
    pub total_distance_m: f64,
    /// Maximum speed reached (m/s).
    pub max_speed_ms: f64,
    /// Mean speed (m/s).
    pub mean_speed_ms: f64,
    /// Peak lateral G (g).
    pub peak_lateral_g: f64,
    /// Peak longitudinal G (g).
    pub peak_longitudinal_g: f64,
    /// Maximum RPM seen.
    pub max_rpm: f64,
    /// Duration of the session (s).
    pub duration_s: f64,
}
impl TelemetryStatistics {
    /// Compute statistics from a `TelemetryRecorder`.
    pub fn from_recorder(rec: &TelemetryRecorder) -> Self {
        if rec.frames.is_empty() {
            return Self {
                total_distance_m: 0.0,
                max_speed_ms: 0.0,
                mean_speed_ms: 0.0,
                peak_lateral_g: 0.0,
                peak_longitudinal_g: 0.0,
                max_rpm: 0.0,
                duration_s: 0.0,
            };
        }
        let max_speed = rec
            .frames
            .iter()
            .map(|f| vec3_len(f.velocity))
            .fold(0.0_f64, f64::max);
        let max_rpm = rec.frames.iter().map(|f| f.rpm).fold(0.0_f64, f64::max);
        let duration = rec
            .frames
            .last()
            .expect("collection should not be empty")
            .time
            - rec
                .frames
                .first()
                .expect("collection should not be empty")
                .time;
        Self {
            total_distance_m: rec.distance_traveled(),
            max_speed_ms: max_speed,
            mean_speed_ms: rec.average_speed(),
            peak_lateral_g: rec.max_lateral_g(),
            peak_longitudinal_g: rec.max_longitudinal_g(),
            max_rpm,
            duration_s: duration,
        }
    }
    /// Produce a human-readable one-line summary.
    pub fn summary_line(&self) -> String {
        format!(
            "dur={:.1}s dist={:.0}m vmax={:.1}km/h vmean={:.1}km/h latG={:.2}g longG={:.2}g rpm={:.0}",
            self.duration_s,
            self.total_distance_m,
            self.max_speed_ms * 3.6,
            self.mean_speed_ms * 3.6,
            self.peak_lateral_g,
            self.peak_longitudinal_g,
            self.max_rpm,
        )
    }
}
/// Describes a track by its total length and a sequence of 3-D checkpoints.
pub struct TrackPosition {
    /// Total track length in metres.
    pub track_length: f64,
    /// Ordered waypoints that approximate the track centre-line.
    pub checkpoints: Vec<[f64; 3]>,
}
impl TrackPosition {
    /// Return the arc-length distance along the track to the closest point on
    /// any segment to `pos`.
    pub fn nearest_track_distance(&self, pos: [f64; 3]) -> f64 {
        if self.checkpoints.len() < 2 {
            return 0.0;
        }
        let mut best_dist_sq = f64::INFINITY;
        let mut best_arc = 0.0_f64;
        let mut arc = 0.0_f64;
        let n = self.checkpoints.len();
        for i in 0..n {
            let a = self.checkpoints[i];
            let b = self.checkpoints[(i + 1) % n];
            let seg_len = vec3_len(vec3_sub(b, a));
            let t = if seg_len > 1e-12 {
                vec3_dot(vec3_sub(pos, a), vec3_sub(b, a)) / (seg_len * seg_len)
            } else {
                0.0
            };
            let t = t.clamp(0.0, 1.0);
            let proj = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            let d_sq = vec3_dist_sq(pos, proj);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_arc = arc + t * seg_len;
            }
            arc += seg_len;
        }
        best_arc
    }
    /// Return `true` when the vehicle crosses the start/finish circle between
    /// two consecutive positions.  The circle is centred on `start_finish`
    /// with the given `radius`.
    pub fn detect_lap_complete(
        &self,
        prev_pos: [f64; 3],
        curr_pos: [f64; 3],
        start_finish: [f64; 3],
        radius: f64,
    ) -> bool {
        let prev_inside = vec3_dist_sq(prev_pos, start_finish) <= radius * radius;
        let curr_inside = vec3_dist_sq(curr_pos, start_finish) <= radius * radius;
        !prev_inside && curr_inside
    }
}
/// Tracks sector and lap times for a circuit with a fixed number of sectors.
pub struct LapTimer {
    pub(super) start_time: Option<f64>,
    pub(super) lap_times: Vec<f64>,
    pub(super) sector_times: Vec<Vec<f64>>,
    pub(super) _num_sectors: usize,
    pub(super) current_sector: usize,
    pub(super) sector_start_time: Option<f64>,
}
impl LapTimer {
    /// Create a timer for a track with `num_sectors` sectors.
    pub fn new(num_sectors: usize) -> Self {
        Self {
            start_time: None,
            lap_times: Vec::new(),
            sector_times: Vec::new(),
            _num_sectors: num_sectors,
            current_sector: 0,
            sector_start_time: None,
        }
    }
    /// Start a new lap at absolute time `time`.
    pub fn start_lap(&mut self, time: f64) {
        self.start_time = Some(time);
        self.sector_start_time = Some(time);
        self.current_sector = 0;
        self.sector_times.push(Vec::new());
    }
    /// Mark the end of the current sector at absolute time `time`.
    /// Automatically advances the sector counter.
    pub fn complete_sector(&mut self, time: f64) {
        if let Some(sector_start) = self.sector_start_time {
            let elapsed = time - sector_start;
            if let Some(row) = self.sector_times.last_mut() {
                row.push(elapsed);
            }
            self.sector_start_time = Some(time);
            self.current_sector += 1;
        }
    }
    /// Close the current lap at absolute time `time`.
    /// Returns the lap time (time since `start_lap`).
    pub fn complete_lap(&mut self, time: f64) -> f64 {
        let lap_time = match self.start_time {
            Some(t) => time - t,
            None => 0.0,
        };
        self.lap_times.push(lap_time);
        self.start_time = None;
        self.sector_start_time = None;
        self.current_sector = 0;
        lap_time
    }
    /// The fastest recorded lap time, or `None` if no laps completed.
    pub fn best_lap(&self) -> Option<f64> {
        self.lap_times.iter().cloned().reduce(f64::min)
    }
    /// The most recently completed lap time, or `None` if no laps completed.
    pub fn last_lap(&self) -> Option<f64> {
        self.lap_times.last().cloned()
    }
}
/// Tracks the temperature trend of all four tires over a rolling window.
///
/// Records up to `max_samples` temperature values and exposes statistics such
/// as the current mean, rate of change, and whether the tire is warming or
/// cooling.
pub struct TireTempTrend {
    /// Historical temperatures for each wheel `[FL, FR, RL, RR]` (°C).
    pub(super) history: [Vec<f64>; 4],
    /// Corresponding timestamps in seconds.
    pub(super) times: Vec<f64>,
    /// Maximum number of samples retained.
    pub(super) max_samples: usize,
}
impl TireTempTrend {
    /// Create a new trend tracker with the given rolling-window size.
    pub fn new(max_samples: usize) -> Self {
        Self {
            history: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            times: Vec::new(),
            max_samples: max_samples.max(2),
        }
    }
    /// Record a temperature snapshot at simulation time `t`.
    ///
    /// `temps` = `[FL, FR, RL, RR]` in °C.
    pub fn record(&mut self, t: f64, temps: [f64; 4]) {
        if self.times.len() >= self.max_samples {
            self.times.remove(0);
            for w in self.history.iter_mut() {
                w.remove(0);
            }
        }
        self.times.push(t);
        for (i, w) in self.history.iter_mut().enumerate() {
            w.push(temps[i]);
        }
    }
    /// Current (most recent) temperature for wheel `wheel` in °C, or `None`
    /// if no samples have been recorded.
    pub fn current(&self, wheel: usize) -> Option<f64> {
        self.history.get(wheel)?.last().cloned()
    }
    /// Mean temperature over the entire window for wheel `wheel`.
    pub fn mean(&self, wheel: usize) -> f64 {
        let h = &self.history[wheel.min(3)];
        if h.is_empty() {
            return 0.0;
        }
        h.iter().sum::<f64>() / h.len() as f64
    }
    /// Estimated rate of temperature change for wheel `wheel` in °C/s.
    ///
    /// Computed as `(latest − oldest) / (t_latest − t_oldest)`.
    /// Returns `0.0` when fewer than two samples are available.
    pub fn rate_of_change(&self, wheel: usize) -> f64 {
        let h = &self.history[wheel.min(3)];
        if h.len() < 2 || self.times.len() < 2 {
            return 0.0;
        }
        let dt = self.times.last().expect("collection should not be empty")
            - self.times.first().expect("collection should not be empty");
        if dt.abs() < 1e-12 {
            return 0.0;
        }
        let dtemp = h.last().expect("collection should not be empty")
            - h.first().expect("collection should not be empty");
        dtemp / dt
    }
    /// Returns `true` when the tire is warming (positive rate of change).
    pub fn is_warming(&self, wheel: usize) -> bool {
        self.rate_of_change(wheel) > 0.0
    }
    /// Maximum temperature reached in the window for wheel `wheel`.
    pub fn peak(&self, wheel: usize) -> f64 {
        self.history[wheel.min(3)]
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum temperature reached in the window for wheel `wheel`.
    pub fn trough(&self, wheel: usize) -> f64 {
        self.history[wheel.min(3)]
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
}
/// Analyzes brake pressure traces to identify braking events and their quality.
pub struct BrakeTraceAnalyzer {
    /// Brake pressure samples (0–1).
    pub(super) pressures: Vec<f64>,
    /// Corresponding timestamps (s).
    pub(super) times: Vec<f64>,
    /// Detection threshold above which a braking event starts.
    pub threshold: f64,
}
impl BrakeTraceAnalyzer {
    /// Create a new analyzer.
    pub fn new(threshold: f64) -> Self {
        Self {
            pressures: Vec::new(),
            times: Vec::new(),
            threshold: threshold.clamp(0.0, 1.0),
        }
    }
    /// Append a brake sample.
    pub fn record(&mut self, t: f64, pressure: f64) {
        self.pressures.push(pressure.clamp(0.0, 1.0));
        self.times.push(t);
    }
    /// Detect braking events as (start_time, end_time, peak_pressure) tuples.
    pub fn detect_events(&self) -> Vec<(f64, f64, f64)> {
        let n = self.pressures.len();
        if n < 2 {
            return Vec::new();
        }
        let mut events = Vec::new();
        let mut in_event = false;
        let mut event_start = 0.0_f64;
        let mut peak = 0.0_f64;
        for i in 0..n {
            let p = self.pressures[i];
            if !in_event && p >= self.threshold {
                in_event = true;
                event_start = self.times[i];
                peak = p;
            } else if in_event {
                if p > peak {
                    peak = p;
                }
                if p < self.threshold || i == n - 1 {
                    events.push((event_start, self.times[i], peak));
                    in_event = false;
                    peak = 0.0;
                }
            }
        }
        events
    }
    /// Number of distinct braking events in the trace.
    pub fn event_count(&self) -> usize {
        self.detect_events().len()
    }
    /// Mean peak brake pressure across all detected events.
    pub fn mean_peak_pressure(&self) -> f64 {
        let events = self.detect_events();
        if events.is_empty() {
            return 0.0;
        }
        events.iter().map(|e| e.2).sum::<f64>() / events.len() as f64
    }
    /// Total time spent braking (above threshold), in seconds.
    pub fn total_braking_time(&self) -> f64 {
        self.detect_events().iter().map(|e| e.1 - e.0).sum()
    }
}
/// Compares two telemetry recordings at equal track distances to find the
/// time gained or lost and the speed difference at each point.
pub struct LapComparison {
    /// Reference (faster / previous) lap frames.
    pub reference: Vec<TelemetryFrame>,
    /// Target (current) lap frames.
    pub target: Vec<TelemetryFrame>,
}
impl LapComparison {
    /// Create a new lap comparison.
    pub fn new(reference: Vec<TelemetryFrame>, target: Vec<TelemetryFrame>) -> Self {
        Self { reference, target }
    }
    /// Time delta at each frame index: `target.time[i] − reference.time[i]`.
    ///
    /// Positive = target is slower at that point.
    pub fn time_deltas(&self) -> Vec<f64> {
        let n = self.reference.len().min(self.target.len());
        (0..n)
            .map(|i| self.target[i].time - self.reference[i].time)
            .collect()
    }
    /// Speed delta at each frame: `|v_target[i]| − |v_ref[i]|` (m/s).
    pub fn speed_deltas(&self) -> Vec<f64> {
        let n = self.reference.len().min(self.target.len());
        (0..n)
            .map(|i| {
                let v_tgt = vec3_len(self.target[i].velocity);
                let v_ref = vec3_len(self.reference[i].velocity);
                v_tgt - v_ref
            })
            .collect()
    }
    /// Throttle delta at each frame: `target.throttle[i] − reference.throttle[i]`.
    pub fn throttle_deltas(&self) -> Vec<f64> {
        let n = self.reference.len().min(self.target.len());
        (0..n)
            .map(|i| self.target[i].throttle - self.reference[i].throttle)
            .collect()
    }
    /// Total lap time difference: `target_lap_time − reference_lap_time`.
    ///
    /// Requires at least one frame in each recording.
    pub fn total_time_delta(&self) -> f64 {
        let ref_dur = match (self.reference.first(), self.reference.last()) {
            (Some(f), Some(l)) => l.time - f.time,
            _ => 0.0,
        };
        let tgt_dur = match (self.target.first(), self.target.last()) {
            (Some(f), Some(l)) => l.time - f.time,
            _ => 0.0,
        };
        tgt_dur - ref_dur
    }
    /// Maximum speed advantage of the target over the reference (positive if target faster).
    pub fn max_speed_advantage(&self) -> f64 {
        self.speed_deltas()
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Number of frames where the target is faster than the reference.
    pub fn frames_faster_count(&self) -> usize {
        self.time_deltas().iter().filter(|&&d| d < 0.0).count()
    }
}
/// A single telemetry sample captured at one simulation step.
#[derive(Debug, Clone)]
pub struct TelemetryFrame {
    /// Simulation time in seconds.
    pub time: f64,
    /// World-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// World-space velocity `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
    /// World-space acceleration `[ax, ay, az]` (m/s²).
    pub acceleration: [f64; 3],
    /// Front-wheel steering angle in radians.
    pub steering_angle: f64,
    /// Throttle position in `[0, 1]`.
    pub throttle: f64,
    /// Brake pressure in `[0, 1]`.
    pub brake: f64,
    /// Current gear (0 = neutral, negative = reverse).
    pub gear: i32,
    /// Engine RPM.
    pub rpm: f64,
    /// Wheel speeds `[FL, FR, RL, RR]` (rad/s).
    pub wheel_speeds: [f64; 4],
    /// Lateral G-force.
    pub lateral_g: f64,
    /// Longitudinal G-force.
    pub longitudinal_g: f64,
}
/// Aggregates multiple lap times, sector times, and pit stop data into a
/// concise session summary.
pub struct SessionSummary {
    /// All recorded lap times (s).
    pub lap_times: Vec<f64>,
    /// All recorded sector times, grouped by lap.
    pub sector_times: Vec<Vec<f64>>,
    /// Total pit time (s).
    pub total_pit_time: f64,
    /// Number of pit stops.
    pub num_pit_stops: usize,
}
impl SessionSummary {
    /// Create an empty session summary.
    pub fn new() -> Self {
        Self {
            lap_times: Vec::new(),
            sector_times: Vec::new(),
            total_pit_time: 0.0,
            num_pit_stops: 0,
        }
    }
    /// Record a completed lap.
    pub fn record_lap(&mut self, lap_time: f64, sectors: Vec<f64>) {
        self.lap_times.push(lap_time);
        self.sector_times.push(sectors);
    }
    /// Set pit stop data.
    pub fn set_pit_data(&mut self, total_pit_time: f64, num_pit_stops: usize) {
        self.total_pit_time = total_pit_time;
        self.num_pit_stops = num_pit_stops;
    }
    /// Best lap time (s), or `None` if no laps completed.
    pub fn best_lap(&self) -> Option<f64> {
        self.lap_times.iter().cloned().reduce(f64::min)
    }
    /// Mean lap time (s).
    pub fn mean_lap(&self) -> f64 {
        if self.lap_times.is_empty() {
            return 0.0;
        }
        self.lap_times.iter().sum::<f64>() / self.lap_times.len() as f64
    }
    /// Standard deviation of lap times (s).
    pub fn lap_time_std(&self) -> f64 {
        let n = self.lap_times.len();
        if n < 2 {
            return 0.0;
        }
        let mean = self.mean_lap();
        let var = self
            .lap_times
            .iter()
            .map(|&t| (t - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        var.sqrt()
    }
    /// Fastest sector times across all laps (one per sector position).
    ///
    /// Returns an empty vector if no sectors were recorded.
    pub fn theoretical_best_lap(&self) -> Vec<f64> {
        if self.sector_times.is_empty() {
            return Vec::new();
        }
        let n_sectors = self.sector_times.iter().map(|s| s.len()).max().unwrap_or(0);
        (0..n_sectors)
            .map(|i| {
                self.sector_times
                    .iter()
                    .filter_map(|row| row.get(i))
                    .cloned()
                    .reduce(f64::min)
                    .unwrap_or(0.0)
            })
            .collect()
    }
    /// Sum of the theoretical best sectors (s).
    pub fn theoretical_best_total(&self) -> f64 {
        self.theoretical_best_lap().iter().sum()
    }
    /// Produce a formatted multi-line session report.
    pub fn report(&self) -> String {
        let best = self
            .best_lap()
            .map(|t| format!("{t:.3}s"))
            .unwrap_or_else(|| "N/A".to_string());
        let mean = self.mean_lap();
        let std = self.lap_time_std();
        let tbl = self.theoretical_best_total();
        format!(
            "Laps: {} | Best: {} | Mean: {:.3}s | Std: {:.3}s | Theo-best: {:.3}s | Pit stops: {} ({:.1}s)",
            self.lap_times.len(),
            best,
            mean,
            std,
            tbl,
            self.num_pit_stops,
            self.total_pit_time,
        )
    }
}
/// Compares two laps sector by sector to find time gains/losses.
pub struct SectorComparison {
    /// Reference (e.g. best) lap sector times.
    pub reference_sectors: Vec<f64>,
    /// Current (or comparison) lap sector times.
    pub current_sectors: Vec<f64>,
}
impl SectorComparison {
    /// Create a new sector comparison.
    pub fn new(reference: Vec<f64>, current: Vec<f64>) -> Self {
        Self {
            reference_sectors: reference,
            current_sectors: current,
        }
    }
    /// Delta time per sector (current − reference, seconds).
    /// Negative = improvement; positive = time lost.
    pub fn sector_deltas(&self) -> Vec<f64> {
        let n = self.reference_sectors.len().min(self.current_sectors.len());
        (0..n)
            .map(|i| self.current_sectors[i] - self.reference_sectors[i])
            .collect()
    }
    /// Cumulative delta at each sector boundary.
    pub fn cumulative_deltas(&self) -> Vec<f64> {
        let deltas = self.sector_deltas();
        let mut cum = 0.0_f64;
        deltas
            .iter()
            .map(|&d| {
                cum += d;
                cum
            })
            .collect()
    }
    /// Total lap time delta (positive = slower than reference).
    pub fn total_delta(&self) -> f64 {
        self.sector_deltas().iter().sum()
    }
    /// Index of the sector with the largest time loss.
    pub fn worst_sector(&self) -> Option<usize> {
        let deltas = self.sector_deltas();
        deltas
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
    /// Index of the sector with the largest time gain (most improvement).
    pub fn best_sector(&self) -> Option<usize> {
        let deltas = self.sector_deltas();
        deltas
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
}
/// Identifies cornering events from lateral G data and computes corner statistics.
pub struct CornerAnalyzer {
    pub(super) lateral_g_history: Vec<f64>,
    pub(super) times: Vec<f64>,
    /// Lateral-G threshold for corner detection.
    pub g_threshold: f64,
}
impl CornerAnalyzer {
    /// Create a new corner analyzer.
    pub fn new(g_threshold: f64) -> Self {
        Self {
            lateral_g_history: Vec::new(),
            times: Vec::new(),
            g_threshold: g_threshold.abs(),
        }
    }
    /// Record a lateral G sample.
    pub fn record(&mut self, t: f64, lateral_g: f64) {
        self.lateral_g_history.push(lateral_g);
        self.times.push(t);
    }
    /// Detect corners as (start_time, end_time, peak_g) tuples.
    pub fn detect_corners(&self) -> Vec<(f64, f64, f64)> {
        let n = self.lateral_g_history.len();
        if n < 2 {
            return Vec::new();
        }
        let mut corners = Vec::new();
        let mut in_corner = false;
        let mut start = 0.0_f64;
        let mut peak_g = 0.0_f64;
        for i in 0..n {
            let g = self.lateral_g_history[i].abs();
            if !in_corner && g >= self.g_threshold {
                in_corner = true;
                start = self.times[i];
                peak_g = g;
            } else if in_corner {
                if g > peak_g {
                    peak_g = g;
                }
                if g < self.g_threshold || i == n - 1 {
                    corners.push((start, self.times[i], peak_g));
                    in_corner = false;
                    peak_g = 0.0;
                }
            }
        }
        corners
    }
    /// Number of detected corners.
    pub fn corner_count(&self) -> usize {
        self.detect_corners().len()
    }
    /// Mean peak lateral G across all corners.
    pub fn mean_corner_g(&self) -> f64 {
        let corners = self.detect_corners();
        if corners.is_empty() {
            return 0.0;
        }
        corners.iter().map(|c| c.2).sum::<f64>() / corners.len() as f64
    }
    /// Total time spent cornering (above threshold), in seconds.
    pub fn total_corner_time(&self) -> f64 {
        self.detect_corners().iter().map(|c| c.1 - c.0).sum()
    }
}
/// Logs Energy Recovery System (ERS) state: deployed energy, harvested energy,
/// battery state-of-charge, and MGU deployment power.
pub struct ErsEnergyLog {
    /// Timestamps (s).
    pub times: Vec<f64>,
    /// ERS deployment power (W) at each step (positive = deploy, negative = harvest).
    pub power: Vec<f64>,
    /// Battery state-of-charge fraction `[0, 1]` at each step.
    pub soc: Vec<f64>,
    /// Maximum battery capacity (J).
    pub capacity_j: f64,
    /// Current stored energy (J).
    pub(super) stored_j: f64,
}
impl ErsEnergyLog {
    /// Create a new ERS logger.
    ///
    /// * `capacity_j`  – battery capacity in Joules
    /// * `initial_soc` – initial state-of-charge (0–1)
    pub fn new(capacity_j: f64, initial_soc: f64) -> Self {
        let soc = initial_soc.clamp(0.0, 1.0);
        Self {
            times: Vec::new(),
            power: Vec::new(),
            soc: Vec::new(),
            capacity_j: capacity_j.max(1.0),
            stored_j: capacity_j * soc,
        }
    }
    /// Step the ERS logger by `dt` seconds with power `p_w` (W).
    ///
    /// Positive `p_w` = deploy (discharge battery); negative = harvest (charge).
    pub fn step(&mut self, t: f64, p_w: f64, dt: f64) {
        let delta_j = p_w * dt;
        self.stored_j = (self.stored_j - delta_j).clamp(0.0, self.capacity_j);
        let soc = self.stored_j / self.capacity_j;
        self.times.push(t);
        self.power.push(p_w);
        self.soc.push(soc);
    }
    /// Current state-of-charge (0–1).
    pub fn current_soc(&self) -> f64 {
        self.stored_j / self.capacity_j
    }
    /// Total energy deployed (J) over the log.
    pub fn total_deployed_j(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let mut total = 0.0f64;
        for i in 1..self.times.len() {
            let dt = self.times[i] - self.times[i - 1];
            let p = self.power[i];
            if p > 0.0 {
                total += p * dt;
            }
        }
        total
    }
    /// Total energy harvested (J) over the log.
    pub fn total_harvested_j(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let mut total = 0.0f64;
        for i in 1..self.times.len() {
            let dt = self.times[i] - self.times[i - 1];
            let p = self.power[i];
            if p < 0.0 {
                total += (-p) * dt;
            }
        }
        total
    }
    /// Peak deployment power (W).
    pub fn peak_deployment_power(&self) -> f64 {
        self.power.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Battery delta over the log (final_soc − initial_soc, in fraction).
    pub fn soc_delta(&self) -> f64 {
        match (self.soc.first(), self.soc.last()) {
            (Some(&first), Some(&last)) => last - first,
            _ => 0.0,
        }
    }
}
/// Replays previously recorded telemetry at a specified playback speed.
pub struct TelemetryReplay {
    /// The frames to replay.
    pub(super) frames: Vec<TelemetryFrame>,
    /// Current playback position (frame index).
    pub(super) cursor: usize,
    /// Playback speed multiplier (1.0 = real-time, 2.0 = 2× speed).
    pub playback_speed: f64,
    /// Simulation time of the last `advance` call.
    pub(super) last_time: f64,
}
impl TelemetryReplay {
    /// Create a new replay session from a set of frames.
    pub fn new(frames: Vec<TelemetryFrame>, playback_speed: f64) -> Self {
        Self {
            frames,
            cursor: 0,
            playback_speed: playback_speed.max(1e-6),
            last_time: 0.0,
        }
    }
    /// Advance replay by `real_dt` seconds of real time.
    ///
    /// Returns the current frame (or the last frame if the replay has ended).
    pub fn advance(&mut self, real_dt: f64) -> Option<&TelemetryFrame> {
        if self.frames.is_empty() {
            return None;
        }
        self.last_time += real_dt * self.playback_speed;
        while self.cursor + 1 < self.frames.len()
            && self.frames[self.cursor + 1].time <= self.last_time
        {
            self.cursor += 1;
        }
        self.frames.get(self.cursor)
    }
    /// Returns `true` when all frames have been replayed.
    pub fn is_finished(&self) -> bool {
        self.cursor + 1 >= self.frames.len()
    }
    /// Reset replay to the beginning.
    pub fn reset(&mut self) {
        self.cursor = 0;
        self.last_time = 0.0;
    }
    /// Total replay duration in recording time.
    pub fn duration(&self) -> f64 {
        match (self.frames.first(), self.frames.last()) {
            (Some(first), Some(last)) => last.time - first.time,
            _ => 0.0,
        }
    }
}
/// Simple fuel consumption model based on throttle position and engine speed.
///
/// The instantaneous consumption rate is computed as:
///
/// `rate = base_rate + throttle_coeff × throttle + rpm_coeff × rpm`
///
/// where `base_rate`, `throttle_coeff`, and `rpm_coeff` are empirical
/// calibration constants.  All rates are in kg/s.
pub struct FuelModel {
    /// Baseline fuel rate at zero throttle and idle RPM (kg/s).
    pub base_rate: f64,
    /// Additional fuel rate per unit throttle (kg/s per unit throttle).
    pub throttle_coeff: f64,
    /// Additional fuel rate per 1000 RPM (kg/s per 1000 rpm).
    pub rpm_coeff: f64,
    /// Current fuel mass remaining in kg.
    pub fuel_kg: f64,
    /// Maximum fuel capacity in kg.
    pub capacity_kg: f64,
}
impl FuelModel {
    /// Create a model suitable for a Formula-style car.
    ///
    /// Typical values: base 0.002 kg/s, 0.04 kg/s per unit throttle,
    /// 0.001 kg/s per 1000 RPM, 100 kg capacity.
    pub fn formula_car() -> Self {
        Self {
            base_rate: 0.002,
            throttle_coeff: 0.04,
            rpm_coeff: 0.001,
            fuel_kg: 100.0,
            capacity_kg: 100.0,
        }
    }
    /// Instantaneous fuel consumption rate in kg/s.
    ///
    /// Returns 0 when fuel is empty.
    pub fn consumption_rate(&self, throttle: f64, rpm: f64) -> f64 {
        if self.fuel_kg <= 0.0 {
            return 0.0;
        }
        let throttle_c = throttle.clamp(0.0, 1.0);
        let rate =
            self.base_rate + self.throttle_coeff * throttle_c + self.rpm_coeff * (rpm / 1000.0);
        rate.max(0.0)
    }
    /// Advance the fuel model by `dt` seconds.
    ///
    /// Deducts the consumed fuel and returns the actual consumption in kg.
    pub fn step(&mut self, throttle: f64, rpm: f64, dt: f64) -> f64 {
        let rate = self.consumption_rate(throttle, rpm);
        let consumed = (rate * dt).min(self.fuel_kg.max(0.0));
        self.fuel_kg = (self.fuel_kg - consumed).max(0.0);
        consumed
    }
    /// State of charge as a fraction in `[0, 1]`.
    pub fn fuel_fraction(&self) -> f64 {
        if self.capacity_kg <= 0.0 {
            return 0.0;
        }
        (self.fuel_kg / self.capacity_kg).clamp(0.0, 1.0)
    }
    /// Estimated remaining range in seconds at the current throttle and RPM.
    ///
    /// Returns `f64::INFINITY` when the consumption rate is effectively zero.
    pub fn remaining_time(&self, throttle: f64, rpm: f64) -> f64 {
        let rate = self.consumption_rate(throttle, rpm);
        if rate < 1e-12 {
            return f64::INFINITY;
        }
        self.fuel_kg / rate
    }
}
/// Tracks energy (fuel + electric) consumption with optional regeneration.
pub struct EnergyConsumption {
    /// Total energy consumed (J or Wh depending on unit choice).
    pub total_consumed: f64,
    /// Total energy recovered (e.g., regen braking).
    pub total_recovered: f64,
    /// Energy per unit distance (J/m) – updated during `record_step`.
    pub energy_per_km: f64,
    /// Internal accumulator for distance.
    pub(super) distance_acc: f64,
    /// Internal accumulator for consumed energy since last km.
    pub(super) energy_acc: f64,
}
impl EnergyConsumption {
    /// Create a new energy consumption tracker.
    pub fn new() -> Self {
        Self {
            total_consumed: 0.0,
            total_recovered: 0.0,
            energy_per_km: 0.0,
            distance_acc: 0.0,
            energy_acc: 0.0,
        }
    }
    /// Record one simulation step.
    ///
    /// `consumed`: energy consumed this step (same unit as the tracker).
    /// `recovered`: energy recovered this step.
    /// `distance`: distance travelled this step (m).
    pub fn record_step(&mut self, consumed: f64, recovered: f64, distance: f64) {
        self.total_consumed += consumed.max(0.0);
        self.total_recovered += recovered.max(0.0);
        self.distance_acc += distance.max(0.0);
        self.energy_acc += consumed.max(0.0);
        if self.distance_acc >= 1000.0 {
            self.energy_per_km = self.energy_acc / (self.distance_acc / 1000.0);
            self.distance_acc = 0.0;
            self.energy_acc = 0.0;
        }
    }
    /// Net energy consumed (consumed − recovered).
    pub fn net_energy(&self) -> f64 {
        self.total_consumed - self.total_recovered
    }
    /// Efficiency: recovered / consumed, in \[0, 1\].  Returns 0 if nothing consumed.
    pub fn recovery_efficiency(&self) -> f64 {
        if self.total_consumed < 1e-15 {
            return 0.0;
        }
        (self.total_recovered / self.total_consumed).clamp(0.0, 1.0)
    }
}
/// A speed-versus-distance trace along a lap.
///
/// Samples vehicle speed at discrete track positions.
pub struct SpeedTrace {
    /// Track distance stamps (m).
    pub distances: Vec<f64>,
    /// Vehicle speeds (m/s) at each distance.
    pub speeds: Vec<f64>,
}
impl SpeedTrace {
    /// Create a new empty trace.
    pub fn new() -> Self {
        Self {
            distances: Vec::new(),
            speeds: Vec::new(),
        }
    }
    /// Append a sample.
    pub fn record(&mut self, distance: f64, speed: f64) {
        self.distances.push(distance.max(0.0));
        self.speeds.push(speed.max(0.0));
    }
    /// Minimum speed (braking points) in the trace.
    pub fn min_speed(&self) -> f64 {
        self.speeds.iter().cloned().fold(f64::INFINITY, f64::min)
    }
    /// Maximum speed (straight-line) in the trace.
    pub fn max_speed(&self) -> f64 {
        self.speeds.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Distance of the speed trap (location of max speed).
    pub fn max_speed_distance(&self) -> Option<f64> {
        self.speeds
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| self.distances[i])
    }
    /// Speed at a given track distance (m/s), using linear interpolation.
    pub fn speed_at(&self, distance: f64) -> f64 {
        let n = self.distances.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 || distance <= self.distances[0] {
            return self.speeds[0];
        }
        if distance >= self.distances[n - 1] {
            return self.speeds[n - 1];
        }
        let idx = self
            .distances
            .partition_point(|&d| d <= distance)
            .saturating_sub(1);
        let idx = idx.min(n - 2);
        let d0 = self.distances[idx];
        let d1 = self.distances[idx + 1];
        let v0 = self.speeds[idx];
        let v1 = self.speeds[idx + 1];
        let dd = d1 - d0;
        if dd.abs() < 1e-12 {
            return v0;
        }
        v0 + (v1 - v0) * (distance - d0) / dd
    }
    /// Delta speed between two traces at the same distance positions.
    /// Returns a new trace of (distance, speed_delta) pairs.
    pub fn delta_versus(&self, other: &SpeedTrace) -> Vec<(f64, f64)> {
        self.distances
            .iter()
            .enumerate()
            .map(|(i, &d)| {
                let v_self = self.speeds[i];
                let v_other = other.speed_at(d);
                (d, v_self - v_other)
            })
            .collect()
    }
}
/// Records 3-D G-force vectors over time.
pub struct GForceRecorder {
    /// Timestamps (s).
    pub times: Vec<f64>,
    /// G-force vectors \[gx, gy, gz\] at each timestamp.
    pub g_vectors: Vec<[f64; 3]>,
    /// Maximum number of samples (ring buffer).
    pub(super) max_samples: usize,
}
impl GForceRecorder {
    /// Create a new recorder.
    pub fn new(max_samples: usize) -> Self {
        Self {
            times: Vec::new(),
            g_vectors: Vec::new(),
            max_samples: max_samples.max(1),
        }
    }
    /// Record a G-force vector at time `t`.
    pub fn record(&mut self, t: f64, g: [f64; 3]) {
        if self.times.len() >= self.max_samples {
            self.times.remove(0);
            self.g_vectors.remove(0);
        }
        self.times.push(t);
        self.g_vectors.push(g);
    }
    /// Peak resultant G (magnitude) across all recorded frames.
    pub fn peak_resultant_g(&self) -> f64 {
        self.g_vectors
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }
    /// Mean resultant G across all frames.
    pub fn mean_resultant_g(&self) -> f64 {
        if self.g_vectors.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .g_vectors
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
            .sum();
        sum / self.g_vectors.len() as f64
    }
    /// Time above threshold G magnitude (seconds).
    pub fn time_above_g(&self, threshold: f64) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let mut total = 0.0_f64;
        for i in 1..self.times.len() {
            let g = (self.g_vectors[i][0].powi(2)
                + self.g_vectors[i][1].powi(2)
                + self.g_vectors[i][2].powi(2))
            .sqrt();
            if g >= threshold {
                total += self.times[i] - self.times[i - 1];
            }
        }
        total
    }
}
/// A single named scalar channel recording values over time.
pub struct TelemetryChannel {
    /// Human-readable name.
    pub name: String,
    /// Sample timestamps in seconds.
    pub times: Vec<f64>,
    /// Sample values (arbitrary units).
    pub values: Vec<f64>,
}
impl TelemetryChannel {
    /// Create an empty channel with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            times: Vec::new(),
            values: Vec::new(),
        }
    }
    /// Append a sample.
    pub fn record(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }
    /// Linear interpolation: returns the interpolated value at `t`.
    /// Clamps to the first/last value outside the recorded range.
    pub fn interpolate_linear(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 || t <= self.times[0] {
            return self.values[0];
        }
        if t >= self.times[n - 1] {
            return self.values[n - 1];
        }
        let idx = self.times.partition_point(|&ts| ts <= t).saturating_sub(1);
        let idx = idx.min(n - 2);
        let t0 = self.times[idx];
        let t1 = self.times[idx + 1];
        let v0 = self.values[idx];
        let v1 = self.values[idx + 1];
        let dt = t1 - t0;
        if dt.abs() < 1e-15 {
            return v0;
        }
        v0 + (v1 - v0) * (t - t0) / dt
    }
    /// Cubic (Catmull-Rom) interpolation at time `t`.
    /// Falls back to linear at the boundaries (first and last interval).
    pub fn interpolate_cubic(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n < 2 {
            return self.interpolate_linear(t);
        }
        if t <= self.times[0] {
            return self.values[0];
        }
        if t >= self.times[n - 1] {
            return self.values[n - 1];
        }
        let idx = self.times.partition_point(|&ts| ts <= t).saturating_sub(1);
        let idx = idx.min(n - 2);
        let i0 = if idx > 0 { idx - 1 } else { 0 };
        let i1 = idx;
        let i2 = idx + 1;
        let i3 = (idx + 2).min(n - 1);
        let t0 = self.times[i1];
        let t1 = self.times[i2];
        let dt = t1 - t0;
        if dt.abs() < 1e-15 {
            return self.values[i1];
        }
        let alpha = (t - t0) / dt;
        let p0 = self.values[i0];
        let p1 = self.values[i1];
        let p2 = self.values[i2];
        let p3 = self.values[i3];
        let a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
        let b = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
        let c = -0.5 * p0 + 0.5 * p2;
        let d = p1;
        ((a * alpha + b) * alpha + c) * alpha + d
    }
    /// Number of samples.
    pub fn len(&self) -> usize {
        self.times.len()
    }
    /// Returns true if no samples have been recorded.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
}
