// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation logging, telemetry, and replay utilities.
//!
//! Provides structured logging of per-step simulation data (energy, timing,
//! body counts), filtering, compression, metrics, snapshotting, and a
//! ring-buffer streaming log.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// SimLogEntry
// ---------------------------------------------------------------------------

/// A single logged snapshot of simulation state at one timestep.
#[derive(Debug, Clone, PartialEq)]
pub struct SimLogEntry {
    /// Simulation timestep index (0-based).
    pub timestep: u64,
    /// Simulation time in seconds.
    pub sim_time: f64,
    /// Wall-clock time in seconds since start.
    pub wall_time: f64,
    /// Total mechanical energy at this step.
    pub energy: f64,
    /// Number of bodies/particles active at this step.
    pub n_bodies: usize,
    /// Arbitrary additional fields (name → value).
    pub custom: HashMap<String, f64>,
}

impl SimLogEntry {
    /// Create a new log entry with no custom fields.
    pub fn new(timestep: u64, sim_time: f64, wall_time: f64, energy: f64, n_bodies: usize) -> Self {
        Self {
            timestep,
            sim_time,
            wall_time,
            energy,
            n_bodies,
            custom: HashMap::new(),
        }
    }

    /// Insert a custom field value.
    pub fn insert(&mut self, key: impl Into<String>, value: f64) {
        self.custom.insert(key.into(), value);
    }

    /// Retrieve a custom field value by name.
    pub fn get(&self, key: &str) -> Option<f64> {
        self.custom.get(key).copied()
    }

    /// Serialise this entry to a CSV row (timestep, sim_time, wall_time,
    /// energy, n_bodies, then sorted custom keys).
    pub fn to_csv_row(&self, extra_keys: &[String]) -> String {
        let mut parts = vec![
            self.timestep.to_string(),
            format!("{:.9}", self.sim_time),
            format!("{:.9}", self.wall_time),
            format!("{:.9}", self.energy),
            self.n_bodies.to_string(),
        ];
        for key in extra_keys {
            if let Some(v) = self.custom.get(key) {
                parts.push(format!("{:.9}", v));
            } else {
                parts.push(String::new());
            }
        }
        parts.join(",")
    }
}

// ---------------------------------------------------------------------------
// SimLogger
// ---------------------------------------------------------------------------

/// Append-only simulation logger with CSV and JSON serialisation.
#[derive(Debug, Default, Clone)]
pub struct SimLogger {
    /// All logged entries in chronological order.
    pub entries: Vec<SimLogEntry>,
    /// How many entries to buffer before auto-flushing (0 = manual only).
    pub flush_interval: usize,
    /// Buffered output lines not yet "flushed" (written out).
    buffer: Vec<String>,
}

impl SimLogger {
    /// Create a new logger.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a logger with an auto-flush interval.
    pub fn with_flush_interval(flush_interval: usize) -> Self {
        Self {
            flush_interval,
            ..Default::default()
        }
    }

    /// Append a new entry.
    pub fn append(&mut self, entry: SimLogEntry) {
        self.entries.push(entry);
    }

    /// Return the number of logged entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if no entries have been logged.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Collect all unique custom field keys across all entries (sorted).
    pub fn all_custom_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = self
            .entries
            .iter()
            .flat_map(|e| e.custom.keys().cloned())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        keys.sort();
        keys
    }

    /// Serialise all entries to a CSV string.
    ///
    /// The first line is a header row.  Custom field columns appear in sorted
    /// key order.
    pub fn to_csv(&self) -> String {
        let keys = self.all_custom_keys();
        let mut lines = Vec::with_capacity(self.entries.len() + 1);
        let mut header = ["timestep", "sim_time", "wall_time", "energy", "n_bodies"]
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        header.extend(keys.iter().cloned());
        lines.push(header.join(","));
        for entry in &self.entries {
            lines.push(entry.to_csv_row(&keys));
        }
        lines.join("\n")
    }

    /// Serialise all entries to a JSON array string.
    pub fn to_json(&self) -> String {
        let items: Vec<String> = self.entries.iter().map(entry_to_json).collect();
        format!("[{}]", items.join(","))
    }

    /// Flush buffered lines and return them (clears the buffer).
    pub fn flush(&mut self) -> Vec<String> {
        std::mem::take(&mut self.buffer)
    }
}

/// Serialise one `SimLogEntry` to a JSON object string.
fn entry_to_json(e: &SimLogEntry) -> String {
    let mut parts = vec![
        format!("\"timestep\":{}", e.timestep),
        format!("\"sim_time\":{:.9}", e.sim_time),
        format!("\"wall_time\":{:.9}", e.wall_time),
        format!("\"energy\":{:.9}", e.energy),
        format!("\"n_bodies\":{}", e.n_bodies),
    ];
    let mut custom_pairs: Vec<(&String, &f64)> = e.custom.iter().collect();
    custom_pairs.sort_by_key(|(k, _)| k.as_str());
    for (k, v) in custom_pairs {
        let esc = k.replace('"', "\\\"");
        parts.push(format!("\"{}\":{:.9}", esc, v));
    }
    format!("{{{}}}", parts.join(","))
}

// ---------------------------------------------------------------------------
// EnergyMonitor
// ---------------------------------------------------------------------------

/// Tracks kinetic and potential energy; detects anomalous jumps.
#[derive(Debug, Clone, Default)]
pub struct EnergyMonitor {
    /// History of (kinetic, potential) energy pairs.
    pub history: Vec<(f64, f64)>,
    /// Relative jump threshold; a change > threshold * prev_energy triggers
    /// an anomaly.
    pub anomaly_threshold: f64,
    /// Indices of steps where an anomaly was detected.
    pub anomaly_steps: Vec<usize>,
}

impl EnergyMonitor {
    /// Create a new energy monitor with the given relative threshold.
    pub fn new(anomaly_threshold: f64) -> Self {
        Self {
            anomaly_threshold,
            ..Default::default()
        }
    }

    /// Record a new (kinetic, potential) energy measurement.
    pub fn record(&mut self, kinetic: f64, potential: f64) {
        let total = kinetic + potential;
        if let Some(&(ke_prev, pe_prev)) = self.history.last() {
            let prev_total = ke_prev + pe_prev;
            if prev_total.abs() > f64::EPSILON {
                let rel_change = (total - prev_total).abs() / prev_total.abs();
                if rel_change > self.anomaly_threshold {
                    self.anomaly_steps.push(self.history.len());
                }
            }
        }
        self.history.push((kinetic, potential));
    }

    /// Total energy at the most recent step.
    pub fn latest_total(&self) -> Option<f64> {
        self.history.last().map(|(ke, pe)| ke + pe)
    }

    /// Mean total energy over all recorded steps.
    pub fn mean_total(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.history.iter().map(|(ke, pe)| ke + pe).sum();
        sum / self.history.len() as f64
    }

    /// Returns `true` if any anomaly was detected.
    pub fn has_anomaly(&self) -> bool {
        !self.anomaly_steps.is_empty()
    }
}

// ---------------------------------------------------------------------------
// PerformanceLog
// ---------------------------------------------------------------------------

/// Timing (in seconds) for the major phases of one simulation step.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StepTiming {
    /// Time spent in broad-phase collision detection.
    pub broadphase: f64,
    /// Time spent in narrow-phase collision detection.
    pub narrowphase: f64,
    /// Time spent in constraint solve.
    pub solve: f64,
    /// Time spent integrating positions/velocities.
    pub integrate: f64,
}

impl StepTiming {
    /// Total time across all phases.
    pub fn total(&self) -> f64 {
        self.broadphase + self.narrowphase + self.solve + self.integrate
    }
}

/// Accumulates per-step timing data.
#[derive(Debug, Clone, Default)]
pub struct PerformanceLog {
    /// Recorded timings in chronological order.
    pub timings: Vec<StepTiming>,
}

impl PerformanceLog {
    /// Create a new performance log.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a timing record.
    pub fn push(&mut self, timing: StepTiming) {
        self.timings.push(timing);
    }

    /// Mean total step time across all recorded steps.
    pub fn mean_total(&self) -> f64 {
        if self.timings.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .timings
            .iter()
            .map(|t| t.total())
            .collect::<Vec<_>>()
            .iter()
            .sum();
        sum / self.timings.len() as f64
    }

    /// Maximum total step time.
    pub fn max_total(&self) -> f64 {
        self.timings
            .iter()
            .map(|t| t.total())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum total step time.
    pub fn min_total(&self) -> f64 {
        self.timings
            .iter()
            .map(|t| t.total())
            .fold(f64::INFINITY, f64::min)
    }
}

// ---------------------------------------------------------------------------
// LogFilter
// ---------------------------------------------------------------------------

/// Filters a slice of `SimLogEntry` by time range or custom field predicates.
#[derive(Debug, Clone, Default)]
pub struct LogFilter {
    /// Optional minimum sim_time (inclusive).
    pub sim_time_min: Option<f64>,
    /// Optional maximum sim_time (inclusive).
    pub sim_time_max: Option<f64>,
    /// Optional minimum energy (inclusive).
    pub energy_min: Option<f64>,
    /// Optional maximum energy (inclusive).
    pub energy_max: Option<f64>,
    /// Optional minimum n_bodies.
    pub n_bodies_min: Option<usize>,
}

impl LogFilter {
    /// Create an empty filter (accepts everything).
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the sim_time range.
    pub fn sim_time_range(mut self, min: f64, max: f64) -> Self {
        self.sim_time_min = Some(min);
        self.sim_time_max = Some(max);
        self
    }

    /// Set the energy range.
    pub fn energy_range(mut self, min: f64, max: f64) -> Self {
        self.energy_min = Some(min);
        self.energy_max = Some(max);
        self
    }

    /// Set minimum n_bodies filter.
    pub fn min_bodies(mut self, n: usize) -> Self {
        self.n_bodies_min = Some(n);
        self
    }

    /// Return all entries from `log` that pass this filter.
    pub fn apply<'a>(&self, log: &'a [SimLogEntry]) -> Vec<&'a SimLogEntry> {
        log.iter()
            .filter(|e| {
                if let Some(min) = self.sim_time_min
                    && e.sim_time < min
                {
                    return false;
                }
                if let Some(max) = self.sim_time_max
                    && e.sim_time > max
                {
                    return false;
                }
                if let Some(min) = self.energy_min
                    && e.energy < min
                {
                    return false;
                }
                if let Some(max) = self.energy_max
                    && e.energy > max
                {
                    return false;
                }
                if let Some(n) = self.n_bodies_min
                    && e.n_bodies < n
                {
                    return false;
                }
                true
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// SimSnapshot
// ---------------------------------------------------------------------------

/// Version identifier for snapshot format.
#[derive(Debug, Clone, PartialEq)]
pub struct SnapshotVersion {
    /// Major version number.
    pub major: u32,
    /// Minor version number.
    pub minor: u32,
}

impl SnapshotVersion {
    /// Current snapshot version.
    pub const CURRENT: SnapshotVersion = SnapshotVersion { major: 1, minor: 0 };
}

impl std::fmt::Display for SnapshotVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}", self.major, self.minor)
    }
}

/// A full serialisable snapshot of the simulation state.
#[derive(Debug, Clone)]
pub struct SimSnapshot {
    /// Format version for forward-compatibility checks.
    pub version: SnapshotVersion,
    /// Simulation time of this snapshot.
    pub sim_time: f64,
    /// Positions of all bodies `[[x,y,z\], ...]`.
    pub positions: Vec<[f64; 3]>,
    /// Velocities of all bodies `[[vx,vy,vz\], ...]`.
    pub velocities: Vec<[f64; 3]>,
    /// Masses of all bodies.
    pub masses: Vec<f64>,
    /// Optional metadata key-value pairs.
    pub metadata: HashMap<String, String>,
}

impl SimSnapshot {
    /// Create a new snapshot.
    pub fn new(
        sim_time: f64,
        positions: Vec<[f64; 3]>,
        velocities: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Self {
        Self {
            version: SnapshotVersion::CURRENT,
            sim_time,
            positions,
            velocities,
            masses,
            metadata: HashMap::new(),
        }
    }

    /// Number of bodies in this snapshot.
    pub fn n_bodies(&self) -> usize {
        self.positions.len()
    }

    /// Serialise to a JSON string.
    pub fn to_json(&self) -> String {
        let pos_json = array3_to_json_arr(&self.positions);
        let vel_json = array3_to_json_arr(&self.velocities);
        let mass_json = format!(
            "[{}]",
            self.masses
                .iter()
                .map(|m| format!("{:.9}", m))
                .collect::<Vec<_>>()
                .join(",")
        );
        let mut meta_parts: Vec<String> = self
            .metadata
            .iter()
            .map(|(k, v)| {
                let ek = k.replace('"', "\\\"");
                let ev = v.replace('"', "\\\"");
                format!("\"{}\":\"{}\"", ek, ev)
            })
            .collect();
        meta_parts.sort();
        let meta_json = format!("{{{}}}", meta_parts.join(","));
        format!(
            "{{\"version\":\"{}\",\"sim_time\":{:.9},\"n_bodies\":{},\"positions\":{},\"velocities\":{},\"masses\":{},\"metadata\":{}}}",
            self.version,
            self.sim_time,
            self.n_bodies(),
            pos_json,
            vel_json,
            mass_json,
            meta_json,
        )
    }
}

fn array3_to_json_arr(data: &[[f64; 3]]) -> String {
    let inner: Vec<String> = data
        .iter()
        .map(|v| format!("[{:.9},{:.9},{:.9}]", v[0], v[1], v[2]))
        .collect();
    format!("[{}]", inner.join(","))
}

// ---------------------------------------------------------------------------
// TelemetryStream — ring-buffer streaming log
// ---------------------------------------------------------------------------

/// Streaming telemetry log that keeps the last `capacity` entries in memory
/// while also accumulating all entries for bulk export.
#[derive(Debug, Clone)]
pub struct TelemetryStream {
    /// Ring buffer holding the last `capacity` entries.
    ring: Vec<Option<SimLogEntry>>,
    /// Next write position in the ring buffer.
    head: usize,
    /// Ring buffer capacity.
    capacity: usize,
    /// Total number of entries ever appended.
    total_count: u64,
    /// All entries ever appended (for full export).
    all_entries: Vec<SimLogEntry>,
}

impl TelemetryStream {
    /// Create a stream with the given ring-buffer capacity.
    pub fn new(capacity: usize) -> Self {
        assert!(capacity > 0, "capacity must be > 0");
        Self {
            ring: vec![None; capacity],
            head: 0,
            capacity,
            total_count: 0,
            all_entries: Vec::new(),
        }
    }

    /// Append a new entry to the stream.
    pub fn push(&mut self, entry: SimLogEntry) {
        self.all_entries.push(entry.clone());
        self.ring[self.head] = Some(entry);
        self.head = (self.head + 1) % self.capacity;
        self.total_count += 1;
    }

    /// Total number of entries ever pushed.
    pub fn total_count(&self) -> u64 {
        self.total_count
    }

    /// Iterate over the entries currently in the ring buffer (may be < capacity).
    pub fn ring_iter(&self) -> impl Iterator<Item = &SimLogEntry> {
        self.ring.iter().filter_map(|e| e.as_ref())
    }

    /// Number of entries currently in the ring buffer.
    pub fn ring_len(&self) -> usize {
        self.ring.iter().filter(|e| e.is_some()).count()
    }

    /// Return all entries ever pushed.
    pub fn all(&self) -> &[SimLogEntry] {
        &self.all_entries
    }
}

// ---------------------------------------------------------------------------
// SimMetrics
// ---------------------------------------------------------------------------

/// Computes statistical summaries over a named field in a log.
#[derive(Debug, Clone)]
pub struct SimMetrics {
    /// Collected values.
    values: Vec<f64>,
}

impl SimMetrics {
    /// Build metrics for the `energy` field of each entry.
    pub fn from_energy(entries: &[SimLogEntry]) -> Self {
        Self {
            values: entries.iter().map(|e| e.energy).collect(),
        }
    }

    /// Build metrics for a named custom field (entries missing the field are
    /// skipped).
    pub fn from_custom(entries: &[SimLogEntry], key: &str) -> Self {
        Self {
            values: entries.iter().filter_map(|e| e.get(key)).collect(),
        }
    }

    /// Build metrics from a raw slice of values.
    pub fn from_values(values: Vec<f64>) -> Self {
        Self { values }
    }

    /// Number of data points.
    pub fn count(&self) -> usize {
        self.values.len()
    }

    /// Mean value (returns `0.0` for empty input).
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().copied().sum::<f64>() / self.values.len() as f64
    }

    /// Population standard deviation.
    pub fn std_dev(&self) -> f64 {
        if self.values.len() < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var =
            self.values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / self.values.len() as f64;
        var.sqrt()
    }

    /// Minimum value.
    pub fn min(&self) -> f64 {
        self.values.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Maximum value.
    pub fn max(&self) -> f64 {
        self.values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Median value (returns `0.0` for empty input).
    pub fn median(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        let mut sorted = self.values.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = sorted.len();
        if n.is_multiple_of(2) {
            (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
        } else {
            sorted[n / 2]
        }
    }
}

// ---------------------------------------------------------------------------
// ReplayLog
// ---------------------------------------------------------------------------

/// Reads back a log and supports time-series reconstruction and interpolation.
#[derive(Debug, Clone)]
pub struct ReplayLog {
    /// Entries sorted by sim_time.
    entries: Vec<SimLogEntry>,
}

impl ReplayLog {
    /// Build a replay log from a vector of entries.  Entries are sorted by
    /// `sim_time`.
    pub fn new(mut entries: Vec<SimLogEntry>) -> Self {
        entries.sort_by(|a, b| {
            a.sim_time
                .partial_cmp(&b.sim_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Self { entries }
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if no entries are stored.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return the entry at index `i`.
    pub fn get(&self, i: usize) -> Option<&SimLogEntry> {
        self.entries.get(i)
    }

    /// Extract the full time series for `energy`.
    pub fn energy_series(&self) -> Vec<(f64, f64)> {
        self.entries
            .iter()
            .map(|e| (e.sim_time, e.energy))
            .collect()
    }

    /// Linearly interpolate the `energy` field at an arbitrary `t`.
    ///
    /// Returns `None` if the log is empty or `t` is outside the covered range.
    pub fn interpolate_energy(&self, t: f64) -> Option<f64> {
        if self.entries.is_empty() {
            return None;
        }
        // Clamp to range
        let first = self
            .entries
            .first()
            .expect("collection should not be empty");
        let last = self.entries.last().expect("collection should not be empty");
        if t <= first.sim_time {
            return Some(first.energy);
        }
        if t >= last.sim_time {
            return Some(last.energy);
        }
        // Binary search for bracketing interval
        let idx = self.entries.partition_point(|e| e.sim_time <= t);
        let e0 = &self.entries[idx - 1];
        let e1 = &self.entries[idx];
        let dt = e1.sim_time - e0.sim_time;
        if dt.abs() < f64::EPSILON {
            return Some(e0.energy);
        }
        let frac = (t - e0.sim_time) / dt;
        Some(e0.energy + frac * (e1.energy - e0.energy))
    }

    /// Iterate over all entries in sim_time order.
    pub fn iter(&self) -> impl Iterator<Item = &SimLogEntry> {
        self.entries.iter()
    }
}

// ---------------------------------------------------------------------------
// LogCompression
// ---------------------------------------------------------------------------

/// Run-length encoded representation of a sequence of f64 values.
#[derive(Debug, Clone, PartialEq)]
pub struct RleRun {
    /// The repeated value.
    pub value: f64,
    /// How many times it repeats.
    pub count: usize,
}

/// Delta-encoded timestamp sequence (stores differences between consecutive
/// values rather than absolute values).
#[derive(Debug, Clone)]
pub struct DeltaTimestamps {
    /// The first (absolute) timestamp.
    pub first: f64,
    /// Differences between consecutive timestamps.
    pub deltas: Vec<f64>,
}

impl DeltaTimestamps {
    /// Encode a slice of timestamps into delta form.
    pub fn encode(timestamps: &[f64]) -> Self {
        if timestamps.is_empty() {
            return Self {
                first: 0.0,
                deltas: Vec::new(),
            };
        }
        let first = timestamps[0];
        let deltas = timestamps.windows(2).map(|w| w[1] - w[0]).collect();
        Self { first, deltas }
    }

    /// Decode back to absolute timestamps.
    pub fn decode(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.deltas.len() + 1);
        out.push(self.first);
        let mut cur = self.first;
        for &d in &self.deltas {
            cur += d;
            out.push(cur);
        }
        out
    }

    /// Number of timestamps represented.
    pub fn len(&self) -> usize {
        self.deltas.len() + 1
    }

    /// Returns `true` if empty (no timestamps).
    pub fn is_empty(&self) -> bool {
        self.deltas.is_empty() && self.first == 0.0
    }
}

/// Run-length encode a sequence of f64 values.
///
/// Consecutive equal values (bitwise) are collapsed into a single `RleRun`.
pub fn rle_encode(values: &[f64]) -> Vec<RleRun> {
    let mut runs = Vec::new();
    let mut iter = values.iter().copied();
    let Some(first) = iter.next() else {
        return runs;
    };
    let mut cur = first;
    let mut count = 1usize;
    for v in iter {
        if v.to_bits() == cur.to_bits() {
            count += 1;
        } else {
            runs.push(RleRun { value: cur, count });
            cur = v;
            count = 1;
        }
    }
    runs.push(RleRun { value: cur, count });
    runs
}

/// Decode a run-length encoded sequence back to flat values.
pub fn rle_decode(runs: &[RleRun]) -> Vec<f64> {
    runs.iter()
        .flat_map(|r| std::iter::repeat_n(r.value, r.count))
        .collect()
}

/// Compression statistics for a log field.
#[derive(Debug, Clone)]
pub struct LogCompression {
    /// Original length.
    pub original_len: usize,
    /// Number of RLE runs.
    pub rle_runs: usize,
    /// Compression ratio (original / runs).
    pub ratio: f64,
}

impl LogCompression {
    /// Analyse compression potential for a sequence of values.
    pub fn analyse(values: &[f64]) -> Self {
        let runs = rle_encode(values);
        let original_len = values.len();
        let rle_runs = runs.len();
        let ratio = if rle_runs == 0 {
            1.0
        } else {
            original_len as f64 / rle_runs as f64
        };
        Self {
            original_len,
            rle_runs,
            ratio,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(step: u64, t: f64, energy: f64) -> SimLogEntry {
        SimLogEntry::new(step, t, t * 0.001, energy, 100)
    }

    // -- SimLogEntry --

    #[test]
    fn test_entry_new_fields() {
        let e = make_entry(5, 0.5, 42.0);
        assert_eq!(e.timestep, 5);
        assert!((e.sim_time - 0.5).abs() < 1e-12);
        assert!((e.energy - 42.0).abs() < 1e-12);
        assert_eq!(e.n_bodies, 100);
        assert!(e.custom.is_empty());
    }

    #[test]
    fn test_entry_insert_and_get() {
        let mut e = make_entry(0, 0.0, 0.0);
        e.insert("temp", 300.0);
        assert!((e.get("temp").unwrap() - 300.0).abs() < 1e-12);
        assert!(e.get("missing").is_none());
    }

    #[test]
    fn test_entry_to_csv_row_no_custom() {
        let e = SimLogEntry::new(1, 0.1, 0.0001, 50.0, 10);
        let row = e.to_csv_row(&[]);
        assert!(row.starts_with("1,"));
        assert!(row.contains("50.0") || row.contains("50."));
    }

    #[test]
    fn test_entry_to_csv_row_with_custom() {
        let mut e = make_entry(2, 1.0, 10.0);
        e.insert("foo", 3.125);
        let keys = vec!["foo".to_string()];
        let row = e.to_csv_row(&keys);
        assert!(row.contains("3.125"));
    }

    #[test]
    fn test_entry_to_csv_row_missing_custom_key() {
        let e = make_entry(0, 0.0, 0.0);
        let keys = vec!["absent".to_string()];
        let row = e.to_csv_row(&keys);
        // Last field should be empty
        assert!(row.ends_with(','));
    }

    // -- SimLogger --

    #[test]
    fn test_logger_append_and_len() {
        let mut logger = SimLogger::new();
        assert!(logger.is_empty());
        logger.append(make_entry(0, 0.0, 100.0));
        assert_eq!(logger.len(), 1);
        assert!(!logger.is_empty());
    }

    #[test]
    fn test_logger_to_csv_header() {
        let mut logger = SimLogger::new();
        logger.append(make_entry(0, 0.0, 1.0));
        let csv = logger.to_csv();
        let first_line = csv.lines().next().unwrap();
        assert!(first_line.contains("timestep"));
        assert!(first_line.contains("energy"));
    }

    #[test]
    fn test_logger_to_csv_row_count() {
        let mut logger = SimLogger::new();
        for i in 0..5 {
            logger.append(make_entry(i, i as f64 * 0.1, 10.0));
        }
        let csv = logger.to_csv();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 6); // 1 header + 5 data
    }

    #[test]
    fn test_logger_to_json_is_array() {
        let mut logger = SimLogger::new();
        logger.append(make_entry(0, 0.0, 5.0));
        let json = logger.to_json();
        assert!(json.starts_with('['));
        assert!(json.ends_with(']'));
    }

    #[test]
    fn test_logger_all_custom_keys_sorted() {
        let mut logger = SimLogger::new();
        let mut e = make_entry(0, 0.0, 0.0);
        e.insert("z_key", 1.0);
        e.insert("a_key", 2.0);
        logger.append(e);
        let keys = logger.all_custom_keys();
        assert_eq!(keys[0], "a_key");
        assert_eq!(keys[1], "z_key");
    }

    #[test]
    fn test_logger_flush_buffer() {
        let mut logger = SimLogger::new();
        let flushed = logger.flush();
        assert!(flushed.is_empty());
    }

    // -- EnergyMonitor --

    #[test]
    fn test_energy_monitor_no_anomaly_steady() {
        let mut mon = EnergyMonitor::new(0.1);
        for _ in 0..10 {
            mon.record(50.0, 50.0); // constant total = 100
        }
        assert!(!mon.has_anomaly());
    }

    #[test]
    fn test_energy_monitor_detects_jump() {
        let mut mon = EnergyMonitor::new(0.1); // 10% threshold
        mon.record(50.0, 50.0); // total = 100
        mon.record(70.0, 80.0); // total = 150 → +50% jump → anomaly
        assert!(mon.has_anomaly());
    }

    #[test]
    fn test_energy_monitor_latest_total() {
        let mut mon = EnergyMonitor::new(1.0);
        mon.record(30.0, 20.0);
        assert!((mon.latest_total().unwrap() - 50.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_mean_total() {
        let mut mon = EnergyMonitor::new(1.0);
        mon.record(10.0, 10.0); // 20
        mon.record(20.0, 20.0); // 40
        assert!((mon.mean_total() - 30.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_empty() {
        let mon = EnergyMonitor::new(0.1);
        assert!(mon.latest_total().is_none());
        assert!((mon.mean_total()).abs() < 1e-12);
        assert!(!mon.has_anomaly());
    }

    // -- PerformanceLog --

    #[test]
    fn test_perf_log_step_timing_total() {
        let t = StepTiming {
            broadphase: 0.001,
            narrowphase: 0.002,
            solve: 0.003,
            integrate: 0.001,
        };
        assert!((t.total() - 0.007).abs() < 1e-12);
    }

    #[test]
    fn test_perf_log_mean_total() {
        let mut log = PerformanceLog::new();
        log.push(StepTiming {
            broadphase: 0.01,
            narrowphase: 0.0,
            solve: 0.0,
            integrate: 0.0,
        });
        log.push(StepTiming {
            broadphase: 0.03,
            narrowphase: 0.0,
            solve: 0.0,
            integrate: 0.0,
        });
        assert!((log.mean_total() - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_perf_log_max_min() {
        let mut log = PerformanceLog::new();
        log.push(StepTiming {
            broadphase: 0.005,
            narrowphase: 0.0,
            solve: 0.0,
            integrate: 0.0,
        });
        log.push(StepTiming {
            broadphase: 0.020,
            narrowphase: 0.0,
            solve: 0.0,
            integrate: 0.0,
        });
        assert!((log.max_total() - 0.020).abs() < 1e-12);
        assert!((log.min_total() - 0.005).abs() < 1e-12);
    }

    #[test]
    fn test_perf_log_empty_mean() {
        let log = PerformanceLog::new();
        assert!((log.mean_total()).abs() < 1e-12);
    }

    // -- LogFilter --

    #[test]
    fn test_filter_sim_time_range() {
        let entries = vec![
            make_entry(0, 0.0, 10.0),
            make_entry(1, 1.0, 10.0),
            make_entry(2, 2.0, 10.0),
        ];
        let filter = LogFilter::new().sim_time_range(0.5, 1.5);
        let filtered = filter.apply(&entries);
        assert_eq!(filtered.len(), 1);
        assert!((filtered[0].sim_time - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_filter_energy_range() {
        let entries = vec![
            make_entry(0, 0.0, 5.0),
            make_entry(1, 1.0, 50.0),
            make_entry(2, 2.0, 500.0),
        ];
        let filter = LogFilter::new().energy_range(10.0, 100.0);
        let filtered = filter.apply(&entries);
        assert_eq!(filtered.len(), 1);
        assert!((filtered[0].energy - 50.0).abs() < 1e-12);
    }

    #[test]
    fn test_filter_min_bodies() {
        let mut entries = vec![make_entry(0, 0.0, 10.0)];
        entries[0].n_bodies = 5;
        let mut e2 = make_entry(1, 1.0, 10.0);
        e2.n_bodies = 200;
        entries.push(e2);
        let filter = LogFilter::new().min_bodies(100);
        let filtered = filter.apply(&entries);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].n_bodies, 200);
    }

    #[test]
    fn test_filter_empty_accepts_all() {
        let entries: Vec<SimLogEntry> = (0..5).map(|i| make_entry(i, i as f64, 1.0)).collect();
        let filter = LogFilter::new();
        assert_eq!(filter.apply(&entries).len(), 5);
    }

    // -- SimSnapshot --

    #[test]
    fn test_snapshot_n_bodies() {
        let snap = SimSnapshot::new(
            1.0,
            vec![[0.0; 3], [1.0; 3]],
            vec![[0.0; 3], [0.0; 3]],
            vec![1.0, 2.0],
        );
        assert_eq!(snap.n_bodies(), 2);
    }

    #[test]
    fn test_snapshot_to_json_contains_version() {
        let snap = SimSnapshot::new(0.5, vec![], vec![], vec![]);
        let json = snap.to_json();
        assert!(json.contains("\"version\""));
        assert!(json.contains("1.0"));
    }

    #[test]
    fn test_snapshot_to_json_contains_n_bodies() {
        let snap = SimSnapshot::new(0.0, vec![[1.0, 2.0, 3.0]], vec![[0.1, 0.2, 0.3]], vec![5.0]);
        let json = snap.to_json();
        assert!(json.contains("\"n_bodies\":1"));
    }

    #[test]
    fn test_snapshot_version_display() {
        assert_eq!(SnapshotVersion::CURRENT.to_string(), "1.0");
    }

    // -- TelemetryStream --

    #[test]
    fn test_telemetry_stream_basic_push() {
        let mut stream = TelemetryStream::new(3);
        stream.push(make_entry(0, 0.0, 1.0));
        assert_eq!(stream.total_count(), 1);
        assert_eq!(stream.ring_len(), 1);
    }

    #[test]
    fn test_telemetry_stream_ring_wraps() {
        let mut stream = TelemetryStream::new(3);
        for i in 0..5u64 {
            stream.push(make_entry(i, i as f64, 1.0));
        }
        assert_eq!(stream.total_count(), 5);
        assert_eq!(stream.ring_len(), 3); // only 3 slots
    }

    #[test]
    fn test_telemetry_stream_all_entries() {
        let mut stream = TelemetryStream::new(2);
        for i in 0..10u64 {
            stream.push(make_entry(i, i as f64, 0.0));
        }
        assert_eq!(stream.all().len(), 10);
    }

    #[test]
    fn test_telemetry_stream_ring_iter() {
        let mut stream = TelemetryStream::new(5);
        stream.push(make_entry(0, 0.0, 42.0));
        let in_ring: Vec<_> = stream.ring_iter().collect();
        assert_eq!(in_ring.len(), 1);
        assert!((in_ring[0].energy - 42.0).abs() < 1e-12);
    }

    // -- SimMetrics --

    #[test]
    fn test_metrics_mean() {
        let m = SimMetrics::from_values(vec![1.0, 2.0, 3.0]);
        assert!((m.mean() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_std_dev() {
        let m = SimMetrics::from_values(vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);
        // Population std ≈ 2.0
        assert!((m.std_dev() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_metrics_min_max() {
        let m = SimMetrics::from_values(vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0]);
        assert!((m.min() - 1.0).abs() < 1e-12);
        assert!((m.max() - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_median_odd() {
        let m = SimMetrics::from_values(vec![5.0, 1.0, 3.0]);
        assert!((m.median() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_median_even() {
        let m = SimMetrics::from_values(vec![1.0, 2.0, 3.0, 4.0]);
        assert!((m.median() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_from_energy() {
        let entries: Vec<SimLogEntry> = vec![make_entry(0, 0.0, 10.0), make_entry(1, 1.0, 20.0)];
        let m = SimMetrics::from_energy(&entries);
        assert_eq!(m.count(), 2);
        assert!((m.mean() - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_from_custom_key() {
        let mut e1 = make_entry(0, 0.0, 0.0);
        e1.insert("temp", 100.0);
        let mut e2 = make_entry(1, 1.0, 0.0);
        e2.insert("temp", 200.0);
        let entries = vec![e1, e2];
        let m = SimMetrics::from_custom(&entries, "temp");
        assert_eq!(m.count(), 2);
        assert!((m.mean() - 150.0).abs() < 1e-12);
    }

    #[test]
    fn test_metrics_empty() {
        let m = SimMetrics::from_values(vec![]);
        assert_eq!(m.count(), 0);
        assert!((m.mean()).abs() < 1e-12);
        assert!((m.std_dev()).abs() < 1e-12);
        assert!((m.median()).abs() < 1e-12);
    }

    // -- ReplayLog --

    #[test]
    fn test_replay_sorted_on_construction() {
        let entries = vec![
            make_entry(2, 2.0, 20.0),
            make_entry(0, 0.0, 0.0),
            make_entry(1, 1.0, 10.0),
        ];
        let replay = ReplayLog::new(entries);
        assert!((replay.get(0).unwrap().sim_time - 0.0).abs() < 1e-12);
        assert!((replay.get(2).unwrap().sim_time - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_energy_series_len() {
        let entries: Vec<_> = (0..5)
            .map(|i| make_entry(i, i as f64, i as f64 * 10.0))
            .collect();
        let replay = ReplayLog::new(entries);
        assert_eq!(replay.energy_series().len(), 5);
    }

    #[test]
    fn test_replay_interpolate_energy_midpoint() {
        let entries = vec![make_entry(0, 0.0, 0.0), make_entry(1, 2.0, 20.0)];
        let replay = ReplayLog::new(entries);
        let e = replay.interpolate_energy(1.0).unwrap();
        assert!((e - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_replay_interpolate_clamp_before_start() {
        let entries = vec![make_entry(0, 1.0, 5.0), make_entry(1, 2.0, 10.0)];
        let replay = ReplayLog::new(entries);
        assert!((replay.interpolate_energy(0.0).unwrap() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_interpolate_clamp_after_end() {
        let entries = vec![make_entry(0, 0.0, 5.0), make_entry(1, 1.0, 10.0)];
        let replay = ReplayLog::new(entries);
        assert!((replay.interpolate_energy(99.0).unwrap() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_replay_empty() {
        let replay = ReplayLog::new(vec![]);
        assert!(replay.is_empty());
        assert!(replay.interpolate_energy(0.5).is_none());
    }

    // -- LogCompression / RLE --

    #[test]
    fn test_rle_encode_constant() {
        let data = vec![3.125; 5];
        let runs = rle_encode(&data);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].count, 5);
        assert!((runs[0].value - 3.125).abs() < 1e-12);
    }

    #[test]
    fn test_rle_decode_roundtrip() {
        let data = vec![1.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        let runs = rle_encode(&data);
        let decoded = rle_decode(&runs);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_rle_encode_all_unique() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let runs = rle_encode(&data);
        assert_eq!(runs.len(), 4);
    }

    #[test]
    fn test_rle_encode_empty() {
        let runs = rle_encode(&[]);
        assert!(runs.is_empty());
    }

    #[test]
    fn test_log_compression_analyse_ratio() {
        let data = vec![42.0f64; 100];
        let comp = LogCompression::analyse(&data);
        assert_eq!(comp.original_len, 100);
        assert_eq!(comp.rle_runs, 1);
        assert!((comp.ratio - 100.0).abs() < 1e-9);
    }

    // -- DeltaTimestamps --

    #[test]
    fn test_delta_timestamps_encode_decode_roundtrip() {
        let ts = vec![0.0, 0.1, 0.2, 0.4, 0.7, 1.1];
        let enc = DeltaTimestamps::encode(&ts);
        let dec = enc.decode();
        for (a, b) in ts.iter().zip(dec.iter()) {
            assert!((a - b).abs() < 1e-10, "a={a} b={b}");
        }
    }

    #[test]
    fn test_delta_timestamps_len() {
        let ts = vec![0.0, 1.0, 2.0, 3.0];
        let enc = DeltaTimestamps::encode(&ts);
        assert_eq!(enc.len(), 4);
    }

    #[test]
    fn test_delta_timestamps_empty() {
        let enc = DeltaTimestamps::encode(&[]);
        assert_eq!(enc.len(), 1); // first=0, no deltas → decode gives [0.0]
    }
}
