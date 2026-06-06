// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SQLite-like in-memory simulation database with time-indexed snapshots,
//! run-length-encoded compression, differential storage, event logging,
//! result aggregation, and CSV/JSON import-export.
//!
//! # Overview
//!
//! - [`SimulationRecord`] — a single simulation run with parameters and metadata.
//! - [`SimulationDatabase`] — a collection of records with CSV/JSON persistence.
//! - [`SnapshotTable`] — time-indexed snapshots with RLE compression and
//!   differential storage.
//! - [`EventLog`] — timestamped event records for a simulation session.
//! - [`ResultAggregator`] — computes min/max/mean/std over snapshot collections.
//! - [`ProvenanceTracker`] — W3C-PROV-style provenance graph.
//! - [`CheckpointManager`] — rolling-window checkpoint storage.
//! - [`ParameterSweep`] — Cartesian-product and Latin-hypercube parameter sweeps.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// SimulationRecord
// ---------------------------------------------------------------------------

/// A single simulation run record storing parameters and metadata.
#[derive(Debug, Clone)]
pub struct SimulationRecord {
    /// Unique identifier for this simulation run.
    pub id: String,
    /// Unix timestamp (seconds since epoch) when the simulation was created.
    pub timestamp: u64,
    /// Numeric simulation parameters (e.g. `"dt"`, `"viscosity"`).
    pub parameters: HashMap<String, f64>,
    /// Arbitrary string metadata (e.g. `"solver"`, `"git_hash"`).
    pub metadata: HashMap<String, String>,
}

impl SimulationRecord {
    /// Construct a new record with the given id and timestamp.
    pub fn new(id: impl Into<String>, timestamp: u64) -> Self {
        Self {
            id: id.into(),
            timestamp,
            parameters: HashMap::new(),
            metadata: HashMap::new(),
        }
    }

    /// Set a numeric parameter, overwriting any existing value.
    pub fn set_param(&mut self, key: impl Into<String>, value: f64) {
        self.parameters.insert(key.into(), value);
    }

    /// Set a string metadata entry.
    pub fn set_meta(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.metadata.insert(key.into(), value.into());
    }
}

// ---------------------------------------------------------------------------
// SimulationDatabase
// ---------------------------------------------------------------------------

/// In-memory database of simulation records with CSV/JSON persistence.
#[derive(Debug, Default)]
pub struct SimulationDatabase {
    /// All stored simulation records.
    pub records: Vec<SimulationRecord>,
    /// Path to the backing CSV file (informational; used by save/load helpers).
    pub file_path: String,
}

impl SimulationDatabase {
    /// Create a new, empty database associated with the given file path.
    pub fn new(file_path: impl Into<String>) -> Self {
        Self {
            records: Vec::new(),
            file_path: file_path.into(),
        }
    }

    /// Append a record to the database.
    pub fn add_record(&mut self, record: SimulationRecord) {
        self.records.push(record);
    }

    /// Find a record by its unique id, returning a reference or `None`.
    pub fn find_by_id(&self, id: &str) -> Option<&SimulationRecord> {
        self.records.iter().find(|r| r.id == id)
    }

    /// Return all records where `param` lies in the inclusive range `[lo, hi]`.
    pub fn query_range(&self, param: &str, lo: f64, hi: f64) -> Vec<&SimulationRecord> {
        self.records
            .iter()
            .filter(|r| r.parameters.get(param).is_some_and(|&v| v >= lo && v <= hi))
            .collect()
    }

    /// Return all records whose timestamp falls in `[t_lo, t_hi]`.
    pub fn query_time_range(&self, t_lo: u64, t_hi: u64) -> Vec<&SimulationRecord> {
        self.records
            .iter()
            .filter(|r| r.timestamp >= t_lo && r.timestamp <= t_hi)
            .collect()
    }

    /// Delete all records matching the given id.  Returns the number of
    /// records removed.
    pub fn delete_by_id(&mut self, id: &str) -> usize {
        let before = self.records.len();
        self.records.retain(|r| r.id != id);
        before - self.records.len()
    }

    /// Serialize all records to a CSV string.
    ///
    /// Format: `id,timestamp,key=value;…,metakey=metaval;…`
    pub fn save_to_csv(&self) -> String {
        let mut out = String::from("id,timestamp,parameters,metadata\n");
        for r in &self.records {
            let params: Vec<String> = r
                .parameters
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect();
            let meta: Vec<String> = r.metadata.iter().map(|(k, v)| format!("{k}={v}")).collect();
            out.push_str(&format!(
                "{},{},{},{}\n",
                r.id,
                r.timestamp,
                params.join(";"),
                meta.join(";")
            ));
        }
        out
    }

    /// Populate the database by parsing a CSV string produced by `save_to_csv`.
    ///
    /// Clears existing records before loading.
    pub fn load_from_csv(&mut self, s: &str) {
        self.records.clear();
        for line in s.lines().skip(1) {
            let parts: Vec<&str> = line.splitn(4, ',').collect();
            if parts.len() < 2 {
                continue;
            }
            let id = parts[0].to_string();
            let timestamp: u64 = parts[1].parse().unwrap_or(0);
            let mut record = SimulationRecord::new(id, timestamp);
            if parts.len() > 2 && !parts[2].is_empty() {
                for pair in parts[2].split(';') {
                    let kv: Vec<&str> = pair.splitn(2, '=').collect();
                    if kv.len() == 2
                        && let Ok(v) = kv[1].parse::<f64>()
                    {
                        record.parameters.insert(kv[0].to_string(), v);
                    }
                }
            }
            if parts.len() > 3 && !parts[3].is_empty() {
                for pair in parts[3].split(';') {
                    let kv: Vec<&str> = pair.splitn(2, '=').collect();
                    if kv.len() == 2 {
                        record.metadata.insert(kv[0].to_string(), kv[1].to_string());
                    }
                }
            }
            self.records.push(record);
        }
    }

    /// Compute `(min, max, mean)` for the named parameter across all records
    /// that contain it. Returns `(0.0, 0.0, 0.0)` if no records match.
    pub fn statistics(&self, param: &str) -> (f64, f64, f64) {
        let values: Vec<f64> = self
            .records
            .iter()
            .filter_map(|r| r.parameters.get(param).copied())
            .collect();
        if values.is_empty() {
            return (0.0, 0.0, 0.0);
        }
        let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        (min, max, mean)
    }

    /// Export all records as a compact JSON string.
    pub fn export_json(&self) -> String {
        let mut out = String::from("[\n");
        for (i, r) in self.records.iter().enumerate() {
            out.push_str("  {\n");
            out.push_str(&format!("    \"id\": \"{}\",\n", r.id));
            out.push_str(&format!("    \"timestamp\": {},\n", r.timestamp));
            out.push_str("    \"parameters\": {");
            let params: Vec<String> = r
                .parameters
                .iter()
                .map(|(k, v)| format!("\"{k}\": {v}"))
                .collect();
            out.push_str(&params.join(", "));
            out.push_str("},\n");
            out.push_str("    \"metadata\": {");
            let meta: Vec<String> = r
                .metadata
                .iter()
                .map(|(k, v)| format!("\"{k}\": \"{v}\""))
                .collect();
            out.push_str(&meta.join(", "));
            out.push_str("}\n");
            if i + 1 < self.records.len() {
                out.push_str("  },\n");
            } else {
                out.push_str("  }\n");
            }
        }
        out.push(']');
        out
    }

    /// Import records from a minimal JSON array string as produced by
    /// `export_json`.  Only the `id` and `timestamp` fields are parsed in
    /// this lightweight implementation; full parameter parsing requires a
    /// proper JSON library.
    pub fn import_json_ids(&mut self, json: &str) {
        // Very lightweight scanner: extract "id": "VALUE" and "timestamp": N pairs.
        for chunk in json.split('{') {
            let id = extract_json_str(chunk, "\"id\"");
            let ts_str = extract_json_number(chunk, "\"timestamp\"");
            if let Some(id) = id {
                let ts: u64 = ts_str.unwrap_or_default().parse().unwrap_or(0);
                self.records.push(SimulationRecord::new(id, ts));
            }
        }
    }

    /// Return the total number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Return `true` if the database holds no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Extract the first string value for `key` from a JSON fragment.
fn extract_json_str(chunk: &str, key: &str) -> Option<String> {
    let pos = chunk.find(key)?;
    let rest = &chunk[pos + key.len()..];
    let colon = rest.find(':')? + 1;
    let rest2 = rest[colon..].trim_start();
    if !rest2.starts_with('"') {
        return None;
    }
    let inner = &rest2[1..];
    let end = inner.find('"')?;
    Some(inner[..end].to_string())
}

/// Extract the first numeric (integer) value for `key` from a JSON fragment.
fn extract_json_number(chunk: &str, key: &str) -> Option<String> {
    let pos = chunk.find(key)?;
    let rest = &chunk[pos + key.len()..];
    let colon = rest.find(':')? + 1;
    let rest2 = rest[colon..].trim_start();
    let end = rest2
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest2.len());
    if end == 0 {
        return None;
    }
    Some(rest2[..end].to_string())
}

// ---------------------------------------------------------------------------
// Run-length encoding helpers
// ---------------------------------------------------------------------------

/// Encode a slice of `f64` values using run-length encoding.
///
/// Returns a `Vec<(f64, usize)>` where each tuple is `(value, count)`.
pub fn rle_encode(data: &[f64]) -> Vec<(f64, usize)> {
    if data.is_empty() {
        return vec![];
    }
    let mut result = Vec::new();
    let mut current = data[0];
    let mut count = 1usize;
    for &v in &data[1..] {
        if (v - current).abs() < f64::EPSILON {
            count += 1;
        } else {
            result.push((current, count));
            current = v;
            count = 1;
        }
    }
    result.push((current, count));
    result
}

/// Decode a run-length-encoded slice back to a flat `Vec`f64`.
pub fn rle_decode(encoded: &[(f64, usize)]) -> Vec<f64> {
    let mut result = Vec::new();
    for &(v, n) in encoded {
        for _ in 0..n {
            result.push(v);
        }
    }
    result
}

/// Compression ratio: `original_len / encoded_len` (number of (value, count) pairs).
///
/// Returns `1.0` if either side is zero.
pub fn rle_compression_ratio(original_len: usize, encoded: &[(f64, usize)]) -> f64 {
    if encoded.is_empty() || original_len == 0 {
        return 1.0;
    }
    original_len as f64 / encoded.len() as f64
}

// ---------------------------------------------------------------------------
// SnapshotTable — time-indexed snapshot storage
// ---------------------------------------------------------------------------

/// A single simulation snapshot: a named collection of scalar fields
/// recorded at a given simulation time.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// Simulation time at which this snapshot was taken.
    pub time: f64,
    /// Named field arrays (e.g. `"velocity_x"`, `"pressure"`).
    pub fields: HashMap<String, Vec<f64>>,
}

impl Snapshot {
    /// Create an empty snapshot at the given time.
    pub fn new(time: f64) -> Self {
        Self {
            time,
            fields: HashMap::new(),
        }
    }

    /// Add or overwrite a named field.
    pub fn set_field(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.fields.insert(name.into(), data);
    }

    /// Return the length of the first field array, or 0 if no fields exist.
    pub fn node_count(&self) -> usize {
        self.fields.values().next().map_or(0, |v| v.len())
    }
}

/// Differential storage entry: stores only the indices and new values
/// that differ from the previous snapshot.
#[derive(Debug, Clone)]
pub struct DiffEntry {
    /// Simulation time this diff applies to.
    pub time: f64,
    /// Field name.
    pub field: String,
    /// Indices of changed values.
    pub indices: Vec<usize>,
    /// New values at those indices.
    pub new_values: Vec<f64>,
}

impl DiffEntry {
    /// Create a new diff entry.
    pub fn new(time: f64, field: impl Into<String>) -> Self {
        Self {
            time,
            field: field.into(),
            indices: Vec::new(),
            new_values: Vec::new(),
        }
    }
}

/// Time-indexed table of simulation snapshots with optional RLE compression
/// and differential storage.
#[derive(Debug, Default)]
pub struct SnapshotTable {
    /// Stored snapshots ordered by insertion time.
    pub snapshots: Vec<Snapshot>,
    /// Differential entries for fields that change sparsely between snapshots.
    pub diffs: Vec<DiffEntry>,
    /// RLE-compressed versions of field arrays, keyed by `"field@time_idx"`.
    pub compressed: HashMap<String, Vec<(f64, usize)>>,
}

impl SnapshotTable {
    /// Create an empty table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a snapshot; snapshots are kept in ascending time order.
    pub fn insert(&mut self, snap: Snapshot) {
        let pos = self.snapshots.partition_point(|s| s.time < snap.time);
        self.snapshots.insert(pos, snap);
    }

    /// Return all snapshots whose time falls in `\[t_lo, t_hi\]`.
    pub fn query_time_range(&self, t_lo: f64, t_hi: f64) -> Vec<&Snapshot> {
        self.snapshots
            .iter()
            .filter(|s| s.time >= t_lo && s.time <= t_hi)
            .collect()
    }

    /// Retrieve a single snapshot whose time is closest to `t`.
    pub fn nearest(&self, t: f64) -> Option<&Snapshot> {
        self.snapshots.iter().min_by(|a, b| {
            (a.time - t)
                .abs()
                .partial_cmp(&(b.time - t).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Compress a named field of the snapshot at index `snap_idx` using RLE.
    ///
    /// The compressed data is stored internally under the key
    /// `"`field`@`snap_idx`"`.
    pub fn compress_field(&mut self, snap_idx: usize, field: &str) {
        if let Some(snap) = self.snapshots.get(snap_idx)
            && let Some(data) = snap.fields.get(field)
        {
            let encoded = rle_encode(data);
            let key = format!("{field}@{snap_idx}");
            self.compressed.insert(key, encoded);
        }
    }

    /// Decompress a previously compressed field, returning the flat `Vec`f64`.
    ///
    /// Returns `None` if the field was not compressed.
    pub fn decompress_field(&self, snap_idx: usize, field: &str) -> Option<Vec<f64>> {
        let key = format!("{field}@{snap_idx}");
        self.compressed.get(&key).map(|enc| rle_decode(enc))
    }

    /// Compute the differential between consecutive snapshots for a named field
    /// and store the result in `self.diffs`.
    ///
    /// Returns the number of changed elements.
    pub fn compute_diff(&mut self, field: &str, snap_idx: usize) -> usize {
        if snap_idx == 0 || snap_idx >= self.snapshots.len() {
            return 0;
        }
        let prev = self.snapshots[snap_idx - 1]
            .fields
            .get(field)
            .cloned()
            .unwrap_or_default();
        let curr = self.snapshots[snap_idx]
            .fields
            .get(field)
            .cloned()
            .unwrap_or_default();
        let time = self.snapshots[snap_idx].time;
        let mut entry = DiffEntry::new(time, field);
        for (i, (&p, &c)) in prev.iter().zip(curr.iter()).enumerate() {
            if (c - p).abs() > f64::EPSILON {
                entry.indices.push(i);
                entry.new_values.push(c);
            }
        }
        let changed = entry.indices.len();
        self.diffs.push(entry);
        changed
    }

    /// Apply a stored diff to a base field array, returning the updated array.
    ///
    /// `base` is modified in-place by overwriting the changed indices.
    pub fn apply_diff(base: &mut [f64], diff: &DiffEntry) {
        for (&idx, &val) in diff.indices.iter().zip(diff.new_values.iter()) {
            if idx < base.len() {
                base[idx] = val;
            }
        }
    }

    /// Return the number of snapshots.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    /// Return `true` if the table is empty.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }
}

// ---------------------------------------------------------------------------
// EventLog — timestamped event logging
// ---------------------------------------------------------------------------

/// Severity level for a logged event.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EventLevel {
    /// Diagnostic trace messages.
    Debug,
    /// General informational messages.
    Info,
    /// Non-fatal warnings.
    Warning,
    /// Recoverable errors.
    Error,
    /// Fatal errors that stop the simulation.
    Critical,
}

impl EventLevel {
    /// Return the level as a static string label.
    pub fn as_str(self) -> &'static str {
        match self {
            EventLevel::Debug => "DEBUG",
            EventLevel::Info => "INFO",
            EventLevel::Warning => "WARNING",
            EventLevel::Error => "ERROR",
            EventLevel::Critical => "CRITICAL",
        }
    }
}

/// A single logged event.
#[derive(Debug, Clone)]
pub struct LogEvent {
    /// Simulation time at which the event occurred.
    pub sim_time: f64,
    /// Wall-clock timestamp (Unix seconds).
    pub wall_time: u64,
    /// Severity level.
    pub level: EventLevel,
    /// Category or subsystem name (e.g. `"solver"`, `"io"`).
    pub category: String,
    /// Human-readable message.
    pub message: String,
}

impl LogEvent {
    /// Construct a new event record.
    pub fn new(
        sim_time: f64,
        wall_time: u64,
        level: EventLevel,
        category: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            sim_time,
            wall_time,
            level,
            category: category.into(),
            message: message.into(),
        }
    }
}

/// Append-only timestamped event log for a simulation session.
#[derive(Debug, Default)]
pub struct EventLog {
    /// All recorded events.
    pub events: Vec<LogEvent>,
}

impl EventLog {
    /// Create an empty event log.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an event.
    pub fn log(&mut self, event: LogEvent) {
        self.events.push(event);
    }

    /// Convenience helper: log a message at `Info` level.
    pub fn info(&mut self, sim_time: f64, category: impl Into<String>, message: impl Into<String>) {
        self.log(LogEvent::new(
            sim_time,
            0,
            EventLevel::Info,
            category,
            message,
        ));
    }

    /// Convenience helper: log a message at `Warning` level.
    pub fn warn(&mut self, sim_time: f64, category: impl Into<String>, message: impl Into<String>) {
        self.log(LogEvent::new(
            sim_time,
            0,
            EventLevel::Warning,
            category,
            message,
        ));
    }

    /// Convenience helper: log a message at `Error` level.
    pub fn error(
        &mut self,
        sim_time: f64,
        category: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.log(LogEvent::new(
            sim_time,
            0,
            EventLevel::Error,
            category,
            message,
        ));
    }

    /// Return all events whose severity is at least `min_level`.
    pub fn filter_level(&self, min_level: EventLevel) -> Vec<&LogEvent> {
        self.events
            .iter()
            .filter(|e| e.level >= min_level)
            .collect()
    }

    /// Return all events for the given category.
    pub fn filter_category<'a>(&'a self, cat: &str) -> Vec<&'a LogEvent> {
        self.events.iter().filter(|e| e.category == cat).collect()
    }

    /// Return all events whose simulation time falls in `[t_lo, t_hi]`.
    pub fn filter_sim_time(&self, t_lo: f64, t_hi: f64) -> Vec<&LogEvent> {
        self.events
            .iter()
            .filter(|e| e.sim_time >= t_lo && e.sim_time <= t_hi)
            .collect()
    }

    /// Serialize the log to a CSV string.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("sim_time,wall_time,level,category,message\n");
        for e in &self.events {
            out.push_str(&format!(
                "{},{},{},{},{}\n",
                e.sim_time,
                e.wall_time,
                e.level.as_str(),
                e.category,
                e.message
            ));
        }
        out
    }

    /// Total number of events recorded.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    /// Return `true` if no events have been logged.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
}

// ---------------------------------------------------------------------------
// ResultAggregator — statistical aggregation over snapshot collections
// ---------------------------------------------------------------------------

/// Summary statistics computed over a sequence of scalar values.
#[derive(Debug, Clone)]
pub struct AggStats {
    /// Number of values in the sample.
    pub count: usize,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Population standard deviation.
    pub std: f64,
    /// Sum of all values.
    pub sum: f64,
}

impl AggStats {
    /// Compute statistics from a slice.  Returns `None` if the slice is empty.
    pub fn from_slice(data: &[f64]) -> Option<Self> {
        if data.is_empty() {
            return None;
        }
        let n = data.len();
        let sum = data.iter().sum::<f64>();
        let mean = sum / n as f64;
        let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let variance = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        let std = variance.sqrt();
        Some(Self {
            count: n,
            min,
            max,
            mean,
            std,
            sum,
        })
    }
}

/// Aggregates statistics over a collection of snapshots for a named field.
#[derive(Debug, Default)]
pub struct ResultAggregator {
    /// Collected scalar samples (one per snapshot time step).
    pub samples: Vec<f64>,
}

impl ResultAggregator {
    /// Create an empty aggregator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add the mean value of `field` from `snapshot` to the sample set.
    pub fn add_snapshot_mean(&mut self, snapshot: &Snapshot, field: &str) {
        if let Some(data) = snapshot.fields.get(field)
            && !data.is_empty()
        {
            let mean = data.iter().sum::<f64>() / data.len() as f64;
            self.samples.push(mean);
        }
    }

    /// Add the maximum absolute value of `field` from `snapshot`.
    pub fn add_snapshot_max_abs(&mut self, snapshot: &Snapshot, field: &str) {
        if let Some(data) = snapshot.fields.get(field)
            && let Some(max_abs) = data.iter().cloned().map(f64::abs).reduce(f64::max)
        {
            self.samples.push(max_abs);
        }
    }

    /// Push a raw scalar sample directly.
    pub fn push(&mut self, v: f64) {
        self.samples.push(v);
    }

    /// Compute aggregate statistics over all pushed samples.
    pub fn compute(&self) -> Option<AggStats> {
        AggStats::from_slice(&self.samples)
    }

    /// Reset all samples.
    pub fn reset(&mut self) {
        self.samples.clear();
    }
}

// ---------------------------------------------------------------------------
// MetadataStore — simulation parameter and provenance metadata
// ---------------------------------------------------------------------------

/// A structured store for simulation-level metadata including git hash,
/// build timestamp, and runtime parameters.
#[derive(Debug, Default, Clone)]
pub struct MetadataStore {
    /// Key-value string pairs (e.g. `"git_hash"`, `"compiler_version"`).
    pub entries: HashMap<String, String>,
}

impl MetadataStore {
    /// Create an empty metadata store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set or overwrite an entry.
    pub fn set(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.entries.insert(key.into(), value.into());
    }

    /// Retrieve an entry by key.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(String::as_str)
    }

    /// Serialize to a newline-delimited `KEY=VALUE` format.
    pub fn to_properties(&self) -> String {
        let mut lines: Vec<String> = self
            .entries
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect();
        lines.sort();
        lines.join("\n")
    }

    /// Parse from a newline-delimited `KEY=VALUE` format.
    pub fn from_properties(s: &str) -> Self {
        let mut store = Self::new();
        for line in s.lines() {
            if let Some(pos) = line.find('=') {
                let k = &line[..pos];
                let v = &line[pos + 1..];
                store.set(k, v);
            }
        }
        store
    }

    /// Merge another store into `self`, overwriting on conflict.
    pub fn merge(&mut self, other: &MetadataStore) {
        for (k, v) in &other.entries {
            self.entries.insert(k.clone(), v.clone());
        }
    }

    /// Return the number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return `true` if the store is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

// ---------------------------------------------------------------------------
// ProvenanceType
// ---------------------------------------------------------------------------

/// Discriminant for provenance graph nodes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProvenanceType {
    /// A computational activity (e.g. a solver step).
    Computation,
    /// A data artifact (e.g. an input file or result dataset).
    Dataset,
    /// A software or human agent responsible for an action.
    Agent,
}

// ---------------------------------------------------------------------------
// ProvenanceNode
// ---------------------------------------------------------------------------

/// A node in the provenance graph.
#[derive(Debug, Clone)]
pub struct ProvenanceNode {
    /// Unique node identifier.
    pub id: String,
    /// Category of this node.
    pub type_: ProvenanceType,
    /// IDs of nodes that are inputs to this node.
    pub inputs: Vec<String>,
    /// IDs of nodes that are outputs of this node.
    pub outputs: Vec<String>,
}

impl ProvenanceNode {
    /// Create a leaf node with no inputs or outputs.
    pub fn new(id: impl Into<String>, type_: ProvenanceType) -> Self {
        Self {
            id: id.into(),
            type_,
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// ProvenanceTracker
// ---------------------------------------------------------------------------

/// Directed acyclic graph for W3C-PROV-style provenance tracking.
#[derive(Debug, Default)]
pub struct ProvenanceTracker {
    /// Ordered list of nodes in the provenance graph.
    pub graph: Vec<ProvenanceNode>,
}

impl ProvenanceTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a computation node that consumes `inputs` and produces `outputs`.
    pub fn add_computation(
        &mut self,
        id: impl Into<String>,
        inputs: Vec<String>,
        outputs: Vec<String>,
    ) {
        let mut node = ProvenanceNode::new(id, ProvenanceType::Computation);
        node.inputs = inputs;
        node.outputs = outputs;
        self.graph.push(node);
    }

    /// Record a data-source (dataset) node.
    pub fn add_data_source(&mut self, id: impl Into<String>, outputs: Vec<String>) {
        let mut node = ProvenanceNode::new(id, ProvenanceType::Dataset);
        node.outputs = outputs;
        self.graph.push(node);
    }

    /// Trace the full ancestry of `target_id`.
    ///
    /// Returns all node IDs that transitively contribute to `target_id`.
    pub fn trace_lineage(&self, target_id: &str) -> Vec<String> {
        let mut visited: Vec<String> = Vec::new();
        let mut queue: Vec<String> = vec![target_id.to_string()];
        while let Some(current) = queue.pop() {
            if visited.contains(&current) {
                continue;
            }
            visited.push(current.clone());
            for node in &self.graph {
                if node.outputs.contains(&current) && !visited.contains(&node.id) {
                    queue.push(node.id.clone());
                }
            }
            if let Some(node) = self.graph.iter().find(|n| n.id == current) {
                for inp in &node.inputs {
                    if !visited.contains(inp) {
                        queue.push(inp.clone());
                    }
                }
            }
        }
        visited
    }

    /// Serialise the provenance graph to a simplified PROV-JSON string.
    pub fn to_prov_json(&self) -> String {
        let mut out = String::from("{\n  \"nodes\": [\n");
        for (i, node) in self.graph.iter().enumerate() {
            let type_str = match node.type_ {
                ProvenanceType::Computation => "Activity",
                ProvenanceType::Dataset => "Entity",
                ProvenanceType::Agent => "Agent",
            };
            out.push_str(&format!(
                "    {{\"id\": \"{}\", \"type\": \"{}\", \"inputs\": [{}], \"outputs\": [{}]}}",
                node.id,
                type_str,
                node.inputs
                    .iter()
                    .map(|s| format!("\"{s}\""))
                    .collect::<Vec<_>>()
                    .join(", "),
                node.outputs
                    .iter()
                    .map(|s| format!("\"{s}\""))
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
            if i + 1 < self.graph.len() {
                out.push_str(",\n");
            } else {
                out.push('\n');
            }
        }
        out.push_str("  ]\n}");
        out
    }
}

// ---------------------------------------------------------------------------
// CheckpointManager
// ---------------------------------------------------------------------------

/// Manages a rolling window of simulation checkpoints stored in memory.
///
/// In a real deployment `base_dir` would name a directory on disk; here the
/// data is kept in an internal `HashMap` so the module has no I/O dependency.
#[derive(Debug)]
pub struct CheckpointManager {
    /// Base directory for checkpoint files (informational, not used for in-memory store).
    pub base_dir: String,
    /// Maximum number of checkpoints to retain.
    pub max_checkpoints: usize,
    store: HashMap<usize, Vec<f64>>,
}

impl CheckpointManager {
    /// Create a new manager for `base_dir` retaining at most `max_checkpoints`.
    pub fn new(base_dir: impl Into<String>, max_checkpoints: usize) -> Self {
        Self {
            base_dir: base_dir.into(),
            max_checkpoints: max_checkpoints.max(1),
            store: HashMap::new(),
        }
    }

    /// Save a checkpoint for `step`, returning the logical path string.
    pub fn save_checkpoint(&mut self, step: usize, data: &[f64]) -> String {
        self.store.insert(step, data.to_vec());
        self.cleanup_old();
        format!("{}/checkpoint_{step:06}.bin", self.base_dir)
    }

    /// Return all checkpoint step numbers in ascending order.
    pub fn list_checkpoints(&self) -> Vec<usize> {
        let mut steps: Vec<usize> = self.store.keys().copied().collect();
        steps.sort_unstable();
        steps
    }

    /// Load the data for `step`, or `None` if no such checkpoint exists.
    pub fn load_checkpoint(&self, step: usize) -> Option<Vec<f64>> {
        self.store.get(&step).cloned()
    }

    /// Remove old checkpoints so that at most `max_checkpoints` are kept.
    pub fn cleanup_old(&mut self) {
        let mut steps = self.list_checkpoints();
        while steps.len() > self.max_checkpoints {
            let oldest = steps.remove(0);
            self.store.remove(&oldest);
        }
    }
}

// ---------------------------------------------------------------------------
// ParameterSweep
// ---------------------------------------------------------------------------

/// A multi-dimensional parameter sweep specification.
///
/// Each entry is a `(name, values)` pair.  The sweep can be evaluated either
/// as a full Cartesian product or via a Latin-hypercube sample.
#[derive(Debug, Default)]
pub struct ParameterSweep {
    /// Ordered list of `(parameter_name, candidate_values)` pairs.
    pub params: Vec<(String, Vec<f64>)>,
}

impl ParameterSweep {
    /// Create an empty sweep.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a parameter with its list of candidate values.
    pub fn add_param(&mut self, name: impl Into<String>, values: Vec<f64>) {
        self.params.push((name.into(), values));
    }

    /// Total number of points in the Cartesian product.
    pub fn count(&self) -> usize {
        if self.params.is_empty() {
            return 0;
        }
        self.params.iter().map(|(_, v)| v.len()).product()
    }

    /// Enumerate every combination of parameter values (Cartesian product).
    pub fn cartesian_product(&self) -> Vec<HashMap<String, f64>> {
        if self.params.is_empty() {
            return vec![];
        }
        let mut result: Vec<HashMap<String, f64>> = vec![HashMap::new()];
        for (name, values) in &self.params {
            let mut next: Vec<HashMap<String, f64>> = Vec::new();
            for existing in &result {
                for &v in values {
                    let mut map = existing.clone();
                    map.insert(name.clone(), v);
                    next.push(map);
                }
            }
            result = next;
        }
        result
    }

    /// Draw `n` samples using a deterministic Latin-hypercube strategy.
    pub fn latin_hypercube_sample(&self, n: usize) -> Vec<HashMap<String, f64>> {
        if n == 0 || self.params.is_empty() {
            return vec![];
        }
        let mut samples: Vec<HashMap<String, f64>> = (0..n).map(|_| HashMap::new()).collect();
        for (dim, (name, values)) in self.params.iter().enumerate() {
            let k = values.len();
            if k == 0 {
                continue;
            }
            let mut perm: Vec<usize> = (0..n).collect();
            let seed: usize = dim
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            for i in (1..n).rev() {
                let j = seed.wrapping_mul(i).wrapping_add(dim) % (i + 1);
                perm.swap(i, j);
            }
            for (i, map) in samples.iter_mut().enumerate() {
                let idx = ((perm[i] * k) / n).min(k - 1);
                let t = (perm[i] as f64 + 0.5) / n as f64;
                let lo = values[idx];
                let hi = if idx + 1 < k { values[idx + 1] } else { lo };
                let local_t = (t * n as f64 - perm[i] as f64).clamp(0.0, 1.0);
                let v = lo + local_t * (hi - lo);
                map.insert(name.clone(), v);
            }
        }
        samples
    }
}

// ---------------------------------------------------------------------------
// QueryBuilder — fluent interface for filtering records
// ---------------------------------------------------------------------------

/// A builder for composing multi-criteria queries over a
/// [`SimulationDatabase`].
#[derive(Debug, Default)]
pub struct QueryBuilder<'a> {
    db: Option<&'a SimulationDatabase>,
    param_filters: Vec<(String, f64, f64)>,
    time_range: Option<(u64, u64)>,
    meta_filter: Option<(String, String)>,
}

impl<'a> QueryBuilder<'a> {
    /// Attach a database to query.
    pub fn from(db: &'a SimulationDatabase) -> Self {
        Self {
            db: Some(db),
            ..Self::default()
        }
    }

    /// Filter by a numeric parameter range `[lo, hi]`.
    pub fn param_range(mut self, name: impl Into<String>, lo: f64, hi: f64) -> Self {
        self.param_filters.push((name.into(), lo, hi));
        self
    }

    /// Filter by wall-clock timestamp range `[t_lo, t_hi]`.
    pub fn time_range(mut self, t_lo: u64, t_hi: u64) -> Self {
        self.time_range = Some((t_lo, t_hi));
        self
    }

    /// Filter by a required metadata key-value pair.
    pub fn meta_eq(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.meta_filter = Some((key.into(), value.into()));
        self
    }

    /// Execute the query and return matching records.
    pub fn execute(&self) -> Vec<&SimulationRecord> {
        let db = match self.db {
            Some(d) => d,
            None => return vec![],
        };
        db.records
            .iter()
            .filter(|r| {
                // Parameter range filters
                for (name, lo, hi) in &self.param_filters {
                    match r.parameters.get(name.as_str()) {
                        Some(&v) if v >= *lo && v <= *hi => {}
                        _ => return false,
                    }
                }
                // Time range filter
                if let Some((t_lo, t_hi)) = self.time_range
                    && (r.timestamp < t_lo || r.timestamp > t_hi)
                {
                    return false;
                }
                // Metadata filter
                if let Some((k, v)) = &self.meta_filter {
                    match r.metadata.get(k.as_str()) {
                        Some(val) if val == v => {}
                        _ => return false,
                    }
                }
                true
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- SimulationRecord ---

    #[test]
    fn test_record_new() {
        let r = SimulationRecord::new("run-001", 1_700_000_000);
        assert_eq!(r.id, "run-001");
        assert_eq!(r.timestamp, 1_700_000_000);
        assert!(r.parameters.is_empty());
        assert!(r.metadata.is_empty());
    }

    #[test]
    fn test_record_set_param() {
        let mut r = SimulationRecord::new("r1", 0);
        r.set_param("dt", 0.01);
        assert_eq!(r.parameters["dt"], 0.01);
    }

    #[test]
    fn test_record_set_meta() {
        let mut r = SimulationRecord::new("r2", 0);
        r.set_meta("solver", "rk4");
        assert_eq!(r.metadata["solver"], "rk4");
    }

    // --- SimulationDatabase ---

    #[test]
    fn test_db_new_empty() {
        let path = std::env::temp_dir().join("test_simdb.csv");
        let path_str = path.to_str().unwrap_or("").to_string();
        let db = SimulationDatabase::new(&path_str);
        assert!(db.records.is_empty());
        assert_eq!(db.file_path, path_str);
    }

    #[test]
    fn test_db_add_and_find() {
        let path = std::env::temp_dir().join("test_simdb_find.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        let r = SimulationRecord::new("abc", 42);
        db.add_record(r);
        assert!(db.find_by_id("abc").is_some());
        assert!(db.find_by_id("xyz").is_none());
    }

    #[test]
    fn test_db_delete_by_id() {
        let path = std::env::temp_dir().join("test_simdb_del.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        db.add_record(SimulationRecord::new("to_delete", 0));
        db.add_record(SimulationRecord::new("keep", 0));
        let removed = db.delete_by_id("to_delete");
        assert_eq!(removed, 1);
        assert!(db.find_by_id("to_delete").is_none());
        assert!(db.find_by_id("keep").is_some());
    }

    #[test]
    fn test_db_query_range_basic() {
        let path = std::env::temp_dir().join("test_simdb_qr.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        let mut r1 = SimulationRecord::new("a", 0);
        r1.set_param("dt", 0.01);
        let mut r2 = SimulationRecord::new("b", 0);
        r2.set_param("dt", 0.1);
        let mut r3 = SimulationRecord::new("c", 0);
        r3.set_param("dt", 1.0);
        db.add_record(r1);
        db.add_record(r2);
        db.add_record(r3);
        let found = db.query_range("dt", 0.005, 0.2);
        assert_eq!(found.len(), 2);
    }

    #[test]
    fn test_db_query_time_range() {
        let path = std::env::temp_dir().join("test_simdb_tr.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        for i in 0_u64..5 {
            db.add_record(SimulationRecord::new(format!("r{i}"), 1000 + i * 100));
        }
        let found = db.query_time_range(1100, 1300);
        assert_eq!(found.len(), 3); // ts 1100, 1200, 1300
    }

    #[test]
    fn test_db_save_and_load_roundtrip() {
        let path = std::env::temp_dir().join("test_simdb_rt.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        let mut r = SimulationRecord::new("run42", 999);
        r.set_param("Re", 1000.0);
        r.set_meta("solver", "rk4");
        db.add_record(r);
        let csv = db.save_to_csv();
        let mut db2 = SimulationDatabase::new(path.to_str().unwrap_or(""));
        db2.load_from_csv(&csv);
        assert_eq!(db2.records.len(), 1);
        assert_eq!(db2.records[0].id, "run42");
        assert_eq!(db2.records[0].timestamp, 999);
        assert!((db2.records[0].parameters["Re"] - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn test_db_statistics_basic() {
        let path = std::env::temp_dir().join("test_simdb_stats.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        for (i, v) in [1.0_f64, 2.0, 3.0, 4.0, 5.0].iter().enumerate() {
            let mut r = SimulationRecord::new(format!("r{i}"), 0);
            r.set_param("x", *v);
            db.add_record(r);
        }
        let (min, max, mean) = db.statistics("x");
        assert!((min - 1.0).abs() < 1e-9);
        assert!((max - 5.0).abs() < 1e-9);
        assert!((mean - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_db_export_json_contains_id() {
        let path = std::env::temp_dir().join("test_simdb_json.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        db.add_record(SimulationRecord::new("sim-1", 0));
        let json = db.export_json();
        assert!(json.contains("\"sim-1\""));
    }

    // --- RLE helpers ---

    #[test]
    fn test_rle_encode_all_same() {
        let data = vec![3.125; 100];
        let enc = rle_encode(&data);
        assert_eq!(enc.len(), 1);
        assert_eq!(enc[0].1, 100);
    }

    #[test]
    fn test_rle_roundtrip() {
        let data = vec![1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0];
        let enc = rle_encode(&data);
        let dec = rle_decode(&enc);
        assert_eq!(dec, data);
    }

    #[test]
    fn test_rle_empty() {
        let enc = rle_encode(&[]);
        assert!(enc.is_empty());
        let dec = rle_decode(&[]);
        assert!(dec.is_empty());
    }

    #[test]
    fn test_rle_compression_ratio() {
        let data = vec![0.0_f64; 50];
        let enc = rle_encode(&data);
        let ratio = rle_compression_ratio(data.len(), &enc);
        assert!(ratio > 1.0);
    }

    #[test]
    fn test_rle_no_compression_all_different() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let enc = rle_encode(&data);
        let ratio = rle_compression_ratio(data.len(), &enc);
        // Each run has length 1, so ratio == 1
        assert!((ratio - 1.0).abs() < 1e-9);
    }

    // --- SnapshotTable ---

    #[test]
    fn test_snapshot_insert_sorted() {
        let mut table = SnapshotTable::new();
        table.insert(Snapshot::new(3.0));
        table.insert(Snapshot::new(1.0));
        table.insert(Snapshot::new(2.0));
        let times: Vec<f64> = table.snapshots.iter().map(|s| s.time).collect();
        assert_eq!(times, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_snapshot_query_time_range() {
        let mut table = SnapshotTable::new();
        for t in [0.0, 1.0, 2.0, 3.0, 4.0] {
            table.insert(Snapshot::new(t));
        }
        let found = table.query_time_range(1.0, 3.0);
        assert_eq!(found.len(), 3);
    }

    #[test]
    fn test_snapshot_nearest() {
        let mut table = SnapshotTable::new();
        for t in [0.0, 1.0, 2.0] {
            table.insert(Snapshot::new(t));
        }
        let near = table.nearest(1.4).unwrap();
        assert!((near.time - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_snapshot_compress_decompress() {
        let mut table = SnapshotTable::new();
        let mut snap = Snapshot::new(0.0);
        snap.set_field("pressure", vec![1.0, 1.0, 1.0, 2.0]);
        table.insert(snap);
        table.compress_field(0, "pressure");
        let dec = table.decompress_field(0, "pressure").unwrap();
        assert_eq!(dec, vec![1.0, 1.0, 1.0, 2.0]);
    }

    #[test]
    fn test_snapshot_diff_changes_counted() {
        let mut table = SnapshotTable::new();
        let mut s0 = Snapshot::new(0.0);
        s0.set_field("u", vec![1.0, 2.0, 3.0]);
        let mut s1 = Snapshot::new(1.0);
        s1.set_field("u", vec![1.0, 9.0, 3.0]); // index 1 changed
        table.insert(s0);
        table.insert(s1);
        let changed = table.compute_diff("u", 1);
        assert_eq!(changed, 1);
        assert_eq!(table.diffs[0].indices, vec![1]);
        assert!((table.diffs[0].new_values[0] - 9.0).abs() < 1e-9);
    }

    #[test]
    fn test_snapshot_apply_diff() {
        let diff = DiffEntry {
            time: 1.0,
            field: "u".to_string(),
            indices: vec![0, 2],
            new_values: vec![10.0, 30.0],
        };
        let mut base = vec![1.0, 2.0, 3.0];
        SnapshotTable::apply_diff(&mut base, &diff);
        assert!((base[0] - 10.0).abs() < 1e-9);
        assert!((base[1] - 2.0).abs() < 1e-9);
        assert!((base[2] - 30.0).abs() < 1e-9);
    }

    // --- EventLog ---

    #[test]
    fn test_event_log_basic() {
        let mut log = EventLog::new();
        log.info(1.0, "solver", "step started");
        assert_eq!(log.len(), 1);
        assert!(!log.is_empty());
    }

    #[test]
    fn test_event_log_filter_level() {
        let mut log = EventLog::new();
        log.info(0.0, "a", "msg");
        log.warn(1.0, "b", "warn");
        log.error(2.0, "c", "err");
        let warns_and_above = log.filter_level(EventLevel::Warning);
        assert_eq!(warns_and_above.len(), 2);
    }

    #[test]
    fn test_event_log_filter_category() {
        let mut log = EventLog::new();
        log.info(0.0, "io", "read file");
        log.info(0.5, "solver", "step 1");
        log.info(1.0, "io", "write file");
        let io_events = log.filter_category("io");
        assert_eq!(io_events.len(), 2);
    }

    #[test]
    fn test_event_log_filter_sim_time() {
        let mut log = EventLog::new();
        for t in [0.0, 1.0, 2.0, 3.0, 4.0_f64] {
            log.info(t, "x", "msg");
        }
        let found = log.filter_sim_time(1.0, 3.0);
        assert_eq!(found.len(), 3);
    }

    #[test]
    fn test_event_log_to_csv() {
        let mut log = EventLog::new();
        log.info(1.0, "solver", "done");
        let csv = log.to_csv();
        assert!(csv.contains("INFO"));
        assert!(csv.contains("solver"));
        assert!(csv.contains("done"));
    }

    #[test]
    fn test_event_level_ordering() {
        assert!(EventLevel::Critical > EventLevel::Error);
        assert!(EventLevel::Error > EventLevel::Warning);
        assert!(EventLevel::Warning > EventLevel::Info);
        assert!(EventLevel::Info > EventLevel::Debug);
    }

    // --- ResultAggregator ---

    #[test]
    fn test_agg_stats_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = AggStats::from_slice(&data).unwrap();
        assert!((stats.min - 1.0).abs() < 1e-9);
        assert!((stats.max - 5.0).abs() < 1e-9);
        assert!((stats.mean - 3.0).abs() < 1e-9);
        // std of [1,2,3,4,5]: variance = 2.0, std = sqrt(2)
        assert!((stats.std - 2_f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_agg_stats_empty() {
        assert!(AggStats::from_slice(&[]).is_none());
    }

    #[test]
    fn test_result_aggregator_push_compute() {
        let mut agg = ResultAggregator::new();
        agg.push(10.0);
        agg.push(20.0);
        agg.push(30.0);
        let stats = agg.compute().unwrap();
        assert!((stats.mean - 20.0).abs() < 1e-9);
    }

    #[test]
    fn test_result_aggregator_add_snapshot_mean() {
        let mut agg = ResultAggregator::new();
        let mut snap = Snapshot::new(0.0);
        snap.set_field("v", vec![2.0, 4.0, 6.0]);
        agg.add_snapshot_mean(&snap, "v");
        let stats = agg.compute().unwrap();
        assert!((stats.mean - 4.0).abs() < 1e-9);
    }

    #[test]
    fn test_result_aggregator_reset() {
        let mut agg = ResultAggregator::new();
        agg.push(1.0);
        agg.reset();
        assert!(agg.compute().is_none());
    }

    // --- MetadataStore ---

    #[test]
    fn test_metadata_store_set_get() {
        let mut store = MetadataStore::new();
        store.set("git_hash", "abc123");
        assert_eq!(store.get("git_hash"), Some("abc123"));
        assert_eq!(store.get("missing"), None);
    }

    #[test]
    fn test_metadata_store_roundtrip() {
        let mut store = MetadataStore::new();
        store.set("a", "1");
        store.set("b", "hello world");
        let props = store.to_properties();
        let loaded = MetadataStore::from_properties(&props);
        assert_eq!(loaded.get("a"), Some("1"));
        assert_eq!(loaded.get("b"), Some("hello world"));
    }

    #[test]
    fn test_metadata_store_merge() {
        let mut a = MetadataStore::new();
        a.set("x", "1");
        let mut b = MetadataStore::new();
        b.set("y", "2");
        b.set("x", "overwritten");
        a.merge(&b);
        assert_eq!(a.get("x"), Some("overwritten"));
        assert_eq!(a.get("y"), Some("2"));
    }

    // --- ProvenanceTracker ---

    #[test]
    fn test_prov_add_computation() {
        let mut tracker = ProvenanceTracker::new();
        tracker.add_computation("comp1", vec!["data_in".into()], vec!["data_out".into()]);
        assert_eq!(tracker.graph.len(), 1);
        assert_eq!(tracker.graph[0].type_, ProvenanceType::Computation);
    }

    #[test]
    fn test_prov_trace_lineage_chain() {
        let mut tracker = ProvenanceTracker::new();
        tracker.add_data_source("raw", vec!["raw".into()]);
        tracker.add_computation("step1", vec!["raw".into()], vec!["processed".into()]);
        tracker.add_computation("step2", vec!["processed".into()], vec!["result".into()]);
        let lineage = tracker.trace_lineage("result");
        assert!(lineage.contains(&"result".to_string()));
        assert!(lineage.contains(&"step2".to_string()));
    }

    // --- CheckpointManager ---

    #[test]
    fn test_checkpoint_save_and_load() {
        let path = std::env::temp_dir().join("test_simdb_checkpoints");
        let mut mgr = CheckpointManager::new(path.to_str().unwrap_or(""), 5);
        mgr.save_checkpoint(0, &[1.0, 2.0, 3.0]);
        let data = mgr.load_checkpoint(0).unwrap();
        assert_eq!(data, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_checkpoint_cleanup_old() {
        let path = std::env::temp_dir().join("test_simdb_ckpt");
        let mut mgr = CheckpointManager::new(path.to_str().unwrap_or(""), 3);
        for step in 0..6_usize {
            mgr.save_checkpoint(step, &[step as f64]);
        }
        assert!(mgr.list_checkpoints().len() <= 3);
    }

    // --- ParameterSweep ---

    #[test]
    fn test_sweep_cartesian_product() {
        let mut sweep = ParameterSweep::new();
        sweep.add_param("dt", vec![0.01, 0.1]);
        sweep.add_param("Re", vec![100.0, 500.0, 1000.0]);
        assert_eq!(sweep.count(), 6);
        let product = sweep.cartesian_product();
        assert_eq!(product.len(), 6);
    }

    #[test]
    fn test_sweep_latin_hypercube() {
        let mut sweep = ParameterSweep::new();
        sweep.add_param("x", vec![0.0, 1.0, 2.0, 3.0]);
        sweep.add_param("y", vec![10.0, 20.0, 30.0]);
        let samples = sweep.latin_hypercube_sample(5);
        assert_eq!(samples.len(), 5);
    }

    // --- QueryBuilder ---

    #[test]
    fn test_query_builder_param_range() {
        let path = std::env::temp_dir().join("test_simdb_qbpr.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        for i in 0..5_usize {
            let mut r = SimulationRecord::new(format!("r{i}"), 0);
            r.set_param("v", i as f64);
            db.add_record(r);
        }
        let binding = QueryBuilder::from(&db).param_range("v", 1.0, 3.0);
        let results = binding.execute();
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_query_builder_time_range() {
        let path = std::env::temp_dir().join("test_simdb_qbtr.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        for i in 0_u64..5 {
            db.add_record(SimulationRecord::new(format!("t{i}"), 1000 + i * 100));
        }
        let binding = QueryBuilder::from(&db).time_range(1100, 1300);
        let results = binding.execute();
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_query_builder_meta_eq() {
        let path = std::env::temp_dir().join("test_simdb_qbme.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        let mut r1 = SimulationRecord::new("a", 0);
        r1.set_meta("solver", "rk4");
        let mut r2 = SimulationRecord::new("b", 0);
        r2.set_meta("solver", "euler");
        db.add_record(r1);
        db.add_record(r2);
        let binding = QueryBuilder::from(&db).meta_eq("solver", "rk4");
        let results = binding.execute();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].id, "a");
    }

    #[test]
    fn test_query_builder_combined_filters() {
        let path = std::env::temp_dir().join("test_simdb_qbcf.csv");
        let mut db = SimulationDatabase::new(path.to_str().unwrap_or(""));
        for i in 0_u64..6 {
            let mut r = SimulationRecord::new(format!("r{i}"), 1000 + i * 100);
            r.set_param("Re", (i as f64) * 100.0);
            r.set_meta("solver", if i % 2 == 0 { "rk4" } else { "euler" });
            db.add_record(r);
        }
        let binding = QueryBuilder::from(&db)
            .param_range("Re", 100.0, 400.0)
            .meta_eq("solver", "euler");
        let results = binding.execute();
        // Re in [100,400]: i=1(100),2(200),3(300),4(400)
        // euler: i=1,3,5 → intersection: i=1(100),3(300)
        assert_eq!(results.len(), 2);
    }
}
