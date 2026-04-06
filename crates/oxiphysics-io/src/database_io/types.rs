//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::collections::{HashMap, VecDeque};

/// A single row in the simulation database, mapping column names to values.
#[derive(Debug, Clone, Default)]
pub struct DbRow {
    /// Map from column name to value.
    pub columns: HashMap<String, DbValue>,
}
impl DbRow {
    /// Create an empty row.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a key-value pair into this row.
    pub fn set(&mut self, key: impl Into<String>, val: impl Into<DbValue>) {
        self.columns.insert(key.into(), val.into());
    }
    /// Get a value by column name.
    pub fn get(&self, key: &str) -> Option<&DbValue> {
        self.columns.get(key)
    }
}
/// Metadata for a single simulation run.
#[derive(Debug, Clone)]
pub struct SimulationMetadata {
    /// Unique simulation identifier.
    pub sim_id: String,
    /// Human-readable description.
    pub description: String,
    /// Creation timestamp (Unix seconds, as f64).
    pub created_at: f64,
    /// Arbitrary key-value parameters.
    pub parameters: HashMap<String, String>,
    /// List of artifact file paths associated with this simulation.
    pub artifacts: Vec<String>,
}
impl SimulationMetadata {
    /// Create a new metadata entry with empty parameters and artifacts.
    pub fn new(sim_id: impl Into<String>, description: impl Into<String>, created_at: f64) -> Self {
        Self {
            sim_id: sim_id.into(),
            description: description.into(),
            created_at,
            parameters: HashMap::new(),
            artifacts: Vec::new(),
        }
    }
    /// Add a simulation parameter.
    pub fn set_param(&mut self, key: impl Into<String>, val: impl Into<String>) {
        self.parameters.insert(key.into(), val.into());
    }
    /// Add an artifact path.
    pub fn add_artifact(&mut self, path: impl Into<String>) {
        self.artifacts.push(path.into());
    }
}
/// Registry of simulation metadata entries (data catalog / provenance store).
#[derive(Debug, Default)]
pub struct DataCatalog {
    pub(super) entries: HashMap<String, SimulationMetadata>,
}
impl DataCatalog {
    /// Create an empty data catalog.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a simulation entry.
    pub fn register(&mut self, meta: SimulationMetadata) {
        self.entries.insert(meta.sim_id.clone(), meta);
    }
    /// Look up a simulation by ID.
    pub fn lookup(&self, sim_id: &str) -> Option<&SimulationMetadata> {
        self.entries.get(sim_id)
    }
    /// Look up a simulation by ID (mutable).
    pub fn lookup_mut(&mut self, sim_id: &str) -> Option<&mut SimulationMetadata> {
        self.entries.get_mut(sim_id)
    }
    /// Remove a simulation by ID.
    ///
    /// Returns `true` if an entry was removed.
    pub fn remove(&mut self, sim_id: &str) -> bool {
        self.entries.remove(sim_id).is_some()
    }
    /// Search for simulations whose description contains `query` (case-insensitive).
    pub fn search_description(&self, query: &str) -> Vec<&SimulationMetadata> {
        let q = query.to_lowercase();
        self.entries
            .values()
            .filter(|m| m.description.to_lowercase().contains(&q))
            .collect()
    }
    /// Return all simulation IDs.
    pub fn all_ids(&self) -> Vec<&str> {
        self.entries.keys().map(|s| s.as_str()).collect()
    }
    /// Number of registered simulations.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` if no simulations are registered.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// Append-only time series with range query and decimation.
#[derive(Debug, Default)]
pub struct TimeSeriesStore {
    /// Internal storage (always sorted by time).
    pub(super) samples: Vec<TsSample>,
    /// Label/name of this time series.
    pub label: String,
}
impl TimeSeriesStore {
    /// Create a new, empty time series with the given label.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            samples: Vec::new(),
            label: label.into(),
        }
    }
    /// Append a sample.  Samples should be added in ascending time order.
    pub fn append(&mut self, time: f64, value: f64) {
        self.samples.push(TsSample::new(time, value));
    }
    /// Number of samples stored.
    pub fn len(&self) -> usize {
        self.samples.len()
    }
    /// Returns `true` if no samples have been recorded.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
    /// Return all samples whose timestamp is in `[t_start, t_end]`.
    pub fn range_query(&self, t_start: f64, t_end: f64) -> Vec<&TsSample> {
        self.samples
            .iter()
            .filter(|s| s.time >= t_start && s.time <= t_end)
            .collect()
    }
    /// Decimate the series by keeping every `n`-th sample.
    ///
    /// Returns a new `TimeSeriesStore` with the decimated data.
    pub fn decimate(&self, n: usize) -> TimeSeriesStore {
        let mut out = TimeSeriesStore::new(format!("{}_dec{n}", self.label));
        if n == 0 {
            return out;
        }
        for (i, s) in self.samples.iter().enumerate() {
            if i % n == 0 {
                out.append(s.time, s.value);
            }
        }
        out
    }
    /// Compute the mean value over a time range `[t_start, t_end]`.
    pub fn mean_in_range(&self, t_start: f64, t_end: f64) -> Option<f64> {
        let vals: Vec<f64> = self
            .range_query(t_start, t_end)
            .iter()
            .map(|s| s.value)
            .collect();
        if vals.is_empty() {
            None
        } else {
            Some(vals.iter().sum::<f64>() / vals.len() as f64)
        }
    }
    /// Return the maximum value in the entire series.
    pub fn max_value(&self) -> Option<f64> {
        self.samples.iter().map(|s| s.value).reduce(f64::max)
    }
    /// Return the minimum value in the entire series.
    pub fn min_value(&self) -> Option<f64> {
        self.samples.iter().map(|s| s.value).reduce(f64::min)
    }
    /// Export the series to a CSV string `"time,value\n..."`.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("time,value\n");
        for s in &self.samples {
            out.push_str(&format!("{},{}\n", s.time, s.value));
        }
        out
    }
    /// Interpolate the value at time `t` using linear interpolation.
    ///
    /// Returns `None` if the series is empty or `t` is out of range.
    pub fn interpolate(&self, t: f64) -> Option<f64> {
        if self.samples.is_empty() {
            return None;
        }
        let pos = self.samples.partition_point(|s| s.time <= t);
        if pos == 0 {
            return Some(self.samples[0].value);
        }
        if pos >= self.samples.len() {
            return Some(
                self.samples
                    .last()
                    .expect("collection should not be empty")
                    .value,
            );
        }
        let lo = &self.samples[pos - 1];
        let hi = &self.samples[pos];
        let dt = hi.time - lo.time;
        if dt < 1e-15 {
            return Some(lo.value);
        }
        let frac = (t - lo.time) / dt;
        Some(lo.value + frac * (hi.value - lo.value))
    }
}
/// Catalog of simulation snapshots supporting lazy-loading semantics.
///
/// Snapshots are stored by ID. The catalog provides range queries by
/// simulation time and tag-based filtering.
#[derive(Debug, Default)]
pub struct SnapshotCatalog {
    pub(super) snapshots: HashMap<String, SnapshotEntry>,
    /// Auto-increment counter for generated IDs.
    pub(super) next_id: u64,
}
impl SnapshotCatalog {
    /// Create an empty catalog.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a snapshot entry and return its id.
    pub fn register(&mut self, entry: SnapshotEntry) -> &str {
        let id = entry.id.clone();
        self.snapshots.insert(id.clone(), entry);
        self.snapshots[&id].id.as_str()
    }
    /// Auto-generate an id and register a snapshot.
    ///
    /// Returns the generated id as an owned `String`.
    pub fn register_auto(&mut self, sim_time: f64, path: impl Into<String>) -> String {
        let id = format!("snap_{:06}", self.next_id);
        self.next_id += 1;
        self.register(SnapshotEntry::new(id.clone(), sim_time, path));
        id
    }
    /// Look up a snapshot by id.
    pub fn get(&self, id: &str) -> Option<&SnapshotEntry> {
        self.snapshots.get(id)
    }
    /// Look up a snapshot mutably by id.
    pub fn get_mut(&mut self, id: &str) -> Option<&mut SnapshotEntry> {
        self.snapshots.get_mut(id)
    }
    /// Remove a snapshot by id.  Returns `true` if removed.
    pub fn remove(&mut self, id: &str) -> bool {
        self.snapshots.remove(id).is_some()
    }
    /// Number of snapshots in the catalog.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }
    /// Returns `true` if the catalog is empty.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }
    /// Return all snapshots whose `sim_time` is in `[t_start, t_end]`.
    pub fn range_query(&self, t_start: f64, t_end: f64) -> Vec<&SnapshotEntry> {
        self.snapshots
            .values()
            .filter(|s| s.sim_time >= t_start && s.sim_time <= t_end)
            .collect()
    }
    /// Return all snapshots that contain `tag`.
    pub fn query_by_tag(&self, tag: &str) -> Vec<&SnapshotEntry> {
        self.snapshots
            .values()
            .filter(|s| s.tags.iter().any(|t| t == tag))
            .collect()
    }
    /// Return ids of all unloaded snapshots.
    pub fn unloaded_ids(&self) -> Vec<&str> {
        self.snapshots
            .values()
            .filter(|s| !s.loaded)
            .map(|s| s.id.as_str())
            .collect()
    }
    /// Mark all snapshots in a time range as loaded.
    pub fn mark_range_loaded(&mut self, t_start: f64, t_end: f64) {
        for s in self.snapshots.values_mut() {
            if s.sim_time >= t_start && s.sim_time <= t_end {
                s.loaded = true;
            }
        }
    }
    /// Compute total file size across all snapshots.
    pub fn total_file_size(&self) -> usize {
        self.snapshots.values().map(|s| s.file_size).sum()
    }
    /// Return snapshots sorted by sim_time (ascending).
    pub fn sorted_by_time(&self) -> Vec<&SnapshotEntry> {
        let mut v: Vec<&SnapshotEntry> = self.snapshots.values().collect();
        v.sort_by(|a, b| {
            a.sim_time
                .partial_cmp(&b.sim_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        v
    }
}
/// A cached simulation result identified by a string key.
#[derive(Debug, Clone)]
pub struct CacheEntry {
    /// Cache key.
    pub key: String,
    /// Serialised or raw result data.
    pub data: Vec<f64>,
    /// Simulation time at which this result was computed.
    pub sim_time: f64,
    /// Cache version tag for invalidation.
    pub version: u64,
}
impl CacheEntry {
    /// Create a new cache entry.
    pub fn new(key: impl Into<String>, data: Vec<f64>, sim_time: f64, version: u64) -> Self {
        Self {
            key: key.into(),
            data,
            sim_time,
            version,
        }
    }
}
/// A single simulation run record with metadata.
///
/// Stores the simulation name, Unix timestamp, and arbitrary string parameters.
#[derive(Debug, Clone)]
pub struct SimulationRecord {
    /// Unique run name / identifier.
    pub name: String,
    /// Unix timestamp (seconds since epoch) as f64.
    pub timestamp: f64,
    /// Arbitrary key-value simulation parameters.
    pub params: HashMap<String, String>,
    /// Optional output file path.
    pub output_path: Option<String>,
}
impl SimulationRecord {
    /// Create a new record with the given name and timestamp.
    pub fn new(name: impl Into<String>, timestamp: f64) -> Self {
        Self {
            name: name.into(),
            timestamp,
            params: HashMap::new(),
            output_path: None,
        }
    }
    /// Set a parameter value.
    pub fn set_param(&mut self, key: impl Into<String>, val: impl Into<String>) {
        self.params.insert(key.into(), val.into());
    }
    /// Get a parameter value by key.
    pub fn get_param(&self, key: &str) -> Option<&str> {
        self.params.get(key).map(|s| s.as_str())
    }
    /// Set the output file path.
    pub fn set_output(&mut self, path: impl Into<String>) {
        self.output_path = Some(path.into());
    }
}
/// Serialize/deserialize `SimulationRecord` to a simple JSON-like text format.
///
/// The format is:
/// ```text
/// {"name":"`n`","timestamp":`t`,"params":{"k1":"v1",...},"output":"`path`"}
/// ```
#[derive(Debug, Clone, Default)]
pub struct DatabaseSerializer;
impl DatabaseSerializer {
    /// Create a new serializer.
    pub fn new() -> Self {
        Self
    }
    /// Serialize a `SimulationRecord` to a JSON-like string.
    pub fn serialize(&self, rec: &SimulationRecord) -> String {
        let params_str: Vec<String> = rec
            .params
            .iter()
            .map(|(k, v)| format!(r#""{k}":"{v}""#))
            .collect();
        let params_json = format!("{{{}}}", params_str.join(","));
        let output_str = match &rec.output_path {
            Some(p) => format!(r#""{p}""#),
            None => "null".to_string(),
        };
        format!(
            r#"{{"name":"{name}","timestamp":{ts},"params":{params},"output":{out}}}"#,
            name = rec.name,
            ts = rec.timestamp,
            params = params_json,
            out = output_str,
        )
    }
    /// Serialize a list of records to a JSON array string.
    pub fn serialize_all(&self, recs: &[&SimulationRecord]) -> String {
        let items: Vec<String> = recs.iter().map(|r| self.serialize(r)).collect();
        format!("[{}]", items.join(","))
    }
    /// Deserialize a single record from a JSON-like string.
    ///
    /// Uses minimal parsing; fields must appear in the order produced by
    /// `serialize`. Returns `None` on parse failure.
    pub fn deserialize(&self, s: &str) -> Option<SimulationRecord> {
        let name = self.extract_string_field(s, "name")?;
        let ts_str = self.extract_raw_field(s, "timestamp")?;
        let timestamp: f64 = ts_str.trim().parse().ok()?;
        let mut rec = SimulationRecord::new(name, timestamp);
        let out_raw = self.extract_raw_field(s, "output").unwrap_or_default();
        let out_raw = out_raw.trim();
        if out_raw != "null" && out_raw.starts_with('"') {
            let cleaned = out_raw.trim_matches('"').to_string();
            rec.set_output(cleaned);
        }
        Some(rec)
    }
    fn extract_string_field(&self, s: &str, field: &str) -> Option<String> {
        let key = format!(r#""{field}":""#);
        let start = s.find(key.as_str())? + key.len();
        let rest = &s[start..];
        let end = rest.find('"')?;
        Some(rest[..end].to_string())
    }
    fn extract_raw_field(&self, s: &str, field: &str) -> Option<String> {
        let key = format!(r#""{field}":"#);
        let start = s.find(key.as_str())? + key.len();
        let rest = &s[start..];
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        Some(rest[..end].to_string())
    }
}
/// A single sample in a time series.
#[derive(Debug, Clone)]
pub struct TsSample {
    /// Timestamp in seconds.
    pub time: f64,
    /// Scalar value at this timestamp.
    pub value: f64,
}
impl TsSample {
    /// Create a new time series sample.
    pub fn new(time: f64, value: f64) -> Self {
        Self { time, value }
    }
}
/// A material entry in the material database.
#[derive(Debug, Clone)]
pub struct MaterialRecord {
    /// Material name/identifier.
    pub name: String,
    /// Mass density in kg/m³.
    pub density: f64,
    /// Young's modulus in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Thermal conductivity in W/(m·K).
    pub thermal_conductivity: f64,
    /// Yield strength in Pa (0 if not applicable).
    pub yield_strength: f64,
    /// Arbitrary tags for fuzzy search.
    pub tags: Vec<String>,
}
impl MaterialRecord {
    /// Create a new material record.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: impl Into<String>,
        density: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
        thermal_conductivity: f64,
        yield_strength: f64,
        tags: Vec<String>,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            youngs_modulus,
            poisson_ratio,
            thermal_conductivity,
            yield_strength,
            tags,
        }
    }
}
/// Column filter for export: keep only these columns (if empty, keep all).
#[derive(Debug, Clone, Default)]
pub struct ExportFilter {
    /// Columns to include (empty = all columns).
    pub columns: Vec<String>,
    /// Optional minimum value for a numeric column named `filter_col`.
    pub min_value: Option<(String, f64)>,
    /// Optional maximum value for a numeric column named `filter_col`.
    pub max_value: Option<(String, f64)>,
}
impl ExportFilter {
    /// Create an empty (pass-all) filter.
    pub fn new() -> Self {
        Self::default()
    }
    /// Check whether a row passes the numeric range filters.
    pub fn row_passes(&self, row: &DbRow) -> bool {
        if let Some((col, lo)) = &self.min_value {
            let v = row
                .get(col)
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::NEG_INFINITY);
            if v < *lo {
                return false;
            }
        }
        if let Some((col, hi)) = &self.max_value {
            let v = row
                .get(col)
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::INFINITY);
            if v > *hi {
                return false;
            }
        }
        true
    }
}
/// Export format selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    /// Comma-separated values.
    Csv,
    /// JSON array of objects.
    Json,
    /// Simplified HDF5-like text representation.
    Hdf5Text,
}
/// Export pipeline: converts `SimulationDatabase` rows to text output.
#[derive(Debug)]
pub struct ExportPipeline {
    /// Output format.
    pub format: ExportFormat,
    /// Row filter.
    pub filter: ExportFilter,
}
impl ExportPipeline {
    /// Create a new export pipeline with the given format and no filter.
    pub fn new(format: ExportFormat) -> Self {
        Self {
            format,
            filter: ExportFilter::new(),
        }
    }
    /// Set the export filter.
    pub fn with_filter(mut self, filter: ExportFilter) -> Self {
        self.filter = filter;
        self
    }
    /// Run the export, returning the serialised output as a `String`.
    pub fn export(&self, db: &SimulationDatabase) -> String {
        let cols = &self.filter.columns;
        let rows: Vec<&DbRow> = db
            .rows
            .iter()
            .filter(|r| self.filter.row_passes(r))
            .collect();
        match self.format {
            ExportFormat::Csv => self.to_csv(&rows, cols),
            ExportFormat::Json => self.to_json(&rows, cols),
            ExportFormat::Hdf5Text => self.to_hdf5_text(&rows, cols, &db.name),
        }
    }
    /// Determine effective column list from rows.
    fn effective_cols(&self, rows: &[&DbRow], cols: &[String]) -> Vec<String> {
        if cols.is_empty() {
            let mut seen: Vec<String> = Vec::new();
            for row in rows {
                for k in row.columns.keys() {
                    if !seen.contains(k) {
                        seen.push(k.clone());
                    }
                }
            }
            seen.sort();
            seen
        } else {
            cols.to_vec()
        }
    }
    fn value_to_str(v: &DbValue) -> String {
        match v {
            DbValue::Int(i) => i.to_string(),
            DbValue::Float(f) => format!("{f}"),
            DbValue::Text(s) => s.clone(),
            DbValue::Bool(b) => b.to_string(),
            DbValue::Null => "".to_string(),
        }
    }
    fn to_csv(&self, rows: &[&DbRow], cols: &[String]) -> String {
        let headers = self.effective_cols(rows, cols);
        let mut out = headers.join(",");
        out.push('\n');
        for row in rows {
            let line: Vec<String> = headers
                .iter()
                .map(|c| row.get(c).map(Self::value_to_str).unwrap_or_default())
                .collect();
            out.push_str(&line.join(","));
            out.push('\n');
        }
        out
    }
    fn to_json(&self, rows: &[&DbRow], cols: &[String]) -> String {
        let headers = self.effective_cols(rows, cols);
        let mut out = String::from("[\n");
        for (ri, row) in rows.iter().enumerate() {
            out.push_str("  {");
            let fields: Vec<String> = headers
                .iter()
                .map(|c| {
                    let val = row
                        .get(c)
                        .map(|v| match v {
                            DbValue::Text(s) => format!(r#""{s}""#),
                            DbValue::Null => "null".to_string(),
                            other => Self::value_to_str(other),
                        })
                        .unwrap_or_else(|| "null".to_string());
                    format!(r#""{c}":{val}"#)
                })
                .collect();
            out.push_str(&fields.join(","));
            out.push('}');
            if ri + 1 < rows.len() {
                out.push(',');
            }
            out.push('\n');
        }
        out.push(']');
        out
    }
    fn to_hdf5_text(&self, rows: &[&DbRow], cols: &[String], table_name: &str) -> String {
        let headers = self.effective_cols(rows, cols);
        let mut out = format!("# HDF5-like text dump of table '{table_name}'\n");
        out.push_str(&format!("# columns: {}\n", headers.join(",")));
        out.push_str(&format!("# rows: {}\n", rows.len()));
        for row in rows {
            let line: Vec<String> = headers
                .iter()
                .map(|c| row.get(c).map(Self::value_to_str).unwrap_or_default())
                .collect();
            out.push_str(&line.join(" "));
            out.push('\n');
        }
        out
    }
}
/// In-memory material property database with fuzzy search and interpolation.
#[derive(Debug, Default)]
pub struct MaterialDatabase {
    pub(super) records: Vec<MaterialRecord>,
}
impl MaterialDatabase {
    /// Create an empty material database.
    pub fn new() -> Self {
        Self::default()
    }
    /// Create a database pre-populated with common engineering materials.
    pub fn with_defaults() -> Self {
        let mut db = Self::new();
        db.insert(MaterialRecord::new(
            "steel_1020",
            7850.0,
            210e9,
            0.29,
            50.0,
            250e6,
            vec!["metal".into(), "steel".into(), "ferrous".into()],
        ));
        db.insert(MaterialRecord::new(
            "aluminium_6061",
            2700.0,
            69e9,
            0.33,
            167.0,
            276e6,
            vec!["metal".into(), "aluminium".into(), "light".into()],
        ));
        db.insert(MaterialRecord::new(
            "copper",
            8960.0,
            110e9,
            0.34,
            385.0,
            70e6,
            vec!["metal".into(), "copper".into(), "conductor".into()],
        ));
        db.insert(MaterialRecord::new(
            "polycarbonate",
            1200.0,
            2.4e9,
            0.37,
            0.2,
            60e6,
            vec!["polymer".into(), "plastic".into(), "transparent".into()],
        ));
        db.insert(MaterialRecord::new(
            "concrete",
            2300.0,
            30e9,
            0.20,
            1.7,
            3e6,
            vec!["composite".into(), "concrete".into(), "brittle".into()],
        ));
        db
    }
    /// Insert a material record.
    pub fn insert(&mut self, record: MaterialRecord) {
        self.records.push(record);
    }
    /// Exact lookup by name.
    pub fn lookup(&self, name: &str) -> Option<&MaterialRecord> {
        self.records.iter().find(|r| r.name == name)
    }
    /// Fuzzy search: return all records whose name or tags contain `query`
    /// (case-insensitive substring match).
    pub fn fuzzy_search(&self, query: &str) -> Vec<&MaterialRecord> {
        let q = query.to_lowercase();
        self.records
            .iter()
            .filter(|r| {
                r.name.to_lowercase().contains(&q)
                    || r.tags.iter().any(|t| t.to_lowercase().contains(&q))
            })
            .collect()
    }
    /// Interpolate material properties between two named materials.
    ///
    /// Returns a new `MaterialRecord` with blended properties at weight `t`
    /// (0 = pure first material, 1 = pure second material).
    ///
    /// Returns `None` if either material is not found.
    pub fn interpolate(&self, name_a: &str, name_b: &str, t: f64) -> Option<MaterialRecord> {
        let a = self.lookup(name_a)?;
        let b = self.lookup(name_b)?;
        let lerp = |va: f64, vb: f64| va + t * (vb - va);
        Some(MaterialRecord::new(
            format!("{name_a}_to_{name_b}_{t:.2}"),
            lerp(a.density, b.density),
            lerp(a.youngs_modulus, b.youngs_modulus),
            lerp(a.poisson_ratio, b.poisson_ratio),
            lerp(a.thermal_conductivity, b.thermal_conductivity),
            lerp(a.yield_strength, b.yield_strength),
            vec![],
        ))
    }
    /// Number of materials in the database.
    pub fn len(&self) -> usize {
        self.records.len()
    }
    /// Returns `true` if no materials are stored.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    /// Return all material names.
    pub fn names(&self) -> Vec<&str> {
        self.records.iter().map(|r| r.name.as_str()).collect()
    }
}
/// LRU cache for expensive simulation results.
///
/// Evicts least-recently-used entries when the capacity limit is reached.
/// Invalidation is by version number: entries with a stale version are evicted.
#[derive(Debug)]
pub struct ResultCache {
    /// Maximum number of entries.
    pub capacity: usize,
    /// Current cache version; entries below this are invalid.
    pub current_version: u64,
    /// LRU queue: front = most recently used.
    pub(super) entries: VecDeque<CacheEntry>,
}
impl ResultCache {
    /// Create a new LRU cache with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            current_version: 0,
            entries: VecDeque::new(),
        }
    }
    /// Insert or update a cache entry.  Evicts LRU entry when over capacity.
    pub fn put(&mut self, entry: CacheEntry) {
        self.entries.retain(|e| e.key != entry.key);
        self.entries.push_front(entry);
        while self.entries.len() > self.capacity {
            self.entries.pop_back();
        }
    }
    /// Retrieve a cached entry by key, moving it to the front (MRU).
    ///
    /// Returns `None` if the key is not found or the entry is stale.
    pub fn get(&mut self, key: &str) -> Option<&CacheEntry> {
        let pos = self.entries.iter().position(|e| e.key == key)?;
        let entry = self.entries.remove(pos)?;
        if entry.version < self.current_version {
            return None;
        }
        self.entries.push_front(entry);
        self.entries.front()
    }
    /// Increment the version, invalidating all current cache entries.
    pub fn invalidate_all(&mut self) {
        self.current_version += 1;
    }
    /// Remove a specific key from the cache.
    pub fn remove(&mut self, key: &str) {
        self.entries.retain(|e| e.key != key);
    }
    /// Number of entries currently in the cache.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Evict all stale entries (version < current_version).
    pub fn evict_stale(&mut self) {
        let ver = self.current_version;
        self.entries.retain(|e| e.version >= ver);
    }
}
/// In-memory CRUD store for `SimulationRecord` entries.
///
/// Supports insert, lookup by name, update, delete, and query filtering.
#[derive(Debug, Default)]
pub struct SimulationRecordDatabase {
    pub(super) records: HashMap<String, SimulationRecord>,
}
impl SimulationRecordDatabase {
    /// Create an empty database.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert or replace a record.
    pub fn insert(&mut self, rec: SimulationRecord) {
        self.records.insert(rec.name.clone(), rec);
    }
    /// Look up a record by name.
    pub fn get(&self, name: &str) -> Option<&SimulationRecord> {
        self.records.get(name)
    }
    /// Look up a record mutably by name.
    pub fn get_mut(&mut self, name: &str) -> Option<&mut SimulationRecord> {
        self.records.get_mut(name)
    }
    /// Delete a record by name.  Returns `true` if the record existed.
    pub fn delete(&mut self, name: &str) -> bool {
        self.records.remove(name).is_some()
    }
    /// Number of records.
    pub fn count(&self) -> usize {
        self.records.len()
    }
    /// Returns `true` if there are no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    /// Query: return all records that match the given `DatabaseQuery`.
    pub fn query(&self, q: &DatabaseQuery) -> Vec<&SimulationRecord> {
        self.records.values().filter(|r| q.matches(r)).collect()
    }
    /// Return all record names.
    pub fn names(&self) -> Vec<&str> {
        self.records.keys().map(|s| s.as_str()).collect()
    }
    /// Clear all records.
    pub fn clear(&mut self) {
        self.records.clear();
    }
}
/// Metadata for a single simulation snapshot.
#[derive(Debug, Clone)]
pub struct SnapshotEntry {
    /// Snapshot identifier (e.g. `"frame_0042"`).
    pub id: String,
    /// Simulation time at which this snapshot was taken (s).
    pub sim_time: f64,
    /// File path or URI for lazy loading (empty = not yet saved).
    pub path: String,
    /// File size in bytes (0 if unknown).
    pub file_size: usize,
    /// Whether this snapshot has been loaded into memory.
    pub loaded: bool,
    /// Optional per-snapshot tags.
    pub tags: Vec<String>,
}
impl SnapshotEntry {
    /// Create a new, unloaded snapshot entry.
    pub fn new(id: impl Into<String>, sim_time: f64, path: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            sim_time,
            path: path.into(),
            file_size: 0,
            loaded: false,
            tags: Vec::new(),
        }
    }
    /// Mark this snapshot as loaded.
    pub fn mark_loaded(&mut self) {
        self.loaded = true;
    }
    /// Add a tag.
    pub fn add_tag(&mut self, tag: impl Into<String>) {
        self.tags.push(tag.into());
    }
}
/// A simple in-memory relational table for simulation data.
///
/// Supports insert, column-based equality filtering, and projection.
#[derive(Debug, Default)]
pub struct SimulationDatabase {
    /// All rows stored in the table.
    pub rows: Vec<DbRow>,
    /// Table name for metadata/export purposes.
    pub name: String,
}
impl SimulationDatabase {
    /// Create a new, named, empty table.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            rows: Vec::new(),
            name: name.into(),
        }
    }
    /// Insert a row into the table.
    pub fn insert(&mut self, row: DbRow) {
        self.rows.push(row);
    }
    /// Count the number of rows.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }
    /// Query rows where column `col` equals `val`.
    pub fn query_eq(&self, col: &str, val: &DbValue) -> Vec<&DbRow> {
        self.rows
            .iter()
            .filter(|r| r.get(col).map(|v| v == val).unwrap_or(false))
            .collect()
    }
    /// Query rows where column `col` (numeric) is in `[lo, hi]`.
    pub fn query_range(&self, col: &str, lo: f64, hi: f64) -> Vec<&DbRow> {
        self.rows
            .iter()
            .filter(|r| {
                r.get(col)
                    .and_then(|v| v.as_f64())
                    .map(|f| f >= lo && f <= hi)
                    .unwrap_or(false)
            })
            .collect()
    }
    /// Return a projected view: for each row, extract the listed columns.
    pub fn project(&self, cols: &[&str]) -> Vec<HashMap<String, DbValue>> {
        self.rows
            .iter()
            .map(|r| {
                cols.iter()
                    .filter_map(|&c| r.get(c).map(|v| (c.to_string(), v.clone())))
                    .collect()
            })
            .collect()
    }
    /// Delete all rows where column `col` equals `val`.
    ///
    /// Returns the number of rows deleted.
    pub fn delete_eq(&mut self, col: &str, val: &DbValue) -> usize {
        let before = self.rows.len();
        self.rows
            .retain(|r| r.get(col).map(|v| v != val).unwrap_or(true));
        before - self.rows.len()
    }
    /// Clear all rows.
    pub fn clear(&mut self) {
        self.rows.clear();
    }
    /// Compute the mean of a numeric column across all rows.
    ///
    /// Returns `None` if no numeric values are found.
    pub fn column_mean(&self, col: &str) -> Option<f64> {
        let vals: Vec<f64> = self
            .rows
            .iter()
            .filter_map(|r| r.get(col).and_then(|v| v.as_f64()))
            .collect();
        if vals.is_empty() {
            None
        } else {
            Some(vals.iter().sum::<f64>() / vals.len() as f64)
        }
    }
}
/// Query predicate for filtering `SimulationRecord` entries.
#[derive(Debug, Clone, Default)]
pub struct DatabaseQuery {
    /// Keep only records whose timestamp >= this value.
    pub time_start: Option<f64>,
    /// Keep only records whose timestamp <= this value.
    pub time_end: Option<f64>,
    /// Keep only records whose name starts with this prefix.
    pub name_prefix: Option<String>,
    /// Keep only records that have param `key` equal to `value`.
    pub param_filter: Option<(String, String)>,
}
impl DatabaseQuery {
    /// Create an empty (pass-all) query.
    pub fn new() -> Self {
        Self::default()
    }
    /// Filter by time range `[t_start, t_end]`.
    pub fn with_time_range(mut self, t_start: f64, t_end: f64) -> Self {
        self.time_start = Some(t_start);
        self.time_end = Some(t_end);
        self
    }
    /// Filter by name prefix.
    pub fn with_name_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.name_prefix = Some(prefix.into());
        self
    }
    /// Filter by a parameter key=value match.
    pub fn with_param(mut self, key: impl Into<String>, val: impl Into<String>) -> Self {
        self.param_filter = Some((key.into(), val.into()));
        self
    }
    /// Check whether a `SimulationRecord` passes this query.
    pub fn matches(&self, rec: &SimulationRecord) -> bool {
        if let Some(t0) = self.time_start
            && rec.timestamp < t0
        {
            return false;
        }
        if let Some(t1) = self.time_end
            && rec.timestamp > t1
        {
            return false;
        }
        if let Some(ref prefix) = self.name_prefix
            && !rec.name.starts_with(prefix.as_str())
        {
            return false;
        }
        if let Some((ref k, ref v)) = self.param_filter {
            match rec.params.get(k.as_str()) {
                Some(pv) if pv == v => {}
                _ => return false,
            }
        }
        true
    }
}
/// A column value in the simulation database.
#[derive(Debug, Clone, PartialEq)]
pub enum DbValue {
    /// Integer value.
    Int(i64),
    /// Floating-point value.
    Float(f64),
    /// Text value.
    Text(String),
    /// Boolean value.
    Bool(bool),
    /// Null / missing value.
    Null,
}
impl DbValue {
    /// Return `Some(f64)` if this value is `Float` or `Int`, else `None`.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            DbValue::Float(v) => Some(*v),
            DbValue::Int(v) => Some(*v as f64),
            _ => None,
        }
    }
    /// Return `Some(&str)` if this value is `Text`, else `None`.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            DbValue::Text(s) => Some(s.as_str()),
            _ => None,
        }
    }
}
