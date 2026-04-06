// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Performance regression testing infrastructure.
//!
//! Provides a lightweight microbenchmark runner, baseline serialization
//! (JSON, no external dependency), and regression detection with
//! configurable thresholds.

use std::fmt;
use std::fs;
use std::io;
use std::path::Path;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

// ─────────────────────────────────────────────────────────────────────────────
// Error type
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during performance-regression workflows.
#[derive(Debug)]
pub enum PerfError {
    /// An I/O error (file read/write).
    Io(io::Error),
    /// JSON parsing / formatting error.
    Parse(String),
}

impl fmt::Display for PerfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PerfError::Io(e) => write!(f, "I/O error: {e}"),
            PerfError::Parse(msg) => write!(f, "parse error: {msg}"),
        }
    }
}

impl std::error::Error for PerfError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            PerfError::Io(e) => Some(e),
            PerfError::Parse(_) => None,
        }
    }
}

impl From<io::Error> for PerfError {
    fn from(e: io::Error) -> Self {
        PerfError::Io(e)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

/// A single benchmark result with timing data.
#[derive(Debug, Clone)]
pub struct BenchmarkResult {
    /// Human-readable benchmark name.
    pub name: String,
    /// Mean elapsed time in nanoseconds.
    pub mean_ns: f64,
    /// Standard deviation of elapsed time in nanoseconds.
    pub std_dev_ns: f64,
    /// Number of iterations executed.
    pub iterations: usize,
    /// Unix timestamp (seconds) when the benchmark was recorded.
    pub timestamp: u64,
}

/// Collection of benchmark results forming a baseline snapshot.
#[derive(Debug, Clone)]
pub struct BenchmarkBaseline {
    /// Individual results in this baseline.
    pub results: Vec<BenchmarkResult>,
    /// Git commit hash (or any identifier) for reproducibility.
    pub commit_hash: String,
    /// Human-readable date string (e.g. "2026-03-30").
    pub date: String,
}

/// Configuration knobs for regression detection.
#[derive(Debug, Clone)]
pub struct RegressionConfig {
    /// A percentage threshold: if a benchmark slows down by more than this
    /// percentage it is flagged as a regression.  For example `10.0` means
    /// a 10 % slowdown triggers an alert.
    pub threshold_pct: f64,
    /// Minimum number of iterations required for a result to be considered
    /// statistically meaningful.
    pub min_iterations: usize,
}

/// The full regression report produced by [`compare_baselines`].
#[derive(Debug, Clone)]
pub struct RegressionReport {
    /// Per-benchmark comparison entries.
    pub comparisons: Vec<BenchmarkComparison>,
    /// `true` if **any** comparison exceeds the configured threshold.
    pub has_regression: bool,
}

/// Side-by-side comparison of one benchmark between baseline and current.
#[derive(Debug, Clone)]
pub struct BenchmarkComparison {
    /// Benchmark name.
    pub name: String,
    /// Baseline mean in nanoseconds.
    pub baseline_ns: f64,
    /// Current mean in nanoseconds.
    pub current_ns: f64,
    /// Percentage change: positive means slower (regression).
    pub change_pct: f64,
    /// Whether this particular comparison constitutes a regression.
    pub is_regression: bool,
}

// ─────────────────────────────────────────────────────────────────────────────
// BenchmarkResult
// ─────────────────────────────────────────────────────────────────────────────

impl BenchmarkResult {
    /// Create a new benchmark result.  The timestamp is set to the current
    /// Unix time (seconds since epoch).
    pub fn new(name: &str, mean_ns: f64, std_dev_ns: f64, iterations: usize) -> Self {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();
        Self {
            name: name.to_string(),
            mean_ns,
            std_dev_ns,
            iterations,
            timestamp,
        }
    }

    /// Create a result with an explicit timestamp (useful for tests).
    pub fn with_timestamp(
        name: &str,
        mean_ns: f64,
        std_dev_ns: f64,
        iterations: usize,
        timestamp: u64,
    ) -> Self {
        Self {
            name: name.to_string(),
            mean_ns,
            std_dev_ns,
            iterations,
            timestamp,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BenchmarkBaseline – construction
// ─────────────────────────────────────────────────────────────────────────────

impl BenchmarkBaseline {
    /// Create a new empty baseline.
    pub fn new(commit_hash: &str, date: &str) -> Self {
        Self {
            results: Vec::new(),
            commit_hash: commit_hash.to_string(),
            date: date.to_string(),
        }
    }

    /// Append a single result.
    pub fn add_result(&mut self, result: BenchmarkResult) {
        self.results.push(result);
    }

    // ── JSON serialization (manual, no serde) ────────────────────────────

    /// Serialize this baseline to a JSON string.
    pub fn to_json(&self) -> String {
        let mut buf = String::with_capacity(512);
        buf.push_str("{\n");
        buf.push_str(&format!(
            "  \"commit_hash\": \"{}\",\n",
            escape_json(&self.commit_hash)
        ));
        buf.push_str(&format!("  \"date\": \"{}\",\n", escape_json(&self.date)));
        buf.push_str("  \"results\": [\n");
        for (i, r) in self.results.iter().enumerate() {
            buf.push_str("    {\n");
            buf.push_str(&format!("      \"name\": \"{}\",\n", escape_json(&r.name)));
            buf.push_str(&format!("      \"mean_ns\": {},\n", format_f64(r.mean_ns)));
            buf.push_str(&format!(
                "      \"std_dev_ns\": {},\n",
                format_f64(r.std_dev_ns)
            ));
            buf.push_str(&format!("      \"iterations\": {},\n", r.iterations));
            buf.push_str(&format!("      \"timestamp\": {}\n", r.timestamp));
            buf.push_str("    }");
            if i + 1 < self.results.len() {
                buf.push(',');
            }
            buf.push('\n');
        }
        buf.push_str("  ]\n");
        buf.push('}');
        buf
    }

    /// Deserialize a baseline from a JSON string previously produced by
    /// [`to_json`](Self::to_json).
    pub fn from_json(s: &str) -> Result<Self, PerfError> {
        let trimmed = s.trim();
        let commit_hash = extract_string_field(trimmed, "commit_hash")?;
        let date = extract_string_field(trimmed, "date")?;

        let results = parse_results_array(trimmed)?;

        Ok(Self {
            results,
            commit_hash,
            date,
        })
    }

    // ── File I/O ─────────────────────────────────────────────────────────

    /// Write the JSON representation to a file, creating parent
    /// directories if necessary.
    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), PerfError> {
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, self.to_json())?;
        Ok(())
    }

    /// Read and parse a baseline from a JSON file.
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> Result<Self, PerfError> {
        let contents = fs::read_to_string(path)?;
        Self::from_json(&contents)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RegressionConfig
// ─────────────────────────────────────────────────────────────────────────────

impl Default for RegressionConfig {
    fn default() -> Self {
        Self {
            threshold_pct: 10.0,
            min_iterations: 100,
        }
    }
}

impl RegressionConfig {
    /// Create a config with custom values.
    pub fn new(threshold_pct: f64, min_iterations: usize) -> Self {
        Self {
            threshold_pct,
            min_iterations,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Comparison logic
// ─────────────────────────────────────────────────────────────────────────────

/// Compare two baselines and produce a regression report.
///
/// For every benchmark that appears in **both** `current` and `baseline`
/// (matched by `name`), the percentage change is computed as:
///
/// ```text
/// change_pct = (current.mean_ns - baseline.mean_ns) / baseline.mean_ns * 100
/// ```
///
/// A positive `change_pct` means the current run is **slower** (potential
/// regression).  If the absolute change exceeds `config.threshold_pct` in
/// the positive direction **and** both results have at least
/// `config.min_iterations`, the comparison is marked as a regression.
pub fn compare_baselines(
    current: &BenchmarkBaseline,
    baseline: &BenchmarkBaseline,
    config: &RegressionConfig,
) -> RegressionReport {
    let mut comparisons = Vec::new();
    let mut has_regression = false;

    for cur in &current.results {
        if let Some(base) = baseline.results.iter().find(|b| b.name == cur.name) {
            let change_pct = if base.mean_ns.abs() < f64::EPSILON {
                0.0
            } else {
                (cur.mean_ns - base.mean_ns) / base.mean_ns * 100.0
            };

            let enough_iterations =
                cur.iterations >= config.min_iterations && base.iterations >= config.min_iterations;

            let is_regression = enough_iterations && change_pct > config.threshold_pct;

            if is_regression {
                has_regression = true;
            }

            comparisons.push(BenchmarkComparison {
                name: cur.name.clone(),
                baseline_ns: base.mean_ns,
                current_ns: cur.mean_ns,
                change_pct,
                is_regression,
            });
        }
    }

    RegressionReport {
        comparisons,
        has_regression,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RegressionReport – human-readable summary
// ─────────────────────────────────────────────────────────────────────────────

impl RegressionReport {
    /// Produce a human-readable, multi-line summary of the report.
    pub fn summary(&self) -> String {
        let mut buf = String::with_capacity(256);
        buf.push_str("=== Performance Regression Report ===\n\n");

        if self.comparisons.is_empty() {
            buf.push_str("No benchmarks compared.\n");
            return buf;
        }

        for cmp in &self.comparisons {
            let direction = if cmp.change_pct > 0.0 {
                "slower"
            } else if cmp.change_pct < 0.0 {
                "faster"
            } else {
                "same"
            };

            let marker = if cmp.is_regression {
                " [REGRESSION]"
            } else {
                ""
            };

            buf.push_str(&format!(
                "  {}: {:.2} ns -> {:.2} ns ({:+.2}% {}){}",
                cmp.name, cmp.baseline_ns, cmp.current_ns, cmp.change_pct, direction, marker,
            ));
            buf.push('\n');
        }

        buf.push('\n');
        if self.has_regression {
            buf.push_str("RESULT: REGRESSION DETECTED\n");
        } else {
            buf.push_str("RESULT: No regression detected.\n");
        }

        buf
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Microbenchmark runner
// ─────────────────────────────────────────────────────────────────────────────

/// Run a closure `iterations` times and return a [`BenchmarkResult`] with
/// mean and standard-deviation timing data.
///
/// The runner performs a short warm-up (10 % of `iterations`, clamped to
/// 1..=1000) before timing begins.
pub fn microbench<F: FnMut()>(name: &str, iterations: usize, mut f: F) -> BenchmarkResult {
    // Warm-up
    let warmup = (iterations / 10).clamp(1, 1000);
    for _ in 0..warmup {
        f();
    }

    // Timed runs – collect individual durations
    let mut durations_ns: Vec<f64> = Vec::with_capacity(iterations);
    for _ in 0..iterations {
        let start = Instant::now();
        f();
        let elapsed = start.elapsed();
        durations_ns.push(elapsed.as_nanos() as f64);
    }

    let n = durations_ns.len() as f64;
    let mean = durations_ns.iter().sum::<f64>() / n;

    let variance = if durations_ns.len() > 1 {
        durations_ns.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / (n - 1.0)
    } else {
        0.0
    };
    let std_dev = variance.sqrt();

    BenchmarkResult::new(name, mean, std_dev, iterations)
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers – JSON
// ─────────────────────────────────────────────────────────────────────────────

fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out
}

fn unescape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(ch) = chars.next() {
        if ch == '\\' {
            match chars.next() {
                Some('"') => out.push('"'),
                Some('\\') => out.push('\\'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('t') => out.push('\t'),
                Some(other) => {
                    out.push('\\');
                    out.push(other);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(ch);
        }
    }
    out
}

/// Format an f64 for JSON output, ensuring it always has a decimal point.
fn format_f64(v: f64) -> String {
    let s = format!("{v}");
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{v}.0")
    }
}

/// Extract a JSON string field value by key (simple, non-nested).
fn extract_string_field(json: &str, key: &str) -> Result<String, PerfError> {
    let pattern = format!("\"{}\"", key);
    let key_pos = json
        .find(&pattern)
        .ok_or_else(|| PerfError::Parse(format!("missing field \"{key}\"")))?;
    let after_key = &json[key_pos + pattern.len()..];
    // skip optional whitespace and the colon
    let after_colon = after_key
        .find(':')
        .map(|i| &after_key[i + 1..])
        .ok_or_else(|| PerfError::Parse(format!("expected ':' after \"{key}\"")))?;
    let after_colon = after_colon.trim_start();
    if !after_colon.starts_with('"') {
        return Err(PerfError::Parse(format!(
            "expected string value for \"{key}\""
        )));
    }
    let value_start = 1; // skip opening quote
    let rest = &after_colon[value_start..];
    let mut end = 0;
    let mut escaped = false;
    for ch in rest.chars() {
        if escaped {
            escaped = false;
            end += ch.len_utf8();
            continue;
        }
        if ch == '\\' {
            escaped = true;
            end += 1;
            continue;
        }
        if ch == '"' {
            break;
        }
        end += ch.len_utf8();
    }
    Ok(unescape_json(&rest[..end]))
}

/// Extract a JSON number (f64) field.
fn extract_f64_field(json: &str, key: &str) -> Result<f64, PerfError> {
    let pattern = format!("\"{}\"", key);
    let key_pos = json
        .find(&pattern)
        .ok_or_else(|| PerfError::Parse(format!("missing field \"{key}\"")))?;
    let after_key = &json[key_pos + pattern.len()..];
    let after_colon = after_key
        .find(':')
        .map(|i| &after_key[i + 1..])
        .ok_or_else(|| PerfError::Parse(format!("expected ':' after \"{key}\"")))?;
    let after_colon = after_colon.trim_start();
    // read until comma, '}', or whitespace
    let end = after_colon
        .find([',', '}', ']', '\n'])
        .unwrap_or(after_colon.len());
    let num_str = after_colon[..end].trim();
    num_str
        .parse::<f64>()
        .map_err(|e| PerfError::Parse(format!("invalid f64 for \"{key}\": {e}")))
}

/// Extract a JSON integer (usize) field.
fn extract_usize_field(json: &str, key: &str) -> Result<usize, PerfError> {
    let pattern = format!("\"{}\"", key);
    let key_pos = json
        .find(&pattern)
        .ok_or_else(|| PerfError::Parse(format!("missing field \"{key}\"")))?;
    let after_key = &json[key_pos + pattern.len()..];
    let after_colon = after_key
        .find(':')
        .map(|i| &after_key[i + 1..])
        .ok_or_else(|| PerfError::Parse(format!("expected ':' after \"{key}\"")))?;
    let after_colon = after_colon.trim_start();
    let end = after_colon
        .find([',', '}', ']', '\n'])
        .unwrap_or(after_colon.len());
    let num_str = after_colon[..end].trim();
    num_str
        .parse::<usize>()
        .map_err(|e| PerfError::Parse(format!("invalid usize for \"{key}\": {e}")))
}

/// Extract a JSON u64 field.
fn extract_u64_field(json: &str, key: &str) -> Result<u64, PerfError> {
    let pattern = format!("\"{}\"", key);
    let key_pos = json
        .find(&pattern)
        .ok_or_else(|| PerfError::Parse(format!("missing field \"{key}\"")))?;
    let after_key = &json[key_pos + pattern.len()..];
    let after_colon = after_key
        .find(':')
        .map(|i| &after_key[i + 1..])
        .ok_or_else(|| PerfError::Parse(format!("expected ':' after \"{key}\"")))?;
    let after_colon = after_colon.trim_start();
    let end = after_colon
        .find([',', '}', ']', '\n'])
        .unwrap_or(after_colon.len());
    let num_str = after_colon[..end].trim();
    num_str
        .parse::<u64>()
        .map_err(|e| PerfError::Parse(format!("invalid u64 for \"{key}\": {e}")))
}

/// Parse the `"results"` array from the top-level JSON object.
fn parse_results_array(json: &str) -> Result<Vec<BenchmarkResult>, PerfError> {
    let results_key = json
        .find("\"results\"")
        .ok_or_else(|| PerfError::Parse("missing \"results\" array".to_string()))?;
    let after_key = &json[results_key + "\"results\"".len()..];
    let bracket_start = after_key
        .find('[')
        .ok_or_else(|| PerfError::Parse("expected '[' after \"results\"".to_string()))?;
    let array_content = &after_key[bracket_start + 1..];

    // Find matching ']'
    let bracket_end = find_matching_bracket(array_content)?;
    let array_str = &array_content[..bracket_end];

    // Split by top-level objects (simple brace matching)
    let objects = split_json_objects(array_str);

    let mut results = Vec::new();
    for obj in &objects {
        let name = extract_string_field(obj, "name")?;
        let mean_ns = extract_f64_field(obj, "mean_ns")?;
        let std_dev_ns = extract_f64_field(obj, "std_dev_ns")?;
        let iterations = extract_usize_field(obj, "iterations")?;
        let timestamp = extract_u64_field(obj, "timestamp")?;
        results.push(BenchmarkResult {
            name,
            mean_ns,
            std_dev_ns,
            iterations,
            timestamp,
        });
    }

    Ok(results)
}

/// Find the position of the matching ']' for content that starts right
/// after an opening '\['.
fn find_matching_bracket(s: &str) -> Result<usize, PerfError> {
    let mut depth: i32 = 0;
    let mut in_string = false;
    let mut escaped = false;
    for (i, ch) in s.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if ch == '\\' && in_string {
            escaped = true;
            continue;
        }
        if ch == '"' {
            in_string = !in_string;
            continue;
        }
        if in_string {
            continue;
        }
        if ch == '[' {
            depth += 1;
        } else if ch == ']' {
            if depth == 0 {
                return Ok(i);
            }
            depth -= 1;
        }
    }
    Err(PerfError::Parse("unmatched '['".to_string()))
}

/// Split a JSON array body into individual object strings.
fn split_json_objects(s: &str) -> Vec<String> {
    let mut objects = Vec::new();
    let mut depth: i32 = 0;
    let mut in_string = false;
    let mut escaped = false;
    let mut start: Option<usize> = None;

    for (i, ch) in s.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if ch == '\\' && in_string {
            escaped = true;
            continue;
        }
        if ch == '"' {
            in_string = !in_string;
            continue;
        }
        if in_string {
            continue;
        }
        if ch == '{' {
            if depth == 0 {
                start = Some(i);
            }
            depth += 1;
        } else if ch == '}' {
            depth -= 1;
            if depth == 0
                && let Some(s_idx) = start
            {
                objects.push(s[s_idx..=i].to_string());
                start = None;
            }
        }
    }
    objects
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialize_deserialize_roundtrip() {
        let mut baseline = BenchmarkBaseline::new("abc123", "2026-03-30");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "broadphase_100",
            1500.5,
            120.3,
            200,
            1711800000,
        ));
        baseline.add_result(BenchmarkResult::with_timestamp(
            "gjk_intersect",
            850.0,
            45.2,
            500,
            1711800001,
        ));

        let json = baseline.to_json();
        let restored =
            BenchmarkBaseline::from_json(&json).expect("round-trip deserialization failed");

        assert_eq!(restored.commit_hash, "abc123");
        assert_eq!(restored.date, "2026-03-30");
        assert_eq!(restored.results.len(), 2);
        assert_eq!(restored.results[0].name, "broadphase_100");
        assert!((restored.results[0].mean_ns - 1500.5).abs() < 1e-9);
        assert!((restored.results[0].std_dev_ns - 120.3).abs() < 1e-9);
        assert_eq!(restored.results[0].iterations, 200);
        assert_eq!(restored.results[0].timestamp, 1711800000);
        assert_eq!(restored.results[1].name, "gjk_intersect");
    }

    #[test]
    fn test_regression_detected() {
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "bench_a", 1000.0, 10.0, 200, 0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        // 15 % slower -> should regress at 10 % threshold
        current.add_result(BenchmarkResult::with_timestamp(
            "bench_a", 1150.0, 12.0, 200, 1,
        ));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        assert!(report.has_regression);
        assert_eq!(report.comparisons.len(), 1);
        assert!(report.comparisons[0].is_regression);
        assert!((report.comparisons[0].change_pct - 15.0).abs() < 1e-9);
    }

    #[test]
    fn test_no_regression_identical_results() {
        let mut baseline = BenchmarkBaseline::new("v1", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "bench_x", 500.0, 5.0, 200, 0,
        ));

        let mut current = BenchmarkBaseline::new("v2", "2026-03-30");
        current.add_result(BenchmarkResult::with_timestamp(
            "bench_x", 500.0, 5.0, 200, 1,
        ));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        assert!(!report.has_regression);
        assert_eq!(report.comparisons.len(), 1);
        assert!(!report.comparisons[0].is_regression);
        assert!(report.comparisons[0].change_pct.abs() < 1e-9);
    }

    #[test]
    fn test_improvement_detected() {
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "bench_fast",
            2000.0,
            20.0,
            300,
            0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        // 25 % faster
        current.add_result(BenchmarkResult::with_timestamp(
            "bench_fast",
            1500.0,
            15.0,
            300,
            1,
        ));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        assert!(!report.has_regression);
        assert_eq!(report.comparisons.len(), 1);
        assert!(!report.comparisons[0].is_regression);
        assert!(report.comparisons[0].change_pct < 0.0); // negative = improvement
    }

    #[test]
    fn test_microbench_returns_valid_data() {
        let mut counter = 0u64;
        let result = microbench("counter_test", 100, || {
            counter += 1;
        });

        assert_eq!(result.name, "counter_test");
        assert_eq!(result.iterations, 100);
        assert!(result.mean_ns >= 0.0);
        assert!(result.std_dev_ns >= 0.0);
        assert!(result.timestamp > 0);
        // The closure should have run iterations + warmup times
        assert!(counter >= 100);
    }

    #[test]
    fn test_report_summary_formatting() {
        let report = RegressionReport {
            has_regression: true,
            comparisons: vec![
                BenchmarkComparison {
                    name: "bench_a".to_string(),
                    baseline_ns: 1000.0,
                    current_ns: 1150.0,
                    change_pct: 15.0,
                    is_regression: true,
                },
                BenchmarkComparison {
                    name: "bench_b".to_string(),
                    baseline_ns: 800.0,
                    current_ns: 750.0,
                    change_pct: -6.25,
                    is_regression: false,
                },
            ],
        };

        let summary = report.summary();
        assert!(summary.contains("Performance Regression Report"));
        assert!(summary.contains("bench_a"));
        assert!(summary.contains("[REGRESSION]"));
        assert!(summary.contains("bench_b"));
        assert!(summary.contains("REGRESSION DETECTED"));
    }

    #[test]
    fn test_report_summary_no_regression() {
        let report = RegressionReport {
            has_regression: false,
            comparisons: vec![BenchmarkComparison {
                name: "ok_bench".to_string(),
                baseline_ns: 100.0,
                current_ns: 105.0,
                change_pct: 5.0,
                is_regression: false,
            }],
        };

        let summary = report.summary();
        assert!(summary.contains("No regression detected"));
        assert!(!summary.contains("[REGRESSION]"));
    }

    #[test]
    fn test_file_roundtrip() {
        let dir = std::env::temp_dir().join("oxiphysics_perf_test");
        let path = dir.join("baseline.json");

        let mut baseline = BenchmarkBaseline::new("deadbeef", "2026-03-30");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "file_bench",
            999.9,
            11.1,
            150,
            1711800000,
        ));

        baseline
            .save_to_file(&path)
            .expect("failed to save baseline");
        let loaded = BenchmarkBaseline::load_from_file(&path).expect("failed to load baseline");

        assert_eq!(loaded.commit_hash, "deadbeef");
        assert_eq!(loaded.results.len(), 1);
        assert_eq!(loaded.results[0].name, "file_bench");
        assert!((loaded.results[0].mean_ns - 999.9).abs() < 1e-9);

        // Cleanup
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_regression_below_min_iterations() {
        // Even though the slowdown is large, insufficient iterations should
        // not flag a regression.
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "small_bench",
            100.0,
            5.0,
            10,
            0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        current.add_result(BenchmarkResult::with_timestamp(
            "small_bench",
            200.0,
            10.0,
            10,
            1,
        ));

        let config = RegressionConfig::default(); // min_iterations = 100
        let report = compare_baselines(&current, &baseline, &config);

        assert!(!report.has_regression);
        assert!(!report.comparisons[0].is_regression);
    }

    #[test]
    fn test_json_escape_roundtrip() {
        let mut baseline = BenchmarkBaseline::new("hash\"with\\special", "2026-03-30");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "bench\tnewline\ntest",
            100.0,
            1.0,
            100,
            0,
        ));

        let json = baseline.to_json();
        let restored = BenchmarkBaseline::from_json(&json).expect("escaped round-trip failed");

        assert_eq!(restored.commit_hash, "hash\"with\\special");
        assert_eq!(restored.results[0].name, "bench\tnewline\ntest");
    }

    #[test]
    fn test_empty_baseline_serialization() {
        let baseline = BenchmarkBaseline::new("empty", "2026-03-30");
        let json = baseline.to_json();
        let restored =
            BenchmarkBaseline::from_json(&json).expect("empty baseline round-trip failed");

        assert_eq!(restored.commit_hash, "empty");
        assert!(restored.results.is_empty());
    }

    #[test]
    fn test_compare_mismatched_names() {
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp("alpha", 100.0, 1.0, 200, 0));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        current.add_result(BenchmarkResult::with_timestamp("beta", 100.0, 1.0, 200, 1));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        // No matching names -> no comparisons
        assert!(report.comparisons.is_empty());
        assert!(!report.has_regression);
    }

    #[test]
    fn test_compare_zero_baseline_mean() {
        // Guard against division by zero.
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "zero_mean",
            0.0,
            0.0,
            200,
            0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        current.add_result(BenchmarkResult::with_timestamp(
            "zero_mean",
            10.0,
            1.0,
            200,
            1,
        ));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        assert_eq!(report.comparisons.len(), 1);
        assert!((report.comparisons[0].change_pct).abs() < 1e-9);
    }

    #[test]
    fn test_multiple_benchmarks_mixed() {
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "fast_one", 100.0, 2.0, 200, 0,
        ));
        baseline.add_result(BenchmarkResult::with_timestamp(
            "slow_one", 500.0, 10.0, 200, 0,
        ));
        baseline.add_result(BenchmarkResult::with_timestamp(
            "stable", 300.0, 5.0, 200, 0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        // Improved
        current.add_result(BenchmarkResult::with_timestamp(
            "fast_one", 80.0, 1.5, 200, 1,
        ));
        // Regressed (20 %)
        current.add_result(BenchmarkResult::with_timestamp(
            "slow_one", 600.0, 12.0, 200, 1,
        ));
        // Within threshold (3.3 %)
        current.add_result(BenchmarkResult::with_timestamp(
            "stable", 310.0, 5.0, 200, 1,
        ));

        let config = RegressionConfig::default();
        let report = compare_baselines(&current, &baseline, &config);

        assert!(report.has_regression);
        assert_eq!(report.comparisons.len(), 3);

        let slow = report
            .comparisons
            .iter()
            .find(|c| c.name == "slow_one")
            .expect("slow_one missing");
        assert!(slow.is_regression);

        let fast = report
            .comparisons
            .iter()
            .find(|c| c.name == "fast_one")
            .expect("fast_one missing");
        assert!(!fast.is_regression);

        let stable = report
            .comparisons
            .iter()
            .find(|c| c.name == "stable")
            .expect("stable missing");
        assert!(!stable.is_regression);
    }

    #[test]
    fn test_custom_threshold() {
        let mut baseline = BenchmarkBaseline::new("old", "2026-01-01");
        baseline.add_result(BenchmarkResult::with_timestamp(
            "tight", 1000.0, 10.0, 200, 0,
        ));

        let mut current = BenchmarkBaseline::new("new", "2026-03-30");
        // 3 % slowdown
        current.add_result(BenchmarkResult::with_timestamp(
            "tight", 1030.0, 11.0, 200, 1,
        ));

        // Default 10 % -> no regression
        let default_report = compare_baselines(&current, &baseline, &RegressionConfig::default());
        assert!(!default_report.has_regression);

        // Strict 2 % -> regression
        let strict = RegressionConfig::new(2.0, 100);
        let strict_report = compare_baselines(&current, &baseline, &strict);
        assert!(strict_report.has_regression);
    }
}
