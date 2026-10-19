// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Material database I/O: CES-style material properties, JSON material libraries.
//!
//! Provides typed property storage, a searchable material database, JSON
//! serialization/deserialization, temperature-dependent property interpolation,
//! and CSV import (Matweb-style).

use std::collections::HashMap;

// ── MaterialProperty ──────────────────────────────────────────────────────────

/// A typed material property value.
///
/// Each variant carries the numeric content; metadata (units, source) lives in
/// the parent `MaterialRecord`.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialProperty {
    /// A single scalar value.
    Scalar(f64),
    /// A min/max range (e.g. yield strength from 200 to 350 MPa).
    Range(f64, f64),
    /// A flat row-major matrix or tensor (e.g. stiffness tensor).
    Matrix(Vec<f64>),
}

impl MaterialProperty {
    /// Return the representative (mid-point or scalar) value for comparisons.
    pub fn representative(&self) -> f64 {
        match self {
            MaterialProperty::Scalar(v) => *v,
            MaterialProperty::Range(lo, hi) => 0.5 * (lo + hi),
            MaterialProperty::Matrix(v) => {
                if v.is_empty() {
                    0.0
                } else {
                    v[0]
                }
            }
        }
    }

    /// Check whether a value falls within this property (scalar exact, range inclusive).
    pub fn contains(&self, value: f64) -> bool {
        match self {
            MaterialProperty::Scalar(v) => (v - value).abs() < 1e-12,
            MaterialProperty::Range(lo, hi) => value >= *lo && value <= *hi,
            MaterialProperty::Matrix(_) => false,
        }
    }
}

// ── MaterialRecord ────────────────────────────────────────────────────────────

/// A single material entry: name, category, and a map of named properties.
///
/// Typical property keys: `"density"`, `"elastic_modulus"`, `"yield_strength"`,
/// `"thermal_conductivity"`, `"poisson_ratio"`.
#[derive(Debug, Clone)]
pub struct MaterialRecord {
    /// Human-readable material name (e.g. `"AISI 304 Stainless Steel"`).
    pub name: String,
    /// Hierarchical category path separated by `>` (e.g. `"Metal>Ferrous>Steel"`).
    pub category: String,
    /// Named properties map.
    pub properties: HashMap<String, MaterialProperty>,
    /// Optional source/reference string.
    pub source: String,
}

impl MaterialRecord {
    /// Create a new empty `MaterialRecord`.
    pub fn new(name: &str, category: &str) -> Self {
        Self {
            name: name.to_string(),
            category: category.to_string(),
            properties: HashMap::new(),
            source: String::new(),
        }
    }

    /// Insert a scalar property.
    pub fn set_scalar(&mut self, key: &str, value: f64) {
        self.properties
            .insert(key.to_string(), MaterialProperty::Scalar(value));
    }

    /// Insert a range property.
    pub fn set_range(&mut self, key: &str, lo: f64, hi: f64) {
        self.properties
            .insert(key.to_string(), MaterialProperty::Range(lo, hi));
    }

    /// Retrieve the representative value of a named property, or `None`.
    pub fn get_value(&self, key: &str) -> Option<f64> {
        self.properties.get(key).map(|p| p.representative())
    }
}

// ── MaterialDatabase ──────────────────────────────────────────────────────────

/// An in-memory collection of `MaterialRecord` entries with search capabilities.
#[derive(Debug, Clone, Default)]
pub struct MaterialDatabase {
    /// All stored records.
    pub records: Vec<MaterialRecord>,
}

impl MaterialDatabase {
    /// Create a new empty database.
    pub fn new() -> Self {
        Self {
            records: Vec::new(),
        }
    }

    /// Add a record to the database.
    pub fn add(&mut self, record: MaterialRecord) {
        self.records.push(record);
    }

    /// Search for records whose `category` starts with the given prefix.
    pub fn by_category(&self, prefix: &str) -> Vec<&MaterialRecord> {
        self.records
            .iter()
            .filter(|r| r.category.starts_with(prefix))
            .collect()
    }

    /// Filter records where a named property's representative value falls in `[lo, hi]`.
    pub fn filter_by_property(&self, key: &str, lo: f64, hi: f64) -> Vec<&MaterialRecord> {
        self.records
            .iter()
            .filter(|r| {
                if let Some(v) = r.get_value(key) {
                    v >= lo && v <= hi
                } else {
                    false
                }
            })
            .collect()
    }

    /// Generate Ashby chart data: pairs of (prop_x, prop_y) for each record that
    /// has both properties.
    pub fn ashby_chart(&self, prop_x: &str, prop_y: &str) -> Vec<(f64, f64, &str)> {
        self.records
            .iter()
            .filter_map(|r| {
                let x = r.get_value(prop_x)?;
                let y = r.get_value(prop_y)?;
                Some((x, y, r.name.as_str()))
            })
            .collect()
    }
}

// ── MaterialCategories ────────────────────────────────────────────────────────

/// Hierarchical material category tree node.
///
/// Each node has a name and optional children, representing the tree:
/// `Metal > Ferrous > Steel > AISI 304`.
#[derive(Debug, Clone)]
pub struct MaterialCategories {
    /// Category name (leaf or node).
    pub name: String,
    /// Child subcategories.
    pub children: Vec<MaterialCategories>,
}

impl MaterialCategories {
    /// Create a leaf category node.
    pub fn leaf(name: &str) -> Self {
        Self {
            name: name.to_string(),
            children: Vec::new(),
        }
    }

    /// Create a node with children.
    pub fn node(name: &str, children: Vec<MaterialCategories>) -> Self {
        Self {
            name: name.to_string(),
            children,
        }
    }

    /// Check whether this node or any descendant has the given name.
    pub fn contains_name(&self, target: &str) -> bool {
        if self.name == target {
            return true;
        }
        self.children.iter().any(|c| c.contains_name(target))
    }

    /// Collect all leaf names under this node.
    pub fn leaf_names(&self) -> Vec<String> {
        if self.children.is_empty() {
            vec![self.name.clone()]
        } else {
            self.children.iter().flat_map(|c| c.leaf_names()).collect()
        }
    }
}

// ── MaterialComparison ────────────────────────────────────────────────────────

/// Compare materials using performance indices and radar chart data.
#[derive(Debug, Clone)]
pub struct MaterialComparison {
    /// Property keys used for comparison.
    pub axes: Vec<String>,
}

impl MaterialComparison {
    /// Create a `MaterialComparison` for the given set of property axes.
    pub fn new(axes: Vec<String>) -> Self {
        Self { axes }
    }

    /// Compute a performance index for a record.
    ///
    /// The index is the ratio of two properties: `numerator / denominator`.
    /// Returns `None` if either property is missing or denominator is zero.
    pub fn performance_index(
        &self,
        record: &MaterialRecord,
        numerator: &str,
        denominator: &str,
    ) -> Option<f64> {
        let num = record.get_value(numerator)?;
        let den = record.get_value(denominator)?;
        if den.abs() < 1e-30 {
            return None;
        }
        Some(num / den)
    }

    /// Rank materials by performance index (descending, highest first).
    pub fn rank_by_index<'a>(
        &self,
        records: &[&'a MaterialRecord],
        numerator: &str,
        denominator: &str,
    ) -> Vec<(&'a MaterialRecord, f64)> {
        let mut scored: Vec<_> = records
            .iter()
            .filter_map(|&r| {
                let idx = self.performance_index(r, numerator, denominator)?;
                Some((r, idx))
            })
            .collect();
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        scored
    }

    /// Produce radar chart data for a single material: normalized values on each axis.
    ///
    /// Each value is divided by `reference_values[i]` for normalization.
    pub fn radar_data(&self, record: &MaterialRecord, reference_values: &[f64]) -> Vec<f64> {
        self.axes
            .iter()
            .enumerate()
            .map(|(i, key)| {
                let val = record.get_value(key).unwrap_or(0.0);
                let ref_val = if i < reference_values.len() {
                    reference_values[i]
                } else {
                    1.0
                };
                if ref_val.abs() < 1e-30 {
                    0.0
                } else {
                    val / ref_val
                }
            })
            .collect()
    }
}

// ── TemperatureDependence ─────────────────────────────────────────────────────

/// Temperature-dependent material property via piecewise linear interpolation.
///
/// Tabulated as `(temperature, value)` pairs sorted by temperature.
#[derive(Debug, Clone)]
pub struct TemperatureDependence {
    /// Sorted (temperature, value) table.
    pub table: Vec<(f64, f64)>,
    /// Property name this table describes.
    pub property_name: String,
}

impl TemperatureDependence {
    /// Create a new `TemperatureDependence` table.
    pub fn new(property_name: &str, mut table: Vec<(f64, f64)>) -> Self {
        table.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self {
            table,
            property_name: property_name.to_string(),
        }
    }

    /// Interpolate property value at the given temperature.
    ///
    /// Clamps to the first/last value outside the tabulated range.
    pub fn value_at(&self, temperature: f64) -> f64 {
        if self.table.is_empty() {
            return 0.0;
        }
        if self.table.len() == 1 {
            return self.table[0].1;
        }
        let first = &self.table[0];
        let last = &self.table[self.table.len() - 1];
        if temperature <= first.0 {
            return first.1;
        }
        if temperature >= last.0 {
            return last.1;
        }
        // Binary search for the bracket
        let mut lo = 0usize;
        let mut hi = self.table.len() - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.table[mid].0 <= temperature {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let (t0, v0) = self.table[lo];
        let (t1, v1) = self.table[hi];
        let u = (temperature - t0) / (t1 - t0);
        v0 + u * (v1 - v0)
    }
}

// ── MaterialFilter ────────────────────────────────────────────────────────────

/// Chainable filter for material databases.
///
/// Filters can be composed: call `filter_category`, `filter_min`, `filter_max`
/// in sequence and collect results with `apply`.
#[derive(Debug, Clone, Default)]
pub struct MaterialFilter {
    /// Optional category prefix filter.
    pub category_prefix: Option<String>,
    /// Per-property minimum value constraints.
    pub min_values: HashMap<String, f64>,
    /// Per-property maximum value constraints.
    pub max_values: HashMap<String, f64>,
}

impl MaterialFilter {
    /// Create an empty filter (passes everything).
    pub fn new() -> Self {
        Self::default()
    }

    /// Restrict to materials whose category starts with `prefix`.
    pub fn filter_category(mut self, prefix: &str) -> Self {
        self.category_prefix = Some(prefix.to_string());
        self
    }

    /// Require property `key` >= `min`.
    pub fn filter_min(mut self, key: &str, min: f64) -> Self {
        self.min_values.insert(key.to_string(), min);
        self
    }

    /// Require property `key` <= `max`.
    pub fn filter_max(mut self, key: &str, max: f64) -> Self {
        self.max_values.insert(key.to_string(), max);
        self
    }

    /// Apply the filter to a database and return matching records.
    pub fn apply<'a>(&self, db: &'a MaterialDatabase) -> Vec<&'a MaterialRecord> {
        db.records
            .iter()
            .filter(|r| {
                if let Some(prefix) = &self.category_prefix
                    && !r.category.starts_with(prefix.as_str())
                {
                    return false;
                }
                for (key, &min_v) in &self.min_values {
                    if let Some(v) = r.get_value(key) {
                        if v < min_v {
                            return false;
                        }
                    } else {
                        return false;
                    }
                }
                for (key, &max_v) in &self.max_values {
                    if let Some(v) = r.get_value(key) {
                        if v > max_v {
                            return false;
                        }
                    } else {
                        return false;
                    }
                }
                true
            })
            .collect()
    }
}

// ── JsonMaterialDb ─────────────────────────────────────────────────────────────

/// JSON serialization and deserialization of a [`MaterialRecord`].
///
/// [`JsonMaterialDb::write_json_material`] writes a small, human-readable JSON
/// object and [`JsonMaterialDb::parse_json_material`] reads it back with the
/// `serde_json` parser, so any valid JSON text is accepted, whatever its
/// whitespace or key order, and strings may hold any Unicode text:
///
/// ```text
/// {
///   "name": "Ti-6Al-4V (α+β)",
///   "category": "Metal>NonFerrous>Titanium",
///   "source": "Matweb",
///   "properties": {
///     "density": {"type": "scalar", "value": 4430.0},
///     "yield_strength": {"type": "range", "lo": 800000000.0, "hi": 1000000000.0},
///     "stiffness": {"type": "matrix", "values": [1.0, 2.0, 3.0]}
///   }
/// }
/// ```
///
/// Numbers are written in the shortest decimal form that parses back to the
/// identical `f64`.  JSON has no spelling for NaN and the infinities, so they
/// are written as the strings `"NaN"`, `"Infinity"` and `"-Infinity"`, which
/// the reader accepts wherever it expects a number.
pub struct JsonMaterialDb;

impl JsonMaterialDb {
    /// Serialize a `MaterialRecord` to a JSON object string.
    ///
    /// The properties are written in key order, so the text of a record is the
    /// same on every run.
    pub fn write_json_material(record: &MaterialRecord) -> String {
        let mut out = String::from("{\n");
        out.push_str(&format!("  \"name\": \"{}\",\n", escape_json(&record.name)));
        out.push_str(&format!(
            "  \"category\": \"{}\",\n",
            escape_json(&record.category)
        ));
        out.push_str(&format!(
            "  \"source\": \"{}\",\n",
            escape_json(&record.source)
        ));
        out.push_str("  \"properties\": {\n");
        let mut entries: Vec<(&String, &MaterialProperty)> = record.properties.iter().collect();
        entries.sort_by(|a, b| a.0.cmp(b.0));
        let last = entries.len().saturating_sub(1);
        for (index, (key, property)) in entries.into_iter().enumerate() {
            let comma = if index < last { "," } else { "" };
            let body = match property {
                MaterialProperty::Scalar(value) => {
                    format!("\"type\": \"scalar\", \"value\": {}", json_number(*value))
                }
                MaterialProperty::Range(lo, hi) => format!(
                    "\"type\": \"range\", \"lo\": {}, \"hi\": {}",
                    json_number(*lo),
                    json_number(*hi)
                ),
                MaterialProperty::Matrix(values) => {
                    let items: Vec<String> = values.iter().map(|x| json_number(*x)).collect();
                    format!("\"type\": \"matrix\", \"values\": [{}]", items.join(", "))
                }
            };
            out.push_str(&format!(
                "    \"{}\": {{{}}}{}\n",
                escape_json(key),
                body,
                comma
            ));
        }
        out.push_str("  }\n");
        out.push('}');
        out
    }

    /// Parse a single material record from a JSON string.
    ///
    /// Reads the format produced by [`JsonMaterialDb::write_json_material`]:
    /// a JSON object with a string `"name"`, optional string `"category"` and
    /// `"source"` (empty when absent) and a `"properties"` object whose members
    /// are `{"type": "scalar", "value": v}`, `{"type": "range", "lo": a, "hi": b}`
    /// or `{"type": "matrix", "values": [..]}`.  Other members are ignored.
    ///
    /// Returns `None` when the text is not valid JSON or does not have that
    /// shape (a missing name, an unknown property type, a missing or
    /// non-numeric value, ...); a record is never returned half-parsed.
    pub fn parse_json_material(json: &str) -> Option<MaterialRecord> {
        let root: serde_json::Value = serde_json::from_str(json).ok()?;
        let fields = root.as_object()?;
        let name = fields.get("name")?.as_str()?;
        let category = optional_text(fields, "category")?;
        let source = optional_text(fields, "source")?;
        let properties = fields.get("properties")?.as_object()?;
        let mut record = MaterialRecord::new(name, category);
        record.source = source.to_string();
        for (key, entry) in properties {
            record
                .properties
                .insert(key.clone(), parse_property(entry)?);
        }
        Some(record)
    }
}

/// A string member that may be absent (`""`); present but not a string is `None`.
fn optional_text<'a>(
    fields: &'a serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> Option<&'a str> {
    match fields.get(key) {
        None => Some(""),
        Some(value) => value.as_str(),
    }
}

/// Decode one member of the `"properties"` object, dispatching on its `"type"`.
fn parse_property(entry: &serde_json::Value) -> Option<MaterialProperty> {
    let fields = entry.as_object()?;
    match fields.get("type")?.as_str()? {
        "scalar" => Some(MaterialProperty::Scalar(json_f64(fields.get("value")?)?)),
        "range" => Some(MaterialProperty::Range(
            json_f64(fields.get("lo")?)?,
            json_f64(fields.get("hi")?)?,
        )),
        "matrix" => {
            let values = fields
                .get("values")?
                .as_array()?
                .iter()
                .map(json_f64)
                .collect::<Option<Vec<f64>>>()?;
            Some(MaterialProperty::Matrix(values))
        }
        _ => None,
    }
}

/// A JSON number, or one of the strings `"NaN"`, `"Infinity"` and `"-Infinity"`.
fn json_f64(value: &serde_json::Value) -> Option<f64> {
    match value {
        serde_json::Value::Number(number) => number.as_f64(),
        serde_json::Value::String(text) => match text.as_str() {
            "NaN" => Some(f64::NAN),
            "Infinity" => Some(f64::INFINITY),
            "-Infinity" => Some(f64::NEG_INFINITY),
            _ => None,
        },
        _ => None,
    }
}

/// Render an `f64` as a JSON token (see [`JsonMaterialDb`] for the spelling of
/// non-finite values).
fn json_number(value: f64) -> String {
    if value.is_nan() {
        "\"NaN\"".to_string()
    } else if value == f64::INFINITY {
        "\"Infinity\"".to_string()
    } else if value == f64::NEG_INFINITY {
        "\"-Infinity\"".to_string()
    } else {
        // `{:?}` is the shortest form that round-trips, and valid JSON for every
        // finite value (never a leading `.`, a trailing `.` or a `+`).
        format!("{value:?}")
    }
}

/// Escape special characters for JSON strings.
///
/// Quotes, backslashes and all control characters are escaped, so the result is
/// valid JSON and decodes back to `s`.
fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000C}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

// ── CsvMaterialDb ─────────────────────────────────────────────────────────────

/// Matweb-style CSV material database reader.
///
/// Expected CSV format: first row is headers, subsequent rows are materials.
/// At minimum, columns `Name` and `Category` must be present.
/// All other columns are treated as scalar properties.
pub struct CsvMaterialDb;

impl CsvMaterialDb {
    /// Parse a Matweb-style CSV string into a `MaterialDatabase`.
    pub fn parse(csv: &str) -> MaterialDatabase {
        let mut db = MaterialDatabase::new();
        let mut lines = csv.lines();

        // Parse header row
        let header_line = match lines.next() {
            Some(h) => h,
            None => return db,
        };
        let headers: Vec<&str> = header_line.split(',').map(|s| s.trim()).collect();
        let name_idx = headers.iter().position(|&h| h.eq_ignore_ascii_case("Name"));
        let cat_idx = headers
            .iter()
            .position(|&h| h.eq_ignore_ascii_case("Category"));

        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            let name = name_idx
                .and_then(|i| fields.get(i))
                .copied()
                .unwrap_or("Unknown");
            let category = cat_idx
                .and_then(|i| fields.get(i))
                .copied()
                .unwrap_or("Uncategorized");
            let mut record = MaterialRecord::new(name, category);

            for (i, header) in headers.iter().enumerate() {
                if Some(i) == name_idx || Some(i) == cat_idx {
                    continue;
                }
                if let Some(val_str) = fields.get(i)
                    && let Ok(val) = val_str.parse::<f64>()
                {
                    record.set_scalar(header, val);
                }
            }
            db.add(record);
        }
        db
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod json_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn make_steel() -> MaterialRecord {
        let mut r = MaterialRecord::new("AISI 304", "Metal>Ferrous>Steel");
        r.set_scalar("density", 8000.0);
        r.set_scalar("elastic_modulus", 200e9);
        r.set_scalar("yield_strength", 215e6);
        r.set_scalar("thermal_conductivity", 16.0);
        r.source = "Matweb".to_string();
        r
    }

    fn make_aluminum() -> MaterialRecord {
        let mut r = MaterialRecord::new("Al 6061", "Metal>NonFerrous>Aluminum");
        r.set_scalar("density", 2700.0);
        r.set_scalar("elastic_modulus", 69e9);
        r.set_scalar("yield_strength", 276e6);
        r.set_scalar("thermal_conductivity", 167.0);
        r
    }

    // --- MaterialProperty ---

    #[test]
    fn test_scalar_representative() {
        let p = MaterialProperty::Scalar(42.0);
        assert!((p.representative() - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_range_representative_midpoint() {
        let p = MaterialProperty::Range(100.0, 200.0);
        assert!((p.representative() - 150.0).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_representative_first() {
        let p = MaterialProperty::Matrix(vec![7.0, 8.0, 9.0]);
        assert!((p.representative() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_contains() {
        let p = MaterialProperty::Scalar(5.0);
        assert!(p.contains(5.0));
        assert!(!p.contains(6.0));
    }

    #[test]
    fn test_range_contains() {
        let p = MaterialProperty::Range(10.0, 20.0);
        assert!(p.contains(15.0));
        assert!(!p.contains(25.0));
    }

    // --- MaterialRecord ---

    #[test]
    fn test_record_get_value() {
        let r = make_steel();
        let d = r.get_value("density").unwrap();
        assert!((d - 8000.0).abs() < 1e-6);
    }

    #[test]
    fn test_record_missing_key() {
        let r = make_steel();
        assert!(r.get_value("nonexistent").is_none());
    }

    // --- MaterialDatabase ---

    #[test]
    fn test_database_by_category() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        let metals = db.by_category("Metal>Ferrous");
        assert_eq!(metals.len(), 1);
        assert_eq!(metals[0].name, "AISI 304");
    }

    #[test]
    fn test_filter_by_property_range() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        // density between 2000 and 5000 → only aluminum
        let results = db.filter_by_property("density", 2000.0, 5000.0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].name, "Al 6061");
    }

    #[test]
    fn test_filter_returns_nothing_outside_range() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        let results = db.filter_by_property("density", 1000.0, 2000.0);
        assert!(results.is_empty());
    }

    #[test]
    fn test_ashby_chart_generates_pairs() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        let pairs = db.ashby_chart("density", "elastic_modulus");
        assert_eq!(pairs.len(), 2);
        for (x, y, _name) in &pairs {
            assert!(*x > 0.0);
            assert!(*y > 0.0);
        }
    }

    // --- JSON round-trip ---

    #[test]
    fn test_json_roundtrip_scalar() {
        let record = make_steel();
        let json = JsonMaterialDb::write_json_material(&record);
        let parsed = JsonMaterialDb::parse_json_material(&json).expect("parse failed");
        assert_eq!(parsed.name, record.name);
        assert_eq!(parsed.category, record.category);
        let d_orig = record.get_value("density").unwrap();
        let d_parsed = parsed.get_value("density").unwrap();
        assert!(
            (d_orig - d_parsed).abs() < 1e-3,
            "orig={d_orig} parsed={d_parsed}"
        );
    }

    #[test]
    fn test_json_roundtrip_source() {
        let record = make_steel();
        let json = JsonMaterialDb::write_json_material(&record);
        let parsed = JsonMaterialDb::parse_json_material(&json).unwrap();
        assert_eq!(parsed.source, "Matweb");
    }

    #[test]
    fn test_json_roundtrip_range_property() {
        let mut r = MaterialRecord::new("Ti-6Al-4V", "Metal>NonFerrous>Titanium");
        r.set_range("yield_strength", 800e6, 1000e6);
        let json = JsonMaterialDb::write_json_material(&r);
        let parsed = JsonMaterialDb::parse_json_material(&json).unwrap();
        let v = parsed.get_value("yield_strength").unwrap();
        assert!((v - 900e6).abs() < 1.0, "v={v}");
    }

    #[test]
    fn test_json_roundtrip_matrix_property() {
        let mut r = MaterialRecord::new("Composite", "Polymer>Composite");
        r.properties.insert(
            "stiffness_tensor".to_string(),
            MaterialProperty::Matrix(vec![1.0, 2.0, 3.0]),
        );
        let json = JsonMaterialDb::write_json_material(&r);
        let parsed = JsonMaterialDb::parse_json_material(&json).unwrap();
        if let Some(MaterialProperty::Matrix(v)) = parsed.properties.get("stiffness_tensor") {
            assert_eq!(v.len(), 3);
            assert!((v[0] - 1.0).abs() < 1e-4);
        } else {
            panic!("matrix property not found");
        }
    }

    #[test]
    fn test_json_roundtrip_non_ascii_strings() {
        let mut r = MaterialRecord::new("Stål 鋼", "Métal>鉄鋼");
        r.source = "Källa: 日本工業規格".to_string();
        r.set_scalar("density", 7850.0);
        let json = JsonMaterialDb::write_json_material(&r);
        let parsed = JsonMaterialDb::parse_json_material(&json).expect("parse failed");
        assert_eq!(parsed.name, "Stål 鋼");
        assert_eq!(parsed.category, "Métal>鉄鋼");
        assert_eq!(parsed.source, "Källa: 日本工業規格");
        assert!((parsed.get_value("density").unwrap() - 7850.0).abs() < 1e-3);
    }

    #[test]
    fn test_json_non_ascii_name_is_not_truncated_or_panicking() {
        // One 3-byte character: the char count (1) is not a byte offset into the text.
        let r = MaterialRecord::new("鋼", "Metal");
        let parsed =
            JsonMaterialDb::parse_json_material(&JsonMaterialDb::write_json_material(&r)).unwrap();
        assert_eq!(parsed.name, "鋼");
        // "Ståhl": 5 chars but 6 bytes; the old code cut it to "Ståh".
        let r = MaterialRecord::new("Ståhl", "Metal");
        let parsed =
            JsonMaterialDb::parse_json_material(&JsonMaterialDb::write_json_material(&r)).unwrap();
        assert_eq!(parsed.name, "Ståhl");
    }

    #[test]
    fn test_json_roundtrip_escapes_and_control_characters() {
        let name = "line1\nline2\t\"quoted\" back\\slash \u{1} \u{8}\u{c}\r end";
        let r = MaterialRecord::new(name, "a/b");
        let json = JsonMaterialDb::write_json_material(&r);
        assert!(
            !json.contains('\u{1}'),
            "control characters must be escaped"
        );
        assert!(json.contains("\\u0001"));
        let parsed = JsonMaterialDb::parse_json_material(&json).unwrap();
        assert_eq!(parsed.name, name);
        assert_eq!(parsed.category, "a/b");
    }

    // --- Performance index ---

    #[test]
    fn test_performance_index_stiffness_per_weight() {
        let cmp =
            MaterialComparison::new(vec!["elastic_modulus".to_string(), "density".to_string()]);
        let steel = make_steel();
        let al = make_aluminum();
        let idx_steel = cmp
            .performance_index(&steel, "elastic_modulus", "density")
            .unwrap();
        let idx_al = cmp
            .performance_index(&al, "elastic_modulus", "density")
            .unwrap();
        // Al should have higher specific stiffness (69e9/2700 ≈ 25.6 MN·m/kg)
        // vs steel (200e9/8000 = 25 MN·m/kg) — approximately equal but Al slightly higher
        assert!(idx_al > 0.0 && idx_steel > 0.0);
    }

    #[test]
    fn test_rank_by_index_order() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        let cmp = MaterialComparison::new(vec!["density".to_string()]);
        let refs: Vec<&MaterialRecord> = db.records.iter().collect();
        let ranked = cmp.rank_by_index(&refs, "elastic_modulus", "density");
        // Both should appear
        assert_eq!(ranked.len(), 2);
        // Highest index first
        assert!(ranked[0].1 >= ranked[1].1);
    }

    #[test]
    fn test_radar_data_normalized() {
        let cmp =
            MaterialComparison::new(vec!["density".to_string(), "elastic_modulus".to_string()]);
        let steel = make_steel();
        let radar = cmp.radar_data(&steel, &[8000.0, 200e9]);
        assert_eq!(radar.len(), 2);
        assert!((radar[0] - 1.0).abs() < 1e-9);
        assert!((radar[1] - 1.0).abs() < 1e-9);
    }

    // --- Temperature interpolation ---

    #[test]
    fn test_temperature_exact_at_table_point() {
        let td = TemperatureDependence::new(
            "yield_strength",
            vec![(20.0, 215e6), (200.0, 185e6), (400.0, 140e6)],
        );
        assert!((td.value_at(20.0) - 215e6).abs() < 1e-3);
        assert!((td.value_at(200.0) - 185e6).abs() < 1e-3);
        assert!((td.value_at(400.0) - 140e6).abs() < 1e-3);
    }

    #[test]
    fn test_temperature_interpolation_midpoint() {
        let td = TemperatureDependence::new("yield_strength", vec![(0.0, 0.0), (100.0, 100.0)]);
        let v = td.value_at(50.0);
        assert!((v - 50.0).abs() < 1e-9, "v={v}");
    }

    #[test]
    fn test_temperature_clamp_below_range() {
        let td = TemperatureDependence::new("E", vec![(100.0, 200e9), (300.0, 180e9)]);
        let v = td.value_at(0.0);
        assert!((v - 200e9).abs() < 1.0, "v={v}");
    }

    #[test]
    fn test_temperature_clamp_above_range() {
        let td = TemperatureDependence::new("E", vec![(100.0, 200e9), (300.0, 180e9)]);
        let v = td.value_at(500.0);
        assert!((v - 180e9).abs() < 1.0, "v={v}");
    }

    // --- MaterialFilter ---

    #[test]
    fn test_filter_category_prefix() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        let flt = MaterialFilter::new().filter_category("Metal>Ferrous");
        let res = flt.apply(&db);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].name, "AISI 304");
    }

    #[test]
    fn test_filter_min_property() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        // density >= 5000 → only steel
        let flt = MaterialFilter::new().filter_min("density", 5000.0);
        let res = flt.apply(&db);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].name, "AISI 304");
    }

    #[test]
    fn test_filter_max_property() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        // density <= 5000 → only aluminum
        let flt = MaterialFilter::new().filter_max("density", 5000.0);
        let res = flt.apply(&db);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].name, "Al 6061");
    }

    #[test]
    fn test_filter_combined() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        db.add(make_aluminum());
        // category Metal, density >= 5000 → steel only
        let flt = MaterialFilter::new()
            .filter_category("Metal")
            .filter_min("density", 5000.0);
        let res = flt.apply(&db);
        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_filter_nothing_passes_impossible_constraint() {
        let mut db = MaterialDatabase::new();
        db.add(make_steel());
        let flt = MaterialFilter::new().filter_min("density", 1e10);
        let res = flt.apply(&db);
        assert!(res.is_empty());
    }

    // --- Category hierarchy ---

    #[test]
    fn test_category_hierarchy_contains_leaf() {
        let tree = MaterialCategories::node(
            "Metal",
            vec![MaterialCategories::node(
                "Ferrous",
                vec![MaterialCategories::leaf("Steel")],
            )],
        );
        assert!(tree.contains_name("Steel"));
        assert!(tree.contains_name("Ferrous"));
        assert!(tree.contains_name("Metal"));
        assert!(!tree.contains_name("Polymer"));
    }

    #[test]
    fn test_category_leaf_names() {
        let tree = MaterialCategories::node(
            "Metal",
            vec![
                MaterialCategories::node(
                    "Ferrous",
                    vec![
                        MaterialCategories::leaf("Steel"),
                        MaterialCategories::leaf("Cast Iron"),
                    ],
                ),
                MaterialCategories::leaf("Aluminum"),
            ],
        );
        let leaves = tree.leaf_names();
        assert_eq!(leaves.len(), 3);
        assert!(leaves.contains(&"Steel".to_string()));
        assert!(leaves.contains(&"Cast Iron".to_string()));
        assert!(leaves.contains(&"Aluminum".to_string()));
    }

    // --- CSV parse ---

    #[test]
    fn test_csv_parse_basic() {
        let csv = "Name,Category,density,elastic_modulus\n\
                   Steel,Metal>Ferrous,7850,210000000000\n\
                   Aluminum,Metal>NonFerrous,2700,69000000000\n";
        let db = CsvMaterialDb::parse(csv);
        assert_eq!(db.records.len(), 2);
        assert_eq!(db.records[0].name, "Steel");
        let d = db.records[0].get_value("density").unwrap();
        assert!((d - 7850.0).abs() < 1e-3, "density={d}");
    }

    #[test]
    fn test_csv_category_mapping() {
        let csv = "Name,Category,density\nAluminum,Metal>NonFerrous,2700\n";
        let db = CsvMaterialDb::parse(csv);
        assert_eq!(db.records[0].category, "Metal>NonFerrous");
    }

    #[test]
    fn test_csv_empty_input() {
        let db = CsvMaterialDb::parse("");
        assert!(db.records.is_empty());
    }

    #[test]
    fn test_csv_skips_empty_lines() {
        let csv = "Name,Category,density\n\nSteel,Metal,7850\n\n";
        let db = CsvMaterialDb::parse(csv);
        assert_eq!(db.records.len(), 1);
    }
}
