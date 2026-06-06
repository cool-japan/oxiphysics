// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Material properties database for physics simulations.
//!
//! Provides a searchable in-memory database of engineering materials with
//! mechanical, thermal, and physical properties. Supports JSON serialisation.

/// A single material entry in the database.
#[derive(Debug, Clone)]
pub struct MaterialEntry {
    /// Material name (e.g. `"steel_1020"`).
    pub name: String,
    /// Mass density in kg/m³.
    pub density: f64,
    /// Young's modulus in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Thermal conductivity in W/(m·K).
    pub thermal_conductivity: f64,
    /// Specific heat capacity in J/(kg·K).
    pub specific_heat: f64,
    /// Coefficient of thermal expansion in 1/K.
    pub thermal_expansion: f64,
}

impl MaterialEntry {
    /// Create a new material entry with all fields.
    pub fn new(
        name: impl Into<String>,
        density: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
        thermal_conductivity: f64,
        specific_heat: f64,
        thermal_expansion: f64,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            youngs_modulus,
            poisson_ratio,
            thermal_conductivity,
            specific_heat,
            thermal_expansion,
        }
    }

    /// Serialise this entry to a JSON-like string fragment (no outer braces).
    fn to_json_fields(&self) -> String {
        format!(
            r#""name":"{name}","density":{density},"youngs_modulus":{ym},"poisson_ratio":{pr},"thermal_conductivity":{tc},"specific_heat":{sh},"thermal_expansion":{te}"#,
            name = self.name,
            density = self.density,
            ym = self.youngs_modulus,
            pr = self.poisson_ratio,
            tc = self.thermal_conductivity,
            sh = self.specific_heat,
            te = self.thermal_expansion,
        )
    }
}

/// An in-memory database of material entries.
#[derive(Debug, Default)]
pub struct MaterialDatabase {
    entries: Vec<MaterialEntry>,
}

impl MaterialDatabase {
    /// Create an empty material database.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a database pre-populated with common engineering materials.
    ///
    /// Includes: steel 1020, aluminium 6061, copper, titanium Ti-6Al-4V,
    /// PTFE, water (liquid, 20 °C), and dry air (20 °C, 1 atm).
    pub fn with_defaults() -> Self {
        let mut db = Self::new();
        // Steel 1020
        db.add_material(MaterialEntry::new(
            "steel_1020",
            7_870.0,
            200.0e9,
            0.29,
            51.9,
            486.0,
            11.7e-6,
        ));
        // Aluminium 6061-T6
        db.add_material(MaterialEntry::new(
            "aluminum_6061",
            2_700.0,
            68.9e9,
            0.33,
            167.0,
            896.0,
            23.6e-6,
        ));
        // Copper (pure, annealed)
        db.add_material(MaterialEntry::new(
            "copper", 8_960.0, 117.0e9, 0.34, 385.0, 385.0, 17.0e-6,
        ));
        // Titanium Ti-6Al-4V
        db.add_material(MaterialEntry::new(
            "titanium_ti6al4v",
            4_430.0,
            113.8e9,
            0.342,
            6.7,
            526.3,
            8.6e-6,
        ));
        // PTFE (Teflon)
        db.add_material(MaterialEntry::new(
            "ptfe", 2_200.0, 0.5e9, 0.46, 0.25, 1_004.0, 135.0e-6,
        ));
        // Water (liquid, 20 °C)
        db.add_material(MaterialEntry::new(
            "water", 998.2, 2.2e9, 0.5, 0.598, 4_182.0, 0.207e-3,
        ));
        // Dry air (20 °C, 1 atm)
        db.add_material(MaterialEntry::new(
            "air", 1.204, 0.0, 0.0, 0.0257, 1_005.0, 3.43e-3,
        ));
        db
    }

    /// Add a material to the database.
    pub fn add_material(&mut self, entry: MaterialEntry) {
        self.entries.push(entry);
    }

    /// Look up a material by name (case-sensitive).
    ///
    /// Returns `Some(&MaterialEntry)` if found, `None` otherwise.
    pub fn get_material(&self, name: &str) -> Option<&MaterialEntry> {
        self.entries.iter().find(|e| e.name == name)
    }

    /// Remove a material by name.
    ///
    /// Returns `true` if a material was removed, `false` if not found.
    pub fn remove_material(&mut self, name: &str) -> bool {
        let before = self.entries.len();
        self.entries.retain(|e| e.name != name);
        self.entries.len() < before
    }

    /// Returns the number of entries in the database.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the database has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return all materials whose density falls within `[min_density, max_density]`.
    pub fn search_by_density(&self, min_density: f64, max_density: f64) -> Vec<&MaterialEntry> {
        self.entries
            .iter()
            .filter(|e| e.density >= min_density && e.density <= max_density)
            .collect()
    }

    /// Return all materials whose Young's modulus falls within `[min_e, max_e]`.
    pub fn search_by_youngs_modulus(&self, min_e: f64, max_e: f64) -> Vec<&MaterialEntry> {
        self.entries
            .iter()
            .filter(|e| e.youngs_modulus >= min_e && e.youngs_modulus <= max_e)
            .collect()
    }

    /// Serialise the entire database to a JSON string.
    ///
    /// The format is a JSON array of objects: `[{"name":..., ...}, ...]`.
    pub fn export_json(&self) -> String {
        let items: Vec<String> = self
            .entries
            .iter()
            .map(|e| format!("{{{}}}", e.to_json_fields()))
            .collect();
        format!("[{}]", items.join(","))
    }

    /// Deserialise a database from a JSON string produced by `export_json`.
    ///
    /// Returns `Err` if parsing fails. This is a minimal hand-rolled parser that
    /// handles the exact format produced by `export_json`.
    pub fn import_json(json: &str) -> Result<Self, String> {
        let mut db = Self::new();
        let trimmed = json.trim();
        if trimmed == "[]" {
            return Ok(db);
        }
        // Strip outer brackets
        let inner = trimmed
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .ok_or("Expected JSON array")?;

        // Split objects: find balanced { } blocks
        let objects = split_json_objects(inner)?;
        for obj in objects {
            let entry = parse_material_json_object(&obj)?;
            db.add_material(entry);
        }
        Ok(db)
    }
}

// ── Private JSON helpers ──────────────────────────────────────────────────────

/// Split a comma-delimited sequence of `{...}` objects into individual strings.
fn split_json_objects(s: &str) -> Result<Vec<String>, String> {
    let mut objects = Vec::new();
    let mut depth = 0i32;
    let mut start = 0usize;

    for (i, c) in s.char_indices() {
        match c {
            '{' => {
                if depth == 0 {
                    start = i;
                }
                depth += 1;
            }
            '}' => {
                depth -= 1;
                if depth == 0 {
                    objects.push(s[start..=i].to_string());
                }
            }
            _ => {}
        }
    }

    if depth != 0 {
        return Err("Unbalanced braces in JSON".into());
    }
    Ok(objects)
}

/// Extract the value of a JSON string field `"key":"value"`.
fn extract_str_field<'a>(obj: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{}\":\"", key);
    let start = obj.find(&needle)? + needle.len();
    let end = obj[start..].find('"')? + start;
    Some(&obj[start..end])
}

/// Extract the value of a JSON numeric field `"key":value`.
fn extract_f64_field(obj: &str, key: &str) -> Option<f64> {
    let needle = format!("\"{}\":", key);
    let start = obj.find(&needle)? + needle.len();
    let rest = &obj[start..];
    // Read until comma, closing brace, or end
    let end = rest.find([',', '}']).unwrap_or(rest.len());
    rest[..end].trim().parse().ok()
}

/// Parse a single JSON object `{...}` into a `MaterialEntry`.
fn parse_material_json_object(obj: &str) -> Result<MaterialEntry, String> {
    let name = extract_str_field(obj, "name")
        .ok_or("missing name")?
        .to_string();
    let density = extract_f64_field(obj, "density").ok_or("missing density")?;
    let youngs_modulus =
        extract_f64_field(obj, "youngs_modulus").ok_or("missing youngs_modulus")?;
    let poisson_ratio = extract_f64_field(obj, "poisson_ratio").ok_or("missing poisson_ratio")?;
    let thermal_conductivity =
        extract_f64_field(obj, "thermal_conductivity").ok_or("missing thermal_conductivity")?;
    let specific_heat = extract_f64_field(obj, "specific_heat").ok_or("missing specific_heat")?;
    let thermal_expansion =
        extract_f64_field(obj, "thermal_expansion").ok_or("missing thermal_expansion")?;

    Ok(MaterialEntry {
        name,
        density,
        youngs_modulus,
        poisson_ratio,
        thermal_conductivity,
        specific_heat,
        thermal_expansion,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_steel() -> MaterialEntry {
        MaterialEntry::new("steel_test", 7_870.0, 200.0e9, 0.29, 51.9, 486.0, 11.7e-6)
    }

    // add / get
    #[test]
    fn test_add_and_get() {
        let mut db = MaterialDatabase::new();
        db.add_material(make_steel());
        assert!(db.get_material("steel_test").is_some());
    }

    #[test]
    fn test_get_missing() {
        let db = MaterialDatabase::new();
        assert!(db.get_material("nonexistent").is_none());
    }

    #[test]
    fn test_len_empty() {
        let db = MaterialDatabase::new();
        assert_eq!(db.len(), 0);
        assert!(db.is_empty());
    }

    #[test]
    fn test_len_after_add() {
        let mut db = MaterialDatabase::new();
        db.add_material(make_steel());
        assert_eq!(db.len(), 1);
        assert!(!db.is_empty());
    }

    // remove
    #[test]
    fn test_remove_existing() {
        let mut db = MaterialDatabase::new();
        db.add_material(make_steel());
        assert!(db.remove_material("steel_test"));
        assert!(db.is_empty());
    }

    #[test]
    fn test_remove_missing() {
        let mut db = MaterialDatabase::new();
        assert!(!db.remove_material("ghost"));
    }

    // defaults
    #[test]
    fn test_defaults_count() {
        let db = MaterialDatabase::with_defaults();
        assert_eq!(db.len(), 7);
    }

    #[test]
    fn test_defaults_steel_exists() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.get_material("steel_1020").is_some());
    }

    #[test]
    fn test_defaults_aluminum_density() {
        let db = MaterialDatabase::with_defaults();
        let al = db.get_material("aluminum_6061").unwrap();
        assert!((al.density - 2_700.0).abs() < 1.0);
    }

    #[test]
    fn test_defaults_copper_youngs() {
        let db = MaterialDatabase::with_defaults();
        let cu = db.get_material("copper").unwrap();
        assert!((cu.youngs_modulus - 117.0e9).abs() < 1e6);
    }

    #[test]
    fn test_defaults_titanium_exists() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.get_material("titanium_ti6al4v").is_some());
    }

    #[test]
    fn test_defaults_ptfe_exists() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.get_material("ptfe").is_some());
    }

    #[test]
    fn test_defaults_water_exists() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.get_material("water").is_some());
    }

    #[test]
    fn test_defaults_air_exists() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.get_material("air").is_some());
    }

    // search_by_density
    #[test]
    fn test_search_by_density_finds_metals() {
        let db = MaterialDatabase::with_defaults();
        let results = db.search_by_density(2_000.0, 9_000.0);
        assert!(results.len() >= 4); // Al, Cu, Ti, steel
    }

    #[test]
    fn test_search_by_density_excludes_air() {
        let db = MaterialDatabase::with_defaults();
        let results = db.search_by_density(2_000.0, 9_000.0);
        assert!(!results.iter().any(|e| e.name == "air"));
    }

    #[test]
    fn test_search_by_density_empty_range() {
        let db = MaterialDatabase::with_defaults();
        let results = db.search_by_density(1e12, 2e12);
        assert!(results.is_empty());
    }

    // search_by_youngs_modulus
    #[test]
    fn test_search_by_youngs_excludes_air() {
        let db = MaterialDatabase::with_defaults();
        let results = db.search_by_youngs_modulus(1.0e9, 300.0e9);
        assert!(!results.iter().any(|e| e.name == "air"));
    }

    #[test]
    fn test_search_by_youngs_finds_steel() {
        let db = MaterialDatabase::with_defaults();
        let results = db.search_by_youngs_modulus(150.0e9, 250.0e9);
        assert!(results.iter().any(|e| e.name == "steel_1020"));
    }

    // JSON roundtrip
    #[test]
    fn test_export_import_roundtrip() {
        let db = MaterialDatabase::with_defaults();
        let json = db.export_json();
        let db2 = MaterialDatabase::import_json(&json).unwrap();
        assert_eq!(db2.len(), db.len());
    }

    #[test]
    fn test_export_import_preserves_density() {
        let db = MaterialDatabase::with_defaults();
        let json = db.export_json();
        let db2 = MaterialDatabase::import_json(&json).unwrap();
        let orig = db.get_material("steel_1020").unwrap();
        let restored = db2.get_material("steel_1020").unwrap();
        assert!((orig.density - restored.density).abs() < 1e-3);
    }

    #[test]
    fn test_export_import_preserves_youngs_modulus() {
        let db = MaterialDatabase::with_defaults();
        let json = db.export_json();
        let db2 = MaterialDatabase::import_json(&json).unwrap();
        let orig = db.get_material("copper").unwrap();
        let restored = db2.get_material("copper").unwrap();
        assert!((orig.youngs_modulus - restored.youngs_modulus).abs() < 1e4);
    }

    #[test]
    fn test_export_import_empty_db() {
        let db = MaterialDatabase::new();
        let json = db.export_json();
        let db2 = MaterialDatabase::import_json(&json).unwrap();
        assert!(db2.is_empty());
    }

    #[test]
    fn test_export_import_single_entry() {
        let mut db = MaterialDatabase::new();
        db.add_material(make_steel());
        let json = db.export_json();
        let db2 = MaterialDatabase::import_json(&json).unwrap();
        assert_eq!(db2.len(), 1);
        let entry = db2.get_material("steel_test").unwrap();
        assert!((entry.poisson_ratio - 0.29).abs() < 1e-9);
    }

    #[test]
    fn test_import_invalid_json() {
        assert!(MaterialDatabase::import_json("{bad json}").is_err());
    }

    // material_entry fields
    #[test]
    fn test_material_entry_fields() {
        let entry = make_steel();
        assert_eq!(entry.name, "steel_test");
        assert!((entry.density - 7870.0).abs() < 1.0);
        assert!((entry.youngs_modulus - 200.0e9).abs() < 1e6);
        assert!((entry.poisson_ratio - 0.29).abs() < 1e-9);
        assert!((entry.thermal_conductivity - 51.9).abs() < 0.01);
        assert!((entry.specific_heat - 486.0).abs() < 0.1);
        assert!((entry.thermal_expansion - 11.7e-6).abs() < 1e-10);
    }
}
