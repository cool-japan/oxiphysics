// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Runtime material interaction table.
//!
//! The `MaterialTable` maps `MaterialId`s to `MaterialDef`s and resolves
//! contact-pair properties (restitution, static and dynamic friction) through
//! configurable `CombineRule`s.  Per-pair overrides let you fine-tune the
//! contact response for specific surface combinations without changing either
//! material globally.
//!
//! > **Scope note**: This module covers *contact physics* — the coefficients
//! > consumed by the collision solver when two surfaces meet.  It is distinct
//! > from the `oxiphysics_materials` crate, which carries bulk physical-property
//! > databases (composites, polymers, aerospace alloys, electrochemistry, etc.).
//!
//! ## Built-in presets
//!
//! The table is pre-populated with six common surfaces:
//! `concrete`, `rubber`, `metal`, `ice`, `wood`, `glass`.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::material_table::{MaterialTable, CombineRule};
//!
//! let mut table = MaterialTable::new();
//! let rubber   = table.id_for_name("rubber").unwrap();
//! let concrete = table.id_for_name("concrete").unwrap();
//!
//! // Rubber-on-concrete: high friction
//! let friction = table.contact_static_friction(rubber, concrete);
//! assert!(friction > 0.5);
//!
//! // Ice-on-ice: override bounce to a specific value
//! let ice = table.id_for_name("ice").unwrap();
//! table.set_pair_restitution(ice, ice, 0.05);
//! assert!((table.contact_restitution(ice, ice) - 0.05).abs() < 1e-9);
//! ```

use std::collections::HashMap;

// ============================================================================
// MaterialId
// ============================================================================

/// Lightweight, copyable handle to a registered material.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct MaterialId(pub u32);

impl std::fmt::Display for MaterialId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MatId({})", self.0)
    }
}

// ============================================================================
// MaterialDef
// ============================================================================

/// Physical surface properties of a single material.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MaterialDef {
    /// Human-readable name (used by [`MaterialTable::id_for_name`]).
    pub name: String,
    /// Bulk density (kg/m³).
    pub density: f64,
    /// Coefficient of restitution (bounciness) in `[0, 1]`.
    pub restitution: f64,
    /// Coefficient of static friction.
    pub static_friction: f64,
    /// Coefficient of dynamic (kinetic) friction.
    pub dynamic_friction: f64,
    /// Linear velocity damping per second.
    pub linear_damping: f64,
    /// Angular velocity damping per second.
    pub angular_damping: f64,
}

impl MaterialDef {
    /// Creates a fully specified material definition.
    pub fn new(
        name: impl Into<String>,
        density: f64,
        restitution: f64,
        static_friction: f64,
        dynamic_friction: f64,
        linear_damping: f64,
        angular_damping: f64,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            restitution,
            static_friction,
            dynamic_friction,
            linear_damping,
            angular_damping,
        }
    }
}

// ============================================================================
// CombineRule
// ============================================================================

/// Rule used to combine two material coefficients into a single contact value.
///
/// Applied when no per-pair override is set.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CombineRule {
    /// Use the smaller of the two values.
    Min,
    /// Use the larger of the two values.
    Max,
    /// Arithmetic mean: `(a + b) / 2`.
    #[default]
    Average,
    /// Geometric mean: `sqrt(a · b)`.
    GeometricMean,
    /// Product: `a · b`.
    Multiply,
}

impl CombineRule {
    /// Applies the rule to two values.
    pub fn combine(self, a: f64, b: f64) -> f64 {
        match self {
            CombineRule::Min => a.min(b),
            CombineRule::Max => a.max(b),
            CombineRule::Average => (a + b) * 0.5,
            CombineRule::GeometricMean => (a * b).max(0.0).sqrt(),
            CombineRule::Multiply => a * b,
        }
    }
}

// ============================================================================
// PairKey (canonicalised so (a,b) == (b,a))
// ============================================================================

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct PairKey(u32, u32);

impl PairKey {
    fn new(a: MaterialId, b: MaterialId) -> Self {
        let lo = a.0.min(b.0);
        let hi = a.0.max(b.0);
        PairKey(lo, hi)
    }
}

// ============================================================================
// PairOverride (internal)
// ============================================================================

#[derive(Clone, Debug, Default)]
struct PairOverride {
    restitution: Option<f64>,
    static_friction: Option<f64>,
    dynamic_friction: Option<f64>,
    /// Per-pair combine rule; overrides table defaults when set.
    combine_rule: Option<CombineRule>,
}

// ============================================================================
// MaterialTable
// ============================================================================

/// Registry of material definitions with contact-pair resolution.
///
/// Pre-populated with built-in presets.  New materials can be registered at
/// runtime and pair-specific overrides applied independently.
#[derive(Debug)]
pub struct MaterialTable {
    materials: Vec<MaterialDef>,
    name_map: HashMap<String, MaterialId>,
    pair_overrides: HashMap<PairKey, PairOverride>,
    /// Default [`CombineRule`] for resolving restitution.
    pub restitution_rule: CombineRule,
    /// Default [`CombineRule`] for resolving friction coefficients.
    pub friction_rule: CombineRule,
}

impl Default for MaterialTable {
    fn default() -> Self {
        Self::new()
    }
}

impl MaterialTable {
    /// Creates a table pre-populated with six built-in surface presets:
    /// `concrete`, `rubber`, `metal`, `ice`, `wood`, `glass`.
    pub fn new() -> Self {
        let mut table = Self {
            materials: Vec::new(),
            name_map: HashMap::new(),
            pair_overrides: HashMap::new(),
            restitution_rule: CombineRule::Average,
            friction_rule: CombineRule::Average,
        };
        //                         name       density  rest   sf     df     ld      ad
        table.register(MaterialDef::new(
            "concrete", 2_300.0, 0.20, 0.70, 0.60, 0.000, 0.000,
        ));
        table.register(MaterialDef::new(
            "rubber", 1_100.0, 0.65, 0.85, 0.75, 0.010, 0.010,
        ));
        table.register(MaterialDef::new(
            "metal", 7_850.0, 0.30, 0.50, 0.45, 0.000, 0.000,
        ));
        table.register(MaterialDef::new(
            "ice", 917.0, 0.10, 0.05, 0.03, 0.000, 0.000,
        ));
        table.register(MaterialDef::new(
            "wood", 700.0, 0.40, 0.55, 0.45, 0.005, 0.005,
        ));
        table.register(MaterialDef::new(
            "glass", 2_500.0, 0.50, 0.40, 0.35, 0.000, 0.000,
        ));
        table
    }

    // -----------------------------------------------------------------------
    // Registration
    // -----------------------------------------------------------------------

    /// Registers a material and returns its [`MaterialId`].
    ///
    /// If a material with the same name already exists, its definition is
    /// **replaced** and the same id is returned.
    pub fn register(&mut self, def: MaterialDef) -> MaterialId {
        if let Some(&id) = self.name_map.get(&def.name) {
            self.materials[id.0 as usize] = def;
            return id;
        }
        let id = MaterialId(self.materials.len() as u32);
        self.name_map.insert(def.name.clone(), id);
        self.materials.push(def);
        id
    }

    /// Convenience wrapper for registering by individual fields.
    pub fn register_named(
        &mut self,
        name: impl Into<String>,
        density: f64,
        restitution: f64,
        static_friction: f64,
        dynamic_friction: f64,
        linear_damping: f64,
        angular_damping: f64,
    ) -> MaterialId {
        self.register(MaterialDef::new(
            name,
            density,
            restitution,
            static_friction,
            dynamic_friction,
            linear_damping,
            angular_damping,
        ))
    }

    // -----------------------------------------------------------------------
    // Lookup
    // -----------------------------------------------------------------------

    /// Returns the definition for `id`, or `None` if unregistered.
    pub fn lookup(&self, id: MaterialId) -> Option<&MaterialDef> {
        self.materials.get(id.0 as usize)
    }

    /// Returns the id for a material name, or `None` if not found.
    pub fn id_for_name(&self, name: &str) -> Option<MaterialId> {
        self.name_map.get(name).copied()
    }

    /// Number of registered materials.
    pub fn material_count(&self) -> usize {
        self.materials.len()
    }

    /// Number of active per-pair overrides.
    pub fn pair_count(&self) -> usize {
        self.pair_overrides.len()
    }

    /// Iterates all registered material definitions.
    pub fn iter_materials(&self) -> impl Iterator<Item = &MaterialDef> {
        self.materials.iter()
    }

    // -----------------------------------------------------------------------
    // Pair overrides
    // -----------------------------------------------------------------------

    /// Sets a per-pair restitution override (bypasses combine rule).
    pub fn set_pair_restitution(&mut self, a: MaterialId, b: MaterialId, value: f64) {
        self.pair_overrides
            .entry(PairKey::new(a, b))
            .or_default()
            .restitution = Some(value);
    }

    /// Sets a per-pair static-friction override.
    pub fn set_pair_static_friction(&mut self, a: MaterialId, b: MaterialId, value: f64) {
        self.pair_overrides
            .entry(PairKey::new(a, b))
            .or_default()
            .static_friction = Some(value);
    }

    /// Sets a per-pair dynamic-friction override.
    pub fn set_pair_dynamic_friction(&mut self, a: MaterialId, b: MaterialId, value: f64) {
        self.pair_overrides
            .entry(PairKey::new(a, b))
            .or_default()
            .dynamic_friction = Some(value);
    }

    /// Sets a per-pair [`CombineRule`] (overrides [`restitution_rule`](Self::restitution_rule)
    /// and [`friction_rule`](Self::friction_rule) for this pair).
    pub fn set_pair_combine_rule(&mut self, a: MaterialId, b: MaterialId, rule: CombineRule) {
        self.pair_overrides
            .entry(PairKey::new(a, b))
            .or_default()
            .combine_rule = Some(rule);
    }

    /// Removes all overrides for a material pair.
    pub fn clear_pair_override(&mut self, a: MaterialId, b: MaterialId) {
        self.pair_overrides.remove(&PairKey::new(a, b));
    }

    // -----------------------------------------------------------------------
    // Contact property resolution
    // -----------------------------------------------------------------------

    /// Resolves the contact restitution for a material pair.
    ///
    /// Precedence:
    /// 1. Per-pair restitution override (if set).
    /// 2. Per-pair [`CombineRule`] applied to both materials' coefficients.
    /// 3. Table-wide [`restitution_rule`](Self::restitution_rule).
    pub fn contact_restitution(&self, a: MaterialId, b: MaterialId) -> f64 {
        let key = PairKey::new(a, b);
        if let Some(r) = self.pair_overrides.get(&key).and_then(|ov| ov.restitution) {
            return r;
        }
        let ra = self.lookup(a).map_or(0.0, |m| m.restitution);
        let rb = self.lookup(b).map_or(0.0, |m| m.restitution);
        let rule = self
            .pair_overrides
            .get(&key)
            .and_then(|ov| ov.combine_rule)
            .unwrap_or(self.restitution_rule);
        rule.combine(ra, rb)
    }

    /// Resolves the contact static friction for a material pair.
    pub fn contact_static_friction(&self, a: MaterialId, b: MaterialId) -> f64 {
        let key = PairKey::new(a, b);
        if let Some(f) = self
            .pair_overrides
            .get(&key)
            .and_then(|ov| ov.static_friction)
        {
            return f;
        }
        let fa = self.lookup(a).map_or(0.5, |m| m.static_friction);
        let fb = self.lookup(b).map_or(0.5, |m| m.static_friction);
        let rule = self
            .pair_overrides
            .get(&key)
            .and_then(|ov| ov.combine_rule)
            .unwrap_or(self.friction_rule);
        rule.combine(fa, fb)
    }

    /// Resolves the contact dynamic friction for a material pair.
    pub fn contact_dynamic_friction(&self, a: MaterialId, b: MaterialId) -> f64 {
        let key = PairKey::new(a, b);
        if let Some(f) = self
            .pair_overrides
            .get(&key)
            .and_then(|ov| ov.dynamic_friction)
        {
            return f;
        }
        let fa = self.lookup(a).map_or(0.4, |m| m.dynamic_friction);
        let fb = self.lookup(b).map_or(0.4, |m| m.dynamic_friction);
        let rule = self
            .pair_overrides
            .get(&key)
            .and_then(|ov| ov.combine_rule)
            .unwrap_or(self.friction_rule);
        rule.combine(fa, fb)
    }

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------

    /// Returns a multiline summary listing all registered materials.
    pub fn summary(&self) -> String {
        let mut s = format!(
            "MaterialTable: {} materials, {} pair overrides\n",
            self.material_count(),
            self.pair_count(),
        );
        for (i, m) in self.materials.iter().enumerate() {
            s.push_str(&format!(
                "  [{i:2}] {:12} | density={:6.0} kg/m³ | rest={:.2} | sf={:.2} df={:.2}\n",
                m.name, m.density, m.restitution, m.static_friction, m.dynamic_friction
            ));
        }
        s
    }
}

impl std::fmt::Display for MaterialTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "MaterialTable {{ {} materials, {} pair overrides }}",
            self.material_count(),
            self.pair_count(),
        )
    }
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_registered() {
        let table = MaterialTable::new();
        assert_eq!(table.material_count(), 6);
        for name in &["concrete", "rubber", "metal", "ice", "wood", "glass"] {
            assert!(table.id_for_name(name).is_some(), "missing preset: {name}");
        }
    }

    #[test]
    fn pair_restitution_override() {
        let mut table = MaterialTable::new();
        let ice = table.id_for_name("ice").unwrap();
        table.set_pair_restitution(ice, ice, 0.05);
        assert!((table.contact_restitution(ice, ice) - 0.05).abs() < 1e-9);
        // Symmetric
        assert!((table.contact_restitution(ice, ice) - 0.05).abs() < 1e-9);
    }

    #[test]
    fn friction_average_combine() {
        let table = MaterialTable::new();
        let concrete = table.id_for_name("concrete").unwrap();
        let ice = table.id_for_name("ice").unwrap();
        // concrete sf = 0.70, ice sf = 0.05 → average = 0.375
        let expected = (0.70 + 0.05) * 0.5;
        let got = table.contact_static_friction(concrete, ice);
        assert!((got - expected).abs() < 1e-9);
    }

    #[test]
    fn register_duplicate_replaces() {
        let mut table = MaterialTable::new();
        let id1 = table.id_for_name("concrete").unwrap();
        let id2 = table.register(MaterialDef::new(
            "concrete", 3000.0, 0.99, 0.99, 0.99, 0.0, 0.0,
        ));
        assert_eq!(id1, id2); // same id
        assert!((table.lookup(id1).unwrap().density - 3000.0).abs() < 1e-9);
    }

    #[test]
    fn combine_rules() {
        assert!((CombineRule::Min.combine(0.3, 0.7) - 0.3).abs() < 1e-9);
        assert!((CombineRule::Max.combine(0.3, 0.7) - 0.7).abs() < 1e-9);
        assert!((CombineRule::Average.combine(0.3, 0.7) - 0.5).abs() < 1e-9);
        assert!((CombineRule::Multiply.combine(0.5, 0.4) - 0.2).abs() < 1e-9);
        let gm = CombineRule::GeometricMean.combine(4.0, 9.0);
        assert!((gm - 6.0).abs() < 1e-9);
    }

    #[test]
    fn clear_override() {
        let mut table = MaterialTable::new();
        let rubber = table.id_for_name("rubber").unwrap();
        let concrete = table.id_for_name("concrete").unwrap();
        table.set_pair_static_friction(rubber, concrete, 0.01);
        assert!((table.contact_static_friction(rubber, concrete) - 0.01).abs() < 1e-9);
        table.clear_pair_override(rubber, concrete);
        // Back to default average
        let expected = (0.85 + 0.70) * 0.5;
        let got = table.contact_static_friction(rubber, concrete);
        assert!((got - expected).abs() < 1e-9);
    }
}
