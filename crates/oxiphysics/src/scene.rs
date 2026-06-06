// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Declarative scene description with JSON round-trip and a fluent builder.
//!
//! A `SceneDescription` is an engine-agnostic specification of a physics
//! world: gravity, simulation parameters, and a list of rigid bodies with
//! shapes, materials, and optional constraints.  It can be serialized to JSON
//! and deserialized back, making it suitable for save files, level editors,
//! and test fixtures.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::scene::{SceneBuilder, SceneMaterial};
//!
//! let scene = SceneBuilder::new("sandbox")
//!     .with_gravity([0.0, -9.81, 0.0])
//!     .with_dt(1.0 / 60.0)
//!     .add_static_box("floor", [0.0, -0.5, 0.0], [10.0, 0.5, 10.0])
//!     .add_sphere("ball", [0.0, 5.0, 0.0], 0.5, false)
//!     .build();
//!
//! assert_eq!(scene.body_count(), 2);
//! let json = scene.to_json().unwrap();
//! let loaded = oxiphysics::scene::SceneDescription::from_json(&json).unwrap();
//! assert_eq!(loaded.name, "sandbox");
//! ```

use serde::{Deserialize, Serialize};

// ============================================================================
// SceneShape
// ============================================================================

/// Collision shape for a scene body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SceneShape {
    /// Sphere with a given radius.
    Sphere {
        /// Sphere radius in metres.
        radius: f64,
    },
    /// Axis-aligned box defined by half-extents.
    Box {
        /// Half-extent in each axis (x, y, z).
        half_extents: [f64; 3],
    },
    /// Capsule (cylinder + two hemisphere end-caps).
    Capsule {
        /// Capsule radius.
        radius: f64,
        /// Half-height of the cylinder portion (excluding end-caps).
        half_height: f64,
    },
    /// Cylinder.
    Cylinder {
        /// Cylinder radius.
        radius: f64,
        /// Half-height.
        half_height: f64,
    },
    /// Infinite half-space defined by an outward normal.
    Plane {
        /// Outward surface normal (should be unit length).
        normal: [f64; 3],
    },
}

impl SceneShape {
    /// Human-readable shape name for diagnostics.
    pub fn kind_name(&self) -> &'static str {
        match self {
            SceneShape::Sphere { .. } => "sphere",
            SceneShape::Box { .. } => "box",
            SceneShape::Capsule { .. } => "capsule",
            SceneShape::Cylinder { .. } => "cylinder",
            SceneShape::Plane { .. } => "plane",
        }
    }

    /// Approximate volume in m³.
    pub fn approx_volume(&self) -> f64 {
        use std::f64::consts::PI;
        match self {
            SceneShape::Sphere { radius } => (4.0 / 3.0) * PI * radius.powi(3),
            SceneShape::Box { half_extents } => {
                8.0 * half_extents[0] * half_extents[1] * half_extents[2]
            }
            SceneShape::Capsule {
                radius,
                half_height,
            } => PI * radius * radius * (half_height * 2.0 + (4.0 / 3.0) * radius),
            SceneShape::Cylinder {
                radius,
                half_height,
            } => PI * radius * radius * half_height * 2.0,
            SceneShape::Plane { .. } => f64::INFINITY,
        }
    }
}

// ============================================================================
// SceneMaterial
// ============================================================================

/// Physical material properties for a scene body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneMaterial {
    /// Mass density in kg/m³ (default 1000 — water).
    pub density: f64,
    /// Coefficient of restitution [0, 1].
    pub restitution: f64,
    /// Coefficient of static friction.
    pub friction: f64,
    /// Coefficient of dynamic friction (defaults to `friction` if not set).
    #[serde(default)]
    pub dynamic_friction: Option<f64>,
    /// Linear damping coefficient (0 = no damping).
    #[serde(default)]
    pub linear_damping: f64,
    /// Angular damping coefficient (0 = no damping).
    #[serde(default)]
    pub angular_damping: f64,
}

impl Default for SceneMaterial {
    fn default() -> Self {
        Self {
            density: 1000.0,
            restitution: 0.3,
            friction: 0.5,
            dynamic_friction: None,
            linear_damping: 0.0,
            angular_damping: 0.0,
        }
    }
}

impl SceneMaterial {
    /// Material for static surfaces (high friction, no bounce).
    pub fn static_surface() -> Self {
        Self {
            density: 0.0,
            restitution: 0.0,
            friction: 0.8,
            ..Default::default()
        }
    }

    /// Rubber-like material (high restitution, high friction).
    pub fn rubber() -> Self {
        Self {
            density: 1200.0,
            restitution: 0.7,
            friction: 0.9,
            ..Default::default()
        }
    }

    /// Metal-like material (low restitution, moderate friction).
    pub fn metal() -> Self {
        Self {
            density: 7800.0,
            restitution: 0.1,
            friction: 0.4,
            ..Default::default()
        }
    }

    /// Ice-like material (very low friction).
    pub fn ice() -> Self {
        Self {
            density: 917.0,
            restitution: 0.1,
            friction: 0.05,
            ..Default::default()
        }
    }

    /// Effective dynamic friction (falls back to `friction`).
    pub fn dynamic_friction(&self) -> f64 {
        self.dynamic_friction.unwrap_or(self.friction)
    }
}

// ============================================================================
// SceneBody
// ============================================================================

/// A single rigid body in the scene.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneBody {
    /// Unique identifier within the scene.
    pub id: String,
    /// Collision shape.
    pub shape: SceneShape,
    /// World-space position [x, y, z] in metres.
    pub position: [f64; 3],
    /// World-space orientation as unit quaternion [qx, qy, qz, qw].
    #[serde(default = "default_quat")]
    pub rotation: [f64; 4],
    /// Linear velocity in m/s.
    #[serde(default)]
    pub velocity: [f64; 3],
    /// Angular velocity in rad/s.
    #[serde(default)]
    pub angular_velocity: [f64; 3],
    /// If `true` the body cannot move (zero inverse mass).
    #[serde(default)]
    pub is_static: bool,
    /// Material properties.
    #[serde(default)]
    pub material: SceneMaterial,
    /// Optional application tags for grouping or filtering.
    #[serde(default)]
    pub tags: Vec<String>,
}

fn default_quat() -> [f64; 4] {
    [0.0, 0.0, 0.0, 1.0]
}

impl SceneBody {
    /// Approximate mass in kg (shape volume × density).
    pub fn approx_mass(&self) -> f64 {
        if self.is_static {
            return f64::INFINITY;
        }
        self.shape.approx_volume() * self.material.density
    }

    /// Returns `true` if this body has the given tag.
    pub fn has_tag(&self, tag: &str) -> bool {
        self.tags.iter().any(|t| t == tag)
    }
}

// ============================================================================
// SceneConstraintKind
// ============================================================================

/// Type of constraint connecting two bodies.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SceneConstraintKind {
    /// Locks all relative degrees of freedom.
    Fixed,
    /// Hinge joint rotating around a single axis.
    Hinge {
        /// Hinge axis in world space (unit vector).
        axis: [f64; 3],
        /// Optional angular limits [lower, upper] in radians.
        #[serde(default)]
        limits: Option<[f64; 2]>,
    },
    /// Ball-and-socket joint.
    Ball {
        /// Optional cone half-angle limit in radians.
        #[serde(default)]
        limit_angle: Option<f64>,
    },
    /// Prismatic (slider) joint along an axis.
    Slider {
        /// Slide axis in world space (unit vector).
        axis: [f64; 3],
        /// Optional linear limits [lower, upper] in metres.
        #[serde(default)]
        limits: Option<[f64; 2]>,
    },
    /// Distance constraint keeping bodies a fixed distance apart.
    Distance {
        /// Target distance in metres.
        target: f64,
        /// Stiffness (spring constant).  `None` = rigid.
        #[serde(default)]
        stiffness: Option<f64>,
    },
}

impl SceneConstraintKind {
    /// Human-readable name.
    pub fn kind_name(&self) -> &'static str {
        match self {
            SceneConstraintKind::Fixed => "fixed",
            SceneConstraintKind::Hinge { .. } => "hinge",
            SceneConstraintKind::Ball { .. } => "ball",
            SceneConstraintKind::Slider { .. } => "slider",
            SceneConstraintKind::Distance { .. } => "distance",
        }
    }
}

// ============================================================================
// SceneConstraint
// ============================================================================

/// A constraint (joint) connecting two bodies.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneConstraint {
    /// Unique identifier.
    pub id: String,
    /// Joint type.
    pub kind: SceneConstraintKind,
    /// ID of the first connected body.
    pub body_a: String,
    /// ID of the second connected body.
    pub body_b: String,
    /// Anchor point in body A's local frame.
    #[serde(default)]
    pub anchor_a: [f64; 3],
    /// Anchor point in body B's local frame.
    #[serde(default)]
    pub anchor_b: [f64; 3],
    /// Maximum force before the constraint breaks (None = unbreakable).
    #[serde(default)]
    pub break_force: Option<f64>,
}

// ============================================================================
// SceneDescription
// ============================================================================

/// A complete, self-contained description of a physics scene.
///
/// Serializes to/from JSON via serde.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SceneDescription {
    /// Format version (current: 1).
    pub version: u32,
    /// Human-readable scene name.
    pub name: String,
    /// Gravitational acceleration [gx, gy, gz] in m/s².
    pub gravity: [f64; 3],
    /// Default simulation time step in seconds.
    pub dt: f64,
    /// Number of solver sub-steps per time step.
    #[serde(default = "default_substeps")]
    pub substeps: u32,
    /// All rigid bodies in the scene.
    #[serde(default)]
    pub bodies: Vec<SceneBody>,
    /// All constraints in the scene.
    #[serde(default)]
    pub constraints: Vec<SceneConstraint>,
}

fn default_substeps() -> u32 {
    4
}

impl SceneDescription {
    /// Create an empty scene with default settings.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            version: 1,
            name: name.into(),
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            substeps: 4,
            bodies: Vec::new(),
            constraints: Vec::new(),
        }
    }

    /// Deserialize from a JSON string.
    pub fn from_json(s: &str) -> Result<Self, String> {
        serde_json::from_str(s).map_err(|e| e.to_string())
    }

    /// Serialize to pretty-printed JSON.
    pub fn to_json(&self) -> Result<String, String> {
        serde_json::to_string_pretty(self).map_err(|e| e.to_string())
    }

    /// Serialize to compact (single-line) JSON.
    pub fn to_json_compact(&self) -> Result<String, String> {
        serde_json::to_string(self).map_err(|e| e.to_string())
    }

    // ------------------------------------------------------------------
    // Body queries
    // ------------------------------------------------------------------

    /// Total number of bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Number of static bodies.
    pub fn static_body_count(&self) -> usize {
        self.bodies.iter().filter(|b| b.is_static).count()
    }

    /// Number of dynamic bodies.
    pub fn dynamic_body_count(&self) -> usize {
        self.bodies.iter().filter(|b| !b.is_static).count()
    }

    /// Look up a body by its ID.
    pub fn body_by_id(&self, id: &str) -> Option<&SceneBody> {
        self.bodies.iter().find(|b| b.id == id)
    }

    /// Returns `true` if a body with the given ID exists.
    pub fn has_body(&self, id: &str) -> bool {
        self.body_by_id(id).is_some()
    }

    /// Collect all body IDs.
    pub fn body_ids(&self) -> Vec<&str> {
        self.bodies.iter().map(|b| b.id.as_str()).collect()
    }

    /// All bodies with the given tag.
    pub fn bodies_with_tag(&self, tag: &str) -> Vec<&SceneBody> {
        self.bodies.iter().filter(|b| b.has_tag(tag)).collect()
    }

    // ------------------------------------------------------------------
    // Constraint queries
    // ------------------------------------------------------------------

    /// Number of constraints.
    pub fn constraint_count(&self) -> usize {
        self.constraints.len()
    }

    /// Look up a constraint by ID.
    pub fn constraint_by_id(&self, id: &str) -> Option<&SceneConstraint> {
        self.constraints.iter().find(|c| c.id == id)
    }

    /// All constraints attached to a given body ID.
    pub fn constraints_for_body(&self, body_id: &str) -> Vec<&SceneConstraint> {
        self.constraints
            .iter()
            .filter(|c| c.body_a == body_id || c.body_b == body_id)
            .collect()
    }

    // ------------------------------------------------------------------
    // Summary
    // ------------------------------------------------------------------

    /// One-line human-readable summary.
    pub fn summary(&self) -> String {
        format!(
            "SceneDescription {{ name={:?}, bodies={} (static={}, dynamic={}), constraints={}, dt={:.4}s }}",
            self.name,
            self.body_count(),
            self.static_body_count(),
            self.dynamic_body_count(),
            self.constraint_count(),
            self.dt,
        )
    }
}

impl std::fmt::Display for SceneDescription {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.summary())
    }
}

// ============================================================================
// SceneBuilder
// ============================================================================

/// Fluent builder for constructing a [`SceneDescription`].
///
/// # Example
///
/// ```rust
/// use oxiphysics::scene::SceneBuilder;
///
/// let scene = SceneBuilder::new("demo")
///     .with_gravity([0.0, -9.81, 0.0])
///     .add_sphere("ball", [0.0, 5.0, 0.0], 0.5, false)
///     .add_static_box("floor", [0.0, 0.0, 0.0], [10.0, 0.1, 10.0])
///     .build();
///
/// assert_eq!(scene.body_count(), 2);
/// ```
pub struct SceneBuilder {
    description: SceneDescription,
}

impl SceneBuilder {
    /// Start a new builder with the given scene name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            description: SceneDescription::new(name),
        }
    }

    /// Set gravitational acceleration.
    pub fn with_gravity(mut self, gravity: [f64; 3]) -> Self {
        self.description.gravity = gravity;
        self
    }

    /// Set the default time step.
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.description.dt = dt;
        self
    }

    /// Set the number of solver sub-steps.
    pub fn with_substeps(mut self, substeps: u32) -> Self {
        self.description.substeps = substeps;
        self
    }

    /// Add a fully specified body.
    pub fn add_body(mut self, body: SceneBody) -> Self {
        self.description.bodies.push(body);
        self
    }

    /// Add a sphere body.
    pub fn add_sphere(
        mut self,
        id: impl Into<String>,
        position: [f64; 3],
        radius: f64,
        is_static: bool,
    ) -> Self {
        self.description.bodies.push(SceneBody {
            id: id.into(),
            shape: SceneShape::Sphere { radius },
            position,
            rotation: default_quat(),
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static,
            material: if is_static {
                SceneMaterial::static_surface()
            } else {
                SceneMaterial::default()
            },
            tags: Vec::new(),
        });
        self
    }

    /// Add a static (immovable) box.
    pub fn add_static_box(
        mut self,
        id: impl Into<String>,
        position: [f64; 3],
        half_extents: [f64; 3],
    ) -> Self {
        self.description.bodies.push(SceneBody {
            id: id.into(),
            shape: SceneShape::Box { half_extents },
            position,
            rotation: default_quat(),
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static: true,
            material: SceneMaterial::static_surface(),
            tags: Vec::new(),
        });
        self
    }

    /// Add a dynamic box.
    pub fn add_box(
        mut self,
        id: impl Into<String>,
        position: [f64; 3],
        half_extents: [f64; 3],
        is_static: bool,
    ) -> Self {
        self.description.bodies.push(SceneBody {
            id: id.into(),
            shape: SceneShape::Box { half_extents },
            position,
            rotation: default_quat(),
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static,
            material: if is_static {
                SceneMaterial::static_surface()
            } else {
                SceneMaterial::default()
            },
            tags: Vec::new(),
        });
        self
    }

    /// Add a capsule body.
    pub fn add_capsule(
        mut self,
        id: impl Into<String>,
        position: [f64; 3],
        radius: f64,
        half_height: f64,
        is_static: bool,
    ) -> Self {
        self.description.bodies.push(SceneBody {
            id: id.into(),
            shape: SceneShape::Capsule {
                radius,
                half_height,
            },
            position,
            rotation: default_quat(),
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static,
            material: if is_static {
                SceneMaterial::static_surface()
            } else {
                SceneMaterial::default()
            },
            tags: Vec::new(),
        });
        self
    }

    /// Add a ground plane (static, facing +Y).
    pub fn add_ground_plane(mut self, id: impl Into<String>) -> Self {
        self.description.bodies.push(SceneBody {
            id: id.into(),
            shape: SceneShape::Plane {
                normal: [0.0, 1.0, 0.0],
            },
            position: [0.0; 3],
            rotation: default_quat(),
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_static: true,
            material: SceneMaterial::static_surface(),
            tags: vec!["ground".to_string()],
        });
        self
    }

    /// Tag the last added body with an additional label.
    pub fn tag_last(mut self, tag: impl Into<String>) -> Self {
        if let Some(b) = self.description.bodies.last_mut() {
            b.tags.push(tag.into());
        }
        self
    }

    /// Set a custom material on the last added body.
    pub fn material_last(mut self, mat: SceneMaterial) -> Self {
        if let Some(b) = self.description.bodies.last_mut() {
            b.material = mat;
        }
        self
    }

    /// Set initial velocity on the last added body.
    pub fn velocity_last(mut self, v: [f64; 3]) -> Self {
        if let Some(b) = self.description.bodies.last_mut() {
            b.velocity = v;
        }
        self
    }

    /// Add a constraint.
    pub fn add_constraint(mut self, constraint: SceneConstraint) -> Self {
        self.description.constraints.push(constraint);
        self
    }

    /// Add a hinge joint between two bodies.
    pub fn add_hinge(
        mut self,
        id: impl Into<String>,
        body_a: impl Into<String>,
        body_b: impl Into<String>,
        axis: [f64; 3],
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        self.description.constraints.push(SceneConstraint {
            id: id.into(),
            kind: SceneConstraintKind::Hinge { axis, limits: None },
            body_a: body_a.into(),
            body_b: body_b.into(),
            anchor_a,
            anchor_b,
            break_force: None,
        });
        self
    }

    /// Add a ball joint between two bodies.
    pub fn add_ball_joint(
        mut self,
        id: impl Into<String>,
        body_a: impl Into<String>,
        body_b: impl Into<String>,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        self.description.constraints.push(SceneConstraint {
            id: id.into(),
            kind: SceneConstraintKind::Ball { limit_angle: None },
            body_a: body_a.into(),
            body_b: body_b.into(),
            anchor_a,
            anchor_b,
            break_force: None,
        });
        self
    }

    /// Finalise and return the [`SceneDescription`].
    pub fn build(self) -> SceneDescription {
        self.description
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builder_round_trip_json() {
        let scene = SceneBuilder::new("test_scene")
            .with_gravity([0.0, -9.81, 0.0])
            .with_dt(1.0 / 60.0)
            .add_sphere("ball", [0.0, 5.0, 0.0], 0.5, false)
            .add_static_box("floor", [0.0, -0.1, 0.0], [10.0, 0.1, 10.0])
            .build();

        assert_eq!(scene.body_count(), 2);
        assert_eq!(scene.static_body_count(), 1);
        assert_eq!(scene.dynamic_body_count(), 1);

        let json = scene.to_json().unwrap();
        let restored = SceneDescription::from_json(&json).unwrap();
        assert_eq!(restored.name, "test_scene");
        assert_eq!(restored.body_count(), 2);
        assert!(restored.has_body("ball"));
        assert!(restored.has_body("floor"));
    }

    #[test]
    fn body_material_presets() {
        let mat = SceneMaterial::rubber();
        assert!(mat.restitution > 0.5);
        let mat2 = SceneMaterial::ice();
        assert!(mat2.friction < 0.1);
    }

    #[test]
    fn constraint_round_trip() {
        let scene = SceneBuilder::new("hinge_test")
            .add_sphere("a", [0.0, 0.0, 0.0], 0.5, false)
            .add_sphere("b", [1.0, 0.0, 0.0], 0.5, false)
            .add_hinge(
                "j1",
                "a",
                "b",
                [0.0, 1.0, 0.0],
                [0.5, 0.0, 0.0],
                [-0.5, 0.0, 0.0],
            )
            .build();

        assert_eq!(scene.constraint_count(), 1);
        let json = scene.to_json().unwrap();
        let restored = SceneDescription::from_json(&json).unwrap();
        assert_eq!(restored.constraint_count(), 1);
        let c = restored.constraint_by_id("j1").unwrap();
        assert_eq!(c.body_a, "a");
        assert_eq!(c.body_b, "b");
    }

    #[test]
    fn tag_filtering() {
        let scene = SceneBuilder::new("tags")
            .add_sphere("s1", [0.0, 0.0, 0.0], 0.5, false)
            .tag_last("projectile")
            .add_sphere("s2", [1.0, 0.0, 0.0], 0.5, false)
            .tag_last("projectile")
            .add_static_box("wall", [5.0, 0.0, 0.0], [0.5, 5.0, 5.0])
            .build();

        assert_eq!(scene.bodies_with_tag("projectile").len(), 2);
        assert_eq!(scene.bodies_with_tag("wall").len(), 0);
    }

    #[test]
    fn shape_volume() {
        let sphere = SceneShape::Sphere { radius: 1.0 };
        let v = sphere.approx_volume();
        assert!((v - 4.0 * std::f64::consts::PI / 3.0).abs() < 1e-6);
    }

    #[test]
    fn scene_display() {
        let scene = SceneDescription::new("empty");
        let s = scene.to_string();
        assert!(s.contains("empty"));
    }

    #[test]
    fn bodies_in_radius_via_constraints_for_body() {
        let scene = SceneBuilder::new("chain")
            .add_sphere("a", [0.0, 0.0, 0.0], 0.3, false)
            .add_sphere("b", [1.0, 0.0, 0.0], 0.3, false)
            .add_sphere("c", [2.0, 0.0, 0.0], 0.3, false)
            .add_ball_joint("j_ab", "a", "b", [0.3, 0.0, 0.0], [-0.3, 0.0, 0.0])
            .add_ball_joint("j_bc", "b", "c", [0.3, 0.0, 0.0], [-0.3, 0.0, 0.0])
            .build();

        assert_eq!(scene.constraints_for_body("b").len(), 2);
        assert_eq!(scene.constraints_for_body("a").len(), 1);
    }
}
