// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Canvas renderer hints for browser-based physics visualizations.
//!
//! This module provides lightweight data structures used to communicate
//! per-frame rendering information from the physics engine to a JavaScript
//! canvas 2D or WebGL renderer. All types are plain-data and serialize cleanly
//! to JSON for `postMessage` transfer or structured clone.
//!
//! ## Typical per-frame flow
//!
//! 1. Call `RendererHints::from_engine(engine)` to collect all hints.
//! 2. Transfer the `RendererHints` JSON blob to the main thread (or worker).
//! 3. The renderer iterates over `bounding_boxes`, `contact_points`, and
//!    `velocity_vectors` to draw debug overlays.

#![allow(missing_docs)]
#![allow(dead_code)]

use crate::engine::WasmPhysicsEngine;
use crate::types::ContactResult;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// AxisAlignedBoundingBox
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box (AABB) for a single body.
///
/// Used by the canvas renderer to draw wireframe extents.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Aabb {
    /// Body handle.
    pub handle: u32,
    /// Minimum corner `[min_x, min_y, min_z]`.
    pub min: [f64; 3],
    /// Maximum corner `[max_x, max_y, max_z]`.
    pub max: [f64; 3],
}

impl Aabb {
    /// Create an `Aabb` from a center point and half-extents.
    pub fn from_center_half_extents(handle: u32, center: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            handle,
            min: [
                center[0] - half[0],
                center[1] - half[1],
                center[2] - half[2],
            ],
            max: [
                center[0] + half[0],
                center[1] + half[1],
                center[2] + half[2],
            ],
        }
    }

    /// Create an `Aabb` from a center point and a uniform radius (sphere AABB).
    pub fn from_sphere(handle: u32, center: [f64; 3], radius: f64) -> Self {
        Self::from_center_half_extents(handle, center, [radius; 3])
    }

    /// Compute the center of this AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Compute the half-extents of this AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }

    /// Returns `true` if the given point `[x, y, z]` is inside (or on) this AABB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Returns `true` if this AABB overlaps with `other`.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Expand this AABB by a margin on all sides.
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            handle: self.handle,
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }
}

// ---------------------------------------------------------------------------
// ContactPointHint
// ---------------------------------------------------------------------------

/// Rendering hint for a single contact point.
///
/// Encodes the world-space position and the contact normal for drawing
/// contact gizmos (arrows, crosses, etc.) on the canvas.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ContactPointHint {
    /// World-space position of the contact.
    pub position: [f64; 3],
    /// Contact normal (unit vector from body B toward body A).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Magnitude of the impulse applied (for scaling the arrow).
    pub impulse: f64,
    /// Handle of body A.
    pub body_a: u32,
    /// Handle of body B.
    pub body_b: u32,
    /// Whether this is a newly detected contact (vs. persisting).
    pub is_new: bool,
}

impl ContactPointHint {
    /// Build a `ContactPointHint` from a [`ContactResult`].
    pub fn from_contact(c: &ContactResult) -> Self {
        // Use midpoint of the two contact points
        let px = (c.point_on_a[0] + c.point_on_b[0]) * 0.5;
        let py = (c.point_on_a[1] + c.point_on_b[1]) * 0.5;
        let pz = (c.point_on_a[2] + c.point_on_b[2]) * 0.5;
        Self {
            position: [px, py, pz],
            normal: c.normal,
            depth: c.depth,
            impulse: c.impulse,
            body_a: c.body_a,
            body_b: c.body_b,
            is_new: c.is_new,
        }
    }
}

// ---------------------------------------------------------------------------
// VelocityVector
// ---------------------------------------------------------------------------

/// Rendering hint for a body velocity vector.
///
/// Used to draw velocity arrows on a debug canvas overlay.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VelocityVector {
    /// Body handle.
    pub handle: u32,
    /// World-space origin of the arrow (body center of mass).
    pub origin: [f64; 3],
    /// Velocity vector `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Speed (magnitude of velocity).
    pub speed: f64,
    /// Whether the body is sleeping (affects rendering color).
    pub is_sleeping: bool,
}

impl VelocityVector {
    /// Create from components.
    pub fn new(handle: u32, origin: [f64; 3], velocity: [f64; 3], is_sleeping: bool) -> Self {
        let speed =
            (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2])
                .sqrt();
        Self {
            handle,
            origin,
            velocity,
            speed,
            is_sleeping,
        }
    }
}

// ---------------------------------------------------------------------------
// BodyRenderHint
// ---------------------------------------------------------------------------

/// Per-body rendering hint: AABB, velocity arrow, sleeping status.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyRenderHint {
    /// Body handle.
    pub handle: u32,
    /// Axis-aligned bounding box (computed from sphere/box collider).
    pub aabb: Option<Aabb>,
    /// Velocity arrow hint.
    pub velocity_vector: VelocityVector,
    /// Whether the body is currently sleeping.
    pub is_sleeping: bool,
    /// Whether this body is static.
    pub is_static: bool,
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
}

// ---------------------------------------------------------------------------
// RendererHints (aggregate per-frame data)
// ---------------------------------------------------------------------------

/// Aggregate per-frame renderer hints collected from a `WasmPhysicsEngine`.
///
/// Serializable to JSON for transfer to the main thread renderer.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::WasmPhysicsEngine;
/// use oxiphysics_wasm::renderer::RendererHints;
///
/// let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
/// let b = engine.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
/// engine.add_sphere_collider(b, 0.5);
/// engine.step(1.0 / 60.0);
///
/// let hints = RendererHints::from_engine(&engine, true, true, true);
/// assert!(!hints.body_hints.is_empty());
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RendererHints {
    /// Per-body rendering information.
    pub body_hints: Vec<BodyRenderHint>,
    /// Contact point hints for this frame.
    pub contact_points: Vec<ContactPointHint>,
    /// Simulation time at which these hints were collected.
    pub sim_time: f64,
    /// Number of active bodies.
    pub active_body_count: u32,
    /// Number of contacts.
    pub contact_count: u32,
    /// Whether AABBs were requested.
    pub include_aabbs: bool,
    /// Whether velocity vectors were requested.
    pub include_velocities: bool,
    /// Whether contact points were requested.
    pub include_contacts: bool,
}

impl RendererHints {
    /// Collect renderer hints from a `WasmPhysicsEngine`.
    ///
    /// - `include_aabbs`: compute AABB from first sphere collider.
    /// - `include_velocities`: build velocity vectors.
    /// - `include_contacts`: copy contact points.
    pub fn from_engine(
        engine: &WasmPhysicsEngine,
        include_aabbs: bool,
        include_velocities: bool,
        include_contacts: bool,
    ) -> Self {
        let handles = engine.get_all_body_handles();
        let mut body_hints = Vec::with_capacity(handles.len());

        for &h in &handles {
            let pos = engine.get_position(h);
            let rot = engine.get_rotation(h);
            let vel = engine.get_velocity(h);
            let sleeping = engine
                .get_body_state(h)
                .map(|s| s.is_sleeping)
                .unwrap_or(false);
            let is_static = engine.body_inv_mass(h) == 0.0 && !engine.body_is_dynamic(h);

            // Build AABB from sphere radius (default radius 0.5 if no collider info)
            let aabb = if include_aabbs {
                Some(Aabb::from_sphere(h, pos, 0.5))
            } else {
                None
            };

            let velocity_vector = VelocityVector::new(h, pos, vel, sleeping);

            body_hints.push(BodyRenderHint {
                handle: h,
                aabb,
                velocity_vector,
                is_sleeping: sleeping,
                is_static,
                position: pos,
                rotation: rot,
            });
        }

        let contact_points = if include_contacts {
            engine
                .get_contacts()
                .iter()
                .map(ContactPointHint::from_contact)
                .collect()
        } else {
            Vec::new()
        };

        Self {
            active_body_count: engine.get_body_count(),
            contact_count: engine.get_contact_count(),
            sim_time: engine.time(),
            body_hints,
            contact_points,
            include_aabbs,
            include_velocities,
            include_contacts,
        }
    }

    /// Serialize to a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }

    /// Deserialize from a JSON string.
    pub fn from_json(json: &str) -> Option<Self> {
        serde_json::from_str(json).ok()
    }

    /// Return all AABB data as a flat `Vec`f64`: `\[handle, min_x, min_y, min_z, max_x, max_y, max_z, ...\]`.
    pub fn aabbs_flat(&self) -> Vec<f64> {
        let mut out = Vec::new();
        for hint in &self.body_hints {
            if let Some(aabb) = &hint.aabb {
                out.push(aabb.handle as f64);
                out.extend_from_slice(&aabb.min);
                out.extend_from_slice(&aabb.max);
            }
        }
        out
    }

    /// Return all contact positions as a flat `Vec`f64`: `[x0, y0, z0, x1, ...]`.
    pub fn contact_positions_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.contact_points.len() * 3);
        for cp in &self.contact_points {
            out.extend_from_slice(&cp.position);
        }
        out
    }

    /// Return all velocity vectors as a flat `Vec`f64`: `\[ox, oy, oz, vx, vy, vz, ...\]`.
    pub fn velocity_vectors_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.body_hints.len() * 6);
        for hint in &self.body_hints {
            out.extend_from_slice(&hint.velocity_vector.origin);
            out.extend_from_slice(&hint.velocity_vector.velocity);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// DebugDrawCommand
// ---------------------------------------------------------------------------

/// A single draw command for a canvas 2D renderer.
///
/// Designed to be consumed by a JavaScript renderer that iterates
/// over a list of commands and executes them in order.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum DebugDrawCommand {
    /// Draw a circle at `center` with `radius` in a given `color` (CSS hex).
    Circle {
        /// Center x (world units).
        cx: f64,
        /// Center y (world units).
        cy: f64,
        /// Radius (world units).
        radius: f64,
        /// CSS color string (e.g. `"#ff0000"`).
        color: String,
        /// Line width.
        line_width: f64,
    },
    /// Draw a line between two points.
    Line {
        /// Start x.
        x0: f64,
        /// Start y.
        y0: f64,
        /// End x.
        x1: f64,
        /// End y.
        y1: f64,
        /// CSS color string.
        color: String,
        /// Line width.
        line_width: f64,
    },
    /// Draw an axis-aligned rectangle.
    Rect {
        /// Left edge.
        x: f64,
        /// Top edge.
        y: f64,
        /// Width.
        width: f64,
        /// Height.
        height: f64,
        /// CSS color string.
        color: String,
        /// Line width (0 = filled).
        line_width: f64,
    },
    /// Draw a text label.
    Text {
        /// X position.
        x: f64,
        /// Y position.
        y: f64,
        /// Label text.
        text: String,
        /// CSS color string.
        color: String,
        /// Font size in pixels.
        font_size: f64,
    },
}

/// Build a list of `DebugDrawCommand` from `RendererHints` for a 2D top-down view (XZ plane).
///
/// Bodies are drawn as circles, contact points as red crosses, velocities as green arrows.
pub fn build_draw_commands(hints: &RendererHints, scale: f64) -> Vec<DebugDrawCommand> {
    let mut cmds = Vec::new();

    for hint in &hints.body_hints {
        let cx = hint.position[0] * scale;
        let cy = hint.position[2] * scale; // top-down: Z maps to canvas Y
        let radius = hint
            .aabb
            .as_ref()
            .map(|a| a.half_extents()[0] * scale)
            .unwrap_or(10.0);

        let color = if hint.is_sleeping {
            "#808080".to_string()
        } else if hint.is_static {
            "#4444ff".to_string()
        } else {
            "#00cc00".to_string()
        };

        cmds.push(DebugDrawCommand::Circle {
            cx,
            cy,
            radius,
            color,
            line_width: 1.5,
        });

        // Velocity arrow (XZ plane)
        let vx = hint.velocity_vector.velocity[0];
        let vz = hint.velocity_vector.velocity[2];
        if vx * vx + vz * vz > 1e-6 {
            cmds.push(DebugDrawCommand::Line {
                x0: cx,
                y0: cy,
                x1: cx + vx * scale * 0.1,
                y1: cy + vz * scale * 0.1,
                color: "#00ffff".to_string(),
                line_width: 1.0,
            });
        }

        // Handle label
        cmds.push(DebugDrawCommand::Text {
            x: cx + 4.0,
            y: cy - 4.0,
            text: format!("{}", hint.handle),
            color: "#ffffff".to_string(),
            font_size: 10.0,
        });
    }

    // Contact points as small crosses
    for cp in &hints.contact_points {
        let cx = cp.position[0] * scale;
        let cy = cp.position[2] * scale;
        let s = 4.0;
        cmds.push(DebugDrawCommand::Line {
            x0: cx - s,
            y0: cy,
            x1: cx + s,
            y1: cy,
            color: "#ff0000".to_string(),
            line_width: 1.5,
        });
        cmds.push(DebugDrawCommand::Line {
            x0: cx,
            y0: cy - s,
            x1: cx,
            y1: cy + s,
            color: "#ff0000".to_string(),
            line_width: 1.5,
        });
    }

    cmds
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WasmPhysicsEngine;
    use crate::renderer::Aabb;

    use crate::renderer::RendererHints;

    // --- Aabb ---

    #[test]
    fn test_aabb_from_sphere() {
        let aabb = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(aabb.min, [-1.0, -1.0, -1.0]);
        assert_eq!(aabb.max, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_aabb_center() {
        let aabb = Aabb::from_sphere(0, [2.0, 3.0, 4.0], 1.0);
        let c = aabb.center();
        assert!((c[0] - 2.0).abs() < 1e-10);
        assert!((c[1] - 3.0).abs() < 1e-10);
        assert!((c[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_contains_point() {
        let aabb = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 2.0);
        assert!(aabb.contains_point([1.0, 0.0, 0.0]));
        assert!(!aabb.contains_point([3.0, 0.0, 0.0]));
    }

    #[test]
    fn test_aabb_overlaps() {
        let a = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        let b = Aabb::from_sphere(1, [1.5, 0.0, 0.0], 1.0);
        assert!(a.overlaps(&b), "should overlap at x=0.5..1.5");
        let c = Aabb::from_sphere(2, [5.0, 0.0, 0.0], 1.0);
        assert!(!a.overlaps(&c), "should not overlap");
    }

    #[test]
    fn test_aabb_expanded() {
        let aabb = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 1.0);
        let exp = aabb.expanded(0.5);
        assert!((exp.min[0] + 1.5).abs() < 1e-10);
        assert!((exp.max[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_half_extents() {
        let aabb = Aabb::from_center_half_extents(0, [0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let he = aabb.half_extents();
        assert!((he[0] - 2.0).abs() < 1e-10);
        assert!((he[1] - 3.0).abs() < 1e-10);
        assert!((he[2] - 4.0).abs() < 1e-10);
    }

    // --- ContactPointHint ---

    #[test]
    fn test_contact_point_hint_from_contact() {
        let mut c = ContactResult::new(0, 1);
        c.point_on_a = [0.0, 0.0, 0.0];
        c.point_on_b = [2.0, 0.0, 0.0];
        c.impulse = 5.0;
        let hint = ContactPointHint::from_contact(&c);
        assert!((hint.position[0] - 1.0).abs() < 1e-10);
        assert!((hint.impulse - 5.0).abs() < 1e-10);
    }

    // --- VelocityVector ---

    #[test]
    fn test_velocity_vector_speed() {
        let v = VelocityVector::new(0, [0.0, 0.0, 0.0], [3.0, 4.0, 0.0], false);
        assert!((v.speed - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_vector_sleeping() {
        let v = VelocityVector::new(1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], true);
        assert!(v.is_sleeping);
        assert!((v.speed).abs() < 1e-10);
    }

    // --- RendererHints ---

    #[test]
    fn test_renderer_hints_from_engine_empty() {
        let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let hints = RendererHints::from_engine(&engine, true, true, true);
        assert_eq!(hints.active_body_count, 0);
        assert!(hints.body_hints.is_empty());
        assert!(hints.contact_points.is_empty());
    }

    #[test]
    fn test_renderer_hints_from_engine_with_bodies() {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        engine.add_static_body(0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, true, false);
        assert_eq!(hints.active_body_count, 2);
        assert_eq!(hints.body_hints.len(), 2);
    }

    #[test]
    fn test_renderer_hints_json_roundtrip() {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
        let hints = RendererHints::from_engine(&engine, true, true, true);
        let json = hints.to_json();
        let back = RendererHints::from_json(&json).expect("should deserialize");
        assert_eq!(back.active_body_count, hints.active_body_count);
    }

    #[test]
    fn test_renderer_hints_aabbs_flat() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let flat = hints.aabbs_flat();
        // 1 body × 7 floats (handle + min[3] + max[3])
        assert_eq!(flat.len(), 7);
    }

    #[test]
    fn test_renderer_hints_contact_positions_flat() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b1 = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0);
        engine.add_sphere_collider(b0, 1.0);
        engine.add_sphere_collider(b1, 1.0);
        engine.step(1.0 / 60.0);
        let hints = RendererHints::from_engine(&engine, false, false, true);
        // flat = 3 floats per contact
        assert_eq!(
            hints.contact_positions_flat().len(),
            hints.contact_points.len() * 3
        );
    }

    #[test]
    fn test_renderer_hints_velocity_vectors_flat() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, false, true, false);
        // 1 body × 6 floats (origin[3] + vel[3])
        assert_eq!(hints.velocity_vectors_flat().len(), 6);
    }

    // --- DebugDrawCommand ---

    #[test]
    fn test_build_draw_commands_no_bodies() {
        let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let hints = RendererHints::from_engine(&engine, true, true, true);
        let cmds = build_draw_commands(&hints, 10.0);
        assert!(cmds.is_empty());
    }

    #[test]
    fn test_build_draw_commands_with_body() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, true, false);
        let cmds = build_draw_commands(&hints, 10.0);
        // At minimum: 1 circle + 1 text label
        assert!(!cmds.is_empty());
        let has_circle = cmds
            .iter()
            .any(|c| matches!(c, DebugDrawCommand::Circle { .. }));
        assert!(has_circle);
    }

    #[test]
    fn test_build_draw_commands_contact_crosses() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b1 = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0);
        engine.add_sphere_collider(b0, 1.0);
        engine.add_sphere_collider(b1, 1.0);
        engine.step(1.0 / 60.0);
        let hints = RendererHints::from_engine(&engine, true, true, true);
        let cmds = build_draw_commands(&hints, 10.0);
        // Contact cross = 2 lines per contact
        let line_count = cmds
            .iter()
            .filter(|c| matches!(c, DebugDrawCommand::Line { .. }))
            .count();
        if hints.contact_count > 0 {
            assert!(line_count >= 2);
        }
    }

    #[test]
    fn test_renderer_hints_no_aabbs() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, false, false, false);
        // No AABBs requested → all body aabb fields should be None
        for h in &hints.body_hints {
            assert!(h.aabb.is_none());
        }
    }
}

// ===========================================================================
// SSAO, instanced mesh, and screen-space AABB extensions
// ===========================================================================

/// Screen Screen Ambient Occlusion (SSAO) configuration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SsaoConfig {
    /// Whether SSAO is enabled.
    pub enabled: bool,
    /// Occlusion intensity in `\[0.0, 1.0\]` (0 = no occlusion, 1 = full).
    pub intensity: f64,
    /// Sampling radius in world units.
    pub radius: f64,
    /// Number of samples per pixel (higher = better quality, lower = faster).
    pub sample_count: u32,
    /// Bilateral blur kernel size (odd integer).
    pub blur_kernel_size: u32,
}

impl SsaoConfig {
    /// Default SSAO configuration (disabled, neutral settings).
    pub fn default_disabled() -> Self {
        Self {
            enabled: false,
            intensity: 0.5,
            radius: 0.5,
            sample_count: 16,
            blur_kernel_size: 3,
        }
    }

    /// Preset for high-quality SSAO.
    pub fn high_quality() -> Self {
        Self {
            enabled: true,
            intensity: 0.8,
            radius: 0.3,
            sample_count: 64,
            blur_kernel_size: 5,
        }
    }

    /// Preset for fast (low-quality) SSAO.
    pub fn fast() -> Self {
        Self {
            enabled: true,
            intensity: 0.5,
            radius: 0.5,
            sample_count: 8,
            blur_kernel_size: 3,
        }
    }

    /// Clamp intensity to `\[0.0, 1.0\]`.
    pub fn with_intensity(mut self, intensity: f64) -> Self {
        self.intensity = intensity.clamp(0.0, 1.0);
        self
    }
}

/// A single instanced mesh entry for GPU instanced rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InstancedMeshEntry {
    /// Body handle this instance corresponds to.
    pub handle: u32,
    /// Position `\[x, y, z\]`.
    pub position: [f64; 3],
    /// Orientation quaternion `\[x, y, z, w\]`.
    pub rotation: [f64; 4],
    /// Non-uniform scale `\[sx, sy, sz\]`.
    pub scale: [f64; 3],
    /// RGBA color `\[r, g, b, a\]` in `\[0.0, 1.0\]`.
    pub color: [f64; 4],
}

impl InstancedMeshEntry {
    /// Create an entry from a body render hint with default white color and unit scale.
    pub fn from_body_hint(hint: &BodyRenderHint) -> Self {
        Self {
            handle: hint.handle,
            position: hint.position,
            rotation: hint.rotation,
            scale: [1.0, 1.0, 1.0],
            color: [1.0, 1.0, 1.0, 1.0],
        }
    }

    /// Return the transform as a flat `\[px, py, pz, rx, ry, rz, rw, sx, sy, sz\]`.
    pub fn to_flat_transform(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(10);
        v.extend_from_slice(&self.position);
        v.extend_from_slice(&self.rotation);
        v.extend_from_slice(&self.scale);
        v
    }
}

/// Screen-space (2D) bounding box of a projected 3D AABB.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScreenSpaceAabb {
    /// Body handle.
    pub handle: u32,
    /// Minimum screen-space pixel `\[x, y\]`.
    pub min_px: [f64; 2],
    /// Maximum screen-space pixel `\[x, y\]`.
    pub max_px: [f64; 2],
    /// Screen-space width in pixels.
    pub width_px: f64,
    /// Screen-space height in pixels.
    pub height_px: f64,
    /// Whether the body is partially or fully behind the camera.
    pub clipped: bool,
}

impl ScreenSpaceAabb {
    /// Compute a screen-space AABB for a body given its world-space AABB
    /// and a simple orthographic camera projection.
    ///
    /// `viewport_size` is `\[width_px, height_px\]`.
    /// `camera_scale` is pixels per world unit.
    /// `camera_offset` is the world-space origin mapped to the screen center.
    pub fn from_aabb_orthographic(
        aabb: &Aabb,
        viewport_size: [f64; 2],
        camera_scale: f64,
        camera_offset: [f64; 2],
    ) -> Self {
        // Project min/max corners using orthographic (XY plane projection)
        let cx = viewport_size[0] * 0.5;
        let cy = viewport_size[1] * 0.5;

        let project = |wx: f64, wy: f64| -> [f64; 2] {
            [
                cx + (wx - camera_offset[0]) * camera_scale,
                cy - (wy - camera_offset[1]) * camera_scale, // Y flipped
            ]
        };

        let corners = [
            project(aabb.min[0], aabb.min[1]),
            project(aabb.min[0], aabb.max[1]),
            project(aabb.max[0], aabb.min[1]),
            project(aabb.max[0], aabb.max[1]),
        ];

        let min_px_x = corners.iter().map(|c| c[0]).fold(f64::INFINITY, f64::min);
        let min_px_y = corners.iter().map(|c| c[1]).fold(f64::INFINITY, f64::min);
        let max_px_x = corners
            .iter()
            .map(|c| c[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let max_px_y = corners
            .iter()
            .map(|c| c[1])
            .fold(f64::NEG_INFINITY, f64::max);

        let clipped = max_px_x < 0.0
            || min_px_x > viewport_size[0]
            || max_px_y < 0.0
            || min_px_y > viewport_size[1];

        Self {
            handle: aabb.handle,
            min_px: [min_px_x, min_px_y],
            max_px: [max_px_x, max_px_y],
            width_px: (max_px_x - min_px_x).max(0.0),
            height_px: (max_px_y - min_px_y).max(0.0),
            clipped,
        }
    }
}

/// WASM renderer with SSAO, instanced mesh, and screen-space AABB support.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WasmRenderer {
    /// Current SSAO configuration.
    ssao: SsaoConfig,
    /// Current instanced mesh entries (one per active body).
    instances: Vec<InstancedMeshEntry>,
    /// Viewport width in pixels.
    viewport_width: f64,
    /// Viewport height in pixels.
    viewport_height: f64,
}

impl WasmRenderer {
    /// Create a new renderer with the given viewport dimensions.
    pub fn new(viewport_width: f64, viewport_height: f64) -> Self {
        Self {
            ssao: SsaoConfig::default_disabled(),
            instances: Vec::new(),
            viewport_width,
            viewport_height,
        }
    }

    // -----------------------------------------------------------------------
    // Ambient occlusion
    // -----------------------------------------------------------------------

    /// Set the SSAO configuration.
    pub fn set_ambient_occlusion(&mut self, config: SsaoConfig) {
        self.ssao = config;
    }

    /// Return the current SSAO configuration.
    pub fn ssao_config(&self) -> &SsaoConfig {
        &self.ssao
    }

    /// Enable or disable SSAO without changing other settings.
    pub fn toggle_ssao(&mut self, enabled: bool) {
        self.ssao.enabled = enabled;
    }

    // -----------------------------------------------------------------------
    // Instanced rendering
    // -----------------------------------------------------------------------

    /// Update the instanced mesh buffer from a set of `RendererHints`.
    ///
    /// One instance entry is created per body hint.
    pub fn update_instanced_mesh(&mut self, hints: &RendererHints) {
        self.instances = hints
            .body_hints
            .iter()
            .map(InstancedMeshEntry::from_body_hint)
            .collect();
    }

    /// Return the instanced mesh entries.
    pub fn get_instances(&self) -> &[InstancedMeshEntry] {
        &self.instances
    }

    /// Return the instance transforms as a flat buffer for GPU upload.
    ///
    /// Format per instance: `\[px, py, pz, rx, ry, rz, rw, sx, sy, sz\]` (10 floats).
    pub fn instance_buffer_flat(&self) -> Vec<f64> {
        self.instances
            .iter()
            .flat_map(|inst| inst.to_flat_transform())
            .collect()
    }

    // -----------------------------------------------------------------------
    // Screen-space bounding box
    // -----------------------------------------------------------------------

    /// Compute screen-space bounding boxes for all body AABBs in the hints,
    /// using an orthographic projection.
    ///
    /// `camera_scale` is pixels per world unit (zoom level).
    /// `camera_offset` is the world position mapped to the screen center.
    pub fn compute_screen_space_bounding_box(
        &self,
        hints: &RendererHints,
        camera_scale: f64,
        camera_offset: [f64; 2],
    ) -> Vec<ScreenSpaceAabb> {
        hints
            .body_hints
            .iter()
            .filter_map(|hint| hint.aabb.as_ref())
            .map(|aabb| {
                ScreenSpaceAabb::from_aabb_orthographic(
                    aabb,
                    [self.viewport_width, self.viewport_height],
                    camera_scale,
                    camera_offset,
                )
            })
            .collect()
    }
}

// ===========================================================================
// Tests for WasmRenderer (SSAO, instanced mesh, screen-space AABB)
// ===========================================================================

#[cfg(test)]
mod wasm_renderer_tests {

    use crate::WasmPhysicsEngine;
    use crate::renderer::Aabb;

    use crate::renderer::RendererHints;
    use crate::renderer::ScreenSpaceAabb;

    use crate::renderer::SsaoConfig;
    use crate::renderer::WasmRenderer;

    // --- SsaoConfig ---

    #[test]
    fn test_ssao_default_disabled() {
        let cfg = SsaoConfig::default_disabled();
        assert!(!cfg.enabled, "default SSAO should be disabled");
    }

    #[test]
    fn test_ssao_high_quality_preset() {
        let cfg = SsaoConfig::high_quality();
        assert!(cfg.enabled);
        assert_eq!(cfg.sample_count, 64);
        assert!(cfg.intensity > 0.5);
    }

    #[test]
    fn test_ssao_intensity_clamped() {
        let cfg = SsaoConfig::default_disabled().with_intensity(2.0);
        assert!(
            (cfg.intensity - 1.0).abs() < 1e-10,
            "intensity clamped to 1.0"
        );
        let cfg2 = SsaoConfig::default_disabled().with_intensity(-1.0);
        assert!(
            (cfg2.intensity - 0.0).abs() < 1e-10,
            "intensity clamped to 0.0"
        );
    }

    #[test]
    fn test_set_ambient_occlusion() {
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.set_ambient_occlusion(SsaoConfig::high_quality());
        assert!(renderer.ssao_config().enabled);
        assert_eq!(renderer.ssao_config().sample_count, 64);
    }

    #[test]
    fn test_toggle_ssao() {
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.toggle_ssao(true);
        assert!(renderer.ssao_config().enabled);
        renderer.toggle_ssao(false);
        assert!(!renderer.ssao_config().enabled);
    }

    // --- Instanced mesh ---

    #[test]
    fn test_update_instanced_mesh_empty() {
        let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        assert!(renderer.get_instances().is_empty());
    }

    #[test]
    fn test_update_instanced_mesh_one_body() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 3.0, 5.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        assert_eq!(renderer.get_instances().len(), 1);
        let inst = &renderer.get_instances()[0];
        assert!((inst.position[0] - 3.0).abs() < 1e-12);
        assert!((inst.position[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_instance_buffer_flat_length() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        engine.add_dynamic_body(2.0, 1.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let mut renderer = WasmRenderer::new(800.0, 600.0);
        renderer.update_instanced_mesh(&hints);
        let flat = renderer.instance_buffer_flat();
        // 2 bodies × 10 floats each
        assert_eq!(flat.len(), 20, "flat buffer should have 20 elements");
    }

    // --- Screen-space AABB ---

    #[test]
    fn test_screen_space_aabb_center_body() {
        // A body at origin with AABB radius 0.5, orthographic scale 100 px/unit,
        // viewport 800×600, camera offset [0,0] → body should appear near center.
        let aabb = Aabb::from_sphere(0, [0.0, 0.0, 0.0], 0.5);
        let ss = ScreenSpaceAabb::from_aabb_orthographic(&aabb, [800.0, 600.0], 100.0, [0.0, 0.0]);
        // Center should be near [400, 300]
        let cx = (ss.min_px[0] + ss.max_px[0]) * 0.5;
        let cy = (ss.min_px[1] + ss.max_px[1]) * 0.5;
        assert!((cx - 400.0).abs() < 1.0, "screen cx = {}", cx);
        assert!((cy - 300.0).abs() < 1.0, "screen cy = {}", cy);
        assert!(!ss.clipped);
    }

    #[test]
    fn test_screen_space_aabb_size() {
        // AABB half=1.0, scale=50 → screen width = 2*1.0*50 = 100 px
        let aabb = Aabb::from_sphere(1, [0.0, 0.0, 0.0], 1.0);
        let ss = ScreenSpaceAabb::from_aabb_orthographic(&aabb, [800.0, 600.0], 50.0, [0.0, 0.0]);
        assert!(
            (ss.width_px - 100.0).abs() < 1.0,
            "width_px = {}",
            ss.width_px
        );
        assert!(
            (ss.height_px - 100.0).abs() < 1.0,
            "height_px = {}",
            ss.height_px
        );
    }

    #[test]
    fn test_compute_screen_space_bounding_box_empty() {
        // No AABBs in hints → empty result
        let engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, false, false, false);
        let renderer = WasmRenderer::new(800.0, 600.0);
        let boxes = renderer.compute_screen_space_bounding_box(&hints, 100.0, [0.0, 0.0]);
        assert!(boxes.is_empty());
    }

    #[test]
    fn test_compute_screen_space_bounding_box_one_body() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let hints = RendererHints::from_engine(&engine, true, false, false);
        let renderer = WasmRenderer::new(800.0, 600.0);
        let boxes = renderer.compute_screen_space_bounding_box(&hints, 100.0, [0.0, 0.0]);
        assert_eq!(boxes.len(), 1);
        assert!(boxes[0].width_px > 0.0);
        assert!(boxes[0].height_px > 0.0);
    }
}

// ===========================================================================
// WebGPU render pipeline descriptor helpers
// ===========================================================================

/// Blend mode for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BlendMode {
    /// No blending — overwrite destination.
    Opaque,
    /// Standard alpha blending.
    AlphaBlend,
    /// Additive blending.
    Additive,
    /// Pre-multiplied alpha blending.
    PremultipliedAlpha,
}

/// Depth test comparison function.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DepthCompare {
    /// Always pass.
    Always,
    /// Pass if new depth < stored depth.
    Less,
    /// Pass if new depth <= stored depth.
    LessEqual,
    /// Pass if new depth > stored depth.
    Greater,
    /// Pass if new depth == stored depth.
    Equal,
    /// Never pass.
    Never,
}

/// A descriptor for a render pipeline (GPU render pass configuration).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RenderPipelineDesc {
    /// Human-readable label.
    pub label: String,
    /// Blend mode.
    pub blend_mode: BlendMode,
    /// Whether depth testing is enabled.
    pub depth_test: bool,
    /// Whether depth writing is enabled.
    pub depth_write: bool,
    /// Depth comparison function.
    pub depth_compare: DepthCompare,
    /// Whether back-face culling is enabled.
    pub backface_culling: bool,
    /// Number of color render targets.
    pub color_target_count: u32,
    /// MSAA sample count (1, 2, 4, or 8).
    pub sample_count: u32,
}

impl RenderPipelineDesc {
    /// Default opaque pipeline (depth test on, writes on, back-face culling on).
    pub fn opaque(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: true,
            depth_write: true,
            depth_compare: DepthCompare::Less,
            backface_culling: true,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Transparent pipeline (alpha blend, depth test on, depth write off).
    pub fn transparent(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::AlphaBlend,
            depth_test: true,
            depth_write: false,
            depth_compare: DepthCompare::Less,
            backface_culling: false,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Shadow-map pipeline (depth-only, no color targets).
    pub fn shadow_map(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: true,
            depth_write: true,
            depth_compare: DepthCompare::Less,
            backface_culling: true,
            color_target_count: 0,
            sample_count: 1,
        }
    }

    /// Post-processing pipeline (no depth test, no depth write, full-screen quad).
    pub fn post_process(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            blend_mode: BlendMode::Opaque,
            depth_test: false,
            depth_write: false,
            depth_compare: DepthCompare::Always,
            backface_culling: false,
            color_target_count: 1,
            sample_count: 1,
        }
    }

    /// Enable MSAA with the given sample count.
    pub fn with_msaa(mut self, samples: u32) -> Self {
        self.sample_count = samples.clamp(1, 8);
        self
    }

    /// Set the number of color targets (for MRT / deferred shading).
    pub fn with_color_targets(mut self, count: u32) -> Self {
        self.color_target_count = count;
        self
    }
}

impl Default for RenderPipelineDesc {
    fn default() -> Self {
        Self::opaque("default")
    }
}

// ===========================================================================
// Instanced rendering data
// ===========================================================================

/// Per-instance data for GPU instanced rendering (transform + color).
///
/// Layout: `\[px, py, pz, qx, qy, qz, qw, sx, sy, sz, r, g, b, a\]` (14 floats).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InstanceData {
    /// World position `\[x, y, z\]`.
    pub position: [f64; 3],
    /// Orientation quaternion `\[x, y, z, w\]`.
    pub rotation: [f64; 4],
    /// Scale `\[sx, sy, sz\]`.
    pub scale: [f64; 3],
    /// RGBA color `\[r, g, b, a\]` in `\[0, 1\]`.
    pub color: [f64; 4],
}

impl InstanceData {
    /// Create with default unit scale and white color.
    pub fn new(position: [f64; 3], rotation: [f64; 4]) -> Self {
        Self {
            position,
            rotation,
            scale: [1.0; 3],
            color: [1.0; 4],
        }
    }

    /// Pack to flat 14-float buffer.
    pub fn to_flat14(&self) -> [f64; 14] {
        [
            self.position[0],
            self.position[1],
            self.position[2],
            self.rotation[0],
            self.rotation[1],
            self.rotation[2],
            self.rotation[3],
            self.scale[0],
            self.scale[1],
            self.scale[2],
            self.color[0],
            self.color[1],
            self.color[2],
            self.color[3],
        ]
    }

    /// Pack many instances into a flat `Vec`f64` (14 floats per instance).
    pub fn pack_flat(instances: &[Self]) -> Vec<f64> {
        instances.iter().flat_map(|i| i.to_flat14()).collect()
    }
}

// ===========================================================================
// Shadow map pass configuration
// ===========================================================================

/// Configuration for a single shadow map render pass.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShadowMapPass {
    /// Unique pass identifier.
    pub id: u32,
    /// Resolution of the shadow map texture (width = height = resolution).
    pub resolution: u32,
    /// Light direction (normalized `[x, y, z]`).
    pub light_dir: [f64; 3],
    /// Orthographic projection half-size (world units from center).
    pub ortho_size: f64,
    /// Near clip distance for the shadow camera.
    pub near: f64,
    /// Far clip distance for the shadow camera.
    pub far: f64,
    /// Depth bias to avoid self-shadowing artefacts.
    pub depth_bias: f64,
}

impl ShadowMapPass {
    /// Create a shadow map pass for a directional light.
    pub fn directional(id: u32, light_dir: [f64; 3], ortho_size: f64) -> Self {
        let len = (light_dir[0] * light_dir[0]
            + light_dir[1] * light_dir[1]
            + light_dir[2] * light_dir[2])
            .sqrt();
        let ld = if len > 1e-15 {
            [light_dir[0] / len, light_dir[1] / len, light_dir[2] / len]
        } else {
            [0.0, -1.0, 0.0]
        };
        Self {
            id,
            resolution: 1024,
            light_dir: ld,
            ortho_size,
            near: 0.1,
            far: 200.0,
            depth_bias: 0.005,
        }
    }

    /// Return the pipeline descriptor for this shadow pass.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::shadow_map(format!("shadow_map_{}", self.id))
    }
}

// ===========================================================================
// Deferred shading G-buffer pass
// ===========================================================================

/// G-buffer layout for deferred shading.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GBufferLayout {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Whether to include a velocity/motion-vector render target.
    pub include_velocity: bool,
    /// Whether to include an emissive render target.
    pub include_emissive: bool,
}

impl GBufferLayout {
    /// Create a standard G-buffer layout.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            include_velocity: false,
            include_emissive: false,
        }
    }

    /// Number of render targets in this G-buffer.
    pub fn render_target_count(&self) -> u32 {
        let mut n = 3u32; // albedo, normals, depth
        if self.include_velocity {
            n += 1;
        }
        if self.include_emissive {
            n += 1;
        }
        n
    }

    /// Return the geometry pass pipeline for filling the G-buffer.
    pub fn geometry_pass_pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::opaque("deferred_geometry")
            .with_color_targets(self.render_target_count())
    }

    /// Return the lighting pass pipeline (full-screen quad reading G-buffer).
    pub fn lighting_pass_pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::post_process("deferred_lighting")
    }
}

// ===========================================================================
// Particle rendering pass
// ===========================================================================

/// Configuration for a point-sprite or mesh particle rendering pass.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParticleRenderPass {
    /// Maximum number of particles supported.
    pub max_particles: u32,
    /// Particle size in world units.
    pub particle_size: f64,
    /// Blend mode (usually `Additive` or `AlphaBlend`).
    pub blend_mode: BlendMode,
    /// Whether to sort particles back-to-front (for correct alpha blending).
    pub sort_back_to_front: bool,
}

impl ParticleRenderPass {
    /// Create a default additive particle pass.
    pub fn additive(max_particles: u32) -> Self {
        Self {
            max_particles,
            particle_size: 0.05,
            blend_mode: BlendMode::Additive,
            sort_back_to_front: false,
        }
    }

    /// Create an alpha-blended particle pass (with sorting).
    pub fn alpha_blended(max_particles: u32) -> Self {
        Self {
            max_particles,
            particle_size: 0.05,
            blend_mode: BlendMode::AlphaBlend,
            sort_back_to_front: true,
        }
    }

    /// Pipeline descriptor for this particle pass.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        let mut p = RenderPipelineDesc::transparent("particles");
        p.blend_mode = self.blend_mode;
        p
    }

    /// Sort particle positions by distance to camera (farthest first).
    pub fn sort_particles_by_depth(positions: &[[f64; 3]], camera_pos: [f64; 3]) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..positions.len()).collect();
        indices.sort_by(|&a, &b| {
            let da = dist_sq(positions[a], camera_pos);
            let db = dist_sq(positions[b], camera_pos);
            db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
        });
        indices
    }
}

#[inline]
fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// ===========================================================================
// Post-processing pass
// ===========================================================================

/// A post-processing effect pass (tone-mapping, bloom, FXAA, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PostProcessPass {
    /// Effect name / shader identifier.
    pub effect: String,
    /// Whether this pass is enabled.
    pub enabled: bool,
    /// Arbitrary float parameters (effect-specific).
    pub params: Vec<f64>,
}

impl PostProcessPass {
    /// Create an FXAA anti-aliasing pass.
    pub fn fxaa() -> Self {
        Self {
            effect: "fxaa".to_string(),
            enabled: true,
            params: vec![0.0833, 0.166, 0.75],
        }
    }

    /// Create a bloom pass with threshold and intensity.
    pub fn bloom(threshold: f64, intensity: f64) -> Self {
        Self {
            effect: "bloom".to_string(),
            enabled: true,
            params: vec![threshold, intensity],
        }
    }

    /// Create a tone-mapping pass (Reinhard).
    pub fn tone_map_reinhard(exposure: f64) -> Self {
        Self {
            effect: "tonemap_reinhard".to_string(),
            enabled: true,
            params: vec![exposure],
        }
    }

    /// Create a vignette effect.
    pub fn vignette(strength: f64, radius: f64) -> Self {
        Self {
            effect: "vignette".to_string(),
            enabled: true,
            params: vec![strength, radius],
        }
    }

    /// Pipeline descriptor.
    pub fn pipeline(&self) -> RenderPipelineDesc {
        RenderPipelineDesc::post_process(&self.effect)
    }
}

/// A post-processing pipeline: an ordered list of effect passes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PostProcessPipeline {
    /// Ordered list of effect passes.
    pub passes: Vec<PostProcessPass>,
}

impl PostProcessPipeline {
    /// Create an empty pipeline.
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    /// Add a pass.
    pub fn add_pass(&mut self, pass: PostProcessPass) {
        self.passes.push(pass);
    }

    /// Count of enabled passes.
    pub fn enabled_count(&self) -> usize {
        self.passes.iter().filter(|p| p.enabled).count()
    }

    /// Common preset: FXAA + bloom + tone-mapping.
    pub fn default_hdr() -> Self {
        let mut p = Self::new();
        p.add_pass(PostProcessPass::bloom(0.8, 0.3));
        p.add_pass(PostProcessPass::tone_map_reinhard(1.0));
        p.add_pass(PostProcessPass::fxaa());
        p
    }
}

impl Default for PostProcessPipeline {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// Render frame graph
// ===========================================================================

/// A simple description of a full render frame graph.
///
/// Lists the shadow, geometry, particle, and post-process passes in order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameGraph {
    /// Shadow map passes (one per shadow-casting light).
    pub shadow_passes: Vec<ShadowMapPass>,
    /// G-buffer layout (for deferred shading).
    pub gbuffer: GBufferLayout,
    /// Particle rendering passes.
    pub particle_passes: Vec<ParticleRenderPass>,
    /// Post-processing pipeline.
    pub post_process: PostProcessPipeline,
}

impl FrameGraph {
    /// Create a minimal forward-rendering frame graph.
    pub fn forward(width: u32, height: u32) -> Self {
        Self {
            shadow_passes: vec![ShadowMapPass::directional(0, [0.3, -0.8, 0.5], 20.0)],
            gbuffer: GBufferLayout::new(width, height),
            particle_passes: Vec::new(),
            post_process: PostProcessPipeline::default_hdr(),
        }
    }

    /// Create a deferred-rendering frame graph with particles.
    pub fn deferred(width: u32, height: u32) -> Self {
        let mut fg = Self::forward(width, height);
        fg.gbuffer.include_velocity = true;
        fg.particle_passes.push(ParticleRenderPass::additive(65536));
        fg
    }

    /// Total number of render passes.
    pub fn total_pass_count(&self) -> usize {
        self.shadow_passes.len()
            + 2 // geometry + lighting
            + self.particle_passes.len()
            + self.post_process.enabled_count()
    }
}

// ===========================================================================
// Tests for new renderer types
// ===========================================================================

#[cfg(test)]
mod renderer_ext_tests {

    use crate::renderer::BlendMode;
    use crate::renderer::FrameGraph;
    use crate::renderer::GBufferLayout;
    use crate::renderer::InstanceData;
    use crate::renderer::ParticleRenderPass;
    use crate::renderer::PostProcessPass;
    use crate::renderer::PostProcessPipeline;
    use crate::renderer::RenderPipelineDesc;

    use crate::renderer::ShadowMapPass;

    // --- RenderPipelineDesc ---

    #[test]
    fn test_pipeline_opaque_defaults() {
        let p = RenderPipelineDesc::opaque("test");
        assert!(p.depth_test);
        assert!(p.depth_write);
        assert_eq!(p.blend_mode, BlendMode::Opaque);
        assert_eq!(p.sample_count, 1);
    }

    #[test]
    fn test_pipeline_transparent_no_depth_write() {
        let p = RenderPipelineDesc::transparent("alpha");
        assert!(!p.depth_write);
        assert_eq!(p.blend_mode, BlendMode::AlphaBlend);
    }

    #[test]
    fn test_pipeline_shadow_map_no_color_targets() {
        let p = RenderPipelineDesc::shadow_map("shadow");
        assert_eq!(p.color_target_count, 0);
    }

    #[test]
    fn test_pipeline_post_process_no_depth() {
        let p = RenderPipelineDesc::post_process("pp");
        assert!(!p.depth_test);
        assert!(!p.depth_write);
    }

    #[test]
    fn test_pipeline_with_msaa() {
        let p = RenderPipelineDesc::opaque("msaa").with_msaa(4);
        assert_eq!(p.sample_count, 4);
    }

    #[test]
    fn test_pipeline_with_color_targets() {
        let p = RenderPipelineDesc::opaque("mrt").with_color_targets(4);
        assert_eq!(p.color_target_count, 4);
    }

    #[test]
    fn test_pipeline_msaa_clamped() {
        let p = RenderPipelineDesc::opaque("msaa").with_msaa(16);
        assert_eq!(p.sample_count, 8); // clamped to 8
    }

    // --- InstanceData ---

    #[test]
    fn test_instance_data_flat14_length() {
        let inst = InstanceData::new([1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 1.0]);
        let flat = inst.to_flat14();
        assert_eq!(flat.len(), 14);
    }

    #[test]
    fn test_instance_data_pack_flat() {
        let instances = vec![
            InstanceData::new([0.0; 3], [0.0, 0.0, 0.0, 1.0]),
            InstanceData::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]),
        ];
        let flat = InstanceData::pack_flat(&instances);
        assert_eq!(flat.len(), 28); // 2 × 14
    }

    #[test]
    fn test_instance_data_position_in_flat() {
        let inst = InstanceData::new([5.0, 6.0, 7.0], [0.0, 0.0, 0.0, 1.0]);
        let flat = inst.to_flat14();
        assert!((flat[0] - 5.0).abs() < 1e-10);
        assert!((flat[1] - 6.0).abs() < 1e-10);
        assert!((flat[2] - 7.0).abs() < 1e-10);
    }

    // --- ShadowMapPass ---

    #[test]
    fn test_shadow_map_pass_directional() {
        let p = ShadowMapPass::directional(0, [0.0, -1.0, 0.0], 10.0);
        assert_eq!(p.id, 0);
        assert_eq!(p.resolution, 1024);
        assert!((p.light_dir[1] + 1.0).abs() < 1e-10); // normalized downward
    }

    #[test]
    fn test_shadow_map_pass_normalizes_direction() {
        let p = ShadowMapPass::directional(1, [3.0, 4.0, 0.0], 5.0);
        let len = (p.light_dir[0] * p.light_dir[0]
            + p.light_dir[1] * p.light_dir[1]
            + p.light_dir[2] * p.light_dir[2])
            .sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_shadow_map_pass_pipeline() {
        let p = ShadowMapPass::directional(0, [0.0, -1.0, 0.0], 10.0);
        let pl = p.pipeline();
        assert_eq!(pl.color_target_count, 0);
    }

    // --- GBufferLayout ---

    #[test]
    fn test_gbuffer_render_target_count_basic() {
        let g = GBufferLayout::new(1920, 1080);
        assert_eq!(g.render_target_count(), 3); // albedo + normals + depth
    }

    #[test]
    fn test_gbuffer_render_target_count_with_extras() {
        let mut g = GBufferLayout::new(1920, 1080);
        g.include_velocity = true;
        g.include_emissive = true;
        assert_eq!(g.render_target_count(), 5);
    }

    #[test]
    fn test_gbuffer_geometry_pass_pipeline() {
        let g = GBufferLayout::new(800, 600);
        let p = g.geometry_pass_pipeline();
        assert!(p.depth_test);
    }

    #[test]
    fn test_gbuffer_lighting_pass_pipeline() {
        let g = GBufferLayout::new(800, 600);
        let p = g.lighting_pass_pipeline();
        assert!(!p.depth_test);
    }

    // --- ParticleRenderPass ---

    #[test]
    fn test_particle_pass_additive() {
        let p = ParticleRenderPass::additive(1000);
        assert_eq!(p.blend_mode, BlendMode::Additive);
        assert!(!p.sort_back_to_front);
    }

    #[test]
    fn test_particle_pass_alpha_blended_sorts() {
        let p = ParticleRenderPass::alpha_blended(1000);
        assert_eq!(p.blend_mode, BlendMode::AlphaBlend);
        assert!(p.sort_back_to_front);
    }

    #[test]
    fn test_particle_sort_by_depth() {
        let positions = vec![[1.0, 0.0, 0.0], [5.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let camera = [0.0, 0.0, 0.0];
        let order = ParticleRenderPass::sort_particles_by_depth(&positions, camera);
        // Farthest first: [5, 2, 1] → indices [1, 2, 0]
        assert_eq!(order[0], 1); // dist 5 is farthest
    }

    // --- PostProcessPass ---

    #[test]
    fn test_post_process_fxaa() {
        let p = PostProcessPass::fxaa();
        assert!(p.enabled);
        assert_eq!(p.effect, "fxaa");
    }

    #[test]
    fn test_post_process_bloom_params() {
        let p = PostProcessPass::bloom(0.9, 0.5);
        assert!((p.params[0] - 0.9).abs() < 1e-10);
        assert!((p.params[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_post_process_tone_map() {
        let p = PostProcessPass::tone_map_reinhard(1.5);
        assert_eq!(p.effect, "tonemap_reinhard");
        assert!((p.params[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_post_process_vignette() {
        let p = PostProcessPass::vignette(0.5, 0.8);
        assert_eq!(p.effect, "vignette");
    }

    // --- PostProcessPipeline ---

    #[test]
    fn test_post_process_pipeline_enabled_count() {
        let p = PostProcessPipeline::default_hdr();
        assert_eq!(p.enabled_count(), 3); // bloom + tonemap + fxaa
    }

    #[test]
    fn test_post_process_pipeline_add_pass() {
        let mut p = PostProcessPipeline::new();
        p.add_pass(PostProcessPass::fxaa());
        assert_eq!(p.passes.len(), 1);
    }

    #[test]
    fn test_post_process_pipeline_empty_enabled_count() {
        let p = PostProcessPipeline::new();
        assert_eq!(p.enabled_count(), 0);
    }

    // --- FrameGraph ---

    #[test]
    fn test_frame_graph_forward_pass_count() {
        let fg = FrameGraph::forward(1920, 1080);
        let total = fg.total_pass_count();
        // 1 shadow + 2 deferred + 0 particles + 3 post-process = 6
        assert_eq!(total, 6);
    }

    #[test]
    fn test_frame_graph_deferred_has_particles() {
        let fg = FrameGraph::deferred(1920, 1080);
        assert!(!fg.particle_passes.is_empty());
    }

    #[test]
    fn test_frame_graph_deferred_has_velocity() {
        let fg = FrameGraph::deferred(1920, 1080);
        assert!(fg.gbuffer.include_velocity);
    }

    #[test]
    fn test_frame_graph_forward_one_shadow() {
        let fg = FrameGraph::forward(800, 600);
        assert_eq!(fg.shadow_passes.len(), 1);
    }

    #[test]
    fn test_render_pipeline_desc_default() {
        let p = RenderPipelineDesc::default();
        assert_eq!(p.label, "default");
        assert_eq!(p.blend_mode, BlendMode::Opaque);
    }

    #[test]
    fn test_post_process_pipeline_default() {
        let p = PostProcessPipeline::default();
        assert_eq!(p.passes.len(), 0);
    }
}
