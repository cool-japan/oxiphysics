// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core renderer hint types: AABB, contact points, velocity vectors,
//! per-body hints, aggregate frame hints, and canvas draw commands.

#![allow(missing_docs)]
#![allow(dead_code)]

use crate::engine::WasmPhysicsEngine;
use crate::types::ContactResult;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Aabb
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
// RendererHints
// ---------------------------------------------------------------------------

/// Aggregate per-frame renderer hints collected from a `WasmPhysicsEngine`.
///
/// Serializable to JSON for transfer to the main thread renderer.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RendererHints {
    /// Per-body rendering information.
    // Vec<BodyRenderHint> is not IntoWasmAbi for pub fields — expose via method.
    #[wasm_bindgen(skip)]
    pub body_hints: Vec<BodyRenderHint>,
    /// Contact point hints for this frame.
    // Vec<ContactPointHint> is not IntoWasmAbi for pub fields.
    #[wasm_bindgen(skip)]
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

    /// Return all AABB data as a flat `Vec<f64>`: `[handle, min_x, min_y, min_z, max_x, max_y, max_z, ...]`.
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

    /// Return all contact positions as a flat `Vec<f64>`: `[x0, y0, z0, x1, ...]`.
    pub fn contact_positions_flat(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.contact_points.len() * 3);
        for cp in &self.contact_points {
            out.extend_from_slice(&cp.position);
        }
        out
    }

    /// Return all velocity vectors as a flat `Vec<f64>`: `[ox, oy, oz, vx, vy, vz, ...]`.
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
// RendererHints — wasm_bindgen JS wrappers
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl RendererHints {
    /// Construct `RendererHints` from a `WasmPhysicsEngine` (JS-friendly).
    #[wasm_bindgen(js_name = "fromEngine")]
    pub fn from_engine_js(
        engine: &WasmPhysicsEngine,
        include_aabbs: bool,
        include_velocities: bool,
        include_contacts: bool,
    ) -> RendererHints {
        RendererHints::from_engine(engine, include_aabbs, include_velocities, include_contacts)
    }

    /// Serialize to a JSON string (JS-friendly).
    #[wasm_bindgen(js_name = "toJson")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Return all AABB data as a flat `Float64Array`-ready `Vec<f64>`.
    ///
    /// Layout: `[handle, min_x, min_y, min_z, max_x, max_y, max_z, ...]`  (7 values per AABB).
    #[wasm_bindgen(js_name = "aabbsFlat")]
    pub fn aabbs_flat_js(&self) -> Vec<f64> {
        self.aabbs_flat()
    }

    /// Return all contact positions as a flat `Float64Array`-ready `Vec<f64>`.
    ///
    /// Layout: `[x0, y0, z0, x1, y1, z1, ...]`  (3 values per contact).
    #[wasm_bindgen(js_name = "contactPositionsFlat")]
    pub fn contact_positions_flat_js(&self) -> Vec<f64> {
        self.contact_positions_flat()
    }

    /// Return all velocity vectors as a flat `Float64Array`-ready `Vec<f64>`.
    ///
    /// Layout: `[ox, oy, oz, vx, vy, vz, ...]`  (6 values per body).
    #[wasm_bindgen(js_name = "velocityVectorsFlat")]
    pub fn velocity_vectors_flat_js(&self) -> Vec<f64> {
        self.velocity_vectors_flat()
    }

    /// Return `body_hints` serialized as a `JsValue` (JSON array).
    #[wasm_bindgen(js_name = "bodyHintsJs")]
    pub fn body_hints_js(&self) -> Result<JsValue, JsValue> {
        use crate::wasm_helpers::to_js_value;
        to_js_value(&self.body_hints)
    }

    /// Return `contact_points` serialized as a `JsValue` (JSON array).
    #[wasm_bindgen(js_name = "contactPointsJs")]
    pub fn contact_points_js(&self) -> Result<JsValue, JsValue> {
        use crate::wasm_helpers::to_js_value;
        to_js_value(&self.contact_points)
    }

    /// Build canvas 2D draw commands from these hints (JS-friendly).
    ///
    /// Returns a `JsValue` containing a JSON array of draw commands.
    #[wasm_bindgen(js_name = "buildDrawCommandsJs")]
    pub fn build_draw_commands_js(&self, scale: f64) -> Result<JsValue, JsValue> {
        use crate::wasm_helpers::to_js_value;
        let cmds = build_draw_commands(self, scale);
        to_js_value(&cmds)
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
        cx: f64,
        cy: f64,
        radius: f64,
        color: String,
        line_width: f64,
    },
    /// Draw a line between two points.
    Line {
        x0: f64,
        y0: f64,
        x1: f64,
        y1: f64,
        color: String,
        line_width: f64,
    },
    /// Draw an axis-aligned rectangle.
    Rect {
        x: f64,
        y: f64,
        width: f64,
        height: f64,
        color: String,
        line_width: f64,
    },
    /// Draw a text label.
    Text {
        x: f64,
        y: f64,
        text: String,
        color: String,
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

        cmds.push(DebugDrawCommand::Text {
            x: cx + 4.0,
            y: cy - 4.0,
            text: format!("{}", hint.handle),
            color: "#ffffff".to_string(),
            font_size: 10.0,
        });
    }

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
        assert_eq!(hints.velocity_vectors_flat().len(), 6);
    }

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
        for h in &hints.body_hints {
            assert!(h.aabb.is_none());
        }
    }

    // ---------------------------------------------------------------------------
    // Integration tests (Slice W6)
    // ---------------------------------------------------------------------------

    /// Verify that RendererHints serializes all body positions correctly via flat API.
    #[test]
    fn test_renderer_hints_velocity_flat_matches_struct() {
        let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
        engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
        let hints = RendererHints::from_engine(&engine, false, true, false);
        let flat = hints.velocity_vectors_flat();
        // origin = position of body = [1, 2, 3]
        assert!((flat[0] - 1.0).abs() < 1e-10);
        assert!((flat[1] - 2.0).abs() < 1e-10);
        assert!((flat[2] - 3.0).abs() < 1e-10);
    }

    /// Verify that DebugDrawCommand serializes to tagged JSON correctly.
    #[test]
    fn test_debug_draw_command_json_tag() {
        let cmd = DebugDrawCommand::Circle {
            cx: 1.0,
            cy: 2.0,
            radius: 5.0,
            color: "#ff0000".to_string(),
            line_width: 1.0,
        };
        let json = serde_json::to_string(&cmd).expect("serialize");
        assert!(
            json.contains("\"type\":\"Circle\""),
            "expected type tag in {json}"
        );
    }
}
