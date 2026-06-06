// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the trigger/sensor volume system.
//!
//! Wraps `oxiphysics::trigger::TriggerWorld` and exposes a JSON-oriented
//! surface suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// WasmBodyEntry — WASM-friendly representation of a physics body
// ---------------------------------------------------------------------------

/// Minimal description of a physics body for overlap testing.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmBodyEntry {
    /// Unique body index (stable per simulation).
    pub index: usize,
    /// World-space position `[x, y, z]`.
    #[wasm_bindgen(skip)]
    pub position: [f64; 3],
    /// Approximate bounding radius (metres).
    pub radius: f64,
    /// Whether the body is in the sleeping state.
    pub is_sleeping: bool,
}

impl WasmBodyEntry {
    /// Create a new body entry.
    pub fn new(index: usize, position: [f64; 3], radius: f64, is_sleeping: bool) -> Self {
        WasmBodyEntry {
            index,
            position,
            radius,
            is_sleeping,
        }
    }
}

// ---------------------------------------------------------------------------
// WasmBodyEntry — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmBodyEntry {
    /// Construct a new body entry from explicit components (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(
        index: usize,
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        is_sleeping: bool,
    ) -> WasmBodyEntry {
        WasmBodyEntry::new(index, [x, y, z], radius, is_sleeping)
    }

    /// Position as a flat `Vec<f64>` of length 3.
    #[wasm_bindgen(js_name = "get_position")]
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Replace the position with `[x, y, z]`.
    #[wasm_bindgen(js_name = "set_position")]
    pub fn set_position_js(&mut self, x: f64, y: f64, z: f64) {
        self.position = [x, y, z];
    }

    /// Serialise to a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// WasmTriggerEvent — serialisable trigger event
// ---------------------------------------------------------------------------

/// A trigger overlap event emitted by [`WasmTriggerWorld::update`].
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmTriggerEvent {
    /// `"enter"`, `"exit"`, or `"stay"`.
    #[wasm_bindgen(skip)]
    pub kind: String,
    /// ID of the trigger volume that fired.
    pub trigger_id: u32,
    /// Index of the body that triggered the event.
    pub body_index: usize,
    /// Simulation step at which the event was recorded.
    pub step: u64,
}

// ---------------------------------------------------------------------------
// WasmTriggerEvent — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmTriggerEvent {
    /// Event kind (`"enter"`, `"exit"`, or `"stay"`).
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }

    /// Replace the event kind string.
    #[wasm_bindgen(setter)]
    pub fn set_kind(&mut self, kind: String) {
        self.kind = kind;
    }

    /// Serialise the event as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// Internal geometry helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum TriggerShape {
    Sphere {
        center: [f64; 3],
        radius: f64,
    },
    Aabb {
        min: [f64; 3],
        max: [f64; 3],
    },
    Capsule {
        start: [f64; 3],
        end: [f64; 3],
        radius: f64,
    },
}

impl TriggerShape {
    fn contains_point_with_radius(&self, p: [f64; 3], r: f64) -> bool {
        match self {
            TriggerShape::Sphere { center, radius } => {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                let threshold = radius + r;
                dist_sq <= threshold * threshold
            }
            TriggerShape::Aabb { min, max } => {
                (p[0] - r) <= max[0]
                    && (p[0] + r) >= min[0]
                    && (p[1] - r) <= max[1]
                    && (p[1] + r) >= min[1]
                    && (p[2] - r) <= max[2]
                    && (p[2] + r) >= min[2]
            }
            TriggerShape::Capsule { start, end, radius } => {
                // Closest point on line segment [start, end] to p, then sphere check.
                let ab = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
                let ap = [p[0] - start[0], p[1] - start[1], p[2] - start[2]];
                let ab_len_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
                let t = if ab_len_sq < 1e-30 {
                    0.0
                } else {
                    let dot = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2];
                    (dot / ab_len_sq).clamp(0.0, 1.0)
                };
                let closest = [
                    start[0] + t * ab[0],
                    start[1] + t * ab[1],
                    start[2] + t * ab[2],
                ];
                let dx = p[0] - closest[0];
                let dy = p[1] - closest[1];
                let dz = p[2] - closest[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                let threshold = radius + r;
                dist_sq <= threshold * threshold
            }
        }
    }
}

#[derive(Clone, Debug)]
struct TriggerVolume {
    id: u32,
    shape: TriggerShape,
    enabled: bool,
}

// ---------------------------------------------------------------------------
// WasmTriggerWorld
// ---------------------------------------------------------------------------

/// WASM wrapper around the trigger/sensor volume system.
///
/// Manages spherical, AABB, and capsule trigger volumes and detects body
/// enter/exit/stay events each simulation step.
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmTriggerWorld {
    volumes: Vec<TriggerVolume>,
    occupied: std::collections::HashMap<usize, Vec<u32>>,
    next_id: u32,
    emit_stay: bool,
}

impl WasmTriggerWorld {
    /// Create a new trigger world with no volumes.
    pub fn new() -> Self {
        WasmTriggerWorld {
            volumes: Vec::new(),
            occupied: std::collections::HashMap::new(),
            next_id: 0,
            emit_stay: false,
        }
    }

    /// Add a spherical trigger volume. Returns the assigned volume ID.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64, _tags: Vec<String>) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Sphere { center, radius },
            enabled: true,
        });
        id
    }

    /// Add an AABB trigger volume. Returns the assigned volume ID.
    pub fn add_aabb(&mut self, min: [f64; 3], max: [f64; 3], _tags: Vec<String>) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Aabb { min, max },
            enabled: true,
        });
        id
    }

    /// Add a capsule trigger volume. Returns the assigned volume ID.
    pub fn add_capsule(
        &mut self,
        start: [f64; 3],
        end: [f64; 3],
        radius: f64,
        _tags: Vec<String>,
    ) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Capsule { start, end, radius },
            enabled: true,
        });
        id
    }

    /// Remove a trigger volume by ID. Also clears any occupancy records for it.
    pub fn remove_volume(&mut self, id: u32) {
        self.volumes.retain(|v| v.id != id);
        for occ in self.occupied.values_mut() {
            occ.retain(|&t| t != id);
        }
        self.occupied.retain(|_, ids| !ids.is_empty());
    }

    /// Enable or disable stay-event emission.
    pub fn set_emit_stay(&mut self, emit: bool) {
        self.emit_stay = emit;
    }

    /// Return a sorted list of volume IDs currently occupied by `body_index`.
    pub fn bodies_in(&self, id: u32) -> Vec<usize> {
        self.occupied
            .iter()
            .filter_map(|(&bi, ids)| if ids.contains(&id) { Some(bi) } else { None })
            .collect()
    }

    /// Number of non-empty occupancy entries across all bodies.
    pub fn active_occupancies(&self) -> usize {
        self.occupied.values().filter(|v| !v.is_empty()).count()
    }

    /// Update trigger occupancy for the given body snapshot and return events.
    ///
    /// `step` is the current simulation tick (monotonically increasing).
    pub fn update(&mut self, step: u64, bodies: &[WasmBodyEntry]) -> Vec<WasmTriggerEvent> {
        let mut events = Vec::new();

        // Step 1: bodies absent this step → Exit events
        let present: std::collections::HashSet<usize> = bodies.iter().map(|b| b.index).collect();
        let absent_exits: Vec<(usize, u32)> = self
            .occupied
            .iter()
            .filter(|&(&bi, _)| !present.contains(&bi))
            .flat_map(|(&bi, ids)| ids.iter().map(move |&tid| (bi, tid)))
            .collect();
        for (body_index, trigger_id) in absent_exits {
            events.push(WasmTriggerEvent {
                kind: "exit".to_string(),
                trigger_id,
                body_index,
                step,
            });
            if let Some(occ) = self.occupied.get_mut(&body_index) {
                occ.retain(|&t| t != trigger_id);
            }
        }
        self.occupied.retain(|_, ids| !ids.is_empty());

        // Step 2: compute overlaps (read-only on volumes)
        struct OverlapTest {
            body_index: usize,
            trigger_id: u32,
            now_inside: bool,
        }
        let tests: Vec<OverlapTest> = self
            .volumes
            .iter()
            .filter(|v| v.enabled)
            .flat_map(|vol| {
                let vid = vol.id;
                let shape = vol.shape.clone();
                bodies.iter().map(move |body| OverlapTest {
                    body_index: body.index,
                    trigger_id: vid,
                    now_inside: shape.contains_point_with_radius(body.position, body.radius),
                })
            })
            .collect();

        // Step 3: apply and emit
        for test in tests {
            let was_inside = self
                .occupied
                .get(&test.body_index)
                .is_some_and(|ids| ids.contains(&test.trigger_id));
            match (was_inside, test.now_inside) {
                (false, true) => {
                    self.occupied
                        .entry(test.body_index)
                        .or_default()
                        .push(test.trigger_id);
                    events.push(WasmTriggerEvent {
                        kind: "enter".to_string(),
                        trigger_id: test.trigger_id,
                        body_index: test.body_index,
                        step,
                    });
                }
                (true, false) => {
                    if let Some(occ) = self.occupied.get_mut(&test.body_index) {
                        occ.retain(|&t| t != test.trigger_id);
                    }
                    events.push(WasmTriggerEvent {
                        kind: "exit".to_string(),
                        trigger_id: test.trigger_id,
                        body_index: test.body_index,
                        step,
                    });
                }
                (true, true) => {
                    if self.emit_stay {
                        events.push(WasmTriggerEvent {
                            kind: "stay".to_string(),
                            trigger_id: test.trigger_id,
                            body_index: test.body_index,
                            step,
                        });
                    }
                }
                (false, false) => {}
            }
        }

        events
    }

    /// Serialise the list of trigger events as a JSON string.
    pub fn events_to_json(events: &[WasmTriggerEvent]) -> String {
        serde_json::to_string(events).unwrap_or_default()
    }

    /// Return the number of volumes registered in this world.
    pub fn volume_count(&self) -> usize {
        self.volumes.len()
    }
}

impl Default for WasmTriggerWorld {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// WasmTriggerWorld — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmTriggerWorld {
    /// Construct an empty trigger world (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmTriggerWorld {
        WasmTriggerWorld::new()
    }

    /// Add a spherical trigger volume. Tags are passed as a JSON string array.
    ///
    /// Returns the assigned volume ID.
    #[wasm_bindgen(js_name = "add_sphere")]
    pub fn add_sphere_js(
        &mut self,
        cx: f64,
        cy: f64,
        cz: f64,
        radius: f64,
        tags_json: &str,
    ) -> Result<u32, JsValue> {
        let tags: Vec<String> = if tags_json.is_empty() {
            Vec::new()
        } else {
            serde_json::from_str(tags_json).map_err(err_to_jsvalue)?
        };
        Ok(self.add_sphere([cx, cy, cz], radius, tags))
    }

    /// Add an AABB trigger volume. Tags are passed as a JSON string array.
    #[wasm_bindgen(js_name = "add_aabb")]
    pub fn add_aabb_js(
        &mut self,
        min_x: f64,
        min_y: f64,
        min_z: f64,
        max_x: f64,
        max_y: f64,
        max_z: f64,
        tags_json: &str,
    ) -> Result<u32, JsValue> {
        let tags: Vec<String> = if tags_json.is_empty() {
            Vec::new()
        } else {
            serde_json::from_str(tags_json).map_err(err_to_jsvalue)?
        };
        Ok(self.add_aabb([min_x, min_y, min_z], [max_x, max_y, max_z], tags))
    }

    /// Add a capsule trigger volume. Tags are passed as a JSON string array.
    #[wasm_bindgen(js_name = "add_capsule")]
    pub fn add_capsule_js(
        &mut self,
        sx: f64,
        sy: f64,
        sz: f64,
        ex: f64,
        ey: f64,
        ez: f64,
        radius: f64,
        tags_json: &str,
    ) -> Result<u32, JsValue> {
        let tags: Vec<String> = if tags_json.is_empty() {
            Vec::new()
        } else {
            serde_json::from_str(tags_json).map_err(err_to_jsvalue)?
        };
        Ok(self.add_capsule([sx, sy, sz], [ex, ey, ez], radius, tags))
    }

    /// Remove a trigger volume by ID.
    #[wasm_bindgen(js_name = "remove_volume")]
    pub fn remove_volume_js(&mut self, id: u32) {
        self.remove_volume(id);
    }

    /// Enable or disable stay-event emission.
    #[wasm_bindgen(js_name = "set_emit_stay")]
    pub fn set_emit_stay_js(&mut self, emit: bool) {
        self.set_emit_stay(emit);
    }

    /// Bodies currently occupying volume `id`, returned as a JSON array.
    #[wasm_bindgen(js_name = "bodies_in")]
    pub fn bodies_in_js(&self, id: u32) -> Result<String, JsValue> {
        let v = self.bodies_in(id);
        serde_json::to_string(&v).map_err(err_to_jsvalue)
    }

    /// Number of non-empty occupancy entries.
    #[wasm_bindgen(js_name = "active_occupancies")]
    pub fn active_occupancies_js(&self) -> usize {
        self.active_occupancies()
    }

    /// Number of registered trigger volumes.
    #[wasm_bindgen(js_name = "volume_count")]
    pub fn volume_count_js(&self) -> usize {
        self.volume_count()
    }

    /// Update occupancy for a JSON-encoded body snapshot and return JSON events.
    ///
    /// `bodies_json` is a JSON array of [`WasmBodyEntry`] objects.
    #[wasm_bindgen(js_name = "update_json")]
    pub fn update_json_js(&mut self, step: u64, bodies_json: &str) -> Result<String, JsValue> {
        let bodies: Vec<WasmBodyEntry> =
            serde_json::from_str(bodies_json).map_err(err_to_jsvalue)?;
        let events = self.update(step, &bodies);
        serde_json::to_string(&events).map_err(err_to_jsvalue)
    }

    /// Update occupancy for a JSON-encoded body snapshot and return events as a `JsValue`.
    #[wasm_bindgen(js_name = "update_js_value")]
    pub fn update_js_value(&mut self, step: u64, bodies_json: &str) -> Result<JsValue, JsValue> {
        let bodies: Vec<WasmBodyEntry> =
            serde_json::from_str(bodies_json).map_err(err_to_jsvalue)?;
        let events = self.update(step, &bodies);
        to_js_value(&events)
    }
}

// ---------------------------------------------------------------------------
// Free wasm-bindgen helpers
// ---------------------------------------------------------------------------

/// Serialise a JSON-encoded list of [`WasmTriggerEvent`] back to a JSON string.
///
/// This mirrors the static [`WasmTriggerWorld::events_to_json`] helper but is
/// callable from JS as a free function.
#[wasm_bindgen(js_name = "trigger_events_to_json")]
pub fn trigger_events_to_json_js(events_json: &str) -> Result<String, JsValue> {
    let events: Vec<WasmTriggerEvent> =
        serde_json::from_str(events_json).map_err(err_to_jsvalue)?;
    serde_json::to_string(&events).map_err(err_to_jsvalue)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn body(index: usize, pos: [f64; 3]) -> WasmBodyEntry {
        WasmBodyEntry::new(index, pos, 0.5, false)
    }

    #[test]
    fn test_triggers_bridge_instantiation() {
        let world = WasmTriggerWorld::new();
        let json = WasmTriggerWorld::events_to_json(&[]);
        assert!(
            !json.is_empty(),
            "events_to_json should return non-empty string"
        );
        assert_eq!(world.volume_count(), 0);
    }

    #[test]
    fn test_triggers_bridge_sphere_enter_exit() {
        let mut world = WasmTriggerWorld::new();
        let id = world.add_sphere([0.0, 0.0, 0.0], 5.0, vec!["zone".into()]);

        let inside = vec![body(0, [0.0, 0.0, 0.0])];
        let ev = world.update(1, &inside);
        assert_eq!(ev.len(), 1);
        assert_eq!(ev[0].kind, "enter");
        assert_eq!(ev[0].body_index, 0);

        // Step 2: body still inside — no events without emit_stay
        let ev2 = world.update(2, &inside);
        assert!(ev2.is_empty());

        // Step 3: body moves outside
        let outside = vec![body(0, [20.0, 0.0, 0.0])];
        let ev3 = world.update(3, &outside);
        assert_eq!(ev3.len(), 1);
        assert_eq!(ev3[0].kind, "exit");
        assert_eq!(world.bodies_in(id).len(), 0);
    }

    #[test]
    fn test_triggers_bridge_stay_events() {
        let mut world = WasmTriggerWorld::new();
        world.set_emit_stay(true);
        world.add_sphere([0.0, 0.0, 0.0], 5.0, vec![]);
        let bodies = vec![body(0, [0.0, 0.0, 0.0])];
        world.update(1, &bodies); // Enter
        let ev = world.update(2, &bodies); // Stay
        assert!(ev.iter().any(|e| e.kind == "stay"));
    }

    #[test]
    fn test_triggers_bridge_aabb() {
        let mut world = WasmTriggerWorld::new();
        world.add_aabb([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], vec![]);
        let ev = world.update(1, &[body(0, [0.0, 0.0, 0.0])]);
        assert!(ev.iter().any(|e| e.kind == "enter"));
        let ev2 = world.update(2, &[body(0, [100.0, 0.0, 0.0])]);
        assert!(ev2.iter().any(|e| e.kind == "exit"));
    }

    #[test]
    fn test_triggers_bridge_remove_volume() {
        let mut world = WasmTriggerWorld::new();
        let id = world.add_sphere([0.0, 0.0, 0.0], 5.0, vec![]);
        world.update(1, &[body(0, [0.0, 0.0, 0.0])]);
        assert_eq!(world.active_occupancies(), 1);
        world.remove_volume(id);
        assert_eq!(world.active_occupancies(), 0);
    }
}
