// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Trigger / sensor volume system.
//!
//! A `TriggerVolume` is a ghost shape (sphere, AABB, or capsule) that detects
//! when physics bodies enter or leave its interior.  Unlike rigid colliders,
//! trigger volumes produce no contact forces — they only fire `TriggerEvent`s.
//!
//! The `TriggerWorld` manages a collection of trigger volumes and updates their
//! occupancy state each simulation step.
//!
//! > **Integration note**: `TriggerEvent` variants mirror the
//! > `PhysicsEvent::BodyEnterRegion` / `PhysicsEvent::BodyLeaveRegion` variants
//! > in [`crate::event_bus`].  You can forward trigger events to an `EventBus`
//! > to unify game logic.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::trigger::{TriggerWorld, BodyEntry};
//!
//! let mut world = TriggerWorld::new();
//! let id = world.add_sphere([0.0, 0.0, 0.0], 5.0, vec!["hazard".into()]);
//!
//! let bodies = vec![
//!     BodyEntry { index: 0, position: [0.0, 0.0, 0.0], radius: 0.5, is_sleeping: false },
//! ];
//! let events = world.update(1, &bodies);
//! // body 0 just entered the sphere trigger
//! assert_eq!(events.len(), 1);
//! assert!(events[0].is_enter());
//! ```

use std::collections::HashMap;
use std::collections::HashSet;

// ============================================================================
// Vec3 helpers (self-contained)
// ============================================================================

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_len2(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ============================================================================
// TriggerShape
// ============================================================================

/// Shape defining the sensing region of a [`TriggerVolume`].
#[derive(Clone, Debug)]
pub enum TriggerShape {
    /// Spherical trigger.
    Sphere {
        /// World-space center.
        center: [f64; 3],
        /// Sphere radius.
        radius: f64,
    },
    /// Axis-aligned bounding-box trigger.
    Aabb {
        /// Minimum corner (inclusive).
        min: [f64; 3],
        /// Maximum corner (inclusive).
        max: [f64; 3],
    },
    /// Capsule trigger (cylinder capped with hemispheres).
    Capsule {
        /// Center of the bottom hemisphere.
        start: [f64; 3],
        /// Center of the top hemisphere.
        end: [f64; 3],
        /// Capsule radius.
        radius: f64,
    },
}

impl TriggerShape {
    /// Returns `true` if the point `p` with bounding-sphere radius `r` overlaps
    /// this trigger shape.
    pub fn contains_point_with_radius(&self, p: [f64; 3], r: f64) -> bool {
        match self {
            TriggerShape::Sphere { center, radius } => {
                let d2 = v3_len2(v3_sub(p, *center));
                let combined = radius + r;
                d2 <= combined * combined
            }
            TriggerShape::Aabb { min, max } => {
                (p[0] >= min[0] - r && p[0] <= max[0] + r)
                    && (p[1] >= min[1] - r && p[1] <= max[1] + r)
                    && (p[2] >= min[2] - r && p[2] <= max[2] + r)
            }
            TriggerShape::Capsule { start, end, radius } => {
                let ab = v3_sub(*end, *start);
                let ap = v3_sub(p, *start);
                let t = (v3_dot(ap, ab) / v3_dot(ab, ab).max(1e-12)).clamp(0.0, 1.0);
                let closest = v3_add(*start, v3_scale(ab, t));
                let d2 = v3_len2(v3_sub(p, closest));
                let combined = radius + r;
                d2 <= combined * combined
            }
        }
    }

    /// Approximate bounding sphere radius (for broad-phase pre-rejection).
    pub fn bounding_radius(&self) -> f64 {
        match self {
            TriggerShape::Sphere { radius, .. } => *radius,
            TriggerShape::Aabb { min, max } => {
                let d = v3_sub(*max, *min);
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() * 0.5
            }
            TriggerShape::Capsule { start, end, radius } => {
                v3_len2(v3_sub(*end, *start)).sqrt() * 0.5 + radius
            }
        }
    }
}

// ============================================================================
// TriggerVolume
// ============================================================================

/// A named ghost shape that detects body overlap without generating contact
/// forces.
#[derive(Clone, Debug)]
pub struct TriggerVolume {
    /// Unique identifier assigned by [`TriggerWorld`].
    pub id: u32,
    /// Geometric shape of the sensor region.
    pub shape: TriggerShape,
    /// Arbitrary string tags (e.g. `"hazard"`, `"checkpoint"`, `"zone_a"`).
    pub tags: Vec<String>,
    /// Whether this trigger is active.  Disabled triggers do not emit events.
    pub enabled: bool,
}

impl TriggerVolume {
    /// Returns `true` if any tag matches `tag`.
    pub fn has_tag(&self, tag: &str) -> bool {
        self.tags.iter().any(|t| t == tag)
    }
}

// ============================================================================
// BodyEntry — caller-supplied body state for one update tick
// ============================================================================

/// Minimal body descriptor required by [`TriggerWorld::update`].
///
/// The caller populates this from their physics state each step.
#[derive(Clone, Debug)]
pub struct BodyEntry {
    /// Stable index identifying this body.
    pub index: usize,
    /// World-space center position.
    pub position: [f64; 3],
    /// Bounding-sphere radius used for overlap testing.
    pub radius: f64,
    /// Whether the body is sleeping.  Sleeping bodies still participate in
    /// trigger tests.
    pub is_sleeping: bool,
}

// ============================================================================
// TriggerEvent
// ============================================================================

/// Events produced by [`TriggerWorld::update`].
#[derive(Clone, Debug)]
pub enum TriggerEvent {
    /// A body entered a trigger volume this step.
    Enter {
        /// The trigger that was entered.
        trigger_id: u32,
        /// The body that entered.
        body_index: usize,
        /// Simulation step at which the event occurred.
        step: u64,
    },
    /// A body exited a trigger volume this step.
    Exit {
        /// The trigger that was exited.
        trigger_id: u32,
        /// The body that exited.
        body_index: usize,
        /// Simulation step at which the event occurred.
        step: u64,
    },
    /// A body is persistently inside a trigger (emitted only when
    /// [`TriggerWorld::set_emit_stay`] is `true`).
    Stay {
        /// The containing trigger.
        trigger_id: u32,
        /// The overlapping body.
        body_index: usize,
        /// Current simulation step.
        step: u64,
    },
}

impl TriggerEvent {
    /// Returns the trigger id.
    pub fn trigger_id(&self) -> u32 {
        match self {
            TriggerEvent::Enter { trigger_id, .. }
            | TriggerEvent::Exit { trigger_id, .. }
            | TriggerEvent::Stay { trigger_id, .. } => *trigger_id,
        }
    }

    /// Returns the body index.
    pub fn body_index(&self) -> usize {
        match self {
            TriggerEvent::Enter { body_index, .. }
            | TriggerEvent::Exit { body_index, .. }
            | TriggerEvent::Stay { body_index, .. } => *body_index,
        }
    }

    /// Returns `true` if this is an `Enter` event.
    pub fn is_enter(&self) -> bool {
        matches!(self, TriggerEvent::Enter { .. })
    }

    /// Returns `true` if this is an `Exit` event.
    pub fn is_exit(&self) -> bool {
        matches!(self, TriggerEvent::Exit { .. })
    }

    /// Returns `true` if this is a `Stay` event.
    pub fn is_stay(&self) -> bool {
        matches!(self, TriggerEvent::Stay { .. })
    }
}

// ============================================================================
// Overlap result (internal scratch type)
// ============================================================================

struct OverlapTest {
    body_index: usize,
    trigger_id: u32,
    now_inside: bool,
}

// ============================================================================
// TriggerWorld
// ============================================================================

/// Manages trigger volumes and tracks body occupancy across simulation steps.
///
/// Call [`update`](TriggerWorld::update) each step with current body positions.
/// The returned events describe transitions (enter / exit) and, optionally,
/// persistent overlaps (stay).
///
/// Occupancy is accumulated across calls — no need to reset between steps.
///
/// ## Example
///
/// ```rust
/// use oxiphysics::trigger::{TriggerWorld, BodyEntry};
///
/// let mut world = TriggerWorld::new();
/// let _id = world.add_aabb([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], vec![]);
///
/// let bodies = vec![BodyEntry { index: 0, position: [0.0, 0.0, 0.0], radius: 0.5, is_sleeping: false }];
/// let ev = world.update(1, &bodies);
/// assert_eq!(ev.len(), 1); // Enter
/// let ev2 = world.update(2, &bodies);
/// assert!(ev2.is_empty()); // Stay suppressed by default
/// ```
#[derive(Debug, Default)]
pub struct TriggerWorld {
    volumes: Vec<TriggerVolume>,
    next_id: u32,
    /// `body_index` → list of trigger ids currently containing that body.
    occupied: HashMap<usize, Vec<u32>>,
    /// When `true`, persistent overlaps emit [`TriggerEvent::Stay`] each step.
    emit_stay: bool,
}

impl TriggerWorld {
    /// Creates an empty [`TriggerWorld`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Enable or disable `Stay` event emission (default: `false`).
    pub fn set_emit_stay(&mut self, emit: bool) {
        self.emit_stay = emit;
    }

    // -----------------------------------------------------------------------
    // Volume registration
    // -----------------------------------------------------------------------

    /// Adds a spherical trigger.  Returns its id.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64, tags: Vec<String>) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Sphere { center, radius },
            tags,
            enabled: true,
        });
        id
    }

    /// Adds an AABB trigger.  Returns its id.
    pub fn add_aabb(&mut self, min: [f64; 3], max: [f64; 3], tags: Vec<String>) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Aabb { min, max },
            tags,
            enabled: true,
        });
        id
    }

    /// Adds a capsule trigger.  Returns its id.
    pub fn add_capsule(
        &mut self,
        start: [f64; 3],
        end: [f64; 3],
        radius: f64,
        tags: Vec<String>,
    ) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.volumes.push(TriggerVolume {
            id,
            shape: TriggerShape::Capsule { start, end, radius },
            tags,
            enabled: true,
        });
        id
    }

    /// Removes a trigger volume by id.  Returns `true` if it existed.
    ///
    /// Any occupancy entries referencing the removed volume are cleared.
    pub fn remove_volume(&mut self, id: u32) -> bool {
        if let Some(pos) = self.volumes.iter().position(|v| v.id == id) {
            self.volumes.swap_remove(pos);
            for occ in self.occupied.values_mut() {
                occ.retain(|&tid| tid != id);
            }
            self.occupied.retain(|_, ids| !ids.is_empty());
            true
        } else {
            false
        }
    }

    /// Enables a previously disabled trigger volume.
    pub fn enable_volume(&mut self, id: u32) {
        if let Some(v) = self.volumes.iter_mut().find(|v| v.id == id) {
            v.enabled = true;
        }
    }

    /// Disables a trigger volume.  Disabled volumes do not emit events.
    pub fn disable_volume(&mut self, id: u32) {
        if let Some(v) = self.volumes.iter_mut().find(|v| v.id == id) {
            v.enabled = false;
        }
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Returns body indices currently inside the given trigger.
    pub fn bodies_in(&self, trigger_id: u32) -> Vec<usize> {
        self.occupied
            .iter()
            .filter_map(|(&bi, ids)| {
                if ids.contains(&trigger_id) {
                    Some(bi)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Returns trigger ids currently containing the given body.
    pub fn triggers_containing(&self, body_index: usize) -> Vec<u32> {
        self.occupied.get(&body_index).cloned().unwrap_or_default()
    }

    /// Number of registered trigger volumes.
    pub fn volume_count(&self) -> usize {
        self.volumes.len()
    }

    /// Total number of active body/trigger occupancy pairs.
    pub fn active_occupancies(&self) -> usize {
        self.occupied.values().map(|v| v.len()).sum()
    }

    /// Returns a reference to the volume with the given id.
    pub fn volume(&self, id: u32) -> Option<&TriggerVolume> {
        self.volumes.iter().find(|v| v.id == id)
    }

    /// Returns all volumes carrying the specified tag.
    pub fn volumes_with_tag(&self, tag: &str) -> Vec<&TriggerVolume> {
        self.volumes.iter().filter(|v| v.has_tag(tag)).collect()
    }

    // -----------------------------------------------------------------------
    // Update
    // -----------------------------------------------------------------------

    /// Advances trigger state one simulation step.
    ///
    /// For each enabled trigger volume tests each body for overlap and compares
    /// against previous occupancy:
    ///
    /// - [`TriggerEvent::Enter`] — body is newly inside the volume.
    /// - [`TriggerEvent::Exit`]  — body has left the volume (or disappeared).
    /// - [`TriggerEvent::Stay`]  — persistent overlap (only when
    ///   [`set_emit_stay`](TriggerWorld::set_emit_stay) is `true`).
    pub fn update(&mut self, step: u64, bodies: &[BodyEntry]) -> Vec<TriggerEvent> {
        let mut events = Vec::new();

        // ── Step 1: handle bodies absent from this step ──────────────────────
        let present: HashSet<usize> = bodies.iter().map(|b| b.index).collect();

        let absent_exits: Vec<(usize, u32)> = self
            .occupied
            .iter()
            .filter(|&(&bi, _)| !present.contains(&bi))
            .flat_map(|(&bi, ids)| ids.iter().map(move |&tid| (bi, tid)))
            .collect();

        for (body_index, trigger_id) in absent_exits {
            events.push(TriggerEvent::Exit {
                trigger_id,
                body_index,
                step,
            });
            if let Some(occ) = self.occupied.get_mut(&body_index) {
                occ.retain(|&t| t != trigger_id);
            }
        }
        self.occupied.retain(|_, ids| !ids.is_empty());

        // ── Step 2: compute overlaps (read-only on self.volumes) ─────────────
        let tests: Vec<OverlapTest> = self
            .volumes
            .iter()
            .filter(|v| v.enabled)
            .flat_map(|vol| {
                let vid = vol.id;
                bodies.iter().map(move |body| OverlapTest {
                    body_index: body.index,
                    trigger_id: vid,
                    now_inside: vol
                        .shape
                        .contains_point_with_radius(body.position, body.radius),
                })
            })
            .collect();

        // ── Step 3: apply results (write to self.occupied) ───────────────────
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
                    events.push(TriggerEvent::Enter {
                        trigger_id: test.trigger_id,
                        body_index: test.body_index,
                        step,
                    });
                }
                (true, false) => {
                    if let Some(occ) = self.occupied.get_mut(&test.body_index) {
                        occ.retain(|&t| t != test.trigger_id);
                    }
                    events.push(TriggerEvent::Exit {
                        trigger_id: test.trigger_id,
                        body_index: test.body_index,
                        step,
                    });
                }
                (true, true) => {
                    if self.emit_stay {
                        events.push(TriggerEvent::Stay {
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
}

// ============================================================================
// Display
// ============================================================================

impl std::fmt::Display for TriggerWorld {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TriggerWorld {{ volumes: {}, occupancies: {} }}",
            self.volumes.len(),
            self.active_occupancies(),
        )
    }
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn body(index: usize, pos: [f64; 3]) -> BodyEntry {
        BodyEntry {
            index,
            position: pos,
            radius: 0.5,
            is_sleeping: false,
        }
    }

    #[test]
    fn sphere_enter_exit() {
        let mut world = TriggerWorld::new();
        let id = world.add_sphere([0.0, 0.0, 0.0], 5.0, vec!["zone".into()]);

        let inside = vec![body(0, [0.0, 0.0, 0.0])];
        let ev = world.update(1, &inside);
        assert_eq!(ev.len(), 1);
        assert!(ev[0].is_enter());
        assert_eq!(ev[0].body_index(), 0);

        // Step 2: still inside, no stay emitted by default
        let ev = world.update(2, &inside);
        assert!(ev.is_empty());

        // Step 3: body moves outside
        let outside = vec![body(0, [20.0, 0.0, 0.0])];
        let ev = world.update(3, &outside);
        assert_eq!(ev.len(), 1);
        assert!(ev[0].is_exit());

        assert_eq!(world.bodies_in(id).len(), 0);
        assert_eq!(world.active_occupancies(), 0);
    }

    #[test]
    fn aabb_trigger() {
        let mut world = TriggerWorld::new();
        world.add_aabb([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], vec![]);

        let ev = world.update(1, &[body(0, [0.0, 0.0, 0.0])]);
        assert!(ev.iter().any(|e| e.is_enter()));
        let ev = world.update(2, &[body(0, [100.0, 0.0, 0.0])]);
        assert!(ev.iter().any(|e| e.is_exit()));
    }

    #[test]
    fn capsule_trigger() {
        let mut world = TriggerWorld::new();
        world.add_capsule([0.0, 0.0, 0.0], [0.0, 10.0, 0.0], 2.0, vec![]);

        let ev = world.update(1, &[body(0, [0.0, 5.0, 0.0])]);
        assert!(ev[0].is_enter());
    }

    #[test]
    fn stay_events_when_enabled() {
        let mut world = TriggerWorld::new();
        world.set_emit_stay(true);
        world.add_sphere([0.0, 0.0, 0.0], 5.0, vec![]);
        let bodies = vec![body(0, [0.0, 0.0, 0.0])];
        world.update(1, &bodies); // Enter
        let ev = world.update(2, &bodies); // Stay
        assert!(ev.iter().any(|e| e.is_stay()));
    }

    #[test]
    fn absent_body_emits_exit() {
        let mut world = TriggerWorld::new();
        world.add_sphere([0.0, 0.0, 0.0], 5.0, vec![]);
        world.update(1, &[body(0, [0.0, 0.0, 0.0])]);
        // Body disappears (not in next frame)
        let ev = world.update(2, &[]);
        assert!(ev.iter().any(|e| e.is_exit()));
    }

    #[test]
    fn remove_volume_clears_occupancy() {
        let mut world = TriggerWorld::new();
        let id = world.add_sphere([0.0, 0.0, 0.0], 5.0, vec![]);
        world.update(1, &[body(0, [0.0, 0.0, 0.0])]);
        assert_eq!(world.active_occupancies(), 1);
        world.remove_volume(id);
        assert_eq!(world.active_occupancies(), 0);
    }
}
