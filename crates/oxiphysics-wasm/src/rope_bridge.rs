// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the rope / chain simulation.
//!
//! Exposes a JSON-oriented surface around Verlet-integrated rope links
//! suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, flatten_vec3s, to_js_value, unflatten_vec3s};

/// Default [`WasmRope::damping`]: 2 % of the carried per-step displacement is
/// removed each Verlet step (a light, air-drag-like loss).
const DEFAULT_DAMPING: f64 = 0.02;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// WasmRopeLink
// ---------------------------------------------------------------------------

/// A single particle in the rope chain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmRopeLink {
    /// Current world-space position.
    pub position: [f64; 3],
    /// Previous world-space position (Verlet).
    pub prev_position: [f64; 3],
    /// Link mass (kg); must be > 0.
    pub mass: f64,
    /// Visual / collision radius (metres).
    pub radius: f64,
    /// Whether this link is pinned in place.
    pub pinned: bool,
}

// ---------------------------------------------------------------------------
// WasmRopeError
// ---------------------------------------------------------------------------

/// Error variants returned by [`WasmRope::new`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmRopeError {
    TooFewLinks,
    ZeroMassLink(usize),
    ZeroSegmentLength,
}

impl std::fmt::Display for WasmRopeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WasmRopeError::TooFewLinks => write!(f, "rope must have at least 2 links"),
            WasmRopeError::ZeroMassLink(i) => write!(f, "link {i} has zero or negative mass"),
            WasmRopeError::ZeroSegmentLength => write!(f, "segment rest-length must be > 0"),
        }
    }
}

// ---------------------------------------------------------------------------
// WasmRope
// ---------------------------------------------------------------------------

/// WASM wrapper for a Verlet-integrated rope / chain.
///
/// Links are connected by distance constraints enforced using Gauss-Seidel
/// projection. Head and tail can optionally be pinned to a world-space point.
///
/// Most fields are private when accessed from JavaScript — use the
/// `#[wasm_bindgen]` accessors for `Vec<f64>` views over arrays and
/// `Option<[f64; 3]>` anchors.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmRope {
    /// Rope particles — private; access from JS via the flat-array accessors.
    pub(crate) links: Vec<WasmRopeLink>,
    /// Rest length between adjacent links (metres).
    pub rest_length: f64,
    /// Gravity vector — private; access from JS via `get_gravity` / `set_gravity`.
    gravity: [f64; 3],
    /// Velocity damping per Verlet step, in `[0, 1]`: the fraction of the
    /// carried per-step displacement (`position - prev_position`) that is
    /// removed each step.  `0` leaves the rope undamped and `1` stops it
    /// after a single step.  Values outside `[0, 1]` are clamped when
    /// stepping, so damping can never amplify the motion.
    pub damping: f64,
    /// Constraint solve iterations — private; access via `iterations()` (`u32`).
    iterations: usize,
    /// Bending stiffness coefficient.
    pub bending_stiffness: f64,
    /// World-space anchor for the head link — private; use the JS accessors.
    anchor_head: Option<[f64; 3]>,
    /// World-space anchor for the tail link — private; use the JS accessors.
    anchor_tail: Option<[f64; 3]>,
}

impl WasmRope {
    /// Construct a rope of `link_count` links hanging downward from `origin`.
    ///
    /// `segment_length` is the rest length between adjacent links (metres).
    /// `link_mass` is the mass of every link (kg).
    pub fn new(
        origin: [f64; 3],
        link_count: usize,
        segment_length: f64,
        link_mass: f64,
    ) -> Result<Self, WasmRopeError> {
        if link_count < 2 {
            return Err(WasmRopeError::TooFewLinks);
        }
        if link_mass <= 0.0 {
            return Err(WasmRopeError::ZeroMassLink(0));
        }
        if segment_length <= 0.0 {
            return Err(WasmRopeError::ZeroSegmentLength);
        }
        let links = (0..link_count)
            .map(|i| {
                let y_offset = -(i as f64) * segment_length;
                let pos = [origin[0], origin[1] + y_offset, origin[2]];
                WasmRopeLink {
                    position: pos,
                    prev_position: pos,
                    mass: link_mass,
                    radius: 0.05,
                    pinned: false,
                }
            })
            .collect();
        Ok(WasmRope {
            links,
            rest_length: segment_length,
            gravity: [0.0, -9.81, 0.0],
            damping: DEFAULT_DAMPING,
            iterations: 8,
            bending_stiffness: 0.0,
            anchor_head: None,
            anchor_tail: None,
        })
    }

    /// Pin the head (index 0) to a fixed world-space position.
    pub fn pin_head(&mut self, pos: [f64; 3]) {
        self.anchor_head = Some(pos);
        if let Some(link) = self.links.first_mut() {
            link.position = pos;
            link.prev_position = pos;
            link.pinned = true;
        }
    }

    /// Pin the tail (last link) to a fixed world-space position.
    pub fn pin_tail(&mut self, pos: [f64; 3]) {
        self.anchor_tail = Some(pos);
        if let Some(link) = self.links.last_mut() {
            link.position = pos;
            link.prev_position = pos;
            link.pinned = true;
        }
    }

    /// Simulate one time step of `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let n = self.links.len();
        if n == 0 {
            return;
        }
        let dt2 = dt * dt;

        // Step 1: damped Verlet integration (skip pinned links).
        //
        // The implicit velocity of a Verlet particle is the displacement
        // `position - prev_position` of the previous step.  Damping removes
        // the fraction `damping` of that carried displacement, so the particle
        // advances by `(1 - damping) * (position - prev_position)` plus the
        // acceleration term `gravity * dt^2`.  With `damping == 0` this is the
        // classic undamped Verlet update and with `damping == 1` the particle
        // keeps no momentum from the previous step.  The motion is therefore
        // never amplified by damping.
        let retention = 1.0 - self.damping.clamp(0.0, 1.0);
        let accel_step = scale(self.gravity, dt2);
        for link in self.links.iter_mut() {
            if link.pinned {
                continue;
            }
            let carried = scale(sub(link.position, link.prev_position), retention);
            link.prev_position = link.position;
            link.position = add(add(link.position, carried), accel_step);
        }

        // Step 2: Gauss-Seidel distance constraint
        for _ in 0..self.iterations {
            for i in 0..n - 1 {
                let pi_pinned = self.links[i].pinned;
                let pj_pinned = self.links[i + 1].pinned;

                let pos_i = self.links[i].position;
                let pos_j = self.links[i + 1].position;
                let mi = self.links[i].mass;
                let mj = self.links[i + 1].mass;

                let delta = sub(pos_j, pos_i);
                let dist = len(delta);
                if dist < 1.0e-12 {
                    continue;
                }

                let correction = (dist - self.rest_length) / dist;
                let wi = if pi_pinned { 0.0 } else { 1.0 / mi };
                let wj = if pj_pinned { 0.0 } else { 1.0 / mj };
                let wsum = wi + wj;
                if wsum < 1.0e-30 {
                    continue;
                }

                if !pi_pinned {
                    self.links[i].position = add(pos_i, scale(delta, correction * (wi / wsum)));
                }
                if !pj_pinned {
                    self.links[i + 1].position = sub(pos_j, scale(delta, correction * (wj / wsum)));
                }
            }
        }

        // Step 3: Re-apply anchors
        if let Some(head_pos) = self.anchor_head
            && let Some(link) = self.links.first_mut()
        {
            link.position = head_pos;
        }
        if let Some(tail_pos) = self.anchor_tail
            && let Some(link) = self.links.last_mut()
        {
            link.position = tail_pos;
        }
    }

    /// Number of links in this rope.
    pub fn link_count(&self) -> usize {
        self.links.len()
    }

    /// Sum of all current segment lengths.
    pub fn total_length(&self) -> f64 {
        self.links
            .windows(2)
            .map(|w| len(sub(w[1].position, w[0].position)))
            .sum()
    }

    /// Serialise the rope state as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// wasm-bindgen JavaScript API
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmRope {
    /// Construct a hanging rope from `(ox, oy, oz)`.
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string if `link_count < 2`, `link_mass <= 0`,
    /// or `segment_length <= 0`.
    #[wasm_bindgen(constructor)]
    pub fn new_js(
        ox: f64,
        oy: f64,
        oz: f64,
        link_count: u32,
        segment_length: f64,
        link_mass: f64,
    ) -> std::result::Result<WasmRope, JsValue> {
        WasmRope::new([ox, oy, oz], link_count as usize, segment_length, link_mass)
            .map_err(err_to_jsvalue)
    }

    /// Pin the head (index 0) to a fixed world-space position.
    #[wasm_bindgen(js_name = "pin_head")]
    pub fn pin_head_js(&mut self, x: f64, y: f64, z: f64) {
        self.pin_head([x, y, z]);
    }

    /// Pin the tail (last link) to a fixed world-space position.
    #[wasm_bindgen(js_name = "pin_tail")]
    pub fn pin_tail_js(&mut self, x: f64, y: f64, z: f64) {
        self.pin_tail([x, y, z]);
    }

    /// Clear the head anchor (allowing the head link to swing freely).
    #[wasm_bindgen(js_name = "unpin_head")]
    pub fn unpin_head_js(&mut self) {
        self.anchor_head = None;
        if let Some(link) = self.links.first_mut() {
            link.pinned = false;
        }
    }

    /// Clear the tail anchor.
    #[wasm_bindgen(js_name = "unpin_tail")]
    pub fn unpin_tail_js(&mut self) {
        self.anchor_tail = None;
        if let Some(link) = self.links.last_mut() {
            link.pinned = false;
        }
    }

    /// Advance the simulation by `dt` seconds.
    #[wasm_bindgen(js_name = "step")]
    pub fn step_js(&mut self, dt: f64) {
        self.step(dt);
    }

    /// Number of links in this rope.
    #[wasm_bindgen(js_name = "link_count")]
    pub fn link_count_js(&self) -> u32 {
        self.link_count() as u32
    }

    /// Sum of all current segment lengths.
    #[wasm_bindgen(js_name = "total_length")]
    pub fn total_length_js(&self) -> f64 {
        self.total_length()
    }

    /// Number of constraint iterations per step.
    #[wasm_bindgen(js_name = "iterations")]
    pub fn iterations_js(&self) -> u32 {
        self.iterations as u32
    }

    /// Set the number of constraint iterations per step (clamped to >= 1).
    #[wasm_bindgen(js_name = "set_iterations")]
    pub fn set_iterations_js(&mut self, iterations: u32) {
        self.iterations = (iterations as usize).max(1);
    }

    /// Gravity vector as `[gx, gy, gz]` (`Float64Array`).
    #[wasm_bindgen(js_name = "get_gravity")]
    pub fn get_gravity_js(&self) -> Vec<f64> {
        self.gravity.to_vec()
    }

    /// Set the gravity vector.
    #[wasm_bindgen(js_name = "set_gravity")]
    pub fn set_gravity_js(&mut self, gx: f64, gy: f64, gz: f64) {
        self.gravity = [gx, gy, gz];
    }

    /// Head anchor as `[x, y, z]` or an empty array when unanchored.
    #[wasm_bindgen(js_name = "get_anchor_head")]
    pub fn get_anchor_head_js(&self) -> Vec<f64> {
        self.anchor_head.map(|a| a.to_vec()).unwrap_or_default()
    }

    /// Tail anchor as `[x, y, z]` or an empty array when unanchored.
    #[wasm_bindgen(js_name = "get_anchor_tail")]
    pub fn get_anchor_tail_js(&self) -> Vec<f64> {
        self.anchor_tail.map(|a| a.to_vec()).unwrap_or_default()
    }

    /// Flat link positions `[x0, y0, z0, x1, y1, z1, ...]` as `Float64Array`.
    #[wasm_bindgen(js_name = "get_positions_flat")]
    pub fn get_positions_flat_js(&self) -> Vec<f64> {
        let positions: Vec<[f64; 3]> = self.links.iter().map(|l| l.position).collect();
        flatten_vec3s(&positions)
    }

    /// Flat previous-step positions used by Verlet integration.
    #[wasm_bindgen(js_name = "get_prev_positions_flat")]
    pub fn get_prev_positions_flat_js(&self) -> Vec<f64> {
        let positions: Vec<[f64; 3]> = self.links.iter().map(|l| l.prev_position).collect();
        flatten_vec3s(&positions)
    }

    /// Per-link masses as `Float64Array`.
    #[wasm_bindgen(js_name = "get_masses")]
    pub fn get_masses_js(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.mass).collect()
    }

    /// Per-link radii as `Float64Array`.
    #[wasm_bindgen(js_name = "get_radii")]
    pub fn get_radii_js(&self) -> Vec<f64> {
        self.links.iter().map(|l| l.radius).collect()
    }

    /// Per-link `pinned` flags packed into a `Uint8Array` (0 / 1).
    #[wasm_bindgen(js_name = "get_pinned_flags")]
    pub fn get_pinned_flags_js(&self) -> Vec<u8> {
        self.links.iter().map(|l| u8::from(l.pinned)).collect()
    }

    /// Replace all link positions from a flat `[x0, y0, z0, ...]` buffer.
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error if `flat.len()` is not divisible by 3 or
    /// does not match the current link count.
    #[wasm_bindgen(js_name = "set_positions_flat")]
    pub fn set_positions_flat_js(&mut self, flat: &[f64]) -> std::result::Result<(), JsValue> {
        let positions = unflatten_vec3s(flat).map_err(err_to_jsvalue)?;
        if positions.len() != self.links.len() {
            return Err(JsValue::from_str(&format!(
                "expected {} positions, got {}",
                self.links.len(),
                positions.len()
            )));
        }
        for (link, pos) in self.links.iter_mut().zip(positions.iter()) {
            link.position = *pos;
        }
        Ok(())
    }

    /// Serialise the rope state as a `JsValue` (plain JS object).
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string if serialisation fails.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> std::result::Result<JsValue, JsValue> {
        to_js_value(self)
    }

    /// Serialise the rope state as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rope_bridge_instantiation() {
        let rope = WasmRope::new([0.0, 5.0, 0.0], 5, 1.0, 1.0).expect("valid rope");
        let json = rope.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_rope_bridge_too_few_links() {
        assert_eq!(
            WasmRope::new([0.0; 3], 1, 1.0, 1.0),
            Err(WasmRopeError::TooFewLinks)
        );
    }

    #[test]
    fn test_rope_bridge_zero_mass() {
        assert_eq!(
            WasmRope::new([0.0; 3], 3, 1.0, 0.0),
            Err(WasmRopeError::ZeroMassLink(0))
        );
    }

    #[test]
    fn test_rope_bridge_step_does_not_panic() {
        let mut rope = WasmRope::new([0.0, 5.0, 0.0], 4, 1.0, 1.0).expect("valid rope");
        for _ in 0..60 {
            rope.step(1.0 / 60.0);
        }
        // Head should have fallen (no anchors)
        assert!(rope.links[0].position[1] < 5.0);
    }

    #[test]
    fn test_rope_bridge_pinned_head_stays() {
        let mut rope = WasmRope::new([0.0, 5.0, 0.0], 4, 1.0, 1.0).expect("valid rope");
        rope.pin_head([0.0, 5.0, 0.0]);
        for _ in 0..200 {
            rope.step(1.0 / 60.0);
        }
        let head = rope.links[0].position;
        assert!(
            (head[1] - 5.0).abs() < 0.01,
            "pinned head y drifted to {}",
            head[1]
        );
    }

    #[test]
    fn test_rope_bridge_serde_roundtrip() {
        let rope = WasmRope::new([1.0, 2.0, 3.0], 4, 0.5, 2.0).expect("valid rope");
        let json = rope.to_json();
        let restored: WasmRope = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(restored.links.len(), rope.links.len());
        assert!((restored.rest_length - rope.rest_length).abs() < 1e-12);
    }

    // ---- damping regression (F5729) --------------------------------------

    const DT: f64 = 1.0 / 60.0;

    /// Weightless, unpinned rope translating along +x at `speed` m/s: every
    /// link carries the same displacement, so the distance constraints stay
    /// satisfied and the rope's motion is governed by damping alone.
    fn drifting_rope(damping: f64, speed: f64) -> WasmRope {
        let mut rope = WasmRope::new([0.0, 5.0, 0.0], 6, 0.25, 1.0).expect("valid rope");
        rope.set_gravity_js(0.0, 0.0, 0.0);
        rope.damping = damping;
        for link in rope.links.iter_mut() {
            link.prev_position[0] = link.position[0] - speed * DT;
        }
        rope
    }

    /// Speed of the head link implied by its Verlet displacement (m/s).
    fn head_speed(rope: &WasmRope) -> f64 {
        let link = &rope.links[0];
        len(sub(link.position, link.prev_position)) / DT
    }

    #[test]
    fn test_rope_default_damping_is_a_small_fraction() {
        let rope = WasmRope::new([0.0; 3], 3, 1.0, 1.0).expect("valid rope");
        assert!((rope.damping - 0.02).abs() < 1e-15, "got {}", rope.damping);
    }

    #[test]
    fn test_rope_zero_damping_preserves_speed() {
        let mut rope = drifting_rope(0.0, 1.0);
        assert!((head_speed(&rope) - 1.0).abs() < 1e-12);
        for _ in 0..600 {
            rope.step(DT);
        }
        assert!(
            (head_speed(&rope) - 1.0).abs() < 1e-9,
            "undamped rope changed speed to {}",
            head_speed(&rope)
        );
        // 600 steps at 1 m/s and 60 Hz travel exactly 10 m along +x.
        assert!((rope.links[0].position[0] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_rope_damping_never_amplifies_speed() {
        // Regression: the old update multiplied the carried velocity by
        // (1 + damping) per step, i.e. ~1.98x per step at the default.
        let mut rope = drifting_rope(0.02, 1.0);
        let initial = head_speed(&rope);
        let mut previous = initial;
        for step in 0..600 {
            rope.step(DT);
            let speed = head_speed(&rope);
            assert!(
                speed <= previous + 1e-12,
                "speed grew at step {step}: {previous} -> {speed}"
            );
            assert!(speed <= initial + 1e-12, "speed exceeded its initial value");
            previous = speed;
        }
        assert!(previous < 1e-4, "rope barely damped: speed {previous}");
    }

    #[test]
    fn test_rope_damping_decays_geometrically() {
        let damping = 0.1;
        let mut rope = drifting_rope(damping, 2.0);
        for _ in 0..10 {
            rope.step(DT);
        }
        let expected = 2.0 * (1.0 - damping).powi(10);
        assert!(
            (head_speed(&rope) - expected).abs() < 1e-9,
            "speed {} != expected {expected}",
            head_speed(&rope)
        );
    }

    #[test]
    fn test_rope_full_damping_stops_after_one_step() {
        let mut rope = drifting_rope(1.0, 3.0);
        rope.step(DT);
        assert!(head_speed(&rope) < 1e-12, "speed {}", head_speed(&rope));
    }

    #[test]
    fn test_rope_out_of_range_damping_is_clamped() {
        // A negative damping used to be an amplifier; it is clamped to 0.
        let mut negative = drifting_rope(-0.5, 1.0);
        for _ in 0..30 {
            negative.step(DT);
        }
        assert!((head_speed(&negative) - 1.0).abs() < 1e-9);
        // A damping above 1 is clamped to 1 (full damping): all carried motion
        // is removed, the rope stays where it was instead of reversing.
        let mut excessive = drifting_rope(2.0, 1.0);
        excessive.step(DT);
        assert!(head_speed(&excessive) < 1e-12);
        let drift = excessive.links[0].position[0];
        assert!(drift.abs() < 1e-12, "rope moved or reversed: x = {drift}");
    }

    #[test]
    fn test_rope_free_fall_reaches_damped_terminal_speed() {
        // v <- (1 - d) v + g dt^2 has the fixed point g dt^2 / d (per step).
        let mut rope = WasmRope::new([0.0, 50.0, 0.0], 5, 0.5, 1.0).expect("valid rope");
        let damping = rope.damping;
        for _ in 0..1200 {
            rope.step(DT);
        }
        let link = &rope.links[2];
        let per_step_drop = link.prev_position[1] - link.position[1];
        let expected = 9.81 * DT * DT / damping;
        assert!(
            (per_step_drop - expected).abs() < 0.01 * expected,
            "per-step drop {per_step_drop} != terminal {expected}"
        );
    }
}
