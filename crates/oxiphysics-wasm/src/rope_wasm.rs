// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `#[wasm_bindgen]` wrapper for [`WasmRope`] — exposes the Verlet rope
//! simulation to JavaScript with a flat-array API surface.

use wasm_bindgen::prelude::*;

use crate::rope_bridge::WasmRope;

/// JavaScript-accessible wrapper around the Verlet-integrated rope simulation.
///
/// ## JavaScript example
///
/// ```js
/// const rope = WasmRopeJs.new(0.0, 5.0, 0.0, 16, 0.3, 1.0);
/// rope.pin_head(0.0, 5.0, 0.0);
/// for (let i = 0; i < 60; i++) rope.step(1 / 60);
/// const positions = rope.all_positions_flat(); // Float64Array [x0,y0,z0, ...]
/// ```
#[wasm_bindgen]
pub struct WasmRopeJs {
    inner: WasmRope,
}

#[wasm_bindgen]
impl WasmRopeJs {
    /// Create a rope of `link_count` links hanging from `(origin_x, origin_y, origin_z)`.
    ///
    /// Returns `null` if parameters are invalid (fewer than 2 links, zero mass, zero length).
    pub fn new(
        origin_x: f64,
        origin_y: f64,
        origin_z: f64,
        link_count: u32,
        segment_length: f64,
        link_mass: f64,
    ) -> Option<WasmRopeJs> {
        WasmRope::new(
            [origin_x, origin_y, origin_z],
            link_count as usize,
            segment_length,
            link_mass,
        )
        .ok()
        .map(|inner| WasmRopeJs { inner })
    }

    /// Pin the head link (index 0) to a fixed world-space position.
    pub fn pin_head(&mut self, x: f64, y: f64, z: f64) {
        self.inner.pin_head([x, y, z]);
    }

    /// Pin the tail link (last) to a fixed world-space position.
    pub fn pin_tail(&mut self, x: f64, y: f64, z: f64) {
        self.inner.pin_tail([x, y, z]);
    }

    /// Apply an impulse velocity to link `index`.
    ///
    /// Adds `(ivx, ivy, ivz)` to the link's velocity by shifting its previous
    /// position backward by the impulse amount times `dt`.
    pub fn apply_impulse(&mut self, index: u32, ivx: f64, ivy: f64, ivz: f64) {
        let idx = index as usize;
        if let Some(link) = self.inner.links.get_mut(idx)
            && !link.pinned
        {
            // Verlet: modifying prev_position changes effective velocity
            link.prev_position[0] -= ivx;
            link.prev_position[1] -= ivy;
            link.prev_position[2] -= ivz;
        }
    }

    /// Advance the rope simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.inner.step(dt);
    }

    /// Number of links in the rope.
    pub fn link_count(&self) -> u32 {
        self.inner.link_count() as u32
    }

    /// Total arc-length of all segments at the current simulation state.
    pub fn total_length(&self) -> f64 {
        self.inner.total_length()
    }

    /// All link positions as a flat `Float64Array` `[x0,y0,z0, x1,y1,z1, ...]`.
    pub fn all_positions_flat(&self) -> Vec<f64> {
        self.inner
            .links
            .iter()
            .flat_map(|l| [l.position[0], l.position[1], l.position[2]])
            .collect()
    }

    /// Position of link `index` as `[x, y, z]`.
    pub fn get_position(&self, index: u32) -> Vec<f64> {
        self.inner
            .links
            .get(index as usize)
            .map(|l| l.position.to_vec())
            .unwrap_or_else(|| vec![0.0; 3])
    }

    /// Set the Verlet damping coefficient (0 = no damping, 1 = fully damped).
    ///
    /// Typical range: 0.90–0.99. Values outside `[0, 1]` are clamped.
    pub fn set_damping(&mut self, damping: f64) {
        self.inner.damping = damping.clamp(0.0, 1.0);
    }

    /// Current Verlet damping coefficient.
    pub fn get_damping(&self) -> f64 {
        self.inner.damping
    }

    /// Serialise the rope state as a JSON string.
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
}
