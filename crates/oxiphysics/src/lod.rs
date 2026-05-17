// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Level-of-Detail (LOD) simulation tier management.
//!
//! Large physics scenes cannot afford full-fidelity simulation for every body.
//! This module assigns bodies to one of four `LodTier`s based on their
//! distance from a focus point (camera, player, area of interest).
//!
//! ## Tiers
//!
//! | Tier | Description | Default substeps |
//! |------|-------------|-----------------|
//! | `Full` | All features, maximum substeps | 4 |
//! | `Reduced` | Simplified collision, fewer substeps | 2 |
//! | `Minimal` | Position extrapolation only | 1 |
//! | `Frozen` | No simulation — body is static | 0 |
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::lod::{LodSystem, LodConfig, LodTier};
//!
//! let config = LodConfig {
//!     full_radius:    50.0,
//!     reduced_radius: 150.0,
//!     minimal_radius: 300.0,
//!     substeps_full:    4,
//!     substeps_reduced: 2,
//!     substeps_minimal: 1,
//! };
//! let mut lod = LodSystem::new(config, [0.0; 3]);
//!
//! let near_body = lod.add_body([10.0, 0.0, 0.0], 0.0);
//! let far_body  = lod.add_body([500.0, 0.0, 0.0], 0.0);
//!
//! let transitions = lod.update();
//! assert_eq!(lod.tier(near_body), Some(LodTier::Full));
//! assert_eq!(lod.tier(far_body), Some(LodTier::Frozen));
//! assert!(!transitions.is_empty());
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// LodTier
// ---------------------------------------------------------------------------

/// Simulation fidelity tier for a body.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LodTier {
    /// Full physics simulation: all features, maximum substeps.
    Full,
    /// Simplified simulation: broadphase only or reduced contact iterations,
    /// fewer substeps.
    Reduced,
    /// Minimal simulation: linear velocity extrapolation only.
    Minimal,
    /// Frozen: no simulation, body is treated as static.
    Frozen,
}

// ---------------------------------------------------------------------------
// LodConfig
// ---------------------------------------------------------------------------

/// Distance thresholds and substep counts for LOD tier assignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LodConfig {
    /// Distance (m) within which bodies receive [`LodTier::Full`] simulation.
    pub full_radius: f64,
    /// Distance (m) within which bodies receive [`LodTier::Reduced`] simulation.
    ///
    /// Must be `>= full_radius`.
    pub reduced_radius: f64,
    /// Distance (m) within which bodies receive [`LodTier::Minimal`] simulation.
    ///
    /// Bodies beyond this distance receive [`LodTier::Frozen`].
    /// Must be `>= reduced_radius`.
    pub minimal_radius: f64,
    /// Substep count for [`LodTier::Full`] bodies.
    pub substeps_full: u32,
    /// Substep count for [`LodTier::Reduced`] bodies.
    pub substeps_reduced: u32,
    /// Substep count for [`LodTier::Minimal`] bodies.
    pub substeps_minimal: u32,
}

impl Default for LodConfig {
    fn default() -> Self {
        Self {
            full_radius: 50.0,
            reduced_radius: 150.0,
            minimal_radius: 300.0,
            substeps_full: 4,
            substeps_reduced: 2,
            substeps_minimal: 1,
        }
    }
}

// ---------------------------------------------------------------------------
// LodBody
// ---------------------------------------------------------------------------

/// A body tracked by the LOD system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LodBody {
    /// Unique body identifier (assigned by [`LodSystem::add_body`]).
    pub id: u32,
    /// Current world-space position of the body (m).
    pub position: [f64; 3],
    /// Current LOD tier (updated by [`LodSystem::update`]).
    pub tier: LodTier,
    /// Priority bias in `[0.0, 1.0]`.
    ///
    /// Higher values shift the effective distance toward zero, biasing the body
    /// toward finer LOD tiers.  A value of `1.0` reduces effective distance by
    /// up to 30%.
    pub priority: f32,
}

// ---------------------------------------------------------------------------
// LodUpdate
// ---------------------------------------------------------------------------

/// A tier transition emitted by [`LodSystem::update`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LodUpdate {
    /// ID of the body that changed tier.
    pub id: u32,
    /// Tier before the update.
    pub old_tier: LodTier,
    /// Tier after the update.
    pub new_tier: LodTier,
}

// ---------------------------------------------------------------------------
// LodSystem
// ---------------------------------------------------------------------------

/// Level-of-detail system that assigns bodies to simulation tiers based on
/// their distance from a focus point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LodSystem {
    /// LOD configuration (distance thresholds and substep counts).
    pub config: LodConfig,
    bodies: Vec<LodBody>,
    next_id: u32,
    /// World-space focus point (camera position, player, area of interest).
    pub focus: [f64; 3],
}

impl LodSystem {
    /// Create a new LOD system.
    pub fn new(config: LodConfig, focus: [f64; 3]) -> Self {
        Self {
            config,
            bodies: Vec::new(),
            next_id: 0,
            focus,
        }
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Register a body and return its ID.
    ///
    /// `priority` is clamped to `[0.0, 1.0]`.  The initial tier is set to
    /// [`LodTier::Full`]; call [`update`][Self::update] to compute the correct
    /// initial tier.
    pub fn add_body(&mut self, position: [f64; 3], priority: f32) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.bodies.push(LodBody {
            id,
            position,
            tier: LodTier::Full,
            priority: priority.clamp(0.0, 1.0),
        });
        id
    }

    /// Unregister a body by ID.  No-op if not found.
    pub fn remove_body(&mut self, id: u32) {
        self.bodies.retain(|b| b.id != id);
    }

    /// Immutable reference to a body entry.
    pub fn body(&self, id: u32) -> Option<&LodBody> {
        self.bodies.iter().find(|b| b.id == id)
    }

    /// Update the world-space position of a body.  No-op if not found.
    pub fn update_position(&mut self, id: u32, position: [f64; 3]) {
        if let Some(b) = self.bodies.iter_mut().find(|b| b.id == id) {
            b.position = position;
        }
    }

    /// Update the priority bias of a body.
    ///
    /// `priority` is clamped to `[0.0, 1.0]`.
    pub fn set_priority(&mut self, id: u32, priority: f32) {
        if let Some(b) = self.bodies.iter_mut().find(|b| b.id == id) {
            b.priority = priority.clamp(0.0, 1.0);
        }
    }

    // -----------------------------------------------------------------------
    // Focus
    // -----------------------------------------------------------------------

    /// Move the focus point.
    ///
    /// Call [`update`][Self::update] afterwards to propagate tier changes.
    pub fn set_focus(&mut self, focus: [f64; 3]) {
        self.focus = focus;
    }

    // -----------------------------------------------------------------------
    // LOD tier computation
    // -----------------------------------------------------------------------

    /// Compute the LOD tier for a body at `dist` from the focus point given a
    /// priority bias.
    fn compute_tier(&self, dist: f64, priority: f32) -> LodTier {
        // High-priority bodies see a reduced effective distance (≤ 30% closer).
        let effective = dist * (1.0 - (priority as f64) * 0.3);
        let cfg = &self.config;
        if effective <= cfg.full_radius {
            LodTier::Full
        } else if effective <= cfg.reduced_radius {
            LodTier::Reduced
        } else if effective <= cfg.minimal_radius {
            LodTier::Minimal
        } else {
            LodTier::Frozen
        }
    }

    /// Recompute LOD tiers for all bodies based on the current focus point.
    ///
    /// Returns [`LodUpdate`] entries for bodies whose tier changed.
    pub fn update(&mut self) -> Vec<LodUpdate> {
        let focus = self.focus;
        let mut transitions = Vec::new();

        for body in &mut self.bodies {
            let dx = body.position[0] - focus[0];
            let dy = body.position[1] - focus[1];
            let dz = body.position[2] - focus[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let new_tier = {
                let effective = dist * (1.0 - (body.priority as f64) * 0.3);
                let cfg = &self.config;
                if effective <= cfg.full_radius {
                    LodTier::Full
                } else if effective <= cfg.reduced_radius {
                    LodTier::Reduced
                } else if effective <= cfg.minimal_radius {
                    LodTier::Minimal
                } else {
                    LodTier::Frozen
                }
            };

            if body.tier != new_tier {
                transitions.push(LodUpdate {
                    id: body.id,
                    old_tier: body.tier,
                    new_tier,
                });
                body.tier = new_tier;
            }
        }

        transitions
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Current LOD tier for a body.  Returns `None` if `id` is not registered.
    pub fn tier(&self, id: u32) -> Option<LodTier> {
        self.bodies.iter().find(|b| b.id == id).map(|b| b.tier)
    }

    /// Substep count recommended for a body given its current tier.
    ///
    /// Returns `0` for [`LodTier::Frozen`] bodies and unregistered IDs.
    pub fn substeps_for(&self, id: u32) -> u32 {
        match self.tier(id) {
            Some(LodTier::Full) => self.config.substeps_full,
            Some(LodTier::Reduced) => self.config.substeps_reduced,
            Some(LodTier::Minimal) => self.config.substeps_minimal,
            Some(LodTier::Frozen) | None => 0,
        }
    }

    /// Collect IDs of all bodies in the given tier.
    pub fn bodies_in_tier(&self, tier: LodTier) -> Vec<u32> {
        self.bodies
            .iter()
            .filter(|b| b.tier == tier)
            .map(|b| b.id)
            .collect()
    }

    /// Number of registered bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Iterate over all body entries.
    pub fn bodies(&self) -> impl Iterator<Item = &LodBody> {
        self.bodies.iter()
    }

    /// Distance from the focus point to a body.  Returns `None` if not found.
    pub fn distance_to_focus(&self, id: u32) -> Option<f64> {
        self.bodies.iter().find(|b| b.id == id).map(|b| {
            let dx = b.position[0] - self.focus[0];
            let dy = b.position[1] - self.focus[1];
            let dz = b.position[2] - self.focus[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        })
    }

    /// Tier distribution summary: `(full, reduced, minimal, frozen)` counts.
    pub fn tier_counts(&self) -> (usize, usize, usize, usize) {
        let mut counts = (0usize, 0usize, 0usize, 0usize);
        for b in &self.bodies {
            match b.tier {
                LodTier::Full => counts.0 += 1,
                LodTier::Reduced => counts.1 += 1,
                LodTier::Minimal => counts.2 += 1,
                LodTier::Frozen => counts.3 += 1,
            }
        }
        counts
    }

    /// Return `true` if the given body would change tier if [`update`][Self::update]
    /// were called right now.
    pub fn would_transition(&self, id: u32) -> bool {
        if let Some(body) = self.bodies.iter().find(|b| b.id == id) {
            let dx = body.position[0] - self.focus[0];
            let dy = body.position[1] - self.focus[1];
            let dz = body.position[2] - self.focus[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            self.compute_tier(dist, body.priority) != body.tier
        } else {
            false
        }
    }
}
