// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rope / chain distance constraints with Verlet integration and Gauss-Seidel projection.
//!
//! ## Types
//!
//! - `RopeLink` — a single particle in the rope chain.
//! - `AnchorKind` — fixed world-space point or dynamic body handle.
//! - `Rope` — the full rope: links, rest length, stiffness, and anchors.
//! - `BodyTransformLookup` — closure type alias for body-handle → world position.
//! - `RopeError` — error variants returned by constructors.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::rope::{AnchorKind, Rope};
//!
//! // Rope of 6 links hanging from [0, 5, 0], segment rest-length 1 m.
//! let mut rope = Rope::new([0.0, 5.0, 0.0], 6, 1.0, 1.0).expect("valid rope");
//! rope.anchor_head = Some(AnchorKind::World([0.0, 5.0, 0.0]));
//!
//! let no_bodies = &(|_: u64| -> Option<[f64; 3]> { None });
//! let dt = 1.0 / 60.0;
//! for _ in 0..200 {
//!     rope.step(dt, no_bodies);
//! }
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Helper math (Vec3 = [f64; 3])
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
// RopeLink
// ---------------------------------------------------------------------------

/// A single particle in the rope chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RopeLink {
    /// Current world-space position.
    pub position: [f64; 3],
    /// Previous world-space position — used by the Verlet integrator.
    pub prev_position: [f64; 3],
    /// Link mass; must be > 0.
    pub mass: f64,
    /// Visual / collision radius (does not affect constraint solving).
    pub radius: f64,
}

// ---------------------------------------------------------------------------
// AnchorKind
// ---------------------------------------------------------------------------

/// Anchor point for the head or tail of the rope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchorKind {
    /// Fixed world-space position.
    World([f64; 3]),
    /// Body handle; caller provides the world-space transform via lookup.
    Body(u64),
}

// ---------------------------------------------------------------------------
// BodyTransformLookup
// ---------------------------------------------------------------------------

/// Minimal body-transform lookup: body-handle → world-space position.
///
/// Implement as a closure: `&|handle: u64| -> Option<[f64; 3]> { … }`.
pub type BodyTransformLookup<'a> = dyn Fn(u64) -> Option<[f64; 3]> + 'a;

// ---------------------------------------------------------------------------
// RopeError
// ---------------------------------------------------------------------------

/// Error type for rope operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RopeError {
    /// A link had mass ≤ 0; contains the offending index.
    ZeroMassLink(usize),
    /// The rope requires at least 2 links.
    TooFewLinks,
}

impl std::fmt::Display for RopeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RopeError::ZeroMassLink(i) => write!(f, "link {i} has non-positive mass"),
            RopeError::TooFewLinks => write!(f, "rope requires at least 2 links"),
        }
    }
}

impl std::error::Error for RopeError {}

// ---------------------------------------------------------------------------
// Rope
// ---------------------------------------------------------------------------

/// A position-based rope / chain of distance constraints.
///
/// Integration uses velocity Verlet; constraints are projected via two-pass
/// Gauss-Seidel.  Optional bending resistance is applied as a midpoint
/// correction on each triplet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rope {
    /// Ordered chain of rope particles.
    pub links: Vec<RopeLink>,
    /// Rest length per segment (uniform).
    pub rest_length: f64,
    /// Bending stiffness in `[0, 1]`; 0 = no bending resistance.
    pub bending_stiffness: f64,
    /// Velocity damping coefficient in `[0, 1)`.
    pub damping: f64,
    /// Per-step gravitational acceleration.
    pub gravity: [f64; 3],
    /// If `Some`, link 0 is pinned to this anchor.
    pub anchor_head: Option<AnchorKind>,
    /// If `Some`, link `N-1` is pinned to this anchor.
    pub anchor_tail: Option<AnchorKind>,
}

impl Rope {
    /// Create a straight rope hanging along −Y from `root`.
    ///
    /// Returns `Err(RopeError::TooFewLinks)` if `num_links < 2`.
    /// Returns `Err(RopeError::ZeroMassLink(0))` if `mass_per_link ≤ 0`.
    pub fn new(
        root: [f64; 3],
        num_links: usize,
        segment_length: f64,
        mass_per_link: f64,
    ) -> Result<Self, RopeError> {
        if num_links < 2 {
            return Err(RopeError::TooFewLinks);
        }
        if mass_per_link <= 0.0 {
            return Err(RopeError::ZeroMassLink(0));
        }

        let links: Vec<RopeLink> = (0..num_links)
            .map(|i| {
                let pos = [root[0], root[1] - i as f64 * segment_length, root[2]];
                RopeLink {
                    position: pos,
                    prev_position: pos,
                    mass: mass_per_link,
                    radius: segment_length * 0.1,
                }
            })
            .collect();

        Ok(Rope {
            links,
            rest_length: segment_length,
            bending_stiffness: 0.0,
            damping: 0.0,
            gravity: [0.0, -9.81, 0.0],
            anchor_head: None,
            anchor_tail: None,
        })
    }

    /// Full constructor from a pre-built link list.
    ///
    /// Returns `Err(RopeError::TooFewLinks)` if fewer than 2 links are provided.
    /// Returns `Err(RopeError::ZeroMassLink(i))` for the first link with mass ≤ 0.
    pub fn from_links(links: Vec<RopeLink>, rest_length: f64) -> Result<Self, RopeError> {
        if links.len() < 2 {
            return Err(RopeError::TooFewLinks);
        }
        for (i, link) in links.iter().enumerate() {
            if link.mass <= 0.0 {
                return Err(RopeError::ZeroMassLink(i));
            }
        }
        Ok(Rope {
            links,
            rest_length,
            bending_stiffness: 0.0,
            damping: 0.0,
            gravity: [0.0, -9.81, 0.0],
            anchor_head: None,
            anchor_tail: None,
        })
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Resolve an anchor to a world-space position, if available.
    fn resolve_anchor(anchor: &AnchorKind, lookup: &BodyTransformLookup<'_>) -> Option<[f64; 3]> {
        match anchor {
            AnchorKind::World(p) => Some(*p),
            AnchorKind::Body(id) => lookup(*id),
        }
    }

    /// Apply head and tail anchor pins (both position and prev_position).
    fn apply_anchors(&mut self, lookup: &BodyTransformLookup<'_>) {
        let n = self.links.len();

        // Clone to avoid borrow conflicts (anchor is part of self, links is another part).
        let head_pos = self
            .anchor_head
            .as_ref()
            .and_then(|a| Self::resolve_anchor(a, lookup));
        let tail_pos = self
            .anchor_tail
            .as_ref()
            .and_then(|a| Self::resolve_anchor(a, lookup));

        if let Some(pos) = head_pos {
            self.links[0].position = pos;
            self.links[0].prev_position = pos;
        }
        if let Some(pos) = tail_pos {
            self.links[n - 1].position = pos;
            self.links[n - 1].prev_position = pos;
        }
    }

    /// Returns (head_pinned, tail_pinned) — whether each end has a resolvable anchor.
    fn pinned_flags(&self, lookup: &BodyTransformLookup<'_>) -> (bool, bool) {
        let head = self
            .anchor_head
            .as_ref()
            .map(|a| Self::resolve_anchor(a, lookup).is_some())
            .unwrap_or(false);
        let tail = self
            .anchor_tail
            .as_ref()
            .map(|a| Self::resolve_anchor(a, lookup).is_some())
            .unwrap_or(false);
        (head, tail)
    }

    // -----------------------------------------------------------------------
    // Step
    // -----------------------------------------------------------------------

    /// Advance the rope by one time step.
    ///
    /// # Steps
    /// 1. Apply anchor overrides (pin).
    /// 2. Velocity Verlet integration for non-pinned links.
    /// 3. Re-apply anchor pins (bodies may have moved).
    /// 4. Two-pass Gauss-Seidel distance constraints.
    /// 5. Bending correction (if `bending_stiffness > 0`).
    /// 6. Final anchor re-pin.
    pub fn step(&mut self, dt: f64, lookup: &BodyTransformLookup<'_>) {
        let n = self.links.len();
        let (head_pinned, tail_pinned) = self.pinned_flags(lookup);

        // Step 1: apply anchors before Verlet so pinned prev_pos is set.
        self.apply_anchors(lookup);

        // Step 2: Velocity Verlet for non-pinned links.
        let grav_acc = self.gravity;
        let damping = self.damping;
        let dt2 = dt * dt;

        for i in 0..n {
            let pinned = (i == 0 && head_pinned) || (i == n - 1 && tail_pinned);
            if pinned {
                continue;
            }
            let link = &mut self.links[i];
            let vel = scale(sub(link.position, link.prev_position), 1.0 - damping);
            let grav = scale(grav_acc, dt2);
            let new_pos = add(add(link.position, vel), grav);
            link.prev_position = link.position;
            link.position = new_pos;
        }

        // Step 3: re-pin anchors (body may have moved).
        self.apply_anchors(lookup);

        // Step 4: two-pass Gauss-Seidel distance constraints.
        for _ in 0..2 {
            for i in 0..n - 1 {
                let pi_pinned = (i == 0 && head_pinned) || (i == n - 1 && tail_pinned);
                // i+1 ranges from 1..n, so i+1==n-1 is the only reachable tail case.
                let pj_pinned = i + 1 == n - 1 && tail_pinned;

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

        // Step 5: bending correction on triplets.
        if self.bending_stiffness > 0.0 {
            let stiffness = self.bending_stiffness;
            for i in 1..n - 1 {
                let pinned = (i == 0 && head_pinned) || (i == n - 1 && tail_pinned);
                if pinned {
                    continue;
                }
                let pa = self.links[i - 1].position;
                let pb = self.links[i].position;
                let pc = self.links[i + 1].position;
                // midpoint of the two neighbours
                let mid = scale(add(pa, pc), 0.5);
                // move the middle link toward the midpoint
                let correction = scale(sub(mid, pb), stiffness);
                self.links[i].position = add(pb, correction);
            }
        }

        // Step 6: final anchor re-pin.
        self.apply_anchors(lookup);
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Iterator over consecutive link-pair positions, suitable for collision queries.
    pub fn iter_segments(&self) -> impl Iterator<Item = ([f64; 3], [f64; 3])> + '_ {
        self.links
            .windows(2)
            .map(|w| (w[0].position, w[1].position))
    }

    /// Sum of all current segment lengths (not rest length).
    pub fn total_length(&self) -> f64 {
        self.links
            .windows(2)
            .map(|w| len(sub(w[1].position, w[0].position)))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn no_lookup() -> impl Fn(u64) -> Option<[f64; 3]> {
        |_: u64| None
    }

    // 1. Free hang: X and Z should not drift (no anchors; whole rope falls straight down).
    #[test]
    fn test_free_hang_straight() {
        // No anchors — rope falls freely under gravity.
        let mut rope = Rope::new([0.0, 5.0, 0.0], 5, 1.0, 1.0).expect("valid rope");
        let lookup = no_lookup();
        let lu: &BodyTransformLookup<'_> = &lookup;
        let dt = 1.0 / 60.0;
        for _ in 0..200 {
            rope.step(dt, lu);
        }
        // Gravity is purely along -Y, rope started on Y-axis:
        // X and Z remain 0 throughout (no lateral force or asymmetry).
        for (i, link) in rope.links.iter().enumerate() {
            assert!(
                link.position[0].abs() < 0.01,
                "link {i} X drifted: {}",
                link.position[0]
            );
            assert!(
                link.position[2].abs() < 0.01,
                "link {i} Z drifted: {}",
                link.position[2]
            );
        }
    }

    // 2. Both ends pinned: all links stay within bounding box.
    #[test]
    fn test_both_ends_pinned() {
        // 6 links: head at [0,5,0], tail at [0,0,0], segment_length 1.
        let mut rope = Rope::new([0.0, 5.0, 0.0], 6, 1.0, 1.0).expect("valid rope");
        rope.anchor_head = Some(AnchorKind::World([0.0, 5.0, 0.0]));
        rope.anchor_tail = Some(AnchorKind::World([0.0, 0.0, 0.0]));

        let lookup = no_lookup();
        let lu: &BodyTransformLookup<'_> = &lookup;
        let dt = 1.0 / 60.0;
        for _ in 0..200 {
            rope.step(dt, lu);
        }
        for (i, link) in rope.links.iter().enumerate() {
            let [x, y, z] = link.position;
            assert!(x.abs() <= 2.0, "link {i} x={x} out of [-2, 2]");
            assert!(
                (-1.0_f64..=5.5).contains(&y),
                "link {i} y={y} out of [-1, 5.5]"
            );
            assert!(z.abs() <= 2.0, "link {i} z={z} out of [-2, 2]");
        }
    }

    // 3. Zero gravity: links at rest should not move.
    #[test]
    fn test_zero_gravity_no_motion() {
        let mut rope = Rope::new([0.0, 0.0, 0.0], 4, 1.0, 1.0).expect("valid rope");
        rope.gravity = [0.0, 0.0, 0.0];

        // Record initial positions.
        let initial: Vec<[f64; 3]> = rope.links.iter().map(|l| l.position).collect();

        let lookup = no_lookup();
        let lu: &BodyTransformLookup<'_> = &lookup;
        let dt = 1.0 / 60.0;
        for _ in 0..100 {
            rope.step(dt, lu);
        }

        let threshold = rope.rest_length * 0.01;
        for (i, (link, &init)) in rope.links.iter().zip(initial.iter()).enumerate() {
            let d = len(sub(link.position, init));
            assert!(
                d < threshold,
                "link {i} moved {d} > threshold {threshold} under zero gravity"
            );
        }
    }

    // 4. Zero-mass rejection.
    #[test]
    fn test_zero_mass_rejection() {
        let result = Rope::new([0.0, 0.0, 0.0], 3, 1.0, 0.0);
        assert!(
            matches!(result, Err(RopeError::ZeroMassLink(_))),
            "Expected ZeroMassLink error, got {result:?}"
        );
    }

    // 5. Too few links.
    #[test]
    fn test_too_few_links() {
        let result = Rope::new([0.0, 0.0, 0.0], 1, 1.0, 1.0);
        assert!(
            matches!(result, Err(RopeError::TooFewLinks)),
            "Expected TooFewLinks error, got {result:?}"
        );
    }

    // 6. Serde round-trip.
    #[test]
    fn test_serde_round_trip() {
        let mut rope = Rope::new([1.0, 2.0, 3.0], 4, 0.5, 2.0).expect("valid rope");
        rope.anchor_head = Some(AnchorKind::World([1.0, 2.0, 3.0]));
        rope.anchor_tail = Some(AnchorKind::Body(99));

        let json = serde_json::to_string(&rope).expect("serialize");
        let restored: Rope = serde_json::from_str(&json).expect("deserialize");

        assert_eq!(restored.links.len(), rope.links.len());
        assert!(
            (restored.rest_length - rope.rest_length).abs() < 1.0e-12,
            "rest_length mismatch"
        );
        assert!(
            matches!(restored.anchor_head, Some(AnchorKind::World(_))),
            "anchor_head not preserved"
        );
        assert!(
            matches!(restored.anchor_tail, Some(AnchorKind::Body(99))),
            "anchor_tail not preserved"
        );
    }
}
