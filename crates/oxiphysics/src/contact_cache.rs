// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Persistent contact pair cache with warm-starting support.
//!
//! Sequential impulse constraint solvers converge significantly faster when
//! initialised with impulse values from the previous simulation step
//! ("warm starting"). This module provides a frame-to-frame contact cache that
//! stores contact manifolds and their associated impulses.
//!
//! ## Life Cycle
//!
//! Each simulation step:
//! 1. Call `ContactCache::begin_step` — ages every entry by one frame.
//! 2. For each active contact pair detected by narrow-phase: call
//!    `ContactCache::update_pair` — inserts or refreshes the entry and resets
//!    its lifetime to 0.
//! 3. Before the impulse solver loop: call `ContactCache::lookup` to retrieve
//!    warm-start impulses.
//! 4. After the solver converges: call `ContactCache::store_impulses` to save
//!    the final impulse values.
//! 5. Call `ContactCache::evict_stale` to remove pairs that disappeared (e.g.
//!    separating bodies).
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::contact_cache::{ContactCache, ContactPoint};
//!
//! let mut cache = ContactCache::new(3);
//!
//! cache.begin_step();
//! let pt = ContactPoint {
//!     pos_a: [0.0, 1.0, 0.0],
//!     pos_b: [0.0, 0.9, 0.0],
//!     normal: [0.0, 1.0, 0.0],
//!     depth: 0.1,
//! };
//! cache.update_pair(0, 1, vec![pt]);
//! cache.evict_stale();
//!
//! let entry = cache.lookup(0, 1).unwrap();
//! assert_eq!(entry.points.len(), 1);
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PairKey — private canonical (a, b) pair helper
// ---------------------------------------------------------------------------

/// A canonicalized body-pair key: min(a, b) < max(a, b).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PairKey(u32, u32);

impl PairKey {
    #[inline]
    fn new(a: u32, b: u32) -> Self {
        if a <= b { PairKey(a, b) } else { PairKey(b, a) }
    }

    #[inline]
    fn matches(&self, entry: &CachedContact) -> bool {
        *self == PairKey::new(entry.body_a, entry.body_b)
    }
}

// ---------------------------------------------------------------------------
// ContactPoint
// ---------------------------------------------------------------------------

/// A single contact point in a collision manifold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactPoint {
    /// World-space position of the contact on body A.
    pub pos_a: [f64; 3],
    /// World-space position of the contact on body B.
    pub pos_b: [f64; 3],
    /// Contact normal, pointing from B toward A.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
}

// ---------------------------------------------------------------------------
// CachedContact
// ---------------------------------------------------------------------------

/// Cached contact data for a body pair, including warm-start impulse values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CachedContact {
    /// Lower body ID (canonical order).
    pub body_a: u32,
    /// Higher body ID (canonical order).
    pub body_b: u32,
    /// Contact points from the last narrow-phase update.
    pub points: Vec<ContactPoint>,
    /// Normal impulse magnitude for each contact point.
    ///
    /// Indexed to match `points`. Populated by [`ContactCache::store_impulses`].
    pub normal_impulse: Vec<f64>,
    /// Tangent (friction) impulse for each contact point: `[u, v]` in the
    /// contact tangent plane.
    pub tangent_impulse: Vec<[f64; 2]>,
    /// Frames since this entry was last refreshed by [`ContactCache::update_pair`].
    pub lifetime: u32,
}

impl CachedContact {
    fn new(a: u32, b: u32, points: Vec<ContactPoint>) -> Self {
        let n = points.len();
        Self {
            body_a: a,
            body_b: b,
            points,
            normal_impulse: vec![0.0; n],
            tangent_impulse: vec![[0.0; 2]; n],
            lifetime: 0,
        }
    }

    /// Warm-start normal impulse slice (matches `points` by index).
    pub fn warm_start_normal(&self) -> &[f64] {
        &self.normal_impulse
    }

    /// Warm-start tangent impulse slice (matches `points` by index).
    pub fn warm_start_tangent(&self) -> &[[f64; 2]] {
        &self.tangent_impulse
    }

    /// Update the stored impulses from solver output.
    ///
    /// Slices are silently truncated or zero-padded to match `points.len()`.
    pub fn update_impulses(&mut self, normal: &[f64], tangent: &[[f64; 2]]) {
        let n = self.points.len();
        self.normal_impulse.resize(n, 0.0);
        self.tangent_impulse.resize(n, [0.0; 2]);
        let nlen = normal.len().min(n);
        self.normal_impulse[..nlen].copy_from_slice(&normal[..nlen]);
        let tlen = tangent.len().min(n);
        self.tangent_impulse[..tlen].copy_from_slice(&tangent[..tlen]);
    }
}

// ---------------------------------------------------------------------------
// ContactCache
// ---------------------------------------------------------------------------

/// Persistent contact pair cache with warm-starting support.
///
/// Internally uses a `Vec` to store entries — contact counts in a typical
/// simulation are bounded (hundreds, not thousands) so linear scan is
/// competitive and avoids HashMap serialisation complexity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactCache {
    entries: Vec<CachedContact>,
    /// Maximum number of frames an entry can go without refresh before
    /// [`evict_stale`][Self::evict_stale] removes it.  Typical value: 2–4.
    pub max_lifetime: u32,
}

impl ContactCache {
    /// Create a new cache with the given maximum entry lifetime in frames.
    pub fn new(max_lifetime: u32) -> Self {
        Self {
            entries: Vec::new(),
            max_lifetime,
        }
    }

    // -----------------------------------------------------------------------
    // Per-step lifecycle
    // -----------------------------------------------------------------------

    /// Increment the lifetime counter for every entry.
    ///
    /// Call this **once at the start** of each simulation step, before
    /// [`update_pair`][Self::update_pair].
    pub fn begin_step(&mut self) {
        for e in &mut self.entries {
            e.lifetime += 1;
        }
    }

    /// Insert or refresh the contact manifold for a body pair.
    ///
    /// - If the pair is **new**, an entry is created with zero impulses.
    /// - If the pair **already exists** and the point count matches, the stored
    ///   impulses are preserved for warm starting.  If the count changed (manifold
    ///   restructured), impulses are reset.
    /// - The entry's `lifetime` is reset to 0.
    pub fn update_pair(&mut self, a: u32, b: u32, points: Vec<ContactPoint>) {
        let key = PairKey::new(a, b);
        let n = points.len();
        if let Some(e) = self.entries.iter_mut().find(|e| key.matches(e)) {
            if e.points.len() != n {
                // Manifold restructured — reset warm-start data.
                e.normal_impulse = vec![0.0; n];
                e.tangent_impulse = vec![[0.0; 2]; n];
            }
            e.points = points;
            e.lifetime = 0;
        } else {
            self.entries.push(CachedContact::new(key.0, key.1, points));
        }
    }

    /// Look up the cached contact data for a body pair (immutable).
    pub fn lookup(&self, a: u32, b: u32) -> Option<&CachedContact> {
        let key = PairKey::new(a, b);
        self.entries.iter().find(|e| key.matches(e))
    }

    /// Look up the cached contact data for a body pair (mutable).
    pub fn lookup_mut(&mut self, a: u32, b: u32) -> Option<&mut CachedContact> {
        let key = PairKey::new(a, b);
        self.entries.iter_mut().find(|e| key.matches(e))
    }

    /// Store final solver impulses for a pair so they can be used as warm-start
    /// values in the next step.
    pub fn store_impulses(&mut self, a: u32, b: u32, normal: &[f64], tangent: &[[f64; 2]]) {
        let key = PairKey::new(a, b);
        if let Some(e) = self.entries.iter_mut().find(|e| key.matches(e)) {
            e.update_impulses(normal, tangent);
        }
    }

    // -----------------------------------------------------------------------
    // Eviction and cleanup
    // -----------------------------------------------------------------------

    /// Remove entries whose `lifetime` exceeds [`max_lifetime`][Self::max_lifetime].
    ///
    /// Returns the number of entries evicted.
    pub fn evict_stale(&mut self) -> usize {
        let max = self.max_lifetime;
        let before = self.entries.len();
        self.entries.retain(|e| e.lifetime <= max);
        before - self.entries.len()
    }

    /// Explicitly remove the entry for a specific body pair.
    pub fn remove(&mut self, a: u32, b: u32) {
        let key = PairKey::new(a, b);
        self.entries.retain(|e| !key.matches(e));
    }

    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Number of contact pairs currently in the cache.
    pub fn entry_count(&self) -> usize {
        self.entries.len()
    }

    /// Return `true` if the cache has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate over all cached contact entries.
    pub fn iter(&self) -> impl Iterator<Item = &CachedContact> {
        self.entries.iter()
    }

    /// Collect pairs whose summed normal impulse exceeds `threshold`.
    ///
    /// Useful for triggering high-impact sound or VFX effects.
    ///
    /// Returns `(body_a, body_b, total_normal_impulse)` triples.
    pub fn pairs_above_impulse_threshold(&self, threshold: f64) -> Vec<(u32, u32, f64)> {
        self.entries
            .iter()
            .filter_map(|e| {
                let sum: f64 = e.normal_impulse.iter().sum();
                if sum > threshold {
                    Some((e.body_a, e.body_b, sum))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Return the body pair with the largest summed normal impulse, if any.
    pub fn highest_impact_pair(&self) -> Option<(u32, u32, f64)> {
        self.entries
            .iter()
            .map(|e| {
                let sum: f64 = e.normal_impulse.iter().sum();
                (e.body_a, e.body_b, sum)
            })
            .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
    }
}
