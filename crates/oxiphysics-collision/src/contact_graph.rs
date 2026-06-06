// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact graph persistence and speculative contacts.
//!
//! Provides persistent contact tracking across frames with warm-starting support,
//! speculative contact prediction for continuous collision avoidance, and a compact
//! LRU contact cache.

use std::collections::HashMap;

// ── ContactKey ───────────────────────────────────────────────────────────────

/// Canonical key for a contact pair. Always stores `body_a < body_b`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ContactKey {
    /// Index of the lower-indexed body.
    pub body_a: usize,
    /// Index of the higher-indexed body.
    pub body_b: usize,
}

impl ContactKey {
    /// Create a canonical `ContactKey`; swaps indices so `body_a < body_b`.
    pub fn new(a: usize, b: usize) -> Self {
        if a <= b {
            Self {
                body_a: a,
                body_b: b,
            }
        } else {
            Self {
                body_a: b,
                body_b: a,
            }
        }
    }
}

// ── PersistedContact ─────────────────────────────────────────────────────────

/// A contact that persists across simulation frames for warm-starting.
#[derive(Debug, Clone)]
pub struct PersistedContact {
    /// The pair of bodies involved.
    pub key: ContactKey,
    /// Contact normal pointing from body_a to body_b.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// World-space contact point.
    pub contact_point: [f64; 3],
    /// Number of consecutive frames this contact has been active.
    pub age: u32,
    /// Accumulated normal impulse for warm starting.
    pub cached_impulse: f64,
    /// Whether this contact was seen in the current frame.
    pub is_active: bool,
}

// ── ContactGraph ─────────────────────────────────────────────────────────────

/// Persistent contact graph: maps body pairs to their latest contact data.
#[derive(Debug, Default)]
pub struct ContactGraph {
    contacts: HashMap<ContactKey, PersistedContact>,
}

impl ContactGraph {
    /// Create an empty `ContactGraph`.
    pub fn new() -> Self {
        Self {
            contacts: HashMap::new(),
        }
    }

    /// Insert or update a contact. Increments `age` on update and marks it active.
    pub fn update(&mut self, key: ContactKey, normal: [f64; 3], depth: f64, point: [f64; 3]) {
        let entry = self.contacts.entry(key).or_insert(PersistedContact {
            key,
            normal,
            depth,
            contact_point: point,
            age: 0,
            cached_impulse: 0.0,
            is_active: true,
        });
        entry.normal = normal;
        entry.depth = depth;
        entry.contact_point = point;
        entry.age += 1;
        entry.is_active = true;
    }

    /// Remove a contact by key.
    pub fn remove(&mut self, key: &ContactKey) {
        self.contacts.remove(key);
    }

    /// Look up a contact.
    pub fn get(&self, key: &ContactKey) -> Option<&PersistedContact> {
        self.contacts.get(key)
    }

    /// Look up a contact mutably.
    pub fn get_mut(&mut self, key: &ContactKey) -> Option<&mut PersistedContact> {
        self.contacts.get_mut(key)
    }

    /// Mark every contact inactive before processing a new frame.
    pub fn mark_inactive_all(&mut self) {
        for c in self.contacts.values_mut() {
            c.is_active = false;
        }
    }

    /// Remove all contacts that are still inactive (not seen this frame).
    pub fn purge_inactive(&mut self) {
        self.contacts.retain(|_, v| v.is_active);
    }

    /// Iterate over active contacts.
    pub fn active_contacts(&self) -> Vec<&PersistedContact> {
        self.contacts.values().filter(|c| c.is_active).collect()
    }

    /// All contacts that involve `body` (active or not).
    pub fn contacts_for_body(&self, body: usize) -> Vec<&PersistedContact> {
        self.contacts
            .values()
            .filter(|c| c.key.body_a == body || c.key.body_b == body)
            .collect()
    }

    /// Total number of contacts stored (active + inactive).
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    /// Returns `true` if there are no contacts stored.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }
}

// ── SpeculativeContact ────────────────────────────────────────────────────────

/// A predicted contact used for speculative collision response.
#[derive(Debug, Clone)]
pub struct SpeculativeContact {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Predicted contact normal (from A toward B).
    pub normal: [f64; 3],
    /// Current gap (positive = separated, negative = penetrating).
    pub separation: f64,
    /// Relative closing velocity along the normal (positive = approaching).
    pub closing_velocity: f64,
}

/// Predict a speculative contact between two spheres within the next `dt` seconds.
///
/// Returns `Some` if the bodies will touch or are already penetrating,
/// or if they are approaching fast enough to close the gap within `dt`.
pub fn speculative_contact(
    pos_a: [f64; 3],
    vel_a: [f64; 3],
    pos_b: [f64; 3],
    vel_b: [f64; 3],
    radius_a: f64,
    radius_b: f64,
    dt: f64,
) -> Option<SpeculativeContact> {
    let dx = pos_b[0] - pos_a[0];
    let dy = pos_b[1] - pos_a[1];
    let dz = pos_b[2] - pos_a[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    let (nx, ny, nz) = if dist > 1e-12 {
        (dx / dist, dy / dist, dz / dist)
    } else {
        (0.0, 1.0, 0.0)
    };

    let separation = dist - (radius_a + radius_b);

    // Relative velocity of B with respect to A, projected onto normal
    let rel_vx = vel_b[0] - vel_a[0];
    let rel_vy = vel_b[1] - vel_a[1];
    let rel_vz = vel_b[2] - vel_a[2];
    // Negative dot = approaching (B moving toward A along normal)
    let rel_v_along_n = rel_vx * nx + rel_vy * ny + rel_vz * nz;
    // closing_velocity > 0 means the gap is shrinking
    let closing_velocity = -rel_v_along_n;

    let will_contact =
        separation <= 0.0 || (closing_velocity > 0.0 && separation < closing_velocity * dt);

    if will_contact {
        Some(SpeculativeContact {
            body_a: 0,
            body_b: 1,
            normal: [nx, ny, nz],
            separation,
            closing_velocity,
        })
    } else {
        None
    }
}

/// Compute the speculative impulse magnitude needed to prevent interpenetration.
///
/// Returns the scalar impulse along the contact normal.
pub fn speculative_impulse(
    contact: &SpeculativeContact,
    inv_mass_a: f64,
    inv_mass_b: f64,
    restitution: f64,
) -> f64 {
    let denom = inv_mass_a + inv_mass_b;
    if denom < 1e-30 {
        return 0.0;
    }
    // closing_velocity > 0 means approaching; impulse must be positive
    let numerator = (1.0 + restitution) * contact.closing_velocity;
    (numerator / denom).max(0.0)
}

// ── ContactCache ──────────────────────────────────────────────────────────────

/// Compact fixed-capacity contact cache with LRU eviction.
///
/// Stores `(ContactKey, cached_normal_impulse)` pairs in insertion order.
/// When full, the oldest entry (front) is evicted.
#[derive(Debug)]
pub struct ContactCache {
    entries: Vec<(ContactKey, f64)>,
    capacity: usize,
}

impl ContactCache {
    /// Create a new cache with the given maximum capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            capacity,
        }
    }

    /// Insert or update a `(key, impulse)` pair. Evicts the oldest entry when full.
    pub fn insert(&mut self, key: ContactKey, impulse: f64) {
        // Update in place if already present
        if let Some(e) = self.entries.iter_mut().find(|e| e.0 == key) {
            e.1 = impulse;
            return;
        }
        // Do nothing if capacity is zero
        if self.capacity == 0 {
            return;
        }
        // Evict oldest (front) if at capacity
        if self.entries.len() >= self.capacity {
            self.entries.remove(0);
        }
        self.entries.push((key, impulse));
    }

    /// Look up a cached impulse by key.
    pub fn lookup(&self, key: &ContactKey) -> Option<f64> {
        self.entries.iter().find(|e| &e.0 == key).map(|e| e.1)
    }

    /// Number of entries currently stored.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` when no entries are stored.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ContactGraph;
    use crate::ContactKey;

    // ── ContactKey tests ──────────────────────────────────────────────────────

    #[test]
    fn contact_key_canonical_order() {
        let k1 = ContactKey::new(3, 1);
        let k2 = ContactKey::new(1, 3);
        assert_eq!(k1, k2);
        assert_eq!(k1.body_a, 1);
        assert_eq!(k1.body_b, 3);
    }

    #[test]
    fn contact_key_equal_indices() {
        let k = ContactKey::new(5, 5);
        assert_eq!(k.body_a, 5);
        assert_eq!(k.body_b, 5);
    }

    #[test]
    fn contact_key_hash_consistent() {
        use std::collections::HashSet;
        let mut s = HashSet::new();
        s.insert(ContactKey::new(0, 1));
        s.insert(ContactKey::new(1, 0));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn contact_key_different_pairs_not_equal() {
        let k1 = ContactKey::new(0, 1);
        let k2 = ContactKey::new(0, 2);
        assert_ne!(k1, k2);
    }

    #[test]
    fn contact_key_zero_indices() {
        let k = ContactKey::new(0, 0);
        assert_eq!(k.body_a, 0);
        assert_eq!(k.body_b, 0);
    }

    // ── ContactGraph tests ────────────────────────────────────────────────────

    #[test]
    fn contact_graph_new_is_empty() {
        let g = ContactGraph::new();
        assert!(g.is_empty());
        assert_eq!(g.len(), 0);
    }

    #[test]
    fn contact_graph_update_inserts() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.1, [0.0, 0.0, 0.0]);
        assert_eq!(g.len(), 1);
    }

    #[test]
    fn contact_graph_update_increments_age() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.1, [0.0, 0.0, 0.0]);
        g.update(key, [0.0, 1.0, 0.0], 0.2, [0.0, 0.0, 0.0]);
        let c = g.get(&key).unwrap();
        assert_eq!(c.age, 2);
    }

    #[test]
    fn contact_graph_update_overwrites_depth() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.5, [0.0, 0.0, 0.0]);
        g.update(key, [0.0, 1.0, 0.0], 1.5, [0.0, 0.0, 0.0]);
        assert!((g.get(&key).unwrap().depth - 1.5).abs() < 1e-12);
    }

    #[test]
    fn contact_graph_remove() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.1, [0.0, 0.0, 0.0]);
        g.remove(&key);
        assert!(g.is_empty());
    }

    #[test]
    fn contact_graph_get_missing_returns_none() {
        let g = ContactGraph::new();
        let key = ContactKey::new(7, 8);
        assert!(g.get(&key).is_none());
    }

    #[test]
    fn contact_graph_get_mut_modifies() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.1, [0.0, 0.0, 0.0]);
        g.get_mut(&key).unwrap().cached_impulse = 3.125;
        assert!((g.get(&key).unwrap().cached_impulse - 3.125).abs() < 1e-12);
    }

    #[test]
    fn contact_graph_mark_inactive_all() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(1, 2), [0.0, 1.0, 0.0], 0.2, [0.0; 3]);
        g.mark_inactive_all();
        assert_eq!(g.active_contacts().len(), 0);
    }

    #[test]
    fn contact_graph_purge_inactive_removes_stale() {
        let mut g = ContactGraph::new();
        let k1 = ContactKey::new(0, 1);
        let k2 = ContactKey::new(1, 2);
        g.update(k1, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(k2, [0.0, 1.0, 0.0], 0.2, [0.0; 3]);
        g.mark_inactive_all();
        // Re-activate only k1
        g.update(k1, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.purge_inactive();
        assert_eq!(g.len(), 1);
        assert!(g.get(&k1).is_some());
        assert!(g.get(&k2).is_none());
    }

    #[test]
    fn contact_graph_active_contacts_count() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(2, 3), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.mark_inactive_all();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        assert_eq!(g.active_contacts().len(), 1);
    }

    #[test]
    fn contact_graph_contacts_for_body() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(0, 2), [0.0, 1.0, 0.0], 0.2, [0.0; 3]);
        g.update(ContactKey::new(3, 4), [0.0, 1.0, 0.0], 0.3, [0.0; 3]);
        let contacts = g.contacts_for_body(0);
        assert_eq!(contacts.len(), 2);
        let contacts_4 = g.contacts_for_body(4);
        assert_eq!(contacts_4.len(), 1);
    }

    #[test]
    fn contact_graph_is_active_on_update() {
        let mut g = ContactGraph::new();
        let key = ContactKey::new(0, 1);
        g.update(key, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        assert!(g.get(&key).unwrap().is_active);
    }

    // ── SpeculativeContact tests ──────────────────────────────────────────────

    #[test]
    fn speculative_contact_approaching_predicts_contact() {
        // Two spheres 1 unit apart (radii 0.4 each → gap = 0.2), approaching at 1 m/s
        let pos_a = [0.0, 0.0, 0.0];
        let vel_a = [0.5, 0.0, 0.0];
        let pos_b = [1.0, 0.0, 0.0];
        let vel_b = [-0.5, 0.0, 0.0];
        let result = speculative_contact(pos_a, vel_a, pos_b, vel_b, 0.4, 0.4, 1.0 / 60.0);
        // closing velocity = 1.0, gap = 0.2, 1.0 * (1/60) ≈ 0.0167 < 0.2 → no prediction this frame
        assert!(result.is_none());
    }

    #[test]
    fn speculative_contact_fast_approach_within_dt() {
        // Gap = 0.01, closing at 10 m/s, dt = 1/60 → 10 * 1/60 ≈ 0.167 > 0.01
        let pos_a = [0.0, 0.0, 0.0];
        let vel_a = [5.0, 0.0, 0.0];
        let pos_b = [1.0, 0.0, 0.0];
        let vel_b = [-5.0, 0.0, 0.0];
        // radii sum = 0.99 → separation = 1.0 - 0.99 = 0.01
        let result = speculative_contact(pos_a, vel_a, pos_b, vel_b, 0.495, 0.495, 1.0 / 60.0);
        assert!(result.is_some());
        let c = result.unwrap();
        assert!(c.closing_velocity > 0.0);
    }

    #[test]
    fn speculative_contact_penetrating_always_some() {
        // Already overlapping (dist < r_a + r_b)
        let pos_a = [0.0, 0.0, 0.0];
        let pos_b = [0.5, 0.0, 0.0];
        let result = speculative_contact(pos_a, [0.0; 3], pos_b, [0.0; 3], 0.4, 0.4, 1.0 / 60.0);
        assert!(result.is_some());
        assert!(result.unwrap().separation < 0.0);
    }

    #[test]
    fn speculative_contact_separating_returns_none() {
        let pos_a = [0.0, 0.0, 0.0];
        let vel_a = [-5.0, 0.0, 0.0];
        let pos_b = [2.0, 0.0, 0.0];
        let vel_b = [5.0, 0.0, 0.0];
        let result = speculative_contact(pos_a, vel_a, pos_b, vel_b, 0.1, 0.1, 1.0 / 60.0);
        assert!(result.is_none());
    }

    #[test]
    fn speculative_impulse_zero_for_static_bodies() {
        let c = SpeculativeContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            separation: -0.05,
            closing_velocity: 1.0,
        };
        let j = speculative_impulse(&c, 0.0, 0.0, 0.5);
        assert!((j).abs() < 1e-12);
    }

    #[test]
    fn speculative_impulse_positive_for_approaching() {
        let c = SpeculativeContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            separation: -0.05,
            closing_velocity: 2.0,
        };
        let j = speculative_impulse(&c, 1.0, 1.0, 0.0);
        assert!(j > 0.0);
    }

    #[test]
    fn speculative_impulse_no_negative() {
        // Separating contact → impulse clamped to 0
        let c = SpeculativeContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            separation: 0.5,
            closing_velocity: -1.0,
        };
        let j = speculative_impulse(&c, 1.0, 1.0, 0.5);
        assert!(j >= 0.0);
    }

    // ── ContactCache tests ────────────────────────────────────────────────────

    #[test]
    fn contact_cache_insert_and_lookup() {
        let mut cache = ContactCache::new(8);
        let key = ContactKey::new(0, 1);
        cache.insert(key, 1.5);
        assert_eq!(cache.lookup(&key), Some(1.5));
    }

    #[test]
    fn contact_cache_lookup_missing_returns_none() {
        let cache = ContactCache::new(8);
        assert_eq!(cache.lookup(&ContactKey::new(9, 10)), None);
    }

    #[test]
    fn contact_cache_update_existing() {
        let mut cache = ContactCache::new(8);
        let key = ContactKey::new(0, 1);
        cache.insert(key, 1.0);
        cache.insert(key, 2.0);
        assert_eq!(cache.lookup(&key), Some(2.0));
        assert_eq!(cache.len(), 1);
    }

    #[test]
    fn contact_cache_lru_eviction() {
        let mut cache = ContactCache::new(3);
        let k0 = ContactKey::new(0, 1);
        let k1 = ContactKey::new(1, 2);
        let k2 = ContactKey::new(2, 3);
        let k3 = ContactKey::new(3, 4);
        cache.insert(k0, 0.0);
        cache.insert(k1, 1.0);
        cache.insert(k2, 2.0);
        // k0 should be evicted
        cache.insert(k3, 3.0);
        assert_eq!(cache.len(), 3);
        assert!(cache.lookup(&k0).is_none());
        assert_eq!(cache.lookup(&k3), Some(3.0));
    }

    #[test]
    fn contact_cache_is_empty() {
        let cache = ContactCache::new(4);
        assert!(cache.is_empty());
    }

    #[test]
    fn contact_cache_zero_capacity_does_not_panic() {
        let mut cache = ContactCache::new(0);
        let key = ContactKey::new(0, 1);
        cache.insert(key, 1.0); // should not panic
        assert!(cache.is_empty());
    }
}

// ── ContactManifold (up to 4 contact points) ─────────────────────────────────

/// A single contact point within a manifold.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// World-space position of the contact point.
    pub position: [f64; 3],
    /// Penetration depth at this contact point.
    pub depth: f64,
    /// Accumulated normal impulse for warm starting.
    pub impulse: f64,
    /// Whether this point was matched to a previous frame's point.
    pub warm_started: bool,
}

impl ContactPoint {
    /// Create a new contact point.
    pub fn new(position: [f64; 3], depth: f64) -> Self {
        Self {
            position,
            depth,
            impulse: 0.0,
            warm_started: false,
        }
    }
}

/// A persistent contact manifold holding up to 4 contact points.
///
/// The manifold is associated with a body pair and stores the contact
/// normal and up to 4 contact points for stable simulation.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Unique key for the body pair.
    pub key: ContactKey,
    /// Contact normal (from body_b toward body_a).
    pub normal: [f64; 3],
    /// Contact points (max 4).
    pub points: Vec<ContactPoint>,
    /// Number of frames this manifold has been active.
    pub age: u32,
    /// Whether this manifold was refreshed this frame.
    pub active: bool,
}

impl ContactManifold {
    /// Create a new empty manifold.
    pub fn new(key: ContactKey, normal: [f64; 3]) -> Self {
        Self {
            key,
            normal,
            points: Vec::with_capacity(4),
            age: 0,
            active: true,
        }
    }

    /// Add a contact point. If already at 4, reduce by removing the point
    /// that minimises the manifold area (contact reduction).
    pub fn add_point(&mut self, point: ContactPoint) {
        if self.points.len() < 4 {
            self.points.push(point);
        } else {
            // Contact reduction: keep the deepest + 3 that maximise area
            self.reduce_contacts(point);
        }
    }

    /// Contact reduction heuristic: replace the point that minimises the
    /// area of the remaining 4-point set.
    fn reduce_contacts(&mut self, new_point: ContactPoint) {
        // Candidate set: existing 4 + new_point = 5 points.
        // For each removal candidate, compute the triangle area of the remaining 4.
        let mut candidates = self.points.clone();
        candidates.push(new_point);

        let mut max_area = -1.0_f64;
        let mut best_set = self.points.clone();

        // Deepest point is always kept.
        let deepest_idx = candidates
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.depth
                    .partial_cmp(&b.depth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Try all combinations of 3 from the remaining 4.
        let remaining: Vec<usize> = (0..candidates.len())
            .filter(|&i| i != deepest_idx)
            .collect();
        for skip in &remaining {
            let subset: Vec<&ContactPoint> = candidates
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != deepest_idx && i != skip)
                .map(|(_, p)| p)
                .collect();
            if subset.len() < 3 {
                continue;
            }
            let area = triangle_area(subset[0].position, subset[1].position, subset[2].position);
            if area > max_area {
                max_area = area;
                let mut new_set = vec![candidates[deepest_idx]];
                new_set.extend(subset.iter().copied().copied());
                best_set = new_set;
            }
        }

        self.points = best_set;
    }

    /// Mark the manifold as inactive (to be purged if not refreshed).
    pub fn mark_inactive(&mut self) {
        self.active = false;
    }

    /// Clear all contact points.
    pub fn clear_points(&mut self) {
        self.points.clear();
    }

    /// Increment the age counter.
    pub fn tick(&mut self) {
        self.age += 1;
    }

    /// Number of contact points.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
}

/// Compute the area of a triangle from 3 world-space points.
fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
}

// ── ManifoldStore ──────────────────────────────────────────────────────────────

/// Stores all persistent contact manifolds.
#[derive(Debug, Default)]
pub struct ManifoldStore {
    manifolds: HashMap<ContactKey, ContactManifold>,
}

impl ManifoldStore {
    /// Create an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Get or create a manifold for the given body pair.
    pub fn get_or_create(&mut self, key: ContactKey, normal: [f64; 3]) -> &mut ContactManifold {
        self.manifolds
            .entry(key)
            .or_insert_with(|| ContactManifold::new(key, normal))
    }

    /// Get a manifold immutably.
    pub fn get(&self, key: &ContactKey) -> Option<&ContactManifold> {
        self.manifolds.get(key)
    }

    /// Mark all manifolds inactive at the start of a frame.
    pub fn mark_all_inactive(&mut self) {
        for m in self.manifolds.values_mut() {
            m.mark_inactive();
        }
    }

    /// Remove all inactive manifolds.
    pub fn purge_inactive(&mut self) {
        self.manifolds.retain(|_, m| m.active);
    }

    /// Tick all manifolds (increment age).
    pub fn tick_all(&mut self) {
        for m in self.manifolds.values_mut() {
            m.tick();
        }
    }

    /// Total number of manifolds.
    pub fn len(&self) -> usize {
        self.manifolds.len()
    }

    /// Whether the store is empty.
    pub fn is_empty(&self) -> bool {
        self.manifolds.is_empty()
    }

    /// Collect all active manifolds.
    pub fn active_manifolds(&self) -> Vec<&ContactManifold> {
        self.manifolds.values().filter(|m| m.active).collect()
    }
}

// ── Contact normal/tangent frame ──────────────────────────────────────────────

/// An orthogonal contact frame: normal + two tangents.
#[derive(Debug, Clone, Copy)]
pub struct ContactFrame {
    /// Contact normal (from B toward A).
    pub normal: [f64; 3],
    /// First tangent (perpendicular to normal).
    pub tangent1: [f64; 3],
    /// Second tangent (perpendicular to both normal and tangent1).
    pub tangent2: [f64; 3],
}

impl ContactFrame {
    /// Build a contact frame from a normal vector.
    ///
    /// Uses the Gram-Schmidt process to find two perpendicular tangents.
    pub fn from_normal(n: [f64; 3]) -> Self {
        // Normalize
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let normal = if len > 1e-12 {
            [n[0] / len, n[1] / len, n[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };

        // Choose a reference vector not parallel to normal
        let ref_vec = if normal[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };

        // Gram-Schmidt: t1 = ref - (ref·n)*n
        let dot_rn = ref_vec[0] * normal[0] + ref_vec[1] * normal[1] + ref_vec[2] * normal[2];
        let t1_raw = [
            ref_vec[0] - dot_rn * normal[0],
            ref_vec[1] - dot_rn * normal[1],
            ref_vec[2] - dot_rn * normal[2],
        ];
        let t1_len = (t1_raw[0] * t1_raw[0] + t1_raw[1] * t1_raw[1] + t1_raw[2] * t1_raw[2]).sqrt();
        let tangent1 = if t1_len > 1e-12 {
            [t1_raw[0] / t1_len, t1_raw[1] / t1_len, t1_raw[2] / t1_len]
        } else {
            [0.0, 0.0, 1.0]
        };

        // t2 = n × t1
        let tangent2 = [
            normal[1] * tangent1[2] - normal[2] * tangent1[1],
            normal[2] * tangent1[0] - normal[0] * tangent1[2],
            normal[0] * tangent1[1] - normal[1] * tangent1[0],
        ];

        Self {
            normal,
            tangent1,
            tangent2,
        }
    }

    /// Project a vector onto the contact frame.
    pub fn project(&self, v: [f64; 3]) -> [f64; 3] {
        let vn = v[0] * self.normal[0] + v[1] * self.normal[1] + v[2] * self.normal[2];
        let vt1 = v[0] * self.tangent1[0] + v[1] * self.tangent1[1] + v[2] * self.tangent1[2];
        let vt2 = v[0] * self.tangent2[0] + v[1] * self.tangent2[1] + v[2] * self.tangent2[2];
        [vn, vt1, vt2]
    }

    /// Unproject from the contact frame back to world space.
    pub fn unproject(&self, v: [f64; 3]) -> [f64; 3] {
        [
            v[0] * self.normal[0] + v[1] * self.tangent1[0] + v[2] * self.tangent2[0],
            v[0] * self.normal[1] + v[1] * self.tangent1[1] + v[2] * self.tangent2[1],
            v[0] * self.normal[2] + v[1] * self.tangent1[2] + v[2] * self.tangent2[2],
        ]
    }
}

// ── Cache invalidation by AABB movement ───────────────────────────────────────

/// Threshold for contact cache invalidation.
///
/// When a body moves by more than `threshold * contact_radius`, cached contacts
/// for that body are invalidated.
pub struct ContactCacheInvalidator {
    /// Previous positions of each body.
    pub prev_positions: Vec<[f64; 3]>,
    /// Movement threshold.
    pub threshold: f64,
}

impl ContactCacheInvalidator {
    /// Create a new invalidator.
    pub fn new(n_bodies: usize, threshold: f64) -> Self {
        Self {
            prev_positions: vec![[0.0; 3]; n_bodies],
            threshold,
        }
    }

    /// Update positions and return list of bodies that moved beyond threshold.
    pub fn update_and_check(
        &mut self,
        new_positions: &[[f64; 3]],
        contact_radius: f64,
    ) -> Vec<usize> {
        let limit = self.threshold * contact_radius;
        let moved: Vec<usize> = new_positions
            .iter()
            .enumerate()
            .filter(|(i, np)| {
                let pp = self.prev_positions[*i];
                let d = [np[0] - pp[0], np[1] - pp[1], np[2] - pp[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() > limit
            })
            .map(|(i, _)| i)
            .collect();
        let copy_len = new_positions.len().min(self.prev_positions.len());
        self.prev_positions[..copy_len].copy_from_slice(&new_positions[..copy_len]);
        moved
    }

    /// Invalidate contacts in a graph for all moved bodies.
    pub fn invalidate_contacts(&self, graph: &mut ContactGraph, moved_bodies: &[usize]) {
        for &body in moved_bodies {
            let keys_to_remove: Vec<ContactKey> = graph
                .contacts_for_body(body)
                .iter()
                .map(|c| c.key)
                .collect();
            for key in keys_to_remove {
                graph.remove(&key);
            }
        }
    }
}

// ── Island Detection (BFS-based connected components) ─────────────────────────

/// A simulation island: a group of bodies connected by active contacts.
///
/// Islands are used to split the constraint solver into independent groups,
/// enabling sleeping and parallel solving.
#[derive(Debug, Clone)]
pub struct ContactIsland {
    /// Indices of bodies belonging to this island.
    pub bodies: Vec<usize>,
    /// Contact keys that form the internal edges of the island.
    pub edges: Vec<ContactKey>,
    /// Whether every body in the island could be put to sleep.
    pub can_sleep: bool,
}

impl ContactIsland {
    /// Number of bodies in the island.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Number of contact edges in the island.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
}

/// Statistics about the current island decomposition.
#[derive(Debug, Clone, Default)]
pub struct IslandStats {
    /// Total number of islands.
    pub island_count: usize,
    /// Average island size (bodies per island).
    pub avg_size: f64,
    /// Maximum island size.
    pub max_size: usize,
    /// Number of islands marked as sleeping.
    pub sleeping_islands: usize,
    /// Total number of bodies assigned to an island.
    pub total_bodies: usize,
}

impl IslandStats {
    /// Compute statistics from a slice of islands.
    pub fn from_islands(islands: &[ContactIsland]) -> Self {
        if islands.is_empty() {
            return Self::default();
        }
        let total: usize = islands.iter().map(|i| i.bodies.len()).sum();
        let max_size = islands.iter().map(|i| i.bodies.len()).max().unwrap_or(0);
        let sleeping = islands.iter().filter(|i| i.can_sleep).count();
        Self {
            island_count: islands.len(),
            avg_size: total as f64 / islands.len() as f64,
            max_size,
            sleeping_islands: sleeping,
            total_bodies: total,
        }
    }
}

/// Detect connected components in the contact graph using BFS.
///
/// Only **active** contacts are traversed.  Each returned [`ContactIsland`]
/// contains the body indices and contact keys of one connected component.
///
/// `n_bodies` is the total number of bodies; body indices must be in
/// `0..n_bodies`.
pub fn detect_islands(graph: &ContactGraph, n_bodies: usize) -> Vec<ContactIsland> {
    // Build adjacency: body → list of (neighbour, key)
    let mut adj: Vec<Vec<(usize, ContactKey)>> = vec![Vec::new(); n_bodies];
    for c in graph.active_contacts() {
        let a = c.key.body_a;
        let b = c.key.body_b;
        if a < n_bodies && b < n_bodies {
            adj[a].push((b, c.key));
            adj[b].push((a, c.key));
        }
    }

    let mut visited = vec![false; n_bodies];
    let mut islands: Vec<ContactIsland> = Vec::new();

    for start in 0..n_bodies {
        if visited[start] || adj[start].is_empty() {
            continue;
        }
        // BFS
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        let mut bodies = Vec::new();
        let mut edges: std::collections::HashSet<ContactKey> = std::collections::HashSet::new();

        while let Some(cur) = queue.pop_front() {
            bodies.push(cur);
            for &(neighbour, key) in &adj[cur] {
                edges.insert(key);
                if !visited[neighbour] {
                    visited[neighbour] = true;
                    queue.push_back(neighbour);
                }
            }
        }

        islands.push(ContactIsland {
            bodies,
            edges: edges.into_iter().collect(),
            can_sleep: false,
        });
    }

    islands
}

/// Propagate sleep state across islands.
///
/// An island may sleep if *all* bodies in it have a velocity magnitude below
/// `sleep_threshold` for at least `min_sleep_frames` consecutive frames.
///
/// `velocities[i]` is the squared speed of body `i`.
/// `sleep_counters[i]` tracks how many frames body `i` has been slow.
///
/// Returns a list of body indices that should be put to sleep this frame.
pub fn propagate_sleep(
    islands: &mut [ContactIsland],
    velocities_sq: &[f64],
    sleep_counters: &mut [u32],
    sleep_threshold_sq: f64,
    min_sleep_frames: u32,
) -> Vec<usize> {
    let mut bodies_to_sleep = Vec::new();

    for island in islands.iter_mut() {
        // Check if all bodies are slow
        let all_slow = island
            .bodies
            .iter()
            .all(|&b| velocities_sq.get(b).copied().unwrap_or(f64::INFINITY) < sleep_threshold_sq);

        if all_slow {
            // First increment all counters, then check if all have reached the threshold
            for &b in &island.bodies {
                if let Some(cnt) = sleep_counters.get_mut(b) {
                    *cnt += 1;
                }
            }
            let ready = island
                .bodies
                .iter()
                .all(|&b| sleep_counters.get(b).copied().unwrap_or(0) >= min_sleep_frames);
            if ready {
                island.can_sleep = true;
                bodies_to_sleep.extend(island.bodies.iter().copied());
            }
        } else {
            // Reset counters for any fast body, wake island
            for &b in &island.bodies {
                if let Some(cnt) = sleep_counters.get_mut(b) {
                    *cnt = 0;
                }
            }
            island.can_sleep = false;
        }
    }

    bodies_to_sleep
}

// ── Contact frequency tracking ────────────────────────────────────────────────

/// Category of a contact based on how long it has been active.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactFrequency {
    /// Contact exists for only 1 frame (just started).
    Transient,
    /// Contact has been active for 2–4 frames.
    ShortLived,
    /// Contact has been active for 5+ frames (persistent warm contact).
    Persistent,
}

impl ContactFrequency {
    /// Classify a contact by its `age` field.
    pub fn classify(age: u32) -> Self {
        match age {
            0..=1 => ContactFrequency::Transient,
            2..=4 => ContactFrequency::ShortLived,
            _ => ContactFrequency::Persistent,
        }
    }
}

impl ContactGraph {
    /// Classify every active contact by frequency and return counts.
    ///
    /// Returns `(transient, short_lived, persistent)`.
    pub fn frequency_counts(&self) -> (usize, usize, usize) {
        let mut t = 0usize;
        let mut s = 0usize;
        let mut p = 0usize;
        for c in self.active_contacts() {
            match ContactFrequency::classify(c.age) {
                ContactFrequency::Transient => t += 1,
                ContactFrequency::ShortLived => s += 1,
                ContactFrequency::Persistent => p += 1,
            }
        }
        (t, s, p)
    }

    /// Return all active contacts classified as persistent (age ≥ 5).
    pub fn persistent_contacts(&self) -> Vec<&PersistedContact> {
        self.active_contacts()
            .into_iter()
            .filter(|c| {
                matches!(
                    ContactFrequency::classify(c.age),
                    ContactFrequency::Persistent
                )
            })
            .collect()
    }
}

// ── Manifold quality scoring ──────────────────────────────────────────────────

/// A quality score for a contact manifold in `[0.0, 1.0]`.
///
/// Higher is better.  The score is based on the depth and age of the contact.
pub fn manifold_quality_score(contact: &PersistedContact) -> f64 {
    // Depth contribution: clamp to [0, 1] (deeper = better, up to 0.5 units)
    let depth_score = (contact.depth / 0.5).clamp(0.0, 1.0);
    // Age contribution: older contacts are more reliable (saturates at age 10)
    let age_score = (contact.age as f64 / 10.0).min(1.0);
    // Warm-start contribution: non-zero impulse indicates solver has converged
    let impulse_score = if contact.cached_impulse.abs() > 1e-6 {
        0.2
    } else {
        0.0
    };
    (0.5 * depth_score + 0.3 * age_score + 0.2 * impulse_score).min(1.0)
}

// ── Bipartite static-dynamic boundary check ───────────────────────────────────

/// Test whether the contact graph is bipartite across a static/dynamic boundary.
///
/// A contact graph is "bipartite-safe" if no two static bodies share a direct
/// contact edge (static–static contacts are usually filtered by the broadphase
/// but this function can verify that invariant).
///
/// `is_static[i]` is `true` when body `i` is a static (infinite mass) body.
///
/// Returns `true` if no static–static edge exists.
pub fn is_static_dynamic_bipartite(graph: &ContactGraph, is_static: &[bool]) -> bool {
    for c in graph.active_contacts() {
        let a_static = is_static.get(c.key.body_a).copied().unwrap_or(false);
        let b_static = is_static.get(c.key.body_b).copied().unwrap_or(false);
        if a_static && b_static {
            return false;
        }
    }
    true
}

// ── Tests for new functionality ────────────────────────────────────────────────

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ── ContactPoint tests ────────────────────────────────────────────────────

    #[test]
    fn contact_point_new() {
        let cp = ContactPoint::new([1.0, 2.0, 3.0], 0.05);
        assert!((cp.depth - 0.05).abs() < 1e-12);
        assert!(!cp.warm_started);
        assert_eq!(cp.impulse, 0.0);
    }

    // ── ContactManifold tests ─────────────────────────────────────────────────

    #[test]
    fn contact_manifold_add_up_to_four() {
        let key = ContactKey::new(0, 1);
        let mut m = ContactManifold::new(key, [0.0, 1.0, 0.0]);
        for i in 0..4 {
            m.add_point(ContactPoint::new([i as f64, 0.0, 0.0], 0.1));
        }
        assert_eq!(m.n_points(), 4);
    }

    #[test]
    fn contact_manifold_reduction_keeps_four() {
        let key = ContactKey::new(0, 1);
        let mut m = ContactManifold::new(key, [0.0, 1.0, 0.0]);
        for i in 0..5 {
            m.add_point(ContactPoint::new(
                [i as f64 * 0.1, 0.0, 0.0],
                (i as f64 + 1.0) * 0.01,
            ));
        }
        assert_eq!(m.n_points(), 4, "Should keep exactly 4 after reduction");
    }

    #[test]
    fn contact_manifold_tick_increments_age() {
        let key = ContactKey::new(0, 1);
        let mut m = ContactManifold::new(key, [0.0, 1.0, 0.0]);
        m.tick();
        m.tick();
        assert_eq!(m.age, 2);
    }

    #[test]
    fn contact_manifold_mark_inactive() {
        let key = ContactKey::new(0, 1);
        let mut m = ContactManifold::new(key, [0.0, 1.0, 0.0]);
        assert!(m.active);
        m.mark_inactive();
        assert!(!m.active);
    }

    // ── ManifoldStore tests ───────────────────────────────────────────────────

    #[test]
    fn manifold_store_get_or_create() {
        let mut store = ManifoldStore::new();
        let key = ContactKey::new(0, 1);
        let m = store.get_or_create(key, [0.0, 1.0, 0.0]);
        m.add_point(ContactPoint::new([0.0; 3], 0.1));
        assert_eq!(store.get(&key).unwrap().n_points(), 1);
    }

    #[test]
    fn manifold_store_purge_inactive() {
        let mut store = ManifoldStore::new();
        let k1 = ContactKey::new(0, 1);
        let k2 = ContactKey::new(2, 3);
        store.get_or_create(k1, [0.0, 1.0, 0.0]);
        store.get_or_create(k2, [0.0, 1.0, 0.0]);
        store.mark_all_inactive();
        // Reactivate k1
        store.get_or_create(k1, [0.0, 1.0, 0.0]).active = true;
        store.purge_inactive();
        assert_eq!(store.len(), 1);
        assert!(store.get(&k1).is_some());
        assert!(store.get(&k2).is_none());
    }

    #[test]
    fn manifold_store_active_count() {
        let mut store = ManifoldStore::new();
        store
            .get_or_create(ContactKey::new(0, 1), [0.0, 1.0, 0.0])
            .active = true;
        store
            .get_or_create(ContactKey::new(2, 3), [0.0, 1.0, 0.0])
            .active = false;
        assert_eq!(store.active_manifolds().len(), 1);
    }

    // ── ContactFrame tests ────────────────────────────────────────────────────

    #[test]
    fn contact_frame_normal_preserved() {
        let f = ContactFrame::from_normal([0.0, 1.0, 0.0]);
        let len =
            (f.normal[0] * f.normal[0] + f.normal[1] * f.normal[1] + f.normal[2] * f.normal[2])
                .sqrt();
        assert!((len - 1.0).abs() < 1e-10);
        assert!((f.normal[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn contact_frame_tangents_perpendicular_to_normal() {
        let f = ContactFrame::from_normal([1.0, 0.0, 0.0]);
        let dt1n =
            f.tangent1[0] * f.normal[0] + f.tangent1[1] * f.normal[1] + f.tangent1[2] * f.normal[2];
        let dt2n =
            f.tangent2[0] * f.normal[0] + f.tangent2[1] * f.normal[1] + f.tangent2[2] * f.normal[2];
        assert!(
            dt1n.abs() < 1e-10,
            "tangent1 not perpendicular to normal: {dt1n}"
        );
        assert!(
            dt2n.abs() < 1e-10,
            "tangent2 not perpendicular to normal: {dt2n}"
        );
    }

    #[test]
    fn contact_frame_project_unproject_roundtrip() {
        let f = ContactFrame::from_normal([0.0, 0.0, 1.0]);
        let v = [1.0, 2.0, 3.0];
        let projected = f.project(v);
        let recovered = f.unproject(projected);
        for k in 0..3 {
            assert!(
                (recovered[k] - v[k]).abs() < 1e-10,
                "component {k}: {} vs {}",
                recovered[k],
                v[k]
            );
        }
    }

    // ── ContactCacheInvalidator tests ─────────────────────────────────────────

    #[test]
    fn invalidator_detects_moved_body() {
        let mut inv = ContactCacheInvalidator::new(3, 0.5);
        let new_pos = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let moved = inv.update_and_check(&new_pos, 1.0);
        // Body 0 moved by 1.0 > 0.5 * 1.0 = 0.5
        assert!(moved.contains(&0), "body 0 should be in moved list");
    }

    #[test]
    fn invalidator_static_body_not_moved() {
        let mut inv = ContactCacheInvalidator::new(2, 0.5);
        inv.prev_positions[0] = [0.1, 0.0, 0.0];
        let new_pos = [[0.1, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let moved = inv.update_and_check(&new_pos, 1.0);
        assert!(
            !moved.contains(&0),
            "static body should not be in moved list"
        );
    }

    #[test]
    fn invalidator_removes_contacts_for_moved_body() {
        let mut inv = ContactCacheInvalidator::new(3, 0.1);
        let mut graph = ContactGraph::new();
        let k01 = ContactKey::new(0, 1);
        let k12 = ContactKey::new(1, 2);
        graph.update(k01, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        graph.update(k12, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);

        // Move body 0 beyond threshold
        let new_pos = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let moved = inv.update_and_check(&new_pos, 1.0);
        inv.invalidate_contacts(&mut graph, &moved);

        // Contact (0,1) should be removed; (1,2) stays
        assert!(graph.get(&k01).is_none(), "contact 0-1 should be removed");
        assert!(graph.get(&k12).is_some(), "contact 1-2 should stay");
    }

    // ── Triangle area helper ──────────────────────────────────────────────────

    #[test]
    fn triangle_area_unit() {
        // Right triangle with legs of length 1
        let area = triangle_area([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((area - 0.5).abs() < 1e-12, "area={area}");
    }

    #[test]
    fn triangle_area_degenerate_zero() {
        // Three collinear points → zero area
        let area = triangle_area([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!(area < 1e-12, "degenerate triangle area={area}");
    }
}

// ── Tests for island detection and related features ───────────────────────────

#[cfg(test)]
mod tests_islands {
    use super::*;

    fn make_graph_chain(n: usize) -> ContactGraph {
        let mut g = ContactGraph::new();
        for i in 0..(n - 1) {
            g.update(ContactKey::new(i, i + 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        }
        g
    }

    // ── detect_islands ────────────────────────────────────────────────────────

    #[test]
    fn island_chain_of_four_is_one_island() {
        let g = make_graph_chain(4);
        let islands = detect_islands(&g, 4);
        assert_eq!(islands.len(), 1, "chain of 4 should be 1 island");
        assert_eq!(islands[0].bodies.len(), 4);
    }

    #[test]
    fn island_two_separate_pairs() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(2, 3), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let islands = detect_islands(&g, 4);
        assert_eq!(
            islands.len(),
            2,
            "two disconnected pairs should give 2 islands"
        );
    }

    #[test]
    fn island_empty_graph_gives_no_islands() {
        let g = ContactGraph::new();
        let islands = detect_islands(&g, 5);
        assert_eq!(islands.len(), 0);
    }

    #[test]
    fn island_single_contact() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.05, [0.0; 3]);
        let islands = detect_islands(&g, 2);
        assert_eq!(islands.len(), 1);
        assert_eq!(islands[0].edges.len(), 1);
    }

    #[test]
    fn island_star_topology() {
        // Body 0 touches bodies 1, 2, 3 → all in one island
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(0, 2), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(0, 3), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let islands = detect_islands(&g, 4);
        assert_eq!(islands.len(), 1);
        assert_eq!(islands[0].bodies.len(), 4);
        assert_eq!(islands[0].edges.len(), 3);
    }

    #[test]
    fn island_inactive_contacts_ignored() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(2, 3), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.mark_inactive_all();
        // Re-activate only (0,1)
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let islands = detect_islands(&g, 4);
        assert_eq!(islands.len(), 1, "only 1 active island");
    }

    // ── IslandStats ───────────────────────────────────────────────────────────

    #[test]
    fn island_stats_empty() {
        let stats = IslandStats::from_islands(&[]);
        assert_eq!(stats.island_count, 0);
        assert_eq!(stats.max_size, 0);
    }

    #[test]
    fn island_stats_two_islands() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(2, 3), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(3, 4), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let islands = detect_islands(&g, 5);
        let stats = IslandStats::from_islands(&islands);
        assert_eq!(stats.island_count, 2);
        assert_eq!(stats.max_size, 3); // island containing 2,3,4
        assert_eq!(stats.total_bodies, 5);
    }

    // ── propagate_sleep ───────────────────────────────────────────────────────

    #[test]
    fn sleep_propagation_slow_bodies_sleep() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let mut islands = detect_islands(&g, 2);

        // Velocities below threshold
        let vel_sq = [0.0001, 0.0001];
        let mut sleep_counters = [0u32; 2];

        // Run for min_sleep_frames = 3 frames
        let mut to_sleep = Vec::new();
        for _ in 0..3 {
            to_sleep = propagate_sleep(&mut islands, &vel_sq, &mut sleep_counters, 0.01, 3);
        }
        assert!(
            !to_sleep.is_empty(),
            "bodies should be put to sleep after 3 slow frames"
        );
    }

    #[test]
    fn sleep_propagation_fast_body_resets_counter() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let mut islands = detect_islands(&g, 2);

        let mut sleep_counters = [5u32; 2]; // already counted
        // Body 0 is fast
        let vel_sq = [1.0, 0.0001];
        let to_sleep = propagate_sleep(&mut islands, &vel_sq, &mut sleep_counters, 0.01, 3);
        // Island won't sleep because body 0 is fast
        assert!(to_sleep.is_empty(), "fast body should prevent sleep");
        assert_eq!(sleep_counters[0], 0, "counter for fast body should reset");
    }

    // ── ContactFrequency ──────────────────────────────────────────────────────

    #[test]
    fn contact_frequency_classify() {
        assert_eq!(ContactFrequency::classify(0), ContactFrequency::Transient);
        assert_eq!(ContactFrequency::classify(1), ContactFrequency::Transient);
        assert_eq!(ContactFrequency::classify(2), ContactFrequency::ShortLived);
        assert_eq!(ContactFrequency::classify(4), ContactFrequency::ShortLived);
        assert_eq!(ContactFrequency::classify(5), ContactFrequency::Persistent);
        assert_eq!(
            ContactFrequency::classify(100),
            ContactFrequency::Persistent
        );
    }

    #[test]
    fn contact_frequency_counts() {
        let mut g = ContactGraph::new();
        let k0 = ContactKey::new(0, 1);
        let k1 = ContactKey::new(1, 2);
        let k2 = ContactKey::new(2, 3);
        // Insert and age k2 manually by calling update repeatedly
        for _ in 0..6 {
            g.update(k2, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        }
        g.update(k0, [0.0, 1.0, 0.0], 0.1, [0.0; 3]); // age=1 → transient
        g.update(k1, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(k1, [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(k1, [0.0, 1.0, 0.0], 0.1, [0.0; 3]); // age=3 → short_lived

        let (t, s, p) = g.frequency_counts();
        assert_eq!(t, 1, "one transient contact");
        assert_eq!(s, 1, "one short-lived contact");
        assert!(p >= 1, "at least one persistent contact");
    }

    #[test]
    fn persistent_contacts_only_old() {
        let mut g = ContactGraph::new();
        let k_new = ContactKey::new(0, 1);
        let k_old = ContactKey::new(1, 2);
        g.update(k_new, [0.0, 1.0, 0.0], 0.1, [0.0; 3]); // age=1
        for _ in 0..6 {
            g.update(k_old, [0.0, 1.0, 0.0], 0.1, [0.0; 3]); // age=6
        }
        let persistent = g.persistent_contacts();
        assert_eq!(persistent.len(), 1);
        assert_eq!(persistent[0].key, k_old);
    }

    // ── manifold_quality_score ────────────────────────────────────────────────

    #[test]
    fn quality_score_deep_old_warm_is_high() {
        let key = ContactKey::new(0, 1);
        let c = PersistedContact {
            key,
            normal: [0.0, 1.0, 0.0],
            depth: 0.5,
            contact_point: [0.0; 3],
            age: 20,
            cached_impulse: 1.0,
            is_active: true,
        };
        let score = manifold_quality_score(&c);
        assert!(
            score > 0.8,
            "deep, old, warm contact should score > 0.8, got {score}"
        );
        assert!(score <= 1.0);
    }

    #[test]
    fn quality_score_shallow_new_is_low() {
        let key = ContactKey::new(0, 1);
        let c = PersistedContact {
            key,
            normal: [0.0, 1.0, 0.0],
            depth: 0.001,
            contact_point: [0.0; 3],
            age: 0,
            cached_impulse: 0.0,
            is_active: true,
        };
        let score = manifold_quality_score(&c);
        assert!(
            score < 0.2,
            "shallow, new, cold contact should score < 0.2, got {score}"
        );
    }

    #[test]
    fn quality_score_in_unit_range() {
        let key = ContactKey::new(0, 1);
        for depth in [0.0, 0.1, 0.5, 1.0, 10.0] {
            for age in [0u32, 5, 20] {
                let c = PersistedContact {
                    key,
                    normal: [0.0, 1.0, 0.0],
                    depth,
                    contact_point: [0.0; 3],
                    age,
                    cached_impulse: 0.0,
                    is_active: true,
                };
                let score = manifold_quality_score(&c);
                assert!((0.0..=1.0).contains(&score), "score out of range: {score}");
            }
        }
    }

    // ── is_static_dynamic_bipartite ───────────────────────────────────────────

    #[test]
    fn bipartite_check_all_dynamic() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        g.update(ContactKey::new(1, 2), [0.0, 1.0, 0.0], 0.1, [0.0; 3]);
        let is_static = [false, false, false];
        assert!(is_static_dynamic_bipartite(&g, &is_static));
    }

    #[test]
    fn bipartite_check_dynamic_on_static() {
        let mut g = ContactGraph::new();
        // Dynamic body 1 rests on static body 0
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.05, [0.0; 3]);
        let is_static = [true, false];
        assert!(is_static_dynamic_bipartite(&g, &is_static));
    }

    #[test]
    fn bipartite_check_static_static_contact_fails() {
        let mut g = ContactGraph::new();
        g.update(ContactKey::new(0, 1), [0.0, 1.0, 0.0], 0.01, [0.0; 3]);
        let is_static = [true, true]; // both static → should fail
        assert!(!is_static_dynamic_bipartite(&g, &is_static));
    }

    #[test]
    fn bipartite_check_empty_graph_is_bipartite() {
        let g = ContactGraph::new();
        let is_static: Vec<bool> = vec![];
        assert!(is_static_dynamic_bipartite(&g, &is_static));
    }
}
