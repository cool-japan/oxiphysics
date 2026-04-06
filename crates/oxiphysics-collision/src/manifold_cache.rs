#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Persistent contact manifold caching for stable simulation.
//!
//! Provides warm-starting of contact impulses across frames by caching
//! contact points and matching them between simulation steps.
#![allow(dead_code)]

use std::collections::HashMap;

// ─── Math helpers ─────────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn len_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    len_sq3(a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l > 1e-10 {
        scale3(a, 1.0 / l)
    } else {
        [0.0, 0.0, 0.0]
    }
}

/// Apply rotation matrix (row-major 3×3) to a vector.
#[inline]
fn mat3_mul_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Transform a local-space point to world space.
/// `transform` = (translation, rotation_matrix)
#[inline]
fn transform_point(transform: ([f64; 3], [[f64; 3]; 3]), local: [f64; 3]) -> [f64; 3] {
    add3(transform.0, mat3_mul_vec(transform.1, local))
}

// ─── ContactPointId ───────────────────────────────────────────────────────────

/// Utility for computing a hash-based ID for contact point matching.
pub struct ContactPointId;

impl ContactPointId {
    /// Compute a stable hash ID from quantized local positions.
    ///
    /// Positions are quantized to a 1 cm grid before hashing to allow
    /// small numerical drift without changing the contact identity.
    pub fn compute(local_a: [f64; 3], local_b: [f64; 3]) -> u64 {
        const QUANT: f64 = 100.0; // 1 cm grid
        let qa = [
            (local_a[0] * QUANT).round() as i64,
            (local_a[1] * QUANT).round() as i64,
            (local_a[2] * QUANT).round() as i64,
        ];
        let qb = [
            (local_b[0] * QUANT).round() as i64,
            (local_b[1] * QUANT).round() as i64,
            (local_b[2] * QUANT).round() as i64,
        ];
        // FNV-1a style hash over the six integers
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for &v in qa.iter().chain(qb.iter()) {
            for b in v.to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01B3);
            }
        }
        h
    }
}

// ─── ContactPoint ─────────────────────────────────────────────────────────────

/// A single contact point with warm-start impulse data.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// Contact point on body A in world space.
    pub world_pos_a: [f64; 3],
    /// Contact point on body B in world space.
    pub world_pos_b: [f64; 3],
    /// Contact point in body A local space.
    pub local_pos_a: [f64; 3],
    /// Contact point in body B local space.
    pub local_pos_b: [f64; 3],
    /// Contact normal pointing from B to A.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Accumulated tangential impulse for warm starting (2 tangent directions).
    pub tangent_impulse: [f64; 2],
    /// Accumulated normal impulse for warm starting.
    pub normal_impulse: f64,
    /// Number of frames this contact has persisted.
    pub lifetime: u32,
    /// Hash-based ID for matching across frames.
    pub id: u64,
}

impl ContactPoint {
    /// Create a new contact point, computing its ID automatically.
    pub fn new(
        world_pos_a: [f64; 3],
        world_pos_b: [f64; 3],
        local_pos_a: [f64; 3],
        local_pos_b: [f64; 3],
        normal: [f64; 3],
        depth: f64,
    ) -> Self {
        let id = ContactPointId::compute(local_pos_a, local_pos_b);
        Self {
            world_pos_a,
            world_pos_b,
            local_pos_a,
            local_pos_b,
            normal,
            depth,
            tangent_impulse: [0.0, 0.0],
            normal_impulse: 0.0,
            lifetime: 0,
            id,
        }
    }
}

// ─── PersistentManifold ───────────────────────────────────────────────────────

/// Persistent contact manifold holding up to 4 contact points with warm-start data.
///
/// Named `PersistentManifold` to distinguish from the simpler `ContactManifold`
/// used in narrow-phase output.
#[derive(Debug, Clone)]
pub struct PersistentManifold {
    /// Handle ID of body A.
    pub body_a: u32,
    /// Handle ID of body B.
    pub body_b: u32,
    /// Contact points (max 4).
    pub points: Vec<ContactPoint>,
    /// Average contact normal.
    pub normal: [f64; 3],
    /// Friction coefficient.
    pub friction: f64,
    /// Restitution coefficient.
    pub restitution: f64,
    /// Whether this manifold was touched this frame.
    pub is_active: bool,
}

/// Distance threshold (squared) below which two contact points are considered
/// the same contact.
const MATCH_DIST_SQ: f64 = 0.01 * 0.01; // 1 cm

/// Separation threshold beyond which a cached contact is considered broken.
const SEPARATION_THRESHOLD: f64 = 0.02; // 2 cm

impl PersistentManifold {
    /// Create a new empty persistent manifold.
    pub fn new(body_a: u32, body_b: u32, friction: f64, restitution: f64) -> Self {
        Self {
            body_a,
            body_b,
            points: Vec::with_capacity(4),
            normal: [0.0, 1.0, 0.0],
            friction,
            restitution,
            is_active: true,
        }
    }

    /// Add a new contact point or update an existing cached one.
    ///
    /// If a point with the same ID (or within the distance threshold) already
    /// exists, the cached warm-start impulses are transferred to the new point.
    /// The manifold is pruned to at most 4 points using area-maximising selection.
    pub fn add_or_update(&mut self, mut new_point: ContactPoint) {
        // Try to find a matching existing point by ID or proximity.
        let match_idx = self.points.iter().position(|p| {
            if p.id == new_point.id {
                return true;
            }
            let d = len_sq3(sub3(p.local_pos_a, new_point.local_pos_a));
            d < MATCH_DIST_SQ
        });

        if let Some(idx) = match_idx {
            // Transfer warm-start impulses.
            new_point.normal_impulse = self.points[idx].normal_impulse;
            new_point.tangent_impulse = self.points[idx].tangent_impulse;
            new_point.lifetime = self.points[idx].lifetime + 1;
            self.points[idx] = new_point;
        } else {
            self.points.push(new_point);
            if self.points.len() > 4 {
                self.points = ManifoldReduction::reduce_to_4_points(&self.points);
            }
        }

        self.recompute_normal();
    }

    /// Remove stale contact points based on world-space re-projection.
    ///
    /// For each cached point the local positions are transformed back to world
    /// space using the current body transforms. The point is removed if the
    /// bodies have separated or if the projected position has drifted too far
    /// from the cached world position.
    pub fn remove_stale(
        &mut self,
        transform_a: ([f64; 3], [[f64; 3]; 3]),
        transform_b: ([f64; 3], [[f64; 3]; 3]),
    ) {
        self.points.retain(|p| {
            let wa = transform_point(transform_a, p.local_pos_a);
            let wb = transform_point(transform_b, p.local_pos_b);

            // Check separation along normal.
            let sep = dot3(sub3(wa, wb), p.normal);
            if sep > SEPARATION_THRESHOLD {
                return false;
            }

            // Check positional drift of world_pos_a.
            if len_sq3(sub3(wa, p.world_pos_a)) > SEPARATION_THRESHOLD * SEPARATION_THRESHOLD {
                return false;
            }

            // Check positional drift of world_pos_b.
            if len_sq3(sub3(wb, p.world_pos_b)) > SEPARATION_THRESHOLD * SEPARATION_THRESHOLD {
                return false;
            }

            true
        });

        self.recompute_normal();
    }

    fn recompute_normal(&mut self) {
        if self.points.is_empty() {
            return;
        }
        let mut avg = [0.0f64; 3];
        for p in &self.points {
            avg = add3(avg, p.normal);
        }
        let n = self.points.len() as f64;
        self.normal = normalize3(scale3(avg, 1.0 / n));
    }
}

// ─── ManifoldReduction ────────────────────────────────────────────────────────

/// Algorithms for reducing a contact set to at most 4 representative points.
pub struct ManifoldReduction;

impl ManifoldReduction {
    /// Reduce a slice of contact points to at most 4 points.
    ///
    /// Selection strategy:
    /// 1. Keep the deepest penetrating point.
    /// 2. Keep the point furthest from point 1.
    /// 3. Keep the point that maximises triangle area with points 1 and 2.
    /// 4. Keep the point that maximises quadrilateral area with points 1–3.
    pub fn reduce_to_4_points(points: &[ContactPoint]) -> Vec<ContactPoint> {
        if points.len() <= 4 {
            return points.to_vec();
        }

        // 1. Deepest point.
        let idx0 = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.depth
                    .partial_cmp(&b.depth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // 2. Furthest from point 0.
        let p0 = points[idx0].world_pos_a;
        let idx1 = points
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx0)
            .max_by(|(_, a), (_, b)| {
                len_sq3(sub3(a.world_pos_a, p0))
                    .partial_cmp(&len_sq3(sub3(b.world_pos_a, p0)))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // 3. Maximise triangle area with p0 and p1.
        let p1 = points[idx1].world_pos_a;
        let idx2 = points
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx0 && *i != idx1)
            .max_by(|(_, a), (_, b)| {
                Self::tri_area_sq(p0, p1, a.world_pos_a)
                    .partial_cmp(&Self::tri_area_sq(p0, p1, b.world_pos_a))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // 4. Maximise quad area.
        let p2 = points[idx2].world_pos_a;
        let maybe_idx3 = points
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx0 && *i != idx1 && *i != idx2)
            .max_by(|(_, a), (_, b)| {
                Self::quad_area_sq(p0, p1, p2, a.world_pos_a)
                    .partial_cmp(&Self::quad_area_sq(p0, p1, p2, b.world_pos_a))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i);

        let mut result = vec![
            points[idx0].clone(),
            points[idx1].clone(),
            points[idx2].clone(),
        ];
        if let Some(idx3) = maybe_idx3 {
            result.push(points[idx3].clone());
        }
        result
    }

    /// Squared magnitude of cross product, proportional to triangle area squared.
    fn tri_area_sq(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        len_sq3(cross3(sub3(b, a), sub3(c, a)))
    }

    /// Approximate quad area: sum of two triangle areas (squared, for comparison only).
    fn quad_area_sq(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
        Self::tri_area_sq(a, b, c) + Self::tri_area_sq(a, c, d)
    }

    /// Compute the approximate area of the convex hull of the contact points
    /// projected onto the contact plane defined by the first point's normal.
    ///
    /// Uses the shoelace formula on the projected 2-D coordinates.
    pub fn contact_area(points: &[ContactPoint]) -> f64 {
        if points.len() < 3 {
            return 0.0;
        }

        let n = points[0].normal;
        let t1 = Self::make_tangent(n);
        let t2 = cross3(n, t1);

        let origin = points[0].world_pos_a;
        let projected: Vec<[f64; 2]> = points
            .iter()
            .map(|p| {
                let d = sub3(p.world_pos_a, origin);
                [dot3(d, t1), dot3(d, t2)]
            })
            .collect();

        shoelace_area(&projected)
    }

    fn make_tangent(n: [f64; 3]) -> [f64; 3] {
        let candidate = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        normalize3(cross3(n, candidate))
    }
}

// ─── Shoelace helper ──────────────────────────────────────────────────────────

fn shoelace_area(pts: &[[f64; 2]]) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pts[i][0] * pts[j][1];
        area -= pts[j][0] * pts[i][1];
    }
    area.abs() * 0.5
}

// ─── ContactManifoldCache ─────────────────────────────────────────────────────

/// World-level cache of persistent contact manifolds.
///
/// Stores one `PersistentManifold` per body pair, keyed by ordered `(min, max)`
/// body-handle pair.
pub struct ContactManifoldCache {
    /// All active manifolds, keyed by ordered body-handle pair.
    pub manifolds: HashMap<(u32, u32), PersistentManifold>,
    /// Maximum number of frames a manifold can persist without being touched.
    pub max_lifetime: u32,
}

impl ContactManifoldCache {
    /// Create a new cache.
    pub fn new(max_lifetime: u32) -> Self {
        Self {
            manifolds: HashMap::new(),
            max_lifetime,
        }
    }

    /// Return the canonical key for a body pair (always min, max).
    #[inline]
    fn key(body_a: u32, body_b: u32) -> (u32, u32) {
        (body_a.min(body_b), body_a.max(body_b))
    }

    /// Get or create the persistent manifold for the given body pair.
    pub fn get_or_create(
        &mut self,
        body_a: u32,
        body_b: u32,
        friction: f64,
        restitution: f64,
    ) -> &mut PersistentManifold {
        let k = Self::key(body_a, body_b);
        self.manifolds
            .entry(k)
            .or_insert_with(|| PersistentManifold::new(body_a, body_b, friction, restitution))
    }

    /// Merge a set of newly detected contact points into the cached manifold.
    ///
    /// Each new point is added-or-updated; existing point lifetimes are
    /// incremented automatically inside `add_or_update`.
    pub fn update_manifold(&mut self, body_a: u32, body_b: u32, new_points: Vec<ContactPoint>) {
        let k = Self::key(body_a, body_b);
        if let Some(manifold) = self.manifolds.get_mut(&k) {
            manifold.is_active = true;
            for pt in new_points {
                manifold.add_or_update(pt);
            }
            for p in &mut manifold.points {
                p.lifetime = p.lifetime.saturating_add(1);
            }
        }
    }

    /// Remove manifolds that are no longer active or have no points.
    pub fn remove_inactive(&mut self) {
        self.manifolds
            .retain(|_, m| m.is_active && !m.points.is_empty());
    }

    /// Mark all manifolds as inactive at the start of a frame.
    ///
    /// Manifolds that are still colliding will be reactivated when
    /// `update_manifold` is called.
    pub fn begin_frame(&mut self) {
        for m in self.manifolds.values_mut() {
            m.is_active = false;
        }
    }
}

// ─── ManifoldPointMatcher ────────────────────────────────────────────────────

/// Utility for matching new contact points to existing cached points.
pub struct ManifoldPointMatcher;

impl ManifoldPointMatcher {
    /// Find the best matching point in `existing` for `new_point`.
    ///
    /// Matching is done by:
    /// 1. Exact ID match (hash equality).
    /// 2. Proximity in local-space A (distance² < threshold).
    ///
    /// Returns the index into `existing` or `None` if no match.
    pub fn find_match(
        existing: &[ContactPoint],
        new_point: &ContactPoint,
        dist_sq_threshold: f64,
    ) -> Option<usize> {
        for (i, p) in existing.iter().enumerate() {
            if p.id == new_point.id {
                return Some(i);
            }
            let d = len_sq3(sub3(p.local_pos_a, new_point.local_pos_a));
            if d < dist_sq_threshold {
                return Some(i);
            }
        }
        None
    }

    /// Match all points in `new_points` against `existing`, returning
    /// pairs `(new_idx, existing_idx)` for each successful match.
    pub fn match_all(
        existing: &[ContactPoint],
        new_points: &[ContactPoint],
        dist_sq_threshold: f64,
    ) -> Vec<(usize, usize)> {
        let mut matches = Vec::new();
        for (ni, np) in new_points.iter().enumerate() {
            if let Some(ei) = Self::find_match(existing, np, dist_sq_threshold) {
                matches.push((ni, ei));
            }
        }
        matches
    }
}

// ─── ManifoldLifetimeManager ─────────────────────────────────────────────────

/// Manages the lifetime of manifolds in a cache: marks active ones,
/// ages inactive ones, and removes stale ones.
pub struct ManifoldLifetimeManager {
    /// Number of frames a manifold can be inactive before removal.
    pub max_inactive_frames: u32,
}

impl ManifoldLifetimeManager {
    /// Create a new manager.
    pub fn new(max_inactive_frames: u32) -> Self {
        Self {
            max_inactive_frames,
        }
    }

    /// Increment lifetime counters for active manifold points; drop ones that
    /// have exceeded `max_lifetime`.
    pub fn age_manifold(&self, manifold: &mut PersistentManifold) {
        manifold
            .points
            .retain(|p| p.lifetime <= self.max_inactive_frames);
    }

    /// Process all manifolds in a cache: age each and remove empty ones.
    pub fn process_cache(&self, cache: &mut ContactManifoldCache) {
        for manifold in cache.manifolds.values_mut() {
            self.age_manifold(manifold);
        }
        cache.manifolds.retain(|_, m| !m.points.is_empty());
    }

    /// Return how many manifolds are active.
    pub fn active_count(&self, cache: &ContactManifoldCache) -> usize {
        cache.manifolds.values().filter(|m| m.is_active).count()
    }
}

// ─── ManifoldCompressor ───────────────────────────────────────────────────────

/// Compresses a manifold by removing redundant contact points while
/// preserving the most physically significant ones.
pub struct ManifoldCompressor;

impl ManifoldCompressor {
    /// Compress `points` to at most `max_count` points.
    ///
    /// Keeps the point with the greatest penetration depth, then selects
    /// remaining points to maximise the covered area (same strategy as
    /// `ManifoldReduction::reduce_to_4_points`).
    pub fn compress(points: &[ContactPoint], max_count: usize) -> Vec<ContactPoint> {
        if points.len() <= max_count {
            return points.to_vec();
        }
        if max_count == 0 {
            return vec![];
        }
        if max_count == 1 {
            // Just the deepest point
            let idx = points
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| {
                    a.depth
                        .partial_cmp(&b.depth)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            return vec![points[idx].clone()];
        }
        // For max_count >= 4 fall through to existing reduction
        ManifoldReduction::reduce_to_4_points(points)
    }

    /// Remove duplicate contact points (same ID).
    pub fn deduplicate(points: &mut Vec<ContactPoint>) {
        let mut seen = std::collections::HashSet::new();
        points.retain(|p| seen.insert(p.id));
    }

    /// Merge two contact sets, preferring warm-start data from the existing set.
    pub fn merge(
        existing: &[ContactPoint],
        incoming: &[ContactPoint],
        max_count: usize,
    ) -> Vec<ContactPoint> {
        let mut merged: Vec<ContactPoint> = existing.to_vec();
        for inc in incoming {
            if let Some(idx) = ManifoldPointMatcher::find_match(&merged, inc, MATCH_DIST_SQ) {
                // Transfer warm-start data
                let mut updated = inc.clone();
                updated.normal_impulse = merged[idx].normal_impulse;
                updated.tangent_impulse = merged[idx].tangent_impulse;
                updated.lifetime = merged[idx].lifetime + 1;
                merged[idx] = updated;
            } else {
                merged.push(inc.clone());
            }
        }
        Self::compress(&merged, max_count)
    }
}

// ─── PersistentManifold extensions ───────────────────────────────────────────

impl PersistentManifold {
    /// Update the manifold with a completely new set of contact points.
    ///
    /// Points are matched against the existing cache to preserve warm-start
    /// impulses, then compressed to at most 4.
    pub fn update_from_new_contacts(&mut self, new_contacts: Vec<ContactPoint>) {
        let merged = ManifoldCompressor::merge(&self.points, &new_contacts, 4);
        self.points = merged;
        self.recompute_normal();
        self.is_active = true;
    }

    /// Extract warm-start data for the solver: returns (normal_impulse, tangent_impulse) per point.
    pub fn warmstart_data(&self) -> Vec<(f64, [f64; 2])> {
        self.points
            .iter()
            .map(|p| (p.normal_impulse, p.tangent_impulse))
            .collect()
    }

    /// Apply solver results back to the warm-start cache.
    pub fn store_solver_impulses(&mut self, impulses: &[(f64, [f64; 2])]) {
        for (p, &(ni, ti)) in self.points.iter_mut().zip(impulses.iter()) {
            p.normal_impulse = ni;
            p.tangent_impulse = ti;
        }
    }

    /// Number of active contact points.
    pub fn contact_count(&self) -> usize {
        self.points.len()
    }

    /// Maximum penetration depth across all contact points.
    pub fn max_depth(&self) -> f64 {
        self.points.iter().map(|p| p.depth).fold(0.0f64, f64::max)
    }

    /// Clamp all accumulated impulses to non-negative (velocity-based solvers).
    pub fn clamp_impulses(&mut self) {
        for p in &mut self.points {
            p.normal_impulse = p.normal_impulse.max(0.0);
        }
    }
}

// ─── ContactManifoldCache extensions ─────────────────────────────────────────

impl ContactManifoldCache {
    /// Total number of contact points across all active manifolds.
    pub fn total_contact_points(&self) -> usize {
        self.manifolds.values().map(|m| m.points.len()).sum()
    }

    /// Number of manifolds in the cache.
    pub fn manifold_count(&self) -> usize {
        self.manifolds.len()
    }

    /// Clear all manifolds.
    pub fn clear(&mut self) {
        self.manifolds.clear();
    }

    /// Collect all warm-start data for a given body pair (if present).
    pub fn get_warmstart(&self, body_a: u32, body_b: u32) -> Option<Vec<(f64, [f64; 2])>> {
        let k = Self::key(body_a, body_b);
        self.manifolds.get(&k).map(|m| m.warmstart_data())
    }

    /// Update all manifolds that belong to a single body (e.g., on body removal).
    pub fn remove_body(&mut self, body_id: u32) {
        self.manifolds
            .retain(|_, m| m.body_a != body_id && m.body_b != body_id);
    }
}

// ─── Contact caching strategies ──────────────────────────────────────────────

/// Strategy for matching contact points across frames.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CachingStrategy {
    /// Match by exact hash ID only.
    IdOnly,
    /// Match by spatial proximity only.
    ProximityOnly,
    /// Match by ID first, then fall back to proximity.
    IdThenProximity,
}

/// Apply a caching strategy to find the best match for a new contact in existing contacts.
pub fn find_match_with_strategy(
    existing: &[ContactPoint],
    new_point: &ContactPoint,
    strategy: CachingStrategy,
    dist_sq_threshold: f64,
) -> Option<usize> {
    match strategy {
        CachingStrategy::IdOnly => existing.iter().position(|p| p.id == new_point.id),
        CachingStrategy::ProximityOnly => existing.iter().enumerate().find_map(|(i, p)| {
            if len_sq3(sub3(p.local_pos_a, new_point.local_pos_a)) < dist_sq_threshold {
                Some(i)
            } else {
                None
            }
        }),
        CachingStrategy::IdThenProximity => {
            ManifoldPointMatcher::find_match(existing, new_point, dist_sq_threshold)
        }
    }
}

// ─── Warm-start data aging ─────────────────────────────────────────────────────

/// Scale warm-start impulses by an age factor.
///
/// Reduces the effective warm-start contribution of older contacts.
/// `age_factor` should be in `(0, 1]` — e.g. `0.95` per frame.
pub fn age_warm_start(point: &mut ContactPoint, age_factor: f64) {
    point.normal_impulse *= age_factor;
    point.tangent_impulse[0] *= age_factor;
    point.tangent_impulse[1] *= age_factor;
}

/// Apply aging to all contact points in a persistent manifold.
pub fn age_manifold_warm_start(manifold: &mut PersistentManifold, age_factor: f64) {
    for p in &mut manifold.points {
        age_warm_start(p, age_factor);
    }
}

/// Apply aging to all manifolds in a cache.
pub fn age_cache_warm_start(cache: &mut ContactManifoldCache, age_factor: f64) {
    for manifold in cache.manifolds.values_mut() {
        age_manifold_warm_start(manifold, age_factor);
    }
}

// ─── Position correction impulses ────────────────────────────────────────────

/// Compute a Baumgarte position correction impulse magnitude.
///
/// `depth` is penetration depth, `beta` is the correction factor (0.1–0.3 typical),
/// `dt` is the time step.  Returns the correction impulse to apply along the normal.
pub fn baumgarte_correction(depth: f64, beta: f64, dt: f64) -> f64 {
    if dt > 1e-12 {
        (beta * depth / dt).max(0.0)
    } else {
        0.0
    }
}

/// Compute a slop-clamped Baumgarte correction.
///
/// `slop` is a small penetration allowance (e.g. 0.005 m) that is not corrected.
pub fn baumgarte_correction_slop(depth: f64, beta: f64, dt: f64, slop: f64) -> f64 {
    baumgarte_correction((depth - slop).max(0.0), beta, dt)
}

/// Apply position correction impulses to all contact points in a manifold.
///
/// Returns the per-point correction magnitudes.
pub fn apply_position_corrections(
    manifold: &PersistentManifold,
    beta: f64,
    dt: f64,
    slop: f64,
) -> Vec<f64> {
    manifold
        .points
        .iter()
        .map(|p| baumgarte_correction_slop(p.depth, beta, dt, slop))
        .collect()
}

// ─── Island-level manifold batching ──────────────────────────────────────────

/// A simulation island: a group of bodies connected by contacts.
#[derive(Debug, Clone)]
pub struct ContactIsland {
    /// Body IDs in this island.
    pub bodies: Vec<u32>,
    /// Manifold keys `(min_id, max_id)` belonging to this island.
    pub manifold_keys: Vec<(u32, u32)>,
}

impl ContactIsland {
    /// Create an empty island.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            manifold_keys: Vec::new(),
        }
    }

    /// Number of bodies in the island.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Number of contacts in the island.
    pub fn contact_count(&self) -> usize {
        self.manifold_keys.len()
    }
}

impl Default for ContactIsland {
    fn default() -> Self {
        Self::new()
    }
}

/// Build contact islands from a manifold cache using union-find.
///
/// Returns a list of islands, each containing body IDs and manifold keys.
pub fn build_contact_islands(cache: &ContactManifoldCache) -> Vec<ContactIsland> {
    // Collect all body IDs
    let mut all_bodies: Vec<u32> = Vec::new();
    for (a, b) in cache.manifolds.keys() {
        if !all_bodies.contains(a) {
            all_bodies.push(*a);
        }
        if !all_bodies.contains(b) {
            all_bodies.push(*b);
        }
    }

    let n = all_bodies.len();
    if n == 0 {
        return Vec::new();
    }

    // Union-Find
    let mut parent: Vec<usize> = (0..n).collect();

    let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path compression
            x = parent[x];
        }
        x
    };

    for (a, b) in cache.manifolds.keys() {
        let Some(ia) = all_bodies.iter().position(|&id| id == *a) else {
            continue;
        };
        let Some(ib) = all_bodies.iter().position(|&id| id == *b) else {
            continue;
        };
        let ra = find(&mut parent, ia);
        let rb = find(&mut parent, ib);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    // Group bodies by root
    let mut island_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut islands: Vec<ContactIsland> = Vec::new();

    for i in 0..n {
        let root = find(&mut parent, i);
        let island_idx = *island_map.entry(root).or_insert_with(|| {
            islands.push(ContactIsland::new());
            islands.len() - 1
        });
        islands[island_idx].bodies.push(all_bodies[i]);
    }

    // Assign manifold keys to islands
    for key in cache.manifolds.keys() {
        let Some(ia) = all_bodies.iter().position(|&id| id == key.0) else {
            continue;
        };
        let root = find(&mut parent, ia);
        let island_idx = island_map[&root];
        islands[island_idx].manifold_keys.push(*key);
    }

    islands
}

// ─── Manifold quality metrics ─────────────────────────────────────────────────

/// Compute quality metrics for a manifold.
#[derive(Debug, Clone)]
pub struct ManifoldMetrics {
    /// Number of contact points.
    pub contact_count: usize,
    /// Maximum penetration depth.
    pub max_depth: f64,
    /// Average penetration depth.
    pub avg_depth: f64,
    /// Maximum pairwise distance between contact points (contact area measure).
    pub spread: f64,
    /// Whether all contact points have valid warm-start data.
    pub is_warm: bool,
}

/// Compute quality metrics for a persistent manifold.
pub fn compute_manifold_metrics(manifold: &PersistentManifold) -> ManifoldMetrics {
    let n = manifold.points.len();
    if n == 0 {
        return ManifoldMetrics {
            contact_count: 0,
            max_depth: 0.0,
            avg_depth: 0.0,
            spread: 0.0,
            is_warm: false,
        };
    }

    let max_depth = manifold
        .points
        .iter()
        .map(|p| p.depth)
        .fold(0.0f64, f64::max);
    let avg_depth = manifold.points.iter().map(|p| p.depth).sum::<f64>() / n as f64;
    let is_warm = manifold
        .points
        .iter()
        .all(|p| p.normal_impulse.abs() > 0.0 || p.lifetime > 0);

    // Max pairwise distance
    let mut spread = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = len3(sub3(
                manifold.points[i].world_pos_a,
                manifold.points[j].world_pos_a,
            ));
            if d > spread {
                spread = d;
            }
        }
    }

    ManifoldMetrics {
        contact_count: n,
        max_depth,
        avg_depth,
        spread,
        is_warm,
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_contact(
        pos_a: [f64; 3],
        pos_b: [f64; 3],
        depth: f64,
        normal_impulse: f64,
    ) -> ContactPoint {
        let mut cp = ContactPoint::new(pos_a, pos_b, pos_a, pos_b, [0.0, 1.0, 0.0], depth);
        cp.normal_impulse = normal_impulse;
        cp
    }

    // ------------------------------------------------------------------
    // PersistentManifold::add_or_update
    // ------------------------------------------------------------------

    #[test]
    fn test_add_same_point_twice_no_overflow() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 1.0);
        m.add_or_update(pt.clone());
        m.add_or_update(pt.clone());
        assert_eq!(m.points.len(), 1, "same point added twice should stay as 1");
    }

    #[test]
    fn test_add_or_update_max_4_points() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let positions: [[f64; 3]; 6] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [0.5, 0.0, 0.5],
            [-0.5, 0.0, -0.5],
        ];
        for &p in &positions {
            let neg_p = scale3(p, -1.0);
            let cp = make_contact(p, neg_p, 0.01, 0.0);
            m.add_or_update(cp);
        }
        assert!(
            m.points.len() <= 4,
            "manifold must not exceed 4 points, got {}",
            m.points.len()
        );
    }

    // ------------------------------------------------------------------
    // ManifoldReduction
    // ------------------------------------------------------------------

    #[test]
    fn test_reduce_to_4_points_with_6_inputs() {
        let positions: [[f64; 3]; 6] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [0.5, 0.0, 0.5],
            [-0.5, 0.0, -0.5],
        ];
        let points: Vec<ContactPoint> = positions
            .iter()
            .map(|&p| ContactPoint::new(p, p, p, p, [0.0, 1.0, 0.0], 0.01))
            .collect();
        let reduced = ManifoldReduction::reduce_to_4_points(&points);
        assert!(
            reduced.len() <= 4,
            "reduce_to_4_points must return <=4 points, got {}",
            reduced.len()
        );
    }

    // ------------------------------------------------------------------
    // ContactManifoldCache::get_or_create
    // ------------------------------------------------------------------

    #[test]
    fn test_get_or_create_new_pair() {
        let mut cache = ContactManifoldCache::new(5);
        let m = cache.get_or_create(1, 2, 0.5, 0.3);
        assert_eq!(m.friction, 0.5);
        assert_eq!(m.restitution, 0.3);
        assert!(m.is_active);
    }

    #[test]
    fn test_get_or_create_ordered_key() {
        let mut cache = ContactManifoldCache::new(5);
        let _ = cache.get_or_create(2, 1, 0.4, 0.2);
        assert!(
            cache.manifolds.contains_key(&(1, 2)),
            "key should be canonicalized to (1,2)"
        );
    }

    // ------------------------------------------------------------------
    // ContactManifoldCache::begin_frame
    // ------------------------------------------------------------------

    #[test]
    fn test_begin_frame_marks_all_inactive() {
        let mut cache = ContactManifoldCache::new(5);
        let _ = cache.get_or_create(0, 1, 0.5, 0.3);
        let _ = cache.get_or_create(2, 3, 0.5, 0.3);
        cache.begin_frame();
        for m in cache.manifolds.values() {
            assert!(
                !m.is_active,
                "manifold should be inactive after begin_frame"
            );
        }
    }

    // ------------------------------------------------------------------
    // Warm start: normal_impulse preserved on update
    // ------------------------------------------------------------------

    #[test]
    fn test_warm_start_impulse_preserved() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        pt.normal_impulse = 42.0;
        m.add_or_update(pt);

        // Second frame: new point at the same location.
        let pt2 = ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        m.add_or_update(pt2);

        assert_eq!(
            m.points[0].normal_impulse, 42.0,
            "normal_impulse should be warm-started from previous frame"
        );
    }

    // ------------------------------------------------------------------
    // contact_area > 0 for 3+ non-collinear points
    // ------------------------------------------------------------------

    #[test]
    fn test_contact_area_nonzero_for_triangle() {
        let pts = vec![
            ContactPoint::new(
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                0.01,
            ),
            ContactPoint::new(
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                0.01,
            ),
            ContactPoint::new(
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                0.01,
            ),
        ];
        let area = ManifoldReduction::contact_area(&pts);
        assert!(
            area > 0.0,
            "contact_area of non-collinear triangle must be > 0, got {}",
            area
        );
    }

    // ------------------------------------------------------------------
    // ManifoldPointMatcher
    // ------------------------------------------------------------------

    #[test]
    fn test_matcher_finds_by_id() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 5.0)];
        // Build a point with the same local positions → same ID
        let new_pt = ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let idx = ManifoldPointMatcher::find_match(&existing, &new_pt, MATCH_DIST_SQ);
        assert_eq!(idx, Some(0), "should match by ID");
    }

    #[test]
    fn test_matcher_finds_by_proximity() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 1.0)];
        // Slightly different position, but within MATCH_DIST_SQ
        let new_pt = ContactPoint::new(
            [0.001, 0.0, 0.0],
            [0.001, -0.01, 0.0],
            [0.001, 0.0, 0.0], // close to [0,0,0]
            [0.001, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let idx = ManifoldPointMatcher::find_match(&existing, &new_pt, MATCH_DIST_SQ);
        assert!(idx.is_some(), "should match by proximity");
    }

    #[test]
    fn test_matcher_no_match_far() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 1.0)];
        let new_pt = ContactPoint::new(
            [5.0, 0.0, 0.0],
            [5.0, -0.01, 0.0],
            [5.0, 0.0, 0.0],
            [5.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let idx = ManifoldPointMatcher::find_match(&existing, &new_pt, MATCH_DIST_SQ);
        assert!(idx.is_none(), "far point should not match");
    }

    // ------------------------------------------------------------------
    // ManifoldLifetimeManager
    // ------------------------------------------------------------------

    #[test]
    fn test_lifetime_manager_ages_out_old_points() {
        let mgr = ManifoldLifetimeManager::new(2);
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 1.0);
        pt.lifetime = 100; // Very old
        m.points.push(pt);
        mgr.age_manifold(&mut m);
        assert!(m.points.is_empty(), "old point should be removed");
    }

    #[test]
    fn test_lifetime_manager_keeps_young_points() {
        let mgr = ManifoldLifetimeManager::new(10);
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 1.0);
        pt.lifetime = 3;
        m.points.push(pt);
        mgr.age_manifold(&mut m);
        assert_eq!(m.points.len(), 1, "young point should be kept");
    }

    // ------------------------------------------------------------------
    // ManifoldCompressor
    // ------------------------------------------------------------------

    #[test]
    fn test_compressor_no_change_small_set() {
        let pts: Vec<ContactPoint> = (0..3)
            .map(|i| make_contact([i as f64, 0.0, 0.0], [i as f64, -0.01, 0.0], 0.01, 0.0))
            .collect();
        let out = ManifoldCompressor::compress(&pts, 4);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn test_compressor_single_point() {
        let pts: Vec<ContactPoint> = vec![
            make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.05, 1.0),
            make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.01, 0.0),
        ];
        let out = ManifoldCompressor::compress(&pts, 1);
        assert_eq!(out.len(), 1);
        // Should keep the deepest
        assert!((out[0].depth - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_compressor_merge_preserves_warmstart() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 7.0)];
        let incoming = vec![ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        )];
        let merged = ManifoldCompressor::merge(&existing, &incoming, 4);
        assert_eq!(merged.len(), 1);
        assert!(
            (merged[0].normal_impulse - 7.0).abs() < 1e-10,
            "warm-start should be 7.0"
        );
    }

    // ------------------------------------------------------------------
    // PersistentManifold extensions
    // ------------------------------------------------------------------

    #[test]
    fn test_update_from_new_contacts() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let contacts = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0)];
        m.update_from_new_contacts(contacts);
        assert_eq!(m.contact_count(), 1);
        assert!(m.is_active);
    }

    #[test]
    fn test_max_depth() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        m.add_or_update(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.02, 0.0));
        m.add_or_update(make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.05, 0.0));
        assert!((m.max_depth() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_warmstart_data_roundtrip() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0);
        pt.normal_impulse = 3.0;
        pt.tangent_impulse = [1.0, 2.0];
        m.points.push(pt);
        let ws = m.warmstart_data();
        assert_eq!(ws.len(), 1);
        assert!((ws[0].0 - 3.0).abs() < 1e-12);
        // Store different values
        m.store_solver_impulses(&[(10.0, [5.0, 6.0])]);
        assert!((m.points[0].normal_impulse - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_clamp_impulses() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0);
        pt.normal_impulse = -3.0;
        m.points.push(pt);
        m.clamp_impulses();
        assert_eq!(m.points[0].normal_impulse, 0.0);
    }

    // ------------------------------------------------------------------
    // ContactManifoldCache extensions
    // ------------------------------------------------------------------

    #[test]
    fn test_cache_total_contact_points() {
        let mut cache = ContactManifoldCache::new(5);
        let m = cache.get_or_create(0, 1, 0.5, 0.3);
        m.add_or_update(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0));
        m.add_or_update(make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.01, 0.0));
        assert_eq!(cache.total_contact_points(), 2);
    }

    #[test]
    fn test_cache_remove_body() {
        let mut cache = ContactManifoldCache::new(5);
        let _ = cache.get_or_create(0, 1, 0.5, 0.3);
        let _ = cache.get_or_create(2, 3, 0.5, 0.3);
        cache.remove_body(0);
        assert_eq!(cache.manifold_count(), 1);
    }

    #[test]
    fn test_cache_clear() {
        let mut cache = ContactManifoldCache::new(5);
        let _ = cache.get_or_create(0, 1, 0.5, 0.3);
        cache.clear();
        assert_eq!(cache.manifold_count(), 0);
    }

    #[test]
    fn test_cache_get_warmstart_none() {
        let cache = ContactManifoldCache::new(5);
        assert!(cache.get_warmstart(0, 1).is_none());
    }

    // ─── CachingStrategy ──────────────────────────────────────────────────────

    #[test]
    fn test_caching_strategy_id_only_finds_match() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0)];
        let new_pt = ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let result =
            find_match_with_strategy(&existing, &new_pt, CachingStrategy::IdOnly, MATCH_DIST_SQ);
        assert_eq!(result, Some(0), "ID-only strategy should match by hash ID");
    }

    #[test]
    fn test_caching_strategy_proximity_only_finds_match() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0)];
        // Slightly different position but within threshold
        let new_pt = ContactPoint::new(
            [0.001, 0.0, 0.0],
            [0.001, -0.01, 0.0],
            [0.001, 0.0, 0.0],
            [0.001, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let result = find_match_with_strategy(
            &existing,
            &new_pt,
            CachingStrategy::ProximityOnly,
            MATCH_DIST_SQ,
        );
        assert!(
            result.is_some(),
            "Proximity strategy should find nearby match"
        );
    }

    #[test]
    fn test_caching_strategy_proximity_only_no_match_far() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0)];
        let new_pt = ContactPoint::new(
            [10.0, 0.0, 0.0],
            [10.0, -0.01, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let result = find_match_with_strategy(
            &existing,
            &new_pt,
            CachingStrategy::ProximityOnly,
            MATCH_DIST_SQ,
        );
        assert!(
            result.is_none(),
            "Far point should not match with proximity strategy"
        );
    }

    #[test]
    fn test_caching_strategy_id_then_proximity_finds_by_id() {
        let existing = vec![make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 3.0)];
        let new_pt = ContactPoint::new(
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -0.01, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        );
        let result = find_match_with_strategy(
            &existing,
            &new_pt,
            CachingStrategy::IdThenProximity,
            MATCH_DIST_SQ,
        );
        assert_eq!(result, Some(0), "IdThenProximity should find by ID");
    }

    // ─── age_warm_start ────────────────────────────────────────────────────────

    #[test]
    fn test_age_warm_start_scales_impulses() {
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 10.0);
        pt.tangent_impulse = [4.0, 2.0];
        age_warm_start(&mut pt, 0.9);
        assert!(
            (pt.normal_impulse - 9.0).abs() < 1e-10,
            "normal_impulse should be 9.0"
        );
        assert!(
            (pt.tangent_impulse[0] - 3.6).abs() < 1e-10,
            "tangent[0] should be 3.6"
        );
        assert!(
            (pt.tangent_impulse[1] - 1.8).abs() < 1e-10,
            "tangent[1] should be 1.8"
        );
    }

    #[test]
    fn test_age_manifold_warm_start() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 5.0);
        pt.tangent_impulse = [2.0, 0.0];
        m.points.push(pt);
        age_manifold_warm_start(&mut m, 0.5);
        assert!(
            (m.points[0].normal_impulse - 2.5).abs() < 1e-10,
            "normal_impulse should be 2.5"
        );
        assert!(
            (m.points[0].tangent_impulse[0] - 1.0).abs() < 1e-10,
            "tangent[0] should be 1.0"
        );
    }

    #[test]
    fn test_age_cache_warm_start() {
        let mut cache = ContactManifoldCache::new(5);
        let m = cache.get_or_create(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 8.0);
        pt.tangent_impulse = [4.0, 0.0];
        m.points.push(pt);
        age_cache_warm_start(&mut cache, 0.5);
        let m2 = cache.manifolds.get(&(0, 1)).unwrap();
        assert!(
            (m2.points[0].normal_impulse - 4.0).abs() < 1e-10,
            "normal_impulse should be 4.0"
        );
    }

    // ─── baumgarte_correction ──────────────────────────────────────────────────

    #[test]
    fn test_baumgarte_correction_basic() {
        let corr = baumgarte_correction(0.1, 0.2, 0.016);
        let expected = 0.2 * 0.1 / 0.016;
        assert!(
            (corr - expected).abs() < 1e-10,
            "Expected {expected}, got {corr}"
        );
    }

    #[test]
    fn test_baumgarte_correction_zero_depth() {
        let corr = baumgarte_correction(0.0, 0.2, 0.016);
        assert_eq!(corr, 0.0, "Zero depth should produce no correction");
    }

    #[test]
    fn test_baumgarte_correction_slop_no_correction_below_slop() {
        let corr = baumgarte_correction_slop(0.003, 0.2, 0.016, 0.005);
        assert_eq!(corr, 0.0, "Depth below slop should produce no correction");
    }

    #[test]
    fn test_baumgarte_correction_slop_correction_above_slop() {
        let corr = baumgarte_correction_slop(0.01, 0.2, 0.016, 0.005);
        let expected = baumgarte_correction(0.005, 0.2, 0.016);
        assert!(
            (corr - expected).abs() < 1e-10,
            "Expected {expected}, got {corr}"
        );
    }

    #[test]
    fn test_apply_position_corrections_returns_per_point() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        m.points
            .push(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0));
        m.points
            .push(make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.02, 0.0));
        let corrections = apply_position_corrections(&m, 0.2, 0.016, 0.005);
        assert_eq!(
            corrections.len(),
            2,
            "Should return one correction per point"
        );
        assert!(
            corrections[0] >= 0.0 && corrections[1] >= 0.0,
            "Corrections should be non-negative"
        );
    }

    // ─── ContactIsland ─────────────────────────────────────────────────────────

    #[test]
    fn test_contact_island_default() {
        let island = ContactIsland::default();
        assert_eq!(island.body_count(), 0);
        assert_eq!(island.contact_count(), 0);
    }

    #[test]
    fn test_build_contact_islands_single_pair() {
        let mut cache = ContactManifoldCache::new(5);
        let m = cache.get_or_create(0, 1, 0.5, 0.3);
        m.add_or_update(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0));
        let islands = build_contact_islands(&cache);
        assert_eq!(islands.len(), 1, "One pair should form one island");
        assert_eq!(islands[0].body_count(), 2, "Island should have 2 bodies");
    }

    #[test]
    fn test_build_contact_islands_two_separate_pairs() {
        let mut cache = ContactManifoldCache::new(5);
        let m1 = cache.get_or_create(0, 1, 0.5, 0.3);
        m1.add_or_update(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0));
        let m2 = cache.get_or_create(2, 3, 0.5, 0.3);
        m2.add_or_update(make_contact([5.0, 0.0, 0.0], [5.0, -0.01, 0.0], 0.01, 0.0));
        let islands = build_contact_islands(&cache);
        // Two disconnected pairs → two islands
        assert_eq!(
            islands.len(),
            2,
            "Two separate pairs should form two islands"
        );
    }

    #[test]
    fn test_build_contact_islands_connected_chain() {
        // Bodies 0-1-2 form a chain: should be one island
        let mut cache = ContactManifoldCache::new(5);
        let m1 = cache.get_or_create(0, 1, 0.5, 0.3);
        m1.add_or_update(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 0.0));
        let m2 = cache.get_or_create(1, 2, 0.5, 0.3);
        m2.add_or_update(make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.01, 0.0));
        let islands = build_contact_islands(&cache);
        assert_eq!(islands.len(), 1, "Connected chain should form one island");
        assert_eq!(
            islands[0].body_count(),
            3,
            "Chain of 3 bodies should have 3 bodies in island"
        );
    }

    // ─── ManifoldMetrics ───────────────────────────────────────────────────────

    #[test]
    fn test_manifold_metrics_empty() {
        let m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let metrics = compute_manifold_metrics(&m);
        assert_eq!(metrics.contact_count, 0);
        assert_eq!(metrics.max_depth, 0.0);
        assert!(!metrics.is_warm);
    }

    #[test]
    fn test_manifold_metrics_single_point() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        m.points
            .push(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.05, 0.0));
        let metrics = compute_manifold_metrics(&m);
        assert_eq!(metrics.contact_count, 1);
        assert!(
            (metrics.max_depth - 0.05).abs() < 1e-10,
            "max_depth should be 0.05"
        );
        assert!(
            (metrics.avg_depth - 0.05).abs() < 1e-10,
            "avg_depth should be 0.05"
        );
        assert_eq!(metrics.spread, 0.0, "Single point has zero spread");
    }

    #[test]
    fn test_manifold_metrics_two_points() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        m.points
            .push(make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.02, 0.0));
        m.points
            .push(make_contact([1.0, 0.0, 0.0], [1.0, -0.01, 0.0], 0.04, 0.0));
        let metrics = compute_manifold_metrics(&m);
        assert_eq!(metrics.contact_count, 2);
        assert!((metrics.max_depth - 0.04).abs() < 1e-10);
        assert!((metrics.avg_depth - 0.03).abs() < 1e-10);
        assert!(
            (metrics.spread - 1.0).abs() < 1e-10,
            "Spread should be 1.0, got {}",
            metrics.spread
        );
    }

    #[test]
    fn test_manifold_metrics_is_warm_with_impulse() {
        let mut m = PersistentManifold::new(0, 1, 0.5, 0.3);
        let mut pt = make_contact([0.0, 0.0, 0.0], [0.0, -0.01, 0.0], 0.01, 5.0);
        pt.normal_impulse = 5.0;
        m.points.push(pt);
        let metrics = compute_manifold_metrics(&m);
        assert!(
            metrics.is_warm,
            "Manifold with non-zero normal_impulse should be warm"
        );
    }
}
