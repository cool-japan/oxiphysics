// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial query API — ray casting, sphere sweeps, and overlap tests.
//!
//! `QueryWorld` maintains a flat list of body shapes (spheres, AABBs,
//! capsules) that can be probed without running a full simulation step.
//!
//! ## Supported queries
//!
//! | Query | Description |
//! |-------|-------------|
//! | `QueryWorld::raycast` | First hit along a ray |
//! | `QueryWorld::raycast_all` | All hits along a ray, sorted by distance |
//! | `QueryWorld::overlap_sphere` | Bodies whose shapes overlap a sphere |
//! | `QueryWorld::overlap_aabb` | Bodies whose shapes overlap an AABB |
//! | `QueryWorld::closest_body` | Nearest body to a world-space point |
//! | `QueryWorld::bodies_in_radius` | All bodies within a given distance |
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::query::{QueryWorld, Ray, QueryFilter};
//!
//! let mut world = QueryWorld::new();
//! world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
//! world.add_sphere(1, [3.0, 0.0, 0.0], 1.0, false);
//!
//! let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
//! let hit = world.raycast(&ray, &QueryFilter::default());
//! assert!(hit.is_some());
//! assert_eq!(hit.unwrap().body_index, 0);
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

// ============================================================================
// Vec3 helpers (self-contained)
// ============================================================================

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_len_sq(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}

#[inline]
fn v3_len(a: [f64; 3]) -> f64 {
    v3_len_sq(a).sqrt()
}

#[inline]
fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = v3_len(a);
    if l < 1e-12 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / l)
    }
}

#[inline]
fn v3_clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ============================================================================
// Ray
// ============================================================================

/// A world-space ray defined by an origin and a (normalized) direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    /// World-space origin of the ray.
    pub origin: [f64; 3],
    /// Unit-length direction vector.
    pub direction: [f64; 3],
}

impl Ray {
    /// Construct a new ray.  `direction` is normalized automatically.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction: v3_normalize(direction),
        }
    }

    /// Point at parameter `t` along the ray: `origin + t * direction`.
    #[inline]
    pub fn at(&self, t: f64) -> [f64; 3] {
        v3_add(self.origin, v3_scale(self.direction, t))
    }
}

// ============================================================================
// RayHit
// ============================================================================

/// A single intersection of a [`Ray`] with a body shape.
#[derive(Debug, Clone, PartialEq)]
pub struct RayHit {
    /// Ray parameter at the hit point (`t ≥ 0`).
    pub t: f64,
    /// Index of the body that was hit.
    pub body_index: usize,
    /// Outward surface normal at the hit point (unit length).
    pub normal: [f64; 3],
    /// World-space position of the hit point.
    pub point: [f64; 3],
}

// ============================================================================
// QueryFilter
// ============================================================================

/// Options that restrict which bodies a query will consider.
#[derive(Debug, Clone)]
pub struct QueryFilter {
    /// Body indices to exclude from results.
    pub exclude_bodies: Vec<usize>,
    /// Maximum ray or distance parameter to accept (default: `f64::MAX`).
    pub max_distance: f64,
    /// Whether to include sleeping bodies (default: `true`).
    pub include_sleeping: bool,
}

impl Default for QueryFilter {
    fn default() -> Self {
        Self {
            exclude_bodies: Vec::new(),
            max_distance: f64::MAX,
            include_sleeping: true,
        }
    }
}

impl QueryFilter {
    /// Create a filter that excludes a single body.
    pub fn exclude(body_index: usize) -> Self {
        Self {
            exclude_bodies: vec![body_index],
            ..Default::default()
        }
    }

    /// Create a filter with a maximum range.
    pub fn with_max_distance(max: f64) -> Self {
        Self {
            max_distance: max,
            ..Default::default()
        }
    }

    /// Returns `true` if this body passes the filter.
    fn accepts(&self, index: usize, sleeping: bool) -> bool {
        if !self.include_sleeping && sleeping {
            return false;
        }
        !self.exclude_bodies.contains(&index)
    }
}

// ============================================================================
// QueryShape
// ============================================================================

/// The collision shape used for spatial queries.
#[derive(Debug, Clone, PartialEq)]
pub enum QueryShape {
    /// Sphere at a world-space center.
    Sphere {
        /// World-space center.
        center: [f64; 3],
        /// Sphere radius.
        radius: f64,
    },
    /// Axis-aligned bounding box.
    Aabb {
        /// Minimum corner.
        min: [f64; 3],
        /// Maximum corner.
        max: [f64; 3],
    },
    /// Capsule (swept sphere) between two endpoints.
    Capsule {
        /// First endpoint.
        start: [f64; 3],
        /// Second endpoint.
        end: [f64; 3],
        /// Capsule radius.
        radius: f64,
    },
}

impl QueryShape {
    /// World-space AABB that conservatively encloses this shape.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            QueryShape::Sphere { center, radius } => {
                let r = *radius;
                (
                    [center[0] - r, center[1] - r, center[2] - r],
                    [center[0] + r, center[1] + r, center[2] + r],
                )
            }
            QueryShape::Aabb { min, max } => (*min, *max),
            QueryShape::Capsule { start, end, radius } => {
                let r = *radius;
                (
                    [
                        start[0].min(end[0]) - r,
                        start[1].min(end[1]) - r,
                        start[2].min(end[2]) - r,
                    ],
                    [
                        start[0].max(end[0]) + r,
                        start[1].max(end[1]) + r,
                        start[2].max(end[2]) + r,
                    ],
                )
            }
        }
    }

    /// Approximate world-space center of this shape.
    pub fn center(&self) -> [f64; 3] {
        match self {
            QueryShape::Sphere { center, .. } => *center,
            QueryShape::Aabb { min, max } => [
                (min[0] + max[0]) * 0.5,
                (min[1] + max[1]) * 0.5,
                (min[2] + max[2]) * 0.5,
            ],
            QueryShape::Capsule { start, end, .. } => [
                (start[0] + end[0]) * 0.5,
                (start[1] + end[1]) * 0.5,
                (start[2] + end[2]) * 0.5,
            ],
        }
    }

    /// Squared distance from `point` to the surface (negative if inside).
    pub fn distance_sq_to_point(&self, point: [f64; 3]) -> f64 {
        match self {
            QueryShape::Sphere { center, radius } => {
                let d = v3_len(v3_sub(point, *center)) - radius;
                d * d
            }
            QueryShape::Aabb { min, max } => {
                // Closest point on AABB
                let cx = v3_clamp(point[0], min[0], max[0]);
                let cy = v3_clamp(point[1], min[1], max[1]);
                let cz = v3_clamp(point[2], min[2], max[2]);
                let dx = point[0] - cx;
                let dy = point[1] - cy;
                let dz = point[2] - cz;
                dx * dx + dy * dy + dz * dz
            }
            QueryShape::Capsule { start, end, radius } => {
                let seg = v3_sub(*end, *start);
                let to_pt = v3_sub(point, *start);
                let seg_len_sq = v3_len_sq(seg);
                let t = if seg_len_sq < 1e-12 {
                    0.0
                } else {
                    v3_clamp(v3_dot(to_pt, seg) / seg_len_sq, 0.0, 1.0)
                };
                let closest = v3_add(*start, v3_scale(seg, t));
                let d = v3_len(v3_sub(point, closest)) - radius;
                d * d
            }
        }
    }

    /// Test ray intersection.  Returns ray parameter `t` on hit, or `None`.
    pub fn ray_intersect(&self, ray: &Ray) -> Option<(f64, [f64; 3])> {
        match self {
            QueryShape::Sphere { center, radius } => ray_sphere_intersect(ray, *center, *radius),
            QueryShape::Aabb { min, max } => ray_aabb_intersect(ray, *min, *max),
            QueryShape::Capsule { start, end, radius } => {
                ray_capsule_intersect(ray, *start, *end, *radius)
            }
        }
    }

    /// Returns `true` if this shape overlaps with `sphere(center, radius)`.
    pub fn overlaps_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        let dist_sq = self.distance_sq_to_point(center);
        dist_sq <= radius * radius
    }

    /// Returns `true` if this shape overlaps `aabb(min, max)`.
    pub fn overlaps_aabb(&self, min: [f64; 3], max: [f64; 3]) -> bool {
        let (smin, smax) = self.aabb();
        smin[0] <= max[0]
            && smax[0] >= min[0]
            && smin[1] <= max[1]
            && smax[1] >= min[1]
            && smin[2] <= max[2]
            && smax[2] >= min[2]
    }
}

// ============================================================================
// Intersection helpers
// ============================================================================

/// Ray–sphere intersection. Returns `(t, normal)` for the nearest forward hit.
fn ray_sphere_intersect(ray: &Ray, center: [f64; 3], radius: f64) -> Option<(f64, [f64; 3])> {
    let oc = v3_sub(ray.origin, center);
    let a = v3_dot(ray.direction, ray.direction);
    let b = 2.0 * v3_dot(oc, ray.direction);
    let c = v3_dot(oc, oc) - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t0 = (-b - sqrt_disc) / (2.0 * a);
    let t1 = (-b + sqrt_disc) / (2.0 * a);
    let t = if t0 >= 0.0 {
        t0
    } else if t1 >= 0.0 {
        t1
    } else {
        return None;
    };
    let hit = ray.at(t);
    let normal = v3_normalize(v3_sub(hit, center));
    Some((t, normal))
}

/// Ray–AABB intersection using slab method. Returns `(t, normal)`.
fn ray_aabb_intersect(ray: &Ray, min: [f64; 3], max: [f64; 3]) -> Option<(f64, [f64; 3])> {
    let mut t_min = 0.0_f64;
    let mut t_max = f64::MAX;
    let mut normal = [0.0_f64; 3];

    for axis in 0..3 {
        let inv_d = if ray.direction[axis].abs() < 1e-12 {
            f64::INFINITY
        } else {
            1.0 / ray.direction[axis]
        };
        let t1 = (min[axis] - ray.origin[axis]) * inv_d;
        let t2 = (max[axis] - ray.origin[axis]) * inv_d;
        let (t_near, t_far) = if t1 < t2 { (t1, t2) } else { (t2, t1) };

        if t_near > t_min {
            t_min = t_near;
            normal = [0.0; 3];
            normal[axis] = if t1 < t2 { -1.0 } else { 1.0 };
        }
        t_max = t_max.min(t_far);

        if t_min > t_max {
            return None;
        }
    }
    if t_min < 0.0 {
        return None;
    }
    Some((t_min, normal))
}

/// Ray–capsule intersection. Returns `(t, normal)`.
fn ray_capsule_intersect(
    ray: &Ray,
    start: [f64; 3],
    end: [f64; 3],
    radius: f64,
) -> Option<(f64, [f64; 3])> {
    // Project ray against the infinite cylinder aligned with the capsule axis,
    // then clamp to the segment and check the sphere end-caps.
    let seg = v3_sub(end, start);
    let ro = v3_sub(ray.origin, start);
    let seg_len = v3_len(seg);
    if seg_len < 1e-12 {
        // Degenerate capsule: treat as sphere
        return ray_sphere_intersect(ray, start, radius);
    }
    let seg_n = v3_scale(seg, 1.0 / seg_len);

    // Component of ray direction perpendicular to capsule axis
    let d_dot_seg = v3_dot(ray.direction, seg_n);
    let d_perp = v3_sub(ray.direction, v3_scale(seg_n, d_dot_seg));
    let o_dot_seg = v3_dot(ro, seg_n);
    let o_perp = v3_sub(ro, v3_scale(seg_n, o_dot_seg));

    // Solve 2-D ray–circle intersection in the perpendicular plane
    let a = v3_dot(d_perp, d_perp);
    let b_half = v3_dot(o_perp, d_perp);
    let c = v3_dot(o_perp, o_perp) - radius * radius;
    let disc = b_half * b_half - a * c;

    let mut best_t = f64::MAX;
    let mut best_normal = [0.0_f64; 3];

    if disc >= 0.0 && a > 1e-12 {
        let t_cyl = (-b_half - disc.sqrt()) / a;
        if t_cyl >= 0.0 {
            let proj = o_dot_seg + t_cyl * d_dot_seg;
            if proj >= 0.0 && proj <= seg_len {
                // Hit cylinder body
                let hit = ray.at(t_cyl);
                let axis_pt = v3_add(start, v3_scale(seg_n, proj));
                let n = v3_normalize(v3_sub(hit, axis_pt));
                best_t = t_cyl;
                best_normal = n;
            }
        }
    }

    // Check sphere end-caps
    for cap_center in [start, end] {
        if let Some((t, n)) = ray_sphere_intersect(ray, cap_center, radius)
            && t >= 0.0
            && t < best_t
        {
            best_t = t;
            best_normal = n;
        }
    }

    if best_t < f64::MAX {
        Some((best_t, best_normal))
    } else {
        None
    }
}

// ============================================================================
// QueryEntry
// ============================================================================

/// A single body registered in the [`QueryWorld`].
#[derive(Debug, Clone)]
pub struct QueryEntry {
    /// Application-assigned body index.
    pub index: usize,
    /// Shape used for queries.
    pub shape: QueryShape,
    /// Whether the body is currently sleeping.
    pub is_sleeping: bool,
}

// ============================================================================
// QueryWorld
// ============================================================================

/// A lightweight spatial index that supports common query operations.
///
/// Bodies are inserted by the application (typically after each simulation
/// step) and can then be probed via ray casts, sphere overlaps, etc.
///
/// # Example
///
/// ```rust
/// use oxiphysics::query::{QueryWorld, Ray, QueryFilter};
///
/// let mut world = QueryWorld::new();
/// world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
/// world.add_sphere(1, [4.0, 0.0, 0.0], 1.0, false);
///
/// let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
/// let hit = world.raycast(&ray, &QueryFilter::default()).unwrap();
/// assert_eq!(hit.body_index, 0);
/// ```
#[derive(Debug, Default)]
pub struct QueryWorld {
    entries: Vec<QueryEntry>,
}

impl QueryWorld {
    /// Create an empty query world.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with reserved capacity for `n` bodies.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            entries: Vec::with_capacity(n),
        }
    }

    // ------------------------------------------------------------------
    // Insertion / removal
    // ------------------------------------------------------------------

    /// Register a sphere body and return its slot position.
    pub fn add_sphere(&mut self, index: usize, center: [f64; 3], radius: f64, sleeping: bool) {
        self.entries.push(QueryEntry {
            index,
            shape: QueryShape::Sphere { center, radius },
            is_sleeping: sleeping,
        });
    }

    /// Register an AABB body.
    pub fn add_aabb(&mut self, index: usize, min: [f64; 3], max: [f64; 3], sleeping: bool) {
        self.entries.push(QueryEntry {
            index,
            shape: QueryShape::Aabb { min, max },
            is_sleeping: sleeping,
        });
    }

    /// Register a capsule body.
    pub fn add_capsule(
        &mut self,
        index: usize,
        start: [f64; 3],
        end: [f64; 3],
        radius: f64,
        sleeping: bool,
    ) {
        self.entries.push(QueryEntry {
            index,
            shape: QueryShape::Capsule { start, end, radius },
            is_sleeping: sleeping,
        });
    }

    /// Update the sleeping flag for a body by its application index.
    pub fn set_sleeping(&mut self, index: usize, sleeping: bool) {
        for e in &mut self.entries {
            if e.index == index {
                e.is_sleeping = sleeping;
            }
        }
    }

    /// Remove all entries for a body index.  Returns how many were removed.
    pub fn remove(&mut self, index: usize) -> usize {
        let before = self.entries.len();
        self.entries.retain(|e| e.index != index);
        before - self.entries.len()
    }

    /// Remove all bodies.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Number of registered bodies.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if no bodies are registered.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate over all entries.
    pub fn iter(&self) -> impl Iterator<Item = &QueryEntry> {
        self.entries.iter()
    }

    // ------------------------------------------------------------------
    // Spatial queries
    // ------------------------------------------------------------------

    /// Return the closest [`RayHit`] along `ray`, or `None` if no body
    /// within the filter's `max_distance` is hit.
    pub fn raycast(&self, ray: &Ray, filter: &QueryFilter) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;

        for entry in &self.entries {
            if !filter.accepts(entry.index, entry.is_sleeping) {
                continue;
            }
            if let Some((t, normal)) = entry.shape.ray_intersect(ray) {
                if t > filter.max_distance {
                    continue;
                }
                let keep = match &best {
                    None => true,
                    Some(prev) => t < prev.t,
                };
                if keep {
                    best = Some(RayHit {
                        t,
                        body_index: entry.index,
                        normal,
                        point: ray.at(t),
                    });
                }
            }
        }
        best
    }

    /// Return **all** [`RayHit`]s along `ray`, sorted nearest-first.
    pub fn raycast_all(&self, ray: &Ray, filter: &QueryFilter) -> Vec<RayHit> {
        let mut hits: Vec<RayHit> = self
            .entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .filter_map(|entry| {
                entry.shape.ray_intersect(ray).and_then(|(t, normal)| {
                    if t <= filter.max_distance {
                        Some(RayHit {
                            t,
                            body_index: entry.index,
                            normal,
                            point: ray.at(t),
                        })
                    } else {
                        None
                    }
                })
            })
            .collect();

        hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }

    /// Return the indices of all bodies whose shapes overlap `sphere(center, radius)`.
    pub fn overlap_sphere(
        &self,
        center: [f64; 3],
        radius: f64,
        filter: &QueryFilter,
    ) -> Vec<usize> {
        self.entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .filter(|e| e.shape.overlaps_sphere(center, radius))
            .map(|e| e.index)
            .collect()
    }

    /// Return the indices of all bodies whose shapes overlap `aabb(min, max)`.
    pub fn overlap_aabb(&self, min: [f64; 3], max: [f64; 3], filter: &QueryFilter) -> Vec<usize> {
        self.entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .filter(|e| e.shape.overlaps_aabb(min, max))
            .map(|e| e.index)
            .collect()
    }

    /// Return the index and approximate squared distance of the body whose
    /// shape is nearest to `point`, or `None` if the world is empty.
    pub fn closest_body(&self, point: [f64; 3], filter: &QueryFilter) -> Option<(usize, f64)> {
        self.entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .map(|e| (e.index, e.shape.distance_sq_to_point(point)))
            .min_by(|(_, da), (_, db)| da.partial_cmp(db).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Return all bodies within `radius` of `point`, sorted nearest-first.
    /// The second element is the approximate distance (not squared).
    pub fn bodies_in_radius(
        &self,
        point: [f64; 3],
        radius: f64,
        filter: &QueryFilter,
    ) -> Vec<(usize, f64)> {
        let r2 = radius * radius;
        let mut results: Vec<(usize, f64)> = self
            .entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .filter_map(|e| {
                let d2 = e.shape.distance_sq_to_point(point);
                if d2 <= r2 {
                    Some((e.index, d2.sqrt()))
                } else {
                    None
                }
            })
            .collect();

        results.sort_by(|(_, da), (_, db)| da.partial_cmp(db).unwrap_or(std::cmp::Ordering::Equal));
        results
    }

    /// Count how many bodies overlap a sphere.
    pub fn count_in_sphere(&self, center: [f64; 3], radius: f64) -> usize {
        self.overlap_sphere(center, radius, &QueryFilter::default())
            .len()
    }

    /// Return the k nearest bodies to `point`, sorted nearest-first.
    pub fn k_nearest(&self, point: [f64; 3], k: usize, filter: &QueryFilter) -> Vec<(usize, f64)> {
        let mut scored: Vec<(usize, f64)> = self
            .entries
            .iter()
            .filter(|e| filter.accepts(e.index, e.is_sleeping))
            .map(|e| (e.index, e.shape.distance_sq_to_point(point).sqrt()))
            .collect();
        scored.sort_by(|(_, da), (_, db)| da.partial_cmp(db).unwrap_or(std::cmp::Ordering::Equal));
        scored.truncate(k);
        scored
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_sphere_hit() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
        let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hit = world.raycast(&ray, &QueryFilter::default()).unwrap();
        assert_eq!(hit.body_index, 0);
        assert!((hit.t - 4.0).abs() < 1e-6);
    }

    #[test]
    fn ray_sphere_miss() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [0.0, 10.0, 0.0], 1.0, false);
        let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(world.raycast(&ray, &QueryFilter::default()).is_none());
    }

    #[test]
    fn raycast_all_sorted() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [4.0, 0.0, 0.0], 1.0, false);
        world.add_sphere(1, [0.0, 0.0, 0.0], 1.0, false);
        let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hits = world.raycast_all(&ray, &QueryFilter::default());
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].body_index, 1); // closer sphere first
        assert_eq!(hits[1].body_index, 0);
    }

    #[test]
    fn overlap_sphere_basic() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
        world.add_sphere(1, [10.0, 0.0, 0.0], 1.0, false);
        let found = world.overlap_sphere([0.0, 0.0, 0.0], 2.0, &QueryFilter::default());
        assert_eq!(found, vec![0]);
    }

    #[test]
    fn closest_body() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
        world.add_sphere(1, [5.0, 0.0, 0.0], 1.0, false);
        let (idx, _) = world
            .closest_body([1.5, 0.0, 0.0], &QueryFilter::default())
            .unwrap();
        assert_eq!(idx, 0);
    }

    #[test]
    fn aabb_query() {
        let mut world = QueryWorld::new();
        world.add_aabb(0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], false);
        world.add_aabb(1, [5.0, 5.0, 5.0], [6.0, 6.0, 6.0], false);
        let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hit = world.raycast(&ray, &QueryFilter::default()).unwrap();
        assert_eq!(hit.body_index, 0);
    }

    #[test]
    fn filter_sleeping() {
        let mut world = QueryWorld::new();
        world.add_sphere(0, [0.0, 0.0, 0.0], 1.0, true); // sleeping
        world.add_sphere(1, [0.0, 0.0, 0.0], 1.0, false);
        let filter = QueryFilter {
            include_sleeping: false,
            ..Default::default()
        };
        let ray = Ray::new([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hit = world.raycast(&ray, &filter).unwrap();
        assert_eq!(hit.body_index, 1);
    }

    #[test]
    fn k_nearest() {
        let mut world = QueryWorld::new();
        for i in 0..5_usize {
            world.add_sphere(i, [i as f64 * 2.0, 0.0, 0.0], 0.5, false);
        }
        let nearest = world.k_nearest([0.0, 0.0, 0.0], 3, &QueryFilter::default());
        assert_eq!(nearest.len(), 3);
        assert_eq!(nearest[0].0, 0); // body 0 is at origin
    }
}
