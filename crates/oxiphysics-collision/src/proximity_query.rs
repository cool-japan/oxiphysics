// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced proximity queries for collision detection.
//!
//! This module provides closest-point and signed-distance-field computations
//! for common geometric primitives, along with CSG SDF operations and a
//! warm-start proximity cache.
//!
//! All positions are represented as `[f64; 3]` (x, y, z) to avoid external
//! algebra dependencies.

// ─────────────────────────────────────────────────────────────────────────────
// Low-level vector helpers (no external dependencies)
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vadd(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vsub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vscale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vcross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vlen_sq(a: [f64; 3]) -> f64 {
    vdot(a, a)
}

#[inline]
fn vlen(a: [f64; 3]) -> f64 {
    vlen_sq(a).sqrt()
}

#[inline]
fn vnormalize(a: [f64; 3]) -> [f64; 3] {
    let l = vlen(a);
    if l < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        vscale(a, 1.0 / l)
    }
}

#[inline]
fn vlerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    vadd(vscale(a, 1.0 - t), vscale(b, t))
}

#[inline]
fn vdist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vlen(vsub(a, b))
}

#[inline]
fn vclamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public data types
// ─────────────────────────────────────────────────────────────────────────────

/// Identifies which geometric feature of a primitive is closest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClosestFeature {
    /// A vertex / corner of the primitive.
    Vertex,
    /// An edge of the primitive.
    Edge,
    /// A face (interior) of the primitive.
    Face,
}

/// Result returned by closest-point queries.
///
/// Contains the closest position, the separation distance, which feature was
/// nearest, and optional barycentric parameters for that feature.
#[derive(Debug, Clone, Copy)]
pub struct ClosestPointResult {
    /// The closest point on the primitive expressed in world space.
    pub position: [f64; 3],
    /// Non-negative Euclidean distance from the query point to `position`.
    pub distance: f64,
    /// Which geometric feature of the primitive is closest.
    pub feature: ClosestFeature,
    /// First barycentric / interpolation parameter (e.g. `t` along an edge).
    pub u_param: f64,
    /// Second barycentric parameter (e.g. `v` in a triangle).
    pub v_param: f64,
}

impl ClosestPointResult {
    fn new(
        position: [f64; 3],
        query_pt: [f64; 3],
        feature: ClosestFeature,
        u_param: f64,
        v_param: f64,
    ) -> Self {
        let distance = vdist(position, query_pt);
        Self {
            position,
            distance,
            feature,
            u_param,
            v_param,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Segment – segment closest points helper result
// ─────────────────────────────────────────────────────────────────────────────

/// Result from a segment–segment closest-point query.
#[derive(Debug, Clone, Copy)]
pub struct SegSegResult {
    /// Closest point on segment AB.
    pub point_a: [f64; 3],
    /// Closest point on segment CD.
    pub point_b: [f64; 3],
    /// Parameter along AB (0 → A, 1 → B).
    pub s: f64,
    /// Parameter along CD (0 → C, 1 → D).
    pub t: f64,
    /// Euclidean distance between the two closest points.
    pub distance: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Point – segment
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest point on segment `[a, b]` to `p`.
///
/// Returns a [`ClosestPointResult`] whose `u_param` is the clamped parameter
/// `t ∈ [0, 1]` such that `closest = a + t * (b - a)`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_segment_closest;
/// let result = point_segment_closest([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]);
/// assert!((result.u_param - 0.5).abs() < 1e-9);
/// ```
pub fn point_segment_closest(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> ClosestPointResult {
    let ab = vsub(b, a);
    let ap = vsub(p, a);
    let len_sq = vlen_sq(ab);

    let (t, feature) = if len_sq < 1e-300 {
        // Degenerate segment – a == b
        (0.0, ClosestFeature::Vertex)
    } else {
        let raw_t = vdot(ap, ab) / len_sq;
        if raw_t <= 0.0 {
            (0.0, ClosestFeature::Vertex)
        } else if raw_t >= 1.0 {
            (1.0, ClosestFeature::Vertex)
        } else {
            (raw_t, ClosestFeature::Edge)
        }
    };

    let pos = vadd(a, vscale(ab, t));
    ClosestPointResult::new(pos, p, feature, t, 0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Point – triangle  (Ericson / Voronoi region method)
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest point on triangle `(v0, v1, v2)` to point `p` using
/// Voronoi-region decomposition (Real-Time Collision Detection, Ericson §5.1.5).
///
/// The returned `u_param` and `v_param` are the barycentric coordinates
/// corresponding to `v1` and `v2` respectively (so `w = 1 - u - v` for `v0`).
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_triangle_closest;
/// // Point directly above centroid
/// let r = point_triangle_closest(
///     [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0],
///     [1.0, 1.0, 2.0]);
/// assert!(r.distance < 2.1);
/// ```
pub fn point_triangle_closest(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    p: [f64; 3],
) -> ClosestPointResult {
    let ab = vsub(v1, v0);
    let ac = vsub(v2, v0);
    let ap = vsub(p, v0);

    let d1 = vdot(ab, ap);
    let d2 = vdot(ac, ap);
    // Vertex region v0
    if d1 <= 0.0 && d2 <= 0.0 {
        return ClosestPointResult::new(v0, p, ClosestFeature::Vertex, 0.0, 0.0);
    }

    let bp = vsub(p, v1);
    let d3 = vdot(ab, bp);
    let d4 = vdot(ac, bp);
    // Vertex region v1
    if d3 >= 0.0 && d4 <= d3 {
        return ClosestPointResult::new(v1, p, ClosestFeature::Vertex, 1.0, 0.0);
    }

    // Edge region v0-v1
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let pos = vadd(v0, vscale(ab, v));
        return ClosestPointResult::new(pos, p, ClosestFeature::Edge, v, 0.0);
    }

    let cp = vsub(p, v2);
    let d5 = vdot(ab, cp);
    let d6 = vdot(ac, cp);
    // Vertex region v2
    if d6 >= 0.0 && d5 <= d6 {
        return ClosestPointResult::new(v2, p, ClosestFeature::Vertex, 0.0, 1.0);
    }

    // Edge region v0-v2
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let pos = vadd(v0, vscale(ac, w));
        return ClosestPointResult::new(pos, p, ClosestFeature::Edge, 0.0, w);
    }

    // Edge region v1-v2
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let pos = vadd(v1, vscale(vsub(v2, v1), w));
        return ClosestPointResult::new(pos, p, ClosestFeature::Edge, 1.0 - w, w);
    }

    // Interior of triangle
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let pos = vadd(v0, vadd(vscale(ab, v), vscale(ac, w)));
    ClosestPointResult::new(pos, p, ClosestFeature::Face, v, w)
}

// ─────────────────────────────────────────────────────────────────────────────
// Segment – segment
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the closest points between segments `AB` and `CD`.
///
/// Handles the degenerate (parallel / point) cases robustly.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::segment_segment_closest;
/// let r = segment_segment_closest(
///     [0.0,0.0,0.0],[1.0,0.0,0.0],
///     [0.5,1.0,0.0],[0.5,-1.0,0.0]);
/// assert!(r.distance < 1e-9);
/// ```
pub fn segment_segment_closest(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> SegSegResult {
    let d1 = vsub(b, a);
    let d2 = vsub(d, c);
    let r = vsub(a, c);

    let e = vdot(d1, d1);
    let f = vdot(d2, d2);

    let (mut s, t);

    const EPS: f64 = 1e-10;

    // Both degenerate
    if e < EPS && f < EPS {
        s = 0.0;
        t = 0.0;
    } else if e < EPS {
        // Segment AB is a point
        s = 0.0;
        t = vclamp(vdot(r, d2) / f, 0.0, 1.0);
    } else {
        let c1 = vdot(d1, r);
        if f < EPS {
            // Segment CD is a point
            t = 0.0;
            s = vclamp(-c1 / e, 0.0, 1.0);
        } else {
            let b_coeff = vdot(d1, d2);
            let denom = e * f - b_coeff * b_coeff;
            // Not parallel
            if denom.abs() > EPS {
                let c2 = vdot(d2, r);
                s = vclamp((b_coeff * c2 - c1 * f) / denom, 0.0, 1.0);
            } else {
                // Parallel – pick s=0
                s = 0.0;
            }
            let c2 = vdot(d2, r);
            let tnom = b_coeff * s + c2;
            if tnom < 0.0 {
                t = 0.0;
                s = vclamp(-c1 / e, 0.0, 1.0);
            } else if tnom > f {
                t = 1.0;
                s = vclamp((b_coeff - c1) / e, 0.0, 1.0);
            } else {
                t = tnom / f;
            }
        }
    }

    let point_a = vadd(a, vscale(d1, s));
    let point_b = vadd(c, vscale(d2, t));
    let distance = vdist(point_a, point_b);
    SegSegResult {
        point_a,
        point_b,
        s,
        t,
        distance,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Point – AABB
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance from point `p` to axis-aligned bounding box defined by
/// `min` and `max` corners.
///
/// Positive when the point is outside the box, negative when inside.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_box_distance;
/// // Inside
/// assert!(point_box_distance([-1.0,-1.0,-1.0],[1.0,1.0,1.0],[0.0,0.0,0.0]) < 0.0);
/// // Outside
/// assert!(point_box_distance([-1.0,-1.0,-1.0],[1.0,1.0,1.0],[2.0,0.0,0.0]) > 0.0);
/// ```
pub fn point_box_distance(min: [f64; 3], max: [f64; 3], p: [f64; 3]) -> f64 {
    // Exterior distance squared
    let mut dist_sq_outside = 0.0_f64;
    for i in 0..3 {
        if p[i] < min[i] {
            let d = min[i] - p[i];
            dist_sq_outside += d * d;
        } else if p[i] > max[i] {
            let d = p[i] - max[i];
            dist_sq_outside += d * d;
        }
    }
    if dist_sq_outside > 0.0 {
        return dist_sq_outside.sqrt();
    }

    // Interior: return negative of deepest penetration
    let mut min_pen = f64::MAX;
    for i in 0..3 {
        let d_lo = p[i] - min[i];
        let d_hi = max[i] - p[i];
        let pen = d_lo.min(d_hi);
        if pen < min_pen {
            min_pen = pen;
        }
    }
    -min_pen
}

// ─────────────────────────────────────────────────────────────────────────────
// Point – sphere
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance from `p` to a sphere with `center` and radius `r`.
///
/// Negative when the point is inside the sphere.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_sphere_distance;
/// assert!((point_sphere_distance([0.0,0.0,0.0], 1.0, [0.0,0.0,0.0]) - (-1.0)).abs() < 1e-9);
/// assert!((point_sphere_distance([0.0,0.0,0.0], 1.0, [2.0,0.0,0.0]) - 1.0).abs() < 1e-9);
/// ```
pub fn point_sphere_distance(center: [f64; 3], r: f64, p: [f64; 3]) -> f64 {
    vdist(p, center) - r
}

// ─────────────────────────────────────────────────────────────────────────────
// Point – capsule
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance from `p` to a capsule defined by segment endpoints `a`, `b`
/// and radius `r`.
///
/// Negative when the point is inside the capsule volume.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_capsule_distance;
/// // Point on the capsule axis → inside if within radius
/// let d = point_capsule_distance([0.0,0.0,0.0],[0.0,1.0,0.0], 0.5, [0.0,0.5,0.0]);
/// assert!(d < 0.0);
/// ```
pub fn point_capsule_distance(a: [f64; 3], b: [f64; 3], r: f64, p: [f64; 3]) -> f64 {
    let res = point_segment_closest(a, b, p);
    res.distance - r
}

// ─────────────────────────────────────────────────────────────────────────────
// Signed Distance Field
// ─────────────────────────────────────────────────────────────────────────────

/// A variant-based signed distance field (SDF) for common primitive shapes.
///
/// Use [`SignedDistanceField::eval`] to obtain the signed distance at any
/// world-space point.
#[derive(Debug, Clone)]
pub enum SignedDistanceField {
    /// A sphere defined by its `center` and radius `r`.
    Sphere {
        /// World-space center of the sphere.
        center: [f64; 3],
        /// Radius of the sphere.
        r: f64,
    },
    /// An axis-aligned box defined by `half`-extents centred at the origin.
    Box {
        /// Half-widths along each axis.
        half: [f64; 3],
    },
    /// A capsule (swept sphere) from point `a` to `b` with radius `r`.
    Capsule {
        /// Start of the capsule's medial axis.
        a: [f64; 3],
        /// End of the capsule's medial axis.
        b: [f64; 3],
        /// Radius of the capsule.
        r: f64,
    },
    /// An upright cylinder centred at the origin, aligned along `axis` (unit
    /// vector), with radius `r` and total height `h`.
    Cylinder {
        /// Unit direction along the cylinder's length.
        axis: [f64; 3],
        /// Radius of the cylinder.
        r: f64,
        /// Total height of the cylinder.
        h: f64,
    },
}

impl SignedDistanceField {
    /// Evaluates the signed distance at world-space point `p`.
    ///
    /// Returns a negative value when `p` is inside the shape.
    ///
    /// # Examples
    /// ```no_run
    /// use oxiphysics_collision::proximity_query::SignedDistanceField;
    /// let sdf = SignedDistanceField::Sphere { center: [0.0;3], r: 1.0 };
    /// assert!((sdf.eval([1.5,0.0,0.0]) - 0.5).abs() < 1e-9);
    /// ```
    pub fn eval(&self, p: [f64; 3]) -> f64 {
        match self {
            SignedDistanceField::Sphere { center, r } => point_sphere_distance(*center, *r, p),
            SignedDistanceField::Box { half } => {
                // Symmetric about origin
                let q = [
                    p[0].abs() - half[0],
                    p[1].abs() - half[1],
                    p[2].abs() - half[2],
                ];
                let q_pos = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
                let outside = vlen(q_pos);
                let inside = q[0].max(q[1]).max(q[2]).min(0.0);
                outside + inside
            }
            SignedDistanceField::Capsule { a, b, r } => point_capsule_distance(*a, *b, *r, p),
            SignedDistanceField::Cylinder { axis, r, h } => {
                // Project onto axis (assumed to pass through origin)
                let ax = vnormalize(*axis);
                let proj = vdot(p, ax);
                let half_h = h / 2.0;
                let radial_pt = vsub(p, vscale(ax, proj));
                let rd = vlen(radial_pt) - r;
                let yd = proj.abs() - half_h;
                // Distance to cylinder
                let outside_r = rd.max(0.0);
                let outside_y = yd.max(0.0);
                let outside = (outside_r * outside_r + outside_y * outside_y).sqrt();
                let inside = rd.max(yd).min(0.0);
                outside + inside
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CSG SDF operations
// ─────────────────────────────────────────────────────────────────────────────

/// Boolean union of two SDF values: the region that is inside *either* shape.
///
/// Returns `min(a, b)`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::sdf_union;
/// assert_eq!(sdf_union(2.0, -1.0), -1.0);
/// ```
#[inline]
pub fn sdf_union(a: f64, b: f64) -> f64 {
    a.min(b)
}

/// Boolean intersection of two SDF values: the region that is inside *both* shapes.
///
/// Returns `max(a, b)`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::sdf_intersect;
/// assert_eq!(sdf_intersect(2.0, -1.0), 2.0);
/// ```
#[inline]
pub fn sdf_intersect(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// Boolean subtraction of SDF `b` from `a`: the region inside `a` but outside `b`.
///
/// Returns `max(a, -b)`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::sdf_subtract;
/// assert!((sdf_subtract(-0.5, -0.3) - 0.3).abs() < 1e-9);
/// ```
#[inline]
pub fn sdf_subtract(a: f64, b: f64) -> f64 {
    a.max(-b)
}

// ─────────────────────────────────────────────────────────────────────────────
// Proximity cache
// ─────────────────────────────────────────────────────────────────────────────

/// Maximum number of entries in a [`ProximityCache`].
pub const PROXIMITY_CACHE_CAPACITY: usize = 64;

/// A cached closest-feature pair between two shapes identified by integer IDs.
#[derive(Debug, Clone, Copy)]
pub struct ProximityCacheEntry {
    /// ID of the first shape.
    pub id_a: u32,
    /// ID of the second shape.
    pub id_b: u32,
    /// Closest point on shape A.
    pub point_a: [f64; 3],
    /// Closest point on shape B.
    pub point_b: [f64; 3],
    /// Distance at the time of caching.
    pub distance: f64,
    /// Which feature of shape A produced the closest point.
    pub feature_a: ClosestFeature,
    /// Which feature of shape B produced the closest point.
    pub feature_b: ClosestFeature,
    /// Generation / frame stamp when this entry was last updated.
    pub stamp: u64,
}

/// LRU-eviction proximity cache holding up to [`PROXIMITY_CACHE_CAPACITY`]
/// closest-feature pairs.
///
/// Useful for warm-starting iterative contact solvers: re-using the previous
/// frame's closest features drastically reduces the number of GJK / distance
/// iterations required.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::{ProximityCache, ClosestFeature};
/// let mut cache = ProximityCache::new();
/// cache.insert(1, 2, [0.0;3], [1.0,0.0,0.0], 1.0,
///              ClosestFeature::Face, ClosestFeature::Face, 0);
/// assert!(cache.get(1, 2).is_some());
/// ```
#[derive(Debug, Clone)]
pub struct ProximityCache {
    entries: Vec<ProximityCacheEntry>,
}

impl Default for ProximityCache {
    fn default() -> Self {
        Self::new()
    }
}

impl ProximityCache {
    /// Creates a new, empty proximity cache.
    pub fn new() -> Self {
        Self {
            entries: Vec::with_capacity(PROXIMITY_CACHE_CAPACITY),
        }
    }

    /// Returns the number of entries currently held.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the cache contains no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Looks up a cached entry for the (unordered) pair `(id_a, id_b)`.
    ///
    /// Returns `None` when the pair has not been cached.
    pub fn get(&self, id_a: u32, id_b: u32) -> Option<&ProximityCacheEntry> {
        let (lo, hi) = if id_a <= id_b {
            (id_a, id_b)
        } else {
            (id_b, id_a)
        };
        self.entries.iter().find(|e| e.id_a == lo && e.id_b == hi)
    }

    /// Inserts or updates a closest-feature pair in the cache.
    ///
    /// When the cache is full the oldest entry (smallest `stamp`) is evicted.
    pub fn insert(
        &mut self,
        id_a: u32,
        id_b: u32,
        point_a: [f64; 3],
        point_b: [f64; 3],
        distance: f64,
        feature_a: ClosestFeature,
        feature_b: ClosestFeature,
        stamp: u64,
    ) {
        let (lo, hi, pa, pb, fa, fb) = if id_a <= id_b {
            (id_a, id_b, point_a, point_b, feature_a, feature_b)
        } else {
            (id_b, id_a, point_b, point_a, feature_b, feature_a)
        };

        // Update in place if already present
        if let Some(e) = self
            .entries
            .iter_mut()
            .find(|e| e.id_a == lo && e.id_b == hi)
        {
            e.point_a = pa;
            e.point_b = pb;
            e.distance = distance;
            e.feature_a = fa;
            e.feature_b = fb;
            e.stamp = stamp;
            return;
        }

        let entry = ProximityCacheEntry {
            id_a: lo,
            id_b: hi,
            point_a: pa,
            point_b: pb,
            distance,
            feature_a: fa,
            feature_b: fb,
            stamp,
        };

        if self.entries.len() < PROXIMITY_CACHE_CAPACITY {
            self.entries.push(entry);
        } else {
            // Evict oldest
            let oldest_idx = self
                .entries
                .iter()
                .enumerate()
                .min_by_key(|(_, e)| e.stamp)
                .map(|(i, _)| i)
                .unwrap_or(0);
            self.entries[oldest_idx] = entry;
        }
    }

    /// Removes the cached entry for pair `(id_a, id_b)` if present.
    pub fn remove(&mut self, id_a: u32, id_b: u32) {
        let (lo, hi) = if id_a <= id_b {
            (id_a, id_b)
        } else {
            (id_b, id_a)
        };
        self.entries.retain(|e| !(e.id_a == lo && e.id_b == hi));
    }

    /// Removes all entries whose stamp is strictly less than `min_stamp`.
    ///
    /// Call once per simulation step with the current frame index to evict
    /// stale entries from objects that have been removed.
    pub fn evict_stale(&mut self, min_stamp: u64) {
        self.entries.retain(|e| e.stamp >= min_stamp);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Triangle – triangle distance
// ─────────────────────────────────────────────────────────────────────────────

/// Result from a triangle–triangle minimum-distance query.
#[derive(Debug, Clone, Copy)]
pub struct TriTriResult {
    /// Closest point on the first triangle.
    pub point_a: [f64; 3],
    /// Closest point on the second triangle.
    pub point_b: [f64; 3],
    /// Minimum Euclidean distance between the two triangles.
    pub distance: f64,
}

/// Point at which the segment `p -> q` pierces the triangle `(v0, v1, v2)`.
///
/// The segment must cross the plane of the triangle (touching it with an end
/// point counts) at a point inside the triangle. A segment that lies in the
/// plane of the triangle never pierces it, and neither does any segment through
/// a zero-area triangle: such configurations are covered by the vertex-face and
/// edge-edge pairs of the distance search, which are exact for them.
fn segment_triangle_hit(
    p: [f64; 3],
    q: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> Option<[f64; 3]> {
    let normal = vcross(vsub(v1, v0), vsub(v2, v0));
    let normal_len_sq = vlen_sq(normal);
    if normal_len_sq <= 0.0 {
        return None;
    }
    // Signed distances (scaled by |normal|) of the end points from the plane.
    let side_p = vdot(normal, vsub(p, v0));
    let side_q = vdot(normal, vsub(q, v0));
    let same_side = (side_p > 0.0 && side_q > 0.0) || (side_p < 0.0 && side_q < 0.0);
    if same_side || (side_p == 0.0 && side_q == 0.0) {
        return None;
    }
    let t = side_p / (side_p - side_q);
    let hit = vadd(p, vscale(vsub(q, p), t));
    // The hit point is inside when it lies on the inner side of all three edges.
    let tolerance = -1e-12 * normal_len_sq;
    let inside = [(v0, v1), (v1, v2), (v2, v0)]
        .iter()
        .all(|&(a, b)| vdot(normal, vcross(vsub(b, a), vsub(hit, a))) >= tolerance);
    inside.then_some(hit)
}

/// Computes the minimum distance between two triangles in 3-D.
///
/// Intersecting triangles (crossing, touching or overlapping) have distance
/// zero: each of the six edges is tested against the other triangle, and the
/// result then carries the same witness point on both triangles, the centroid
/// of the piercing points (a point of the intersection, which is convex).
/// Otherwise the minimum is attained at a feature pair, and an exhaustive search
/// over vertex–triangle and edge–edge pairs finds it, following the approach
/// described in Real-Time Collision Detection §5.2. Coplanar triangles that
/// overlap are handled by that search as well (a vertex inside the other
/// triangle, or two crossing edges, give distance zero).
///
/// The triangles are defined by their three vertices `a0/a1/a2` and
/// `b0/b1/b2`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::triangle_triangle_distance;
/// let r = triangle_triangle_distance(
///     [0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],
///     [0.0,0.0,2.0],[1.0,0.0,2.0],[0.0,1.0,2.0]);
/// assert!((r.distance - 2.0).abs() < 1e-9);
/// ```
pub fn triangle_triangle_distance(
    a0: [f64; 3],
    a1: [f64; 3],
    a2: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
    b2: [f64; 3],
) -> TriTriResult {
    // Intersecting triangles: an edge of one pierces the other.
    let mut hit_sum = [0.0_f64; 3];
    let mut hit_count = 0.0_f64;
    for (p, q) in [(a0, a1), (a1, a2), (a2, a0)] {
        if let Some(hit) = segment_triangle_hit(p, q, b0, b1, b2) {
            hit_sum = vadd(hit_sum, hit);
            hit_count += 1.0;
        }
    }
    for (p, q) in [(b0, b1), (b1, b2), (b2, b0)] {
        if let Some(hit) = segment_triangle_hit(p, q, a0, a1, a2) {
            hit_sum = vadd(hit_sum, hit);
            hit_count += 1.0;
        }
    }
    if hit_count > 0.0 {
        let witness = vscale(hit_sum, 1.0 / hit_count);
        return TriTriResult {
            point_a: witness,
            point_b: witness,
            distance: 0.0,
        };
    }

    let mut best_dist = f64::MAX;
    let mut best_pa = a0;
    let mut best_pb = b0;

    // Helper: update best from a closest-point result against a fixed "other" point
    macro_rules! try_pt_tri {
        ($pt:expr, $q0:expr, $q1:expr, $q2:expr, $is_a:expr) => {{
            let r = point_triangle_closest($q0, $q1, $q2, $pt);
            if r.distance < best_dist {
                best_dist = r.distance;
                if $is_a {
                    best_pa = $pt;
                    best_pb = r.position;
                } else {
                    best_pa = r.position;
                    best_pb = $pt;
                }
            }
        }};
    }

    // Vertices of A against triangle B
    try_pt_tri!(a0, b0, b1, b2, true);
    try_pt_tri!(a1, b0, b1, b2, true);
    try_pt_tri!(a2, b0, b1, b2, true);

    // Vertices of B against triangle A
    try_pt_tri!(b0, a0, a1, a2, false);
    try_pt_tri!(b1, a0, a1, a2, false);
    try_pt_tri!(b2, a0, a1, a2, false);

    // Edge–edge pairs (3×3 = 9)
    let edges_a = [(a0, a1), (a1, a2), (a2, a0)];
    let edges_b = [(b0, b1), (b1, b2), (b2, b0)];
    for &(ea, eb) in &edges_a {
        for &(ec, ed) in &edges_b {
            let r = segment_segment_closest(ea, eb, ec, ed);
            if r.distance < best_dist {
                best_dist = r.distance;
                best_pa = r.point_a;
                best_pb = r.point_b;
            }
        }
    }

    TriTriResult {
        point_a: best_pa,
        point_b: best_pb,
        distance: best_dist,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Smooth-min (k-factor) CSG helper (bonus, not in spec but useful)
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth union of two SDF values using an exponential smooth-min kernel.
///
/// `k` controls the blending radius; larger `k` → more aggressive smoothing.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::sdf_smooth_union;
/// let v = sdf_smooth_union(0.0, 0.0, 32.0);
/// assert!(v <= 0.0);
/// ```
pub fn sdf_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    let h = (0.5 + 0.5 * (b - a) / k).clamp(0.0, 1.0);
    vlerp([b, 0.0, 0.0], [a, 0.0, 0.0], h)[0] - k * h * (1.0 - h)
}

// ─────────────────────────────────────────────────────────────────────────────
// Closest point on AABB
// ─────────────────────────────────────────────────────────────────────────────

/// Returns the closest point on an AABB (defined by `min` and `max`) to `p`.
///
/// When `p` is inside the box the closest point is `p` itself and the distance
/// is negative (use [`point_box_distance`] for the signed scalar).
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_box_closest;
/// let cp = point_box_closest([-1.0,-1.0,-1.0],[1.0,1.0,1.0],[3.0,0.0,0.0]);
/// assert!((cp[0] - 1.0).abs() < 1e-9);
/// ```
pub fn point_box_closest(min: [f64; 3], max: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        p[0].clamp(min[0], max[0]),
        p[1].clamp(min[1], max[1]),
        p[2].clamp(min[2], max[2]),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Point-in-triangle test (barycentric method)
// ─────────────────────────────────────────────────────────────────────────────

/// Tests whether point `p` (assumed co-planar with the triangle) lies inside
/// triangle `(v0, v1, v2)` using the barycentric method.
///
/// Returns `(u, v)` with `p = v0 + u * (v1 - v0) + v * (v2 - v0)` (so `u` and
/// `v` are the barycentric weights of `v1` and `v2`) if `p` lies inside or on
/// the boundary of the triangle, otherwise `None`. The coordinates are exact
/// for triangles of any size and orientation; degenerate (zero-area) triangles
/// return `None`.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_in_triangle_barycentric;
/// let r = point_in_triangle_barycentric(
///     [0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.25,0.25,0.0]);
/// assert!(r.is_some());
/// ```
pub fn point_in_triangle_barycentric(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    p: [f64; 3],
) -> Option<(f64, f64)> {
    let ab = vsub(v1, v0);
    let ac = vsub(v2, v0);
    let ap = vsub(p, v0);
    let n = vcross(ab, ac);
    // For a co-planar `p = v0 + u*ab + v*ac` the triple products satisfy
    // `(ap x ac) . n = u * |n|^2` and `(ab x ap) . n = v * |n|^2`, so the
    // normalisation must divide by `|n|^2` (not `|n|`) for the coordinates to
    // be independent of the triangle's size.
    let nn = vdot(n, n);
    if nn < 1e-300 {
        return None;
    }
    let u = vdot(vcross(ap, ac), n) / nn;
    let v = vdot(vcross(ab, ap), n) / nn;
    if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
        Some((u, v))
    } else {
        None
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ray – sphere / AABB distance helpers (complementary)
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance from point `p` to an infinite plane defined by `normal`
/// (unit vector) and scalar offset `d` (plane equation: dot(n, x) = d).
///
/// Positive on the side the normal points toward, negative on the other.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_plane_distance;
/// let d = point_plane_distance([0.0,1.0,0.0], 0.0, [0.0,2.0,0.0]);
/// assert!((d - 2.0).abs() < 1e-9);
/// ```
pub fn point_plane_distance(normal: [f64; 3], offset_d: f64, p: [f64; 3]) -> f64 {
    vdot(normal, p) - offset_d
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: squared distance point to segment
// ─────────────────────────────────────────────────────────────────────────────

/// Returns the squared distance from `p` to segment `[a, b]`.
///
/// This avoids a square-root and is suitable for comparisons.
///
/// # Examples
/// ```no_run
/// use oxiphysics_collision::proximity_query::point_segment_dist_sq;
/// assert!((point_segment_dist_sq([0.0,0.0,0.0],[2.0,0.0,0.0],[1.0,1.0,0.0]) - 1.0).abs() < 1e-9);
/// ```
pub fn point_segment_dist_sq(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> f64 {
    let res = point_segment_closest(a, b, p);
    res.distance * res.distance
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── point_segment_closest ────────────────────────────────────────────────

    #[test]
    fn test_point_segment_at_start() {
        let r = point_segment_closest([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!((r.u_param - 0.0).abs() < EPS, "t should be 0 at start");
        assert!((r.distance - 1.0).abs() < EPS);
        assert_eq!(r.feature, ClosestFeature::Vertex);
    }

    #[test]
    fn test_point_segment_at_end() {
        let r = point_segment_closest([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((r.u_param - 1.0).abs() < EPS, "t should be 1 at end");
        assert!((r.distance - 1.0).abs() < EPS);
        assert_eq!(r.feature, ClosestFeature::Vertex);
    }

    #[test]
    fn test_point_segment_midpoint() {
        let r = point_segment_closest([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]);
        assert!((r.u_param - 0.5).abs() < EPS);
        assert!((r.distance - 1.0).abs() < EPS);
        assert_eq!(r.feature, ClosestFeature::Edge);
    }

    #[test]
    fn test_point_segment_degenerate() {
        // Both endpoints coincide
        let r = point_segment_closest([3.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((r.distance - 3.0).abs() < EPS);
    }

    // ── point_triangle_closest ──────────────────────────────────────────────

    #[test]
    fn test_point_triangle_vertex_v0() {
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        );
        assert_eq!(r.feature, ClosestFeature::Vertex);
        assert!((r.position[0] - 0.0).abs() < EPS && (r.position[1] - 0.0).abs() < EPS);
    }

    #[test]
    fn test_point_triangle_vertex_v1() {
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
        );
        assert_eq!(r.feature, ClosestFeature::Vertex);
        assert!((r.position[0] - 2.0).abs() < EPS);
    }

    #[test]
    fn test_point_triangle_vertex_v2() {
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 4.0, 0.0],
        );
        assert_eq!(r.feature, ClosestFeature::Vertex);
        assert!((r.position[1] - 2.0).abs() < EPS);
    }

    #[test]
    fn test_point_triangle_edge_v0v1() {
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, -1.0, 0.0],
        );
        assert_eq!(r.feature, ClosestFeature::Edge);
        assert!((r.position[1] - 0.0).abs() < EPS);
    }

    #[test]
    fn test_point_triangle_face_interior() {
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert_eq!(r.feature, ClosestFeature::Face);
        assert!((r.distance).abs() < EPS);
    }

    #[test]
    fn test_point_triangle_above_face() {
        // Point directly above centroid
        let r = point_triangle_closest(
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [1.0, 1.0, 5.0],
        );
        assert!((r.distance - 5.0).abs() < EPS);
    }

    // ── segment_segment_closest ─────────────────────────────────────────────

    #[test]
    fn test_segseg_crossing() {
        // Two crossing segments in XY plane
        let r = segment_segment_closest(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!(r.distance < EPS, "crossing segments: dist={}", r.distance);
    }

    #[test]
    fn test_segseg_parallel() {
        // Parallel segments, no crossing
        let r = segment_segment_closest(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        );
        assert!(
            (r.distance - 1.0).abs() < EPS,
            "parallel dist should be 1.0"
        );
    }

    #[test]
    fn test_segseg_skew() {
        // Skew segments
        let r = segment_segment_closest(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 1.0],
            [0.5, 1.0, -1.0],
        );
        // Closest should be about 1.0 (vertical offset)
        assert!((r.distance - 1.0).abs() < EPS);
    }

    #[test]
    fn test_segseg_endpoint_to_endpoint() {
        let r = segment_segment_closest(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        );
        assert!((r.distance - 1.0).abs() < EPS);
        assert!((r.s - 1.0).abs() < EPS);
        assert!((r.t - 0.0).abs() < EPS);
    }

    // ── point_box_distance ──────────────────────────────────────────────────

    #[test]
    fn test_point_box_inside() {
        let d = point_box_distance([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(d < 0.0, "origin should be inside unit box");
        assert!((d - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_point_box_outside_face() {
        let d = point_box_distance([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [2.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < EPS);
    }

    #[test]
    fn test_point_box_outside_corner() {
        let d = point_box_distance([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let expected = (3.0_f64).sqrt();
        assert!((d - expected).abs() < EPS);
    }

    #[test]
    fn test_point_box_on_surface() {
        let d = point_box_distance([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [1.0, 0.0, 0.0]);
        assert!(
            d.abs() < EPS,
            "point on surface should have zero signed dist"
        );
    }

    // ── point_sphere_distance ───────────────────────────────────────────────

    #[test]
    fn test_point_sphere_inside() {
        let d = point_sphere_distance([0.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);
        assert!((d - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_point_sphere_outside() {
        let d = point_sphere_distance([0.0, 0.0, 0.0], 1.0, [3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < EPS);
    }

    #[test]
    fn test_point_sphere_on_surface() {
        let d = point_sphere_distance([1.0, 2.0, 3.0], 1.5, [2.5, 2.0, 3.0]);
        assert!(d.abs() < EPS);
    }

    // ── point_capsule_distance ──────────────────────────────────────────────

    #[test]
    fn test_point_capsule_inside_cylinder_portion() {
        let d = point_capsule_distance([0.0, 0.0, 0.0], [0.0, 4.0, 0.0], 1.0, [0.0, 2.0, 0.0]);
        assert!((d - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_point_capsule_outside_end_cap() {
        let d = point_capsule_distance([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5, [0.0, 3.0, 0.0]);
        // Nearest segment point is (0,1,0), dist=2.0, minus radius 0.5 → 1.5
        assert!((d - 1.5).abs() < EPS);
    }

    // ── SignedDistanceField ─────────────────────────────────────────────────

    #[test]
    fn test_sdf_sphere_inside() {
        let sdf = SignedDistanceField::Sphere {
            center: [0.0; 3],
            r: 2.0,
        };
        assert!(sdf.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_sdf_sphere_outside() {
        let sdf = SignedDistanceField::Sphere {
            center: [0.0; 3],
            r: 1.0,
        };
        let d = sdf.eval([2.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < EPS);
    }

    #[test]
    fn test_sdf_box_center() {
        let sdf = SignedDistanceField::Box {
            half: [1.0, 1.0, 1.0],
        };
        let d = sdf.eval([0.0, 0.0, 0.0]);
        assert!(d < 0.0);
        assert!((d - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_sdf_box_corner() {
        let sdf = SignedDistanceField::Box {
            half: [1.0, 1.0, 1.0],
        };
        let d = sdf.eval([2.0, 2.0, 2.0]);
        let expected = (3.0_f64).sqrt();
        assert!((d - expected).abs() < EPS);
    }

    #[test]
    fn test_sdf_box_face() {
        let sdf = SignedDistanceField::Box {
            half: [1.0, 1.0, 1.0],
        };
        let d = sdf.eval([3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < EPS);
    }

    #[test]
    fn test_sdf_capsule_inside() {
        let sdf = SignedDistanceField::Capsule {
            a: [0.0, 0.0, 0.0],
            b: [0.0, 2.0, 0.0],
            r: 1.0,
        };
        assert!(sdf.eval([0.0, 1.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_sdf_capsule_outside() {
        let sdf = SignedDistanceField::Capsule {
            a: [0.0, 0.0, 0.0],
            b: [0.0, 2.0, 0.0],
            r: 0.5,
        };
        let d = sdf.eval([0.0, 5.0, 0.0]);
        assert!((d - 2.5).abs() < EPS);
    }

    #[test]
    fn test_sdf_cylinder_inside() {
        let sdf = SignedDistanceField::Cylinder {
            axis: [0.0, 1.0, 0.0],
            r: 1.0,
            h: 4.0,
        };
        assert!(sdf.eval([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_sdf_cylinder_outside_radially() {
        let sdf = SignedDistanceField::Cylinder {
            axis: [0.0, 1.0, 0.0],
            r: 1.0,
            h: 4.0,
        };
        let d = sdf.eval([3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < EPS);
    }

    // ── CSG operations ──────────────────────────────────────────────────────

    #[test]
    fn test_csg_union_picks_smaller() {
        assert!((sdf_union(1.0, -0.5) - (-0.5)).abs() < EPS);
        assert!((sdf_union(-0.3, 0.7) - (-0.3)).abs() < EPS);
    }

    #[test]
    fn test_csg_intersect_picks_larger() {
        assert!((sdf_intersect(-1.0, 0.5) - 0.5).abs() < EPS);
        assert!((sdf_intersect(-0.5, -0.3) - (-0.3)).abs() < EPS);
    }

    #[test]
    fn test_csg_subtract() {
        // subtract positive region: should shrink inside
        assert!((sdf_subtract(-0.5, -0.3) - 0.3).abs() < EPS);
        assert!((sdf_subtract(-1.0, 2.0) - (-1.0)).abs() < EPS);
    }

    #[test]
    fn test_csg_two_spheres_union() {
        let s1 = SignedDistanceField::Sphere {
            center: [-2.0, 0.0, 0.0],
            r: 1.5,
        };
        let s2 = SignedDistanceField::Sphere {
            center: [2.0, 0.0, 0.0],
            r: 1.5,
        };
        // Point between: inside union if inside either
        let p = [0.0, 0.0, 0.0];
        let combined = sdf_union(s1.eval(p), s2.eval(p));
        // both spheres are 2 units away, radius 1.5 → outside, distance 0.5
        assert!((combined - 0.5).abs() < EPS);
    }

    // ── ProximityCache ──────────────────────────────────────────────────────

    #[test]
    fn test_proximity_cache_insert_and_get() {
        let mut cache = ProximityCache::new();
        cache.insert(
            1,
            2,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            ClosestFeature::Face,
            ClosestFeature::Face,
            0,
        );
        let e = cache.get(1, 2).expect("entry should exist");
        assert!((e.distance - 1.0).abs() < EPS);
    }

    #[test]
    fn test_proximity_cache_order_independent() {
        let mut cache = ProximityCache::new();
        cache.insert(
            5,
            3,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            ClosestFeature::Edge,
            ClosestFeature::Vertex,
            0,
        );
        assert!(cache.get(3, 5).is_some());
        assert!(cache.get(5, 3).is_some());
    }

    #[test]
    fn test_proximity_cache_remove() {
        let mut cache = ProximityCache::new();
        cache.insert(
            1,
            2,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            ClosestFeature::Face,
            ClosestFeature::Face,
            0,
        );
        cache.remove(1, 2);
        assert!(cache.get(1, 2).is_none());
    }

    #[test]
    fn test_proximity_cache_evict_stale() {
        let mut cache = ProximityCache::new();
        cache.insert(
            1,
            2,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            ClosestFeature::Face,
            ClosestFeature::Face,
            0,
        );
        cache.insert(
            3,
            4,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            ClosestFeature::Face,
            ClosestFeature::Face,
            10,
        );
        cache.evict_stale(5);
        assert!(cache.get(1, 2).is_none());
        assert!(cache.get(3, 4).is_some());
    }

    #[test]
    fn test_proximity_cache_capacity_evicts_oldest() {
        let mut cache = ProximityCache::new();
        for i in 0..PROXIMITY_CACHE_CAPACITY as u32 {
            cache.insert(
                i * 2,
                i * 2 + 1,
                [0.0; 3],
                [1.0, 0.0, 0.0],
                1.0,
                ClosestFeature::Face,
                ClosestFeature::Face,
                i as u64,
            );
        }
        assert_eq!(cache.len(), PROXIMITY_CACHE_CAPACITY);
        // Insert one more: oldest (stamp=0) should be evicted
        cache.insert(
            200,
            201,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            ClosestFeature::Face,
            ClosestFeature::Face,
            1000,
        );
        assert_eq!(cache.len(), PROXIMITY_CACHE_CAPACITY);
        assert!(cache.get(0, 1).is_none());
    }

    // ── triangle_triangle_distance ──────────────────────────────────────────

    #[test]
    fn test_tri_tri_parallel_offset() {
        let r = triangle_triangle_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 3.0],
            [1.0, 0.0, 3.0],
            [0.0, 1.0, 3.0],
        );
        assert!((r.distance - 3.0).abs() < EPS);
    }

    #[test]
    fn test_tri_tri_touching() {
        let r = triangle_triangle_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        );
        assert!(r.distance < EPS, "touching at origin, dist={}", r.distance);
    }

    // ── point_box_closest ───────────────────────────────────────────────────

    #[test]
    fn test_point_box_closest_outside() {
        let cp = point_box_closest([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [3.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < EPS);
        assert!((cp[1] - 0.0).abs() < EPS);
    }

    // ── sdf_smooth_union ────────────────────────────────────────────────────

    #[test]
    fn test_sdf_smooth_union_equals_min_when_far_apart() {
        // When values are far apart, smooth union ≈ regular union
        let a = -10.0_f64;
        let b = 10.0_f64;
        let su = sdf_smooth_union(a, b, 0.1);
        let u = sdf_union(a, b);
        assert!((su - u).abs() < 0.01, "smooth_union far apart should ≈ min");
    }

    // ── point_in_triangle_barycentric ───────────────────────────────────────

    #[test]
    fn test_point_in_triangle_barycentric_unit_triangle() {
        let r = point_in_triangle_barycentric(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.25, 0.25, 0.0],
        )
        .expect("interior point of the unit triangle");
        assert!((r.0 - 0.25).abs() < EPS, "u={}", r.0);
        assert!((r.1 - 0.25).abs() < EPS, "v={}", r.1);
    }

    #[test]
    fn test_point_in_triangle_barycentric_legs_two_triangle() {
        // |n| = 4: the buggy `/ |n|` normalisation returned u = 3 here.
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [2.0, 0.0, 0.0];
        let v2 = [0.0, 2.0, 0.0];
        let (u, v) = point_in_triangle_barycentric(v0, v1, v2, [1.5, 0.4, 0.0])
            .expect("(1.5, 0.4) lies inside the legs-2 triangle");
        assert!((u - 0.75).abs() < EPS, "u={u}");
        assert!((v - 0.2).abs() < EPS, "v={v}");
        assert!(
            point_in_triangle_barycentric(v0, v1, v2, [2.0, 2.0, 0.0]).is_none(),
            "(2, 2) is outside the triangle"
        );
        assert!(
            point_in_triangle_barycentric(v0, v1, v2, [-0.1, 0.5, 0.0]).is_none(),
            "(-0.1, 0.5) is outside the triangle"
        );
    }

    #[test]
    fn test_point_in_triangle_barycentric_is_scale_invariant() {
        for scale in [1.0e-3, 1.0, 1.0e3] {
            let v0 = [0.0, 0.0, 0.0];
            let v1 = [3.0 * scale, 0.0, 0.0];
            let v2 = [0.0, 3.0 * scale, 0.0];
            let centroid = [scale, scale, 0.0];
            let (u, v) = point_in_triangle_barycentric(v0, v1, v2, centroid)
                .unwrap_or_else(|| panic!("centroid must be inside at scale {scale}"));
            assert!((u - 1.0 / 3.0).abs() < 1e-9, "scale {scale}: u={u}");
            assert!((v - 1.0 / 3.0).abs() < 1e-9, "scale {scale}: v={v}");
        }
    }

    #[test]
    fn test_point_in_triangle_barycentric_tilted_triangle() {
        // Scalene triangle lying in the oblique plane x + y + z = 1.
        let v0 = [1.0, 0.0, 0.0];
        let v1 = [0.0, 2.0, -1.0];
        let v2 = [0.0, 0.0, 1.0];
        let (u, v) = (0.3, 0.45);
        let p = [
            v0[0] + u * (v1[0] - v0[0]) + v * (v2[0] - v0[0]),
            v0[1] + u * (v1[1] - v0[1]) + v * (v2[1] - v0[1]),
            v0[2] + u * (v1[2] - v0[2]) + v * (v2[2] - v0[2]),
        ];
        let (ru, rv) = point_in_triangle_barycentric(v0, v1, v2, p).expect("interior point");
        assert!((ru - u).abs() < 1e-12, "u={ru}");
        assert!((rv - v).abs() < 1e-12, "v={rv}");
    }

    #[test]
    fn test_point_in_triangle_barycentric_degenerate_triangle() {
        let r = point_in_triangle_barycentric(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [0.5, 0.5, 0.5],
        );
        assert!(r.is_none(), "collinear vertices span no triangle");
    }

    // ── triangle_triangle_distance: intersecting triangles ──────────────────

    type Tri = [[f64; 3]; 3];

    fn tri_distance(a: &Tri, b: &Tri) -> TriTriResult {
        triangle_triangle_distance(a[0], a[1], a[2], b[0], b[1], b[2])
    }

    /// Deterministic generator in `[0, 1)`.
    struct Lcg(u64);

    impl Lcg {
        fn next_unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        fn point(&mut self, extent: f64) -> [f64; 3] {
            [
                (self.next_unit() * 2.0 - 1.0) * extent,
                (self.next_unit() * 2.0 - 1.0) * extent,
                (self.next_unit() * 2.0 - 1.0) * extent,
            ]
        }
    }

    /// Distance between two triangles by alternating projections between the
    /// two convex sets: a monotone upper bound on the distance that converges
    /// to it. Independent of the feature-pair search under test.
    fn alternating_projection_distance(a: &Tri, b: &Tri) -> f64 {
        let mut p = vscale(vadd(vadd(a[0], a[1]), a[2]), 1.0 / 3.0);
        let mut best = f64::MAX;
        for _ in 0..50_000 {
            let q = point_triangle_closest(b[0], b[1], b[2], p).position;
            let next = point_triangle_closest(a[0], a[1], a[2], q).position;
            let d = vdist(next, q);
            if d > best - 1e-15 {
                best = best.min(d);
                break;
            }
            best = d;
            p = next;
        }
        best
    }

    #[test]
    fn test_tri_tri_piercing_edge_has_zero_distance() {
        // B stands upright in the plane x = 0 and its edges cross the plane of A
        // at (0, -0.5, 0) and (0, 0, 0), both inside A. No vertex of either
        // triangle lies on the other and no two edges meet, so the feature-pair
        // search alone would report a distance of about one.
        let a: Tri = [[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [0.0, 2.0, 0.0]];
        let b: Tri = [[0.0, -0.5, -1.0], [0.0, -0.5, 1.0], [0.0, 0.5, 1.0]];
        for (first, second) in [(&a, &b), (&b, &a)] {
            let r = tri_distance(first, second);
            assert_eq!(r.distance, 0.0, "pierced triangle");
            assert_eq!(r.point_a, r.point_b, "one witness for both triangles");
            // The witness lies on the segment where B crosses A.
            let w = r.point_a;
            assert!(w[0].abs() < EPS && w[2].abs() < EPS, "witness {w:?}");
            assert!(w[1] >= -0.5 - EPS && w[1] <= EPS, "witness {w:?}");
        }
    }

    #[test]
    fn test_tri_tri_crossed_x_shape_has_zero_distance() {
        // Two triangles crossing like an X, edges through each other's interior.
        let a: Tri = [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]];
        let b: Tri = [[-0.3, 0.0, -1.0], [0.3, 0.0, -1.0], [0.0, 0.0, 1.0]];
        let r = tri_distance(&a, &b);
        assert_eq!(r.distance, 0.0);
        assert!(r.point_a[2].abs() < EPS, "witness is in the plane of A");
        assert!(r.point_a[1].abs() < EPS, "witness is in the plane of B");
    }

    #[test]
    fn test_tri_tri_sharing_a_vertex_has_zero_distance() {
        let a: Tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let b: Tri = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]];
        let r = tri_distance(&a, &b);
        assert_eq!(r.distance, 0.0);
        assert!(vlen(r.point_a) < EPS);
    }

    #[test]
    fn test_tri_tri_coplanar_overlap_and_gap() {
        let a: Tri = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let overlapping: Tri = [[1.0, 0.5, 0.0], [3.0, 0.5, 0.0], [1.0, 2.5, 0.0]];
        assert!(tri_distance(&a, &overlapping).distance < EPS);
        let inside: Tri = [[0.2, 0.2, 0.0], [0.6, 0.2, 0.0], [0.2, 0.6, 0.0]];
        assert!(tri_distance(&a, &inside).distance < EPS);
        let apart: Tri = [[4.0, 0.0, 0.0], [5.0, 0.0, 0.0], [4.0, 1.0, 0.0]];
        let r = tri_distance(&a, &apart);
        assert!((r.distance - 2.0).abs() < EPS, "gap {}", r.distance);
        assert!((r.point_a[0] - 2.0).abs() < EPS && (r.point_b[0] - 4.0).abs() < EPS);
    }

    #[test]
    fn test_tri_tri_separated_triangle_above_the_plane_keeps_its_distance() {
        // B passes above A and does not reach its plane: a segment test that
        // ignored the end-point sides would report an intersection.
        let a: Tri = [[-2.0, -2.0, 0.0], [2.0, -2.0, 0.0], [0.0, 2.0, 0.0]];
        let b: Tri = [[0.0, -0.5, 0.25], [0.0, -0.5, 1.0], [0.0, 0.5, 1.0]];
        let r = tri_distance(&a, &b);
        assert!((r.distance - 0.25).abs() < EPS, "distance {}", r.distance);
        // B touches the plane of A outside of A: still a positive distance.
        let outside: Tri = [[5.0, 0.0, -1.0], [5.0, 0.0, 1.0], [5.0, 1.0, 1.0]];
        assert!(tri_distance(&a, &outside).distance > 2.9);
    }

    #[test]
    fn test_tri_tri_matches_alternating_projections_on_random_pairs() {
        let mut rng = Lcg(0x7121_9A55_0000_0042);
        let mut intersecting = 0;
        for case in 0..300 {
            let a: Tri = [rng.point(1.0), rng.point(1.0), rng.point(1.0)];
            let shift = rng.point(0.6);
            let b: Tri = [
                vadd(rng.point(1.0), shift),
                vadd(rng.point(1.0), shift),
                vadd(rng.point(1.0), shift),
            ];
            let exact = tri_distance(&a, &b);
            let reference = alternating_projection_distance(&a, &b);
            assert!(
                exact.distance <= reference + 1e-9,
                "case {case}: exact {} exceeds the reference upper bound {reference}",
                exact.distance
            );
            if exact.distance == 0.0 {
                intersecting += 1;
                assert!(
                    reference < 1e-3,
                    "case {case}: reported as intersecting but the triangles are {reference} apart"
                );
                assert_eq!(exact.point_a, exact.point_b);
            } else {
                assert!(
                    (exact.distance - vdist(exact.point_a, exact.point_b)).abs() < 1e-9,
                    "case {case}: witnesses must realise the distance"
                );
                assert!(
                    (reference - exact.distance).abs() < 1e-6,
                    "case {case}: exact {} vs reference {reference}",
                    exact.distance
                );
            }
        }
        assert!(
            (20..=200).contains(&intersecting),
            "the sample must contain both intersecting and separated pairs, got {intersecting} intersecting"
        );
    }
}
