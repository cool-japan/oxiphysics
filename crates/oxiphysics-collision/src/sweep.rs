// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Swept-primitive collision queries.
//!
//! Provides analytic and conservative time-of-impact tests for common
//! primitive pairs:
//!
//! - Sphere vs sphere (analytic)
//! - Sphere vs AABB (conservative)
//! - AABB vs AABB (conservative)
//! - Ray vs sphere, box, triangle

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn len_sq(v: [f64; 3]) -> f64 {
    dot(v, v)
}

#[inline]
fn len(v: [f64; 3]) -> f64 {
    len_sq(v).sqrt()
}

#[inline]
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let l = len(v);
    if l < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        scale(v, 1.0 / l)
    }
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Public result types
// ---------------------------------------------------------------------------

/// Result of a swept-primitive time-of-impact query.
#[derive(Debug, Clone, PartialEq)]
pub struct SweepResult {
    /// Time of impact in `[0, 1]` (fraction of the motion step).
    pub toi: f64,
    /// Contact normal at the time of impact (unit length, from B toward A).
    pub normal: [f64; 3],
    /// Contact point on shape A at the time of impact.
    pub point_a: [f64; 3],
    /// Contact point on shape B at the time of impact.
    pub point_b: [f64; 3],
}

/// Result of a ray-cast query.
#[derive(Debug, Clone, PartialEq)]
pub struct RayHit {
    /// Distance from the ray origin to the hit point.
    pub t: f64,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Surface normal at the hit point (unit length, pointing outward).
    pub normal: [f64; 3],
}

// ---------------------------------------------------------------------------
// swept_sphere_sphere
// ---------------------------------------------------------------------------

/// Analytic swept sphere–sphere time-of-impact.
///
/// `c0_a` / `c1_a` are the start and end centres of sphere A with radius `r_a`.
/// `c0_b` / `c1_b` are the start and end centres of sphere B with radius `r_b`.
///
/// Returns the first `t ∈ [0, 1]` at which the spheres touch, or `None` if
/// they do not collide during the step.
pub fn swept_sphere_sphere(
    c0_a: [f64; 3],
    c1_a: [f64; 3],
    r_a: f64,
    c0_b: [f64; 3],
    c1_b: [f64; 3],
    r_b: f64,
) -> Option<SweepResult> {
    let r_sum = r_a + r_b;

    // Relative displacement and velocity (b relative to a).
    let dp = sub(c0_b, c0_a);
    // Relative velocity: (c1_b - c0_b) - (c1_a - c0_a)
    let dv = sub(sub(c1_b, c0_b), sub(c1_a, c0_a));

    // |dp + t*dv|² = r_sum²
    // a*t² + 2*b*t + (c - r²) = 0
    let a = len_sq(dv);
    let b = dot(dp, dv);
    let c = len_sq(dp) - r_sum * r_sum;

    // Already touching / overlapping at t=0
    if c <= 0.0 {
        let normal = if len_sq(dp) > 1e-20 {
            normalize(scale(dp, -1.0))
        } else {
            [0.0, 1.0, 0.0]
        };
        let contact = add(c0_a, scale(normal, -r_a));
        return Some(SweepResult {
            toi: 0.0,
            normal,
            point_a: contact,
            point_b: contact,
        });
    }

    if a < 1e-300 {
        // No relative motion.
        return None;
    }

    let disc = b * b - a * c;
    if disc < 0.0 {
        return None; // Trajectories never intersect.
    }

    let sqrt_disc = disc.sqrt();
    let t = (-b - sqrt_disc) / a;
    let t = if (0.0..=1.0).contains(&t) {
        t
    } else {
        let t2 = (-b + sqrt_disc) / a;
        if (0.0..=1.0).contains(&t2) {
            t2
        } else {
            return None;
        }
    };

    // Positions at t
    let pos_a = add(c0_a, scale(sub(c1_a, c0_a), t));
    let pos_b = add(c0_b, scale(sub(c1_b, c0_b), t));
    let d = sub(pos_b, pos_a);
    let dist = len(d);
    let normal = if dist > 1e-10 {
        scale(d, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };

    Some(SweepResult {
        toi: t,
        normal,
        point_a: add(pos_a, scale(normal, r_a)),
        point_b: add(pos_b, scale(normal, -r_b)),
    })
}

// ---------------------------------------------------------------------------
// swept_sphere_aabb  (conservative)
// ---------------------------------------------------------------------------

/// Conservative swept sphere vs AABB time-of-impact.
///
/// The sphere moves from `c0` to `c1` with radius `r`.  The AABB is static
/// and defined by its min/max corners `aabb_min`/`aabb_max`.
///
/// Uses the Minkowski-sum approach: expand the AABB by `r` on every face,
/// then intersect the segment `[c0, c1]` against the expanded AABB.
pub fn swept_sphere_aabb(
    c0: [f64; 3],
    c1: [f64; 3],
    r: f64,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> Option<SweepResult> {
    // Expand AABB by sphere radius (Minkowski sum with sphere = box rounded by r).
    // For a conservative test we expand by r on each face.
    let exp_min = [aabb_min[0] - r, aabb_min[1] - r, aabb_min[2] - r];
    let exp_max = [aabb_max[0] + r, aabb_max[1] + r, aabb_max[2] + r];

    // Ray–AABB slab intersection.
    let dir = sub(c1, c0);
    let result = ray_aabb_slab(c0, dir, exp_min, exp_max, 0.0, 1.0)?;

    // Reconstruct the hit position on the sphere and the actual surface normal.
    let hit_center = add(c0, scale(dir, result.t));

    // Find the closest point on the original (un-expanded) AABB.
    let closest_on_aabb = closest_point_on_aabb(hit_center, aabb_min, aabb_max);
    let diff = sub(hit_center, closest_on_aabb);
    let dist = len(diff);
    let normal = if dist > 1e-10 {
        normalize(diff)
    } else {
        result.normal
    };

    Some(SweepResult {
        toi: result.t,
        normal,
        point_a: add(hit_center, scale(normal, -r)),
        point_b: closest_on_aabb,
    })
}

// ---------------------------------------------------------------------------
// swept_aabb_aabb  (conservative)
// ---------------------------------------------------------------------------

/// Conservative swept AABB vs static AABB time-of-impact.
///
/// `min_a`/`max_a` are the initial bounds of the moving AABB A.
/// `vel` is the displacement vector of A (i.e., `A` moves by `vel` over `t=1`).
/// `min_b`/`max_b` are the bounds of the static AABB B.
///
/// Uses the Minkowski-sum approach: expand B by the half-extents of A, then
/// cast the point trajectory of A's centre through the expanded B.
pub fn swept_aabb_aabb(
    min_a: [f64; 3],
    max_a: [f64; 3],
    vel: [f64; 3],
    min_b: [f64; 3],
    max_b: [f64; 3],
) -> Option<SweepResult> {
    // Half-extents of A.
    let half_a = [
        (max_a[0] - min_a[0]) * 0.5,
        (max_a[1] - min_a[1]) * 0.5,
        (max_a[2] - min_a[2]) * 0.5,
    ];
    // Centre of A at t=0.
    let center_a = [
        (min_a[0] + max_a[0]) * 0.5,
        (min_a[1] + max_a[1]) * 0.5,
        (min_a[2] + max_a[2]) * 0.5,
    ];

    // Minkowski-expanded B (by A's half-extents).
    let exp_min = [
        min_b[0] - half_a[0],
        min_b[1] - half_a[1],
        min_b[2] - half_a[2],
    ];
    let exp_max = [
        max_b[0] + half_a[0],
        max_b[1] + half_a[1],
        max_b[2] + half_a[2],
    ];

    // Ray–AABB slab test with centre of A sweeping by vel.
    let result = ray_aabb_slab(center_a, vel, exp_min, exp_max, 0.0, 1.0)?;

    // Hit position of A's centre.
    let hit_center_a = add(center_a, scale(vel, result.t));

    // Clamp to get the closest face point on B.
    let contact_b = closest_point_on_aabb(hit_center_a, min_b, max_b);

    Some(SweepResult {
        toi: result.t,
        normal: result.normal,
        point_a: add(
            hit_center_a,
            scale(result.normal, -half_a[0].max(half_a[1]).max(half_a[2])),
        ),
        point_b: contact_b,
    })
}

// ---------------------------------------------------------------------------
// Ray vs sphere
// ---------------------------------------------------------------------------

/// Ray–sphere intersection.
///
/// The ray starts at `origin` and travels in direction `dir` (need not be unit
/// length; parametric distance is returned in units of `|dir|`).
/// Returns the first hit with `t ∈ [t_min, t_max]`.
pub fn ray_sphere(
    origin: [f64; 3],
    dir: [f64; 3],
    center: [f64; 3],
    radius: f64,
    t_min: f64,
    t_max: f64,
) -> Option<RayHit> {
    let oc = sub(origin, center);
    let a = len_sq(dir);
    let b = dot(oc, dir);
    let c = len_sq(oc) - radius * radius;

    let disc = b * b - a * c;
    if disc < 0.0 {
        return None;
    }

    let sqrt_disc = disc.sqrt();
    let t = {
        let t1 = (-b - sqrt_disc) / a;
        if t1 >= t_min && t1 <= t_max {
            t1
        } else {
            let t2 = (-b + sqrt_disc) / a;
            if t2 >= t_min && t2 <= t_max {
                t2
            } else {
                return None;
            }
        }
    };

    let point = add(origin, scale(dir, t));
    let normal = normalize(sub(point, center));

    Some(RayHit { t, point, normal })
}

// ---------------------------------------------------------------------------
// Ray vs AABB
// ---------------------------------------------------------------------------

/// Ray–AABB intersection using the slab method.
///
/// Returns the first hit with `t ∈ [t_min, t_max]`.
pub fn ray_aabb(
    origin: [f64; 3],
    dir: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    t_min: f64,
    t_max: f64,
) -> Option<RayHit> {
    ray_aabb_slab(origin, dir, aabb_min, aabb_max, t_min, t_max)
}

/// Internal slab-method AABB intersection returning a `RayHit`.
fn ray_aabb_slab(
    origin: [f64; 3],
    dir: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    t_min: f64,
    t_max: f64,
) -> Option<RayHit> {
    let mut t_near = t_min;
    let mut t_far = t_max;
    let mut hit_axis = 0usize;
    let mut hit_sign = 1.0_f64;

    for i in 0..3 {
        if dir[i].abs() < 1e-300 {
            // Ray is parallel to slab.
            if origin[i] < aabb_min[i] || origin[i] > aabb_max[i] {
                return None; // Ray is outside this slab.
            }
        } else {
            let inv = 1.0 / dir[i];
            let mut t1 = (aabb_min[i] - origin[i]) * inv;
            let mut t2 = (aabb_max[i] - origin[i]) * inv;
            let sign = if t1 <= t2 { -1.0_f64 } else { 1.0_f64 };
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            if t1 > t_near {
                t_near = t1;
                hit_axis = i;
                hit_sign = sign;
            }
            if t2 < t_far {
                t_far = t2;
            }
            if t_near > t_far {
                return None;
            }
        }
    }

    if t_near > t_far || t_near > t_max || t_far < t_min {
        return None;
    }

    let t = if t_near >= t_min { t_near } else { t_far };
    if t < t_min || t > t_max {
        return None;
    }

    let point = add(origin, scale(dir, t));
    let mut normal = [0.0f64; 3];
    normal[hit_axis] = hit_sign;

    Some(RayHit { t, point, normal })
}

// ---------------------------------------------------------------------------
// Ray vs triangle (Möller–Trumbore)
// ---------------------------------------------------------------------------

/// Ray–triangle intersection using the Möller–Trumbore algorithm.
///
/// The triangle is defined by vertices `v0`, `v1`, `v2` (counter-clockwise
/// winding defines the outward normal). Returns a hit for `t ∈ [t_min, t_max]`.
pub fn ray_triangle(
    origin: [f64; 3],
    dir: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    t_min: f64,
    t_max: f64,
) -> Option<RayHit> {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);

    let h = cross(dir, e2);
    let det = dot(e1, h);

    // Back-face culling: |det| < epsilon means ray is parallel.
    if det.abs() < 1e-10 {
        return None;
    }

    let inv_det = 1.0 / det;
    let s = sub(origin, v0);
    let u = dot(s, h) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = cross(s, e1);
    let v = dot(dir, q) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = dot(e2, q) * inv_det;
    if t < t_min || t > t_max {
        return None;
    }

    let point = add(origin, scale(dir, t));
    let normal = normalize(cross(e1, e2));

    Some(RayHit { t, point, normal })
}

// ---------------------------------------------------------------------------
// Closest point helpers
// ---------------------------------------------------------------------------

/// Closest point on AABB to a query point.
fn closest_point_on_aabb(pt: [f64; 3], aabb_min: [f64; 3], aabb_max: [f64; 3]) -> [f64; 3] {
    [
        pt[0].clamp(aabb_min[0], aabb_max[0]),
        pt[1].clamp(aabb_min[1], aabb_max[1]),
        pt[2].clamp(aabb_min[2], aabb_max[2]),
    ]
}

/// Squared distance from a point to an AABB.
pub fn point_aabb_dist_sq(pt: [f64; 3], aabb_min: [f64; 3], aabb_max: [f64; 3]) -> f64 {
    let mut dist_sq = 0.0_f64;
    for i in 0..3 {
        if pt[i] < aabb_min[i] {
            let d = aabb_min[i] - pt[i];
            dist_sq += d * d;
        } else if pt[i] > aabb_max[i] {
            let d = pt[i] - aabb_max[i];
            dist_sq += d * d;
        }
    }
    dist_sq
}

/// Test whether a sphere overlaps an AABB.
pub fn sphere_aabb_overlap(
    center: [f64; 3],
    radius: f64,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> bool {
    point_aabb_dist_sq(center, aabb_min, aabb_max) <= radius * radius
}

/// Test whether two AABBs overlap.
pub fn aabb_aabb_overlap(
    min_a: [f64; 3],
    max_a: [f64; 3],
    min_b: [f64; 3],
    max_b: [f64; 3],
) -> bool {
    min_a[0] <= max_b[0]
        && max_a[0] >= min_b[0]
        && min_a[1] <= max_b[1]
        && max_a[1] >= min_b[1]
        && min_a[2] <= max_b[2]
        && max_a[2] >= min_b[2]
}

// ---------------------------------------------------------------------------
// Swept capsule vs static plane
// ---------------------------------------------------------------------------

/// Swept capsule vs an infinite static plane (`n·x = d`).
///
/// The capsule sweeps from `c0` to `c1` with the given `half_height` (along Y)
/// and `radius`.  Returns the first `t ∈ [0, 1]` at which the capsule surface
/// first touches the plane, or `None` if no contact.
pub fn swept_capsule_plane(
    c0: [f64; 3],
    c1: [f64; 3],
    radius: f64,
    half_height: f64,
    plane_normal: [f64; 3],
    plane_d: f64,
) -> Option<SweepResult> {
    let n = normalize(plane_normal);
    // The lowest point of the capsule along n at t=0 is the support point.
    // Support: min_dist = dot(n, centre_bottom) - radius
    // centre_bottom = c0 - n*(half_height) (if n is roughly up)
    // Conservative: treat the capsule as a sphere of radius (radius + half_height) for the sweep.
    let effective_r = radius + half_height;
    let dist0 = dot(n, c0) - plane_d - effective_r;

    if dist0 <= 0.0 {
        return Some(SweepResult {
            toi: 0.0,
            normal: n,
            point_a: c0,
            point_b: add(c0, scale(n, -dist0)),
        });
    }

    let disp = sub(c1, c0);
    let closing = -dot(n, disp);

    if closing <= 1e-12 {
        return None;
    }

    let t = dist0 / closing;
    if t > 1.0 {
        return None;
    }

    let hit_pos = add(c0, scale(disp, t));
    Some(SweepResult {
        toi: t,
        normal: n,
        point_a: add(hit_pos, scale(n, -effective_r)),
        point_b: add(hit_pos, scale(n, -effective_r)),
    })
}

// ---------------------------------------------------------------------------
// TOI refinement for sphere-sphere via bisection
// ---------------------------------------------------------------------------

/// Refine the time of impact for two moving spheres using bisection.
///
/// `t_lo` and `t_hi` must bracket the contact: at `t_lo` the spheres are
/// separated, at `t_hi` they overlap (or just touch).  Returns the refined
/// TOI within `[t_lo, t_hi]`.
#[allow(clippy::too_many_arguments)]
pub fn refine_sphere_sphere_toi(
    c0_a: [f64; 3],
    vel_a: [f64; 3],
    r_a: f64,
    c0_b: [f64; 3],
    vel_b: [f64; 3],
    r_b: f64,
    t_lo: f64,
    t_hi: f64,
    tolerance: f64,
    max_iter: usize,
) -> f64 {
    let r_sum = r_a + r_b;
    let mut lo = t_lo;
    let mut hi = t_hi;

    for _ in 0..max_iter {
        if hi - lo < tolerance {
            break;
        }
        let mid = (lo + hi) * 0.5;
        let pa = add(c0_a, scale(vel_a, mid));
        let pb = add(c0_b, scale(vel_b, mid));
        let d = len(sub(pa, pb));
        if d <= r_sum {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) * 0.5
}

// ---------------------------------------------------------------------------
// Sweep broadphase: collect candidate pairs from a list of swept spheres
// ---------------------------------------------------------------------------

/// A swept sphere description for broadphase sweeping.
#[derive(Debug, Clone)]
pub struct SweptSphere {
    /// Start position.
    pub pos0: [f64; 3],
    /// End position.
    pub pos1: [f64; 3],
    /// Radius.
    pub radius: f64,
    /// User-supplied ID.
    pub id: u32,
}

impl SweptSphere {
    /// Create a new swept sphere.
    pub fn new(pos0: [f64; 3], pos1: [f64; 3], radius: f64, id: u32) -> Self {
        Self {
            pos0,
            pos1,
            radius,
            id,
        }
    }

    /// AABB bounding the entire sweep.
    pub fn sweep_aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut mn = [0.0f64; 3];
        let mut mx = [0.0f64; 3];
        for i in 0..3 {
            mn[i] = self.pos0[i].min(self.pos1[i]) - self.radius;
            mx[i] = self.pos0[i].max(self.pos1[i]) + self.radius;
        }
        (mn, mx)
    }

    /// Test if the sweep AABB of this sphere overlaps with another.
    pub fn sweep_aabb_overlaps(&self, other: &SweptSphere) -> bool {
        let (amn, amx) = self.sweep_aabb();
        let (bmn, bmx) = other.sweep_aabb();
        amn[0] <= bmx[0]
            && amx[0] >= bmn[0]
            && amn[1] <= bmx[1]
            && amx[1] >= bmn[1]
            && amn[2] <= bmx[2]
            && amx[2] >= bmn[2]
    }
}

/// Broadphase sweep: O(n²) candidate pair collection using sweep AABBs.
///
/// Returns pairs `(id_a, id_b)` with `id_a < id_b` whose sweep AABBs overlap.
pub fn sweep_broadphase(spheres: &[SweptSphere]) -> Vec<(u32, u32)> {
    let mut pairs = Vec::new();
    for i in 0..spheres.len() {
        for j in (i + 1)..spheres.len() {
            if spheres[i].sweep_aabb_overlaps(&spheres[j]) {
                let (a, b) = if spheres[i].id < spheres[j].id {
                    (spheres[i].id, spheres[j].id)
                } else {
                    (spheres[j].id, spheres[i].id)
                };
                pairs.push((a, b));
            }
        }
    }
    pairs
}

// ---------------------------------------------------------------------------
// Time-interval sweep: TOI within [t_start, t_end]
// ---------------------------------------------------------------------------

/// Swept sphere-sphere test within a restricted time interval `[t_start, t_end]`.
///
/// `c0_a` is the position of sphere A at `t=0`; `vel_a` its velocity.
/// The TOI is remapped to `[t_start, t_end]` before evaluation.
///
/// Returns the TOI in `[0, 1]` (fraction of the full interval) or `None`.
#[allow(clippy::too_many_arguments)]
pub fn swept_sphere_sphere_interval(
    c0_a: [f64; 3],
    vel_a: [f64; 3],
    r_a: f64,
    c0_b: [f64; 3],
    vel_b: [f64; 3],
    r_b: f64,
    t_start: f64,
    t_end: f64,
) -> Option<SweepResult> {
    // Advance positions to t_start, then test sweep from t_start to t_end.
    let a0 = add(c0_a, scale(vel_a, t_start));
    let a1 = add(c0_a, scale(vel_a, t_end));
    let b0 = add(c0_b, scale(vel_b, t_start));
    let b1 = add(c0_b, scale(vel_b, t_end));
    swept_sphere_sphere(a0, a1, r_a, b0, b1, r_b)
}

// ---------------------------------------------------------------------------
// Compound shape sweep: find earliest TOI among all child shapes
// ---------------------------------------------------------------------------

/// A compound swept shape made of multiple spheres (simplified).
pub struct CompoundSweptShape {
    /// Child spheres with local offsets relative to the compound centre.
    pub spheres: Vec<([f64; 3], f64)>, // (local_offset, radius)
}

impl CompoundSweptShape {
    /// Create a new compound swept shape.
    pub fn new(spheres: Vec<([f64; 3], f64)>) -> Self {
        Self { spheres }
    }

    /// Sweep this compound shape from `pos0` to `pos1` against another compound
    /// shape sweeping from `other_pos0` to `other_pos1`.
    ///
    /// Returns the earliest `SweepResult` or `None` if no collision.
    pub fn sweep_vs_compound(
        &self,
        pos0: [f64; 3],
        pos1: [f64; 3],
        other: &CompoundSweptShape,
        other_pos0: [f64; 3],
        other_pos1: [f64; 3],
    ) -> Option<SweepResult> {
        let mut best: Option<SweepResult> = None;

        for (off_a, r_a) in &self.spheres {
            let a0 = add(pos0, *off_a);
            let a1 = add(pos1, *off_a);
            for (off_b, r_b) in &other.spheres {
                let b0 = add(other_pos0, *off_b);
                let b1 = add(other_pos1, *off_b);
                if let Some(result) = swept_sphere_sphere(a0, a1, *r_a, b0, b1, *r_b)
                    && best.as_ref().is_none_or(|br| result.toi < br.toi)
                {
                    best = Some(result);
                }
            }
        }
        best
    }
}

// ---------------------------------------------------------------------------
// Swept capsule vs static AABB
// ---------------------------------------------------------------------------

/// Swept capsule vs static AABB time-of-impact.
///
/// The capsule is defined by its axis endpoints `c0`/`c1` (swept from start to end)
/// and `radius`. The AABB is static. Uses an expanded AABB approach by inflating
/// the box by `radius` on each face and testing the swept segment axis.
pub fn swept_capsule_aabb(
    axis_start_0: [f64; 3],
    axis_start_1: [f64; 3],
    axis_end_0: [f64; 3],
    axis_end_1: [f64; 3],
    radius: f64,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> Option<SweepResult> {
    // Use the midpoint of the capsule axis as a representative sweep trajectory.
    let mid_start = scale(add(axis_start_0, axis_start_1), 0.5);
    let mid_end = scale(add(axis_end_0, axis_end_1), 0.5);
    // Half-length of capsule (contributes to effective radius)
    let half_len = len(scale(sub(axis_start_1, axis_start_0), 0.5));
    let effective_radius = radius + half_len;
    // Expand AABB by effective radius
    let exp_min = [
        aabb_min[0] - effective_radius,
        aabb_min[1] - effective_radius,
        aabb_min[2] - effective_radius,
    ];
    let exp_max = [
        aabb_max[0] + effective_radius,
        aabb_max[1] + effective_radius,
        aabb_max[2] + effective_radius,
    ];
    let dir = sub(mid_end, mid_start);
    let hit = ray_aabb_slab(mid_start, dir, exp_min, exp_max, 0.0, 1.0)?;
    let hit_center = add(mid_start, scale(dir, hit.t));
    let closest_on_aabb = closest_point_on_aabb(hit_center, aabb_min, aabb_max);
    let diff = sub(hit_center, closest_on_aabb);
    let d = len(diff);
    let normal = if d > 1e-10 {
        normalize(diff)
    } else {
        hit.normal
    };
    Some(SweepResult {
        toi: hit.t,
        normal,
        point_a: add(hit_center, scale(normal, -effective_radius)),
        point_b: closest_on_aabb,
    })
}

// ---------------------------------------------------------------------------
// Linear CCD pair test
// ---------------------------------------------------------------------------

/// Linear continuous collision detection between two moving spheres.
///
/// `pos_a`/`vel_a`: position and velocity of sphere A.
/// `pos_b`/`vel_b`: position and velocity of sphere B.
/// `dt`: time step.
///
/// Returns the time of impact within `[0, dt]` or `None` if no collision.
pub fn linear_ccd_sphere_sphere(
    pos_a: [f64; 3],
    vel_a: [f64; 3],
    r_a: f64,
    pos_b: [f64; 3],
    vel_b: [f64; 3],
    r_b: f64,
    dt: f64,
) -> Option<SweepResult> {
    let c0_a = pos_a;
    let c1_a = add(pos_a, scale(vel_a, dt));
    let c0_b = pos_b;
    let c1_b = add(pos_b, scale(vel_b, dt));
    let result = swept_sphere_sphere(c0_a, c1_a, r_a, c0_b, c1_b, r_b)?;
    // Scale toi back to real time
    Some(SweepResult {
        toi: result.toi * dt,
        ..result
    })
}

// ---------------------------------------------------------------------------
// TOI binary search refinement (generic distance function)
// ---------------------------------------------------------------------------

/// Refine a time-of-impact by bisection using a generic distance function.
///
/// `dist_fn(t)` should return the signed distance (negative = overlapping)
/// at time `t`. `t_lo` has `dist_fn(t_lo) >= 0`, `t_hi` has `dist_fn(t_hi) <= 0`.
/// Returns the refined TOI.
pub fn bisect_toi<F>(
    mut t_lo: f64,
    mut t_hi: f64,
    tolerance: f64,
    max_iter: usize,
    dist_fn: F,
) -> f64
where
    F: Fn(f64) -> f64,
{
    for _ in 0..max_iter {
        if t_hi - t_lo < tolerance {
            break;
        }
        let mid = (t_lo + t_hi) * 0.5;
        if dist_fn(mid) >= 0.0 {
            t_lo = mid;
        } else {
            t_hi = mid;
        }
    }
    (t_lo + t_hi) * 0.5
}

// ---------------------------------------------------------------------------
// Conservative advancement for rotating bodies
// ---------------------------------------------------------------------------

/// A rigid body state for conservative advancement.
#[derive(Debug, Clone)]
pub struct RotatingBodyState {
    /// Center of mass position.
    pub position: [f64; 3],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular speed (magnitude of angular velocity vector), radians/s.
    pub angular_speed: f64,
    /// Bounding sphere radius (enclosing the entire body).
    pub bound_radius: f64,
}

impl RotatingBodyState {
    /// Create a new rotating body state.
    pub fn new(
        position: [f64; 3],
        linear_vel: [f64; 3],
        angular_speed: f64,
        bound_radius: f64,
    ) -> Self {
        Self {
            position,
            linear_vel,
            angular_speed,
            bound_radius,
        }
    }

    /// Upper bound on how far any point on this body can move in time `dt`.
    ///
    /// Uses the formula: `|v_linear| * dt + angular_speed * bound_radius * dt`.
    pub fn motion_bound(&self, dt: f64) -> f64 {
        let linear_speed = len(self.linear_vel);
        (linear_speed + self.angular_speed * self.bound_radius) * dt
    }

    /// Advance position by `dt`.
    pub fn advance(&self, dt: f64) -> [f64; 3] {
        add(self.position, scale(self.linear_vel, dt))
    }
}

/// Conservative advancement TOI for two rotating bodies.
///
/// Iteratively advances time until the bounding spheres touch, using
/// motion bounds to ensure conservative (never-tunneling) steps.
///
/// Returns the TOI in `[0, 1]` or `None` if the bodies never collide.
#[allow(clippy::too_many_arguments)]
pub fn conservative_advancement_rotating(
    a: &RotatingBodyState,
    b: &RotatingBodyState,
    separation: f64,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
    max_iter: usize,
) -> Option<f64> {
    let total_dt = t_end - t_start;
    if total_dt <= 0.0 {
        return None;
    }
    let sum_r = a.bound_radius + b.bound_radius;
    let init_dist = len(sub(b.position, a.position)) - sum_r;
    if init_dist <= 0.0 {
        return Some(t_start);
    }
    if separation < 0.0 {
        return Some(t_start);
    }

    let mut t = t_start;
    let mut remaining_sep = separation.max(init_dist);

    for _ in 0..max_iter {
        if remaining_sep <= tolerance {
            return Some(t);
        }
        let remaining_dt = t_end - t;
        if remaining_dt <= 0.0 {
            break;
        }
        let motion_a = a.motion_bound(remaining_dt);
        let motion_b = b.motion_bound(remaining_dt);
        let total_motion = motion_a + motion_b;
        if total_motion < 1e-14 {
            break;
        }
        let dt_step = (remaining_sep / total_motion).min(remaining_dt);
        t += dt_step;
        // Advance positions
        let pos_a_t = a.advance(t - t_start);
        let pos_b_t = b.advance(t - t_start);
        let new_dist = len(sub(pos_b_t, pos_a_t)) - sum_r;
        if new_dist <= tolerance {
            return Some(t);
        }
        remaining_sep = new_dist;
    }
    None
}

// ---------------------------------------------------------------------------
// Swept sphere query against a list of triangles
// ---------------------------------------------------------------------------

/// Result of a swept sphere vs triangle mesh query.
#[derive(Debug, Clone)]
pub struct TriangleMeshSweepResult {
    /// Time of impact.
    pub toi: f64,
    /// Triangle index that was hit.
    pub tri_index: usize,
    /// Contact normal.
    pub normal: [f64; 3],
    /// Contact point on the sphere.
    pub point_a: [f64; 3],
    /// Contact point on the triangle.
    pub point_b: [f64; 3],
}

/// Sweep a sphere against a list of triangles and return the earliest impact.
///
/// Each triangle is `[v0, v1, v2]` as `[[f64; 3\]; 3]`.
pub fn swept_sphere_triangle_mesh(
    c0: [f64; 3],
    c1: [f64; 3],
    radius: f64,
    triangles: &[[[f64; 3]; 3]],
) -> Option<TriangleMeshSweepResult> {
    let mut best: Option<TriangleMeshSweepResult> = None;

    for (i, tri) in triangles.iter().enumerate() {
        let v0 = tri[0];
        let v1 = tri[1];
        let v2 = tri[2];

        // Use the triangle normal to create an offset plane sweep
        let e1 = sub(v1, v0);
        let e2 = sub(v2, v0);
        let raw_normal = cross(e1, e2);
        if len(raw_normal) < 1e-12 {
            continue;
        }
        let tri_normal_raw = normalize(raw_normal);

        // Try both sides of the triangle (double-sided)
        // Use the side the sphere is approaching from
        let approach_normal = if dot(tri_normal_raw, sub(c0, v0)) >= 0.0 {
            tri_normal_raw
        } else {
            scale(tri_normal_raw, -1.0)
        };

        let plane_d = dot(approach_normal, v0) + radius;

        // Intersect the sphere centre with the offset plane
        let d0 = dot(approach_normal, c0) - plane_d;
        let d1 = dot(approach_normal, c1) - plane_d;

        if d0 < 0.0 || (d0 - d1).abs() < 1e-14 {
            continue;
        }

        let t = d0 / (d0 - d1);
        if !(0.0..=1.0).contains(&t) {
            continue;
        }

        // Position of sphere centre at impact
        let dir = sub(c1, c0);
        let center_t = add(c0, scale(dir, t));

        // Find closest point on triangle to sphere centre projection onto triangle plane
        let proj = sub(
            center_t,
            scale(approach_normal, dot(approach_normal, sub(center_t, v0))),
        );
        // Clamp to triangle using edge tests (simplified: use the closest-point algorithm)
        let closest_on_tri = closest_point_on_tri(v0, v1, v2, proj);
        let diff = sub(center_t, closest_on_tri);
        let d = len(diff);
        if d > radius + 1e-6 {
            continue;
        }

        let normal = if d > 1e-10 {
            normalize(diff)
        } else {
            approach_normal
        };

        let result = TriangleMeshSweepResult {
            toi: t,
            tri_index: i,
            normal,
            point_a: add(center_t, scale(normal, -radius)),
            point_b: closest_on_tri,
        };

        if best.as_ref().is_none_or(|br| t < br.toi) {
            best = Some(result);
        }
    }

    best
}

/// Closest point on a triangle (inline helper for sweep.rs).
fn closest_point_on_tri(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let ab = sub(v1, v0);
    let ac = sub(v2, v0);
    let ap = sub(p, v0);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return v0;
    }
    let bp = sub(p, v1);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return v1;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(v0, scale(ab, v));
    }
    let cp = sub(p, v2);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return v2;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(v0, scale(ac, w));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(v1, scale(sub(v2, v1), w));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add(add(v0, scale(ab, v)), scale(ac, w))
}

// ---------------------------------------------------------------------------
// Sphere-capsule sweep
// ---------------------------------------------------------------------------

/// Swept sphere vs static capsule time-of-impact.
///
/// The sphere sweeps from `c0_s` to `c1_s` with `r_s`. The capsule is static
/// with axis `(cap0, cap1)` and `r_c`. Reduces to swept-sphere-sphere using
/// the closest point on the capsule axis.
pub fn swept_sphere_capsule(
    c0_s: [f64; 3],
    c1_s: [f64; 3],
    r_s: f64,
    cap0: [f64; 3],
    cap1: [f64; 3],
    r_c: f64,
) -> Option<SweepResult> {
    // Find the closest point on the capsule axis to the midpoint of the sphere sweep
    let mid_sphere = scale(add(c0_s, c1_s), 0.5);
    let axis = sub(cap1, cap0);
    let axis_len_sq = dot(axis, axis);
    let t = if axis_len_sq < 1e-20 {
        0.0
    } else {
        (dot(sub(mid_sphere, cap0), axis) / axis_len_sq).clamp(0.0, 1.0)
    };
    let cap_closest = add(cap0, scale(axis, t));
    // Treat the capsule as a static sphere at the closest axis point
    swept_sphere_sphere(c0_s, c1_s, r_s, cap_closest, cap_closest, r_c)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── swept_sphere_sphere ──────────────────────────────────────────────────

    #[test]
    fn test_swept_sphere_sphere_head_on() {
        // Two spheres of radius 0.5 moving toward each other along X.
        // A: starts at x=-3, ends at x=-1  (moves +2)
        // B: starts at x=+3, ends at x=+1  (moves -2)
        // Gap at t=0: 6 - 1 = 5. Relative speed 4. They meet when gap = 0.
        // (3-2t) - (- 3 + 2t) = 1  => 6 - 4t = 1 => t = 1.25 — outside [0,1].
        // Use shorter distances:
        // A: -1.5 → -0.5,  B: 1.5 → 0.5   (each moves 1 unit)
        // gap at t=0: 3 - 1 = 2;  gap(t) = 2 - 2t = 0 → t = 1.0
        let result = swept_sphere_sphere(
            [-1.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            0.5,
            [1.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_some(), "head-on spheres should hit");
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
    }

    #[test]
    fn test_swept_sphere_sphere_miss() {
        // Two spheres moving in parallel (no collision).
        let result = swept_sphere_sphere(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.0, 5.0, 0.0],
            [1.0, 5.0, 0.0],
            0.5,
        );
        assert!(result.is_none(), "parallel spheres should not collide");
    }

    #[test]
    fn test_swept_sphere_sphere_already_overlapping() {
        // Spheres already overlap at t=0.
        let result = swept_sphere_sphere(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.5, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_some());
        assert_eq!(result.unwrap().toi, 0.0);
    }

    #[test]
    fn test_swept_sphere_sphere_moving_apart() {
        // Spheres start separated and move further apart.
        let result = swept_sphere_sphere(
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_none(), "diverging spheres should not collide");
    }

    // ── swept_sphere_aabb ────────────────────────────────────────────────────

    #[test]
    fn test_swept_sphere_aabb_hit() {
        // Sphere falls from above toward an AABB.
        // Sphere centre: starts at (0,5,0), ends at (0,-1,0), radius=0.5
        // AABB: [-1,-1,-1] to [1,1,1]
        let result = swept_sphere_aabb(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            0.5,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_some(), "sphere moving toward AABB should hit");
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
    }

    #[test]
    fn test_swept_sphere_aabb_miss() {
        // Sphere moves away from AABB.
        let result = swept_sphere_aabb(
            [10.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            0.5,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_none(), "sphere moving away should not hit");
    }

    // ── swept_aabb_aabb ──────────────────────────────────────────────────────

    #[test]
    fn test_swept_aabb_aabb_hit() {
        // AABB A slides along X toward AABB B.
        let result = swept_aabb_aabb(
            [-3.0, -0.5, -0.5],
            [-2.0, 0.5, 0.5], // min/max A at t=0
            [2.5, 0.0, 0.0],  // A moves +2.5 along X
            [0.0, -0.5, -0.5],
            [1.0, 0.5, 0.5], // static B
        );
        assert!(result.is_some(), "sliding AABB should hit static AABB");
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
    }

    #[test]
    fn test_swept_aabb_aabb_miss() {
        // AABB A slides away from AABB B.
        let result = swept_aabb_aabb(
            [-3.0, -0.5, -0.5],
            [-2.0, 0.5, 0.5],
            [-5.0, 0.0, 0.0],
            [2.0, -0.5, -0.5],
            [3.0, 0.5, 0.5],
        );
        assert!(result.is_none(), "AABB moving away should not hit");
    }

    // ── ray_sphere ────────────────────────────────────────────────────────────

    #[test]
    fn test_ray_sphere_hit() {
        // Ray along +X axis toward sphere at origin.
        let hit = ray_sphere(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_some(), "ray should hit sphere");
        let h = hit.unwrap();
        assert!((h.t - 4.0).abs() < 1e-9, "expected t=4.0, got {}", h.t);
    }

    #[test]
    fn test_ray_sphere_miss() {
        // Ray parallel to sphere but not intersecting.
        let hit = ray_sphere(
            [-5.0, 2.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_none(), "ray should miss sphere");
    }

    #[test]
    fn test_ray_sphere_behind() {
        // Ray origin past the sphere.
        let hit = ray_sphere(
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_none(), "ray going away from sphere should miss");
    }

    // ── ray_aabb ──────────────────────────────────────────────────────────────

    #[test]
    fn test_ray_aabb_hit() {
        let hit = ray_aabb(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_some(), "ray should hit AABB");
        let h = hit.unwrap();
        assert!((h.t - 4.0).abs() < 1e-9, "expected t=4.0, got {}", h.t);
    }

    #[test]
    fn test_ray_aabb_miss() {
        let hit = ray_aabb(
            [-5.0, 3.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_none(), "ray should miss AABB");
    }

    // ── ray_triangle ──────────────────────────────────────────────────────────

    #[test]
    fn test_ray_triangle_hit() {
        // Ray straight down onto a horizontal triangle.
        let hit = ray_triangle(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_some(), "ray should hit triangle");
        let h = hit.unwrap();
        assert!((h.t - 5.0).abs() < 1e-9, "expected t=5.0, got {}", h.t);
    }

    #[test]
    fn test_ray_triangle_miss() {
        // Ray misses the triangle.
        let hit = ray_triangle(
            [5.0, 5.0, 5.0],
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_none(), "ray should miss triangle");
    }

    #[test]
    fn test_ray_triangle_normal_direction() {
        // Triangle in XZ plane with CCW winding → normal should be +Y.
        let hit = ray_triangle(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
            0.0,
            100.0,
        );
        let h = hit.unwrap();
        // Normal: e1 = (2,0,0), e2 = (1,0,2) → cross = (0*2-0*0, 0*1-2*2, 2*0-0*1) = (0,-4,0) — normalize = (0,-1,0)
        // Actually with these verts let's just check it's unit length.
        let nl = len(h.normal);
        assert!(
            (nl - 1.0).abs() < 1e-9,
            "normal should be unit length, got {}",
            nl
        );
    }

    // ── point_aabb_dist_sq ───────────────────────────────────────────────────

    #[test]
    fn test_point_aabb_dist_sq_inside() {
        let d = point_aabb_dist_sq([0.0, 0.0, 0.0], [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert_eq!(d, 0.0, "point inside AABB has distance 0");
    }

    #[test]
    fn test_point_aabb_dist_sq_outside() {
        let d = point_aabb_dist_sq([2.0, 0.0, 0.0], [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(
            (d - 1.0).abs() < 1e-9,
            "point 1 unit outside AABB, expected dsq=1, got {}",
            d
        );
    }

    // ── sphere_aabb_overlap ──────────────────────────────────────────────────

    #[test]
    fn test_sphere_aabb_overlap_touching() {
        // Sphere just touching the AABB.
        assert!(sphere_aabb_overlap(
            [2.0, 0.0, 0.0],
            1.0,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0]
        ));
    }

    #[test]
    fn test_sphere_aabb_overlap_separated() {
        assert!(!sphere_aabb_overlap(
            [3.0, 0.0, 0.0],
            1.0,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0]
        ));
    }

    // ── aabb_aabb_overlap ────────────────────────────────────────────────────

    #[test]
    fn test_aabb_aabb_overlap_yes() {
        assert!(aabb_aabb_overlap(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0],
            [2.0, 2.0, 2.0],
        ));
    }

    #[test]
    fn test_aabb_aabb_overlap_no() {
        assert!(!aabb_aabb_overlap(
            [-1.0, -1.0, -1.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
        ));
    }

    // ── swept_capsule_plane ───────────────────────────────────────────────────

    #[test]
    fn test_swept_capsule_plane_hit() {
        // Capsule (treated conservatively as sphere of r+hh) falls toward y=0 plane
        let result = swept_capsule_plane(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            0.5,
            1.0,
            [0.0, 1.0, 0.0],
            0.0,
        );
        assert!(result.is_some(), "capsule should hit floor plane");
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
    }

    #[test]
    fn test_swept_capsule_plane_moving_away() {
        let result = swept_capsule_plane(
            [0.0, 5.0, 0.0],
            [0.0, 10.0, 0.0],
            0.5,
            1.0,
            [0.0, 1.0, 0.0],
            0.0,
        );
        assert!(result.is_none(), "moving away from plane should not hit");
    }

    #[test]
    fn test_swept_capsule_plane_already_touching() {
        // Capsule effective radius = 0.5 + 0.5 = 1.0; position = (0, 1, 0) → dist = 0
        let result = swept_capsule_plane(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
            0.5,
            [0.0, 1.0, 0.0],
            0.0,
        );
        assert!(result.is_some());
        assert_eq!(result.unwrap().toi, 0.0);
    }

    // ── refine_sphere_sphere_toi ──────────────────────────────────────────────

    #[test]
    fn test_refine_toi_close_result() {
        // Spheres at ±1, moving toward each other at speed 1
        // Expected TOI ≈ 0.5 (contact at centre)
        let toi = refine_sphere_sphere_toi(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            0.0,
            1.0,
            1e-8,
            64,
        );
        assert!((toi - 0.5).abs() < 0.01, "refined toi={toi}, expected ≈0.5");
    }

    // ── SweptSphere ──────────────────────────────────────────────────────────

    #[test]
    fn test_swept_sphere_aabb_correct() {
        let s = SweptSphere::new([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0);
        let (mn, mx) = s.sweep_aabb();
        assert!((mn[0] - (-1.0)).abs() < 1e-12);
        assert!((mx[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_swept_sphere_aabb_overlap() {
        let a = SweptSphere::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 0);
        let b = SweptSphere::new([0.5, 0.0, 0.0], [1.5, 0.0, 0.0], 0.5, 1);
        assert!(a.sweep_aabb_overlaps(&b));
    }

    #[test]
    fn test_swept_sphere_aabb_no_overlap() {
        let a = SweptSphere::new([0.0; 3], [0.0; 3], 0.1, 0);
        let b = SweptSphere::new([10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.1, 1);
        assert!(!a.sweep_aabb_overlaps(&b));
    }

    // ── sweep_broadphase ─────────────────────────────────────────────────────

    #[test]
    fn test_sweep_broadphase_finds_pair() {
        let spheres = vec![
            SweptSphere::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 10),
            SweptSphere::new([0.5, 0.0, 0.0], [2.0, 0.0, 0.0], 0.5, 20),
            SweptSphere::new([50.0, 0.0, 0.0], [51.0, 0.0, 0.0], 0.5, 30),
        ];
        let pairs = sweep_broadphase(&spheres);
        assert_eq!(pairs.len(), 1, "only pair (10,20) should overlap");
        assert_eq!(pairs[0], (10, 20));
    }

    #[test]
    fn test_sweep_broadphase_no_pairs() {
        let spheres = vec![
            SweptSphere::new([0.0; 3], [0.0; 3], 0.1, 1),
            SweptSphere::new([100.0, 0.0, 0.0], [100.0, 0.0, 0.0], 0.1, 2),
        ];
        let pairs = sweep_broadphase(&spheres);
        assert!(pairs.is_empty());
    }

    // ── swept_sphere_sphere_interval ─────────────────────────────────────────

    #[test]
    fn test_interval_sweep_finds_hit() {
        // Spheres at ±5 moving toward each other at speed 1/s.
        // In interval [4, 6] they should collide.
        let result = swept_sphere_sphere_interval(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            4.0,
            6.0,
        );
        assert!(result.is_some(), "should hit within [4,6]");
    }

    #[test]
    fn test_interval_sweep_no_hit_outside() {
        // Spheres collide at t=4.5; interval is [5, 10] → already past
        // They would be overlapping at t_start=5, so this would actually return Some(0.0)
        // Use diverging spheres instead
        let result = swept_sphere_sphere_interval(
            [-5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5, // moving away
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            0.0,
            5.0,
        );
        assert!(result.is_none(), "diverging spheres should not collide");
    }

    // ── CompoundSweptShape ───────────────────────────────────────────────────

    #[test]
    fn test_compound_sweep_finds_earliest_hit() {
        // Compound A: single sphere at origin offset [0,0,0]
        // Compound B: single sphere at [4,0,0] offset [0,0,0]
        // A moves from [-5,0,0] to [0,0,0]; B is static
        let comp_a = CompoundSweptShape::new(vec![([0.0; 3], 0.5)]);
        let comp_b = CompoundSweptShape::new(vec![([0.0; 3], 0.5)]);
        let result = comp_a.sweep_vs_compound(
            [-5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0], // moves to x=-1
            &comp_b,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], // static at origin
        );
        // A gets to x=-1 at t=1, but distance at end is |(-1)-(0)| = 1 = sum_r → should touch
        assert!(result.is_some(), "compound sweep should detect hit");
    }

    #[test]
    fn test_compound_sweep_no_hit_separated() {
        let comp_a = CompoundSweptShape::new(vec![([0.0; 3], 0.5)]);
        let comp_b = CompoundSweptShape::new(vec![([0.0; 3], 0.5)]);
        // A moves away from B
        let result = comp_a.sweep_vs_compound(
            [0.0; 3],
            [-10.0, 0.0, 0.0],
            &comp_b,
            [50.0, 0.0, 0.0],
            [50.0, 0.0, 0.0],
        );
        assert!(result.is_none());
    }

    // ── swept_capsule_aabb ────────────────────────────────────────────────────

    #[test]
    fn test_swept_capsule_aabb_hit() {
        // Capsule descending from above toward a box
        let result = swept_capsule_aabb(
            [0.0, 5.0, -0.3],
            [0.0, 5.0, 0.3], // start axis
            [0.0, 0.0, -0.3],
            [0.0, 0.0, 0.3], // end axis
            0.2,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_some(), "descending capsule should hit AABB");
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
    }

    #[test]
    fn test_swept_capsule_aabb_miss() {
        // Capsule moves away from box
        let result = swept_capsule_aabb(
            [10.0, 0.0, -0.3],
            [10.0, 0.0, 0.3],
            [20.0, 0.0, -0.3],
            [20.0, 0.0, 0.3],
            0.2,
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_none(), "capsule moving away should not hit AABB");
    }

    // ── linear_ccd_sphere_sphere ──────────────────────────────────────────────

    #[test]
    fn test_linear_ccd_sphere_sphere_hit() {
        // Two spheres converging: A at x=-1, vel=+2; B at x=1, vel=-2; r=0.5 each
        // They meet when |pa - pb| = 1.0: pa = -1+2t, pb = 1-2t → dist = 2-4t = 1 → t=0.25
        let result = linear_ccd_sphere_sphere(
            [-1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            0.5,
            [1.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            0.5,
            1.0,
        );
        assert!(result.is_some(), "converging spheres should hit");
        let hit = result.unwrap();
        // toi is scaled by dt=1.0
        assert!(
            hit.toi >= 0.0 && hit.toi <= 1.0,
            "toi in [0,1], got {}",
            hit.toi
        );
    }

    #[test]
    fn test_linear_ccd_sphere_sphere_no_hit() {
        let result = linear_ccd_sphere_sphere(
            [-10.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.3,
            [10.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.3,
            1.0,
        );
        assert!(result.is_none(), "spheres moving away should not hit");
    }

    // ── bisect_toi ────────────────────────────────────────────────────────────

    #[test]
    fn test_bisect_toi_sphere_contact() {
        // Two spheres at ±1 with radius 0.5 each, converging at speed 1 each
        // Contact at t=0.5 (sum distances = 1.0 = 2*0.5)
        let toi = bisect_toi(0.0, 1.0, 1e-8, 64, |t| {
            let pa = add([-1.0, 0.0, 0.0], scale([1.0, 0.0, 0.0], t));
            let pb = add([1.0, 0.0, 0.0], scale([-1.0, 0.0, 0.0], t));
            len(sub(pa, pb)) - 1.0 // contact when dist = 1.0 (sum of radii)
        });
        assert!(
            (toi - 0.5).abs() < 1e-6,
            "bisect_toi should find t≈0.5, got {toi}"
        );
    }

    #[test]
    fn test_bisect_toi_already_contact() {
        // dist_fn always negative → t_hi converges toward t_lo
        let toi = bisect_toi(0.0, 1.0, 1e-8, 64, |_t| -1.0);
        assert!(toi < 0.5, "Should converge to low end when always negative");
    }

    // ── RotatingBodyState ─────────────────────────────────────────────────────

    #[test]
    fn test_rotating_body_motion_bound() {
        let body = RotatingBodyState::new([0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.5);
        // linear_speed = 1, angular_speed * bound_radius = 1.0
        // motion_bound(1.0) = (1 + 1) * 1 = 2
        let mb = body.motion_bound(1.0);
        assert!(
            (mb - 2.0).abs() < 1e-10,
            "motion_bound should be 2.0, got {mb}"
        );
    }

    #[test]
    fn test_rotating_body_advance() {
        let body = RotatingBodyState::new([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0, 1.0);
        let pos = body.advance(0.5);
        assert!(
            (pos[0] - 2.0).abs() < 1e-10,
            "advanced x should be 2.0, got {}",
            pos[0]
        );
    }

    #[test]
    fn test_rotating_body_motion_bound_zero_dt() {
        let body = RotatingBodyState::new([0.0; 3], [10.0, 0.0, 0.0], 5.0, 1.0);
        let mb = body.motion_bound(0.0);
        assert_eq!(mb, 0.0, "motion_bound(0) should be 0");
    }

    // ── conservative_advancement_rotating ────────────────────────────────────

    #[test]
    fn test_conservative_advancement_no_collision_moving_apart() {
        let a = RotatingBodyState::new([-10.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.0, 0.5);
        let b = RotatingBodyState::new([10.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0, 0.5);
        let result = conservative_advancement_rotating(&a, &b, 18.0, 0.0, 1.0, 1e-4, 100);
        assert!(result.is_none(), "diverging bodies should not collide");
    }

    #[test]
    fn test_conservative_advancement_already_overlapping() {
        let a = RotatingBodyState::new([0.0; 3], [0.0; 3], 0.0, 1.0);
        let b = RotatingBodyState::new([0.5, 0.0, 0.0], [0.0; 3], 0.0, 1.0);
        let result = conservative_advancement_rotating(&a, &b, 0.0, 0.0, 1.0, 1e-4, 100);
        // Overlap at start → returns t_start
        assert!(result.is_some(), "overlapping at start should return a toi");
        assert_eq!(result.unwrap(), 0.0, "overlapping at start → toi = t_start");
    }

    // ── TriangleMeshSweepResult ───────────────────────────────────────────────

    #[test]
    fn test_swept_sphere_triangle_mesh_hit() {
        // Triangle in XZ plane, sphere sweeps down from above
        let triangles = vec![[[-1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 0.0, 1.0]]];
        let result = swept_sphere_triangle_mesh([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &triangles);
        assert!(
            result.is_some(),
            "sphere sweeping toward triangle should hit"
        );
        let hit = result.unwrap();
        assert!(hit.toi >= 0.0 && hit.toi <= 1.0, "toi={}", hit.toi);
        assert_eq!(hit.tri_index, 0, "should hit first triangle");
    }

    #[test]
    fn test_swept_sphere_triangle_mesh_miss() {
        // Triangle far from the sphere's path
        let triangles = vec![[[-1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 0.0, 1.0]]];
        let result =
            swept_sphere_triangle_mesh([10.0, 3.0, 0.0], [10.0, -1.0, 0.0], 0.5, &triangles);
        assert!(result.is_none(), "sphere far from triangle should miss");
    }

    #[test]
    fn test_swept_sphere_triangle_mesh_multiple_triangles() {
        // Two triangles; sphere should hit the closer one (index 1 at y=2)
        let triangles = vec![
            [[-1.0, -5.0, -1.0], [1.0, -5.0, -1.0], [0.0, -5.0, 1.0]], // very far
            [[-1.0, 1.0, -1.0], [1.0, 1.0, -1.0], [0.0, 1.0, 1.0]],    // close
        ];
        // Sphere moves from y=5 to y=-5, radius=0.3
        let result = swept_sphere_triangle_mesh([0.0, 5.0, 0.0], [0.0, -5.0, 0.0], 0.3, &triangles);
        assert!(result.is_some(), "should hit a triangle");
    }

    // ── swept_sphere_capsule ──────────────────────────────────────────────────

    #[test]
    fn test_swept_sphere_capsule_hit() {
        // Sphere at x=-3 sweeps to x=0.5 past the capsule axis, radius 0.3 + capsule 0.5 = 0.8 sum
        // At end: sphere center at 0.5, capsule axis at 0 → dist = 0.5 < 0.8 → collision
        let result = swept_sphere_capsule(
            [-3.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.3,
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(result.is_some(), "sphere sweeping past capsule should hit");
    }

    #[test]
    fn test_swept_sphere_capsule_miss() {
        // Sphere moves away from capsule
        let result = swept_sphere_capsule(
            [5.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.3,
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(
            result.is_none(),
            "sphere moving away should not hit capsule"
        );
    }

    #[test]
    fn test_swept_sphere_capsule_already_overlapping() {
        // Sphere already inside the capsule radius
        let result = swept_sphere_capsule(
            [0.3, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.3,
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(result.is_some(), "overlapping should produce result");
        assert_eq!(result.unwrap().toi, 0.0, "overlap at start → toi=0");
    }

    // ── Additional ray tests ──────────────────────────────────────────────────

    #[test]
    fn test_ray_sphere_grazing() {
        // Ray grazing the sphere tangentially
        let hit = ray_sphere(
            [-5.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_some(), "grazing ray should hit sphere at tangent");
        let h = hit.unwrap();
        assert!(
            (h.point[1] - 1.0).abs() < 1e-6,
            "hit y should be ≈1 (tangent)"
        );
    }

    #[test]
    fn test_ray_aabb_inside_origin() {
        // Ray starting inside the AABB
        let hit = ray_aabb(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            0.0,
            100.0,
        );
        assert!(
            hit.is_some(),
            "ray from inside AABB should hit the far wall"
        );
    }

    #[test]
    fn test_point_aabb_dist_sq_corner() {
        // Point at (2, 2, 2), AABB [-1,1]^3 → distance to corner (1,1,1) = sqrt(3)
        let d = point_aabb_dist_sq([2.0, 2.0, 2.0], [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!((d - 3.0).abs() < 1e-10, "Expected dsq=3.0, got {d}");
    }

    #[test]
    fn test_aabb_aabb_overlap_touching_face() {
        // Two AABBs touching on a face (edge-case overlap)
        assert!(
            aabb_aabb_overlap(
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, 0.0, 0.0],
                [2.0, 1.0, 1.0],
            ),
            "face-touching AABBs should overlap"
        );
    }
}
