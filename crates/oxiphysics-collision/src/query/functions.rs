// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    AllPairsRayResult, NearestShape, PointQueryResult, QueryBox, QueryCapsule, QueryShapeRef,
    QuerySphere, RayCastResult, SegmentDistResult, SphereCastResult, ToiResult, ViewFrustum,
};

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / n)
    }
}
pub(super) fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if v > hi {
        hi
    } else {
        v
    }
}
/// Apply a 4×4 homogeneous transform to a point.
pub fn transform_point(p: [f64; 3], t: [[f64; 4]; 4]) -> [f64; 3] {
    [
        t[0][0] * p[0] + t[0][1] * p[1] + t[0][2] * p[2] + t[0][3],
        t[1][0] * p[0] + t[1][1] * p[1] + t[1][2] * p[2] + t[1][3],
        t[2][0] * p[0] + t[2][1] * p[1] + t[2][2] * p[2] + t[2][3],
    ]
}
/// Invert a rigid-body 4×4 transform (rotation + translation only).
pub(super) fn invert_transform(t: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let r = [
        [t[0][0], t[1][0], t[2][0]],
        [t[0][1], t[1][1], t[2][1]],
        [t[0][2], t[1][2], t[2][2]],
    ];
    let tr = [t[0][3], t[1][3], t[2][3]];
    let neg_rt = [
        -(r[0][0] * tr[0] + r[0][1] * tr[1] + r[0][2] * tr[2]),
        -(r[1][0] * tr[0] + r[1][1] * tr[1] + r[1][2] * tr[2]),
        -(r[2][0] * tr[0] + r[2][1] * tr[1] + r[2][2] * tr[2]),
    ];
    [
        [r[0][0], r[0][1], r[0][2], neg_rt[0]],
        [r[1][0], r[1][1], r[1][2], neg_rt[1]],
        [r[2][0], r[2][1], r[2][2], neg_rt[2]],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
/// Transform a ray by the *inverse* of `t` so the query can be done in local space.
pub fn transform_ray(origin: [f64; 3], dir: [f64; 3], t: [[f64; 4]; 4]) -> ([f64; 3], [f64; 3]) {
    let inv = invert_transform(t);
    let new_origin = transform_point(origin, inv);
    let new_dir = [
        inv[0][0] * dir[0] + inv[0][1] * dir[1] + inv[0][2] * dir[2],
        inv[1][0] * dir[0] + inv[1][1] * dir[1] + inv[1][2] * dir[2],
        inv[2][0] * dir[0] + inv[2][1] * dir[1] + inv[2][2] * dir[2],
    ];
    (new_origin, new_dir)
}
/// Common shape query interface.
pub trait ShapeQuery {
    /// Cast a ray from `ray_origin` in direction `ray_dir` (need not be normalised)
    /// up to parameter `max_toi`. Returns `None` when there is no hit.
    fn ray_cast(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_toi: f64,
    ) -> Option<RayCastResult>;
    /// Find the closest point on the shape to `point` and return it together
    /// with the signed distance and an inside flag.
    fn point_query(&self, point: [f64; 3]) -> PointQueryResult;
    /// Axis-aligned bounding box as (min, max) corner pair.
    fn aabb(&self) -> ([f64; 3], [f64; 3]);
}
/// Closest point on segment \[p0, p1\] to a given point.
/// Returns `(point, t)` where t ∈ \[0,1\].
pub fn closest_point_on_segment(p: [f64; 3], p0: [f64; 3], p1: [f64; 3]) -> ([f64; 3], f64) {
    let ab = sub(p1, p0);
    let ap = sub(p, p0);
    let len_sq = dot(ab, ab);
    if len_sq < 1e-30 {
        return (p0, 0.0);
    }
    let t = clamp(dot(ap, ab) / len_sq, 0.0, 1.0);
    let point = add(p0, scale(ab, t));
    (point, t)
}
/// Signed distance from point `p` to the segment \[p0, p1\] (unsigned, always ≥ 0).
pub fn dist_point_to_segment(p: [f64; 3], p0: [f64; 3], p1: [f64; 3]) -> f64 {
    let (cp, _) = closest_point_on_segment(p, p0, p1);
    norm(sub(p, cp))
}
/// Closest point on triangle (a, b, c) to point `p`.
/// Returns the closest point as a `[f64; 3]`.
pub fn closest_point_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let cp_pt = sub(p, c);
    let d5 = dot(ab, cp_pt);
    let d6 = dot(ac, cp_pt);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(a, scale(ab, v));
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(a, scale(ac, w));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add(a, add(scale(ab, v), scale(ac, w)))
}
/// Signed distance from point `p` to the triangle (a, b, c) surface.
/// Always returns a non-negative distance (unsigned).
pub fn dist_point_to_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let cp = closest_point_on_triangle(p, a, b, c);
    norm(sub(p, cp))
}
/// Test whether point `p` is inside a convex hull defined by a set of half-space
/// planes.  Each plane is `[a, b, c, d]` where `ax + by + cz + d >= 0` is "inside".
///
/// Returns `true` if `p` satisfies all half-spaces.
pub fn point_in_convex_hull(p: [f64; 3], planes: &[[f64; 4]]) -> bool {
    for plane in planes {
        let n = [plane[0], plane[1], plane[2]];
        let d = plane[3];
        if dot(p, n) + d < -1e-10 {
            return false;
        }
    }
    true
}
/// Test whether point `p` is inside a convex polygon mesh (list of triangles).
/// Computes the half-spaces from each triangle and checks all.
///
/// `triangles` is a list of `(a, b, c)` triangles; the mesh must be convex
/// and outward-facing normals.
pub fn point_in_convex_mesh(p: [f64; 3], triangles: &[([f64; 3], [f64; 3], [f64; 3])]) -> bool {
    for &(a, b, c) in triangles {
        let ab = sub(b, a);
        let ac = sub(c, a);
        let n = normalize(cross3(ab, ac));
        let d = dot(n, a);
        if dot(n, p) - d > 1e-8 {
            return false;
        }
    }
    true
}
/// Signed distance from point `p` to a convex sphere (centered at `center`, radius `r`).
/// Negative if inside.
pub fn signed_dist_to_sphere(p: [f64; 3], center: [f64; 3], radius: f64) -> f64 {
    norm(sub(p, center)) - radius
}
/// Signed distance from point `p` to an axis-aligned box centered at `box_center`
/// with half-extents `he`.  Negative if inside.
pub fn signed_dist_to_box(p: [f64; 3], box_center: [f64; 3], he: [f64; 3]) -> f64 {
    let q = [
        (p[0] - box_center[0]).abs() - he[0],
        (p[1] - box_center[1]).abs() - he[1],
        (p[2] - box_center[2]).abs() - he[2],
    ];
    let outside_dist = norm([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
    let inside_dist = q[0].max(q[1]).max(q[2]).min(0.0);
    outside_dist + inside_dist
}
/// Signed distance from point `p` to a capsule with segment \[p0, p1\] and radius r.
pub fn signed_dist_to_capsule(p: [f64; 3], p0: [f64; 3], p1: [f64; 3], r: f64) -> f64 {
    let (cp, _) = closest_point_on_segment(p, p0, p1);
    norm(sub(p, cp)) - r
}
/// Ray vs sphere intersection returning the first hit time.
/// Returns `Some(t)` if hit, `None` otherwise.
pub fn ray_vs_sphere(
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    center: [f64; 3],
    radius: f64,
    max_toi: f64,
) -> Option<f64> {
    let sphere = QuerySphere { center, radius };
    sphere.ray_cast(ray_origin, ray_dir, max_toi).map(|r| r.toi)
}
/// Ray vs axis-aligned box intersection returning the first hit time.
pub fn ray_vs_box(
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    box_center: [f64; 3],
    half_extents: [f64; 3],
    max_toi: f64,
) -> Option<f64> {
    let b = QueryBox {
        center: box_center,
        half_extents,
    };
    b.ray_cast(ray_origin, ray_dir, max_toi).map(|r| r.toi)
}
/// Ray vs triangle intersection (Möller-Trumbore algorithm).
/// Returns `Some(t)` if the ray hits the triangle front face.
pub fn ray_vs_triangle(
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    max_toi: f64,
) -> Option<f64> {
    let edge1 = sub(v1, v0);
    let edge2 = sub(v2, v0);
    let h = cross3(ray_dir, edge2);
    let a = dot(edge1, h);
    if a.abs() < 1e-10 {
        return None;
    }
    let f = 1.0 / a;
    let s = sub(ray_origin, v0);
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross3(s, edge1);
    let v = f * dot(ray_dir, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot(edge2, q);
    if t >= 0.0 && t <= max_toi {
        Some(t)
    } else {
        None
    }
}
/// Moving sphere sweep vs static sphere.
///
/// Returns `Some(t)` in `[0, max_t]` where the sphere first touches the static sphere.
pub fn sphere_sweep_vs_sphere(
    start: [f64; 3],
    velocity: [f64; 3],
    moving_radius: f64,
    static_center: [f64; 3],
    static_radius: f64,
    max_t: f64,
) -> Option<f64> {
    let combined_r = moving_radius + static_radius;
    let dx = start[0] - static_center[0];
    let dy = start[1] - static_center[1];
    let dz = start[2] - static_center[2];
    let a = dot(velocity, velocity);
    let b = 2.0 * (dx * velocity[0] + dy * velocity[1] + dz * velocity[2]);
    let c = dx * dx + dy * dy + dz * dz - combined_r * combined_r;
    if a < 1e-30 {
        return if c <= 0.0 { Some(0.0) } else { None };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    if t1 >= 0.0 && t1 <= max_t {
        Some(t1)
    } else if t2 >= 0.0 && t2 <= max_t {
        Some(t2)
    } else {
        None
    }
}
/// Moving sphere sweep vs static AABB.
///
/// Uses the Minkowski sum: expand box by sphere radius and cast a ray.
pub fn sphere_sweep_vs_box(
    start: [f64; 3],
    velocity: [f64; 3],
    sphere_radius: f64,
    box_center: [f64; 3],
    box_half_extents: [f64; 3],
    max_t: f64,
) -> Option<f64> {
    let expanded_he = [
        box_half_extents[0] + sphere_radius,
        box_half_extents[1] + sphere_radius,
        box_half_extents[2] + sphere_radius,
    ];
    let b = QueryBox {
        center: box_center,
        half_extents: expanded_he,
    };
    b.ray_cast(start, velocity, max_t).map(|r| r.toi)
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Test whether two axis-aligned boxes overlap (no contact data).
///
/// Returns `true` if the AABB `[min_a, max_a]` overlaps `[min_b, max_b]`.
pub fn aabb_overlap(min_a: [f64; 3], max_a: [f64; 3], min_b: [f64; 3], max_b: [f64; 3]) -> bool {
    min_a[0] <= max_b[0]
        && max_a[0] >= min_b[0]
        && min_a[1] <= max_b[1]
        && max_a[1] >= min_b[1]
        && min_a[2] <= max_b[2]
        && max_a[2] >= min_b[2]
}
/// Test whether two spheres overlap (no contact data).
pub fn sphere_sphere_overlap(
    center_a: [f64; 3],
    radius_a: f64,
    center_b: [f64; 3],
    radius_b: f64,
) -> bool {
    let d = sub(center_a, center_b);
    let dist_sq = dot(d, d);
    let r_sum = radius_a + radius_b;
    dist_sq <= r_sum * r_sum
}
/// Test whether a sphere overlaps an axis-aligned box.
///
/// Uses the standard closest-point projection.
pub fn sphere_aabb_overlap(
    center: [f64; 3],
    radius: f64,
    box_center: [f64; 3],
    box_half_extents: [f64; 3],
) -> bool {
    let lo = sub(box_center, box_half_extents);
    let hi = add(box_center, box_half_extents);
    let closest = [
        clamp(center[0], lo[0], hi[0]),
        clamp(center[1], lo[1], hi[1]),
        clamp(center[2], lo[2], hi[2]),
    ];
    let d = sub(center, closest);
    dot(d, d) <= radius * radius
}
/// Test whether two capsules overlap (no contact data).
///
/// Computes the closest distance between the two central segments and compares
/// it against the sum of radii.
pub fn capsule_capsule_overlap(
    p0_a: [f64; 3],
    p1_a: [f64; 3],
    radius_a: f64,
    p0_b: [f64; 3],
    p1_b: [f64; 3],
    radius_b: f64,
) -> bool {
    let da = sub(p1_a, p0_a);
    let db = sub(p1_b, p0_b);
    let r = sub(p0_a, p0_b);
    let a = dot(da, da);
    let e = dot(db, db);
    let f = dot(db, r);
    let (s, t) = if a <= 1e-15 && e <= 1e-15 {
        (0.0, 0.0)
    } else if a <= 1e-15 {
        (0.0, clamp(f / e, 0.0, 1.0))
    } else {
        let c = dot(da, r);
        if e <= 1e-15 {
            (clamp(-c / a, 0.0, 1.0), 0.0)
        } else {
            let b = dot(da, db);
            let denom = a * e - b * b;
            let s = if denom.abs() > 1e-15 {
                clamp((b * f - c * e) / denom, 0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            let (t_clamped, s_final) = if t < 0.0 {
                (0.0, clamp(-c / a, 0.0, 1.0))
            } else if t > 1.0 {
                (1.0, clamp((b - c) / a, 0.0, 1.0))
            } else {
                (t, s)
            };
            (s_final, t_clamped)
        }
    };
    let closest_a = add(p0_a, scale(da, s));
    let closest_b = add(p0_b, scale(db, t));
    let d = sub(closest_a, closest_b);
    let r_sum = radius_a + radius_b;
    dot(d, d) <= r_sum * r_sum
}
/// Cast a ray against a list of shapes (brute-force all-pairs).
///
/// Returns a sorted list of hits (nearest first).
pub fn ray_cast_all_pairs(
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    max_toi: f64,
    shapes: &[QueryShapeRef<'_>],
) -> Vec<AllPairsRayResult> {
    let mut results: Vec<AllPairsRayResult> = Vec::new();
    for (idx, shape) in shapes.iter().enumerate() {
        let hit = match shape {
            QueryShapeRef::Sphere(s) => s.ray_cast(ray_origin, ray_dir, max_toi),
            QueryShapeRef::Box(b) => b.ray_cast(ray_origin, ray_dir, max_toi),
            QueryShapeRef::Capsule(c) => c.ray_cast(ray_origin, ray_dir, max_toi),
        };
        if let Some(r) = hit {
            results.push(AllPairsRayResult {
                shape_index: idx,
                toi: r.toi,
                normal: r.normal,
                point: r.point,
            });
        }
    }
    results.sort_by(|a, b| {
        a.toi
            .partial_cmp(&b.toi)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results
}
/// Find all shapes that contain the given point (brute-force all-pairs).
///
/// Returns a list of shape indices whose interiors contain `point`.
pub fn point_containment_all_pairs(point: [f64; 3], shapes: &[QueryShapeRef<'_>]) -> Vec<usize> {
    let mut inside = Vec::new();
    for (idx, shape) in shapes.iter().enumerate() {
        let r = match shape {
            QueryShapeRef::Sphere(s) => s.point_query(point),
            QueryShapeRef::Box(b) => b.point_query(point),
            QueryShapeRef::Capsule(c) => c.point_query(point),
        };
        if r.inside {
            inside.push(idx);
        }
    }
    inside
}
/// Query shapes potentially visible from the given frustum (brute-force).
///
/// Returns a list of shape indices that are not fully culled by the frustum.
pub fn frustum_query(frustum: &ViewFrustum, shapes: &[QueryShapeRef<'_>]) -> Vec<usize> {
    let mut visible = Vec::new();
    for (idx, shape) in shapes.iter().enumerate() {
        let visible_flag = match shape {
            QueryShapeRef::Sphere(s) => frustum.test_sphere(s.center, s.radius),
            QueryShapeRef::Box(b) => {
                let (lo, hi) = b.aabb();
                frustum.test_aabb(lo, hi)
            }
            QueryShapeRef::Capsule(c) => {
                let (lo, hi) = c.aabb();
                frustum.test_aabb(lo, hi)
            }
        };
        if visible_flag {
            visible.push(idx);
        }
    }
    visible
}
/// Find the K nearest shapes to `query_point`, sorted by ascending distance.
///
/// Shapes that contain the point have a negative distance and appear first.
pub fn k_nearest_shapes(
    query_point: [f64; 3],
    shapes: &[QueryShapeRef<'_>],
    k: usize,
) -> Vec<NearestShape> {
    let mut all: Vec<NearestShape> = shapes
        .iter()
        .enumerate()
        .map(|(idx, shape)| {
            let r = match shape {
                QueryShapeRef::Sphere(s) => s.point_query(query_point),
                QueryShapeRef::Box(b) => b.point_query(query_point),
                QueryShapeRef::Capsule(c) => c.point_query(query_point),
            };
            NearestShape {
                index: idx,
                distance: r.distance,
                closest_point: r.point,
            }
        })
        .collect();
    all.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    all.truncate(k);
    all
}
/// Cast a moving sphere of radius `sphere_radius` from `start` along `velocity`
/// against a scene of static shapes.  Returns the first hit within `max_t`.
///
/// Implemented by expanding each shape by `sphere_radius` (Minkowski sum) and
/// casting a ray from `start` in the velocity direction.
pub fn sphere_cast_scene(
    start: [f64; 3],
    velocity: [f64; 3],
    sphere_radius: f64,
    max_t: f64,
    shapes: &[QueryShapeRef<'_>],
) -> Option<SphereCastResult> {
    let mut best_toi = max_t + 1.0;
    let mut best_result: Option<SphereCastResult> = None;
    for (idx, shape) in shapes.iter().enumerate() {
        match shape {
            QueryShapeRef::Sphere(s) => {
                let expanded = QuerySphere {
                    center: s.center,
                    radius: s.radius + sphere_radius,
                };
                if let Some(r) = expanded.ray_cast(start, velocity, max_t)
                    && r.toi < best_toi
                {
                    best_toi = r.toi;
                    best_result = Some(SphereCastResult {
                        shape_index: idx,
                        toi: r.toi,
                        normal: r.normal,
                    });
                }
            }
            QueryShapeRef::Box(b) => {
                let expanded = QueryBox {
                    center: b.center,
                    half_extents: [
                        b.half_extents[0] + sphere_radius,
                        b.half_extents[1] + sphere_radius,
                        b.half_extents[2] + sphere_radius,
                    ],
                };
                if let Some(r) = expanded.ray_cast(start, velocity, max_t)
                    && r.toi < best_toi
                {
                    best_toi = r.toi;
                    best_result = Some(SphereCastResult {
                        shape_index: idx,
                        toi: r.toi,
                        normal: r.normal,
                    });
                }
            }
            QueryShapeRef::Capsule(c) => {
                let expanded = QueryCapsule {
                    p0: c.p0,
                    p1: c.p1,
                    radius: c.radius + sphere_radius,
                };
                if let Some(r) = expanded.ray_cast(start, velocity, max_t)
                    && r.toi < best_toi
                {
                    best_toi = r.toi;
                    best_result = Some(SphereCastResult {
                        shape_index: idx,
                        toi: r.toi,
                        normal: r.normal,
                    });
                }
            }
        }
    }
    best_result
}
/// Compute the time of impact for two linearly-moving spheres.
///
/// Both spheres move at constant velocity over the interval `[0, max_t]`.
/// Returns `None` if no collision occurs within `max_t`.
pub fn compute_time_of_impact_spheres(
    center_a: [f64; 3],
    vel_a: [f64; 3],
    radius_a: f64,
    center_b: [f64; 3],
    vel_b: [f64; 3],
    radius_b: f64,
    max_t: f64,
) -> Option<ToiResult> {
    let rel_vel = sub(vel_a, vel_b);
    let d = sub(center_a, center_b);
    let combined_r = radius_a + radius_b;
    let a = dot(rel_vel, rel_vel);
    let b = 2.0 * dot(d, rel_vel);
    let c = dot(d, d) - combined_r * combined_r;
    if a < 1e-30 {
        return if c <= 0.0 {
            let n = if norm(d) > 1e-15 {
                normalize(d)
            } else {
                [1.0, 0.0, 0.0]
            };
            let contact = add(center_b, scale(n, radius_b));
            Some(ToiResult {
                toi: 0.0,
                normal: n,
                point: contact,
            })
        } else {
            None
        };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    let toi = if t1 >= 0.0 && t1 <= max_t {
        t1
    } else if t2 >= 0.0 && t2 <= max_t {
        t2
    } else {
        return None;
    };
    let pos_a = add(center_a, scale(vel_a, toi));
    let pos_b = add(center_b, scale(vel_b, toi));
    let diff = sub(pos_a, pos_b);
    let n = if norm(diff) > 1e-15 {
        normalize(diff)
    } else {
        [1.0, 0.0, 0.0]
    };
    let point = add(pos_b, scale(n, radius_b));
    Some(ToiResult {
        toi,
        normal: n,
        point,
    })
}
/// Compute the time of impact for two linearly-moving AABBs.
///
/// Uses the separating axis theorem on all three axes.  Returns `None` if the
/// AABBs do not collide within `[0, max_t]`.
pub fn compute_time_of_impact_aabbs(
    min_a: [f64; 3],
    max_a: [f64; 3],
    vel_a: [f64; 3],
    min_b: [f64; 3],
    max_b: [f64; 3],
    vel_b: [f64; 3],
    max_t: f64,
) -> Option<f64> {
    let rel_vel = sub(vel_a, vel_b);
    let mut t_enter = 0.0_f64;
    let mut t_exit = max_t;
    for i in 0..3 {
        let v = rel_vel[i];
        let d0 = min_b[i] - max_a[i];
        let d1 = max_b[i] - min_a[i];
        if v.abs() < 1e-15 {
            if min_a[i] > max_b[i] || min_b[i] > max_a[i] {
                return None;
            }
        } else {
            let (ta, tb) = if v > 0.0 {
                (d0 / v, d1 / v)
            } else {
                (d1 / v, d0 / v)
            };
            t_enter = t_enter.max(ta);
            t_exit = t_exit.min(tb);
            if t_enter > t_exit {
                return None;
            }
        }
    }
    if t_enter < 0.0 || t_enter > max_t {
        return None;
    }
    Some(t_enter)
}
/// Project a point onto the surface of a convex shape represented by a set of
/// support vertices (convex hull vertices).
///
/// Returns the closest point on the convex hull surface to `point`.
///
/// **Algorithm**: finds the closest vertex, then checks each edge and face
/// (triangle fan from centroid) to refine.  For a full convex hull this is a
/// simplified GJK-like approach; here we use the closest-point-to-triangle
/// test on the convex polygon formed by adjacent vertices.
pub fn project_point_convex(point: [f64; 3], vertices: &[[f64; 3]]) -> [f64; 3] {
    if vertices.is_empty() {
        return point;
    }
    if vertices.len() == 1 {
        return vertices[0];
    }
    if vertices.len() == 2 {
        let (cp, _) = closest_point_on_segment(point, vertices[0], vertices[1]);
        return cp;
    }
    let n = vertices.len() as f64;
    let mut centroid = [0.0_f64; 3];
    for v in vertices {
        centroid[0] += v[0];
        centroid[1] += v[1];
        centroid[2] += v[2];
    }
    centroid[0] /= n;
    centroid[1] /= n;
    centroid[2] /= n;
    let mut best_dist_sq = f64::INFINITY;
    let mut best_pt = vertices[0];
    let nv = vertices.len();
    for i in 0..nv {
        let a = vertices[i];
        let b = vertices[(i + 1) % nv];
        let cp = closest_point_on_triangle_local(point, centroid, a, b);
        let d = sub(point, cp);
        let dsq = dot(d, d);
        if dsq < best_dist_sq {
            best_dist_sq = dsq;
            best_pt = cp;
        }
    }
    best_pt
}
/// Closest point on a triangle `(v0, v1, v2)` to `point` (internal helper).
pub(super) fn closest_point_on_triangle_local(
    point: [f64; 3],
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
) -> [f64; 3] {
    let ab = sub(v1, v0);
    let ac = sub(v2, v0);
    let ap = sub(point, v0);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return v0;
    }
    let bp = sub(point, v1);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return v1;
    }
    let cp2 = sub(point, v2);
    let d5 = dot(ab, cp2);
    let d6 = dot(ac, cp2);
    if d6 >= 0.0 && d5 <= d6 {
        return v2;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(v0, scale(ab, v));
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
/// Compute the minimum distance (and closest points) between two line segments.
///
/// Segments are defined by their endpoints: `(p0_a, p1_a)` and `(p0_b, p1_b)`.
pub fn distance_between_segments(
    p0_a: [f64; 3],
    p1_a: [f64; 3],
    p0_b: [f64; 3],
    p1_b: [f64; 3],
) -> SegmentDistResult {
    let da = sub(p1_a, p0_a);
    let db = sub(p1_b, p0_b);
    let r = sub(p0_a, p0_b);
    let a = dot(da, da);
    let e = dot(db, db);
    let f = dot(db, r);
    let (s, t) = if a <= 1e-15 {
        let t = clamp(f / e.max(1e-15), 0.0, 1.0);
        (0.0, t)
    } else {
        let c = dot(da, r);
        if e <= 1e-15 {
            let s = clamp(-c / a, 0.0, 1.0);
            (s, 0.0)
        } else {
            let b = dot(da, db);
            let denom = a * e - b * b;
            let s = if denom.abs() > 1e-15 {
                clamp((b * f - c * e) / denom, 0.0, 1.0)
            } else {
                0.0
            };
            let t_raw = (b * s + f) / e;
            let t = clamp(t_raw, 0.0, 1.0);
            let s = clamp((-c + b * t) / a, 0.0, 1.0);
            (s, t)
        }
    };
    let point_a = add(p0_a, scale(da, s));
    let point_b = add(p0_b, scale(db, t));
    let diff = sub(point_a, point_b);
    let distance = norm(diff);
    SegmentDistResult {
        point_a,
        point_b,
        distance,
    }
}
/// Closest point on a sphere surface to `point`.
pub fn closest_point_on_sphere(point: [f64; 3], center: [f64; 3], radius: f64) -> [f64; 3] {
    let s = QuerySphere { center, radius };
    s.point_query(point).point
}
/// Closest point on an AABB surface to `point`.
pub fn closest_point_on_box(
    point: [f64; 3],
    box_center: [f64; 3],
    half_extents: [f64; 3],
) -> [f64; 3] {
    let b = QueryBox {
        center: box_center,
        half_extents,
    };
    b.point_query(point).point
}
/// Closest point on a capsule surface to `point`.
pub fn closest_point_on_capsule(
    point: [f64; 3],
    p0: [f64; 3],
    p1: [f64; 3],
    radius: f64,
) -> [f64; 3] {
    let c = QueryCapsule { p0, p1, radius };
    c.point_query(point).point
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::CompoundQuery;
    use crate::query::QueryBox;

    use crate::query::QuerySphere;

    #[test]
    fn sphere_ray_cast_hits() {
        let sphere = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let result = sphere.ray_cast([3.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some(), "expected a hit");
        let r = result.unwrap();
        assert!((r.toi - 2.0).abs() < 1e-10, "toi={}", r.toi);
    }
    #[test]
    fn sphere_ray_cast_miss() {
        let sphere = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let result = sphere.ray_cast([3.0, 2.0, 0.0], [-1.0, 0.0, 0.0], 100.0);
        assert!(result.is_none(), "expected a miss");
    }
    #[test]
    fn sphere_point_query_outside() {
        let sphere = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let r = sphere.point_query([2.0, 0.0, 0.0]);
        assert!(!r.inside);
        assert!((r.distance - 1.0).abs() < 1e-10, "distance={}", r.distance);
        assert!((r.point[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn box_ray_cast_hits_face() {
        let b = QueryBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let result = b.ray_cast([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some(), "expected a hit");
        let r = result.unwrap();
        assert!((r.toi - 4.0).abs() < 1e-10, "toi={}", r.toi);
        assert!(
            (r.normal[0] - (-1.0)).abs() < 1e-10,
            "normal={:?}",
            r.normal
        );
    }
    #[test]
    fn box_point_query_inside() {
        let b = QueryBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [2.0, 2.0, 2.0],
        };
        let r = b.point_query([0.0, 0.0, 0.0]);
        assert!(r.inside, "expected inside");
        assert!(r.distance < 0.0, "signed distance should be negative");
    }
    #[test]
    fn box_point_query_outside() {
        let b = QueryBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let r = b.point_query([3.0, 0.0, 0.0]);
        assert!(!r.inside);
        assert!((r.distance - 2.0).abs() < 1e-10, "distance={}", r.distance);
    }
    #[test]
    fn capsule_point_query_between_endpoints() {
        let cap = QueryCapsule {
            p0: [0.0, 0.0, -2.0],
            p1: [0.0, 0.0, 2.0],
            radius: 1.0,
        };
        let r = cap.point_query([2.0, 0.0, 0.0]);
        assert!(!r.inside);
        assert!((r.distance - 1.0).abs() < 1e-10, "distance={}", r.distance);
        assert!((r.point[0] - 1.0).abs() < 1e-10, "closest_x={}", r.point[0]);
        assert!(r.point[2].abs() < 1e-10, "closest_z={}", r.point[2]);
    }
    #[test]
    fn capsule_ray_cast_hits() {
        let cap = QueryCapsule {
            p0: [0.0, -1.0, 0.0],
            p1: [0.0, 1.0, 0.0],
            radius: 1.0,
        };
        let result = cap.ray_cast([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some(), "expected a hit");
        let r = result.unwrap();
        assert!((r.toi - 4.0).abs() < 1e-10, "toi={}", r.toi);
    }
    #[test]
    fn transform_point_identity() {
        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let p = [1.0, 2.0, 3.0];
        let q = transform_point(p, identity);
        assert!((q[0] - 1.0).abs() < 1e-10);
        assert!((q[1] - 2.0).abs() < 1e-10);
        assert!((q[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn transform_point_translation() {
        let t = [
            [1.0, 0.0, 0.0, 3.0],
            [0.0, 1.0, 0.0, 5.0],
            [0.0, 0.0, 1.0, 7.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let p = [1.0, 0.0, 0.0];
        let q = transform_point(p, t);
        assert!((q[0] - 4.0).abs() < 1e-10);
        assert!((q[1] - 5.0).abs() < 1e-10);
        assert!((q[2] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn closest_on_segment_midpoint() {
        let (cp, t) = closest_point_on_segment([1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((t - 0.5).abs() < 1e-10, "t={t}");
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
    }
    #[test]
    fn closest_on_segment_clamps_start() {
        let (cp, t) = closest_point_on_segment([-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((t).abs() < 1e-10, "t should be 0, got {t}");
        assert_eq!(cp, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn closest_on_segment_clamps_end() {
        let (cp, t) = closest_point_on_segment([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((t - 1.0).abs() < 1e-10, "t should be 1, got {t}");
        assert!((cp[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn closest_on_triangle_interior() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 0.0, 0.0];
        let c = [1.5, 2.0, 0.0];
        let p = [1.5, 0.6, 2.0];
        let cp = closest_point_on_triangle(p, a, b, c);
        assert!(cp[2].abs() < 1e-8, "z={}", cp[2]);
    }
    #[test]
    fn closest_on_triangle_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [10.0, 0.0, 0.0];
        let c = [0.0, 10.0, 0.0];
        let p = [-2.0, -2.0, 0.0];
        let cp = closest_point_on_triangle(p, a, b, c);
        assert!((cp[0]).abs() < 1e-10);
        assert!((cp[1]).abs() < 1e-10);
    }
    #[test]
    fn point_in_unit_box_hull() {
        let planes: Vec<[f64; 4]> = vec![
            [1.0, 0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0, 1.0],
        ];
        assert!(point_in_convex_hull([0.5, 0.5, 0.5], &planes));
        assert!(!point_in_convex_hull([2.0, 0.5, 0.5], &planes));
        assert!(!point_in_convex_hull([-0.1, 0.5, 0.5], &planes));
    }
    #[test]
    fn signed_dist_sphere_outside() {
        let d = signed_dist_to_sphere([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        assert!((d - 2.0).abs() < 1e-10, "d={d}");
    }
    #[test]
    fn signed_dist_sphere_inside() {
        let d = signed_dist_to_sphere([0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        assert!(d < 0.0, "should be negative inside, d={d}");
        assert!((d - (-0.5)).abs() < 1e-10, "d={d}");
    }
    #[test]
    fn signed_dist_box_outside() {
        let d = signed_dist_to_box([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((d - 2.0).abs() < 1e-10, "d={d}");
    }
    #[test]
    fn signed_dist_box_inside() {
        let d = signed_dist_to_box([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(d < 0.0, "center should be inside, d={d}");
    }
    #[test]
    fn signed_dist_capsule_outside() {
        let d = signed_dist_to_capsule([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!((d - 1.5).abs() < 1e-10, "d={d}");
    }
    #[test]
    fn signed_dist_capsule_inside() {
        let d = signed_dist_to_capsule([0.1, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!(d < 0.0, "should be inside capsule, d={d}");
    }
    #[test]
    fn ray_vs_sphere_hit() {
        let t = ray_vs_sphere(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            100.0,
        );
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn ray_vs_sphere_miss() {
        let t = ray_vs_sphere(
            [-5.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            100.0,
        );
        assert!(t.is_none());
    }
    #[test]
    fn ray_vs_box_hit() {
        let t = ray_vs_box(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            100.0,
        );
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn ray_vs_triangle_hit() {
        let t = ray_vs_triangle(
            [0.5, 0.5, 5.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            100.0,
        );
        assert!(t.is_some(), "ray straight down should hit triangle");
        assert!((t.unwrap() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn ray_vs_triangle_miss_parallel() {
        let t = ray_vs_triangle(
            [0.5, 0.5, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            100.0,
        );
        assert!(t.is_none());
    }
    #[test]
    fn sphere_sweep_vs_sphere_hit() {
        let t = sphere_sweep_vs_sphere(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            0.5,
            10.0,
        );
        assert!(t.is_some(), "sweep should hit");
        assert!((t.unwrap() - 4.0).abs() < 1e-10, "t={}", t.unwrap());
    }
    #[test]
    fn sphere_sweep_vs_sphere_miss() {
        let t = sphere_sweep_vs_sphere(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.0, 5.0, 0.0],
            0.5,
            10.0,
        );
        assert!(t.is_none());
    }
    #[test]
    fn sphere_sweep_vs_box_hit() {
        let t = sphere_sweep_vs_box(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            100.0,
        );
        assert!(t.is_some(), "sweep should hit expanded box");
    }
    #[test]
    fn compound_ray_cast_finds_closest() {
        let mut compound = CompoundQuery::new();
        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        compound.add_shape(
            Box::new(QuerySphere {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
            }),
            identity,
        );
        compound.add_shape(
            Box::new(QuerySphere {
                center: [10.0, 0.0, 0.0],
                radius: 1.0,
            }),
            identity,
        );
        let result = compound.ray_cast_any([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!((r.toi - 4.0).abs() < 1e-10, "toi={}", r.toi);
    }
}
#[cfg(test)]
mod tests_extended {
    use super::*;

    use crate::query::QueryBox;
    use crate::query::QueryShapeRef;
    use crate::query::QuerySphere;
    use crate::query::ViewFrustum;
    use crate::query::aabb_overlap;
    use crate::query::capsule_capsule_overlap;
    use crate::query::closest_point_on_capsule;
    use crate::query::compute_time_of_impact_aabbs;
    use crate::query::compute_time_of_impact_spheres;
    use crate::query::distance_between_segments;
    use crate::query::frustum_query;
    use crate::query::k_nearest_shapes;
    use crate::query::norm;
    use crate::query::point_containment_all_pairs;
    use crate::query::project_point_convex;
    use crate::query::ray_cast_all_pairs;
    use crate::query::sphere_aabb_overlap;
    use crate::query::sphere_cast_scene;
    use crate::query::sphere_sphere_overlap;
    use crate::query::sub;
    #[test]
    fn aabb_overlap_true_when_touching() {
        assert!(aabb_overlap(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 1.0],
        ));
    }
    #[test]
    fn aabb_overlap_false_when_separated() {
        assert!(!aabb_overlap(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 1.0],
        ));
    }
    #[test]
    fn sphere_sphere_overlap_touching() {
        assert!(sphere_sphere_overlap(
            [0.0, 0.0, 0.0],
            1.0,
            [2.0, 0.0, 0.0],
            1.0
        ));
    }
    #[test]
    fn sphere_sphere_overlap_separated() {
        assert!(!sphere_sphere_overlap(
            [0.0, 0.0, 0.0],
            0.4,
            [2.0, 0.0, 0.0],
            0.4
        ));
    }
    #[test]
    fn sphere_aabb_overlap_inside() {
        assert!(sphere_aabb_overlap(
            [0.0, 0.0, 0.0],
            0.5,
            [0.0, 0.0, 0.0],
            [2.0, 2.0, 2.0]
        ));
    }
    #[test]
    fn sphere_aabb_overlap_separated() {
        assert!(!sphere_aabb_overlap(
            [5.0, 0.0, 0.0],
            0.5,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0]
        ));
    }
    #[test]
    fn sphere_aabb_overlap_touching_face() {
        assert!(sphere_aabb_overlap(
            [1.5, 0.0, 0.0],
            0.5,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0]
        ));
    }
    #[test]
    fn capsule_capsule_overlap_crossing() {
        assert!(capsule_capsule_overlap(
            [-2.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            0.1,
            [0.0, -2.0, 0.0],
            [0.0, 2.0, 0.0],
            0.1,
        ));
    }
    #[test]
    fn capsule_capsule_overlap_separated() {
        assert!(!capsule_capsule_overlap(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.1,
            [-1.0, 5.0, 0.0],
            [1.0, 5.0, 0.0],
            0.1,
        ));
    }
    #[test]
    fn all_pairs_ray_cast_hits_all() {
        let s1 = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let s2 = QuerySphere {
            center: [5.0, 0.0, 0.0],
            radius: 1.0,
        };
        let shapes = vec![QueryShapeRef::Sphere(&s1), QueryShapeRef::Sphere(&s2)];
        let hits = ray_cast_all_pairs([-10.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, &shapes);
        assert_eq!(hits.len(), 2, "both spheres should be hit");
        assert!(hits[0].toi < hits[1].toi);
    }
    #[test]
    fn all_pairs_ray_cast_misses_all() {
        let s1 = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
        };
        let shapes = vec![QueryShapeRef::Sphere(&s1)];
        let hits = ray_cast_all_pairs([0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 100.0, &shapes);
        assert!(hits.is_empty(), "ray above sphere should miss");
    }
    #[test]
    fn all_pairs_point_containment() {
        let s = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 2.0,
        };
        let b = QueryBox {
            center: [10.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let shapes = vec![QueryShapeRef::Sphere(&s), QueryShapeRef::Box(&b)];
        let inside = point_containment_all_pairs([0.5, 0.0, 0.0], &shapes);
        assert!(inside.contains(&0), "point inside sphere");
        assert!(!inside.contains(&1), "point not inside box");
    }
    #[test]
    fn all_pairs_point_containment_both() {
        let s = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let b = QueryBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [5.0, 5.0, 5.0],
        };
        let shapes = vec![QueryShapeRef::Sphere(&s), QueryShapeRef::Box(&b)];
        let inside = point_containment_all_pairs([0.5, 0.0, 0.0], &shapes);
        assert_eq!(inside.len(), 2, "point inside both shapes");
    }
    fn make_unit_frustum() -> ViewFrustum {
        ViewFrustum::new([
            [1.0, 0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ])
    }
    #[test]
    fn frustum_test_sphere_inside() {
        let f = make_unit_frustum();
        assert!(
            f.test_sphere([0.0, 0.0, 0.0], 0.1),
            "sphere at origin should be visible"
        );
    }
    #[test]
    fn frustum_test_sphere_outside() {
        let f = make_unit_frustum();
        assert!(
            !f.test_sphere([5.0, 0.0, 0.0], 0.1),
            "sphere far away should be culled"
        );
    }
    #[test]
    fn frustum_test_sphere_straddling() {
        let f = make_unit_frustum();
        assert!(
            f.test_sphere([1.5, 0.0, 0.0], 1.0),
            "straddling sphere should not be culled"
        );
    }
    #[test]
    fn frustum_test_aabb_inside() {
        let f = make_unit_frustum();
        assert!(f.test_aabb([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]));
    }
    #[test]
    fn frustum_test_aabb_outside() {
        let f = make_unit_frustum();
        assert!(!f.test_aabb([3.0, 3.0, 3.0], [4.0, 4.0, 4.0]));
    }
    #[test]
    fn frustum_query_filters_correctly() {
        let f = make_unit_frustum();
        let s_in = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.1,
        };
        let s_out = QuerySphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.1,
        };
        let b_in = QueryBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![
            QueryShapeRef::Sphere(&s_in),
            QueryShapeRef::Sphere(&s_out),
            QueryShapeRef::Box(&b_in),
        ];
        let visible = frustum_query(&f, &shapes);
        assert!(visible.contains(&0), "in-frustum sphere should be visible");
        assert!(
            !visible.contains(&1),
            "out-of-frustum sphere should be culled"
        );
        assert!(visible.contains(&2), "in-frustum box should be visible");
    }
    #[test]
    fn k_nearest_returns_closest_first() {
        let s_near = QuerySphere {
            center: [1.0, 0.0, 0.0],
            radius: 0.1,
        };
        let s_far = QuerySphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.1,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![
            QueryShapeRef::Sphere(&s_far),
            QueryShapeRef::Sphere(&s_near),
        ];
        let nearest = k_nearest_shapes([0.0, 0.0, 0.0], &shapes, 2);
        assert_eq!(nearest.len(), 2);
        assert!(
            nearest[0].distance < nearest[1].distance,
            "nearest should be first"
        );
        assert_eq!(nearest[0].index, 1, "index 1 (near sphere) should be first");
    }
    #[test]
    fn k_nearest_k_larger_than_shapes() {
        let s = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Sphere(&s)];
        let nearest = k_nearest_shapes([5.0, 0.0, 0.0], &shapes, 100);
        assert_eq!(nearest.len(), 1, "should only return what exists");
    }
    #[test]
    fn k_nearest_inside_shape_is_negative_distance() {
        let s = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 5.0,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Sphere(&s)];
        let nearest = k_nearest_shapes([0.0, 0.0, 0.0], &shapes, 1);
        assert!(
            nearest[0].distance < 0.0,
            "inside sphere → negative distance"
        );
    }
    #[test]
    fn k_nearest_zero_k_returns_empty() {
        let s = QuerySphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Sphere(&s)];
        let nearest = k_nearest_shapes([5.0, 0.0, 0.0], &shapes, 0);
        assert!(nearest.is_empty());
    }
    #[test]
    fn sphere_cast_hits_sphere() {
        let target = QuerySphere {
            center: [5.0, 0.0, 0.0],
            radius: 1.0,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Sphere(&target)];
        let result = sphere_cast_scene([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 10.0, &shapes);
        assert!(result.is_some(), "sphere cast should hit expanded sphere");
        let r = result.unwrap();
        assert_eq!(r.shape_index, 0);
        assert!(r.toi > 0.0 && r.toi <= 10.0);
    }
    #[test]
    fn sphere_cast_hits_box() {
        let target = QueryBox {
            center: [5.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Box(&target)];
        let result = sphere_cast_scene([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 10.0, &shapes);
        assert!(result.is_some(), "sphere cast should hit box");
    }
    #[test]
    fn sphere_cast_misses() {
        let target = QuerySphere {
            center: [0.0, 5.0, 0.0],
            radius: 0.1,
        };
        let shapes: Vec<QueryShapeRef<'_>> = vec![QueryShapeRef::Sphere(&target)];
        let result = sphere_cast_scene([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 100.0, &shapes);
        assert!(result.is_none(), "sphere cast perpendicular should miss");
    }
    #[test]
    fn sphere_cast_finds_nearest_of_two() {
        let near = QuerySphere {
            center: [3.0, 0.0, 0.0],
            radius: 0.5,
        };
        let far = QuerySphere {
            center: [8.0, 0.0, 0.0],
            radius: 0.5,
        };
        let shapes: Vec<QueryShapeRef<'_>> =
            vec![QueryShapeRef::Sphere(&far), QueryShapeRef::Sphere(&near)];
        let result = sphere_cast_scene([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 100.0, &shapes);
        assert!(result.is_some());
        let r = result.unwrap();
        assert_eq!(
            r.shape_index, 1,
            "should hit the nearer sphere first (index 1)"
        );
    }
    #[test]
    fn closest_point_sphere_outside() {
        let p = closest_point_on_sphere([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        assert!((p[0] - 1.0).abs() < 1e-10, "closest point should be at x=1");
    }
    #[test]
    fn closest_point_box_outside() {
        let p = closest_point_on_box([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn closest_point_capsule_side() {
        let p = closest_point_on_capsule([2.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!(
            (p[0] - 0.5).abs() < 1e-10,
            "closest x should be at capsule surface"
        );
        assert!(
            p[1].abs() < 1e-10,
            "closest point should be at midpoint height"
        );
    }
    #[test]
    fn toi_spheres_approaching_head_on() {
        let result = compute_time_of_impact_spheres(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            [5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            1.0,
            100.0,
        );
        assert!(result.is_some(), "approaching spheres must collide");
        let r = result.unwrap();
        assert!(
            (r.toi - 4.0).abs() < 1e-10,
            "toi should be 4.0, got {}",
            r.toi
        );
    }
    #[test]
    fn toi_spheres_no_collision_diverging() {
        let result = compute_time_of_impact_spheres(
            [-5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            1.0,
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            100.0,
        );
        assert!(result.is_none(), "diverging spheres should not collide");
    }
    #[test]
    fn toi_aabbs_approaching() {
        let result = compute_time_of_impact_aabbs(
            [-6.0, -1.0, -1.0],
            [-4.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [4.0, -1.0, -1.0],
            [6.0, 1.0, 1.0],
            [-1.0, 0.0, 0.0],
            100.0,
        );
        assert!(result.is_some(), "approaching boxes must collide");
        let toi = result.unwrap();
        assert!((toi - 4.0).abs() < 1e-10, "toi should be 4.0, got {}", toi);
    }
    #[test]
    fn toi_aabbs_no_collision_miss() {
        let result = compute_time_of_impact_aabbs(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 0.0],
            [5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0],
            [0.0, 1.0, 0.0],
            100.0,
        );
        assert!(
            result.is_none(),
            "parallel-moving separated boxes should not collide"
        );
    }
    #[test]
    fn project_point_convex_square_outside() {
        let verts = vec![
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let p = project_point_convex([3.0, 0.0, 0.0], &verts);
        let d = norm(sub(p, [1.0, 0.0, 0.0]));
        assert!(d < 0.5, "closest point should be near (1,0,0), got {:?}", p);
    }
    #[test]
    fn project_point_convex_single_vertex() {
        let p = project_point_convex([5.0, 5.0, 5.0], &[[1.0, 2.0, 3.0]]);
        assert_eq!(
            p,
            [1.0, 2.0, 3.0],
            "single vertex: project returns that vertex"
        );
    }
    #[test]
    fn segment_distance_perpendicular_segments() {
        let r = distance_between_segments(
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [2.0, -2.0, 0.0],
            [2.0, 2.0, 0.0],
        );
        assert!(
            r.distance < 1e-10,
            "crossing segments: distance should be ~0, got {}",
            r.distance
        );
    }
    #[test]
    fn segment_distance_parallel_segments() {
        let r = distance_between_segments(
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [4.0, 3.0, 0.0],
        );
        assert!(
            (r.distance - 3.0).abs() < 1e-10,
            "parallel segments 3 apart: distance should be 3, got {}",
            r.distance
        );
    }
    #[test]
    fn segment_distance_disjoint_segments() {
        let r = distance_between_segments(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
        );
        assert!(
            (r.distance - 4.0).abs() < 1e-10,
            "disjoint segments: distance should be 4, got {}",
            r.distance
        );
    }
    #[test]
    fn segment_distance_touching_at_endpoint() {
        let r = distance_between_segments(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        );
        assert!(
            r.distance < 1e-10,
            "touching segments: distance must be 0, got {}",
            r.distance
        );
    }
}
