//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    Aabb, BatchRayResult, Capsule, ConvexMesh, Heightfield, Obb, Ray, RayHit, Sphere,
};

#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn len_sq(v: [f64; 3]) -> f64 {
    dot(v, v)
}
#[inline]
pub(super) fn len(v: [f64; 3]) -> f64 {
    len_sq(v).sqrt()
}
#[inline]
pub(super) fn normalize(v: [f64; 3]) -> [f64; 3] {
    let l = len(v);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        scale(v, 1.0 / l)
    }
}
#[inline]
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn neg(v: [f64; 3]) -> [f64; 3] {
    [-v[0], -v[1], -v[2]]
}
#[inline]
pub(super) fn reflect_vec(d: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    sub(d, scale(n, 2.0 * dot(d, n)))
}
/// Apply a rotation matrix (row-major 3×3) to a vector.
#[inline]
pub(super) fn mat3_mul_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}
/// Transpose a 3×3 rotation matrix.
#[inline]
pub(super) fn mat3_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Intersect a ray with a sphere.
///
/// Uses the geometric / discriminant method.  Returns the nearest hit
/// within `[ray.t_min, ray.t_max]`, or `None`.
pub fn ray_sphere(ray: &Ray, sphere: &Sphere) -> Option<RayHit> {
    let oc = sub(ray.origin, sphere.centre);
    let a = len_sq(ray.direction);
    let half_b = dot(oc, ray.direction);
    let c = len_sq(oc) - sphere.radius * sphere.radius;
    let disc = half_b * half_b - a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_d = disc.sqrt();
    let t0 = (-half_b - sqrt_d) / a;
    let t1 = (-half_b + sqrt_d) / a;
    let t = if ray.valid_t(t0) {
        t0
    } else if ray.valid_t(t1) {
        t1
    } else {
        return None;
    };
    let point = ray.at(t);
    let normal = sphere.normal_at(point);
    let uv = sphere.uv_at(point);
    Some(RayHit::new(t, point, normal, uv, 0))
}
/// Intersect a ray with an AABB using the slab method.
///
/// Returns the entry and exit `t` values, or `None` if the ray misses.
pub fn ray_aabb(ray: &Ray, aabb: &Aabb) -> Option<(f64, f64)> {
    let inv_dir = [
        1.0 / ray.direction[0],
        1.0 / ray.direction[1],
        1.0 / ray.direction[2],
    ];
    let mut t_min = ray.t_min;
    let mut t_max = ray.t_max;
    for i in 0..3 {
        let t1 = (aabb.min[i] - ray.origin[i]) * inv_dir[i];
        let t2 = (aabb.max[i] - ray.origin[i]) * inv_dir[i];
        let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        t_min = t_min.max(ta);
        t_max = t_max.min(tb);
        if t_min > t_max {
            return None;
        }
    }
    Some((t_min, t_max))
}
/// Intersect a ray with an AABB and return a full [`RayHit`] with normal.
///
/// The normal points outward from the face that was first hit.
pub fn ray_aabb_hit(ray: &Ray, aabb: &Aabb, prim_id: usize) -> Option<RayHit> {
    let (t_entry, _t_exit) = ray_aabb(ray, aabb)?;
    if !ray.valid_t(t_entry) {
        return None;
    }
    let point = ray.at(t_entry);
    let c = aabb.centre();
    let h = aabb.half_extents();
    let local = sub(point, c);
    let mut normal = [0.0f64; 3];
    let mut min_diff = f64::INFINITY;
    for i in 0..3 {
        let diff = (local[i].abs() - h[i]).abs();
        if diff < min_diff {
            min_diff = diff;
            normal = [0.0; 3];
            normal[i] = if local[i] >= 0.0 { 1.0 } else { -1.0 };
        }
    }
    Some(RayHit::new(t_entry, point, normal, [0.0, 0.0], prim_id))
}
/// Intersect a ray with an OBB.
///
/// Transforms the ray into OBB local space and applies the slab method.
pub fn ray_obb(ray: &Ray, obb: &Obb) -> Option<RayHit> {
    let d = sub(ray.origin, obb.centre);
    let axes_t = mat3_transpose(&obb.axes);
    let local_origin = mat3_mul_vec(&axes_t, d);
    let local_dir = mat3_mul_vec(&axes_t, ray.direction);
    let local_ray = Ray {
        origin: local_origin,
        direction: local_dir,
        t_min: ray.t_min,
        t_max: ray.t_max,
    };
    let local_aabb = Aabb::new(neg(obb.half_extents), obb.half_extents);
    let (t_entry, _t_exit) = ray_aabb(&local_ray, &local_aabb)?;
    if !ray.valid_t(t_entry) {
        return None;
    }
    let local_point = local_ray.at(t_entry);
    let world_point = add(obb.centre, mat3_mul_vec(&obb.axes, local_point));
    let mut normal_local = [0.0f64; 3];
    let mut min_diff = f64::INFINITY;
    for i in 0..3 {
        let diff = (local_point[i].abs() - obb.half_extents[i]).abs();
        if diff < min_diff {
            min_diff = diff;
            normal_local = [0.0; 3];
            normal_local[i] = if local_point[i] >= 0.0 { 1.0 } else { -1.0 };
        }
    }
    let world_normal = normalize(mat3_mul_vec(&obb.axes, normal_local));
    Some(RayHit::new(
        t_entry,
        world_point,
        world_normal,
        [0.0, 0.0],
        0,
    ))
}
/// Intersect a ray with a triangle using the Möller-Trumbore algorithm.
///
/// # Arguments
/// * `ray`  – The ray.
/// * `v0, v1, v2` – Triangle vertices.
/// * `prim_id`    – Index of this primitive for the hit record.
///
/// Returns `Some(RayHit)` with barycentric UVs `(u, v)` if the ray hits,
/// `None` otherwise.
pub fn ray_triangle(
    ray: &Ray,
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    prim_id: usize,
) -> Option<RayHit> {
    let edge1 = sub(v1, v0);
    let edge2 = sub(v2, v0);
    let h = cross(ray.direction, edge2);
    let a = dot(edge1, h);
    pub(super) const EPS: f64 = 1e-10;
    if a.abs() < EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = sub(ray.origin, v0);
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, edge1);
    let v = f * dot(ray.direction, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot(edge2, q);
    if !ray.valid_t(t) {
        return None;
    }
    let point = ray.at(t);
    let normal = normalize(cross(edge1, edge2));
    Some(RayHit::new(t, point, normal, [u, v], prim_id))
}
/// Intersect a ray with a triangle, respecting back-face culling.
///
/// Returns `None` if the ray hits the back face.
pub fn ray_triangle_cull(
    ray: &Ray,
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    prim_id: usize,
) -> Option<RayHit> {
    let hit = ray_triangle(ray, v0, v1, v2, prim_id)?;
    if dot(hit.normal, ray.direction) >= 0.0 {
        None
    } else {
        Some(hit)
    }
}
/// Intersect a ray with a capsule.
///
/// The implementation tests the infinite cylinder and clamps to the
/// hemispherical endcaps.
pub fn ray_capsule(ray: &Ray, capsule: &Capsule) -> Option<RayHit> {
    let ab = capsule.axis();
    let ao = sub(ray.origin, capsule.a);
    let ab_len_sq = len_sq(ab);
    if ab_len_sq < 1e-300 {
        return ray_sphere(ray, &Sphere::new(capsule.a, capsule.radius));
    }
    let d = ray.direction;
    let ab_dot_d = dot(ab, d);
    let ab_dot_ao = dot(ab, ao);
    let a_coef = len_sq(d) - ab_dot_d * ab_dot_d / ab_len_sq;
    let b_coef = dot(ao, d) - ab_dot_ao * ab_dot_d / ab_len_sq;
    let c_coef = len_sq(ao) - ab_dot_ao * ab_dot_ao / ab_len_sq - capsule.radius * capsule.radius;
    let mut best_t = f64::INFINITY;
    if a_coef.abs() > 1e-12 {
        let disc = b_coef * b_coef - a_coef * c_coef;
        if disc >= 0.0 {
            let sqrt_d = disc.sqrt();
            for sign in [-1.0f64, 1.0] {
                let t = (-b_coef + sign * sqrt_d) / a_coef;
                if !ray.valid_t(t) {
                    continue;
                }
                let q = add(ray.origin, scale(d, t));
                let proj = dot(sub(q, capsule.a), ab) / ab_len_sq;
                if (0.0..=1.0).contains(&proj) && t < best_t {
                    best_t = t;
                }
            }
        }
    }
    for cap_centre in [capsule.a, capsule.b] {
        let cap_sphere = Sphere::new(cap_centre, capsule.radius);
        if let Some(hit) = ray_sphere(ray, &cap_sphere)
            && hit.t < best_t
        {
            best_t = hit.t;
        }
    }
    if best_t.is_finite() {
        let point = ray.at(best_t);
        let proj_t = dot(sub(point, capsule.a), ab) / ab_len_sq;
        let seg_pt = add(capsule.a, scale(ab, proj_t.clamp(0.0, 1.0)));
        let normal = normalize(sub(point, seg_pt));
        Some(RayHit::new(best_t, point, normal, [0.0, 0.0], 0))
    } else {
        None
    }
}
/// Intersect a ray with a convex mesh by ray-marching along the ray and
/// testing the AABB first, then performing binary search on the penetration
/// depth (simplified approach for a convex hull).
///
/// For correctness, the mesh faces must be provided as triangles.
/// This function tests all triangles of a triangle-soup convex mesh.
///
/// # Arguments
/// * `ray`       – The ray.
/// * `mesh`      – Convex mesh.
/// * `triangles` – Index triples `(i, j, k)` referencing `mesh.vertices`.
pub fn ray_convex_mesh(
    ray: &Ray,
    mesh: &ConvexMesh,
    triangles: &[(usize, usize, usize)],
) -> Option<RayHit> {
    let aabb = mesh.aabb();
    ray_aabb(ray, &aabb)?;
    let mut best: Option<RayHit> = None;
    for &(i, j, k) in triangles {
        if i >= mesh.vertices.len() || j >= mesh.vertices.len() || k >= mesh.vertices.len() {
            continue;
        }
        let v0 = mesh.vertices[i];
        let v1 = mesh.vertices[j];
        let v2 = mesh.vertices[k];
        if let Some(hit) = ray_triangle(ray, v0, v1, v2, i) {
            best = Some(match best {
                None => hit,
                Some(prev) => {
                    if hit.t < prev.t {
                        hit
                    } else {
                        prev
                    }
                }
            });
        }
    }
    best
}
/// Intersect a ray with a heightfield by stepping along the ray in XZ and
/// checking whether the ray drops below the terrain surface.
///
/// Returns `None` if no intersection is found within `max_steps` steps.
#[allow(clippy::too_many_arguments)]
pub fn ray_heightfield(ray: &Ray, hf: &Heightfield, max_steps: usize) -> Option<RayHit> {
    let step = hf.cell_size * 0.5;
    let mut t = ray.t_min;
    let mut prev_above = true;
    let mut prev_t = t;
    for _step_i in 0..max_steps {
        let p = ray.at(t);
        let x = p[0];
        let z = p[2];
        if x < 0.0 || z < 0.0 || x > hf.width() || z > hf.depth() {
            t += step;
            continue;
        }
        let terrain_y = hf.height_at(x, z);
        let above = p[1] >= terrain_y;
        if !above && prev_above {
            let mut ta = prev_t;
            let mut tb = t;
            for _ in 0..20 {
                let tm = 0.5 * (ta + tb);
                let pm = ray.at(tm);
                let ym = hf.height_at(pm[0], pm[2]);
                if pm[1] >= ym {
                    ta = tm;
                } else {
                    tb = tm;
                }
            }
            let t_hit = 0.5 * (ta + tb);
            let p_hit = ray.at(t_hit);
            let eps = hf.cell_size * 0.01;
            let dh_dx = (hf.height_at(p_hit[0] + eps, p_hit[2])
                - hf.height_at(p_hit[0] - eps, p_hit[2]))
                / (2.0 * eps);
            let dh_dz = (hf.height_at(p_hit[0], p_hit[2] + eps)
                - hf.height_at(p_hit[0], p_hit[2] - eps))
                / (2.0 * eps);
            let normal = normalize([-dh_dx, 1.0, -dh_dz]);
            return Some(RayHit::new(
                t_hit,
                p_hit,
                normal,
                [p_hit[0] / hf.width(), p_hit[2] / hf.depth()],
                0,
            ));
        }
        prev_above = above;
        prev_t = t;
        t += step;
        if t > ray.t_max {
            break;
        }
    }
    None
}
/// Cast multiple rays against a list of spheres, returning the closest hit
/// per ray.
///
/// This is the brute-force O(N·M) reference implementation.
pub fn batch_ray_spheres(rays: &[Ray], spheres: &[Sphere]) -> BatchRayResult {
    let hits = rays
        .iter()
        .map(|ray| {
            let mut best: Option<RayHit> = None;
            for sphere in spheres {
                if let Some(hit) = ray_sphere(ray, sphere) {
                    best = Some(match best {
                        None => hit,
                        Some(prev) => {
                            if hit.t < prev.t {
                                hit
                            } else {
                                prev
                            }
                        }
                    });
                }
            }
            best
        })
        .collect();
    BatchRayResult { hits }
}
/// Cast multiple rays against a list of triangles (triangle mesh).
///
/// Returns the closest hit per ray across all triangles.
pub fn batch_ray_triangles(
    rays: &[Ray],
    vertices: &[[f64; 3]],
    indices: &[(usize, usize, usize)],
) -> BatchRayResult {
    let hits = rays
        .iter()
        .map(|ray| {
            let mut best: Option<RayHit> = None;
            for (fi, &(i, j, k)) in indices.iter().enumerate() {
                if i >= vertices.len() || j >= vertices.len() || k >= vertices.len() {
                    continue;
                }
                if let Some(hit) = ray_triangle(ray, vertices[i], vertices[j], vertices[k], fi) {
                    best = Some(match best {
                        None => hit,
                        Some(prev) => {
                            if hit.t < prev.t {
                                hit
                            } else {
                                prev
                            }
                        }
                    });
                }
            }
            best
        })
        .collect();
    BatchRayResult { hits }
}
/// Unproject a screen-space pixel into a world-space ray for object picking.
///
/// # Arguments
/// * `screen_x, screen_y` – Pixel coordinates.
/// * `screen_width, screen_height` – Screen dimensions.
/// * `fov_y_rad`          – Vertical field of view in radians.
/// * `near`               – Near clip plane distance.
/// * `cam_pos`            – Camera position in world space.
/// * `cam_forward`        – Camera forward direction (normalised).
/// * `cam_up`             – Camera up direction (normalised).
#[allow(clippy::too_many_arguments)]
pub fn pick_ray(
    screen_x: f64,
    screen_y: f64,
    screen_width: f64,
    screen_height: f64,
    fov_y_rad: f64,
    _near: f64,
    cam_pos: [f64; 3],
    cam_forward: [f64; 3],
    cam_up: [f64; 3],
) -> Ray {
    let aspect = screen_width / screen_height.max(1e-10);
    let half_h = (fov_y_rad * 0.5).tan();
    let half_w = aspect * half_h;
    let ndc_x = (screen_x / screen_width.max(1e-10)) * 2.0 - 1.0;
    let ndc_y = 1.0 - (screen_y / screen_height.max(1e-10)) * 2.0;
    let cam_right = normalize(cross(cam_forward, cam_up));
    let cam_up_ortho = normalize(cross(cam_right, cam_forward));
    let dir = normalize(add(
        add(cam_forward, scale(cam_right, ndc_x * half_w)),
        scale(cam_up_ortho, ndc_y * half_h),
    ));
    Ray::new(cam_pos, dir)
}
/// Test whether the path between `origin` and `light_pos` is occluded by
/// any sphere in `occluders`.
///
/// Returns `true` if the point is in shadow (a sphere blocks the path).
pub fn is_in_shadow(origin: [f64; 3], light_pos: [f64; 3], occluders: &[Sphere]) -> bool {
    let to_light = sub(light_pos, origin);
    let dist = len(to_light);
    if dist < 1e-12 {
        return false;
    }
    let shadow_ray = Ray {
        origin,
        direction: normalize(to_light),
        t_min: 1e-4,
        t_max: dist - 1e-4,
    };
    occluders
        .iter()
        .any(|s| ray_sphere(&shadow_ray, s).is_some())
}
/// Cast a shadow ray toward an area light and estimate the shadow attenuation
/// by sampling multiple points on the light disc.
///
/// Returns the fraction of light samples that are not occluded (`0.0` = fully
/// in shadow, `1.0` = fully lit).
pub fn soft_shadow_factor(
    origin: [f64; 3],
    light_centre: [f64; 3],
    light_radius: f64,
    occluders: &[Sphere],
    n_samples: usize,
    seed: u64,
) -> f64 {
    let mut rng_state = seed.max(1);
    let mut lcg = || -> f64 {
        rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((rng_state >> 11) as f64) * (1.0 / (1u64 << 53) as f64)
    };
    let to_light = normalize(sub(light_centre, origin));
    let up = if to_light[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let right = normalize(cross(to_light, up));
    let up2 = cross(to_light, right);
    let mut lit = 0usize;
    for _ in 0..n_samples {
        let u = lcg() * 2.0 - 1.0;
        let v = lcg() * 2.0 - 1.0;
        let offset = add(scale(right, u * light_radius), scale(up2, v * light_radius));
        let sample = add(light_centre, offset);
        if !is_in_shadow(origin, sample, occluders) {
            lit += 1;
        }
    }
    lit as f64 / n_samples.max(1) as f64
}
/// Compute the reflected ray direction.
///
/// `d_reflect = d - 2(d·n)n` where `n` is the outward unit normal.
pub fn reflect_ray(incident: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    normalize(reflect_vec(incident, normal))
}
/// Compute the refracted ray direction using Snell's law.
///
/// # Arguments
/// * `incident`  – Incident ray direction (unit vector, pointing toward surface).
/// * `normal`    – Surface normal (unit vector, outward).
/// * `n1`        – Index of refraction of the incident medium.
/// * `n2`        – Index of refraction of the transmitted medium.
///
/// Returns `None` on total internal reflection.
pub fn refract_ray(incident: [f64; 3], normal: [f64; 3], n1: f64, n2: f64) -> Option<[f64; 3]> {
    let n_ratio = n1 / n2.max(1e-300);
    let cos_i = -dot(incident, normal);
    let sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        return None;
    }
    let cos_t = (1.0 - sin2_t).sqrt();
    let refracted = add(
        scale(incident, n_ratio),
        scale(normal, n_ratio * cos_i - cos_t),
    );
    Some(normalize(refracted))
}
/// Schlick's approximation for the Fresnel reflectance.
///
/// # Arguments
/// * `cos_theta` – Cosine of the angle of incidence.
/// * `n1`        – Refractive index of incident medium.
/// * `n2`        – Refractive index of transmitted medium.
pub fn schlick_reflectance(cos_theta: f64, n1: f64, n2: f64) -> f64 {
    let r0 = ((n1 - n2) / (n1 + n2)).powi(2);
    r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
}
/// Compute the closest point on a ray to a query point `q`.
///
/// Returns `(t, closest_point)` where `t ≥ ray.t_min`.
pub fn ray_closest_point(ray: &Ray, q: [f64; 3]) -> (f64, [f64; 3]) {
    let t = dot(sub(q, ray.origin), ray.direction).max(ray.t_min);
    let t_clamped = t.min(ray.t_max);
    (t_clamped, ray.at(t_clamped))
}
/// Squared distance from a point `q` to the nearest point on the ray.
pub fn ray_point_dist_sq(ray: &Ray, q: [f64; 3]) -> f64 {
    let (_t, closest) = ray_closest_point(ray, q);
    len_sq(sub(q, closest))
}
/// Compute the distance between two rays (shortest segment between infinite lines,
/// clamped to the valid `[t_min, t_max]` range of each ray).
///
/// Returns `(t1, t2, distance)`.
pub fn ray_ray_distance(r1: &Ray, r2: &Ray) -> (f64, f64, f64) {
    let w = sub(r1.origin, r2.origin);
    let a = dot(r1.direction, r1.direction);
    let b = dot(r1.direction, r2.direction);
    let c = dot(r2.direction, r2.direction);
    let d = dot(r1.direction, w);
    let e = dot(r2.direction, w);
    let denom = a * c - b * b;
    let (t1, t2) = if denom.abs() < 1e-12 {
        let t1 = r1.t_min;
        let t2 = (dot(sub(r1.origin, r2.origin), r2.direction) / c.max(1e-300))
            .clamp(r2.t_min, r2.t_max);
        (t1, t2)
    } else {
        let t1 = ((b * e - c * d) / denom).clamp(r1.t_min, r1.t_max);
        let t2 = ((a * e - b * d) / denom).clamp(r2.t_min, r2.t_max);
        (t1, t2)
    };
    let dist = len(sub(r1.at(t1), r2.at(t2)));
    (t1, t2, dist)
}
/// Closest distance from a point to a line segment.
///
/// Returns `(t, point_on_segment, distance)`.
pub fn segment_point_distance(a: [f64; 3], b: [f64; 3], q: [f64; 3]) -> (f64, [f64; 3], f64) {
    let ab = sub(b, a);
    let len_ab_sq = len_sq(ab);
    if len_ab_sq < 1e-300 {
        return (0.0, a, len(sub(q, a)));
    }
    let t = (dot(sub(q, a), ab) / len_ab_sq).clamp(0.0, 1.0);
    let p = add(a, scale(ab, t));
    (t, p, len(sub(q, p)))
}
