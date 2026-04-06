// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU ray tracing with f32 precision (CPU mock).
//!
//! Provides BVH nodes, triangle/sphere/AABB intersection tests,
//! full scene traversal, ambient occlusion sampling, and a simple
//! orthographic renderer — all running on the CPU as a mock GPU backend.

// ── Core types ────────────────────────────────────────────────────────────────

/// A ray defined by an origin and a direction (f32 precision).
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f32; 3],
    /// Ray direction (should be normalised for physically correct distances).
    pub direction: [f32; 3],
}

impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f32; 3], direction: [f32; 3]) -> Self {
        Self { origin, direction }
    }

    /// Evaluate the ray at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: f32) -> [f32; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}

/// Information about a ray–surface intersection.
#[derive(Debug, Clone, Copy)]
pub struct HitRecord {
    /// Ray parameter at the hit point.
    pub t: f32,
    /// World-space hit point.
    pub point: [f32; 3],
    /// Surface normal at the hit point (outward facing).
    pub normal: [f32; 3],
    /// Index into the material table.
    pub material_id: u32,
}

/// A BVH (bounding-volume hierarchy) node.
#[derive(Debug, Clone, Copy)]
pub struct BvhNode {
    /// Minimum corner of the axis-aligned bounding box.
    pub aabb_min: [f32; 3],
    /// Maximum corner of the axis-aligned bounding box.
    pub aabb_max: [f32; 3],
    /// Index of the left child (used when `is_leaf == false`).
    pub left: u32,
    /// Index of the right child (used when `is_leaf == false`).
    pub right: u32,
    /// True if this node stores a triangle directly.
    pub is_leaf: bool,
    /// Triangle index stored by this leaf node.
    pub tri_idx: u32,
}

/// A triangle in GPU-friendly format.
#[derive(Debug, Clone, Copy)]
pub struct GpuTriangle {
    /// First vertex.
    pub v0: [f32; 3],
    /// Second vertex.
    pub v1: [f32; 3],
    /// Third vertex.
    pub v2: [f32; 3],
    /// Pre-computed face normal.
    pub normal: [f32; 3],
    /// Index into the material table.
    pub material_id: u32,
}

// ── Vector helpers (f32) ──────────────────────────────────────────────────────

#[inline]
fn dot3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn sub3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn normalize3f(v: [f32; 3]) -> [f32; 3] {
    let len = dot3f(v, v).sqrt();
    if len < 1e-10 {
        return [0.0; 3];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

// ── Intersection tests ────────────────────────────────────────────────────────

/// Test ray–sphere intersection.
///
/// Returns the smallest positive `t`, or `None` if there is no hit.
pub fn ray_sphere_intersect(ray: &Ray, center: [f32; 3], radius: f32) -> Option<f32> {
    let oc = sub3f(ray.origin, center);
    let a = dot3f(ray.direction, ray.direction);
    let half_b = dot3f(oc, ray.direction);
    let c = dot3f(oc, oc) - radius * radius;
    let discriminant = half_b * half_b - a * c;
    if discriminant < 0.0 {
        return None;
    }
    let sqrt_d = discriminant.sqrt();
    let t1 = (-half_b - sqrt_d) / a;
    if t1 > 1e-4 {
        return Some(t1);
    }
    let t2 = (-half_b + sqrt_d) / a;
    if t2 > 1e-4 { Some(t2) } else { None }
}

/// Test ray–triangle intersection using the Möller–Trumbore algorithm.
///
/// Returns `Some(t)` on a valid hit, or `None`.
pub fn ray_triangle_intersect(ray: &Ray, tri: &GpuTriangle) -> Option<f32> {
    let edge1 = sub3f(tri.v1, tri.v0);
    let edge2 = sub3f(tri.v2, tri.v0);
    let h = cross3f(ray.direction, edge2);
    let a = dot3f(edge1, h);
    if a.abs() < 1e-8 {
        return None; // Ray is parallel to triangle
    }
    let f = 1.0 / a;
    let s = sub3f(ray.origin, tri.v0);
    let u = f * dot3f(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross3f(s, edge1);
    let v = f * dot3f(ray.direction, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot3f(edge2, q);
    if t > 1e-4 { Some(t) } else { None }
}

/// Test ray–AABB intersection using the slab method.
///
/// Returns `true` if the ray hits the box `[aabb_min, aabb_max]`.
pub fn ray_aabb_intersect(ray: &Ray, aabb_min: [f32; 3], aabb_max: [f32; 3]) -> bool {
    let mut t_min = 0.0_f32;
    let mut t_max = f32::MAX;
    for i in 0..3 {
        let inv_d = 1.0 / ray.direction[i];
        let t0 = (aabb_min[i] - ray.origin[i]) * inv_d;
        let t1 = (aabb_max[i] - ray.origin[i]) * inv_d;
        let (t_near, t_far) = if inv_d >= 0.0 { (t0, t1) } else { (t1, t0) };
        t_min = t_min.max(t_near);
        t_max = t_max.min(t_far);
        if t_max < t_min {
            return false;
        }
    }
    t_max >= 0.0
}

// ── Scene traversal ───────────────────────────────────────────────────────────

/// Find the closest triangle hit by `ray` in `triangles`.
///
/// Returns `Some(HitRecord)` if any triangle is hit, `None` otherwise.
pub fn trace_ray(ray: &Ray, triangles: &[GpuTriangle]) -> Option<HitRecord> {
    let mut best_t = f32::MAX;
    let mut best_hit: Option<HitRecord> = None;

    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(ray, tri)
            && t < best_t
        {
            best_t = t;
            let point = ray.at(t);
            best_hit = Some(HitRecord {
                t,
                point,
                normal: tri.normal,
                material_id: tri.material_id,
            });
        }
    }
    best_hit
}

// ── Rendering ─────────────────────────────────────────────────────────────────

/// Render a scene with orthographic projection, returning one RGB pixel per
/// sample.
///
/// * `camera_pos` — origin of all primary rays (ortho camera).
/// * `nx`, `ny` — image dimensions in pixels.
/// * Returns a flat `Vec` of `[r, g, b]` triples (row-major).
///
/// Each pixel fires one ray in the `-Z` direction; hit pixels get the
/// absolute-value of the triangle normal as colour, miss pixels return black.
pub fn gpu_render_pixels(
    triangles: &[GpuTriangle],
    camera_pos: [f32; 3],
    nx: usize,
    ny: usize,
) -> Vec<[f32; 3]> {
    let mut pixels = Vec::with_capacity(nx * ny);
    for row in 0..ny {
        for col in 0..nx {
            let u = (col as f32 + 0.5) / nx as f32 * 2.0 - 1.0;
            let v = (row as f32 + 0.5) / ny as f32 * 2.0 - 1.0;
            let ray = Ray::new(
                [camera_pos[0] + u, camera_pos[1] + v, camera_pos[2]],
                [0.0, 0.0, -1.0],
            );
            let colour = match trace_ray(&ray, triangles) {
                Some(hit) => [
                    hit.normal[0].abs(),
                    hit.normal[1].abs(),
                    hit.normal[2].abs(),
                ],
                None => [0.0, 0.0, 0.0],
            };
            pixels.push(colour);
        }
    }
    pixels
}

// ── Ambient occlusion ─────────────────────────────────────────────────────────

/// Estimate ambient occlusion at a hit point by sampling hemisphere directions.
///
/// `n_samples` rays are fired in directions cosine-weighted around `hit.normal`.
/// Returns the fraction of rays that are *not* occluded (1.0 = fully lit).
pub fn ambient_occlusion_sample(
    hit: &HitRecord,
    triangles: &[GpuTriangle],
    n_samples: usize,
) -> f32 {
    use rand::RngExt;
    if n_samples == 0 {
        return 1.0;
    }
    let mut rng = rand::rng();
    let mut unoccluded = 0usize;

    let n = normalize3f(hit.normal);
    // Build a local tangent frame
    let up = if n[0].abs() < 0.9 {
        [1.0_f32, 0.0, 0.0]
    } else {
        [0.0_f32, 1.0, 0.0]
    };
    let tangent = normalize3f(cross3f(n, up));
    let bitangent = cross3f(n, tangent);

    for _ in 0..n_samples {
        // Sample hemisphere using cosine-weighting
        let r1: f32 = rng.random_range(0.0_f32..1.0_f32);
        let r2: f32 = rng.random_range(0.0_f32..1.0_f32);
        let phi = 2.0 * std::f32::consts::PI * r1;
        let cos_theta = r2.sqrt();
        let sin_theta = (1.0_f32 - cos_theta * cos_theta).sqrt();
        let lx = sin_theta * phi.cos();
        let ly = sin_theta * phi.sin();
        let lz = cos_theta;
        let dir = [
            lx * tangent[0] + ly * bitangent[0] + lz * n[0],
            lx * tangent[1] + ly * bitangent[1] + lz * n[1],
            lx * tangent[2] + ly * bitangent[2] + lz * n[2],
        ];
        let ao_ray = Ray::new(hit.point, dir);
        if trace_ray(&ao_ray, triangles).is_none() {
            unoccluded += 1;
        }
    }
    unoccluded as f32 / n_samples as f32
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> GpuTriangle {
        GpuTriangle {
            v0: [0.0, 0.0, -1.0],
            v1: [1.0, 0.0, -1.0],
            v2: [0.0, 1.0, -1.0],
            normal: [0.0, 0.0, 1.0],
            material_id: 0,
        }
    }

    fn centered_ray() -> Ray {
        Ray::new([0.25, 0.25, 0.0], [0.0, 0.0, -1.0])
    }

    #[test]
    fn test_ray_at() {
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = r.at(2.0);
        assert!((p[0] - 2.0).abs() < 1e-6);
        assert!(p[1].abs() < 1e-6);
        assert!(p[2].abs() < 1e-6);
    }

    #[test]
    fn test_ray_sphere_hit() {
        let r = Ray::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        let t = ray_sphere_intersect(&r, [0.0, 0.0, 0.0], 1.0);
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-4);
    }

    #[test]
    fn test_ray_sphere_miss() {
        let r = Ray::new([5.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(ray_sphere_intersect(&r, [0.0, 0.0, 0.0], 1.0).is_none());
    }

    #[test]
    fn test_ray_sphere_inside() {
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let t = ray_sphere_intersect(&r, [0.0, 0.0, 0.0], 2.0);
        assert!(t.is_some());
    }

    #[test]
    fn test_ray_triangle_hit() {
        let tri = unit_triangle();
        let r = centered_ray();
        let t = ray_triangle_intersect(&r, &tri);
        assert!(t.is_some());
        assert!((t.unwrap() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_ray_triangle_miss_outside() {
        let tri = unit_triangle();
        let r = Ray::new([2.0, 2.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(ray_triangle_intersect(&r, &tri).is_none());
    }

    #[test]
    fn test_ray_triangle_parallel() {
        let tri = unit_triangle();
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(ray_triangle_intersect(&r, &tri).is_none());
    }

    #[test]
    fn test_ray_aabb_hit_direct() {
        let r = Ray::new([0.0, 0.0, 2.0], [0.0, 0.0, -1.0]);
        assert!(ray_aabb_intersect(&r, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]));
    }

    #[test]
    fn test_ray_aabb_miss() {
        let r = Ray::new([5.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(!ray_aabb_intersect(&r, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]));
    }

    #[test]
    fn test_ray_aabb_from_inside() {
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(ray_aabb_intersect(&r, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]));
    }

    #[test]
    fn test_trace_ray_hit() {
        let tri = unit_triangle();
        let r = centered_ray();
        let hit = trace_ray(&r, &[tri]);
        assert!(hit.is_some());
    }

    #[test]
    fn test_trace_ray_miss() {
        let tri = unit_triangle();
        let r = Ray::new([5.0, 5.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(trace_ray(&r, &[tri]).is_none());
    }

    #[test]
    fn test_trace_ray_closest() {
        // Two triangles at z=-1 and z=-2; ray should hit the nearer one
        let tri1 = unit_triangle(); // z=-1
        let tri2 = GpuTriangle {
            v0: [0.0, 0.0, -2.0],
            v1: [1.0, 0.0, -2.0],
            v2: [0.0, 1.0, -2.0],
            normal: [0.0, 0.0, 1.0],
            material_id: 1,
        };
        let r = centered_ray();
        let hit = trace_ray(&r, &[tri1, tri2]).unwrap();
        assert_eq!(hit.material_id, 0);
    }

    #[test]
    fn test_trace_ray_empty_scene() {
        let r = centered_ray();
        assert!(trace_ray(&r, &[]).is_none());
    }

    #[test]
    fn test_hit_record_point() {
        let tri = unit_triangle();
        let r = centered_ray();
        let hit = trace_ray(&r, &[tri]).unwrap();
        assert!((hit.point[2] - (-1.0)).abs() < 1e-4);
    }

    #[test]
    fn test_hit_record_normal() {
        let tri = unit_triangle();
        let r = centered_ray();
        let hit = trace_ray(&r, &[tri]).unwrap();
        assert!((hit.normal[2] - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_gpu_render_pixels_count() {
        let pixels = gpu_render_pixels(&[], [0.0, 0.0, 5.0], 4, 4);
        assert_eq!(pixels.len(), 16);
    }

    #[test]
    fn test_gpu_render_pixels_miss_black() {
        // Empty scene: all pixels should be black
        let pixels = gpu_render_pixels(&[], [0.0, 0.0, 5.0], 2, 2);
        for p in &pixels {
            assert_eq!(*p, [0.0, 0.0, 0.0]);
        }
    }

    #[test]
    fn test_gpu_render_pixels_hit_coloured() {
        // Single triangle centred at the viewport
        let tri = GpuTriangle {
            v0: [-2.0, -2.0, -1.0],
            v1: [2.0, -2.0, -1.0],
            v2: [0.0, 2.0, -1.0],
            normal: [0.0, 0.0, 1.0],
            material_id: 0,
        };
        let pixels = gpu_render_pixels(&[tri], [0.0, 0.0, 0.0], 3, 3);
        // Middle pixel should be non-black
        let has_hit = pixels.iter().any(|p| p[2] > 0.5);
        assert!(has_hit);
    }

    #[test]
    fn test_ambient_occlusion_zero_samples() {
        let hit = HitRecord {
            t: 1.0,
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            material_id: 0,
        };
        let ao = ambient_occlusion_sample(&hit, &[], 0);
        assert!((ao - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_ambient_occlusion_empty_scene() {
        let hit = HitRecord {
            t: 1.0,
            point: [0.0, 1.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            material_id: 0,
        };
        let ao = ambient_occlusion_sample(&hit, &[], 32);
        assert!((ao - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_ambient_occlusion_range() {
        let hit = HitRecord {
            t: 1.0,
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
            material_id: 0,
        };
        let ao = ambient_occlusion_sample(&hit, &[], 16);
        assert!((0.0..=1.0).contains(&ao));
    }

    #[test]
    fn test_bvh_node_fields() {
        let node = BvhNode {
            aabb_min: [-1.0, -1.0, -1.0],
            aabb_max: [1.0, 1.0, 1.0],
            left: 0,
            right: 1,
            is_leaf: true,
            tri_idx: 42,
        };
        assert_eq!(node.tri_idx, 42);
        assert!(node.is_leaf);
    }

    #[test]
    fn test_gpu_triangle_fields() {
        let tri = unit_triangle();
        assert_eq!(tri.material_id, 0);
        assert!((tri.normal[2] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_ray_sphere_tangent() {
        // Ray that just grazes the sphere
        let r = Ray::new([1.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        // sphere radius 1 at origin: tangent at x=1
        let t = ray_sphere_intersect(&r, [0.0, 0.0, 0.0], 1.0);
        assert!(t.is_some());
    }

    #[test]
    fn test_render_1x1_empty() {
        let pixels = gpu_render_pixels(&[], [0.0, 0.0, 1.0], 1, 1);
        assert_eq!(pixels.len(), 1);
        assert_eq!(pixels[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_ray_aabb_negative_direction() {
        let r = Ray::new([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!(ray_aabb_intersect(&r, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]));
    }
}
