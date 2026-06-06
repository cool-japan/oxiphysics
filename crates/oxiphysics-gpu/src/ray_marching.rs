// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ray marching and Signed Distance Field (SDF) rendering utilities.
//!
//! Provides a CPU-side implementation of common ray marching primitives:
//! SDFs for spheres, boxes, capsules, and tori; boolean CSG operators;
//! smooth blending; ambient occlusion; soft shadows; and a simple renderer
//! that produces depth and normal buffers.

// ── Vector helpers ────────────────────────────────────────────────────────────

/// Add two 3-D vectors component-wise.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors component-wise.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by a scalar.
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length of a 3-D vector.
#[inline]
fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalize a 3-D vector; returns zero vector when length < 1e-15.
#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = length3(v);
    if len < 1e-15 {
        return [0.0; 3];
    }
    scale3(v, 1.0 / len)
}

/// Clamp a value to \[lo, hi\].
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// Component-wise max of two 3-D vectors.
#[inline]
fn max3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])]
}

/// Component-wise absolute value of a 3-D vector.
#[cfg(test)]
#[inline]
fn abs3(v: [f64; 3]) -> [f64; 3] {
    [v[0].abs(), v[1].abs(), v[2].abs()]
}

// ── Ray ───────────────────────────────────────────────────────────────────────

/// A ray defined by an origin point and a (normalized) direction vector.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f64; 3],
    /// Ray direction (should be normalized).
    pub direction: [f64; 3],
}

impl Ray {
    /// Construct a new ray from origin and direction.
    ///
    /// The direction is normalized automatically.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction: normalize3(direction),
        }
    }

    /// Evaluate the ray at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}

// ── SDF trait ─────────────────────────────────────────────────────────────────

/// A Signed Distance Field: maps a 3-D point to the signed distance to the
/// nearest surface (negative inside, positive outside).
pub trait Sdf {
    /// Return the signed distance from point `p` to the nearest surface.
    fn distance(&self, p: [f64; 3]) -> f64;

    /// Return the surface normal at point `p` via central-difference gradient.
    ///
    /// The default implementation uses a small finite-difference epsilon of
    /// `1e-5`.
    fn normal(&self, p: [f64; 3]) -> [f64; 3] {
        let e = 1e-5;
        let dx = self.distance([p[0] + e, p[1], p[2]]) - self.distance([p[0] - e, p[1], p[2]]);
        let dy = self.distance([p[0], p[1] + e, p[2]]) - self.distance([p[0], p[1] - e, p[2]]);
        let dz = self.distance([p[0], p[1], p[2] + e]) - self.distance([p[0], p[1], p[2] - e]);
        normalize3([dx, dy, dz])
    }
}

// ── SphereSdf ─────────────────────────────────────────────────────────────────

/// SDF for a sphere.
#[derive(Debug, Clone, Copy)]
pub struct SphereSdf {
    /// Center of the sphere in world space.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}

impl SphereSdf {
    /// Construct a new `SphereSdf`.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
}

impl Sdf for SphereSdf {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let d = sub3(p, self.center);
        length3(d) - self.radius
    }
}

// ── BoxSdf ────────────────────────────────────────────────────────────────────

/// SDF for an axis-aligned box centred at the origin.
#[derive(Debug, Clone, Copy)]
pub struct BoxSdf {
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}

impl BoxSdf {
    /// Construct a new `BoxSdf`.
    pub fn new(half_extents: [f64; 3]) -> Self {
        Self { half_extents }
    }
}

impl Sdf for BoxSdf {
    fn distance(&self, p: [f64; 3]) -> f64 {
        // q = abs(p) - half_extents
        let q = [
            p[0].abs() - self.half_extents[0],
            p[1].abs() - self.half_extents[1],
            p[2].abs() - self.half_extents[2],
        ];
        let outside = length3(max3(q, [0.0; 3]));
        let inside = q[0].max(q[1]).max(q[2]).min(0.0);
        outside + inside
    }
}

// ── CapsuleSdf ────────────────────────────────────────────────────────────────

/// SDF for a capsule (cylinder with hemispherical caps).
#[derive(Debug, Clone, Copy)]
pub struct CapsuleSdf {
    /// One endpoint of the capsule axis.
    pub a: [f64; 3],
    /// Other endpoint of the capsule axis.
    pub b: [f64; 3],
    /// Radius of the capsule.
    pub r: f64,
}

impl CapsuleSdf {
    /// Construct a new `CapsuleSdf`.
    pub fn new(a: [f64; 3], b: [f64; 3], r: f64) -> Self {
        Self { a, b, r }
    }
}

impl Sdf for CapsuleSdf {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let pa = sub3(p, self.a);
        let ba = sub3(self.b, self.a);
        let h = clamp(dot3(pa, ba) / dot3(ba, ba), 0.0, 1.0);
        let closest = sub3(pa, scale3(ba, h));
        length3(closest) - self.r
    }
}

// ── TorusSdf ──────────────────────────────────────────────────────────────────

/// SDF for a torus lying in the XZ plane, centred at the origin.
#[derive(Debug, Clone, Copy)]
pub struct TorusSdf {
    /// Major radius (from origin to centre of tube).
    pub r_major: f64,
    /// Minor radius (radius of the tube).
    pub r_minor: f64,
}

impl TorusSdf {
    /// Construct a new `TorusSdf`.
    pub fn new(r_major: f64, r_minor: f64) -> Self {
        Self { r_major, r_minor }
    }
}

impl Sdf for TorusSdf {
    fn distance(&self, p: [f64; 3]) -> f64 {
        // Distance to ring in XZ plane
        let q_x = (p[0] * p[0] + p[2] * p[2]).sqrt() - self.r_major;
        let q_y = p[1];
        (q_x * q_x + q_y * q_y).sqrt() - self.r_minor
    }
}

// ── CSG operators ─────────────────────────────────────────────────────────────

/// SDF union: take the minimum of two distances.
///
/// The result is the SDF of the shape that is the union of the two objects.
pub fn sdf_union(a: f64, b: f64) -> f64 {
    a.min(b)
}

/// SDF intersection: take the maximum of two distances.
///
/// The result is the SDF of the shape that is the intersection of the two
/// objects.
pub fn sdf_intersection(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// SDF subtraction: subtract shape B from shape A.
///
/// Returns the SDF of `A \ B`.
pub fn sdf_subtraction(a: f64, b: f64) -> f64 {
    a.max(-b)
}

/// Polynomial smooth union of two SDF values with blend radius `k`.
///
/// When `k == 0` this degenerates to a hard union.  Larger `k` values produce
/// a smoother blend at the boundary between the two shapes.
pub fn sdf_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    if k < 1e-15 {
        return sdf_union(a, b);
    }
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    let mix = b + h * (a - b);
    mix - k * h * (1.0 - h)
}

// ── RayMarchResult ────────────────────────────────────────────────────────────

/// Result of a single ray march step.
#[derive(Debug, Clone, Copy)]
pub struct RayMarchResult {
    /// Whether the ray hit a surface.
    pub hit: bool,
    /// Ray parameter `t` at the hit point (or `max_dist` on miss).
    pub t: f64,
    /// Number of sphere-marching steps taken.
    pub steps: usize,
    /// Surface normal at the hit point (zero vector on miss).
    pub normal: [f64; 3],
}

impl RayMarchResult {
    /// Convenience constructor for a miss result.
    pub fn miss(t: f64, steps: usize) -> Self {
        Self {
            hit: false,
            t,
            steps,
            normal: [0.0; 3],
        }
    }

    /// Convenience constructor for a hit result.
    pub fn hit(t: f64, steps: usize, normal: [f64; 3]) -> Self {
        Self {
            hit: true,
            t,
            steps,
            normal,
        }
    }
}

// ── ray_march ─────────────────────────────────────────────────────────────────

/// Perform sphere-marching along `ray` against the provided SDF.
///
/// # Parameters
/// - `ray` – the ray to march.
/// - `sdf` – the scene SDF.
/// - `max_steps` – maximum iteration count (e.g. 128).
/// - `max_dist` – maximum travel distance before declaring a miss.
/// - `eps` – surface threshold: when `|d| < eps` a hit is declared.
pub fn ray_march(
    ray: &Ray,
    sdf: &dyn Sdf,
    max_steps: usize,
    max_dist: f64,
    eps: f64,
) -> RayMarchResult {
    let mut t = 0.0_f64;
    for step in 0..max_steps {
        let p = ray.at(t);
        let d = sdf.distance(p);
        if d.abs() < eps {
            let normal = sdf.normal(p);
            return RayMarchResult::hit(t, step + 1, normal);
        }
        t += d;
        if t >= max_dist {
            return RayMarchResult::miss(max_dist, step + 1);
        }
    }
    RayMarchResult::miss(t, max_steps)
}

// ── AmbientOcclusion ──────────────────────────────────────────────────────────

/// Ambient occlusion estimator using SDF step-based sampling.
#[derive(Debug, Clone, Copy)]
pub struct AmbientOcclusion {
    /// Number of sampling steps along the normal.
    pub samples: usize,
    /// Maximum AO sampling radius.
    pub radius: f64,
    /// Falloff exponent (higher = quicker falloff).
    pub falloff: f64,
}

impl AmbientOcclusion {
    /// Construct a new `AmbientOcclusion` estimator.
    pub fn new(samples: usize, radius: f64, falloff: f64) -> Self {
        Self {
            samples,
            radius,
            falloff,
        }
    }

    /// Estimate ambient occlusion at position `pos` with surface normal
    /// `normal` against the given SDF.
    ///
    /// Returns a value in `[0, 1]` where 1.0 means fully unoccluded.
    pub fn compute(&self, pos: [f64; 3], normal: [f64; 3], sdf: &dyn Sdf) -> f64 {
        if self.samples == 0 {
            return 1.0;
        }
        let mut occ = 0.0_f64;
        for i in 1..=self.samples {
            let t = self.radius * (i as f64) / (self.samples as f64);
            let sample_pos = add3(pos, scale3(normal, t));
            let d = sdf.distance(sample_pos);
            occ += (t - d).max(0.0) / t.powf(self.falloff);
        }
        (1.0 - occ / (self.samples as f64)).clamp(0.0, 1.0)
    }
}

// ── SoftShadow ────────────────────────────────────────────────────────────────

/// Penumbra-aware soft shadow estimator.
#[derive(Debug, Clone, Copy)]
pub struct SoftShadow {
    /// Direction toward the light source (will be normalized).
    pub light_dir: [f64; 3],
    /// Sharpness of the shadow penumbra; higher `k` → harder shadow.
    pub k: f64,
}

impl SoftShadow {
    /// Construct a new `SoftShadow`.
    pub fn new(light_dir: [f64; 3], k: f64) -> Self {
        Self {
            light_dir: normalize3(light_dir),
            k,
        }
    }

    /// Compute soft shadow factor at surface point `pos`.
    ///
    /// Returns a value in `[0, 1]`: 0.0 = fully in shadow, 1.0 = fully lit.
    pub fn compute(&self, pos: [f64; 3], sdf: &dyn Sdf, max_dist: f64) -> f64 {
        let eps = 1e-4;
        let mut res = 1.0_f64;
        let mut t = eps;
        while t < max_dist {
            let p = add3(pos, scale3(self.light_dir, t));
            let d = sdf.distance(p);
            if d < eps {
                return 0.0;
            }
            res = res.min(self.k * d / t);
            t += d;
        }
        res.clamp(0.0, 1.0)
    }
}

// ── RayMarchRenderer ──────────────────────────────────────────────────────────

/// A simple camera/renderer that generates rays and produces image buffers via
/// ray marching.
#[derive(Debug, Clone)]
pub struct RayMarchRenderer {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Vertical field of view in radians.
    pub fov: f64,
    /// Camera position in world space.
    pub camera_pos: [f64; 3],
    /// Point the camera looks at.
    pub camera_target: [f64; 3],
}

impl RayMarchRenderer {
    /// Construct a new renderer.
    pub fn new(
        width: usize,
        height: usize,
        fov: f64,
        camera_pos: [f64; 3],
        camera_target: [f64; 3],
    ) -> Self {
        Self {
            width,
            height,
            fov,
            camera_pos,
            camera_target,
        }
    }

    /// Build the camera coordinate frame (right, up, forward).
    fn camera_basis(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let world_up = [0.0_f64, 1.0, 0.0];
        let forward = normalize3(sub3(self.camera_target, self.camera_pos));
        let right = normalize3(cross3_fn(forward, world_up));
        let up = cross3_fn(right, forward);
        (right, up, forward)
    }

    /// Generate the primary ray for pixel `(px, py)`.
    ///
    /// Pixel coordinates are measured from the top-left corner.
    pub fn generate_ray(&self, px: usize, py: usize) -> Ray {
        let (right, up, forward) = self.camera_basis();
        let aspect = self.width as f64 / self.height.max(1) as f64;
        let half_h = (self.fov * 0.5).tan();
        let half_w = half_h * aspect;

        // NDC → [-1, 1], y flipped
        let u = (2.0 * (px as f64 + 0.5) / self.width as f64 - 1.0) * half_w;
        let v = (1.0 - 2.0 * (py as f64 + 0.5) / self.height as f64) * half_h;

        let dir = add3(add3(forward, scale3(right, u)), scale3(up, v));
        Ray::new(self.camera_pos, dir)
    }

    /// Render the scene as a depth buffer.
    ///
    /// Returns a flat `width × height` `Vec`f64` where each element is the
    /// `t` value at the nearest surface hit, or `f64::INFINITY` on a miss.
    pub fn render_depth(&self, sdf: &dyn Sdf) -> Vec<f64> {
        let n = self.width * self.height;
        let mut depth = Vec::with_capacity(n);
        for py in 0..self.height {
            for px in 0..self.width {
                let ray = self.generate_ray(px, py);
                let result = ray_march(&ray, sdf, 256, 1000.0, 1e-5);
                depth.push(if result.hit { result.t } else { f64::INFINITY });
            }
        }
        depth
    }

    /// Render the scene as a normal buffer.
    ///
    /// Returns a flat `width × height` `Vec<\[f64; 3\]>` where each element is
    /// the surface normal at the hit point, or `\[0,0,0\]` on a miss.
    pub fn render_normals(&self, sdf: &dyn Sdf) -> Vec<[f64; 3]> {
        let n = self.width * self.height;
        let mut normals = Vec::with_capacity(n);
        for py in 0..self.height {
            for px in 0..self.width {
                let ray = self.generate_ray(px, py);
                let result = ray_march(&ray, sdf, 256, 1000.0, 1e-5);
                normals.push(result.normal);
            }
        }
        normals
    }
}

/// Internal cross product helper (avoids collision with public function name).
#[inline]
fn cross3_fn(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ── Subsurface scattering ─────────────────────────────────────────────────────

/// Approximate subsurface scattering transmittance.
///
/// Models the fraction of light that penetrates to `depth` below the surface,
/// using a simple exponential Beer-Lambert model modified by the absorption
/// coefficient.
///
/// # Parameters
/// - `depth` – penetration depth (metres, or consistent units).
/// - `scatter_coeff` – scattering coefficient (1/m).
/// - `absorption` – absorption coefficient (1/m).
///
/// Returns a value in `\[0, 1\]`.
pub fn subsurface_scattering_approx(depth: f64, scatter_coeff: f64, absorption: f64) -> f64 {
    let extinction = scatter_coeff + absorption;
    (-extinction * depth.max(0.0)).exp()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── Ray tests ────────────────────────────────────────────────────────

    #[test]
    fn ray_at_origin() {
        let ray = Ray::new([0.0; 3], [0.0, 0.0, 1.0]);
        let p = ray.at(0.0);
        assert_eq!(p, [0.0; 3]);
    }

    #[test]
    fn ray_at_t_positive() {
        let ray = Ray::new([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let p = ray.at(5.0);
        assert!((p[0] - 1.0).abs() < 1e-12);
        assert!((p[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn ray_direction_normalized() {
        let ray = Ray::new([0.0; 3], [3.0, 4.0, 0.0]);
        let len = length3(ray.direction);
        assert!((len - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ray_at_negative_t() {
        let ray = Ray::new([0.0; 3], [1.0, 0.0, 0.0]);
        let p = ray.at(-2.0);
        assert!((p[0] + 2.0).abs() < 1e-12);
    }

    // ── SphereSdf tests ──────────────────────────────────────────────────

    #[test]
    fn sphere_distance_outside() {
        let s = SphereSdf::new([0.0; 3], 1.0);
        let d = s.distance([2.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_distance_inside() {
        let s = SphereSdf::new([0.0; 3], 1.0);
        let d = s.distance([0.0; 3]);
        assert!((d + 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_distance_surface() {
        let s = SphereSdf::new([0.0; 3], 1.0);
        let d = s.distance([1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn sphere_normal_at_x_axis() {
        let s = SphereSdf::new([0.0; 3], 1.0);
        let n = s.normal([1.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-4);
        assert!(n[1].abs() < 1e-4);
        assert!(n[2].abs() < 1e-4);
    }

    #[test]
    fn sphere_translated_center() {
        let s = SphereSdf::new([5.0, 0.0, 0.0], 2.0);
        let d = s.distance([7.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);
    }

    // ── BoxSdf tests ─────────────────────────────────────────────────────

    #[test]
    fn box_distance_outside_along_x() {
        let b = BoxSdf::new([1.0; 3]);
        let d = b.distance([3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-12);
    }

    #[test]
    fn box_distance_inside() {
        let b = BoxSdf::new([1.0; 3]);
        let d = b.distance([0.0; 3]);
        assert!(d < 0.0);
    }

    #[test]
    fn box_distance_on_face() {
        let b = BoxSdf::new([1.0, 1.0, 1.0]);
        let d = b.distance([1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn box_normal_x_face() {
        let b = BoxSdf::new([1.0; 3]);
        let n = b.normal([1.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-3);
    }

    // ── CapsuleSdf tests ─────────────────────────────────────────────────

    #[test]
    fn capsule_distance_at_midpoint() {
        let c = CapsuleSdf::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let d = c.distance([0.5, 0.0, 0.0]);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn capsule_distance_outside() {
        let c = CapsuleSdf::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let d = c.distance([2.0, 0.0, 0.0]);
        assert!((d - 1.5).abs() < 1e-10);
    }

    #[test]
    fn capsule_distance_at_cap() {
        let c = CapsuleSdf::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        // Point directly above top cap
        let d = c.distance([0.0, 2.0, 0.0]);
        assert!((d - 0.5).abs() < 1e-10);
    }

    // ── TorusSdf tests ───────────────────────────────────────────────────

    #[test]
    fn torus_distance_on_ring() {
        // Point on the major ring at (r_major, 0, 0) displaced by r_minor in y
        let t = TorusSdf::new(2.0, 0.5);
        let d = t.distance([2.0, 0.5, 0.0]);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn torus_distance_centre_negative() {
        let t = TorusSdf::new(2.0, 0.5);
        // Origin is inside if r_major < r_minor, otherwise outside
        let d = t.distance([0.0; 3]);
        // distance = sqrt((0 - r_major)^2 + 0^2) - r_minor = r_major - r_minor
        assert!((d - 1.5).abs() < 1e-10);
    }

    #[test]
    fn torus_distance_outside_far() {
        let t = TorusSdf::new(1.0, 0.2);
        let d = t.distance([10.0, 0.0, 0.0]);
        assert!(d > 0.0);
    }

    // ── CSG operator tests ───────────────────────────────────────────────

    #[test]
    fn sdf_union_picks_min() {
        assert_eq!(sdf_union(3.0, 1.0), 1.0);
        assert_eq!(sdf_union(-1.0, 2.0), -1.0);
    }

    #[test]
    fn sdf_intersection_picks_max() {
        assert_eq!(sdf_intersection(3.0, 1.0), 3.0);
    }

    #[test]
    fn sdf_subtraction_correctness() {
        // subtract b from a: max(a, -b)
        assert_eq!(sdf_subtraction(1.0, -2.0), 2.0);
        assert_eq!(sdf_subtraction(1.0, 3.0), 1.0);
    }

    #[test]
    fn sdf_smooth_union_degenerate() {
        // k=0 → hard union
        let su = sdf_smooth_union(3.0, 1.0, 0.0);
        assert!((su - sdf_union(3.0, 1.0)).abs() < 1e-10);
    }

    #[test]
    fn sdf_smooth_union_blends() {
        let su = sdf_smooth_union(0.0, 0.0, 1.0);
        // At equal distances with large k, should be <= 0
        assert!(su <= 0.0 + 1e-10);
    }

    #[test]
    fn sdf_smooth_union_between_values() {
        // result should be <= min(a, b)
        let a = 2.0_f64;
        let b = 3.0_f64;
        let su = sdf_smooth_union(a, b, 0.5);
        assert!(su <= a + 1e-10);
    }

    // ── RayMarchResult tests ─────────────────────────────────────────────

    #[test]
    fn ray_march_result_miss() {
        let r = RayMarchResult::miss(100.0, 10);
        assert!(!r.hit);
        assert!((r.t - 100.0).abs() < 1e-12);
    }

    #[test]
    fn ray_march_result_hit() {
        let r = RayMarchResult::hit(5.0, 20, [1.0, 0.0, 0.0]);
        assert!(r.hit);
        assert!((r.normal[0] - 1.0).abs() < 1e-12);
    }

    // ── ray_march function tests ─────────────────────────────────────────

    #[test]
    fn ray_march_hits_sphere() {
        let sphere = SphereSdf::new([0.0, 0.0, 5.0], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let result = ray_march(&ray, &sphere, 256, 100.0, 1e-5);
        assert!(result.hit);
        assert!((result.t - 4.0).abs() < 1e-3);
    }

    #[test]
    fn ray_march_misses_sphere() {
        let sphere = SphereSdf::new([0.0, 0.0, 5.0], 1.0);
        let ray = Ray::new([10.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let result = ray_march(&ray, &sphere, 256, 100.0, 1e-5);
        assert!(!result.hit);
    }

    #[test]
    fn ray_march_steps_bounded() {
        let sphere = SphereSdf::new([0.0, 0.0, 5.0], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let result = ray_march(&ray, &sphere, 128, 100.0, 1e-5);
        assert!(result.steps <= 128);
    }

    #[test]
    fn ray_march_hit_normal_approx_minus_z() {
        let sphere = SphereSdf::new([0.0, 0.0, 5.0], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let result = ray_march(&ray, &sphere, 256, 100.0, 1e-5);
        // Normal at front face should point toward camera (-z direction)
        assert!(result.hit);
        assert!(result.normal[2] < 0.0);
    }

    // ── AmbientOcclusion tests ───────────────────────────────────────────

    #[test]
    fn ao_open_space_returns_one() {
        // A point far from any surface should be fully unoccluded
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let ao = AmbientOcclusion::new(8, 1.0, 1.0);
        let pos = [100.0, 0.0, 0.0];
        let normal = [1.0, 0.0, 0.0];
        let v = ao.compute(pos, normal, &sphere);
        assert!(v > 0.5);
    }

    #[test]
    fn ao_samples_zero_returns_one() {
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let ao = AmbientOcclusion::new(0, 1.0, 1.0);
        let v = ao.compute([0.0; 3], [0.0, 1.0, 0.0], &sphere);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ao_result_clamped() {
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let ao = AmbientOcclusion::new(4, 0.5, 1.0);
        let pos = [1.0, 0.0, 0.0];
        let normal = [1.0, 0.0, 0.0];
        let v = ao.compute(pos, normal, &sphere);
        assert!((0.0..=1.0).contains(&v));
    }

    // ── SoftShadow tests ─────────────────────────────────────────────────

    #[test]
    fn soft_shadow_lit_no_occluder() {
        // Point with no occluder between it and the light direction
        let sphere = SphereSdf::new([0.0, 0.0, -100.0], 0.1);
        let ss = SoftShadow::new([0.0, 1.0, 0.0], 8.0);
        let pos = [0.0, 0.0, 0.0];
        let v = ss.compute(pos, &sphere, 50.0);
        assert!(v > 0.9);
    }

    #[test]
    fn soft_shadow_shadow_with_occluder() {
        // Position the sphere directly in the shadow ray path
        let sphere = SphereSdf::new([0.0, 5.0, 0.0], 1.0);
        let ss = SoftShadow::new([0.0, 1.0, 0.0], 8.0);
        let pos = [0.0, 0.0, 0.0];
        let v = ss.compute(pos, &sphere, 20.0);
        assert!(v < 0.5);
    }

    #[test]
    fn soft_shadow_result_in_range() {
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let ss = SoftShadow::new([1.0, 1.0, 1.0], 4.0);
        let v = ss.compute([3.0, 0.0, 0.0], &sphere, 20.0);
        assert!((0.0..=1.0).contains(&v));
    }

    // ── RayMarchRenderer tests ───────────────────────────────────────────

    #[test]
    fn renderer_depth_buffer_size() {
        let renderer = RayMarchRenderer::new(4, 4, PI / 3.0, [0.0, 0.0, -5.0], [0.0, 0.0, 0.0]);
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let depth = renderer.render_depth(&sphere);
        assert_eq!(depth.len(), 16);
    }

    #[test]
    fn renderer_normals_buffer_size() {
        let renderer = RayMarchRenderer::new(2, 2, PI / 3.0, [0.0, 0.0, -5.0], [0.0, 0.0, 0.0]);
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let normals = renderer.render_normals(&sphere);
        assert_eq!(normals.len(), 4);
    }

    #[test]
    fn renderer_center_ray_hits_sphere() {
        let renderer = RayMarchRenderer::new(11, 11, PI / 3.0, [0.0, 0.0, -5.0], [0.0, 0.0, 0.0]);
        let sphere = SphereSdf::new([0.0; 3], 1.0);
        let depth = renderer.render_depth(&sphere);
        // Center pixel should hit the sphere
        let center_t = depth[5 * 11 + 5];
        assert!(center_t < f64::INFINITY);
    }

    #[test]
    fn renderer_corner_ray_misses_small_sphere() {
        let renderer = RayMarchRenderer::new(11, 11, PI / 3.0, [0.0, 0.0, -5.0], [0.0, 0.0, 0.0]);
        let sphere = SphereSdf::new([0.0; 3], 0.01);
        let depth = renderer.render_depth(&sphere);
        // Corner pixel (0,0) should miss the tiny sphere
        assert_eq!(depth[0], f64::INFINITY);
    }

    #[test]
    fn renderer_generate_ray_center() {
        let renderer = RayMarchRenderer::new(11, 11, PI / 3.0, [0.0, 0.0, -5.0], [0.0, 0.0, 0.0]);
        let ray = renderer.generate_ray(5, 5);
        // Center ray should point roughly along +z
        assert!(ray.direction[2] > 0.0);
    }

    // ── subsurface_scattering_approx tests ───────────────────────────────

    #[test]
    fn sss_zero_depth_returns_one() {
        let v = subsurface_scattering_approx(0.0, 1.0, 1.0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sss_large_depth_near_zero() {
        let v = subsurface_scattering_approx(1000.0, 1.0, 0.0);
        assert!(v < 1e-10);
    }

    #[test]
    fn sss_negative_depth_clamped() {
        let v = subsurface_scattering_approx(-5.0, 1.0, 1.0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sss_result_in_range() {
        let v = subsurface_scattering_approx(0.5, 2.0, 1.0);
        assert!((0.0..=1.0).contains(&v));
    }

    #[test]
    fn sss_monotone_decreasing() {
        let v1 = subsurface_scattering_approx(1.0, 1.0, 0.5);
        let v2 = subsurface_scattering_approx(2.0, 1.0, 0.5);
        assert!(v1 > v2);
    }

    // ── Vector helper internal tests ─────────────────────────────────────

    #[test]
    fn test_add3() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_sub3() {
        let r = sub3([4.0, 5.0, 6.0], [1.0, 2.0, 3.0]);
        assert_eq!(r, [3.0, 3.0, 3.0]);
    }

    #[test]
    fn test_scale3() {
        let r = scale3([1.0, 2.0, 3.0], 2.0);
        assert_eq!(r, [2.0, 4.0, 6.0]);
    }

    #[test]
    fn test_abs3() {
        let r = abs3([-1.0, 2.0, -3.0]);
        assert_eq!(r, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_max3() {
        let r = max3([1.0, -1.0, 2.0], [0.0, 0.0, 0.0]);
        assert_eq!(r, [1.0, 0.0, 2.0]);
    }

    #[test]
    fn test_normalize3_zero() {
        let n = normalize3([0.0; 3]);
        assert_eq!(n, [0.0; 3]);
    }

    #[test]
    fn test_clamp() {
        assert_eq!(clamp(-1.0, 0.0, 1.0), 0.0);
        assert_eq!(clamp(0.5, 0.0, 1.0), 0.5);
        assert_eq!(clamp(2.0, 0.0, 1.0), 1.0);
    }
}
