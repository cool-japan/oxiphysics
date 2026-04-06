// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU path tracer using GPU-style algorithms.
//!
//! Implements a Monte Carlo path tracer with Lambertian, Metal, and Dielectric
//! materials, sphere and triangle primitives, and a progressive pixel buffer.

#![allow(dead_code)]

use rand::Rng;

use rand::RngExt;
// ── Vector helpers (f32, 3D) ─────────────────────────────────────────────────

#[inline]
fn vadd(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vsub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vmul(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vmul3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

#[inline]
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn length(v: [f32; 3]) -> f32 {
    dot(v, v).sqrt()
}

#[inline]
fn normalize(v: [f32; 3]) -> [f32; 3] {
    let l = length(v);
    if l < 1e-8 {
        return [0.0; 3];
    }
    vmul(v, 1.0 / l)
}

#[inline]
fn reflect(d: [f32; 3], n: [f32; 3]) -> [f32; 3] {
    vsub(d, vmul(n, 2.0 * dot(d, n)))
}

fn refract(uv: [f32; 3], n: [f32; 3], ni_over_nt: f32) -> Option<[f32; 3]> {
    let cos_theta = (-dot(uv, n)).min(1.0);
    let r_out_perp = vmul(vadd(uv, vmul(n, cos_theta)), ni_over_nt);
    let r_out_parallel_len2 = (1.0 - dot(r_out_perp, r_out_perp)).abs();
    let r_out_parallel = vmul(n, -(r_out_parallel_len2.sqrt()));
    Some(vadd(r_out_perp, r_out_parallel))
}

fn schlick(cosine: f32, ref_idx: f32) -> f32 {
    let r0 = ((1.0 - ref_idx) / (1.0 + ref_idx)).powi(2);
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}

fn random_in_unit_sphere(rng: &mut impl Rng) -> [f32; 3] {
    loop {
        let v = [
            rng.random_range(-1.0f32..1.0),
            rng.random_range(-1.0f32..1.0),
            rng.random_range(-1.0f32..1.0),
        ];
        if dot(v, v) < 1.0 {
            return v;
        }
    }
}

fn random_unit_vector(rng: &mut impl Rng) -> [f32; 3] {
    normalize(random_in_unit_sphere(rng))
}

// ── Ray ──────────────────────────────────────────────────────────────────────

/// A ray with an origin and a direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f32; 3],
    /// Ray direction (need not be normalised, but usually is).
    pub direction: [f32; 3],
}

impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f32; 3], direction: [f32; 3]) -> Self {
        Self { origin, direction }
    }

    /// Evaluate the ray at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: f32) -> [f32; 3] {
        vadd(self.origin, vmul(self.direction, t))
    }
}

// ── Material ─────────────────────────────────────────────────────────────────

/// The scattering type of a material.
#[derive(Debug, Clone, Copy)]
pub enum MaterialType {
    /// Perfectly diffuse (Lambertian) scattering.
    Lambertian,
    /// Metallic reflection with a roughness fuzz factor in \[0,1\].
    Metal(f32),
    /// Dielectric (glass-like) material with index of refraction.
    Dielectric(f32),
}

/// A surface material with albedo colour and scattering type.
#[derive(Debug, Clone, Copy)]
pub struct Material {
    /// Base colour of the material (RGB in \[0,1\]).
    pub albedo: [f32; 3],
    /// Scattering behaviour.
    pub kind: MaterialType,
}

impl Material {
    /// Create a Lambertian (diffuse) material.
    pub fn lambertian(albedo: [f32; 3]) -> Self {
        Self {
            albedo,
            kind: MaterialType::Lambertian,
        }
    }

    /// Create a metallic material.
    pub fn metal(albedo: [f32; 3], fuzz: f32) -> Self {
        Self {
            albedo,
            kind: MaterialType::Metal(fuzz.clamp(0.0, 1.0)),
        }
    }

    /// Create a dielectric (glass) material.
    pub fn dielectric(ior: f32) -> Self {
        Self {
            albedo: [1.0; 3],
            kind: MaterialType::Dielectric(ior),
        }
    }

    /// Scatter a ray off this material surface.
    ///
    /// Returns `Some((scattered_ray, attenuation))` or `None` if absorbed.
    pub fn scatter(
        &self,
        ray: &Ray,
        hit: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<(Ray, [f32; 3])> {
        match self.kind {
            MaterialType::Lambertian => {
                let target = vadd(vadd(hit.point, hit.normal), random_unit_vector(rng));
                let scattered = Ray::new(hit.point, vsub(target, hit.point));
                Some((scattered, self.albedo))
            }
            MaterialType::Metal(fuzz) => {
                let reflected = reflect(normalize(ray.direction), hit.normal);
                let fuzzed = vadd(reflected, vmul(random_in_unit_sphere(rng), fuzz));
                if dot(fuzzed, hit.normal) > 0.0 {
                    Some((Ray::new(hit.point, fuzzed), self.albedo))
                } else {
                    None
                }
            }
            MaterialType::Dielectric(ior) => {
                let attenuation = [1.0f32; 3];
                let refraction_ratio = if hit.front_face { 1.0 / ior } else { ior };
                let unit_dir = normalize(ray.direction);
                let cos_theta = (-dot(unit_dir, hit.normal)).min(1.0);
                let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
                let cannot_refract = refraction_ratio * sin_theta > 1.0;
                let scattered_dir = if cannot_refract
                    || schlick(cos_theta, refraction_ratio) > rng.random::<f32>()
                {
                    reflect(unit_dir, hit.normal)
                } else {
                    refract(unit_dir, hit.normal, refraction_ratio)
                        .unwrap_or_else(|| reflect(unit_dir, hit.normal))
                };
                Some((Ray::new(hit.point, scattered_dir), attenuation))
            }
        }
    }
}

// ── HitRecord ────────────────────────────────────────────────────────────────

/// Record of a ray–surface intersection.
#[derive(Debug, Clone, Copy)]
pub struct HitRecord {
    /// Ray parameter at intersection.
    pub t: f32,
    /// World-space hit point.
    pub point: [f32; 3],
    /// Outward-facing surface normal (normalised).
    pub normal: [f32; 3],
    /// Index into the scene's material list.
    pub material_index: usize,
    /// True when the ray hits the front face.
    pub front_face: bool,
}

impl HitRecord {
    fn new(t: f32, point: [f32; 3], outward_normal: [f32; 3], ray: &Ray, mat: usize) -> Self {
        let front_face = dot(ray.direction, outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            vmul(outward_normal, -1.0)
        };
        Self {
            t,
            point,
            normal,
            material_index: mat,
            front_face,
        }
    }
}

// ── Sphere ───────────────────────────────────────────────────────────────────

/// A sphere primitive.
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    /// Centre of the sphere.
    pub center: [f32; 3],
    /// Radius of the sphere.
    pub radius: f32,
    /// Index into the scene's material list.
    pub material_index: usize,
}

impl Sphere {
    /// Create a new sphere.
    pub fn new(center: [f32; 3], radius: f32, material_index: usize) -> Self {
        Self {
            center,
            radius,
            material_index,
        }
    }

    /// Test ray–sphere intersection in the interval `(t_min, t_max)`.
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let oc = vsub(ray.origin, self.center);
        let a = dot(ray.direction, ray.direction);
        let half_b = dot(oc, ray.direction);
        let c = dot(oc, oc) - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }
        let sqrt_d = discriminant.sqrt();
        let mut root = (-half_b - sqrt_d) / a;
        if root < t_min || root > t_max {
            root = (-half_b + sqrt_d) / a;
            if root < t_min || root > t_max {
                return None;
            }
        }
        let point = ray.at(root);
        let outward_normal = vmul(vsub(point, self.center), 1.0 / self.radius);
        Some(HitRecord::new(
            root,
            point,
            outward_normal,
            ray,
            self.material_index,
        ))
    }
}

// ── Triangle ─────────────────────────────────────────────────────────────────

/// A triangle primitive.
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// Vertex A.
    pub v0: [f32; 3],
    /// Vertex B.
    pub v1: [f32; 3],
    /// Vertex C.
    pub v2: [f32; 3],
    /// Precomputed geometric normal (normalised).
    pub normal: [f32; 3],
    /// Index into the scene's material list.
    pub material_index: usize,
}

impl Triangle {
    /// Create a new triangle, computing the normal automatically.
    pub fn new(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3], material_index: usize) -> Self {
        let edge1 = vsub(v1, v0);
        let edge2 = vsub(v2, v0);
        let normal = normalize(cross(edge1, edge2));
        Self {
            v0,
            v1,
            v2,
            normal,
            material_index,
        }
    }

    /// Möller–Trumbore ray–triangle intersection.
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        const EPSILON: f32 = 1e-7;
        let edge1 = vsub(self.v1, self.v0);
        let edge2 = vsub(self.v2, self.v0);
        let h = cross(ray.direction, edge2);
        let a = dot(edge1, h);
        if a.abs() < EPSILON {
            return None; // parallel
        }
        let f = 1.0 / a;
        let s = vsub(ray.origin, self.v0);
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
        if t < t_min || t > t_max {
            return None;
        }
        let point = ray.at(t);
        Some(HitRecord::new(
            t,
            point,
            self.normal,
            ray,
            self.material_index,
        ))
    }
}

// ── Light ────────────────────────────────────────────────────────────────────

/// A point light source.
#[derive(Debug, Clone, Copy)]
pub struct PointLight {
    /// Light position.
    pub position: [f32; 3],
    /// Light colour/intensity (RGB).
    pub color: [f32; 3],
    /// Light intensity multiplier.
    pub intensity: f32,
}

impl PointLight {
    /// Create a new point light.
    pub fn new(position: [f32; 3], color: [f32; 3], intensity: f32) -> Self {
        Self {
            position,
            color,
            intensity,
        }
    }
}

// ── PathTracerScene ──────────────────────────────────────────────────────────

/// A scene containing geometry, materials, and lights.
#[derive(Debug, Clone, Default)]
pub struct PathTracerScene {
    /// Sphere primitives in the scene.
    pub spheres: Vec<Sphere>,
    /// Triangle primitives in the scene.
    pub triangles: Vec<Triangle>,
    /// Point lights in the scene.
    pub lights: Vec<PointLight>,
    /// Material list shared by all primitives.
    pub materials: Vec<Material>,
    /// Background (sky) gradient top colour.
    pub sky_top: [f32; 3],
    /// Background (sky) gradient bottom colour.
    pub sky_bottom: [f32; 3],
}

impl PathTracerScene {
    /// Create a new empty scene.
    pub fn new() -> Self {
        Self {
            sky_top: [0.5, 0.7, 1.0],
            sky_bottom: [1.0, 1.0, 1.0],
            ..Default::default()
        }
    }

    /// Add a material and return its index.
    pub fn add_material(&mut self, mat: Material) -> usize {
        let idx = self.materials.len();
        self.materials.push(mat);
        idx
    }

    /// Add a sphere to the scene.
    pub fn add_sphere(&mut self, sphere: Sphere) {
        self.spheres.push(sphere);
    }

    /// Add a triangle to the scene.
    pub fn add_triangle(&mut self, triangle: Triangle) {
        self.triangles.push(triangle);
    }

    /// Add a point light to the scene.
    pub fn add_light(&mut self, light: PointLight) {
        self.lights.push(light);
    }

    /// Find the closest intersection along a ray.
    pub fn hit_scene(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut closest: Option<HitRecord> = None;
        let mut t_closest = t_max;
        for sphere in &self.spheres {
            if let Some(rec) = sphere.hit(ray, t_min, t_closest) {
                t_closest = rec.t;
                closest = Some(rec);
            }
        }
        for tri in &self.triangles {
            if let Some(rec) = tri.hit(ray, t_min, t_closest) {
                t_closest = rec.t;
                closest = Some(rec);
            }
        }
        closest
    }

    /// Sky background colour for a given ray direction.
    fn sky_color(&self, ray: &Ray) -> [f32; 3] {
        let unit = normalize(ray.direction);
        let t = 0.5 * (unit[1] + 1.0);
        let a = self.sky_bottom;
        let b = self.sky_top;
        [
            a[0] * (1.0 - t) + b[0] * t,
            a[1] * (1.0 - t) + b[1] * t,
            a[2] * (1.0 - t) + b[2] * t,
        ]
    }

    /// Trace a ray through the scene with Monte Carlo path tracing.
    ///
    /// Returns the estimated radiance (RGB) for the ray.
    pub fn trace(&self, ray: &Ray, max_depth: usize, rng: &mut impl Rng) -> [f32; 3] {
        if max_depth == 0 {
            return [0.0; 3];
        }
        if let Some(hit) = self.hit_scene(ray, 1e-4, f32::INFINITY) {
            let mat = &self.materials[hit.material_index];
            if let Some((scattered, attenuation)) = mat.scatter(ray, &hit, rng) {
                let incoming = self.trace(&scattered, max_depth - 1, rng);
                vmul3(attenuation, incoming)
            } else {
                [0.0; 3]
            }
        } else {
            self.sky_color(ray)
        }
    }
}

// ── PathTracerBuffer ─────────────────────────────────────────────────────────

/// A pixel accumulation buffer for progressive rendering.
#[derive(Debug, Clone)]
pub struct PathTracerBuffer {
    /// Buffer width in pixels.
    pub width: usize,
    /// Buffer height in pixels.
    pub height: usize,
    /// Accumulated colour per pixel (RGB, floating point).
    pub accumulator: Vec<[f32; 3]>,
    /// Number of samples accumulated per pixel.
    pub sample_count: Vec<u32>,
}

impl PathTracerBuffer {
    /// Create a new buffer of `width × height` pixels, all zeroed.
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            accumulator: vec![[0.0; 3]; n],
            sample_count: vec![0; n],
        }
    }

    /// Add a colour sample to pixel `(x, y)`.
    pub fn add_sample(&mut self, x: usize, y: usize, color: [f32; 3]) {
        let idx = y * self.width + x;
        let acc = &mut self.accumulator[idx];
        acc[0] += color[0];
        acc[1] += color[1];
        acc[2] += color[2];
        self.sample_count[idx] += 1;
    }

    /// Get the averaged colour at pixel `(x, y)`.
    pub fn get_pixel(&self, x: usize, y: usize) -> [f32; 3] {
        let idx = y * self.width + x;
        let n = self.sample_count[idx] as f32;
        if n == 0.0 {
            return [0.0; 3];
        }
        let acc = self.accumulator[idx];
        [acc[0] / n, acc[1] / n, acc[2] / n]
    }

    /// Return a gamma-corrected (gamma=2) u8 RGB image row-major.
    pub fn to_rgb8(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.width * self.height * 3);
        for y in 0..self.height {
            for x in 0..self.width {
                let c = self.get_pixel(x, y);
                for ch in c.iter() {
                    let linear = ch.clamp(0.0, 1.0);
                    let gamma = linear.sqrt(); // gamma = 2
                    out.push((gamma * 255.999) as u8);
                }
            }
        }
        out
    }

    /// Reset all accumulated samples.
    pub fn clear(&mut self) {
        for acc in &mut self.accumulator {
            *acc = [0.0; 3];
        }
        for s in &mut self.sample_count {
            *s = 0;
        }
    }

    /// Total number of samples accumulated across all pixels.
    pub fn total_samples(&self) -> u64 {
        self.sample_count.iter().map(|&s| s as u64).sum()
    }
}

// ── Camera ───────────────────────────────────────────────────────────────────

/// A simple pinhole camera.
#[derive(Debug, Clone)]
pub struct Camera {
    /// Camera origin.
    pub origin: [f32; 3],
    lower_left_corner: [f32; 3],
    horizontal: [f32; 3],
    vertical: [f32; 3],
    lens_radius: f32,
    u: [f32; 3],
    v: [f32; 3],
}

impl Camera {
    /// Construct a camera.
    ///
    /// * `look_from` — eye position
    /// * `look_at` — target position
    /// * `vup` — world up vector
    /// * `vfov` — vertical field of view in degrees
    /// * `aspect_ratio` — image width / height
    /// * `aperture` — lens aperture (0 = pinhole)
    /// * `focus_dist` — focus distance
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        look_from: [f32; 3],
        look_at: [f32; 3],
        vup: [f32; 3],
        vfov: f32,
        aspect_ratio: f32,
        aperture: f32,
        focus_dist: f32,
    ) -> Self {
        let theta = vfov.to_radians();
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        let w = normalize(vsub(look_from, look_at));
        let u = normalize(cross(vup, w));
        let v = cross(w, u);

        let horizontal = vmul(u, viewport_width * focus_dist);
        let vertical = vmul(v, viewport_height * focus_dist);
        let lower_left_corner = vsub(
            vsub(vsub(look_from, vmul(horizontal, 0.5)), vmul(vertical, 0.5)),
            vmul(w, focus_dist),
        );

        Self {
            origin: look_from,
            lower_left_corner,
            horizontal,
            vertical,
            lens_radius: aperture / 2.0,
            u,
            v,
        }
    }

    /// Generate a ray for pixel coordinates `(s, t)` in \[0,1\]×\[0,1\].
    pub fn get_ray(&self, s: f32, t: f32, rng: &mut impl Rng) -> Ray {
        let rd = vmul(self.random_in_unit_disk(rng), self.lens_radius);
        let offset = vadd(vmul(self.u, rd[0]), vmul(self.v, rd[1]));
        let dir = vsub(
            vadd(
                vadd(self.lower_left_corner, vmul(self.horizontal, s)),
                vmul(self.vertical, t),
            ),
            vadd(self.origin, offset),
        );
        Ray::new(vadd(self.origin, offset), dir)
    }

    fn random_in_unit_disk(&self, rng: &mut impl Rng) -> [f32; 3] {
        loop {
            let p = [
                rng.random_range(-1.0f32..1.0),
                rng.random_range(-1.0f32..1.0),
                0.0,
            ];
            if dot(p, p) < 1.0 {
                return p;
            }
        }
    }
}

// ── PathTracerRenderer ───────────────────────────────────────────────────────

/// Renderer that progressively renders a scene into a `PathTracerBuffer`.
#[derive(Debug, Clone)]
pub struct PathTracerRenderer {
    /// Scene to render.
    pub scene: PathTracerScene,
    /// Camera.
    pub camera: Camera,
    /// Maximum ray bounce depth.
    pub max_depth: usize,
    /// Samples per pixel per call to `render_pass`.
    pub samples_per_pass: usize,
}

impl PathTracerRenderer {
    /// Create a new renderer.
    pub fn new(
        scene: PathTracerScene,
        camera: Camera,
        max_depth: usize,
        samples_per_pass: usize,
    ) -> Self {
        Self {
            scene,
            camera,
            max_depth,
            samples_per_pass,
        }
    }

    /// Render one progressive pass into `buffer`.
    ///
    /// Each pixel receives `self.samples_per_pass` new samples.
    pub fn render_pass(&self, buffer: &mut PathTracerBuffer) {
        let w = buffer.width;
        let h = buffer.height;
        let mut rng = rand::rng();
        for y in 0..h {
            for x in 0..w {
                let mut color = [0.0f32; 3];
                for _ in 0..self.samples_per_pass {
                    let u = (x as f32 + rng.random::<f32>()) / (w - 1) as f32;
                    let v = (y as f32 + rng.random::<f32>()) / (h - 1) as f32;
                    let ray = self.camera.get_ray(u, v, &mut rng);
                    let c = self.scene.trace(&ray, self.max_depth, &mut rng);
                    color[0] += c[0];
                    color[1] += c[1];
                    color[2] += c[2];
                }
                let inv = 1.0 / self.samples_per_pass as f32;
                buffer.add_sample(x, y, vmul(color, inv));
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_rng() -> impl Rng {
        rand::rng()
    }

    // ── Ray tests ────────────────────────────────────────────────────────

    #[test]
    fn test_ray_at_origin() {
        let r = Ray::new([0.0; 3], [1.0, 0.0, 0.0]);
        let p = r.at(0.0);
        assert_eq!(p, [0.0; 3]);
    }

    #[test]
    fn test_ray_at_t() {
        let r = Ray::new([1.0, 2.0, 3.0], [1.0, 0.0, 0.0]);
        let p = r.at(3.0);
        assert!((p[0] - 4.0).abs() < 1e-6);
        assert!((p[1] - 2.0).abs() < 1e-6);
        assert!((p[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_ray_at_negative_t() {
        let r = Ray::new([0.0; 3], [0.0, 1.0, 0.0]);
        let p = r.at(-2.0);
        assert!((p[1] - (-2.0)).abs() < 1e-6);
    }

    // ── Sphere tests ─────────────────────────────────────────────────────

    #[test]
    fn test_sphere_hit_center() {
        let s = Sphere::new([0.0, 0.0, -1.0], 0.5, 0);
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        let hit = s.hit(&r, 0.001, f32::INFINITY);
        assert!(hit.is_some());
        let rec = hit.unwrap();
        assert!(rec.t > 0.4 && rec.t < 0.6);
    }

    #[test]
    fn test_sphere_miss() {
        let s = Sphere::new([0.0, 0.0, -1.0], 0.5, 0);
        let r = Ray::new([0.0; 3], [0.0, 1.0, 0.0]);
        assert!(s.hit(&r, 0.001, f32::INFINITY).is_none());
    }

    #[test]
    fn test_sphere_hit_from_inside() {
        let s = Sphere::new([0.0; 3], 1.0, 0);
        let r = Ray::new([0.0; 3], [1.0, 0.0, 0.0]);
        let hit = s.hit(&r, 0.001, f32::INFINITY);
        assert!(hit.is_some());
        let rec = hit.unwrap();
        assert!(!rec.front_face);
    }

    #[test]
    fn test_sphere_normal_outward() {
        let s = Sphere::new([0.0; 3], 1.0, 0);
        let r = Ray::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        let hit = s.hit(&r, 0.001, f32::INFINITY).unwrap();
        assert!(hit.front_face);
        assert!((hit.normal[2] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_sphere_t_range_cull() {
        let s = Sphere::new([0.0, 0.0, -1.0], 0.5, 0);
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        // t_max < actual hit
        assert!(s.hit(&r, 0.001, 0.1).is_none());
    }

    // ── Triangle tests ───────────────────────────────────────────────────

    #[test]
    fn test_triangle_hit() {
        let tri = Triangle::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0);
        let r = Ray::new([0.0, 0.3, 1.0], [0.0, 0.0, -1.0]);
        let hit = tri.hit(&r, 0.001, f32::INFINITY);
        assert!(hit.is_some());
    }

    #[test]
    fn test_triangle_miss_outside() {
        let tri = Triangle::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0);
        let r = Ray::new([5.0, 5.0, 1.0], [0.0, 0.0, -1.0]);
        assert!(tri.hit(&r, 0.001, f32::INFINITY).is_none());
    }

    #[test]
    fn test_triangle_miss_parallel() {
        let tri = Triangle::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0);
        // Ray parallel to triangle plane
        let r = Ray::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        assert!(tri.hit(&r, 0.001, f32::INFINITY).is_none());
    }

    #[test]
    fn test_triangle_normal_direction() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0);
        // Normal should point in +Z for counter-clockwise winding
        assert!(tri.normal[2].abs() > 0.9);
    }

    #[test]
    fn test_triangle_hit_at_vertex() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0);
        // Ray aimed near edge midpoint
        let r = Ray::new([0.5, 0.5, 1.0], [0.0, 0.0, -1.0]);
        assert!(tri.hit(&r, 0.001, f32::INFINITY).is_some());
    }

    // ── Material tests ───────────────────────────────────────────────────

    #[test]
    fn test_lambertian_scatter() {
        let mat = Material::lambertian([0.8, 0.3, 0.3]);
        let ray = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        let hit = HitRecord {
            t: 1.0,
            point: [0.0, 0.0, -1.0],
            normal: [0.0, 0.0, 1.0],
            material_index: 0,
            front_face: true,
        };
        let mut rng = make_rng();
        let result = mat.scatter(&ray, &hit, &mut rng);
        assert!(result.is_some());
        let (_scattered, attenuation) = result.unwrap();
        assert!((attenuation[0] - 0.8).abs() < 1e-6);
    }

    #[test]
    fn test_metal_scatter() {
        let mat = Material::metal([0.8, 0.8, 0.8], 0.0);
        let ray = Ray::new([0.0; 3], normalize([1.0, -1.0, 0.0]));
        let hit = HitRecord {
            t: 1.0,
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            material_index: 0,
            front_face: true,
        };
        let mut rng = make_rng();
        let result = mat.scatter(&ray, &hit, &mut rng);
        assert!(result.is_some());
        let (scattered, _attenuation) = result.unwrap();
        // Reflected ray should go up (positive y)
        assert!(scattered.direction[1] > 0.0);
    }

    #[test]
    fn test_dielectric_scatter() {
        let mat = Material::dielectric(1.5);
        let ray = Ray::new([0.0, 0.0, 1.0], normalize([0.0, 0.0, -1.0]));
        let hit = HitRecord {
            t: 1.0,
            point: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            material_index: 0,
            front_face: true,
        };
        let mut rng = make_rng();
        let result = mat.scatter(&ray, &hit, &mut rng);
        assert!(result.is_some());
        let (_s, attn) = result.unwrap();
        assert!((attn[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_metal_fuzz_clamped() {
        let mat = Material::metal([1.0; 3], 5.0);
        if let MaterialType::Metal(f) = mat.kind {
            assert!(f <= 1.0);
        } else {
            panic!("expected Metal");
        }
    }

    // ── Scene tests ──────────────────────────────────────────────────────

    #[test]
    fn test_scene_add_material() {
        let mut scene = PathTracerScene::new();
        let idx = scene.add_material(Material::lambertian([1.0; 3]));
        assert_eq!(idx, 0);
        let idx2 = scene.add_material(Material::lambertian([0.5; 3]));
        assert_eq!(idx2, 1);
    }

    #[test]
    fn test_scene_hit_sphere() {
        let mut scene = PathTracerScene::new();
        let m = scene.add_material(Material::lambertian([0.5; 3]));
        scene.add_sphere(Sphere::new([0.0, 0.0, -1.0], 0.5, m));
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        assert!(scene.hit_scene(&r, 0.001, f32::INFINITY).is_some());
    }

    #[test]
    fn test_scene_miss() {
        let scene = PathTracerScene::new();
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        assert!(scene.hit_scene(&r, 0.001, f32::INFINITY).is_none());
    }

    #[test]
    fn test_scene_sky_color_up() {
        let scene = PathTracerScene::new();
        let r = Ray::new([0.0; 3], [0.0, 1.0, 0.0]);
        let c = scene.sky_color(&r);
        // Should be close to sky_top
        assert!(c[2] > 0.9);
    }

    #[test]
    fn test_scene_trace_no_hit() {
        let scene = PathTracerScene::new();
        let r = Ray::new([0.0; 3], [0.0, 1.0, 0.0]);
        let mut rng = make_rng();
        let c = scene.trace(&r, 5, &mut rng);
        // Should be sky colour
        assert!(c[2] > 0.0);
    }

    #[test]
    fn test_scene_trace_depth_zero() {
        let mut scene = PathTracerScene::new();
        let m = scene.add_material(Material::lambertian([0.5; 3]));
        scene.add_sphere(Sphere::new([0.0, 0.0, -1.0], 0.5, m));
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        let mut rng = make_rng();
        let c = scene.trace(&r, 0, &mut rng);
        assert_eq!(c, [0.0; 3]);
    }

    #[test]
    fn test_scene_closest_hit() {
        let mut scene = PathTracerScene::new();
        let m = scene.add_material(Material::lambertian([0.5; 3]));
        scene.add_sphere(Sphere::new([0.0, 0.0, -2.0], 0.5, m));
        scene.add_sphere(Sphere::new([0.0, 0.0, -1.0], 0.5, m));
        let r = Ray::new([0.0; 3], [0.0, 0.0, -1.0]);
        let hit = scene.hit_scene(&r, 0.001, f32::INFINITY).unwrap();
        // Closer sphere is at z=-1, so t ~ 0.5
        assert!(hit.t < 1.0);
    }

    // ── PathTracerBuffer tests ───────────────────────────────────────────

    #[test]
    fn test_buffer_new() {
        let buf = PathTracerBuffer::new(4, 4);
        assert_eq!(buf.width, 4);
        assert_eq!(buf.height, 4);
        assert_eq!(buf.total_samples(), 0);
    }

    #[test]
    fn test_buffer_add_and_get() {
        let mut buf = PathTracerBuffer::new(4, 4);
        buf.add_sample(1, 2, [0.6, 0.4, 0.2]);
        buf.add_sample(1, 2, [0.4, 0.6, 0.8]);
        let p = buf.get_pixel(1, 2);
        assert!((p[0] - 0.5).abs() < 1e-5);
        assert!((p[1] - 0.5).abs() < 1e-5);
        assert!((p[2] - 0.5).abs() < 1e-5);
    }

    #[test]
    fn test_buffer_zero_samples() {
        let buf = PathTracerBuffer::new(4, 4);
        let p = buf.get_pixel(0, 0);
        assert_eq!(p, [0.0; 3]);
    }

    #[test]
    fn test_buffer_to_rgb8_white() {
        let mut buf = PathTracerBuffer::new(1, 1);
        buf.add_sample(0, 0, [1.0; 3]);
        let rgb = buf.to_rgb8();
        assert_eq!(rgb.len(), 3);
        assert_eq!(rgb[0], 255);
    }

    #[test]
    fn test_buffer_to_rgb8_black() {
        let mut buf = PathTracerBuffer::new(1, 1);
        buf.add_sample(0, 0, [0.0; 3]);
        let rgb = buf.to_rgb8();
        assert_eq!(rgb[0], 0);
    }

    #[test]
    fn test_buffer_total_samples() {
        let mut buf = PathTracerBuffer::new(2, 2);
        buf.add_sample(0, 0, [1.0; 3]);
        buf.add_sample(0, 0, [1.0; 3]);
        buf.add_sample(1, 1, [0.5; 3]);
        assert_eq!(buf.total_samples(), 3);
    }

    #[test]
    fn test_buffer_clear() {
        let mut buf = PathTracerBuffer::new(2, 2);
        buf.add_sample(0, 0, [1.0; 3]);
        buf.clear();
        assert_eq!(buf.total_samples(), 0);
        assert_eq!(buf.get_pixel(0, 0), [0.0; 3]);
    }

    #[test]
    fn test_buffer_size() {
        let buf = PathTracerBuffer::new(8, 6);
        assert_eq!(buf.accumulator.len(), 48);
        assert_eq!(buf.sample_count.len(), 48);
    }

    // ── Camera tests ─────────────────────────────────────────────────────

    #[test]
    fn test_camera_get_ray_center() {
        let cam = Camera::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            90.0,
            1.0,
            0.0,
            1.0,
        );
        let mut rng = make_rng();
        let ray = cam.get_ray(0.5, 0.5, &mut rng);
        // Center ray should point roughly in -Z
        let d = normalize(ray.direction);
        assert!(d[2] < -0.9);
    }

    #[test]
    fn test_camera_origin() {
        let cam = Camera::new(
            [1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            60.0,
            1.5,
            0.0,
            1.0,
        );
        assert!((cam.origin[0] - 1.0).abs() < 1e-5);
    }

    // ── Renderer integration test ────────────────────────────────────────

    #[test]
    fn test_renderer_render_pass_small() {
        let mut scene = PathTracerScene::new();
        let m = scene.add_material(Material::lambertian([0.7, 0.3, 0.5]));
        scene.add_sphere(Sphere::new([0.0, 0.0, -1.0], 0.5, m));
        let cam = Camera::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            90.0,
            1.0,
            0.0,
            1.0,
        );
        let renderer = PathTracerRenderer::new(scene, cam, 3, 2);
        let mut buf = PathTracerBuffer::new(4, 4);
        renderer.render_pass(&mut buf);
        assert!(buf.total_samples() > 0);
        // Each pixel gets exactly one accumulated sample entry (averaged internally)
        assert_eq!(buf.total_samples(), (4 * 4) as u64);
    }

    #[test]
    fn test_renderer_rgb8_output_valid() {
        let mut scene = PathTracerScene::new();
        let m = scene.add_material(Material::lambertian([0.5; 3]));
        scene.add_sphere(Sphere::new([0.0, 0.0, -1.0], 0.5, m));
        let cam = Camera::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            90.0,
            1.0,
            0.0,
            1.0,
        );
        let renderer = PathTracerRenderer::new(scene, cam, 2, 1);
        let mut buf = PathTracerBuffer::new(8, 8);
        renderer.render_pass(&mut buf);
        let rgb = buf.to_rgb8();
        assert_eq!(rgb.len(), 8 * 8 * 3);
        // All values valid u8
        for &v in &rgb {
            let _ = v; // just confirm they exist
        }
    }

    // ── Vector helper tests ──────────────────────────────────────────────

    #[test]
    fn test_vadd() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = vadd(a, b);
        assert_eq!(c, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_vsub() {
        let a = [3.0, 2.0, 1.0];
        let b = [1.0, 1.0, 1.0];
        assert_eq!(vsub(a, b), [2.0, 1.0, 0.0]);
    }

    #[test]
    fn test_dot_orthogonal() {
        assert!((dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])).abs() < 1e-7);
    }

    #[test]
    fn test_cross_unit_vectors() {
        let k = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((k[2] - 1.0).abs() < 1e-7);
    }

    #[test]
    fn test_normalize_length() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize(v);
        let l = length(n);
        assert!((l - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_reflect_normal_incidence() {
        let d = [0.0, -1.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let r = reflect(d, n);
        assert!((r[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_schlick_zero_angle() {
        // cosine=0 → schlick = r0 + (1-r0)*(1-0)^5 = 1.0
        let s = schlick(0.0, 1.5);
        assert!((s - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_schlick_grazing() {
        // cosine=1 → schlick = r0 + (1-r0)*(1-1)^5 = r0
        let s = schlick(1.0, 1.5);
        let r0 = ((1.0 - 1.5f32) / (1.0 + 1.5)).powi(2);
        assert!((s - r0).abs() < 1e-5);
    }
}
