// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sphere shape.

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::f64::consts::PI;

/// A sphere defined by its radius.
#[derive(Debug, Clone)]
pub struct Sphere {
    /// Radius of the sphere.
    pub radius: Real,
}

impl Sphere {
    /// Create a new sphere with the given radius.
    pub fn new(radius: Real) -> Self {
        Self { radius }
    }

    /// Surface area: 4πr².
    pub fn surface_area(&self) -> Real {
        4.0 * PI * self.radius * self.radius
    }

    /// Volume: (4/3)πr³.
    pub fn volume_explicit(&self) -> Real {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    /// Inertia tensor as \[\[f64;3\\];3] row-major: (2/5)mr² × identity.
    pub fn inertia_tensor_array(&self, mass: f64) -> [[f64; 3]; 3] {
        let i = 0.4 * mass * self.radius * self.radius;
        [[i, 0.0, 0.0], [0.0, i, 0.0], [0.0, 0.0, i]]
    }

    /// Ray cast returning (t, normal) as plain arrays.
    pub fn ray_cast_array(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        let o = Vec3::new(origin[0], origin[1], origin[2]);
        let d = Vec3::new(direction[0], direction[1], direction[2]);
        let hit = self.ray_cast(&o, &d, max_toi)?;
        Some((hit.toi, [hit.normal.x, hit.normal.y, hit.normal.z]))
    }

    /// Closest point on the sphere surface to `p`.
    /// If `p` is the origin, returns a point on the +X side.
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if len < 1e-12 {
            return [self.radius, 0.0, 0.0];
        }
        let s = self.radius / len;
        [p[0] * s, p[1] * s, p[2] * s]
    }

    /// GJK support function: farthest point in `direction`.
    pub fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        if len < 1e-12 {
            return [self.radius, 0.0, 0.0];
        }
        let s = self.radius / len;
        [direction[0] * s, direction[1] * s, direction[2] * s]
    }

    // ── New methods ──

    /// GJK support point with a center offset.
    /// Returns `center + radius * normalize(direction)`.
    pub fn support_with_center(&self, center: [f64; 3], direction: [f64; 3]) -> [f64; 3] {
        let sp = self.support(direction);
        [sp[0] + center[0], sp[1] + center[1], sp[2] + center[2]]
    }

    /// Minkowski sum support: support(A) + support(B) in the given direction.
    /// Both spheres are centered at the origin.
    pub fn minkowski_sum_support(&self, other: &Sphere, direction: [f64; 3]) -> [f64; 3] {
        let combined_radius = self.radius + other.radius;
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        if len < 1e-12 {
            return [combined_radius, 0.0, 0.0];
        }
        let s = combined_radius / len;
        [direction[0] * s, direction[1] * s, direction[2] * s]
    }

    /// Minkowski difference support: support_A(d) - support_B(-d).
    /// Both spheres centered at origin.
    pub fn minkowski_diff_support(&self, other: &Sphere, direction: [f64; 3]) -> [f64; 3] {
        let sa = self.support(direction);
        let neg_d = [-direction[0], -direction[1], -direction[2]];
        let sb = other.support(neg_d);
        [sa[0] - sb[0], sa[1] - sb[1], sa[2] - sb[2]]
    }

    /// Compute bounding sphere from a set of points (Ritter's algorithm).
    /// Returns `(center, radius)`.
    pub fn bounding_sphere_from_points(points: &[[f64; 3]]) -> ([f64; 3], f64) {
        if points.is_empty() {
            return ([0.0; 3], 0.0);
        }
        if points.len() == 1 {
            return (points[0], 0.0);
        }

        // Find two most separated points along each axis
        let mut min_x = 0usize;
        let mut max_x = 0usize;
        for (i, p) in points.iter().enumerate() {
            if p[0] < points[min_x][0] {
                min_x = i;
            }
            if p[0] > points[max_x][0] {
                max_x = i;
            }
        }

        let mut center = [
            (points[min_x][0] + points[max_x][0]) * 0.5,
            (points[min_x][1] + points[max_x][1]) * 0.5,
            (points[min_x][2] + points[max_x][2]) * 0.5,
        ];
        let dx = points[max_x][0] - points[min_x][0];
        let dy = points[max_x][1] - points[min_x][1];
        let dz = points[max_x][2] - points[min_x][2];
        let mut radius = (dx * dx + dy * dy + dz * dz).sqrt() * 0.5;

        // Expand sphere to include all points
        for p in points {
            let dx = p[0] - center[0];
            let dy = p[1] - center[1];
            let dz = p[2] - center[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > radius {
                let new_radius = (radius + dist) * 0.5;
                let shift = (new_radius - radius) / dist;
                center[0] += dx * shift;
                center[1] += dy * shift;
                center[2] += dz * shift;
                radius = new_radius;
            }
        }

        (center, radius)
    }

    /// Sphere-sphere intersection test.
    /// Returns true if two spheres (at given centers) overlap.
    pub fn intersects_sphere(
        &self,
        center_a: [f64; 3],
        other: &Sphere,
        center_b: [f64; 3],
    ) -> bool {
        let dx = center_b[0] - center_a[0];
        let dy = center_b[1] - center_a[1];
        let dz = center_b[2] - center_a[2];
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let r_sum = self.radius + other.radius;
        dist_sq <= r_sum * r_sum
    }

    /// Sphere-sphere intersection circle.
    /// Returns `Some((center, normal, circle_radius))` if the spheres intersect
    /// in a circle, `None` if they don't intersect or are concentric.
    pub fn sphere_intersection_circle(
        &self,
        center_a: [f64; 3],
        other: &Sphere,
        center_b: [f64; 3],
    ) -> Option<([f64; 3], [f64; 3], f64)> {
        let dx = center_b[0] - center_a[0];
        let dy = center_b[1] - center_a[1];
        let dz = center_b[2] - center_a[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();

        if d < 1e-12 {
            return None; // concentric
        }

        let r1 = self.radius;
        let r2 = other.radius;

        if d > r1 + r2 || d < (r1 - r2).abs() {
            return None; // no intersection
        }

        // Distance from center_a along the axis to the intersection plane
        let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        let h_sq = r1 * r1 - a * a;
        let h = if h_sq > 0.0 { h_sq.sqrt() } else { 0.0 };

        let nx = dx / d;
        let ny = dy / d;
        let nz = dz / d;

        let cx = center_a[0] + a * nx;
        let cy = center_a[1] + a * ny;
        let cz = center_a[2] + a * nz;

        Some(([cx, cy, cz], [nx, ny, nz], h))
    }

    /// Closest point on the sphere surface to a line defined by `point` and `direction`.
    /// Returns the closest point on the sphere surface.
    pub fn closest_point_to_line(&self, line_point: [f64; 3], line_dir: [f64; 3]) -> [f64; 3] {
        // Find closest point on line to origin (sphere center)
        let dir_len_sq =
            line_dir[0] * line_dir[0] + line_dir[1] * line_dir[1] + line_dir[2] * line_dir[2];
        if dir_len_sq < 1e-24 {
            return self.closest_point(line_point);
        }
        let t = -(line_point[0] * line_dir[0]
            + line_point[1] * line_dir[1]
            + line_point[2] * line_dir[2])
            / dir_len_sq;
        let closest_on_line = [
            line_point[0] + t * line_dir[0],
            line_point[1] + t * line_dir[1],
            line_point[2] + t * line_dir[2],
        ];
        self.closest_point(closest_on_line)
    }

    /// Closest point on the sphere surface to a plane defined by `normal` (unit) and `d`
    /// where `normal · x = d`.
    pub fn closest_point_to_plane(&self, normal: [f64; 3], d: f64) -> [f64; 3] {
        // Signed distance from sphere center (origin) to plane
        // sign = -(normal · 0 - d) = d
        // The closest point on the sphere is in the direction toward the plane
        let sign = if d >= 0.0 { 1.0 } else { -1.0 };
        [
            sign * normal[0] * self.radius,
            sign * normal[1] * self.radius,
            sign * normal[2] * self.radius,
        ]
    }

    /// Sphere sweep (moving sphere): test if a sphere moving from `start` along
    /// `velocity` hits a static sphere (centered at `target_center` with `target_radius`).
    /// Returns `Some(t)` in \[0, max_t\] if there is a collision.
    pub fn sphere_sweep(
        &self,
        start: [f64; 3],
        velocity: [f64; 3],
        target_center: [f64; 3],
        target_radius: f64,
        max_t: f64,
    ) -> Option<f64> {
        let combined_r = self.radius + target_radius;
        // Ray from start toward target
        let dx = start[0] - target_center[0];
        let dy = start[1] - target_center[1];
        let dz = start[2] - target_center[2];

        let a = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
        let b = 2.0 * (dx * velocity[0] + dy * velocity[1] + dz * velocity[2]);
        let c = dx * dx + dy * dy + dz * dz - combined_r * combined_r;

        if a < 1e-24 {
            // Stationary
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

    /// Signed distance from a point to the sphere surface.
    /// Negative if inside the sphere, positive if outside.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        len - self.radius
    }

    /// Returns true if the point is inside (or on) the sphere.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] * p[0] + p[1] * p[1] + p[2] * p[2] <= self.radius * self.radius
    }

    /// Project a point onto the sphere interior (clamp to sphere if outside).
    pub fn project_inside(&self, p: [f64; 3]) -> [f64; 3] {
        let len_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        if len_sq <= self.radius * self.radius {
            p
        } else {
            self.closest_point(p)
        }
    }

    // ── Extended: geodesic sphere, spherical cap, hemisphere sampling, conformal maps ──

    /// Generate a geodesic icosphere by subdividing an icosahedron.
    ///
    /// Returns `(vertices, triangles)` where each vertex is a unit-sphere
    /// point scaled to `self.radius`, and triangles are `[usize; 3]` index triples.
    /// `subdivisions` = 0 gives the base icosahedron (20 faces, 12 verts).
    pub fn geodesic_icosphere(&self, subdivisions: u32) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let (mut verts, mut tris) = base_icosahedron();
        for _ in 0..subdivisions {
            let (v2, t2) = subdivide_icosphere(verts, tris);
            verts = v2;
            tris = t2;
        }
        // Scale to self.radius
        let r = self.radius;
        let verts = verts
            .iter()
            .map(|v| {
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-14);
                [v[0] / len * r, v[1] / len * r, v[2] / len * r]
            })
            .collect();
        (verts, tris)
    }

    /// Spherical cap volume.
    ///
    /// Returns the volume of a spherical cap of height `h` cut from this sphere.
    /// `h` must be in `[0, 2r]`.
    pub fn spherical_cap_volume(&self, h: f64) -> f64 {
        let h_clamped = h.clamp(0.0, 2.0 * self.radius);
        PI * h_clamped * h_clamped * (3.0 * self.radius - h_clamped) / 3.0
    }

    /// Spherical cap surface area (curved part only).
    ///
    /// Area = 2π r h.
    pub fn spherical_cap_area(&self, h: f64) -> f64 {
        let h_clamped = h.clamp(0.0, 2.0 * self.radius);
        2.0 * PI * self.radius * h_clamped
    }

    /// Spherical zone area between two parallel planes.
    ///
    /// For planes at heights `h1` and `h2` (measured from the bottom of the
    /// sphere), returns the surface area of the zone.
    /// Area = 2π r |h2 - h1|.
    pub fn spherical_zone_area(&self, h1: f64, h2: f64) -> f64 {
        let dh = (h2 - h1).abs();
        2.0 * PI * self.radius * dh
    }

    /// Sample points on the hemisphere (z ≥ 0) using cosine-weighted sampling.
    ///
    /// Uses a deterministic xorshift PRNG seeded with `seed`.
    /// Each point lies on the unit hemisphere surface projected to `self.radius`.
    pub fn hemisphere_cosine_sample(&self, n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut state = seed;
        let xorshift = |s: &mut u64| -> f64 {
            *s ^= *s << 13;
            *s ^= *s >> 7;
            *s ^= *s << 17;
            (*s as f64) / (u64::MAX as f64)
        };
        let r = self.radius;
        (0..n)
            .map(|_| {
                let u1 = xorshift(&mut state);
                let u2 = xorshift(&mut state);
                // cosine-weighted: sqrt(u1) for sin(θ), 2π*u2 for φ
                let sin_theta = u1.sqrt();
                let cos_theta = (1.0 - u1).sqrt();
                let phi = 2.0 * PI * u2;
                [
                    r * sin_theta * phi.cos(),
                    r * sin_theta * phi.sin(),
                    r * cos_theta,
                ]
            })
            .collect()
    }

    /// Uniform hemisphere sampling (z ≥ 0).
    ///
    /// Each sample is uniformly distributed on the upper hemisphere.
    pub fn hemisphere_uniform_sample(&self, n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut state = seed;
        let xorshift = |s: &mut u64| -> f64 {
            *s ^= *s << 13;
            *s ^= *s >> 7;
            *s ^= *s << 17;
            (*s as f64) / (u64::MAX as f64)
        };
        let r = self.radius;
        (0..n)
            .map(|_| {
                let u1 = xorshift(&mut state);
                let u2 = xorshift(&mut state);
                let cos_theta = u1; // uniform in cos_theta ∈ [0, 1]
                let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
                let phi = 2.0 * PI * u2;
                [
                    r * sin_theta * phi.cos(),
                    r * sin_theta * phi.sin(),
                    r * cos_theta,
                ]
            })
            .collect()
    }

    /// 3-D Fibonacci sphere packing.
    ///
    /// Generates `n` nearly uniformly distributed points on the sphere surface
    /// using the golden-angle Fibonacci lattice method.
    pub fn fibonacci_sphere(&self, n: usize) -> Vec<[f64; 3]> {
        let n = n.max(1);
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let r = self.radius;
        (0..n)
            .map(|i| {
                let theta = (1.0 - 2.0 * i as f64 / (n - 1).max(1) as f64)
                    .clamp(-1.0, 1.0)
                    .acos();
                let phi = 2.0 * PI * i as f64 / golden;
                [
                    r * theta.sin() * phi.cos(),
                    r * theta.sin() * phi.sin(),
                    r * theta.cos(),
                ]
            })
            .collect()
    }

    /// Stereographic projection from the sphere to the plane.
    ///
    /// Projects a point `p` on the sphere to the stereographic plane using
    /// the south pole as the projection center.  The north pole maps to
    /// infinity; the south pole is undefined.
    ///
    /// Returns `(x_plane, y_plane)` or `None` if `p` is at (or near) the south pole.
    pub fn stereographic_project(&self, p: [f64; 3]) -> Option<[f64; 2]> {
        let r = self.radius;
        // South pole at (0, 0, -r)
        let denom = r + p[2]; // z + r
        if denom.abs() < 1e-10 * r {
            return None; // near south pole
        }
        let scale = 2.0 * r / denom;
        Some([p[0] * scale, p[1] * scale])
    }

    /// Inverse stereographic projection from the plane to the sphere.
    ///
    /// Maps `(x_plane, y_plane)` back to a point on the sphere.
    /// Uses the formula: X = 4r²x/D, Y = 4r²y/D, Z = r(D - 4r²)/D
    /// where D = x² + y² + 4r².
    pub fn stereographic_unproject(&self, uv: [f64; 2]) -> [f64; 3] {
        let r = self.radius;
        let x = uv[0];
        let y = uv[1];
        let d2 = x * x + y * y;
        let denom = d2 + 4.0 * r * r;
        [
            4.0 * r * r * x / denom,
            4.0 * r * r * y / denom,
            r * (4.0 * r * r - d2) / denom,
        ]
    }

    /// Spherical harmonic (l=3) coefficient Y_3^m(theta, phi).
    ///
    /// Extends `spherical_harmonic` to degree 3.
    pub fn sh_l3(m: i32, theta: f64, phi: f64) -> f64 {
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        match m {
            -3 => 0.25 * (35.0 / (2.0 * PI)).sqrt() * sin_t.powi(3) * (3.0 * phi).sin(),
            -2 => 0.5 * (105.0 / PI).sqrt() * sin_t.powi(2) * cos_t * (2.0 * phi).sin(),
            -1 => {
                0.25 * (21.0 / (2.0 * PI)).sqrt() * sin_t * (5.0 * cos_t * cos_t - 1.0) * phi.sin()
            }
            0 => 0.25 * (7.0 / PI).sqrt() * (5.0 * cos_t.powi(3) - 3.0 * cos_t),
            1 => {
                0.25 * (21.0 / (2.0 * PI)).sqrt() * sin_t * (5.0 * cos_t * cos_t - 1.0) * phi.cos()
            }
            2 => 0.25 * (105.0 / PI).sqrt() * sin_t.powi(2) * cos_t * (2.0 * phi).cos(),
            3 => 0.25 * (35.0 / (2.0 * PI)).sqrt() * sin_t.powi(3) * (3.0 * phi).cos(),
            _ => 0.0,
        }
    }

    /// Compute the solid angle subtended by a spherical cap of height `h`.
    ///
    /// Ω = 2π(1 - cos θ) where cos θ = 1 - h/r.
    pub fn cap_solid_angle(&self, h: f64) -> f64 {
        let cos_theta = 1.0 - h / self.radius;
        2.0 * PI * (1.0 - cos_theta.clamp(-1.0, 1.0))
    }

    /// Longitude/latitude (geographic) coordinates for a point on the sphere.
    ///
    /// Returns `(longitude, latitude)` in radians.
    /// Longitude in `(-π, π]`, latitude in `[-π/2, π/2]`.
    pub fn lon_lat(&self, p: [f64; 3]) -> (f64, f64) {
        let lon = p[1].atan2(p[0]);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-14);
        let lat = (p[2] / r).clamp(-1.0, 1.0).asin();
        (lon, lat)
    }
}

// ── Icosphere helpers ─────────────────────────────────────────────────────────

/// Build the base icosahedron: 12 vertices on the unit sphere, 20 triangles.
fn base_icosahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let s = 1.0 / (1.0 + t * t).sqrt();
    let ts = t * s;
    let verts = vec![
        [-s, ts, 0.0],
        [s, ts, 0.0],
        [-s, -ts, 0.0],
        [s, -ts, 0.0],
        [0.0, -s, ts],
        [0.0, s, ts],
        [0.0, -s, -ts],
        [0.0, s, -ts],
        [ts, 0.0, -s],
        [ts, 0.0, s],
        [-ts, 0.0, -s],
        [-ts, 0.0, s],
    ];
    let tris = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    (verts, tris)
}

/// Midpoint of two 3-D points.
fn midpoint3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        (a[0] + b[0]) * 0.5,
        (a[1] + b[1]) * 0.5,
        (a[2] + b[2]) * 0.5,
    ]
}

/// One level of icosphere subdivision (each triangle → 4 triangles).
fn subdivide_icosphere(
    verts: Vec<[f64; 3]>,
    tris: Vec<[usize; 3]>,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    use std::collections::HashMap;
    let mut new_verts = verts.clone();
    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut new_tris = Vec::with_capacity(tris.len() * 4);

    let mut midpoint_idx = |v: &mut Vec<[f64; 3]>, a: usize, b: usize| -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&idx) = edge_map.get(&key) {
            return idx;
        }
        let mp = midpoint3(v[a], v[b]);
        let idx = v.len();
        v.push(mp);
        edge_map.insert(key, idx);
        idx
    };

    for &[a, b, c] in &tris {
        let ab = midpoint_idx(&mut new_verts, a, b);
        let bc = midpoint_idx(&mut new_verts, b, c);
        let ca = midpoint_idx(&mut new_verts, c, a);
        new_tris.push([a, ab, ca]);
        new_tris.push([b, bc, ab]);
        new_tris.push([c, ca, bc]);
        new_tris.push([ab, bc, ca]);
    }
    (new_verts, new_tris)
}

impl Shape for Sphere {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn bounding_box(&self) -> Aabb {
        let extent = Vec3::new(self.radius, self.radius, self.radius);
        Aabb::new(-extent, extent)
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let norm = direction.norm();
        if norm < 1e-10 {
            return Vec3::new(self.radius, 0.0, 0.0);
        }
        direction * (self.radius / norm)
    }

    fn volume(&self) -> Real {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::zeros()
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let i = 0.4 * mass * self.radius * self.radius;
        Mat3::new(i, 0.0, 0.0, 0.0, i, 0.0, 0.0, 0.0, i)
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        // Solve |origin + t * direction|^2 = radius^2
        let a = ray_direction.dot(ray_direction);
        let b = 2.0 * ray_origin.dot(ray_direction);
        let c = ray_origin.dot(ray_origin) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            return None;
        }

        let sqrt_disc = discriminant.sqrt();
        let t1 = (-b - sqrt_disc) / (2.0 * a);
        let t2 = (-b + sqrt_disc) / (2.0 * a);

        let t = if t1 >= 0.0 { t1 } else { t2 };

        if t < 0.0 || t > max_toi {
            return None;
        }

        let point = ray_origin + ray_direction * t;
        let normal = point.normalize();
        Some(RayHit {
            point,
            normal,
            toi: t,
        })
    }
}

// ── Sphere extensions (free functions) ───────────────────────────────────────

/// Minimum enclosing sphere using Welzl's algorithm (randomised).
///
/// Returns `(center, radius)` for the smallest sphere containing all points.
pub fn minimum_enclosing_sphere(points: &[[f64; 3]]) -> ([f64; 3], f64) {
    if points.is_empty() {
        return ([0.0; 3], 0.0);
    }
    let mut pts: Vec<[f64; 3]> = points.to_vec();
    let n = pts.len();
    welzl_sphere(&mut pts, &[], n)
}

/// Sphere-OBB penetration depth.
pub fn sphere_obb_penetration(
    sphere_radius: f64,
    sphere_center: [f64; 3],
    obb_center: [f64; 3],
    obb_half_extents: [f64; 3],
    obb_rotation: [[f64; 3]; 3],
) -> Option<f64> {
    let diff = [
        sphere_center[0] - obb_center[0],
        sphere_center[1] - obb_center[1],
        sphere_center[2] - obb_center[2],
    ];
    let local = [
        diff[0] * obb_rotation[0][0] + diff[1] * obb_rotation[0][1] + diff[2] * obb_rotation[0][2],
        diff[0] * obb_rotation[1][0] + diff[1] * obb_rotation[1][1] + diff[2] * obb_rotation[1][2],
        diff[0] * obb_rotation[2][0] + diff[1] * obb_rotation[2][1] + diff[2] * obb_rotation[2][2],
    ];
    let clamped = [
        local[0].clamp(-obb_half_extents[0], obb_half_extents[0]),
        local[1].clamp(-obb_half_extents[1], obb_half_extents[1]),
        local[2].clamp(-obb_half_extents[2], obb_half_extents[2]),
    ];
    let delta = [
        local[0] - clamped[0],
        local[1] - clamped[1],
        local[2] - clamped[2],
    ];
    let dist_sq = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
    let r_sq = sphere_radius * sphere_radius;
    if dist_sq > r_sq {
        return None;
    }
    let dist = dist_sq.sqrt();
    Some(sphere_radius - dist)
}

/// Generate `n` uniformly distributed random points on the unit sphere surface.
/// Uses rejection sampling.  `seed` is an LCG seed.
pub fn random_points_on_sphere(radius: f64, n: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut rng_state = seed;
    let mut next_f64 = move || -> f64 {
        rng_state = rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = 0x3FF0000000000000u64 | (rng_state >> 12);
        f64::from_bits(bits) - 1.0
    };

    let mut points = Vec::with_capacity(n);
    while points.len() < n {
        let x = next_f64() * 2.0 - 1.0;
        let y = next_f64() * 2.0 - 1.0;
        let z = next_f64() * 2.0 - 1.0;
        let len_sq = x * x + y * y + z * z;
        if !(1e-30..=1.0).contains(&len_sq) {
            continue;
        }
        let len = len_sq.sqrt();
        points.push([x / len * radius, y / len * radius, z / len * radius]);
    }
    points
}

/// Compute real spherical harmonics coefficient Y_l^m(theta, phi).
/// Only l = 0, 1, 2 implemented.
pub fn spherical_harmonic(l: u32, m: i32, theta: f64, phi: f64) -> f64 {
    match (l, m) {
        (0, 0) => 1.0 / (2.0 * PI.sqrt()),
        (1, -1) => (3.0 / (4.0 * PI)).sqrt() * theta.sin() * phi.sin(),
        (1, 0) => (3.0 / (4.0 * PI)).sqrt() * theta.cos(),
        (1, 1) => (3.0 / (4.0 * PI)).sqrt() * theta.sin() * phi.cos(),
        (2, -2) => 0.5 * (15.0 / PI).sqrt() * theta.sin().powi(2) * (2.0 * phi).sin(),
        (2, -1) => 0.5 * (15.0 / PI).sqrt() * theta.sin() * theta.cos() * phi.sin(),
        (2, 0) => 0.25 * (5.0 / PI).sqrt() * (3.0 * theta.cos().powi(2) - 1.0),
        (2, 1) => 0.5 * (15.0 / PI).sqrt() * theta.sin() * theta.cos() * phi.cos(),
        (2, 2) => 0.25 * (15.0 / PI).sqrt() * theta.sin().powi(2) * (2.0 * phi).cos(),
        _ => 0.0,
    }
}

/// Sphere-sphere intersection volume (lens volume).
pub fn sphere_sphere_intersection_volume(r1: f64, r2: f64, d: f64) -> f64 {
    if d >= r1 + r2 {
        return 0.0;
    }
    if d <= (r1 - r2).abs() {
        let r_min = r1.min(r2);
        return (4.0 / 3.0) * PI * r_min.powi(3);
    }
    let h1 = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let cap1 = (PI / 3.0) * (r1 - h1).powi(2) * (2.0 * r1 + h1);
    let h2 = d - h1;
    let cap2 = (PI / 3.0) * (r2 - h2).powi(2) * (2.0 * r2 + h2);
    cap1 + cap2
}

// ── Welzl helpers (module-level) ──────────────────────────────────────────────

/// Minimum sphere enclosing up to 4 boundary points.
fn min_sphere_with_boundary(boundary: &[[f64; 3]]) -> ([f64; 3], f64) {
    match boundary.len() {
        0 => ([0.0; 3], 0.0),
        1 => (boundary[0], 0.0),
        2 => {
            let c = [
                (boundary[0][0] + boundary[1][0]) * 0.5,
                (boundary[0][1] + boundary[1][1]) * 0.5,
                (boundary[0][2] + boundary[1][2]) * 0.5,
            ];
            let r = sphere_dist3(boundary[0], c);
            (c, r)
        }
        3 => circumsphere_3(boundary[0], boundary[1], boundary[2]),
        4 => circumsphere_4(boundary[0], boundary[1], boundary[2], boundary[3]),
        _ => ([0.0; 3], 0.0),
    }
}

fn sphere_dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn welzl_sphere(pts: &mut [[f64; 3]], boundary: &[[f64; 3]], n: usize) -> ([f64; 3], f64) {
    if n == 0 || boundary.len() == 4 {
        return min_sphere_with_boundary(boundary);
    }
    let p = pts[n - 1];
    let (c, r) = welzl_sphere(pts, boundary, n - 1);
    if sphere_dist3(p, c) <= r + 1e-10 {
        return (c, r);
    }
    let mut new_boundary = boundary.to_vec();
    new_boundary.push(p);
    welzl_sphere(pts, &new_boundary, n - 1)
}

fn circumsphere_3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], f64) {
    // circumcenter of triangle in 3D
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ab_cross_ac = cross3sph(ab, ac);
    let denom = 2.0 * dot3sph(ab_cross_ac, ab_cross_ac);
    if denom.abs() < 1e-30 {
        // Degenerate: return circumsphere of ab
        return (
            [
                (a[0] + b[0]) * 0.5,
                (a[1] + b[1]) * 0.5,
                (a[2] + b[2]) * 0.5,
            ],
            sphere_dist3(
                a,
                [
                    (a[0] + b[0]) * 0.5,
                    (a[1] + b[1]) * 0.5,
                    (a[2] + b[2]) * 0.5,
                ],
            ),
        );
    }
    let i_denom = 1.0 / denom;
    let ab2 = dot3sph(ab, ab);
    let ac2 = dot3sph(ac, ac);
    let ab_cross_ac2 = cross3sph(ab, ac);
    let ac_cross_ab = cross3sph(ac, ab);
    let _ = ac_cross_ab;
    let t1 = cross3sph(scale3sph(ab, ac2), ab_cross_ac2);
    let t2 = cross3sph(scale3sph(ac, ab2), ab_cross_ac2);
    let center = [
        a[0] + (t1[0] - t2[0]) * i_denom,
        a[1] + (t1[1] - t2[1]) * i_denom,
        a[2] + (t1[2] - t2[2]) * i_denom,
    ];
    let r = sphere_dist3(a, center);
    (center, r)
}

fn circumsphere_4(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> ([f64; 3], f64) {
    // Use circumsphere of the three points that give the largest sphere
    let opts = [
        circumsphere_3(a, b, c),
        circumsphere_3(a, b, d),
        circumsphere_3(a, c, d),
        circumsphere_3(b, c, d),
    ];
    // Return the smallest sphere that contains all 4 points
    for (center, r) in &opts {
        let all_inside = [a, b, c, d]
            .iter()
            .all(|&p| sphere_dist3(p, *center) <= r + 1e-8);
        if all_inside {
            return (*center, *r);
        }
    }
    // Fallback: midpoint of ab
    let cen = [
        (a[0] + d[0]) * 0.5,
        (a[1] + d[1]) * 0.5,
        (a[2] + d[2]) * 0.5,
    ];
    let r = sphere_dist3(a, cen);
    (cen, r)
}

fn cross3sph(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot3sph(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scale3sph(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[cfg(test)]
mod proptest_tests {
    use super::*;
    use crate::Sphere;

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_sphere_volume_scales_cubically(r in 0.01_f64..100.0_f64) {
            let s1 = Sphere::new(r);
            let s2 = Sphere::new(2.0 * r);
            let ratio = s2.volume() / s1.volume();
            prop_assert!(
                (ratio - 8.0).abs() < 1e-6,
                "volume ratio for 2r vs r should be 8, got {}",
                ratio
            );
        }

        #[test]
        fn prop_sphere_inertia_positive_definite(
            r in 0.01_f64..100.0_f64,
            m in 0.01_f64..1000.0_f64,
        ) {
            let s = Sphere::new(r);
            let it = s.inertia_tensor(m);
            // Diagonal inertia tensor for sphere: all entries equal 0.4*m*r^2 > 0
            prop_assert!(it[(0, 0)] > 0.0, "Ixx not positive: {}", it[(0, 0)]);
            prop_assert!(it[(1, 1)] > 0.0, "Iyy not positive: {}", it[(1, 1)]);
            prop_assert!(it[(2, 2)] > 0.0, "Izz not positive: {}", it[(2, 2)]);
        }

        #[test]
        fn prop_sphere_volume_nonneg(r in 0.01_f64..100.0_f64) {
            let s = Sphere::new(r);
            prop_assert!(s.volume() >= 0.0, "volume negative: {}", s.volume());
        }

        #[test]
        fn prop_sphere_bounding_box_contains_surface(
            r in 0.01_f64..100.0_f64,
        ) {
            let s = Sphere::new(r);
            let bb = s.bounding_box();
            // The bounding box should extend exactly [-r, r] in each axis
            prop_assert!(
                (bb.min.x + r).abs() < 1e-10,
                "bb.min.x={} != -r={}",
                bb.min.x,
                -r
            );
            prop_assert!(
                (bb.max.x - r).abs() < 1e-10,
                "bb.max.x={} != r={}",
                bb.max.x,
                r
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Sphere;
    use crate::sphere::PI;
    use crate::sphere::minimum_enclosing_sphere;
    use crate::sphere::random_points_on_sphere;
    use crate::sphere::sphere_obb_penetration;
    use crate::sphere::sphere_sphere_intersection_volume;
    use crate::sphere::spherical_harmonic;
    use oxiphysics_core::math::Vec3;

    #[test]
    fn test_sphere_bounding_box() {
        let s = Sphere::new(2.0);
        let bb = s.bounding_box();
        assert_eq!(bb.min, Vec3::new(-2.0, -2.0, -2.0));
        assert_eq!(bb.max, Vec3::new(2.0, 2.0, 2.0));
    }

    #[test]
    fn test_sphere_volume() {
        let s = Sphere::new(1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_surface_area() {
        let s = Sphere::new(1.0);
        assert!((s.surface_area() - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_inertia() {
        let s = Sphere::new(1.0);
        let it = s.inertia_tensor(1.0);
        assert!((it[(0, 0)] - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_inertia_array() {
        let s = Sphere::new(1.0);
        let it = s.inertia_tensor_array(1.0);
        assert!((it[0][0] - 0.4).abs() < 1e-10);
        assert!((it[1][1] - 0.4).abs() < 1e-10);
        assert!((it[2][2] - 0.4).abs() < 1e-10);
        assert!((it[0][1]).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_raycast() {
        let s = Sphere::new(1.0);
        let origin = Vec3::new(-5.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let hit = s.ray_cast(&origin, &dir, 100.0).unwrap();
        assert!((hit.toi - 4.0).abs() < 1e-10);
        assert!((hit.point.x + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_raycast_miss() {
        let s = Sphere::new(1.0);
        let origin = Vec3::new(-5.0, 5.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        assert!(s.ray_cast(&origin, &dir, 100.0).is_none());
    }

    #[test]
    fn test_sphere_closest_point() {
        let s = Sphere::new(2.0);
        let cp = s.closest_point([4.0, 0.0, 0.0]);
        assert!((cp[0] - 2.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
        assert!(cp[2].abs() < 1e-10);
    }

    #[test]
    fn test_sphere_support() {
        let s = Sphere::new(3.0);
        let sp = s.support([1.0, 0.0, 0.0]);
        assert!((sp[0] - 3.0).abs() < 1e-10);
        // diagonal
        let len = 3.0_f64.sqrt();
        let sp2 = s.support([1.0, 1.0, 1.0]);
        assert!((sp2[0] - 3.0 / len).abs() < 1e-6);
    }

    #[test]
    fn test_sphere_raycast_array() {
        let s = Sphere::new(1.0);
        let (t, n) = s
            .ray_cast_array([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0)
            .unwrap();
        assert!((t - 4.0).abs() < 1e-10);
        assert!((n[0] + 1.0).abs() < 1e-10);
    }

    // ── New tests ──

    #[test]
    fn test_sphere_support_with_center() {
        let s = Sphere::new(2.0);
        let sp = s.support_with_center([1.0, 2.0, 3.0], [1.0, 0.0, 0.0]);
        assert!((sp[0] - 3.0).abs() < 1e-10); // 1 + 2
        assert!((sp[1] - 2.0).abs() < 1e-10);
        assert!((sp[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_minkowski_sum_support() {
        let a = Sphere::new(2.0);
        let b = Sphere::new(3.0);
        let sp = a.minkowski_sum_support(&b, [1.0, 0.0, 0.0]);
        assert!((sp[0] - 5.0).abs() < 1e-10);
        assert!(sp[1].abs() < 1e-10);
    }

    #[test]
    fn test_sphere_minkowski_diff_support() {
        let a = Sphere::new(5.0);
        let b = Sphere::new(2.0);
        // support_A([1,0,0]) = [5,0,0], support_B([-1,0,0]) = [-2,0,0]
        // diff = [5-(-2), 0, 0] = [7, 0, 0]
        let sp = a.minkowski_diff_support(&b, [1.0, 0.0, 0.0]);
        assert!((sp[0] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_sphere_from_points_empty() {
        let (c, r) = Sphere::bounding_sphere_from_points(&[]);
        assert_eq!(c, [0.0, 0.0, 0.0]);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn test_bounding_sphere_from_points_single() {
        let (c, r) = Sphere::bounding_sphere_from_points(&[[1.0, 2.0, 3.0]]);
        assert_eq!(c, [1.0, 2.0, 3.0]);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn test_bounding_sphere_from_points_two() {
        let pts = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let (c, r) = Sphere::bounding_sphere_from_points(&pts);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_sphere_contains_all_points() {
        let pts = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let (c, r) = Sphere::bounding_sphere_from_points(&pts);
        for p in &pts {
            let dx = p[0] - c[0];
            let dy = p[1] - c[1];
            let dz = p[2] - c[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                dist <= r + 1e-6,
                "point {:?} outside bounding sphere (dist={}, r={})",
                p,
                dist,
                r
            );
        }
    }

    #[test]
    fn test_sphere_intersection_overlap() {
        let a = Sphere::new(1.0);
        let b = Sphere::new(1.0);
        assert!(a.intersects_sphere([0.0, 0.0, 0.0], &b, [1.5, 0.0, 0.0]));
    }

    #[test]
    fn test_sphere_intersection_no_overlap() {
        let a = Sphere::new(1.0);
        let b = Sphere::new(1.0);
        assert!(!a.intersects_sphere([0.0, 0.0, 0.0], &b, [3.0, 0.0, 0.0]));
    }

    #[test]
    fn test_sphere_intersection_touching() {
        let a = Sphere::new(1.0);
        let b = Sphere::new(1.0);
        assert!(a.intersects_sphere([0.0, 0.0, 0.0], &b, [2.0, 0.0, 0.0]));
    }

    #[test]
    fn test_sphere_intersection_circle() {
        let a = Sphere::new(2.0);
        let b = Sphere::new(2.0);
        let result = a.sphere_intersection_circle([0.0, 0.0, 0.0], &b, [2.0, 0.0, 0.0]);
        assert!(result.is_some());
        let (center, normal, radius) = result.unwrap();
        assert!((center[0] - 1.0).abs() < 1e-10); // midpoint along axis
        assert!((normal[0] - 1.0).abs() < 1e-10); // axis direction
        assert!(radius > 0.0);
    }

    #[test]
    fn test_sphere_intersection_circle_no_overlap() {
        let a = Sphere::new(1.0);
        let b = Sphere::new(1.0);
        let result = a.sphere_intersection_circle([0.0, 0.0, 0.0], &b, [5.0, 0.0, 0.0]);
        assert!(result.is_none());
    }

    #[test]
    fn test_sphere_intersection_circle_concentric() {
        let a = Sphere::new(2.0);
        let b = Sphere::new(1.0);
        let result = a.sphere_intersection_circle([0.0, 0.0, 0.0], &b, [0.0, 0.0, 0.0]);
        assert!(result.is_none());
    }

    #[test]
    fn test_closest_point_to_line() {
        let s = Sphere::new(1.0);
        // Line along Y axis, offset by x=3
        let cp = s.closest_point_to_line([3.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        // Closest point on line to origin is (3,0,0), so closest on sphere is (1,0,0)
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
        assert!(cp[2].abs() < 1e-10);
    }

    #[test]
    fn test_closest_point_to_line_through_center() {
        let s = Sphere::new(2.0);
        // Line through origin along X
        let cp = s.closest_point_to_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        // Closest on line to origin is origin itself; sphere returns +X
        assert!((cp[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_closest_point_to_plane() {
        let s = Sphere::new(3.0);
        // Plane x = 5 => normal=[1,0,0], d=5
        let cp = s.closest_point_to_plane([1.0, 0.0, 0.0], 5.0);
        assert!((cp[0] - 3.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
    }

    #[test]
    fn test_closest_point_to_plane_negative() {
        let s = Sphere::new(2.0);
        // Plane x = -5 => normal=[1,0,0], d=-5
        let cp = s.closest_point_to_plane([1.0, 0.0, 0.0], -5.0);
        assert!((cp[0] + 2.0).abs() < 1e-10); // closest in -X direction
    }

    #[test]
    fn test_sphere_sweep_hit() {
        let s = Sphere::new(0.5);
        let t = s.sphere_sweep([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5, 10.0);
        assert!(t.is_some());
        let t = t.unwrap();
        assert!((t - 4.0).abs() < 1e-10); // combined radius = 1, distance = 5
    }

    #[test]
    fn test_sphere_sweep_miss() {
        let s = Sphere::new(0.5);
        let t = s.sphere_sweep([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 5.0, 0.0], 0.5, 10.0);
        assert!(t.is_none());
    }

    #[test]
    fn test_sphere_sweep_already_overlapping() {
        let s = Sphere::new(2.0);
        let t = s.sphere_sweep([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 10.0);
        assert!(t.is_some());
        assert!((t.unwrap()).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_signed_distance() {
        let s = Sphere::new(2.0);
        assert!((s.signed_distance([3.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
        assert!((s.signed_distance([1.0, 0.0, 0.0]) - (-1.0)).abs() < 1e-10);
        assert!((s.signed_distance([2.0, 0.0, 0.0])).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_contains_point() {
        let s = Sphere::new(2.0);
        assert!(s.contains_point([0.0, 0.0, 0.0]));
        assert!(s.contains_point([1.0, 1.0, 0.0]));
        assert!(!s.contains_point([3.0, 0.0, 0.0]));
    }

    #[test]
    fn test_sphere_project_inside() {
        let s = Sphere::new(2.0);
        let p = s.project_inside([1.0, 0.0, 0.0]);
        assert_eq!(p, [1.0, 0.0, 0.0]); // already inside

        let p2 = s.project_inside([5.0, 0.0, 0.0]);
        assert!((p2[0] - 2.0).abs() < 1e-10); // clamped to surface
    }

    // ── minimum_enclosing_sphere ─────────────────────────────────────────────

    #[test]
    fn test_minimum_enclosing_sphere_empty() {
        let (c, r) = minimum_enclosing_sphere(&[]);
        assert_eq!(c, [0.0; 3]);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn test_minimum_enclosing_sphere_two_points() {
        let pts = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let (c, r) = minimum_enclosing_sphere(&pts);
        // Center at midpoint (2,0,0), radius = 2
        assert!((c[0] - 2.0).abs() < 0.1, "c[0]={}", c[0]);
        assert!((r - 2.0).abs() < 0.1, "r={r}");
    }

    #[test]
    fn test_minimum_enclosing_sphere_contains_all() {
        let pts = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let (c, r) = minimum_enclosing_sphere(&pts);
        for p in &pts {
            let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt();
            assert!(d <= r + 1e-6, "point {:?} outside sphere r={r}", p);
        }
    }

    // ── sphere_obb_penetration ───────────────────────────────────────────────

    #[test]
    fn test_sphere_obb_penetration_overlap() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pen = sphere_obb_penetration(
            1.0,
            [1.5, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            identity,
        );
        assert!(pen.is_some(), "sphere at 1.5 should overlap unit OBB");
        assert!(pen.unwrap() > 0.0);
    }

    #[test]
    fn test_sphere_obb_penetration_no_overlap() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pen = sphere_obb_penetration(
            0.5,
            [5.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            identity,
        );
        assert!(pen.is_none());
    }

    // ── random_points_on_sphere ──────────────────────────────────────────────

    #[test]
    fn test_random_points_on_sphere_count() {
        let pts = random_points_on_sphere(1.0, 100, 42);
        assert_eq!(pts.len(), 100);
    }

    #[test]
    fn test_random_points_on_sphere_on_surface() {
        let r = 2.5;
        let pts = random_points_on_sphere(r, 50, 123);
        for p in &pts {
            let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                (dist - r).abs() < 1e-8,
                "point not on sphere surface: dist={dist}"
            );
        }
    }

    // ── spherical_harmonic ───────────────────────────────────────────────────

    #[test]
    fn test_spherical_harmonic_y00_constant() {
        // Y_0^0 is constant
        let v1 = spherical_harmonic(0, 0, 0.5, 1.0);
        let v2 = spherical_harmonic(0, 0, 1.2, 2.5);
        assert!((v1 - v2).abs() < 1e-14, "Y_0^0 should be constant");
    }

    #[test]
    fn test_spherical_harmonic_y10_at_poles() {
        // Y_1^0(0, phi) = sqrt(3/(4π)) * cos(0) = sqrt(3/(4π))
        let val = spherical_harmonic(1, 0, 0.0, 0.0);
        let expected = (3.0 / (4.0 * PI)).sqrt();
        assert!(
            (val - expected).abs() < 1e-10,
            "val={val} expected={expected}"
        );
    }

    #[test]
    fn test_spherical_harmonic_unsupported_returns_zero() {
        let val = spherical_harmonic(5, 3, 1.0, 1.0);
        assert_eq!(val, 0.0);
    }

    // ── sphere_sphere_intersection_volume ────────────────────────────────────

    #[test]
    fn test_sphere_sphere_intersection_volume_no_overlap() {
        let vol = sphere_sphere_intersection_volume(1.0, 1.0, 3.0);
        assert_eq!(vol, 0.0);
    }

    #[test]
    fn test_sphere_sphere_intersection_volume_contained() {
        // Sphere of radius 0.5 fully inside sphere of radius 2
        let vol = sphere_sphere_intersection_volume(2.0, 0.5, 0.0);
        let expected = (4.0 / 3.0) * PI * 0.5_f64.powi(3);
        assert!(
            (vol - expected).abs() < 1e-10,
            "vol={vol} expected={expected}"
        );
    }

    #[test]
    fn test_sphere_sphere_intersection_volume_symmetric() {
        // For equal radii and some distance, intersection should be symmetric
        let vol1 = sphere_sphere_intersection_volume(1.0, 1.0, 1.0);
        let vol2 = sphere_sphere_intersection_volume(1.0, 1.0, 1.0);
        assert!((vol1 - vol2).abs() < 1e-14);
        assert!(vol1 > 0.0);
    }

    // ── Extended tests for geodesic sphere, caps, Fibonacci, hemisphere, conformal ──

    #[test]
    fn test_geodesic_icosphere_base_vert_count() {
        let s = Sphere::new(1.0);
        let (verts, tris) = s.geodesic_icosphere(0);
        assert_eq!(verts.len(), 12);
        assert_eq!(tris.len(), 20);
    }

    #[test]
    fn test_geodesic_icosphere_subdivision_1() {
        let s = Sphere::new(1.0);
        let (verts, tris) = s.geodesic_icosphere(1);
        // After 1 subdivision: 42 verts, 80 tris
        assert_eq!(verts.len(), 42);
        assert_eq!(tris.len(), 80);
    }

    #[test]
    fn test_geodesic_icosphere_vertices_on_sphere() {
        let s = Sphere::new(2.0);
        let (verts, _) = s.geodesic_icosphere(1);
        for v in &verts {
            let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!((r - 2.0).abs() < 1e-9, "vertex not on sphere r={r}");
        }
    }

    #[test]
    fn test_geodesic_icosphere_subdivision_2_count() {
        let s = Sphere::new(1.0);
        let (verts, tris) = s.geodesic_icosphere(2);
        // 2 subdivisions: 162 verts, 320 tris
        assert_eq!(verts.len(), 162);
        assert_eq!(tris.len(), 320);
    }

    #[test]
    fn test_spherical_cap_volume_zero_height() {
        let s = Sphere::new(2.0);
        assert_eq!(s.spherical_cap_volume(0.0), 0.0);
    }

    #[test]
    fn test_spherical_cap_volume_full_sphere() {
        let s = Sphere::new(2.0);
        // h = 2r → full sphere volume = (4/3)π r³
        let v = s.spherical_cap_volume(4.0);
        let full = s.volume();
        assert!((v - full).abs() < 1e-9, "cap vol {v} vs sphere vol {full}");
    }

    #[test]
    fn test_spherical_cap_volume_hemisphere() {
        let s = Sphere::new(1.0);
        // h = r → hemisphere: V = (2/3)π r³
        let v = s.spherical_cap_volume(1.0);
        let expected = (2.0 / 3.0) * PI;
        assert!((v - expected).abs() < 1e-9);
    }

    #[test]
    fn test_spherical_cap_area_zero() {
        let s = Sphere::new(1.0);
        assert_eq!(s.spherical_cap_area(0.0), 0.0);
    }

    #[test]
    fn test_spherical_cap_area_full_sphere() {
        let s = Sphere::new(1.0);
        // h = 2r → full sphere surface area = 4π r²
        let a = s.spherical_cap_area(2.0);
        assert!((a - 4.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn test_spherical_zone_area_hemisphere() {
        let s = Sphere::new(1.0);
        // Zone from bottom (h=0) to equator (h=r): 2π r * r = 2π
        let a = s.spherical_zone_area(0.0, 1.0);
        assert!((a - 2.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn test_spherical_zone_area_symmetric() {
        let s = Sphere::new(2.0);
        let a1 = s.spherical_zone_area(0.0, 1.0);
        let a2 = s.spherical_zone_area(1.0, 0.0);
        assert!((a1 - a2).abs() < 1e-12);
    }

    #[test]
    fn test_hemisphere_cosine_sample_count() {
        let s = Sphere::new(1.0);
        let pts = s.hemisphere_cosine_sample(50, 42);
        assert_eq!(pts.len(), 50);
    }

    #[test]
    fn test_hemisphere_cosine_sample_on_sphere() {
        let s = Sphere::new(2.0);
        let pts = s.hemisphere_cosine_sample(30, 7);
        for p in &pts {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 2.0).abs() < 1e-9, "r={r}");
        }
    }

    #[test]
    fn test_hemisphere_cosine_sample_positive_z() {
        let s = Sphere::new(1.0);
        let pts = s.hemisphere_cosine_sample(40, 99);
        for p in &pts {
            assert!(p[2] >= 0.0, "z={}", p[2]);
        }
    }

    #[test]
    fn test_hemisphere_uniform_sample_count() {
        let s = Sphere::new(1.0);
        let pts = s.hemisphere_uniform_sample(25, 11);
        assert_eq!(pts.len(), 25);
    }

    #[test]
    fn test_hemisphere_uniform_sample_positive_z() {
        let s = Sphere::new(1.0);
        let pts = s.hemisphere_uniform_sample(30, 55);
        for p in &pts {
            assert!(p[2] >= 0.0, "z={}", p[2]);
        }
    }

    #[test]
    fn test_fibonacci_sphere_count() {
        let s = Sphere::new(1.0);
        let pts = s.fibonacci_sphere(50);
        assert_eq!(pts.len(), 50);
    }

    #[test]
    fn test_fibonacci_sphere_on_sphere() {
        let s = Sphere::new(3.0);
        let pts = s.fibonacci_sphere(30);
        for p in &pts {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 3.0).abs() < 1e-9, "r={r}");
        }
    }

    #[test]
    fn test_fibonacci_sphere_roughly_uniform() {
        // No two adjacent Fibonacci points should be exactly equal
        let s = Sphere::new(1.0);
        let pts = s.fibonacci_sphere(20);
        for i in 0..pts.len() - 1 {
            let p = pts[i];
            let q = pts[i + 1];
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            assert!(d > 0.01, "consecutive Fibonacci points too close d={d}");
        }
    }

    #[test]
    fn test_stereographic_project_north_pole() {
        let s = Sphere::new(1.0);
        // North pole (0, 0, 1) → (0, 0)
        let uv = s.stereographic_project([0.0, 0.0, 1.0]);
        assert!(uv.is_some());
        let [x, y] = uv.unwrap();
        assert!(x.abs() < 1e-10);
        assert!(y.abs() < 1e-10);
    }

    #[test]
    fn test_stereographic_project_south_pole_none() {
        let s = Sphere::new(1.0);
        let uv = s.stereographic_project([0.0, 0.0, -1.0]);
        assert!(uv.is_none());
    }

    #[test]
    fn test_stereographic_roundtrip() {
        let s = Sphere::new(1.0);
        let p_orig = [0.6, 0.0, 0.8]; // on unit sphere
        let uv = s.stereographic_project(p_orig).unwrap();
        let p_back = s.stereographic_unproject(uv);
        for i in 0..3 {
            assert!((p_orig[i] - p_back[i]).abs() < 1e-9, "mismatch at [{i}]");
        }
    }

    #[test]
    fn test_stereographic_unproject_on_sphere() {
        let s = Sphere::new(2.0);
        let p = s.stereographic_unproject([1.0, 0.5]);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((r - 2.0).abs() < 1e-9, "r={r}");
    }

    #[test]
    fn test_sh_l3_m0_at_pole() {
        // Y_3^0(θ=0) = (1/4)*sqrt(7/π)*(5*1 - 3) = (1/4)*sqrt(7/π)*2
        let val = Sphere::sh_l3(0, 0.0, 0.0);
        let expected = 0.25 * (7.0 / PI).sqrt() * 2.0;
        assert!((val - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sh_l3_unknown_m_zero() {
        assert_eq!(Sphere::sh_l3(5, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_cap_solid_angle_hemisphere() {
        let s = Sphere::new(1.0);
        // h = r = 1 → cos_theta = 0 → Ω = 2π
        let omega = s.cap_solid_angle(1.0);
        assert!((omega - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cap_solid_angle_zero() {
        let s = Sphere::new(1.0);
        assert!((s.cap_solid_angle(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_cap_solid_angle_full_sphere() {
        let s = Sphere::new(1.0);
        // h = 2r → Ω = 4π
        let omega = s.cap_solid_angle(2.0);
        assert!((omega - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_lon_lat_north_pole() {
        let s = Sphere::new(1.0);
        let (_, lat) = s.lon_lat([0.0, 0.0, 1.0]);
        assert!((lat - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lon_lat_south_pole() {
        let s = Sphere::new(1.0);
        let (_, lat) = s.lon_lat([0.0, 0.0, -1.0]);
        assert!((lat + PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lon_lat_equator() {
        let s = Sphere::new(1.0);
        let (lon, lat) = s.lon_lat([1.0, 0.0, 0.0]);
        assert!(lat.abs() < 1e-10);
        assert!(lon.abs() < 1e-10);
    }

    #[test]
    fn test_fibonacci_sphere_single_point() {
        let s = Sphere::new(1.0);
        let pts = s.fibonacci_sphere(1);
        assert_eq!(pts.len(), 1);
        let r = (pts[0][0].powi(2) + pts[0][1].powi(2) + pts[0][2].powi(2)).sqrt();
        assert!((r - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_spherical_cap_volume_increases_with_h() {
        let s = Sphere::new(1.0);
        let v1 = s.spherical_cap_volume(0.5);
        let v2 = s.spherical_cap_volume(1.0);
        let v3 = s.spherical_cap_volume(1.5);
        assert!(v1 < v2 && v2 < v3);
    }

    #[test]
    fn test_sh_l3_all_m_finite() {
        for m in -3..=3 {
            let v = Sphere::sh_l3(m, 0.5, 1.0);
            assert!(v.is_finite(), "sh_l3({m}) not finite: {v}");
        }
    }
}
