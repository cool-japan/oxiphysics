// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cylinder shape (Y-axis aligned).

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::f64::consts::PI;

/// A cylinder defined by radius and half-height along the Y axis.
#[derive(Debug, Clone)]
pub struct Cylinder {
    /// Radius of the cylinder.
    pub radius: Real,
    /// Half-height along the Y axis.
    pub half_height: Real,
}

impl Cylinder {
    /// Create a new cylinder.
    pub fn new(radius: Real, half_height: Real) -> Self {
        Self {
            radius,
            half_height,
        }
    }

    /// Surface area: 2πr² + 2πrh (two caps + lateral surface).
    pub fn surface_area(&self) -> Real {
        let h = 2.0 * self.half_height;
        2.0 * PI * self.radius * self.radius + 2.0 * PI * self.radius * h
    }

    /// Volume: πr²h (full height = 2 * half_height).
    pub fn volume_explicit(&self) -> Real {
        PI * self.radius * self.radius * 2.0 * self.half_height
    }

    /// Inertia tensor as \[\[f64;3\\];3] row-major.
    /// Ixx = Izz = m(3r²+h²)/12, Iyy = mr²/2 (Y is symmetry axis).
    pub fn inertia_tensor_array(&self, mass: f64) -> [[f64; 3]; 3] {
        let r2 = self.radius * self.radius;
        let h = 2.0 * self.half_height;
        let h2 = h * h;
        let iy = 0.5 * mass * r2;
        let ixz = mass * (3.0 * r2 + h2) / 12.0;
        [[ixz, 0.0, 0.0], [0.0, iy, 0.0], [0.0, 0.0, ixz]]
    }

    /// Ray cast returning (t, normal) as plain arrays.
    /// Returns `None` if no intersection within `max_toi`.
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

    /// Closest point on (or inside) the cylinder to point `p`.
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        let px = p[0];
        let py = p[1];
        let pz = p[2];

        // Clamp Y to cylinder height range
        let cy = py.clamp(-self.half_height, self.half_height);

        // Project onto XZ and clamp to disk radius
        let xz_len = (px * px + pz * pz).sqrt();
        let (cx, cz) = if xz_len <= self.radius {
            (px, pz)
        } else {
            let s = self.radius / xz_len;
            (px * s, pz * s)
        };

        // Now find whether the closest point is on the cap or side
        let inside_xz = xz_len <= self.radius;
        let inside_y = py.abs() <= self.half_height;

        if inside_xz && inside_y {
            // Point is inside — find nearest surface
            let dist_to_top = self.half_height - py;
            let dist_to_bot = py + self.half_height;
            let dist_to_side = self.radius - xz_len;
            let min_d = dist_to_top.min(dist_to_bot).min(dist_to_side);

            if min_d == dist_to_side {
                let s = self.radius / xz_len.max(1e-30);
                [px * s, py, pz * s]
            } else if dist_to_top <= dist_to_bot {
                [px, self.half_height, pz]
            } else {
                [px, -self.half_height, pz]
            }
        } else {
            [cx, cy, cz]
        }
    }

    /// Returns true if point `p` is strictly inside the cylinder.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let xz2 = p[0] * p[0] + p[2] * p[2];
        xz2 <= self.radius * self.radius && p[1].abs() <= self.half_height
    }

    // ── New methods ──

    /// GJK support function using plain arrays.
    /// Returns the farthest point on the cylinder in the given direction.
    pub fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let xz_len = (direction[0] * direction[0] + direction[2] * direction[2]).sqrt();
        let (sx, sz) = if xz_len > 1e-10 {
            (
                self.radius * direction[0] / xz_len,
                self.radius * direction[2] / xz_len,
            )
        } else {
            (self.radius, 0.0)
        };
        let sy = self.half_height.copysign(direction[1]);
        [sx, sy, sz]
    }

    /// Signed distance from a point to the cylinder surface.
    /// Negative if inside, positive if outside.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        let xz_len = (p[0] * p[0] + p[2] * p[2]).sqrt();
        let dist_side = xz_len - self.radius;
        let dist_cap = p[1].abs() - self.half_height;

        if dist_side <= 0.0 && dist_cap <= 0.0 {
            // Inside: return the larger (less negative) distance
            dist_side.max(dist_cap)
        } else if dist_side > 0.0 && dist_cap > 0.0 {
            // Outside corner: Euclidean distance to edge
            (dist_side * dist_side + dist_cap * dist_cap).sqrt()
        } else {
            // Outside one dimension, inside the other
            dist_side.max(dist_cap)
        }
    }

    /// Closest point on the cylinder to a line segment from `a` to `b`.
    /// Returns the closest point on the cylinder surface.
    pub fn closest_point_to_segment(&self, a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        // Sample several points on the segment and find the one whose
        // closest-point on cylinder is nearest.
        let n_samples = 16;
        let mut best_dist_sq = f64::INFINITY;
        let mut best_cp = [0.0; 3];

        for i in 0..=n_samples {
            let t = i as f64 / n_samples as f64;
            let seg_pt = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            let cp = self.closest_point(seg_pt);
            let dx = cp[0] - seg_pt[0];
            let dy = cp[1] - seg_pt[1];
            let dz = cp[2] - seg_pt[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_dist_sq {
                best_dist_sq = d2;
                best_cp = cp;
            }
        }
        best_cp
    }

    /// Ray cast against only the top cap (y = +half_height).
    /// Returns `Some((t, normal))` if hit.
    pub fn ray_cast_top_cap(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        if direction[1].abs() < 1e-12 {
            return None;
        }
        let t = (self.half_height - origin[1]) / direction[1];
        if t < 0.0 || t > max_toi {
            return None;
        }
        let px = origin[0] + t * direction[0];
        let pz = origin[2] + t * direction[2];
        if px * px + pz * pz <= self.radius * self.radius {
            Some((t, [0.0, 1.0, 0.0]))
        } else {
            None
        }
    }

    /// Ray cast against only the bottom cap (y = -half_height).
    /// Returns `Some((t, normal))` if hit.
    pub fn ray_cast_bottom_cap(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        if direction[1].abs() < 1e-12 {
            return None;
        }
        let t = (-self.half_height - origin[1]) / direction[1];
        if t < 0.0 || t > max_toi {
            return None;
        }
        let px = origin[0] + t * direction[0];
        let pz = origin[2] + t * direction[2];
        if px * px + pz * pz <= self.radius * self.radius {
            Some((t, [0.0, -1.0, 0.0]))
        } else {
            None
        }
    }

    /// Ray cast against only the lateral (curved) surface.
    /// Returns `Some((t, normal))` if hit.
    pub fn ray_cast_lateral(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        let a = direction[0] * direction[0] + direction[2] * direction[2];
        if a < 1e-12 {
            return None;
        }
        let b = 2.0 * (origin[0] * direction[0] + origin[2] * direction[2]);
        let c = origin[0] * origin[0] + origin[2] * origin[2] - self.radius * self.radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        for sign in [-1.0, 1.0] {
            let t = (-b + sign * sqrt_disc) / (2.0 * a);
            if t >= 0.0 && t <= max_toi {
                let py = origin[1] + t * direction[1];
                if py.abs() <= self.half_height {
                    let px = origin[0] + t * direction[0];
                    let pz = origin[2] + t * direction[2];
                    let xz_len = (px * px + pz * pz).sqrt();
                    let nx = if xz_len > 1e-12 { px / xz_len } else { 1.0 };
                    let nz = if xz_len > 1e-12 { pz / xz_len } else { 0.0 };
                    return Some((t, [nx, 0.0, nz]));
                }
            }
        }
        None
    }

    /// Compute the lateral surface area (excluding caps).
    pub fn lateral_surface_area(&self) -> f64 {
        2.0 * PI * self.radius * 2.0 * self.half_height
    }

    /// Compute the cap area (one cap).
    pub fn cap_area(&self) -> f64 {
        PI * self.radius * self.radius
    }

    /// Project a point onto the cylinder axis (the Y-axis segment).
    /// Returns the clamped Y coordinate.
    pub fn project_on_axis(&self, p: [f64; 3]) -> f64 {
        p[1].clamp(-self.half_height, self.half_height)
    }

    /// Distance from a point to the cylinder axis (the Y axis segment).
    pub fn distance_to_axis(&self, p: [f64; 3]) -> f64 {
        let clamped_y = p[1].clamp(-self.half_height, self.half_height);
        let dy = p[1] - clamped_y;
        let xz2 = p[0] * p[0] + p[2] * p[2];
        (xz2 + dy * dy).sqrt()
    }

    // ── New expanded methods ──

    /// Cylinder vs infinite plane intersection test.
    ///
    /// Plane: `dot(plane_normal, x) = plane_d`.
    /// Returns `true` if any part of the finite cylinder intersects the plane.
    pub fn intersects_plane(&self, plane_normal: [f64; 3], plane_d: f64) -> bool {
        // The most positive and negative signed distances are attained at the
        // four extreme points: (±r_proj, ±hh, 0) projected onto the plane normal.
        // XZ projection contribution: r * sqrt(nx² + nz²)
        let xz_proj =
            (plane_normal[0] * plane_normal[0] + plane_normal[2] * plane_normal[2]).sqrt();
        let y_contrib_top = plane_normal[1] * self.half_height;
        let y_contrib_bot = -plane_normal[1] * self.half_height;
        let r_contrib = self.radius * xz_proj;

        let max_sd = y_contrib_top.max(y_contrib_bot) + r_contrib;
        let min_sd = y_contrib_top.min(y_contrib_bot) - r_contrib;

        min_sd <= plane_d && plane_d <= max_sd
    }

    /// Cylinder SDF (Signed Distance Function).
    ///
    /// Identical to `signed_distance` but exposed under the SDF naming convention.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        self.signed_distance(p)
    }

    /// Cylinder-sphere intersection test.
    ///
    /// Returns `true` if the sphere (center + radius) overlaps the finite cylinder.
    pub fn intersects_sphere(&self, sphere_center: [f64; 3], sphere_radius: f64) -> bool {
        let cp = self.closest_point(sphere_center);
        let dx = sphere_center[0] - cp[0];
        let dy = sphere_center[1] - cp[1];
        let dz = sphere_center[2] - cp[2];
        dx * dx + dy * dy + dz * dz <= sphere_radius * sphere_radius
    }

    /// Test a finite cylinder against an infinite cylinder (same Y axis, different radius).
    ///
    /// An infinite cylinder has the same Y axis but extends infinitely in ±Y.
    /// They intersect if the radii overlap (i.e., `|r_finite - r_infinite| ≤ lateral dist`).
    /// Since both are Y-aligned and centered at origin, this reduces to a 2-D circle test.
    pub fn intersects_infinite_cylinder(&self, infinite_radius: f64, center_xz: [f64; 2]) -> bool {
        // Distance between the two Y-axis lines in XZ
        let dist = (center_xz[0] * center_xz[0] + center_xz[1] * center_xz[1]).sqrt();
        // Overlap when distance < sum of radii
        dist < self.radius + infinite_radius
    }

    /// Generate `n` uniformly distributed random points on the cylinder surface
    /// using a deterministic xorshift64 PRNG.
    pub fn random_surface_points(&self, n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut points = Vec::with_capacity(n);
        let lat = self.lateral_surface_area();
        let cap = self.cap_area();
        let total = lat + 2.0 * cap;
        let p_lat = lat / total;
        let p_top = cap / total;

        let mut state = seed;
        let next = |s: &mut u64| -> f64 {
            *s ^= *s << 13;
            *s ^= *s >> 7;
            *s ^= *s << 17;
            (*s as f64) / (u64::MAX as f64)
        };

        while points.len() < n {
            let u = next(&mut state);
            let v = next(&mut state);
            let w = next(&mut state);

            if u < p_lat {
                let theta = v * 2.0 * PI;
                let y = (w * 2.0 - 1.0) * self.half_height;
                points.push([self.radius * theta.cos(), y, self.radius * theta.sin()]);
            } else if u < p_lat + p_top {
                // Top cap: uniform disk
                let r = self.radius * v.sqrt();
                let theta = w * 2.0 * PI;
                points.push([r * theta.cos(), self.half_height, r * theta.sin()]);
            } else {
                // Bottom cap
                let r = self.radius * v.sqrt();
                let theta = w * 2.0 * PI;
                points.push([r * theta.cos(), -self.half_height, r * theta.sin()]);
            }
        }
        points
    }

    /// Cylinder support function (plain array, same as `support` but for trait-like usage).
    pub fn support_array(&self, direction: [f64; 3]) -> [f64; 3] {
        self.support(direction)
    }

    /// Full height of the cylinder (2 × half_height).
    pub fn full_height(&self) -> f64 {
        2.0 * self.half_height
    }

    /// Radius of gyration about the Y (symmetry) axis.
    ///
    /// k = sqrt(Iy / m) = sqrt(r²/2) = r / sqrt(2).
    pub fn radius_of_gyration_y(&self) -> f64 {
        self.radius / 2.0_f64.sqrt()
    }

    /// Radius of gyration about the X (transverse) axis.
    ///
    /// k = sqrt(Ix / m) = sqrt((3r² + h²) / 12).
    pub fn radius_of_gyration_x(&self, _mass: f64) -> f64 {
        let r2 = self.radius * self.radius;
        let h = 2.0 * self.half_height;
        let h2 = h * h;
        ((3.0 * r2 + h2) / 12.0).sqrt()
    }

    // ── Extended geometry: hollow, frustum, oblique, ruled surface, UV cap ──

    /// Volume of the swept region when the cylinder moves by translation `delta`.
    ///
    /// For a cylinder swept along a vector `delta` (not along its own axis),
    /// this returns the approximate swept volume: V_cyl + A_cap * |delta|.
    pub fn volume_swept(&self, delta: [f64; 3]) -> f64 {
        let dist = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
        self.volume() + self.cap_area() * dist
    }

    /// UV mapping for the cylinder cap at `+half_height`.
    ///
    /// Maps a 3-D point `p` on the cap disk to UV in `[0, 1]²`.
    /// Uses a radial projection: `u = 0.5 + x/(2r)`, `v = 0.5 + z/(2r)`.
    pub fn cap_uv(&self, p: [f64; 3]) -> [f64; 2] {
        let u = 0.5 + p[0] / (2.0 * self.radius);
        let v = 0.5 + p[2] / (2.0 * self.radius);
        [u.clamp(0.0, 1.0), v.clamp(0.0, 1.0)]
    }

    /// UV mapping for the lateral (rolled-out) surface.
    ///
    /// `u` is the normalised circumferential angle (azimuth), `v` is the
    /// normalised height. Both in `[0, 1]`.
    pub fn lateral_uv(&self, p: [f64; 3]) -> [f64; 2] {
        let theta = p[2].atan2(p[0]); // in (-π, π]
        let u = ((theta / (2.0 * PI)) + 1.0) % 1.0;
        let v = (p[1] + self.half_height) / (2.0 * self.half_height);
        [u, v.clamp(0.0, 1.0)]
    }

    /// Hollow cylinder check: ring-shaped cross-section.
    ///
    /// Returns `true` if `p` is inside the hollow cylinder with inner radius
    /// `inner_r` (a tube with outer radius = `self.radius`).
    pub fn contains_point_hollow(&self, p: [f64; 3], inner_r: f64) -> bool {
        let xz2 = p[0] * p[0] + p[2] * p[2];
        let r_inner = inner_r.max(0.0).min(self.radius);
        xz2 <= self.radius * self.radius
            && xz2 >= r_inner * r_inner
            && p[1].abs() <= self.half_height
    }

    /// Volume of the hollow cylinder (annular cylinder).
    ///
    /// Volume = π(R² - r²) * h where R = outer, r = inner, h = full height.
    pub fn volume_hollow(&self, inner_radius: f64) -> f64 {
        let r = inner_radius.max(0.0).min(self.radius);
        PI * (self.radius * self.radius - r * r) * 2.0 * self.half_height
    }

    /// Frustum (truncated cone) slant height given top and bottom radii.
    ///
    /// `r_top` and `r_bottom` are radii at `+half_height` and `-half_height`.
    /// Returns the slant height of the frustum.
    pub fn frustum_slant_height(&self, r_top: f64, r_bottom: f64) -> f64 {
        let h = 2.0 * self.half_height;
        let dr = (r_top - r_bottom).abs();
        (h * h + dr * dr).sqrt()
    }

    /// Lateral area of a frustum with given top and bottom radii.
    ///
    /// Lateral area = π * (r_top + r_bottom) * slant_height.
    pub fn frustum_lateral_area(&self, r_top: f64, r_bottom: f64) -> f64 {
        let slant = self.frustum_slant_height(r_top, r_bottom);
        PI * (r_top + r_bottom) * slant
    }

    /// Volume of a frustum with given top and bottom radii.
    ///
    /// V = π h (r_top² + r_top*r_bottom + r_bottom²) / 3.
    pub fn frustum_volume(&self, r_top: f64, r_bottom: f64) -> f64 {
        let h = 2.0 * self.half_height;
        PI * h * (r_top * r_top + r_top * r_bottom + r_bottom * r_bottom) / 3.0
    }

    /// Oblique cylinder end-cap area at a given tilt.
    ///
    /// When the top cap is tilted at `tilt_angle` radians from horizontal,
    /// the effective cap area becomes `π r² / cos(tilt_angle)`.
    pub fn oblique_cap_area(&self, tilt_angle: f64) -> f64 {
        let cos_a = tilt_angle.cos().abs().max(1e-12);
        PI * self.radius * self.radius / cos_a
    }

    /// Ruled surface interpolation: return a point on the ruled surface
    /// between two guide curves (bottom circle at `-half_height` and top
    /// circle at `+half_height`).
    ///
    /// `u ∈ [0, 2π)` is the azimuth, `t ∈ [0, 1\]` is the height parameter.
    /// For a straight cylinder both circles have the same radius, so the
    /// ruled surface is just the lateral surface.
    pub fn ruled_surface_point(&self, u: f64, t: f64) -> [f64; 3] {
        let y = -self.half_height + t * 2.0 * self.half_height;
        [self.radius * u.cos(), y, self.radius * u.sin()]
    }

    /// Stack two cylinders and return the combined volume.
    ///
    /// If the second cylinder has the same radius but different height it is
    /// stacked end-to-end on top of this cylinder.  Returns the sum of their
    /// volumes.
    pub fn stacked_volume(&self, other: &Cylinder) -> f64 {
        self.volume() + other.volume()
    }

    /// Stack two cylinders and return the combined bounding half-height.
    ///
    /// Returns the total half-height of the axis-aligned bounding box of the
    /// two cylinders stacked along the Y axis.
    pub fn stacked_half_height(&self, other: &Cylinder) -> f64 {
        self.half_height + other.half_height
    }

    /// Combined bounding radius when two cylinders are stacked.
    ///
    /// Returns the maximum of the two radii.
    pub fn stacked_bounding_radius(&self, other: &Cylinder) -> f64 {
        self.radius.max(other.radius)
    }

    /// Surface normal at a lateral point `p` (unit outward normal in XZ).
    ///
    /// Ignores the Y component; suitable for points on the curved surface.
    pub fn lateral_normal_at(&self, p: [f64; 3]) -> [f64; 3] {
        let len = (p[0] * p[0] + p[2] * p[2]).sqrt();
        if len < 1e-12 {
            return [1.0, 0.0, 0.0];
        }
        [p[0] / len, 0.0, p[2] / len]
    }

    /// Cylinder SDF gradient (approximate unit normal direction).
    ///
    /// Returns the gradient of the SDF at `p`, i.e. the direction of
    /// fastest increase.
    pub fn sdf_gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let eps = 1e-5;
        let sdf0 = self.sdf(p);
        let gx = (self.sdf([p[0] + eps, p[1], p[2]]) - sdf0) / eps;
        let gy = (self.sdf([p[0], p[1] + eps, p[2]]) - sdf0) / eps;
        let gz = (self.sdf([p[0], p[1], p[2] + eps]) - sdf0) / eps;
        let len = (gx * gx + gy * gy + gz * gz).sqrt().max(1e-14);
        [gx / len, gy / len, gz / len]
    }

    /// Count the number of times a ray intersects the infinite cylinder
    /// (lateral surface only, ignoring caps).
    ///
    /// Returns 0, 1, or 2.
    pub fn lateral_ray_intersection_count(&self, origin: [f64; 3], direction: [f64; 3]) -> usize {
        let a = direction[0] * direction[0] + direction[2] * direction[2];
        if a < 1e-14 {
            return 0;
        }
        let b = 2.0 * (origin[0] * direction[0] + origin[2] * direction[2]);
        let c = origin[0] * origin[0] + origin[2] * origin[2] - self.radius * self.radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return 0;
        }
        if disc < 1e-14 {
            return 1;
        }
        2
    }

    /// Generate `n` points uniformly distributed along the top cap edge (rim).
    ///
    /// Returns points on the circle of radius `r` at `y = +half_height`.
    pub fn top_rim_points(&self, n: usize) -> Vec<[f64; 3]> {
        let n = n.max(2);
        (0..n)
            .map(|i| {
                let t = 2.0 * PI * i as f64 / n as f64;
                [
                    self.radius * t.cos(),
                    self.half_height,
                    self.radius * t.sin(),
                ]
            })
            .collect()
    }
}

impl Shape for Cylinder {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn bounding_box(&self) -> Aabb {
        Aabb::new(
            Vec3::new(-self.radius, -self.half_height, -self.radius),
            Vec3::new(self.radius, self.half_height, self.radius),
        )
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let xz_len = (direction.x * direction.x + direction.z * direction.z).sqrt();
        let (sx, sz) = if xz_len > 1e-10 {
            (
                self.radius * direction.x / xz_len,
                self.radius * direction.z / xz_len,
            )
        } else {
            (self.radius, 0.0)
        };
        let sy = self.half_height.copysign(direction.y);
        Vec3::new(sx, sy, sz)
    }

    fn volume(&self) -> Real {
        PI * self.radius * self.radius * 2.0 * self.half_height
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::zeros()
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let r2 = self.radius * self.radius;
        let h2 = (2.0 * self.half_height).powi(2);
        let iy = 0.5 * mass * r2;
        let ixz = mass * (3.0 * r2 + h2) / 12.0;
        Mat3::new(ixz, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ixz)
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;

        // Infinite cylinder in XZ
        let a = ray_direction.x.powi(2) + ray_direction.z.powi(2);
        let b = 2.0 * (ray_origin.x * ray_direction.x + ray_origin.z * ray_direction.z);
        let c = ray_origin.x.powi(2) + ray_origin.z.powi(2) - self.radius.powi(2);

        if a > 1e-12 {
            let disc = b * b - 4.0 * a * c;
            if disc >= 0.0 {
                let sqrt_disc = disc.sqrt();
                for sign in [-1.0, 1.0] {
                    let t = (-b + sign * sqrt_disc) / (2.0 * a);
                    if t >= 0.0 && t <= max_toi {
                        let p = ray_origin + ray_direction * t;
                        if p.y.abs() <= self.half_height {
                            let normal = Vec3::new(p.x, 0.0, p.z).normalize();
                            if best.as_ref().is_none_or(|prev| t < prev.toi) {
                                best = Some(RayHit {
                                    point: p,
                                    normal,
                                    toi: t,
                                });
                            }
                        }
                    }
                }
            }
        }

        // Top and bottom caps
        if ray_direction.y.abs() > 1e-12 {
            for &cap_y in &[self.half_height, -self.half_height] {
                let t = (cap_y - ray_origin.y) / ray_direction.y;
                if t >= 0.0 && t <= max_toi {
                    let p = ray_origin + ray_direction * t;
                    if p.x.powi(2) + p.z.powi(2) <= self.radius.powi(2) {
                        let normal = Vec3::new(0.0, cap_y.signum(), 0.0);
                        if best.as_ref().is_none_or(|prev| t < prev.toi) {
                            best = Some(RayHit {
                                point: p,
                                normal,
                                toi: t,
                            });
                        }
                    }
                }
            }
        }

        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cylinder_volume() {
        let c = Cylinder::new(1.0, 1.0);
        let expected = PI * 2.0;
        assert!((c.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_surface_area() {
        // r=1, h=2 (half_height=1): 2π*1 + 2π*1*2 = 2π + 4π = 6π
        let c = Cylinder::new(1.0, 1.0);
        let expected = 6.0 * PI;
        assert!((c.surface_area() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_inertia_symmetry() {
        let c = Cylinder::new(1.5, 2.0);
        let it = c.inertia_tensor(3.0);
        // I_xx == I_zz (symmetry around Y axis)
        let diff = (it[(0, 0)] - it[(2, 2)]).abs();
        assert!(
            diff < 1e-10,
            "I_xx != I_zz: {} vs {}",
            it[(0, 0)],
            it[(2, 2)]
        );
        // Off-diagonal should be zero
        assert!(it[(0, 1)].abs() < 1e-10);
        assert!(it[(1, 0)].abs() < 1e-10);
        assert!(it[(0, 2)].abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_inertia_array() {
        let c = Cylinder::new(1.0, 0.5);
        let it = c.inertia_tensor_array(2.0);
        // Iyy = mr²/2 = 2*1/2 = 1.0
        assert!((it[1][1] - 1.0).abs() < 1e-10);
        // Ixx = m(3r²+h²)/12 = 2*(3+1)/12 = 8/12 = 2/3
        let h = 1.0_f64; // full height
        let ixx = 2.0 * (3.0 + h * h) / 12.0;
        assert!((it[0][0] - ixx).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_raycast() {
        let c = Cylinder::new(1.0, 2.0);
        let origin = Vec3::new(-5.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let hit = c.ray_cast(&origin, &dir, 100.0).unwrap();
        assert!((hit.toi - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_raycast_array() {
        let c = Cylinder::new(1.0, 2.0);
        let (t, n) = c
            .ray_cast_array([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0)
            .unwrap();
        assert!((t - 4.0).abs() < 1e-10);
        assert!((n[0] + 1.0).abs() < 1e-10); // normal points -X
    }

    #[test]
    fn test_cylinder_raycast_cap() {
        let c = Cylinder::new(1.0, 2.0);
        let origin = Vec3::new(0.0, 10.0, 0.0);
        let dir = Vec3::new(0.0, -1.0, 0.0);
        let hit = c.ray_cast(&origin, &dir, 100.0).unwrap();
        assert!((hit.toi - 8.0).abs() < 1e-10);
        assert!((hit.normal.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_contains_point() {
        let c = Cylinder::new(1.0, 2.0);
        assert!(c.contains_point([0.0, 0.0, 0.0]));
        assert!(c.contains_point([0.9, 1.9, 0.0]));
        assert!(!c.contains_point([0.0, 2.1, 0.0]));
        assert!(!c.contains_point([1.1, 0.0, 0.0]));
    }

    #[test]
    fn test_cylinder_closest_point_outside_side() {
        let c = Cylinder::new(1.0, 2.0);
        let cp = c.closest_point([3.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!((cp[1]).abs() < 1e-10);
        assert!((cp[2]).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_closest_point_above_cap() {
        let c = Cylinder::new(1.0, 2.0);
        let cp = c.closest_point([0.0, 5.0, 0.0]);
        assert!((cp[1] - 2.0).abs() < 1e-10);
    }

    // ── New tests ──

    #[test]
    fn test_cylinder_support_array() {
        let c = Cylinder::new(2.0, 3.0);
        let sp = c.support([1.0, 0.0, 0.0]);
        assert!((sp[0] - 2.0).abs() < 1e-10);
        assert!((sp[1] - 3.0).abs() < 1e-10); // positive Y for positive-ish Y (copysign(0) = +)
    }

    #[test]
    fn test_cylinder_support_negative_y() {
        let c = Cylinder::new(1.0, 2.0);
        let sp = c.support([0.0, -1.0, 0.0]);
        assert!((sp[1] + 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_support_diagonal() {
        let c = Cylinder::new(1.0, 2.0);
        let sp = c.support([1.0, 1.0, 0.0]);
        assert!((sp[0] - 1.0).abs() < 1e-10); // full radius in XZ direction
        assert!((sp[1] - 2.0).abs() < 1e-10); // top cap
    }

    #[test]
    fn test_cylinder_signed_distance_inside() {
        let c = Cylinder::new(2.0, 3.0);
        let d = c.signed_distance([0.0, 0.0, 0.0]);
        assert!(d < 0.0); // center is inside
        // Distance to nearest surface = min(2, 3) = 2 (radial)
        assert!((d + 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_signed_distance_outside_side() {
        let c = Cylinder::new(1.0, 2.0);
        let d = c.signed_distance([3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-10); // 3 - 1 = 2
    }

    #[test]
    fn test_cylinder_signed_distance_outside_cap() {
        let c = Cylinder::new(1.0, 2.0);
        let d = c.signed_distance([0.0, 4.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-10); // 4 - 2 = 2
    }

    #[test]
    fn test_cylinder_signed_distance_outside_corner() {
        let c = Cylinder::new(1.0, 1.0);
        // Outside both radially and vertically
        let d = c.signed_distance([2.0, 2.0, 0.0]);
        // dist_side=1, dist_cap=1, corner distance = sqrt(2)
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_ray_cast_top_cap() {
        let c = Cylinder::new(1.0, 2.0);
        let result = c.ray_cast_top_cap([0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, n) = result.unwrap();
        assert!((t - 3.0).abs() < 1e-10);
        assert!((n[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_ray_cast_top_cap_miss() {
        let c = Cylinder::new(1.0, 2.0);
        // Ray misses the cap disk
        let result = c.ray_cast_top_cap([3.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_cylinder_ray_cast_bottom_cap() {
        let c = Cylinder::new(1.0, 2.0);
        let result = c.ray_cast_bottom_cap([0.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, n) = result.unwrap();
        assert!((t - 3.0).abs() < 1e-10);
        assert!((n[1] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_ray_cast_lateral() {
        let c = Cylinder::new(1.0, 2.0);
        let result = c.ray_cast_lateral([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, n) = result.unwrap();
        assert!((t - 4.0).abs() < 1e-10);
        assert!((n[0] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_ray_cast_lateral_miss_height() {
        let c = Cylinder::new(1.0, 1.0);
        // Ray at y=5, which is above the cylinder
        let result = c.ray_cast_lateral([-5.0, 5.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_cylinder_lateral_surface_area() {
        let c = Cylinder::new(1.0, 1.0);
        // 2π * r * h = 2π * 1 * 2 = 4π
        assert!((c.lateral_surface_area() - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_cap_area() {
        let c = Cylinder::new(2.0, 1.0);
        assert!((c.cap_area() - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_surface_area_decomposition() {
        let c = Cylinder::new(1.5, 2.5);
        let total = c.surface_area();
        let decomposed = c.lateral_surface_area() + 2.0 * c.cap_area();
        assert!((total - decomposed).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_project_on_axis() {
        let c = Cylinder::new(1.0, 2.0);
        assert!((c.project_on_axis([0.0, 5.0, 0.0]) - 2.0).abs() < 1e-10);
        assert!((c.project_on_axis([0.0, -5.0, 0.0]) + 2.0).abs() < 1e-10);
        assert!((c.project_on_axis([3.0, 1.0, 0.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_distance_to_axis() {
        let c = Cylinder::new(1.0, 2.0);
        // Point on axis
        assert!(c.distance_to_axis([0.0, 0.0, 0.0]).abs() < 1e-10);
        // Point at (3,0,0), on the midplane
        assert!((c.distance_to_axis([3.0, 0.0, 0.0]) - 3.0).abs() < 1e-10);
        // Point above cap, on Y axis
        assert!((c.distance_to_axis([0.0, 5.0, 0.0]) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_closest_point_to_segment() {
        let c = Cylinder::new(1.0, 2.0);
        // Segment far away in X
        let cp = c.closest_point_to_segment([5.0, 0.0, 0.0], [5.0, 1.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_volume_explicit_matches() {
        let c = Cylinder::new(1.5, 2.5);
        assert!((c.volume_explicit() - c.volume()).abs() < 1e-10);
    }

    // ── Expanded tests ──

    #[test]
    fn test_cylinder_intersects_plane_through_center() {
        let c = Cylinder::new(1.0, 2.0);
        // Horizontal plane at y=0: normal=(0,1,0), d=0
        assert!(c.intersects_plane([0.0, 1.0, 0.0], 0.0));
    }

    #[test]
    fn test_cylinder_intersects_plane_above() {
        let c = Cylinder::new(1.0, 2.0);
        // Plane far above: normal=(0,1,0), d=10 → cylinder max y contribution = 2, so no hit
        assert!(!c.intersects_plane([0.0, 1.0, 0.0], 10.0));
    }

    #[test]
    fn test_cylinder_intersects_plane_lateral() {
        let c = Cylinder::new(2.0, 1.0);
        // Vertical plane at x=1: normal=(1,0,0), d=1 → radius=2, should intersect
        assert!(c.intersects_plane([1.0, 0.0, 0.0], 1.0));
        // Plane at x=3: outside radius
        assert!(!c.intersects_plane([1.0, 0.0, 0.0], 3.0));
    }

    #[test]
    fn test_cylinder_sdf_matches_signed_distance() {
        let c = Cylinder::new(1.0, 2.0);
        let p = [3.0, 0.0, 0.0];
        assert!((c.sdf(p) - c.signed_distance(p)).abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_intersects_sphere_yes() {
        let c = Cylinder::new(1.0, 2.0);
        assert!(c.intersects_sphere([0.0, 0.0, 0.0], 0.5));
        assert!(c.intersects_sphere([1.5, 0.0, 0.0], 0.6)); // slightly outside, sphere big enough
    }

    #[test]
    fn test_cylinder_intersects_sphere_no() {
        let c = Cylinder::new(1.0, 2.0);
        assert!(!c.intersects_sphere([5.0, 0.0, 0.0], 0.5));
    }

    #[test]
    fn test_cylinder_intersects_infinite_cylinder_overlap() {
        let c = Cylinder::new(1.0, 2.0);
        // Infinite cylinder centered at (0,0) with radius 1 → total 2, distance 0 → overlap
        assert!(c.intersects_infinite_cylinder(1.0, [0.0, 0.0]));
    }

    #[test]
    fn test_cylinder_intersects_infinite_cylinder_no_overlap() {
        let c = Cylinder::new(1.0, 2.0);
        // Infinite cylinder centered far away
        assert!(!c.intersects_infinite_cylinder(0.5, [5.0, 0.0]));
    }

    #[test]
    fn test_cylinder_random_surface_points_count() {
        let c = Cylinder::new(1.0, 2.0);
        let pts = c.random_surface_points(50, 12345);
        assert_eq!(pts.len(), 50);
    }

    #[test]
    fn test_cylinder_random_surface_points_on_surface() {
        let c = Cylinder::new(1.0, 2.0);
        let pts = c.random_surface_points(80, 99);
        for p in &pts {
            let sdf = c.sdf(*p);
            assert!(sdf.abs() < 0.01, "point {:?} has sdf={sdf}", p);
        }
    }

    #[test]
    fn test_cylinder_support_array_matches_support() {
        let c = Cylinder::new(2.0, 3.0);
        let d = [1.0, 1.0, 0.0];
        let a = c.support_array(d);
        let b = c.support(d);
        assert_eq!(a, b);
    }

    #[test]
    fn test_cylinder_full_height() {
        let c = Cylinder::new(1.0, 3.0);
        assert!((c.full_height() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_radius_of_gyration_y() {
        let c = Cylinder::new(2.0, 1.0);
        // k_y = r/sqrt(2) = 2/sqrt(2) = sqrt(2)
        assert!((c.radius_of_gyration_y() - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_radius_of_gyration_x() {
        let c = Cylinder::new(1.0, 1.0); // r=1, h=2
        // k_x = sqrt((3 + 4)/12) = sqrt(7/12)
        let expected = (7.0_f64 / 12.0).sqrt();
        assert!((c.radius_of_gyration_x(1.0) - expected).abs() < 1e-12);
    }

    // ── Extended tests for new methods ──────────────────────────────────────

    #[test]
    fn test_cylinder_volume_swept_zero_delta() {
        let c = Cylinder::new(1.0, 2.0);
        // Delta = 0 → swept volume = cylinder volume
        let swept = c.volume_swept([0.0, 0.0, 0.0]);
        assert!((swept - c.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_volume_swept_positive_delta() {
        let c = Cylinder::new(1.0, 1.0);
        let swept = c.volume_swept([1.0, 0.0, 0.0]);
        // Additional contribution: π*r²*1 = π
        assert!(swept > c.volume());
        assert!((swept - c.volume() - PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_cap_uv_center() {
        let c = Cylinder::new(2.0, 1.0);
        let uv = c.cap_uv([0.0, 2.0, 0.0]);
        assert!((uv[0] - 0.5).abs() < 1e-10);
        assert!((uv[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_cap_uv_edge() {
        let c = Cylinder::new(2.0, 1.0);
        // Point at (r, h, 0) → u = 0.5 + 2/4 = 0.5 + 0.5 = 1.0 → clamped to 1
        let uv = c.cap_uv([2.0, 1.0, 0.0]);
        assert!((uv[0] - 1.0).abs() < 1e-10);
        assert!((uv[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_lateral_uv_bottom_left() {
        let c = Cylinder::new(1.0, 1.0);
        // At θ=0 (x=1, z=0), y=-half_height: u≈0, v=0
        let uv = c.lateral_uv([1.0, -1.0, 0.0]);
        assert!(uv[1].abs() < 1e-10, "v={}", uv[1]);
    }

    #[test]
    fn test_cylinder_lateral_uv_top() {
        let c = Cylinder::new(1.0, 1.0);
        let uv = c.lateral_uv([1.0, 1.0, 0.0]);
        assert!((uv[1] - 1.0).abs() < 1e-10, "v={}", uv[1]);
    }

    #[test]
    fn test_cylinder_hollow_inside_ring() {
        let c = Cylinder::new(2.0, 1.0);
        // Point at r=1.5 (between inner=1 and outer=2)
        assert!(c.contains_point_hollow([1.5, 0.0, 0.0], 1.0));
    }

    #[test]
    fn test_cylinder_hollow_outside_outer() {
        let c = Cylinder::new(2.0, 1.0);
        assert!(!c.contains_point_hollow([3.0, 0.0, 0.0], 1.0));
    }

    #[test]
    fn test_cylinder_hollow_inside_inner_hole() {
        let c = Cylinder::new(2.0, 1.0);
        // Point at r=0.5, inside the hole
        assert!(!c.contains_point_hollow([0.5, 0.0, 0.0], 1.0));
    }

    #[test]
    fn test_cylinder_volume_hollow_full_solid() {
        // Inner radius = 0 → volume should equal solid volume
        let c = Cylinder::new(1.0, 1.0);
        assert!((c.volume_hollow(0.0) - c.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_volume_hollow_positive() {
        let c = Cylinder::new(2.0, 1.0);
        let v = c.volume_hollow(1.0);
        assert!(v > 0.0);
        // π(4-1)*2 = 6π
        assert!((v - 6.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_frustum_slant_height_equal_radii() {
        let c = Cylinder::new(1.0, 1.0);
        // Equal radii → slant = full height = 2
        let sl = c.frustum_slant_height(1.0, 1.0);
        assert!((sl - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_frustum_slant_height_cone() {
        let c = Cylinder::new(1.0, 1.0);
        // r_top=0, r_bottom=1: slant = sqrt(4+1) = sqrt(5)
        let sl = c.frustum_slant_height(0.0, 1.0);
        assert!((sl - 5.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_frustum_lateral_area() {
        let c = Cylinder::new(1.0, 1.0); // h=2
        // Equal radii r=1: lateral area = 2πr*h = 2π*1*2 = 4π
        let area = c.frustum_lateral_area(1.0, 1.0);
        assert!((area - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_frustum_volume_cylinder() {
        let c = Cylinder::new(1.0, 1.0);
        // Both radii equal to cylinder radius → frustum volume = cylinder volume
        let fv = c.frustum_volume(1.0, 1.0);
        assert!((fv - c.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_frustum_volume_cone() {
        let c = Cylinder::new(1.0, 1.0);
        // r_top=0, r_bottom=1: frustum volume = cone volume = πr²h/3 = π*1*2/3
        let fv = c.frustum_volume(0.0, 1.0);
        let expected = PI * 1.0 * 2.0 / 3.0;
        assert!((fv - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_oblique_cap_area_zero_tilt() {
        let c = Cylinder::new(1.0, 1.0);
        // Zero tilt: cap area = π r²
        assert!((c.oblique_cap_area(0.0) - PI).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_oblique_cap_area_larger_at_tilt() {
        let c = Cylinder::new(1.0, 1.0);
        assert!(c.oblique_cap_area(0.5) > c.oblique_cap_area(0.0));
    }

    #[test]
    fn test_cylinder_ruled_surface_point_on_surface() {
        let c = Cylinder::new(1.0, 2.0);
        // Check a few points
        for i in 0..8 {
            let u = 2.0 * PI * i as f64 / 8.0;
            for j in 0..=4 {
                let t = j as f64 / 4.0;
                let p = c.ruled_surface_point(u, t);
                let sdf = c.sdf(p);
                assert!(sdf.abs() < 1e-9, "ruled sdf={sdf}");
            }
        }
    }

    #[test]
    fn test_cylinder_stacked_volume() {
        let c1 = Cylinder::new(1.0, 1.0);
        let c2 = Cylinder::new(1.0, 2.0);
        let combined = c1.stacked_volume(&c2);
        assert!((combined - c1.volume() - c2.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_stacked_half_height() {
        let c1 = Cylinder::new(1.0, 1.0);
        let c2 = Cylinder::new(1.0, 2.0);
        assert!((c1.stacked_half_height(&c2) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_stacked_bounding_radius() {
        let c1 = Cylinder::new(1.0, 1.0);
        let c2 = Cylinder::new(3.0, 0.5);
        assert!((c1.stacked_bounding_radius(&c2) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_lateral_normal_at_x_axis() {
        let c = Cylinder::new(1.0, 2.0);
        let n = c.lateral_normal_at([1.0, 0.5, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-10);
        assert!(n[1].abs() < 1e-12);
        assert!(n[2].abs() < 1e-12);
    }

    #[test]
    fn test_cylinder_lateral_normal_unit() {
        let c = Cylinder::new(1.0, 2.0);
        let n = c.lateral_normal_at([0.6, 1.0, 0.8]);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_cylinder_sdf_gradient_outside_points_outward() {
        let c = Cylinder::new(1.0, 2.0);
        let g = c.sdf_gradient([3.0, 0.0, 0.0]);
        // Gradient should point outward: positive x component
        assert!(g[0] > 0.5, "gx={}", g[0]);
    }

    #[test]
    fn test_cylinder_lateral_ray_intersection_count_two() {
        let c = Cylinder::new(1.0, 2.0);
        let count = c.lateral_ray_intersection_count([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(count, 2);
    }

    #[test]
    fn test_cylinder_lateral_ray_intersection_count_miss() {
        let c = Cylinder::new(1.0, 2.0);
        let count = c.lateral_ray_intersection_count([0.0, 5.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(count, 0);
    }

    #[test]
    fn test_cylinder_top_rim_points_count() {
        let c = Cylinder::new(1.0, 2.0);
        let pts = c.top_rim_points(16);
        assert_eq!(pts.len(), 16);
    }

    #[test]
    fn test_cylinder_top_rim_points_at_top() {
        let c = Cylinder::new(1.0, 2.0);
        let pts = c.top_rim_points(8);
        for p in &pts {
            assert!((p[1] - c.half_height).abs() < 1e-12);
            let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
            assert!((xz - c.radius).abs() < 1e-9);
        }
    }
}
