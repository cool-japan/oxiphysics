// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cone shape (Y-axis aligned, apex at top).

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::f64::consts::PI;

/// A cone defined by radius and half-height along the Y axis.
/// Base is at y = -half_height, apex at y = +half_height.
#[derive(Debug, Clone)]
pub struct Cone {
    /// Radius of the base.
    pub radius: Real,
    /// Half-height along the Y axis.
    pub half_height: Real,
}

impl Cone {
    /// Create a new cone.
    pub fn new(radius: Real, half_height: Real) -> Self {
        Self {
            radius,
            half_height,
        }
    }

    /// Volume: πr²h/3 (h = full height = 2*half_height).
    pub fn volume_explicit(&self) -> Real {
        let h = 2.0 * self.half_height;
        PI * self.radius * self.radius * h / 3.0
    }

    /// Surface area: base disk + lateral surface = πr² + πr*slant.
    /// slant = sqrt(r² + h²) where h = full height.
    pub fn surface_area(&self) -> Real {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        let slant = (r * r + h * h).sqrt();
        PI * r * r + PI * r * slant
    }

    /// Inertia tensor as \[\[f64;3\\];3] row-major.
    pub fn inertia_tensor_array(&self, mass: f64) -> [[f64; 3]; 3] {
        let r2 = self.radius * self.radius;
        let h2 = (2.0 * self.half_height).powi(2);
        let iy = 3.0 * mass * r2 / 10.0;
        let ixz = mass * (3.0 * r2 + 2.0 * h2) / 20.0;
        [[ixz, 0.0, 0.0], [0.0, iy, 0.0], [0.0, 0.0, ixz]]
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

    /// Closest point on (or inside) the cone to point `p`.
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        let px = p[0];
        let py = p[1];
        let pz = p[2];

        // Clamp Y to cone range
        let cy = py.clamp(-self.half_height, self.half_height);

        // Radius at height cy: from apex (half_height) to base (-half_height)
        let h = 2.0 * self.half_height;
        let r_at_cy = self.radius * (self.half_height - cy) / h;

        // XZ distance
        let xz_len = (px * px + pz * pz).sqrt();

        if xz_len <= r_at_cy {
            // Project onto closest surface
            // Distance to base cap
            let dist_to_base = py - (-self.half_height);
            // Distance to lateral surface (approximate)
            let dist_to_side = r_at_cy - xz_len;
            if dist_to_base <= dist_to_side {
                [px, -self.half_height, pz]
            } else {
                // Clamp to cone surface at given Y
                if xz_len < 1e-12 {
                    [r_at_cy, cy, 0.0]
                } else {
                    let s = r_at_cy / xz_len;
                    [px * s, cy, pz * s]
                }
            }
        } else {
            // Outside XZ — project to cone side at clamped Y
            if xz_len < 1e-12 {
                [r_at_cy, cy, 0.0]
            } else {
                let s = r_at_cy / xz_len;
                [px * s, cy, pz * s]
            }
        }
    }

    // ── New methods ──

    /// GJK support function using plain arrays.
    /// Returns the farthest point on the cone in the given direction.
    pub fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let h = 2.0 * self.half_height;
        let xz_len = (direction[0] * direction[0] + direction[2] * direction[2]).sqrt();

        let apex_dot = direction[1] * self.half_height;
        let base_center_dot = -direction[1] * self.half_height;
        let base_rim_dot = base_center_dot + self.radius * xz_len;

        if apex_dot >= base_rim_dot {
            [0.0, self.half_height, 0.0]
        } else if xz_len > 1e-10 {
            let _sin_angle = self.radius / (self.radius * self.radius + h * h).sqrt();
            [
                self.radius * direction[0] / xz_len,
                -self.half_height,
                self.radius * direction[2] / xz_len,
            ]
        } else {
            [self.radius, -self.half_height, 0.0]
        }
    }

    /// Slant height of the cone: sqrt(r² + h²).
    pub fn slant_height(&self) -> f64 {
        let h = 2.0 * self.half_height;
        (self.radius * self.radius + h * h).sqrt()
    }

    /// Half angle of the cone in radians: atan(r / h).
    pub fn half_angle(&self) -> f64 {
        let h = 2.0 * self.half_height;
        (self.radius / h).atan()
    }

    /// Lateral (side) surface area only: πr * slant.
    pub fn lateral_surface_area(&self) -> f64 {
        PI * self.radius * self.slant_height()
    }

    /// Base area: πr².
    pub fn base_area(&self) -> f64 {
        PI * self.radius * self.radius
    }

    /// Returns true if `p` is inside the cone.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        if p[1] < -self.half_height || p[1] > self.half_height {
            return false;
        }
        let h = 2.0 * self.half_height;
        let r_at_y = self.radius * (self.half_height - p[1]) / h;
        let xz2 = p[0] * p[0] + p[2] * p[2];
        xz2 <= r_at_y * r_at_y
    }

    /// Signed distance from a point to the cone surface.
    /// Negative inside, positive outside.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        let cp = self.closest_point(p);
        let dx = p[0] - cp[0];
        let dy = p[1] - cp[1];
        let dz = p[2] - cp[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if self.contains_point(p) { -dist } else { dist }
    }

    /// Radius at a given y coordinate (clamped to cone range).
    pub fn radius_at_y(&self, y: f64) -> f64 {
        let clamped = y.clamp(-self.half_height, self.half_height);
        let h = 2.0 * self.half_height;
        self.radius * (self.half_height - clamped) / h
    }

    /// Create a truncated cone (frustum) by cutting at a given height ratio.
    /// `cut_ratio` is in \[0, 1\] where 0 = base, 1 = apex.
    /// Returns `(top_radius, bottom_radius, half_height)`.
    pub fn truncated_cone(&self, cut_ratio: f64) -> (f64, f64, f64) {
        let ratio = cut_ratio.clamp(0.0, 1.0);
        let top_r = self.radius * (1.0 - ratio);
        let bottom_r = self.radius;
        let hh = self.half_height * ratio;
        (top_r, bottom_r, hh)
    }

    /// Volume of a truncated cone (frustum).
    /// `cut_ratio` is in \[0, 1\].
    pub fn truncated_cone_volume(&self, cut_ratio: f64) -> f64 {
        let (r1, r2, _) = self.truncated_cone(cut_ratio);
        let h = 2.0 * self.half_height * cut_ratio.clamp(0.0, 1.0);
        PI * h * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0
    }

    /// Ray cast against only the base cap at y = -half_height.
    pub fn ray_cast_base_cap(
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

    /// Closest point on the cone's lateral surface to a point.
    /// This projects onto the side surface only (ignoring base cap).
    pub fn closest_point_on_lateral(&self, p: [f64; 3]) -> [f64; 3] {
        let xz_len = (p[0] * p[0] + p[2] * p[2]).sqrt();
        let h = 2.0 * self.half_height;

        let side_len = self.slant_height();
        let side_dx = self.radius / side_len;
        let side_dy = -h / side_len;

        let vx = xz_len;
        let vy = p[1] - self.half_height;

        let t = (vx * side_dx + vy * side_dy).clamp(0.0, side_len);
        let proj_xz = t * side_dx;
        let proj_y = self.half_height + t * side_dy;

        if xz_len < 1e-12 {
            [proj_xz, proj_y, 0.0]
        } else {
            let s = proj_xz / xz_len;
            [p[0] * s, proj_y, p[2] * s]
        }
    }

    // ── New expanded methods ──

    /// Analytic ray-cone intersection (same as `ray_cast` but returns plain arrays).
    /// Identical to `ray_cast_array` but additionally returns a boolean indicating
    /// whether the hit is on the lateral surface (`true`) or base cap (`false`).
    pub fn ray_cone_analytic(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3], bool)> {
        let o = Vec3::new(origin[0], origin[1], origin[2]);
        let d = Vec3::new(direction[0], direction[1], direction[2]);

        let h = 2.0 * self.half_height;
        let k = self.radius / h;
        let k2 = k * k;
        let apex_y = self.half_height;

        let fy = apex_y - o.y;
        let a = d.x * d.x + d.z * d.z - k2 * d.y * d.y;
        let b = 2.0 * (o.x * d.x + o.z * d.z + k2 * fy * d.y);
        let c_val = o.x * o.x + o.z * o.z - k2 * fy * fy;

        let mut best_t = f64::INFINITY;
        let mut best_is_lateral = true;

        if a.abs() > 1e-12 {
            let disc = b * b - 4.0 * a * c_val;
            if disc >= 0.0 {
                let sqrt_disc = disc.sqrt();
                for t in [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)] {
                    if t >= 0.0 && t <= max_toi {
                        let py = o.y + t * d.y;
                        if py >= -self.half_height && py <= self.half_height && t < best_t {
                            best_t = t;
                            best_is_lateral = true;
                        }
                    }
                }
            }
        }

        // Base cap
        if d.y.abs() > 1e-12 {
            let t = (-self.half_height - o.y) / d.y;
            if t >= 0.0 && t <= max_toi {
                let px = o.x + t * d.x;
                let pz = o.z + t * d.z;
                if px * px + pz * pz <= self.radius * self.radius && t < best_t {
                    best_t = t;
                    best_is_lateral = false;
                }
            }
        }

        if best_t.is_infinite() {
            return None;
        }

        let hit_pt = [
            origin[0] + best_t * direction[0],
            origin[1] + best_t * direction[1],
            origin[2] + best_t * direction[2],
        ];
        let normal = if best_is_lateral {
            let xz_len = (hit_pt[0] * hit_pt[0] + hit_pt[2] * hit_pt[2])
                .sqrt()
                .max(1e-30);
            let slope = self.radius / h;
            let nx = hit_pt[0] / xz_len;
            let nz = hit_pt[2] / xz_len;
            let len = (nx * nx + slope * slope + nz * nz).sqrt();
            [nx / len, slope / len, nz / len]
        } else {
            [0.0, -1.0, 0.0]
        };

        Some((best_t, normal, best_is_lateral))
    }

    /// Cone-plane intersection.
    ///
    /// Tests whether the cone intersects a plane defined by `plane_normal` (unit)
    /// and `plane_d` (signed distance from origin: dot(n, x) = plane_d).
    ///
    /// Returns `true` if any part of the cone (including base disk) is on or
    /// crosses the plane.
    pub fn intersects_plane(&self, plane_normal: [f64; 3], plane_d: f64) -> bool {
        // Evaluate the signed distance of the apex and the base rim extremes.
        let apex = [0.0_f64, self.half_height, 0.0];
        let apex_sd =
            plane_normal[0] * apex[0] + plane_normal[1] * apex[1] + plane_normal[2] * apex[2]
                - plane_d;

        // Base center
        let base_center_sd = -plane_normal[1] * self.half_height - plane_d;

        // Farthest base rim point in plane_normal direction
        let xz_proj =
            (plane_normal[0] * plane_normal[0] + plane_normal[2] * plane_normal[2]).sqrt();
        let rim_sd = base_center_sd + self.radius * xz_proj;
        let rim_sd_neg = base_center_sd - self.radius * xz_proj;

        // Intersection if any signed distance changes sign
        let signs = [apex_sd, rim_sd, rim_sd_neg];
        let any_pos = signs.iter().any(|&s| s >= 0.0);
        let any_neg = signs.iter().any(|&s| s <= 0.0);
        any_pos && any_neg
    }

    /// Compute a frustum (truncated cone) from this cone.
    ///
    /// Cuts the cone at heights `y_bottom` and `y_top` (in cone-local space,
    /// where -half_height ≤ y_bottom < y_top ≤ half_height).
    ///
    /// Returns `(bottom_radius, top_radius, frustum_half_height)`.
    pub fn frustum_from_cone(&self, y_bottom: f64, y_top: f64) -> (f64, f64, f64) {
        let yb = y_bottom.clamp(-self.half_height, self.half_height);
        let yt = y_top.clamp(yb, self.half_height);
        let r_bottom = self.radius_at_y(yb);
        let r_top = self.radius_at_y(yt);
        let frustum_hh = (yt - yb) * 0.5;
        (r_bottom, r_top, frustum_hh)
    }

    /// Cone SDF (Signed Distance Function).
    ///
    /// Returns negative inside the cone, positive outside.
    /// Same as `signed_distance` but exposed under the SDF naming convention.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        self.signed_distance(p)
    }

    /// Cone-sphere intersection test.
    ///
    /// Returns `true` if the sphere (center + radius) overlaps the cone.
    /// Uses the closest-point approach: if the closest point on the cone to
    /// the sphere center is within `sphere_radius`, they intersect.
    pub fn intersects_sphere(&self, sphere_center: [f64; 3], sphere_radius: f64) -> bool {
        let cp = self.closest_point(sphere_center);
        let dx = sphere_center[0] - cp[0];
        let dy = sphere_center[1] - cp[1];
        let dz = sphere_center[2] - cp[2];
        let dist_sq = dx * dx + dy * dy + dz * dz;
        dist_sq <= sphere_radius * sphere_radius
    }

    /// Lateral surface area of a frustum derived from this cone.
    ///
    /// `y_bottom` and `y_top` are in cone-local Y coordinates.
    /// Lateral area = π * (r_bottom + r_top) * slant_height_frustum.
    pub fn frustum_lateral_area(&self, y_bottom: f64, y_top: f64) -> f64 {
        let (r_bottom, r_top, _hh) = self.frustum_from_cone(y_bottom, y_top);
        let dy = (y_top - y_bottom).abs();
        let dr = (r_bottom - r_top).abs();
        let slant = (dy * dy + dr * dr).sqrt();
        PI * (r_bottom + r_top) * slant
    }

    // ── Extended geometry: double cone, oblique, unfolding, bounding cone ──

    /// Apex angle of the cone (full angle at the apex), in radians.
    ///
    /// The full apex angle is `2 * atan(r / h)` where h is the full height.
    pub fn apex_angle(&self) -> f64 {
        2.0 * self.half_angle()
    }

    /// Create a double cone (bicone) by reflecting this cone across its base.
    ///
    /// Returns the parameters of the upper half (this cone) and the lower half
    /// (same dimensions, mirrored).  Both have apex at `±half_height`.
    pub fn double_cone_params(&self) -> (f64, f64) {
        (self.radius, self.half_height)
    }

    /// Returns `true` if `p` is inside the double cone (both upper and lower).
    ///
    /// The double cone is symmetric about y=0.  Upper half: apex at
    /// `+half_height`, lower half: apex at `-half_height`.
    pub fn contains_point_double_cone(&self, p: [f64; 3]) -> bool {
        let h = 2.0 * self.half_height;
        let xz2 = p[0] * p[0] + p[2] * p[2];
        let y_abs = p[1].abs();
        if y_abs > self.half_height {
            return false;
        }
        // At y_abs, radius = r * y_abs / half_height
        let r_at = self.radius * y_abs / self.half_height;
        xz2 <= r_at * r_at + 1e-14 * h
    }

    /// Lateral surface unfolding: convert `(u, t)` to a point in the unrolled
    /// lateral surface plane.
    ///
    /// In the unrolled sector the slant height becomes the radial distance and
    /// the azimuth `u` is scaled by `sin(half_angle)`.  Returns `(x_unrolled, y_unrolled)`.
    pub fn unfold_lateral(&self, u: f64, t_slant: f64) -> [f64; 2] {
        // When unrolled, the cone becomes a sector with radius = slant_height
        // and arc angle = 2π * sin(half_angle).
        let sin_alpha = self.radius / self.slant_height().max(1e-14);
        let phi = u * sin_alpha;
        [t_slant * phi.cos(), t_slant * phi.sin()]
    }

    /// Bounding cone for a set of points: returns the smallest half-angle cone
    /// with its apex at the origin and axis along +Y that contains all points.
    ///
    /// Returns `(axis_half_angle_radians, axis_elevation)` where the axis is
    /// fixed to +Y.
    pub fn bounding_cone_of_points(points: &[[f64; 3]]) -> (f64, f64) {
        if points.is_empty() {
            return (0.0, 0.0);
        }
        let mut max_half_angle: f64 = 0.0;
        let mut max_elev: f64 = 0.0;
        for &p in points {
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if len < 1e-14 {
                continue;
            }
            // Half-angle from +Y axis
            let cos_a = (p[1] / len).clamp(-1.0, 1.0);
            let half_angle = cos_a.acos();
            let elev = p[1] / len;
            if half_angle > max_half_angle {
                max_half_angle = half_angle;
                max_elev = elev;
            }
        }
        (max_half_angle, max_elev)
    }

    /// Volume swept when the cone translates by vector `delta`.
    ///
    /// Approximate: V_cone + base_area * |delta|.
    pub fn volume_swept(&self, delta: [f64; 3]) -> f64 {
        let dist = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
        self.volume() + self.base_area() * dist
    }

    /// Cone SDF gradient (approximate unit normal direction).
    ///
    /// Returns the gradient of the SDF at `p`, suitable as the collision normal.
    pub fn sdf_gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let eps = 1e-5;
        let s0 = self.sdf(p);
        let gx = (self.sdf([p[0] + eps, p[1], p[2]]) - s0) / eps;
        let gy = (self.sdf([p[0], p[1] + eps, p[2]]) - s0) / eps;
        let gz = (self.sdf([p[0], p[1], p[2] + eps]) - s0) / eps;
        let len = (gx * gx + gy * gy + gz * gz).sqrt().max(1e-14);
        [gx / len, gy / len, gz / len]
    }

    /// Ray-cone intersection count (lateral + base cap).
    ///
    /// Returns the number of valid positive-t intersections within `[0, max_toi]`.
    pub fn ray_intersection_count(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> usize {
        let mut count = 0usize;
        let h = 2.0 * self.half_height;
        let k = self.radius / h;
        let k2 = k * k;
        let apex_y = self.half_height;
        let fy = apex_y - origin[1];
        let a = direction[0] * direction[0] + direction[2] * direction[2]
            - k2 * direction[1] * direction[1];
        let b =
            2.0 * (origin[0] * direction[0] + origin[2] * direction[2] + k2 * fy * direction[1]);
        let c = origin[0] * origin[0] + origin[2] * origin[2] - k2 * fy * fy;
        if a.abs() > 1e-12 {
            let disc = b * b - 4.0 * a * c;
            if disc >= 0.0 {
                let sq = disc.sqrt();
                for t in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                    if t >= 0.0 && t <= max_toi {
                        let py = origin[1] + t * direction[1];
                        if py >= -self.half_height && py <= self.half_height {
                            count += 1;
                        }
                    }
                }
            }
        }
        // Base cap
        if direction[1].abs() > 1e-12 {
            let t = (-self.half_height - origin[1]) / direction[1];
            if t >= 0.0 && t <= max_toi {
                let px = origin[0] + t * direction[0];
                let pz = origin[2] + t * direction[2];
                if px * px + pz * pz <= self.radius * self.radius {
                    count += 1;
                }
            }
        }
        count
    }

    /// Compute the solid angle subtended by the cone at the apex.
    ///
    /// For a cone with half-angle `α`, the solid angle is `2π(1 - cos α)`.
    pub fn solid_angle(&self) -> f64 {
        let alpha = self.half_angle();
        2.0 * PI * (1.0 - alpha.cos())
    }

    /// Generate points on the base circle of the cone.
    ///
    /// Returns `n` equally-spaced points on the base rim at `y = -half_height`.
    pub fn base_rim_points(&self, n: usize) -> Vec<[f64; 3]> {
        let n = n.max(3);
        (0..n)
            .map(|i| {
                let t = 2.0 * PI * i as f64 / n as f64;
                [
                    self.radius * t.cos(),
                    -self.half_height,
                    self.radius * t.sin(),
                ]
            })
            .collect()
    }

    /// Test whether a sphere at `center` with `sphere_radius` is entirely
    /// inside the cone (a stricter test than `intersects_sphere`).
    pub fn sphere_inside_cone(&self, center: [f64; 3], sphere_radius: f64) -> bool {
        // The sphere is inside if all points on the sphere are inside the cone.
        // Approximate: check the sphere center is inside and the distance to
        // the cone surface is at least sphere_radius.
        let sd = self.sdf(center);
        sd + sphere_radius <= 0.0
    }

    /// Sample a random point on the cone's lateral surface using a deterministic PRNG.
    pub fn random_lateral_surface_point(&self, seed: u64) -> [f64; 3] {
        // Use a simple xorshift64 to get reproducible results
        let mut state = seed;
        let xorshift = |s: &mut u64| -> f64 {
            *s ^= *s << 13;
            *s ^= *s >> 7;
            *s ^= *s << 17;
            (*s as f64) / (u64::MAX as f64)
        };

        // Sample t along slant height proportional to circumference (inverse CDF)
        // Circumference at height t: 2π * r(t), where r(t) = r*(1 - t/h)
        // CDF of area up to t: t - t²/(2h), so sample with u = t/h - (t/h)²/2
        // Invert: t/h = 1 - sqrt(1-u) (for top-aligned cone where apex has r=0)
        let u = xorshift(&mut state);
        let frac_t = 1.0 - (1.0 - u).sqrt(); // normalized height fraction from apex
        let y = self.half_height - frac_t * 2.0 * self.half_height;
        let r_at_y = self.radius_at_y(y);
        let theta = xorshift(&mut state) * 2.0 * PI;
        [r_at_y * theta.cos(), y, r_at_y * theta.sin()]
    }
}

impl Shape for Cone {
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
        let h = 2.0 * self.half_height;
        let xz_len = (direction.x * direction.x + direction.z * direction.z).sqrt();

        // Check if apex is more extreme
        let apex_dot = direction.y * self.half_height;
        let base_center_dot = -direction.y * self.half_height;
        let base_rim_dot = base_center_dot + self.radius * xz_len;

        if apex_dot >= base_rim_dot {
            Vec3::new(0.0, self.half_height, 0.0)
        } else if xz_len > 1e-10 {
            let sin_angle = self.radius / (self.radius * self.radius + h * h).sqrt();
            let _ = sin_angle; // used for more precise support, but rim suffices
            Vec3::new(
                self.radius * direction.x / xz_len,
                -self.half_height,
                self.radius * direction.z / xz_len,
            )
        } else {
            Vec3::new(self.radius, -self.half_height, 0.0)
        }
    }

    fn volume(&self) -> Real {
        let h = 2.0 * self.half_height;
        PI * self.radius * self.radius * h / 3.0
    }

    fn center_of_mass(&self) -> Vec3 {
        // Center of mass is at 1/4 height from the base
        Vec3::new(0.0, -self.half_height * 0.5, 0.0)
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let r2 = self.radius * self.radius;
        // Full height h = 2 * half_height
        let h2 = (2.0 * self.half_height).powi(2);
        // I_zz (spin axis) = 3*m*r²/10
        let iy = 3.0 * mass * r2 / 10.0;
        // I_xx = I_yy = m*(3r² + 2h²)/20  (h = full height)
        let ixz = mass * (3.0 * r2 + 2.0 * h2) / 20.0;
        Mat3::new(ixz, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ixz)
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        // Cone: apex at y=+half_height, base at y=-half_height, radius r.
        // For a point P on the cone surface: r_at_y = radius * (half_height - y) / (2*half_height)
        // Constraint: x² + z² = [r * (h - y) / (2h)]²
        // Let k = r / (2*half_height), apex_y = half_height
        // (ox + t*dx)² + (oz + t*dz)² = [k*(apex_y - (oy + t*dy))]²
        let h = 2.0 * self.half_height;
        let k = self.radius / h;
        let k2 = k * k;
        let apex_y = self.half_height;

        let ox = ray_origin.x;
        let oy = ray_origin.y;
        let oz = ray_origin.z;
        let dx = ray_direction.x;
        let dy = ray_direction.y;
        let dz = ray_direction.z;

        // (dx²+dz² - k²*dy²) t² + 2*(ox*dx+oz*dz + k²*(apex_y-oy)*dy) t +
        //   (ox²+oz² - k²*(apex_y-oy)²) = 0
        let fy = apex_y - oy;
        let a = dx * dx + dz * dz - k2 * dy * dy;
        let b = 2.0 * (ox * dx + oz * dz + k2 * fy * dy);
        let c_val = ox * ox + oz * oz - k2 * fy * fy;

        let mut best: Option<RayHit> = None;

        let try_t = |t: Real, best: &mut Option<RayHit>| {
            if t < 0.0 || t > max_toi {
                return;
            }
            let p = ray_origin + ray_direction * t;
            if p.y < -self.half_height || p.y > self.half_height {
                return;
            }
            // Normal: outward in XZ, with slope component
            let xz_len = (p.x * p.x + p.z * p.z).sqrt();
            let normal = if xz_len > 1e-12 {
                let slope = self.radius / h;
                Vec3::new(p.x / xz_len, slope, p.z / xz_len).normalize()
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            if best.as_ref().is_none_or(|prev| t < prev.toi) {
                *best = Some(RayHit {
                    point: p,
                    normal,
                    toi: t,
                });
            }
        };

        if a.abs() > 1e-12 {
            let disc = b * b - 4.0 * a * c_val;
            if disc >= 0.0 {
                let sqrt_disc = disc.sqrt();
                try_t((-b - sqrt_disc) / (2.0 * a), &mut best);
                try_t((-b + sqrt_disc) / (2.0 * a), &mut best);
            }
        } else if b.abs() > 1e-12 {
            try_t(-c_val / b, &mut best);
        }

        // Base cap at y = -half_height
        if ray_direction.y.abs() > 1e-12 {
            let t = (-self.half_height - ray_origin.y) / ray_direction.y;
            if t >= 0.0 && t <= max_toi {
                let p = ray_origin + ray_direction * t;
                if p.x * p.x + p.z * p.z <= self.radius * self.radius {
                    let normal = Vec3::new(0.0, -1.0, 0.0);
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

        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cone_volume() {
        let c = Cone::new(1.0, 1.5);
        let expected = PI * 1.0 * 3.0 / 3.0;
        assert!((c.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_surface_area() {
        // r=1, h=2 (half_height=1): base=π, slant=sqrt(1+4)=sqrt(5), lateral=π*sqrt(5)
        let c = Cone::new(1.0, 1.0);
        let slant = (1.0_f64 + 4.0_f64).sqrt();
        let expected = PI + PI * slant;
        assert!((c.surface_area() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_support_apex() {
        let c = Cone::new(1.0, 2.0);
        let sp = c.support_point(&Vec3::new(0.0, 1.0, 0.0));
        assert!((sp.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_inertia_symmetry() {
        let c = Cone::new(1.5, 2.0);
        let it = c.inertia_tensor(3.0);
        // I_xx == I_zz (both are ixz), I_yy is spin axis
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
    }

    #[test]
    fn test_cone_inertia_array() {
        let c = Cone::new(1.0, 1.0);
        let it = c.inertia_tensor_array(1.0);
        assert!(it[0][0] > 0.0);
        assert!(it[1][1] > 0.0);
        assert!((it[0][0] - it[2][2]).abs() < 1e-10);
    }

    #[test]
    fn test_cone_raycast_side() {
        // Ray along +X hitting the side of the cone
        let c = Cone::new(1.0, 2.0); // apex at y=2, base at y=-2, radius 1
        // Ray from far left at y=0 (midpoint height), should intersect side
        let origin = Vec3::new(-5.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        // At y=0, the cone radius is r*(h - y)/(2h) = 1*(2-0)/4 = 0.5
        let hit = c.ray_cast(&origin, &dir, 100.0);
        assert!(hit.is_some(), "expected a hit on cone side");
        let hit = hit.unwrap();
        assert!(hit.toi > 0.0);
        assert!(
            (hit.point.x + 0.5).abs() < 1e-6,
            "expected hit at x=-0.5, got x={}",
            hit.point.x
        );
    }

    #[test]
    fn test_cone_raycast_base_cap() {
        // Ray from below hitting the base cap
        let c = Cone::new(1.0, 1.0); // apex at y=1, base at y=-1
        let origin = Vec3::new(0.0, -5.0, 0.0);
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let hit = c.ray_cast(&origin, &dir, 100.0);
        assert!(hit.is_some(), "expected a hit on cone base");
        let hit = hit.unwrap();
        assert!(
            (hit.toi - 4.0).abs() < 1e-6,
            "expected toi=4, got {}",
            hit.toi
        );
        assert!((hit.normal.y + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_closest_point_outside() {
        let c = Cone::new(1.0, 1.0); // apex y=1, base y=-1, r=1
        // Point far in X at y=0: cone radius at y=0 is 0.5
        let cp = c.closest_point([5.0, 0.0, 0.0]);
        // Should project to cone side at y=0, x=0.5
        assert!(
            (cp[0] - 0.5).abs() < 1e-6,
            "expected cp[0]=0.5, got {}",
            cp[0]
        );
        assert!((cp[1]).abs() < 1e-6);
    }

    // ── New tests ──

    #[test]
    fn test_cone_support_array_apex() {
        let c = Cone::new(1.0, 2.0);
        let sp = c.support([0.0, 1.0, 0.0]);
        assert!((sp[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_support_array_base_rim() {
        let c = Cone::new(1.0, 2.0);
        let sp = c.support([1.0, -1.0, 0.0]);
        assert!((sp[0] - 1.0).abs() < 1e-10);
        assert!((sp[1] + 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_slant_height() {
        let c = Cone::new(3.0, 2.0); // h=4
        let expected = (9.0 + 16.0_f64).sqrt(); // 5
        assert!((c.slant_height() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_half_angle() {
        let c = Cone::new(1.0, 1.0); // r=1, h=2
        let expected = (0.5_f64).atan(); // atan(r/h) = atan(0.5)
        assert!((c.half_angle() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_lateral_surface_area() {
        let c = Cone::new(1.0, 1.0);
        let slant = c.slant_height();
        let expected = PI * 1.0 * slant;
        assert!((c.lateral_surface_area() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_base_area() {
        let c = Cone::new(2.0, 1.0);
        assert!((c.base_area() - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_cone_surface_area_decomposition() {
        let c = Cone::new(1.5, 2.0);
        let total = c.surface_area();
        let decomposed = c.lateral_surface_area() + c.base_area();
        assert!((total - decomposed).abs() < 1e-10);
    }

    #[test]
    fn test_cone_contains_point() {
        let c = Cone::new(1.0, 1.0); // apex y=1, base y=-1, r=1
        assert!(c.contains_point([0.0, 0.0, 0.0])); // center, r_at_0 = 0.5
        assert!(c.contains_point([0.0, -1.0, 0.0])); // base center
        assert!(!c.contains_point([0.0, 2.0, 0.0])); // above apex
        assert!(!c.contains_point([1.0, 0.0, 0.0])); // at y=0, r=0.5, x=1 > 0.5
    }

    #[test]
    fn test_cone_contains_point_base_rim() {
        let c = Cone::new(1.0, 1.0);
        // At y=-1, r_at_y = 1.0, so (1,0) should be on boundary
        assert!(c.contains_point([1.0, -1.0, 0.0]));
        assert!(!c.contains_point([1.1, -1.0, 0.0]));
    }

    #[test]
    fn test_cone_signed_distance_inside() {
        let c = Cone::new(2.0, 2.0);
        let d = c.signed_distance([0.0, -2.0, 0.0]); // base center
        assert!(d <= 0.0);
    }

    #[test]
    fn test_cone_signed_distance_outside() {
        let c = Cone::new(1.0, 1.0);
        let d = c.signed_distance([5.0, 0.0, 0.0]);
        assert!(d > 0.0);
    }

    #[test]
    fn test_cone_radius_at_y() {
        let c = Cone::new(2.0, 2.0); // h=4
        assert!((c.radius_at_y(-2.0) - 2.0).abs() < 1e-10); // base
        assert!((c.radius_at_y(2.0)).abs() < 1e-10); // apex
        assert!((c.radius_at_y(0.0) - 1.0).abs() < 1e-10); // midpoint
    }

    #[test]
    fn test_cone_truncated_cone() {
        let c = Cone::new(2.0, 2.0);
        let (top_r, bot_r, hh) = c.truncated_cone(0.5);
        assert!((bot_r - 2.0).abs() < 1e-10);
        assert!((top_r - 1.0).abs() < 1e-10); // 2 * (1 - 0.5)
        assert!((hh - 1.0).abs() < 1e-10); // 2 * 0.5
    }

    #[test]
    fn test_cone_truncated_cone_full() {
        let c = Cone::new(2.0, 2.0);
        let (top_r, _, _) = c.truncated_cone(1.0);
        assert!(top_r.abs() < 1e-10); // full cone, top radius = 0
    }

    #[test]
    fn test_cone_truncated_cone_volume() {
        let c = Cone::new(1.0, 1.0);
        // Full cone volume = π * 1 * 2 / 3
        let full_vol = c.volume();
        // Truncated at ratio=1.0 should give the full cone volume
        let trunc_vol = c.truncated_cone_volume(1.0);
        assert!((trunc_vol - full_vol).abs() < 1e-10);
    }

    #[test]
    fn test_cone_ray_cast_base_cap_array() {
        let c = Cone::new(1.0, 1.0);
        let result = c.ray_cast_base_cap([0.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, n) = result.unwrap();
        assert!((t - 4.0).abs() < 1e-10);
        assert!((n[1] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_ray_cast_base_cap_miss() {
        let c = Cone::new(1.0, 1.0);
        let result = c.ray_cast_base_cap([5.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_cone_closest_point_on_lateral() {
        let c = Cone::new(1.0, 1.0);
        // Point far out along X at y=0
        let cp = c.closest_point_on_lateral([5.0, 0.0, 0.0]);
        // At y=0 on the side, r should be about 0.5
        // But closest_point_on_lateral projects onto the actual side line
        assert!(cp[0] > 0.0);
        // The point should be on the cone side
        let r = (cp[0] * cp[0] + cp[2] * cp[2]).sqrt();
        let expected_r = c.radius_at_y(cp[1]);
        assert!((r - expected_r).abs() < 0.1);
    }

    #[test]
    fn test_cone_closest_point_on_lateral_apex() {
        let c = Cone::new(1.0, 1.0);
        // Point above the apex
        let cp = c.closest_point_on_lateral([0.0, 5.0, 0.0]);
        // Should project to apex
        assert!((cp[1] - 1.0).abs() < 1e-10);
        assert!(cp[0].abs() < 1e-10);
    }

    #[test]
    fn test_cone_volume_explicit_matches() {
        let c = Cone::new(1.5, 2.5);
        assert!((c.volume_explicit() - c.volume()).abs() < 1e-10);
    }

    // ── Expanded tests ──

    #[test]
    fn test_cone_ray_analytic_lateral() {
        let c = Cone::new(1.0, 2.0);
        let result = c.ray_cone_analytic([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some(), "should hit lateral surface");
        let (t, n, is_lateral) = result.unwrap();
        assert!(t > 0.0);
        assert!(is_lateral, "should be lateral hit");
        assert!(
            n[1] > 0.0,
            "normal should have positive Y component for lateral cone"
        );
    }

    #[test]
    fn test_cone_ray_analytic_base() {
        let c = Cone::new(1.0, 1.0);
        let result = c.ray_cone_analytic([0.0, -5.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, n, is_lateral) = result.unwrap();
        assert!((t - 4.0).abs() < 1e-6, "toi should be 4, got {t}");
        assert!(!is_lateral);
        assert!((n[1] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_plane_intersection_horizontal_through_base() {
        let c = Cone::new(1.0, 1.0);
        // Horizontal plane at y = -1 (exactly the base) – normal is (0,1,0), d=-1
        let hit = c.intersects_plane([0.0, 1.0, 0.0], -1.0);
        assert!(hit, "plane through base should intersect cone");
    }

    #[test]
    fn test_cone_plane_intersection_above_apex() {
        let c = Cone::new(1.0, 1.0);
        // Plane far above apex: normal (0,1,0), d=5 → all cone points have y·1 < 5
        let hit = c.intersects_plane([0.0, 1.0, 0.0], 5.0);
        assert!(!hit, "plane far above should not intersect cone");
    }

    #[test]
    fn test_cone_frustum_from_cone() {
        let c = Cone::new(2.0, 2.0); // apex y=2, base y=-2, r=2, h=4
        // Cut from y=-1 to y=1
        let (r_bot, r_top, hh) = c.frustum_from_cone(-1.0, 1.0);
        // r_at(-1) = 2*(2-(-1))/4 = 2*3/4 = 1.5
        assert!((r_bot - 1.5).abs() < 1e-10, "r_bot={r_bot}");
        // r_at(1) = 2*(2-1)/4 = 0.5
        assert!((r_top - 0.5).abs() < 1e-10, "r_top={r_top}");
        assert!((hh - 1.0).abs() < 1e-10, "hh={hh}");
    }

    #[test]
    fn test_cone_sdf_matches_signed_distance() {
        let c = Cone::new(1.0, 1.0);
        let p = [3.0, 0.0, 0.0];
        assert!((c.sdf(p) - c.signed_distance(p)).abs() < 1e-12);
    }

    #[test]
    fn test_cone_intersects_sphere_overlap() {
        let c = Cone::new(1.0, 1.0);
        // Large sphere centered at origin should overlap the cone
        assert!(c.intersects_sphere([0.0, 0.0, 0.0], 5.0));
    }

    #[test]
    fn test_cone_intersects_sphere_no_overlap() {
        let c = Cone::new(1.0, 1.0);
        // Small sphere far away
        assert!(!c.intersects_sphere([10.0, 0.0, 0.0], 0.1));
    }

    #[test]
    fn test_cone_frustum_lateral_area_full_cone() {
        let c = Cone::new(1.0, 1.0); // half_height=1, radius=1
        // Full lateral area as frustum from -1 to 1 (base to apex)
        let area = c.frustum_lateral_area(-1.0, 1.0);
        let expected = c.lateral_surface_area();
        // The frustum from base to apex should match full lateral area
        assert!(
            (area - expected).abs() < 1e-9,
            "frustum area {area} vs lateral {expected}"
        );
    }

    #[test]
    fn test_cone_random_lateral_surface_point_on_surface() {
        let c = Cone::new(1.0, 2.0);
        for seed in [1u64, 42, 99, 1337, 65535] {
            let p = c.random_lateral_surface_point(seed);
            // Point should be within cone Y range
            assert!(
                p[1] >= -c.half_height - 1e-9 && p[1] <= c.half_height + 1e-9,
                "y={} out of range",
                p[1]
            );
            // XZ radius should match cone radius at that Y
            let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
            let r_at = c.radius_at_y(p[1]);
            assert!((xz - r_at).abs() < 1e-9, "xz={xz}, r_at={r_at}");
        }
    }

    #[test]
    fn test_cone_ray_analytic_miss() {
        let c = Cone::new(1.0, 1.0);
        // Ray going away from cone
        let result = c.ray_cone_analytic([0.0, 10.0, 0.0], [0.0, 1.0, 0.0], 5.0);
        assert!(result.is_none(), "ray going away should miss");
    }

    // ── Extended tests for new methods ──────────────────────────────────────

    #[test]
    fn test_cone_apex_angle() {
        let c = Cone::new(1.0, 1.0); // r=1, h=2
        // apex angle = 2 * atan(1/2)
        let expected = 2.0 * (0.5_f64).atan();
        assert!((c.apex_angle() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_apex_angle_90_degrees() {
        // r=h → half_angle = 45° → apex_angle = 90°
        let c = Cone::new(1.0, 0.5); // h = 1.0
        // half_angle = atan(1/1) = π/4 → apex = π/2
        let expected = PI / 2.0;
        assert!((c.apex_angle() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_double_cone_params() {
        let c = Cone::new(1.5, 2.0);
        let (r, hh) = c.double_cone_params();
        assert!((r - 1.5).abs() < 1e-12);
        assert!((hh - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_cone_contains_point_double_cone_at_origin() {
        // Double cone: at y=0 the radius is 0, so origin is on the waist
        let c = Cone::new(1.0, 1.0);
        assert!(c.contains_point_double_cone([0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_cone_contains_point_double_cone_above_apex() {
        let c = Cone::new(1.0, 1.0);
        assert!(!c.contains_point_double_cone([0.0, 1.5, 0.0]));
    }

    #[test]
    fn test_cone_contains_point_double_cone_on_surface() {
        let c = Cone::new(1.0, 1.0);
        // At |y| = half_height, r = radius
        assert!(c.contains_point_double_cone([1.0, 1.0, 0.0]));
    }

    #[test]
    fn test_cone_unfold_lateral_origin() {
        let c = Cone::new(1.0, 1.0);
        // t_slant=0 (apex) should map to origin
        let [x, y] = c.unfold_lateral(0.0, 0.0);
        assert!(x.abs() < 1e-12);
        assert!(y.abs() < 1e-12);
    }

    #[test]
    fn test_cone_unfold_lateral_nonzero() {
        let c = Cone::new(1.0, 1.0);
        let [x, y] = c.unfold_lateral(0.0, 1.0);
        let len = (x * x + y * y).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "len={len}");
    }

    #[test]
    fn test_cone_bounding_cone_of_points_empty() {
        let (ha, _) = Cone::bounding_cone_of_points(&[]);
        assert_eq!(ha, 0.0);
    }

    #[test]
    fn test_cone_bounding_cone_of_points_single() {
        let pts = [[1.0, 1.0, 0.0]];
        let (ha, _) = Cone::bounding_cone_of_points(&pts);
        let expected = PI / 4.0; // 45° from +Y for (1,1,0) normalized
        assert!((ha - expected).abs() < 1e-9, "ha={ha} expected={expected}");
    }

    #[test]
    fn test_cone_bounding_cone_contains_all_points() {
        let pts = [[0.1, 1.0, 0.0], [0.0, 1.0, 0.1], [0.0, 1.0, -0.1]];
        let (max_ha, _) = Cone::bounding_cone_of_points(&pts);
        // All points within the bounding cone
        for &p in &pts {
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let cos_a = (p[1] / len).clamp(-1.0, 1.0);
            let ha = cos_a.acos();
            assert!(ha <= max_ha + 1e-12, "ha={ha} > max={max_ha}");
        }
    }

    #[test]
    fn test_cone_volume_swept_zero_delta() {
        let c = Cone::new(1.0, 1.0);
        let swept = c.volume_swept([0.0, 0.0, 0.0]);
        assert!((swept - c.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_cone_volume_swept_positive_delta() {
        let c = Cone::new(1.0, 1.0);
        let swept = c.volume_swept([1.0, 0.0, 0.0]);
        assert!(swept > c.volume());
        // Extra: base_area * 1 = π
        assert!((swept - c.volume() - PI).abs() < 1e-10);
    }

    #[test]
    fn test_cone_sdf_gradient_outside_outward() {
        let c = Cone::new(1.0, 1.0);
        let g = c.sdf_gradient([5.0, 0.0, 0.0]);
        // Should point away from cone
        assert!(g[0] > 0.0, "gx={}", g[0]);
    }

    #[test]
    fn test_cone_sdf_gradient_unit_length() {
        let c = Cone::new(1.0, 1.0);
        let g = c.sdf_gradient([3.0, 0.0, 0.0]);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-6, "gradient not unit: len={len}");
    }

    #[test]
    fn test_cone_ray_intersection_count_lateral() {
        let c = Cone::new(1.0, 2.0);
        // Ray along X at y=0: should hit lateral surface once on entry
        let count = c.ray_intersection_count([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(count >= 1, "expected >=1 hit, got {count}");
    }

    #[test]
    fn test_cone_ray_intersection_count_miss() {
        let c = Cone::new(1.0, 1.0);
        let count = c.ray_intersection_count([0.0, 10.0, 0.0], [0.0, 1.0, 0.0], 5.0);
        assert_eq!(count, 0);
    }

    #[test]
    fn test_cone_solid_angle_full_sphere() {
        // half_angle ≈ π → solid_angle ≈ 4π (full sphere)
        // half_angle = π/2 → solid_angle = 2π (hemisphere)
        let c = Cone::new(1.0, 0.5); // h=1, r=1 → half_angle = atan(1) = π/4
        let sa = c.solid_angle();
        let expected = 2.0 * PI * (1.0 - (PI / 4.0).cos());
        assert!((sa - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cone_solid_angle_positive() {
        let c = Cone::new(1.0, 1.0);
        assert!(c.solid_angle() > 0.0);
    }

    #[test]
    fn test_cone_base_rim_points_count() {
        let c = Cone::new(1.0, 2.0);
        let pts = c.base_rim_points(12);
        assert_eq!(pts.len(), 12);
    }

    #[test]
    fn test_cone_base_rim_points_on_base() {
        let c = Cone::new(1.5, 2.0);
        let pts = c.base_rim_points(16);
        for p in &pts {
            assert!((p[1] + c.half_height).abs() < 1e-12);
            let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
            assert!((xz - c.radius).abs() < 1e-9);
        }
    }

    #[test]
    fn test_cone_sphere_inside_cone_yes() {
        let c = Cone::new(2.0, 2.0);
        // Small sphere deep inside the base
        assert!(c.sphere_inside_cone([0.0, -1.9, 0.0], 0.05));
    }

    #[test]
    fn test_cone_sphere_inside_cone_no() {
        let c = Cone::new(1.0, 1.0);
        // Large sphere centered outside
        assert!(!c.sphere_inside_cone([5.0, 0.0, 0.0], 1.0));
    }

    #[test]
    fn test_cone_unfold_lateral_angle_scaling() {
        let c = Cone::new(1.0, 1.0);
        // Two points at same slant height but different u should differ in angle
        let p1 = c.unfold_lateral(0.0, 2.0);
        let p2 = c.unfold_lateral(PI / 2.0, 2.0);
        let d = ((p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2)).sqrt();
        assert!(d > 0.1, "unrolled points should differ");
    }

    #[test]
    fn test_cone_double_cone_waist_radius_zero() {
        let c = Cone::new(2.0, 1.0);
        // At y=0, double cone has zero radius (at the center)
        // contains_point_double_cone([0,0,0]) should be true (on axis)
        assert!(c.contains_point_double_cone([0.0, 0.0, 0.0]));
    }
}
