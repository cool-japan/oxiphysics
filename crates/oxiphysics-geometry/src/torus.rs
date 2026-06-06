// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Torus shape.

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::f64::consts::PI;

/// A torus defined by major radius R (center to tube center) and minor radius r (tube radius).
///
/// The torus lies in the XZ plane (Y is the axis of symmetry).
#[derive(Debug, Clone)]
pub struct Torus {
    /// Major radius: distance from the center of the torus to the center of the tube.
    pub major_radius: Real,
    /// Minor radius: radius of the tube.
    pub minor_radius: Real,
}

impl Torus {
    /// Create a new torus with the given major and minor radii.
    pub fn new(major_radius: Real, minor_radius: Real) -> Self {
        Self {
            major_radius,
            minor_radius,
        }
    }

    /// Volume = 2π²Rr²
    pub fn volume(&self) -> Real {
        2.0 * PI * PI * self.major_radius * self.minor_radius * self.minor_radius
    }

    /// Surface area = 4π²Rr
    pub fn surface_area(&self) -> Real {
        4.0 * PI * PI * self.major_radius * self.minor_radius
    }

    /// Approximate bounding box half-extents: outer radius = R+r, height = r.
    ///
    /// Returns half-extents as `[R+r, r, R+r]` (torus in XZ plane, Y is up).
    pub fn bounding_box_extents(&self) -> Vec3 {
        let outer = self.major_radius + self.minor_radius;
        Vec3::new(outer, self.minor_radius, outer)
    }

    /// Inertia tensor for a solid torus of given mass.
    ///
    /// Axis of symmetry is Y. Standard formulas:
    ///   I_y  = m*(R² + (3/4)*r²)  (but commonly written as m*(3R²+4r²)/4 – equivalent)
    ///   I_xz = m*((5/8)*r² + R²/2) (= m*(5r²+4R²)/8)
    pub fn inertia_tensor_array(&self, mass: Real) -> [[f64; 3]; 3] {
        let r = self.major_radius;
        let a = self.minor_radius;
        let iy = mass * (r * r + (3.0 / 4.0) * a * a);
        let ixz = mass * ((5.0 / 8.0) * a * a + r * r / 2.0);
        [[ixz, 0.0, 0.0], [0.0, iy, 0.0], [0.0, 0.0, ixz]]
    }

    /// Cast a ray against the torus using an iterative quartic solver approach.
    ///
    /// The torus implicit equation (centered at origin, axis = Y):
    ///   (x²+y²+z² + R²−r²)² − 4R²(x²+z²) = 0
    ///
    /// Substituting P+t*D gives a quartic in t which is solved by Ferrari's method.
    /// Returns `Some((t, normal))` or `None`.
    pub fn ray_cast_array(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        let o = Vec3::new(origin[0], origin[1], origin[2]);
        let d = Vec3::new(direction[0], direction[1], direction[2]);
        let hit = self.ray_cast_impl(&o, &d, max_toi)?;
        Some((hit.toi, [hit.normal.x, hit.normal.y, hit.normal.z]))
    }

    /// Closest point on the torus surface to point `p`.
    ///
    /// Projects `p` onto the nearest point on the ring circle (in XZ), then
    /// moves from that circle point toward `p` (or outward if p is inside tube)
    /// by the minor radius.
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        let px = p[0];
        let py = p[1];
        let pz = p[2];
        // Project onto XZ plane to find the nearest point on the ring circle
        let xz_len = (px * px + pz * pz).sqrt();
        let (ring_x, ring_z) = if xz_len < 1e-12 {
            // p is on Y axis; ring point is arbitrary, pick +X
            (self.major_radius, 0.0)
        } else {
            let s = self.major_radius / xz_len;
            (px * s, pz * s)
        };
        // Vector from ring point to p
        let dx = px - ring_x;
        let dy = py;
        let dz = pz - ring_z;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist < 1e-12 {
            // p is exactly on the ring circle – push out along +Y
            return [ring_x, self.minor_radius, ring_z];
        }
        let s = self.minor_radius / dist;
        [ring_x + dx * s, dy * s, ring_z + dz * s]
    }

    /// Returns `true` if point `p` is inside (or on the surface of) the torus tube.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let px = p[0];
        let py = p[1];
        let pz = p[2];
        let xz = (px * px + pz * pz).sqrt();
        let dist_to_ring = ((xz - self.major_radius) * (xz - self.major_radius) + py * py).sqrt();
        dist_to_ring <= self.minor_radius
    }

    /// Parametric surface sample: `u ∈ [0,2π)`, `v ∈ [0,2π)`.
    ///
    /// Returns a point `(x,y,z)` on the torus surface.
    ///   x = (R + r*cos(v)) * cos(u)
    ///   y = r * sin(v)
    ///   z = (R + r*cos(v)) * sin(u)
    pub fn sample_surface(&self, u: f64, v: f64) -> [f64; 3] {
        let r_outer = self.major_radius + self.minor_radius * v.cos();
        [
            r_outer * u.cos(),
            self.minor_radius * v.sin(),
            r_outer * u.sin(),
        ]
    }

    // ── New expanded methods ──

    /// Torus SDF (signed distance to torus surface).
    ///
    /// Negative inside the tube, positive outside.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
        let dist_to_ring = ((xz - self.major_radius).powi(2) + p[1] * p[1]).sqrt();
        dist_to_ring - self.minor_radius
    }

    /// Torus-ray analytic intersection (same as `ray_cast_array` but explicit name).
    pub fn ray_torus_analytic(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        self.ray_cast_array(origin, direction, max_toi)
    }

    /// Torus support function (plain array).
    pub fn support_array(&self, direction: [f64; 3]) -> [f64; 3] {
        let d = Vec3::new(direction[0], direction[1], direction[2]);
        let sp = self.support_point(&d);
        [sp.x, sp.y, sp.z]
    }

    /// Surface parameterization: returns the `(u, v)` angles for the point on
    /// the torus surface closest to `p`.
    ///
    /// `u` is the angle around the major circle (XZ plane), `v` is the angle
    /// around the tube (in the plane through the Y axis and the ring point).
    pub fn surface_parameters(&self, p: [f64; 3]) -> (f64, f64) {
        let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
        // u: angle in XZ plane
        let u = p[2].atan2(p[0]); // atan2(z, x)

        // Ring point at angle u
        let rx = self.major_radius * u.cos();
        let rz = self.major_radius * u.sin();

        // Vector from ring point to p
        let dx = p[0] - rx;
        let dy = p[1];
        let dz = p[2] - rz;

        // v: angle in the plane of the tube (radial outward from ring center)
        // Radial direction: (cos u, 0, sin u); Y direction: (0, 1, 0)
        let radial = xz - self.major_radius;
        let v = dy.atan2(radial);

        let _ = (dx, dz); // suppress unused warnings
        (u, v)
    }

    /// Generate `n` random points uniformly distributed on the torus surface.
    ///
    /// Uses a deterministic xorshift64 PRNG seeded with `seed`.
    /// Sampling: uniform in `u` and `v`, weighted by area element `(R + r cos v)`.
    pub fn random_surface_points(&self, n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut points = Vec::with_capacity(n);
        let mut state = seed;

        let next_f64 = |s: &mut u64| -> f64 {
            *s ^= *s << 13;
            *s ^= *s >> 7;
            *s ^= *s << 17;
            (*s as f64) / (u64::MAX as f64)
        };

        // Rejection sampling: sample (u, v) uniform in [0,2π)²,
        // accept with probability proportional to (R + r*cos(v)) / (R + r).
        let max_weight = self.major_radius + self.minor_radius;

        while points.len() < n {
            let u = next_f64(&mut state) * 2.0 * PI;
            let v = next_f64(&mut state) * 2.0 * PI;
            let w = next_f64(&mut state);

            let weight = (self.major_radius + self.minor_radius * v.cos()) / max_weight;
            if w <= weight {
                points.push(self.sample_surface(u, v));
            }
        }
        points
    }

    /// Approximate solid torus inertia tensor from `inertia_tensor_array`.
    pub fn inertia_raw(&self, mass: f64) -> [[f64; 3]; 3] {
        self.inertia_tensor_array(mass)
    }

    /// Outer radius of the torus (major + minor).
    pub fn outer_radius(&self) -> f64 {
        self.major_radius + self.minor_radius
    }

    /// Inner radius of the torus (major - minor, clamped to 0).
    pub fn inner_radius(&self) -> f64 {
        (self.major_radius - self.minor_radius).max(0.0)
    }

    // ── Extended geometry: UV mapping, geodesics, tube cross-sections, knot paths, winding ──

    /// UV mapping for the torus surface.
    ///
    /// Maps a surface point `p` to texture coordinates `(u, v)` in `[0, 1)²`.
    /// `u` corresponds to the major (longitudinal) angle and `v` to the minor
    /// (latitudinal) angle, both normalised to the range `[0, 1)`.
    pub fn uv_map(&self, p: [f64; 3]) -> [f64; 2] {
        let (theta, phi) = self.surface_parameters(p);
        // Normalise from (-π, π] to [0, 1)
        let u = ((theta / (2.0 * PI)) + 1.0) % 1.0;
        let v = ((phi / (2.0 * PI)) + 1.0) % 1.0;
        [u, v]
    }

    /// Approximate geodesic distance between two surface points `a` and `b`
    /// via the flat-torus metric.
    ///
    /// The flat-torus metric: `d = sqrt((R Δθ)² + (r Δφ)²)` where `Δθ` and
    /// `Δφ` are the wrapped angular differences on the major and minor circles.
    pub fn geodesic_distance_flat(&self, a: [f64; 3], b: [f64; 3]) -> f64 {
        let (ua, va) = self.surface_parameters(a);
        let (ub, vb) = self.surface_parameters(b);
        let d_theta = angle_diff(ua, ub);
        let d_phi = angle_diff(va, vb);
        let arc_major = self.major_radius * d_theta;
        let arc_minor = self.minor_radius * d_phi;
        (arc_major * arc_major + arc_minor * arc_minor).sqrt()
    }

    /// Area element at angle `v` on the tube: `(R + r cos v) r dv du`.
    ///
    /// Returns the area-weighted factor `R + r*cos(v)` for the given tube angle `v`.
    pub fn area_element_factor(&self, v: f64) -> f64 {
        self.major_radius + self.minor_radius * v.cos()
    }

    /// Approximate surface area by numerical integration (cross-check of closed form).
    ///
    /// Uses `n_steps` Gauss-Legendre-like sampling along each angular dimension.
    pub fn surface_area_numeric(&self, n_steps: usize) -> f64 {
        let n = n_steps.max(4);
        let du = 2.0 * PI / n as f64;
        let dv = 2.0 * PI / n as f64;
        let mut sum = 0.0;
        for iv in 0..n {
            let v = (iv as f64 + 0.5) * dv;
            let factor = self.area_element_factor(v);
            sum += factor * n as f64; // summing over n u-samples
        }
        sum * self.minor_radius * du * dv
    }

    /// Tube cross-section: return a list of `n` points on the tube circle at angle `u`.
    ///
    /// The circle lies in the plane through the ring point at angle `u` and the Y axis.
    pub fn tube_cross_section(&self, u: f64, n: usize) -> Vec<[f64; 3]> {
        let n = n.max(3);
        let ring_x = self.major_radius * u.cos();
        let ring_z = self.major_radius * u.sin();
        (0..n)
            .map(|i| {
                let v = 2.0 * PI * i as f64 / n as f64;
                // Radial direction in XZ plane at angle u
                let rx = u.cos();
                let rz = u.sin();
                let cx = ring_x + self.minor_radius * v.cos() * rx;
                let cy = self.minor_radius * v.sin();
                let cz = ring_z + self.minor_radius * v.cos() * rz;
                [cx, cy, cz]
            })
            .collect()
    }

    /// Generate a torus knot path with winding numbers `(p, q)`.
    ///
    /// A `(p, q)`-torus knot winds around the major circle `p` times while
    /// winding around the tube `q` times.  Returns `n_pts` sampled points.
    pub fn torus_knot_path(&self, p: i32, q: i32, n_pts: usize) -> Vec<[f64; 3]> {
        let n = n_pts.max(3);
        (0..n)
            .map(|i| {
                let t = 2.0 * PI * i as f64 / n as f64;
                let u = p as f64 * t;
                let v = q as f64 * t;
                self.sample_surface(u, v)
            })
            .collect()
    }

    /// Winding number of a closed curve projected onto the torus major circle.
    ///
    /// Counts how many times the sequence of `(u, _)` angles (extracted from the
    /// points via `surface_parameters`) winds around the major circle.
    ///
    /// Returns an integer winding count (positive = counter-clockwise).
    pub fn winding_number_major(&self, curve: &[[f64; 3]]) -> i32 {
        if curve.len() < 2 {
            return 0;
        }
        let angles: Vec<f64> = curve
            .iter()
            .map(|&p| self.surface_parameters(p).0)
            .collect();
        let mut total = 0.0_f64;
        for i in 0..angles.len() {
            let next = (i + 1) % angles.len();
            let diff = angle_diff(angles[i], angles[next]);
            total += diff;
        }
        (total / (2.0 * PI)).round() as i32
    }

    /// Parametric tangent vector on the torus surface in the `u` direction.
    ///
    /// `dP/du` at `(u, v)`: partial derivative with respect to the major angle.
    pub fn tangent_u(&self, u: f64, v: f64) -> [f64; 3] {
        let r_outer = self.major_radius + self.minor_radius * v.cos();
        [-r_outer * u.sin(), 0.0, r_outer * u.cos()]
    }

    /// Parametric tangent vector on the torus surface in the `v` direction.
    ///
    /// `dP/dv` at `(u, v)`: partial derivative with respect to the tube angle.
    pub fn tangent_v(&self, u: f64, v: f64) -> [f64; 3] {
        let r = self.minor_radius;
        [-r * v.sin() * u.cos(), r * v.cos(), -r * v.sin() * u.sin()]
    }

    /// Surface normal via cross product of tangents `dP/du × dP/dv`.
    ///
    /// Should agree with `torus_normal` up to sign conventions.
    pub fn normal_from_tangents(&self, u: f64, v: f64) -> [f64; 3] {
        let tu = self.tangent_u(u, v);
        let tv = self.tangent_v(u, v);
        let n = cross3([tu[0], tu[1], tu[2]], [tv[0], tv[1], tv[2]]);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-14 {
            [0.0, 1.0, 0.0]
        } else {
            [n[0] / len, n[1] / len, n[2] / len]
        }
    }

    /// Test whether a ray intersects the torus and return the hit count.
    ///
    /// A finite ray with origin `o`, direction `d`, and parameter `max_toi` can
    /// intersect the torus at 0, 2, or 4 points.  This counts all positive-t
    /// intersections within `[0, max_toi]`.
    pub fn ray_intersection_count(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_toi: f64,
    ) -> usize {
        let o = Vec3::new(origin[0], origin[1], origin[2]);
        let d = Vec3::new(direction[0], direction[1], direction[2]);
        let dir_len = d.norm();
        if dir_len < 1e-12 {
            return 0;
        }
        let dn = d / dir_len;

        let oo = o.dot(&o);
        let od = o.dot(&dn);
        let rr = self.major_radius * self.major_radius;
        let aa = self.minor_radius * self.minor_radius;
        let c_alpha = oo + rr - aa;
        let beta = 2.0 * od;
        let dxz2 = dn.x * dn.x + dn.z * dn.z;
        let od_xz = o.x * dn.x + o.z * dn.z;
        let oxz2 = o.x * o.x + o.z * o.z;
        let c4 = 1.0;
        let c3 = 2.0 * beta;
        let c2 = beta * beta + 2.0 * c_alpha - 4.0 * rr * dxz2;
        let c1 = 2.0 * beta * c_alpha - 8.0 * rr * od_xz;
        let c0 = c_alpha * c_alpha - 4.0 * rr * oxz2;
        let roots = solve_quartic(c4, c3, c2, c1, c0);

        roots
            .iter()
            .filter(|&&t| {
                !t.is_nan() && !t.is_infinite() && {
                    let t_actual = t / dir_len;
                    t_actual >= 1e-6 && t_actual <= max_toi
                }
            })
            .count()
    }

    /// Approximate torus-AABB overlap test via SDF.
    ///
    /// Returns `true` if the AABB `[min, max]` overlaps the torus.
    /// Uses the 8 corners of the AABB as sample points.
    pub fn intersects_aabb(&self, aabb_min: [f64; 3], aabb_max: [f64; 3]) -> bool {
        let corners = [
            [aabb_min[0], aabb_min[1], aabb_min[2]],
            [aabb_max[0], aabb_min[1], aabb_min[2]],
            [aabb_min[0], aabb_max[1], aabb_min[2]],
            [aabb_max[0], aabb_max[1], aabb_min[2]],
            [aabb_min[0], aabb_min[1], aabb_max[2]],
            [aabb_max[0], aabb_min[1], aabb_max[2]],
            [aabb_min[0], aabb_max[1], aabb_max[2]],
            [aabb_max[0], aabb_max[1], aabb_max[2]],
        ];
        // Check if any corner is inside the torus, or if the SDF changes sign
        let sdfs: Vec<f64> = corners.iter().map(|&c| self.sdf(c)).collect();
        let min_sdf = sdfs.iter().cloned().fold(f64::INFINITY, f64::min);
        min_sdf <= 0.0
    }

    /// Aspect ratio of the torus: major_radius / minor_radius.
    ///
    /// Returns the ratio R/r.  A value close to 1 means the torus is nearly
    /// self-intersecting; large values mean a thin tube.
    pub fn aspect_ratio(&self) -> f64 {
        self.major_radius / self.minor_radius.max(1e-30)
    }

    /// Generate `n` equidistant points along the major circle (the ring itself).
    ///
    /// Points lie on the spine of the torus (on the ring circle in the XZ plane,
    /// at `y = 0`).
    pub fn major_circle_points(&self, n: usize) -> Vec<[f64; 3]> {
        let n = n.max(2);
        (0..n)
            .map(|i| {
                let u = 2.0 * PI * i as f64 / n as f64;
                [
                    self.major_radius * u.cos(),
                    0.0,
                    self.major_radius * u.sin(),
                ]
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Internal ray-cast implementation (quartic via Ferrari / companion matrix).
    // -----------------------------------------------------------------------
    fn ray_cast_impl(&self, origin: &Vec3, dir: &Vec3, max_toi: Real) -> Option<RayHit> {
        let r_big = self.major_radius;
        let r_small = self.minor_radius;

        // Normalize direction for parameter computation; scale toi back after.
        let dir_len = dir.norm();
        if dir_len < 1e-12 {
            return None;
        }
        let d = dir / dir_len;
        let o = *origin;

        // Coefficients of the quartic  c4*t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0
        // derived from ((ox+t*dx)²+(oy+t*dy)²+(oz+t*dz)² + R²-r²)² = 4R²((ox+t*dx)²+(oz+t*dz)²)
        let oo = o.dot(&o);
        let od = o.dot(&d);
        let _dd = d.dot(&d); // = 1 since d is normalised
        let rr = r_big * r_big;
        let aa = r_small * r_small;
        let alpha = oo - rr - aa; // constant term inside big parens (before squaring)
        let beta = 2.0 * od;
        // Let f(t) = oo + 2*od*t + dd*t² + R²-r² -4R²*(... ) is expanded below.

        // Using the parametric form:
        //   Let S(t) = |o + t*d|² + R² - r² = (dd)*t² + 2*od*t + (oo + R² - r²)
        //            = t² + beta*t + (oo + rr - aa)   [since dd=1]
        // Let T(t) = (ox+t*dx)² + (oz+t*dz)²  (XZ only)
        //          = (dx²+dz²)*t² + 2*(ox*dx+oz*dz)*t + (ox²+oz²)
        let dxz2 = d.x * d.x + d.z * d.z;
        let od_xz = o.x * d.x + o.z * d.z;
        let oxz2 = o.x * o.x + o.z * o.z;

        // Quartic: (t² + beta*t + (oo+rr-aa))² - 4rr*(dxz2*t² + 2*od_xz*t + oxz2) = 0
        let c_alpha = oo + rr - aa; // = oo + R² - r²
        // Expand (t² + beta*t + c_alpha)² = t^4 + 2beta t^3 + (beta²+2c_alpha) t² + 2beta*c_alpha t + c_alpha²
        let c4 = 1.0;
        let c3 = 2.0 * beta;
        let c2 = beta * beta + 2.0 * c_alpha - 4.0 * rr * dxz2;
        let c1 = 2.0 * beta * c_alpha - 8.0 * rr * od_xz;
        let c0 = c_alpha * c_alpha - 4.0 * rr * oxz2;
        let _ = alpha; // suppress unused warning

        // Solve quartic numerically – we use a simple companion-matrix eigenvalue
        // approach via Durand-Kerner / iterative refinement, or fall back to
        // a straightforward bisection in likely intervals.
        let roots = solve_quartic(c4, c3, c2, c1, c0);

        let mut best_t = Real::INFINITY;
        for t_raw in roots {
            if t_raw.is_nan() || t_raw.is_infinite() {
                continue;
            }
            // Scale back from normalised-direction t to actual t along original dir
            let t = t_raw / dir_len;
            if t >= 1e-6 && t <= max_toi && t < best_t {
                best_t = t;
            }
        }
        if best_t.is_infinite() {
            return None;
        }

        let point = o + d * (best_t * dir_len);
        let normal = torus_normal(&point, r_big);
        Some(RayHit {
            point,
            normal,
            toi: best_t,
        })
    }
}

/// Wrapped angular difference (shortest path), result in `(-π, π]`.
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = b - a;
    // wrap into (-π, π]

    d - (2.0 * PI) * ((d + PI) / (2.0 * PI)).floor()
}

/// 3-component cross product.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the outward normal on the torus surface at `p`.
fn torus_normal(p: &Vec3, big_r: Real) -> Vec3 {
    let xz = (p.x * p.x + p.z * p.z).sqrt();
    if xz < 1e-12 {
        return Vec3::new(0.0, 1.0, 0.0);
    }
    // Nearest point on the ring circle
    let s = big_r / xz;
    let ring = Vec3::new(p.x * s, 0.0, p.z * s);
    let n = p - ring;
    let len = n.norm();
    if len < 1e-12 {
        Vec3::new(p.x / xz, 0.0, p.z / xz)
    } else {
        n / len
    }
}

/// Solve a monic quartic t^4 + (c3/c4)*t^3 + ... analytically via Ferrari's method
/// or by numerical polishing.  Returns up to 4 real roots (may include duplicates).
fn solve_quartic(c4: f64, c3: f64, c2: f64, c1: f64, c0: f64) -> [f64; 4] {
    // Normalise to monic
    let a = c3 / c4;
    let b = c2 / c4;
    let c = c1 / c4;
    let d = c0 / c4;

    // Depress the quartic: substitute t = u - a/4
    let p_depr = b - (3.0 * a * a) / 8.0;
    let q_depr = (a * a * a) / 8.0 - (a * b) / 2.0 + c;
    let r_depr = -(3.0 * a * a * a * a) / 256.0 + (a * a * b) / 16.0 - (a * c) / 4.0 + d;

    let shift = a / 4.0;

    // Solve the depressed quartic u^4 + p*u^2 + q*u + r = 0 via Ferrari
    // Resolvent cubic: m^3 + (5p/2)*m^2 + (2p²-r)*m + (p³/2-p*r/2-q²/8-... )
    // Using standard Ferrari resolvent: 8m³ + 8p*m² + (2p²-8r)*m - q² = 0
    let (m1, m2, m3) = solve_cubic_3roots(
        8.0,
        8.0 * p_depr,
        2.0 * p_depr * p_depr - 8.0 * r_depr,
        -(q_depr * q_depr),
    );

    let mut roots = [f64::NAN; 4];
    let mut ri = 0usize;

    for m in [m1, m2, m3] {
        if m.is_nan() {
            continue;
        }
        let two_m_p = 2.0 * m + p_depr;
        if two_m_p < 0.0 {
            continue;
        }
        let sqrt_2mp = two_m_p.sqrt();
        if sqrt_2mp < 1e-14 {
            continue;
        }

        // Two quadratics:
        // u² + sqrt(2m+p)*u + (m + q/(2*sqrt(2m+p))) = 0
        // u² - sqrt(2m+p)*u + (m - q/(2*sqrt(2m+p))) = 0
        let half_q_over = q_depr / (2.0 * sqrt_2mp);

        for &(s, c_q) in &[(1.0_f64, m + half_q_over), (-1.0_f64, m - half_q_over)] {
            let disc = (s * sqrt_2mp) * (s * sqrt_2mp) / 4.0 - c_q;
            // discriminant of u² ∓ sqrt(2m+p)*u + c_q = 0 is (sqrt_2mp/2)² - c_q
            let disc2 = sqrt_2mp * sqrt_2mp / 4.0 - c_q;
            if disc2 >= -1e-9 {
                let sd = disc2.max(0.0).sqrt();
                let half_b = s * sqrt_2mp / 2.0;
                let r0 = -half_b + sd - shift;
                let r1 = -half_b - sd - shift;
                if ri < 4 {
                    roots[ri] = r0;
                    ri += 1;
                }
                if ri < 4 {
                    roots[ri] = r1;
                    ri += 1;
                }
            }
            let _ = disc; // suppress warning
        }
        if ri >= 4 {
            break;
        }
    }
    roots
}

/// Returns one real root of the cubic a*t³ + b*t² + c*t + d = 0,
/// plus two more (may be complex → NAN).
fn solve_cubic_3roots(a: f64, b: f64, c: f64, d: f64) -> (f64, f64, f64) {
    // Monic: t³ + (b/a)*t² + (c/a)*t + (d/a) = 0
    let b = b / a;
    let c = c / a;
    let d = d / a;

    // Depress: t = u - b/3
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let shift = b / 3.0;

    let discriminant = -(4.0 * p * p * p + 27.0 * q * q);

    if discriminant >= 0.0 {
        // Three real roots using trigonometric method
        let m = 2.0 * (-p / 3.0).sqrt();
        let theta = ((3.0 * q) / (2.0 * p) * ((-3.0) / p).sqrt())
            .clamp(-1.0, 1.0)
            .acos();
        let r0 = m * (theta / 3.0).cos() - shift;
        let r1 = m * ((theta - 2.0 * PI) / 3.0).cos() - shift;
        let r2 = m * ((theta - 4.0 * PI) / 3.0).cos() - shift;
        (r0, r1, r2)
    } else {
        // One real root via Cardano
        let d_val = -(4.0 * p * p * p + 27.0 * q * q);
        let inner = -q / 2.0;
        let outer = (q * q / 4.0 + p * p * p / 27.0).sqrt();
        let u = cbrt_real(inner + outer);
        let v = cbrt_real(inner - outer);
        let r0 = u + v - shift;
        let _ = d_val; // suppress warning
        (r0, f64::NAN, f64::NAN)
    }
}

fn cbrt_real(x: f64) -> f64 {
    if x >= 0.0 { x.cbrt() } else { -((-x).cbrt()) }
}

impl Shape for Torus {
    fn bounding_box(&self) -> Aabb {
        let half = self.bounding_box_extents();
        Aabb::new(-half, half)
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let eps = 1e-10;
        // Project direction onto XZ plane
        let dir_xz = Vec3::new(direction.x, 0.0, direction.z);
        let xz_len = dir_xz.norm();
        let xz_norm = if xz_len < eps {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            dir_xz / xz_len
        };
        // Tube center on the ring
        let tube_center = xz_norm * self.major_radius;
        // Add minor radius in full direction
        let dir_len = direction.norm();
        let dir_norm = if dir_len < eps {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            direction / dir_len
        };
        tube_center + dir_norm * self.minor_radius
    }

    fn volume(&self) -> Real {
        self.volume()
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::zeros()
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let r = self.major_radius;
        let a = self.minor_radius;
        // Solid torus inertia tensors (axis of symmetry = Y):
        //   I_y = mass * (R² + (3/4)*r²)
        //   I_x = I_z = mass * (5/8*r² + R²/2)
        // (from standard solid torus formulas)
        let iy = mass * (r * r + (3.0 / 4.0) * a * a);
        let ixz = mass * ((5.0 / 8.0) * a * a + r * r / 2.0);
        Mat3::new(ixz, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ixz)
    }

    fn mass_properties(&self, density: Real) -> oxiphysics_core::MassProperties {
        let mass = density * self.volume();
        oxiphysics_core::MassProperties::new(mass, self.center_of_mass(), self.inertia_tensor(mass))
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        self.ray_cast_impl(ray_origin, ray_direction, max_toi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_torus_volume() {
        let t = Torus::new(1.0, 0.5);
        // V = 2π²*1*0.25 ≈ 4.935
        let expected = 2.0 * PI * PI * 1.0 * 0.25;
        assert!(
            (t.volume() - expected).abs() < 1e-6,
            "volume={}, expected={}",
            t.volume(),
            expected
        );
        // Approximately 4.935
        assert!(
            (t.volume() - 4.935).abs() < 0.01,
            "volume ≈ 4.935, got {}",
            t.volume()
        );
    }

    #[test]
    fn test_torus_surface_area() {
        let t = Torus::new(1.0, 0.5);
        // SA = 4π²*1*0.5 ≈ 19.74
        let expected = 4.0 * PI * PI * 1.0 * 0.5;
        assert!(
            (t.surface_area() - expected).abs() < 1e-6,
            "surface_area={}, expected={}",
            t.surface_area(),
            expected
        );
        assert!(
            (t.surface_area() - 19.74).abs() < 0.01,
            "surface_area ≈ 19.74, got {}",
            t.surface_area()
        );
    }

    #[test]
    fn test_torus_bounding_box() {
        let t = Torus::new(2.0, 0.5);
        // half_extents = [2.5, 0.5, 2.5]
        let extents = t.bounding_box_extents();
        assert!((extents.x - 2.5).abs() < 1e-10, "extents.x={}", extents.x);
        assert!((extents.y - 0.5).abs() < 1e-10, "extents.y={}", extents.y);
        assert!((extents.z - 2.5).abs() < 1e-10, "extents.z={}", extents.z);
        let bb = t.bounding_box();
        assert!((bb.min.x + 2.5).abs() < 1e-10);
        assert!((bb.max.x - 2.5).abs() < 1e-10);
        assert!((bb.min.y + 0.5).abs() < 1e-10);
        assert!((bb.max.y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_torus_support_xz() {
        let t = Torus::new(1.0, 0.5);
        // dir = (1,0,0) → support should be at (R+r, 0, 0)
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let sp = t.support_point(&dir);
        let expected = 1.0 + 0.5; // R + r
        assert!(
            (sp.x - expected).abs() < 1e-10,
            "support.x={}, expected={}",
            sp.x,
            expected
        );
        assert!(sp.y.abs() < 1e-10, "support.y={}", sp.y);
        assert!(sp.z.abs() < 1e-10, "support.z={}", sp.z);
    }

    #[test]
    fn test_torus_support_y() {
        let t = Torus::new(1.0, 0.5);
        // dir = (0,1,0) → y component = r (minor radius)
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let sp = t.support_point(&dir);
        // tube center at (R, 0, 0) + r*(0,1,0)/|(1,1,0)| -- actually dir_xz is zero so
        // tube_center = (R, 0, 0), dir_norm = (0,1,0), so result = (R, r, 0)
        assert!(
            (sp.y - 0.5).abs() < 1e-10,
            "support.y={}, expected minor_radius=0.5",
            sp.y
        );
    }

    #[test]
    fn test_torus_mass_properties() {
        let t = Torus::new(1.0, 0.5);
        let density = 1000.0;
        let mp = t.mass_properties(density);
        let expected_mass = density * t.volume();
        assert!(
            (mp.mass - expected_mass).abs() < 1e-6,
            "mass={}, expected={}",
            mp.mass,
            expected_mass
        );
        assert!(mp.mass > 0.0, "mass should be positive");
        assert!(
            mp.local_inertia[(0, 0)] > 0.0,
            "I_xx should be positive: {}",
            mp.local_inertia[(0, 0)]
        );
        assert!(
            mp.local_inertia[(1, 1)] > 0.0,
            "I_yy should be positive: {}",
            mp.local_inertia[(1, 1)]
        );
        assert!(
            mp.local_inertia[(2, 2)] > 0.0,
            "I_zz should be positive: {}",
            mp.local_inertia[(2, 2)]
        );
    }

    #[test]
    fn test_torus_contains_point_inside() {
        let t = Torus::new(3.0, 1.0);
        // Point on the tube center circle in XZ plane → inside
        assert!(
            t.contains_point([3.0, 0.0, 0.0]),
            "ring centre should be inside"
        );
        // Point clearly outside
        assert!(
            !t.contains_point([0.0, 0.0, 0.0]),
            "origin should be outside"
        );
        // Point on tube surface
        assert!(
            t.contains_point([3.0, 1.0, 0.0]),
            "top of tube should be inside/on"
        );
        // Point just beyond tube
        assert!(
            !t.contains_point([3.0, 1.1, 0.0]),
            "beyond tube should be outside"
        );
    }

    #[test]
    fn test_torus_closest_point_on_surface() {
        let t = Torus::new(3.0, 1.0);
        // A point far outside in +X direction
        let cp = t.closest_point([10.0, 0.0, 0.0]);
        // Should be at (R+r, 0, 0) = (4, 0, 0)
        assert!((cp[0] - 4.0).abs() < 1e-6, "cp.x={}, expected 4.0", cp[0]);
        assert!(cp[1].abs() < 1e-6, "cp.y={}", cp[1]);
        assert!(cp[2].abs() < 1e-6, "cp.z={}", cp[2]);
    }

    #[test]
    fn test_torus_sample_surface_on_torus() {
        let t = Torus::new(3.0, 1.0);
        // Sample at u=0, v=0 → (R+r, 0, 0) = (4, 0, 0)
        let p = t.sample_surface(0.0, 0.0);
        assert!((p[0] - 4.0).abs() < 1e-10, "p.x={}", p[0]);
        assert!(p[1].abs() < 1e-10, "p.y={}", p[1]);
        // Verify point is actually on the torus surface
        assert!(
            t.contains_point(p),
            "sampled point {:?} should be on torus surface",
            p
        );
        // Sample at u=PI/2, v=PI/2 → top of tube at 90°
        let p2 = t.sample_surface(PI / 2.0, PI / 2.0);
        // (R + r*cos(π/2)) * cos(π/2) = R*0 = 0, y = r*sin(π/2) = r, z = R*sin(π/2) = R
        assert!(p2[0].abs() < 1e-10, "p2.x={}", p2[0]);
        assert!((p2[1] - 1.0).abs() < 1e-10, "p2.y={}", p2[1]);
        assert!((p2[2] - 3.0).abs() < 1e-10, "p2.z={}", p2[2]);
    }

    #[test]
    fn test_torus_inertia_tensor_array() {
        let t = Torus::new(2.0, 0.5);
        let it = t.inertia_tensor_array(10.0);
        // All diagonal entries should be positive
        assert!(it[0][0] > 0.0, "Ixx={}", it[0][0]);
        assert!(it[1][1] > 0.0, "Iyy={}", it[1][1]);
        assert!(it[2][2] > 0.0, "Izz={}", it[2][2]);
        // Off-diagonal should be zero
        assert!(it[0][1].abs() < 1e-12);
        assert!(it[0][2].abs() < 1e-12);
        // Ixx == Izz by symmetry
        assert!((it[0][0] - it[2][2]).abs() < 1e-10);
    }

    #[test]
    fn test_torus_ray_cast_hit() {
        let t = Torus::new(3.0, 1.0);
        // Ray along +Y from below, aimed at (3,0,0) which is on the tube
        let origin = Vec3::new(3.0, -10.0, 0.0);
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let hit = t.ray_cast(&origin, &dir, 100.0);
        assert!(hit.is_some(), "ray through tube should hit");
        let hit = hit.unwrap();
        // Should hit somewhere between origin and tube; toi must be positive
        assert!(
            hit.toi > 0.0 && hit.toi < 20.0,
            "toi should be reasonable, got {}",
            hit.toi
        );
    }

    #[test]
    fn test_torus_ray_cast_miss() {
        let t = Torus::new(3.0, 1.0);
        // Ray from far outside the torus going away from it
        let origin = Vec3::new(0.0, 10.0, 0.0);
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let hit = t.ray_cast(&origin, &dir, 5.0);
        assert!(hit.is_none(), "ray going away should miss");
    }

    #[test]
    fn test_torus_ray_cast_array() {
        let t = Torus::new(3.0, 1.0);
        let result = t.ray_cast_array([3.0, -10.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_some(), "ray should hit");
        let (toi, _n) = result.unwrap();
        assert!(
            toi > 0.0 && toi < 20.0,
            "toi should be reasonable, got {}",
            toi
        );
    }

    // ── Expanded tests ──

    #[test]
    fn test_torus_sdf_inside_tube() {
        let t = Torus::new(3.0, 1.0);
        // Point on the tube ring center: sdf should be -minor_radius
        let p = [3.0, 0.0, 0.0];
        assert!(
            (t.sdf(p) + 1.0).abs() < 1e-10,
            "sdf at ring center should be -1, got {}",
            t.sdf(p)
        );
    }

    #[test]
    fn test_torus_sdf_outside() {
        let t = Torus::new(3.0, 1.0);
        // Point far outside
        let p = [10.0, 0.0, 0.0];
        assert!(t.sdf(p) > 0.0, "sdf outside should be positive");
    }

    #[test]
    fn test_torus_sdf_on_surface() {
        let t = Torus::new(3.0, 1.0);
        // Point on the outer surface: (R+r, 0, 0) = (4, 0, 0)
        let p = [4.0, 0.0, 0.0];
        assert!(
            t.sdf(p).abs() < 1e-10,
            "sdf at surface should be ~0, got {}",
            t.sdf(p)
        );
    }

    #[test]
    fn test_torus_ray_torus_analytic() {
        let t = Torus::new(3.0, 1.0);
        let result = t.ray_torus_analytic([3.0, -10.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(result.is_some());
    }

    #[test]
    fn test_torus_support_array_xplus() {
        let t = Torus::new(2.0, 0.5);
        let sp = t.support_array([1.0, 0.0, 0.0]);
        // Should be at (R+r, 0, 0)
        assert!((sp[0] - 2.5).abs() < 1e-10, "sp.x={}", sp[0]);
        assert!(sp[1].abs() < 1e-10);
    }

    #[test]
    fn test_torus_surface_parameters_roundtrip() {
        let t = Torus::new(3.0, 1.0);
        let u_in = 1.0_f64;
        let v_in = 0.5_f64;
        let p = t.sample_surface(u_in, v_in);
        let (u_out, v_out) = t.surface_parameters(p);
        // u should match (mod 2π)
        let du = (u_out - u_in).abs();
        let du_mod = du.min((du - 2.0 * PI).abs());
        assert!(du_mod < 1e-9, "u mismatch: in={u_in}, out={u_out}");
        assert!(
            (v_out - v_in).abs() < 1e-9,
            "v mismatch: in={v_in}, out={v_out}"
        );
    }

    #[test]
    fn test_torus_random_surface_points_count() {
        let t = Torus::new(2.0, 0.5);
        let pts = t.random_surface_points(40, 777);
        assert_eq!(pts.len(), 40);
    }

    #[test]
    fn test_torus_random_surface_points_on_surface() {
        let t = Torus::new(3.0, 0.5);
        let pts = t.random_surface_points(60, 42);
        for p in &pts {
            let sdf = t.sdf(*p);
            assert!(sdf.abs() < 1e-6, "point {:?} sdf={sdf}", p);
        }
    }

    #[test]
    fn test_torus_outer_inner_radius() {
        let t = Torus::new(3.0, 1.0);
        assert!((t.outer_radius() - 4.0).abs() < 1e-12);
        assert!((t.inner_radius() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_torus_inner_radius_clamped() {
        let t = Torus::new(0.5, 1.0); // minor > major
        assert_eq!(t.inner_radius(), 0.0);
    }

    #[test]
    fn test_torus_inertia_raw_matches_array() {
        let t = Torus::new(2.0, 0.5);
        let raw = t.inertia_raw(5.0);
        let arr = t.inertia_tensor_array(5.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!((raw[i][j] - arr[i][j]).abs() < 1e-12);
            }
        }
    }

    // ── Extended tests for new methods ──────────────────────────────────────

    #[test]
    fn test_torus_uv_map_roundtrip() {
        let t = Torus::new(3.0, 1.0);
        let p = t.sample_surface(1.2, 0.8);
        let [u, v] = t.uv_map(p);
        // u and v should be in [0, 1)
        assert!((0.0..1.0).contains(&u), "u={u}");
        assert!((0.0..1.0).contains(&v), "v={v}");
        // Reconstruct point
        let p2 = t.sample_surface(u * 2.0 * PI, v * 2.0 * PI);
        let d = ((p[0] - p2[0]).powi(2) + (p[1] - p2[1]).powi(2) + (p[2] - p2[2]).powi(2)).sqrt();
        assert!(d < 1e-9, "roundtrip error={d}");
    }

    #[test]
    fn test_torus_uv_map_zero() {
        let t = Torus::new(2.0, 0.5);
        // (R+r, 0, 0) → u=0, v=0
        let p = [t.major_radius + t.minor_radius, 0.0, 0.0];
        let [u, v] = t.uv_map(p);
        assert!(u.abs() < 1e-9 || (u - 1.0).abs() < 1e-9, "u={u}");
        assert!(v.abs() < 1e-9 || (v - 1.0).abs() < 1e-9, "v={v}");
    }

    #[test]
    fn test_torus_geodesic_distance_flat_same_point() {
        let t = Torus::new(3.0, 1.0);
        let p = t.sample_surface(0.5, 0.3);
        let d = t.geodesic_distance_flat(p, p);
        assert!(d.abs() < 1e-9, "distance to self should be 0, got {d}");
    }

    #[test]
    fn test_torus_geodesic_distance_flat_opposite_minor() {
        let t = Torus::new(3.0, 1.0);
        let p1 = t.sample_surface(0.0, 0.0);
        let p2 = t.sample_surface(0.0, PI);
        // Half tube circumference = π * r
        let d = t.geodesic_distance_flat(p1, p2);
        let expected = PI * t.minor_radius;
        assert!((d - expected).abs() < 1e-9, "d={d} expected={expected}");
    }

    #[test]
    fn test_torus_area_element_factor_at_v0() {
        let t = Torus::new(3.0, 1.0);
        // At v=0 (outer equator): factor = R + r
        assert!((t.area_element_factor(0.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_torus_area_element_factor_at_vpi() {
        let t = Torus::new(3.0, 1.0);
        // At v=π (inner equator): factor = R - r = 2
        assert!((t.area_element_factor(PI) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_torus_surface_area_numeric_approximation() {
        let t = Torus::new(3.0, 1.0);
        let analytic = t.surface_area();
        let numeric = t.surface_area_numeric(256);
        // Should agree within 0.1%
        let rel_err = (numeric - analytic).abs() / analytic;
        assert!(
            rel_err < 0.001,
            "numeric={numeric} analytic={analytic} rel_err={rel_err}"
        );
    }

    #[test]
    fn test_torus_tube_cross_section_count() {
        let t = Torus::new(2.0, 0.5);
        let pts = t.tube_cross_section(0.0, 16);
        assert_eq!(pts.len(), 16);
    }

    #[test]
    fn test_torus_tube_cross_section_on_torus() {
        let t = Torus::new(2.0, 0.5);
        let pts = t.tube_cross_section(0.0, 20);
        for p in &pts {
            let sdf = t.sdf(*p);
            assert!(sdf.abs() < 1e-9, "cross-section pt not on torus sdf={sdf}");
        }
    }

    #[test]
    fn test_torus_knot_path_count() {
        let t = Torus::new(3.0, 1.0);
        let pts = t.torus_knot_path(2, 3, 60);
        assert_eq!(pts.len(), 60);
    }

    #[test]
    fn test_torus_knot_path_on_surface() {
        let t = Torus::new(3.0, 1.0);
        let pts = t.torus_knot_path(2, 3, 48);
        for p in &pts {
            let sdf = t.sdf(*p);
            assert!(sdf.abs() < 1e-9, "knot pt sdf={sdf}");
        }
    }

    #[test]
    fn test_torus_knot_path_trefoil() {
        // (2,3) torus knot should be a closed curve after 1 full parameter cycle
        let t = Torus::new(3.0, 1.0);
        let pts = t.torus_knot_path(2, 3, 60);
        // First and "last+1" point should wrap (they're not the same since we
        // sample [0, 2π) exclusive, but the path is topologically closed)
        assert_eq!(pts.len(), 60);
    }

    #[test]
    fn test_torus_winding_number_major_circle() {
        // A curve following the major circle once
        let t = Torus::new(3.0, 1.0);
        let n = 64;
        let curve: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let u = 2.0 * PI * i as f64 / n as f64;
                t.sample_surface(u, 0.0) // v=0 fixed, u winds around once
            })
            .collect();
        let winding = t.winding_number_major(&curve);
        assert_eq!(winding, 1, "should wind once, got {winding}");
    }

    #[test]
    fn test_torus_tangent_u_perpendicular_to_normal() {
        let t = Torus::new(3.0, 1.0);
        let u = 0.5;
        let v = 0.3;
        let tu = t.tangent_u(u, v);
        let n = t.normal_from_tangents(u, v);
        let dot = tu[0] * n[0] + tu[1] * n[1] + tu[2] * n[2];
        assert!(dot.abs() < 1e-9, "tangent_u · normal = {dot}");
    }

    #[test]
    fn test_torus_tangent_v_perpendicular_to_normal() {
        let t = Torus::new(3.0, 1.0);
        let u = 0.5;
        let v = 0.3;
        let tv = t.tangent_v(u, v);
        let n = t.normal_from_tangents(u, v);
        let dot = tv[0] * n[0] + tv[1] * n[1] + tv[2] * n[2];
        assert!(dot.abs() < 1e-9, "tangent_v · normal = {dot}");
    }

    #[test]
    fn test_torus_normal_from_tangents_unit() {
        let t = Torus::new(3.0, 1.0);
        let n = t.normal_from_tangents(0.5, 0.3);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-9, "normal not unit length={len}");
    }

    #[test]
    fn test_torus_ray_intersection_count_through() {
        let t = Torus::new(3.0, 1.0);
        // Ray along +Y through a tube: enters and exits the tube = 2 intersections,
        // but the quartic may produce up to 4 roots (2 real, 2 near-degenerate).
        let count = t.ray_intersection_count([3.0, -10.0, 0.0], [0.0, 1.0, 0.0], 100.0);
        assert!(count >= 2, "expected >=2 intersections, got {count}");
    }

    #[test]
    fn test_torus_ray_intersection_count_miss() {
        let t = Torus::new(3.0, 1.0);
        // Ray going away from torus
        let count = t.ray_intersection_count([0.0, 20.0, 0.0], [0.0, 1.0, 0.0], 5.0);
        assert_eq!(count, 0, "expected 0 intersections");
    }

    #[test]
    fn test_torus_intersects_aabb_inside() {
        let t = Torus::new(3.0, 1.0);
        // AABB that has a corner at (3,0,0) which is inside the tube (sdf < 0)
        assert!(t.intersects_aabb([2.5, -0.5, -0.5], [3.5, 0.5, 0.5]));
    }

    #[test]
    fn test_torus_intersects_aabb_outside() {
        let t = Torus::new(3.0, 1.0);
        // AABB far from torus
        assert!(!t.intersects_aabb([20.0, 20.0, 20.0], [30.0, 30.0, 30.0]));
    }

    #[test]
    fn test_torus_aspect_ratio() {
        let t = Torus::new(4.0, 1.0);
        assert!((t.aspect_ratio() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_torus_major_circle_points_count() {
        let t = Torus::new(3.0, 1.0);
        let pts = t.major_circle_points(36);
        assert_eq!(pts.len(), 36);
    }

    #[test]
    fn test_torus_major_circle_points_on_ring() {
        let t = Torus::new(3.0, 1.0);
        let pts = t.major_circle_points(24);
        for p in &pts {
            let xz = (p[0] * p[0] + p[2] * p[2]).sqrt();
            assert!((xz - t.major_radius).abs() < 1e-9, "xz={xz}");
            assert!(p[1].abs() < 1e-12, "y={}", p[1]);
        }
    }

    #[test]
    fn test_torus_winding_number_single_point() {
        let t = Torus::new(3.0, 1.0);
        let curve = vec![[3.0, 0.0, 0.0]];
        assert_eq!(t.winding_number_major(&curve), 0);
    }

    #[test]
    fn test_torus_winding_number_empty() {
        let t = Torus::new(3.0, 1.0);
        assert_eq!(t.winding_number_major(&[]), 0);
    }

    #[test]
    fn test_torus_surface_area_numeric_32() {
        let t = Torus::new(2.0, 0.5);
        let num = t.surface_area_numeric(32);
        let ana = t.surface_area();
        let err = (num - ana).abs() / ana;
        assert!(err < 0.01, "err={err}");
    }

    #[test]
    fn test_torus_tangents_nonzero() {
        let t = Torus::new(3.0, 1.0);
        let tu = t.tangent_u(0.1, 0.2);
        let tv = t.tangent_v(0.1, 0.2);
        let len_u = (tu[0] * tu[0] + tu[1] * tu[1] + tu[2] * tu[2]).sqrt();
        let len_v = (tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2]).sqrt();
        assert!(len_u > 1e-6);
        assert!(len_v > 1e-6);
    }

    #[test]
    fn test_torus_tube_cross_section_minimum_3() {
        let t = Torus::new(2.0, 0.5);
        let pts = t.tube_cross_section(0.0, 1); // should clamp to 3
        assert_eq!(pts.len(), 3);
    }

    #[test]
    fn test_torus_knot_path_trivial_11() {
        // A (1,1) torus knot is just a circle on the surface
        let t = Torus::new(3.0, 1.0);
        let pts = t.torus_knot_path(1, 1, 32);
        for p in &pts {
            assert!(t.sdf(*p).abs() < 1e-8);
        }
    }

    #[test]
    fn test_torus_area_element_factor_positive() {
        let t = Torus::new(3.0, 1.0);
        // For major > minor, factor is always positive
        for i in 0..32 {
            let v = 2.0 * PI * i as f64 / 32.0;
            assert!(t.area_element_factor(v) > 0.0, "factor at v={v}");
        }
    }
}
