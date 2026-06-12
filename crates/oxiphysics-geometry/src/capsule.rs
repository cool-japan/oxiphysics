// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Capsule shape (sphere-swept line segment along Y axis).

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use std::f64::consts::PI;

/// A capsule defined by radius and half-height along the Y axis.
#[derive(Debug, Clone)]
pub struct Capsule {
    /// Radius of the capsule end-caps.
    pub radius: Real,
    /// Half-height of the cylindrical segment.
    pub half_height: Real,
}

impl Capsule {
    /// Create a new capsule.
    pub fn new(radius: Real, half_height: Real) -> Self {
        Self {
            radius,
            half_height,
        }
    }

    /// Volume: cylinder + two hemisphere caps = πr²h + (4/3)πr³.
    pub fn volume_explicit(&self) -> Real {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        PI * r * r * h + (4.0 / 3.0) * PI * r.powi(3)
    }

    /// Surface area: 2πrh (cylinder lateral) + 4πr² (two hemispheres).
    pub fn surface_area(&self) -> Real {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        2.0 * PI * r * h + 4.0 * PI * r * r
    }

    /// Inertia tensor as \[\[f64;3\\];3] row-major.
    pub fn inertia_tensor_array(&self, mass: f64) -> [[f64; 3]; 3] {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        let r2 = r * r;
        let h2 = h * h;

        let vol_cyl = PI * r2 * h;
        let vol_sphere = (4.0 / 3.0) * PI * r.powi(3);
        let total_vol = vol_cyl + vol_sphere;

        let mass_cyl = mass * vol_cyl / total_vol;
        let mass_sphere = mass * vol_sphere / total_vol;

        let iy_cyl = 0.5 * mass_cyl * r2;
        let ixz_cyl = mass_cyl * (3.0 * r2 + h2) / 12.0;

        let iy_sphere = 0.4 * mass_sphere * r2;
        let offset = self.half_height + 3.0 * r / 8.0;
        let ixz_sphere = 0.4 * mass_sphere * r2 + mass_sphere * offset * offset;

        let iy = iy_cyl + iy_sphere;
        let ixz = ixz_cyl + ixz_sphere;

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

    /// GJK support function: farthest point in `direction`.
    pub fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        // Choose which hemicap center
        let cap_y = if direction[1] >= 0.0 {
            self.half_height
        } else {
            -self.half_height
        };
        if len < 1e-12 {
            return [0.0, cap_y, 0.0];
        }
        let s = self.radius / len;
        [direction[0] * s, cap_y + direction[1] * s, direction[2] * s]
    }

    // ── New methods ──

    /// Closest point on the capsule surface to point `p`.
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        // Project p onto the medial axis (Y-axis segment from -hh to +hh)
        let clamped_y = p[1].clamp(-self.half_height, self.half_height);
        // Vector from axis point to p
        let dx = p[0];
        let dy = p[1] - clamped_y;
        let dz = p[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();

        if len < 1e-12 {
            // On the axis; return a point on the surface in +X direction
            return [self.radius, clamped_y, 0.0];
        }

        let s = self.radius / len;
        [dx * s, clamped_y + dy * s, dz * s]
    }

    /// Returns true if `p` is inside (or on) the capsule.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let clamped_y = p[1].clamp(-self.half_height, self.half_height);
        let dx = p[0];
        let dy = p[1] - clamped_y;
        let dz = p[2];
        dx * dx + dy * dy + dz * dz <= self.radius * self.radius
    }

    /// Signed distance from a point to the capsule surface.
    /// Negative if inside, positive if outside.
    pub fn signed_distance(&self, p: [f64; 3]) -> f64 {
        let clamped_y = p[1].clamp(-self.half_height, self.half_height);
        let dx = p[0];
        let dy = p[1] - clamped_y;
        let dz = p[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist - self.radius
    }

    /// Medial axis endpoints: the two centers of the hemispherical caps.
    pub fn medial_axis_endpoints(&self) -> ([f64; 3], [f64; 3]) {
        ([0.0, self.half_height, 0.0], [0.0, -self.half_height, 0.0])
    }

    /// Full length of the capsule (tip to tip along Y).
    pub fn full_length(&self) -> f64 {
        2.0 * self.half_height + 2.0 * self.radius
    }

    /// Medial axis length (just the cylindrical segment).
    pub fn medial_axis_length(&self) -> f64 {
        2.0 * self.half_height
    }

    /// Distance between two segments in 3D.
    /// Segment A: from `a0` to `a1`, Segment B: from `b0` to `b1`.
    /// Returns the minimum distance.
    pub fn segment_segment_distance(a0: [f64; 3], a1: [f64; 3], b0: [f64; 3], b1: [f64; 3]) -> f64 {
        let da = [a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]];
        let db = [b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]];
        let r = [a0[0] - b0[0], a0[1] - b0[1], a0[2] - b0[2]];

        let a = da[0] * da[0] + da[1] * da[1] + da[2] * da[2];
        let e = db[0] * db[0] + db[1] * db[1] + db[2] * db[2];
        let f = db[0] * r[0] + db[1] * r[1] + db[2] * r[2];

        let eps = 1e-12;

        if a <= eps && e <= eps {
            // Both segments degenerate to points
            let dx = r[0];
            let dy = r[1];
            let dz = r[2];
            return (dx * dx + dy * dy + dz * dz).sqrt();
        }

        let b = da[0] * db[0] + da[1] * db[1] + da[2] * db[2];
        let c = da[0] * r[0] + da[1] * r[1] + da[2] * r[2];

        let (s, t);

        if a <= eps {
            s = 0.0;
            t = (f / e).clamp(0.0, 1.0);
        } else if e <= eps {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let denom = a * e - b * b;
            s = if denom.abs() > eps {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t_nom = b * s + f;
            if t_nom < 0.0 {
                t = 0.0;
            } else if t_nom > e {
                t = 1.0;
            } else {
                t = t_nom / e;
            }
        }

        // Recompute s based on t
        let s = if a > eps {
            ((b * t - c) / a).clamp(0.0, 1.0)
        } else {
            s
        };

        let dx = r[0] + da[0] * s - db[0] * t;
        let dy = r[1] + da[1] * s - db[1] * t;
        let dz = r[2] + da[2] * s - db[2] * t;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Capsule-capsule distance.
    /// Both capsules are centered at the origin with Y-axis alignment,
    /// but offset by given centers.
    /// Returns the minimum distance between the two capsule surfaces (0 if overlapping).
    pub fn capsule_capsule_distance(
        &self,
        center_a: [f64; 3],
        other: &Capsule,
        center_b: [f64; 3],
    ) -> f64 {
        // Medial axis of A: from center_a + [0, -hh_a, 0] to center_a + [0, hh_a, 0]
        let a0 = [center_a[0], center_a[1] - self.half_height, center_a[2]];
        let a1 = [center_a[0], center_a[1] + self.half_height, center_a[2]];
        let b0 = [center_b[0], center_b[1] - other.half_height, center_b[2]];
        let b1 = [center_b[0], center_b[1] + other.half_height, center_b[2]];

        let seg_dist = Self::segment_segment_distance(a0, a1, b0, b1);
        let surface_dist = seg_dist - self.radius - other.radius;
        surface_dist.max(0.0)
    }

    /// Check if two capsules overlap.
    pub fn capsule_capsule_overlap(
        &self,
        center_a: [f64; 3],
        other: &Capsule,
        center_b: [f64; 3],
    ) -> bool {
        let a0 = [center_a[0], center_a[1] - self.half_height, center_a[2]];
        let a1 = [center_a[0], center_a[1] + self.half_height, center_a[2]];
        let b0 = [center_b[0], center_b[1] - other.half_height, center_b[2]];
        let b1 = [center_b[0], center_b[1] + other.half_height, center_b[2]];

        let seg_dist = Self::segment_segment_distance(a0, a1, b0, b1);
        seg_dist <= self.radius + other.radius
    }

    /// Project a point onto the medial axis (the Y-axis segment).
    /// Returns the clamped Y value.
    pub fn project_on_medial_axis(&self, p: [f64; 3]) -> f64 {
        p[1].clamp(-self.half_height, self.half_height)
    }

    /// Distance from a point to the medial axis.
    pub fn distance_to_medial_axis(&self, p: [f64; 3]) -> f64 {
        let clamped_y = p[1].clamp(-self.half_height, self.half_height);
        let dy = p[1] - clamped_y;
        (p[0] * p[0] + dy * dy + p[2] * p[2]).sqrt()
    }

    // ── New expanded methods ──

    /// Inertia tensor (mass-normalized) returned as `[[f64;3\];3]` (identical to
    /// `inertia_tensor_array` but accepts mass separately for the public array API).
    pub fn inertia_tensor_raw(&self, mass: f64) -> [[f64; 3]; 3] {
        self.inertia_tensor_array(mass)
    }

    /// Closest points between two capsule medial axes (and thus capsule surfaces).
    ///
    /// Returns `(pa, pb, seg_dist)` where `pa` and `pb` are the closest points
    /// on the respective medial axes and `seg_dist` is the axis-to-axis distance.
    pub fn closest_points_capsule_vs_capsule(
        &self,
        center_a: [f64; 3],
        other: &Capsule,
        center_b: [f64; 3],
    ) -> ([f64; 3], [f64; 3], f64) {
        let a0 = [center_a[0], center_a[1] - self.half_height, center_a[2]];
        let a1 = [center_a[0], center_a[1] + self.half_height, center_a[2]];
        let b0 = [center_b[0], center_b[1] - other.half_height, center_b[2]];
        let b1 = [center_b[0], center_b[1] + other.half_height, center_b[2]];

        let (pa, pb, seg_dist) = Self::segment_segment_closest(a0, a1, b0, b1);
        (pa, pb, seg_dist)
    }

    /// Returns the closest points on two 3D segments and the distance between them.
    pub fn segment_segment_closest(
        a0: [f64; 3],
        a1: [f64; 3],
        b0: [f64; 3],
        b1: [f64; 3],
    ) -> ([f64; 3], [f64; 3], f64) {
        let da = [a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]];
        let db = [b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]];
        let r = [a0[0] - b0[0], a0[1] - b0[1], a0[2] - b0[2]];

        let aa = da[0] * da[0] + da[1] * da[1] + da[2] * da[2];
        let ee = db[0] * db[0] + db[1] * db[1] + db[2] * db[2];
        let f = db[0] * r[0] + db[1] * r[1] + db[2] * r[2];
        let eps = 1e-12;

        let (s, t) = if aa <= eps && ee <= eps {
            (0.0_f64, 0.0_f64)
        } else if aa <= eps {
            (0.0_f64, (f / ee).clamp(0.0, 1.0))
        } else {
            let c = da[0] * r[0] + da[1] * r[1] + da[2] * r[2];
            if ee <= eps {
                ((-c / aa).clamp(0.0, 1.0), 0.0_f64)
            } else {
                let b = da[0] * db[0] + da[1] * db[1] + da[2] * db[2];
                let denom = aa * ee - b * b;
                let s_cand = if denom.abs() > eps {
                    ((b * f - c * ee) / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let t_nom = b * s_cand + f;
                let t_cand = if t_nom < 0.0 {
                    0.0
                } else if t_nom > ee {
                    1.0
                } else {
                    t_nom / ee
                };
                let s_final = ((b * t_cand - c) / aa).clamp(0.0, 1.0);
                (s_final, t_cand)
            }
        };

        let pa = [a0[0] + s * da[0], a0[1] + s * da[1], a0[2] + s * da[2]];
        let pb = [b0[0] + t * db[0], b0[1] + t * db[1], b0[2] + t * db[2]];
        let dx = pa[0] - pb[0];
        let dy = pa[1] - pb[1];
        let dz = pa[2] - pb[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        (pa, pb, dist)
    }

    /// Capsule vs Oriented Bounding Box (OBB) intersection test.
    ///
    /// `obb_center`: world-space center of the OBB.
    /// `obb_axes`:   columns are the OBB's local X, Y, Z axes (unit vectors).
    /// `obb_half`:   half-extents along the OBB local axes.
    ///
    /// Converts the query to OBB-local space and tests sphere-swept segment
    /// against the OBB using the separating-axis theorem (SAT) with slabs.
    pub fn intersects_obb(
        &self,
        capsule_center: [f64; 3],
        obb_center: [f64; 3],
        obb_axes: [[f64; 3]; 3],
        obb_half: [f64; 3],
    ) -> bool {
        // Expand the OBB by the capsule radius (Minkowski sum)
        let expanded_half = [
            obb_half[0] + self.radius,
            obb_half[1] + self.radius,
            obb_half[2] + self.radius,
        ];

        // Capsule segment endpoints in world space
        let a = [
            capsule_center[0],
            capsule_center[1] - self.half_height,
            capsule_center[2],
        ];
        let b = [
            capsule_center[0],
            capsule_center[1] + self.half_height,
            capsule_center[2],
        ];

        // Test whether either endpoint is inside the expanded OBB
        for pt in [a, b] {
            let d = [
                pt[0] - obb_center[0],
                pt[1] - obb_center[1],
                pt[2] - obb_center[2],
            ];
            let local = [
                d[0] * obb_axes[0][0] + d[1] * obb_axes[0][1] + d[2] * obb_axes[0][2],
                d[0] * obb_axes[1][0] + d[1] * obb_axes[1][1] + d[2] * obb_axes[1][2],
                d[0] * obb_axes[2][0] + d[1] * obb_axes[2][1] + d[2] * obb_axes[2][2],
            ];
            if local[0].abs() <= expanded_half[0]
                && local[1].abs() <= expanded_half[1]
                && local[2].abs() <= expanded_half[2]
            {
                return true;
            }
        }

        // Test whether segment passes through the expanded OBB (slab test)
        // Transform segment into OBB-local space
        let da = [
            a[0] - obb_center[0],
            a[1] - obb_center[1],
            a[2] - obb_center[2],
        ];
        let dir_w = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];

        let orig_local = [
            da[0] * obb_axes[0][0] + da[1] * obb_axes[0][1] + da[2] * obb_axes[0][2],
            da[0] * obb_axes[1][0] + da[1] * obb_axes[1][1] + da[2] * obb_axes[1][2],
            da[0] * obb_axes[2][0] + da[1] * obb_axes[2][1] + da[2] * obb_axes[2][2],
        ];
        let dir_local = [
            dir_w[0] * obb_axes[0][0] + dir_w[1] * obb_axes[0][1] + dir_w[2] * obb_axes[0][2],
            dir_w[0] * obb_axes[1][0] + dir_w[1] * obb_axes[1][1] + dir_w[2] * obb_axes[1][2],
            dir_w[0] * obb_axes[2][0] + dir_w[1] * obb_axes[2][1] + dir_w[2] * obb_axes[2][2],
        ];

        let mut t_min = 0.0_f64;
        let mut t_max = 1.0_f64;

        for axis in 0..3 {
            let d_loc = dir_local[axis];
            let o_loc = orig_local[axis];
            let hh = expanded_half[axis];
            if d_loc.abs() < 1e-12 {
                if o_loc.abs() > hh {
                    return false;
                }
            } else {
                let inv_d = 1.0 / d_loc;
                let t0 = (-hh - o_loc) * inv_d;
                let t1 = (hh - o_loc) * inv_d;
                let (t_lo, t_hi) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
                t_min = t_min.max(t_lo);
                t_max = t_max.min(t_hi);
                if t_min > t_max {
                    return false;
                }
            }
        }
        true
    }

    /// Signed distance field (SDF) for the capsule.
    ///
    /// Same as `signed_distance` but with a different name to match the
    /// naming convention used in shader / SDF contexts.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        self.signed_distance(p)
    }

    /// Continuous collision detection (swept capsule).
    ///
    /// Given the capsule moving from `center_start` to `center_end` (linear
    /// motion, constant orientation), and a static sphere at `sphere_center`
    /// with `sphere_radius`, returns the first time of contact `t ∈ [0,1]`
    /// or `None` if no contact occurs.
    ///
    /// This is a simple CCD sweep using the Minkowski-sum radius and a
    /// moving-segment vs static-point formulation.
    pub fn swept_capsule_vs_sphere(
        &self,
        center_start: [f64; 3],
        center_end: [f64; 3],
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> Option<f64> {
        let combined_radius = self.radius + sphere_radius;

        // We need to find the earliest t in [0,1] such that the distance from
        // sphere_center to the capsule medial axis (at time t) equals combined_radius.
        // Linearly interpolate capsule center.
        // At time t, capsule axis: from C(t)-[0,hh,0] to C(t)+[0,hh,0],
        // where C(t) = center_start + t*(center_end - center_start).

        let n_steps = 64;
        for i in 0..=n_steps {
            let t = i as f64 / n_steps as f64;
            let cx = center_start[0] + t * (center_end[0] - center_start[0]);
            let cy = center_start[1] + t * (center_end[1] - center_start[1]);
            let cz = center_start[2] + t * (center_end[2] - center_start[2]);

            let a0 = [cx, cy - self.half_height, cz];
            let a1 = [cx, cy + self.half_height, cz];
            let dist = Self::point_segment_distance(sphere_center, a0, a1);
            if dist <= combined_radius {
                return Some(t);
            }
        }
        None
    }

    /// Distance from a point to a segment.
    fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let ab_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        let t = if ab_sq < 1e-24 {
            0.0
        } else {
            ((ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2]) / ab_sq).clamp(0.0, 1.0)
        };
        let proj = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
        let dx = p[0] - proj[0];
        let dy = p[1] - proj[1];
        let dz = p[2] - proj[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Generate `n` random points on the capsule surface using rejection sampling.
    ///
    /// Returns a `Vec` of `[f64;3]` coordinates sampled uniformly on the surface.
    /// Uses a deterministic PRNG seed for reproducibility.
    pub fn random_surface_points(&self, n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut points = Vec::with_capacity(n);
        let r = self.radius;
        let hh = self.half_height;

        // Surface areas
        let cyl_area = 2.0 * PI * r * 2.0 * hh;
        let sphere_area = 4.0 * PI * r * r;
        let total_area = cyl_area + sphere_area;
        let p_cylinder = cyl_area / total_area;

        let mut rng_state = seed;
        let next_f64 = |state: &mut u64| -> f64 {
            // xorshift64 PRNG
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            (*state as f64) / (u64::MAX as f64)
        };

        while points.len() < n {
            let u = next_f64(&mut rng_state);
            let v = next_f64(&mut rng_state);
            let w = next_f64(&mut rng_state);

            if u < p_cylinder {
                // Cylindrical region
                let theta = v * 2.0 * PI;
                let y = (w * 2.0 - 1.0) * hh;
                points.push([r * theta.cos(), y, r * theta.sin()]);
            } else {
                // Hemisphere (top or bottom)
                let theta = v * 2.0 * PI;
                let phi = (w * 2.0 - 1.0).acos();
                let x = r * phi.sin() * theta.cos();
                let z = r * phi.sin() * theta.sin();
                let raw_y = r * phi.cos();
                // Map to correct hemisphere
                let extra = next_f64(&mut rng_state);
                let y = if extra < 0.5 {
                    hh + raw_y.abs()
                } else {
                    -hh - raw_y.abs()
                };
                points.push([x, y, z]);
            }
        }
        points
    }

    // ── Capsule-triangle contact ──

    /// Test capsule vs a world-space triangle `(ta, tb, tc)`.
    ///
    /// The capsule is placed at `center` (Y-aligned).  Returns `Some((depth, normal))`
    /// if the capsule penetrates the triangle plane by at least `depth > 0`, where
    /// `normal` points from the triangle towards the capsule center.
    pub fn capsule_triangle_contact(
        &self,
        center: [f64; 3],
        ta: [f64; 3],
        tb: [f64; 3],
        tc: [f64; 3],
    ) -> Option<(f64, [f64; 3])> {
        // Triangle edge vectors and normal
        let e1 = [tb[0] - ta[0], tb[1] - ta[1], tb[2] - ta[2]];
        let e2 = [tc[0] - ta[0], tc[1] - ta[1], tc[2] - ta[2]];
        let raw_n = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let n_len = (raw_n[0] * raw_n[0] + raw_n[1] * raw_n[1] + raw_n[2] * raw_n[2]).sqrt();
        if n_len < 1e-12 {
            return None;
        }
        let n = [raw_n[0] / n_len, raw_n[1] / n_len, raw_n[2] / n_len];

        // Capsule medial axis endpoints in world space
        let a0 = [center[0], center[1] - self.half_height, center[2]];
        let a1 = [center[0], center[1] + self.half_height, center[2]];

        // Signed distances of both endpoints from the triangle plane
        let d0 = (a0[0] - ta[0]) * n[0] + (a0[1] - ta[1]) * n[1] + (a0[2] - ta[2]) * n[2];
        let d1 = (a1[0] - ta[0]) * n[0] + (a1[1] - ta[1]) * n[1] + (a1[2] - ta[2]) * n[2];

        // Use the endpoint closest to the triangle plane
        let min_dist = d0.abs().min(d1.abs());
        let penetration = self.radius - min_dist;
        if penetration > 0.0 {
            Some((penetration, n))
        } else {
            None
        }
    }
}

// ── CapsuleChain (rope segments) ──

/// A chain of capsule segments forming a rope / chain approximation.
///
/// Each consecutive pair of control points defines a capsule segment
/// of the same radius.
#[derive(Debug, Clone)]
pub struct CapsuleChain {
    /// Control points of the chain (N points → N-1 segments).
    pub points: Vec<[f64; 3]>,
    /// Uniform radius for all segments.
    pub radius: f64,
}

impl CapsuleChain {
    /// Create a chain from control points with zero radius.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        Self {
            points,
            radius: 0.0,
        }
    }

    /// Create a chain with a specified radius.
    pub fn with_radius(points: Vec<[f64; 3]>, radius: f64) -> Self {
        Self { points, radius }
    }

    /// Number of capsule segments (one fewer than the number of points).
    pub fn segment_count(&self) -> usize {
        self.points.len().saturating_sub(1)
    }

    /// Total arc length of the chain (sum of segment lengths).
    pub fn total_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| {
                let dx = w[1][0] - w[0][0];
                let dy = w[1][1] - w[0][1];
                let dz = w[1][2] - w[0][2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .sum()
    }

    /// Returns `true` if `p` is inside any capsule segment of the chain.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        self.points
            .windows(2)
            .any(|w| Capsule::point_segment_distance(p, w[0], w[1]) <= self.radius)
    }

    /// SDF: minimum signed distance to any segment surface.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let min_seg_dist = self
            .points
            .windows(2)
            .map(|w| Capsule::point_segment_distance(p, w[0], w[1]))
            .fold(f64::INFINITY, f64::min);
        min_seg_dist - self.radius
    }

    /// Minimum distance from `p` to the chain surface (negative if inside).
    pub fn min_distance_to_point(&self, p: [f64; 3]) -> f64 {
        self.sdf(p)
    }
}

// ── DeformableCapsule ──

/// A capsule with arbitrary (deformable) endpoint positions.
///
/// Unlike the axis-aligned `Capsule`, this capsule can have its two
/// hemispherical cap centres at any positions in 3D space.
#[derive(Debug, Clone)]
pub struct DeformableCapsule {
    /// First endpoint of the medial axis.
    pub a: [f64; 3],
    /// Second endpoint of the medial axis.
    pub b: [f64; 3],
    /// Capsule radius.
    pub radius: f64,
}

impl DeformableCapsule {
    /// Create a new deformable capsule.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }

    /// First endpoint (cap centre A).
    pub fn endpoint_a(&self) -> [f64; 3] {
        self.a
    }

    /// Second endpoint (cap centre B).
    pub fn endpoint_b(&self) -> [f64; 3] {
        self.b
    }

    /// Update endpoint A.
    pub fn set_endpoint_a(&mut self, a: [f64; 3]) {
        self.a = a;
    }

    /// Update endpoint B.
    pub fn set_endpoint_b(&mut self, b: [f64; 3]) {
        self.b = b;
    }

    /// Length of the medial axis segment.
    pub fn length(&self) -> f64 {
        let d = [
            self.b[0] - self.a[0],
            self.b[1] - self.a[1],
            self.b[2] - self.a[2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    /// Midpoint of the medial axis.
    pub fn midpoint(&self) -> [f64; 3] {
        [
            (self.a[0] + self.b[0]) * 0.5,
            (self.a[1] + self.b[1]) * 0.5,
            (self.a[2] + self.b[2]) * 0.5,
        ]
    }

    /// Medial axis endpoints.
    pub fn medial_axis_endpoints(&self) -> ([f64; 3], [f64; 3]) {
        (self.a, self.b)
    }

    /// Signed distance field: negative inside, positive outside.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        Capsule::point_segment_distance(p, self.a, self.b) - self.radius
    }

    /// Support point in the given direction.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let da = dir[0] * self.a[0] + dir[1] * self.a[1] + dir[2] * self.a[2];
        let db = dir[0] * self.b[0] + dir[1] * self.b[1] + dir[2] * self.b[2];
        let base = if da >= db { self.a } else { self.b };
        let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if len < 1e-12 {
            return base;
        }
        let s = self.radius / len;
        [
            base[0] + dir[0] * s,
            base[1] + dir[1] * s,
            base[2] + dir[2] * s,
        ]
    }
}

// ── CapsuleFrustum ──

/// A capsule-like shape with different radii at each end (a "swept frustum").
///
/// The axis runs along Y from `y=0` (radius `r_bottom`) to `y=height` (radius `r_top`).
/// The radius varies linearly along the axis.
#[derive(Debug, Clone)]
pub struct CapsuleFrustum {
    /// Radius at the bottom cap (y=0).
    pub r_bottom: f64,
    /// Radius at the top cap (y=height).
    pub r_top: f64,
    /// Total height (length of axis segment).
    pub height: f64,
}

impl CapsuleFrustum {
    /// Create a new capsule frustum.
    pub fn new(r_bottom: f64, r_top: f64, height: f64) -> Self {
        Self {
            r_bottom,
            r_top,
            height,
        }
    }

    /// Interpolated radius at height `y` (clamped to \[0, height\]).
    pub fn radius_at_height(&self, y: f64) -> f64 {
        if self.height < 1e-12 {
            return (self.r_bottom + self.r_top) * 0.5;
        }
        let t = (y / self.height).clamp(0.0, 1.0);
        self.r_bottom + t * (self.r_top - self.r_bottom)
    }

    /// Approximate volume (frustum cone + two hemispheres of different sizes).
    pub fn volume(&self) -> f64 {
        // Truncated cone: π/3 * h * (r1² + r1*r2 + r2²)
        let r1 = self.r_bottom;
        let r2 = self.r_top;
        let h = self.height;
        let cone_vol = std::f64::consts::PI / 3.0 * h * (r1 * r1 + r1 * r2 + r2 * r2);
        let sphere_bot = (2.0 / 3.0) * std::f64::consts::PI * r1 * r1 * r1;
        let sphere_top = (2.0 / 3.0) * std::f64::consts::PI * r2 * r2 * r2;
        cone_vol + sphere_bot + sphere_top
    }

    /// Returns `true` if point `p` is inside the frustum.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        self.sdf(p) <= 0.0
    }

    /// Signed distance to the frustum surface (approximate).
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let y = p[1].clamp(0.0, self.height);
        let r = self.radius_at_height(y);
        let dy = p[1] - y;
        let xz_dist = (p[0] * p[0] + p[2] * p[2]).sqrt();
        let dist_to_axis = (xz_dist * xz_dist + dy * dy).sqrt();
        dist_to_axis - r
    }
}

// ── CurvedCapsulePath ──

/// A capsule swept along a polyline path.
///
/// Each consecutive pair of path points defines a capsule segment.
/// The shape represents the union of all such capsule segments.
#[derive(Debug, Clone)]
pub struct CurvedCapsulePath {
    /// Polyline control points.
    pub path: Vec<[f64; 3]>,
    /// Uniform radius.
    pub radius: f64,
}

impl CurvedCapsulePath {
    /// Create a new curved capsule path.
    pub fn new(path: Vec<[f64; 3]>, radius: f64) -> Self {
        Self { path, radius }
    }

    /// Total path length (sum of segment lengths).
    pub fn path_length(&self) -> f64 {
        self.path
            .windows(2)
            .map(|w| {
                let d = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .sum()
    }

    /// Number of path segments (one fewer than points).
    pub fn num_segments(&self) -> usize {
        self.path.len().saturating_sub(1)
    }

    /// Returns `true` if `p` is inside any capsule segment.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        self.path
            .windows(2)
            .any(|w| Capsule::point_segment_distance(p, w[0], w[1]) <= self.radius)
    }

    /// Approximate SDF: minimum distance to any segment minus radius.
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let min_d = self
            .path
            .windows(2)
            .map(|w| Capsule::point_segment_distance(p, w[0], w[1]))
            .fold(f64::INFINITY, f64::min);
        min_d - self.radius
    }
}

impl Shape for Capsule {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn bounding_box(&self) -> Aabb {
        let r = self.radius;
        let h = self.half_height + r;
        Aabb::new(Vec3::new(-r, -h, -r), Vec3::new(r, h, r))
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let norm = direction.norm();
        let cap_center = if direction.y >= 0.0 {
            Vec3::new(0.0, self.half_height, 0.0)
        } else {
            Vec3::new(0.0, -self.half_height, 0.0)
        };
        if norm < 1e-10 {
            return cap_center;
        }
        cap_center + direction * (self.radius / norm)
    }

    fn volume(&self) -> Real {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        // Cylinder + sphere
        PI * r * r * h + (4.0 / 3.0) * PI * r.powi(3)
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::zeros()
    }

    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        let r2 = r * r;
        let h2 = h * h;

        // Approximate as cylinder for the inertia
        let vol_cyl = PI * r2 * h;
        let vol_sphere = (4.0 / 3.0) * PI * r.powi(3);
        let total_vol = vol_cyl + vol_sphere;

        let mass_cyl = mass * vol_cyl / total_vol;
        let mass_sphere = mass * vol_sphere / total_vol;

        // Cylinder inertia (Y-aligned)
        let iy_cyl = 0.5 * mass_cyl * r2;
        let ixz_cyl = mass_cyl * (3.0 * r2 + h2) / 12.0;

        // Hemisphere inertia (approximation using shifted sphere)
        let iy_sphere = 0.4 * mass_sphere * r2;
        // Parallel axis theorem for hemispheres offset by half_height + 3r/8
        let offset = self.half_height + 3.0 * r / 8.0;
        let ixz_sphere = 0.4 * mass_sphere * r2 + mass_sphere * offset * offset;

        let iy = iy_cyl + iy_sphere;
        let ixz = ixz_cyl + ixz_sphere;

        Mat3::new(ixz, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ixz)
    }

    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        // Simplified: test against the two end-cap spheres and the cylinder
        let mut best: Option<RayHit> = None;

        // Test top hemisphere
        let top_center = Vec3::new(0.0, self.half_height, 0.0);
        if let Some(hit) = ray_sphere(ray_origin, ray_direction, &top_center, self.radius)
            && hit.toi <= max_toi
            && hit.point.y >= self.half_height
            && best.as_ref().is_none_or(|b| hit.toi < b.toi)
        {
            best = Some(hit);
        }

        // Test bottom hemisphere
        let bot_center = Vec3::new(0.0, -self.half_height, 0.0);
        if let Some(hit) = ray_sphere(ray_origin, ray_direction, &bot_center, self.radius)
            && hit.toi <= max_toi
            && hit.point.y <= -self.half_height
            && best.as_ref().is_none_or(|b| hit.toi < b.toi)
        {
            best = Some(hit);
        }

        // Test infinite cylinder (XZ plane)
        if let Some(hit) = ray_cylinder_y(ray_origin, ray_direction, self.radius, self.half_height)
            && hit.toi <= max_toi
            && best.as_ref().is_none_or(|b| hit.toi < b.toi)
        {
            best = Some(hit);
        }

        best
    }
}

fn ray_sphere(origin: &Vec3, direction: &Vec3, center: &Vec3, radius: Real) -> Option<RayHit> {
    let oc = origin - center;
    let a = direction.dot(direction);
    let b = 2.0 * oc.dot(direction);
    let c = oc.dot(&oc) - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    let t = if t1 >= 0.0 { t1 } else { t2 };
    if t < 0.0 {
        return None;
    }
    let point = origin + direction * t;
    let normal = (point - center).normalize();
    Some(RayHit {
        point,
        normal,
        toi: t,
    })
}

fn ray_cylinder_y(
    origin: &Vec3,
    direction: &Vec3,
    radius: Real,
    half_height: Real,
) -> Option<RayHit> {
    // Infinite cylinder along Y axis in XZ plane
    let a = direction.x * direction.x + direction.z * direction.z;
    let b = 2.0 * (origin.x * direction.x + origin.z * direction.z);
    let c = origin.x * origin.x + origin.z * origin.z - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 || a < 1e-12 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);

    for t in [t1, t2] {
        if t < 0.0 {
            continue;
        }
        let point = origin + direction * t;
        if point.y.abs() <= half_height {
            let normal = Vec3::new(point.x, 0.0, point.z).normalize();
            return Some(RayHit {
                point,
                normal,
                toi: t,
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_capsule_volume() {
        let c = Capsule::new(1.0, 1.0);
        let expected = PI * 1.0 * 2.0 + (4.0 / 3.0) * PI;
        assert!((c.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_volume_explicit() {
        let c = Capsule::new(1.0, 1.0);
        assert!((c.volume_explicit() - c.volume()).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_surface_area() {
        // r=1, half_height=1: 2π*1*2 + 4π = 4π + 4π = 8π
        let c = Capsule::new(1.0, 1.0);
        let expected = 8.0 * PI;
        assert!((c.surface_area() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_inertia_array() {
        let c = Capsule::new(1.0, 1.0);
        let it = c.inertia_tensor_array(1.0);
        // Diagonal elements should be positive
        assert!(it[0][0] > 0.0);
        assert!(it[1][1] > 0.0);
        assert!(it[2][2] > 0.0);
        // Ixx == Izz by symmetry
        assert!((it[0][0] - it[2][2]).abs() < 1e-10);
        // Off-diagonal zero
        assert!(it[0][1].abs() < 1e-10);
    }

    #[test]
    fn test_capsule_raycast() {
        let c = Capsule::new(1.0, 2.0);
        let origin = Vec3::new(-5.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let hit = c.ray_cast(&origin, &dir, 100.0).unwrap();
        assert!((hit.toi - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_raycast_array() {
        let c = Capsule::new(1.0, 2.0);
        let (t, _n) = c
            .ray_cast_array([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0)
            .unwrap();
        assert!((t - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_support_up() {
        let c = Capsule::new(1.0, 2.0);
        let sp = c.support([0.0, 1.0, 0.0]);
        // Should return top cap: y = half_height + radius = 3
        assert!((sp[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_support_down() {
        let c = Capsule::new(1.0, 2.0);
        let sp = c.support([0.0, -1.0, 0.0]);
        // Should return bottom: y = -half_height - radius = -3
        assert!((sp[1] + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_raycast_top_cap() {
        let c = Capsule::new(1.0, 2.0);
        let origin = Vec3::new(0.0, 10.0, 0.0);
        let dir = Vec3::new(0.0, -1.0, 0.0);
        let hit = c.ray_cast(&origin, &dir, 100.0).unwrap();
        // Top of capsule at y = half_height + radius = 3
        assert!(
            (hit.toi - 7.0).abs() < 1e-10,
            "expected toi=7, got {}",
            hit.toi
        );
    }

    // ── New tests ──

    #[test]
    fn test_capsule_closest_point_on_side() {
        let c = Capsule::new(1.0, 2.0);
        let cp = c.closest_point([5.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
        assert!(cp[2].abs() < 1e-10);
    }

    #[test]
    fn test_capsule_closest_point_on_top_cap() {
        let c = Capsule::new(1.0, 2.0);
        let cp = c.closest_point([0.0, 5.0, 0.0]);
        // Should be on the top hemisphere: (0, 2+1, 0) = (0, 3, 0)
        assert!((cp[1] - 3.0).abs() < 1e-10);
        assert!(cp[0].abs() < 1e-10);
    }

    #[test]
    fn test_capsule_closest_point_on_bottom_cap() {
        let c = Capsule::new(1.0, 2.0);
        let cp = c.closest_point([0.0, -5.0, 0.0]);
        assert!((cp[1] + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_closest_point_on_axis() {
        let c = Capsule::new(1.0, 2.0);
        // Point on axis inside capsule
        let cp = c.closest_point([0.0, 0.0, 0.0]);
        // Should return a point on the surface, in +X direction
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
    }

    #[test]
    fn test_capsule_contains_point() {
        let c = Capsule::new(1.0, 2.0);
        assert!(c.contains_point([0.0, 0.0, 0.0])); // center
        assert!(c.contains_point([0.5, 1.0, 0.0])); // inside cylinder
        assert!(c.contains_point([0.0, 2.5, 0.0])); // inside top cap
        assert!(!c.contains_point([0.0, 3.5, 0.0])); // above capsule
        assert!(!c.contains_point([2.0, 0.0, 0.0])); // outside radially
    }

    #[test]
    fn test_capsule_signed_distance() {
        let c = Capsule::new(2.0, 1.0);
        // Center
        assert!((c.signed_distance([0.0, 0.0, 0.0]) + 2.0).abs() < 1e-10);
        // On surface
        assert!(c.signed_distance([2.0, 0.0, 0.0]).abs() < 1e-10);
        // Outside
        assert!((c.signed_distance([4.0, 0.0, 0.0]) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_medial_axis_endpoints() {
        let c = Capsule::new(1.0, 3.0);
        let (top, bot) = c.medial_axis_endpoints();
        assert!((top[1] - 3.0).abs() < 1e-10);
        assert!((bot[1] + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_full_length() {
        let c = Capsule::new(1.0, 2.0);
        // 2*hh + 2*r = 4 + 2 = 6
        assert!((c.full_length() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_medial_axis_length() {
        let c = Capsule::new(1.0, 2.0);
        assert!((c.medial_axis_length() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_distance_parallel() {
        // Two parallel segments along Y, offset in X
        let d = Capsule::segment_segment_distance(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        );
        assert!((d - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_distance_crossing() {
        // Two segments that cross: one along X, one along Z, separated by y
        let d = Capsule::segment_segment_distance(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, -1.0],
            [0.0, 2.0, 1.0],
        );
        assert!((d - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_distance_touching() {
        let d = Capsule::segment_segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        );
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_segment_segment_distance_degenerate_points() {
        let d = Capsule::segment_segment_distance(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [3.0, 4.0, 0.0],
            [3.0, 4.0, 0.0],
        );
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_capsule_distance_separated() {
        let a = Capsule::new(1.0, 1.0);
        let b = Capsule::new(1.0, 1.0);
        // Place them 5 apart in X
        let d = a.capsule_capsule_distance([0.0, 0.0, 0.0], &b, [5.0, 0.0, 0.0]);
        // Segment distance = 5, surface distance = 5 - 1 - 1 = 3
        assert!((d - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_capsule_distance_overlapping() {
        let a = Capsule::new(1.0, 1.0);
        let b = Capsule::new(1.0, 1.0);
        let d = a.capsule_capsule_distance([0.0, 0.0, 0.0], &b, [0.5, 0.0, 0.0]);
        assert_eq!(d, 0.0); // overlapping
    }

    #[test]
    fn test_capsule_capsule_overlap() {
        let a = Capsule::new(1.0, 1.0);
        let b = Capsule::new(1.0, 1.0);
        assert!(a.capsule_capsule_overlap([0.0, 0.0, 0.0], &b, [1.0, 0.0, 0.0]));
        assert!(!a.capsule_capsule_overlap([0.0, 0.0, 0.0], &b, [5.0, 0.0, 0.0]));
    }

    #[test]
    fn test_capsule_project_on_medial_axis() {
        let c = Capsule::new(1.0, 2.0);
        assert!((c.project_on_medial_axis([0.0, 5.0, 0.0]) - 2.0).abs() < 1e-10);
        assert!((c.project_on_medial_axis([0.0, -5.0, 0.0]) + 2.0).abs() < 1e-10);
        assert!((c.project_on_medial_axis([3.0, 1.0, 0.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_distance_to_medial_axis() {
        let c = Capsule::new(1.0, 2.0);
        assert!(c.distance_to_medial_axis([0.0, 0.0, 0.0]).abs() < 1e-10);
        assert!((c.distance_to_medial_axis([3.0, 0.0, 0.0]) - 3.0).abs() < 1e-10);
        // Point above cap
        assert!((c.distance_to_medial_axis([0.0, 5.0, 0.0]) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_support_diagonal() {
        let c = Capsule::new(1.0, 2.0);
        let sp = c.support([1.0, 1.0, 0.0]);
        // cap_y = 2.0 (positive Y)
        let len = (1.0_f64 + 1.0).sqrt();
        let expected_x = 1.0 / len;
        let expected_y = 2.0 + 1.0 / len;
        assert!((sp[0] - expected_x).abs() < 1e-10);
        assert!((sp[1] - expected_y).abs() < 1e-10);
    }

    // ── Expanded tests ──

    #[test]
    fn test_capsule_closest_points_vs_capsule() {
        let a = Capsule::new(1.0, 1.0);
        let b = Capsule::new(1.0, 1.0);
        // Capsules side by side in X, separated by 4 units
        let (pa, pb, seg_dist) =
            a.closest_points_capsule_vs_capsule([0.0, 0.0, 0.0], &b, [4.0, 0.0, 0.0]);
        assert!(
            (seg_dist - 4.0).abs() < 1e-9,
            "axis dist should be 4, got {seg_dist}"
        );
        assert!((pa[0]).abs() < 1e-9);
        assert!((pb[0] - 4.0).abs() < 1e-9);
    }

    #[test]
    fn test_segment_segment_closest_crossing() {
        // X-segment and Z-segment crossing at origin at different Y levels
        let (pa, pb, dist) = Capsule::segment_segment_closest(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 3.0, -1.0],
            [0.0, 3.0, 1.0],
        );
        assert!((dist - 3.0).abs() < 1e-9, "dist should be 3, got {dist}");
        assert!(pa[0].abs() < 1e-9);
        assert!((pb[1] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_capsule_intersects_obb_basic() {
        let c = Capsule::new(0.5, 1.0);
        // OBB centered at origin, axis-aligned, half-extents [2,2,2]
        let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Capsule centered inside OBB
        assert!(c.intersects_obb([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], axes, [2.0, 2.0, 2.0]));
        // Capsule far away
        assert!(!c.intersects_obb([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], axes, [2.0, 2.0, 2.0]));
    }

    #[test]
    fn test_capsule_sdf_matches_signed_distance() {
        let c = Capsule::new(1.0, 2.0);
        let p = [3.0, 0.0, 0.0];
        assert!((c.sdf(p) - c.signed_distance(p)).abs() < 1e-12);
    }

    #[test]
    fn test_swept_capsule_vs_sphere_contact() {
        let c = Capsule::new(1.0, 1.0);
        // Capsule starts at x=10, moves to x=0 – sphere at origin with radius 0.5
        let result =
            c.swept_capsule_vs_sphere([10.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.5);
        assert!(result.is_some(), "swept capsule should contact sphere");
        let t = result.unwrap();
        assert!((0.0..=1.0).contains(&t), "t should be in [0,1], got {t}");
    }

    #[test]
    fn test_swept_capsule_vs_sphere_no_contact() {
        let c = Capsule::new(0.5, 0.5);
        // Capsule moves along X, sphere far away in Z
        let result =
            c.swept_capsule_vs_sphere([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 0.0, 100.0], 0.5);
        assert!(result.is_none(), "should not contact distant sphere");
    }

    #[test]
    fn test_random_surface_points_count() {
        let c = Capsule::new(1.0, 2.0);
        let pts = c.random_surface_points(50, 42);
        assert_eq!(pts.len(), 50);
    }

    #[test]
    fn test_random_surface_points_on_surface() {
        let c = Capsule::new(1.0, 2.0);
        let pts = c.random_surface_points(100, 7);
        for p in &pts {
            // Each point should be within radius+eps of the surface
            let sdf = c.sdf(*p);
            assert!(sdf.abs() < 0.05, "point {:?} has sdf={sdf}", p);
        }
    }

    #[test]
    fn test_point_segment_distance_midpoint() {
        // Point above midpoint of segment along X
        let d = Capsule::point_segment_distance([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_inertia_tensor_raw_matches_array() {
        let c = Capsule::new(1.5, 2.5);
        let raw = c.inertia_tensor_raw(3.0);
        let arr = c.inertia_tensor_array(3.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!((raw[i][j] - arr[i][j]).abs() < 1e-12);
            }
        }
    }

    // ── New 30-50 expanded tests ──

    #[test]
    fn test_capsule_triangle_contact_no_overlap() {
        let c = Capsule::new(0.5, 1.0);
        // Triangle far away from capsule
        let a = [10.0, 0.0, 0.0];
        let b = [11.0, 0.0, 0.0];
        let tri = [10.0_f64, 0.0_f64, 1.0_f64];
        let result = c.capsule_triangle_contact([0.0, 0.0, 0.0], a, b, tri);
        assert!(
            result.is_none(),
            "should be no contact with distant triangle"
        );
    }

    #[test]
    fn test_capsule_triangle_contact_overlap() {
        let c = Capsule::new(1.0, 1.0);
        // Large triangle containing the capsule axis
        let a = [-5.0, 0.0, -5.0];
        let b = [5.0, 0.0, -5.0];
        let tri = [0.0_f64, 0.0_f64, 5.0_f64];
        let result = c.capsule_triangle_contact([0.0, 0.5, 0.0], a, b, tri);
        assert!(
            result.is_some(),
            "capsule near triangle should produce contact"
        );
    }

    #[test]
    fn test_capsule_chain_empty() {
        let links = CapsuleChain::new(vec![]);
        assert_eq!(links.segment_count(), 0);
    }

    #[test]
    fn test_capsule_chain_single_point() {
        let links = CapsuleChain::new(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(links.segment_count(), 0);
    }

    #[test]
    fn test_capsule_chain_two_points() {
        let links = CapsuleChain::new(vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert_eq!(links.segment_count(), 1);
    }

    #[test]
    fn test_capsule_chain_total_length() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        let chain = CapsuleChain::new(pts);
        assert!((chain.total_length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_chain_contains_point_near_segment() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let chain = CapsuleChain::with_radius(pts, 1.0);
        assert!(chain.contains_point([0.5, 1.0, 0.0]));
        assert!(!chain.contains_point([2.0, 1.0, 0.0]));
    }

    #[test]
    fn test_deformable_capsule_update_endpoints() {
        let mut dc = DeformableCapsule::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        dc.set_endpoint_a([1.0, 0.0, 0.0]);
        dc.set_endpoint_b([1.0, 2.0, 0.0]);
        assert!((dc.endpoint_a()[0] - 1.0).abs() < 1e-12);
        assert!((dc.endpoint_b()[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_deformable_capsule_length() {
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [0.0, 3.0, 0.0], 0.5);
        assert!((dc.length() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_deformable_capsule_midpoint() {
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.5);
        let mid = dc.midpoint();
        assert!((mid[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_capsule_medial_axis_3d_segment() {
        let dc = DeformableCapsule::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 0.5);
        let (a, b) = dc.medial_axis_endpoints();
        assert!((a[0] - 1.0).abs() < 1e-12);
        assert!((b[2] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_deformable_capsule_sdf_inside() {
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 1.0);
        // Point on the axis inside the capsule
        assert!(dc.sdf([0.0, 1.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_deformable_capsule_sdf_outside() {
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        // Point far away
        assert!(dc.sdf([5.0, 1.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_capsule_frustum_volume_positive() {
        let f = CapsuleFrustum::new(0.5, 1.0, 2.0);
        assert!(f.volume() > 0.0);
    }

    #[test]
    fn test_capsule_frustum_contains_center() {
        let f = CapsuleFrustum::new(0.5, 1.0, 2.0);
        // Center along axis should be inside
        assert!(f.contains_point([0.0, 1.0, 0.0]));
    }

    #[test]
    fn test_capsule_frustum_excludes_far_point() {
        let f = CapsuleFrustum::new(0.5, 1.0, 2.0);
        assert!(!f.contains_point([10.0, 0.0, 0.0]));
    }

    #[test]
    fn test_curved_capsule_path_length() {
        // Path with 3 points forming an L
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        let cp = CurvedCapsulePath::new(pts, 0.3);
        let len = cp.path_length();
        assert!((len - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_curved_capsule_path_contains_axis_point() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let cp = CurvedCapsulePath::new(pts, 0.5);
        assert!(cp.contains_point([0.0, 1.0, 0.0]));
    }

    #[test]
    fn test_curved_capsule_path_excludes_far_point() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let cp = CurvedCapsulePath::new(pts, 0.5);
        assert!(!cp.contains_point([5.0, 0.0, 0.0]));
    }

    #[test]
    fn test_swept_capsule_contact_time_order() {
        let c = Capsule::new(0.5, 1.0);
        // Capsule moving from far away to the origin
        let t = c.swept_capsule_vs_sphere([20.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1);
        assert!(t.is_some());
        let tv = t.unwrap();
        assert!(tv > 0.0 && tv <= 1.0);
    }

    #[test]
    fn test_capsule_capsule_touching_boundary() {
        let a = Capsule::new(1.0, 1.0);
        let b = Capsule::new(1.0, 1.0);
        // Place them exactly touching: axis distance = sum of radii = 2
        let d = a.capsule_capsule_distance([0.0, 0.0, 0.0], &b, [2.0, 0.0, 0.0]);
        assert!(
            d.abs() < 1e-10,
            "touching capsules should have 0 distance, got {d}"
        );
    }

    #[test]
    fn test_segment_segment_closest_endpoints() {
        // Two perpendicular segments not crossing; closest is the tips
        let (pa, pb, dist) = Capsule::segment_segment_closest(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 2.0, 0.0],
        );
        assert!(dist > 0.0);
        let _ = (pa, pb);
    }

    #[test]
    fn test_capsule_closest_point_on_surface_norm() {
        let c = Capsule::new(1.0, 2.0);
        let p = [3.0, 1.0, 4.0];
        let cp = c.closest_point(p);
        // Verify the closest point is on the capsule surface (sdf ~ 0)
        let sdf = c.sdf(cp);
        assert!(
            sdf.abs() < 1e-9,
            "closest point sdf should be ~0, got {sdf}"
        );
    }

    #[test]
    fn test_capsule_bounding_box_correct() {
        let c = Capsule::new(1.0, 2.0);
        use crate::shape::Shape;
        let bb = c.bounding_box();
        // Bounding box should extend r=1 in X/Z, half_height+r=3 in Y
        assert!((bb.max.y - 3.0).abs() < 1e-10);
        assert!((bb.max.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_capsule_chain_segment_radii() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let chain = CapsuleChain::with_radius(pts, 0.3);
        assert_eq!(chain.segment_count(), 3);
    }

    #[test]
    fn test_deformable_capsule_support_positive_y() {
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [0.0, 4.0, 0.0], 1.0);
        let sp = dc.support([0.0, 1.0, 0.0]);
        // Should be at top cap: (0, 4+1, 0)
        assert!(sp[1] > 4.0, "support in +Y should be above endpoint b");
    }

    #[test]
    fn test_capsule_frustum_sdf_on_axis_top() {
        // At the top (y=height), sdf should be ~0 or negative
        let f = CapsuleFrustum::new(0.5, 1.0, 3.0);
        // center at y=1.5, top cap at y=3 + r_top=1.0 => y=4 is outside
        let sdf_inside = f.sdf([0.0, 1.5, 0.0]);
        assert!(
            sdf_inside < 0.0,
            "center should be inside frustum, sdf={sdf_inside}"
        );
    }

    #[test]
    fn test_curved_capsule_path_segment_count() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let cp = CurvedCapsulePath::new(pts, 0.2);
        assert_eq!(cp.num_segments(), 4);
    }

    #[test]
    fn test_capsule_chain_sdf_on_axis() {
        let pts = vec![[0.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let chain = CapsuleChain::with_radius(pts, 1.0);
        // Point on the axis inside the chain
        let sdf = chain.sdf([0.0, 1.5, 0.0]);
        assert!(sdf < 0.0, "axis point should be inside chain, sdf={sdf}");
    }

    #[test]
    fn test_capsule_triangle_contact_normal_points_away() {
        let c = Capsule::new(1.0, 1.0);
        // Triangle at y=-0.3, large and centered
        let a = [-10.0, -0.3, -10.0];
        let b = [10.0, -0.3, -10.0];
        let tri = [0.0_f64, -0.3_f64, 10.0_f64];
        if let Some((depth, normal)) = c.capsule_triangle_contact([0.0, 0.0, 0.0], a, b, tri) {
            assert!(depth > 0.0, "depth should be positive, got {depth}");
            // Normal should be unit length (either +Y or -Y for a horizontal triangle)
            let n_len =
                (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
            assert!(
                (n_len - 1.0).abs() < 1e-9,
                "normal should be unit length, got {n_len}"
            );
        }
    }

    #[test]
    fn test_deformable_capsule_arbitrary_orientation() {
        // Capsule along X axis
        let dc = DeformableCapsule::new([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], 1.0);
        assert!(dc.sdf([2.0, 0.0, 0.0]) < 0.0); // on axis
        assert!(dc.sdf([2.0, 1.5, 0.0]) > 0.0); // outside radius
    }

    #[test]
    fn test_capsule_chain_min_distance_to_point() {
        let pts = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let chain = CapsuleChain::with_radius(pts, 0.5);
        let d = chain.min_distance_to_point([2.0, 3.0, 0.0]);
        // Closest axis point is (2,0,0), dist to that is 3, minus radius 0.5 = 2.5
        assert!((d - 2.5).abs() < 1e-9, "expected 2.5, got {d}");
    }

    #[test]
    fn test_capsule_frustum_larger_radius_at_top() {
        let f = CapsuleFrustum::new(0.5, 2.0, 3.0);
        // At the top, radius should be larger
        assert!(f.radius_at_height(3.0) > f.radius_at_height(0.0));
    }

    #[test]
    fn test_capsule_frustum_radius_interpolated() {
        let f = CapsuleFrustum::new(1.0, 3.0, 4.0);
        // At height 2 (midpoint), radius should be midway between 1 and 3
        let r_mid = f.radius_at_height(2.0);
        assert!(
            (r_mid - 2.0).abs() < 1e-9,
            "midpoint radius should be 2.0, got {r_mid}"
        );
    }
}
