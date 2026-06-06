//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// A rotational sweep (lathe) — rotate a 2-D profile around the Y axis.
///
/// The profile is specified as `(r, y)` pairs where `r` is the radial distance
/// from the Y axis.  The resulting solid of revolution is discretised into
/// `segments` triangular strips.
#[derive(Debug, Clone)]
pub struct RotationalSweep {
    /// Profile: `(radius, height)` pairs.
    pub profile: Vec<[f64; 2]>,
    /// Number of angular segments (resolution).
    pub segments: usize,
}
impl RotationalSweep {
    /// Create a new rotational sweep.
    pub fn new(profile: Vec<[f64; 2]>, segments: usize) -> Self {
        Self {
            profile,
            segments: segments.max(3),
        }
    }
    /// Conservative AABB.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        if self.profile.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let max_r = self
            .profile
            .iter()
            .map(|p| p[0].abs())
            .fold(0.0f64, f64::max);
        let min_y = self
            .profile
            .iter()
            .map(|p| p[1])
            .fold(f64::INFINITY, f64::min);
        let max_y = self
            .profile
            .iter()
            .map(|p| p[1])
            .fold(f64::NEG_INFINITY, f64::max);
        ([-max_r, min_y, -max_r], [max_r, max_y, max_r])
    }
    /// Approximate volume via the disk/washer method (trapezoidal integration).
    pub fn volume(&self) -> f64 {
        if self.profile.len() < 2 {
            return 0.0;
        }
        let mut vol = 0.0f64;
        for i in 0..self.profile.len() - 1 {
            let [r0, y0] = self.profile[i];
            let [r1, y1] = self.profile[i + 1];
            let dy = (y1 - y0).abs();
            let avg_area = std::f64::consts::PI * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
            vol += avg_area * dy;
        }
        vol
    }
    /// Approximate surface area (lateral only, no caps).
    pub fn lateral_surface_area(&self) -> f64 {
        if self.profile.len() < 2 {
            return 0.0;
        }
        let mut area = 0.0f64;
        for i in 0..self.profile.len() - 1 {
            let [r0, y0] = self.profile[i];
            let [r1, y1] = self.profile[i + 1];
            let dr = r1 - r0;
            let dy = y1 - y0;
            let slant = (dr * dr + dy * dy).sqrt();
            area += std::f64::consts::PI * (r0 + r1) * slant;
        }
        area
    }
    /// Generate vertex positions for the mesh.
    ///
    /// Returns a `Vec` of `[f64; 3]` positions.
    pub fn vertices(&self) -> Vec<[f64; 3]> {
        let segs = self.segments;
        let mut verts = Vec::with_capacity(self.profile.len() * segs);
        for &[r, y] in &self.profile {
            for s in 0..segs {
                let angle = 2.0 * std::f64::consts::PI * s as f64 / segs as f64;
                verts.push([r * angle.cos(), y, r * angle.sin()]);
            }
        }
        verts
    }
}
/// Result of a linear cast query.
#[derive(Debug, Clone)]
pub struct LinearCastResult {
    /// Time of impact in `[0, 1]`.
    pub toi: f64,
    /// Approximate contact point at time of impact.
    pub contact_point: [f64; 3],
    /// Approximate contact normal (from B towards A).
    pub normal: [f64; 3],
}
/// A sphere swept along a straight line segment — effectively a capsule.
#[derive(Debug, Clone)]
pub struct SweptSphere {
    /// World-space centre at the start of the sweep.
    pub center_start: [f64; 3],
    /// World-space centre at the end of the sweep.
    pub center_end: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}
impl SweptSphere {
    /// Create a new `SweptSphere`.
    pub fn new(center_start: [f64; 3], center_end: [f64; 3], radius: f64) -> Self {
        Self {
            center_start,
            center_end,
            radius,
        }
    }
    /// Compute a tight axis-aligned bounding box enclosing the entire sweep.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let r = self.radius;
        let min = [
            self.center_start[0].min(self.center_end[0]) - r,
            self.center_start[1].min(self.center_end[1]) - r,
            self.center_start[2].min(self.center_end[2]) - r,
        ];
        let max = [
            self.center_start[0].max(self.center_end[0]) + r,
            self.center_start[1].max(self.center_end[1]) + r,
            self.center_start[2].max(self.center_end[2]) + r,
        ];
        (min, max)
    }
    /// Length of the sweep path.
    pub fn sweep_length(&self) -> f64 {
        len3(sub3(self.center_end, self.center_start))
    }
    /// Centre position at parametric time `t` in `[0, 1]`.
    pub fn center_at(&self, t: f64) -> [f64; 3] {
        lerp3(self.center_start, self.center_end, t)
    }
    /// Sweep direction (unnormalized).
    pub fn direction(&self) -> [f64; 3] {
        sub3(self.center_end, self.center_start)
    }
    /// Surface area of the swept volume (capsule surface area).
    pub fn surface_area(&self) -> f64 {
        let l = self.sweep_length();
        2.0 * std::f64::consts::PI * self.radius * l
            + 4.0 * std::f64::consts::PI * self.radius * self.radius
    }
    /// Volume of the swept volume (capsule volume).
    pub fn volume(&self) -> f64 {
        let l = self.sweep_length();
        let r = self.radius;
        std::f64::consts::PI * r * r * l + (4.0 / 3.0) * std::f64::consts::PI * r * r * r
    }
    /// Ray vs swept sphere (capsule) intersection.
    ///
    /// Returns the smallest non-negative *t* such that the ray
    /// `ray_origin + t * ray_dir` touches the capsule surface, or `None` on
    /// miss.
    ///
    /// The capsule is the Minkowski sum of the line segment
    /// `[center_start, center_end]` with a ball of `radius`.
    pub fn ray_intersect(&self, ray_origin: [f64; 3], ray_dir: [f64; 3]) -> Option<f64> {
        let pa = self.center_start;
        let pb = self.center_end;
        let r = self.radius;
        let d = sub3(pb, pa);
        let ro = sub3(ray_origin, pa);
        let dd = dot3(d, d);
        let rd = dot3(ray_dir, d);
        let ro_d = dot3(ro, d);
        let ro_ro = dot3(ro, ro);
        let rd_rd = dot3(ray_dir, ray_dir);
        let ro_rd = dot3(ro, ray_dir);
        let a = rd_rd - rd * rd / dd;
        let b = 2.0 * (ro_rd - ro_d * rd / dd);
        let c = ro_ro - ro_d * ro_d / dd - r * r;
        let mut t_min = f64::INFINITY;
        let disc = b * b - 4.0 * a * c;
        if disc >= 0.0 && a.abs() > 1e-14 {
            let sq = disc.sqrt();
            for &sign in &[-1.0_f64, 1.0_f64] {
                let t = (-b + sign * sq) / (2.0 * a);
                if t >= 0.0 {
                    let proj = (ro_d + t * rd) / dd;
                    if (0.0..=1.0).contains(&proj) && t < t_min {
                        t_min = t;
                    }
                }
            }
        }
        {
            let qa = rd_rd;
            let qb = 2.0 * ro_rd;
            let qc = ro_ro - r * r;
            let disc_a = qb * qb - 4.0 * qa * qc;
            if disc_a >= 0.0 {
                let sq = disc_a.sqrt();
                for &sign in &[-1.0_f64, 1.0_f64] {
                    let t = (-qb + sign * sq) / (2.0 * qa);
                    if t >= 0.0 {
                        let proj = (ro_d + t * rd) / dd;
                        if proj <= 0.0 && t < t_min {
                            t_min = t;
                        }
                    }
                }
            }
        }
        {
            let ro_b = sub3(ray_origin, pb);
            let ro_b_ro_b = dot3(ro_b, ro_b);
            let ro_b_rd = dot3(ro_b, ray_dir);
            let qa = rd_rd;
            let qb = 2.0 * ro_b_rd;
            let qc = ro_b_ro_b - r * r;
            let disc_b = qb * qb - 4.0 * qa * qc;
            if disc_b >= 0.0 {
                let sq = disc_b.sqrt();
                for &sign in &[-1.0_f64, 1.0_f64] {
                    let t = (-qb + sign * sq) / (2.0 * qa);
                    if t >= 0.0 {
                        let proj = (ro_d + t * rd) / dd;
                        if proj >= 1.0 && t < t_min {
                            t_min = t;
                        }
                    }
                }
            }
        }
        if t_min.is_finite() { Some(t_min) } else { None }
    }
}
/// An oriented bounding box swept linearly along a displacement vector.
///
/// The box is described by its centre, three orientation axes (rows of a 3×3
/// rotation matrix), and half-extents along each local axis.  The swept volume
/// is conservative: AABB of start + end OBB.
#[derive(Debug, Clone)]
pub struct SweptObb {
    /// OBB centre at the start.
    pub center_start: [f64; 3],
    /// OBB orientation: three unit axis vectors.
    pub axes: [[f64; 3]; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
    /// Displacement vector (center_end = center_start + displacement).
    pub displacement: [f64; 3],
}
impl SweptObb {
    /// Create a new `SweptObb`.
    pub fn new(
        center_start: [f64; 3],
        axes: [[f64; 3]; 3],
        half_extents: [f64; 3],
        displacement: [f64; 3],
    ) -> Self {
        Self {
            center_start,
            axes,
            half_extents,
            displacement,
        }
    }
    /// Centre at the end of the sweep.
    pub fn center_end(&self) -> [f64; 3] {
        add3(self.center_start, self.displacement)
    }
    /// Conservative AABB enclosing the OBB at both start and end.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut world_min = [f64::INFINITY; 3];
        let mut world_max = [f64::NEG_INFINITY; 3];
        for &center in &[self.center_start, self.center_end()] {
            for k in 0..3 {
                let mut r = 0.0f64;
                for j in 0..3 {
                    r += self.half_extents[j] * self.axes[j][k].abs();
                }
                if center[k] - r < world_min[k] {
                    world_min[k] = center[k] - r;
                }
                if center[k] + r > world_max[k] {
                    world_max[k] = center[k] + r;
                }
            }
        }
        (world_min, world_max)
    }
    /// Volume of the OBB (static).
    pub fn volume(&self) -> f64 {
        8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
    }
    /// Support point of the static OBB in direction `dir`.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let mut result = self.center_start;
        for j in 0..3 {
            let s = if dot3(self.axes[j], dir) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            result = add3(result, scale3(self.axes[j], s * self.half_extents[j]));
        }
        result
    }
}
/// A capsule swept along a straight line defined by start and end positions.
///
/// The capsule itself is axis-aligned along the Y axis with given `radius`
/// and `half_height`. The swept volume is the Minkowski sum of the capsule
/// with the line segment.
#[derive(Debug, Clone)]
pub struct SweptCapsule {
    /// World-space position of the capsule center at start.
    pub position_start: [f64; 3],
    /// World-space position of the capsule center at end.
    pub position_end: [f64; 3],
    /// Capsule radius.
    pub radius: f64,
    /// Capsule half-height (distance from centre to each hemisphere centre).
    pub half_height: f64,
}
impl SweptCapsule {
    /// Create a new `SweptCapsule`.
    pub fn new(
        position_start: [f64; 3],
        position_end: [f64; 3],
        radius: f64,
        half_height: f64,
    ) -> Self {
        Self {
            position_start,
            position_end,
            radius,
            half_height,
        }
    }
    /// Conservative AABB enclosing the swept capsule.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let r = self.radius;
        let h = self.half_height;
        let expand = [r, h + r, r];
        let mut world_min = [f64::INFINITY; 3];
        let mut world_max = [f64::NEG_INFINITY; 3];
        for &pos in &[self.position_start, self.position_end] {
            for k in 0..3 {
                let lo = pos[k] - expand[k];
                let hi = pos[k] + expand[k];
                if lo < world_min[k] {
                    world_min[k] = lo;
                }
                if hi > world_max[k] {
                    world_max[k] = hi;
                }
            }
        }
        (world_min, world_max)
    }
    /// Position at parametric time `t`.
    pub fn position_at(&self, t: f64) -> [f64; 3] {
        lerp3(self.position_start, self.position_end, t)
    }
    /// Volume of the capsule (static, not swept).
    pub fn capsule_volume(&self) -> f64 {
        let r = self.radius;
        let h = 2.0 * self.half_height;
        std::f64::consts::PI * r * r * h + (4.0 / 3.0) * std::f64::consts::PI * r * r * r
    }
    /// Time of impact against a static sphere using conservative advancement.
    ///
    /// Returns `t` in `[0, 1]` or `None`.
    pub fn toi_vs_sphere(&self, sphere_center: [f64; 3], sphere_radius: f64) -> Option<f64> {
        let effective_radius = self.radius + self.half_height;
        let combined = effective_radius + sphere_radius;
        let vel = sub3(self.position_end, self.position_start);
        let rel = sub3(self.position_start, sphere_center);
        let a = dot3(vel, vel);
        let b = 2.0 * dot3(rel, vel);
        let c = dot3(rel, rel) - combined * combined;
        if a < 1e-14 {
            return if c <= 0.0 { Some(0.0) } else { None };
        }
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        let t1 = (-b - sq) / (2.0 * a);
        let t2 = (-b + sq) / (2.0 * a);
        let t = if t1 >= 0.0 { t1 } else { t2 };
        if (0.0..=1.0).contains(&t) {
            Some(t)
        } else {
            None
        }
    }
}
/// A box swept along a path defined by a start and end transform.
///
/// The conservative AABB encloses the box at both the start and end poses.
#[derive(Debug, Clone)]
pub struct SweptBox {
    /// Homogeneous transform (row-major 4x4) at the start of the sweep.
    pub transform_start: [[f64; 4]; 4],
    /// Homogeneous transform (row-major 4x4) at the end of the sweep.
    pub transform_end: [[f64; 4]; 4],
    /// Half-extents of the box along its local axes.
    pub half_extents: [f64; 3],
}
impl SweptBox {
    /// Create a new `SweptBox`.
    pub fn new(
        transform_start: [[f64; 4]; 4],
        transform_end: [[f64; 4]; 4],
        half_extents: [f64; 3],
    ) -> Self {
        Self {
            transform_start,
            transform_end,
            half_extents,
        }
    }
    /// Conservative AABB enclosing the box at both the start and end poses.
    ///
    /// All 8 corners of the box are transformed into world space at each
    /// endpoint and the result is the union of the two world-space AABBs.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let hx = self.half_extents[0];
        let hy = self.half_extents[1];
        let hz = self.half_extents[2];
        let corners: [[f64; 3]; 8] = [
            [-hx, -hy, -hz],
            [hx, -hy, -hz],
            [-hx, hy, -hz],
            [hx, hy, -hz],
            [-hx, -hy, hz],
            [hx, -hy, hz],
            [-hx, hy, hz],
            [hx, hy, hz],
        ];
        let mut world_min = [f64::INFINITY; 3];
        let mut world_max = [f64::NEG_INFINITY; 3];
        for &m in &[self.transform_start, self.transform_end] {
            for &lc in &corners {
                let wc = transform_point(m, lc);
                for k in 0..3 {
                    if wc[k] < world_min[k] {
                        world_min[k] = wc[k];
                    }
                    if wc[k] > world_max[k] {
                        world_max[k] = wc[k];
                    }
                }
            }
        }
        (world_min, world_max)
    }
    /// Conservative AABB including `n_samples` intermediate poses.
    pub fn aabb_sampled(&self, n_samples: usize) -> ([f64; 3], [f64; 3]) {
        let hx = self.half_extents[0];
        let hy = self.half_extents[1];
        let hz = self.half_extents[2];
        let corners: [[f64; 3]; 8] = [
            [-hx, -hy, -hz],
            [hx, -hy, -hz],
            [-hx, hy, -hz],
            [hx, hy, -hz],
            [-hx, -hy, hz],
            [hx, -hy, hz],
            [-hx, hy, hz],
            [hx, hy, hz],
        ];
        let mut world_min = [f64::INFINITY; 3];
        let mut world_max = [f64::NEG_INFINITY; 3];
        let steps = n_samples.max(2);
        for i in 0..steps {
            let t = i as f64 / (steps - 1) as f64;
            let m = lerp_matrix(self.transform_start, self.transform_end, t);
            for &lc in &corners {
                let wc = transform_point(m, lc);
                for k in 0..3 {
                    if wc[k] < world_min[k] {
                        world_min[k] = wc[k];
                    }
                    if wc[k] > world_max[k] {
                        world_max[k] = wc[k];
                    }
                }
            }
        }
        (world_min, world_max)
    }
    /// Volume of the box (static).
    pub fn box_volume(&self) -> f64 {
        8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
    }
    /// Translation at the start pose.
    pub fn start_translation(&self) -> [f64; 3] {
        [
            self.transform_start[0][3],
            self.transform_start[1][3],
            self.transform_start[2][3],
        ]
    }
    /// Translation at the end pose.
    pub fn end_translation(&self) -> [f64; 3] {
        [
            self.transform_end[0][3],
            self.transform_end[1][3],
            self.transform_end[2][3],
        ]
    }
    /// Displacement vector from start to end translation.
    pub fn displacement(&self) -> [f64; 3] {
        sub3(self.end_translation(), self.start_translation())
    }
}
/// A shape extruded linearly along a vector — a prism swept volume.
///
/// The base polygon is defined in the XY plane and the extrusion direction
/// is given by `sweep_vec`.
#[derive(Debug, Clone)]
pub struct LinearExtrusion {
    /// 2-D profile points (X, Y components).
    pub profile: Vec<[f64; 2]>,
    /// Extrusion vector in 3-D space.
    pub sweep_vec: [f64; 3],
}
impl LinearExtrusion {
    /// Create a new linear extrusion.
    pub fn new(profile: Vec<[f64; 2]>, sweep_vec: [f64; 3]) -> Self {
        Self { profile, sweep_vec }
    }
    /// Conservative AABB enclosing the extruded prism.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        if self.profile.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for &[x, y] in &self.profile {
            mn[0] = mn[0].min(x);
            mx[0] = mx[0].max(x);
            mn[1] = mn[1].min(y);
            mx[1] = mx[1].max(y);
            mn[2] = mn[2].min(0.0);
            mx[2] = mx[2].max(0.0);
            let sx = x + self.sweep_vec[0];
            let sy = y + self.sweep_vec[1];
            let sz = self.sweep_vec[2];
            mn[0] = mn[0].min(sx);
            mx[0] = mx[0].max(sx);
            mn[1] = mn[1].min(sy);
            mx[1] = mx[1].max(sy);
            mn[2] = mn[2].min(sz);
            mx[2] = mx[2].max(sz);
        }
        (mn, mx)
    }
    /// Approximate volume: profile area × extrusion length.
    pub fn volume(&self) -> f64 {
        let area = self.profile_area();
        let len = len3(self.sweep_vec);
        area * len
    }
    /// Area of the 2-D profile using the shoelace formula.
    pub fn profile_area(&self) -> f64 {
        let n = self.profile.len();
        if n < 3 {
            return 0.0;
        }
        let mut signed = 0.0f64;
        for i in 0..n {
            let [x0, y0] = self.profile[i];
            let [x1, y1] = self.profile[(i + 1) % n];
            signed += x0 * y1 - x1 * y0;
        }
        (signed * 0.5).abs()
    }
    /// Perimeter of the 2-D profile.
    pub fn profile_perimeter(&self) -> f64 {
        let n = self.profile.len();
        if n < 2 {
            return 0.0;
        }
        let mut perim = 0.0f64;
        for i in 0..n {
            let [x0, y0] = self.profile[i];
            let [x1, y1] = self.profile[(i + 1) % n];
            let dx = x1 - x0;
            let dy = y1 - y0;
            perim += (dx * dx + dy * dy).sqrt();
        }
        perim
    }
    /// Surface area of the extruded prism.
    ///
    /// = 2 * profile_area + perimeter * extrusion_length
    pub fn surface_area(&self) -> f64 {
        let area = self.profile_area();
        let perim = self.profile_perimeter();
        let len = len3(self.sweep_vec);
        2.0 * area + perim * len
    }
}
/// A shape defined by its start and end AABB — the swept union.
///
/// This represents the conservative volume swept by an AABB moving linearly
/// from `start` to `end` (i.e., the union of the two AABBs and everything in
/// between along each axis).
#[derive(Debug, Clone)]
pub struct SweptAabb {
    /// AABB at the start: `(min, max)`.
    pub start_min: [f64; 3],
    /// AABB at the start: max corner.
    pub start_max: [f64; 3],
    /// AABB at the end: min corner.
    pub end_min: [f64; 3],
    /// AABB at the end: max corner.
    pub end_max: [f64; 3],
}
impl SweptAabb {
    /// Create a new `SweptAabb`.
    pub fn new(
        start_min: [f64; 3],
        start_max: [f64; 3],
        end_min: [f64; 3],
        end_max: [f64; 3],
    ) -> Self {
        Self {
            start_min,
            start_max,
            end_min,
            end_max,
        }
    }
    /// Conservative AABB covering the entire sweep.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for k in 0..3 {
            mn[k] = self.start_min[k].min(self.end_min[k]);
            mx[k] = self.start_max[k].max(self.end_max[k]);
        }
        (mn, mx)
    }
    /// Test whether a static point is inside the swept volume AABB.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        let (mn, mx) = self.aabb();
        (0..3).all(|k| p[k] >= mn[k] && p[k] <= mx[k])
    }
    /// Displacement of the AABB centre from start to end.
    pub fn displacement(&self) -> [f64; 3] {
        let start_center = scale3(add3(self.start_min, self.start_max), 0.5);
        let end_center = scale3(add3(self.end_min, self.end_max), 0.5);
        sub3(end_center, start_center)
    }
}
/// Result of a CCD time-of-impact query.
#[derive(Debug, Clone)]
pub struct CcdResult {
    /// Normalised time of impact in `[0, 1]`.
    pub toi: f64,
    /// Approximate contact point.
    pub contact_point: [f64; 3],
}
