//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use crate::math::{Quat, Real, Vec3};

use super::functions::{quat_mul, quat_rotate};

/// Builder for PhysicsConfig with a fluent API.
#[derive(Debug, Clone)]
pub struct PhysicsConfigBuilder {
    pub(super) config: PhysicsConfig,
}
impl PhysicsConfigBuilder {
    /// Start building a PhysicsConfig from defaults.
    pub fn new() -> Self {
        Self {
            config: PhysicsConfig::default(),
        }
    }
    /// Set the gravity vector.
    pub fn gravity(mut self, g: Vec3) -> Self {
        self.config.gravity = g;
        self
    }
    /// Set the number of solver iterations.
    pub fn solver_iterations(mut self, n: u32) -> Self {
        self.config.solver_iterations = n;
        self
    }
    /// Set the linear sleep threshold.
    pub fn linear_sleep_threshold(mut self, t: Real) -> Self {
        self.config.linear_sleep_threshold = t;
        self
    }
    /// Set the angular sleep threshold.
    pub fn angular_sleep_threshold(mut self, t: Real) -> Self {
        self.config.angular_sleep_threshold = t;
        self
    }
    /// Set the time before sleep.
    pub fn time_before_sleep(mut self, t: Real) -> Self {
        self.config.time_before_sleep = t;
        self
    }
    /// Enable or disable CCD.
    pub fn ccd_enabled(mut self, enabled: bool) -> Self {
        self.config.ccd_enabled = enabled;
        self
    }
    /// Build the PhysicsConfig.
    pub fn build(self) -> PhysicsConfig {
        self.config
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// An infinite plane defined by normal and offset.
#[derive(Debug, Clone)]
pub struct Plane {
    /// Outward-facing unit normal.
    pub normal: Vec3,
    /// Signed distance from the origin to the plane (along the normal).
    pub offset: Real,
}
impl Plane {
    /// Create a plane from normal and offset. Normal is normalised internally.
    pub fn new(normal: Vec3, offset: Real) -> Self {
        let n = normal
            .try_normalize(1e-12)
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
        Self { normal: n, offset }
    }
    /// Create a plane passing through `point` with the given normal.
    pub fn from_point_normal(point: Vec3, normal: Vec3) -> Self {
        let n = normal
            .try_normalize(1e-12)
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
        let offset = n.dot(&point);
        Self { normal: n, offset }
    }
    /// Signed distance from `point` to the plane.
    pub fn signed_distance(&self, point: &Vec3) -> Real {
        self.normal.dot(point) - self.offset
    }
    /// Project a point onto the plane (closest point on the plane).
    pub fn project_point(&self, point: &Vec3) -> Vec3 {
        *point - self.normal * self.signed_distance(point)
    }
    /// Ray–plane intersection. Returns the t-value or `None` if parallel.
    pub fn ray_intersect(&self, ray: &Ray) -> Option<Real> {
        let denom = self.normal.dot(&ray.direction);
        if denom.abs() < 1e-12 {
            return None;
        }
        let t = (self.offset - self.normal.dot(&ray.origin)) / denom;
        if t >= 0.0 && t <= ray.max_t {
            Some(t)
        } else {
            None
        }
    }
    /// Test whether two spheres straddle the plane.
    pub fn sphere_straddles(&self, sphere: &Sphere) -> bool {
        self.signed_distance(&sphere.center).abs() <= sphere.radius
    }
}
/// A sphere in 3D space.
#[derive(Debug, Clone)]
pub struct Sphere {
    /// Center of the sphere.
    pub center: Vec3,
    /// Radius of the sphere.
    pub radius: Real,
}
impl Sphere {
    /// Create a new sphere.
    pub fn new(center: Vec3, radius: Real) -> Self {
        Self { center, radius }
    }
    /// Axis-aligned bounding box enclosing this sphere.
    pub fn aabb(&self) -> Aabb {
        let r = Vec3::new(self.radius, self.radius, self.radius);
        Aabb {
            min: self.center - r,
            max: self.center + r,
        }
    }
    /// Test if `point` is inside or on the surface of the sphere.
    pub fn contains_point(&self, point: &Vec3) -> bool {
        (point - self.center).norm_squared() <= self.radius * self.radius
    }
    /// Test if this sphere overlaps another sphere.
    pub fn overlaps_sphere(&self, other: &Sphere) -> bool {
        let r_sum = self.radius + other.radius;
        (self.center - other.center).norm_squared() <= r_sum * r_sum
    }
    /// Test if this sphere overlaps an AABB.
    pub fn overlaps_aabb(&self, aabb: &Aabb) -> bool {
        aabb.sq_distance_to_point(&self.center) <= self.radius * self.radius
    }
    /// Ray–sphere intersection. Returns the entry t-value or `None`.
    pub fn ray_intersect(&self, ray: &Ray) -> Option<Real> {
        let oc = ray.origin - self.center;
        let b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;
        let disc = b * b - c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t = -b - sqrt_disc;
        if t >= 0.0 && t <= ray.max_t {
            return Some(t);
        }
        let t2 = -b + sqrt_disc;
        if t2 >= 0.0 && t2 <= ray.max_t {
            return Some(t2);
        }
        None
    }
    /// Volume of the sphere.
    pub fn volume(&self) -> Real {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius.powi(3)
    }
    /// Surface area of the sphere.
    pub fn surface_area(&self) -> Real {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }
    /// Minimum enclosing sphere of a set of points (Ritter's algorithm).
    pub fn bounding_sphere(points: &[Vec3]) -> Option<Sphere> {
        if points.is_empty() {
            return None;
        }
        if points.len() == 1 {
            return Some(Sphere::new(points[0], 0.0));
        }
        let mut center = points[0];
        let mut radius = 0.0f64;
        for &p in points.iter().skip(1) {
            let d = (p - center).norm();
            if d > radius {
                radius = (radius + d) * 0.5;
                center = center + (p - center) * ((d - radius) / d + 0.5);
            }
        }
        for &p in points {
            let d = (p - center).norm();
            if d > radius {
                radius = d;
            }
        }
        Some(Sphere::new(center, radius))
    }
}
/// Axis-aligned bounding box for broad-phase collision detection.
#[derive(Debug, Clone)]
pub struct Aabb {
    /// Minimum corner of the bounding box.
    pub min: Vec3,
    /// Maximum corner of the bounding box.
    pub max: Vec3,
}
impl Aabb {
    /// Create a new AABB from min and max corners.
    pub fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }
    /// Check whether this AABB intersects another.
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }
    /// Merge this AABB with another, returning the smallest AABB containing both.
    pub fn merge(&self, other: &Aabb) -> Self {
        Self {
            min: self.min.inf(&other.min),
            max: self.max.sup(&other.max),
        }
    }
    /// Check if a point is inside the AABB.
    pub fn contains_point(&self, point: &Vec3) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }
    /// Expand the AABB by a margin in all directions.
    pub fn expand(&self, margin: Real) -> Self {
        let v = Vec3::new(margin, margin, margin);
        Self {
            min: self.min - v,
            max: self.max + v,
        }
    }
    /// Center of the AABB.
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }
    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> Vec3 {
        (self.max - self.min) * 0.5
    }
    /// Surface area of the AABB.
    pub fn surface_area(&self) -> Real {
        let d = self.max - self.min;
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }
    /// Volume of the AABB.
    pub fn volume(&self) -> Real {
        let d = self.max - self.min;
        d.x * d.y * d.z
    }
}
impl Aabb {
    /// Create an AABB from a center point and half-extents.
    pub fn from_center_half_extents(center: Vec3, half_extents: Vec3) -> Self {
        Self {
            min: center - half_extents,
            max: center + half_extents,
        }
    }
    /// Create an AABB enclosing a set of points.
    pub fn from_points(points: &[Vec3]) -> Option<Self> {
        if points.is_empty() {
            return None;
        }
        let mut min = points[0];
        let mut max = points[0];
        for p in points.iter().skip(1) {
            min = min.inf(p);
            max = max.sup(p);
        }
        Some(Self { min, max })
    }
    /// Longest axis (0=x, 1=y, 2=z).
    pub fn longest_axis(&self) -> usize {
        let d = self.max - self.min;
        if d.x >= d.y && d.x >= d.z {
            0
        } else if d.y >= d.z {
            1
        } else {
            2
        }
    }
    /// Diagonal length of the AABB.
    pub fn diagonal(&self) -> f64 {
        (self.max - self.min).norm()
    }
}
impl Aabb {
    /// Test whether `ray` intersects this AABB using the slab method.
    ///
    /// Returns `Some(hit)` with the entry/exit t values if the ray hits,
    /// or `None` if it misses.  Handles degenerate (zero-direction) rays.
    pub fn ray_intersect(&self, ray: &Ray) -> Option<RayAabbHit> {
        let mut t_min = 0.0_f64;
        let mut t_max = ray.max_t;
        let origin = [ray.origin.x, ray.origin.y, ray.origin.z];
        let dir = [ray.direction.x, ray.direction.y, ray.direction.z];
        let bmin = [self.min.x, self.min.y, self.min.z];
        let bmax = [self.max.x, self.max.y, self.max.z];
        for i in 0..3 {
            if dir[i].abs() < 1e-12 {
                if origin[i] < bmin[i] || origin[i] > bmax[i] {
                    return None;
                }
            } else {
                let inv_d = 1.0 / dir[i];
                let t0 = (bmin[i] - origin[i]) * inv_d;
                let t1 = (bmax[i] - origin[i]) * inv_d;
                let (t0, t1) = if inv_d < 0.0 { (t1, t0) } else { (t0, t1) };
                t_min = t_min.max(t0);
                t_max = t_max.min(t1);
                if t_max < t_min {
                    return None;
                }
            }
        }
        Some(RayAabbHit { t_min, t_max })
    }
    /// Test ray against AABB and return just the entry t (simpler API).
    pub fn ray_hit_t(&self, ray: &Ray) -> Option<Real> {
        self.ray_intersect(ray).map(|h| h.t_min)
    }
    /// Compute the AABB of the swept volume as the box moves from `transform_a`
    /// to `transform_b` (uses the merge of two oriented boxes).
    pub fn swept_aabb(half_extents: Vec3, from: &Transform, to: &Transform) -> Aabb {
        let aabb_from = Aabb::oriented_box_aabb(half_extents, from);
        let aabb_to = Aabb::oriented_box_aabb(half_extents, to);
        aabb_from.merge(&aabb_to)
    }
    /// Compute the AABB enclosing an oriented box (defined by half-extents and a transform).
    pub fn oriented_box_aabb(half_extents: Vec3, transform: &Transform) -> Aabb {
        let rot = transform.rotation.to_rotation_matrix();
        let m = rot.matrix();
        let wx = m[(0, 0)].abs() * half_extents.x
            + m[(0, 1)].abs() * half_extents.y
            + m[(0, 2)].abs() * half_extents.z;
        let wy = m[(1, 0)].abs() * half_extents.x
            + m[(1, 1)].abs() * half_extents.y
            + m[(1, 2)].abs() * half_extents.z;
        let wz = m[(2, 0)].abs() * half_extents.x
            + m[(2, 1)].abs() * half_extents.y
            + m[(2, 2)].abs() * half_extents.z;
        let half = Vec3::new(wx, wy, wz);
        Aabb {
            min: transform.position - half,
            max: transform.position + half,
        }
    }
    /// Compute the intersection (overlap) AABB of two boxes.
    ///
    /// Returns `None` if the boxes do not overlap.
    pub fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        let min = Vec3::new(
            self.min.x.max(other.min.x),
            self.min.y.max(other.min.y),
            self.min.z.max(other.min.z),
        );
        let max = Vec3::new(
            self.max.x.min(other.max.x),
            self.max.y.min(other.max.y),
            self.max.z.min(other.max.z),
        );
        if min.x <= max.x && min.y <= max.y && min.z <= max.z {
            Some(Aabb { min, max })
        } else {
            None
        }
    }
    /// Compute the squared distance from a point to the nearest point on the AABB surface.
    ///
    /// Returns 0 if the point is inside the AABB.
    pub fn sq_distance_to_point(&self, point: &Vec3) -> Real {
        let mut dist_sq = 0.0;
        for (p, lo, hi) in [
            (point.x, self.min.x, self.max.x),
            (point.y, self.min.y, self.max.y),
            (point.z, self.min.z, self.max.z),
        ] {
            if p < lo {
                dist_sq += (lo - p) * (lo - p);
            } else if p > hi {
                dist_sq += (p - hi) * (p - hi);
            }
        }
        dist_sq
    }
    /// Return the closest point on (or inside) the AABB to `point`.
    pub fn closest_point(&self, point: &Vec3) -> Vec3 {
        Vec3::new(
            point.x.clamp(self.min.x, self.max.x),
            point.y.clamp(self.min.y, self.max.y),
            point.z.clamp(self.min.z, self.max.z),
        )
    }
}
/// Transform representing position and orientation in 3D space.
#[derive(Debug, Clone)]
pub struct Transform {
    /// Position in world space.
    pub position: Vec3,
    /// Orientation as a unit quaternion.
    pub rotation: Quat,
}
impl Transform {
    /// Create a new transform with the given position and identity rotation.
    pub fn from_position(position: Vec3) -> Self {
        Self {
            position,
            rotation: Quat::identity(),
        }
    }
    /// Create a new transform with the given position and rotation.
    pub fn new(position: Vec3, rotation: Quat) -> Self {
        Self { position, rotation }
    }
    /// Transform a point from local space to world space.
    pub fn transform_point(&self, point: &Vec3) -> Vec3 {
        self.rotation * point + self.position
    }
    /// Transform a vector (direction) from local space to world space.
    pub fn transform_vector(&self, vector: &Vec3) -> Vec3 {
        self.rotation * vector
    }
    /// Compute the inverse transform.
    pub fn inverse(&self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            position: inv_rot * (-self.position),
            rotation: inv_rot,
        }
    }
    /// Compose this transform with another (self * other).
    pub fn compose(&self, other: &Transform) -> Self {
        Self {
            position: self.rotation * other.position + self.position,
            rotation: self.rotation * other.rotation,
        }
    }
}
impl Transform {
    /// Convert to a 4x4 transformation matrix (column-major).
    ///
    /// The matrix is `[R | t; 0 0 0 1]` where R is the rotation matrix.
    pub fn to_matrix4(&self) -> [[f64; 4]; 4] {
        let r = self.rotation.to_rotation_matrix();
        let m = r.matrix();
        [
            [m[(0, 0)], m[(1, 0)], m[(2, 0)], 0.0],
            [m[(0, 1)], m[(1, 1)], m[(2, 1)], 0.0],
            [m[(0, 2)], m[(1, 2)], m[(2, 2)], 0.0],
            [self.position.x, self.position.y, self.position.z, 1.0],
        ]
    }
    /// Extract Euler angles (roll, pitch, yaw) from the rotation.
    pub fn euler_angles(&self) -> (f64, f64, f64) {
        self.rotation.euler_angles()
    }
    /// Create a transform from an axis-angle rotation and a position.
    pub fn from_axis_angle(position: Vec3, axis: Vec3, angle: f64) -> Self {
        let axis_unit = nalgebra::Unit::new_normalize(axis);
        let rotation = Quat::from_axis_angle(&axis_unit, angle);
        Self { position, rotation }
    }
    /// Linear interpolation between two transforms.
    pub fn lerp(&self, other: &Transform, t: f64) -> Transform {
        let pos = self.position + (other.position - self.position) * t;
        let rot = self.rotation.slerp(&other.rotation, t);
        Transform {
            position: pos,
            rotation: rot,
        }
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// Result of a ray–AABB intersection test.
#[derive(Debug, Clone, Copy)]
pub struct RayAabbHit {
    /// Entry parameter (t_min).
    pub t_min: Real,
    /// Exit parameter (t_max).
    pub t_max: Real,
}
/// A ray in 3D space defined by an origin and a unit-direction vector.
#[derive(Debug, Clone)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: Vec3,
    /// Unit direction vector.
    pub direction: Vec3,
    /// Maximum ray length (use `f64::MAX` for infinite rays).
    pub max_t: Real,
}
impl Ray {
    /// Create a new ray. `direction` is normalised internally.
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        let dir_norm = direction
            .try_normalize(1e-12)
            .unwrap_or_else(|| Vec3::new(0.0, 0.0, 1.0));
        Self {
            origin,
            direction: dir_norm,
            max_t: Real::MAX,
        }
    }
    /// Create a ray with an explicit maximum extent.
    pub fn new_bounded(origin: Vec3, direction: Vec3, max_t: Real) -> Self {
        let dir_norm = direction
            .try_normalize(1e-12)
            .unwrap_or_else(|| Vec3::new(0.0, 0.0, 1.0));
        Self {
            origin,
            direction: dir_norm,
            max_t,
        }
    }
    /// Evaluate the point at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: Real) -> Vec3 {
        self.origin + self.direction * t
    }
    /// Create a ray from `from` toward `to`, with length = distance(from, to).
    pub fn from_to(from: Vec3, to: Vec3) -> Self {
        let d = to - from;
        let len = d.norm();
        let dir = if len > 1e-12 {
            d / len
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        Self {
            origin: from,
            direction: dir,
            max_t: len,
        }
    }
}
/// A triangle in 3D space.
#[derive(Debug, Clone)]
pub struct Triangle {
    /// First vertex.
    pub a: Vec3,
    /// Second vertex.
    pub b: Vec3,
    /// Third vertex.
    pub c: Vec3,
}
impl Triangle {
    /// Create a new triangle.
    pub fn new(a: Vec3, b: Vec3, c: Vec3) -> Self {
        Self { a, b, c }
    }
    /// Outward-facing normal (not normalised).
    pub fn normal_unnorm(&self) -> Vec3 {
        (self.b - self.a).cross(&(self.c - self.a))
    }
    /// Unit normal.
    pub fn normal(&self) -> Vec3 {
        self.normal_unnorm()
            .try_normalize(1e-12)
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0))
    }
    /// Area of the triangle.
    pub fn area(&self) -> Real {
        self.normal_unnorm().norm() * 0.5
    }
    /// Centroid of the triangle.
    pub fn centroid(&self) -> Vec3 {
        (self.a + self.b + self.c) / 3.0
    }
    /// Möller–Trumbore ray–triangle intersection.
    ///
    /// Returns `Some(t)` if the ray hits the front face, otherwise `None`.
    pub fn ray_intersect(&self, ray: &Ray) -> Option<Real> {
        let edge1 = self.b - self.a;
        let edge2 = self.c - self.a;
        let h = ray.direction.cross(&edge2);
        let det = edge1.dot(&h);
        if det.abs() < 1e-12 {
            return None;
        }
        let inv_det = 1.0 / det;
        let s = ray.origin - self.a;
        let u = inv_det * s.dot(&h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = s.cross(&edge1);
        let v = inv_det * ray.direction.dot(&q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = inv_det * edge2.dot(&q);
        if t >= 0.0 && t <= ray.max_t {
            Some(t)
        } else {
            None
        }
    }
    /// Return the barycentric coordinates (u, v, w) of a point projected
    /// onto the triangle plane.
    pub fn barycentric(&self, point: &Vec3) -> (Real, Real, Real) {
        let v0 = self.b - self.a;
        let v1 = self.c - self.a;
        let v2 = *point - self.a;
        let d00 = v0.dot(&v0);
        let d01 = v0.dot(&v1);
        let d11 = v1.dot(&v1);
        let d20 = v2.dot(&v0);
        let d21 = v2.dot(&v1);
        let denom = d00 * d11 - d01 * d01;
        if denom.abs() < 1e-20 {
            return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        }
        let v = (d11 * d20 - d01 * d21) / denom;
        let w = (d00 * d21 - d01 * d20) / denom;
        let u = 1.0 - v - w;
        (u, v, w)
    }
    /// Axis-aligned bounding box of the triangle.
    pub fn aabb(&self) -> Aabb {
        Aabb {
            min: self.a.inf(&self.b).inf(&self.c),
            max: self.a.sup(&self.b).sup(&self.c),
        }
    }
}
/// A capsule: a line segment swept by a sphere.
#[derive(Debug, Clone)]
pub struct Capsule {
    /// First endpoint of the central axis.
    pub a: Vec3,
    /// Second endpoint of the central axis.
    pub b: Vec3,
    /// Radius of the capsule.
    pub radius: Real,
}
impl Capsule {
    /// Create a new capsule.
    pub fn new(a: Vec3, b: Vec3, radius: Real) -> Self {
        Self { a, b, radius }
    }
    /// Create an upright capsule centered at `center` with given height and radius.
    pub fn upright(center: Vec3, half_height: Real, radius: Real) -> Self {
        Self {
            a: Vec3::new(center.x, center.y - half_height, center.z),
            b: Vec3::new(center.x, center.y + half_height, center.z),
            radius,
        }
    }
    /// Axis-aligned bounding box of the capsule.
    pub fn aabb(&self) -> Aabb {
        let r = Vec3::new(self.radius, self.radius, self.radius);
        let min = self.a.inf(&self.b) - r;
        let max = self.a.sup(&self.b) + r;
        Aabb { min, max }
    }
    /// Squared distance from a point to the capsule's central axis segment.
    pub fn sq_dist_point_to_axis(&self, point: &Vec3) -> Real {
        let ab = self.b - self.a;
        let ap = *point - self.a;
        let t = ap.dot(&ab) / ab.dot(&ab).max(1e-20);
        let t_clamped = t.clamp(0.0, 1.0);
        let closest = self.a + ab * t_clamped;
        (*point - closest).norm_squared()
    }
    /// Test if a point is inside the capsule.
    pub fn contains_point(&self, point: &Vec3) -> bool {
        self.sq_dist_point_to_axis(point) <= self.radius * self.radius
    }
    /// Volume of the capsule (cylinder + two hemispheres).
    pub fn volume(&self) -> Real {
        let h = (self.b - self.a).norm();
        std::f64::consts::PI * self.radius * self.radius * (h + (4.0 / 3.0) * self.radius)
    }
}
/// A 3-D transform stored as position + quaternion (xyzw) + uniform scale,
/// using plain arrays to avoid nalgebra dependencies outside `oxiphysics-core`.
#[derive(Debug, Clone)]
pub struct Transform3D {
    /// Translation (world-space position).
    pub position: [f64; 3],
    /// Rotation quaternion stored as `[x, y, z, w]` (unit quaternion).
    pub rotation: [f64; 4],
    /// Uniform scale factor.
    pub scale: f64,
}
impl Transform3D {
    /// Create an identity transform (no translation, no rotation, scale = 1).
    pub fn identity() -> Self {
        Self::default()
    }
    /// Create a transform with position only (identity rotation, scale = 1).
    pub fn from_position(position: [f64; 3]) -> Self {
        Self {
            position,
            ..Self::default()
        }
    }
    /// Apply the transform to a point: `scale * rotate(point) + position`.
    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let rotated = quat_rotate(self.rotation, point);
        [
            self.scale * rotated[0] + self.position[0],
            self.scale * rotated[1] + self.position[1],
            self.scale * rotated[2] + self.position[2],
        ]
    }
    /// Compute the inverse transform.
    ///
    /// For a transform `T = (p, q, s)`, the inverse is `(-q^-1(p)/s, q^-1, 1/s)`.
    pub fn inverse(&self) -> Self {
        let inv_scale = 1.0 / self.scale;
        let inv_rot = [
            -self.rotation[0],
            -self.rotation[1],
            -self.rotation[2],
            self.rotation[3],
        ];
        let neg_p = [
            -self.position[0] * inv_scale,
            -self.position[1] * inv_scale,
            -self.position[2] * inv_scale,
        ];
        let inv_pos = quat_rotate(inv_rot, neg_p);
        Self {
            position: inv_pos,
            rotation: inv_rot,
            scale: inv_scale,
        }
    }
    /// Compose this transform with `other` (self applied after other):
    /// `result = self ∘ other`.
    pub fn compose(&self, other: &Transform3D) -> Self {
        let new_rot = quat_mul(self.rotation, other.rotation);
        let scaled_other_pos = [
            other.position[0] * self.scale,
            other.position[1] * self.scale,
            other.position[2] * self.scale,
        ];
        let rotated = quat_rotate(self.rotation, scaled_other_pos);
        let new_pos = [
            self.position[0] + rotated[0],
            self.position[1] + rotated[1],
            self.position[2] + rotated[2],
        ];
        Self {
            position: new_pos,
            rotation: new_rot,
            scale: self.scale * other.scale,
        }
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// Global physics configuration.
#[derive(Debug, Clone)]
pub struct PhysicsConfig {
    /// Gravity vector (default: -9.81 on Y).
    pub gravity: Vec3,
    /// Number of constraint solver iterations.
    pub solver_iterations: u32,
    /// Linear velocity threshold below which bodies may sleep.
    pub linear_sleep_threshold: Real,
    /// Angular velocity threshold below which bodies may sleep.
    pub angular_sleep_threshold: Real,
    /// Time a body must be below thresholds before sleeping.
    pub time_before_sleep: Real,
    /// Enable continuous collision detection (CCD) to prevent tunneling.
    pub ccd_enabled: bool,
}
impl PhysicsConfig {
    /// Validate the configuration. Returns an error message if invalid.
    pub fn validate(&self) -> Result<(), String> {
        if self.solver_iterations == 0 {
            return Err("solver_iterations must be > 0".to_string());
        }
        if self.linear_sleep_threshold < 0.0 {
            return Err("linear_sleep_threshold must be >= 0".to_string());
        }
        if self.angular_sleep_threshold < 0.0 {
            return Err("angular_sleep_threshold must be >= 0".to_string());
        }
        if self.time_before_sleep < 0.0 {
            return Err("time_before_sleep must be >= 0".to_string());
        }
        Ok(())
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// Pre-built configurations for common physics scenarios.
pub struct DefaultConfigs;
impl DefaultConfigs {
    /// Standard Earth-gravity rigid body config.
    pub fn rigid_body() -> PhysicsConfig {
        PhysicsConfig::default()
    }
    /// Soft-body simulation config (more iterations, lower sleep thresholds).
    pub fn soft_body() -> PhysicsConfig {
        PhysicsConfig {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 20,
            linear_sleep_threshold: 0.001,
            angular_sleep_threshold: 0.001,
            time_before_sleep: 1.0,
            ccd_enabled: false,
        }
    }
    /// Fluid-like simulation (zero gravity, many iterations).
    pub fn fluid() -> PhysicsConfig {
        PhysicsConfig {
            gravity: Vec3::zeros(),
            solver_iterations: 50,
            linear_sleep_threshold: 0.0,
            angular_sleep_threshold: 0.0,
            time_before_sleep: f64::MAX,
            ccd_enabled: false,
        }
    }
    /// Space simulation (no gravity, CCD enabled).
    pub fn space() -> PhysicsConfig {
        PhysicsConfig {
            gravity: Vec3::zeros(),
            solver_iterations: 8,
            linear_sleep_threshold: 0.001,
            angular_sleep_threshold: 0.001,
            time_before_sleep: 2.0,
            ccd_enabled: true,
        }
    }
    /// Moon gravity simulation.
    pub fn moon() -> PhysicsConfig {
        PhysicsConfig {
            gravity: Vec3::new(0.0, -1.625, 0.0),
            ..PhysicsConfig::default()
        }
    }
    /// Mars gravity simulation.
    pub fn mars() -> PhysicsConfig {
        PhysicsConfig {
            gravity: Vec3::new(0.0, -3.72, 0.0),
            ..PhysicsConfig::default()
        }
    }
}
/// Handle to a collider, with generation counter for safe reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ColliderHandle {
    /// Index into the collider storage.
    pub index: u32,
    /// Generation counter for detecting stale handles.
    pub generation: u32,
}
impl ColliderHandle {
    /// Create a new collider handle.
    pub fn new(index: u32, generation: u32) -> Self {
        Self { index, generation }
    }
}
impl ColliderHandle {
    /// Check if this is a valid handle.
    pub fn is_valid(&self) -> bool {
        self.generation > 0 || self.index > 0
    }
    /// Create a null handle.
    pub fn null() -> Self {
        Self {
            index: 0,
            generation: 0,
        }
    }
}
/// A 3-D axis-aligned bounding box with richer query helpers.
///
/// Internally identical to `Aabb` but offers additional bulk/intersection APIs.
#[derive(Debug, Clone)]
pub struct BoundingBox3D {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl BoundingBox3D {
    /// Construct from explicit min/max corner arrays.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Merge a list of `BoundingBox3D`s into the smallest enclosing box.
    ///
    /// Returns `None` when the list is empty.
    pub fn merge_many(boxes: &[BoundingBox3D]) -> Option<Self> {
        if boxes.is_empty() {
            return None;
        }
        let mut min = boxes[0].min;
        let mut max = boxes[0].max;
        for b in boxes.iter().skip(1) {
            for i in 0..3 {
                if b.min[i] < min[i] {
                    min[i] = b.min[i];
                }
                if b.max[i] > max[i] {
                    max[i] = b.max[i];
                }
            }
        }
        Some(Self { min, max })
    }
    /// Compute the intersection (overlap) of two bounding boxes.
    ///
    /// Returns `None` if the boxes do not overlap.
    pub fn intersection(&self, other: &BoundingBox3D) -> Option<BoundingBox3D> {
        let mut min = [0.0f64; 3];
        let mut max = [0.0f64; 3];
        for i in 0..3 {
            min[i] = self.min[i].max(other.min[i]);
            max[i] = self.max[i].min(other.max[i]);
            if min[i] > max[i] {
                return None;
            }
        }
        Some(BoundingBox3D { min, max })
    }
    /// Volume of the bounding box.
    pub fn volume(&self) -> f64 {
        (0..3)
            .map(|i| (self.max[i] - self.min[i]).max(0.0))
            .product()
    }
    /// Whether a point is inside (inclusive) the box.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] <= self.max[i])
    }
}
/// A discrete simulation time step.
#[derive(Debug, Clone, Copy)]
pub struct TimeStep {
    /// Duration of this time step in seconds.
    pub dt: Real,
}
impl TimeStep {
    /// Create a new time step with the given duration.
    pub fn new(dt: Real) -> Self {
        Self { dt }
    }
}
impl TimeStep {
    /// Validate the time step. Returns error if dt <= 0.
    pub fn validate(&self) -> Result<(), String> {
        if self.dt <= 0.0 {
            return Err("dt must be positive".to_string());
        }
        Ok(())
    }
    /// Frequency corresponding to this time step: 1/dt.
    pub fn frequency(&self) -> f64 {
        1.0 / self.dt
    }
}
/// A ray in 3-D space stored with plain `[f64; 3]` arrays (no nalgebra).
#[derive(Debug, Clone)]
pub struct Ray3D {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Unit direction.
    pub direction: [f64; 3],
}
impl Ray3D {
    /// Create a new ray. `direction` is normalised internally.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        let dir = if len > 1e-12 {
            [direction[0] / len, direction[1] / len, direction[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            origin,
            direction: dir,
        }
    }
    /// Return the point along the ray at parameter `t`: `origin + t * direction`.
    pub fn parameter_at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}
/// An infinite plane stored with plain `[f64; 3]` arrays (no nalgebra).
///
/// Equation: `dot(normal, p) = offset`.
#[derive(Debug, Clone)]
pub struct Plane3D {
    /// Unit normal.
    pub normal: [f64; 3],
    /// Plane offset along the normal from the origin.
    pub offset: f64,
}
impl Plane3D {
    /// Create a plane from a unit normal and offset.  The normal is normalised.
    pub fn new(normal: [f64; 3], offset: f64) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let n = if len > 1e-12 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self { normal: n, offset }
    }
    /// Signed distance from `point` to the plane.
    pub fn signed_distance(&self, point: [f64; 3]) -> f64 {
        self.normal[0] * point[0] + self.normal[1] * point[1] + self.normal[2] * point[2]
            - self.offset
    }
    /// Orthogonal projection of `point` onto the plane.
    pub fn project_point(&self, point: [f64; 3]) -> [f64; 3] {
        let d = self.signed_distance(point);
        [
            point[0] - d * self.normal[0],
            point[1] - d * self.normal[1],
            point[2] - d * self.normal[2],
        ]
    }
}
