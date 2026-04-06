// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core data types: transforms, bounding volumes, handles, configs.

pub use crate::math::{Mat3, Quat, Real, Vec3};

/// Transform representing position and orientation in 3D space.
#[derive(Debug, Clone)]
pub struct Transform {
    /// Position in world space.
    pub position: Vec3,
    /// Orientation as a unit quaternion.
    pub rotation: Quat,
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            position: Vec3::zeros(),
            rotation: Quat::identity(),
        }
    }
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

/// Mass properties of a rigid body.
#[derive(Debug, Clone)]
pub struct MassProperties {
    /// Mass in kilograms.
    pub mass: Real,
    /// Center of mass in local space.
    pub center_of_mass: Vec3,
    /// Local inertia tensor (3x3 matrix).
    pub local_inertia: Mat3,
}

impl MassProperties {
    /// Create new mass properties.
    pub fn new(mass: Real, center_of_mass: Vec3, local_inertia: Mat3) -> Self {
        Self {
            mass,
            center_of_mass,
            local_inertia,
        }
    }

    /// Mass properties for a point mass.
    pub fn point_mass(mass: Real) -> Self {
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::zeros(),
        }
    }

    /// Inverse mass (0 for infinite/static mass).
    pub fn inverse_mass(&self) -> Real {
        if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        }
    }

    /// Inverse inertia tensor (zero matrix for infinite/static inertia).
    pub fn inverse_inertia(&self) -> Mat3 {
        if self.mass > 0.0 {
            self.local_inertia.try_inverse().unwrap_or_else(Mat3::zeros)
        } else {
            Mat3::zeros()
        }
    }
}

/// Handle to a rigid body, with generation counter for safe reuse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyHandle {
    /// Index into the body storage.
    pub index: u32,
    /// Generation counter for detecting stale handles.
    pub generation: u32,
}

impl BodyHandle {
    /// Create a new body handle.
    pub fn new(index: u32, generation: u32) -> Self {
        Self { index, generation }
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

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 8,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            ccd_enabled: true,
        }
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

#[cfg(test)]
mod prop_tests {

    use crate::Aabb;

    use crate::Quat;

    use crate::Transform;

    use crate::Vec3;
    use proptest::prelude::*;

    fn coord_strategy() -> impl Strategy<Value = f64> {
        -100.0_f64..100.0_f64
    }

    fn vec3_strategy() -> impl Strategy<Value = Vec3> {
        (coord_strategy(), coord_strategy(), coord_strategy())
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }

    fn positive_coord_strategy() -> impl Strategy<Value = f64> {
        0.01_f64..100.0_f64
    }

    fn positive_vec3_strategy() -> impl Strategy<Value = Vec3> {
        (
            positive_coord_strategy(),
            positive_coord_strategy(),
            positive_coord_strategy(),
        )
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }

    fn transform_strategy() -> impl Strategy<Value = Transform> {
        (
            vec3_strategy(),
            coord_strategy(),
            coord_strategy(),
            coord_strategy(),
            coord_strategy(),
        )
            .prop_map(|(pos, qi, qj, qk, qw)| {
                // Build a unit quaternion from raw values by normalizing
                let raw = nalgebra::Quaternion::new(qw, qi, qj, qk);
                let norm = raw.norm();
                let quat = if norm < 1e-10 {
                    Quat::identity()
                } else {
                    Quat::from_quaternion(nalgebra::Quaternion::new(
                        raw.w / norm,
                        raw.i / norm,
                        raw.j / norm,
                        raw.k / norm,
                    ))
                };
                Transform::new(pos, quat)
            })
    }

    proptest! {
        #[test]
        fn prop_transform_compose_associative(
            ta in transform_strategy(),
            tb in transform_strategy(),
            tc in transform_strategy(),
        ) {
            let lhs = ta.compose(&tb).compose(&tc);
            let rhs = ta.compose(&tb.compose(&tc));
            let pos_diff = (lhs.position - rhs.position).norm();
            prop_assert!(pos_diff < 1e-6, "positions differ by {}", pos_diff);
        }

        #[test]
        fn prop_transform_inverse(ta in transform_strategy()) {
            let inv = ta.inverse();
            let composed = ta.compose(&inv);
            let pos_diff = composed.position.norm();
            prop_assert!(pos_diff < 1e-6, "T * T^-1 position not zero: norm={}", pos_diff);
            let rot_diff = (composed.rotation * Quat::identity().inverse()).angle();
            prop_assert!(rot_diff < 1e-6, "T * T^-1 rotation not identity: angle={}", rot_diff);
        }

        #[test]
        fn prop_aabb_merge_contains_both(
            min_a in vec3_strategy(),
            ext_a in positive_vec3_strategy(),
            min_b in vec3_strategy(),
            ext_b in positive_vec3_strategy(),
        ) {
            let a = Aabb::new(min_a, min_a + ext_a);
            let b = Aabb::new(min_b, min_b + ext_b);
            let merged = a.merge(&b);
            prop_assert!(merged.contains_point(&a.min), "merged does not contain a.min");
            prop_assert!(merged.contains_point(&a.max), "merged does not contain a.max");
            prop_assert!(merged.contains_point(&b.min), "merged does not contain b.min");
            prop_assert!(merged.contains_point(&b.max), "merged does not contain b.max");
        }

        #[test]
        fn prop_aabb_surface_area_positive(
            min in vec3_strategy(),
            ext in positive_vec3_strategy(),
        ) {
            let aabb = Aabb::new(min, min + ext);
            let sa = aabb.surface_area();
            prop_assert!(sa > 0.0, "surface area not positive: {}", sa);
        }

        #[test]
        fn prop_vec3_dot_commutative(
            a in vec3_strategy(),
            b in vec3_strategy(),
        ) {
            let ab = a.dot(&b);
            let ba = b.dot(&a);
            prop_assert!((ab - ba).abs() < 1e-6, "dot not commutative: {} vs {}", ab, ba);
        }

        #[test]
        fn prop_vec3_cross_anti_commutative(
            a in vec3_strategy(),
            b in vec3_strategy(),
        ) {
            let ab = a.cross(&b);
            let ba = b.cross(&a);
            let diff = (ab + ba).norm();
            prop_assert!(diff < 1e-6, "cross not anti-commutative: diff norm={}", diff);
        }
    }
}

// ─── Type Conversions ────────────────────────────────────────────────────────

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

// ─── AABB Conversion ─────────────────────────────────────────────────────────

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

// ─── Type Validation ─────────────────────────────────────────────────────────

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

// ─── Builder Patterns ────────────────────────────────────────────────────────

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

impl Default for PhysicsConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Default Configs for Physics Domains ─────────────────────────────────────

/// Pre-built configurations for common physics scenarios.
pub struct DefaultConfigs;

#[allow(dead_code)]
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

// ─── MassProperties Extensions ──────────────────────────────────────────────

impl MassProperties {
    /// Mass properties for a uniform-density sphere.
    pub fn sphere(mass: Real, radius: Real) -> Self {
        let i_scalar = (2.0 / 5.0) * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(i_scalar, 0.0, 0.0, 0.0, i_scalar, 0.0, 0.0, 0.0, i_scalar),
        }
    }

    /// Mass properties for a uniform-density box with dimensions (wx, wy, wz).
    pub fn cuboid(mass: Real, wx: Real, wy: Real, wz: Real) -> Self {
        let c = mass / 12.0;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(
                c * (wy * wy + wz * wz),
                0.0,
                0.0,
                0.0,
                c * (wx * wx + wz * wz),
                0.0,
                0.0,
                0.0,
                c * (wx * wx + wy * wy),
            ),
        }
    }

    /// Mass properties for a uniform-density cylinder (axis along Y).
    pub fn cylinder(mass: Real, radius: Real, height: Real) -> Self {
        let ix = mass * (3.0 * radius * radius + height * height) / 12.0;
        let iy = 0.5 * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }

    /// Check if this represents a static (infinite mass) body.
    pub fn is_static(&self) -> bool {
        self.mass <= 0.0
    }

    /// Combine two mass properties (additive).
    pub fn combine(&self, other: &MassProperties) -> MassProperties {
        MassProperties {
            mass: self.mass + other.mass,
            center_of_mass: (self.center_of_mass * self.mass + other.center_of_mass * other.mass)
                / (self.mass + other.mass),
            local_inertia: self.local_inertia + other.local_inertia,
        }
    }
}

// ─── Additional Handle Operations ───────────────────────────────────────────

impl BodyHandle {
    /// Check if this is a valid (non-null) handle.
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

#[cfg(test)]
mod tests {

    use crate::Aabb;
    use crate::BodyHandle;

    use crate::ColliderHandle;
    use crate::DefaultConfigs;
    use crate::MassProperties;
    use crate::PhysicsConfig;
    use crate::PhysicsConfigBuilder;

    use crate::Quat;

    use crate::TimeStep;
    use crate::Transform;

    use crate::Vec3;

    #[test]
    fn test_aabb_intersection() {
        let a = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(0.5, 0.5, 0.5), Vec3::new(1.5, 1.5, 1.5));
        let c = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));

        assert!(a.intersects(&b));
        assert!(!a.intersects(&c));
    }

    #[test]
    fn test_aabb_merge() {
        let a = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        let merged = a.merge(&b);
        assert_eq!(merged.min, Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(merged.max, Vec3::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_aabb_contains_point() {
        let aabb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        assert!(aabb.contains_point(&Vec3::new(0.5, 0.5, 0.5)));
        assert!(!aabb.contains_point(&Vec3::new(2.0, 0.5, 0.5)));
    }

    #[test]
    fn test_aabb_surface_area() {
        let aabb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 3.0));
        let sa = aabb.surface_area();
        // 2*(1*2 + 2*3 + 3*1) = 2*(2+6+3) = 22
        assert!((sa - 22.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_compose() {
        // T1 translates by (1,0,0), T2 translates by (0,2,0)
        let t1 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.0, 2.0, 0.0));
        let composed = t1.compose(&t2);
        // Applying composed to origin should give (1,2,0)
        let p = composed.transform_point(&Vec3::zeros());
        assert!((p.x - 1.0).abs() < 1e-10);
        assert!((p.y - 2.0).abs() < 1e-10);
        assert!(p.z.abs() < 1e-10);
    }

    #[test]
    fn test_transform_inverse() {
        let t = Transform::from_position(Vec3::new(1.0, 2.0, 3.0));
        let inv = t.inverse();
        let composed = t.compose(&inv);
        assert!(composed.position.norm() < 1e-10);
        // rotation should be identity: angle of (composed.rotation * identity^-1) ~ 0
        let angle = (composed.rotation * Quat::identity().inverse()).angle();
        assert!(angle < 1e-10);
    }

    #[test]
    fn test_transform_point() {
        let t = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let p = t.transform_point(&Vec3::new(1.0, 0.0, 0.0));
        assert!((p.x - 11.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_vector() {
        // With identity rotation, transform_vector should not add translation
        let t = Transform::from_position(Vec3::new(10.0, 5.0, 3.0));
        let v = Vec3::new(1.0, 0.0, 0.0);
        let tv = t.transform_vector(&v);
        // translation must not affect vectors
        assert!((tv.x - 1.0).abs() < 1e-10);
        assert!(tv.y.abs() < 1e-10);
        assert!(tv.z.abs() < 1e-10);
    }

    #[test]
    fn test_transform_default() {
        let t = Transform::default();
        assert_eq!(t.position, Vec3::zeros());
    }

    #[test]
    fn test_mass_properties() {
        let mp = MassProperties::point_mass(2.0);
        assert!((mp.inverse_mass() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_body_handle() {
        let h = BodyHandle::new(0, 0);
        assert_eq!(h.index, 0);
        assert_eq!(h.generation, 0);
    }

    #[test]
    fn test_physics_config_default() {
        let cfg = PhysicsConfig::default();
        assert!((cfg.gravity.y + 9.81).abs() < 1e-10);
    }

    // ─── New tests for expanded types ───────────────────────────────────

    #[test]
    fn test_transform_to_matrix4_identity() {
        let t = Transform::default();
        let m = t.to_matrix4();
        // Diagonal should be 1.0
        assert!((m[0][0] - 1.0).abs() < 1e-10);
        assert!((m[1][1] - 1.0).abs() < 1e-10);
        assert!((m[2][2] - 1.0).abs() < 1e-10);
        assert!((m[3][3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_euler_angles_identity() {
        let t = Transform::default();
        let (roll, pitch, yaw) = t.euler_angles();
        assert!(roll.abs() < 1e-10);
        assert!(pitch.abs() < 1e-10);
        assert!(yaw.abs() < 1e-10);
    }

    #[test]
    fn test_transform_from_axis_angle() {
        let t = Transform::from_axis_angle(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::PI / 2.0,
        );
        // After 90° rotation about Z, [1,0,0] → [0,1,0]
        let v = t.transform_vector(&Vec3::new(1.0, 0.0, 0.0));
        assert!((v.y - 1.0).abs() < 1e-10, "v = {:?}", v);
    }

    #[test]
    fn test_transform_lerp() {
        let t1 = Transform::from_position(Vec3::zeros());
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let mid = t1.lerp(&t2, 0.5);
        assert!((mid.position.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_from_center_half_extents() {
        let aabb =
            Aabb::from_center_half_extents(Vec3::new(5.0, 5.0, 5.0), Vec3::new(1.0, 2.0, 3.0));
        assert!((aabb.min.x - 4.0).abs() < 1e-10);
        assert!((aabb.max.x - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_from_points() {
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, -2.0, -3.0),
        ];
        let aabb = Aabb::from_points(&points).unwrap();
        assert!((aabb.min.x + 1.0).abs() < 1e-10);
        assert!((aabb.max.z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_from_points_empty() {
        let result = Aabb::from_points(&[]);
        assert!(result.is_none());
    }

    #[test]
    fn test_aabb_longest_axis() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 3.0, 2.0));
        assert_eq!(aabb.longest_axis(), 1); // y is longest
    }

    #[test]
    fn test_aabb_diagonal() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(3.0, 4.0, 0.0));
        assert!((aabb.diagonal() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_physics_config_validate_valid() {
        let cfg = PhysicsConfig::default();
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_physics_config_validate_invalid() {
        let cfg = PhysicsConfig {
            solver_iterations: 0,
            ..PhysicsConfig::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_timestep_validate() {
        let ts = TimeStep::new(0.01);
        assert!(ts.validate().is_ok());
        let ts_bad = TimeStep::new(-0.01);
        assert!(ts_bad.validate().is_err());
    }

    #[test]
    fn test_timestep_frequency() {
        let ts = TimeStep::new(0.01);
        assert!((ts.frequency() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_builder_pattern() {
        let cfg = PhysicsConfigBuilder::new()
            .gravity(Vec3::new(0.0, -1.625, 0.0))
            .solver_iterations(16)
            .ccd_enabled(false)
            .build();
        assert!((cfg.gravity.y + 1.625).abs() < 1e-10);
        assert_eq!(cfg.solver_iterations, 16);
        assert!(!cfg.ccd_enabled);
    }

    #[test]
    fn test_builder_default() {
        let builder = PhysicsConfigBuilder::default();
        let cfg = builder.build();
        assert!((cfg.gravity.y + 9.81).abs() < 1e-10);
    }

    #[test]
    fn test_default_configs_rigid_body() {
        let cfg = DefaultConfigs::rigid_body();
        assert!(cfg.ccd_enabled);
    }

    #[test]
    fn test_default_configs_soft_body() {
        let cfg = DefaultConfigs::soft_body();
        assert_eq!(cfg.solver_iterations, 20);
    }

    #[test]
    fn test_default_configs_fluid() {
        let cfg = DefaultConfigs::fluid();
        assert_eq!(cfg.solver_iterations, 50);
        assert!(cfg.gravity.norm() < 1e-10);
    }

    #[test]
    fn test_default_configs_space() {
        let cfg = DefaultConfigs::space();
        assert!(cfg.gravity.norm() < 1e-10);
        assert!(cfg.ccd_enabled);
    }

    #[test]
    fn test_default_configs_moon() {
        let cfg = DefaultConfigs::moon();
        assert!((cfg.gravity.y + 1.625).abs() < 1e-10);
    }

    #[test]
    fn test_default_configs_mars() {
        let cfg = DefaultConfigs::mars();
        assert!((cfg.gravity.y + 3.72).abs() < 1e-10);
    }

    #[test]
    fn test_mass_properties_sphere() {
        let mp = MassProperties::sphere(10.0, 1.0);
        assert!((mp.mass - 10.0).abs() < 1e-10);
        // I = 2/5 * m * r^2 = 4.0
        let i = mp.local_inertia[(0, 0)];
        assert!((i - 4.0).abs() < 1e-10, "I_sphere = {i}");
    }

    #[test]
    fn test_mass_properties_cuboid() {
        let mp = MassProperties::cuboid(12.0, 1.0, 2.0, 3.0);
        assert!((mp.mass - 12.0).abs() < 1e-10);
        // Ixx = m/12 * (wy^2 + wz^2) = 12/12 * (4+9) = 13
        let ixx = mp.local_inertia[(0, 0)];
        assert!((ixx - 13.0).abs() < 1e-10, "Ixx = {ixx}");
    }

    #[test]
    fn test_mass_properties_cylinder() {
        let mp = MassProperties::cylinder(6.0, 1.0, 2.0);
        assert!((mp.mass - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_mass_properties_is_static() {
        let mp = MassProperties::point_mass(0.0);
        assert!(mp.is_static());
        let mp2 = MassProperties::point_mass(1.0);
        assert!(!mp2.is_static());
    }

    #[test]
    fn test_mass_properties_combine() {
        let mp1 = MassProperties::point_mass(2.0);
        let mp2 = MassProperties::point_mass(3.0);
        let combined = mp1.combine(&mp2);
        assert!((combined.mass - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_body_handle_null() {
        let h = BodyHandle::null();
        assert!(!h.is_valid());
        let h2 = BodyHandle::new(1, 0);
        assert!(h2.is_valid());
    }

    #[test]
    fn test_collider_handle_null() {
        let h = ColliderHandle::null();
        assert!(!h.is_valid());
    }
}

// ─── Ray / Segment types ─────────────────────────────────────────────────────

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
            .unwrap_or_else(Vec3::z_axis_unit);
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
            .unwrap_or_else(Vec3::z_axis_unit);
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

// ─── Ray–AABB intersection ────────────────────────────────────────────────────

/// Result of a ray–AABB intersection test.
#[derive(Debug, Clone, Copy)]
pub struct RayAabbHit {
    /// Entry parameter (t_min).
    pub t_min: Real,
    /// Exit parameter (t_max).
    pub t_max: Real,
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
                // Ray is parallel to slab — check if origin is inside
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
        // Project each column of the rotation matrix onto the world axes.
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

// ─── Sphere primitive ─────────────────────────────────────────────────────────

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
        // Initial approximate sphere
        let mut center = points[0];
        let mut radius = 0.0f64;

        for &p in points.iter().skip(1) {
            let d = (p - center).norm();
            if d > radius {
                // Expand sphere to include p
                radius = (radius + d) * 0.5;
                center = center + (p - center) * ((d - radius) / d + 0.5);
            }
        }

        // Second pass: make sure all points are inside
        for &p in points {
            let d = (p - center).norm();
            if d > radius {
                radius = d;
            }
        }

        Some(Sphere::new(center, radius))
    }
}

// ─── Capsule primitive ────────────────────────────────────────────────────────

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

// ─── Plane primitive ──────────────────────────────────────────────────────────

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

// ─── Triangle primitive ───────────────────────────────────────────────────────

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
            return None; // Ray is parallel to triangle
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

// ─── Vec3 extension helpers ───────────────────────────────────────────────────

/// Extension trait providing a unit-Z-axis fallback (used internally).
trait Vec3ZAxisUnit {
    fn z_axis_unit() -> Vec3;
}

impl Vec3ZAxisUnit for Vec3 {
    fn z_axis_unit() -> Vec3 {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

// ─── Additional MassProperties shapes ────────────────────────────────────────

impl MassProperties {
    /// Mass properties for a uniform-density cone (apex at top, base at bottom,
    /// axis along +Y). Height = `height`, base radius = `radius`.
    pub fn cone(mass: Real, radius: Real, height: Real) -> Self {
        let ix = mass * (3.0 * radius * radius / 20.0 + height * height * 3.0 / 80.0);
        let iy = 3.0 * mass * radius * radius / 10.0;
        Self {
            mass,
            // Center of mass is at 1/4 height from base.
            center_of_mass: Vec3::new(0.0, height * 0.25, 0.0),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }

    /// Mass properties for a hollow spherical shell.
    pub fn spherical_shell(mass: Real, radius: Real) -> Self {
        let i_scalar = (2.0 / 3.0) * mass * radius * radius;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(i_scalar, 0.0, 0.0, 0.0, i_scalar, 0.0, 0.0, 0.0, i_scalar),
        }
    }

    /// Mass properties for a uniform-density capsule (cylinder + 2 hemispheres,
    /// axis along Y).
    pub fn capsule(mass: Real, radius: Real, half_height: Real) -> Self {
        // Decompose into a cylinder and two hemispheres.
        let h = 2.0 * half_height;
        let vol_cyl = std::f64::consts::PI * radius * radius * h;
        let vol_hemi = (2.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        let vol_total = vol_cyl + 2.0 * vol_hemi;
        let density = mass / vol_total;
        let m_cyl = density * vol_cyl;
        let m_hemi = density * vol_hemi;

        // Cylinder inertia
        let ix_cyl = m_cyl * (3.0 * radius * radius + h * h) / 12.0;
        let iy_cyl = 0.5 * m_cyl * radius * radius;

        // Hemisphere inertia about axis through own centroid, then parallel-axis to capsule center
        let ix_hemi = 83.0 * m_hemi * radius * radius / 320.0;
        let iy_hemi = 2.0 * m_hemi * radius * radius / 5.0;
        // Distance from hemisphere centroid to capsule axis center
        let d = half_height + 3.0 * radius / 8.0;
        let ix_hemi_total = ix_hemi + m_hemi * d * d;

        let ix = ix_cyl + 2.0 * ix_hemi_total;
        let iy = iy_cyl + 2.0 * iy_hemi;
        Self {
            mass,
            center_of_mass: Vec3::zeros(),
            local_inertia: Mat3::new(ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, ix),
        }
    }

    /// Parallel-axis theorem: shift the inertia tensor by offset `d` from the
    /// center of mass.
    pub fn shifted_inertia(&self, d: Vec3) -> Mat3 {
        let m = self.mass;
        let d2 = d.dot(&d);
        // I_new = I_cm + m * (|d|² E - d ⊗ d)
        let outer = Mat3::new(
            d.x * d.x,
            d.x * d.y,
            d.x * d.z,
            d.y * d.x,
            d.y * d.y,
            d.y * d.z,
            d.z * d.x,
            d.z * d.y,
            d.z * d.z,
        );
        self.local_inertia + m * (Mat3::identity() * d2 - outer)
    }
}

// ─── BoundingBox3D (alias-style extensions on Aabb) ───────────────────────────

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
    #[allow(dead_code)]
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Merge a list of `BoundingBox3D`s into the smallest enclosing box.
    ///
    /// Returns `None` when the list is empty.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn volume(&self) -> f64 {
        (0..3)
            .map(|i| (self.max[i] - self.min[i]).max(0.0))
            .product()
    }

    /// Whether a point is inside (inclusive) the box.
    #[allow(dead_code)]
    pub fn contains(&self, p: [f64; 3]) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] <= self.max[i])
    }
}

// ─── Ray3D ────────────────────────────────────────────────────────────────────

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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn parameter_at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}

// ─── Plane3D ──────────────────────────────────────────────────────────────────

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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn signed_distance(&self, point: [f64; 3]) -> f64 {
        self.normal[0] * point[0] + self.normal[1] * point[1] + self.normal[2] * point[2]
            - self.offset
    }

    /// Orthogonal projection of `point` onto the plane.
    #[allow(dead_code)]
    pub fn project_point(&self, point: [f64; 3]) -> [f64; 3] {
        let d = self.signed_distance(point);
        [
            point[0] - d * self.normal[0],
            point[1] - d * self.normal[1],
            point[2] - d * self.normal[2],
        ]
    }
}

// ─── Transform3D ─────────────────────────────────────────────────────────────

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

impl Default for Transform3D {
    fn default() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0], // identity
            scale: 1.0,
        }
    }
}

/// Multiply two quaternions `a * b` (Hamilton product), both stored as `[x,y,z,w]`.
#[allow(dead_code)]
fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let (ax, ay, az, aw) = (a[0], a[1], a[2], a[3]);
    let (bx, by, bz, bw) = (b[0], b[1], b[2], b[3]);
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}

/// Rotate a vector `v` by quaternion `q` (`[x,y,z,w]`) using `q * v * q^-1`.
#[allow(dead_code)]
fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    // Pure quaternion representation of v
    let pv: [f64; 4] = [v[0], v[1], v[2], 0.0];
    // Conjugate of q
    let qc: [f64; 4] = [-q[0], -q[1], -q[2], q[3]];
    let tmp = quat_mul(q, pv);
    let res = quat_mul(tmp, qc);
    [res[0], res[1], res[2]]
}

impl Transform3D {
    /// Create an identity transform (no translation, no rotation, scale = 1).
    #[allow(dead_code)]
    pub fn identity() -> Self {
        Self::default()
    }

    /// Create a transform with position only (identity rotation, scale = 1).
    #[allow(dead_code)]
    pub fn from_position(position: [f64; 3]) -> Self {
        Self {
            position,
            ..Self::default()
        }
    }

    /// Apply the transform to a point: `scale * rotate(point) + position`.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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

// ─── Additional tests ─────────────────────────────────────────────────────────

#[cfg(test)]
mod types_expanded_tests {

    use crate::Aabb;

    use crate::BoundingBox3D;
    use crate::Capsule;

    use crate::MassProperties;

    use crate::Plane;
    use crate::Plane3D;

    use crate::Ray;
    use crate::Ray3D;
    use crate::Sphere;

    use crate::Transform;
    use crate::Transform3D;
    use crate::Triangle;
    use crate::Vec3;
    use std::f64::consts::PI;

    // ── Ray ───────────────────────────────────────────────────────────────────

    #[test]
    fn test_ray_at() {
        let ray = Ray::new(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let p = ray.at(3.0);
        assert!((p.x - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_from_to() {
        let ray = Ray::from_to(Vec3::zeros(), Vec3::new(5.0, 0.0, 0.0));
        assert!((ray.max_t - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_normalizes_direction() {
        let ray = Ray::new(Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0));
        assert!((ray.direction.norm() - 1.0).abs() < 1e-12);
    }

    // ── AABB ray intersection ─────────────────────────────────────────────────

    #[test]
    fn test_aabb_ray_hit() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let hit = aabb.ray_intersect(&ray);
        assert!(hit.is_some(), "ray should hit unit cube");
        let h = hit.unwrap();
        assert!((h.t_min - 4.0).abs() < 1e-10, "t_min={}", h.t_min);
    }

    #[test]
    fn test_aabb_ray_miss() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(-5.0, 5.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(aabb.ray_intersect(&ray).is_none());
    }

    #[test]
    fn test_aabb_ray_inside() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let hit = aabb.ray_intersect(&ray);
        assert!(hit.is_some());
    }

    #[test]
    fn test_aabb_intersection_overlapping() {
        let a = Aabb::new(Vec3::zeros(), Vec3::new(2.0, 2.0, 2.0));
        let b = Aabb::new(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
        let inter = a.intersection(&b).unwrap();
        assert!((inter.min.x - 1.0).abs() < 1e-10);
        assert!((inter.max.x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_intersection_disjoint() {
        let a = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        assert!(a.intersection(&b).is_none());
    }

    #[test]
    fn test_aabb_sq_distance_outside() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(2.0, 0.5, 0.5); // 1 unit outside on X
        assert!((aabb.sq_distance_to_point(&p) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_sq_distance_inside_is_zero() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(0.5, 0.5, 0.5);
        assert_eq!(aabb.sq_distance_to_point(&p), 0.0);
    }

    #[test]
    fn test_aabb_closest_point() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(2.0, 0.5, 0.5);
        let cp = aabb.closest_point(&p);
        assert!((cp.x - 1.0).abs() < 1e-10);
        assert!((cp.y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_aabb_oriented_box() {
        let t = Transform::default();
        let he = Vec3::new(1.0, 2.0, 3.0);
        let aabb = Aabb::oriented_box_aabb(he, &t);
        assert!((aabb.max.x - 1.0).abs() < 1e-10);
        assert!((aabb.max.y - 2.0).abs() < 1e-10);
        assert!((aabb.max.z - 3.0).abs() < 1e-10);
    }

    // ── Sphere ────────────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_contains_point() {
        let s = Sphere::new(Vec3::zeros(), 2.0);
        assert!(s.contains_point(&Vec3::new(1.0, 0.0, 0.0)));
        assert!(!s.contains_point(&Vec3::new(3.0, 0.0, 0.0)));
    }

    #[test]
    fn test_sphere_overlaps_sphere() {
        let a = Sphere::new(Vec3::zeros(), 1.0);
        let b = Sphere::new(Vec3::new(1.5, 0.0, 0.0), 1.0);
        let c = Sphere::new(Vec3::new(5.0, 0.0, 0.0), 1.0);
        assert!(a.overlaps_sphere(&b));
        assert!(!a.overlaps_sphere(&c));
    }

    #[test]
    fn test_sphere_volume() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let expected = 4.0 / 3.0 * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_ray_intersect_hit() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let ray = Ray::new(Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let t = s.ray_intersect(&ray);
        assert!(t.is_some(), "ray should hit sphere");
        assert!((t.unwrap() - 4.0).abs() < 1e-10, "t={}", t.unwrap());
    }

    #[test]
    fn test_sphere_ray_intersect_miss() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let ray = Ray::new(Vec3::new(-5.0, 2.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(s.ray_intersect(&ray).is_none());
    }

    #[test]
    fn test_sphere_bounding_sphere() {
        let pts = vec![Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let s = Sphere::bounding_sphere(&pts).unwrap();
        assert!(s.radius >= 1.0 - 1e-9);
    }

    #[test]
    fn test_sphere_bounding_sphere_single_point() {
        let pts = vec![Vec3::new(3.0, 4.0, 0.0)];
        let s = Sphere::bounding_sphere(&pts).unwrap();
        assert_eq!(s.radius, 0.0);
    }

    // ── Capsule ───────────────────────────────────────────────────────────────

    #[test]
    fn test_capsule_contains_point_on_axis() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        assert!(c.contains_point(&Vec3::new(0.0, 0.0, 0.0)));
    }

    #[test]
    fn test_capsule_volume() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        let v = c.volume();
        assert!(v > 0.0, "volume should be positive: {v}");
    }

    #[test]
    fn test_capsule_aabb_upright() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        // height = 2; radius = 0.5; total Y-extent = 2.0 + 1.0 (caps) = 3.0
        let aabb = c.aabb();
        assert!(aabb.max.y >= 1.0, "y max should be >= 1.0: {}", aabb.max.y);
    }

    // ── Plane ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_plane_signed_distance() {
        let p = Plane::from_point_normal(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        let above = Vec3::new(0.0, 5.0, 0.0);
        let below = Vec3::new(0.0, -1.0, 0.0);
        assert!(p.signed_distance(&above) > 0.0);
        assert!(p.signed_distance(&below) < 0.0);
    }

    #[test]
    fn test_plane_project_point() {
        let p = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let projected = p.project_point(&Vec3::new(3.0, 5.0, 2.0));
        assert!(
            projected.y.abs() < 1e-10,
            "projected y should be 0: {}",
            projected.y
        );
    }

    #[test]
    fn test_plane_ray_intersect() {
        let plane = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let ray = Ray::new(Vec3::new(0.0, 5.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let t = plane.ray_intersect(&ray);
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-10, "t={}", t.unwrap());
    }

    #[test]
    fn test_plane_ray_parallel_miss() {
        let plane = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let ray = Ray::new(Vec3::new(0.0, 5.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(plane.ray_intersect(&ray).is_none());
    }

    // ── Triangle ──────────────────────────────────────────────────────────────

    #[test]
    fn test_triangle_area() {
        // Right triangle with legs 3 and 4 → area = 6
        let t = Triangle::new(
            Vec3::zeros(),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.0, 4.0, 0.0),
        );
        assert!((t.area() - 6.0).abs() < 1e-10, "area={}", t.area());
    }

    #[test]
    fn test_triangle_normal() {
        let t = Triangle::new(
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let n = t.normal();
        // Normal should be (0, 0, 1)
        assert!((n.z - 1.0).abs() < 1e-10, "n={:?}", n);
    }

    #[test]
    fn test_triangle_ray_hit() {
        let tri = Triangle::new(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let ray = Ray::new(Vec3::new(0.0, 0.3, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let t = tri.ray_intersect(&ray);
        assert!(t.is_some(), "ray should hit triangle");
    }

    #[test]
    fn test_triangle_ray_miss() {
        let tri = Triangle::new(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let ray = Ray::new(Vec3::new(5.0, 5.0, -1.0), Vec3::new(0.0, 0.0, 1.0));
        assert!(tri.ray_intersect(&ray).is_none());
    }

    #[test]
    fn test_triangle_barycentric_centroid() {
        let tri = Triangle::new(
            Vec3::zeros(),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.0, 3.0, 0.0),
        );
        let centroid = tri.centroid();
        let (u, v, w) = tri.barycentric(&centroid);
        assert!((u - 1.0 / 3.0).abs() < 1e-10, "u={}", u);
        assert!((v - 1.0 / 3.0).abs() < 1e-10, "v={}", v);
        assert!((w - 1.0 / 3.0).abs() < 1e-10, "w={}", w);
    }

    // ── MassProperties shapes ─────────────────────────────────────────────────

    #[test]
    fn test_mass_properties_cone() {
        let mp = MassProperties::cone(3.0, 1.0, 2.0);
        assert!((mp.mass - 3.0).abs() < 1e-10);
        assert!(mp.local_inertia[(0, 0)] > 0.0);
    }

    #[test]
    fn test_mass_properties_spherical_shell() {
        let mp = MassProperties::spherical_shell(1.0, 1.0);
        // I = 2/3 * 1 * 1² = 0.666...
        assert!((mp.local_inertia[(0, 0)] - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_mass_properties_capsule() {
        let mp = MassProperties::capsule(5.0, 0.5, 1.0);
        assert!((mp.mass - 5.0).abs() < 1e-10);
        assert!(mp.local_inertia[(0, 0)] > 0.0);
    }

    #[test]
    fn test_mass_properties_shifted_inertia() {
        let mp = MassProperties::sphere(1.0, 1.0);
        // Shift by [1, 0, 0]: Ixx_new = Ixx + m*(0) = 0.4; Iyy_new = Iyy + m*1 = 1.4
        let shifted = mp.shifted_inertia(Vec3::new(1.0, 0.0, 0.0));
        assert!(
            (shifted[(0, 0)] - 0.4).abs() < 1e-10,
            "Ixx={}",
            shifted[(0, 0)]
        );
        assert!(
            (shifted[(1, 1)] - 1.4).abs() < 1e-10,
            "Iyy={}",
            shifted[(1, 1)]
        );
    }

    // ── BoundingBox3D ─────────────────────────────────────────────────────────

    #[test]
    fn test_bounding_box3d_merge_many_empty() {
        let result = BoundingBox3D::merge_many(&[]);
        assert!(result.is_none(), "merge of empty list should be None");
    }

    #[test]
    fn test_bounding_box3d_merge_many_single() {
        let b = BoundingBox3D::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let merged = BoundingBox3D::merge_many(&[b]).unwrap();
        assert_eq!(merged.min, [1.0, 2.0, 3.0]);
        assert_eq!(merged.max, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_bounding_box3d_merge_many_two() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = BoundingBox3D::new([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0]);
        let merged = BoundingBox3D::merge_many(&[a, b]).unwrap();
        assert_eq!(merged.min, [-1.0, -1.0, -1.0]);
        assert_eq!(merged.max, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_bounding_box3d_intersection_overlap() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = BoundingBox3D::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let inter = a.intersection(&b).unwrap();
        assert_eq!(inter.min, [1.0, 1.0, 1.0]);
        assert_eq!(inter.max, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_bounding_box3d_intersection_no_overlap() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = BoundingBox3D::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(
            a.intersection(&b).is_none(),
            "non-overlapping boxes should return None"
        );
    }

    // ── Ray3D ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_ray3d_parameter_at_origin() {
        let ray = Ray3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = ray.parameter_at(0.0);
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_ray3d_parameter_at_positive_t() {
        let ray = Ray3D::new([1.0, 2.0, 3.0], [0.0, 0.0, 1.0]);
        let p = ray.parameter_at(5.0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 8.0).abs() < 1e-10);
    }

    // ── Plane3D ───────────────────────────────────────────────────────────────

    #[test]
    fn test_plane3d_project_point_on_plane() {
        // XY plane at z=0: project (1,2,0) → should stay (1,2,0)
        let plane = Plane3D::new([0.0, 0.0, 1.0], 0.0);
        let proj = plane.project_point([1.0, 2.0, 0.0]);
        assert!((proj[0] - 1.0).abs() < 1e-10);
        assert!((proj[1] - 2.0).abs() < 1e-10);
        assert!(proj[2].abs() < 1e-10);
    }

    #[test]
    fn test_plane3d_project_point_above_plane() {
        // XY plane at z=0: project (1,2,3) → (1,2,0)
        let plane = Plane3D::new([0.0, 0.0, 1.0], 0.0);
        let proj = plane.project_point([1.0, 2.0, 3.0]);
        assert!((proj[0] - 1.0).abs() < 1e-10);
        assert!((proj[1] - 2.0).abs() < 1e-10);
        assert!(proj[2].abs() < 1e-10);
    }

    // ── Transform3D ───────────────────────────────────────────────────────────

    #[test]
    fn test_transform3d_identity_apply() {
        let t = Transform3D::identity();
        let p = t.apply([3.0, 4.0, 5.0]);
        assert!((p[0] - 3.0).abs() < 1e-10);
        assert!((p[1] - 4.0).abs() < 1e-10);
        assert!((p[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform3d_translation_apply() {
        let mut t = Transform3D::identity();
        t.position = [1.0, 2.0, 3.0];
        let p = t.apply([0.0, 0.0, 0.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform3d_inverse_cancels() {
        let mut t = Transform3D::identity();
        t.position = [1.0, 2.0, 3.0];
        let inv = t.inverse();
        let composed = t.compose(&inv);
        // Compose should give identity-like transform
        assert!(
            composed.position[0].abs() < 1e-9,
            "x={}",
            composed.position[0]
        );
        assert!(
            composed.position[1].abs() < 1e-9,
            "y={}",
            composed.position[1]
        );
        assert!(
            composed.position[2].abs() < 1e-9,
            "z={}",
            composed.position[2]
        );
    }
}
