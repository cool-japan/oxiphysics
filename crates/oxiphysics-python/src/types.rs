// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python-friendly wrapper types with serde support.
//!
//! All types in this module are designed for easy FFI/Python interop:
//! - No lifetimes
//! - Public fields only
//! - Serialize/Deserialize support
//! - Arrays instead of nalgebra types where applicable

use oxiphysics_core::math::Vec3;
use oxiphysics_core::{Aabb, Transform};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PyVec3
// ---------------------------------------------------------------------------

/// A 3-component vector suitable for Python interop.
///
/// Uses `f64` components and maps directly to/from `oxiphysics_core::math::Vec3`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct PyVec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl PyVec3 {
    /// Create a new `PyVec3`.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Dot product with another vector.
    pub fn dot(&self, other: &PyVec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product with another vector.
    pub fn cross(&self, other: &PyVec3) -> PyVec3 {
        PyVec3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Squared length.
    pub fn length_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Euclidean length.
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    /// Normalized (unit) vector; returns zero vector if length is near zero.
    pub fn normalized(&self) -> PyVec3 {
        let len = self.length();
        if len < 1e-12 {
            PyVec3::zero()
        } else {
            PyVec3::new(self.x / len, self.y / len, self.z / len)
        }
    }

    /// Scale all components by `s`.
    pub fn scale(&self, s: f64) -> PyVec3 {
        PyVec3::new(self.x * s, self.y * s, self.z * s)
    }

    /// Add another vector.
    pub fn add(&self, other: &PyVec3) -> PyVec3 {
        PyVec3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    /// Subtract another vector.
    pub fn sub(&self, other: &PyVec3) -> PyVec3 {
        PyVec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    /// Convert to a plain `[f64; 3]` array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Create from a `[f64; 3]` array.
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }
}

impl From<Vec3> for PyVec3 {
    fn from(v: Vec3) -> Self {
        Self {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

impl From<PyVec3> for Vec3 {
    fn from(v: PyVec3) -> Self {
        Vec3::new(v.x, v.y, v.z)
    }
}

impl From<[f64; 3]> for PyVec3 {
    fn from(arr: [f64; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }
}

impl From<PyVec3> for [f64; 3] {
    fn from(v: PyVec3) -> Self {
        [v.x, v.y, v.z]
    }
}

// ---------------------------------------------------------------------------
// PyTransform
// ---------------------------------------------------------------------------

/// A transform (position + quaternion rotation) suitable for Python interop.
///
/// Rotation is stored as `[x, y, z, w]` quaternion components.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTransform {
    /// Position in world space.
    pub position: PyVec3,
    /// Rotation as a quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
}

impl PyTransform {
    /// Create an identity transform (origin, no rotation).
    pub fn identity() -> Self {
        Self {
            position: PyVec3::zero(),
            rotation: [0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Compose two transforms: `self` applied first, then `other`.
    ///
    /// Position is summed; rotation is the quaternion product `other * self`.
    pub fn compose(&self, other: &PyTransform) -> PyTransform {
        let (ax, ay, az, aw) = (
            self.rotation[0],
            self.rotation[1],
            self.rotation[2],
            self.rotation[3],
        );
        let (bx, by, bz, bw) = (
            other.rotation[0],
            other.rotation[1],
            other.rotation[2],
            other.rotation[3],
        );
        let rx = bw * ax + bx * aw + by * az - bz * ay;
        let ry = bw * ay - bx * az + by * aw + bz * ax;
        let rz = bw * az + bx * ay - by * ax + bz * aw;
        let rw = bw * aw - bx * ax - by * ay - bz * az;
        PyTransform {
            position: PyVec3::new(
                self.position.x + other.position.x,
                self.position.y + other.position.y,
                self.position.z + other.position.z,
            ),
            rotation: [rx, ry, rz, rw],
        }
    }

    /// Invert the transform.
    pub fn inverse(&self) -> PyTransform {
        // Invert quaternion: conjugate (negated xyz, same w) / norm^2
        // For unit quaternion conjugate == inverse
        let inv_rot = [
            -self.rotation[0],
            -self.rotation[1],
            -self.rotation[2],
            self.rotation[3],
        ];
        // Inverse position: -inv_rot * position
        let px = self.position.x;
        let py = self.position.y;
        let pz = self.position.z;
        let (qx, qy, qz, qw) = (inv_rot[0], inv_rot[1], inv_rot[2], inv_rot[3]);
        // Rotate position by inv_rot
        let t_x = 2.0 * (qy * pz - qz * py);
        let t_y = 2.0 * (qz * px - qx * pz);
        let t_z = 2.0 * (qx * py - qy * px);
        let rx = px + qw * t_x + (qy * t_z - qz * t_y);
        let ry = py + qw * t_y + (qz * t_x - qx * t_z);
        let rz = pz + qw * t_z + (qx * t_y - qy * t_x);
        PyTransform {
            position: PyVec3::new(-rx, -ry, -rz),
            rotation: inv_rot,
        }
    }
}

impl Default for PyTransform {
    fn default() -> Self {
        Self::identity()
    }
}

impl From<Transform> for PyTransform {
    fn from(t: Transform) -> Self {
        let q = t.rotation;
        Self {
            position: PyVec3::from(t.position),
            rotation: [q.i, q.j, q.k, q.w],
        }
    }
}

impl From<PyTransform> for Transform {
    fn from(t: PyTransform) -> Self {
        let pos: Vec3 = t.position.into();
        let q = nalgebra::UnitQuaternion::from_quaternion(nalgebra::Quaternion::new(
            t.rotation[3],
            t.rotation[0],
            t.rotation[1],
            t.rotation[2],
        ));
        Transform::new(pos, q)
    }
}

// ---------------------------------------------------------------------------
// PyAabb
// ---------------------------------------------------------------------------

/// An axis-aligned bounding box suitable for Python interop.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyAabb {
    /// Minimum corner.
    pub min: PyVec3,
    /// Maximum corner.
    pub max: PyVec3,
}

impl PyAabb {
    /// Create a new AABB from min and max corners.
    pub fn new(min: PyVec3, max: PyVec3) -> Self {
        Self { min, max }
    }

    /// Center of the AABB.
    pub fn center(&self) -> PyVec3 {
        PyVec3::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
            (self.min.z + self.max.z) * 0.5,
        )
    }

    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> PyVec3 {
        PyVec3::new(
            (self.max.x - self.min.x) * 0.5,
            (self.max.y - self.min.y) * 0.5,
            (self.max.z - self.min.z) * 0.5,
        )
    }

    /// Whether the given point is inside this AABB.
    pub fn contains_point(&self, p: PyVec3) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// Whether this AABB intersects another.
    pub fn intersects(&self, other: &PyAabb) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }
}

impl From<Aabb> for PyAabb {
    fn from(a: Aabb) -> Self {
        Self {
            min: PyVec3::from(a.min),
            max: PyVec3::from(a.max),
        }
    }
}

impl From<PyAabb> for Aabb {
    fn from(a: PyAabb) -> Self {
        Aabb::new(a.min.into(), a.max.into())
    }
}

// ---------------------------------------------------------------------------
// PyColliderShape
// ---------------------------------------------------------------------------

/// Shape variant for a rigid body collider.
///
/// Only primitive shapes are supported; they map to underlying collision shapes.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum PyColliderShape {
    /// A sphere with the given radius.
    Sphere {
        /// Radius in meters.
        radius: f64,
    },
    /// An axis-aligned box with given half-extents.
    Box {
        /// Half-extents `[hx, hy, hz]` in meters.
        half_extents: [f64; 3],
    },
    /// An infinite plane defined by a normal and distance from origin.
    Plane {
        /// Outward normal (should be unit length).
        normal: [f64; 3],
        /// Signed distance from origin along the normal.
        distance: f64,
    },
    /// A capsule (cylinder + two hemispheres) along the Y axis.
    Capsule {
        /// Radius of the cylinder and end hemispheres.
        radius: f64,
        /// Half-height of the cylindrical portion (total height = 2 * half_height + 2 * radius).
        half_height: f64,
    },
    /// A cylinder aligned with the Y axis.
    Cylinder {
        /// Radius of the cylinder.
        radius: f64,
        /// Half-height of the cylinder.
        half_height: f64,
    },
    /// A cone aligned with the Y axis, apex pointing up.
    Cone {
        /// Radius of the base.
        radius: f64,
        /// Half-height of the cone.
        half_height: f64,
    },
}

impl PyColliderShape {
    /// Create a sphere shape.
    pub fn sphere(radius: f64) -> Self {
        PyColliderShape::Sphere { radius }
    }

    /// Create a box shape from half-extents.
    pub fn box_shape(hx: f64, hy: f64, hz: f64) -> Self {
        PyColliderShape::Box {
            half_extents: [hx, hy, hz],
        }
    }

    /// Create a plane shape.
    pub fn plane(normal: [f64; 3], distance: f64) -> Self {
        PyColliderShape::Plane { normal, distance }
    }

    /// Create a capsule shape.
    pub fn capsule(radius: f64, half_height: f64) -> Self {
        PyColliderShape::Capsule {
            radius,
            half_height,
        }
    }

    /// Compute an approximate volume for mass/inertia estimation.
    pub fn approximate_volume(&self) -> f64 {
        const PI: f64 = std::f64::consts::PI;
        match self {
            PyColliderShape::Sphere { radius } => (4.0 / 3.0) * PI * radius * radius * radius,
            PyColliderShape::Box { half_extents } => {
                8.0 * half_extents[0] * half_extents[1] * half_extents[2]
            }
            PyColliderShape::Plane { .. } => f64::INFINITY,
            PyColliderShape::Capsule {
                radius,
                half_height,
            } => PI * radius * radius * (2.0 * half_height + (4.0 / 3.0) * radius),
            PyColliderShape::Cylinder {
                radius,
                half_height,
            } => PI * radius * radius * 2.0 * half_height,
            PyColliderShape::Cone {
                radius,
                half_height,
            } => (1.0 / 3.0) * PI * radius * radius * 2.0 * half_height,
        }
    }

    /// Whether this shape is a plane (infinite extents).
    pub fn is_infinite(&self) -> bool {
        matches!(self, PyColliderShape::Plane { .. })
    }
}

// ---------------------------------------------------------------------------
// PyRigidBodyConfig
// ---------------------------------------------------------------------------

/// Configuration for creating a new rigid body.
///
/// Passed to `PyPhysicsWorld::add_rigid_body`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyRigidBodyConfig {
    /// Mass in kilograms. Use 0.0 or `f64::INFINITY` for a static/kinematic body.
    pub mass: f64,
    /// Initial position `[x, y, z]`.
    pub position: [f64; 3],
    /// Initial linear velocity `[vx, vy, vz]`. Defaults to zero.
    pub velocity: [f64; 3],
    /// Initial orientation as quaternion `[x, y, z, w]`. Defaults to identity.
    pub orientation: [f64; 4],
    /// Initial angular velocity `[wx, wy, wz]`. Defaults to zero.
    pub angular_velocity: [f64; 3],
    /// Collider shapes attached to this body.
    pub shapes: Vec<PyColliderShape>,
    /// Friction coefficient (0 = frictionless, 1 = high friction).
    pub friction: f64,
    /// Restitution/bounciness (0 = perfectly inelastic, 1 = perfectly elastic).
    pub restitution: f64,
    /// Whether this body is static (immovable, infinite mass).
    pub is_static: bool,
    /// Whether this body is kinematic (moved manually, not by forces).
    pub is_kinematic: bool,
    /// Whether this body can go to sleep when motion is negligible.
    pub can_sleep: bool,
    /// Linear damping coefficient (drag).
    pub linear_damping: f64,
    /// Angular damping coefficient (rotational drag).
    pub angular_damping: f64,
    /// Optional user-defined tag/name for identification.
    pub tag: Option<String>,
}

impl PyRigidBodyConfig {
    /// Create a dynamic body at the given position with default settings.
    pub fn dynamic(mass: f64, position: [f64; 3]) -> Self {
        Self {
            mass,
            position,
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            shapes: vec![],
            friction: 0.5,
            restitution: 0.3,
            is_static: false,
            is_kinematic: false,
            can_sleep: true,
            linear_damping: 0.0,
            angular_damping: 0.0,
            tag: None,
        }
    }

    /// Create a static body at the given position with default settings.
    pub fn static_body(position: [f64; 3]) -> Self {
        Self {
            mass: 0.0,
            position,
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            shapes: vec![],
            friction: 0.5,
            restitution: 0.3,
            is_static: true,
            is_kinematic: false,
            can_sleep: false,
            linear_damping: 0.0,
            angular_damping: 0.0,
            tag: None,
        }
    }

    /// Add a collider shape to this body config and return `self` for chaining.
    pub fn with_shape(mut self, shape: PyColliderShape) -> Self {
        self.shapes.push(shape);
        self
    }

    /// Set friction and return `self` for chaining.
    pub fn with_friction(mut self, friction: f64) -> Self {
        self.friction = friction;
        self
    }

    /// Set restitution and return `self` for chaining.
    pub fn with_restitution(mut self, restitution: f64) -> Self {
        self.restitution = restitution;
        self
    }

    /// Set linear damping and return `self` for chaining.
    pub fn with_linear_damping(mut self, damping: f64) -> Self {
        self.linear_damping = damping;
        self
    }

    /// Set angular damping and return `self` for chaining.
    pub fn with_angular_damping(mut self, damping: f64) -> Self {
        self.angular_damping = damping;
        self
    }

    /// Set tag and return `self` for chaining.
    pub fn with_tag(mut self, tag: impl Into<String>) -> Self {
        self.tag = Some(tag.into());
        self
    }

    /// Effective inverse mass (0 for static/infinite).
    pub fn inverse_mass(&self) -> f64 {
        if self.is_static || self.mass <= 0.0 {
            0.0
        } else {
            1.0 / self.mass
        }
    }
}

impl Default for PyRigidBodyConfig {
    fn default() -> Self {
        Self::dynamic(1.0, [0.0; 3])
    }
}

// ---------------------------------------------------------------------------
// PyContactResult
// ---------------------------------------------------------------------------

/// Result of a contact/collision between two bodies.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyContactResult {
    /// Handle of the first body involved.
    pub body_a: u32,
    /// Handle of the second body involved.
    pub body_b: u32,
    /// Contact point in world space.
    pub contact_point: [f64; 3],
    /// Contact normal (from body_a to body_b).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Combined friction coefficient for this contact.
    pub friction: f64,
    /// Combined restitution coefficient for this contact.
    pub restitution: f64,
    /// Impulse applied to resolve this contact.
    pub impulse: f64,
}

impl PyContactResult {
    /// Create a new contact result.
    pub fn new(
        body_a: u32,
        body_b: u32,
        contact_point: [f64; 3],
        normal: [f64; 3],
        depth: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            contact_point,
            normal,
            depth,
            friction: 0.5,
            restitution: 0.3,
            impulse: 0.0,
        }
    }

    /// Whether this contact represents a genuine collision (positive depth).
    pub fn is_colliding(&self) -> bool {
        self.depth > 0.0
    }

    /// Relative velocity of body_a with respect to body_b projected onto the normal.
    /// A negative value means separating; positive means approaching.
    pub fn separating_velocity(&self, vel_a: [f64; 3], vel_b: [f64; 3]) -> f64 {
        let rel = [
            vel_a[0] - vel_b[0],
            vel_a[1] - vel_b[1],
            vel_a[2] - vel_b[2],
        ];
        rel[0] * self.normal[0] + rel[1] * self.normal[1] + rel[2] * self.normal[2]
    }
}

// ---------------------------------------------------------------------------
// PySimConfig
// ---------------------------------------------------------------------------

/// Top-level simulation configuration.
///
/// Passed to `PyPhysicsWorld::new` to configure global simulation parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySimConfig {
    /// Gravity vector `[gx, gy, gz]`.
    pub gravity: [f64; 3],
    /// Number of constraint solver iterations per step.
    pub solver_iterations: u32,
    /// Linear velocity threshold for sleep.
    pub linear_sleep_threshold: f64,
    /// Angular velocity threshold for sleep.
    pub angular_sleep_threshold: f64,
    /// Time (in seconds) a body must be below thresholds before sleeping.
    pub time_before_sleep: f64,
    /// Enable continuous collision detection.
    pub ccd_enabled: bool,
    /// Coefficient of restitution mixing method: "average", "min", or "max".
    pub restitution_mixing: String,
    /// Friction mixing method: "average", "min", or "max".
    pub friction_mixing: String,
    /// Maximum allowed penetration depth before a contact is ignored.
    pub max_penetration: f64,
    /// Baumgarte stabilization factor (0..1, typical 0.1–0.3).
    pub baumgarte_factor: f64,
    /// Whether sleeping is enabled globally.
    pub sleep_enabled: bool,
    /// Maximum substeps for CCD.
    pub ccd_max_substeps: u32,
}

impl PySimConfig {
    /// Standard Earth-gravity configuration.
    pub fn earth_gravity() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            ..Self::default()
        }
    }

    /// Zero-gravity (space) configuration.
    pub fn zero_gravity() -> Self {
        Self {
            gravity: [0.0, 0.0, 0.0],
            ..Self::default()
        }
    }

    /// Moon-gravity configuration (~1/6 of Earth).
    pub fn moon_gravity() -> Self {
        Self {
            gravity: [0.0, -1.62, 0.0],
            ..Self::default()
        }
    }
}

impl Default for PySimConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            solver_iterations: 8,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            ccd_enabled: true,
            restitution_mixing: "average".to_string(),
            friction_mixing: "average".to_string(),
            max_penetration: 0.001,
            baumgarte_factor: 0.2,
            sleep_enabled: true,
            ccd_max_substeps: 10,
        }
    }
}

// ---------------------------------------------------------------------------
// Legacy types kept for backward compatibility
// ---------------------------------------------------------------------------

/// Description for creating a rigid body from Python (legacy).
///
/// Prefer `PyRigidBodyConfig` for new code.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyRigidBodyDesc {
    /// Mass in kilograms.
    pub mass: f64,
    /// Initial position.
    pub position: PyVec3,
    /// Whether this body is static (immovable).
    pub is_static: bool,
}

/// Description for creating a collider from Python (legacy).
///
/// Prefer `PyColliderShape` for new code.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyColliderDesc {
    /// Shape type name (e.g. "sphere", "box").
    pub shape_type: String,
    /// Half-extents for box shapes.
    pub half_extents: Option<PyVec3>,
    /// Radius for sphere shapes.
    pub radius: Option<f64>,
    /// Friction coefficient.
    pub friction: f64,
    /// Restitution (bounciness) coefficient.
    pub restitution: f64,
}

impl PyColliderDesc {
    /// Create a sphere collider description with default material properties.
    pub fn sphere(radius: f64) -> Self {
        Self {
            shape_type: "sphere".to_string(),
            half_extents: None,
            radius: Some(radius),
            friction: 0.5,
            restitution: 0.3,
        }
    }

    /// Create a box collider description with default material properties.
    pub fn box_shape(half_extents: PyVec3) -> Self {
        Self {
            shape_type: "box".to_string(),
            half_extents: Some(half_extents),
            radius: None,
            friction: 0.5,
            restitution: 0.3,
        }
    }
}

/// Material properties for Python interop.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMaterial {
    /// Material name.
    pub name: String,
    /// Density in kg/m^3.
    pub density: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// Restitution coefficient.
    pub restitution: f64,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pyvec3_roundtrip() {
        let original = PyVec3::new(1.0, 2.0, 3.0);
        let v3: Vec3 = original.into();
        let back = PyVec3::from(v3);
        assert_eq!(original, back);
    }

    #[test]
    fn test_pytransform_conversion() {
        let t = Transform::default();
        let py_t = PyTransform::from(t);
        assert!((py_t.position.x).abs() < 1e-10);
        assert!((py_t.rotation[3] - 1.0).abs() < 1e-10); // w=1 for identity

        let back: Transform = py_t.into();
        assert!(back.position.norm() < 1e-10);
    }

    #[test]
    fn test_py_vec3_arithmetic() {
        let a = PyVec3::new(1.0, 2.0, 3.0);
        let b = PyVec3::new(4.0, 5.0, 6.0);

        // add
        let sum = a.add(&b);
        assert!((sum.x - 5.0).abs() < 1e-10);
        assert!((sum.y - 7.0).abs() < 1e-10);
        assert!((sum.z - 9.0).abs() < 1e-10);

        // sub
        let diff = b.sub(&a);
        assert!((diff.x - 3.0).abs() < 1e-10);
        assert!((diff.y - 3.0).abs() < 1e-10);
        assert!((diff.z - 3.0).abs() < 1e-10);

        // scale
        let scaled = a.scale(2.0);
        assert!((scaled.x - 2.0).abs() < 1e-10);
        assert!((scaled.y - 4.0).abs() < 1e-10);
        assert!((scaled.z - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_py_transform_compose() {
        let t1 = PyTransform {
            position: PyVec3::new(1.0, 0.0, 0.0),
            rotation: [0.0, 0.0, 0.0, 1.0],
        };
        let t2 = PyTransform {
            position: PyVec3::new(2.0, 0.0, 0.0),
            rotation: [0.0, 0.0, 0.0, 1.0],
        };
        let composed = t1.compose(&t2);
        assert!((composed.position.x - 3.0).abs() < 1e-10);
        assert!((composed.position.y).abs() < 1e-10);
        assert!((composed.position.z).abs() < 1e-10);
        // identity * identity = identity quaternion
        assert!((composed.rotation[3] - 1.0).abs() < 1e-10);
        assert!((composed.rotation[0]).abs() < 1e-10);
    }

    #[test]
    fn test_py_collider_desc_sphere() {
        let desc = PyColliderDesc::sphere(0.5);
        assert_eq!(desc.shape_type, "sphere");
        assert!((desc.radius.unwrap() - 0.5).abs() < 1e-10);
        assert!(desc.half_extents.is_none());

        let json = serde_json::to_string(&desc).expect("serialize");
        assert!(json.contains("0.5"));
        assert!(json.contains("sphere"));
    }

    #[test]
    fn test_py_json_roundtrip() {
        let original = PyVec3::new(1.5, -2.3, 0.7);
        let json = serde_json::to_string(&original).expect("serialize PyVec3");
        let restored: PyVec3 = serde_json::from_str(&json).expect("deserialize PyVec3");
        assert!((restored.x - original.x).abs() < 1e-10);
        assert!((restored.y - original.y).abs() < 1e-10);
        assert!((restored.z - original.z).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_sphere_volume() {
        let s = PyColliderShape::sphere(1.0);
        let vol = s.approximate_volume();
        let expected = (4.0 / 3.0) * std::f64::consts::PI;
        assert!((vol - expected).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_box_volume() {
        let b = PyColliderShape::box_shape(1.0, 2.0, 3.0);
        let vol = b.approximate_volume();
        // 8 * 1 * 2 * 3 = 48
        assert!((vol - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_plane_is_infinite() {
        let p = PyColliderShape::plane([0.0, 1.0, 0.0], 0.0);
        assert!(p.is_infinite());
    }

    #[test]
    fn test_rigid_body_config_inverse_mass() {
        let cfg = PyRigidBodyConfig::dynamic(2.0, [0.0; 3]);
        assert!((cfg.inverse_mass() - 0.5).abs() < 1e-10);

        let static_cfg = PyRigidBodyConfig::static_body([0.0; 3]);
        assert!((static_cfg.inverse_mass()).abs() < 1e-10);
    }

    #[test]
    fn test_rigid_body_config_builder() {
        let cfg = PyRigidBodyConfig::dynamic(5.0, [1.0, 2.0, 3.0])
            .with_friction(0.8)
            .with_restitution(0.1)
            .with_tag("test_body")
            .with_shape(PyColliderShape::sphere(0.5));

        assert!((cfg.friction - 0.8).abs() < 1e-10);
        assert!((cfg.restitution - 0.1).abs() < 1e-10);
        assert_eq!(cfg.tag.as_deref(), Some("test_body"));
        assert_eq!(cfg.shapes.len(), 1);
    }

    #[test]
    fn test_sim_config_default() {
        let cfg = PySimConfig::default();
        assert!((cfg.gravity[1] + 9.81).abs() < 1e-10);
        assert_eq!(cfg.solver_iterations, 8);
    }

    #[test]
    fn test_sim_config_moon_gravity() {
        let cfg = PySimConfig::moon_gravity();
        assert!((cfg.gravity[1] + 1.62).abs() < 1e-10);
    }

    #[test]
    fn test_contact_result_is_colliding() {
        let contact = PyContactResult::new(0, 1, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.05);
        assert!(contact.is_colliding());

        let no_contact = PyContactResult::new(0, 1, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], -0.01);
        assert!(!no_contact.is_colliding());
    }

    #[test]
    fn test_pyvec3_dot_cross() {
        let a = PyVec3::new(1.0, 0.0, 0.0);
        let b = PyVec3::new(0.0, 1.0, 0.0);
        assert!((a.dot(&b)).abs() < 1e-10);

        let c = a.cross(&b);
        // (1,0,0) x (0,1,0) = (0,0,1)
        assert!(c.x.abs() < 1e-10);
        assert!(c.y.abs() < 1e-10);
        assert!((c.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pyvec3_normalize() {
        let v = PyVec3::new(3.0, 0.0, 4.0);
        let n = v.normalized();
        assert!((n.length() - 1.0).abs() < 1e-10);

        let zero = PyVec3::zero();
        let nz = zero.normalized();
        assert!(nz.length() < 1e-10);
    }

    #[test]
    fn test_pyaabb_contains() {
        let aabb = PyAabb::new(PyVec3::new(-1.0, -1.0, -1.0), PyVec3::new(1.0, 1.0, 1.0));
        assert!(aabb.contains_point(PyVec3::new(0.0, 0.0, 0.0)));
        assert!(!aabb.contains_point(PyVec3::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn test_pyaabb_intersects() {
        let a = PyAabb::new(PyVec3::new(0.0, 0.0, 0.0), PyVec3::new(2.0, 2.0, 2.0));
        let b = PyAabb::new(PyVec3::new(1.0, 1.0, 1.0), PyVec3::new(3.0, 3.0, 3.0));
        let c = PyAabb::new(PyVec3::new(5.0, 5.0, 5.0), PyVec3::new(6.0, 6.0, 6.0));
        assert!(a.intersects(&b));
        assert!(!a.intersects(&c));
    }

    #[test]
    fn test_rigid_body_config_json_roundtrip() {
        let cfg = PyRigidBodyConfig::dynamic(3.0, [1.0, 2.0, 3.0])
            .with_shape(PyColliderShape::sphere(0.5))
            .with_friction(0.7);
        let json = serde_json::to_string(&cfg).expect("serialize");
        let restored: PyRigidBodyConfig = serde_json::from_str(&json).expect("deserialize");
        assert!((restored.mass - 3.0).abs() < 1e-10);
        assert_eq!(restored.shapes.len(), 1);
        assert!((restored.friction - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_contact_separating_velocity() {
        let contact = PyContactResult::new(0, 1, [0.0; 3], [0.0, 1.0, 0.0], 0.01);
        // Bodies moving apart: vel_a = (0,1,0), vel_b = (0,-1,0)
        let sv = contact.separating_velocity([0.0, 1.0, 0.0], [0.0, -1.0, 0.0]);
        // relative = (0,2,0), dot with normal (0,1,0) = 2
        assert!((sv - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_shape_capsule_volume() {
        let c = PyColliderShape::capsule(1.0, 2.0);
        let vol = c.approximate_volume();
        let expected = std::f64::consts::PI * 1.0 * 1.0 * (4.0 + (4.0 / 3.0) * 1.0);
        assert!((vol - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sim_config_zero_gravity() {
        let cfg = PySimConfig::zero_gravity();
        for &g in &cfg.gravity {
            assert!(g.abs() < 1e-10);
        }
    }

    #[test]
    fn test_collider_shape_serde_roundtrip() {
        let shape = PyColliderShape::box_shape(0.5, 1.0, 1.5);
        let json = serde_json::to_string(&shape).expect("serialize shape");
        let back: PyColliderShape = serde_json::from_str(&json).expect("deserialize shape");
        assert_eq!(shape, back);
    }
}
