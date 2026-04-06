// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WASM-friendly types using flat f64 arrays for JavaScript interop.
//!
//! All types in this module are designed to cross the WASM boundary efficiently.
//! They avoid heap-allocated types where possible and use flat arrays instead
//! of nalgebra types, making them trivially serializable and transferable.
//!
//! ## Design Principles
//!
//! - No nalgebra dependencies: all math uses raw `[f64; N]` arrays.
//! - All fields are `pub` for direct JavaScript access via wasm-bindgen.
//! - All types derive `Clone`, `Debug`, and `serde` traits.
//! - Comprehensive doc comments on every type and field.

#![allow(missing_docs)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Vec3
// ---------------------------------------------------------------------------

/// A 3D vector represented as a flat `[x, y, z]` array for WASM interop.
///
/// This type avoids any nalgebra dependency and is directly usable at the
/// WASM boundary. All arithmetic methods are provided as plain functions on
/// the `[f64; 3]` values.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::Vec3Wasm;
///
/// let v = Vec3Wasm::new(1.0, 2.0, 3.0);
/// assert!((v.length() - f64::sqrt(14.0)).abs() < 1e-10);
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Vec3Wasm {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3Wasm {
    /// Create a new `Vec3Wasm` from components.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Create the zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the X axis.
    pub fn unit_x() -> Self {
        Self {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the Y axis.
    pub fn unit_y() -> Self {
        Self {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    }

    /// Create a unit vector along the Z axis.
    pub fn unit_z() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    }

    /// Return the inner `[x, y, z]` array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Create from a `[f64; 3]` array.
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }

    /// Compute the squared length (L2 norm squared).
    pub fn length_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Compute the length (L2 norm).
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }

    /// Return a normalized copy, or the zero vector if length is zero.
    pub fn normalized(&self) -> Self {
        let len = self.length();
        if len < 1e-15 {
            Self::zero()
        } else {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
            }
        }
    }

    /// Dot product with another vector.
    pub fn dot(&self, other: &Vec3Wasm) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product with another vector.
    pub fn cross(&self, other: &Vec3Wasm) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Component-wise addition.
    pub fn add(&self, other: &Vec3Wasm) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    /// Component-wise subtraction.
    pub fn sub(&self, other: &Vec3Wasm) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    /// Linear interpolation towards `other` by factor `t` (clamped to \[0, 1\]).
    pub fn lerp(&self, other: &Vec3Wasm, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        self.add(&other.sub(self).scale(t))
    }

    /// Distance to another vector.
    pub fn distance_to(&self, other: &Vec3Wasm) -> f64 {
        self.sub(other).length()
    }
}

impl Default for Vec3Wasm {
    fn default() -> Self {
        Self::zero()
    }
}

// ---------------------------------------------------------------------------
// Quaternion (rotation)
// ---------------------------------------------------------------------------

/// A unit quaternion for rotations, stored as `[x, y, z, w]`.
///
/// Follows the convention where `w` is the scalar part and `(x, y, z)` is
/// the vector part. The identity quaternion is `[0, 0, 0, 1]`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct QuatWasm {
    /// X (imaginary i) component.
    pub x: f64,
    /// Y (imaginary j) component.
    pub y: f64,
    /// Z (imaginary k) component.
    pub z: f64,
    /// W (real) component.
    pub w: f64,
}

impl QuatWasm {
    /// Create a quaternion from components `(x, y, z, w)`.
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    /// The identity quaternion `(0, 0, 0, 1)`.
    pub fn identity() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    /// Return the inner `[x, y, z, w]` array.
    pub fn to_array(&self) -> [f64; 4] {
        [self.x, self.y, self.z, self.w]
    }

    /// Create from a `[f64; 4]` array `[x, y, z, w]`.
    pub fn from_array(arr: [f64; 4]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
            w: arr[3],
        }
    }

    /// Normalize the quaternion (so it has unit length).
    pub fn normalized(&self) -> Self {
        let len = (self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w).sqrt();
        if len < 1e-15 {
            Self::identity()
        } else {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
                w: self.w / len,
            }
        }
    }

    /// Conjugate of this quaternion (`-x, -y, -z, w`).
    pub fn conjugate(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    /// Quaternion multiplication (Hamilton product).
    pub fn mul(&self, other: &QuatWasm) -> Self {
        Self {
            x: self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            y: self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            z: self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
            w: self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        }
    }

    /// Rotate a 3D vector by this quaternion.
    pub fn rotate_vec(&self, v: &Vec3Wasm) -> Vec3Wasm {
        // p' = q * p * q*
        let qv = Vec3Wasm::new(self.x, self.y, self.z);
        let uv = qv.cross(v);
        let uuv = qv.cross(&uv);
        let scaled_uv = uv.scale(2.0 * self.w);
        let scaled_uuv = uuv.scale(2.0);
        v.add(&scaled_uv).add(&scaled_uuv)
    }

    /// Spherical linear interpolation (SLERP) toward `other` by factor `t`.
    pub fn slerp(&self, other: &QuatWasm, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        let mut dot = self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w;

        // If dot negative, negate one quaternion to take the shorter arc
        let other_adj = if dot < 0.0 {
            dot = -dot;
            QuatWasm::new(-other.x, -other.y, -other.z, -other.w)
        } else {
            *other
        };

        if dot > 0.9995 {
            // Nearly identical: linear interpolation
            let r = QuatWasm {
                x: self.x + t * (other_adj.x - self.x),
                y: self.y + t * (other_adj.y - self.y),
                z: self.z + t * (other_adj.z - self.z),
                w: self.w + t * (other_adj.w - self.w),
            };
            return r.normalized();
        }

        let theta_0 = dot.acos();
        let theta = theta_0 * t;
        let sin_theta = theta.sin();
        let sin_theta_0 = theta_0.sin();

        let s1 = theta.cos() - dot * sin_theta / sin_theta_0;
        let s2 = sin_theta / sin_theta_0;
        QuatWasm {
            x: s1 * self.x + s2 * other_adj.x,
            y: s1 * self.y + s2 * other_adj.y,
            z: s1 * self.z + s2 * other_adj.z,
            w: s1 * self.w + s2 * other_adj.w,
        }
    }
}

impl Default for QuatWasm {
    fn default() -> Self {
        Self::identity()
    }
}

// ---------------------------------------------------------------------------
// Transform
// ---------------------------------------------------------------------------

/// A rigid body transform: position + orientation.
///
/// Stored as a `Vec3Wasm` position and a `QuatWasm` orientation to avoid
/// the overhead of a 4x4 matrix while still supporting full 3D transformations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct TransformWasm {
    /// Position in world space `[x, y, z]`.
    pub position: Vec3Wasm,
    /// Orientation as a unit quaternion `[x, y, z, w]`.
    pub rotation: QuatWasm,
}

impl TransformWasm {
    /// Create a transform from position and rotation.
    pub fn new(position: Vec3Wasm, rotation: QuatWasm) -> Self {
        Self { position, rotation }
    }

    /// Create a transform at the given position with identity rotation.
    pub fn from_position(x: f64, y: f64, z: f64) -> Self {
        Self {
            position: Vec3Wasm::new(x, y, z),
            rotation: QuatWasm::identity(),
        }
    }

    /// The identity transform (origin, no rotation).
    pub fn identity() -> Self {
        Self {
            position: Vec3Wasm::zero(),
            rotation: QuatWasm::identity(),
        }
    }

    /// Encode this transform as a flat `[px, py, pz, qx, qy, qz, qw]` array.
    pub fn to_array7(&self) -> [f64; 7] {
        [
            self.position.x,
            self.position.y,
            self.position.z,
            self.rotation.x,
            self.rotation.y,
            self.rotation.z,
            self.rotation.w,
        ]
    }

    /// Decode from a flat `[px, py, pz, qx, qy, qz, qw]` array.
    pub fn from_array7(arr: [f64; 7]) -> Self {
        Self {
            position: Vec3Wasm::new(arr[0], arr[1], arr[2]),
            rotation: QuatWasm::new(arr[3], arr[4], arr[5], arr[6]),
        }
    }

    /// Encode this transform as a column-major 4x4 matrix `[f64; 16]`.
    pub fn to_matrix4(&self) -> [f64; 16] {
        let q = &self.rotation;
        let (x, y, z, w) = (q.x, q.y, q.z, q.w);
        let (tx, ty, tz) = (self.position.x, self.position.y, self.position.z);

        // Column-major layout
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y + w * z),
            2.0 * (x * z - w * y),
            0.0,
            2.0 * (x * y - w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z + w * x),
            0.0,
            2.0 * (x * z + w * y),
            2.0 * (y * z - w * x),
            1.0 - 2.0 * (x * x + y * y),
            0.0,
            tx,
            ty,
            tz,
            1.0,
        ]
    }

    /// Transform a point from local space to world space.
    pub fn transform_point(&self, p: &Vec3Wasm) -> Vec3Wasm {
        self.rotation.rotate_vec(p).add(&self.position)
    }

    /// Transform a direction vector (ignores translation).
    pub fn transform_vector(&self, v: &Vec3Wasm) -> Vec3Wasm {
        self.rotation.rotate_vec(v)
    }

    /// Compute the inverse transform.
    pub fn inverse(&self) -> Self {
        let inv_rot = self.rotation.conjugate();
        let neg_pos = self.position.scale(-1.0);
        let inv_pos = inv_rot.rotate_vec(&neg_pos);
        Self {
            position: inv_pos,
            rotation: inv_rot,
        }
    }
}

impl Default for TransformWasm {
    fn default() -> Self {
        Self::identity()
    }
}

// ---------------------------------------------------------------------------
// RigidBodyConfig
// ---------------------------------------------------------------------------

/// Configuration for creating a new rigid body.
///
/// All fields have sensible defaults via the `Default` implementation.
/// Static bodies (mass = 0.0) are immovable and do not respond to forces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RigidBodyConfig {
    /// Mass in kilograms. Use `0.0` for a static/kinematic body.
    pub mass: f64,
    /// Initial position `[x, y, z]` in world space.
    pub position: [f64; 3],
    /// Initial orientation as a quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
    /// Initial linear velocity `[x, y, z]`.
    pub linear_velocity: [f64; 3],
    /// Initial angular velocity `[wx, wy, wz]` in radians/second.
    pub angular_velocity: [f64; 3],
    /// Linear damping coefficient (0 = no damping, 1 = fully damped).
    pub linear_damping: f64,
    /// Angular damping coefficient (0 = no damping, 1 = fully damped).
    pub angular_damping: f64,
    /// Restitution (bounciness) coefficient in the range \[0, 1\].
    /// 0 = perfectly inelastic, 1 = perfectly elastic.
    pub restitution: f64,
    /// Friction coefficient (Coulomb model). Typical values: 0.3–0.8.
    pub friction: f64,
    /// Whether this body can go to sleep when it has been at rest.
    pub can_sleep: bool,
    /// Whether this is a sensor (triggers contact events but exerts no forces).
    pub is_sensor: bool,
    /// Whether the body is kinematic (position-driven, not force-driven).
    pub is_kinematic: bool,
    /// Gravity scale multiplier (1.0 = normal gravity, 0.0 = no gravity).
    pub gravity_scale: f64,
}

impl Default for RigidBodyConfig {
    fn default() -> Self {
        Self {
            mass: 1.0,
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            linear_damping: 0.01,
            angular_damping: 0.05,
            restitution: 0.3,
            friction: 0.5,
            can_sleep: true,
            is_sensor: false,
            is_kinematic: false,
            gravity_scale: 1.0,
        }
    }
}

impl RigidBodyConfig {
    /// Create a config with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a static body config at the origin.
    pub fn static_body() -> Self {
        Self {
            mass: 0.0,
            can_sleep: false,
            ..Default::default()
        }
    }

    /// Set the mass.
    pub fn with_mass(mut self, mass: f64) -> Self {
        self.mass = mass;
        self
    }

    /// Set the initial position.
    pub fn with_position(mut self, x: f64, y: f64, z: f64) -> Self {
        self.position = [x, y, z];
        self
    }

    /// Set the initial linear velocity.
    pub fn with_linear_velocity(mut self, vx: f64, vy: f64, vz: f64) -> Self {
        self.linear_velocity = [vx, vy, vz];
        self
    }

    /// Set the restitution coefficient.
    pub fn with_restitution(mut self, r: f64) -> Self {
        self.restitution = r.clamp(0.0, 1.0);
        self
    }

    /// Set the friction coefficient.
    pub fn with_friction(mut self, f: f64) -> Self {
        self.friction = f.max(0.0);
        self
    }

    /// Mark this body as a sensor (no collision response).
    pub fn as_sensor(mut self) -> Self {
        self.is_sensor = true;
        self
    }

    /// Mark this body as kinematic (position-controlled).
    pub fn as_kinematic(mut self) -> Self {
        self.is_kinematic = true;
        self.mass = 0.0;
        self
    }

    /// Returns `true` if this body is static (mass == 0 and not kinematic).
    pub fn is_static(&self) -> bool {
        self.mass == 0.0 && !self.is_kinematic
    }
}

// ---------------------------------------------------------------------------
// ColliderConfig
// ---------------------------------------------------------------------------

/// Shape type for a collider.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColliderShapeType {
    /// Sphere shape, parameterized by radius.
    Sphere,
    /// Axis-aligned box shape, parameterized by half-extents.
    Box,
    /// Infinite static plane, defined by normal + offset.
    Plane,
    /// Capsule shape (cylinder capped with hemispheres).
    Capsule,
    /// Cylinder shape.
    Cylinder,
    /// Convex hull shape (user-provided vertex list).
    ConvexHull,
}

/// Configuration for creating a new collider attached to a body.
///
/// A collider defines the shape used for collision detection.
/// One body can have multiple colliders (compound shapes).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColliderConfig {
    /// Shape type of this collider.
    pub shape_type: ColliderShapeType,
    /// Radius (used for `Sphere`, `Capsule`, `Cylinder`).
    pub radius: f64,
    /// Half-extents `[hx, hy, hz]` (used for `Box`).
    pub half_extents: [f64; 3],
    /// Height (used for `Capsule` and `Cylinder`).
    pub height: f64,
    /// Plane normal `[nx, ny, nz]` (used for `Plane`).
    pub plane_normal: [f64; 3],
    /// Plane offset along normal (used for `Plane`).
    pub plane_offset: f64,
    /// Local transform of this collider relative to the body origin.
    pub local_position: [f64; 3],
    /// Local rotation as a quaternion `[x, y, z, w]`.
    pub local_rotation: [f64; 4],
    /// Density used to compute mass properties (kg/m³).
    pub density: f64,
    /// Whether this collider is a sensor only (no collision response).
    pub is_sensor: bool,
    /// Collision group bitmask (bodies only collide if groups match).
    pub collision_group: u32,
    /// Collision mask bitmask (which groups this collider interacts with).
    pub collision_mask: u32,
}

impl Default for ColliderConfig {
    fn default() -> Self {
        Self {
            shape_type: ColliderShapeType::Sphere,
            radius: 0.5,
            half_extents: [0.5; 3],
            height: 1.0,
            plane_normal: [0.0, 1.0, 0.0],
            plane_offset: 0.0,
            local_position: [0.0; 3],
            local_rotation: [0.0, 0.0, 0.0, 1.0],
            density: 1000.0,
            is_sensor: false,
            collision_group: 0xFFFF_FFFF,
            collision_mask: 0xFFFF_FFFF,
        }
    }
}

impl ColliderConfig {
    /// Create a sphere collider with the given radius.
    pub fn sphere(radius: f64) -> Self {
        Self {
            shape_type: ColliderShapeType::Sphere,
            radius,
            ..Default::default()
        }
    }

    /// Create a box collider with the given half-extents.
    pub fn cuboid(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            shape_type: ColliderShapeType::Box,
            half_extents: [hx, hy, hz],
            ..Default::default()
        }
    }

    /// Create an infinite static plane collider.
    ///
    /// The plane's equation is `normal · x = offset`.
    pub fn plane(nx: f64, ny: f64, nz: f64, offset: f64) -> Self {
        Self {
            shape_type: ColliderShapeType::Plane,
            plane_normal: [nx, ny, nz],
            plane_offset: offset,
            ..Default::default()
        }
    }

    /// Create a capsule collider.
    pub fn capsule(radius: f64, height: f64) -> Self {
        Self {
            shape_type: ColliderShapeType::Capsule,
            radius,
            height,
            ..Default::default()
        }
    }

    /// Create a cylinder collider.
    pub fn cylinder(radius: f64, height: f64) -> Self {
        Self {
            shape_type: ColliderShapeType::Cylinder,
            radius,
            height,
            ..Default::default()
        }
    }

    /// Compute the approximate volume of this collider shape.
    pub fn volume(&self) -> f64 {
        use std::f64::consts::PI;
        match self.shape_type {
            ColliderShapeType::Sphere => (4.0 / 3.0) * PI * self.radius.powi(3),
            ColliderShapeType::Box => {
                8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
            }
            ColliderShapeType::Capsule => {
                PI * self.radius * self.radius * self.height
                    + (4.0 / 3.0) * PI * self.radius.powi(3)
            }
            ColliderShapeType::Cylinder => PI * self.radius * self.radius * self.height,
            ColliderShapeType::Plane | ColliderShapeType::ConvexHull => 0.0,
        }
    }

    /// Compute the approximate mass of this collider given `density`.
    pub fn mass(&self) -> f64 {
        self.volume() * self.density
    }
}

// ---------------------------------------------------------------------------
// ContactResult
// ---------------------------------------------------------------------------

/// The result of a contact query between two bodies.
///
/// Contains the full contact manifold information needed to resolve
/// the collision and trigger physics events on the JavaScript side.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactResult {
    /// Handle of the first body involved in the contact.
    pub body_a: u32,
    /// Handle of the second body involved in the contact.
    pub body_b: u32,
    /// Contact point in world space on body A's surface.
    pub point_on_a: [f64; 3],
    /// Contact point in world space on body B's surface.
    pub point_on_b: [f64; 3],
    /// Contact normal pointing from body B toward body A.
    pub normal: [f64; 3],
    /// Penetration depth (positive means overlapping).
    pub depth: f64,
    /// Relative velocity at the contact point along the normal.
    pub relative_velocity: f64,
    /// Impulse applied to resolve this contact (magnitude).
    pub impulse: f64,
    /// Whether this is a new contact (vs. a persisting one from the previous frame).
    pub is_new: bool,
    /// Friction impulse applied at this contact.
    pub friction_impulse: f64,
}

impl ContactResult {
    /// Create a new `ContactResult` with default values.
    pub fn new(body_a: u32, body_b: u32) -> Self {
        Self {
            body_a,
            body_b,
            point_on_a: [0.0; 3],
            point_on_b: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.0,
            relative_velocity: 0.0,
            impulse: 0.0,
            is_new: true,
            friction_impulse: 0.0,
        }
    }

    /// Return `true` if the contact is penetrating (depth > threshold).
    pub fn is_penetrating(&self) -> bool {
        self.depth > 1e-6
    }

    /// Return the mid-point between the two contact points.
    pub fn contact_midpoint(&self) -> [f64; 3] {
        [
            (self.point_on_a[0] + self.point_on_b[0]) * 0.5,
            (self.point_on_a[1] + self.point_on_b[1]) * 0.5,
            (self.point_on_a[2] + self.point_on_b[2]) * 0.5,
        ]
    }
}

// ---------------------------------------------------------------------------
// SimulationConfig
// ---------------------------------------------------------------------------

/// Global simulation configuration passed when creating the engine.
///
/// This covers gravity, solver settings, time stepping, and sleeping parameters.
/// Changing these at runtime requires creating a new engine (or calling `reset`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Gravity vector `[gx, gy, gz]` in m/s². Default: `[0, -9.81, 0]`.
    pub gravity: [f64; 3],
    /// Fixed simulation time step in seconds. Default: `1/60`.
    pub fixed_dt: f64,
    /// Maximum number of substeps per `step()` call.
    pub max_substeps: u32,
    /// Number of constraint solver iterations. Higher = more accurate.
    pub solver_iterations: u32,
    /// Linear velocity below which a body may go to sleep (m/s).
    pub linear_sleep_threshold: f64,
    /// Angular velocity below which a body may go to sleep (rad/s).
    pub angular_sleep_threshold: f64,
    /// Time a body must stay below sleep thresholds before sleeping (s).
    pub time_before_sleep: f64,
    /// Enable body sleeping to save CPU when bodies come to rest.
    pub sleeping_enabled: bool,
    /// Enable continuous collision detection (CCD) to prevent tunneling.
    pub ccd_enabled: bool,
    /// Broadphase margin added to AABBs (m).
    pub broadphase_margin: f64,
    /// Number of position-correction iterations (for PBD / XPBD).
    pub position_iterations: u32,
    /// Enable warm starting of constraints (reuse impulses from last frame).
    pub warm_starting: bool,
    /// Contact slop distance: penetration allowed before correction (m).
    pub contact_slop: f64,
    /// Maximum allowed penetration depth before full correction (m).
    pub max_penetration_correction: f64,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            fixed_dt: 1.0 / 60.0,
            max_substeps: 4,
            solver_iterations: 8,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            sleeping_enabled: true,
            ccd_enabled: false,
            broadphase_margin: 0.02,
            position_iterations: 4,
            warm_starting: true,
            contact_slop: 0.005,
            max_penetration_correction: 0.2,
        }
    }
}

impl SimulationConfig {
    /// Create a config with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a high-accuracy config (more solver iterations, smaller time step).
    pub fn high_accuracy() -> Self {
        Self {
            fixed_dt: 1.0 / 120.0,
            max_substeps: 8,
            solver_iterations: 16,
            position_iterations: 8,
            ..Default::default()
        }
    }

    /// Create a performance config (fewer iterations, larger time step).
    pub fn performance() -> Self {
        Self {
            fixed_dt: 1.0 / 30.0,
            max_substeps: 2,
            solver_iterations: 4,
            position_iterations: 2,
            ccd_enabled: false,
            ..Default::default()
        }
    }

    /// Set gravity.
    pub fn with_gravity(mut self, gx: f64, gy: f64, gz: f64) -> Self {
        self.gravity = [gx, gy, gz];
        self
    }

    /// Set the fixed time step.
    pub fn with_fixed_dt(mut self, dt: f64) -> Self {
        self.fixed_dt = dt;
        self
    }

    /// Set the number of solver iterations.
    pub fn with_solver_iterations(mut self, n: u32) -> Self {
        self.solver_iterations = n;
        self
    }

    /// Enable or disable CCD.
    pub fn with_ccd(mut self, enabled: bool) -> Self {
        self.ccd_enabled = enabled;
        self
    }
}

// ---------------------------------------------------------------------------
// BodyState
// ---------------------------------------------------------------------------

/// Full runtime state snapshot of a single rigid body.
///
/// This is returned by `WasmPhysicsEngine::get_body_state()` and can be
/// serialized to JSON for web worker messaging.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyState {
    /// Body handle (index).
    pub handle: u32,
    /// Current position `[x, y, z]`.
    pub position: [f64; 3],
    /// Current orientation as a quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
    /// Current linear velocity `[vx, vy, vz]`.
    pub linear_velocity: [f64; 3],
    /// Current angular velocity `[wx, wy, wz]`.
    pub angular_velocity: [f64; 3],
    /// Whether the body is currently sleeping.
    pub is_sleeping: bool,
    /// Whether the body is active (not removed).
    pub is_active: bool,
    /// Kinetic energy (approximation: 0.5 * m * v²).
    pub kinetic_energy: f64,
}

impl BodyState {
    /// Create a new `BodyState` with default zero values.
    pub fn new(handle: u32) -> Self {
        Self {
            handle,
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            is_sleeping: false,
            is_active: true,
            kinetic_energy: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// RaycastResult
// ---------------------------------------------------------------------------

/// Result of a raycast query.
///
/// If no body was hit, `hit` is `false` and other fields are undefined.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RaycastResult {
    /// Whether the ray hit any body.
    pub hit: bool,
    /// Handle of the body that was hit (only valid when `hit == true`).
    pub body_handle: u32,
    /// World-space hit point `[x, y, z]`.
    pub point: [f64; 3],
    /// Surface normal at the hit point `[nx, ny, nz]`.
    pub normal: [f64; 3],
    /// Distance from ray origin to hit point.
    pub distance: f64,
}

impl RaycastResult {
    /// Create a "no hit" result.
    pub fn no_hit() -> Self {
        Self {
            hit: false,
            body_handle: u32::MAX,
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            distance: f64::INFINITY,
        }
    }
}

// ---------------------------------------------------------------------------
// DebugInfo
// ---------------------------------------------------------------------------

/// Debug and performance information from the last simulation step.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DebugInfo {
    /// Number of broadphase pairs tested this step.
    pub broadphase_pairs: u32,
    /// Number of narrowphase tests performed this step.
    pub narrowphase_tests: u32,
    /// Number of active contacts in the manifold cache.
    pub active_contacts: u32,
    /// Number of sleeping bodies.
    pub sleeping_bodies: u32,
    /// Number of solver iterations actually performed.
    pub solver_iterations_performed: u32,
    /// Wall-clock time for the integration step (microseconds).
    pub integration_time_us: u64,
    /// Wall-clock time for the collision detection step (microseconds).
    pub collision_time_us: u64,
    /// Wall-clock time for the constraint solver step (microseconds).
    pub solver_time_us: u64,
    /// Total wall-clock time for the simulation step (microseconds).
    pub total_time_us: u64,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- Vec3Wasm ---

    #[test]
    fn test_vec3_zero() {
        let v = Vec3Wasm::zero();
        assert_eq!(v.length(), 0.0);
    }

    #[test]
    fn test_vec3_length() {
        let v = Vec3Wasm::new(3.0, 4.0, 0.0);
        assert!((v.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_normalized() {
        let v = Vec3Wasm::new(0.0, 3.0, 0.0);
        let n = v.normalized();
        assert!((n.y - 1.0).abs() < 1e-10);
        assert!(n.x.abs() < 1e-10);
    }

    #[test]
    fn test_vec3_dot() {
        let a = Vec3Wasm::new(1.0, 0.0, 0.0);
        let b = Vec3Wasm::new(0.0, 1.0, 0.0);
        assert!((a.dot(&b)).abs() < 1e-10);
        assert!((a.dot(&a) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_cross() {
        let x = Vec3Wasm::unit_x();
        let y = Vec3Wasm::unit_y();
        let z = x.cross(&y);
        assert!((z.z - 1.0).abs() < 1e-10);
        assert!(z.x.abs() < 1e-10);
        assert!(z.y.abs() < 1e-10);
    }

    #[test]
    fn test_vec3_lerp() {
        let a = Vec3Wasm::new(0.0, 0.0, 0.0);
        let b = Vec3Wasm::new(10.0, 0.0, 0.0);
        let mid = a.lerp(&b, 0.5);
        assert!((mid.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_array_roundtrip() {
        let arr = [1.0, 2.0, 3.0];
        let v = Vec3Wasm::from_array(arr);
        assert_eq!(v.to_array(), arr);
    }

    // --- QuatWasm ---

    #[test]
    fn test_quat_identity_rotate() {
        let q = QuatWasm::identity();
        let v = Vec3Wasm::new(1.0, 2.0, 3.0);
        let rotated = q.rotate_vec(&v);
        assert!((rotated.x - v.x).abs() < 1e-10);
        assert!((rotated.y - v.y).abs() < 1e-10);
        assert!((rotated.z - v.z).abs() < 1e-10);
    }

    #[test]
    fn test_quat_mul_identity() {
        let q = QuatWasm::new(0.1, 0.2, 0.3, 0.9).normalized();
        let i = QuatWasm::identity();
        let qi = q.mul(&i);
        assert!((qi.x - q.x).abs() < 1e-10);
        assert!((qi.w - q.w).abs() < 1e-10);
    }

    #[test]
    fn test_quat_slerp_endpoints() {
        let a = QuatWasm::identity();
        let b = QuatWasm::new(0.0, 1.0, 0.0, 0.0).normalized();
        let r0 = a.slerp(&b, 0.0);
        let r1 = a.slerp(&b, 1.0);
        assert!((r0.w - 1.0).abs() < 1e-6);
        assert!((r1.w).abs() < 1e-6);
    }

    // --- TransformWasm ---

    #[test]
    fn test_transform_identity() {
        let t = TransformWasm::identity();
        let p = Vec3Wasm::new(1.0, 2.0, 3.0);
        let out = t.transform_point(&p);
        assert!((out.x - p.x).abs() < 1e-10);
        assert!((out.y - p.y).abs() < 1e-10);
        assert!((out.z - p.z).abs() < 1e-10);
    }

    #[test]
    fn test_transform_translation() {
        let t = TransformWasm::from_position(5.0, 0.0, 0.0);
        let p = Vec3Wasm::zero();
        let out = t.transform_point(&p);
        assert!((out.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_inverse() {
        let t = TransformWasm::from_position(3.0, -2.0, 1.0);
        let inv = t.inverse();
        let p = Vec3Wasm::zero();
        let roundtrip = TransformWasm::from_position(
            t.inverse().position.x,
            t.inverse().position.y,
            t.inverse().position.z,
        );
        // inv of translation should be negative translation
        let out = inv.transform_point(&p);
        let original = t.transform_point(&out);
        let _ = roundtrip;
        // original should be near zero
        assert!(original.length() < 1e-6, "roundtrip error: {:?}", original);
    }

    #[test]
    fn test_transform_matrix4() {
        let t = TransformWasm::identity();
        let m = t.to_matrix4();
        // Diagonal elements should be 1 for identity
        assert!((m[0] - 1.0).abs() < 1e-10); // col0[0]
        assert!((m[5] - 1.0).abs() < 1e-10); // col1[1]
        assert!((m[10] - 1.0).abs() < 1e-10); // col2[2]
        assert!((m[15] - 1.0).abs() < 1e-10); // col3[3]
    }

    #[test]
    fn test_transform_array7_roundtrip() {
        let t = TransformWasm::from_position(1.0, 2.0, 3.0);
        let arr = t.to_array7();
        let t2 = TransformWasm::from_array7(arr);
        assert!((t2.position.x - 1.0).abs() < 1e-10);
        assert!((t2.position.y - 2.0).abs() < 1e-10);
        assert!((t2.position.z - 3.0).abs() < 1e-10);
    }

    // --- RigidBodyConfig ---

    #[test]
    fn test_rigid_body_config_default() {
        let cfg = RigidBodyConfig::default();
        assert!((cfg.mass - 1.0).abs() < 1e-10);
        assert!(!cfg.is_static());
    }

    #[test]
    fn test_rigid_body_config_static() {
        let cfg = RigidBodyConfig::static_body();
        assert!(cfg.is_static());
    }

    #[test]
    fn test_rigid_body_config_builder() {
        let cfg = RigidBodyConfig::new()
            .with_mass(5.0)
            .with_position(1.0, 2.0, 3.0)
            .with_restitution(0.8)
            .with_friction(0.4);
        assert!((cfg.mass - 5.0).abs() < 1e-10);
        assert!((cfg.position[0] - 1.0).abs() < 1e-10);
        assert!((cfg.restitution - 0.8).abs() < 1e-10);
        assert!((cfg.friction - 0.4).abs() < 1e-10);
    }

    // --- ColliderConfig ---

    #[test]
    fn test_collider_sphere_volume() {
        let c = ColliderConfig::sphere(1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((c.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_collider_box_volume() {
        let c = ColliderConfig::cuboid(1.0, 2.0, 3.0);
        assert!((c.volume() - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_plane() {
        let c = ColliderConfig::plane(0.0, 1.0, 0.0, 0.0);
        assert_eq!(c.shape_type, ColliderShapeType::Plane);
        assert!((c.plane_normal[1] - 1.0).abs() < 1e-10);
    }

    // --- ContactResult ---

    #[test]
    fn test_contact_result_new() {
        let cr = ContactResult::new(0, 1);
        assert_eq!(cr.body_a, 0);
        assert_eq!(cr.body_b, 1);
        assert!(cr.is_new);
    }

    #[test]
    fn test_contact_midpoint() {
        let mut cr = ContactResult::new(0, 1);
        cr.point_on_a = [0.0, 0.0, 0.0];
        cr.point_on_b = [2.0, 0.0, 0.0];
        let mid = cr.contact_midpoint();
        assert!((mid[0] - 1.0).abs() < 1e-10);
    }

    // --- SimulationConfig ---

    #[test]
    fn test_simulation_config_default() {
        let cfg = SimulationConfig::default();
        assert!((cfg.gravity[1] + 9.81).abs() < 1e-10);
        assert!((cfg.fixed_dt - 1.0 / 60.0).abs() < 1e-10);
    }

    #[test]
    fn test_simulation_config_high_accuracy() {
        let cfg = SimulationConfig::high_accuracy();
        assert!(cfg.solver_iterations > SimulationConfig::default().solver_iterations);
    }

    #[test]
    fn test_simulation_config_builder() {
        let cfg = SimulationConfig::new()
            .with_gravity(0.0, -1.62, 0.0) // Moon gravity
            .with_solver_iterations(16)
            .with_ccd(true);
        assert!((cfg.gravity[1] + 1.62).abs() < 1e-10);
        assert_eq!(cfg.solver_iterations, 16);
        assert!(cfg.ccd_enabled);
    }

    // --- BodyState ---

    #[test]
    fn test_body_state_new() {
        let s = BodyState::new(42);
        assert_eq!(s.handle, 42);
        assert!(s.is_active);
        assert!(!s.is_sleeping);
    }

    // --- RaycastResult ---

    #[test]
    fn test_raycast_result_no_hit() {
        let r = RaycastResult::no_hit();
        assert!(!r.hit);
        assert_eq!(r.distance, f64::INFINITY);
    }

    // --- DebugInfo ---

    #[test]
    fn test_debug_info_default() {
        let d = DebugInfo::default();
        assert_eq!(d.active_contacts, 0);
        assert_eq!(d.total_time_us, 0);
    }
}
