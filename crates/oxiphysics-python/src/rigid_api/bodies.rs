// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body and collider types for Python bindings.

use super::{quat_identity, vec3_dot, vec3_length, vec3_scale};
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// PyShapeType
// ---------------------------------------------------------------------------

/// Collider shape enumeration.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PyShapeType {
    /// Sphere shape.
    Sphere = 0,
    /// Axis-aligned box shape.
    Box = 1,
    /// Capsule shape (cylinder with hemispherical caps).
    Capsule = 2,
    /// Convex hull mesh.
    ConvexHull = 3,
    /// Triangle mesh (static only).
    TriMesh = 4,
}

// ---------------------------------------------------------------------------
// PyCollider
// ---------------------------------------------------------------------------

/// A collider attached to a rigid body.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyCollider {
    /// Shape type.
    #[pyo3(get, set)]
    pub shape: PyShapeType,
    /// Half-extents `[hx, hy, hz]` for Box; `[radius, half_height, 0]` for Capsule.
    #[pyo3(get, set)]
    pub half_extents: [f64; 3],
    /// Radius (used for Sphere and Capsule).
    #[pyo3(get, set)]
    pub radius: f64,
    /// Friction coefficient.
    #[pyo3(get, set)]
    pub friction: f64,
    /// Restitution coefficient.
    #[pyo3(get, set)]
    pub restitution: f64,
    /// Local offset `[x, y, z]` from the body centre.
    #[pyo3(get, set)]
    pub offset: [f64; 3],
    /// Whether this collider acts as a sensor (no impulse response).
    #[pyo3(get, set)]
    pub is_sensor: bool,
}

#[pymethods]
impl PyCollider {
    /// Create a sphere collider.
    #[new]
    #[pyo3(signature = (radius = 0.5, friction = 0.5, restitution = 0.3))]
    pub fn new(radius: f64, friction: f64, restitution: f64) -> Self {
        Self {
            shape: PyShapeType::Sphere,
            half_extents: [radius, radius, radius],
            radius,
            friction,
            restitution,
            offset: [0.0; 3],
            is_sensor: false,
        }
    }

    /// Create a box collider with the given half-extents.
    #[staticmethod]
    pub fn box_collider(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            shape: PyShapeType::Box,
            half_extents: [hx, hy, hz],
            radius: hx.min(hy).min(hz),
            friction: 0.5,
            restitution: 0.3,
            offset: [0.0; 3],
            is_sensor: false,
        }
    }

    /// Create a capsule collider.
    #[staticmethod]
    pub fn capsule_collider(radius: f64, half_height: f64) -> Self {
        Self {
            shape: PyShapeType::Capsule,
            half_extents: [radius, half_height, 0.0],
            radius,
            friction: 0.5,
            restitution: 0.3,
            offset: [0.0; 3],
            is_sensor: false,
        }
    }

    /// Volume estimate for this collider shape.
    pub fn volume(&self) -> f64 {
        match self.shape {
            PyShapeType::Sphere => (4.0 / 3.0) * std::f64::consts::PI * self.radius.powi(3),
            PyShapeType::Box => {
                8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
            }
            PyShapeType::Capsule => {
                let r = self.radius;
                let h = self.half_extents[1] * 2.0;
                std::f64::consts::PI * r * r * h + (4.0 / 3.0) * std::f64::consts::PI * r * r * r
            }
            _ => 0.0,
        }
    }

    fn __repr__(&self) -> String {
        format!("PyCollider(shape={:?}, radius={})", self.shape, self.radius)
    }
}

// ---------------------------------------------------------------------------
// PyContact
// ---------------------------------------------------------------------------

/// Contact information returned after a collision step.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyContact {
    /// Handle of body A.
    #[pyo3(get)]
    pub body_a: u32,
    /// Handle of body B.
    #[pyo3(get)]
    pub body_b: u32,
    /// World-space contact point.
    #[pyo3(get)]
    pub contact_point: [f64; 3],
    /// Contact normal (pointing from B toward A).
    #[pyo3(get)]
    pub normal: [f64; 3],
    /// Penetration depth.
    #[pyo3(get)]
    pub depth: f64,
    /// Impulse magnitude applied during resolution.
    #[pyo3(get)]
    pub impulse: f64,
}

#[pymethods]
impl PyContact {
    #[new]
    pub fn new(
        body_a: u32,
        body_b: u32,
        contact_point: [f64; 3],
        normal: [f64; 3],
        depth: f64,
        impulse: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            contact_point,
            normal,
            depth,
            impulse,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyContact(a={}, b={}, depth={:.4})",
            self.body_a, self.body_b, self.depth
        )
    }
}

// ---------------------------------------------------------------------------
// PyRayResult
// ---------------------------------------------------------------------------

/// Result of a ray-cast query.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyRayResult {
    /// Handle of the body that was hit.
    #[pyo3(get)]
    pub handle: u32,
    /// Distance along the ray to the hit point.
    #[pyo3(get)]
    pub distance: f64,
    /// World-space hit point.
    #[pyo3(get)]
    pub hit_point: [f64; 3],
    /// Surface normal at the hit point.
    #[pyo3(get)]
    pub normal: [f64; 3],
}

#[pymethods]
impl PyRayResult {
    #[new]
    pub fn new(handle: u32, distance: f64, hit_point: [f64; 3], normal: [f64; 3]) -> Self {
        Self {
            handle,
            distance,
            hit_point,
            normal,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyRayResult(handle={}, distance={:.4})",
            self.handle, self.distance
        )
    }
}

// ---------------------------------------------------------------------------
// PyRigidBody
// ---------------------------------------------------------------------------

/// A single rigid body in the physics world.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyRigidBody {
    /// Unique body handle.
    #[pyo3(get)]
    pub handle: u32,
    /// Mass in kilograms.
    #[pyo3(get, set)]
    pub mass: f64,
    /// World-space position `[x, y, z]`.
    #[pyo3(get, set)]
    pub position: [f64; 3],
    /// Linear velocity `[vx, vy, vz]`.
    #[pyo3(get, set)]
    pub velocity: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    #[pyo3(get, set)]
    pub orientation: [f64; 4],
    /// Angular velocity `[wx, wy, wz]`.
    #[pyo3(get, set)]
    pub angular_velocity: [f64; 3],
    /// Whether this body is static (immovable).
    #[pyo3(get, set)]
    pub is_static: bool,
    /// Whether this body is kinematic.
    #[pyo3(get, set)]
    pub is_kinematic: bool,
    /// Whether sleep is enabled for this body.
    #[pyo3(get, set)]
    pub can_sleep: bool,
    /// Friction coefficient.
    #[pyo3(get, set)]
    pub friction: f64,
    /// Restitution coefficient.
    #[pyo3(get, set)]
    pub restitution: f64,
    /// Linear damping coefficient.
    #[pyo3(get, set)]
    pub linear_damping: f64,
    /// Angular damping coefficient.
    #[pyo3(get, set)]
    pub angular_damping: f64,
    /// Optional string tag.
    #[pyo3(get, set)]
    pub tag: Option<String>,
    /// Attached colliders.
    #[pyo3(get)]
    pub colliders: Vec<PyCollider>,
}

#[pymethods]
impl PyRigidBody {
    /// Create a new dynamic rigid body.
    #[new]
    #[pyo3(signature = (handle = 0, mass = 1.0, position = [0.0, 0.0, 0.0]))]
    pub fn new(handle: u32, mass: f64, position: [f64; 3]) -> Self {
        Self {
            handle,
            mass,
            position,
            velocity: [0.0; 3],
            orientation: quat_identity(),
            angular_velocity: [0.0; 3],
            is_static: false,
            is_kinematic: false,
            can_sleep: true,
            friction: 0.5,
            restitution: 0.3,
            linear_damping: 0.0,
            angular_damping: 0.0,
            tag: None,
            colliders: Vec::new(),
        }
    }

    /// Create a static (immovable) body.
    #[staticmethod]
    pub fn static_body(handle: u32, position: [f64; 3]) -> Self {
        let mut body = Self::new(handle, 0.0, position);
        body.is_static = true;
        body
    }

    /// Add a collider to this body.
    pub fn add_collider(&mut self, collider: PyCollider) {
        self.colliders.push(collider);
    }

    /// Compute kinetic energy (½mv²).
    pub fn kinetic_energy(&self) -> f64 {
        if self.is_static || self.mass <= 0.0 {
            return 0.0;
        }
        let v2 = vec3_dot(self.velocity, self.velocity);
        0.5 * self.mass * v2
    }

    /// Speed (magnitude of velocity vector).
    pub fn speed(&self) -> f64 {
        vec3_length(self.velocity)
    }

    fn __repr__(&self) -> String {
        format!(
            "PyRigidBody(handle={}, mass={}, pos={:?})",
            self.handle, self.mass, self.position
        )
    }
}

// ---------------------------------------------------------------------------
// PyRigidBodySet
// ---------------------------------------------------------------------------

/// A collection of rigid bodies.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct PyRigidBodySet {
    pub(crate) bodies: Vec<PyRigidBody>,
    pub(crate) next_handle: u32,
}

#[pymethods]
impl PyRigidBodySet {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a body to the set. Returns its assigned handle.
    pub fn add(&mut self, mut body: PyRigidBody) -> u32 {
        let h = self.next_handle;
        body.handle = h;
        self.bodies.push(body);
        self.next_handle += 1;
        h
    }

    /// Remove a body by handle. Returns `True` if found and removed.
    pub fn remove(&mut self, handle: u32) -> bool {
        let before = self.bodies.len();
        self.bodies.retain(|b| b.handle != handle);
        self.bodies.len() < before
    }

    /// Get a body by handle, or `None`.
    pub fn get(&self, handle: u32) -> Option<PyRigidBody> {
        self.bodies.iter().find(|b| b.handle == handle).cloned()
    }

    /// Return all bodies as a list.
    pub fn all_bodies(&self) -> Vec<PyRigidBody> {
        self.bodies.clone()
    }

    /// Number of bodies in the set.
    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    /// Apply an impulse to a body.
    pub fn apply_impulse(&mut self, handle: u32, impulse: [f64; 3]) {
        if let Some(body) = self.bodies.iter_mut().find(|b| b.handle == handle)
            && !body.is_static
            && !body.is_kinematic
            && body.mass > 0.0
        {
            let inv_m = 1.0 / body.mass;
            body.velocity = [
                body.velocity[0] + impulse[0] * inv_m,
                body.velocity[1] + impulse[1] * inv_m,
                body.velocity[2] + impulse[2] * inv_m,
            ];
        }
    }

    /// Return handles of all bodies within `radius` of `center`.
    pub fn bodies_in_radius(&self, center: [f64; 3], radius: f64) -> Vec<u32> {
        let r2 = radius * radius;
        self.bodies
            .iter()
            .filter(|b| {
                let dp = vec3_scale(
                    [
                        b.position[0] - center[0],
                        b.position[1] - center[1],
                        b.position[2] - center[2],
                    ],
                    1.0,
                );
                vec3_dot(dp, dp) <= r2
            })
            .map(|b| b.handle)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!("PyRigidBodySet(len={})", self.bodies.len())
    }
}
