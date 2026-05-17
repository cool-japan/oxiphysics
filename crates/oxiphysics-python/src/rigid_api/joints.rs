// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint types and physics pipeline for Python bindings.

use super::bodies::{PyContact, PyRigidBodySet};
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// PyJointType
// ---------------------------------------------------------------------------

/// Type of joint between two rigid bodies.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PyJointType {
    /// Ball-socket: rotation in all directions, no translation.
    BallSocket = 0,
    /// Hinge: rotation around one axis only.
    Hinge = 1,
    /// Slider: translation along one axis only.
    Slider = 2,
    /// Fixed: no relative motion.
    Fixed = 3,
    /// Distance: maintains a fixed or bounded distance.
    Distance = 4,
}

// ---------------------------------------------------------------------------
// PyJoint
// ---------------------------------------------------------------------------

/// A constraint joint between two rigid bodies.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyJoint {
    /// Unique joint handle.
    #[pyo3(get)]
    pub handle: u32,
    /// Type of joint.
    #[pyo3(get)]
    pub joint_type: PyJointType,
    /// Handle of body A.
    #[pyo3(get)]
    pub body_a: u32,
    /// Handle of body B.
    #[pyo3(get)]
    pub body_b: u32,
    /// Anchor in world space (or body-A local if `local_anchors = true`).
    #[pyo3(get, set)]
    pub anchor_a: [f64; 3],
    /// Anchor on body B in world space (or body-B local).
    #[pyo3(get, set)]
    pub anchor_b: [f64; 3],
    /// Hinge / slider axis in world space.
    #[pyo3(get, set)]
    pub axis: [f64; 3],
    /// Target distance for Distance joints.
    #[pyo3(get, set)]
    pub target_distance: f64,
    /// Joint stiffness `[0, 1]`.
    #[pyo3(get, set)]
    pub stiffness: f64,
    /// Whether the joint is currently active.
    #[pyo3(get, set)]
    pub enabled: bool,
}

#[pymethods]
impl PyJoint {
    /// Create a ball-socket joint at the given pivot point.
    #[new]
    pub fn new(handle: u32, body_a: u32, body_b: u32, pivot: [f64; 3]) -> Self {
        Self {
            handle,
            joint_type: PyJointType::BallSocket,
            body_a,
            body_b,
            anchor_a: pivot,
            anchor_b: pivot,
            axis: [0.0, 1.0, 0.0],
            target_distance: 0.0,
            stiffness: 0.8,
            enabled: true,
        }
    }

    /// Create a hinge joint.
    #[staticmethod]
    pub fn hinge(handle: u32, body_a: u32, body_b: u32, pivot: [f64; 3], axis: [f64; 3]) -> Self {
        Self {
            handle,
            joint_type: PyJointType::Hinge,
            body_a,
            body_b,
            anchor_a: pivot,
            anchor_b: pivot,
            axis,
            target_distance: 0.0,
            stiffness: 0.8,
            enabled: true,
        }
    }

    /// Create a distance joint.
    #[staticmethod]
    pub fn distance(
        handle: u32,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        target: f64,
    ) -> Self {
        Self {
            handle,
            joint_type: PyJointType::Distance,
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            axis: [0.0, 1.0, 0.0],
            target_distance: target,
            stiffness: 0.8,
            enabled: true,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PyJoint(handle={}, type={:?}, a={}, b={})",
            self.handle, self.joint_type, self.body_a, self.body_b
        )
    }
}

// ---------------------------------------------------------------------------
// PyJointSet
// ---------------------------------------------------------------------------

/// A collection of joints.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct PyJointSet {
    joints: Vec<PyJoint>,
    next_handle: u32,
}

#[pymethods]
impl PyJointSet {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a joint to the set. Returns its handle.
    pub fn add(&mut self, mut joint: PyJoint) -> u32 {
        let h = self.next_handle;
        joint.handle = h;
        self.joints.push(joint);
        self.next_handle += 1;
        h
    }

    /// Remove a joint by handle. Returns `True` if found.
    pub fn remove(&mut self, handle: u32) -> bool {
        let before = self.joints.len();
        self.joints.retain(|j| j.handle != handle);
        self.joints.len() < before
    }

    /// Get a joint by handle, or `None`.
    pub fn get(&self, handle: u32) -> Option<PyJoint> {
        self.joints.iter().find(|j| j.handle == handle).cloned()
    }

    /// Return all joints.
    pub fn all_joints(&self) -> Vec<PyJoint> {
        self.joints.clone()
    }

    /// Number of joints.
    pub fn len(&self) -> usize {
        self.joints.len()
    }

    /// Whether the joint set is empty.
    pub fn is_empty(&self) -> bool {
        self.joints.is_empty()
    }

    fn __repr__(&self) -> String {
        format!("PyJointSet(len={})", self.joints.len())
    }
}

// ---------------------------------------------------------------------------
// PyPipelineConfig
// ---------------------------------------------------------------------------

/// Configuration for the physics pipeline.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyPipelineConfig {
    /// Gravity vector `[gx, gy, gz]`.
    #[pyo3(get, set)]
    pub gravity: [f64; 3],
    /// Maximum number of substeps per `step()` call.
    #[pyo3(get, set)]
    pub max_substeps: u32,
    /// Enable sleeping.
    #[pyo3(get, set)]
    pub sleep_enabled: bool,
    /// Linear sleep threshold.
    #[pyo3(get, set)]
    pub linear_sleep_threshold: f64,
    /// Angular sleep threshold.
    #[pyo3(get, set)]
    pub angular_sleep_threshold: f64,
    /// Baumgarte stabilisation factor.
    #[pyo3(get, set)]
    pub baumgarte_factor: f64,
}

#[pymethods]
impl PyPipelineConfig {
    #[new]
    #[pyo3(signature = (gravity = [0.0, -9.81, 0.0]))]
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            gravity,
            max_substeps: 1,
            sleep_enabled: true,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            baumgarte_factor: 0.2,
        }
    }

    fn __repr__(&self) -> String {
        format!("PyPipelineConfig(gravity={:?})", self.gravity)
    }
}

impl Default for PyPipelineConfig {
    fn default() -> Self {
        Self::new([0.0, -9.81, 0.0])
    }
}

// ---------------------------------------------------------------------------
// PyWorld
// ---------------------------------------------------------------------------

/// A self-contained physics world combining rigid bodies and joints.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyWorld {
    /// Body storage.
    #[pyo3(get)]
    pub bodies: PyRigidBodySet,
    /// Joint storage.
    #[pyo3(get)]
    pub joints: PyJointSet,
    /// Pipeline configuration.
    #[pyo3(get)]
    pub config: PyPipelineConfig,
    /// Accumulated simulation time.
    #[pyo3(get)]
    pub time: f64,
    /// Contacts from the most recent step.
    contacts: Vec<PyContact>,
}

#[pymethods]
impl PyWorld {
    #[new]
    #[pyo3(signature = (gravity = [0.0, -9.81, 0.0]))]
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            bodies: PyRigidBodySet::new(),
            joints: PyJointSet::new(),
            config: PyPipelineConfig::new(gravity),
            time: 0.0,
            contacts: Vec::new(),
        }
    }

    /// Get contacts from the most recent step.
    pub fn get_contacts(&self) -> Vec<PyContact> {
        self.contacts.clone()
    }

    /// Number of bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Number of joints.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }

    /// Accumulated simulation time in seconds.
    pub fn elapsed_time(&self) -> f64 {
        self.time
    }

    fn __repr__(&self) -> String {
        format!(
            "PyWorld(bodies={}, joints={}, time={:.3})",
            self.bodies.len(),
            self.joints.len(),
            self.time
        )
    }
}

// ---------------------------------------------------------------------------
// PyPhysicsPipeline
// ---------------------------------------------------------------------------

/// A stateless step executor for `PyWorld`.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct PyPhysicsPipeline {
    /// Number of steps executed.
    #[pyo3(get)]
    pub step_count: u64,
}

#[pymethods]
impl PyPhysicsPipeline {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Step the world forward by `dt` seconds using `substeps` substeps.
    ///
    /// Applies gravity, integrates velocities and positions, and resolves joints
    /// using a simple soft-constraint correction.
    pub fn step(&mut self, world: &mut PyWorld, dt: f64, substeps: u32) {
        let substeps = substeps.max(1);
        let sub_dt = dt / substeps as f64;

        for _ in 0..substeps {
            self.step_once(world, sub_dt);
        }
    }

    /// Return total number of steps executed.
    pub fn total_steps(&self) -> u64 {
        self.step_count
    }

    fn __repr__(&self) -> String {
        format!("PyPhysicsPipeline(steps={})", self.step_count)
    }
}

impl PyPhysicsPipeline {
    fn step_once(&mut self, world: &mut PyWorld, dt: f64) {
        let g = world.config.gravity;

        // Integrate all dynamic bodies
        for body in world.bodies.bodies.iter_mut() {
            if body.is_static || body.is_kinematic || body.mass <= 0.0 {
                continue;
            }
            let inv_m = 1.0 / body.mass;

            // Gravity
            body.velocity[0] += g[0] * dt;
            body.velocity[1] += g[1] * dt;
            body.velocity[2] += g[2] * dt;

            // Linear damping
            let ld = (1.0_f64 - body.linear_damping * dt).max(0.0_f64);
            body.velocity[0] *= ld;
            body.velocity[1] *= ld;
            body.velocity[2] *= ld;

            // Integrate position
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;

            let _ = inv_m; // used implicitly above
        }

        // Soft joint resolution
        let joints: Vec<PyJoint> = world.joints.joints.clone();
        for joint in &joints {
            if !joint.enabled {
                continue;
            }
            let pa = world
                .bodies
                .bodies
                .iter()
                .find(|b| b.handle == joint.body_a)
                .map(|b| b.position);
            let pb = world
                .bodies
                .bodies
                .iter()
                .find(|b| b.handle == joint.body_b)
                .map(|b| b.position);

            let (pa, pb) = match (pa, pb) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };

            let diff = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
            let dist = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();

            if dist < 1e-12 {
                continue;
            }

            let target = joint.target_distance;
            let error = dist - target;
            if error.abs() < 1e-8 {
                continue;
            }

            let n = [diff[0] / dist, diff[1] / dist, diff[2] / dist];
            let correction = error * joint.stiffness * 0.5;

            if let Some(ba) = world
                .bodies
                .bodies
                .iter_mut()
                .find(|b| b.handle == joint.body_a && !b.is_static)
            {
                ba.position[0] += n[0] * correction;
                ba.position[1] += n[1] * correction;
                ba.position[2] += n[2] * correction;
            }
            if let Some(bb) = world
                .bodies
                .bodies
                .iter_mut()
                .find(|b| b.handle == joint.body_b && !b.is_static)
            {
                bb.position[0] -= n[0] * correction;
                bb.position[1] -= n[1] * correction;
                bb.position[2] -= n[2] * correction;
            }
        }

        world.time += dt;
        self.step_count += 1;
    }
}
