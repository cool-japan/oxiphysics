// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python-facing `#[pyclass]` wrapper types for the OxiPhysics engine.
//!
//! These types form the public Python API surface and are registered in the
//! `oxiphysics` extension module defined in `lib.rs`.
//!
//! ## Quick start
//!
//! ```python
//! import oxiphysics
//!
//! world = oxiphysics.PhysicsWorld(gy=-9.81)
//!
//! cfg = oxiphysics.RigidBodyConfig.dynamic(mass=1.0, y=10.0)
//! cfg.add_sphere(radius=0.5)
//! ball = world.add_rigid_body(cfg)
//!
//! for _ in range(60):
//!     world.step(1.0 / 60.0)
//!
//! x, y, z = world.get_position(ball)
//! print(f"y = {y:.3f}")
//!
//! # NumPy
//! import numpy as np
//! pos = np.array(world.all_positions_flat()).reshape(-1, 3)
//! ```

use pyo3::prelude::*;

use crate::types::{PyColliderShape, PyContactResult, PyRigidBodyConfig, PySimConfig};
use crate::world_api::PyPhysicsWorld;

// ============================================================================
// SimConfig
// ============================================================================

/// Simulation-wide configuration (gravity, solver settings, sleep thresholds).
///
/// Create via static factory methods or the constructor:
///
/// ```python
/// cfg = oxiphysics.SimConfig.earth_gravity()
/// cfg.solver_iterations = 16
/// ```
#[pyclass(name = "SimConfig", skip_from_py_object)]
#[derive(Clone)]
pub struct SimConfig {
    pub(crate) inner: PySimConfig,
}

#[pymethods]
impl SimConfig {
    /// Create a `SimConfig` with the given gravity components.
    #[new]
    #[pyo3(signature = (gx=0.0, gy=-9.81, gz=0.0))]
    fn new(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            inner: PySimConfig {
                gravity: [gx, gy, gz],
                ..PySimConfig::default()
            },
        }
    }

    /// Standard Earth gravity (`[0, -9.81, 0]`).
    #[staticmethod]
    fn earth_gravity() -> Self {
        Self {
            inner: PySimConfig::earth_gravity(),
        }
    }

    /// Zero gravity (space / micro-gravity).
    #[staticmethod]
    fn zero_gravity() -> Self {
        Self {
            inner: PySimConfig::zero_gravity(),
        }
    }

    /// Moon gravity (~1/6 of Earth, `[0, -1.62, 0]`).
    #[staticmethod]
    fn moon_gravity() -> Self {
        Self {
            inner: PySimConfig::moon_gravity(),
        }
    }

    /// Gravity vector as `(gx, gy, gz)`.
    #[getter]
    fn gravity(&self) -> (f64, f64, f64) {
        let g = self.inner.gravity;
        (g[0], g[1], g[2])
    }

    /// Set gravity from a `(gx, gy, gz)` tuple or list.
    #[setter]
    fn set_gravity(&mut self, g: (f64, f64, f64)) {
        self.inner.gravity = [g.0, g.1, g.2];
    }

    /// Constraint-solver iterations per step (default 8, increase for stability).
    #[getter]
    fn solver_iterations(&self) -> u32 {
        self.inner.solver_iterations
    }

    #[setter]
    fn set_solver_iterations(&mut self, v: u32) {
        self.inner.solver_iterations = v;
    }

    /// Whether sleeping is globally enabled.
    #[getter]
    fn sleep_enabled(&self) -> bool {
        self.inner.sleep_enabled
    }

    #[setter]
    fn set_sleep_enabled(&mut self, v: bool) {
        self.inner.sleep_enabled = v;
    }

    /// Linear velocity threshold for sleep detection.
    #[getter]
    fn linear_sleep_threshold(&self) -> f64 {
        self.inner.linear_sleep_threshold
    }

    #[setter]
    fn set_linear_sleep_threshold(&mut self, v: f64) {
        self.inner.linear_sleep_threshold = v;
    }

    /// Angular velocity threshold for sleep detection.
    #[getter]
    fn angular_sleep_threshold(&self) -> f64 {
        self.inner.angular_sleep_threshold
    }

    #[setter]
    fn set_angular_sleep_threshold(&mut self, v: f64) {
        self.inner.angular_sleep_threshold = v;
    }

    /// Whether continuous collision detection (CCD) is enabled.
    #[getter]
    fn ccd_enabled(&self) -> bool {
        self.inner.ccd_enabled
    }

    #[setter]
    fn set_ccd_enabled(&mut self, v: bool) {
        self.inner.ccd_enabled = v;
    }

    /// Baumgarte stabilisation factor (0–1, typical 0.1–0.3).
    #[getter]
    fn baumgarte_factor(&self) -> f64 {
        self.inner.baumgarte_factor
    }

    #[setter]
    fn set_baumgarte_factor(&mut self, v: f64) {
        self.inner.baumgarte_factor = v;
    }

    fn __repr__(&self) -> String {
        let g = self.inner.gravity;
        format!(
            "SimConfig(gravity=({:.3}, {:.3}, {:.3}), solver_iterations={})",
            g[0], g[1], g[2], self.inner.solver_iterations,
        )
    }
}

// ============================================================================
// RigidBodyConfig
// ============================================================================

/// Configuration for creating a rigid body.
///
/// ```python
/// # Dynamic sphere
/// cfg = oxiphysics.RigidBodyConfig.dynamic(mass=1.0, x=0.0, y=10.0, z=0.0)
/// cfg.add_sphere(radius=0.5)
/// cfg.set_restitution(0.7)
///
/// # Static ground plane
/// ground = oxiphysics.RigidBodyConfig.static_body(y=-0.5)
/// ground.add_box(hx=50.0, hy=0.5, hz=50.0)
/// ```
#[pyclass(name = "RigidBodyConfig", skip_from_py_object)]
#[derive(Clone)]
pub struct RigidBodyConfig {
    pub(crate) inner: PyRigidBodyConfig,
}

#[pymethods]
impl RigidBodyConfig {
    /// Create a dynamic (physics-simulated) body at `(x, y, z)`.
    #[staticmethod]
    #[pyo3(signature = (mass, x=0.0, y=0.0, z=0.0))]
    fn dynamic(mass: f64, x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: PyRigidBodyConfig::dynamic(mass, [x, y, z]),
        }
    }

    /// Create a static (immovable) body at `(x, y, z)`.
    #[staticmethod]
    #[pyo3(signature = (x=0.0, y=0.0, z=0.0))]
    fn static_body(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: PyRigidBodyConfig::static_body([x, y, z]),
        }
    }

    /// Attach a sphere collider with the given `radius`.
    fn add_sphere(&mut self, radius: f64) -> PyResult<()> {
        if radius <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "radius must be positive",
            ));
        }
        self.inner.shapes.push(PyColliderShape::sphere(radius));
        Ok(())
    }

    /// Attach an axis-aligned box collider (`hx`, `hy`, `hz` are half-extents).
    fn add_box(&mut self, hx: f64, hy: f64, hz: f64) -> PyResult<()> {
        if hx <= 0.0 || hy <= 0.0 || hz <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "half-extents must be positive",
            ));
        }
        self.inner
            .shapes
            .push(PyColliderShape::box_shape(hx, hy, hz));
        Ok(())
    }

    /// Attach a capsule collider (cylinder + two hemispheres, Y-axis aligned).
    fn add_capsule(&mut self, radius: f64, half_height: f64) -> PyResult<()> {
        if radius <= 0.0 || half_height <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "radius and half_height must be positive",
            ));
        }
        self.inner
            .shapes
            .push(PyColliderShape::capsule(radius, half_height));
        Ok(())
    }

    /// Set friction coefficient (0 = frictionless, 1 = high friction).
    fn set_friction(&mut self, friction: f64) {
        self.inner.friction = friction;
    }

    /// Set restitution / bounciness (0 = inelastic, 1 = perfectly elastic).
    fn set_restitution(&mut self, restitution: f64) {
        self.inner.restitution = restitution;
    }

    /// Set linear damping (air drag on translation).
    fn set_linear_damping(&mut self, damping: f64) {
        self.inner.linear_damping = damping;
    }

    /// Set angular damping (air drag on rotation).
    fn set_angular_damping(&mut self, damping: f64) {
        self.inner.angular_damping = damping;
    }

    /// Attach a string tag for identification.
    fn set_tag(&mut self, tag: String) {
        self.inner.tag = Some(tag);
    }

    /// Enable or disable sleep for this body.
    fn set_can_sleep(&mut self, can_sleep: bool) {
        self.inner.can_sleep = can_sleep;
    }

    /// Set the initial linear velocity `(vx, vy, vz)`.
    fn set_velocity(&mut self, vx: f64, vy: f64, vz: f64) {
        self.inner.velocity = [vx, vy, vz];
    }

    /// Set the initial orientation as a unit quaternion `(qx, qy, qz, qw)`.
    fn set_orientation(&mut self, qx: f64, qy: f64, qz: f64, qw: f64) {
        self.inner.orientation = [qx, qy, qz, qw];
    }

    /// Body mass in kilograms.
    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass
    }

    /// Initial position `(x, y, z)`.
    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        let p = self.inner.position;
        (p[0], p[1], p[2])
    }

    /// Friction coefficient.
    #[getter]
    fn friction(&self) -> f64 {
        self.inner.friction
    }

    /// Restitution coefficient.
    #[getter]
    fn restitution(&self) -> f64 {
        self.inner.restitution
    }

    /// `True` if this is a static (immovable) body.
    #[getter]
    fn is_static(&self) -> bool {
        self.inner.is_static
    }

    /// Number of collider shapes attached.
    #[getter]
    fn shape_count(&self) -> usize {
        self.inner.shapes.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "RigidBodyConfig(mass={}, position={:?}, shapes={}, is_static={})",
            self.inner.mass,
            self.inner.position,
            self.inner.shapes.len(),
            self.inner.is_static,
        )
    }
}

// ============================================================================
// ContactResult
// ============================================================================

/// A contact / collision event from the most recent simulation step.
///
/// Retrieve via `PhysicsWorld.get_contacts()`.
#[pyclass(name = "ContactResult", skip_from_py_object)]
#[derive(Clone)]
pub struct ContactResult {
    pub(crate) inner: PyContactResult,
}

#[pymethods]
impl ContactResult {
    /// Handle of the first body involved.
    #[getter]
    fn body_a(&self) -> u32 {
        self.inner.body_a
    }

    /// Handle of the second body involved.
    #[getter]
    fn body_b(&self) -> u32 {
        self.inner.body_b
    }

    /// Contact point in world space `(x, y, z)`.
    #[getter]
    fn contact_point(&self) -> (f64, f64, f64) {
        let p = self.inner.contact_point;
        (p[0], p[1], p[2])
    }

    /// Contact normal from body_a toward body_b `(nx, ny, nz)`.
    #[getter]
    fn normal(&self) -> (f64, f64, f64) {
        let n = self.inner.normal;
        (n[0], n[1], n[2])
    }

    /// Penetration depth (positive = overlapping).
    #[getter]
    fn depth(&self) -> f64 {
        self.inner.depth
    }

    /// Impulse applied to resolve this contact.
    #[getter]
    fn impulse(&self) -> f64 {
        self.inner.impulse
    }

    /// Combined friction coefficient.
    #[getter]
    fn friction(&self) -> f64 {
        self.inner.friction
    }

    /// Combined restitution coefficient.
    #[getter]
    fn restitution(&self) -> f64 {
        self.inner.restitution
    }

    /// `True` if this contact represents genuine penetration (depth > 0).
    fn is_colliding(&self) -> bool {
        self.inner.is_colliding()
    }

    fn __repr__(&self) -> String {
        format!(
            "ContactResult(bodies=({}, {}), depth={:.4}, impulse={:.4})",
            self.inner.body_a, self.inner.body_b, self.inner.depth, self.inner.impulse,
        )
    }
}

// ============================================================================
// PhysicsWorld
// ============================================================================

/// Self-contained physics simulation world.
///
/// # Quick start
///
/// ```python
/// import oxiphysics
///
/// world = oxiphysics.PhysicsWorld(gy=-9.81)
///
/// # Add a sphere that falls
/// cfg = oxiphysics.RigidBodyConfig.dynamic(mass=1.0, y=10.0)
/// cfg.add_sphere(radius=0.5)
/// ball = world.add_rigid_body(cfg)
///
/// # Add a static ground box
/// gcfg = oxiphysics.RigidBodyConfig.static_body()
/// gcfg.add_box(hx=50.0, hy=0.1, hz=50.0)
/// world.add_rigid_body(gcfg)
///
/// # Simulate 1 s at 60 Hz
/// for _ in range(60):
///     world.step(1.0 / 60.0)
///
/// x, y, z = world.get_position(ball)
/// print(f"Ball at y={y:.3f}")
///
/// # NumPy bulk access
/// import numpy as np
/// positions = np.array(world.all_positions_flat()).reshape(-1, 3)
/// ```
#[pyclass(name = "PhysicsWorld")]
pub struct PhysicsWorld {
    inner: PyPhysicsWorld,
}

#[pymethods]
impl PhysicsWorld {
    /// Create a world with the given gravity components (default: Earth gravity).
    #[new]
    #[pyo3(signature = (gx=0.0, gy=-9.81, gz=0.0))]
    fn new(gx: f64, gy: f64, gz: f64) -> Self {
        Self {
            inner: PyPhysicsWorld::new_with_config(PySimConfig {
                gravity: [gx, gy, gz],
                ..PySimConfig::default()
            }),
        }
    }

    /// Create a world with Earth gravity `(0, -9.81, 0)`.
    #[staticmethod]
    fn with_earth_gravity() -> Self {
        Self {
            inner: PyPhysicsWorld::new_with_config(PySimConfig::earth_gravity()),
        }
    }

    /// Create a zero-gravity world.
    #[staticmethod]
    fn with_zero_gravity() -> Self {
        Self {
            inner: PyPhysicsWorld::new_with_config(PySimConfig::zero_gravity()),
        }
    }

    /// Create a world from a `SimConfig` object.
    #[staticmethod]
    fn from_config(config: PyRef<'_, SimConfig>) -> Self {
        Self {
            inner: PyPhysicsWorld::new_with_config(config.inner.clone()),
        }
    }

    // -----------------------------------------------------------------------
    // Body management
    // -----------------------------------------------------------------------

    /// Add a body described by `config` and return its integer handle.
    fn add_rigid_body(&mut self, config: PyRef<'_, RigidBodyConfig>) -> u32 {
        self.inner.add_rigid_body(config.inner.clone())
    }

    /// Shorthand: add a dynamic sphere at `(x, y, z)`.
    fn add_sphere_body(&mut self, mass: f64, x: f64, y: f64, z: f64, radius: f64) -> u32 {
        let cfg =
            PyRigidBodyConfig::dynamic(mass, [x, y, z]).with_shape(PyColliderShape::sphere(radius));
        self.inner.add_rigid_body(cfg)
    }

    /// Shorthand: add a dynamic box body (half-extents `hx, hy, hz`).
    fn add_box_body(
        &mut self,
        mass: f64,
        x: f64,
        y: f64,
        z: f64,
        hx: f64,
        hy: f64,
        hz: f64,
    ) -> u32 {
        let cfg = PyRigidBodyConfig::dynamic(mass, [x, y, z])
            .with_shape(PyColliderShape::box_shape(hx, hy, hz));
        self.inner.add_rigid_body(cfg)
    }

    /// Shorthand: add a dynamic capsule body.
    fn add_capsule_body(
        &mut self,
        mass: f64,
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        half_height: f64,
    ) -> u32 {
        let cfg = PyRigidBodyConfig::dynamic(mass, [x, y, z])
            .with_shape(PyColliderShape::capsule(radius, half_height));
        self.inner.add_rigid_body(cfg)
    }

    /// Shorthand: add a static body at `(x, y, z)` (no colliders attached).
    fn add_static_body(&mut self, x: f64, y: f64, z: f64) -> u32 {
        self.inner
            .add_rigid_body(PyRigidBodyConfig::static_body([x, y, z]))
    }

    /// Remove a body by handle. Returns `True` if it existed.
    fn remove_body(&mut self, handle: u32) -> bool {
        self.inner.remove_body(handle)
    }

    // -----------------------------------------------------------------------
    // Simulation control
    // -----------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds (one step).
    ///
    /// ## asyncio integration
    ///
    /// To step the simulation without blocking an asyncio event loop, run it in
    /// a thread-pool via `asyncio.to_thread` (Python 3.9+) or
    /// `loop.run_in_executor`:
    ///
    /// ```python
    /// import asyncio, oxiphysics
    ///
    /// async def simulate(world: oxiphysics.PhysicsWorld, frames: int = 60):
    ///     for _ in range(frames):
    ///         await asyncio.to_thread(world.step, 1.0 / 60.0)
    ///
    /// # Or with an explicit executor:
    /// async def simulate_executor(world, frames=60):
    ///     loop = asyncio.get_running_loop()
    ///     import concurrent.futures
    ///     pool = concurrent.futures.ThreadPoolExecutor()
    ///     for _ in range(frames):
    ///         await loop.run_in_executor(pool, world.step, 1.0 / 60.0)
    /// ```
    fn step(&mut self, dt: f64) {
        self.inner.step(dt);
    }

    /// Advance the simulation by `dt` seconds using `substeps` equal sub-steps.
    ///
    /// Useful for fast-moving objects that would otherwise tunnel through geometry.
    ///
    /// ## asyncio integration
    ///
    /// ```python
    /// import asyncio, oxiphysics
    ///
    /// async def tick(world: oxiphysics.PhysicsWorld):
    ///     await asyncio.to_thread(world.step_substeps, 1.0 / 60.0, 4)
    /// ```
    fn step_substeps(&mut self, dt: f64, substeps: u32) {
        self.inner.step_substeps(dt, substeps);
    }

    /// Remove all bodies, contacts and constraints; reset simulation time to 0.
    fn reset(&mut self) {
        self.inner.reset();
    }

    // -----------------------------------------------------------------------
    // Body state – queries
    // -----------------------------------------------------------------------

    /// Position `(x, y, z)` of `handle`, or `None` if invalid.
    fn get_position(&self, handle: u32) -> Option<(f64, f64, f64)> {
        self.inner.get_position(handle).map(|p| (p[0], p[1], p[2]))
    }

    /// Linear velocity `(vx, vy, vz)`, or `None`.
    fn get_velocity(&self, handle: u32) -> Option<(f64, f64, f64)> {
        self.inner.get_velocity(handle).map(|v| (v[0], v[1], v[2]))
    }

    /// Orientation quaternion `(qx, qy, qz, qw)`, or `None`.
    fn get_orientation(&self, handle: u32) -> Option<(f64, f64, f64, f64)> {
        self.inner
            .get_orientation(handle)
            .map(|q| (q[0], q[1], q[2], q[3]))
    }

    /// Angular velocity `(wx, wy, wz)` in rad/s, or `None`.
    fn get_angular_velocity(&self, handle: u32) -> Option<(f64, f64, f64)> {
        self.inner
            .get_angular_velocity(handle)
            .map(|v| (v[0], v[1], v[2]))
    }

    /// Optional string tag of `handle`, or `None`.
    fn get_tag(&self, handle: u32) -> Option<String> {
        self.inner.get_tag(handle)
    }

    // -----------------------------------------------------------------------
    // Body state – mutations
    // -----------------------------------------------------------------------

    /// Teleport the body to `(x, y, z)` (bypasses physics).
    fn set_position(&mut self, handle: u32, x: f64, y: f64, z: f64) {
        self.inner.set_position(handle, [x, y, z]);
    }

    /// Set linear velocity to `(vx, vy, vz)`.
    fn set_velocity(&mut self, handle: u32, vx: f64, vy: f64, vz: f64) {
        self.inner.set_velocity(handle, [vx, vy, vz]);
    }

    /// Set orientation from quaternion `(qx, qy, qz, qw)`.
    fn set_orientation(&mut self, handle: u32, qx: f64, qy: f64, qz: f64, qw: f64) {
        self.inner.set_orientation(handle, [qx, qy, qz, qw]);
    }

    /// Set angular velocity `(wx, wy, wz)` in rad/s.
    fn set_angular_velocity(&mut self, handle: u32, wx: f64, wy: f64, wz: f64) {
        self.inner.set_angular_velocity(handle, [wx, wy, wz]);
    }

    // -----------------------------------------------------------------------
    // Force / impulse
    // -----------------------------------------------------------------------

    /// Apply a force `(fx, fy, fz)` at the body's centre of mass.
    ///
    /// The force accumulates until the next `step()`.
    fn apply_force(&mut self, handle: u32, fx: f64, fy: f64, fz: f64) {
        self.inner.apply_force(handle, [fx, fy, fz], None);
    }

    /// Apply a force `(fx, fy, fz)` at world-space point `(px, py, pz)`.
    fn apply_force_at(
        &mut self,
        handle: u32,
        fx: f64,
        fy: f64,
        fz: f64,
        px: f64,
        py: f64,
        pz: f64,
    ) {
        self.inner
            .apply_force(handle, [fx, fy, fz], Some([px, py, pz]));
    }

    /// Apply an instantaneous linear impulse `(ix, iy, iz)`.
    fn apply_impulse(&mut self, handle: u32, ix: f64, iy: f64, iz: f64) {
        self.inner.apply_impulse(handle, [ix, iy, iz], None);
    }

    /// Apply a torque `(tx, ty, tz)` (accumulates until next step).
    fn apply_torque(&mut self, handle: u32, tx: f64, ty: f64, tz: f64) {
        self.inner.apply_torque(handle, [tx, ty, tz]);
    }

    // -----------------------------------------------------------------------
    // Sleep
    // -----------------------------------------------------------------------

    /// `True` if `handle` is currently sleeping.
    fn is_sleeping(&self, handle: u32) -> bool {
        self.inner.is_sleeping(handle)
    }

    /// Wake a sleeping body.
    fn wake_body(&mut self, handle: u32) {
        self.inner.wake_body(handle);
    }

    /// Force a body into sleep.
    fn sleep_body(&mut self, handle: u32) {
        self.inner.sleep_body(handle);
    }

    // -----------------------------------------------------------------------
    // Global configuration
    // -----------------------------------------------------------------------

    /// Set the global gravity vector.
    fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.inner.set_gravity([gx, gy, gz]);
    }

    /// Current gravity as `(gx, gy, gz)`.
    fn gravity(&self) -> (f64, f64, f64) {
        let g = self.inner.gravity();
        (g[0], g[1], g[2])
    }

    // -----------------------------------------------------------------------
    // Statistics
    // -----------------------------------------------------------------------

    /// Number of active bodies.
    fn body_count(&self) -> usize {
        self.inner.body_count()
    }

    /// Number of bodies currently sleeping.
    fn sleeping_count(&self) -> usize {
        self.inner.sleeping_count()
    }

    /// Number of contacts detected in the most recent step.
    fn contact_count(&self) -> usize {
        self.inner.contact_count()
    }

    /// Total accumulated simulation time in seconds.
    fn time(&self) -> f64 {
        self.inner.time()
    }

    // -----------------------------------------------------------------------
    // Bulk queries
    // -----------------------------------------------------------------------

    /// Handles of all active bodies.
    fn active_handles(&self) -> Vec<u32> {
        self.inner.active_handles()
    }

    /// All body positions as a flat list `[x0, y0, z0, x1, y1, z1, ...]`.
    ///
    /// Ideal for NumPy:
    ///
    /// ```python
    /// import numpy as np
    /// pos = np.array(world.all_positions_flat()).reshape(-1, 3)
    /// ```
    fn all_positions_flat(&self) -> Vec<f64> {
        self.inner
            .all_positions()
            .into_iter()
            .flat_map(|p: [f64; 3]| p)
            .collect()
    }

    /// All body velocities as a flat list `[vx0, vy0, vz0, vx1, vy1, vz1, ...]`.
    ///
    /// ```python
    /// import numpy as np
    /// vel = np.array(world.all_velocities_flat()).reshape(-1, 3)
    /// ```
    fn all_velocities_flat(&self) -> Vec<f64> {
        self.inner
            .all_velocities()
            .into_iter()
            .flat_map(|v: [f64; 3]| v)
            .collect()
    }

    /// Contacts from the most recent step as a list of `ContactResult` objects.
    fn get_contacts(&self) -> Vec<ContactResult> {
        self.inner
            .get_contacts()
            .into_iter()
            .map(|c| ContactResult { inner: c })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Spatial queries
    // -----------------------------------------------------------------------

    /// Handles of bodies whose centre lies within the given AABB.
    fn bodies_in_aabb(
        &self,
        xmin: f64,
        ymin: f64,
        zmin: f64,
        xmax: f64,
        ymax: f64,
        zmax: f64,
    ) -> Vec<u32> {
        self.inner
            .bodies_in_aabb([xmin, ymin, zmin, xmax, ymax, zmax])
    }

    /// Cast a ray from `(ox, oy, oz)` in direction `(dx, dy, dz)`.
    ///
    /// Returns `(handle, distance)` of the nearest sphere hit, or `None`.
    fn raycast(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_dist: f64,
    ) -> Option<(u32, f64)> {
        self.inner.raycast([ox, oy, oz], [dx, dy, dz], max_dist)
    }

    fn __repr__(&self) -> String {
        let g = self.inner.gravity();
        format!(
            "PhysicsWorld(bodies={}, time={:.3}, gravity=({:.2}, {:.2}, {:.2}))",
            self.inner.body_count(),
            self.inner.time(),
            g[0],
            g[1],
            g[2],
        )
    }
}
