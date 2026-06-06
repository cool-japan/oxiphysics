// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Energy, Gravity Field, Contact List, and Inertia Extensions.

use super::PyPhysicsWorld;
use pyo3::prelude::*;

// ===========================================================================
// Energy, Gravity Field, Contact List, and Inertia Extensions
// ===========================================================================

/// Contact pair returned by `get_contact_list`.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, PartialEq)]
pub struct ContactPair {
    /// Handle of body A.
    #[pyo3(get, set)]
    pub body_a: u32,
    /// Handle of body B.
    #[pyo3(get, set)]
    pub body_b: u32,
    /// World-space contact point.
    #[pyo3(get, set)]
    pub contact_point: [f64; 3],
    /// Contact normal (from B toward A).
    #[pyo3(get, set)]
    pub normal: [f64; 3],
    /// Penetration depth.
    #[pyo3(get, set)]
    pub depth: f64,
    /// Impulse magnitude applied.
    #[pyo3(get, set)]
    pub impulse: f64,
}

/// Inertia tensor (3×3 matrix stored as 9 elements, row-major).
#[pyclass(from_py_object)]
#[derive(Debug, Clone, PartialEq)]
pub struct InertiaTensor {
    /// Row-major 3×3 inertia tensor elements.
    #[pyo3(get, set)]
    pub elements: [f64; 9],
}

#[pymethods]
impl InertiaTensor {
    /// Create from diagonal elements (assumes principal axes aligned).
    #[staticmethod]
    pub fn from_diagonal(ix: f64, iy: f64, iz: f64) -> Self {
        let elements = [ix, 0.0, 0.0, 0.0, iy, 0.0, 0.0, 0.0, iz];
        Self { elements }
    }

    /// Return diagonal elements `[Ixx, Iyy, Izz]`.
    pub fn diagonal(&self) -> [f64; 3] {
        [self.elements[0], self.elements[4], self.elements[8]]
    }

    /// Trace (sum of diagonal).
    pub fn trace(&self) -> f64 {
        self.elements[0] + self.elements[4] + self.elements[8]
    }
}

/// `PyRigidBody` is a lightweight handle-based wrapper that exposes
/// body-level computations (moment of inertia, etc.) without owning state.
#[pyclass]
pub struct PyRigidBody {
    /// Mass of the body (kg).
    #[pyo3(get, set)]
    pub mass: f64,
    /// Half-extents `[hx, hy, hz]` for box-shaped bodies.
    #[pyo3(get, set)]
    pub half_extents: [f64; 3],
}

#[pymethods]
impl PyRigidBody {
    /// Create a new `PyRigidBody` with the given mass and box half-extents.
    #[new]
    pub fn new(mass: f64, half_extents: [f64; 3]) -> Self {
        Self { mass, half_extents }
    }

    /// Compute the solid box inertia tensor for this body.
    ///
    /// For a solid box with mass `m` and half-extents `(hx, hy, hz)`:
    /// - `Ixx = m/3 * (hy² + hz²)`
    /// - `Iyy = m/3 * (hx² + hz²)`
    /// - `Izz = m/3 * (hx² + hy²)`
    ///
    /// Returns an `InertiaTensor` with the principal-axis diagonal.
    pub fn compute_moment_of_inertia_box(&self) -> InertiaTensor {
        let m = self.mass;
        let hx = self.half_extents[0];
        let hy = self.half_extents[1];
        let hz = self.half_extents[2];
        // Full-extent dimensions: lx=2*hx, ly=2*hy, lz=2*hz
        // I_xx = m/12 * (ly² + lz²) = m/12 * (4hy² + 4hz²) = m/3 * (hy² + hz²)
        let ixx = (m / 3.0) * (hy * hy + hz * hz);
        let iyy = (m / 3.0) * (hx * hx + hz * hz);
        let izz = (m / 3.0) * (hx * hx + hy * hy);
        InertiaTensor::from_diagonal(ixx, iyy, izz)
    }

    /// Compute moment of inertia for a solid sphere.
    ///
    /// `I = 2/5 * m * r²` for each axis (isotropic).
    pub fn compute_moment_of_inertia_sphere(&self, radius: f64) -> InertiaTensor {
        let i = 0.4 * self.mass * radius * radius;
        InertiaTensor::from_diagonal(i, i, i)
    }
}

impl PyPhysicsWorld {
    // -----------------------------------------------------------------------
    // Total energy
    // -----------------------------------------------------------------------

    /// Compute total mechanical energy: kinetic energy + gravitational potential energy.
    ///
    /// KE = Σ ½ m v²  (summed over dynamic bodies)
    /// PE = Σ m * |g| * h  (height measured along gravity direction; uses Y component of gravity)
    ///
    /// Returns `(kinetic, potential, total)`.
    pub fn compute_total_energy(&self) -> (f64, f64, f64) {
        let g = self.config.gravity;
        let g_mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        // Gravity direction unit vector (pointing downward)
        let g_dir = if g_mag > 1e-10 {
            [g[0] / g_mag, g[1] / g_mag, g[2] / g_mag]
        } else {
            [0.0, -1.0, 0.0]
        };

        let mut ke = 0.0_f64;
        let mut pe = 0.0_f64;

        for slot in &self.slots {
            let body = match slot.body.as_ref() {
                Some(b) if !b.is_static && !b.is_kinematic => b,
                _ => continue,
            };

            let v = body.velocity;
            let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            ke += 0.5 * body.mass * v2;

            // Height above origin along negative-gravity direction
            // h = -dot(position, g_dir)  (positive when above ground)
            let h = -(body.position[0] * g_dir[0]
                + body.position[1] * g_dir[1]
                + body.position[2] * g_dir[2]);
            pe += body.mass * g_mag * h;
        }

        (ke, pe, ke + pe)
    }

    // -----------------------------------------------------------------------
    // Spatially varying gravity
    // -----------------------------------------------------------------------

    /// Apply a spatially varying gravity field to all dynamic bodies for one step.
    ///
    /// The field is specified as a closure `f(position) -> [f64; 3]` mapping
    /// a body's world-space position to a gravitational acceleration vector.
    ///
    /// This method accumulates forces based on the field and does **not**
    /// step the simulation; call `step()` afterwards.
    ///
    /// # Example
    ///
    /// Radial gravity toward the origin:
    /// ```ignore
    /// world.apply_gravity_field(|p| {
    ///     let r2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + 1e-6;
    ///     let r = r2.sqrt();
    ///     [-9.81 * p[0] / r, -9.81 * p[1] / r, -9.81 * p[2] / r]
    /// });
    /// ```
    pub fn apply_gravity_field<F>(&mut self, field: F)
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        for slot in &mut self.slots {
            let body = match slot.body.as_mut() {
                Some(b) if !b.is_static && !b.is_kinematic => b,
                _ => continue,
            };
            let g_local = field(body.position);
            // F = m * g
            body.accumulated_force[0] += body.mass * g_local[0];
            body.accumulated_force[1] += body.mass * g_local[1];
            body.accumulated_force[2] += body.mass * g_local[2];
        }
    }

    // -----------------------------------------------------------------------
    // Contact list
    // -----------------------------------------------------------------------

    /// Return the list of contact pairs from the most recent simulation step.
    ///
    /// Each entry contains the two body handles and contact geometry.
    pub fn get_contact_list(&self) -> Vec<ContactPair> {
        self.contacts
            .iter()
            .map(|c| ContactPair {
                body_a: c.body_a,
                body_b: c.body_b,
                contact_point: c.contact_point,
                normal: c.normal,
                depth: c.depth,
                impulse: c.impulse,
            })
            .collect()
    }
}
