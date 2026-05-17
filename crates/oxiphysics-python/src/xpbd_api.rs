// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics xpbd module.
//!
//! Exposes the Extended Position-Based Dynamics integrator to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::xpbd::{XpbdConstraint, XpbdSolver};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyXpbdSolver
// ─────────────────────────────────────────────────────────────────────────────

/// XPBD solver with particles and constraints.
#[pyclass(name = "XpbdSolver")]
pub struct PyXpbdSolver {
    inner: XpbdSolver,
}

impl Default for PyXpbdSolver {
    fn default() -> Self {
        Self {
            inner: XpbdSolver::new(),
        }
    }
}

#[pymethods]
impl PyXpbdSolver {
    /// Create a new empty XPBD solver with default gravity (0, -9.81, 0).
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set gravity.
    pub fn set_gravity(&mut self, gx: f64, gy: f64, gz: f64) {
        self.inner.gravity = [gx, gy, gz];
    }

    /// Set number of substeps.
    pub fn set_substeps(&mut self, substeps: usize) {
        self.inner.substeps = substeps;
    }

    /// Add a particle at `pos` with inverse mass `inv_mass` (0 = pinned). Returns particle index.
    pub fn add_particle(&mut self, pos: [f64; 3], inv_mass: f64) -> usize {
        self.inner.add_particle(pos, inv_mass)
    }

    /// Pin a particle (set inv_mass = 0).
    pub fn pin_particle(&mut self, idx: usize) {
        self.inner.pin_particle(idx);
    }

    /// Add a distance constraint between particles `a` and `b`.
    pub fn add_distance_constraint(
        &mut self,
        a: usize,
        b: usize,
        rest_length: f64,
        compliance: f64,
    ) -> usize {
        self.inner.add_constraint(XpbdConstraint::Distance {
            a,
            b,
            rest_length,
            compliance,
        })
    }

    /// Add an angle constraint at apex `b` between particles `a` and `c`.
    pub fn add_angle_constraint(
        &mut self,
        a: usize,
        b: usize,
        c: usize,
        rest_angle: f64,
        compliance: f64,
    ) -> usize {
        self.inner.add_constraint(XpbdConstraint::Angle {
            a,
            b,
            c,
            rest_angle,
            compliance,
        })
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.inner.step(dt);
    }

    /// Get the position of particle `idx`, or None if out of bounds.
    pub fn particle_position(&self, idx: usize) -> Option<[f64; 3]> {
        self.inner.particles.get(idx).map(|p| p.pos)
    }

    /// Estimate the velocity of particle `idx` given outer timestep `dt`.
    pub fn particle_velocity(&self, idx: usize, dt: f64) -> [f64; 3] {
        self.inner.particle_velocity(idx, dt)
    }

    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.inner.particles.len()
    }

    /// Number of constraints.
    pub fn constraint_count(&self) -> usize {
        self.inner.constraints.len()
    }

    /// Get all particle positions as a flat list `[x0, y0, z0, x1, y1, z1, ...]`.
    pub fn all_positions(&self) -> Vec<f64> {
        self.inner
            .particles
            .iter()
            .flat_map(|p| [p.pos[0], p.pos[1], p.pos[2]])
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyXpbdSolver>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xpbd_solver_instantiation() {
        let mut solver = PyXpbdSolver::new();
        assert_eq!(solver.particle_count(), 0);
        let i0 = solver.add_particle([0.0, 5.0, 0.0], 1.0);
        let i1 = solver.add_particle([1.0, 5.0, 0.0], 1.0);
        assert_eq!(solver.particle_count(), 2);
        let _cid = solver.add_distance_constraint(i0, i1, 1.0, 0.0);
        assert_eq!(solver.constraint_count(), 1);
    }

    #[test]
    fn test_xpbd_solver_step() {
        let mut solver = PyXpbdSolver::new();
        // Pin a particle at origin, add a free particle above
        let i0 = solver.add_particle([0.0, 0.0, 0.0], 0.0); // pinned
        let i1 = solver.add_particle([0.0, 1.0, 0.0], 1.0); // free
        solver.add_distance_constraint(i0, i1, 1.0, 0.0);
        solver.step(1.0 / 60.0);
        // Position should change but constraint should be maintained
        let p0 = solver.particle_position(i0).expect("p0 should exist");
        let p1 = solver.particle_position(i1).expect("p1 should exist");
        // Pinned particle should not move
        assert!((p0[0]).abs() < 1e-9);
        // Free particle should have moved due to gravity
        assert!(p1[1] < 1.0);
    }
}
