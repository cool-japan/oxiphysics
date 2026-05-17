// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics interpolator module.
//!
//! Exposes smooth damp, exponential decay, spring followers, and lerp utilities.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::interpolator::{SpringFollower, SpringFollower3};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// Free functions
// ─────────────────────────────────────────────────────────────────────────────

/// Linear interpolation between two f64 values.
#[pyfunction]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    oxiphysics::interpolator::lerp(a, b, t)
}

/// Linear interpolation between two Vec3 values.
#[pyfunction]
pub fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    oxiphysics::interpolator::lerp3(a, b, t)
}

/// Exponential decay toward target (scalar). Returns new value.
#[pyfunction]
pub fn exp_decay(current: f64, target: f64, rate: f64, dt: f64) -> f64 {
    oxiphysics::interpolator::exp_decay(current, target, rate, dt)
}

/// Exponential decay toward target (Vec3). Returns new value.
#[pyfunction]
pub fn exp_decay3(current: [f64; 3], target: [f64; 3], rate: f64, dt: f64) -> [f64; 3] {
    oxiphysics::interpolator::exp_decay3(current, target, rate, dt)
}

/// Remap a value from one range to another.
#[pyfunction]
pub fn remap(value: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    oxiphysics::interpolator::remap(value, in_min, in_max, out_min, out_max)
}

/// Cubic Hermite smoothstep.
#[pyfunction]
pub fn smoothstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    oxiphysics::interpolator::smoothstep(edge0, edge1, x)
}

// ─────────────────────────────────────────────────────────────────────────────
// PySpringFollower
// ─────────────────────────────────────────────────────────────────────────────

/// A critically-damped spring follower (scalar).
#[pyclass(name = "SpringFollower")]
pub struct PySpringFollower {
    inner: SpringFollower,
    /// Current target the spring is moving toward.
    target: f64,
}

#[pymethods]
impl PySpringFollower {
    /// Create a spring follower at `initial_value` with `stiffness`.
    #[new]
    pub fn new(initial_value: f64, stiffness: f64) -> Self {
        let mut inner = SpringFollower::new(stiffness);
        inner.set_value(initial_value);
        Self {
            inner,
            target: initial_value,
        }
    }

    /// Create with explicit damping ratio.
    #[staticmethod]
    pub fn with_damping(initial_value: f64, stiffness: f64, damping: f64) -> Self {
        let mut inner = SpringFollower::with_damping(stiffness, damping);
        inner.set_value(initial_value);
        Self {
            inner,
            target: initial_value,
        }
    }

    /// Set the target value.
    pub fn set_target(&mut self, target: f64) {
        self.target = target;
    }

    /// Set the current value directly (teleport).
    pub fn set_value(&mut self, value: f64) {
        self.inner.set_value(value);
    }

    /// Advance by dt and return the new value.
    pub fn step(&mut self, dt: f64) -> f64 {
        self.inner.update(self.target, dt)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySpringFollower3
// ─────────────────────────────────────────────────────────────────────────────

/// A critically-damped spring follower (Vec3).
#[pyclass(name = "SpringFollower3")]
pub struct PySpringFollower3 {
    inner: SpringFollower3,
    /// Current target the spring is moving toward.
    target: [f64; 3],
}

#[pymethods]
impl PySpringFollower3 {
    /// Create a Vec3 spring follower at the given initial value.
    #[new]
    pub fn new(initial_value: [f64; 3], stiffness: f64) -> Self {
        let mut inner = SpringFollower3::new(stiffness);
        inner.set_value(initial_value);
        Self {
            inner,
            target: initial_value,
        }
    }

    /// Set the target Vec3.
    pub fn set_target(&mut self, target: [f64; 3]) {
        self.target = target;
    }

    /// Set the current value directly.
    pub fn set_value(&mut self, value: [f64; 3]) {
        self.inner.set_value(value);
    }

    /// Advance by dt and return the new Vec3 value.
    pub fn step(&mut self, dt: f64) -> [f64; 3] {
        self.inner.update(self.target, dt)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySpringFollower>()?;
    m.add_class::<PySpringFollower3>()?;
    m.add_function(wrap_pyfunction!(lerp, m)?)?;
    m.add_function(wrap_pyfunction!(lerp3, m)?)?;
    m.add_function(wrap_pyfunction!(exp_decay, m)?)?;
    m.add_function(wrap_pyfunction!(exp_decay3, m)?)?;
    m.add_function(wrap_pyfunction!(remap, m)?)?;
    m.add_function(wrap_pyfunction!(smoothstep, m)?)?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lerp() {
        assert!((lerp(0.0, 10.0, 0.5) - 5.0).abs() < 1e-9);
        assert!((lerp(0.0, 10.0, 0.0) - 0.0).abs() < 1e-9);
        assert!((lerp(0.0, 10.0, 1.0) - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_spring_follower() {
        let mut spring = PySpringFollower::new(0.0, 20.0);
        spring.set_target(10.0);
        let v = spring.step(0.1);
        // Should move toward target
        assert!(v > 0.0 && v < 10.0);
    }

    #[test]
    fn test_spring_follower3() {
        let mut spring = PySpringFollower3::new([0.0, 0.0, 0.0], 20.0);
        spring.set_target([10.0, 0.0, 0.0]);
        let v = spring.step(0.1);
        assert!(v[0] > 0.0);
    }
}
