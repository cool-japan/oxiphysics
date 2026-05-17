// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics noise module.
//!
//! Exposes procedural 3D value noise and fractal Brownian motion (fBm)
//! to Python via PyO3.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::noise::{FractalNoise, ValueNoise3D};
use pyo3::prelude::*;

#[cfg(feature = "numpy-bridge")]
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
#[cfg(feature = "numpy-bridge")]
use pyo3::Bound;

// ─────────────────────────────────────────────────────────────────────────────
// PyValueNoise3D
// ─────────────────────────────────────────────────────────────────────────────

/// Single-octave 3D value noise with quintic interpolation.
///
/// Outputs are in `[-1.0, 1.0]` and continuous everywhere.
#[pyclass(name = "ValueNoise3D")]
pub struct PyValueNoise3D {
    inner: ValueNoise3D,
}

#[pymethods]
impl PyValueNoise3D {
    /// Create a new ValueNoise3D with the given seed.
    #[new]
    pub fn new(seed: u32) -> Self {
        Self {
            inner: ValueNoise3D::new(seed),
        }
    }

    /// Sample noise at (x, y, z). Returns a value in [-1.0, 1.0].
    pub fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        self.inner.sample(x, y, z)
    }

    /// Returns the seed used to initialise this noise instance.
    pub fn seed(&self) -> u32 {
        self.inner.seed
    }

    /// Sample a regular grid of noise values and return them as a
    /// `numpy.ndarray` of shape `(nz, ny, nx)` and dtype `float64`.
    ///
    /// # Arguments
    /// * `origin` – `(x0, y0, z0)` world-space position of the `[0,0,0]` cell.
    /// * `step`   – `(dx, dy, dz)` spacing between adjacent cells.
    /// * `shape`  – `(nz, ny, nx)` number of cells along each axis.
    ///
    /// The outermost index is Z so the array is C-contiguous in the
    /// standard (Z, Y, X) layout used by NumPy / SciPy.
    #[cfg(feature = "numpy-bridge")]
    #[pyo3(signature = (origin, step, shape))]
    pub fn sample_grid_to_numpy<'py>(
        &self,
        py: Python<'py>,
        origin: (f64, f64, f64),
        step: (f64, f64, f64),
        shape: (usize, usize, usize),
    ) -> PyResult<Bound<'py, numpy::PyArray3<f64>>> {
        let (nz, ny, nx) = shape;
        let total = nz * ny * nx;
        let mut data = Vec::with_capacity(total);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let x = origin.0 + ix as f64 * step.0;
                    let y = origin.1 + iy as f64 * step.1;
                    let z = origin.2 + iz as f64 * step.2;
                    data.push(self.inner.sample(x, y, z));
                }
            }
        }
        let flat: Bound<'py, PyArray1<f64>> = data.into_pyarray(py);
        flat.reshape([nz, ny, nx])
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyFractalNoise
// ─────────────────────────────────────────────────────────────────────────────

/// Fractal Brownian motion (fBm) noise.
///
/// Sums multiple octaves of value noise with persistence / lacunarity scaling.
#[pyclass(name = "FractalNoise")]
pub struct PyFractalNoise {
    inner: FractalNoise,
}

#[pymethods]
impl PyFractalNoise {
    /// Create fractal noise with the given seed and number of octaves.
    #[new]
    pub fn new(seed: u32, octaves: u32) -> Self {
        Self {
            inner: FractalNoise::new(seed, octaves),
        }
    }

    /// Returns the seed.
    pub fn seed(&self) -> u32 {
        self.inner.seed()
    }

    /// Fractional Brownian motion at (x, y, z).
    pub fn fbm(&self, x: f64, y: f64, z: f64) -> f64 {
        self.inner.sample(x, y, z)
    }

    /// Turbulence (absolute fBm) at (x, y, z).
    pub fn turbulence(&self, x: f64, y: f64, z: f64) -> f64 {
        self.inner.turbulence(x, y, z)
    }

    /// Ridged multifractal at (x, y, z).
    pub fn ridged(&self, x: f64, y: f64, z: f64) -> f64 {
        self.inner.ridged(x, y, z)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyValueNoise3D>()?;
    m.add_class::<PyFractalNoise>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_value_noise_instantiation() {
        let n = PyValueNoise3D::new(42);
        assert_eq!(n.seed(), 42);
        let v = n.sample(1.0, 2.0, 3.0);
        assert!((-1.0..=1.0).contains(&v));
    }

    #[test]
    fn test_fractal_noise_instantiation() {
        let f = PyFractalNoise::new(123, 4);
        assert_eq!(f.seed(), 123);
        let v = f.turbulence(1.0, 2.0, 3.0);
        assert!(v >= 0.0);
    }
}
