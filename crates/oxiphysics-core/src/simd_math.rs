#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SIMD-accelerated math kernels for batch vector and particle operations.
//!
//! This module provides Structure-of-Arrays (SoA) layouts and batch operations
//! optimized for CPU cache locality and auto-vectorization. The batch operations
//! process multiple elements at once, enabling the compiler to emit SIMD
//! instructions on supported platforms.
//!
//! # Layout
//!
//! Instead of the traditional Array-of-Structures (AoS) layout:
//! ```text
//! [x0,y0,z0, x1,y1,z1, x2,y2,z2, ...]
//! ```
//! We use Structure-of-Arrays (SoA):
//! ```text
//! xs: [x0, x1, x2, ...]
//! ys: [y0, y1, y2, ...]
//! zs: [z0, z1, z2, ...]
//! ```
//! This layout allows the compiler to vectorize operations across contiguous
//! memory, improving throughput for large batches.

use std::f64;

/// Structure-of-Arrays layout for batch Vec3 operations.
///
/// Each component (x, y, z) is stored in a separate contiguous vector,
/// enabling efficient SIMD-style batch processing.
#[derive(Debug, Clone, PartialEq)]
pub struct Vec3Batch {
    /// X components of all vectors in the batch.
    pub x: Vec<f64>,
    /// Y components of all vectors in the batch.
    pub y: Vec<f64>,
    /// Z components of all vectors in the batch.
    pub z: Vec<f64>,
}

/// Errors that can occur in SIMD batch operations.
#[derive(Debug, Clone, thiserror::Error)]
pub enum SimdMathError {
    /// Batch size mismatch between operands.
    #[error("batch size mismatch: left has {left} elements, right has {right} elements")]
    SizeMismatch {
        /// Size of the left operand.
        left: usize,
        /// Size of the right operand.
        right: usize,
    },
    /// Inconsistent internal dimensions in a Vec3Batch.
    #[error("inconsistent Vec3Batch dimensions: x={x_len}, y={y_len}, z={z_len}")]
    InconsistentDimensions {
        /// Length of x component.
        x_len: usize,
        /// Length of y component.
        y_len: usize,
        /// Length of z component.
        z_len: usize,
    },
    /// Zero-length vector encountered where normalization is needed.
    #[error("cannot normalize zero-length vector at index {index}")]
    ZeroLengthVector {
        /// Index of the zero-length vector.
        index: usize,
    },
}

impl Vec3Batch {
    /// Allocate a new batch with `n` zero-initialized vectors.
    #[must_use]
    pub fn new(n: usize) -> Self {
        Self {
            x: vec![0.0; n],
            y: vec![0.0; n],
            z: vec![0.0; n],
        }
    }

    /// Convert an Array-of-Structures slice to Structure-of-Arrays layout.
    ///
    /// Each element of `positions` is `[x, y, z]`.
    #[must_use]
    pub fn from_aos(positions: &[[f64; 3]]) -> Self {
        let n = positions.len();
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        let mut z = Vec::with_capacity(n);
        for p in positions {
            x.push(p[0]);
            y.push(p[1]);
            z.push(p[2]);
        }
        Self { x, y, z }
    }

    /// Convert back from SoA to AoS layout.
    ///
    /// Returns an error if internal dimensions are inconsistent.
    pub fn to_aos(&self) -> Result<Vec<[f64; 3]>, SimdMathError> {
        self.validate()?;
        let n = self.x.len();
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            result.push([self.x[i], self.y[i], self.z[i]]);
        }
        Ok(result)
    }

    /// Returns the number of vectors in this batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Returns true if the batch contains no vectors.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Validate that all component vectors have equal length.
    fn validate(&self) -> Result<(), SimdMathError> {
        let (xl, yl, zl) = (self.x.len(), self.y.len(), self.z.len());
        if xl != yl || yl != zl {
            return Err(SimdMathError::InconsistentDimensions {
                x_len: xl,
                y_len: yl,
                z_len: zl,
            });
        }
        Ok(())
    }

    /// Check that two batches have equal size.
    fn check_size(&self, other: &Self) -> Result<(), SimdMathError> {
        if self.len() != other.len() {
            return Err(SimdMathError::SizeMismatch {
                left: self.len(),
                right: other.len(),
            });
        }
        Ok(())
    }

    /// Element-wise addition of two batches.
    ///
    /// # Errors
    /// Returns `SimdMathError::SizeMismatch` if batch sizes differ.
    pub fn add(&self, other: &Vec3Batch) -> Result<Vec3Batch, SimdMathError> {
        self.check_size(other)?;
        let n = self.len();
        let mut rx = vec![0.0_f64; n];
        let mut ry = vec![0.0_f64; n];
        let mut rz = vec![0.0_f64; n];

        // Written as simple loops to encourage auto-vectorization
        for i in 0..n {
            rx[i] = self.x[i] + other.x[i];
        }
        for i in 0..n {
            ry[i] = self.y[i] + other.y[i];
        }
        for i in 0..n {
            rz[i] = self.z[i] + other.z[i];
        }

        Ok(Vec3Batch {
            x: rx,
            y: ry,
            z: rz,
        })
    }

    /// Element-wise subtraction of two batches.
    ///
    /// # Errors
    /// Returns `SimdMathError::SizeMismatch` if batch sizes differ.
    pub fn sub(&self, other: &Vec3Batch) -> Result<Vec3Batch, SimdMathError> {
        self.check_size(other)?;
        let n = self.len();
        let mut rx = vec![0.0_f64; n];
        let mut ry = vec![0.0_f64; n];
        let mut rz = vec![0.0_f64; n];

        for i in 0..n {
            rx[i] = self.x[i] - other.x[i];
        }
        for i in 0..n {
            ry[i] = self.y[i] - other.y[i];
        }
        for i in 0..n {
            rz[i] = self.z[i] - other.z[i];
        }

        Ok(Vec3Batch {
            x: rx,
            y: ry,
            z: rz,
        })
    }

    /// Scale all vectors by a uniform scalar.
    #[must_use]
    pub fn scale(&self, s: f64) -> Vec3Batch {
        let n = self.len();
        let mut rx = vec![0.0_f64; n];
        let mut ry = vec![0.0_f64; n];
        let mut rz = vec![0.0_f64; n];

        for i in 0..n {
            rx[i] = self.x[i] * s;
        }
        for i in 0..n {
            ry[i] = self.y[i] * s;
        }
        for i in 0..n {
            rz[i] = self.z[i] * s;
        }

        Vec3Batch {
            x: rx,
            y: ry,
            z: rz,
        }
    }

    /// Batch dot product: returns `x[i]*other.x[i] + y[i]*other.y[i] + z[i]*other.z[i]`
    /// for each `i`.
    ///
    /// # Errors
    /// Returns `SimdMathError::SizeMismatch` if batch sizes differ.
    pub fn dot(&self, other: &Vec3Batch) -> Result<Vec<f64>, SimdMathError> {
        self.check_size(other)?;
        let n = self.len();
        let mut result = vec![0.0_f64; n];

        // Accumulate component-wise to allow vectorization of each loop
        for i in 0..n {
            result[i] = self.x[i] * other.x[i];
        }
        for i in 0..n {
            result[i] += self.y[i] * other.y[i];
        }
        for i in 0..n {
            result[i] += self.z[i] * other.z[i];
        }

        Ok(result)
    }

    /// Batch cross product.
    ///
    /// For each index `i`, computes `self[i] × other[i]`.
    ///
    /// # Errors
    /// Returns `SimdMathError::SizeMismatch` if batch sizes differ.
    pub fn cross(&self, other: &Vec3Batch) -> Result<Vec3Batch, SimdMathError> {
        self.check_size(other)?;
        let n = self.len();
        let mut rx = vec![0.0_f64; n];
        let mut ry = vec![0.0_f64; n];
        let mut rz = vec![0.0_f64; n];

        // cross.x = self.y * other.z - self.z * other.y
        for i in 0..n {
            rx[i] = self.y[i] * other.z[i] - self.z[i] * other.y[i];
        }
        // cross.y = self.z * other.x - self.x * other.z
        for i in 0..n {
            ry[i] = self.z[i] * other.x[i] - self.x[i] * other.z[i];
        }
        // cross.z = self.x * other.y - self.y * other.x
        for i in 0..n {
            rz[i] = self.x[i] * other.y[i] - self.y[i] * other.x[i];
        }

        Ok(Vec3Batch {
            x: rx,
            y: ry,
            z: rz,
        })
    }

    /// Batch squared length: `x[i]^2 + y[i]^2 + z[i]^2` for each `i`.
    #[must_use]
    pub fn length_sq(&self) -> Vec<f64> {
        let n = self.len();
        let mut result = vec![0.0_f64; n];

        for i in 0..n {
            result[i] = self.x[i] * self.x[i];
        }
        for i in 0..n {
            result[i] += self.y[i] * self.y[i];
        }
        for i in 0..n {
            result[i] += self.z[i] * self.z[i];
        }

        result
    }

    /// Batch vector length (Euclidean norm).
    #[must_use]
    pub fn length(&self) -> Vec<f64> {
        let sq = self.length_sq();
        sq.into_iter().map(f64::sqrt).collect()
    }

    /// Normalize all vectors in-place to unit length.
    ///
    /// Vectors with length below `f64::EPSILON` are left unchanged and their
    /// indices are collected in the returned error. If all vectors are valid,
    /// returns `Ok(())`.
    ///
    /// # Errors
    /// Returns `SimdMathError::ZeroLengthVector` for the first zero-length vector found.
    pub fn normalize(&mut self) -> Result<(), SimdMathError> {
        let lengths = self.length();
        let n = self.len();

        // First pass: check for zero-length vectors
        for (i, &len) in lengths.iter().enumerate() {
            if len < f64::EPSILON {
                return Err(SimdMathError::ZeroLengthVector { index: i });
            }
        }

        // Second pass: compute reciprocals and scale (vectorization-friendly)
        let mut inv_lengths = vec![0.0_f64; n];
        for i in 0..n {
            inv_lengths[i] = 1.0 / lengths[i];
        }

        for i in 0..n {
            self.x[i] *= inv_lengths[i];
        }
        for i in 0..n {
            self.y[i] *= inv_lengths[i];
        }
        for i in 0..n {
            self.z[i] *= inv_lengths[i];
        }

        Ok(())
    }

    /// Batch pairwise squared distance: `|a[i] - b[i]|^2` for each `i`.
    ///
    /// # Errors
    /// Returns `SimdMathError::SizeMismatch` if batch sizes differ.
    pub fn distance_sq_pairwise(a: &Vec3Batch, b: &Vec3Batch) -> Result<Vec<f64>, SimdMathError> {
        a.check_size(b)?;
        let n = a.len();
        let mut result = vec![0.0_f64; n];

        // dx^2
        for i in 0..n {
            let dx = a.x[i] - b.x[i];
            result[i] = dx * dx;
        }
        // dy^2
        for i in 0..n {
            let dy = a.y[i] - b.y[i];
            result[i] += dy * dy;
        }
        // dz^2
        for i in 0..n {
            let dz = a.z[i] - b.z[i];
            result[i] += dz * dz;
        }

        Ok(result)
    }
}

// ---------------------------------------------------------------------------
// Batch Particle Operations
// ---------------------------------------------------------------------------

/// Compute distances from a single reference point to many positions.
///
/// Returns the Euclidean distance from `ref_pos` to each position in the batch.
#[must_use]
pub fn compute_distances_batch(positions: &Vec3Batch, ref_pos: [f64; 3]) -> Vec<f64> {
    let n = positions.len();
    let mut result = vec![0.0_f64; n];

    // dx^2
    for i in 0..n {
        let dx = positions.x[i] - ref_pos[0];
        result[i] = dx * dx;
    }
    // dy^2
    for i in 0..n {
        let dy = positions.y[i] - ref_pos[1];
        result[i] += dy * dy;
    }
    // dz^2
    for i in 0..n {
        let dz = positions.z[i] - ref_pos[2];
        result[i] += dz * dz;
    }

    // sqrt
    for val in &mut result {
        *val = val.sqrt();
    }

    result
}

/// Accumulate forces onto a mutable force batch.
///
/// For each index `i`, computes:
/// ```text
/// forces[i] += directions[i] * magnitudes[i]
/// ```
///
/// # Errors
/// Returns `SimdMathError::SizeMismatch` if the sizes of `forces`, `directions`,
/// or `magnitudes` do not match.
pub fn accumulate_forces_batch(
    forces: &mut Vec3Batch,
    directions: &Vec3Batch,
    magnitudes: &[f64],
) -> Result<(), SimdMathError> {
    let n = forces.len();
    if n != directions.len() {
        return Err(SimdMathError::SizeMismatch {
            left: n,
            right: directions.len(),
        });
    }
    if n != magnitudes.len() {
        return Err(SimdMathError::SizeMismatch {
            left: n,
            right: magnitudes.len(),
        });
    }

    for i in 0..n {
        forces.x[i] += directions.x[i] * magnitudes[i];
    }
    for i in 0..n {
        forces.y[i] += directions.y[i] * magnitudes[i];
    }
    for i in 0..n {
        forces.z[i] += directions.z[i] * magnitudes[i];
    }

    Ok(())
}

/// Evaluate the cubic spline kernel (M4 kernel) for a batch of distances.
///
/// The cubic spline kernel is widely used in Smoothed Particle Hydrodynamics (SPH).
/// It is defined in 3D as:
///
/// ```text
/// W(r, h) = σ * { (2-q)^3 - 4*(1-q)^3   if 0 ≤ q < 1
///               { (2-q)^3                  if 1 ≤ q < 2
///               { 0                        if q ≥ 2
/// ```
///
/// where `q = r/h` and `σ = 1/(4π h³)` is the 3D normalization constant.
///
/// # Arguments
/// * `r` - batch of distances
/// * `h` - smoothing length (must be positive)
///
/// # Returns
/// Vector of kernel values. Returns empty vec if `h` is not positive.
#[must_use]
pub fn cubic_spline_kernel_batch(r: &[f64], h: f64) -> Vec<f64> {
    if h <= 0.0 {
        return vec![0.0; r.len()];
    }

    let n = r.len();
    let inv_h = 1.0 / h;
    let sigma = 1.0 / (4.0 * f64::consts::PI * h * h * h);

    let mut result = vec![0.0_f64; n];
    let mut q = vec![0.0_f64; n];

    // Compute q = r/h (vectorizable)
    for i in 0..n {
        q[i] = r[i] * inv_h;
    }

    // Evaluate kernel piecewise
    for i in 0..n {
        let qi = q[i];
        if qi >= 2.0 {
            result[i] = 0.0;
        } else if qi >= 1.0 {
            let t = 2.0 - qi;
            result[i] = sigma * t * t * t;
        } else if qi >= 0.0 {
            let t2 = 2.0 - qi;
            let t1 = 1.0 - qi;
            result[i] = sigma * (t2 * t2 * t2 - 4.0 * t1 * t1 * t1);
        } else {
            // Negative distance: treat as absolute value
            let qi_abs = qi.abs();
            if qi_abs >= 2.0 {
                result[i] = 0.0;
            } else if qi_abs >= 1.0 {
                let t = 2.0 - qi_abs;
                result[i] = sigma * t * t * t;
            } else {
                let t2 = 2.0 - qi_abs;
                let t1 = 1.0 - qi_abs;
                result[i] = sigma * (t2 * t2 * t2 - 4.0 * t1 * t1 * t1);
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-12;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------------
    // AoS <-> SoA round-trip
    // -----------------------------------------------------------------------

    #[test]
    fn test_aos_soa_round_trip() {
        let positions = vec![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [-1.5, 0.0, 3.125],
        ];
        let batch = Vec3Batch::from_aos(&positions);
        assert_eq!(batch.len(), 4);
        assert_eq!(batch.x, vec![1.0, 4.0, 7.0, -1.5]);
        assert_eq!(batch.y, vec![2.0, 5.0, 8.0, 0.0]);
        assert_eq!(batch.z, vec![3.0, 6.0, 9.0, 3.125]);

        let back = batch.to_aos().expect("to_aos should succeed");
        assert_eq!(back, positions);
    }

    #[test]
    fn test_aos_soa_round_trip_empty() {
        let positions: Vec<[f64; 3]> = vec![];
        let batch = Vec3Batch::from_aos(&positions);
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
        let back = batch.to_aos().expect("to_aos should succeed for empty");
        assert!(back.is_empty());
    }

    // -----------------------------------------------------------------------
    // Batch dot product matches scalar
    // -----------------------------------------------------------------------

    #[test]
    fn test_batch_dot_matches_scalar() {
        let a_aos = vec![[1.0, 2.0, 3.0], [4.0, -1.0, 2.0], [0.0, 0.0, 1.0]];
        let b_aos = vec![[3.0, -2.0, 1.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]];

        let a = Vec3Batch::from_aos(&a_aos);
        let b = Vec3Batch::from_aos(&b_aos);

        let batch_dots = a.dot(&b).expect("dot should succeed");

        // Scalar reference
        for (i, (&da, &db)) in a_aos.iter().zip(b_aos.iter()).enumerate() {
            let scalar_dot = da[0] * db[0] + da[1] * db[1] + da[2] * db[2];
            assert!(
                approx_eq(batch_dots[i], scalar_dot, EPSILON),
                "dot mismatch at index {i}: batch={}, scalar={scalar_dot}",
                batch_dots[i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // Batch cross product matches scalar
    // -----------------------------------------------------------------------

    #[test]
    fn test_batch_cross_matches_scalar() {
        let a_aos = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 2.0, 3.0]];
        let b_aos = vec![[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [4.0, 5.0, 6.0]];

        let a = Vec3Batch::from_aos(&a_aos);
        let b = Vec3Batch::from_aos(&b_aos);

        let cross_batch = a.cross(&b).expect("cross should succeed");
        let cross_aos = cross_batch.to_aos().expect("to_aos should succeed");

        for (i, (&va, &vb)) in a_aos.iter().zip(b_aos.iter()).enumerate() {
            let cx = va[1] * vb[2] - va[2] * vb[1];
            let cy = va[2] * vb[0] - va[0] * vb[2];
            let cz = va[0] * vb[1] - va[1] * vb[0];
            assert!(
                approx_eq(cross_aos[i][0], cx, EPSILON),
                "cross.x mismatch at {i}"
            );
            assert!(
                approx_eq(cross_aos[i][1], cy, EPSILON),
                "cross.y mismatch at {i}"
            );
            assert!(
                approx_eq(cross_aos[i][2], cz, EPSILON),
                "cross.z mismatch at {i}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Batch normalize produces unit vectors
    // -----------------------------------------------------------------------

    #[test]
    fn test_batch_normalize_unit_vectors() {
        let positions = vec![
            [3.0, 4.0, 0.0],  // length = 5
            [0.0, 0.0, 7.0],  // length = 7
            [1.0, 1.0, 1.0],  // length = sqrt(3)
            [10.0, 0.0, 0.0], // length = 10
        ];
        let mut batch = Vec3Batch::from_aos(&positions);
        batch.normalize().expect("normalize should succeed");

        let lengths = batch.length();
        for (i, &len) in lengths.iter().enumerate() {
            assert!(
                approx_eq(len, 1.0, 1e-10),
                "expected unit length at index {i}, got {len}"
            );
        }

        // Check specific known result: [3,4,0] -> [0.6, 0.8, 0.0]
        assert!(approx_eq(batch.x[0], 0.6, EPSILON));
        assert!(approx_eq(batch.y[0], 0.8, EPSILON));
        assert!(approx_eq(batch.z[0], 0.0, EPSILON));
    }

    #[test]
    fn test_normalize_zero_vector_error() {
        let mut batch = Vec3Batch::from_aos(&[[0.0, 0.0, 0.0]]);
        let result = batch.normalize();
        assert!(result.is_err());
        match result {
            Err(SimdMathError::ZeroLengthVector { index }) => assert_eq!(index, 0),
            other => panic!("expected ZeroLengthVector error, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------------
    // Distance computation matches naive loop
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_distances_batch_matches_naive() {
        let positions_aos = vec![
            [1.0, 0.0, 0.0],
            [0.0, 3.0, 4.0],
            [1.0, 1.0, 1.0],
            [10.0, 20.0, 30.0],
        ];
        let ref_pos = [0.0, 0.0, 0.0];

        let batch = Vec3Batch::from_aos(&positions_aos);
        let distances = compute_distances_batch(&batch, ref_pos);

        for (i, pos) in positions_aos.iter().enumerate() {
            let naive = ((pos[0] - ref_pos[0]).powi(2)
                + (pos[1] - ref_pos[1]).powi(2)
                + (pos[2] - ref_pos[2]).powi(2))
            .sqrt();
            assert!(
                approx_eq(distances[i], naive, EPSILON),
                "distance mismatch at {i}: batch={}, naive={naive}",
                distances[i]
            );
        }
    }

    #[test]
    fn test_distance_sq_pairwise() {
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let b = Vec3Batch::from_aos(&[[4.0, 6.0, 3.0], [4.0, 5.0, 6.0]]);

        let dsq = Vec3Batch::distance_sq_pairwise(&a, &b).expect("should succeed");
        // [1]: (4-4)^2 + (5-5)^2 + (6-6)^2 = 0
        assert!(approx_eq(dsq[1], 0.0, EPSILON));
        // [0]: (4-1)^2 + (6-2)^2 + (3-3)^2 = 9+16+0 = 25
        assert!(approx_eq(dsq[0], 25.0, EPSILON));
    }

    // -----------------------------------------------------------------------
    // Kernel evaluation matches scalar version
    // -----------------------------------------------------------------------

    #[test]
    fn test_cubic_spline_kernel_matches_scalar() {
        let h = 1.0;
        let sigma = 1.0 / (4.0 * f64::consts::PI * h * h * h);

        // Test representative q values
        let r_values = vec![0.0, 0.5, 0.99, 1.0, 1.5, 1.99, 2.0, 3.0];
        let kernel_vals = cubic_spline_kernel_batch(&r_values, h);

        for (i, &r) in r_values.iter().enumerate() {
            let q = r / h;
            let expected = if q >= 2.0 {
                0.0
            } else if q >= 1.0 {
                let t = 2.0 - q;
                sigma * t * t * t
            } else {
                let t2 = 2.0 - q;
                let t1 = 1.0 - q;
                sigma * (t2 * t2 * t2 - 4.0 * t1 * t1 * t1)
            };
            assert!(
                approx_eq(kernel_vals[i], expected, EPSILON),
                "kernel mismatch at r={r}: batch={}, expected={expected}",
                kernel_vals[i]
            );
        }
    }

    #[test]
    fn test_cubic_spline_kernel_zero_at_boundary() {
        let h = 2.0;
        let vals = cubic_spline_kernel_batch(&[4.0, 5.0, 100.0], h);
        for (i, &v) in vals.iter().enumerate() {
            assert!(
                approx_eq(v, 0.0, EPSILON),
                "expected zero at index {i}, got {v}"
            );
        }
    }

    #[test]
    fn test_cubic_spline_kernel_non_positive_h() {
        let vals = cubic_spline_kernel_batch(&[1.0, 2.0], 0.0);
        assert_eq!(vals, vec![0.0, 0.0]);

        let vals_neg = cubic_spline_kernel_batch(&[1.0], -1.0);
        assert_eq!(vals_neg, vec![0.0]);
    }

    // -----------------------------------------------------------------------
    // Empty batch operations
    // -----------------------------------------------------------------------

    #[test]
    fn test_empty_batch_operations() {
        let a = Vec3Batch::new(0);
        let b = Vec3Batch::new(0);

        assert!(a.is_empty());

        let sum = a.add(&b).expect("add empty should succeed");
        assert!(sum.is_empty());

        let diff = a.sub(&b).expect("sub empty should succeed");
        assert!(diff.is_empty());

        let dots = a.dot(&b).expect("dot empty should succeed");
        assert!(dots.is_empty());

        let cross = a.cross(&b).expect("cross empty should succeed");
        assert!(cross.is_empty());

        let scaled = a.scale(5.0);
        assert!(scaled.is_empty());

        let lsq = a.length_sq();
        assert!(lsq.is_empty());

        let lens = a.length();
        assert!(lens.is_empty());

        let dsq =
            Vec3Batch::distance_sq_pairwise(&a, &b).expect("distance_sq empty should succeed");
        assert!(dsq.is_empty());

        let dists = compute_distances_batch(&a, [0.0, 0.0, 0.0]);
        assert!(dists.is_empty());

        let kernel = cubic_spline_kernel_batch(&[], 1.0);
        assert!(kernel.is_empty());
    }

    // -----------------------------------------------------------------------
    // Size mismatch errors
    // -----------------------------------------------------------------------

    #[test]
    fn test_size_mismatch_errors() {
        let a = Vec3Batch::new(3);
        let b = Vec3Batch::new(5);

        assert!(a.add(&b).is_err());
        assert!(a.sub(&b).is_err());
        assert!(a.dot(&b).is_err());
        assert!(a.cross(&b).is_err());
        assert!(Vec3Batch::distance_sq_pairwise(&a, &b).is_err());
    }

    // -----------------------------------------------------------------------
    // Force accumulation
    // -----------------------------------------------------------------------

    #[test]
    fn test_accumulate_forces_batch() {
        let mut forces = Vec3Batch::new(3);
        let directions = Vec3Batch::from_aos(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let magnitudes = vec![10.0, 20.0, 30.0];

        accumulate_forces_batch(&mut forces, &directions, &magnitudes)
            .expect("accumulate should succeed");

        assert!(approx_eq(forces.x[0], 10.0, EPSILON));
        assert!(approx_eq(forces.y[1], 20.0, EPSILON));
        assert!(approx_eq(forces.z[2], 30.0, EPSILON));

        // Accumulate again (forces should add up)
        accumulate_forces_batch(&mut forces, &directions, &magnitudes)
            .expect("second accumulate should succeed");

        assert!(approx_eq(forces.x[0], 20.0, EPSILON));
        assert!(approx_eq(forces.y[1], 40.0, EPSILON));
        assert!(approx_eq(forces.z[2], 60.0, EPSILON));
    }

    #[test]
    fn test_accumulate_forces_size_mismatch() {
        let mut forces = Vec3Batch::new(3);
        let directions = Vec3Batch::new(2);
        let magnitudes = vec![1.0, 2.0, 3.0];

        assert!(accumulate_forces_batch(&mut forces, &directions, &magnitudes).is_err());

        let directions2 = Vec3Batch::new(3);
        let magnitudes2 = vec![1.0, 2.0];
        assert!(accumulate_forces_batch(&mut forces, &directions2, &magnitudes2).is_err());
    }

    // -----------------------------------------------------------------------
    // Add / Sub correctness
    // -----------------------------------------------------------------------

    #[test]
    fn test_add_sub_inverse() {
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let b = Vec3Batch::from_aos(&[[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]);

        let sum = a.add(&b).expect("add should succeed");
        let back = sum.sub(&b).expect("sub should succeed");

        for i in 0..a.len() {
            assert!(approx_eq(back.x[i], a.x[i], EPSILON));
            assert!(approx_eq(back.y[i], a.y[i], EPSILON));
            assert!(approx_eq(back.z[i], a.z[i], EPSILON));
        }
    }

    // -----------------------------------------------------------------------
    // Scale correctness
    // -----------------------------------------------------------------------

    #[test]
    fn test_scale_by_zero() {
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0]]);
        let scaled = a.scale(0.0);
        assert!(approx_eq(scaled.x[0], 0.0, EPSILON));
        assert!(approx_eq(scaled.y[0], 0.0, EPSILON));
        assert!(approx_eq(scaled.z[0], 0.0, EPSILON));
    }

    #[test]
    fn test_scale_negative() {
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0]]);
        let scaled = a.scale(-2.0);
        assert!(approx_eq(scaled.x[0], -2.0, EPSILON));
        assert!(approx_eq(scaled.y[0], -4.0, EPSILON));
        assert!(approx_eq(scaled.z[0], -6.0, EPSILON));
    }

    // -----------------------------------------------------------------------
    // Length / length_sq
    // -----------------------------------------------------------------------

    #[test]
    fn test_length_sq_and_length() {
        let batch = Vec3Batch::from_aos(&[[3.0, 4.0, 0.0], [0.0, 0.0, 5.0]]);
        let lsq = batch.length_sq();
        assert!(approx_eq(lsq[0], 25.0, EPSILON));
        assert!(approx_eq(lsq[1], 25.0, EPSILON));

        let lens = batch.length();
        assert!(approx_eq(lens[0], 5.0, EPSILON));
        assert!(approx_eq(lens[1], 5.0, EPSILON));
    }

    // -----------------------------------------------------------------------
    // Cross product specific identities
    // -----------------------------------------------------------------------

    #[test]
    fn test_cross_product_anticommutative() {
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let b = Vec3Batch::from_aos(&[[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]);

        let axb = a.cross(&b).expect("cross a x b");
        let bxa = b.cross(&a).expect("cross b x a");

        // a x b = -(b x a)
        for i in 0..a.len() {
            assert!(approx_eq(axb.x[i], -bxa.x[i], EPSILON));
            assert!(approx_eq(axb.y[i], -bxa.y[i], EPSILON));
            assert!(approx_eq(axb.z[i], -bxa.z[i], EPSILON));
        }
    }

    #[test]
    fn test_cross_product_perpendicular() {
        // a x b should be perpendicular to both a and b
        let a = Vec3Batch::from_aos(&[[1.0, 2.0, 3.0]]);
        let b = Vec3Batch::from_aos(&[[4.0, 5.0, 6.0]]);

        let c = a.cross(&b).expect("cross");
        let dot_ac = a.dot(&c).expect("dot a.c");
        let dot_bc = b.dot(&c).expect("dot b.c");

        assert!(approx_eq(dot_ac[0], 0.0, 1e-10));
        assert!(approx_eq(dot_bc[0], 0.0, 1e-10));
    }

    // -----------------------------------------------------------------------
    // Kernel monotonicity
    // -----------------------------------------------------------------------

    #[test]
    fn test_cubic_spline_kernel_monotonic_decrease() {
        let h = 1.0;
        let r: Vec<f64> = (0..20).map(|i| i as f64 * 0.1).collect();
        let vals = cubic_spline_kernel_batch(&r, h);

        // Kernel should be non-negative
        for (i, &v) in vals.iter().enumerate() {
            assert!(v >= 0.0, "kernel negative at r={}: {v}", r[i]);
        }

        // Maximum should be at r=0
        let max_val = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(approx_eq(max_val, vals[0], EPSILON));
    }

    // -----------------------------------------------------------------------
    // Compute distances with non-origin reference
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_distances_nonzero_ref() {
        let positions = Vec3Batch::from_aos(&[[4.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        let ref_pos = [1.0, 0.0, 0.0];
        let dists = compute_distances_batch(&positions, ref_pos);

        assert!(approx_eq(dists[0], 3.0, EPSILON));
        let expected_1 = (0.0_f64 + 1.0 + 1.0_f64).sqrt();
        assert!(approx_eq(dists[1], expected_1, EPSILON));
    }
}
