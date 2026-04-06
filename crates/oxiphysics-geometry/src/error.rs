// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-geometry
//!
//! This module provides structured error handling for geometry operations,
//! including mesh validation, ray-cast diagnostics, parameter checking,
//! and I/O errors.  All public types implement `std::error::Error` via
//! `thiserror`.

use thiserror::Error;

// ── Core error type ──────────────────────────────────────────────────────────

/// Main error type for the geometry module.
#[derive(Debug, Error)]
pub enum Error {
    /// Generic catch-all error with a human-readable message.
    #[error("{0}")]
    General(String),

    /// A required parameter was out of its valid range.
    #[error("parameter '{name}' = {value} is out of range [{min}, {max}]")]
    OutOfRange {
        /// Name of the parameter.
        name: &'static str,
        /// Actual value supplied.
        value: f64,
        /// Inclusive lower bound.
        min: f64,
        /// Inclusive upper bound.
        max: f64,
    },

    /// A mesh is topologically or geometrically invalid.
    #[error("mesh validation failed: {reason}")]
    InvalidMesh {
        /// Human-readable description of the problem.
        reason: String,
    },

    /// A buffer had the wrong length.
    #[error("buffer length mismatch: expected {expected}, got {actual}")]
    LengthMismatch {
        /// Expected length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },

    /// Numerical computation did not converge.
    #[error(
        "numerical computation '{operation}' did not converge after {iterations} iterations (residual {residual:.3e})"
    )]
    ConvergenceFailure {
        /// Name of the operation that failed.
        operation: &'static str,
        /// Number of iterations attempted.
        iterations: usize,
        /// Final residual value.
        residual: f64,
    },

    /// An index was out of bounds for the given container.
    #[error("index out of bounds: {index} >= {len}")]
    IndexOutOfBounds {
        /// The invalid index.
        index: usize,
        /// The length of the container.
        len: usize,
    },

    /// A shape requires at least a minimum number of vertices/points.
    #[error("shape requires at least {required} vertices/points, got {actual}")]
    TooFewPoints {
        /// Minimum required.
        required: usize,
        /// Actual count.
        actual: usize,
    },

    /// A degenerate geometry was encountered (zero area, zero volume, etc.).
    #[error("degenerate geometry: {details}")]
    DegenerateGeometry {
        /// Description of the degenerate case.
        details: String,
    },

    /// Two arrays that must have equal length do not.
    #[error("array dimension mismatch: '{lhs}' has {lhs_len} elements, '{rhs}' has {rhs_len}")]
    DimensionMismatch {
        /// Name of the first array.
        lhs: &'static str,
        /// Length of the first array.
        lhs_len: usize,
        /// Name of the second array.
        rhs: &'static str,
        /// Length of the second array.
        rhs_len: usize,
    },

    /// A requested feature or algorithm is not supported for this input.
    #[error("unsupported operation '{operation}': {reason}")]
    Unsupported {
        /// Name of the operation.
        operation: &'static str,
        /// Reason it is unsupported.
        reason: String,
    },

    /// An I/O error occurred during geometry serialization or deserialization.
    #[error("I/O error in '{context}': {message}")]
    Io {
        /// The operation context (e.g. "serialize heightfield").
        context: &'static str,
        /// Human-readable message.
        message: String,
    },
}

/// Result type alias for geometry operations.
pub type Result<T> = std::result::Result<T, Error>;

// ── Convenience constructors ─────────────────────────────────────────────────

impl Error {
    /// Create a `General` error from any string-like value.
    #[allow(dead_code)]
    pub fn general(msg: impl Into<String>) -> Self {
        Self::General(msg.into())
    }

    /// Create an `OutOfRange` error.
    #[allow(dead_code)]
    pub fn out_of_range(name: &'static str, value: f64, min: f64, max: f64) -> Self {
        Self::OutOfRange {
            name,
            value,
            min,
            max,
        }
    }

    /// Create an `InvalidMesh` error.
    #[allow(dead_code)]
    pub fn invalid_mesh(reason: impl Into<String>) -> Self {
        Self::InvalidMesh {
            reason: reason.into(),
        }
    }

    /// Create a `LengthMismatch` error.
    #[allow(dead_code)]
    pub fn length_mismatch(expected: usize, actual: usize) -> Self {
        Self::LengthMismatch { expected, actual }
    }

    /// Create a `ConvergenceFailure` error.
    #[allow(dead_code)]
    pub fn convergence_failure(operation: &'static str, iterations: usize, residual: f64) -> Self {
        Self::ConvergenceFailure {
            operation,
            iterations,
            residual,
        }
    }

    /// Create an `IndexOutOfBounds` error.
    #[allow(dead_code)]
    pub fn index_out_of_bounds(index: usize, len: usize) -> Self {
        Self::IndexOutOfBounds { index, len }
    }

    /// Create a `TooFewPoints` error.
    #[allow(dead_code)]
    pub fn too_few_points(required: usize, actual: usize) -> Self {
        Self::TooFewPoints { required, actual }
    }

    /// Create a `DegenerateGeometry` error.
    #[allow(dead_code)]
    pub fn degenerate_geometry(details: impl Into<String>) -> Self {
        Self::DegenerateGeometry {
            details: details.into(),
        }
    }

    /// Create a `DimensionMismatch` error.
    #[allow(dead_code)]
    pub fn dimension_mismatch(
        lhs: &'static str,
        lhs_len: usize,
        rhs: &'static str,
        rhs_len: usize,
    ) -> Self {
        Self::DimensionMismatch {
            lhs,
            lhs_len,
            rhs,
            rhs_len,
        }
    }

    /// Create an `Unsupported` error.
    #[allow(dead_code)]
    pub fn unsupported(operation: &'static str, reason: impl Into<String>) -> Self {
        Self::Unsupported {
            operation,
            reason: reason.into(),
        }
    }

    /// Create an `Io` error.
    #[allow(dead_code)]
    pub fn io(context: &'static str, message: impl Into<String>) -> Self {
        Self::Io {
            context,
            message: message.into(),
        }
    }

    /// Returns `true` if this is a `General` error.
    #[allow(dead_code)]
    pub fn is_general(&self) -> bool {
        matches!(self, Self::General(_))
    }

    /// Returns `true` if this is a `LengthMismatch` error.
    #[allow(dead_code)]
    pub fn is_length_mismatch(&self) -> bool {
        matches!(self, Self::LengthMismatch { .. })
    }

    /// Returns `true` if this is a convergence failure.
    #[allow(dead_code)]
    pub fn is_convergence_failure(&self) -> bool {
        matches!(self, Self::ConvergenceFailure { .. })
    }

    /// Returns `true` if this is an index-out-of-bounds error.
    #[allow(dead_code)]
    pub fn is_index_out_of_bounds(&self) -> bool {
        matches!(self, Self::IndexOutOfBounds { .. })
    }

    /// Returns `true` if this error indicates degenerate geometry.
    #[allow(dead_code)]
    pub fn is_degenerate(&self) -> bool {
        matches!(self, Self::DegenerateGeometry { .. })
    }
}

// ── Validation helpers ───────────────────────────────────────────────────────

/// Assert that `value` lies in `[min, max]`, returning `Err(Error::OutOfRange)`
/// if it does not.
#[allow(dead_code)]
pub fn check_range(name: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value >= min && value <= max {
        Ok(())
    } else {
        Err(Error::out_of_range(name, value, min, max))
    }
}

/// Assert that `len == expected`, returning `Err(Error::LengthMismatch)` if not.
#[allow(dead_code)]
pub fn check_len(expected: usize, actual: usize) -> Result<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(Error::length_mismatch(expected, actual))
    }
}

/// Assert that `index < len`, returning `Err(Error::IndexOutOfBounds)` if not.
#[allow(dead_code)]
pub fn check_index(index: usize, len: usize) -> Result<()> {
    if index < len {
        Ok(())
    } else {
        Err(Error::index_out_of_bounds(index, len))
    }
}

/// Assert that `count >= required`, returning `Err(Error::TooFewPoints)` if not.
#[allow(dead_code)]
pub fn check_min_points(required: usize, actual: usize) -> Result<()> {
    if actual >= required {
        Ok(())
    } else {
        Err(Error::too_few_points(required, actual))
    }
}

/// Assert that `lhs_len == rhs_len`, returning `Err(Error::DimensionMismatch)`.
#[allow(dead_code)]
pub fn check_dim_match(
    lhs: &'static str,
    lhs_len: usize,
    rhs: &'static str,
    rhs_len: usize,
) -> Result<()> {
    if lhs_len == rhs_len {
        Ok(())
    } else {
        Err(Error::dimension_mismatch(lhs, lhs_len, rhs, rhs_len))
    }
}

/// Assert that a value is strictly positive, returning `Err(Error::OutOfRange)`.
#[allow(dead_code)]
pub fn check_positive(name: &'static str, value: f64) -> Result<()> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(Error::out_of_range(name, value, f64::EPSILON, f64::MAX))
    }
}

/// Assert that a value is non-negative, returning `Err(Error::OutOfRange)`.
#[allow(dead_code)]
pub fn check_non_negative(name: &'static str, value: f64) -> Result<()> {
    if value >= 0.0 {
        Ok(())
    } else {
        Err(Error::out_of_range(name, value, 0.0, f64::MAX))
    }
}

/// Assert that a value is finite (not NaN, not infinite).
#[allow(dead_code)]
pub fn check_finite(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::general(format!(
            "parameter '{name}' is not finite: {value}"
        )))
    }
}

/// Check all values in a slice are finite.
#[allow(dead_code)]
pub fn check_finite_slice(name: &'static str, values: &[f64]) -> Result<()> {
    for (i, &v) in values.iter().enumerate() {
        if !v.is_finite() {
            return Err(Error::general(format!(
                "parameter '{name}[{i}]' is not finite: {v}"
            )));
        }
    }
    Ok(())
}

// ── Mesh validation helpers ──────────────────────────────────────────────────

/// Validate that a triangle mesh has consistent vertex/index arrays.
///
/// Checks:
/// 1. `vertices` is non-empty.
/// 2. All indices in `triangles` are less than `vertices.len()`.
/// 3. No triangle has repeated vertex indices (degenerate triangles).
#[allow(dead_code)]
pub fn validate_mesh(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Result<()> {
    check_min_points(1, vertices.len())?;
    for (i, tri) in triangles.iter().enumerate() {
        for &idx in tri {
            if idx >= vertices.len() {
                return Err(Error::InvalidMesh {
                    reason: format!(
                        "triangle {i}: index {idx} >= vertex count {}",
                        vertices.len()
                    ),
                });
            }
        }
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
            return Err(Error::DegenerateGeometry {
                details: format!(
                    "triangle {i} has repeated indices: [{}, {}, {}]",
                    tri[0], tri[1], tri[2]
                ),
            });
        }
    }
    Ok(())
}

/// Validate a height-field descriptor.
///
/// Checks that `rows >= 2`, `cols >= 2`, `scale_x > 0`, `scale_z > 0`, and
/// `heights.len() == rows * cols`.
#[allow(dead_code)]
pub fn validate_heightfield(
    heights: &[f64],
    rows: usize,
    cols: usize,
    scale_x: f64,
    scale_z: f64,
) -> Result<()> {
    if rows < 2 {
        return Err(Error::TooFewPoints {
            required: 2,
            actual: rows,
        });
    }
    if cols < 2 {
        return Err(Error::TooFewPoints {
            required: 2,
            actual: cols,
        });
    }
    check_positive("scale_x", scale_x)?;
    check_positive("scale_z", scale_z)?;
    check_len(rows * cols, heights.len())?;
    check_finite_slice("heights", heights)?;
    Ok(())
}

/// Validate that a ray direction is non-zero and finite.
#[allow(dead_code)]
pub fn validate_ray_dir(dir: [f64; 3]) -> Result<()> {
    check_finite_slice("ray_dir", &dir)?;
    let len_sq = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
    if len_sq < 1e-30 {
        return Err(Error::DegenerateGeometry {
            details: "ray direction is zero (or near-zero)".into(),
        });
    }
    Ok(())
}

/// Validate a point cloud: non-empty, all coordinates finite.
#[allow(dead_code)]
pub fn validate_point_cloud(points: &[[f64; 3]]) -> Result<()> {
    check_min_points(1, points.len())?;
    for (i, p) in points.iter().enumerate() {
        if !p[0].is_finite() || !p[1].is_finite() || !p[2].is_finite() {
            return Err(Error::general(format!(
                "point[{i}] contains non-finite coordinate: {p:?}"
            )));
        }
    }
    Ok(())
}

// ── Error context helpers ────────────────────────────────────────────────────

/// Attach a context string to a `Result`, wrapping the error in a new
/// `General` message that includes the original error's display text.
#[allow(dead_code)]
pub trait WithContext<T> {
    /// Wrap any error with an additional context prefix.
    fn with_context(self, ctx: &str) -> Result<T>;
}

impl<T, E: std::fmt::Display> WithContext<T> for std::result::Result<T, E> {
    fn with_context(self, ctx: &str) -> Result<T> {
        self.map_err(|e| Error::General(format!("{ctx}: {e}")))
    }
}

// ── Iterative-solver convergence tracker ────────────────────────────────────

/// Tracks residual progress for an iterative solver and raises
/// `Error::ConvergenceFailure` when the iteration limit is exceeded.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConvergenceTracker {
    operation: &'static str,
    max_iterations: usize,
    tolerance: f64,
    current_iteration: usize,
    last_residual: f64,
}

#[allow(dead_code)]
impl ConvergenceTracker {
    /// Create a new tracker.
    ///
    /// # Arguments
    /// - `operation` — human-readable name of the algorithm.
    /// - `max_iterations` — maximum allowed iterations before failure.
    /// - `tolerance` — convergence criterion (residual < tolerance ⟹ converged).
    pub fn new(operation: &'static str, max_iterations: usize, tolerance: f64) -> Self {
        Self {
            operation,
            max_iterations,
            tolerance,
            current_iteration: 0,
            last_residual: f64::INFINITY,
        }
    }

    /// Record a residual for the current iteration and advance the counter.
    ///
    /// Returns `Ok(true)` if converged, `Ok(false)` if still iterating,
    /// or `Err(Error::ConvergenceFailure)` if the limit was reached.
    pub fn update(&mut self, residual: f64) -> Result<bool> {
        self.last_residual = residual;
        self.current_iteration += 1;
        if residual < self.tolerance {
            return Ok(true);
        }
        if self.current_iteration >= self.max_iterations {
            return Err(Error::convergence_failure(
                self.operation,
                self.current_iteration,
                residual,
            ));
        }
        Ok(false)
    }

    /// Current iteration count.
    pub fn iterations(&self) -> usize {
        self.current_iteration
    }

    /// Last recorded residual.
    pub fn residual(&self) -> f64 {
        self.last_residual
    }

    /// Reset the tracker for reuse.
    pub fn reset(&mut self) {
        self.current_iteration = 0;
        self.last_residual = f64::INFINITY;
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Error construction ────────────────────────────────────────────────────

    #[test]
    fn test_general_error_display() {
        let e = Error::general("something went wrong");
        assert!(e.to_string().contains("something went wrong"));
    }

    #[test]
    fn test_out_of_range_display() {
        let e = Error::out_of_range("radius", -1.0, 0.0, 100.0);
        let s = e.to_string();
        assert!(s.contains("radius"), "should mention parameter name");
        assert!(s.contains("-1"), "should mention actual value");
    }

    #[test]
    fn test_length_mismatch_display() {
        let e = Error::length_mismatch(10, 7);
        let s = e.to_string();
        assert!(s.contains("10"));
        assert!(s.contains("7"));
    }

    #[test]
    fn test_convergence_failure_display() {
        let e = Error::convergence_failure("smoothing", 100, 0.001);
        let s = e.to_string();
        assert!(s.contains("smoothing"));
        assert!(s.contains("100"));
    }

    #[test]
    fn test_index_out_of_bounds_display() {
        let e = Error::index_out_of_bounds(5, 3);
        let s = e.to_string();
        assert!(s.contains('5'));
        assert!(s.contains('3'));
    }

    #[test]
    fn test_too_few_points_display() {
        let e = Error::too_few_points(3, 1);
        let s = e.to_string();
        assert!(s.contains('3'));
        assert!(s.contains('1'));
    }

    #[test]
    fn test_degenerate_geometry_display() {
        let e = Error::degenerate_geometry("zero area triangle");
        assert!(e.to_string().contains("zero area triangle"));
    }

    #[test]
    fn test_dimension_mismatch_display() {
        let e = Error::dimension_mismatch("positions", 10, "normals", 8);
        let s = e.to_string();
        assert!(s.contains("positions"));
        assert!(s.contains("normals"));
    }

    #[test]
    fn test_unsupported_display() {
        let e = Error::unsupported("CSG union", "non-manifold mesh");
        let s = e.to_string();
        assert!(s.contains("CSG union"));
        assert!(s.contains("non-manifold mesh"));
    }

    #[test]
    fn test_io_error_display() {
        let e = Error::io("serialize heightfield", "disk full");
        let s = e.to_string();
        assert!(s.contains("serialize heightfield"));
        assert!(s.contains("disk full"));
    }

    // ── is_* predicates ───────────────────────────────────────────────────────

    #[test]
    fn test_is_general() {
        assert!(Error::general("x").is_general());
        assert!(!Error::length_mismatch(1, 2).is_general());
    }

    #[test]
    fn test_is_length_mismatch() {
        assert!(Error::length_mismatch(3, 4).is_length_mismatch());
        assert!(!Error::general("x").is_length_mismatch());
    }

    #[test]
    fn test_is_convergence_failure() {
        assert!(Error::convergence_failure("op", 10, 0.1).is_convergence_failure());
        assert!(!Error::general("x").is_convergence_failure());
    }

    #[test]
    fn test_is_degenerate() {
        assert!(Error::degenerate_geometry("zero vol").is_degenerate());
        assert!(!Error::general("x").is_degenerate());
    }

    // ── Validation helpers ────────────────────────────────────────────────────

    #[test]
    fn test_check_range_ok() {
        assert!(check_range("r", 5.0, 0.0, 10.0).is_ok());
    }

    #[test]
    fn test_check_range_below_min() {
        let r = check_range("r", -1.0, 0.0, 10.0);
        assert!(r.is_err());
        assert!(matches!(r.unwrap_err(), Error::OutOfRange { .. }));
    }

    #[test]
    fn test_check_range_above_max() {
        let r = check_range("r", 11.0, 0.0, 10.0);
        assert!(r.is_err());
    }

    #[test]
    fn test_check_len_ok() {
        assert!(check_len(5, 5).is_ok());
    }

    #[test]
    fn test_check_len_mismatch() {
        assert!(check_len(5, 4).is_err());
    }

    #[test]
    fn test_check_index_ok() {
        assert!(check_index(0, 1).is_ok());
        assert!(check_index(4, 5).is_ok());
    }

    #[test]
    fn test_check_index_equal_to_len_fails() {
        assert!(check_index(5, 5).is_err());
    }

    #[test]
    fn test_check_min_points_ok() {
        assert!(check_min_points(3, 3).is_ok());
        assert!(check_min_points(3, 10).is_ok());
    }

    #[test]
    fn test_check_min_points_too_few() {
        assert!(check_min_points(4, 2).is_err());
    }

    #[test]
    fn test_check_positive_ok() {
        assert!(check_positive("s", 0.001).is_ok());
    }

    #[test]
    fn test_check_positive_zero_fails() {
        assert!(check_positive("s", 0.0).is_err());
    }

    #[test]
    fn test_check_non_negative_ok() {
        assert!(check_non_negative("v", 0.0).is_ok());
        assert!(check_non_negative("v", 1.5).is_ok());
    }

    #[test]
    fn test_check_non_negative_negative_fails() {
        assert!(check_non_negative("v", -0.1).is_err());
    }

    #[test]
    fn test_check_finite_ok() {
        assert!(check_finite("x", 3.125).is_ok());
    }

    #[test]
    fn test_check_finite_nan_fails() {
        assert!(check_finite("x", f64::NAN).is_err());
    }

    #[test]
    fn test_check_finite_inf_fails() {
        assert!(check_finite("x", f64::INFINITY).is_err());
    }

    #[test]
    fn test_check_finite_slice_ok() {
        assert!(check_finite_slice("pts", &[1.0, 2.0, 3.0]).is_ok());
    }

    #[test]
    fn test_check_finite_slice_nan_fails() {
        assert!(check_finite_slice("pts", &[1.0, f64::NAN, 3.0]).is_err());
    }

    #[test]
    fn test_check_dim_match_ok() {
        assert!(check_dim_match("pos", 5, "nrm", 5).is_ok());
    }

    #[test]
    fn test_check_dim_match_fail() {
        let r = check_dim_match("pos", 5, "nrm", 3);
        assert!(r.is_err());
        assert!(matches!(r.unwrap_err(), Error::DimensionMismatch { .. }));
    }

    // ── Mesh validation ───────────────────────────────────────────────────────

    #[test]
    fn test_validate_mesh_ok() {
        let verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = [[0, 1, 2]];
        assert!(validate_mesh(&verts, &tris).is_ok());
    }

    #[test]
    fn test_validate_mesh_index_out_of_range() {
        let verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let tris = [[0, 1, 5]]; // index 5 >= len 2
        assert!(validate_mesh(&verts, &tris).is_err());
    }

    #[test]
    fn test_validate_mesh_degenerate_triangle() {
        let verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = [[0, 0, 2]]; // repeated index
        let r = validate_mesh(&verts, &tris);
        assert!(r.is_err());
        assert!(matches!(r.unwrap_err(), Error::DegenerateGeometry { .. }));
    }

    #[test]
    fn test_validate_mesh_empty_vertices() {
        let r = validate_mesh(&[], &[[0, 1, 2]]);
        assert!(r.is_err());
    }

    // ── HeightField validation ────────────────────────────────────────────────

    #[test]
    fn test_validate_heightfield_ok() {
        let heights = vec![0.0f64; 4 * 4];
        assert!(validate_heightfield(&heights, 4, 4, 1.0, 1.0).is_ok());
    }

    #[test]
    fn test_validate_heightfield_too_few_rows() {
        let heights = vec![0.0f64; 4];
        assert!(validate_heightfield(&heights, 1, 4, 1.0, 1.0).is_err());
    }

    #[test]
    fn test_validate_heightfield_bad_scale() {
        let heights = vec![0.0f64; 4 * 4];
        assert!(validate_heightfield(&heights, 4, 4, 0.0, 1.0).is_err());
    }

    #[test]
    fn test_validate_heightfield_len_mismatch() {
        let heights = vec![0.0f64; 10]; // wrong length
        assert!(validate_heightfield(&heights, 4, 4, 1.0, 1.0).is_err());
    }

    #[test]
    fn test_validate_heightfield_nan_height() {
        let mut heights = vec![0.0f64; 4 * 4];
        heights[5] = f64::NAN;
        assert!(validate_heightfield(&heights, 4, 4, 1.0, 1.0).is_err());
    }

    // ── Ray validation ────────────────────────────────────────────────────────

    #[test]
    fn test_validate_ray_dir_ok() {
        assert!(validate_ray_dir([0.0, -1.0, 0.0]).is_ok());
    }

    #[test]
    fn test_validate_ray_dir_zero_fails() {
        assert!(validate_ray_dir([0.0, 0.0, 0.0]).is_err());
    }

    #[test]
    fn test_validate_ray_dir_nan_fails() {
        assert!(validate_ray_dir([f64::NAN, 0.0, 0.0]).is_err());
    }

    // ── Point cloud validation ────────────────────────────────────────────────

    #[test]
    fn test_validate_point_cloud_ok() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        assert!(validate_point_cloud(&pts).is_ok());
    }

    #[test]
    fn test_validate_point_cloud_empty_fails() {
        assert!(validate_point_cloud(&[]).is_err());
    }

    #[test]
    fn test_validate_point_cloud_nan_fails() {
        let pts = [[f64::NAN, 0.0, 0.0]];
        assert!(validate_point_cloud(&pts).is_err());
    }

    // ── WithContext ───────────────────────────────────────────────────────────

    #[test]
    fn test_with_context_ok_passes_through() {
        let r: std::result::Result<i32, &str> = Ok(42);
        let r2: Result<i32> = r.with_context("test");
        assert_eq!(r2.unwrap(), 42);
    }

    #[test]
    fn test_with_context_wraps_error() {
        let r: std::result::Result<i32, &str> = Err("original error");
        let r2: Result<i32> = r.with_context("loading mesh");
        let e = r2.unwrap_err();
        let s = e.to_string();
        assert!(s.contains("loading mesh"), "context missing: {s}");
        assert!(s.contains("original error"), "original missing: {s}");
    }

    // ── ConvergenceTracker ────────────────────────────────────────────────────

    #[test]
    fn test_convergence_tracker_converges() {
        let mut tracker = ConvergenceTracker::new("test_op", 100, 1e-6);
        // First update with small residual should converge
        let result = tracker.update(1e-8);
        assert!(result.is_ok());
        assert!(result.unwrap(), "should report converged");
        assert_eq!(tracker.iterations(), 1);
    }

    #[test]
    fn test_convergence_tracker_not_yet_converged() {
        let mut tracker = ConvergenceTracker::new("test_op", 100, 1e-6);
        let result = tracker.update(0.5);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "should not report converged yet");
    }

    #[test]
    fn test_convergence_tracker_failure() {
        let mut tracker = ConvergenceTracker::new("mesh_smooth", 3, 1e-10);
        let _ = tracker.update(1.0);
        let _ = tracker.update(0.5);
        let r = tracker.update(0.3); // 3rd update exceeds limit
        assert!(r.is_err());
        assert!(matches!(r.unwrap_err(), Error::ConvergenceFailure { .. }));
    }

    #[test]
    fn test_convergence_tracker_residual_tracked() {
        let mut tracker = ConvergenceTracker::new("op", 100, 1e-6);
        let _ = tracker.update(0.7);
        assert!((tracker.residual() - 0.7).abs() < 1e-12);
    }

    #[test]
    fn test_convergence_tracker_reset() {
        let mut tracker = ConvergenceTracker::new("op", 100, 1e-6);
        let _ = tracker.update(0.5);
        tracker.reset();
        assert_eq!(tracker.iterations(), 0);
        assert!(tracker.residual().is_infinite());
    }

    #[test]
    fn test_convergence_tracker_iterations_counted() {
        let mut tracker = ConvergenceTracker::new("op", 100, 1e-6);
        for _ in 0..5 {
            let _ = tracker.update(1.0);
        }
        assert_eq!(tracker.iterations(), 5);
    }
}
