// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics ik module.
//!
//! Exposes FABRIK and 2-bone analytic IK to Python.

use oxiphysics::ik::{IkChain, IkSolver};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyIkChain
// ─────────────────────────────────────────────────────────────────────────────

/// An IK chain (root position + joint segments).
#[pyclass(name = "IkChain")]
pub struct PyIkChain {
    inner: IkChain,
}

#[pymethods]
impl PyIkChain {
    /// Create a new IK chain from the given root position and segments.
    ///
    /// `segments` is a list of `(length, min_angle_deg, max_angle_deg)`.
    #[new]
    pub fn new(root: [f64; 3], segments: Vec<(f64, f64, f64)>) -> Self {
        Self {
            inner: IkChain::new(root, &segments),
        }
    }

    /// Total maximum reach of the chain.
    pub fn total_reach(&self) -> f64 {
        self.inner.total_reach()
    }

    /// Number of segments in the chain.
    pub fn segment_count(&self) -> usize {
        self.inner.joints.len()
    }

    /// Compute world-space joint positions from root and local offsets.
    fn joint_world_positions(&self) -> Vec<[f64; 3]> {
        let mut positions = Vec::with_capacity(self.inner.joints.len() + 1);
        positions.push(self.inner.root_world);
        let mut prev = self.inner.root_world;
        for joint in &self.inner.joints {
            let next = [
                prev[0] + joint.local_offset[0],
                prev[1] + joint.local_offset[1],
                prev[2] + joint.local_offset[2],
            ];
            positions.push(next);
            prev = next;
        }
        positions
    }

    /// Solve FABRIK IK to reach the target and return the joint positions as JSON.
    pub fn solve_fabrik_json(&mut self, target: [f64; 3]) -> PyResult<String> {
        let solver = IkSolver::Fabrik {
            max_iterations: 20,
            tolerance: 1e-4,
        };
        let report = solver.solve(&mut self.inner, target);
        let positions = self.joint_world_positions();
        let result = serde_json::json!({
            "converged": report.converged,
            "iterations": report.iterations,
            "residual": report.residual,
            "positions": positions,
        });
        serde_json::to_string(&result)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Solve 2-bone analytic IK and return joint positions as JSON.
    pub fn solve_two_bone_json(&mut self, target: [f64; 3]) -> PyResult<String> {
        let report = IkSolver::TwoBone.solve(&mut self.inner, target);
        let positions = self.joint_world_positions();
        let result = serde_json::json!({
            "converged": report.converged,
            "iterations": report.iterations,
            "residual": report.residual,
            "positions": positions,
        });
        serde_json::to_string(&result)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyIkChain>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ik_chain_instantiation() {
        let chain = PyIkChain::new(
            [0.0, 0.0, 0.0],
            vec![(1.0, -90.0, 90.0), (1.0, -90.0, 90.0)],
        );
        assert_eq!(chain.segment_count(), 2);
        assert!((chain.total_reach() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_ik_fabrik_solve() {
        let mut chain = PyIkChain::new(
            [0.0, 0.0, 0.0],
            vec![(1.0, -90.0, 90.0), (1.0, -90.0, 90.0)],
        );
        let json = chain
            .solve_fabrik_json([1.0, 1.0, 0.0])
            .expect("solve failed");
        assert!(json.contains("positions"));
    }
}
