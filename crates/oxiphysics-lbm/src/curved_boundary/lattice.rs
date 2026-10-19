// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Velocity-set descriptor used by the curved boundary schemes.

use crate::error::{Error, Result};
use crate::lattice::{
    CS2, D2Q9_OPPOSITES, D2Q9_VELOCITIES, D2Q9_WEIGHTS, D3Q19_OPPOSITES, D3Q19_VELOCITIES,
    D3Q19_WEIGHTS, D3Q27_OPPOSITES, D3Q27_VELOCITIES, D3Q27_WEIGHTS,
};

/// Discrete velocity set (velocities, weights and opposite-direction table).
///
/// Two-dimensional sets store `c_z = 0` so that every scheme can work with
/// three-component vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvedLattice {
    velocities: Vec<[i32; 3]>,
    weights: Vec<f64>,
    opposites: Vec<usize>,
    dimensions: usize,
}

impl CurvedLattice {
    /// Build a velocity set after checking its internal consistency.
    ///
    /// Requirements: equal table lengths, weights summing to one, the opposite
    /// table being an involution with `c_opp(i) = −c_i`, and `dimensions` in
    /// `{2, 3}` (2-D sets must have `c_z = 0`).
    pub fn new(
        velocities: Vec<[i32; 3]>,
        weights: Vec<f64>,
        opposites: Vec<usize>,
        dimensions: usize,
    ) -> Result<Self> {
        let n = velocities.len();
        if n == 0 || weights.len() != n || opposites.len() != n {
            return Err(Error::ComponentCountMismatch {
                expected: n,
                found: weights.len().min(opposites.len()),
                field: "curved lattice tables".to_string(),
            });
        }
        if !(2..=3).contains(&dimensions) {
            return Err(Error::UnsupportedVelocitySet {
                name: format!("{dimensions}-dimensional velocity set"),
            });
        }
        let weight_sum: f64 = weights.iter().sum();
        if (weight_sum - 1.0).abs() > 1e-12 || weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(Error::UnsupportedVelocitySet {
                name: format!("weights must be positive and sum to 1 (sum = {weight_sum})"),
            });
        }
        for (i, (&opp, c)) in opposites.iter().zip(velocities.iter()).enumerate() {
            let c_opp = velocities
                .get(opp)
                .ok_or_else(|| Error::UnsupportedVelocitySet {
                    name: format!("opposite of direction {i} out of range ({opp})"),
                })?;
            let back = opposites.get(opp).copied();
            if back != Some(i) || (0..3).any(|d| c_opp[d] != -c[d]) {
                return Err(Error::UnsupportedVelocitySet {
                    name: format!("inconsistent opposite table at direction {i}"),
                });
            }
            if dimensions == 2 && c[2] != 0 {
                return Err(Error::UnsupportedVelocitySet {
                    name: format!("2-D velocity set has c_z != 0 at direction {i}"),
                });
            }
        }
        Ok(Self {
            velocities,
            weights,
            opposites,
            dimensions,
        })
    }

    /// The D2Q9 set (order: rest, E, N, W, S, NE, NW, SW, SE).
    pub fn d2q9() -> Self {
        Self {
            velocities: D2Q9_VELOCITIES.iter().map(|c| [c[0], c[1], 0]).collect(),
            weights: D2Q9_WEIGHTS.to_vec(),
            opposites: D2Q9_OPPOSITES.to_vec(),
            dimensions: 2,
        }
    }

    /// The D3Q19 set (same ordering as [`crate::lattice::D3Q19_VELOCITIES`]).
    pub fn d3q19() -> Self {
        Self {
            velocities: D3Q19_VELOCITIES.to_vec(),
            weights: D3Q19_WEIGHTS.to_vec(),
            opposites: D3Q19_OPPOSITES.to_vec(),
            dimensions: 3,
        }
    }

    /// The D3Q27 set (same ordering as [`crate::lattice::D3Q27_VELOCITIES`]).
    pub fn d3q27() -> Self {
        Self {
            velocities: D3Q27_VELOCITIES.to_vec(),
            weights: D3Q27_WEIGHTS.to_vec(),
            opposites: D3Q27_OPPOSITES.to_vec(),
            dimensions: 3,
        }
    }

    /// Number of discrete velocities `Q`.
    pub fn num_directions(&self) -> usize {
        self.velocities.len()
    }

    /// Spatial dimension (2 or 3).
    pub fn dimensions(&self) -> usize {
        self.dimensions
    }

    /// Integer velocity vectors.
    pub fn velocities(&self) -> &[[i32; 3]] {
        &self.velocities
    }

    /// Lattice weights.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Opposite-direction table.
    pub fn opposites(&self) -> &[usize] {
        &self.opposites
    }

    /// Velocity `c_i` as floating-point vector, `None` if `i` is out of range.
    pub fn velocity_f64(&self, i: usize) -> Option<[f64; 3]> {
        self.velocities
            .get(i)
            .map(|c| [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])])
    }

    /// Weight `w_i`, `None` if `i` is out of range.
    pub fn weight(&self, i: usize) -> Option<f64> {
        self.weights.get(i).copied()
    }

    /// Opposite direction `ī`, `None` if `i` is out of range.
    pub fn opposite(&self, i: usize) -> Option<usize> {
        self.opposites.get(i).copied()
    }

    /// Index of the direction whose velocity equals `c`, if any.
    pub fn direction_of(&self, c: [i32; 3]) -> Option<usize> {
        self.velocities.iter().position(|v| *v == c)
    }

    /// Zeroth and first moments `(ρ, Σ c_i f_i)` of a node's populations.
    ///
    /// Returns an error if `f` does not have exactly `Q` entries.
    pub fn moments(&self, f: &[f64]) -> Result<(f64, [f64; 3])> {
        if f.len() != self.velocities.len() {
            return Err(Error::ComponentCountMismatch {
                expected: self.velocities.len(),
                found: f.len(),
                field: "node populations".to_string(),
            });
        }
        let mut rho = 0.0;
        let mut j = [0.0; 3];
        for (fi, c) in f.iter().zip(self.velocities.iter()) {
            rho += fi;
            for d in 0..3 {
                j[d] += f64::from(c[d]) * fi;
            }
        }
        Ok((rho, j))
    }

    /// Second-order equilibrium `f_i^eq(ρ, u)`; returns `0` for an invalid `i`.
    pub fn equilibrium(&self, i: usize, rho: f64, u: [f64; 3]) -> f64 {
        match (self.velocity_f64(i), self.weight(i)) {
            (Some(c), Some(w)) => {
                let cu = dot3(c, u);
                let uu = dot3(u, u);
                w * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - uu / (2.0 * CS2))
            }
            _ => 0.0,
        }
    }
}

/// Three-component dot product.
pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
