// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Equilibrium (PDE) residual of a stress field on a regular grid.
//!
//! The static equilibrium equation of a continuum with body force density `f`
//! is `div(sigma) + f = 0`.  For a stress field sampled on a regular grid of
//! spacing `h` (flat index `i + nx * (j + ny * k)`, x fastest) the divergence
//! at an interior node is evaluated with second-order central differences,
//!
//! ```text
//! (div sigma)_i = sum_j [ sigma_ij(x + h e_j) - sigma_ij(x - h e_j) ] / (2 h),
//! ```
//!
//! and the physics-informed loss is the mean squared residual
//! `L_pde = mean over interior nodes of || div sigma + f ||^2`.
//!
//! A grid axis with a single node (`n = 1`) describes a field that does not
//! vary along that axis (plane problems); such an axis contributes no
//! derivative and all of its nodes are interior along it.  An axis needs at
//! least three nodes otherwise.

use std::ops::Range;

use super::PhysicsInformedLoss;
use crate::{Error, Result};

/// Voigt index of the stress component `sigma_ij` (`[xx, yy, zz, xy, xz, yz]`).
const VOIGT_INDEX: [[usize; 3]; 3] = [[0, 3, 4], [3, 1, 5], [4, 5, 2]];

/// Voigt index of the stress component `sigma_(row, col)`.
fn voigt_index(row: usize, col: usize) -> usize {
    VOIGT_INDEX[row][col]
}

/// Validated grid description shared by the PDE operators.
struct Grid {
    dims: [usize; 3],
    /// Interior node ranges along each axis.
    interior: [Range<usize>; 3],
    /// Grid spacing `h`.
    spacing: f64,
}

impl Grid {
    /// Validate `dims`, the stress field length and the grid spacing.
    fn new(node_count: usize, dims: [usize; 3], spacing: f64) -> Result<Self> {
        let expected = dims[0]
            .checked_mul(dims[1])
            .and_then(|n| n.checked_mul(dims[2]))
            .ok_or_else(|| {
                Error::InvalidParameter(format!("grid dimensions {dims:?} overflow usize"))
            })?;
        if expected == 0 || node_count != expected {
            return Err(Error::InvalidParameter(format!(
                "stress field has {node_count} nodes but grid {dims:?} requires {expected} (and at least one node per axis)"
            )));
        }
        if !(spacing.is_finite() && spacing > 0.0) {
            return Err(Error::InvalidParameter(format!(
                "grid spacing (fd_step) must be finite and positive, got {spacing}"
            )));
        }
        let mut interior = [0..0, 0..0, 0..0];
        for (axis, range) in interior.iter_mut().enumerate() {
            *range = match dims[axis] {
                1 => 0..1,
                2 => {
                    return Err(Error::InvalidParameter(format!(
                        "grid {dims:?} has no interior node along axis {axis}: use 1 or at least 3 nodes"
                    )));
                }
                n => 1..n - 1,
            };
        }
        Ok(Self {
            dims,
            interior,
            spacing,
        })
    }

    fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.dims[0] * (j + self.dims[1] * k)
    }

    /// Number of interior nodes.
    fn interior_count(&self) -> usize {
        self.interior.iter().map(|r| r.len()).product()
    }

    /// Visit every interior node `(i, j, k)` in x-fastest order.
    fn for_each_interior(&self, mut visit: impl FnMut([usize; 3])) {
        for k in self.interior[2].clone() {
            for j in self.interior[1].clone() {
                for i in self.interior[0].clone() {
                    visit([i, j, k]);
                }
            }
        }
    }

    /// Flat indices of the neighbours of `node` along `axis`, or `None` if the
    /// field does not vary along that axis.
    fn neighbours(&self, node: [usize; 3], axis: usize) -> Option<(usize, usize)> {
        if self.dims[axis] < 3 {
            return None;
        }
        let mut plus = node;
        let mut minus = node;
        plus[axis] += 1;
        minus[axis] -= 1;
        Some((
            self.index(plus[0], plus[1], plus[2]),
            self.index(minus[0], minus[1], minus[2]),
        ))
    }

    /// Equilibrium residual `div(sigma) + f` at an interior node.
    fn residual(&self, stress: &[[f64; 6]], node: [usize; 3], body_force: [f64; 3]) -> [f64; 3] {
        let mut residual = body_force;
        for axis in 0..3 {
            if let Some((plus, minus)) = self.neighbours(node, axis) {
                for (comp, r) in residual.iter_mut().enumerate() {
                    let voigt = voigt_index(comp, axis);
                    *r += (stress[plus][voigt] - stress[minus][voigt]) / (2.0 * self.spacing);
                }
            }
        }
        residual
    }
}

impl PhysicsInformedLoss {
    /// Equilibrium residuals `div(sigma) + f` at the interior nodes of a
    /// regular grid, in x-fastest order of the interior nodes.
    ///
    /// `stress` holds one Voigt stress `[xx, yy, zz, xy, xz, yz]` per grid
    /// node (`dims = [nx, ny, nz]`, flat index `i + nx * (j + ny * k)`), the
    /// grid spacing is [`fd_step`](PhysicsInformedLoss::fd_step) and
    /// `body_force` is the body force density `f`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `stress.len()` does not equal
    /// `nx * ny * nz`, an axis has exactly two nodes (no interior node), or
    /// `fd_step` is not finite and positive.
    pub fn equilibrium_residuals(
        &self,
        stress: &[[f64; 6]],
        dims: [usize; 3],
        body_force: [f64; 3],
    ) -> Result<Vec<[f64; 3]>> {
        let grid = Grid::new(stress.len(), dims, self.fd_step)?;
        let mut residuals = Vec::with_capacity(grid.interior_count());
        grid.for_each_interior(|node| residuals.push(grid.residual(stress, node, body_force)));
        Ok(residuals)
    }

    /// Physics-informed equilibrium loss
    /// `L_pde = mean over interior nodes of || div(sigma) + f ||^2`.
    ///
    /// A constant stress field with `body_force = 0` has `L_pde = 0`.
    ///
    /// # Errors
    ///
    /// See [`equilibrium_residuals`](Self::equilibrium_residuals).
    pub fn pde_loss(
        &self,
        stress: &[[f64; 6]],
        dims: [usize; 3],
        body_force: [f64; 3],
    ) -> Result<f64> {
        let grid = Grid::new(stress.len(), dims, self.fd_step)?;
        let mut sum = 0.0;
        grid.for_each_interior(|node| {
            let r = grid.residual(stress, node, body_force);
            sum += r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        });
        Ok(sum / grid.interior_count() as f64)
    }

    /// Gradient `d L_pde / d sigma` of [`pde_loss`](Self::pde_loss) with
    /// respect to every Voigt stress component of every grid node.
    ///
    /// A shear component `sigma_xy` enters the residual through both
    /// `sigma_xy` and `sigma_yx`; the returned gradient is the derivative with
    /// respect to the single stored value.  Multiplied by
    /// [`lambda_pde`](PhysicsInformedLoss::lambda_pde) it is the PDE term of
    /// the loss gradient that [`NeuralPinnSolver::train_step`](super::NeuralPinnSolver::train_step)
    /// backpropagates.
    ///
    /// # Errors
    ///
    /// See [`equilibrium_residuals`](Self::equilibrium_residuals).
    pub fn pde_gradient(
        &self,
        stress: &[[f64; 6]],
        dims: [usize; 3],
        body_force: [f64; 3],
    ) -> Result<Vec<[f64; 6]>> {
        let grid = Grid::new(stress.len(), dims, self.fd_step)?;
        let mut gradient = vec![[0.0f64; 6]; stress.len()];
        // d L / d r = 2 r / N_interior and d r / d sigma_plus/minus = +-1 / (2 h).
        let scale = 2.0 / grid.interior_count() as f64 / (2.0 * grid.spacing);
        grid.for_each_interior(|node| {
            let r = grid.residual(stress, node, body_force);
            for axis in 0..3 {
                if let Some((plus, minus)) = grid.neighbours(node, axis) {
                    for (comp, &r_comp) in r.iter().enumerate() {
                        let voigt = voigt_index(comp, axis);
                        gradient[plus][voigt] += scale * r_comp;
                        gradient[minus][voigt] -= scale * r_comp;
                    }
                }
            }
        });
        Ok(gradient)
    }

    /// Total physics-informed loss `L = L_data + lambda_pde * L_pde` of a
    /// predicted stress field on a regular grid.
    ///
    /// `L_data` is the mean squared error against `target` over all nodes and
    /// stress components, `L_pde` the equilibrium residual of `stress` (see
    /// [`pde_loss`](Self::pde_loss)) for the body force density `body_force`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidParameter`] if `target` and `stress` differ in
    /// length, or for the reasons listed at
    /// [`equilibrium_residuals`](Self::equilibrium_residuals).
    pub fn compute_field(
        &self,
        stress: &[[f64; 6]],
        dims: [usize; 3],
        target: &[[f64; 6]],
        body_force: [f64; 3],
    ) -> Result<f64> {
        if stress.len() != target.len() {
            return Err(Error::InvalidParameter(format!(
                "target field has {} nodes but the predicted field has {}",
                target.len(),
                stress.len()
            )));
        }
        let pde = self.pde_loss(stress, dims, body_force)?;
        let flat_pred: Vec<f64> = stress.iter().flatten().copied().collect();
        let flat_target: Vec<f64> = target.iter().flatten().copied().collect();
        Ok(self.data_loss(&flat_pred, &flat_target) + self.lambda_pde * pde)
    }
}
