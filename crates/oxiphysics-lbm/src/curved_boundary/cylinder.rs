// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rotating cylinder used as a moving curved wall.

use super::CurvedBcMethod;
use super::lattice::CurvedLattice;
use super::set::{CurvedBoundarySet, LatticeDomain};
use crate::error::Result;

/// A cylinder rotating about its axis (parallel to z) with angular velocity
/// `omega_rot` (radians per lattice time step, positive = counter-clockwise).
#[derive(Debug, Clone, PartialEq)]
pub struct RotatingCylinder {
    /// Centre of the cylinder in the xy-plane (lattice units).
    pub center: [f64; 2],
    /// Cylinder radius (lattice units).
    pub radius: f64,
    /// Angular velocity (rad per lattice time step).
    pub omega_rot: f64,
}

impl RotatingCylinder {
    /// Create a new `RotatingCylinder`.
    pub fn new(center: [f64; 2], radius: f64, omega_rot: f64) -> Self {
        Self {
            center,
            radius,
            omega_rot,
        }
    }

    /// Rigid-body velocity `ω ẑ × (x − x_c)` at `pos`:
    ///
    /// ```text
    /// v_x = -omega_rot * (y - cy)
    /// v_y =  omega_rot * (x - cx)
    /// ```
    pub fn wall_velocity(&self, pos: [f64; 2]) -> [f64; 3] {
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        [-self.omega_rot * dy, self.omega_rot * dx, 0.0]
    }

    /// Whether `pos` lies inside (or on) the cylinder.
    pub fn contains(&self, pos: [f64; 2]) -> bool {
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        dx * dx + dy * dy <= self.radius * self.radius
    }

    /// Signed distance to the cylinder surface: positive in the fluid
    /// (outside), negative inside.
    pub fn signed_distance(&self, pos: [f64; 2]) -> f64 {
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        (dx * dx + dy * dy).sqrt() - self.radius
    }

    /// Parameter `t ∈ (0, 1]` of the first intersection of the segment
    /// `from + t·dir` with the cylinder surface, `None` if there is none.
    ///
    /// Passing a lattice link vector as `dir` yields the link fraction `q`.
    pub fn ray_q(&self, from: [f64; 2], dir: [f64; 2]) -> Option<f64> {
        let ox = from[0] - self.center[0];
        let oy = from[1] - self.center[1];
        let a = dir[0] * dir[0] + dir[1] * dir[1];
        if a <= 0.0 {
            return None;
        }
        let b = 2.0 * (ox * dir[0] + oy * dir[1]);
        let c = ox * ox + oy * oy - self.radius * self.radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let t = (-b - disc.sqrt()) / (2.0 * a);
        if t > 0.0 && t <= 1.0 { Some(t) } else { None }
    }

    /// Wall-cut D2Q9 links of the cylinder on an `nx × ny` lattice with unit
    /// spacing, using the linear Bouzidi scheme.
    ///
    /// Fluid nodes lie outside the cylinder, solid nodes inside; every link
    /// carries the rigid-body wall velocity at its intersection point.
    ///
    /// # Errors
    /// Returns an error if `nx` or `ny` is zero.
    pub fn boundary_set(&self, nx: usize, ny: usize) -> Result<CurvedBoundarySet> {
        let domain = LatticeDomain::new(nx, ny, 1, 1.0)?;
        Ok(self.boundary_set_with(&domain, &CurvedLattice::d2q9(), CurvedBcMethod::Bouzidi))
    }

    /// Wall-cut links of the cylinder on an arbitrary domain and velocity set.
    pub fn boundary_set_with(
        &self,
        domain: &LatticeDomain,
        lattice: &CurvedLattice,
        method: CurvedBcMethod,
    ) -> CurvedBoundarySet {
        let sdf = |pos: [f64; 3]| self.signed_distance([pos[0], pos[1]]);
        let wall = |pos: [f64; 3]| self.wall_velocity([pos[0], pos[1]]);
        CurvedBoundarySet::from_geometry_with(&sdf, &wall, domain, lattice, method)
    }
}
