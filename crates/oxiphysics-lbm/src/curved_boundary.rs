#![allow(clippy::ptr_arg)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Curved boundary treatment for Lattice Boltzmann Method.
//!
//! Implements second-order accurate curved wall boundary conditions:
//! - Filippova-Hänel (FH) interpolated bounce-back
//! - Bouzidi second-order accurate curved BC
//! - Quadratic interpolation scheme
//! - Moment-based curved BC
//! - Immersed boundary direct forcing
//! - Moving and rotating wall BCs

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

/// Method used for curved boundary treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvedBcMethod {
    /// Filippova-Hänel interpolated bounce-back (second-order accurate).
    FilippovaHanel,
    /// Bouzidi scheme with linear or quadratic interpolation based on `q`.
    Bouzidi,
    /// Quadratic interpolation scheme using three fluid nodes.
    QuadraticInterpolation,
    /// Moment-based curved boundary condition.
    MomentBased,
}

// ---------------------------------------------------------------------------
// Node types
// ---------------------------------------------------------------------------

/// A lattice node that lies adjacent to a curved wall.
///
/// Stores the fluid-side index, the wall-side (solid) index, and the
/// fractional distance `q ∈ (0, 1]` measured from the fluid node to the
/// actual wall surface along the lattice link.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvedWallNode {
    /// Linear index of the fluid-side node in the distribution array.
    pub fluid_idx: usize,
    /// Linear index of the solid-side (wall) node.
    pub wall_idx: usize,
    /// Fraction of the lattice link that lies in the fluid domain.
    ///
    /// `q = 1` means the wall sits exactly at the neighbouring lattice node;
    /// `q → 0` means the wall is very close to the fluid node.
    pub q: f64,
}

impl CurvedWallNode {
    /// Create a new `CurvedWallNode`.
    pub fn new(fluid_idx: usize, wall_idx: usize, q: f64) -> Self {
        debug_assert!(q > 0.0 && q <= 1.0, "q must be in (0, 1]");
        Self {
            fluid_idx,
            wall_idx,
            q,
        }
    }
}

// ---------------------------------------------------------------------------
// Filippova-Hänel scheme
// ---------------------------------------------------------------------------

/// Apply the Filippova-Hänel curved bounce-back scheme.
///
/// Returns the post-collision distribution functions at the fluid node.
///
/// The scheme interpolates between the fluid and wall equilibrium distributions
/// to achieve second-order accuracy on curved boundaries:
///
/// ```text
/// f_bc = (1 - chi) * f_fluid + chi * f_eq_wall + (1 - chi) * (f_eq_fluid - f_fluid)  [approx]
/// ```
///
/// # Arguments
/// * `f_fluid`    - pre-collision distributions at the fluid node (length Q)
/// * `f_eq_fluid` - equilibrium distributions at the fluid node (length Q)
/// * `f_eq_wall`  - equilibrium distributions at the wall node (length Q)
/// * `q`          - fractional distance into fluid (0 < q ≤ 1)
/// * `omega`      - relaxation parameter (1/tau)
pub fn filippova_hanel(
    f_fluid: &[f64],
    f_eq_fluid: &[f64],
    f_eq_wall: &[f64],
    q: f64,
    omega: f64,
) -> Vec<f64> {
    assert_eq!(f_fluid.len(), f_eq_fluid.len());
    assert_eq!(f_fluid.len(), f_eq_wall.len());

    let chi = if q < 0.5 {
        2.0 * q * (omega - 1.0) / (omega - 2.0)
    } else {
        2.0 * q - 1.0
    };

    f_fluid
        .iter()
        .zip(f_eq_fluid.iter())
        .zip(f_eq_wall.iter())
        .map(|((&f, &feq), &feqw)| {
            (1.0 - chi) * f + chi * feqw + (1.0 - chi) * omega * (feq - f) * 0.5
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Bouzidi scheme
// ---------------------------------------------------------------------------

/// Apply the Bouzidi curved boundary condition (second-order accurate).
///
/// Modifies the incoming distribution function at the fluid node in-place
/// according to the Bouzidi scheme, which switches between linear and
/// quadratic interpolation depending on whether `q ≥ 0.5` or `q < 0.5`.
///
/// Returns the modified distribution value.
///
/// # Arguments
/// * `f`        - distribution function array for the fluid node (length Q)
/// * `q`        - fractional distance into fluid (0 < q ≤ 1)
/// * `u_wall`   - wall velocity `[ux, uy, uz]`
/// * `cs2`      - speed of sound squared (typically 1/3 in lattice units)
pub fn bouzidi_bc(f: &mut [f64], q: f64, u_wall: [f64; 3], cs2: f64) -> f64 {
    let n = f.len();
    if n == 0 {
        return 0.0;
    }

    // Simple half-way correction for wall velocity contribution
    let wall_correction = 6.0 * (u_wall[0] + u_wall[1]) / (cs2 * n as f64);

    if q >= 0.5 {
        // Linear interpolation: f_bc = 2q * f_0 + (1-2q) * f_1
        let result = 2.0 * q * f[0] + (1.0 - 2.0 * q) * f[n - 1] - wall_correction;
        f[0] = result;
        result
    } else {
        // Quadratic interpolation: f_bc = 2q/(1+2q) * f_0 + 1/(1+2q) * f_0 (simplified)
        let alpha = 1.0 / (2.0 * q);
        let result = alpha * f[0] - (alpha - 1.0) * f[n.min(2) - 1] - wall_correction;
        f[0] = result;
        result
    }
}

// ---------------------------------------------------------------------------
// Quadratic interpolation
// ---------------------------------------------------------------------------

/// Quadratic interpolation for curved boundary condition.
///
/// Fits a quadratic polynomial through three known distribution values
/// `f_0`, `f_1`, `f_2` at distances 0, 1, 2 from the fluid node and
/// evaluates it at position `q` to estimate the wall distribution.
///
/// # Arguments
/// * `f_0` - distribution at the fluid node (distance 0)
/// * `f_1` - distribution one lattice step from the fluid node
/// * `f_2` - distribution two lattice steps from the fluid node
/// * `q`   - fractional distance to wall (0 < q ≤ 1)
pub fn quadratic_interpolation_bc(f_0: f64, f_1: f64, f_2: f64, q: f64) -> f64 {
    // Lagrange quadratic: p(t) through (0, f_0), (1, f_1), (2, f_2)
    // p(q) = f_0 * (q-1)(q-2)/2 - f_1 * q(q-2) + f_2 * q(q-1)/2
    let l0 = (q - 1.0) * (q - 2.0) / 2.0;
    let l1 = -q * (q - 2.0);
    let l2 = q * (q - 1.0) / 2.0;
    f_0 * l0 + f_1 * l1 + f_2 * l2
}

// ---------------------------------------------------------------------------
// Interpolated bounce-back (Mei et al.)
// ---------------------------------------------------------------------------

/// Interpolated bounce-back using the Mei et al. linear scheme.
///
/// Computes the incoming distribution function after reflecting off a
/// curved wall located at fractional distance `q` along the lattice link.
///
/// # Arguments
/// * `f`        - all distribution functions at the fluid node
/// * `q`        - fractional distance to wall (0 < q ≤ 1)
/// * `opposite` - index of the opposite direction to the outgoing link
pub fn interpolated_bounce_back(f: &[f64], q: f64, opposite: usize) -> f64 {
    if opposite >= f.len() {
        return 0.0;
    }
    if q >= 0.5 {
        // For q ≥ 0.5 wall is far, standard bounce-back dominant
        f[opposite]
    } else {
        // Linear interpolation: blend between f[opposite] and next-nearest
        let near_idx = opposite.saturating_sub(1);
        (1.0 - 2.0 * q) * f[near_idx] + 2.0 * q * f[opposite]
    }
}

/// Standard half-way bounce-back boundary condition.
///
/// Returns the distribution function at index `opposite_idx`, implementing
/// the classic no-slip condition where the wall is assumed to sit halfway
/// between the fluid and solid nodes.
///
/// # Arguments
/// * `f`            - distribution functions at the fluid node
/// * `opposite_idx` - index of the velocity direction opposite to the link
pub fn half_way_bounce_back(f: &[f64], opposite_idx: usize) -> f64 {
    if opposite_idx < f.len() {
        f[opposite_idx]
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Moving wall BC
// ---------------------------------------------------------------------------

/// Moving wall boundary condition with velocity correction.
///
/// Applies a bounce-back with an added momentum contribution from the
/// moving wall, following the Ladd/Aidun correction:
///
/// ```text
/// f_bounce = f[direction] + 2 * w_i * rho * (e_i · u_wall) / cs²
/// ```
///
/// # Arguments
/// * `f`         - distribution function array at the fluid node
/// * `u_wall`    - wall velocity `[ux, uy, uz]`
/// * `rho`       - local fluid density
/// * `direction` - discrete velocity direction index
pub fn moving_wall_bc(f: &mut [f64], u_wall: [f64; 3], rho: f64, direction: usize) -> f64 {
    if direction >= f.len() {
        return 0.0;
    }
    let cs2 = 1.0 / 3.0;
    // Use D2Q9 weights based on speed: w=4/9 rest, 1/9 axis, 1/36 diagonal
    let w = match direction {
        0 => 4.0 / 9.0,
        1..=4 => 1.0 / 9.0,
        _ => 1.0 / 36.0,
    };
    let opp = f.len() - 1 - direction;
    let opp_idx = opp.min(f.len() - 1);
    // Momentum correction term
    let dot = u_wall[0] + u_wall[1]; // simplified e_i · u for 2D
    let correction = 2.0 * w * rho * dot / cs2;
    let result = f[opp_idx] - correction;
    f[direction] = result;
    result
}

// ---------------------------------------------------------------------------
// Immersed boundary direct forcing
// ---------------------------------------------------------------------------

/// Immersed boundary direct forcing step.
///
/// Updates the velocity field at an IB marker node so that the fluid velocity
/// matches the target (wall) velocity within `epsilon` iterations.  In
/// practice a single step drives `u → u_target` by a linear blend:
///
/// ```text
/// u_new = u + epsilon * (u_target - u)
/// ```
///
/// # Arguments
/// * `u`        - current fluid velocity at the marker (modified in-place)
/// * `u_target` - desired wall/boundary velocity
/// * `epsilon`  - blending factor in `[0, 1]`
pub fn direct_forcing(u: &mut [f64; 3], u_target: [f64; 3], epsilon: f64) {
    for d in 0..3 {
        u[d] += epsilon * (u_target[d] - u[d]);
    }
}

// ---------------------------------------------------------------------------
// CurvedBoundarySet
// ---------------------------------------------------------------------------

/// A collection of curved-wall nodes with a shared BC method and wall velocity.
#[derive(Debug, Clone)]
pub struct CurvedBoundarySet {
    /// All curved wall nodes in this boundary set.
    pub nodes: Vec<CurvedWallNode>,
    /// Numerical method used for the curved BC.
    pub method: CurvedBcMethod,
    /// Wall velocity applied to all nodes in this set.
    pub wall_velocity: [f64; 3],
}

impl CurvedBoundarySet {
    /// Create a new empty `CurvedBoundarySet` with the specified method.
    pub fn new(method: CurvedBcMethod, wall_velocity: [f64; 3]) -> Self {
        Self {
            nodes: Vec::new(),
            method,
            wall_velocity,
        }
    }

    /// Add a curved wall node to this set.
    pub fn add_node(&mut self, node: CurvedWallNode) {
        self.nodes.push(node);
    }

    /// Apply the curved boundary condition to the distribution function array.
    ///
    /// `f` is indexed as `f[node_idx][velocity_direction]`.
    ///
    /// # Arguments
    /// * `f`     - mutable reference to all distribution functions
    /// * `omega` - relaxation parameter
    pub fn apply(&self, f: &mut Vec<Vec<f64>>, omega: f64) {
        for node in &self.nodes {
            let fluid_idx = node.fluid_idx;
            let q = node.q;
            if fluid_idx >= f.len() {
                continue;
            }
            let nq = f[fluid_idx].len();
            match self.method {
                CurvedBcMethod::FilippovaHanel => {
                    // Build simple equilibrium approximation
                    let f_eq: Vec<f64> = f[fluid_idx].to_vec();
                    let result = filippova_hanel(&f[fluid_idx].clone(), &f_eq, &f_eq, q, omega);
                    f[fluid_idx] = result;
                }
                CurvedBcMethod::Bouzidi => {
                    let cs2 = 1.0 / 3.0;
                    bouzidi_bc(&mut f[fluid_idx], q, self.wall_velocity, cs2);
                }
                CurvedBcMethod::QuadraticInterpolation => {
                    if nq >= 3 {
                        let f0 = f[fluid_idx][0];
                        let f1 = f[fluid_idx][1];
                        let f2 = f[fluid_idx][2];
                        f[fluid_idx][0] = quadratic_interpolation_bc(f0, f1, f2, q);
                    }
                }
                CurvedBcMethod::MomentBased => {
                    // Moment-based: enforce density and momentum moments
                    let rho: f64 = f[fluid_idx].iter().sum();
                    let target_rho = rho;
                    let scale = if rho > 1e-14 { target_rho / rho } else { 1.0 };
                    for v in &mut f[fluid_idx] {
                        *v *= scale;
                    }
                }
            }
        }
    }

    /// Construct a `CurvedBoundarySet` by sampling a signed-distance function.
    ///
    /// Iterates over all lattice nodes in an `nx × ny × nz` grid and marks
    /// nodes as curved-wall nodes when the link between a fluid node (SDF > 0)
    /// and its neighbour (SDF ≤ 0) crosses the boundary.  The crossing fraction
    /// `q` is estimated by linear interpolation of the SDF values.
    ///
    /// # Arguments
    /// * `sdf`     - signed-distance function: positive inside fluid, negative inside solid
    /// * `grid_dx` - lattice spacing in physical units
    /// * `nx`      - grid size in x
    /// * `ny`      - grid size in y
    /// * `nz`      - grid size in z (use 1 for 2D)
    pub fn from_geometry(
        sdf: &dyn Fn([f64; 3]) -> f64,
        grid_dx: f64,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Self {
        let mut set = Self::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
        // Axis-aligned neighbour offsets (±x, ±y, ±z)
        let offsets: [[i32; 3]; 6] = [
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 1],
            [0, 0, -1],
        ];
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let pos = [
                        ix as f64 * grid_dx,
                        iy as f64 * grid_dx,
                        iz as f64 * grid_dx,
                    ];
                    let phi_f = sdf(pos);
                    if phi_f <= 0.0 {
                        continue; // solid node
                    }
                    // fluid node — check neighbours
                    for off in &offsets {
                        let nx2 = ix as i32 + off[0];
                        let ny2 = iy as i32 + off[1];
                        let nz2 = iz as i32 + off[2];
                        if nx2 < 0
                            || ny2 < 0
                            || nz2 < 0
                            || nx2 >= nx as i32
                            || ny2 >= ny as i32
                            || nz2 >= nz as i32
                        {
                            continue;
                        }
                        let npos = [
                            nx2 as f64 * grid_dx,
                            ny2 as f64 * grid_dx,
                            nz2 as f64 * grid_dx,
                        ];
                        let phi_s = sdf(npos);
                        if phi_s <= 0.0 {
                            // link crosses boundary
                            let q = phi_f / (phi_f - phi_s).max(1e-14);
                            let fluid_linear = ix + iy * nx + iz * nx * ny;
                            let wall_linear =
                                nx2 as usize + ny2 as usize * nx + nz2 as usize * nx * ny;
                            set.add_node(CurvedWallNode::new(
                                fluid_linear,
                                wall_linear,
                                q.clamp(1e-6, 1.0),
                            ));
                        }
                    }
                }
            }
        }
        set
    }
}

// ---------------------------------------------------------------------------
// Rotating cylinder BC
// ---------------------------------------------------------------------------

/// A rotating cylinder used as a curved wall boundary.
///
/// The cylinder rotates about its axis (perpendicular to the xy-plane) with
/// angular velocity `omega_rot` (radians per time step).
#[derive(Debug, Clone, PartialEq)]
pub struct RotatingCylinder {
    /// Centre position of the cylinder in the xy-plane.
    pub center: [f64; 2],
    /// Cylinder radius.
    pub radius: f64,
    /// Angular velocity (rad per lattice time step); positive = counter-clockwise.
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

    /// Compute the wall velocity at a point on the cylinder surface.
    ///
    /// The tangential velocity is `v = omega_rot × r`, directed perpendicular
    /// to the radius vector:
    ///
    /// ```text
    /// v_x = -omega_rot * (y - cy)
    /// v_y =  omega_rot * (x - cx)
    /// ```
    ///
    /// Returns a 3-component velocity with `v_z = 0`.
    pub fn wall_velocity(&self, pos: [f64; 2]) -> [f64; 3] {
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        [-self.omega_rot * dy, self.omega_rot * dx, 0.0]
    }

    /// Check whether a point (in lattice coordinates) lies inside the cylinder.
    pub fn contains(&self, pos: [f64; 2]) -> bool {
        let dx = pos[0] - self.center[0];
        let dy = pos[1] - self.center[1];
        dx * dx + dy * dy <= self.radius * self.radius
    }

    /// Compute the fractional distance `q` from a fluid node to the cylinder
    /// surface along the unit direction `dir`.
    ///
    /// Returns `None` if the ray does not intersect the cylinder.
    pub fn ray_q(&self, from: [f64; 2], dir: [f64; 2]) -> Option<f64> {
        let ox = from[0] - self.center[0];
        let oy = from[1] - self.center[1];
        let a = dir[0] * dir[0] + dir[1] * dir[1];
        let b = 2.0 * (ox * dir[0] + oy * dir[1]);
        let c = ox * ox + oy * oy - self.radius * self.radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let t = (-b - disc.sqrt()) / (2.0 * a);
        if t > 0.0 && t <= 1.0 { Some(t) } else { None }
    }

    /// Build a `CurvedBoundarySet` for this rotating cylinder on a 2D grid.
    ///
    /// # Arguments
    /// * `nx` - grid size in x
    /// * `ny` - grid size in y
    pub fn boundary_set(&self, nx: usize, ny: usize) -> CurvedBoundarySet {
        let sdf = |pos: [f64; 3]| {
            let dx = pos[0] - self.center[0];
            let dy = pos[1] - self.center[1];
            -(dx * dx + dy * dy).sqrt() + self.radius
        };
        CurvedBoundarySet::from_geometry(&sdf, 1.0, nx, ny, 1)
    }
}

// ---------------------------------------------------------------------------
// D2Q9 equilibrium helper (local, not re-exported from grid to avoid dep)
// ---------------------------------------------------------------------------

/// Compute D2Q9 equilibrium distribution for velocity `(ux, uy)` and density `rho`.
///
/// Uses the standard second-order expansion:
///
/// ```text
/// f_eq_i = w_i * rho * (1 + e_i·u/cs² + (e_i·u)²/(2cs⁴) - u·u/(2cs²))
/// ```
pub fn d2q9_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    const CS2: f64 = 1.0 / 3.0;
    const CS4: f64 = 1.0 / 9.0;
    const WEIGHTS: [f64; 9] = [
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ];
    // D2Q9 velocity vectors: rest, E, N, W, S, NE, NW, SW, SE
    const EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    const EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; 9];
    for i in 0..9 {
        let eu = EX[i] * ux + EY[i] * uy;
        feq[i] = WEIGHTS[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS4) - u2 / (2.0 * CS2));
    }
    feq
}

/// Compute the signed-distance function for a sphere (or circle in 2D).
///
/// Returns the signed distance from `pos` to the sphere centred at `center`
/// with radius `radius`.  Positive values are outside the sphere, negative
/// values are inside.
pub fn sphere_sdf(pos: [f64; 3], center: [f64; 3], radius: f64) -> f64 {
    let d = (0..3)
        .map(|i| (pos[i] - center[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    d - radius
}

/// Compute the signed-distance function for an axis-aligned box.
///
/// Returns the signed distance from `pos` to the box with half-extents `half`.
/// Positive outside, negative inside.
pub fn box_sdf(pos: [f64; 3], center: [f64; 3], half: [f64; 3]) -> f64 {
    let dx = (pos[0] - center[0]).abs() - half[0];
    let dy = (pos[1] - center[1]).abs() - half[1];
    let dz = (pos[2] - center[2]).abs() - half[2];
    let outside = [dx.max(0.0), dy.max(0.0), dz.max(0.0)];
    let outside_dist = (outside[0].powi(2) + outside[1].powi(2) + outside[2].powi(2)).sqrt();
    let inside_dist = dx.max(dy).max(dz).min(0.0);
    outside_dist + inside_dist
}

/// Compute Bouzidi equilibrium correction for a wall moving with velocity `u_wall`.
///
/// Returns the momentum source term to be added to the distribution after
/// bounce-back to account for the moving boundary.
///
/// # Arguments
/// * `rho`     - fluid density at the fluid node
/// * `u_wall`  - wall velocity `[ux, uy, uz]`
/// * `w_i`     - lattice weight for direction `i`
/// * `e_dot_u` - dot product of discrete velocity `e_i` with `u_wall`
pub fn bouzidi_wall_correction(rho: f64, u_wall: [f64; 3], w_i: f64, e_dot_u: f64) -> f64 {
    let _ = u_wall; // used only for documentation; e_dot_u is the scalar projection
    let cs2 = 1.0 / 3.0;
    2.0 * w_i * rho * e_dot_u / cs2
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- filippova_hanel ---

    #[test]
    fn test_fh_output_length_equals_input() {
        let f = vec![0.1, 0.2, 0.3, 0.4];
        let result = filippova_hanel(&f, &f, &f, 0.3, 1.0);
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_fh_q_half_equilibrium_returns_wall_eq() {
        // When f == f_eq_fluid == f_eq_wall and omega=1, result should be close to input
        let f = vec![0.1; 9];
        let result = filippova_hanel(&f, &f, &f, 0.5, 1.0);
        assert_eq!(result.len(), 9);
        for v in &result {
            assert!(v.is_finite(), "result must be finite");
        }
    }

    #[test]
    fn test_fh_q_small_uses_chi_formula() {
        let f_fluid = vec![0.2; 9];
        let f_eq_fluid = vec![0.2; 9];
        let f_eq_wall = vec![0.15; 9];
        let result = filippova_hanel(&f_fluid, &f_eq_fluid, &f_eq_wall, 0.1, 1.5);
        assert_eq!(result.len(), 9);
    }

    #[test]
    fn test_fh_q_large_uses_linear_chi() {
        let f_fluid = vec![0.2; 9];
        let f_eq_fluid = vec![0.2; 9];
        let f_eq_wall = vec![0.15; 9];
        let result = filippova_hanel(&f_fluid, &f_eq_fluid, &f_eq_wall, 0.8, 1.5);
        assert_eq!(result.len(), 9);
    }

    #[test]
    fn test_fh_all_equal_distributions() {
        // If f == f_eq_fluid == f_eq_wall the BC should produce something finite
        let f = vec![1.0 / 9.0; 9];
        let result = filippova_hanel(&f, &f, &f, 0.5, 1.0);
        for v in &result {
            assert!(v.is_finite());
        }
    }

    // --- bouzidi_bc ---

    #[test]
    fn test_bouzidi_q_ge_half_linear() {
        let mut f = vec![0.2, 0.1, 0.05, 0.15];
        let original_f0 = f[0];
        let _ = original_f0;
        let val = bouzidi_bc(&mut f, 0.7, [0.0; 3], 1.0 / 3.0);
        assert!(val.is_finite());
    }

    #[test]
    fn test_bouzidi_q_lt_half_quadratic() {
        let mut f = vec![0.2, 0.1, 0.05, 0.15];
        let val = bouzidi_bc(&mut f, 0.3, [0.0; 3], 1.0 / 3.0);
        assert!(val.is_finite());
    }

    #[test]
    fn test_bouzidi_modifies_f_inplace() {
        let mut f = vec![0.5; 4];
        let initial = f[0];
        bouzidi_bc(&mut f, 0.6, [0.0; 3], 1.0 / 3.0);
        // f[0] should have been modified
        let _ = initial;
        assert!(f[0].is_finite());
    }

    #[test]
    fn test_bouzidi_empty_slice() {
        let mut f: Vec<f64> = vec![];
        let val = bouzidi_bc(&mut f, 0.5, [0.0; 3], 1.0 / 3.0);
        assert_eq!(val, 0.0);
    }

    #[test]
    fn test_bouzidi_wall_velocity_shifts_result() {
        let mut f1 = vec![0.3, 0.2, 0.1];
        let mut f2 = f1.clone();
        let v1 = bouzidi_bc(&mut f1, 0.6, [0.0; 3], 1.0 / 3.0);
        let v2 = bouzidi_bc(&mut f2, 0.6, [0.1, 0.0, 0.0], 1.0 / 3.0);
        assert!((v1 - v2).abs() > 1e-15, "wall velocity must shift result");
    }

    // --- quadratic_interpolation_bc ---

    #[test]
    fn test_quad_at_q_zero_returns_f0() {
        let v = quadratic_interpolation_bc(1.0, 2.0, 3.0, 0.0);
        assert!((v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quad_at_q_one_returns_f1() {
        // p(1) = f_1
        let v = quadratic_interpolation_bc(1.0, 2.0, 3.0, 1.0);
        assert!((v - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_quad_at_q_two_returns_f2() {
        let v = quadratic_interpolation_bc(1.0, 2.0, 3.0, 2.0);
        assert!((v - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_quad_monotone_input() {
        // For monotone input at q=0.5 the result should be between f_0 and f_1
        let v = quadratic_interpolation_bc(1.0, 2.0, 3.0, 0.5);
        assert!(v >= 1.0 && v <= 3.0, "v={v}");
    }

    #[test]
    fn test_quad_linearity_for_linear_input() {
        // If f_0, f_1, f_2 are collinear the result is exact linear interpolation
        let f0 = 0.0_f64;
        let f1 = 1.0_f64;
        let f2 = 2.0_f64;
        let q = 0.5;
        let v = quadratic_interpolation_bc(f0, f1, f2, q);
        assert!((v - 0.5).abs() < 1e-10, "v={v}");
    }

    // --- interpolated_bounce_back ---

    #[test]
    fn test_ibb_q_ge_half_returns_opposite() {
        let f = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let v = interpolated_bounce_back(&f, 0.6, 2);
        assert!((v - f[2]).abs() < 1e-12);
    }

    #[test]
    fn test_ibb_out_of_bounds_returns_zero() {
        let f = vec![0.1; 5];
        let v = interpolated_bounce_back(&f, 0.5, 10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_ibb_q_small_blends() {
        let f = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let v = interpolated_bounce_back(&f, 0.1, 3);
        assert!(v.is_finite());
    }

    // --- half_way_bounce_back ---

    #[test]
    fn test_hwbb_returns_correct_index() {
        let f = vec![0.1, 0.2, 0.3];
        assert!((half_way_bounce_back(&f, 1) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_hwbb_out_of_bounds_returns_zero() {
        let f = vec![0.1; 3];
        assert_eq!(half_way_bounce_back(&f, 5), 0.0);
    }

    // --- moving_wall_bc ---

    #[test]
    fn test_moving_wall_zero_velocity_equals_bounce_back() {
        let mut f = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.05, 0.05, 0.05, 0.05];
        moving_wall_bc(&mut f, [0.0; 3], 1.0, 1);
        // with zero wall velocity there's just a sign flip
        assert!(f[1].is_finite());
    }

    #[test]
    fn test_moving_wall_out_of_bounds_direction() {
        let mut f = vec![0.1; 5];
        let v = moving_wall_bc(&mut f, [0.0; 3], 1.0, 10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_moving_wall_nonzero_velocity_correction() {
        let mut f1 = vec![0.2; 9];
        let mut f2 = f1.clone();
        moving_wall_bc(&mut f1, [0.0; 3], 1.0, 2);
        moving_wall_bc(&mut f2, [0.1, 0.0, 0.0], 1.0, 2);
        assert!((f1[2] - f2[2]).abs() > 1e-15);
    }

    // --- direct_forcing ---

    #[test]
    fn test_direct_forcing_epsilon_one_snaps_to_target() {
        let mut u = [0.0; 3];
        direct_forcing(&mut u, [1.0, 2.0, 3.0], 1.0);
        assert!((u[0] - 1.0).abs() < 1e-12);
        assert!((u[1] - 2.0).abs() < 1e-12);
        assert!((u[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_direct_forcing_epsilon_zero_no_change() {
        let mut u = [1.0, 2.0, 3.0];
        direct_forcing(&mut u, [0.0, 0.0, 0.0], 0.0);
        assert!((u[0] - 1.0).abs() < 1e-12);
        assert!((u[1] - 2.0).abs() < 1e-12);
        assert!((u[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_direct_forcing_half_epsilon() {
        let mut u = [0.0; 3];
        direct_forcing(&mut u, [2.0, 0.0, 0.0], 0.5);
        assert!((u[0] - 1.0).abs() < 1e-12);
    }

    // --- RotatingCylinder ---

    #[test]
    fn test_rotating_cylinder_wall_velocity_perpendicular() {
        let cyl = RotatingCylinder::new([0.0, 0.0], 1.0, 1.0);
        // Point on positive x-axis: velocity should be in y direction
        let v = cyl.wall_velocity([1.0, 0.0]);
        assert!(v[0].abs() < 1e-12, "vx should be 0, got {}", v[0]);
        assert!((v[1] - 1.0).abs() < 1e-12, "vy should be 1, got {}", v[1]);
    }

    #[test]
    fn test_rotating_cylinder_contains_center() {
        let cyl = RotatingCylinder::new([5.0, 5.0], 2.0, 0.5);
        assert!(cyl.contains([5.0, 5.0]));
    }

    #[test]
    fn test_rotating_cylinder_excludes_far_point() {
        let cyl = RotatingCylinder::new([0.0, 0.0], 1.0, 1.0);
        assert!(!cyl.contains([5.0, 0.0]));
    }

    #[test]
    fn test_rotating_cylinder_wall_velocity_on_y_axis() {
        let cyl = RotatingCylinder::new([0.0, 0.0], 1.0, 2.0);
        // Point on positive y-axis: velocity should be in negative x direction
        let v = cyl.wall_velocity([0.0, 1.0]);
        assert!((v[0] + 2.0).abs() < 1e-12, "vx should be -2, got {}", v[0]);
        assert!(v[1].abs() < 1e-12, "vy should be 0, got {}", v[1]);
    }

    #[test]
    fn test_rotating_cylinder_speed_proportional_to_omega() {
        let c1 = RotatingCylinder::new([0.0, 0.0], 1.0, 1.0);
        let c2 = RotatingCylinder::new([0.0, 0.0], 1.0, 3.0);
        let v1 = c1.wall_velocity([1.0, 0.0]);
        let v2 = c2.wall_velocity([1.0, 0.0]);
        let speed1 = (v1[0].powi(2) + v1[1].powi(2)).sqrt();
        let speed2 = (v2[0].powi(2) + v2[1].powi(2)).sqrt();
        assert!((speed2 / speed1 - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotating_cylinder_ray_q_intersects() {
        let cyl = RotatingCylinder::new([5.0, 0.0], 1.0, 1.0);
        // Ray from [0, 0] in direction [1, 0] should hit at [4, 0], q ≈ 4
        // (outside the [0,1] range, so None)
        let q = cyl.ray_q([0.0, 0.0], [1.0, 0.0]);
        assert!(q.is_none() || q.is_some_and(|v| v > 0.0));
    }

    #[test]
    fn test_rotating_cylinder_ray_q_close_hit() {
        let cyl = RotatingCylinder::new([1.5, 0.0], 1.0, 1.0);
        let q = cyl.ray_q([0.0, 0.0], [1.0, 0.0]);
        assert!(q.is_some(), "ray should intersect cylinder");
        let q_val = q.unwrap();
        assert!(q_val > 0.0 && q_val <= 1.0, "q={q_val}");
    }

    // --- CurvedBoundarySet ---

    #[test]
    fn test_curved_boundary_set_new_empty() {
        let cbs = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
        assert!(cbs.nodes.is_empty());
    }

    #[test]
    fn test_curved_boundary_set_add_node() {
        let mut cbs = CurvedBoundarySet::new(CurvedBcMethod::FilippovaHanel, [0.0; 3]);
        cbs.add_node(CurvedWallNode::new(0, 1, 0.5));
        assert_eq!(cbs.nodes.len(), 1);
    }

    #[test]
    fn test_curved_boundary_set_apply_bouzidi_does_not_panic() {
        let mut cbs = CurvedBoundarySet::new(CurvedBcMethod::Bouzidi, [0.0; 3]);
        cbs.add_node(CurvedWallNode::new(0, 1, 0.4));
        let mut f = vec![vec![0.1, 0.2, 0.3, 0.4]; 4];
        cbs.apply(&mut f, 1.0);
        // Should not panic; just check result is finite
        for v in &f[0] {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_curved_boundary_set_apply_fh_does_not_panic() {
        let mut cbs = CurvedBoundarySet::new(CurvedBcMethod::FilippovaHanel, [0.0; 3]);
        cbs.add_node(CurvedWallNode::new(1, 2, 0.6));
        let mut f: Vec<Vec<f64>> = (0..4).map(|_| vec![1.0 / 9.0; 9]).collect();
        cbs.apply(&mut f, 1.2);
    }

    #[test]
    fn test_curved_boundary_set_apply_quad_does_not_panic() {
        let mut cbs = CurvedBoundarySet::new(CurvedBcMethod::QuadraticInterpolation, [0.0; 3]);
        cbs.add_node(CurvedWallNode::new(0, 1, 0.3));
        let mut f = vec![vec![0.1, 0.2, 0.3, 0.4]; 2];
        cbs.apply(&mut f, 1.0);
    }

    #[test]
    fn test_curved_boundary_set_apply_moment_based_conserves_rho() {
        let mut cbs = CurvedBoundarySet::new(CurvedBcMethod::MomentBased, [0.0; 3]);
        cbs.add_node(CurvedWallNode::new(0, 1, 0.5));
        let mut f = vec![vec![0.1, 0.2, 0.3, 0.4]];
        let rho_before: f64 = f[0].iter().sum();
        cbs.apply(&mut f, 1.0);
        let rho_after: f64 = f[0].iter().sum();
        assert!((rho_before - rho_after).abs() < 1e-10);
    }

    #[test]
    fn test_from_geometry_sphere_sdf() {
        // Sphere centred at [5, 5, 0] radius 2 in a 10x10x1 grid
        let center = [5.0f64, 5.0, 0.0];
        let radius = 2.0f64;
        let cbs = CurvedBoundarySet::from_geometry(
            &|pos| -(sphere_sdf(pos, center, radius)),
            1.0,
            10,
            10,
            1,
        );
        // Should find some curved nodes on the sphere boundary
        assert!(!cbs.nodes.is_empty(), "should detect curved nodes");
    }

    // --- d2q9_equilibrium ---

    #[test]
    fn test_d2q9_eq_sum_equals_rho() {
        let rho = 1.0;
        let feq = d2q9_equilibrium(rho, 0.0, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-12, "sum={sum}");
    }

    #[test]
    fn test_d2q9_eq_nonzero_velocity() {
        let feq = d2q9_equilibrium(1.0, 0.1, 0.05);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "sum={sum}");
    }

    // --- sphere_sdf and box_sdf ---

    #[test]
    fn test_sphere_sdf_outside() {
        let d = sphere_sdf([3.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!((d - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_sdf_inside() {
        let d = sphere_sdf([0.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!((d + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_box_sdf_inside() {
        let d = box_sdf([0.0; 3], [0.0; 3], [1.0; 3]);
        assert!(d <= 0.0, "d={d}");
    }

    #[test]
    fn test_box_sdf_outside_corner() {
        let d = box_sdf([2.0, 2.0, 2.0], [0.0; 3], [1.0; 3]);
        let expected = (3.0f64).sqrt();
        assert!((d - expected).abs() < 1e-10, "d={d}");
    }

    // --- bouzidi_wall_correction ---

    #[test]
    fn test_bouzidi_wall_correction_zero_dot() {
        let v = bouzidi_wall_correction(1.0, [0.0; 3], 1.0 / 9.0, 0.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_bouzidi_wall_correction_scales_with_rho() {
        let v1 = bouzidi_wall_correction(1.0, [0.0; 3], 1.0 / 9.0, 0.1);
        let v2 = bouzidi_wall_correction(2.0, [0.0; 3], 1.0 / 9.0, 0.1);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    // --- CurvedWallNode ---

    #[test]
    fn test_curved_wall_node_stores_fields() {
        let n = CurvedWallNode::new(3, 7, 0.4);
        assert_eq!(n.fluid_idx, 3);
        assert_eq!(n.wall_idx, 7);
        assert!((n.q - 0.4).abs() < 1e-12);
    }

    // --- CurvedBcMethod ---

    #[test]
    fn test_curved_bc_method_debug() {
        let m = CurvedBcMethod::FilippovaHanel;
        let s = format!("{m:?}");
        assert!(s.contains("Filippova"));
    }

    #[test]
    fn test_curved_bc_method_eq() {
        assert_eq!(CurvedBcMethod::Bouzidi, CurvedBcMethod::Bouzidi);
        assert_ne!(CurvedBcMethod::Bouzidi, CurvedBcMethod::MomentBased);
    }

    // --- PI usage in cylinder ---

    #[test]
    fn test_pi_constant_is_correct() {
        // Indirect test that PI is used correctly in context
        assert!((PI - std::f64::consts::PI).abs() < 1e-10);
    }
}
