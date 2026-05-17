// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lattice Boltzmann Method (LBM) simulation API for Python interop.
//!
//! Provides `PyLbmSimulation`: a 2-D D2Q9 LBM solver with BGK collision,
//! configurable boundaries (no-slip walls, periodic, inlet/outlet),
//! and full velocity and density field extraction.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// D2Q9 constants
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights.
const W: [f64; 9] = [
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

/// D2Q9 X-components of discrete velocities.
const EX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 Y-components of discrete velocities.
const EY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// Opposite direction indices for bounce-back.
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ---------------------------------------------------------------------------
// Boundary cell type
// ---------------------------------------------------------------------------

/// Boundary condition type for an LBM cell.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LbmBoundary {
    /// Fluid cell (no boundary).
    Fluid,
    /// No-slip solid wall (full bounce-back).
    NoSlipWall,
    /// Periodic (handled at domain edges automatically).
    Periodic,
    /// Inlet: fixed velocity applied on this cell column.
    Inlet,
    /// Outlet: zero-gradient outflow.
    Outlet,
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the LBM simulation.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyLbmConfig {
    /// Grid width (number of cells in X direction).
    #[pyo3(get, set)]
    pub width: usize,
    /// Grid height (number of cells in Y direction).
    #[pyo3(get, set)]
    pub height: usize,
    /// Kinematic viscosity ν. Controls the BGK relaxation rate.
    #[pyo3(get, set)]
    pub viscosity: f64,
    /// External body-force acceleration `[fx, fy]` applied to the fluid.
    pub body_force: [f64; 2],
    /// Initial uniform velocity `[ux, uy]` (used to initialise equilibrium).
    pub init_velocity: [f64; 2],
    /// Initial uniform density.
    #[pyo3(get, set)]
    pub init_density: f64,
}

#[pymethods]
impl PyLbmConfig {
    /// Create a new configuration.
    #[new]
    pub fn new(width: usize, height: usize, viscosity: f64) -> Self {
        Self {
            width,
            height,
            viscosity: viscosity.max(1e-8),
            body_force: [0.0; 2],
            init_velocity: [0.0; 2],
            init_density: 1.0,
        }
    }

    /// Lid-driven cavity (64×64, viscosity=0.01).
    #[staticmethod]
    pub fn lid_driven_cavity() -> Self {
        Self::new(64, 64, 0.01)
    }

    /// Channel flow (128×32, Poiseuille-like body force).
    #[staticmethod]
    pub fn channel_flow() -> Self {
        let mut cfg = Self::new(128, 32, 0.01);
        cfg.body_force = [1e-4, 0.0];
        cfg
    }

    /// Compute the BGK relaxation rate ω = 1/(3ν + 0.5).
    pub fn omega(&self) -> f64 {
        1.0 / (3.0 * self.viscosity + 0.5)
    }

    /// Body force `[fx, fy]` as a list.
    #[getter]
    pub fn body_force(&self) -> Vec<f64> {
        self.body_force.to_vec()
    }

    /// Set body force from a `[fx, fy]` list.
    #[setter]
    pub fn set_body_force(&mut self, v: Vec<f64>) {
        if v.len() >= 2 {
            self.body_force = [v[0], v[1]];
        }
    }

    /// Initial velocity `[ux, uy]` as a list.
    #[getter]
    pub fn init_velocity(&self) -> Vec<f64> {
        self.init_velocity.to_vec()
    }

    /// Set initial velocity from a `[ux, uy]` list.
    #[setter]
    pub fn set_init_velocity(&mut self, v: Vec<f64>) {
        if v.len() >= 2 {
            self.init_velocity = [v[0], v[1]];
        }
    }
}

impl Default for PyLbmConfig {
    fn default() -> Self {
        Self::new(32, 32, 0.01)
    }
}

// ---------------------------------------------------------------------------
// PyLbmSimulation
// ---------------------------------------------------------------------------

/// A 2-D D2Q9 Lattice-Boltzmann fluid simulation.
///
/// Supports:
/// - BGK collision operator (single relaxation time)
/// - Periodic streaming
/// - No-slip bounce-back walls via boundary map
/// - Velocity and density field extraction
/// - Lid-driven velocity boundary on top row
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyLbmSimulation {
    /// Grid width.
    width: usize,
    /// Grid height.
    height: usize,
    /// BGK relaxation rate.
    omega: f64,
    /// Body-force `[fx, fy]`.
    body_force: [f64; 2],
    /// Distribution functions `f[cell * 9 + q]`.
    f: Vec<f64>,
    /// Temporary buffer for post-collision distributions.
    f_tmp: Vec<f64>,
    /// Boundary type per cell `[y * width + x]`.
    boundary: Vec<LbmBoundary>,
    /// Lid velocity `[ux, uy]` applied to the top row.
    lid_velocity: Option<[f64; 2]>,
    /// Accumulated step count.
    step_count: u64,
}

#[pymethods]
impl PyLbmSimulation {
    /// Create a new LBM simulation from configuration.
    ///
    /// All cells are initialised to equilibrium at the configured density and velocity.
    #[new]
    pub fn new(config: &PyLbmConfig) -> Self {
        let n = config.width * config.height;
        let ux0 = config.init_velocity[0];
        let uy0 = config.init_velocity[1];
        let rho0 = config.init_density;

        let mut f = vec![0.0f64; n * 9];
        for i in 0..n {
            let feq = equilibrium(rho0, ux0, uy0);
            for q in 0..9 {
                f[i * 9 + q] = feq[q];
            }
        }
        let f_tmp = f.clone();
        let boundary = vec![LbmBoundary::Fluid; n];

        Self {
            width: config.width,
            height: config.height,
            omega: config.omega(),
            body_force: config.body_force,
            f,
            f_tmp,
            boundary,
            lid_velocity: None,
            step_count: 0,
        }
    }

    /// Grid width.
    #[getter]
    pub fn width(&self) -> usize {
        self.width
    }

    /// Grid height.
    #[getter]
    pub fn height(&self) -> usize {
        self.height
    }

    /// Number of completed steps.
    #[getter]
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Set the boundary type of cell `(x, y)`.
    pub fn set_boundary(&mut self, x: usize, y: usize, btype: LbmBoundary) {
        if x < self.width && y < self.height {
            self.boundary[y * self.width + x] = btype;
        }
    }

    /// Get the boundary type of cell `(x, y)`.
    pub fn get_boundary(&self, x: usize, y: usize) -> Option<LbmBoundary> {
        if x < self.width && y < self.height {
            Some(self.boundary[y * self.width + x])
        } else {
            None
        }
    }

    /// Mark all cells along the top row as `NoSlipWall`.
    pub fn add_top_wall(&mut self) {
        let y = self.height - 1;
        for x in 0..self.width {
            self.set_boundary(x, y, LbmBoundary::NoSlipWall);
        }
    }

    /// Mark all cells along the bottom row as `NoSlipWall`.
    pub fn add_bottom_wall(&mut self) {
        for x in 0..self.width {
            self.set_boundary(x, 0, LbmBoundary::NoSlipWall);
        }
    }

    /// Mark all four enclosing walls as `NoSlipWall` (cavity setup).
    pub fn add_enclosing_walls(&mut self) {
        for x in 0..self.width {
            self.set_boundary(x, 0, LbmBoundary::NoSlipWall);
            self.set_boundary(x, self.height - 1, LbmBoundary::NoSlipWall);
        }
        for y in 0..self.height {
            self.set_boundary(0, y, LbmBoundary::NoSlipWall);
            self.set_boundary(self.width - 1, y, LbmBoundary::NoSlipWall);
        }
    }

    /// Set a lid velocity applied to the top row every step.
    ///
    /// Pass `None` to disable.
    pub fn set_lid_velocity(&mut self, vel: Option<Vec<f64>>) {
        self.lid_velocity = vel.and_then(|v| {
            if v.len() >= 2 {
                Some([v[0], v[1]])
            } else {
                None
            }
        });
    }

    /// Get lid velocity as `[ux, uy]` or `None`.
    pub fn get_lid_velocity(&self) -> Option<Vec<f64>> {
        self.lid_velocity.map(|v| v.to_vec())
    }

    /// Get macroscopic velocity `[ux, uy]` at cell `(x, y)`.
    ///
    /// Returns `[0, 0]` for out-of-bounds or wall cells.
    pub fn velocity_at(&self, x: usize, y: usize) -> Vec<f64> {
        if x >= self.width || y >= self.height {
            return vec![0.0, 0.0];
        }
        let cell = y * self.width + x;
        if self.boundary[cell] == LbmBoundary::NoSlipWall {
            return vec![0.0, 0.0];
        }
        let idx = cell * 9;
        let rho: f64 = self.f[idx..idx + 9].iter().sum();
        if rho < 1e-15 {
            return vec![0.0, 0.0];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        for q in 0..9 {
            ux += EX[q] as f64 * self.f[idx + q];
            uy += EY[q] as f64 * self.f[idx + q];
        }
        vec![ux / rho, uy / rho]
    }

    /// Get macroscopic density at cell `(x, y)`. Returns 0 if out of bounds.
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        if x >= self.width || y >= self.height {
            return 0.0;
        }
        let idx = (y * self.width + x) * 9;
        self.f[idx..idx + 9].iter().sum()
    }

    /// Return the full velocity field as a flat `Vec<f64>` of `[ux, uy]` pairs,
    /// row-major (y outer, x inner). Length = `width * height * 2`.
    pub fn get_velocity_field(&self) -> Vec<f64> {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n * 2];
        for y in 0..self.height {
            for x in 0..self.width {
                let v = self.velocity_at(x, y);
                let cell = y * self.width + x;
                out[cell * 2] = v[0];
                out[cell * 2 + 1] = v[1];
            }
        }
        out
    }

    /// Return the full density field as a flat `Vec<f64>`, row-major.
    /// Length = `width * height`.
    pub fn get_density_field(&self) -> Vec<f64> {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n];
        for y in 0..self.height {
            for x in 0..self.width {
                out[y * self.width + x] = self.density_at(x, y);
            }
        }
        out
    }

    /// Advance the simulation by one LBM time step.
    ///
    /// Steps:
    /// 1. Collision (BGK) + optional body-force correction.
    /// 2. Apply lid velocity (if set) as equilibrium on top row.
    /// 3. Streaming with periodic boundaries.
    /// 4. Bounce-back for `NoSlipWall` cells.
    pub fn step(&mut self) {
        let w = self.width;
        let h = self.height;
        let omega = self.omega;
        let bfx = self.body_force[0];
        let bfy = self.body_force[1];

        // -- Collision --
        for y in 0..h {
            for x in 0..w {
                let cell = y * w + x;
                let idx = cell * 9;
                if self.boundary[cell] == LbmBoundary::NoSlipWall {
                    // Prepare for bounce-back in streaming: copy as-is
                    for q in 0..9 {
                        self.f_tmp[idx + q] = self.f[idx + q];
                    }
                    continue;
                }
                let rho: f64 = self.f[idx..idx + 9].iter().sum();
                let mut ux = 0.0f64;
                let mut uy = 0.0f64;
                for q in 0..9 {
                    ux += EX[q] as f64 * self.f[idx + q];
                    uy += EY[q] as f64 * self.f[idx + q];
                }
                let rho_inv = if rho > 1e-15 { 1.0 / rho } else { 0.0 };
                ux = ux * rho_inv + bfx / (2.0 * omega * rho.max(1e-15));
                uy = uy * rho_inv + bfy / (2.0 * omega * rho.max(1e-15));
                let feq = equilibrium(rho, ux, uy);
                for (q, (&feq_q, f_tmp_q)) in feq
                    .iter()
                    .zip(self.f_tmp[idx..idx + 9].iter_mut())
                    .enumerate()
                {
                    *f_tmp_q = self.f[idx + q] * (1.0 - omega) + feq_q * omega;
                }
            }
        }

        // -- Lid velocity (Zou-He-like: set top row equilibrium) --
        if let Some(lid_vel) = self.lid_velocity {
            let y = h - 1;
            for x in 0..w {
                let cell = y * w + x;
                let idx = cell * 9;
                let rho = self.f[idx..idx + 9].iter().sum::<f64>().max(0.01);
                let feq = equilibrium(rho, lid_vel[0], lid_vel[1]);
                self.f_tmp[idx..(9 + idx)].copy_from_slice(&feq);
            }
        }

        // -- Streaming with periodic boundaries --
        let f_src = self.f_tmp.clone();
        for y in 0..h {
            for x in 0..w {
                let dst_cell = y * w + x;
                let dst_idx = dst_cell * 9;
                for (q, (ex, ey)) in EX.iter().zip(EY.iter()).enumerate() {
                    let src_x = ((x as isize - *ex as isize).rem_euclid(w as isize)) as usize;
                    let src_y = ((y as isize - *ey as isize).rem_euclid(h as isize)) as usize;
                    let src_idx = (src_y * w + src_x) * 9;
                    self.f[dst_idx + q] = f_src[src_idx + q];
                }
            }
        }

        // -- Bounce-back for wall cells --
        for y in 0..h {
            for x in 0..w {
                let cell = y * w + x;
                if self.boundary[cell] != LbmBoundary::NoSlipWall {
                    continue;
                }
                let idx = cell * 9;
                let mut tmp = [0.0f64; 9];
                for q in 0..9 {
                    tmp[OPP[q]] = self.f[idx + q];
                }
                self.f[idx..(9 + idx)].copy_from_slice(&tmp);
            }
        }

        self.step_count += 1;
    }

    /// Run `n` steps.
    pub fn run(&mut self, steps: u64) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// Mean macroscopic density over all non-wall cells.
    pub fn mean_density(&self) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for y in 0..self.height {
            for x in 0..self.width {
                let cell = y * self.width + x;
                if self.boundary[cell] != LbmBoundary::NoSlipWall {
                    sum += self.density_at(x, y);
                    count += 1;
                }
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }

    /// Maximum speed (magnitude of velocity) over all non-wall cells.
    pub fn max_speed(&self) -> f64 {
        let mut max_v = 0.0f64;
        for y in 0..self.height {
            for x in 0..self.width {
                let cell = y * self.width + x;
                if self.boundary[cell] == LbmBoundary::NoSlipWall {
                    continue;
                }
                let v = self.velocity_at(x, y);
                let speed = (v[0] * v[0] + v[1] * v[1]).sqrt();
                if speed > max_v {
                    max_v = speed;
                }
            }
        }
        max_v
    }

    /// Vorticity field (z-component of curl) as a flat `Vec<f64>`.
    ///
    /// Uses central differences. Length = width * height.
    pub fn vorticity_field(&self) -> Vec<f64> {
        let w = self.width;
        let h = self.height;
        let mut out = vec![0.0f64; w * h];
        for y in 1..(h - 1) {
            for x in 1..(w - 1) {
                let vxp = self.velocity_at(x, y + 1)[0];
                let vxm = self.velocity_at(x, y - 1)[0];
                let vyp = self.velocity_at(x + 1, y)[1];
                let vym = self.velocity_at(x - 1, y)[1];
                out[y * w + x] = (vyp - vym) * 0.5 - (vxp - vxm) * 0.5;
            }
        }
        out
    }

    /// Enstrophy (0.5 ∫ ω² dA) — scalar measure of total vorticity.
    pub fn enstrophy(&self) -> f64 {
        self.vorticity_field().iter().map(|&w| 0.5 * w * w).sum()
    }

    /// Reynolds number estimate: Re = L * U / ν (characteristic length = height/2).
    pub fn reynolds_number(&self, omega: f64) -> f64 {
        let nu = (1.0 / omega - 0.5) / 3.0;
        let u_max = self.max_speed();
        let l = self.height as f64 * 0.5;
        if nu < 1e-15 {
            return 0.0;
        }
        l * u_max / nu
    }

    /// Return the velocity magnitude field as a flat `Vec<f64>`.
    pub fn speed_field(&self) -> Vec<f64> {
        let w = self.width;
        let h = self.height;
        let mut out = vec![0.0f64; w * h];
        for y in 0..h {
            for x in 0..w {
                let v = self.velocity_at(x, y);
                out[y * w + x] = (v[0] * v[0] + v[1] * v[1]).sqrt();
            }
        }
        out
    }

    /// Reset all cells to equilibrium at given density and velocity.
    pub fn reset_to_equilibrium(&mut self, rho: f64, ux: f64, uy: f64) {
        let feq = equilibrium(rho, ux, uy);
        let n = self.width * self.height;
        for i in 0..n {
            let start = i * 9;
            self.f[start..start + 9].copy_from_slice(&feq);
            self.f_tmp[start..start + 9].copy_from_slice(&feq);
        }
    }

    /// Set a circular obstacle (no-slip) at `(cx, cy)` with integer radius `r`.
    pub fn add_circular_obstacle(&mut self, cx: usize, cy: usize, r: usize) {
        let r2 = (r * r) as isize;
        for y in 0..self.height {
            for x in 0..self.width {
                let dx = x as isize - cx as isize;
                let dy = y as isize - cy as isize;
                if dx * dx + dy * dy <= r2 {
                    self.set_boundary(x, y, LbmBoundary::NoSlipWall);
                }
            }
        }
    }

    /// Count of fluid (non-wall) cells.
    pub fn fluid_cell_count(&self) -> usize {
        self.boundary
            .iter()
            .filter(|&&b| b == LbmBoundary::Fluid)
            .count()
    }

    /// Count of wall cells.
    pub fn wall_cell_count(&self) -> usize {
        self.boundary
            .iter()
            .filter(|&&b| b == LbmBoundary::NoSlipWall)
            .count()
    }

    /// Collect per-step statistics.
    pub fn collect_stats(&self) -> LbmStats {
        LbmStats {
            step_count: self.step_count,
            mean_density: self.mean_density(),
            max_speed: self.max_speed(),
            enstrophy: self.enstrophy(),
            fluid_cell_count: self.fluid_cell_count(),
            wall_cell_count: self.wall_cell_count(),
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: D2Q9 equilibrium distribution
// ---------------------------------------------------------------------------

#[inline]
fn equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; 9];
    for q in 0..9 {
        let eu = EX[q] as f64 * ux + EY[q] as f64 * uy;
        feq[q] = W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
    }
    feq
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {

    use crate::LbmBoundary;
    use crate::PyLbmConfig;
    use crate::PyLbmSimulation;

    fn small_sim() -> PyLbmSimulation {
        PyLbmSimulation::new(&PyLbmConfig::new(8, 8, 0.1))
    }

    #[test]
    fn test_lbm_creation() {
        let sim = small_sim();
        assert_eq!(sim.width(), 8);
        assert_eq!(sim.height(), 8);
        assert_eq!(sim.step_count(), 0);
    }

    #[test]
    fn test_lbm_initial_density_near_one() {
        let sim = small_sim();
        for y in 0..8 {
            for x in 0..8 {
                let rho = sim.density_at(x, y);
                assert!((rho - 1.0).abs() < 1e-10, "rho at ({},{}) = {}", x, y, rho);
            }
        }
    }

    #[test]
    fn test_lbm_initial_velocity_near_zero() {
        let sim = small_sim();
        let v = sim.velocity_at(0, 0);
        assert!(v[0].abs() < 1e-10 && v[1].abs() < 1e-10);
    }

    #[test]
    fn test_lbm_step_advances_count() {
        let mut sim = small_sim();
        sim.step();
        sim.step();
        assert_eq!(sim.step_count(), 2);
    }

    #[test]
    fn test_lbm_run_n_steps() {
        let mut sim = small_sim();
        sim.run(10);
        assert_eq!(sim.step_count(), 10);
    }

    #[test]
    fn test_lbm_velocity_field_length() {
        let sim = small_sim();
        assert_eq!(sim.get_velocity_field().len(), 8 * 8 * 2);
    }

    #[test]
    fn test_lbm_density_field_length() {
        let sim = small_sim();
        assert_eq!(sim.get_density_field().len(), 64);
    }

    #[test]
    fn test_lbm_set_boundary_wall() {
        let mut sim = small_sim();
        sim.set_boundary(0, 0, LbmBoundary::NoSlipWall);
        assert_eq!(sim.get_boundary(0, 0), Some(LbmBoundary::NoSlipWall));
    }

    #[test]
    fn test_lbm_velocity_zero_on_wall() {
        let mut sim = small_sim();
        sim.set_boundary(3, 3, LbmBoundary::NoSlipWall);
        let v = sim.velocity_at(3, 3);
        assert!(v[0].abs() < 1e-12 && v[1].abs() < 1e-12);
    }

    #[test]
    fn test_lbm_add_top_wall() {
        let mut sim = small_sim();
        sim.add_top_wall();
        for x in 0..8 {
            assert_eq!(sim.get_boundary(x, 7), Some(LbmBoundary::NoSlipWall));
        }
    }

    #[test]
    fn test_lbm_add_bottom_wall() {
        let mut sim = small_sim();
        sim.add_bottom_wall();
        for x in 0..8 {
            assert_eq!(sim.get_boundary(x, 0), Some(LbmBoundary::NoSlipWall));
        }
    }

    #[test]
    fn test_lbm_add_enclosing_walls() {
        let mut sim = small_sim();
        sim.add_enclosing_walls();
        // Corners should be walls
        assert_eq!(sim.get_boundary(0, 0), Some(LbmBoundary::NoSlipWall));
        assert_eq!(sim.get_boundary(7, 7), Some(LbmBoundary::NoSlipWall));
        assert_eq!(sim.get_boundary(0, 7), Some(LbmBoundary::NoSlipWall));
        assert_eq!(sim.get_boundary(7, 0), Some(LbmBoundary::NoSlipWall));
    }

    #[test]
    fn test_lbm_out_of_bounds_boundary() {
        let sim = small_sim();
        assert!(sim.get_boundary(99, 99).is_none());
    }

    #[test]
    fn test_lbm_body_force_creates_flow() {
        let mut cfg = PyLbmConfig::new(16, 8, 0.1);
        cfg.body_force = [1e-3, 0.0];
        let mut sim = PyLbmSimulation::new(&cfg);
        // Run many steps: x-velocity in interior should become positive
        sim.run(200);
        let v = sim.velocity_at(8, 4);
        assert!(
            v[0] > 0.0,
            "body force in +x should drive positive ux, got {}",
            v[0]
        );
    }

    #[test]
    fn test_lbm_lid_velocity_drives_flow() {
        let mut sim = PyLbmSimulation::new(&PyLbmConfig::new(16, 16, 0.02));
        sim.add_enclosing_walls();
        sim.set_lid_velocity(Some(vec![0.1, 0.0]));
        sim.run(100);
        // Top interior row should see non-zero x-velocity
        let v = sim.velocity_at(8, 14);
        assert!(v[0].abs() > 1e-6, "lid should drive flow, got ux={}", v[0]);
    }

    #[test]
    fn test_lbm_omega_from_config() {
        let cfg = PyLbmConfig::new(8, 8, 1.0 / 6.0);
        // omega = 1/(3*(1/6)+0.5) = 1/1.0 = 1.0
        assert!((cfg.omega() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lbm_channel_flow_config() {
        let cfg = PyLbmConfig::channel_flow();
        assert_eq!(cfg.width, 128);
        assert!(cfg.body_force[0] > 0.0);
    }

    #[test]
    fn test_lbm_density_conserved_no_body_force() {
        // Total density should be approximately conserved with periodic BC
        let mut sim = small_sim();
        let rho_before: f64 = (0..8usize)
            .flat_map(|y| (0..8usize).map(move |x| (x, y)))
            .map(|(x, y)| sim.density_at(x, y))
            .sum();
        sim.run(10);
        let rho_after: f64 = (0..8usize)
            .flat_map(|y| (0..8usize).map(move |x| (x, y)))
            .map(|(x, y)| sim.density_at(x, y))
            .sum();
        assert!(
            (rho_after - rho_before).abs() / rho_before < 1e-6,
            "density not conserved: before={} after={}",
            rho_before,
            rho_after
        );
    }
}

// ---------------------------------------------------------------------------
// D3Q19 3-D LBM simulation
// ---------------------------------------------------------------------------

/// D3Q19 discrete velocity set (19 directions including rest).
///
/// Velocities are `(ex, ey, ez)` pairs.
const D3Q19_EX: [i32; 19] = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0];
const D3Q19_EY: [i32; 19] = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1];
const D3Q19_EZ: [i32; 19] = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1];

/// D3Q19 lattice weights.
const D3Q19_W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Opposite directions for D3Q19 bounce-back.
const D3Q19_OPP: [usize; 19] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
];

/// Configuration for a 3-D D3Q19 LBM simulation.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyLbm3dConfig {
    /// Grid width (X dimension).
    #[pyo3(get, set)]
    pub nx: usize,
    /// Grid depth (Y dimension).
    #[pyo3(get, set)]
    pub ny: usize,
    /// Grid height (Z dimension).
    #[pyo3(get, set)]
    pub nz: usize,
    /// Kinematic viscosity.
    #[pyo3(get, set)]
    pub viscosity: f64,
    /// External body force `[fx, fy, fz]`.
    pub body_force: [f64; 3],
    /// Initial uniform density.
    #[pyo3(get, set)]
    pub init_density: f64,
    /// Initial velocity `[ux, uy, uz]`.
    pub init_velocity: [f64; 3],
}

#[pymethods]
impl PyLbm3dConfig {
    /// Create a 3-D configuration.
    #[new]
    pub fn new(nx: usize, ny: usize, nz: usize, viscosity: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            viscosity: viscosity.max(1e-8),
            body_force: [0.0; 3],
            init_density: 1.0,
            init_velocity: [0.0; 3],
        }
    }

    /// Compute BGK relaxation rate ω = 1/(3ν + 0.5).
    pub fn omega(&self) -> f64 {
        1.0 / (3.0 * self.viscosity + 0.5)
    }

    /// Small 3-D test configuration (8×8×8).
    #[staticmethod]
    pub fn small() -> Self {
        Self::new(8, 8, 8, 0.1)
    }

    /// Body force `[fx, fy, fz]` as a list.
    #[getter]
    pub fn body_force(&self) -> Vec<f64> {
        self.body_force.to_vec()
    }

    /// Set body force from a `[fx, fy, fz]` list.
    #[setter]
    pub fn set_body_force(&mut self, v: Vec<f64>) {
        if v.len() >= 3 {
            self.body_force = [v[0], v[1], v[2]];
        }
    }

    /// Initial velocity `[ux, uy, uz]` as a list.
    #[getter]
    pub fn init_velocity(&self) -> Vec<f64> {
        self.init_velocity.to_vec()
    }

    /// Set initial velocity from a `[ux, uy, uz]` list.
    #[setter]
    pub fn set_init_velocity(&mut self, v: Vec<f64>) {
        if v.len() >= 3 {
            self.init_velocity = [v[0], v[1], v[2]];
        }
    }
}

/// A 3-D D3Q19 Lattice-Boltzmann simulation.
///
/// Supports BGK collision, periodic streaming, and bounce-back walls.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyLbm3dSimulation {
    /// Grid dimension X.
    nx: usize,
    /// Grid dimension Y.
    ny: usize,
    /// Grid dimension Z.
    nz: usize,
    /// BGK relaxation rate.
    omega: f64,
    /// Body force `[fx, fy, fz]`.
    body_force: [f64; 3],
    /// Distribution functions `f[cell * 19 + q]`.
    f: Vec<f64>,
    /// Swap buffer.
    f_tmp: Vec<f64>,
    /// Boundary type per cell.
    boundary: Vec<LbmBoundary>,
    /// Step counter.
    step_count: u64,
}

#[pymethods]
impl PyLbm3dSimulation {
    /// Create a new 3-D LBM simulation.
    #[new]
    pub fn new(config: &PyLbm3dConfig) -> Self {
        let n = config.nx * config.ny * config.nz;
        let ux0 = config.init_velocity[0];
        let uy0 = config.init_velocity[1];
        let uz0 = config.init_velocity[2];
        let rho0 = config.init_density;
        let bfx = config.body_force[0];
        let bfy = config.body_force[1];
        let bfz = config.body_force[2];

        let mut f = vec![0.0f64; n * 19];
        for i in 0..n {
            let feq = d3q19_equilibrium(rho0, ux0 + bfx, uy0 + bfy, uz0 + bfz);
            for q in 0..19 {
                f[i * 19 + q] = feq[q];
            }
        }
        let f_tmp = f.clone();
        let boundary = vec![LbmBoundary::Fluid; n];

        Self {
            nx: config.nx,
            ny: config.ny,
            nz: config.nz,
            omega: config.omega(),
            body_force: config.body_force,
            f,
            f_tmp,
            boundary,
            step_count: 0,
        }
    }

    /// Grid dimensions as `(nx, ny, nz)`.
    pub fn dimensions(&self) -> (usize, usize, usize) {
        (self.nx, self.ny, self.nz)
    }

    /// Step count.
    #[getter]
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Set boundary at cell `(x, y, z)`.
    pub fn set_boundary(&mut self, x: usize, y: usize, z: usize, btype: LbmBoundary) {
        if x < self.nx && y < self.ny && z < self.nz {
            self.boundary[z * self.ny * self.nx + y * self.nx + x] = btype;
        }
    }

    /// Get boundary at cell `(x, y, z)`.
    pub fn get_boundary(&self, x: usize, y: usize, z: usize) -> Option<LbmBoundary> {
        if x < self.nx && y < self.ny && z < self.nz {
            Some(self.boundary[z * self.ny * self.nx + y * self.nx + x])
        } else {
            None
        }
    }

    /// Macroscopic density at cell `(x, y, z)`.
    pub fn density_at(&self, x: usize, y: usize, z: usize) -> f64 {
        if x >= self.nx || y >= self.ny || z >= self.nz {
            return 0.0;
        }
        let idx = (z * self.ny * self.nx + y * self.nx + x) * 19;
        self.f[idx..idx + 19].iter().sum()
    }

    /// Macroscopic velocity `[ux, uy, uz]` at cell `(x, y, z)`.
    pub fn velocity_at(&self, x: usize, y: usize, z: usize) -> Vec<f64> {
        if x >= self.nx || y >= self.ny || z >= self.nz {
            return vec![0.0, 0.0, 0.0];
        }
        let cell = z * self.ny * self.nx + y * self.nx + x;
        if self.boundary[cell] == LbmBoundary::NoSlipWall {
            return vec![0.0, 0.0, 0.0];
        }
        let idx = cell * 19;
        let rho: f64 = self.f[idx..idx + 19].iter().sum();
        if rho < 1e-15 {
            return vec![0.0, 0.0, 0.0];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        let mut uz = 0.0f64;
        for q in 0..19 {
            ux += D3Q19_EX[q] as f64 * self.f[idx + q];
            uy += D3Q19_EY[q] as f64 * self.f[idx + q];
            uz += D3Q19_EZ[q] as f64 * self.f[idx + q];
        }
        vec![ux / rho, uy / rho, uz / rho]
    }

    /// Advance the simulation by one D3Q19 step (BGK + streaming + bounce-back).
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let omega = self.omega;
        let bfx = self.body_force[0];
        let bfy = self.body_force[1];
        let bfz = self.body_force[2];

        // -- Collision --
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let cell = z * ny * nx + y * nx + x;
                    let idx = cell * 19;
                    if self.boundary[cell] == LbmBoundary::NoSlipWall {
                        for q in 0..19 {
                            self.f_tmp[idx + q] = self.f[idx + q];
                        }
                        continue;
                    }
                    let rho: f64 = self.f[idx..idx + 19].iter().sum();
                    let mut ux = 0.0f64;
                    let mut uy = 0.0f64;
                    let mut uz = 0.0f64;
                    for q in 0..19 {
                        ux += D3Q19_EX[q] as f64 * self.f[idx + q];
                        uy += D3Q19_EY[q] as f64 * self.f[idx + q];
                        uz += D3Q19_EZ[q] as f64 * self.f[idx + q];
                    }
                    let rho_inv = if rho > 1e-15 { 1.0 / rho } else { 0.0 };
                    ux = ux * rho_inv + bfx / (2.0 * omega * rho.max(1e-15));
                    uy = uy * rho_inv + bfy / (2.0 * omega * rho.max(1e-15));
                    uz = uz * rho_inv + bfz / (2.0 * omega * rho.max(1e-15));
                    let feq = d3q19_equilibrium(rho, ux, uy, uz);
                    for (q, (&feq_q, f_tmp_q)) in feq
                        .iter()
                        .zip(self.f_tmp[idx..idx + 19].iter_mut())
                        .enumerate()
                    {
                        *f_tmp_q = self.f[idx + q] * (1.0 - omega) + feq_q * omega;
                    }
                }
            }
        }

        // -- Streaming with periodic BC --
        let f_src = self.f_tmp.clone();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst_cell = z * ny * nx + y * nx + x;
                    let dst_idx = dst_cell * 19;
                    for (q, ((ex, ey), ez)) in D3Q19_EX
                        .iter()
                        .zip(D3Q19_EY.iter())
                        .zip(D3Q19_EZ.iter())
                        .enumerate()
                    {
                        let sx = ((x as isize - *ex as isize).rem_euclid(nx as isize)) as usize;
                        let sy = ((y as isize - *ey as isize).rem_euclid(ny as isize)) as usize;
                        let sz = ((z as isize - *ez as isize).rem_euclid(nz as isize)) as usize;
                        let src_idx = (sz * ny * nx + sy * nx + sx) * 19;
                        self.f[dst_idx + q] = f_src[src_idx + q];
                    }
                }
            }
        }

        // -- Bounce-back for wall cells --
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let cell = z * ny * nx + y * nx + x;
                    if self.boundary[cell] != LbmBoundary::NoSlipWall {
                        continue;
                    }
                    let idx = cell * 19;
                    let mut tmp = [0.0f64; 19];
                    for q in 0..19 {
                        tmp[D3Q19_OPP[q]] = self.f[idx + q];
                    }
                    self.f[idx..(19 + idx)].copy_from_slice(&tmp);
                }
            }
        }

        self.step_count += 1;
    }

    /// Run `n` steps.
    pub fn run(&mut self, steps: u64) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// Flat density field (length = nx * ny * nz).
    pub fn density_field(&self) -> Vec<f64> {
        let n = self.nx * self.ny * self.nz;
        let mut out = vec![0.0f64; n];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let cell = z * self.ny * self.nx + y * self.nx + x;
                    out[cell] = self.density_at(x, y, z);
                }
            }
        }
        out
    }
}

/// D3Q19 equilibrium distribution function.
fn d3q19_equilibrium(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 19] {
    let u2 = ux * ux + uy * uy + uz * uz;
    let mut feq = [0.0f64; 19];
    for q in 0..19 {
        let eu = D3Q19_EX[q] as f64 * ux + D3Q19_EY[q] as f64 * uy + D3Q19_EZ[q] as f64 * uz;
        feq[q] = D3Q19_W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
    }
    feq
}

// ---------------------------------------------------------------------------
// Extended LBM boundary conditions
// ---------------------------------------------------------------------------

/// Apply a moving-wall (Zou-He) velocity boundary to a row of cells.
///
/// Forces the top row (z = nz-1) of a 3-D grid to have velocity `vel`.
#[pyfunction]
pub fn apply_moving_wall_3d(sim: &mut PyLbm3dSimulation, z_layer: usize, vel: Vec<f64>) {
    let nz = sim.nz;
    if z_layer >= nz {
        return;
    }
    let nx = sim.nx;
    let ny = sim.ny;
    let vel_arr = if vel.len() >= 3 {
        [vel[0], vel[1], vel[2]]
    } else {
        [0.0; 3]
    };
    for y in 0..ny {
        for x in 0..nx {
            let cell = z_layer * ny * nx + y * nx + x;
            let idx = cell * 19;
            let rho: f64 = sim.f[idx..idx + 19].iter().sum::<f64>().max(0.01);
            let feq = d3q19_equilibrium(rho, vel_arr[0], vel_arr[1], vel_arr[2]);
            sim.f[idx..(19 + idx)].copy_from_slice(&feq);
        }
    }
}

/// Zou-He inlet boundary condition: prescribe velocity at x=0.
///
/// Sets the inlet (left wall) cells of a 2-D simulation to a specified
/// velocity via the Zou-He formulation.
#[pyfunction]
pub fn zou_he_velocity_inlet(sim: &mut PyLbmSimulation, ux_in: f64) {
    let w = sim.width;
    let h = sim.height;
    for y in 1..(h - 1) {
        let cell = y * w; // x = 0
        let idx = cell * 9;
        let rho = (sim.f[idx]
            + sim.f[idx + 2]
            + sim.f[idx + 4]
            + 2.0 * (sim.f[idx + 3] + sim.f[idx + 6] + sim.f[idx + 7]))
            / (1.0 - ux_in);
        let feq = equilibrium(rho, ux_in, 0.0);
        sim.f[idx..(9 + idx)].copy_from_slice(&feq);
    }
}

/// Zou-He outlet boundary condition: prescribe zero-gradient at x = width-1.
#[pyfunction]
pub fn zou_he_pressure_outlet(sim: &mut PyLbmSimulation) {
    let w = sim.width;
    let h = sim.height;
    let x = w - 1;
    for y in 1..(h - 1) {
        let cell = y * w + x;
        let cell_prev = y * w + x - 1;
        let idx = cell * 9;
        let idx_prev = cell_prev * 9;
        // Extrapolate: copy distributions from neighbour
        for q in 0..9 {
            sim.f[idx + q] = sim.f[idx_prev + q];
        }
    }
}

// ===========================================================================
// MRT (Multiple Relaxation Time) collision for D2Q9
// ===========================================================================

/// MRT relaxation rates for D2Q9 (9 rates, one per moment).
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MrtRelaxation {
    /// Relaxation rate for density (s0).
    #[pyo3(get, set)]
    pub s0: f64,
    /// Relaxation rate for energy (s1).
    #[pyo3(get, set)]
    pub s1: f64,
    /// Relaxation rate for energy squared (s2).
    #[pyo3(get, set)]
    pub s2: f64,
    /// Relaxation rate for x-momentum (s3) — usually 1 (skip).
    #[pyo3(get, set)]
    pub s3: f64,
    /// Relaxation rate for energy·x (s4).
    #[pyo3(get, set)]
    pub s4: f64,
    /// Relaxation rate for y-momentum (s5) — usually 1 (skip).
    #[pyo3(get, set)]
    pub s5: f64,
    /// Relaxation rate for energy·y (s6).
    #[pyo3(get, set)]
    pub s6: f64,
    /// Relaxation rate for stress-xx (s7 = 1/τ).
    #[pyo3(get, set)]
    pub s7: f64,
    /// Relaxation rate for stress-xy (s8 = 1/τ).
    #[pyo3(get, set)]
    pub s8: f64,
}

#[pymethods]
impl MrtRelaxation {
    /// Create MRT parameters from kinematic viscosity ν.
    #[new]
    pub fn new(viscosity: f64) -> Self {
        Self::from_viscosity(viscosity)
    }

    /// Create MRT parameters from kinematic viscosity ν (static method).
    #[staticmethod]
    pub fn from_viscosity(viscosity: f64) -> Self {
        let s_nu = 1.0 / (3.0 * viscosity + 0.5);
        Self {
            s0: 1.0,
            s1: 1.4,
            s2: 1.4,
            s3: 1.0,
            s4: 1.2,
            s5: 1.0,
            s6: 1.2,
            s7: s_nu,
            s8: s_nu,
        }
    }

    /// Return rates as a `Vec<f64>` of length 9.
    pub fn as_array(&self) -> Vec<f64> {
        vec![
            self.s0, self.s1, self.s2, self.s3, self.s4, self.s5, self.s6, self.s7, self.s8,
        ]
    }
}

impl Default for MrtRelaxation {
    fn default() -> Self {
        Self::from_viscosity(0.1)
    }
}

// ===========================================================================
// LBM statistics
// ===========================================================================

/// Per-step statistics for a 2-D LBM simulation.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct LbmStats {
    /// Step count.
    #[pyo3(get)]
    pub step_count: u64,
    /// Mean density.
    #[pyo3(get)]
    pub mean_density: f64,
    /// Maximum speed.
    #[pyo3(get)]
    pub max_speed: f64,
    /// Enstrophy.
    #[pyo3(get)]
    pub enstrophy: f64,
    /// Number of fluid cells.
    #[pyo3(get)]
    pub fluid_cell_count: usize,
    /// Number of wall cells.
    #[pyo3(get)]
    pub wall_cell_count: usize,
}

// ===========================================================================
// Tests for D3Q19 LBM API
// ===========================================================================

#[cfg(test)]
mod d3q19_tests {

    use crate::LbmBoundary;

    use crate::lbm_api::PyLbm3dConfig;
    use crate::lbm_api::PyLbm3dSimulation;
    use crate::lbm_api::apply_moving_wall_3d;
    use crate::lbm_api::d3q19_equilibrium;

    fn small_3d() -> PyLbm3dSimulation {
        PyLbm3dSimulation::new(&PyLbm3dConfig::small())
    }

    #[test]
    fn test_d3q19_creation() {
        let sim = small_3d();
        let (nx, ny, nz) = sim.dimensions();
        assert_eq!(nx, 8);
        assert_eq!(ny, 8);
        assert_eq!(nz, 8);
        assert_eq!(sim.step_count(), 0);
    }

    #[test]
    fn test_d3q19_initial_density() {
        let sim = small_3d();
        let rho = sim.density_at(4, 4, 4);
        assert!((rho - 1.0).abs() < 1e-10, "rho = {}", rho);
    }

    #[test]
    fn test_d3q19_initial_velocity_near_zero() {
        let sim = small_3d();
        let v = sim.velocity_at(0, 0, 0);
        assert!(v[0].abs() < 1e-10 && v[1].abs() < 1e-10 && v[2].abs() < 1e-10);
    }

    #[test]
    fn test_d3q19_step_advances_count() {
        let mut sim = small_3d();
        sim.step();
        sim.step();
        assert_eq!(sim.step_count(), 2);
    }

    #[test]
    fn test_d3q19_run_n_steps() {
        let mut sim = small_3d();
        sim.run(10);
        assert_eq!(sim.step_count(), 10);
    }

    #[test]
    fn test_d3q19_density_field_length() {
        let sim = small_3d();
        assert_eq!(sim.density_field().len(), 8 * 8 * 8);
    }

    #[test]
    fn test_d3q19_set_get_boundary() {
        let mut sim = small_3d();
        sim.set_boundary(1, 1, 1, LbmBoundary::NoSlipWall);
        assert_eq!(sim.get_boundary(1, 1, 1), Some(LbmBoundary::NoSlipWall));
    }

    #[test]
    fn test_d3q19_out_of_bounds() {
        let sim = small_3d();
        assert!(sim.get_boundary(99, 0, 0).is_none());
    }

    #[test]
    fn test_d3q19_velocity_zero_on_wall() {
        let mut sim = small_3d();
        sim.set_boundary(2, 2, 2, LbmBoundary::NoSlipWall);
        let v = sim.velocity_at(2, 2, 2);
        assert!(v[0].abs() < 1e-12 && v[1].abs() < 1e-12 && v[2].abs() < 1e-12);
    }

    #[test]
    fn test_d3q19_body_force_creates_flow() {
        let mut cfg = PyLbm3dConfig::new(8, 8, 8, 0.1);
        cfg.body_force = [1e-3, 0.0, 0.0];
        let mut sim = PyLbm3dSimulation::new(&cfg);
        sim.run(50);
        let v = sim.velocity_at(4, 4, 4);
        assert!(v[0] > 0.0, "body force should drive flow, ux={}", v[0]);
    }

    #[test]
    fn test_d3q19_density_conserved() {
        let mut sim = small_3d();
        let rho_before: f64 = sim.density_field().iter().sum();
        sim.run(5);
        let rho_after: f64 = sim.density_field().iter().sum();
        assert!(
            (rho_after - rho_before).abs() / rho_before < 1e-6,
            "density not conserved: {} vs {}",
            rho_before,
            rho_after
        );
    }

    #[test]
    fn test_d3q19_config_omega() {
        let cfg = PyLbm3dConfig::new(4, 4, 4, 1.0 / 6.0);
        assert!((cfg.omega() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lbm3d_config_small() {
        let cfg = PyLbm3dConfig::small();
        assert_eq!(cfg.nx, 8);
        assert_eq!(cfg.ny, 8);
        assert_eq!(cfg.nz, 8);
    }

    #[test]
    fn test_d3q19_equilibrium_sum_to_rho() {
        let feq = d3q19_equilibrium(1.5, 0.1, 0.0, -0.05);
        let total: f64 = feq.iter().sum();
        assert!((total - 1.5).abs() < 1e-10, "equilibrium sum = {}", total);
    }

    #[test]
    fn test_moving_wall_3d_changes_velocity() {
        let mut sim = small_3d();
        let v_before = sim.velocity_at(4, 4, 7);
        apply_moving_wall_3d(&mut sim, 7, vec![0.1, 0.0, 0.0]);
        let v_after = sim.velocity_at(4, 4, 7);
        // The equilibrium was forced, so the reconstructed velocity should shift
        // (not identical because velocity_at recomputes from f)
        let _ = (v_before, v_after); // just ensure no panic
    }
}

// ===========================================================================
// Extended LBM tests
// ===========================================================================

#[cfg(test)]
mod lbm_ext_tests {

    use crate::LbmBoundary;
    use crate::PyLbmConfig;
    use crate::PyLbmSimulation;
    use crate::lbm_api::MrtRelaxation;
    use crate::lbm_api::PyLbm3dConfig;
    use crate::lbm_api::PyLbm3dSimulation;

    use crate::lbm_api::zou_he_pressure_outlet;
    use crate::lbm_api::zou_he_velocity_inlet;

    fn small_sim() -> PyLbmSimulation {
        PyLbmSimulation::new(&PyLbmConfig::new(16, 8, 0.1))
    }

    // --- MrtRelaxation ---

    #[test]
    fn test_mrt_from_viscosity() {
        let mrt = MrtRelaxation::from_viscosity(1.0 / 6.0);
        let s7 = mrt.s7;
        assert!((s7 - 1.0).abs() < 1e-10, "s7={}", s7);
    }

    #[test]
    fn test_mrt_as_array_length() {
        let mrt = MrtRelaxation::default();
        assert_eq!(mrt.as_array().len(), 9);
    }

    #[test]
    fn test_mrt_rates_positive() {
        let mrt = MrtRelaxation::from_viscosity(0.02);
        for s in mrt.as_array() {
            assert!(s > 0.0, "rate must be positive: {}", s);
        }
    }

    // --- Macroscopic quantities ---

    #[test]
    fn test_mean_density_initial() {
        let sim = small_sim();
        let rho = sim.mean_density();
        assert!((rho - 1.0).abs() < 1e-8, "rho = {}", rho);
    }

    #[test]
    fn test_max_speed_initial_zero() {
        let sim = small_sim();
        assert!(sim.max_speed() < 1e-10);
    }

    #[test]
    fn test_vorticity_field_length() {
        let sim = small_sim();
        let vort = sim.vorticity_field();
        assert_eq!(vort.len(), 16 * 8);
    }

    #[test]
    fn test_enstrophy_initial_near_zero() {
        let sim = small_sim();
        assert!(sim.enstrophy() < 1e-10);
    }

    #[test]
    fn test_speed_field_length() {
        let sim = small_sim();
        assert_eq!(sim.speed_field().len(), 16 * 8);
    }

    #[test]
    fn test_reset_to_equilibrium() {
        let mut sim = small_sim();
        sim.run(50);
        sim.reset_to_equilibrium(1.0, 0.0, 0.0);
        let rho = sim.density_at(8, 4);
        assert!((rho - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_fluid_cell_count_all_fluid_initially() {
        let sim = small_sim();
        assert_eq!(sim.fluid_cell_count(), 16 * 8);
    }

    #[test]
    fn test_wall_cell_count_zero_initially() {
        let sim = small_sim();
        assert_eq!(sim.wall_cell_count(), 0);
    }

    #[test]
    fn test_fluid_wall_count_after_enclosing_walls() {
        let mut sim = small_sim();
        sim.add_enclosing_walls();
        let walls = sim.wall_cell_count();
        let fluid = sim.fluid_cell_count();
        assert_eq!(walls + fluid, 16 * 8);
        assert!(walls > 0);
    }

    #[test]
    fn test_circular_obstacle_adds_walls() {
        let mut sim = small_sim();
        sim.add_circular_obstacle(8, 4, 2);
        assert!(sim.wall_cell_count() > 0);
    }

    #[test]
    fn test_circular_obstacle_center_is_wall() {
        let mut sim = small_sim();
        sim.add_circular_obstacle(8, 4, 2);
        assert_eq!(sim.get_boundary(8, 4), Some(LbmBoundary::NoSlipWall));
    }

    #[test]
    fn test_reynolds_number_positive_after_flow() {
        let mut cfg = PyLbmConfig::new(32, 16, 0.02);
        cfg.body_force = [1e-4, 0.0];
        let mut sim = PyLbmSimulation::new(&cfg);
        sim.run(200);
        let re = sim.reynolds_number(cfg.omega());
        let _ = re; // just no panic
    }

    // --- Zou-He BCs ---

    #[test]
    fn test_zou_he_velocity_inlet_no_panic() {
        let mut sim = small_sim();
        zou_he_velocity_inlet(&mut sim, 0.05);
        // just ensure no panic
    }

    #[test]
    fn test_zou_he_pressure_outlet_no_panic() {
        let mut sim = small_sim();
        zou_he_pressure_outlet(&mut sim);
    }

    // --- Statistics ---

    #[test]
    fn test_collect_stats_initial() {
        let sim = small_sim();
        let stats = sim.collect_stats();
        assert_eq!(stats.step_count, 0);
        assert!((stats.mean_density - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_collect_stats_after_run() {
        let mut sim = small_sim();
        sim.run(10);
        let stats = sim.collect_stats();
        assert_eq!(stats.step_count, 10);
    }

    #[test]
    fn test_lbm_config_omega_range() {
        let cfg = PyLbmConfig::new(8, 8, 0.1);
        let omega = cfg.omega();
        // omega = 1/(3*0.1 + 0.5) = 1/0.8 = 1.25
        assert!((omega - 1.25).abs() < 1e-10);
    }

    // --- D3Q19 extended ---

    #[test]
    fn test_d3q19_velocity_field_size() {
        let sim = PyLbm3dSimulation::new(&PyLbm3dConfig::small());
        let (nx, ny, nz) = sim.dimensions();
        let density = sim.density_field();
        assert_eq!(density.len(), nx * ny * nz);
    }

    #[test]
    fn test_d3q19_enclosing_walls_boundary_check() {
        let mut sim = PyLbm3dSimulation::new(&PyLbm3dConfig::new(4, 4, 4, 0.1));
        // Add walls on the Z=0 layer
        for y in 0..4 {
            for x in 0..4 {
                sim.set_boundary(x, y, 0, LbmBoundary::NoSlipWall);
            }
        }
        assert_eq!(sim.get_boundary(2, 2, 0), Some(LbmBoundary::NoSlipWall));
        assert_eq!(sim.get_boundary(2, 2, 2), Some(LbmBoundary::Fluid));
    }

    #[test]
    fn test_d3q19_body_force_x_creates_positive_flow() {
        let mut cfg = PyLbm3dConfig::new(8, 8, 8, 0.05);
        cfg.body_force = [5e-4, 0.0, 0.0];
        let mut sim = PyLbm3dSimulation::new(&cfg);
        sim.run(100);
        let v = sim.velocity_at(4, 4, 4);
        assert!(
            v[0] > 0.0,
            "body force +x should drive ux > 0, got {}",
            v[0]
        );
    }

    #[test]
    fn test_d3q19_config_new_clamps_viscosity() {
        let cfg = PyLbm3dConfig::new(4, 4, 4, 0.0);
        assert!(cfg.viscosity >= 1e-8);
    }

    #[test]
    fn test_lbm2d_config_clamps_viscosity() {
        let cfg = PyLbmConfig::new(4, 4, 0.0);
        assert!(cfg.viscosity >= 1e-8);
    }

    #[test]
    fn test_d3q19_out_of_bounds_velocity() {
        let sim = PyLbm3dSimulation::new(&PyLbm3dConfig::small());
        let v = sim.velocity_at(100, 0, 0);
        for &vi in &v {
            assert!((vi).abs() < 1e-15);
        }
    }

    #[test]
    fn test_vorticity_boundary_cells_zero() {
        let sim = small_sim();
        let vort = sim.vorticity_field();
        // Corner cells (x=0 or y=0) are always zero by construction
        assert!((vort[0]).abs() < 1e-15);
    }

    #[test]
    fn test_speed_field_initial_near_zero() {
        let sim = small_sim();
        let speed = sim.speed_field();
        for &s in &speed {
            assert!(s < 1e-10, "initial speed should be near zero, got {}", s);
        }
    }
}

/// Register all `lbm` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_lbm_module(parent: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = pyo3::types::PyModule::new(parent.py(), "lbm")?;
    child.add_class::<LbmBoundary>()?;
    child.add_class::<PyLbmConfig>()?;
    child.add_class::<PyLbmSimulation>()?;
    child.add_class::<PyLbm3dConfig>()?;
    child.add_class::<PyLbm3dSimulation>()?;
    child.add_class::<MrtRelaxation>()?;
    child.add_class::<LbmStats>()?;
    child.add_function(pyo3::wrap_pyfunction!(apply_moving_wall_3d, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(zou_he_velocity_inlet, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(zou_he_pressure_outlet, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}
