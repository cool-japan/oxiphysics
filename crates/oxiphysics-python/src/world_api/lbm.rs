// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM (Lattice Boltzmann Method) Simulation.

#![allow(missing_docs)]

// ===========================================================================
// LBM (Lattice Boltzmann Method) Simulation
// ===========================================================================

/// Configuration for a 2-D D2Q9 Lattice-Boltzmann simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyLbmConfig {
    /// Grid width (number of cells in X direction).
    pub width: usize,
    /// Grid height (number of cells in Y direction).
    pub height: usize,
    /// Kinematic viscosity (nu). Controls the relaxation rate.
    pub viscosity: f64,
    /// External body-force acceleration `[fx, fy]` applied to the fluid.
    pub body_force: [f64; 2],
}

impl PyLbmConfig {
    /// Create a new LBM configuration.
    pub fn new(width: usize, height: usize, viscosity: f64) -> Self {
        Self {
            width,
            height,
            viscosity: viscosity.max(1e-6),
            body_force: [0.0; 2],
        }
    }

    /// Create a default lid-driven cavity configuration (64×64).
    pub fn lid_driven_cavity() -> Self {
        Self::new(64, 64, 0.01)
    }

    /// Compute the relaxation parameter omega from viscosity.
    ///
    /// `omega = 1 / (3 * nu + 0.5)`
    pub fn omega(&self) -> f64 {
        1.0 / (3.0 * self.viscosity + 0.5)
    }
}

/// D2Q9 Lattice-Boltzmann grid simulation.
///
/// Stores distribution functions for all 9 velocity directions at every cell.
/// Uses the BGK collision operator with a single relaxation time.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyLbmGrid {
    /// Number of cells in X.
    width: usize,
    /// Number of cells in Y.
    height: usize,
    /// Relaxation rate omega.
    omega: f64,
    /// External body force `[fx, fy]`.
    body_force: [f64; 2],
    /// Distribution functions, indexed `[cell_idx * 9 + q]`.
    f: Vec<f64>,
    /// Post-collision / swap buffer.
    f_tmp: Vec<f64>,
    /// Accumulated simulation steps.
    step_count: u64,
}

// D2Q9 weights and discrete velocities
const D2Q9_W: [f64; 9] = [
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
// ex, ey for each of the 9 directions
const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

impl PyLbmGrid {
    /// Create a new LBM grid from configuration. All cells initialised to
    /// equilibrium with unit density and zero velocity.
    pub fn new(config: &PyLbmConfig) -> Self {
        let n = config.width * config.height;
        let mut f = vec![0.0f64; n * 9];
        // Initialise each cell to equilibrium at rho=1, u=(0,0)
        for i in 0..n {
            for q in 0..9 {
                f[i * 9 + q] = D2Q9_W[q];
            }
        }
        Self {
            width: config.width,
            height: config.height,
            omega: config.omega(),
            body_force: config.body_force,
            f: f.clone(),
            f_tmp: f,
            step_count: 0,
        }
    }

    /// Return the grid width.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Return the grid height.
    pub fn height(&self) -> usize {
        self.height
    }

    /// Return the number of completed steps.
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Advance the simulation by one LBM time step.
    ///
    /// Performs: collision (BGK) + streaming + periodic bounce-back walls.
    pub fn step(&mut self) {
        let w = self.width;
        let h = self.height;
        let omega = self.omega;
        let bfx = self.body_force[0];
        let bfy = self.body_force[1];

        // --- Collision + compute macroscopic fields ---
        for y in 0..h {
            for x in 0..w {
                let idx = (y * w + x) * 9;
                // Density
                let rho: f64 = self.f[idx..idx + 9].iter().sum();
                // Velocity
                let mut ux = 0.0f64;
                let mut uy = 0.0f64;
                for q in 0..9 {
                    ux += D2Q9_EX[q] * self.f[idx + q];
                    uy += D2Q9_EY[q] * self.f[idx + q];
                }
                let rho_inv = if rho > 1e-15 { 1.0 / rho } else { 0.0 };
                ux = ux * rho_inv + bfx * 0.5 / omega;
                uy = uy * rho_inv + bfy * 0.5 / omega;
                let u2 = ux * ux + uy * uy;
                // BGK collision → post-collision stored in f_tmp
                for q in 0..9 {
                    let eu = D2Q9_EX[q] * ux + D2Q9_EY[q] * uy;
                    let feq = D2Q9_W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
                    self.f_tmp[idx + q] = self.f[idx + q] * (1.0 - omega) + feq * omega;
                }
            }
        }

        // --- Streaming with periodic boundaries ---
        let f_src = self.f_tmp.clone();
        for y in 0..h {
            for x in 0..w {
                let dst_idx = (y * w + x) * 9;
                #[allow(clippy::manual_memcpy)]
                for q in 0..9 {
                    let src_x =
                        ((x as isize - D2Q9_EX[q] as isize).rem_euclid(w as isize)) as usize;
                    let src_y =
                        ((y as isize - D2Q9_EY[q] as isize).rem_euclid(h as isize)) as usize;
                    let src_idx = (src_y * w + src_x) * 9;
                    self.f[dst_idx + q] = f_src[src_idx + q];
                }
            }
        }

        self.step_count += 1;
    }

    /// Get the macroscopic velocity `[ux, uy]` at cell `(x, y)`.
    ///
    /// Returns `[0.0, 0.0]` if coordinates are out of bounds.
    pub fn velocity_at(&self, x: usize, y: usize) -> [f64; 2] {
        if x >= self.width || y >= self.height {
            return [0.0; 2];
        }
        let idx = (y * self.width + x) * 9;
        let rho: f64 = self.f[idx..idx + 9].iter().sum();
        if rho < 1e-15 {
            return [0.0; 2];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        for q in 0..9 {
            ux += D2Q9_EX[q] * self.f[idx + q];
            uy += D2Q9_EY[q] * self.f[idx + q];
        }
        [ux / rho, uy / rho]
    }

    /// Get the macroscopic density at cell `(x, y)`.
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        if x >= self.width || y >= self.height {
            return 0.0;
        }
        let idx = (y * self.width + x) * 9;
        self.f[idx..idx + 9].iter().sum()
    }

    /// Return the entire velocity field as a flat `Vec`f64` of `\[ux, uy\]` pairs,
    /// row-major order (x-major).
    pub fn velocity_field(&self) -> Vec<f64> {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n * 2];
        for y in 0..self.height {
            for x in 0..self.width {
                let cell = y * self.width + x;
                let v = self.velocity_at(x, y);
                out[cell * 2] = v[0];
                out[cell * 2 + 1] = v[1];
            }
        }
        out
    }

    /// Return the entire density field as a flat `Vec`f64`, row-major.
    pub fn density_field(&self) -> Vec<f64> {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n];
        for y in 0..self.height {
            for x in 0..self.width {
                out[y * self.width + x] = self.density_at(x, y);
            }
        }
        out
    }
}
