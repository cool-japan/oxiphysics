// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Concentration fields and basic reaction functions.

// ---------------------------------------------------------------------------
// ConcentrationField
// ---------------------------------------------------------------------------

/// 2D concentration field with explicit Fickian diffusion and periodic BC.
pub struct ConcentrationField {
    /// Domain width (number of cells in x).
    pub nx: usize,
    /// Domain height (number of cells in y).
    pub ny: usize,
    /// Concentration data, row-major: `data[j * nx + i]`.
    pub data: Vec<f64>,
    /// Molecular diffusivity D.
    pub diffusivity: f64,
}

impl ConcentrationField {
    /// Create a new field initialised to a uniform concentration `init`.
    pub fn new(nx: usize, ny: usize, diffusivity: f64, init: f64) -> Self {
        Self {
            nx,
            ny,
            data: vec![init; nx * ny],
            diffusivity,
        }
    }

    /// Get the concentration at grid cell (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.nx + i]
    }

    /// Set the concentration at grid cell (i, j).
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[j * self.nx + i] = v;
    }

    /// Advance the concentration field by one explicit diffusion step.
    ///
    /// Uses the standard second-order finite difference Laplacian with
    /// periodic boundary conditions:
    ///
    /// C_new(i,j) = C(i,j) + D * dt / dx^2 * (C(i+1,j) + C(i-1,j)
    ///                                        + C(i,j+1) + C(i,j-1)
    ///                                        - 4*C(i,j))
    pub fn diffusion_step(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let coeff = self.diffusivity * dt / (dx * dx);
        let old = self.data.clone();

        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;

                let laplacian =
                    old[j * nx + ip] + old[j * nx + im] + old[jp * nx + i] + old[jm * nx + i]
                        - 4.0 * old[j * nx + i];

                self.data[j * nx + i] = old[j * nx + i] + coeff * laplacian;
            }
        }
    }

    /// Return the sum of all concentration values (total mass / dx^2).
    pub fn total_mass(&self) -> f64 {
        self.data.iter().sum()
    }

    /// Return the maximum concentration in the field.
    pub fn max_concentration(&self) -> f64 {
        self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Return the minimum concentration in the field.
    pub fn min_concentration(&self) -> f64 {
        self.data.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Apply a floor so that all concentrations are at least `min_val`.
    pub fn clamp_min(&mut self, min_val: f64) {
        for v in self.data.iter_mut() {
            if *v < min_val {
                *v = min_val;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Reaction functions
// ---------------------------------------------------------------------------

/// Apply a first-order decay reaction: C(t+dt) = C(t) * exp(-k * dt).
///
/// This is the exact (operator-split) solution to dC/dt = -k*C.
pub fn first_order_reaction(c: &mut ConcentrationField, rate_constant: f64, dt: f64) {
    let factor = (-rate_constant * dt).exp();
    for v in c.data.iter_mut() {
        *v *= factor;
    }
}

/// Apply a bimolecular reaction A + B -> products using explicit Euler.
///
/// delta_A = delta_B = -rate * A * B * dt  (applied pointwise, same grid).
///
/// Both `a` and `b` must have the same dimensions.
pub fn bimolecular_reaction(
    a: &mut ConcentrationField,
    b: &mut ConcentrationField,
    rate: f64,
    dt: f64,
) {
    debug_assert_eq!(a.nx, b.nx);
    debug_assert_eq!(a.ny, b.ny);
    let n = a.nx * a.ny;
    for k in 0..n {
        let delta = rate * a.data[k] * b.data[k] * dt;
        a.data[k] -= delta;
        b.data[k] -= delta;
    }
}

/// Apply a second-order reaction A + B -> C using explicit Euler.
///
/// delta_A = delta_B = -rate * A * B * dt, delta_C = +rate * A * B * dt.
pub fn bimolecular_reaction_with_product(
    a: &mut ConcentrationField,
    b: &mut ConcentrationField,
    c: &mut ConcentrationField,
    rate: f64,
    dt: f64,
) {
    debug_assert_eq!(a.nx, b.nx);
    debug_assert_eq!(a.ny, b.ny);
    debug_assert_eq!(a.nx, c.nx);
    debug_assert_eq!(a.ny, c.ny);
    let n = a.nx * a.ny;
    for k in 0..n {
        let delta = rate * a.data[k] * b.data[k] * dt;
        a.data[k] -= delta;
        b.data[k] -= delta;
        c.data[k] += delta;
    }
}

/// Apply a reversible reaction A <-> B with forward rate kf and reverse rate kr.
///
/// dA/dt = -kf*A + kr*B,  dB/dt = +kf*A - kr*B
pub fn reversible_reaction(
    a: &mut ConcentrationField,
    b: &mut ConcentrationField,
    kf: f64,
    kr: f64,
    dt: f64,
) {
    debug_assert_eq!(a.nx, b.nx);
    debug_assert_eq!(a.ny, b.ny);
    let n = a.nx * a.ny;
    for k in 0..n {
        let forward = kf * a.data[k] * dt;
        let reverse = kr * b.data[k] * dt;
        a.data[k] += -forward + reverse;
        b.data[k] += forward - reverse;
    }
}
