// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Catalytic surfaces, species boundary conditions, reaction front tracking,
//! and heat release coupling.

use super::concentration::ConcentrationField;

/// A catalytic surface that enhances reaction rates at specific grid cells.
///
/// Cells marked as catalytic have their reaction rate multiplied by
/// `enhancement_factor`.
#[derive(Debug, Clone)]
pub struct CatalyticSurface {
    /// Grid mask: true = catalytic cell.
    pub mask: Vec<bool>,
    /// Enhancement factor for catalytic cells.
    pub enhancement_factor: f64,
    /// Adsorption rate for species uptake onto surface.
    pub adsorption_rate: f64,
    /// Desorption rate for species release from surface.
    pub desorption_rate: f64,
    /// Surface coverage (fraction 0..1) at each catalytic cell.
    pub coverage: Vec<f64>,
}

impl CatalyticSurface {
    /// Create a catalytic surface for a domain of size nx*ny.
    pub fn new(nx: usize, ny: usize, enhancement_factor: f64) -> Self {
        let n = nx * ny;
        Self {
            mask: vec![false; n],
            enhancement_factor,
            adsorption_rate: 0.0,
            desorption_rate: 0.0,
            coverage: vec![0.0; n],
        }
    }

    /// Set adsorption/desorption rates.
    pub fn with_adsorption(mut self, ads_rate: f64, des_rate: f64) -> Self {
        self.adsorption_rate = ads_rate;
        self.desorption_rate = des_rate;
        self
    }

    /// Mark a row of cells as catalytic (e.g., a wall at j=wall_j).
    pub fn set_catalytic_row(&mut self, nx: usize, j: usize) {
        for i in 0..nx {
            self.mask[j * nx + i] = true;
        }
    }

    /// Get the effective rate at a given cell index.
    ///
    /// Returns `base_rate * enhancement_factor` for catalytic cells,
    /// `base_rate` otherwise.
    pub fn effective_rate(&self, cell: usize, base_rate: f64) -> f64 {
        if self.mask[cell] {
            base_rate * self.enhancement_factor
        } else {
            base_rate
        }
    }

    /// Update surface coverage using Langmuir kinetics.
    ///
    /// d(theta)/dt = k_ads * c * (1 - theta) - k_des * theta
    pub fn update_coverage(&mut self, concentration: &[f64], dt: f64) {
        for (idx, is_cat) in self.mask.iter().enumerate() {
            if *is_cat {
                let theta = self.coverage[idx];
                let c = concentration[idx];
                let d_theta =
                    self.adsorption_rate * c * (1.0 - theta) - self.desorption_rate * theta;
                self.coverage[idx] = (theta + d_theta * dt).clamp(0.0, 1.0);
            }
        }
    }

    /// Number of catalytic cells.
    pub fn num_catalytic_cells(&self) -> usize {
        self.mask.iter().filter(|&&v| v).count()
    }
}

// ---------------------------------------------------------------------------
// SpeciesBoundaryCondition
// ---------------------------------------------------------------------------

/// Boundary condition types for species transport.
#[derive(Debug, Clone, Copy)]
pub enum SpeciesBcType {
    /// Fixed concentration (Dirichlet).
    FixedConcentration(f64),
    /// Zero gradient (Neumann).
    ZeroFlux,
    /// Convective outflow: dC/dn = 0 (extrapolation).
    ConvectiveOutflow,
}

/// Boundary condition applied to one edge of the 2D domain.
#[derive(Debug, Clone, Copy)]
pub enum BoundaryEdge {
    /// Left edge (i=0).
    Left,
    /// Right edge (i=nx-1).
    Right,
    /// Bottom edge (j=0).
    Bottom,
    /// Top edge (j=ny-1).
    Top,
}

/// A species boundary condition specification.
#[derive(Debug, Clone, Copy)]
pub struct SpeciesBoundaryCondition {
    /// Which edge to apply this BC on.
    pub edge: BoundaryEdge,
    /// The type of BC.
    pub bc_type: SpeciesBcType,
}

impl SpeciesBoundaryCondition {
    /// Create a new species BC.
    pub fn new(edge: BoundaryEdge, bc_type: SpeciesBcType) -> Self {
        Self { edge, bc_type }
    }

    /// Apply this boundary condition to a concentration field.
    pub fn apply(&self, field: &mut ConcentrationField) {
        let nx = field.nx;
        let ny = field.ny;
        match self.edge {
            BoundaryEdge::Left => {
                for j in 0..ny {
                    match self.bc_type {
                        SpeciesBcType::FixedConcentration(c) => field.set(0, j, c),
                        SpeciesBcType::ZeroFlux => {
                            let val = field.get(1, j);
                            field.set(0, j, val);
                        }
                        SpeciesBcType::ConvectiveOutflow => {
                            let val = field.get(1, j);
                            field.set(0, j, val);
                        }
                    }
                }
            }
            BoundaryEdge::Right => {
                for j in 0..ny {
                    match self.bc_type {
                        SpeciesBcType::FixedConcentration(c) => field.set(nx - 1, j, c),
                        SpeciesBcType::ZeroFlux => {
                            let val = field.get(nx - 2, j);
                            field.set(nx - 1, j, val);
                        }
                        SpeciesBcType::ConvectiveOutflow => {
                            let val = field.get(nx - 2, j);
                            field.set(nx - 1, j, val);
                        }
                    }
                }
            }
            BoundaryEdge::Bottom => {
                for i in 0..nx {
                    match self.bc_type {
                        SpeciesBcType::FixedConcentration(c) => field.set(i, 0, c),
                        SpeciesBcType::ZeroFlux => {
                            let val = field.get(i, 1);
                            field.set(i, 0, val);
                        }
                        SpeciesBcType::ConvectiveOutflow => {
                            let val = field.get(i, 1);
                            field.set(i, 0, val);
                        }
                    }
                }
            }
            BoundaryEdge::Top => {
                for i in 0..nx {
                    match self.bc_type {
                        SpeciesBcType::FixedConcentration(c) => field.set(i, ny - 1, c),
                        SpeciesBcType::ZeroFlux => {
                            let val = field.get(i, ny - 2);
                            field.set(i, ny - 1, val);
                        }
                        SpeciesBcType::ConvectiveOutflow => {
                            let val = field.get(i, ny - 2);
                            field.set(i, ny - 1, val);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ReactionFrontTracker
// ---------------------------------------------------------------------------

/// Tracks the position of a reaction front in a 2D domain.
///
/// The front is defined as the iso-contour where concentration equals a
/// threshold value. The tracker records the average x-position of this
/// front over time.
#[derive(Debug, Clone)]
pub struct ReactionFrontTracker {
    /// Concentration threshold defining the front.
    pub threshold: f64,
    /// History of (time, average_x_position) pairs.
    pub history: Vec<(f64, f64)>,
}

impl ReactionFrontTracker {
    /// Create a new front tracker with the given threshold.
    pub fn new(threshold: f64) -> Self {
        Self {
            threshold,
            history: Vec::new(),
        }
    }

    /// Locate the front and record its position at the given time.
    ///
    /// The front position is the average x-coordinate of cells where
    /// the concentration crosses the threshold.
    pub fn track(&mut self, field: &ConcentrationField, time: f64) {
        let nx = field.nx;
        let ny = field.ny;
        let mut sum_x = 0.0;
        let mut count = 0usize;

        for j in 0..ny {
            for i in 0..(nx - 1) {
                let c_here = field.get(i, j);
                let c_next = field.get(i + 1, j);
                // Check if threshold is between c_here and c_next
                if (c_here - self.threshold) * (c_next - self.threshold) <= 0.0 {
                    // Linear interpolation for crossing position
                    let dc = c_next - c_here;
                    if dc.abs() > 1e-20 {
                        let frac = (self.threshold - c_here) / dc;
                        sum_x += i as f64 + frac;
                    } else {
                        sum_x += i as f64 + 0.5;
                    }
                    count += 1;
                }
            }
        }

        let avg_x = if count > 0 {
            sum_x / count as f64
        } else {
            f64::NAN
        };
        self.history.push((time, avg_x));
    }

    /// Estimate front speed from the last two recorded positions.
    pub fn front_speed(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        let n = self.history.len();
        let (t1, x1) = self.history[n - 2];
        let (t2, x2) = self.history[n - 1];
        let dt = t2 - t1;
        if dt.abs() < 1e-20 {
            return 0.0;
        }
        (x2 - x1) / dt
    }

    /// Number of recorded front positions.
    pub fn num_records(&self) -> usize {
        self.history.len()
    }
}

// ---------------------------------------------------------------------------
// HeatRelease
// ---------------------------------------------------------------------------

/// Coupling between reaction progress and a temperature field.
///
/// When a reaction consumes reactant at rate R, heat is released:
/// dT/dt += heat_of_reaction * R
#[derive(Debug, Clone)]
pub struct HeatRelease {
    /// Heat of reaction (dimensionless or in lattice units).
    /// Positive = exothermic.
    pub heat_of_reaction: f64,
    /// Temperature field (same grid as concentration).
    pub temperature: Vec<f64>,
    /// Thermal diffusivity for temperature diffusion.
    pub thermal_diffusivity: f64,
}

impl HeatRelease {
    /// Create a heat release model for an nx*ny domain at initial temperature T0.
    pub fn new(nx: usize, ny: usize, heat_of_reaction: f64, t0: f64) -> Self {
        Self {
            heat_of_reaction,
            temperature: vec![t0; nx * ny],
            thermal_diffusivity: 0.01,
        }
    }

    /// Set thermal diffusivity.
    pub fn with_thermal_diffusivity(mut self, alpha: f64) -> Self {
        self.thermal_diffusivity = alpha;
        self
    }

    /// Apply heat release from a first-order reaction consuming species c.
    ///
    /// For each cell: dT = heat_of_reaction * k * c * dt
    pub fn apply_first_order(&mut self, concentration: &[f64], rate_constant: f64, dt: f64) {
        for (idx, &c) in concentration.iter().enumerate() {
            self.temperature[idx] += self.heat_of_reaction * rate_constant * c * dt;
        }
    }

    /// Diffuse the temperature field (periodic BC).
    pub fn diffuse_temperature(&mut self, nx: usize, ny: usize, dt: f64, dx: f64) {
        let coeff = self.thermal_diffusivity * dt / (dx * dx);
        let old = self.temperature.clone();
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                let laplacian =
                    old[j * nx + ip] + old[j * nx + im] + old[jp * nx + i] + old[jm * nx + i]
                        - 4.0 * old[j * nx + i];
                self.temperature[j * nx + i] = old[j * nx + i] + coeff * laplacian;
            }
        }
    }

    /// Mean temperature across all cells.
    pub fn mean_temperature(&self) -> f64 {
        if self.temperature.is_empty() {
            return 0.0;
        }
        self.temperature.iter().sum::<f64>() / self.temperature.len() as f64
    }

    /// Maximum temperature.
    pub fn max_temperature(&self) -> f64 {
        self.temperature
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
