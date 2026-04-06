// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Configuration types and state for molecular dynamics simulation.
//!
//! Contains the physical constants, ensemble definitions, configuration
//! parameters, and simulation state types.

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant in reduced / amu · Å² · ps⁻² units (kJ mol⁻¹ K⁻¹).
#[allow(dead_code)]
pub const KB_REDUCED: f64 = 8.314_462_618e-3; // kJ mol⁻¹ K⁻¹

// ---------------------------------------------------------------------------
// Ensemble enum
// ---------------------------------------------------------------------------

/// Statistical-mechanical ensemble for the MD run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum Ensemble {
    /// Micro-canonical: constant N, V, E.
    NVE,
    /// Canonical: constant N, V, T (thermostat active).
    NVT,
    /// Isothermal-isobaric: constant N, P, T (thermostat + barostat active).
    NPT,
}

// ---------------------------------------------------------------------------
// MdConfig
// ---------------------------------------------------------------------------

/// Configuration parameters for a self-contained MD run.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MdConfig {
    /// Time step (ps).
    pub dt: f64,
    /// Total number of MD steps.
    pub n_steps: u64,
    /// Target temperature (K).
    pub temperature: f64,
    /// Target pressure (bar).
    pub pressure: f64,
    /// Orthorhombic simulation box lengths \[Lx, Ly, Lz\] (Å).
    pub box_lengths: [f64; 3],
    /// Apply periodic boundary conditions.
    pub pbc: bool,
    /// Statistical ensemble.
    pub ensemble: Ensemble,
}

impl Default for MdConfig {
    fn default() -> Self {
        Self {
            dt: 0.002,
            n_steps: 1000,
            temperature: 300.0,
            pressure: 1.0,
            box_lengths: [30.0, 30.0, 30.0],
            pbc: true,
            ensemble: Ensemble::NVT,
        }
    }
}

// ---------------------------------------------------------------------------
// MdState
// ---------------------------------------------------------------------------

/// Dynamical state of a self-contained MD system.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MdState {
    /// Atom positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities (Å ps⁻¹).
    pub velocities: Vec<[f64; 3]>,
    /// Forces on atoms (kJ mol⁻¹ Å⁻¹).
    pub forces: Vec<[f64; 3]>,
    /// Atom masses (amu).
    pub masses: Vec<f64>,
    /// Current MD step index.
    pub step: u64,
    /// Simulated time (ps).
    pub time: f64,
    /// Instantaneous temperature (K).
    pub temperature: f64,
    /// Instantaneous pressure (bar).
    pub pressure: f64,
    /// Potential energy (kJ mol⁻¹).
    pub potential_energy: f64,
    /// Kinetic energy (kJ mol⁻¹).
    pub kinetic_energy: f64,
}

impl MdState {
    /// Create a new state from atom data (forces initialised to zero).
    #[allow(dead_code)]
    pub fn new(positions: Vec<[f64; 3]>, velocities: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        let n = positions.len();
        assert_eq!(velocities.len(), n, "velocity count mismatch");
        assert_eq!(masses.len(), n, "mass count mismatch");
        Self {
            forces: vec![[0.0; 3]; n],
            positions,
            velocities,
            masses,
            step: 0,
            time: 0.0,
            temperature: 0.0,
            pressure: 0.0,
            potential_energy: 0.0,
            kinetic_energy: 0.0,
        }
    }

    /// Number of atoms.
    #[allow(dead_code)]
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }
}
