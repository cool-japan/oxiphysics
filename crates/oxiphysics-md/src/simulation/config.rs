// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Configuration types and state for molecular dynamics simulation.
//!
//! Contains the physical constants, ensemble definitions, configuration
//! parameters, and simulation state types.
//!
//! # Unit system
//!
//! The self-contained drivers of this module ([`MdSim`](super::MdSim),
//! [`MdDriver`](super::MdDriver), [`MdSimulationPlain`](super::MdSimulationPlain))
//! work in the common biomolecular MD unit system:
//!
//! | quantity    | unit                         |
//! |-------------|------------------------------|
//! | length      | Å                            |
//! | time        | ps                           |
//! | mass        | amu                          |
//! | energy      | kJ mol⁻¹                     |
//! | force       | kJ mol⁻¹ Å⁻¹                 |
//! | velocity    | Å ps⁻¹                       |
//! | temperature | K                            |
//! | pressure    | bar                          |
//!
//! The mechanical combination amu · Å² · ps⁻² is **not** the energy unit
//! kJ mol⁻¹.  With 1 amu = 1.660 539 × 10⁻²⁷ kg one amu · Å² · ps⁻² is
//! 1.660 539 × 10⁻²³ J per particle, i.e. 10 J mol⁻¹:
//!
//! ```text
//! 1 amu · Å² · ps⁻²  =  0.01 kJ mol⁻¹          (KJMOL_PER_AMU_A2_PS2)
//! 1 kJ mol⁻¹         =  100 amu · Å² · ps⁻²    (ACCEL_PER_FORCE_OVER_MASS)
//! ```
//!
//! Consequently a force `F` (kJ mol⁻¹ Å⁻¹) on a mass `m` (amu) produces the
//! acceleration `100 · F / m` (Å ps⁻²), and a velocity `v` (Å ps⁻¹) carries the
//! kinetic energy `½ m v² · 0.01` (kJ mol⁻¹).  The helpers
//! [`atom_kinetic_energy`], [`velocity_increment`] and
//! [`thermal_velocity_sigma`] apply these factors; every integrator,
//! thermostat and velocity initialiser of the module goes through them, so
//! that time steps, velocities and temperatures are physical (a spring with
//! `k = 100 kJ mol⁻¹ Å⁻²` and `m = 1 amu` oscillates with the period
//! `0.2π sqrt(m / k) ps ≈ 0.0628 ps`).

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant in kJ mol⁻¹ K⁻¹.
///
/// This is the gas constant `R` expressed per mole, the value matching the
/// kJ mol⁻¹ energy unit of the module (see the [module documentation](self)).
/// Multiplying it by a temperature gives an energy in kJ mol⁻¹; to obtain a
/// velocity scale from `k_B T / m` use [`thermal_velocity_sigma`], which
/// includes the conversion [`ACCEL_PER_FORCE_OVER_MASS`].
pub const KB_REDUCED: f64 = 8.314_462_618e-3; // kJ mol⁻¹ K⁻¹

/// Energy in kJ mol⁻¹ of one amu · Å² · ps⁻²: `1 amu Å² ps⁻² = 0.01 kJ mol⁻¹`.
///
/// Multiply `½ m v²` (amu Å² ps⁻²) by this factor to obtain a kinetic energy
/// in kJ mol⁻¹.
pub const KJMOL_PER_AMU_A2_PS2: f64 = 0.01;

/// Acceleration in Å ps⁻² produced by a force of 1 kJ mol⁻¹ Å⁻¹ on a mass of
/// 1 amu: `1 kJ mol⁻¹ = 100 amu Å² ps⁻²`.
///
/// Multiply `F / m` (kJ mol⁻¹ Å⁻¹ amu⁻¹) by this factor to obtain an
/// acceleration in Å ps⁻².  It is the reciprocal of [`KJMOL_PER_AMU_A2_PS2`].
pub const ACCEL_PER_FORCE_OVER_MASS: f64 = 100.0;

/// Kinetic energy `½ m |v|²` (kJ mol⁻¹) of one atom of mass `mass` (amu) moving
/// with `velocity` (Å ps⁻¹).
#[inline]
pub fn atom_kinetic_energy(velocity: &[f64; 3], mass: f64) -> f64 {
    let speed_sq =
        velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
    0.5 * mass * speed_sq * KJMOL_PER_AMU_A2_PS2
}

/// Velocity change `Δv = (F / m) · dt` (Å ps⁻¹) of an atom of mass `mass` (amu)
/// under the constant `force` (kJ mol⁻¹ Å⁻¹) acting for `dt` (ps).
///
/// Includes the unit factor [`ACCEL_PER_FORCE_OVER_MASS`].
#[inline]
pub fn velocity_increment(force: &[f64; 3], mass: f64, dt: f64) -> [f64; 3] {
    let factor = ACCEL_PER_FORCE_OVER_MASS * dt / mass;
    [force[0] * factor, force[1] * factor, force[2] * factor]
}

/// Standard deviation (Å ps⁻¹) of one Cartesian velocity component of an atom
/// of mass `mass` (amu) in thermal equilibrium at `temperature` (K):
/// `σ = sqrt(k_B T / m)` with the unit factor [`ACCEL_PER_FORCE_OVER_MASS`].
///
/// For example argon (39.948 amu) at 120 K has `σ ≈ 1.58 Å ps⁻¹`.
#[inline]
pub fn thermal_velocity_sigma(temperature: f64, mass: f64) -> f64 {
    (ACCEL_PER_FORCE_OVER_MASS * KB_REDUCED * temperature / mass).sqrt()
}

// ---------------------------------------------------------------------------
// Ensemble enum
// ---------------------------------------------------------------------------

/// Statistical-mechanical ensemble for the MD run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Advance every velocity by `F / m · dt` using the stored forces.
    ///
    /// `dt` is the duration of the kick (ps; pass `0.5 * dt` for a velocity-Verlet
    /// half step).  The forces are in kJ mol⁻¹ Å⁻¹, the masses in amu and the
    /// velocities in Å ps⁻¹, so the unit factor
    /// [`ACCEL_PER_FORCE_OVER_MASS`] is included (see [`velocity_increment`]).
    ///
    /// # Panics
    ///
    /// Panics if `forces` or `masses` hold fewer entries than `velocities`.
    pub fn kick_velocities(&mut self, dt: f64) {
        for i in 0..self.n_atoms() {
            let dv = velocity_increment(&self.forces[i], self.masses[i], dt);
            for (v, d) in self.velocities[i].iter_mut().zip(dv.iter()) {
                *v += d;
            }
        }
    }
}
