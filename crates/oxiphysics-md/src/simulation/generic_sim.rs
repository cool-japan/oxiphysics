// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Generic trait-based MD simulation driver and energy records.

use crate::atom::AtomSet;
use crate::barostat::Barostat;
use crate::forcefield::ForceField;
use crate::integrator::Integrator;
use crate::neighbor::PeriodicBox;
use crate::thermostat::Thermostat;

// ---------------------------------------------------------------------------
// EnergyRecord
// ---------------------------------------------------------------------------

/// Tracks energy components during a simulation.
#[derive(Debug, Clone, Copy, Default)]
pub struct EnergyRecord {
    /// Kinetic energy.
    pub kinetic: f64,
    /// Potential energy.
    pub potential: f64,
    /// Total energy (kinetic + potential).
    pub total: f64,
}

// ---------------------------------------------------------------------------
// Generic MdSimulation (trait-based, original code)
// ---------------------------------------------------------------------------

/// High-level MD simulation driver.
///
/// Combines atoms, force field, integrator, thermostat, barostat,
/// and periodic box into a single simulation object.
pub struct MdSimulation<F, I, T, B>
where
    F: ForceField,
    I: Integrator,
    T: Thermostat,
    B: Barostat,
{
    /// The atom set being simulated.
    pub atoms: AtomSet,
    /// The force field.
    pub force_field: F,
    /// The time integrator.
    pub integrator: I,
    /// The thermostat.
    pub thermostat: T,
    /// The barostat.
    pub barostat: B,
    /// Periodic simulation box.
    pub pbox: PeriodicBox,
    /// Target temperature (for thermostat).
    pub target_temp: f64,
    /// Target pressure (for barostat).
    pub target_pressure: f64,
    /// Boltzmann constant in the chosen unit system.
    pub boltzmann_k: f64,
    /// Current step number.
    pub step_count: u64,
    /// Latest energy record.
    pub energy: EnergyRecord,
}

impl<F, I, T, B> MdSimulation<F, I, T, B>
where
    F: ForceField,
    I: Integrator,
    T: Thermostat,
    B: Barostat,
{
    /// Create a new MD simulation.
    pub fn new(
        atoms: AtomSet,
        force_field: F,
        integrator: I,
        thermostat: T,
        barostat: B,
        pbox: PeriodicBox,
        target_temp: f64,
        target_pressure: f64,
        boltzmann_k: f64,
    ) -> Self {
        Self {
            atoms,
            force_field,
            integrator,
            thermostat,
            barostat,
            pbox,
            target_temp,
            target_pressure,
            boltzmann_k,
            step_count: 0,
            energy: EnergyRecord::default(),
        }
    }

    /// Perform a single MD time step.
    pub fn step(&mut self, dt: f64) {
        let ff = &self.force_field;
        let pbox = &self.pbox;
        let mut potential_energy = 0.0;

        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.clear_forces();
            potential_energy = ff.compute_forces(atoms, pbox);
        };

        self.integrator.step(&mut self.atoms, dt, &mut force_fn);

        self.thermostat
            .apply(&mut self.atoms, self.target_temp, dt, self.boltzmann_k);

        self.barostat
            .apply(&mut self.atoms, &mut self.pbox, self.target_pressure, dt);

        for pos in &mut self.atoms.positions {
            *pos = self.pbox.wrap_position(pos);
        }

        let ke = self.atoms.kinetic_energy();
        self.energy = EnergyRecord {
            kinetic: ke,
            potential: potential_energy,
            total: ke + potential_energy,
        };

        self.step_count += 1;
    }

    /// Run the simulation for `n_steps` with time step `dt`.
    ///
    /// Optionally collects energy records at the given interval.
    /// Returns the collected energy records.
    pub fn run(
        &mut self,
        n_steps: u64,
        dt: f64,
        record_interval: Option<u64>,
    ) -> Vec<EnergyRecord> {
        let interval = record_interval.unwrap_or(0);
        let mut records = Vec::new();

        self.atoms.clear_forces();
        let _pe = self.force_field.compute_forces(&mut self.atoms, &self.pbox);

        for i in 0..n_steps {
            self.step(dt);
            if interval > 0 && (i + 1) % interval == 0 {
                records.push(self.energy);
            }
        }

        records
    }
}
