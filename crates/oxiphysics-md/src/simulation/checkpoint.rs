// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Checkpoint, profiling, and monitoring types for MD simulation.
//!
//! Contains SimulationCheckpoint, SimulationProfile, NveMonitor, NptBarostat,
//! and additional MdSim methods for NPT and diagnostics.

use super::config::{Ensemble, MdConfig, MdState};
use super::core_sim::MdSim;

// ---------------------------------------------------------------------------
// SimulationCheckpoint
// ---------------------------------------------------------------------------

/// Snapshot of an MD simulation state for checkpointing / restart.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SimulationCheckpoint {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle forces at checkpoint.
    pub forces: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Step number at checkpoint.
    pub step: u64,
    /// Simulation time at checkpoint.
    pub time: f64,
    /// Kinetic energy at checkpoint.
    pub kinetic_energy: f64,
    /// Potential energy at checkpoint.
    pub potential_energy: f64,
}

// ---------------------------------------------------------------------------
// SimulationProfile — profiling support
// ---------------------------------------------------------------------------

/// Performance profile recorded during a simulation run.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SimulationProfile {
    /// Total number of integration steps taken.
    pub total_steps: u64,
    /// Wall-clock time in seconds.
    pub wall_time_s: f64,
    /// Steps per second throughput.
    pub steps_per_second: f64,
    /// Average kinetic energy over the run.
    pub avg_kinetic_energy: f64,
    /// Average temperature over the run.
    pub avg_temperature: f64,
}

// ---------------------------------------------------------------------------
// NveMonitor — NVE energy monitoring
// ---------------------------------------------------------------------------

/// Monitoring utility for NVE (micro-canonical) simulations.
///
/// Tracks total energy drift and momentum conservation.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct NveMonitor {
    /// Initial total energy (KE + PE).
    pub e0: f64,
    /// Maximum energy drift observed.
    pub max_drift: f64,
    /// Initial total momentum magnitude.
    pub p0: f64,
}

#[allow(dead_code)]
impl NveMonitor {
    /// Initialise the monitor from a simulation state.
    pub fn init(&mut self, sim: &MdSim) {
        self.e0 = sim.compute_kinetic_energy() + sim.state.potential_energy;
        let p = sim.total_momentum();
        self.p0 = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        self.max_drift = 0.0;
    }

    /// Update the monitor after a step.
    pub fn update(&mut self, sim: &MdSim) {
        let e_now = sim.compute_kinetic_energy() + sim.state.potential_energy;
        let drift = (e_now - self.e0).abs();
        if drift > self.max_drift {
            self.max_drift = drift;
        }
    }

    /// Relative energy drift: |E(t) - E(0)| / |E(0)|.
    pub fn relative_drift(&self) -> f64 {
        if self.e0.abs() < 1e-30 {
            return 0.0;
        }
        self.max_drift / self.e0.abs()
    }
}

// ---------------------------------------------------------------------------
// NptBarostat — simple Berendsen pressure control for MdSim
// ---------------------------------------------------------------------------

/// Simple isotropic Berendsen barostat for NPT simulations.
///
/// Rescales the box and positions by a factor μ = 1 - β·dt/τ·(P_target - P).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NptBarostat {
    /// Target pressure (bar).
    pub target_pressure: f64,
    /// Coupling time constant (ps).
    pub tau_p: f64,
    /// Isothermal compressibility (bar⁻¹).
    pub compressibility: f64,
}

#[allow(dead_code)]
impl NptBarostat {
    /// Create a new NPT Berendsen barostat.
    pub fn new(target_pressure: f64, tau_p: f64, compressibility: f64) -> Self {
        Self {
            target_pressure,
            tau_p,
            compressibility,
        }
    }

    /// Compute the linear scaling factor μ for the box.
    pub fn scale_factor(&self, p_current: f64, dt: f64) -> f64 {
        let mu3 = 1.0 - self.compressibility * dt / self.tau_p * (self.target_pressure - p_current);
        mu3.cbrt()
    }

    /// Apply pressure coupling to the simulation: rescale box and positions.
    pub fn apply(&self, sim: &mut MdSim, virial: f64, dt: f64) {
        let p_current = sim.compute_pressure(virial);
        let mu = self.scale_factor(p_current, dt);
        // Rescale positions and box lengths.
        for pos in sim.state.positions.iter_mut() {
            pos[0] *= mu;
            pos[1] *= mu;
            pos[2] *= mu;
        }
        for a in 0..3 {
            sim.config.box_lengths[a] *= mu;
        }
    }
}

// ---------------------------------------------------------------------------
// MdSim extensions: checkpointing, profiling, NPT
// ---------------------------------------------------------------------------

impl MdSim {
    /// Save a checkpoint of the current simulation state.
    #[allow(dead_code)]
    pub fn checkpoint(&self) -> SimulationCheckpoint {
        SimulationCheckpoint {
            positions: self.state.positions.clone(),
            velocities: self.state.velocities.clone(),
            forces: self.state.forces.clone(),
            masses: self.state.masses.clone(),
            step: self.state.step,
            time: self.state.time,
            kinetic_energy: self.state.kinetic_energy,
            potential_energy: self.state.potential_energy,
        }
    }

    /// Restore a simulation from a checkpoint.
    #[allow(dead_code)]
    pub fn restore_checkpoint(config: MdConfig, ckpt: SimulationCheckpoint) -> Self {
        let mut state = MdState::new(ckpt.positions, ckpt.velocities, ckpt.masses);
        state.forces = ckpt.forces;
        state.step = ckpt.step;
        state.time = ckpt.time;
        state.kinetic_energy = ckpt.kinetic_energy;
        state.potential_energy = ckpt.potential_energy;
        Self { config, state }
    }

    /// Run for `n_steps` with forces and collect a performance profile.
    ///
    /// Records the number of steps, wall-clock time, and average energy/temperature.
    #[allow(dead_code)]
    pub fn run_with_profiling(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        record_every: u64,
    ) -> SimulationProfile {
        use std::time::Instant;

        let t_start = Instant::now();
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let n_steps = self.config.n_steps;
        let mut sum_ke = 0.0f64;
        let mut sum_t = 0.0f64;
        let mut n_records = 0u64;

        for s in 0..n_steps {
            self.velocity_verlet_step(&forces_fn);

            if self.config.ensemble == Ensemble::NVT || self.config.ensemble == Ensemble::NPT {
                self.apply_velocity_rescaling(self.config.temperature);
            }

            let ke = self.compute_kinetic_energy();
            self.state.kinetic_energy = ke;
            self.state.temperature = self.compute_temperature();
            self.state.step += 1;
            self.state.time += self.config.dt;

            if record_every > 0 && (s + 1) % record_every == 0 {
                sum_ke += ke;
                sum_t += self.state.temperature;
                n_records += 1;
            }
        }

        let wall = t_start.elapsed().as_secs_f64();
        let (avg_ke, avg_t) = if n_records > 0 {
            (sum_ke / n_records as f64, sum_t / n_records as f64)
        } else {
            (self.state.kinetic_energy, self.state.temperature)
        };

        SimulationProfile {
            total_steps: n_steps,
            wall_time_s: wall,
            steps_per_second: if wall > 1e-12 {
                n_steps as f64 / wall
            } else {
                f64::INFINITY
            },
            avg_kinetic_energy: avg_ke,
            avg_temperature: avg_t,
        }
    }

    /// Run an NPT simulation using the Berendsen barostat and thermostat.
    ///
    /// `forces_fn` should also return a scalar virial (sum r·F).
    /// For simplicity this version uses a zero virial (ideal gas).
    #[allow(dead_code)]
    pub fn run_npt(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        barostat: &NptBarostat,
        tau_t: f64,
        record_every: u64,
    ) -> Vec<(f64, f64, f64)> {
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let mut records = Vec::new();
        let n_steps = self.config.n_steps;

        for s in 0..n_steps {
            self.velocity_verlet_step(&forces_fn);
            self.apply_berendsen_thermostat(self.config.temperature, tau_t);
            barostat.apply(self, 0.0, self.config.dt); // zero virial for simplicity

            let ke = self.compute_kinetic_energy();
            self.state.kinetic_energy = ke;
            self.state.temperature = self.compute_temperature();
            self.state.step += 1;
            self.state.time += self.config.dt;

            if record_every > 0 && (s + 1) % record_every == 0 {
                records.push((ke, self.state.potential_energy, self.state.temperature));
            }
        }
        records
    }
}
