#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core MD simulation driver (MdSim) and Berendsen thermostat.

use super::config::{Ensemble, KB_REDUCED, MdConfig, MdState};

// ---------------------------------------------------------------------------
// MdSim — self-contained MD driver
// ---------------------------------------------------------------------------

/// Self-contained MD simulation driver using plain arrays (no nalgebra).
///
/// Operates on [`MdState`] with an explicit velocity-Verlet integration loop
/// and optional velocity-rescaling thermostat.
#[allow(dead_code)]
pub struct MdSim {
    /// Run configuration.
    pub config: MdConfig,
    /// Current dynamical state.
    pub state: MdState,
}

impl MdSim {
    /// Create a new simulation.
    #[allow(dead_code)]
    pub fn new(config: MdConfig, state: MdState) -> Self {
        Self { config, state }
    }

    // -----------------------------------------------------------------------
    // Kinetic energy
    // -----------------------------------------------------------------------

    /// Compute kinetic energy (kJ mol⁻¹) from current velocities and masses.
    ///
    /// KE = ½ Σᵢ mᵢ |vᵢ|²
    #[allow(dead_code)]
    pub fn compute_kinetic_energy(&self) -> f64 {
        self.state
            .velocities
            .iter()
            .zip(self.state.masses.iter())
            .map(|(v, &m)| {
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * m * v2
            })
            .sum()
    }

    // -----------------------------------------------------------------------
    // Temperature
    // -----------------------------------------------------------------------

    /// Compute instantaneous temperature (K) from kinetic energy.
    ///
    /// T = 2 KE / (3 N k_B)
    #[allow(dead_code)]
    pub fn compute_temperature(&self) -> f64 {
        let n = self.state.n_atoms() as f64;
        if n == 0.0 {
            return 0.0;
        }
        let ke = self.compute_kinetic_energy();
        2.0 * ke / (3.0 * n * KB_REDUCED)
    }

    // -----------------------------------------------------------------------
    // Pressure
    // -----------------------------------------------------------------------

    /// Compute instantaneous pressure (bar) using ideal + virial contributions.
    ///
    /// P = (N k_B T + W) / (3 V)
    ///
    /// where W is the virial (kJ mol⁻¹) and V is the box volume (Å³).
    /// Result converted to bar: 1 kJ mol⁻¹ Å⁻³ ≈ 16.6054 bar.
    #[allow(dead_code)]
    pub fn compute_pressure(&self, virial: f64) -> f64 {
        let n = self.state.n_atoms() as f64;
        let v =
            self.config.box_lengths[0] * self.config.box_lengths[1] * self.config.box_lengths[2];
        if v == 0.0 {
            return 0.0;
        }
        let t = self.compute_temperature();
        let ideal = n * KB_REDUCED * t;
        // conversion factor: 1 kJ mol⁻¹ Å⁻³ = 16605.4 bar → divide by 3V
        let p_kj_per_ang3 = (ideal + virial) / (3.0 * v);
        p_kj_per_ang3 * 16_605.4
    }

    // -----------------------------------------------------------------------
    // Velocity-Verlet integrator (higher-order version with forces_fn)
    // -----------------------------------------------------------------------

    /// Perform one velocity-Verlet step with an external force function.
    ///
    /// `forces_fn(positions, box_lengths) -> forces`
    #[allow(dead_code)]
    pub fn velocity_verlet_step(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
    ) {
        let dt = self.config.dt;
        let n = self.state.n_atoms();

        // Half-step velocity update: v(t + dt/2) = v(t) + F(t)/(2m) * dt
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt * self.state.forces[i][a] * inv_m;
            }
        }

        // Full-step position update: x(t + dt) = x(t) + v(t + dt/2) * dt
        for i in 0..n {
            for a in 0..3 {
                self.state.positions[i][a] += dt * self.state.velocities[i][a];
            }
        }

        // PBC wrap
        if self.config.pbc {
            for i in 0..n {
                for a in 0..3 {
                    let l = self.config.box_lengths[a];
                    self.state.positions[i][a] = self.state.positions[i][a].rem_euclid(l);
                }
            }
        }

        // Compute new forces
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        // Second half-step velocity update
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt * self.state.forces[i][a] * inv_m;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Velocity-rescaling thermostat
    // -----------------------------------------------------------------------

    /// Apply simple velocity-rescaling thermostat to hit `target_T` (K).
    ///
    /// Scales all velocities by √(T_target / T_current).
    #[allow(dead_code)]
    pub fn apply_velocity_rescaling(&mut self, target_t: f64) {
        let t_current = self.compute_temperature();
        if t_current < 1e-10 {
            return;
        }
        let scale = (target_t / t_current).sqrt();
        for v in self.state.velocities.iter_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }

    // -----------------------------------------------------------------------
    // Single combined step
    // -----------------------------------------------------------------------

    /// Perform one full MD step: integrate, then apply ensemble controls.
    ///
    /// Uses a zero-force placeholder; real simulations should call
    /// `velocity_verlet_step` directly with a physical force function.
    #[allow(dead_code)]
    pub fn step(&mut self) {
        let zero_forces = |pos: &[[f64; 3]], _box: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
        self.velocity_verlet_step(zero_forces);

        match self.config.ensemble {
            Ensemble::NVT | Ensemble::NPT => {
                self.apply_velocity_rescaling(self.config.temperature);
            }
            Ensemble::NVE => {}
        }

        let ke = self.compute_kinetic_energy();
        let t = self.compute_temperature();
        self.state.kinetic_energy = ke;
        self.state.temperature = t;
        self.state.step += 1;
        self.state.time += self.config.dt;
    }

    /// Run for `n_steps` steps, recording (KE, PE, T) every `record_every` steps.
    #[allow(dead_code)]
    pub fn run_with_forces(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        record_every: u64,
    ) -> Vec<(f64, f64, f64)> {
        // Initial forces
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let mut records = Vec::new();
        let n_steps = self.config.n_steps;
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
                records.push((ke, self.state.potential_energy, self.state.temperature));
            }
        }
        records
    }
}

// ---------------------------------------------------------------------------
// Berendsen thermostat (self-contained, no trait dep)
// ---------------------------------------------------------------------------

/// Berendsen velocity-rescaling thermostat for the self-contained MD driver.
///
/// Rescales velocities toward `target_t` with a coupling constant `tau`:
///
/// ```text
/// scale = sqrt(1 + dt/tau * (T_target / T_current - 1))
/// ```
#[allow(dead_code)]
pub struct BerendsenRescaler {
    /// Coupling time constant (ps).
    pub tau: f64,
}

impl BerendsenRescaler {
    /// Create a new Berendsen rescaler.
    #[allow(dead_code)]
    pub fn new(tau: f64) -> Self {
        assert!(tau > 0.0, "tau must be positive");
        Self { tau }
    }

    /// Compute the velocity scaling factor for this time step.
    #[allow(dead_code)]
    pub fn scale_factor(&self, t_current: f64, target_t: f64, dt: f64) -> f64 {
        if t_current < 1e-10 {
            return 1.0;
        }
        let ratio = target_t / t_current;
        (1.0 + dt / self.tau * (ratio - 1.0)).max(0.0).sqrt()
    }
}

impl MdSim {
    /// Apply Berendsen velocity rescaling with coupling time `tau` (ps).
    #[allow(dead_code)]
    pub fn apply_berendsen_thermostat(&mut self, target_t: f64, tau: f64) {
        let t_current = self.compute_temperature();
        let rescaler = BerendsenRescaler::new(tau);
        let scale = rescaler.scale_factor(t_current, target_t, self.config.dt);
        for v in self.state.velocities.iter_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }

    /// Remove centre-of-mass velocity from the system.
    #[allow(dead_code)]
    pub fn remove_com_velocity(&mut self) {
        let n = self.state.n_atoms();
        if n == 0 {
            return;
        }
        let mut total_mass = 0.0;
        let mut com_vel = [0.0; 3];
        for i in 0..n {
            let m = self.state.masses[i];
            total_mass += m;
            for a in 0..3 {
                com_vel[a] += m * self.state.velocities[i][a];
            }
        }
        if total_mass > 0.0 {
            for a in 0..3 {
                com_vel[a] /= total_mass;
            }
            for i in 0..n {
                for a in 0..3 {
                    self.state.velocities[i][a] -= com_vel[a];
                }
            }
        }
    }

    /// Compute total momentum \[px, py, pz\] (amu Å ps⁻¹).
    #[allow(dead_code)]
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut p = [0.0; 3];
        for i in 0..self.state.n_atoms() {
            let m = self.state.masses[i];
            for a in 0..3 {
                p[a] += m * self.state.velocities[i][a];
            }
        }
        p
    }

    /// Compute total energy (kinetic + potential) from stored state (kJ mol⁻¹).
    #[allow(dead_code)]
    pub fn total_energy(&self) -> f64 {
        self.state.kinetic_energy + self.state.potential_energy
    }

    /// Perform a full MD step using an external force function and Berendsen
    /// thermostat (NVT) or simple velocity-Verlet (NVE).
    #[allow(dead_code)]
    pub fn step_with_forces(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        tau: f64,
    ) {
        self.velocity_verlet_step(forces_fn);
        match self.config.ensemble {
            Ensemble::NVT | Ensemble::NPT => {
                self.apply_berendsen_thermostat(self.config.temperature, tau);
            }
            Ensemble::NVE => {}
        }
        self.state.kinetic_energy = self.compute_kinetic_energy();
        self.state.temperature = self.compute_temperature();
        self.state.step += 1;
        self.state.time += self.config.dt;
    }
}
