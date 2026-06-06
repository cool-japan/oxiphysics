// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Analysis and restart utilities for MD simulation.
//!
//! Contains EquilibrationDetector, RestartData, RunStatistics,
//! MdSim analysis methods, and velocity initialization functions.

use super::config::{Ensemble, KB_REDUCED, MdState};
use super::core_sim::MdSim;
use super::plain_sim::{ObservableRecord, compute_msd};

// ---------------------------------------------------------------------------
// EquilibrationDetector
// ---------------------------------------------------------------------------

/// Detects equilibration by monitoring variance of total energy over a window.
///
/// When the relative variance of energy falls below a threshold for a
/// sufficient number of consecutive windows, the system is considered
/// equilibrated.
#[derive(Debug, Clone)]
pub struct EquilibrationDetector {
    /// Window size (number of energy samples per window).
    pub window: usize,
    /// Relative variance threshold below which the system is equilibrated.
    pub threshold: f64,
    /// Number of consecutive windows below threshold before declaring equilibration.
    pub required_windows: usize,
    energy_buffer: Vec<f64>,
    consecutive_ok: usize,
    /// Whether the system has been detected as equilibrated.
    pub equilibrated: bool,
    /// The simulation step at which equilibration was first detected, if any.
    pub equilibration_step: Option<u64>,
}

impl EquilibrationDetector {
    /// Create a new equilibration detector.
    ///
    /// # Arguments
    /// * `window`           - samples per evaluation window.
    /// * `threshold`        - relative variance (Var(E) / `E`²) below which the
    ///   system is considered equilibrated.
    /// * `required_windows` - consecutive windows needed to declare equilibration.
    pub fn new(window: usize, threshold: f64, required_windows: usize) -> Self {
        assert!(window > 1, "window must be > 1");
        assert!(threshold > 0.0, "threshold must be positive");
        assert!(required_windows >= 1, "required_windows must be >= 1");
        Self {
            window,
            threshold,
            required_windows,
            energy_buffer: Vec::with_capacity(window),
            consecutive_ok: 0,
            equilibrated: false,
            equilibration_step: None,
        }
    }

    /// Feed a new total-energy sample.  Returns `true` if the system just
    /// became equilibrated with this sample.
    pub fn update(&mut self, total_energy: f64, step: u64) -> bool {
        if self.equilibrated {
            return false;
        }
        self.energy_buffer.push(total_energy);
        if self.energy_buffer.len() < self.window {
            return false;
        }
        // Evaluate variance over the window
        let n = self.energy_buffer.len();
        let mean = self.energy_buffer.iter().sum::<f64>() / n as f64;
        let var = self
            .energy_buffer
            .iter()
            .map(|e| (e - mean) * (e - mean))
            .sum::<f64>()
            / (n as f64 - 1.0);
        let rel_var = if mean.abs() > 1e-30 {
            var / (mean * mean)
        } else {
            var
        };

        // Drain window for the next evaluation
        self.energy_buffer.clear();

        if rel_var < self.threshold {
            self.consecutive_ok += 1;
        } else {
            self.consecutive_ok = 0;
        }

        if self.consecutive_ok >= self.required_windows {
            self.equilibrated = true;
            self.equilibration_step = Some(step);
            return true;
        }
        false
    }

    /// Reset the detector state.
    pub fn reset(&mut self) {
        self.energy_buffer.clear();
        self.consecutive_ok = 0;
        self.equilibrated = false;
        self.equilibration_step = None;
    }
}

// ---------------------------------------------------------------------------
// RestartData — serialization to/from Vec<f64>
// ---------------------------------------------------------------------------

/// Restart data for serializing an MD simulation state to a flat buffer.
///
/// Layout:
/// ```text
/// [n_atoms (as f64), step, time, ke, pe,
///  positions (3*n), velocities (3*n), forces (3*n), masses (n)]
/// ```
#[derive(Debug, Clone)]
pub struct RestartData {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Step at which the restart was saved.
    pub step: u64,
    /// Simulation time (ps).
    pub time: f64,
    /// Kinetic energy (kJ mol⁻¹).
    pub kinetic_energy: f64,
    /// Potential energy (kJ mol⁻¹).
    pub potential_energy: f64,
    /// Atom positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities (Å ps⁻¹).
    pub velocities: Vec<[f64; 3]>,
    /// Atom forces (kJ mol⁻¹ Å⁻¹).
    pub forces: Vec<[f64; 3]>,
    /// Atom masses (amu).
    pub masses: Vec<f64>,
}

impl RestartData {
    /// Serialise the restart data into a flat `Vec`f64` buffer.
    ///
    /// The buffer can later be passed back to [`RestartData::from_buffer`].
    pub fn to_buffer(&self) -> Vec<f64> {
        let n = self.n_atoms;
        let mut buf = Vec::with_capacity(5 + 10 * n);
        buf.push(n as f64);
        buf.push(self.step as f64);
        buf.push(self.time);
        buf.push(self.kinetic_energy);
        buf.push(self.potential_energy);
        for pos in &self.positions {
            buf.push(pos[0]);
            buf.push(pos[1]);
            buf.push(pos[2]);
        }
        for vel in &self.velocities {
            buf.push(vel[0]);
            buf.push(vel[1]);
            buf.push(vel[2]);
        }
        for frc in &self.forces {
            buf.push(frc[0]);
            buf.push(frc[1]);
            buf.push(frc[2]);
        }
        for &m in &self.masses {
            buf.push(m);
        }
        buf
    }

    /// Deserialise a restart buffer produced by `to_buffer`.
    ///
    /// Returns `None` if the buffer is malformed or too short.
    pub fn from_buffer(buf: &[f64]) -> Option<Self> {
        if buf.len() < 5 {
            return None;
        }
        let n = buf[0] as usize;
        let step = buf[1] as u64;
        let time = buf[2];
        let ke = buf[3];
        let pe = buf[4];
        let expected = 5 + 10 * n;
        if buf.len() < expected {
            return None;
        }
        let mut idx = 5usize;
        let mut positions = Vec::with_capacity(n);
        for _ in 0..n {
            positions.push([buf[idx], buf[idx + 1], buf[idx + 2]]);
            idx += 3;
        }
        let mut velocities = Vec::with_capacity(n);
        for _ in 0..n {
            velocities.push([buf[idx], buf[idx + 1], buf[idx + 2]]);
            idx += 3;
        }
        let mut forces = Vec::with_capacity(n);
        for _ in 0..n {
            forces.push([buf[idx], buf[idx + 1], buf[idx + 2]]);
            idx += 3;
        }
        let masses: Vec<f64> = buf[idx..idx + n].to_vec();
        Some(Self {
            n_atoms: n,
            step,
            time,
            kinetic_energy: ke,
            potential_energy: pe,
            positions,
            velocities,
            forces,
            masses,
        })
    }

    /// Build a [`RestartData`] from an [`MdState`].
    pub fn from_state(state: &MdState) -> Self {
        let n = state.n_atoms();
        Self {
            n_atoms: n,
            step: state.step,
            time: state.time,
            kinetic_energy: state.kinetic_energy,
            potential_energy: state.potential_energy,
            positions: state.positions.clone(),
            velocities: state.velocities.clone(),
            forces: state.forces.clone(),
            masses: state.masses.clone(),
        }
    }

    /// Restore an [`MdState`] from restart data.
    pub fn to_state(&self) -> MdState {
        let mut s = MdState::new(
            self.positions.clone(),
            self.velocities.clone(),
            self.masses.clone(),
        );
        s.forces = self.forces.clone();
        s.step = self.step;
        s.time = self.time;
        s.kinetic_energy = self.kinetic_energy;
        s.potential_energy = self.potential_energy;
        s
    }
}

// ---------------------------------------------------------------------------
// RunStatistics — comprehensive run stats
// ---------------------------------------------------------------------------

/// Comprehensive statistics collected during an MD run.
#[derive(Debug, Clone, Default)]
pub struct RunStatistics {
    /// Total number of steps completed.
    pub total_steps: u64,
    /// Total simulated time (ps).
    pub total_time: f64,
    /// Wall-clock elapsed time (seconds).
    pub wall_time_s: f64,
    /// Throughput: steps per second.
    pub steps_per_second: f64,
    /// Mean temperature (K) over the run.
    pub mean_temperature: f64,
    /// Standard deviation of temperature (K).
    pub std_temperature: f64,
    /// Mean total energy (kJ mol⁻¹).
    pub mean_total_energy: f64,
    /// Maximum energy drift |E(t) - E(0)| (kJ mol⁻¹).
    pub max_energy_drift: f64,
    /// Relative energy drift |ΔE| / |E₀|.
    pub relative_energy_drift: f64,
    /// Mean pressure (bar) over the run.
    pub mean_pressure: f64,
}

impl RunStatistics {
    /// Build from an [`ObservableRecord`] and timing info.
    pub fn from_record(record: &ObservableRecord, wall_time_s: f64) -> Self {
        let n = record.len();
        let total_steps = n as u64;
        let total_time = record.time.last().copied().unwrap_or(0.0);

        let mean_t = record.mean_temperature();
        let std_t = {
            if n < 2 {
                0.0
            } else {
                let var = record
                    .temperature
                    .iter()
                    .map(|t| (t - mean_t) * (t - mean_t))
                    .sum::<f64>()
                    / (n as f64 - 1.0);
                var.sqrt()
            }
        };

        let mean_e = record.mean_total_energy();
        let max_drift = record.max_energy_drift();
        let rel_drift = record.relative_energy_drift();
        let mean_p = if record.pressure.is_empty() {
            0.0
        } else {
            record.pressure.iter().sum::<f64>() / record.pressure.len() as f64
        };

        RunStatistics {
            total_steps,
            total_time,
            wall_time_s,
            steps_per_second: if wall_time_s > 1e-12 {
                total_steps as f64 / wall_time_s
            } else {
                f64::INFINITY
            },
            mean_temperature: mean_t,
            std_temperature: std_t,
            mean_total_energy: mean_e,
            max_energy_drift: max_drift,
            relative_energy_drift: rel_drift,
            mean_pressure: mean_p,
        }
    }

    /// Returns `true` if the energy drift is below `tol` (relative).
    pub fn energy_conserved(&self, tol: f64) -> bool {
        self.relative_energy_drift <= tol
    }
}

// ---------------------------------------------------------------------------
// MdSim extensions: MSD tracking, observable recording, full run driver
// ---------------------------------------------------------------------------

impl MdSim {
    /// Compute instantaneous MSD relative to reference positions.
    pub fn msd_from_reference(&self, reference: &[[f64; 3]]) -> f64 {
        compute_msd(&self.state.positions, reference)
    }

    /// Run a full MD simulation and return [`ObservableRecord`] and [`RunStatistics`].
    ///
    /// Combines velocity-Verlet integration, optional thermostat, observable
    /// recording, and wall-clock profiling into one convenience method.
    ///
    /// # Arguments
    /// * `forces_fn`      - external force function `(positions, box) -> forces`.
    /// * `record_every`   - record observables every this many steps (0 = never).
    /// * `virial_fn`      - optional function to compute the scalar virial for pressure.
    pub fn run_full(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        record_every: u64,
        virial_fn: Option<impl Fn(&[[f64; 3]], &[f64; 3]) -> f64>,
    ) -> (ObservableRecord, RunStatistics) {
        use std::time::Instant;

        let ref_positions = self.state.positions.clone();
        let t_start = Instant::now();

        // Initial force evaluation
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let n_steps = self.config.n_steps;
        let mut record = ObservableRecord::new();

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
                let virial = if let Some(ref vfn) = virial_fn {
                    vfn(&self.state.positions, &self.config.box_lengths)
                } else {
                    0.0
                };
                let p = self.compute_pressure(virial);
                let msd = compute_msd(&self.state.positions, &ref_positions);
                record.push(
                    self.state.time,
                    ke,
                    self.state.potential_energy,
                    self.state.temperature,
                    p,
                    msd,
                );
            }
        }

        let wall = t_start.elapsed().as_secs_f64();
        let stats = RunStatistics::from_record(&record, wall);
        (record, stats)
    }

    /// Detect equilibration using an [`EquilibrationDetector`].
    ///
    /// Runs for at most `max_steps` steps or until equilibration is detected,
    /// returning the number of steps taken.
    pub fn run_until_equilibrated(
        &mut self,
        forces_fn: impl Fn(&[[f64; 3]], &[f64; 3]) -> Vec<[f64; 3]>,
        detector: &mut EquilibrationDetector,
        max_steps: u64,
    ) -> u64 {
        self.state.forces = forces_fn(&self.state.positions, &self.config.box_lengths);

        let mut steps_taken = 0u64;
        for _ in 0..max_steps {
            self.velocity_verlet_step(&forces_fn);

            if self.config.ensemble == Ensemble::NVT || self.config.ensemble == Ensemble::NPT {
                self.apply_velocity_rescaling(self.config.temperature);
            }

            let ke = self.compute_kinetic_energy();
            self.state.kinetic_energy = ke;
            self.state.temperature = self.compute_temperature();
            self.state.step += 1;
            self.state.time += self.config.dt;
            steps_taken += 1;

            let total_e = ke + self.state.potential_energy;
            if detector.update(total_e, self.state.step) {
                break;
            }
        }
        steps_taken
    }

    /// Save the current state as a [`RestartData`] serialized to a `Vec`f64`.
    pub fn save_restart(&self) -> Vec<f64> {
        RestartData::from_state(&self.state).to_buffer()
    }

    /// Restore simulation state from a restart buffer produced by `save_restart`.
    ///
    /// Returns `true` on success, `false` if the buffer was malformed.
    pub fn load_restart(&mut self, buf: &[f64]) -> bool {
        if let Some(rd) = RestartData::from_buffer(buf) {
            self.state = rd.to_state();
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// Velocity initialization utilities
// ---------------------------------------------------------------------------

/// Assign Maxwell-Boltzmann velocities to atoms at temperature `T` (K).
///
/// Draws velocities from a Gaussian distribution with
/// σ_α = sqrt(k_B T / m) per degree of freedom.
///
/// Removes the centre-of-mass drift after assignment.
/// Uses the Box-Muller transform for Gaussian sampling.
///
/// # Arguments
/// * `masses`      - atom masses (amu).
/// * `temperature` - target temperature (K).
/// * `seed_offset` - deterministic offset to initialise the PRNG stream.
pub fn maxwell_boltzmann_velocities(
    masses: &[f64],
    temperature: f64,
    seed_offset: u64,
) -> Vec<[f64; 3]> {
    let n = masses.len();
    let mut velocities = vec![[0.0f64; 3]; n];

    // Simple deterministic LCG PRNG (no external crate needed).
    let mut rng_state = 6_364_136_223_846_793_005u64
        .wrapping_add(seed_offset)
        .wrapping_add(1_442_695_040_888_963_407);

    let mut next_rand = move || -> f64 {
        // LCG step
        rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // Convert to uniform [0, 1)
        (rng_state >> 11) as f64 / (1u64 << 53) as f64
    };

    // Box-Muller transform
    let mut next_gauss = move || -> f64 {
        loop {
            let u1 = next_rand();
            let u2 = next_rand();
            if u1 > 1e-20 {
                return (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            }
        }
    };

    for i in 0..n {
        let sigma = (KB_REDUCED * temperature / masses[i]).sqrt();
        for v in &mut velocities[i] {
            *v = sigma * next_gauss();
        }
    }

    // Remove COM velocity
    let mut total_mass = 0.0;
    let mut com_v = [0.0; 3];
    for i in 0..n {
        total_mass += masses[i];
        for a in 0..3 {
            com_v[a] += masses[i] * velocities[i][a];
        }
    }
    if total_mass > 1e-20 {
        for v in &mut com_v {
            *v /= total_mass;
        }
        for vel in velocities.iter_mut().take(n) {
            for (v, &cv) in vel.iter_mut().zip(com_v.iter()) {
                *v -= cv;
            }
        }
    }

    velocities
}

/// Rescale velocities to match a target temperature exactly.
///
/// Computes the current temperature from the given velocities and masses,
/// then scales all velocities by √(T_target / T_current).
pub fn rescale_to_temperature(velocities: &mut [[f64; 3]], masses: &[f64], target_t: f64) {
    let n = velocities.len();
    if n == 0 {
        return;
    }
    let ke: f64 = velocities
        .iter()
        .zip(masses.iter())
        .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum();
    let t_current = 2.0 * ke / (3.0 * n as f64 * KB_REDUCED);
    if t_current < 1e-10 {
        return;
    }
    let scale = (target_t / t_current).sqrt();
    for v in velocities.iter_mut() {
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;
    }
}
