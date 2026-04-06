#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Plain (non-generic) MD simulation and observable tracking.
//!
//! Contains MdSimulationPlain, ObservableRecord, and displacement functions.

use super::config::KB_REDUCED;
use super::forces::{Atom, BarostatType, ThermostatType};

// ---------------------------------------------------------------------------
// MdSimulationPlain — self-contained simulation on Vec<Atom>
// ---------------------------------------------------------------------------

/// Self-contained MD simulation driver operating on a `Vec`Atom`.
///
/// Provides velocity-Verlet integration, PBC wrapping, temperature and
/// kinetic energy calculation, and virial pressure.  Optional `thermostat`
/// and `barostat` fields select ensemble control without external trait objects.
#[allow(dead_code)]
pub struct MdSimulationPlain {
    /// All atoms in the simulation.
    pub atoms: Vec<Atom>,
    /// Orthorhombic box lengths [Lx, Ly, Lz] (Å).
    pub box_lengths: [f64; 3],
    /// Number of completed integration steps.
    pub step_count: u64,
    /// Simulated time (ps).
    pub time: f64,
    /// Optional thermostat for NVT / NPT.
    pub thermostat: Option<ThermostatType>,
    /// Optional barostat for NPT.
    pub barostat: Option<BarostatType>,
    /// Berendsen thermostat coupling time (ps).
    pub tau_t: f64,
    /// Target temperature for thermostat (K).
    pub target_temp: f64,
}

#[allow(dead_code)]
impl MdSimulationPlain {
    /// Create a new simulation in the NVE ensemble.
    pub fn new(atoms: Vec<Atom>, box_lengths: [f64; 3]) -> Self {
        Self {
            atoms,
            box_lengths,
            step_count: 0,
            time: 0.0,
            thermostat: None,
            barostat: None,
            tau_t: 0.1,
            target_temp: 300.0,
        }
    }

    /// Number of atoms.
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    // -----------------------------------------------------------------------
    // PBC wrapping
    // -----------------------------------------------------------------------

    /// Wrap all atom positions into the primary box [0, L) using PBC.
    pub fn apply_pbc(&mut self) {
        for atom in &mut self.atoms {
            for a in 0..3 {
                let l = self.box_lengths[a];
                atom.position[a] = atom.position[a].rem_euclid(l);
            }
        }
    }

    /// Minimum image displacement dr → [-L/2, L/2) for an orthorhombic box.
    #[inline]
    fn minimum_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let mut r = dr;
        for a in 0..3 {
            let l = self.box_lengths[a];
            r[a] -= l * (r[a] / l).round();
        }
        r
    }

    // -----------------------------------------------------------------------
    // Energy and temperature
    // -----------------------------------------------------------------------

    /// Kinetic energy (kJ mol⁻¹): KE = ½ Σ mᵢ |vᵢ|².
    pub fn kinetic_energy(&self) -> f64 {
        self.atoms.iter().map(|a| a.kinetic_energy()).sum()
    }

    /// Instantaneous temperature (K): T = 2 KE / (3 N k_B).
    pub fn temperature(&self) -> f64 {
        let n = self.n_atoms() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (3.0 * n * KB_REDUCED)
    }

    /// Instantaneous pressure (bar) using the virial theorem.
    ///
    /// P = (N k_B T + W / 3) / V, converted to bar.
    ///
    /// `virial` is the scalar W = Σ rᵢ · Fᵢ (kJ mol⁻¹).
    pub fn pressure_virial(&self, virial: f64) -> f64 {
        let n = self.n_atoms() as f64;
        let v = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        if v < 1e-30 {
            return 0.0;
        }
        let t = self.temperature();
        let ideal = n * KB_REDUCED * t;
        let p_kj_ang3 = (ideal + virial / 3.0) / v;
        // 1 kJ mol⁻¹ Å⁻³ = 16 605.4 bar
        p_kj_ang3 * 16_605.4
    }

    // -----------------------------------------------------------------------
    // Velocity-Verlet integration
    // -----------------------------------------------------------------------

    /// Perform one velocity-Verlet step with a zero-force placeholder.
    ///
    /// All force components are assumed to be pre-stored in `atom.force`.
    /// After the step, forces are left unchanged (caller must update them
    /// before the next step for a real physical simulation).
    pub fn step(&mut self, dt: f64) {
        let n = self.n_atoms();

        // Half-step velocity: v += F/(2m) * dt
        for i in 0..n {
            let inv_m = 1.0 / self.atoms[i].mass;
            for a in 0..3 {
                self.atoms[i].velocity[a] += 0.5 * dt * self.atoms[i].force[a] * inv_m;
            }
        }

        // Full-step position: x += v * dt
        for i in 0..n {
            for a in 0..3 {
                self.atoms[i].position[a] += dt * self.atoms[i].velocity[a];
            }
        }

        // PBC wrap
        self.apply_pbc();

        // NOTE: a real simulation would compute new forces here.
        // For free-particle / demo integration the forces stay zero.

        // Second half-step velocity: v += F/(2m) * dt
        for i in 0..n {
            let inv_m = 1.0 / self.atoms[i].mass;
            for a in 0..3 {
                self.atoms[i].velocity[a] += 0.5 * dt * self.atoms[i].force[a] * inv_m;
            }
        }

        // Ensemble control
        match self.thermostat {
            Some(ThermostatType::Berendsen) | Some(ThermostatType::Csvr) => {
                self.apply_berendsen_rescale(self.target_temp, dt);
            }
            Some(ThermostatType::NoseHoover) => {
                // Simplified: rescale (full Nosé-Hoover would require chain variables)
                self.apply_berendsen_rescale(self.target_temp, dt);
            }
            None => {}
        }

        self.step_count += 1;
        self.time += dt;
    }

    /// Apply Berendsen velocity rescaling toward `target_t` (K) with coupling `tau_t`.
    fn apply_berendsen_rescale(&mut self, target_t: f64, dt: f64) {
        let t_current = self.temperature();
        if t_current < 1e-10 {
            return;
        }
        let ratio = target_t / t_current;
        let scale = (1.0 + dt / self.tau_t * (ratio - 1.0)).max(0.0).sqrt();
        for atom in &mut self.atoms {
            for a in 0..3 {
                atom.velocity[a] *= scale;
            }
        }
    }

    /// Run for `n_steps` steps of size `dt`, returning (KE, T) at each step.
    pub fn run(&mut self, n_steps: u64, dt: f64) -> Vec<(f64, f64)> {
        let mut records = Vec::with_capacity(n_steps as usize);
        for _ in 0..n_steps {
            self.step(dt);
            records.push((self.kinetic_energy(), self.temperature()));
        }
        records
    }

    /// Compute the virial scalar W = Σᵢ Σⱼ>ᵢ (r_ij · F_ij) using stored forces.
    ///
    /// Uses the scalar contribution of position dotted with force per atom:
    /// W ≈ Σᵢ rᵢ · Fᵢ (approximate; exact virial requires pair decomposition).
    pub fn virial_approx(&self) -> f64 {
        self.atoms
            .iter()
            .map(|a| {
                a.position[0] * a.force[0] + a.position[1] * a.force[1] + a.position[2] * a.force[2]
            })
            .sum()
    }

    /// Remove centre-of-mass velocity.
    pub fn remove_com_velocity(&mut self) {
        let n = self.n_atoms();
        if n == 0 {
            return;
        }
        let mut total_mass = 0.0_f64;
        let mut com_v = [0.0_f64; 3];
        for atom in &self.atoms {
            total_mass += atom.mass;
            for a in 0..3 {
                com_v[a] += atom.mass * atom.velocity[a];
            }
        }
        if total_mass > 1e-20 {
            for a in 0..3 {
                com_v[a] /= total_mass;
            }
            for atom in &mut self.atoms {
                for a in 0..3 {
                    atom.velocity[a] -= com_v[a];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ObservableRecord — observable history for analysis
// ---------------------------------------------------------------------------

/// Time-series record of MD observables collected during a run.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ObservableRecord {
    /// Simulation time (ps) at each recorded frame.
    pub time: Vec<f64>,
    /// Kinetic energy (kJ mol⁻¹) per frame.
    pub kinetic_energy: Vec<f64>,
    /// Potential energy (kJ mol⁻¹) per frame.
    pub potential_energy: Vec<f64>,
    /// Total energy (kJ mol⁻¹) per frame.
    pub total_energy: Vec<f64>,
    /// Instantaneous temperature (K) per frame.
    pub temperature: Vec<f64>,
    /// Instantaneous pressure (bar) per frame.
    pub pressure: Vec<f64>,
    /// Mean squared displacement (Å²) per frame.
    pub msd: Vec<f64>,
}

#[allow(dead_code)]
impl ObservableRecord {
    /// Create a new empty record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a frame.
    #[allow(clippy::too_many_arguments)]
    pub fn push(&mut self, time: f64, ke: f64, pe: f64, temp: f64, pressure: f64, msd: f64) {
        self.time.push(time);
        self.kinetic_energy.push(ke);
        self.potential_energy.push(pe);
        self.total_energy.push(ke + pe);
        self.temperature.push(temp);
        self.pressure.push(pressure);
        self.msd.push(msd);
    }

    /// Number of recorded frames.
    pub fn len(&self) -> usize {
        self.time.len()
    }

    /// Is the record empty?
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// Mean temperature over all frames.
    pub fn mean_temperature(&self) -> f64 {
        if self.temperature.is_empty() {
            return 0.0;
        }
        self.temperature.iter().sum::<f64>() / self.temperature.len() as f64
    }

    /// Mean total energy over all frames.
    pub fn mean_total_energy(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        self.total_energy.iter().sum::<f64>() / self.total_energy.len() as f64
    }

    /// Variance of total energy over all frames.
    ///
    /// Var(E) = `E²` - `E`²
    pub fn energy_variance(&self) -> f64 {
        let n = self.total_energy.len();
        if n < 2 {
            return 0.0;
        }
        let mean = self.mean_total_energy();
        let sum_sq: f64 = self
            .total_energy
            .iter()
            .map(|e| (e - mean) * (e - mean))
            .sum();
        sum_sq / (n as f64 - 1.0)
    }

    /// Maximum absolute energy drift from frame 0.
    pub fn max_energy_drift(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        let e0 = self.total_energy[0];
        self.total_energy
            .iter()
            .map(|&e| (e - e0).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Relative energy drift: max_drift / |E₀|.
    pub fn relative_energy_drift(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        let e0 = self.total_energy[0].abs();
        if e0 < 1e-30 {
            return 0.0;
        }
        self.max_energy_drift() / e0
    }
}

// ---------------------------------------------------------------------------
// MSD — Mean Squared Displacement
// ---------------------------------------------------------------------------

/// Compute the mean squared displacement (Å²) between two sets of positions.
///
/// MSD = (1/N) Σᵢ |rᵢ(t) - rᵢ(0)|²
///
/// Note: this uses unwrapped positions (no PBC folding) for correct MSD.
#[allow(dead_code)]
pub fn compute_msd(positions_now: &[[f64; 3]], positions_ref: &[[f64; 3]]) -> f64 {
    let n = positions_now.len().min(positions_ref.len());
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        for a in 0..3 {
            let d = positions_now[i][a] - positions_ref[i][a];
            sum += d * d;
        }
    }
    sum / n as f64
}

/// Compute per-atom squared displacements (Å²).
///
/// Returns a vector of |rᵢ(t) - rᵢ(0)|² for each atom.
#[allow(dead_code)]
pub fn per_atom_displacement_sq(
    positions_now: &[[f64; 3]],
    positions_ref: &[[f64; 3]],
) -> Vec<f64> {
    let n = positions_now.len().min(positions_ref.len());
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let mut d2 = 0.0;
        for a in 0..3 {
            let d = positions_now[i][a] - positions_ref[i][a];
            d2 += d * d;
        }
        result.push(d2);
    }
    result
}
