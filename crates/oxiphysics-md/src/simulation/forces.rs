// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Force computation utilities, MD driver, and atom types.
//!
//! Contains standalone force functions, the MdDriver high-level runner,
//! thermostat/barostat type enums, and the standalone Atom type.

use super::config::{KB_REDUCED, MdState};

// ---------------------------------------------------------------------------
// Standalone PBC / pair helpers
// ---------------------------------------------------------------------------

/// Apply the minimum image convention for a single orthorhombic box.
///
/// Returns the displacement vector mapped to the principal image \[-L/2, L/2).
pub fn apply_minimum_image(r: [f64; 3], box_len: f64) -> [f64; 3] {
    [
        r[0] - box_len * (r[0] / box_len).round(),
        r[1] - box_len * (r[1] / box_len).round(),
        r[2] - box_len * (r[2] / box_len).round(),
    ]
}

/// Squared distance between two positions (no PBC).
pub fn pair_distance_sq(ri: [f64; 3], rj: [f64; 3]) -> f64 {
    let dx = rj[0] - ri[0];
    let dy = rj[1] - ri[1];
    let dz = rj[2] - ri[2];
    dx * dx + dy * dy + dz * dz
}

// ---------------------------------------------------------------------------
// compute_lj_forces — full O(N²) Lennard-Jones force loop with PBC
// ---------------------------------------------------------------------------

/// Compute Lennard-Jones forces and accumulate into `state.forces`.
///
/// Uses the 12-6 LJ potential with a spherical cutoff and minimum image
/// convention for a cubic periodic box.
///
/// Force magnitudes are in kJ mol⁻¹ Å⁻¹ (if masses are amu and distances Å).
/// epsilon: well depth (kJ mol⁻¹), sigma: finite distance (Å),
/// cutoff: pair cutoff (Å), box_len: cubic box side (Å).
pub fn compute_lj_forces(state: &mut MdState, epsilon: f64, sigma: f64, cutoff: f64, box_len: f64) {
    let n = state.n_atoms();
    // Zero forces first
    for f in state.forces.iter_mut() {
        *f = [0.0; 3];
    }

    let cutoff2 = cutoff * cutoff;
    let sigma2 = sigma * sigma;

    for i in 0..n {
        for j in (i + 1)..n {
            let dr_raw = [
                state.positions[j][0] - state.positions[i][0],
                state.positions[j][1] - state.positions[i][1],
                state.positions[j][2] - state.positions[i][2],
            ];
            let dr = apply_minimum_image(dr_raw, box_len);
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r2 >= cutoff2 || r2 < 1e-20 {
                continue;
            }

            // sr2 = (sigma/r)^2
            let sr2 = sigma2 / r2;
            let sr6 = sr2 * sr2 * sr2;
            let sr12 = sr6 * sr6;

            // F_i = (24 eps / r) * [(sigma/r)^6 - 2*(sigma/r)^12] * r_hat_ij
            // f_over_r = 24 * epsilon * (sr6 - 2*sr12) / r^2
            // (positive at r > r_min = repulsive for atom i toward j)
            let f_over_r = 24.0 * epsilon * (sr6 - 2.0 * sr12) / r2;

            let (top, bot) = state.forces.split_at_mut(j);
            for ((&dra, fi), fj) in dr.iter().zip(top[i].iter_mut()).zip(bot[0].iter_mut()) {
                let f_a = f_over_r * dra;
                *fi += f_a;
                *fj -= f_a;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MdDriver — alternative MD driver struct
// ---------------------------------------------------------------------------

/// Standalone MD driver operating on [`MdState`].
///
/// Uses velocity-Verlet integration and optional velocity-rescaling thermostat.
/// No dependencies on external crates.
pub struct MdDriver {
    /// Dynamical state.
    pub state: MdState,
    /// Time step (ps).
    pub dt: f64,
    /// Number of completed integration steps.
    pub step_count: u64,
    /// Apply thermostat every this many steps (0 = never).
    pub thermostat_freq: usize,
}

impl MdDriver {
    /// Create a new [`MdDriver`].
    pub fn new(state: MdState, dt: f64) -> Self {
        Self {
            state,
            dt,
            step_count: 0,
            thermostat_freq: 0,
        }
    }

    /// Kinetic energy (kJ mol⁻¹) of the current state.
    pub fn kinetic_energy(&self) -> f64 {
        self.state
            .velocities
            .iter()
            .zip(self.state.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Instantaneous temperature (K) via equipartition.
    ///
    /// T = 2 KE / (3 N k_B)
    pub fn temperature(&self) -> f64 {
        let n = self.state.n_atoms() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (3.0 * n * KB_REDUCED)
    }

    /// Perform one velocity-Verlet integration step.
    ///
    /// `force_fn` receives `&mut MdState` and must fill `state.forces`.
    pub fn velocity_verlet_step<F: Fn(&mut MdState)>(&mut self, force_fn: F) {
        let dt = self.dt;
        let n = self.state.n_atoms();

        // Half-step velocity: v += F/(2m) * dt
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt * self.state.forces[i][a] * inv_m;
            }
        }

        // Full-step position: x += v * dt
        for i in 0..n {
            for a in 0..3 {
                self.state.positions[i][a] += dt * self.state.velocities[i][a];
            }
        }

        // Compute new forces
        force_fn(&mut self.state);

        // Second half-step velocity: v += F/(2m) * dt
        for i in 0..n {
            let inv_m = 1.0 / self.state.masses[i];
            for a in 0..3 {
                self.state.velocities[i][a] += 0.5 * dt * self.state.forces[i][a] * inv_m;
            }
        }

        self.step_count += 1;
    }

    /// Velocity-rescaling thermostat: scale velocities to hit `t_target` (K).
    pub fn rescale_velocities_to_temperature(&mut self, t_target: f64) {
        let t_current = self.temperature();
        if t_current < 1e-10 {
            return;
        }
        let scale = (t_target / t_current).sqrt();
        for v in self.state.velocities.iter_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }

    /// Remove centre-of-mass velocity drift.
    pub fn remove_com_velocity(&mut self) {
        let n = self.state.n_atoms();
        if n == 0 {
            return;
        }
        let mut total_mass = 0.0f64;
        let mut com_v = [0.0f64; 3];
        for (&m, vel) in self
            .state
            .masses
            .iter()
            .zip(self.state.velocities.iter())
            .take(n)
        {
            total_mass += m;
            for (cv, &v) in com_v.iter_mut().zip(vel.iter()) {
                *cv += m * v;
            }
        }
        if total_mass > 1e-20 {
            for v in &mut com_v {
                *v /= total_mass;
            }
            for vel in self.state.velocities[..n].iter_mut() {
                for (v, &cv) in vel.iter_mut().zip(com_v.iter()) {
                    *v -= cv;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ThermostatType / BarostatType — plain enums for self-contained MD
// ---------------------------------------------------------------------------

/// Thermostat algorithm for NVT / NPT self-contained simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermostatType {
    /// Berendsen velocity-rescaling thermostat.
    Berendsen,
    /// Nosé-Hoover chain thermostat.
    NoseHoover,
    /// Canonical sampling through velocity rescaling (CSVR).
    Csvr,
}

/// Barostat algorithm for NPT self-contained simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarostatType {
    /// Berendsen isotropic pressure coupling.
    Berendsen,
    /// Monte-Carlo barostat (volume moves).
    MonteCarlo,
    /// Parrinello-Rahman fully flexible cell barostat.
    Parrinello,
}

// ---------------------------------------------------------------------------
// Atom — simple self-contained atom representation
// ---------------------------------------------------------------------------

/// A single atom with position, velocity, force, mass, and charge.
///
/// Uses plain `[f64;3]` arrays for all vector quantities.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Position (Å).
    pub position: [f64; 3],
    /// Velocity (Å ps⁻¹).
    pub velocity: [f64; 3],
    /// Force (kJ mol⁻¹ Å⁻¹).
    pub force: [f64; 3],
    /// Mass (amu).
    pub mass: f64,
    /// Partial charge (e).
    pub charge: f64,
}

impl Atom {
    /// Create a new atom.
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, charge: f64) -> Self {
        Self {
            position,
            velocity,
            force: [0.0; 3],
            mass,
            charge,
        }
    }

    /// Kinetic energy of this atom (kJ mol⁻¹): ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.velocity[0] * self.velocity[0]
            + self.velocity[1] * self.velocity[1]
            + self.velocity[2] * self.velocity[2];
        0.5 * self.mass * v2
    }
}
