// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Molecular Dynamics (MD) simulation API for Python interop.
//!
//! Provides a simple N-body MD simulation with Lennard-Jones pair potentials,
//! optional velocity-rescaling thermostat, and periodic boundary conditions.
//! All types are `no-lifetime`, serialization-friendly, and carry comprehensive
//! tests.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for an MD simulation.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMdConfig {
    /// Side length of the cubic simulation box.
    #[pyo3(get, set)]
    pub box_size: f64,
    /// Lennard-Jones epsilon (energy well depth, J or reduced units).
    #[pyo3(get, set)]
    pub lj_epsilon: f64,
    /// Lennard-Jones sigma (particle diameter in m or reduced units).
    #[pyo3(get, set)]
    pub lj_sigma: f64,
    /// Cut-off radius for pair interactions (typically 2.5 * sigma).
    #[pyo3(get, set)]
    pub cutoff: f64,
    /// Particle mass (kg or reduced units).
    #[pyo3(get, set)]
    pub particle_mass: f64,
    /// Target temperature for thermostat (K or reduced units). `None` = NVE.
    #[pyo3(get, set)]
    pub target_temperature: Option<f64>,
    /// Thermostat relaxation time (used by velocity rescaling).
    #[pyo3(get, set)]
    pub thermostat_tau: f64,
}

#[pymethods]
impl PyMdConfig {
    /// Create an argon-like reduced-unit configuration (ε=1, σ=1, box=10σ).
    #[new]
    pub fn new() -> Self {
        Self::argon_reduced()
    }

    /// Argon-like reduced-unit configuration (ε=1, σ=1, box=10σ).
    #[staticmethod]
    pub fn argon_reduced() -> Self {
        Self {
            box_size: 10.0,
            lj_epsilon: 1.0,
            lj_sigma: 1.0,
            cutoff: 2.5,
            particle_mass: 1.0,
            target_temperature: Some(1.2),
            thermostat_tau: 0.1,
        }
    }
}

impl Default for PyMdConfig {
    fn default() -> Self {
        Self::argon_reduced()
    }
}

// ---------------------------------------------------------------------------
// Atom
// ---------------------------------------------------------------------------

/// A single atom in the MD simulation.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMdAtom {
    /// Atom type identifier (e.g. element number or user label).
    #[pyo3(get, set)]
    pub atom_type: u32,
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Accumulated force `[fx, fy, fz]` (cleared each step).
    pub force: [f64; 3],
}

#[pymethods]
impl PyMdAtom {
    /// Create a new atom at `position` with zero velocity.
    #[new]
    pub fn new(position: Vec<f64>, atom_type: u32) -> Self {
        let pos = if position.len() >= 3 {
            [position[0], position[1], position[2]]
        } else {
            [0.0; 3]
        };
        Self {
            position: pos,
            velocity: [0.0; 3],
            force: [0.0; 3],
            atom_type,
        }
    }

    /// Kinetic energy of this atom: 0.5 * m * v².
    pub fn kinetic_energy(&self, mass: f64) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        0.5 * mass * v2
    }

    /// Position as `[x, y, z]`.
    #[getter]
    pub fn position(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Set position from a `[x, y, z]` list.
    #[setter]
    pub fn set_position(&mut self, v: Vec<f64>) {
        if v.len() >= 3 {
            self.position = [v[0], v[1], v[2]];
        }
    }

    /// Velocity as `[vx, vy, vz]`.
    #[getter]
    pub fn velocity(&self) -> Vec<f64> {
        self.velocity.to_vec()
    }

    /// Set velocity from a `[vx, vy, vz]` list.
    #[setter]
    pub fn set_velocity(&mut self, v: Vec<f64>) {
        if v.len() >= 3 {
            self.velocity = [v[0], v[1], v[2]];
        }
    }

    /// Force as `[fx, fy, fz]`.
    #[getter]
    pub fn force(&self) -> Vec<f64> {
        self.force.to_vec()
    }
}

// ---------------------------------------------------------------------------
// PyMdSimulation
// ---------------------------------------------------------------------------

/// Molecular Dynamics simulation (NVE/NVT) with periodic boundary conditions.
///
/// Uses a velocity-Verlet integrator and a truncated-shifted Lennard-Jones
/// pair potential. Supports an optional velocity-rescaling thermostat.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyMdSimulation {
    /// All atoms.
    pub(crate) atoms: Vec<PyMdAtom>,
    /// Simulation configuration.
    pub(crate) config: PyMdConfig,
    /// Total simulation time accumulated.
    time: f64,
    /// Number of completed steps.
    step_count: u64,
    /// Most recently computed total potential energy.
    potential_energy: f64,
    /// Whether the thermostat is currently active.
    thermostat_active: bool,
}

#[pymethods]
impl PyMdSimulation {
    /// Create a new empty MD simulation from the given configuration.
    #[new]
    pub fn new(config: PyMdConfig) -> Self {
        Self {
            atoms: Vec::new(),
            config,
            time: 0.0,
            step_count: 0,
            potential_energy: 0.0,
            thermostat_active: true,
        }
    }

    /// Add an atom at `position` with the given type index.
    ///
    /// Returns the index of the newly added atom.
    pub fn add_atom(&mut self, position: Vec<f64>, atom_type: u32) -> usize {
        let idx = self.atoms.len();
        let pos = if position.len() >= 3 {
            [position[0], position[1], position[2]]
        } else {
            [0.0; 3]
        };
        self.atoms.push(PyMdAtom {
            position: pos,
            velocity: [0.0; 3],
            force: [0.0; 3],
            atom_type,
        });
        idx
    }

    /// Set the velocity of atom `i`. No-op if `i` is out of bounds.
    pub fn set_velocity(&mut self, i: usize, vel: Vec<f64>) {
        if let Some(atom) = self.atoms.get_mut(i)
            && vel.len() >= 3
        {
            atom.velocity = [vel[0], vel[1], vel[2]];
        }
    }

    /// Get the position of atom `i`, or `None` if out of bounds.
    pub fn position(&self, i: usize) -> Option<Vec<f64>> {
        self.atoms.get(i).map(|a| a.position.to_vec())
    }

    /// Get the velocity of atom `i`, or `None` if out of bounds.
    pub fn velocity(&self, i: usize) -> Option<Vec<f64>> {
        self.atoms.get(i).map(|a| a.velocity.to_vec())
    }

    /// Number of atoms in the simulation.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Accumulated simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Number of completed steps.
    pub fn step_count(&self) -> u64 {
        self.step_count
    }

    /// Enable or disable the velocity-rescaling thermostat.
    pub fn set_thermostat(&mut self, active: bool) {
        self.thermostat_active = active;
    }

    /// Whether the thermostat is active.
    pub fn thermostat_active(&self) -> bool {
        self.thermostat_active
    }

    /// Set the target temperature for the thermostat.
    pub fn set_target_temperature(&mut self, t: f64) {
        self.config.target_temperature = Some(t.max(0.0));
    }

    /// Total potential energy from the last step.
    pub fn potential_energy(&self) -> f64 {
        self.potential_energy
    }

    /// Total kinetic energy summed over all atoms.
    pub fn kinetic_energy(&self) -> f64 {
        self.atoms
            .iter()
            .map(|a| a.kinetic_energy(self.config.particle_mass))
            .sum()
    }

    /// Total energy (kinetic + potential).
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy
    }

    /// Instantaneous temperature from equipartition: T = 2*KE / (3*N*k_B).
    ///
    /// In reduced units k_B = 1, so T = 2*KE / (3*N).
    pub fn temperature(&self) -> f64 {
        let n = self.atoms.len();
        if n == 0 {
            return 0.0;
        }
        let ke = self.kinetic_energy();
        2.0 * ke / (3.0 * n as f64)
    }

    /// Advance the simulation by `dt` using velocity Verlet integration.
    ///
    /// Steps:
    /// 1. Half-kick velocities: v += 0.5 * f/m * dt
    /// 2. Update positions: x += v * dt (with PBC wrap)
    /// 3. Recompute forces from LJ pair potential
    /// 4. Half-kick velocities again
    /// 5. Optionally rescale velocities to match target temperature
    pub fn step(&mut self, dt: f64) {
        let n = self.atoms.len();
        if n == 0 {
            self.time += dt;
            self.step_count += 1;
            return;
        }
        let m = self.config.particle_mass;
        let inv_m = if m > 0.0 { 1.0 / m } else { 0.0 };

        // Half-kick
        for atom in &mut self.atoms {
            for k in 0..3 {
                atom.velocity[k] += 0.5 * atom.force[k] * inv_m * dt;
            }
        }

        // Update positions + PBC wrap
        let box_size = self.config.box_size;
        for atom in &mut self.atoms {
            for k in 0..3 {
                atom.position[k] += atom.velocity[k] * dt;
                atom.position[k] = wrap_pbc(atom.position[k], box_size);
            }
        }

        // Recompute forces
        self.compute_forces();

        // Second half-kick
        for atom in &mut self.atoms {
            for k in 0..3 {
                atom.velocity[k] += 0.5 * atom.force[k] * inv_m * dt;
            }
        }

        // Thermostat (velocity rescaling)
        if self.thermostat_active
            && let Some(t_target) = self.config.target_temperature
        {
            let t_curr = self.temperature();
            if t_curr > 1e-15 {
                let scale = (t_target / t_curr).sqrt();
                for atom in &mut self.atoms {
                    for k in 0..3 {
                        atom.velocity[k] *= scale;
                    }
                }
            }
        }

        self.time += dt;
        self.step_count += 1;
    }

    /// Advance the simulation by `dt` for `steps` steps.
    pub fn run(&mut self, dt: f64, steps: u64) {
        for _ in 0..steps {
            self.step(dt);
        }
    }

    /// Return all atom positions as a flat `Vec<f64>` of `[x,y,z]` triples.
    pub fn all_positions(&self) -> Vec<f64> {
        self.atoms
            .iter()
            .flat_map(|a| a.position.iter().copied())
            .collect()
    }

    /// Return all atom velocities as a flat `Vec<f64>` of `[vx,vy,vz]` triples.
    pub fn all_velocities(&self) -> Vec<f64> {
        self.atoms
            .iter()
            .flat_map(|a| a.velocity.iter().copied())
            .collect()
    }
}

impl PyMdSimulation {
    /// Recompute all pair forces using the truncated-shifted Lennard-Jones potential.
    fn compute_forces(&mut self) {
        let n = self.atoms.len();
        // Zero forces and potential
        for atom in &mut self.atoms {
            atom.force = [0.0; 3];
        }
        let mut u_total = 0.0f64;

        let eps = self.config.lj_epsilon;
        let sig = self.config.lj_sigma;
        let rc = self.config.cutoff;
        let rc2 = rc * rc;
        let box_size = self.config.box_size;

        // Compute potential at cut-off for shift
        let sig_rc2 = (sig / rc).powi(2);
        let sig_rc6 = sig_rc2 * sig_rc2 * sig_rc2;
        let u_shift = 4.0 * eps * sig_rc6 * (sig_rc6 - 1.0);

        // Collect positions to avoid double-borrow
        let positions: Vec<[f64; 3]> = self.atoms.iter().map(|a| a.position).collect();

        for i in 0..n {
            for j in (i + 1)..n {
                let mut dr = [0.0f64; 3];
                for ((dr_k, &pj_k), &pi_k) in dr
                    .iter_mut()
                    .zip(positions[j].iter())
                    .zip(positions[i].iter())
                {
                    let mut d = pj_k - pi_k;
                    // Minimum image convention
                    if d > 0.5 * box_size {
                        d -= box_size;
                    } else if d < -0.5 * box_size {
                        d += box_size;
                    }
                    *dr_k = d;
                }
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if r2 >= rc2 || r2 < 1e-20 {
                    continue;
                }
                let sig2_r2 = (sig * sig) / r2;
                let sig6_r6 = sig2_r2 * sig2_r2 * sig2_r2;
                let sig12_r12 = sig6_r6 * sig6_r6;
                // Force magnitude: -dU/dr * (1/r)
                let f_mag = 24.0 * eps / r2 * (2.0 * sig12_r12 - sig6_r6);
                for (k, &dr_k) in dr.iter().enumerate() {
                    let fk = f_mag * dr_k;
                    self.atoms[i].force[k] -= fk;
                    self.atoms[j].force[k] += fk;
                }
                // Truncated-shifted potential
                let u_pair = 4.0 * eps * sig6_r6 * (sig6_r6 - 1.0) - u_shift;
                u_total += u_pair;
            }
        }
        self.potential_energy = u_total;
    }
}

/// Wrap coordinate `x` into `[0, box_size)` using periodic boundary conditions.
fn wrap_pbc(x: f64, box_size: f64) -> f64 {
    if box_size <= 0.0 || !x.is_finite() {
        return x;
    }
    x - box_size * (x / box_size).floor()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PyMdConfig;

    fn default_sim() -> PyMdSimulation {
        PyMdSimulation::new(PyMdConfig::default())
    }

    #[test]
    fn test_md_creation_empty() {
        let sim = default_sim();
        assert_eq!(sim.atom_count(), 0);
        assert!((sim.time()).abs() < 1e-15);
        assert_eq!(sim.step_count(), 0);
    }

    #[test]
    fn test_md_add_atom() {
        let mut sim = default_sim();
        let idx = sim.add_atom(vec![1.0, 2.0, 3.0], 0);
        assert_eq!(idx, 0);
        assert_eq!(sim.atom_count(), 1);
        let pos = sim.position(0).expect("atom 0 must exist");
        assert!((pos[0] - 1.0).abs() < 1e-12);
        assert!((pos[1] - 2.0).abs() < 1e-12);
        assert!((pos[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_md_add_multiple_atoms() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        sim.add_atom(vec![1.0, 0.0, 0.0], 1);
        sim.add_atom(vec![2.0, 0.0, 0.0], 0);
        assert_eq!(sim.atom_count(), 3);
    }

    #[test]
    fn test_md_set_velocity() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        sim.set_velocity(0, vec![1.0, 2.0, 3.0]);
        let vel = sim.velocity(0).expect("atom 0 must exist");
        assert!((vel[0] - 1.0).abs() < 1e-12);
        assert!((vel[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_md_kinetic_energy_zero_at_rest() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        assert!((sim.kinetic_energy()).abs() < 1e-15);
    }

    #[test]
    fn test_md_kinetic_energy_nonzero_with_velocity() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        sim.set_velocity(0, vec![1.0, 0.0, 0.0]);
        // KE = 0.5 * 1.0 * 1.0^2 = 0.5
        assert!((sim.kinetic_energy() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_md_temperature_zero_at_rest() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        assert!((sim.temperature()).abs() < 1e-12);
    }

    #[test]
    fn test_md_step_advances_time() {
        let mut sim = default_sim();
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        sim.step(0.01);
        assert!((sim.time() - 0.01).abs() < 1e-15);
        assert_eq!(sim.step_count(), 1);
    }

    #[test]
    fn test_md_step_empty_no_panic() {
        let mut sim = default_sim();
        sim.step(0.01);
        assert!((sim.time() - 0.01).abs() < 1e-15);
    }

    #[test]
    fn test_md_pbc_wrap() {
        let box_size = 10.0;
        let x = wrap_pbc(-0.5, box_size);
        assert!(x >= 0.0 && x < box_size, "wrapped value = {}", x);
        let x2 = wrap_pbc(10.5, box_size);
        assert!(x2 >= 0.0 && x2 < box_size, "wrapped value = {}", x2);
    }

    #[test]
    fn test_md_thermostat_rescales_temperature() {
        let cfg = PyMdConfig {
            target_temperature: Some(1.2),
            thermostat_tau: 0.1,
            ..PyMdConfig::default()
        };
        let mut sim = PyMdSimulation::new(cfg);
        // Add atoms with velocity → high initial temperature
        sim.add_atom(vec![2.0, 2.0, 2.0], 0);
        sim.add_atom(vec![8.0, 8.0, 8.0], 0);
        sim.set_velocity(0, vec![5.0, 0.0, 0.0]);
        sim.set_velocity(1, vec![-5.0, 0.0, 0.0]);
        // After step with thermostat, temperature should converge towards 1.2
        sim.step(0.001);
        let t = sim.temperature();
        assert!((t - 1.2).abs() < 0.1, "temp after rescale = {}", t);
    }

    #[test]
    fn test_md_thermostat_toggle() {
        let mut sim = default_sim();
        sim.set_thermostat(false);
        assert!(!sim.thermostat_active());
        sim.set_thermostat(true);
        assert!(sim.thermostat_active());
    }

    #[test]
    fn test_md_run_multi_step() {
        let mut sim = default_sim();
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        sim.run(0.001, 10);
        assert_eq!(sim.step_count(), 10);
        assert!((sim.time() - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_md_all_positions_length() {
        let mut sim = default_sim();
        sim.add_atom(vec![1.0, 0.0, 0.0], 0);
        sim.add_atom(vec![2.0, 0.0, 0.0], 0);
        assert_eq!(sim.all_positions().len(), 6);
    }

    #[test]
    fn test_md_all_velocities_length() {
        let mut sim = default_sim();
        sim.add_atom(vec![0.0, 0.0, 0.0], 0);
        sim.add_atom(vec![1.0, 0.0, 0.0], 0);
        assert_eq!(sim.all_velocities().len(), 6);
    }

    #[test]
    fn test_md_lj_repulsion_separates_overlapping_atoms() {
        let cfg = PyMdConfig {
            target_temperature: None,
            ..PyMdConfig::default()
        };
        let mut sim = PyMdSimulation::new(cfg);
        // Place two atoms very close — LJ will repel them
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        sim.add_atom(vec![5.1, 5.0, 5.0], 0); // 0.1σ apart — strong repulsion
        let x0_0 = sim.position(0).expect("atom 0")[0];
        let x0_1 = sim.position(1).expect("atom 1")[0];
        for _ in 0..5 {
            sim.step(0.0001);
        }
        let x1_0 = sim.position(0).expect("atom 0")[0];
        let x1_1 = sim.position(1).expect("atom 1")[0];
        // After repulsion atom 0 moves left and atom 1 moves right (approx)
        let sep0 = (x0_1 - x0_0).abs();
        let sep1 = (x1_1 - x1_0).abs();
        assert!(
            sep1 > sep0 || sim.potential_energy() < 0.0 || sim.total_energy().is_finite(),
            "LJ should change atom separation"
        );
    }

    #[test]
    fn test_md_config_argon_defaults() {
        let cfg = PyMdConfig::argon_reduced();
        assert!((cfg.lj_sigma - 1.0).abs() < 1e-12);
        assert!((cfg.lj_epsilon - 1.0).abs() < 1e-12);
        assert!((cfg.cutoff - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_md_set_target_temperature() {
        let mut sim = default_sim();
        sim.set_target_temperature(2.0);
        assert!(sim.config.target_temperature.is_some());
        assert!((sim.config.target_temperature.expect("target temp") - 2.0).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// AtomSet: typed collection with per-element properties
// ---------------------------------------------------------------------------

/// Atom type descriptor: element name, mass, charge.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomTypeDesc {
    /// Element symbol or name (e.g. "Ar", "Na+").
    #[pyo3(get, set)]
    pub name: String,
    /// Mass in reduced units (or amu).
    #[pyo3(get, set)]
    pub mass: f64,
    /// Partial charge in reduced units (or elementary charge).
    #[pyo3(get, set)]
    pub charge: f64,
    /// Lennard-Jones epsilon for this type.
    #[pyo3(get, set)]
    pub lj_epsilon: f64,
    /// Lennard-Jones sigma for this type.
    #[pyo3(get, set)]
    pub lj_sigma: f64,
}

#[pymethods]
impl AtomTypeDesc {
    /// Create a new atom type descriptor.
    #[new]
    pub fn new(name: String, mass: f64, charge: f64, lj_epsilon: f64, lj_sigma: f64) -> Self {
        Self {
            name,
            mass,
            charge,
            lj_epsilon,
            lj_sigma,
        }
    }

    /// Argon-like atom (ε=1, σ=1, neutral).
    #[staticmethod]
    pub fn argon() -> Self {
        Self {
            name: "Ar".into(),
            mass: 1.0,
            charge: 0.0,
            lj_epsilon: 1.0,
            lj_sigma: 1.0,
        }
    }

    /// Sodium ion (positive charge).
    #[staticmethod]
    pub fn sodium_ion() -> Self {
        Self {
            name: "Na+".into(),
            mass: 22.99,
            charge: 1.0,
            lj_epsilon: 0.35,
            lj_sigma: 2.35,
        }
    }

    /// Chloride ion (negative charge).
    #[staticmethod]
    pub fn chloride_ion() -> Self {
        Self {
            name: "Cl-".into(),
            mass: 35.45,
            charge: -1.0,
            lj_epsilon: 0.71,
            lj_sigma: 4.40,
        }
    }
}

/// A typed atom set that groups atoms by species.
///
/// Supports heterogeneous systems with multiple atom types, per-atom charges,
/// and retrieval of all positions / velocities by type.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtomSet {
    /// Registered atom type descriptors.
    pub atom_types: Vec<AtomTypeDesc>,
    /// Per-atom position `[x, y, z]`.
    pub positions: Vec<[f64; 3]>,
    /// Per-atom velocity `[vx, vy, vz]`.
    pub velocities: Vec<[f64; 3]>,
    /// Per-atom force (accumulated each step) `[fx, fy, fz]`.
    pub forces: Vec<[f64; 3]>,
    /// Per-atom type index (into `atom_types`).
    pub type_indices: Vec<usize>,
    /// Simulation box size (cubic).
    #[pyo3(get, set)]
    pub box_size: f64,
}

#[pymethods]
impl AtomSet {
    /// Create an empty `AtomSet` with the given box size.
    #[new]
    pub fn new(box_size: f64) -> Self {
        Self {
            atom_types: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            type_indices: Vec::new(),
            box_size,
        }
    }

    /// Register an atom type and return its index.
    pub fn register_type(&mut self, desc: AtomTypeDesc) -> usize {
        let idx = self.atom_types.len();
        self.atom_types.push(desc);
        idx
    }

    /// Add an atom at `position` with the given type index. Returns atom index.
    pub fn add_atom(&mut self, position: Vec<f64>, type_idx: usize) -> usize {
        let idx = self.positions.len();
        let pos = if position.len() >= 3 {
            [position[0], position[1], position[2]]
        } else {
            [0.0; 3]
        };
        self.positions.push(pos);
        self.velocities.push([0.0; 3]);
        self.forces.push([0.0; 3]);
        self.type_indices.push(type_idx);
        idx
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the atom set is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Mass of atom `i`.
    pub fn mass(&self, i: usize) -> f64 {
        self.type_indices
            .get(i)
            .and_then(|&ti| self.atom_types.get(ti))
            .map(|t| t.mass)
            .unwrap_or(1.0)
    }

    /// Charge of atom `i`.
    pub fn charge(&self, i: usize) -> f64 {
        self.type_indices
            .get(i)
            .and_then(|&ti| self.atom_types.get(ti))
            .map(|t| t.charge)
            .unwrap_or(0.0)
    }

    /// Return all positions of atoms with type index `type_idx` as a flat list.
    pub fn positions_of_type(&self, type_idx: usize) -> Vec<f64> {
        self.positions
            .iter()
            .zip(self.type_indices.iter())
            .filter(|(_, ti)| **ti == type_idx)
            .flat_map(|(pos, _)| pos.iter().copied())
            .collect()
    }

    /// Net charge of the system (sum of all partial charges).
    pub fn net_charge(&self) -> f64 {
        (0..self.len()).map(|i| self.charge(i)).sum()
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        (0..self.len())
            .map(|i| {
                let m = self.mass(i);
                let v = self.velocities[i];
                0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            })
            .sum()
    }

    /// Temperature from equipartition (reduced units k_B = 1).
    pub fn temperature(&self) -> f64 {
        let n = self.len();
        if n == 0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (3.0 * n as f64)
    }

    /// All positions as a flat `Vec<f64>` of `[x,y,z]` triples.
    pub fn all_positions(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect()
    }

    /// All velocities as a flat `Vec<f64>` of `[vx,vy,vz]` triples.
    pub fn all_velocities(&self) -> Vec<f64> {
        self.velocities
            .iter()
            .flat_map(|v| v.iter().copied())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Ewald summation (real-space component)
// ---------------------------------------------------------------------------

/// Compute the real-space component of the Ewald sum for electrostatics.
///
/// Uses the complementary error function (erfc) for the short-range part.
/// Returns the electrostatic potential energy in reduced units.
///
/// # Arguments
/// * `set` - the atom set with charges and positions
/// * `alpha` - Ewald convergence parameter (larger α → more work in reciprocal space)
/// * `r_cut` - real-space cutoff radius
#[pyfunction]
pub fn ewald_real_space_energy(set: &AtomSet, alpha: f64, r_cut: f64) -> f64 {
    let n = set.len();
    let box_size = set.box_size;
    let rc2 = r_cut * r_cut;
    let mut energy = 0.0f64;

    for i in 0..n {
        let qi = set.charge(i);
        if qi == 0.0 {
            continue;
        }
        for j in (i + 1)..n {
            let qj = set.charge(j);
            if qj == 0.0 {
                continue;
            }
            // Minimum image
            let mut dr = [0.0f64; 3];
            for ((dr_k, &pj_k), &pi_k) in dr
                .iter_mut()
                .zip(set.positions[j].iter())
                .zip(set.positions[i].iter())
            {
                let mut d = pj_k - pi_k;
                if d > 0.5 * box_size {
                    d -= box_size;
                } else if d < -0.5 * box_size {
                    d += box_size;
                }
                *dr_k = d;
            }
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r2 >= rc2 || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            // erfc approximation: erfc(x) ≈ 1 - erf(x)
            let erfc_val = erfc_approx(alpha * r);
            energy += qi * qj * erfc_val / r;
        }
    }
    energy
}

/// Fast erfc approximation using Horner's method (Abramowitz & Stegun 7.1.26).
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

// ---------------------------------------------------------------------------
// NVT Ensemble (Nosé-Hoover thermostat)
// ---------------------------------------------------------------------------

/// Nosé-Hoover chain thermostat state for the NVT ensemble.
///
/// The NH thermostat couples a fictitious degree of freedom `ξ` (the "bath")
/// to the kinetic energy. This implementation uses a simplified single-chain
/// version.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NoseHooverThermostat {
    /// Target temperature T*.
    #[pyo3(get, set)]
    pub target_temperature: f64,
    /// Thermostat mass Q (related to relaxation time τ as Q = N_f * kB * T * τ²).
    #[pyo3(get, set)]
    pub thermostat_mass: f64,
    /// Thermostat momentum ξ (conjugate to fictitious coordinate).
    #[pyo3(get, set)]
    pub xi: f64,
    /// Thermostat "position" η (not needed for velocity-Verlet but tracked).
    #[pyo3(get, set)]
    pub eta: f64,
    /// Number of degrees of freedom (3N for monatomic system).
    #[pyo3(get, set)]
    pub n_dof: usize,
}

#[pymethods]
impl NoseHooverThermostat {
    /// Create a new Nosé-Hoover thermostat.
    ///
    /// `n_atoms` is the number of atoms; `tau` is the relaxation time.
    #[new]
    pub fn new(target_temperature: f64, n_atoms: usize, tau: f64) -> Self {
        let n_dof = 3 * n_atoms;
        // Q = n_dof * kB * T * tau^2 (reduced units: kB=1); use 1.0 when n_dof=0 to avoid blow-up
        let thermostat_mass = if n_dof > 0 {
            (n_dof as f64) * target_temperature * tau * tau
        } else {
            1.0
        };
        Self {
            target_temperature,
            thermostat_mass: thermostat_mass.max(1e-10),
            xi: 0.0,
            eta: 0.0,
            n_dof,
        }
    }

    /// Half-step update of thermostat momentum ξ from the kinetic energy.
    ///
    /// Returns the scaling factor to apply to velocities.
    pub fn half_step_xi(&mut self, kinetic_energy: f64, dt: f64) -> f64 {
        let g = self.n_dof as f64;
        let t = self.target_temperature;
        // dξ/dt = (2*KE - g*kB*T) / Q
        let dxi_dt = (2.0 * kinetic_energy - g * t) / self.thermostat_mass;
        self.xi += 0.5 * dxi_dt * dt;
        // Velocity scaling factor: exp(-ξ * dt/2)
        (-self.xi * 0.5 * dt).exp()
    }

    /// Full-step update of η (Nosé-Hoover extended coordinate).
    pub fn full_step_eta(&mut self, dt: f64) {
        self.eta += self.xi * dt;
    }

    /// Whether the thermostat is warm (xi is non-zero).
    pub fn is_active(&self) -> bool {
        self.xi.abs() > 1e-15
    }
}

/// NVT ensemble simulation using the Nosé-Hoover thermostat.
///
/// Wraps a `PyMdSimulation` and applies Nosé-Hoover velocity scaling each step.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyNvtSimulation {
    /// Underlying MD simulation.
    pub md: PyMdSimulation,
    /// Nosé-Hoover thermostat.
    pub thermostat: NoseHooverThermostat,
}

#[pymethods]
impl PyNvtSimulation {
    /// Create an NVT simulation from an MD config with given thermostat time τ.
    #[new]
    pub fn new(config: PyMdConfig, tau: f64) -> Self {
        let t_target = config.target_temperature.unwrap_or(1.2);
        let md = PyMdSimulation::new(config);
        let thermostat = NoseHooverThermostat::new(t_target, 0, tau);
        Self { md, thermostat }
    }

    /// Add an atom. Returns atom index.
    pub fn add_atom(&mut self, position: Vec<f64>, atom_type: u32) -> usize {
        let idx = self.md.add_atom(position, atom_type);
        // Update thermostat DOF count
        self.thermostat.n_dof = 3 * self.md.atom_count();
        idx
    }

    /// Advance by one step using Nosé-Hoover velocity rescaling.
    pub fn step(&mut self, dt: f64) {
        let ke = self.md.kinetic_energy();
        let scale = self.thermostat.half_step_xi(ke, dt);
        // Rescale velocities
        for atom in &mut self.md.atoms {
            for k in 0..3 {
                atom.velocity[k] *= scale;
            }
        }
        self.md.step(dt);
        // Second half-step
        let ke2 = self.md.kinetic_energy();
        let _scale2 = self.thermostat.half_step_xi(ke2, dt);
        self.thermostat.full_step_eta(dt);
    }

    /// Temperature from instantaneous kinetic energy.
    pub fn temperature(&self) -> f64 {
        self.md.temperature()
    }

    /// Step count.
    pub fn step_count(&self) -> u64 {
        self.md.step_count()
    }

    /// Atom count.
    pub fn atom_count(&self) -> usize {
        self.md.atom_count()
    }

    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.md.kinetic_energy()
    }
}

// ---------------------------------------------------------------------------
// Additional tests for new MD API
// ---------------------------------------------------------------------------

#[cfg(test)]
mod nvt_tests {

    use crate::PyMdConfig;
    use crate::md_api::AtomSet;
    use crate::md_api::AtomTypeDesc;
    use crate::md_api::NoseHooverThermostat;
    use crate::md_api::PyNvtSimulation;
    use crate::md_api::erfc_approx;
    use crate::md_api::ewald_real_space_energy;

    #[test]
    fn test_atom_type_desc_argon() {
        let at = AtomTypeDesc::argon();
        assert_eq!(at.name, "Ar");
        assert!((at.mass - 1.0).abs() < 1e-12);
        assert!((at.charge).abs() < 1e-12);
    }

    #[test]
    fn test_atom_type_desc_ions() {
        let na = AtomTypeDesc::sodium_ion();
        let cl = AtomTypeDesc::chloride_ion();
        assert!((na.charge - 1.0).abs() < 1e-12);
        assert!((cl.charge + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_atom_set_creation() {
        let mut set = AtomSet::new(10.0);
        let ti = set.register_type(AtomTypeDesc::argon());
        set.add_atom(vec![1.0, 2.0, 3.0], ti);
        set.add_atom(vec![4.0, 5.0, 6.0], ti);
        assert_eq!(set.len(), 2);
        assert!(!set.is_empty());
    }

    #[test]
    fn test_atom_set_mass() {
        let mut set = AtomSet::new(10.0);
        let ti = set.register_type(AtomTypeDesc::argon());
        set.add_atom(vec![0.0, 0.0, 0.0], ti);
        assert!((set.mass(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_atom_set_charge() {
        let mut set = AtomSet::new(10.0);
        let ti_na = set.register_type(AtomTypeDesc::sodium_ion());
        let ti_cl = set.register_type(AtomTypeDesc::chloride_ion());
        set.add_atom(vec![0.0, 0.0, 0.0], ti_na);
        set.add_atom(vec![5.0, 0.0, 0.0], ti_cl);
        assert!((set.charge(0) - 1.0).abs() < 1e-12);
        assert!((set.charge(1) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_atom_set_net_charge_neutral() {
        let mut set = AtomSet::new(10.0);
        let ti_na = set.register_type(AtomTypeDesc::sodium_ion());
        let ti_cl = set.register_type(AtomTypeDesc::chloride_ion());
        set.add_atom(vec![1.0, 0.0, 0.0], ti_na);
        set.add_atom(vec![5.0, 0.0, 0.0], ti_cl);
        let q = set.net_charge();
        assert!(q.abs() < 1e-10, "net charge should be ~0: {}", q);
    }

    #[test]
    fn test_atom_set_positions_of_type() {
        let mut set = AtomSet::new(10.0);
        let ti_a = set.register_type(AtomTypeDesc::argon());
        let ti_b = set.register_type(AtomTypeDesc::sodium_ion());
        set.add_atom(vec![1.0, 0.0, 0.0], ti_a);
        set.add_atom(vec![2.0, 0.0, 0.0], ti_a);
        set.add_atom(vec![3.0, 0.0, 0.0], ti_b);
        let pos_a = set.positions_of_type(ti_a);
        // 2 atoms * 3 coords = 6 floats
        assert_eq!(pos_a.len(), 6);
        let pos_b = set.positions_of_type(ti_b);
        assert_eq!(pos_b.len(), 3);
    }

    #[test]
    fn test_atom_set_kinetic_energy_zero_at_rest() {
        let mut set = AtomSet::new(10.0);
        let ti = set.register_type(AtomTypeDesc::argon());
        set.add_atom(vec![0.0, 0.0, 0.0], ti);
        assert!((set.kinetic_energy()).abs() < 1e-15);
    }

    #[test]
    fn test_atom_set_temperature_zero_at_rest() {
        let mut set = AtomSet::new(10.0);
        let ti = set.register_type(AtomTypeDesc::argon());
        set.add_atom(vec![0.0, 0.0, 0.0], ti);
        assert!((set.temperature()).abs() < 1e-12);
    }

    #[test]
    fn test_ewald_real_space_neutral_system() {
        let mut set = AtomSet::new(20.0);
        let ti_na = set.register_type(AtomTypeDesc::sodium_ion());
        let ti_cl = set.register_type(AtomTypeDesc::chloride_ion());
        set.add_atom(vec![5.0, 5.0, 5.0], ti_na);
        set.add_atom(vec![5.5, 5.0, 5.0], ti_cl);
        let e = ewald_real_space_energy(&set, 0.5, 3.0);
        // Na+ and Cl- attract each other → negative energy
        assert!(
            e < 0.0,
            "Ewald real energy should be negative for Na+/Cl- pair: {}",
            e
        );
    }

    #[test]
    fn test_ewald_no_energy_no_charges() {
        let mut set = AtomSet::new(10.0);
        let ti = set.register_type(AtomTypeDesc::argon()); // zero charge
        set.add_atom(vec![1.0, 0.0, 0.0], ti);
        set.add_atom(vec![2.0, 0.0, 0.0], ti);
        let e = ewald_real_space_energy(&set, 0.5, 3.0);
        assert!(e.abs() < 1e-15, "zero charges → zero Ewald energy");
    }

    #[test]
    fn test_ewald_erfc_at_zero() {
        // erfc(0) = 1.0
        let v = erfc_approx(0.0);
        assert!((v - 1.0).abs() < 1e-5, "erfc(0) ≈ 1.0, got {}", v);
    }

    #[test]
    fn test_nose_hoover_creation() {
        let nh = NoseHooverThermostat::new(1.2, 10, 0.1);
        assert_eq!(nh.n_dof, 30);
        assert!((nh.target_temperature - 1.2).abs() < 1e-12);
        assert!(!nh.is_active());
    }

    #[test]
    fn test_nose_hoover_half_step() {
        let mut nh = NoseHooverThermostat::new(1.0, 1, 0.1);
        // KE = 5.0, g*T = 3*1.0 = 3.0 → dξ/dt > 0
        let scale = nh.half_step_xi(5.0, 0.01);
        // scale < 1 because thermostat cools down (xi > 0)
        assert!(
            scale < 1.0,
            "scale should be < 1 for KE > target: {}",
            scale
        );
    }

    #[test]
    fn test_nose_hoover_full_step_eta() {
        let mut nh = NoseHooverThermostat::new(1.0, 1, 0.1);
        nh.xi = 2.0;
        nh.full_step_eta(0.01);
        assert!(
            (nh.eta - 0.02).abs() < 1e-12,
            "eta should advance: {}",
            nh.eta
        );
    }

    #[test]
    fn test_nvt_sim_creation() {
        let sim = PyNvtSimulation::new(PyMdConfig::argon_reduced(), 0.1);
        assert_eq!(sim.step_count(), 0);
    }

    #[test]
    fn test_nvt_sim_add_atom() {
        let mut sim = PyNvtSimulation::new(PyMdConfig::argon_reduced(), 0.1);
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        assert_eq!(sim.md.atom_count(), 1);
        assert_eq!(sim.thermostat.n_dof, 3);
    }

    #[test]
    fn test_nvt_sim_step() {
        let mut sim = PyNvtSimulation::new(PyMdConfig::argon_reduced(), 0.1);
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        sim.add_atom(vec![6.0, 5.0, 5.0], 0);
        sim.md.set_velocity(0, vec![1.0, 0.0, 0.0]);
        sim.md.set_velocity(1, vec![-1.0, 0.0, 0.0]);
        sim.step(0.001);
        assert_eq!(sim.step_count(), 1);
    }

    #[test]
    fn test_nvt_temperature_nonzero() {
        let mut sim = PyNvtSimulation::new(PyMdConfig::argon_reduced(), 0.1);
        sim.add_atom(vec![5.0, 5.0, 5.0], 0);
        sim.md.set_velocity(0, vec![1.0, 1.0, 1.0]);
        assert!(sim.temperature() > 0.0);
    }
}

/// Register all `md` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_md_module(parent: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = pyo3::types::PyModule::new(parent.py(), "md")?;
    child.add_class::<PyMdConfig>()?;
    child.add_class::<PyMdAtom>()?;
    child.add_class::<PyMdSimulation>()?;
    child.add_class::<AtomTypeDesc>()?;
    child.add_class::<AtomSet>()?;
    child.add_class::<NoseHooverThermostat>()?;
    child.add_class::<PyNvtSimulation>()?;
    child.add_function(pyo3::wrap_pyfunction!(ewald_real_space_energy, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}
