// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared fixtures of the physics regression tests of the self-contained
//! drivers: a shifted-force Lennard-Jones force routine that reports the pair
//! virial and the potential energy, FCC lattices and argon-like fluids.

use super::analysis::{maxwell_boltzmann_velocities, rescale_to_temperature};
use super::forces::Atom;
use super::plain_sim::ForceEvaluation;

/// Lennard-Jones well depth of argon (kJ mol⁻¹); `ε / k_B = 119.8 K`.
pub(super) const ARGON_EPSILON: f64 = 0.996_07;
/// Lennard-Jones diameter of argon (Å).
pub(super) const ARGON_SIGMA: f64 = 3.405;
/// Mass of argon (amu).
pub(super) const ARGON_MASS: f64 = 39.948;

/// Forces on every atom (kJ mol⁻¹ Å⁻¹).
pub(super) type Forces = Vec<[f64; 3]>;

/// Forces, pair virial and potential energy of a pair-interacting system.
pub(super) struct PairResult {
    /// Force on every atom (kJ mol⁻¹ Å⁻¹).
    pub forces: Vec<[f64; 3]>,
    /// Pair virial `W = Σ_{i<j} r_ij · F_ij` (kJ mol⁻¹).
    pub virial: f64,
    /// Potential energy (kJ mol⁻¹).
    pub potential_energy: f64,
}

/// Shifted-force Lennard-Jones interaction, `u(r) = u_LJ(r) − u_LJ(r_c) −
/// (r − r_c) u_LJ'(r_c)` for `r < r_c`: energy and force both vanish at the
/// cutoff, so the dynamics conserves energy without cutoff noise.
#[derive(Clone, Copy)]
pub(super) struct ShiftedForceLj {
    /// Well depth (kJ mol⁻¹).
    pub epsilon: f64,
    /// Diameter (Å).
    pub sigma: f64,
    /// Cutoff radius (Å); must not exceed half the shortest box edge.
    pub cutoff: f64,
}

impl ShiftedForceLj {
    /// Argon-like interaction with the cutoff given in units of `σ`.
    pub fn argon(cutoff_in_sigma: f64) -> Self {
        Self {
            epsilon: ARGON_EPSILON,
            sigma: ARGON_SIGMA,
            cutoff: cutoff_in_sigma * ARGON_SIGMA,
        }
    }

    /// Plain 12-6 energy `u_LJ(r)`.
    fn energy_lj(&self, r: f64) -> f64 {
        let sr6 = (self.sigma / r).powi(6);
        4.0 * self.epsilon * (sr6 * sr6 - sr6)
    }

    /// Plain 12-6 repulsive force magnitude `−u_LJ'(r)`.
    fn force_lj(&self, r: f64) -> f64 {
        let sr6 = (self.sigma / r).powi(6);
        24.0 * self.epsilon * (2.0 * sr6 * sr6 - sr6) / r
    }

    /// Evaluate the system with minimum-image separations in the orthorhombic
    /// box `box_lengths` (the cutoff has to be below half of every edge).
    pub fn evaluate(&self, positions: &[[f64; 3]], box_lengths: &[f64; 3]) -> PairResult {
        let n = positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut virial = 0.0;
        let mut potential_energy = 0.0;
        let cutoff_sq = self.cutoff * self.cutoff;
        let u_cut = self.energy_lj(self.cutoff);
        let f_cut = self.force_lj(self.cutoff);
        for i in 0..n {
            for j in (i + 1)..n {
                let mut d = [0.0_f64; 3];
                for (a, da) in d.iter_mut().enumerate() {
                    let raw = positions[i][a] - positions[j][a];
                    *da = raw - box_lengths[a] * (raw / box_lengths[a]).round();
                }
                let r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if r_sq >= cutoff_sq || r_sq < 1e-24 {
                    continue;
                }
                let r = r_sq.sqrt();
                let magnitude = self.force_lj(r) - f_cut;
                potential_energy += self.energy_lj(r) - u_cut + f_cut * (r - self.cutoff);
                virial += magnitude * r;
                let scale = magnitude / r;
                for a in 0..3 {
                    forces[i][a] += scale * d[a];
                    forces[j][a] -= scale * d[a];
                }
            }
        }
        PairResult {
            forces,
            virial,
            potential_energy,
        }
    }

    /// Force routine of [`MdSim`](super::MdSim) drivers: `(forces, virial)`.
    pub fn virial_fn(&self) -> impl Fn(&[[f64; 3]], &[f64; 3]) -> (Forces, f64) + use<> {
        let lj = *self;
        move |positions, box_lengths| {
            let result = lj.evaluate(positions, box_lengths);
            (result.forces, result.virial)
        }
    }

    /// Force routine of [`MdSimulationPlain`](super::MdSimulationPlain)
    /// supplying forces, virial and potential energy.
    pub fn plain_fn(&self) -> impl FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation + use<> {
        let lj = *self;
        move |atoms, box_lengths| {
            let positions: Vec<[f64; 3]> = atoms.iter().map(|atom| atom.position).collect();
            let result = lj.evaluate(&positions, box_lengths);
            ForceEvaluation::new(result.forces)
                .with_virial(result.virial)
                .with_potential_energy(result.potential_energy)
        }
    }
}

/// Positions of an FCC crystal of `cells`³ cubic cells with the cubic lattice
/// constant `lattice_constant` (Å), and the edges of the periodic box.
pub(super) fn fcc_lattice(cells: usize, lattice_constant: f64) -> (Vec<[f64; 3]>, [f64; 3]) {
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut positions = Vec::with_capacity(4 * cells * cells * cells);
    for ix in 0..cells {
        for iy in 0..cells {
            for iz in 0..cells {
                for b in &basis {
                    positions.push([
                        (ix as f64 + b[0]) * lattice_constant,
                        (iy as f64 + b[1]) * lattice_constant,
                        (iz as f64 + b[2]) * lattice_constant,
                    ]);
                }
            }
        }
    }
    let edge = cells as f64 * lattice_constant;
    (positions, [edge; 3])
}

/// Argon-like fluid on an FCC lattice with `cells`³ cells at the reduced number
/// density `rho * = N σ³ / V`, with Maxwell-Boltzmann velocities at exactly
/// `temperature` (K) and zero total momentum.  Returns the atoms and the box.
pub(super) fn argon_fcc_fluid(
    cells: usize,
    reduced_density: f64,
    temperature: f64,
    seed: u64,
) -> (Vec<Atom>, [f64; 3]) {
    let lattice_constant = (4.0 / reduced_density).cbrt() * ARGON_SIGMA;
    let (positions, box_lengths) = fcc_lattice(cells, lattice_constant);
    let masses = vec![ARGON_MASS; positions.len()];
    let mut velocities = maxwell_boltzmann_velocities(&masses, temperature, seed);
    rescale_to_temperature(&mut velocities, &masses, temperature);
    let atoms = positions
        .into_iter()
        .zip(velocities)
        .map(|(position, velocity)| Atom::new(position, velocity, ARGON_MASS, 0.0))
        .collect();
    (atoms, box_lengths)
}

/// Reduced number density `N σ³ / V` of `n_atoms` argon atoms in `box_lengths`.
pub(super) fn reduced_density(n_atoms: usize, box_lengths: &[f64; 3]) -> f64 {
    n_atoms as f64 * ARGON_SIGMA.powi(3) / (box_lengths[0] * box_lengths[1] * box_lengths[2])
}
