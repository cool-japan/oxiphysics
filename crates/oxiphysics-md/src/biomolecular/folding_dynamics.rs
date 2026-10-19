// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gō-model Langevin dynamics of [`ProteinFolding`].
//!
//! The beads of the topology move on the structure-based (Gō) potential energy
//! surface of Clementi, Nymeyer & Onuchic (J. Mol. Biol. 298, 937, 2000), whose
//! global minimum is the reference ("native") structure:
//!
//! ```text
//! V = Σ_bonds ½ k_b (r − r0)²
//!   + Σ_native contacts ε [ 5 (r0/r)¹² − 6 (r0/r)¹⁰ ]
//!   + Σ_non-native pairs ε (σ/r)¹²,
//! ```
//!
//! with `r0` the separation in the reference structure, the stiffness
//! `k_b = 200 ε / Å²` (the `K_r = 100 ε` of Clementi *et al.* in the `K (r − r0)²`
//! convention), the contact well depth `ε` = [`ProteinFolding::contact_epsilon`]
//! and the non-native hard-core diameter `σ = 4 Å`.  Pairs joined by a bond, pairs
//! of consecutive atoms and the two ends of a bond angle carry no non-native
//! repulsion, and a native contact that is also a bond is not counted twice.
//! The angle and dihedral terms of the full Clementi model are not included.
//!
//! The beads are propagated by the BAOAB splitting of underdamped Langevin
//! dynamics (Leimkuhler & Matthews 2013) at the temperature
//! [`ProteinFolding::temperature`] with friction [`ProteinFolding::gamma`].
//!
//! # Units
//! Lengths in Å, time in ps, energy in kcal mol⁻¹, mass in amu, so velocities are in
//! Å/ps and forces in kcal mol⁻¹ Å⁻¹; one kcal mol⁻¹ is
//! [`KCAL_MOL_TO_AMU_A2_PS2`] = 418.4 amu·Å²/ps².

use super::types::ProteinFolding;
use crate::units::{K_B_AMU_A2_PS2, KCAL_MOL_TO_AMU_A2_PS2};
use oxiphysics_core::random::{Xoshiro256, box_muller};

#[cfg(test)]
mod tests;

/// Default bead mass (amu): the mean mass of an amino-acid residue.
pub const DEFAULT_BEAD_MASS: f64 = 110.0;

/// Default native-contact well depth `ε` (kcal mol⁻¹).
pub const DEFAULT_CONTACT_EPSILON: f64 = 1.0;

/// Default Langevin friction coefficient `γ` (ps⁻¹).
pub const DEFAULT_FRICTION: f64 = 1.0;

/// Default seed of the noise generator.
pub const DEFAULT_FOLDING_SEED: u64 = 0x601D_C0DE_2000_0001;

/// Bond force constant in units of `ε / Å²` (`V = ½ k (r − r0)²`).
pub const BOND_STIFFNESS_PER_EPSILON: f64 = 200.0;

/// Hard-core diameter `σ` (Å) of the non-native repulsion `ε (σ/r)¹²`.
pub const NON_NATIVE_SIGMA: f64 = 4.0;

/// A native contact counts as formed while its separation stays below this
/// multiple of its native separation (the standard definition of `Q`).
pub const CONTACT_FORMED_FACTOR: f64 = 1.2;

/// Separations below this value (Å) are treated as coincident: the direction of the
/// pair force is undefined there.
const COINCIDENT_DISTANCE: f64 = 1.0e-9;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Add `scale · d` to `forces[a]` and subtract it from `forces[b]`.
fn apply_pair_force(forces: &mut [[f64; 3]], a: usize, b: usize, d: [f64; 3], scale: f64) {
    for axis in 0..3 {
        forces[a][axis] += scale * d[axis];
        forces[b][axis] -= scale * d[axis];
    }
}

/// Sorted per-atom lists of the partners `j > i` of every atom `i` in a pair class.
struct PairTable {
    partners: Vec<Vec<usize>>,
}

impl PairTable {
    fn new(n_atoms: usize) -> Self {
        Self {
            partners: vec![Vec::new(); n_atoms],
        }
    }

    fn insert(&mut self, a: usize, b: usize) {
        let (lo, hi) = (a.min(b), a.max(b));
        if lo != hi
            && let Some(list) = self.partners.get_mut(lo)
        {
            list.push(hi);
        }
    }

    fn finish(&mut self) {
        for list in &mut self.partners {
            list.sort_unstable();
            list.dedup();
        }
    }

    /// Whether the pair `(i, j)` with `i < j` is in the table.
    fn contains(&self, i: usize, j: usize) -> bool {
        self.partners
            .get(i)
            .is_some_and(|list| list.binary_search(&j).is_ok())
    }
}

/// The pair classes that decide which interaction a bead pair carries.
struct PairClasses {
    /// Pairs joined by a bond.
    bonded: PairTable,
    /// Pairs without non-native repulsion: bonded pairs, consecutive atoms, the
    /// two ends of a bond angle and the native contacts.
    no_repulsion: PairTable,
}

impl ProteinFolding {
    /// Replace the noise seed: the generator state becomes the one of a
    /// `Xoshiro256` seeded with `seed`.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.rng_state = Xoshiro256::new(seed).state;
        self
    }

    /// Force constant `k_b = 200 ε / Å²` (kcal mol⁻¹ Å⁻²) of the harmonic bonds.
    pub fn bond_force_constant(&self) -> f64 {
        BOND_STIFFNESS_PER_EPSILON * self.contact_epsilon
    }

    fn pair_classes(&self) -> PairClasses {
        let n_atoms = self.topology.positions.len();
        let mut bonded = PairTable::new(n_atoms);
        let mut no_repulsion = PairTable::new(n_atoms);
        let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
        for &[a, b] in &self.topology.bonds {
            if a < n_atoms && b < n_atoms && a != b {
                bonded.insert(a, b);
                no_repulsion.insert(a, b);
                neighbours[a].push(b);
                neighbours[b].push(a);
            }
        }
        // 1-3 pairs: the two neighbours of a common atom.
        for around in &neighbours {
            for (position, &first) in around.iter().enumerate() {
                for &second in &around[position + 1..] {
                    no_repulsion.insert(first, second);
                }
            }
        }
        // Consecutive atoms and the native contacts.
        for atom in 1..n_atoms {
            no_repulsion.insert(atom - 1, atom);
        }
        for &[i, j] in &self.native_contacts {
            if i < n_atoms && j < n_atoms {
                no_repulsion.insert(i, j);
            }
        }
        bonded.finish();
        no_repulsion.finish();
        PairClasses {
            bonded,
            no_repulsion,
        }
    }

    /// Potential energy (kcal mol⁻¹) and forces (kcal mol⁻¹ Å⁻¹) of the current
    /// bead positions on the Gō surface.
    fn evaluate(&self) -> (f64, Vec<[f64; 3]>) {
        let positions = &self.topology.positions;
        let reference = &self.ref_positions;
        let n_atoms = positions.len();
        let epsilon = self.contact_epsilon;
        let k_bond = self.bond_force_constant();
        let classes = self.pair_classes();
        let mut forces = vec![[0.0_f64; 3]; n_atoms];
        let mut energy = 0.0_f64;

        // Harmonic bonds about their native lengths.
        for &[a, b] in &self.topology.bonds {
            if a == b {
                continue;
            }
            let (Some(&x_a), Some(&x_b), Some(&ref_a), Some(&ref_b)) = (
                positions.get(a),
                positions.get(b),
                reference.get(a),
                reference.get(b),
            ) else {
                continue;
            };
            let d = sub(x_a, x_b);
            let r = norm(d);
            if r < COINCIDENT_DISTANCE {
                continue;
            }
            let stretch = r - norm(sub(ref_a, ref_b));
            energy += 0.5 * k_bond * stretch * stretch;
            apply_pair_force(&mut forces, a, b, d, -k_bond * stretch / r);
        }

        // Native contacts: 10-12 potential with its minimum -ε at the native distance.
        for &[i, j] in &self.native_contacts {
            let (lo, hi) = (i.min(j), i.max(j));
            if lo == hi || classes.bonded.contains(lo, hi) {
                continue;
            }
            let (Some(&x_i), Some(&x_j), Some(&ref_i), Some(&ref_j)) = (
                positions.get(i),
                positions.get(j),
                reference.get(i),
                reference.get(j),
            ) else {
                continue;
            };
            let d = sub(x_i, x_j);
            let r = norm(d);
            if r < COINCIDENT_DISTANCE {
                continue;
            }
            let q = norm(sub(ref_i, ref_j)) / r;
            let q10 = q.powi(10);
            let q12 = q10 * q * q;
            energy += epsilon * (5.0 * q12 - 6.0 * q10);
            apply_pair_force(&mut forces, i, j, d, 60.0 * epsilon * (q12 - q10) / (r * r));
        }

        // Non-native pairs repel with ε (σ/r)¹².
        for i in 0..n_atoms {
            for j in (i + 1)..n_atoms {
                if classes.no_repulsion.contains(i, j) {
                    continue;
                }
                let d = sub(positions[i], positions[j]);
                let r = norm(d);
                if r < COINCIDENT_DISTANCE {
                    continue;
                }
                let repulsion = epsilon * (NON_NATIVE_SIGMA / r).powi(12);
                energy += repulsion;
                apply_pair_force(&mut forces, i, j, d, 12.0 * repulsion / (r * r));
            }
        }
        (energy, forces)
    }

    /// Potential energy of the current conformation on the Gō surface (kcal mol⁻¹).
    ///
    /// The native structure is the global minimum; its energy is `-ε` times the
    /// number of native contacts that are not bonds, plus a negligible non-native
    /// repulsion.
    pub fn potential_energy(&self) -> f64 {
        self.evaluate().0
    }

    /// Forces on the beads in the current conformation (kcal mol⁻¹ Å⁻¹): the negative
    /// gradient of [`ProteinFolding::potential_energy`].
    pub fn forces(&self) -> Vec<[f64; 3]> {
        self.evaluate().1
    }

    /// Kinetic energy of the bead velocities (kcal mol⁻¹).
    pub fn kinetic_energy(&self) -> f64 {
        let sum_v_sq: f64 = self
            .velocities
            .iter()
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum();
        0.5 * self.bead_mass * sum_v_sq / KCAL_MOL_TO_AMU_A2_PS2
    }

    /// Instantaneous kinetic temperature `2 KE / (3 N k_B)` (K); `0.0` without beads.
    pub fn kinetic_temperature(&self) -> f64 {
        let n_atoms = self.velocities.len();
        if n_atoms == 0 {
            return 0.0;
        }
        let kinetic_amu = self.kinetic_energy() * KCAL_MOL_TO_AMU_A2_PS2;
        2.0 * kinetic_amu / (3.0 * n_atoms as f64 * K_B_AMU_A2_PS2)
    }

    /// The Xoshiro256 generator continuing from the stored state; an all-zero state
    /// (which would make the generator output zeros forever) is re-seeded.
    fn noise_generator(&self) -> Xoshiro256 {
        if self.rng_state == [0; 4] {
            Xoshiro256::new(DEFAULT_FOLDING_SEED)
        } else {
            Xoshiro256 {
                state: self.rng_state,
            }
        }
    }

    /// Draw the velocities of all beads from the Maxwell-Boltzmann distribution at
    /// [`ProteinFolding::temperature`] (`σ_v = √(k_B T / m)`).
    pub fn thermalize_velocities(&mut self) {
        let n_atoms = self.topology.positions.len();
        let sigma_v = self.thermal_velocity_scale();
        let mut rng = self.noise_generator();
        let mut velocities = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let (gx, gy) = box_muller(&mut rng);
            let (gz, _) = box_muller(&mut rng);
            velocities.push([sigma_v * gx, sigma_v * gy, sigma_v * gz]);
        }
        self.rng_state = rng.state;
        self.velocities = velocities;
    }

    /// `√(k_B T / m)` in Å/ps at the current temperature.
    fn thermal_velocity_scale(&self) -> f64 {
        (K_B_AMU_A2_PS2 * self.temperature.max(0.0) / self.bead_mass).sqrt()
    }

    /// One BAOAB Langevin step of length `dt` (ps) at temperature `T` with friction `γ`.
    ///
    /// ```text
    /// B: v += ½ Δt F(x)/m        A: x += ½ Δt v
    /// O: v  = e^{-γΔt} v + √((1 − e^{-2γΔt}) k_B T / m) ξ
    /// A: x += ½ Δt v            B: v += ½ Δt F(x)/m
    /// ```
    fn langevin_step(&mut self) {
        let dt = self.dt;
        let half = 0.5 * dt;
        // Acceleration (Å/ps²) per unit force (kcal mol⁻¹ Å⁻¹).
        let acceleration = KCAL_MOL_TO_AMU_A2_PS2 / self.bead_mass;

        let (_, forces) = self.evaluate();
        for ((velocity, position), force) in self
            .velocities
            .iter_mut()
            .zip(self.topology.positions.iter_mut())
            .zip(forces.iter())
        {
            for ((v, x), f) in velocity.iter_mut().zip(position.iter_mut()).zip(force) {
                *v += half * acceleration * f;
                *x += half * *v;
            }
        }

        let decay = (-self.gamma.max(0.0) * dt).exp();
        let sigma_v = (1.0 - decay * decay).sqrt() * self.thermal_velocity_scale();
        let mut rng = self.noise_generator();
        for (velocity, position) in self
            .velocities
            .iter_mut()
            .zip(self.topology.positions.iter_mut())
        {
            let (gx, gy) = box_muller(&mut rng);
            let (gz, _) = box_muller(&mut rng);
            for ((v, x), noise) in velocity
                .iter_mut()
                .zip(position.iter_mut())
                .zip([gx, gy, gz])
            {
                *v = decay * *v + sigma_v * noise;
                *x += half * *v;
            }
        }
        self.rng_state = rng.state;

        let (_, forces) = self.evaluate();
        for (velocity, force) in self.velocities.iter_mut().zip(forces.iter()) {
            for (v, f) in velocity.iter_mut().zip(force) {
                *v += half * acceleration * f;
            }
        }
    }

    /// Advance the simulation by one Langevin MD step and record the native-contact
    /// fraction `Q`.
    ///
    /// The beads move on the Gō surface (see the [module documentation](self)) with a
    /// BAOAB Langevin integrator at [`ProteinFolding::temperature`] (which
    /// [`ProteinFolding::quench`] can change between steps), friction
    /// [`ProteinFolding::gamma`] and time step [`ProteinFolding::dt`] (ps).  The
    /// velocities and noise generator state are kept in the struct, so a run is
    /// reproducible from the seed.  A velocity list of the wrong length is reset
    /// to rest.
    ///
    /// Nothing moves (but `Q` is still recorded) for a non-positive or non-finite
    /// `dt` or `bead_mass`.
    pub fn step(&mut self) {
        let n_atoms = self.topology.positions.len();
        if self.velocities.len() != n_atoms {
            self.velocities = vec![[0.0; 3]; n_atoms];
        }
        let valid = self.dt.is_finite()
            && self.dt > 0.0
            && self.bead_mass.is_finite()
            && self.bead_mass > 0.0;
        if valid && n_atoms > 0 {
            self.langevin_step();
        }
        self.step_count += 1;
        let q = self.contact_fraction_q();
        self.q_history.push(q);
    }

    /// Run `n_steps` Langevin steps, recording `Q` after each.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }
}
