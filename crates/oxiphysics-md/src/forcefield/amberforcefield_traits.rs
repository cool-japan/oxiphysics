//! # AmberForceField - Trait Implementations
//!
//! This module contains trait implementations for `AmberForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! The `ForceField` implementation evaluates the non-bonded interaction with the
//! AMBER exclusion rules: 1-2 and 1-3 neighbours are skipped, 1-4 neighbours are
//! scaled, and the partial charges interact through a reaction-field Coulomb term.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::{CellList, PeriodicBox};

use super::functions::ForceField;
use super::types::AmberForceField;

impl Default for AmberForceField {
    fn default() -> Self {
        Self::new()
    }
}

/// Coefficients of the reaction-field Coulomb interaction.
///
/// For unit charges the pair energy is `1/r + k_rf r² − c_rf` and the repulsive
/// force `1/r² − 2 k_rf r` for `r` below the cutoff `r_c`, and zero beyond.  The
/// continuum outside the cutoff has the relative permittivity `ε_rf`, and with
/// `ρ = 1/ε_rf`
///
/// ```text
/// k_rf = (1 − ρ) / ((2 + ρ) r_c³),      c_rf = 3 / ((2 + ρ) r_c),
/// ```
///
/// which makes the energy vanish at `r_c`.  For `ε_rf → ∞` (`ρ = 0`) the force
/// vanishes there too.  For finite `ε_rf` these are the coefficients of
/// [`crate::electrostatics::ReactionFieldElectrostatics`].
#[derive(Debug, Clone, Copy)]
struct ReactionField {
    cutoff: f64,
    k_rf: f64,
    c_rf: f64,
}

impl ReactionField {
    /// Coefficients for the permittivity `epsilon_rf`; a non-finite or
    /// non-positive value selects the conducting boundary `ε_rf = ∞`.
    fn new(epsilon_rf: f64, cutoff: f64) -> Self {
        let rho = if epsilon_rf.is_finite() && epsilon_rf > 0.0 {
            1.0 / epsilon_rf
        } else {
            0.0
        };
        Self {
            cutoff,
            k_rf: (1.0 - rho) / ((2.0 + rho) * cutoff.powi(3)),
            c_rf: 3.0 / ((2.0 + rho) * cutoff),
        }
    }

    /// Energy and repulsive-positive force of two unit charges at distance `r`.
    fn unit_energy_and_force(&self, r: f64) -> (f64, f64) {
        (
            1.0 / r + self.k_rf * r * r - self.c_rf,
            1.0 / (r * r) - 2.0 * self.k_rf * r,
        )
    }
}

impl AmberForceField {
    /// Accumulate the non-bonded forces and return `(energy, virial)`.
    ///
    /// Every pair within the cutoff that is not excluded interacts through the
    /// pair potential registered for its atom types and through the
    /// reaction-field Coulomb interaction of its partial charges.  The 1-4
    /// pairs have both energies (and forces) scaled by
    /// [`AmberForceField::scaling_14`].  With the minimum-image separation
    /// `dr = r_j − r_i` and the repulsive-positive force magnitude `f = −dE/dr`
    /// the force on `j` is `f·dr/|dr|` and the pair adds `dr·F_j = f·|dr|` to
    /// the virial, so the virial depends on separations only.
    fn accumulate_nonbonded(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        let pair_cutoff = self.pair.max_cutoff();
        let coulomb_cutoff = if self.coulomb_cutoff > 0.0 {
            self.coulomb_cutoff
        } else {
            pair_cutoff
        };
        let electrostatics = self.coulomb_k != 0.0
            && coulomb_cutoff > 0.0
            && atoms.charges.iter().any(|&charge| charge != 0.0);
        let cutoff = if electrostatics {
            pair_cutoff.max(coulomb_cutoff)
        } else {
            pair_cutoff
        };
        if cutoff <= 0.0 || atoms.is_empty() {
            return (0.0, 0.0);
        }
        let reaction_field =
            electrostatics.then(|| ReactionField::new(self.reaction_field_epsilon, coulomb_cutoff));
        let (exclusions, pairs_14) = self.effective_topology();
        let cell_list = CellList::build(&atoms.positions, pbox, cutoff);
        let pairs = cell_list.neighbor_pairs(&atoms.positions, pbox);
        let mut energy = 0.0;
        let mut virial = 0.0;
        for (i, j, dr, dist) in pairs {
            if dist <= 1e-15 || exclusions.is_excluded(i, j) {
                continue;
            }
            let one_four = pairs_14.contains(&(i.min(j), i.max(j)));
            let mut pair_energy = 0.0;
            let mut pair_force = 0.0;
            if let Some(potential) = self.pair.find_potential(atoms.types[i], atoms.types[j])
                && dist < potential.cutoff()
            {
                let (e, f) = (potential.energy(dist), potential.force(dist));
                if one_four {
                    pair_energy += self.scaling_14.scale_lj(e);
                    pair_force += self.scaling_14.scale_lj(f);
                } else {
                    pair_energy += e;
                    pair_force += f;
                }
            }
            if let Some(field) = &reaction_field
                && dist < field.cutoff
            {
                let charge_product = self.coulomb_k * atoms.charges[i] * atoms.charges[j];
                if charge_product != 0.0 {
                    let (unit_energy, unit_force) = field.unit_energy_and_force(dist);
                    let (e, f) = (charge_product * unit_energy, charge_product * unit_force);
                    if one_four {
                        pair_energy += self.scaling_14.scale_coulomb(e);
                        pair_force += self.scaling_14.scale_coulomb(f);
                    } else {
                        pair_energy += e;
                        pair_force += f;
                    }
                }
            }
            let f_vec = dr * (pair_force / dist);
            virial += dr.dot(&f_vec);
            atoms.forces[i] -= f_vec;
            atoms.forces[j] += f_vec;
            energy += pair_energy;
        }
        (energy, virial)
    }
}

impl ForceField for AmberForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        energy += self.accumulate_nonbonded(atoms, pbox).0;
        energy += self.bonded.compute_forces(atoms, pbox);
        energy += self.angle.compute_forces(atoms, pbox);
        energy += self.dihedral.compute_forces(atoms, pbox);
        energy
    }

    /// Sum of the component virials: the non-bonded and bond terms use their
    /// analytic pair sums, the angle and dihedral terms the energy-based
    /// default of [`ForceField::compute_forces_virial`].
    fn compute_forces_virial(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        let (e_nonbonded, w_nonbonded) = self.accumulate_nonbonded(atoms, pbox);
        let (e_bond, w_bond) = self.bonded.compute_forces_virial(atoms, pbox);
        let (e_angle, w_angle) = self.angle.compute_forces_virial(atoms, pbox);
        let (e_dihedral, w_dihedral) = self.dihedral.compute_forces_virial(atoms, pbox);
        (
            e_nonbonded + e_bond + e_angle + e_dihedral,
            w_nonbonded + w_bond + w_angle + w_dihedral,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electrostatics::{COULOMB_K, ReactionFieldElectrostatics};
    use crate::forcefield::virial::scaling_virial;
    use crate::forcefield::{OneFourScaling, PairForceField};
    use crate::neighbor::ExclusionList;
    use crate::potential::{LennardJones, Potential};
    use oxiphysics_core::math::Vec3;
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    const LJ_EPSILON: f64 = 0.4;
    const LJ_SIGMA: f64 = 3.4;
    const CUTOFF: f64 = 12.0;
    const BOND_K: f64 = 1500.0;
    const BOND_R0: f64 = 1.54;
    const ANGLE_K: f64 = 250.0;
    const BOX: f64 = 50.0;
    const CHARGES: [f64; 6] = [0.30, -0.20, -0.10, 0.25, -0.15, 0.10];
    /// Central-difference step for the numerical gradients.
    const FD_STEP: f64 = 1.0e-5;

    fn tetrahedral_angle() -> f64 {
        112.0_f64.to_radians()
    }

    /// Position of the atom `d` bonded to `c` at distance `bond`, with the bond
    /// angle `b-c-d` and the dihedral angle `a-b-c-d` (natural extension of
    /// reference frames).
    fn place_atom(a: Vec3, b: Vec3, c: Vec3, bond: f64, angle: f64, dihedral: f64) -> Vec3 {
        let bc = (c - b).normalize();
        let normal = (b - a).cross(&bc).normalize();
        let in_plane = normal.cross(&bc);
        c + bc * (-bond * angle.cos())
            + in_plane * (bond * angle.sin() * dihedral.cos())
            + normal * (bond * angle.sin() * dihedral.sin())
    }

    /// United-atom alkyl chain around the box centre; every atom is displaced
    /// by a random vector of at most `jitter` per axis so that all bonded
    /// terms exert non-zero forces.
    fn chain_positions(n_atoms: usize, jitter: f64, seed: u64) -> Vec<Vec3> {
        let centre = Vec3::new(BOX / 2.0, BOX / 2.0, BOX / 2.0);
        let theta = tetrahedral_angle();
        let mut positions = vec![centre, centre + Vec3::new(BOND_R0, 0.0, 0.0)];
        positions
            .push(positions[1] + Vec3::new(-BOND_R0 * theta.cos(), BOND_R0 * theta.sin(), 0.0));
        let dihedrals = [65.0_f64, 190.0, -70.0, 175.0, 60.0];
        for k in 3..n_atoms {
            let next = place_atom(
                positions[k - 3],
                positions[k - 2],
                positions[k - 1],
                BOND_R0,
                theta,
                dihedrals[(k - 3) % dihedrals.len()].to_radians(),
            );
            positions.push(next);
        }
        let mut rng = StdRng::seed_from_u64(seed);
        positions
            .into_iter()
            .map(|p| {
                p + Vec3::new(
                    rng.random_range(-jitter..jitter),
                    rng.random_range(-jitter..jitter),
                    rng.random_range(-jitter..jitter),
                )
            })
            .collect()
    }

    fn chain_atoms(positions: &[Vec3]) -> AtomSet {
        let mut atoms = AtomSet::new();
        for (i, p) in positions.iter().enumerate() {
            atoms.add_atom(*p, Vec3::zeros(), 12.0, CHARGES[i % CHARGES.len()], 0);
        }
        atoms
    }

    /// Chain with bonds, angles and dihedrals along the backbone and a
    /// Lennard-Jones interaction of all atoms; the topology is not built yet.
    fn chain_force_field(n_atoms: usize) -> AmberForceField {
        let mut ff = AmberForceField::new();
        ff.pair
            .add_interaction(0, 0, LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF));
        for i in 0..n_atoms - 1 {
            ff.bonded.add_bond(i, i + 1, BOND_K, BOND_R0);
        }
        for i in 0..n_atoms.saturating_sub(2) {
            ff.angle
                .add_angle(i, i + 1, i + 2, ANGLE_K, tetrahedral_angle());
        }
        for i in 0..n_atoms.saturating_sub(3) {
            ff.dihedral
                .add_dihedral(i, i + 1, i + 2, i + 3, 2.0, 0.0, 3);
        }
        ff
    }

    /// Sum of the bond, angle and dihedral energies.
    fn bonded_energy(ff: &AmberForceField, atoms: &AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut scratch = atoms.clone();
        scratch.clear_forces();
        ff.bonded.compute_forces(&mut scratch, pbox)
            + ff.angle.compute_forces(&mut scratch, pbox)
            + ff.dihedral.compute_forces(&mut scratch, pbox)
    }

    /// Conducting-boundary reaction-field energy of two charges at `r`.
    fn reaction_field_energy(q_i: f64, q_j: f64, r: f64) -> f64 {
        COULOMB_K * q_i * q_j * (1.0 / r + r * r / (2.0 * CUTOFF.powi(3)) - 1.5 / CUTOFF)
    }

    fn pair_distance(atoms: &AtomSet, i: usize, j: usize) -> f64 {
        (atoms.positions[j] - atoms.positions[i]).norm()
    }

    /// Independent non-bonded reference for a chain: atoms up to two bonds
    /// apart do not interact, atoms three bonds apart interact with the AMBER
    /// 1-4 scale factors and all others in full.
    fn reference_nonbonded_energy(atoms: &AtomSet) -> f64 {
        let lj = LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF);
        let mut energy = 0.0;
        for i in 0..atoms.len() {
            for j in (i + 1)..atoms.len() {
                let separation = j - i;
                if separation <= 2 {
                    continue;
                }
                let r = pair_distance(atoms, i, j);
                let (lj_scale, coulomb_scale) = if separation == 3 {
                    (0.5, 1.0 / 1.2)
                } else {
                    (1.0, 1.0)
                };
                energy += lj_scale * lj.energy(r)
                    + coulomb_scale * reaction_field_energy(atoms.charges[i], atoms.charges[j], r);
            }
        }
        energy
    }

    fn total_energy(ff: &AmberForceField, positions: &[Vec3], pbox: &PeriodicBox) -> f64 {
        let mut atoms = chain_atoms(positions);
        ff.compute_forces(&mut atoms, pbox)
    }

    fn numeric_forces(ff: &AmberForceField, positions: &[Vec3], pbox: &PeriodicBox) -> Vec<Vec3> {
        let mut out = vec![Vec3::zeros(); positions.len()];
        for atom in 0..positions.len() {
            for axis in 0..3 {
                let mut plus = positions.to_vec();
                let mut minus = positions.to_vec();
                plus[atom][axis] += FD_STEP;
                minus[atom][axis] -= FD_STEP;
                out[atom][axis] = -(total_energy(ff, &plus, pbox) - total_energy(ff, &minus, pbox))
                    / (2.0 * FD_STEP);
            }
        }
        out
    }

    fn two_atoms(separation: f64, q_i: f64, q_j: f64) -> AtomSet {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(20.0, 20.0, 20.0), Vec3::zeros(), 12.0, q_i, 0);
        atoms.add_atom(
            Vec3::new(20.0 + separation, 20.0, 20.0),
            Vec3::zeros(),
            12.0,
            q_j,
            0,
        );
        atoms
    }

    #[test]
    fn bonded_atoms_do_not_feel_the_lennard_jones_potential() {
        let atoms = two_atoms(BOND_R0, 0.0, 0.0);
        let pbox = PeriodicBox::cubic(BOX);
        let mut amber = AmberForceField::new();
        amber
            .pair
            .add_interaction(0, 0, LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF));
        amber.bonded.add_bond(0, 1, BOND_K, BOND_R0);

        // The bond sits at its equilibrium length: with the pair excluded
        // there is neither energy nor force.  This holds with and without an
        // explicit build_topology() call.
        for build in [true, false] {
            let mut ff = amber.clone();
            if build {
                ff.build_topology();
            }
            let mut a = atoms.clone();
            let energy = ff.compute_forces(&mut a, &pbox);
            assert!(energy.abs() < 1.0e-12, "build = {build}: energy {energy}");
            // The positions carry rounding noise of about 1e-15 A, so the
            // bond (k = 1500) feels a force of order 1e-12; the Lennard-Jones
            // force this test guards against is of order 1e5.
            assert!(
                a.forces.iter().all(|f| f.norm() < 1.0e-9),
                "build = {build}"
            );
        }

        // The plain pair force field puts the same atoms deep into the
        // repulsive core, which is what the exclusion removes.
        let mut plain = PairForceField::new();
        plain.add_interaction(0, 0, LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF));
        let mut b = atoms.clone();
        let repulsion = plain.compute_forces(&mut b, &pbox);
        assert!(repulsion > 1.0e3, "plain LJ energy {repulsion}");
    }

    #[test]
    fn butane_energy_is_the_bonded_terms_plus_the_scaled_one_four_pair() {
        let positions = chain_positions(4, 0.04, 11);
        let mut atoms = chain_atoms(&positions);
        let mut ff = chain_force_field(4);
        ff.build_topology();
        let pbox = PeriodicBox::cubic(BOX);
        let total = ff.compute_forces(&mut atoms, &pbox);

        // Only atoms 0 and 3 are neither bonded nor 1-3 neighbours.
        let r03 = pair_distance(&atoms, 0, 3);
        let lj = LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF);
        let one_four = 0.5 * lj.energy(r03)
            + reaction_field_energy(atoms.charges[0], atoms.charges[3], r03) / 1.2;
        let expected = bonded_energy(&ff, &atoms, &pbox) + one_four;
        assert!(
            (total - expected).abs() < 1.0e-9 * (1.0 + expected.abs()),
            "total {total} vs bonded + 1-4 pair {expected}"
        );
        assert!(
            (one_four - reference_nonbonded_energy(&atoms)).abs() < 1.0e-12,
            "reference helper disagrees with the explicit 1-4 energy"
        );
        assert!(one_four.abs() > 1.0e-3, "the 1-4 pair must contribute");
    }

    #[test]
    fn pentane_distinguishes_excluded_scaled_and_full_pairs() {
        let positions = chain_positions(5, 0.04, 12);
        let mut atoms = chain_atoms(&positions);
        let mut ff = chain_force_field(5);
        ff.build_topology();
        let pbox = PeriodicBox::cubic(BOX);
        let total = ff.compute_forces(&mut atoms, &pbox);
        let expected = bonded_energy(&ff, &atoms, &pbox) + reference_nonbonded_energy(&atoms);
        assert!(
            (total - expected).abs() < 1.0e-9 * (1.0 + expected.abs()),
            "total {total} vs reference {expected}"
        );
    }

    #[test]
    fn forces_equal_minus_the_energy_gradient() {
        let positions = chain_positions(5, 0.05, 13);
        let mut atoms = chain_atoms(&positions);
        let mut ff = chain_force_field(5);
        ff.build_topology();
        let pbox = PeriodicBox::cubic(BOX);
        ff.compute_forces(&mut atoms, &pbox);
        let numeric = numeric_forces(&ff, &positions, &pbox);
        for (atom, (analytic, fd)) in atoms.forces.iter().zip(&numeric).enumerate() {
            for axis in 0..3 {
                let (a, n) = (analytic[axis], fd[axis]);
                assert!(
                    (a - n).abs() < 1.0e-5 * (1.0 + n.abs()),
                    "atom {atom} axis {axis}: analytic {a} vs finite difference {n}"
                );
            }
        }
        let largest = atoms.forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
        assert!(
            largest > 1.0,
            "the test geometry must be strained, |F| = {largest}"
        );
    }

    #[test]
    fn net_force_and_torque_vanish() {
        let positions = chain_positions(6, 0.05, 14);
        let mut atoms = chain_atoms(&positions);
        let mut ff = chain_force_field(6);
        ff.build_topology();
        ff.compute_forces(&mut atoms, &PeriodicBox::cubic(BOX));
        let net: Vec3 = atoms.forces.iter().fold(Vec3::zeros(), |acc, f| acc + f);
        assert!(net.norm() < 1.0e-9, "net force {net:?}");
        let torque: Vec3 = atoms
            .positions
            .iter()
            .zip(&atoms.forces)
            .fold(Vec3::zeros(), |acc, (r, f)| acc + r.cross(f));
        assert!(torque.norm() < 1.0e-7, "net torque {torque:?}");
    }

    #[test]
    fn virial_matches_the_dilation_derivative_of_the_energy() {
        let positions = chain_positions(5, 0.05, 15);
        let atoms = chain_atoms(&positions);
        let mut ff = chain_force_field(5);
        ff.build_topology();
        let pbox = PeriodicBox::cubic(BOX);
        let mut scratch = atoms.clone();
        let (energy, virial) = ff.compute_forces_virial(&mut scratch, &pbox);
        let mut reference = atoms.clone();
        assert!((energy - ff.compute_forces(&mut reference, &pbox)).abs() < 1.0e-12);
        let numeric = scaling_virial(&ff, &atoms, &pbox);
        assert!(
            (virial - numeric).abs() < 1.0e-6 * (1.0 + virial.abs()),
            "analytic virial {virial} vs numeric {numeric}"
        );
        assert!(virial.abs() > 1.0e-3, "the virial must be non-trivial");
    }

    #[test]
    fn coulomb_interaction_matches_the_reaction_field_closed_forms() {
        let pbox = PeriodicBox::cubic(BOX);
        let (q_i, q_j) = (0.6, -0.4);
        let rc = 9.0;
        let r = 4.5;

        // Conducting boundary: E = K q q (1/r + r²/(2 rc³) − 3/(2 rc)).
        let mut ff = AmberForceField::new();
        ff.coulomb_cutoff = rc;
        let mut atoms = two_atoms(r, q_i, q_j);
        let energy = ff.compute_forces(&mut atoms, &pbox);
        let expected = COULOMB_K * q_i * q_j * (1.0 / r + r * r / (2.0 * rc.powi(3)) - 1.5 / rc);
        assert!(
            (energy - expected).abs() < 1.0e-9 * expected.abs(),
            "{energy} vs {expected}"
        );
        // Opposite charges attract: the force on atom j points towards atom i.
        let f_mag = COULOMB_K * q_i * q_j * (1.0 / (r * r) - r / rc.powi(3));
        assert!(f_mag < 0.0);
        assert!((atoms.forces[1].x - f_mag).abs() < 1.0e-9 * f_mag.abs());
        assert!((atoms.forces[0].x + f_mag).abs() < 1.0e-9 * f_mag.abs());
        assert!(atoms.forces[1].y.abs() < 1.0e-12 && atoms.forces[1].z.abs() < 1.0e-12);

        // Finite permittivity: the crate's reaction-field model is the reference.
        let mut finite = AmberForceField::new();
        finite.coulomb_cutoff = rc;
        finite.reaction_field_epsilon = 78.5;
        let reference = ReactionFieldElectrostatics::new(78.5, rc);
        let mut atoms = two_atoms(r, q_i, q_j);
        let energy = finite.compute_forces(&mut atoms, &pbox);
        let expected = reference.energy(q_i, q_j, r);
        assert!(
            (energy - expected).abs() < 1.0e-9 * expected.abs(),
            "{energy} vs {expected}"
        );
        let expected_force = reference.force_magnitude(q_i, q_j, r);
        assert!((atoms.forces[1].x - expected_force).abs() < 1.0e-9 * expected_force.abs());

        // Energy (and, for the conducting boundary, force) vanish at the cutoff.
        let mut atoms = two_atoms(rc - 1.0e-6, q_i, q_j);
        let energy = ff.compute_forces(&mut atoms, &pbox);
        assert!(energy.abs() < 1.0e-6, "energy at the cutoff {energy}");
        assert!(atoms.forces[1].norm() < 1.0e-5, "force at the cutoff");
        // Beyond the cutoff there is nothing.
        let mut atoms = two_atoms(rc + 0.5, q_i, q_j);
        assert_eq!(ff.compute_forces(&mut atoms, &pbox), 0.0);
        assert!(atoms.forces.iter().all(|f| f.norm() == 0.0));
    }

    #[test]
    fn coulomb_constant_zero_switches_the_electrostatics_off() {
        let positions = chain_positions(5, 0.04, 16);
        let pbox = PeriodicBox::cubic(BOX);
        let mut charged = chain_force_field(5);
        charged.coulomb_k = 0.0;
        charged.build_topology();
        let mut neutral = charged.clone();
        neutral.coulomb_k = COULOMB_K;
        let mut with_charges = chain_atoms(&positions);
        let mut without_charges = chain_atoms(&positions);
        without_charges.charges.iter_mut().for_each(|q| *q = 0.0);
        let a = charged.compute_forces(&mut with_charges, &pbox);
        let b = neutral.compute_forces(&mut without_charges, &pbox);
        assert!((a - b).abs() < 1.0e-12, "{a} vs {b}");
    }

    #[test]
    fn one_four_scaling_acts_on_the_one_four_pair_only() {
        let positions = chain_positions(5, 0.04, 17);
        let pbox = PeriodicBox::cubic(BOX);
        let energy_with = |scaling: OneFourScaling| {
            let mut ff = chain_force_field(5);
            ff.scaling_14 = scaling;
            ff.build_topology();
            let mut atoms = chain_atoms(&positions);
            ff.compute_forces(&mut atoms, &pbox)
        };
        let atoms = chain_atoms(&positions);
        let lj = LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF);
        let one_four_pairs = [(0, 3), (1, 4)];
        let full_lj: f64 = one_four_pairs
            .iter()
            .map(|&(i, j)| lj.energy(pair_distance(&atoms, i, j)))
            .sum();
        let full_coulomb: f64 = one_four_pairs
            .iter()
            .map(|&(i, j)| {
                reaction_field_energy(
                    atoms.charges[i],
                    atoms.charges[j],
                    pair_distance(&atoms, i, j),
                )
            })
            .sum();
        let off = energy_with(OneFourScaling {
            lj_scale: 0.0,
            coulomb_scale: 0.0,
        });
        let amber = energy_with(OneFourScaling::amber());
        let full = energy_with(OneFourScaling::full());
        assert!(
            (amber - off - (0.5 * full_lj + full_coulomb / 1.2)).abs() < 1.0e-9,
            "AMBER scaling: {amber} - {off}"
        );
        assert!(
            (full - off - (full_lj + full_coulomb)).abs() < 1.0e-9,
            "unscaled 1-4: {full} - {off}"
        );
    }

    #[test]
    fn manual_exclusions_remove_a_pair_and_survive_build_topology() {
        let positions = chain_positions(5, 0.04, 18);
        let pbox = PeriodicBox::cubic(BOX);
        let mut ff = chain_force_field(5);
        ff.exclusions.add_12(0, 4);
        ff.build_topology();
        assert!(ff.exclusions.is_excluded(0, 4));
        let mut atoms = chain_atoms(&positions);
        let total = ff.compute_forces(&mut atoms, &pbox);
        // Remove the (0, 4) pair, the only fully interacting one, from the reference.
        let lj = LennardJones::new(LJ_EPSILON, LJ_SIGMA, CUTOFF);
        let r04 = pair_distance(&atoms, 0, 4);
        let pair_04 =
            lj.energy(r04) + reaction_field_energy(atoms.charges[0], atoms.charges[4], r04);
        let expected =
            bonded_energy(&ff, &atoms, &pbox) + reference_nonbonded_energy(&atoms) - pair_04;
        assert!((total - expected).abs() < 1.0e-9 * (1.0 + expected.abs()));
    }

    #[test]
    fn forgetting_build_topology_gives_the_same_forces() {
        let positions = chain_positions(6, 0.04, 19);
        let pbox = PeriodicBox::cubic(BOX);
        let mut built = chain_force_field(6);
        built.build_topology();
        let unbuilt = chain_force_field(6);
        assert!(!unbuilt.topology_is_current());
        let mut a = chain_atoms(&positions);
        let mut b = chain_atoms(&positions);
        let ea = built.compute_forces(&mut a, &pbox);
        let eb = unbuilt.compute_forces(&mut b, &pbox);
        assert!((ea - eb).abs() < 1.0e-12, "{ea} vs {eb}");
        for (fa, fb) in a.forces.iter().zip(&b.forces) {
            assert!((fa - fb).norm() < 1.0e-12);
        }
    }

    #[test]
    fn molecule_across_the_periodic_boundary_has_the_same_energy_and_forces() {
        let positions = chain_positions(5, 0.04, 20);
        let pbox = PeriodicBox::cubic(BOX);
        let mut ff = chain_force_field(5);
        ff.build_topology();
        let mut centred = chain_atoms(&positions);
        let e_centred = ff.compute_forces(&mut centred, &pbox);

        // Shift by half a box along x and wrap: the chain straddles x = 0 / L.
        let wrapped: Vec<Vec3> = positions
            .iter()
            .map(|p| Vec3::new((p.x + BOX / 2.0).rem_euclid(BOX), p.y, p.z))
            .collect();
        assert!(
            wrapped.iter().any(|p| p.x < 5.0) && wrapped.iter().any(|p| p.x > BOX - 5.0),
            "the shifted molecule must straddle the boundary"
        );
        let mut straddling = chain_atoms(&wrapped);
        let e_straddling = ff.compute_forces(&mut straddling, &pbox);
        assert!(
            (e_centred - e_straddling).abs() < 1.0e-9 * (1.0 + e_centred.abs()),
            "{e_centred} vs {e_straddling}"
        );
        for (a, b) in centred.forces.iter().zip(&straddling.forces) {
            assert!((a - b).norm() < 1.0e-8 * (1.0 + a.norm()), "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn empty_and_unparameterised_systems_are_harmless() {
        let pbox = PeriodicBox::cubic(BOX);
        let ff = AmberForceField::new();
        let mut empty = AtomSet::new();
        assert_eq!(ff.compute_forces(&mut empty, &pbox), 0.0);
        // Charges without any cutoff (no pair potential, no Coulomb cutoff).
        let mut atoms = two_atoms(3.0, 1.0, -1.0);
        assert_eq!(ff.compute_forces(&mut atoms, &pbox), 0.0);
        let mut atoms = two_atoms(3.0, 1.0, -1.0);
        assert_eq!(ff.compute_forces_virial(&mut atoms, &pbox), (0.0, 0.0));
        assert!(ExclusionList::new().is_empty());
    }
}
