//! # DihedralForceField - Trait Implementations
//!
//! This module contains trait implementations for `DihedralForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::{PeriodicBox, distance_pbc};

use super::functions::ForceField;
use super::types::DihedralForceField;

impl Default for DihedralForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for DihedralForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        for dih in &self.dihedrals {
            let (r_ij, _) = distance_pbc(&atoms.positions[dih.i], &atoms.positions[dih.j], pbox);
            let (r_jk, _) = distance_pbc(&atoms.positions[dih.j], &atoms.positions[dih.k], pbox);
            let (r_kl, _) = distance_pbc(&atoms.positions[dih.k], &atoms.positions[dih.l], pbox);
            let n1 = r_ij.cross(&r_jk);
            let n2 = r_jk.cross(&r_kl);
            let n1_norm = n1.norm();
            let n2_norm = n2.norm();
            if n1_norm < 1e-15 || n2_norm < 1e-15 {
                continue;
            }
            let n1_hat = n1 / n1_norm;
            let n2_hat = n2 / n2_norm;
            let cos_phi = n1_hat.dot(&n2_hat).clamp(-1.0, 1.0);
            let m1 = n1_hat.cross(&n2_hat);
            let rjk_hat = r_jk / r_jk.norm().max(1e-15);
            let sin_phi = m1.dot(&rjk_hat);
            let phi = sin_phi.atan2(cos_phi);
            let n_f = dih.n as f64;
            energy += dih.v_n * (1.0 + (n_f * phi - dih.gamma).cos());
            let d_energy_dphi = -dih.v_n * n_f * (n_f * phi - dih.gamma).sin();
            // Blondel-Karplus gradient of the signed torsion with b1 = r_ij,
            // b2 = r_jk, b3 = r_kl, L = |b2|, n1 = b1 x b2, n2 = b2 x b3,
            // x = b1.b2 / L^2 and y = b3.b2 / L^2:
            //   dphi/dr_i = -L n1 / |n1|^2,
            //   dphi/dr_l = +L n2 / |n2|^2,
            //   dphi/dr_j = -(1 + x) dphi/dr_i + y dphi/dr_l,
            //   dphi/dr_k = -(dphi/dr_i + dphi/dr_j + dphi/dr_l).
            // The forces F = -(dE/dphi) dphi/dr are therefore
            //   F_i = +(dE/dphi) L n1 / |n1|^2,
            //   F_l = -(dE/dphi) L n2 / |n2|^2,
            //   F_j = -(1 + x) F_i + y F_l,
            //   F_k = -(F_i + F_j + F_l)   (translation invariance).
            let rjk_len = r_jk.norm().max(1e-15);
            let fi = n1 * (d_energy_dphi * rjk_len / (n1_norm * n1_norm));
            let fl = n2 * (-d_energy_dphi * rjk_len / (n2_norm * n2_norm));
            let x_proj = r_ij.dot(&r_jk) / (rjk_len * rjk_len);
            let y_proj = r_kl.dot(&r_jk) / (rjk_len * rjk_len);
            let fj = fi * (-1.0 - x_proj) + fl * y_proj;
            let fk = -(fi + fj + fl);
            atoms.forces[dih.i] += fi;
            atoms.forces[dih.j] += fj;
            atoms.forces[dih.k] += fk;
            atoms.forces[dih.l] += fl;
        }
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::math::Vec3;
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    /// Central-difference step used for the numerical gradient.
    const FD_STEP: f64 = 1.0e-6;

    fn atom_set(positions: &[Vec3]) -> AtomSet {
        let mut atoms = AtomSet::new();
        for p in positions {
            atoms.add_atom(*p, Vec3::zeros(), 12.0, 0.0, 0);
        }
        atoms
    }

    /// Potential energy returned by `compute_forces` at `positions`.
    fn energy_at(ff: &DihedralForceField, positions: &[Vec3], pbox: &PeriodicBox) -> f64 {
        let mut atoms = atom_set(positions);
        ff.compute_forces(&mut atoms, pbox)
    }

    /// Analytic forces accumulated by `compute_forces` (starting from zero).
    fn analytic_forces(
        ff: &DihedralForceField,
        positions: &[Vec3],
        pbox: &PeriodicBox,
    ) -> Vec<Vec3> {
        let mut atoms = atom_set(positions);
        ff.compute_forces(&mut atoms, pbox);
        atoms.forces
    }

    /// `-dE/dx` by central differences for every Cartesian component.
    fn numeric_forces(
        ff: &DihedralForceField,
        positions: &[Vec3],
        pbox: &PeriodicBox,
    ) -> Vec<Vec3> {
        let mut out = vec![Vec3::zeros(); positions.len()];
        for atom in 0..positions.len() {
            for axis in 0..3 {
                let mut plus = positions.to_vec();
                let mut minus = positions.to_vec();
                plus[atom][axis] += FD_STEP;
                minus[atom][axis] -= FD_STEP;
                let e_plus = energy_at(ff, &plus, pbox);
                let e_minus = energy_at(ff, &minus, pbox);
                out[atom][axis] = -(e_plus - e_minus) / (2.0 * FD_STEP);
            }
        }
        out
    }

    fn assert_forces_match(analytic: &[Vec3], numeric: &[Vec3], context: &str) {
        for (atom, (f_ana, f_num)) in analytic.iter().zip(numeric).enumerate() {
            for axis in 0..3 {
                let a = f_ana[axis];
                let n = f_num[axis];
                assert!(
                    (a - n).abs() < 1.0e-6 * (1.0 + n.abs()),
                    "{context}: atom {atom} axis {axis}: analytic {a} vs finite difference {n}"
                );
            }
        }
    }

    /// Four non-planar atoms i-j-k-l (a bent chain with a general torsion).
    fn random_chain(rng: &mut StdRng, origin: Vec3) -> Vec<Vec3> {
        loop {
            let mut unit = || {
                let v = Vec3::new(
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                );
                if v.norm() < 0.3 {
                    None
                } else {
                    Some(v / v.norm())
                }
            };
            let (Some(a), Some(b), Some(c)) = (unit(), unit(), unit()) else {
                continue;
            };
            let rj = origin;
            let rk = rj + b * 1.5;
            let ri = rj + a * 1.4;
            let rl = rk + c * 1.6;
            let n1 = (rj - ri).cross(&(rk - rj));
            let n2 = (rk - rj).cross(&(rl - rk));
            // Reject near-collinear triples, where the torsion is ill defined.
            if n1.norm() > 0.5 && n2.norm() > 0.5 {
                return vec![ri, rj, rk, rl];
            }
        }
    }

    #[test]
    fn dihedral_forces_equal_minus_gradient_for_all_multiplicities_and_phases() {
        let pbox = PeriodicBox::cubic(60.0);
        let mut rng = StdRng::seed_from_u64(3449);
        for trial in 0..8 {
            let positions = random_chain(&mut rng, Vec3::new(30.0, 30.0, 30.0));
            for n in 1..=3_u32 {
                for gamma in [0.0, std::f64::consts::FRAC_PI_3] {
                    let mut ff = DihedralForceField::new();
                    ff.add_dihedral(0, 1, 2, 3, 6.3, gamma, n);
                    let analytic = analytic_forces(&ff, &positions, &pbox);
                    let numeric = numeric_forces(&ff, &positions, &pbox);
                    assert_forces_match(
                        &analytic,
                        &numeric,
                        &format!("trial {trial}, n = {n}, gamma = {gamma}"),
                    );
                }
            }
        }
    }

    #[test]
    fn dihedral_forces_have_zero_net_force_and_torque() {
        let pbox = PeriodicBox::cubic(60.0);
        let mut rng = StdRng::seed_from_u64(91);
        let positions = random_chain(&mut rng, Vec3::new(30.0, 30.0, 30.0));
        let mut ff = DihedralForceField::new();
        ff.add_dihedral(0, 1, 2, 3, 4.0, 0.7, 2);
        let forces = analytic_forces(&ff, &positions, &pbox);
        let net_force: Vec3 = forces.iter().fold(Vec3::zeros(), |acc, f| acc + f);
        assert!(net_force.norm() < 1.0e-10, "net force {net_force:?}");
        // The energy is rotation invariant, so the net torque about any point vanishes.
        let net_torque: Vec3 = positions
            .iter()
            .zip(&forces)
            .fold(Vec3::zeros(), |acc, (r, f)| acc + r.cross(f));
        assert!(net_torque.norm() < 1.0e-9, "net torque {net_torque:?}");
    }

    #[test]
    fn dihedral_forces_are_a_descent_direction() {
        let pbox = PeriodicBox::cubic(60.0);
        let mut rng = StdRng::seed_from_u64(5);
        let positions = random_chain(&mut rng, Vec3::new(30.0, 30.0, 30.0));
        let mut ff = DihedralForceField::new();
        ff.add_dihedral(0, 1, 2, 3, 5.0, 0.0, 3);
        let forces = analytic_forces(&ff, &positions, &pbox);
        let e0 = energy_at(&ff, &positions, &pbox);
        let stepped: Vec<Vec3> = positions
            .iter()
            .zip(&forces)
            .map(|(p, f)| p + f * 1.0e-5)
            .collect();
        let e1 = energy_at(&ff, &stepped, &pbox);
        assert!(e1 < e0, "energy rose along the force: {e0} -> {e1}");
    }

    #[test]
    fn dihedral_planar_trans_minimum_has_zero_force_and_a_restoring_tilt_force() {
        // Planar zig-zag (trans) chain: phi = pi.  With V (1 + cos(phi)), n = 1 and
        // gamma = 0 this is the energy minimum, so every force vanishes.
        let pbox = PeriodicBox::cubic(30.0);
        let positions = vec![
            Vec3::new(10.0, 11.0, 10.0),
            Vec3::new(10.0, 10.0, 10.0),
            Vec3::new(11.5, 10.0, 10.0),
            Vec3::new(11.5, 9.0, 10.0),
        ];
        let mut ff = DihedralForceField::new();
        ff.add_dihedral(0, 1, 2, 3, 3.0, 0.0, 1);
        assert!(energy_at(&ff, &positions, &pbox).abs() < 1.0e-12);
        for f in &analytic_forces(&ff, &positions, &pbox) {
            assert!(f.norm() < 1.0e-9, "force at the torsional minimum: {f:?}");
        }

        // Lifting the end atom out of the plane raises the energy, so the force on
        // it must point back towards the plane (-z), and it must equal -grad E.
        let mut tilted = positions.clone();
        tilted[3].z += 0.3;
        let analytic = analytic_forces(&ff, &tilted, &pbox);
        let numeric = numeric_forces(&ff, &tilted, &pbox);
        assert_forces_match(&analytic, &numeric, "tilted trans");
        assert!(
            analytic[3].z < 0.0,
            "restoring z-force on atom l: {:?}",
            analytic[3]
        );
    }

    #[test]
    fn dihedral_forces_are_consistent_across_the_periodic_boundary() {
        // The chain straddles the x = 0 / x = L face of a small box.
        let pbox = PeriodicBox::cubic(12.0);
        let positions = vec![
            Vec3::new(0.3, 5.0, 5.0),
            Vec3::new(11.3, 5.6, 5.2),
            Vec3::new(11.9, 6.9, 4.6),
            Vec3::new(0.9, 7.4, 5.8),
        ];
        let mut ff = DihedralForceField::new();
        ff.add_dihedral(0, 1, 2, 3, 2.5, std::f64::consts::FRAC_PI_3, 2);
        let analytic = analytic_forces(&ff, &positions, &pbox);
        let numeric = numeric_forces(&ff, &positions, &pbox);
        assert_forces_match(&analytic, &numeric, "periodic");
    }
}
