//! # AngleForceField - Trait Implementations
//!
//! This module contains trait implementations for `AngleForceField`.
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
use super::types::AngleForceField;

impl Default for AngleForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for AngleForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        for angle in &self.angles {
            let (rji, dji) =
                distance_pbc(&atoms.positions[angle.j], &atoms.positions[angle.i], pbox);
            let (rjk, djk) =
                distance_pbc(&atoms.positions[angle.j], &atoms.positions[angle.k], pbox);
            if dji < 1e-15 || djk < 1e-15 {
                continue;
            }
            let rji_n = rji / dji;
            let rjk_n = rjk / djk;
            let cos_theta = rji_n.dot(&rjk_n).clamp(-1.0, 1.0);
            let theta = cos_theta.acos();
            let delta = theta - angle.theta0;
            energy += 0.5 * angle.k_theta * delta * delta;
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-15);
            // dE/dθ for E = ½ k_θ (θ − θ₀)².
            let d_energy = angle.k_theta * delta;
            // With a = r_i − r_j and b = r_k − r_j, cos θ = â·b̂ and
            //   ∂θ/∂r_i = −(b̂ − cos θ â) / (|a| sin θ),
            // so the force F_i = −dE/dr_i = −(dE/dθ) ∂θ/∂r_i points along
            // +(b̂ − cos θ â): it moves atom i towards the side of atom k when
            // θ > θ₀ and away from it when θ < θ₀, restoring θ₀.  F_k follows
            // by symmetry and F_j balances the two (translation invariance).
            let fi = (rjk_n - rji_n * cos_theta) / (dji * sin_theta) * d_energy;
            let fk = (rji_n - rjk_n * cos_theta) / (djk * sin_theta) * d_energy;
            atoms.forces[angle.i] += fi;
            atoms.forces[angle.k] += fk;
            atoms.forces[angle.j] -= fi + fk;
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
    fn energy_at(ff: &AngleForceField, positions: &[Vec3], pbox: &PeriodicBox) -> f64 {
        let mut atoms = atom_set(positions);
        ff.compute_forces(&mut atoms, pbox)
    }

    /// Analytic forces accumulated by `compute_forces` (starting from zero).
    fn analytic_forces(ff: &AngleForceField, positions: &[Vec3], pbox: &PeriodicBox) -> Vec<Vec3> {
        let mut atoms = atom_set(positions);
        ff.compute_forces(&mut atoms, pbox);
        atoms.forces
    }

    /// `-dE/dx` by central differences for every Cartesian component.
    fn numeric_forces(ff: &AngleForceField, positions: &[Vec3], pbox: &PeriodicBox) -> Vec<Vec3> {
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

    /// Bond angle (radians) at the central atom `j` of the triple `(i, j, k)`.
    fn angle_of(positions: &[Vec3], pbox: &PeriodicBox) -> f64 {
        let (rji, dji) = distance_pbc(&positions[1], &positions[0], pbox);
        let (rjk, djk) = distance_pbc(&positions[1], &positions[2], pbox);
        (rji.dot(&rjk) / (dji * djk)).clamp(-1.0, 1.0).acos()
    }

    /// Three non-collinear atoms around a central atom inside a large box.
    fn random_triple(rng: &mut StdRng, centre: Vec3) -> Vec<Vec3> {
        loop {
            let mut arm = || {
                Vec3::new(
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                )
            };
            let a = arm();
            let b = arm();
            if a.norm() < 0.4 || b.norm() < 0.4 {
                continue;
            }
            let sin_theta = (1.0 - (a.dot(&b) / (a.norm() * b.norm())).powi(2)).sqrt();
            if sin_theta > 0.3 {
                return vec![centre + a * 1.6, centre, centre + b * 1.4];
            }
        }
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

    #[test]
    fn angle_forces_equal_minus_gradient_above_and_below_equilibrium() {
        let pbox = PeriodicBox::cubic(60.0);
        let mut rng = StdRng::seed_from_u64(3448);
        for trial in 0..20 {
            let positions = random_triple(&mut rng, Vec3::new(30.0, 30.0, 30.0));
            let theta = angle_of(&positions, &pbox);
            // theta0 below theta: the angle is too wide (delta > 0); above: too narrow.
            for (label, theta0) in [("wide", theta - 0.35), ("narrow", theta + 0.35)] {
                let mut ff = AngleForceField::new();
                ff.add_angle(0, 1, 2, 85.0, theta0);
                let analytic = analytic_forces(&ff, &positions, &pbox);
                let numeric = numeric_forces(&ff, &positions, &pbox);
                assert_forces_match(&analytic, &numeric, &format!("trial {trial} ({label})"));
            }
        }
    }

    #[test]
    fn angle_forces_sum_to_zero_and_are_a_descent_direction() {
        let pbox = PeriodicBox::cubic(60.0);
        let mut rng = StdRng::seed_from_u64(77);
        let positions = random_triple(&mut rng, Vec3::new(30.0, 30.0, 30.0));
        let theta = angle_of(&positions, &pbox);
        let mut ff = AngleForceField::new();
        ff.add_angle(0, 1, 2, 120.0, theta - 0.5);
        let forces = analytic_forces(&ff, &positions, &pbox);
        let net = forces[0] + forces[1] + forces[2];
        assert!(net.norm() < 1.0e-12, "net force {net:?}");

        // A small step along the force must lower the energy; with the reversed
        // sign (the original bug) it raised it.
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
    fn angle_restoring_force_closes_a_wide_angle_and_opens_a_narrow_one() {
        let pbox = PeriodicBox::cubic(40.0);
        // 90-degree geometry: i on +x, k on +y of the central atom j.
        let positions = vec![
            Vec3::new(21.0, 20.0, 20.0),
            Vec3::new(20.0, 20.0, 20.0),
            Vec3::new(20.0, 21.0, 20.0),
        ];
        // Equilibrium 60 degrees < 90: the angle is too wide, so atom i is pushed
        // towards +y (towards atom k) and atom k towards +x.
        let mut wide = AngleForceField::new();
        wide.add_angle(0, 1, 2, 50.0, std::f64::consts::FRAC_PI_3);
        let f = analytic_forces(&wide, &positions, &pbox);
        assert!(f[0].y > 0.0 && f[2].x > 0.0, "{:?} {:?}", f[0], f[2]);
        // |F| = k_theta * |theta - theta0| / r for a right angle with unit arms.
        let expected = 50.0 * (std::f64::consts::FRAC_PI_2 - std::f64::consts::FRAC_PI_3);
        assert!((f[0].y - expected).abs() < 1.0e-9, "{}", f[0].y);

        // Equilibrium 120 degrees > 90: the angle is too narrow, forces reverse.
        let mut narrow = AngleForceField::new();
        narrow.add_angle(0, 1, 2, 50.0, 2.0 * std::f64::consts::FRAC_PI_3);
        let g = analytic_forces(&narrow, &positions, &pbox);
        assert!(g[0].y < 0.0 && g[2].x < 0.0, "{:?} {:?}", g[0], g[2]);
    }

    #[test]
    fn angle_forces_are_consistent_across_the_periodic_boundary() {
        // The triple straddles the x = 0 / x = L face of a small box, so the
        // minimum-image vectors differ from the raw coordinate differences.
        let pbox = PeriodicBox::cubic(12.0);
        let positions = vec![
            Vec3::new(0.4, 5.0, 5.0),
            Vec3::new(11.2, 5.3, 5.1),
            Vec3::new(11.5, 6.6, 4.4),
        ];
        let theta = angle_of(&positions, &pbox);
        let mut ff = AngleForceField::new();
        ff.add_angle(0, 1, 2, 70.0, theta + 0.4);
        let analytic = analytic_forces(&ff, &positions, &pbox);
        let numeric = numeric_forces(&ff, &positions, &pbox);
        assert_forces_match(&analytic, &numeric, "periodic");
    }
}
