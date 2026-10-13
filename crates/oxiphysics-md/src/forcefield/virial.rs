//! Translation-invariant scalar virial of a [`ForceField`].
//!
//! The pressure of a periodic system follows from the virial theorem,
//! `P = (2·KE + W) / (3·V)`, where `W` is the *pair* virial
//! `W = Σ_{pairs} r_ij · F_ij` evaluated with minimum-image separations `r_ij`.
//!
//! The "atomic" virial `Σ_i r_i · F_i` of wrapped coordinates is **not** a valid
//! substitute under periodic boundary conditions: a pair that interacts across
//! a cell face contributes the spurious term `L · F_ij` (the box length times
//! the pair force) because the two wrapped positions differ by a lattice
//! vector.  The pair virial depends only on separations and is therefore
//! invariant under translation of the whole system.
//!
//! [`ForceField::compute_forces_virial`] returns this quantity.  Force fields
//! with an analytic pair decomposition (pair and bond terms) override it;
//! everything else falls back on [`scaling_virial`], which obtains the virial
//! from the energy alone.

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::ForceField;

/// Relative amount `h` by which [`scaling_virial`] dilates the system
/// (`λ = 1 ± h`) for its central finite difference.
///
/// The truncation error of the central difference is `O(h²)` times the third
/// `λ`-derivative of the energy and the round-off error is `O(ε_mach / h)`
/// times the energy scale, so `10⁻⁶` balances both at about `10⁻¹⁰` relative
/// accuracy for smooth interactions.
pub const VIRIAL_SCALING_STEP: f64 = 1.0e-6;

/// Scalar virial `W = Σ r · F` of `ff`, obtained from the energy alone.
///
/// Dilating every coordinate *and* the periodic box by a factor `λ` scales each
/// minimum-image separation by exactly `λ`, so the energy `U(λ)` of the
/// dilated system satisfies
///
/// ```text
/// dU/dλ |_{λ=1} = Σ_{pairs} ∇U · r_ij = −W.
/// ```
///
/// The derivative is evaluated with a central difference of step
/// [`VIRIAL_SCALING_STEP`] on a scratch copy of `atoms` (the input is not
/// modified).  Because only energies enter, the result is valid for **any**
/// translation-invariant force field, including many-body terms (angles,
/// dihedrals), and handles periodic boundaries exactly.  Energies that are
/// invariant under uniform dilation (pure angle and torsion terms) correctly
/// yield a zero virial.
///
/// The cost is two extra energy evaluations.  A potential with a hard cutoff
/// that is not shifted to zero has a discontinuous energy; a pair that
/// crosses the cutoff within the finite-difference window perturbs the result.
/// Force fields in that situation should override
/// [`ForceField::compute_forces_virial`] with an analytic pair sum.
pub fn scaling_virial<F: ForceField + ?Sized>(ff: &F, atoms: &AtomSet, pbox: &PeriodicBox) -> f64 {
    let dilated_energy = |lambda: f64| -> f64 {
        let mut scratch = atoms.clone();
        scratch.clear_forces();
        for position in &mut scratch.positions {
            *position *= lambda;
        }
        let dilated_box = PeriodicBox::new(pbox.dims * lambda);
        ff.compute_forces(&mut scratch, &dilated_box)
    };
    let h = VIRIAL_SCALING_STEP;
    let e_plus = dilated_energy(1.0 + h);
    let e_minus = dilated_energy(1.0 - h);
    -(e_plus - e_minus) / (2.0 * h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{
        AngleForceField, BondedForceField, DihedralForceField, PairForceField,
    };
    use crate::potential::{LennardJones, Potential};
    use oxiphysics_core::math::Vec3;

    fn lj_pair_atoms(x0: f64, x1: f64) -> AtomSet {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(x0, 5.0, 5.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(x1, 5.0, 5.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms
    }

    fn lj_field() -> PairForceField {
        let mut ff = PairForceField::new();
        ff.add_interaction(0, 0, LennardJones::new(1.0, 1.0, 2.5));
        ff
    }

    /// Wrapper that hides the analytic override so the numeric default runs.
    struct NumericOnly<'a, F: ForceField>(&'a F);

    impl<F: ForceField> ForceField for NumericOnly<'_, F> {
        fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
            self.0.compute_forces(atoms, pbox)
        }
    }

    #[test]
    fn scaling_virial_of_lj_pair_matches_f_times_r() {
        let ff = lj_field();
        let pbox = PeriodicBox::cubic(10.0);
        let atoms = lj_pair_atoms(4.6, 5.5);
        let r = 0.9_f64;
        let lj = LennardJones::new(1.0, 1.0, 2.5);
        let expected = lj.force(r) * r;
        let numeric = scaling_virial(&ff, &atoms, &pbox);
        assert!(
            (numeric - expected).abs() <= 1e-6 * expected.abs(),
            "numeric virial {numeric} vs analytic f*r {expected}"
        );
    }

    #[test]
    fn scaling_virial_does_not_modify_the_input() {
        let ff = lj_field();
        let pbox = PeriodicBox::cubic(10.0);
        let mut atoms = lj_pair_atoms(4.6, 5.5);
        let _ = ff.compute_forces(&mut atoms, &pbox);
        let positions_before = atoms.positions.clone();
        let forces_before = atoms.forces.clone();
        let _ = scaling_virial(&ff, &atoms, &pbox);
        assert_eq!(atoms.positions, positions_before);
        assert_eq!(atoms.forces, forces_before);
    }

    #[test]
    fn default_compute_forces_virial_agrees_with_pair_override() {
        let ff = lj_field();
        let pbox = PeriodicBox::cubic(10.0);
        let mut atoms = AtomSet::new();
        let coordinates = [
            [1.0, 1.0, 1.0],
            [1.9, 1.2, 0.8],
            [2.6, 2.1, 1.9],
            [9.4, 9.8, 1.3],
            [0.3, 0.1, 1.1],
            [5.0, 5.0, 5.0],
            [5.8, 5.3, 4.6],
        ];
        for c in coordinates {
            atoms.add_atom(Vec3::new(c[0], c[1], c[2]), Vec3::zeros(), 1.0, 0.0, 0);
        }
        let mut analytic_atoms = atoms.clone();
        let (e_analytic, w_analytic) = ff.compute_forces_virial(&mut analytic_atoms, &pbox);
        let mut numeric_atoms = atoms.clone();
        let (e_numeric, w_numeric) =
            NumericOnly(&ff).compute_forces_virial(&mut numeric_atoms, &pbox);
        assert!(
            (e_analytic - e_numeric).abs() < 1e-12,
            "energies {e_analytic} vs {e_numeric}"
        );
        assert!(
            (w_analytic - w_numeric).abs() <= 1e-6 * w_analytic.abs().max(1.0),
            "analytic virial {w_analytic} vs numeric {w_numeric}"
        );
    }

    #[test]
    fn bonded_virial_override_matches_numeric_default() {
        let mut bonded = BondedForceField::new();
        bonded.add_bond(0, 1, 100.0, 1.0);
        bonded.add_bond(1, 2, 80.0, 1.2);
        let pbox = PeriodicBox::cubic(6.0);
        let mut atoms = AtomSet::new();
        // The molecule straddles the periodic boundary in x.
        for c in [[5.6, 3.0, 3.0], [0.4, 3.3, 3.1], [1.5, 2.9, 3.6]] {
            atoms.add_atom(Vec3::new(c[0], c[1], c[2]), Vec3::zeros(), 1.0, 0.0, 0);
        }
        let mut a = atoms.clone();
        let (e_a, w_a) = bonded.compute_forces_virial(&mut a, &pbox);
        let mut b = atoms.clone();
        let (e_b, w_b) = NumericOnly(&bonded).compute_forces_virial(&mut b, &pbox);
        assert!((e_a - e_b).abs() < 1e-12);
        assert!(
            (w_a - w_b).abs() <= 1e-6 * w_a.abs().max(1.0),
            "analytic bond virial {w_a} vs numeric {w_b}"
        );
        assert!(
            w_a.abs() > 1e-3,
            "test geometry must give a non-trivial virial"
        );
    }

    #[test]
    fn angle_and_dihedral_terms_have_zero_virial() {
        // Angle and torsion energies depend only on angles, which are
        // invariant under uniform dilation: W = 0.
        let pbox = PeriodicBox::cubic(20.0);
        let mut atoms = AtomSet::new();
        for c in [
            [10.0, 10.0, 10.0],
            [11.2, 10.4, 10.1],
            [11.9, 11.5, 10.3],
            [13.0, 11.7, 11.0],
        ] {
            atoms.add_atom(Vec3::new(c[0], c[1], c[2]), Vec3::zeros(), 1.0, 0.0, 0);
        }
        let mut angles = AngleForceField::new();
        angles.add_angle(0, 1, 2, 50.0, 1.9);
        let mut dihedrals = DihedralForceField::new();
        dihedrals.add_dihedral(0, 1, 2, 3, 2.0, 0.3, 2);
        let mut a = atoms.clone();
        let (e_angle, w_angle) = angles.compute_forces_virial(&mut a, &pbox);
        assert!(e_angle > 0.0);
        assert!(w_angle.abs() < 1e-6, "angle virial {w_angle}");
        let mut d = atoms.clone();
        let (_, w_dihedral) = dihedrals.compute_forces_virial(&mut d, &pbox);
        assert!(w_dihedral.abs() < 1e-6, "dihedral virial {w_dihedral}");
    }
}
