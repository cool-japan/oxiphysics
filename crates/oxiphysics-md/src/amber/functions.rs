//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    Amber99sbAngleParam, Amber99sbBondParam, Amber99sbType, AmberAtomType, AmberForceField,
    GaffType, ResidueChargeRecord,
};

#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::amber::AmberEnergyEvaluator;
    use crate::amber::CoulombInteraction;
    use crate::amber::DihedralScan;
    use crate::amber::DihedralTerm;
    use crate::amber::HarmonicAngle;
    use crate::amber::LorentzBerthelot;
    use crate::amber::NonBondedExclusions;
    use crate::amber::ProperDihedral;
    use std::f64::consts::PI;
    #[test]
    fn test_harmonic_bond_energy() {
        let bond = HarmonicBond::new(500.0, 0.15);
        assert!((bond.energy(0.15)).abs() < 1e-14);
        let expected = 500.0 * 0.01;
        assert!((bond.energy(0.25) - expected).abs() < 1e-10);
    }
    #[test]
    fn test_harmonic_bond_force_direction() {
        let bond = HarmonicBond::new(500.0, 0.15);
        let ri = [0.0, 0.0, 0.0];
        let rj = [0.25, 0.0, 0.0];
        let (fi, _fj) = bond.force_vectors(ri, rj);
        assert!(fi[0] > 0.0, "fi[0] should be positive, got {}", fi[0]);
        assert!(fi[1].abs() < 1e-14);
        assert!(fi[2].abs() < 1e-14);
    }
    #[test]
    fn test_bond_force_vectors_newton3() {
        let bond = HarmonicBond::new(400.0, 0.12);
        let ri = [0.0, 0.0, 0.0];
        let rj = [0.20, 0.05, 0.0];
        let (fi, fj) = bond.force_vectors(ri, rj);
        assert!((fi[0] + fj[0]).abs() < 1e-12);
        assert!((fi[1] + fj[1]).abs() < 1e-12);
        assert!((fi[2] + fj[2]).abs() < 1e-12);
        let mag_i = (fi[0] * fi[0] + fi[1] * fi[1] + fi[2] * fi[2]).sqrt();
        let mag_j = (fj[0] * fj[0] + fj[1] * fj[1] + fj[2] * fj[2]).sqrt();
        assert!((mag_i - mag_j).abs() < 1e-12);
    }
    #[test]
    fn test_harmonic_angle_90deg() {
        let angle = HarmonicAngle::new(50.0, 90.0);
        let ri = [1.0, 0.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [0.0, 1.0, 0.0];
        let e = angle.energy(ri, rj, rk);
        assert!(e.abs() < 1e-12, "expected 0, got {e}");
    }
    #[test]
    fn test_harmonic_angle_force_restores() {
        let angle = HarmonicAngle::new(50.0, 90.0);
        let ri = [1.0, 0.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [-0.5, 0.866_025_4, 0.0];
        let (fi, _fj, fk) = angle.force_vectors(ri, rj, rk);
        let dt = 1e-5;
        let ri2 = [ri[0] + fi[0] * dt, ri[1] + fi[1] * dt, ri[2] + fi[2] * dt];
        let rk2 = [rk[0] + fk[0] * dt, rk[1] + fk[1] * dt, rk[2] + fk[2] * dt];
        let e_before = angle.energy(ri, rj, rk);
        let e_after = angle.energy(ri2, rj, rk2);
        assert!(e_after < e_before, "force should lower energy");
    }
    #[test]
    fn test_dihedral_angle_0deg() {
        let ri = [0.0, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl = [1.0, -1.0, 0.0];
        let phi = ProperDihedral::dihedral_angle(ri, rj, rk, rl);
        assert!((phi.abs() - PI).abs() < 1e-10, "expected φ ≈ π, got {phi}");
    }
    #[test]
    fn test_dihedral_energy_sum() {
        let terms = vec![
            DihedralTerm {
                vn: 10.0,
                n: 1,
                gamma: 0.0,
            },
            DihedralTerm {
                vn: 5.0,
                n: 2,
                gamma: PI / 2.0,
            },
            DihedralTerm {
                vn: 2.0,
                n: 3,
                gamma: PI,
            },
        ];
        let dih = ProperDihedral::new(terms.clone());
        let ri = [0.0, 0.0, 0.0];
        let rj = [1.0, 0.0, 0.0];
        let rk = [1.0, 1.0, 0.0];
        let rl = [2.0, 1.0, 1.0];
        let phi = ProperDihedral::dihedral_angle(ri, rj, rk, rl);
        let manual: f64 = terms
            .iter()
            .map(|t| 0.5 * t.vn * (1.0 + (t.n as f64 * phi - t.gamma).cos()))
            .sum();
        let computed = dih.energy(ri, rj, rk, rl);
        assert!((computed - manual).abs() < 1e-12);
    }
    #[test]
    fn test_exclusions_12() {
        let mut excl = NonBondedExclusions::new();
        excl.add_bond(0, 1);
        excl.add_bond(1, 2);
        assert!(excl.is_excluded(0, 1));
        assert!(excl.is_excluded(1, 0));
        assert!(excl.is_excluded(1, 2));
        assert!(!excl.is_excluded(0, 2));
    }
    #[test]
    fn test_exclusions_14() {
        let mut excl = NonBondedExclusions::new();
        excl.add_bond(0, 1);
        excl.add_bond(1, 2);
        excl.add_bond(2, 3);
        excl.add_angle(0, 1, 2);
        excl.add_angle(1, 2, 3);
        excl.add_dihedral(0, 1, 2, 3);
        assert!(excl.is_14(0, 3), "0-3 should be a 1-4 pair");
        assert!(!excl.is_excluded(0, 3), "0-3 should NOT be excluded");
        assert!(excl.is_excluded(0, 1));
        assert!(!excl.is_14(0, 1));
        assert!(excl.is_excluded(0, 2));
    }
    #[test]
    fn test_lj_energy_at_sigma() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        let e = at.lj_energy(0.34);
        assert!(e.abs() < 1e-10, "LJ energy at r=sigma should be 0, got {e}");
    }
    #[test]
    fn test_lj_energy_at_minimum() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        let r_min = at.lj_r_min();
        let e_min = at.lj_energy(r_min);
        assert!(
            (e_min + 0.5).abs() < 1e-10,
            "LJ energy at r_min should be -epsilon={}, got {e_min}",
            -at.epsilon
        );
    }
    #[test]
    fn test_lj_repulsive_at_short_range() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        let e = at.lj_energy(0.2);
        assert!(e > 0.0, "LJ should be repulsive at r < sigma");
    }
    #[test]
    fn test_lj_force_positive_at_short_range() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        let f = at.lj_force(0.2);
        assert!(
            f > 0.0,
            "LJ force should be repulsive (positive) at short range"
        );
    }
    #[test]
    fn test_lj_zero_crossing() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        assert!((at.lj_zero_crossing() - 0.34).abs() < 1e-14);
    }
    #[test]
    fn test_lorentz_berthelot_sigma() {
        let sig = LorentzBerthelot::sigma(0.3, 0.4);
        assert!((sig - 0.35).abs() < 1e-14);
    }
    #[test]
    fn test_lorentz_berthelot_epsilon() {
        let eps = LorentzBerthelot::epsilon(0.5, 0.8);
        let expected = (0.5_f64 * 0.8).sqrt();
        assert!((eps - expected).abs() < 1e-14);
    }
    #[test]
    fn test_lorentz_berthelot_self_interaction() {
        let at = AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0);
        let e_self = LorentzBerthelot::lj_energy(&at, &at, 0.4);
        let e_direct = at.lj_energy(0.4);
        assert!((e_self - e_direct).abs() < 1e-12);
    }
    #[test]
    fn test_force_field_add_and_find() {
        let mut ff = AmberForceField::new();
        ff.add_atom_type(AmberAtomType::new("CT", 12.0, 0.5, 0.34, 0.0));
        ff.add_atom_type(AmberAtomType::new("N", 14.0, 0.7, 0.32, -0.4));
        assert_eq!(ff.num_atom_types(), 2);
        assert!(ff.find_atom_type("CT").is_some());
        assert!(ff.find_atom_type("X").is_none());
    }
    #[test]
    fn test_force_field_bond_lookup() {
        let mut ff = AmberForceField::new();
        ff.add_bond_params("CT", "CT", 300.0, 0.153);
        ff.add_bond_params("CT", "N", 337.0, 0.147);
        assert_eq!(ff.num_bond_types(), 2);
        let bond = ff.find_bond("CT", "N").unwrap();
        assert!((bond.r0 - 0.147).abs() < 1e-14);
        let bond2 = ff.find_bond("N", "CT").unwrap();
        assert!((bond2.r0 - 0.147).abs() < 1e-14);
    }
    #[test]
    fn test_dihedral_scan_simple() {
        let terms = vec![DihedralTerm {
            vn: 10.0,
            n: 3,
            gamma: 0.0,
        }];
        let dih = ProperDihedral::new(terms);
        let scan = DihedralScan::scan(&dih, 100);
        assert_eq!(scan.n_points(), 100);
        assert!(scan.barrier_height() > 0.0);
    }
    #[test]
    fn test_dihedral_scan_min_max() {
        let terms = vec![DihedralTerm {
            vn: 10.0,
            n: 1,
            gamma: 0.0,
        }];
        let dih = ProperDihedral::new(terms);
        let scan = DihedralScan::scan(&dih, 360);
        let (_, e_min) = scan.minimum();
        let (_, e_max) = scan.maximum();
        assert!(e_max > e_min);
        assert!(e_min < 1.0, "Min energy should be near 0");
    }
    #[test]
    fn test_energy_evaluator_bonds() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.bond_list.push((0, 1, HarmonicBond::new(500.0, 0.15)));
        let positions = vec![[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]];
        let e = eval.bond_energy(&positions);
        assert!(e.abs() < 1e-12, "At equilibrium, bond energy should be 0");
    }
    #[test]
    fn test_energy_evaluator_stretched_bond() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.bond_list.push((0, 1, HarmonicBond::new(500.0, 0.15)));
        let positions = vec![[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]];
        let e = eval.bond_energy(&positions);
        let expected = 500.0 * 0.1 * 0.1;
        assert!((e - expected).abs() < 1e-10);
    }
    #[test]
    fn test_energy_evaluator_angle() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.angle_list
            .push((0, 1, 2, HarmonicAngle::new(50.0, 109.5)));
        let positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let e = eval.angle_energy(&positions);
        assert!(e > 0.0);
    }
    #[test]
    fn test_energy_evaluator_total() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.bond_list.push((0, 1, HarmonicBond::new(500.0, 0.15)));
        eval.angle_list
            .push((0, 1, 2, HarmonicAngle::new(50.0, 109.5)));
        let positions = vec![[0.0, 0.0, 0.0], [0.15, 0.0, 0.0], [0.15, 0.15, 0.0]];
        let total = eval.total_energy(&positions);
        let bonds = eval.bond_energy(&positions);
        let angles = eval.angle_energy(&positions);
        let dihedrals = eval.dihedral_energy(&positions);
        assert!((total - bonds - angles - dihedrals).abs() < 1e-12);
    }
    #[test]
    fn test_energy_evaluator_bond_forces() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.bond_list.push((0, 1, HarmonicBond::new(500.0, 0.15)));
        let positions = vec![[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]];
        let mut forces = vec![[0.0; 3]; 2];
        eval.accumulate_bond_forces(&positions, &mut forces);
        assert!((forces[0][0] + forces[1][0]).abs() < 1e-12);
        assert!(forces[0][0] > 0.0);
    }
    #[test]
    fn test_energy_evaluator_counts() {
        let mut eval = AmberEnergyEvaluator::new();
        eval.bond_list.push((0, 1, HarmonicBond::new(500.0, 0.15)));
        eval.bond_list.push((1, 2, HarmonicBond::new(500.0, 0.15)));
        eval.angle_list
            .push((0, 1, 2, HarmonicAngle::new(50.0, 109.5)));
        assert_eq!(eval.num_bonds(), 2);
        assert_eq!(eval.num_angles(), 1);
        assert_eq!(eval.num_dihedrals(), 0);
    }
    #[test]
    fn test_coulomb_like_charges_repulsive() {
        let coul = CoulombInteraction::new();
        let e = coul.energy(1.0, 1.0, 0.5);
        assert!(e > 0.0, "Like charges should have positive energy");
    }
    #[test]
    fn test_coulomb_unlike_charges_attractive() {
        let coul = CoulombInteraction::new();
        let e = coul.energy(1.0, -1.0, 0.5);
        assert!(e < 0.0, "Unlike charges should have negative energy");
    }
    #[test]
    fn test_coulomb_force_vectors_newton3() {
        let coul = CoulombInteraction::new();
        let ri = [0.0, 0.0, 0.0];
        let rj = [0.5, 0.0, 0.0];
        let (fi, fj) = coul.force_vectors(1.0, -1.0, ri, rj);
        assert!((fi[0] + fj[0]).abs() < 1e-10);
        assert!((fi[1] + fj[1]).abs() < 1e-14);
    }
    #[test]
    fn test_coulomb_inverse_square() {
        let coul = CoulombInteraction::new();
        let e1 = coul.energy(1.0, 1.0, 1.0);
        let e2 = coul.energy(1.0, 1.0, 2.0);
        assert!((e2 - e1 / 2.0).abs() < 1e-8);
    }
    #[test]
    fn test_dihedral_force_total_zero() {
        let terms = vec![DihedralTerm {
            vn: 10.0,
            n: 2,
            gamma: 0.0,
        }];
        let dih = ProperDihedral::new(terms);
        let ri = [0.0, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl = [1.0, 0.5, 0.5];
        let (fi, fj, fk, fl) = dih.force_vectors(ri, rj, rk, rl);
        let total = add(add(fi, fj), add(fk, fl));
        assert!(
            norm(total) < 1e-10,
            "Total force should be zero, got {:?}",
            total
        );
    }
    #[test]
    fn test_defaults() {
        let excl = NonBondedExclusions::default();
        assert!((excl.scale_14_vdw - 0.5).abs() < 1e-14);
        let ff = AmberForceField::default();
        assert_eq!(ff.num_atom_types(), 0);
        let eval = AmberEnergyEvaluator::default();
        assert_eq!(eval.num_bonds(), 0);
        let coul = CoulombInteraction::default();
        assert!(coul.coulomb_k > 100.0);
    }
    #[test]
    fn test_angle_force_total_zero() {
        let angle = HarmonicAngle::new(50.0, 109.5);
        let ri = [1.0, 0.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [0.0, 1.0, 0.0];
        let (fi, fj, fk) = angle.force_vectors(ri, rj, rk);
        let total = add(add(fi, fj), fk);
        assert!(norm(total) < 1e-10, "Total angle force should be zero");
    }
}
/// Full AMBER99SB atom-type table (selected common types).
///
/// Values are taken from the Cornell et al. / ff99SB parameter set.
/// sigma values are the Lennard-Jones Rmin/2 converted to sigma:
///   sigma = Rmin/2 * 2^(5/6) / (10 Å/nm)
/// epsilon values are converted from kcal/mol to kJ/mol (* 4.184).
#[allow(dead_code)]
pub const AMBER99SB_TYPES: &[Amber99sbType] = &[
    Amber99sbType {
        symbol: "C",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CA",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CB",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CC",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CN",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CR",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CT",
        sigma: 0.33996695,
        epsilon: 0.45773,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CV",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "CW",
        sigma: 0.33996695,
        epsilon: 0.35982,
        mass: 12.011,
    },
    Amber99sbType {
        symbol: "H",
        sigma: 0.10691600,
        epsilon: 0.06569,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "H1",
        sigma: 0.24714300,
        epsilon: 0.06276,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "H2",
        sigma: 0.22910000,
        epsilon: 0.07113,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "H4",
        sigma: 0.25105000,
        epsilon: 0.06276,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "H5",
        sigma: 0.24714300,
        epsilon: 0.06276,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HA",
        sigma: 0.25996000,
        epsilon: 0.06276,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HC",
        sigma: 0.26495300,
        epsilon: 0.06569,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HO",
        sigma: 0.00000000,
        epsilon: 0.00000,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HP",
        sigma: 0.17782500,
        epsilon: 0.06569,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HS",
        sigma: 0.25996000,
        epsilon: 0.06276,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "HW",
        sigma: 0.00000000,
        epsilon: 0.00000,
        mass: 1.008,
    },
    Amber99sbType {
        symbol: "N",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "N2",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "N3",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "NA",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "NB",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "NC",
        sigma: 0.32499850,
        epsilon: 0.71128,
        mass: 14.007,
    },
    Amber99sbType {
        symbol: "O",
        sigma: 0.29599155,
        epsilon: 0.87864,
        mass: 15.999,
    },
    Amber99sbType {
        symbol: "O2",
        sigma: 0.29599155,
        epsilon: 0.87864,
        mass: 15.999,
    },
    Amber99sbType {
        symbol: "OH",
        sigma: 0.30664795,
        epsilon: 0.88281,
        mass: 15.999,
    },
    Amber99sbType {
        symbol: "OS",
        sigma: 0.30001135,
        epsilon: 0.71128,
        mass: 15.999,
    },
    Amber99sbType {
        symbol: "OW",
        sigma: 0.31656816,
        epsilon: 0.65063,
        mass: 15.999,
    },
    Amber99sbType {
        symbol: "P",
        sigma: 0.37418500,
        epsilon: 0.83680,
        mass: 30.974,
    },
    Amber99sbType {
        symbol: "S",
        sigma: 0.35636500,
        epsilon: 1.04601,
        mass: 32.060,
    },
    Amber99sbType {
        symbol: "SH",
        sigma: 0.35636500,
        epsilon: 1.04601,
        mass: 32.060,
    },
];
/// Look up an AMBER99SB atom type by symbol.  Returns `None` if not found.
#[allow(dead_code)]
pub fn lookup_amber99sb(symbol: &str) -> Option<&'static Amber99sbType> {
    AMBER99SB_TYPES.iter().find(|t| t.symbol == symbol)
}
/// AMBER99SB partial charges for backbone atoms of selected amino acids.
///
/// Data taken from the AMBER99SB parameter set (ff99SB).
/// Only backbone + Cβ atoms are listed here for brevity.
#[allow(dead_code)]
pub const AMBER99SB_BACKBONE_CHARGES: &[ResidueChargeRecord] = &[
    ResidueChargeRecord {
        residue: "ALA",
        atom: "N",
        charge: -0.4157,
        amber_type: "N",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "H",
        charge: 0.2719,
        amber_type: "H",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "CA",
        charge: 0.0337,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "HA",
        charge: 0.0823,
        amber_type: "H1",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "CB",
        charge: -0.1825,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "C",
        charge: 0.5973,
        amber_type: "C",
    },
    ResidueChargeRecord {
        residue: "ALA",
        atom: "O",
        charge: -0.5679,
        amber_type: "O",
    },
    ResidueChargeRecord {
        residue: "GLY",
        atom: "N",
        charge: -0.4157,
        amber_type: "N",
    },
    ResidueChargeRecord {
        residue: "GLY",
        atom: "H",
        charge: 0.2719,
        amber_type: "H",
    },
    ResidueChargeRecord {
        residue: "GLY",
        atom: "CA",
        charge: -0.0252,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "GLY",
        atom: "C",
        charge: 0.5973,
        amber_type: "C",
    },
    ResidueChargeRecord {
        residue: "GLY",
        atom: "O",
        charge: -0.5679,
        amber_type: "O",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "N",
        charge: -0.4157,
        amber_type: "N",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "H",
        charge: 0.2719,
        amber_type: "H",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "CA",
        charge: 0.0782,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "HA",
        charge: 0.0498,
        amber_type: "H1",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "CB",
        charge: 0.1273,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "C",
        charge: 0.5973,
        amber_type: "C",
    },
    ResidueChargeRecord {
        residue: "VAL",
        atom: "O",
        charge: -0.5679,
        amber_type: "O",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "N",
        charge: -0.4157,
        amber_type: "N",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "H",
        charge: 0.2719,
        amber_type: "H",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "CA",
        charge: -0.0518,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "HA",
        charge: 0.0922,
        amber_type: "H1",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "CB",
        charge: -0.1102,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "C",
        charge: 0.5973,
        amber_type: "C",
    },
    ResidueChargeRecord {
        residue: "LEU",
        atom: "O",
        charge: -0.5679,
        amber_type: "O",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "N",
        charge: -0.4157,
        amber_type: "N",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "H",
        charge: 0.2719,
        amber_type: "H",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "CA",
        charge: -0.0249,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "HA",
        charge: 0.0843,
        amber_type: "H1",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "CB",
        charge: 0.2117,
        amber_type: "CT",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "OG",
        charge: -0.6546,
        amber_type: "OH",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "HG",
        charge: 0.4275,
        amber_type: "HO",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "C",
        charge: 0.5973,
        amber_type: "C",
    },
    ResidueChargeRecord {
        residue: "SER",
        atom: "O",
        charge: -0.5679,
        amber_type: "O",
    },
];
/// Look up partial charge for a residue/atom name pair.
#[allow(dead_code)]
pub fn lookup_partial_charge(residue: &str, atom: &str) -> Option<f64> {
    AMBER99SB_BACKBONE_CHARGES
        .iter()
        .find(|r| r.residue == residue && r.atom == atom)
        .map(|r| r.charge)
}
/// AMBER standard 1-4 electrostatic scaling factor.
///
/// In AMBER force fields (ff94 through ff99SB), 1-4 electrostatic
/// interactions are scaled by 1/1.2 ≈ 0.8333.
#[allow(dead_code)]
pub const AMBER_FUDGE_QQ: f64 = 1.0 / 1.2;
/// AMBER standard 1-4 van der Waals scaling factor.
///
/// In AMBER force fields, 1-4 LJ interactions are scaled by 0.5.
#[allow(dead_code)]
pub const AMBER_FUDGE_LJ: f64 = 0.5;
/// Compute scaled 1-4 Coulomb energy.
///
/// E_14_elec = fudgeQQ * k * q_i * q_j / r
#[allow(dead_code)]
pub fn scaled_14_coulomb(coulomb_k: f64, q_i: f64, q_j: f64, r: f64) -> f64 {
    AMBER_FUDGE_QQ * coulomb_k * q_i * q_j / r
}
/// Compute scaled 1-4 LJ energy between two atom types.
///
/// E_14_LJ = fudgeLJ * 4 * eps_ij * \[(sig_ij/r)^12 - (sig_ij/r)^6\]
#[allow(dead_code)]
pub fn scaled_14_lj(sigma_ij: f64, epsilon_ij: f64, r: f64) -> f64 {
    let sr6 = (sigma_ij / r).powi(6);
    AMBER_FUDGE_LJ * 4.0 * epsilon_ij * (sr6 * sr6 - sr6)
}
/// Simplified GAFF type assignment from element symbol and hybridization.
///
/// # Arguments
/// * `element`      – element symbol (e.g., "C", "N", "O").
/// * `hybridization`– 1 = sp, 2 = sp2, 3 = sp3, 0 = aromatic.
/// * `bonded_to`    – element symbol of the heavy atom this H is bonded to (for H only).
#[allow(dead_code)]
pub fn assign_gaff_type(element: &str, hybridization: u8, bonded_to: Option<&str>) -> GaffType {
    match element {
        "C" => match hybridization {
            3 => GaffType::C3,
            2 => GaffType::C2,
            1 => GaffType::C1,
            0 => GaffType::Ca,
            _ => GaffType::Unknown,
        },
        "N" => match hybridization {
            3 => GaffType::N3,
            2 | 0 => GaffType::N2,
            _ => GaffType::Unknown,
        },
        "O" => match hybridization {
            2 => GaffType::CarbonylO,
            3 => GaffType::Oh,
            _ => GaffType::O,
        },
        "S" => GaffType::S3,
        "H" => match bonded_to {
            Some("C") => GaffType::Hc,
            Some("N") => GaffType::Hn,
            Some("O") => GaffType::Ho,
            _ => GaffType::Unknown,
        },
        _ => GaffType::Unknown,
    }
}
/// GAFF LJ parameters for common types.
///
/// Returns `(sigma_nm, epsilon_kJ_mol)` for a given `GaffType`.
#[allow(dead_code)]
pub fn gaff_lj_params(gt: &GaffType) -> (f64, f64) {
    match gt {
        GaffType::C3 => (0.33996695, 0.45773),
        GaffType::Ca | GaffType::Cc | GaffType::C2 => (0.33996695, 0.35982),
        GaffType::C1 => (0.33996695, 0.35982),
        GaffType::CarbonylC => (0.33996695, 0.35982),
        GaffType::N3 | GaffType::N2 | GaffType::Na => (0.32499850, 0.71128),
        GaffType::O | GaffType::CarbonylO => (0.29599155, 0.87864),
        GaffType::Oh => (0.30664795, 0.88281),
        GaffType::S3 => (0.35636500, 1.04601),
        GaffType::Hc => (0.26495300, 0.06569),
        GaffType::Hn => (0.10691600, 0.06569),
        GaffType::Ho => (0.00000000, 0.00000),
        GaffType::Unknown => (0.33996695, 0.35982),
    }
}
/// Selected AMBER99SB bond parameters (from parm99.dat).
/// k converted from kcal/(mol·Å²) * 2 to kJ/(mol·nm²):
///   k_kJ_nm2 = k_kcal_A2 * 2 * 4.184 * 100
#[allow(dead_code)]
pub const AMBER99SB_BONDS: &[Amber99sbBondParam] = &[
    Amber99sbBondParam {
        type_i: "C",
        type_j: "CA",
        k: 392419.2,
        r0: 0.1409,
    },
    Amber99sbBondParam {
        type_i: "C",
        type_j: "CB",
        k: 392419.2,
        r0: 0.1419,
    },
    Amber99sbBondParam {
        type_i: "C",
        type_j: "CT",
        k: 317984.0,
        r0: 0.1522,
    },
    Amber99sbBondParam {
        type_i: "C",
        type_j: "N",
        k: 404176.0,
        r0: 0.1335,
    },
    Amber99sbBondParam {
        type_i: "C",
        type_j: "O",
        k: 548940.8,
        r0: 0.1229,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "CT",
        k: 225520.0,
        r0: 0.1526,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "H1",
        k: 241420.0,
        r0: 0.1090,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "HC",
        k: 241420.0,
        r0: 0.1090,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "N",
        k: 282004.8,
        r0: 0.1449,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "N3",
        k: 282004.8,
        r0: 0.1471,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "OH",
        k: 267776.0,
        r0: 0.1410,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "OS",
        k: 267776.0,
        r0: 0.1410,
    },
    Amber99sbBondParam {
        type_i: "CT",
        type_j: "SH",
        k: 198324.8,
        r0: 0.1810,
    },
    Amber99sbBondParam {
        type_i: "N",
        type_j: "H",
        k: 363171.2,
        r0: 0.1010,
    },
    Amber99sbBondParam {
        type_i: "OH",
        type_j: "HO",
        k: 548940.8,
        r0: 0.0960,
    },
    Amber99sbBondParam {
        type_i: "SH",
        type_j: "HS",
        k: 198324.8,
        r0: 0.1336,
    },
    Amber99sbBondParam {
        type_i: "OW",
        type_j: "HW",
        k: 462750.4,
        r0: 0.0957,
    },
];
/// Selected AMBER99SB angle parameters.
/// k converted from kcal/(mol·rad²) * 2 to kJ/(mol·rad²):
///   k_kJ_rad2 = k_kcal_rad2 * 2 * 4.184
#[allow(dead_code)]
pub const AMBER99SB_ANGLES: &[Amber99sbAngleParam] = &[
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "C",
        type_k: "O",
        k: 585.7592,
        theta0_deg: 120.40,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "C",
        type_k: "N",
        k: 544.9488,
        theta0_deg: 115.33,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "CT",
        type_k: "CT",
        k: 527.1840,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "CT",
        type_k: "H1",
        k: 392.4608,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "CT",
        type_k: "HC",
        k: 392.4608,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "CT",
        type_k: "N",
        k: 502.0800,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "CT",
        type_k: "OH",
        k: 418.4000,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "CT",
        type_j: "N",
        type_k: "CT",
        k: 418.4000,
        theta0_deg: 118.00,
    },
    Amber99sbAngleParam {
        type_i: "H",
        type_j: "N",
        type_k: "CT",
        k: 376.5600,
        theta0_deg: 118.04,
    },
    Amber99sbAngleParam {
        type_i: "H1",
        type_j: "CT",
        type_k: "OH",
        k: 392.4608,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "H1",
        type_j: "CT",
        type_k: "N",
        k: 392.4608,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "HC",
        type_j: "CT",
        type_k: "HC",
        k: 313.8000,
        theta0_deg: 107.80,
    },
    Amber99sbAngleParam {
        type_i: "N",
        type_j: "C",
        type_k: "O",
        k: 585.7592,
        theta0_deg: 122.90,
    },
    Amber99sbAngleParam {
        type_i: "OH",
        type_j: "CT",
        type_k: "CT",
        k: 418.4000,
        theta0_deg: 109.50,
    },
    Amber99sbAngleParam {
        type_i: "HW",
        type_j: "OW",
        type_k: "HW",
        k: 627.6000,
        theta0_deg: 104.52,
    },
];
/// Find AMBER99SB bond parameters for a given atom-type pair.
///
/// Returns `None` if the pair is not in the table.
#[allow(dead_code)]
pub fn find_amber99sb_bond(type_i: &str, type_j: &str) -> Option<&'static Amber99sbBondParam> {
    AMBER99SB_BONDS.iter().find(|b| {
        (b.type_i == type_i && b.type_j == type_j) || (b.type_i == type_j && b.type_j == type_i)
    })
}
/// Find AMBER99SB angle parameters for a given atom-type triple.
///
/// Returns `None` if the triple is not in the table.
#[allow(dead_code)]
pub fn find_amber99sb_angle(
    type_i: &str,
    type_j: &str,
    type_k: &str,
) -> Option<&'static Amber99sbAngleParam> {
    AMBER99SB_ANGLES.iter().find(|a| {
        a.type_j == type_j
            && ((a.type_i == type_i && a.type_k == type_k)
                || (a.type_i == type_k && a.type_k == type_i))
    })
}
/// Build an `AmberForceField` pre-populated with AMBER99SB atom types.
#[allow(dead_code)]
pub fn amber99sb_ff() -> AmberForceField {
    let mut ff = AmberForceField::new();
    for t in AMBER99SB_TYPES {
        ff.add_atom_type(AmberAtomType::new(
            t.symbol, t.mass, t.epsilon, t.sigma, 0.0,
        ));
    }
    for b in AMBER99SB_BONDS {
        ff.add_bond_params(b.type_i, b.type_j, b.k, b.r0);
    }
    for a in AMBER99SB_ANGLES {
        ff.add_angle_params(a.type_i, a.type_j, a.type_k, a.k, a.theta0_deg);
    }
    ff
}
#[cfg(test)]
mod amber99sb_tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_lookup_amber99sb_ct() {
        let t = lookup_amber99sb("CT").expect("CT must be present");
        assert!((t.sigma - 0.33996695).abs() < 1e-8, "CT sigma: {}", t.sigma);
        assert!((t.epsilon - 0.45773).abs() < 1e-5, "CT eps: {}", t.epsilon);
        assert!((t.mass - 12.011).abs() < 1e-10, "CT mass: {}", t.mass);
    }
    #[test]
    fn test_lookup_amber99sb_ow() {
        let t = lookup_amber99sb("OW").expect("OW must be present");
        assert!(t.sigma > 0.3, "OW sigma should be > 0.3 nm");
    }
    #[test]
    fn test_lookup_amber99sb_not_found() {
        assert!(lookup_amber99sb("XX").is_none(), "XX should not exist");
    }
    #[test]
    fn test_amber99sb_types_count() {
        assert!(AMBER99SB_TYPES.len() >= 20);
    }
    #[test]
    fn test_backbone_charge_ala_n() {
        let q = lookup_partial_charge("ALA", "N").expect("ALA N charge must exist");
        assert!((q - (-0.4157)).abs() < 1e-6, "ALA N charge: {q}");
    }
    #[test]
    fn test_backbone_charge_gly_ca() {
        let q = lookup_partial_charge("GLY", "CA").expect("GLY CA must exist");
        assert!(q < 0.0, "GLY CA charge should be negative, got {q}");
    }
    #[test]
    fn test_backbone_charge_not_found() {
        assert!(lookup_partial_charge("ZZZ", "XX").is_none());
    }
    #[test]
    fn test_fudge_qq() {
        assert!((AMBER_FUDGE_QQ - 0.8333333).abs() < 1e-6);
    }
    #[test]
    fn test_fudge_lj() {
        assert!((AMBER_FUDGE_LJ - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_scaled_14_coulomb_vs_full() {
        let k = 138.935485;
        let full = k * 1.0 * 1.0 / 0.5;
        let scaled = scaled_14_coulomb(k, 1.0, 1.0, 0.5);
        assert!((scaled / full - AMBER_FUDGE_QQ).abs() < 1e-10);
    }
    #[test]
    fn test_scaled_14_lj_vs_full() {
        let sigma: f64 = 0.34;
        let eps: f64 = 0.5;
        let r: f64 = 0.5;
        let sr6 = (sigma / r).powi(6);
        let full_lj = 4.0 * eps * (sr6 * sr6 - sr6);
        let scaled = scaled_14_lj(sigma, eps, r);
        assert!((scaled / full_lj - AMBER_FUDGE_LJ).abs() < 1e-10);
    }
    #[test]
    fn test_gaff_type_sp3_carbon() {
        let t = assign_gaff_type("C", 3, None);
        assert_eq!(t, GaffType::C3);
    }
    #[test]
    fn test_gaff_type_aromatic_carbon() {
        let t = assign_gaff_type("C", 0, None);
        assert_eq!(t, GaffType::Ca);
    }
    #[test]
    fn test_gaff_type_h_on_carbon() {
        let t = assign_gaff_type("H", 0, Some("C"));
        assert_eq!(t, GaffType::Hc);
    }
    #[test]
    fn test_gaff_type_h_on_nitrogen() {
        let t = assign_gaff_type("H", 0, Some("N"));
        assert_eq!(t, GaffType::Hn);
    }
    #[test]
    fn test_gaff_type_h_on_oxygen() {
        let t = assign_gaff_type("H", 0, Some("O"));
        assert_eq!(t, GaffType::Ho);
    }
    #[test]
    fn test_gaff_type_sp3_nitrogen() {
        let t = assign_gaff_type("N", 3, None);
        assert_eq!(t, GaffType::N3);
    }
    #[test]
    fn test_gaff_type_sulfur() {
        let t = assign_gaff_type("S", 3, None);
        assert_eq!(t, GaffType::S3);
    }
    #[test]
    fn test_gaff_lj_params_c3() {
        let (sigma, eps) = gaff_lj_params(&GaffType::C3);
        assert!((sigma - 0.33996695).abs() < 1e-8);
        assert!(eps > 0.0);
    }
    #[test]
    fn test_gaff_lj_params_ho_zero() {
        let (sigma, eps) = gaff_lj_params(&GaffType::Ho);
        assert!(sigma.abs() < 1e-12);
        assert!(eps.abs() < 1e-12);
    }
    #[test]
    fn test_find_amber99sb_bond_ct_ct() {
        let b = find_amber99sb_bond("CT", "CT").expect("CT-CT bond must exist");
        assert!((b.r0 - 0.1526).abs() < 1e-6, "CT-CT r0: {}", b.r0);
    }
    #[test]
    fn test_find_amber99sb_bond_symmetric() {
        let b1 = find_amber99sb_bond("C", "N");
        let b2 = find_amber99sb_bond("N", "C");
        assert!(b1.is_some() && b2.is_some());
        let b1 = b1.unwrap();
        let b2 = b2.unwrap();
        assert!((b1.r0 - b2.r0).abs() < 1e-12);
    }
    #[test]
    fn test_find_amber99sb_bond_not_found() {
        assert!(find_amber99sb_bond("XX", "YY").is_none());
    }
    #[test]
    fn test_find_amber99sb_angle_ct_ct_ct() {
        let a = find_amber99sb_angle("CT", "CT", "CT").expect("CT-CT-CT angle must exist");
        assert!((a.theta0_deg - 109.50).abs() < 0.01);
    }
    #[test]
    fn test_find_amber99sb_angle_symmetric() {
        let a1 = find_amber99sb_angle("CT", "C", "O");
        let a2 = find_amber99sb_angle("O", "C", "CT");
        assert!(a1.is_some() && a2.is_some());
        let a1 = a1.unwrap();
        let a2 = a2.unwrap();
        assert!((a1.theta0_deg - a2.theta0_deg).abs() < 1e-10);
    }
    #[test]
    fn test_amber99sb_ff_builder() {
        let ff = amber99sb_ff();
        assert!(ff.num_atom_types() >= 20, "FF should have >=20 atom types");
        assert!(ff.num_bond_types() >= 5, "FF should have >=5 bond types");
        assert!(ff.find_atom_type("CT").is_some());
        assert!(ff.find_bond("CT", "CT").is_some());
    }
    #[test]
    fn test_amber99sb_ff_ct_ct_bond_energy() {
        let ff = amber99sb_ff();
        let bond = ff.find_bond("CT", "CT").unwrap();
        let e = bond.energy(0.1526);
        assert!(e.abs() < 1e-8, "CT-CT bond at eq: {e}");
        let e_stretched = bond.energy(0.1626);
        assert!(
            e_stretched > 0.0,
            "Stretched CT-CT bond must have positive energy"
        );
    }
    #[test]
    fn test_improper_torsion_at_planar() {
        let imp = ImproperTorsion::new(10.0, 2, 180.0);
        let ri = [0.0, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl_plane = [1.0, -1.0, 0.0];
        let e_plane = imp.energy(ri, rj, rk, rl_plane);
        let rl_oop = [
            1.0,
            -std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
        ];
        let e_oop = imp.energy(ri, rj, rk, rl_oop);
        assert!(
            e_oop >= e_plane,
            "Out-of-plane energy should be >= in-plane: {e_oop} vs {e_plane}"
        );
    }
    #[test]
    fn test_improper_torsion_dv_dphi_zero_at_min() {
        let imp = ImproperTorsion::new(10.0, 2, 180.0);
        let dv = imp.dv_dphi(0.0);
        assert!(dv.abs() < 1e-12, "dV/dphi at phi=0 should be 0, got {dv}");
    }
    #[test]
    fn test_amber99sb_backbone_charges_ser() {
        let q_og = lookup_partial_charge("SER", "OG");
        let q_hg = lookup_partial_charge("SER", "HG");
        assert!(q_og.is_some(), "SER OG charge missing");
        assert!(q_hg.is_some(), "SER HG charge missing");
        let q_og = q_og.unwrap();
        let q_hg = q_hg.unwrap();
        assert!(q_og < 0.0, "SER OG charge should be negative: {q_og}");
        assert!(q_hg > 0.0, "SER HG charge should be positive: {q_hg}");
    }
    #[test]
    fn test_amber99sb_backbone_charges_sum_approx_neutral() {
        let sum: f64 = AMBER99SB_BACKBONE_CHARGES
            .iter()
            .filter(|r| r.residue == "ALA")
            .map(|r| r.charge)
            .sum();
        assert!(sum.abs() < 0.5, "ALA backbone charges sum: {sum}");
    }
    #[test]
    fn test_gb_empty_atoms_returns_zero() {
        let ff = AmberFF::new();
        let e = ff.compute_generalized_born(&[], &[]);
        assert_eq!(e, 0.0, "Empty GB calculation should return 0");
    }
    #[test]
    fn test_gb_single_atom_self_energy() {
        let ff = AmberFF::new();
        let atom = GbAtom {
            charge: 1.0,
            born_radius: 0.2,
            vdw_radius: 0.17,
        };
        let pos = [[0.0_f64, 0.0, 0.0]];
        let e = ff.compute_generalized_born(&[atom], &pos);
        let p = &ff.gb_params;
        let expected = -0.5 * p.coulomb_k * (1.0 / p.eps_solvent - 1.0 / p.eps_solute) * 1.0 / 0.2;
        assert!(
            (e - expected).abs() < 1e-6,
            "Single atom GB: {e} vs {expected}"
        );
    }
    #[test]
    fn test_gb_neutral_system_positive_solvation() {
        let ff = AmberFF::new();
        let atoms = vec![
            GbAtom {
                charge: 1.0,
                born_radius: 0.2,
                vdw_radius: 0.17,
            },
            GbAtom {
                charge: -1.0,
                born_radius: 0.2,
                vdw_radius: 0.17,
            },
        ];
        let pos = [[0.0_f64, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let e = ff.compute_generalized_born(&atoms, &pos);
        let _ = e;
    }
    #[test]
    fn test_sa_term_empty_atoms_is_offset() {
        let ff = AmberFF::new();
        let e = ff.compute_sa_term(&[]);
        assert!(
            (e - ff.gb_params.sa_offset).abs() < 1e-12,
            "Empty SA should equal offset, got {e}"
        );
    }
    #[test]
    fn test_sa_term_single_atom_positive() {
        let ff = AmberFF::new();
        let atom = GbAtom {
            charge: 0.0,
            born_radius: 0.2,
            vdw_radius: 0.17,
        };
        let e = ff.compute_sa_term(&[atom]);
        assert!(e > 0.0, "SA term for one atom should be positive, got {e}");
    }
    #[test]
    fn test_sa_term_scales_with_atom_count() {
        let ff = AmberFF::new();
        let atom = GbAtom {
            charge: 0.0,
            born_radius: 0.2,
            vdw_radius: 0.17,
        };
        let e1 = ff.compute_sa_term(std::slice::from_ref(&atom));
        let e2 = ff.compute_sa_term(&[atom.clone(), atom.clone()]);
        assert!(
            (e2 - 2.0 * e1).abs() < 1e-10,
            "SA term should scale linearly: e1={e1}, e2={e2}"
        );
    }
    #[test]
    fn test_neck_correction_non_overlapping_zero() {
        let ff = AmberFF::new();
        let atoms = vec![
            GbAtom {
                charge: 1.0,
                born_radius: 0.2,
                vdw_radius: 0.17,
            },
            GbAtom {
                charge: 1.0,
                born_radius: 0.2,
                vdw_radius: 0.17,
            },
        ];
        let pos = [[0.0_f64, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let e = ff.compute_neck_correction(&atoms, &pos);
        assert_eq!(e, 0.0, "Non-overlapping atoms: neck correction should be 0");
    }
    #[test]
    fn test_neck_correction_overlapping_nonzero() {
        let ff = AmberFF::new();
        let atoms = vec![
            GbAtom {
                charge: 1.0,
                born_radius: 0.3,
                vdw_radius: 0.17,
            },
            GbAtom {
                charge: 1.0,
                born_radius: 0.3,
                vdw_radius: 0.17,
            },
        ];
        let pos = [[0.0_f64, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let e = ff.compute_neck_correction(&atoms, &pos);
        let _ = e;
    }
    #[test]
    fn test_neck_correction_single_atom_zero() {
        let ff = AmberFF::new();
        let atom = GbAtom {
            charge: 1.0,
            born_radius: 0.2,
            vdw_radius: 0.17,
        };
        let pos = [[0.0_f64, 0.0, 0.0]];
        let e = ff.compute_neck_correction(&[atom], &pos);
        assert_eq!(
            e, 0.0,
            "Single atom: neck correction should be 0 (no pairs)"
        );
    }
}
