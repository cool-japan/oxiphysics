//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;
use crate::potential::Potential;

use super::types::{
    AngleTerm, Bond, CgenffForceField, DihedralTerm, GromosForceField, ImproperDihedral,
    IntraEnergyDecomposition, LjAtomType, MixingRule, UreyBradleyTerm, ValidationError,
    ValidationSeverity,
};

/// Trait for force fields that compute forces on atoms.
pub trait ForceField: Send + Sync {
    /// Compute forces on all atoms and return the potential energy.
    ///
    /// This should **add** to `atoms.forces` (not clear them first)
    /// so that multiple force fields can be composed.
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64;
}
/// Helper trait to allow cloning boxed potentials.
pub trait PotentialClone: Potential {
    /// Clone into a boxed trait object.
    fn clone_box(&self) -> Box<dyn PotentialClone>;
}
impl<T: Potential + Clone + 'static> PotentialClone for T {
    fn clone_box(&self) -> Box<dyn PotentialClone> {
        Box::new(self.clone())
    }
}
impl Clone for Box<dyn PotentialClone> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
/// Apply mixing rules to get cross-type LJ parameters.
///
/// Returns `(sigma_ij, epsilon_ij)`.
pub fn apply_mixing_rule(rule: MixingRule, ti: &LjAtomType, tj: &LjAtomType) -> (f64, f64) {
    match rule {
        MixingRule::LorentzBerthelot => {
            let sigma = (ti.sigma + tj.sigma) / 2.0;
            let epsilon = (ti.epsilon * tj.epsilon).sqrt();
            (sigma, epsilon)
        }
        MixingRule::Geometric => {
            let sigma = (ti.sigma * tj.sigma).sqrt();
            let epsilon = (ti.epsilon * tj.epsilon).sqrt();
            (sigma, epsilon)
        }
        MixingRule::Arithmetic => {
            let sigma = (ti.sigma + tj.sigma) / 2.0;
            let epsilon = (ti.epsilon + tj.epsilon) / 2.0;
            (sigma, epsilon)
        }
        MixingRule::WaldmanHagler => {
            let si6 = ti.sigma.powi(6);
            let sj6 = tj.sigma.powi(6);
            let sigma = ((si6 + sj6) / 2.0).powf(1.0 / 6.0);
            let epsilon = if si6 + sj6 > 1e-30 {
                2.0 * (ti.epsilon * tj.epsilon).sqrt() * ti.sigma.powi(3) * tj.sigma.powi(3)
                    / (si6 + sj6)
            } else {
                0.0
            };
            (sigma, epsilon)
        }
    }
}
/// Generate all cross-type LJ parameters from a set of atom types.
pub fn generate_cross_parameters(
    types: &[LjAtomType],
    rule: MixingRule,
) -> Vec<(u32, u32, f64, f64)> {
    let mut params = Vec::new();
    for i in 0..types.len() {
        for j in i..types.len() {
            let (sigma, epsilon) = apply_mixing_rule(rule, &types[i], &types[j]);
            params.push((types[i].type_id, types[j].type_id, sigma, epsilon));
        }
    }
    params
}
/// Validate basic force field parameters.
///
/// Checks for:
/// - Negative spring constants
/// - Zero or negative masses
/// - Bond indices out of range
/// - Unreasonable equilibrium distances
///
/// Returns a list of validation errors/warnings.
pub fn validate_force_field(
    n_atoms: usize,
    bonds: &[Bond],
    angles: &[AngleTerm],
    dihedrals: &[DihedralTerm],
) -> Vec<ValidationError> {
    let mut errors = Vec::new();
    for (idx, bond) in bonds.iter().enumerate() {
        if bond.i >= n_atoms || bond.j >= n_atoms {
            errors.push(ValidationError {
                message: format!(
                    "Bond {idx}: atom index out of range (i={}, j={}, n_atoms={n_atoms})",
                    bond.i, bond.j
                ),
                severity: ValidationSeverity::Error,
            });
        }
        if bond.i == bond.j {
            errors.push(ValidationError {
                message: format!("Bond {idx}: self-bond (i == j == {})", bond.i),
                severity: ValidationSeverity::Error,
            });
        }
        if bond.k < 0.0 {
            errors.push(ValidationError {
                message: format!("Bond {idx}: negative spring constant k={}", bond.k),
                severity: ValidationSeverity::Error,
            });
        }
        if bond.r0 <= 0.0 {
            errors.push(ValidationError {
                message: format!(
                    "Bond {idx}: non-positive equilibrium distance r0={}",
                    bond.r0
                ),
                severity: ValidationSeverity::Warning,
            });
        }
    }
    for (idx, angle) in angles.iter().enumerate() {
        if angle.i >= n_atoms || angle.j >= n_atoms || angle.k >= n_atoms {
            errors.push(ValidationError {
                message: format!(
                    "Angle {idx}: atom index out of range (i={}, j={}, k={}, n_atoms={n_atoms})",
                    angle.i, angle.j, angle.k
                ),
                severity: ValidationSeverity::Error,
            });
        }
        if angle.k_theta < 0.0 {
            errors.push(ValidationError {
                message: format!(
                    "Angle {idx}: negative force constant k_theta={}",
                    angle.k_theta
                ),
                severity: ValidationSeverity::Error,
            });
        }
    }
    for (idx, dih) in dihedrals.iter().enumerate() {
        if dih.i >= n_atoms || dih.j >= n_atoms || dih.k >= n_atoms || dih.l >= n_atoms {
            errors
                .push(ValidationError {
                    message: format!(
                        "Dihedral {idx}: atom index out of range (i={}, j={}, k={}, l={}, n_atoms={n_atoms})",
                        dih.i, dih.j, dih.k, dih.l
                    ),
                    severity: ValidationSeverity::Error,
                });
        }
    }
    errors
}
/// Check if force field has any critical errors.
pub fn has_critical_errors(errors: &[ValidationError]) -> bool {
    errors
        .iter()
        .any(|e| e.severity == ValidationSeverity::Error)
}
/// Count warnings in validation results.
pub fn count_warnings(errors: &[ValidationError]) -> usize {
    errors
        .iter()
        .filter(|e| e.severity == ValidationSeverity::Warning)
        .count()
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::potential::LennardJones;
    use oxiphysics_core::math::Vec3;

    #[test]
    fn test_pair_force_field_two_atoms() {
        let mut atoms = AtomSet::new();
        let sigma = 1.0;
        let r = 2.0_f64.powf(1.0 / 6.0) * sigma;
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(r, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let mut ff = PairForceField::new();
        ff.add_interaction(0, 0, LennardJones::new(1.0, sigma, 5.0));
        let pbox = PeriodicBox::cubic(10.0);
        let energy = ff.compute_forces(&mut atoms, &pbox);
        assert!(atoms.forces[0].norm() < 1e-8);
        assert!(atoms.forces[1].norm() < 1e-8);
        assert!(energy < 0.0);
    }
    #[test]
    fn test_bonded_force_field() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(2.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let mut ff = BondedForceField::new();
        ff.add_bond(0, 1, 100.0, 1.5);
        let pbox = PeriodicBox::cubic(10.0);
        let energy = ff.compute_forces(&mut atoms, &pbox);
        assert!((energy - 12.5).abs() < 1e-10);
        assert!(atoms.forces[0].x > 0.0);
        assert!(atoms.forces[1].x < 0.0);
    }
    #[test]
    fn test_lorentz_berthelot_mixing() {
        let ti = LjAtomType {
            type_id: 0,
            sigma: 3.0,
            epsilon: 1.0,
        };
        let tj = LjAtomType {
            type_id: 1,
            sigma: 4.0,
            epsilon: 4.0,
        };
        let (s, e) = apply_mixing_rule(MixingRule::LorentzBerthelot, &ti, &tj);
        assert!((s - 3.5).abs() < 1e-12, "sigma LB: {s}");
        assert!((e - 2.0).abs() < 1e-12, "epsilon LB: {e}");
    }
    #[test]
    fn test_geometric_mixing() {
        let ti = LjAtomType {
            type_id: 0,
            sigma: 2.0,
            epsilon: 2.0,
        };
        let tj = LjAtomType {
            type_id: 1,
            sigma: 8.0,
            epsilon: 8.0,
        };
        let (s, e) = apply_mixing_rule(MixingRule::Geometric, &ti, &tj);
        assert!((s - 4.0).abs() < 1e-12, "sigma geo: {s}");
        assert!((e - 4.0).abs() < 1e-12, "epsilon geo: {e}");
    }
    #[test]
    fn test_arithmetic_mixing() {
        let ti = LjAtomType {
            type_id: 0,
            sigma: 2.0,
            epsilon: 1.0,
        };
        let tj = LjAtomType {
            type_id: 1,
            sigma: 4.0,
            epsilon: 3.0,
        };
        let (s, e) = apply_mixing_rule(MixingRule::Arithmetic, &ti, &tj);
        assert!((s - 3.0).abs() < 1e-12, "sigma arith: {s}");
        assert!((e - 2.0).abs() < 1e-12, "epsilon arith: {e}");
    }
    #[test]
    fn test_waldman_hagler_mixing() {
        let ti = LjAtomType {
            type_id: 0,
            sigma: 3.0,
            epsilon: 1.0,
        };
        let tj = LjAtomType {
            type_id: 1,
            sigma: 3.0,
            epsilon: 1.0,
        };
        let (s, e) = apply_mixing_rule(MixingRule::WaldmanHagler, &ti, &tj);
        assert!((s - 3.0).abs() < 1e-12, "sigma WH same type: {s}");
        assert!((e - 1.0).abs() < 1e-12, "epsilon WH same type: {e}");
    }
    #[test]
    fn test_self_mixing_identity() {
        let t = LjAtomType {
            type_id: 0,
            sigma: 3.5,
            epsilon: 0.5,
        };
        for rule in &[
            MixingRule::LorentzBerthelot,
            MixingRule::Geometric,
            MixingRule::Arithmetic,
        ] {
            let (s, e) = apply_mixing_rule(*rule, &t, &t);
            assert!((s - 3.5).abs() < 1e-12, "self-mixing sigma {rule:?}: {s}");
            assert!((e - 0.5).abs() < 1e-12, "self-mixing epsilon {rule:?}: {e}");
        }
    }
    #[test]
    fn test_generate_cross_parameters() {
        let types = vec![
            LjAtomType {
                type_id: 0,
                sigma: 3.0,
                epsilon: 1.0,
            },
            LjAtomType {
                type_id: 1,
                sigma: 4.0,
                epsilon: 2.0,
            },
        ];
        let params = generate_cross_parameters(&types, MixingRule::LorentzBerthelot);
        assert_eq!(params.len(), 3, "should have 3 cross pairs");
    }
    #[test]
    fn test_14_scaling_amber() {
        let s = OneFourScaling::amber();
        assert!((s.lj_scale - 0.5).abs() < 1e-12);
        assert!((s.coulomb_scale - 1.0 / 1.2).abs() < 1e-12);
    }
    #[test]
    fn test_14_scaling_opls() {
        let s = OneFourScaling::opls();
        assert!((s.lj_scale - 0.5).abs() < 1e-12);
        assert!((s.coulomb_scale - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_14_scale_energy() {
        let s = OneFourScaling::amber();
        let e = s.scale_lj(10.0);
        assert!((e - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_opls_creation() {
        let mut opls = OplsForceField::new();
        opls.add_atom_type("opls_135", "C", -0.18, 0.35, 0.276, 12.011);
        opls.add_atom_type("opls_140", "H", 0.06, 0.25, 0.126, 1.008);
        assert_eq!(opls.n_atom_types(), 2);
        assert!(opls.find_atom_type("opls_135").is_some());
        assert!(opls.find_atom_type("opls_999").is_none());
    }
    #[test]
    fn test_opls_dihedral_energy_zero() {
        let e = OplsForceField::opls_dihedral_energy(0.0, 2.0, 1.0, 3.0, 0.5);
        assert!((e - 5.0).abs() < 1e-12, "OPLS dih at phi=0: {e}");
    }
    #[test]
    fn test_opls_dihedral_energy_pi() {
        let e = OplsForceField::opls_dihedral_energy(std::f64::consts::PI, 2.0, 1.0, 3.0, 0.5);
        assert!(e.abs() < 1e-12, "OPLS dih at phi=pi: {e}");
    }
    #[test]
    fn test_charmm_creation() {
        let mut charmm = CharmmForceField::new();
        charmm.add_atom_type("CT1", "C", -0.09, 2.275, -0.02, 12.011);
        assert_eq!(charmm.n_atom_types(), 1);
    }
    #[test]
    fn test_charmm_rmin_to_sigma() {
        let sigma = CharmmForceField::rmin_half_to_sigma(1.0);
        let expected = 2.0 / 2.0_f64.powf(1.0 / 6.0);
        assert!(
            (sigma - expected).abs() < 1e-12,
            "sigma: {sigma}, expected {expected}"
        );
    }
    #[test]
    fn test_charmm_urey_bradley_energy() {
        let e = CharmmForceField::urey_bradley_energy(2.5, 100.0, 2.0);
        assert!((e - 25.0).abs() < 1e-12, "UB energy: {e}");
    }
    #[test]
    fn test_charmm_improper_energy() {
        let e = CharmmForceField::improper_energy(0.1, 50.0, 0.0);
        assert!((e - 0.5).abs() < 1e-12, "improper energy: {e}");
    }
    #[test]
    fn test_validate_valid_ff() {
        let bonds = vec![Bond {
            i: 0,
            j: 1,
            k: 100.0,
            r0: 1.5,
        }];
        let angles = vec![];
        let dihedrals = vec![];
        let errors = validate_force_field(5, &bonds, &angles, &dihedrals);
        assert!(
            !has_critical_errors(&errors),
            "Valid FF should have no critical errors"
        );
    }
    #[test]
    fn test_validate_out_of_range() {
        let bonds = vec![Bond {
            i: 0,
            j: 10,
            k: 100.0,
            r0: 1.5,
        }];
        let errors = validate_force_field(5, &bonds, &[], &[]);
        assert!(
            has_critical_errors(&errors),
            "Should detect out-of-range index"
        );
    }
    #[test]
    fn test_validate_self_bond() {
        let bonds = vec![Bond {
            i: 2,
            j: 2,
            k: 100.0,
            r0: 1.5,
        }];
        let errors = validate_force_field(5, &bonds, &[], &[]);
        assert!(has_critical_errors(&errors), "Should detect self-bond");
    }
    #[test]
    fn test_validate_negative_k() {
        let bonds = vec![Bond {
            i: 0,
            j: 1,
            k: -100.0,
            r0: 1.5,
        }];
        let errors = validate_force_field(5, &bonds, &[], &[]);
        assert!(has_critical_errors(&errors), "Should detect negative k");
    }
    #[test]
    fn test_validate_nonpositive_r0() {
        let bonds = vec![Bond {
            i: 0,
            j: 1,
            k: 100.0,
            r0: 0.0,
        }];
        let errors = validate_force_field(5, &bonds, &[], &[]);
        let warnings = count_warnings(&errors);
        assert!(warnings > 0, "Should warn about non-positive r0");
    }
    #[test]
    fn test_validate_angle_out_of_range() {
        let angles = vec![AngleTerm {
            i: 0,
            j: 1,
            k: 100,
            k_theta: 50.0,
            theta0: 1.9,
        }];
        let errors = validate_force_field(5, &[], &angles, &[]);
        assert!(
            has_critical_errors(&errors),
            "Should detect angle atom out of range"
        );
    }
    #[test]
    fn test_validate_dihedral_out_of_range() {
        let dihedrals = vec![DihedralTerm {
            i: 0,
            j: 1,
            k: 2,
            l: 99,
            v_n: 1.0,
            gamma: 0.0,
            n: 1,
        }];
        let errors = validate_force_field(5, &[], &[], &dihedrals);
        assert!(
            has_critical_errors(&errors),
            "Should detect dihedral atom out of range"
        );
    }
}
/// Compute intramolecular energy decomposition for a CHARMM force field setup.
pub fn decompose_intramolecular_energy(
    atoms: &crate::atom::AtomSet,
    bonds: &[Bond],
    angles: &[AngleTerm],
    dihedrals: &[DihedralTerm],
    impropers: &[ImproperDihedral],
    ub_terms: &[UreyBradleyTerm],
    pbox: &crate::neighbor::PeriodicBox,
) -> IntraEnergyDecomposition {
    let mut decomp = IntraEnergyDecomposition::default();
    for bond in bonds {
        let (_, dist) =
            crate::neighbor::distance_pbc(&atoms.positions[bond.i], &atoms.positions[bond.j], pbox);
        let delta = dist - bond.r0;
        decomp.bond_energy += 0.5 * bond.k * delta * delta;
    }
    for angle in angles {
        let (rji, dji) = crate::neighbor::distance_pbc(
            &atoms.positions[angle.j],
            &atoms.positions[angle.i],
            pbox,
        );
        let (rjk, djk) = crate::neighbor::distance_pbc(
            &atoms.positions[angle.j],
            &atoms.positions[angle.k],
            pbox,
        );
        if dji > 1e-15 && djk > 1e-15 {
            let cos_t = (rji / dji).dot(&(rjk / djk)).clamp(-1.0, 1.0);
            let theta = cos_t.acos();
            let delta = theta - angle.theta0;
            decomp.angle_energy += 0.5 * angle.k_theta * delta * delta;
        }
    }
    for dih in dihedrals {
        let (r_ij, _) =
            crate::neighbor::distance_pbc(&atoms.positions[dih.i], &atoms.positions[dih.j], pbox);
        let (r_jk, _) =
            crate::neighbor::distance_pbc(&atoms.positions[dih.j], &atoms.positions[dih.k], pbox);
        let (r_kl, _) =
            crate::neighbor::distance_pbc(&atoms.positions[dih.k], &atoms.positions[dih.l], pbox);
        let n1 = r_ij.cross(&r_jk);
        let n2 = r_jk.cross(&r_kl);
        let n1n = n1.norm();
        let n2n = n2.norm();
        if n1n > 1e-15 && n2n > 1e-15 {
            let cos_phi = (n1 / n1n).dot(&(n2 / n2n)).clamp(-1.0, 1.0);
            let m1 = (n1 / n1n).cross(&(n2 / n2n));
            let rjk_hat = r_jk / r_jk.norm().max(1e-15);
            let phi = m1.dot(&rjk_hat).atan2(cos_phi);
            let nf = dih.n as f64;
            decomp.dihedral_energy += dih.v_n * (1.0 + (nf * phi - dih.gamma).cos());
        }
    }
    for imp in impropers {
        let (r_ij, _) =
            crate::neighbor::distance_pbc(&atoms.positions[imp.i], &atoms.positions[imp.j], pbox);
        let (r_ik, _) =
            crate::neighbor::distance_pbc(&atoms.positions[imp.i], &atoms.positions[imp.k], pbox);
        let (r_il, _) =
            crate::neighbor::distance_pbc(&atoms.positions[imp.i], &atoms.positions[imp.l], pbox);
        let n_plane = r_ij.cross(&r_ik);
        let nn = n_plane.norm();
        if nn > 1e-15 {
            let n_hat = n_plane / nn;
            let sin_xi = n_hat.dot(&r_il) / r_il.norm().max(1e-15);
            let xi = sin_xi.clamp(-1.0, 1.0).asin();
            let delta = xi - imp.phi0;
            decomp.improper_energy += imp.k_imp * delta * delta;
        }
    }
    for ub in ub_terms {
        let (_, dist) =
            crate::neighbor::distance_pbc(&atoms.positions[ub.i], &atoms.positions[ub.k], pbox);
        let delta = dist - ub.r0;
        decomp.urey_bradley_energy += ub.k_ub * delta * delta;
    }
    decomp
}
/// Validate GROMOS force field parameters.
pub fn validate_gromos(ff: &GromosForceField, n_atoms: usize) -> Vec<ValidationError> {
    let mut errors = Vec::new();
    for (idx, at) in ff.atom_types.iter().enumerate() {
        if at.c6 < 0.0 {
            errors.push(ValidationError {
                message: format!(
                    "GROMOS atom type {idx} ({}): negative C6={}",
                    at.name, at.c6
                ),
                severity: ValidationSeverity::Error,
            });
        }
        if at.c12 < 0.0 {
            errors.push(ValidationError {
                message: format!(
                    "GROMOS atom type {idx} ({}): negative C12={}",
                    at.name, at.c12
                ),
                severity: ValidationSeverity::Error,
            });
        }
        if at.mass <= 0.0 {
            errors.push(ValidationError {
                message: format!(
                    "GROMOS atom type {idx} ({}): non-positive mass={}",
                    at.name, at.mass
                ),
                severity: ValidationSeverity::Warning,
            });
        }
    }
    let base = validate_force_field(n_atoms, &ff.bonds, &ff.angles, &ff.dihedrals);
    errors.extend(base);
    errors
}
/// Validate CGenFF force field parameters.
pub fn validate_cgenff(ff: &CgenffForceField, n_atoms: usize) -> Vec<ValidationError> {
    let mut errors = Vec::new();
    for at in &ff.atom_types {
        if at.charge_penalty > ff.charge_penalty_threshold {
            errors.push(ValidationError {
                message: format!(
                    "CGenFF atom type {} has high charge penalty {:.1} (threshold {:.1})",
                    at.name, at.charge_penalty, ff.charge_penalty_threshold
                ),
                severity: ValidationSeverity::Warning,
            });
        }
        if at.lj_penalty > ff.lj_penalty_threshold {
            errors.push(ValidationError {
                message: format!(
                    "CGenFF atom type {} has high LJ penalty {:.1} (threshold {:.1})",
                    at.name, at.lj_penalty, ff.lj_penalty_threshold
                ),
                severity: ValidationSeverity::Warning,
            });
        }
    }
    let base = validate_force_field(n_atoms, &ff.bonds, &ff.angles, &[]);
    errors.extend(base);
    errors
}
/// Compute the dihedral (torsion) angle φ for atoms I-J-K-L (radians, −π…+π).
///
/// Standard definition: angle between the I-J-K plane and the J-K-L plane.
pub fn dihedral_angle(ri: [f64; 3], rj: [f64; 3], rk: [f64; 3], rl: [f64; 3]) -> f64 {
    let b1 = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let b2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
    let b3 = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];
    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let n1_len = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
    let n2_len = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
    if n1_len < 1e-15 || n2_len < 1e-15 {
        return 0.0;
    }
    let cos_phi = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (n1_len * n2_len);
    let cos_phi = cos_phi.clamp(-1.0, 1.0);
    let b2_hat = [b2[0] / n1_len, b2[1] / n1_len, b2[2] / n1_len];
    let m1 = cross3(n1, n2);
    let sign = if m1[0] * b2_hat[0] + m1[1] * b2_hat[1] + m1[2] * b2_hat[2] < 0.0 {
        -1.0
    } else {
        1.0
    };
    sign * cos_phi.acos()
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests_extended {
    use super::super::types::*;
    use super::*;
    use std::f64::consts::PI;
    #[test]
    fn test_gromos_creation() {
        let mut ff = GromosForceField::new();
        ff.add_atom_type("CH2", 14.027, 0.0, 0.0088, 1.23e-4);
        assert_eq!(ff.n_atom_types(), 1);
        assert!(ff.find_atom_type("CH2").is_some());
        assert!(ff.find_atom_type("CH3").is_none());
    }
    #[test]
    fn test_gromos_lj_energy_at_minimum() {
        let c6: f64 = 0.0088;
        let c12: f64 = 1.23e-4;
        let r_min = (2.0 * c12 / c6).powf(1.0 / 6.0);
        let e = GromosForceField::gromos_lj_energy(r_min, c6, c12);
        let e_min = -c6 * c6 / (4.0 * c12);
        assert!((e - e_min).abs() < 1e-10, "e={e}, e_min={e_min}");
    }
    #[test]
    fn test_gromos_lj_force_at_minimum_zero() {
        let c6: f64 = 0.0088;
        let c12: f64 = 1.23e-4;
        let r_min = (2.0 * c12 / c6).powf(1.0 / 6.0);
        let f = GromosForceField::gromos_lj_force(r_min, c6, c12);
        assert!(f.abs() < 1e-8, "Force at minimum should be ~0: {f}");
    }
    #[test]
    fn test_gromos_c6c12_roundtrip() {
        let sigma = 0.35;
        let epsilon = 0.276;
        let (c6, c12) = GromosForceField::sigma_epsilon_to_c6c12(sigma, epsilon);
        let (s2, e2) = GromosForceField::c6c12_to_sigma_epsilon(c6, c12);
        assert!((s2 - sigma).abs() < 1e-10, "sigma roundtrip: {s2}");
        assert!((e2 - epsilon).abs() < 1e-10, "epsilon roundtrip: {e2}");
    }
    #[test]
    fn test_gromos_validate_valid() {
        let ff = GromosForceField::new();
        let errors = validate_gromos(&ff, 10);
        assert!(!has_critical_errors(&errors));
    }
    #[test]
    fn test_gromos_validate_negative_c6() {
        let mut ff = GromosForceField::new();
        ff.add_atom_type("BAD", 12.0, 0.0, -1.0, 1e-4);
        let errors = validate_gromos(&ff, 10);
        assert!(has_critical_errors(&errors));
    }
    #[test]
    fn test_gromos_validate_negative_c12() {
        let mut ff = GromosForceField::new();
        ff.add_atom_type("BAD", 12.0, 0.0, 1e-3, -1e-4);
        let errors = validate_gromos(&ff, 10);
        assert!(has_critical_errors(&errors));
    }
    #[test]
    fn test_cgenff_creation() {
        let mut ff = CgenffForceField::new();
        ff.add_atom_type("CG2R61", "C", -0.115, 1.9924, -0.07, 12.011);
        assert_eq!(ff.n_atom_types(), 1);
        assert!(ff.find_atom_type("CG2R61").is_some());
    }
    #[test]
    fn test_cgenff_high_penalty_detection() {
        let mut ff = CgenffForceField::new();
        ff.add_atom_type("CUSTOM", "C", 0.1, 1.9, -0.05, 12.0);
        ff.atom_types[0].charge_penalty = 50.0;
        assert!(ff.has_high_penalty_types());
    }
    #[test]
    fn test_cgenff_no_penalty_ok() {
        let mut ff = CgenffForceField::new();
        ff.add_atom_type("KNOWN", "C", -0.1, 1.9, -0.07, 12.0);
        assert!(!ff.has_high_penalty_types());
    }
    #[test]
    fn test_cgenff_torsion_energy_zero_phase() {
        let e = CgenffForceField::cgenff_torsion_energy(0.0, &[(1.5, 1, 0.0)]);
        assert!((e - 3.0).abs() < 1e-12, "e={e}");
    }
    #[test]
    fn test_cgenff_torsion_energy_multiple_terms() {
        let terms = vec![(1.0, 1, 0.0), (0.5, 2, 0.0), (0.25, 3, 0.0)];
        let e = CgenffForceField::cgenff_torsion_energy(0.0, &terms);
        let expected = 2.0 * (1.0 + 0.5 + 0.25);
        assert!((e - expected).abs() < 1e-12, "e={e} expected={expected}");
    }
    #[test]
    fn test_cgenff_validate_valid() {
        let ff = CgenffForceField::new();
        let errors = validate_cgenff(&ff, 10);
        assert!(!has_critical_errors(&errors));
    }
    #[test]
    fn test_intra_decomp_default_zero() {
        let d = IntraEnergyDecomposition::default();
        assert_eq!(d.total(), 0.0);
        assert_eq!(d.bonded_total(), 0.0);
        assert_eq!(d.nonbonded14_total(), 0.0);
    }
    #[test]
    fn test_intra_decomp_total() {
        let d = IntraEnergyDecomposition {
            bond_energy: 1.0,
            angle_energy: 2.0,
            dihedral_energy: 3.0,
            improper_energy: 0.5,
            urey_bradley_energy: 0.25,
            lj14_energy: 0.1,
            coulomb14_energy: 0.15,
        };
        let expected = 1.0 + 2.0 + 3.0 + 0.5 + 0.25 + 0.1 + 0.15;
        assert!((d.total() - expected).abs() < 1e-12, "total={}", d.total());
    }
    #[test]
    fn test_intra_decomp_bonded_total() {
        let d = IntraEnergyDecomposition {
            bond_energy: 1.0,
            angle_energy: 2.0,
            dihedral_energy: 3.0,
            improper_energy: 0.5,
            urey_bradley_energy: 0.25,
            lj14_energy: 10.0,
            coulomb14_energy: 10.0,
        };
        assert!((d.bonded_total() - 6.75).abs() < 1e-12);
        assert!((d.nonbonded14_total() - 20.0).abs() < 1e-12);
    }
    #[test]
    fn test_intra_decomp_bond_fraction() {
        let d = IntraEnergyDecomposition {
            bond_energy: 5.0,
            angle_energy: 5.0,
            ..IntraEnergyDecomposition::default()
        };
        let frac = d.bond_fraction();
        assert!((frac - 0.5).abs() < 1e-10, "frac={frac}");
    }
    #[test]
    fn test_intra_decomp_zero_total_fraction() {
        let d = IntraEnergyDecomposition::default();
        assert_eq!(d.bond_fraction(), 0.0);
    }
    #[test]
    fn test_generate_cross_parameters_single_type() {
        let types = vec![LjAtomType {
            type_id: 0,
            sigma: 3.0,
            epsilon: 1.0,
        }];
        let params = generate_cross_parameters(&types, MixingRule::Geometric);
        assert_eq!(params.len(), 1, "Single type → 1 pair (self)");
        assert_eq!(params[0].0, 0);
        assert_eq!(params[0].1, 0);
    }
    #[test]
    fn test_generate_cross_parameters_three_types() {
        let types = vec![
            LjAtomType {
                type_id: 0,
                sigma: 3.0,
                epsilon: 1.0,
            },
            LjAtomType {
                type_id: 1,
                sigma: 3.5,
                epsilon: 0.5,
            },
            LjAtomType {
                type_id: 2,
                sigma: 4.0,
                epsilon: 2.0,
            },
        ];
        let params = generate_cross_parameters(&types, MixingRule::LorentzBerthelot);
        assert_eq!(params.len(), 6, "3 types → 6 cross pairs");
    }
    #[test]
    fn test_opls_dihedral_periodicity() {
        use std::f64::consts::PI;
        let e1 = OplsForceField::opls_dihedral_energy(0.5, 1.0, 0.5, 0.25, 0.1);
        let e2 = OplsForceField::opls_dihedral_energy(0.5 + 2.0 * PI, 1.0, 0.5, 0.25, 0.1);
        assert!((e1 - e2).abs() < 1e-10, "OPLS dihedral should be periodic");
    }
    #[test]
    fn test_charmm_urey_bradley_at_equilibrium() {
        let e = CharmmForceField::urey_bradley_energy(2.0, 100.0, 2.0);
        assert_eq!(e, 0.0, "At equilibrium, UB energy should be 0");
    }
    #[test]
    fn test_charmm_improper_at_equilibrium() {
        let e = CharmmForceField::improper_energy(0.5, 50.0, 0.5);
        assert_eq!(e, 0.0, "At equilibrium, improper energy should be 0");
    }
    #[test]
    fn test_forcefield_improper_at_equilibrium_zero() {
        let ri = [0.0_f64, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl = [1.0, -1.0, 0.0];
        let phi0 = dihedral_angle(ri, rj, rk, rl);
        let e = Forcefield::compute_improper_dihedral(ri, rj, rk, rl, 10.0, phi0);
        assert!(
            e.abs() < 1e-12,
            "At phi0, improper energy should be 0, got {e}"
        );
    }
    #[test]
    fn test_forcefield_improper_positive_when_distorted() {
        let ri = [0.0_f64, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl = [
            1.0,
            -std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
        ];
        let e = Forcefield::compute_improper_dihedral(ri, rj, rk, rl, 10.0, 0.0);
        assert!(
            e > 0.0,
            "Out-of-plane improper should have positive energy, got {e}"
        );
    }
    #[test]
    fn test_forcefield_improper_scales_with_k() {
        let ri = [0.0_f64, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        let rl = [
            1.0,
            -std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
        ];
        let e1 = Forcefield::compute_improper_dihedral(ri, rj, rk, rl, 10.0, 0.0);
        let e2 = Forcefield::compute_improper_dihedral(ri, rj, rk, rl, 20.0, 0.0);
        assert!(
            (e2 - 2.0 * e1).abs() < 1e-10,
            "Improper energy should scale linearly with k"
        );
    }
    #[test]
    fn test_forcefield_ub_at_equilibrium_zero() {
        let r1 = [0.0_f64, 0.0, 0.0];
        let r3 = [2.0, 0.0, 0.0];
        let e = Forcefield::compute_urey_bradley(r1, r3, 100.0, 2.0);
        assert_eq!(e, 0.0, "At r0, UB energy should be 0");
    }
    #[test]
    fn test_forcefield_ub_positive_when_stretched() {
        let r1 = [0.0_f64, 0.0, 0.0];
        let r3 = [2.5, 0.0, 0.0];
        let e = Forcefield::compute_urey_bradley(r1, r3, 100.0, 2.0);
        assert!(
            (e - 100.0 * 0.25).abs() < 1e-10,
            "UB energy should be k*(dr)^2 = 25.0, got {e}"
        );
    }
    #[test]
    fn test_forcefield_ub_symmetric_around_r0() {
        let r1 = [0.0_f64, 0.0, 0.0];
        let r3_plus = [2.3, 0.0, 0.0];
        let r3_minus = [1.7, 0.0, 0.0];
        let e_plus = Forcefield::compute_urey_bradley(r1, r3_plus, 100.0, 2.0);
        let e_minus = Forcefield::compute_urey_bradley(r1, r3_minus, 100.0, 2.0);
        assert!(
            (e_plus - e_minus).abs() < 1e-10,
            "UB energy should be symmetric: e+={e_plus}, e-={e_minus}"
        );
    }
    #[test]
    fn test_forcefield_cmap_uniform_grid_constant() {
        let n = 4usize;
        let c = 3.125_f64;
        let grid = vec![c; n * n];
        let e = Forcefield::compute_cmap_correction(0.0, 0.0, &grid, n);
        assert!(
            (e - c).abs() < 1e-12,
            "Uniform CMAP grid should return {c}, got {e}"
        );
    }
    #[test]
    fn test_forcefield_cmap_empty_grid_returns_zero() {
        let e = Forcefield::compute_cmap_correction(0.0, 0.0, &[], 0);
        assert_eq!(e, 0.0, "Empty CMAP grid should return 0");
    }
    #[test]
    fn test_forcefield_cmap_wrong_size_returns_zero() {
        let grid = vec![1.0_f64; 5];
        let e = Forcefield::compute_cmap_correction(0.0, 0.0, &grid, 4);
        assert_eq!(e, 0.0, "Wrong-size CMAP grid should return 0");
    }
    #[test]
    fn test_forcefield_cmap_bilinear_bounds() {
        let n = 3usize;
        let mut grid = vec![0.0_f64; n * n];
        grid[0] = 10.0;
        let e = Forcefield::compute_cmap_correction(-PI + 1e-10, -PI + 1e-10, &grid, n);
        assert!(
            e > 0.0,
            "CMAP near grid[0][0]=10 should return positive value, got {e}"
        );
    }
}
