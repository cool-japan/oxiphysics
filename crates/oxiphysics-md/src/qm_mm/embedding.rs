// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared helpers for electrostatic embedding of a QM region in external
//! (MM) point charges, and unit conversions used by every QM engine.
//!
//! Positions in a [`QmRegion`] are in Å, charges in units of e, energies in
//! kcal/mol and forces in kcal/(mol·Å).  The QM engines work in atomic units, so
//! the external electrostatic potential at a nucleus is reported in Hartree/e:
//!
//! φ_A = Σ_J Q_J / R_AJ   (R_AJ in Bohr)

use super::{QmError, QmRegion};

/// Ångström → Bohr.
pub(crate) const ANG_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;
/// Hartree → kcal/mol.
pub(crate) const HARTREE_TO_KCAL: f64 = 627.509_474_063_1;
/// Hartree/Bohr → kcal/(mol·Å).
pub(crate) const HARTREE_PER_BOHR_TO_KCAL_PER_ANG: f64 = HARTREE_TO_KCAL * ANG_TO_BOHR;
/// Minimum allowed distance (Å) between an external charge and a QM nucleus.
const MIN_CHARGE_DISTANCE_ANG: f64 = 1e-6;

/// Validate the geometric input of `region` and reset all of its outputs.
///
/// Called at the start of every engine so that a failed calculation never leaves
/// the energy / forces of a previous step behind.  Returns the external
/// electrostatic potential φ_A (Hartree/e) at every QM nucleus.
pub(crate) fn begin_run(region: &mut QmRegion) -> Result<Vec<f64>, QmError> {
    let n = region.atomic_numbers.len();
    region.energy = 0.0;
    region.embedding_energy = 0.0;
    region.forces = vec![[0.0; 3]; n];
    region.partial_charges = vec![0.0; n];
    region.external_charge_forces = vec![[0.0; 3]; region.external_charges.len()];
    if n == 0 {
        return Err(QmError::EmptyRegion);
    }
    if region.positions.len() != n {
        return Err(QmError::DimensionMismatch {
            what: "QM positions vs atomic numbers",
            expected: n,
            found: region.positions.len(),
        });
    }
    external_potential(&region.positions, &region.external_charges)
}

/// External electrostatic potential φ_A = Σ_J Q_J / R_AJ (Hartree/e) at each QM
/// nucleus.
///
/// Fails if a charge is non-finite or (numerically) coincides with a nucleus.
pub(crate) fn external_potential(
    positions_ang: &[[f64; 3]],
    external_charges: &[([f64; 3], f64)],
) -> Result<Vec<f64>, QmError> {
    let mut phi = vec![0.0_f64; positions_ang.len()];
    for (index, (pos_j, q_j)) in external_charges.iter().enumerate() {
        if !q_j.is_finite() || pos_j.iter().any(|x| !x.is_finite()) {
            return Err(QmError::InvalidExternalCharge { index });
        }
        for (phi_a, pos_a) in phi.iter_mut().zip(positions_ang.iter()) {
            let r_ang = distance(pos_a, pos_j);
            if r_ang < MIN_CHARGE_DISTANCE_ANG {
                return Err(QmError::InvalidExternalCharge { index });
            }
            *phi_a += q_j / (r_ang * ANG_TO_BOHR);
        }
    }
    Ok(phi)
}

/// Forces (kcal/(mol·Å)) on the external charges from the derivatives of the
/// QM energy with respect to the nuclear potentials, dE/dφ_A (Hartree per
/// Hartree/e, i.e. units of e).
///
/// Valid whenever the external charges enter the QM energy only through
/// φ_A (PM3 and SCC-DFTB):
///
/// F_J = −Σ_A (∂E/∂φ_A) ∂φ_A/∂R_J,   ∂φ_A/∂R_J = Q_J (R_A − R_J)/R_AJ³.
pub(crate) fn forces_from_potential_derivatives(
    positions_ang: &[[f64; 3]],
    external_charges: &[([f64; 3], f64)],
    de_dphi: &[f64],
) -> Vec<[f64; 3]> {
    external_charges
        .iter()
        .map(|(pos_j, q_j)| {
            let mut f = [0.0_f64; 3];
            for (pos_a, &dedphi) in positions_ang.iter().zip(de_dphi.iter()) {
                let d_bohr = [
                    (pos_a[0] - pos_j[0]) * ANG_TO_BOHR,
                    (pos_a[1] - pos_j[1]) * ANG_TO_BOHR,
                    (pos_a[2] - pos_j[2]) * ANG_TO_BOHR,
                ];
                let r2 = d_bohr[0] * d_bohr[0] + d_bohr[1] * d_bohr[1] + d_bohr[2] * d_bohr[2];
                let r3 = r2 * r2.sqrt();
                for (f_k, d_k) in f.iter_mut().zip(d_bohr.iter()) {
                    *f_k -= dedphi * q_j * d_k / r3;
                }
            }
            [
                f[0] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
                f[1] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
                f[2] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
            ]
        })
        .collect()
}

/// Euclidean distance between two points.
pub(crate) fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Central finite-difference forces (kcal/(mol·Å)) on the QM nuclei of an energy
/// function `energy(positions) -> Hartree`, with step `delta_ang` (Å).
///
/// Any failure of a displaced single point is propagated instead of being
/// replaced by a default energy.
pub(crate) fn finite_difference_forces<F>(
    positions_ang: &[[f64; 3]],
    delta_ang: f64,
    mut energy: F,
) -> Result<Vec<[f64; 3]>, QmError>
where
    F: FnMut(&[[f64; 3]]) -> Result<f64, QmError>,
{
    let mut forces = vec![[0.0_f64; 3]; positions_ang.len()];
    let mut displaced = positions_ang.to_vec();
    for (i, f_i) in forces.iter_mut().enumerate() {
        for (k, f_ik) in f_i.iter_mut().enumerate() {
            let original = displaced[i][k];
            displaced[i][k] = original + delta_ang;
            let e_plus = energy(&displaced)?;
            displaced[i][k] = original - delta_ang;
            let e_minus = energy(&displaced)?;
            displaced[i][k] = original;
            *f_ik = -(e_plus - e_minus) / (2.0 * delta_ang) * HARTREE_TO_KCAL;
        }
    }
    Ok(forces)
}
