// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin wrapper delegating `QmRegion::run()` for `QmMethod::DftB3lyp` to the
//! restricted Kohn–Sham B3LYP/STO-3G solver in `crate::quantum_chemistry`
//! (grid-based V_xc with B88 + LYP + VWN(RPA) and 20 % exact exchange inside
//! the SCF).
//!
//! External point charges are embedded in the core Hamiltonian shared with
//! the HF solver; the embedding energy and the forces on the charges use the
//! Kohn-Sham density with the Hellmann–Feynman expression of
//! [`super::hf`], which is exact for a variational KS energy.

use super::QmRegion;
use super::embedding::{HARTREE_TO_KCAL, begin_run, finite_difference_forces};
use super::error::QmError;
use super::hf::{build_solver, embedding_terms, mulliken_charges};
use crate::quantum_chemistry::{DensityFunctionalTheory, ScfStatus, XcFunctional};

/// Functional used for `QmMethod::DftB3lyp`.
const FUNCTIONAL: XcFunctional = XcFunctional::B3lyp;

/// Finite-difference step for the nuclear forces (Å).
const FD_STEP_ANG: f64 = 0.001;

/// Converged KS-B3LYP energy (Hartree) of a geometry in the external field.
fn dft_energy(region: &QmRegion, positions_ang: &[[f64; 3]]) -> Result<f64, QmError> {
    let solver = build_solver(
        &region.atomic_numbers,
        positions_ang,
        region.charge,
        region.multiplicity,
        &region.external_charges,
    )?;
    let result = DensityFunctionalTheory::new(solver, FUNCTIONAL).run_ks_scf()?;
    if result.status != ScfStatus::Converged {
        return Err(QmError::ScfNotConverged);
    }
    Ok(result.energy)
}

/// Run a restricted Kohn–Sham B3LYP/STO-3G calculation on `region`
/// (used for `QmMethod::DftB3lyp`).
///
/// Populates `region.energy`, `region.embedding_energy` (kcal/mol),
/// `region.forces`, `region.external_charge_forces` (kcal/mol/Å) and
/// `region.partial_charges` (Mulliken, e).
///
/// # Errors
/// Empty or inconsistent region, element without an STO-3G basis, open shell,
/// singular overlap matrix, invalid external charge, or SCF non-convergence
/// (also at a displaced geometry).  The outputs of `region` are then left cleared.
pub fn run_dft_b3lyp(region: &mut QmRegion) -> Result<(), QmError> {
    begin_run(region)?;
    let solver = build_solver(
        &region.atomic_numbers,
        &region.positions,
        region.charge,
        region.multiplicity,
        &region.external_charges,
    )?;
    let result = DensityFunctionalTheory::new(solver.clone(), FUNCTIONAL).run_ks_scf()?;
    if result.status != ScfStatus::Converged {
        return Err(QmError::ScfNotConverged);
    }
    let (e_emb, external_forces) = embedding_terms(&solver, &result.density_matrix);
    let charges = mulliken_charges(&solver, &result.density_matrix);
    let forces = finite_difference_forces(&region.positions, FD_STEP_ANG, |pos| {
        dft_energy(region, pos)
    })?;

    region.energy = (result.energy - e_emb) * HARTREE_TO_KCAL;
    region.embedding_energy = e_emb * HARTREE_TO_KCAL;
    region.forces = forces;
    region.partial_charges = charges;
    region.external_charge_forces = external_forces;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::{QmMethod, QmRegion};

    fn h2() -> QmRegion {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::DftB3lyp, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        region
    }

    #[test]
    fn h2_dft_energy_negative() {
        let mut region = h2();
        run_dft_b3lyp(&mut region).expect("B3LYP H2");
        assert!(
            region.energy < 0.0,
            "H2 B3LYP energy should be negative: {}",
            region.energy
        );
        assert!(region.energy.is_finite(), "H2 B3LYP energy must be finite");
    }

    #[test]
    fn unsupported_element_is_an_error() {
        let mut region = QmRegion::new(vec![0], QmMethod::DftB3lyp, 0, 1);
        region.atomic_numbers = vec![36];
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(
            run_dft_b3lyp(&mut region),
            Err(QmError::UnsupportedElement(36))
        );
    }

    #[test]
    fn external_charges_polarise_the_density() {
        let mut region = h2();
        region.external_charges = vec![([3.0, 0.0, 0.0], 1.0)];
        run_dft_b3lyp(&mut region).expect("embedded B3LYP H2");
        let q = &region.partial_charges;
        assert!(q[1] < -1e-3 && q[0] > 1e-3, "charges {q:?}");
        assert!(region.embedding_energy < 0.0);
        assert_eq!(region.external_charge_forces.len(), 1);
        // The positive charge is attracted towards the polarised molecule.
        assert!(region.external_charge_forces[0][0] < 0.0);
    }
}
