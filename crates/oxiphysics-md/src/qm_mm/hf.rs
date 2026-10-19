// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin wrapper delegating `QmRegion::run()` for `QmMethod::Hf` to
//! the restricted Hartree-Fock solver in `crate::quantum_chemistry`.
//!
//! Positions are given in Å (the `QmRegion` convention) and converted to Bohr by
//! [`HartreeFockSolver::from_atoms`].  Energy is returned in kcal/mol; forces in
//! kcal/(mol·Å).
//!
//! External point charges are embedded in the one-electron Hamiltonian
//! ([`HartreeFockSolver::add_external_point_charges`]).  Because the SCF energy
//! is variational and the basis does not depend on the charge positions, the
//! forces on the external charges follow from the Hellmann–Feynman theorem:
//!
//! dE/dR_J = Tr(P ∂V^J/∂R_J) + Σ_A Z_A Q_J ∂(1/R_AJ)/∂R_J.

use super::QmRegion;
use super::embedding::{
    ANG_TO_BOHR, HARTREE_PER_BOHR_TO_KCAL_PER_ANG, HARTREE_TO_KCAL, begin_run,
    finite_difference_forces,
};
use super::error::QmError;
use crate::quantum_chemistry::{HartreeFockSolver, ScfStatus};
use oxiphysics_core::numerical_methods::cholesky;

/// Finite-difference step for the nuclear forces (Å).
const FD_STEP_ANG: f64 = 0.001;
/// Finite-difference step (Bohr) for the derivative of the external-potential
/// integrals with respect to a charge position.
const FD_STEP_CHARGE_BOHR: f64 = 1e-4;

/// Convert `(position Å, charge)` pairs to `(position Bohr, charge)`.
pub(crate) fn charges_to_bohr(charges: &[([f64; 3], f64)]) -> Vec<([f64; 3], f64)> {
    charges
        .iter()
        .map(|(p, q)| {
            (
                [p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR],
                *q,
            )
        })
        .collect()
}

/// First element of `atomic_numbers` without an STO-3G basis.
fn unsupported_element_error(atomic_numbers: &[u32]) -> QmError {
    atomic_numbers
        .iter()
        .find(|&&z| HartreeFockSolver::from_atoms(&[z], &[[0.0; 3]], 0).is_none())
        .map_or(QmError::EmptyRegion, |&z| QmError::UnsupportedElement(z))
}

/// Build a validated closed-shell STO-3G solver with the external charges
/// (Å, e) embedded in its core Hamiltonian.
pub(crate) fn build_solver(
    atomic_numbers: &[u32],
    positions_ang: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    external_charges_ang: &[([f64; 3], f64)],
) -> Result<HartreeFockSolver, QmError> {
    if multiplicity != 1 {
        return Err(QmError::OpenShellNotSupported);
    }
    let mut solver = HartreeFockSolver::from_atoms(atomic_numbers, positions_ang, charge)
        .ok_or_else(|| unsupported_element_error(atomic_numbers))?;
    if solver.n_electrons % 2 != 0 {
        return Err(QmError::OpenShellNotSupported);
    }
    if solver.n_electrons > 2 * solver.n_basis {
        return Err(QmError::InvalidElectronCount {
            electrons: solver.n_electrons as i64,
            n_basis: solver.n_basis,
        });
    }
    // The Roothaan-Hall problem needs a positive-definite overlap matrix.
    if cholesky(&solver.overlap).is_none() {
        return Err(QmError::EigensolverFailed);
    }
    solver.add_external_point_charges(&charges_to_bohr(external_charges_ang))?;
    Ok(solver)
}

/// Electrostatic embedding energy (Hartree) and forces on the external charges
/// (kcal/(mol·Å)) for a converged density `density` of `solver`.
///
/// E_emb = Tr(P V_ext) + Σ_A Σ_J Z_A Q_J / R_AJ.
pub(crate) fn embedding_terms(
    solver: &HartreeFockSolver,
    density: &[Vec<f64>],
) -> (f64, Vec<[f64; 3]>) {
    let trace = |v: &[Vec<f64>]| -> f64 {
        density
            .iter()
            .zip(v.iter())
            .map(|(p_row, v_row)| {
                p_row
                    .iter()
                    .zip(v_row.iter())
                    .map(|(p, x)| p * x)
                    .sum::<f64>()
            })
            .sum()
    };
    let mut e_emb = trace(&solver.external_potential_matrix(&solver.external_charges));
    let mut forces = Vec::with_capacity(solver.external_charges.len());
    for &(pos_j, q_j) in &solver.external_charges {
        let mut grad = [0.0_f64; 3];
        for (z_a, r_a) in solver
            .atomic_charges
            .iter()
            .zip(solver.atomic_positions.iter())
        {
            let d = [r_a[0] - pos_j[0], r_a[1] - pos_j[1], r_a[2] - pos_j[2]];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            e_emb += z_a * q_j / r;
            for (g, d_k) in grad.iter_mut().zip(d.iter()) {
                *g += z_a * q_j * d_k / (r * r * r);
            }
        }
        for (k, g) in grad.iter_mut().enumerate() {
            let mut plus = pos_j;
            plus[k] += FD_STEP_CHARGE_BOHR;
            let mut minus = pos_j;
            minus[k] -= FD_STEP_CHARGE_BOHR;
            let v_plus = solver.external_potential_matrix(&[(plus, q_j)]);
            let v_minus = solver.external_potential_matrix(&[(minus, q_j)]);
            *g += (trace(&v_plus) - trace(&v_minus)) / (2.0 * FD_STEP_CHARGE_BOHR);
        }
        forces.push([
            -grad[0] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
            -grad[1] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
            -grad[2] * HARTREE_PER_BOHR_TO_KCAL_PER_ANG,
        ]);
    }
    (e_emb, forces)
}

/// Mulliken net charges of every atom of `solver` for density `density`.
pub(crate) fn mulliken_charges(solver: &HartreeFockSolver, density: &[Vec<f64>]) -> Vec<f64> {
    solver
        .atomic_charges
        .iter()
        .enumerate()
        .map(|(a, &z)| {
            let indices: Vec<usize> = solver
                .basis_atoms
                .iter()
                .enumerate()
                .filter(|&(_, &atom)| atom == a)
                .map(|(mu, _)| mu)
                .collect();
            solver.mulliken_charge(density, z, &indices)
        })
        .collect()
}

/// Converged HF total energy (Hartree) of a geometry in the external field.
fn hf_energy(region: &QmRegion, positions_ang: &[[f64; 3]]) -> Result<f64, QmError> {
    let solver = build_solver(
        &region.atomic_numbers,
        positions_ang,
        region.charge,
        region.multiplicity,
        &region.external_charges,
    )?;
    let result = solver.run_scf();
    if result.status != ScfStatus::Converged {
        return Err(QmError::ScfNotConverged);
    }
    Ok(result.energy)
}

/// Run a restricted Hartree-Fock/STO-3G calculation on `region`.
///
/// Populates `region.energy`, `region.embedding_energy` (kcal/mol),
/// `region.forces`, `region.external_charge_forces` (kcal/mol/Å) and
/// `region.partial_charges` (Mulliken, e).  Nuclear forces are central finite
/// differences (±0.001 Å) of the SCF energy — O(6N) SCF calls, suitable only
/// for small QM regions.
///
/// # Errors
/// Empty or inconsistent region, element without an STO-3G basis, open shell,
/// singular overlap matrix, invalid external charge, or SCF non-convergence
/// (also at a displaced geometry).  The outputs of `region` are then left cleared.
pub fn run_hf(region: &mut QmRegion) -> Result<(), QmError> {
    begin_run(region)?;
    let solver = build_solver(
        &region.atomic_numbers,
        &region.positions,
        region.charge,
        region.multiplicity,
        &region.external_charges,
    )?;
    let result = solver.run_scf();
    if result.status != ScfStatus::Converged {
        return Err(QmError::ScfNotConverged);
    }
    let (e_emb, external_forces) = embedding_terms(&solver, &result.density_matrix);
    let charges = mulliken_charges(&solver, &result.density_matrix);
    let forces =
        finite_difference_forces(&region.positions, FD_STEP_ANG, |pos| hf_energy(region, pos))?;

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
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Hf, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        region
    }

    #[test]
    fn h2_hf_energy_negative() {
        let mut region = h2();
        run_hf(&mut region).expect("HF H2");
        assert!(
            region.energy < 0.0,
            "H2 HF energy should be negative: {}",
            region.energy
        );
        // E_HF(H2, 0.74 Å) ≈ −1.1167 Eh.
        assert!((region.energy / HARTREE_TO_KCAL + 1.1167).abs() < 1e-3);
        assert!((region.forces[0][0] + region.forces[1][0]).abs() < 1e-3);
    }

    #[test]
    fn unsupported_element_is_an_error() {
        let mut region = QmRegion::new(vec![0], QmMethod::Hf, 0, 1);
        region.atomic_numbers = vec![36]; // Kr — beyond the STO-3G table (H–Ar)
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(run_hf(&mut region), Err(QmError::UnsupportedElement(36)));
        assert_eq!(region.energy, 0.0);
    }

    #[test]
    fn open_shell_is_an_error() {
        let mut region = QmRegion::new(vec![0], QmMethod::Hf, 0, 1);
        region.atomic_numbers = vec![1];
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(run_hf(&mut region), Err(QmError::OpenShellNotSupported));
    }

    #[test]
    fn h2_is_polarised_by_external_charge() {
        let mut vac = h2();
        run_hf(&mut vac).expect("vacuum");
        let mut emb = h2();
        emb.external_charges = vec![([3.0, 0.0, 0.0], 1.0)];
        run_hf(&mut emb).expect("embedded");
        let q = &emb.partial_charges;
        assert!(q[1] < -1e-3 && q[0] > 1e-3, "charges {q:?}");
        assert!(emb.dipole_moment()[0] < -1e-3);
        assert!(emb.embedding_energy < 0.0);
        // The vacuum H₂ has no charge separation at all.
        assert!(vac.partial_charges.iter().all(|q| q.abs() < 1e-8));
    }

    #[test]
    fn zero_external_charge_reproduces_vacuum() {
        let mut vac = h2();
        run_hf(&mut vac).expect("vacuum");
        let mut emb = h2();
        emb.external_charges = vec![([0.3, 2.0, 0.0], 0.0)];
        run_hf(&mut emb).expect("embedded");
        assert!((vac.energy - emb.energy).abs() < 1e-8);
        assert!(emb.embedding_energy.abs() < 1e-12);
    }

    #[test]
    fn external_charge_force_matches_finite_difference() {
        let q_pos = [0.5, 1.8, -0.7];
        let mut region = h2();
        region.external_charges = vec![(q_pos, -0.7)];
        run_hf(&mut region).expect("embedded H2");
        let f = region.external_charge_forces[0];
        let h = 1e-3;
        for k in 0..3 {
            let mut e = [0.0; 2];
            for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
                let mut r = h2();
                let mut p = q_pos;
                p[k] += sign * h;
                r.external_charges = vec![(p, -0.7)];
                run_hf(&mut r).expect("displaced charge");
                *slot = r.total_energy_in_field();
            }
            let fd = -(e[0] - e[1]) / (2.0 * h);
            assert!(
                (fd - f[k]).abs() < 2e-3 + 1e-3 * fd.abs(),
                "component {k}: Hellmann-Feynman {} vs FD {fd}",
                f[k]
            );
        }
        for k in 0..3 {
            let total: f64 = region.forces.iter().map(|g| g[k]).sum::<f64>() + f[k];
            assert!(total.abs() < 5e-3, "net force component {k}: {total}");
        }
    }
}
