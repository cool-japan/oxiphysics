// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin wrapper delegating `QmRegion::run()` for `QmMethod::DftB3lyp` to the
//! Kohn-Sham DFT solver in `crate::quantum_chemistry`.
//!
//! **Note:** `DftB3lyp` is implemented at the LDA/VWN level, not true B3LYP.
//! B3LYP requires exact-exchange HF integrals woven into the KS equations;
//! that is deferred to a future implementation pass.  The module documentation
//! clearly states the approximation level.

use super::QmRegion;
use crate::quantum_chemistry::{DensityFunctionalTheory, HartreeFockSolver, XcFunctional};

/// Run a Kohn-Sham LDA/VWN calculation on `region` (used for `QmMethod::DftB3lyp`).
///
/// **Approximation level:** LDA (Slater exchange + VWN correlation), not B3LYP.
/// Populates `region.energy` (kcal/mol) and `region.forces` (kcal/mol/Å).
pub fn run_dft_lda(region: &mut QmRegion) {
    let Some(hf_solver) =
        HartreeFockSolver::from_atoms(&region.atomic_numbers, &region.positions, region.charge)
    else {
        return;
    };
    let dft = DensityFunctionalTheory::new(hf_solver, XcFunctional::Lda);
    let result = dft.run_ks_scf();
    region.energy = result.energy * 627.5094;
    region.forces = dft_numerical_forces(&region.atomic_numbers, &region.positions, region.charge);
}

/// Central-difference numerical forces for the KS-DFT energy.
///
/// O(6N) SCF calls — suitable only for small QM regions.
fn dft_numerical_forces(
    atomic_numbers: &[u32],
    positions: &[[f64; 3]],
    charge: i32,
) -> Vec<[f64; 3]> {
    const DELTA: f64 = 0.001;
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for (i, force) in forces.iter_mut().enumerate() {
        for (k, f) in force.iter_mut().enumerate() {
            let e_plus = dft_energy_point(atomic_numbers, positions, charge, i, k, DELTA);
            let e_minus = dft_energy_point(atomic_numbers, positions, charge, i, k, -DELTA);
            *f = -(e_plus - e_minus) / (2.0 * DELTA);
        }
    }
    forces
}

fn dft_energy_point(
    atomic_numbers: &[u32],
    positions: &[[f64; 3]],
    charge: i32,
    atom: usize,
    coord: usize,
    delta: f64,
) -> f64 {
    let mut pos = positions.to_vec();
    pos[atom][coord] += delta;
    let solver = HartreeFockSolver::from_atoms(atomic_numbers, &pos, charge);
    solver
        .map(|s| {
            DensityFunctionalTheory::new(s, XcFunctional::Lda)
                .run_ks_scf()
                .energy
                * 627.5094
        })
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::{QmMethod, QmRegion};

    #[test]
    fn h2_dft_energy_negative() {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::DftB3lyp, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_dft_lda(&mut region);
        assert!(
            region.energy < 0.0,
            "H2 LDA energy should be negative: {}",
            region.energy
        );
        assert!(region.energy.is_finite(), "H2 LDA energy must be finite");
    }

    #[test]
    fn unsupported_element_no_panic() {
        let mut region = QmRegion::new(vec![0], QmMethod::DftB3lyp, 0, 1);
        region.atomic_numbers = vec![10];
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        run_dft_lda(&mut region);
    }
}
