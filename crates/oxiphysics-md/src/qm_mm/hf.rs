// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin wrapper delegating `QmRegion::run()` for `QmMethod::Hf` to
//! the restricted Hartree-Fock solver in `crate::quantum_chemistry`.
//!
//! Positions are converted from Å (the `QmRegion` convention) to Bohr
//! internally.  Energy is returned in kcal/mol; forces in kcal/(mol·Å).

use super::QmRegion;
use crate::quantum_chemistry::HartreeFockSolver;

/// Run a restricted Hartree-Fock/STO-3G calculation on `region`.
///
/// Populates `region.energy` (kcal/mol) and `region.forces` (kcal/mol/Å).
/// Returns silently without modifying the region if any atom is unsupported
/// by the STO-3G basis (H, C, N, O, F, S are supported).
pub fn run_hf(region: &mut QmRegion) {
    let Some(solver) =
        HartreeFockSolver::from_atoms(&region.atomic_numbers, &region.positions, region.charge)
    else {
        return;
    };
    let result = solver.run_scf();
    // Hartree → kcal/mol
    region.energy = result.energy * 627.5094;
    // Numerical gradient via central differences
    region.forces = hf_numerical_forces(&region.atomic_numbers, &region.positions, region.charge);
}

/// Central-difference numerical forces for the HF energy.
///
/// O(6N) SCF calls — suitable only for small QM regions.
fn hf_numerical_forces(
    atomic_numbers: &[u32],
    positions: &[[f64; 3]],
    charge: i32,
) -> Vec<[f64; 3]> {
    const DELTA: f64 = 0.001; // Å
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for (i, force) in forces.iter_mut().enumerate() {
        for (k, f) in force.iter_mut().enumerate() {
            let e_plus = hf_energy_point(atomic_numbers, positions, charge, i, k, DELTA);
            let e_minus = hf_energy_point(atomic_numbers, positions, charge, i, k, -DELTA);
            *f = -(e_plus - e_minus) / (2.0 * DELTA);
        }
    }
    forces
}

fn hf_energy_point(
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
    solver.map(|s| s.run_scf().energy * 627.5094).unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qm_mm::{QmMethod, QmRegion};

    #[test]
    fn h2_hf_energy_negative() {
        let mut region = QmRegion::new(vec![0, 1], QmMethod::Hf, 0, 1);
        region.atomic_numbers = vec![1, 1];
        region.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        run_hf(&mut region);
        assert!(
            region.energy < 0.0,
            "H2 HF energy should be negative: {}",
            region.energy
        );
        assert!(region.energy.is_finite(), "H2 HF energy must be finite");
    }

    #[test]
    fn unsupported_element_no_panic() {
        let mut region = QmRegion::new(vec![0], QmMethod::Hf, 0, 1);
        region.atomic_numbers = vec![10]; // Ne — not supported
        region.set_positions(vec![[0.0, 0.0, 0.0]]);
        run_hf(&mut region); // must not panic; energy stays 0
    }
}
