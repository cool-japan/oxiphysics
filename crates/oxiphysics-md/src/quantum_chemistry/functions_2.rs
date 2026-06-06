//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::types::BasisFunction;
#[cfg(test)]
use super::types_2::HartreeFockSolver;

#[cfg(test)]
mod hf_tests {
    use super::*;
    #[test]
    fn h2_scf_energy_negative() {
        let solver =
            HartreeFockSolver::from_atoms(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0)
                .expect("H2 solver creation");
        let result = solver.run_scf();
        assert!(
            result.energy < 0.0,
            "H2 total energy should be negative: {}",
            result.energy
        );
    }
    #[test]
    fn h2_orbital_energies_sorted() {
        let solver =
            HartreeFockSolver::from_atoms(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0)
                .expect("H2");
        let result = solver.run_scf();
        let first = result.orbital_energies[0];
        let last = result.orbital_energies.last().copied().unwrap_or(0.0);
        assert!(
            first <= last,
            "orbital energies should be sorted ascending: first={first}, last={last}"
        );
    }
    #[test]
    fn h2_mulliken_neutral() {
        let solver =
            HartreeFockSolver::from_atoms(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0)
                .expect("H2");
        let result = solver.run_scf();
        let q0 = solver.mulliken_charge(&result.density_matrix, 1.0, &[0]);
        let q1 = solver.mulliken_charge(&result.density_matrix, 1.0, &[1]);
        assert!(
            (q0 + q1).abs() < 0.2,
            "total Mulliken charge not neutral: q0={q0}, q1={q1}"
        );
    }
    #[test]
    fn nuclear_attraction_sign() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.24);
        let v = bf.nuclear_attraction_with(&bf, 1.0, [0.0, 0.0, 0.0]);
        assert!(v < 0.0, "nuclear attraction should be negative: {v}");
    }
    #[test]
    fn eri_positive_coulomb() {
        let solver =
            HartreeFockSolver::from_atoms(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0)
                .expect("H2");
        let basis = vec![
            BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.24),
            BasisFunction::sto3g_s([1.398, 0.0, 0.0], 1.24),
        ];
        let eri = solver.compute_eris(&basis);
        let idx = solver.eri_index(0, 0, 0, 0);
        assert!(
            eri[idx] > 0.0,
            "(0,0|0,0) self-Coulomb must be positive: {}",
            eri[idx]
        );
    }
}
