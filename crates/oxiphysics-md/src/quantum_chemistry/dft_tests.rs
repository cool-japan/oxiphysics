// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the all-electron STO-3G Hartree–Fock and the grid-based
//! Kohn–Sham solvers against PySCF 2.14 (libxc) reference energies.

use super::error::QuantumChemistryError;
use super::functions::BOHR_TO_ANGSTROM;
use super::kohn_sham::xc_contribution;
use super::molecular_grid::{MolecularGrid, MolecularGridSettings};
use super::sto3g_basis::sto3g_atom_basis;
use super::types::BasisFunction;
use super::types_2::{BasisSetType, DensityFunctionalTheory, HartreeFockSolver};
use super::types_3::{ScfStatus, XcFunctional};
use super::xc_functional::evaluate_xc;

/// Crawford's H₂O geometry (bohr): O, H, H.
const WATER: [(u32, [f64; 3]); 3] = [
    (8, [0.0, -0.143_225_816_552, 0.0]),
    (1, [1.638_036_840_407, 1.136_548_822_547, 0.0]),
    (1, [-1.638_036_840_407, 1.136_548_822_547, 0.0]),
];

/// STO-3G solver for atoms given in bohr.
fn sto3g_solver_bohr(atoms: &[(u32, [f64; 3])]) -> HartreeFockSolver {
    let mut basis = Vec::new();
    let mut basis_atoms = Vec::new();
    for (a, (z, pos)) in atoms.iter().enumerate() {
        let shells = sto3g_atom_basis(*z, *pos).expect("STO-3G element");
        basis_atoms.extend(std::iter::repeat_n(a, shells.len()));
        basis.extend(shells);
    }
    let charges: Vec<f64> = atoms.iter().map(|(z, _)| f64::from(*z)).collect();
    let n_electrons = atoms.iter().map(|(z, _)| *z as usize).sum();
    let mut solver = HartreeFockSolver::from_basis(
        basis,
        basis_atoms,
        charges,
        atoms.iter().map(|(_, p)| *p).collect(),
        n_electrons,
        BasisSetType::Sto3G,
    )
    .expect("solver");
    solver.energy_threshold = 1e-10;
    solver
}

/// He with 14 uncontracted even-tempered s functions α_k = 0.05·2.5^k.
fn helium_even_tempered() -> HartreeFockSolver {
    let basis: Vec<BasisFunction> = (0..14)
        .map(|k| BasisFunction::contracted([0.0; 3], [0, 0, 0], &[0.05 * 2.5_f64.powi(k)], &[1.0]))
        .collect();
    let n = basis.len();
    let mut solver = HartreeFockSolver::from_basis(
        basis,
        vec![0; n],
        vec![2.0],
        vec![[0.0; 3]],
        2,
        BasisSetType::Custom,
    )
    .expect("He solver");
    solver.energy_threshold = 1e-10;
    solver
}

fn h2() -> HartreeFockSolver {
    sto3g_solver_bohr(&[(1, [0.0; 3]), (1, [0.0, 0.0, 1.4])])
}

fn ks_energy(solver: HartreeFockSolver, xc: XcFunctional) -> f64 {
    let n_electrons = solver.n_electrons as f64;
    let result = DensityFunctionalTheory::new(solver, xc)
        .run_ks_scf()
        .expect("KS-SCF");
    assert_eq!(
        result.status,
        ScfStatus::Converged,
        "{xc:?} did not converge"
    );
    assert!(
        (result.integrated_electrons - n_electrons).abs() < 1e-5,
        "∫ρ = {}",
        result.integrated_electrons
    );
    result.energy
}

#[test]
fn water_rhf_sto3g_matches_crawford_reference() {
    // Crawford programming project #3 / PySCF: −74.942079928192 Eh.
    let solver = sto3g_solver_bohr(&WATER);
    assert_eq!(solver.n_basis, 7);
    assert_eq!(solver.n_electrons, 10);
    assert!((solver.nuclear_repulsion() - 8.002_367_061_810_77).abs() < 1e-9);
    let result = solver.run_scf();
    assert_eq!(result.status, ScfStatus::Converged);
    assert!(
        (result.energy - (-74.942_079_928_192_4)).abs() < 1e-6,
        "E(RHF) = {}",
        result.energy
    );
}

#[test]
fn water_from_atoms_at_szabo_geometry() {
    // R(OH) = 1.809 bohr, ∠HOH = 104.52°: Szabo & Ostlund −74.963, PySCF −74.962940028.
    let r = 1.809 * BOHR_TO_ANGSTROM;
    let theta = 104.52_f64.to_radians();
    let solver = HartreeFockSolver::from_atoms(
        &[8, 1, 1],
        &[
            [0.0, 0.0, 0.0],
            [r, 0.0, 0.0],
            [r * theta.cos(), r * theta.sin(), 0.0],
        ],
        0,
    )
    .expect("H2O");
    assert_eq!(solver.n_electrons, 10, "all electrons, not valence only");
    let mut solver = solver;
    solver.energy_threshold = 1e-10;
    let result = solver.run_scf();
    assert_eq!(result.status, ScfStatus::Converged);
    assert!(
        (result.energy - (-74.962_940_028_256_74)).abs() < 1e-5,
        "E(RHF) = {}",
        result.energy
    );
}

#[test]
fn hydrogen_sulfide_rhf_uses_the_3sp_shell() {
    // PySCF RHF/STO-3G: −394.305545941049 Eh.
    let solver = sto3g_solver_bohr(&[(16, [0.0; 3]), (1, [0.0, 0.0, 2.5]), (1, [2.4, 0.0, -0.6])]);
    assert_eq!(solver.n_basis, 11);
    let result = solver.run_scf();
    assert_eq!(result.status, ScfStatus::Converged);
    assert!(
        (result.energy - (-394.305_545_941_049)).abs() < 1e-5,
        "E(RHF) = {}",
        result.energy
    );
}

#[test]
fn helium_hf_limit_with_even_tempered_basis() {
    // PySCF: −2.8616763918548 (HF limit −2.86168).
    let result = helium_even_tempered().run_scf();
    assert_eq!(result.status, ScfStatus::Converged);
    assert!((result.energy - (-2.861_676_391_854_796)).abs() < 1e-7);
}

#[test]
fn helium_lda_and_pbe_match_pyscf() {
    // PySCF (200×1202 grid): LDA (Slater+VWN5) −2.8348326042, PBE −2.8929308947.
    // (NIST LDA basis-set limit: −2.834836.)
    let lda = ks_energy(helium_even_tempered(), XcFunctional::Lda);
    assert!((lda - (-2.834_832_604_225_022_5)).abs() < 1e-5, "LDA {lda}");
    let pbe = ks_energy(helium_even_tempered(), XcFunctional::Pbe);
    assert!((pbe - (-2.892_930_894_652_247)).abs() < 1e-5, "PBE {pbe}");
}

#[test]
fn helium_b3lyp_and_tpss_match_pyscf() {
    let b3lyp = ks_energy(helium_even_tempered(), XcFunctional::B3lyp);
    assert!(
        (b3lyp - (-2.915_214_879_665_825_3)).abs() < 1e-5,
        "B3LYP {b3lyp}"
    );
    let tpss = ks_energy(helium_even_tempered(), XcFunctional::Tpss);
    assert!(
        (tpss - (-2.909_653_510_889_225)).abs() < 1e-5,
        "TPSS {tpss}"
    );
}

#[test]
fn hydrogen_molecule_ks_energies_match_pyscf() {
    let cases = [
        (XcFunctional::Lda, -1.121_200_704_284_351_3),
        (XcFunctional::Pbe, -1.152_064_373_123_530_5),
        (XcFunctional::B3lyp, -1.165_400_928_155_918_3),
        (XcFunctional::Tpss, -1.166_754_551_603_303),
    ];
    for (xc, reference) in cases {
        let e = ks_energy(h2(), xc);
        assert!((e - reference).abs() < 1e-5, "{xc:?}: {e} vs {reference}");
    }
}

#[test]
fn water_ks_energies_match_pyscf() {
    let cases = [
        (XcFunctional::Lda, -74.733_272_086_109_04),
        (XcFunctional::Pbe, -75.229_898_449_471_46),
        (XcFunctional::B3lyp, -75.312_291_140_686_31),
        (XcFunctional::Tpss, -75.331_663_461_672_66),
    ];
    for (xc, reference) in cases {
        let e = ks_energy(sto3g_solver_bohr(&WATER), xc);
        assert!((e - reference).abs() < 2e-5, "{xc:?}: {e} vs {reference}");
    }
}

#[test]
fn ks_result_components_are_consistent() {
    let solver = h2();
    let t = solver.kinetic_matrix();
    let dft = DensityFunctionalTheory::new(solver, XcFunctional::Lda);
    let r = dft.run_ks_scf().expect("KS");
    // T = Tr(P T) is positive and the virial ratio −V/T is close to 2 near equilibrium.
    let tr: f64 = r
        .density_matrix
        .iter()
        .zip(t.iter())
        .map(|(p, t)| p.iter().zip(t.iter()).map(|(a, b)| a * b).sum::<f64>())
        .sum();
    assert!((r.kinetic_energy - tr).abs() < 1e-12);
    assert!(r.kinetic_energy > 0.0);
    let virial = -(r.energy - r.kinetic_energy) / r.kinetic_energy;
    assert!((virial - 2.0).abs() < 0.2, "virial ratio {virial}");
    assert!(r.xc_energy < 0.0);
    // Orbital energies ascending, MO coefficients one row per orbital.
    assert!(r.orbital_energies.windows(2).all(|w| w[0] <= w[1]));
    assert_eq!(r.mo_coefficients.len(), 2);
    assert!(r.iterations >= 2);
}

#[test]
fn hybrid_exact_exchange_enters_the_scf() {
    // Changing a₀ for B3LYP must change the energy (exact exchange is in the Fock matrix).
    let base = DensityFunctionalTheory::new(h2(), XcFunctional::B3lyp)
        .run_ks_scf()
        .expect("B3LYP");
    let mut no_hf = DensityFunctionalTheory::new(h2(), XcFunctional::B3lyp);
    no_hf.alpha_x = 0.0;
    let pure = no_hf.run_ks_scf().expect("B3LYP a0=0");
    assert!((base.energy - pure.energy).abs() > 1e-3);
}

#[test]
fn xc_potential_matrix_is_the_density_derivative_of_exc() {
    // dE_xc/dλ at P + λD equals Tr(V_xc D) for a symmetric perturbation D.
    let solver = sto3g_solver_bohr(&WATER);
    let hf = solver.run_scf();
    let grid = MolecularGrid::new(
        &solver.atomic_charges,
        &solver.atomic_positions,
        MolecularGridSettings {
            radial_points: 40,
            polar_points: 8,
        },
    )
    .expect("grid");
    let n = solver.n_basis;
    // Symmetric in (i, j) by construction.
    let mut d = vec![vec![0.0_f64; n]; n];
    for (i, row) in d.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            *x = 0.01 * (((i * 7 + j * 7 + i * j) % 5) as f64 - 2.0);
        }
    }
    for xc in [
        XcFunctional::Lda,
        XcFunctional::Pbe,
        XcFunctional::B3lyp,
        XcFunctional::Tpss,
    ] {
        let a0 = if xc == XcFunctional::B3lyp { 0.2 } else { 0.0 };
        let p = &hf.density_matrix;
        let shifted = |h: f64| -> Vec<Vec<f64>> {
            p.iter()
                .zip(d.iter())
                .map(|(pr, dr)| pr.iter().zip(dr.iter()).map(|(x, y)| x + h * y).collect())
                .collect()
        };
        let h = 1e-5;
        let e_plus = xc_contribution(&grid, &solver.basis_functions, &shifted(h), xc, a0).energy;
        let e_minus = xc_contribution(&grid, &solver.basis_functions, &shifted(-h), xc, a0).energy;
        let v = xc_contribution(&grid, &solver.basis_functions, p, xc, a0).matrix;
        let analytic: f64 = v
            .iter()
            .zip(d.iter())
            .map(|(vr, dr)| vr.iter().zip(dr.iter()).map(|(x, y)| x * y).sum::<f64>())
            .sum();
        let numeric = (e_plus - e_minus) / (2.0 * h);
        assert!(
            (analytic - numeric).abs() < 1e-6 * analytic.abs().max(1.0),
            "{xc:?}: analytic {analytic}, numeric {numeric}"
        );
    }
}

#[test]
fn ks_requires_basis_functions_and_closed_shell() {
    // Integrals supplied directly (no basis functions): no grid can be evaluated.
    let bare = HartreeFockSolver::new(
        1,
        2,
        BasisSetType::Sto3G,
        vec![vec![1.0]],
        vec![vec![-1.0]],
        vec![0.5],
        vec![1.0],
        vec![[0.0; 3]],
    );
    assert_eq!(
        DensityFunctionalTheory::new(bare, XcFunctional::Lda)
            .run_ks_scf()
            .err(),
        Some(QuantumChemistryError::MissingData("basis functions"))
    );
    let mut odd = h2();
    odd.n_electrons = 1;
    assert_eq!(
        DensityFunctionalTheory::new(odd, XcFunctional::Lda)
            .run_ks_scf()
            .err(),
        Some(QuantumChemistryError::OpenShellReference(1))
    );
}

#[test]
fn pseudopotential_shifts_the_total_energy() {
    let plain = DensityFunctionalTheory::new(h2(), XcFunctional::Lda)
        .run_ks_scf()
        .expect("LDA");
    let mut shifted = DensityFunctionalTheory::new(h2(), XcFunctional::Lda);
    shifted.set_pseudopotential(-0.1);
    let r = shifted.run_ks_scf().expect("LDA + PP");
    assert!((r.energy - plain.energy + 0.1).abs() < 1e-9);
}

#[test]
fn vwn5_correlation_matches_libxc() {
    for (rho, eps) in [
        (0.1, -0.053_397_289_185_949_81),
        (1.0, -0.071_592_612_306_790_65),
        (0.01, -0.037_645_190_262_171_42),
    ] {
        let got = DensityFunctionalTheory::vwn_correlation_energy_density(rho);
        assert!((got - eps).abs() < 1e-12, "ρ = {rho}: {got} vs {eps}");
    }
}

#[test]
fn uniform_gas_xc_energy_uses_the_selected_functional() {
    let rho = 0.3;
    let volume = 2.0;
    let lda = DensityFunctionalTheory::new(h2(), XcFunctional::Lda);
    let expected_lda = (DensityFunctionalTheory::lda_exchange_energy_density(rho)
        + DensityFunctionalTheory::vwn_correlation_energy_density(rho))
        * rho
        * volume;
    assert!((lda.xc_energy_uniform(rho, volume) - expected_lda).abs() < 1e-12);
    // PBE reduces to Slater + PW92 for a uniform gas (not VWN).
    let pbe = DensityFunctionalTheory::new(h2(), XcFunctional::Pbe);
    let pbe_e = pbe.xc_energy_uniform(rho, volume);
    let pbe_point = evaluate_xc(XcFunctional::Pbe, 0.0, rho, 0.0, 0.0);
    assert!((pbe_e - pbe_point.energy_density * volume).abs() < 1e-12);
    assert!((pbe_e - expected_lda).abs() > 1e-6);
    assert!((pbe_e - expected_lda).abs() < 1e-3);
    // B3LYP keeps LYP correlation; its exact exchange counts as Slater exchange.
    let b3 = DensityFunctionalTheory::new(h2(), XcFunctional::B3lyp);
    assert!((b3.xc_energy_uniform(rho, volume) - expected_lda).abs() > 1e-3);
}
