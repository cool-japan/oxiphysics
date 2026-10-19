// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the grid Kohn-Sham DFT driver.

use super::dft_md::{smeared_coulomb, smeared_coulomb_derivative_over_u};
use super::real_space_grid::erf;
use super::{DensityFunctionalMd, QuantumChemistryMdError, UniformGrid, XcFunctional};

const PROTON_MASS: f64 = 1_836.152_673;

/// Gaussian density exp(−α|r − c|²) normalised on the grid to `n_el` electrons.
fn gaussian_density(grid: &UniformGrid, centre: [f64; 3], alpha: f64, n_el: f64) -> Vec<f64> {
    let raw: Vec<f64> = grid
        .points()
        .iter()
        .map(|p| {
            let d2 = (p[0] - centre[0]).powi(2)
                + (p[1] - centre[1]).powi(2)
                + (p[2] - centre[2]).powi(2);
            (-alpha * d2).exp()
        })
        .collect();
    let total: f64 = raw.iter().sum::<f64>() * grid.dv();
    raw.iter().map(|x| x * n_el / total).collect()
}

fn helium(functional: XcFunctional) -> DensityFunctionalMd {
    let grid = UniformGrid::centered([0.0; 3], 0.25, [20, 20, 20]).expect("grid");
    let density = gaussian_density(&grid, [0.0; 3], 1.5, 2.0);
    DensityFunctionalMd::new(functional, grid, density)
        .expect("dft")
        .with_ions(vec![2.0], vec![[0.0; 3]], vec![4.0 * PROTON_MASS])
        .expect("ions")
}

fn hydrogen_molecule(functional: XcFunctional, half_bond: f64) -> DensityFunctionalMd {
    let grid = UniformGrid::centered([0.0; 3], 0.25, [24, 20, 20]).expect("grid");
    let density: Vec<f64> = gaussian_density(&grid, [-half_bond, 0.0, 0.0], 1.0, 1.0)
        .iter()
        .zip(gaussian_density(&grid, [half_bond, 0.0, 0.0], 1.0, 1.0).iter())
        .map(|(a, b)| a + b)
        .collect();
    DensityFunctionalMd::new(functional, grid, density)
        .expect("dft")
        .with_ions(
            vec![1.0, 1.0],
            vec![[-half_bond, 0.0, 0.0], [half_bond, 0.0, 0.0]],
            vec![PROTON_MASS, PROTON_MASS],
        )
        .expect("ions")
}

#[test]
fn test_smeared_coulomb_derivative_is_consistent() {
    let sigma = 0.4;
    for &u in &[1e-4_f64, 3e-3, 0.05, 0.3, 0.9, 2.5, 8.0] {
        let h = (1e-4 * u).max(1e-6);
        let fd = (smeared_coulomb(u + h, sigma) - smeared_coulomb(u - h, sigma)) / (2.0 * h);
        let analytic = smeared_coulomb_derivative_over_u(u, sigma) * u;
        assert!(
            (fd - analytic).abs() < 1e-6 * analytic.abs().max(1e-3),
            "u {u}: {fd} vs {analytic}"
        );
    }
    // Far from the Gaussian nucleus the potential is the bare Coulomb 1/u.
    assert!((smeared_coulomb(10.0, sigma) - 0.1).abs() < 1e-15);
    assert!(
        (smeared_coulomb(0.0, sigma) - 2.0 / (sigma * std::f64::consts::PI.sqrt())).abs() < 1e-15
    );
}

#[test]
fn test_dft_symmetric_density_gives_zero_force() {
    let grid = UniformGrid::centered([0.0; 3], 0.3, [16, 16, 16]).expect("grid");
    let density = gaussian_density(&grid, [0.0; 3], 1.0, 2.0);
    let dft = DensityFunctionalMd::new(XcFunctional::Lda, grid, density)
        .expect("dft")
        .with_ions(vec![2.0], vec![[0.0; 3]], vec![1.0])
        .expect("ions");
    let f = dft.forces().expect("forces");
    for c in f[0] {
        assert!(c.abs() < 1e-10, "{c}");
    }
}

#[test]
fn test_dft_off_centre_density_attracts_nucleus_with_correct_magnitude() {
    // Gaussian nucleus (width σ, charge Z) and Gaussian electron cloud
    // (exponent α, N electrons) at distance d: E = −Z N erf(d/s)/d with
    // s = √(σ² + 1/α); the nucleus is pulled towards the cloud.
    let grid = UniformGrid::centered([0.0; 3], 0.2, [28, 24, 24]).expect("grid");
    let (alpha, d, z, n_el) = (2.0_f64, 1.0_f64, 2.0_f64, 2.0_f64);
    let density = gaussian_density(&grid, [d, 0.0, 0.0], alpha, n_el);
    let mut dft = DensityFunctionalMd::new(XcFunctional::Lda, grid, density)
        .expect("dft")
        .with_ions(vec![z], vec![[0.0; 3]], vec![1.0])
        .expect("ions");
    dft.nuclear_smearing = 0.3;
    let f = dft.forces().expect("forces");
    let s = (0.3_f64.powi(2) + 1.0 / alpha).sqrt();
    let phi_prime = (2.0 / (s * std::f64::consts::PI.sqrt()) * (-(d / s).powi(2)).exp() * d
        - erf(d / s))
        / (d * d);
    let expected = -z * n_el * phi_prime;
    assert!(f[0][0] > 0.0);
    assert!(
        ((f[0][0] - expected) / expected).abs() < 5e-3,
        "{} vs {expected}",
        f[0][0]
    );
    assert!(f[0][1].abs() < 1e-10 && f[0][2].abs() < 1e-10);
}

#[test]
fn test_dft_dipole_uses_volume_element_and_electron_charge() {
    let grid = UniformGrid::centered([0.0; 3], 0.25, [24, 16, 16]).expect("grid");
    let density = gaussian_density(&grid, [1.0, 0.0, 0.0], 2.0, 2.0);
    let dft = DensityFunctionalMd::new(XcFunctional::Lda, grid, density)
        .expect("dft")
        .with_ions(vec![2.0], vec![[1.0, 0.0, 0.0]], vec![1.0])
        .expect("ions");
    let d = dft.dipole_moment();
    // −2 e × 1 bohr = −5.0835 D.
    assert!((d[0] + 2.0 * 2.541_746_473).abs() < 1e-3, "{d:?}");
    assert!(d[1].abs() < 1e-10 && d[2].abs() < 1e-10);
    // A neutral atom centred on its own electron cloud has no total dipole.
    let total = dft.total_dipole_moment();
    assert!(total[0].abs() < 1e-3, "{total:?}");
    assert!((dft.integrated_density() - 2.0).abs() < 1e-12);
}

#[test]
fn test_dft_rejects_inconsistent_inputs() {
    let grid = UniformGrid::centered([0.0; 3], 0.3, [8, 8, 8]).expect("grid");
    let err = DensityFunctionalMd::new(XcFunctional::Lda, grid.clone(), vec![0.0; 10]);
    assert!(matches!(
        err,
        Err(QuantumChemistryMdError::DimensionMismatch { .. })
    ));
    let mut hydrogen_atom = DensityFunctionalMd::new(XcFunctional::Lda, grid, vec![0.0; 512])
        .expect("dft")
        .with_ions(vec![1.0], vec![[0.0; 3]], vec![1.0])
        .expect("ions");
    assert_eq!(
        hydrogen_atom.kohn_sham_step().err(),
        Some(QuantumChemistryMdError::OpenShell(1))
    );
    hydrogen_atom.mixing = 1.5;
    assert!(hydrogen_atom.forces().is_err());
}

#[test]
fn test_dft_xc_potential_is_functional_derivative_on_grid() {
    for functional in [XcFunctional::Lda, XcFunctional::Pbe, XcFunctional::B3lyp] {
        let grid = UniformGrid::centered([0.0; 3], 0.3, [12, 12, 12]).expect("grid");
        let density = gaussian_density(&grid, [0.1, -0.05, 0.0], 0.8, 2.0);
        let dft = DensityFunctionalMd::new(functional.clone(), grid.clone(), density.clone())
            .expect("dft");
        let base = dft.xc_on_grid(&density);
        let dv = grid.dv();
        for &idx in &[
            grid.index(6, 6, 6),
            grid.index(4, 7, 5),
            grid.index(8, 3, 6),
        ] {
            let delta = 1e-5 * density[idx];
            let mut plus = density.clone();
            plus[idx] += delta;
            let mut minus = density.clone();
            minus[idx] -= delta;
            let e_plus = dft.xc_on_grid(&plus);
            let e_minus = dft.xc_on_grid(&minus);
            let fd = ((e_plus.exchange + e_plus.correlation)
                - (e_minus.exchange + e_minus.correlation))
                / (2.0 * delta * dv);
            let v = base.potential[idx];
            assert!(
                (fd - v).abs() < 1e-6 * v.abs().max(1.0),
                "{functional:?} idx {idx}: {fd} vs {v}"
            );
        }
    }
}

#[test]
fn test_dft_helium_lda_scf() {
    let mut he = helium(XcFunctional::Lda);
    let report = he.run_scf().expect("scf converges");
    let e = &report.energies;
    assert!((he.integrated_density() - 2.0).abs() < 1e-9);
    assert!(e.kinetic > 0.0 && e.external < 0.0 && e.hartree > 0.0);
    assert!(e.exchange < 0.0 && e.correlation < 0.0);
    assert_eq!(e.exact_exchange, 0.0);
    assert_eq!(e.ion_ion, 0.0);
    let sum = e.kinetic + e.external + e.hartree + e.exchange + e.correlation;
    assert!((sum - e.total).abs() < 1e-12);
    // The 1s level is bound and the atom is bound relative to He²⁺ + 2e⁻.
    assert!(report.eigenvalues[0] < -0.25, "{:?}", report.eigenvalues);
    // A Gaussian nucleus (σ = 1.5 h = 0.375 bohr) binds less than the point
    // nucleus of the all-electron LDA atom (−2.834 Ha).
    assert!(e.total < -1.9 && e.total > -2.834, "E = {}", e.total);
    // Spherical symmetry on a symmetric grid: no net force.
    let f = he.forces().expect("forces");
    for c in f[0] {
        assert!(c.abs() < 1e-8, "{c}");
    }
}

#[test]
fn test_dft_helium_b3lyp_exact_exchange_is_minus_half_hartree() {
    let mut he = helium(XcFunctional::B3lyp);
    let report = he.run_scf().expect("scf converges");
    // One doubly occupied orbital: E_x^HF = −(φφ|φφ) = −E_H/2.
    let e_x_hf = he.exact_exchange_energy().expect("orbitals");
    let e_h = he.hartree_energy();
    assert!((e_x_hf + 0.5 * e_h).abs() < 1e-9 * e_h, "{e_x_hf} vs {e_h}");
    assert!((report.energies.exact_exchange - 0.2 * e_x_hf).abs() < 1e-6);
    assert!(report.energies.exchange < 0.0 && report.energies.correlation < 0.0);
}

#[test]
fn test_dft_helium_pbe_scf_converges() {
    let mut lda = helium(XcFunctional::Lda);
    let mut pbe = helium(XcFunctional::Pbe);
    let e_lda = lda.run_scf().expect("lda").energies;
    let e_pbe = pbe.run_scf().expect("pbe").energies;
    // The GGA exchange enhancement makes E_x more negative than LDA exchange
    // and the PBE correlation is much smaller in magnitude than LDA's.
    assert!(e_pbe.exchange < e_lda.exchange);
    assert!(e_pbe.correlation > e_lda.correlation);
    assert!((pbe.integrated_density() - 2.0).abs() < 1e-9);
}

#[test]
fn test_dft_hellmann_feynman_force_matches_energy_derivative() {
    let mut h2 = hydrogen_molecule(XcFunctional::Lda, 0.8);
    h2.scf_tolerance = 1e-9;
    h2.eigensolver_tolerance = 1e-9;
    h2.mixing = 0.5;
    h2.run_scf().expect("scf");
    let f = h2.forces().expect("forces");
    let delta = 0.02;
    let mut energies = [0.0; 2];
    for (slot, sign) in energies.iter_mut().zip([1.0, -1.0]) {
        let mut displaced = h2.clone();
        displaced.ionic_positions[1][0] += sign * delta;
        *slot = displaced.run_scf().expect("scf").energies.total;
    }
    let fd_force = -(energies[0] - energies[1]) / (2.0 * delta);
    assert!(
        (f[1][0] - fd_force).abs() < 2e-4,
        "HF force {} vs −dE/dR {fd_force}",
        f[1][0]
    );
}

#[test]
fn test_dft_h2_bond_restoring_forces() {
    let mut stretched = hydrogen_molecule(XcFunctional::Lda, 1.2);
    stretched.run_scf().expect("scf");
    let f = stretched.forces().expect("forces");
    assert!(f[1][0] < 0.0 && f[0][0] > 0.0, "{f:?}");
    let mut compressed = hydrogen_molecule(XcFunctional::Lda, 0.5);
    compressed.run_scf().expect("scf");
    let f = compressed.forces().expect("forces");
    assert!(f[1][0] > 0.0 && f[0][0] < 0.0, "{f:?}");
}

#[test]
fn test_dft_bomd_step_conserves_energy() {
    let mut h2 = hydrogen_molecule(XcFunctional::Lda, 0.9);
    h2.scf_tolerance = 1e-9;
    h2.eigensolver_tolerance = 1e-9;
    h2.mixing = 0.5;
    let e0 = h2.run_scf().expect("scf").energies.total + h2.ionic_kinetic_energy();
    let x0 = h2.ionic_positions[1][0];
    let mut e1 = e0;
    for _ in 0..3 {
        let report = h2.step(10.0).expect("step");
        e1 = report.energies.total + h2.ionic_kinetic_energy();
    }
    assert!((h2.ionic_positions[1][0] - x0).abs() > 1e-4);
    assert!((e1 - e0).abs() < 2e-4, "drift {}", e1 - e0);
    assert!(h2.step(-1.0).is_err());
}

#[test]
fn test_dft_helium_energy_decreases_towards_point_nucleus_limit() {
    let mut wide = helium(XcFunctional::Lda);
    let e_wide = wide.run_scf().expect("scf").energies.total;
    let mut narrow = helium(XcFunctional::Lda);
    narrow.nuclear_smearing = narrow.grid.spacing;
    let e_narrow = narrow.run_scf().expect("scf").energies.total;
    assert!(e_narrow < e_wide - 0.1, "{e_narrow} vs {e_wide}");
    assert!(e_narrow > -2.834, "{e_narrow}");
}

#[test]
fn test_dft_pulay_and_linear_mixing_reach_the_same_ground_state() {
    let mut pulay = helium(XcFunctional::Lda);
    pulay.scf_tolerance = 1e-8;
    let r_pulay = pulay.run_scf().expect("pulay scf");
    let mut linear = helium(XcFunctional::Lda);
    linear.scf_tolerance = 1e-8;
    linear.pulay_history = 0;
    let r_linear = linear.run_scf().expect("linear scf");
    assert!((r_pulay.energies.total - r_linear.energies.total).abs() < 1e-7);
    assert!(
        r_pulay.steps < r_linear.steps,
        "{} vs {}",
        r_pulay.steps,
        r_linear.steps
    );
}
