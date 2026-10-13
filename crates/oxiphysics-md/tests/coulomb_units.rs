// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit-consistency regression tests for the Å-based Coulomb constant.
//!
//! Every Å-based electrostatics routine of the crate (direct Coulomb, Ewald,
//! PME, Born/GB solvation) must share the single constant
//! `e² N_A / (4π ε₀) = 1389.354576 kJ·Å/(mol·e²)`.  The nm-valued GROMACS
//! constant 138.935 is ten times too small for Å distances; the Madelung
//! ratio tests divide by the constant and therefore cannot see such an
//! error, so these tests pin the *absolute* energy scale.

use oxiphysics_md::electrostatics::{
    BornSolvation, COULOMB_K, coulomb_energy, coulomb_force_scalar, pme,
};
use oxiphysics_md::ewald::{self, EwaldParams, EwaldSummation};
use oxiphysics_md::ionic_liquid_md::WolfParams;

/// Madelung constant of the rock-salt structure.
const NACL_MADELUNG: f64 = 1.747_564_594_6;

/// `e² N_A / (4π ε₀)` from CODATA 2018 SI values, converted to kJ·Å/mol.
fn coulomb_constant_from_si() -> f64 {
    let elementary_charge = 1.602_176_634e-19_f64; // C
    let avogadro = 6.022_140_76e23_f64; // 1/mol
    let vacuum_permittivity = 8.854_187_812_8e-12_f64; // F/m
    // J·m/mol -> kJ·Å/mol is a factor 1e-3 · 1e10 = 1e7.
    elementary_charge * elementary_charge * avogadro
        / (4.0 * std::f64::consts::PI * vacuum_permittivity)
        * 1.0e7
}

/// Rock-salt checkerboard of `n_per_side`³ ions with nearest-neighbour
/// distance `d_nn`, returned as `(positions, charges, box_edge)`.
fn rock_salt(n_per_side: usize, d_nn: f64) -> (Vec<[f64; 3]>, Vec<f64>, f64) {
    let mut positions = Vec::with_capacity(n_per_side.pow(3));
    let mut charges = Vec::with_capacity(n_per_side.pow(3));
    for ix in 0..n_per_side {
        for iy in 0..n_per_side {
            for iz in 0..n_per_side {
                positions.push([ix as f64 * d_nn, iy as f64 * d_nn, iz as f64 * d_nn]);
                charges.push(if (ix + iy + iz) % 2 == 0 { 1.0 } else { -1.0 });
            }
        }
    }
    (positions, charges, n_per_side as f64 * d_nn)
}

#[test]
fn coulomb_constant_matches_si_derivation() {
    let reference = coulomb_constant_from_si();
    assert!(
        ((COULOMB_K - reference) / reference).abs() < 1e-7,
        "COULOMB_K = {COULOMB_K}, SI derivation = {reference}"
    );
}

#[test]
fn ewald_pme_and_direct_coulomb_share_one_constant() {
    assert_eq!(ewald::COULOMB_K, COULOMB_K);
    assert_eq!(pme::COULOMB_K, COULOMB_K);
    // The direct Coulomb energy and force are expressed through the same value.
    assert!((coulomb_energy(1.0, 1.0, 1.0) - COULOMB_K).abs() < 1e-9);
    assert!((coulomb_force_scalar(1.0, 1.0, 1.0) - COULOMB_K).abs() < 1e-9);
}

#[test]
fn ewald_nacl_lattice_energy_per_ion_pair_is_861_kj_per_mol() {
    // NaCl: d_nn = 2.82 Å.  E_pair = -M K / d_nn = -860.99 kJ/mol; the old
    // nm-valued constant would give -86.1 kJ/mol.
    let d_nn = 2.82_f64;
    let (positions, charges, box_edge) = rock_salt(4, d_nn);
    let n_pairs = (positions.len() / 2) as f64;

    // Minimum-image cutoff just below L/2 with erfc(alpha r_cut) ~ 1e-12.
    let r_cut = 0.5 * box_edge - 0.04;
    let alpha = 5.0 / r_cut;
    let params = EwaldParams {
        alpha,
        r_cutoff: r_cut,
        k_cutoff: 2.0 * std::f64::consts::PI * alpha,
        epsilon_r: 1.0,
    };
    let ewald = EwaldSummation::with_k_max(params, [box_edge; 3], 12);
    let e_per_pair = ewald.total_energy(&positions, &charges) / n_pairs;

    let expected = -NACL_MADELUNG * COULOMB_K / d_nn;
    assert!(
        (expected + 860.988).abs() < 0.01,
        "reference energy per ion pair {expected}"
    );
    assert!(
        ((e_per_pair - expected) / expected).abs() < 5e-3,
        "Ewald energy per ion pair {e_per_pair:.3} kJ/mol, expected {expected:.3} kJ/mol"
    );
}

#[test]
fn born_solvation_energy_of_a_monovalent_cation_is_408_kj_per_mol() {
    // -(K / 2a)(1 - 1/eps) for q = +1, a = 1.68 Å in water (eps = 78.5).
    let born = BornSolvation::new(78.5);
    let delta_g = born.solvation_energy(1.0, 1.68);
    assert!(
        (delta_g + 408.230_890).abs() < 1e-3,
        "Born solvation energy = {delta_g} kJ/mol"
    );
}

#[test]
fn ionic_liquid_wolf_summation_uses_the_angstrom_constant() {
    // alpha = 0.2 / Å, r_cut = 12 Å: a +1/-1 pair at 4 Å has
    // K [erfc(0.8)/4 - erfc(2.4)/12] with the sign of -1, i.e. -89.4986 kJ/mol; the
    // nm-valued constant 138.935 would give only -8.95 kJ/mol.
    let wolf = WolfParams::default_il();
    assert_eq!(wolf.coulomb_conv, COULOMB_K);
    let energy = wolf.energy(1.0, -1.0, 4.0);
    // The module's erfc approximation is accurate to 1.5e-7 absolute.
    assert!(
        (energy + 89.498_585_379_604_1).abs() < 1e-3,
        "Wolf pair energy = {energy} kJ/mol"
    );
}
