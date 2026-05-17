// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase 21.6 validation: NaCl Madelung constant via Ewald summation.
//!
//! Computes the total electrostatic energy of a 4x4x4 simple-cubic checkerboard
//! ion lattice (64 ions, alternating +/- charges, unit spacing) under 3D
//! periodic boundary conditions using the Ewald sum implemented in
//! `oxiphysics_md::ewald`.  The energy per ion pair, divided by the Coulomb
//! constant and the nearest-neighbour distance, yields -M where M is the
//! dimensionless Madelung constant of the NaCl rock-salt lattice.
//!
//! Reference value: M = 1.7475645946 (Madelung constant for NaCl structure).

#![allow(missing_docs)]

#[path = "regression_harness.rs"]
mod harness;

use oxiphysics::md::ewald::{EwaldParams, EwaldSummation};

/// Coulomb constant used by the Ewald implementation (kJ*angstrom*mol^-1*e^-2).
///
/// Mirrors `oxiphysics_md::ewald::params::COULOMB_K` (not re-exported as a pub
/// constant at the module root, so we redefine it locally with the identical
/// value; used only to convert the measured energy into dimensionless units).
const COULOMB_K_LOCAL: f64 = 138.935;

/// Build the 4x4x4 simple-cubic rock-salt ion lattice.
///
/// - Simple-cubic spacing `a0 = 1` (reduced units).
/// - 4 ions along each axis -> N = 64 ions.
/// - Periodic cubic box of edge `L = 4 * a0 = 4`.
/// - Charge `q = (-1)^(ix + iy + iz)`.
///
/// The checkerboard is image-consistent: the box edge is an even multiple of
/// the simple-cubic spacing, so translating by L preserves the checkerboard
/// sign pattern.  Nearest-neighbour distance is `a0 = 1`.
fn build_nacl_lattice(a0: f64, n_per_side: usize) -> (Vec<[f64; 3]>, Vec<f64>, f64) {
    let mut positions: Vec<[f64; 3]> = Vec::with_capacity(n_per_side.pow(3));
    let mut charges: Vec<f64> = Vec::with_capacity(n_per_side.pow(3));
    for ix in 0..n_per_side {
        for iy in 0..n_per_side {
            for iz in 0..n_per_side {
                let x = ix as f64 * a0;
                let y = iy as f64 * a0;
                let z = iz as f64 * a0;
                let parity = (ix + iy + iz) & 1;
                let q = if parity == 0 { 1.0 } else { -1.0 };
                positions.push([x, y, z]);
                charges.push(q);
            }
        }
    }
    let box_len = n_per_side as f64 * a0;
    (positions, charges, box_len)
}

#[test]
fn test_nacl_madelung_constant() {
    // Lattice parameters.
    let a0: f64 = 1.0;
    let n_per_side: usize = 4;
    let (positions, charges, box_len) = build_nacl_lattice(a0, n_per_side);
    let n_ions = positions.len();
    assert_eq!(n_ions, 64, "expected 64 ions in 4x4x4 lattice");

    // Charge-neutrality sanity check.
    let q_net: f64 = charges.iter().sum();
    assert!(
        q_net.abs() < 1e-12,
        "NaCl lattice must be charge-neutral, got net charge {q_net}"
    );

    // Ewald parameters.  Box edge L = 4, so the minimum-image cutoff is L/2 = 2.
    // Pick r_cut strictly below L/2 so all in-range pairs are captured.
    // alpha ~ 5/r_cut gives erfc(alpha*r_cut) ~ 1.5e-12 in real space.
    // k_max = 10 gives exp(-(pi*k/(alpha*L))^2) ~ 1e-4 in reciprocal space --
    // comfortably below the 5e-3 Madelung tolerance.
    let r_cut = 1.9_f64;
    let alpha = 5.0 / r_cut;
    let k_max = 10_i32;

    let params = EwaldParams {
        alpha,
        r_cutoff: r_cut,
        k_cutoff: 2.0 * std::f64::consts::PI * alpha,
        epsilon_r: 1.0,
    };
    let ewald = EwaldSummation::with_k_max(params, [box_len, box_len, box_len], k_max);

    // Total Ewald energy (real + reciprocal + self), in kJ/mol.
    let e_total = ewald.total_energy(&positions, &charges);
    assert!(e_total.is_finite(), "Ewald total energy must be finite");
    assert!(
        e_total < 0.0,
        "NaCl lattice energy must be negative (attractive), got {e_total}"
    );

    // Madelung constant via energy per ion pair:
    //   E_total = - M * COULOMB_K * (N/2) / a0
    //   => -M = (E_total / (N/2)) * a0 / COULOMB_K
    //
    // Baseline `ewald_madelung` stores the signed quantity -M = -1.7475645946.
    let n_pairs = (n_ions as f64) / 2.0;
    let e_per_pair = e_total / n_pairs;
    let measured_neg_m = e_per_pair * a0 / COULOMB_K_LOCAL;

    let baseline = harness::load_baseline("ewald_madelung").expect("load ewald_madelung baseline");
    eprintln!(
        "[ewald_madelung] measured = {measured_neg_m:.10}, expected = {:.10}, |diff| = {:.3e}",
        baseline.expected,
        (measured_neg_m - baseline.expected).abs(),
    );
    assert_close!(measured_neg_m, baseline);
}
