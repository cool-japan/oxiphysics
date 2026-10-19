// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the crack-tip stress-intensity estimate of
//! [`FractureMesh::propagate_crack`].
//!
//! The estimate has to come from the elastic energy stored in the tensile
//! bonds at the tip.  The net force on the tip node is the wrong quantity: in
//! static equilibrium it vanishes however heavily the crack is loaded (so a
//! loaded crack never propagated), and a force divided by `sqrt(pi r)` has the
//! unit N/sqrt(m), not Pa sqrt(m).

use super::{CrackTip, FractureMaterial, FractureMesh, notched_tension_specimen};

const YOUNG: f64 = 70.0e9;
const K_IC: f64 = 0.75e6;
const DX: f64 = 0.1;

fn glass() -> FractureMaterial {
    FractureMaterial::new(YOUNG, 0.23, 2500.0, 50.0e6, K_IC)
}

/// A 5 x 5 grid without gravity and a crack tip on its interior centre node.
fn grid_with_interior_tip() -> (FractureMesh, usize) {
    let mut mesh = FractureMesh::new_grid(5, 5, DX, glass());
    mesh.gravity = [0.0; 3];
    let tip_node = 2 * 5 + 2;
    mesh.crack_tips
        .push(CrackTip::new(tip_node, [1.0, 0.0, 0.0]));
    (mesh, tip_node)
}

/// Uniform affine strain `strain` along y (a homogeneous tension state).
fn stretch_y(mesh: &mut FractureMesh, strain: f64) {
    for node in &mut mesh.nodes {
        node.pos[1] *= 1.0 + strain;
    }
}

/// Closed-form estimate for the interior node of a y-stretched `new_grid`.
///
/// The node carries 2 horizontal bonds (no extension), 2 vertical bonds and 4
/// diagonal bonds of stiffness `E dx`.  Only tensile bonds store energy that
/// counts; the new crack surface is the mean bond length times a unit
/// thickness.
fn expected_k_i(strain: f64) -> f64 {
    let stiffness = YOUNG * DX;
    let ext_vertical = strain * DX;
    let ext_diagonal = DX * ((1.0 + (1.0 + strain).powi(2)).sqrt() - 2.0_f64.sqrt());
    let tensile = |ext: f64| {
        if ext > 0.0 {
            0.5 * stiffness * ext * ext
        } else {
            0.0
        }
    };
    let energy = 2.0 * tensile(ext_vertical) + 4.0 * tensile(ext_diagonal);
    let mean_bond_length = DX * (4.0 + 4.0 * 2.0_f64.sqrt()) / 8.0;
    (YOUNG * energy / mean_bond_length).sqrt()
}

fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|c| c * c).sum::<f64>().sqrt()
}

#[test]
fn stretched_grid_propagates_although_the_net_force_on_the_tip_vanishes() {
    let strain = 2.0e-4;
    let (mut mesh, tip_node) = grid_with_interior_tip();
    stretch_y(&mut mesh, strain);
    assert!(
        strain < mesh.material.critical_stretch(),
        "the test must stay below the bond-break threshold"
    );

    // The interior node of a homogeneously strained lattice is in equilibrium,
    // exactly the situation in which the old force-based estimate was zero.
    let forces = mesh.compute_bond_forces();
    let bond_force_scale = YOUNG * DX * strain * DX;
    assert!(
        bond_force_scale > 1.0e3,
        "the bonds themselves carry a large force: {bond_force_scale}"
    );
    assert!(
        norm(forces[tip_node]) < 1.0e-9 * bond_force_scale,
        "net force on the tip node {} must vanish",
        norm(forces[tip_node])
    );

    mesh.propagate_crack(1.0e-6);

    let k_i = mesh.crack_tips[0].stress_intensity;
    let expected = expected_k_i(strain);
    assert!(
        (k_i - expected).abs() < 1.0e-9 * expected,
        "K_I = {k_i} but the stored tensile bond energy gives {expected}"
    );
    assert!(k_i > K_IC, "K_I = {k_i} must exceed K_Ic = {K_IC}");
    // The bond ahead of the tip (towards +x) broke and the tip moved onto it.
    assert_eq!(mesh.broken_bond_count(), 1);
    assert_eq!(mesh.crack_tips[0].node_idx, tip_node + 1);
    assert!((mesh.crack_length() - DX).abs() < 1.0e-12);
}

/// Rounding of the node positions makes a bond's current length differ from its
/// rest length in the last ulp, so an unloaded lattice shows a stored energy of
/// ~1e-34 J and a K_I of ~1e-6 Pa sqrt(m): this is numerical noise, twelve
/// orders of magnitude below K_Ic, not a loading of the tip.
fn assert_unloaded(k_i: f64) {
    assert!(
        (0.0..1.0e-9 * K_IC).contains(&k_i),
        "K_I = {k_i} must vanish up to rounding (K_Ic = {K_IC})"
    );
}

#[test]
fn unloaded_grid_does_not_propagate() {
    let (mut mesh, tip_node) = grid_with_interior_tip();
    mesh.propagate_crack(1.0e-6);
    assert_unloaded(mesh.crack_tips[0].stress_intensity);
    assert_eq!(mesh.broken_bond_count(), 0);
    assert_eq!(mesh.crack_tips[0].node_idx, tip_node);
    assert_eq!(mesh.crack_length(), 0.0);
}

#[test]
fn subcritical_stretch_does_not_propagate() {
    let strain = 1.0e-5;
    let (mut mesh, tip_node) = grid_with_interior_tip();
    stretch_y(&mut mesh, strain);
    mesh.propagate_crack(1.0e-6);

    let k_i = mesh.crack_tips[0].stress_intensity;
    let expected = expected_k_i(strain);
    assert!(
        expected < K_IC,
        "test set-up: K = {expected} is subcritical"
    );
    assert!((k_i - expected).abs() < 1.0e-9 * expected);
    assert!(k_i > 0.0, "a stretched tip has a positive K_I");
    assert_eq!(mesh.broken_bond_count(), 0);
    assert_eq!(mesh.crack_tips[0].node_idx, tip_node);
}

#[test]
fn compression_does_not_drive_a_crack() {
    let (mut mesh, tip_node) = grid_with_interior_tip();
    stretch_y(&mut mesh, -2.0e-4);
    mesh.propagate_crack(1.0e-6);
    assert_unloaded(mesh.crack_tips[0].stress_intensity);
    assert_eq!(mesh.broken_bond_count(), 0);
    assert_eq!(mesh.crack_tips[0].node_idx, tip_node);
}

#[test]
fn stress_intensity_is_linear_in_strain_with_pa_sqrt_m_scale() {
    let k_of = |strain: f64| {
        let (mut mesh, _) = grid_with_interior_tip();
        stretch_y(&mut mesh, strain);
        mesh.propagate_crack(1.0e-6);
        mesh.crack_tips[0].stress_intensity
    };
    let k1 = k_of(1.0e-5);
    let k2 = k_of(2.0e-5);
    assert!(
        (k2 / k1 - 2.0).abs() < 1.0e-3,
        "K_I = sqrt(E G) with G ~ strain^2 must be linear in the strain: {k1} -> {k2}"
    );
    // sigma sqrt(pi r) with sigma = E strain and r = dx is the expected order
    // of magnitude of a Pa sqrt(m) quantity at this strain.
    let sigma_sqrt_pi_r = YOUNG * 1.0e-5 * (std::f64::consts::PI * DX).sqrt();
    assert!(
        k1 > 0.1 * sigma_sqrt_pi_r && k1 < 10.0 * sigma_sqrt_pi_r,
        "K_I = {k1} is not of the order of sigma sqrt(pi r) = {sigma_sqrt_pi_r}"
    );
}

#[test]
fn notched_specimen_under_tension_advances_its_crack() {
    let mut mesh = notched_tension_specimen(1.0, 0.5, 0.2, DX, glass());
    mesh.gravity = [0.0; 3];
    let intact_before = mesh.intact_bond_count();
    let tip_before = mesh.crack_tips[0].node_idx;
    stretch_y(&mut mesh, 2.0e-4);

    mesh.propagate_crack(1.0e-6);

    assert!(
        mesh.crack_tips[0].stress_intensity > K_IC,
        "K_I = {} must exceed K_Ic",
        mesh.crack_tips[0].stress_intensity
    );
    assert_eq!(mesh.intact_bond_count(), intact_before - 1);
    assert_ne!(mesh.crack_tips[0].node_idx, tip_before);
    assert!(mesh.crack_tips[0].total_length > 0.0);
}
