// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `dipole_gradient_force` must be the exact negative gradient of the
//! charge/dipole pair energy
//!
//! ```text
//! U_ij = q_i q_j / r + q_i (mu_j.r) / r^3 - q_j (mu_i.r) / r^3 + U_dd(mu_i, mu_j, r),
//! r = r_i - r_j,
//! ```
//!
//! with the dipoles held fixed.  The reference energy is written out here and the
//! force is checked against central finite differences of it.

use oxiphysics_md::polarizable_md::{
    PolarizableAtom, dipole_dipole_interaction, dipole_gradient_force,
};

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn atom(position: [f64; 3], charge: f64, dipole: [f64; 3]) -> PolarizableAtom {
    let mut a = PolarizableAtom::new(position, charge, 1.0, 1.0);
    a.dipole = dipole;
    a
}

/// Pair energy of the charge-dipole model, `r = r_a - r_b`.
fn pair_energy(a: &PolarizableAtom, b: &PolarizableAtom) -> f64 {
    let r_vec = sub(a.position, b.position);
    let r = dot(r_vec, r_vec).sqrt();
    let r3 = r * r * r;
    let charge_charge = a.charge * b.charge / r;
    let charge_dipole = (a.charge * dot(b.dipole, r_vec) - b.charge * dot(a.dipole, r_vec)) / r3;
    charge_charge + charge_dipole + dipole_dipole_interaction(a.dipole, b.dipole, r_vec)
}

fn total_energy(atoms: &[PolarizableAtom]) -> f64 {
    let mut energy = 0.0;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            energy += pair_energy(&atoms[i], &atoms[j]);
        }
    }
    energy
}

/// `-dU/dx_{i,d}` by central finite differences with step `h`.
fn numerical_force(atoms: &[PolarizableAtom], i: usize, d: usize, h: f64) -> f64 {
    let mut plus = atoms.to_vec();
    plus[i].position[d] += h;
    let mut minus = atoms.to_vec();
    minus[i].position[d] -= h;
    -(total_energy(&plus) - total_energy(&minus)) / (2.0 * h)
}

fn generic_cluster() -> Vec<PolarizableAtom> {
    vec![
        atom([0.0, 0.0, 0.0], 0.8, [0.31, -0.22, 0.17]),
        atom([2.3, 0.4, -0.3], -0.5, [-0.12, 0.40, 0.26]),
        atom([0.7, 2.6, 0.9], 0.3, [0.28, 0.09, -0.35]),
        atom([-1.4, 1.1, 2.4], -0.6, [-0.20, -0.30, 0.15]),
    ]
}

#[test]
fn force_is_the_negative_gradient_of_the_pair_energy() {
    let atoms = generic_cluster();
    let mut worst = 0.0_f64;
    let mut scale = 0.0_f64;
    for i in 0..atoms.len() {
        let force = dipole_gradient_force(&atoms, i);
        for (d, &f) in force.iter().enumerate() {
            let fd = numerical_force(&atoms, i, d, 1.0e-5);
            worst = worst.max((f - fd).abs());
            scale = scale.max(fd.abs());
        }
    }
    assert!(
        scale > 1.0e-2,
        "the test forces must be non-trivial: {scale}"
    );
    assert!(
        worst < 1.0e-6 * scale,
        "max |F - (-dU/dx)| = {worst:e} for force scale {scale:e}"
    );
}

#[test]
fn dipole_dipole_part_alone_matches_the_gradient_of_its_energy() {
    // Neutral atoms: only the (mu_i . grad) E_j term is left.
    let mut atoms = generic_cluster();
    for a in &mut atoms {
        a.charge = 0.0;
    }
    let mut scale = 0.0_f64;
    let mut worst = 0.0_f64;
    for i in 0..atoms.len() {
        let force = dipole_gradient_force(&atoms, i);
        for (d, &f) in force.iter().enumerate() {
            let fd = numerical_force(&atoms, i, d, 1.0e-5);
            worst = worst.max((f - fd).abs());
            scale = scale.max(fd.abs());
        }
    }
    assert!(scale > 1.0e-3, "{scale}");
    assert!(worst < 1.0e-6 * scale, "{worst:e} vs {scale:e}");
}

#[test]
fn charge_dipole_part_alone_matches_the_gradient_of_its_energy() {
    // Charged atoms, but only atom 0 and atom 2 carry dipoles.
    let mut atoms = generic_cluster();
    atoms[1].dipole = [0.0; 3];
    atoms[3].dipole = [0.0; 3];
    let mut scale = 0.0_f64;
    let mut worst = 0.0_f64;
    for i in 0..atoms.len() {
        let force = dipole_gradient_force(&atoms, i);
        for (d, &f) in force.iter().enumerate() {
            let fd = numerical_force(&atoms, i, d, 1.0e-5);
            worst = worst.max((f - fd).abs());
            scale = scale.max(fd.abs());
        }
    }
    assert!(worst < 1.0e-6 * scale, "{worst:e} vs {scale:e}");
}

#[test]
fn forces_obey_newtons_third_law() {
    let atoms = generic_cluster();
    let mut total = [0.0_f64; 3];
    let mut scale = 0.0_f64;
    for i in 0..atoms.len() {
        let force = dipole_gradient_force(&atoms, i);
        for d in 0..3 {
            total[d] += force[d];
            scale = scale.max(force[d].abs());
        }
    }
    for t in total {
        assert!(t.abs() < 1.0e-12 * scale.max(1.0), "net force {t:e}");
    }
}

#[test]
fn collinear_head_to_tail_dipoles_attract_with_six_mu_squared_over_d_to_the_fourth() {
    let (mu, d) = (0.4, 3.0);
    let atoms = vec![
        atom([d, 0.0, 0.0], 0.0, [mu, 0.0, 0.0]),
        atom([0.0, 0.0, 0.0], 0.0, [mu, 0.0, 0.0]),
    ];
    // U = -2 mu^2 / d^3 and r = r_0 - r_1 = +d x: F_0 = -6 mu^2 / d^4 x (towards atom 1).
    let expected = -6.0 * mu * mu / d.powi(4);
    let f0 = dipole_gradient_force(&atoms, 0);
    let f1 = dipole_gradient_force(&atoms, 1);
    assert!((f0[0] - expected).abs() < 1.0e-14, "{f0:?}");
    assert!((f1[0] + expected).abs() < 1.0e-14, "{f1:?}");
    assert!(f0[1].abs() < 1.0e-14 && f0[2].abs() < 1.0e-14);
}

#[test]
fn point_charges_reduce_to_the_coulomb_force() {
    let atoms = vec![
        atom([0.0, 0.0, 0.0], 1.5, [0.0; 3]),
        atom([0.0, 2.0, 0.0], 2.0, [0.0; 3]),
    ];
    // Like charges repel: atom 1 is pushed to +y with q1 q2 / r^2.
    let f1 = dipole_gradient_force(&atoms, 1);
    assert!((f1[1] - 1.5 * 2.0 / 4.0).abs() < 1.0e-14, "{f1:?}");
    let f0 = dipole_gradient_force(&atoms, 0);
    assert!((f0[1] + 1.5 * 2.0 / 4.0).abs() < 1.0e-14, "{f0:?}");
}

#[test]
fn dipole_in_the_field_of_a_point_charge_feels_the_field_gradient() {
    // A dipole mu on atom 0 near the charge q on atom 1 (r = r_0 - r_1 = d x):
    // U = -q (mu.r) / r^3, F_0 = q [mu / r^3 - 3 (mu.r) r / r^5].
    let (q, d) = (1.2, 2.5);
    let mu = [0.3, 0.2, -0.1];
    let atoms = vec![atom([d, 0.0, 0.0], 0.0, mu), atom([0.0; 3], q, [0.0; 3])];
    let f = dipole_gradient_force(&atoms, 0);
    let expected = [
        q * (mu[0] / d.powi(3) - 3.0 * mu[0] * d * d / d.powi(5)),
        q * mu[1] / d.powi(3),
        q * mu[2] / d.powi(3),
    ];
    for axis in 0..3 {
        assert!(
            (f[axis] - expected[axis]).abs() < 1.0e-14,
            "axis {axis}: {f:?} vs {expected:?}"
        );
    }
}

#[test]
fn out_of_range_index_and_coincident_atoms_give_zero_force() {
    let atoms = generic_cluster();
    assert_eq!(dipole_gradient_force(&atoms, atoms.len()), [0.0; 3]);
    let overlapping = vec![
        atom([1.0, 1.0, 1.0], 1.0, [0.1, 0.2, 0.3]),
        atom([1.0, 1.0, 1.0], -1.0, [0.3, 0.2, 0.1]),
    ];
    assert_eq!(dipole_gradient_force(&overlapping, 0), [0.0; 3]);
    assert_eq!(dipole_gradient_force(&[], 0), [0.0; 3]);
}
