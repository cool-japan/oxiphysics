// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for ply-level recovery and first-ply failure.
//!
//! * `ply_stresses` rotated the laminate strain with the *stress*
//!   transformation, which doubled the shear coupling for every off-axis ply
//!   (`ε₁ += 2cs·γ_xy` instead of `cs·γ_xy`); the old tests used 0° plies only.
//! * `first_ply_failure` invented the ply strengths from the elastic moduli
//!   (`Xt = 0.01·E1`, ...), looked at the mid-plane only and returned the first
//!   stacking index with `FI ≥ 1` instead of the ply with the smallest strength
//!   ratio.

use oxiphysics_fem::Error;
use oxiphysics_fem::composite_fem::{
    Laminate, Ply, PlyStrength, first_ply_failure, laminate_abd_matrix, laminate_midplane_strains,
    ply_stresses, reduced_stiffness, rotation_matrix_2d, strain_transform_2d,
    transformed_stiffness, tsai_wu_failure, tsai_wu_strength_ratio,
};

const E1: f64 = 140.0e9;
const E2: f64 = 10.0e9;
const G12: f64 = 5.0e9;
const NU12: f64 = 0.3;
const THICKNESS: f64 = 2.0e-4;

fn ply(angle_deg: f64) -> Ply {
    Ply::new(THICKNESS, angle_deg, E1, E2, G12, NU12)
}

fn strength() -> PlyStrength {
    PlyStrength {
        xt: 1500.0e6,
        xc: 1200.0e6,
        yt: 50.0e6,
        yc: 200.0e6,
        s12: 70.0e6,
    }
}

fn laminate(angles: &[f64]) -> Laminate {
    let mut lam = Laminate::new();
    for &angle in angles {
        lam.add_ply(ply(angle));
    }
    lam
}

fn mat_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[test]
fn strain_transform_is_the_inverse_transpose_of_the_stress_transform() {
    for angle_deg in [-75.0_f64, -45.0, 0.0, 15.0, 30.0, 45.0, 60.0, 90.0, 137.0] {
        let theta = angle_deg.to_radians();
        let strain = strain_transform_2d(theta);
        let stress = rotation_matrix_2d(theta);
        // T_ε · Tᵀ = I ⇔ T_ε = (T⁻¹)ᵀ.
        for (i, strain_row) in strain.iter().enumerate() {
            for (j, stress_row) in stress.iter().enumerate() {
                let product: f64 = strain_row.iter().zip(stress_row).map(|(a, b)| a * b).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product - expected).abs() < 1e-14,
                    "{angle_deg}°: (T_ε·Tᵀ)[{i}][{j}] = {product}"
                );
            }
        }
    }
}

/// A 45° ply loaded by `γ_xy = 10⁻³` alone: `ε₁ = cs·γ = 5·10⁻⁴`,
/// `ε₂ = −cs·γ`, `γ₁₂ = (c² − s²)·γ = 0`.  The stress transformation applied to
/// the strain gave `ε₁ = 2cs·γ = 10⁻³`, twice the correct value.
#[test]
fn shear_strain_on_a_45_degree_ply_splits_into_equal_and_opposite_normal_strains() {
    let eps_mat = mat_vec(
        &strain_transform_2d(45.0_f64.to_radians()),
        [0.0, 0.0, 1.0e-3],
    );
    assert!((eps_mat[0] - 5.0e-4).abs() < 1e-16, "{eps_mat:?}");
    assert!((eps_mat[1] + 5.0e-4).abs() < 1e-16, "{eps_mat:?}");
    assert!(eps_mat[2].abs() < 1e-16, "{eps_mat:?}");

    let p = ply(45.0);
    let q = reduced_stiffness(&p);
    let sigma = ply_stresses(&p, [0.0, 0.0, 1.0e-3], 0.0, [0.0; 3]);
    let expected = [
        q[0][0] * 5.0e-4 - q[0][1] * 5.0e-4,
        q[1][0] * 5.0e-4 - q[1][1] * 5.0e-4,
        0.0,
    ];
    for i in 0..3 {
        assert!(
            (sigma[i] - expected[i]).abs() <= 1e-9 * expected[0].abs(),
            "σ[{i}] = {}, expected {}",
            sigma[i],
            expected[i]
        );
    }
}

/// `T⁻¹ · σ_mat` must equal `Q̄ · ε_global` for every angle: the recovered
/// ply stress is consistent with the transformed stiffness that built the ABD
/// matrix.
#[test]
fn recovered_ply_stress_is_consistent_with_the_transformed_stiffness() {
    let eps = [8.0e-4, -3.0e-4, 6.0e-4];
    let kappa = [0.2, -0.1, 0.05];
    let z = 0.7 * THICKNESS;
    for angle_deg in [-60.0_f64, -45.0, -20.0, 0.0, 10.0, 30.0, 45.0, 75.0, 90.0] {
        let p = ply(angle_deg);
        let sigma_mat = ply_stresses(&p, eps, z, kappa);
        let t_inverse = rotation_matrix_2d(-angle_deg.to_radians());
        let sigma_global = mat_vec(&t_inverse, sigma_mat);
        let strain_global = [
            eps[0] + z * kappa[0],
            eps[1] + z * kappa[1],
            eps[2] + z * kappa[2],
        ];
        let expected = mat_vec(&transformed_stiffness(&p), strain_global);
        for i in 0..3 {
            assert!(
                (sigma_global[i] - expected[i]).abs() <= 1e-9 * expected[0].abs().max(1.0e5),
                "{angle_deg}°, component {i}: {} vs {}",
                sigma_global[i],
                expected[i]
            );
        }
    }
}

#[test]
fn plies_rotated_by_ninety_degrees_swap_the_material_axes() {
    // A 90° ply under ε_x sees it as a transverse strain.
    let sigma = ply_stresses(&ply(90.0), [1.0e-3, 0.0, 0.0], 0.0, [0.0; 3]);
    let q = reduced_stiffness(&ply(0.0));
    assert!(
        (sigma[0] - q[0][1] * 1.0e-3).abs() <= 1e-9 * q[1][1] * 1.0e-3,
        "{sigma:?}"
    );
    assert!(
        (sigma[1] - q[1][1] * 1.0e-3).abs() <= 1e-9 * q[1][1] * 1.0e-3,
        "{sigma:?}"
    );
    assert!(sigma[2].abs() <= 1e-9 * q[1][1] * 1.0e-3, "{sigma:?}");
}

// ---------------------------------------------------------------------------
// Strength ratio
// ---------------------------------------------------------------------------

#[test]
fn strength_ratio_of_simple_stress_states() {
    let s = strength();
    let ratio = |stress: [f64; 3]| {
        tsai_wu_strength_ratio(stress, &s)
            .expect("valid strengths")
            .expect("finite ratio")
    };
    assert!((ratio([s.xt / 2.0, 0.0, 0.0]) - 2.0).abs() < 1e-12);
    assert!((ratio([-s.xc / 4.0, 0.0, 0.0]) - 4.0).abs() < 1e-12);
    assert!((ratio([0.0, s.yt / 3.0, 0.0]) - 3.0).abs() < 1e-12);
    assert!((ratio([0.0, -s.yc / 5.0, 0.0]) - 5.0).abs() < 1e-12);
    assert!((ratio([0.0, 0.0, s.s12 / 2.0]) - 2.0).abs() < 1e-12);
    assert!((ratio([0.0, 0.0, -s.s12 / 2.0]) - 2.0).abs() < 1e-12);
    assert!(
        tsai_wu_strength_ratio([0.0; 3], &s)
            .expect("valid strengths")
            .is_none()
    );
}

#[test]
fn scaling_a_biaxial_stress_by_the_ratio_gives_a_unit_tsai_wu_index() {
    let s = strength();
    let f1 = 1.0 / s.xt - 1.0 / s.xc;
    let f2 = 1.0 / s.yt - 1.0 / s.yc;
    let f11 = 1.0 / (s.xt * s.xc);
    let f22 = 1.0 / (s.yt * s.yc);
    let f12 = -0.5 * (f11 * f22).sqrt();
    for stress in [
        [300.0e6, 12.0e6, 20.0e6],
        [-400.0e6, -50.0e6, 10.0e6],
        [200.0e6, -60.0e6, -25.0e6],
        [-100.0e6, 15.0e6, 30.0e6],
    ] {
        let r = tsai_wu_strength_ratio(stress, &s)
            .expect("valid strengths")
            .expect("finite ratio");
        let scaled = [r * stress[0], r * stress[1], r * stress[2]];
        let index = tsai_wu_failure(scaled, f1, f2, f12, f11, f22, s.s12);
        assert!(
            (index - 1.0).abs() < 1e-12,
            "stress {stress:?}: F(R·σ) = {index}"
        );
    }
}

#[test]
fn invalid_strengths_and_stresses_are_rejected() {
    for bad in [
        PlyStrength {
            xt: 0.0,
            ..strength()
        },
        PlyStrength {
            xc: -1.0,
            ..strength()
        },
        PlyStrength {
            yt: f64::NAN,
            ..strength()
        },
        PlyStrength {
            yc: f64::INFINITY,
            ..strength()
        },
        PlyStrength {
            s12: 0.0,
            ..strength()
        },
    ] {
        assert!(
            matches!(bad.validate(), Err(Error::InvalidInput(_))),
            "{bad:?}"
        );
        assert!(matches!(
            tsai_wu_strength_ratio([1.0e6, 0.0, 0.0], &bad),
            Err(Error::InvalidInput(_))
        ));
    }
    assert!(PlyStrength::new(1.0, 1.0, 1.0, 1.0, 1.0).is_ok());
    assert!(PlyStrength::new(1.0, 1.0, -1.0, 1.0, 1.0).is_err());
    assert!(matches!(
        tsai_wu_strength_ratio([f64::NAN, 0.0, 0.0], &strength()),
        Err(Error::InvalidInput(_))
    ));
}

// ---------------------------------------------------------------------------
// First-ply failure
// ---------------------------------------------------------------------------

/// A single 0° ply under `N_x` carries `σ₁ = N_x/t`, so the strength ratio is
/// `R = X·t/|N_x|`: `R·N_x = Xᵀ·t` in tension and `Xᶜ·t` in compression.
#[test]
fn single_zero_degree_ply_fails_at_the_longitudinal_strength() {
    let lam = laminate(&[0.0]);
    let s = [strength()];
    let nx = 1000.0;
    let (index, ratio) = first_ply_failure(&lam, &s, [nx, 0.0, 0.0, 0.0, 0.0, 0.0])
        .expect("valid input")
        .expect("loaded ply has a finite ratio");
    assert_eq!(index, 0);
    let expected = strength().xt * THICKNESS / nx;
    assert!(
        (ratio - expected).abs() <= 1e-9 * expected,
        "R = {ratio}, expected {expected}"
    );

    let (_, compression) = first_ply_failure(&lam, &s, [-nx, 0.0, 0.0, 0.0, 0.0, 0.0])
        .expect("valid input")
        .expect("finite ratio");
    let expected = strength().xc * THICKNESS / nx;
    assert!(
        (compression - expected).abs() <= 1e-9 * expected,
        "R = {compression}, expected {expected}"
    );
}

#[test]
fn single_ply_transverse_and_shear_resultants() {
    let lam = laminate(&[0.0]);
    let s = [strength()];
    let load = 500.0;
    let (_, transverse) = first_ply_failure(&lam, &s, [0.0, load, 0.0, 0.0, 0.0, 0.0])
        .expect("valid input")
        .expect("finite ratio");
    let expected = strength().yt * THICKNESS / load;
    assert!(
        (transverse - expected).abs() <= 1e-9 * expected,
        "{transverse} vs {expected}"
    );

    let (_, shear) = first_ply_failure(&lam, &s, [0.0, 0.0, -load, 0.0, 0.0, 0.0])
        .expect("valid input")
        .expect("finite ratio");
    let expected = strength().s12 * THICKNESS / load;
    assert!(
        (shear - expected).abs() <= 1e-9 * expected,
        "{shear} vs {expected}"
    );
}

/// Pure bending of a single 0° ply: `σ₁ = 6·M_x/t²` on the surfaces.  The
/// compressive side (`Xᶜ < Xᵀ`) governs, which a mid-plane-only check cannot see.
#[test]
fn surfaces_not_the_mid_plane_govern_bending_failure() {
    let lam = laminate(&[0.0]);
    let s = [strength()];
    let moment = 0.8;
    let (index, ratio) = first_ply_failure(&lam, &s, [0.0, 0.0, 0.0, moment, 0.0, 0.0])
        .expect("valid input")
        .expect("finite ratio");
    assert_eq!(index, 0);
    let expected = strength().xc * THICKNESS * THICKNESS / (6.0 * moment);
    assert!(
        (ratio - expected).abs() <= 1e-9 * expected,
        "{ratio} vs {expected}"
    );
}

/// The ply with the smallest strength ratio is reported, not the first one in
/// the stacking order that exceeds 1.  The transverse tensile strength is far
/// below the longitudinal one, so under `N_x` the 90° plies of a [0/90/90/0]
/// laminate fail first; scaling the loads by `R` brings exactly that ply (and
/// no other) to `FI = 1`.
#[test]
fn the_weakest_ply_is_found_and_reaches_a_unit_failure_index() {
    let lam = laminate(&[0.0, 90.0, 90.0, 0.0]);
    let strengths = vec![strength(); 4];
    let loads = [800.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let (index, ratio) = first_ply_failure(&lam, &strengths, loads)
        .expect("valid input")
        .expect("finite ratio");
    // Plies 1 and 2 are mirror images: either may be reported.
    assert!(
        index == 1 || index == 2,
        "weakest ply is a 90° ply, got {index}"
    );

    // Independent evaluation of every ply at the scaled load.
    let s = strength();
    let f1 = 1.0 / s.xt - 1.0 / s.xc;
    let f2 = 1.0 / s.yt - 1.0 / s.yc;
    let f11 = 1.0 / (s.xt * s.xc);
    let f22 = 1.0 / (s.yt * s.yc);
    let f12 = -0.5 * (f11 * f22).sqrt();
    let abd = laminate_abd_matrix(&lam);
    let sol = laminate_midplane_strains(&abd, loads.map(|l| l * ratio));
    let eps0 = [sol[0], sol[1], sol[2]];
    let kappa = [sol[3], sol[4], sol[5]];
    let mut z = -lam.total_thickness() / 2.0;
    let mut ply_index = [0.0_f64; 4];
    for (k, p) in lam.plies.iter().enumerate() {
        for surface in [z, z + p.thickness] {
            let stress = ply_stresses(p, eps0, surface, kappa);
            let fi = tsai_wu_failure(stress, f1, f2, f12, f11, f22, s.s12);
            ply_index[k] = ply_index[k].max(fi);
        }
        z += p.thickness;
    }
    assert!(
        (ply_index[index] - 1.0).abs() < 1e-9,
        "FI of the reported ply at R·loads = {}",
        ply_index[index]
    );
    for (k, &fi) in ply_index.iter().enumerate() {
        assert!(
            fi <= 1.0 + 1e-9,
            "ply {k} exceeds the criterion first: FI = {fi}"
        );
    }
    // The 0° plies are far from failure at that load.
    assert!(ply_index[0] < 0.5 && ply_index[3] < 0.5, "{ply_index:?}");

    // The response is linear: halving the loads doubles R.
    let (_, half) = first_ply_failure(&lam, &strengths, loads.map(|l| 0.5 * l))
        .expect("valid input")
        .expect("finite ratio");
    assert!((half - 2.0 * ratio).abs() <= 1e-9 * ratio);
}

/// A balanced [45/−45]_s laminate under pure `N_xy` has `A₁₆ = A₂₆ = B = 0`, so
/// `ε_x = ε_y = 0` and `γ_xy = N_xy/A₆₆`.  The material-frame strains of a
/// ±45° ply are then `ε₁ = ±γ/2`, `ε₂ = ∓γ/2`, `γ₁₂ = 0` (strain transformation
/// `T_ε`).  The stress transformation gave `±γ`, i.e. stresses twice as large
/// and a strength ratio half as large.
#[test]
fn off_axis_plies_use_the_corrected_strain_transformation_in_the_failure_check() {
    let lam = laminate(&[45.0, -45.0, -45.0, 45.0]);
    let strengths = vec![strength(); 4];
    let nxy = 600.0;
    let loads = [0.0, 0.0, nxy, 0.0, 0.0, 0.0];
    let (index, ratio) = first_ply_failure(&lam, &strengths, loads)
        .expect("valid input")
        .expect("finite ratio");

    let abd = laminate_abd_matrix(&lam);
    assert!(abd[0][2].abs() < 1e-6 * abd[2][2] && abd[1][2].abs() < 1e-6 * abd[2][2]);
    let gamma = nxy / abd[2][2];
    let q = reduced_stiffness(&ply(0.0));
    let sigma_plus = [
        (q[0][0] - q[0][1]) * gamma / 2.0,
        (q[0][1] - q[1][1]) * gamma / 2.0,
        0.0,
    ];
    let sigma_minus = [-sigma_plus[0], -sigma_plus[1], 0.0];
    let ratio_of = |stress: [f64; 3]| {
        tsai_wu_strength_ratio(stress, &strength())
            .expect("valid strengths")
            .expect("finite ratio")
    };
    let (ratio_plus, ratio_minus) = (ratio_of(sigma_plus), ratio_of(sigma_minus));
    let expected = ratio_plus.min(ratio_minus);
    assert!(
        (ratio - expected).abs() <= 1e-9 * expected,
        "R = {ratio}, expected {expected} (+45°: {ratio_plus}, −45°: {ratio_minus})"
    );
    // The reported ply is the one carrying the smaller ratio.
    let reported_angle = lam.plies[index].angle_deg;
    if ratio_plus < ratio_minus {
        assert_eq!(reported_angle, 45.0);
    } else {
        assert_eq!(reported_angle, -45.0);
    }
}

#[test]
fn no_failure_is_reported_when_nothing_is_loaded_or_there_are_no_plies() {
    let lam = laminate(&[0.0, 90.0]);
    let strengths = vec![strength(); 2];
    assert!(
        first_ply_failure(&lam, &strengths, [0.0; 6])
            .expect("valid input")
            .is_none()
    );
    assert!(
        first_ply_failure(&Laminate::new(), &[], [1.0e3; 6])
            .expect("valid input")
            .is_none()
    );
}

#[test]
fn invalid_input_is_reported() {
    let lam = laminate(&[0.0, 90.0]);
    let loads = [1.0e3, 0.0, 0.0, 0.0, 0.0, 0.0];
    // Strength count mismatch (both directions).
    assert!(matches!(
        first_ply_failure(&lam, &[strength()], loads),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        first_ply_failure(&Laminate::new(), &[strength()], loads),
        Err(Error::InvalidInput(_))
    ));
    // Bad strength and non-finite loads.
    let strengths = vec![strength(); 2];
    let mut bad = strengths.clone();
    bad[1].yt = -1.0;
    assert!(matches!(
        first_ply_failure(&lam, &bad, loads),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        first_ply_failure(&lam, &strengths, [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0]),
        Err(Error::InvalidInput(_))
    ));
    // Non-physical plies: zero thickness, ν12·ν21 ≥ 1.
    let mut thin = laminate(&[0.0]);
    thin.plies[0].thickness = 0.0;
    assert!(matches!(
        first_ply_failure(&thin, &[strength()], loads),
        Err(Error::InvalidInput(_))
    ));
    let mut impossible = laminate(&[0.0]);
    impossible.plies[0].nu12 = 5.0;
    assert!(matches!(
        first_ply_failure(&impossible, &[strength()], loads),
        Err(Error::InvalidInput(_))
    ));
}
