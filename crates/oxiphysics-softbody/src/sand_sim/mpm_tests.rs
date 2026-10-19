// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the MLS-MPM sand kernel.
//!
//! They target the failure classes a hand-rolled MPM typically has: weights
//! that do not form a partition of unity (mass loss and spurious velocity
//! gradients), a decomposition that is not a decomposition, a return mapping
//! that drops the rotation or lands off the yield surface, a P2G/G2P pair
//! that leaks momentum, and a missing `Fᵀ` in the force.

use super::SandParticle;
use super::mpm::{
    IDENTITY, Mat3, Material, SparseGrid, Stencil, StepConstants, Vec3, decompose, floor_contact,
    is_active, kirchhoff_stress, mat_det, mat_from_col_major, mat_mul, mat_to_col_major,
    mat_transpose, mat_vec, project_deformation, project_log_strain, quadratic_weights,
    stable_time_step, substep, sym_eigen3,
};
use super::test_support::{Lcg, assert_close, assert_vec_close, max_abs_diff, rotation, sandwich};

fn symmetric(r: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| r[i][j] + r[j][i]))
}

fn default_material(friction_angle_deg: f64) -> Material {
    Material::new(3.537e5, 0.3, friction_angle_deg).expect("valid material")
}

fn frobenius(a: &Mat3) -> f64 {
    a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

fn sorted(mut v: [f64; 3]) -> [f64; 3] {
    v.sort_by(|a, b| a.total_cmp(b));
    v
}

// ---------------------------------------------------------------------------
// Symmetric eigen-decomposition
// ---------------------------------------------------------------------------

#[test]
fn eigen_reconstructs_random_symmetric_matrices() {
    let mut rng = Lcg::new(1);
    for _ in 0..400 {
        let a = symmetric(&rng.mat3(-3.0, 3.0));
        let (lambda, v) = sym_eigen3(&a);
        let vtv = mat_mul(&mat_transpose(&v), &v);
        assert!(
            max_abs_diff(&vtv, &IDENTITY) < 1e-13,
            "V is not orthonormal"
        );
        let back = sandwich(&v, lambda, &mat_transpose(&v));
        assert!(max_abs_diff(&back, &a) < 1e-12, "reconstruction error");
        assert_close(
            lambda.iter().sum(),
            a[0][0] + a[1][1] + a[2][2],
            1e-13,
            "trace",
        );
        assert_close(lambda.iter().product(), mat_det(&a), 1e-11, "determinant");
    }
}

#[test]
fn eigen_handles_repeated_and_degenerate_spectra() {
    // Repeated eigenvalue in a rotated basis: A = Q diag(2, 2, 5) Qᵀ.
    let q = rotation([1.0, 2.0, 3.0], 0.7);
    let a = sandwich(&q, [2.0, 2.0, 5.0], &mat_transpose(&q));
    let (lambda, v) = sym_eigen3(&a);
    let sorted_lambda = sorted(lambda);
    assert_vec_close(sorted_lambda, [2.0, 2.0, 5.0], 1e-13, "repeated spectrum");
    assert!(max_abs_diff(&sandwich(&v, lambda, &mat_transpose(&v)), &a) < 1e-13);

    // Already diagonal, identity, zero and rank one.
    for (m, expected) in [
        (
            [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]],
            [1.0, 2.0, 3.0],
        ),
        (IDENTITY, [1.0, 1.0, 1.0]),
        ([[0.0; 3]; 3], [0.0, 0.0, 0.0]),
        (
            [[1.0, 2.0, 2.0], [2.0, 4.0, 4.0], [2.0, 4.0, 4.0]],
            [0.0, 0.0, 9.0],
        ),
    ] {
        let (lambda, v) = sym_eigen3(&m);
        assert_vec_close(sorted(lambda), expected, 1e-13, "spectrum");
        assert!(max_abs_diff(&sandwich(&v, lambda, &mat_transpose(&v)), &m) < 1e-13);
    }

    // Off-diagonal entries far below the diagonal must not stall or corrupt.
    let nearly_diagonal = [[2.0, 1e-20, 0.0], [1e-20, 3.0, 1e-21], [0.0, 1e-21, 4.0]];
    let (lambda, _) = sym_eigen3(&nearly_diagonal);
    assert_vec_close(sorted(lambda), [2.0, 3.0, 4.0], 1e-14, "nearly diagonal");

    // Scale invariance: 1e8 and 1e-8 times a generic matrix.
    let base = symmetric(&Lcg::new(5).mat3(-1.0, 1.0));
    for scale in [1.0e8, 1.0e-8] {
        let a: Mat3 = std::array::from_fn(|i| std::array::from_fn(|j| base[i][j] * scale));
        let (lambda, v) = sym_eigen3(&a);
        let back = sandwich(&v, lambda, &mat_transpose(&v));
        assert!(max_abs_diff(&back, &a) < 1e-12 * scale, "scale {scale}");
    }
}

// ---------------------------------------------------------------------------
// Quadratic B-spline stencil
// ---------------------------------------------------------------------------

#[test]
fn weights_form_a_partition_of_unity_with_zero_first_moment() {
    // The defect class of an off-by-one stencil: weights summing to < 1 and a
    // non-zero first moment (which turns a uniform velocity into a velocity
    // gradient).
    let mut rng = Lcg::new(2);
    for step in 0..=1000 {
        let fx = 0.5 + f64::from(step) / 1000.0 * (1.5 - 0.5 - 1e-12);
        let w = quadratic_weights(fx);
        assert!(w.iter().all(|&x| x >= 0.0), "negative weight at {fx}");
        assert_close(w.iter().sum(), 1.0, 1e-14, "sum of weights");
        let first: f64 = w
            .iter()
            .enumerate()
            .map(|(k, wk)| wk * (k as f64 - fx))
            .sum();
        assert!(first.abs() < 1e-14, "first moment {first} at fx={fx}");
        let second: f64 = w
            .iter()
            .enumerate()
            .map(|(k, wk)| wk * (k as f64 - fx).powi(2))
            .sum();
        assert_close(second, 0.25, 1e-14, "second moment = 1/4 cell²");
    }
    // Random positions on grids of different resolution.
    for h in [0.05, 1.0, 0.0123] {
        for _ in 0..300 {
            let pos = rng.vec3(-5.0, 5.0);
            let stencil = Stencil::new(pos, 1.0 / h).expect("finite position");
            let mut total = 0.0;
            let mut moment = [0.0; 3];
            let mut second = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    for c in 0..3 {
                        let (w, offset) = stencil.node(a, b, c);
                        assert!(w >= 0.0);
                        total += w;
                        for ((m, o), row) in moment.iter_mut().zip(offset).zip(second.iter_mut()) {
                            *m += w * o;
                            for (entry, o2) in row.iter_mut().zip(offset) {
                                *entry += w * o * o2;
                            }
                        }
                        // The node really sits `offset` cells from the particle.
                        let node = [
                            (stencil.base[0] + a as i64) as f64 * h,
                            (stencil.base[1] + b as i64) as f64 * h,
                            (stencil.base[2] + c as i64) as f64 * h,
                        ];
                        for ((n, p), o) in node.iter().zip(pos).zip(offset) {
                            assert!((n - p - o * h).abs() < 1e-12, "node position");
                        }
                    }
                }
            }
            assert_close(total, 1.0, 1e-13, "27 weights");
            assert_vec_close(moment, [0.0; 3], 1e-13, "first moment");
            let quarter_identity = [[0.25, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.25]];
            assert!(
                max_abs_diff(&second, &quarter_identity) < 1e-13,
                "second moment"
            );
            for (b, x) in stencil.base.iter().zip(pos) {
                assert_eq!(*b, (x / h - 0.5).floor() as i64);
            }
            assert!(stencil.fx.iter().all(|f| (0.5..1.5).contains(f)));
        }
    }
}

#[test]
fn stencil_rejects_positions_it_cannot_represent() {
    let inv_h = 20.0;
    assert!(Stencil::new([f64::NAN, 0.0, 0.0], inv_h).is_none());
    assert!(Stencil::new([0.0, f64::INFINITY, 0.0], inv_h).is_none());
    assert!(Stencil::new([0.0, 0.0, f64::NEG_INFINITY], inv_h).is_none());
    // More than 2^20 cells from the origin.
    assert!(Stencil::new([1.0e6, 0.0, 0.0], inv_h).is_none());
    assert!(Stencil::new([-1.0e6, 0.0, 0.0], inv_h).is_none());
    assert!(Stencil::new([1.0e4, -1.0e4, 5.0e4], inv_h).is_some());
}

#[test]
fn sparse_grid_registers_each_node_once() {
    let h = 0.1;
    let stencil_at = |x: f64, y: f64, z: f64| Stencil::new([x, y, z], 1.0 / h).expect("finite");
    let mut grid = SparseGrid::with_capacity(8);
    let first = grid.register(&stencil_at(0.31, 0.42, 0.53));
    assert_eq!(grid.len(), 27);
    // Same cell again (cached block) and a second particle in the same cell.
    assert_eq!(grid.register(&stencil_at(0.31, 0.42, 0.53)), first);
    assert_eq!(grid.register(&stencil_at(0.33, 0.44, 0.52)), first);
    assert_eq!(grid.len(), 27);
    // One cell further along x shares 18 of the 27 nodes ...
    let shifted = grid.register(&stencil_at(0.41, 0.42, 0.53));
    assert_eq!(grid.len(), 36);
    assert_ne!(shifted, first);
    // ... and returning to the first block still gives the first ids.
    assert_eq!(grid.register(&stencil_at(0.31, 0.42, 0.53)), first);
    assert_eq!(grid.len(), 36);
    // Ids index the node arrays consistently.
    assert!(first.iter().chain(&shifted).all(|&id| id < grid.len()));
    assert_eq!(grid.cells.len(), grid.mass.len());
    assert_eq!(grid.cells.len(), grid.vel.len());
}

#[test]
fn sparse_grid_keys_are_unique_over_negative_and_positive_cells() {
    // Stencils with bases −3..=3 on every axis cover the node block −3..=5.
    let mut grid = SparseGrid::with_capacity(16);
    let h = 0.5;
    for ix in -3i64..=3 {
        for iy in -3i64..=3 {
            for iz in -3i64..=3 {
                // fx = 1.0 puts the particle exactly on cell `base + 1`.
                let pos = [
                    (ix as f64 + 1.0) * h,
                    (iy as f64 + 1.0) * h,
                    (iz as f64 + 1.0) * h,
                ];
                let stencil = Stencil::new(pos, 1.0 / h).expect("finite");
                assert_eq!(stencil.base, [ix, iy, iz]);
                let _ = grid.register(&stencil);
            }
        }
    }
    assert_eq!(grid.len(), 9 * 9 * 9);
    let mut cells = grid.cells.clone();
    cells.sort_unstable();
    cells.dedup();
    assert_eq!(cells.len(), 9 * 9 * 9, "cells must be distinct");
}

// ---------------------------------------------------------------------------
// Material constants
// ---------------------------------------------------------------------------

#[test]
fn material_constants_match_the_closed_forms() {
    let m = Material::new(3.537e5, 0.3, 30.0).expect("valid");
    assert_close(m.mu, 3.537e5 / 2.6, 1e-13, "mu");
    assert_close(m.lambda, 3.537e5 * 0.3 / (1.3 * 0.4), 1e-13, "lambda");
    assert_close(
        m.p_wave_modulus(),
        m.lambda + 2.0 * m.mu,
        1e-15,
        "P-wave modulus",
    );
    // sin 30° = 1/2, so alpha = sqrt(2/3) * 2 * 0.5 / 2.5.
    assert_close(m.alpha, (2.0f64 / 3.0).sqrt() * 0.4, 1e-14, "alpha(30°)");
    assert_eq!(Material::new(1.0e5, 0.2, 0.0).map(|m| m.alpha), Some(0.0));

    let invalid = [
        (0.0, 0.3, 30.0),
        (-1.0, 0.3, 30.0),
        (f64::NAN, 0.3, 30.0),
        (f64::INFINITY, 0.3, 30.0),
        (1.0e5, 0.5, 30.0),
        (1.0e5, -0.1, 30.0),
        (1.0e5, f64::NAN, 30.0),
        (1.0e5, 0.3, 90.0),
        (1.0e5, 0.3, -1.0),
        (1.0e5, 0.3, f64::NAN),
    ];
    for (e, nu, phi) in invalid {
        assert!(
            Material::new(e, nu, phi).is_none(),
            "E={e} nu={nu} phi={phi}"
        );
    }
}

// ---------------------------------------------------------------------------
// Hencky stress
// ---------------------------------------------------------------------------

#[test]
fn stress_is_zero_for_rigid_motions_and_follows_linear_elasticity_for_small_strain() {
    let mat = default_material(30.0);
    let zero = kirchhoff_stress(&IDENTITY, &mat);
    assert!(frobenius(&zero) == 0.0);
    // A pure rotation carries no stress.
    let q = rotation([0.3, -1.0, 0.5], 1.1);
    assert!(frobenius(&kirchhoff_stress(&q, &mat)) < 1e-9);
    // Small uniaxial stretch: τ_xx = (λ + 2μ) e, τ_yy = τ_zz = λ e.
    let e = 1.0e-7;
    let tau = kirchhoff_stress(
        &[[1.0 + e, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &mat,
    );
    assert_close(tau[0][0] / e, mat.lambda + 2.0 * mat.mu, 1e-6, "tau_xx");
    assert_close(tau[1][1] / e, mat.lambda, 1e-6, "tau_yy");
    assert_close(tau[2][2] / e, mat.lambda, 1e-6, "tau_zz");
    assert!(tau[0][1].abs() < 1e-3 && tau[1][2].abs() < 1e-3);
}

#[test]
fn stress_is_symmetric_and_objective() {
    // τ is symmetric (zero net torque on the grid: the "missing Fᵀ" class) and
    // transforms as τ(R F) = R τ(F) Rᵀ under a rigid rotation R.
    let mat = default_material(30.0);
    let mut rng = Lcg::new(11);
    for _ in 0..200 {
        let f: Mat3 = std::array::from_fn(|i| {
            std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 } + rng.range(-0.25, 0.25))
        });
        let tau = kirchhoff_stress(&f, &mat);
        let asym: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .map(|(i, j)| (tau[i][j] - tau[j][i]).abs())
            .fold(0.0, f64::max);
        assert!(asym < 1e-8 * (1.0 + frobenius(&tau)), "asymmetry {asym}");
        let r = rotation(rng.vec3(-1.0, 1.0), rng.range(-3.0, 3.0));
        let rotated = kirchhoff_stress(&mat_mul(&r, &f), &mat);
        let expected = mat_mul(&mat_mul(&r, &tau), &mat_transpose(&r));
        assert!(
            max_abs_diff(&rotated, &expected) < 1e-8 * (1.0 + frobenius(&tau)),
            "objectivity"
        );
    }
}

#[test]
fn decomposition_recovers_the_principal_stretches() {
    let mut rng = Lcg::new(3);
    for _ in 0..200 {
        let q = rotation(rng.vec3(-1.0, 1.0), rng.range(-3.0, 3.0));
        let r2 = rotation(rng.vec3(-1.0, 1.0), rng.range(-3.0, 3.0));
        let stretches = [
            rng.range(0.4, 1.8),
            rng.range(0.4, 1.8),
            rng.range(0.4, 1.8),
        ];
        let f = sandwich(&q, stretches, &mat_transpose(&r2));
        let dec = decompose(&f);
        assert_vec_close(sorted(dec.sigma), sorted(stretches), 1e-12, "sigma");
        for (s, e) in dec.sigma.iter().zip(dec.eps) {
            assert_close(e, s.ln(), 1e-14, "eps = ln sigma");
        }
        // F Fᵀ = U Σ² Uᵀ.
        let b = mat_mul(&f, &mat_transpose(&f));
        let sigma_sq = dec.sigma.map(|s| s * s);
        let back = sandwich(&dec.u, sigma_sq, &mat_transpose(&dec.u));
        assert!(max_abs_diff(&back, &b) < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Return mapping (Klár et al. 2016)
// ---------------------------------------------------------------------------

/// `K = (3λ + 2μ) α / (2μ)`: slope of the cone in (|tr ε|, ‖ε̂‖) space.
fn cone_slope(mat: &Material) -> f64 {
    (3.0 * mat.lambda + 2.0 * mat.mu) / (2.0 * mat.mu) * mat.alpha
}

fn deviator_norm(eps: [f64; 3]) -> f64 {
    let mean = (eps[0] + eps[1] + eps[2]) / 3.0;
    eps.iter().map(|e| (e - mean).powi(2)).sum::<f64>().sqrt()
}

#[test]
fn tension_keeps_the_rotation_and_drops_the_stretch() {
    // The defect class "returns the identity and loses the rotation".
    let mat = default_material(30.0);
    let q = rotation([1.0, 2.0, 3.0], 0.7);
    let r2 = rotation([-2.0, 1.0, 0.5], 0.4);
    let f = sandwich(&q, [1.2, 1.1, 1.05], &mat_transpose(&r2));
    let (projected, plastic_volume) = project_deformation(&f, &mat);
    let rotation_part = mat_mul(&q, &mat_transpose(&r2));
    assert!(
        max_abs_diff(&projected, &rotation_part) < 1e-12,
        "F' must be U Vᵀ"
    );
    assert_close(mat_det(&projected), 1.0, 1e-12, "det");
    assert_close(
        plastic_volume,
        (1.2f64 * 1.1 * 1.05).ln(),
        1e-12,
        "dilation",
    );
    assert!(
        frobenius(&kirchhoff_stress(&projected, &mat)) < 1e-9,
        "stress free"
    );
    // A tension state with some compression direction is still a cut-off.
    let mixed = sandwich(&q, [1.4, 0.95, 1.02], &mat_transpose(&r2));
    let (mixed_projected, mixed_volume) = project_deformation(&mixed, &mat);
    assert!(max_abs_diff(&mixed_projected, &rotation_part) < 1e-12);
    assert!(mixed_volume > 0.0);
}

#[test]
fn elastic_states_are_returned_unchanged() {
    let mat = default_material(30.0);
    let q = rotation([0.2, 0.9, -0.4], 1.3);
    let r2 = rotation([1.0, 0.0, 1.0], -0.6);
    // Pure hydrostatic compression, mild shear inside the cone, no strain.
    for stretches in [[0.9, 0.9, 0.9], [0.8, 0.82, 0.78], [1.0, 1.0, 1.0]] {
        let f = sandwich(&q, stretches, &mat_transpose(&r2));
        let (projected, plastic_volume) = project_deformation(&f, &mat);
        assert_eq!(projected, f, "elastic state must be bit-for-bit unchanged");
        assert_eq!(plastic_volume, 0.0);
    }
}

#[test]
fn over_sheared_compression_lands_on_the_cone_and_keeps_trace_and_direction() {
    for phi in [15.0, 30.0, 45.0, 60.0] {
        let mat = default_material(phi);
        let k = cone_slope(&mat);
        let eps = [-0.30, 0.05, 0.0];
        assert!(deviator_norm(eps) > k * (-(eps[0] + eps[1] + eps[2])) || phi >= 45.0);
        let proj = project_log_strain(eps, &mat);
        let trace: f64 = eps.iter().sum();
        if proj.yielded {
            // tr ε preserved, ‖ε̂‖ = K |tr ε| (f = 0), ε̂ direction preserved.
            assert_close(proj.strain.iter().sum(), trace, 1e-14, "trace");
            assert_close(
                deviator_norm(proj.strain),
                k * trace.abs(),
                1e-12,
                "on the cone",
            );
            let mean = trace / 3.0;
            let dir_in: Vec3 = eps.map(|e| e - mean);
            let dir_out: Vec3 = proj.strain.map(|e| e - mean);
            let cos = dir_in.iter().zip(dir_out).map(|(a, b)| a * b).sum::<f64>()
                / (deviator_norm(eps) * deviator_norm(proj.strain));
            assert_close(cos, 1.0, 1e-12, "direction");
            assert_eq!(proj.plastic_volume, 0.0, "the flow rule is isochoric");
        } else {
            assert_eq!(proj.strain, eps);
        }
    }
    // Idempotent: projecting the projected state changes nothing.
    let mat = default_material(30.0);
    let once = project_log_strain([-0.2, 0.15, 0.05], &mat);
    assert!(once.yielded);
    let twice = project_log_strain(once.strain, &mat);
    for (a, b) in once.strain.iter().zip(twice.strain) {
        assert_close(*a, b, 1e-12, "idempotent");
    }
}

#[test]
fn projection_of_a_deformation_gradient_realises_the_projected_strain() {
    // F' = U exp(H) Vᵀ must have principal log strains H and keep the principal
    // directions of F.
    let mat = default_material(30.0);
    let mut rng = Lcg::new(17);
    let mut yielded = 0;
    for _ in 0..200 {
        let q = rotation(rng.vec3(-1.0, 1.0), rng.range(-3.0, 3.0));
        let r2 = rotation(rng.vec3(-1.0, 1.0), rng.range(-3.0, 3.0));
        let stretches = [
            rng.range(0.6, 1.25),
            rng.range(0.6, 1.25),
            rng.range(0.6, 1.25),
        ];
        let f = sandwich(&q, stretches, &mat_transpose(&r2));
        let (projected, plastic_volume) = project_deformation(&f, &mat);
        let eps = stretches.map(f64::ln);
        let expected = project_log_strain(eps, &mat);
        if !expected.yielded {
            assert_eq!(projected, f);
            continue;
        }
        yielded += 1;
        assert_close(
            plastic_volume,
            expected.plastic_volume,
            1e-12,
            "plastic volume",
        );
        let h_stretches = expected.strain.map(f64::exp);
        let target = sandwich(&q, h_stretches, &mat_transpose(&r2));
        assert!(
            max_abs_diff(&projected, &target) < 1e-10,
            "F' = U exp(H) Vᵀ"
        );
    }
    assert!(yielded > 30, "too few yielding samples: {yielded}");
}

#[test]
fn triaxial_compression_at_yield_mobilises_exactly_the_friction_angle() {
    // The cone is calibrated to Mohr-Coulomb in triaxial compression:
    // (σ_a − σ_l) / (σ_a + σ_l) = sin φ with compressive principal stresses.
    for phi in [15.0, 30.0, 45.0] {
        let mat = default_material(phi);
        for (e_axial, e_lateral) in [(-0.30, 0.10), (-0.20, 0.06), (-0.12, 0.04)] {
            let f = [
                [f64::exp(e_axial), 0.0, 0.0],
                [0.0, f64::exp(e_lateral), 0.0],
                [0.0, 0.0, f64::exp(e_lateral)],
            ];
            let (projected, _) = project_deformation(&f, &mat);
            assert_ne!(
                projected, f,
                "the trial state must be over-sheared (phi={phi})"
            );
            let tau = kirchhoff_stress(&projected, &mat);
            // Principal axes of this state are the coordinate axes.
            let (sigma_a, sigma_l) = (-tau[0][0], -tau[1][1]);
            assert_close(tau[1][1], tau[2][2], 1e-9, "lateral symmetry");
            assert!(sigma_a > 0.0 && sigma_l > 0.0, "compressive state");
            let mobilised = (sigma_a - sigma_l) / (sigma_a + sigma_l);
            assert_close(
                mobilised,
                phi.to_radians().sin(),
                1e-9,
                &format!("mobilised friction at phi={phi}"),
            );
        }
        // Inside the cone the same state is elastic and mobilises less.
        let mild = [
            [f64::exp(-0.05), 0.0, 0.0],
            [0.0, f64::exp(-0.04), 0.0],
            [0.0, 0.0, f64::exp(-0.04)],
        ];
        let (unchanged, _) = project_deformation(&mild, &mat);
        assert_eq!(unchanged, mild);
        let tau = kirchhoff_stress(&mild, &mat);
        let mobilised = (tau[1][1] - tau[0][0]) / (-tau[0][0] - tau[1][1]);
        assert!(mobilised < phi.to_radians().sin());
    }
}

#[test]
fn zero_friction_sand_loses_all_shear_resistance() {
    let mat = default_material(0.0);
    let proj = project_log_strain([-0.1, 0.02, 0.02], &mat);
    assert!(proj.yielded);
    let mean = (-0.1 + 0.02 + 0.02) / 3.0;
    for e in proj.strain {
        assert_close(e, mean, 1e-14, "hydrostatic after projection");
    }
}

// ---------------------------------------------------------------------------
// Floor contact
// ---------------------------------------------------------------------------

#[test]
fn floor_contact_clamps_normal_velocity_and_applies_coulomb_friction() {
    // Separating nodes are untouched.
    let mut v = [1.0, 0.5, -2.0];
    floor_contact(&mut v, 0.5);
    assert_eq!(v, [1.0, 0.5, -2.0]);
    // Frictionless: only the normal component goes.
    let mut v = [1.0, -2.0, 0.5];
    floor_contact(&mut v, 0.0);
    assert_eq!(v, [1.0, 0.0, 0.5]);
    // Slip: tangential speed 3 loses µ·|Δv_n| = 0.5·2 = 1 along its direction.
    let mut v = [3.0, -2.0, 0.0];
    floor_contact(&mut v, 0.5);
    assert_vec_close(v, [2.0, 0.0, 0.0], 1e-15, "slip");
    let mut v = [0.0, -2.0, -4.0];
    floor_contact(&mut v, 0.5);
    assert_vec_close(v, [0.0, 0.0, -3.0], 1e-15, "slip along z");
    // Stick: friction impulse exceeds the tangential momentum.
    let mut v = [0.5, -2.0, 0.5];
    floor_contact(&mut v, 1.0);
    assert_eq!(v, [0.0, 0.0, 0.0]);
    // Sticky floor.
    let mut v = [7.0, -1.0e-9, 3.0];
    floor_contact(&mut v, f64::INFINITY);
    assert_eq!(v, [0.0, 0.0, 0.0]);
    // Exactly resting nodes get no spurious friction (and no NaN from ∞ · 0).
    let mut v = [1.0, 0.0, 0.0];
    floor_contact(&mut v, f64::INFINITY);
    assert_eq!(v, [1.0, 0.0, 0.0]);
}

// ---------------------------------------------------------------------------
// Particle / grid transfers
// ---------------------------------------------------------------------------

fn free_constants(h: f64, dt: f64) -> StepConstants {
    StepConstants {
        h,
        dt,
        gravity: [0.0; 3],
        material: default_material(30.0),
        floor_friction: 1.0,
        plane_strain_z: None,
    }
}

/// Random stress-free particles well above the floor with random velocities and
/// affine matrices.
fn flying_cloud(n: usize, seed: u64) -> Vec<SandParticle> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|_| {
            let position = [
                rng.range(-0.3, 0.3),
                rng.range(0.8, 1.4),
                rng.range(-0.3, 0.3),
            ];
            let mut p = SandParticle::new(
                position,
                rng.range(0.5e-3, 2.0e-3),
                rng.range(1.0e-6, 2.0e-6),
            );
            p.velocity = rng.vec3(-2.0, 2.0);
            p.affine_velocity = mat_to_col_major(&rng.mat3(-3.0, 3.0));
            p
        })
        .collect()
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Linear momentum Σ m v.
fn momentum(ps: &[SandParticle]) -> Vec3 {
    ps.iter().fold([0.0; 3], |acc, p| {
        [
            acc[0] + p.mass * p.velocity[0],
            acc[1] + p.mass * p.velocity[1],
            acc[2] + p.mass * p.velocity[2],
        ]
    })
}

/// APIC angular momentum Σ m (x × v + s), `s_a = (h²/4) ε_abc C_cb`.
fn angular_momentum(ps: &[SandParticle], h: f64) -> Vec3 {
    ps.iter().fold([0.0; 3], |acc, p| {
        let c = mat_from_col_major(&p.affine_velocity);
        let spin = [
            0.25 * h * h * (c[2][1] - c[1][2]),
            0.25 * h * h * (c[0][2] - c[2][0]),
            0.25 * h * h * (c[1][0] - c[0][1]),
        ];
        let l = cross(p.position, p.velocity);
        [
            acc[0] + p.mass * (l[0] + spin[0]),
            acc[1] + p.mass * (l[1] + spin[1]),
            acc[2] + p.mass * (l[2] + spin[2]),
        ]
    })
}

#[test]
fn transfers_conserve_mass_linear_and_angular_momentum() {
    let h = 0.05;
    let constants = free_constants(h, 1.0e-4);
    let mut ps = flying_cloud(150, 21);
    let p0 = momentum(&ps);
    let l0 = angular_momentum(&ps, h);
    let mass0: f64 = ps.iter().map(|p| p.mass).sum();
    substep(&mut ps, &constants);
    let p1 = momentum(&ps);
    let l1 = angular_momentum(&ps, h);
    assert_vec_close(p1, p0, 1e-12, "linear momentum");
    assert_vec_close(l1, l0, 1e-12, "angular momentum (APIC)");
    assert_close(ps.iter().map(|p| p.mass).sum(), mass0, 1e-15, "mass");
    // F = I means no stress in the first sub-step: particles only translate.
    assert!(ps.iter().all(|p| p.position.iter().all(|x| x.is_finite())));
    // A second sub-step keeps conserving (stress now acts, gravity is off):
    // linear momentum stays exact, since the internal forces are balanced.
    substep(&mut ps, &constants);
    assert_vec_close(momentum(&ps), p0, 1e-9, "linear momentum, second step");
}

#[test]
fn affine_velocity_fields_survive_a_round_trip_exactly() {
    // v_p = v0 + A x_p with C_p = A is reproduced exactly by P2G → G2P: the
    // test of the weights' moments and of D⁻¹ = 4/h².
    let h = 0.04;
    let constants = free_constants(h, 1.0e-4);
    let mut rng = Lcg::new(33);
    let v0 = [0.3, -0.7, 1.1];
    let a: Mat3 = rng.mat3(-2.0, 2.0);
    let mut ps = flying_cloud(120, 34);
    for p in &mut ps {
        let ax = mat_vec(&a, p.position);
        p.velocity = [v0[0] + ax[0], v0[1] + ax[1], v0[2] + ax[2]];
        p.affine_velocity = mat_to_col_major(&a);
    }
    let before: Vec<Vec3> = ps.iter().map(|p| p.position).collect();
    substep(&mut ps, &constants);
    for (p, x0) in ps.iter().zip(&before) {
        let ax = mat_vec(&a, *x0);
        assert_vec_close(
            p.velocity,
            [v0[0] + ax[0], v0[1] + ax[1], v0[2] + ax[2]],
            1e-11,
            "affine velocity",
        );
        let c = mat_from_col_major(&p.affine_velocity);
        assert!(max_abs_diff(&c, &a) < 1e-10, "affine matrix C");
    }
}

#[test]
fn deformation_gradient_is_advanced_by_the_velocity_gradient() {
    // After one sub-step F' = (I + dt C') F with C' = ∇v of the affine field.
    // Compression with a little shear stays inside the cone, so the update is
    // purely elastic and F' = I + dt A exactly. (Regression for a wrongly scaled
    // update, which changes how fast stress builds up.)
    let h = 0.04;
    let dt = 1.0e-3;
    let constants = free_constants(h, dt);
    let a: Mat3 = [[-0.5, 0.2, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, -0.5]];
    let mut ps = flying_cloud(60, 41);
    for p in &mut ps {
        p.velocity = mat_vec(&a, p.position);
        p.affine_velocity = mat_to_col_major(&a);
    }
    substep(&mut ps, &constants);
    let expected: Mat3 =
        std::array::from_fn(|i| std::array::from_fn(|j| IDENTITY[i][j] + dt * a[i][j]));
    for p in &ps {
        let f = mat_from_col_major(&p.deformation_gradient);
        assert!(max_abs_diff(&f, &expected) < 1e-12, "F' = I + dt A");
        assert_eq!(p.plastic_strain_vol, 0.0);
    }
}

#[test]
fn a_block_in_free_fall_accelerates_at_g_without_damping() {
    // The defect class of a stencil that loses weight: free fall would be
    // damped to a terminal velocity.
    let h = 0.05;
    let g = [0.0, -9.81, 0.0];
    let dt = 1.0e-3;
    let mut constants = free_constants(h, dt);
    constants.gravity = g;
    let mut ps: Vec<SandParticle> = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let position = [
                    0.2 + 0.027 * i as f64,
                    1.0 + 0.027 * j as f64,
                    -0.3 + 0.027 * k as f64,
                ];
                ps.push(SandParticle::new(position, 1.0e-3, 6.0e-7));
            }
        }
    }
    let start = ps.clone();
    let steps = 100;
    for _ in 0..steps {
        substep(&mut ps, &constants);
    }
    let n = steps as f64;
    for (p, p0) in ps.iter().zip(&start) {
        assert_vec_close(p.velocity, [0.0, g[1] * dt * n, 0.0], 1e-12, "velocity");
        // Semi-implicit Euler: y_n = y_0 + g dt² n (n + 1) / 2.
        assert_close(
            p.position[1],
            p0.position[1] + g[1] * dt * dt * n * (n + 1.0) / 2.0,
            1e-12,
            "height",
        );
        assert_close(p.position[0], p0.position[0], 1e-14, "x");
        // Pure translation: no velocity gradient and no strain.
        assert!(p.affine_velocity.iter().all(|c| c.abs() < 1e-9));
        let f = mat_from_col_major(&p.deformation_gradient);
        assert!(max_abs_diff(&f, &IDENTITY) < 1e-9, "no strain in free fall");
        assert_eq!(p.plastic_strain_vol, 0.0);
    }
}

#[test]
fn substep_is_bitwise_deterministic() {
    let constants = free_constants(0.05, 2.0e-4);
    let mut first = flying_cloud(200, 5);
    let mut second = first.clone();
    for _ in 0..5 {
        substep(&mut first, &constants);
        substep(&mut second, &constants);
    }
    assert_eq!(first, second);
}

#[test]
fn contact_with_the_floor_stops_a_falling_particle_and_removes_normal_momentum() {
    // A block falling onto the floor must stop (not tunnel): after contact the
    // lowest particles have no downward velocity and stay on or above y = 0.
    let h = 0.05;
    let mut constants = free_constants(h, 5.0e-4);
    constants.gravity = [0.0, -9.81, 0.0];
    let mut ps: Vec<SandParticle> = Vec::new();
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                let mut p = SandParticle::new(
                    [0.025 * i as f64, 0.1 + 0.025 * j as f64, 0.025 * k as f64],
                    1600.0 * 1.5625e-5,
                    1.5625e-5,
                );
                p.velocity = [0.0, -2.0, 0.0];
                ps.push(p);
            }
        }
    }
    for _ in 0..400 {
        substep(&mut ps, &constants);
    }
    assert!(
        ps.iter().all(|p| p.position[1] >= 0.0),
        "tunnelled through the floor"
    );
    let lowest = ps
        .iter()
        .map(|p| p.position[1])
        .fold(f64::INFINITY, f64::min);
    assert!(
        lowest < 0.05,
        "block should have reached the floor: {lowest}"
    );
    let kinetic: f64 = ps.iter().map(SandParticle::kinetic_energy).sum();
    let initial = ps.iter().map(|p| 0.5 * p.mass * 4.0).sum::<f64>();
    assert!(
        kinetic < 0.5 * initial,
        "impact must dissipate energy: {kinetic} vs {initial}"
    );
}

// ---------------------------------------------------------------------------
// Activity and time-step limit
// ---------------------------------------------------------------------------

#[test]
fn only_well_formed_particles_are_active() {
    let inv_h = 20.0;
    let good = SandParticle::new([0.1, 0.2, 0.3], 1.0e-3, 1.0e-6);
    assert!(is_active(&good, inv_h));
    type Mutation = Box<dyn Fn(&mut SandParticle)>;
    let mutations: Vec<(&str, Mutation)> = vec![
        ("zero mass", Box::new(|p| p.mass = 0.0)),
        ("negative mass", Box::new(|p| p.mass = -1.0)),
        ("NaN mass", Box::new(|p| p.mass = f64::NAN)),
        ("zero volume", Box::new(|p| p.volume = 0.0)),
        ("negative volume", Box::new(|p| p.volume = -1.0e-6)),
        ("NaN position", Box::new(|p| p.position[1] = f64::NAN)),
        (
            "infinite position",
            Box::new(|p| p.position[0] = f64::INFINITY),
        ),
        ("far position", Box::new(|p| p.position[2] = 1.0e9)),
        ("NaN velocity", Box::new(|p| p.velocity[2] = f64::NAN)),
        ("NaN F", Box::new(|p| p.deformation_gradient[4] = f64::NAN)),
        (
            "infinite C",
            Box::new(|p| p.affine_velocity[7] = f64::INFINITY),
        ),
        (
            "NaN plastic strain",
            Box::new(|p| p.plastic_strain_vol = f64::NAN),
        ),
    ];
    for (what, mutate) in mutations {
        let mut p = good.clone();
        mutate(&mut p);
        assert!(!is_active(&p, inv_h), "{what} must be inert");
    }
}

#[test]
fn stable_time_step_follows_the_wave_speed_and_the_fastest_particle() {
    let mat = default_material(30.0);
    let h = 0.05;
    let cfl = 0.3;
    let mut ps = vec![
        SandParticle::new([0.0, 1.0, 0.0], 1.6e-3, 1.0e-6),
        SandParticle::new([0.1, 1.0, 0.0], 3.2e-3, 1.0e-6),
    ];
    let c_min_density = (mat.p_wave_modulus() / 1600.0).sqrt();
    assert_close(
        stable_time_step(&ps, h, cfl, &mat),
        cfl * h / c_min_density,
        1e-13,
        "at rest",
    );
    ps[1].velocity = [3.0, 4.0, 0.0];
    assert_close(
        stable_time_step(&ps, h, cfl, &mat),
        cfl * h / (c_min_density + 5.0),
        1e-13,
        "with a 5 m/s particle",
    );
    // Inert particles do not count; with none left the step is unbounded.
    ps[0].mass = f64::NAN;
    ps[1].volume = 0.0;
    assert!(stable_time_step(&ps, h, cfl, &mat).is_infinite());
    assert!(stable_time_step(&[], h, cfl, &mat).is_infinite());
}
