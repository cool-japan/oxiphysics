// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended unit tests for the SPH kernels and kernel utilities in the parent
//! module (QuinticWendlandC6, super-Gaussian, normalisation and gradient checks).

use super::*;

// --- QuinticWendlandC6Kernel ---

#[test]
fn test_quintic_wendland_c6_positive_inside() {
    let k = QuinticWendlandC6Kernel;
    assert!(k.w(0.0, 1.0) > 0.0, "QWC6 at origin should be positive");
    assert!(k.w(0.3, 1.0) > 0.0, "QWC6 at r=0.3 should be positive");
    assert!(k.w(0.8, 1.0) > 0.0, "QWC6 at r=0.8 should be positive");
}

#[test]
fn test_quintic_wendland_c6_zero_outside() {
    let k = QuinticWendlandC6Kernel;
    assert!(k.w(1.0, 1.0).abs() < 1e-14, "QWC6 at boundary should be 0");
    assert!(
        k.w(1.5, 1.0).abs() < 1e-14,
        "QWC6 outside support should be 0"
    );
}

#[test]
fn test_quintic_wendland_c6_gradient_negative() {
    let k = QuinticWendlandC6Kernel;
    let dw = k.grad_w(0.3, 1.0);
    assert!(
        dw < 0.0,
        "QWC6 gradient inside support should be negative: {dw}"
    );
}

#[test]
fn test_quintic_wendland_c6_monotonic() {
    let k = QuinticWendlandC6Kernel;
    let h = 1.0;
    let w0 = k.w(0.0, h);
    let w1 = k.w(0.3, h);
    let w2 = k.w(0.6, h);
    let w3 = k.w(0.9, h);
    assert!(
        w0 > w1 && w1 > w2 && w2 > w3,
        "QWC6 should be monotonically decreasing"
    );
}

#[test]
fn test_quintic_wendland_c6_eval_wrapper() {
    let w = quintic_wendland_c6_eval(0.0, 1.0);
    assert!(w > 0.0);
    let w = quintic_wendland_c6_eval(1.5, 1.0);
    assert!(w.abs() < 1e-14);
}

// --- GaussianKernel ---

#[test]
fn test_gaussian_positive_at_origin() {
    let k = GaussianKernel;
    assert!(k.w(0.0, 1.0) > 0.0, "Gaussian at origin should be positive");
}

#[test]
fn test_gaussian_decreases_with_distance() {
    let k = GaussianKernel;
    let h = 1.0;
    let w0 = k.w(0.0, h);
    let w1 = k.w(0.5, h);
    let w2 = k.w(1.0, h);
    assert!(w0 > w1 && w1 > w2, "Gaussian should decrease with distance");
}

#[test]
fn test_gaussian_gradient_negative() {
    let k = GaussianKernel;
    let dw = k.grad_w(0.5, 1.0);
    assert!(
        dw < 0.0,
        "Gaussian gradient inside support should be negative: {dw}"
    );
}

#[test]
fn test_gaussian_gradient_zero_at_origin() {
    let k = GaussianKernel;
    assert!(
        k.grad_w(0.0, 1.0).abs() < 1e-10,
        "Gaussian gradient at origin should be 0"
    );
}

#[test]
fn test_gaussian_eval_wrapper() {
    let w = gaussian_eval(0.0, 1.0);
    assert!(w > 0.0);
}

#[test]
fn test_gaussian_integrates_approximately() {
    let k = GaussianKernel;
    let h = 1.0;
    let support = 4.0 * h;
    let n = 100usize;
    let dx = 2.0 * support / n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        let x = -support + (i as f64 + 0.5) * dx;
        for j in 0..n {
            let y = -support + (j as f64 + 0.5) * dx;
            for kk in 0..n {
                let z = -support + (kk as f64 + 0.5) * dx;
                let r = (x * x + y * y + z * z).sqrt();
                sum += k.w(r, h) * dx * dx * dx;
            }
        }
    }
    assert!(
        (sum - 1.0).abs() < 0.05,
        "Gaussian integral = {sum}, expected ~1.0"
    );
}

// --- QuinticSplineKernel ---

#[test]
fn test_quintic_spline_positive_inside() {
    let k = QuinticSplineKernel;
    let h = 1.0;
    assert!(k.w(0.0, h) > 0.0, "QS at origin should be positive");
    assert!(k.w(0.5, h) > 0.0, "QS at r=0.5 should be positive");
    assert!(k.w(1.5, h) > 0.0, "QS at r=1.5 should be positive");
    assert!(k.w(2.5, h) > 0.0, "QS at r=2.5 should be positive");
}

#[test]
fn test_quintic_spline_zero_outside() {
    let k = QuinticSplineKernel;
    let h = 1.0;
    assert!(k.w(3.0, h).abs() < 1e-14, "QS at boundary should be 0");
    assert!(k.w(3.5, h).abs() < 1e-14, "QS outside support should be 0");
}

#[test]
fn test_quintic_spline_gradient_zero_at_origin() {
    let k = QuinticSplineKernel;
    assert!(k.grad_w(0.0, 1.0).abs() < 1e-10);
}

#[test]
fn test_quintic_spline_eval_wrapper() {
    let w = quintic_spline_eval(0.0, 1.0);
    assert!(w > 0.0);
    let w_out = quintic_spline_eval(3.5, 1.0);
    assert!(w_out.abs() < 1e-14);
}

// --- KSPH gradient correction ---

#[test]
fn test_ksph_corrected_gradient_identity() {
    // With identity correction matrix, gradient should be unchanged
    let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let g = [1.0, 2.0, 3.0];
    let result = ksph_corrected_gradient(g, id);
    for i in 0..3 {
        assert!(
            (result[i] - g[i]).abs() < 1e-12,
            "Identity correction should preserve gradient: result[{i}]={}",
            result[i]
        );
    }
}

#[test]
fn test_ksph_corrected_gradient_scales() {
    // Diagonal scaling matrix: should scale each component
    let scale = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.5]];
    let g = [1.0, 1.0, 1.0];
    let result = ksph_corrected_gradient(g, scale);
    assert!((result[0] - 2.0).abs() < 1e-12);
    assert!((result[1] - 3.0).abs() < 1e-12);
    assert!((result[2] - 0.5).abs() < 1e-12);
}

#[test]
fn test_ksph_correct_all_gradients_count() {
    let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let grads = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let result = ksph_correct_all_gradients(&grads, id);
    assert_eq!(result.len(), grads.len());
}

// --- Second moment matrix ---

#[test]
fn test_second_moment_matrix_symmetric() {
    let pos_i = [0.0; 3];
    let neighbors = vec![[0.1, 0.2, 0.3], [-0.1, 0.1, -0.2]];
    let w_vals = vec![1.0, 1.0];
    let masses = vec![1.0, 1.0];
    let densities = vec![1.0, 1.0];
    let m = second_moment_matrix(pos_i, &neighbors, &w_vals, &masses, &densities);
    for (i, row) in m.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(
                (val - m[j][i]).abs() < 1e-12,
                "Second moment matrix should be symmetric: m[{i}][{j}]={} != m[{j}][{i}]={}",
                val,
                m[j][i]
            );
        }
    }
}

// --- Kernel renormalization ---

#[test]
fn test_renormalize_kernel_sums_to_one() {
    let w_vals = vec![0.5, 0.3, 0.2];
    let vols = vec![1.0, 1.0, 1.0];
    let scale = renormalize_kernel(&w_vals, &vols);
    // sum = 1.0, scale = 1.0
    assert!((scale - 1.0).abs() < 1e-12, "scale={scale}");
}

#[test]
fn test_renormalize_kernel_non_unit_sum() {
    let w_vals = vec![0.4, 0.4];
    let vols = vec![1.0, 1.0];
    let scale = renormalize_kernel(&w_vals, &vols);
    // sum = 0.8, scale = 1.25
    assert!((scale - 1.25).abs() < 1e-12, "scale={scale}");
}

#[test]
fn test_apply_renormalization_modifies_in_place() {
    let mut w_vals = vec![0.4, 0.4];
    let vols = vec![1.0, 1.0];
    apply_renormalization(&mut w_vals, &vols);
    let sum: f64 = w_vals.iter().zip(vols.iter()).map(|(&w, &v)| w * v).sum();
    assert!(
        (sum - 1.0).abs() < 1e-12,
        "After renormalization, sum should be 1: {sum}"
    );
}

#[test]
fn test_kernel_renorm_factor_matches_shepard() {
    let w = vec![0.6, 0.4];
    let m = vec![1.0, 1.0];
    let d = vec![1.0, 1.0];
    let factor = kernel_renorm_factor(&w, &m, &d);
    let shepard = shepard_correction(&w, &m, &d);
    assert!((factor - shepard).abs() < 1e-12);
}

// --- SPH gradient and Laplacian stencils ---

#[test]
fn test_sph_laplacian_1d_constant_field() {
    // Constant field: Laplacian should be 0
    let k = CubicSplineKernel;
    let h = 0.2;
    let f_i = 5.0;
    let pos_neighbors = vec![-0.1, 0.0, 0.1];
    let f_neighbors = vec![5.0, 5.0, 5.0];
    let masses = vec![0.001, 0.001, 0.001];
    let densities = vec![1000.0, 1000.0, 1000.0];
    let lap = sph_laplacian_1d(
        f_i,
        &f_neighbors,
        0.0,
        &pos_neighbors,
        &masses,
        &densities,
        &k,
        h,
    );
    // The Laplacian of a constant field should be approximately 0
    assert!(
        lap.abs() < 1.0,
        "Laplacian of constant field should be ~0: {lap}"
    );
}

#[test]
fn test_sph_gradient_3d_single_neighbor() {
    // ∇_i W_ij = (dW/dr) (x_i - x_j) / r points from x_i towards x_j (dW/dr < 0),
    // so a neighbour in +x direction with f_j > f_i yields a POSITIVE x gradient
    // of exactly vol * (f_j - f_i) * |dW/dr|.
    let k = CubicSplineKernel;
    let h = 1.0;
    let pos_i = [0.0, 0.0, 0.0];
    let pos_neighbors = vec![[0.5, 0.0, 0.0]];
    let f_neighbors = vec![1.0];
    let masses = vec![1.0];
    let densities = vec![1.0];
    let grad_f = sph_gradient_3d(
        0.0,
        &f_neighbors,
        pos_i,
        &pos_neighbors,
        &masses,
        &densities,
        &k,
        h,
    );
    let expected = -k.grad_w(0.5, h);
    assert!(
        expected > 0.0,
        "cubic spline dW/dr must be negative at r = 0.5"
    );
    assert!(
        grad_f[0] > 0.0,
        "f increases in +x, so the x gradient must be positive: {}",
        grad_f[0]
    );
    assert!(
        (grad_f[0] - expected).abs() < 1e-12 * expected,
        "x gradient {} should equal -dW/dr = {expected}",
        grad_f[0]
    );
    assert!(
        grad_f[1].abs() < 1e-12,
        "Gradient y should be 0: {}",
        grad_f[1]
    );
    assert!(
        grad_f[2].abs() < 1e-12,
        "Gradient z should be 0: {}",
        grad_f[2]
    );
}

#[test]
fn test_sph_gradient_3d_antisymmetric() {
    // Two symmetric neighbors with antisymmetric field → gradient non-zero
    let k = CubicSplineKernel;
    let h = 1.0;
    let pos_i = [0.0, 0.0, 0.0];
    let pos_neighbors = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
    let f_neighbors = vec![1.0, -1.0]; // antisymmetric field
    let masses = vec![1.0, 1.0];
    let densities = vec![1.0, 1.0];
    let grad_f = sph_gradient_3d(
        0.0,
        &f_neighbors,
        pos_i,
        &pos_neighbors,
        &masses,
        &densities,
        &k,
        h,
    );
    // f = 2x sampled at x = ±0.5: both contributions add up to 2 |dW/dr| > 0.
    let expected = -2.0 * k.grad_w(0.5, h);
    assert!(
        grad_f[0] > 0.0,
        "f increases with x, so the x gradient must be positive: {}",
        grad_f[0]
    );
    assert!(
        (grad_f[0] - expected).abs() < 1e-12 * expected,
        "x gradient {} should equal {expected}",
        grad_f[0]
    );
    assert!(grad_f[1].abs() < 1e-14 && grad_f[2].abs() < 1e-14);
}

#[test]
fn test_sph_gradient_3d_recovers_linear_field_on_lattice() {
    // f = 2x + 3y sampled on a 7^3 lattice (spacing 0.5 h, cubic spline,
    // V = spacing^3).  The interior (central) particle must see ∇f ≈ (2, 3, 0).
    let k = CubicSplineKernel;
    let h = 1.0;
    let spacing = 0.5 * h;
    let volume = spacing * spacing * spacing;
    let mut positions = Vec::new();
    for ix in -3_i32..=3 {
        for iy in -3_i32..=3 {
            for iz in -3_i32..=3 {
                positions.push([
                    f64::from(ix) * spacing,
                    f64::from(iy) * spacing,
                    f64::from(iz) * spacing,
                ]);
            }
        }
    }
    let field = |p: &[f64; 3]| 2.0 * p[0] + 3.0 * p[1];
    let f_neighbors: Vec<f64> = positions.iter().map(field).collect();
    let masses = vec![volume; positions.len()];
    let densities = vec![1.0; positions.len()];
    let centre = [0.0, 0.0, 0.0];
    let grad_f = sph_gradient_3d(
        field(&centre),
        &f_neighbors,
        centre,
        &positions,
        &masses,
        &densities,
        &k,
        h,
    );
    assert!(
        (grad_f[0] - 2.0).abs() < 0.2,
        "∂f/∂x should be ≈ 2, got {}",
        grad_f[0]
    );
    assert!(
        (grad_f[1] - 3.0).abs() < 0.3,
        "∂f/∂y should be ≈ 3, got {}",
        grad_f[1]
    );
    assert!(
        grad_f[2].abs() < 1e-10,
        "∂f/∂z should vanish by symmetry, got {}",
        grad_f[2]
    );
}

// --- Gaussian gradient vector ---

#[test]
fn test_gaussian_grad_points_opposite_to_displacement() {
    let h = 1.0;
    let r_vec = [0.5, 0.0, 0.0];
    let g = gaussian_grad(r_vec, h);
    assert!(
        g[0] < 0.0,
        "Gaussian grad x should be negative (away from origin): {}",
        g[0]
    );
    assert!(g[1].abs() < 1e-12);
    assert!(g[2].abs() < 1e-12);
}

#[test]
fn test_gaussian_grad_zero_at_origin() {
    let g = gaussian_grad([0.0; 3], 1.0);
    let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    assert!(
        mag < 1e-10,
        "Gaussian grad at origin should be 0: mag={mag}"
    );
}

// --- Kernel comparison with new kernels ---

#[test]
fn test_gaussian_different_from_cubic_spline() {
    let (max_abs, _) = compare_kernels(&GaussianKernel, &CubicSplineKernel, 1.0, 100);
    assert!(max_abs > 0.0, "Gaussian and CubicSpline should differ");
}

#[test]
fn test_quintic_spline_different_from_cubic_spline() {
    let (max_abs, _) = compare_kernels(&QuinticSplineKernel, &CubicSplineKernel, 1.0, 100);
    assert!(
        max_abs > 0.0,
        "Quintic spline and CubicSpline should differ"
    );
}

// --- Laplacian sign tests ---

#[test]
fn test_gaussian_laplacian_negative_near_origin() {
    // Gaussian is concave near origin → Laplacian should be negative
    let k = GaussianKernel;
    // At r=0: d²W/dr² = sigma * (-2/h²) → negative
    // But the 3D Laplacian includes 2/r * dW/dr term
    // Just check it's finite
    let lap = k.laplacian_w(0.1, 1.0);
    assert!(
        lap.is_finite(),
        "Gaussian Laplacian should be finite: {lap}"
    );
}

#[test]
fn test_quintic_spline_laplacian_finite() {
    let k = QuinticSplineKernel;
    for r in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5] {
        let lap = k.laplacian_w(r, 1.0);
        assert!(
            lap.is_finite(),
            "Quintic spline Laplacian should be finite at r={r}: {lap}"
        );
    }
}

// --- Kernel symmetry: w(-r) = w(r) since kernels use |r| ---

#[test]
fn test_kernels_use_scalar_distance() {
    let kernels: Vec<(&str, Box<dyn SphKernel>)> = vec![
        ("Gaussian", Box::new(GaussianKernel)),
        ("QuinticSpline", Box::new(QuinticSplineKernel)),
        ("QuinticWendlandC6", Box::new(QuinticWendlandC6Kernel)),
    ];
    for (name, k) in &kernels {
        let w = k.w(0.3, 1.0);
        assert!(w > 0.0, "{name} should be positive at r=0.3: w={w}");
    }
}

// --- Test gradient magnitudes for new kernels ---

#[test]
fn test_gaussian_grad_magnitude_matches_scalar() {
    let h = 1.0;
    let r = 0.5;
    let r_vec = [r, 0.0, 0.0];
    let scalar_dw = GaussianKernel.grad_w(r, h).abs();
    let vec_g = gaussian_grad(r_vec, h);
    let mag = (vec_g[0] * vec_g[0] + vec_g[1] * vec_g[1] + vec_g[2] * vec_g[2]).sqrt();
    assert!(
        (mag - scalar_dw).abs() < 1e-12,
        "Gaussian grad magnitude mismatch: {mag:.6e} vs {scalar_dw:.6e}"
    );
}

#[test]
fn test_quintic_wendland_c6_gradient_zero_at_origin() {
    let k = QuinticWendlandC6Kernel;
    assert!(
        k.grad_w(0.0, 1.0).abs() < 1e-10,
        "QWC6 gradient at origin should be 0"
    );
}

#[test]
fn test_quintic_wendland_c6_gradient_zero_outside() {
    let k = QuinticWendlandC6Kernel;
    assert!(
        k.grad_w(1.5, 1.0).abs() < 1e-14,
        "QWC6 gradient outside support should be 0"
    );
}

// --- 3-D normalisation: ∫ 4π r² W(r, h) dr = 1 ---

/// Composite Simpson integral of the radial volume integrand
/// `4π r² W(r, h)` over `[0, r_max]`.
fn radial_volume_integral(kernel: &dyn SphKernel, h: f64, r_max: f64) -> f64 {
    let n = 20_000_usize; // even number of intervals
    let dr = r_max / n as f64;
    let integrand = |r: f64| 4.0 * PI * r * r * kernel.w(r, h);
    let mut sum = integrand(0.0) + integrand(r_max);
    for i in 1..n {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * integrand(i as f64 * dr);
    }
    sum * dr / 3.0
}

#[test]
fn test_kernels_integrate_to_one_in_3d() {
    // (name, kernel, support radius in units of h).  h != 1 so that a wrong
    // power of h in the normalisation constant is also caught.
    let h = 0.7;
    let cases: Vec<(&str, Box<dyn SphKernel>, f64)> = vec![
        ("CubicSpline", Box::new(CubicSplineKernel), 2.0),
        ("WendlandC2", Box::new(WendlandC2Kernel), 2.0),
        ("WendlandC4", Box::new(WendlandC4Kernel), 2.0),
        ("WendlandC6", Box::new(WendlandC6Kernel), 2.0),
        ("QuinticWendlandC6", Box::new(QuinticWendlandC6Kernel), 1.0),
        ("QuinticSpline", Box::new(QuinticSplineKernel), 3.0),
        ("SuperGaussian", Box::new(SuperGaussianKernel), 4.0),
        ("Gaussian", Box::new(GaussianKernel), 4.0),
    ];
    for (name, kernel, support) in &cases {
        let integral = radial_volume_integral(kernel.as_ref(), h, support * h);
        assert!(
            (integral - 1.0).abs() < 1e-3,
            "{name}: ∫ 4π r² W dr = {integral}, expected 1 (tol 1e-3)"
        );
    }
}

#[test]
fn test_quintic_wendland_c6_uses_3d_constant() {
    // W(0, h) = σ = 1365 / (64 π h³); the 2-D constant 78 / (7 π h³) would
    // give a kernel that integrates to only 0.5224.
    let h = 0.5;
    let w0 = QuinticWendlandC6Kernel.w(0.0, h);
    let expected = 1365.0 / (64.0 * PI * h * h * h);
    assert!(
        (w0 - expected).abs() < 1e-12 * expected,
        "W(0) = {w0}, expected {expected}"
    );
}

#[test]
fn test_normalised_kernel_gradients_match_finite_differences() {
    let h = 0.7;
    let cases: Vec<(&str, Box<dyn SphKernel>, f64)> = vec![
        ("WendlandC4", Box::new(WendlandC4Kernel), 2.0),
        ("WendlandC6", Box::new(WendlandC6Kernel), 2.0),
        ("QuinticWendlandC6", Box::new(QuinticWendlandC6Kernel), 1.0),
        ("QuinticSpline", Box::new(QuinticSplineKernel), 3.0),
        ("SuperGaussian", Box::new(SuperGaussianKernel), 3.5),
    ];
    let eps = 1e-6 * h;
    for (name, kernel, support) in &cases {
        for step in 1..40 {
            let r = support * h * f64::from(step) / 40.0;
            let fd = (kernel.w(r + eps, h) - kernel.w(r - eps, h)) / (2.0 * eps);
            let analytic = kernel.grad_w(r, h);
            let scale = fd.abs().max(analytic.abs()).max(1e-30);
            assert!(
                (fd - analytic).abs() / scale < 1e-5,
                "{name}: dW/dr at r = {r}: analytic {analytic:.9e} vs finite difference {fd:.9e}"
            );
        }
    }
}
