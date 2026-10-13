// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the scientific visualization primitives.

use super::*;

// ---- TensorGlyph ----

#[test]
fn tensor_glyph_ellipsoid_axes_match_eigenvalues() {
    let evals = [3.0, 2.0, 1.0];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
    let axes = g.ellipsoid_axes();
    assert!((axes[0] - 3.0).abs() < 1e-10);
    assert!((axes[1] - 2.0).abs() < 1e-10);
    assert!((axes[2] - 1.0).abs() < 1e-10);
}

#[test]
fn tensor_glyph_scale_applied_to_axes() {
    let evals = [2.0, 1.0, 0.5];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 2.0);
    let axes = g.ellipsoid_axes();
    assert!((axes[0] - 4.0).abs() < 1e-10);
}

#[test]
fn tensor_glyph_isotropic_fa_is_zero() {
    let evals = [1.0, 1.0, 1.0];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
    assert!(g.fractional_anisotropy() < 1e-10);
}

#[test]
fn tensor_glyph_anisotropic_fa_near_one() {
    let evals = [100.0, 0.01, 0.01];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
    let fa = g.fractional_anisotropy();
    assert!(
        fa > 0.9,
        "FA={fa:.6} should be close to 1 for very anisotropic tensor"
    );
}

#[test]
fn tensor_glyph_surface_points_count() {
    let evals = [2.0, 1.5, 1.0];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
    let pts = g.ellipsoid_surface_points(16, 8);
    assert_eq!(pts.len(), 16 * 8);
}

#[test]
fn build_tensor_glyphs_from_diagonal_tensor() {
    let tensors = vec![[4.0, 1.0, 0.25, 0.0, 0.0, 0.0]];
    let positions = vec![[0.0, 0.0, 0.0]];
    let glyphs = build_tensor_glyphs(&tensors, &positions, GlyphShape::Ellipsoid, 1.0);
    assert_eq!(glyphs.len(), 1);
    // All eigenvalues should be real and match the diagonal entries approximately
    let evals = glyphs[0].eigenvalues;
    let sorted: Vec<f64> = {
        let mut v = evals.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    };
    assert!(
        (sorted[2] - 4.0).abs() < 0.1,
        "largest eigenval≈4: {:.6}",
        sorted[2]
    );
}

/// Assert `T·e_k = λ_k·e_k`, unit length and mutual orthogonality for the
/// eigenpairs returned for the Voigt tensor `[Txx, Tyy, Tzz, Txy, Txz, Tyz]`.
fn assert_eigenpairs(voigt: [f64; 6], evals: [f64; 3], evecs: [f64; 9]) {
    let [txx, tyy, tzz, txy, txz, tyz] = voigt;
    let t = [[txx, txy, txz], [txy, tyy, tyz], [txz, tyz, tzz]];
    let row = |k: usize| [evecs[3 * k], evecs[3 * k + 1], evecs[3 * k + 2]];
    for (k, &eval) in evals.iter().enumerate() {
        let e = row(k);
        let norm = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
        assert!((norm - 1.0).abs() < 1e-9, "eigenvector {k} has norm {norm}");
        for i in 0..3 {
            let te: f64 = (0..3).map(|j| t[i][j] * e[j]).sum();
            assert!(
                (te - eval * e[i]).abs() < 1e-9,
                "T·e{k} component {i}: {te} vs λ·e = {}",
                eval * e[i]
            );
        }
        for l in (k + 1)..3 {
            let f = row(l);
            let dot = e[0] * f[0] + e[1] * f[1] + e[2] * f[2];
            assert!(
                dot.abs() < 1e-9,
                "eigenvectors {k} and {l} not orthogonal: {dot}"
            );
        }
    }
}

#[test]
fn eigen3_sym_returns_sorted_eigenvalues_and_row_eigenvectors() {
    // [[1,1,0],[1,1,0],[0,0,0]]: eigenvalues 2, 0, 0 with e1 = (1,1,0)/sqrt(2).
    let voigt = [1.0, 1.0, 0.0, 1.0, 0.0, 0.0];
    let (evals, evecs) = eigen3_sym(voigt);
    assert!((evals[0] - 2.0).abs() < 1e-9, "λ1 = {}", evals[0]);
    assert!(evals[1].abs() < 1e-9 && evals[2].abs() < 1e-9, "{evals:?}");
    let s = std::f64::consts::FRAC_1_SQRT_2;
    assert!((evecs[0].abs() - s).abs() < 1e-9, "e1x = {}", evecs[0]);
    assert!((evecs[1].abs() - s).abs() < 1e-9, "e1y = {}", evecs[1]);
    assert!(evecs[2].abs() < 1e-9, "e1z = {}", evecs[2]);
    assert!(
        evecs[0] * evecs[1] > 0.0,
        "e1 must be along (1,1,0), not (1,-1,0)"
    );
    assert_eigenpairs(voigt, evals, evecs);
}

#[test]
fn eigen3_sym_sorts_a_diagonal_tensor_descending() {
    // Diagonal entries 1, 5, 3 used to come back in storage order.
    let voigt = [1.0, 5.0, 3.0, 0.0, 0.0, 0.0];
    let (evals, evecs) = eigen3_sym(voigt);
    assert!((evals[0] - 5.0).abs() < 1e-12);
    assert!((evals[1] - 3.0).abs() < 1e-12);
    assert!((evals[2] - 1.0).abs() < 1e-12);
    // Largest eigenvalue belongs to the y axis, then z, then x.
    assert!(
        (evecs[1].abs() - 1.0).abs() < 1e-12,
        "e1 = {:?}",
        &evecs[0..3]
    );
    assert!(
        (evecs[5].abs() - 1.0).abs() < 1e-12,
        "e2 = {:?}",
        &evecs[3..6]
    );
    assert!(
        (evecs[6].abs() - 1.0).abs() < 1e-12,
        "e3 = {:?}",
        &evecs[6..9]
    );
    assert_eigenpairs(voigt, evals, evecs);
}

#[test]
fn eigen3_sym_distinguishes_xz_and_yz_shear() {
    // The scientific_viz Voigt order is [Txx, Tyy, Tzz, Txy, Txz, Tyz].
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let txz = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let (evals, evecs) = eigen3_sym(txz);
    assert!((evals[0] - 1.0).abs() < 1e-9 && evals[1].abs() < 1e-9);
    assert!((evals[2] + 1.0).abs() < 1e-9, "{evals:?}");
    assert!((evecs[0].abs() - s).abs() < 1e-9 && evecs[1].abs() < 1e-9);
    assert!((evecs[2].abs() - s).abs() < 1e-9, "e1 = {:?}", &evecs[0..3]);
    assert_eigenpairs(txz, evals, evecs);

    let tyz = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let (evals, evecs) = eigen3_sym(tyz);
    assert!((evals[0] - 1.0).abs() < 1e-9);
    assert!(evecs[0].abs() < 1e-9);
    assert!((evecs[1].abs() - s).abs() < 1e-9, "e1 = {:?}", &evecs[0..3]);
    assert!((evecs[2].abs() - s).abs() < 1e-9, "e1 = {:?}", &evecs[0..3]);
    assert_eigenpairs(tyz, evals, evecs);
}

#[test]
fn eigen3_sym_general_tensor_satisfies_the_eigen_equation() {
    let voigt = [4.0, 3.0, 2.0, 1.0, 0.5, -0.7];
    let (evals, evecs) = eigen3_sym(voigt);
    assert!(evals[0] >= evals[1] && evals[1] >= evals[2], "{evals:?}");
    assert!(
        (evals.iter().sum::<f64>() - 9.0).abs() < 1e-9,
        "trace preserved"
    );
    assert_eigenpairs(voigt, evals, evecs);
}

#[test]
fn tensor_glyph_surface_points_follow_the_eigenvectors() {
    // Semi-axes 3, 1, 0.5 along world y, x and z respectively.
    let evals = [3.0, 1.0, 0.5];
    let evecs = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let centre = [1.0, 2.0, 3.0];
    let g = TensorGlyph::new(centre, evals, evecs, GlyphShape::Ellipsoid, 1.0);
    let pts = g.ellipsoid_surface_points(64, 32);
    assert_eq!(pts.len(), 64 * 32);
    let (mut y_reach, mut x_reach, mut z_reach) = (0.0_f64, 0.0_f64, 0.0_f64);
    for p in &pts {
        let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
        // Implicit ellipsoid equation in the eigenvector frame.
        let q = (d[1] / 3.0).powi(2) + (d[0] / 1.0).powi(2) + (d[2] / 0.5).powi(2);
        assert!((q - 1.0).abs() < 1e-9, "point {p:?} off the surface: {q}");
        y_reach = y_reach.max(d[1].abs());
        x_reach = x_reach.max(d[0].abs());
        z_reach = z_reach.max(d[2].abs());
    }
    assert!(
        (y_reach - 3.0).abs() < 0.02,
        "long axis must lie along y: {y_reach}"
    );
    assert!(
        (x_reach - 1.0).abs() < 0.02,
        "middle axis must lie along x: {x_reach}"
    );
    assert!(
        (z_reach - 0.5).abs() < 0.02,
        "short axis must lie along z: {z_reach}"
    );
}

#[test]
fn built_glyph_surface_is_elongated_along_the_major_eigenvector() {
    // [[2,1,0],[1,2,0],[0,0,1]]: eigenvalues 3, 1, 1 with e1 = (1,1,0)/sqrt(2).
    let glyphs = build_tensor_glyphs(
        &[[2.0, 2.0, 1.0, 1.0, 0.0, 0.0]],
        &[[0.0, 0.0, 0.0]],
        GlyphShape::Ellipsoid,
        1.0,
    );
    let g = &glyphs[0];
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let mut along_major = 0.0_f64;
    let mut along_minor = 0.0_f64;
    for p in g.ellipsoid_surface_points(72, 36) {
        along_major = along_major.max(((p[0] + p[1]) * s).abs());
        along_minor = along_minor.max(((p[0] - p[1]) * s).abs());
    }
    assert!(
        (along_major - 3.0).abs() < 0.02,
        "reach along (1,1,0): {along_major}"
    );
    assert!(
        (along_minor - 1.0).abs() < 0.02,
        "reach along (1,-1,0): {along_minor}"
    );
}

#[test]
fn tensor_glyph_linear_anisotropy_extreme() {
    let evals = [10.0, 10.0, 0.001];
    let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
    assert!(g.planar_anisotropy() > 0.5);
}

// ---- PhaseDiagramViz ----

#[test]
fn phase_diagram_clausius_clapeyron_points_count() {
    let pts =
        PhaseDiagramViz::clausius_clapeyron_line([300.0, 400.0], 50, 40_000.0, 373.0, 101325.0);
    assert_eq!(pts.len(), 50);
}

#[test]
fn phase_diagram_clausius_clapeyron_pressure_increases() {
    let pts =
        PhaseDiagramViz::clausius_clapeyron_line([300.0, 500.0], 10, 40_000.0, 373.0, 101325.0);
    // Vapor pressure should increase with temperature
    assert!(
        pts[9][1] > pts[0][1],
        "vapor pressure should increase with T"
    );
}

#[test]
fn phase_diagram_tie_line_at_correct_temperature() {
    let mut pd = PhaseDiagramViz::new();
    let temps = vec![400.0, 450.0, 500.0];
    pd.compute_tie_lines(&temps, |t| 0.3 + 0.001 * t, |t| 0.1 + 0.0005 * t);
    assert_eq!(pd.tie_lines.len(), 3);
    assert!((pd.tie_lines[1].temperature - 450.0).abs() < 1e-10);
}

#[test]
fn phase_diagram_triple_and_critical_point() {
    let mut pd = PhaseDiagramViz::new();
    pd.set_triple_point(273.16, 611.73);
    pd.set_critical_point(647.1, 22.064e6);
    assert!(pd.triple_point.is_some());
    assert!(pd.critical_point.is_some());
    let tp = pd.triple_point.unwrap();
    assert!((tp.temperature - 273.16).abs() < 1e-6);
}

#[test]
fn phase_diagram_add_region_and_line() {
    let mut pd = PhaseDiagramViz::new();
    pd.add_region(PhaseRegion {
        name: "Liquid".to_string(),
        color: [0.2, 0.5, 1.0, 0.4],
        boundary: vec![[200.0, 1e3], [300.0, 1e5]],
    });
    pd.add_coexistence_line(vec![[273.0, 1e5], [373.0, 2e5]]);
    assert_eq!(pd.regions.len(), 1);
    assert_eq!(pd.coexistence_lines.len(), 1);
}

// ---- ErrorBandPlot ----

#[test]
fn error_band_mean_within_band() {
    let x = vec![0.0, 1.0, 2.0, 3.0];
    let samples: Vec<Vec<f64>> = (0..100)
        .map(|i| {
            x.iter()
                .map(|&xi| xi * xi + (i as f64 - 50.0) * 0.01)
                .collect()
        })
        .collect();
    let band = ErrorBandPlot::from_samples(x, &samples, 2.0);
    assert!(band.mean_within_band(), "mean should lie within ±2σ band");
}

#[test]
fn error_band_upper_greater_than_lower() {
    let x = vec![0.0, 1.0, 2.0];
    let samples: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64; 3]).collect();
    let band = ErrorBandPlot::from_samples(x, &samples, 1.0);
    for (lo, hi) in band.lower.iter().zip(band.upper.iter()) {
        assert!(hi >= lo, "upper >= lower must hold");
    }
}

#[test]
fn error_band_empty_samples() {
    let x = vec![0.0, 1.0];
    let band = ErrorBandPlot::from_samples(x, &[], 1.0);
    assert_eq!(band.mean.len(), 2);
    assert!(band.mean_within_band());
}

#[test]
fn error_band_bootstrap_ci_structure() {
    let x = vec![0.0, 1.0, 2.0];
    let data: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64 * 0.1; 3]).collect();
    let band = ErrorBandPlot::bootstrap_ci(x, &data, 0.95, 50, 1);
    assert_eq!(band.x.len(), 3);
    assert!(band.mean_within_band());
}

/// `rows × cols` observations `j + U(-1, 1)` drawn from a seeded generator.
fn noisy_rows(seed: u64, rows: usize, cols: usize) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..rows)
        .map(|_| {
            (0..cols)
                .map(|j| j as f64 + rng.random_range(-1.0..1.0))
                .collect()
        })
        .collect()
}

#[test]
fn error_band_bootstrap_ci_is_a_real_interval() {
    // Rotating the rows (the old "resampling") always reproduced the
    // sample mean, giving lower == upper == mean.  A genuine bootstrap of
    // the mean of U(-1, 1) noise has half-width ≈ 1.96·σ/√m.
    let rows = 400;
    let data = noisy_rows(1, rows, 3);
    let band = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0, 2.0], &data, 0.95, 2000, 42);
    let expected = 1.96 * (1.0_f64 / 3.0).sqrt() / (rows as f64).sqrt();
    for j in 0..3 {
        assert!(
            band.lower[j] < band.mean[j] && band.mean[j] < band.upper[j],
            "column {j}: [{}, {}] must strictly contain the mean {}",
            band.lower[j],
            band.upper[j],
            band.mean[j]
        );
        let half_width = 0.5 * (band.upper[j] - band.lower[j]);
        assert!(
            half_width > 0.7 * expected && half_width < 1.3 * expected,
            "column {j}: half-width {half_width}, expected about {expected}"
        );
    }
    assert!(band.mean_within_band());
}

#[test]
fn error_band_bootstrap_ci_is_deterministic_per_seed() {
    let data = noisy_rows(5, 60, 4);
    let x = vec![0.0, 1.0, 2.0, 3.0];
    let a = ErrorBandPlot::bootstrap_ci(x.clone(), &data, 0.9, 300, 7);
    let b = ErrorBandPlot::bootstrap_ci(x.clone(), &data, 0.9, 300, 7);
    let c = ErrorBandPlot::bootstrap_ci(x, &data, 0.9, 300, 8);
    assert_eq!(a.lower, b.lower);
    assert_eq!(a.upper, b.upper);
    assert_ne!(
        a.lower, c.lower,
        "a different seed must resample differently"
    );
}

#[test]
fn error_band_bootstrap_ci_zero_replicates_collapse_to_the_mean() {
    let data = noisy_rows(3, 10, 3);
    let band = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0, 2.0], &data, 0.95, 0, 0);
    assert_eq!(band.lower, band.mean);
    assert_eq!(band.upper, band.mean);
}

#[test]
fn error_band_bootstrap_ci_wider_level_gives_a_wider_band() {
    let data = noisy_rows(9, 100, 2);
    let x = vec![0.0, 1.0];
    let narrow = ErrorBandPlot::bootstrap_ci(x.clone(), &data, 0.5, 1000, 3);
    let wide = ErrorBandPlot::bootstrap_ci(x, &data, 0.99, 1000, 3);
    for j in 0..2 {
        assert!(
            wide.upper[j] - wide.lower[j] > narrow.upper[j] - narrow.lower[j],
            "column {j}: 99% band must be wider than the 50% band"
        );
    }
}

#[test]
fn error_band_bootstrap_ci_constant_data_has_zero_width_only_without_variance() {
    let data: Vec<Vec<f64>> = (0..30).map(|_| vec![2.5, -1.0]).collect();
    let band = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0], &data, 0.95, 100, 11);
    assert_eq!(band.mean, vec![2.5, -1.0]);
    assert_eq!(band.lower, band.mean);
    assert_eq!(band.upper, band.mean);
}

#[test]
fn error_band_bootstrap_ci_ragged_rows_and_odd_levels_do_not_panic() {
    let data = vec![vec![1.0, 2.0, 3.0], vec![2.0, 3.0], vec![3.0]];
    let band = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0, 2.0], &data, 0.95, 200, 5);
    assert!((band.mean[0] - 2.0).abs() < 1e-12);
    assert!((band.mean[1] - 2.5).abs() < 1e-12);
    assert!(
        (band.mean[2] - 3.0).abs() < 1e-12,
        "only one row has column 2"
    );
    assert_eq!(band.lower[2], 3.0);
    assert_eq!(band.upper[2], 3.0);
    for j in 0..3 {
        assert!(band.lower[j] <= band.upper[j]);
    }
    // Out-of-range and NaN confidence levels, and more columns than data.
    let clamped = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0], &data, 5.0, 50, 1);
    assert!(clamped.lower[0] <= clamped.upper[0]);
    let nan_level = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0], &data, f64::NAN, 50, 1);
    assert_eq!(nan_level.lower, nan_level.mean);
    let wide = ErrorBandPlot::bootstrap_ci(vec![0.0; 5], &data, 0.95, 50, 1);
    assert_eq!(wide.mean.len(), 5);
    assert_eq!(wide.mean[4], 0.0);
}

#[test]
fn error_band_bootstrap_ci_empty_inputs() {
    let band = ErrorBandPlot::bootstrap_ci(vec![0.0, 1.0], &[], 0.95, 100, 0);
    assert_eq!(band.mean, vec![0.0, 0.0]);
    assert!(band.mean_within_band());
    let band = ErrorBandPlot::bootstrap_ci(Vec::new(), &[vec![1.0]], 0.95, 100, 0);
    assert!(band.x.is_empty() && band.mean.is_empty());
}

#[test]
fn quantile_sorted_interpolates_linearly() {
    let sorted = [1.0, 2.0, 3.0, 4.0];
    assert_eq!(quantile_sorted(&sorted, 0.0), Some(1.0));
    assert_eq!(quantile_sorted(&sorted, 1.0), Some(4.0));
    assert_eq!(quantile_sorted(&sorted, 0.5), Some(2.5));
    assert!((quantile_sorted(&sorted, 1.0 / 3.0).unwrap_or(f64::NAN) - 2.0).abs() < 1e-12);
    assert_eq!(quantile_sorted(&[7.0], 0.3), Some(7.0));
    assert_eq!(quantile_sorted(&[], 0.5), None);
    // Out-of-range probabilities are clamped.
    assert_eq!(quantile_sorted(&sorted, -1.0), Some(1.0));
    assert_eq!(quantile_sorted(&sorted, 2.0), Some(4.0));
}

// ---- BifurcationDiagram ----

#[test]
fn bifurcation_logistic_period1_for_small_r() {
    // r < 3: logistic map has a single stable fixed point
    let diag = BifurcationDiagram::logistic_map([2.5, 2.9], 5, 500, 100);
    let stable = diag.stable_points();
    assert!(
        !stable.is_empty(),
        "should have stable points for r ∈ [2.5, 2.9]"
    );
}

#[test]
fn bifurcation_logistic_chaotic_for_large_r() {
    // r ≈ 4.0: logistic map is fully chaotic
    let diag = BifurcationDiagram::logistic_map([3.9, 4.0], 5, 500, 200);
    let chaotic = diag.chaotic_points();
    assert!(!chaotic.is_empty(), "should have chaotic points for r ≈ 4");
}

#[test]
fn bifurcation_diagram_points_in_unit_interval() {
    let diag = BifurcationDiagram::logistic_map([2.5, 4.0], 10, 100, 50);
    for p in &diag.points {
        assert!(p.x >= 0.0 && p.x <= 1.0, "x={:.6} not in [0,1]", p.x);
    }
}

// ---- FractalViz ----

#[test]
fn mandelbrot_origin_does_not_escape() {
    let f = FractalViz::mandelbrot(10, 10, 256);
    let v = f.escape_value(0.0, 0.0);
    // Origin z=0 maps to c=0, stays bounded → max_iter
    assert!(
        (v - 256.0).abs() < 10.0,
        "origin escape value={v:.6} should be near max_iter"
    );
}

#[test]
fn mandelbrot_far_point_escapes_fast() {
    let f = FractalViz::mandelbrot(10, 10, 256);
    let v = f.escape_value(10.0, 10.0);
    assert!(v < 10.0, "far-away point should escape quickly, got {v:.6}");
}

#[test]
fn julia_set_render_dimensions() {
    let f = FractalViz::julia(8, 6, 64, -0.7, 0.27);
    let buf = f.render(-1.5, 1.5, -1.0, 1.0);
    assert_eq!(buf.len(), 48);
}

#[test]
fn mandelbrot_render_dimensions() {
    let f = FractalViz::mandelbrot(4, 4, 32);
    let buf = f.render(-2.0, 1.0, -1.5, 1.5);
    assert_eq!(buf.len(), 16);
}

#[test]
fn fractal_orbit_trap_positive() {
    let f = FractalViz {
        width: 4,
        height: 4,
        max_iter: 64,
        escape_radius: 2.0,
        coloring: FractalColoring::OrbitTrap,
        julia_c: None,
    };
    let v = f.orbit_trap_value(0.0, 0.0);
    assert!(v >= 0.0, "orbit trap should be non-negative: {v:.6}");
}

#[test]
fn fractal_escape_time_coloring() {
    let mut f = FractalViz::mandelbrot(4, 4, 100);
    f.coloring = FractalColoring::EscapeTime;
    let v = f.escape_value(2.5, 0.0);
    assert!(v < 100.0);
}

// ---- NetworkViz ----

#[test]
fn network_layout_nodes_on_canvas() {
    let mut net = NetworkViz::new(100.0, 100.0);
    for i in 0..5 {
        net.add_node(i, 50.0, 50.0);
    }
    net.add_edge(0, 1, 1.0);
    net.add_edge(1, 2, 1.0);
    let k = (100.0 * 100.0 / 5.0_f64).sqrt();
    net.layout_fruchterman_reingold(50, k);
    for node in &net.nodes {
        assert!(node.position[0] >= 0.0 && node.position[0] <= 100.0);
        assert!(node.position[1] >= 0.0 && node.position[1] <= 100.0);
    }
}

#[test]
fn network_degree_updated_on_add_edge() {
    let mut net = NetworkViz::new(100.0, 100.0);
    net.add_node(0, 10.0, 10.0);
    net.add_node(1, 90.0, 90.0);
    net.add_edge(0, 1, 1.0);
    assert_eq!(net.nodes[0].degree, 1);
    assert_eq!(net.nodes[1].degree, 1);
}

#[test]
fn network_color_by_degree_changes_colors() {
    let mut net = NetworkViz::new(100.0, 100.0);
    for i in 0..3 {
        net.add_node(i, 50.0, 50.0);
    }
    net.add_edge(0, 1, 1.0);
    net.add_edge(0, 2, 1.0);
    let orig_color = net.nodes[0].color;
    net.color_by_degree();
    // Node 0 has degree 2 vs nodes 1,2 with degree 1 → different color
    assert!(
        (net.nodes[0].color[0] - orig_color[0]).abs() > 0.01
            || (net.nodes[1].color[0] - orig_color[0]).abs() > 0.01
    );
}

// ---- PhasePortrait ----

#[test]
fn phase_portrait_fixed_point_at_zero_velocity() {
    // f(x,y) = (-x, -y): origin is a stable node
    // Integrate long enough with a small step to reach dist < 0.1
    let f = |x: f64, y: f64| (-x, -y);
    let seeds = [[-0.5_f64, 0.5], [0.5, -0.5]];
    let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 5, &seeds, 0.05, 200);
    // Trajectories should approach origin (e^{-10} ≈ 4.5e-5 << 0.1)
    for traj in &pp.trajectories {
        let last = traj.states.last().unwrap();
        let dist = (last[0] * last[0] + last[1] * last[1]).sqrt();
        assert!(
            dist < 0.1,
            "trajectory should converge to origin: dist={dist:.6}"
        );
    }
}

#[test]
fn phase_portrait_vector_field_count() {
    let f = |x: f64, y: f64| (y, -x);
    let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 8, &[], 0.01, 10);
    assert_eq!(pp.vector_field.len(), 64);
}

#[test]
fn phase_portrait_trajectory_length() {
    let f = |_x: f64, _y: f64| (0.0, 0.0);
    let seeds = [[0.5_f64, 0.5]];
    let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 4, &seeds, 0.1, 20);
    assert_eq!(pp.trajectories[0].states.len(), 21);
}

#[test]
fn phase_portrait_classify_stable_node() {
    let f = |x: f64, y: f64| (-2.0 * x, -3.0 * y);
    let mut pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 4, &[], 0.01, 10);
    pp.classify_fixed_points(&[[0.0, 0.0]], &f, 1e-5);
    assert_eq!(pp.fixed_points.len(), 1);
    let fp = pp.fixed_points[0];
    assert!(
        fp.stability == FixedPointType::StableNode || fp.stability == FixedPointType::StableFocus,
        "origin of (-2x,-3y) should be stable: {:?}",
        fp.stability
    );
}

fn classify_at(f: &impl Fn(f64, f64) -> (f64, f64), point: [f64; 2]) -> Option<FixedPointType> {
    let mut pp = PhasePortrait::compute(f, [-1.0, 1.0], [-1.0, 1.0], 2, &[], 0.01, 1);
    pp.classify_fixed_points(&[point], f, 1e-5);
    pp.fixed_points.first().map(|fp| fp.stability)
}

#[test]
fn phase_portrait_classifies_a_linear_saddle() {
    // Jacobian [[1, 3], [1, 2]]: trace 3, det -1.  Halving ∂f/∂y gave
    // det +0.5 and an "unstable node".
    let f = |x: f64, y: f64| (x + 3.0 * y, x + 2.0 * y);
    assert_eq!(classify_at(&f, [0.0, 0.0]), Some(FixedPointType::Saddle));
    // The Jacobian of a linear field is the same everywhere.
    assert_eq!(classify_at(&f, [0.3, -0.2]), Some(FixedPointType::Saddle));
}

#[test]
fn phase_portrait_classifies_a_stable_focus() {
    // Jacobian [[0, 1], [-1, -1.5]]: trace -1.5, det 1, discriminant < 0.
    let f = |x: f64, y: f64| (y, -x - 1.5 * y);
    assert_eq!(
        classify_at(&f, [0.0, 0.0]),
        Some(FixedPointType::StableFocus)
    );
}

#[test]
fn phase_portrait_classifies_every_fixed_point_type() {
    // Unstable node: Jacobian diag(2, 3).
    let node = |x: f64, y: f64| (2.0 * x, 3.0 * y);
    assert_eq!(
        classify_at(&node, [0.0, 0.0]),
        Some(FixedPointType::UnstableNode)
    );
    // Stable node: diag(-2, -3).
    let stable = |x: f64, y: f64| (-2.0 * x, -3.0 * y);
    assert_eq!(
        classify_at(&stable, [0.0, 0.0]),
        Some(FixedPointType::StableNode)
    );
    // Unstable focus: [[0.5, -1], [1, 0.5]].
    let focus = |x: f64, y: f64| (0.5 * x - y, x + 0.5 * y);
    assert_eq!(
        classify_at(&focus, [0.0, 0.0]),
        Some(FixedPointType::UnstableFocus)
    );
    // Centre: [[0, 1], [-1, 0]].
    let centre = |x: f64, y: f64| (y, -x);
    assert_eq!(
        classify_at(&centre, [0.0, 0.0]),
        Some(FixedPointType::Center)
    );
}

#[test]
fn phase_portrait_centre_test_is_relative_to_the_eigenvalue_scale() {
    // Weakly damped oscillators at very small and very large scale: the
    // trace is 0.1 of sqrt(det) in both, so both are foci, never centres.
    for scale in [1e-9, 1.0, 1e6] {
        let f = move |x: f64, y: f64| (scale * y, scale * (-x - 0.1 * y));
        assert_eq!(
            classify_at(&f, [0.0, 0.0]),
            Some(FixedPointType::StableFocus),
            "scale {scale}"
        );
    }
    // An undamped oscillator is a centre at any scale.
    for scale in [1e-9, 1.0, 1e6] {
        let f = move |x: f64, y: f64| (scale * y, -scale * x);
        assert_eq!(
            classify_at(&f, [0.0, 0.0]),
            Some(FixedPointType::Center),
            "scale {scale}"
        );
    }
}

#[test]
fn phase_portrait_classification_of_nonlinear_field_uses_the_local_jacobian() {
    // Damped pendulum x' = v, v' = -sin(x) - 0.3 v: the downward position
    // (x = 0) is a stable focus, the inverted one (x = pi) a saddle.
    let f = |x: f64, v: f64| (v, -x.sin() - 0.3 * v);
    assert_eq!(
        classify_at(&f, [0.0, 0.0]),
        Some(FixedPointType::StableFocus)
    );
    assert_eq!(classify_at(&f, [PI, 0.0]), Some(FixedPointType::Saddle));
}

#[test]
fn phase_portrait_skips_unclassifiable_candidates_and_defaults_the_step() {
    let f = |x: f64, y: f64| if x > 0.5 { (f64::NAN, 0.0) } else { (-x, -y) };
    let mut pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 2, &[], 0.01, 1);
    // The stencil of the second candidate straddles the NaN region.
    pp.classify_fixed_points(&[[0.0, 0.0], [0.5, 0.0]], &f, 1e-5);
    assert_eq!(pp.fixed_points.len(), 1);
    assert_eq!(pp.fixed_points[0].location, [0.0, 0.0]);
    // eps = 0 (or non-finite) falls back to a sensible default step.
    let g = |x: f64, y: f64| (-2.0 * x, -3.0 * y);
    let mut pp = PhasePortrait::compute(&g, [-1.0, 1.0], [-1.0, 1.0], 2, &[], 0.01, 1);
    pp.classify_fixed_points(&[[0.0, 0.0]], &g, 0.0);
    pp.classify_fixed_points(&[[0.0, 0.0]], &g, f64::NAN);
    assert_eq!(pp.fixed_points.len(), 2);
    for fp in &pp.fixed_points {
        assert_eq!(fp.stability, FixedPointType::StableNode);
    }
}

// ---- Hodograph ----

#[test]
fn hodograph_projectile_linear_vy() {
    let h = Hodograph::projectile(10.0, 20.0, 9.81, 0.1, 10);
    // vx stays constant
    for v in &h.velocities {
        assert!((v[0] - 10.0).abs() < 1e-10, "vx should be constant");
    }
    // vy decreases linearly
    assert!(h.velocities.last().unwrap()[1] < h.velocities[0][1]);
}

#[test]
fn hodograph_arc_length_positive() {
    let h = Hodograph::projectile(5.0, 10.0, 9.81, 0.1, 20);
    assert!(h.arc_length() > 0.0);
}

#[test]
fn hodograph_mean_speed_correct() {
    let mut h = Hodograph::new();
    h.push(0.0, 3.0, 4.0); // speed = 5
    h.push(1.0, 0.0, 5.0); // speed = 5
    let ms = h.mean_speed();
    assert!((ms - 5.0).abs() < 1e-10, "mean speed={ms:.6}");
}

// ---- LissajousFigures ----

#[test]
fn lissajous_point_count() {
    let fig = LissajousFigure::new(1.0, 2.0, 0.0, 1.0, 1.0, [0.0, 2.0 * PI], 200);
    assert_eq!(fig.points.len(), 200);
}

#[test]
fn lissajous_amplitude_bounds() {
    let amp_x = 3.0;
    let amp_y = 2.0;
    let fig = LissajousFigure::new(1.0, 2.0, PI / 4.0, amp_x, amp_y, [0.0, 4.0 * PI], 500);
    for &[x, y] in &fig.points {
        assert!(
            x.abs() <= amp_x + 1e-10,
            "x={x:.6} exceeds amplitude {amp_x}"
        );
        assert!(
            y.abs() <= amp_y + 1e-10,
            "y={y:.6} exceeds amplitude {amp_y}"
        );
    }
}

#[test]
fn lissajous_integer_ratio_is_closed() {
    // 1:2 ratio → closed figure
    let fig = LissajousFigure::new(1.0, 2.0, 0.0, 1.0, 1.0, [0.0, 2.0 * PI], 100);
    assert!(
        fig.is_closed(1e-6),
        "1:2 frequency ratio should be a closed curve"
    );
}

#[test]
fn lissajous_bounding_box_within_amplitude() {
    let fig = LissajousFigure::new(3.0, 2.0, 0.5, 2.0, 1.5, [0.0, 10.0 * PI], 1000);
    let (xb, yb) = fig.bounding_box();
    assert!(xb[0] >= -2.0 - 1e-10 && xb[1] <= 2.0 + 1e-10);
    assert!(yb[0] >= -1.5 - 1e-10 && yb[1] <= 1.5 + 1e-10);
}

#[test]
fn lissajous_circle_at_equal_freqs_and_phase() {
    // freq_x = freq_y, phase = π/2 → circle
    let fig = LissajousFigure::new(1.0, 1.0, PI / 2.0, 1.0, 1.0, [0.0, 2.0 * PI], 360);
    // All points should lie on unit circle
    for &[x, y] in &fig.points {
        let r = (x * x + y * y).sqrt();
        assert!((r - 1.0).abs() < 0.01, "r={r:.6} should be 1 for a circle");
    }
}
