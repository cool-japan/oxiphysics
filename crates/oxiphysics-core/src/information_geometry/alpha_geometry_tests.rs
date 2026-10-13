// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use super::*;

const PARAMS: [f64; 2] = [0.3, 1.2];

#[test]
fn test_alpha_zero_is_levi_civita() {
    let manifold = StatisticalManifold::gaussian_1d();
    let ag = AlphaGeometry::new(0.0, manifold.clone());
    let a = ag.alpha_connection(&PARAMS).expect("alpha connection");
    let lc = manifold.christoffel_symbols(&PARAMS).expect("christoffel");
    for k in 0..2 {
        for i in 0..2 {
            for j in 0..2 {
                assert!((a[k][i][j] - lc[k][i][j]).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn test_dual_connections_average_to_levi_civita() {
    // (Gamma^(alpha) + Gamma^(-alpha)) / 2 = Gamma^(0), and they differ by alpha g^-1 T.
    let manifold = StatisticalManifold::gaussian_1d();
    let ag = AlphaGeometry::new(0.6, manifold.clone());
    let a = ag.alpha_connection(&PARAMS).expect("alpha");
    let d = ag.dual_connection(&PARAMS).expect("dual");
    let lc = manifold.christoffel_symbols(&PARAMS).expect("lc");
    let sigma = PARAMS[1];
    for k in 0..2 {
        for i in 0..2 {
            for j in 0..2 {
                assert!((0.5 * (a[k][i][j] + d[k][i][j]) - lc[k][i][j]).abs() < 1e-7);
            }
        }
    }
    // Gamma^(-a) - Gamma^(a) = a g^{-1} T: e.g. the sigma-component of the
    // (sigma, sigma) entry is a * (sigma^2 / 2) * 8 / sigma^3 = 4 a / sigma.
    let diff = d[1][1][1] - a[1][1][1];
    assert!((diff - 4.0 * 0.6 / sigma).abs() < 1e-6, "diff = {diff}");
    // And the mu-component of (mu, sigma) is a * sigma^2 * 2 / sigma^3 = 2 a / sigma.
    let diff = d[0][0][1] - a[0][0][1];
    assert!((diff - 2.0 * 0.6 / sigma).abs() < 1e-6, "diff = {diff}");
}

#[test]
fn test_levi_civita_curvature_is_minus_one_half() {
    // The Gaussian manifold with the Fisher metric has constant curvature -1/2.
    let ag = AlphaGeometry::new(0.0, StatisticalManifold::gaussian_1d());
    let k = ag.sectional_curvature(&PARAMS, 0, 1).expect("curvature");
    assert!((k + 0.5).abs() < 1e-4, "K = {k}");
    let r = ag.curvature_tensor(&PARAMS).expect("tensor");
    // R^mu_{sigma mu sigma} = -1/sigma^2.
    let expected = -1.0 / (PARAMS[1] * PARAMS[1]);
    assert!(
        (r[0][1][0][1] - expected).abs() < 1e-4,
        "R = {}",
        r[0][1][0][1]
    );
    // Antisymmetry in the last two indices.
    assert!((r[0][1][0][1] + r[0][1][1][0]).abs() < 1e-4);
}

#[test]
fn test_exponential_family_is_dually_flat() {
    // Exponential families are flat for alpha = +1 and alpha = -1, in any
    // coordinates, so the whole curvature tensor vanishes.
    for alpha in [1.0, -1.0] {
        let ag = AlphaGeometry::new(alpha, StatisticalManifold::gaussian_1d());
        let r = ag.curvature_tensor(&PARAMS).expect("tensor");
        let max_abs = r
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(max_abs < 1e-4, "alpha = {alpha}: max |R| = {max_abs}");
    }
}

#[test]
fn test_alpha_curvature_scales_with_one_minus_alpha_squared() {
    // For the Gaussian family the alpha-curvature is -(1 - alpha^2) / 2.
    let ag = AlphaGeometry::new(0.5, StatisticalManifold::gaussian_1d());
    let k = ag.sectional_curvature(&PARAMS, 0, 1).expect("curvature");
    assert!((k + 0.375).abs() < 1e-4, "K = {k}");
    assert!(ag.sectional_curvature(&PARAMS, 0, 0).is_err());
    assert_eq!(ag.manifold().dim(), 2);
}
