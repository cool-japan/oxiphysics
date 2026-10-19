// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use super::*;

const MU: f64 = 0.3;
const SIGMA: f64 = 1.2;

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: expected {expected}, got {actual}"
    );
}

/// Discrete exponential family `p(x; θ) ∝ exp(θ x)` on `x ∈ {0, 1, 2, 3}`.
fn discrete_exponential_family() -> StatisticalManifold {
    let nodes = (0..4).map(|j| (vec![j as f64], 1.0)).collect();
    StatisticalManifold::new(
        1,
        |x: &[f64], theta: &[f64]| {
            let log_partition = (0..4)
                .map(|j| (theta[0] * j as f64).exp())
                .sum::<f64>()
                .ln();
            theta[0] * x[0] - log_partition
        },
        nodes,
    )
    .expect("valid manifold")
}

#[test]
fn test_gaussian_fisher_metric_matches_closed_form() {
    let m = StatisticalManifold::gaussian_1d();
    assert_eq!(m.dim(), 2);
    let g = m.fisher_metric(&[MU, SIGMA]).expect("metric");
    let s2 = SIGMA * SIGMA;
    assert_close(g[0][0], 1.0 / s2, 1e-5, "g_mu_mu");
    assert_close(g[1][1], 2.0 / s2, 1e-5, "g_sigma_sigma");
    assert_close(g[0][1], 0.0, 1e-5, "g_mu_sigma");
    assert_close(g[1][0], 0.0, 1e-5, "g_sigma_mu");
}

#[test]
fn test_gaussian_fisher_metric_at_standard_normal() {
    // The removed proxy returned g_sigma_sigma = -1 here; the true value is 2.
    let m = StatisticalManifold::gaussian_1d();
    let g = m.fisher_metric(&[0.0, 1.0]).expect("metric");
    assert_close(g[0][0], 1.0, 1e-6, "g_mu_mu");
    assert_close(g[1][1], 2.0, 1e-6, "g_sigma_sigma");
}

#[test]
fn test_gaussian_skewness_tensor_matches_closed_form() {
    let m = StatisticalManifold::gaussian_1d();
    let t = m.skewness_tensor(&[MU, SIGMA]).expect("skewness");
    let s3 = SIGMA.powi(3);
    // T_{mu mu sigma} = E[z^2 (z^2 - 1)] / sigma^3 = 2 / sigma^3,
    // T_{sigma sigma sigma} = E[(z^2 - 1)^3] / sigma^3 = 8 / sigma^3.
    assert_close(t[0][0][1], 2.0 / s3, 1e-6, "T_mmS");
    assert_close(t[0][1][0], 2.0 / s3, 1e-6, "T_mSm");
    assert_close(t[1][0][0], 2.0 / s3, 1e-6, "T_Smm");
    assert_close(t[1][1][1], 8.0 / s3, 1e-6, "T_SSS");
    assert_close(t[0][0][0], 0.0, 1e-6, "T_mmm");
    assert_close(t[0][1][1], 0.0, 1e-6, "T_mSS");
}

#[test]
fn test_gaussian_christoffel_symbols_match_closed_form() {
    let m = StatisticalManifold::gaussian_1d();
    let gamma = m.christoffel_symbols(&[MU, SIGMA]).expect("christoffel");
    // Gamma^mu_{mu sigma} = -1/sigma, Gamma^sigma_{mu mu} = 1/(2 sigma),
    // Gamma^sigma_{sigma sigma} = -1/sigma, all others vanish.
    assert_close(gamma[0][0][1], -1.0 / SIGMA, 1e-6, "G^m_mS");
    assert_close(gamma[0][1][0], -1.0 / SIGMA, 1e-6, "G^m_Sm");
    assert_close(gamma[1][0][0], 0.5 / SIGMA, 1e-6, "G^S_mm");
    assert_close(gamma[1][1][1], -1.0 / SIGMA, 1e-6, "G^S_SS");
    assert_close(gamma[0][0][0], 0.0, 1e-6, "G^m_mm");
    assert_close(gamma[0][1][1], 0.0, 1e-6, "G^m_SS");
    assert_close(gamma[1][0][1], 0.0, 1e-6, "G^S_mS");
}

#[test]
fn test_alpha_connections_satisfy_amari_relation() {
    // Gamma^(alpha) = Gamma^(0) - (alpha / 2) T (first kind).
    let m = StatisticalManifold::gaussian_1d();
    let params = [MU, SIGMA];
    let g0 = m.alpha_christoffel_first_kind(&params, 0.0).expect("lc");
    let ga = m.alpha_christoffel_first_kind(&params, 0.7).expect("alpha");
    let t = m.skewness_tensor(&params).expect("skewness");
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                assert_close(
                    ga[i][j][k],
                    g0[i][j][k] - 0.35 * t[i][j][k],
                    1e-7,
                    "Gamma^(0.7)",
                );
            }
        }
    }
}

#[test]
fn test_e_connection_vanishes_in_natural_parameters() {
    // For an exponential family in natural coordinates the e-connection
    // (alpha = 1) vanishes and the m-connection equals the skewness tensor.
    let m = discrete_exponential_family();
    let theta = [0.4];
    let e = m.alpha_christoffel_first_kind(&theta, 1.0).expect("e");
    assert_close(e[0][0][0], 0.0, 1e-7, "e-connection");
    let mix = m.alpha_christoffel_first_kind(&theta, -1.0).expect("m");
    let t = m.skewness_tensor(&theta).expect("skewness");
    assert_close(mix[0][0][0], t[0][0][0], 1e-7, "m-connection");
    // The metric equals Var_theta[x].
    let weights: Vec<f64> = (0..4).map(|j| (theta[0] * j as f64).exp()).collect();
    let z: f64 = weights.iter().sum();
    let mean: f64 = weights
        .iter()
        .enumerate()
        .map(|(j, w)| j as f64 * w)
        .sum::<f64>()
        / z;
    let var: f64 = weights
        .iter()
        .enumerate()
        .map(|(j, w)| (j as f64 - mean).powi(2) * w)
        .sum::<f64>()
        / z;
    let g = m.fisher_metric(&theta).expect("metric");
    assert_close(g[0][0], var, 1e-9, "Var_theta[x]");
}

#[test]
fn test_exponential_map_matches_closed_form_geodesic() {
    // Along the sigma axis the Fisher-Rao geodesic is sigma(t) = sigma0 exp(t v / sigma0).
    let m = StatisticalManifold::gaussian_1d();
    let end = m
        .exponential_map(&[0.0, 1.0], &[0.0, 0.5], 1.0)
        .expect("exp map");
    assert_close(end[0], 0.0, 1e-9, "mu stays fixed");
    assert_close(end[1], 0.5f64.exp(), 1e-6, "sigma(1)");
}

#[test]
fn test_geodesic_hits_endpoints_and_distance_matches_closed_form() {
    let m = StatisticalManifold::gaussian_1d();
    let p = [0.0, 1.0];
    let q = [1.0, 1.5];
    // Shooting once, then sampling the geodesic through the exponential map.
    let v = m.log_map(&p, &q).expect("log map");
    let end = m.exponential_map(&p, &v, 1.0).expect("gamma(1)");
    for d in 0..2 {
        assert_close(end[d], q[d], 1e-6, "gamma(1)");
    }
    let expected = crate::information_theory::fisher_rao_gaussian(p[0], p[1], q[0], q[1]);
    let g = m.fisher_metric(&p).expect("metric");
    let speed2: f64 = (0..2)
        .map(|i| (0..2).map(|j| g[i][j] * v[i] * v[j]).sum::<f64>())
        .sum();
    assert_close(speed2.sqrt(), expected, 1e-5, "Fisher-Rao distance");
    // Geodesic midpoint is equidistant from both ends.
    let mid = m.geodesic(&p, &q, 0.5).expect("midpoint");
    let d1 = crate::information_theory::fisher_rao_gaussian(p[0], p[1], mid[0], mid[1]);
    let d2 = crate::information_theory::fisher_rao_gaussian(mid[0], mid[1], q[0], q[1]);
    assert_close(d1, 0.5 * expected, 1e-5, "first half");
    assert_close(d2, 0.5 * expected, 1e-5, "second half");
    // The geodesic bends toward larger sigma (it is not the straight chord).
    assert!(mid[1] > 0.5 * (p[1] + q[1]) + 1e-3, "midpoint {mid:?}");
}

#[test]
fn test_geodesic_distance_along_sigma_axis() {
    // d(N(0, 1), N(0, 2)) = sqrt(2) ln 2.  The shooting iterates reach
    // sigma ~ e, beyond the default [-12, 12] grid, so use a wider one.
    let m = StatisticalManifold::gaussian_1d_on_grid(-30.0, 30.0, 601).expect("grid");
    let d = m
        .geodesic_distance(&[0.0, 1.0], &[0.0, 2.0])
        .expect("distance");
    assert_close(d, std::f64::consts::SQRT_2 * 2.0f64.ln(), 1e-5, "distance");
}

#[test]
fn test_invalid_inputs_are_rejected() {
    let m = StatisticalManifold::gaussian_1d();
    assert!(matches!(
        m.fisher_metric(&[0.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    assert!(matches!(
        m.fisher_metric(&[0.0, -1.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    assert!(m.fisher_metric(&[0.0, f64::NAN]).is_err());
    // A third, unused parameter makes the metric singular.
    let redundant = StatisticalManifold::new(
        3,
        |x: &[f64], theta: &[f64]| gaussian_log_density(x, &theta[..2]),
        trapezoid_nodes(-12.0, 12.0, 481),
    )
    .expect("valid manifold");
    assert!(matches!(
        redundant.christoffel_symbols(&[0.0, 1.0, 0.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    // Unnormalised density is detected via the quadrature mass.
    let unnormalised = StatisticalManifold::new(
        2,
        |x: &[f64], theta: &[f64]| gaussian_log_density(x, theta) + 1.0,
        trapezoid_nodes(-12.0, 12.0, 481),
    )
    .expect("valid manifold");
    assert!(unnormalised.fisher_metric(&[0.0, 1.0]).is_err());
}

#[test]
fn test_constructor_validation() {
    let ld = |_: &[f64], _: &[f64]| 0.0;
    assert!(StatisticalManifold::new(0, ld, vec![(vec![0.0], 1.0)]).is_err());
    assert!(StatisticalManifold::new(1, ld, vec![]).is_err());
    assert!(StatisticalManifold::new(1, ld, vec![(vec![0.0], -1.0)]).is_err());
    assert!(StatisticalManifold::new(1, ld, vec![(vec![0.0], 0.0)]).is_err());
    assert!(StatisticalManifold::new(1, ld, vec![(vec![f64::NAN], 1.0)]).is_err());
    assert!(StatisticalManifold::gaussian_1d_on_grid(1.0, -1.0, 100).is_err());
    assert!(StatisticalManifold::gaussian_1d_on_grid(-1.0, 1.0, 2).is_err());
    let wide = StatisticalManifold::gaussian_1d_on_grid(-60.0, 60.0, 2401).expect("grid");
    let g = wide.fisher_metric(&[5.0, 4.0]).expect("metric");
    assert_close(g[1][1], 2.0 / 16.0, 1e-6, "g_sigma_sigma on wide grid");
    let dbg = format!("{wide:?}");
    assert!(dbg.contains("StatisticalManifold") && dbg.contains("2401"));
}

#[test]
fn test_gaussian_fisher_metric_at_sigma_two() {
    // N(0, 2^2): g = diag(1 / sigma^2, 2 / sigma^2) = diag(1/4, 1/2).
    let m = StatisticalManifold::gaussian_1d();
    let g = m.fisher_metric(&[0.0, 2.0]).expect("metric");
    assert_close(g[0][0], 0.25, 1e-6, "g_mu_mu");
    assert_close(g[1][1], 0.5, 1e-6, "g_sigma_sigma");
    assert_close(g[0][1], 0.0, 1e-6, "g_mu_sigma");
}

#[test]
fn test_poisson_fisher_metric_equals_rate() {
    // Poisson family in natural coordinates: psi(theta) = e^theta, so the Fisher
    // metric is psi''(theta) = e^theta.
    let nodes: Vec<(Vec<f64>, f64)> = (0..=60).map(|j| (vec![j as f64], 1.0)).collect();
    let m = StatisticalManifold::new(
        1,
        |x: &[f64], theta: &[f64]| {
            let count = x[0];
            let log_factorial: f64 = (2..=(count as u32)).map(|i| f64::from(i).ln()).sum();
            theta[0] * count - theta[0].exp() - log_factorial
        },
        nodes,
    )
    .expect("valid manifold");
    for theta in [-0.5, 0.0, 0.7, 1.3] {
        let g = m.fisher_metric(&[theta]).expect("metric");
        assert_close(g[0][0], theta.exp(), 1e-6, "g = e^theta");
    }
}
