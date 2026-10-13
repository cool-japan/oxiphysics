// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use super::*;

fn stat_x(x: &[f64]) -> f64 {
    x.first().copied().unwrap_or(0.0)
}

fn stat_x2(x: &[f64]) -> f64 {
    let v = x.first().copied().unwrap_or(0.0);
    v * v
}

fn stat_const(_: &[f64]) -> f64 {
    1.0
}

/// Unnormalised `p_j ∝ exp(Σ_a θ_a T_a(j))`.
fn family_member(theta: &[f64], stats: &[fn(&[f64]) -> f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|j| {
            let x = [j as f64];
            theta
                .iter()
                .zip(stats)
                .map(|(t, s)| t * s(&x))
                .sum::<f64>()
                .exp()
        })
        .collect()
}

fn kl(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .filter(|(x, _)| **x > 0.0)
        .map(|(x, y)| x * (x / y).ln())
        .sum()
}

fn normalised(p: &[f64]) -> Vec<f64> {
    let s: f64 = p.iter().sum();
    p.iter().map(|v| v / s).collect()
}

#[test]
fn test_projection_recovers_family_member() {
    let stats: Vec<fn(&[f64]) -> f64> = vec![stat_x];
    let ip = InformationProjection::new(stats.clone());
    let p = family_member(&[0.5], &stats, 4);
    let theta = ip.project(&p).expect("forward projection");
    assert_eq!(theta.len(), 1);
    assert!((theta[0] - 0.5).abs() < 1e-9, "forward theta = {theta:?}");
    let theta_rev = ip
        .reverse_kl_projection(&p, &[0.0])
        .expect("reverse projection");
    assert!(
        (theta_rev[0] - 0.5).abs() < 1e-9,
        "reverse theta = {theta_rev:?}"
    );
}

#[test]
fn test_projection_recovers_two_parameter_member() {
    let stats: Vec<fn(&[f64]) -> f64> = vec![stat_x, stat_x2];
    let ip = InformationProjection::new(stats.clone());
    let truth = [0.3, -0.2];
    let p = family_member(&truth, &stats, 6);
    let fwd = ip.project(&p).expect("forward");
    let rev = ip.reverse_kl_projection(&p, &[0.0, 0.0]).expect("reverse");
    for a in 0..2 {
        assert!((fwd[a] - truth[a]).abs() < 1e-8, "fwd = {fwd:?}");
        assert!((rev[a] - truth[a]).abs() < 1e-8, "rev = {rev:?}");
    }
    let q = ip.distribution(&fwd, 6).expect("distribution");
    let pn = normalised(&p);
    for (qj, pj) in q.iter().zip(&pn) {
        assert!((qj - pj).abs() < 1e-10);
    }
}

#[test]
fn test_forward_and_reverse_projection_of_out_of_family_distribution() {
    let stats: Vec<fn(&[f64]) -> f64> = vec![stat_x];
    let ip = InformationProjection::new(stats);
    // Asymmetric bimodal distribution: not a member of exp(theta x).
    let p = vec![0.45, 0.05, 0.05, 0.1, 0.35];
    let fwd = ip.project(&p).expect("forward");
    let rev = ip.reverse_kl_projection(&p, &[0.0]).expect("reverse");
    let q_fwd = ip.distribution(&fwd, p.len()).expect("q_fwd");
    let q_rev = ip.distribution(&rev, p.len()).expect("q_rev");
    // Forward projection matches the mean.
    let mean = |d: &[f64]| d.iter().enumerate().map(|(j, v)| j as f64 * v).sum::<f64>();
    assert!((mean(&q_fwd) - mean(&p)).abs() < 1e-10);
    // Each projection is optimal for its own divergence.
    assert!(kl(&p, &q_fwd) <= kl(&p, &q_rev) + 1e-12);
    assert!(kl(&q_rev, &p) <= kl(&q_fwd, &p) + 1e-12);
    // Reverse projection is stationary for KL(q_theta || p).
    let h = 1e-5;
    let f = |t: f64| {
        let q = ip.distribution(&[t], p.len()).expect("q");
        kl(&q, &p)
    };
    let grad = (f(rev[0] + h) - f(rev[0] - h)) / (2.0 * h);
    assert!(grad.abs() < 1e-8, "reverse KL gradient = {grad}");
    assert!(f(rev[0]) < f(rev[0] + 0.05) && f(rev[0]) < f(rev[0] - 0.05));
    // The two projections genuinely differ for this p.
    assert!((fwd[0] - rev[0]).abs() > 1e-3);
}

#[test]
fn test_forward_projection_with_zero_mass_in_interior() {
    let ip = InformationProjection::new(vec![stat_x]);
    // Mean 1.3 lies strictly inside [0, 3] although p_0 = p_3 = 0.
    let p = vec![0.0, 0.7, 0.3, 0.0];
    let theta = ip.project(&p).expect("forward");
    let q = ip.distribution(&theta, 4).expect("q");
    let mean: f64 = q.iter().enumerate().map(|(j, v)| j as f64 * v).sum();
    assert!((mean - 1.3).abs() < 1e-10);
}

#[test]
fn test_forward_projection_boundary_moment_is_error() {
    let ip = InformationProjection::new(vec![stat_x]);
    // E_p[x] = 0 is the minimum of the support: theta would be -infinity.
    assert!(matches!(
        ip.project(&[1.0, 0.0, 0.0, 0.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    assert!(matches!(
        ip.project(&[0.0, 0.0, 0.0, 2.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
}

#[test]
fn test_reverse_projection_rejects_zero_probability() {
    let ip = InformationProjection::new(vec![stat_x]);
    assert!(matches!(
        ip.reverse_kl_projection(&[0.5, 0.0, 0.5], &[0.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    assert!(matches!(
        ip.reverse_kl_projection(&[0.5, 0.5], &[0.0, 1.0]),
        Err(PhysicsError::InvalidInput { .. })
    ));
}

#[test]
fn test_invalid_inputs_are_rejected() {
    let ip = InformationProjection::new(vec![stat_x]);
    assert!(ip.project(&[]).is_err());
    assert!(ip.project(&[-0.1, 1.1]).is_err());
    assert!(ip.project(&[0.0, 0.0]).is_err());
    assert!(ip.project(&[f64::NAN, 1.0]).is_err());
    // Constant statistic: natural parameter not identifiable.
    let degenerate = InformationProjection::new(vec![stat_const]);
    assert!(matches!(
        degenerate.project(&[0.2, 0.3, 0.5]),
        Err(PhysicsError::InvalidInput { .. })
    ));
    // Two statistics on a two-point support are affinely dependent.
    let too_many = InformationProjection::new(vec![stat_x, stat_x2]);
    assert!(too_many.project(&[0.4, 0.6]).is_err());
    let empty = InformationProjection::new(vec![]);
    assert!(empty.project(&[0.5, 0.5]).is_err());
}

#[test]
fn test_unnormalised_input_is_normalised() {
    let ip = InformationProjection::new(vec![stat_x]);
    let p = vec![1.0, 2.0, 3.0, 4.0];
    let scaled: Vec<f64> = p.iter().map(|v| 10.0 * v).collect();
    let a = ip.project(&p).expect("p");
    let b = ip.project(&scaled).expect("scaled");
    assert!((a[0] - b[0]).abs() < 1e-12);
}
