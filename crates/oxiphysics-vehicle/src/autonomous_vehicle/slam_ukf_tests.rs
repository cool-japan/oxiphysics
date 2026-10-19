// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the FastSLAM 1.0 particle filter
//! ([`ParticleFilterSlam`]) and the unscented Kalman filter ([`UkfState`]).

use super::functions::wrap_angle;
use super::types::*;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use std::f64::consts::{PI, TAU};

// ===========================================================================
// UKF
// ===========================================================================

const CORRELATED_COV: [f64; 9] = [2.0, 0.8, 0.3, 0.8, 1.5, -0.4, 0.3, -0.4, 0.7];

fn ukf_with(mean: [f64; 3], cov: [f64; 9], alpha: f64, kappa: f64) -> UkfState {
    let mut ukf = UkfState::new();
    ukf.mean = mean;
    ukf.cov = cov;
    ukf.alpha = alpha;
    ukf.kappa = kappa;
    ukf
}

fn diagonal(a: f64, b: f64, c: f64) -> [f64; 9] {
    [a, 0.0, 0.0, 0.0, b, 0.0, 0.0, 0.0, c]
}

/// Element-wise comparison of two equally long slices.
fn assert_all_close(actual: &[f64], expected: &[f64], tol: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length mismatch");
    for (k, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!((a - e).abs() < tol, "{what}[{k}]: {a} vs {e}");
    }
}

/// Textbook linear Kalman update for the measurement `z = [x, y] + noise`.
fn reference_kalman_update(
    mean: [f64; 3],
    cov: [f64; 9],
    z: [f64; 2],
    r: f64,
) -> ([f64; 3], [f64; 9]) {
    let s = [cov[0] + r, cov[1], cov[3], cov[4] + r];
    let det = s[0] * s[3] - s[1] * s[2];
    let s_inv = [s[3] / det, -s[1] / det, -s[2] / det, s[0] / det];
    let mut gain = [[0.0; 2]; 3];
    for (i, row) in gain.iter_mut().enumerate() {
        let (a, b) = (cov[i * 3], cov[i * 3 + 1]);
        row[0] = a * s_inv[0] + b * s_inv[2];
        row[1] = a * s_inv[1] + b * s_inv[3];
    }
    let innovation = [z[0] - mean[0], z[1] - mean[1]];
    let mut new_mean = mean;
    let mut new_cov = cov;
    for i in 0..3 {
        new_mean[i] += gain[i][0] * innovation[0] + gain[i][1] * innovation[1];
        for j in 0..3 {
            // P' = P - K (H P)
            new_cov[i * 3 + j] -= gain[i][0] * cov[j] + gain[i][1] * cov[3 + j];
        }
    }
    (new_mean, new_cov)
}

#[test]
fn ukf_sigma_points_reproduce_mean_and_covariance_including_correlations() {
    for (alpha, kappa) in [(1e-3, 0.0), (0.5, 1.0), (1.0, 0.0)] {
        let ukf = ukf_with([1.0, -2.0, 0.5], CORRELATED_COV, alpha, kappa);
        let points = ukf.sigma_points().expect("valid covariance");
        let (wm, wc) = ukf.weights().expect("valid parameters");
        assert!((wm.iter().sum::<f64>() - 1.0).abs() < 1e-8, "alpha {alpha}");
        let mut mean = [0.0; 3];
        let mut cov = [0.0; 9];
        for ((point, wm_i), wc_i) in points.iter().zip(wm.iter()).zip(wc.iter()) {
            for k in 0..3 {
                mean[k] += wm_i * point[k];
            }
            for r in 0..3 {
                for c in 0..3 {
                    cov[r * 3 + c] += wc_i * (point[r] - ukf.mean[r]) * (point[c] - ukf.mean[c]);
                }
            }
        }
        assert_all_close(&mean, &ukf.mean, 1e-7, &format!("alpha {alpha}: mean"));
        assert_all_close(&cov, &CORRELATED_COV, 1e-7, &format!("alpha {alpha}: cov"));
    }
}

#[test]
fn ukf_weights_follow_the_scaled_unscented_transform() {
    let ukf = ukf_with([0.0; 3], diagonal(1.0, 1.0, 1.0), 0.5, 1.0);
    let (wm, wc) = ukf.weights().expect("valid parameters");
    let n = 3.0;
    let lambda = 0.5 * 0.5 * (n + 1.0) - n;
    assert!((wm[0] - lambda / (n + lambda)).abs() < 1e-12);
    assert!((wc[0] - (lambda / (n + lambda) + 1.0 - 0.25 + 2.0)).abs() < 1e-12);
    for i in 1..7 {
        assert!((wm[i] - 0.5 / (n + lambda)).abs() < 1e-12);
        assert!((wc[i] - 0.5 / (n + lambda)).abs() < 1e-12);
    }
}

#[test]
fn ukf_sigma_points_honour_off_diagonal_covariance() {
    // Strongly correlated x/y: the sigma points must lie along the diagonal,
    // which the former "diagonal sqrt only" construction could not represent.
    let cov = [1.0, 0.9, 0.0, 0.9, 1.0, 0.0, 0.0, 0.0, 0.1];
    let ukf = ukf_with([0.0; 3], cov, 1.0, 0.0);
    let points = ukf.sigma_points().expect("valid covariance");
    for point in &points[1..] {
        let cross = point[0] * point[1];
        assert!(
            cross >= -1e-12,
            "x and y must not have opposite signs: {point:?}"
        );
    }
    // sqrt(3) * (1, 0.9): the first Cholesky column lies along the diagonal.
    assert!(
        points.iter().any(|p| p[0] * p[1] > 2.0),
        "no sigma point along the correlated direction: {points:?}"
    );
}

#[test]
fn ukf_sigma_points_reject_invalid_covariances_and_parameters() {
    let mut ukf = UkfState::new();
    ukf.cov = diagonal(1.0, -1.0, 1.0);
    assert!(ukf.sigma_points().is_err());
    ukf.cov = [1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(ukf.sigma_points().is_err(), "indefinite covariance");
    ukf.cov = diagonal(f64::NAN, 1.0, 1.0);
    assert!(ukf.sigma_points().is_err());

    let mut ukf = UkfState::new();
    ukf.alpha = 0.0;
    assert!(ukf.sigma_points().is_err());
    assert!(ukf.weights().is_err());
    ukf.alpha = 1.0;
    ukf.kappa = -3.0;
    assert!(ukf.weights().is_err(), "n + kappa = 0 gives a zero spread");
}

#[test]
fn ukf_sigma_points_of_a_certain_state_collapse_to_the_mean() {
    let mean = [3.0, -1.0, 0.25];
    let ukf = ukf_with(mean, [0.0; 9], 1e-3, 0.0);
    for point in ukf.sigma_points().expect("zero covariance is valid") {
        assert_eq!(point, mean);
    }
    // Uncertainty in x only: only the x coordinate varies.
    let ukf = ukf_with(mean, diagonal(4.0, 0.0, 0.0), 1.0, 0.0);
    for point in ukf
        .sigma_points()
        .expect("semi-definite covariance is valid")
    {
        assert_eq!(point[1], mean[1]);
        assert_eq!(point[2], mean[2]);
    }
}

#[test]
fn ukf_predict_follows_dead_reckoning_for_straight_motion() {
    for yaw0 in [0.0, PI / 6.0, -2.0] {
        let mut ukf = ukf_with([0.0, 0.0, yaw0], diagonal(1e-10, 1e-10, 1e-10), 1e-3, 0.0);
        for _ in 0..10 {
            ukf.predict(0.1, 10.0, 0.0, [0.0; 3])
                .expect("valid prediction");
        }
        assert!(
            (ukf.mean[0] - 10.0 * yaw0.cos()).abs() < 1e-6,
            "yaw0 {yaw0}: x = {}",
            ukf.mean[0]
        );
        assert!((ukf.mean[1] - 10.0 * yaw0.sin()).abs() < 1e-6);
        assert!((ukf.mean[2] - yaw0).abs() < 1e-9);
    }
}

#[test]
fn ukf_predict_turning_matches_euler_integration_of_the_unicycle() {
    let (v, omega, dt) = (5.0, 0.4, 0.05);
    let mut ukf = ukf_with([1.0, -2.0, 0.3], diagonal(1e-12, 1e-12, 1e-12), 1e-3, 0.0);
    let mut expected: [f64; 3] = [1.0, -2.0, 0.3];
    for _ in 0..40 {
        ukf.predict(dt, v, omega, [0.0; 3])
            .expect("valid prediction");
        expected = [
            expected[0] + v * expected[2].cos() * dt,
            expected[1] + v * expected[2].sin() * dt,
            expected[2] + omega * dt,
        ];
    }
    assert_all_close(&ukf.mean, &expected, 1e-6, "state");
}

#[test]
fn ukf_predict_captures_second_order_heading_uncertainty() {
    // y ~ N(0, 0.5), yaw ~ N(0, 0.04): a linearised (EKF) prediction would give
    // x = v dt = 10 and var(y) = 0.5 + v² dt² var(yaw) only to first order in the
    // mean.  The unscented transform sees E[cos(yaw)] = 1 - var/2 = 0.98.
    let mut ukf = ukf_with([0.0; 3], diagonal(1e-9, 0.5, 0.04), 1e-3, 0.0);
    ukf.predict(1.0, 10.0, 0.0, [0.0; 3])
        .expect("valid prediction");
    assert!(
        (ukf.mean[0] - 9.8).abs() < 1e-5,
        "x = {} (EKF would give 10.0)",
        ukf.mean[0]
    );
    assert!(ukf.mean[0] < 9.81);
    assert!(ukf.mean[1].abs() < 1e-9 && ukf.mean[2].abs() < 1e-9);
    // var(y') = var(y) + v² dt² var(yaw)
    assert!((ukf.cov[4] - 4.5).abs() < 1e-4, "var(y) = {}", ukf.cov[4]);
    // cov(y', yaw') = v dt var(yaw) = 0.4
    assert!(
        (ukf.cov[5] - 0.4).abs() < 1e-5,
        "cov(y, yaw) = {}",
        ukf.cov[5]
    );
    assert!(
        (ukf.cov[7] - ukf.cov[5]).abs() < 1e-12,
        "covariance symmetric"
    );
    assert!((ukf.cov[8] - 0.04).abs() < 1e-8);
}

#[test]
fn ukf_predict_adds_process_noise_to_the_diagonal() {
    let mut ukf = ukf_with([0.0; 3], diagonal(1e-9, 1e-9, 1e-9), 1e-3, 0.0);
    ukf.predict(0.1, 0.0, 0.0, [0.01, 0.02, 0.03])
        .expect("valid prediction");
    assert!((ukf.cov[0] - 0.01).abs() < 1e-6);
    assert!((ukf.cov[4] - 0.02).abs() < 1e-6);
    assert!((ukf.cov[8] - 0.03).abs() < 1e-6);
}

#[test]
fn ukf_predict_wraps_yaw_across_the_branch_cut() {
    // alpha = 1 spreads the sigma points ±0.17 rad, so some of them cross pi.
    let mut ukf = ukf_with([0.0, 0.0, 3.0], diagonal(1e-6, 1e-6, 0.01), 1.0, 0.0);
    ukf.predict(1.0, 0.0, 0.2, [0.0; 3])
        .expect("valid prediction");
    let expected = wrap_angle(3.2);
    assert!(
        (ukf.mean[2] - expected).abs() < 1e-9,
        "yaw = {}",
        ukf.mean[2]
    );
    assert!(ukf.mean[2] < 0.0);
    assert!(
        (ukf.cov[8] - 0.01).abs() < 1e-6,
        "yaw variance must not pick up a 2π jump: {}",
        ukf.cov[8]
    );
    let mut ukf = ukf_with([0.0, 0.0, -3.1], diagonal(1e-6, 1e-6, 0.01), 1.0, 0.0);
    ukf.predict(1.0, 0.0, -0.2, [0.0; 3])
        .expect("valid prediction");
    assert!((ukf.mean[2] - wrap_angle(-3.3)).abs() < 1e-9);
    assert!((ukf.cov[8] - 0.01).abs() < 1e-6);
}

#[test]
fn ukf_update_matches_the_linear_kalman_filter_for_position_measurements() {
    let prior_mean = [0.3, -0.2, 0.1];
    let prior_cov = [1.0, 0.4, 0.2, 0.4, 2.0, -0.3, 0.2, -0.3, 0.5];
    let (z, r) = ([1.0, -0.5], 0.5);
    let (ref_mean, ref_cov) = reference_kalman_update(prior_mean, prior_cov, z, r);
    for (alpha, kappa) in [(1e-3, 0.0), (0.5, 1.0), (1.0, 0.0)] {
        let mut ukf = ukf_with(prior_mean, prior_cov, alpha, kappa);
        ukf.update(z, r).expect("valid update");
        assert_all_close(&ukf.mean, &ref_mean, 1e-7, &format!("alpha {alpha}: mean"));
        assert_all_close(&ukf.cov, &ref_cov, 1e-7, &format!("alpha {alpha}: cov"));
        // The unmeasured yaw is corrected through the cross-covariance.
        assert!((ukf.mean[2] - prior_mean[2]).abs() > 1e-3);
    }
}

#[test]
fn ukf_update_shrinks_uncertainty_and_keeps_the_covariance_symmetric() {
    let mut ukf = ukf_with([0.0; 3], CORRELATED_COV, 0.5, 0.0);
    let before = ukf.cov;
    ukf.update([0.4, -0.3], 0.25).expect("valid update");
    for i in 0..3 {
        for j in 0..3 {
            assert!((ukf.cov[i * 3 + j] - ukf.cov[j * 3 + i]).abs() < 1e-12);
        }
    }
    assert!(ukf.cov[0] < before[0] && ukf.cov[4] < before[4]);
    assert!(
        ukf.cov[8] < before[8],
        "yaw variance shrinks via x/y correlation"
    );
    assert!(
        ukf.sigma_points().is_ok(),
        "posterior stays positive semi-definite"
    );
}

#[test]
fn ukf_position_measurement_corrects_heading_after_a_prediction() {
    // After a prediction with heading uncertainty, lateral position and heading
    // are positively correlated; a measurement to the left of the predicted
    // position must therefore rotate the heading estimate to the left.
    let mut ukf = ukf_with([0.0; 3], diagonal(0.01, 0.01, 0.04), 1e-3, 0.0);
    ukf.predict(1.0, 10.0, 0.0, [0.0; 3])
        .expect("valid prediction");
    assert!(ukf.cov[5] > 0.3);
    let prior_yaw_var = ukf.cov[8];
    let (ref_mean, ref_cov) = reference_kalman_update(ukf.mean, ukf.cov, [9.8, 0.5], 0.25);
    ukf.update([9.8, 0.5], 0.25).expect("valid update");
    assert!(ukf.mean[2] > 0.02, "yaw = {}", ukf.mean[2]);
    assert!(ukf.cov[8] < prior_yaw_var);
    assert_all_close(&ukf.mean, &ref_mean, 1e-6, "mean");
    assert_all_close(&ukf.cov, &ref_cov, 1e-6, "cov");
}

#[test]
fn ukf_tracks_a_moving_vehicle_from_noisy_positions() {
    let mut rng = StdRng::seed_from_u64(11);
    let mut ukf = ukf_with([0.0, 0.0, 0.0], diagonal(1.0, 1.0, 0.25), 0.5, 0.0);
    let (v, omega, dt) = (8.0, 0.1, 0.1);
    let mut truth: [f64; 3] = [0.0, 0.0, 0.0];
    for _ in 0..100 {
        truth = [
            truth[0] + v * truth[2].cos() * dt,
            truth[1] + v * truth[2].sin() * dt,
            truth[2] + omega * dt,
        ];
        ukf.predict(dt, v, omega, [1e-4, 1e-4, 1e-5])
            .expect("valid prediction");
        let noisy = [
            truth[0] + rng.random_range(-0.3..0.3),
            truth[1] + rng.random_range(-0.3..0.3),
        ];
        ukf.update(noisy, 0.03).expect("valid update");
    }
    assert!((ukf.mean[0] - truth[0]).abs() < 0.3, "x error");
    assert!((ukf.mean[1] - truth[1]).abs() < 0.3, "y error");
    assert!(
        wrap_angle(ukf.mean[2] - truth[2]).abs() < 0.1,
        "heading error {}",
        wrap_angle(ukf.mean[2] - truth[2])
    );
}

#[test]
fn ukf_rejects_invalid_arguments_without_touching_the_state() {
    let mut ukf = ukf_with([1.0, 2.0, 0.3], CORRELATED_COV, 0.5, 0.0);
    let (mean, cov) = (ukf.mean, ukf.cov);
    assert!(ukf.predict(-0.1, 1.0, 0.0, [0.0; 3]).is_err());
    assert!(ukf.predict(f64::NAN, 1.0, 0.0, [0.0; 3]).is_err());
    assert!(ukf.predict(0.1, f64::INFINITY, 0.0, [0.0; 3]).is_err());
    assert!(ukf.predict(0.1, 1.0, f64::NAN, [0.0; 3]).is_err());
    assert!(ukf.predict(0.1, 1.0, 0.0, [-1.0, 0.0, 0.0]).is_err());
    assert!(ukf.predict(0.1, 1.0, 0.0, [0.0, f64::NAN, 0.0]).is_err());
    assert!(ukf.update([f64::NAN, 0.0], 1.0).is_err());
    assert!(ukf.update([0.0, 0.0], -1.0).is_err());
    assert!(ukf.update([0.0, 0.0], f64::NAN).is_err());
    assert_eq!(ukf.mean, mean);
    assert_eq!(ukf.cov, cov);

    // A non-positive-semi-definite covariance cannot be propagated or updated.
    ukf.cov = diagonal(1.0, -1.0, 1.0);
    assert!(ukf.predict(0.1, 1.0, 0.0, [0.0; 3]).is_err());
    assert!(ukf.update([0.0, 0.0], 1.0).is_err());
    assert_eq!(ukf.mean, mean);

    // No uncertainty and a perfect sensor: the innovation covariance is singular.
    let mut ukf = ukf_with([0.0; 3], [0.0; 9], 1.0, 0.0);
    assert!(ukf.update([1.0, 1.0], 0.0).is_err());
    assert_eq!(ukf.mean, [0.0; 3]);
}

// ===========================================================================
// FastSLAM
// ===========================================================================

const SIGMA_R: f64 = 0.1;
const SIGMA_B: f64 = 0.05;

/// Range and bearing (relative to the heading) of `landmark` seen from `pose`.
fn range_bearing(pose: (f64, f64, f64), landmark: Vec2) -> (f64, f64) {
    let dx = landmark.x - pose.0;
    let dy = landmark.y - pose.1;
    (dx.hypot(dy), wrap_angle(dy.atan2(dx) - pose.2))
}

/// Noise-free motion with the model used by `SlamParticle::propagate`.
fn advance(pose: (f64, f64, f64), delta_s: f64, delta_theta: f64) -> (f64, f64, f64) {
    let yaw = pose.2 + delta_theta;
    (
        pose.0 + delta_s * yaw.cos(),
        pose.1 + delta_s * yaw.sin(),
        yaw,
    )
}

fn scan(pose: (f64, f64, f64), landmarks: &[Vec2]) -> Vec<(usize, f64, f64)> {
    landmarks
        .iter()
        .enumerate()
        .map(|(id, &lm)| {
            let (r, b) = range_bearing(pose, lm);
            (id, r, b)
        })
        .collect()
}

fn five_landmarks() -> Vec<Vec2> {
    vec![
        Vec2::new(6.0, 4.0),
        Vec2::new(10.0, -3.0),
        Vec2::new(14.0, 5.0),
        Vec2::new(8.0, 8.0),
        Vec2::new(3.0, -5.0),
    ]
}

#[test]
fn slam_first_sighting_initialises_a_landmark_with_the_inverse_jacobian_covariance() {
    let mut particle = SlamParticle::new(1.0, 2.0, 0.5, 0.25);
    particle.update_weight(&[(3, 4.0, 0.2)], SIGMA_R, SIGMA_B);
    assert_eq!(particle.landmarks.len(), 1);
    let lm = particle.landmarks[0];
    let phi: f64 = 0.5 + 0.2;
    assert_eq!(lm.id, 3);
    assert!((lm.position.x - (1.0 + 4.0 * phi.cos())).abs() < 1e-12);
    assert!((lm.position.y - (2.0 + 4.0 * phi.sin())).abs() < 1e-12);
    // cov = sigma_r² u uᵀ + (r sigma_b)² v vᵀ with u along, v across the ray.
    let (u, v) = ([phi.cos(), phi.sin()], [-phi.sin(), phi.cos()]);
    let across = (4.0 * SIGMA_B).powi(2);
    let expected = [
        SIGMA_R * SIGMA_R * u[0] * u[0] + across * v[0] * v[0],
        SIGMA_R * SIGMA_R * u[0] * u[1] + across * v[0] * v[1],
        SIGMA_R * SIGMA_R * u[1] * u[0] + across * v[1] * v[0],
        SIGMA_R * SIGMA_R * u[1] * u[1] + across * v[1] * v[1],
    ];
    assert_all_close(&lm.cov, &expected, 1e-12, "landmark covariance");
    assert_eq!(
        particle.weight, 0.25,
        "a first sighting carries no information"
    );
}

#[test]
fn slam_second_sighting_applies_the_ekf_update_and_the_innovation_likelihood() {
    // Robot at the origin looking at a landmark 10 m ahead: G = diag(1, 0.1).
    let mut particle = SlamParticle::new(0.0, 0.0, 0.0, 1.0);
    particle.observe(&[(0, 10.0, 0.0)], SIGMA_R, SIGMA_B);
    let first = particle.landmarks[0];
    assert!((first.cov[0] - SIGMA_R * SIGMA_R).abs() < 1e-12);
    assert!((first.cov[3] - 100.0 * SIGMA_B * SIGMA_B).abs() < 1e-12);

    let log_likelihood = particle.observe(&[(0, 10.2, 0.001)], SIGMA_R, SIGMA_B);
    // S = diag(2 sigma_r², 2 sigma_b²), nu = (0.2, 0.001), K = diag(0.5, 5).
    let s = [2.0 * SIGMA_R * SIGMA_R, 2.0 * SIGMA_B * SIGMA_B];
    let nu = [0.2, 0.001];
    let expected_ll =
        -TAU.ln() - 0.5 * (s[0] * s[1]).ln() - 0.5 * (nu[0] * nu[0] / s[0] + nu[1] * nu[1] / s[1]);
    assert!(
        (log_likelihood - expected_ll).abs() < 1e-9,
        "log-likelihood {log_likelihood} vs {expected_ll}"
    );
    let lm = particle.landmarks[0];
    assert_eq!(particle.landmarks.len(), 1);
    assert!((lm.position.x - 10.1).abs() < 1e-9, "x = {}", lm.position.x);
    assert!(
        (lm.position.y - 0.005).abs() < 1e-9,
        "y = {}",
        lm.position.y
    );
    // Joseph-form covariance: exactly half of the prior for this geometry.
    assert!((lm.cov[0] - 0.5 * first.cov[0]).abs() < 1e-12);
    assert!((lm.cov[3] - 0.5 * first.cov[3]).abs() < 1e-12);
    assert!(lm.cov[1].abs() < 1e-12 && lm.cov[2].abs() < 1e-12);
}

#[test]
fn slam_update_weight_is_no_longer_a_no_op() {
    let mut particle = SlamParticle::new(0.0, 0.0, 0.0, 0.5);
    particle.update_weight(&[(0, 10.0, 0.0)], SIGMA_R, SIGMA_B);
    assert_eq!(particle.weight, 0.5);
    particle.update_weight(&[(0, 10.2, 0.001)], SIGMA_R, SIGMA_B);
    let s = [2.0 * SIGMA_R * SIGMA_R, 2.0 * SIGMA_B * SIGMA_B];
    let ll = -TAU.ln()
        - 0.5 * (s[0] * s[1]).ln()
        - 0.5 * (0.2_f64.powi(2) / s[0] + 0.001_f64.powi(2) / s[1]);
    assert!((particle.weight - 0.5 * ll.exp()).abs() < 1e-9);
    assert_ne!(particle.weight, 0.5);
}

#[test]
fn slam_particles_with_a_wrong_pose_lose_weight() {
    let landmark = Vec2::new(10.0, 0.0);
    let mut slam = ParticleFilterSlam::new(2, 0.0, 0.0, 0.1).with_bearing_sigma(0.05);
    for p in &mut slam.particles {
        p.landmarks.push(Landmark::known(0, landmark));
    }
    slam.particles[1].pose = (1.0, 0.0, 0.0);
    let (range, bearing) = range_bearing((0.0, 0.0, 0.0), landmark);
    slam.update(&[(0, range, bearing)]);
    let (w_good, w_bad) = (slam.particles[0].weight, slam.particles[1].weight);
    assert!(w_good > 0.999_999, "good hypothesis weight {w_good}");
    assert!(w_bad < 1e-15, "bad hypothesis weight {w_bad}");
    assert!((w_good + w_bad - 1.0).abs() < 1e-12);
    assert_eq!(slam.best_particle().map(|p| p.pose), Some((0.0, 0.0, 0.0)));
}

#[test]
fn slam_weight_scaling_survives_scans_that_underflow_every_likelihood() {
    // Every particle is ~900 sigma away from the measurement: the plain
    // likelihoods are all exp(-4e5) = 0 in double precision.
    let mut slam = ParticleFilterSlam::new(3, 0.0, 0.0, 0.1);
    for (i, p) in slam.particles.iter_mut().enumerate() {
        p.landmarks.push(Landmark::known(0, Vec2::new(10.0, 0.0)));
        p.pose = (0.5 * i as f64, 0.0, 0.0);
    }
    slam.update(&[(0, 100.0, 0.0)]);
    let weights: Vec<f64> = slam.particles.iter().map(|p| p.weight).collect();
    assert!(weights.iter().all(|w| w.is_finite()));
    assert!((weights.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!((weights[0] - 1.0).abs() < 1e-12, "weights {weights:?}");
}

#[test]
fn slam_degenerate_weights_reset_to_uniform_before_resampling() {
    for bad in [0.0, f64::NAN, 1e-310] {
        let mut slam = ParticleFilterSlam::new(4, 0.0, 0.0, 0.1);
        for (i, p) in slam.particles.iter_mut().enumerate() {
            p.weight = bad;
            p.pose.0 = i as f64;
        }
        slam.normalise_weights();
        assert!(
            slam.particles
                .iter()
                .all(|p| (p.weight - 0.25).abs() < 1e-12),
            "bad weight {bad}"
        );
        for p in &mut slam.particles {
            p.weight = bad;
        }
        let mut rng = StdRng::seed_from_u64(5);
        slam.resample_with_rng(&mut rng);
        assert_eq!(slam.particles.len(), 4);
        assert!(
            slam.particles
                .iter()
                .all(|p| (p.weight - 0.25).abs() < 1e-12)
        );
        // A vanishing total used to pin every sample on the last particle.
        let distinct: std::collections::HashSet<u64> =
            slam.particles.iter().map(|p| p.pose.0.to_bits()).collect();
        assert!(
            distinct.len() > 1,
            "bad weight {bad}: degenerate resampling"
        );
    }
}

#[test]
fn slam_systematic_resampling_copies_particles_in_proportion_to_their_weights() {
    let mut slam = ParticleFilterSlam::new(100, 0.0, 0.0, 0.1);
    for (i, p) in slam.particles.iter_mut().enumerate() {
        p.pose.0 = i as f64;
        p.weight = if i == 7 { 0.5 } else { 0.5 / 99.0 };
    }
    let mut rng = StdRng::seed_from_u64(99);
    slam.resample_with_rng(&mut rng);
    assert_eq!(slam.particles.len(), 100);
    let copies = slam.particles.iter().filter(|p| p.pose.0 == 7.0).count();
    assert!(
        (49..=51).contains(&copies),
        "copies of particle 7: {copies}"
    );
    assert!(
        slam.particles
            .iter()
            .all(|p| (p.weight - 0.01).abs() < 1e-12)
    );
    assert!((slam.effective_sample_size() - 100.0).abs() < 1e-9);
}

#[test]
fn slam_bearing_innovation_is_wrapped_across_the_branch_cut() {
    let mut particle = SlamParticle::new(0.0, 0.0, 0.0, 1.0);
    particle.observe(&[(0, 10.0, PI)], SIGMA_R, SIGMA_B);
    // The same direction expressed as -pi + 0.001 differs by 0.001 rad, not 2 pi.
    let ll = particle.observe(&[(0, 10.0, -PI + 0.001)], SIGMA_R, SIGMA_B);
    let lm = particle.landmarks[0];
    assert!((lm.position.x + 10.0).abs() < 1e-6, "x = {}", lm.position.x);
    assert!(
        (lm.position.y + 0.005).abs() < 1e-6,
        "y = {}",
        lm.position.y
    );
    assert!(ll > -10.0, "a 0.001 rad innovation is likely: ll = {ll}");
}

#[test]
fn slam_skips_invalid_and_degenerate_observations() {
    let mut particle = SlamParticle::new(3.0, 4.0, 0.3, 1.0);
    let ll = particle.observe(
        &[
            (0, f64::NAN, 0.0),
            (1, -2.0, 0.0),
            (2, 5.0, f64::INFINITY),
            (3, f64::INFINITY, 0.0),
        ],
        SIGMA_R,
        SIGMA_B,
    );
    assert_eq!(ll, 0.0);
    assert!(particle.landmarks.is_empty());
    // A mapped landmark exactly at the sensor cannot be linearised.
    particle
        .landmarks
        .push(Landmark::known(9, Vec2::new(3.0, 4.0)));
    let ll = particle.observe(&[(9, 0.0, 0.0)], SIGMA_R, SIGMA_B);
    assert_eq!(ll, 0.0);
    assert_eq!(particle.landmarks[0].position, Vec2::new(3.0, 4.0));
}

#[test]
fn slam_zero_motion_noise_is_deterministic_and_does_not_panic() {
    let mut slam = ParticleFilterSlam::new(5, 0.0, 0.0, 0.1);
    let mut rng = StdRng::seed_from_u64(1);
    slam.predict_with_rng(1.0, 0.0, &mut rng);
    for p in &slam.particles {
        assert_eq!(p.pose, (1.0, 0.0, 0.0));
    }
    slam.predict(0.5, 0.0);
    for p in &slam.particles {
        assert!((p.pose.0 - 1.5).abs() < 1e-12);
    }
}

#[test]
fn slam_seeded_runs_are_reproducible() {
    let run = |seed: u64| {
        let mut slam = ParticleFilterSlam::new(20, 0.05, 0.02, 0.1);
        let mut rng = StdRng::seed_from_u64(seed);
        for _ in 0..10 {
            slam.predict_with_rng(0.3, 0.05, &mut rng);
        }
        slam.particles.iter().map(|p| p.pose).collect::<Vec<_>>()
    };
    assert_eq!(run(21), run(21));
    assert_ne!(run(21), run(22));
}

#[test]
fn slam_noise_free_odometry_recovers_the_exact_map() {
    let landmarks = five_landmarks();
    let mut slam = ParticleFilterSlam::new(10, 0.0, 0.0, 0.05).with_bearing_sigma(0.01);
    let mut rng = StdRng::seed_from_u64(3);
    let mut truth = (0.0, 0.0, 0.0);
    slam.update(&scan(truth, &landmarks));
    for _ in 0..30 {
        truth = advance(truth, 0.4, 0.03);
        slam.predict_with_rng(0.4, 0.03, &mut rng);
        slam.update(&scan(truth, &landmarks));
    }
    for (id, lm) in landmarks.iter().enumerate() {
        let estimate = slam.landmark_estimate(id).expect("landmark was observed");
        assert!(
            estimate.dist(*lm) < 1e-9,
            "landmark {id}: {estimate:?} vs {lm:?}"
        );
    }
    let (x, y, yaw) = slam.mean_pose();
    assert!((x - truth.0).abs() < 1e-9 && (y - truth.1).abs() < 1e-9);
    assert!((yaw - truth.2).abs() < 1e-9);
    assert!((slam.effective_sample_size() - 10.0).abs() < 1e-6);
}

#[test]
fn slam_builds_a_map_of_five_landmarks_from_known_odometry() {
    let landmarks = five_landmarks();
    for seed in [1_u64, 2, 3] {
        let n = 150;
        let mut slam = ParticleFilterSlam::new(n, 0.01, 0.004, 0.05).with_bearing_sigma(0.01);
        let mut rng = StdRng::seed_from_u64(seed);
        let noisy_scan = |pose: (f64, f64, f64), rng: &mut StdRng| -> Vec<(usize, f64, f64)> {
            scan(pose, &landmarks)
                .into_iter()
                .map(|(id, r, b)| {
                    (
                        id,
                        r + rng.random_range(-0.06..0.06),
                        b + rng.random_range(-0.012..0.012),
                    )
                })
                .collect()
        };
        let mut truth = (0.0, 0.0, 0.0);
        slam.update(&noisy_scan(truth, &mut rng));
        for _ in 0..50 {
            truth = advance(truth, 0.4, 0.03);
            slam.predict_with_rng(0.4, 0.03, &mut rng);
            slam.update(&noisy_scan(truth, &mut rng));
            if slam.effective_sample_size() < 0.5 * n as f64 {
                slam.resample_with_rng(&mut rng);
            }
        }
        for (id, lm) in landmarks.iter().enumerate() {
            let estimate = slam.landmark_estimate(id).expect("landmark was observed");
            let error = estimate.dist(*lm);
            assert!(error < 0.2, "seed {seed}, landmark {id}: error {error} m");
        }
        let (x, y, _) = slam.mean_pose();
        let pose_error = ((x - truth.0).powi(2) + (y - truth.1).powi(2)).sqrt();
        assert!(pose_error < 0.2, "seed {seed}: pose error {pose_error} m");
        // The landmark covariances of the best particle shrank with repeated sightings.
        let best = slam.best_particle().expect("particles exist");
        assert_eq!(best.landmarks.len(), 5);
        for lm in &best.landmarks {
            assert!(lm.cov[0] < 0.01 && lm.cov[3] < 0.01, "cov {:?}", lm.cov);
        }
    }
}

#[test]
fn slam_landmark_lookup_and_best_particle_handle_empty_maps() {
    let slam = ParticleFilterSlam::new(3, 0.0, 0.0, 0.1);
    assert!(slam.landmark_estimate(0).is_none());
    assert!(slam.best_particle().is_some());
    assert!((slam.effective_sample_size() - 3.0).abs() < 1e-12);
    let empty = ParticleFilterSlam::new(0, 0.0, 0.0, 0.1);
    assert!(empty.best_particle().is_none());
    assert!(empty.landmark_estimate(0).is_none());
}
