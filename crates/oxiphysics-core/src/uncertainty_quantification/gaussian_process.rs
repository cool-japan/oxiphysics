// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact Gaussian-process regression with a squared-exponential kernel.
//!
//! Fitting factorises the covariance matrix `K + σ_n² I` by Cholesky
//! decomposition, `K + σ_n² I = L Lᵀ`. The posterior at a new point `x*` is
//!
//! * mean `k*ᵀ α` with `α = (K + σ_n² I)⁻¹ y`, and
//! * variance `k(x*, x*) − ‖L⁻¹ k*‖²`,
//!
//! where `k*` holds the kernel values between `x*` and the training inputs.

use crate::error::{PhysicsError, PhysicsResult};
use nalgebra::{DMatrix, DVector};

/// Relative diagonal jitters (multiples of `σ_f²`) tried, in order, when the
/// covariance matrix is not numerically positive definite. A first attempt
/// without any jitter is always made.
const JITTER_LADDER: [f64; 3] = [1e-10, 1e-8, 1e-6];

/// Squared-exponential (RBF) kernel: k(x, x') = σ² exp(-‖x−x'‖² / (2ℓ²)).
///
/// # Arguments
/// * `x` - First input vector.
/// * `xp` - Second input vector.
/// * `sigma_f` - Signal standard deviation σ_f.
/// * `length_scale` - Length scale ℓ.
pub fn rbf_kernel(x: &[f64], xp: &[f64], sigma_f: f64, length_scale: f64) -> f64 {
    let sq_dist: f64 = x.iter().zip(xp.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    sigma_f * sigma_f * (-sq_dist / (2.0 * length_scale * length_scale)).exp()
}

/// A fitted Gaussian Process regression model.
#[derive(Debug, Clone)]
pub struct GaussianProcessSurrogate {
    /// Training input points (row-major: n_train × n_dims).
    pub x_train: Vec<Vec<f64>>,
    /// Training output values.
    pub y_train: Vec<f64>,
    /// Alpha vector = (K + σ²_n I)^{-1} y  (pre-computed for prediction).
    pub alpha: Vec<f64>,
    /// Signal standard deviation σ_f.
    pub sigma_f: f64,
    /// Length scale ℓ.
    pub length_scale: f64,
    /// Noise standard deviation σ_n.
    pub noise_std: f64,
    /// Lower Cholesky factor `L` of `K + σ_n² I + jitter·I` (`n_train × n_train`).
    pub chol_l: DMatrix<f64>,
    /// Diagonal jitter that had to be added to make the covariance matrix
    /// positive definite; `0.0` when the plain matrix factorised.
    pub jitter: f64,
}

impl GaussianProcessSurrogate {
    /// A surrogate without data: its posterior is the prior, with mean `0` and
    /// variance `σ_f²` everywhere.
    fn prior(sigma_f: f64, length_scale: f64, noise_std: f64) -> Self {
        Self {
            x_train: Vec::new(),
            y_train: Vec::new(),
            alpha: Vec::new(),
            sigma_f,
            length_scale,
            noise_std,
            chol_l: DMatrix::zeros(0, 0),
            jitter: 0.0,
        }
    }

    /// Predict the posterior mean and variance of the latent function at `x_star`.
    ///
    /// The variance is `k(x*, x*) − ‖L⁻¹ k*‖²` clamped at zero; it describes the
    /// uncertainty of the underlying function and does **not** include the
    /// observation noise (add `noise_std²` for the variance of a new noisy
    /// measurement). It approaches `σ_f²` far away from the data and drops to
    /// roughly `σ_n²` or below at a training input.
    pub fn predict(&self, x_star: &[f64]) -> (f64, f64) {
        let k_ss = rbf_kernel(x_star, x_star, self.sigma_f, self.length_scale);
        let n = self.x_train.len();
        if n == 0 || self.alpha.len() != n || self.chol_l.nrows() != n || self.chol_l.ncols() != n {
            return (0.0, k_ss);
        }
        let k_star = DVector::from_iterator(
            n,
            self.x_train
                .iter()
                .map(|xi| rbf_kernel(xi, x_star, self.sigma_f, self.length_scale)),
        );
        let mean_pred: f64 = k_star
            .iter()
            .zip(self.alpha.iter())
            .map(|(k, a)| k * a)
            .sum();
        // v = L⁻¹ k*, so ‖v‖² = k*ᵀ (K + σ²I)⁻¹ k*.
        let var_pred = match self.chol_l.solve_lower_triangular(&k_star) {
            Some(v) => (k_ss - v.dot(&v)).max(0.0),
            // A zero pivot can only appear in a hand-edited factor; the prior
            // variance is the safe upper bound.
            None => k_ss,
        };
        (mean_pred, var_pred)
    }
}

/// Reject non-finite or ill-shaped inputs of a GP fit.
fn validate_gp_inputs(
    x_train: &[Vec<f64>],
    y_train: &[f64],
    sigma_f: f64,
    length_scale: f64,
    noise_std: f64,
) -> PhysicsResult<()> {
    let invalid = |field: &str, reason: String| PhysicsError::InvalidInput {
        field: field.to_string(),
        reason,
    };
    if x_train.len() != y_train.len() {
        return Err(invalid(
            "y_train",
            format!(
                "expected one output per training input ({}), got {}",
                x_train.len(),
                y_train.len()
            ),
        ));
    }
    if !sigma_f.is_finite() || sigma_f <= 0.0 {
        return Err(invalid(
            "sigma_f",
            format!("must be finite and positive, got {sigma_f}"),
        ));
    }
    if !length_scale.is_finite() || length_scale <= 0.0 {
        return Err(invalid(
            "length_scale",
            format!("must be finite and positive, got {length_scale}"),
        ));
    }
    if !noise_std.is_finite() || noise_std < 0.0 {
        return Err(invalid(
            "noise_std",
            format!("must be finite and non-negative, got {noise_std}"),
        ));
    }
    if y_train.iter().any(|y| !y.is_finite()) {
        return Err(invalid("y_train", "all outputs must be finite".to_string()));
    }
    let dims = x_train.first().map_or(0, Vec::len);
    if x_train.iter().any(|x| x.len() != dims) {
        return Err(invalid(
            "x_train",
            format!("all training inputs must have {dims} coordinates"),
        ));
    }
    if x_train.iter().flatten().any(|c| !c.is_finite()) {
        return Err(invalid(
            "x_train",
            "all coordinates must be finite".to_string(),
        ));
    }
    Ok(())
}

/// Fit a Gaussian Process surrogate, reporting failures.
///
/// Builds `K + σ_n² I`, factorises it by Cholesky decomposition and solves
/// `(K + σ_n² I) α = y`. If the matrix is not numerically positive definite
/// (typically because training inputs nearly coincide and the noise is tiny)
/// the factorisation is retried with a diagonal jitter of `σ_f²·{1e-10, 1e-8,
/// 1e-6}`; the jitter that was needed is recorded in
/// [`GaussianProcessSurrogate::jitter`].
///
/// # Arguments
/// * `x_train` - Training inputs (one vector per sample, all of equal length).
/// * `y_train` - Training outputs.
/// * `sigma_f` - Signal standard deviation (must be positive).
/// * `length_scale` - RBF length scale (must be positive).
/// * `noise_std` - Observation noise standard deviation (non-negative).
///
/// # Errors
///
/// [`PhysicsError::InvalidInput`] if the inputs are not finite or have
/// inconsistent shapes, if a hyper-parameter is out of range, or if the
/// covariance matrix is not positive definite even with the largest jitter.
pub fn try_fit_gp(
    x_train: Vec<Vec<f64>>,
    y_train: Vec<f64>,
    sigma_f: f64,
    length_scale: f64,
    noise_std: f64,
) -> PhysicsResult<GaussianProcessSurrogate> {
    validate_gp_inputs(&x_train, &y_train, sigma_f, length_scale, noise_std)?;
    let n = x_train.len();
    if n == 0 {
        return Ok(GaussianProcessSurrogate::prior(
            sigma_f,
            length_scale,
            noise_std,
        ));
    }
    let mut covariance = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let value = rbf_kernel(&x_train[i], &x_train[j], sigma_f, length_scale);
            covariance[(i, j)] = value;
            covariance[(j, i)] = value;
        }
        covariance[(i, i)] += noise_std * noise_std;
    }
    let targets = DVector::from_column_slice(&y_train);
    let signal_variance = sigma_f * sigma_f;
    let jitters = std::iter::once(0.0).chain(JITTER_LADDER.iter().map(|rel| rel * signal_variance));
    for jitter in jitters {
        let mut attempt = covariance.clone();
        if jitter > 0.0 {
            for i in 0..n {
                attempt[(i, i)] += jitter;
            }
        }
        let Some(cholesky) = attempt.cholesky() else {
            continue;
        };
        let alpha = cholesky.solve(&targets);
        let chol_l = cholesky.l();
        if alpha.iter().any(|a| !a.is_finite()) || chol_l.iter().any(|v| !v.is_finite()) {
            continue;
        }
        return Ok(GaussianProcessSurrogate {
            x_train,
            y_train,
            alpha: alpha.iter().copied().collect(),
            sigma_f,
            length_scale,
            noise_std,
            chol_l,
            jitter,
        });
    }
    Err(PhysicsError::InvalidInput {
        field: "x_train".to_string(),
        reason: "the covariance matrix K + noise_std² I is not positive definite even with \
                 diagonal jitter; training inputs (nearly) coincide while noise_std is too small"
            .to_string(),
    })
}

/// Fit a Gaussian Process surrogate using exact GP regression.
///
/// Same as [`try_fit_gp`], but infallible: if the data cannot be fitted (see
/// the errors of [`try_fit_gp`]) the result is the *prior* surrogate, which has
/// no training data and predicts mean `0` and variance `σ_f²` everywhere.
/// Use [`try_fit_gp`] to find out why a fit failed.
pub fn fit_gp(
    x_train: Vec<Vec<f64>>,
    y_train: Vec<f64>,
    sigma_f: f64,
    length_scale: f64,
    noise_std: f64,
) -> GaussianProcessSurrogate {
    try_fit_gp(x_train, y_train, sigma_f, length_scale, noise_std)
        .unwrap_or_else(|_| GaussianProcessSurrogate::prior(sigma_f, length_scale, noise_std))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic uniform generator in [0, 1).
    struct Uniform(u64);

    impl Uniform {
        fn next(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn quadratic_data() -> (Vec<Vec<f64>>, Vec<f64>) {
        let x: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = x.iter().map(|v| v[0] * v[0]).collect();
        (x, y)
    }

    fn invalid_field(error: &PhysicsError) -> &str {
        match error {
            PhysicsError::InvalidInput { field, .. } => field.as_str(),
            other => panic!("expected InvalidInput, got {other:?}"),
        }
    }

    #[test]
    fn variance_collapses_at_training_points_and_the_mean_interpolates() {
        let (x, y) = quadratic_data();
        let gp = try_fit_gp(x.clone(), y.clone(), 2.0, 1.5, 1e-4).expect("well-posed fit");
        assert_eq!(gp.jitter, 0.0);
        for (xi, yi) in x.iter().zip(y.iter()) {
            let (mean, variance) = gp.predict(xi);
            assert!((mean - yi).abs() < 1e-3, "mean {mean} vs {yi}");
            assert!(variance < 1e-6, "variance {variance} at a training point");
        }
    }

    #[test]
    fn variance_returns_to_the_prior_far_from_the_data() {
        let (x, y) = quadratic_data();
        let gp = try_fit_gp(x, y, 2.0, 1.5, 1e-4).expect("well-posed fit");
        // 20 length scales away from the last training input.
        let (mean, variance) = gp.predict(&[4.0 + 20.0 * 1.5]);
        assert!(mean.abs() < 1e-9, "mean {mean}");
        assert!((variance - 4.0).abs() < 1e-9, "variance {variance}");
    }

    #[test]
    fn variance_grows_between_training_points() {
        let (x, y) = quadratic_data();
        let gp = try_fit_gp(x, y, 2.0, 1.5, 1e-4).expect("well-posed fit");
        let (_, at_data) = gp.predict(&[2.0]);
        let (_, between) = gp.predict(&[2.5]);
        let (_, outside) = gp.predict(&[7.0]);
        assert!(at_data < between && between < outside);
        assert!(outside <= 4.0 + 1e-12);
    }

    #[test]
    fn single_observation_matches_the_closed_form_posterior() {
        // n = 1: K + σ_n² = σ_f² + σ_n², k* = σ_f² exp(-d² / 2ℓ²). The old
        // "diagonal approximation" ignored σ_n and gave a different variance.
        let (sigma_f, length_scale, noise_std) = (1.5, 2.0, 1.0);
        let gp = try_fit_gp(vec![vec![0.0]], vec![3.0], sigma_f, length_scale, noise_std)
            .expect("well-posed fit");
        let distance = 1.0_f64;
        let k_star =
            sigma_f * sigma_f * (-distance * distance / (2.0 * length_scale * length_scale)).exp();
        let denominator = sigma_f * sigma_f + noise_std * noise_std;
        let (mean, variance) = gp.predict(&[distance]);
        assert!((mean - k_star * 3.0 / denominator).abs() < 1e-12);
        let expected = sigma_f * sigma_f - k_star * k_star / denominator;
        assert!(
            (variance - expected).abs() < 1e-12,
            "{variance} vs {expected}"
        );
    }

    #[test]
    fn prediction_matches_the_explicit_matrix_inverse_formula() {
        let mut rng = Uniform(0xdead_beef_1234_5678);
        let n = 6;
        let x: Vec<Vec<f64>> = (0..n)
            .map(|_| vec![4.0 * rng.next(), 4.0 * rng.next()])
            .collect();
        let y: Vec<f64> = x.iter().map(|p| (p[0] * 1.3).sin() + 0.5 * p[1]).collect();
        let (sigma_f, length_scale, noise_std) = (1.2, 1.1, 0.2);
        let gp = try_fit_gp(x.clone(), y.clone(), sigma_f, length_scale, noise_std)
            .expect("well-posed fit");

        let mut k = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                k[(i, j)] = rbf_kernel(&x[i], &x[j], sigma_f, length_scale);
            }
            k[(i, i)] += noise_std * noise_std;
        }
        let k_inverse = k.clone().try_inverse().expect("K is positive definite");
        let targets = DVector::from_column_slice(&y);
        for _ in 0..20 {
            let probe = [4.0 * rng.next(), 4.0 * rng.next()];
            let k_star = DVector::from_iterator(
                n,
                x.iter()
                    .map(|xi| rbf_kernel(xi, &probe, sigma_f, length_scale)),
            );
            let expected_mean = (k_star.transpose() * &k_inverse * &targets)[(0, 0)];
            let expected_variance = rbf_kernel(&probe, &probe, sigma_f, length_scale)
                - (k_star.transpose() * &k_inverse * &k_star)[(0, 0)];
            let (mean, variance) = gp.predict(&probe);
            assert!(
                (mean - expected_mean).abs() < 1e-9,
                "{mean} vs {expected_mean}"
            );
            assert!(
                (variance - expected_variance).abs() < 1e-9,
                "{variance} vs {expected_variance}"
            );
            assert!((0.0..=sigma_f * sigma_f + 1e-12).contains(&variance));
        }
    }

    #[test]
    fn duplicated_inputs_without_noise_are_rescued_by_jitter() {
        let x = vec![vec![0.0], vec![0.0], vec![1.0], vec![2.0]];
        let y = vec![1.0, 1.0, 2.0, 0.5];
        let gp = try_fit_gp(x, y, 1.0, 1.0, 0.0).expect("jitter must make the fit possible");
        assert!(gp.jitter > 0.0 && gp.jitter <= 1e-6, "jitter {}", gp.jitter);
        let (mean, variance) = gp.predict(&[0.0]);
        assert!(mean.is_finite() && (mean - 1.0).abs() < 1e-3, "mean {mean}");
        assert!((0.0..1e-3).contains(&variance), "variance {variance}");
    }

    #[test]
    fn invalid_inputs_are_reported_by_try_fit_gp() {
        let (x, y) = quadratic_data();
        let fit = |x: Vec<Vec<f64>>, y: Vec<f64>, sf: f64, ls: f64, ns: f64| {
            try_fit_gp(x, y, sf, ls, ns).expect_err("must be rejected")
        };
        assert_eq!(
            invalid_field(&fit(x.clone(), y[..4].to_vec(), 1.0, 1.0, 0.1)),
            "y_train"
        );
        assert_eq!(
            invalid_field(&fit(x.clone(), y.clone(), 0.0, 1.0, 0.1)),
            "sigma_f"
        );
        assert_eq!(
            invalid_field(&fit(x.clone(), y.clone(), f64::NAN, 1.0, 0.1)),
            "sigma_f"
        );
        assert_eq!(
            invalid_field(&fit(x.clone(), y.clone(), 1.0, -1.0, 0.1)),
            "length_scale"
        );
        assert_eq!(
            invalid_field(&fit(x.clone(), y.clone(), 1.0, 1.0, -0.1)),
            "noise_std"
        );
        assert_eq!(
            invalid_field(&fit(x.clone(), y.clone(), 1.0, 1.0, f64::NAN)),
            "noise_std"
        );
        let mut bad_y = y.clone();
        bad_y[2] = f64::INFINITY;
        assert_eq!(
            invalid_field(&fit(x.clone(), bad_y, 1.0, 1.0, 0.1)),
            "y_train"
        );
        let mut bad_x = x.clone();
        bad_x[1][0] = f64::NAN;
        assert_eq!(
            invalid_field(&fit(bad_x, y.clone(), 1.0, 1.0, 0.1)),
            "x_train"
        );
        let mut ragged = x;
        ragged[3].push(0.0);
        assert_eq!(invalid_field(&fit(ragged, y, 1.0, 1.0, 0.1)), "x_train");
    }

    #[test]
    fn a_covariance_matrix_that_cannot_be_factorised_is_an_error() {
        // σ_f² underflows to zero, so K is the zero matrix and no jitter helps.
        let (x, y) = quadratic_data();
        let error = try_fit_gp(x, y, 1e-200, 1.0, 0.0).expect_err("K is not positive definite");
        assert_eq!(invalid_field(&error), "x_train");
    }

    #[test]
    fn fit_gp_degrades_to_the_prior_when_the_data_cannot_be_fitted() {
        let (x, y) = quadratic_data();
        let gp = fit_gp(x, y, 2.0, 1.5, f64::NAN);
        assert!(gp.x_train.is_empty() && gp.alpha.is_empty());
        let (mean, variance) = gp.predict(&[1.0]);
        assert_eq!(mean, 0.0);
        assert!((variance - 4.0).abs() < 1e-12);
    }

    #[test]
    fn empty_training_set_predicts_the_prior() {
        let gp =
            try_fit_gp(Vec::new(), Vec::new(), 3.0, 0.5, 0.1).expect("an empty fit is the prior");
        let (mean, variance) = gp.predict(&[0.3, 0.7]);
        assert_eq!(mean, 0.0);
        assert!((variance - 9.0).abs() < 1e-12);
    }

    #[test]
    fn fit_gp_matches_try_fit_gp_on_valid_data() {
        let (x, y) = quadratic_data();
        let direct = try_fit_gp(x.clone(), y.clone(), 2.0, 1.5, 0.01).expect("well-posed fit");
        let wrapped = fit_gp(x, y, 2.0, 1.5, 0.01);
        assert_eq!(direct.alpha, wrapped.alpha);
        assert_eq!(direct.chol_l, wrapped.chol_l);
        assert_eq!(direct.predict(&[2.2]), wrapped.predict(&[2.2]));
    }
}
