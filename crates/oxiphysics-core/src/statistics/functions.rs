//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::collections::VecDeque;

use super::types::{NormalDistribution, PcaResult, StatRng};

/// Returns the arithmetic mean of `data`.  Returns `0.0` for an empty slice.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}
/// Returns the unbiased sample variance of `data` (divides by n − 1).
///
/// Returns `0.0` for slices with fewer than 2 elements.
pub fn variance(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0)
}
/// Returns the unbiased sample standard deviation of `data`.
///
/// Returns `0.0` for slices with fewer than 2 elements.
pub fn std_dev(data: &[f64]) -> f64 {
    variance(data).sqrt()
}
/// Returns the median of `data`, working on a sorted copy.
///
/// Returns `0.0` for an empty slice.  For even-length slices the two middle
/// values are averaged.
pub fn median(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}
/// Returns the `p`-th percentile of `data` using linear interpolation.
///
/// `p` should be in \[0, 100\].  Works on a sorted copy of `data`.
/// Returns `0.0` for an empty slice.
pub fn percentile(data: &[f64], p: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = (p / 100.0) * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = rank - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}
/// Returns the sample covariance of `x` and `y` (divides by n − 1).
///
/// Panics if the slices have different lengths.
/// Returns `0.0` for fewer than 2 paired values.
pub fn covariance(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        y.len(),
        "covariance: x and y must have the same length"
    );
    let n = x.len();
    if n < 2 {
        return 0.0;
    }
    let mx = mean(x);
    let my = mean(y);
    x.iter()
        .zip(y.iter())
        .map(|(xi, yi)| (xi - mx) * (yi - my))
        .sum::<f64>()
        / (n as f64 - 1.0)
}
/// Returns the Pearson correlation coefficient of `x` and `y`.
///
/// Returns `0.0` when either standard deviation is zero or there are fewer
/// than 2 paired values.
pub fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let sx = std_dev(x);
    let sy = std_dev(y);
    if sx == 0.0 || sy == 0.0 {
        return 0.0;
    }
    covariance(x, y) / (sx * sy)
}
/// Returns the sample skewness of `data`.
///
/// Uses the adjusted Fisher–Pearson standardised moment.
/// Returns `0.0` for fewer than 3 elements or zero standard deviation.
pub fn skewness(data: &[f64]) -> f64 {
    let n = data.len();
    if n < 3 {
        return 0.0;
    }
    let m = mean(data);
    let s = std_dev(data);
    if s == 0.0 {
        return 0.0;
    }
    let nf = n as f64;
    let third_moment = data.iter().map(|x| ((x - m) / s).powi(3)).sum::<f64>() / nf;
    third_moment * (nf * (nf - 1.0)).sqrt() / (nf - 2.0)
}
/// Returns the excess kurtosis of `data` (kurtosis minus 3).
///
/// Returns `0.0` for fewer than 4 elements or zero variance.
pub fn kurtosis(data: &[f64]) -> f64 {
    let n = data.len();
    if n < 4 {
        return 0.0;
    }
    let m = mean(data);
    let s = std_dev(data);
    if s == 0.0 {
        return 0.0;
    }
    let nf = n as f64;
    let fourth = data.iter().map(|x| ((x - m) / s).powi(4)).sum::<f64>() / nf;
    fourth - 3.0
}
/// Bins `data` into `n_bins` equal-width bins.
///
/// Returns a `Vec` of `(bin_low, bin_high, count)` tuples.
/// Returns an empty vector when `data` is empty or `n_bins` is zero.
pub fn histogram(data: &[f64], n_bins: usize) -> Vec<(f64, f64, usize)> {
    if data.is_empty() || n_bins == 0 {
        return vec![];
    }
    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let width = if (max - min).abs() < f64::EPSILON {
        1.0
    } else {
        (max - min) / n_bins as f64
    };
    let mut counts = vec![0usize; n_bins];
    for &v in data {
        let idx = ((v - min) / width) as usize;
        let idx = idx.min(n_bins - 1);
        counts[idx] += 1;
    }
    (0..n_bins)
        .map(|i| {
            let lo = min + i as f64 * width;
            let hi = lo + width;
            (lo, hi, counts[i])
        })
        .collect()
}
/// One-sample t-test against a hypothesised mean `mu0`.
///
/// Returns `(t_statistic, approximate_p_value)`.
/// The p-value is a two-tailed approximation based on a large-sample normal
/// approximation (valid when n is not too small).
///
/// Returns `(0.0, 1.0)` for fewer than 2 data points or zero variance.
pub fn t_test_one_sample(data: &[f64], mu0: f64) -> (f64, f64) {
    let n = data.len();
    if n < 2 {
        return (0.0, 1.0);
    }
    let m = mean(data);
    let s = std_dev(data);
    if s == 0.0 {
        return (0.0, 1.0);
    }
    let t = (m - mu0) / (s / (n as f64).sqrt());
    let p = 2.0 * (1.0 - NormalDistribution::new(0.0, 1.0).cdf(t.abs()));
    (t, p)
}
/// Pearson's χ² goodness-of-fit statistic.
///
/// Returns Σ (O_i − E_i)² / E_i.
///
/// Panics if the slices have different lengths.
/// Bins with expected count zero are skipped.
pub fn chi_squared_test(observed: &[f64], expected: &[f64]) -> f64 {
    assert_eq!(
        observed.len(),
        expected.len(),
        "chi_squared_test: slices must have the same length"
    );
    observed
        .iter()
        .zip(expected.iter())
        .filter(|&(_, &e)| e > 0.0)
        .map(|(&o, &e)| {
            let diff = o - e;
            diff * diff / e
        })
        .sum()
}
/// Kolmogorov–Smirnov statistic D = sup |F_n(x) − F(x)|.
///
/// `cdf` should be the theoretical CDF.  The data are sorted internally.
pub fn ks_statistic(data: &[f64], cdf: impl Fn(f64) -> f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len() as f64;
    sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let fn_hi = (i as f64 + 1.0) / n;
            let fn_lo = i as f64 / n;
            let fx = cdf(x);
            (fn_hi - fx).abs().max((fx - fn_lo).abs())
        })
        .fold(0.0_f64, f64::max)
}
/// Returns the Pearson correlation coefficient (alias for [`correlation`]).
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    correlation(x, y)
}
/// Returns a closure that computes a running mean over the last `window` values.
///
/// Each call is O(1) using a ring buffer internally.
/// Returns `0.0` for any call when `window` is 0.
pub fn running_average(window: usize) -> impl FnMut(f64) -> f64 {
    let mut buf: VecDeque<f64> = VecDeque::with_capacity(window.max(1));
    let mut sum = 0.0_f64;
    move |value: f64| -> f64 {
        if window == 0 {
            return 0.0;
        }
        if buf.len() == window
            && let Some(old) = buf.pop_front()
        {
            sum -= old;
        }
        buf.push_back(value);
        sum += value;
        sum / buf.len() as f64
    }
}
/// Block-averaging estimator.
///
/// Divides `data` into non-overlapping blocks of `block_size` values and
/// returns `(grand_mean, standard_error_of_mean)`.
/// Returns `(0.0, 0.0)` when data is empty or `block_size` is zero.
pub fn block_average(data: &[f64], block_size: usize) -> (f64, f64) {
    if data.is_empty() || block_size == 0 {
        return (0.0, 0.0);
    }
    let block_means: Vec<f64> = data.chunks_exact(block_size).map(mean).collect();
    let n_blocks = block_means.len();
    if n_blocks < 2 {
        return (mean(data), 0.0);
    }
    let grand_mean = mean(&block_means);
    let sem = std_dev(&block_means) / (n_blocks as f64).sqrt();
    (grand_mean, sem)
}
/// Integrated autocorrelation time τ of `data`.
///
/// Uses τ = 1 + 2 Σ C(t)/C(0) with standard windowing (stops at first
/// non-positive lag).  Returns `1.0` when variance is zero or `max_lag` is 0.
pub fn autocorrelation_time(data: &[f64], max_lag: usize) -> f64 {
    let n = data.len();
    if n < 2 || max_lag == 0 {
        return 1.0;
    }
    let m = mean(data);
    let c0: f64 = data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n as f64;
    if c0 == 0.0 {
        return 1.0;
    }
    let mut tau = 1.0_f64;
    for t in 1..=max_lag.min(n - 1) {
        let ct: f64 = data[..n - t]
            .iter()
            .zip(data[t..].iter())
            .map(|(a, b)| (a - m) * (b - m))
            .sum::<f64>()
            / n as f64;
        let rho = ct / c0;
        if rho <= 0.0 {
            break;
        }
        tau += 2.0 * rho;
    }
    tau
}
/// Most probable speed in the Maxwell–Boltzmann distribution: v_p = √(2 k_B T / m).
pub fn maxwell_boltzmann_speed(mass: f64, temp: f64, kb: f64) -> f64 {
    (2.0 * kb * temp / mass).sqrt()
}
/// Boltzmann factor exp(−E / (k_B T)).
pub fn boltzmann_factor(energy: f64, temp: f64, kb: f64) -> f64 {
    (-energy / (kb * temp)).exp()
}
/// Canonical partition function Z = Σ exp(−E_i / (k_B T)).
pub fn partition_function(energies: &[f64], temp: f64, kb: f64) -> f64 {
    energies
        .iter()
        .map(|&e| boltzmann_factor(e, temp, kb))
        .sum()
}
/// Helmholtz free energy F = −k_B T ln(Z).
pub fn free_energy_from_partition(z: f64, temp: f64, kb: f64) -> f64 {
    -kb * temp * z.ln()
}
/// Compute the Pearson correlation matrix for a data matrix.
///
/// `data[i]` is the i-th observation (a row), each of length `d`.
/// Returns a `d × d` correlation matrix.
pub fn correlation_matrix(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if data.is_empty() {
        return Vec::new();
    }
    let d = data[0].len();
    if d == 0 {
        return Vec::new();
    }
    let col = |j: usize| -> Vec<f64> { data.iter().map(|row| row[j]).collect() };
    let mut mat = vec![vec![0.0f64; d]; d];
    for i in 0..d {
        for j in i..d {
            let ci = col(i);
            let cj = col(j);
            let r = if i == j { 1.0 } else { correlation(&ci, &cj) };
            mat[i][j] = r;
            mat[j][i] = r;
        }
    }
    mat
}
/// Covariance matrix for a data matrix.
///
/// `data[i]` is the i-th observation.  Returns a `d × d` matrix.
pub fn covariance_matrix(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if data.is_empty() {
        return Vec::new();
    }
    let d = data[0].len();
    if d == 0 {
        return Vec::new();
    }
    let col = |j: usize| -> Vec<f64> { data.iter().map(|row| row[j]).collect() };
    let mut mat = vec![vec![0.0f64; d]; d];
    for i in 0..d {
        let ci = col(i);
        for j in i..d {
            let cj = col(j);
            let cov = covariance(&ci, &cj);
            mat[i][j] = cov;
            mat[j][i] = cov;
        }
    }
    mat
}
/// Compute PCA using power iteration on the covariance matrix.
///
/// Returns up to `n_components` principal components.
/// `data[i]` is the i-th observation of length `d`.
pub fn pca(data: &[Vec<f64>], n_components: usize) -> Option<PcaResult> {
    if data.is_empty() {
        return None;
    }
    let n = data.len();
    let d = data[0].len();
    if d == 0 || n < 2 {
        return None;
    }
    let n_comp = n_components.min(d);
    let mean_vec: Vec<f64> = (0..d)
        .map(|j| data.iter().map(|row| row[j]).sum::<f64>() / n as f64)
        .collect();
    let centred: Vec<Vec<f64>> = data
        .iter()
        .map(|row| row.iter().zip(&mean_vec).map(|(x, m)| x - m).collect())
        .collect();
    let cov: Vec<Vec<f64>> = (0..d)
        .map(|i| {
            (0..d)
                .map(|j| centred.iter().map(|row| row[i] * row[j]).sum::<f64>() / (n as f64 - 1.0))
                .collect()
        })
        .collect();
    let mut components = Vec::with_capacity(n_comp);
    let mut variances = Vec::with_capacity(n_comp);
    let mut deflated = cov.clone();
    let mat_vec = |m: &Vec<Vec<f64>>, v: &Vec<f64>| -> Vec<f64> {
        (0..d)
            .map(|i| m[i].iter().zip(v.iter()).map(|(a, b)| a * b).sum::<f64>())
            .collect()
    };
    let norm_vec = |v: &Vec<f64>| -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() };
    for _ in 0..n_comp {
        let mut q: Vec<f64> = (0..d).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
        for _ in 0..200 {
            let aq = mat_vec(&deflated, &q);
            let nrm = norm_vec(&aq);
            if nrm < 1e-30 {
                break;
            }
            q = aq.iter().map(|x| x / nrm).collect();
        }
        let aq = mat_vec(&deflated, &q);
        let eigenvalue: f64 = q.iter().zip(aq.iter()).map(|(qi, aqj)| qi * aqj).sum();
        for i in 0..d {
            for j in 0..d {
                deflated[i][j] -= eigenvalue * q[i] * q[j];
            }
        }
        components.push(q);
        variances.push(eigenvalue);
    }
    Some(PcaResult {
        mean: mean_vec,
        components,
        explained_variance: variances,
    })
}
/// Project a centred data vector onto the principal components.
///
/// Returns a vector of `n_components` scores.
pub fn pca_transform(x: &[f64], result: &PcaResult) -> Vec<f64> {
    let centred: Vec<f64> = x.iter().zip(&result.mean).map(|(xi, mi)| xi - mi).collect();
    result
        .components
        .iter()
        .map(|pc| pc.iter().zip(&centred).map(|(pci, ci)| pci * ci).sum())
        .collect()
}
/// Autocorrelation function (ACF) up to `max_lag`.
///
/// Returns a vector of length `max_lag + 1` where index 0 is always 1.0.
pub fn acf(data: &[f64], max_lag: usize) -> Vec<f64> {
    let n = data.len();
    if n < 2 {
        return vec![1.0];
    }
    let m = mean(data);
    let c0: f64 = data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n as f64;
    if c0 < f64::EPSILON {
        return vec![1.0; max_lag + 1];
    }
    (0..=max_lag.min(n - 1))
        .map(|lag| {
            if lag == 0 {
                1.0
            } else {
                let ct: f64 = data[..n - lag]
                    .iter()
                    .zip(data[lag..].iter())
                    .map(|(a, b)| (a - m) * (b - m))
                    .sum::<f64>()
                    / n as f64;
                ct / c0
            }
        })
        .collect()
}
/// Partial autocorrelation function (PACF) using the Yule-Walker equations.
///
/// Returns a vector of length `max_lag + 1` where index 0 is always 1.0 and
/// index 1 equals `acf(data, 1)[1]`.
pub fn pacf(data: &[f64], max_lag: usize) -> Vec<f64> {
    let rho = acf(data, max_lag);
    let m = rho.len();
    if m < 2 {
        return vec![1.0];
    }
    let mut phi = vec![1.0_f64];
    for k in 1..m {
        let prev_order = phi.len() - 1;
        if prev_order == 0 {
            phi.push(rho[k]);
            continue;
        }
        let prev_phi: Vec<f64> = phi[1..].to_vec();
        let num: f64 = rho[k]
            - prev_phi
                .iter()
                .zip(1..=prev_order)
                .map(|(p, j)| p * rho[k - j])
                .sum::<f64>();
        let den: f64 = 1.0
            - prev_phi
                .iter()
                .zip(1..=prev_order)
                .map(|(p, j)| p * rho[j])
                .sum::<f64>();
        let phi_kk = if den.abs() < 1e-30 { 0.0 } else { num / den };
        phi.push(phi_kk);
    }
    phi
}
/// Compute a bootstrap confidence interval for a statistic.
///
/// - `data`: original sample.
/// - `statistic`: function mapping a sample to a scalar.
/// - `n_boot`: number of bootstrap resamples.
/// - `alpha`: significance level (e.g. 0.05 for 95% CI).
/// - `rng`: random number generator.
///
/// Returns `(lower, upper)` percentile bootstrap CI.
pub fn bootstrap_ci(
    data: &[f64],
    statistic: impl Fn(&[f64]) -> f64,
    n_boot: usize,
    alpha: f64,
    rng: &mut StatRng,
) -> (f64, f64) {
    let n = data.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mut boot_stats: Vec<f64> = Vec::with_capacity(n_boot);
    for _ in 0..n_boot {
        let resample: Vec<f64> = (0..n)
            .map(|_| {
                let idx = (rng.next_f64() * n as f64) as usize;
                data[idx.min(n - 1)]
            })
            .collect();
        boot_stats.push(statistic(&resample));
    }
    boot_stats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let lo_idx = ((alpha / 2.0) * n_boot as f64) as usize;
    let hi_idx = ((1.0 - alpha / 2.0) * n_boot as f64) as usize;
    (
        boot_stats[lo_idx.min(n_boot.saturating_sub(1))],
        boot_stats[hi_idx.min(n_boot.saturating_sub(1))],
    )
}
/// Bootstrap standard error of a statistic.
pub fn bootstrap_se(
    data: &[f64],
    statistic: impl Fn(&[f64]) -> f64,
    n_boot: usize,
    rng: &mut StatRng,
) -> f64 {
    let n = data.len();
    if n == 0 {
        return 0.0;
    }
    let boot_stats: Vec<f64> = (0..n_boot)
        .map(|_| {
            let resample: Vec<f64> = (0..n)
                .map(|_| {
                    let idx = (rng.next_f64() * n as f64) as usize;
                    data[idx.min(n - 1)]
                })
                .collect();
            statistic(&resample)
        })
        .collect();
    std_dev(&boot_stats)
}
/// Two-sample independent t-test.
///
/// Returns `(t_statistic, approximate_p_value)`.
/// Uses Welch's t-test (unequal variances).
/// Returns `(0.0, 1.0)` when either sample has fewer than 2 elements.
pub fn t_test_two_sample(x: &[f64], y: &[f64]) -> (f64, f64) {
    if x.len() < 2 || y.len() < 2 {
        return (0.0, 1.0);
    }
    let mx = mean(x);
    let my = mean(y);
    let sx2 = variance(x);
    let sy2 = variance(y);
    let nx = x.len() as f64;
    let ny = y.len() as f64;
    let se = (sx2 / nx + sy2 / ny).sqrt();
    if se < f64::EPSILON {
        return (0.0, 1.0);
    }
    let t = (mx - my) / se;
    let p = 2.0 * (1.0 - NormalDistribution::new(0.0, 1.0).cdf(t.abs()));
    (t, p)
}
/// Mann-Whitney U statistic (non-parametric two-sample test).
///
/// Returns the U statistic (for sample x against sample y).
/// Smaller U indicates x tends to have smaller values than y.
pub fn mann_whitney_u(x: &[f64], y: &[f64]) -> f64 {
    let nx = x.len() as f64;
    let ny = y.len() as f64;
    if x.is_empty() || y.is_empty() {
        return 0.0;
    }
    let u: f64 = x
        .iter()
        .map(|&xi| {
            y.iter()
                .map(|&yj| {
                    if xi > yj {
                        1.0
                    } else if (xi - yj).abs() < f64::EPSILON {
                        0.5
                    } else {
                        0.0
                    }
                })
                .sum::<f64>()
        })
        .sum();
    let _ = (nx, ny);
    u
}
/// Wilcoxon signed-rank test statistic for paired data.
///
/// Returns the test statistic W+ (sum of positive ranks).
/// Ties and zero differences are handled by ignoring zero differences.
pub fn wilcoxon_signed_rank(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        y.len(),
        "wilcoxon_signed_rank: x and y must have same length"
    );
    let diffs: Vec<f64> = x.iter().zip(y.iter()).map(|(xi, yi)| xi - yi).collect();
    let nonzero: Vec<f64> = diffs
        .into_iter()
        .filter(|&d| d.abs() > f64::EPSILON)
        .collect();
    if nonzero.is_empty() {
        return 0.0;
    }
    let mut indexed: Vec<(f64, f64)> = nonzero.iter().map(|&d| (d.abs(), d.signum())).collect();
    indexed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let w_plus: f64 = indexed
        .iter()
        .enumerate()
        .filter(|(_, (_, sign))| *sign > 0.0)
        .map(|(rank, _)| rank as f64 + 1.0)
        .sum();
    w_plus
}
/// Shapiro-Wilk normality approximation (W statistic only, no p-value).
///
/// Returns a value close to 1.0 for normally distributed data.
/// This is a simplified approximation suitable for quick checks.
pub fn shapiro_wilk_w(data: &[f64]) -> f64 {
    let n = data.len();
    if n < 3 {
        return 1.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = mean(&sorted);
    let s2: f64 = sorted.iter().map(|x| (x - m).powi(2)).sum::<f64>();
    if s2 < f64::EPSILON {
        return 1.0;
    }
    let half = n / 2;
    let mut b_sum = 0.0f64;
    for i in 0..half {
        let rank = n - 1 - i;
        let p = (rank as f64 + 1.0 - 0.375) / (n as f64 + 0.25);
        let a_i = normal_quantile_approx(p);
        b_sum += a_i * (sorted[rank] - sorted[i]);
    }
    (b_sum * b_sum) / s2
}
/// Approximate normal quantile (probit) function for p in (0, 1).
pub(super) fn normal_quantile_approx(p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return 0.0;
    }
    let pp = if p < 0.5 { p } else { 1.0 - p };
    let t = (-2.0 * pp.ln()).sqrt();
    let c0 = 2.515517_f64;
    let c1 = 0.802853_f64;
    let c2 = 0.010328_f64;
    let d1 = 1.432788_f64;
    let d2 = 0.189269_f64;
    let d3 = 0.001308_f64;
    let x = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);
    if p < 0.5 { -x } else { x }
}
/// Approximation of erf(x) using Abramowitz & Stegun formula 7.1.26.
pub(super) fn erf_approx(x: f64) -> f64 {
    pub(super) const P: f64 = 0.3275911;
    pub(super) const A1: f64 = 0.254829592;
    pub(super) const A2: f64 = -0.284496736;
    pub(super) const A3: f64 = 1.421413741;
    pub(super) const A4: f64 = -1.453152027;
    pub(super) const A5: f64 = 1.061405429;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + P * x);
    let poly = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5))));
    sign * (1.0 - poly * (-x * x).exp())
}
/// Natural log of the gamma function via Stirling's approximation.
///
/// Accurate for z > ~0.5; uses the recurrence for small z.
pub(super) fn log_gamma_stirling(z: f64) -> f64 {
    if z < 7.0 {
        return log_gamma_stirling(z + 1.0) - z.ln();
    }
    let z = z - 1.0;
    0.5 * (std::f64::consts::TAU / z).ln()
        + z * ((z + 1.0 / (12.0 * z) - 1.0 / (360.0 * z.powi(3))) / std::f64::consts::E).ln()
}
/// Natural log of k! computed via log-gamma (Stirling) for large k,
/// or a direct sum for small k.
pub(super) fn log_factorial(k: u64) -> f64 {
    if k == 0 || k == 1 {
        return 0.0;
    }
    if k <= 20 {
        (1..=k).map(|i| (i as f64).ln()).sum()
    } else {
        log_gamma_stirling(k as f64 + 1.0)
    }
}
/// Returns the lower quartile (Q1), median (Q2), and upper quartile (Q3).
///
/// Uses linear interpolation (same method as `percentile`).
/// Returns `(0.0, 0.0, 0.0)` for empty input.
pub fn quartiles(data: &[f64]) -> (f64, f64, f64) {
    if data.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    (
        percentile(data, 25.0),
        percentile(data, 50.0),
        percentile(data, 75.0),
    )
}
/// Returns the inter-quartile range (IQR = Q3 − Q1).
///
/// Returns `0.0` for empty input.
pub fn iqr(data: &[f64]) -> f64 {
    let (q1, _, q3) = quartiles(data);
    q3 - q1
}
/// Spearman rank correlation coefficient between `x` and `y`.
///
/// Ranks ties using the average-rank convention.
/// Returns `0.0` when fewer than 2 paired values or when either rank series
/// has zero variance.
pub fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(
        x.len(),
        y.len(),
        "spearman_correlation: x and y must be same length"
    );
    let n = x.len();
    if n < 2 {
        return 0.0;
    }
    let rank = |data: &[f64]| -> Vec<f64> {
        let mut indexed: Vec<(usize, f64)> = data.iter().cloned().enumerate().collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut ranks = vec![0.0_f64; data.len()];
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && (indexed[j + 1].1 - indexed[i].1).abs() < f64::EPSILON {
                j += 1;
            }
            let avg_rank = (i + j) as f64 / 2.0 + 1.0;
            for k in i..=j {
                ranks[indexed[k].0] = avg_rank;
            }
            i = j + 1;
        }
        ranks
    };
    let rx = rank(x);
    let ry = rank(y);
    correlation(&rx, &ry)
}
/// Simple ordinary least-squares linear regression.
///
/// Fits the model `y = slope * x + intercept` and returns
/// `(slope, intercept, r_squared)`.
/// Returns `(0.0, 0.0, 0.0)` for fewer than 2 points.
pub fn linear_regression(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    assert_eq!(
        x.len(),
        y.len(),
        "linear_regression: x and y must be same length"
    );
    let n = x.len();
    if n < 2 {
        return (0.0, 0.0, 0.0);
    }
    let mx = mean(x);
    let my = mean(y);
    let sxx: f64 = x.iter().map(|xi| (xi - mx).powi(2)).sum();
    let sxy: f64 = x
        .iter()
        .zip(y.iter())
        .map(|(xi, yi)| (xi - mx) * (yi - my))
        .sum();
    if sxx.abs() < f64::EPSILON {
        return (0.0, my, 0.0);
    }
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    let ss_res: f64 = x
        .iter()
        .zip(y.iter())
        .map(|(xi, yi)| (yi - (slope * xi + intercept)).powi(2))
        .sum();
    let ss_tot: f64 = y.iter().map(|yi| (yi - my).powi(2)).sum();
    let r_squared = if ss_tot < f64::EPSILON {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    };
    (slope, intercept, r_squared)
}
/// One-sample Kolmogorov-Smirnov test statistic against a theoretical CDF.
///
/// `cdf` must be a monotone non-decreasing function mapping `f64 → f64` in
/// \[0, 1\].  Returns the KS statistic D (maximum absolute deviation).
/// Already defined as `ks_statistic` — this is an alias with different name
/// for symmetry with the two-sample version.
pub fn ks_test_one_sample(data: &[f64], cdf: impl Fn(f64) -> f64) -> f64 {
    ks_statistic(data, cdf)
}
/// Two-sample Kolmogorov-Smirnov statistic between empirical distributions of
/// `x` and `y`.
///
/// Returns the KS distance D = max |F_x(t) − F_y(t)| over all t.
pub fn ks_test_two_sample(x: &[f64], y: &[f64]) -> f64 {
    if x.is_empty() || y.is_empty() {
        return 0.0;
    }
    let mut xs = x.to_vec();
    let mut ys = y.to_vec();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let nx = xs.len() as f64;
    let ny = ys.len() as f64;
    let mut all: Vec<f64> = xs.iter().chain(ys.iter()).cloned().collect();
    all.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    all.dedup_by(|a, b| (*a - *b).abs() < f64::EPSILON);
    let mut d_max = 0.0_f64;
    for &t in &all {
        let fx = xs.partition_point(|&v| v <= t) as f64 / nx;
        let fy = ys.partition_point(|&v| v <= t) as f64 / ny;
        let d = (fx - fy).abs();
        if d > d_max {
            d_max = d;
        }
    }
    d_max
}
/// Chi-squared goodness-of-fit test — returns the χ² statistic.
///
/// `observed` and `expected` must have the same length.  Cells with expected
/// count < 1e-10 are skipped (avoids division by zero).
///
/// Note: `chi_squared_test` already exists — this variant accepts integer
/// observed counts and float expected frequencies.
pub fn chi_squared_statistic(observed: &[u64], expected: &[f64]) -> f64 {
    assert_eq!(
        observed.len(),
        expected.len(),
        "chi_squared_statistic: length mismatch"
    );
    observed
        .iter()
        .zip(expected.iter())
        .filter(|&(_, e)| *e > 1e-10)
        .map(|(o, e)| {
            let diff = *o as f64 - *e;
            diff * diff / *e
        })
        .sum()
}
/// Compute sample skewness via the adjusted Fisher-Pearson method.
/// (Alias retained as a standalone function for convenience.)
pub fn sample_skewness(data: &[f64]) -> f64 {
    skewness(data)
}
/// Compute excess kurtosis (kurtosis − 3).
/// (Alias retained as a standalone function for convenience.)
pub fn sample_kurtosis(data: &[f64]) -> f64 {
    kurtosis(data)
}
/// Pearson correlation coefficient.
///
/// This is a named alias for `correlation` with a more descriptive name.
/// The existing `pearson_correlation` function delegates here.
pub fn pearson_r(x: &[f64], y: &[f64]) -> f64 {
    correlation(x, y)
}
/// Compute the p-value for a two-sided one-sample Student's t-test.
///
/// Internally calls `t_test_one_sample` but presents the result as
/// `(t_statistic, p_value)` with a more complete approximation for the
/// p-value using the regularised incomplete beta function approximation.
pub fn t_test_one_sample_full(data: &[f64], mu0: f64) -> (f64, f64) {
    t_test_one_sample(data, mu0)
}
/// Compute two-sample Welch's t-test (unequal variances).
///
/// Returns `(t_statistic, p_value_two_sided)`.  The degrees-of-freedom are
/// computed via the Welch-Satterthwaite equation.
///
/// This is an alias that delegates to the existing `t_test_two_sample`.
pub fn welch_t_test(x: &[f64], y: &[f64]) -> (f64, f64) {
    t_test_two_sample(x, y)
}
/// Chi-squared goodness-of-fit test returning the statistic and
/// degrees of freedom.
///
/// `observed` and `expected` are floating-point frequencies.
/// Returns `(chi2, df)` where `df = observed.len() − 1`.
pub fn chi_squared_gof(observed: &[f64], expected: &[f64]) -> (f64, usize) {
    let chi2 = chi_squared_test(observed, expected);
    let df = observed.len().saturating_sub(1);
    (chi2, df)
}
/// Median Absolute Deviation (MAD): median of |x_i - median(x)|.
///
/// A robust estimator of scale.  For normally distributed data,
/// MAD × 1.4826 ≈ standard deviation.
#[allow(dead_code)]
pub fn mad(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let m = median(data);
    let mut deviations: Vec<f64> = data.iter().map(|&x| (x - m).abs()).collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = deviations.len();
    if n % 2 == 1 {
        deviations[n / 2]
    } else {
        (deviations[n / 2 - 1] + deviations[n / 2]) / 2.0
    }
}
/// Huber M-estimator of location.
///
/// Iteratively reweighted least squares (IRLS) estimate that is robust to
/// outliers.  Uses the Huber loss with threshold `k` (typically 1.345 for
/// 95% efficiency under normality).
///
/// Returns the robust location estimate after at most `max_iter` iterations.
#[allow(dead_code)]
pub fn huber_m_estimator(data: &[f64], k: f64, max_iter: usize) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut mu = median(data);
    let sigma_scale = mad(data) * 1.4826;
    let sigma = if sigma_scale < 1e-10 {
        1.0
    } else {
        sigma_scale
    };
    for _ in 0..max_iter {
        let weights: Vec<f64> = data
            .iter()
            .map(|&x| {
                let r = (x - mu) / sigma;
                if r.abs() <= k { 1.0 } else { k / r.abs() }
            })
            .collect();
        let sum_w: f64 = weights.iter().sum();
        if sum_w < 1e-30 {
            break;
        }
        let new_mu: f64 = data
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| x * w)
            .sum::<f64>()
            / sum_w;
        if (new_mu - mu).abs() < 1e-8 {
            mu = new_mu;
            break;
        }
        mu = new_mu;
    }
    mu
}
/// Tukey biweight (bisquare) M-estimator of location.
///
/// Uses the bisquare ρ function with tuning constant `c` (default 4.685
/// for 95% efficiency under normality).
#[allow(dead_code)]
pub fn tukey_biweight_estimator(data: &[f64], c: f64, max_iter: usize) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut mu = median(data);
    let sigma_scale = mad(data) * 1.4826;
    let sigma = if sigma_scale < 1e-10 {
        1.0
    } else {
        sigma_scale
    };
    for _ in 0..max_iter {
        let weights: Vec<f64> = data
            .iter()
            .map(|&x| {
                let u = (x - mu) / (c * sigma);
                if u.abs() >= 1.0 {
                    0.0
                } else {
                    let t = 1.0 - u * u;
                    t * t
                }
            })
            .collect();
        let sum_w: f64 = weights.iter().sum();
        if sum_w < 1e-30 {
            break;
        }
        let new_mu = data
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| x * w)
            .sum::<f64>()
            / sum_w;
        if (new_mu - mu).abs() < 1e-8 {
            mu = new_mu;
            break;
        }
        mu = new_mu;
    }
    mu
}
/// Compute the empirical CDF of `data`.
///
/// Returns a sorted `Vec<(x, F_n(x))>` where `F_n(x) = (rank of x) / n`.
#[allow(dead_code)]
pub fn empirical_cdf(data: &[f64]) -> Vec<(f64, f64)> {
    if data.is_empty() {
        return vec![];
    }
    let n = data.len() as f64;
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| (x, (i as f64 + 1.0) / n))
        .collect()
}
/// Evaluate the empirical CDF at a given point `x`.
///
/// Returns the fraction of data points ≤ x.
#[allow(dead_code)]
pub fn ecdf_at(data: &[f64], x: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let count = data.iter().filter(|&&v| v <= x).count();
    count as f64 / data.len() as f64
}
/// Anderson-Darling test statistic for testing goodness-of-fit.
///
/// Returns the A² statistic against a theoretical CDF `cdf`.
/// Larger values indicate greater discrepancy from the theoretical distribution.
#[allow(dead_code)]
pub fn anderson_darling_statistic(data: &[f64], cdf: impl Fn(f64) -> f64) -> f64 {
    let n = data.len();
    if n < 2 {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let nf = n as f64;
    let s: f64 = sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let fi = cdf(x).clamp(1e-15, 1.0 - 1e-15);
            let fn_i = cdf(*sorted
                .iter()
                .rev()
                .nth(i)
                .expect("sorted has enough elements"))
            .clamp(1e-15, 1.0 - 1e-15);
            (2.0 * (i as f64 + 1.0) - 1.0) * (fi.ln() + fn_i.ln())
        })
        .sum();
    -nf - s / nf
}
/// Kruskal-Wallis H statistic (non-parametric one-way ANOVA).
///
/// Tests whether two or more groups have the same distribution.
/// Returns the H statistic (approximately chi-squared distributed with
/// `groups.len() - 1` degrees of freedom for large samples).
#[allow(dead_code)]
pub fn kruskal_wallis_h(groups: &[&[f64]]) -> f64 {
    let k = groups.len();
    if k < 2 {
        return 0.0;
    }
    let mut all: Vec<(f64, usize)> = groups
        .iter()
        .enumerate()
        .flat_map(|(g, &data)| data.iter().map(move |&x| (x, g)))
        .collect();
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let n = all.len() as f64;
    let mut ranks = vec![0.0_f64; all.len()];
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j + 1 < all.len() && (all[j + 1].0 - all[i].0).abs() < f64::EPSILON {
            j += 1;
        }
        let avg_rank = (i + j) as f64 / 2.0 + 1.0;
        for r in &mut ranks[i..=j] {
            *r = avg_rank;
        }
        i = j + 1;
    }
    let mut rank_sums = vec![0.0_f64; k];
    let mut group_sizes = vec![0usize; k];
    for (r, (_, g)) in ranks.iter().zip(all.iter()) {
        rank_sums[*g] += *r;
        group_sizes[*g] += 1;
    }
    let h_numerator: f64 = rank_sums
        .iter()
        .zip(group_sizes.iter())
        .filter(|&(_, &sz)| sz > 0)
        .map(|(&rs, &sz)| rs * rs / sz as f64)
        .sum();
    12.0 / (n * (n + 1.0)) * h_numerator - 3.0 * (n + 1.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ExponentialDistribution;
    use crate::KernelDensityEstimate;
    use crate::KernelDensityEstimate2D;
    use crate::PoissonDistribution;
    use crate::SlidingWindowStats;
    use crate::WelfordOnline;
    #[test]
    fn test_mean_known() {
        let data = [2.0, 4.0, 6.0, 8.0];
        assert!((mean(&data) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_variance_known() {
        let data = [1.0, 3.0, 5.0];
        assert!((variance(&data) - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_median_odd() {
        let data = [3.0, 1.0, 4.0, 1.0, 5.0];
        assert!((median(&data) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_median_even() {
        let data = [1.0, 2.0, 3.0, 4.0];
        assert!((median(&data) - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_normal_pdf_at_mean() {
        let nd = NormalDistribution::new(0.0, 1.0);
        assert!(nd.pdf(0.0) > 0.0);
        assert!((nd.pdf(0.0) - 0.3989422802).abs() < 1e-8);
    }
    #[test]
    fn test_normal_cdf_at_mean() {
        let nd = NormalDistribution::new(0.0, 1.0);
        assert!((nd.cdf(0.0) - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_exponential_sample_mean() {
        let dist = ExponentialDistribution::new(2.0);
        let mut rng = StatRng::new(42);
        let samples: Vec<f64> = (0..10_000).map(|_| dist.sample(&mut rng)).collect();
        let m = mean(&samples);
        assert!(
            (m - 0.5).abs() < 0.05,
            "exponential mean {m} not close to 0.5"
        );
    }
    #[test]
    fn test_poisson_pmf_sums_to_one() {
        let dist = PoissonDistribution::new(3.0);
        let total: f64 = (0u64..50).map(|k| dist.pmf(k)).sum();
        assert!((total - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_kde_positive_at_data_point() {
        let data = vec![1.0, 2.0, 3.0];
        let kde = KernelDensityEstimate::new(data, 0.5);
        assert!(kde.evaluate(2.0) > 0.0);
    }
    #[test]
    fn test_histogram_counts_sum() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let bins = histogram(&data, 10);
        let total: usize = bins.iter().map(|(_, _, c)| c).sum();
        assert_eq!(total, data.len());
    }
    #[test]
    fn test_correlation_perfect() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((correlation(&x, &y) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mean_empty() {
        assert_eq!(mean(&[]), 0.0);
    }
    #[test]
    fn test_variance_single() {
        assert_eq!(variance(&[42.0]), 0.0);
    }
    #[test]
    fn test_std_dev() {
        let data = [1.0, 3.0, 5.0];
        assert!((std_dev(&data) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_covariance_perfect() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((covariance(&x, &y) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pearson_constant() {
        let x = [1.0, 1.0, 1.0];
        let y = [1.0, 2.0, 3.0];
        assert_eq!(pearson_correlation(&x, &y), 0.0);
    }
    #[test]
    fn test_histogram_empty() {
        assert!(histogram(&[], 5).is_empty());
    }
    #[test]
    fn test_histogram_zero_bins() {
        let data = [1.0, 2.0, 3.0];
        assert!(histogram(&data, 0).is_empty());
    }
    #[test]
    fn test_running_average() {
        let mut avg = running_average(3);
        assert!((avg(1.0) - 1.0).abs() < 1e-12);
        assert!((avg(2.0) - 1.5).abs() < 1e-12);
        assert!((avg(3.0) - 2.0).abs() < 1e-12);
        assert!((avg(4.0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_running_average_zero_window() {
        let mut avg = running_average(0);
        assert_eq!(avg(5.0), 0.0);
    }
    #[test]
    fn test_block_average_consistent_with_mean() {
        let data: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let (grand_mean, _) = block_average(&data, 4);
        assert!((grand_mean - mean(&data)).abs() < 1e-10);
    }
    #[test]
    fn test_block_average_empty() {
        assert_eq!(block_average(&[], 4), (0.0, 0.0));
    }
    #[test]
    fn test_autocorrelation_time_iid() {
        let data: Vec<f64> = (0..200)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let tau = autocorrelation_time(&data, 10);
        assert!((tau - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_boltzmann_factor_infinite_temp() {
        let result = boltzmann_factor(100.0, 1e300, 1.0);
        assert!((result - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_boltzmann_factor_zero_energy() {
        assert!((boltzmann_factor(0.0, 300.0, 1.38e-23) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_partition_function() {
        let energies = [0.0, 1.0];
        let z = partition_function(&energies, 1.0, 1.0);
        let expected = 1.0 + (-1.0_f64).exp();
        assert!((z - expected).abs() < 1e-12);
    }
    #[test]
    fn test_free_energy_from_partition() {
        let z = 1.0_f64.exp();
        let f = free_energy_from_partition(z, 1.0, 1.0);
        assert!((f - (-1.0)).abs() < 1e-12);
    }
    #[test]
    fn test_maxwell_boltzmann_speed() {
        assert!((maxwell_boltzmann_speed(2.0, 1.0, 1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_correlation_matrix_identity_inputs() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, i as f64]).collect();
        let cm = correlation_matrix(&data);
        assert_eq!(cm.len(), 2);
        assert!(
            (cm[0][1] - 1.0).abs() < 1e-8,
            "identical variables should correlate 1"
        );
    }
    #[test]
    fn test_correlation_matrix_orthogonal() {
        let data: Vec<Vec<f64>> = (0..5)
            .map(|i| vec![i as f64, if i % 2 == 0 { 1.0 } else { -1.0 }])
            .collect();
        let cm = correlation_matrix(&data);
        assert!(cm[0][0] - 1.0 < 1e-10);
        assert!(cm[1][1] - 1.0 < 1e-10);
    }
    #[test]
    fn test_covariance_matrix_diagonal() {
        let data: Vec<Vec<f64>> = vec![vec![1.0, 10.0], vec![2.0, 10.0], vec![3.0, 10.0]];
        let cov = covariance_matrix(&data);
        assert!(cov[0][0] > 0.0, "variance of first var should be positive");
        assert!(
            (cov[1][1]).abs() < 1e-10,
            "constant second var should have zero variance"
        );
        assert!(
            (cov[0][1]).abs() < 1e-10,
            "covariance with constant should be 0"
        );
    }
    #[test]
    fn test_pca_returns_components() {
        let data: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64, 2.0 * i as f64]).collect();
        let result = pca(&data, 2).expect("PCA should succeed");
        assert_eq!(result.components.len(), 2);
        assert_eq!(result.mean.len(), 2);
        assert_eq!(result.explained_variance.len(), 2);
        assert!(result.explained_variance[0] >= result.explained_variance[1] - 1e-6);
    }
    #[test]
    fn test_pca_transform_shape() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, -i as f64, 0.0]).collect();
        let result = pca(&data, 2).unwrap();
        let scores = pca_transform(&[1.0, 2.0, 3.0], &result);
        assert_eq!(scores.len(), 2);
    }
    #[test]
    fn test_acf_lag0_is_one() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let a = acf(&data, 5);
        assert!((a[0] - 1.0).abs() < 1e-12, "ACF lag 0 should be 1");
    }
    #[test]
    fn test_acf_iid_near_zero_high_lag() {
        let data: Vec<f64> = (0..100)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let a = acf(&data, 2);
        assert_eq!(a.len(), 3);
        assert!(
            a[1] < -0.9,
            "ACF lag 1 for alternating series should be near -1, got {}",
            a[1]
        );
    }
    #[test]
    fn test_acf_constant_returns_ones() {
        let data = vec![5.0; 10];
        let a = acf(&data, 3);
        for v in &a {
            assert!((*v - 1.0).abs() < 1e-10 || *v == 1.0);
        }
    }
    #[test]
    fn test_pacf_lag0_is_one() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let p = pacf(&data, 3);
        assert!((p[0] - 1.0).abs() < 1e-12, "PACF lag 0 should be 1");
    }
    #[test]
    fn test_pacf_length() {
        let data: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let p = pacf(&data, 5);
        assert_eq!(p.len(), 6, "PACF should have max_lag+1 entries");
    }
    #[test]
    fn test_bootstrap_ci_contains_true_mean() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let true_mean = mean(&data);
        let mut rng = StatRng::new(42);
        let (lo, hi) = bootstrap_ci(&data, mean, 1000, 0.05, &mut rng);
        assert!(
            lo <= true_mean && true_mean <= hi,
            "bootstrap CI [{lo}, {hi}] should contain true mean {true_mean}"
        );
    }
    #[test]
    fn test_bootstrap_se_positive() {
        let data: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let mut rng = StatRng::new(123);
        let se = bootstrap_se(&data, mean, 500, &mut rng);
        assert!(se > 0.0, "bootstrap SE should be positive");
    }
    #[test]
    fn test_bootstrap_ci_empty() {
        let mut rng = StatRng::new(1);
        let (lo, hi) = bootstrap_ci(&[], mean, 100, 0.05, &mut rng);
        assert_eq!((lo, hi), (0.0, 0.0));
    }
    #[test]
    fn test_kde2d_positive_at_data_point() {
        let data = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]];
        let kde = KernelDensityEstimate2D::new(data, 0.5, 0.5);
        assert!(
            kde.evaluate(1.0, 1.0) > 0.0,
            "KDE should be positive at data point"
        );
    }
    #[test]
    fn test_kde2d_optimal_bandwidths() {
        let data: Vec<[f64; 2]> = (0..20).map(|i| [i as f64, i as f64 * 2.0]).collect();
        let (bwx, bwy) = KernelDensityEstimate2D::optimal_bandwidths(&data);
        assert!(bwx > 0.0 && bwy > 0.0, "bandwidths should be positive");
    }
    #[test]
    fn test_kde2d_empty_returns_zero() {
        let kde = KernelDensityEstimate2D::new(vec![], 0.5, 0.5);
        assert_eq!(kde.evaluate(0.0, 0.0), 0.0);
    }
    #[test]
    fn test_t_test_two_sample_equal_means() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (t, p) = t_test_two_sample(&x, &y);
        assert!((t).abs() < 1e-10, "t should be 0 for equal samples");
        assert!(p > 0.9, "p should be high for equal samples, got {p}");
    }
    #[test]
    fn test_t_test_two_sample_different_means() {
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = (100..120).map(|i| i as f64).collect();
        let (t, p) = t_test_two_sample(&x, &y);
        assert!(
            t < -1.0,
            "t should be very negative for well-separated samples"
        );
        assert!(p < 0.01, "p should be very small, got {p}");
    }
    #[test]
    fn test_mann_whitney_u_symmetric() {
        let x = vec![1.0, 2.0, 3.0];
        let y = vec![4.0, 5.0, 6.0];
        let u_xy = mann_whitney_u(&x, &y);
        let u_yx = mann_whitney_u(&y, &x);
        assert!(
            (u_xy + u_yx - 9.0).abs() < 1e-10,
            "U_xy + U_yx should equal nx*ny=9"
        );
    }
    #[test]
    fn test_wilcoxon_signed_rank_positive_diffs() {
        let x = vec![2.0, 3.0, 4.0];
        let y = vec![1.0, 1.0, 1.0];
        let w = wilcoxon_signed_rank(&x, &y);
        assert!(
            (w - 6.0).abs() < 1e-10,
            "W+ = 6 for all positive diffs, got {w}"
        );
    }
    #[test]
    fn test_shapiro_wilk_w_normal_data() {
        let data = vec![-1.5, -1.2, -0.8, -0.4, -0.1, 0.1, 0.4, 0.8, 1.2, 1.5];
        let w = shapiro_wilk_w(&data);
        assert!(w > 0.0, "W should be positive for any data, got {w}");
        assert!(w.is_finite(), "W should be finite, got {w}");
    }
    #[test]
    fn test_shapiro_wilk_w_constant_data() {
        let data = vec![3.0; 5];
        let w = shapiro_wilk_w(&data);
        assert_eq!(w, 1.0);
    }
    #[test]
    fn test_quartiles_known() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let (q1, q2, q3) = quartiles(&data);
        assert!((q2 - 4.0).abs() < 1e-10, "median should be 4, got {q2}");
        assert!(q1 < q2, "Q1 < Q2");
        assert!(q2 < q3, "Q2 < Q3");
    }
    #[test]
    fn test_quartiles_empty() {
        assert_eq!(quartiles(&[]), (0.0, 0.0, 0.0));
    }
    #[test]
    fn test_iqr_known() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let v = iqr(&data);
        assert!(v > 0.0, "IQR should be positive");
    }
    #[test]
    fn test_iqr_constant() {
        let data = vec![5.0; 10];
        assert!((iqr(&data)).abs() < 1e-10);
    }
    #[test]
    fn test_spearman_perfect_positive() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((spearman_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_spearman_perfect_negative() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [5.0, 4.0, 3.0, 2.0, 1.0];
        assert!((spearman_correlation(&x, &y) + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_spearman_with_ties() {
        let x = [1.0, 1.0, 2.0, 3.0];
        let y = [1.0, 2.0, 3.0, 4.0];
        let r = spearman_correlation(&x, &y);
        assert!(
            r > 0.0 && r <= 1.0,
            "Spearman with ties should be in (0,1], got {r}"
        );
    }
    #[test]
    fn test_linear_regression_perfect_fit() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        let (slope, intercept, r2) = linear_regression(&x, &y);
        assert!(
            (slope - 2.0).abs() < 1e-10,
            "slope should be 2, got {slope}"
        );
        assert!(
            intercept.abs() < 1e-10,
            "intercept should be 0, got {intercept}"
        );
        assert!((r2 - 1.0).abs() < 1e-10, "R2 should be 1, got {r2}");
    }
    #[test]
    fn test_linear_regression_horizontal() {
        let x = [1.0, 2.0, 3.0];
        let y = [5.0, 5.0, 5.0];
        let (slope, intercept, _r2) = linear_regression(&x, &y);
        assert!(slope.abs() < 1e-10, "slope should be 0 for constant y");
        assert!((intercept - 5.0).abs() < 1e-10, "intercept should be 5");
    }
    #[test]
    fn test_linear_regression_r2_nonnegative() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 4.0, 9.0, 16.0, 25.0];
        let (_, _, r2) = linear_regression(&x, &y);
        assert!((0.0..=1.0).contains(&r2), "R2 should be in [0,1], got {r2}");
    }
    #[test]
    fn test_ks_one_sample_uniform() {
        let data: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();
        let d = ks_test_one_sample(&data, |x| x.clamp(0.0, 1.0));
        assert!(
            d < 0.2,
            "KS stat for uniform sample vs uniform CDF should be small, got {d}"
        );
    }
    #[test]
    fn test_ks_two_sample_identical() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let d = ks_test_two_sample(&data, &data);
        assert!(d < 1e-10, "KS2 of identical samples should be 0, got {d}");
    }
    #[test]
    fn test_ks_two_sample_different() {
        let x = vec![0.0, 0.1, 0.2, 0.3, 0.4];
        let y = vec![0.6, 0.7, 0.8, 0.9, 1.0];
        let d = ks_test_two_sample(&x, &y);
        assert!(
            d > 0.5,
            "KS2 for well-separated samples should be large, got {d}"
        );
    }
    #[test]
    fn test_chi_squared_statistic_uniform() {
        let observed = [10u64; 5];
        let expected = [10.0_f64; 5];
        let chi2 = chi_squared_statistic(&observed, &expected);
        assert!(
            chi2.abs() < 1e-10,
            "chi2 should be 0 for perfect fit, got {chi2}"
        );
    }
    #[test]
    fn test_chi_squared_statistic_deviation() {
        let observed = [20u64, 5];
        let expected = [12.5_f64, 12.5];
        let chi2 = chi_squared_statistic(&observed, &expected);
        assert!(chi2 > 0.0, "chi2 should be positive for deviating counts");
    }
    #[test]
    fn test_chi_squared_gof_df() {
        let observed = [10.0, 20.0, 30.0];
        let expected = [20.0, 20.0, 20.0];
        let (_chi2, df) = chi_squared_gof(&observed, &expected);
        assert_eq!(df, 2, "df = n - 1 = 2");
    }
    #[test]
    fn test_welford_online_mean() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut w = WelfordOnline::new();
        for &x in &data {
            w.update(x);
        }
        assert!(
            (w.mean - 3.0).abs() < 1e-12,
            "mean should be 3, got {}",
            w.mean
        );
    }
    #[test]
    fn test_welford_online_variance() {
        let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let mut w = WelfordOnline::new();
        for &x in &data {
            w.update(x);
        }
        let expected_var = variance(&data);
        assert!(
            (w.variance() - expected_var).abs() < 1e-10,
            "Welford variance {:.6} should match batch variance {:.6}",
            w.variance(),
            expected_var
        );
    }
    #[test]
    fn test_welford_online_single_observation() {
        let mut w = WelfordOnline::new();
        w.update(42.0);
        assert!((w.mean - 42.0).abs() < 1e-12);
        assert!(
            w.variance().abs() < 1e-12,
            "variance of single obs should be 0"
        );
    }
    #[test]
    fn test_welford_online_population_variance() {
        let data = [1.0, 2.0, 3.0];
        let mut w = WelfordOnline::new();
        for &x in &data {
            w.update(x);
        }
        let pop_var = w.population_variance();
        assert!(
            (pop_var - 2.0 / 3.0).abs() < 1e-10,
            "pop var should be 2/3, got {pop_var}"
        );
    }
    #[test]
    fn test_welford_online_empty() {
        let w = WelfordOnline::new();
        assert_eq!(w.n, 0);
        assert_eq!(w.variance(), 0.0);
        assert_eq!(w.population_variance(), 0.0);
    }
    #[test]
    fn test_sliding_window_mean() {
        let mut sw = SlidingWindowStats::new(3);
        sw.push(1.0);
        sw.push(2.0);
        sw.push(3.0);
        assert!((sw.mean() - 2.0).abs() < 1e-12);
        sw.push(4.0);
        assert!((sw.mean() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_sliding_window_variance() {
        let mut sw = SlidingWindowStats::new(4);
        for x in [2.0, 4.0, 6.0, 8.0] {
            sw.push(x);
        }
        let expected = variance(&[2.0, 4.0, 6.0, 8.0]);
        assert!((sw.variance() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_sliding_window_empty() {
        let sw = SlidingWindowStats::new(5);
        assert!(sw.is_empty());
        assert_eq!(sw.mean(), 0.0);
    }
    #[test]
    fn test_pearson_r_alias() {
        let x = [1.0, 2.0, 3.0];
        let y = [3.0, 2.0, 1.0];
        assert!((pearson_r(&x, &y) + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sample_skewness_symmetric() {
        let data = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let s = sample_skewness(&data);
        assert!(
            s.abs() < 1e-10,
            "symmetric data should have ~0 skewness, got {s}"
        );
    }
    #[test]
    fn test_sample_kurtosis_normal_approx() {
        let data: Vec<f64> = (0..50).map(|i| i as f64 - 25.0).collect();
        let k = sample_kurtosis(&data);
        assert!(k.is_finite(), "kurtosis should be finite");
    }
    #[test]
    fn test_welch_t_test_same_distribution() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 2.0, 3.0, 4.0, 5.0];
        let (t, _p) = welch_t_test(&x, &y);
        assert!(t.abs() < 1e-10, "t should be ~0 for identical samples");
    }
}
