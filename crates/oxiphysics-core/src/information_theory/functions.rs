//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// x * log2(x), defined as 0 when x <= 0.
#[inline]
pub fn xlog2x(x: f64) -> f64 {
    if x <= 0.0 { 0.0 } else { x * x.log2() }
}
/// x * ln(x), defined as 0 when x <= 0.
#[inline]
pub fn xlnx(x: f64) -> f64 {
    if x <= 0.0 { 0.0 } else { x * x.ln() }
}
/// Marginal distribution over rows (sum_j p(i,j)).
pub fn marginal_x(joint: &[Vec<f64>]) -> Vec<f64> {
    joint.iter().map(|row| row.iter().sum()).collect()
}
/// Marginal distribution over columns (sum_i p(i,j)).
pub fn marginal_y(joint: &[Vec<f64>]) -> Vec<f64> {
    if joint.is_empty() {
        return vec![];
    }
    let ncols = joint[0].len();
    let mut py = vec![0.0_f64; ncols];
    for row in joint {
        for (j, &v) in row.iter().enumerate() {
            py[j] += v;
        }
    }
    py
}
/// Normalise a probability vector in-place. Returns false if sum is zero.
pub fn normalise_inplace(v: &mut [f64]) -> bool {
    let s: f64 = v.iter().sum();
    if s <= 0.0 {
        return false;
    }
    for x in v.iter_mut() {
        *x /= s;
    }
    true
}
/// Softmax of a slice (for numerical stability).
pub fn softmax(logits: &[f64]) -> Vec<f64> {
    if logits.is_empty() {
        return vec![];
    }
    let max_val = logits.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = logits.iter().map(|&x| (x - max_val).exp()).collect();
    let s: f64 = exps.iter().sum();
    exps.iter().map(|&e| e / s).collect()
}
/// Shannon entropy H(X) = -sum p * log2(p), in bits.
///
/// Returns `0.0` for an empty slice.
pub fn entropy_bits(probs: &[f64]) -> f64 {
    -probs.iter().map(|&p| xlog2x(p)).sum::<f64>()
}
/// Shannon entropy H(X) = -sum p * ln(p), in nats.
///
/// Returns `0.0` for an empty slice.
pub fn entropy_nats(probs: &[f64]) -> f64 {
    -probs.iter().map(|&p| xlnx(p)).sum::<f64>()
}
/// Joint entropy H(X,Y) = -sum_{i,j} p(i,j) * log2(p(i,j)), in bits.
///
/// `joint[i][j]` is the probability of outcome (i, j).
pub fn joint_entropy(joint: &[Vec<f64>]) -> f64 {
    -joint
        .iter()
        .flat_map(|row| row.iter())
        .map(|&p| xlog2x(p))
        .sum::<f64>()
}
/// Conditional entropy H(Y|X) = H(X,Y) - H(X), in bits.
///
/// `joint[i][j]` is the probability of outcome (i, j).
pub fn conditional_entropy(joint: &[Vec<f64>]) -> f64 {
    let px = marginal_x(joint);
    let hx = entropy_bits(&px);
    let hxy = joint_entropy(joint);
    hxy - hx
}
/// Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y), in bits.
///
/// `joint[i][j]` is the probability of outcome (i, j).
pub fn mutual_information(joint: &[Vec<f64>]) -> f64 {
    let px = marginal_x(joint);
    let py = marginal_y(joint);
    entropy_bits(&px) + entropy_bits(&py) - joint_entropy(joint)
}
/// Normalised mutual information NMI = I(X;Y) / sqrt(H(X) * H(Y)).
///
/// Returns `0.0` when either marginal entropy is zero.
pub fn normalised_mutual_information(joint: &[Vec<f64>]) -> f64 {
    let px = marginal_x(joint);
    let py = marginal_y(joint);
    let hx = entropy_bits(&px);
    let hy = entropy_bits(&py);
    if hx <= 0.0 || hy <= 0.0 {
        return 0.0;
    }
    mutual_information(joint) / (hx * hy).sqrt()
}
/// Variation of information VI(X;Y) = H(X,Y) - I(X;Y), in bits.
///
/// Measures the amount of information lost and gained between two variables.
pub fn variation_of_information(joint: &[Vec<f64>]) -> f64 {
    let hxy = joint_entropy(joint);
    let mi = mutual_information(joint);
    hxy - mi
}
/// Interaction information (three-variable) from joint marginals.
///
/// I(X;Y;Z) = I(X;Y) - I(X;Y|Z)
/// Approximated from pairwise joints and the triple joint.
/// `joint_xyz[i][j*nz + k]` is P(x=i, y=j, z=k).
pub fn interaction_information_3(
    joint_xy: &[Vec<f64>],
    joint_xz: &[Vec<f64>],
    joint_yz: &[Vec<f64>],
    joint_xyz: &[Vec<f64>],
) -> f64 {
    let px = marginal_x(joint_xy);
    let py = marginal_y(joint_xy);
    let pz = marginal_y(joint_xz);
    -entropy_bits(&px) - entropy_bits(&py) - entropy_bits(&pz)
        + joint_entropy(joint_xy)
        + joint_entropy(joint_xz)
        + joint_entropy(joint_yz)
        - joint_entropy(joint_xyz)
}
/// KL divergence D_KL(P || Q) = sum p * ln(p/q), in nats.
///
/// Returns `f64::INFINITY` if q\[i\] = 0 and p\[i\] > 0.
pub fn kl_divergence(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "kl_divergence: length mismatch");
    p.iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            if pi <= 0.0 {
                0.0
            } else if qi <= 0.0 {
                f64::INFINITY
            } else {
                pi * (pi / qi).ln()
            }
        })
        .sum()
}
/// KL divergence in bits: D_KL(P || Q) = sum p * log2(p/q).
///
/// Returns `f64::INFINITY` if q\[i\] = 0 and p\[i\] > 0.
pub fn kl_divergence_bits(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "kl_divergence_bits: length mismatch");
    p.iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            if pi <= 0.0 {
                0.0
            } else if qi <= 0.0 {
                f64::INFINITY
            } else {
                pi * (pi / qi).log2()
            }
        })
        .sum()
}
/// Jensen-Shannon divergence JS(P || Q) = (KL(P||M) + KL(Q||M)) / 2
/// where M = (P + Q) / 2.  Always in \[0, ln 2\].
pub fn js_divergence(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "js_divergence: length mismatch");
    let m: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (pi + qi) / 2.0)
        .collect();
    (kl_divergence(p, &m) + kl_divergence(q, &m)) / 2.0
}
/// Cross-entropy H(P, Q) = -sum p * log2(q), in bits.
///
/// Returns `f64::INFINITY` if q\[i\] = 0 and p\[i\] > 0.
pub fn cross_entropy(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "cross_entropy: length mismatch");
    -p.iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            if pi <= 0.0 {
                0.0
            } else if qi <= 0.0 {
                f64::INFINITY
            } else {
                pi * qi.log2()
            }
        })
        .sum::<f64>()
}
/// Total variation distance TV(P, Q) = 1/2 sum |p - q|.  Always in \[0, 1\].
pub fn total_variation_distance(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "tv_distance: length mismatch");
    0.5 * p
        .iter()
        .zip(q.iter())
        .map(|(&a, &b)| (a - b).abs())
        .sum::<f64>()
}
/// Hellinger distance H^2(P, Q) = 1 - sum sqrt(p*q).
///
/// The squared Hellinger distance, always in \[0, 1\].
pub fn hellinger_distance_sq(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "hellinger: length mismatch");
    let bc: f64 = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (pi * qi).sqrt())
        .sum();
    (1.0 - bc).max(0.0)
}
/// Bhattacharyya distance D_B(P, Q) = -ln(sum sqrt(p*q)).
///
/// Returns `f64::INFINITY` if Bhattacharyya coefficient is zero.
pub fn bhattacharyya_distance(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "bhattacharyya: length mismatch");
    let bc: f64 = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (pi * qi).sqrt())
        .sum();
    if bc <= 0.0 { f64::INFINITY } else { -bc.ln() }
}
/// Chi-squared divergence chi2(P || Q) = sum (p-q)^2 / q.
///
/// Returns `f64::INFINITY` if q\[i\] = 0 and p\[i\] > 0.
pub fn chi_squared_divergence(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "chi_squared: length mismatch");
    p.iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| {
            if qi <= 0.0 {
                if pi > 0.0 { f64::INFINITY } else { 0.0 }
            } else {
                (pi - qi).powi(2) / qi
            }
        })
        .sum()
}
/// General f-divergence D_f(P || Q) = sum q * f(p/q).
///
/// `f_fn` is a convex function with f(1) = 0.
/// Skips terms where q = 0.
pub fn f_divergence(p: &[f64], q: &[f64], f_fn: impl Fn(f64) -> f64) -> f64 {
    assert_eq!(p.len(), q.len(), "f_divergence: length mismatch");
    p.iter()
        .zip(q.iter())
        .map(
            |(&pi, &qi)| {
                if qi <= 0.0 { 0.0 } else { qi * f_fn(pi / qi) }
            },
        )
        .sum()
}
/// Renyi entropy H_alpha(X) = log2(sum p^alpha) / (1 - alpha), in bits.
///
/// Degenerates to Shannon entropy as alpha -> 1.
pub fn renyi_entropy(probs: &[f64], alpha: f64) -> f64 {
    if (alpha - 1.0).abs() < 1e-12 {
        return entropy_bits(probs);
    }
    let sum: f64 = probs.iter().map(|&p| p.powf(alpha)).sum();
    sum.log2() / (1.0 - alpha)
}
/// Min-entropy H_inf(X) = -log2(max p).
///
/// Returns `0.0` for an empty distribution.
pub fn min_entropy(probs: &[f64]) -> f64 {
    let p_max = probs.iter().cloned().fold(0.0_f64, f64::max);
    if p_max <= 0.0 { 0.0 } else { -p_max.log2() }
}
/// Tsallis entropy S_q(X) = (1 - sum p^q) / (q - 1).
///
/// Degenerates to Shannon entropy in the limit q -> 1.
pub fn tsallis_entropy(probs: &[f64], q: f64) -> f64 {
    if (q - 1.0).abs() < 1e-12 {
        return entropy_nats(probs);
    }
    let sum: f64 = probs.iter().map(|&p| p.powf(q)).sum();
    (1.0 - sum) / (q - 1.0)
}
/// Collision entropy H_2(X) = -log2(sum p^2), in bits.
///
/// Special case of Renyi entropy with alpha = 2.
pub fn collision_entropy(probs: &[f64]) -> f64 {
    let sum_sq: f64 = probs.iter().map(|&p| p * p).sum();
    if sum_sq <= 0.0 { 0.0 } else { -sum_sq.log2() }
}
/// Hartley entropy H_0(X) = log2(|support|), in bits.
///
/// Counts the number of symbols with positive probability.
pub fn hartley_entropy(probs: &[f64]) -> f64 {
    let count = probs.iter().filter(|&&p| p > 0.0).count() as f64;
    if count <= 0.0 { 0.0 } else { count.log2() }
}
/// Differential entropy of a Gaussian N(mu, sigma^2) = 1/2 ln(2*pi*e * sigma^2), in nats.
///
/// `sigma` is the standard deviation.
pub fn gaussian_differential_entropy(sigma: f64) -> f64 {
    0.5 * (2.0 * PI * std::f64::consts::E * sigma * sigma).ln()
}
/// Differential entropy of a uniform distribution on \[a, b\] = ln(b - a).
///
/// Returns `f64::NEG_INFINITY` if a >= b.
pub fn uniform_differential_entropy(a: f64, b: f64) -> f64 {
    let width = b - a;
    if width <= 0.0 {
        f64::NEG_INFINITY
    } else {
        width.ln()
    }
}
/// Differential entropy of an exponential distribution with rate lambda = 1 - ln(lambda).
///
/// `rate` must be positive.
pub fn exponential_differential_entropy(rate: f64) -> f64 {
    1.0 - rate.ln()
}
/// Differential entropy of a Laplace distribution with scale b: h = 1 + ln(2b).
///
/// `scale` must be positive.
pub fn laplace_differential_entropy(scale: f64) -> f64 {
    1.0 + (2.0 * scale).ln()
}
/// Differential entropy of a multivariate Gaussian with covariance determinant.
///
/// h = d/2 * ln(2*pi*e) + 1/2 * ln(det_cov), where d is the dimension.
pub fn multivariate_gaussian_entropy(dimension: usize, det_cov: f64) -> f64 {
    let d = dimension as f64;
    0.5 * d * (2.0 * PI * std::f64::consts::E).ln() + 0.5 * det_cov.ln()
}
/// Differential entropy of a gamma distribution Gamma(alpha, beta).
///
/// h = alpha - ln(beta) + ln(Gamma(alpha)) + (1-alpha)*psi(alpha)
/// where psi is the digamma function. Uses Stirling approximation for large alpha.
pub fn gamma_differential_entropy(alpha: f64, beta: f64) -> f64 {
    let ln_gamma = if alpha > 10.0 {
        (alpha - 0.5) * alpha.ln() - alpha + 0.5 * (2.0 * PI).ln()
    } else {
        let mut val = 1.0_f64;
        let mut a = alpha;
        while a < 10.0 {
            val *= a;
            a += 1.0;
        }
        let lg_a = (a - 0.5) * a.ln() - a + 0.5 * (2.0 * PI).ln();
        lg_a - val.ln()
    };
    let psi_alpha = if alpha > 5.0 {
        alpha.ln() - 0.5 / alpha
    } else {
        let mut psi = 0.0_f64;
        let mut a = alpha;
        while a < 5.0 {
            psi -= 1.0 / a;
            a += 1.0;
        }
        psi + a.ln() - 0.5 / a
    };
    alpha - beta.ln() + ln_gamma + (1.0 - alpha) * psi_alpha
}
/// Maximum entropy distribution on {0, ..., n-1} subject to a mean constraint.
///
/// Solves for the distribution p_i proportional to exp(-lambda * i) that achieves
/// the specified mean. Returns the probability vector.
///
/// `n` is the number of symbols, `target_mean` is the desired E\[X\].
pub fn max_entropy_mean_constraint(n: usize, target_mean: f64) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![1.0];
    }
    let mut lo = -20.0_f64;
    let mut hi = 20.0_f64;
    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        let probs = softmax(&(0..n).map(|i| -mid * i as f64).collect::<Vec<_>>());
        let mean: f64 = probs.iter().enumerate().map(|(i, &p)| i as f64 * p).sum();
        if mean < target_mean {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let lambda = (lo + hi) / 2.0;
    softmax(&(0..n).map(|i| -lambda * i as f64).collect::<Vec<_>>())
}
/// Maximum entropy distribution subject to mean and variance constraints.
///
/// Finds p_i proportional to exp(-lambda1 * i - lambda2 * i^2) matching
/// the desired mean and variance. Uses gradient descent on the Lagrange multipliers.
pub fn max_entropy_mean_variance_constraint(
    n: usize,
    target_mean: f64,
    target_var: f64,
) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![1.0];
    }
    let mut lambda1 = 0.0_f64;
    let mut lambda2 = 0.0_f64;
    let lr = 0.01_f64;
    for _ in 0..500 {
        let logits: Vec<f64> = (0..n)
            .map(|i| -lambda1 * i as f64 - lambda2 * (i as f64).powi(2))
            .collect();
        let probs = softmax(&logits);
        let mean: f64 = probs.iter().enumerate().map(|(i, &p)| i as f64 * p).sum();
        let var: f64 = probs
            .iter()
            .enumerate()
            .map(|(i, &p)| (i as f64 - mean).powi(2) * p)
            .sum();
        lambda1 += lr * (mean - target_mean);
        lambda2 += lr * (var - target_var);
    }
    let logits: Vec<f64> = (0..n)
        .map(|i| -lambda1 * i as f64 - lambda2 * (i as f64).powi(2))
        .collect();
    softmax(&logits)
}
/// Maximum entropy distribution on {0, ..., n-1} (uniform).
///
/// The uniform distribution maximises entropy with no constraints.
pub fn max_entropy_uniform(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    vec![1.0 / n as f64; n]
}
/// Maximum entropy rate for a Gaussian channel with power constraint P and noise N.
///
/// The capacity-achieving distribution is Gaussian, giving C = 1/2 log2(1 + P/N).
pub fn max_entropy_gaussian_channel(power: f64, noise: f64) -> f64 {
    0.5 * (1.0 + power / noise).log2()
}
/// Transfer entropy T(X->Y) with lag `lag` estimated from time-series data.
///
/// Uses a simple histogram with 8 equally-spaced bins over \[min, max\].
/// Returns 0 when any marginal is too small to be meaningful.
pub fn transfer_entropy(source: &[f64], target: &[f64], lag: usize) -> f64 {
    let n = source.len().min(target.len());
    if n <= lag + 1 {
        return 0.0;
    }
    let bins = 8_usize;
    let all_vals: Vec<f64> = source[..n]
        .iter()
        .chain(target[..n].iter())
        .copied()
        .collect();
    let lo = all_vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = all_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if (hi - lo).abs() < 1e-12 {
        return 0.0;
    }
    let bin_of = |v: f64| -> usize {
        let idx = ((v - lo) / (hi - lo) * bins as f64) as usize;
        idx.min(bins - 1)
    };
    let mut p3 = vec![0.0_f64; bins * bins * bins];
    let mut p2_ty = vec![0.0_f64; bins * bins];
    let mut p2_src_ty = vec![0.0_f64; bins * bins];
    let mut p1_ty = vec![0.0_f64; bins];
    let count = (n - lag) as f64;
    for t in lag..n {
        let yt = bin_of(target[t]);
        let yt_lag = bin_of(target[t - lag]);
        let xt_lag = bin_of(source[t - lag]);
        p3[yt * bins * bins + yt_lag * bins + xt_lag] += 1.0 / count;
        p2_ty[yt * bins + yt_lag] += 1.0 / count;
        p2_src_ty[yt_lag * bins + xt_lag] += 1.0 / count;
        p1_ty[yt_lag] += 1.0 / count;
    }
    let mut te = 0.0_f64;
    for yt in 0..bins {
        for yt_l in 0..bins {
            for xt_l in 0..bins {
                let joint = p3[yt * bins * bins + yt_l * bins + xt_l];
                if joint <= 0.0 {
                    continue;
                }
                let a = p2_ty[yt * bins + yt_l];
                let b = p2_src_ty[yt_l * bins + xt_l];
                let c = p1_ty[yt_l];
                if a > 0.0 && b > 0.0 && c > 0.0 {
                    te += joint * (joint * c / (a * b)).ln();
                }
            }
        }
    }
    te.max(0.0)
}
/// Conditional transfer entropy with a conditioning variable.
///
/// T(X->Y|Z) estimated from time-series with lag and 8-bin histogram.
pub fn conditional_transfer_entropy(
    source: &[f64],
    target: &[f64],
    cond: &[f64],
    lag: usize,
) -> f64 {
    let n = source.len().min(target.len()).min(cond.len());
    if n <= lag + 1 {
        return 0.0;
    }
    let te_xy = transfer_entropy(source, target, lag);
    let te_zy = transfer_entropy(cond, target, lag);
    (te_xy - te_zy).max(0.0)
}
/// Shannon-Hartley channel capacity C = B * log2(1 + S/N), in bits/s.
///
/// `bandwidth_hz` is the channel bandwidth; `snr` is the signal-to-noise ratio
/// (linear, not dB).
pub fn shannon_hartley_capacity(bandwidth_hz: f64, snr: f64) -> f64 {
    bandwidth_hz * (1.0 + snr).log2()
}
/// Channel capacity of a discrete memoryless channel via Blahut-Arimoto
/// (fixed 200 iterations).
///
/// `transition[i][j]` = P(output j | input i).  Returns capacity in bits.
pub fn channel_capacity_blahut(transition: &[Vec<f64>]) -> f64 {
    let n_in = transition.len();
    if n_in == 0 {
        return 0.0;
    }
    let n_out = transition[0].len();
    if n_out == 0 {
        return 0.0;
    }
    let mut q = vec![1.0 / n_in as f64; n_in];
    for _ in 0..200 {
        let mut py = vec![0.0_f64; n_out];
        for (i, qi) in q.iter().enumerate() {
            for (j, &tij) in transition[i].iter().enumerate() {
                py[j] += qi * tij;
            }
        }
        let mut c = vec![0.0_f64; n_in];
        for (ci, ti) in c.iter_mut().zip(transition.iter()) {
            let mut s = 0.0_f64;
            for (&pij, &py_j) in ti.iter().zip(py.iter()) {
                if pij > 0.0 && py_j > 0.0 {
                    s += pij * (pij / py_j).ln();
                }
            }
            *ci = s;
        }
        let c_max = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut new_q = vec![0.0_f64; n_in];
        let mut sum = 0.0_f64;
        for ((nq, &qi), &ci) in new_q.iter_mut().zip(q.iter()).zip(c.iter()) {
            *nq = qi * (ci - c_max).exp();
            sum += *nq;
        }
        for qi in &mut new_q {
            *qi /= sum;
        }
        q = new_q;
    }
    let mut py = vec![0.0_f64; n_out];
    for (i, qi) in q.iter().enumerate() {
        for (j, &tij) in transition[i].iter().enumerate() {
            py[j] += qi * tij;
        }
    }
    let mut cap = 0.0_f64;
    for (&qi, ti) in q.iter().zip(transition.iter()) {
        if qi <= 0.0 {
            continue;
        }
        for (&pij, &py_j) in ti.iter().zip(py.iter()) {
            if pij > 0.0 && py_j > 0.0 {
                cap += qi * pij * (pij / py_j).log2();
            }
        }
    }
    cap.max(0.0)
}
/// Capacity of a binary erasure channel BEC(epsilon).
///
/// C = 1 - epsilon.
pub fn bec_capacity(epsilon: f64) -> f64 {
    (1.0 - epsilon).max(0.0)
}
/// Capacity of a binary symmetric channel BSC(p).
///
/// C = 1 - H_b(p), where H_b is the binary entropy function.
pub fn bsc_capacity(p: f64) -> f64 {
    let hb = -xlog2x(p) - xlog2x(1.0 - p);
    (1.0 - hb).max(0.0)
}
/// Capacity of an AWGN channel with power constraint P and noise variance N.
///
/// C = 1/2 log2(1 + P/N).
pub fn awgn_capacity(power: f64, noise: f64) -> f64 {
    0.5 * (1.0 + power / noise).log2()
}
/// Sphere-packing (Hamming) upper bound on channel capacity.
///
/// For a binary code of length n and minimum distance d:
/// Rate <= 1 - H_b(d/(2n)) for large n.
pub fn sphere_packing_bound(n: usize, d: usize) -> f64 {
    let ratio = d as f64 / (2.0 * n as f64);
    let hb = -xlog2x(ratio) - xlog2x(1.0 - ratio);
    (1.0 - hb).max(0.0)
}
/// Simulate Huffman code lengths for a probability distribution.
///
/// Returns a vector of code lengths (in bits) for each symbol, in the same
/// order as `probs`.  Uses a min-heap priority queue simulation.
pub fn huffman_lengths(probs: &[f64]) -> Vec<usize> {
    let n = probs.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![1];
    }
    #[derive(Clone)]
    pub(super) struct Node {
        prob: f64,
        leaves: Vec<usize>,
    }
    let mut heap: Vec<Node> = probs
        .iter()
        .enumerate()
        .map(|(i, &p)| Node {
            prob: p,
            leaves: vec![i],
        })
        .collect();
    let mut lengths = vec![0_usize; n];
    while heap.len() > 1 {
        let (i1, _) = heap
            .iter()
            .enumerate()
            .min_by(|a, b| {
                a.1.prob
                    .partial_cmp(&b.1.prob)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("heap has at least 2 elements");
        let n1 = heap.remove(i1);
        let (i2, _) = heap
            .iter()
            .enumerate()
            .min_by(|a, b| {
                a.1.prob
                    .partial_cmp(&b.1.prob)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("heap has at least 2 elements");
        let n2 = heap.remove(i2);
        for &l in &n1.leaves {
            lengths[l] += 1;
        }
        for &l in &n2.leaves {
            lengths[l] += 1;
        }
        let combined_leaves: Vec<usize> = n1.leaves.into_iter().chain(n2.leaves).collect();
        heap.push(Node {
            prob: n1.prob + n2.prob,
            leaves: combined_leaves,
        });
    }
    lengths
}
/// Expected code length E\[L\] = sum p_i * l_i for a Huffman code.
pub fn huffman_expected_length(probs: &[f64]) -> f64 {
    let lengths = huffman_lengths(probs);
    probs
        .iter()
        .zip(lengths.iter())
        .map(|(&p, &l)| p * l as f64)
        .sum()
}
/// Rate-distortion for binary source with Hamming distortion.
///
/// R(D) = 1 - H_b(D) for D in \[0, 0.5\].
pub fn binary_rate_distortion(d: f64) -> f64 {
    let d = d.clamp(0.0, 0.5);
    let hd = -xlog2x(d) - xlog2x(1.0 - d);
    (1.0 - hd).max(0.0)
}
/// Gaussian rate-distortion R(D) = 1/2 * log2(sigma^2 / D), for D <= sigma^2.
///
/// `sigma_sq` is the source variance.  Returns 0 when D >= sigma^2.
pub fn gaussian_rate_distortion(sigma_sq: f64, d: f64) -> f64 {
    if d >= sigma_sq || d <= 0.0 {
        return 0.0;
    }
    0.5 * (sigma_sq / d).log2()
}
/// Rate-distortion for a discrete uniform source with alphabet size m.
///
/// Under Hamming distortion: R(D) = log2(m) - H_b(D) - D*log2(m-1) for small D.
/// Simplified to the parametric form.
pub fn uniform_rate_distortion(m: usize, d: f64) -> f64 {
    if m <= 1 {
        return 0.0;
    }
    let d = d.clamp(0.0, 1.0 - 1.0 / m as f64);
    let hd = -xlog2x(d) - xlog2x(1.0 - d);
    let r = (m as f64).log2() - hd - d * ((m - 1) as f64).log2();
    r.max(0.0)
}
/// Distortion-rate function for Gaussian source: D(R) = sigma^2 * 2^(-2R).
pub fn gaussian_distortion_at_rate(sigma_sq: f64, rate: f64) -> f64 {
    sigma_sq * 2.0_f64.powf(-2.0 * rate)
}
/// Kraft inequality: checks that lengths satisfy sum 2^(-l_i) <= 1.
pub fn kraft_inequality(lengths: &[usize]) -> bool {
    let sum: f64 = lengths.iter().map(|&l| 2.0_f64.powi(-(l as i32))).sum();
    sum <= 1.0 + 1e-9
}
/// Redundancy of a code = E\[L\] - H(X), in bits.
///
/// Positive redundancy means the code is longer than the entropy lower bound.
pub fn code_redundancy(probs: &[f64], lengths: &[usize]) -> f64 {
    assert_eq!(
        probs.len(),
        lengths.len(),
        "code_redundancy: length mismatch"
    );
    let expected_len: f64 = probs
        .iter()
        .zip(lengths.iter())
        .map(|(&p, &l)| p * l as f64)
        .sum();
    expected_len - entropy_bits(probs)
}
/// Shannon-Fano-Elias coding: theoretical code length for symbol i is
/// ceil(-log2(p_i)) + 1 bits.
pub fn shannon_fano_lengths(probs: &[f64]) -> Vec<usize> {
    probs
        .iter()
        .map(|&p| {
            if p <= 0.0 {
                0
            } else {
                (-p.log2()).ceil() as usize + 1
            }
        })
        .collect()
}
/// Akaike Information Criterion AIC = 2k - 2 ln L.
///
/// `k` = number of parameters, `log_likelihood` = ln L.
pub fn aic(k: usize, log_likelihood: f64) -> f64 {
    2.0 * k as f64 - 2.0 * log_likelihood
}
/// Corrected AIC for small samples: AICc = AIC + 2k(k+1)/(n-k-1).
pub fn aicc(k: usize, n: usize, log_likelihood: f64) -> f64 {
    let base = aic(k, log_likelihood);
    let kf = k as f64;
    let nf = n as f64;
    let denom = nf - kf - 1.0;
    if denom <= 0.0 {
        base
    } else {
        base + 2.0 * kf * (kf + 1.0) / denom
    }
}
/// Bayesian Information Criterion BIC = k ln(n) - 2 ln L.
pub fn bic(k: usize, n: usize, log_likelihood: f64) -> f64 {
    k as f64 * (n as f64).ln() - 2.0 * log_likelihood
}
/// Select the best model index by lowest AIC.
///
/// Returns `None` if `models` is empty.
pub fn select_by_aic(models: &[(usize, f64)]) -> Option<usize> {
    models
        .iter()
        .enumerate()
        .min_by(|a, b| {
            let aic_a = aic(a.1.0, a.1.1);
            let aic_b = aic(b.1.0, b.1.1);
            aic_a
                .partial_cmp(&aic_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
}
/// Minimum description length (MDL) criterion: MDL = -ln L + k/2 * ln(n).
///
/// Equivalent to BIC/2 up to a constant.
pub fn mdl(k: usize, n: usize, log_likelihood: f64) -> f64 {
    -log_likelihood + 0.5 * k as f64 * (n as f64).ln()
}
/// Fisher information matrix for a categorical distribution with `k` categories.
///
/// The Fisher metric on the probability simplex is diagonal with entries
/// g_{ii}(p) = 1/p_i.  Returns a `k x k` matrix (row-major) with 1/p_i on the
/// diagonal and 0 elsewhere.
pub fn fisher_metric_categorical(probs: &[f64]) -> Vec<f64> {
    let k = probs.len();
    let mut mat = vec![0.0_f64; k * k];
    for (i, &p) in probs.iter().enumerate() {
        if p > 0.0 {
            mat[i * k + i] = 1.0 / p;
        }
    }
    mat
}
/// Fisher information for a Gaussian model N(mu, sigma^2) w.r.t. mu.
///
/// I(mu) = 1 / sigma^2.
pub fn fisher_information_gaussian_mean(sigma: f64) -> f64 {
    if sigma <= 0.0 {
        f64::INFINITY
    } else {
        1.0 / (sigma * sigma)
    }
}
/// Fisher information for a Bernoulli(p) model.
///
/// I(p) = 1 / (p(1-p)).
pub fn fisher_information_bernoulli(p: f64) -> f64 {
    let denom = p * (1.0 - p);
    if denom <= 0.0 {
        f64::INFINITY
    } else {
        1.0 / denom
    }
}
/// Fisher information for a Poisson(lambda) model.
///
/// I(lambda) = 1 / lambda.
pub fn fisher_information_poisson(lambda: f64) -> f64 {
    if lambda <= 0.0 {
        f64::INFINITY
    } else {
        1.0 / lambda
    }
}
/// Fisher information for an exponential(lambda) model.
///
/// I(lambda) = 1 / lambda^2.
pub fn fisher_information_exponential(lambda: f64) -> f64 {
    if lambda <= 0.0 {
        f64::INFINITY
    } else {
        1.0 / (lambda * lambda)
    }
}
/// Fisher information matrix for a Gaussian N(mu, sigma^2) w.r.t. (mu, sigma).
///
/// Returns a 2x2 matrix (row-major):
/// \[\[1/sigma^2, 0\], \[0, 2/sigma^2\]\]
pub fn fisher_matrix_gaussian(sigma: f64) -> [f64; 4] {
    let s2 = sigma * sigma;
    [1.0 / s2, 0.0, 0.0, 2.0 / s2]
}
/// Numerical Fisher information from a parametric log-likelihood.
///
/// Given a function `log_p(theta, x)` and sample data, computes
/// I(theta) = E\[(d/dtheta log p)^2\] via numerical differentiation.
pub fn fisher_information_numerical(
    log_p: impl Fn(f64, f64) -> f64,
    theta: f64,
    data: &[f64],
    _h: f64,
) -> f64 {
    let h = if _h > 0.0 { _h } else { 1e-6 };
    let n = data.len() as f64;
    if n <= 0.0 {
        return 0.0;
    }
    let mut sum_sq = 0.0_f64;
    for &x in data {
        let dlog = (log_p(theta + h, x) - log_p(theta - h, x)) / (2.0 * h);
        sum_sq += dlog * dlog;
    }
    sum_sq / n
}
/// Natural gradient step: theta_new = theta - eta * G(theta)^{-1} * grad_L.
///
/// For a diagonal Fisher matrix (categorical), G^{-1}_{ii} = p_i.
pub fn natural_gradient_step_categorical(
    probs: &[f64],
    grad: &[f64],
    learning_rate: f64,
) -> Vec<f64> {
    assert_eq!(
        probs.len(),
        grad.len(),
        "natural_gradient_step: length mismatch"
    );
    probs
        .iter()
        .zip(grad.iter())
        .map(|(&p, &g)| p - learning_rate * p * g)
        .collect()
}
/// Fisher-Rao distance between two categorical distributions.
///
/// On the probability simplex, the geodesic distance under the Fisher metric is
/// d(P, Q) = 2 * arccos(sum sqrt(p_i * q_i)).
pub fn fisher_rao_distance(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "fisher_rao_distance: length mismatch");
    let bc: f64 = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| (pi * qi).sqrt())
        .sum();
    2.0 * bc.clamp(-1.0, 1.0).acos()
}
/// Fisher-Rao distance between two univariate Gaussians N(mu1, s1^2) and N(mu2, s2^2).
///
/// d = sqrt(2) * |arctan((mu2-mu1) / (sqrt(2) * s2)) - arctan((mu2-mu1) / (sqrt(2) * s1))|
/// (simplified form for the 1D case on the (mu, sigma) manifold).
///
/// Uses the exact formula: d^2 = 2*ln((s1^2+s2^2+(mu1-mu2)^2) / (2*s1*s2)).
pub fn fisher_rao_gaussian(mu1: f64, sigma1: f64, mu2: f64, sigma2: f64) -> f64 {
    let s1 = sigma1.abs();
    let s2 = sigma2.abs();
    if s1 <= 0.0 || s2 <= 0.0 {
        return f64::INFINITY;
    }
    let num = s1 * s1 + s2 * s2 + (mu1 - mu2).powi(2);
    let den = 2.0 * s1 * s2;
    (2.0 * (num / den).ln()).sqrt()
}
/// Christoffel symbols of the first kind for the Fisher-Rao metric on the simplex.
///
/// For the categorical manifold: Gamma_{ij,k} = -delta_{ijk} / (2 * p_i^2)
/// Returns the value for indices (i, j, k) given the probability vector.
pub fn christoffel_first_kind_simplex(probs: &[f64], i: usize, j: usize, k: usize) -> f64 {
    if i == j && j == k && probs[i] > 0.0 {
        -1.0 / (2.0 * probs[i] * probs[i])
    } else {
        0.0
    }
}
/// Exponential map on the probability simplex at point p in direction v.
///
/// Uses the spherical representation: p_i = xi_i^2, xi_i = sqrt(p_i).
/// The exponential map in xi-coordinates is just the geodesic on the sphere.
pub fn exponential_map_simplex(p: &[f64], v: &[f64], t: f64) -> Vec<f64> {
    assert_eq!(p.len(), v.len(), "exponential_map: length mismatch");
    let n = p.len();
    let xi: Vec<f64> = p.iter().map(|&pi| pi.sqrt()).collect();
    let mut tang: Vec<f64> = (0..n)
        .map(|i| {
            if p[i] > 0.0 {
                v[i] / (2.0 * p[i].sqrt())
            } else {
                0.0
            }
        })
        .collect();
    let dot: f64 = xi.iter().zip(tang.iter()).map(|(&x, &t_)| x * t_).sum();
    for i in 0..n {
        tang[i] -= dot * xi[i];
    }
    let tang_norm: f64 = tang.iter().map(|&t_| t_ * t_).sum::<f64>().sqrt();
    if tang_norm < 1e-15 {
        return p.to_vec();
    }
    let angle = tang_norm * t;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let mut result = vec![0.0_f64; n];
    for i in 0..n {
        let new_xi = cos_a * xi[i] + sin_a * tang[i] / tang_norm;
        result[i] = new_xi * new_xi;
    }
    let _ = normalise_inplace(&mut result);
    result
}
/// Parallel transport on the probability simplex along a geodesic.
///
/// Transports vector v from point p toward point q by parameter t.
/// Uses the spherical representation.
pub fn parallel_transport_simplex(p: &[f64], q: &[f64], v: &[f64], t: f64) -> Vec<f64> {
    let n = p.len();
    assert_eq!(n, q.len());
    assert_eq!(n, v.len());
    let d = fisher_rao_distance(p, q);
    if d < 1e-15 {
        return v.to_vec();
    }
    let mut result = vec![0.0_f64; n];
    for i in 0..n {
        let pi_t = (1.0 - t) * p[i] + t * q[i];
        if p[i] > 0.0 && pi_t > 0.0 {
            result[i] = v[i] * (pi_t / p[i]).sqrt();
        }
    }
    result
}
