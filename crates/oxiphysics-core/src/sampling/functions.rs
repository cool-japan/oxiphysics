//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{HaltonSequence, LatinHypercube, Lcg, Sobol, SobolSequence};
use rand::RngExt;

/// Sample a point uniformly on the surface of the unit sphere (Marsaglia 1972).
///
/// Returns `[x, y, z]` with `x² + y² + z² = 1`.
pub fn sample_sphere_surface(rng: &mut Lcg) -> [f64; 3] {
    loop {
        let x = rng.next_f64_range(-1.0, 1.0);
        let y = rng.next_f64_range(-1.0, 1.0);
        let s = x * x + y * y;
        if s >= 1.0 {
            continue;
        }
        let r = (1.0 - s).sqrt();
        return [2.0 * x * r, 2.0 * y * r, 1.0 - 2.0 * s];
    }
}
/// Sample a point uniformly inside the closed unit ball (rejection sampling).
///
/// Returns `[x, y, z]` with `x² + y² + z² ≤ 1`.
pub fn sample_sphere_volume(rng: &mut Lcg) -> [f64; 3] {
    loop {
        let x = rng.next_f64_range(-1.0, 1.0);
        let y = rng.next_f64_range(-1.0, 1.0);
        let z = rng.next_f64_range(-1.0, 1.0);
        if x * x + y * y + z * z <= 1.0 {
            return [x, y, z];
        }
    }
}
/// Sample a 3-D Gaussian vector with zero mean and isotropic standard deviation
/// `sigma`.
pub fn sample_gaussian_3d(rng: &mut Lcg, sigma: f64) -> [f64; 3] {
    let (z0, z1) = rng.next_normal_pair();
    let z2 = rng.next_normal();
    [sigma * z0, sigma * z1, sigma * z2]
}
/// Sample a 3-D velocity vector from the Maxwell–Boltzmann distribution.
///
/// Each velocity component is drawn from `N(0, sqrt(k_B T / m))`.
///
/// # Arguments
/// * `mass` – particle mass
/// * `temp` – temperature in Kelvin
/// * `kb`   – Boltzmann constant (use `1.380649e-23` for SI units)
pub fn sample_maxwell_boltzmann(rng: &mut Lcg, mass: f64, temp: f64, kb: f64) -> [f64; 3] {
    let sigma = (kb * temp / mass).sqrt();
    sample_gaussian_3d(rng, sigma)
}
/// Stratified 1-D sampling.
///
/// Divides `[0, 1)` into `n` equal strata of width `1/n`.  One uniform random
/// point is placed inside each stratum, so the returned vector has exactly `n`
/// elements and every stratum `[k/n, (k+1)/n)` is represented.
pub fn stratified_sample_1d(n: usize, rng: &mut Lcg) -> Vec<f64> {
    let inv_n = 1.0 / n as f64;
    (0..n)
        .map(|k| (k as f64 + rng.next_f64()) * inv_n)
        .collect()
}
/// Trait for quasi-random sequences.
pub trait QuasiRandomSequence {
    /// Return the next point in the sequence.
    fn next(&mut self) -> Vec<f64>;
}
/// Generate 2D Poisson disk (blue noise) samples using Bridson's algorithm.
///
/// Returns at most `n` sample points in the unit square \[0,1\]² with minimum
/// separation distance `r`. The actual number of generated points may be less
/// than `n` if the domain cannot accommodate more.
pub fn blue_noise_2d(n: usize, r: f64, rng: &mut Lcg) -> Vec<[f64; 2]> {
    let mut samples: Vec<[f64; 2]> = Vec::with_capacity(n);
    let mut active: Vec<usize> = Vec::new();
    let r2 = r * r;
    let k = 30;
    let first = [rng.next_f64(), rng.next_f64()];
    samples.push(first);
    active.push(0);
    while !active.is_empty() && samples.len() < n {
        let idx = (rng.next_u64() as usize) % active.len();
        let parent_idx = active[idx];
        let parent = samples[parent_idx];
        let mut found = false;
        for _ in 0..k {
            let angle = rng.next_f64_range(0.0, 2.0 * std::f64::consts::PI);
            let radius = rng.next_f64_range(r, 2.0 * r);
            let candidate = [
                parent[0] + radius * angle.cos(),
                parent[1] + radius * angle.sin(),
            ];
            if candidate[0] < 0.0
                || candidate[0] >= 1.0
                || candidate[1] < 0.0
                || candidate[1] >= 1.0
            {
                continue;
            }
            let valid = samples.iter().all(|&s| {
                let dx = s[0] - candidate[0];
                let dy = s[1] - candidate[1];
                dx * dx + dy * dy >= r2
            });
            if valid {
                let new_idx = samples.len();
                samples.push(candidate);
                active.push(new_idx);
                found = true;
                break;
            }
        }
        if !found {
            active.remove(idx);
        }
    }
    samples
}
/// Generate `n_samples` stratified samples in `n_dims` dimensions.
///
/// Each dimension is independently stratified into `n_samples` bins;
/// samples are uniformly jittered within each bin.  This is a simplified
/// (independent-per-dimension) version with no joint stratification.
pub fn stratified_sample_nd(n_samples: usize, n_dims: usize, rng: &mut Lcg) -> Vec<Vec<f64>> {
    let inv_n = 1.0 / n_samples as f64;
    let mut result: Vec<Vec<f64>> = (0..n_samples).map(|_| Vec::with_capacity(n_dims)).collect();
    for _ in 0..n_dims {
        let mut strata: Vec<usize> = (0..n_samples).collect();
        for i in (1..n_samples).rev() {
            let j = (rng.next_u64() as usize) % (i + 1);
            strata.swap(i, j);
        }
        for (res, &s) in result.iter_mut().zip(strata.iter()) {
            let val = (s as f64 + rng.next_f64()) * inv_n;
            res.push(val);
        }
    }
    result
}
/// Systematic resampling for particle filters.
///
/// Given `particles` and their (unnormalized) `weights`, draw `n_out` particles
/// with replacement proportional to weights.
pub fn systematic_resample(
    particles: &[f64],
    weights: &[f64],
    n_out: usize,
    rng: &mut Lcg,
) -> Vec<f64> {
    let n = particles.len().min(weights.len());
    if n == 0 || n_out == 0 {
        return Vec::new();
    }
    let total: f64 = weights[..n].iter().sum();
    let norm: Vec<f64> = weights[..n].iter().map(|&w| w / total).collect();
    let mut cumsum = Vec::with_capacity(n);
    let mut acc = 0.0;
    for &w in &norm {
        acc += w;
        cumsum.push(acc);
    }
    let u0 = rng.next_f64() / n_out as f64;
    let step = 1.0 / n_out as f64;
    let mut result = Vec::with_capacity(n_out);
    let mut j = 0usize;
    for i in 0..n_out {
        let u = u0 + i as f64 * step;
        while j < n - 1 && cumsum[j] < u {
            j += 1;
        }
        result.push(particles[j]);
    }
    result
}
/// Generate `n` multi-dimensional Halton points using the first `n_dims` primes as bases.
///
/// Returns `n` points, each of dimension `n_dims`, in `[0,1)^{n_dims}`.
pub fn halton_multivariate(n: usize, n_dims: usize) -> Vec<Vec<f64>> {
    pub(super) const PRIMES: [u32; 10] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    let d = n_dims.min(PRIMES.len());
    (1..=n)
        .map(|i| {
            (0..d)
                .map(|k| HaltonSequence::van_der_corput(i as u32, PRIMES[k]))
                .collect()
        })
        .collect()
}
/// Latin Hypercube Sample: `n` points in `d` dimensions, each in `[0,1)`.
///
/// Uses `rand::Rng` for shuffling and jittering.
pub fn latin_hypercube_sample(n: usize, d: usize, rng: &mut impl rand::Rng) -> Vec<Vec<f64>> {
    let inv_n = 1.0 / n as f64;
    let mut result: Vec<Vec<f64>> = (0..n).map(|_| Vec::with_capacity(d)).collect();
    for _ in 0..d {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            perm.swap(i, j);
        }
        for i in 0..n {
            let jitter: f64 = rng.random();
            let val = (perm[i] as f64 + jitter) * inv_n;
            result[i].push(val);
        }
    }
    result
}
/// Sobol' quasi-random sequence: `n` points in `d` dimensions (up to 3 dims).
///
/// Uses the same direction numbers as [`SobolSequence`].
pub fn sobol_sequence(n: usize, d: usize) -> Vec<Vec<f64>> {
    SobolSequence::new(d).sample(n)
}
/// Stratified (jittered) 1-D samples using `rand::Rng`.
///
/// Divides `[0,1)` into `n` strata; one jittered sample per stratum.
pub fn stratified_sample_1d_rng(n: usize, rng: &mut impl rand::Rng) -> Vec<f64> {
    let inv_n = 1.0 / n as f64;
    (0..n)
        .map(|k| {
            let jitter: f64 = rng.random();
            (k as f64 + jitter) * inv_n
        })
        .collect()
}
/// Discrete inverse-CDF (importance) sampling.
///
/// Draws `n` indices from `0..pdf.len()` with probability proportional to
/// `pdf[i]`.  Returns indices (with replacement).
///
/// # Panics
/// Panics if `pdf` is empty or all weights are zero.
pub fn importance_sample(pdf: &[f64], n: usize, rng: &mut impl rand::Rng) -> Vec<usize> {
    assert!(!pdf.is_empty(), "pdf must not be empty");
    let total: f64 = pdf.iter().sum();
    assert!(total > 0.0, "pdf weights must sum to a positive value");
    let mut cdf = Vec::with_capacity(pdf.len());
    let mut acc = 0.0;
    for &w in pdf {
        acc += w / total;
        cdf.push(acc);
    }
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }
    let mut samples = Vec::with_capacity(n);
    for _ in 0..n {
        let u: f64 = rng.random();
        let idx = cdf.partition_point(|&c| c < u);
        samples.push(idx.min(pdf.len() - 1));
    }
    samples
}
/// 2-D rejection sampling.
///
/// Draws `n` uniform samples from `[x_range.0, x_range.1) × [y_range.0, y_range.1)`
/// and accepts `(x, y)` when `u < f(x, y)` where `u ~ U[0, max_f]`.
///
/// `f` must be non-negative; the function evaluates an internal upper bound
/// by sampling a small grid.
pub fn rejection_sample_2d<F>(
    f: F,
    x_range: (f64, f64),
    y_range: (f64, f64),
    n: usize,
    rng: &mut impl rand::Rng,
) -> Vec<[f64; 2]>
where
    F: Fn(f64, f64) -> f64,
{
    let mut max_f = 1e-30_f64;
    for ix in 0..20 {
        for iy in 0..20 {
            let x = x_range.0 + (ix as f64 + 0.5) / 20.0 * (x_range.1 - x_range.0);
            let y = y_range.0 + (iy as f64 + 0.5) / 20.0 * (y_range.1 - y_range.0);
            max_f = max_f.max(f(x, y));
        }
    }
    if max_f <= 0.0 {
        return Vec::new();
    }
    let mut samples = Vec::with_capacity(n);
    while samples.len() < n {
        let x: f64 = x_range.0 + rng.random::<f64>() * (x_range.1 - x_range.0);
        let y: f64 = y_range.0 + rng.random::<f64>() * (y_range.1 - y_range.0);
        let u: f64 = rng.random::<f64>() * max_f;
        if u <= f(x, y) {
            samples.push([x, y]);
        }
    }
    samples
}
/// Van der Corput / Halton sequence in a given `base`.  Returns `n` values in `[0,1)`.
pub fn halton_sequence(n: usize, base: u32) -> Vec<f64> {
    HaltonSequence::sample(n, base)
}
/// Direction numbers for dimensions 1–10 of the Sobol sequence.
///
/// These are the standard Joe & Kuo direction numbers for `s ≤ 10` dimensions.
/// Dimension 0 is the trivial sequence `1<<(BITS-1-i)`.
/// Dimensions 1–9 use the published primitive polynomials and initialisation
/// values from Joe & Kuo (2003).
pub(super) fn sobol10_direction_numbers(dim: usize) -> Vec<u32> {
    pub(super) const BITS: usize = 32;
    match dim {
        0 => (0..BITS).map(|i| 1u32 << (BITS - 1 - i)).collect(),
        1 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            for i in 1..BITS {
                v[i] = v[i - 1] ^ (v[i - 1] >> 1);
            }
            v
        }
        2 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 1u32 << (BITS - 2);
            for i in 2..BITS {
                v[i] = v[i - 2] ^ (v[i - 2] >> 2) ^ v[i - 1] ^ (v[i - 1] >> 1);
            }
            v
        }
        3 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 1u32 << (BITS - 2);
            v[2] = (3u32 << (BITS - 3)).wrapping_shr(0);
            v[2] = 3u32 << (BITS - 3);
            for i in 3..BITS {
                v[i] = v[i - 3] ^ (v[i - 3] >> 3) ^ v[i - 1] ^ (v[i - 1] >> 1);
            }
            v
        }
        4 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 3u32 << (BITS - 3);
            v[2] = 7u32 << (BITS - 4);
            for i in 3..BITS {
                v[i] = v[i - 3]
                    ^ (v[i - 3] >> 3)
                    ^ (v[i - 2] << 1)
                    ^ (v[i - 2] >> 1)
                    ^ v[i - 1]
                    ^ (v[i - 1] >> 1);
            }
            v
        }
        5 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 1u32 << (BITS - 2);
            v[2] = 5u32 << (BITS - 4);
            v[3] = 13u32 << (BITS - 5);
            for i in 4..BITS {
                v[i] = v[i - 4] ^ (v[i - 4] >> 4) ^ v[i - 1] ^ (v[i - 1] >> 1);
            }
            v
        }
        6 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 3u32 << (BITS - 3);
            v[2] = 13u32 << (BITS - 5);
            v[3] = 5u32 << (BITS - 4);
            for i in 4..BITS {
                v[i] = v[i - 4]
                    ^ (v[i - 4] >> 4)
                    ^ (v[i - 3] << 1)
                    ^ (v[i - 3] >> 2)
                    ^ v[i - 1]
                    ^ (v[i - 1] >> 1);
            }
            v
        }
        7 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 1u32 << (BITS - 2);
            v[2] = 5u32 << (BITS - 4);
            v[3] = 11u32 << (BITS - 5);
            v[4] = 19u32 << (BITS - 6);
            for i in 5..BITS {
                v[i] = v[i - 5]
                    ^ (v[i - 5] >> 5)
                    ^ v[i - 2]
                    ^ (v[i - 2] >> 2)
                    ^ v[i - 1]
                    ^ (v[i - 1] >> 1);
            }
            v
        }
        8 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 3u32 << (BITS - 3);
            v[2] = 7u32 << (BITS - 4);
            v[3] = 5u32 << (BITS - 4);
            v[4] = 21u32 << (BITS - 6);
            for i in 5..BITS {
                v[i] = v[i - 5]
                    ^ (v[i - 5] >> 5)
                    ^ (v[i - 3] << 1)
                    ^ (v[i - 3] >> 2)
                    ^ v[i - 1]
                    ^ (v[i - 1] >> 1);
            }
            v
        }
        9 => {
            let mut v = vec![0u32; BITS];
            v[0] = 1u32 << (BITS - 1);
            v[1] = 1u32 << (BITS - 2);
            v[2] = 1u32 << (BITS - 3);
            v[3] = 7u32 << (BITS - 5);
            v[4] = 13u32 << (BITS - 6);
            v[5] = 37u32 << (BITS - 7);
            for i in 6..BITS {
                v[i] = v[i - 6] ^ (v[i - 6] >> 6) ^ v[i - 1] ^ (v[i - 1] >> 1);
            }
            v
        }
        _ => (0..BITS).map(|i| 1u32 << (BITS - 1 - i)).collect(),
    }
}
/// Generate `n` Sobol points in 10 dimensions.
///
/// Each point is a `Vec`f64` of length 10 with values in `\[0, 1)`.
pub fn sobol_sequence_10d(n: usize) -> Vec<Vec<f64>> {
    pub(super) const BITS: usize = 32;
    pub(super) const N_DIMS: usize = 10;
    let dirs: Vec<Vec<u32>> = (0..N_DIMS).map(sobol10_direction_numbers).collect();
    (0..n)
        .map(|idx| {
            let gray = idx ^ (idx >> 1);
            (0..N_DIMS)
                .map(|d| {
                    let x = (0..BITS)
                        .filter(|&i| (gray >> i) & 1 == 1)
                        .fold(0u32, |acc, i| acc ^ dirs[d][i]);
                    x as f64 / (1u64 << BITS) as f64
                })
                .collect()
        })
        .collect()
}
/// Monte Carlo integration of `f` over `\[a, b\]` using `n` uniform samples.
///
/// Returns `(estimate, standard_error)`.
///
/// The estimate converges as O(1/√n) to ∫_a^b f(x) dx.
pub fn monte_carlo_integrate(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    n: usize,
    rng: &mut Lcg,
) -> (f64, f64) {
    let width = b - a;
    let mut sum = 0.0_f64;
    let mut sum2 = 0.0_f64;
    for _ in 0..n {
        let x = rng.next_f64_range(a, b);
        let v = f(x);
        sum += v;
        sum2 += v * v;
    }
    let n_f = n as f64;
    let mean = sum / n_f;
    let variance = (sum2 / n_f - mean * mean).max(0.0);
    let stderr = (variance / n_f).sqrt() * width;
    (mean * width, stderr)
}
/// Monte Carlo integration of `f` over the unit hypercube `\[0,1)^d` in `d` dimensions.
///
/// Returns `(estimate, standard_error)`.
pub fn monte_carlo_integrate_nd(
    f: impl Fn(&[f64]) -> f64,
    d: usize,
    n: usize,
    rng: &mut Lcg,
) -> (f64, f64) {
    let mut sum = 0.0_f64;
    let mut sum2 = 0.0_f64;
    let mut pt = vec![0.0_f64; d];
    for _ in 0..n {
        for x in pt.iter_mut() {
            *x = rng.next_f64();
        }
        let v = f(&pt);
        sum += v;
        sum2 += v * v;
    }
    let n_f = n as f64;
    let mean = sum / n_f;
    let variance = (sum2 / n_f - mean * mean).max(0.0);
    (mean, (variance / n_f).sqrt())
}
/// Quasi-Monte Carlo integration using the Halton sequence.
///
/// Integrates `f` over `\[a, b\]` using `n` Halton (van der Corput base-2) points.
/// Typically O(log(n)/n) convergence for smooth integrands.
pub fn qmc_integrate_halton(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let width = b - a;
    let sum: f64 = (1..=n)
        .map(|i| {
            let t = HaltonSequence::van_der_corput(i as u32, 2);
            f(a + t * width)
        })
        .sum();
    sum / n as f64 * width
}
/// Quasi-Monte Carlo integration using the 1-D Sobol sequence.
///
/// Integrates `f` over `\[a, b\]` using `n` Sobol points.
pub fn qmc_integrate_sobol(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let width = b - a;
    let mut sobol = Sobol::new_1d();
    let _ = sobol.next_1d();
    let sum: f64 = (0..n).map(|_| f(a + sobol.next_1d() * width)).sum();
    sum / n as f64 * width
}
/// Control-variate Monte Carlo estimator.
///
/// Estimates `E\[f(X)\]` using control variate `g` whose mean `E\[g(X)\]` is known.
///
/// Returns the reduced-variance estimate.
pub fn mc_control_variate(
    f: impl Fn(f64) -> f64,
    g: impl Fn(f64) -> f64,
    g_mean: f64,
    a: f64,
    b: f64,
    n: usize,
    rng: &mut Lcg,
) -> f64 {
    let width = b - a;
    let mut f_vals = Vec::with_capacity(n);
    let mut g_vals = Vec::with_capacity(n);
    for _ in 0..n {
        let x = rng.next_f64_range(a, b);
        f_vals.push(f(x));
        g_vals.push(g(x));
    }
    let n_f = n as f64;
    let f_mean = f_vals.iter().sum::<f64>() / n_f;
    let g_sample_mean = g_vals.iter().sum::<f64>() / n_f;
    let mut cov_fg = 0.0_f64;
    let mut var_g = 0.0_f64;
    for (&fi, &gi) in f_vals.iter().zip(g_vals.iter()) {
        cov_fg += (fi - f_mean) * (gi - g_sample_mean);
        var_g += (gi - g_sample_mean) * (gi - g_sample_mean);
    }
    let c = if var_g > 1e-30 { cov_fg / var_g } else { 0.0 };
    (f_mean - c * (g_sample_mean - g_mean)) * width
}
/// Antithetic variates Monte Carlo estimator.
///
/// For each sample `u`, also evaluates at `1-u`, halving the number of
/// independent evaluations needed for the same variance reduction.
pub fn mc_antithetic(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize, rng: &mut Lcg) -> f64 {
    let width = b - a;
    let m = n.div_ceil(2);
    let sum: f64 = (0..m)
        .map(|_| {
            let u = rng.next_f64();
            let x1 = a + u * width;
            let x2 = a + (1.0 - u) * width;
            (f(x1) + f(x2)) * 0.5
        })
        .sum();
    sum / m as f64 * width
}
/// Compute unnormalized importance weights for a set of samples.
///
/// Given samples drawn from `proposal(x)`, returns the importance weight
/// `target(x) / proposal(x)` for each sample.
pub fn importance_weights(
    samples: &[f64],
    target: impl Fn(f64) -> f64,
    proposal: impl Fn(f64) -> f64,
) -> Vec<f64> {
    samples
        .iter()
        .map(|&x| {
            let p = proposal(x);
            if p > 1e-300 { target(x) / p } else { 0.0 }
        })
        .collect()
}
/// Self-normalized importance sampling estimator for `E_target\[h(X)\]`.
///
/// Samples are drawn from `proposal`; returns the weighted mean of `h`.
pub fn self_normalized_is(
    samples: &[f64],
    h: impl Fn(f64) -> f64,
    target: impl Fn(f64) -> f64,
    proposal: impl Fn(f64) -> f64,
) -> f64 {
    let mut sum_w = 0.0_f64;
    let mut sum_wh = 0.0_f64;
    for &x in samples {
        let p = proposal(x);
        if p < 1e-300 {
            continue;
        }
        let w = target(x) / p;
        sum_w += w;
        sum_wh += w * h(x);
    }
    if sum_w < 1e-300 { 0.0 } else { sum_wh / sum_w }
}
/// Effective sample size (ESS) from normalized importance weights.
///
/// ESS = (Σ wᵢ)² / Σ wᵢ²
pub fn effective_sample_size(weights: &[f64]) -> f64 {
    let sum_w: f64 = weights.iter().sum();
    if sum_w < 1e-300 {
        return 0.0;
    }
    let sum_w2: f64 = weights.iter().map(|&w| (w / sum_w) * (w / sum_w)).sum();
    if sum_w2 < 1e-300 {
        weights.len() as f64
    } else {
        1.0 / sum_w2
    }
}
/// Generate `n_samples` stratified samples over the unit hypercube `\[0,1)^d`.
///
/// This is the Latin Hypercube design: for each dimension, the `n_samples`
/// strata are shuffled independently so that every stratum is visited exactly
/// once per dimension.  Identical to [`LatinHypercube::sample`] but as a
/// free function.
pub fn stratified_unit_hypercube(n_samples: usize, d: usize, rng: &mut Lcg) -> Vec<Vec<f64>> {
    LatinHypercube::new(n_samples, d).sample(rng)
}
/// Latin Hypercube Sampling with maximin distance optimization.
///
/// Generates `n_candidates` LHS designs and returns the one whose minimum
/// inter-point distance (maximin criterion) is largest.
///
/// # Arguments
/// * `n_samples`    – number of points per design
/// * `n_dims`       – number of dimensions
/// * `n_candidates` – how many candidate designs to evaluate
/// * `rng`          – random number generator
pub fn lhs_maximin(
    n_samples: usize,
    n_dims: usize,
    n_candidates: usize,
    rng: &mut Lcg,
) -> Vec<Vec<f64>> {
    let lhs = LatinHypercube::new(n_samples, n_dims);
    let mut best: Vec<Vec<f64>> = lhs.sample(rng);
    let mut best_min_dist = min_pairwise_dist(&best);
    for _ in 1..n_candidates {
        let candidate = lhs.sample(rng);
        let d = min_pairwise_dist(&candidate);
        if d > best_min_dist {
            best_min_dist = d;
            best = candidate;
        }
    }
    best
}
/// Compute the minimum pairwise Euclidean distance among a set of points.
pub(super) fn min_pairwise_dist(pts: &[Vec<f64>]) -> f64 {
    let n = pts.len();
    if n < 2 {
        return f64::INFINITY;
    }
    let mut min_d2 = f64::INFINITY;
    for i in 0..n {
        for j in (i + 1)..n {
            let d2: f64 = pts[i]
                .iter()
                .zip(pts[j].iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum();
            if d2 < min_d2 {
                min_d2 = d2;
            }
        }
    }
    min_d2.sqrt()
}
/// Bootstrap resample `data` with replacement, returning `n` resampled datasets.
///
/// Each dataset is a vector of length `data.len()` drawn uniformly with
/// replacement from `data`.
pub fn bootstrap_resample(data: &[f64], n_resamples: usize, rng: &mut Lcg) -> Vec<Vec<f64>> {
    let m = data.len();
    (0..n_resamples)
        .map(|_| {
            (0..m)
                .map(|_| {
                    let idx = (rng.next_u64() as usize) % m;
                    data[idx]
                })
                .collect()
        })
        .collect()
}
/// Bootstrap estimate of the mean and its 95% confidence interval.
///
/// Returns `(mean, ci_low, ci_high)`.
pub fn bootstrap_mean_ci(data: &[f64], n_resamples: usize, rng: &mut Lcg) -> (f64, f64, f64) {
    let resamples = bootstrap_resample(data, n_resamples, rng);
    let mut means: Vec<f64> = resamples
        .iter()
        .map(|s| s.iter().sum::<f64>() / s.len() as f64)
        .collect();
    means.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    let lo_idx = (0.025 * n_resamples as f64) as usize;
    let hi_idx = ((0.975 * n_resamples as f64) as usize).min(n_resamples - 1);
    (mean, means[lo_idx], means[hi_idx])
}
/// Draw `n` samples from an arbitrary 1-D PDF using rejection sampling.
///
/// `pdf` must be non-negative; `envelope` must satisfy `pdf(x) ≤ envelope`
/// for all `x ∈ \[lo, hi)`.
///
/// # Arguments
/// * `pdf`      – target density function
/// * `lo`, `hi` – support of the distribution
/// * `envelope` – upper bound of `pdf` over the support
/// * `n`        – number of samples to draw
/// * `rng`      – LCG random source
pub fn rejection_sample_1d(
    pdf: impl Fn(f64) -> f64,
    lo: f64,
    hi: f64,
    envelope: f64,
    n: usize,
    rng: &mut Lcg,
) -> Vec<f64> {
    let mut samples = Vec::with_capacity(n);
    while samples.len() < n {
        let x = rng.next_f64_range(lo, hi);
        let u = rng.next_f64() * envelope;
        if u < pdf(x) {
            samples.push(x);
        }
    }
    samples
}
/// Compute IS estimate of `E\[h(X)\]` where `X ~ target` but samples are drawn
/// from `proposal`.
///
/// Returns `(estimate, effective_sample_size)`.
///
/// `weights\[i\] = target(x_i) / proposal(x_i)`.
pub fn importance_sampling_estimate(
    samples: &[f64],
    h: impl Fn(f64) -> f64,
    target: impl Fn(f64) -> f64,
    proposal: impl Fn(f64) -> f64,
) -> (f64, f64) {
    let mut sum_w = 0.0_f64;
    let mut sum_wh = 0.0_f64;
    let mut sum_w2 = 0.0_f64;
    for &x in samples {
        let p = proposal(x);
        if p < 1e-300 {
            continue;
        }
        let w = target(x) / p;
        sum_w += w;
        sum_wh += w * h(x);
        sum_w2 += (w / 1.0).powi(2);
    }
    let estimate = if sum_w > 1e-300 { sum_wh / sum_w } else { 0.0 };
    let ess = if sum_w2 > 1e-300 {
        (sum_w * sum_w) / sum_w2
    } else {
        0.0
    };
    (estimate, ess)
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::sampling::GibbsSampler;

    use crate::sampling::ImportanceSampler;
    use crate::sampling::Lcg;
    use crate::sampling::MetropolisHastings;
    use crate::sampling::StratifiedSampler;

    /// next_f64 must lie in [0, 1).
    #[test]
    fn test_next_f64_range() {
        let mut rng = Lcg::new(42);
        for _ in 0..10_000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v), "out of range: {v}");
        }
    }
    /// Sample 1 000 standard normals; empirical mean must be close to 0.
    #[test]
    fn test_next_normal_mean() {
        let mut rng = Lcg::new(123);
        let n = 1_000;
        let mean: f64 = (0..n).map(|_| rng.next_normal()).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.15, "mean {mean} too far from 0");
    }
    /// All points on the sphere surface must have |v| ≈ 1.
    #[test]
    fn test_sample_sphere_surface_norm() {
        let mut rng = Lcg::new(7);
        for _ in 0..1_000 {
            let [x, y, z] = sample_sphere_surface(&mut rng);
            let r = (x * x + y * y + z * z).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "sphere surface radius = {r}");
        }
    }
    /// All points inside the ball must have |v| ≤ 1.
    #[test]
    fn test_sample_sphere_volume_inside() {
        let mut rng = Lcg::new(99);
        for _ in 0..1_000 {
            let [x, y, z] = sample_sphere_volume(&mut rng);
            let r2 = x * x + y * y + z * z;
            assert!(r2 <= 1.0 + 1e-12, "point outside unit ball: r² = {r2}");
        }
    }
    /// Stratified samples: every stratum `\[k/n, (k+1)/n)` must be occupied.
    #[test]
    fn test_stratified_sample_1d_coverage() {
        let mut rng = Lcg::new(55);
        let n = 20;
        let samples = stratified_sample_1d(n, &mut rng);
        assert_eq!(samples.len(), n);
        let inv_n = 1.0 / n as f64;
        for (k, &s) in samples.iter().enumerate() {
            let lo = k as f64 * inv_n;
            let hi = lo + inv_n;
            assert!(
                s >= lo && s < hi,
                "stratum {k}: sample {s} not in [{lo}, {hi})"
            );
        }
    }
    #[test]
    fn test_lhs_sample_count() {
        let mut rng = Lcg::new(100);
        let lhs = LatinHypercube::new(10, 3);
        let samples = lhs.sample(&mut rng);
        assert_eq!(samples.len(), 10);
        for s in &samples {
            assert_eq!(s.len(), 3);
        }
    }
    #[test]
    fn test_lhs_covers_all_strata() {
        let mut rng = Lcg::new(200);
        let n = 8;
        let lhs = LatinHypercube::new(n, 2);
        let samples = lhs.sample(&mut rng);
        let inv_n = 1.0 / n as f64;
        for d in 0..2 {
            let mut strata = vec![false; n];
            for s in &samples {
                let stratum = (s[d] / inv_n) as usize;
                let stratum = stratum.min(n - 1);
                strata[stratum] = true;
            }
            for (k, &occupied) in strata.iter().enumerate() {
                assert!(occupied, "dim {d} stratum {k} not occupied");
            }
        }
    }
    #[test]
    fn test_stratified_2d_count() {
        let mut rng = Lcg::new(300);
        let ss = StratifiedSampler::new(5);
        let samples = ss.sample_2d(&mut rng);
        assert_eq!(samples.len(), 25);
    }
    #[test]
    fn test_stratified_2d_coverage() {
        let mut rng = Lcg::new(400);
        let n = 4;
        let ss = StratifiedSampler::new(n);
        let samples = ss.sample_2d(&mut rng);
        let inv_n = 1.0 / n as f64;
        for (idx, s) in samples.iter().enumerate() {
            let ix = idx % n;
            let iy = idx / n;
            let x_lo = ix as f64 * inv_n;
            let y_lo = iy as f64 * inv_n;
            assert!(s[0] >= x_lo && s[0] < x_lo + inv_n, "x out of stratum");
            assert!(s[1] >= y_lo && s[1] < y_lo + inv_n, "y out of stratum");
        }
    }
    #[test]
    fn test_halton_in_unit_interval() {
        let samples = HaltonSequence::sample(100, 2);
        for &s in &samples {
            assert!((0.0..1.0).contains(&s), "Halton sample {s} out of [0,1)");
        }
    }
    #[test]
    fn test_halton_base2_first_values() {
        let samples = HaltonSequence::sample(4, 2);
        assert!((samples[0] - 0.5).abs() < 1e-12);
        assert!((samples[1] - 0.25).abs() < 1e-12);
        assert!((samples[2] - 0.75).abs() < 1e-12);
        assert!((samples[3] - 0.125).abs() < 1e-12);
    }
    #[test]
    fn test_sobol_multi_dim_count() {
        let sobol = SobolSequence::new(3);
        let samples = sobol.sample(16);
        assert_eq!(samples.len(), 16);
        for s in &samples {
            assert_eq!(s.len(), 3);
        }
    }
    #[test]
    fn test_sobol_in_unit_interval() {
        let sobol = SobolSequence::new(2);
        let samples = sobol.sample(64);
        for s in &samples {
            for &v in s {
                assert!((0.0..1.0).contains(&v), "Sobol value {v} out of [0,1)");
            }
        }
    }
    #[test]
    fn test_importance_sampler_generates_samples() {
        pub(super) fn uniform_pdf(_x: f64) -> f64 {
            1.0
        }
        let sampler = ImportanceSampler::new(uniform_pdf, 1.0);
        let mut rng = Lcg::new(500);
        let samples = sampler.sample(50, &mut rng);
        assert_eq!(samples.len(), 50);
        for &s in &samples {
            assert!((0.0..1.0).contains(&s));
        }
    }
    #[test]
    fn test_importance_sampler_biased_pdf() {
        pub(super) fn biased_pdf(x: f64) -> f64 {
            if x >= 0.5 { 2.0 } else { 0.0 }
        }
        let sampler = ImportanceSampler::new(biased_pdf, 2.0);
        let mut rng = Lcg::new(600);
        let samples = sampler.sample(100, &mut rng);
        for &s in &samples {
            assert!(s >= 0.5, "all samples should be >= 0.5, got {s}");
        }
    }
    #[test]
    fn test_metropolis_hastings_stays_in_range() {
        pub(super) fn log_target(x: f64) -> f64 {
            -0.5 * x * x
        }
        let mut mcmc = MetropolisHastings::new(0.0, 1.0, 42);
        let samples = mcmc.sample(500, log_target);
        let n_outliers = samples.iter().filter(|&&x| x.abs() > 5.0).count();
        assert!(n_outliers < 20, "{} outliers > 5σ", n_outliers);
    }
    #[test]
    fn test_metropolis_hastings_mean_near_zero() {
        pub(super) fn log_target(x: f64) -> f64 {
            -0.5 * x * x
        }
        let mut mcmc = MetropolisHastings::new(0.0, 0.5, 99);
        let samples = mcmc.sample(2000, log_target);
        let mean: f64 = samples[500..].iter().sum::<f64>() / (samples.len() - 500) as f64;
        assert!(mean.abs() < 0.5, "mean={}", mean);
    }
    #[test]
    fn test_gibbs_sampler_basic() {
        let mut gibbs = GibbsSampler::new(2, 42);
        let samples = gibbs.sample(
            200,
            &[
                Box::new(|_, _: &[f64]| 0.0_f64),
                Box::new(|_, _: &[f64]| 0.0_f64),
            ],
            &[1.0, 1.0],
        );
        assert_eq!(samples.len(), 200);
        for s in &samples {
            assert_eq!(s.len(), 2);
        }
    }
    #[test]
    fn test_blue_noise_count() {
        let mut rng = Lcg::new(777);
        let samples = blue_noise_2d(16, 0.1, &mut rng);
        assert_eq!(samples.len(), 16);
    }
    #[test]
    fn test_blue_noise_min_distance() {
        let mut rng = Lcg::new(888);
        let r = 0.15;
        let n = 20;
        let samples = blue_noise_2d(n, r, &mut rng);
        for i in 0..samples.len() {
            for j in (i + 1)..samples.len() {
                let dx = samples[i][0] - samples[j][0];
                let dy = samples[i][1] - samples[j][1];
                let d = (dx * dx + dy * dy).sqrt();
                assert!(d >= r - 1e-9, "samples {} and {} too close: d={}", i, j, d);
            }
        }
    }
    #[test]
    fn test_stratified_nd_count() {
        let mut rng = Lcg::new(999);
        let samples = stratified_sample_nd(3, 4, &mut rng);
        assert_eq!(samples.len(), 3);
        for s in &samples {
            assert_eq!(s.len(), 4);
        }
    }
    #[test]
    fn test_stratified_nd_in_range() {
        let mut rng = Lcg::new(1001);
        let samples = stratified_sample_nd(10, 3, &mut rng);
        for s in &samples {
            for &v in s {
                assert!((0.0..=1.0).contains(&v), "value {} out of [0,1]", v);
            }
        }
    }
    #[test]
    fn test_sobol_1d_low_discrepancy() {
        let mut sobol = Sobol::new_1d();
        let vals: Vec<f64> = (0..16).map(|_| sobol.next_1d()).collect();
        for &v in &vals {
            assert!((0.0..1.0).contains(&v), "Sobol val {} out of range", v);
        }
    }
    #[test]
    fn test_lhs_rand_count_and_dims() {
        let mut rng = rand::rng();
        let samples = latin_hypercube_sample(12, 4, &mut rng);
        assert_eq!(samples.len(), 12);
        for s in &samples {
            assert_eq!(s.len(), 4);
        }
    }
    #[test]
    fn test_sobol_sequence_equidistributed_dim0() {
        let samples = sobol_sequence(16, 1);
        assert_eq!(samples.len(), 16);
        let mut sorted: Vec<f64> = samples.iter().map(|s| s[0]).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for &v in &sorted {
            assert!((0.0..1.0).contains(&v), "sobol val {} out of range", v);
        }
    }
    #[test]
    fn test_stratified_sample_1d_rng_covers_unit_interval() {
        let mut rng = rand::rng();
        let n = 10;
        let samples = stratified_sample_1d_rng(n, &mut rng);
        assert_eq!(samples.len(), n);
        for &s in &samples {
            assert!((0.0..=1.0).contains(&s), "value {} out of [0,1]", s);
        }
    }
    #[test]
    fn test_halton_sequence_in_unit_interval() {
        let samples = halton_sequence(50, 2);
        assert_eq!(samples.len(), 50);
        let sum: f64 = samples.iter().sum();
        assert!(
            sum > 0.0 && sum < 50.0,
            "Halton sum {} out of expected range",
            sum
        );
        for &s in &samples {
            assert!((0.0..1.0).contains(&s), "Halton value {} out of [0,1)", s);
        }
    }
    #[test]
    fn test_importance_sample_count() {
        let mut rng = rand::rng();
        let pdf = vec![1.0, 2.0, 3.0, 4.0];
        let samples = importance_sample(&pdf, 100, &mut rng);
        assert_eq!(samples.len(), 100);
        for &idx in &samples {
            assert!(idx < pdf.len(), "index {} out of range", idx);
        }
    }
    #[test]
    fn test_importance_sample_biased_toward_last() {
        let mut rng = rand::rng();
        let pdf = vec![0.001, 0.001, 0.001, 100.0];
        let samples = importance_sample(&pdf, 200, &mut rng);
        let n_last = samples.iter().filter(|&&i| i == 3).count();
        assert!(
            n_last > 150,
            "expected mostly index 3, got {} out of 200",
            n_last
        );
    }
    #[test]
    fn test_rejection_sample_2d_count() {
        let mut rng = rand::rng();
        let samples = rejection_sample_2d(|_x, _y| 1.0, (0.0, 1.0), (0.0, 1.0), 20, &mut rng);
        assert_eq!(samples.len(), 20);
        for s in &samples {
            assert!(s[0] >= 0.0 && s[0] < 1.0, "x={} out of range", s[0]);
            assert!(s[1] >= 0.0 && s[1] < 1.0, "y={} out of range", s[1]);
        }
    }
    #[test]
    fn test_weighted_resample_count() {
        let particles: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let weights: Vec<f64> = vec![1.0; 10];
        let mut rng = Lcg::new(42);
        let resampled = systematic_resample(&particles, &weights, 10, &mut rng);
        assert_eq!(resampled.len(), 10);
    }
    #[test]
    fn test_weighted_resample_biased() {
        let particles: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let mut weights = vec![0.001; 10];
        weights[5] = 100.0;
        let mut rng = Lcg::new(42);
        let resampled = systematic_resample(&particles, &weights, 20, &mut rng);
        let n_fives = resampled
            .iter()
            .filter(|&&v| (v - 5.0).abs() < 1e-9)
            .count();
        assert!(
            n_fives > 15,
            "most resampled should be 5.0, got {}",
            n_fives
        );
    }
    #[test]
    fn test_sobol10_in_unit_interval() {
        let pts = sobol_sequence_10d(64);
        assert_eq!(pts.len(), 64);
        for p in &pts {
            assert_eq!(p.len(), 10);
            for &v in p {
                assert!((0.0..1.0).contains(&v), "Sobol10 value {} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_sobol10_dim0_matches_1d() {
        let pts10 = sobol_sequence_10d(16);
        let pts1 = sobol_sequence(16, 1);
        for (p10, p1) in pts10.iter().zip(pts1.iter()) {
            assert!((p10[0] - p1[0]).abs() < 1e-14, "dim0 mismatch");
        }
    }
    #[test]
    fn test_monte_carlo_integrate_constant() {
        let mut rng = Lcg::new(42);
        let (est, err) = monte_carlo_integrate(|_| 1.0, 0.0, 1.0, 10_000, &mut rng);
        assert!((est - 1.0).abs() < 0.01, "estimate={}", est);
        assert!(err < 0.01, "error={}", err);
    }
    #[test]
    fn test_monte_carlo_integrate_linear() {
        let mut rng = Lcg::new(7);
        let (est, _err) = monte_carlo_integrate(|x| x, 0.0, 1.0, 50_000, &mut rng);
        assert!((est - 0.5).abs() < 0.01, "estimate={}", est);
    }
    #[test]
    fn test_monte_carlo_integrate_quadratic() {
        let mut rng = Lcg::new(13);
        let (est, _err) = monte_carlo_integrate(|x| x * x, 0.0, 1.0, 100_000, &mut rng);
        assert!((est - 1.0 / 3.0).abs() < 0.01, "estimate={}", est);
    }
    #[test]
    fn test_qmc_halton_vs_mc_variance() {
        let qmc_est = qmc_integrate_halton(|x| x * x, 0.0, 1.0, 1024);
        assert!(
            (qmc_est - 1.0 / 3.0).abs() < 0.005,
            "QMC estimate={}",
            qmc_est
        );
    }
    #[test]
    fn test_qmc_sobol_vs_mc_variance() {
        let est = qmc_integrate_sobol(|x| x, 0.0, 1.0, 1024);
        assert!((est - 0.5).abs() < 0.005, "Sobol QMC estimate={}", est);
    }
    #[test]
    fn test_importance_weight_estimator_uniform() {
        let samples: Vec<f64> = (0..10).map(|i| (i as f64 + 0.5) / 10.0).collect();
        let weights = importance_weights(&samples, |_| 1.0, |_| 1.0);
        for &w in &weights {
            assert!((w - 1.0).abs() < 1e-12, "weight={}", w);
        }
    }
    #[test]
    fn test_importance_weight_estimator_ratio() {
        let samples = vec![0.25, 0.5, 0.75];
        let weights = importance_weights(&samples, |x| 2.0 * x, |_| 1.0);
        assert!((weights[0] - 0.5).abs() < 1e-12);
        assert!((weights[1] - 1.0).abs() < 1e-12);
        assert!((weights[2] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_stratified_hypercube_all_in_range() {
        let mut rng = Lcg::new(2025);
        let pts = stratified_unit_hypercube(8, 5, &mut rng);
        assert_eq!(pts.len(), 8);
        for p in &pts {
            assert_eq!(p.len(), 5);
            for &v in p {
                assert!((0.0..1.0).contains(&v), "v={} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_stratified_hypercube_coverage() {
        let mut rng = Lcg::new(3030);
        let n = 6;
        let pts = stratified_unit_hypercube(n, 2, &mut rng);
        let inv = 1.0 / n as f64;
        for d in 0..2 {
            let mut occupied = vec![false; n];
            for p in &pts {
                let s = (p[d] / inv) as usize;
                occupied[s.min(n - 1)] = true;
            }
            for (k, &hit) in occupied.iter().enumerate() {
                assert!(hit, "dim {d} stratum {k} not covered");
            }
        }
    }
}
#[cfg(test)]
mod extended_sampling_tests {
    use crate::sampling::GaussianKde;

    use crate::sampling::HammersleySequence;

    use crate::sampling::Lcg;

    use crate::sampling::bootstrap_mean_ci;
    use crate::sampling::bootstrap_resample;
    use crate::sampling::importance_sampling_estimate;
    use crate::sampling::lhs_maximin;
    use crate::sampling::min_pairwise_dist;
    use crate::sampling::rejection_sample_1d;
    #[test]
    fn test_lhs_maximin_returns_correct_shape() {
        let mut rng = Lcg::new(1111);
        let pts = lhs_maximin(8, 3, 5, &mut rng);
        assert_eq!(pts.len(), 8);
        for p in &pts {
            assert_eq!(p.len(), 3);
            for &v in p {
                assert!((0.0..1.0).contains(&v), "v={} out of [0,1)", v);
            }
        }
    }
    #[test]
    fn test_lhs_maximin_better_than_single_lhs() {
        let mut rng = Lcg::new(2222);
        let optimized = lhs_maximin(10, 2, 20, &mut rng);
        let d_opt = min_pairwise_dist(&optimized);
        assert!(d_opt > 0.0, "min dist should be positive");
    }
    #[test]
    fn test_lhs_maximin_single_candidate_equals_lhs() {
        let mut rng = Lcg::new(42);
        let pts = lhs_maximin(5, 2, 1, &mut rng);
        let inv = 1.0 / 5.0_f64;
        for d in 0..2 {
            let mut occupied = [false; 5];
            for p in &pts {
                let s = (p[d] / inv) as usize;
                occupied[s.min(4)] = true;
            }
            for (k, &hit) in occupied.iter().enumerate() {
                assert!(
                    hit,
                    "dim {d} stratum {k} not covered in single-candidate LHS"
                );
            }
        }
    }
    #[test]
    fn test_hammersley_shape() {
        let pts = HammersleySequence::sample(16, 3);
        assert_eq!(pts.len(), 16);
        for p in &pts {
            assert_eq!(p.len(), 3);
            for &v in p {
                assert!((0.0..=1.0).contains(&v), "v={} out of [0,1]", v);
            }
        }
    }
    #[test]
    fn test_hammersley_first_dim_uniform_grid() {
        let n = 8;
        let pts = HammersleySequence::sample(n, 2);
        for (i, p) in pts.iter().enumerate() {
            let expected = i as f64 / n as f64;
            assert!((p[0] - expected).abs() < 1e-12, "dim0[{}]={}", i, p[0]);
        }
    }
    #[test]
    fn test_hammersley_second_dim_halton_base2() {
        let pts = HammersleySequence::sample(4, 2);
        let expected = [0.5, 0.25, 0.75, 0.125];
        for (i, &exp) in expected.iter().enumerate() {
            assert!(
                (pts[i][1] - exp).abs() < 1e-12,
                "dim1[{}]={} expected={}",
                i,
                pts[i][1],
                exp
            );
        }
    }
    #[test]
    fn test_bootstrap_resample_shape() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let mut rng = Lcg::new(9999);
        let resamples = bootstrap_resample(&data, 50, &mut rng);
        assert_eq!(resamples.len(), 50);
        for r in &resamples {
            assert_eq!(r.len(), 10);
        }
    }
    #[test]
    fn test_bootstrap_resample_values_from_data() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut rng = Lcg::new(7777);
        let resamples = bootstrap_resample(&data, 20, &mut rng);
        for r in &resamples {
            for &v in r {
                assert!(
                    data.contains(&v),
                    "bootstrap value {} not in original data",
                    v
                );
            }
        }
    }
    #[test]
    fn test_bootstrap_mean_ci_contains_true_mean() {
        let mut rng = Lcg::new(5555);
        let data: Vec<f64> = (0..200).map(|_| rng.next_f64()).collect();
        let (mean, ci_lo, ci_hi) = bootstrap_mean_ci(&data, 500, &mut rng);
        assert!(ci_lo <= mean, "ci_lo={} > mean={}", ci_lo, mean);
        assert!(ci_hi >= mean, "ci_hi={} < mean={}", ci_hi, mean);
        assert!(
            ci_lo < 0.5 && ci_hi > 0.5,
            "CI [{}, {}] does not contain 0.5",
            ci_lo,
            ci_hi
        );
    }
    #[test]
    fn test_kde_positive_density() {
        let data = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let kde = GaussianKde::new(data, 0.5);
        for &x in &[0.0_f64, 0.5, 1.0, 1.5, 2.0] {
            let d = kde.evaluate(x);
            assert!(d > 0.0, "KDE density at {} should be positive: {}", x, d);
        }
    }
    #[test]
    fn test_kde_integrates_to_one() {
        let data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let kde = GaussianKde::new(data, 0.5);
        let n = 10_000;
        let lo = -5.0_f64;
        let hi = 9.0_f64;
        let width = (hi - lo) / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let x = lo + (i as f64 + 0.5) * width;
                kde.evaluate(x) * width
            })
            .sum();
        assert!((integral - 1.0).abs() < 0.05, "KDE integral={}", integral);
    }
    #[test]
    fn test_kde_silverman_bandwidth_positive() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let h = GaussianKde::silverman_bandwidth(&data);
        assert!(h > 0.0, "bandwidth should be positive: {}", h);
    }
    #[test]
    fn test_kde_auto_bandwidth() {
        let data = vec![0.0, 1.0, 2.0, 3.0];
        let kde = GaussianKde::new(data, -1.0);
        assert!(kde.bandwidth > 0.0, "auto bandwidth should be positive");
    }
    #[test]
    fn test_kde_evaluate_batch() {
        let data = vec![0.0, 1.0, 2.0];
        let kde = GaussianKde::new(data, 0.3);
        let queries = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let batch = kde.evaluate_batch(&queries);
        assert_eq!(batch.len(), 5);
        for &v in &batch {
            assert!(v > 0.0, "KDE batch value non-positive: {}", v);
        }
        let kde2 = GaussianKde::new(vec![0.0, 1.0, 2.0], 0.3);
        for (i, &q) in queries.iter().enumerate() {
            assert!((batch[i] - kde2.evaluate(q)).abs() < 1e-14);
        }
    }
    #[test]
    fn test_rejection_sample_1d_uniform() {
        let mut rng = Lcg::new(1234);
        let samples = rejection_sample_1d(|_| 1.0, 0.0, 1.0, 1.0, 50, &mut rng);
        assert_eq!(samples.len(), 50);
        for &s in &samples {
            assert!((0.0..1.0).contains(&s), "sample {} out of [0,1)", s);
        }
    }
    #[test]
    fn test_rejection_sample_1d_upper_half_only() {
        let mut rng = Lcg::new(5678);
        let samples = rejection_sample_1d(
            |x| if x >= 0.5 { 2.0 } else { 0.0 },
            0.0,
            1.0,
            2.0,
            30,
            &mut rng,
        );
        assert_eq!(samples.len(), 30);
        for &s in &samples {
            assert!(s >= 0.5, "sample {} should be >= 0.5", s);
        }
    }
    #[test]
    fn test_rejection_sample_1d_range_correct() {
        let mut rng = Lcg::new(4321);
        let samples = rejection_sample_1d(|x| 3.0 * x * x, 0.0, 1.0, 3.0, 100, &mut rng);
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!((mean - 0.75).abs() < 0.15, "mean={}", mean);
    }
    #[test]
    fn test_importance_sampling_estimate_uniform() {
        let samples: Vec<f64> = (0..100).map(|i| (i as f64 + 0.5) / 100.0).collect();
        let (est, _ess) = importance_sampling_estimate(&samples, |x| x, |_| 1.0, |_| 1.0);
        assert!((est - 0.5).abs() < 0.01, "estimate={}", est);
    }
    #[test]
    fn test_importance_sampling_estimate_ess_equals_n_uniform() {
        let samples: Vec<f64> = (0..50).map(|i| i as f64 / 50.0).collect();
        let (_est, ess) = importance_sampling_estimate(&samples, |x| x, |_| 1.0, |_| 1.0);
        assert!((ess - 50.0).abs() < 0.01, "ESS={}", ess);
    }
    #[test]
    fn test_importance_sampling_estimate_reweighted() {
        let mut rng = Lcg::new(8888);
        let samples: Vec<f64> = (0..500).map(|_| rng.next_f64()).collect();
        let (est, ess) = importance_sampling_estimate(&samples, |x| x, |x| 2.0 * x, |_| 1.0);
        assert!(
            (est - 2.0 / 3.0).abs() < 0.05,
            "IS estimate={} expected~2/3",
            est
        );
        assert!(ess > 0.0, "ESS should be positive: {}", ess);
    }
}
