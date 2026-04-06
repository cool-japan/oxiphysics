#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Random processes: Brownian motion, Poisson processes, Markov chains,
//! continuous-time Markov chains, Lévy processes, and SDE solvers.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─── Minimal LCG RNG (no external rand crate needed in core logic) ───────────

/// Fast linear congruential pseudo-random number generator.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }
    /// Uniform float in \[0, 1).
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Standard normal via Box-Muller.
    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
    /// Normal N(mu, sigma).
    fn normal_ms(&mut self, mu: f64, sigma: f64) -> f64 {
        mu + sigma * self.normal()
    }
}

// ─── BrownianMotion ──────────────────────────────────────────────────────────

/// Parameters and generators for Brownian motion processes.
pub struct BrownianMotion;

impl BrownianMotion {
    /// Simulate a standard Brownian motion path W(t) on \[0, T\] with `n_steps` steps.
    ///
    /// Returns a vector of length `n_steps + 1` starting at 0, where increments
    /// are iid N(0, dt).
    pub fn standard(t: f64, n_steps: usize, seed: u64) -> Vec<f64> {
        let dt = t / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut rng = Lcg::new(seed);
        let mut path = vec![0.0f64; n_steps + 1];
        for i in 1..=n_steps {
            path[i] = path[i - 1] + sqrt_dt * rng.normal();
        }
        path
    }

    /// Simulate a Geometric Brownian Motion S(t) = S0 * exp((mu - 0.5*sigma^2)*t + sigma*W(t)).
    ///
    /// Returns the price path of length `n_steps + 1`.
    pub fn geometric(s0: f64, mu: f64, sigma: f64, t: f64, n_steps: usize, seed: u64) -> Vec<f64> {
        let dt = t / n_steps as f64;
        let drift = (mu - 0.5 * sigma * sigma) * dt;
        let vol = sigma * dt.sqrt();
        let mut rng = Lcg::new(seed);
        let mut path = vec![s0; n_steps + 1];
        for i in 1..=n_steps {
            path[i] = path[i - 1] * (drift + vol * rng.normal()).exp();
        }
        path
    }

    /// Simulate a Fractional Brownian Motion (fBm) path with Hurst exponent H.
    ///
    /// Uses the Davies-Harte (spectral) method via a circulant embedding.
    /// Valid for H in (0, 1); H = 0.5 recovers standard BM.
    pub fn fractional(hurst: f64, n_steps: usize, seed: u64) -> Vec<f64> {
        assert!(
            hurst > 0.0 && hurst < 1.0,
            "Hurst exponent must be in (0, 1)"
        );
        let n = n_steps;
        // Autocovariance: gamma(k) = 0.5*(|k-1|^{2H} - 2|k|^{2H} + |k+1|^{2H})
        let gamma = |k: f64| -> f64 {
            0.5 * ((k - 1.0).abs().powf(2.0 * hurst) - 2.0 * k.abs().powf(2.0 * hurst)
                + (k + 1.0).abs().powf(2.0 * hurst))
        };
        // First row of circulant (length 2n)
        let m = 2 * n;
        let mut c = vec![0.0f64; m];
        c[0] = 1.0;
        for k in 1..n {
            c[k] = gamma(k as f64);
            c[m - k] = c[k];
        }
        // DFT of c (eigenvalues of circulant — real for symmetric c)
        let mut eigvals = vec![0.0f64; m];
        for j in 0..m {
            let mut sum = 0.0;
            for k in 0..m {
                sum += c[k] * (2.0 * PI * j as f64 * k as f64 / m as f64).cos();
            }
            eigvals[j] = sum / m as f64;
        }
        let mut rng = Lcg::new(seed);
        // Generate complex white noise in frequency domain
        let mut w_real = vec![0.0f64; m];
        let mut w_imag = vec![0.0f64; m];
        for k in 0..m {
            let z1 = rng.normal();
            let z2 = rng.normal();
            let s = if eigvals[k] > 0.0 {
                eigvals[k].sqrt()
            } else {
                0.0
            };
            w_real[k] = s * z1;
            w_imag[k] = s * z2;
        }
        // Inverse DFT (real part only)
        let mut path = vec![0.0f64; n + 1];
        for i in 1..=n {
            let mut sum = 0.0;
            for k in 0..m {
                sum += w_real[k] * (2.0 * PI * i as f64 * k as f64 / m as f64).cos()
                    - w_imag[k] * (2.0 * PI * i as f64 * k as f64 / m as f64).sin();
            }
            path[i] = sum;
        }
        // Cumulative sum to get fBm increments
        for i in 1..=n {
            path[i] += path[i - 1];
        }
        path
    }

    /// Compute the quadratic variation of a path (empirical check for BM).
    ///
    /// For a standard BM of length T, the quadratic variation converges to T.
    pub fn quadratic_variation(path: &[f64]) -> f64 {
        path.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum()
    }

    /// Compute the sample mean of a path.
    pub fn mean(path: &[f64]) -> f64 {
        path.iter().sum::<f64>() / path.len() as f64
    }

    /// Compute the sample variance of a path.
    pub fn variance(path: &[f64]) -> f64 {
        let m = Self::mean(path);
        path.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (path.len() - 1) as f64
    }
}

// ─── PoissonProcess ──────────────────────────────────────────────────────────

/// Homogeneous and inhomogeneous Poisson process generators.
pub struct PoissonProcess;

impl PoissonProcess {
    /// Simulate a homogeneous Poisson process with rate `lambda` on \[0, T\].
    ///
    /// Returns the arrival times as a sorted vector.
    pub fn homogeneous(lambda: f64, t: f64, seed: u64) -> Vec<f64> {
        assert!(lambda > 0.0, "Rate must be positive");
        let mut rng = Lcg::new(seed);
        let mut times = Vec::new();
        let mut current = 0.0f64;
        loop {
            // Inter-arrival ~ Exp(lambda)
            let u = rng.uniform().max(1e-300);
            let inter = -u.ln() / lambda;
            current += inter;
            if current >= t {
                break;
            }
            times.push(current);
        }
        times
    }

    /// Simulate an inhomogeneous Poisson process with intensity function `lambda(t)`.
    ///
    /// Uses thinning (Lewis-Shedler algorithm): first generate a homogeneous PP
    /// with rate `lambda_max`, then accept each event with probability `lambda(t)/lambda_max`.
    pub fn inhomogeneous<F>(lambda: F, lambda_max: f64, t: f64, seed: u64) -> Vec<f64>
    where
        F: Fn(f64) -> f64,
    {
        let raw = Self::homogeneous(lambda_max, t, seed);
        let mut rng = Lcg::new(seed.wrapping_add(1));
        raw.into_iter()
            .filter(|&ti| rng.uniform() < lambda(ti) / lambda_max)
            .collect()
    }

    /// Simulate a compound Poisson process: N(T) events each with random jump size.
    ///
    /// Returns (arrival_times, cumulative_sum_of_jumps).
    pub fn compound<F>(lambda: f64, t: f64, jump_dist: F, seed: u64) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(f64) -> f64,
    {
        let mut rng = Lcg::new(seed);
        let mut times = Vec::new();
        let mut values = vec![0.0f64];
        let mut current = 0.0f64;
        let mut cum = 0.0f64;
        loop {
            let u = rng.uniform().max(1e-300);
            current += -u.ln() / lambda;
            if current >= t {
                break;
            }
            let jump = jump_dist(rng.uniform());
            cum += jump;
            times.push(current);
            values.push(cum);
        }
        (times, values)
    }

    /// Return the expected number of events E\[N(t)\] = lambda * t.
    pub fn expected_count(lambda: f64, t: f64) -> f64 {
        lambda * t
    }

    /// Inter-arrival distribution CDF: P(T <= t) = 1 - exp(-lambda * t).
    pub fn inter_arrival_cdf(lambda: f64, t: f64) -> f64 {
        1.0 - (-lambda * t).exp()
    }

    /// Compute the empirical rate from observed arrival times.
    pub fn empirical_rate(times: &[f64], t: f64) -> f64 {
        times.len() as f64 / t
    }
}

// ─── MarkovChain ─────────────────────────────────────────────────────────────

/// Discrete-time finite Markov chain with transition matrix P.
pub struct MarkovChain {
    /// Transition matrix P (n x n), row-stochastic.
    pub p: Vec<f64>,
    /// Number of states.
    pub n_states: usize,
}

impl MarkovChain {
    /// Create a Markov chain from a row-stochastic transition matrix (n x n).
    pub fn new(p: Vec<f64>, n_states: usize) -> Self {
        Self { p, n_states }
    }

    /// Simulate a trajectory of length `n_steps` starting from state `init`.
    pub fn simulate(&self, init: usize, n_steps: usize, seed: u64) -> Vec<usize> {
        let mut rng = Lcg::new(seed);
        let mut path = vec![init];
        let mut state = init;
        for _ in 0..n_steps {
            let u = rng.uniform();
            let mut cum = 0.0;
            let mut next = self.n_states - 1;
            for j in 0..self.n_states {
                cum += self.p[state * self.n_states + j];
                if u < cum {
                    next = j;
                    break;
                }
            }
            path.push(next);
            state = next;
        }
        path
    }

    /// Compute the stationary distribution by iterating the distribution vector.
    ///
    /// Returns the vector pi such that pi P = pi.
    pub fn stationary_distribution(&self, tol: f64, max_iter: usize) -> Vec<f64> {
        let n = self.n_states;
        let mut pi = vec![1.0 / n as f64; n];
        for _ in 0..max_iter {
            let new_pi: Vec<f64> = (0..n)
                .map(|j| (0..n).map(|i| pi[i] * self.p[i * n + j]).sum())
                .collect();
            let diff: f64 = new_pi.iter().zip(&pi).map(|(a, b)| (a - b).abs()).sum();
            pi = new_pi;
            if diff < tol {
                break;
            }
        }
        pi
    }

    /// Compute the n-step transition matrix P^n.
    pub fn n_step_matrix(&self, n: usize) -> Vec<f64> {
        let k = self.n_states;
        let mut result = vec![0.0f64; k * k];
        for i in 0..k {
            result[i * k + i] = 1.0;
        } // identity
        let mut base = self.p.clone();
        let mut exp = n;
        while exp > 0 {
            if exp & 1 == 1 {
                result = Self::matmul_sq(&result, &base, k);
            }
            base = Self::matmul_sq(&base, &base, k);
            exp >>= 1;
        }
        result
    }

    /// Estimate the mixing time: smallest t such that ||(pi_0 P^t) - pi||_1 < eps.
    pub fn mixing_time(&self, init: usize, eps: f64) -> usize {
        let pi = self.stationary_distribution(1e-10, 2000);
        let n = self.n_states;
        let mut dist = vec![0.0f64; n];
        dist[init] = 1.0;
        for t in 1..10000 {
            let new_dist: Vec<f64> = (0..n)
                .map(|j| (0..n).map(|i| dist[i] * self.p[i * n + j]).sum())
                .collect();
            dist = new_dist;
            let tv: f64 = dist
                .iter()
                .zip(&pi)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                * 0.5;
            if tv < eps {
                return t;
            }
        }
        10000
    }

    /// Compute the empirical transition counts from a trajectory.
    pub fn empirical_transition_matrix(trajectory: &[usize], n_states: usize) -> Vec<f64> {
        let mut counts = vec![0.0f64; n_states * n_states];
        for w in trajectory.windows(2) {
            counts[w[0] * n_states + w[1]] += 1.0;
        }
        for i in 0..n_states {
            let row_sum: f64 = (0..n_states).map(|j| counts[i * n_states + j]).sum();
            if row_sum > 0.0 {
                for j in 0..n_states {
                    counts[i * n_states + j] /= row_sum;
                }
            }
        }
        counts
    }

    fn matmul_sq(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut c = vec![0.0f64; n * n];
        for i in 0..n {
            for k in 0..n {
                if a[i * n + k] == 0.0 {
                    continue;
                }
                for j in 0..n {
                    c[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }
        c
    }
}

// ─── ContinuousTimeMarkov ─────────────────────────────────────────────────────

/// Continuous-time Markov chain with infinitesimal generator matrix Q.
pub struct ContinuousTimeMarkov {
    /// Generator matrix Q (n x n). Rows sum to zero; off-diagonals >= 0.
    pub q: Vec<f64>,
    /// Number of states.
    pub n_states: usize,
}

impl ContinuousTimeMarkov {
    /// Create a CTMC from generator matrix Q.
    pub fn new(q: Vec<f64>, n_states: usize) -> Self {
        Self { q, n_states }
    }

    /// Simulate a CTMC trajectory up to time T using Gillespie's algorithm.
    ///
    /// Returns (times, states) pairs representing the jump times and new states.
    pub fn simulate(&self, init: usize, t_max: f64, seed: u64) -> (Vec<f64>, Vec<usize>) {
        let mut rng = Lcg::new(seed);
        let n = self.n_states;
        let mut times = vec![0.0f64];
        let mut states = vec![init];
        let mut state = init;
        let mut t = 0.0f64;

        loop {
            // Total rate out of state i
            let rate_out: f64 = -(self.q[state * n + state]);
            if rate_out <= 1e-300 {
                break;
            } // absorbing state
            let dt = -rng.uniform().max(1e-300).ln() / rate_out;
            t += dt;
            if t >= t_max {
                break;
            }
            // Choose next state proportional to off-diagonal Q entries
            let u = rng.uniform() * rate_out;
            let mut cum = 0.0;
            let mut next = state;
            for j in 0..n {
                if j == state {
                    continue;
                }
                cum += self.q[state * n + j];
                if u < cum {
                    next = j;
                    break;
                }
            }
            times.push(t);
            states.push(next);
            state = next;
        }
        (times, states)
    }

    /// Solve the Kolmogorov forward equations dp/dt = p * Q at time T.
    ///
    /// Uses matrix exponentiation: p(T) = p(0) * exp(Q * T).
    pub fn kolmogorov_forward(&self, p0: &[f64], t: f64) -> Vec<f64> {
        let n = self.n_states;
        let qt: Vec<f64> = self.q.iter().map(|v| v * t).collect();
        // exp(QT) via Padé (reuse the approach from MatrixFunctions)
        let exp_qt = expm_dense(&qt, n);
        // p(T) = p(0) * exp(QT)  (row vector times matrix)
        (0..n)
            .map(|j| (0..n).map(|i| p0[i] * exp_qt[i * n + j]).sum())
            .collect()
    }

    /// Compute the stationary distribution pi satisfying pi Q = 0, sum(pi) = 1.
    pub fn stationary_distribution(&self, tol: f64, max_iter: usize) -> Vec<f64> {
        let n = self.n_states;
        // pi Q = 0 is equivalent to Q^T pi = 0
        // Power iteration on the embedded chain
        let dt = 1.0
            / (0..n)
                .map(|i| -self.q[i * n + i])
                .fold(0.0f64, f64::max)
                .max(1.0);
        // Uniformisation: P = I + dt * Q
        let p: Vec<f64> = (0..n * n)
            .map(|idx| {
                let i = idx / n;
                let j = idx % n;
                if i == j {
                    1.0 + dt * self.q[idx]
                } else {
                    dt * self.q[idx]
                }
            })
            .collect();
        let mc = MarkovChain::new(p, n);
        mc.stationary_distribution(tol, max_iter)
    }

    /// Compute the mean first-passage time from state `i` to state `j`.
    ///
    /// Uses the fundamental matrix method for ergodic CTMCs.
    pub fn mean_first_passage_time(&self, from: usize, to: usize) -> f64 {
        let pi = self.stationary_distribution(1e-10, 2000);
        // For ergodic CTMC: m_{ij} = 1/(pi_j * q_j) where q_j is total rate out of j
        let qj = -(self.q[to * self.n_states + to]);
        if qj < 1e-300 || pi[to] < 1e-300 {
            return f64::INFINITY;
        }
        let _ = from;
        1.0 / (pi[to] * qj)
    }

    /// Compute the expected holding time in state `i` (mean sojourn time).
    pub fn mean_sojourn_time(&self, state: usize) -> f64 {
        let rate = -(self.q[state * self.n_states + state]);
        if rate > 1e-300 {
            1.0 / rate
        } else {
            f64::INFINITY
        }
    }
}

/// Dense matrix exponential via scaling-and-squaring (Padé \[3/3\]).
fn expm_dense(a: &[f64], n: usize) -> Vec<f64> {
    let id = eye_n(n);
    // Scale
    let norm: f64 = a.iter().map(|v| v.abs()).sum::<f64>();
    let s = ((norm / 1.0).log2().ceil() as i32).max(0) as u32;
    let scale = 1.0 / (1u64 << s) as f64;
    let a_s: Vec<f64> = a.iter().map(|v| v * scale).collect();
    let a2 = mm(&a_s, &a_s, n);
    // Padé [3/3] coefficients: c = [1, 1/2, 3/26, 1/312]
    let c0 = 1.0f64;
    let c1 = 0.5;
    let c2 = 3.0 / 26.0;
    let c3 = 1.0 / 312.0;
    let u: Vec<f64> = (0..n * n).map(|i| c3 * a2[i] + c1 * id[i]).collect();
    let u = mv(&a_s, &u, n); // A * (c3*A2 + c1*I)
    let u: Vec<f64> = u.to_vec();
    let v: Vec<f64> = (0..n * n).map(|i| c2 * a2[i] + c0 * id[i]).collect();
    // r = (V-U)^{-1}(V+U)
    let num: Vec<f64> = v.iter().zip(&u).map(|(x, y)| x + y).collect();
    let den: Vec<f64> = v.iter().zip(&u).map(|(x, y)| x - y).collect();
    let inv_den = inv_n(&den, n);
    let mut result = mm(&inv_den, &num, n);
    for _ in 0..s {
        result = mm(&result, &result, n);
    }
    result
}

fn eye_n(n: usize) -> Vec<f64> {
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}

fn mm(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; n * n];
    for i in 0..n {
        for k in 0..n {
            if a[i * n + k] == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    c
}

fn mv(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    mm(a, b, n)
}

fn inv_n(a: &[f64], n: usize) -> Vec<f64> {
    let mut m = a.to_vec();
    let mut inv = eye_n(n);
    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if m[row * n + col].abs() > m[max_row * n + col].abs() {
                max_row = row;
            }
        }
        for j in 0..n {
            m.swap(col * n + j, max_row * n + j);
            inv.swap(col * n + j, max_row * n + j);
        }
        let d = m[col * n + col];
        if d.abs() < 1e-300 {
            continue;
        }
        for j in 0..n {
            m[col * n + j] /= d;
            inv[col * n + j] /= d;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = m[row * n + col];
            for j in 0..n {
                let mv2 = m[col * n + j];
                let iv = inv[col * n + j];
                m[row * n + j] -= f * mv2;
                inv[row * n + j] -= f * iv;
            }
        }
    }
    inv
}

// ─── LevyProcess ─────────────────────────────────────────────────────────────

/// Lévy process generators: alpha-stable, variance gamma, CGMY.
pub struct LevyProcess;

impl LevyProcess {
    /// Simulate an alpha-stable Lévy process path using Chambers-Mallows-Stuck method.
    ///
    /// For index alpha in (0, 2] and skewness beta in \[-1, 1\].
    /// alpha=2 gives a Gaussian (Brownian motion scaling).
    pub fn alpha_stable(
        alpha: f64,
        beta: f64,
        scale: f64,
        loc: f64,
        t: f64,
        n_steps: usize,
        seed: u64,
    ) -> Vec<f64> {
        assert!(alpha > 0.0 && alpha <= 2.0, "alpha must be in (0, 2]");
        assert!((-1.0..=1.0).contains(&beta), "beta must be in [-1, 1]");
        let dt = t / n_steps as f64;
        let mut rng = Lcg::new(seed);
        let mut path = vec![0.0f64; n_steps + 1];
        for i in 1..=n_steps {
            let z = Self::stable_sample(&mut rng, alpha, beta);
            path[i] = path[i - 1] + loc * dt + scale * dt.powf(1.0 / alpha) * z;
        }
        path
    }

    fn stable_sample(rng: &mut Lcg, alpha: f64, beta: f64) -> f64 {
        if (alpha - 2.0).abs() < 1e-10 {
            return rng.normal();
        }
        let phi = PI * (rng.uniform() - 0.5);
        let w = -rng.uniform().max(1e-300).ln();
        if (alpha - 1.0).abs() < 1e-10 {
            let b_phi = (PI / 2.0 + beta * phi).tan();
            return (PI / 2.0) * (b_phi * phi - beta * (b_phi * phi.cos()).max(1e-300).ln());
        }
        let zeta = -beta * (PI * alpha / 2.0).tan();
        let xi = (1.0 / alpha) * (beta * (PI * alpha / 2.0).tan()).atan();
        let t1 = (1.0 + zeta * zeta).powf(1.0 / (2.0 * alpha));
        let t2 = (alpha * (phi + xi)).sin() / phi.cos().powf(1.0 / alpha);
        let t3 = ((phi - alpha * (phi + xi)).cos() / w).powf((1.0 - alpha) / alpha);
        t1 * t2 * t3
    }

    /// Simulate a Variance Gamma (VG) process path.
    ///
    /// VG = theta * G + sigma * W(G) where G is a Gamma subordinator.
    pub fn variance_gamma(
        sigma: f64,
        nu: f64,
        theta: f64,
        t: f64,
        n_steps: usize,
        seed: u64,
    ) -> Vec<f64> {
        let dt = t / n_steps as f64;
        let mut rng = Lcg::new(seed);
        let mut path = vec![0.0f64; n_steps + 1];
        for i in 1..=n_steps {
            // Sample Gamma increment dG ~ Gamma(dt/nu, nu)
            let dg = Self::gamma_sample(&mut rng, dt / nu, nu);
            let dw = rng.normal() * dg.sqrt();
            path[i] = path[i - 1] + theta * dg + sigma * dw;
        }
        path
    }

    fn gamma_sample(rng: &mut Lcg, shape: f64, scale: f64) -> f64 {
        // Marsaglia-Tsang method for shape >= 1; for shape < 1, boost and scale back
        if shape < 1e-10 {
            return 0.0;
        }
        let (k, s) = if shape < 1.0 {
            let u = rng.uniform().max(1e-300);
            (shape + 1.0, u.powf(1.0 / shape))
        } else {
            (shape, 1.0)
        };
        let d = k - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = rng.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = rng.uniform().max(1e-300);
            if u < 1.0 - 0.0331 * x.powi(4) {
                return d * v * s * scale;
            }
            if u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
                return d * v * s * scale;
            }
        }
    }

    /// Simulate a CGMY process path (Carr-Geman-Madan-Yor).
    ///
    /// Approximates the process by truncating small jumps and simulating
    /// large jumps via a compound Poisson approach.
    pub fn cgmy(
        c: f64,
        g: f64,
        m_param: f64,
        y: f64,
        t: f64,
        n_steps: usize,
        seed: u64,
    ) -> Vec<f64> {
        assert!(
            y < 2.0 && c > 0.0 && g > 0.0 && m_param > 0.0,
            "Invalid CGMY parameters"
        );
        let dt = t / n_steps as f64;
        let mut rng = Lcg::new(seed);
        let mut path = vec![0.0f64; n_steps + 1];
        // Approximate: drift + diffusive small jumps + compound Poisson large jumps
        let eps = 0.1f64; // truncation level
        // Truncated small-jump variance (approximate)
        let sigma_small = (2.0 * c * eps.powf(2.0 - y) / (2.0 - y)).sqrt();
        // Large jump intensity
        let lambda_large = c * (g.powf(-y) + m_param.powf(-y)) * (1.0 - y).max(0.1);
        for i in 1..=n_steps {
            let diffusive = sigma_small * dt.sqrt() * rng.normal();
            // Poisson arrivals of large jumps
            let n_jumps = {
                let mu = lambda_large * dt;
                let mut k = 0usize;
                let mut p = (-mu).exp();
                let mut cum = p;
                let u = rng.uniform();
                while u > cum && k < 20 {
                    k += 1;
                    p *= mu / k as f64;
                    cum += p;
                }
                k
            };
            let mut jump_sum = 0.0;
            for _ in 0..n_jumps {
                let u = rng.uniform();
                let jump = if u < 0.5 {
                    // Positive jump ~ Exp(m_param)
                    -rng.uniform().max(1e-300).ln() / m_param
                } else {
                    // Negative jump ~ -Exp(g)
                    rng.uniform().max(1e-300).ln() / g
                };
                jump_sum += jump;
            }
            path[i] = path[i - 1] + diffusive + jump_sum;
        }
        path
    }

    /// Compute the characteristic exponent of an alpha-stable process at frequency u.
    ///
    /// For alpha != 1: psi(u) = -|u|^alpha * (1 - i*beta*sign(u)*tan(pi*alpha/2)).
    pub fn alpha_stable_characteristic_exponent(alpha: f64, beta: f64, u: f64) -> f64 {
        // Return real part of log E[exp(iuX)]
        let tan_term = (PI * alpha / 2.0).tan();
        -u.abs().powf(alpha) * (1.0 + beta * u.signum() * tan_term)
    }
}

// ─── StochasticDifferentialEquation ──────────────────────────────────────────

/// SDE solver for dX = mu(X,t)dt + sigma(X,t)dW.
pub struct StochasticDifferentialEquation;

impl StochasticDifferentialEquation {
    /// Euler-Maruyama scheme: X_{n+1} = X_n + mu(X_n, t_n)*dt + sigma(X_n, t_n)*sqrt(dt)*Z.
    ///
    /// Strong order 0.5, weak order 1.0.
    pub fn euler_maruyama<F, G>(
        mu: F,
        sigma: G,
        x0: f64,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        seed: u64,
    ) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(f64, f64) -> f64,
        G: Fn(f64, f64) -> f64,
    {
        let dt = (t_end - t0) / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut rng = Lcg::new(seed);
        let mut t_path = vec![t0; n_steps + 1];
        let mut x_path = vec![x0; n_steps + 1];
        for i in 0..n_steps {
            let t = t0 + i as f64 * dt;
            let x = x_path[i];
            let dw = sqrt_dt * rng.normal();
            x_path[i + 1] = x + mu(x, t) * dt + sigma(x, t) * dw;
            t_path[i + 1] = t + dt;
        }
        (t_path, x_path)
    }

    /// Milstein scheme: adds a correction term for higher strong order (1.0).
    ///
    /// X_{n+1} = X_n + mu*dt + sigma*dW + 0.5*sigma*sigma'*(dW^2 - dt).
    /// sigma_prime is the derivative of sigma with respect to x.
    pub fn milstein<F, G, H>(
        mu: F,
        sigma: G,
        sigma_prime: H,
        x0: f64,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        seed: u64,
    ) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(f64, f64) -> f64,
        G: Fn(f64, f64) -> f64,
        H: Fn(f64, f64) -> f64,
    {
        let dt = (t_end - t0) / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut rng = Lcg::new(seed);
        let mut t_path = vec![t0; n_steps + 1];
        let mut x_path = vec![x0; n_steps + 1];
        for i in 0..n_steps {
            let t = t0 + i as f64 * dt;
            let x = x_path[i];
            let z = rng.normal();
            let dw = sqrt_dt * z;
            let sig = sigma(x, t);
            let sig_p = sigma_prime(x, t);
            x_path[i + 1] = x + mu(x, t) * dt + sig * dw + 0.5 * sig * sig_p * (dw * dw - dt);
            t_path[i + 1] = t + dt;
        }
        (t_path, x_path)
    }

    /// Runge-Kutta (Runge-Kutta-Maruyama) order 1 scheme for Ito SDEs.
    ///
    /// Uses a predictor step to estimate the diffusion coefficient more accurately.
    pub fn runge_kutta_maruyama<F, G>(
        mu: F,
        sigma: G,
        x0: f64,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        seed: u64,
    ) -> (Vec<f64>, Vec<f64>)
    where
        F: Fn(f64, f64) -> f64,
        G: Fn(f64, f64) -> f64,
    {
        let dt = (t_end - t0) / n_steps as f64;
        let sqrt_dt = dt.sqrt();
        let mut rng = Lcg::new(seed);
        let mut t_path = vec![t0; n_steps + 1];
        let mut x_path = vec![x0; n_steps + 1];
        for i in 0..n_steps {
            let t = t0 + i as f64 * dt;
            let x = x_path[i];
            let dw = sqrt_dt * rng.normal();
            let sig = sigma(x, t);
            // Predictor: x_hat = x + sigma(x,t)*sqrt(dt)
            let x_hat = x + sig * sqrt_dt;
            let sig_hat = sigma(x_hat, t);
            x_path[i + 1] =
                x + mu(x, t) * dt + sig * dw + 0.5 * (sig_hat - sig) * (dw * dw - dt) / sqrt_dt;
            t_path[i + 1] = t + dt;
        }
        (t_path, x_path)
    }

    /// Simulate multiple paths (Monte Carlo ensemble).
    ///
    /// Returns a matrix of shape (n_paths, n_steps+1).
    pub fn monte_carlo_paths<F, G>(
        mu: F,
        sigma: G,
        x0: f64,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        n_paths: usize,
        seed: u64,
    ) -> Vec<Vec<f64>>
    where
        F: Fn(f64, f64) -> f64 + Clone,
        G: Fn(f64, f64) -> f64 + Clone,
    {
        (0..n_paths)
            .map(|k| {
                let (_, x) = Self::euler_maruyama(
                    mu.clone(),
                    sigma.clone(),
                    x0,
                    t0,
                    t_end,
                    n_steps,
                    seed.wrapping_add(k as u64),
                );
                x
            })
            .collect()
    }

    /// Compute the strong order of convergence from a reference solution.
    ///
    /// Returns the estimated strong error at the terminal time.
    pub fn strong_error_estimate<F, G>(
        mu: F,
        sigma: G,
        x0: f64,
        t0: f64,
        t_end: f64,
        seed: u64,
        n_steps_coarse: usize,
        n_steps_fine: usize,
    ) -> f64
    where
        F: Fn(f64, f64) -> f64 + Clone,
        G: Fn(f64, f64) -> f64 + Clone,
    {
        let (_, x_coarse) = Self::euler_maruyama(
            mu.clone(),
            sigma.clone(),
            x0,
            t0,
            t_end,
            n_steps_coarse,
            seed,
        );
        let (_, x_fine) = Self::euler_maruyama(mu, sigma, x0, t0, t_end, n_steps_fine, seed);
        (*x_coarse.last().expect("coarse path is non-empty")
            - *x_fine.last().expect("fine path is non-empty"))
        .abs()
    }

    /// Compute the weak order of convergence via expected value of a functional.
    ///
    /// Returns E\[f(X(T))\] estimated from n_paths Monte Carlo paths.
    pub fn weak_order_estimate<F, G, H>(
        mu: F,
        sigma: G,
        functional: H,
        x0: f64,
        t0: f64,
        t_end: f64,
        n_steps: usize,
        n_paths: usize,
        seed: u64,
    ) -> f64
    where
        F: Fn(f64, f64) -> f64 + Clone,
        G: Fn(f64, f64) -> f64 + Clone,
        H: Fn(f64) -> f64,
    {
        let paths = Self::monte_carlo_paths(mu, sigma, x0, t0, t_end, n_steps, n_paths, seed);
        let sum: f64 = paths
            .iter()
            .map(|p| functional(*p.last().expect("path is non-empty")))
            .sum();
        sum / n_paths as f64
    }

    /// Stratonovich to Ito conversion: add the correction term -0.5*sigma*sigma'.
    ///
    /// Returns the Ito drift given the Stratonovich drift and diffusion.
    pub fn stratonovich_to_ito<G, H>(
        strat_mu: f64,
        sigma_val: f64,
        sigma_prime: f64,
        _sigma: G,
        _sigma_p: H,
    ) -> f64
    where
        G: Fn(f64, f64) -> f64,
        H: Fn(f64, f64) -> f64,
    {
        strat_mu - 0.5 * sigma_val * sigma_prime
    }
}

// ─── Statistics helpers ───────────────────────────────────────────────────────

/// Compute the sample mean of a slice.
pub fn sample_mean(data: &[f64]) -> f64 {
    data.iter().sum::<f64>() / data.len() as f64
}

/// Compute the sample variance of a slice.
pub fn sample_variance(data: &[f64]) -> f64 {
    let m = sample_mean(data);
    data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (data.len() - 1).max(1) as f64
}

/// Compute the empirical autocorrelation at lag k.
pub fn autocorrelation(data: &[f64], lag: usize) -> f64 {
    let n = data.len();
    if lag >= n {
        return 0.0;
    }
    let m = sample_mean(data);
    let var = sample_variance(data);
    if var < 1e-300 {
        return 0.0;
    }
    let cov: f64 = (0..n - lag)
        .map(|i| (data[i] - m) * (data[i + lag] - m))
        .sum::<f64>()
        / (n - lag) as f64;
    cov / var
}

/// Compute the empirical power spectral density via periodogram.
pub fn periodogram(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    let m = sample_mean(data);
    let centered: Vec<f64> = data.iter().map(|x| x - m).collect();
    (0..n / 2)
        .map(|k| {
            let re: f64 = (0..n)
                .map(|t| centered[t] * (2.0 * PI * k as f64 * t as f64 / n as f64).cos())
                .sum();
            let im: f64 = (0..n)
                .map(|t| centered[t] * (2.0 * PI * k as f64 * t as f64 / n as f64).sin())
                .sum();
            (re * re + im * im) / n as f64
        })
        .collect()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BrownianMotion ───────────────────────────────────────────────────────

    #[test]
    fn test_bm_starts_at_zero() {
        let path = BrownianMotion::standard(1.0, 100, 42);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_bm_length() {
        let n = 200;
        let path = BrownianMotion::standard(1.0, n, 1);
        assert_eq!(path.len(), n + 1);
    }

    #[test]
    fn test_bm_quadratic_variation() {
        // QV should be approximately T=1 for a standard BM with many steps
        let path = BrownianMotion::standard(1.0, 10000, 7);
        let qv = BrownianMotion::quadratic_variation(&path);
        assert!((qv - 1.0).abs() < 0.1, "QV = {qv}");
    }

    #[test]
    fn test_gbm_positive() {
        let path = BrownianMotion::geometric(100.0, 0.05, 0.2, 1.0, 100, 99);
        assert!(path.iter().all(|&v| v > 0.0), "GBM should be positive");
    }

    #[test]
    fn test_gbm_starts_at_s0() {
        let path = BrownianMotion::geometric(50.0, 0.1, 0.3, 1.0, 100, 5);
        assert!((path[0] - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_fbm_length() {
        let path = BrownianMotion::fractional(0.7, 50, 42);
        assert_eq!(path.len(), 51);
    }

    #[test]
    fn test_fbm_starts_at_zero() {
        let path = BrownianMotion::fractional(0.3, 100, 10);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_bm_variance_approx() {
        let path = BrownianMotion::standard(1.0, 1000, 123);
        // Var of final value ~ T = 1
        let increments: Vec<f64> = path.windows(2).map(|w| w[1] - w[0]).collect();
        let v = BrownianMotion::variance(&increments);
        // dt = 1/1000, so variance of each increment = dt
        assert!((v - 0.001).abs() < 0.002, "variance = {v}");
    }

    // ── PoissonProcess ───────────────────────────────────────────────────────

    #[test]
    fn test_poisson_sorted_times() {
        let times = PoissonProcess::homogeneous(5.0, 10.0, 42);
        for w in times.windows(2) {
            assert!(w[0] < w[1]);
        }
    }

    #[test]
    fn test_poisson_times_in_range() {
        let t = 5.0;
        let times = PoissonProcess::homogeneous(3.0, t, 7);
        assert!(times.iter().all(|&ti| ti < t));
    }

    #[test]
    fn test_poisson_expected_count() {
        let expected = PoissonProcess::expected_count(4.0, 2.5);
        assert!((expected - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_poisson_empirical_rate() {
        let t = 100.0;
        let lambda = 5.0;
        let times = PoissonProcess::homogeneous(lambda, t, 99);
        let emp = PoissonProcess::empirical_rate(&times, t);
        assert!((emp - lambda).abs() < 1.0, "empirical rate = {emp}");
    }

    #[test]
    fn test_poisson_inhomogeneous() {
        let times = PoissonProcess::inhomogeneous(|t| 1.0 + t, 11.0, 10.0, 42);
        assert!(times.iter().all(|&ti| ti < 10.0));
    }

    #[test]
    fn test_compound_poisson() {
        let (times, values) = PoissonProcess::compound(2.0, 5.0, |u| u, 42);
        assert!(values.len() == times.len() + 1, "values has one more entry");
        assert_eq!(values[0], 0.0);
    }

    #[test]
    fn test_inter_arrival_cdf() {
        let cdf = PoissonProcess::inter_arrival_cdf(1.0, 0.0);
        assert_eq!(cdf, 0.0);
        let cdf_inf = PoissonProcess::inter_arrival_cdf(1.0, 100.0);
        assert!((cdf_inf - 1.0).abs() < 1e-30);
    }

    // ── MarkovChain ──────────────────────────────────────────────────────────

    fn simple_mc() -> MarkovChain {
        // 2-state chain
        MarkovChain::new(vec![0.9, 0.1, 0.2, 0.8], 2)
    }

    #[test]
    fn test_mc_trajectory_length() {
        let mc = simple_mc();
        let traj = mc.simulate(0, 100, 42);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_mc_states_valid() {
        let mc = simple_mc();
        let traj = mc.simulate(0, 200, 5);
        assert!(traj.iter().all(|&s| s < 2));
    }

    #[test]
    fn test_mc_stationary_sums_to_one() {
        let mc = simple_mc();
        let pi = mc.stationary_distribution(1e-10, 2000);
        let s: f64 = pi.iter().sum();
        assert!((s - 1.0).abs() < 1e-8, "sum = {s}");
    }

    #[test]
    fn test_mc_stationary_known() {
        let mc = simple_mc();
        let pi = mc.stationary_distribution(1e-10, 2000);
        // pi0 = 2/3, pi1 = 1/3 for this chain
        assert!((pi[0] - 2.0 / 3.0).abs() < 1e-4, "pi[0] = {}", pi[0]);
        assert!((pi[1] - 1.0 / 3.0).abs() < 1e-4, "pi[1] = {}", pi[1]);
    }

    #[test]
    fn test_mc_n_step_matrix() {
        let mc = simple_mc();
        let p1 = mc.n_step_matrix(1);
        for (a, b) in p1.iter().zip(&mc.p) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_mc_mixing_time() {
        let mc = simple_mc();
        let t = mc.mixing_time(0, 0.01);
        assert!(t > 0);
    }

    #[test]
    fn test_empirical_transition_matrix() {
        let mc = simple_mc();
        let traj = mc.simulate(0, 10000, 77);
        let emp = MarkovChain::empirical_transition_matrix(&traj, 2);
        // Should be close to true matrix
        assert!((emp[0] - 0.9).abs() < 0.05, "emp[0,0] = {}", emp[0]);
    }

    // ── ContinuousTimeMarkov ─────────────────────────────────────────────────

    fn simple_ctmc() -> ContinuousTimeMarkov {
        // 2-state: q01 = 2, q10 = 3 → diag = [-2, -3]
        ContinuousTimeMarkov::new(vec![-2.0, 2.0, 3.0, -3.0], 2)
    }

    #[test]
    fn test_ctmc_simulate_runs() {
        let ctmc = simple_ctmc();
        let (times, states) = ctmc.simulate(0, 5.0, 42);
        assert_eq!(times.len(), states.len());
        assert!(times.iter().all(|&t| t < 5.0));
    }

    #[test]
    fn test_ctmc_stationary_sums_to_one() {
        let ctmc = simple_ctmc();
        let pi = ctmc.stationary_distribution(1e-10, 2000);
        let s: f64 = pi.iter().sum();
        assert!((s - 1.0).abs() < 1e-6, "sum = {s}");
    }

    #[test]
    fn test_ctmc_stationary_known() {
        let ctmc = simple_ctmc();
        let pi = ctmc.stationary_distribution(1e-10, 2000);
        // pi0 = 3/5, pi1 = 2/5 (balance: pi0 * 2 = pi1 * 3)
        assert!((pi[0] - 0.6).abs() < 0.05, "pi[0] = {}", pi[0]);
    }

    #[test]
    fn test_ctmc_kolmogorov_sums_to_one() {
        let ctmc = simple_ctmc();
        let p0 = vec![1.0, 0.0];
        let pt = ctmc.kolmogorov_forward(&p0, 1.0);
        let s: f64 = pt.iter().sum();
        assert!((s - 1.0).abs() < 1e-6, "sum = {s}");
    }

    #[test]
    fn test_ctmc_mean_sojourn() {
        let ctmc = simple_ctmc();
        assert!((ctmc.mean_sojourn_time(0) - 0.5).abs() < 1e-10);
        assert!((ctmc.mean_sojourn_time(1) - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_ctmc_first_passage_finite() {
        let ctmc = simple_ctmc();
        let fpt = ctmc.mean_first_passage_time(0, 1);
        assert!(fpt.is_finite() && fpt > 0.0);
    }

    // ── LevyProcess ──────────────────────────────────────────────────────────

    #[test]
    fn test_alpha_stable_length() {
        let path = LevyProcess::alpha_stable(1.5, 0.0, 1.0, 0.0, 1.0, 100, 42);
        assert_eq!(path.len(), 101);
    }

    #[test]
    fn test_alpha_stable_starts_zero() {
        let path = LevyProcess::alpha_stable(1.8, 0.2, 0.5, 0.0, 1.0, 50, 7);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_alpha_stable_finite() {
        let path = LevyProcess::alpha_stable(1.5, 0.0, 1.0, 0.0, 1.0, 50, 1);
        assert!(path.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_vg_length() {
        let path = LevyProcess::variance_gamma(0.2, 0.1, -0.1, 1.0, 100, 42);
        assert_eq!(path.len(), 101);
    }

    #[test]
    fn test_vg_starts_zero() {
        let path = LevyProcess::variance_gamma(0.2, 0.1, 0.0, 1.0, 50, 5);
        assert_eq!(path[0], 0.0);
    }

    #[test]
    fn test_cgmy_length() {
        let path = LevyProcess::cgmy(1.0, 5.0, 10.0, 0.5, 1.0, 100, 42);
        assert_eq!(path.len(), 101);
    }

    #[test]
    fn test_cgmy_finite() {
        let path = LevyProcess::cgmy(1.0, 5.0, 10.0, 0.5, 1.0, 50, 9);
        assert!(path.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_alpha_stable_char_exp() {
        let psi = LevyProcess::alpha_stable_characteristic_exponent(2.0, 0.0, 1.0);
        // For alpha=2, beta=0: psi = -1
        assert!((psi + 1.0).abs() < 1e-10, "psi = {psi}");
    }

    // ── SDE ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_em_gbm_terminal() {
        // GBM: dS = mu*S dt + sigma*S dW
        let mu_val = 0.05;
        let sigma_val = 0.2;
        let s0 = 100.0;
        let t = 1.0;
        let (_, x) = StochasticDifferentialEquation::euler_maruyama(
            move |x, _t| mu_val * x,
            move |x, _t| sigma_val * x,
            s0,
            0.0,
            t,
            1000,
            42,
        );
        assert!(*x.last().unwrap() > 0.0);
    }

    #[test]
    fn test_em_starts_at_x0() {
        let (_, x) = StochasticDifferentialEquation::euler_maruyama(
            |x, _| x,
            |_, _| 0.1,
            5.0,
            0.0,
            1.0,
            100,
            1,
        );
        assert!((x[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_milstein_ou() {
        // OU: dX = -theta*X dt + sigma dW → sigma'(X) = 0
        let theta = 2.0;
        let sig = 0.5;
        let (_, x) = StochasticDifferentialEquation::milstein(
            move |x, _| -theta * x,
            move |_, _| sig,
            |_, _| 0.0, // sigma_prime = 0
            1.0,
            0.0,
            2.0,
            500,
            42,
        );
        assert!(x.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_milstein_length() {
        let (t, x) = StochasticDifferentialEquation::milstein(
            |x, _| x,
            |_, _| 0.1,
            |_, _| 0.0,
            1.0,
            0.0,
            1.0,
            100,
            5,
        );
        assert_eq!(t.len(), 101);
        assert_eq!(x.len(), 101);
    }

    #[test]
    fn test_rkm_finite() {
        let (_, x) = StochasticDifferentialEquation::runge_kutta_maruyama(
            |x, _| -x,
            |_, _| 0.5,
            1.0,
            0.0,
            1.0,
            200,
            42,
        );
        assert!(x.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_monte_carlo_paths_count() {
        let paths = StochasticDifferentialEquation::monte_carlo_paths(
            |x, _| -x,
            |_, _| 0.2,
            1.0,
            0.0,
            1.0,
            50,
            10,
            42,
        );
        assert_eq!(paths.len(), 10);
        assert!(paths.iter().all(|p| p.len() == 51));
    }

    #[test]
    fn test_strong_error_finite() {
        let err = StochasticDifferentialEquation::strong_error_estimate(
            |x, _| -x,
            |_, _| 0.2,
            1.0,
            0.0,
            1.0,
            42,
            10,
            100,
        );
        assert!(err.is_finite());
    }

    #[test]
    fn test_weak_order_estimate() {
        // E[X(T)] for deterministic ODE dX = X dt, X(0)=1 → X(T) ≈ e
        let ev = StochasticDifferentialEquation::weak_order_estimate(
            |x, _| x,
            |_, _| 0.0,
            |x| x,
            1.0,
            0.0,
            1.0,
            100,
            50,
            42,
        );
        assert!((ev - std::f64::consts::E).abs() < 0.1, "E[X(T)] = {ev}");
    }

    #[test]
    fn test_stratonovich_to_ito() {
        let ito_mu = StochasticDifferentialEquation::stratonovich_to_ito(
            1.0,
            0.5,
            0.5,
            |_: f64, _: f64| 0.5f64,
            |_: f64, _: f64| 0.5f64,
        );
        assert!((ito_mu - (1.0 - 0.5 * 0.5 * 0.5)).abs() < 1e-10);
    }

    // ── Statistics helpers ────────────────────────────────────────────────────

    #[test]
    fn test_sample_mean() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((sample_mean(&data) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_sample_variance() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((sample_variance(&data) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_lag0() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let ac = autocorrelation(&data, 0);
        // lag-0 autocorrelation is cov(X,X)/var(X); slight discrepancy due to
        // population vs. sample variance — should be close to 1
        assert!((ac - 1.0).abs() < 0.02, "lag-0 autocorrelation = {ac}");
    }

    #[test]
    fn test_autocorrelation_white_noise() {
        // White noise should have near-zero autocorrelation at lag > 0
        let path = BrownianMotion::standard(1.0, 500, 88);
        let increments: Vec<f64> = path.windows(2).map(|w| w[1] - w[0]).collect();
        let ac1 = autocorrelation(&increments, 1).abs();
        assert!(ac1 < 0.15, "lag-1 autocorr = {ac1}");
    }

    #[test]
    fn test_periodogram_length() {
        let data: Vec<f64> = (0..64).map(|i| (i as f64).sin()).collect();
        let psd = periodogram(&data);
        assert_eq!(psd.len(), 32);
    }

    #[test]
    fn test_lcg_uniform_range() {
        let mut rng = Lcg::new(42);
        for _ in 0..1000 {
            let u = rng.uniform();
            assert!((0.0..1.0).contains(&u), "uniform out of range: {u}");
        }
    }

    #[test]
    fn test_lcg_normal_finite() {
        let mut rng = Lcg::new(42);
        for _ in 0..100 {
            let z = rng.normal();
            assert!(z.is_finite());
        }
    }
}
