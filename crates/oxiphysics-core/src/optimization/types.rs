//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{cholesky_lower, forward_solve_lower, simulated_annealing};
use rand::RngExt;

/// Facade that bundles several optimization algorithms under a single struct.
///
/// Each method delegates to the free-standing implementations in this module.
pub struct Optimizer;
impl Optimizer {
    /// Create a new `Optimizer` facade.
    pub fn new() -> Self {
        Self
    }
    /// Global minimisation via Differential Evolution (DE/rand/1/bin).
    ///
    /// Initialises a population in `bounds`, runs `max_gens` generations,
    /// and returns the best solution found.
    ///
    /// # Arguments
    /// * `f`        - Objective function (minimised).
    /// * `bounds`   - Per-dimension `(lo, hi)` constraints.
    /// * `pop_size` - Number of individuals (≥ 4).
    /// * `f_weight` - Differential weight F ∈ (0, 2].
    /// * `cr`       - Crossover probability ∈ \[0, 1\].
    /// * `max_gens` - Number of generations.
    pub fn differential_evolution<F>(
        &self,
        f: F,
        bounds: &[(f64, f64)],
        pop_size: usize,
        f_weight: f64,
        cr: f64,
        max_gens: u32,
    ) -> OptResult
    where
        F: Fn(&[f64]) -> f64,
    {
        let dim = bounds.len();
        let mut rng = rand::rng();
        let mut de = DifferentialEvolution::new(pop_size.max(4), dim, bounds, &mut rng);
        de.f_weight = f_weight;
        de.cr = cr;
        for _ in 0..max_gens {
            de.step(&f, &mut rng);
        }
        let best = de.best(&f);
        let f_val = f(&best);
        OptResult {
            x: best,
            f_val,
            n_iter: max_gens,
            converged: true,
        }
    }
    /// Global minimisation via Simulated Annealing.
    ///
    /// Delegates to [`simulated_annealing`].
    pub fn simulated_annealing<F>(
        &self,
        f: F,
        x0: Vec<f64>,
        initial_temp: f64,
        cooling_rate: f64,
        step_size: f64,
        max_iter: u32,
    ) -> OptResult
    where
        F: Fn(&[f64]) -> f64,
    {
        simulated_annealing(f, x0, initial_temp, cooling_rate, step_size, max_iter)
    }
    /// One CMA-ES step: update mean, covariance, and step size.
    ///
    /// Delegates to [`CmaEs::step`].
    pub fn cmaes_step<F>(&self, cmaes: &mut CmaEs, f: &F) -> f64
    where
        F: Fn(&[f64]) -> f64,
    {
        cmaes.step(f)
    }
}
/// A single particle in PSO.
#[derive(Clone)]
pub(crate) struct Particle {
    pub(crate) pos: Vec<f64>,
    pub(crate) vel: Vec<f64>,
    pub(crate) best_pos: Vec<f64>,
    pub(crate) best_val: f64,
}
/// Differential Evolution optimizer.
///
/// Uses DE/rand/1/bin strategy: for each individual, a trial vector is formed
/// from a random base vector plus `scale_f * (r1 - r2)` and crossed with the target.
pub struct DifferentialEvolution {
    /// Population of candidate vectors.
    pub pop: Vec<Vec<f64>>,
    /// Differential weight (`F` in the DE literature) ∈ (0, 2].
    pub f_weight: f64,
    /// Crossover probability (`CR` in the DE literature) ∈ \[0, 1\].
    pub cr: f64,
}
impl DifferentialEvolution {
    /// Initialise a population of `n` individuals in `dim` dimensions,
    /// each component drawn uniformly from `bounds[i] = (lo, hi)`.
    pub fn new(n: usize, dim: usize, bounds: &[(f64, f64)], rng: &mut impl rand::Rng) -> Self {
        let pop: Vec<Vec<f64>> = (0..n)
            .map(|_| {
                (0..dim)
                    .map(|d| {
                        let (lo, hi) = bounds[d];
                        lo + rng.random::<f64>() * (hi - lo)
                    })
                    .collect()
            })
            .collect();
        Self {
            pop,
            f_weight: 0.8,
            cr: 0.9,
        }
    }
    /// Perform one generation of DE.
    ///
    /// Each individual is challenged by a trial vector formed by mutation and
    /// binomial crossover.
    pub fn step<E: Fn(&[f64]) -> f64>(&mut self, fitness: E, rng: &mut impl rand::Rng) {
        let n = self.pop.len();
        if n < 4 {
            return;
        }
        let dim = self.pop[0].len();
        let mut new_pop = self.pop.clone();
        for (i, new_entry) in new_pop.iter_mut().enumerate() {
            let mut indices = Vec::with_capacity(3);
            while indices.len() < 3 {
                let r = (rng.random::<f64>() * n as f64) as usize % n;
                if r != i && !indices.contains(&r) {
                    indices.push(r);
                }
            }
            let (a, b, c) = (indices[0], indices[1], indices[2]);
            let mutant: Vec<f64> = (0..dim)
                .map(|d| self.pop[a][d] + self.f_weight * (self.pop[b][d] - self.pop[c][d]))
                .collect();
            let j_rand = (rng.random::<f64>() * dim as f64) as usize % dim;
            let trial: Vec<f64> = (0..dim)
                .map(|d| {
                    if d == j_rand || rng.random::<f64>() < self.cr {
                        mutant[d]
                    } else {
                        self.pop[i][d]
                    }
                })
                .collect();
            if fitness(&trial) <= fitness(&self.pop[i]) {
                *new_entry = trial;
            }
        }
        self.pop = new_pop;
    }
    /// Return the best individual in the current population.
    pub fn best<E: Fn(&[f64]) -> f64>(&self, fitness: E) -> Vec<f64> {
        self.pop
            .iter()
            .min_by(|a, b| {
                fitness(a)
                    .partial_cmp(&fitness(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .cloned()
            .unwrap_or_default()
    }
}
/// Nelder-Mead simplex optimizer (struct-based, complement to the free function).
///
/// Maintains the simplex state so that you can inspect or mutate it between
/// calls.
pub struct NelderMead {
    /// Simplex vertices (n+1 vectors of length n).
    pub simplex: Vec<Vec<f64>>,
    /// Dimension of the search space.
    pub n: usize,
}
impl NelderMead {
    /// Initialise a regular simplex around `x0`.
    ///
    /// Each vertex is `x0` with one coordinate shifted by `step`.
    pub fn new(x0: Vec<f64>, step: f64) -> Self {
        let n = x0.len();
        let mut simplex: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
        simplex.push(x0.clone());
        for i in 0..n {
            let mut v = x0.clone();
            v[i] += step;
            simplex.push(v);
        }
        Self { simplex, n }
    }
    /// Run the Nelder-Mead algorithm until convergence or `max_iter` is reached.
    ///
    /// Returns the best vertex found.
    pub fn minimize<F>(&mut self, f: F, max_iter: usize, tol: f64) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = self.n;
        let mut fvals: Vec<f64> = self.simplex.iter().map(|v| f(v)).collect();
        let alpha = 1.0_f64;
        let gamma = 2.0_f64;
        let rho = 0.5_f64;
        let sigma = 0.5_f64;
        for _ in 0..max_iter {
            let mut order: Vec<usize> = (0..=n).collect();
            order.sort_by(|&a, &b| {
                fvals[a]
                    .partial_cmp(&fvals[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            self.simplex = order.iter().map(|&i| self.simplex[i].clone()).collect();
            fvals = order.iter().map(|&i| fvals[i]).collect();
            let fmean = fvals.iter().sum::<f64>() / (n + 1) as f64;
            let fstd: f64 =
                (fvals.iter().map(|fi| (fi - fmean).powi(2)).sum::<f64>() / (n + 1) as f64).sqrt();
            if fstd < tol {
                break;
            }
            let mut centroid = vec![0.0f64; n];
            for v in &self.simplex[..n] {
                for (c, vi) in centroid.iter_mut().zip(v.iter()) {
                    *c += vi;
                }
            }
            for c in centroid.iter_mut() {
                *c /= n as f64;
            }
            let xr: Vec<f64> = centroid
                .iter()
                .zip(self.simplex[n].iter())
                .map(|(c, w)| c + alpha * (c - w))
                .collect();
            let fr = f(&xr);
            if fr < fvals[0] {
                let xe: Vec<f64> = centroid
                    .iter()
                    .zip(xr.iter())
                    .map(|(c, r)| c + gamma * (r - c))
                    .collect();
                let fe = f(&xe);
                if fe < fr {
                    self.simplex[n] = xe;
                    fvals[n] = fe;
                } else {
                    self.simplex[n] = xr;
                    fvals[n] = fr;
                }
            } else if fr < fvals[n - 1] {
                self.simplex[n] = xr;
                fvals[n] = fr;
            } else {
                let xc: Vec<f64> = centroid
                    .iter()
                    .zip(self.simplex[n].iter())
                    .map(|(c, w)| c + rho * (w - c))
                    .collect();
                let fc = f(&xc);
                if fc < fvals[n] {
                    self.simplex[n] = xc;
                    fvals[n] = fc;
                } else {
                    let x0b = self.simplex[0].clone();
                    for (simplex_entry, fval) in
                        self.simplex[1..=n].iter_mut().zip(fvals[1..=n].iter_mut())
                    {
                        *simplex_entry = x0b
                            .iter()
                            .zip(simplex_entry.iter())
                            .map(|(b, v)| b + sigma * (v - b))
                            .collect();
                        *fval = f(simplex_entry);
                    }
                }
            }
        }
        self.simplex[0].clone()
    }
}
/// CMA-ES: Covariance Matrix Adaptation Evolution Strategy.
///
/// A state-of-the-art derivative-free optimizer for moderately-sized problems
/// (typically `n ≤ 100` dimensions).  Maintains a multivariate Gaussian
/// distribution over the search space, adapting the mean, covariance matrix,
/// and step size at each generation.
///
/// ## Reference
/// Hansen, N. (2016). *The CMA Evolution Strategy: A Tutorial*.
/// arXiv:1604.00772.
pub struct CmaEs {
    /// Current mean (best estimate of the optimum).
    pub mean: Vec<f64>,
    /// Step size σ.
    pub sigma: f64,
    /// Population size λ.
    pub lam: usize,
    /// Number of parents μ = λ / 2.
    pub mu: usize,
    /// Recombination weights (length μ, sum to 1).
    pub weights: Vec<f64>,
    /// Covariance matrix (flattened row-major, n×n).
    pub cov: Vec<f64>,
    /// Evolution path p_c (for rank-one update).
    pub p_c: Vec<f64>,
    /// Evolution path p_σ (for step-size control).
    pub p_sigma: Vec<f64>,
    /// Cumulation factor c_c.
    pub c_c: f64,
    /// Step-size cumulation c_σ.
    pub c_sigma: f64,
    /// Learning rate for rank-one update c_1.
    pub c1: f64,
    /// Learning rate for rank-μ update c_mu.
    pub c_mu: f64,
    /// Damping for step-size control d_σ.
    pub d_sigma: f64,
    /// Expected value ||N(0,I)|| ≈ sqrt(n).
    pub chi_n: f64,
    /// Generation counter.
    pub generation: u64,
}
impl CmaEs {
    /// Initialise CMA-ES at `mean` with step size `sigma`.
    ///
    /// Uses the default Hansen (2016) parameter settings.
    pub fn new(mean: Vec<f64>, sigma: f64) -> Self {
        let n = mean.len();
        assert!(n > 0, "CMA-ES requires at least 1 dimension");
        let lam = 4 + (3.0 * (n as f64).ln()).floor() as usize;
        let mu = lam / 2;
        let raw: Vec<f64> = (0..mu)
            .map(|i| ((mu as f64 + 0.5).ln() - (i as f64 + 1.0).ln()).max(0.0))
            .collect();
        let w_sum: f64 = raw.iter().sum();
        let weights: Vec<f64> = raw.iter().map(|&w| w / w_sum).collect();
        let mu_eff: f64 = 1.0 / weights.iter().map(|&w| w * w).sum::<f64>();
        let c_c = (4.0 + mu_eff / n as f64) / (n as f64 + 4.0 + 2.0 * mu_eff / n as f64);
        let c_sigma = (mu_eff + 2.0) / (n as f64 + mu_eff + 5.0);
        let c1 = 2.0 / ((n as f64 + 1.3).powi(2) + mu_eff);
        let c_mu = (2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((n as f64 + 2.0).powi(2) + mu_eff))
            .min(1.0 - c1);
        let d_sigma =
            1.0 + 2.0 * (0.0_f64.max((mu_eff - 1.0) / (n as f64 + 1.0) - 1.0)).sqrt() + c_sigma;
        let chi_n =
            (n as f64).sqrt() * (1.0 - 1.0 / (4.0 * n as f64) + 1.0 / (21.0 * (n as f64).powi(2)));
        let mut cov = vec![0.0_f64; n * n];
        for i in 0..n {
            cov[i * n + i] = 1.0;
        }
        Self {
            mean,
            sigma,
            lam,
            mu,
            weights,
            cov,
            p_c: vec![0.0; n],
            p_sigma: vec![0.0; n],
            c_c,
            c_sigma,
            c1,
            c_mu,
            d_sigma,
            chi_n,
            generation: 0,
        }
    }
    /// Perform one CMA-ES generation.
    ///
    /// Samples `λ` candidate solutions from the current multivariate Gaussian
    /// `N(mean, σ² * C)`, evaluates the objective, ranks the candidates,
    /// and updates the mean, covariance paths, and step size.
    ///
    /// Returns the best objective value in this generation.
    pub fn step<F>(&mut self, f: &F) -> f64
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = self.mean.len();
        let mut rng = rand::rng();
        let sqrt_c = cholesky_lower(&self.cov, n);
        let mut samples: Vec<Vec<f64>> = (0..self.lam)
            .map(|_| {
                let z: Vec<f64> = (0..n)
                    .map(|_| {
                        let u1: f64 = rng.random::<f64>().max(f64::EPSILON);
                        let u2: f64 = rng.random::<f64>();
                        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
                    })
                    .collect();
                let mut y = self.mean.clone();
                for i in 0..n {
                    for j in 0..=i {
                        y[i] += self.sigma * sqrt_c[i * n + j] * z[j];
                    }
                }
                y
            })
            .collect();
        let mut fitvals: Vec<(f64, usize)> =
            samples.iter().enumerate().map(|(i, x)| (f(x), i)).collect();
        fitvals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let best_f = fitvals[0].0;
        let old_mean = self.mean.clone();
        let mut new_mean = vec![0.0_f64; n];
        for (rank, &(_, idx)) in fitvals.iter().take(self.mu).enumerate() {
            let w = self.weights[rank];
            for d in 0..n {
                new_mean[d] += w * samples[idx][d];
            }
        }
        self.mean = new_mean;
        let mean_delta: Vec<f64> = (0..n)
            .map(|i| (self.mean[i] - old_mean[i]) / self.sigma)
            .collect();
        let l_inv_delta = forward_solve_lower(&sqrt_c, &mean_delta, n);
        let c_s_comp = (self.c_sigma * (2.0 - self.c_sigma)).sqrt();
        let mu_eff: f64 = 1.0 / self.weights.iter().map(|&w| w * w).sum::<f64>();
        let mu_eff_sqrt = mu_eff.sqrt();
        for (ps, &ld) in self.p_sigma.iter_mut().zip(l_inv_delta.iter()) {
            *ps = (1.0 - self.c_sigma) * *ps + c_s_comp * mu_eff_sqrt * ld;
        }
        let p_sigma_norm = self.p_sigma.iter().map(|v| v * v).sum::<f64>().sqrt();
        self.sigma *= ((self.c_sigma / self.d_sigma) * (p_sigma_norm / self.chi_n - 1.0)).exp();
        let h_sigma = if p_sigma_norm < (1.4 + 2.0 / (n as f64 + 1.0)) * self.chi_n {
            1.0_f64
        } else {
            0.0_f64
        };
        let c_c_comp = (self.c_c * (2.0 - self.c_c)).sqrt();
        for (pc, &md) in self.p_c.iter_mut().zip(mean_delta.iter()) {
            *pc = (1.0 - self.c_c) * *pc + h_sigma * c_c_comp * mu_eff_sqrt * md;
        }
        let mut rank_mu_sum = vec![0.0_f64; n * n];
        for (rank, &(_, idx)) in fitvals.iter().take(self.mu).enumerate() {
            let w = self.weights[rank];
            let delta: Vec<f64> = (0..n)
                .map(|i| (samples[idx][i] - old_mean[i]) / self.sigma)
                .collect();
            for i in 0..n {
                for j in 0..n {
                    rank_mu_sum[i * n + j] += w * delta[i] * delta[j];
                }
            }
        }
        let delta_h = (1.0 - h_sigma) * self.c_c * (2.0 - self.c_c);
        for i in 0..n {
            for j in 0..n {
                let rank1 = self.p_c[i] * self.p_c[j];
                self.cov[i * n + j] = (1.0 - self.c1 - self.c_mu) * self.cov[i * n + j]
                    + self.c1 * (rank1 + delta_h * self.cov[i * n + j])
                    + self.c_mu * rank_mu_sum[i * n + j];
            }
        }
        for i in 0..n {
            for j in 0..i {
                let avg = (self.cov[i * n + j] + self.cov[j * n + i]) / 2.0;
                self.cov[i * n + j] = avg;
                self.cov[j * n + i] = avg;
            }
        }
        samples.sort_by(|a, b| f(a).partial_cmp(&f(b)).unwrap_or(std::cmp::Ordering::Equal));
        self.generation += 1;
        best_f
    }
}
/// Result of an iterative optimization run.
#[derive(Debug, Clone)]
pub struct OptResult {
    /// The best parameter vector found.
    pub x: Vec<f64>,
    /// The objective function value at `x`.
    pub f_val: f64,
    /// Number of iterations performed.
    pub n_iter: u32,
    /// Whether the algorithm converged within the given tolerance.
    pub converged: bool,
}
/// Simple trust-region optimizer using a dogleg step.
///
/// Maintains a trust-region radius `delta` and adjusts it based on the ratio
/// of actual to predicted reduction.
pub struct TrustRegion {
    /// Current trust-region radius.
    pub delta: f64,
    /// Acceptance threshold η ∈ (0, 0.25).
    pub eta: f64,
}
impl TrustRegion {
    /// Create a new `TrustRegion` with given radius and acceptance threshold.
    pub fn new(delta: f64, eta: f64) -> Self {
        Self { delta, eta }
    }
    /// Perform one trust-region step.
    ///
    /// Given the current point `x`, computes the gradient and a Cauchy/dogleg
    /// step within the trust region.  Updates `self.delta` based on the
    /// actual-to-predicted reduction ratio.
    ///
    /// Returns the new point (may equal `x` if step is rejected).
    pub fn step<F, G>(&mut self, f: F, grad: G, x: &[f64]) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
        G: Fn(&[f64]) -> Vec<f64>,
    {
        let n = x.len();
        let fx = f(x);
        let g = grad(x);
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < 1e-15 {
            return x.to_vec();
        }
        let cauchy_scale = self.delta / gnorm;
        let p: Vec<f64> = g.iter().map(|gi| -cauchy_scale * gi).collect();
        let p_norm: f64 = p.iter().map(|v| v * v).sum::<f64>().sqrt();
        let p = if p_norm > self.delta {
            p.iter()
                .map(|pi| pi * self.delta / p_norm)
                .collect::<Vec<_>>()
        } else {
            p
        };
        let x_new: Vec<f64> = (0..n).map(|i| x[i] + p[i]).collect();
        let fx_new = f(&x_new);
        let pred_red: f64 = -g.iter().zip(p.iter()).map(|(gi, pi)| gi * pi).sum::<f64>();
        let actual_red = fx - fx_new;
        let rho = if pred_red.abs() < 1e-15 {
            1.0
        } else {
            actual_red / pred_red
        };
        if rho < 0.25 {
            self.delta *= 0.25;
        } else if rho > 0.75 {
            self.delta = (2.0 * self.delta).min(1e6);
        }
        if rho > self.eta { x_new } else { x.to_vec() }
    }
}
