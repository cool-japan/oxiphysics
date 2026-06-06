//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Thermodynamic integration with statistical error propagation.
pub struct ThermoIntegration {
    /// Lambda grid points in \[0, 1\].
    pub lambdas: Vec<f64>,
    /// Mean ⟨dU/dλ⟩ at each lambda point.
    pub du_dl: Vec<f64>,
    /// Standard error of ⟨dU/dλ⟩ at each lambda point.
    pub du_dl_err: Vec<f64>,
}
impl ThermoIntegration {
    /// Create a new `ThermoIntegration` from raw sample arrays at each lambda.
    ///
    /// # Arguments
    /// * `lambdas`  – lambda grid (ascending, in \[0,1\]).
    /// * `samples`  – `samples[i]` = raw dU/dλ samples at lambda\[i\].
    pub fn from_samples(lambdas: Vec<f64>, samples: Vec<Vec<f64>>) -> Self {
        let n = lambdas.len().min(samples.len());
        let mut du_dl = Vec::with_capacity(n);
        let mut du_dl_err = Vec::with_capacity(n);
        for s in samples.iter().take(n) {
            let (mean, stderr) = block_mean_stderr(s);
            du_dl.push(mean);
            du_dl_err.push(stderr);
        }
        ThermoIntegration {
            lambdas,
            du_dl,
            du_dl_err,
        }
    }
    /// Compute the TI free energy estimate using the trapezoidal rule.
    pub fn compute_free_energy(&self) -> f64 {
        thermodynamic_integration(&self.lambdas, &self.du_dl)
    }
    /// Compute the statistical error estimate of the TI integral.
    ///
    /// Propagates the standard errors through the trapezoidal rule in
    /// quadrature:
    ///
    /// ```text
    /// σ²(ΔF) = Σ_i [ ½(σᵢ + σᵢ₊₁) Δλ ]²
    /// ```
    pub fn compute_error_estimate(&self) -> f64 {
        ti_uncertainty(&self.lambdas, &self.du_dl_err)
    }
    /// Number of lambda windows.
    pub fn n_windows(&self) -> usize {
        self.lambdas.len()
    }
}
/// Extended umbrella-sampling window descriptor.
pub struct UmbrellaWindow {
    /// Equilibrium position of the harmonic bias (same units as ξ).
    pub xi_ref: f64,
    /// Spring constant (energy / length²).
    pub k_spring: f64,
    /// Accumulated ξ samples.
    pub samples: Vec<f64>,
    /// Target temperature (K).
    pub temperature: f64,
}
impl UmbrellaWindow {
    /// Create a new umbrella-sampling window.
    pub fn new(xi_ref: f64, k_spring: f64, temperature: f64) -> Self {
        Self {
            xi_ref,
            k_spring,
            samples: Vec::new(),
            temperature,
        }
    }
    /// Add a sample.
    pub fn push(&mut self, xi: f64) {
        self.samples.push(xi);
    }
    /// Compute the harmonic bias potential (energy).
    pub fn bias(&self, xi: f64) -> f64 {
        0.5 * self.k_spring * (xi - self.xi_ref).powi(2)
    }
    /// Mean ξ from samples.
    pub fn mean_xi(&self) -> Option<f64> {
        let n = self.samples.len();
        if n == 0 {
            return None;
        }
        Some(self.samples.iter().sum::<f64>() / n as f64)
    }
    /// Variance of ξ samples.
    pub fn var_xi(&self) -> Option<f64> {
        let n = self.samples.len();
        if n < 2 {
            return None;
        }
        let mean = self.mean_xi()?;
        Some(
            self.samples
                .iter()
                .map(|&x| (x - mean).powi(2))
                .sum::<f64>()
                / (n - 1) as f64,
        )
    }
    /// Effective spring constant from fluctuations: k_eff = kT / ⟨δξ²⟩.
    pub fn effective_spring_constant(&self, kt: f64) -> Option<f64> {
        let var = self.var_xi()?;
        if var < 1e-300 {
            return None;
        }
        Some(kt / var)
    }
}
/// Two-dimensional WHAM window for a 2-D PMF calculation.
pub struct WhamWindow2D {
    /// Reference positions of the bias (xi1_ref, xi2_ref).
    pub xi_ref: [f64; 2],
    /// Spring constants (k1, k2).
    pub k_bias: [f64; 2],
    /// Samples as (xi1, xi2) pairs.
    pub samples: Vec<[f64; 2]>,
    /// Current free energy offset estimate (dimensionless).
    pub f_estimate: f64,
}
impl WhamWindow2D {
    /// Create a new 2D WHAM window.
    pub fn new(xi_ref: [f64; 2], k_bias: [f64; 2]) -> Self {
        Self {
            xi_ref,
            k_bias,
            samples: Vec::new(),
            f_estimate: 0.0,
        }
    }
    /// Harmonic bias energy at 2-D position `xi` (dimensionless).
    pub fn bias_energy_2d(&self, xi: [f64; 2], kt: f64) -> f64 {
        let e1 = 0.5 * self.k_bias[0] * (xi[0] - self.xi_ref[0]).powi(2);
        let e2 = 0.5 * self.k_bias[1] * (xi[1] - self.xi_ref[1]).powi(2);
        (e1 + e2) / kt
    }
    /// Add a sample.
    pub fn add_sample_2d(&mut self, xi: [f64; 2]) {
        self.samples.push(xi);
    }
}
/// A single simulation state (λ-window or umbrella window) for MBAR.
///
/// Stores the reduced potential energies of the samples from this state
/// evaluated at every other state (the "u_kn" matrix).
#[derive(Debug, Clone)]
pub struct MbarState {
    /// State index.
    pub index: usize,
    /// Number of samples from this state.
    pub n_samples: usize,
    /// Reduced potential energy u_l(x_kn): potential of each sample from THIS state
    /// evaluated at EVERY state l.  Shape: n_states × n_samples.
    pub u_kl: Vec<Vec<f64>>,
}
impl MbarState {
    /// Create a new MBAR state.
    pub fn new(index: usize) -> Self {
        Self {
            index,
            n_samples: 0,
            u_kl: Vec::new(),
        }
    }
}
/// Kullback-Leibler and Jensen-Shannon divergence helper.
///
/// Both `p` and `q` are assumed to be normalized probability vectors
/// (each sums to 1, all entries ≥ 0).
pub struct KlDivergence;
impl KlDivergence {
    /// Kullback-Leibler divergence D_KL(P ‖ Q) = Σ p_i · ln(p_i / q_i).
    ///
    /// Terms where p_i = 0 are skipped (0 · ln(0) = 0 by convention).
    /// Returns `f64::INFINITY` if any q_i = 0 while p_i > 0.
    pub fn kl_divergence(p: &[f64], q: &[f64]) -> f64 {
        let n = p.len().min(q.len());
        let mut sum = 0.0;
        for i in 0..n {
            if p[i] <= 0.0 {
                continue;
            }
            if q[i] <= 0.0 {
                return f64::INFINITY;
            }
            sum += p[i] * (p[i] / q[i]).ln();
        }
        sum
    }
    /// Jensen-Shannon divergence.
    ///
    /// JS(P ‖ Q) = ½ D_KL(P ‖ M) + ½ D_KL(Q ‖ M),  M = ½(P + Q).
    ///
    /// Always finite and symmetric; values in \[0, ln 2\].
    pub fn js_divergence(p: &[f64], q: &[f64]) -> f64 {
        let n = p.len().min(q.len());
        let m: Vec<f64> = (0..n).map(|i| 0.5 * (p[i] + q[i])).collect();
        0.5 * Self::kl_divergence(&p[..n], &m) + 0.5 * Self::kl_divergence(&q[..n], &m)
    }
}
/// A single λ-window in an alchemical FEP calculation.
///
/// Stores the λ value, the ΔU samples collected at this λ, and optional
/// overlap statistics with the adjacent window.
#[derive(Debug, Clone)]
pub struct LambdaWindow {
    /// λ value (0 = fully coupled, 1 = fully decoupled, or vice-versa).
    pub lambda: f64,
    /// ΔU samples from state λ evaluated at state λ+1.
    pub du_forward: Vec<f64>,
    /// ΔU samples from state λ evaluated at state λ-1.
    pub du_backward: Vec<f64>,
}
impl LambdaWindow {
    /// Create a new empty λ-window.
    pub fn new(lambda: f64) -> Self {
        Self {
            lambda,
            du_forward: Vec::new(),
            du_backward: Vec::new(),
        }
    }
    /// Add a forward ΔU sample.
    pub fn add_forward(&mut self, du: f64) {
        self.du_forward.push(du);
    }
    /// Add a backward ΔU sample.
    pub fn add_backward(&mut self, du: f64) {
        self.du_backward.push(du);
    }
    /// Zwanzig estimate for this window (forward direction).
    pub fn zwanzig_forward(&self, kt: f64) -> f64 {
        fep_zwanzig(&self.du_forward, kt)
    }
    /// BAR estimate between this window and the next.
    pub fn bar_estimate(&self, next_window: &LambdaWindow, kt: f64) -> f64 {
        fep_bar(&self.du_forward, &next_window.du_backward, kt, 1e-10)
    }
    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.du_forward.len().max(self.du_backward.len())
    }
}
/// Staged alchemical transformation result.
#[derive(Debug, Clone)]
pub struct AlchemicalResult {
    /// Free energy of decharging stage.
    pub decharge_df: f64,
    /// Free energy of van-der-Waals decoupling stage.
    pub vdw_df: f64,
    /// Total ΔF = decharge + vdW.
    pub total_df: f64,
}
/// MBAR estimator with bootstrap uncertainty quantification.
pub struct Mbar {
    /// K × N_total reduced potential matrix.
    pub u_kn: Vec<Vec<f64>>,
    /// Number of samples from each state.
    pub n_k: Vec<usize>,
    /// MBAR free energies (normalised so f\[0\] = 0).
    pub free_energies: Vec<f64>,
}
impl Mbar {
    /// Create and solve a new MBAR instance.
    ///
    /// # Arguments
    /// * `u_kn`     – K × N_total reduced potential matrix.
    /// * `n_k`      – number of samples from each state.
    /// * `tol`      – convergence tolerance.
    /// * `max_iter` – maximum MBAR iterations.
    pub fn new(u_kn: Vec<Vec<f64>>, n_k: Vec<usize>, tol: f64, max_iter: usize) -> Self {
        let free_energies = mbar_free_energies(&u_kn, &n_k, tol, max_iter);
        Self {
            u_kn,
            n_k,
            free_energies,
        }
    }
    /// Estimate the uncertainty of the MBAR free energies via bootstrapping.
    ///
    /// Resamples the columns of `u_kn` (with replacement) `n_bootstrap` times,
    /// re-runs MBAR, and returns the standard deviation of the bootstrap
    /// free energy estimates for each state.
    ///
    /// # Arguments
    /// * `n_bootstrap` – number of bootstrap replicates.
    /// * `seed_offset`  – added to the replicate index to vary random draws.
    pub fn compute_uncertainty(&self, n_bootstrap: usize, seed_offset: u64) -> Vec<f64> {
        let k = self.free_energies.len();
        if k == 0 || n_bootstrap == 0 {
            return vec![0.0; k];
        }
        let n_total: usize = self.n_k.iter().sum();
        if n_total == 0 {
            return vec![0.0; k];
        }
        let tol = 1e-6;
        let max_iter = 1000;
        let mut boot_estimates: Vec<Vec<f64>> = Vec::with_capacity(n_bootstrap);
        for rep in 0..n_bootstrap {
            let mut lcg: u64 = (rep as u64)
                .wrapping_add(seed_offset)
                .wrapping_add(1)
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let mut boot_u: Vec<Vec<f64>> = (0..k).map(|_| Vec::with_capacity(n_total)).collect();
            for _ in 0..n_total {
                lcg = lcg
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let idx = (lcg >> 33) as usize % n_total;
                for (ki, bu) in boot_u.iter_mut().enumerate().take(k) {
                    bu.push(self.u_kn[ki].get(idx).copied().unwrap_or(0.0));
                }
            }
            let f_boot = mbar_free_energies(&boot_u, &self.n_k, tol, max_iter);
            boot_estimates.push(f_boot);
        }
        (0..k)
            .map(|ki| {
                let vals: Vec<f64> = boot_estimates
                    .iter()
                    .filter_map(|v| v.get(ki).copied())
                    .collect();
                let n = vals.len();
                if n < 2 {
                    return 0.0;
                }
                let mean = vals.iter().sum::<f64>() / n as f64;
                let var = vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
                var.sqrt()
            })
            .collect()
    }
}
/// A single umbrella-sampling window used in WHAM.
///
/// The window contains a harmonic bias potential applied at reference position
/// `xi_ref` with spring constant `k_bias`.
#[derive(Debug, Clone)]
pub struct WhamWindow {
    /// Reference position of the bias potential (reaction-coordinate units).
    pub xi_ref: f64,
    /// Spring constant of the harmonic bias (energy/length²).
    pub k_bias: f64,
    /// Samples of the collective variable (unbiased positions) from this window.
    pub samples: Vec<f64>,
    /// Current estimate of the free energy offset `f_i` (dimensionless, kT units).
    pub f_estimate: f64,
}
impl WhamWindow {
    /// Create a new WHAM window.
    pub fn new(xi_ref: f64, k_bias: f64) -> Self {
        Self {
            xi_ref,
            k_bias,
            samples: Vec::new(),
            f_estimate: 0.0,
        }
    }
    /// Harmonic bias energy at position `xi` (dimensionless / kT).
    pub fn bias_energy(&self, xi: f64, kt: f64) -> f64 {
        0.5 * self.k_bias * (xi - self.xi_ref).powi(2) / kt
    }
    /// Add a sample to this window.
    pub fn add_sample(&mut self, xi: f64) {
        self.samples.push(xi);
    }
}
/// Free energy profile (PMF) along a collective variable.
pub struct FreeEnergyProfile {
    /// Collective variable values (bin centers or sample points).
    pub xi: Vec<f64>,
    /// Free energy F(ξ) values (in units consistent with kT).
    pub pmf: Vec<f64>,
}
impl FreeEnergyProfile {
    /// Create a `FreeEnergyProfile` from (xi, pmf) data.
    pub fn new(xi: Vec<f64>, pmf: Vec<f64>) -> Self {
        Self { xi, pmf }
    }
    /// Compute the barrier height from the PMF.
    ///
    /// The barrier height is:
    ///
    /// ```text
    /// ΔF‡ = F_max - F_min (over the range [xi_start, xi_end])
    /// ```
    ///
    /// # Arguments
    /// * `xi_start` – lower bound of the collective variable range to search.
    /// * `xi_end`   – upper bound.
    ///
    /// Returns `None` if no data points lie in the range or PMF is empty.
    pub fn compute_barrier(&self, xi_start: f64, xi_end: f64) -> Option<f64> {
        let values: Vec<f64> = self
            .xi
            .iter()
            .zip(self.pmf.iter())
            .filter(|(x, _)| **x >= xi_start && **x <= xi_end)
            .map(|(_, f)| *f)
            .collect();
        if values.is_empty() {
            return None;
        }
        let f_min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let f_max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Some(f_max - f_min)
    }
    /// Find the position of the PMF maximum (transition state estimate).
    pub fn transition_state(&self) -> Option<f64> {
        let (max_idx, _) = self
            .pmf
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        self.xi.get(max_idx).copied()
    }
    /// Find the position of the PMF minimum (reactant/product well).
    pub fn minimum_position(&self) -> Option<f64> {
        let (min_idx, _) = self
            .pmf
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        self.xi.get(min_idx).copied()
    }
    /// Number of data points.
    pub fn len(&self) -> usize {
        self.xi.len()
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.xi.is_empty()
    }
}
