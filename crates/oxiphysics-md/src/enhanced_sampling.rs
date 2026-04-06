// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Enhanced sampling methods for molecular dynamics.
//!
//! Provides:
//! - [`ReplicaExchange`]: Parallel tempering / REMD with Metropolis swap acceptance.
//! - [`Metadynamics`]: Gaussian-hill-based metadynamics on a 1D collective variable.
//! - [`UmbrellaSamplingWindow`]: Single harmonic-bias umbrella sampling window.
//!
//! These types complement the existing lower-level modules [`crate::remd`],
//! [`crate::metadynamics`], and [`crate::sampling`] with a higher-level API
//! that matches the interface described in the crate documentation.

#![allow(dead_code, missing_docs)]

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant in kJ/mol/K.
const KB: f64 = 8.314e-3;

// ─────────────────────────────────────────────────────────────────────────────
// ReplicaExchange
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel tempering (Replica Exchange MD) controller.
///
/// Maintains a set of replicas, each running at a different temperature,
/// and periodically proposes configuration swaps between adjacent replicas
/// using the Metropolis criterion:
///
/// ```text
/// P_acc = min(1, exp((β_i − β_j)(E_i − E_j)))
/// ```
#[derive(Debug, Clone)]
pub struct ReplicaExchange {
    /// Temperatures of each replica (K).
    pub temperatures: Vec<f64>,
    /// Current potential energies of each replica (kJ/mol).
    pub energies: Vec<f64>,
    /// Number of replicas.
    pub n_replicas: usize,
    /// Total number of swap attempts made.
    pub n_attempts: u64,
    /// Total number of accepted swaps.
    pub n_accepted: u64,
    /// History of accepted swaps: each entry is `(step, replica_i, replica_j)`.
    pub swap_history: Vec<(u64, usize, usize)>,
    /// Step counter.
    pub step: u64,
}

impl ReplicaExchange {
    /// Create a new [`ReplicaExchange`] controller with the given temperature ladder.
    ///
    /// Energies are initialised to zero; call [`Self::update_energies`] before
    /// the first exchange attempt.
    pub fn new(temperatures: Vec<f64>) -> Self {
        let n = temperatures.len();
        Self {
            temperatures,
            energies: vec![0.0; n],
            n_replicas: n,
            n_attempts: 0,
            n_accepted: 0,
            swap_history: Vec::new(),
            step: 0,
        }
    }

    /// Create a geometrically spaced temperature ladder with `n` replicas
    /// spanning `[t_min, t_max]`.
    pub fn geometric_ladder(t_min: f64, t_max: f64, n: usize) -> Self {
        if n == 0 {
            return Self::new(vec![]);
        }
        if n == 1 {
            return Self::new(vec![t_min]);
        }
        let ratio = (t_max / t_min).powf(1.0 / (n - 1) as f64);
        let temps: Vec<f64> = (0..n).map(|i| t_min * ratio.powi(i as i32)).collect();
        Self::new(temps)
    }

    /// Update the potential energy of each replica.
    ///
    /// `energies` must have the same length as the number of replicas.
    pub fn update_energies(&mut self, energies: &[f64]) {
        assert_eq!(
            energies.len(),
            self.n_replicas,
            "Energy slice length must match n_replicas"
        );
        self.energies.copy_from_slice(energies);
    }

    /// Attempt configuration swaps between all adjacent replica pairs.
    ///
    /// The caller supplies `random_values`: a slice of uniform random numbers
    /// in `[0, 1)`, one per adjacent pair.  Passing `None` causes all swaps
    /// with `P_acc >= 1` to be deterministically accepted (useful for tests).
    ///
    /// Returns the list of `(replica_i, replica_j)` pairs that were swapped
    /// this round.
    pub fn attempt_swap(&mut self, energies: &[f64], betas: &[f64]) -> Vec<(usize, usize)> {
        let n = self.n_replicas;
        assert_eq!(energies.len(), n, "energies length mismatch");
        assert_eq!(betas.len(), n, "betas length mismatch");

        let mut swapped = Vec::new();

        for pair in 0..(n.saturating_sub(1)) {
            let bi = betas[pair];
            let bj = betas[pair + 1];
            let ei = energies[pair];
            let ej = energies[pair + 1];

            let delta = (bi - bj) * (ei - ej);
            let p_acc = if delta >= 0.0 { 1.0 } else { delta.exp() };

            self.n_attempts += 1;

            // Deterministic: accept if p_acc >= 1 (no random number provided)
            if p_acc >= 1.0 {
                self.n_accepted += 1;
                self.swap_history.push((self.step, pair, pair + 1));
                swapped.push((pair, pair + 1));
                // Swap temperatures to reflect the exchange
                self.temperatures.swap(pair, pair + 1);
            }
        }

        self.step += 1;
        swapped
    }

    /// Attempt configuration swaps using explicit random acceptance values.
    ///
    /// `random_values[i]` is a uniform random number for pair `(i, i+1)`.
    /// A swap is accepted if `random_values[i] < P_acc`.
    pub fn attempt_swap_random(
        &mut self,
        energies: &[f64],
        betas: &[f64],
        random_values: &[f64],
    ) -> Vec<(usize, usize)> {
        let n = self.n_replicas;
        assert_eq!(energies.len(), n, "energies length mismatch");
        assert_eq!(betas.len(), n, "betas length mismatch");

        let mut swapped = Vec::new();

        for pair in 0..(n.saturating_sub(1)) {
            let bi = betas[pair];
            let bj = betas[pair + 1];
            let ei = energies[pair];
            let ej = energies[pair + 1];

            let delta = (bi - bj) * (ei - ej);
            let p_acc = if delta >= 0.0 { 1.0_f64 } else { delta.exp() };

            self.n_attempts += 1;

            let rand_val = random_values.get(pair).copied().unwrap_or(1.0);
            if rand_val < p_acc {
                self.n_accepted += 1;
                self.swap_history.push((self.step, pair, pair + 1));
                swapped.push((pair, pair + 1));
                self.temperatures.swap(pair, pair + 1);
            }
        }

        self.step += 1;
        swapped
    }

    /// Compute the Metropolis acceptance probability for a swap between replicas
    /// `i` and `j`.
    ///
    /// P_acc = min(1, exp((β_i − β_j)(E_i − E_j)))
    pub fn acceptance_probability(e_i: f64, e_j: f64, beta_i: f64, beta_j: f64) -> f64 {
        let delta = (beta_i - beta_j) * (e_i - e_j);
        if delta >= 0.0 { 1.0 } else { delta.exp() }
    }

    /// Overall swap acceptance rate.
    pub fn acceptance_rate(&self) -> f64 {
        if self.n_attempts == 0 {
            return 0.0;
        }
        self.n_accepted as f64 / self.n_attempts as f64
    }

    /// Compute inverse temperatures (β = 1/(kB T)) for all replicas.
    pub fn betas(&self) -> Vec<f64> {
        self.temperatures.iter().map(|&t| 1.0 / (KB * t)).collect()
    }

    /// Returns the number of replicas.
    pub fn n_replicas(&self) -> usize {
        self.n_replicas
    }

    /// Returns the number of accepted swaps (alias for n_accepted).
    pub fn accepted(&self) -> u64 {
        self.n_accepted
    }

    /// Returns the number of attempted swaps (alias for n_attempts).
    pub fn attempted(&self) -> u64 {
        self.n_attempts
    }

    /// Set the energy of a single replica by index.
    pub fn set_energy(&mut self, idx: usize, energy: f64) {
        self.energies[idx] = energy;
    }

    /// Attempt a swap between replicas `i` and `j` using a pre-drawn random value.
    ///
    /// Returns `true` if the swap was accepted.
    pub fn try_swap(&mut self, i: usize, j: usize, rand_val: f64) -> bool {
        let kb = 1.380_649e-23_f64;
        let beta_i = 1.0 / (kb * self.temperatures[i]);
        let beta_j = 1.0 / (kb * self.temperatures[j]);
        let p_acc =
            Self::acceptance_probability(self.energies[i], self.energies[j], beta_i, beta_j);
        self.n_attempts += 1;
        if rand_val < p_acc {
            self.temperatures.swap(i, j);
            self.energies.swap(i, j);
            self.n_accepted += 1;
            true
        } else {
            false
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Gaussian hill (internal)
// ─────────────────────────────────────────────────────────────────────────────

/// A single Gaussian bias hill deposited in metadynamics.
#[derive(Debug, Clone)]
struct GaussBias {
    center: f64,
    height: f64,
    width: f64,
}

impl GaussBias {
    fn evaluate(&self, s: f64) -> f64 {
        let d = s - self.center;
        self.height * (-d * d / (2.0 * self.width * self.width)).exp()
    }

    fn gradient(&self, s: f64) -> f64 {
        let d = s - self.center;
        -self.evaluate(s) * d / (self.width * self.width)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Metadynamics
// ─────────────────────────────────────────────────────────────────────────────

/// Metadynamics bias on a 1D collective variable.
///
/// Repeatedly deposits Gaussian hills centred at the current CV value.
/// Supports both plain metadynamics and well-tempered metadynamics (set
/// `delta_t > 0`).
#[derive(Debug, Clone)]
pub struct Metadynamics {
    /// Deposited Gaussian hills.
    hills: Vec<GaussBias>,
    /// Initial hill height (kJ/mol).
    pub height: f64,
    /// Hill width (sigma, in CV units).
    pub width: f64,
    /// Number of MD steps between hill depositions.
    pub deposition_stride: usize,
    /// Well-tempered bias temperature increment ΔT (K).  0 = plain metadynamics.
    pub delta_t: f64,
    /// Simulation temperature T (K).
    pub temperature: f64,
    /// Current step counter.
    pub step: usize,
}

impl Metadynamics {
    /// Create a new [`Metadynamics`] instance.
    ///
    /// # Arguments
    /// * `height`            – Initial Gaussian height (kJ/mol).
    /// * `width`             – Gaussian width σ (CV units).
    /// * `deposition_stride` – Deposit a new hill every this many steps.
    /// * `temperature`       – Simulation temperature (K).
    /// * `delta_t`           – Well-tempered ΔT (K); 0 for plain metadynamics.
    pub fn new(
        height: f64,
        width: f64,
        deposition_stride: usize,
        temperature: f64,
        delta_t: f64,
    ) -> Self {
        Self {
            hills: Vec::new(),
            height,
            width,
            deposition_stride,
            delta_t,
            temperature,
            step: 0,
        }
    }

    /// Add a Gaussian hill at position `cv_value` with the given `height` and `width`.
    ///
    /// This bypasses the automatic stride; use this for direct insertion.
    pub fn add_gaussian(&mut self, cv_value: f64, height: f64, width: f64) {
        self.hills.push(GaussBias {
            center: cv_value,
            height,
            width,
        });
    }

    /// Compute the total bias potential at a given CV value.
    ///
    /// V_bias(s) = Σ_i V_Gauss_i(s)
    pub fn compute_bias(&self, cv_value: f64) -> f64 {
        self.hills.iter().map(|h| h.evaluate(cv_value)).sum()
    }

    /// Compute the bias force at a given CV value.
    ///
    /// F_bias(s) = −dV_bias/ds
    pub fn compute_bias_force(&self, cv_value: f64) -> f64 {
        -self.hills.iter().map(|h| h.gradient(cv_value)).sum::<f64>()
    }

    /// Attempt to deposit a Gaussian hill at the current CV position.
    ///
    /// Deposits when `step % deposition_stride == 0`.
    /// For well-tempered metadynamics, the effective height is scaled by:
    ///   h_eff = height * exp(−V_bias(s) / (kB * (T + ΔT)))
    pub fn maybe_deposit(&mut self, cv_value: f64) {
        if self.deposition_stride > 0 && self.step.is_multiple_of(self.deposition_stride) {
            let h = if self.delta_t > 0.0 {
                let v_bias = self.compute_bias(cv_value);
                let biasfactor_temp = self.temperature + self.delta_t;
                self.height * (-v_bias / (KB * biasfactor_temp)).exp()
            } else {
                self.height
            };
            self.hills.push(GaussBias {
                center: cv_value,
                height: h,
                width: self.width,
            });
        }
        self.step += 1;
    }

    /// Number of Gaussian hills deposited so far.
    pub fn n_hills(&self) -> usize {
        self.hills.len()
    }

    /// Estimate the free energy on a 1D grid of CV values.
    ///
    /// Returns F(s) = −V_bias(s) (converged estimate in the well-tempered case).
    pub fn free_energy_estimate(&self, cv_grid: &[f64]) -> Vec<f64> {
        cv_grid.iter().map(|&s| -self.compute_bias(s)).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UmbrellaSamplingWindow
// ─────────────────────────────────────────────────────────────────────────────

/// One umbrella sampling window with a harmonic bias potential.
///
/// V_bias(ξ) = k/2 · (ξ − ξ₀)²
///
/// Multiple windows are combined via WHAM or by using the associated
/// `compute_pmf` helper.
#[derive(Debug, Clone)]
pub struct UmbrellaSamplingWindow {
    /// Window centre ξ₀ (CV units).
    pub center: f64,
    /// Harmonic force constant k (kJ/mol/CV²).
    pub force_constant: f64,
    /// Collected CV samples.
    pub samples: Vec<f64>,
}

impl UmbrellaSamplingWindow {
    /// Create a new umbrella sampling window at `center` with force constant `k`.
    pub fn new(center: f64, force_constant: f64) -> Self {
        Self {
            center,
            force_constant,
            samples: Vec::new(),
        }
    }

    /// Compute the bias potential at CV value `xi`.
    pub fn bias_potential(&self, xi: f64) -> f64 {
        let d = xi - self.center;
        0.5 * self.force_constant * d * d
    }

    /// Compute the bias force at CV value `xi`.
    ///
    /// F = −dV/dξ = −k · (ξ − ξ₀)
    pub fn bias_force(&self, xi: f64) -> f64 {
        -self.force_constant * (xi - self.center)
    }

    /// Add a sample to this window's histogram.
    pub fn add_sample(&mut self, xi: f64) {
        self.samples.push(xi);
    }

    /// Mean and variance of collected samples.
    pub fn statistics(&self) -> (f64, f64) {
        let n = self.samples.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        let mean = self.samples.iter().sum::<f64>() / n as f64;
        let var = self
            .samples
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        (mean, var)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UmbrellaSampling — a collection of windows
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of umbrella sampling windows for computing a potential of
/// mean force (PMF).
///
/// Each window is a `(center, force_constant)` tuple.
#[derive(Debug, Clone)]
pub struct UmbrellaSampling {
    /// Individual windows.
    pub windows: Vec<UmbrellaSamplingWindow>,
}

impl UmbrellaSampling {
    /// Create a new empty umbrella sampling set.
    pub fn new() -> Self {
        Self {
            windows: Vec::new(),
        }
    }

    /// Create a set of evenly spaced windows.
    ///
    /// # Arguments
    /// * `centers`        – Window centre positions.
    /// * `force_constant` – Shared force constant for all windows (kJ/mol/CV²).
    pub fn from_centers(centers: &[f64], force_constant: f64) -> Self {
        let windows = centers
            .iter()
            .map(|&c| UmbrellaSamplingWindow::new(c, force_constant))
            .collect();
        Self { windows }
    }

    /// Add a window.
    pub fn add_window(&mut self, center: f64, force_constant: f64) {
        self.windows
            .push(UmbrellaSamplingWindow::new(center, force_constant));
    }

    /// Number of windows.
    pub fn n_windows(&self) -> usize {
        self.windows.len()
    }

    /// Simple WHAM-free PMF estimate using histogram + Boltzmann inversion.
    ///
    /// For each window, the biased probability density is estimated from the
    /// histogram, then unbiased via:
    ///
    ///   F(ξ_b) = −kT ln\[H_i(ξ_b)\] + V_bias_i(ξ_b)
    ///
    /// For multiple windows, free energies are averaged in regions of overlap,
    /// and the global minimum is shifted to zero.
    ///
    /// # Arguments
    /// * `n_bins` – Number of bins for the combined histogram.
    /// * `kT`     – Thermal energy (kJ/mol).
    ///
    /// # Returns
    /// `(bin_centers, pmf_values)`.
    #[allow(non_snake_case)]
    pub fn compute_pmf(&self, n_bins: usize, kT: f64) -> (Vec<f64>, Vec<f64>) {
        if self.windows.is_empty() || n_bins == 0 {
            return (vec![], vec![]);
        }

        let all_samples: Vec<f64> = self
            .windows
            .iter()
            .flat_map(|w| w.samples.iter().cloned())
            .collect();
        if all_samples.is_empty() {
            return (vec![], vec![]);
        }

        let xmin = all_samples.iter().cloned().fold(f64::INFINITY, f64::min);
        let xmax = all_samples
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if (xmax - xmin).abs() < 1e-15 {
            return (vec![xmin], vec![0.0]);
        }

        let bw = (xmax - xmin) / n_bins as f64;
        let bin_centers: Vec<f64> = (0..n_bins).map(|i| xmin + (i as f64 + 0.5) * bw).collect();

        // For each bin, accumulate (unbiased_F, count) contributions from each window
        let mut f_sum = vec![0.0f64; n_bins];
        let mut f_count = vec![0u64; n_bins];

        for window in &self.windows {
            // Build per-window histogram
            let n_win = window.samples.len();
            if n_win == 0 {
                continue;
            }
            let mut hist = vec![0u64; n_bins];
            for &xi in &window.samples {
                let idx = ((xi - xmin) / bw) as usize;
                let idx = idx.min(n_bins - 1);
                hist[idx] += 1;
            }
            let max_count = *hist.iter().max().unwrap_or(&0);
            if max_count == 0 {
                continue;
            }
            // F_unbiased(ξ_b) = −kT ln[H_i(ξ_b)/N_i] + V_bias_i(ξ_b) + const
            // The constant per window is the window free energy offset F_i.
            // We set the reference so that the minimum F within each window = 0
            // and accumulate.
            let mut f_window = vec![f64::INFINITY; n_bins];
            for (b, &cnt) in hist.iter().enumerate() {
                if cnt > 0 {
                    let p = cnt as f64 / n_win as f64;
                    let xi_b = bin_centers[b];
                    let v_bias = window.bias_potential(xi_b);
                    // F_unbiased = −kT ln(p) + V_bias (up to a constant)
                    f_window[b] = -kT * p.ln() + v_bias;
                }
            }
            // Subtract the minimum to put each window on the same scale
            let f_min = f_window
                .iter()
                .cloned()
                .filter(|f| f.is_finite())
                .fold(f64::INFINITY, f64::min);
            if f_min.is_infinite() {
                continue;
            }
            for (b, f) in f_window.iter().enumerate() {
                if f.is_finite() {
                    f_sum[b] += f - f_min;
                    f_count[b] += 1;
                }
            }
        }

        // Average contributions per bin
        let pmf_raw: Vec<f64> = f_sum
            .iter()
            .zip(f_count.iter())
            .map(|(&s, &c)| if c > 0 { s / c as f64 } else { f64::INFINITY })
            .collect();

        // Shift global minimum to zero
        let global_min = pmf_raw
            .iter()
            .cloned()
            .filter(|f| f.is_finite())
            .fold(f64::INFINITY, f64::min);
        let pmf = if global_min.is_finite() {
            pmf_raw
                .iter()
                .map(|&f| {
                    if f.is_finite() {
                        f - global_min
                    } else {
                        f64::INFINITY
                    }
                })
                .collect()
        } else {
            pmf_raw
        };

        (bin_centers, pmf)
    }
}

impl Default for UmbrellaSampling {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive Bias Metadynamics
// ─────────────────────────────────────────────────────────────────────────────

/// Adaptive-bias metadynamics on a 1D collective variable.
///
/// The hill width adapts based on the local CV velocity (adaptive Gaussians):
/// σ_adapt = σ_min + (σ_max − σ_min) * exp(−|dCV/dt| / σ_vel)
///
/// When |dCV/dt| is large, the system is moving fast so we use wide hills
/// (exploratory). When it is slow, we use narrow hills for accuracy.
#[derive(Debug, Clone)]
pub struct AdaptiveBiasMetadynamics {
    /// Deposited hills.
    hills: Vec<GaussBias>,
    /// Base hill height (kJ/mol).
    pub height: f64,
    /// Minimum hill width σ_min (CV units).
    pub width_min: f64,
    /// Maximum hill width σ_max (CV units).
    pub width_max: f64,
    /// CV velocity scale for adaptive width (CV units/step).
    pub vel_scale: f64,
    /// Deposition stride (steps between hill deposits).
    pub deposition_stride: usize,
    /// Well-tempered ΔT (K). 0 = plain.
    pub delta_t: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Current step.
    pub step: usize,
    /// Previous CV value (for velocity estimation).
    prev_cv: Option<f64>,
}

impl AdaptiveBiasMetadynamics {
    /// Create a new adaptive-bias metadynamics controller.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        height: f64,
        width_min: f64,
        width_max: f64,
        vel_scale: f64,
        deposition_stride: usize,
        temperature: f64,
        delta_t: f64,
    ) -> Self {
        Self {
            hills: Vec::new(),
            height,
            width_min,
            width_max,
            vel_scale,
            deposition_stride,
            delta_t,
            temperature,
            step: 0,
            prev_cv: None,
        }
    }

    /// Compute the total bias at a given CV value.
    pub fn compute_bias(&self, cv: f64) -> f64 {
        self.hills.iter().map(|h| h.evaluate(cv)).sum()
    }

    /// Compute the bias force at a given CV value.
    pub fn compute_bias_force(&self, cv: f64) -> f64 {
        -self.hills.iter().map(|h| h.gradient(cv)).sum::<f64>()
    }

    /// Attempt to deposit a Gaussian hill.
    ///
    /// Computes the adaptive width from the CV velocity (|cv − prev_cv| / dt_steps).
    pub fn maybe_deposit(&mut self, cv: f64) {
        let cv_vel = if let Some(prev) = self.prev_cv {
            (cv - prev).abs()
        } else {
            0.0
        };

        if self.deposition_stride > 0 && self.step.is_multiple_of(self.deposition_stride) {
            // Adaptive width
            let sigma = if cv_vel > 0.0 {
                self.width_min
                    + (self.width_max - self.width_min) * (-cv_vel / self.vel_scale).exp()
            } else {
                self.width_max
            };

            let h = if self.delta_t > 0.0 {
                let v_bias = self.compute_bias(cv);
                self.height * (-v_bias / (KB * (self.temperature + self.delta_t))).exp()
            } else {
                self.height
            };

            self.hills.push(GaussBias {
                center: cv,
                height: h,
                width: sigma,
            });
        }
        self.prev_cv = Some(cv);
        self.step += 1;
    }

    /// Number of hills deposited.
    pub fn n_hills(&self) -> usize {
        self.hills.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Well-Tempered ABF (Adaptive Biasing Force)
// ─────────────────────────────────────────────────────────────────────────────

/// Adaptive biasing force (ABF) method on a 1D CV.
///
/// Accumulates the mean force ⟨F·∂s/∂r⟩ in bins and subtracts it as
/// an adaptive bias after reaching a minimum sample threshold.
#[derive(Debug, Clone)]
pub struct AdaptiveBiasingForce {
    /// CV bin centres.
    pub bins: Vec<f64>,
    /// Accumulated force sum per bin.
    force_sum: Vec<f64>,
    /// Sample count per bin.
    counts: Vec<usize>,
    /// Minimum sample count before applying the bias.
    pub min_samples: usize,
    /// Well-tempered scale factor (0 = full ABF, 1.0 = no bias).
    pub well_tempered_scale: f64,
}

impl AdaptiveBiasingForce {
    /// Create a new ABF on a uniform grid from `cv_min` to `cv_max` with `n_bins` bins.
    pub fn new(
        cv_min: f64,
        cv_max: f64,
        n_bins: usize,
        min_samples: usize,
        well_tempered_scale: f64,
    ) -> Self {
        let bw = (cv_max - cv_min) / n_bins as f64;
        let bins: Vec<f64> = (0..n_bins)
            .map(|i| cv_min + (i as f64 + 0.5) * bw)
            .collect();
        Self {
            force_sum: vec![0.0; n_bins],
            counts: vec![0; n_bins],
            bins,
            min_samples,
            well_tempered_scale,
        }
    }

    fn bin_index(&self, cv: f64) -> Option<usize> {
        let n = self.bins.len();
        if n == 0 {
            return None;
        }
        let bw = if n > 1 {
            self.bins[1] - self.bins[0]
        } else {
            1.0
        };
        let cv_min = self.bins[0] - 0.5 * bw;
        let idx = ((cv - cv_min) / bw) as isize;
        if idx >= 0 && (idx as usize) < n {
            Some(idx as usize)
        } else {
            None
        }
    }

    /// Record an instantaneous force sample at the given CV value.
    pub fn accumulate(&mut self, cv: f64, force: f64) {
        if let Some(b) = self.bin_index(cv) {
            self.force_sum[b] += force;
            self.counts[b] += 1;
        }
    }

    /// Return the adaptive biasing force at the given CV position.
    ///
    /// Returns 0 if the bin has fewer than `min_samples` samples.
    pub fn bias_force(&self, cv: f64) -> f64 {
        if let Some(b) = self.bin_index(cv)
            && self.counts[b] >= self.min_samples
        {
            let mean_f = self.force_sum[b] / self.counts[b] as f64;
            return -mean_f * (1.0 - self.well_tempered_scale);
        }
        0.0
    }

    /// Mean force estimate at bin `b`.
    pub fn mean_force_at_bin(&self, b: usize) -> f64 {
        if self.counts[b] == 0 {
            0.0
        } else {
            self.force_sum[b] / self.counts[b] as f64
        }
    }

    /// Integrate the mean force to obtain a free energy profile (PMF).
    ///
    /// Uses the trapezoidal rule with 0 at the first bin.
    pub fn free_energy_profile(&self) -> Vec<f64> {
        let n = self.bins.len();
        if n < 2 {
            return vec![0.0; n];
        }
        let bw = self.bins[1] - self.bins[0];
        let mut fe = vec![0.0f64; n];
        for i in 1..n {
            let mf = 0.5 * (self.mean_force_at_bin(i - 1) + self.mean_force_at_bin(i));
            fe[i] = fe[i - 1] - mf * bw;
        }
        fe
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OPES (On-the-fly Probability Enhanced Sampling) — basic stub
// ─────────────────────────────────────────────────────────────────────────────

/// Basic OPES (On-the-fly Probability Enhanced Sampling) method.
///
/// Reference: Invernizzi & Parrinello, J. Phys. Chem. Lett. 11, 2731 (2020).
///
/// OPES estimates the probability density of the CV on-the-fly and
/// constructs a bias potential:
///   V(s) = (kT) * biasfactor * ln\[p(s) / p_0\]
///
/// where p(s) is the running KDE estimate of the biased distribution and
/// biasfactor = (T + ΔT) / T.
///
/// This implementation uses a simple kernel density estimate (KDE) with a
/// Gaussian kernel.
#[derive(Debug, Clone)]
pub struct Opes {
    /// Gaussian kernel bandwidth for KDE (CV units).
    pub bandwidth: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Bias factor = (T + ΔT) / T. Should be > 1.
    pub biasfactor: f64,
    /// Collected CV samples for KDE.
    samples: Vec<f64>,
    /// Deposition stride.
    pub deposition_stride: usize,
    /// Current step.
    step: usize,
}

impl Opes {
    /// Create a new OPES controller.
    pub fn new(
        bandwidth: f64,
        temperature: f64,
        biasfactor: f64,
        deposition_stride: usize,
    ) -> Self {
        Self {
            bandwidth,
            temperature,
            biasfactor,
            samples: Vec::new(),
            deposition_stride,
            step: 0,
        }
    }

    /// Record the current CV value (adds it to the kernel density estimate).
    pub fn observe(&mut self, cv: f64) {
        if self.deposition_stride > 0 && self.step.is_multiple_of(self.deposition_stride) {
            self.samples.push(cv);
        }
        self.step += 1;
    }

    /// Estimate the KDE probability density at `cv`.
    pub fn kde_density(&self, cv: f64) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 1e-30;
        }
        let bw = self.bandwidth;
        let norm = 1.0 / (bw * (2.0 * std::f64::consts::PI).sqrt() * n as f64);
        let sum: f64 = self
            .samples
            .iter()
            .map(|&s| {
                let z = (cv - s) / bw;
                (-0.5 * z * z).exp()
            })
            .sum();
        (norm * sum).max(1e-30)
    }

    /// Compute the OPES bias potential at `cv`.
    ///
    /// V(s) = kT * biasfactor * ln(p(s) / p_target)
    /// where p_target is estimated from the minimum observed density.
    pub fn bias_potential(&self, cv: f64) -> f64 {
        let kbt = KB * self.temperature;
        let p = self.kde_density(cv);
        // Use minimum density as reference (ensures non-negative bias in explored regions)
        let p_min = if self.samples.is_empty() {
            1e-30
        } else {
            self.samples
                .iter()
                .map(|&s| self.kde_density(s))
                .fold(f64::INFINITY, f64::min)
        };
        kbt * self.biasfactor * (p / p_min.max(1e-30)).ln()
    }

    /// Number of KDE samples collected.
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Funnel Metadynamics
// ─────────────────────────────────────────────────────────────────────────────

/// Funnel metadynamics restraint.
///
/// In funnel metadynamics a funnel-shaped restraint is applied to the CV
/// to confine the simulation in the binding pocket region and accelerate
/// ligand binding/unbinding sampling.
///
/// The funnel consists of:
/// - A **cylindrical** region (radius `r_cyl`) for small CV values (s < `s_cone`).
/// - A **conical** region that opens linearly for larger CV values.
///
/// The restraint energy (kJ/mol) is:
///   V_fun = k_fun * max(0, r - r_fun(s))²
///
/// where `r` is the perpendicular distance from the funnel axis and
/// `r_fun(s)` is the funnel radius at axial coordinate `s`.
#[derive(Debug, Clone)]
pub struct FunnelRestraint {
    /// Force constant for the funnel wall (kJ/mol/Å²).
    pub k_fun: f64,
    /// Radius of the cylindrical stem (Å).
    pub r_cyl: f64,
    /// CV value at which the funnel transitions from cylinder to cone.
    pub s_cone: f64,
    /// Half-angle of the cone (radians).
    pub half_angle: f64,
}

impl FunnelRestraint {
    /// Create a new funnel restraint.
    pub fn new(k_fun: f64, r_cyl: f64, s_cone: f64, half_angle_deg: f64) -> Self {
        Self {
            k_fun,
            r_cyl,
            s_cone,
            half_angle: half_angle_deg.to_radians(),
        }
    }

    /// Funnel radius at axial coordinate `s`.
    pub fn funnel_radius(&self, s: f64) -> f64 {
        if s <= self.s_cone {
            self.r_cyl
        } else {
            self.r_cyl + (s - self.s_cone) * self.half_angle.tan()
        }
    }

    /// Restraint energy at (s, r_perp).
    pub fn energy(&self, s: f64, r_perp: f64) -> f64 {
        let r_fun = self.funnel_radius(s);
        let excess = r_perp - r_fun;
        if excess > 0.0 {
            self.k_fun * excess * excess
        } else {
            0.0
        }
    }

    /// Radial restraint force at (s, r_perp).
    ///
    /// Positive force pushes inward (restoring).
    pub fn radial_force(&self, s: f64, r_perp: f64) -> f64 {
        let r_fun = self.funnel_radius(s);
        let excess = r_perp - r_fun;
        if excess > 0.0 {
            -2.0 * self.k_fun * excess
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Conformational Sampling via Simulated Annealing
// ─────────────────────────────────────────────────────────────────────────────

/// Simple simulated annealing schedule for conformational sampling.
///
/// Provides a monotonically decreasing temperature from `t_start` to `t_end`
/// over `n_steps` steps using either a linear or geometric cooling schedule.
#[derive(Debug, Clone)]
pub struct SimulatedAnnealingSchedule {
    /// Start temperature (K).
    pub t_start: f64,
    /// End temperature (K).
    pub t_end: f64,
    /// Total number of annealing steps.
    pub n_steps: usize,
    /// Cooling mode.
    pub mode: AnnealingMode,
}

/// Cooling schedule mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnealingMode {
    /// T decreases linearly.
    Linear,
    /// T decreases geometrically (exponential cooling).
    Geometric,
}

impl SimulatedAnnealingSchedule {
    /// Create a new annealing schedule.
    pub fn new(t_start: f64, t_end: f64, n_steps: usize, mode: AnnealingMode) -> Self {
        Self {
            t_start,
            t_end,
            n_steps,
            mode,
        }
    }

    /// Temperature at step `i`.
    pub fn temperature_at(&self, i: usize) -> f64 {
        if self.n_steps == 0 || i >= self.n_steps {
            return self.t_end;
        }
        let frac = i as f64 / (self.n_steps - 1).max(1) as f64;
        match self.mode {
            AnnealingMode::Linear => self.t_start + frac * (self.t_end - self.t_start),
            AnnealingMode::Geometric => self.t_start * (self.t_end / self.t_start).powf(frac),
        }
    }

    /// Metropolis acceptance probability for a move with energy change `delta_e`.
    pub fn accept_probability(&self, delta_e: f64, step: usize) -> f64 {
        let t = self.temperature_at(step);
        if delta_e <= 0.0 {
            1.0
        } else {
            (-delta_e / (KB * t)).exp()
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ReplicaExchange ───────────────────────────────────────────────────────

    #[test]
    fn test_replica_exchange_new() {
        let re = ReplicaExchange::new(vec![300.0, 400.0, 500.0]);
        assert_eq!(re.n_replicas, 3);
        assert_eq!(re.energies, vec![0.0, 0.0, 0.0]);
        assert_eq!(re.n_attempts, 0);
        assert_eq!(re.n_accepted, 0);
    }

    #[test]
    fn test_replica_exchange_geometric_ladder() {
        let re = ReplicaExchange::geometric_ladder(300.0, 600.0, 4);
        assert_eq!(re.n_replicas, 4);
        assert!((re.temperatures[0] - 300.0).abs() < 1e-6);
        assert!((re.temperatures[3] - 600.0).abs() < 1e-6);
        // Geometric → constant ratio
        let r0 = re.temperatures[1] / re.temperatures[0];
        let r1 = re.temperatures[2] / re.temperatures[1];
        assert!((r0 - r1).abs() < 1e-8, "Geometric ratio must be constant");
    }

    #[test]
    fn test_acceptance_probability_zero_delta() {
        // When energies are equal or beta difference × energy difference > 0 → p = 1
        let p = ReplicaExchange::acceptance_probability(
            -100.0,
            -100.0,
            1.0 / (KB * 300.0),
            1.0 / (KB * 400.0),
        );
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_acceptance_probability_negative_delta() {
        let bi = 1.0 / (KB * 300.0);
        let bj = 1.0 / (KB * 400.0);
        // ei > ej and bi > bj → delta < 0 → p < 1
        let p = ReplicaExchange::acceptance_probability(-50.0, -200.0, bj, bi);
        // (bj - bi) < 0, (ei - ej) = -50 - (-200) = 150 > 0 → delta < 0 → p < 1
        assert!(
            p > 0.0 && p <= 1.0,
            "Probability must be in (0, 1], got {p}"
        );
    }

    #[test]
    fn test_attempt_swap_high_energy_replica_swaps_deterministically() {
        // For P_acc = 1 we need delta = (bi - bj)(Ei - Ej) >= 0.
        // With temps = [300, 600]: beta0 = 1/(kB*300) > beta1 = 1/(kB*600)
        //   so bi - bj > 0.
        // For delta >= 0 we need Ei - Ej >= 0, i.e. E0 >= E1.
        // Set E0 = -50 (higher energy), E1 = -200 (lower energy):
        //   delta = (positive)(positive) = positive → P_acc = 1 → deterministic accept.
        let mut re = ReplicaExchange::new(vec![300.0, 600.0]);
        let betas = re.betas();
        let energies = vec![-50.0, -200.0]; // E0 > E1
        let swapped = re.attempt_swap(&energies, &betas);
        assert_eq!(swapped.len(), 1, "Should accept exactly one swap");
        assert_eq!(swapped[0], (0, 1));
    }

    #[test]
    fn test_attempt_swap_updates_step() {
        let mut re = ReplicaExchange::new(vec![300.0, 400.0]);
        let betas = re.betas();
        let energies = vec![0.0, 0.0];
        re.attempt_swap(&energies, &betas);
        assert_eq!(re.step, 1);
    }

    #[test]
    fn test_attempt_swap_random_rejection() {
        let mut re = ReplicaExchange::new(vec![300.0, 400.0]);
        let betas = vec![1.0 / (KB * 300.0), 1.0 / (KB * 400.0)];
        // P_acc < 1 scenario: bi > bj (300 K has higher beta), E0 < E1
        // delta = (bi - bj)(E0 - E1) = (positive)(negative) < 0 → p_acc < 1
        // Use E0 very negative, E1 near zero → E0 - E1 very negative → p_acc ≈ 0
        let energies = vec![-1000.0, 0.0];
        let random_vals = vec![0.9999]; // near 1 → rejects unless p_acc > 0.9999
        let swapped = re.attempt_swap_random(&energies, &betas, &random_vals);
        // p_acc = exp((bi - bj)(-1000)) ≈ exp(-very_large) ≈ 0 → 0.9999 > p_acc → reject
        assert_eq!(
            swapped.len(),
            0,
            "Should reject swap when random value > P_acc"
        );
    }

    #[test]
    fn test_acceptance_rate_tracking() {
        let mut re = ReplicaExchange::new(vec![300.0, 600.0]);
        // Use attempt_swap_random with random=0.0 to always accept
        // (p_acc >= 0.0 is always true)
        for _ in 0..5 {
            let betas = re.betas(); // recalculate after each swap (temps change)
            let energies = vec![-50.0, -200.0]; // E0 > E1 with beta0 > beta1 → delta > 0 → p_acc=1
            re.attempt_swap(&energies, &betas);
        }
        let rate = re.acceptance_rate();
        assert!(
            rate > 0.0 && rate <= 1.0,
            "Acceptance rate must be in [0,1], got {rate}"
        );
        assert_eq!(re.n_attempts, 5);
    }

    #[test]
    fn test_betas_are_inverse_kt() {
        let temps = vec![300.0, 400.0, 500.0];
        let re = ReplicaExchange::new(temps.clone());
        let betas = re.betas();
        for (i, (&t, &b)) in temps.iter().zip(betas.iter()).enumerate() {
            let expected = 1.0 / (KB * t);
            assert!((b - expected).abs() < 1e-10, "beta[{i}] mismatch");
        }
    }

    // ── Metadynamics ──────────────────────────────────────────────────────────

    #[test]
    fn test_metadynamics_bias_zero_before_first_hill() {
        let meta = Metadynamics::new(1.0, 0.5, 10, 300.0, 0.0);
        let bias = meta.compute_bias(0.0);
        assert!(
            bias.abs() < 1e-12,
            "Bias must be zero before any hills, got {bias}"
        );
    }

    #[test]
    fn test_metadynamics_add_gaussian_increases_bias() {
        let mut meta = Metadynamics::new(1.0, 0.5, 10, 300.0, 0.0);
        let b0 = meta.compute_bias(0.0);
        meta.add_gaussian(0.0, 1.0, 0.5);
        let b1 = meta.compute_bias(0.0);
        assert!(b1 > b0, "Bias must increase after adding a Gaussian hill");
    }

    #[test]
    fn test_metadynamics_add_gaussian_peaks_at_center() {
        let mut meta = Metadynamics::new(1.0, 0.5, 1, 300.0, 0.0);
        meta.add_gaussian(2.0, 1.0, 0.5);
        let at_center = meta.compute_bias(2.0);
        let off_center = meta.compute_bias(3.0);
        assert!(at_center > off_center, "Bias must peak at hill center");
        assert!(
            (at_center - 1.0).abs() < 1e-10,
            "Bias at center must equal hill height"
        );
    }

    #[test]
    fn test_metadynamics_bias_force_zero_at_hill_center() {
        let mut meta = Metadynamics::new(1.0, 0.5, 1, 300.0, 0.0);
        meta.add_gaussian(0.0, 1.0, 0.5);
        let force = meta.compute_bias_force(0.0);
        assert!(
            force.abs() < 1e-10,
            "Bias force at hill center must be zero, got {force}"
        );
    }

    #[test]
    fn test_metadynamics_maybe_deposit_at_stride() {
        let mut meta = Metadynamics::new(1.0, 0.5, 3, 300.0, 0.0);
        // step=0 → deposit
        meta.maybe_deposit(0.0);
        assert_eq!(meta.n_hills(), 1, "Should deposit at step 0");
        // step=1 → no deposit
        meta.maybe_deposit(0.0);
        assert_eq!(meta.n_hills(), 1, "Should not deposit at step 1");
        // step=2 → no deposit
        meta.maybe_deposit(0.0);
        assert_eq!(meta.n_hills(), 1, "Should not deposit at step 2");
        // step=3 → deposit
        meta.maybe_deposit(0.0);
        assert_eq!(meta.n_hills(), 2, "Should deposit at step 3");
    }

    #[test]
    fn test_metadynamics_well_tempered_height_decreases() {
        let mut meta = Metadynamics::new(1.0, 0.5, 1, 300.0, 300.0);
        // deposit 5 hills at the same point; heights should decrease
        for _ in 0..5 {
            meta.maybe_deposit(0.0);
        }
        assert!(meta.n_hills() == 5, "Should have 5 hills");
        // Check compute_bias grows but slows (first hill height > later hills)
        // We can't access internal hills directly, but bias grows less per step
        let b_before = meta.compute_bias(0.0);
        meta.maybe_deposit(0.0);
        let b_after = meta.compute_bias(0.0);
        // The increment should be < 1.0 (initial height)
        assert!(
            b_after - b_before < 1.0,
            "Well-tempered: increment should be < initial height"
        );
    }

    #[test]
    fn test_metadynamics_free_energy_estimate_is_negated_bias() {
        let mut meta = Metadynamics::new(2.0, 0.5, 1, 300.0, 0.0);
        meta.maybe_deposit(1.0);
        let grid = vec![0.0, 1.0, 2.0];
        let fe = meta.free_energy_estimate(&grid);
        for &cv in &grid {
            let bias = meta.compute_bias(cv);
            let fe_at = fe[grid.iter().position(|&x| (x - cv).abs() < 1e-12).unwrap()];
            assert!((fe_at + bias).abs() < 1e-10, "FE at cv={cv} must be -bias");
        }
    }

    // ── UmbrellaSamplingWindow ────────────────────────────────────────────────

    #[test]
    fn test_umbrella_window_bias_zero_at_center() {
        let w = UmbrellaSamplingWindow::new(2.5, 100.0);
        let v = w.bias_potential(2.5);
        assert!(v.abs() < 1e-12, "Bias at center must be zero, got {v}");
    }

    #[test]
    fn test_umbrella_window_force_zero_at_center() {
        let w = UmbrellaSamplingWindow::new(2.5, 100.0);
        let f = w.bias_force(2.5);
        assert!(f.abs() < 1e-12, "Force at center must be zero, got {f}");
    }

    #[test]
    fn test_umbrella_window_force_restoring() {
        let w = UmbrellaSamplingWindow::new(1.0, 50.0);
        // xi > xi0 → force < 0 (pulls back)
        assert!(w.bias_force(2.0) < 0.0);
        // xi < xi0 → force > 0 (pushes back)
        assert!(w.bias_force(0.0) > 0.0);
    }

    #[test]
    fn test_umbrella_window_statistics() {
        let mut w = UmbrellaSamplingWindow::new(3.0, 1.0);
        for x in [1.0, 2.0, 3.0, 4.0, 5.0] {
            w.add_sample(x);
        }
        let (mean, var) = w.statistics();
        assert!((mean - 3.0).abs() < 1e-10, "Mean should be 3, got {mean}");
        assert!((var - 2.0).abs() < 1e-10, "Variance should be 2, got {var}");
    }

    // ── UmbrellaSampling (collection) ─────────────────────────────────────────

    #[test]
    fn test_umbrella_sampling_from_centers() {
        let centers = vec![0.0, 1.0, 2.0, 3.0];
        let us = UmbrellaSampling::from_centers(&centers, 50.0);
        assert_eq!(us.n_windows(), 4);
        assert!((us.windows[2].center - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_umbrella_sampling_add_window() {
        let mut us = UmbrellaSampling::new();
        us.add_window(1.5, 100.0);
        us.add_window(2.5, 100.0);
        assert_eq!(us.n_windows(), 2);
    }

    #[test]
    fn test_umbrella_sampling_compute_pmf_returns_correct_n_bins() {
        let centers = vec![0.0, 1.0, 2.0];
        let mut us = UmbrellaSampling::from_centers(&centers, 50.0);
        // Populate windows with some samples
        for (i, w) in us.windows.iter_mut().enumerate() {
            for j in 0..20 {
                w.add_sample(i as f64 + j as f64 * 0.05 - 0.5);
            }
        }
        let (centers_out, pmf) = us.compute_pmf(20, 2.479); // kT ≈ 300K in kJ/mol
        assert_eq!(centers_out.len(), 20);
        assert_eq!(pmf.len(), 20);
    }

    #[test]
    fn test_umbrella_sampling_compute_pmf_minimum_at_highest_density() {
        // Single window at xi0=1.0 with many samples near center
        let mut us = UmbrellaSampling::new();
        us.add_window(1.0, 200.0);
        // Add samples tightly around center
        for j in 0..100 {
            us.windows[0].add_sample(1.0 + (j as f64 - 50.0) * 0.01);
        }
        let (bin_centers, pmf) = us.compute_pmf(20, 2.479);

        // Find the finite minimum
        let min_idx = pmf
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v.is_finite())
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .expect("Should find finite minimum in PMF");

        let min_pos = bin_centers[min_idx];
        assert!(
            (min_pos - 1.0).abs() < 0.3,
            "PMF minimum at {min_pos} should be near window center 1.0"
        );
    }

    #[test]
    fn test_replica_exchange_acceptance_probability_range() {
        // For any (ei, ej, beta_i, beta_j), probability must be in [0, 1]
        let bi = 1.0 / (KB * 300.0);
        let bj = 1.0 / (KB * 400.0);
        for &ei in &[-200.0, -100.0, 0.0, 100.0] {
            for &ej in &[-200.0, -100.0, 0.0, 100.0] {
                let p = ReplicaExchange::acceptance_probability(ei, ej, bi, bj);
                assert!(
                    (0.0..=1.0 + 1e-12).contains(&p),
                    "P_acc={p} must be in [0,1] for ei={ei}, ej={ej}"
                );
            }
        }
    }

    // ── AdaptiveBiasMetadynamics ───────────────────────────────────────────

    #[test]
    fn test_adaptive_meta_no_bias_before_first_deposit() {
        let meta = AdaptiveBiasMetadynamics::new(1.0, 0.2, 0.8, 0.5, 10, 300.0, 0.0);
        let bias = meta.compute_bias(0.0);
        assert!(
            bias.abs() < 1e-12,
            "Adaptive meta bias should be 0 initially"
        );
    }

    #[test]
    fn test_adaptive_meta_deposits_at_stride() {
        let mut meta = AdaptiveBiasMetadynamics::new(1.0, 0.2, 0.8, 0.5, 3, 300.0, 0.0);
        meta.maybe_deposit(0.0); // step 0 → deposit
        assert_eq!(meta.n_hills(), 1);
        meta.maybe_deposit(0.1); // step 1 → no
        meta.maybe_deposit(0.2); // step 2 → no
        meta.maybe_deposit(0.3); // step 3 → deposit
        assert_eq!(meta.n_hills(), 2);
    }

    #[test]
    fn test_adaptive_meta_bias_increases() {
        let mut meta = AdaptiveBiasMetadynamics::new(1.0, 0.3, 0.6, 0.5, 1, 300.0, 0.0);
        meta.maybe_deposit(0.0);
        let b0 = meta.compute_bias(0.0);
        meta.maybe_deposit(0.0);
        let b1 = meta.compute_bias(0.0);
        assert!(b1 > b0, "Bias should increase after more deposits");
    }

    #[test]
    fn test_adaptive_meta_bias_force_finite() {
        let mut meta = AdaptiveBiasMetadynamics::new(1.0, 0.2, 0.8, 0.5, 1, 300.0, 0.0);
        meta.maybe_deposit(0.0);
        let f = meta.compute_bias_force(0.5);
        assert!(f.is_finite(), "Bias force should be finite");
    }

    // ── AdaptiveBiasingForce ───────────────────────────────────────────────

    #[test]
    fn test_abf_zero_before_min_samples() {
        let abf = AdaptiveBiasingForce::new(-1.0, 1.0, 10, 5, 0.0);
        // No samples accumulated → bias force should be 0
        let f = abf.bias_force(0.0);
        assert!(
            f.abs() < 1e-12,
            "ABF bias should be 0 before min_samples, got {f}"
        );
    }

    #[test]
    fn test_abf_accumulates_and_computes_mean_force() {
        let mut abf = AdaptiveBiasingForce::new(0.0, 1.0, 5, 3, 0.0);
        // Add 4 force samples of 10.0 at cv=0.1 (bin 0)
        for _ in 0..4 {
            abf.accumulate(0.1, 10.0);
        }
        let mf = abf.mean_force_at_bin(0);
        assert!(
            (mf - 10.0).abs() < 1e-10,
            "Mean force should be 10, got {mf}"
        );
        let f = abf.bias_force(0.1);
        assert!(
            (f - (-10.0)).abs() < 1e-10,
            "ABF bias force should negate mean force, got {f}"
        );
    }

    #[test]
    fn test_abf_free_energy_profile_length() {
        let abf = AdaptiveBiasingForce::new(0.0, 2.0, 10, 1, 0.0);
        let fe = abf.free_energy_profile();
        assert_eq!(fe.len(), 10, "FE profile should have n_bins entries");
    }

    #[test]
    fn test_abf_free_energy_starts_at_zero() {
        let abf = AdaptiveBiasingForce::new(0.0, 1.0, 5, 1, 0.0);
        let fe = abf.free_energy_profile();
        assert!(fe[0].abs() < 1e-12, "FE profile should start at 0");
    }

    // ── OPES ──────────────────────────────────────────────────────────────

    #[test]
    fn test_opes_zero_density_before_samples() {
        let opes = Opes::new(0.2, 300.0, 5.0, 1);
        let p = opes.kde_density(0.0);
        assert!(
            p > 0.0,
            "KDE density should be > 0 even before samples (floor)"
        );
    }

    #[test]
    fn test_opes_density_peaks_at_observed_point() {
        let mut opes = Opes::new(0.2, 300.0, 5.0, 1);
        opes.observe(1.0);
        let p_at = opes.kde_density(1.0);
        let p_far = opes.kde_density(5.0);
        assert!(p_at > p_far, "KDE density should peak at observed CV");
    }

    #[test]
    fn test_opes_n_samples_increases_at_stride() {
        let mut opes = Opes::new(0.2, 300.0, 5.0, 2);
        opes.observe(0.0); // step 0 → sample
        opes.observe(1.0); // step 1 → no
        opes.observe(2.0); // step 2 → sample
        assert_eq!(opes.n_samples(), 2);
    }

    #[test]
    fn test_opes_bias_potential_finite() {
        let mut opes = Opes::new(0.3, 300.0, 5.0, 1);
        for i in 0..20 {
            opes.observe(i as f64 * 0.1);
        }
        let v = opes.bias_potential(0.5);
        assert!(v.is_finite(), "OPES bias potential should be finite");
    }

    // ── FunnelRestraint ───────────────────────────────────────────────────

    #[test]
    fn test_funnel_radius_cylindrical() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        // s < s_cone → cylindrical radius
        assert!((fun.funnel_radius(3.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_funnel_radius_increases_beyond_cone() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        let r5 = fun.funnel_radius(5.0);
        let r10 = fun.funnel_radius(10.0);
        assert!(r10 > r5, "Funnel radius should increase beyond s_cone");
    }

    #[test]
    fn test_funnel_energy_zero_inside() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        // r_perp = 1.0 < r_cyl = 2.0 → energy = 0
        let e = fun.energy(3.0, 1.0);
        assert!(e.abs() < 1e-12, "Energy inside funnel should be 0, got {e}");
    }

    #[test]
    fn test_funnel_energy_positive_outside() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        let e = fun.energy(3.0, 3.0); // r_perp=3 > r_cyl=2 → penalty
        assert!(e > 0.0, "Energy outside funnel should be positive, got {e}");
    }

    #[test]
    fn test_funnel_radial_force_zero_inside() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        let f = fun.radial_force(3.0, 1.0);
        assert!(f.abs() < 1e-12, "Force inside funnel should be 0");
    }

    #[test]
    fn test_funnel_radial_force_restoring() {
        let fun = FunnelRestraint::new(100.0, 2.0, 5.0, 30.0);
        let f = fun.radial_force(3.0, 3.0);
        assert!(
            f < 0.0,
            "Restoring radial force should be negative (pushes inward)"
        );
    }

    // ── SimulatedAnnealingSchedule ─────────────────────────────────────────

    #[test]
    fn test_annealing_linear_endpoints() {
        let sched = SimulatedAnnealingSchedule::new(1000.0, 100.0, 100, AnnealingMode::Linear);
        assert!((sched.temperature_at(0) - 1000.0).abs() < 1e-6);
        assert!((sched.temperature_at(99) - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_annealing_linear_monotone_decreasing() {
        let sched = SimulatedAnnealingSchedule::new(1000.0, 100.0, 10, AnnealingMode::Linear);
        let temps: Vec<f64> = (0..10).map(|i| sched.temperature_at(i)).collect();
        for w in temps.windows(2) {
            assert!(w[1] <= w[0], "Temperature should decrease monotonically");
        }
    }

    #[test]
    fn test_annealing_geometric_endpoints() {
        let sched = SimulatedAnnealingSchedule::new(1000.0, 10.0, 100, AnnealingMode::Geometric);
        assert!((sched.temperature_at(0) - 1000.0).abs() < 1e-6);
        assert!((sched.temperature_at(99) - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_annealing_accept_probability_downhill() {
        let sched = SimulatedAnnealingSchedule::new(1000.0, 100.0, 10, AnnealingMode::Linear);
        let p = sched.accept_probability(-5.0, 0);
        assert!(
            (p - 1.0).abs() < 1e-12,
            "Downhill move should always be accepted"
        );
    }

    #[test]
    fn test_annealing_accept_probability_uphill_decreases_with_step() {
        let sched = SimulatedAnnealingSchedule::new(1000.0, 10.0, 100, AnnealingMode::Geometric);
        let p_hot = sched.accept_probability(10.0, 0);
        let p_cold = sched.accept_probability(10.0, 90);
        assert!(
            p_hot > p_cold,
            "Hot acceptance should be higher than cold for uphill move"
        );
    }
}
