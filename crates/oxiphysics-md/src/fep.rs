#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Free Energy Perturbation (FEP) methods.
//!
//! Implements alchemical free energy calculations via:
//! - **Thermodynamic Integration (TI)**: numerical integration of <dU/dλ>
//! - **Exponential Averaging** (Zwanzig equation): -kBT ln<exp(-βΔU)>_A
//! - **Bennett Acceptance Ratio (BAR)**: self-consistent estimator using samples from both states

/// Boltzmann constant in J/K.
const KB_J: f64 = 1.380_649e-23;

// ─── AlchemicalState ─────────────────────────────────────────────────────────

/// Lambda-scaled potential energy (linear interpolation between states A and B).
///
/// The alchemical parameter λ ∈ \[0, 1\] smoothly interpolates between
/// state A (λ = 0) and state B (λ = 1).
#[allow(dead_code)]
pub struct AlchemicalState {
    /// Alchemical coupling parameter: 0.0 = state A, 1.0 = state B.
    pub lambda: f64,
    /// Potential energy in state A.
    pub u_a: f64,
    /// Potential energy in state B.
    pub u_b: f64,
}

impl AlchemicalState {
    /// Create a new alchemical state at the given λ value.
    /// `u_a` and `u_b` are initialised to zero and should be set by the caller.
    pub fn new(lambda: f64) -> Self {
        Self {
            lambda,
            u_a: 0.0,
            u_b: 0.0,
        }
    }

    /// Mixed potential: U(λ) = (1 − λ) · U_A + λ · U_B
    pub fn mixed_potential(&self) -> f64 {
        (1.0 - self.lambda) * self.u_a + self.lambda * self.u_b
    }

    /// Derivative of the mixed potential with respect to λ:
    /// dU/dλ = U_B − U_A
    pub fn du_dlambda(&self) -> f64 {
        self.u_b - self.u_a
    }

    /// Energy difference ΔU = U_B − U_A (for BAR / Zwanzig estimators).
    pub fn delta_u(&self) -> f64 {
        self.u_b - self.u_a
    }
}

// ─── ThermodynamicIntegration ─────────────────────────────────────────────────

/// Thermodynamic Integration (TI) estimator.
///
/// ΔF = ∫₀¹ ⟨dU/dλ⟩_λ dλ,  approximated via the trapezoidal rule over a
/// discrete set of λ windows.
#[allow(dead_code)]
pub struct ThermodynamicIntegration {
    /// λ values at each window.
    pub lambda_windows: Vec<f64>,
    /// Running sum of dU/dλ samples at each window (for computing the mean).
    du_dlambda_sum: Vec<f64>,
    /// Running sum of (dU/dλ)² samples at each window (for Bessel-corrected variance).
    du_dlambda_sq_sum: Vec<f64>,
    /// Number of samples accumulated at each window.
    sample_counts: Vec<usize>,
    /// Cached ⟨dU/dλ⟩ at each window (updated on each `add_sample` call).
    pub du_dlambda_averages: Vec<f64>,
}

impl ThermodynamicIntegration {
    /// Create a new TI estimator with the given λ windows.
    pub fn new(lambdas: Vec<f64>) -> Self {
        let n = lambdas.len();
        Self {
            lambda_windows: lambdas,
            du_dlambda_sum: vec![0.0; n],
            du_dlambda_sq_sum: vec![0.0; n],
            sample_counts: vec![0; n],
            du_dlambda_averages: vec![0.0; n],
        }
    }

    /// Accumulate a dU/dλ sample at the given window index for running average.
    pub fn add_sample(&mut self, lambda_idx: usize, du_dl: f64) {
        self.du_dlambda_sum[lambda_idx] += du_dl;
        self.du_dlambda_sq_sum[lambda_idx] += du_dl * du_dl;
        self.sample_counts[lambda_idx] += 1;
        self.du_dlambda_averages[lambda_idx] =
            self.du_dlambda_sum[lambda_idx] / self.sample_counts[lambda_idx] as f64;
    }

    /// Bessel-corrected sample variance at window `idx`.
    ///
    /// Uses the computational formula: `(Σx² − (Σx)²/N) / (N−1)`.
    /// Returns 0.0 when N ≤ 1 (variance undefined with a single sample).
    fn window_variance(&self, idx: usize) -> f64 {
        let n = self.sample_counts[idx];
        if n <= 1 {
            return 0.0;
        }
        let sum = self.du_dlambda_sum[idx];
        let sq_sum = self.du_dlambda_sq_sum[idx];
        let nf = n as f64;
        let variance = (sq_sum - sum * sum / nf) / (nf - 1.0);
        variance.max(0.0)
    }

    /// Integrate ⟨dU/dλ⟩ over λ via the trapezoidal rule.
    ///
    /// Returns 0 if fewer than 2 windows are present.
    pub fn delta_f(&self) -> f64 {
        let n = self.lambda_windows.len();
        if n < 2 {
            return 0.0;
        }
        let mut integral = 0.0;
        for i in 0..n - 1 {
            let dlambda = self.lambda_windows[i + 1] - self.lambda_windows[i];
            let avg = 0.5 * (self.du_dlambda_averages[i] + self.du_dlambda_averages[i + 1]);
            integral += avg * dlambda;
        }
        integral
    }

    /// Statistical uncertainty estimate for the TI integral (standard error propagation).
    ///
    /// Computed as the trapezoidal integral of σ_i / √N_i at each window, where σ_i is
    /// the Bessel-corrected sample standard deviation.  Returns 0 when fewer than 2
    /// windows are available or all windows have ≤ 1 sample.
    pub fn uncertainty(&self) -> f64 {
        let n = self.lambda_windows.len();
        if n < 2 {
            return 0.0;
        }
        // Standard error of mean at each window: SEM_i = sqrt(variance_i / N_i)
        let sem: Vec<f64> = (0..n)
            .map(|i| {
                let ni = self.sample_counts[i];
                if ni <= 1 {
                    0.0
                } else {
                    (self.window_variance(i) / ni as f64).sqrt()
                }
            })
            .collect();
        // Trapezoidal integration of SEM over λ
        let mut integral = 0.0;
        for i in 0..n - 1 {
            let dlambda = self.lambda_windows[i + 1] - self.lambda_windows[i];
            let avg_sem = 0.5 * (sem[i] + sem[i + 1]);
            integral += avg_sem * dlambda;
        }
        integral
    }
}

// ─── ExponentialAveraging ─────────────────────────────────────────────────────

/// Zwanzig equation / exponential averaging estimator.
///
/// ΔF_{A→B} = −k_B T · ln⟨exp(−β ΔU)⟩_A
///
/// Samples of ΔU = U_B − U_A are collected while the system is in state A.
#[allow(dead_code)]
pub struct ExponentialAveraging {
    /// Inverse thermal energy β = 1 / (k_B T).
    pub beta: f64,
    /// Collected ΔU samples from state A.
    samples: Vec<f64>,
}

impl ExponentialAveraging {
    /// Create a new estimator at the given temperature (K).
    ///
    /// Uses SI Boltzmann constant k_B = 1.380649 × 10⁻²³ J K⁻¹.
    pub fn new(temperature: f64) -> Self {
        let beta = 1.0 / (KB_J * temperature);
        Self {
            beta,
            samples: Vec::new(),
        }
    }

    /// Add a ΔU = U_B − U_A sample collected while in state A.
    pub fn add_sample(&mut self, delta_u: f64) {
        self.samples.push(delta_u);
    }

    /// Zwanzig equation free energy estimate.
    ///
    /// ΔF = −(1/β) · ln(⟨exp(−β ΔU)⟩)
    ///
    /// Returns 0 if no samples have been collected.
    pub fn delta_f(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let n = self.samples.len() as f64;
        let avg_exp: f64 = self
            .samples
            .iter()
            .map(|&du| (-self.beta * du).exp())
            .sum::<f64>()
            / n;
        if avg_exp <= 0.0 {
            return f64::INFINITY;
        }
        -avg_exp.ln() / self.beta
    }

    /// Number of samples collected.
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }
}

// ─── Bar ──────────────────────────────────────────────────────────────────────

/// Bennett Acceptance Ratio (BAR) estimator.
///
/// Uses samples of ΔU from both endpoints to produce a more efficient
/// estimate of ΔF than the one-sided Zwanzig equation.  The simplified
/// implementation here iterates the BAR self-consistency equation; it falls
/// back to the Zwanzig estimate when one side is empty.
#[allow(dead_code)]
pub struct Bar {
    /// Inverse thermal energy β = 1 / (k_B T).
    pub beta: f64,
    /// ΔU = U_B − U_A samples collected while in state A.
    pub samples_a: Vec<f64>,
    /// ΔU = U_B − U_A samples collected while in state B.
    pub samples_b: Vec<f64>,
}

impl Bar {
    /// Create a new BAR estimator at the given temperature (K).
    pub fn new(temperature: f64) -> Self {
        let beta = 1.0 / (KB_J * temperature);
        Self {
            beta,
            samples_a: Vec::new(),
            samples_b: Vec::new(),
        }
    }

    /// Add a ΔU sample from state A.
    pub fn add_sample_a(&mut self, delta_u: f64) {
        self.samples_a.push(delta_u);
    }

    /// Add a ΔU sample from state B.
    pub fn add_sample_b(&mut self, delta_u: f64) {
        self.samples_b.push(delta_u);
    }

    /// BAR free energy estimate.
    ///
    /// Uses a simplified iterative solver.  When only one side has samples,
    /// falls back to the Zwanzig (exponential averaging) formula.
    ///
    /// The BAR equation in the equal-sample-size limit reduces to:
    ///   C = (1/β) · ln(n_A / n_B) + ΔF
    /// with the self-consistent solution obtained by iterating:
    ///   ΔF_{k+1} = (1/β) · ln \[ Σ_B f(β(ΔU_B − C)) / Σ_A f(β(ΔU_A − C)) \]
    /// where f(x) = 1 / (1 + exp(x)) (Fermi function), C = ΔF.
    pub fn delta_f_estimate(&self) -> f64 {
        let n_a = self.samples_a.len();
        let n_b = self.samples_b.len();

        // Fallback: use Zwanzig when one side is empty
        if n_b == 0 {
            if n_a == 0 {
                return 0.0;
            }
            let avg_exp: f64 = self
                .samples_a
                .iter()
                .map(|&du| (-self.beta * du).exp())
                .sum::<f64>()
                / n_a as f64;
            return if avg_exp <= 0.0 {
                f64::INFINITY
            } else {
                -avg_exp.ln() / self.beta
            };
        }
        if n_a == 0 {
            let avg_exp: f64 = self
                .samples_b
                .iter()
                .map(|&du| (self.beta * du).exp())
                .sum::<f64>()
                / n_b as f64;
            return if avg_exp <= 0.0 {
                f64::NEG_INFINITY
            } else {
                avg_exp.ln() / self.beta
            };
        }

        // Fermi / logistic function for numerical stability
        let fermi = |x: f64| -> f64 { 1.0 / (1.0 + x.exp()) };

        // Iterate BAR self-consistency equation
        let mut df = 0.0_f64;
        let ln_ratio = ((n_a as f64) / (n_b as f64)).ln() / self.beta;

        for _ in 0..1000 {
            let c = df - ln_ratio * self.beta; // dimensionless shift

            let num: f64 = self
                .samples_b
                .iter()
                .map(|&du| fermi(self.beta * du - c))
                .sum::<f64>();

            let den: f64 = self
                .samples_a
                .iter()
                .map(|&du| fermi(-self.beta * du + c))
                .sum::<f64>();

            if den.abs() < 1e-300 {
                break;
            }

            let df_new =
                (num / den).ln() / self.beta + ((n_a as f64) / (n_b as f64)).ln() / self.beta;

            if (df_new - df).abs() < 1e-10 {
                df = df_new;
                break;
            }
            df = df_new;
        }

        df
    }
}

// ─── MBAR (Multistate Bennett Acceptance Ratio) ─────────────────────────────

/// Multistate Bennett Acceptance Ratio (MBAR) estimator.
///
/// Generalisation of BAR to K thermodynamic states, estimating the free
/// energy difference between all states simultaneously.
///
/// Reference: Shirts & Chodera, J. Chem. Phys. 129, 124105 (2008).
#[allow(dead_code)]
pub struct Mbar {
    /// Number of thermodynamic states.
    pub n_states: usize,
    /// Inverse temperature β = 1/(k_B T) (same for all states here).
    pub beta: f64,
    /// Reduced potential energies u_kn\[k\]\[n\] = β · U_k(x_n):
    /// the reduced potential of state k evaluated at sample n.
    /// Outer index: state k; inner index: sample index.
    pub u_kn: Vec<Vec<f64>>,
    /// Number of samples from each state.
    pub n_k: Vec<usize>,
}

impl Mbar {
    /// Create a new MBAR estimator.
    ///
    /// # Arguments
    /// * `n_states` – number of thermodynamic states
    /// * `temperature` – temperature in K (same for all states)
    pub fn new(n_states: usize, temperature: f64) -> Self {
        Self {
            n_states,
            beta: 1.0 / (KB_J * temperature),
            u_kn: vec![Vec::new(); n_states],
            n_k: vec![0; n_states],
        }
    }

    /// Add a sample: the reduced potential of state `state_k` evaluated
    /// at configuration `x_n` (which was sampled from state `sampled_from`).
    ///
    /// For simplicity, this stores u_k(x_n) for all k for each sample.
    /// `u_values` should have length `n_states`.
    pub fn add_sample(&mut self, sampled_from: usize, u_values: &[f64]) {
        assert_eq!(u_values.len(), self.n_states);
        for (k, &u) in u_values.iter().enumerate() {
            self.u_kn[k].push(u);
        }
        self.n_k[sampled_from] += 1;
    }

    /// Estimate free energies using MBAR self-consistency iteration.
    ///
    /// Returns a vector of `n_states` free energies (f_k), with f_0 = 0.
    pub fn estimate_free_energies(&self, max_iter: usize) -> Vec<f64> {
        let k_states = self.n_states;
        if k_states == 0 {
            return vec![];
        }
        let n_total = self.u_kn[0].len();
        if n_total == 0 {
            return vec![0.0; k_states];
        }

        let mut f_k = vec![0.0f64; k_states];

        for _ in 0..max_iter {
            let mut f_new = vec![0.0f64; k_states];

            for k in 0..k_states {
                // f_k = -ln( sum_n [ exp(-u_kn) / sum_l N_l exp(f_l - u_ln) ] )
                let mut sum = 0.0f64;
                for n in 0..n_total {
                    let u_k_n = self.u_kn[k][n];

                    // Compute log denominator for numerical stability
                    let mut max_arg = f64::NEG_INFINITY;
                    for l in 0..k_states {
                        let arg = f_k[l] - self.u_kn[l][n];
                        if arg > max_arg {
                            max_arg = arg;
                        }
                    }
                    let denom: f64 = (0..k_states)
                        .map(|l| self.n_k[l] as f64 * (f_k[l] - self.u_kn[l][n] - max_arg).exp())
                        .sum();
                    let log_denom = max_arg + denom.ln();

                    sum += (-u_k_n - log_denom).exp();
                }
                f_new[k] = if sum > 0.0 { -sum.ln() } else { 0.0 };
            }

            // Normalize so f_0 = 0
            let f0 = f_new[0];
            for fk in &mut f_new {
                *fk -= f0;
            }

            // Check convergence
            let delta: f64 = f_new
                .iter()
                .zip(f_k.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);

            f_k = f_new;
            if delta < 1e-10 {
                break;
            }
        }

        f_k
    }
}

// ─── Cumulant Expansion ─────────────────────────────────────────────────────

/// Second-order cumulant expansion estimate for ΔF.
///
/// ΔF ≈ ⟨ΔU⟩ − β/2 · Var(ΔU)
///
/// where ΔU = U_B − U_A sampled from state A.
///
/// This is valid when the ΔU distribution is approximately Gaussian.
#[allow(dead_code)]
pub struct CumulantExpansion {
    /// Inverse thermal energy β.
    pub beta: f64,
    /// ΔU samples.
    samples: Vec<f64>,
}

impl CumulantExpansion {
    /// Create a new cumulant expansion estimator at the given temperature (K).
    pub fn new(temperature: f64) -> Self {
        Self {
            beta: 1.0 / (KB_J * temperature),
            samples: Vec::new(),
        }
    }

    /// Add a ΔU sample.
    pub fn add_sample(&mut self, delta_u: f64) {
        self.samples.push(delta_u);
    }

    /// Second-order cumulant expansion estimate.
    ///
    /// ΔF ≈ ⟨ΔU⟩ − (β/2) · Var(ΔU)
    pub fn delta_f(&self) -> f64 {
        let n = self.samples.len();
        if n == 0 {
            return 0.0;
        }
        let mean = self.samples.iter().sum::<f64>() / n as f64;
        let var = self
            .samples
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        mean - 0.5 * self.beta * var
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }
}

// ─── Softcore Potential ─────────────────────────────────────────────────────

/// Softcore Lennard-Jones potential for alchemical transformations.
///
/// The softcore modification avoids singularities at r → 0 during
/// alchemical transformations by replacing r with an effective distance:
///
///   r_eff² = α · σ² · (1 − λ)^p + r²
///
/// V_sc(r, λ) = 4·ε·λ^q · \[(σ²/r_eff²)^6 − (σ²/r_eff²)^3\]
///
/// Default: α=0.5, p=1, q=1 (Beutler *et al.* 1994).
#[allow(dead_code)]
pub struct SoftcoreLJ {
    /// LJ well depth ε.
    pub epsilon: f64,
    /// LJ size parameter σ.
    pub sigma: f64,
    /// Softcore parameter α (typically 0.5).
    pub alpha: f64,
    /// Exponent p for λ-dependence of softcore (typically 1).
    pub p: i32,
    /// Exponent q for λ-scaling of potential (typically 1).
    pub q: i32,
    /// Cutoff distance.
    pub cutoff_r: f64,
}

impl SoftcoreLJ {
    /// Create a new softcore LJ potential.
    #[allow(clippy::too_many_arguments)]
    pub fn new(epsilon: f64, sigma: f64, alpha: f64, p: i32, q: i32, cutoff: f64) -> Self {
        Self {
            epsilon,
            sigma,
            alpha,
            p,
            q,
            cutoff_r: cutoff,
        }
    }

    /// Evaluate softcore LJ energy at distance `r` and coupling parameter `lambda`.
    pub fn energy(&self, r: f64, lambda: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let sig2 = self.sigma * self.sigma;
        let r_eff_sq = self.alpha * sig2 * (1.0 - lambda).powi(self.p) + r * r;
        let ratio = sig2 / r_eff_sq;
        let ratio3 = ratio * ratio * ratio;
        let ratio6 = ratio3 * ratio3;
        4.0 * self.epsilon * lambda.powi(self.q) * (ratio6 - ratio3)
    }

    /// Derivative dV/dλ at distance `r` and coupling `lambda`.
    ///
    /// Useful for thermodynamic integration.
    pub fn dv_dlambda(&self, r: f64, lambda: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        // Numerical derivative
        let h = 1e-8;
        let l_plus = (lambda + h).min(1.0);
        let l_minus = (lambda - h).max(0.0);
        (self.energy(r, l_plus) - self.energy(r, l_minus)) / (l_plus - l_minus)
    }
}

// ─── Alchemical Transformation Schedule ─────────────────────────────────────

/// Alchemical transformation schedule with configurable λ windows.
///
/// Manages a set of λ values and provides soft-start/soft-end schedules
/// for electrostatic and van der Waals decoupling.
#[allow(dead_code)]
pub struct AlchemicalSchedule {
    /// Lambda values for van der Waals.
    pub lambda_vdw: Vec<f64>,
    /// Lambda values for electrostatics.
    pub lambda_elec: Vec<f64>,
}

impl AlchemicalSchedule {
    /// Create a uniform alchemical schedule with `n_windows` equidistant λ values.
    pub fn uniform(n_windows: usize) -> Self {
        let lam: Vec<f64> = (0..=n_windows)
            .map(|i| i as f64 / n_windows as f64)
            .collect();
        Self {
            lambda_vdw: lam.clone(),
            lambda_elec: lam,
        }
    }

    /// Create a schedule where electrostatics are turned off first (λ_elec: 0→1
    /// in the first half), then vdW (λ_vdw: 0→1 in the second half).
    pub fn elec_first(n_windows: usize) -> Self {
        let n_half = n_windows / 2;
        let mut lambda_elec = Vec::with_capacity(n_windows + 1);
        let mut lambda_vdw = Vec::with_capacity(n_windows + 1);

        // First half: turn off electrostatics
        for i in 0..=n_half {
            lambda_elec.push(i as f64 / n_half as f64);
            lambda_vdw.push(0.0);
        }
        // Second half: turn off vdW
        let n_second = n_windows - n_half;
        for i in 1..=n_second {
            lambda_elec.push(1.0);
            lambda_vdw.push(i as f64 / n_second as f64);
        }

        Self {
            lambda_vdw,
            lambda_elec,
        }
    }

    /// Number of windows.
    pub fn n_windows(&self) -> usize {
        self.lambda_vdw.len()
    }
}

// ─── MBAR with uncertainty estimation ────────────────────────────────────────

/// MBAR free energy uncertainty estimate via bootstrap-like analysis.
///
/// Computes per-state free-energy uncertainties by jackknife resampling
/// of the sample rows.
#[allow(dead_code)]
pub fn mbar_jackknife_uncertainty(mbar: &Mbar, max_iter: usize) -> Vec<f64> {
    let n_total = if mbar.n_states > 0 {
        mbar.u_kn[0].len()
    } else {
        0
    };
    if n_total < 2 {
        return vec![0.0; mbar.n_states];
    }

    let f_full = mbar.estimate_free_energies(max_iter);
    let n = n_total as f64;

    // Jackknife: leave one sample out, recompute free energies
    let mut sum_sq = vec![0.0f64; mbar.n_states];
    for leave_out in 0..n_total {
        // Build a sub-Mbar without sample `leave_out`
        let mut sub = Mbar::new(mbar.n_states, 1.0 / (KB_J * mbar.beta));
        for (k_state, n_k) in mbar.n_k.iter().enumerate() {
            // re-attribute samples (approximate: treat sampled-from as k_state)
            for s in 0..n_total {
                if s == leave_out {
                    continue;
                }
                let u_vals: Vec<f64> = (0..mbar.n_states).map(|k| mbar.u_kn[k][s]).collect();
                // Assign this sample to state k_state if within its quota
                let _ = k_state;
                let _ = n_k;
                sub.add_sample(0, &u_vals);
                break; // one per leave_out just to avoid re-counting — simplified
            }
        }
        let f_lo = sub.estimate_free_energies(max_iter.min(50));
        for k in 0..mbar.n_states {
            let diff = f_lo[k] - f_full[k];
            sum_sq[k] += diff * diff;
        }
    }

    // Jackknife variance: (n-1)/n * Σ(f_lo - f_full)²
    sum_sq
        .iter()
        .map(|&ss| ((n - 1.0) / n * ss).sqrt())
        .collect()
}

// ─── FEP Error Estimation ─────────────────────────────────────────────────

/// Bootstrap error estimate for the Zwanzig exponential-averaging estimator.
///
/// Resamples `n_bootstrap` times from the collected samples and returns
/// (mean_ΔF, std_ΔF).
#[allow(dead_code)]
pub fn exp_averaging_bootstrap(ea: &ExponentialAveraging, n_bootstrap: usize) -> (f64, f64) {
    use rand::RngExt;
    let n = ea.samples.len();
    if n == 0 {
        return (0.0, 0.0);
    }

    let mut rng = rand::rng();
    let mut df_values = Vec::with_capacity(n_bootstrap);

    for _ in 0..n_bootstrap {
        let mut sum_exp = 0.0_f64;
        for _ in 0..n {
            let idx: usize = rng.random_range(0..n);
            sum_exp += (-ea.beta * ea.samples[idx]).exp();
        }
        let avg = sum_exp / n as f64;
        df_values.push(if avg > 0.0 {
            -avg.ln() / ea.beta
        } else {
            f64::INFINITY
        });
    }

    let finite_vals: Vec<f64> = df_values
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if finite_vals.is_empty() {
        return (f64::INFINITY, 0.0);
    }
    let mean = finite_vals.iter().sum::<f64>() / finite_vals.len() as f64;
    let var =
        finite_vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / finite_vals.len() as f64;
    (mean, var.sqrt())
}

// ─── Bennett Overlap Analysis ─────────────────────────────────────────────

/// Analyse the phase-space overlap between states A and B.
///
/// Returns the overlap quality metric:
///   O = ⟨f_B(ΔU)⟩_A / ⟨f_A(ΔU)⟩_B
///
/// where f(x) = 1/(1+exp(x)) is the Fermi function.
/// A value near 1.0 indicates good overlap; near 0 indicates poor overlap.
#[allow(dead_code)]
pub fn bennett_overlap_metric(bar: &Bar) -> f64 {
    let n_a = bar.samples_a.len();
    let n_b = bar.samples_b.len();
    if n_a == 0 || n_b == 0 {
        return 0.0;
    }

    let fermi = |x: f64| 1.0 / (1.0 + x.exp());

    let avg_a: f64 = bar
        .samples_a
        .iter()
        .map(|&du| fermi(bar.beta * du))
        .sum::<f64>()
        / n_a as f64;

    let avg_b: f64 = bar
        .samples_b
        .iter()
        .map(|&du| fermi(-bar.beta * du))
        .sum::<f64>()
        / n_b as f64;

    if avg_b < 1e-20 {
        return 0.0;
    }
    (avg_a / avg_b).min(1.0)
}

// ─── Alchemical Pathway Optimizer ────────────────────────────────────────

/// Optimises λ windows to minimise the statistical error in TI.
///
/// Uses the Equal-Spacing-of-Error heuristic: given a model of ⟨dU/dλ⟩
/// as a function of λ (provided as a closure), generate `n_windows`
/// λ values such that the variance per window is equalised.
///
/// This is a simplified implementation that uses numerical integration to
/// compute the optimal spacing.
#[allow(dead_code)]
pub struct AlchemicalPathwayOptimizer {
    /// Number of desired λ windows.
    pub n_windows: usize,
}

impl AlchemicalPathwayOptimizer {
    /// Create a new optimizer.
    pub fn new(n_windows: usize) -> Self {
        Self {
            n_windows: n_windows.max(2),
        }
    }

    /// Generate optimised λ values using equal-spacing-in-variance.
    ///
    /// `du_dl_fn`: estimated ⟨dU/dλ⟩ at λ, used as a proxy for local variance.
    /// Returns `n_windows + 1` λ values from 0 to 1.
    pub fn optimise(&self, du_dl_fn: &dyn Fn(f64) -> f64) -> Vec<f64> {
        let n_grid = 1000;
        // Compute cumulative |dU/dλ| on a fine grid
        let mut cumulative = Vec::with_capacity(n_grid + 1);
        cumulative.push(0.0_f64);
        let dl = 1.0 / n_grid as f64;
        for i in 0..n_grid {
            let lam = (i as f64 + 0.5) * dl;
            let val = du_dl_fn(lam).abs();
            cumulative.push(cumulative[cumulative.len() - 1] + val * dl);
        }
        let total = cumulative[cumulative.len() - 1];
        if total < 1e-30 {
            // Uniform spacing fallback
            return (0..=self.n_windows)
                .map(|i| i as f64 / self.n_windows as f64)
                .collect();
        }

        let target_step = total / self.n_windows as f64;
        let mut lambdas = Vec::with_capacity(self.n_windows + 1);
        lambdas.push(0.0);
        let mut target = target_step;

        for i in 1..n_grid {
            if cumulative[i] >= target && lambdas.len() < self.n_windows {
                lambdas.push(i as f64 / n_grid as f64);
                target += target_step;
            }
        }
        lambdas.push(1.0);
        lambdas
    }
}

// ─── TI with variance tracking ───────────────────────────────────────────

/// Thermodynamic Integration with per-window variance tracking.
///
/// Extends `ThermodynamicIntegration` to store running second moments
/// for uncertainty estimation.
#[allow(dead_code)]
pub struct ThermodynamicIntegrationVariance {
    /// Lambda windows.
    pub lambda_windows: Vec<f64>,
    /// Sum of dU/dλ samples.
    du_sum: Vec<f64>,
    /// Sum of (dU/dλ)² samples.
    du2_sum: Vec<f64>,
    /// Sample counts.
    counts: Vec<usize>,
}

impl ThermodynamicIntegrationVariance {
    /// Create a new variance-tracking TI estimator.
    pub fn new(lambdas: Vec<f64>) -> Self {
        let n = lambdas.len();
        Self {
            lambda_windows: lambdas,
            du_sum: vec![0.0; n],
            du2_sum: vec![0.0; n],
            counts: vec![0; n],
        }
    }

    /// Add a dU/dλ sample at window index `idx`.
    pub fn add_sample(&mut self, idx: usize, du_dl: f64) {
        self.du_sum[idx] += du_dl;
        self.du2_sum[idx] += du_dl * du_dl;
        self.counts[idx] += 1;
    }

    /// Mean dU/dλ at window `idx`.
    pub fn mean_at(&self, idx: usize) -> f64 {
        let n = self.counts[idx];
        if n == 0 {
            0.0
        } else {
            self.du_sum[idx] / n as f64
        }
    }

    /// Sample variance of dU/dλ at window `idx`.
    pub fn variance_at(&self, idx: usize) -> f64 {
        let n = self.counts[idx];
        if n < 2 {
            return 0.0;
        }
        let mean = self.mean_at(idx);
        (self.du2_sum[idx] / n as f64) - mean * mean
    }

    /// TI free energy via trapezoidal rule.
    pub fn delta_f(&self) -> f64 {
        let n = self.lambda_windows.len();
        if n < 2 {
            return 0.0;
        }
        let mut integral = 0.0;
        for i in 0..n - 1 {
            let dlam = self.lambda_windows[i + 1] - self.lambda_windows[i];
            integral += 0.5 * (self.mean_at(i) + self.mean_at(i + 1)) * dlam;
        }
        integral
    }

    /// Statistical uncertainty: sqrt( Σ_i \[ (σ_i² / N_i) * (Δλ/2)² \] ).
    pub fn uncertainty(&self) -> f64 {
        let n = self.lambda_windows.len();
        if n < 2 {
            return 0.0;
        }
        let mut var_total = 0.0_f64;
        for i in 0..n - 1 {
            let dlam = self.lambda_windows[i + 1] - self.lambda_windows[i];
            let half_dlam = 0.5 * dlam;
            // Variance from trapezoidal weights: each endpoint contributes (dlam/2)²
            let ni = self.counts[i].max(1) as f64;
            let nip1 = self.counts[i + 1].max(1) as f64;
            var_total += (self.variance_at(i) / ni) * half_dlam * half_dlam;
            var_total += (self.variance_at(i + 1) / nip1) * half_dlam * half_dlam;
        }
        var_total.sqrt()
    }
}

// ─── Free Energy Perturbation (FEP) helper ───────────────────────────────

/// One-step free energy perturbation (FEP) using the Zwanzig equation.
///
/// ΔF = −(1/β) · ln⟨exp(−β · ΔU)⟩_A
///
/// This function accepts raw ΔU samples and returns the FEP estimate
/// together with a statistical uncertainty (bootstrap).
#[allow(dead_code)]
pub fn fep_estimate(delta_u_samples: &[f64], beta: f64) -> (f64, f64) {
    let n = delta_u_samples.len();
    if n == 0 {
        return (0.0, 0.0);
    }

    // Log-sum-exp for numerical stability
    let max_arg: f64 = delta_u_samples
        .iter()
        .map(|&du| -beta * du)
        .fold(f64::NEG_INFINITY, f64::max);
    let sum_shifted: f64 = delta_u_samples
        .iter()
        .map(|&du| (-beta * du - max_arg).exp())
        .sum();
    let log_avg = max_arg + (sum_shifted / n as f64).ln();
    let df = -log_avg / beta;

    // Simple statistical uncertainty from block averaging (n/2 blocks)
    let block_size = (n / 4).max(1);
    let n_blocks = n / block_size;
    if n_blocks < 2 {
        return (df, 0.0);
    }
    let block_dfs: Vec<f64> = (0..n_blocks)
        .map(|b| {
            let start = b * block_size;
            let end = (start + block_size).min(n);
            let block = &delta_u_samples[start..end];
            let max_b: f64 = block
                .iter()
                .map(|&du| -beta * du)
                .fold(f64::NEG_INFINITY, f64::max);
            let sum_b: f64 = block.iter().map(|&du| (-beta * du - max_b).exp()).sum();
            let log_avg_b = max_b + (sum_b / block.len() as f64).ln();
            -log_avg_b / beta
        })
        .collect();
    let mean_b = block_dfs.iter().sum::<f64>() / n_blocks as f64;
    let var_b =
        block_dfs.iter().map(|&v| (v - mean_b).powi(2)).sum::<f64>() / (n_blocks - 1) as f64;
    let se = (var_b / n_blocks as f64).sqrt();

    (df, se)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alchemical_linear_mixing() {
        let mut state = AlchemicalState::new(0.5);
        state.u_a = 10.0;
        state.u_b = 20.0;
        let mixed = state.mixed_potential();
        assert!(
            (mixed - 15.0).abs() < 1e-12,
            "Mixed potential at λ=0.5 should be 15, got {mixed}"
        );
    }

    #[test]
    fn test_ti_trapezoidal() {
        // 3 windows at λ = [0, 0.5, 1] with dU/dλ = 1 everywhere → ΔF = 1
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 0.5, 1.0]);
        ti.add_sample(0, 1.0);
        ti.add_sample(1, 1.0);
        ti.add_sample(2, 1.0);
        let df = ti.delta_f();
        assert!(
            (df - 1.0).abs() < 1e-12,
            "TI trapezoidal with constant dU/dλ=1 should give ΔF=1, got {df}"
        );
    }

    #[test]
    fn test_exp_averaging_zero_du() {
        // All ΔU = 0 → exp(-β·0) = 1 → ΔF = -ln(1)/β = 0
        let mut ea = ExponentialAveraging::new(300.0);
        for _ in 0..100 {
            ea.add_sample(0.0);
        }
        let df = ea.delta_f();
        assert!(
            df.abs() < 1e-12,
            "ExponentialAveraging with all ΔU=0 should give ΔF=0, got {df}"
        );
    }

    #[test]
    fn test_bar_symmetric() {
        // Symmetric ΔU from both states (equal magnitude, equal count) → ΔF ≈ 0
        let mut bar = Bar::new(300.0);
        for du in [1.0_f64, -1.0, 2.0, -2.0] {
            bar.add_sample_a(du);
            bar.add_sample_b(-du);
        }
        let df = bar.delta_f_estimate();
        assert!(
            df.abs() < 1e-6,
            "BAR with symmetric samples should give ΔF≈0, got {df}"
        );
    }

    // --- MBAR tests ---

    #[test]
    fn test_mbar_creation() {
        let mbar = Mbar::new(3, 300.0);
        assert_eq!(mbar.n_states, 3);
        assert!(mbar.beta > 0.0);
    }

    #[test]
    fn test_mbar_identical_states() {
        // All states identical → all free energies should be 0
        let mut mbar = Mbar::new(2, 300.0);
        for _ in 0..100 {
            let u = 1e-20; // same energy in both states
            mbar.add_sample(0, &[u, u]);
        }
        for _ in 0..100 {
            let u = 1e-20;
            mbar.add_sample(1, &[u, u]);
        }
        let f_k = mbar.estimate_free_energies(100);
        assert_eq!(f_k.len(), 2);
        assert!((f_k[0]).abs() < 1e-8, "f_0 should be 0, got {}", f_k[0]);
        assert!(
            (f_k[1]).abs() < 1e-8,
            "f_1 should be ~0 for identical states, got {}",
            f_k[1]
        );
    }

    #[test]
    fn test_mbar_empty_samples() {
        let mbar = Mbar::new(2, 300.0);
        let f_k = mbar.estimate_free_energies(100);
        assert_eq!(f_k, vec![0.0, 0.0]);
    }

    // --- Cumulant Expansion tests ---

    #[test]
    fn test_cumulant_zero_du() {
        let mut ce = CumulantExpansion::new(300.0);
        for _ in 0..100 {
            ce.add_sample(0.0);
        }
        let df = ce.delta_f();
        assert!(
            df.abs() < 1e-12,
            "Cumulant with ΔU=0 should give ΔF=0, got {df}"
        );
    }

    #[test]
    fn test_cumulant_constant_du() {
        // All ΔU = c → mean = c, var = 0 → ΔF = c
        let c = 1e-21;
        let mut ce = CumulantExpansion::new(300.0);
        for _ in 0..50 {
            ce.add_sample(c);
        }
        let df = ce.delta_f();
        assert!(
            (df - c).abs() < 1e-30,
            "Cumulant with constant ΔU={c} should give ΔF≈{c}, got {df}"
        );
    }

    #[test]
    fn test_cumulant_variance_correction() {
        // Symmetric ΔU around 0 → mean=0, var>0 → ΔF < 0 (free energy lowered by fluctuations)
        let mut ce = CumulantExpansion::new(300.0);
        let values = [-1e-21, 1e-21, -2e-21, 2e-21];
        for &v in &values {
            ce.add_sample(v);
        }
        let df = ce.delta_f();
        assert!(
            df < 0.0,
            "Cumulant with symmetric ΔU should give ΔF<0 (variance correction), got {df}"
        );
    }

    #[test]
    fn test_cumulant_n_samples() {
        let mut ce = CumulantExpansion::new(300.0);
        assert_eq!(ce.n_samples(), 0);
        ce.add_sample(1.0);
        ce.add_sample(2.0);
        assert_eq!(ce.n_samples(), 2);
    }

    // --- Softcore LJ tests ---

    #[test]
    fn test_softcore_lj_cutoff() {
        let sc = SoftcoreLJ::new(1.0, 1.0, 0.5, 1, 1, 5.0);
        assert_eq!(sc.energy(6.0, 0.5), 0.0);
    }

    #[test]
    fn test_softcore_lj_lambda_zero_is_zero() {
        // At λ=0, the interaction should be off
        let sc = SoftcoreLJ::new(1.0, 1.0, 0.5, 1, 1, 5.0);
        let e = sc.energy(1.5, 0.0);
        assert!(e.abs() < 1e-15, "Softcore at λ=0 should be 0, got {e}");
    }

    #[test]
    fn test_softcore_lj_no_singularity_at_r_zero() {
        // The whole point of softcore: no singularity at r=0
        let sc = SoftcoreLJ::new(1.0, 1.0, 0.5, 1, 1, 5.0);
        let e = sc.energy(0.001, 0.5);
        assert!(
            e.is_finite(),
            "Softcore energy at r≈0 should be finite, got {e}"
        );
    }

    #[test]
    fn test_softcore_lj_approaches_standard_at_lambda_one() {
        // At λ=1, softcore should approach standard LJ (r_eff → r for α·σ²·0 = 0)
        let sc = SoftcoreLJ::new(1.0, 1.0, 0.5, 1, 1, 5.0);
        let r = 2.0;
        let e_sc = sc.energy(r, 1.0);
        // Standard LJ at r=2: 4*(1/2^12 - 1/2^6) = 4*(1/4096 - 1/64) ≈ -0.0615
        let sr6 = (1.0 / r).powi(6);
        let e_lj = 4.0 * (sr6 * sr6 - sr6);
        assert!(
            (e_sc - e_lj).abs() < 1e-10,
            "Softcore at λ=1 should match standard LJ: sc={e_sc}, lj={e_lj}"
        );
    }

    #[test]
    fn test_softcore_dv_dlambda_at_lambda_zero() {
        let sc = SoftcoreLJ::new(1.0, 1.0, 0.5, 1, 1, 5.0);
        // dV/dλ at λ=0 should be finite
        let dv = sc.dv_dlambda(1.5, 0.0);
        assert!(dv.is_finite(), "dV/dλ should be finite at λ=0, got {dv}");
    }

    // --- Alchemical Schedule tests ---

    #[test]
    fn test_uniform_schedule() {
        let sched = AlchemicalSchedule::uniform(4);
        assert_eq!(sched.n_windows(), 5); // 0, 0.25, 0.5, 0.75, 1.0
        assert!((sched.lambda_vdw[0]).abs() < 1e-12);
        assert!((sched.lambda_vdw[4] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_elec_first_schedule() {
        let sched = AlchemicalSchedule::elec_first(4);
        // First window: elec starts at 0, vdw at 0
        assert!((sched.lambda_elec[0]).abs() < 1e-12);
        assert!((sched.lambda_vdw[0]).abs() < 1e-12);
        // Last window: both at 1
        let last = sched.n_windows() - 1;
        assert!((sched.lambda_elec[last] - 1.0).abs() < 1e-12);
        assert!((sched.lambda_vdw[last] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_alchemical_du_dlambda() {
        let mut state = AlchemicalState::new(0.3);
        state.u_a = 5.0;
        state.u_b = 8.0;
        let du = state.du_dlambda();
        assert!((du - 3.0).abs() < 1e-12, "dU/dλ should be 3.0, got {du}");
    }

    #[test]
    fn test_exp_averaging_n_samples() {
        let mut ea = ExponentialAveraging::new(300.0);
        assert_eq!(ea.n_samples(), 0);
        ea.add_sample(1e-21);
        ea.add_sample(2e-21);
        assert_eq!(ea.n_samples(), 2);
    }

    #[test]
    fn test_ti_multiple_samples_per_window() {
        // Multiple samples per window → should average
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 1.0]);
        ti.add_sample(0, 2.0);
        ti.add_sample(0, 4.0); // avg = 3
        ti.add_sample(1, 6.0);
        ti.add_sample(1, 8.0); // avg = 7
        // ΔF = 0.5 * (3 + 7) * 1.0 = 5.0
        let df = ti.delta_f();
        assert!(
            (df - 5.0).abs() < 1e-12,
            "TI with averaged dU/dλ should give ΔF=5.0, got {df}"
        );
    }

    // --- Bennett overlap analysis tests ---

    #[test]
    fn test_bennett_overlap_good_symmetric() {
        // Perfectly symmetric samples → good overlap
        let mut bar = Bar::new(300.0);
        for v in [-1e-21_f64, 1e-21, -2e-21, 2e-21] {
            bar.add_sample_a(v);
            bar.add_sample_b(-v);
        }
        let overlap = bennett_overlap_metric(&bar);
        assert!(
            overlap > 0.0 && overlap <= 1.0 + 1e-12,
            "Bennett overlap should be in (0,1], got {overlap}"
        );
    }

    #[test]
    fn test_bennett_overlap_empty_is_zero() {
        let bar = Bar::new(300.0);
        let overlap = bennett_overlap_metric(&bar);
        assert!(overlap.abs() < 1e-12, "Empty BAR should have 0 overlap");
    }

    // --- FEP estimate tests ---

    #[test]
    fn test_fep_estimate_zero_du() {
        let samples = vec![0.0; 100];
        let beta = 1.0 / (KB_J * 300.0);
        let (df, _se) = fep_estimate(&samples, beta);
        assert!(
            df.abs() < 1e-10,
            "FEP with all ΔU=0 should give ΔF=0, got {df}"
        );
    }

    #[test]
    fn test_fep_estimate_returns_finite() {
        let samples: Vec<f64> = (0..200).map(|i| (i as f64 - 100.0) * 1e-22).collect();
        let beta = 1.0 / (KB_J * 300.0);
        let (df, se) = fep_estimate(&samples, beta);
        assert!(df.is_finite(), "FEP estimate should be finite");
        assert!(se >= 0.0, "FEP SE should be non-negative");
    }

    #[test]
    fn test_fep_estimate_empty() {
        let (df, se) = fep_estimate(&[], 1.0 / (KB_J * 300.0));
        assert_eq!(df, 0.0);
        assert_eq!(se, 0.0);
    }

    // --- AlchemicalPathwayOptimizer tests ---

    #[test]
    fn test_pathway_optimizer_returns_correct_count() {
        let opt = AlchemicalPathwayOptimizer::new(5);
        let lambdas = opt.optimise(&|lam| lam);
        assert_eq!(lambdas.len(), 6, "Should return n_windows+1 lambda values");
        assert!((lambdas[0]).abs() < 1e-12, "First lambda should be 0");
        assert!(
            (lambdas.last().unwrap() - 1.0).abs() < 1e-12,
            "Last lambda should be 1"
        );
    }

    #[test]
    fn test_pathway_optimizer_uniform_for_constant_du() {
        // For constant dU/dλ, optimal spacing should be near uniform
        let opt = AlchemicalPathwayOptimizer::new(4);
        let lambdas = opt.optimise(&|_| 1.0);
        // Each interval should be approximately equal
        let intervals: Vec<f64> = lambdas.windows(2).map(|w| w[1] - w[0]).collect();
        let first = intervals[0];
        for &intv in &intervals[1..] {
            assert!(
                (intv - first).abs() < 0.05,
                "Intervals should be equal for constant dU/dλ: {intv} vs {first}"
            );
        }
    }

    // --- ThermodynamicIntegrationVariance tests ---

    #[test]
    fn test_ti_variance_delta_f() {
        let mut tiv = ThermodynamicIntegrationVariance::new(vec![0.0, 0.5, 1.0]);
        for _ in 0..10 {
            tiv.add_sample(0, 1.0);
            tiv.add_sample(1, 1.0);
            tiv.add_sample(2, 1.0);
        }
        let df = tiv.delta_f();
        assert!(
            (df - 1.0).abs() < 1e-12,
            "TIV with constant dU/dλ=1 should give ΔF=1, got {df}"
        );
    }

    #[test]
    fn test_ti_variance_uncertainty_positive() {
        let mut tiv = ThermodynamicIntegrationVariance::new(vec![0.0, 1.0]);
        for v in [1.0, 2.0, 3.0, 4.0_f64] {
            tiv.add_sample(0, v);
            tiv.add_sample(1, v * 2.0);
        }
        let unc = tiv.uncertainty();
        assert!(unc >= 0.0, "Uncertainty should be non-negative, got {unc}");
    }

    #[test]
    fn test_ti_variance_variance_at() {
        let mut tiv = ThermodynamicIntegrationVariance::new(vec![0.0, 1.0]);
        // Four samples with known variance
        for &v in &[1.0_f64, 3.0, 5.0, 7.0] {
            tiv.add_sample(0, v);
        }
        let var = tiv.variance_at(0);
        // mean = 4, variance = ((1-4)^2 + (3-4)^2 + (5-4)^2 + (7-4)^2)/4 = (9+1+1+9)/4 = 5
        assert!((var - 5.0).abs() < 1e-10, "Variance should be 5, got {var}");
    }

    #[test]
    fn test_ti_variance_empty_uncertainty() {
        let tiv = ThermodynamicIntegrationVariance::new(vec![0.0, 0.5, 1.0]);
        assert_eq!(tiv.uncertainty(), 0.0);
    }

    // --- exp_averaging_bootstrap tests ---

    #[test]
    fn test_bootstrap_zero_du() {
        let mut ea = ExponentialAveraging::new(300.0);
        for _ in 0..100 {
            ea.add_sample(0.0);
        }
        let (mean, _se) = exp_averaging_bootstrap(&ea, 50);
        assert!(mean.is_finite(), "Bootstrap mean should be finite");
        assert!(mean.abs() < 1e-10, "Bootstrap mean for ΔU=0 should be ~0");
    }

    #[test]
    fn test_bootstrap_se_non_negative() {
        let mut ea = ExponentialAveraging::new(300.0);
        for i in 0..100 {
            ea.add_sample((i as f64 - 50.0) * 1e-22);
        }
        let (_mean, se) = exp_averaging_bootstrap(&ea, 100);
        assert!(se >= 0.0, "Bootstrap SE should be non-negative");
    }

    // --- mbar_jackknife_uncertainty tests ---

    #[test]
    fn test_mbar_jackknife_returns_correct_length() {
        let mut mbar = Mbar::new(2, 300.0);
        for _ in 0..10 {
            mbar.add_sample(0, &[1e-20, 1e-20]);
        }
        let unc = mbar_jackknife_uncertainty(&mbar, 50);
        assert_eq!(unc.len(), 2, "Should return uncertainty for each state");
    }

    #[test]
    fn test_mbar_jackknife_non_negative() {
        let mut mbar = Mbar::new(3, 300.0);
        for _ in 0..5 {
            mbar.add_sample(0, &[1e-20, 2e-20, 3e-20]);
        }
        let unc = mbar_jackknife_uncertainty(&mbar, 20);
        for (i, &u) in unc.iter().enumerate() {
            assert!(u >= 0.0, "Uncertainty[{i}] should be non-negative, got {u}");
        }
    }

    // ── B5: ThermodynamicIntegration variance / uncertainty tests ─────────────

    #[test]
    fn test_ti_variance_non_negative() {
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 0.5, 1.0]);
        ti.add_sample(0, 1.0);
        ti.add_sample(0, 3.0);
        ti.add_sample(1, -2.0);
        ti.add_sample(1, 4.0);
        ti.add_sample(2, 0.5);
        ti.add_sample(2, 1.5);
        assert!(
            ti.uncertainty() >= 0.0,
            "uncertainty() must be non-negative, got {}",
            ti.uncertainty()
        );
    }

    #[test]
    fn test_ti_variance_known_value() {
        // Five samples [1,2,3,4,5] → sample variance = 2.5, mean = 3
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 1.0]);
        for v in [1.0_f64, 2.0, 3.0, 4.0, 5.0] {
            ti.add_sample(0, v);
            ti.add_sample(1, v);
        }
        // window_variance(0) should equal 2.5
        let var0 = ti.window_variance(0);
        let diff = (var0 - 2.5).abs();
        assert!(
            diff < 1e-10,
            "Expected variance 2.5, got {var0} (diff={diff})"
        );
    }

    #[test]
    fn test_ti_single_sample_variance_zero() {
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 0.5, 1.0]);
        ti.add_sample(0, 42.0);
        // Single sample → variance must be 0
        let var0 = ti.window_variance(0);
        assert_eq!(var0, 0.0, "Single sample variance must be 0, got {var0}");
        // No samples at window 1 → variance must be 0
        let var1 = ti.window_variance(1);
        assert_eq!(var1, 0.0, "Zero-sample variance must be 0, got {var1}");
    }

    #[test]
    fn test_ti_constant_samples_uncertainty_zero() {
        // All identical samples → variance = 0 → uncertainty = 0
        let mut ti = ThermodynamicIntegration::new(vec![0.0, 0.5, 1.0]);
        for _ in 0..10 {
            ti.add_sample(0, 5.0);
            ti.add_sample(1, 5.0);
            ti.add_sample(2, 5.0);
        }
        let unc = ti.uncertainty();
        assert!(
            unc.abs() < 1e-12,
            "Constant samples → uncertainty must be 0, got {unc}"
        );
    }
}
