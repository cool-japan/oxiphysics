//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::error::MdError;
use rand::rngs::StdRng;
use rand::{Rng, RngExt, SeedableRng};

/// Metropolis acceptance test for a move whose reduced cost is `delta`.
///
/// `delta` is the change of the dimensionless energy `β E`, so `exp(-delta)` is
/// the Boltzmann factor, and the move is accepted with probability
/// `min(1, exp(-delta))`.  Exactly one uniform random number is consumed per
/// call regardless of `delta`, which keeps seeded streams aligned between runs
/// with different bias landscapes.
fn metropolis_accept<R: Rng + ?Sized>(delta: f64, rng: &mut R) -> bool {
    let u: f64 = rng.random();
    u < (-delta).exp().min(1.0)
}
/// Check that `index` addresses one of `n_replicas` replicas.
fn check_replica_index(index: usize, n_replicas: usize) -> Result<(), MdError> {
    if index < n_replicas {
        Ok(())
    } else {
        Err(MdError::General(format!(
            "replica index {index} out of range for {n_replicas} replicas"
        )))
    }
}
/// Inverse thermal energy `1 / (k_B T)` in mol/kJ, rejecting non-positive or
/// non-finite temperatures.
fn inverse_thermal_energy(temperature: f64) -> Result<f64, MdError> {
    if temperature.is_finite() && temperature > 0.0 {
        Ok(1.0 / (KB * temperature))
    } else {
        Err(MdError::General(format!(
            "replica temperature must be positive and finite, got {temperature}"
        )))
    }
}
/// State for a metadynamics (or well-tempered metadynamics) simulation.
pub struct MetadynamicsState {
    /// Deposited Gaussian hills.
    pub hills: Vec<GaussianHill>,
    /// Metadynamics parameters.
    pub params: MetadynamicsParams,
    /// Current simulation step count.
    pub step_count: usize,
}
impl MetadynamicsState {
    /// Create a new metadynamics state.
    pub fn new(params: MetadynamicsParams) -> Self {
        Self {
            hills: Vec::new(),
            params,
            step_count: 0,
        }
    }
    /// Compute total bias potential at `cv_values`:
    /// V_bias = sum of all Gaussian hills at current CV values
    pub fn bias_potential(&self, cv_values: &[f64]) -> f64 {
        self.hills.iter().map(|h| h.evaluate(cv_values)).sum()
    }
    /// Compute bias force (negative gradient of bias potential) in CV space.
    /// F_bias_s_i = -dV_bias/ds_i
    pub fn bias_force(&self, cv_values: &[f64]) -> Vec<f64> {
        let n = cv_values.len();
        let mut force = vec![0.0f64; n];
        for hill in &self.hills {
            let grad = hill.gradient_wrt_cv(cv_values);
            for (i, g) in grad.iter().enumerate() {
                force[i] -= g;
            }
        }
        force
    }
    /// Conditionally deposit a new hill, then increment step count.
    ///
    /// Deposits when `step_count % deposition_stride == 0`.
    ///
    /// For well-tempered metadynamics (`delta_t > 0`), the height is scaled by:
    ///   height *= exp(-V_bias(s) / (k_B * delta_T))
    ///   where delta_T = delta_t * temperature
    pub fn maybe_deposit(&mut self, cv_values: &[f64]) {
        if self.params.deposition_stride > 0
            && self
                .step_count
                .is_multiple_of(self.params.deposition_stride)
        {
            let mut h = self.params.height;
            if self.params.delta_t > 0.0 {
                let v_bias = self.bias_potential(cv_values);
                let delta_temp = self.params.delta_t * self.params.temperature;
                h *= (-v_bias / (KB * delta_temp)).exp();
            }
            let hill = GaussianHill::new(cv_values.to_vec(), h, self.params.widths.clone());
            self.hills.push(hill);
        }
        self.step_count += 1;
    }
    /// Estimate the free energy on a 1D grid.
    ///
    /// In well-tempered metadynamics, F(s) ≈ -V_bias(s) as hills converge.
    /// `bins` is a slice of 1D CV grid points (each a `Vec`f64` of length 1).
    pub fn free_energy_estimate(&self, bins: &[Vec<f64>]) -> Vec<f64> {
        bins.iter().map(|s| -self.bias_potential(s)).collect()
    }
}
/// OPES state (on-the-fly probability enhanced sampling): a reweighted kernel
/// density estimate of the *unbiased* CV distribution `P(s)` and the bias built
/// from it,
///
/// ```text
/// V(s) = (1 - 1/gamma) * k_B*T * ln( P_hat(s) / Z + epsilon )
/// ```
///
/// which drives the CV towards the well-tempered target `P(s)^(1/gamma)`.
///
/// The samples `s_k` are drawn from the *biased* ensemble, so every kernel is
/// weighted by the reweighting factor `w_k = exp(V(s_k) / k_B*T)` evaluated with
/// the bias that was active when the sample was taken.  The kernel height stored
/// in [`OpesState::kernels`] is that weight; without it the estimate converges
/// to the biased distribution and the bias flattens instead of reaching the
/// target.  `Z` is the average of the (un-normalised) density estimate over the
/// kernel centres, so `P_hat(s) / Z` is the density relative to its mean over
/// the explored region, and `epsilon = exp(-barrier / (k_B*T * (1 - 1/gamma)))`
/// bounds the bias from below by `-barrier`.
///
/// Reference: Invernizzi & Parrinello, J. Phys. Chem. Lett. 11, 2731 (2020).
#[derive(Debug, Clone)]
pub struct OpesState {
    /// Gaussian kernels forming the density estimate.  The height of each kernel
    /// is its reweighting factor `w_k = exp(V(s_k) / k_B*T)`.
    pub kernels: Vec<GaussianHill>,
    /// Normalisation constant Z: the mean of [`OpesState::kernel_density`] over
    /// the kernel centres.  Kept up to date by [`OpesState::update`]; call
    /// [`OpesState::recompute_normalisation`] after editing `kernels` by hand.
    pub z: f64,
    /// Bias factor gamma.  `gamma = 1` is unbiased sampling and no bias is
    /// applied for `gamma <= 1`; the target is `P(s)^(1/gamma)`, so `gamma -> inf`
    /// aims at a uniform distribution of the CV.
    pub gamma: f64,
    /// Simulation temperature (K).
    pub temperature: f64,
    /// Step counter.
    pub step: usize,
    /// Deposition stride.
    pub stride: usize,
    /// Kernel bandwidth: Gaussian standard deviation per CV dimension.
    pub bandwidth: Vec<f64>,
    /// Free-energy barrier `ΔE` (kJ/mol) the sampling should be able to cross.
    /// It sets the regularisation `epsilon` (the bias never drops below
    /// `-barrier`).  [`OpesState::new`] uses `gamma * k_B*T`, the convention in
    /// which the bias factor is `barrier / k_B*T`.
    pub barrier: f64,
    /// `sum_k kernel_density(centre_k)`, the numerator of `z`.
    centre_density_sum: f64,
}
impl OpesState {
    /// Create a new OPES state.
    ///
    /// The barrier is set to `gamma * k_B*T`; use [`OpesState::with_barrier`] to
    /// choose it independently of the bias factor.
    pub fn new(gamma: f64, temperature: f64, stride: usize, bandwidth: Vec<f64>) -> Self {
        Self {
            kernels: Vec::new(),
            z: 0.0,
            gamma,
            temperature,
            step: 0,
            stride,
            bandwidth,
            barrier: gamma * KB * temperature,
            centre_density_sum: 0.0,
        }
    }
    /// Set the barrier `ΔE` (kJ/mol) that defines the regularisation `epsilon`.
    ///
    /// # Errors
    /// [`MdError::General`] if `barrier` is NaN or negative (`+inf` is allowed
    /// and gives `epsilon = 0`).
    pub fn with_barrier(mut self, barrier: f64) -> Result<Self, MdError> {
        if barrier.is_nan() || barrier < 0.0 {
            return Err(MdError::General(format!(
                "OPES barrier must be non-negative, got {barrier}"
            )));
        }
        self.barrier = barrier;
        Ok(self)
    }
    /// Thermal energy `k_B*T` in kJ/mol.
    fn thermal_energy(&self) -> f64 {
        KB * self.temperature
    }
    /// Prefactor `1 - 1/gamma` of the bias.
    fn bias_prefactor(&self) -> f64 {
        1.0 - 1.0 / self.gamma
    }
    /// Whether a bias is applied: at least one kernel, `gamma > 1` and a usable
    /// temperature.
    fn is_biasing(&self) -> bool {
        let kt = self.thermal_energy();
        self.z >= 1e-30 && self.gamma > 1.0 && kt > 0.0 && kt.is_finite()
    }
    /// Regularisation `epsilon = exp(-barrier / (k_B*T * (1 - 1/gamma)))`.
    ///
    /// It is the value of `P_hat / Z` below which the bias stops decreasing:
    /// far from every kernel `V = (1 - 1/gamma) k_B*T ln(epsilon) = -barrier`.
    pub fn epsilon(&self) -> f64 {
        let scale = self.bias_prefactor() * self.thermal_energy();
        if scale > 0.0 && scale.is_finite() {
            (-self.barrier / scale).exp().max(0.0)
        } else {
            0.0
        }
    }
    /// `P_hat(s) / Z + epsilon`, the argument of the logarithm in the bias.
    fn density_ratio(&self, cv: &[f64]) -> f64 {
        self.kernel_density(cv) / self.z + self.epsilon()
    }
    /// Evaluate the (un-normalised) kernel density estimate at `cv`.
    ///
    /// Because the kernels are weighted by their reweighting factors this is
    /// proportional to the estimate of the unbiased distribution `P(s)`.
    pub fn kernel_density(&self, cv: &[f64]) -> f64 {
        self.kernels.iter().map(|k| k.evaluate(cv)).sum()
    }
    /// Compute the OPES bias potential at `cv`.
    ///
    /// `V(s) = (1 - 1/gamma) * k_B*T * ln(P_hat(s) / Z + epsilon)`, see
    /// [`OpesState::epsilon`] for the regularisation.  Zero before the first
    /// kernel is deposited and for `gamma <= 1`.
    pub fn bias_potential(&self, cv: &[f64]) -> f64 {
        if !self.is_biasing() {
            return 0.0;
        }
        let ratio = self.density_ratio(cv).max(f64::MIN_POSITIVE);
        self.bias_prefactor() * self.thermal_energy() * ratio.ln()
    }
    /// Reweighting factor `exp(V(cv) / k_B*T)` of a sample taken at `cv` in the
    /// current biased ensemble (1 while no bias is applied).
    pub fn reweighting_factor(&self, cv: &[f64]) -> f64 {
        if !self.is_biasing() {
            return 1.0;
        }
        // exp(beta * V) = (P_hat / Z + epsilon)^(1 - 1/gamma)
        let ratio = self.density_ratio(cv).max(f64::MIN_POSITIVE);
        ratio.powf(self.bias_prefactor())
    }
    /// Update the OPES state at the current CV value.
    ///
    /// Every `stride` steps a kernel centred at `cv` is deposited with the
    /// weight `exp(V(cv) / k_B*T)` of the bias that produced this sample, and
    /// the normalisation `Z` is updated.
    pub fn update(&mut self, cv: &[f64]) {
        if self.step.is_multiple_of(self.stride) {
            let weight = self.reweighting_factor(cv);
            self.push_weighted_kernel(cv, weight);
        }
        self.step += 1;
    }
    /// Append a kernel of the given weight and update `Z` incrementally.
    ///
    /// With `S = sum_k P_hat(c_k)` over the kernel centres, adding the kernel
    /// `n` changes `S` by `P_hat_old(c_n) + w_n + w_n * sum_k G_n(c_k)`
    /// (the new kernel's own value, the old density at the new centre and the
    /// new kernel's tail at the old centres), where `G_n` is the unit Gaussian.
    fn push_weighted_kernel(&mut self, cv: &[f64], weight: f64) {
        let kernel = GaussianHill::new(cv.to_vec(), weight, self.bandwidth.clone());
        let unit_kernel = GaussianHill::new(cv.to_vec(), 1.0, self.bandwidth.clone());
        let old_density_at_new_centre = self.kernel_density(cv);
        let new_kernel_at_old_centres: f64 = self
            .kernels
            .iter()
            .map(|k| unit_kernel.evaluate(&k.center))
            .sum();
        self.centre_density_sum +=
            old_density_at_new_centre + weight * (1.0 + new_kernel_at_old_centres);
        self.kernels.push(kernel);
        self.z = self.centre_density_sum / self.kernels.len() as f64;
    }
    /// Recompute `Z` from scratch (`O(n²)` in the number of kernels).
    ///
    /// Needed only after `kernels` was modified directly.
    pub fn recompute_normalisation(&mut self) {
        if self.kernels.is_empty() {
            self.centre_density_sum = 0.0;
            self.z = 0.0;
            return;
        }
        self.centre_density_sum = self
            .kernels
            .iter()
            .map(|k| self.kernel_density(&k.center))
            .sum();
        self.z = self.centre_density_sum / self.kernels.len() as f64;
    }
    /// Merge nearby kernels (see [`compress_opes_kernels`]) and refresh `Z`.
    ///
    /// Merged kernels carry the sum of the weights of the originals, so the
    /// density estimate is preserved up to the smoothing of the merge.
    pub fn compress(&mut self, compression_threshold: f64) {
        compress_opes_kernels(&mut self.kernels, compression_threshold);
        self.recompute_normalisation();
    }
    /// Free-energy estimate `F(s) = -k_B*T ln(P_hat(s) / P_hat_max)` of the
    /// *unbiased* system at the given CV points (kJ/mol, lowest density
    /// maximum shifted to zero; points without density are `+inf`).
    pub fn free_energy_estimate(&self, points: &[Vec<f64>]) -> Vec<f64> {
        let densities: Vec<f64> = points.iter().map(|s| self.kernel_density(s)).collect();
        let peak = densities.iter().copied().fold(0.0_f64, f64::max);
        if peak <= 0.0 {
            return vec![0.0; points.len()];
        }
        let kt = self.thermal_energy();
        densities
            .iter()
            .map(|&d| {
                if d > 0.0 {
                    -kt * (d / peak).ln()
                } else {
                    f64::INFINITY
                }
            })
            .collect()
    }
}
/// Coordination number CV: sum_j sigma(r_ij) where sigma is a switching function.
///
/// sigma(r) = (1 - (r/r_0)^n) / (1 - (r/r_0)^m)
///
/// This is the standard PLUMED-style coordination number.
pub struct CoordinationNumberCV {
    /// Reference atom index.
    pub atom_ref: usize,
    /// Indices of neighbour atoms.
    pub neighbours: Vec<usize>,
    /// Reference distance r_0.
    pub r_0: f64,
    /// Numerator exponent n (default: 6).
    pub n: u32,
    /// Denominator exponent m (default: 12).
    pub m: u32,
    /// Name.
    pub name: String,
}
impl CoordinationNumberCV {
    /// Create a new coordination number CV.
    pub fn new(
        atom_ref: usize,
        neighbours: Vec<usize>,
        r_0: f64,
        n: u32,
        m: u32,
        name: impl Into<String>,
    ) -> Self {
        Self {
            atom_ref,
            neighbours,
            r_0,
            n,
            m,
            name: name.into(),
        }
    }
    /// Evaluate the switching function sigma(r).
    ///
    /// The point `r = r_0` is a removable singularity (`0/0`) with the finite
    /// limit `sigma(r_0) = n / m`, which is returned there.
    pub fn switching_function(&self, r: f64) -> f64 {
        self.switching_value_and_derivative(r / self.r_0).0
    }
    /// Derivative `d sigma / d r` of the switching function.
    ///
    /// At the removable singularity `r = r_0` the limit
    /// `d sigma / d r = n (n - m) / (2 m r_0)` is returned.
    pub fn switching_derivative(&self, r: f64) -> f64 {
        self.switching_value_and_derivative(r / self.r_0).1 / self.r_0
    }
    /// `(sigma, d sigma / dx)` as functions of the reduced distance `x = r / r_0`.
    ///
    /// Away from `x = 1`
    ///
    /// ```text
    /// sigma  = (1 - x^n) / (1 - x^m)
    /// sigma' = [-n x^(n-1) (1 - x^m) + m x^(m-1) (1 - x^n)] / (1 - x^m)^2
    /// ```
    ///
    /// Both the numerator and the denominator vanish at `x = 1`, so within
    /// [`SWITCHING_SERIES_HALF_WIDTH`] of it the Taylor series of the ratio of the
    /// two geometric sums is used instead of the cancelling closed form.  The
    /// function is degenerate (and returns zero) for `n = 0` (numerator
    /// identically zero) and `m = 0`; for `n = m` it is identically one.
    fn switching_value_and_derivative(&self, x: f64) -> (f64, f64) {
        if self.n == 0 || self.m == 0 {
            return (0.0, 0.0);
        }
        if self.n == self.m {
            return (1.0, 0.0);
        }
        let n = f64::from(self.n);
        let m = f64::from(self.m);
        let eps = x - 1.0;
        if eps.abs() < SWITCHING_SERIES_HALF_WIDTH {
            let coefficients = removable_singularity_series(n, m);
            let prefactor = n / m;
            let mut value = 0.0;
            let mut slope = 0.0;
            for (order, c) in coefficients.iter().enumerate().rev() {
                value = value * eps + c;
                if order > 0 {
                    slope = slope * eps + c * order as f64;
                }
            }
            return (prefactor * value, prefactor * slope);
        }
        let xn = x.powi(self.n as i32);
        let xm = x.powi(self.m as i32);
        if !xn.is_finite() || !xm.is_finite() {
            // x^m overflowed: far beyond r_0 the switching function has decayed to 0.
            return (0.0, 0.0);
        }
        let xn_1 = x.powi(self.n as i32 - 1);
        let xm_1 = x.powi(self.m as i32 - 1);
        let denominator = 1.0 - xm;
        let value = (1.0 - xn) / denominator;
        let slope = (-n * xn_1 * denominator + m * xm_1 * (1.0 - xn)) / (denominator * denominator);
        (value, slope)
    }
}
/// Half width of the interval around `x = r / r_0 = 1` in which the switching
/// function is evaluated from its Taylor series.
pub(super) const SWITCHING_SERIES_HALF_WIDTH: f64 = 1.0e-4;
/// Taylor coefficients `c_0..c_4` in `eps = x - 1` of `m * sigma(x) / n`.
///
/// With `P_k(eps) = ((1 + eps)^k - 1) / (k eps) = sum_j C(k-1, j) eps^j / (j + 1)`
/// the switching function is `sigma = (n / m) P_n(eps) / P_m(eps)`; the ratio of
/// the two power series is formed by long division.  `c_0 = 1`,
/// `c_1 = (n - m) / 2`.
fn removable_singularity_series(n: f64, m: f64) -> [f64; 5] {
    let geometric_coefficient = |k: f64, order: usize| -> f64 {
        let mut binomial = 1.0;
        for i in 0..order {
            binomial *= (k - 1.0 - i as f64) / (i as f64 + 1.0);
        }
        binomial / (order as f64 + 1.0)
    };
    let mut c = [0.0_f64; 5];
    c[0] = 1.0;
    for order in 1..c.len() {
        let mut coefficient = geometric_coefficient(n, order);
        for i in 1..=order {
            coefficient -= geometric_coefficient(m, i) * c[order - i];
        }
        c[order] = coefficient;
    }
    c
}
/// CV measuring the end-to-end distance of a polymer chain.
///
/// Defined as the distance between the first and last atoms in the chain.
pub struct EndToEndDistanceCV {
    /// Index of the first atom.
    pub first: usize,
    /// Index of the last atom.
    pub last: usize,
    /// Name.
    pub name: String,
}
impl EndToEndDistanceCV {
    /// Create a new end-to-end distance CV.
    pub fn new(first: usize, last: usize, name: impl Into<String>) -> Self {
        Self {
            first,
            last,
            name: name.into(),
        }
    }
}
/// Parallel tempering metadynamics: multiple replicas at different temperatures,
/// each running well-tempered metadynamics with exchange attempts between replicas.
pub struct PtMetaD {
    /// Replica states (one per temperature).
    pub replicas: Vec<MetadynamicsState>,
    /// Temperatures for each replica (K).
    pub temperatures: Vec<f64>,
    /// Number of accepted exchanges.
    pub n_accepted: usize,
    /// Total exchange attempts.
    pub n_attempts: usize,
}
impl PtMetaD {
    /// Create a PTMetaD simulation with given temperatures.
    ///
    /// Each replica gets identical metadynamics parameters except temperature.
    pub fn new(
        temperatures: Vec<f64>,
        heights: Vec<f64>,
        widths_per_replica: Vec<Vec<f64>>,
        stride: usize,
        delta_t: f64,
    ) -> Self {
        assert_eq!(
            temperatures.len(),
            heights.len(),
            "need one height per replica"
        );
        assert_eq!(
            temperatures.len(),
            widths_per_replica.len(),
            "need widths per replica"
        );
        let replicas = temperatures
            .iter()
            .enumerate()
            .map(|(i, &temp)| {
                let params = MetadynamicsParams::new(
                    heights[i],
                    widths_per_replica[i].clone(),
                    stride,
                    temp,
                    delta_t,
                );
                MetadynamicsState::new(params)
            })
            .collect();
        Self {
            replicas,
            temperatures,
            n_accepted: 0,
            n_attempts: 0,
        }
    }
    /// Number of replicas.
    pub fn n_replicas(&self) -> usize {
        self.replicas.len()
    }
    /// Reduced cost `Δ` of exchanging the configurations of the adjacent
    /// replicas `i` and `j = i + 1`.
    ///
    /// Replica `k` samples the weight `exp(-β_k [U(x) + V_k(s(x))])`, where `U` is
    /// the potential energy, `V_k` the metadynamics bias of that replica and
    /// `s(x)` the collective variable.  Swapping the configuration `x_i` of
    /// replica `i` (CV value `cv_i`, potential energy `potential_energy_i`) with
    /// the configuration `x_j` of replica `j` (`cv_j`, `potential_energy_j`)
    /// changes the weight by
    ///
    /// ```text
    /// Δ = (β_i − β_j) (U_j − U_i) + β_i [V_i(s_j) − V_i(s_i)] + β_j [V_j(s_i) − V_j(s_j)]
    /// ```
    ///
    /// and the exchange has to be accepted with probability `min(1, exp(-Δ))`.
    /// Potential energies are in kJ/mol, like the hill heights.
    ///
    /// # Errors
    /// [`MdError::General`] if `i + 1` is not a valid replica, a replica
    /// temperature is not positive and finite, or `Δ` is NaN.
    pub fn exchange_delta(
        &self,
        i: usize,
        cv_i: &[f64],
        cv_j: &[f64],
        potential_energy_i: f64,
        potential_energy_j: f64,
    ) -> Result<f64, MdError> {
        let j = i.saturating_add(1);
        check_replica_index(j, self.replicas.len())?;
        let temperature_of = |k: usize| {
            self.temperatures
                .get(k)
                .copied()
                .ok_or_else(|| MdError::General(format!("no temperature stored for replica {k}")))
        };
        let beta_i = inverse_thermal_energy(temperature_of(i)?)?;
        let beta_j = inverse_thermal_energy(temperature_of(j)?)?;
        let v_i_at_i = self.replicas[i].bias_potential(cv_i);
        let v_i_at_j = self.replicas[i].bias_potential(cv_j);
        let v_j_at_i = self.replicas[j].bias_potential(cv_i);
        let v_j_at_j = self.replicas[j].bias_potential(cv_j);
        let delta = (beta_i - beta_j) * (potential_energy_j - potential_energy_i)
            + beta_i * (v_i_at_j - v_i_at_i)
            + beta_j * (v_j_at_i - v_j_at_j);
        if delta.is_nan() {
            return Err(MdError::General(
                "replica-exchange cost is NaN (check the potential energies and CV values)".into(),
            ));
        }
        Ok(delta)
    }
    /// Attempt a replica exchange between adjacent replicas `i` and `i + 1`.
    ///
    /// The cost `Δ` of [`PtMetaD::exchange_delta`] is evaluated and the exchange
    /// is accepted with the Metropolis probability `min(1, exp(-Δ))`, drawing one
    /// uniform number from `rng`.  The attempt (and acceptance) counters are
    /// updated; it is up to the caller to swap the configurations of the two
    /// replicas when `Ok(true)` is returned.
    ///
    /// # Errors
    /// See [`PtMetaD::exchange_delta`]; nothing is counted when an error occurs.
    pub fn try_exchange_adjacent<R: Rng + ?Sized>(
        &mut self,
        i: usize,
        cv_i: &[f64],
        cv_j: &[f64],
        potential_energy_i: f64,
        potential_energy_j: f64,
        rng: &mut R,
    ) -> Result<bool, MdError> {
        let delta = self.exchange_delta(i, cv_i, cv_j, potential_energy_i, potential_energy_j)?;
        let accept = metropolis_accept(delta, rng);
        self.n_attempts += 1;
        if accept {
            self.n_accepted += 1;
        }
        Ok(accept)
    }
    /// Acceptance rate.
    pub fn acceptance_rate(&self) -> f64 {
        if self.n_attempts == 0 {
            return 0.0;
        }
        self.n_accepted as f64 / self.n_attempts as f64
    }
    /// Total hills across all replicas.
    pub fn total_hills(&self) -> usize {
        self.replicas.iter().map(|r| r.hills.len()).sum()
    }
}
/// Local elevation (flooding) state: adds Gaussian repulsors to flatten the
/// energy landscape around known local minima.
///
/// Reference: Huber, Torda, van Gunsteren, J. Comput.-Aided Mol. Des. 8, 695 (1994).
pub struct LocalElevation {
    /// Deposited Gaussian repulsors (same structure as metadynamics hills).
    pub repulsors: Vec<GaussianHill>,
    /// Maximum allowed flooding height (prevents runaway bias).
    pub max_height: f64,
    /// Decay rate: heights decrease by factor (1-decay) per step.
    pub decay_rate: f64,
}
impl LocalElevation {
    /// Create a new local elevation state.
    pub fn new(max_height: f64, decay_rate: f64) -> Self {
        Self {
            repulsors: Vec::new(),
            max_height,
            decay_rate,
        }
    }
    /// Deposit a new repulsor at `cv_values` with `height` (clamped to max_height).
    pub fn deposit(&mut self, cv_values: &[f64], height: f64, widths: Vec<f64>) {
        let h = height.min(self.max_height);
        self.repulsors
            .push(GaussianHill::new(cv_values.to_vec(), h, widths));
    }
    /// Decay all repulsor heights by factor (1 - decay_rate).
    pub fn decay(&mut self) {
        let factor = 1.0 - self.decay_rate;
        for r in &mut self.repulsors {
            r.height *= factor;
        }
        self.repulsors.retain(|r| r.height > 1e-10);
    }
    /// Total flooding potential at `cv_values`.
    pub fn flooding_potential(&self, cv_values: &[f64]) -> f64 {
        self.repulsors.iter().map(|r| r.evaluate(cv_values)).sum()
    }
}
/// Parameters for a metadynamics simulation.
pub struct MetadynamicsParams {
    /// Initial hill height (kJ/mol).
    pub height: f64,
    /// Gaussian widths per CV dimension (one per CV).
    pub widths: Vec<f64>,
    /// Number of steps between hill depositions.
    pub deposition_stride: usize,
    /// Simulation temperature (K).
    pub temperature: f64,
    /// Well-tempered bias factor (0 = plain metadynamics, > 0 = well-tempered).
    pub delta_t: f64,
}
impl MetadynamicsParams {
    /// Create new metadynamics parameters.
    pub fn new(
        height: f64,
        widths: Vec<f64>,
        deposition_stride: usize,
        temperature: f64,
        delta_t: f64,
    ) -> Self {
        Self {
            height,
            widths,
            deposition_stride,
            temperature,
            delta_t,
        }
    }
}
/// Multiple-walker metadynamics: combines hills from several parallel walkers.
///
/// Each walker deposits hills independently; between deposition steps they share
/// all hills so that each walker feels the combined bias from the entire ensemble.
///
/// The struct keeps a shared pool of every hill that has been synchronised so
/// far together with, per walker, how many entries of its hill list are already
/// part of that pool.  A synchronisation therefore only moves the hills
/// deposited since the previous one, and repeated synchronisations without new
/// deposits are no-ops.
pub struct MultipleWalkersMetadynamics {
    /// Per-walker metadynamics states.
    pub walkers: Vec<MetadynamicsState>,
    /// Every hill deposited by any walker, in synchronisation order.
    pool: Vec<GaussianHill>,
    /// Per walker, the length of its hill list at the last synchronisation
    /// (the leading entries of `walkers[w].hills` that already equal the pool).
    synced_len: Vec<usize>,
}
impl MultipleWalkersMetadynamics {
    /// Create `n_walkers` independent metadynamics walkers with identical parameters.
    pub fn new(n_walkers: usize, params_fn: impl Fn() -> MetadynamicsParams) -> Self {
        let walkers = (0..n_walkers)
            .map(|_| MetadynamicsState::new(params_fn()))
            .collect();
        Self {
            walkers,
            pool: Vec::new(),
            synced_len: vec![0; n_walkers],
        }
    }
    /// Synchronise all walkers: share hills across the entire ensemble.
    ///
    /// The hills each walker deposited since the previous synchronisation are
    /// appended to the shared pool, and then every walker's hill list is set to
    /// the pool, so each walker holds exactly one copy of every hill ever
    /// deposited.  Calling this again without new deposits changes nothing (the
    /// bias is never counted twice).
    pub fn synchronise(&mut self) {
        // Walkers appended after construction have not contributed anything yet.
        self.synced_len.resize(self.walkers.len(), 0);
        for (walker, &synced) in self.walkers.iter().zip(self.synced_len.iter()) {
            let first_new = synced.min(walker.hills.len());
            self.pool.extend(walker.hills[first_new..].iter().cloned());
        }
        for (walker, synced) in self.walkers.iter_mut().zip(self.synced_len.iter_mut()) {
            walker.hills = self.pool.clone();
            *synced = self.pool.len();
        }
    }
    /// The shared pool: every hill that has been synchronised so far, once.
    pub fn shared_hills(&self) -> &[GaussianHill] {
        &self.pool
    }
    /// Total number of hills held by all walkers.
    ///
    /// Right after a synchronisation every walker holds the whole pool, so this
    /// equals `n_walkers() * shared_hills().len()`.
    pub fn total_hills(&self) -> usize {
        self.walkers.iter().map(|w| w.hills.len()).sum()
    }
    /// Number of walkers.
    pub fn n_walkers(&self) -> usize {
        self.walkers.len()
    }
}
/// CV measuring the (unweighted) radius of gyration for a subset of atoms.
pub struct RadiusOfGyrationCV {
    /// Indices of atoms to include in the Rg calculation.
    pub atom_indices: Vec<usize>,
    /// Name of this CV.
    pub name: String,
}
impl RadiusOfGyrationCV {
    /// Create a new radius of gyration CV.
    pub fn new(atom_indices: Vec<usize>, name: impl Into<String>) -> Self {
        Self {
            atom_indices,
            name: name.into(),
        }
    }
}
/// A Gaussian hill deposited in metadynamics.
#[derive(Debug, Clone)]
pub struct GaussianHill {
    /// Center in CV space (one value per CV dimension).
    pub center: Vec<f64>,
    /// Height of the Gaussian (kJ/mol).
    pub height: f64,
    /// Width (sigma) per CV dimension.
    pub widths: Vec<f64>,
}
impl GaussianHill {
    /// Create a new Gaussian hill.
    pub fn new(center: Vec<f64>, height: f64, widths: Vec<f64>) -> Self {
        Self {
            center,
            height,
            widths,
        }
    }
    /// Evaluate the Gaussian at `cv_values`:
    /// V = height * exp(-sum_i (s_i - center_i)^2 / (2*width_i^2))
    pub fn evaluate(&self, cv_values: &[f64]) -> f64 {
        let exponent: f64 = cv_values
            .iter()
            .zip(self.center.iter())
            .zip(self.widths.iter())
            .map(|((&s, &c), &w)| {
                let d = s - c;
                d * d / (2.0 * w * w)
            })
            .sum();
        self.height * (-exponent).exp()
    }
    /// Gradient of the Gaussian with respect to CV values:
    /// dV/ds_i = -V * (s_i - center_i) / width_i^2
    pub fn gradient_wrt_cv(&self, cv_values: &[f64]) -> Vec<f64> {
        let v = self.evaluate(cv_values);
        cv_values
            .iter()
            .zip(self.center.iter())
            .zip(self.widths.iter())
            .map(|((&s, &c), &w)| -v * (s - c) / (w * w))
            .collect()
    }
}
/// Steered molecular dynamics with harmonic bias on a CV.
pub struct SteeredMD {
    /// Index of the CV to steer.
    pub cv_idx: usize,
    /// Initial target value of the CV.
    pub target_value: f64,
    /// Spring constant k (kJ/mol/unit^2).
    pub spring_constant: f64,
    /// Pulling velocity v (unit/step).
    pub pull_velocity: f64,
}
impl SteeredMD {
    /// Create a new steered MD protocol.
    pub fn new(cv_idx: usize, target_value: f64, spring_constant: f64, pull_velocity: f64) -> Self {
        Self {
            cv_idx,
            target_value,
            spring_constant,
            pull_velocity,
        }
    }
    /// Compute the harmonic bias potential at the current step.
    ///
    /// The moving target is: target(step) = target_value + pull_velocity * step
    /// U = 0.5 * k * (cv - target)^2
    pub fn harmonic_bias(&self, cv_value: f64, step: usize) -> f64 {
        let target = self.target_value + self.pull_velocity * step as f64;
        let delta = cv_value - target;
        0.5 * self.spring_constant * delta * delta
    }
    /// Compute Jarzynski-style work along a CV trajectory.
    ///
    /// W = sum_t F(t) * delta_cv(t)
    /// where F(t) = -dU/dcv = -k*(cv(t) - target(t)) and delta_cv = cv(t+1) - cv(t).
    pub fn work_done(&self, cv_trajectory: &[f64], initial: f64) -> f64 {
        if cv_trajectory.len() < 2 {
            return 0.0;
        }
        let mut work = 0.0;
        let mut target = initial;
        for t in 0..cv_trajectory.len() - 1 {
            let cv = cv_trajectory[t];
            let force = -self.spring_constant * (cv - target);
            let dcv = cv_trajectory[t + 1] - cv_trajectory[t];
            work += force * dcv;
            target += self.pull_velocity;
        }
        work
    }
}
/// Bias-exchange metadynamics (BE-META): multiple replicas each biased along
/// a different CV, with periodic exchange between replicas.
///
/// Each replica carries an independent metadynamics state for its own bias.
/// Exchanges are attempted between replicas with the Metropolis criterion,
/// driven by a seeded random number generator (see
/// [`BiasExchangeMetadynamics::with_seed`]) or by a generator supplied by the
/// caller (see [`BiasExchangeMetadynamics::try_exchange_with_rng`]).
pub struct BiasExchangeMetadynamics {
    /// Per-replica metadynamics states.
    pub replicas: Vec<MetadynamicsState>,
    /// Temperature of all replicas (K).
    pub temperature: f64,
    /// Number of accepted exchanges.
    pub n_accepted: usize,
    /// Total exchange attempts.
    pub n_attempts: usize,
    /// Random stream behind [`BiasExchangeMetadynamics::try_exchange`].
    rng: StdRng,
}
impl BiasExchangeMetadynamics {
    /// Seed of the exchange random stream used by [`BiasExchangeMetadynamics::new`].
    pub const DEFAULT_SEED: u64 = 0x00BE_5EED_0001_0001;
    /// Create `n_replicas` replicas with given parameters.
    ///
    /// The exchange random stream is seeded with
    /// [`BiasExchangeMetadynamics::DEFAULT_SEED`], so runs are reproducible;
    /// use [`BiasExchangeMetadynamics::with_seed`] for independent streams.
    pub fn new(
        n_replicas: usize,
        temperature: f64,
        params_fn: impl Fn(usize) -> MetadynamicsParams,
    ) -> Self {
        Self::with_seed(n_replicas, temperature, params_fn, Self::DEFAULT_SEED)
    }
    /// Like [`BiasExchangeMetadynamics::new`] with an explicit seed for the
    /// random stream that decides the exchange attempts.
    pub fn with_seed(
        n_replicas: usize,
        temperature: f64,
        params_fn: impl Fn(usize) -> MetadynamicsParams,
        seed: u64,
    ) -> Self {
        let replicas = (0..n_replicas)
            .map(|i| MetadynamicsState::new(params_fn(i)))
            .collect();
        Self {
            replicas,
            temperature,
            n_accepted: 0,
            n_attempts: 0,
            rng: StdRng::seed_from_u64(seed),
        }
    }
    /// Reduced cost `Δ` of swapping the configurations of replicas `i` and `j`.
    ///
    /// With `V_k` the bias of replica `k`, `cv_i` the collective variable of the
    /// configuration currently held by replica `i` and `cv_j` that of replica
    /// `j`,
    ///
    /// ```text
    /// Δ = β [ V_i(cv_j) + V_j(cv_i) − V_i(cv_i) − V_j(cv_j) ]
    /// ```
    ///
    /// (all replicas share the temperature, so the potential energies cancel).
    /// Both CV values are taken in the CV space seen by both replica biases; a
    /// bias is evaluated at each of them.  The exchange has to be accepted with
    /// probability `min(1, exp(-Δ))`.
    ///
    /// # Errors
    /// [`MdError::General`] if `i` or `j` is not a replica, if `i == j`, if the
    /// temperature is not positive and finite, or if `Δ` is NaN.
    pub fn exchange_delta(
        &self,
        i: usize,
        j: usize,
        cv_i: &[f64],
        cv_j: &[f64],
    ) -> Result<f64, MdError> {
        check_replica_index(i, self.replicas.len())?;
        check_replica_index(j, self.replicas.len())?;
        if i == j {
            return Err(MdError::General(format!(
                "cannot exchange replica {i} with itself"
            )));
        }
        let beta = inverse_thermal_energy(self.temperature)?;
        let v_i_at_i = self.replicas[i].bias_potential(cv_i);
        let v_i_at_j = self.replicas[i].bias_potential(cv_j);
        let v_j_at_j = self.replicas[j].bias_potential(cv_j);
        let v_j_at_i = self.replicas[j].bias_potential(cv_i);
        let delta = beta * ((v_i_at_j + v_j_at_i) - (v_i_at_i + v_j_at_j));
        if delta.is_nan() {
            return Err(MdError::General(
                "bias-exchange cost is NaN (check the CV values)".into(),
            ));
        }
        Ok(delta)
    }
    /// Attempt a replica exchange between replicas `i` and `j`.
    ///
    /// Evaluates [`BiasExchangeMetadynamics::exchange_delta`] and accepts with
    /// the Metropolis probability `min(1, exp(-Δ))`, drawing one number from the
    /// internal seeded stream.  Returns `Ok(true)` if the exchange was accepted;
    /// the caller then swaps the configurations.
    ///
    /// # Errors
    /// See [`BiasExchangeMetadynamics::exchange_delta`]; nothing is counted when
    /// an error occurs.
    pub fn try_exchange(
        &mut self,
        i: usize,
        j: usize,
        cv_i: &[f64],
        cv_j: &[f64],
    ) -> Result<bool, MdError> {
        let delta = self.exchange_delta(i, j, cv_i, cv_j)?;
        let accept = metropolis_accept(delta, &mut self.rng);
        self.record_attempt(accept);
        Ok(accept)
    }
    /// Like [`BiasExchangeMetadynamics::try_exchange`], but the acceptance draw
    /// comes from the caller's generator `rng` (the internal stream is untouched).
    pub fn try_exchange_with_rng<R: Rng + ?Sized>(
        &mut self,
        i: usize,
        j: usize,
        cv_i: &[f64],
        cv_j: &[f64],
        rng: &mut R,
    ) -> Result<bool, MdError> {
        let delta = self.exchange_delta(i, j, cv_i, cv_j)?;
        let accept = metropolis_accept(delta, rng);
        self.record_attempt(accept);
        Ok(accept)
    }
    fn record_attempt(&mut self, accepted: bool) {
        self.n_attempts += 1;
        if accepted {
            self.n_accepted += 1;
        }
    }
    /// Overall acceptance rate.
    pub fn acceptance_rate(&self) -> f64 {
        if self.n_attempts == 0 {
            return 0.0;
        }
        self.n_accepted as f64 / self.n_attempts as f64
    }
}
/// Adaptive Gaussian hill whose width is determined locally from CV fluctuations.
///
/// The width adapts so that hills are deposited with widths proportional to
/// the short-time diffusion coefficient of the CV, preventing over-filling of
/// visited regions.
#[derive(Debug, Clone)]
pub struct AdaptiveGaussianHill {
    /// Hill center in CV space.
    pub center: Vec<f64>,
    /// Hill height (kJ/mol).
    pub height: f64,
    /// Adaptive width per dimension, estimated from CV fluctuations.
    pub widths: Vec<f64>,
    /// Step at which this hill was deposited.
    pub step: usize,
}
impl AdaptiveGaussianHill {
    /// Create a new adaptive hill.
    pub fn new(center: Vec<f64>, height: f64, widths: Vec<f64>, step: usize) -> Self {
        Self {
            center,
            height,
            widths,
            step,
        }
    }
    /// Evaluate the Gaussian at `cv_values`.
    pub fn evaluate(&self, cv_values: &[f64]) -> f64 {
        let exponent: f64 = cv_values
            .iter()
            .zip(self.center.iter())
            .zip(self.widths.iter())
            .map(|((&s, &c), &w)| {
                let w2 = (w * w).max(1e-30);
                let d = s - c;
                d * d / (2.0 * w2)
            })
            .sum();
        self.height * (-exponent).exp()
    }
}
/// CV measuring the distance between two atoms.
pub struct DistanceCV {
    /// Index of the first atom.
    pub atom_i: usize,
    /// Index of the second atom.
    pub atom_j: usize,
    /// Name of this CV.
    pub name: String,
}
impl DistanceCV {
    /// Create a new distance CV between atoms `i` and `j`.
    pub fn new(atom_i: usize, atom_j: usize, name: impl Into<String>) -> Self {
        Self {
            atom_i,
            atom_j,
            name: name.into(),
        }
    }
}
/// CV measuring the angle at atom_j formed by atoms i-j-k (in radians).
pub struct AngleCV {
    /// Index of the first atom.
    pub atom_i: usize,
    /// Index of the vertex atom (angle measured here).
    pub atom_j: usize,
    /// Index of the third atom.
    pub atom_k: usize,
    /// Name of this CV.
    pub name: String,
}
impl AngleCV {
    /// Create a new angle CV for atoms i-j-k.
    pub fn new(atom_i: usize, atom_j: usize, atom_k: usize, name: impl Into<String>) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            name: name.into(),
        }
    }
}
/// Funnel metadynamics configurator: builds a `FunnelMetadynamics` from
/// standard binding/unbinding simulation parameters.
pub struct FunnelConfig {
    /// Radius of the funnel entrance (Å or nm, consistent with CV units).
    pub entrance_radius: f64,
    /// Half-angle of the funnel cone (radians).
    pub half_angle: f64,
    /// Height of the funnel (extent along binding axis, same units as CV).
    pub height: f64,
    /// Harmonic wall constant k (kJ/mol/unit²).
    pub wall_constant: f64,
}
impl FunnelConfig {
    /// Create a funnel configuration.
    pub fn new(entrance_radius: f64, half_angle: f64, height: f64, wall_constant: f64) -> Self {
        Self {
            entrance_radius,
            half_angle,
            height,
            wall_constant,
        }
    }
    /// Maximum allowed radius at axial position `z` along the funnel.
    pub fn max_radius(&self, z: f64) -> f64 {
        if z < 0.0 {
            return self.entrance_radius;
        }
        if z > self.height {
            return 0.01;
        }
        (self.entrance_radius + z * self.half_angle.tan()).max(0.01)
    }
    /// Harmonic wall potential energy for a ligand at (z, r_perp).
    pub fn wall_potential(&self, z: f64, r_perp: f64) -> f64 {
        let r_max = self.max_radius(z);
        if r_perp > r_max {
            let delta = r_perp - r_max;
            0.5 * self.wall_constant * delta * delta
        } else {
            0.0
        }
    }
    /// Volume correction factor for the funnel (approximate).
    pub fn volume_correction(&self) -> f64 {
        std::f64::consts::PI * self.entrance_radius * self.entrance_radius * self.height / 3.0
    }
}
/// Funnel-restrained metadynamics for binding free energy calculations.
pub struct FunnelMetadynamics {
    /// Underlying metadynamics state.
    pub state: MetadynamicsState,
    /// Width of the funnel (Angstroms or nm, consistent with CV units).
    pub funnel_width: f64,
    /// Slope of the funnel (restraint stiffens along the funnel axis).
    pub funnel_slope: f64,
}
impl FunnelMetadynamics {
    /// Create a new funnel metadynamics simulation.
    pub fn new(state: MetadynamicsState, funnel_width: f64, funnel_slope: f64) -> Self {
        Self {
            state,
            funnel_width,
            funnel_slope,
        }
    }
    /// Compute the funnel potential energy as a harmonic restraint.
    ///
    /// The allowed radial distance at a given `projection` along the axis is:
    ///   r_max = funnel_width + funnel_slope * projection
    ///
    /// If `dist_from_axis > r_max`, a harmonic penalty is applied:
    ///   U = 0.5 * k * (dist_from_axis - r_max)^2
    ///
    /// where k is derived from the slope (units: kJ/mol/unit^2).
    pub fn funnel_potential(&self, dist_from_axis: f64, projection: f64) -> f64 {
        let r_max = self.funnel_width + self.funnel_slope * projection.abs();
        if dist_from_axis > r_max {
            let delta = dist_from_axis - r_max;
            0.5 * self.funnel_slope * delta * delta
        } else {
            0.0
        }
    }
}
/// CV time-series monitor that tracks values and computes running statistics.
pub struct CvMonitor {
    /// Name of the CV being monitored.
    pub name: String,
    /// Recorded CV values over time.
    pub values: Vec<f64>,
    /// Time points corresponding to each CV value.
    pub times: Vec<f64>,
}
impl CvMonitor {
    /// Create a new CV monitor.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            values: Vec::new(),
            times: Vec::new(),
        }
    }
    /// Record a CV value at a given time.
    pub fn record(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }
    /// Running mean of the CV.
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }
    /// Running standard deviation of the CV.
    pub fn std_dev(&self) -> f64 {
        let n = self.values.len();
        if n < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var = self.values.iter().map(|&v| (v - m).powi(2)).sum::<f64>() / (n - 1) as f64;
        var.sqrt()
    }
    /// Returns the minimum and maximum CV values observed.
    pub fn range(&self) -> (f64, f64) {
        if self.values.is_empty() {
            return (0.0, 0.0);
        }
        let lo = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    }
    /// Number of recorded data points.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Returns true if no data has been recorded.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}
/// CV measuring the dihedral angle for atoms i-j-k-l (in radians).
///
/// Gradient is computed via finite differences for simplicity.
pub struct DihedralCV {
    /// Index of atom i.
    pub i: usize,
    /// Index of atom j.
    pub j: usize,
    /// Index of atom k.
    pub k: usize,
    /// Index of atom l.
    pub l: usize,
    /// Name of this CV.
    pub name: String,
}
impl DihedralCV {
    /// Create a new dihedral CV for atoms i-j-k-l.
    pub fn new(i: usize, j: usize, k: usize, l: usize, name: impl Into<String>) -> Self {
        Self {
            i,
            j,
            k,
            l,
            name: name.into(),
        }
    }
    pub(super) fn compute_dihedral(
        positions: &[[f64; 3]],
        i: usize,
        j: usize,
        k: usize,
        l: usize,
    ) -> f64 {
        let ri = positions[i];
        let rj = positions[j];
        let rk = positions[k];
        let rl = positions[l];
        let b1 = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
        let b2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
        let b3 = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = cross(n1, b2);
        let n1_norm = norm3(n1);
        let n2_norm = norm3(n2);
        let m1_norm = norm3(m1);
        if n1_norm < 1e-15 || n2_norm < 1e-15 || m1_norm < 1e-15 {
            return 0.0;
        }
        let x = dot3(n1, n2) / (n1_norm * n2_norm);
        let y = dot3(m1, n2) / (m1_norm * n2_norm);
        x.clamp(-1.0, 1.0).acos().copysign(y)
    }
}
