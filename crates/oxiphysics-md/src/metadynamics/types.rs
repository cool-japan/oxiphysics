//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

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
/// OPES state: maintains a running estimate of the probability distribution
/// p(s) of the CV and computes a bias V(s) = (1 - 1/gamma) * k_B*T * log(p(s)/p_0).
///
/// Reference: Invernizzi & Parrinello, J. Phys. Chem. Lett. 11, 2731 (2020).
#[allow(dead_code)]
pub struct OpesState {
    /// Gaussian kernels forming the density estimate.
    pub kernels: Vec<GaussianHill>,
    /// Normalisation constant Z (running sum of kernel weights).
    pub z: f64,
    /// Bias factor gamma (> 1; use gamma=∞ for unbiased target).
    pub gamma: f64,
    /// Simulation temperature (K).
    pub temperature: f64,
    /// Step counter.
    pub step: usize,
    /// Deposition stride.
    pub stride: usize,
    /// Initial kernel height (determines bandwidth in probability space).
    pub bandwidth: Vec<f64>,
}
impl OpesState {
    /// Create a new OPES state.
    #[allow(dead_code)]
    pub fn new(gamma: f64, temperature: f64, stride: usize, bandwidth: Vec<f64>) -> Self {
        Self {
            kernels: Vec::new(),
            z: 0.0,
            gamma,
            temperature,
            step: 0,
            stride,
            bandwidth,
        }
    }
    /// Evaluate the (un-normalised) kernel density estimate at `cv`.
    #[allow(dead_code)]
    pub fn kernel_density(&self, cv: &[f64]) -> f64 {
        self.kernels.iter().map(|k| k.evaluate(cv)).sum()
    }
    /// Compute the OPES bias potential at `cv`.
    ///
    /// V(s) = (1 - 1/gamma) * k_B*T * ln(p_hat(s) / Z + epsilon)
    /// where epsilon prevents log(0).
    #[allow(dead_code)]
    pub fn bias_potential(&self, cv: &[f64]) -> f64 {
        if self.z < 1e-30 || self.gamma <= 1.0 {
            return 0.0;
        }
        let p_hat = self.kernel_density(cv);
        let factor = 1.0 - 1.0 / self.gamma;
        let log_arg = (p_hat / self.z).max(1e-300);
        factor * KB * self.temperature * log_arg.ln()
    }
    /// Update the OPES state at the current CV value.
    ///
    /// Adds a new kernel every `stride` steps and updates the normalisation Z.
    #[allow(dead_code)]
    pub fn update(&mut self, cv: &[f64]) {
        if self.step.is_multiple_of(self.stride) {
            let kernel = GaussianHill::new(cv.to_vec(), 1.0, self.bandwidth.clone());
            self.z += 1.0;
            self.kernels.push(kernel);
        }
        self.step += 1;
    }
}
/// Coordination number CV: sum_j sigma(r_ij) where sigma is a switching function.
///
/// sigma(r) = (1 - (r/r_0)^n) / (1 - (r/r_0)^m)
///
/// This is the standard PLUMED-style coordination number.
#[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn switching_function(&self, r: f64) -> f64 {
        let x = r / self.r_0;
        let xn = x.powi(self.n as i32);
        let xm = x.powi(self.m as i32);
        if (1.0 - xm).abs() < 1e-15 {
            return 0.0;
        }
        (1.0 - xn) / (1.0 - xm)
    }
}
/// CV measuring the end-to-end distance of a polymer chain.
///
/// Defined as the distance between the first and last atoms in the chain.
#[allow(dead_code)]
pub struct EndToEndDistanceCV {
    /// Index of the first atom.
    pub first: usize,
    /// Index of the last atom.
    pub last: usize,
    /// Name.
    pub name: String,
}
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
    /// Attempt a replica exchange between adjacent replicas `i` and `i+1`.
    ///
    /// Exchange probability:
    ///   P = min(1, exp(Δ))
    ///   Δ = (beta_i - beta_{i+1}) * (V_i(s_i) - V_i(s_{i+1}))
    ///
    /// Returns true if exchange was accepted (uses deterministic criterion here).
    pub fn try_exchange_adjacent(&mut self, i: usize, cv_i: &[f64], cv_j: &[f64]) -> bool {
        assert!(i + 1 < self.replicas.len(), "replica index out of range");
        self.n_attempts += 1;
        let beta_i = 1.0 / (KB * self.temperatures[i]);
        let beta_j = 1.0 / (KB * self.temperatures[i + 1]);
        let v_i_at_i = self.replicas[i].bias_potential(cv_i);
        let v_i_at_j = self.replicas[i].bias_potential(cv_j);
        let delta = (beta_j - beta_i) * (v_i_at_i - v_i_at_j);
        let accept = delta >= 0.0 || (-delta).exp() > 0.5;
        if accept {
            self.n_accepted += 1;
        }
        accept
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(max_height: f64, decay_rate: f64) -> Self {
        Self {
            repulsors: Vec::new(),
            max_height,
            decay_rate,
        }
    }
    /// Deposit a new repulsor at `cv_values` with `height` (clamped to max_height).
    #[allow(dead_code)]
    pub fn deposit(&mut self, cv_values: &[f64], height: f64, widths: Vec<f64>) {
        let h = height.min(self.max_height);
        self.repulsors
            .push(GaussianHill::new(cv_values.to_vec(), h, widths));
    }
    /// Decay all repulsor heights by factor (1 - decay_rate).
    #[allow(dead_code)]
    pub fn decay(&mut self) {
        let factor = 1.0 - self.decay_rate;
        for r in &mut self.repulsors {
            r.height *= factor;
        }
        self.repulsors.retain(|r| r.height > 1e-10);
    }
    /// Total flooding potential at `cv_values`.
    #[allow(dead_code)]
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
pub struct MultipleWalkersMetadynamics {
    /// Per-walker metadynamics states.
    pub walkers: Vec<MetadynamicsState>,
}
impl MultipleWalkersMetadynamics {
    /// Create `n_walkers` independent metadynamics walkers with identical parameters.
    pub fn new(n_walkers: usize, params_fn: impl Fn() -> MetadynamicsParams) -> Self {
        let walkers = (0..n_walkers)
            .map(|_| MetadynamicsState::new(params_fn()))
            .collect();
        Self { walkers }
    }
    /// Synchronise all walkers: share hills across the entire ensemble.
    ///
    /// Each walker receives a copy of every hill deposited by any other walker.
    /// This merges hills from all walkers into a single shared pool for all.
    pub fn synchronise(&mut self) {
        let all_hills: Vec<GaussianHill> = self
            .walkers
            .iter()
            .flat_map(|w| {
                w.hills
                    .iter()
                    .map(|h| GaussianHill::new(h.center.clone(), h.height, h.widths.clone()))
            })
            .collect();
        for walker in &mut self.walkers {
            walker.hills.clear();
            for h in &all_hills {
                walker.hills.push(GaussianHill::new(
                    h.center.clone(),
                    h.height,
                    h.widths.clone(),
                ));
            }
        }
    }
    /// Total number of deposited hills across all walkers (before synchronisation).
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
/// Each replica carries an independent metadynamics state for its own CV.
/// Exchanges are attempted between adjacent replicas using Metropolis criterion.
#[allow(dead_code)]
pub struct BiasExchangeMetadynamics {
    /// Per-replica metadynamics states.
    pub replicas: Vec<MetadynamicsState>,
    /// Temperature of all replicas (K).
    pub temperature: f64,
    /// Number of accepted exchanges.
    pub n_accepted: usize,
    /// Total exchange attempts.
    pub n_attempts: usize,
}
impl BiasExchangeMetadynamics {
    /// Create `n_replicas` replicas with given parameters.
    #[allow(dead_code)]
    pub fn new(
        n_replicas: usize,
        temperature: f64,
        params_fn: impl Fn(usize) -> MetadynamicsParams,
    ) -> Self {
        let replicas = (0..n_replicas)
            .map(|i| MetadynamicsState::new(params_fn(i)))
            .collect();
        Self {
            replicas,
            temperature,
            n_accepted: 0,
            n_attempts: 0,
        }
    }
    /// Attempt a replica exchange between replicas `i` and `j`.
    ///
    /// Uses the Metropolis criterion:
    /// Δ = β * (V_bias_j(cv_i) - V_bias_j(cv_j) + V_bias_i(cv_j) - V_bias_i(cv_i))
    ///
    /// Returns `true` if the exchange was accepted.
    #[allow(dead_code)]
    pub fn try_exchange(&mut self, i: usize, j: usize, cv_i: &[f64], cv_j: &[f64]) -> bool {
        self.n_attempts += 1;
        let beta = 1.0 / (KB * self.temperature);
        let v_i_at_i = self.replicas[i].bias_potential(cv_i);
        let v_i_at_j = self.replicas[i].bias_potential(cv_j);
        let v_j_at_j = self.replicas[j].bias_potential(cv_j);
        let v_j_at_i = self.replicas[j].bias_potential(cv_i);
        let delta = beta * ((v_i_at_j + v_j_at_i) - (v_i_at_i + v_j_at_j));
        let accept = if delta <= 0.0 {
            true
        } else {
            (-delta).exp() > 0.5
        };
        if accept {
            self.n_accepted += 1;
        }
        accept
    }
    /// Overall acceptance rate.
    #[allow(dead_code)]
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(center: Vec<f64>, height: f64, widths: Vec<f64>, step: usize) -> Self {
        Self {
            center,
            height,
            widths,
            step,
        }
    }
    /// Evaluate the Gaussian at `cv_values`.
    #[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            values: Vec::new(),
            times: Vec::new(),
        }
    }
    /// Record a CV value at a given time.
    #[allow(dead_code)]
    pub fn record(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }
    /// Running mean of the CV.
    #[allow(dead_code)]
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }
    /// Running standard deviation of the CV.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Returns true if no data has been recorded.
    #[allow(dead_code)]
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
