#![allow(clippy::too_many_arguments)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Replica Exchange Molecular Dynamics (REMD / Parallel Tempering).
//!
//! REMD runs multiple independent MD replicas at different temperatures and
//! periodically attempts to swap configurations between adjacent replicas
//! using a Metropolis criterion.  This greatly accelerates the exploration of
//! rugged free-energy landscapes.
//!
//! # Metropolis exchange criterion
//!
//! Replica `i` (inverse temperature β_i = 1/(k_B T_i)) and replica `j`
//! (β_j) attempt a swap with acceptance probability:
//!
//! ```text
//! P_acc = min(1, exp((β_i − β_j)(E_i − E_j)))
//! ```

/// State of one REMD replica.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RemdReplica {
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Current potential energy (kJ mol⁻¹).
    pub potential_energy: f64,
    /// Unique replica identifier.
    pub replica_id: usize,
    /// Number of accepted exchanges this replica has participated in.
    pub exchange_count: u32,
}

impl RemdReplica {
    /// Create a new [`RemdReplica`].
    pub fn new(temperature: f64, potential_energy: f64, replica_id: usize) -> Self {
        Self {
            temperature,
            potential_energy,
            replica_id,
            exchange_count: 0,
        }
    }
}

// ─── RemdConfig ──────────────────────────────────────────────────────────────

/// Configuration for a REMD run: temperature ladder and Boltzmann constant.
#[derive(Debug, Clone)]
pub struct RemdConfig {
    /// Temperature ladder (K).
    pub temperatures: Vec<f64>,
    /// Boltzmann constant (kJ mol⁻¹ K⁻¹).
    pub kb: f64,
}

impl RemdConfig {
    /// Build a geometrically spaced temperature ladder.
    ///
    /// # Arguments
    /// * `t_min`      – Lowest temperature (K).
    /// * `t_max`      – Highest temperature (K).
    /// * `n_replicas` – Total number of replicas (including endpoints).
    pub fn new(t_min: f64, t_max: f64, n_replicas: usize) -> Self {
        Self {
            temperatures: Self::geometric_temperatures(t_min, t_max, n_replicas),
            kb: 8.314e-3,
        }
    }

    /// Compute geometrically spaced temperatures between `t_min` and `t_max`.
    pub fn geometric_temperatures(t_min: f64, t_max: f64, n: usize) -> Vec<f64> {
        if n == 1 {
            return vec![t_min];
        }
        let ratio = (t_max / t_min).powf(1.0 / (n - 1) as f64);
        (0..n).map(|i| t_min * ratio.powi(i as i32)).collect()
    }
}

// ─── metropolis_exchange_probability ─────────────────────────────────────────

/// Metropolis exchange probability between two replicas.
///
/// Returns `min(1, exp((β_i − β_j)(E_i − E_j)))`.
///
/// # Arguments
/// * `e_i`    – Potential energy of replica i.
/// * `e_j`    – Potential energy of replica j.
/// * `beta_i` – Inverse temperature of replica i (1 / (k_B T_i)).
/// * `beta_j` – Inverse temperature of replica j.
pub fn metropolis_exchange_probability(e_i: f64, e_j: f64, beta_i: f64, beta_j: f64) -> f64 {
    let delta = (beta_i - beta_j) * (e_i - e_j);
    if delta >= 0.0 { 1.0 } else { delta.exp() }
}

// ─── attempt_exchange ────────────────────────────────────────────────────────

/// Attempt a configuration exchange between two replicas.
///
/// If `accept` is true, the temperatures of the two replicas are swapped and
/// both `exchange_count` fields are incremented.
///
/// # Arguments
/// * `replica_i` – First replica (mutably borrowed).
/// * `replica_j` – Second replica (mutably borrowed).
/// * `kb`        – Boltzmann constant (kJ mol⁻¹ K⁻¹).
/// * `accept`    – Whether to accept the exchange.
pub fn attempt_exchange(
    replica_i: &mut RemdReplica,
    replica_j: &mut RemdReplica,
    _kb: f64,
    accept: bool,
) {
    if accept {
        std::mem::swap(&mut replica_i.temperature, &mut replica_j.temperature);
        replica_i.exchange_count += 1;
        replica_j.exchange_count += 1;
    }
}

// ─── RemdExchangeLog ─────────────────────────────────────────────────────────

/// Log of all exchange attempts for post-analysis.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct RemdExchangeLog {
    /// Each entry is `(replica_i_id, replica_j_id, probability, accepted)`.
    pub exchanges: Vec<(usize, usize, f64, bool)>,
}

impl RemdExchangeLog {
    /// Create an empty log.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record one exchange attempt.
    pub fn record(&mut self, i: usize, j: usize, prob: f64, accepted: bool) {
        self.exchanges.push((i, j, prob, accepted));
    }

    /// Overall acceptance rate across all logged attempts.
    ///
    /// Returns `0.0` if no attempts have been logged.
    pub fn acceptance_rate(&self) -> f64 {
        if self.exchanges.is_empty() {
            return 0.0;
        }
        let accepted = self.exchanges.iter().filter(|e| e.3).count();
        accepted as f64 / self.exchanges.len() as f64
    }

    /// Acceptance rate for the specific pair `(i, j)`.
    ///
    /// Returns `0.0` if no attempts for that pair exist.
    pub fn acceptance_rate_for_pair(&self, i: usize, j: usize) -> f64 {
        let pair_entries: Vec<_> = self
            .exchanges
            .iter()
            .filter(|e| (e.0 == i && e.1 == j) || (e.0 == j && e.1 == i))
            .collect();
        if pair_entries.is_empty() {
            return 0.0;
        }
        let accepted = pair_entries.iter().filter(|e| e.3).count();
        accepted as f64 / pair_entries.len() as f64
    }
}

// ─── optimal_exchange_spacing ─────────────────────────────────────────────────

/// Estimate the optimal fractional temperature spacing `ΔT/T` for a target
/// exchange acceptance rate.
///
/// Simplified formula: `ΔT/T ≈ sqrt(2 / (n_dof * Cv)) * f_target`
///
/// # Arguments
/// * `target_acceptance` – Desired acceptance rate (e.g. 0.2 for 20 %).
/// * `heat_capacity`     – Dimensionless heat capacity per degree of freedom.
/// * `n_dof`             – Number of degrees of freedom.
pub fn optimal_exchange_spacing(target_acceptance: f64, heat_capacity: f64, n_dof: usize) -> f64 {
    (2.0 / (n_dof as f64 * heat_capacity)).sqrt() * target_acceptance
}

// ─── RemdTrajectory ──────────────────────────────────────────────────────────

/// Trajectory of a full REMD simulation.
#[derive(Debug)]
pub struct RemdTrajectory {
    /// Current replicas (ordered by index, not temperature).
    pub replicas: Vec<RemdReplica>,
    /// Number of completed steps.
    pub n_steps: u32,
    /// Exchange log.
    pub log: RemdExchangeLog,
    /// Per-replica temperature history (replica index → list of temperatures).
    temperature_history: Vec<Vec<f64>>,
}

impl RemdTrajectory {
    /// Construct a new [`RemdTrajectory`] from a [`RemdConfig`].
    pub fn new(config: &RemdConfig) -> Self {
        let replicas: Vec<RemdReplica> = config
            .temperatures
            .iter()
            .enumerate()
            .map(|(id, &t)| RemdReplica::new(t, 0.0, id))
            .collect();
        let n = replicas.len();
        Self {
            replicas,
            n_steps: 0,
            log: RemdExchangeLog::new(),
            temperature_history: vec![Vec::new(); n],
        }
    }

    /// Advance the trajectory by one exchange round.
    ///
    /// Updates all replica potential energies from `energies`, then attempts
    /// exchanges between every adjacent pair using the corresponding element
    /// of `random_bits` to decide acceptance.
    ///
    /// # Arguments
    /// * `energies`     – New potential energies, one per replica (same order).
    /// * `random_bits`  – One boolean per adjacent pair: `true` means accept
    ///   regardless of probability (for deterministic testing).
    pub fn step_with_energies(&mut self, energies: &[f64], random_bits: &[bool]) {
        assert_eq!(energies.len(), self.replicas.len());

        // Record current temperatures before updates
        for (i, rep) in self.replicas.iter().enumerate() {
            self.temperature_history[i].push(rep.temperature);
        }

        // Update energies
        for (rep, &e) in self.replicas.iter_mut().zip(energies.iter()) {
            rep.potential_energy = e;
        }

        // Attempt adjacent exchanges
        let n_pairs = self.replicas.len().saturating_sub(1);
        let kb = 8.314e-3_f64;
        for pair in 0..n_pairs {
            let beta_i = 1.0 / (kb * self.replicas[pair].temperature);
            let beta_j = 1.0 / (kb * self.replicas[pair + 1].temperature);
            let e_i = self.replicas[pair].potential_energy;
            let e_j = self.replicas[pair + 1].potential_energy;
            let prob = metropolis_exchange_probability(e_i, e_j, beta_i, beta_j);

            let accept = if pair < random_bits.len() {
                random_bits[pair] && prob > 0.0
            } else {
                prob >= 1.0
            };

            let id_i = self.replicas[pair].replica_id;
            let id_j = self.replicas[pair + 1].replica_id;
            self.log.record(id_i, id_j, prob, accept);

            let (left, right) = self.replicas.split_at_mut(pair + 1);
            attempt_exchange(&mut left[pair], &mut right[0], kb, accept);
        }

        self.n_steps += 1;
    }

    /// Current temperature ladder (temperatures in replica order).
    pub fn temperature_ladder(&self) -> Vec<f64> {
        self.replicas.iter().map(|r| r.temperature).collect()
    }

    /// Estimate round-trip time (in steps) for a replica to travel from the
    /// lowest to the highest temperature and back.
    ///
    /// Returns `None` if the replica has not yet completed a full round trip.
    pub fn round_trip_time(&self, replica_id: usize) -> Option<u32> {
        if replica_id >= self.temperature_history.len() {
            return None;
        }
        let history = &self.temperature_history[replica_id];
        if history.is_empty() {
            return None;
        }

        let t_min = self
            .replicas
            .iter()
            .map(|r| r.temperature)
            .fold(f64::INFINITY, f64::min);
        let t_max = self
            .replicas
            .iter()
            .map(|r| r.temperature)
            .fold(f64::NEG_INFINITY, f64::max);

        // Also check original config temperatures from history
        let hist_min = history.iter().cloned().fold(f64::INFINITY, f64::min);
        let hist_max = history.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let effective_min = t_min.min(hist_min);
        let effective_max = t_max.max(hist_max);

        // Find a round trip: min → max → min (or max → min → max)
        let mut at_min = history[0] <= effective_min + 1e-6;
        let mut at_max = history[0] >= effective_max - 1e-6;
        let mut phase = if at_min {
            1i32
        } else if at_max {
            -1i32
        } else {
            0i32
        };
        let mut start_step = 0usize;

        for (step, &t) in history.iter().enumerate() {
            at_min = t <= effective_min + 1e-6;
            at_max = t >= effective_max - 1e-6;

            match phase {
                0 => {
                    if at_min {
                        phase = 1;
                        start_step = step;
                    } else if at_max {
                        phase = -1;
                        start_step = step;
                    }
                }
                1 if at_max => {
                    phase = 2;
                }
                2 if at_min => {
                    return Some((step - start_step) as u32);
                }
                -1 if at_min => {
                    phase = -2;
                }
                -2 if at_max => {
                    return Some((step - start_step) as u32);
                }
                _ => {}
            }
        }
        None
    }
}

// ─── demux_trajectory ────────────────────────────────────────────────────────

/// Extract the temperature history for a given `replica_id` across all frames.
///
/// Returns a `Vec`f64` of length equal to the number of steps recorded so far.
pub fn demux_trajectory(traj: &RemdTrajectory, replica_id: usize) -> Vec<f64> {
    if replica_id >= traj.temperature_history.len() {
        return Vec::new();
    }
    traj.temperature_history[replica_id].clone()
}

// ─── weighted_histogram_analysis ─────────────────────────────────────────────

/// Simplified WHAM-like density of states estimation.
///
/// Takes energy histograms from multiple temperatures and returns a normalised
/// (sums to 1) density-of-states vector over `n_bins` bins.
///
/// # Arguments
/// * `energies`     – One inner `Vec`f64` per temperature, listing sampled energies.
/// * `temperatures` – Temperatures corresponding to each inner vector (K).
/// * `kb`           – Boltzmann constant (kJ mol⁻¹ K⁻¹).
/// * `n_bins`       – Number of histogram bins.
pub fn weighted_histogram_analysis(
    energies: &[Vec<f64>],
    temperatures: &[f64],
    kb: f64,
    n_bins: usize,
) -> Vec<f64> {
    if energies.is_empty() || n_bins == 0 {
        return vec![0.0; n_bins];
    }

    // Find global energy range
    let all_energies: Vec<f64> = energies.iter().flatten().cloned().collect();
    if all_energies.is_empty() {
        return vec![0.0; n_bins];
    }
    let e_min = all_energies.iter().cloned().fold(f64::INFINITY, f64::min);
    let e_max = all_energies
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);

    if (e_max - e_min).abs() < 1e-12 {
        // All energies identical — put everything in first bin
        let mut dos = vec![0.0; n_bins];
        dos[0] = 1.0;
        return dos;
    }

    let bin_width = (e_max - e_min) / n_bins as f64;

    // Weighted histogram: weight each sample by 1/Z_k where Z_k is estimated
    // from the Boltzmann factor at the reference (lowest) temperature.
    let t_ref = temperatures.iter().cloned().fold(f64::INFINITY, f64::min);
    let beta_ref = 1.0 / (kb * t_ref);

    let mut dos = vec![0.0f64; n_bins];

    for (temp_energies, &temp) in energies.iter().zip(temperatures.iter()) {
        let beta = 1.0 / (kb * temp);
        for &e in temp_energies {
            let bin = ((e - e_min) / bin_width).floor() as usize;
            let bin = bin.min(n_bins - 1);
            // Reweight to reference temperature
            let weight = ((beta_ref - beta) * e).exp();
            dos[bin] += weight;
        }
    }

    // Normalise
    let total: f64 = dos.iter().sum();
    if total > 0.0 {
        for d in dos.iter_mut() {
            *d /= total;
        }
    }

    dos
}

// ─── RemdAnalysis ─────────────────────────────────────────────────────────────

/// Post-simulation analysis utilities for REMD trajectories.
pub struct RemdAnalysis;

impl RemdAnalysis {
    /// Compute the mean potential energy per replica.
    ///
    /// `energy_history[replica][step]` → mean over steps.
    #[allow(dead_code)]
    pub fn mean_energies(energy_history: &[Vec<f64>]) -> Vec<f64> {
        energy_history
            .iter()
            .map(|hist| {
                if hist.is_empty() {
                    return 0.0;
                }
                hist.iter().sum::<f64>() / hist.len() as f64
            })
            .collect()
    }

    /// Compute the variance of potential energy per replica.
    #[allow(dead_code)]
    pub fn energy_variance(energy_history: &[Vec<f64>]) -> Vec<f64> {
        energy_history
            .iter()
            .map(|hist| {
                if hist.len() < 2 {
                    return 0.0;
                }
                let n = hist.len() as f64;
                let mean = hist.iter().sum::<f64>() / n;
                hist.iter().map(|&e| (e - mean) * (e - mean)).sum::<f64>() / n
            })
            .collect()
    }

    /// Exchange matrix: `matrix[i][j]` = total number of swaps between replicas i and j.
    #[allow(dead_code)]
    pub fn exchange_matrix(log: &RemdExchangeLog, n_replicas: usize) -> Vec<Vec<u32>> {
        let mut mat = vec![vec![0u32; n_replicas]; n_replicas];
        for &(i, j, _, accepted) in &log.exchanges {
            if accepted && i < n_replicas && j < n_replicas {
                mat[i][j] += 1;
                mat[j][i] += 1;
            }
        }
        mat
    }

    /// Total number of accepted exchanges.
    #[allow(dead_code)]
    pub fn total_exchanges(log: &RemdExchangeLog) -> usize {
        log.exchanges.iter().filter(|e| e.3).count()
    }

    /// Optimal temperature ratio for target acceptance rate.
    ///
    /// Based on the empirical formula: T_{i+1}/T_i ≈ 1 + sqrt(2/N_dof)/sqrt(k_B*T/sigma_E)
    #[allow(dead_code)]
    pub fn optimal_temperature_ratio(n_dof: usize, target_rate: f64) -> f64 {
        1.0 + target_rate * (2.0 / n_dof as f64).sqrt()
    }

    /// Build an exponential temperature ladder from T_min to T_max.
    ///
    /// Same as `RemdConfig::geometric_temperatures` but returns the ratio.
    #[allow(dead_code)]
    pub fn exponential_spacing(t_min: f64, t_max: f64, n: usize) -> (Vec<f64>, f64) {
        let temps = RemdConfig::geometric_temperatures(t_min, t_max, n);
        let ratio = if n > 1 {
            (t_max / t_min).powf(1.0 / (n - 1) as f64)
        } else {
            1.0
        };
        (temps, ratio)
    }

    /// Compute the effective sample size (ESS) for a set of energies.
    ///
    /// ESS ≈ n² / Σ(w_i)² where w_i = exp(-β E_i) / Z.
    #[allow(dead_code)]
    pub fn effective_sample_size(energies: &[f64], beta: f64) -> f64 {
        if energies.is_empty() {
            return 0.0;
        }
        let e_min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
        let weights: Vec<f64> = energies
            .iter()
            .map(|&e| (-(beta * (e - e_min))).exp())
            .collect();
        let sum_w: f64 = weights.iter().sum();
        if sum_w < 1e-30 {
            return 0.0;
        }
        let sum_w2: f64 = weights.iter().map(|&w| (w / sum_w) * (w / sum_w)).sum();
        if sum_w2 < 1e-30 {
            return energies.len() as f64;
        }
        1.0 / sum_w2
    }
}

// ─── Potential energy vs replica statistics ───────────────────────────────────

/// Per-replica energy statistics for analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ReplicaEnergyStats {
    /// Replica index.
    pub replica_id: usize,
    /// Temperature (K).
    pub temperature: f64,
    /// Mean potential energy (kJ/mol).
    pub mean_energy: f64,
    /// Energy standard deviation.
    pub std_energy: f64,
    /// Number of samples.
    pub n_samples: usize,
}

impl ReplicaEnergyStats {
    /// Compute statistics from a history of energies for one replica.
    #[allow(dead_code)]
    pub fn from_history(replica_id: usize, temperature: f64, energies: &[f64]) -> Self {
        let n = energies.len();
        if n == 0 {
            return Self {
                replica_id,
                temperature,
                mean_energy: 0.0,
                std_energy: 0.0,
                n_samples: 0,
            };
        }
        let mean = energies.iter().sum::<f64>() / n as f64;
        let var = energies
            .iter()
            .map(|&e| (e - mean) * (e - mean))
            .sum::<f64>()
            / n as f64;
        Self {
            replica_id,
            temperature,
            mean_energy: mean,
            std_energy: var.sqrt(),
            n_samples: n,
        }
    }

    /// Heat capacity estimate from energy fluctuations: Cv = `δE²` / (k_B T²).
    #[allow(dead_code)]
    pub fn heat_capacity(&self, kb: f64) -> f64 {
        if self.temperature < 1e-10 {
            return 0.0;
        }
        self.std_energy * self.std_energy / (kb * self.temperature * self.temperature)
    }
}

// ─── Metropolis criterion variants ────────────────────────────────────────────

/// Generalised exchange criterion that accounts for pressure-volume work
/// (NPT ensemble REMD).
///
/// Δ = (β_i - β_j)(E_i - E_j) + (β_i P_i - β_j P_j)(V_i - V_j)
#[allow(dead_code)]
pub fn npt_exchange_probability(
    e_i: f64,
    e_j: f64,
    beta_i: f64,
    beta_j: f64,
    p_i: f64,
    p_j: f64,
    v_i: f64,
    v_j: f64,
) -> f64 {
    let delta = (beta_i - beta_j) * (e_i - e_j) + (beta_i * p_i - beta_j * p_j) * (v_i - v_j);
    if delta >= 0.0 { 1.0 } else { delta.exp() }
}

/// Temperature swap acceptance rate estimator based on energy overlap.
///
/// Uses the analytic approximation for Gaussian-distributed energies:
/// P_acc ≈ erfc(|β_i - β_j| * σ_E / √2) where σ_E is the RMS energy fluctuation.
#[allow(dead_code)]
pub fn expected_acceptance_rate(beta_i: f64, beta_j: f64, sigma_e: f64) -> f64 {
    let x = (beta_i - beta_j).abs() * sigma_e / std::f64::consts::SQRT_2;
    // Approximate erfc(x) for small x
    if x < 1e-10 {
        return 1.0;
    }
    let erfc_approx = 1.0 - erf_approx(x);
    erfc_approx.max(0.0)
}

/// Abramowitz & Stegun approximation for erf(x), valid for x >= 0.
fn erf_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    1.0 - poly * (-x * x).exp()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── RemdReplica ──────────────────────────────────────────────────────────

    #[test]
    fn test_replica_new() {
        let r = RemdReplica::new(300.0, -100.0, 0);
        assert_eq!(r.temperature, 300.0);
        assert_eq!(r.potential_energy, -100.0);
        assert_eq!(r.replica_id, 0);
        assert_eq!(r.exchange_count, 0);
    }

    #[test]
    fn test_replica_fields() {
        let r = RemdReplica::new(500.0, 42.0, 3);
        assert_eq!(r.replica_id, 3);
        assert_eq!(r.exchange_count, 0);
    }

    // ── RemdConfig ───────────────────────────────────────────────────────────

    #[test]
    fn test_config_new_endpoints() {
        let config = RemdConfig::new(300.0, 600.0, 4);
        assert_eq!(config.temperatures.len(), 4);
        assert!((config.temperatures[0] - 300.0).abs() < 1e-6);
        assert!((config.temperatures[3] - 600.0).abs() < 1e-6);
    }

    #[test]
    fn test_config_kb() {
        let config = RemdConfig::new(300.0, 600.0, 4);
        assert!((config.kb - 8.314e-3).abs() < 1e-10);
    }

    #[test]
    fn test_geometric_temperatures_monotonic() {
        let temps = RemdConfig::geometric_temperatures(300.0, 600.0, 5);
        for w in temps.windows(2) {
            assert!(w[1] > w[0], "temperatures must be strictly increasing");
        }
    }

    #[test]
    fn test_geometric_temperatures_ratio_constant() {
        let temps = RemdConfig::geometric_temperatures(300.0, 900.0, 4);
        let ratios: Vec<f64> = temps.windows(2).map(|w| w[1] / w[0]).collect();
        let r0 = ratios[0];
        for &r in &ratios[1..] {
            assert!(
                (r - r0).abs() < 1e-8,
                "constant ratio expected, got {r0} vs {r}"
            );
        }
    }

    #[test]
    fn test_geometric_temperatures_single() {
        let temps = RemdConfig::geometric_temperatures(300.0, 900.0, 1);
        assert_eq!(temps.len(), 1);
        assert!((temps[0] - 300.0).abs() < 1e-6);
    }

    #[test]
    fn test_geometric_temperatures_two() {
        let temps = RemdConfig::geometric_temperatures(300.0, 600.0, 2);
        assert_eq!(temps.len(), 2);
        assert!((temps[0] - 300.0).abs() < 1e-6);
        assert!((temps[1] - 600.0).abs() < 1e-6);
    }

    // ── metropolis_exchange_probability ──────────────────────────────────────

    #[test]
    fn test_metropolis_delta_positive_gives_one() {
        // beta_i > beta_j, e_i > e_j → delta > 0 → prob = 1
        let beta_i = 1.0 / (8.314e-3 * 300.0);
        let beta_j = 1.0 / (8.314e-3 * 600.0);
        let p = metropolis_exchange_probability(-50.0, -200.0, beta_i, beta_j);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metropolis_delta_negative_less_than_one() {
        let beta_i = 1.0 / (8.314e-3 * 300.0);
        let beta_j = 1.0 / (8.314e-3 * 600.0);
        let p = metropolis_exchange_probability(-200.0, -50.0, beta_i, beta_j);
        assert!(p < 1.0 && p > 0.0);
    }

    #[test]
    fn test_metropolis_equal_energies() {
        let beta_i = 1.0 / (8.314e-3 * 300.0);
        let beta_j = 1.0 / (8.314e-3 * 600.0);
        let p = metropolis_exchange_probability(-100.0, -100.0, beta_i, beta_j);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metropolis_in_range() {
        let beta_i = 1.0 / (8.314e-3 * 300.0);
        let beta_j = 1.0 / (8.314e-3 * 600.0);
        for e_i in [-200.0f64, -100.0, 0.0, 100.0] {
            for e_j in [-200.0f64, -100.0, 0.0, 100.0] {
                let p = metropolis_exchange_probability(e_i, e_j, beta_i, beta_j);
                assert!((0.0..=1.0).contains(&p), "prob out of range: {p}");
            }
        }
    }

    // ── attempt_exchange ─────────────────────────────────────────────────────

    #[test]
    fn test_attempt_exchange_accepted_swaps_temperature() {
        let mut r0 = RemdReplica::new(300.0, -100.0, 0);
        let mut r1 = RemdReplica::new(600.0, -80.0, 1);
        attempt_exchange(&mut r0, &mut r1, 8.314e-3, true);
        assert!((r0.temperature - 600.0).abs() < 1e-6);
        assert!((r1.temperature - 300.0).abs() < 1e-6);
    }

    #[test]
    fn test_attempt_exchange_accepted_increments_counts() {
        let mut r0 = RemdReplica::new(300.0, -100.0, 0);
        let mut r1 = RemdReplica::new(600.0, -80.0, 1);
        attempt_exchange(&mut r0, &mut r1, 8.314e-3, true);
        assert_eq!(r0.exchange_count, 1);
        assert_eq!(r1.exchange_count, 1);
    }

    #[test]
    fn test_attempt_exchange_rejected_no_swap() {
        let mut r0 = RemdReplica::new(300.0, -100.0, 0);
        let mut r1 = RemdReplica::new(600.0, -80.0, 1);
        attempt_exchange(&mut r0, &mut r1, 8.314e-3, false);
        assert!((r0.temperature - 300.0).abs() < 1e-6);
        assert!((r1.temperature - 600.0).abs() < 1e-6);
        assert_eq!(r0.exchange_count, 0);
        assert_eq!(r1.exchange_count, 0);
    }

    // ── RemdExchangeLog ──────────────────────────────────────────────────────

    #[test]
    fn test_log_empty() {
        let log = RemdExchangeLog::new();
        assert_eq!(log.acceptance_rate(), 0.0);
    }

    #[test]
    fn test_log_record_and_rate() {
        let mut log = RemdExchangeLog::new();
        log.record(0, 1, 1.0, true);
        log.record(0, 1, 0.5, true);
        log.record(0, 1, 0.1, false);
        log.record(0, 1, 0.1, false);
        let rate = log.acceptance_rate();
        assert!((rate - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_log_acceptance_rate_for_pair() {
        let mut log = RemdExchangeLog::new();
        log.record(0, 1, 1.0, true);
        log.record(1, 2, 0.5, false);
        log.record(0, 1, 0.8, true);
        assert!((log.acceptance_rate_for_pair(0, 1) - 1.0).abs() < 1e-10);
        assert!((log.acceptance_rate_for_pair(1, 2) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_log_acceptance_rate_for_missing_pair() {
        let log = RemdExchangeLog::new();
        assert_eq!(log.acceptance_rate_for_pair(5, 6), 0.0);
    }

    // ── optimal_exchange_spacing ─────────────────────────────────────────────

    #[test]
    fn test_optimal_spacing_positive() {
        let spacing = optimal_exchange_spacing(0.2, 1.5, 100);
        assert!(spacing > 0.0);
    }

    #[test]
    fn test_optimal_spacing_scales_with_target() {
        let s1 = optimal_exchange_spacing(0.1, 1.5, 100);
        let s2 = optimal_exchange_spacing(0.2, 1.5, 100);
        assert!((s2 / s1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_optimal_spacing_decreases_with_ndof() {
        let s_small = optimal_exchange_spacing(0.2, 1.5, 10);
        let s_large = optimal_exchange_spacing(0.2, 1.5, 1000);
        assert!(s_large < s_small);
    }

    // ── RemdTrajectory ───────────────────────────────────────────────────────

    #[test]
    fn test_trajectory_new() {
        let config = RemdConfig::new(300.0, 600.0, 4);
        let traj = RemdTrajectory::new(&config);
        assert_eq!(traj.replicas.len(), 4);
        assert_eq!(traj.n_steps, 0);
    }

    #[test]
    fn test_trajectory_temperature_ladder() {
        let config = RemdConfig::new(300.0, 600.0, 3);
        let traj = RemdTrajectory::new(&config);
        let ladder = traj.temperature_ladder();
        assert_eq!(ladder.len(), 3);
        assert!((ladder[0] - 300.0).abs() < 1e-6);
        assert!((ladder[2] - 600.0).abs() < 1e-6);
    }

    #[test]
    fn test_trajectory_step_increments_n_steps() {
        let config = RemdConfig::new(300.0, 600.0, 3);
        let mut traj = RemdTrajectory::new(&config);
        let energies = vec![-100.0, -80.0, -60.0];
        let bits = vec![false, false];
        traj.step_with_energies(&energies, &bits);
        assert_eq!(traj.n_steps, 1);
    }

    #[test]
    fn test_trajectory_step_updates_energies() {
        let config = RemdConfig::new(300.0, 600.0, 3);
        let mut traj = RemdTrajectory::new(&config);
        let energies = vec![-111.0, -222.0, -333.0];
        let bits = vec![false, false];
        traj.step_with_energies(&energies, &bits);
        // Energies updated (no swap since bits = false)
        assert!((traj.replicas[0].potential_energy - (-111.0)).abs() < 1e-6);
    }

    #[test]
    fn test_demux_trajectory() {
        let config = RemdConfig::new(300.0, 600.0, 2);
        let mut traj = RemdTrajectory::new(&config);
        for _ in 0..5 {
            traj.step_with_energies(&[-100.0, -80.0], &[false]);
        }
        let hist = demux_trajectory(&traj, 0);
        assert_eq!(hist.len(), 5);
    }

    #[test]
    fn test_demux_trajectory_invalid_id() {
        let config = RemdConfig::new(300.0, 600.0, 2);
        let traj = RemdTrajectory::new(&config);
        let hist = demux_trajectory(&traj, 99);
        assert!(hist.is_empty());
    }

    #[test]
    fn test_round_trip_none_initially() {
        let config = RemdConfig::new(300.0, 600.0, 2);
        let traj = RemdTrajectory::new(&config);
        assert!(traj.round_trip_time(0).is_none());
    }

    // ── weighted_histogram_analysis ──────────────────────────────────────────

    #[test]
    fn test_wham_empty() {
        let dos = weighted_histogram_analysis(&[], &[], 8.314e-3, 10);
        assert_eq!(dos.len(), 10);
        assert!(dos.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn test_wham_normalised() {
        let energies = vec![vec![-200.0, -180.0, -160.0], vec![-150.0, -140.0, -130.0]];
        let temperatures = vec![300.0, 600.0];
        let dos = weighted_histogram_analysis(&energies, &temperatures, 8.314e-3, 10);
        let total: f64 = dos.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-6,
            "DOS should sum to 1, got {total}"
        );
    }

    #[test]
    fn test_wham_n_bins() {
        let energies = vec![vec![-100.0, -80.0, -60.0]];
        let temperatures = vec![300.0];
        let dos = weighted_histogram_analysis(&energies, &temperatures, 8.314e-3, 7);
        assert_eq!(dos.len(), 7);
    }

    #[test]
    fn test_wham_uniform_energies() {
        // All energies at the same value → single bin gets all weight
        let energies = vec![vec![-100.0; 10]];
        let temperatures = vec![300.0];
        let dos = weighted_histogram_analysis(&energies, &temperatures, 8.314e-3, 5);
        let total: f64 = dos.iter().sum();
        assert!((total - 1.0).abs() < 1e-6);
    }

    // ── RemdAnalysis tests ───────────────────────────────────────────────────

    #[test]
    fn test_remd_mean_energies() {
        let history = vec![vec![-100.0, -200.0], vec![-80.0, -80.0]];
        let means = RemdAnalysis::mean_energies(&history);
        assert!((means[0] - (-150.0)).abs() < 1e-10);
        assert!((means[1] - (-80.0)).abs() < 1e-10);
    }

    #[test]
    fn test_remd_energy_variance() {
        let history = vec![vec![-100.0, -100.0]];
        let var = RemdAnalysis::energy_variance(&history);
        assert!(var[0].abs() < 1e-10, "zero variance for equal energies");
    }

    #[test]
    fn test_remd_exchange_matrix_symmetric() {
        let mut log = RemdExchangeLog::new();
        log.record(0, 1, 0.9, true);
        log.record(1, 2, 0.7, true);
        let mat = RemdAnalysis::exchange_matrix(&log, 3);
        assert_eq!(mat[0][1], mat[1][0]);
        assert_eq!(mat[1][2], mat[2][1]);
    }

    #[test]
    fn test_remd_total_exchanges() {
        let mut log = RemdExchangeLog::new();
        log.record(0, 1, 0.9, true);
        log.record(0, 1, 0.5, false);
        log.record(1, 2, 0.7, true);
        assert_eq!(RemdAnalysis::total_exchanges(&log), 2);
    }

    #[test]
    fn test_exponential_spacing_ratio() {
        let (temps, ratio) = RemdAnalysis::exponential_spacing(300.0, 1200.0, 5);
        assert_eq!(temps.len(), 5);
        // Ratio should be (1200/300)^(1/4) = 4^0.25 ≈ 1.414
        assert!((ratio - 4.0_f64.powf(0.25)).abs() < 1e-8, "ratio = {ratio}");
    }

    #[test]
    fn test_optimal_temperature_ratio_positive() {
        let r = RemdAnalysis::optimal_temperature_ratio(100, 0.2);
        assert!(r > 1.0, "ratio should be > 1, got {r}");
    }

    #[test]
    fn test_effective_sample_size_uniform() {
        // Equal weights → ESS = n
        let energies = vec![0.0; 10];
        let ess = RemdAnalysis::effective_sample_size(&energies, 0.001);
        assert!((ess - 10.0).abs() < 1e-6, "ESS for uniform = {ess}");
    }

    // ── ReplicaEnergyStats tests ─────────────────────────────────────────────

    #[test]
    fn test_replica_stats_mean() {
        let energies = vec![-100.0, -200.0, -150.0];
        let stats = ReplicaEnergyStats::from_history(0, 300.0, &energies);
        assert!((stats.mean_energy - (-150.0)).abs() < 1e-10);
    }

    #[test]
    fn test_replica_stats_std() {
        let energies = vec![0.0, 0.0, 0.0];
        let stats = ReplicaEnergyStats::from_history(0, 300.0, &energies);
        assert!(stats.std_energy.abs() < 1e-10);
    }

    #[test]
    fn test_replica_stats_empty() {
        let stats = ReplicaEnergyStats::from_history(0, 300.0, &[]);
        assert_eq!(stats.n_samples, 0);
    }

    #[test]
    fn test_replica_heat_capacity_positive() {
        let energies: Vec<f64> = (0..20).map(|i| -100.0 + i as f64 * 5.0).collect();
        let stats = ReplicaEnergyStats::from_history(0, 300.0, &energies);
        let cv = stats.heat_capacity(8.314e-3);
        assert!(cv > 0.0, "heat capacity should be positive, got {cv}");
    }

    // ── npt_exchange_probability tests ───────────────────────────────────────

    #[test]
    fn test_npt_exchange_prob_in_range() {
        let kb = 8.314e-3;
        let beta_i = 1.0 / (kb * 300.0);
        let beta_j = 1.0 / (kb * 600.0);
        let p = npt_exchange_probability(-100.0, -80.0, beta_i, beta_j, 1e5, 1e5, 1.0, 1.0);
        assert!((0.0..=1.0).contains(&p), "prob = {p}");
    }

    // ── expected_acceptance_rate tests ───────────────────────────────────────

    #[test]
    fn test_expected_acceptance_rate_equal_temps() {
        // Same temperature → delta_beta = 0, rate ≈ 1
        let rate = expected_acceptance_rate(0.001, 0.001, 10.0);
        assert!((rate - 1.0).abs() < 1e-6, "rate at equal temps = {rate}");
    }

    #[test]
    fn test_expected_acceptance_rate_decreases_with_delta_beta() {
        let r1 = expected_acceptance_rate(0.001, 0.0015, 10.0);
        let r2 = expected_acceptance_rate(0.001, 0.003, 10.0);
        assert!(r1 > r2, "larger temp gap should give lower acceptance rate");
    }
}

// ─── Hamiltonian REMD ─────────────────────────────────────────────────────────

/// A single replica in Hamiltonian REMD (H-REMD).
///
/// In H-REMD replicas run at the same temperature but with different Hamiltonians
/// (e.g. different force-field parameters, alchemical lambda states, or solute
/// scaling). Exchange is attempted between adjacent lambda states.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HremdReplica {
    /// Lambda parameter (0 = reference, 1 = target Hamiltonian).
    pub lambda: f64,
    /// Potential energy at this lambda.
    pub energy: f64,
    /// Potential energy evaluated at the *other* end-state (needed for exchange).
    pub energy_at_neighbour_lambda: f64,
    /// Simulation temperature (K) – identical for all H-REMD replicas.
    pub temperature: f64,
}

impl HremdReplica {
    /// Create a new H-REMD replica.
    #[allow(dead_code)]
    pub fn new(lambda: f64, temperature: f64) -> Self {
        Self {
            lambda,
            energy: 0.0,
            energy_at_neighbour_lambda: 0.0,
            temperature,
        }
    }
}

/// Metropolis exchange criterion for Hamiltonian REMD.
///
/// Replicas `i` (lambda_i, E_i(lambda_i), E_i(lambda_j)) and
///          `j` (lambda_j, E_j(lambda_j), E_j(lambda_i)):
///
/// Δ = β * \[ (E_i(lambda_j) + E_j(lambda_i)) − (E_i(lambda_i) + E_j(lambda_j)) \]
///
/// In the simplified form where only one cross-evaluation is stored:
/// Δ = β * (E_i(lambda_j) - E_i(lambda_i))   (standard for solute scaling)
#[allow(dead_code)]
pub fn hremd_exchange_probability(e_i_own: f64, e_i_neighbour: f64, beta: f64) -> f64 {
    let delta = beta * (e_i_neighbour - e_i_own);
    if delta <= 0.0 { 1.0 } else { (-delta).exp() }
}

// ─── REST2 (Replica Exchange with Solute Tempering 2) ─────────────────────────

/// REST2 effective Hamiltonian scaling for solute atoms.
///
/// Scales solute-solute interactions by `beta_m / beta_0` and
/// solute-solvent interactions by `sqrt(beta_m / beta_0)`,
/// where `beta_m` is the effective inverse temperature of the replica and
/// `beta_0` is the physical temperature beta.
///
/// # Arguments
/// * `e_ss`  – Solute-solute interaction energy.
/// * `e_sv`  – Solute-solvent interaction energy.
/// * `e_vv`  – Solvent-solvent interaction energy (unscaled).
/// * `scale` – `beta_m / beta_0` (> 1 for heated solute replicas).
#[allow(dead_code)]
pub fn rest2_total_energy(e_ss: f64, e_sv: f64, e_vv: f64, scale: f64) -> f64 {
    let scaled_ss = e_ss * scale;
    let scaled_sv = e_sv * scale.sqrt();
    scaled_ss + scaled_sv + e_vv
}

/// Compute the REST2 exchange probability between two replicas.
///
/// Based on the Hamiltonian REMD criterion applied to REST2 energies.
#[allow(dead_code)]
pub fn rest2_exchange_probability(
    e_ss_i: f64,
    e_sv_i: f64,
    e_vv_i: f64,
    scale_i: f64,
    e_ss_j: f64,
    e_sv_j: f64,
    e_vv_j: f64,
    scale_j: f64,
    beta: f64,
) -> f64 {
    let u_i_at_i = rest2_total_energy(e_ss_i, e_sv_i, e_vv_i, scale_i);
    let u_j_at_j = rest2_total_energy(e_ss_j, e_sv_j, e_vv_j, scale_j);
    // After swap: replica i takes scale_j, replica j takes scale_i
    let u_i_at_j = rest2_total_energy(e_ss_i, e_sv_i, e_vv_i, scale_j);
    let u_j_at_i = rest2_total_energy(e_ss_j, e_sv_j, e_vv_j, scale_i);

    let delta = beta * ((u_i_at_j + u_j_at_i) - (u_i_at_i + u_j_at_j));
    if delta <= 0.0 { 1.0 } else { (-delta).exp() }
}

// ─── REMD temperature optimiser ───────────────────────────────────────────────

/// Optimise a temperature ladder to achieve a target acceptance rate.
///
/// Uses the iterative algorithm:
/// 1. Start with the geometric ladder.
/// 2. Estimate acceptance from energy fluctuations at each temperature.
/// 3. Rescale gaps to meet the target acceptance rate.
///
/// Returns the adjusted temperature ladder.
#[allow(dead_code)]
pub fn optimise_temperature_ladder(
    t_min: f64,
    t_max: f64,
    n_replicas: usize,
    energy_sigma: f64,
    target_acceptance: f64,
    kb: f64,
    max_iter: usize,
) -> Vec<f64> {
    if n_replicas <= 1 {
        return vec![t_min];
    }
    let mut temps = RemdConfig::geometric_temperatures(t_min, t_max, n_replicas);
    for _ in 0..max_iter {
        let mut converged = true;
        for i in 0..n_replicas - 1 {
            let beta_i = 1.0 / (kb * temps[i]);
            let beta_j = 1.0 / (kb * temps[i + 1]);
            let rate = expected_acceptance_rate(beta_i, beta_j, energy_sigma);
            if (rate - target_acceptance).abs() > 0.01 {
                // Adjust gap by scaling the ratio
                let ratio = temps[i + 1] / temps[i];
                let new_ratio = if rate > target_acceptance {
                    ratio * 1.05 // wider gap
                } else {
                    ratio * 0.95 // narrower gap
                };
                temps[i + 1] = (temps[i] * new_ratio).min(t_max * 2.0);
                converged = false;
            }
        }
        if converged {
            break;
        }
    }
    temps
}

// ─── Multi-dimensional exchange statistics ─────────────────────────────────────

/// Pair-wise exchange statistics matrix.
///
/// `ExchangeMatrix[i][j]` = (n_attempts, n_accepted) for the pair (i,j).
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct ExchangeMatrix {
    /// Number of replicas.
    pub n: usize,
    /// Attempt counts: `attempts[i][j]`.
    pub attempts: Vec<Vec<u32>>,
    /// Accept counts: `accepted[i][j]`.
    pub accepted: Vec<Vec<u32>>,
}

impl ExchangeMatrix {
    /// Create a new exchange matrix for `n` replicas.
    #[allow(dead_code)]
    pub fn new(n: usize) -> Self {
        Self {
            n,
            attempts: vec![vec![0u32; n]; n],
            accepted: vec![vec![0u32; n]; n],
        }
    }

    /// Record an exchange attempt between replicas `i` and `j`.
    #[allow(dead_code)]
    pub fn record(&mut self, i: usize, j: usize, accepted: bool) {
        if i < self.n && j < self.n {
            self.attempts[i][j] += 1;
            self.attempts[j][i] += 1;
            if accepted {
                self.accepted[i][j] += 1;
                self.accepted[j][i] += 1;
            }
        }
    }

    /// Acceptance rate for pair (i, j).
    #[allow(dead_code)]
    pub fn rate(&self, i: usize, j: usize) -> f64 {
        if i >= self.n || j >= self.n {
            return 0.0;
        }
        let att = self.attempts[i][j];
        if att == 0 {
            return 0.0;
        }
        self.accepted[i][j] as f64 / att as f64
    }

    /// Global acceptance rate across all pairs.
    #[allow(dead_code)]
    pub fn global_rate(&self) -> f64 {
        let mut total_att = 0u64;
        let mut total_acc = 0u64;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                total_att += self.attempts[i][j] as u64;
                total_acc += self.accepted[i][j] as u64;
            }
        }
        if total_att == 0 {
            return 0.0;
        }
        total_acc as f64 / total_att as f64
    }
}

// ─── Replica ladder diagnostics ───────────────────────────────────────────────

/// Compute diffusivity of replicas in temperature space.
///
/// Estimates the mean squared displacement per step in replica-index space
/// from a temperature history.
///
/// `temp_history[step]` = temperature observed at each step for one replica.
/// `temps` = the complete ordered temperature ladder.
///
/// Returns mean squared displacement in units of (replica index)^2 per step.
#[allow(dead_code)]
pub fn replica_diffusivity(temp_history: &[f64], temps: &[f64]) -> f64 {
    if temp_history.len() < 2 || temps.is_empty() {
        return 0.0;
    }
    let n_temps = temps.len() as f64;
    // Map temperature to fractional index
    let to_index = |t: f64| -> f64 {
        let t_min = temps.iter().cloned().fold(f64::INFINITY, f64::min);
        let t_max = temps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if (t_max - t_min).abs() < 1e-12 {
            return 0.0;
        }
        (t - t_min) / (t_max - t_min) * (n_temps - 1.0)
    };
    let indices: Vec<f64> = temp_history.iter().map(|&t| to_index(t)).collect();
    let n = indices.len() - 1;
    if n == 0 {
        return 0.0;
    }
    let msd: f64 = indices
        .windows(2)
        .map(|w| (w[1] - w[0]).powi(2))
        .sum::<f64>()
        / n as f64;
    msd
}

/// Estimate the ergodic measure for a REMD simulation.
///
/// Returns a value in \[0, 1\]: 0 means no mixing (all replicas stuck),
/// 1 means perfect mixing (every replica visits every temperature equally).
#[allow(dead_code)]
pub fn ergodic_measure(exchange_matrix: &ExchangeMatrix) -> f64 {
    let n = exchange_matrix.n;
    if n <= 1 {
        return 1.0;
    }
    let mut total_rate = 0.0;
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            total_rate += exchange_matrix.rate(i, j);
            count += 1;
        }
    }
    if count == 0 {
        return 0.0;
    }
    total_rate / count as f64
}

// ─── REMD with Solute Scaling (SS-REMD) ──────────────────────────────────────

/// Solute-scaling REMD: replicas differ by a scaling of solute torsion potential.
///
/// At lambda=0: solute torsion barrier is removed (free diffusion of dihedrals)
/// At lambda=1: full force field (standard MD)
///
/// Exchange probability uses the H-REMD criterion at a single temperature.
#[allow(dead_code)]
pub fn ss_remd_exchange_probability(
    torsion_i: f64,
    torsion_j: f64,
    lambda_i: f64,
    lambda_j: f64,
    beta: f64,
) -> f64 {
    // Energy difference due to lambda scaling
    // U_torsion(lambda) = lambda * E_torsion
    let u_i_at_i = lambda_i * torsion_i;
    let u_j_at_j = lambda_j * torsion_j;
    let u_i_at_j = lambda_j * torsion_i;
    let u_j_at_i = lambda_i * torsion_j;
    let delta = beta * ((u_i_at_j + u_j_at_i) - (u_i_at_i + u_j_at_j));
    if delta <= 0.0 { 1.0 } else { (-delta).exp() }
}

// ─── Tests for new H-REMD / REST2 / diagnostics ────────────────────────────────

#[cfg(test)]
mod tests_hremd {
    use super::*;

    // ── HremdReplica ──────────────────────────────────────────────────────────

    #[test]
    fn test_hremd_replica_new() {
        let r = HremdReplica::new(0.5, 300.0);
        assert!((r.lambda - 0.5).abs() < 1e-12);
        assert!((r.temperature - 300.0).abs() < 1e-12);
        assert_eq!(r.energy, 0.0);
    }

    #[test]
    fn test_hremd_exchange_probability_same_energy() {
        // Same energy at own and neighbour → full exchange
        let p = hremd_exchange_probability(-100.0, -100.0, 1.0);
        assert!((p - 1.0).abs() < 1e-12, "Same energy → prob=1, got {p}");
    }

    #[test]
    fn test_hremd_exchange_probability_lower_at_neighbour() {
        // Lower energy at neighbour → should accept (prob = 1)
        let p = hremd_exchange_probability(-100.0, -200.0, 0.001);
        assert!(
            (p - 1.0).abs() < 1e-12,
            "Lower E at neighbour → accept, got {p}"
        );
    }

    #[test]
    fn test_hremd_exchange_probability_range() {
        for delta in [-100.0f64, -10.0, 0.0, 10.0, 100.0] {
            let p = hremd_exchange_probability(-100.0, -100.0 + delta, 0.001);
            assert!(
                (0.0..=1.0).contains(&p),
                "prob out of range for delta={delta}: {p}"
            );
        }
    }

    // ── REST2 ────────────────────────────────────────────────────────────────

    #[test]
    fn test_rest2_total_energy_scale_one() {
        // scale=1 → energy unchanged
        let e = rest2_total_energy(-50.0, -30.0, -20.0, 1.0);
        assert!(
            (e - (-100.0)).abs() < 1e-10,
            "scale=1 → full energy, got {e}"
        );
    }

    #[test]
    fn test_rest2_total_energy_scale_zero_vv_only() {
        // scale=0 → only solvent-solvent contributes
        let e = rest2_total_energy(-50.0, -30.0, -20.0, 0.0);
        assert!(
            (e - (-20.0)).abs() < 1e-10,
            "scale=0 → only e_vv survives, got {e}"
        );
    }

    #[test]
    fn test_rest2_exchange_probability_identical_replicas() {
        let p =
            rest2_exchange_probability(-50.0, -30.0, -20.0, 1.0, -50.0, -30.0, -20.0, 1.0, 0.001);
        assert!(
            (p - 1.0).abs() < 1e-12,
            "Identical replicas → prob=1, got {p}"
        );
    }

    #[test]
    fn test_rest2_exchange_probability_in_range() {
        let p =
            rest2_exchange_probability(-50.0, -30.0, -20.0, 2.0, -40.0, -20.0, -15.0, 1.0, 0.001);
        assert!((0.0..=1.0).contains(&p), "REST2 prob out of range: {p}");
    }

    // ── Temperature ladder optimiser ──────────────────────────────────────────

    #[test]
    fn test_optimise_ladder_length() {
        let temps = optimise_temperature_ladder(300.0, 600.0, 4, 20.0, 0.2, 8.314e-3, 5);
        assert_eq!(temps.len(), 4, "Optimised ladder should have 4 replicas");
    }

    #[test]
    fn test_optimise_ladder_starts_at_t_min() {
        let temps = optimise_temperature_ladder(300.0, 600.0, 4, 20.0, 0.2, 8.314e-3, 5);
        assert!((temps[0] - 300.0).abs() < 1e-6);
    }

    #[test]
    fn test_optimise_ladder_single_replica() {
        let temps = optimise_temperature_ladder(300.0, 600.0, 1, 20.0, 0.2, 8.314e-3, 5);
        assert_eq!(temps.len(), 1);
        assert!((temps[0] - 300.0).abs() < 1e-6);
    }

    // ── ExchangeMatrix ────────────────────────────────────────────────────────

    #[test]
    fn test_exchange_matrix_new() {
        let mat = ExchangeMatrix::new(3);
        assert_eq!(mat.n, 3);
        assert_eq!(mat.rate(0, 1), 0.0);
    }

    #[test]
    fn test_exchange_matrix_record_accepted() {
        let mut mat = ExchangeMatrix::new(3);
        mat.record(0, 1, true);
        assert_eq!(mat.attempts[0][1], 1);
        assert_eq!(mat.accepted[0][1], 1);
        // Symmetric
        assert_eq!(mat.attempts[1][0], 1);
        assert_eq!(mat.accepted[1][0], 1);
    }

    #[test]
    fn test_exchange_matrix_record_rejected() {
        let mut mat = ExchangeMatrix::new(3);
        mat.record(0, 1, false);
        assert_eq!(mat.attempts[0][1], 1);
        assert_eq!(mat.accepted[0][1], 0);
    }

    #[test]
    fn test_exchange_matrix_rate() {
        let mut mat = ExchangeMatrix::new(3);
        mat.record(0, 1, true);
        mat.record(0, 1, true);
        mat.record(0, 1, false);
        let rate = mat.rate(0, 1);
        assert!((rate - 2.0 / 3.0).abs() < 1e-10, "rate = {rate}");
    }

    #[test]
    fn test_exchange_matrix_global_rate() {
        let mut mat = ExchangeMatrix::new(2);
        mat.record(0, 1, true);
        mat.record(0, 1, false);
        let g = mat.global_rate();
        assert!((g - 0.5).abs() < 1e-10, "global rate = {g}");
    }

    #[test]
    fn test_exchange_matrix_out_of_range() {
        let mat = ExchangeMatrix::new(2);
        assert_eq!(mat.rate(5, 6), 0.0); // out of bounds → 0
    }

    // ── Replica diffusivity ────────────────────────────────────────────────────

    #[test]
    fn test_replica_diffusivity_empty() {
        let d = replica_diffusivity(&[], &[300.0, 400.0]);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_replica_diffusivity_single_step() {
        let d = replica_diffusivity(&[300.0], &[300.0, 600.0]);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_replica_diffusivity_constant() {
        // Replica always stays at same temperature → diffusivity = 0
        let hist = vec![300.0; 10];
        let temps = vec![300.0, 450.0, 600.0];
        let d = replica_diffusivity(&hist, &temps);
        assert!(d.abs() < 1e-10, "No movement → diffusivity = 0, got {d}");
    }

    #[test]
    fn test_replica_diffusivity_alternating() {
        // Replica bounces between min and max each step → high diffusivity
        let hist = vec![300.0, 600.0, 300.0, 600.0, 300.0];
        let temps = vec![300.0, 600.0];
        let d = replica_diffusivity(&hist, &temps);
        assert!(
            d > 0.0,
            "Alternating replica should have positive diffusivity, got {d}"
        );
    }

    // ── Ergodic measure ────────────────────────────────────────────────────────

    #[test]
    fn test_ergodic_measure_no_exchanges() {
        let mat = ExchangeMatrix::new(3);
        let e = ergodic_measure(&mat);
        assert_eq!(e, 0.0, "No exchanges → ergodic measure = 0");
    }

    #[test]
    fn test_ergodic_measure_all_accepted() {
        let mut mat = ExchangeMatrix::new(3);
        for (i, j) in [(0, 1), (1, 2), (0, 2)] {
            mat.record(i, j, true);
        }
        let e = ergodic_measure(&mat);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "All accepted → ergodic measure = 1, got {e}"
        );
    }

    #[test]
    fn test_ergodic_measure_single_replica() {
        let mat = ExchangeMatrix::new(1);
        let e = ergodic_measure(&mat);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "Single replica → trivially ergodic"
        );
    }

    // ── SS-REMD ───────────────────────────────────────────────────────────────

    #[test]
    fn test_ss_remd_exchange_probability_same_lambda() {
        // Both replicas at same lambda → exchange probability = 1
        let p = ss_remd_exchange_probability(10.0, 10.0, 0.5, 0.5, 0.001);
        assert!((p - 1.0).abs() < 1e-12, "Same lambda → prob=1, got {p}");
    }

    #[test]
    fn test_ss_remd_exchange_probability_in_range() {
        let p = ss_remd_exchange_probability(20.0, 5.0, 0.8, 0.2, 0.001);
        assert!((0.0..=1.0).contains(&p), "SS-REMD prob out of range: {p}");
    }

    // ── Additional Metropolis tests ───────────────────────────────────────────

    #[test]
    fn test_metropolis_symmetry() {
        // P(i→j) may differ from P(j→i) in general, but swapping beta and E should work
        let beta_i = 1.0 / (8.314e-3 * 300.0);
        let beta_j = 1.0 / (8.314e-3 * 600.0);
        let e_i = -150.0;
        let e_j = -100.0;
        let p_fwd = metropolis_exchange_probability(e_i, e_j, beta_i, beta_j);
        let p_rev = metropolis_exchange_probability(e_j, e_i, beta_j, beta_i);
        // Both should be in [0,1]
        assert!((0.0..=1.0).contains(&p_fwd));
        assert!((0.0..=1.0).contains(&p_rev));
    }

    #[test]
    fn test_npt_exchange_prob_same_pv() {
        // When PV terms are equal, reduces to NVT criterion
        let kb = 8.314e-3;
        let beta_i = 1.0 / (kb * 300.0);
        let beta_j = 1.0 / (kb * 600.0);
        let p_nvt = metropolis_exchange_probability(-100.0, -80.0, beta_i, beta_j);
        let p_npt = npt_exchange_probability(
            -100.0, -80.0, beta_i, beta_j, 1.0, 1.0, // equal PV
            0.0, 0.0,
        );
        assert!(
            (p_nvt - p_npt).abs() < 1e-10,
            "NPT with PV=0 should equal NVT"
        );
    }

    #[test]
    fn test_wham_two_temperatures_positive_entries() {
        let energies = vec![vec![-300.0, -290.0, -280.0], vec![-200.0, -190.0, -180.0]];
        let temperatures = vec![300.0, 600.0];
        let dos = weighted_histogram_analysis(&energies, &temperatures, 8.314e-3, 5);
        assert!(
            dos.iter().any(|&d| d > 0.0),
            "DOS should have positive entries"
        );
    }
}
