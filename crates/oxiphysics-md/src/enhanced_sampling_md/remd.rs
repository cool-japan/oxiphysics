// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Replica exchange MD (REMD) controller: temperature REMD and Hamiltonian
//! REMD with the Metropolis swap criterion.
//!
//! Split out of [`super`] to keep every source file below 2000 lines; the
//! public items are re-exported from `crate::enhanced_sampling_md`.

use rand::Rng;
use rand::RngExt;

use super::KB_KJMOL;
use crate::error::MdError;

// ─────────────────────────────────────────────────────────────────────────────
// RemdController
// ─────────────────────────────────────────────────────────────────────────────

/// Swap scheduling strategy for REMD.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwapSchedule {
    /// Try neighbours in order 0↔1, 2↔3, … then 1↔2, 3↔4, … alternating.
    EvenOdd,
    /// Random pair selection.
    Random,
}

/// Replica Exchange Molecular Dynamics controller.
///
/// Supports temperature REMD and Hamiltonian REMD, with configurable swap
/// scheduling and detailed acceptance tracking per pair.
///
/// # Swap criterion
///
/// Slot `k` runs at inverse temperature `beta_k = 1 / (kB T_k)` with the
/// Hamiltonian `U_k`, and currently holds the configuration `x_k`.  Exchanging
/// the configurations of slots `i` and `j` is accepted with the Metropolis
/// probability `min(1, exp(delta))` with
///
/// * temperature REMD (all slots share one Hamiltonian, `lambda == None`):
///   `delta = (beta_i - beta_j) (E_i - E_j)`,
/// * Hamiltonian REMD (`lambda == Some`):
///   `delta = -beta_i [U_i(x_j) - U_i(x_i)] - beta_j [U_j(x_i) - U_j(x_j)]`,
///   which needs the energy of every configuration under every Hamiltonian,
///   supplied through [`RemdController::set_cross_energies`].
///
/// Swaps that lower the combined weighted energy (`delta >= 0`) are always
/// accepted.
#[derive(Debug, Clone)]
pub struct RemdController {
    /// Temperature of each replica (K).
    pub temperatures: Vec<f64>,
    /// Potential energy of each replica (kJ/mol): the energy of the
    /// configuration held by slot `k` under the Hamiltonian of slot `k`.
    ///
    /// In Hamiltonian REMD the swap criterion reads only
    /// [`RemdController::cross_energies`]; `energies` is then kept equal to the
    /// diagonal of that matrix by [`RemdController::set_cross_energies`] and by
    /// accepted swaps.
    pub energies: Vec<f64>,
    /// Mapping: `replica_index[k]` is the configuration currently at slot `k`.
    pub replica_index: Vec<usize>,
    /// Swap schedule.
    pub schedule: SwapSchedule,
    /// Per-pair acceptance counters: `(attempts, accepted)`.
    pub pair_stats: Vec<(u64, u64)>,
    /// Global step counter.
    pub step: u64,
    /// Swap attempt interval (MD steps between swap attempts).
    pub swap_interval: u64,
    /// Hamiltonian REMD lambda parameters (one per replica). `None` → temperature REMD.
    ///
    /// The values label the Hamiltonians `U_lambda`; the controller itself never
    /// evaluates them, the caller computes the energies from them and passes
    /// the result via [`RemdController::set_cross_energies`].
    pub lambda: Option<Vec<f64>>,
    /// Hamiltonian REMD cross-energy matrix (kJ/mol): `cross_energies[k][l]` is
    /// `U_{lambda_k}(x_l)`, the energy of the configuration currently held by
    /// slot `l` evaluated with the Hamiltonian of slot `k`.
    ///
    /// `None` until supplied.  Only used when `lambda` is `Some`; swap proposals
    /// are rejected (and not counted) while it is missing.
    pub cross_energies: Option<Vec<Vec<f64>>>,
}

impl RemdController {
    /// Create a temperature REMD controller with a geometric ladder.
    pub fn new_geometric(t_min: f64, t_max: f64, n_replicas: usize, swap_interval: u64) -> Self {
        let mut temps = Vec::with_capacity(n_replicas);
        if n_replicas <= 1 {
            temps.push(t_min);
        } else {
            let ratio = (t_max / t_min).powf(1.0 / (n_replicas as f64 - 1.0));
            for i in 0..n_replicas {
                temps.push(t_min * ratio.powi(i as i32));
            }
        }
        let n = temps.len();
        Self {
            temperatures: temps,
            energies: vec![0.0; n],
            replica_index: (0..n).collect(),
            schedule: SwapSchedule::EvenOdd,
            pair_stats: vec![(0, 0); n.saturating_sub(1)],
            step: 0,
            swap_interval,
            lambda: None,
            cross_energies: None,
        }
    }

    /// Create a Hamiltonian REMD controller.
    ///
    /// All replicas run at `temperature`; replica `k` uses the Hamiltonian
    /// labelled `lambdas[k]`.  Call [`RemdController::set_cross_energies`] with
    /// the freshly computed cross-energy matrix before every swap attempt:
    /// without it every proposal is rejected.
    pub fn new_hamiltonian(temperature: f64, lambdas: Vec<f64>, swap_interval: u64) -> Self {
        let n = lambdas.len();
        Self {
            temperatures: vec![temperature; n],
            energies: vec![0.0; n],
            replica_index: (0..n).collect(),
            schedule: SwapSchedule::EvenOdd,
            pair_stats: vec![(0, 0); n.saturating_sub(1)],
            step: 0,
            swap_interval,
            lambda: Some(lambdas),
            cross_energies: None,
        }
    }

    /// Number of replicas.
    pub fn n_replicas(&self) -> usize {
        self.temperatures.len()
    }

    /// Update energies from the simulation.
    ///
    /// In temperature REMD these energies decide the swaps.  In Hamiltonian REMD
    /// they are informational only; the swaps use
    /// [`RemdController::set_cross_energies`].
    pub fn set_energies(&mut self, energies: &[f64]) {
        let n = energies.len().min(self.energies.len());
        self.energies[..n].copy_from_slice(&energies[..n]);
    }

    /// Supply the Hamiltonian-REMD cross-energy matrix.
    ///
    /// `u[k][l]` must be `U_{lambda_k}(x_l)` in kJ/mol: the energy of the
    /// configuration currently held by slot `l` under the Hamiltonian of slot
    /// `k`.  Call this after every MD segment, before
    /// [`RemdController::attempt_swaps`].  The diagonal also refreshes
    /// [`RemdController::energies`].
    ///
    /// # Errors
    /// Returns an error if this is a temperature-REMD controller, if `u` is not
    /// `n_replicas x n_replicas`, or if an entry is NaN or `-inf`.  The
    /// controller is left unchanged in that case.  An entry of `+inf` (a
    /// forbidden overlap) is allowed and makes the corresponding swap
    /// impossible.
    pub fn set_cross_energies(&mut self, u: Vec<Vec<f64>>) -> Result<(), MdError> {
        if self.lambda.is_none() {
            return Err(MdError::General(
                "cross-energies are only used by Hamiltonian REMD; \
                 create the controller with RemdController::new_hamiltonian"
                    .to_string(),
            ));
        }
        let n = self.n_replicas();
        if u.len() != n || u.iter().any(|row| row.len() != n) {
            return Err(MdError::General(format!(
                "cross-energy matrix must be {n} x {n}, got {} rows of lengths {:?}",
                u.len(),
                u.iter().map(Vec::len).collect::<Vec<_>>()
            )));
        }
        if u.iter()
            .flatten()
            .any(|x| x.is_nan() || *x == f64::NEG_INFINITY)
        {
            return Err(MdError::General(
                "cross-energy matrix entries must not be NaN or -inf".to_string(),
            ));
        }
        for (k, row) in u.iter().enumerate() {
            self.energies[k] = row[k];
        }
        self.cross_energies = Some(u);
        Ok(())
    }

    /// Natural logarithm `delta` of the Metropolis ratio for exchanging the
    /// configurations of slots `i` and `j` (see the type-level documentation).
    /// The swap is accepted with probability `min(1, exp(delta))`.
    ///
    /// Returns `None` if a slot index is out of range or, in Hamiltonian REMD,
    /// if no well-formed cross-energy matrix has been supplied.
    pub fn swap_log_ratio(&self, i: usize, j: usize) -> Option<f64> {
        let beta_i = 1.0 / (KB_KJMOL * *self.temperatures.get(i)?);
        let beta_j = 1.0 / (KB_KJMOL * *self.temperatures.get(j)?);
        if self.lambda.is_some() {
            let u = self.cross_energies.as_ref()?;
            let u_ii = *u.get(i)?.get(i)?;
            let u_ij = *u.get(i)?.get(j)?;
            let u_ji = *u.get(j)?.get(i)?;
            let u_jj = *u.get(j)?.get(j)?;
            Some(-beta_i * (u_ij - u_ii) - beta_j * (u_ji - u_jj))
        } else {
            let e_i = *self.energies.get(i)?;
            let e_j = *self.energies.get(j)?;
            Some((beta_i - beta_j) * (e_i - e_j))
        }
    }

    /// Attempt swaps according to the current schedule.
    /// Returns a list of pairs that were actually swapped.
    pub fn attempt_swaps(&mut self) -> Vec<(usize, usize)> {
        let n = self.n_replicas();
        if n < 2 {
            return Vec::new();
        }
        let mut rng = rand::rng();
        let mut swapped = Vec::new();

        match self.schedule {
            SwapSchedule::EvenOdd => {
                let offset = (self.step % 2) as usize;
                let mut i = offset;
                while i + 1 < n {
                    if self.try_swap(i, i + 1, &mut rng) {
                        swapped.push((i, i + 1));
                    }
                    i += 2;
                }
            }
            SwapSchedule::Random => {
                let i = rng.random_range(0..n - 1);
                if self.try_swap(i, i + 1, &mut rng) {
                    swapped.push((i, i + 1));
                }
            }
        }
        self.step += 1;
        swapped
    }

    /// Try to swap replicas `i` and `j` using the Metropolis criterion
    /// `min(1, exp(delta))` of [`RemdController::swap_log_ratio`].
    ///
    /// A proposal whose ratio cannot be evaluated (Hamiltonian REMD without
    /// cross-energies) is rejected and not counted in the pair statistics.
    fn try_swap(&mut self, i: usize, j: usize, rng: &mut impl Rng) -> bool {
        let Some(delta) = self.swap_log_ratio(i, j) else {
            return false;
        };
        let pair_idx = i.min(j);
        let Some(stats) = self.pair_stats.get_mut(pair_idx) else {
            return false;
        };
        stats.0 += 1;
        let accept = delta >= 0.0 || rng.random_range(0.0..1.0) < delta.exp();
        if accept {
            stats.1 += 1;
            self.replica_index.swap(i, j);
            self.exchange_configurations(i, j);
        }
        accept
    }

    /// Bookkeeping after an accepted swap: slots `i` and `j` exchange their
    /// configurations, so the energy data follows them.
    ///
    /// Temperature REMD swaps the stored energies.  Hamiltonian REMD swaps the
    /// columns of the cross-energy matrix (column `l` belongs to the
    /// configuration of slot `l`) and re-reads the diagonal, because the
    /// configuration that moved to slot `i` now feels the Hamiltonian of slot
    /// `i`.
    fn exchange_configurations(&mut self, i: usize, j: usize) {
        match self.cross_energies.as_mut() {
            Some(u) if self.lambda.is_some() => {
                for row in u.iter_mut() {
                    row.swap(i, j);
                }
                for k in [i, j] {
                    if let Some(&diagonal) = u.get(k).and_then(|row| row.get(k)) {
                        self.energies[k] = diagonal;
                    }
                }
            }
            _ => self.energies.swap(i, j),
        }
    }

    /// Acceptance ratio for the pair `(i, i+1)`.
    pub fn acceptance_ratio(&self, pair_index: usize) -> f64 {
        let (att, acc) = self.pair_stats[pair_index];
        if att == 0 {
            0.0
        } else {
            acc as f64 / att as f64
        }
    }

    /// Overall acceptance ratio across all pairs.
    pub fn overall_acceptance(&self) -> f64 {
        let (total_att, total_acc): (u64, u64) = self
            .pair_stats
            .iter()
            .fold((0, 0), |(a, b), &(x, y)| (a + x, b + y));
        if total_att == 0 {
            0.0
        } else {
            total_acc as f64 / total_att as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── RemdController ──────────────────────────────────────────────────

    #[test]
    fn test_remd_geometric_ladder() {
        let remd = RemdController::new_geometric(300.0, 600.0, 4, 100);
        assert_eq!(remd.n_replicas(), 4);
        assert!((remd.temperatures[0] - 300.0).abs() < 1e-6);
        assert!((remd.temperatures[3] - 600.0).abs() < 1.0);
        // Check geometric spacing
        let ratio = remd.temperatures[1] / remd.temperatures[0];
        let ratio2 = remd.temperatures[2] / remd.temperatures[1];
        assert!((ratio - ratio2).abs() < 1e-6);
    }

    #[test]
    fn test_remd_swap_always_accept_downhill() {
        let mut remd = RemdController::new_geometric(300.0, 400.0, 2, 1);
        // The cold replica (slot 0) holds the higher energy: moving the
        // low-energy configuration to the cold replica is always favourable.
        remd.set_energies(&[100.0, 50.0]);
        let swapped = remd.attempt_swaps();
        assert_eq!(swapped, vec![(0, 1)]);
        assert_eq!(remd.pair_stats[0], (1, 1));
        assert_eq!(remd.energies, vec![50.0, 100.0]);
        assert_eq!(remd.replica_index, vec![1, 0]);
    }

    /// Reset a two-replica controller to the given slot energies.
    fn two_replica_remd(t_cold: f64, t_hot: f64, e_cold: f64, e_hot: f64) -> RemdController {
        let mut remd = RemdController::new_geometric(t_cold, t_hot, 2, 1);
        remd.set_energies(&[e_cold, e_hot]);
        remd
    }

    #[test]
    fn test_remd_metropolis_direction_is_not_inverted() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let mut rng = StdRng::seed_from_u64(3419);

        // Cold replica at -100 kJ/mol, hot replica at -200 kJ/mol: delta =
        // (beta_c - beta_h)(E_c - E_h) = +10, so every attempt is accepted.
        let mut accepted = 0;
        for _ in 0..10_000 {
            let mut remd = two_replica_remd(300.0, 400.0, -100.0, -200.0);
            assert!(remd.swap_log_ratio(0, 1).unwrap() > 9.0);
            if remd.try_swap(0, 1, &mut rng) {
                accepted += 1;
            }
        }
        assert_eq!(accepted, 10_000, "favourable swaps must always be accepted");

        // Reversed energies: delta = -10, acceptance exp(-10) = 4.5e-5.
        let mut accepted = 0;
        for _ in 0..10_000 {
            let mut remd = two_replica_remd(300.0, 400.0, -200.0, -100.0);
            assert!(remd.swap_log_ratio(0, 1).unwrap() < -9.0);
            if remd.try_swap(0, 1, &mut rng) {
                accepted += 1;
            }
        }
        assert!(
            accepted < 10,
            "unfavourable swaps accepted {accepted}/10000 times (expected < 0.1 %)"
        );
    }

    #[test]
    fn test_remd_acceptance_frequency_matches_exp_delta() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let mut rng = StdRng::seed_from_u64(1_000_003);
        // beta(300 K) - beta(330 K) ~ 0.0364 mol/kJ; E_c - E_h = -19.2 kJ/mol
        // gives delta ~ -0.70 and an acceptance of exp(delta) ~ 0.50.
        let trials = 100_000;
        let probe = two_replica_remd(300.0, 330.0, -50.0, -30.8);
        let delta = probe.swap_log_ratio(0, 1).unwrap();
        let expected = delta.exp();
        assert!(delta < -0.6 && delta > -0.8, "delta = {delta}");
        let mut accepted = 0_u32;
        for _ in 0..trials {
            let mut remd = two_replica_remd(300.0, 330.0, -50.0, -30.8);
            if remd.try_swap(0, 1, &mut rng) {
                accepted += 1;
            }
        }
        let frequency = f64::from(accepted) / f64::from(trials);
        assert!(
            (frequency - expected).abs() < 0.01,
            "acceptance {frequency} vs exp(delta) = {expected}"
        );
    }

    #[test]
    fn test_remd_swap_log_ratio_matches_the_closed_form() {
        let remd = two_replica_remd(300.0, 400.0, -100.0, -200.0);
        let beta_c = 1.0 / (KB_KJMOL * 300.0);
        let beta_h = 1.0 / (KB_KJMOL * 400.0);
        let expected = (beta_c - beta_h) * (-100.0 - -200.0);
        let delta = remd.swap_log_ratio(0, 1).unwrap();
        assert!((delta - expected).abs() < 1e-12, "{delta} vs {expected}");
        // The reverse proposal has the opposite log ratio.
        let reverse = remd.swap_log_ratio(1, 0).unwrap();
        assert!((delta + reverse).abs() < 1e-12);
        assert!(remd.swap_log_ratio(0, 5).is_none());
    }

    #[test]
    fn test_remd_hamiltonian() {
        let remd = RemdController::new_hamiltonian(300.0, vec![0.0, 0.5, 1.0], 50);
        assert_eq!(remd.n_replicas(), 3);
        assert!(remd.lambda.is_some());
        assert!(remd.cross_energies.is_none());
    }

    /// Two harmonic Hamiltonians `U_k(x) = kappa_k x^2 / 2` (kappa = 1 and 10
    /// kJ mol^-1 per unit^2) and the configurations held by slots 0 and 1.
    fn harmonic_cross_energies(x0: f64, x1: f64) -> Vec<Vec<f64>> {
        let kappa = [1.0_f64, 10.0];
        let x = [x0, x1];
        (0..2)
            .map(|k| (0..2).map(|l| 0.5 * kappa[k] * x[l] * x[l]).collect())
            .collect()
    }

    #[test]
    fn test_hremd_swap_log_ratio_uses_the_cross_energies() {
        let mut remd = RemdController::new_hamiltonian(300.0, vec![0.0, 1.0], 1);
        assert!(remd.swap_log_ratio(0, 1).is_none(), "no cross-energies yet");
        remd.set_cross_energies(harmonic_cross_energies(0.7, 0.3))
            .unwrap();
        // delta = -beta [(u01 + u10) - (u00 + u11)] with u00 = 0.245, u01 = 0.045,
        // u10 = 2.45, u11 = 0.45.
        let beta = 1.0 / (KB_KJMOL * 300.0);
        let expected = -beta * ((0.045 + 2.45) - (0.245 + 0.45));
        let delta = remd.swap_log_ratio(0, 1).unwrap();
        assert!((delta - expected).abs() < 1e-12, "{delta} vs {expected}");
        assert!(delta < 0.0 && delta.is_finite());
        // The diagonal refreshed the per-slot energies.
        assert!((remd.energies[0] - 0.245).abs() < 1e-12);
        assert!((remd.energies[1] - 0.45).abs() < 1e-12);
    }

    #[test]
    fn test_hremd_acceptance_frequency_matches_exp_delta() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let mut rng = StdRng::seed_from_u64(3418);
        let trials = 100_000;
        let u = harmonic_cross_energies(0.7, 0.3);
        let mut remd = RemdController::new_hamiltonian(300.0, vec![0.0, 1.0], 1);
        remd.set_cross_energies(u.clone()).unwrap();
        let expected = remd.swap_log_ratio(0, 1).unwrap().exp();
        assert!(expected > 0.3 && expected < 0.7, "exp(delta) = {expected}");

        let mut accepted = 0_u32;
        for _ in 0..trials {
            // Restore the pre-swap configuration assignment for every trial.
            remd.set_cross_energies(u.clone()).unwrap();
            if remd.try_swap(0, 1, &mut rng) {
                accepted += 1;
            }
        }
        let frequency = f64::from(accepted) / f64::from(trials);
        assert!(
            (frequency - expected).abs() < 0.01,
            "acceptance {frequency} vs exp(delta) = {expected}"
        );
        // Without the lambda-dependent criterion every swap would be accepted
        // (beta_i == beta_j gives delta = 0 from the temperatures alone).
        assert!(frequency < 0.9);
    }

    #[test]
    fn test_hremd_without_cross_energies_rejects_without_counting() {
        let mut remd = RemdController::new_hamiltonian(300.0, vec![0.0, 0.5, 1.0], 1);
        for _ in 0..4 {
            assert!(remd.attempt_swaps().is_empty());
        }
        assert!(remd.pair_stats.iter().all(|&s| s == (0, 0)));
        assert_eq!(remd.replica_index, vec![0, 1, 2]);
    }

    #[test]
    fn test_hremd_accepted_swap_exchanges_configurations_consistently() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let mut rng = StdRng::seed_from_u64(7);
        let mut remd = RemdController::new_hamiltonian(300.0, vec![0.0, 1.0], 1);
        // Strongly favourable: cross energies far above the diagonal are
        // the *un*favourable case, so make the diagonal the high one.
        let u = vec![vec![5.0, 1.0], vec![2.0, 6.0]];
        remd.set_cross_energies(u).unwrap();
        let before = remd.swap_log_ratio(0, 1).unwrap();
        assert!(before > 0.0, "delta = {before}");
        assert!(
            remd.try_swap(0, 1, &mut rng),
            "delta > 0 is always accepted"
        );
        assert_eq!(remd.replica_index, vec![1, 0]);
        // Columns exchanged: slot 0 now holds the configuration with u(0, .) = 1.0.
        let swapped = remd.cross_energies.as_ref().unwrap();
        assert_eq!(swapped[0], vec![1.0, 5.0]);
        assert_eq!(swapped[1], vec![6.0, 2.0]);
        assert_eq!(remd.energies, vec![1.0, 2.0]);
        // Detailed balance: the reverse move has exactly the opposite log ratio.
        let after = remd.swap_log_ratio(0, 1).unwrap();
        assert!((before + after).abs() < 1e-12, "{before} vs {after}");
        assert_eq!(remd.pair_stats[0], (1, 1));
    }

    #[test]
    fn test_hremd_set_cross_energies_validates_its_input() {
        let mut remd = RemdController::new_hamiltonian(300.0, vec![0.0, 0.5, 1.0], 1);
        // Wrong number of rows, wrong row length, NaN and -inf entries.
        assert!(remd.set_cross_energies(vec![vec![0.0; 3]; 2]).is_err());
        assert!(
            remd.set_cross_energies(vec![vec![0.0; 3], vec![0.0; 2], vec![0.0; 3]])
                .is_err()
        );
        let mut nan = vec![vec![0.0; 3]; 3];
        nan[1][2] = f64::NAN;
        assert!(remd.set_cross_energies(nan).is_err());
        let mut neg_inf = vec![vec![0.0; 3]; 3];
        neg_inf[0][1] = f64::NEG_INFINITY;
        assert!(remd.set_cross_energies(neg_inf).is_err());
        assert!(remd.cross_energies.is_none(), "failed calls change nothing");
        // +inf marks a forbidden overlap and makes the swap impossible.
        let mut forbidden = vec![vec![1.0; 3]; 3];
        forbidden[0][1] = f64::INFINITY;
        remd.set_cross_energies(forbidden).unwrap();
        assert_eq!(remd.swap_log_ratio(0, 1), Some(f64::NEG_INFINITY));
        // Temperature REMD controllers have no use for cross-energies.
        let mut trem = RemdController::new_geometric(300.0, 400.0, 2, 1);
        assert!(trem.set_cross_energies(vec![vec![0.0; 2]; 2]).is_err());
    }

    #[test]
    fn test_remd_acceptance_initially_zero() {
        let remd = RemdController::new_geometric(300.0, 400.0, 3, 10);
        assert!((remd.overall_acceptance() - 0.0).abs() < 1e-12);
    }
}
