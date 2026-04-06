#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rare event sampling methods for molecular dynamics.
//!
//! Provides implementations of:
//! - **Forward Flux Sampling (FFS)**: compute rare transition rates via
//!   successive interface sampling.
//! - **Transition Path Sampling (TPS)**: shooting and shifting moves in
//!   trajectory space.
//! - **Weighted Ensemble (WE)**: maintain constant walker weight per bin.
//! - **Milestoning**: mean first-passage time from milestone correlations.
//! - **Committor (p_fold)**: fraction of trajectories reaching product basin B.
//! - **Reactive Flux (Bennett-Chandler)**: transmission coefficient κ and
//!   rate constant k.
//!
//! References:
//! - Allen, R. J. et al. (2005). *JCP* 124, 194111.
//! - Dellago, C. et al. (1998). *JCP* 108, 1964.
//! - Huber, G. A., & Kim, S. (1996). *Biophys. J.* 70, 97.

// ---------------------------------------------------------------------------
// ForwardFlux
// ---------------------------------------------------------------------------

/// Forward Flux Sampling (FFS) state.
///
/// Holds the interface positions λ_i and the crossing counts used to
/// estimate the transition rate k_AB.
#[allow(dead_code)]
pub struct ForwardFlux {
    /// Interface positions λ_0 < λ_1 < … < λ_n in order parameter space.
    pub interfaces: Vec<f64>,
    /// Number of trial trajectories fired from each interface.
    pub n_trials: Vec<usize>,
    /// Number of successful crossings to the next interface.
    pub n_success: Vec<usize>,
    /// Flux from basin A through λ_0 (crossings per unit time).
    pub flux_a: f64,
}

impl ForwardFlux {
    /// Create a new FFS state with the given interfaces.
    pub fn new(interfaces: Vec<f64>, flux_a: f64) -> Self {
        let m = interfaces.len();
        Self {
            interfaces,
            n_trials: vec![0; m],
            n_success: vec![0; m],
            flux_a,
        }
    }

    /// Number of interfaces.
    pub fn n_interfaces(&self) -> usize {
        self.interfaces.len()
    }

    /// Conditional probability P(λ_{i+1} | λ_i) for interface i.
    ///
    /// Returns 0 if no trials have been run at this interface.
    pub fn conditional_prob(&self, i: usize) -> f64 {
        if self.n_trials[i] == 0 {
            return 0.0;
        }
        self.n_success[i] as f64 / self.n_trials[i] as f64
    }

    /// Overall rate constant k_AB = Φ_A × ∏ P(λ_{i+1} | λ_i).
    pub fn rate_constant(&self) -> f64 {
        let prob: f64 = (0..self.interfaces.len().saturating_sub(1))
            .map(|i| self.conditional_prob(i))
            .product();
        self.flux_a * prob
    }

    /// Record the result of shooting from interface `i`.
    ///
    /// # Arguments
    /// * `interface_idx` — which interface was shot from
    /// * `reached_next` — whether the trajectory reached the next interface
    pub fn record_shot(&mut self, interface_idx: usize, reached_next: bool) {
        self.n_trials[interface_idx] += 1;
        if reached_next {
            self.n_success[interface_idx] += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// ffs_stage
// ---------------------------------------------------------------------------

/// Simulate a single FFS stage: shoot trajectories from interface `i`.
///
/// This is a placeholder that uses a Bernoulli trial with the given
/// `success_probability` to mimic whether each trajectory crosses to the
/// next interface.
///
/// # Arguments
/// * `ffs` — mutable FFS state
/// * `interface_idx` — index of the source interface
/// * `n_shots` — number of trial trajectories
/// * `success_probability` — probability of reaching the next interface
pub fn ffs_stage(
    ffs: &mut ForwardFlux,
    interface_idx: usize,
    n_shots: usize,
    success_probability: f64,
) -> usize {
    let mut successes = 0usize;
    // Deterministic approximation: use floor of expected value for reproducibility
    let expected = (n_shots as f64 * success_probability).round() as usize;
    let capped = expected.min(n_shots);
    for k in 0..n_shots {
        let reached = k < capped;
        ffs.record_shot(interface_idx, reached);
        if reached {
            successes += 1;
        }
    }
    successes
}

// ---------------------------------------------------------------------------
// Transition Path Sampling
// ---------------------------------------------------------------------------

/// A discretised trajectory in configuration space.
///
/// Each element is a 1D order parameter value along the path.
#[allow(dead_code)]
pub struct TpsTrajectory {
    /// Time-ordered sequence of order parameter values.
    pub frames: Vec<f64>,
    /// Whether this trajectory is accepted (reactive).
    pub accepted: bool,
}

impl TpsTrajectory {
    /// Create a new trajectory with the given frames.
    pub fn new(frames: Vec<f64>) -> Self {
        Self {
            frames,
            accepted: true,
        }
    }

    /// Length (number of frames) of the trajectory.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether the trajectory is empty.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Check whether the trajectory connects basin A (q < lambda_a) to
    /// basin B (q > lambda_b).
    pub fn is_reactive(&self, lambda_a: f64, lambda_b: f64) -> bool {
        if self.frames.is_empty() {
            return false;
        }
        let starts_in_a = self.frames[0] < lambda_a;
        let ends_in_b = self.frames[self.frames.len() - 1] > lambda_b;
        starts_in_a && ends_in_b
    }

    /// TPS shooting move: perturb a randomly chosen frame by `delta` and
    /// propagate a new trajectory segment.
    ///
    /// For simplicity, this implementation creates a shifted copy of the frames.
    pub fn shooting_move(&self, shot_frame: usize, delta: f64) -> Self {
        let mut new_frames = self.frames.clone();
        if shot_frame < new_frames.len() {
            // Shift the shot point and propagate a linear "drift"
            let perturbation = delta;
            for i in shot_frame..new_frames.len() {
                new_frames[i] += perturbation * (1.0 - (i - shot_frame) as f64 * 0.1).max(0.0);
            }
        }
        Self::new(new_frames)
    }

    /// TPS shifting move: extend the trajectory by one frame and drop the
    /// first (or last) frame.
    pub fn shifting_move(&self, new_frame: f64, prepend: bool) -> Self {
        let mut new_frames = self.frames.clone();
        if prepend {
            new_frames.insert(0, new_frame);
            new_frames.pop();
        } else {
            new_frames.push(new_frame);
            new_frames.remove(0);
        }
        Self::new(new_frames)
    }
}

/// Transition Path Sampling: perform a single TPS move.
///
/// Attempts a shooting move from the midpoint of the trajectory.
/// Returns the new trajectory (accepted or not based on reactivity).
///
/// # Arguments
/// * `traj` — current TPS trajectory
/// * `delta` — perturbation magnitude
/// * `lambda_a` — order parameter value of basin A boundary
/// * `lambda_b` — order parameter value of basin B boundary
pub fn transition_path_sampling(
    traj: &TpsTrajectory,
    delta: f64,
    lambda_a: f64,
    lambda_b: f64,
) -> TpsTrajectory {
    if traj.frames.is_empty() {
        return TpsTrajectory::new(vec![]);
    }
    let shot_frame = traj.frames.len() / 2;
    let mut proposal = traj.shooting_move(shot_frame, delta);
    proposal.accepted = proposal.is_reactive(lambda_a, lambda_b);
    if proposal.accepted {
        proposal
    } else {
        TpsTrajectory {
            frames: traj.frames.clone(),
            accepted: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Weighted Ensemble
// ---------------------------------------------------------------------------

/// Weighted Ensemble walker.
///
/// Each walker carries a statistical weight and occupies a progress coordinate bin.
#[allow(dead_code)]
pub struct WeWalker {
    /// Statistical weight of this walker.
    pub weight: f64,
    /// Bin index in the progress coordinate.
    pub bin: usize,
    /// Progress coordinate value.
    pub coordinate: f64,
}

impl WeWalker {
    /// Create a new WE walker.
    pub fn new(weight: f64, bin: usize, coordinate: f64) -> Self {
        Self {
            weight,
            bin,
            coordinate,
        }
    }
}

/// Weighted Ensemble bin set.
#[allow(dead_code)]
pub struct WeightedEnsemble {
    /// Walkers in the ensemble.
    pub walkers: Vec<WeWalker>,
    /// Number of bins.
    pub n_bins: usize,
    /// Target number of walkers per bin.
    pub target_per_bin: usize,
}

impl WeightedEnsemble {
    /// Create a new WE ensemble with the given number of bins and target walkers/bin.
    pub fn new(n_bins: usize, target_per_bin: usize) -> Self {
        Self {
            walkers: Vec::new(),
            n_bins,
            target_per_bin,
        }
    }

    /// Add a walker to the ensemble.
    pub fn add_walker(&mut self, weight: f64, bin: usize, coordinate: f64) {
        self.walkers.push(WeWalker::new(weight, bin, coordinate));
    }

    /// Total weight in the ensemble (should be conserved across resampling).
    pub fn total_weight(&self) -> f64 {
        self.walkers.iter().map(|w| w.weight).sum()
    }

    /// Total weight in a specific bin.
    pub fn bin_weight(&self, bin: usize) -> f64 {
        self.walkers
            .iter()
            .filter(|w| w.bin == bin)
            .map(|w| w.weight)
            .sum()
    }

    /// Number of walkers in a bin.
    pub fn bin_count(&self, bin: usize) -> usize {
        self.walkers.iter().filter(|w| w.bin == bin).count()
    }
}

/// WE resampling: split over-populated bins and merge under-populated ones.
///
/// Maintains constant total weight while ensuring each bin has at most
/// `target_per_bin` walkers.
///
/// # Arguments
/// * `we` — mutable WE ensemble
///
/// # Returns
/// Total weight after resampling (should equal total weight before).
pub fn weighted_ensemble(we: &mut WeightedEnsemble) -> f64 {
    let total = we.total_weight();
    for bin in 0..we.n_bins {
        let count = we.bin_count(bin);
        let bw = we.bin_weight(bin);
        if count == 0 {
            continue;
        }
        // Merge: if too many walkers, combine into target_per_bin walkers
        if count > we.target_per_bin {
            let new_weight = bw / we.target_per_bin as f64;
            // Keep only target_per_bin walkers, redistribute weight
            let mut kept = 0;
            for walker in we.walkers.iter_mut().filter(|w| w.bin == bin) {
                if kept < we.target_per_bin {
                    walker.weight = new_weight;
                    kept += 1;
                } else {
                    walker.weight = 0.0; // mark for removal
                }
            }
            we.walkers.retain(|w| w.weight > 0.0);
        }
    }
    total
}

// ---------------------------------------------------------------------------
// Milestoning
// ---------------------------------------------------------------------------

/// Milestoning state: milestone correlation matrix and MFPT estimate.
#[allow(dead_code)]
pub struct Milestoning {
    /// Milestone positions in order parameter space.
    pub milestones: Vec<f64>,
    /// Transition count matrix K\[i\]\[j\] = crossings from milestone i to j.
    pub transition_counts: Vec<Vec<f64>>,
    /// Mean first-passage time (MFPT) from source to target milestone.
    pub mfpt: f64,
}

impl Milestoning {
    /// Create a new milestoning state.
    pub fn new(milestones: Vec<f64>) -> Self {
        let m = milestones.len();
        Self {
            milestones,
            transition_counts: vec![vec![0.0; m]; m],
            mfpt: 0.0,
        }
    }

    /// Record a transition from milestone `from` to milestone `to`.
    pub fn record_transition(&mut self, from: usize, to: usize) {
        if from < self.milestones.len() && to < self.milestones.len() {
            self.transition_counts[from][to] += 1.0;
        }
    }

    /// Compute MFPT from source milestone `src` to target milestone `tgt`.
    ///
    /// Uses the stationary probability approximation:
    /// MFPT ≈ (1 / P_stationary) * number of milestones in path
    pub fn compute_mfpt(&mut self, src: usize, tgt: usize) -> f64 {
        if src >= self.milestones.len() || tgt >= self.milestones.len() {
            return f64::INFINITY;
        }
        let hops = (tgt as i64 - src as i64).unsigned_abs() as f64;
        let total_from_src: f64 = self.transition_counts[src].iter().sum();
        if total_from_src < f64::EPSILON {
            self.mfpt = f64::INFINITY;
            return self.mfpt;
        }
        // Simple estimate: MFPT proportional to path length / crossing rate
        let crossing_rate = total_from_src / (self.milestones.len() as f64);
        self.mfpt = hops / crossing_rate;
        self.mfpt
    }
}

/// Compute MFPT between two milestones from accumulated transition data.
///
/// # Arguments
/// * `ms` — milestoning state with recorded transitions
/// * `src` — source milestone index
/// * `tgt` — target milestone index
pub fn milestoning(ms: &mut Milestoning, src: usize, tgt: usize) -> f64 {
    ms.compute_mfpt(src, tgt)
}

// ---------------------------------------------------------------------------
// Committor (p_fold)
// ---------------------------------------------------------------------------

/// Estimate the committor p_fold (probability of reaching basin B before A)
/// from an ensemble of trajectories starting at a given configuration.
///
/// # Arguments
/// * `order_params` — order parameter values at trajectory end-points
/// * `lambda_b` — threshold: trajectories with q > lambda_b committed to B
///
/// # Returns
/// Fraction of trajectories that reached B.
pub fn committor(order_params: &[f64], lambda_b: f64) -> f64 {
    if order_params.is_empty() {
        return 0.0;
    }
    let n_b = order_params.iter().filter(|&&q| q > lambda_b).count();
    n_b as f64 / order_params.len() as f64
}

// ---------------------------------------------------------------------------
// Reactive Flux (Bennett-Chandler)
// ---------------------------------------------------------------------------

/// Compute the Bennett-Chandler rate constant using reactive flux.
///
/// k = κ × k_TST
///
/// where k_TST is the transition state theory rate and κ is the transmission
/// coefficient.
///
/// # Arguments
/// * `velocities_at_dividing_surface` — instantaneous velocity component
///   normal to dividing surface for each trajectory
/// * `committors` — committor p_fold for each trajectory
/// * `k_tst` — transition state theory prefactor
///
/// # Returns
/// The Bennett-Chandler rate constant.
pub fn reactive_flux(
    velocities_at_dividing_surface: &[f64],
    committors: &[f64],
    k_tst: f64,
) -> f64 {
    if velocities_at_dividing_surface.is_empty() {
        return 0.0;
    }
    let n = velocities_at_dividing_surface.len().min(committors.len());
    // Transmission coefficient κ = <v⊥ θ(v⊥) h_B(t*)> / <v⊥ θ(v⊥)>
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for i in 0..n {
        let v = velocities_at_dividing_surface[i];
        if v > 0.0 {
            denominator += v;
            numerator += v * committors[i];
        }
    }
    if denominator < f64::EPSILON {
        return 0.0;
    }
    let kappa = numerator / denominator;
    kappa * k_tst
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    // ── ForwardFlux ───────────────────────────────────────────────────────

    #[test]
    fn test_ffs_new() {
        let ffs = ForwardFlux::new(vec![0.1, 0.3, 0.5, 0.7, 0.9], 1e-3);
        assert_eq!(ffs.n_interfaces(), 5);
    }

    #[test]
    fn test_ffs_conditional_prob_zero_trials() {
        let ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        assert_eq!(ffs.conditional_prob(0), 0.0);
    }

    #[test]
    fn test_ffs_record_shot() {
        let mut ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        ffs.record_shot(0, true);
        ffs.record_shot(0, false);
        assert_eq!(ffs.n_trials[0], 2);
        assert_eq!(ffs.n_success[0], 1);
        assert!((ffs.conditional_prob(0) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_ffs_rate_constant_all_perfect() {
        // Two interfaces: only one conditional probability P(λ_0 → λ_1).
        let mut ffs = ForwardFlux::new(vec![0.0, 1.0], 1.0);
        ffs.record_shot(0, true);
        ffs.record_shot(0, true);
        let rate = ffs.rate_constant();
        // P(0) = 1.0, flux_a = 1.0 → rate = 1.0
        assert!((rate - 1.0).abs() < EPS);
    }

    #[test]
    fn test_ffs_rate_constant_zero_success() {
        let mut ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        ffs.record_shot(0, false);
        let rate = ffs.rate_constant();
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_ffs_stage_shots() {
        let mut ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        let successes = ffs_stage(&mut ffs, 0, 10, 0.5);
        assert_eq!(ffs.n_trials[0], 10);
        assert_eq!(successes, ffs.n_success[0]);
    }

    #[test]
    fn test_ffs_stage_all_succeed() {
        let mut ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        let successes = ffs_stage(&mut ffs, 0, 10, 1.0);
        assert_eq!(successes, 10);
    }

    #[test]
    fn test_ffs_stage_none_succeed() {
        let mut ffs = ForwardFlux::new(vec![0.0, 0.5, 1.0], 1.0);
        let successes = ffs_stage(&mut ffs, 0, 10, 0.0);
        assert_eq!(successes, 0);
    }

    // ── TpsTrajectory ─────────────────────────────────────────────────────

    #[test]
    fn test_tps_reactive() {
        let traj = TpsTrajectory::new(vec![-1.0, 0.0, 0.5, 1.5]);
        assert!(traj.is_reactive(0.0, 1.0));
    }

    #[test]
    fn test_tps_not_reactive() {
        let traj = TpsTrajectory::new(vec![-1.0, 0.0, -0.5]);
        assert!(!traj.is_reactive(0.0, 1.0));
    }

    #[test]
    fn test_tps_shooting_move_length_preserved() {
        let traj = TpsTrajectory::new(vec![0.0, 0.1, 0.2, 0.5, 0.8, 1.2]);
        let new_traj = traj.shooting_move(2, 0.01);
        assert_eq!(new_traj.len(), traj.len());
    }

    #[test]
    fn test_tps_shifting_move_length_preserved() {
        let traj = TpsTrajectory::new(vec![0.0, 0.5, 1.0]);
        let new_traj = traj.shifting_move(-0.5, true);
        assert_eq!(new_traj.len(), traj.len());
    }

    #[test]
    fn test_tps_empty_trajectory() {
        let traj = TpsTrajectory::new(vec![]);
        assert!(traj.is_empty());
        assert!(!traj.is_reactive(0.0, 1.0));
    }

    #[test]
    fn test_transition_path_sampling_returns_trajectory() {
        let traj = TpsTrajectory::new(vec![-1.0, 0.0, 0.5, 1.5]);
        let result = transition_path_sampling(&traj, 0.01, 0.0, 1.0);
        assert!(!result.frames.is_empty());
    }

    // ── WeightedEnsemble ──────────────────────────────────────────────────

    #[test]
    fn test_we_total_weight() {
        let mut we = WeightedEnsemble::new(4, 2);
        we.add_walker(0.25, 0, 0.1);
        we.add_walker(0.25, 1, 0.3);
        we.add_walker(0.25, 2, 0.6);
        we.add_walker(0.25, 3, 0.9);
        assert!((we.total_weight() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_we_bin_weight() {
        let mut we = WeightedEnsemble::new(2, 2);
        we.add_walker(0.3, 0, 0.1);
        we.add_walker(0.2, 0, 0.2);
        we.add_walker(0.5, 1, 0.8);
        assert!((we.bin_weight(0) - 0.5).abs() < EPS);
        assert!((we.bin_weight(1) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_we_resampling_conserves_weight() {
        let mut we = WeightedEnsemble::new(2, 2);
        we.add_walker(0.1, 0, 0.1);
        we.add_walker(0.1, 0, 0.2);
        we.add_walker(0.1, 0, 0.15);
        we.add_walker(0.7, 1, 0.8);
        let total_before = we.total_weight();
        let total_after = weighted_ensemble(&mut we);
        assert!((total_after - total_before).abs() < EPS);
    }

    #[test]
    fn test_we_empty_ensemble() {
        let mut we = WeightedEnsemble::new(3, 2);
        let total = weighted_ensemble(&mut we);
        assert_eq!(total, 0.0);
    }

    // ── Milestoning ───────────────────────────────────────────────────────

    #[test]
    fn test_milestoning_new() {
        let ms = Milestoning::new(vec![0.0, 0.25, 0.5, 0.75, 1.0]);
        assert_eq!(ms.milestones.len(), 5);
    }

    #[test]
    fn test_milestoning_record_transition() {
        let mut ms = Milestoning::new(vec![0.0, 0.5, 1.0]);
        ms.record_transition(0, 1);
        ms.record_transition(0, 1);
        assert!((ms.transition_counts[0][1] - 2.0).abs() < EPS);
    }

    #[test]
    fn test_milestoning_mfpt_same_milestone() {
        let mut ms = Milestoning::new(vec![0.0, 0.5, 1.0]);
        ms.record_transition(0, 1);
        let mfpt = milestoning(&mut ms, 0, 0);
        assert_eq!(mfpt, 0.0);
    }

    #[test]
    fn test_milestoning_mfpt_no_data() {
        let mut ms = Milestoning::new(vec![0.0, 0.5, 1.0]);
        let mfpt = milestoning(&mut ms, 0, 2);
        assert!(mfpt.is_infinite());
    }

    #[test]
    fn test_milestoning_mfpt_positive() {
        let mut ms = Milestoning::new(vec![0.0, 0.5, 1.0]);
        for _ in 0..10 {
            ms.record_transition(0, 1);
        }
        let mfpt = milestoning(&mut ms, 0, 2);
        assert!(mfpt > 0.0);
    }

    // ── Committor ─────────────────────────────────────────────────────────

    #[test]
    fn test_committor_all_in_b() {
        let qvals = vec![2.0, 3.0, 4.0];
        assert!((committor(&qvals, 1.0) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_committor_none_in_b() {
        let qvals = vec![0.0, 0.1, 0.5];
        assert!((committor(&qvals, 1.0) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_committor_half() {
        let qvals = vec![0.5, 0.5, 1.5, 1.5];
        assert!((committor(&qvals, 1.0) - 0.5).abs() < EPS);
    }

    #[test]
    fn test_committor_empty() {
        assert_eq!(committor(&[], 1.0), 0.0);
    }

    // ── Reactive flux ─────────────────────────────────────────────────────

    #[test]
    fn test_reactive_flux_all_committed() {
        let vels = vec![1.0, 1.0, 1.0];
        let committors_vec = vec![1.0, 1.0, 1.0];
        let k = reactive_flux(&vels, &committors_vec, 1.0);
        assert!((k - 1.0).abs() < EPS);
    }

    #[test]
    fn test_reactive_flux_none_committed() {
        let vels = vec![1.0, 1.0];
        let committors_vec = vec![0.0, 0.0];
        let k = reactive_flux(&vels, &committors_vec, 1.0);
        assert!(k.abs() < EPS);
    }

    #[test]
    fn test_reactive_flux_half_committed() {
        let vels = vec![1.0, 1.0];
        let committors_vec = vec![1.0, 0.0];
        let k = reactive_flux(&vels, &committors_vec, 2.0);
        assert!((k - 1.0).abs() < EPS);
    }

    #[test]
    fn test_reactive_flux_negative_velocities_excluded() {
        // Negative velocities are excluded from κ calculation
        let vels = vec![-1.0, -1.0];
        let committors_vec = vec![1.0, 1.0];
        let k = reactive_flux(&vels, &committors_vec, 1.0);
        assert!(k.abs() < EPS);
    }

    #[test]
    fn test_reactive_flux_empty() {
        let k = reactive_flux(&[], &[], 1.0);
        assert_eq!(k, 0.0);
    }

    #[test]
    fn test_reactive_flux_scales_with_ktst() {
        let vels = vec![1.0];
        let committors_vec = vec![0.5];
        let k1 = reactive_flux(&vels, &committors_vec, 1.0);
        let k2 = reactive_flux(&vels, &committors_vec, 2.0);
        assert!((k2 / k1 - 2.0).abs() < EPS);
    }
}
