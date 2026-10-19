//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{GaussianHill, MetadynamicsState};
use crate::free_energy::bar::log_mean_fermi;

/// Boltzmann constant in kJ/mol/K
pub(super) const KB: f64 = 8.314e-3;
/// Trait for a collective variable (CV) in metadynamics.
pub trait CollectiveVariable {
    /// Compute the current CV value from atomic positions.
    fn value(&self, positions: &[[f64; 3]]) -> f64;
    /// Compute the gradient dCV/dr for each atom (returns Vec of length n_atoms).
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]>;
    /// Human-readable name for this CV.
    fn name(&self) -> &str;
}
/// Reconstruct a free energy surface (FES) on a 1D grid from metadynamics hills.
///
/// The FES estimate is `F(s) = -V_bias(s)` (valid as hills converge for
/// plain metadynamics; for well-tempered, multiply by (T + ΔT)/ΔT).
///
/// # Arguments
/// * `state`     – converged (or ongoing) metadynamics state.
/// * `s_min`     – lower bound of the CV grid.
/// * `s_max`     – upper bound of the CV grid.
/// * `n_points`  – number of grid points.
/// * `gamma`     – bias factor `(T + ΔT)/T` for well-tempered (use 1.0 for plain meta).
///
/// Returns a vector of `(cv_value, fes_value)` pairs.
pub fn reconstruct_fes_1d(
    state: &MetadynamicsState,
    s_min: f64,
    s_max: f64,
    n_points: usize,
    gamma: f64,
) -> Vec<(f64, f64)> {
    if n_points == 0 || s_min >= s_max {
        return Vec::new();
    }
    let ds = (s_max - s_min) / (n_points - 1).max(1) as f64;
    let fes: Vec<(f64, f64)> = (0..n_points)
        .map(|i| {
            let s = s_min + i as f64 * ds;
            let v_bias = state.bias_potential(&[s]);
            (s, -gamma * v_bias)
        })
        .collect();
    if let Some(&(_, min_val)) = fes
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    {
        fes.into_iter().map(|(s, f)| (s, f - min_val)).collect()
    } else {
        fes
    }
}
/// Check convergence of a metadynamics simulation by monitoring the change in
/// the free energy surface over successive block intervals.
///
/// Splits the hill history into `n_blocks` equal blocks, estimates the FES
/// from each block, and returns the maximum point-wise change between the
/// last two blocks.  A value smaller than `tolerance` indicates convergence.
pub fn check_convergence(
    state: &MetadynamicsState,
    s_min: f64,
    s_max: f64,
    n_grid: usize,
    n_blocks: usize,
    tolerance: f64,
) -> bool {
    let n_hills = state.hills.len();
    if n_hills < 2 * n_blocks {
        return false;
    }
    let block_size = n_hills / n_blocks;
    let block1_start = (n_blocks - 2) * block_size;
    let block1_end = (n_blocks - 1) * block_size;
    let block2_start = block1_end;
    let fes1 = fes_from_hill_range(state, block1_start, block1_end, s_min, s_max, n_grid);
    let fes2 = fes_from_hill_range(state, block2_start, n_hills, s_min, s_max, n_grid);
    let max_diff = fes1
        .iter()
        .zip(fes2.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    max_diff < tolerance
}
/// Compute FES contribution from a slice of hills (indexed by \[start, end)).
pub(super) fn fes_from_hill_range(
    state: &MetadynamicsState,
    start: usize,
    end: usize,
    s_min: f64,
    s_max: f64,
    n_grid: usize,
) -> Vec<f64> {
    if n_grid == 0 || s_min >= s_max {
        return Vec::new();
    }
    let ds = (s_max - s_min) / (n_grid - 1).max(1) as f64;
    (0..n_grid)
        .map(|i| {
            let s = s_min + i as f64 * ds;
            let v: f64 = state.hills[start..end.min(state.hills.len())]
                .iter()
                .map(|h| h.evaluate(&[s]))
                .sum();
            -v
        })
        .collect()
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_gaussian_hill_evaluate_peaks_at_center() {
        let hill = GaussianHill::new(vec![1.0, 2.0], 3.0, vec![0.5, 0.5]);
        let at_center = hill.evaluate(&[1.0, 2.0]);
        assert!(
            (at_center - 3.0).abs() < 1e-12,
            "Expected 3.0 at center, got {at_center}"
        );
        let away = hill.evaluate(&[2.0, 3.0]);
        assert!(away < at_center, "Value should be smaller away from center");
    }
    #[test]
    fn test_gaussian_hill_gradient_zero_at_center() {
        let hill = GaussianHill::new(vec![0.0], 1.0, vec![0.3]);
        let grad = hill.gradient_wrt_cv(&[0.0]);
        assert_eq!(grad.len(), 1);
        assert!(
            grad[0].abs() < 1e-12,
            "Gradient at center should be zero, got {}",
            grad[0]
        );
    }
    #[test]
    fn test_metadynamics_bias_increases_after_depositing() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        let v0 = state.bias_potential(&cv);
        assert_eq!(v0, 0.0, "No hills yet");
        state.maybe_deposit(&cv);
        let v1 = state.bias_potential(&cv);
        assert!(v1 > v0, "Bias should increase after depositing: v1={v1}");
    }
    #[test]
    fn test_well_tempered_height_decreases_in_visited_regions() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 1.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        for _ in 0..5 {
            state.maybe_deposit(&cv);
        }
        let first_height = state.hills[0].height;
        let last_height = state.hills[state.hills.len() - 1].height;
        assert!(
            last_height < first_height,
            "Well-tempered: last hill height {last_height} should be < first {first_height}"
        );
    }
    #[test]
    fn test_distance_cv_known_positions() {
        let cv = DistanceCV::new(0, 1, "d01");
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        let d = cv.value(&positions);
        assert!((d - 5.0).abs() < 1e-12, "Expected distance 5.0, got {d}");
    }
    #[test]
    fn test_steered_md_harmonic_bias_zero_at_target() {
        let smd = SteeredMD::new(0, 2.0, 100.0, 0.1);
        let bias = smd.harmonic_bias(2.0, 0);
        assert!(
            bias.abs() < 1e-12,
            "Bias at current target should be zero, got {bias}"
        );
        let bias5 = smd.harmonic_bias(2.5, 5);
        assert!(
            bias5.abs() < 1e-12,
            "Bias at step-5 target should be zero, got {bias5}"
        );
    }
    #[test]
    fn test_maybe_deposit_only_at_correct_stride() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 3, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        state.maybe_deposit(&cv);
        assert_eq!(state.hills.len(), 1, "Should deposit at step 0");
        assert_eq!(state.step_count, 1);
        state.maybe_deposit(&cv);
        assert_eq!(state.hills.len(), 1, "Should not deposit at step 1");
        assert_eq!(state.step_count, 2);
        state.maybe_deposit(&cv);
        assert_eq!(state.hills.len(), 1, "Should not deposit at step 2");
        assert_eq!(state.step_count, 3);
        state.maybe_deposit(&cv);
        assert_eq!(state.hills.len(), 2, "Should deposit at step 3");
        assert_eq!(state.step_count, 4);
    }
    #[test]
    fn test_distance_cv_gradient_direction() {
        let cv = DistanceCV::new(0, 1, "d");
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let grad = cv.gradient(&positions);
        assert!(
            grad[0][0] < 0.0,
            "grad[0].x should be negative, got {}",
            grad[0][0]
        );
        assert!(
            grad[1][0] > 0.0,
            "grad[1].x should be positive, got {}",
            grad[1][0]
        );
    }
    #[test]
    fn test_free_energy_estimate_returns_negated_bias() {
        let params = MetadynamicsParams::new(2.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        state.maybe_deposit(&cv);
        let bins = vec![vec![0.0f64]];
        let fe = state.free_energy_estimate(&bins);
        let bias = state.bias_potential(&[0.0]);
        assert!(
            (fe[0] + bias).abs() < 1e-12,
            "Free energy should be -bias_potential"
        );
    }
    #[test]
    fn test_multiple_walkers_n_walkers() {
        let mw = MultipleWalkersMetadynamics::new(3, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        assert_eq!(mw.n_walkers(), 3);
    }
    #[test]
    fn test_multiple_walkers_deposit_and_sync() {
        let mut mw = MultipleWalkersMetadynamics::new(2, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        mw.walkers[0].maybe_deposit(&[0.0]);
        mw.walkers[1].maybe_deposit(&[2.0]);
        let pre_sync_total = mw.total_hills();
        assert_eq!(pre_sync_total, 2, "2 hills total before sync");
        mw.synchronise();
        for w in &mw.walkers {
            assert_eq!(
                w.hills.len(),
                2,
                "Each walker should have 2 hills after sync"
            );
        }
    }
    #[test]
    fn test_multiple_walkers_repeated_sync_does_not_multiply_hills() {
        let mut mw = MultipleWalkersMetadynamics::new(3, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        for (k, walker) in mw.walkers.iter_mut().enumerate() {
            walker.maybe_deposit(&[k as f64]);
        }
        mw.synchronise();
        mw.synchronise();
        for walker in &mw.walkers {
            assert_eq!(
                walker.hills.len(),
                3,
                "a second sync must not duplicate hills"
            );
        }
        assert_eq!(mw.shared_hills().len(), 3);
        // One more deposit on walker 0 only.
        mw.walkers[0].maybe_deposit(&[5.0]);
        mw.synchronise();
        for walker in &mw.walkers {
            assert_eq!(walker.hills.len(), 4);
        }
        mw.synchronise();
        mw.synchronise();
        for walker in &mw.walkers {
            assert_eq!(walker.hills.len(), 4);
        }
        assert_eq!(mw.total_hills(), 12);
        assert_eq!(mw.shared_hills().len(), 4);
    }
    #[test]
    fn test_multiple_walkers_bias_counts_every_hill_exactly_once() {
        let mut mw = MultipleWalkersMetadynamics::new(2, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        mw.walkers[0].maybe_deposit(&[0.0]);
        mw.walkers[1].maybe_deposit(&[1.0]);
        for _ in 0..3 {
            mw.synchronise();
        }
        mw.walkers[1].maybe_deposit(&[2.0]);
        for _ in 0..3 {
            mw.synchronise();
        }
        let expected: f64 = [0.0, 1.0, 2.0]
            .iter()
            .map(|&c| GaussianHill::new(vec![c], 1.0, vec![0.5]).evaluate(&[0.8]))
            .sum();
        for walker in &mw.walkers {
            let bias = walker.bias_potential(&[0.8]);
            assert!(
                (bias - expected).abs() < 1e-12,
                "bias {bias} vs the three distinct hills {expected}"
            );
        }
    }
    #[test]
    fn test_multiple_walkers_simultaneous_deposits_are_all_shared() {
        let mut mw = MultipleWalkersMetadynamics::new(2, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        mw.walkers[0].maybe_deposit(&[0.0]);
        mw.walkers[1].maybe_deposit(&[1.0]);
        mw.synchronise();
        // Both walkers deposit again between two synchronisations.
        mw.walkers[0].maybe_deposit(&[2.0]);
        mw.walkers[1].maybe_deposit(&[3.0]);
        mw.walkers[1].maybe_deposit(&[4.0]);
        mw.synchronise();
        for walker in &mw.walkers {
            assert_eq!(walker.hills.len(), 5);
            let mut centers: Vec<f64> = walker.hills.iter().map(|h| h.center[0]).collect();
            centers.sort_by(|a, b| a.total_cmp(b));
            assert_eq!(centers, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        }
    }
    #[test]
    fn test_multiple_walkers_late_walker_and_empty_sync() {
        let mut mw = MultipleWalkersMetadynamics::new(1, || {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        mw.synchronise();
        assert!(mw.shared_hills().is_empty());
        mw.walkers[0].maybe_deposit(&[0.0]);
        mw.synchronise();
        // A walker added after construction brings its own hills into the pool.
        let mut late =
            MetadynamicsState::new(MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0));
        late.maybe_deposit(&[7.0]);
        mw.walkers.push(late);
        mw.synchronise();
        assert_eq!(mw.n_walkers(), 2);
        for walker in &mw.walkers {
            assert_eq!(walker.hills.len(), 2);
        }
        assert_eq!(mw.shared_hills().len(), 2);
    }
    #[test]
    fn test_gaussian_hill_is_clone_and_debug() {
        let hill = GaussianHill::new(vec![0.5, 1.5], 2.0, vec![0.1, 0.2]);
        let copy = hill.clone();
        assert_eq!(copy.center, hill.center);
        assert_eq!(copy.height, hill.height);
        assert_eq!(copy.widths, hill.widths);
        assert!(format!("{hill:?}").contains("GaussianHill"));
    }
    #[test]
    fn test_reconstruct_fes_1d_length() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        for _ in 0..5 {
            state.maybe_deposit(&cv);
        }
        let fes = reconstruct_fes_1d(&state, -5.0, 5.0, 50, 1.0);
        assert_eq!(fes.len(), 50, "FES should have 50 grid points");
    }
    #[test]
    fn test_reconstruct_fes_1d_minimum_zero() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        state.maybe_deposit(&[0.0]);
        let fes = reconstruct_fes_1d(&state, -5.0, 5.0, 50, 1.0);
        let min_val = fes.iter().map(|&(_, f)| f).fold(f64::INFINITY, f64::min);
        assert!(
            min_val.abs() < 1e-10,
            "Minimum FES should be 0 after shift, got {min_val}"
        );
    }
    #[test]
    fn test_cv_monitor_mean() {
        let mut mon = CvMonitor::new("test_cv");
        for i in 0..5 {
            mon.record(i as f64, i as f64);
        }
        assert!(
            (mon.mean() - 2.0).abs() < 1e-12,
            "mean should be 2.0, got {}",
            mon.mean()
        );
    }
    #[test]
    fn test_cv_monitor_std_dev() {
        let mut mon = CvMonitor::new("cv");
        mon.record(0.0, 0.0);
        mon.record(1.0, 0.0);
        assert!(mon.std_dev().abs() < 1e-12);
    }
    #[test]
    fn test_cv_monitor_range() {
        let mut mon = CvMonitor::new("cv");
        mon.record(0.0, 3.0);
        mon.record(1.0, 7.0);
        mon.record(2.0, 1.0);
        let (lo, hi) = mon.range();
        assert!((lo - 1.0).abs() < 1e-12, "min should be 1.0, got {lo}");
        assert!((hi - 7.0).abs() < 1e-12, "max should be 7.0, got {hi}");
    }
    #[test]
    fn test_cv_monitor_empty() {
        let mon = CvMonitor::new("empty_cv");
        assert!(mon.is_empty());
        assert_eq!(mon.mean(), 0.0);
    }
    #[test]
    fn test_check_convergence_not_converged_few_hills() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let state = MetadynamicsState::new(params);
        let converged = check_convergence(&state, -5.0, 5.0, 20, 4, 0.1);
        assert!(!converged, "Should not be converged with no hills");
    }
    #[test]
    fn test_check_convergence_returns_bool() {
        let params = MetadynamicsParams::new(0.1, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        let cv = [0.0f64];
        for _ in 0..20 {
            state.maybe_deposit(&cv);
        }
        let _converged = check_convergence(&state, -3.0, 3.0, 20, 4, 1.0);
    }
}
/// Estimate adaptive widths from recent CV trajectory fluctuations.
///
/// The width for each CV dimension is set to `scale * std(recent_cvs)`.
///
/// # Arguments
/// * `cv_history` – Slice of recent CV vectors (each of length `n_dim`).
/// * `scale`      – Scale factor (typically 1.0 to 2.0).
/// * `min_width`  – Minimum width to avoid numerical issues.
pub fn adaptive_width_from_history(
    cv_history: &[Vec<f64>],
    scale: f64,
    min_width: f64,
) -> Vec<f64> {
    if cv_history.is_empty() {
        return vec![];
    }
    let n_dim = cv_history[0].len();
    let n = cv_history.len() as f64;
    (0..n_dim)
        .map(|d| {
            let mean = cv_history.iter().map(|cv| cv[d]).sum::<f64>() / n;
            let var = cv_history
                .iter()
                .map(|cv| (cv[d] - mean).powi(2))
                .sum::<f64>()
                / n;
            (scale * var.sqrt()).max(min_width)
        })
        .collect()
}
/// Compute free energy difference between two states using exponential averaging.
///
/// ΔF = -k_B*T * ln(<exp(-β*(U_1 - U_0))>_0)
///
/// # Arguments
/// * `delta_u` – Energy differences U_1 - U_0 sampled from state 0.
/// * `temperature` – Simulation temperature (K).
pub fn free_energy_perturbation(delta_u: &[f64], temperature: f64) -> f64 {
    if delta_u.is_empty() {
        return 0.0;
    }
    let beta = 1.0 / (KB * temperature);
    let n = delta_u.len() as f64;
    let du_min = delta_u.iter().cloned().fold(f64::INFINITY, f64::min);
    let sum: f64 = delta_u
        .iter()
        .map(|&du| (-beta * (du - du_min)).exp())
        .sum();
    -temperature * KB * (sum / n).ln() + du_min
}
/// Bennett acceptance ratio (BAR) single-iteration estimate.
///
/// ΔF = C + k_B*T * ln( <f(β(U_0 - U_1 + C))>_1 / <f(β(U_1 - U_0 - C))>_0 )
/// where f(x) = 1 / (1 + exp(x)) is the Fermi function, β = 1/(k_B*T) and C is
/// the trial shift `c_guess`.
///
/// This is Bennett's identity `exp(-βΔF) = <f(β(ΔU - C))>_0 / <f(-β(ΔU - C))>_1 * exp(-βC)`,
/// which holds exactly (in expectation) for every constant C, so a single
/// evaluation already returns a consistent estimate.  The variance is smallest
/// for `C = ΔF + k_B*T*ln(n_1/n_0)`; for equal sample counts feeding the result
/// back in as the new `c_guess` converges to the self-consistent BAR solution.
/// The two Fermi averages are evaluated in log space, so energy differences of
/// thousands of `k_B*T` do not overflow.
///
/// # Arguments
/// * `du_fwd` – U_1 - U_0 sampled from state 0.
/// * `du_rev` – U_0 - U_1 sampled from state 1.
/// * `temperature` – Simulation temperature (K).
/// * `c_guess` – Initial estimate of ΔF (kJ/mol).
///
/// Returns `c_guess` when either sample set is empty or when one of the Fermi
/// averages vanishes identically (no overlap between the two ensembles).
pub fn bar_estimate(du_fwd: &[f64], du_rev: &[f64], temperature: f64, c_guess: f64) -> f64 {
    if du_fwd.is_empty() || du_rev.is_empty() {
        return c_guess;
    }
    let beta = 1.0 / (KB * temperature);
    let ln_mean0 = log_mean_fermi(du_fwd, beta, -beta * c_guess);
    let ln_mean1 = log_mean_fermi(du_rev, beta, beta * c_guess);
    if !(ln_mean0.is_finite() && ln_mean1.is_finite()) {
        return c_guess;
    }
    temperature * KB * (ln_mean1 - ln_mean0) + c_guess
}
#[cfg(test)]
mod tests_extended {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_adaptive_hill_evaluate_at_center() {
        let hill = AdaptiveGaussianHill::new(vec![0.0], 2.0, vec![0.5], 0);
        let v = hill.evaluate(&[0.0]);
        assert!(
            (v - 2.0).abs() < 1e-12,
            "Adaptive hill at center = height, got {v}"
        );
    }
    #[test]
    fn test_adaptive_hill_decays_away_from_center() {
        let hill = AdaptiveGaussianHill::new(vec![0.0], 1.0, vec![0.3], 0);
        let v_center = hill.evaluate(&[0.0]);
        let v_far = hill.evaluate(&[2.0]);
        assert!(
            v_far < v_center,
            "Adaptive hill should decay away from center"
        );
    }
    #[test]
    fn test_adaptive_width_from_history_single_point() {
        let history = vec![vec![1.0f64]];
        let widths = adaptive_width_from_history(&history, 1.0, 0.01);
        assert!(
            (widths[0] - 0.01).abs() < 1e-12,
            "min_width = 0.01, got {}",
            widths[0]
        );
    }
    #[test]
    fn test_adaptive_width_from_history_two_points() {
        let history = vec![vec![0.0f64], vec![2.0f64]];
        let widths = adaptive_width_from_history(&history, 1.0, 0.001);
        assert!(widths[0] > 0.001, "Width should be > min_width");
    }
    #[test]
    fn test_adaptive_width_empty_history() {
        let widths = adaptive_width_from_history(&[], 1.0, 0.01);
        assert!(widths.is_empty(), "Empty history → empty widths");
    }
    #[test]
    fn test_opes_initial_bias_zero() {
        let opes = OpesState::new(10.0, 300.0, 1, vec![0.5]);
        let v = opes.bias_potential(&[0.0]);
        assert_eq!(v, 0.0, "No kernels → zero bias");
    }
    #[test]
    fn test_opes_update_adds_kernel() {
        let mut opes = OpesState::new(10.0, 300.0, 1, vec![0.5]);
        opes.update(&[0.0]);
        assert_eq!(opes.kernels.len(), 1, "Should have 1 kernel after update");
    }
    #[test]
    fn test_opes_stride_deposit_only_at_stride() {
        let mut opes = OpesState::new(10.0, 300.0, 3, vec![0.5]);
        opes.update(&[0.0]);
        opes.update(&[1.0]);
        opes.update(&[2.0]);
        opes.update(&[3.0]);
        assert_eq!(
            opes.kernels.len(),
            2,
            "Should have 2 kernels after 4 updates with stride 3"
        );
    }
    #[test]
    fn test_opes_kernel_density_positive_after_deposit() {
        let mut opes = OpesState::new(10.0, 300.0, 1, vec![0.5]);
        opes.update(&[0.0]);
        let d = opes.kernel_density(&[0.0]);
        assert!(
            d > 0.0,
            "Kernel density should be positive after deposit, got {d}"
        );
    }
    #[test]
    fn test_be_meta_new() {
        let be = BiasExchangeMetadynamics::new(3, 300.0, |_| {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        assert_eq!(be.replicas.len(), 3);
    }
    #[test]
    fn test_be_meta_acceptance_rate_initially_zero() {
        let be = BiasExchangeMetadynamics::new(2, 300.0, |_| {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        assert_eq!(be.acceptance_rate(), 0.0);
    }
    #[test]
    fn test_be_meta_try_exchange_counts_the_attempt() {
        let mut be = BiasExchangeMetadynamics::new(2, 300.0, |_| {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        let accepted = be.try_exchange(0, 1, &[0.0], &[1.0]).unwrap();
        // Without any hills both biases vanish, so Delta = 0 and the move is
        // accepted with probability one.
        assert!(accepted);
        assert_eq!(be.n_attempts, 1);
        assert_eq!(be.n_accepted, 1);
    }
    /// Two replicas at 300 K whose exchange cost for `cv_i = [100]`,
    /// `cv_j = [0]` is exactly `delta`: replica 0 carries one hill of height
    /// `delta * k_B T` at the origin.
    fn bias_exchange_with_cost(delta: f64, seed: u64) -> BiasExchangeMetadynamics {
        let temperature = 300.0;
        let height = delta * KB * temperature;
        let mut be = BiasExchangeMetadynamics::with_seed(
            2,
            temperature,
            |_| MetadynamicsParams::new(height, vec![0.5], 1, temperature, 0.0),
            seed,
        );
        be.replicas[0].maybe_deposit(&[0.0]);
        be
    }
    #[test]
    fn test_be_meta_exchange_delta_matches_the_bias_difference() {
        let temperature = 350.0;
        let mut be = BiasExchangeMetadynamics::new(3, temperature, |k| {
            MetadynamicsParams::new(1.5 + k as f64, vec![0.4], 1, temperature, 0.0)
        });
        be.replicas[0].maybe_deposit(&[0.1]);
        be.replicas[2].maybe_deposit(&[0.9]);
        let (cv_i, cv_j) = ([0.3], [0.7]);
        let beta = 1.0 / (KB * temperature);
        let hill_0 = GaussianHill::new(vec![0.1], 1.5, vec![0.4]);
        let hill_2 = GaussianHill::new(vec![0.9], 3.5, vec![0.4]);
        let expected = beta
            * ((hill_0.evaluate(&cv_j) + hill_2.evaluate(&cv_i))
                - (hill_0.evaluate(&cv_i) + hill_2.evaluate(&cv_j)));
        let delta = be.exchange_delta(0, 2, &cv_i, &cv_j).unwrap();
        assert!((delta - expected).abs() < 1e-12, "{delta} vs {expected}");
        // Swapping the roles of the two replicas gives the same cost.
        let swapped = be.exchange_delta(2, 0, &cv_j, &cv_i).unwrap();
        assert!((delta - swapped).abs() < 1e-12);
    }
    #[test]
    fn test_be_meta_acceptance_rate_follows_boltzmann_factor() {
        // Delta = 1: acceptance probability exp(-1) = 0.3679.  The old rule
        // accepted deterministically whenever Delta < ln 2.
        let mut be = bias_exchange_with_cost(1.0, 7);
        let delta = be.exchange_delta(0, 1, &[100.0], &[0.0]).unwrap();
        assert!((delta - 1.0).abs() < 1e-12, "Delta = {delta}");
        for _ in 0..10_000 {
            be.try_exchange(0, 1, &[100.0], &[0.0]).unwrap();
        }
        assert_eq!(be.n_attempts, 10_000);
        let rate = be.acceptance_rate();
        assert!(
            (rate - (-1.0f64).exp()).abs() < 0.02,
            "acceptance rate {rate}, expected {}",
            (-1.0f64).exp()
        );
        // Delta = 0.5 < ln 2 used to be accepted every single time.
        let mut easy = bias_exchange_with_cost(0.5, 8);
        for _ in 0..10_000 {
            easy.try_exchange(0, 1, &[100.0], &[0.0]).unwrap();
        }
        let easy_rate = easy.acceptance_rate();
        assert!(
            (easy_rate - (-0.5f64).exp()).abs() < 0.02 && easy_rate < 0.95,
            "acceptance rate {easy_rate}, expected {}",
            (-0.5f64).exp()
        );
    }
    #[test]
    fn test_be_meta_favourable_exchange_is_always_accepted() {
        // Reversed roles: the cost is -2, so every attempt is accepted.
        let mut be = bias_exchange_with_cost(2.0, 3);
        let delta = be.exchange_delta(0, 1, &[0.0], &[100.0]).unwrap();
        assert!((delta + 2.0).abs() < 1e-12);
        for _ in 0..1_000 {
            assert!(be.try_exchange(0, 1, &[0.0], &[100.0]).unwrap());
        }
        assert_eq!(be.acceptance_rate(), 1.0);
    }
    #[test]
    fn test_be_meta_unfavourable_exchange_is_almost_never_accepted() {
        let mut be = bias_exchange_with_cost(8.0, 5);
        for _ in 0..10_000 {
            be.try_exchange(0, 1, &[100.0], &[0.0]).unwrap();
        }
        // exp(-8) = 3.4e-4: a handful of acceptances at most.
        assert!(be.acceptance_rate() < 0.003, "{}", be.acceptance_rate());
    }
    #[test]
    fn test_be_meta_streams_are_reproducible_and_seed_dependent() {
        let run = |seed: u64| -> Vec<bool> {
            let mut be = bias_exchange_with_cost(1.0, seed);
            (0..64)
                .map(|_| be.try_exchange(0, 1, &[100.0], &[0.0]).unwrap())
                .collect()
        };
        assert_eq!(run(11), run(11));
        assert_ne!(run(11), run(12));
        // `new` is `with_seed(DEFAULT_SEED)`.
        let params = |_: usize| MetadynamicsParams::new(2.0, vec![0.5], 1, 300.0, 0.0);
        let mut a = BiasExchangeMetadynamics::new(2, 300.0, params);
        let mut b = BiasExchangeMetadynamics::with_seed(
            2,
            300.0,
            params,
            BiasExchangeMetadynamics::DEFAULT_SEED,
        );
        for be in [&mut a, &mut b] {
            be.replicas[0].maybe_deposit(&[0.0]);
        }
        for _ in 0..32 {
            assert_eq!(
                a.try_exchange(0, 1, &[9.0], &[0.0]).unwrap(),
                b.try_exchange(0, 1, &[9.0], &[0.0]).unwrap()
            );
        }
    }
    #[test]
    fn test_be_meta_external_rng_drives_the_decision() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let mut be = bias_exchange_with_cost(1.0, 1);
        let mut rng = StdRng::seed_from_u64(99);
        for _ in 0..10_000 {
            be.try_exchange_with_rng(0, 1, &[100.0], &[0.0], &mut rng)
                .unwrap();
        }
        let rate = be.acceptance_rate();
        assert!((rate - (-1.0f64).exp()).abs() < 0.02, "rate {rate}");
        // Reusing the same seed reproduces the decisions of the first run.
        let mut again = bias_exchange_with_cost(1.0, 2);
        let mut rng_again = StdRng::seed_from_u64(99);
        for _ in 0..10_000 {
            again
                .try_exchange_with_rng(0, 1, &[100.0], &[0.0], &mut rng_again)
                .unwrap();
        }
        assert_eq!(again.n_accepted, be.n_accepted);
    }
    #[test]
    fn test_be_meta_invalid_exchange_requests_are_errors() {
        use crate::MdError;
        let mut be = BiasExchangeMetadynamics::new(2, 300.0, |_| {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        assert!(matches!(
            be.try_exchange(0, 2, &[0.0], &[1.0]),
            Err(MdError::General(_))
        ));
        assert!(matches!(
            be.try_exchange(5, 0, &[0.0], &[1.0]),
            Err(MdError::General(_))
        ));
        assert!(matches!(
            be.try_exchange(1, 1, &[0.0], &[1.0]),
            Err(MdError::General(_))
        ));
        assert_eq!(be.n_attempts, 0, "failed requests are not counted");
        be.temperature = 0.0;
        assert!(be.try_exchange(0, 1, &[0.0], &[1.0]).is_err());
        be.temperature = f64::NAN;
        assert!(be.exchange_delta(0, 1, &[0.0], &[1.0]).is_err());
    }
    #[test]
    fn test_coordination_number_at_zero_distance() {
        let cv = CoordinationNumberCV::new(0, vec![1], 1.0, 6, 12, "cn");
        let positions = vec![[0.0; 3]; 2];
        // sigma(0) = (1 - 0) / (1 - 0) = 1: a coincident neighbour counts fully.
        let v = cv.value(&positions);
        assert!((v - 1.0).abs() < 1e-15, "CN at r = 0 is 1, got {v}");
    }
    #[test]
    fn test_coordination_number_at_large_distance() {
        let cv = CoordinationNumberCV::new(0, vec![1], 1.0, 6, 12, "cn");
        let positions = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let v = cv.value(&positions);
        assert!(v.abs() < 1e-6, "CN at large r should be ≈0, got {v}");
    }
    #[test]
    fn test_coordination_number_at_r0() {
        let cv = CoordinationNumberCV::new(0, vec![1], 1.0, 6, 12, "cn");
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let v = cv.value(&positions);
        assert!(
            (v - 0.5).abs() < 1e-12,
            "CN at r=r0 is the limit n/m = 0.5, got {v}"
        );
    }
    #[test]
    fn test_switching_function_limit_at_r0_is_n_over_m() {
        for &(n, m) in &[(6u32, 12u32), (8, 16), (4, 12), (6, 9), (3, 5)] {
            let cv = CoordinationNumberCV::new(0, vec![], 1.7, n, m, "cn");
            let limit = f64::from(n) / f64::from(m);
            let at_r0 = cv.switching_function(1.7);
            assert!((at_r0 - limit).abs() < 1e-12, "({n},{m}): {at_r0}");
            for &side in &[1.0 - 1e-7, 1.0 + 1e-7, 1.0 - 3e-5, 1.0 + 3e-5] {
                let v = cv.switching_function(1.7 * side);
                assert!(
                    (v - limit).abs() < 1e-4,
                    "({n},{m}) at x = {side}: {v} vs limit {limit}"
                );
            }
        }
    }
    #[test]
    fn test_switching_function_is_smooth_across_the_series_band() {
        let cv = CoordinationNumberCV::new(0, vec![], 1.0, 6, 12, "cn");
        let edge = SWITCHING_SERIES_HALF_WIDTH;
        for &sign in &[-1.0, 1.0] {
            let inside = 1.0 + sign * edge * 0.999_999;
            let outside = 1.0 + sign * edge * 1.000_001;
            let v_in = cv.switching_function(inside);
            let v_out = cv.switching_function(outside);
            assert!((v_in - v_out).abs() < 1e-9, "value jump {v_in} vs {v_out}");
            let d_in = cv.switching_derivative(inside);
            let d_out = cv.switching_derivative(outside);
            assert!(
                (d_in - d_out).abs() < 1e-7,
                "derivative jump {d_in} vs {d_out}"
            );
        }
    }
    #[test]
    fn test_switching_function_matches_closed_form_for_half_exponent() {
        // For m = 2n the ratio simplifies to 1 / (1 + x^n), which has no
        // singularity at x = 1 and so is an independent reference everywhere.
        let cv = CoordinationNumberCV::new(0, vec![], 2.0, 6, 12, "cn");
        for &r in &[0.1, 0.9, 1.9999, 2.0, 2.0001, 2.5, 4.0, 40.0] {
            let x: f64 = r / 2.0;
            let exact = 1.0 / (1.0 + x.powi(6));
            let exact_slope = -6.0 * x.powi(5) / (1.0 + x.powi(6)).powi(2) / 2.0;
            assert!(
                (cv.switching_function(r) - exact).abs() < 1e-11,
                "value at r = {r}"
            );
            assert!(
                (cv.switching_derivative(r) - exact_slope).abs() < 1e-9,
                "slope at r = {r}"
            );
        }
        // At the singular point: sigma = 1/2 and d sigma/dr = -n/(4 r_0) = -0.75.
        assert!((cv.switching_derivative(2.0) + 0.75).abs() < 1e-12);
    }
    #[test]
    fn test_switching_function_derivative_matches_finite_difference() {
        let cv = CoordinationNumberCV::new(0, vec![], 1.3, 8, 14, "cn");
        let h = 1e-6;
        for &r in &[0.2, 0.9, 1.29, 1.3, 1.31, 1.9, 3.0] {
            let fd = (cv.switching_function(r + h) - cv.switching_function(r - h)) / (2.0 * h);
            let analytic = cv.switching_derivative(r);
            assert!(
                (fd - analytic).abs() < 1e-6 * analytic.abs().max(1.0),
                "r = {r}: analytic {analytic} vs finite difference {fd}"
            );
        }
    }
    #[test]
    fn test_switching_function_far_field_and_degenerate_exponents() {
        let cv = CoordinationNumberCV::new(0, vec![], 1.0, 6, 12, "cn");
        assert_eq!(cv.switching_function(1e60), 0.0);
        assert_eq!(cv.switching_derivative(1e60), 0.0);
        assert!((cv.switching_function(0.0) - 1.0).abs() < 1e-15);
        let flat = CoordinationNumberCV::new(0, vec![], 1.0, 4, 4, "cn");
        assert_eq!(flat.switching_function(0.7), 1.0);
        assert_eq!(flat.switching_derivative(0.7), 0.0);
        let zero_n = CoordinationNumberCV::new(0, vec![], 1.0, 0, 12, "cn");
        assert_eq!(zero_n.switching_function(0.7), 0.0);
    }
    #[test]
    fn test_coordination_gradient_matches_finite_difference_of_value() {
        let cv = CoordinationNumberCV::new(0, vec![1, 2, 3], 1.2, 6, 12, "cn");
        let positions = vec![
            [0.1, -0.2, 0.3],
            [1.0, 0.4, -0.1],
            [-0.7, 0.9, 0.5],
            [0.3, -1.1, 0.8],
            [5.0, 5.0, 5.0],
        ];
        let grad = cv.gradient(&positions);
        assert_eq!(grad.len(), positions.len());
        let h = 1e-6;
        for atom in 0..positions.len() {
            for dim in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                plus[atom][dim] += h;
                minus[atom][dim] -= h;
                let fd = (cv.value(&plus) - cv.value(&minus)) / (2.0 * h);
                assert!(
                    (grad[atom][dim] - fd).abs() < 1e-6,
                    "atom {atom} dim {dim}: analytic {} vs finite difference {fd}",
                    grad[atom][dim]
                );
            }
        }
        // The spectator atom does not influence the coordination number.
        assert_eq!(grad[4], [0.0; 3]);
    }
    #[test]
    fn test_coordination_gradient_sums_to_zero() {
        let cv = CoordinationNumberCV::new(1, vec![0, 2, 3], 1.5, 6, 12, "cn");
        let positions = vec![
            [0.4, 0.1, -0.3],
            [0.0, 0.0, 0.0],
            [1.1, 0.7, 0.2],
            [-0.9, 0.8, 1.0],
        ];
        let grad = cv.gradient(&positions);
        for dim in 0..3 {
            let total: f64 = grad.iter().map(|g| g[dim]).sum();
            assert!(total.abs() < 1e-12, "dim {dim}: sum of gradients = {total}");
        }
        // The reference atom carries minus the sum over the neighbours, which
        // the old per-neighbour finite difference scaled by N_nb.
        for (dim, reference) in grad[1].iter().enumerate() {
            let neighbours: f64 = [0usize, 2, 3].iter().map(|&j| grad[j][dim]).sum();
            assert!((reference + neighbours).abs() < 1e-12);
        }
    }
    #[test]
    fn test_coordination_gradient_at_r0_uses_analytic_limit() {
        // sigma'(r_0) = n (n - m) / (2 m r_0) = -1.5 / r_0 for n = 6, m = 12.
        let r_0 = 1.25;
        let cv = CoordinationNumberCV::new(0, vec![1], r_0, 6, 12, "cn");
        let positions = vec![[0.0, 0.0, 0.0], [r_0, 0.0, 0.0]];
        let grad = cv.gradient(&positions);
        assert!((grad[1][0] + 1.5 / r_0).abs() < 1e-12, "{:?}", grad[1]);
        assert!((grad[0][0] - 1.5 / r_0).abs() < 1e-12, "{:?}", grad[0]);
        assert!(grad[1][1].abs() < 1e-15 && grad[1][2].abs() < 1e-15);
    }
    #[test]
    fn test_coordination_gradient_skips_coincident_neighbour() {
        let cv = CoordinationNumberCV::new(0, vec![1, 2], 1.0, 6, 12, "cn");
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.8, 0.0, 0.0]];
        let grad = cv.gradient(&positions);
        assert!(grad.iter().all(|g| g.iter().all(|c| c.is_finite())));
        assert_eq!(grad[1], [0.0; 3]);
        assert!(grad[2][0] < 0.0, "approaching shell edge lowers s");
    }
    #[test]
    fn test_switching_function_at_small_r() {
        let cv = CoordinationNumberCV::new(0, vec![], 1.0, 6, 12, "cn");
        let sigma = cv.switching_function(0.1);
        assert!(
            sigma > 0.9,
            "sigma at r << r_0 should be near 1, got {sigma}"
        );
    }
    #[test]
    fn test_fep_zero_for_empty() {
        let df = free_energy_perturbation(&[], 300.0);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_fep_zero_for_zero_delta_u() {
        let du = vec![0.0; 10];
        let df = free_energy_perturbation(&du, 300.0);
        assert!(df.abs() < 1e-10, "ΔF for ΔU=0 = 0, got {df}");
    }
    #[test]
    fn test_fep_negative_for_negative_delta_u() {
        let du = vec![-10.0; 100];
        let df = free_energy_perturbation(&du, 300.0);
        assert!(df < 0.0, "ΔF for negative ΔU should be negative, got {df}");
    }
    #[test]
    fn test_bar_empty_returns_c_guess() {
        let df = bar_estimate(&[], &[1.0, 2.0], 300.0, 5.0);
        assert_eq!(df, 5.0, "Empty forward → return c_guess");
    }
    #[test]
    fn test_bar_symmetric_gives_zero() {
        let du = vec![0.0f64; 50];
        let df = bar_estimate(&du, &du, 300.0, 0.0);
        assert!(
            df.abs() < 1.0,
            "Symmetric BAR should give small ΔF, got {df}"
        );
    }
    #[test]
    fn test_bar_estimate_is_exact_for_any_trial_shift() {
        // U_1 = U_0 + 2: Bennett's identity returns 2 for every C.
        let du_fwd = vec![2.0f64; 50];
        let du_rev = vec![-2.0f64; 50];
        for &c_guess in &[-5.0, -3.0, 0.0, 2.0, 5.0, 7.0] {
            let df = bar_estimate(&du_fwd, &du_rev, 300.0, c_guess);
            assert!(
                (df - 2.0).abs() < 1e-9,
                "c_guess = {c_guess}: expected 2, got {df}"
            );
        }
    }
    #[test]
    fn test_bar_estimate_is_exact_for_unbalanced_sample_counts() {
        let du_fwd = vec![2.0f64; 50];
        let du_rev = vec![-2.0f64; 10];
        for &c_guess in &[0.0, 1.0, 4.0] {
            let df = bar_estimate(&du_fwd, &du_rev, 300.0, c_guess);
            assert!(
                (df - 2.0).abs() < 1e-9,
                "c_guess = {c_guess}: expected 2, got {df}"
            );
        }
    }
    #[test]
    fn test_bar_estimate_survives_thousands_of_kt() {
        // 5000 kJ/mol is ~2000 kT at 300 K: a naive Fermi function overflows.
        let du_fwd = vec![5000.0f64; 10];
        let du_rev = vec![-5000.0f64; 10];
        let df = bar_estimate(&du_fwd, &du_rev, 300.0, 0.0);
        assert!((df - 5000.0).abs() < 1e-6, "expected 5000, got {df}");
    }
    #[test]
    fn test_bar_estimate_gaussian_work_recovers_crooks_free_energy() {
        use crate::free_energy::test_support::standard_normal;
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let temperature = 300.0;
        let kt = KB * temperature;
        let (mu, sigma) = (10.0, 5.0);
        let exact = mu - sigma * sigma / (2.0 * kt);
        let mut rng = StdRng::seed_from_u64(42);
        let n = 20_000;
        let du_fwd: Vec<f64> = (0..n)
            .map(|_| mu + sigma * standard_normal(&mut rng))
            .collect();
        let du_rev: Vec<f64> = (0..n)
            .map(|_| -(mu - sigma * sigma / kt + sigma * standard_normal(&mut rng)))
            .collect();
        let df = bar_estimate(&du_fwd, &du_rev, temperature, exact);
        assert!(
            (df - exact).abs() < 0.15,
            "BAR = {df} kJ/mol, Crooks value = {exact} kJ/mol"
        );
    }
    #[test]
    fn test_local_elevation_deposit_increases_potential() {
        let mut le = LocalElevation::new(100.0, 0.01);
        let v0 = le.flooding_potential(&[0.0]);
        le.deposit(&[0.0], 1.0, vec![0.5]);
        let v1 = le.flooding_potential(&[0.0]);
        assert!(
            v1 > v0,
            "Depositing repulsor should increase flooding potential"
        );
    }
    #[test]
    fn test_local_elevation_height_clamped() {
        let mut le = LocalElevation::new(1.0, 0.01);
        le.deposit(&[0.0], 100.0, vec![0.5]);
        assert!(
            le.repulsors[0].height <= 1.0,
            "Height should be clamped to max_height"
        );
    }
    #[test]
    fn test_local_elevation_decay_reduces_height() {
        let mut le = LocalElevation::new(10.0, 0.1);
        le.deposit(&[0.0], 5.0, vec![0.5]);
        let h0 = le.repulsors[0].height;
        le.decay();
        let h1 = le.repulsors[0].height;
        assert!(h1 < h0, "Decay should reduce repulsor height: {h0} → {h1}");
    }
    #[test]
    fn test_local_elevation_decay_removes_negligible() {
        let mut le = LocalElevation::new(10.0, 0.9999999);
        le.deposit(&[0.0], 1e-8, vec![0.5]);
        le.decay();
        let _ = le.repulsors.len();
    }
    #[test]
    fn test_bias_force_nonzero_away_from_center() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        state.maybe_deposit(&[0.0]);
        let force = state.bias_force(&[1.0]);
        assert!(!force.is_empty());
        assert!(force[0].is_finite(), "Bias force should be finite");
    }
    #[test]
    fn test_gaussian_hill_height_positive() {
        let hill = GaussianHill::new(vec![0.0], 2.5, vec![0.3]);
        assert!(hill.evaluate(&[0.0]) > 0.0);
        assert!((hill.evaluate(&[0.0]) - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_steered_md_work_sign_on_extension() {
        let smd = SteeredMD::new(0, 1.0, 100.0, 0.1);
        let cv_traj = vec![1.0, 1.1, 1.2, 1.3, 1.4];
        let work = smd.work_done(&cv_traj, 1.0);
        assert!(work.is_finite(), "Work should be finite");
    }
    #[test]
    fn test_steered_md_no_work_at_equilibrium() {
        let smd = SteeredMD::new(0, 1.0, 100.0, 0.0);
        let cv_traj = vec![1.0; 5];
        let work = smd.work_done(&cv_traj, 1.0);
        assert!(work.abs() < 1e-12, "No displacement → no work, got {work}");
    }
    #[test]
    fn test_funnel_potential_zero_inside_funnel() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let state = MetadynamicsState::new(params);
        let funnel = FunnelMetadynamics::new(state, 1.0, 0.5);
        let u = funnel.funnel_potential(0.5, 1.0);
        assert_eq!(u, 0.0, "Inside funnel: potential = 0");
    }
    #[test]
    fn test_funnel_potential_positive_outside_funnel() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let state = MetadynamicsState::new(params);
        let funnel = FunnelMetadynamics::new(state, 0.5, 0.1);
        let u = funnel.funnel_potential(2.0, 0.5);
        assert!(
            u > 0.0,
            "Outside funnel: potential should be positive, got {u}"
        );
    }
    #[test]
    fn test_rg_cv_value_positive() {
        let cv = RadiusOfGyrationCV {
            atom_indices: vec![0, 1, 2, 3],
            name: "Rg".into(),
        };
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let rg = cv.value(&positions);
        assert!(rg > 0.0, "Rg should be positive, got {rg}");
    }
    #[test]
    fn test_angle_cv_known_angle() {
        let cv = AngleCV::new(0, 1, 2, "theta");
        let positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let theta = cv.value(&positions);
        assert!(
            (theta - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "angle = {theta}"
        );
    }
}
/// OPES compression: merges nearby kernels to control the number of kernels.
///
/// When two kernels are closer than `compression_threshold * bandwidth`,
/// they are merged into a single kernel with combined height.
pub fn compress_opes_kernels(kernels: &mut Vec<GaussianHill>, compression_threshold: f64) {
    if kernels.len() < 2 {
        return;
    }
    let mut merged = true;
    while merged {
        merged = false;
        let n = kernels.len();
        let mut to_remove = vec![false; n];
        for i in 0..n {
            if to_remove[i] {
                continue;
            }
            for j in (i + 1)..n {
                if to_remove[j] {
                    continue;
                }
                let dist_sq: f64 = kernels[i]
                    .center
                    .iter()
                    .zip(kernels[j].center.iter())
                    .zip(kernels[i].widths.iter())
                    .map(|((&ci, &cj), &w)| {
                        let w2 = (w * w).max(1e-30);
                        (ci - cj).powi(2) / w2
                    })
                    .sum();
                if dist_sq < compression_threshold * compression_threshold {
                    let h_sum = kernels[i].height + kernels[j].height;
                    if h_sum < 1e-30 {
                        continue;
                    }
                    let dim = kernels[i].center.len();
                    let mut new_center = vec![0.0; dim];
                    for (d, nc) in new_center.iter_mut().enumerate() {
                        *nc = (kernels[i].height * kernels[i].center[d]
                            + kernels[j].height * kernels[j].center[d])
                            / h_sum;
                    }
                    let new_widths = kernels[i].widths.clone();
                    kernels[i] = GaussianHill::new(new_center, h_sum, new_widths);
                    to_remove[j] = true;
                    merged = true;
                }
            }
        }
        let mut idx = n;
        while idx > 0 {
            idx -= 1;
            if to_remove[idx] {
                kernels.remove(idx);
            }
        }
    }
}
/// Free energy `F(s) = -k_B*T ln(p_hat(s) / p_hat_max)` of the unbiased system
/// at each grid point, estimated from OPES kernels (so that
/// `ln p(s) = -beta * F(s)` up to a constant).
///
/// The kernel heights must be the OPES reweighting factors, as deposited by
/// [`OpesState::update`](crate::metadynamics::OpesState::update); only then
/// does the kernel sum estimate the unbiased distribution instead of the
/// biased one.
pub fn opes_target_free_energy(
    kernels: &[GaussianHill],
    grid: &[f64],
    temperature: f64,
) -> Vec<f64> {
    if kernels.is_empty() || grid.is_empty() {
        return vec![0.0; grid.len()];
    }
    let p: Vec<f64> = grid
        .iter()
        .map(|&s| {
            kernels
                .iter()
                .map(|k| k.evaluate(&[s]))
                .sum::<f64>()
                .max(1e-300)
        })
        .collect();
    let p_max = p.iter().cloned().fold(0.0_f64, f64::max).max(1e-300);
    p.iter()
        .map(|&pi| {
            let log_ratio = (pi / p_max).ln();
            -KB * temperature * log_ratio
        })
        .collect()
}
/// Estimate free energy perturbation ΔF from metadynamics bias potential.
///
/// Uses the relation: ΔF(A→B) = -kT * ln(∫_B exp(V_bias/kT) / ∫_A exp(V_bias/kT))
///
/// Numerically evaluates the integrals on a 1D grid.
pub fn fep_from_metadynamics(
    state: &MetadynamicsState,
    s_min_a: f64,
    s_max_a: f64,
    s_min_b: f64,
    s_max_b: f64,
    n_grid: usize,
    temperature: f64,
) -> f64 {
    if n_grid == 0 {
        return 0.0;
    }
    let beta = 1.0 / (KB * temperature);
    let integrate = |s_lo: f64, s_hi: f64| -> f64 {
        if s_lo >= s_hi || n_grid < 2 {
            return 0.0;
        }
        let ds = (s_hi - s_lo) / (n_grid - 1) as f64;
        let sum: f64 = (0..n_grid)
            .map(|i| {
                let s = s_lo + i as f64 * ds;
                let v = state.bias_potential(&[s]);
                (beta * v).exp()
            })
            .sum();
        sum * ds
    };
    let z_a = integrate(s_min_a, s_max_a);
    let z_b = integrate(s_min_b, s_max_b);
    if z_a < 1e-300 || z_b < 1e-300 {
        return 0.0;
    }
    -KB * temperature * (z_b / z_a).ln()
}
/// Build a histogram of CV values from a trajectory.
///
/// Returns (bin_centers, counts).
pub fn build_cv_histogram(
    cv_trajectory: &[f64],
    s_min: f64,
    s_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    if n_bins == 0 || s_min >= s_max {
        return (Vec::new(), Vec::new());
    }
    let bin_width = (s_max - s_min) / n_bins as f64;
    let mut counts = vec![0usize; n_bins];
    for &s in cv_trajectory {
        if s < s_min || s >= s_max {
            continue;
        }
        let idx = ((s - s_min) / bin_width) as usize;
        let idx = idx.min(n_bins - 1);
        counts[idx] += 1;
    }
    let centers: Vec<f64> = (0..n_bins)
        .map(|i| s_min + (i as f64 + 0.5) * bin_width)
        .collect();
    (centers, counts)
}
/// Free energy from histogram: F(s) = -kT * ln(P(s)).
///
/// Returns free energy at each bin center, shifted so minimum is 0.
pub fn free_energy_from_histogram(counts: &[usize], temperature: f64) -> Vec<f64> {
    let n = counts.iter().sum::<usize>();
    if n == 0 {
        return vec![0.0; counts.len()];
    }
    let nf = n as f64;
    let mut fe: Vec<f64> = counts
        .iter()
        .map(|&c| {
            if c == 0 {
                f64::INFINITY
            } else {
                -KB * temperature * (c as f64 / nf).ln()
            }
        })
        .collect();
    let min_fe = fe
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    for v in fe.iter_mut() {
        if v.is_finite() {
            *v -= min_fe;
        }
    }
    fe
}
#[cfg(test)]
mod tests_ptmetad_extended {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_ptmetad_n_replicas() {
        let temps = vec![300.0, 400.0, 500.0];
        let heights = vec![1.0; 3];
        let widths = vec![vec![0.5]; 3];
        let pt = PtMetaD::new(temps, heights, widths, 1, 1.0);
        assert_eq!(pt.n_replicas(), 3);
    }
    #[test]
    fn test_ptmetad_exchange_attempt_counted() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let temps = vec![300.0, 400.0];
        let heights = vec![1.0; 2];
        let widths = vec![vec![0.5]; 2];
        let mut pt = PtMetaD::new(temps, heights, widths, 1, 1.0);
        let mut rng = StdRng::seed_from_u64(1);
        let accepted = pt
            .try_exchange_adjacent(0, &[0.0], &[1.0], -50.0, -50.0, &mut rng)
            .unwrap();
        assert_eq!(pt.n_attempts, 1);
        assert_eq!(pt.n_accepted, usize::from(accepted));
    }
    /// Two replicas (300 K and 330 K, no hills) and the potential-energy gap
    /// `U_j - U_i` for which the exchange cost is exactly `delta`.
    fn ptmetad_with_cost(delta: f64) -> (PtMetaD, f64) {
        let temps = vec![300.0, 330.0];
        let pt = PtMetaD::new(temps, vec![1.0; 2], vec![vec![0.5]; 2], 1, 0.0);
        let beta_gap = 1.0 / (KB * 300.0) - 1.0 / (KB * 330.0);
        (pt, delta / beta_gap)
    }
    #[test]
    fn test_ptmetad_exchange_delta_includes_both_biases_and_energies() {
        let temps = vec![300.0, 360.0];
        let mut pt = PtMetaD::new(temps, vec![2.0, 3.0], vec![vec![0.5], vec![0.6]], 1, 0.0);
        pt.replicas[0].maybe_deposit(&[0.0]);
        pt.replicas[1].maybe_deposit(&[1.0]);
        let (cv_i, cv_j) = ([0.2], [0.9]);
        let (u_i, u_j) = (-120.0, -100.0);
        let beta_i = 1.0 / (KB * 300.0);
        let beta_j = 1.0 / (KB * 360.0);
        let hill_0 = GaussianHill::new(vec![0.0], 2.0, vec![0.5]);
        let hill_1 = GaussianHill::new(vec![1.0], 3.0, vec![0.6]);
        let expected = (beta_i - beta_j) * (u_j - u_i)
            + beta_i * (hill_0.evaluate(&cv_j) - hill_0.evaluate(&cv_i))
            + beta_j * (hill_1.evaluate(&cv_i) - hill_1.evaluate(&cv_j));
        let delta = pt.exchange_delta(0, &cv_i, &cv_j, u_i, u_j).unwrap();
        assert!((delta - expected).abs() < 1e-10, "{delta} vs {expected}");
        // The partner bias is not optional: dropping replica 1's hill changes Delta.
        let mut without_partner_bias = PtMetaD::new(
            vec![300.0, 360.0],
            vec![2.0, 3.0],
            vec![vec![0.5], vec![0.6]],
            1,
            0.0,
        );
        without_partner_bias.replicas[0].maybe_deposit(&[0.0]);
        let partial = without_partner_bias
            .exchange_delta(0, &cv_i, &cv_j, u_i, u_j)
            .unwrap();
        assert!((partial - delta).abs() > 1e-3);
    }
    #[test]
    fn test_ptmetad_acceptance_rate_follows_boltzmann_factor() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        for &(delta, tolerance) in &[(1.0, 0.01), (0.5, 0.01), (3.0, 0.005)] {
            let (mut pt, gap) = ptmetad_with_cost(delta);
            let computed = pt.exchange_delta(0, &[0.0], &[0.0], 0.0, gap).unwrap();
            assert!((computed - delta).abs() < 1e-9, "Delta = {computed}");
            let mut rng = StdRng::seed_from_u64(2025);
            for _ in 0..100_000 {
                pt.try_exchange_adjacent(0, &[0.0], &[0.0], 0.0, gap, &mut rng)
                    .unwrap();
            }
            let expected = (-delta).exp();
            let rate = pt.acceptance_rate();
            assert!(
                (rate - expected).abs() < tolerance,
                "Delta = {delta}: acceptance {rate}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_ptmetad_downhill_exchange_is_always_accepted_and_uphill_is_not() {
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let (mut pt, gap) = ptmetad_with_cost(2.0);
        let mut rng = StdRng::seed_from_u64(4);
        // Reversed energy gap: Delta = -2, every attempt succeeds.
        for _ in 0..1_000 {
            assert!(
                pt.try_exchange_adjacent(0, &[0.0], &[0.0], gap, 0.0, &mut rng)
                    .unwrap()
            );
        }
        // Delta = +2 used to be accepted unconditionally by the broken test
        // `delta >= 0 || exp(-delta) > 0.5`; now it is the exception.
        let mut uphill = 0;
        for _ in 0..1_000 {
            if pt
                .try_exchange_adjacent(0, &[0.0], &[0.0], 0.0, gap, &mut rng)
                .unwrap()
            {
                uphill += 1;
            }
        }
        assert!(uphill < 250, "{uphill} of 1000 uphill moves accepted");
    }
    #[test]
    fn test_ptmetad_invalid_replica_is_an_error_not_a_panic() {
        use crate::MdError;
        use rand::SeedableRng;
        use rand::rngs::StdRng;
        let (mut pt, _) = ptmetad_with_cost(1.0);
        let mut rng = StdRng::seed_from_u64(0);
        let err = pt.try_exchange_adjacent(1, &[0.0], &[0.0], 0.0, 0.0, &mut rng);
        assert!(matches!(err, Err(MdError::General(_))));
        let err = pt.try_exchange_adjacent(usize::MAX, &[0.0], &[0.0], 0.0, 0.0, &mut rng);
        assert!(matches!(err, Err(MdError::General(_))));
        assert_eq!(pt.n_attempts, 0);
        assert!(
            pt.exchange_delta(0, &[0.0], &[0.0], f64::NAN, 0.0).is_err(),
            "NaN energies must be rejected"
        );
        pt.temperatures[1] = -5.0;
        assert!(pt.exchange_delta(0, &[0.0], &[0.0], 0.0, 0.0).is_err());
    }
    #[test]
    fn test_ptmetad_total_hills_accumulate() {
        let temps = vec![300.0, 400.0];
        let heights = vec![1.0; 2];
        let widths = vec![vec![0.5]; 2];
        let mut pt = PtMetaD::new(temps, heights, widths, 1, 0.0);
        pt.replicas[0].maybe_deposit(&[0.0]);
        pt.replicas[1].maybe_deposit(&[1.0]);
        assert_eq!(pt.total_hills(), 2);
    }
    #[test]
    fn test_ptmetad_acceptance_rate_zero_initially() {
        let temps = vec![300.0, 400.0];
        let heights = vec![1.0; 2];
        let widths = vec![vec![0.5]; 2];
        let pt = PtMetaD::new(temps, heights, widths, 1, 0.0);
        assert_eq!(pt.acceptance_rate(), 0.0);
    }
    #[test]
    fn test_funnel_wall_potential_zero_inside() {
        let fc = FunnelConfig::new(2.0, 0.2, 10.0, 100.0);
        let u = fc.wall_potential(2.0, 1.0);
        assert_eq!(u, 0.0, "Inside funnel: potential = 0");
    }
    #[test]
    fn test_funnel_wall_potential_positive_outside() {
        let fc = FunnelConfig::new(0.5, 0.1, 5.0, 100.0);
        let u = fc.wall_potential(0.0, 3.0);
        assert!(u > 0.0, "Outside funnel: potential > 0, got {u}");
    }
    #[test]
    fn test_funnel_max_radius_increases_with_z() {
        let fc = FunnelConfig::new(1.0, 0.3, 10.0, 50.0);
        let r1 = fc.max_radius(1.0);
        let r2 = fc.max_radius(5.0);
        assert!(r2 > r1, "radius increases with z: {r1} → {r2}");
    }
    #[test]
    fn test_funnel_volume_correction_positive() {
        let fc = FunnelConfig::new(2.0, 0.2, 5.0, 100.0);
        assert!(fc.volume_correction() > 0.0);
    }
    #[test]
    fn test_compress_kernels_reduces_count() {
        let mut kernels = vec![
            GaussianHill::new(vec![0.0], 1.0, vec![0.5]),
            GaussianHill::new(vec![0.01], 1.0, vec![0.5]),
        ];
        let n_before = kernels.len();
        compress_opes_kernels(&mut kernels, 0.5);
        assert!(
            kernels.len() <= n_before,
            "compression should not increase count"
        );
    }
    #[test]
    fn test_compress_kernels_far_apart_unchanged() {
        let mut kernels = vec![
            GaussianHill::new(vec![0.0], 1.0, vec![0.5]),
            GaussianHill::new(vec![10.0], 1.0, vec![0.5]),
        ];
        compress_opes_kernels(&mut kernels, 0.1);
        assert_eq!(kernels.len(), 2, "far-apart kernels should not merge");
    }
    #[test]
    fn test_opes_fes_length() {
        let kernels = vec![GaussianHill::new(vec![0.0], 1.0, vec![0.5])];
        let grid: Vec<f64> = (0..10).map(|i| i as f64 * 0.5 - 2.0).collect();
        let fes = opes_target_free_energy(&kernels, &grid, 300.0);
        assert_eq!(fes.len(), grid.len());
    }
    #[test]
    fn test_opes_fes_nonnegative() {
        let kernels = vec![GaussianHill::new(vec![0.0], 2.0, vec![1.0])];
        let grid: Vec<f64> = (0..20).map(|i| i as f64 * 0.5 - 5.0).collect();
        let fes = opes_target_free_energy(&kernels, &grid, 300.0);
        let min_val = fes.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(min_val >= 0.0, "FES minimum should be ≥ 0, got {min_val}");
    }
    #[test]
    fn test_fep_meta_returns_finite() {
        let params = MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0);
        let mut state = MetadynamicsState::new(params);
        for _ in 0..5 {
            state.maybe_deposit(&[0.0]);
        }
        let df = fep_from_metadynamics(&state, -2.0, -0.5, 0.5, 2.0, 20, 300.0);
        assert!(df.is_finite(), "FEP from metadynamics = {df}");
    }
    #[test]
    fn test_fep_meta_zero_for_same_region() {
        let params = MetadynamicsParams::new(0.0, vec![0.5], 1, 300.0, 0.0);
        let state = MetadynamicsState::new(params);
        let df = fep_from_metadynamics(&state, -1.0, 0.0, 0.0, 1.0, 10, 300.0);
        assert!(df.abs() < 1e-10, "symmetric flat FES → ΔF ≈ 0, got {df}");
    }
    #[test]
    fn test_e2e_cv_value() {
        let cv = EndToEndDistanceCV::new(0, 3, "e2e");
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let d = cv.value(&positions);
        assert!((d - 3.0).abs() < 1e-12, "e2e distance = 3, got {d}");
    }
    #[test]
    fn test_e2e_cv_gradient_nonzero() {
        let cv = EndToEndDistanceCV::new(0, 1, "e2e");
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        let grad = cv.gradient(&positions);
        let gnorm = grad
            .iter()
            .map(|g| g.iter().map(|&v| v * v).sum::<f64>())
            .sum::<f64>()
            .sqrt();
        assert!(gnorm > 0.0, "gradient should be nonzero");
    }
    #[test]
    fn test_histogram_bin_count() {
        let traj: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let (centers, counts) = build_cv_histogram(&traj, 0.0, 10.0, 10);
        assert_eq!(centers.len(), 10);
        assert_eq!(counts.len(), 10);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 100);
    }
    #[test]
    fn test_histogram_fes_minimum_zero() {
        let counts = vec![10usize, 20, 30, 20, 10];
        let fes = free_energy_from_histogram(&counts, 300.0);
        let min = fes.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            min.abs() < 1e-10,
            "min FES should be 0 after shift, got {min}"
        );
    }
    #[test]
    fn test_histogram_fes_length() {
        let counts = vec![5usize; 8];
        let fes = free_energy_from_histogram(&counts, 300.0);
        assert_eq!(fes.len(), 8);
    }
    #[test]
    fn test_histogram_zero_count_gives_inf() {
        let counts = vec![0usize, 5, 10];
        let fes = free_energy_from_histogram(&counts, 300.0);
        assert!(fes[0].is_infinite(), "zero count → infinite FES");
    }
    #[test]
    fn test_histogram_uniform_gives_flat_fes() {
        let counts = vec![10usize; 5];
        let fes = free_energy_from_histogram(&counts, 300.0);
        for &f in &fes {
            assert!(f.abs() < 1e-10, "uniform histogram → flat FES = 0, got {f}");
        }
    }
}
