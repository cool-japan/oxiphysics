//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{GaussianHill, MetadynamicsState};

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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
/// ΔF ≈ k_B*T * ln(<f(U_0 - U_1 + C)>_1 / <f(U_1 - U_0 - C)>_0) + C
/// where f(x) = 1 / (1 + exp(x*beta)) is the Fermi function and C is an initial guess.
///
/// # Arguments
/// * `du_fwd` – U_1 - U_0 sampled from state 0.
/// * `du_rev` – U_0 - U_1 sampled from state 1.
/// * `temperature` – Simulation temperature (K).
/// * `c_guess` – Initial estimate of ΔF (kJ/mol).
#[allow(dead_code)]
pub fn bar_estimate(du_fwd: &[f64], du_rev: &[f64], temperature: f64, c_guess: f64) -> f64 {
    if du_fwd.is_empty() || du_rev.is_empty() {
        return c_guess;
    }
    let beta = 1.0 / (KB * temperature);
    let fermi = |x: f64| 1.0 / (1.0 + (beta * x).exp());
    let n0 = du_fwd.len() as f64;
    let n1 = du_rev.len() as f64;
    let mean0: f64 = du_fwd.iter().map(|&du| fermi(du - c_guess)).sum::<f64>() / n0;
    let mean1: f64 = du_rev.iter().map(|&du| fermi(-du + c_guess)).sum::<f64>() / n1;
    if mean1.abs() < 1e-30 {
        return c_guess;
    }
    temperature * KB * (mean0 / mean1).ln() + c_guess
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
    fn test_be_meta_try_exchange_returns_bool() {
        let mut be = BiasExchangeMetadynamics::new(2, 300.0, |_| {
            MetadynamicsParams::new(1.0, vec![0.5], 1, 300.0, 0.0)
        });
        let _accepted = be.try_exchange(0, 1, &[0.0], &[1.0]);
        assert_eq!(be.n_attempts, 1);
    }
    #[test]
    fn test_coordination_number_at_zero_distance() {
        let cv = CoordinationNumberCV::new(0, vec![1], 1.0, 6, 12, "cn");
        let positions = vec![[0.0; 3]; 2];
        let _v = cv.value(&positions);
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
        assert!(v.is_finite(), "CN at r=r0 should be finite, got {v}");
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
#[allow(dead_code)]
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
                    for d in 0..dim {
                        new_center[d] = (kernels[i].height * kernels[i].center[d]
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
/// OPES target distribution: log p_tgt(s) = -beta * F(s).
///
/// Returns estimated free energy at each grid point using kernels.
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
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
        let temps = vec![300.0, 400.0];
        let heights = vec![1.0; 2];
        let widths = vec![vec![0.5]; 2];
        let mut pt = PtMetaD::new(temps, heights, widths, 1, 1.0);
        pt.try_exchange_adjacent(0, &[0.0], &[1.0]);
        assert_eq!(pt.n_attempts, 1);
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
