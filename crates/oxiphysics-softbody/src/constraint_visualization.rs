#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constraint visualisation utilities for soft-body simulations.
//!
//! Provides debug statistics (violation histograms, RMS, max violation, colour
//! mapping) and per-particle energy density maps that help diagnose solver
//! convergence and constraint quality during development.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// ConstraintDebugInfo
// ---------------------------------------------------------------------------

/// Diagnostic information for a single constraint at one frame.
#[derive(Debug, Clone)]
pub struct ConstraintDebugInfo {
    /// Index of the first particle.
    pub particle_a: usize,
    /// Index of the second particle.
    pub particle_b: usize,
    /// Rest length of the constraint (m).
    pub rest_length: f64,
    /// Current distance between the two particles (m).
    pub current_length: f64,
    /// Signed violation: current_length − rest_length (m).
    pub violation: f64,
    /// Constraint stiffness (N/m or compliance inverse).
    pub stiffness: f64,
}

impl ConstraintDebugInfo {
    /// Construct a [`ConstraintDebugInfo`] from raw values.
    pub fn new(
        particle_a: usize,
        particle_b: usize,
        rest_length: f64,
        current_length: f64,
        stiffness: f64,
    ) -> Self {
        Self {
            particle_a,
            particle_b,
            rest_length,
            current_length,
            violation: current_length - rest_length,
            stiffness,
        }
    }
}

// ---------------------------------------------------------------------------
// ConstraintVisualization
// ---------------------------------------------------------------------------

/// Collection of constraint debug data for a single simulation frame.
#[derive(Debug, Clone, Default)]
pub struct ConstraintVisualization {
    /// Per-constraint debug records.
    pub infos: Vec<ConstraintDebugInfo>,
    /// Frame index (0-based).
    pub frame: usize,
}

impl ConstraintVisualization {
    /// Create an empty [`ConstraintVisualization`] for the given frame.
    pub fn new(frame: usize) -> Self {
        Self {
            infos: Vec::new(),
            frame,
        }
    }

    /// Push a single debug record into this frame.
    pub fn push(&mut self, info: ConstraintDebugInfo) {
        self.infos.push(info);
    }
}

// ---------------------------------------------------------------------------
// compute_constraint_violations
// ---------------------------------------------------------------------------

/// Compute signed violations (current_length − rest_length) for each constraint.
///
/// # Arguments
/// * `positions`    – particle positions, each `[x, y, z]` (m).
/// * `constraints`  – list of `(particle_a, particle_b, rest_length)` tuples.
///
/// Returns a `Vec`f64` of the same length as `constraints`.
pub fn compute_constraint_violations(
    positions: &[[f64; 3]],
    constraints: &[(usize, usize, f64)],
) -> Vec<f64> {
    constraints
        .iter()
        .map(|&(a, b, rest)| {
            let pa = positions[a];
            let pb = positions[b];
            let dx = pa[0] - pb[0];
            let dy = pa[1] - pb[1];
            let dz = pa[2] - pb[2];
            let cur = (dx * dx + dy * dy + dz * dz).sqrt();
            cur - rest
        })
        .collect()
}

// ---------------------------------------------------------------------------
// violation_histogram
// ---------------------------------------------------------------------------

/// Build a histogram of absolute violation magnitudes.
///
/// The histogram spans `\[0, max_abs_violation\]`, divided into `n_bins` equal bins.
/// Returns a `Vec`usize` of bin counts with length `n_bins`.
/// Returns an empty vector when `violations` is empty or `n_bins` is 0.
pub fn violation_histogram(violations: &[f64], n_bins: usize) -> Vec<usize> {
    if violations.is_empty() || n_bins == 0 {
        return vec![0; n_bins];
    }
    let abs_max = violations.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if abs_max < 1e-30 {
        let mut bins = vec![0usize; n_bins];
        bins[0] = violations.len();
        return bins;
    }
    let mut bins = vec![0usize; n_bins];
    for &v in violations {
        let frac = v.abs() / abs_max;
        let idx = ((frac * n_bins as f64) as usize).min(n_bins - 1);
        bins[idx] += 1;
    }
    bins
}

// ---------------------------------------------------------------------------
// max_violation
// ---------------------------------------------------------------------------

/// Return the maximum absolute violation (m).
///
/// Returns 0.0 for an empty slice.
pub fn max_violation(violations: &[f64]) -> f64 {
    violations.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
}

// ---------------------------------------------------------------------------
// rms_violation
// ---------------------------------------------------------------------------

/// Root-mean-square violation magnitude (m).
///
/// Returns 0.0 for an empty slice.
pub fn rms_violation(violations: &[f64]) -> f64 {
    if violations.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = violations.iter().map(|v| v * v).sum();
    (sum_sq / violations.len() as f64).sqrt()
}

// ---------------------------------------------------------------------------
// color_by_violation
// ---------------------------------------------------------------------------

/// Map a violation magnitude to a green-to-red RGB colour.
///
/// Returns `[r, g, b]` where each channel is in `[0.0, 1.0]`.
///
/// * `violation = 0`        → green `[0, 1, 0]`
/// * `violation = max_viol` → red   `[1, 0, 0]`
/// * In between            → linear interpolation through yellow.
///
/// # Arguments
/// * `violation`  – absolute violation magnitude (m).
/// * `max_viol`   – reference maximum violation used for normalisation.
pub fn color_by_violation(violation: f64, max_viol: f64) -> [f64; 3] {
    if max_viol <= 0.0 {
        return [0.0, 1.0, 0.0];
    }
    let t = (violation.abs() / max_viol).clamp(0.0, 1.0);
    // green [0,1,0] → yellow [1,1,0] → red [1,0,0]
    let r = (2.0 * t).min(1.0);
    let g = (2.0 * (1.0 - t)).min(1.0);
    [r, g, 0.0]
}

// ---------------------------------------------------------------------------
// constraint_force_magnitude
// ---------------------------------------------------------------------------

/// Estimate the magnitude of the constraint correction force.
///
/// `|F| = stiffness × |violation|`
///
/// # Arguments
/// * `violation`  – signed or unsigned violation (m).
/// * `stiffness`  – constraint stiffness (N/m).
///
/// Returns force magnitude (N).
pub fn constraint_force_magnitude(violation: f64, stiffness: f64) -> f64 {
    stiffness * violation.abs()
}

// ---------------------------------------------------------------------------
// energy_density_map
// ---------------------------------------------------------------------------

/// Per-particle potential energy density from spring-like constraints.
///
/// Each constraint contributes `0.5 · k · violation²` split equally between its
/// two endpoint particles.
///
/// # Arguments
/// * `positions`   – particle positions `[x, y, z]` (m).
/// * `constraints` – list of `(particle_a, particle_b, rest_length)`.
/// * `stiffness`   – uniform stiffness k (N/m) for all constraints.
///
/// Returns a `Vec`f64` of per-particle energy (J) with length `positions.len()`.
pub fn energy_density_map(
    positions: &[[f64; 3]],
    constraints: &[(usize, usize, f64)],
    stiffness: f64,
) -> Vec<f64> {
    let mut energy = vec![0.0_f64; positions.len()];
    for &(a, b, rest) in constraints {
        let pa = positions[a];
        let pb = positions[b];
        let dx = pa[0] - pb[0];
        let dy = pa[1] - pb[1];
        let dz = pa[2] - pb[2];
        let cur = (dx * dx + dy * dy + dz * dz).sqrt();
        let viol = cur - rest;
        let e_half = 0.25 * stiffness * viol * viol; // half split to each particle
        energy[a] += e_half;
        energy[b] += e_half;
    }
    energy
}

// ---------------------------------------------------------------------------
// convergence_metric
// ---------------------------------------------------------------------------

/// Compute per-frame RMS violation from a history of violation vectors.
///
/// # Arguments
/// * `violations_history` – slice of per-frame violation vectors.
///
/// Returns a `Vec`f64` of RMS values, one per frame.
pub fn convergence_metric(violations_history: &[Vec<f64>]) -> Vec<f64> {
    violations_history
        .iter()
        .map(|frame_violations| rms_violation(frame_violations))
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a simple two-particle line
    fn two_particle_setup() -> (Vec<[f64; 3]>, Vec<(usize, usize, f64)>) {
        let positions = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let constraints = vec![(0usize, 1usize, 1.0_f64)];
        (positions, constraints)
    }

    // --- compute_constraint_violations ---

    #[test]
    fn test_violation_exact() {
        let (pos, cons) = two_particle_setup();
        let viols = compute_constraint_violations(&pos, &cons);
        assert_eq!(viols.len(), 1);
        assert!((viols[0] - 0.5).abs() < 1e-12, "violation = 0.5");
    }

    #[test]
    fn test_violation_at_rest() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0)];
        let viols = compute_constraint_violations(&positions, &constraints);
        assert!(viols[0].abs() < 1e-12, "at rest → zero violation");
    }

    #[test]
    fn test_violation_compressed() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0)];
        let viols = compute_constraint_violations(&positions, &constraints);
        assert!(
            (viols[0] - (-0.5)).abs() < 1e-12,
            "compressed → negative violation"
        );
    }

    #[test]
    fn test_violation_empty_constraints() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let constraints: Vec<(usize, usize, f64)> = vec![];
        let viols = compute_constraint_violations(&positions, &constraints);
        assert!(viols.is_empty());
    }

    #[test]
    fn test_violation_multiple_constraints() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0), (1, 2, 1.0)];
        let viols = compute_constraint_violations(&positions, &constraints);
        assert_eq!(viols.len(), 2);
        assert!((viols[0] - 1.0).abs() < 1e-12);
        assert!((viols[1] - 1.0).abs() < 1e-12);
    }

    // --- violation_histogram ---

    #[test]
    fn test_histogram_all_same() {
        let viols = vec![1.0, 1.0, 1.0, 1.0];
        let hist = violation_histogram(&viols, 4);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 4, "all 4 violations should appear in histogram");
    }

    #[test]
    fn test_histogram_length() {
        let viols = vec![0.1, 0.5, 0.9];
        let hist = violation_histogram(&viols, 10);
        assert_eq!(hist.len(), 10);
    }

    #[test]
    fn test_histogram_empty() {
        let hist = violation_histogram(&[], 5);
        assert_eq!(hist.len(), 5);
        assert!(hist.iter().all(|&c| c == 0));
    }

    #[test]
    fn test_histogram_zero_bins() {
        let hist = violation_histogram(&[1.0, 2.0], 0);
        assert!(hist.is_empty());
    }

    #[test]
    fn test_histogram_count_preserved() {
        let viols: Vec<f64> = (0..20).map(|i| i as f64 * 0.1).collect();
        let hist = violation_histogram(&viols, 5);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 20);
    }

    // --- max_violation ---

    #[test]
    fn test_max_violation_basic() {
        let v = vec![0.1, -0.5, 0.3, -0.2];
        assert!((max_violation(&v) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_max_violation_empty() {
        assert_eq!(max_violation(&[]), 0.0);
    }

    #[test]
    fn test_max_violation_all_negative() {
        let v = vec![-0.3, -0.7, -0.1];
        assert!((max_violation(&v) - 0.7).abs() < 1e-12);
    }

    // --- rms_violation ---

    #[test]
    fn test_rms_violation_basic() {
        let v = vec![3.0, 4.0];
        // rms = sqrt((9+16)/2) = sqrt(12.5)
        let expected = (12.5_f64).sqrt();
        assert!((rms_violation(&v) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_rms_violation_empty() {
        assert_eq!(rms_violation(&[]), 0.0);
    }

    #[test]
    fn test_rms_violation_zero() {
        let v = vec![0.0, 0.0, 0.0];
        assert_eq!(rms_violation(&v), 0.0);
    }

    #[test]
    fn test_rms_violation_non_negative() {
        let v = vec![-1.0, -2.0, -3.0];
        assert!(rms_violation(&v) > 0.0);
    }

    // --- color_by_violation ---

    #[test]
    fn test_color_zero_violation_is_green() {
        let c = color_by_violation(0.0, 1.0);
        assert!(c[1] > 0.99 && c[0] < 0.01, "zero violation → green");
    }

    #[test]
    fn test_color_max_violation_is_red() {
        let c = color_by_violation(1.0, 1.0);
        assert!(c[0] > 0.99 && c[1] < 0.01, "max violation → red");
    }

    #[test]
    fn test_color_channels_in_range() {
        for v in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let c = color_by_violation(v, 1.0);
            for ch in c {
                assert!((0.0..=1.0).contains(&ch), "channel out of range");
            }
        }
    }

    #[test]
    fn test_color_zero_max_viol() {
        let c = color_by_violation(0.5, 0.0);
        // Should return green (safe fallback)
        assert!(c[1] > 0.99);
    }

    // --- constraint_force_magnitude ---

    #[test]
    fn test_force_magnitude_basic() {
        let f = constraint_force_magnitude(0.1, 1000.0);
        assert!((f - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_force_magnitude_negative_violation() {
        let f = constraint_force_magnitude(-0.2, 500.0);
        assert!(
            (f - 100.0).abs() < 1e-10,
            "absolute value of violation used"
        );
    }

    #[test]
    fn test_force_magnitude_zero() {
        assert_eq!(constraint_force_magnitude(0.0, 1000.0), 0.0);
    }

    // --- energy_density_map ---

    #[test]
    fn test_energy_density_length() {
        let positions = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let constraints = vec![(0usize, 1usize, 1.0_f64)];
        let e = energy_density_map(&positions, &constraints, 1000.0);
        assert_eq!(e.len(), 2);
    }

    #[test]
    fn test_energy_density_symmetric() {
        let positions = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0)];
        let e = energy_density_map(&positions, &constraints, 1000.0);
        assert!((e[0] - e[1]).abs() < 1e-12, "energy split equally");
    }

    #[test]
    fn test_energy_density_at_rest_zero() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0)];
        let e = energy_density_map(&positions, &constraints, 1000.0);
        assert!(e[0].abs() < 1e-12 && e[1].abs() < 1e-12);
    }

    #[test]
    fn test_energy_density_positive() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let constraints = vec![(0, 1, 1.0)];
        let e = energy_density_map(&positions, &constraints, 1000.0);
        assert!(e[0] > 0.0 && e[1] > 0.0);
    }

    // --- convergence_metric ---

    #[test]
    fn test_convergence_metric_decreasing() {
        let history = vec![
            vec![1.0, 1.0, 1.0],
            vec![0.5, 0.5, 0.5],
            vec![0.1, 0.1, 0.1],
        ];
        let c = convergence_metric(&history);
        assert_eq!(c.len(), 3);
        assert!(c[0] > c[1] && c[1] > c[2], "should decrease");
    }

    #[test]
    fn test_convergence_metric_empty_history() {
        let c = convergence_metric(&[]);
        assert!(c.is_empty());
    }

    #[test]
    fn test_convergence_metric_single_frame() {
        let history = vec![vec![3.0, 4.0]];
        let c = convergence_metric(&history);
        assert_eq!(c.len(), 1);
        let expected = (12.5_f64).sqrt();
        assert!((c[0] - expected).abs() < 1e-12);
    }

    // --- ConstraintVisualization ---

    #[test]
    fn test_constraint_visualization_push() {
        let mut cv = ConstraintVisualization::new(0);
        let info = ConstraintDebugInfo::new(0, 1, 1.0, 1.5, 1000.0);
        cv.push(info);
        assert_eq!(cv.infos.len(), 1);
        assert!((cv.infos[0].violation - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_constraint_debug_info_violation_sign() {
        let info = ConstraintDebugInfo::new(0, 1, 1.0, 0.5, 1000.0);
        assert!(
            (info.violation - (-0.5)).abs() < 1e-12,
            "compressed → negative violation"
        );
    }
}
