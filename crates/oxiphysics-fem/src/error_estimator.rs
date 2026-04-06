#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! A posteriori error estimation for adaptive FEM.
//!
//! Implements the Zienkiewicz-Zhu (ZZ) patch recovery technique, element-wise
//! error indicators, global error norms, and mesh refinement strategies
//! including Dörfler (bulk) marking.

/// Strategy for adaptive mesh refinement.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum RefineStrategy {
    /// Refine all elements uniformly.
    Uniform,
    /// Refine the top fraction of elements by error indicator.
    FixedFraction(f64),
    /// Dörfler (bulk) criterion: mark fewest elements covering given fraction of total error squared.
    Dorflfer(f64),
}

/// Zienkiewicz-Zhu patch recovery error estimator.
///
/// Stores recovered (superconvergent) nodal stress values obtained by
/// least-squares patch fitting, and element error indicators.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ZzErrorEstimator {
    /// Recovered nodal stress components. Shape: `[n_nodes][n_stress_components]`.
    pub recovered_stress: Vec<Vec<f64>>,
    /// Element error indicators (one per element).
    pub element_indicators: Vec<f64>,
    /// Global L2 error norm.
    pub global_error: f64,
    /// Number of stress components (e.g. 3 for 2D: σ_xx, σ_yy, σ_xy).
    pub n_stress_components: usize,
}

impl ZzErrorEstimator {
    /// Create a new ZZ error estimator.
    ///
    /// # Arguments
    /// * `n_nodes` - total number of mesh nodes
    /// * `n_elements` - total number of mesh elements
    /// * `n_stress_components` - number of stress components (3 for 2D, 6 for 3D)
    pub fn new(n_nodes: usize, n_elements: usize, n_stress_components: usize) -> Self {
        Self {
            recovered_stress: vec![vec![0.0; n_stress_components]; n_nodes],
            element_indicators: vec![0.0; n_elements],
            global_error: 0.0,
            n_stress_components,
        }
    }

    /// Number of nodes in the estimator.
    pub fn n_nodes(&self) -> usize {
        self.recovered_stress.len()
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.element_indicators.len()
    }
}

/// Superconvergent stress recovery at nodes using ZZ patch averaging.
///
/// For each node, averages the raw (Gauss-point) stress values from all
/// surrounding elements in its patch. The recovered values have higher
/// accuracy than the raw element stresses.
///
/// # Arguments
/// * `node_patches` - for each node, list of element indices in its patch
/// * `raw_stress` - raw element stress values, shape `[n_elements][n_stress_components]`
/// * `n_nodes` - total node count
/// * `n_stress_components` - stress component count
///
/// Returns recovered nodal stress array, shape `[n_nodes][n_stress_components]`.
pub fn stress_patch_recovery(
    node_patches: &[Vec<usize>],
    raw_stress: &[Vec<f64>],
    n_nodes: usize,
    n_stress_components: usize,
) -> Vec<Vec<f64>> {
    assert_eq!(node_patches.len(), n_nodes);
    let mut recovered = vec![vec![0.0f64; n_stress_components]; n_nodes];
    for (node, patch) in node_patches.iter().enumerate() {
        if patch.is_empty() {
            continue;
        }
        let n_patch = patch.len() as f64;
        for &elem in patch {
            for c in 0..n_stress_components {
                recovered[node][c] += raw_stress[elem][c] / n_patch;
            }
        }
    }
    recovered
}

/// Compute local element error indicator using ZZ recovery.
///
/// The element error indicator `eta_e` is the L2 norm of the difference
/// between recovered nodal stress (interpolated to element) and the raw
/// element stress:
/// `eta_e^2 = volume_e * || sigma_recovered - sigma_raw ||^2`
///
/// # Arguments
/// * `elem_idx` - element index
/// * `raw_stress` - raw element stress, shape `[n_stress_components]`
/// * `recovered_stress_at_elem` - recovered stress interpolated to element centroid
/// * `elem_volume` - element volume (or area in 2D)
///
/// Returns the element error indicator `eta_e`.
pub fn element_error_indicator(
    _elem_idx: usize,
    raw_stress: &[f64],
    recovered_stress_at_elem: &[f64],
    elem_volume: f64,
) -> f64 {
    assert_eq!(raw_stress.len(), recovered_stress_at_elem.len());
    let sq_diff: f64 = raw_stress
        .iter()
        .zip(recovered_stress_at_elem.iter())
        .map(|(&r, &rec)| (rec - r).powi(2))
        .sum();
    (elem_volume * sq_diff).sqrt()
}

/// Compute the global L2 error norm from element indicators.
///
/// `||e||_{L2} = sqrt( sum_e eta_e^2 )`
///
/// # Arguments
/// * `indicators` - element error indicators (one per element)
pub fn global_error_norm(indicators: &[f64]) -> f64 {
    let sum_sq: f64 = indicators.iter().map(|&e| e * e).sum();
    sum_sq.sqrt()
}

/// Mark elements for refinement using the Dörfler (bulk) criterion.
///
/// Selects the minimal subset of elements such that the sum of their squared
/// error indicators is at least `theta` times the total squared error.
/// Returns a boolean mask: `true` means the element is marked for refinement.
///
/// # Arguments
/// * `indicators` - element error indicators (one per element)
/// * `theta` - bulk criterion parameter in `(0, 1]`, typically 0.5–0.9
pub fn dorflfer_marking(indicators: &[f64], theta: f64) -> Vec<bool> {
    let n = indicators.len();
    let total_sq: f64 = indicators.iter().map(|&e| e * e).sum();
    let target = theta * total_sq;

    // Sort element indices by indicator descending
    let mut sorted_idx: Vec<usize> = (0..n).collect();
    sorted_idx.sort_by(|&a, &b| {
        indicators[b]
            .partial_cmp(&indicators[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut marked = vec![false; n];
    let mut accumulated = 0.0;
    for idx in sorted_idx {
        if accumulated >= target {
            break;
        }
        marked[idx] = true;
        accumulated += indicators[idx] * indicators[idx];
    }
    marked
}

/// Apply a refinement strategy to produce a marked element set.
///
/// Returns a boolean vector indicating which elements to refine.
///
/// # Arguments
/// * `indicators` - element error indicators
/// * `strategy` - refinement strategy
pub fn apply_refine_strategy(indicators: &[f64], strategy: &RefineStrategy) -> Vec<bool> {
    match strategy {
        RefineStrategy::Uniform => vec![true; indicators.len()],
        RefineStrategy::FixedFraction(frac) => {
            let n = indicators.len();
            let n_mark = ((frac * n as f64).ceil() as usize).min(n);
            let mut sorted_idx: Vec<usize> = (0..n).collect();
            sorted_idx.sort_by(|&a, &b| {
                indicators[b]
                    .partial_cmp(&indicators[a])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut marked = vec![false; n];
            for &idx in &sorted_idx[..n_mark] {
                marked[idx] = true;
            }
            marked
        }
        RefineStrategy::Dorflfer(theta) => dorflfer_marking(indicators, *theta),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── ZzErrorEstimator ───────────────────────────────────────────────────

    #[test]
    fn test_zz_new_sizes() {
        let est = ZzErrorEstimator::new(10, 8, 3);
        assert_eq!(est.n_nodes(), 10);
        assert_eq!(est.n_elements(), 8);
    }

    #[test]
    fn test_zz_recovered_stress_initialized_zero() {
        let est = ZzErrorEstimator::new(5, 4, 3);
        for row in &est.recovered_stress {
            for &v in row {
                assert!(v.abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_zz_element_indicators_initialized_zero() {
        let est = ZzErrorEstimator::new(5, 4, 3);
        for &ind in &est.element_indicators {
            assert!(ind.abs() < 1e-15);
        }
    }

    #[test]
    fn test_zz_global_error_initialized_zero() {
        let est = ZzErrorEstimator::new(5, 4, 3);
        assert!(est.global_error.abs() < 1e-15);
    }

    #[test]
    fn test_zz_3d_components() {
        let est = ZzErrorEstimator::new(4, 2, 6);
        assert_eq!(est.n_stress_components, 6);
        assert_eq!(est.recovered_stress[0].len(), 6);
    }

    // ── stress_patch_recovery ──────────────────────────────────────────────

    #[test]
    fn test_patch_recovery_single_patch_elem() {
        let patches = vec![vec![0usize]];
        let raw = vec![vec![1.0, 2.0, 3.0]];
        let rec = stress_patch_recovery(&patches, &raw, 1, 3);
        assert_eq!(rec.len(), 1);
        assert!((rec[0][0] - 1.0).abs() < 1e-14);
        assert!((rec[0][1] - 2.0).abs() < 1e-14);
        assert!((rec[0][2] - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_patch_recovery_average_two_elems() {
        let patches = vec![vec![0, 1]];
        let raw = vec![vec![2.0, 4.0], vec![4.0, 0.0]];
        let rec = stress_patch_recovery(&patches, &raw, 1, 2);
        assert!((rec[0][0] - 3.0).abs() < 1e-14);
        assert!((rec[0][1] - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_patch_recovery_empty_patch() {
        let patches = vec![vec![]];
        let raw: Vec<Vec<f64>> = vec![];
        let rec = stress_patch_recovery(&patches, &raw, 1, 2);
        assert!((rec[0][0]).abs() < 1e-15);
        assert!((rec[0][1]).abs() < 1e-15);
    }

    #[test]
    fn test_patch_recovery_multiple_nodes() {
        let patches = vec![vec![0], vec![1], vec![0, 1]];
        let raw = vec![vec![10.0], vec![20.0]];
        let rec = stress_patch_recovery(&patches, &raw, 3, 1);
        assert!((rec[0][0] - 10.0).abs() < 1e-14);
        assert!((rec[1][0] - 20.0).abs() < 1e-14);
        assert!((rec[2][0] - 15.0).abs() < 1e-14);
    }

    #[test]
    fn test_patch_recovery_zero_raw() {
        let patches = vec![vec![0, 1, 2]];
        let raw = vec![vec![0.0], vec![0.0], vec![0.0]];
        let rec = stress_patch_recovery(&patches, &raw, 1, 1);
        assert!(rec[0][0].abs() < 1e-14);
    }

    // ── element_error_indicator ────────────────────────────────────────────

    #[test]
    fn test_element_error_zero_difference() {
        let raw = vec![1.0, 2.0, 3.0];
        let rec = vec![1.0, 2.0, 3.0];
        let eta = element_error_indicator(0, &raw, &rec, 1.0);
        assert!(eta.abs() < 1e-14);
    }

    #[test]
    fn test_element_error_nonzero() {
        let raw = vec![0.0];
        let rec = vec![1.0];
        let eta = element_error_indicator(0, &raw, &rec, 1.0);
        assert!((eta - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_element_error_scales_with_volume() {
        let raw = vec![0.0];
        let rec = vec![1.0];
        let eta1 = element_error_indicator(0, &raw, &rec, 1.0);
        let eta2 = element_error_indicator(0, &raw, &rec, 4.0);
        assert!((eta2 / eta1 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_error_non_negative() {
        let raw = vec![5.0, -3.0, 1.0];
        let rec = vec![2.0, 1.0, -4.0];
        let eta = element_error_indicator(0, &raw, &rec, 0.5);
        assert!(eta >= 0.0);
    }

    #[test]
    fn test_element_error_3_components() {
        // (3-1)^2 + (4-2)^2 + (5-3)^2 = 4+4+4=12, *vol=1 => sqrt(12)
        let raw = vec![1.0, 2.0, 3.0];
        let rec = vec![3.0, 4.0, 5.0];
        let eta = element_error_indicator(0, &raw, &rec, 1.0);
        assert!((eta - 12.0f64.sqrt()).abs() < 1e-12);
    }

    // ── global_error_norm ──────────────────────────────────────────────────

    #[test]
    fn test_global_error_norm_zero() {
        let ind = vec![0.0, 0.0, 0.0];
        assert!(global_error_norm(&ind).abs() < 1e-14);
    }

    #[test]
    fn test_global_error_norm_single() {
        let ind = vec![3.0];
        assert!((global_error_norm(&ind) - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_global_error_norm_pythagorean() {
        let ind = vec![3.0, 4.0];
        assert!((global_error_norm(&ind) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_global_error_norm_all_ones() {
        let ind = vec![1.0; 9];
        assert!((global_error_norm(&ind) - 3.0).abs() < 1e-12);
    }

    // ── dorflfer_marking ───────────────────────────────────────────────────

    #[test]
    fn test_dorfler_theta_one_marks_all() {
        let ind = vec![1.0, 2.0, 3.0, 0.5];
        let marked = dorflfer_marking(&ind, 1.0);
        assert!(marked.iter().all(|&m| m));
    }

    #[test]
    fn test_dorfler_marks_largest_first() {
        let ind = vec![1.0, 10.0, 2.0];
        // theta=0.5: need 50% of total^2 = 0.5*(1+100+4)=52.5; largest=10, 10^2=100 > 52.5
        let marked = dorflfer_marking(&ind, 0.5);
        assert!(marked[1]); // largest must be marked
    }

    #[test]
    fn test_dorfler_count_reasonable() {
        let ind = vec![1.0, 1.0, 1.0, 1.0];
        // equal indicators, theta=0.5: need 50% of 4 = 2; marking 2 suffices
        let marked = dorflfer_marking(&ind, 0.5);
        let n_marked: usize = marked.iter().filter(|&&m| m).count();
        assert!(n_marked >= 1 && n_marked <= 4);
    }

    #[test]
    fn test_dorfler_returns_correct_length() {
        let ind = vec![3.0, 1.0, 4.0, 1.0, 5.0];
        let marked = dorflfer_marking(&ind, 0.7);
        assert_eq!(marked.len(), 5);
    }

    #[test]
    fn test_dorfler_all_zero_indicators() {
        let ind = vec![0.0, 0.0, 0.0];
        // total_sq = 0, target = 0, accumulated starts >= 0 immediately => no elements marked
        let marked = dorflfer_marking(&ind, 0.5);
        assert_eq!(marked.len(), 3);
    }

    // ── RefineStrategy / apply_refine_strategy ─────────────────────────────

    #[test]
    fn test_uniform_marks_all() {
        let ind = vec![1.0, 0.5, 2.0];
        let marked = apply_refine_strategy(&ind, &RefineStrategy::Uniform);
        assert!(marked.iter().all(|&m| m));
    }

    #[test]
    fn test_fixed_fraction_marks_correct_count() {
        let ind = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let marked = apply_refine_strategy(&ind, &RefineStrategy::FixedFraction(0.4));
        let n_marked: usize = marked.iter().filter(|&&m| m).count();
        assert_eq!(n_marked, 2); // ceil(0.4*5)=2
    }

    #[test]
    fn test_fixed_fraction_marks_largest() {
        let ind = vec![1.0, 5.0, 2.0, 4.0];
        let marked = apply_refine_strategy(&ind, &RefineStrategy::FixedFraction(0.25));
        // ceil(0.25*4)=1, should mark largest=idx 1
        assert!(marked[1]);
    }

    #[test]
    fn test_dorfler_strategy_wrapper() {
        let ind = vec![1.0, 10.0, 2.0];
        let marked = apply_refine_strategy(&ind, &RefineStrategy::Dorflfer(0.5));
        assert!(marked[1]);
    }

    #[test]
    fn test_refine_strategy_clone() {
        let s = RefineStrategy::FixedFraction(0.3);
        let _s2 = s.clone();
    }

    #[test]
    fn test_refine_strategy_eq() {
        assert_eq!(RefineStrategy::Uniform, RefineStrategy::Uniform);
        assert_ne!(
            RefineStrategy::FixedFraction(0.3),
            RefineStrategy::FixedFraction(0.7)
        );
    }

    #[test]
    fn test_fixed_fraction_full() {
        let ind = vec![1.0, 2.0];
        let marked = apply_refine_strategy(&ind, &RefineStrategy::FixedFraction(1.0));
        assert!(marked.iter().all(|&m| m));
    }
}
