// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ewald parameters, erfc approximation, optimizer, accuracy budget, and configuration.

// ---------------------------------------------------------------------------
// Physical constant
// ---------------------------------------------------------------------------

/// Coulomb constant in GROMACS-compatible units: kJ*angstrom*mol^-1*e^-2.
///
/// Equivalent to k_e / kBT at 300 K scaled so that energy is in kJ/mol
/// and distances are in angstrom.  Value: 138.935 kJ*angstrom/(mol*e^2).
pub const COULOMB_K: f64 = 138.935;

// ---------------------------------------------------------------------------
// erfc approximation (Abramowitz & Stegun 7.1.26)
// ---------------------------------------------------------------------------

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
///
/// Maximum absolute error ~ 1.5e-7.
pub fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

// ---------------------------------------------------------------------------
// EwaldParams
// ---------------------------------------------------------------------------

/// Parameters controlling the Ewald splitting between real and reciprocal space.
#[derive(Debug, Clone)]
pub struct EwaldParams {
    /// Splitting parameter alpha (angstrom^-1).
    ///
    /// Larger alpha -> shorter-range real-space sum, more work in reciprocal space.
    pub alpha: f64,
    /// Real-space cutoff distance (angstrom).
    pub r_cutoff: f64,
    /// Reciprocal-space cutoff |k| (angstrom^-1).
    pub k_cutoff: f64,
    /// Relative dielectric constant (dimensionless).
    pub epsilon_r: f64,
}

/// Maximum number of k-vectors per dimension for reciprocal-space sums.
pub const DEFAULT_K_MAX: i32 = 5;

impl EwaldParams {
    /// Create [`EwaldParams`] with common defaults given a real-space cutoff.
    ///
    /// Uses `alpha = 5 / r_cutoff` and `k_cutoff = 2*pi*alpha`.
    pub fn new(r_cutoff: f64) -> Self {
        let alpha = 5.0 / r_cutoff;
        let k_cutoff = 2.0 * std::f64::consts::PI * alpha;
        Self {
            alpha,
            r_cutoff,
            k_cutoff,
            epsilon_r: 1.0,
        }
    }

    /// Compute optimal alpha for a given real-space cutoff and desired accuracy.
    ///
    /// Uses the heuristic `alpha = sqrt(-ln(accuracy)) / r_cut`.
    pub fn optimal_alpha(r_cut: f64, accuracy: f64) -> f64 {
        assert!(r_cut > 0.0, "r_cut must be positive");
        assert!(
            accuracy > 0.0 && accuracy < 1.0,
            "accuracy must be in (0,1)"
        );
        (-accuracy.ln()).sqrt() / r_cut
    }

    /// Compute optimal parameters (alpha, k_max) for given accuracy and box size.
    ///
    /// Returns `(alpha, k_max)` where k_max is the number of k-vectors
    /// per dimension needed to achieve the desired accuracy.
    pub fn optimize(r_cut: f64, box_length: f64, accuracy: f64) -> (f64, i32) {
        let alpha = Self::optimal_alpha(r_cut, accuracy);
        // k_max such that exp(-(pi*k_max/(alpha*L))^2) < accuracy
        let ratio = alpha * box_length / std::f64::consts::PI;
        let k_max = (ratio * (-accuracy.ln()).sqrt()).ceil() as i32;
        let k_max = k_max.max(1);
        (alpha, k_max)
    }

    /// Set the dielectric constant.
    pub fn with_epsilon_r(mut self, eps: f64) -> Self {
        self.epsilon_r = eps;
        self
    }

    /// Ewald real-space pair energy (kJ mol^-1) for charges `q1` and `q2`
    /// separated by distance `r` (angstrom).
    ///
    /// ```text
    /// E = COULOMB_K / epsilon_r * q1*q2 * erfc(alpha*r) / r
    /// ```
    pub fn real_space_energy(&self, q1: f64, q2: f64, r: f64) -> f64 {
        if r >= self.r_cutoff || r <= 0.0 {
            return 0.0;
        }
        (COULOMB_K / self.epsilon_r) * q1 * q2 * erfc_approx(self.alpha * r) / r
    }

    /// Ewald real-space pair force magnitude (kJ mol^-1 angstrom^-1) for charges
    /// `q1` and `q2` separated by distance `r` (angstrom).
    ///
    /// Positive value -> repulsive (pushes particles apart).
    /// Negative value -> attractive (pulls particles together).
    ///
    /// Derived from `F = -dE/dr`:
    /// ```text
    /// F = COULOMB_K/epsilon_r * q1*q2 * [erfc(alpha*r)/r^2 + 2*alpha*exp(-alpha^2*r^2)/(sqrt(pi)*r)]
    /// ```
    pub fn real_space_force_mag(&self, q1: f64, q2: f64, r: f64) -> f64 {
        if r >= self.r_cutoff || r <= 0.0 {
            return 0.0;
        }
        let ar = self.alpha * r;
        let erfc_val = erfc_approx(ar);
        let exp_val = (-ar * ar).exp();
        let prefactor = COULOMB_K / self.epsilon_r * q1 * q2;
        // F = -dE/dr; erfc term + Gaussian correction
        prefactor
            * (erfc_val / (r * r) + 2.0 * self.alpha * exp_val / (std::f64::consts::PI.sqrt() * r))
    }

    /// Estimate the real-space error for the Ewald sum.
    ///
    /// ```text
    /// Delta_F_real ~ q^2 * erfc(alpha * r_cut) / r_cut^2
    /// ```
    pub fn real_space_error_estimate(&self, charge_sum_sq: f64) -> f64 {
        let erfc_val = erfc_approx(self.alpha * self.r_cutoff);
        charge_sum_sq * erfc_val / (self.r_cutoff * self.r_cutoff)
    }

    /// Estimate the reciprocal-space error for the Ewald sum.
    ///
    /// ```text
    /// Delta_F_recip ~ q^2 * alpha * exp(-(pi*k_max/(alpha*L))^2)
    /// ```
    pub fn reciprocal_space_error_estimate(
        &self,
        charge_sum_sq: f64,
        box_length: f64,
        k_max: i32,
    ) -> f64 {
        let arg = std::f64::consts::PI * k_max as f64 / (self.alpha * box_length);
        charge_sum_sq * self.alpha * (-arg * arg).exp()
    }
}

// ---------------------------------------------------------------------------
// EwaldSumConfig — compact configuration struct
// ---------------------------------------------------------------------------

/// Compact configuration for a complete Ewald sum.
///
/// Bundles the splitting parameter, k-vector cutoff, and box volume
/// into a single struct for easy parameter passing.
#[derive(Debug, Clone)]
pub struct EwaldSumConfig {
    /// Ewald splitting parameter alpha (angstrom^-1).
    pub alpha: f64,
    /// Maximum k-vector index per dimension (number of k-vectors per axis = 2*k_max+1).
    pub k_max: usize,
    /// Box volume (angstrom^3).
    pub volume: f64,
    /// Cubic box side length (angstrom; volume = box_len^3 assumed).
    pub box_len: f64,
}

impl EwaldSumConfig {
    /// Create a new config for a cubic box.
    ///
    /// `box_len` is the side length in angstrom; volume is computed automatically.
    pub fn new(alpha: f64, k_max: usize, box_len: f64) -> Self {
        assert!(alpha > 0.0, "alpha must be positive");
        assert!(box_len > 0.0, "box_len must be positive");
        Self {
            alpha,
            k_max,
            volume: box_len * box_len * box_len,
            box_len,
        }
    }

    /// Create from explicit volume (non-cubic boxes).
    pub fn with_volume(alpha: f64, k_max: usize, box_len: f64, volume: f64) -> Self {
        assert!(alpha > 0.0, "alpha must be positive");
        assert!(box_len > 0.0, "box_len must be positive");
        assert!(volume > 0.0, "volume must be positive");
        Self {
            alpha,
            k_max,
            volume,
            box_len,
        }
    }
}

// ---------------------------------------------------------------------------
// EwaldAccuracyBudget — choose alpha and k_max for a given error budget
// ---------------------------------------------------------------------------

/// Estimates and manages Ewald accuracy parameters for a given error budget.
///
/// Given a target force error, selects the optimal combination of `alpha`
/// and `k_max` for a cubic box, and provides estimates of the resulting
/// real-space and reciprocal-space errors.
#[derive(Debug, Clone)]
pub struct EwaldAccuracyBudget {
    /// Optimal Ewald alpha (angstrom⁻¹).
    pub alpha: f64,
    /// Number of k-vectors per dimension.
    pub k_max: i32,
    /// Estimated real-space truncation error.
    pub real_error: f64,
    /// Estimated reciprocal-space truncation error.
    pub recip_error: f64,
}

impl EwaldAccuracyBudget {
    /// Compute an accuracy budget for a given system.
    ///
    /// # Arguments
    /// * `r_cut`       - real-space cutoff (Å).
    /// * `box_len`     - cubic box side length (Å).
    /// * `charge_sq`   - sum of squared charges Σ qᵢ² (e²).
    /// * `target_err`  - desired force error (kJ mol⁻¹ Å⁻¹).
    pub fn compute(r_cut: f64, box_len: f64, charge_sq: f64, target_err: f64) -> Self {
        let (alpha, k_max) = EwaldParams::optimize(r_cut, box_len, target_err.max(1e-15));
        let params = EwaldParams {
            alpha,
            r_cutoff: r_cut,
            k_cutoff: 2.0 * std::f64::consts::PI * alpha,
            epsilon_r: 1.0,
        };
        let real_error = params.real_space_error_estimate(charge_sq);
        let recip_error = params.reciprocal_space_error_estimate(charge_sq, box_len, k_max);
        Self {
            alpha,
            k_max,
            real_error,
            recip_error,
        }
    }

    /// Total estimated error (upper bound).
    pub fn total_error(&self) -> f64 {
        self.real_error + self.recip_error
    }

    /// Whether the budget meets the target error.
    pub fn meets_target(&self, target: f64) -> bool {
        self.total_error() <= target
    }
}

// ---------------------------------------------------------------------------
// ErfTable — lookup table for erf/erfc (performance optimization demo)
// ---------------------------------------------------------------------------

/// Precomputed erfc lookup table for fast evaluation.
///
/// Tabulates `erfc(x)` for x ∈ \[0, x_max\] with `n` equal-spaced points.
/// Linear interpolation is used for values between table entries.
#[derive(Debug, Clone)]
pub struct ErfcTable {
    /// Maximum x value in the table.
    pub x_max: f64,
    /// Spacing between table entries.
    pub dx: f64,
    /// Tabulated erfc values.
    pub table: Vec<f64>,
}

impl ErfcTable {
    /// Build a new erfc lookup table.
    ///
    /// # Arguments
    /// * `x_max` - maximum argument (beyond this, erfc is treated as 0).
    /// * `n`     - number of table points (must be >= 2).
    pub fn new(x_max: f64, n: usize) -> Self {
        assert!(n >= 2, "table must have at least 2 entries");
        assert!(x_max > 0.0, "x_max must be positive");
        let dx = x_max / (n as f64 - 1.0);
        let table: Vec<f64> = (0..n).map(|i| erfc_approx(i as f64 * dx)).collect();
        Self { x_max, dx, table }
    }

    /// Look up erfc(x) using linear interpolation.
    ///
    /// Returns 1.0 for x ≤ 0, 0.0 for x > x_max.
    pub fn eval(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 1.0;
        }
        if x >= self.x_max {
            return 0.0;
        }
        let idx_f = x / self.dx;
        let idx = idx_f as usize;
        let frac = idx_f - idx as f64;
        let lo = self.table[idx];
        let hi = if idx + 1 < self.table.len() {
            self.table[idx + 1]
        } else {
            0.0
        };
        lo + frac * (hi - lo)
    }

    /// Number of table entries.
    pub fn size(&self) -> usize {
        self.table.len()
    }
}

// ---------------------------------------------------------------------------
// AlphaOptimizer — automatically choose alpha given N, V, rc, accuracy
// ---------------------------------------------------------------------------

/// Finds the optimal Ewald splitting parameter `alpha` for a given system size
/// and desired force accuracy.
///
/// The optimal alpha minimises the total computation cost while keeping
/// the real-space truncation error below `accuracy`.
#[derive(Debug, Clone, Copy)]
pub struct AlphaOptimizer {
    /// Real-space cutoff (Å).
    pub r_cut: f64,
    /// Desired relative force accuracy.
    pub accuracy: f64,
    /// Box length (Å), assumed cubic.
    pub box_len: f64,
}

impl AlphaOptimizer {
    /// Create an optimizer with given cutoff, accuracy, and box size.
    pub fn new(r_cut: f64, accuracy: f64, box_len: f64) -> Self {
        Self {
            r_cut,
            accuracy,
            box_len,
        }
    }

    /// Compute the optimal alpha using the heuristic:
    /// `alpha = sqrt(-ln(accuracy)) / r_cut`.
    pub fn optimal_alpha(&self) -> f64 {
        (-self.accuracy.ln()).max(0.0_f64).sqrt() / self.r_cut
    }

    /// Estimate the required `k_max` for a given `alpha` and box length.
    ///
    /// `k_max = ceil(alpha * box_len * sqrt(-ln(accuracy)) / π)`.
    pub fn required_k_max(&self, alpha: f64) -> i32 {
        let arg = (-self.accuracy.ln()).max(0.0_f64).sqrt();
        let k_max_f = alpha * self.box_len * arg / std::f64::consts::PI;
        (k_max_f.ceil() as i32).max(1)
    }

    /// Build an [`EwaldParams`] struct with optimal parameters.
    pub fn build_params(&self) -> EwaldParams {
        let alpha = self.optimal_alpha();
        let k_max = self.required_k_max(alpha);
        EwaldParams {
            alpha,
            r_cutoff: self.r_cut,
            k_cutoff: 2.0_f64 * std::f64::consts::PI * k_max as f64 / self.box_len,
            epsilon_r: 1.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erfc_approx() {
        // erfc(0) = 1
        let v0 = erfc_approx(0.0);
        assert!((v0 - 1.0).abs() < 1e-6, "erfc(0) = {v0}");

        // erfc(large) ~ 0
        let v_large = erfc_approx(10.0);
        assert!(v_large < 1e-6, "erfc(10) should be ~0, got {v_large}");

        // erfc(-x) = 2 - erfc(x)
        let x = 1.5;
        let lhs = erfc_approx(-x);
        let rhs = 2.0 - erfc_approx(x);
        assert!((lhs - rhs).abs() < 1e-10, "erfc symmetry: {lhs} vs {rhs}");
    }

    #[test]
    fn test_ewald_sum_erfc_approx_at_zero() {
        let v = erfc_approx(0.0);
        assert!(
            (v - 1.0).abs() < 1e-6,
            "erfc_approx(0) should be 1.0, got {v}"
        );
    }

    #[test]
    fn test_ewald_sum_config_volume() {
        let cfg = EwaldSumConfig::new(0.5, 5, 20.0);
        assert!(
            (cfg.volume - 8000.0).abs() < 1e-8,
            "volume = {}, expected 8000",
            cfg.volume
        );
    }

    #[test]
    fn test_ewald_sum_config_with_volume() {
        let cfg = EwaldSumConfig::with_volume(0.3, 3, 15.0, 1000.0);
        assert!((cfg.volume - 1000.0).abs() < 1e-8);
        assert!((cfg.alpha - 0.3).abs() < 1e-14);
        assert_eq!(cfg.k_max, 3);
    }

    #[test]
    fn test_ewald_accuracy_budget_alpha_positive() {
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 2.0, 1e-4);
        assert!(budget.alpha > 0.0, "alpha should be positive");
        assert!(budget.k_max >= 1, "k_max should be >= 1");
    }

    #[test]
    fn test_ewald_accuracy_budget_errors_nonneg() {
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 2.0, 1e-4);
        assert!(budget.real_error >= 0.0);
        assert!(budget.recip_error >= 0.0);
        assert!(budget.total_error() >= 0.0);
    }

    #[test]
    fn test_ewald_accuracy_budget_meets_lenient_target() {
        // With a very lenient target, the budget should meet it
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 2.0, 1e-2);
        assert!(budget.meets_target(1.0), "should meet target of 1.0");
    }

    #[test]
    fn test_erfc_table_at_zero() {
        let table = ErfcTable::new(8.0, 1000);
        let v = table.eval(0.0);
        assert!(
            (v - 1.0).abs() < 1e-6,
            "erfc(0) from table should be ~1.0, got {v}"
        );
    }

    #[test]
    fn test_erfc_table_at_large_x() {
        let table = ErfcTable::new(8.0, 1000);
        let v = table.eval(7.9);
        assert!(v < 1e-4, "erfc(7.9) should be ~0, got {v}");
        let v_beyond = table.eval(9.0);
        assert_eq!(v_beyond, 0.0, "beyond x_max should return 0");
    }

    #[test]
    fn test_erfc_table_matches_approx() {
        let table = ErfcTable::new(5.0, 10000);
        // Check a few reference values
        for &x in &[0.5, 1.0, 1.5, 2.0, 3.0] {
            let t_val = table.eval(x);
            let ref_val = erfc_approx(x);
            let err = (t_val - ref_val).abs();
            assert!(
                err < 1e-3,
                "table erfc({x}) = {t_val} vs approx {ref_val}, err = {err}"
            );
        }
    }

    #[test]
    fn test_erfc_table_size() {
        let table = ErfcTable::new(6.0, 200);
        assert_eq!(table.size(), 200);
    }

    #[test]
    fn test_erfc_table_monotone_decreasing() {
        let table = ErfcTable::new(6.0, 100);
        for i in 1..table.table.len() {
            assert!(
                table.table[i] <= table.table[i - 1] + 1e-14,
                "erfc table should be monotone decreasing at index {i}"
            );
        }
    }

    #[test]
    fn test_ewald_accuracy_budget_alpha_positive_v2() {
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 4.0, 1e-4);
        assert!(budget.alpha > 0.0);
    }

    #[test]
    fn test_ewald_accuracy_budget_total_error_finite_v2() {
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 4.0, 1e-4);
        assert!(budget.total_error().is_finite());
    }

    #[test]
    fn test_ewald_accuracy_budget_k_max_positive() {
        let budget = EwaldAccuracyBudget::compute(10.0, 30.0, 4.0, 1e-4);
        assert!(budget.k_max >= 1);
    }

    #[test]
    fn test_ewald_accuracy_budget_meets_target_loose() {
        // Use a target in (0, 1) as required by optimal_alpha
        let budget = EwaldAccuracyBudget::compute(8.0, 25.0, 1.0, 0.5);
        // meets_target with a very high threshold should always be true
        assert!(budget.meets_target(1e10));
    }

    #[test]
    fn test_alpha_optimizer_optimal_alpha_positive() {
        let opt = AlphaOptimizer::new(10.0, 1e-4, 30.0);
        let alpha = opt.optimal_alpha();
        assert!(alpha > 0.0, "alpha must be positive: {alpha}");
    }

    #[test]
    fn test_alpha_optimizer_required_k_max_positive() {
        let opt = AlphaOptimizer::new(10.0, 1e-4, 30.0);
        let alpha = opt.optimal_alpha();
        let k_max = opt.required_k_max(alpha);
        assert!(k_max >= 1, "k_max must be at least 1: {k_max}");
    }

    #[test]
    fn test_alpha_optimizer_build_params_reasonable_alpha() {
        let opt = AlphaOptimizer::new(10.0, 1e-5, 30.0);
        let params = opt.build_params();
        assert!(
            params.alpha > 0.0 && params.alpha < 2.0,
            "alpha={} should be in (0,2)",
            params.alpha
        );
    }

    #[test]
    fn test_alpha_optimizer_larger_accuracy_gives_smaller_alpha() {
        let opt_tight = AlphaOptimizer::new(10.0, 1e-8, 30.0);
        let opt_loose = AlphaOptimizer::new(10.0, 1e-2, 30.0);
        // tighter accuracy → larger alpha
        assert!(
            opt_tight.optimal_alpha() > opt_loose.optimal_alpha(),
            "tight accuracy should give larger alpha"
        );
    }
}
