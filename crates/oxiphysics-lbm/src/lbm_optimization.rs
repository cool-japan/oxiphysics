// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM optimization, parameter tuning, and performance analysis.
//!
//! Provides tools for:
//! - Stability analysis via von Neumann criterion for BGK
//! - TRT magic parameter optimization
//! - Convergence monitoring via L2 norm of velocity residuals
//! - Automatic omega tuning via bisection for a target Reynolds number
//! - Performance metrics (MLUPS)
//! - Richardson extrapolation error estimation

// ---------------------------------------------------------------------------
// LbmParameters
// ---------------------------------------------------------------------------

/// LBM simulation parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LbmParameters {
    /// BGK relaxation frequency omega (must be in (0, 2) for stability).
    pub omega: f64,
    /// Number of lattice nodes in x direction.
    pub nx: usize,
    /// Number of lattice nodes in y direction.
    pub ny: usize,
    /// Simulation timestep (lattice units, typically 1).
    pub timestep: f64,
    /// Physical length scale (meters per lattice unit).
    pub dx: f64,
    /// Physical time scale (seconds per lattice timestep).
    pub dt: f64,
}

impl LbmParameters {
    /// Create new LBM parameters with default values.
    ///
    /// # Arguments
    /// - `omega`: relaxation frequency
    /// - `nx`, `ny`: lattice dimensions
    pub fn new(omega: f64, nx: usize, ny: usize) -> Self {
        Self {
            omega,
            nx,
            ny,
            timestep: 1.0,
            dx: 1.0,
            dt: 1.0,
        }
    }

    /// Kinematic viscosity in lattice units: nu = cs^2 * (1/omega - 0.5).
    pub fn viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega - 0.5)
    }

    /// Total number of lattice nodes.
    pub fn total_nodes(&self) -> usize {
        self.nx * self.ny
    }
}

// ---------------------------------------------------------------------------
// Stability analysis
// ---------------------------------------------------------------------------

/// Check von Neumann stability for BGK: omega must be in (0, 2).
///
/// Returns `true` if omega is in the stable range (0, 2), exclusive.
#[allow(dead_code)]
pub fn stability_analysis(omega: f64) -> bool {
    omega > 0.0 && omega < 2.0
}

/// Compute the stability margin: distance from omega to the nearest boundary.
///
/// Returns `min(omega, 2.0 - omega)` for omega in (0, 2), else a negative value.
#[allow(dead_code)]
pub fn stability_margin(omega: f64) -> f64 {
    if !stability_analysis(omega) {
        return f64::NEG_INFINITY;
    }
    omega.min(2.0 - omega)
}

// ---------------------------------------------------------------------------
// TRT magic parameter optimization
// ---------------------------------------------------------------------------

/// Compute the TRT "magic parameter" Lambda = tau_plus * tau_minus.
///
/// The optimal value Lambda = 3/16 eliminates the viscosity-dependent
/// numerical slip on no-slip boundaries.
///
/// # Arguments
/// - `omega_plus`: symmetric relaxation frequency (controls viscosity)
///
/// Returns `omega_minus` such that `tau_plus * tau_minus = magic`.
#[allow(dead_code)]
pub fn optimal_relaxation(omega_plus: f64, magic: f64) -> f64 {
    let tau_plus = 1.0 / omega_plus;
    // tau_minus = magic / tau_plus
    let tau_minus = magic / tau_plus;
    1.0 / tau_minus
}

/// Standard TRT magic parameter Lambda = 3/16.
pub const TRT_MAGIC: f64 = 3.0 / 16.0;

// ---------------------------------------------------------------------------
// Convergence criterion
// ---------------------------------------------------------------------------

/// Compute L2 norm of velocity residual between two timesteps.
///
/// # Arguments
/// - `ux_old`, `uy_old`: velocity components at previous step (length n)
/// - `ux_new`, `uy_new`: velocity components at current step (length n)
///
/// Returns the L2 norm normalized by the number of nodes.
#[allow(dead_code)]
pub fn convergence_criterion(
    ux_old: &[f64],
    uy_old: &[f64],
    ux_new: &[f64],
    uy_new: &[f64],
) -> f64 {
    assert_eq!(ux_old.len(), uy_old.len());
    assert_eq!(ux_new.len(), uy_new.len());
    assert_eq!(ux_old.len(), ux_new.len());
    let n = ux_old.len();
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = ux_old
        .iter()
        .zip(uy_old)
        .zip(ux_new.iter().zip(uy_new))
        .map(|((ox, oy), (nx, ny))| {
            let dx = nx - ox;
            let dy = ny - oy;
            dx * dx + dy * dy
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}

// ---------------------------------------------------------------------------
// Auto-tune omega
// ---------------------------------------------------------------------------

/// Auto-tune omega for a target Reynolds number using bisection.
///
/// The Reynolds number is Re = U * L / nu, where nu = cs^2 * (1/omega - 0.5).
///
/// # Arguments
/// - `target_re`: desired Reynolds number
/// - `u_char`: characteristic velocity (lattice units)
/// - `l_char`: characteristic length (lattice units)
/// - `tol`: tolerance for bisection convergence
///
/// Returns the omega value achieving the target Re, clamped to (0.01, 1.99).
#[allow(dead_code)]
pub fn auto_tune_omega(target_re: f64, u_char: f64, l_char: f64, tol: f64) -> f64 {
    // nu = cs2 * (1/omega - 0.5), Re = u * l / nu
    // => 1/omega - 0.5 = u*l / (Re * cs2)
    // => omega = 1 / (u*l / (Re * cs2) + 0.5)
    let cs2 = 1.0 / 3.0;
    let nu_target = u_char * l_char / target_re;
    let tau = nu_target / cs2 + 0.5;
    let omega = 1.0 / tau;

    // Bisection to refine if needed (for generality)
    let mut lo = 0.01_f64;
    let mut hi = 1.99_f64;
    let mut mid = omega.clamp(lo, hi);

    for _ in 0..100 {
        mid = (lo + hi) / 2.0;
        let nu_mid = cs2 * (1.0 / mid - 0.5);
        let re_mid = u_char * l_char / nu_mid.max(1e-30);
        if (re_mid - target_re).abs() < tol {
            break;
        }
        // re increases as nu decreases, nu decreases as omega increases
        if re_mid < target_re {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    mid
}

// ---------------------------------------------------------------------------
// Performance metrics
// ---------------------------------------------------------------------------

/// Compute MLUPS (million lattice updates per second).
///
/// # Arguments
/// - `n_nodes`: total lattice nodes updated
/// - `n_steps`: number of timesteps executed
/// - `elapsed_secs`: wall-clock time in seconds
#[allow(dead_code)]
pub fn performance_metrics(n_nodes: usize, n_steps: usize, elapsed_secs: f64) -> f64 {
    if elapsed_secs <= 0.0 {
        return 0.0;
    }
    (n_nodes as f64 * n_steps as f64) / (elapsed_secs * 1.0e6)
}

// ---------------------------------------------------------------------------
// Grid refinement error (Richardson extrapolation)
// ---------------------------------------------------------------------------

/// Estimate the Richardson extrapolation error between two grid levels.
///
/// Given coarse (`q_coarse`) and fine (`q_fine`) solutions and a refinement
/// ratio `r` (typically 2), the observed order `p`, returns the estimated
/// error of the fine solution:
///
/// `error = (q_fine - q_coarse) / (r^p - 1)`
///
/// # Arguments
/// - `q_coarse`: quantity of interest on coarse grid
/// - `q_fine`:   quantity of interest on fine grid
/// - `r`:        grid refinement ratio (e.g. 2.0)
/// - `p`:        observed convergence order (e.g. 2.0 for second-order)
#[allow(dead_code)]
pub fn grid_refinement_error(q_coarse: f64, q_fine: f64, r: f64, p: f64) -> f64 {
    let factor = r.powf(p) - 1.0;
    if factor.abs() < 1e-30 {
        return 0.0;
    }
    (q_fine - q_coarse) / factor
}

/// Estimate the convergence order from three grid levels.
///
/// Uses Richardson's formula:
/// `p = ln((q_coarse - q_medium) / (q_medium - q_fine)) / ln(r)`
///
/// # Arguments
/// - `q_coarse`, `q_medium`, `q_fine`: quantity at three grid levels
/// - `r`: refinement ratio between consecutive levels
#[allow(dead_code)]
pub fn convergence_order(q_coarse: f64, q_medium: f64, q_fine: f64, r: f64) -> f64 {
    let num = (q_coarse - q_medium).abs();
    let den = (q_medium - q_fine).abs();
    if den < 1e-30 || num < 1e-30 {
        return 0.0;
    }
    (num / den).ln() / r.ln()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- LbmParameters ---

    #[test]
    fn test_params_new() {
        let p = LbmParameters::new(1.0, 64, 32);
        assert_eq!(p.omega, 1.0);
        assert_eq!(p.nx, 64);
        assert_eq!(p.ny, 32);
        assert_eq!(p.timestep, 1.0);
    }

    #[test]
    fn test_params_total_nodes() {
        let p = LbmParameters::new(1.0, 8, 4);
        assert_eq!(p.total_nodes(), 32);
    }

    #[test]
    fn test_params_viscosity_omega_1() {
        let p = LbmParameters::new(1.0, 10, 10);
        // nu = (1/3) * (1/1.0 - 0.5) = (1/3) * 0.5 = 1/6
        let expected = (1.0 / 3.0) * 0.5;
        assert!((p.viscosity() - expected).abs() < 1e-14);
    }

    #[test]
    fn test_params_viscosity_omega_2() {
        let p = LbmParameters::new(2.0, 10, 10);
        // nu = (1/3) * (0.5 - 0.5) = 0
        let expected = 0.0;
        assert!((p.viscosity() - expected).abs() < 1e-14);
    }

    #[test]
    fn test_params_viscosity_small_omega() {
        let p = LbmParameters::new(0.1, 10, 10);
        // nu = (1/3) * (10 - 0.5) = (1/3) * 9.5
        let expected = (1.0 / 3.0) * 9.5;
        assert!((p.viscosity() - expected).abs() < 1e-12);
    }

    // --- stability_analysis ---

    #[test]
    fn test_stability_in_range() {
        assert!(stability_analysis(1.0));
        assert!(stability_analysis(0.5));
        assert!(stability_analysis(1.9));
        assert!(stability_analysis(0.01));
    }

    #[test]
    fn test_stability_out_of_range() {
        assert!(!stability_analysis(0.0));
        assert!(!stability_analysis(2.0));
        assert!(!stability_analysis(-0.1));
        assert!(!stability_analysis(2.1));
    }

    #[test]
    fn test_stability_margin_at_center() {
        // omega = 1.0: margin = min(1.0, 1.0) = 1.0
        let m = stability_margin(1.0);
        assert!((m - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_stability_margin_near_boundary() {
        let m = stability_margin(0.1);
        assert!((m - 0.1).abs() < 1e-14);
        let m2 = stability_margin(1.9);
        assert!((m2 - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_stability_margin_unstable() {
        assert_eq!(stability_margin(0.0), f64::NEG_INFINITY);
        assert_eq!(stability_margin(2.0), f64::NEG_INFINITY);
    }

    // --- optimal_relaxation (TRT) ---

    #[test]
    fn test_trt_magic_constant() {
        assert!((TRT_MAGIC - 3.0 / 16.0).abs() < 1e-15);
    }

    #[test]
    fn test_optimal_relaxation_magic_3_16() {
        let omega_plus = 1.0;
        let omega_minus = optimal_relaxation(omega_plus, TRT_MAGIC);
        let tau_plus = 1.0 / omega_plus;
        let tau_minus = 1.0 / omega_minus;
        assert!((tau_plus * tau_minus - TRT_MAGIC).abs() < 1e-12);
    }

    #[test]
    fn test_optimal_relaxation_roundtrip() {
        let omega_plus = 1.5;
        let omega_minus = optimal_relaxation(omega_plus, TRT_MAGIC);
        let tau_plus = 1.0 / omega_plus;
        let tau_minus = 1.0 / omega_minus;
        assert!((tau_plus * tau_minus - TRT_MAGIC).abs() < 1e-12);
    }

    #[test]
    fn test_optimal_relaxation_custom_magic() {
        let magic = 0.25;
        let omega_plus = 1.2;
        let omega_minus = optimal_relaxation(omega_plus, magic);
        let tau_plus = 1.0 / omega_plus;
        let tau_minus = 1.0 / omega_minus;
        assert!((tau_plus * tau_minus - magic).abs() < 1e-12);
    }

    // --- convergence_criterion ---

    #[test]
    fn test_convergence_zero_residual() {
        let u = vec![0.1, 0.2, 0.3];
        let v = vec![0.0, 0.1, 0.2];
        let res = convergence_criterion(&u, &v, &u, &v);
        assert!(res < 1e-15);
    }

    #[test]
    fn test_convergence_known_residual() {
        // delta_ux = [1,0], delta_uy = [0,0]
        // sum_sq = 1, n = 2 => L2 = sqrt(0.5)
        let ux_old = vec![0.0, 0.0];
        let uy_old = vec![0.0, 0.0];
        let ux_new = vec![1.0, 0.0];
        let uy_new = vec![0.0, 0.0];
        let res = convergence_criterion(&ux_old, &uy_old, &ux_new, &uy_new);
        let expected = (0.5_f64).sqrt();
        assert!((res - expected).abs() < 1e-14);
    }

    #[test]
    fn test_convergence_empty() {
        let res = convergence_criterion(&[], &[], &[], &[]);
        assert_eq!(res, 0.0);
    }

    // --- auto_tune_omega ---

    #[test]
    fn test_auto_tune_low_re() {
        let omega = auto_tune_omega(10.0, 0.1, 100.0, 1e-6);
        assert!(stability_analysis(omega));
        // Check that resulting Re matches
        let cs2 = 1.0 / 3.0;
        let nu = cs2 * (1.0 / omega - 0.5);
        let re = 0.1 * 100.0 / nu;
        assert!((re - 10.0).abs() < 0.01);
    }

    #[test]
    fn test_auto_tune_high_re() {
        let omega = auto_tune_omega(1000.0, 0.1, 100.0, 1e-6);
        assert!(stability_analysis(omega));
        let cs2 = 1.0 / 3.0;
        let nu = cs2 * (1.0 / omega - 0.5);
        let re = 0.1 * 100.0 / nu;
        assert!((re - 1000.0).abs() < 1.0);
    }

    #[test]
    fn test_auto_tune_in_stable_range() {
        for &re in &[50.0, 200.0, 500.0] {
            let omega = auto_tune_omega(re, 0.05, 100.0, 1e-8);
            assert!(
                stability_analysis(omega),
                "omega={omega} not stable for Re={re}"
            );
        }
    }

    // --- performance_metrics ---

    #[test]
    fn test_performance_metrics_basic() {
        // 1M nodes, 100 steps, 1 second => 100 MLUPS
        let mlups = performance_metrics(1_000_000, 100, 1.0);
        assert!((mlups - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_performance_metrics_zero_time() {
        let mlups = performance_metrics(1_000_000, 100, 0.0);
        assert_eq!(mlups, 0.0);
    }

    #[test]
    fn test_performance_metrics_small() {
        // 100 nodes, 1 step, 1 second => 100/1e6 = 0.0001 MLUPS
        let mlups = performance_metrics(100, 1, 1.0);
        assert!((mlups - 1e-4).abs() < 1e-10);
    }

    // --- grid_refinement_error ---

    #[test]
    fn test_richardson_error_second_order() {
        // For 2nd order: error = (q_fine - q_coarse) / (4 - 1) = (q_fine - q_coarse) / 3
        let q_c = 1.0;
        let q_f = 1.1;
        let err = grid_refinement_error(q_c, q_f, 2.0, 2.0);
        let expected = (1.1 - 1.0) / 3.0;
        assert!((err - expected).abs() < 1e-14);
    }

    #[test]
    fn test_richardson_error_zero() {
        let err = grid_refinement_error(1.5, 1.5, 2.0, 2.0);
        assert_eq!(err, 0.0);
    }

    #[test]
    fn test_richardson_r1_zero() {
        // r=1.0 => r^p - 1 = 0 => return 0
        let err = grid_refinement_error(1.0, 2.0, 1.0, 2.0);
        assert_eq!(err, 0.0);
    }

    // --- convergence_order ---

    #[test]
    fn test_convergence_order_second_order() {
        // For exact 2nd order: q_coarse - q_medium = 4 * (q_medium - q_fine) with r=2
        // p = ln(4) / ln(2) = 2
        let q_f = 1.0;
        let q_m = 1.001;
        let q_c = 1.001 + 4.0 * 0.001; // = 1.005
        let p = convergence_order(q_c, q_m, q_f, 2.0);
        assert!((p - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_convergence_order_first_order() {
        // p=1: q_c - q_m = 2 * (q_m - q_f)
        let q_f = 1.0;
        let q_m = 1.001;
        let q_c = 1.001 + 2.0 * 0.001; // = 1.003
        let p = convergence_order(q_c, q_m, q_f, 2.0);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_convergence_order_zero_denominator() {
        // q_m == q_f => den == 0 => return 0
        let p = convergence_order(1.1, 1.0, 1.0, 2.0);
        assert_eq!(p, 0.0);
    }
}
