// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PME auto-tuning: closed-form Kolafa-Perram / Deserno-Holm error estimates and
//! optimal (α, grid) selection for a target RMS force error.
//!
//! Every estimate is an RMS force error per particle in kJ mol⁻¹ Å⁻¹ (the Coulomb
//! constant `K` = [`COULOMB_K`] is part of the estimate), so a target error in the
//! same unit can be compared with it directly.
//!
//! # Which Deserno-Holm order for which PME
//!
//! [`pme_reciprocal_rms_force_error`] is the Deserno-Holm estimate for a charge
//! assignment of order `P` whose mesh field is obtained by *spectral* (`ik`)
//! differentiation.  The PME of this crate ([`super::pme::pme_reciprocal_forces`])
//! differentiates the order-`P` B-spline analytically instead, so the mesh force is
//! interpolated with the order-`P − 1` derivative spline and its error follows the
//! Deserno-Holm estimate of order `P − 1` (measured: within 20 % over meshes of
//! 16 to 64 points, while the order-`P` estimate is 5 to 25 times too small).
//! [`PmeAutoTuner`] therefore evaluates the estimate at order `spline_order − 1`.
//!
//! # References
//!
//! - Kolafa & Perram (1992): real-space RMS force error estimate.
//! - Deserno & Holm, J. Chem. Phys. 109, 7694 (1998): reciprocal-space PME-mesh error.
//! - Essmann et al., J. Chem. Phys. 103, 8577 (1995): smooth PME formulation.

use super::coulomb::COULOMB_K;
use std::f64::consts::PI;
use thiserror::Error;

/// Highest charge-assignment order with a tabulated Deserno-Holm error series.
pub const MAX_ESTIMATED_SPLINE_ORDER: usize = 7;

/// `true` for a finite, strictly positive number.
fn is_positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

// ---------------------------------------------------------------------------
// Real-space RMS force error estimate (Kolafa-Perram 1992)
// ---------------------------------------------------------------------------

/// Kolafa-Perram (1992) real-space RMS force error estimate for PME.
///
/// Computes the per-atom RMS force error from truncating the real-space
/// Ewald sum at cutoff `rc` for a system with Ewald splitting parameter `α`.
///
/// Formula:
/// ```text
/// ΔF_real ≈ 2 K Σqᵢ² / √(N · rc · V) · exp(−α²·rc²)
/// ```
/// with `K` = [`COULOMB_K`] (kJ mol⁻¹ Å e⁻²), so the result is in kJ mol⁻¹ Å⁻¹.
///
/// Returns `NaN` for `n = 0` or a cutoff or volume that is not finite and positive.
///
/// # Arguments
/// * `q_sq_sum` — Σ qᵢ² (sum of squared charges, e²)
/// * `n`        — number of particles
/// * `rc`       — real-space cutoff (Å)
/// * `volume`   — box volume (Å³)
/// * `alpha`    — Ewald splitting parameter (Å⁻¹)
pub fn pme_real_space_rms_force_error(
    q_sq_sum: f64,
    n: usize,
    rc: f64,
    volume: f64,
    alpha: f64,
) -> f64 {
    if n == 0 || !is_positive(rc) || !is_positive(volume) {
        return f64::NAN;
    }
    let scale = (n as f64 * rc * volume).sqrt();
    2.0 * COULOMB_K * q_sq_sum / scale * (-alpha * alpha * rc * rc).exp()
}

// ---------------------------------------------------------------------------
// Reciprocal-space RMS force error estimate (Deserno-Holm 1998)
// ---------------------------------------------------------------------------

/// Coefficients `a_m^{(P)}`, `m = 0 .. P-1`, of the Deserno-Holm mesh-error series
/// for the charge-assignment order `P` (Deserno & Holm 1998, Appendix; also tabulated
/// in the LAMMPS PPPM error estimator).  `None` outside `1 ..= 7`.
fn deserno_holm_coefficients(order: usize) -> Option<&'static [f64]> {
    const P1: [f64; 1] = [2.0 / 3.0];
    const P2: [f64; 2] = [1.0 / 50.0, 5.0 / 294.0];
    const P3: [f64; 3] = [1.0 / 588.0, 7.0 / 1440.0, 21.0 / 3872.0];
    const P4: [f64; 4] = [
        1.0 / 4320.0,
        3.0 / 1936.0,
        7601.0 / 2271360.0,
        143.0 / 28800.0,
    ];
    const P5: [f64; 5] = [
        1.0 / 23232.0,
        7601.0 / 13628160.0,
        143.0 / 69120.0,
        517231.0 / 106536960.0,
        106640677.0 / 11737571328.0,
    ];
    const P6: [f64; 6] = [
        691.0 / 68140800.0,
        13.0 / 57600.0,
        47021.0 / 35512320.0,
        9694607.0 / 2095994880.0,
        733191589.0 / 59609088000.0,
        326190917.0 / 11700633600.0,
    ];
    const P7: [f64; 7] = [
        1.0 / 345600.0,
        3617.0 / 35512320.0,
        745739.0 / 838397952.0,
        56399353.0 / 12773376000.0,
        25091609.0 / 1560084480.0,
        1755948832039.0 / 36229939200000.0,
        4887769399.0 / 37838389248.0,
    ];
    match order {
        1 => Some(&P1),
        2 => Some(&P2),
        3 => Some(&P3),
        4 => Some(&P4),
        5 => Some(&P5),
        6 => Some(&P6),
        7 => Some(&P7),
        _ => None,
    }
}

/// Deserno-Holm (1998, JCP 109, 7694) PME-mesh reciprocal-space RMS force error
/// estimate.
///
/// For a mesh with spacing `h = L/M` along an axis of length `L` and a charge
/// assignment of order `P` with spectral (`ik`) differentiation of the mesh
/// potential,
/// ```text
/// ΔF_axis ≈ K Σqᵢ² · (h α)^P · √( α L √(2π) · Σ_{m<P} a_m^{(P)} (h α)^{2m} / N ) / L²
/// ```
/// with the tabulated coefficients `a_m^{(P)}` (`P = 1 ..= 7`) and `K` =
/// [`COULOMB_K`], so the result is in kJ mol⁻¹ Å⁻¹.  The three axes are combined as
/// the RMS `√((ΔF_x² + ΔF_y² + ΔF_z²) / 3)`.  The estimate grows monotonically with
/// the mesh spacing `h`.
///
/// Returns `NaN` for a spline order outside `1 ..= 7`, `n = 0`, an empty grid
/// dimension or a box length that is not finite and positive.
///
/// The series was checked against an independent P3M implementation with the optimal
/// influence function for `P = 1 ..= 7`.  The analytic-gradient PME of this crate
/// corresponds to the order `P − 1` of this estimate, see
/// [`PmeAutoTuner::reciprocal_error`].
///
/// # Arguments
/// * `q_sq_sum`    — Σ qᵢ² (e²)
/// * `n`           — number of particles
/// * `alpha`       — Ewald splitting parameter (Å⁻¹)
/// * `grid`        — mesh dimensions [mx, my, mz]
/// * `box_lengths` — box side lengths [Lx, Ly, Lz] (Å)
/// * `spline_order`— B-spline order `P` (typically 4)
pub fn pme_reciprocal_rms_force_error(
    q_sq_sum: f64,
    n: usize,
    alpha: f64,
    grid: [usize; 3],
    box_lengths: [f64; 3],
    spline_order: usize,
) -> f64 {
    let Some(coefficients) = deserno_holm_coefficients(spline_order) else {
        return f64::NAN;
    };
    if n == 0
        || grid.contains(&0)
        || !alpha.is_finite()
        || box_lengths.iter().any(|length| !is_positive(*length))
    {
        return f64::NAN;
    }

    let order = spline_order as i32;
    let mut sum_of_squares = 0.0_f64;
    for (&length, &cells) in box_lengths.iter().zip(grid.iter()) {
        let h_alpha = length / cells as f64 * alpha;
        // Σ_m a_m x^m with x = (h α)², evaluated by Horner's rule.
        let series = coefficients
            .iter()
            .rev()
            .fold(0.0_f64, |acc, a| acc * h_alpha * h_alpha + a);
        let mesh_factor =
            h_alpha.powi(order) * (alpha * length * (2.0 * PI).sqrt() * series / n as f64).sqrt();
        let axis_error = COULOMB_K * q_sq_sum * mesh_factor / (length * length);
        sum_of_squares += axis_error * axis_error;
    }
    (sum_of_squares / 3.0).sqrt()
}

// ---------------------------------------------------------------------------
// Error types
// ---------------------------------------------------------------------------

/// Errors that can occur during PME parameter auto-tuning.
#[derive(Debug, Error)]
pub enum PmeTuningError {
    /// Real-space error is irreducible — cannot find alpha achieving the target.
    #[error("Failed to find valid alpha for target error {0}: real-space error irreducible")]
    AlphaNotFound(f64),
    /// Reciprocal error requires a grid larger than the configured maximum.
    #[error("Grid too large: reciprocal error not achievable within max grid {0}")]
    GridTooLarge(usize),
    /// An input parameter was invalid.
    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),
}

// ---------------------------------------------------------------------------
// PmeParams — tuned parameter set
// ---------------------------------------------------------------------------

/// A fully-specified set of PME parameters (α, grid, r_cut).
///
/// Can be constructed manually or via [`PmeParams::auto_tune`].
#[derive(Debug, Clone)]
pub struct PmeParams {
    /// Ewald splitting parameter α (Å⁻¹).
    pub alpha: f64,
    /// PME grid dimensions [mx, my, mz].
    pub grid: [usize; 3],
    /// Real-space cutoff (Å).
    pub r_cut: f64,
}

impl PmeParams {
    /// Construct PME parameters manually.
    pub fn new(alpha: f64, grid: [usize; 3], r_cut: f64) -> Self {
        Self { alpha, grid, r_cut }
    }

    /// Auto-select (α, grid) to achieve `target_rms_force_error` with minimal grid cost.
    ///
    /// # Arguments
    /// * `charges`                — per-atom charge array (e), length N
    /// * `box_lengths`            — orthorhombic box sides [Lx, Ly, Lz] (Å)
    /// * `target_rms_force_error` — desired RMS force error (kJ mol⁻¹ Å⁻¹)
    /// * `r_cut`                  — real-space cutoff (Å)
    ///
    /// # Errors
    /// Returns [`PmeTuningError`] if no valid parameter set is found.
    pub fn auto_tune(
        charges: &[f64],
        box_lengths: [f64; 3],
        target_rms_force_error: f64,
        r_cut: f64,
    ) -> Result<Self, PmeTuningError> {
        PmeAutoTuner::new(target_rms_force_error, r_cut)?.tune(charges, box_lengths)
    }
}

// ---------------------------------------------------------------------------
// PmeAutoTuner
// ---------------------------------------------------------------------------

/// Auto-tuner that finds optimal PME parameters (α, grid) for a given target
/// RMS force error, using the Kolafa-Perram + Deserno-Holm error estimates.
///
/// The error budget is split equally between real-space and reciprocal-space:
/// `target_real = target_recip = target / √2`.
///
/// Alpha is found analytically from the KP real-space estimate: it is the smallest
/// splitting parameter whose real-space error is within budget, which also gives the
/// cheapest mesh.  The grid is the smallest FFT-friendly grid (sizes with the prime
/// factors 2, 3, 5 only) whose reciprocal error estimate falls below budget; the
/// estimate is the Deserno-Holm formula evaluated for the analytic-gradient PME of
/// this crate (see [`PmeAutoTuner::reciprocal_error`]).
pub struct PmeAutoTuner {
    /// Target RMS force error (kJ mol⁻¹ Å⁻¹).
    pub target_rms_force_error: f64,
    /// Real-space cutoff (Å).
    pub r_cut: f64,
    /// B-spline interpolation order of the PME (default 4, supported `2 ..= 8`).
    pub spline_order: usize,
}

impl PmeAutoTuner {
    /// Create a new auto-tuner.
    ///
    /// # Arguments
    /// * `target_rms_force_error` — desired RMS force error (must be > 0)
    /// * `r_cut`                  — real-space cutoff (must be > 0)
    pub fn new(target_rms_force_error: f64, r_cut: f64) -> Result<Self, PmeTuningError> {
        if !is_positive(target_rms_force_error) {
            return Err(PmeTuningError::InvalidParameter(
                "target error must be positive".to_owned(),
            ));
        }
        if !is_positive(r_cut) {
            return Err(PmeTuningError::InvalidParameter(
                "r_cut must be positive".to_owned(),
            ));
        }
        Ok(Self {
            target_rms_force_error,
            r_cut,
            spline_order: 4,
        })
    }

    /// Estimated RMS force error (kJ mol⁻¹ Å⁻¹) of the reciprocal-space part of the
    /// PME of this crate for the tuner's B-spline order.
    ///
    /// The PME forces are obtained by differentiating the order-`P` B-spline of the
    /// charge assignment analytically, so the mesh force is interpolated with the
    /// order-`P − 1` derivative spline and the error follows the Deserno-Holm estimate
    /// ([`pme_reciprocal_rms_force_error`]) of order `P − 1`.  Against exact Ewald
    /// forces this predicts the measured error of `pme_reciprocal_forces` within
    /// about 20 % for cubic splines (`P = 4`), where the order-`P` estimate would be
    /// 5 to 25 times too small.
    ///
    /// Returns `NaN` for a `spline_order` outside `2 ..= 8` and for invalid arguments.
    pub fn reciprocal_error(
        &self,
        q_sq_sum: f64,
        n: usize,
        alpha: f64,
        grid: [usize; 3],
        box_lengths: [f64; 3],
    ) -> f64 {
        pme_reciprocal_rms_force_error(
            q_sq_sum,
            n,
            alpha,
            grid,
            box_lengths,
            self.spline_order.saturating_sub(1),
        )
    }

    /// Find optimal (α, grid) for the given charge array and box.
    ///
    /// # Arguments
    /// * `charges`     — per-atom charge array (e), length N
    /// * `box_lengths` — orthorhombic box sides [Lx, Ly, Lz] (Å)
    ///
    /// # Errors
    /// * [`PmeTuningError::InvalidParameter`] for an empty or non-finite charge array,
    ///   a box edge that is not finite and positive, a non-positive target or cutoff, or
    ///   a `spline_order` outside `2 ..= 8`.
    /// * [`PmeTuningError::AlphaNotFound`] if no splitting parameter meets the budget.
    /// * [`PmeTuningError::GridTooLarge`] if the target needs more than 512 cells along
    ///   the longest box edge.
    pub fn tune(
        &self,
        charges: &[f64],
        box_lengths: [f64; 3],
    ) -> Result<PmeParams, PmeTuningError> {
        let supported = 2..=MAX_ESTIMATED_SPLINE_ORDER + 1;
        if !supported.contains(&self.spline_order) {
            return Err(PmeTuningError::InvalidParameter(format!(
                "spline order {} is not supported ({}..={})",
                self.spline_order,
                supported.start(),
                supported.end()
            )));
        }
        if !is_positive(self.target_rms_force_error) || !is_positive(self.r_cut) {
            return Err(PmeTuningError::InvalidParameter(
                "target error and r_cut must be positive".to_owned(),
            ));
        }
        let n = charges.len();
        if n == 0 {
            return Err(PmeTuningError::InvalidParameter("no particles".to_owned()));
        }
        if charges.iter().any(|q| !q.is_finite()) {
            return Err(PmeTuningError::InvalidParameter(
                "charges must be finite".to_owned(),
            ));
        }
        if box_lengths.iter().any(|length| !is_positive(*length)) {
            return Err(PmeTuningError::InvalidParameter(
                "box edge lengths must be positive".to_owned(),
            ));
        }
        let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
        let q_sq_sum: f64 = charges.iter().map(|q| q * q).sum();

        // Equally split error budget between real and reciprocal contributions.
        // By Pythagoras (uncorrelated): target² ≥ target_real² + target_recip².
        // Optimal equal split: target_real = target_recip = target / √2.
        let target_each = self.target_rms_force_error / 2_f64.sqrt();

        // ---- Find alpha from real-space KP estimate ----
        // ΔF_real = 2·K·Σq²/√(N·rc·V)·exp(-α²·rc²) = target_each
        // => exp(-α²·rc²) = target_each / prefactor
        // => α = √(-ln(ratio)) / rc
        let prefactor = pme_real_space_rms_force_error(q_sq_sum, n, self.r_cut, volume, 0.0);

        let alpha = if q_sq_sum < 1e-30 {
            // No charge, use a small default alpha
            0.3 / self.r_cut
        } else {
            let ratio = target_each / prefactor;
            if ratio >= 1.0 {
                // Real-space error already below budget at α→0; use small alpha
                0.1 / self.r_cut
            } else if ratio <= 0.0 || !ratio.is_finite() {
                return Err(PmeTuningError::AlphaNotFound(self.target_rms_force_error));
            } else {
                let ln_ratio = ratio.ln(); // ln_ratio < 0 since ratio < 1
                (-ln_ratio / (self.r_cut * self.r_cut)).sqrt()
            }
        };

        // ---- Find minimal grid such that reciprocal error ≤ target_each ----
        // The grid of "resolution" k scales every dimension by the box aspect ratio and
        // rounds up to a "good" FFT size (only prime factors 2, 3, 5).
        let k_min = 4_usize;
        let k_max_search = 512_usize;
        let l_max = box_lengths[0].max(box_lengths[1]).max(box_lengths[2]);
        let grid_for = |k: usize| -> [usize; 3] {
            let size =
                |length: f64| next_good_grid_size((k as f64 * length / l_max).ceil() as usize);
            [
                size(box_lengths[0]),
                size(box_lengths[1]),
                size(box_lengths[2]),
            ]
        };
        let within_budget = |k: usize| {
            self.reciprocal_error(q_sq_sum, n, alpha, grid_for(k), box_lengths) <= target_each
        };

        // Coarse search: increase k by ~25% (at least 2) until the budget is met.
        let mut k_fail = k_min - 1;
        let mut k_pass = None;
        let mut k = k_min;
        while k <= k_max_search {
            if within_budget(k) {
                k_pass = Some(k);
                break;
            }
            k_fail = k;
            k = (k * 5 / 4).max(k + 2);
        }
        let k_pass = k_pass.ok_or(PmeTuningError::GridTooLarge(k_max_search))?;

        // Refine: the estimated error never increases with k, so bisect for the
        // smallest k in (k_fail, k_pass] that still meets the budget.
        let (mut below, mut above) = (k_fail, k_pass);
        while above - below > 1 {
            let middle = below + (above - below) / 2;
            if within_budget(middle) {
                above = middle;
            } else {
                below = middle;
            }
        }

        Ok(PmeParams {
            alpha,
            grid: grid_for(above),
            r_cut: self.r_cut,
        })
    }
}

// ---------------------------------------------------------------------------
// Grid-size utilities
// ---------------------------------------------------------------------------

/// Return the smallest integer ≥ `n` (and ≥ 4, even) whose prime factors are
/// only 2, 3, and 5 — a "smooth" or "regular" FFT size.
pub fn next_good_grid_size(n: usize) -> usize {
    // Enforce minimum of 4 and round up to at least the next even number.
    let start = n.max(4);
    let mut candidate = if start.is_multiple_of(2) {
        start
    } else {
        start + 1
    };
    loop {
        if is_good_size(candidate) {
            return candidate;
        }
        candidate += 2; // stay even
    }
}

/// Return `true` if all prime factors of `n` are in {2, 3, 5}.
fn is_good_size(mut n: usize) -> bool {
    if n == 0 {
        return false;
    }
    for p in [2_usize, 3, 5] {
        while n.is_multiple_of(p) {
            n /= p;
        }
    }
    n == 1
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electrostatics::pme::pme_reciprocal_forces;
    use oxiphysics_core::random::Xoshiro256;

    // Helper: build alternating ±1 charge array of length n
    fn alternating_charges(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 })
            .collect()
    }

    /// SPC/E water charges (O -0.8476 e, H +0.4238 e) of `n_waters` molecules.
    fn spce_charges(n_waters: usize) -> Vec<f64> {
        (0..n_waters)
            .flat_map(|_| [-0.8476_f64, 0.4238, 0.4238])
            .collect()
    }

    /// 1000 SPC/E waters in a 30 Å cube: N = 3000 atoms, Σq² = 1077.64 e².
    const WATER_BOX: [f64; 3] = [30.0; 3];
    const WATER_VOLUME: f64 = 27_000.0;
    const WATER_ATOMS: usize = 3000;

    fn water_q_sq_sum() -> f64 {
        spce_charges(1000).iter().map(|q| q * q).sum()
    }

    // -----------------------------------------------------------------------
    // Test 1: next_good_grid_size correctness
    // -----------------------------------------------------------------------
    #[test]
    fn test_next_good_grid_size_correctness() {
        assert_eq!(next_good_grid_size(4), 4);
        assert_eq!(next_good_grid_size(5), 6);
        assert_eq!(next_good_grid_size(7), 8);
        assert_eq!(next_good_grid_size(10), 10);
        assert_eq!(next_good_grid_size(11), 12);
        assert_eq!(next_good_grid_size(13), 16);
        // Additional spot-checks
        assert_eq!(next_good_grid_size(1), 4); // minimum is 4
        assert_eq!(next_good_grid_size(4), 4);
        assert_eq!(next_good_grid_size(24), 24); // 24 = 2³·3 → good
        assert_eq!(next_good_grid_size(25), 30); // 25 = 5²; but 26=2·13 not good, 28=4·7 not good, 30=2·3·5 good
        assert_eq!(next_good_grid_size(32), 32); // 32 = 2⁵ → good
    }

    // -----------------------------------------------------------------------
    // Test 2: real-space error is the Kolafa-Perram estimate in kJ mol⁻¹ Å⁻¹
    // -----------------------------------------------------------------------
    #[test]
    fn test_real_space_error_matches_kolafa_perram_in_kj_per_mol_per_angstrom() {
        let q2 = water_q_sq_sum();
        let (rc, alpha) = (10.0_f64, 0.31_f64);
        let err = pme_real_space_rms_force_error(q2, WATER_ATOMS, rc, WATER_VOLUME, alpha);
        let expected = 2.0 * COULOMB_K * q2 / (3000.0 * 10.0 * 27_000.0_f64).sqrt()
            * (-(alpha * rc).powi(2)).exp();
        assert!(
            ((err - expected) / expected).abs() < 1e-12,
            "err = {err}, expected {expected}"
        );
        // Order of magnitude: ~7e-3 kJ/mol/Å at α = 0.31 Å⁻¹, rc = 10 Å.
        assert!(err > 6.0e-3 && err < 8.0e-3, "err = {err}");
        // Without any screening the truncation error is the prefactor, 105 kJ/mol/Å
        // (the old √Σq² / no-K formula returned ~0.2 here).
        let unscreened = pme_real_space_rms_force_error(q2, WATER_ATOMS, rc, WATER_VOLUME, 0.0);
        assert!((unscreened - 105.214).abs() < 0.01, "{unscreened}");
    }

    #[test]
    fn test_real_space_error_decreases_with_alpha_and_cutoff() {
        let q2 = water_q_sq_sum();
        let mut previous = f64::INFINITY;
        for step in 1..=20 {
            let alpha = 0.05 * f64::from(step);
            let err = pme_real_space_rms_force_error(q2, WATER_ATOMS, 10.0, WATER_VOLUME, alpha);
            assert!(err < previous, "alpha = {alpha}: {err} !< {previous}");
            previous = err;
        }
        let short = pme_real_space_rms_force_error(q2, WATER_ATOMS, 8.0, WATER_VOLUME, 0.3);
        let long = pme_real_space_rms_force_error(q2, WATER_ATOMS, 12.0, WATER_VOLUME, 0.3);
        assert!(long < short);
    }

    #[test]
    fn test_error_estimates_reject_degenerate_inputs() {
        assert!(pme_real_space_rms_force_error(10.0, 0, 8.0, 1000.0, 0.3).is_nan());
        assert!(pme_real_space_rms_force_error(10.0, 10, 0.0, 1000.0, 0.3).is_nan());
        assert!(pme_real_space_rms_force_error(10.0, 10, 8.0, -1.0, 0.3).is_nan());
        let box_lengths = [10.0; 3];
        assert!(pme_reciprocal_rms_force_error(10.0, 10, 0.3, [16; 3], box_lengths, 0).is_nan());
        assert!(pme_reciprocal_rms_force_error(10.0, 10, 0.3, [16; 3], box_lengths, 8).is_nan());
        assert!(pme_reciprocal_rms_force_error(10.0, 0, 0.3, [16; 3], box_lengths, 4).is_nan());
        assert!(
            pme_reciprocal_rms_force_error(10.0, 10, 0.3, [0, 16, 16], box_lengths, 4).is_nan()
        );
        assert!(
            pme_reciprocal_rms_force_error(10.0, 10, 0.3, [16; 3], [10.0, -10.0, 10.0], 4).is_nan()
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: reciprocal error follows the Deserno-Holm formula
    // -----------------------------------------------------------------------
    #[test]
    fn test_reciprocal_error_matches_the_deserno_holm_formula() {
        // Direct evaluation for the cubic water box, h = L / M, with the p = 4 and p = 6
        // coefficient series written out with their literal fractions.
        let q2 = water_q_sq_sum();
        let alpha = 0.31_f64;
        let l = 30.0_f64;
        let reference = |m: usize, order: i32, coefficients: &[f64]| -> f64 {
            let h = l / m as f64;
            let x = h * alpha;
            let sum: f64 = coefficients
                .iter()
                .enumerate()
                .map(|(i, a)| a * x.powi(2 * i as i32))
                .sum();
            COULOMB_K * q2 * x.powi(order) * (alpha * l * (2.0 * PI).sqrt() * sum / 3000.0).sqrt()
                / (l * l)
        };
        let p4 = [
            1.0 / 4320.0,
            3.0 / 1936.0,
            7601.0 / 2271360.0,
            143.0 / 28800.0,
        ];
        let p6 = [
            691.0 / 68140800.0,
            13.0 / 57600.0,
            47021.0 / 35512320.0,
            9694607.0 / 2095994880.0,
            733191589.0 / 59609088000.0,
            326190917.0 / 11700633600.0,
        ];
        for m in [24_usize, 32, 48, 64] {
            let e4 = pme_reciprocal_rms_force_error(q2, 3000, alpha, [m; 3], WATER_BOX, 4);
            let e6 = pme_reciprocal_rms_force_error(q2, 3000, alpha, [m; 3], WATER_BOX, 6);
            assert!(
                ((e4 - reference(m, 4, &p4)) / e4).abs() < 1e-12,
                "p = 4, M = {m}"
            );
            assert!(
                ((e6 - reference(m, 6, &p6)) / e6).abs() < 1e-12,
                "p = 6, M = {m}"
            );
        }
        // Frozen values of the estimate at α = 0.31 Å⁻¹ (kJ mol⁻¹ Å⁻¹).
        let e32 = pme_reciprocal_rms_force_error(q2, 3000, alpha, [32; 3], WATER_BOX, 4);
        let e48 = pme_reciprocal_rms_force_error(q2, 3000, alpha, [48; 3], WATER_BOX, 4);
        assert!((e32 - 0.020_639_743).abs() < 1e-8, "M = 32: {e32}");
        assert!((e48 - 0.003_547_111).abs() < 1e-8, "M = 48: {e48}");
    }

    #[test]
    fn test_reciprocal_error_is_the_rms_of_the_three_axes() {
        // An anisotropic mesh: the estimate is sqrt((fx² + fy² + fz²) / 3) with each axis
        // evaluated for its own box edge and spacing.
        let q2 = 40.0;
        let (n, alpha) = (200, 0.35);
        let box_lengths = [20.0, 30.0, 40.0];
        let grid = [24, 32, 48];
        let combined = pme_reciprocal_rms_force_error(q2, n, alpha, grid, box_lengths, 4);
        let mut sum_sq = 0.0;
        for axis in 0..3 {
            // A cube whose edge and cell count equal this axis' reproduces its error.
            let single = pme_reciprocal_rms_force_error(
                q2,
                n,
                alpha,
                [grid[axis]; 3],
                [box_lengths[axis]; 3],
                4,
            );
            sum_sq += single * single;
        }
        let expected = (sum_sq / 3.0).sqrt();
        assert!(((combined - expected) / expected).abs() < 1e-12);
    }

    #[test]
    fn test_reciprocal_error_increases_monotonically_with_the_mesh_spacing() {
        let q2 = water_q_sq_sum();
        for order in 2..=7 {
            let mut previous = 0.0_f64;
            for cells in [96_usize, 64, 48, 40, 32, 24, 16, 12, 8] {
                let err =
                    pme_reciprocal_rms_force_error(q2, 3000, 0.31, [cells; 3], WATER_BOX, order);
                assert!(
                    err > previous && err.is_finite(),
                    "order {order}, M = {cells}: {err} !> {previous}"
                );
                previous = err;
            }
        }
    }

    #[test]
    fn test_reciprocal_error_decreases_with_the_spline_order_for_fine_meshes() {
        let q2 = water_q_sq_sum();
        let mut previous = f64::INFINITY;
        for order in 3..=7 {
            let err = pme_reciprocal_rms_force_error(q2, 3000, 0.31, [32; 3], WATER_BOX, order);
            assert!(err < previous, "order {order}: {err} !< {previous}");
            previous = err;
        }
    }

    // -----------------------------------------------------------------------
    // Test 4: the tuner on 1000 SPC/E waters (30 Å box, rc = 10 Å, 1e-2 kJ/mol/Å)
    // -----------------------------------------------------------------------
    #[test]
    fn test_auto_tune_one_thousand_spce_waters() {
        let charges = spce_charges(1000);
        let (target, r_cut) = (1.0e-2, 10.0);
        let tuner = PmeAutoTuner::new(target, r_cut).expect("valid tuner");
        let params = tuner
            .tune(&charges, WATER_BOX)
            .expect("auto_tune should succeed");
        // The old estimators picked α = 0.01 Å⁻¹ and a 4³ grid for this system.
        assert!(
            params.alpha > 0.25 && params.alpha < 0.35,
            "alpha = {}",
            params.alpha
        );
        assert!(
            params.grid.iter().product::<usize>() > 4 * 4 * 4,
            "grid {:?}",
            params.grid
        );
        assert!(
            params.grid.iter().all(|&cells| cells >= 24),
            "{:?}",
            params.grid
        );
        // α = 0.3100 Å⁻¹ and the smallest FFT-friendly mesh within budget (80³ misses it).
        assert!(
            (params.alpha - 0.309_963_556).abs() < 1e-8,
            "{}",
            params.alpha
        );
        assert_eq!(params.grid, [90, 90, 90]);
        assert_eq!(params.r_cut, r_cut);

        let q2 = water_q_sq_sum();
        let target_each = target / 2_f64.sqrt();
        let err_real =
            pme_real_space_rms_force_error(q2, WATER_ATOMS, r_cut, WATER_VOLUME, params.alpha);
        let err_recip =
            tuner.reciprocal_error(q2, WATER_ATOMS, params.alpha, params.grid, WATER_BOX);
        assert!(err_real <= target_each * (1.0 + 1e-9), "real {err_real}");
        assert!(err_recip <= target_each, "reciprocal {err_recip}");
        assert!((err_real * err_real + err_recip * err_recip).sqrt() <= target);

        // The grid is minimal among FFT-friendly sizes: the next smaller one misses the
        // budget (82, 84, 86, 88 all have a prime factor above 5).
        let smaller = tuner.reciprocal_error(q2, WATER_ATOMS, params.alpha, [80; 3], WATER_BOX);
        assert!(smaller > target_each, "an 80³ mesh gives {smaller}");
    }

    // -----------------------------------------------------------------------
    // Test 5: Auto-tune produces params achieving target
    // -----------------------------------------------------------------------
    #[test]
    fn test_auto_tune_achieves_target() {
        // Small ionic system: N=32, alternating ±1, cubic box L=10 Å
        let n = 32_usize;
        let box_len = 10.0_f64;
        let box_lengths = [box_len; 3];
        let charges = alternating_charges(n);
        let target = 1e-4_f64;
        let r_cut = 5.0_f64;

        let tuner = PmeAutoTuner::new(target, r_cut).expect("valid tuner");
        let params = tuner
            .tune(&charges, box_lengths)
            .expect("auto_tune should succeed");
        let via_params = PmeParams::auto_tune(&charges, box_lengths, target, r_cut)
            .expect("auto_tune should succeed");
        assert_eq!(via_params.grid, params.grid);

        assert!(
            params.alpha > 0.0,
            "alpha must be positive, got {}",
            params.alpha
        );
        assert!(
            params.grid.iter().all(|&cells| cells >= 4),
            "{:?}",
            params.grid
        );

        // Verify the tuner's own error estimates are within target
        // (this is what the tuner guarantees; the comparison with actual PME forces is
        //  done in `test_reciprocal_estimate_tracks_the_measured_pme_error`).
        let q_sq_sum: f64 = charges.iter().map(|q| q * q).sum();
        let volume = box_len * box_len * box_len;

        let err_real = pme_real_space_rms_force_error(q_sq_sum, n, r_cut, volume, params.alpha);
        let err_recip = tuner.reciprocal_error(q_sq_sum, n, params.alpha, params.grid, box_lengths);
        // Each component must be ≤ target / √2
        let target_each = target / 2_f64.sqrt();
        assert!(
            err_real <= target_each * 1.01, // 1% tolerance for floating-point rounding
            "real-space error estimate {err_real:.3e} exceeds budget {target_each:.3e}"
        );
        assert!(
            err_recip <= target_each * 1.01,
            "reciprocal error estimate {err_recip:.3e} exceeds budget {target_each:.3e}"
        );

        // Combined (Pythagorean) error should be ≤ target
        let total_estimated = (err_real * err_real + err_recip * err_recip).sqrt();
        assert!(
            total_estimated <= target * 1.05, // 5% tolerance
            "total estimated error {total_estimated:.3e} exceeds target {target:.3e}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 6: tighter targets cost more, looser targets cost less
    // -----------------------------------------------------------------------
    #[test]
    fn test_tuned_parameters_follow_the_target_accuracy() {
        let charges = spce_charges(1000);
        let tune = |target: f64| {
            PmeParams::auto_tune(&charges, WATER_BOX, target, 10.0).expect("auto_tune failed")
        };
        let loose = tune(1.0e-1);
        let medium = tune(1.0e-2);
        let tight = tune(1.0e-3);
        let cells = |p: &PmeParams| p.grid.iter().product::<usize>();
        assert!(loose.alpha < medium.alpha && medium.alpha < tight.alpha);
        assert!(
            cells(&loose) < cells(&medium),
            "{:?} vs {:?}",
            loose.grid,
            medium.grid
        );
        assert!(
            cells(&medium) < cells(&tight),
            "{:?} vs {:?}",
            medium.grid,
            tight.grid
        );
    }

    // -----------------------------------------------------------------------
    // Test 7: Error handling — invalid inputs
    // -----------------------------------------------------------------------
    #[test]
    fn test_auto_tune_error_handling() {
        // Zero target error
        assert!(
            PmeAutoTuner::new(0.0, 5.0).is_err(),
            "zero target error should fail"
        );
        assert!(PmeAutoTuner::new(f64::NAN, 5.0).is_err());
        // Negative r_cut
        assert!(
            PmeAutoTuner::new(1e-4, -1.0).is_err(),
            "negative r_cut should fail"
        );
        // Empty charge array
        let tuner = PmeAutoTuner::new(1e-4, 5.0).expect("valid tuner");
        assert!(
            tuner.tune(&[], [10.0, 10.0, 10.0]).is_err(),
            "empty charges should fail"
        );
        // Degenerate boxes and non-finite charges
        let charges = alternating_charges(8);
        assert!(tuner.tune(&charges, [10.0, 0.0, 10.0]).is_err());
        assert!(tuner.tune(&charges, [-10.0, -10.0, 10.0]).is_err());
        let mut bad = alternating_charges(8);
        bad[3] = f64::NAN;
        assert!(tuner.tune(&bad, [10.0; 3]).is_err());
    }

    #[test]
    fn test_tune_validates_the_spline_order() {
        let charges = alternating_charges(8);
        for order in [0_usize, 1, 9, 12] {
            let mut tuner = PmeAutoTuner::new(1e-1, 4.0).expect("valid tuner");
            tuner.spline_order = order;
            assert!(
                matches!(
                    tuner.tune(&charges, [10.0; 3]),
                    Err(PmeTuningError::InvalidParameter(_))
                ),
                "spline order {order} must be rejected"
            );
            assert!(
                tuner
                    .reciprocal_error(8.0, 8, 0.3, [16; 3], [10.0; 3])
                    .is_nan()
            );
        }
        // Linear and quadratic splines (effective Deserno-Holm orders 1 and 2) converge too
        // slowly for the 512-cell limit at this accuracy; every higher order tunes, and a
        // higher order needs no finer mesh.
        let charges = spce_charges(1000);
        let mut tuner = PmeAutoTuner::new(1e-2, 10.0).expect("valid tuner");
        for order in [2_usize, 3] {
            tuner.spline_order = order;
            assert!(
                matches!(
                    tuner.tune(&charges, WATER_BOX),
                    Err(PmeTuningError::GridTooLarge(512))
                ),
                "spline order {order}"
            );
        }
        let mut previous_cells = usize::MAX;
        for order in 4..=8 {
            tuner.spline_order = order;
            let params = tuner.tune(&charges, WATER_BOX).expect("tunes");
            let cells = params.grid.iter().product::<usize>();
            assert!(
                cells <= previous_cells,
                "order {order}: {cells} > {previous_cells}"
            );
            previous_cells = cells;
        }
    }

    #[test]
    fn test_uncharged_system_gets_the_minimal_grid() {
        let tuner = PmeAutoTuner::new(1e-2, 8.0).expect("valid tuner");
        let params = tuner.tune(&[0.0; 5], [12.0; 3]).expect("tunes");
        assert_eq!(params.grid, [4, 4, 4]);
        assert!((params.alpha - 0.3 / 8.0).abs() < 1e-15);
    }

    // -----------------------------------------------------------------------
    // Test 8: is_good_size / next_good_grid_size edge cases
    // -----------------------------------------------------------------------
    #[test]
    fn test_good_size_properties() {
        // All results of next_good_grid_size must be even and have only factors 2,3,5
        for n in [1, 3, 6, 7, 11, 17, 23, 31, 37, 41, 47, 100, 127] {
            let g = next_good_grid_size(n);
            assert!(g >= 4, "result {g} should be >= 4");
            assert!(g.is_multiple_of(2), "result {g} should be even");
            assert!(
                is_good_size(g),
                "result {g} (from n={n}) should be a smooth number"
            );
            assert!(g >= n, "result {g} should be >= input {n}");
        }
    }

    // -----------------------------------------------------------------------
    // Test 9: the reciprocal estimate tracks the measured PME force error
    // -----------------------------------------------------------------------
    #[test]
    fn test_reciprocal_estimate_tracks_the_measured_pme_error() {
        // 100 random neutral ±1 charges in a 20 Å cube. The mesh error of the reciprocal
        // forces at grid M is measured against a much finer mesh of the same α.  The PME
        // of this crate differentiates the cubic B-spline analytically, so its error is
        // that of the Deserno-Holm estimate one order lower (the tuner's estimate); the
        // plain order-4 estimate, valid for ik differentiation, is several times too small.
        let n = 100;
        let length = 20.0_f64;
        let box_lengths = [length; 3];
        let alpha = 0.35_f64;
        let mut rng = Xoshiro256::new(20_241_019);
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                [
                    length * rng.next_f64(),
                    length * rng.next_f64(),
                    length * rng.next_f64(),
                ]
            })
            .collect();
        let charges = alternating_charges(n);
        let q2: f64 = charges.iter().map(|q| q * q).sum();
        let tuner = PmeAutoTuner::new(1e-2, 8.0).expect("valid tuner");

        let reference = pme_reciprocal_forces(&positions, &charges, box_lengths, alpha, [72; 3]);
        let mut previous_measured = f64::INFINITY;
        for cells in [24_usize, 32, 48] {
            let forces =
                pme_reciprocal_forces(&positions, &charges, box_lengths, alpha, [cells; 3]);
            let mean_sq_error = forces
                .iter()
                .zip(reference.iter())
                .map(|(f, r)| (f[0] - r[0]).powi(2) + (f[1] - r[1]).powi(2) + (f[2] - r[2]).powi(2))
                .sum::<f64>()
                / n as f64;
            let measured = mean_sq_error.sqrt();
            let estimated = tuner.reciprocal_error(q2, n, alpha, [cells; 3], box_lengths);
            let ratio = estimated / measured;
            assert!(
                measured < previous_measured,
                "M = {cells}: measured error {measured} did not shrink"
            );
            previous_measured = measured;
            assert!(
                (0.5..2.0).contains(&ratio),
                "M = {cells}: estimate {estimated:.3e} vs measured {measured:.3e} (ratio {ratio:.2})"
            );
            let plain = pme_reciprocal_rms_force_error(q2, n, alpha, [cells; 3], box_lengths, 4);
            assert!(
                plain < measured / 3.0,
                "M = {cells}: the ik estimate {plain:.3e} should underestimate {measured:.3e}"
            );
        }
    }
}
