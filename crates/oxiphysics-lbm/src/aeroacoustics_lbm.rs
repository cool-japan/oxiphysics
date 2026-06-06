// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM aeroacoustics module.
//!
//! Implements acoustic extraction from LBM simulations including:
//! - Sound pressure level (SPL) computation
//! - Far-field pressure projection with 1/r decay
//! - Ffowcs Williams-Hawkings (FW-H) monopole analogy
//! - Strouhal number and Aeolian tone frequency
//! - Curle's analogy approximation for vortex sound power
//! - Ring-buffer pressure history tracking at monitor points

use std::f64::consts::PI;

// ============================================================================
// Parameter structs
// ============================================================================

/// Parameters for LBM aeroacoustics simulations.
#[derive(Debug, Clone)]
pub struct AeroacousticParams {
    /// Reference Mach number M = U_ref / c_s.
    pub mach_ref: f64,
    /// Lattice speed of sound c_s (typically 1/√3 in LBM units).
    pub c_s: f64,
    /// Kinematic viscosity ν (lattice units).
    pub viscosity: f64,
    /// Indices \[ix, iy\] of acoustic monitor points.
    pub monitor_points: Vec<[usize; 2]>,
}

impl AeroacousticParams {
    /// Create a new `AeroacousticParams` with sensible LBM defaults.
    pub fn new(mach_ref: f64, viscosity: f64, monitor_points: Vec<[usize; 2]>) -> Self {
        Self {
            mach_ref,
            c_s: 1.0 / 3.0_f64.sqrt(),
            viscosity,
            monitor_points,
        }
    }
}

// ============================================================================
// Aeroacoustic LBM state
// ============================================================================

/// LBM aeroacoustics simulation state.
///
/// Stores distribution functions and a pressure history at user-defined
/// monitor points for post-processing acoustic quantities.
#[derive(Debug, Clone)]
pub struct AeroacousticLBM {
    /// Number of lattice nodes in x direction.
    pub nx: usize,
    /// Number of lattice nodes in y direction.
    pub ny: usize,
    /// D2Q9 distribution functions, length nx * ny * 9.
    pub f: Vec<f64>,
    /// Per-monitor-point pressure history: `pressure_history[m][t]`.
    pub pressure_history: Vec<Vec<f64>>,
}

impl AeroacousticLBM {
    /// Create a new uniform-rest `AeroacousticLBM` grid.
    ///
    /// All distributions are initialised to the D2Q9 equilibrium at rest
    /// (rho=1, u=0), so w_0=4/9 and the 4 axis weights=1/9 and
    /// 4 diagonal weights=1/36.
    pub fn new(nx: usize, ny: usize, n_monitor: usize) -> Self {
        let d2q9_weights: [f64; 9] = [
            4.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
        ];
        let mut f = vec![0.0_f64; nx * ny * 9];
        for node in 0..(nx * ny) {
            for (alpha, &w) in d2q9_weights.iter().enumerate() {
                f[node * 9 + alpha] = w;
            }
        }
        Self {
            nx,
            ny,
            f,
            pressure_history: vec![Vec::new(); n_monitor],
        }
    }

    /// Record instantaneous pressure at each monitor point into `pressure_history`.
    ///
    /// Pressure in LBM is p = c_s² * rho where rho = Σ f_α.
    pub fn record_pressure(&mut self, params: &AeroacousticParams) {
        let cs2 = params.c_s * params.c_s;
        for (m, &[ix, iy]) in params.monitor_points.iter().enumerate() {
            let node = iy * self.nx + ix;
            let rho: f64 = (0..9).map(|a| self.f[node * 9 + a]).sum();
            let p = cs2 * rho;
            self.pressure_history[m].push(p);
        }
    }
}

// ============================================================================
// Acoustic quantities
// ============================================================================

/// Instantaneous acoustic pressure fluctuation p' = p_inst − p_ref.
pub fn acoustic_pressure_fluctuation(p_inst: f64, p_ref: f64) -> f64 {
    p_inst - p_ref
}

/// Sound pressure level SPL = 20 log₁₀(p_rms / p_ref) \[dB\].
///
/// Returns −∞ dB (i.e. `f64::NEG_INFINITY`) when `p_rms` ≤ 0 or `p_ref` ≤ 0.
pub fn sound_pressure_level(p_rms: f64, p_ref: f64) -> f64 {
    if p_rms <= 0.0 || p_ref <= 0.0 {
        return f64::NEG_INFINITY;
    }
    20.0 * (p_rms / p_ref).log10()
}

/// Root-mean-square pressure from a time series of pressure samples.
pub fn rms_pressure(history: &[f64]) -> f64 {
    if history.is_empty() {
        return 0.0;
    }
    let mean: f64 = history.iter().sum::<f64>() / history.len() as f64;
    let variance: f64 =
        history.iter().map(|&p| (p - mean).powi(2)).sum::<f64>() / history.len() as f64;
    variance.sqrt()
}

/// Project near-field pressure to far field using 1/r cylindrical/spherical decay.
///
/// `near_field` holds pressure samples at radius `r_near`.  Each sample is
/// scaled by r_near / r_far (1/r amplitude decay).
pub fn far_field_pressure(near_field: &[f64], r_near: f64, r_far: f64) -> f64 {
    if r_far <= 0.0 || r_near <= 0.0 || near_field.is_empty() {
        return 0.0;
    }
    let scale = r_near / r_far;
    let rms = rms_pressure(near_field);
    rms * scale
}

/// Simplified Ffowcs Williams-Hawkings monopole far-field pressure.
///
/// p'(r, t) ≈ (1 / (4πc)) * Q(t_ret) / r
///
/// where Q is the monopole source strength, r is observer distance,
/// c is the speed of sound, and t_ret is the retarded time.
pub fn ffowcs_williams_hawkings_monopole(q: f64, r: f64, c: f64, _t_ret: f64) -> f64 {
    if r <= 0.0 || c <= 0.0 {
        return 0.0;
    }
    q / (4.0 * PI * c * r)
}

/// Strouhal number St = f * L / U.
///
/// `freq` is the shedding frequency, `length` is the characteristic length
/// (e.g. cylinder diameter), and `velocity` is the free-stream speed.
pub fn strouhal_number(freq: f64, length: f64, velocity: f64) -> f64 {
    if velocity == 0.0 {
        return 0.0;
    }
    freq * length / velocity
}

/// Aeolian tone frequency predicted from Strouhal number.
///
/// f = St * U / D
pub fn aeolian_tone_frequency(velocity: f64, diameter: f64, strouhal: f64) -> f64 {
    if diameter == 0.0 {
        return 0.0;
    }
    strouhal * velocity / diameter
}

/// Vortex sound power approximation (Curle's analogy).
///
/// W ≈ ρ * Γ² * U / (2π)
///
/// `circulation` is Γ, `velocity` is the convection speed, `density` is ρ.
pub fn vortex_sound_power(circulation: f64, velocity: f64, density: f64) -> f64 {
    density * circulation * circulation * velocity / (2.0 * PI)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- AeroacousticLBM construction -----------------------------------------

    #[test]
    fn test_new_size() {
        let lbm = AeroacousticLBM::new(4, 4, 2);
        assert_eq!(lbm.nx, 4);
        assert_eq!(lbm.ny, 4);
        assert_eq!(lbm.f.len(), 4 * 4 * 9);
        assert_eq!(lbm.pressure_history.len(), 2);
    }

    #[test]
    fn test_new_weights_sum_to_one() {
        let lbm = AeroacousticLBM::new(3, 3, 0);
        for node in 0..9 {
            let sum: f64 = (0..9).map(|a| lbm.f[node * 9 + a]).sum();
            assert!((sum - 1.0).abs() < 1e-12, "weights sum = {sum}");
        }
    }

    #[test]
    fn test_params_new() {
        let p = AeroacousticParams::new(0.1, 0.01, vec![[1, 1]]);
        assert!((p.c_s - 1.0 / 3.0_f64.sqrt()).abs() < 1e-12);
        assert_eq!(p.monitor_points.len(), 1);
    }

    #[test]
    fn test_record_pressure_length() {
        let params = AeroacousticParams::new(0.1, 0.01, vec![[0, 0], [1, 0]]);
        let mut lbm = AeroacousticLBM::new(4, 4, 2);
        lbm.record_pressure(&params);
        lbm.record_pressure(&params);
        assert_eq!(lbm.pressure_history[0].len(), 2);
        assert_eq!(lbm.pressure_history[1].len(), 2);
    }

    #[test]
    fn test_record_pressure_positive() {
        let params = AeroacousticParams::new(0.1, 0.01, vec![[0, 0]]);
        let mut lbm = AeroacousticLBM::new(4, 4, 1);
        lbm.record_pressure(&params);
        assert!(lbm.pressure_history[0][0] > 0.0);
    }

    // -- acoustic_pressure_fluctuation ----------------------------------------

    #[test]
    fn test_pressure_fluctuation_positive() {
        assert!((acoustic_pressure_fluctuation(1.05, 1.0) - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_fluctuation_zero() {
        assert_eq!(acoustic_pressure_fluctuation(1.0, 1.0), 0.0);
    }

    #[test]
    fn test_pressure_fluctuation_negative() {
        let pf = acoustic_pressure_fluctuation(0.9, 1.0);
        assert!((pf + 0.1).abs() < 1e-12);
    }

    // -- sound_pressure_level -------------------------------------------------

    #[test]
    fn test_spl_reference() {
        // p_rms == p_ref  =>  0 dB
        let spl = sound_pressure_level(1.0, 1.0);
        assert!((spl).abs() < 1e-10);
    }

    #[test]
    fn test_spl_ten_times() {
        // p_rms = 10 * p_ref  =>  20 dB
        let spl = sound_pressure_level(10.0, 1.0);
        assert!((spl - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_spl_zero_p_rms() {
        let spl = sound_pressure_level(0.0, 1.0);
        assert!(spl.is_infinite() && spl < 0.0);
    }

    #[test]
    fn test_spl_zero_p_ref() {
        let spl = sound_pressure_level(1.0, 0.0);
        assert!(spl.is_infinite() && spl < 0.0);
    }

    // -- rms_pressure ---------------------------------------------------------

    #[test]
    fn test_rms_empty() {
        assert_eq!(rms_pressure(&[]), 0.0);
    }

    #[test]
    fn test_rms_constant() {
        // Constant signal has zero RMS fluctuation
        let h = vec![1.0_f64; 100];
        assert!(rms_pressure(&h) < 1e-12);
    }

    #[test]
    fn test_rms_symmetric() {
        // +1 / -1 alternating around 0: variance = 1, rms = 1
        let h = vec![1.0, -1.0, 1.0, -1.0];
        let r = rms_pressure(&h);
        assert!((r - 1.0).abs() < 1e-10, "rms={r}");
    }

    #[test]
    fn test_rms_single() {
        let h = vec![5.0];
        assert!(rms_pressure(&h) < 1e-12);
    }

    // -- far_field_pressure ---------------------------------------------------

    #[test]
    fn test_far_field_decay() {
        let nf = vec![1.0, -1.0, 1.0, -1.0];
        let p_near = far_field_pressure(&nf, 1.0, 1.0);
        let p_far = far_field_pressure(&nf, 1.0, 2.0);
        assert!((p_near - 2.0 * p_far).abs() < 1e-10);
    }

    #[test]
    fn test_far_field_zero_radius() {
        let nf = vec![1.0];
        assert_eq!(far_field_pressure(&nf, 0.0, 1.0), 0.0);
    }

    #[test]
    fn test_far_field_empty() {
        assert_eq!(far_field_pressure(&[], 1.0, 2.0), 0.0);
    }

    // -- ffowcs_williams_hawkings_monopole ------------------------------------

    #[test]
    fn test_fwh_basic() {
        let p = ffowcs_williams_hawkings_monopole(4.0 * PI, 1.0, 1.0, 0.0);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fwh_zero_r() {
        assert_eq!(ffowcs_williams_hawkings_monopole(1.0, 0.0, 1.0, 0.0), 0.0);
    }

    #[test]
    fn test_fwh_zero_c() {
        assert_eq!(ffowcs_williams_hawkings_monopole(1.0, 1.0, 0.0, 0.0), 0.0);
    }

    #[test]
    fn test_fwh_inverse_r() {
        let p1 = ffowcs_williams_hawkings_monopole(1.0, 1.0, 1.0, 0.0);
        let p2 = ffowcs_williams_hawkings_monopole(1.0, 2.0, 1.0, 0.0);
        assert!((p1 - 2.0 * p2).abs() < 1e-12);
    }

    // -- strouhal_number ------------------------------------------------------

    #[test]
    fn test_strouhal_basic() {
        let st = strouhal_number(100.0, 0.01, 2.0);
        assert!((st - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_strouhal_zero_velocity() {
        assert_eq!(strouhal_number(100.0, 0.01, 0.0), 0.0);
    }

    #[test]
    fn test_strouhal_cylinder() {
        // Typical cylinder shedding: St ≈ 0.2
        let st = strouhal_number(4.0, 1.0, 20.0);
        assert!((st - 0.2).abs() < 1e-10);
    }

    // -- aeolian_tone_frequency -----------------------------------------------

    #[test]
    fn test_aeolian_basic() {
        let f = aeolian_tone_frequency(10.0, 1.0, 0.2);
        assert!((f - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_aeolian_zero_diameter() {
        assert_eq!(aeolian_tone_frequency(10.0, 0.0, 0.2), 0.0);
    }

    #[test]
    fn test_aeolian_proportional_to_velocity() {
        let f1 = aeolian_tone_frequency(10.0, 0.1, 0.2);
        let f2 = aeolian_tone_frequency(20.0, 0.1, 0.2);
        assert!((f2 - 2.0 * f1).abs() < 1e-10);
    }

    // -- vortex_sound_power ---------------------------------------------------

    #[test]
    fn test_vortex_power_positive() {
        let w = vortex_sound_power(1.0, 1.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_vortex_power_zero_circulation() {
        assert_eq!(vortex_sound_power(0.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_vortex_power_quadratic_circulation() {
        let w1 = vortex_sound_power(1.0, 1.0, 1.0);
        let w2 = vortex_sound_power(2.0, 1.0, 1.0);
        assert!((w2 - 4.0 * w1).abs() < 1e-10);
    }

    #[test]
    fn test_vortex_power_formula() {
        let gamma = 2.0;
        let u = 3.0;
        let rho = 1.5;
        let expected = rho * gamma * gamma * u / (2.0 * PI);
        assert!((vortex_sound_power(gamma, u, rho) - expected).abs() < 1e-12);
    }
}
