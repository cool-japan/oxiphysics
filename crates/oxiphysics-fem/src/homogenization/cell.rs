// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HomogenizationCell: effective thermal and mechanical property computation.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// HomogenizationCell: effective thermal and mechanical properties
// ---------------------------------------------------------------------------

/// A periodic homogenization cell for computing effective material properties.
///
/// Contains multiple phases with known volume fractions, thermal conductivities,
/// and elastic moduli. Used to compute effective properties of composite materials.
#[allow(dead_code)]
pub struct HomogenizationCell {
    /// Young's moduli of each phase (Pa).
    pub youngs_moduli: Vec<f64>,
    /// Poisson's ratios of each phase.
    pub poisson_ratios: Vec<f64>,
    /// Volume fractions of each phase (must sum to 1).
    pub volume_fractions: Vec<f64>,
    /// Thermal conductivities of each phase (W/m·K).
    pub thermal_conductivities: Vec<f64>,
}

impl HomogenizationCell {
    /// Create a new `HomogenizationCell`.
    ///
    /// # Panics
    /// Panics if the input vectors have different lengths.
    pub fn new(
        youngs_moduli: Vec<f64>,
        poisson_ratios: Vec<f64>,
        volume_fractions: Vec<f64>,
        thermal_conductivities: Vec<f64>,
    ) -> Self {
        assert_eq!(youngs_moduli.len(), poisson_ratios.len());
        assert_eq!(youngs_moduli.len(), volume_fractions.len());
        assert_eq!(youngs_moduli.len(), thermal_conductivities.len());
        Self {
            youngs_moduli,
            poisson_ratios,
            volume_fractions,
            thermal_conductivities,
        }
    }

    /// Compute effective thermal conductivity using the Voigt (upper) and
    /// Reuss (lower) bounds, returning the Hill average.
    ///
    /// Voigt bound (parallel): k_V = Σ f_i * k_i
    /// Reuss bound (series):   k_R = 1 / (Σ f_i / k_i)
    /// Hill average:           k_H = (k_V + k_R) / 2
    ///
    /// Returns `(k_voigt, k_reuss, k_hill)`.
    pub fn compute_effective_thermal_conductivity(&self) -> (f64, f64, f64) {
        let k_voigt: f64 = self
            .volume_fractions
            .iter()
            .zip(self.thermal_conductivities.iter())
            .map(|(&f, &k)| f * k)
            .sum();

        let reuss_sum: f64 = self
            .volume_fractions
            .iter()
            .zip(self.thermal_conductivities.iter())
            .map(|(&f, &k)| if k > 1e-30 { f / k } else { f64::INFINITY })
            .sum();

        let k_reuss = if reuss_sum > 1e-30 {
            1.0 / reuss_sum
        } else {
            0.0
        };
        let k_hill = 0.5 * (k_voigt + k_reuss);

        (k_voigt, k_reuss, k_hill)
    }

    /// Compute effective bulk modulus using Voigt-Reuss-Hill bounds.
    ///
    /// Bulk modulus for each phase: K_i = E_i / (3 * (1 - 2 * nu_i))
    ///
    /// Voigt bound: K_V = Σ f_i * K_i
    /// Reuss bound: K_R = 1 / (Σ f_i / K_i)
    /// Hill average: K_H = (K_V + K_R) / 2
    ///
    /// Returns `(k_voigt, k_reuss, k_hill)`.
    pub fn compute_effective_bulk_modulus(&self) -> (f64, f64, f64) {
        let bulk_moduli: Vec<f64> = self
            .youngs_moduli
            .iter()
            .zip(self.poisson_ratios.iter())
            .map(|(&e, &nu)| {
                let denom = 3.0 * (1.0 - 2.0 * nu);
                if denom.abs() > 1e-30 { e / denom } else { 0.0 }
            })
            .collect();

        let k_voigt: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .map(|(&f, &k)| f * k)
            .sum();

        let reuss_sum: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .map(|(&f, &k)| if k > 1e-30 { f / k } else { f64::INFINITY })
            .sum();

        let k_reuss = if reuss_sum > 1e-30 {
            1.0 / reuss_sum
        } else {
            0.0
        };
        let k_hill = 0.5 * (k_voigt + k_reuss);

        (k_voigt, k_reuss, k_hill)
    }

    /// Compute Hill (Voigt-Reuss-Hill) bounds for effective Young's modulus,
    /// shear modulus, and bulk modulus.
    ///
    /// For a multi-phase composite, this returns:
    /// - Voigt upper bounds for E, G, K
    /// - Reuss lower bounds for E, G, K
    /// - Hill average (E_V + E_R)/2, (G_V + G_R)/2, (K_V + K_R)/2
    ///
    /// Returns `(e_voigt, e_reuss, e_hill, g_voigt, g_reuss, g_hill, k_voigt, k_reuss, k_hill)`.
    #[allow(clippy::type_complexity)]
    pub fn compute_hill_bounds(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let shear_moduli: Vec<f64> = self
            .youngs_moduli
            .iter()
            .zip(self.poisson_ratios.iter())
            .map(|(&e, &nu)| {
                let denom = 2.0 * (1.0 + nu);
                if denom.abs() > 1e-30 { e / denom } else { 0.0 }
            })
            .collect();

        let bulk_moduli: Vec<f64> = self
            .youngs_moduli
            .iter()
            .zip(self.poisson_ratios.iter())
            .map(|(&e, &nu)| {
                let denom = 3.0 * (1.0 - 2.0 * nu);
                if denom.abs() > 1e-30 { e / denom } else { 0.0 }
            })
            .collect();

        // Voigt bounds
        let e_voigt: f64 = self
            .volume_fractions
            .iter()
            .zip(self.youngs_moduli.iter())
            .map(|(&f, &e)| f * e)
            .sum();
        let g_voigt: f64 = self
            .volume_fractions
            .iter()
            .zip(shear_moduli.iter())
            .map(|(&f, &g)| f * g)
            .sum();
        let k_voigt: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .map(|(&f, &k)| f * k)
            .sum();

        // Reuss bounds
        let e_reuss_sum: f64 = self
            .volume_fractions
            .iter()
            .zip(self.youngs_moduli.iter())
            .map(|(&f, &e)| if e > 1e-30 { f / e } else { f64::INFINITY })
            .sum();
        let g_reuss_sum: f64 = self
            .volume_fractions
            .iter()
            .zip(shear_moduli.iter())
            .map(|(&f, &g)| if g > 1e-30 { f / g } else { f64::INFINITY })
            .sum();
        let k_reuss_sum: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .map(|(&f, &k)| if k > 1e-30 { f / k } else { f64::INFINITY })
            .sum();

        let e_reuss = if e_reuss_sum > 1e-30 {
            1.0 / e_reuss_sum
        } else {
            0.0
        };
        let g_reuss = if g_reuss_sum > 1e-30 {
            1.0 / g_reuss_sum
        } else {
            0.0
        };
        let k_reuss = if k_reuss_sum > 1e-30 {
            1.0 / k_reuss_sum
        } else {
            0.0
        };

        let e_hill = 0.5 * (e_voigt + e_reuss);
        let g_hill = 0.5 * (g_voigt + g_reuss);
        let k_hill = 0.5 * (k_voigt + k_reuss);

        (
            e_voigt, e_reuss, e_hill, g_voigt, g_reuss, g_hill, k_voigt, k_reuss, k_hill,
        )
    }

    /// Compute effective thermal expansion coefficient using Turner's model.
    ///
    /// alpha_eff = Σ (f_i * K_i * alpha_i) / Σ (f_i * K_i)
    ///
    /// where K_i is the bulk modulus of phase i.
    pub fn compute_effective_cte(&self, alphas: &[f64]) -> f64 {
        assert_eq!(alphas.len(), self.youngs_moduli.len());
        let bulk_moduli: Vec<f64> = self
            .youngs_moduli
            .iter()
            .zip(self.poisson_ratios.iter())
            .map(|(&e, &nu)| {
                let denom = 3.0 * (1.0 - 2.0 * nu);
                if denom.abs() > 1e-30 { e / denom } else { 0.0 }
            })
            .collect();

        let numerator: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .zip(alphas.iter())
            .map(|((&f, &k), &a)| f * k * a)
            .sum();
        let denominator: f64 = self
            .volume_fractions
            .iter()
            .zip(bulk_moduli.iter())
            .map(|(&f, &k)| f * k)
            .sum();

        if denominator.abs() < 1e-30 {
            0.0
        } else {
            numerator / denominator
        }
    }

    /// Number of phases in this homogenization cell.
    pub fn num_phases(&self) -> usize {
        self.youngs_moduli.len()
    }
}

#[cfg(test)]
mod tests_homogenization_cell {
    use super::*;

    fn two_phase_cell() -> HomogenizationCell {
        // Steel (60%) + Aluminum (40%)
        HomogenizationCell::new(
            vec![200e9, 70e9],
            vec![0.3, 0.33],
            vec![0.6, 0.4],
            vec![50.0, 200.0], // W/m·K
        )
    }

    // ── Thermal conductivity tests ───────────────────────────────────────

    #[test]
    fn test_thermal_conductivity_voigt_ge_reuss() {
        let cell = two_phase_cell();
        let (k_v, k_r, _k_h) = cell.compute_effective_thermal_conductivity();
        assert!(
            k_v >= k_r,
            "Voigt bound must be >= Reuss bound: k_V={}, k_R={}",
            k_v,
            k_r
        );
    }

    #[test]
    fn test_thermal_conductivity_hill_between_bounds() {
        let cell = two_phase_cell();
        let (k_v, k_r, k_h) = cell.compute_effective_thermal_conductivity();
        assert!(
            k_h >= k_r && k_h <= k_v,
            "Hill average must be between Reuss and Voigt: k_R={}, k_H={}, k_V={}",
            k_r,
            k_h,
            k_v
        );
    }

    #[test]
    fn test_thermal_conductivity_single_phase_exact() {
        let cell = HomogenizationCell::new(vec![200e9], vec![0.3], vec![1.0], vec![50.0]);
        let (k_v, k_r, k_h) = cell.compute_effective_thermal_conductivity();
        assert!((k_v - 50.0).abs() < 1e-10, "Voigt single phase: {}", k_v);
        assert!((k_r - 50.0).abs() < 1e-10, "Reuss single phase: {}", k_r);
        assert!((k_h - 50.0).abs() < 1e-10, "Hill single phase: {}", k_h);
    }

    #[test]
    fn test_thermal_conductivity_voigt_is_weighted_average() {
        let cell = two_phase_cell();
        let (k_v, _k_r, _k_h) = cell.compute_effective_thermal_conductivity();
        let expected = 0.6 * 50.0 + 0.4 * 200.0;
        assert!(
            (k_v - expected).abs() < 1e-9,
            "Voigt k = {}, expected {}",
            k_v,
            expected
        );
    }

    // ── Bulk modulus tests ───────────────────────────────────────────────

    #[test]
    fn test_bulk_modulus_voigt_ge_reuss() {
        let cell = two_phase_cell();
        let (k_v, k_r, _k_h) = cell.compute_effective_bulk_modulus();
        assert!(
            k_v >= k_r,
            "Voigt K must be >= Reuss K: K_V={}, K_R={}",
            k_v,
            k_r
        );
    }

    #[test]
    fn test_bulk_modulus_hill_between_bounds() {
        let cell = two_phase_cell();
        let (k_v, k_r, k_h) = cell.compute_effective_bulk_modulus();
        assert!(
            k_h >= k_r && k_h <= k_v,
            "Hill K must be between bounds: K_R={:.3e}, K_H={:.3e}, K_V={:.3e}",
            k_r,
            k_h,
            k_v
        );
    }

    #[test]
    fn test_bulk_modulus_single_phase_exact() {
        let e = 200e9;
        let nu = 0.3;
        let k_expected = e / (3.0 * (1.0 - 2.0 * nu));
        let cell = HomogenizationCell::new(vec![e], vec![nu], vec![1.0], vec![50.0]);
        let (k_v, k_r, k_h) = cell.compute_effective_bulk_modulus();
        assert!(
            (k_v - k_expected).abs() / k_expected < 1e-9,
            "Voigt K: {} vs {}",
            k_v,
            k_expected
        );
        assert!(
            (k_r - k_expected).abs() / k_expected < 1e-9,
            "Reuss K: {} vs {}",
            k_r,
            k_expected
        );
        assert!(
            (k_h - k_expected).abs() / k_expected < 1e-9,
            "Hill K: {} vs {}",
            k_h,
            k_expected
        );
    }

    // ── Hill bounds (full tensor) tests ──────────────────────────────────

    #[test]
    fn test_hill_bounds_e_voigt_ge_e_reuss() {
        let cell = two_phase_cell();
        let (e_v, e_r, _e_h, _gv, _gr, _gh, _kv, _kr, _kh) = cell.compute_hill_bounds();
        assert!(e_v >= e_r, "E_Voigt must be >= E_Reuss: {} vs {}", e_v, e_r);
    }

    #[test]
    fn test_hill_bounds_e_hill_between() {
        let cell = two_phase_cell();
        let (e_v, e_r, e_h, _gv, _gr, _gh, _kv, _kr, _kh) = cell.compute_hill_bounds();
        assert!(
            e_h >= e_r && e_h <= e_v,
            "E_Hill={} must be between E_R={} and E_V={}",
            e_h,
            e_r,
            e_v
        );
    }

    #[test]
    fn test_hill_bounds_single_phase_collapses() {
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let cell = HomogenizationCell::new(vec![e], vec![nu], vec![1.0], vec![50.0]);
        let (e_v, e_r, e_h, _gv, _gr, _gh, _kv, _kr, _kh) = cell.compute_hill_bounds();
        assert!((e_v - e).abs() / e < 1e-9, "Single phase E_Voigt: {}", e_v);
        assert!((e_r - e).abs() / e < 1e-9, "Single phase E_Reuss: {}", e_r);
        assert!((e_h - e).abs() / e < 1e-9, "Single phase E_Hill: {}", e_h);
    }
}
