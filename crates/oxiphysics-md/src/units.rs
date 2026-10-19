// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit-conversion constants shared by the MD modules.
//!
//! Several modules work in the "metal"-like unit system of atomistic MD codes:
//! masses in amu (g mol⁻¹), lengths in Å and times in ps, so that velocities are
//! in Å/ps.  The derived energy unit `amu·Å²/ps²` is `10 J mol⁻¹` exactly
//! (`1 g/mol · (100 m/s)² = 10 J/mol`), which fixes the conversion factors below.

/// Energy of one `amu·(Å/ps)²` expressed in electron-volts.
///
/// 1 Å/ps = 100 m/s, so 1 amu·(Å/ps)² = 1.66054e-27 kg · 1e4 m²/s² = 1.66054e-23 J
/// = 1.036427e-4 eV (taking 1 amu = 1 g/mol exactly, as MD codes do).  This is
/// the `mvv2e` factor of the "metal" unit system.
pub const AMU_A2_PS2_TO_EV: f64 = 1.036_426_965_626_239_8e-4;

/// Boltzmann constant in `amu·Å²/ps²` per kelvin.
///
/// `k_B = R / 10` with `R = N_A k_B = 8.314 462 618 153 24 J mol⁻¹ K⁻¹` (exact,
/// 2019 SI), because one `amu·Å²/ps²` is `10 J mol⁻¹`.  Equivalently
/// `8.617333e-5 eV/K` divided by [`AMU_A2_PS2_TO_EV`].
pub const K_B_AMU_A2_PS2: f64 = 0.831_446_261_815_324;

/// Number of `amu·Å²/ps²` in one kcal mol⁻¹ (`4184 J mol⁻¹ / 10 J mol⁻¹`).
///
/// Multiplying a force in kcal mol⁻¹ Å⁻¹ by this factor and dividing by a mass in
/// amu gives an acceleration in Å/ps².
pub const KCAL_MOL_TO_AMU_A2_PS2: f64 = 418.4;

/// Boltzmann constant in kcal mol⁻¹ K⁻¹ (`R / 4184`).
pub const K_B_KCAL_MOL_K: f64 = 1.987_204_258_640_832e-3;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_b_is_the_gas_constant_over_ten() {
        // R = N_A k_B = 6.02214076e23 * 1.380649e-23 = 8.31446261815324 J/(mol K).
        let gas_constant = 6.022_140_76e23_f64 * 1.380_649e-23_f64;
        assert!((K_B_AMU_A2_PS2 - gas_constant / 10.0).abs() < 1e-12);
    }

    #[test]
    fn k_b_in_ev_is_consistent_with_the_energy_unit() {
        let k_b_ev = K_B_AMU_A2_PS2 * AMU_A2_PS2_TO_EV;
        assert!((k_b_ev - 8.617_333_262e-5).abs() < 1e-13, "{k_b_ev}");
    }

    #[test]
    fn kcal_per_mole_constants_agree() {
        // k_B in kcal/mol/K times the conversion factor is k_B in amu Å²/ps²/K.
        let converted = K_B_KCAL_MOL_K * KCAL_MOL_TO_AMU_A2_PS2;
        assert!((converted - K_B_AMU_A2_PS2).abs() < 1e-12, "{converted}");
        // 1 kcal/mol = 4184 J/mol = 418.4 * 10 J/mol.
        assert_eq!(KCAL_MOL_TO_AMU_A2_PS2 * 10.0, 4184.0);
        assert!((K_B_KCAL_MOL_K - 0.001_987_204).abs() < 1e-9);
    }
}
