// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum materials properties: band structure, Fermi statistics, carrier
//! transport, thermoelectric effects, and optical response.
//!
//! All energies are in electron-volts (eV) unless stated otherwise.
//! Temperature is in Kelvin.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

/// Boltzmann constant in eV/K.
pub const K_B_EV: f64 = 8.617_333_262e-5;
/// Boltzmann constant in J/K.
pub const K_B_J: f64 = 1.380_649e-23;
/// Electron charge in Coulombs.
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Free electron mass in kg.
pub const M_ELECTRON: f64 = 9.109_383_701_5e-31;
/// Reduced Planck constant in J·s.
pub const HBAR: f64 = 1.054_571_817e-34;

// ─────────────────────────────────────────────────────────────────────────────
// Core struct
// ─────────────────────────────────────────────────────────────────────────────

/// Quantum material descriptor capturing electronic structure parameters.
///
/// All energy quantities are in eV.
#[derive(Debug, Clone)]
pub struct QuantumMaterial {
    /// Band gap energy (eV). 0.0 for metals.
    pub band_gap: f64,
    /// Fermi energy measured from valence-band maximum (eV).
    pub fermi_energy: f64,
    /// Effective carrier mass as a fraction of the free electron mass (dimensionless).
    pub effective_mass: f64,
    /// Carrier mobility at 300 K in m²/(V·s).
    pub mobility_300k: f64,
    /// Optical scattering time (Drude relaxation time) in seconds.
    pub scattering_time: f64,
    /// Material name.
    pub name: String,
}

impl QuantumMaterial {
    /// Create a new `QuantumMaterial`.
    ///
    /// # Arguments
    /// * `band_gap` — band gap in eV
    /// * `fermi_energy` — Fermi energy in eV
    /// * `effective_mass` — effective mass ratio m*/m_e
    /// * `mobility_300k` — mobility at 300 K in m²/(V·s)
    /// * `scattering_time` — Drude relaxation time in seconds
    /// * `name` — material identifier
    pub fn new(
        band_gap: f64,
        fermi_energy: f64,
        effective_mass: f64,
        mobility_300k: f64,
        scattering_time: f64,
        name: impl Into<String>,
    ) -> Self {
        Self {
            band_gap,
            fermi_energy,
            effective_mass,
            mobility_300k,
            scattering_time,
            name: name.into(),
        }
    }

    /// Silicon preset (indirect band gap semiconductor).
    pub fn silicon() -> Self {
        Self::new(1.12, 0.56, 0.26, 0.135, 3.0e-13, "Silicon")
    }

    /// GaAs preset (direct band gap III-V semiconductor).
    pub fn gaas() -> Self {
        Self::new(1.42, 0.71, 0.063, 0.85, 1.0e-13, "GaAs")
    }

    /// Copper preset (metal, effectively zero band gap).
    pub fn copper() -> Self {
        Self::new(0.0, 7.04, 1.0, 3.2e-3, 2.5e-14, "Copper")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Statistical functions
// ─────────────────────────────────────────────────────────────────────────────

/// Fermi-Dirac distribution function f(E, E_F, T).
///
/// Returns the probability that a state at energy `energy` (eV) is occupied
/// at temperature `temp_k` (K) with Fermi level `fermi_energy` (eV).
///
/// At T → 0 K this becomes a Heaviside step: f = 1 for E < E_F, f = 0 for E > E_F.
pub fn fermi_dirac(energy: f64, fermi_energy: f64, temp_k: f64) -> f64 {
    if temp_k <= 0.0 {
        if energy < fermi_energy {
            1.0
        } else if (energy - fermi_energy).abs() < 1e-12 {
            0.5
        } else {
            0.0
        }
    } else {
        let kbt = K_B_EV * temp_k;
        let x = (energy - fermi_energy) / kbt;
        // Numerically stable form to avoid overflow.
        if x > 500.0 {
            0.0
        } else if x < -500.0 {
            1.0
        } else {
            1.0 / (1.0 + x.exp())
        }
    }
}

/// Parabolic-band density of states (DOS) for a 3D semiconductor (states/eV/m³).
///
/// D(E) = (1/2π²) × (2 m* / ℏ²)^(3/2) × sqrt(E - E_c) for E > E_c
///
/// # Arguments
/// * `energy` — energy of interest (eV), relative to the band edge
/// * `effective_mass_ratio` — m*/m_e (dimensionless)
///
/// Returns zero for `energy` ≤ 0 (below the band edge).
pub fn density_of_states(energy: f64, effective_mass_ratio: f64) -> f64 {
    if energy <= 0.0 {
        return 0.0;
    }
    let m_star = effective_mass_ratio * M_ELECTRON;
    let energy_j = energy * E_CHARGE; // eV → J
    // D(E) = (1/(2π²)) * (2 m* / ħ²)^(3/2) * sqrt(E)
    let prefactor = (1.0 / (2.0 * PI * PI)) * (2.0 * m_star / (HBAR * HBAR)).powf(1.5);
    let dos_per_j = prefactor * energy_j.sqrt();
    // Convert states/J/m³ → states/eV/m³
    dos_per_j * E_CHARGE
}

/// Electron carrier concentration n (m⁻³) using the Maxwell-Boltzmann approximation
/// (non-degenerate semiconductor).
///
/// n = N_c × exp(-(E_c - E_F) / k_B T)
///
/// where N_c = 2 (2π m* k_B T / h²)^(3/2) is the effective density of states.
///
/// # Arguments
/// * `band_gap` — band gap E_g (eV)
/// * `fermi_energy` — Fermi level E_F measured from valence-band maximum (eV)
/// * `effective_mass_ratio` — m*/m_e
/// * `temp_k` — temperature (K)
pub fn carrier_concentration(
    band_gap: f64,
    fermi_energy: f64,
    effective_mass_ratio: f64,
    temp_k: f64,
) -> f64 {
    if temp_k <= 0.0 {
        return 0.0;
    }
    let m_star = effective_mass_ratio * M_ELECTRON;
    let kbt_j = K_B_J * temp_k;
    let kbt_ev = K_B_EV * temp_k;
    // Effective DOS at conduction band edge.
    let h = 2.0 * PI * HBAR;
    let n_c = 2.0 * (2.0 * PI * m_star * kbt_j / (h * h)).powf(1.5);
    // Conduction band edge relative to valence band maximum.
    let e_c = band_gap;
    // Activation exponent.
    let delta = (e_c - fermi_energy) / kbt_ev;
    if delta > 500.0 {
        return 0.0;
    }
    n_c * (-delta).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Transport properties
// ─────────────────────────────────────────────────────────────────────────────

/// Phonon-scattering mobility as a function of temperature.
///
/// Uses the empirical power-law model: μ(T) = μ₃₀₀ × (300 / T)^α
/// where α ≈ 1.5 for acoustic-phonon scattering (Boltzmann limit).
///
/// # Arguments
/// * `mobility_300k` — mobility at 300 K (m²/(V·s))
/// * `temp_k` — temperature (K)
/// * `alpha` — scattering exponent (typically 1.5 for acoustic phonons)
///
/// Returns mobility in m²/(V·s).
pub fn mobility_temperature(mobility_300k: f64, temp_k: f64, alpha: f64) -> f64 {
    if temp_k <= 0.0 {
        return 0.0;
    }
    mobility_300k * (300.0 / temp_k).powf(alpha)
}

/// Seebeck coefficient (thermopower) S in V/K for a parabolic-band semiconductor.
///
/// S = -(k_B / e) × \[(E_c - E_F)/(k_B T) + (r + 5/2)\]
///
/// where r is the scattering exponent (r = -1/2 for acoustic phonons).
///
/// # Arguments
/// * `band_gap` — E_g (eV)
/// * `fermi_energy` — E_F from valence-band maximum (eV)
/// * `temp_k` — temperature (K)
/// * `r_scatter` — scattering exponent (default: -0.5)
pub fn seebeck_coefficient(band_gap: f64, fermi_energy: f64, temp_k: f64, r_scatter: f64) -> f64 {
    if temp_k <= 0.0 {
        return 0.0;
    }
    let kbt_ev = K_B_EV * temp_k;
    let e_c = band_gap;
    let eta = (e_c - fermi_energy) / kbt_ev; // reduced Fermi energy
    // S = -(k_B/e) * [eta + (r + 5/2)]

    -(K_B_EV / 1.0) * (eta + r_scatter + 2.5) // units: eV/K (numerically same as V/K since e cancels)
}

/// Hall coefficient R_H (m³/C) for an n-type semiconductor.
///
/// R_H = -r_H / (n × e)
///
/// where r_H is the Hall scattering factor (≈ 1 for acoustic phonons) and n is
/// the electron concentration.
///
/// # Arguments
/// * `carrier_density` — electron concentration (m⁻³)
/// * `hall_factor` — Hall scattering factor r_H (dimensionless, typically ≈ 1)
///
/// Returns R_H in m³/C (negative for electrons, positive for holes).
pub fn hall_coefficient(carrier_density: f64, hall_factor: f64) -> f64 {
    if carrier_density <= 0.0 {
        return 0.0;
    }
    -hall_factor / (carrier_density * E_CHARGE)
}

// ─────────────────────────────────────────────────────────────────────────────
// Optical properties
// ─────────────────────────────────────────────────────────────────────────────

/// Drude-model optical conductivity σ(ω) (S/m) at angular frequency ω (rad/s).
///
/// σ(ω) = σ₀ / (1 + (ω τ)²)
///
/// where σ₀ = n e² τ / m* is the DC conductivity.
///
/// # Arguments
/// * `carrier_density` — carrier concentration (m⁻³)
/// * `effective_mass_ratio` — m*/m_e
/// * `scattering_time` — Drude relaxation time τ (s)
/// * `omega` — angular frequency (rad/s)
///
/// Returns the real part of the optical conductivity in S/m.
pub fn optical_conductivity(
    carrier_density: f64,
    effective_mass_ratio: f64,
    scattering_time: f64,
    omega: f64,
) -> f64 {
    let m_star = effective_mass_ratio * M_ELECTRON;
    let sigma0 = carrier_density * E_CHARGE * E_CHARGE * scattering_time / m_star;
    let omega_tau = omega * scattering_time;
    sigma0 / (1.0 + omega_tau * omega_tau)
}

/// Plasma frequency ω_p (rad/s) of a free-carrier system.
///
/// ω_p = sqrt(n e² / (ε₀ m*))
///
/// # Arguments
/// * `carrier_density` — n (m⁻³)
/// * `effective_mass_ratio` — m*/m_e
/// * `epsilon_r` — relative permittivity of the host lattice
pub fn plasma_frequency(carrier_density: f64, effective_mass_ratio: f64, epsilon_r: f64) -> f64 {
    let epsilon_0 = 8.854_187_817e-12_f64; // F/m
    let m_star = effective_mass_ratio * M_ELECTRON;
    (carrier_density * E_CHARGE * E_CHARGE / (epsilon_0 * epsilon_r * m_star)).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // 1. Fermi-Dirac at T=0: below Fermi level → 1.
    #[test]
    fn test_fermi_dirac_t0_below() {
        let f = fermi_dirac(0.5, 1.0, 0.0);
        assert!(
            (f - 1.0).abs() < EPS,
            "FD at T=0, E<E_F should be 1, got {f}"
        );
    }

    // 2. Fermi-Dirac at T=0: above Fermi level → 0.
    #[test]
    fn test_fermi_dirac_t0_above() {
        let f = fermi_dirac(1.5, 1.0, 0.0);
        assert!(f.abs() < EPS, "FD at T=0, E>E_F should be 0, got {f}");
    }

    // 3. Fermi-Dirac at T=0: exactly at Fermi level → 0.5.
    #[test]
    fn test_fermi_dirac_t0_at_ef() {
        let f = fermi_dirac(1.0, 1.0, 0.0);
        assert!(
            (f - 0.5).abs() < EPS,
            "FD at T=0, E=E_F should be 0.5, got {f}"
        );
    }

    // 4. Fermi-Dirac at finite T: at Fermi level → exactly 0.5.
    #[test]
    fn test_fermi_dirac_at_fermi_level_finite_t() {
        let f = fermi_dirac(1.0, 1.0, 300.0);
        assert!(
            (f - 0.5).abs() < EPS,
            "FD at E=E_F for any T should be 0.5, got {f}"
        );
    }

    // 5. Fermi-Dirac: high T thermal smearing, f at E_F+5kT ≈ 0.0067.
    #[test]
    fn test_fermi_dirac_thermal_tail() {
        let kbt = K_B_EV * 1000.0;
        let ef = 1.0;
        let e = ef + 5.0 * kbt;
        let f = fermi_dirac(e, ef, 1000.0);
        let expected = 1.0 / (1.0 + 5.0_f64.exp());
        assert!(
            (f - expected).abs() < 1e-12,
            "FD thermal tail mismatch: {f} vs {expected}"
        );
    }

    // 6. Fermi-Dirac: very negative exponent → 1.
    #[test]
    fn test_fermi_dirac_large_negative() {
        let f = fermi_dirac(-10.0, 1.0, 0.01); // Very cold, way below E_F
        assert!(f > 0.99, "FD way below E_F at low T should be ~1, got {f}");
    }

    // 7. DOS is zero at band edge (E=0).
    #[test]
    fn test_dos_zero_at_band_edge() {
        let d = density_of_states(0.0, 0.26);
        assert!(d.abs() < EPS, "DOS at band edge should be 0, got {d}");
    }

    // 8. DOS increases with energy (parabolic band).
    #[test]
    fn test_dos_increases_with_energy() {
        let d1 = density_of_states(0.1, 0.26);
        let d2 = density_of_states(0.5, 0.26);
        assert!(
            d2 > d1,
            "DOS should increase with energy in parabolic band, d1={d1}, d2={d2}"
        );
    }

    // 9. DOS scales with effective mass (heavier mass → more states).
    #[test]
    fn test_dos_scales_with_mass() {
        let d_light = density_of_states(0.3, 0.1);
        let d_heavy = density_of_states(0.3, 1.0);
        assert!(
            d_heavy > d_light,
            "Heavier effective mass should give larger DOS"
        );
    }

    // 10. Carrier concentration increases with temperature.
    #[test]
    fn test_carrier_concentration_increases_with_temp() {
        let n_low = carrier_concentration(1.12, 0.56, 0.26, 200.0);
        let n_high = carrier_concentration(1.12, 0.56, 0.26, 400.0);
        assert!(
            n_high > n_low,
            "Carrier conc should increase with T: n_low={n_low}, n_high={n_high}"
        );
    }

    // 11. Carrier concentration at T=0 is zero.
    #[test]
    fn test_carrier_concentration_t0() {
        let n = carrier_concentration(1.12, 0.56, 0.26, 0.0);
        assert!(n.abs() < EPS, "Carrier conc at T=0 should be 0, got {n}");
    }

    // 12. Carrier concentration decreases with larger band gap.
    #[test]
    fn test_carrier_concentration_band_gap_effect() {
        let n_small_gap = carrier_concentration(0.5, 0.25, 0.26, 300.0);
        let n_large_gap = carrier_concentration(2.0, 1.0, 0.26, 300.0);
        assert!(
            n_small_gap > n_large_gap,
            "Smaller band gap should yield higher carrier density"
        );
    }

    // 13. Mobility decreases with temperature (alpha > 0).
    #[test]
    fn test_mobility_decreases_with_temp() {
        let mu_low = mobility_temperature(0.135, 200.0, 1.5);
        let mu_high = mobility_temperature(0.135, 400.0, 1.5);
        assert!(
            mu_high < mu_low,
            "Mobility should decrease with temperature"
        );
    }

    // 14. Mobility at 300 K equals the given reference value.
    #[test]
    fn test_mobility_at_300k() {
        let mu = mobility_temperature(0.135, 300.0, 1.5);
        assert!(
            (mu - 0.135).abs() < 1e-12,
            "Mobility at 300K should be 0.135, got {mu}"
        );
    }

    // 15. Mobility at T=0 returns 0.
    #[test]
    fn test_mobility_t0() {
        let mu = mobility_temperature(0.135, 0.0, 1.5);
        assert!(mu.abs() < EPS, "Mobility at T=0 should be 0, got {mu}");
    }

    // 16. Seebeck coefficient is non-zero for a typical semiconductor.
    #[test]
    fn test_seebeck_nonzero() {
        let s = seebeck_coefficient(1.12, 0.56, 300.0, -0.5);
        assert!(
            s.abs() > 1e-6,
            "Seebeck coefficient should be non-zero, got {s}"
        );
    }

    // 17. Seebeck coefficient at T=0 returns 0.
    #[test]
    fn test_seebeck_t0() {
        let s = seebeck_coefficient(1.12, 0.56, 0.0, -0.5);
        assert!(s.abs() < EPS, "Seebeck at T=0 should be 0, got {s}");
    }

    // 18. Seebeck sign: for E_F below conduction band, S < 0 (n-type).
    #[test]
    fn test_seebeck_n_type_negative() {
        // E_F (0.56 eV) < E_c (band_gap = 1.12 eV): n-type → negative Seebeck
        let s = seebeck_coefficient(1.12, 0.56, 300.0, -0.5);
        assert!(s < 0.0, "n-type Seebeck should be negative, got {s}");
    }

    // 19. Hall coefficient is negative for electrons.
    #[test]
    fn test_hall_coefficient_negative_for_electrons() {
        let rh = hall_coefficient(1e22, 1.0);
        assert!(
            rh < 0.0,
            "Hall coefficient for n-type should be negative, got {rh}"
        );
    }

    // 20. Hall coefficient magnitude: R_H = 1 / (n e) for r_H = 1.
    #[test]
    fn test_hall_coefficient_magnitude() {
        let n = 1e22_f64;
        let rh = hall_coefficient(n, 1.0);
        let expected = -1.0 / (n * E_CHARGE);
        assert!(
            (rh - expected).abs() < EPS.abs(),
            "Hall coefficient mismatch: {rh} vs {expected}"
        );
    }

    // 21. Hall coefficient for zero carrier density returns 0.
    #[test]
    fn test_hall_coefficient_zero_density() {
        let rh = hall_coefficient(0.0, 1.0);
        assert!(
            rh.abs() < EPS,
            "Hall coefficient for zero density should be 0, got {rh}"
        );
    }

    // 22. Optical conductivity at DC (ω=0) equals σ₀ = n e² τ / m*.
    #[test]
    fn test_optical_conductivity_dc() {
        let n = 8.49e28_f64; // copper-like carrier density
        let m_ratio = 1.0;
        let tau = 2.5e-14;
        let m_star = m_ratio * M_ELECTRON;
        let sigma0 = n * E_CHARGE * E_CHARGE * tau / m_star;
        let sigma = optical_conductivity(n, m_ratio, tau, 0.0);
        assert!(
            (sigma - sigma0).abs() / sigma0 < 1e-10,
            "DC conductivity mismatch"
        );
    }

    // 23. Optical conductivity decreases at high frequency.
    #[test]
    fn test_optical_conductivity_decreases_at_high_freq() {
        let n = 8.49e28_f64;
        let sigma_dc = optical_conductivity(n, 1.0, 2.5e-14, 0.0);
        let sigma_hf = optical_conductivity(n, 1.0, 2.5e-14, 1e15);
        assert!(
            sigma_hf < sigma_dc,
            "Conductivity should decrease at high frequency"
        );
    }

    // 24. Plasma frequency is positive for finite carrier density.
    #[test]
    fn test_plasma_frequency_positive() {
        let wp = plasma_frequency(8.49e28, 1.0, 1.0);
        assert!(wp > 0.0, "Plasma frequency should be positive, got {wp}");
    }

    // 25. QuantumMaterial presets can be constructed without panicking.
    #[test]
    fn test_presets_construction() {
        let si = QuantumMaterial::silicon();
        let gaas = QuantumMaterial::gaas();
        let cu = QuantumMaterial::copper();
        assert!(si.band_gap > 0.0);
        assert!(gaas.band_gap > 0.0);
        assert!(cu.band_gap == 0.0);
    }

    // 26. Higher effective mass raises the carrier concentration (heavier N_c).
    #[test]
    fn test_carrier_concentration_mass_effect() {
        let n_light = carrier_concentration(1.12, 0.56, 0.1, 300.0);
        let n_heavy = carrier_concentration(1.12, 0.56, 1.0, 300.0);
        assert!(
            n_heavy > n_light,
            "Heavier effective mass should increase N_c and hence n"
        );
    }

    // 27. Fermi-Dirac is monotonically decreasing in energy at finite T.
    #[test]
    fn test_fermi_dirac_monotone() {
        let ef = 1.0;
        let t = 300.0;
        let f1 = fermi_dirac(0.5, ef, t);
        let f2 = fermi_dirac(1.0, ef, t);
        let f3 = fermi_dirac(1.5, ef, t);
        assert!(
            f1 > f2 && f2 > f3,
            "FD must be decreasing in energy: f1={f1} f2={f2} f3={f3}"
        );
    }

    // 28. DOS has correct sqrt(E) dependence.
    #[test]
    fn test_dos_sqrt_scaling() {
        let d1 = density_of_states(1.0, 0.26);
        let d4 = density_of_states(4.0, 0.26);
        let ratio = d4 / d1;
        // sqrt(4)/sqrt(1) = 2
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "DOS should scale as sqrt(E): ratio={ratio}"
        );
    }

    // 29. Optical conductivity at ω = 1/τ equals σ₀/2 (half-power point).
    #[test]
    fn test_optical_conductivity_half_power() {
        let n = 1e28_f64;
        let tau = 1e-14;
        let omega = 1.0 / tau; // ω τ = 1
        let sigma = optical_conductivity(n, 1.0, tau, omega);
        let sigma0 = optical_conductivity(n, 1.0, tau, 0.0);
        assert!(
            (sigma - sigma0 / 2.0).abs() / sigma0 < 1e-12,
            "At ωτ=1, σ should be σ₀/2"
        );
    }

    // 30. Silicon seebeck magnitude is in physically reasonable range (> 100 µV/K).
    #[test]
    fn test_seebeck_silicon_reasonable() {
        // Typical silicon at moderate doping has |S| ~ 300-900 µV/K
        let s = seebeck_coefficient(1.12, 0.56, 300.0, -0.5);
        assert!(
            s.abs() > 1e-4,
            "Silicon Seebeck should be > 100 µV/K, got {s} V/K"
        );
    }
}
