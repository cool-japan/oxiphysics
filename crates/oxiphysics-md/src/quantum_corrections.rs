// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum corrections for molecular dynamics simulations.
//!
//! This module provides quantum-mechanical corrections to classical MD observables:
//!
//! - [`WignerKirkwoodExpansion`] — ℏ²/12mkT correction to the classical partition function
//! - [`QuantumThermalEnergy`] — Bose-Einstein / Fermi-Dirac occupation + zero-point energy
//! - [`PathIntegralEstimator`] — centroid virial estimator for PIMD kinetic energy
//! - [`TunnellingRate`] — WKB tunnelling probability
//! - [`ZeroPointEnergy`] — harmonic ZPE: E₀ = ℏω/2; anharmonic correction
//! - [`IsotopeEffect`] — kinetic isotope effect from mass-dependent ZPE difference
//! - [`QuantumCorrection`] — Feynman-Hibbs effective potential: V_eff = V + ℏ²∇²V/24mkT
//! - [`DebyeModel`] — Debye heat capacity Cv(T)
//! - [`EinsteinModel`] — Einstein heat capacity model
//! - [`NuclearQuantumEffect`] — proton tunnelling fraction, H/D isotope fractionation

use crate::solid_state_md::debye_cv_function;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Reduced Planck constant (J·s).
pub const HBAR: f64 = 1.054_571_817e-34;

/// Boltzmann constant (J K⁻¹).
pub const KB: f64 = 1.380_649e-23;

/// Proton mass (kg).
pub const M_PROTON: f64 = 1.672_621_923e-27;

/// Deuteron mass (kg).
pub const M_DEUTERON: f64 = 3.344_493_695e-27;

/// Planck constant h (J·s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;

// ---------------------------------------------------------------------------
// WignerKirkwoodExpansion
// ---------------------------------------------------------------------------

/// Wigner-Kirkwood quantum correction to the classical partition function.
///
/// The leading correction to the free energy is:
///   ΔF = −ℏ²/(12mkT) ⟨∇²U⟩
///
/// See: Wigner (1932), Kirkwood (1933).
#[derive(Debug, Clone)]
pub struct WignerKirkwoodExpansion {
    /// Particle mass (kg).
    pub mass: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl WignerKirkwoodExpansion {
    /// Create a new WK expansion for the given mass and temperature.
    pub fn new(mass: f64, temperature: f64) -> Self {
        Self { mass, temperature }
    }

    /// Thermal de Broglie wavelength: Λ = ℏ√(2π/mkBT)  (m).
    pub fn thermal_de_broglie(&self) -> f64 {
        HBAR * (2.0 * PI / (self.mass * KB * self.temperature)).sqrt()
    }

    /// WK quantum correction prefactor: ℏ²/(12 m kB T).
    pub fn wk_prefactor(&self) -> f64 {
        HBAR * HBAR / (12.0 * self.mass * KB * self.temperature)
    }

    /// Free energy correction: ΔF = prefactor × ⟨∇²U⟩.
    pub fn free_energy_correction(&self, mean_laplacian_u: f64) -> f64 {
        -self.wk_prefactor() * mean_laplacian_u
    }

    /// Quantum correction to internal energy: ΔU = ℏ²/(12mkBT) × ⟨(∇U)²/kBT − ∇²U⟩.
    pub fn internal_energy_correction(&self, mean_grad_u_sq: f64, mean_laplacian_u: f64) -> f64 {
        let pf = self.wk_prefactor();
        pf * (mean_grad_u_sq / (KB * self.temperature) - mean_laplacian_u)
    }

    /// Ratio of quantum to classical partition function (first order):
    ///   Z_q / Z_cl ≈ 1 − β² ℏ²/(12m) ⟨∇²U⟩
    pub fn partition_function_ratio(&self, mean_laplacian_u: f64) -> f64 {
        let beta = 1.0 / (KB * self.temperature);
        1.0 - beta * beta * HBAR * HBAR / (12.0 * self.mass) * mean_laplacian_u
    }
}

// ---------------------------------------------------------------------------
// QuantumThermalEnergy
// ---------------------------------------------------------------------------

/// Quantum thermal energy including zero-point contributions.
///
/// Handles both Bose-Einstein (bosons) and Fermi-Dirac (fermions) statistics.
#[derive(Debug, Clone)]
pub struct QuantumThermalEnergy {
    /// Angular frequency ω (rad s⁻¹).
    pub omega: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl QuantumThermalEnergy {
    /// Create a new quantum thermal energy calculator.
    pub fn new(omega: f64, temperature: f64) -> Self {
        Self { omega, temperature }
    }

    /// Zero-point energy: E₀ = ℏω/2  (J).
    pub fn zero_point_energy(&self) -> f64 {
        0.5 * HBAR * self.omega
    }

    /// Dimensionless parameter x = ℏω/(kBT).
    pub fn x_param(&self) -> f64 {
        HBAR * self.omega / (KB * self.temperature)
    }

    /// Bose-Einstein mean occupation number: n̄ = 1/(e^x − 1).
    pub fn bose_einstein_occupation(&self) -> f64 {
        let x = self.x_param();
        if x > 700.0 {
            return 0.0;
        }
        1.0 / (x.exp() - 1.0)
    }

    /// Fermi-Dirac occupation at energy ℏω with chemical potential μ=0:
    ///   n̄ = 1/(e^x + 1).
    pub fn fermi_dirac_occupation(&self) -> f64 {
        let x = self.x_param();
        1.0 / (x.exp() + 1.0)
    }

    /// Total quantum harmonic oscillator energy (boson):
    ///   E = ℏω (n̄ + 1/2) = ℏω/2 · coth(x/2).
    pub fn boson_energy(&self) -> f64 {
        let x = self.x_param();
        if x < 1e-10 {
            return KB * self.temperature; // classical limit
        }
        HBAR * self.omega * (0.5 * x).tanh().recip() * 0.5
    }

    /// Classical equipartition energy: E_cl = kBT (one mode).
    pub fn classical_energy(&self) -> f64 {
        KB * self.temperature
    }

    /// Quantum correction to classical energy: ΔE = E_boson − E_cl.
    pub fn quantum_correction_energy(&self) -> f64 {
        self.boson_energy() - self.classical_energy()
    }
}

// ---------------------------------------------------------------------------
// PathIntegralEstimator
// ---------------------------------------------------------------------------

/// Centroid virial estimator for kinetic energy from path-integral MD (PIMD).
///
/// Reference: Herman, Bruskin, Berne (1982); Yamamoto (2005).
#[derive(Debug, Clone)]
pub struct PathIntegralEstimator {
    /// Number of ring-polymer beads P.
    pub n_beads: usize,
    /// Particle mass (kg).
    pub mass: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl PathIntegralEstimator {
    /// Create a new path-integral estimator.
    pub fn new(n_beads: usize, mass: f64, temperature: f64) -> Self {
        Self {
            n_beads,
            mass,
            temperature,
        }
    }

    /// Centroid position from bead positions `r` (1-D, length = n_beads).
    pub fn centroid(&self, r: &[f64]) -> f64 {
        r.iter().sum::<f64>() / self.n_beads as f64
    }

    /// Primitive (thermodynamic) kinetic energy estimator (1-D):
    ///
    /// ```text
    /// T_prim = P / (2β) − Σₖ ½ m ω_P² (rₖ₊₁ − rₖ)²
    /// ```
    ///
    /// where `ω_P = √P / (ℏβ)` and the ring is closed (`r_P ≡ r_0`).
    pub fn primitive_kinetic_energy(&self, r: &[f64]) -> f64 {
        let p = self.n_beads as f64;
        let beta = 1.0 / (KB * self.temperature);
        let omega_p = p.sqrt() / (HBAR * beta);
        let spring_const = self.mass * omega_p * omega_p;

        let classical_term = p / (2.0 * beta); // = P kBT / 2
        let mut spring_energy = 0.0;
        for k in 0..self.n_beads {
            let rk = r[k];
            let rk1 = r[(k + 1) % self.n_beads];
            spring_energy += 0.5 * spring_const * (rk1 - rk).powi(2);
        }
        classical_term - spring_energy
    }

    /// Centroid virial kinetic energy estimator (1-D):
    ///
    /// ```text
    /// T_cv = kB T / 2 + (1 / 2P) Σₖ (rₖ − r_c) · ∂V/∂rₖ
    ///      = kB T / 2 − (1 / 2P) Σₖ (rₖ − r_c) · fₖ
    /// ```
    ///
    /// where `r_c` is the centroid and `forces[k]` is the physical force
    /// `fₖ = −∂V/∂rₖ` (the negative gradient of the full potential, not the
    /// ring-polymer spring force) acting on bead `k`.  For a harmonic potential
    /// `V = mω²r²/2` the force is `fₖ = −mω²rₖ`, so the correction is the
    /// non-negative term `(mω²/2P) Σₖ (rₖ − r_c)²` and `T_cv ≥ kB T / 2`, as
    /// required for a quantum kinetic energy.
    ///
    /// Beads without a matching entry in `forces` (and vice versa) are ignored.
    pub fn centroid_virial_kinetic_energy(&self, r: &[f64], forces: &[f64]) -> f64 {
        let rc = self.centroid(r);
        let p = self.n_beads as f64;
        let virial_sum: f64 = r
            .iter()
            .zip(forces.iter())
            .take(self.n_beads)
            .map(|(&rk, &fk)| (rk - rc) * fk)
            .sum();
        0.5 * KB * self.temperature - virial_sum / (2.0 * p)
    }

    /// Quantum kinetic energy estimate from bead spread:
    ///   T_q ≈ P * ℏ² / (4 * m * β) * ⟨(r - r_c)²⟩  (order-of-magnitude).
    pub fn bead_spread_estimator(&self, r: &[f64]) -> f64 {
        let rc = self.centroid(r);
        let mean_sq_dev = r.iter().map(|&ri| (ri - rc).powi(2)).sum::<f64>() / self.n_beads as f64;
        let beta = 1.0 / (KB * self.temperature);
        self.n_beads as f64 * HBAR * HBAR / (4.0 * self.mass * beta) * mean_sq_dev
    }
}

// ---------------------------------------------------------------------------
// TunnellingRate
// ---------------------------------------------------------------------------

/// WKB tunnelling probability and rate through a 1-D potential barrier.
///
/// WKB transmission coefficient:
///   T_WKB = exp(−2 ∫√(2m(V(x)−E)) dx / ℏ)
#[derive(Debug, Clone)]
pub struct TunnellingRate {
    /// Particle mass (kg).
    pub mass: f64,
    /// Barrier height (J).
    pub barrier_height: f64,
    /// Barrier width a (m).
    pub barrier_width: f64,
    /// Incident particle energy E (J).
    pub energy: f64,
}

impl TunnellingRate {
    /// Create a tunnelling rate calculator for a rectangular barrier.
    pub fn new(mass: f64, barrier_height: f64, barrier_width: f64, energy: f64) -> Self {
        Self {
            mass,
            barrier_height,
            barrier_width,
            energy,
        }
    }

    /// WKB exponent κ = √(2m(V−E)) / ℏ  (m⁻¹).
    pub fn wkb_decay_constant(&self) -> f64 {
        let dv = self.barrier_height - self.energy;
        if dv <= 0.0 {
            return 0.0;
        }
        (2.0 * self.mass * dv).sqrt() / HBAR
    }

    /// WKB transmission coefficient for a rectangular barrier.
    pub fn transmission_coefficient(&self) -> f64 {
        let kappa = self.wkb_decay_constant();
        (-2.0 * kappa * self.barrier_width).exp()
    }

    /// Attempt frequency (Boltzmann; classical): ν = kBT / h.
    pub fn attempt_frequency(&self, temperature: f64) -> f64 {
        KB * temperature / H_PLANCK
    }

    /// Tunnelling rate: k = ν × T_WKB  (s⁻¹).
    pub fn tunnelling_rate(&self, temperature: f64) -> f64 {
        self.attempt_frequency(temperature) * self.transmission_coefficient()
    }

    /// Crossover temperature below which tunnelling dominates over classical hopping.
    ///
    /// T_c = ℏ ω_b / (2π kB)  where ω_b = √(2V/m a²).
    pub fn crossover_temperature(&self) -> f64 {
        let omega_b = (2.0 * self.barrier_height / (self.mass * self.barrier_width.powi(2))).sqrt();
        HBAR * omega_b / (2.0 * PI * KB)
    }
}

// ---------------------------------------------------------------------------
// ZeroPointEnergy
// ---------------------------------------------------------------------------

/// Zero-point energy calculator for harmonic and anharmonically-corrected oscillators.
#[derive(Debug, Clone)]
pub struct ZeroPointEnergy {
    /// Angular frequency ω (rad s⁻¹).
    pub omega: f64,
    /// Particle mass (kg).
    pub mass: f64,
}

impl ZeroPointEnergy {
    /// Create a ZPE calculator.
    pub fn new(omega: f64, mass: f64) -> Self {
        Self { omega, mass }
    }

    /// Harmonic zero-point energy: E₀ = ℏω/2  (J).
    pub fn harmonic_zpe(&self) -> f64 {
        0.5 * HBAR * self.omega
    }

    /// Convert ZPE to wavenumbers (cm⁻¹): ν̃ = ω / (2π c) in cm⁻¹.
    pub fn wavenumber_zpe(&self) -> f64 {
        let c_cm = 2.997_924_58e10; // speed of light in cm/s
        self.omega / (2.0 * PI * c_cm)
    }

    /// Anharmonic ZPE correction (Morse oscillator):
    ///   E₀ = ℏω/2 · (1 − x_e/2)
    ///
    /// where `anharmonicity` xₑ = ℏω/(4Dₑ) and Dₑ is the well depth (J).
    pub fn morse_zpe(&self, well_depth: f64) -> f64 {
        let x_e = HBAR * self.omega / (4.0 * well_depth);
        self.harmonic_zpe() * (1.0 - x_e * 0.5)
    }

    /// ZPE for a 3-D harmonic oscillator with three distinct frequencies.
    pub fn zpe_3d(omega_x: f64, omega_y: f64, omega_z: f64) -> f64 {
        0.5 * HBAR * (omega_x + omega_y + omega_z)
    }

    /// Relative ZPE change on mass substitution m → m' (harmonic approximation).
    ///
    /// ZPE ∝ 1/√m, so ΔE/E = (√(m/m') − 1).
    pub fn mass_substitution_zpe_change(&self, new_mass: f64) -> f64 {
        self.harmonic_zpe() * ((self.mass / new_mass).sqrt() - 1.0)
    }
}

// ---------------------------------------------------------------------------
// IsotopeEffect
// ---------------------------------------------------------------------------

/// Kinetic isotope effect (KIE) from mass-dependent ZPE differences.
///
/// KIE = exp(ΔZPE_reactant / kBT) for the primary KIE.
#[derive(Debug, Clone)]
pub struct IsotopeEffect {
    /// Angular frequency of the bond in the reactant (rad s⁻¹).
    pub omega_reactant: f64,
    /// Angular frequency of the bond in the transition state (rad s⁻¹).
    pub omega_ts: f64,
    /// Light isotope mass (kg), e.g. proton.
    pub mass_light: f64,
    /// Heavy isotope mass (kg), e.g. deuteron.
    pub mass_heavy: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl IsotopeEffect {
    /// Create an isotope effect calculator (default: H/D).
    pub fn new_hd(omega_reactant: f64, omega_ts: f64, temperature: f64) -> Self {
        Self {
            omega_reactant,
            omega_ts,
            mass_light: M_PROTON,
            mass_heavy: M_DEUTERON,
            temperature,
        }
    }

    /// Create with explicit masses.
    pub fn new(
        omega_reactant: f64,
        omega_ts: f64,
        mass_light: f64,
        mass_heavy: f64,
        temperature: f64,
    ) -> Self {
        Self {
            omega_reactant,
            omega_ts,
            mass_light,
            mass_heavy,
            temperature,
        }
    }

    /// ZPE difference: ΔZPE = ZPE_H − ZPE_D (reactant).
    pub fn delta_zpe_reactant(&self) -> f64 {
        let zpe_light = 0.5 * HBAR * self.omega_reactant;
        let zpe_heavy =
            0.5 * HBAR * self.omega_reactant * (self.mass_light / self.mass_heavy).sqrt();
        zpe_light - zpe_heavy
    }

    /// ZPE difference at the transition state.
    pub fn delta_zpe_ts(&self) -> f64 {
        let zpe_light = 0.5 * HBAR * self.omega_ts;
        let zpe_heavy = 0.5 * HBAR * self.omega_ts * (self.mass_light / self.mass_heavy).sqrt();
        zpe_light - zpe_heavy
    }

    /// Primary kinetic isotope effect: KIE = exp((ΔZPE_reactant − ΔZPE_ts) / kBT).
    pub fn primary_kie(&self) -> f64 {
        let delta = self.delta_zpe_reactant() - self.delta_zpe_ts();
        (delta / (KB * self.temperature)).exp()
    }

    /// Swain-Schaad exponent relating H/D and H/T KIEs: exponent ≈ 1.44 (semiclassical).
    pub fn swain_schaad_exponent() -> f64 {
        // log(KIE_HT) / log(KIE_HD) ≈ 1.44 semiclassically
        let r = (M_PROTON / (3.0 * M_PROTON)).sqrt(); // H/T mass ratio proxy
        (1.0 - r).ln() / (1.0 - (M_PROTON / M_DEUTERON).sqrt()).ln()
    }
}

// ---------------------------------------------------------------------------
// QuantumCorrection (Feynman-Hibbs)
// ---------------------------------------------------------------------------

/// Feynman-Hibbs effective potential quantum correction.
///
/// V_eff(r) = V(r) + ℏ²/(24mkBT) · ∇²V(r)
///
/// Valid at high temperature where quantum effects are small.
#[derive(Debug, Clone)]
pub struct QuantumCorrection {
    /// Particle mass (kg).
    pub mass: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl QuantumCorrection {
    /// Create a Feynman-Hibbs corrector.
    pub fn new(mass: f64, temperature: f64) -> Self {
        Self { mass, temperature }
    }

    /// Feynman-Hibbs prefactor: α = ℏ²/(24mkBT).
    pub fn fh_prefactor(&self) -> f64 {
        HBAR * HBAR / (24.0 * self.mass * KB * self.temperature)
    }

    /// Feynman-Hibbs effective potential: V_eff = V + α ∇²V.
    pub fn effective_potential(&self, v: f64, laplacian_v: f64) -> f64 {
        v + self.fh_prefactor() * laplacian_v
    }

    /// Feynman-Hibbs effective force: F_eff = F − α ∇(∇²V).
    ///
    /// `grad_laplacian_v` = ∂(∇²V)/∂r.
    pub fn effective_force(&self, force: f64, grad_laplacian_v: f64) -> f64 {
        force - self.fh_prefactor() * grad_laplacian_v
    }

    /// Correction energy for a LJ pair at distance r.
    ///
    /// For V_LJ = 4ε((σ/r)¹² − (σ/r)⁶), the Laplacian in 3-D is:
    ///   ∇²V = d²V/dr² + (2/r) dV/dr
    pub fn lj_correction(&self, r: f64, epsilon_lj: f64, sigma_lj: f64) -> f64 {
        let sr = sigma_lj / r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        // dV/dr = 4ε(-12σ¹²/r¹³ + 6σ⁶/r⁷)
        let dvdr = 4.0 * epsilon_lj * (-12.0 * sr12 / r + 6.0 * sr6 / r);
        // d²V/dr² = 4ε(12·13σ¹²/r¹⁴ - 6·7σ⁶/r⁸)
        let d2vdr2 =
            4.0 * epsilon_lj * (12.0 * 13.0 * sr12 / r.powi(2) - 6.0 * 7.0 * sr6 / r.powi(2));
        let laplacian_v = d2vdr2 + 2.0 / r * dvdr;
        self.fh_prefactor() * laplacian_v
    }
}

// ---------------------------------------------------------------------------
// DebyeModel
// ---------------------------------------------------------------------------

/// Debye heat capacity model.
///
/// Cv(T) = 9 N kB (T/θD)³ ∫₀^(θD/T) x⁴eˣ/(eˣ−1)² dx
#[derive(Debug, Clone)]
pub struct DebyeModel {
    /// Debye temperature θ_D (K).
    pub debye_temperature: f64,
    /// Number of atoms N.
    pub n_atoms: usize,
}

impl DebyeModel {
    /// Create a Debye model.
    pub fn new(debye_temperature: f64, n_atoms: usize) -> Self {
        Self {
            debye_temperature,
            n_atoms,
        }
    }

    /// Debye *energy* function `D₃(x) = (3/x³) ∫₀ˣ t³/(eᵗ − 1) dt` at `x = θ_D/T`.
    ///
    /// `D₃` gives the thermal part of the internal energy,
    /// `U_th = 3 N k_B T · D₃(θ_D/T)`.  It is **not** the heat-capacity function:
    /// `D₃(1) = 0.674` whereas `C_V(θ_D)/(3 N k_B) = 0.952`.
    ///
    /// Integrating the heat-capacity integral by parts gives the exact identity
    ///
    /// ```text
    /// C_V / (3 N k_B) = 4 D₃(x) − 3x / (eˣ − 1),
    /// ```
    ///
    /// so `D₃ = [c(x) + 3x/(eˣ − 1)] / 4`, where `c(x)` is the accurate (~1e-13)
    /// heat-capacity function [`crate::solid_state_md::debye_cv_function`].
    /// Both terms are non-negative, so no cancellation occurs at any `x`.
    /// `D₃(0) = 1` (classical limit) and `D₃(x) → π⁴/(5x³)` for `x → ∞`.
    fn debye_integral(x_d: f64) -> f64 {
        let cv_function = debye_cv_function(x_d);
        // x / (eˣ − 1) → 1 for x → 0 and → 0 for x → ∞.
        let bose_term = if x_d == 0.0 { 1.0 } else { x_d / x_d.exp_m1() };
        0.25 * (cv_function + 3.0 * bose_term)
    }

    /// Heat capacity at constant volume `C_V(T)` (J K⁻¹).
    ///
    /// ```text
    /// C_V = 9 N k_B (T/θ_D)³ ∫₀^{θ_D/T} x⁴ eˣ / (eˣ − 1)² dx
    ///     = 3 N k_B [4 D₃(θ_D/T) − 3 (θ_D/T) / (e^{θ_D/T} − 1)]
    /// ```
    ///
    /// It follows the Debye `T³` law [`DebyeModel::low_temperature_cv`] for
    /// `T ≪ θ_D`, equals `2.855 N k_B` at `T = θ_D` and approaches the
    /// Dulong-Petit limit [`DebyeModel::dulong_petit_limit`] from below for
    /// `T ≫ θ_D`.  The integral is evaluated by the shared, overflow-free
    /// quadrature in [`crate::solid_state_md::debye_cv_function`].
    ///
    /// Returns `0.0` for `T ≤ 0` (all modes frozen out) and `NaN` for a `NaN`
    /// temperature.
    pub fn heat_capacity(&self, temperature: f64) -> f64 {
        if temperature.is_nan() {
            return f64::NAN;
        }
        if temperature <= 0.0 {
            return 0.0;
        }
        3.0 * self.n_atoms as f64 * KB * debye_cv_function(self.debye_temperature / temperature)
    }

    /// Thermal vibrational energy above the zero-point level (J):
    /// `U_th(T) = 3 N k_B T · D₃(θ_D/T)`.
    ///
    /// It tends to the classical `3 N k_B T` for `T ≫ θ_D` and to
    /// `(3π⁴/5) N k_B T⁴/θ_D³` for `T ≪ θ_D`; its temperature derivative is
    /// [`DebyeModel::heat_capacity`].  Returns `0.0` for `T ≤ 0` and `NaN` for a
    /// `NaN` temperature.
    pub fn thermal_energy(&self, temperature: f64) -> f64 {
        if temperature.is_nan() {
            return f64::NAN;
        }
        if temperature <= 0.0 {
            return 0.0;
        }
        3.0 * self.n_atoms as f64
            * KB
            * temperature
            * Self::debye_integral(self.debye_temperature / temperature)
    }

    /// Total vibrational energy `U(T) = E₀ + U_th(T)` (J), the sum of the
    /// zero-point energy [`DebyeModel::zero_point_energy`] and the thermal energy
    /// [`DebyeModel::thermal_energy`].
    pub fn internal_energy(&self, temperature: f64) -> f64 {
        self.zero_point_energy() + self.thermal_energy(temperature)
    }

    /// High-temperature (Dulong-Petit) limit: Cv → 3 N kB.
    pub fn dulong_petit_limit(&self) -> f64 {
        3.0 * self.n_atoms as f64 * KB
    }

    /// Debye T³ low-temperature heat capacity: Cv ≈ (12π⁴/5) N kB (T/θD)³.
    pub fn low_temperature_cv(&self, temperature: f64) -> f64 {
        let ratio = temperature / self.debye_temperature;
        12.0 * PI.powi(4) / 5.0 * self.n_atoms as f64 * KB * ratio.powi(3)
    }

    /// Total zero-point energy of the `N` atoms: `E₀ = (9/8) N k_B θ_D`.
    pub fn zero_point_energy(&self) -> f64 {
        9.0 / 8.0 * KB * self.debye_temperature * self.n_atoms as f64
    }
}

// ---------------------------------------------------------------------------
// EinsteinModel
// ---------------------------------------------------------------------------

/// Einstein heat capacity model.
///
/// Cv(T) = 3 N kB (x_E)² e^x_E / (e^x_E − 1)²
/// where x_E = θ_E / T = ℏω_E / (kBT).
#[derive(Debug, Clone)]
pub struct EinsteinModel {
    /// Einstein temperature θ_E (K).
    pub einstein_temperature: f64,
    /// Number of atoms.
    pub n_atoms: usize,
}

impl EinsteinModel {
    /// Create an Einstein model.
    pub fn new(einstein_temperature: f64, n_atoms: usize) -> Self {
        Self {
            einstein_temperature,
            n_atoms,
        }
    }

    /// Dimensionless parameter x_E = θ_E / T.
    pub fn x_param(&self, temperature: f64) -> f64 {
        self.einstein_temperature / temperature
    }

    /// Heat capacity at constant volume (J K⁻¹).
    pub fn heat_capacity(&self, temperature: f64) -> f64 {
        let x = self.x_param(temperature);
        let ex = x.exp();
        let denom = ex - 1.0;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        3.0 * self.n_atoms as f64 * KB * x * x * ex / (denom * denom)
    }

    /// Zero-point energy: E_0 = (3/2) N kB θ_E = (3/2) N ℏω_E.
    pub fn zero_point_energy(&self) -> f64 {
        1.5 * self.n_atoms as f64 * KB * self.einstein_temperature
    }

    /// High-temperature limit → Dulong-Petit: 3 N kB.
    pub fn high_temperature_limit(&self) -> f64 {
        3.0 * self.n_atoms as f64 * KB
    }

    /// Mean total energy: E = 3 N kB θ_E (n̄ + 1/2) (3 modes per atom).
    ///
    /// Consistent with the 3-D Einstein model and `zero_point_energy`.
    pub fn mean_energy(&self, temperature: f64) -> f64 {
        let x = self.x_param(temperature);
        let n_bar = if x > 700.0 {
            0.0
        } else {
            1.0 / (x.exp() - 1.0)
        };
        3.0 * self.n_atoms as f64 * KB * self.einstein_temperature * (n_bar + 0.5)
    }
}

// ---------------------------------------------------------------------------
// NuclearQuantumEffect
// ---------------------------------------------------------------------------

/// Nuclear quantum effects: proton tunnelling and H/D isotope fractionation.
#[derive(Debug, Clone)]
pub struct NuclearQuantumEffect {
    /// Angular frequency of the O-H (or N-H) stretch (rad s⁻¹).
    pub omega_oh: f64,
    /// Barrier height for proton transfer (J).
    pub barrier_height: f64,
    /// Transfer distance (m).
    pub transfer_distance: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl NuclearQuantumEffect {
    /// Create a nuclear quantum effect calculator.
    pub fn new(
        omega_oh: f64,
        barrier_height: f64,
        transfer_distance: f64,
        temperature: f64,
    ) -> Self {
        Self {
            omega_oh,
            barrier_height,
            transfer_distance,
            temperature,
        }
    }

    /// Proton tunnelling WKB transmission coefficient.
    pub fn proton_tunnelling_coefficient(&self) -> f64 {
        let kappa = (2.0 * M_PROTON * self.barrier_height).sqrt() / HBAR;
        (-2.0 * kappa * self.transfer_distance).exp()
    }

    /// Deuteron tunnelling coefficient (heavier → less tunnelling).
    pub fn deuteron_tunnelling_coefficient(&self) -> f64 {
        let kappa = (2.0 * M_DEUTERON * self.barrier_height).sqrt() / HBAR;
        (-2.0 * kappa * self.transfer_distance).exp()
    }

    /// Tunnelling enhancement factor relative to classical rate.
    ///
    /// Uses Bell correction: Qt = (u/2) / sin(u/2)  where u = ℏω_b / (kBT).
    pub fn bell_tunnelling_correction(&self) -> f64 {
        let omega_b =
            (2.0 * self.barrier_height / (M_PROTON * self.transfer_distance.powi(2))).sqrt();
        let u = HBAR * omega_b / (KB * self.temperature);
        if u.abs() < 1e-6 {
            return 1.0;
        }
        let half_u = u * 0.5;
        half_u / half_u.sin().max(1e-30)
    }

    /// H/D isotope fractionation factor:
    ///   α_H/D = (T_H / T_D) × exp((ZPE_D − ZPE_H) / kBT)
    pub fn hd_fractionation_factor(&self) -> f64 {
        let t_h = self.proton_tunnelling_coefficient();
        let t_d = self.deuteron_tunnelling_coefficient();
        let zpe_h = 0.5 * HBAR * self.omega_oh;
        let zpe_d = 0.5 * HBAR * self.omega_oh * (M_PROTON / M_DEUTERON).sqrt();
        let zpe_factor = ((zpe_d - zpe_h) / (KB * self.temperature)).exp();
        if t_d.abs() < 1e-300 {
            return f64::INFINITY;
        }
        (t_h / t_d) * zpe_factor
    }

    /// Quantum delocalization (bead spread) for a proton at temperature T:
    ///   Δr ≈ √(ℏ / (m ω)) = thermal de Broglie / √(2π)  (order of magnitude).
    pub fn quantum_delocalization(&self) -> f64 {
        (HBAR / (M_PROTON * self.omega_oh)).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- WignerKirkwoodExpansion ----

    #[test]
    fn test_wk_thermal_de_broglie_positive() {
        let wk = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        assert!(wk.thermal_de_broglie() > 0.0);
    }

    #[test]
    fn test_wk_prefactor_decreases_with_temperature() {
        let wk1 = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        let wk2 = WignerKirkwoodExpansion::new(M_PROTON, 600.0);
        assert!(wk1.wk_prefactor() > wk2.wk_prefactor());
    }

    #[test]
    fn test_wk_prefactor_decreases_with_mass() {
        let wk_h = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        let wk_d = WignerKirkwoodExpansion::new(M_DEUTERON, 300.0);
        assert!(wk_h.wk_prefactor() > wk_d.wk_prefactor());
    }

    #[test]
    fn test_wk_free_energy_correction_zero_laplacian() {
        let wk = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        assert_eq!(wk.free_energy_correction(0.0), 0.0);
    }

    #[test]
    fn test_wk_partition_ratio_near_unity() {
        let wk = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        // Very small Laplacian → ratio ≈ 1
        let ratio = wk.partition_function_ratio(1e-20);
        assert!((ratio - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_wk_internal_energy_correction() {
        let wk = WignerKirkwoodExpansion::new(M_PROTON, 300.0);
        // Zero gradient and zero Laplacian → correction is exactly zero
        let du = wk.internal_energy_correction(0.0, 0.0);
        assert_eq!(du, 0.0);
    }

    // ---- QuantumThermalEnergy ----

    #[test]
    fn test_qte_zero_point_energy_positive() {
        let q = QuantumThermalEnergy::new(1e13, 300.0);
        assert!(q.zero_point_energy() > 0.0);
    }

    #[test]
    fn test_qte_boson_energy_ge_zpe() {
        let q = QuantumThermalEnergy::new(1e13, 300.0);
        assert!(q.boson_energy() >= q.zero_point_energy() - 1e-40);
    }

    #[test]
    fn test_qte_high_temperature_classical_limit() {
        // At very high T, boson energy → kBT (classical)
        let omega = 1e10; // very low frequency
        let t = 10000.0;
        let q = QuantumThermalEnergy::new(omega, t);
        let ratio = q.boson_energy() / q.classical_energy();
        assert!((ratio - 1.0).abs() < 0.05);
    }

    #[test]
    fn test_qte_bose_einstein_occupation_positive() {
        let q = QuantumThermalEnergy::new(1e13, 300.0);
        assert!(q.bose_einstein_occupation() >= 0.0);
    }

    #[test]
    fn test_qte_fermi_dirac_between_zero_and_one() {
        let q = QuantumThermalEnergy::new(1e13, 300.0);
        let n = q.fermi_dirac_occupation();
        assert!(n > 0.0 && n < 1.0);
    }

    #[test]
    fn test_qte_x_param_scaling() {
        let q1 = QuantumThermalEnergy::new(1e13, 300.0);
        let q2 = QuantumThermalEnergy::new(1e13, 600.0);
        assert!((q1.x_param() / q2.x_param() - 2.0).abs() < 1e-10);
    }

    // ---- PathIntegralEstimator ----

    #[test]
    fn test_pie_centroid_uniform() {
        let pie = PathIntegralEstimator::new(16, M_PROTON, 300.0);
        let r = vec![1.0; 16];
        assert!((pie.centroid(&r) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_pie_centroid_linear() {
        let pie = PathIntegralEstimator::new(4, M_PROTON, 300.0);
        let r = vec![0.0, 1.0, 2.0, 3.0];
        assert!((pie.centroid(&r) - 1.5).abs() < 1e-14);
    }

    #[test]
    fn test_pie_primitive_ke_positive_low_spread() {
        let pie = PathIntegralEstimator::new(8, M_PROTON, 300.0);
        // Beads near centroid → spring energy ≈ 0 → T_prim ≈ P kBT/2
        let r = vec![1.0; 8];
        let t_prim = pie.primitive_kinetic_energy(&r);
        let classical = 4.0 * KB * 300.0; // P kBT / 2 = 8 * kBT / 2
        assert!((t_prim - classical).abs() / classical < 0.01);
    }

    #[test]
    fn test_pie_virial_estimator_uniform_force() {
        let pie = PathIntegralEstimator::new(4, M_PROTON, 300.0);
        let r = vec![0.0; 4]; // all at centroid
        let forces = vec![1e-10; 4];
        let t_vir = pie.centroid_virial_kinetic_energy(&r, &forces);
        // With r_k = r_c, virial term = 0 → T_vir = kBT/2
        assert!((t_vir - 0.5 * KB * 300.0).abs() < 1e-30);
    }

    #[test]
    fn test_pie_bead_spread_estimator_zero_for_uniform() {
        let pie = PathIntegralEstimator::new(8, M_PROTON, 300.0);
        let r = vec![0.5; 8]; // all at same position
        let spread = pie.bead_spread_estimator(&r);
        assert!(spread.abs() < 1e-50);
    }

    /// Pre-fix centroid-virial formula `kBT/2 + (1/2P) Σ (rₖ − r_c)·fₖ`, kept
    /// here only to pin down the sign convention of the corrected estimator.
    fn legacy_virial_with_plus_sign(pie: &PathIntegralEstimator, r: &[f64], f: &[f64]) -> f64 {
        let rc = pie.centroid(r);
        let sum: f64 = r.iter().zip(f).map(|(&rk, &fk)| (rk - rc) * fk).sum();
        0.5 * KB * pie.temperature + sum / (2.0 * pie.n_beads as f64)
    }

    #[test]
    fn test_pie_virial_estimator_harmonic_beads_exceeds_classical() {
        // Spread-out beads in V = m w² r²/2 with the physical force f_k = −m w² r_k:
        // T_cv = kBT/2 + (m w²/2P) Σ (r_k − r_c)², strictly above kBT/2.
        let mass = M_PROTON;
        let omega = 3.0e14;
        let temperature = 300.0;
        let pie = PathIntegralEstimator::new(4, mass, temperature);
        let r = [1.0e-11, 3.0e-11, -2.0e-11, 0.5e-11];
        let forces: Vec<f64> = r.iter().map(|&x| -mass * omega * omega * x).collect();
        let t_cv = pie.centroid_virial_kinetic_energy(&r, &forces);
        let rc = pie.centroid(&r);
        let spread: f64 = r.iter().map(|&x| (x - rc).powi(2)).sum();
        let expected = 0.5 * KB * temperature + mass * omega * omega / (2.0 * 4.0) * spread;
        assert!(
            t_cv > 0.5 * KB * temperature,
            "T_cv = {t_cv} must exceed kBT/2 = {}",
            0.5 * KB * temperature
        );
        assert!(
            (t_cv - expected).abs() <= 1e-12 * expected.abs(),
            "T_cv = {t_cv}, expected {expected}"
        );
    }

    #[test]
    fn test_pie_virial_estimator_sign_convention_is_minus_force() {
        // The corrected estimator is the old '+' formula evaluated with the force
        // sign flipped, and a gradient (+dV/dr) input reproduces the old value.
        let mass = M_PROTON;
        let omega = 2.0e14;
        let pie = PathIntegralEstimator::new(6, mass, 250.0);
        let r = [2.0e-11, -1.0e-11, 0.4e-11, 1.5e-11, -2.2e-11, 0.9e-11];
        let force: Vec<f64> = r.iter().map(|&x| -mass * omega * omega * x).collect();
        let neg_force: Vec<f64> = force.iter().map(|&f| -f).collect();
        let corrected = pie.centroid_virial_kinetic_energy(&r, &force);
        let legacy_on_flipped = legacy_virial_with_plus_sign(&pie, &r, &neg_force);
        assert!(
            (corrected - legacy_on_flipped).abs() <= 1e-12 * corrected.abs(),
            "corrected = {corrected}, legacy(−f) = {legacy_on_flipped}"
        );
        // Feeding +dV/dr (instead of the force) lands on the old, wrong value
        // below the classical kinetic energy.
        let wrong = pie.centroid_virial_kinetic_energy(&r, &neg_force);
        let legacy_on_force = legacy_virial_with_plus_sign(&pie, &r, &force);
        assert!((wrong - legacy_on_force).abs() <= 1e-12 * wrong.abs());
        assert!(wrong < 0.5 * KB * 250.0);
    }

    /// Bead positions of a 1-D ring polymer in `V = m w² r²/2` whose real normal-mode
    /// amplitudes all equal their canonical rms value.
    ///
    /// The primitive and centroid-virial estimators are quadratic forms that are
    /// diagonal in the normal modes, so this configuration returns exactly their
    /// canonical average.  The mode variance is
    /// `1/(β m [ω_P² 4 sin²(π j/P) + ω²/P])` with `ω_P² = P/(βħ)²`.
    fn harmonic_ring_polymer_snapshot(
        p: usize,
        mass: f64,
        omega: f64,
        temperature: f64,
    ) -> Vec<f64> {
        let beta = 1.0 / (KB * temperature);
        let omega_p_sq = p as f64 / (beta * HBAR).powi(2);
        let variance = |j: usize| {
            let s = (PI * j as f64 / p as f64).sin();
            1.0 / (beta * mass * (omega_p_sq * 4.0 * s * s + omega * omega / p as f64))
        };
        let norm = (2.0 / p as f64).sqrt();
        (0..p)
            .map(|k| {
                let mut x = 0.0;
                for j in 1..p / 2 {
                    let angle = 2.0 * PI * (j * k) as f64 / p as f64;
                    x += norm * variance(j).sqrt() * (angle.cos() + angle.sin());
                }
                let alternating = if k % 2 == 0 { 1.0 } else { -1.0 };
                x + alternating * variance(p / 2).sqrt() / (p as f64).sqrt()
            })
            .collect()
    }

    #[test]
    fn test_pie_estimators_reproduce_harmonic_oscillator_kinetic_energy() {
        // beta*hbar*omega = 4: <T> = (hbar w / 4) coth(beta hbar w / 2) per degree
        // of freedom, i.e. 1.0373 * hbar w / 4 against the classical kBT/2 = hbar w / 8.
        let mass = M_PROTON;
        let omega = 2.0e14;
        let temperature = HBAR * omega / (4.0 * KB);
        let p = 32;
        let pie = PathIntegralEstimator::new(p, mass, temperature);
        let r = harmonic_ring_polymer_snapshot(p, mass, omega, temperature);
        let forces: Vec<f64> = r.iter().map(|&x| -mass * omega * omega * x).collect();
        let t_prim = pie.primitive_kinetic_energy(&r);
        let t_cv = pie.centroid_virial_kinetic_energy(&r, &forces);
        let exact = 0.25 * HBAR * omega / (0.5 * HBAR * omega / (KB * temperature)).tanh();
        assert!(
            (t_cv - exact).abs() / exact < 0.01,
            "centroid virial {t_cv} vs exact {exact}"
        );
        assert!(
            (t_prim - exact).abs() / exact < 0.01,
            "primitive {t_prim} vs exact {exact}"
        );
        // Both are unbiased estimators of the same finite-P canonical average.
        assert!(
            (t_cv - t_prim).abs() <= 1e-9 * exact,
            "virial {t_cv} and primitive {t_prim} must coincide on average"
        );
        // The quantum kinetic energy is well above the classical value.
        assert!(t_cv > 1.9 * 0.5 * KB * temperature);
    }

    // ---- TunnellingRate ----

    #[test]
    fn test_tunnelling_wkb_decay_positive() {
        let tr = TunnellingRate::new(M_PROTON, 5e-21, 1e-10, 1e-21);
        assert!(tr.wkb_decay_constant() > 0.0);
    }

    #[test]
    fn test_tunnelling_coefficient_between_zero_one() {
        let tr = TunnellingRate::new(M_PROTON, 5e-21, 1e-10, 1e-21);
        let t = tr.transmission_coefficient();
        assert!((0.0..=1.0).contains(&t));
    }

    #[test]
    fn test_tunnelling_coefficient_zero_above_barrier() {
        // E >= V → kappa = 0 → T = 1
        let tr = TunnellingRate::new(M_PROTON, 5e-21, 1e-10, 6e-21);
        assert_eq!(tr.wkb_decay_constant(), 0.0);
        assert_eq!(tr.transmission_coefficient(), 1.0);
    }

    #[test]
    fn test_tunnelling_rate_positive() {
        let tr = TunnellingRate::new(M_PROTON, 5e-21, 5e-11, 1e-21);
        assert!(tr.tunnelling_rate(300.0) > 0.0);
    }

    #[test]
    fn test_tunnelling_crossover_temperature_positive() {
        let tr = TunnellingRate::new(M_PROTON, 5e-21, 1e-10, 1e-21);
        assert!(tr.crossover_temperature() > 0.0);
    }

    #[test]
    fn test_tunnelling_heavier_particle_lower_rate() {
        let tr_h = TunnellingRate::new(M_PROTON, 5e-21, 1e-10, 1e-21);
        let tr_d = TunnellingRate::new(M_DEUTERON, 5e-21, 1e-10, 1e-21);
        assert!(tr_h.transmission_coefficient() > tr_d.transmission_coefficient());
    }

    // ---- ZeroPointEnergy ----

    #[test]
    fn test_zpe_harmonic_positive() {
        let zpe = ZeroPointEnergy::new(1e14, M_PROTON);
        assert!(zpe.harmonic_zpe() > 0.0);
    }

    #[test]
    fn test_zpe_wavenumber_positive() {
        let zpe = ZeroPointEnergy::new(1e14, M_PROTON);
        assert!(zpe.wavenumber_zpe() > 0.0);
    }

    #[test]
    fn test_zpe_morse_less_than_harmonic() {
        let zpe = ZeroPointEnergy::new(1e14, M_PROTON);
        let well_depth = 1e-18; // large well → small anharmonicity
        let morse = zpe.morse_zpe(well_depth);
        assert!(morse < zpe.harmonic_zpe());
    }

    #[test]
    fn test_zpe_3d_sum() {
        let zpe = ZeroPointEnergy::zpe_3d(1e13, 2e13, 3e13);
        assert!((zpe - 0.5 * HBAR * 6e13).abs() < 1e-40);
    }

    #[test]
    fn test_zpe_mass_substitution_heavier_lower() {
        let zpe = ZeroPointEnergy::new(1e14, M_PROTON);
        let delta = zpe.mass_substitution_zpe_change(M_DEUTERON);
        // ZPE ∝ 1/√m; heavier D → lower ZPE, so delta = ZPE*(√(m_H/m_D)-1) < 0
        assert!(delta < 0.0);
    }

    // ---- IsotopeEffect ----

    #[test]
    fn test_kie_delta_zpe_positive() {
        let ie = IsotopeEffect::new_hd(1e14, 1e14, 300.0);
        assert!(ie.delta_zpe_reactant() > 0.0);
    }

    #[test]
    fn test_kie_primary_greater_than_one() {
        // Standard H/D KIE for a normal reaction should be > 1
        let ie = IsotopeEffect::new_hd(3.5e14, 2.0e14, 298.0);
        let kie = ie.primary_kie();
        assert!(kie > 1.0);
    }

    #[test]
    fn test_kie_increases_at_lower_temperature() {
        let ie1 = IsotopeEffect::new_hd(3.5e14, 2.0e14, 300.0);
        let ie2 = IsotopeEffect::new_hd(3.5e14, 2.0e14, 150.0);
        assert!(ie2.primary_kie() > ie1.primary_kie());
    }

    // ---- QuantumCorrection ----

    #[test]
    fn test_fh_prefactor_positive() {
        let qc = QuantumCorrection::new(M_PROTON, 300.0);
        assert!(qc.fh_prefactor() > 0.0);
    }

    #[test]
    fn test_fh_effective_potential_adds_correction() {
        let qc = QuantumCorrection::new(M_PROTON, 300.0);
        let v = 1.0e-20;
        let lap = 1.0e30;
        let v_eff = qc.effective_potential(v, lap);
        assert!(v_eff > v);
    }

    #[test]
    fn test_fh_effective_force() {
        let qc = QuantumCorrection::new(M_PROTON, 300.0);
        let f = -1e-10;
        let glap = 1e40;
        let f_eff = qc.effective_force(f, glap);
        // Gradient correction reduces magnitude
        assert!(f_eff.abs() != f.abs());
    }

    #[test]
    fn test_fh_lj_correction_finite() {
        let qc = QuantumCorrection::new(M_PROTON, 300.0);
        let corr = qc.lj_correction(3.5e-10, 1e-21, 3.4e-10);
        assert!(corr.is_finite());
    }

    // ---- DebyeModel ----

    #[test]
    fn test_debye_high_temperature_dulong_petit() {
        let dm = DebyeModel::new(300.0, 100);
        let cv_high = dm.heat_capacity(10000.0);
        let dp = dm.dulong_petit_limit();
        assert!((cv_high / dp - 1.0).abs() < 0.02);
    }

    #[test]
    fn test_debye_low_temperature_t3() {
        let dm = DebyeModel::new(300.0, 100);
        let cv1 = dm.low_temperature_cv(10.0);
        let cv2 = dm.low_temperature_cv(20.0);
        // Cv ∝ T³ → cv2 = 8 × cv1
        assert!((cv2 / cv1 - 8.0).abs() < 0.01);
    }

    #[test]
    fn test_debye_zero_point_energy_positive() {
        let dm = DebyeModel::new(400.0, 10);
        assert!(dm.zero_point_energy() > 0.0);
    }

    #[test]
    fn test_debye_cv_monotone_increasing() {
        let dm = DebyeModel::new(300.0, 50);
        let cv1 = dm.heat_capacity(50.0);
        let cv2 = dm.heat_capacity(100.0);
        let cv3 = dm.heat_capacity(500.0);
        assert!(cv1 < cv2 && cv2 < cv3);
    }

    #[test]
    fn test_debye_cv_at_debye_temperature_is_2p855_nk() {
        // C_V(T = theta_D) = 3 N kB [4 D3(1) − 3/(e − 1)] = 2.8551964 N kB.
        // The energy function D3 alone would give 3 D3(1) = 2.0232 N kB.
        let n_atoms = 100;
        let dm = DebyeModel::new(300.0, n_atoms);
        let nk = n_atoms as f64 * KB;
        let cv = dm.heat_capacity(300.0);
        assert!(
            (cv / nk - 2.855_196_407).abs() < 1e-8,
            "C_V(theta_D) = {} N kB",
            cv / nk
        );
        let energy_function_value = 3.0 * 0.674_415_564_078 * nk;
        assert!((cv - energy_function_value).abs() > 0.5 * nk);
    }

    #[test]
    fn test_debye_cv_follows_t3_law_well_below_debye_temperature() {
        let dm = DebyeModel::new(400.0, 25);
        for fraction in [50.0, 100.0, 200.0] {
            let t = dm.debye_temperature / fraction;
            let ratio = dm.heat_capacity(t) / dm.low_temperature_cv(t);
            assert!(
                (ratio - 1.0).abs() < 1e-6,
                "C_V / C_T3 at theta/{fraction} = {ratio}"
            );
        }
        // At theta/10 the T³ law is still within a few percent, but no better.
        let t = dm.debye_temperature / 10.0;
        let ratio = dm.heat_capacity(t) / dm.low_temperature_cv(t);
        assert!((ratio - 0.972_971_339_6).abs() < 1e-8, "ratio = {ratio}");
    }

    #[test]
    fn test_debye_cv_approaches_dulong_petit_from_below() {
        let dm = DebyeModel::new(300.0, 64);
        let dulong_petit = dm.dulong_petit_limit();
        let cv_10_theta = dm.heat_capacity(3000.0);
        // C_V/3NkB = 1 − x²/20 + ... with x = 0.1.
        assert!(cv_10_theta < dulong_petit);
        assert!(
            (cv_10_theta / dulong_petit - 0.999_500_178_5).abs() < 1e-9,
            "C_V(10 theta)/3NkB = {}",
            cv_10_theta / dulong_petit
        );
        let mut previous = 0.0;
        for t in [
            5.0, 20.0, 60.0, 150.0, 300.0, 600.0, 1500.0, 6000.0, 60000.0,
        ] {
            let cv = dm.heat_capacity(t);
            assert!(cv > previous, "C_V must increase with T (T = {t})");
            assert!(cv <= dulong_petit * (1.0 + 1e-12));
            previous = cv;
        }
    }

    #[test]
    fn test_debye_cv_degenerate_temperatures() {
        let dm = DebyeModel::new(300.0, 10);
        assert_eq!(dm.heat_capacity(0.0), 0.0);
        assert_eq!(dm.heat_capacity(-12.0), 0.0);
        assert!(dm.heat_capacity(f64::NAN).is_nan());
        // A vanishing Debye temperature means no freeze-out at all.
        let classical = DebyeModel::new(0.0, 10);
        assert!((classical.heat_capacity(300.0) - classical.dulong_petit_limit()).abs() < 1e-30);
    }

    #[test]
    fn test_debye_thermal_energy_derivative_is_heat_capacity() {
        // C_V = dU/dT ties the energy function D3 to the Cv function.
        let dm = DebyeModel::new(350.0, 40);
        for t in [20.0, 90.0, 350.0, 1000.0] {
            let h = 1e-3 * t;
            let derivative = (dm.internal_energy(t + h) - dm.internal_energy(t - h)) / (2.0 * h);
            let cv = dm.heat_capacity(t);
            assert!(
                (derivative - cv).abs() < 1e-5 * cv,
                "dU/dT = {derivative}, C_V = {cv} at T = {t}"
            );
        }
    }

    #[test]
    fn test_debye_thermal_energy_limits() {
        let dm = DebyeModel::new(300.0, 12);
        let nk = 12.0 * KB;
        // High temperature: equipartition 3 N kB T (minus the small x/8 correction).
        let t_hot = 1.0e5;
        assert!((dm.thermal_energy(t_hot) / (3.0 * nk * t_hot) - 1.0).abs() < 2e-3);
        // Low temperature: (3 pi^4 / 5) N kB T^4 / theta^3.
        let t_cold = 3.0_f64;
        let t4_law = 3.0 * PI.powi(4) / 5.0 * nk * t_cold.powi(4) / 300.0_f64.powi(3);
        assert!((dm.thermal_energy(t_cold) / t4_law - 1.0).abs() < 1e-9);
        // Zero temperature leaves only the zero-point energy.
        assert_eq!(dm.thermal_energy(0.0), 0.0);
        assert_eq!(dm.internal_energy(0.0), dm.zero_point_energy());
        assert!(dm.thermal_energy(f64::NAN).is_nan());
    }

    // ---- EinsteinModel ----

    #[test]
    fn test_einstein_high_temp_dulong_petit() {
        let em = EinsteinModel::new(300.0, 100);
        let cv = em.heat_capacity(10000.0);
        let dp = em.high_temperature_limit();
        assert!((cv / dp - 1.0).abs() < 0.02);
    }

    #[test]
    fn test_einstein_zpe_positive() {
        let em = EinsteinModel::new(500.0, 10);
        assert!(em.zero_point_energy() > 0.0);
    }

    #[test]
    fn test_einstein_mean_energy_above_zpe() {
        let em = EinsteinModel::new(500.0, 10);
        let e_mean = em.mean_energy(300.0);
        let zpe = em.zero_point_energy();
        assert!(e_mean >= zpe - 1e-30);
    }

    #[test]
    fn test_einstein_cv_increasing_with_temperature() {
        let em = EinsteinModel::new(300.0, 50);
        let cv1 = em.heat_capacity(100.0);
        let cv2 = em.heat_capacity(500.0);
        assert!(cv1 < cv2);
    }

    // ---- NuclearQuantumEffect ----

    #[test]
    fn test_nqe_proton_tunnelling_positive() {
        let nqe = NuclearQuantumEffect::new(3.5e14, 5e-21, 1e-10, 300.0);
        assert!(nqe.proton_tunnelling_coefficient() > 0.0);
    }

    #[test]
    fn test_nqe_proton_tunnels_more_than_deuteron() {
        let nqe = NuclearQuantumEffect::new(3.5e14, 5e-21, 5e-11, 300.0);
        assert!(nqe.proton_tunnelling_coefficient() > nqe.deuteron_tunnelling_coefficient());
    }

    #[test]
    fn test_nqe_bell_correction_ge_one() {
        let nqe = NuclearQuantumEffect::new(3.5e14, 5e-21, 1e-10, 300.0);
        let qt = nqe.bell_tunnelling_correction();
        assert!(qt >= 1.0);
    }

    #[test]
    fn test_nqe_fractionation_factor_positive() {
        let nqe = NuclearQuantumEffect::new(3.5e14, 5e-21, 5e-11, 300.0);
        assert!(nqe.hd_fractionation_factor() > 0.0);
    }

    #[test]
    fn test_nqe_delocalization_positive() {
        let nqe = NuclearQuantumEffect::new(3.5e14, 5e-21, 1e-10, 300.0);
        assert!(nqe.quantum_delocalization() > 0.0);
    }
}
