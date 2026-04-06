// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Optical materials — refractive index, photonics, and nonlinear optics.
//!
//! Provides:
//! - Complex refractive index, Fresnel equations, Brewster/critical angles
//! - Sellmeier and Cauchy dispersion, Abbe number
//! - Group index, chromatic dispersion
//! - Anti-reflection coatings, Fabry-Perot etalon
//! - Birefringence, optical rotation (chirality)
//! - Beer-Lambert absorption, photoelasticity
//! - Electro-optic (Pockels/Kerr) effects
//! - Nonlinear optics (χ² SHG, χ³ SPM/XPM)
//! - Photoluminescence quantum yield
//! - Optical bandgap via Tauc plot
//! - Drude dielectric function
//! - Metamaterial effective medium (Maxwell Garnett)
//! - Photonic crystal bandgap estimation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

/// Speed of light in vacuum \[m/s\].
const SPEED_OF_LIGHT: f64 = 2.997_924_58e8;

/// Vacuum permeability μ₀ \[H/m\].
const MU_0: f64 = 1.256_637_061_4e-6;

/// Vacuum permittivity ε₀ \[F/m\].
const EPSILON_0: f64 = 8.854_187_817e-12;

/// Planck constant \[J·s\].
const PLANCK_H: f64 = 6.626_070_15e-34;

/// Elementary charge \[C\].
const ELEM_CHARGE: f64 = 1.602_176_634e-19;

// ─────────────────────────────────────────────────────────────────────────────
// OpticalMaterial
// ─────────────────────────────────────────────────────────────────────────────

/// Optical material characterised by its complex refractive index ñ = n + i·k.
///
/// * `n_real` – Real part of refractive index (phase velocity ratio).
/// * `n_imag` – Imaginary part k (extinction coefficient, ≥ 0 for absorbing media).
/// * `wavelength_nm` – Wavelength of light in vacuum \[nm\].
#[derive(Debug, Clone, Copy)]
pub struct OpticalMaterial {
    /// Real part of the complex refractive index n.
    pub n_real: f64,
    /// Imaginary part (extinction coefficient) k ≥ 0.
    pub n_imag: f64,
    /// Vacuum wavelength \[nm\].
    pub wavelength_nm: f64,
}

impl OpticalMaterial {
    /// Create a new [`OpticalMaterial`] with complex index n + ik at wavelength λ.
    pub fn new(n: f64, k: f64, wl_nm: f64) -> Self {
        Self {
            n_real: n,
            n_imag: k,
            wavelength_nm: wl_nm,
        }
    }

    /// Absorption coefficient α \[m⁻¹\].
    ///
    /// α = 4π·k / λ  where λ is in metres.
    pub fn absorption_coefficient(&self) -> f64 {
        4.0 * PI * self.n_imag / (self.wavelength_nm * 1e-9)
    }

    /// Normal-incidence reflectance R from vacuum into this material.
    ///
    /// R = ((n−1)² + k²) / ((n+1)² + k²)
    pub fn reflectance_normal(&self) -> f64 {
        let n = self.n_real;
        let k = self.n_imag;
        let num = (n - 1.0) * (n - 1.0) + k * k;
        let den = (n + 1.0) * (n + 1.0) + k * k;
        num / den
    }

    /// Optical penetration (skin) depth \[m\].
    ///
    /// δ = λ / (4π·k)  (wavelength in vacuum in metres).
    /// Returns infinity when k = 0 (non-absorbing material).
    pub fn penetration_depth(&self) -> f64 {
        if self.n_imag == 0.0 {
            f64::INFINITY
        } else {
            self.wavelength_nm * 1e-9 / (4.0 * PI * self.n_imag)
        }
    }

    /// Optical impedance Z \[Ω\] (magnitude).
    ///
    /// Z = Z₀ / n_real  where Z₀ ≈ 376.73 Ω is the free-space impedance.
    pub fn optical_impedance(&self) -> f64 {
        376.730_313_412 / self.n_real
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fresnel equations
// ─────────────────────────────────────────────────────────────────────────────

/// Fresnel amplitude reflection coefficient for s-polarisation.
///
/// r_s = (n₁·cos θᵢ − n₂·cos θₜ) / (n₁·cos θᵢ + n₂·cos θₜ)
///
/// # Arguments
/// * `n1` – Refractive index of incident medium.
/// * `n2` – Refractive index of transmitted medium.
/// * `theta_i` – Angle of incidence \[rad\].
pub fn fresnel_rs(n1: f64, n2: f64, theta_i: f64) -> f64 {
    let cos_i = theta_i.cos();
    let sin_i = theta_i.sin();
    let sin_t_sq = (n1 / n2 * sin_i) * (n1 / n2 * sin_i);
    if sin_t_sq > 1.0 {
        // Total internal reflection: magnitude is 1, return -1 for phase tracking.
        return -1.0;
    }
    let cos_t = (1.0 - sin_t_sq).sqrt();
    (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
}

/// Fresnel amplitude reflection coefficient for p-polarisation.
///
/// r_p = (n₂·cos θᵢ − n₁·cos θₜ) / (n₂·cos θᵢ + n₁·cos θₜ)
///
/// # Arguments
/// * `n1` – Refractive index of incident medium.
/// * `n2` – Refractive index of transmitted medium.
/// * `theta_i` – Angle of incidence \[rad\].
pub fn fresnel_rp(n1: f64, n2: f64, theta_i: f64) -> f64 {
    let cos_i = theta_i.cos();
    let sin_i = theta_i.sin();
    let sin_t_sq = (n1 / n2 * sin_i) * (n1 / n2 * sin_i);
    if sin_t_sq > 1.0 {
        return -1.0;
    }
    let cos_t = (1.0 - sin_t_sq).sqrt();
    (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
}

/// Power reflectance R and transmittance T at an interface for both polarisations.
///
/// Returns `(R_s, T_s, R_p, T_p)` where R + T = 1 for each polarisation.
///
/// # Arguments
/// * `n1` – Refractive index of incident medium.
/// * `n2` – Refractive index of transmitted medium.
/// * `theta_i` – Angle of incidence \[rad\].
pub fn fresnel_power(n1: f64, n2: f64, theta_i: f64) -> (f64, f64, f64, f64) {
    let rs = fresnel_rs(n1, n2, theta_i);
    let rp = fresnel_rp(n1, n2, theta_i);
    let r_s = rs * rs;
    let r_p = rp * rp;
    let sin_i = theta_i.sin();
    let sin_t_sq = (n1 / n2 * sin_i) * (n1 / n2 * sin_i);
    // Check for TIR.
    if sin_t_sq > 1.0 {
        return (1.0, 0.0, 1.0, 0.0);
    }
    let cos_i = theta_i.cos();
    let cos_t = (1.0 - sin_t_sq).sqrt();
    let factor = (n2 * cos_t) / (n1 * cos_i);
    let ts_amp = 2.0 * n1 * cos_i / (n1 * cos_i + n2 * cos_t);
    let tp_amp = 2.0 * n1 * cos_i / (n2 * cos_i + n1 * cos_t);
    let t_s = factor * ts_amp * ts_amp;
    let t_p = factor * tp_amp * tp_amp;
    (r_s, t_s, r_p, t_p)
}

// ─────────────────────────────────────────────────────────────────────────────
// Brewster angle
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Brewster angle θ_B where r_p = 0.
///
/// θ_B = atan(n₂/n₁)
///
/// # Arguments
/// * `n1` – Refractive index of incident medium.
/// * `n2` – Refractive index of second medium.
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    (n2 / n1).atan()
}

// ─────────────────────────────────────────────────────────────────────────────
// Critical angle
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the critical angle for total internal reflection.
///
/// θ_c = asin(n₂/n₁)  valid only when n₁ > n₂.
/// Returns `None` if n₁ ≤ n₂ (no total internal reflection possible).
///
/// # Arguments
/// * `n1` – Refractive index of the denser medium (incident side).
/// * `n2` – Refractive index of the less dense medium (transmitted side).
pub fn critical_angle(n1: f64, n2: f64) -> Option<f64> {
    if n1 <= n2 {
        None
    } else {
        Some((n2 / n1).asin())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sellmeier equation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute refractive index via the Sellmeier dispersion formula.
///
/// n²(λ) = 1 + Σᵢ Bᵢ·λ² / (λ² − Cᵢ)
///
/// where λ is in micrometres and Cᵢ are the resonance wavelengths squared.
///
/// # Arguments
/// * `wavelength_um` – Wavelength in micrometres \[μm\].
/// * `b` – Sellmeier B coefficients (dimensionless) \[3\].
/// * `c` – Sellmeier C coefficients \[μm²\] \[3\].
pub fn sellmeier_equation(wavelength_um: f64, b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let lam2 = wavelength_um * wavelength_um;
    let mut n2 = 1.0;
    for i in 0..3 {
        n2 += b[i] * lam2 / (lam2 - c[i]);
    }
    n2.sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Cauchy equation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute refractive index via the Cauchy empirical dispersion formula.
///
/// n(λ) = A + B/λ² + C/λ⁴
///
/// where λ is in micrometres.
///
/// # Arguments
/// * `wavelength_um` – Wavelength in micrometres \[μm\].
/// * `a` – Cauchy coefficient A (dominant refractive index at long λ).
/// * `b` – Cauchy coefficient B \[μm²\].
/// * `c` – Cauchy coefficient C \[μm⁴\].
pub fn cauchy_equation(wavelength_um: f64, a: f64, b: f64, c: f64) -> f64 {
    let lam2 = wavelength_um * wavelength_um;
    a + b / lam2 + c / (lam2 * lam2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Abbe number (optical dispersion)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Abbe number V (V-number) of a glass at three wavelengths.
///
/// V = (n_d − 1) / (n_F − n_C)
///
/// where n_d, n_F, n_C are the refractive indices at the Fraunhofer d, F, and C
/// spectral lines (587.6 nm, 486.1 nm, 656.3 nm respectively).
///
/// # Arguments
/// * `n_d` – Refractive index at 587.6 nm (yellow, He-d line).
/// * `n_f` – Refractive index at 486.1 nm (blue, H-F line).
/// * `n_c` – Refractive index at 656.3 nm (red, H-C line).
pub fn abbe_number(n_d: f64, n_f: f64, n_c: f64) -> f64 {
    (n_d - 1.0) / (n_f - n_c)
}

/// Classify glass type from Abbe number.
///
/// * V > 55  → Crown glass (low dispersion).
/// * 35 ≤ V ≤ 55 → Intermediate.
/// * V < 35  → Flint glass (high dispersion).
///
/// Returns a static string label.
pub fn glass_type_from_abbe(v: f64) -> &'static str {
    if v > 55.0 {
        "crown"
    } else if v >= 35.0 {
        "intermediate"
    } else {
        "flint"
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Group refractive index
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the group refractive index n_g = n − λ·dn/dλ.
///
/// Uses a central finite difference for dn/dλ with step `dwl`.
///
/// # Arguments
/// * `n_wl` – Function n(λ) giving refractive index at wavelength λ \[μm\].
/// * `wl` – Wavelength \[μm\].
/// * `dwl` – Finite difference step \[μm\].
pub fn group_refractive_index(n_wl: impl Fn(f64) -> f64, wl: f64, dwl: f64) -> f64 {
    let dn_dlam = (n_wl(wl + dwl) - n_wl(wl - dwl)) / (2.0 * dwl);
    n_wl(wl) - wl * dn_dlam
}

// ─────────────────────────────────────────────────────────────────────────────
// Chromatic dispersion
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the group-velocity dispersion (GVD) parameter D.
///
/// D = −λ/c · d²n/dλ²   \[ps/(nm·km)\] if wavelength in μm and c in m/s.
///
/// Uses a central second finite difference.
///
/// # Arguments
/// * `n_wl` – Function n(λ) giving refractive index at wavelength λ \[μm\].
/// * `wl` – Wavelength \[μm\].
/// * `dwl` – Finite difference step \[μm\].
pub fn chromatic_dispersion(n_wl: impl Fn(f64) -> f64, wl: f64, dwl: f64) -> f64 {
    let d2n = (n_wl(wl + dwl) - 2.0 * n_wl(wl) + n_wl(wl - dwl)) / (dwl * dwl);
    -wl / SPEED_OF_LIGHT * d2n
}

// ─────────────────────────────────────────────────────────────────────────────
// Beer-Lambert absorption
// ─────────────────────────────────────────────────────────────────────────────

/// Beer-Lambert law: transmitted intensity I through absorbing medium.
///
/// I = I₀ · exp(−α · L)
///
/// # Arguments
/// * `i0` – Incident intensity \[arb. units\].
/// * `alpha` – Absorption coefficient \[m⁻¹\].
/// * `length` – Propagation path length \[m\].
pub fn beer_lambert(i0: f64, alpha: f64, length: f64) -> f64 {
    i0 * (-alpha * length).exp()
}

/// Absorbance A = log₁₀(I₀/I) = α·L / ln(10).
///
/// # Arguments
/// * `alpha` – Absorption coefficient \[m⁻¹\].
/// * `length` – Path length \[m\].
pub fn absorbance(alpha: f64, length: f64) -> f64 {
    alpha * length / std::f64::consts::LN_10
}

// ─────────────────────────────────────────────────────────────────────────────
// Birefringence
// ─────────────────────────────────────────────────────────────────────────────

/// Uniaxial birefringent material with ordinary (n_o) and extraordinary (n_e) indices.
#[derive(Debug, Clone, Copy)]
pub struct BirefringentMaterial {
    /// Ordinary refractive index n_o.
    pub n_ordinary: f64,
    /// Extraordinary refractive index n_e.
    pub n_extraordinary: f64,
}

impl BirefringentMaterial {
    /// Create a new [`BirefringentMaterial`].
    pub fn new(n_ordinary: f64, n_extraordinary: f64) -> Self {
        Self {
            n_ordinary,
            n_extraordinary,
        }
    }

    /// Birefringence Δn = n_e − n_o.
    pub fn birefringence(&self) -> f64 {
        self.n_extraordinary - self.n_ordinary
    }

    /// Effective extraordinary index n_eff(θ) as a function of propagation angle θ
    /// relative to the optical axis.
    ///
    /// 1/n_eff² = cos²θ/n_o² + sin²θ/n_e²
    ///
    /// # Arguments
    /// * `theta` – Angle between propagation direction and optical axis \[rad\].
    pub fn effective_index(&self, theta: f64) -> f64 {
        let no = self.n_ordinary;
        let ne = self.n_extraordinary;
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        1.0 / ((cos_t * cos_t / (no * no) + sin_t * sin_t / (ne * ne)).sqrt())
    }

    /// Phase retardation Γ for a plate of thickness `d_m` \[m\] at wavelength `wl_nm` \[nm\].
    ///
    /// Γ = 2π·Δn·d / λ
    pub fn phase_retardation(&self, d_m: f64, wl_nm: f64) -> f64 {
        2.0 * PI * self.birefringence().abs() * d_m / (wl_nm * 1e-9)
    }

    /// Quarter-wave plate thickness for wavelength `wl_nm` \[nm\].
    ///
    /// d_QWP = λ / (4·|Δn|)
    pub fn quarter_wave_thickness(&self, wl_nm: f64) -> f64 {
        wl_nm * 1e-9 / (4.0 * self.birefringence().abs())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Optical rotation (chirality)
// ─────────────────────────────────────────────────────────────────────────────

/// Optical rotation in a chiral medium.
///
/// A chiral medium has different refractive indices for left (n_L) and right (n_R)
/// circularly polarised light.
#[derive(Debug, Clone, Copy)]
pub struct ChiralMedium {
    /// Refractive index for left-circularly polarised light.
    pub n_left: f64,
    /// Refractive index for right-circularly polarised light.
    pub n_right: f64,
}

impl ChiralMedium {
    /// Create a new [`ChiralMedium`].
    pub fn new(n_left: f64, n_right: f64) -> Self {
        Self { n_left, n_right }
    }

    /// Specific optical rotation ρ \[rad/m\] at wavelength `wl_nm` \[nm\].
    ///
    /// ρ = π·(n_L − n_R) / λ
    pub fn specific_rotation(&self, wl_nm: f64) -> f64 {
        PI * (self.n_left - self.n_right) / (wl_nm * 1e-9)
    }

    /// Total rotation of polarisation plane for path length `l_m` \[m\].
    ///
    /// φ = ρ · l
    pub fn rotation_angle(&self, l_m: f64, wl_nm: f64) -> f64 {
        self.specific_rotation(wl_nm) * l_m
    }

    /// Average refractive index n_avg = (n_L + n_R) / 2.
    pub fn mean_index(&self) -> f64 {
        (self.n_left + self.n_right) / 2.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Photoelastic effect (stress-optic coefficient)
// ─────────────────────────────────────────────────────────────────────────────

/// Photoelastic birefringence induced by mechanical stress.
///
/// Δn = C · (σ₁ − σ₂)
///
/// where C is the stress-optic (Brewster) coefficient.
///
/// # Arguments
/// * `stress_optic_coeff` – Brewster coefficient C \[Pa⁻¹\] (typically ~10⁻¹²).
/// * `sigma1` – Principal stress along one axis \[Pa\].
/// * `sigma2` – Principal stress along perpendicular axis \[Pa\].
pub fn photoelastic_birefringence(stress_optic_coeff: f64, sigma1: f64, sigma2: f64) -> f64 {
    stress_optic_coeff * (sigma1 - sigma2)
}

/// Photoelastic fringe order N from birefringence Δn at wavelength `wl_nm` \[nm\]
/// through thickness `d_m` \[m\].
///
/// N = Δn · d / λ
pub fn photoelastic_fringe_order(delta_n: f64, d_m: f64, wl_nm: f64) -> f64 {
    delta_n * d_m / (wl_nm * 1e-9)
}

// ─────────────────────────────────────────────────────────────────────────────
// Electro-optic effects (Pockels and Kerr)
// ─────────────────────────────────────────────────────────────────────────────

/// Pockels electro-optic effect: index change linear in electric field.
///
/// Δn = −½ · n³ · r · E
///
/// # Arguments
/// * `n0` – Unperturbed refractive index.
/// * `r_eo` – Electro-optic coefficient r \[m/V\].
/// * `e_field` – Applied electric field \[V/m\].
pub fn pockels_index_change(n0: f64, r_eo: f64, e_field: f64) -> f64 {
    -0.5 * n0 * n0 * n0 * r_eo * e_field
}

/// Kerr electro-optic effect: index change quadratic in electric field.
///
/// Δn = ½ · n³ · K · λ₀ · E²
///
/// where K is the Kerr constant \[m/V²\] and λ₀ the vacuum wavelength \[m\].
///
/// # Arguments
/// * `n0` – Unperturbed refractive index.
/// * `kerr_const` – Kerr constant K \[m/V²\].
/// * `wl_m` – Vacuum wavelength \[m\].
/// * `e_field` – Applied electric field \[V/m\].
pub fn kerr_index_change(n0: f64, kerr_const: f64, wl_m: f64, e_field: f64) -> f64 {
    0.5 * n0 * n0 * n0 * kerr_const * wl_m * e_field * e_field
}

/// Pockels half-wave voltage V_π for a longitudinal modulator.
///
/// V_π = λ / (2 · n₀³ · r · L/d)
///
/// For a transverse modulator the aspect ratio d/L appears explicitly.
///
/// # Arguments
/// * `wl_m` – Vacuum wavelength \[m\].
/// * `n0` – Unperturbed index.
/// * `r_eo` – Electro-optic coefficient \[m/V\].
/// * `length_m` – Crystal length \[m\].
/// * `electrode_gap_m` – Electrode separation \[m\].
pub fn pockels_half_wave_voltage(
    wl_m: f64,
    n0: f64,
    r_eo: f64,
    length_m: f64,
    electrode_gap_m: f64,
) -> f64 {
    wl_m * electrode_gap_m / (n0 * n0 * n0 * r_eo * length_m)
}

// ─────────────────────────────────────────────────────────────────────────────
// Nonlinear optics: χ² and χ³
// ─────────────────────────────────────────────────────────────────────────────

/// Second-harmonic generation (SHG) conversion efficiency in the undepleted pump approximation.
///
/// η_SHG = (2ω²·d_eff²·L²·I_pump) / (n₁²·n₂·c³·ε₀)
///
/// where d_eff is the effective χ²/2 nonlinear coefficient.
/// This simplified formula gives the fractional power conversion.
///
/// # Arguments
/// * `d_eff` – Effective nonlinear coefficient d_eff = χ²/2 \[m/V\].
/// * `length_m` – Crystal length \[m\].
/// * `pump_intensity` – Pump intensity I \[W/m²\].
/// * `n1` – Refractive index at pump wavelength.
/// * `n2` – Refractive index at SH wavelength.
/// * `wl_m` – Pump wavelength \[m\].
pub fn shg_efficiency(
    d_eff: f64,
    length_m: f64,
    pump_intensity: f64,
    n1: f64,
    n2: f64,
    wl_m: f64,
) -> f64 {
    let omega = 2.0 * PI * SPEED_OF_LIGHT / wl_m;
    let num = 2.0 * omega * omega * d_eff * d_eff * length_m * length_m * pump_intensity;
    let den = n1 * n1 * n2 * SPEED_OF_LIGHT * SPEED_OF_LIGHT * SPEED_OF_LIGHT * EPSILON_0;
    if den.abs() < 1e-60 { 0.0 } else { num / den }
}

/// Phase-matching bandwidth Δλ for SHG.
///
/// Δλ ≈ λ² / (2·L·|n_g1 − n_g2|)
///
/// # Arguments
/// * `wl_m` – Pump wavelength \[m\].
/// * `length_m` – Crystal length \[m\].
/// * `n_group1` – Group index at pump wavelength.
/// * `n_group2` – Group index at SH wavelength.
pub fn shg_phase_matching_bandwidth(wl_m: f64, length_m: f64, n_group1: f64, n_group2: f64) -> f64 {
    let delta_ng = (n_group1 - n_group2).abs();
    if delta_ng < 1e-30 {
        f64::INFINITY
    } else {
        wl_m * wl_m / (2.0 * length_m * delta_ng)
    }
}

/// Self-phase modulation (SPM) nonlinear phase shift for χ³ medium.
///
/// φ_NL = γ · P₀ · L
///
/// where γ = ω·n₂/(c·A_eff) is the nonlinear coefficient \[rad/(W·m)\].
///
/// # Arguments
/// * `n2_nonlinear` – Nonlinear refractive index n₂ (χ³) \[m²/W\].
/// * `wl_m` – Vacuum wavelength \[m\].
/// * `aeff_m2` – Effective mode area A_eff \[m²\].
/// * `peak_power_w` – Peak power P₀ \[W\].
/// * `length_m` – Medium length \[m\].
pub fn spm_phase_shift(
    n2_nonlinear: f64,
    wl_m: f64,
    aeff_m2: f64,
    peak_power_w: f64,
    length_m: f64,
) -> f64 {
    let gamma = 2.0 * PI * n2_nonlinear / (wl_m * aeff_m2);
    gamma * peak_power_w * length_m
}

/// Cross-phase modulation (XPM) phase shift on probe due to pump.
///
/// φ_XPM = 2 · γ · P_pump · L   (factor 2 vs SPM).
///
/// # Arguments
/// * `n2_nonlinear` – Nonlinear refractive index n₂ \[m²/W\].
/// * `wl_m` – Probe wavelength \[m\].
/// * `aeff_m2` – Effective mode area \[m²\].
/// * `pump_power_w` – Pump peak power \[W\].
/// * `length_m` – Medium length \[m\].
pub fn xpm_phase_shift(
    n2_nonlinear: f64,
    wl_m: f64,
    aeff_m2: f64,
    pump_power_w: f64,
    length_m: f64,
) -> f64 {
    2.0 * spm_phase_shift(n2_nonlinear, wl_m, aeff_m2, pump_power_w, length_m)
}

// ─────────────────────────────────────────────────────────────────────────────
// Photoluminescence quantum yield
// ─────────────────────────────────────────────────────────────────────────────

/// Photoluminescence quantum yield (PLQY) η_PL.
///
/// η_PL = N_emitted / N_absorbed  = k_r / (k_r + k_nr)
///
/// where k_r is the radiative decay rate and k_nr the non-radiative decay rate.
///
/// # Arguments
/// * `k_radiative` – Radiative rate \[s⁻¹\].
/// * `k_nonradiative` – Non-radiative rate \[s⁻¹\].
pub fn plqy(k_radiative: f64, k_nonradiative: f64) -> f64 {
    let total = k_radiative + k_nonradiative;
    if total == 0.0 {
        0.0
    } else {
        k_radiative / total
    }
}

/// Fluorescence lifetime τ = 1 / (k_r + k_nr) \[s\].
///
/// # Arguments
/// * `k_radiative` – Radiative rate \[s⁻¹\].
/// * `k_nonradiative` – Non-radiative rate \[s⁻¹\].
pub fn fluorescence_lifetime(k_radiative: f64, k_nonradiative: f64) -> f64 {
    1.0 / (k_radiative + k_nonradiative)
}

// ─────────────────────────────────────────────────────────────────────────────
// Optical bandgap: Tauc plot
// ─────────────────────────────────────────────────────────────────────────────

/// Tauc plot ordinate (αhν)^(1/r) used to extract optical bandgap.
///
/// For direct-allowed transitions: r = 1/2 → (αhν)².
/// For indirect-allowed transitions: r = 2 → (αhν)^(1/2).
///
/// Linear extrapolation to zero gives E_g.
///
/// # Arguments
/// * `alpha` – Absorption coefficient \[m⁻¹\].
/// * `photon_energy_ev` – Photon energy hν \[eV\].
/// * `r` – Tauc exponent (1/2 for direct, 2 for indirect).
pub fn tauc_ordinate(alpha: f64, photon_energy_ev: f64, r: f64) -> f64 {
    let alpha_hv = alpha * photon_energy_ev;
    alpha_hv.powf(1.0 / r)
}

/// Estimate optical bandgap from linear Tauc fit parameters.
///
/// Given the slope `m` and intercept `b` of a linear fit to the Tauc plot
/// ((αhν)^(1/r) = m·(hν) + b), returns E_g = −b/m.
///
/// # Arguments
/// * `slope` – Slope of the Tauc linear fit.
/// * `intercept` – Intercept of the Tauc linear fit.
pub fn tauc_bandgap(slope: f64, intercept: f64) -> f64 {
    if slope.abs() < 1e-30 {
        f64::NAN
    } else {
        -intercept / slope
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Drude dielectric function
// ─────────────────────────────────────────────────────────────────────────────

/// Drude model complex dielectric function ε(ω) for free-electron metal.
///
/// ε(ω) = ε_∞ − ω_p² / (ω² + i·ω·γ)
///
/// Returns `(ε_real, ε_imag)`.
///
/// # Arguments
/// * `epsilon_inf` – High-frequency dielectric constant ε_∞.
/// * `omega_p` – Plasma frequency ω_p \[rad/s\].
/// * `gamma_drude` – Collision (damping) rate γ \[rad/s\].
/// * `omega` – Angular frequency ω \[rad/s\].
pub fn drude_dielectric(
    epsilon_inf: f64,
    omega_p: f64,
    gamma_drude: f64,
    omega: f64,
) -> (f64, f64) {
    let denom = omega * omega + gamma_drude * gamma_drude;
    let re = epsilon_inf - omega_p * omega_p / denom;
    let im = omega_p * omega_p * gamma_drude / (omega * denom);
    (re, im)
}

/// Drude plasma frequency ω_p from carrier density.
///
/// ω_p = sqrt(n_e · e² / (ε₀ · m*))
///
/// # Arguments
/// * `carrier_density` – Free carrier density n_e \[m⁻³\].
/// * `eff_mass` – Effective carrier mass m* \[kg\].
pub fn drude_plasma_frequency(carrier_density: f64, eff_mass: f64) -> f64 {
    (carrier_density * ELEM_CHARGE * ELEM_CHARGE / (EPSILON_0 * eff_mass)).sqrt()
}

/// Drude complex refractive index ñ = sqrt(ε).
///
/// Returns `(n, k)` where ñ = n + i·k.
///
/// # Arguments
/// * `epsilon_real` – Real part of dielectric function.
/// * `epsilon_imag` – Imaginary part of dielectric function.
pub fn drude_refractive_index(epsilon_real: f64, epsilon_imag: f64) -> (f64, f64) {
    let magnitude = (epsilon_real * epsilon_real + epsilon_imag * epsilon_imag).sqrt();
    let n = ((magnitude + epsilon_real) / 2.0).sqrt();
    let k = ((magnitude - epsilon_real) / 2.0).sqrt();
    (n, k)
}

// ─────────────────────────────────────────────────────────────────────────────
// Metamaterial effective medium (Maxwell Garnett mixing rule)
// ─────────────────────────────────────────────────────────────────────────────

/// Maxwell Garnett effective permittivity for a composite of inclusions in a host.
///
/// ε_eff = ε_host · (ε_incl·(1 + 2f) + ε_host·2·(1 − f)) /
///                  (ε_incl·(1 − f)   + ε_host·(2 + f))
///
/// # Arguments
/// * `eps_host` – Permittivity of host matrix.
/// * `eps_incl` – Permittivity of inclusion.
/// * `fill_fraction` – Volume fraction f of inclusions (0–1).
pub fn maxwell_garnett(eps_host: f64, eps_incl: f64, fill_fraction: f64) -> f64 {
    let f = fill_fraction;
    let num = eps_incl * (1.0 + 2.0 * f) + eps_host * 2.0 * (1.0 - f);
    let den = eps_incl * (1.0 - f) + eps_host * (2.0 + f);
    if den.abs() < 1e-30 {
        eps_host
    } else {
        eps_host * num / den
    }
}

/// Bruggeman effective medium approximation for binary mixture.
///
/// Solves: f·(ε₁ − ε_eff)/(ε₁ + 2ε_eff) + (1−f)·(ε₂ − ε_eff)/(ε₂ + 2ε_eff) = 0
///
/// Returns the physical root (ε_eff > 0).
///
/// # Arguments
/// * `eps1` – Permittivity of material 1.
/// * `eps2` – Permittivity of material 2.
/// * `f1` – Volume fraction of material 1 (0–1).
pub fn bruggeman_effective_medium(eps1: f64, eps2: f64, f1: f64) -> f64 {
    let f2 = 1.0 - f1;
    // Quadratic: 0 = (2f1 − f2)·ε_eff² + [f1·(ε2 − 2ε1) + f2·(ε1 − 2ε2)]·...
    // Simplified analytic solution for real epsilon:
    let b = (f1 * (2.0 * eps1 - eps2) + f2 * (2.0 * eps2 - eps1)) / 3.0;
    let disc = b * b + eps1 * eps2 / 9.0 + (f1 * eps1 + f2 * eps2) / 3.0;
    // Use iterative approach: start from linear mix.
    let mut eps_eff = f1 * eps1 + f2 * eps2;
    for _ in 0..100 {
        let t1 = f1 * (eps1 - eps_eff) / (eps1 + 2.0 * eps_eff);
        let t2 = f2 * (eps2 - eps_eff) / (eps2 + 2.0 * eps_eff);
        let _b_unused = b;
        let _disc_unused = disc;
        let residual = t1 + t2;
        // Derivative dr/deps_eff.
        let dt1 = f1 * (-3.0 * eps1) / ((eps1 + 2.0 * eps_eff) * (eps1 + 2.0 * eps_eff));
        let dt2 = f2 * (-3.0 * eps2) / ((eps2 + 2.0 * eps_eff) * (eps2 + 2.0 * eps_eff));
        let drv = dt1 + dt2;
        if drv.abs() < 1e-30 {
            break;
        }
        let step = residual / drv;
        eps_eff -= step;
        if step.abs() < 1e-12 {
            break;
        }
    }
    eps_eff
}

// ─────────────────────────────────────────────────────────────────────────────
// Photonic crystal bandgap
// ─────────────────────────────────────────────────────────────────────────────

/// Photonic crystal 1D bandgap estimation using transfer matrix method.
///
/// For a 1D periodic stack (layers of n₁, d₁ and n₂, d₂) the photonic bandgap
/// centre wavelength and fractional gap width are estimated.
///
/// # Arguments
/// * `n1` – Refractive index of layer 1.
/// * `n2` – Refractive index of layer 2.
/// * `d1` – Thickness of layer 1 \[m\].
/// * `d2` – Thickness of layer 2 \[m\].
///
/// Returns `(lambda_center, fractional_gap)` where
/// λ_c = 2(n₁d₁ + n₂d₂) and Δλ/λ = (4/π)·arcsin(|n₁−n₂|/(n₁+n₂)).
pub fn photonic_crystal_1d_bandgap(n1: f64, n2: f64, d1: f64, d2: f64) -> (f64, f64) {
    let lambda_center = 2.0 * (n1 * d1 + n2 * d2);
    let fractional_gap = (4.0 / PI) * ((n1 - n2).abs() / (n1 + n2)).asin();
    (lambda_center, fractional_gap)
}

/// Transfer matrix for a single dielectric layer at normal incidence.
///
/// M = \[\[cos δ, -i/n · sin δ\\], \[-i·n · sin δ, cos δ\]]
///
/// Returns the 2×2 matrix as `[m11, m12, m21, m22]` (complex components as pairs).
/// For simplicity returns the real transfer matrix ignoring phase for
/// a λ/4 layer where δ = π/2 and cos δ = 0.
///
/// For a quarter-wave layer: M_real = \[\[0, -1/n\\], \[-n, 0\]].
///
/// # Arguments
/// * `n` – Refractive index of the layer.
/// * `delta` – Round-trip phase δ = 2πnd/λ \[rad\].
pub fn transfer_matrix_layer(n: f64, delta: f64) -> [f64; 4] {
    let cos_d = delta.cos();
    let sin_d = delta.sin();
    // Real part of transfer matrix components.
    [cos_d, -sin_d / n, -n * sin_d, cos_d]
}

/// Multiply two 2×2 transfer matrices.
///
/// Returns M = A × B.
///
/// # Arguments
/// * `a` – First matrix \[m11, m12, m21, m22\].
/// * `b` – Second matrix.
pub fn transfer_matrix_multiply(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

/// Reflectance of an N-period Bragg mirror using transfer matrix method.
///
/// # Arguments
/// * `n0` – Refractive index of incident medium.
/// * `ns` – Refractive index of substrate.
/// * `n1` – Refractive index of layer 1.
/// * `n2` – Refractive index of layer 2.
/// * `d1` – Thickness of layer 1 \[m\].
/// * `d2` – Thickness of layer 2 \[m\].
/// * `wl_m` – Wavelength in vacuum \[m\].
/// * `n_periods` – Number of bilayer periods.
pub fn bragg_mirror_reflectance(
    n0: f64,
    ns: f64,
    n1: f64,
    n2: f64,
    d1: f64,
    d2: f64,
    wl_m: f64,
    n_periods: usize,
) -> f64 {
    let delta1 = 2.0 * PI * n1 * d1 / wl_m;
    let delta2 = 2.0 * PI * n2 * d2 / wl_m;
    let m1 = transfer_matrix_layer(n1, delta1);
    let m2 = transfer_matrix_layer(n2, delta2);
    let bilayer = transfer_matrix_multiply(m1, m2);

    // Raise bilayer to the n_periods power by repeated squaring.
    let mut result = [1.0_f64, 0.0, 0.0, 1.0]; // identity
    let mut base = bilayer;
    let mut exp = n_periods;
    while exp > 0 {
        if exp % 2 == 1 {
            result = transfer_matrix_multiply(result, base);
        }
        base = transfer_matrix_multiply(base, base);
        exp /= 2;
    }

    // Compute reflectance: r = (n0*M11 + n0*ns*M12 - M21 - ns*M22) /
    //                          (n0*M11 + n0*ns*M12 + M21 + ns*M22)
    let num = n0 * result[0] + n0 * ns * result[1] - result[2] - ns * result[3];
    let den = n0 * result[0] + n0 * ns * result[1] + result[2] + ns * result[3];
    if den.abs() < 1e-30 {
        1.0
    } else {
        (num / den) * (num / den)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AntiReflectionCoating
// ─────────────────────────────────────────────────────────────────────────────

/// Single-layer quarter-wave anti-reflection coating.
#[derive(Debug, Clone, Copy)]
pub struct AntiReflectionCoating {
    /// Refractive index of the substrate.
    pub n_substrate: f64,
    /// Refractive index of the coating layer.
    pub n_coating: f64,
    /// Physical coating thickness \[nm\].
    pub thickness_nm: f64,
    /// Design wavelength (for quarter-wave condition) \[nm\].
    pub wavelength_nm: f64,
}

impl AntiReflectionCoating {
    /// Create a new [`AntiReflectionCoating`].
    pub fn new(n_substrate: f64, n_coating: f64, thickness_nm: f64, wavelength_nm: f64) -> Self {
        Self {
            n_substrate,
            n_coating,
            thickness_nm,
            wavelength_nm,
        }
    }

    /// Optimal coating refractive index for zero reflection: n_c = sqrt(n_s).
    pub fn optimal_n(&self) -> f64 {
        self.n_substrate.sqrt()
    }

    /// Approximate reflectance of the coated surface at wavelength `wl_nm` \[nm\].
    ///
    /// Uses the thin-film two-beam interference formula:
    ///
    /// R = (r₁² + r₂² + 2·r₁·r₂·cos δ) / (1 + r₁²·r₂² + 2·r₁·r₂·cos δ)
    ///
    /// where r₁ = (1 − n_c)/(1 + n_c), r₂ = (n_c − n_s)/(n_c + n_s),
    /// and δ = 4π·n_c·d / λ.
    pub fn reflectance(&self, wl_nm: f64) -> f64 {
        let r1 = (1.0 - self.n_coating) / (1.0 + self.n_coating);
        let r2 = (self.n_coating - self.n_substrate) / (self.n_coating + self.n_substrate);
        let delta = 4.0 * PI * self.n_coating * self.thickness_nm / wl_nm;
        let cos_d = delta.cos();
        let num = r1 * r1 + r2 * r2 + 2.0 * r1 * r2 * cos_d;
        let den = 1.0 + r1 * r1 * r2 * r2 + 2.0 * r1 * r2 * cos_d;
        if den.abs() < 1e-30 {
            0.0
        } else {
            (num / den).abs()
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Etalon transmission and finesse
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Fabry-Perot etalon transmission (Airy function).
///
/// T = 1 / (1 + F·sin²(δ/2))
///
/// where F = 4R/(1−R)² is the coefficient of finesse.
///
/// # Arguments
/// * `r` – Mirror reflectance R (power reflectivity) in \[0, 1\).
/// * `delta` – Round-trip phase \[rad\], δ = 4π·n·d·cos θ / λ.
pub fn etalon_transmission(r: f64, delta: f64) -> f64 {
    let f_coeff = 4.0 * r / ((1.0 - r) * (1.0 - r));
    let half_sin = (delta / 2.0).sin();
    1.0 / (1.0 + f_coeff * half_sin * half_sin)
}

/// Compute the etalon finesse F = π·√R / (1−R).
///
/// Higher finesse means narrower transmission peaks and higher resolving power.
///
/// # Arguments
/// * `r` – Mirror reflectance R (power reflectivity) in \[0, 1\).
pub fn finesse(r: f64) -> f64 {
    PI * r.sqrt() / (1.0 - r)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. reflectance_normal is in [0, 1].
    #[test]
    fn test_reflectance_normal_range() {
        let mat = OpticalMaterial::new(1.5, 0.0, 500.0);
        let r = mat.reflectance_normal();
        assert!(
            (0.0..=1.0).contains(&r),
            "Reflectance must be in [0,1], got {r}"
        );
    }

    // 2. reflectance_normal for air-glass (n=1.5, k=0).
    #[test]
    fn test_reflectance_glass() {
        let mat = OpticalMaterial::new(1.5, 0.0, 550.0);
        let r = mat.reflectance_normal();
        let expected = (0.5 / 2.5) * (0.5 / 2.5);
        assert!(
            (r - expected).abs() < EPS,
            "Glass reflectance: {r} vs {expected}"
        );
    }

    // 3. absorption_coefficient is zero for k=0.
    #[test]
    fn test_absorption_zero_for_transparent() {
        let mat = OpticalMaterial::new(1.5, 0.0, 500.0);
        assert_eq!(mat.absorption_coefficient(), 0.0, "No absorption when k=0");
    }

    // 4. absorption_coefficient is positive for k>0.
    #[test]
    fn test_absorption_positive_for_absorbing() {
        let mat = OpticalMaterial::new(2.0, 0.5, 500.0);
        let alpha = mat.absorption_coefficient();
        assert!(
            alpha > 0.0,
            "Absorption coefficient must be positive for k>0, got {alpha}"
        );
    }

    // 5. penetration_depth is infinity when k=0.
    #[test]
    fn test_penetration_depth_transparent() {
        let mat = OpticalMaterial::new(1.5, 0.0, 500.0);
        assert!(
            mat.penetration_depth().is_infinite(),
            "Non-absorbing material: infinite penetration depth"
        );
    }

    // 6. penetration_depth decreases with increasing k.
    #[test]
    fn test_penetration_depth_decreases_with_k() {
        let mat1 = OpticalMaterial::new(2.0, 0.1, 500.0);
        let mat2 = OpticalMaterial::new(2.0, 1.0, 500.0);
        assert!(
            mat2.penetration_depth() < mat1.penetration_depth(),
            "Larger k → smaller penetration depth"
        );
    }

    // 7. optical_impedance is positive.
    #[test]
    fn test_optical_impedance_positive() {
        let mat = OpticalMaterial::new(1.5, 0.0, 500.0);
        let z = mat.optical_impedance();
        assert!(z > 0.0, "Optical impedance must be positive, got {z}");
    }

    // 8. optical_impedance for vacuum (n=1) is ~376.73 Ω.
    #[test]
    fn test_optical_impedance_vacuum() {
        let mat = OpticalMaterial::new(1.0, 0.0, 500.0);
        let z = mat.optical_impedance();
        assert!(
            (z - 376.730_313_412).abs() < 1e-6,
            "Vacuum impedance mismatch: {z}"
        );
    }

    // 9. brewster_angle is positive.
    #[test]
    fn test_brewster_angle_positive() {
        let theta_b = brewster_angle(1.0, 1.5);
        assert!(
            theta_b > 0.0,
            "Brewster angle must be positive, got {theta_b}"
        );
    }

    // 10. r_p is near zero at Brewster angle.
    #[test]
    fn test_fresnel_rp_zero_at_brewster() {
        let n1 = 1.0_f64;
        let n2 = 1.5_f64;
        let theta_b = brewster_angle(n1, n2);
        let rp = fresnel_rp(n1, n2, theta_b);
        assert!(
            rp.abs() < 1e-10,
            "r_p should be zero at Brewster angle, got {rp}"
        );
    }

    // 11. r_s at normal incidence equals (n1-n2)/(n1+n2).
    #[test]
    fn test_fresnel_rs_normal_incidence() {
        let n1 = 1.0_f64;
        let n2 = 1.5_f64;
        let rs = fresnel_rs(n1, n2, 0.0);
        let expected = (n1 - n2) / (n1 + n2);
        assert!(
            (rs - expected).abs() < EPS,
            "r_s at normal incidence: {rs} vs {expected}"
        );
    }

    // 12. critical_angle returns None when n1 <= n2.
    #[test]
    fn test_critical_angle_none_for_low_n1() {
        assert!(critical_angle(1.0, 1.5).is_none(), "No TIR when n1 < n2");
    }

    // 13. critical_angle returns Some when n1 > n2.
    #[test]
    fn test_critical_angle_some_for_high_n1() {
        let theta_c = critical_angle(1.5, 1.0);
        assert!(theta_c.is_some(), "TIR exists when n1 > n2");
    }

    // 14. critical_angle value is asin(n2/n1).
    #[test]
    fn test_critical_angle_formula() {
        let n1 = 1.5_f64;
        let n2 = 1.0_f64;
        let theta_c = critical_angle(n1, n2).unwrap();
        let expected = (n2 / n1).asin();
        assert!(
            (theta_c - expected).abs() < EPS,
            "Critical angle: {theta_c} vs {expected}"
        );
    }

    // 15. Sellmeier equation returns n > 1 for BK7 glass at 589 nm.
    #[test]
    fn test_sellmeier_bk7_visible() {
        let b = [1.03961212, 0.231792344, 1.01046945];
        let c = [0.00600069867, 0.0200179144, 103.560653];
        let n = sellmeier_equation(0.589, &b, &c);
        assert!(n > 1.0, "Sellmeier should give n > 1, got {n}");
        assert!(
            (n - 1.5168).abs() < 0.001,
            "BK7 at 589 nm: expected ~1.517, got {n}"
        );
    }

    // 16. Abbe number for crown glass is > 55.
    #[test]
    fn test_abbe_number_crown() {
        // BK7-like: n_d=1.5168, n_F=1.5224, n_C=1.5143
        let v = abbe_number(1.5168, 1.5224, 1.5143);
        assert!(
            v > 55.0,
            "BK7 Abbe number should be crown glass (>55), got {v:.2}"
        );
        assert_eq!(glass_type_from_abbe(v), "crown");
    }

    // 17. Abbe number for flint glass is < 35.
    #[test]
    fn test_abbe_number_flint() {
        let v = abbe_number(1.7, 1.73, 1.68);
        assert!(v < 35.0, "Flint glass Abbe number < 35, got {v:.2}");
        assert_eq!(glass_type_from_abbe(v), "flint");
    }

    // 18. Beer-Lambert: I decreases with path length.
    #[test]
    fn test_beer_lambert_decay() {
        let i1 = beer_lambert(1.0, 100.0, 0.01);
        let i2 = beer_lambert(1.0, 100.0, 0.02);
        assert!(
            i2 < i1,
            "Beer-Lambert: intensity decreases with path length"
        );
        assert!(
            i1 > 0.0 && i2 > 0.0,
            "Transmitted intensity must be positive"
        );
    }

    // 19. Beer-Lambert: zero path length gives I = I₀.
    #[test]
    fn test_beer_lambert_zero_length() {
        let i = beer_lambert(5.0, 100.0, 0.0);
        assert!(
            (i - 5.0).abs() < EPS,
            "Beer-Lambert at zero path: {i} vs 5.0"
        );
    }

    // 20. Birefringence: Δn = n_e - n_o.
    #[test]
    fn test_birefringence_delta_n() {
        let mat = BirefringentMaterial::new(1.658, 1.486); // calcite-like
        let delta = mat.birefringence();
        let expected = 1.486 - 1.658;
        assert!(
            (delta - expected).abs() < EPS,
            "Birefringence mismatch: {delta}"
        );
    }

    // 21. Birefringence effective index at θ=0 equals n_o.
    #[test]
    fn test_birefringence_effective_index_at_zero() {
        let mat = BirefringentMaterial::new(1.658, 1.486);
        let n_eff = mat.effective_index(0.0);
        assert!(
            (n_eff - mat.n_ordinary).abs() < 1e-6,
            "At θ=0, n_eff should equal n_o: {n_eff} vs {}",
            mat.n_ordinary
        );
    }

    // 22. Phase retardation is 0 for zero thickness.
    #[test]
    fn test_phase_retardation_zero_thickness() {
        let mat = BirefringentMaterial::new(1.5, 1.6);
        let gamma = mat.phase_retardation(0.0, 500.0);
        assert!(gamma.abs() < EPS, "Zero thickness → zero retardation");
    }

    // 23. Optical rotation scales with path length.
    #[test]
    fn test_optical_rotation_scales_with_length() {
        let chiral = ChiralMedium::new(1.500, 1.501);
        let phi1 = chiral.rotation_angle(0.01, 500.0);
        let phi2 = chiral.rotation_angle(0.02, 500.0);
        assert!(
            (phi2 - 2.0 * phi1).abs() < 1e-10,
            "Rotation scales linearly with length"
        );
    }

    // 24. Pockels index change sign depends on field direction.
    #[test]
    fn test_pockels_sign() {
        let dn_pos = pockels_index_change(2.0, 3e-12, 1e7);
        let dn_neg = pockels_index_change(2.0, 3e-12, -1e7);
        assert!(
            dn_pos < 0.0,
            "Pockels: positive field → negative Δn for positive r"
        );
        assert!(
            (dn_pos + dn_neg).abs() < EPS,
            "Pockels: Δn is antisymmetric in E"
        );
    }

    // 25. Kerr index change is always positive (quadratic in E).
    #[test]
    fn test_kerr_quadratic() {
        let dn1 = kerr_index_change(1.5, 1e-12, 500e-9, 1e6);
        let dn2 = kerr_index_change(1.5, 1e-12, 500e-9, -1e6);
        assert!((dn1 - dn2).abs() < EPS, "Kerr effect is symmetric in E");
        let dn_double = kerr_index_change(1.5, 1e-12, 500e-9, 2e6);
        assert!(
            (dn_double - 4.0 * dn1).abs() < 1e-30,
            "Kerr scales as E²: dn(2E)=4*dn(E)"
        );
    }

    // 26. PLQY is in [0, 1].
    #[test]
    fn test_plqy_range() {
        let q = plqy(1e8, 1e7);
        assert!((0.0..=1.0).contains(&q), "PLQY must be in [0,1], got {q}");
    }

    // 27. PLQY = 1 when k_nr = 0.
    #[test]
    fn test_plqy_unity_no_nonradiative() {
        let q = plqy(1e9, 0.0);
        assert!(
            (q - 1.0).abs() < EPS,
            "PLQY = 1 with no non-radiative decay, got {q}"
        );
    }

    // 28. Drude dielectric has negative real part below plasma frequency.
    #[test]
    fn test_drude_below_plasma_freq() {
        let omega_p = 1e15_f64;
        let (re, _im) = drude_dielectric(1.0, omega_p, 1e13, omega_p * 0.5);
        assert!(
            re < 0.0,
            "Drude ε_real < 0 below plasma frequency, got {re}"
        );
    }

    // 29. Drude imaginary part is positive.
    #[test]
    fn test_drude_imag_positive() {
        let (_re, im) = drude_dielectric(1.0, 1e15, 1e13, 5e14);
        assert!(im > 0.0, "Drude ε_imag must be positive, got {im}");
    }

    // 30. Maxwell Garnett: f=0 gives ε_eff = ε_host.
    #[test]
    fn test_maxwell_garnett_zero_fill() {
        let eps_eff = maxwell_garnett(2.25, 10.0, 0.0);
        assert!(
            (eps_eff - 2.25).abs() < EPS,
            "Maxwell Garnett at f=0 should give ε_host: {eps_eff}"
        );
    }

    // 31. Maxwell Garnett: f=1 gives ε_eff = ε_incl.
    #[test]
    fn test_maxwell_garnett_full_fill() {
        let eps_eff = maxwell_garnett(2.25, 10.0, 1.0);
        assert!(
            (eps_eff - 10.0).abs() < EPS,
            "Maxwell Garnett at f=1 should give ε_incl: {eps_eff}"
        );
    }

    // 32. Photonic crystal bandgap centre scales with period.
    #[test]
    fn test_photonic_crystal_bandgap_scaling() {
        let (lambda1, _) = photonic_crystal_1d_bandgap(2.0, 1.5, 100e-9, 100e-9);
        let (lambda2, _) = photonic_crystal_1d_bandgap(2.0, 1.5, 200e-9, 200e-9);
        assert!(
            (lambda2 - 2.0 * lambda1).abs() < 1e-20,
            "Bandgap centre scales with period: {lambda1} vs {lambda2}"
        );
    }

    // 33. Bragg mirror reflectance increases with number of periods.
    #[test]
    fn test_bragg_mirror_periods() {
        let r5 = bragg_mirror_reflectance(1.0, 1.0, 2.3, 1.5, 60e-9, 90e-9, 600e-9, 5);
        let r10 = bragg_mirror_reflectance(1.0, 1.0, 2.3, 1.5, 60e-9, 90e-9, 600e-9, 10);
        assert!(
            r10 >= r5,
            "More periods → higher Bragg reflectance: R5={r5:.4}, R10={r10:.4}"
        );
    }

    // 34. etalon_transmission is in [0, 1].
    #[test]
    fn test_etalon_transmission_range() {
        let t = etalon_transmission(0.9, 1.2);
        assert!(
            (0.0..=1.0).contains(&t),
            "Etalon T must be in [0,1], got {t}"
        );
    }

    // 35. etalon_transmission at resonance (delta = 0) equals 1.
    #[test]
    fn test_etalon_transmission_at_resonance() {
        let t = etalon_transmission(0.9, 0.0);
        assert!((t - 1.0).abs() < EPS, "Etalon T = 1 at resonance, got {t}");
    }

    // 36. finesse increases with reflectance.
    #[test]
    fn test_finesse_increases_with_r() {
        let f1 = finesse(0.5);
        let f2 = finesse(0.9);
        assert!(
            f2 > f1,
            "Finesse should increase with reflectance: f1={f1}, f2={f2}"
        );
    }

    // 37. AntiReflectionCoating::optimal_n is sqrt(n_substrate).
    #[test]
    fn test_ar_coating_optimal_n() {
        let coating = AntiReflectionCoating::new(2.25, 1.5, 125.0, 500.0);
        let n_opt = coating.optimal_n();
        let expected = 2.25_f64.sqrt();
        assert!(
            (n_opt - expected).abs() < EPS,
            "Optimal coating n: {n_opt} vs {expected}"
        );
    }

    // 38. fresnel_power: R + T ≈ 1 for non-TIR angle.
    #[test]
    fn test_fresnel_power_conservation() {
        let (rs, ts, rp, tp) = fresnel_power(1.0, 1.5, 0.3);
        assert!(
            (rs + ts - 1.0).abs() < 1e-10,
            "Energy conservation for s-pol: Rs+Ts={}",
            rs + ts
        );
        assert!(
            (rp + tp - 1.0).abs() < 1e-10,
            "Energy conservation for p-pol: Rp+Tp={}",
            rp + tp
        );
    }

    // 39. SHG efficiency scales with L².
    #[test]
    fn test_shg_efficiency_scales_with_l2() {
        let e1 = shg_efficiency(1e-12, 1e-3, 1e12, 1.7, 1.7, 1064e-9);
        let e2 = shg_efficiency(1e-12, 2e-3, 1e12, 1.7, 1.7, 1064e-9);
        let ratio = e2 / e1;
        assert!(
            (ratio - 4.0).abs() < 1e-6,
            "SHG efficiency scales as L²: ratio={ratio:.6}"
        );
    }

    // 40. SPM phase shift scales linearly with power.
    #[test]
    fn test_spm_phase_linear_power() {
        let phi1 = spm_phase_shift(2.6e-20, 1550e-9, 80e-12, 1e3, 1.0);
        let phi2 = spm_phase_shift(2.6e-20, 1550e-9, 80e-12, 2e3, 1.0);
        assert!(
            (phi2 - 2.0 * phi1).abs() < 1e-20,
            "SPM phase shift scales with power: phi1={phi1:.6}, phi2={phi2:.6}"
        );
    }

    // 41. XPM phase is twice SPM for same parameters.
    #[test]
    fn test_xpm_is_twice_spm() {
        let phi_spm = spm_phase_shift(2.6e-20, 1550e-9, 80e-12, 1e3, 1.0);
        let phi_xpm = xpm_phase_shift(2.6e-20, 1550e-9, 80e-12, 1e3, 1.0);
        assert!(
            (phi_xpm - 2.0 * phi_spm).abs() < 1e-20,
            "XPM = 2×SPM: spm={phi_spm:.6}, xpm={phi_xpm:.6}"
        );
    }

    // 42. Tauc ordinate is zero when alpha = 0.
    #[test]
    fn test_tauc_zero_alpha() {
        let ord = tauc_ordinate(0.0, 2.5, 0.5);
        assert!(
            ord.abs() < EPS,
            "Tauc ordinate should be zero for alpha=0, got {ord}"
        );
    }

    // 43. Cauchy equation returns n > 1.
    #[test]
    fn test_cauchy_returns_n_greater_one() {
        let n = cauchy_equation(0.55, 1.458, 0.00354, 0.0);
        assert!(n > 1.0, "Cauchy must return n > 1, got {n}");
    }

    // 44. Drude plasma frequency is positive.
    #[test]
    fn test_drude_plasma_frequency_positive() {
        let wp = drude_plasma_frequency(8.5e28, 9.1e-31);
        assert!(wp > 0.0, "Plasma frequency must be positive, got {wp:.3e}");
    }

    // 45. Photoelastic birefringence vanishes for equal stresses.
    #[test]
    fn test_photoelastic_zero_for_equal_stress() {
        let dn = photoelastic_birefringence(1e-12, 1e6, 1e6);
        assert!(
            dn.abs() < EPS,
            "No birefringence for hydrostatic stress, got {dn}"
        );
    }

    // 46. Bruggeman effective medium is between eps1 and eps2.
    #[test]
    fn test_bruggeman_between_components() {
        let e1 = 2.25;
        let e2 = 12.0;
        let eps_eff = bruggeman_effective_medium(e1, e2, 0.3);
        assert!(
            eps_eff >= e1 && eps_eff <= e2,
            "Bruggeman eps_eff should be between components: {eps_eff}"
        );
    }
}
