// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Semiconductor material models.
//!
//! Covers band gap (direct/indirect), carrier concentration (intrinsic/extrinsic),
//! mobility models (Caughey-Thomas), resistivity, Hall effect, Fermi level,
//! Debye length, PN junction characteristics, Schottky barrier, photovoltaic
//! (Shockley-Queisser limit), and thermoelectric effects (Seebeck/Peltier/ZT).
//!
//! All quantities are in SI units unless otherwise noted.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J/K).
const K_B: f64 = 1.380649e-23;
/// Elementary charge (C).
const Q_E: f64 = 1.602176634e-19;
/// Electron rest mass (kg).
const M_E: f64 = 9.1093837015e-31;
/// Planck constant (J·s).
const H_PLANCK: f64 = 6.62607015e-34;
/// Reduced Planck constant (J·s).
const H_BAR: f64 = 1.054571817e-34;
/// Vacuum permittivity (F/m).
const EPSILON_0: f64 = 8.854187817e-12;
/// Speed of light (m/s).
const C_LIGHT: f64 = 2.99792458e8;
/// Stefan-Boltzmann constant (W/m²·K⁴).
const SIGMA_SB: f64 = 5.670374419e-8;

// ---------------------------------------------------------------------------
// Band gap
// ---------------------------------------------------------------------------

/// Type of band gap transition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BandGapType {
    /// Direct band gap (e.g., GaAs).
    Direct,
    /// Indirect band gap (e.g., Si).
    Indirect,
}

/// Band gap model with Varshni temperature dependence.
///
/// E_g(T) = E_g(0) - α T² / (T + β)
#[derive(Clone, Debug)]
pub struct BandGapModel {
    /// Band gap at T = 0 K (eV).
    pub eg0_ev: f64,
    /// Varshni α parameter (eV/K).
    pub alpha: f64,
    /// Varshni β parameter (K).
    pub beta: f64,
    /// Band gap type.
    pub gap_type: BandGapType,
}

impl BandGapModel {
    /// Create a new band gap model.
    pub fn new(eg0_ev: f64, alpha: f64, beta: f64, gap_type: BandGapType) -> Self {
        Self {
            eg0_ev,
            alpha,
            beta,
            gap_type,
        }
    }

    /// Silicon band gap model (Varshni parameters).
    pub fn silicon() -> Self {
        Self::new(1.166, 4.73e-4, 636.0, BandGapType::Indirect)
    }

    /// Gallium arsenide band gap model.
    pub fn gaas() -> Self {
        Self::new(1.519, 5.405e-4, 204.0, BandGapType::Direct)
    }

    /// Germanium band gap model.
    pub fn germanium() -> Self {
        Self::new(0.7437, 4.774e-4, 235.0, BandGapType::Indirect)
    }

    /// Band gap at temperature `t_kelvin` (eV).
    pub fn band_gap_ev(&self, t_kelvin: f64) -> f64 {
        self.eg0_ev - self.alpha * t_kelvin * t_kelvin / (t_kelvin + self.beta)
    }

    /// Band gap at temperature `t_kelvin` (Joules).
    pub fn band_gap_j(&self, t_kelvin: f64) -> f64 {
        self.band_gap_ev(t_kelvin) * Q_E
    }

    /// Returns whether this is a direct gap semiconductor.
    pub fn is_direct(&self) -> bool {
        self.gap_type == BandGapType::Direct
    }
}

// ---------------------------------------------------------------------------
// Effective density of states
// ---------------------------------------------------------------------------

/// Compute the effective density of states in the conduction band (1/m³).
///
/// N_C = 2 (2π m_e* k_B T / h²)^(3/2)
pub fn effective_dos_conduction(m_eff_ratio: f64, t_kelvin: f64) -> f64 {
    let m_eff = m_eff_ratio * M_E;
    2.0 * (2.0 * PI * m_eff * K_B * t_kelvin / (H_PLANCK * H_PLANCK)).powf(1.5)
}

/// Compute the effective density of states in the valence band (1/m³).
///
/// N_V = 2 (2π m_h* k_B T / h²)^(3/2)
pub fn effective_dos_valence(m_hole_ratio: f64, t_kelvin: f64) -> f64 {
    let m_hole = m_hole_ratio * M_E;
    2.0 * (2.0 * PI * m_hole * K_B * t_kelvin / (H_PLANCK * H_PLANCK)).powf(1.5)
}

// ---------------------------------------------------------------------------
// Intrinsic carrier concentration
// ---------------------------------------------------------------------------

/// Intrinsic carrier concentration n_i (1/m³).
///
/// n_i = sqrt(N_C * N_V) * exp(-E_g / (2 k_B T))
pub fn intrinsic_carrier_concentration(nc: f64, nv: f64, eg_j: f64, t_kelvin: f64) -> f64 {
    (nc * nv).sqrt() * (-eg_j / (2.0 * K_B * t_kelvin)).exp()
}

/// Convenience: intrinsic carrier concentration for silicon at given T.
pub fn intrinsic_carrier_si(t_kelvin: f64) -> f64 {
    let bg = BandGapModel::silicon();
    let eg_j = bg.band_gap_j(t_kelvin);
    // Effective mass ratios for Si: m_de* ≈ 1.08, m_dh* ≈ 0.81
    let nc = effective_dos_conduction(1.08, t_kelvin);
    let nv = effective_dos_valence(0.81, t_kelvin);
    intrinsic_carrier_concentration(nc, nv, eg_j, t_kelvin)
}

// ---------------------------------------------------------------------------
// Extrinsic carrier concentration
// ---------------------------------------------------------------------------

/// Extrinsic carrier concentrations for an n-type semiconductor.
///
/// Returns (n, p) where n is electron concentration, p is hole concentration.
/// Assumes complete ionisation: n ≈ N_D, p = n_i² / n.
pub fn extrinsic_n_type(n_d: f64, ni: f64) -> (f64, f64) {
    let n = n_d;
    let p = ni * ni / n;
    (n, p)
}

/// Extrinsic carrier concentrations for a p-type semiconductor.
///
/// Returns (n, p) where p ≈ N_A, n = n_i² / p.
pub fn extrinsic_p_type(n_a: f64, ni: f64) -> (f64, f64) {
    let p = n_a;
    let n = ni * ni / p;
    (n, p)
}

/// Compensated semiconductor carrier concentration.
///
/// For N_D > N_A: n ≈ N_D - N_A, p = n_i² / n.
/// For N_A > N_D: p ≈ N_A - N_D, n = n_i² / p.
pub fn extrinsic_compensated(n_d: f64, n_a: f64, ni: f64) -> (f64, f64) {
    if n_d > n_a {
        let n = n_d - n_a;
        let p = ni * ni / n;
        (n, p)
    } else {
        let p = n_a - n_d;
        let n = ni * ni / p;
        (n, p)
    }
}

// ---------------------------------------------------------------------------
// Fermi level
// ---------------------------------------------------------------------------

/// Intrinsic Fermi level position relative to conduction band edge (eV).
///
/// E_i = (E_c + E_v)/2 + (k_B T / 2) ln(N_V / N_C)
/// Here we return E_F - E_v for intrinsic: E_g/2 + (k_B T / 2) ln(N_V / N_C).
pub fn intrinsic_fermi_level_ev(eg_ev: f64, nc: f64, nv: f64, t_kelvin: f64) -> f64 {
    let kt_ev = K_B * t_kelvin / Q_E;
    eg_ev / 2.0 + kt_ev / 2.0 * (nv / nc).ln()
}

/// Fermi level for n-type semiconductor (eV above valence band edge).
///
/// E_F - E_v = E_g + k_B T ln(n / N_C) ... but measured from E_v:
/// E_F = E_c + k_B T ln(n / N_C) = (E_v + E_g) + k_B T ln(n / N_C)
/// Relative to E_v: E_F - E_v = E_g + k_B T ln(n / N_C).
pub fn fermi_level_n_type_ev(eg_ev: f64, n: f64, nc: f64, t_kelvin: f64) -> f64 {
    let kt_ev = K_B * t_kelvin / Q_E;
    eg_ev + kt_ev * (n / nc).ln()
}

/// Fermi level for p-type semiconductor (eV above valence band edge).
///
/// E_F - E_v = -k_B T ln(p / N_V).
pub fn fermi_level_p_type_ev(p: f64, nv: f64, t_kelvin: f64) -> f64 {
    let kt_ev = K_B * t_kelvin / Q_E;
    -kt_ev * (p / nv).ln()
}

// ---------------------------------------------------------------------------
// Mobility models (Caughey-Thomas)
// ---------------------------------------------------------------------------

/// Caughey-Thomas mobility model parameters.
///
/// μ = μ_min + (μ_max - μ_min) / (1 + (N / N_ref)^α)
#[derive(Clone, Debug)]
pub struct CaugheyThomasMobility {
    /// Minimum mobility (m²/V·s).
    pub mu_min: f64,
    /// Maximum (lattice) mobility (m²/V·s).
    pub mu_max: f64,
    /// Reference doping concentration (1/m³).
    pub n_ref: f64,
    /// Exponent α.
    pub alpha: f64,
}

impl CaugheyThomasMobility {
    /// Create a new Caughey-Thomas model.
    pub fn new(mu_min: f64, mu_max: f64, n_ref: f64, alpha: f64) -> Self {
        Self {
            mu_min,
            mu_max,
            n_ref,
            alpha,
        }
    }

    /// Electron mobility model for silicon at 300 K.
    ///
    /// Parameters from Caughey & Thomas (1967) converted to SI.
    pub fn si_electron() -> Self {
        Self::new(
            0.0065, // 65 cm²/V·s → 0.0065 m²/V·s
            0.1414, // 1414 cm²/V·s
            9.2e22, // 9.2e16 cm⁻³ → 9.2e22 m⁻³
            0.711,
        )
    }

    /// Hole mobility model for silicon at 300 K.
    pub fn si_hole() -> Self {
        Self::new(
            0.00474, // 47.4 cm²/V·s
            0.0471,  // 471 cm²/V·s
            2.23e23, // 2.23e17 cm⁻³ → 2.23e23 m⁻³
            0.719,
        )
    }

    /// Compute mobility at doping concentration `n_doping` (1/m³).
    ///
    /// Returns mobility in m²/V·s.
    pub fn mobility(&self, n_doping: f64) -> f64 {
        self.mu_min + (self.mu_max - self.mu_min) / (1.0 + (n_doping / self.n_ref).powf(self.alpha))
    }
}

/// High-field saturation velocity model.
///
/// μ_eff = μ_low / (1 + (μ_low E / v_sat))
/// where `e_field` is the electric field (V/m) and `v_sat` is saturation
/// velocity (m/s).
pub fn high_field_mobility(mu_low: f64, e_field: f64, v_sat: f64) -> f64 {
    mu_low / (1.0 + mu_low * e_field.abs() / v_sat)
}

/// Drift velocity under electric field (m/s).
pub fn drift_velocity(mobility: f64, e_field: f64) -> f64 {
    mobility * e_field
}

// ---------------------------------------------------------------------------
// Resistivity and conductivity
// ---------------------------------------------------------------------------

/// Electrical conductivity σ = q (n μ_n + p μ_p) (S/m).
pub fn conductivity(n: f64, mu_n: f64, p: f64, mu_p: f64) -> f64 {
    Q_E * (n * mu_n + p * mu_p)
}

/// Resistivity ρ = 1/σ (Ω·m).
pub fn resistivity(n: f64, mu_n: f64, p: f64, mu_p: f64) -> f64 {
    1.0 / conductivity(n, mu_n, p, mu_p)
}

/// Sheet resistance for a thin layer of thickness `t` (Ω/□).
pub fn sheet_resistance(rho: f64, thickness: f64) -> f64 {
    rho / thickness
}

// ---------------------------------------------------------------------------
// Hall effect
// ---------------------------------------------------------------------------

/// Hall coefficient R_H for a predominantly n-type semiconductor (m³/C).
///
/// R_H = -1 / (n q) for electrons.
pub fn hall_coefficient_n_type(n: f64) -> f64 {
    -1.0 / (n * Q_E)
}

/// Hall coefficient R_H for a predominantly p-type semiconductor (m³/C).
///
/// R_H = +1 / (p q) for holes.
pub fn hall_coefficient_p_type(p: f64) -> f64 {
    1.0 / (p * Q_E)
}

/// Hall voltage for a sample of thickness `d` carrying current `i_current`
/// in magnetic field `b_field` (T).
///
/// V_H = R_H * I * B / d.
pub fn hall_voltage(r_h: f64, i_current: f64, b_field: f64, thickness: f64) -> f64 {
    r_h * i_current * b_field / thickness
}

/// Hall mobility μ_H = |R_H| σ (m²/V·s).
pub fn hall_mobility(r_h: f64, sigma: f64) -> f64 {
    r_h.abs() * sigma
}

/// Determine carrier type and concentration from Hall coefficient.
///
/// Returns (carrier_concentration, is_n_type).
pub fn hall_carrier_analysis(r_h: f64) -> (f64, bool) {
    let is_n_type = r_h < 0.0;
    let concentration = 1.0 / (r_h.abs() * Q_E);
    (concentration, is_n_type)
}

// ---------------------------------------------------------------------------
// Debye length
// ---------------------------------------------------------------------------

/// Extrinsic Debye length (m).
///
/// L_D = sqrt(ε_s ε_0 k_B T / (q² n))
/// where `eps_r` is relative permittivity and `n_carrier` is the majority
/// carrier concentration (1/m³).
pub fn debye_length(eps_r: f64, t_kelvin: f64, n_carrier: f64) -> f64 {
    (eps_r * EPSILON_0 * K_B * t_kelvin / (Q_E * Q_E * n_carrier)).sqrt()
}

/// Intrinsic Debye length using intrinsic carrier concentration.
pub fn intrinsic_debye_length(eps_r: f64, t_kelvin: f64, ni: f64) -> f64 {
    debye_length(eps_r, t_kelvin, ni)
}

// ---------------------------------------------------------------------------
// PN junction
// ---------------------------------------------------------------------------

/// PN junction parameters.
#[derive(Clone, Debug)]
pub struct PnJunction {
    /// Donor concentration in n-region (1/m³).
    pub n_d: f64,
    /// Acceptor concentration in p-region (1/m³).
    pub n_a: f64,
    /// Intrinsic carrier concentration (1/m³).
    pub ni: f64,
    /// Relative permittivity of the semiconductor.
    pub eps_r: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Junction area (m²).
    pub area: f64,
}

impl PnJunction {
    /// Create a new PN junction.
    pub fn new(n_d: f64, n_a: f64, ni: f64, eps_r: f64, temperature: f64, area: f64) -> Self {
        Self {
            n_d,
            n_a,
            ni,
            eps_r,
            temperature,
            area,
        }
    }

    /// Typical silicon PN junction at 300 K.
    pub fn silicon_default() -> Self {
        let ni = intrinsic_carrier_si(300.0);
        Self::new(1e22, 1e22, ni, 11.7, 300.0, 1e-6)
    }

    /// Thermal voltage V_T = k_B T / q (V).
    pub fn thermal_voltage(&self) -> f64 {
        K_B * self.temperature / Q_E
    }

    /// Built-in potential V_bi (V).
    ///
    /// V_bi = (k_B T / q) ln(N_A N_D / n_i²)
    pub fn built_in_potential(&self) -> f64 {
        self.thermal_voltage() * (self.n_a * self.n_d / (self.ni * self.ni)).ln()
    }

    /// Depletion width W at reverse bias `v_r` (V, positive for reverse).
    ///
    /// W = sqrt(2 ε_s ε_0 (V_bi + V_R) (1/N_A + 1/N_D) / q)
    pub fn depletion_width(&self, v_r: f64) -> f64 {
        let vbi = self.built_in_potential();
        let eps = self.eps_r * EPSILON_0;
        (2.0 * eps * (vbi + v_r) * (1.0 / self.n_a + 1.0 / self.n_d) / Q_E).sqrt()
    }

    /// Depletion width on the n-side (m).
    pub fn depletion_width_n(&self, v_r: f64) -> f64 {
        let w = self.depletion_width(v_r);
        w * self.n_a / (self.n_a + self.n_d)
    }

    /// Depletion width on the p-side (m).
    pub fn depletion_width_p(&self, v_r: f64) -> f64 {
        let w = self.depletion_width(v_r);
        w * self.n_d / (self.n_a + self.n_d)
    }

    /// Maximum electric field at the junction (V/m).
    pub fn max_electric_field(&self, v_r: f64) -> f64 {
        let w = self.depletion_width(v_r);
        let vbi = self.built_in_potential();
        2.0 * (vbi + v_r) / w
    }

    /// Junction capacitance per unit area (F/m²).
    ///
    /// C_j = ε_s ε_0 / W
    pub fn junction_capacitance_per_area(&self, v_r: f64) -> f64 {
        let eps = self.eps_r * EPSILON_0;
        let w = self.depletion_width(v_r);
        eps / w
    }

    /// Total junction capacitance (F).
    pub fn junction_capacitance(&self, v_r: f64) -> f64 {
        self.junction_capacitance_per_area(v_r) * self.area
    }

    /// Diode saturation current I_0 (A).
    ///
    /// I_0 = A q n_i² (D_n/(L_n N_A) + D_p/(L_p N_D))
    /// Simplified using typical Si diffusion coefficients.
    pub fn saturation_current(&self, d_n: f64, l_n: f64, d_p: f64, l_p: f64) -> f64 {
        self.area * Q_E * self.ni * self.ni * (d_n / (l_n * self.n_a) + d_p / (l_p * self.n_d))
    }

    /// Ideal diode current I(V) = I_0 (exp(V/V_T) - 1) (A).
    pub fn diode_current(&self, voltage: f64, i0: f64) -> f64 {
        i0 * ((voltage / self.thermal_voltage()).exp() - 1.0)
    }

    /// Diode current with ideality factor n.
    pub fn diode_current_ideality(&self, voltage: f64, i0: f64, ideality: f64) -> f64 {
        i0 * ((voltage / (ideality * self.thermal_voltage())).exp() - 1.0)
    }

    /// Breakdown voltage estimate (V) using critical field `e_crit` (V/m).
    ///
    /// V_BD ≈ ε_s ε_0 E_crit² / (2 q N_B) where N_B = min(N_A, N_D).
    pub fn breakdown_voltage(&self, e_crit: f64) -> f64 {
        let n_b = self.n_a.min(self.n_d);
        let eps = self.eps_r * EPSILON_0;
        eps * e_crit * e_crit / (2.0 * Q_E * n_b)
    }
}

// ---------------------------------------------------------------------------
// Schottky barrier
// ---------------------------------------------------------------------------

/// Schottky barrier diode model.
#[derive(Clone, Debug)]
pub struct SchottkyBarrier {
    /// Metal work function (eV).
    pub phi_m_ev: f64,
    /// Semiconductor electron affinity (eV).
    pub chi_s_ev: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Effective Richardson constant (A/m²·K²).
    pub a_star: f64,
    /// Contact area (m²).
    pub area: f64,
    /// Semiconductor relative permittivity.
    pub eps_r: f64,
    /// Doping concentration (1/m³).
    pub n_doping: f64,
}

impl SchottkyBarrier {
    /// Create a new Schottky barrier model.
    pub fn new(
        phi_m_ev: f64,
        chi_s_ev: f64,
        temperature: f64,
        a_star: f64,
        area: f64,
        eps_r: f64,
        n_doping: f64,
    ) -> Self {
        Self {
            phi_m_ev,
            chi_s_ev,
            temperature,
            a_star,
            area,
            eps_r,
            n_doping,
        }
    }

    /// Barrier height φ_B = φ_m - χ_s (eV).
    pub fn barrier_height_ev(&self) -> f64 {
        self.phi_m_ev - self.chi_s_ev
    }

    /// Thermal voltage (V).
    pub fn thermal_voltage(&self) -> f64 {
        K_B * self.temperature / Q_E
    }

    /// Saturation current density J_0 = A* T² exp(-q φ_B / k_B T) (A/m²).
    pub fn saturation_current_density(&self) -> f64 {
        let phi_b_j = self.barrier_height_ev() * Q_E;
        self.a_star
            * self.temperature
            * self.temperature
            * (-phi_b_j / (K_B * self.temperature)).exp()
    }

    /// Saturation current I_0 = J_0 * A (A).
    pub fn saturation_current(&self) -> f64 {
        self.saturation_current_density() * self.area
    }

    /// Forward current at voltage `v` (A).
    ///
    /// I = I_0 (exp(V / V_T) - 1)
    pub fn current(&self, voltage: f64) -> f64 {
        let i0 = self.saturation_current();
        i0 * ((voltage / self.thermal_voltage()).exp() - 1.0)
    }

    /// Schottky depletion width at reverse bias `v_r` (m).
    pub fn depletion_width(&self, v_r: f64) -> f64 {
        let phi_b = self.barrier_height_ev() * Q_E / Q_E; // in V
        let eps = self.eps_r * EPSILON_0;
        (2.0 * eps * (phi_b + v_r) / (Q_E * self.n_doping)).sqrt()
    }

    /// Schottky barrier capacitance per unit area (F/m²).
    pub fn capacitance_per_area(&self, v_r: f64) -> f64 {
        let eps = self.eps_r * EPSILON_0;
        eps / self.depletion_width(v_r)
    }

    /// Image-force barrier lowering Δφ (eV).
    ///
    /// Δφ = sqrt(q E_max / (4π ε_s ε_0))
    pub fn barrier_lowering_ev(&self, e_max: f64) -> f64 {
        let eps = self.eps_r * EPSILON_0;
        (Q_E * e_max / (4.0 * PI * eps)).sqrt() / Q_E
    }
}

// ---------------------------------------------------------------------------
// Photovoltaic: Shockley-Queisser limit
// ---------------------------------------------------------------------------

/// Shockley-Queisser single-junction photovoltaic model.
#[derive(Clone, Debug)]
pub struct ShockleyQueisser {
    /// Band gap (eV).
    pub eg_ev: f64,
    /// Cell temperature (K).
    pub temperature: f64,
    /// Sun temperature (K), default 5778 K.
    pub t_sun: f64,
    /// Concentration factor (1 for 1-sun).
    pub concentration: f64,
}

impl ShockleyQueisser {
    /// Create a new Shockley-Queisser model.
    pub fn new(eg_ev: f64, temperature: f64) -> Self {
        Self {
            eg_ev,
            temperature,
            t_sun: 5778.0,
            concentration: 1.0,
        }
    }

    /// With concentration ratio.
    pub fn with_concentration(mut self, c: f64) -> Self {
        self.concentration = c;
        self
    }

    /// Thermal voltage (V).
    pub fn thermal_voltage(&self) -> f64 {
        K_B * self.temperature / Q_E
    }

    /// Approximate short-circuit current density (A/m²).
    ///
    /// Integrates Planck spectrum above E_g. Uses simplified analytical form:
    /// J_sc ≈ q * (2π / (c² h³)) * (k_B T_sun)³ * Σ_{n≥0} ...
    /// We use the blackbody photon flux approximation.
    pub fn short_circuit_current_density(&self) -> f64 {
        let eg_j = self.eg_ev * Q_E;
        let kt_sun = K_B * self.t_sun;
        let x = eg_j / kt_sun;
        // Approximate photon flux above E_g using analytical integral.
        // Φ ≈ (2π/(c²h³)) (k_B T_sun)³ [x² + 2x + 2] exp(-x)
        let prefactor = 2.0 * PI / (C_LIGHT * C_LIGHT * H_PLANCK * H_PLANCK * H_PLANCK);
        let flux = prefactor * kt_sun * kt_sun * kt_sun * (x * x + 2.0 * x + 2.0) * (-x).exp();
        // Scale by geometric factor (sun solid angle / 2π) ≈ 2.18e-5 for 1 sun
        let geometric = 2.18e-5 * self.concentration;
        Q_E * flux * geometric
    }

    /// Open-circuit voltage (V).
    ///
    /// V_oc ≈ E_g/q - (k_B T / q) ln(J_0 / J_sc)
    /// Simplified: V_oc ≈ E_g/q - k_B T/q * ln(...)
    pub fn open_circuit_voltage(&self) -> f64 {
        let jsc = self.short_circuit_current_density();
        let vt = self.thermal_voltage();
        // Radiative saturation current density
        let j0 = self.radiative_saturation_current();
        if j0 <= 0.0 || jsc <= 0.0 {
            return 0.0;
        }
        vt * (jsc / j0 + 1.0).ln()
    }

    /// Radiative saturation current density J_0 (A/m²).
    pub fn radiative_saturation_current(&self) -> f64 {
        let eg_j = self.eg_ev * Q_E;
        let kt = K_B * self.temperature;
        let x = eg_j / kt;
        let prefactor = 2.0 * PI * Q_E / (C_LIGHT * C_LIGHT * H_PLANCK * H_PLANCK * H_PLANCK);
        prefactor * kt * kt * kt * (x * x + 2.0 * x + 2.0) * (-x).exp()
    }

    /// Fill factor estimate (empirical formula by Green).
    ///
    /// FF ≈ (v_oc - ln(v_oc + 0.72)) / (v_oc + 1)
    /// where v_oc = V_oc / V_T (normalised).
    pub fn fill_factor(&self) -> f64 {
        let voc = self.open_circuit_voltage();
        let vt = self.thermal_voltage();
        let v_norm = voc / vt;
        if v_norm <= 0.0 {
            return 0.0;
        }
        (v_norm - (v_norm + 0.72).ln()) / (v_norm + 1.0)
    }

    /// Maximum power conversion efficiency η.
    ///
    /// η = J_sc * V_oc * FF / P_in
    /// where P_in = σ T_sun⁴ * (geometric factor).
    pub fn efficiency(&self) -> f64 {
        let jsc = self.short_circuit_current_density();
        let voc = self.open_circuit_voltage();
        let ff = self.fill_factor();
        let geometric = 2.18e-5 * self.concentration;
        let p_in = SIGMA_SB * self.t_sun.powi(4) * geometric;
        if p_in <= 0.0 {
            return 0.0;
        }
        jsc * voc * ff / p_in
    }

    /// Maximum power density (W/m²).
    pub fn max_power_density(&self) -> f64 {
        self.short_circuit_current_density() * self.open_circuit_voltage() * self.fill_factor()
    }
}

// ---------------------------------------------------------------------------
// Thermoelectric effects
// ---------------------------------------------------------------------------

/// Thermoelectric material parameters.
#[derive(Clone, Debug)]
pub struct ThermoelectricMaterial {
    /// Seebeck coefficient S (V/K).
    pub seebeck: f64,
    /// Electrical conductivity σ (S/m).
    pub sigma: f64,
    /// Thermal conductivity κ (W/m·K).
    pub kappa: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl ThermoelectricMaterial {
    /// Create a new thermoelectric material.
    pub fn new(seebeck: f64, sigma: f64, kappa: f64, temperature: f64) -> Self {
        Self {
            seebeck,
            sigma,
            kappa,
            temperature,
        }
    }

    /// Typical Bi₂Te₃ at 300 K.
    pub fn bi2te3() -> Self {
        Self::new(200e-6, 1.0e5, 1.5, 300.0)
    }

    /// Thermoelectric figure of merit ZT = S² σ T / κ.
    pub fn zt(&self) -> f64 {
        self.seebeck * self.seebeck * self.sigma * self.temperature / self.kappa
    }

    /// Power factor S² σ (W/m·K²).
    pub fn power_factor(&self) -> f64 {
        self.seebeck * self.seebeck * self.sigma
    }

    /// Peltier coefficient Π = S T (V).
    pub fn peltier_coefficient(&self) -> f64 {
        self.seebeck * self.temperature
    }

    /// Thomson coefficient μ_T = T dS/dT.
    ///
    /// Approximation: given Seebeck at two temperatures.
    pub fn thomson_coefficient(s1: f64, t1: f64, s2: f64, t2: f64) -> f64 {
        let ds_dt = (s2 - s1) / (t2 - t1);
        let t_avg = (t1 + t2) / 2.0;
        t_avg * ds_dt
    }

    /// Seebeck voltage for temperature difference ΔT (V).
    pub fn seebeck_voltage(&self, delta_t: f64) -> f64 {
        self.seebeck * delta_t
    }

    /// Peltier heat flux Q = Π I = S T I (W) for current `i_current` (A).
    pub fn peltier_heat(&self, i_current: f64) -> f64 {
        self.peltier_coefficient() * i_current
    }

    /// Maximum cooling ΔT for a thermoelectric cooler.
    ///
    /// ΔT_max = ZT_c² / 2 * T_c  (simplified).
    /// More precisely: ΔT_max ≈ 0.5 * ZT * T_c for ZT not too large.
    pub fn max_cooling_delta_t(&self, t_cold: f64) -> f64 {
        0.5 * self.zt() * t_cold / (1.0 + 0.5 * self.zt())
    }

    /// Maximum COP of thermoelectric cooler.
    ///
    /// COP_max = (T_c / ΔT) * ((1 + ZT_m)^0.5 - T_h/T_c) / ((1 + ZT_m)^0.5 + 1)
    pub fn max_cop(&self, t_hot: f64, t_cold: f64) -> f64 {
        let zt_m = self.zt();
        let gamma = (1.0 + zt_m).sqrt();
        let dt = t_hot - t_cold;
        if dt <= 0.0 {
            return f64::INFINITY;
        }
        (t_cold / dt) * (gamma - t_hot / t_cold) / (gamma + 1.0)
    }

    /// Thermoelectric generator efficiency.
    ///
    /// η = (ΔT/T_h) * (√(1+ZT) - 1) / (√(1+ZT) + T_c/T_h)
    pub fn generator_efficiency(&self, t_hot: f64, t_cold: f64) -> f64 {
        let dt = t_hot - t_cold;
        if t_hot <= 0.0 {
            return 0.0;
        }
        let carnot = dt / t_hot;
        let zt = self.zt();
        let gamma = (1.0 + zt).sqrt();
        carnot * (gamma - 1.0) / (gamma + t_cold / t_hot)
    }

    /// Electrical resistivity ρ = 1/σ (Ω·m).
    pub fn resistivity(&self) -> f64 {
        1.0 / self.sigma
    }

    /// Lorenz number estimate L = κ / (σ T) (W·Ω/K²).
    pub fn lorenz_number(&self) -> f64 {
        self.kappa / (self.sigma * self.temperature)
    }
}

/// Thermoelectric module (n-p couple) power output.
///
/// Returns (power, voltage, current) for given ΔT and load resistance.
pub fn thermoelectric_module_power(
    seebeck_np: f64,
    r_internal: f64,
    r_load: f64,
    delta_t: f64,
) -> (f64, f64, f64) {
    let v_oc = seebeck_np * delta_t;
    let current = v_oc / (r_internal + r_load);
    let voltage = current * r_load;
    let power = current * voltage;
    (power, voltage, current)
}

/// Maximum thermoelectric module power (at matched load R_L = R_int).
pub fn thermoelectric_max_power(seebeck_np: f64, r_internal: f64, delta_t: f64) -> f64 {
    let v_oc = seebeck_np * delta_t;
    v_oc * v_oc / (4.0 * r_internal)
}

// ---------------------------------------------------------------------------
// Additional utility functions
// ---------------------------------------------------------------------------

/// Einstein relation: D = μ k_B T / q (m²/s).
pub fn einstein_diffusivity(mobility: f64, t_kelvin: f64) -> f64 {
    mobility * K_B * t_kelvin / Q_E
}

/// Diffusion length L = √(D τ) (m).
pub fn diffusion_length(diffusivity: f64, lifetime: f64) -> f64 {
    (diffusivity * lifetime).sqrt()
}

/// Generation rate for photon absorption (1/m³·s).
///
/// G = α Φ exp(-α x) where α is absorption coefficient (1/m),
/// Φ is photon flux (1/m²·s), and x is depth (m).
pub fn generation_rate(alpha: f64, photon_flux: f64, depth: f64) -> f64 {
    alpha * photon_flux * (-alpha * depth).exp()
}

/// Recombination rate (Shockley-Read-Hall) for a single trap level.
///
/// R_SRH = (n p - n_i²) / (τ_p (n + n_1) + τ_n (p + p_1))
pub fn recombination_srh(n: f64, p: f64, ni: f64, tau_n: f64, tau_p: f64, n1: f64, p1: f64) -> f64 {
    let numerator = n * p - ni * ni;
    let denominator = tau_p * (n + n1) + tau_n * (p + p1);
    if denominator.abs() < 1e-30 {
        return 0.0;
    }
    numerator / denominator
}

/// Auger recombination rate (1/m³·s).
///
/// R_Auger = C_n n² p + C_p n p²
pub fn recombination_auger(n: f64, p: f64, c_n: f64, c_p: f64) -> f64 {
    c_n * n * n * p + c_p * n * p * p
}

/// Radiative recombination rate (1/m³·s).
///
/// R_rad = B (n p - n_i²)
pub fn recombination_radiative(n: f64, p: f64, ni: f64, b_coeff: f64) -> f64 {
    b_coeff * (n * p - ni * ni)
}

/// Total recombination lifetime combining SRH, Auger, and radiative.
///
/// 1/τ_eff = 1/τ_SRH + 1/τ_Auger + 1/τ_rad
pub fn effective_lifetime(tau_srh: f64, tau_auger: f64, tau_rad: f64) -> f64 {
    1.0 / (1.0 / tau_srh + 1.0 / tau_auger + 1.0 / tau_rad)
}

/// Absorption coefficient for direct band gap (simplified parabolic model).
///
/// α = A √(E - E_g)   for E > E_g
pub fn absorption_direct(photon_energy_ev: f64, eg_ev: f64, a_coeff: f64) -> f64 {
    if photon_energy_ev <= eg_ev {
        0.0
    } else {
        a_coeff * (photon_energy_ev - eg_ev).sqrt()
    }
}

/// Absorption coefficient for indirect band gap (simplified model).
///
/// α = B (E - E_g ± E_phonon)² / E  (phonon-assisted transitions).
pub fn absorption_indirect(
    photon_energy_ev: f64,
    eg_ev: f64,
    e_phonon_ev: f64,
    b_coeff: f64,
) -> f64 {
    // Phonon absorption process
    let e_abs = photon_energy_ev - eg_ev + e_phonon_ev;
    // Phonon emission process
    let e_emi = photon_energy_ev - eg_ev - e_phonon_ev;
    let mut alpha = 0.0;
    if e_abs > 0.0 {
        alpha += b_coeff * e_abs * e_abs / photon_energy_ev;
    }
    if e_emi > 0.0 {
        alpha += b_coeff * e_emi * e_emi / photon_energy_ev;
    }
    alpha
}

/// Built-in potential for a heterojunction (eV).
///
/// V_bi = (χ₁ - χ₂) + (E_g2 - E_g1) + k_B T/q ln(N_A N_D / (n_i1 n_i2))
pub fn heterojunction_built_in(
    chi1_ev: f64,
    chi2_ev: f64,
    eg1_ev: f64,
    eg2_ev: f64,
    n_a: f64,
    n_d: f64,
    ni1: f64,
    ni2: f64,
    t_kelvin: f64,
) -> f64 {
    let kt_ev = K_B * t_kelvin / Q_E;
    (chi1_ev - chi2_ev) + (eg2_ev - eg1_ev) + kt_ev * (n_a * n_d / (ni1 * ni2)).ln()
}

/// Metal-semiconductor contact: Ohmic or Schottky determination.
///
/// Returns true if the contact is Schottky (rectifying) for n-type.
/// Schottky condition: φ_m > χ_s for n-type.
pub fn is_schottky_n_type(phi_m_ev: f64, chi_s_ev: f64) -> bool {
    phi_m_ev > chi_s_ev
}

/// Metal-semiconductor contact: Ohmic or Schottky for p-type.
///
/// Schottky condition: φ_m < χ_s + E_g for p-type.
pub fn is_schottky_p_type(phi_m_ev: f64, chi_s_ev: f64, eg_ev: f64) -> bool {
    phi_m_ev < chi_s_ev + eg_ev
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // 1. Band gap: Si at 300 K ≈ 1.12 eV
    #[test]
    fn test_silicon_band_gap_300k() {
        let bg = BandGapModel::silicon();
        let eg = bg.band_gap_ev(300.0);
        assert!(
            (eg - 1.12).abs() < 0.02,
            "Si band gap at 300K should be ~1.12 eV, got {:.6}",
            eg
        );
    }

    // 2. Band gap: GaAs at 0 K
    #[test]
    fn test_gaas_band_gap_0k() {
        let bg = BandGapModel::gaas();
        let eg = bg.band_gap_ev(0.0);
        assert!(
            (eg - 1.519).abs() < TOL,
            "GaAs band gap at 0K should be 1.519 eV, got {:.6}",
            eg
        );
    }

    // 3. Band gap type identification
    #[test]
    fn test_band_gap_type() {
        let si = BandGapModel::silicon();
        let gaas = BandGapModel::gaas();
        assert!(!si.is_direct());
        assert!(gaas.is_direct());
    }

    // 4. Effective DOS positive and increases with temperature
    #[test]
    fn test_effective_dos_temperature_dependence() {
        let nc_300 = effective_dos_conduction(1.08, 300.0);
        let nc_400 = effective_dos_conduction(1.08, 400.0);
        assert!(nc_300 > 0.0);
        assert!(nc_400 > nc_300, "DOS should increase with temperature");
    }

    // 5. Intrinsic carrier concentration: Si at 300 K ≈ 1e16 m⁻³ (1e10 cm⁻³)
    #[test]
    fn test_intrinsic_carrier_si() {
        let ni = intrinsic_carrier_si(300.0);
        // ni for Si at 300 K is about 1.0e10 cm⁻³ = 1.0e16 m⁻³
        assert!(
            ni > 1e14 && ni < 1e18,
            "Si ni at 300K should be ~1e16 m⁻³, got {:.6e}",
            ni
        );
    }

    // 6. Extrinsic n-type: n ≈ N_D, p << n
    #[test]
    fn test_extrinsic_n_type() {
        let ni = 1e16;
        let n_d = 1e22;
        let (n, p) = extrinsic_n_type(n_d, ni);
        assert!((n - n_d).abs() < 1.0);
        assert!(p < n * 1e-6, "Minority carrier p should be << n");
    }

    // 7. Extrinsic p-type: p ≈ N_A, n << p
    #[test]
    fn test_extrinsic_p_type() {
        let ni = 1e16;
        let n_a = 1e22;
        let (n, p) = extrinsic_p_type(n_a, ni);
        assert!((p - n_a).abs() < 1.0);
        assert!(n < p * 1e-6, "Minority carrier n should be << p");
    }

    // 8. Caughey-Thomas: low doping → high mobility, high doping → low mobility
    #[test]
    fn test_caughey_thomas_mobility() {
        let ct = CaugheyThomasMobility::si_electron();
        let mu_low = ct.mobility(1e18); // low doping
        let mu_high = ct.mobility(1e26); // very high doping
        assert!(
            mu_low > mu_high,
            "Low doping mobility {:.6} should exceed high doping {:.6}",
            mu_low,
            mu_high
        );
        assert!(mu_low > 0.0);
        assert!(mu_high > ct.mu_min - 1e-10);
    }

    // 9. Conductivity is positive
    #[test]
    fn test_conductivity_positive() {
        let sigma = conductivity(1e22, 0.14, 1e10, 0.05);
        assert!(
            sigma > 0.0,
            "Conductivity must be positive, got {:.6}",
            sigma
        );
    }

    // 10. Resistivity = 1/σ
    #[test]
    fn test_resistivity_inverse_conductivity() {
        let sigma = conductivity(1e22, 0.14, 1e10, 0.05);
        let rho = resistivity(1e22, 0.14, 1e10, 0.05);
        assert!(
            (rho - 1.0 / sigma).abs() < 1e-20,
            "Resistivity should be inverse of conductivity"
        );
    }

    // 11. Hall coefficient sign: n-type → negative
    #[test]
    fn test_hall_coefficient_sign() {
        let rh_n = hall_coefficient_n_type(1e22);
        let rh_p = hall_coefficient_p_type(1e22);
        assert!(rh_n < 0.0, "n-type Hall coefficient should be negative");
        assert!(rh_p > 0.0, "p-type Hall coefficient should be positive");
    }

    // 12. Hall carrier analysis roundtrip
    #[test]
    fn test_hall_carrier_analysis() {
        let n = 1e22;
        let rh = hall_coefficient_n_type(n);
        let (n_recovered, is_n) = hall_carrier_analysis(rh);
        assert!(is_n, "Should identify as n-type");
        assert!(
            (n_recovered - n).abs() / n < 1e-6,
            "Recovered concentration should match"
        );
    }

    // 13. Debye length positive and decreases with higher carrier concentration
    #[test]
    fn test_debye_length() {
        let ld1 = debye_length(11.7, 300.0, 1e20);
        let ld2 = debye_length(11.7, 300.0, 1e22);
        assert!(ld1 > 0.0);
        assert!(
            ld1 > ld2,
            "Higher carrier → shorter Debye length: {:.6e} vs {:.6e}",
            ld1,
            ld2
        );
    }

    // 14. PN junction built-in potential > 0
    #[test]
    fn test_pn_junction_built_in() {
        let pn = PnJunction::silicon_default();
        let vbi = pn.built_in_potential();
        assert!(
            vbi > 0.0 && vbi < 1.5,
            "Si PN junction Vbi should be 0.5-0.9 V, got {:.6}",
            vbi
        );
    }

    // 15. PN junction depletion width increases with reverse bias
    #[test]
    fn test_pn_depletion_width_bias() {
        let pn = PnJunction::silicon_default();
        let w0 = pn.depletion_width(0.0);
        let w5 = pn.depletion_width(5.0);
        assert!(
            w5 > w0,
            "Depletion width should increase with reverse bias: {:.6e} > {:.6e}",
            w5,
            w0
        );
    }

    // 16. PN junction capacitance decreases with reverse bias
    #[test]
    fn test_pn_capacitance_bias() {
        let pn = PnJunction::silicon_default();
        let c0 = pn.junction_capacitance(0.0);
        let c5 = pn.junction_capacitance(5.0);
        assert!(c0 > c5, "Capacitance should decrease with reverse bias");
    }

    // 17. Diode current: forward > 0, reverse ≈ -I_0
    #[test]
    fn test_diode_current() {
        let pn = PnJunction::silicon_default();
        let i0 = 1e-12;
        let i_fwd = pn.diode_current(0.6, i0);
        let i_rev = pn.diode_current(-5.0, i0);
        assert!(i_fwd > 0.0, "Forward current must be positive");
        assert!(
            (i_rev + i0).abs() < i0 * 0.01,
            "Reverse current should approach -I0"
        );
    }

    // 18. Schottky barrier height
    #[test]
    fn test_schottky_barrier_height() {
        let sb = SchottkyBarrier::new(4.8, 4.05, 300.0, 1.12e6, 1e-6, 11.7, 1e22);
        let phi_b = sb.barrier_height_ev();
        assert!(
            (phi_b - 0.75).abs() < 0.01,
            "Barrier height should be ~0.75 eV, got {:.6}",
            phi_b
        );
    }

    // 19. Schottky: forward current increases exponentially
    #[test]
    fn test_schottky_current_exponential() {
        let sb = SchottkyBarrier::new(4.8, 4.05, 300.0, 1.12e6, 1e-6, 11.7, 1e22);
        let i1 = sb.current(0.3);
        let i2 = sb.current(0.4);
        assert!(
            i2 > i1 * 10.0,
            "100 mV increase should give ~50x more current"
        );
    }

    // 20. Shockley-Queisser: Jsc > 0 for reasonable band gap
    #[test]
    fn test_sq_short_circuit_current() {
        let sq = ShockleyQueisser::new(1.4, 300.0);
        let jsc = sq.short_circuit_current_density();
        assert!(jsc > 0.0, "Jsc must be positive, got {:.6e}", jsc);
    }

    // 21. Shockley-Queisser: efficiency between 0 and 1
    #[test]
    fn test_sq_efficiency_bounds() {
        let sq = ShockleyQueisser::new(1.34, 300.0);
        let eta = sq.efficiency();
        assert!(
            eta > 0.0 && eta < 1.0,
            "Efficiency should be between 0 and 1, got {:.6}",
            eta
        );
    }

    // 22. Thermoelectric ZT for Bi2Te3 ≈ 0.8
    #[test]
    fn test_thermoelectric_zt() {
        let te = ThermoelectricMaterial::bi2te3();
        let zt = te.zt();
        assert!(
            zt > 0.5 && zt < 1.5,
            "Bi2Te3 ZT at 300K should be ~0.8, got {:.6}",
            zt
        );
    }

    // 23. Seebeck voltage proportional to ΔT
    #[test]
    fn test_seebeck_voltage_linear() {
        let te = ThermoelectricMaterial::bi2te3();
        let v1 = te.seebeck_voltage(10.0);
        let v2 = te.seebeck_voltage(20.0);
        assert!(
            (v2 - 2.0 * v1).abs() < 1e-15,
            "Seebeck voltage should be linear in ΔT"
        );
    }

    // 24. Peltier coefficient Π = S T
    #[test]
    fn test_peltier_coefficient() {
        let te = ThermoelectricMaterial::bi2te3();
        let pi = te.peltier_coefficient();
        let expected = te.seebeck * te.temperature;
        assert!(
            (pi - expected).abs() < 1e-15,
            "Peltier coeff should equal S*T"
        );
    }

    // 25. Einstein relation: D = μ k_B T / q
    #[test]
    fn test_einstein_relation() {
        let mu = 0.14; // electron mobility in Si (m²/V·s)
        let t = 300.0;
        let d = einstein_diffusivity(mu, t);
        let expected = mu * K_B * t / Q_E;
        assert!(
            (d - expected).abs() < 1e-15,
            "Einstein relation mismatch: got {:.6e}, expected {:.6e}",
            d,
            expected
        );
    }

    // 26. Diffusion length positive
    #[test]
    fn test_diffusion_length() {
        let d = 36e-4; // typical Si electron diffusivity (m²/s)
        let tau = 1e-6; // lifetime (s)
        let l = diffusion_length(d, tau);
        assert!(l > 0.0, "Diffusion length must be positive");
        let expected = (d * tau).sqrt();
        assert!((l - expected).abs() < 1e-15, "Diffusion length mismatch");
    }

    // 27. Thermoelectric generator efficiency bounded by Carnot
    #[test]
    fn test_teg_efficiency_vs_carnot() {
        let te = ThermoelectricMaterial::bi2te3();
        let t_hot = 500.0;
        let t_cold = 300.0;
        let eta = te.generator_efficiency(t_hot, t_cold);
        let carnot = (t_hot - t_cold) / t_hot;
        assert!(
            eta > 0.0 && eta < carnot,
            "TEG efficiency {:.6} should be between 0 and Carnot {:.6}",
            eta,
            carnot
        );
    }

    // 28. Thermoelectric max power at matched load
    #[test]
    fn test_thermoelectric_max_power() {
        let s_np = 400e-6; // combined Seebeck (V/K)
        let r_int = 1.0; // internal resistance (Ω)
        let dt = 100.0;
        let p_max = thermoelectric_max_power(s_np, r_int, dt);
        let expected = (s_np * dt).powi(2) / (4.0 * r_int);
        assert!((p_max - expected).abs() < 1e-15, "Max power mismatch");
    }

    // 29. SRH recombination: equilibrium → R = 0
    #[test]
    fn test_srh_equilibrium() {
        let ni = 1e16;
        let r = recombination_srh(ni, ni, ni, 1e-6, 1e-6, ni, ni);
        assert!(
            r.abs() < 1e-6,
            "SRH recombination should be ~0 at equilibrium, got {:.6e}",
            r
        );
    }

    // 30. Compensated doping: N_D = N_A → intrinsic behavior
    #[test]
    fn test_compensated_doping() {
        let ni = 1e16;
        let n_d = 1e22;
        // Slight excess n-type
        let (n, p) = extrinsic_compensated(n_d, n_d - 1e18, ni);
        assert!(n > p, "Slight n-type excess should give n > p");
    }
}
