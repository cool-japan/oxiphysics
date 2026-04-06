// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aerospace structural material models.
//!
//! Covers temperature-dependent titanium alloys (Ti-6Al-4V), nickel
//! superalloys (Inconel 718 creep), carbon-fibre reinforced polymer (CFRP)
//! with classical laminate theory, ceramic matrix composites (SiC/SiC CMC),
//! thermal barrier coatings (7YSZ zirconia), ablative materials (char
//! formation), metal foams (Gibson–Ashby), shape-memory alloys (NiTi
//! Nitinol), high-entropy alloys (HEA), and reentry vehicle thermal
//! protection systems (TPS).
//!
//! All computations use plain `f64` — no external linear-algebra dependency.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 1. Titanium alloy Ti-6Al-4V — temperature-dependent properties
// ---------------------------------------------------------------------------

/// Temperature-dependent mechanical and thermal properties of Ti-6Al-4V.
///
/// Data fit to MMPDS/MIL-HDBK-5 curves valid from room temperature to
/// ~850 °C.  Above the beta-transus (~995 °C) the alloy softens rapidly;
/// this model clamps at 950 °C.
#[derive(Debug, Clone)]
pub struct Ti6Al4V {
    /// Reference temperature for property evaluation (°C).
    pub temperature_c: f64,
}

impl Ti6Al4V {
    /// Construct with a given temperature in °C.
    pub fn new(temperature_c: f64) -> Self {
        Self {
            temperature_c: temperature_c.min(950.0),
        }
    }

    /// Young's modulus (GPa) as a function of temperature.
    ///
    /// Fitted linear-quadratic curve through MMPDS data points.
    pub fn youngs_modulus_gpa(&self) -> f64 {
        let t = self.temperature_c;
        // E(T) ≈ 113.8 − 0.051·T − 1.2e-5·T² GPa
        113.8 - 0.051 * t - 1.2e-5 * t * t
    }

    /// Yield strength (MPa, 0.2 % proof) as a function of temperature.
    pub fn yield_strength_mpa(&self) -> f64 {
        let t = self.temperature_c;
        // σ_y(T) ≈ 950 − 0.55·T GPa (decreasing above RT)
        (950.0 - 0.55 * t).max(100.0)
    }

    /// Ultimate tensile strength (MPa).
    pub fn uts_mpa(&self) -> f64 {
        self.yield_strength_mpa() * 1.10 // ~10 % above yield across temperature range
    }

    /// Density (kg/m³) — effectively temperature-independent for engineering
    /// purposes.
    pub fn density_kg_m3(&self) -> f64 {
        4430.0
    }

    /// Thermal conductivity (W/m·K).
    pub fn thermal_conductivity_w_mk(&self) -> f64 {
        let t = self.temperature_c;
        // k(T) ≈ 6.7 + 0.0103·T
        6.7 + 0.0103 * t
    }

    /// Specific heat capacity (J/kg·K).
    pub fn specific_heat_j_kgk(&self) -> f64 {
        let t = self.temperature_c;
        // cp(T) ≈ 526 + 0.13·T J/(kg·K)
        526.0 + 0.13 * t
    }

    /// Coefficient of thermal expansion (µm/m·K, i.e. 10⁻⁶ /K).
    pub fn cte_per_k(&self) -> f64 {
        let t = self.temperature_c;
        (8.6 + 0.0014 * t) * 1e-6
    }

    /// Thermal diffusivity α = k / (ρ·cₚ) in m²/s.
    pub fn thermal_diffusivity_m2_s(&self) -> f64 {
        let k = self.thermal_conductivity_w_mk();
        let rho = self.density_kg_m3();
        let cp = self.specific_heat_j_kgk();
        k / (rho * cp)
    }

    /// Poisson's ratio (essentially constant with temperature).
    pub fn poisson_ratio(&self) -> f64 {
        0.342
    }
}

// ---------------------------------------------------------------------------
// 2. Nickel superalloy Inconel 718 — creep model
// ---------------------------------------------------------------------------

/// Inconel 718 material properties and Norton power-law creep model.
///
/// Creep data reference: Ludwigson & Berger (1969); temperature range
/// 540–760 °C.
#[derive(Debug, Clone)]
pub struct Inconel718 {
    /// Service temperature (°C).
    pub temperature_c: f64,
}

impl Inconel718 {
    /// Create with service temperature.
    pub fn new(temperature_c: f64) -> Self {
        Self { temperature_c }
    }

    /// Young's modulus (GPa).
    pub fn youngs_modulus_gpa(&self) -> f64 {
        let t = self.temperature_c;
        200.0 - 0.08 * t
    }

    /// 0.2 % yield strength (MPa).
    pub fn yield_strength_mpa(&self) -> f64 {
        let t = self.temperature_c;
        (1185.0 - 0.85 * t).max(200.0)
    }

    /// Density (kg/m³).
    pub fn density_kg_m3(&self) -> f64 {
        8190.0
    }

    /// Norton power-law creep rate (s⁻¹) at a given stress (MPa).
    ///
    /// `ε̇_cr = A · σⁿ · exp(−Q / RT)`
    ///
    /// Constants fit to Inconel 718 creep data:
    /// - A = 3.0×10⁻²⁸ (MPa⁻ⁿ · s⁻¹)
    /// - n = 4.8
    /// - Q = 310 kJ/mol (activation energy)
    pub fn norton_creep_rate_per_s(&self, stress_mpa: f64) -> f64 {
        const A: f64 = 3.0e-28;
        const N: f64 = 4.8;
        const Q: f64 = 310_000.0; // J/mol
        const R: f64 = 8.314; // J/(mol·K)
        let t_k = self.temperature_c + 273.15;
        A * stress_mpa.powf(N) * (-Q / (R * t_k)).exp()
    }

    /// Creep rupture life estimate (hours) using the Larson–Miller parameter.
    ///
    /// LMP = T(K) × (log₁₀(t_r) + C) where C ≈ 20 for Inconel 718.
    /// This is the inverse: given LMP, solve for t_r in hours.
    pub fn rupture_life_hours(&self, stress_mpa: f64) -> f64 {
        const C: f64 = 20.0;
        // LMP master curve (simplified fit):  LMP = 29000 - 22·σ  (σ in MPa)
        let lmp = 29_000.0 - 22.0 * stress_mpa;
        let t_k = self.temperature_c + 273.15;
        let log_tr = lmp / t_k - C;
        10.0_f64.powf(log_tr)
    }

    /// Thermal conductivity (W/m·K).
    pub fn thermal_conductivity_w_mk(&self) -> f64 {
        let t = self.temperature_c;
        11.4 + 0.013 * t
    }
}

// ---------------------------------------------------------------------------
// 3. CFRP — Classical Laminate Theory (CLT)
// ---------------------------------------------------------------------------

/// A single unidirectional CFRP ply for CLT analysis.
#[derive(Debug, Clone, Copy)]
pub struct CfrpPly {
    /// Ply orientation angle (degrees, measured from laminate x-axis).
    pub angle_deg: f64,
    /// Ply thickness (mm).
    pub thickness_mm: f64,
    /// Longitudinal modulus E₁ (GPa).
    pub e1_gpa: f64,
    /// Transverse modulus E₂ (GPa).
    pub e2_gpa: f64,
    /// In-plane shear modulus G₁₂ (GPa).
    pub g12_gpa: f64,
    /// Major Poisson's ratio ν₁₂.
    pub nu12: f64,
}

impl CfrpPly {
    /// Standard IM7/8552 carbon/epoxy ply.
    pub fn im7_8552(angle_deg: f64) -> Self {
        Self {
            angle_deg,
            thickness_mm: 0.125,
            e1_gpa: 161.0,
            e2_gpa: 11.38,
            g12_gpa: 5.17,
            nu12: 0.32,
        }
    }

    /// Reduced stiffness components Q₁₁, Q₂₂, Q₁₂, Q₆₆ in the ply
    /// material axes (GPa).
    pub fn reduced_stiffness_gpa(&self) -> [f64; 4] {
        let e1 = self.e1_gpa;
        let e2 = self.e2_gpa;
        let g12 = self.g12_gpa;
        let nu12 = self.nu12;
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        let q11 = e1 / denom;
        let q22 = e2 / denom;
        let q12 = nu12 * e2 / denom;
        let q66 = g12;
        [q11, q22, q12, q66]
    }

    /// Transformed reduced stiffness Q̄ components \[Q̄₁₁, Q̄₂₂, Q̄₁₂, Q̄₆₆,
    /// Q̄₁₆, Q̄₂₆] (GPa) in the laminate coordinate frame.
    pub fn transformed_stiffness_gpa(&self) -> [f64; 6] {
        let [q11, q22, q12, q66] = self.reduced_stiffness_gpa();
        let theta = self.angle_deg.to_radians();
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let c2s2 = c2 * s2;
        let q11b = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * s4;
        let q22b = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * c4;
        let q12b = (q11 + q22 - 4.0 * q66) * c2s2 + q12 * (c4 + s4);
        let q66b = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2s2 + q66 * (c2 - s2).powi(2);
        let q16b = (q11 - q12 - 2.0 * q66) * c * c2 * s - (q22 - q12 - 2.0 * q66) * s * s2 * c;
        let q26b = (q11 - q12 - 2.0 * q66) * s * s2 * c - (q22 - q12 - 2.0 * q66) * c * c2 * s;
        [q11b, q22b, q12b, q66b, q16b, q26b]
    }
}

/// Classical laminate theory analysis result for a symmetric balanced
/// laminate.
#[derive(Debug, Clone)]
pub struct CltResult {
    /// In-plane stiffness matrix A (N/mm): \[A₁₁, A₂₂, A₁₂, A₆₆, A₁₆, A₂₆\].
    pub a_matrix: [f64; 6],
    /// Total laminate thickness (mm).
    pub total_thickness_mm: f64,
    /// Effective in-plane modulus Ex (GPa).
    pub ex_gpa: f64,
    /// Effective in-plane modulus Ey (GPa).
    pub ey_gpa: f64,
    /// Effective in-plane shear modulus Gxy (GPa).
    pub gxy_gpa: f64,
}

/// Compute CLT A-matrix and effective properties for a list of plies.
pub fn clt_analysis(plies: &[CfrpPly]) -> CltResult {
    let mut a = [0.0f64; 6];
    let mut total_t = 0.0f64;
    for ply in plies {
        let q_bar = ply.transformed_stiffness_gpa();
        let t = ply.thickness_mm;
        total_t += t;
        for i in 0..6 {
            a[i] += q_bar[i] * t; // A_ij = Σ Q̄_ij * t_k (GPa·mm = kN/mm effectively)
        }
    }
    // Convert A from GPa·mm to N/mm by ×1000
    let a_nmm: [f64; 6] = a.map(|v| v * 1000.0);
    let h = total_t;
    // Ex = (A11*A22 - A12²) / (A22*h)
    let det = a_nmm[0] * a_nmm[1] - a_nmm[2] * a_nmm[2];
    let ex = if a_nmm[1] * h > 1e-12 {
        det / (a_nmm[1] * h) / 1000.0
    } else {
        0.0
    }; // GPa
    let ey = if a_nmm[0] * h > 1e-12 {
        det / (a_nmm[0] * h) / 1000.0
    } else {
        0.0
    };
    let gxy = a_nmm[3] / h / 1000.0;
    CltResult {
        a_matrix: a_nmm,
        total_thickness_mm: h,
        ex_gpa: ex,
        ey_gpa: ey,
        gxy_gpa: gxy,
    }
}

// ---------------------------------------------------------------------------
// 4. Ceramic matrix composite — SiC/SiC
// ---------------------------------------------------------------------------

/// Mechanical and thermophysical properties of a SiC/SiC 2-D woven CMC.
///
/// Property ranges are representative of Hi-Nicalon Type S / CVI-SiC matrix
/// with BN interphase, as reported in UEET programme data.
#[derive(Debug, Clone)]
pub struct SicSicCmc {
    /// Service temperature (°C).
    pub temperature_c: f64,
    /// Fibre volume fraction (0–1).
    pub fibre_volume_fraction: f64,
}

impl SicSicCmc {
    /// Construct with temperature and fibre volume fraction.
    pub fn new(temperature_c: f64, fibre_volume_fraction: f64) -> Self {
        Self {
            temperature_c,
            fibre_volume_fraction: fibre_volume_fraction.clamp(0.35, 0.55),
        }
    }

    /// In-plane Young's modulus (GPa) using rule-of-mixtures.
    pub fn youngs_modulus_gpa(&self) -> f64 {
        let vf = self.fibre_volume_fraction;
        let e_fibre = 380.0; // Hi-Nicalon Type S, GPa
        let e_matrix = 350.0; // CVI SiC matrix, GPa
        // For a 2-D woven fabric, modulus ≈ 0.5·E_longitudinal (average over ±0°/90°)
        0.5 * (e_fibre * vf + e_matrix * (1.0 - vf))
    }

    /// Proportional limit stress (MPa) — onset of matrix cracking.
    pub fn proportional_limit_mpa(&self) -> f64 {
        150.0 - 0.08 * self.temperature_c
    }

    /// Ultimate tensile strength (MPa).
    pub fn uts_mpa(&self) -> f64 {
        (230.0 - 0.10 * self.temperature_c).max(50.0)
    }

    /// Interlaminar shear strength (MPa).
    pub fn ilss_mpa(&self) -> f64 {
        40.0 - 0.02 * self.temperature_c
    }

    /// Density (kg/m³).
    pub fn density_kg_m3(&self) -> f64 {
        2700.0
    }

    /// Thermal conductivity (W/m·K).
    pub fn thermal_conductivity_w_mk(&self) -> f64 {
        let t = self.temperature_c;
        // Decreasing with temperature (phonon scattering)
        (18.0 * (1.0 - t / 1600.0)).max(3.0)
    }

    /// Maximum use temperature for load-bearing applications (°C).
    pub fn max_use_temperature_c(&self) -> f64 {
        1200.0
    }
}

// ---------------------------------------------------------------------------
// 5. Thermal barrier coating (TBC) — 7YSZ zirconia
// ---------------------------------------------------------------------------

/// Thermal barrier coating properties for yttria-stabilised zirconia (7YSZ).
///
/// Represents an EB-PVD columnar TBC as used on high-pressure turbine blades.
#[derive(Debug, Clone)]
pub struct TbcYsz {
    /// Coating thickness (µm).
    pub thickness_um: f64,
    /// Service gas-side surface temperature (°C).
    pub surface_temperature_c: f64,
    /// Substrate (bond coat) temperature (°C).
    pub substrate_temperature_c: f64,
}

impl TbcYsz {
    /// Create a TBC descriptor.
    pub fn new(
        thickness_um: f64,
        surface_temperature_c: f64,
        substrate_temperature_c: f64,
    ) -> Self {
        Self {
            thickness_um,
            surface_temperature_c,
            substrate_temperature_c,
        }
    }

    /// Mean temperature of the coating (°C).
    pub fn mean_temperature_c(&self) -> f64 {
        (self.surface_temperature_c + self.substrate_temperature_c) * 0.5
    }

    /// Thermal conductivity of 7YSZ (W/m·K).
    ///
    /// EB-PVD coatings have k ≈ 1.5–2.0 W/m·K; sintering reduces k above
    /// ~1200 °C.
    pub fn thermal_conductivity_w_mk(&self) -> f64 {
        let t = self.mean_temperature_c();
        if t < 1200.0 {
            1.9 - 2.0e-4 * t
        } else {
            // Sintering-induced densification increases k
            1.7 + 1.5e-4 * (t - 1200.0)
        }
    }

    /// Temperature drop across the TBC (°C).
    ///
    /// ΔT = q'' · t / k, where q'' is the heat flux (W/m²).
    pub fn temperature_drop_c(&self, heat_flux_w_m2: f64) -> f64 {
        let k = self.thermal_conductivity_w_mk();
        let t_m = self.thickness_um * 1e-6; // convert µm to m
        heat_flux_w_m2 * t_m / k
    }

    /// Thermal cycling life estimate (cycles) using a simplified Paris-law
    /// spallation model.
    ///
    /// N_f ≈ C / (ΔT_cycle)^m, with C = 5×10⁶, m = 2.5 (calibrated to
    /// standard furnace-cycle data for 120 µm EB-PVD TBC).
    pub fn spallation_life_cycles(&self, delta_t_cycle_c: f64) -> f64 {
        const C: f64 = 5.0e6;
        const M: f64 = 2.5;
        C / delta_t_cycle_c.powf(M)
    }

    /// Coefficient of thermal expansion mismatch strain with Ni superalloy
    /// substrate (7YSZ α ≈ 11 µm/m·K vs substrate α ≈ 14 µm/m·K).
    pub fn cte_mismatch_strain(&self, delta_t_c: f64) -> f64 {
        let alpha_ysz = 11.0e-6;
        let alpha_substrate = 14.0e-6;
        (alpha_substrate - alpha_ysz) * delta_t_c
    }
}

// ---------------------------------------------------------------------------
// 6. Ablative material — char formation and recession
// ---------------------------------------------------------------------------

/// Ablative heat shield material descriptor.
///
/// Models phenolic-impregnated carbon ablator (PICA-like) charring and
/// surface recession during atmospheric reentry.
#[derive(Debug, Clone)]
pub struct AblativeMaterial {
    /// Virgin material density (kg/m³).
    pub virgin_density_kg_m3: f64,
    /// Char density (kg/m³).
    pub char_density_kg_m3: f64,
    /// Pyrolysis onset temperature (°C).
    pub pyrolysis_onset_c: f64,
    /// Pyrolysis completion temperature (°C).
    pub pyrolysis_complete_c: f64,
    /// Effective heat of pyrolysis (kJ/kg).
    pub heat_of_pyrolysis_kj_kg: f64,
    /// Char specific heat (J/kg·K).
    pub char_cp_j_kgk: f64,
    /// Virgin specific heat (J/kg·K).
    pub virgin_cp_j_kgk: f64,
}

impl AblativeMaterial {
    /// PICA-representative ablator (Phenolic Impregnated Carbon Ablator).
    pub fn pica() -> Self {
        Self {
            virgin_density_kg_m3: 256.0,
            char_density_kg_m3: 182.0,
            pyrolysis_onset_c: 400.0,
            pyrolysis_complete_c: 850.0,
            heat_of_pyrolysis_kj_kg: 1750.0,
            char_cp_j_kgk: 1500.0,
            virgin_cp_j_kgk: 1200.0,
        }
    }

    /// Local char fraction β (0 = fully virgin, 1 = fully charred) as a
    /// function of local temperature, using a linear pyrolysis model.
    pub fn char_fraction(&self, temperature_c: f64) -> f64 {
        if temperature_c <= self.pyrolysis_onset_c {
            0.0
        } else if temperature_c >= self.pyrolysis_complete_c {
            1.0
        } else {
            (temperature_c - self.pyrolysis_onset_c)
                / (self.pyrolysis_complete_c - self.pyrolysis_onset_c)
        }
    }

    /// Local bulk density (kg/m³) as a function of char fraction.
    pub fn local_density_kg_m3(&self, beta: f64) -> f64 {
        self.virgin_density_kg_m3 * (1.0 - beta) + self.char_density_kg_m3 * beta
    }

    /// Blowing correction factor B' (dimensionless) for surface energy
    /// balance in a laminar boundary layer (simplified Sutton–Graves).
    ///
    /// B' = ṁ_w / (ρ_e · u_e · C_H) where ṁ_w is the ablation mass flux.
    /// Uses a simplified correlation: B' ≈ 0.5 · (T_wall / T_recovery)^0.5
    pub fn blowing_correction(&self, t_wall_c: f64, t_recovery_c: f64) -> f64 {
        let t_wall = t_wall_c + 273.15;
        let t_rec = (t_recovery_c + 273.15).max(1.0);
        0.5 * (t_wall / t_rec).sqrt()
    }

    /// Surface recession rate (mm/s) as a function of heat flux (MW/m²) and
    /// wall temperature (°C), using the Arrhenius char oxidation kinetics.
    ///
    /// ṡ = A_ox · exp(−E_a / RT) · q'' / (ρ_c · H_v)
    pub fn recession_rate_mm_s(&self, heat_flux_mw_m2: f64, wall_temperature_c: f64) -> f64 {
        const A_OX: f64 = 2.5e4; // pre-exponential factor
        const E_A: f64 = 1.8e5; // J/mol activation energy
        const R: f64 = 8.314;
        let t_k = wall_temperature_c + 273.15;
        let kinetic = A_OX * (-E_A / (R * t_k)).exp();
        let hv = self.heat_of_pyrolysis_kj_kg * 1000.0; // J/kg
        let rho_c = self.char_density_kg_m3;
        let q = heat_flux_mw_m2 * 1e6; // W/m²
        kinetic * q / (rho_c * hv) * 1000.0 // mm/s
    }
}

// ---------------------------------------------------------------------------
// 7. Metal foam — Gibson–Ashby model
// ---------------------------------------------------------------------------

/// Metal foam descriptor using the Gibson–Ashby cellular solid model.
#[derive(Debug, Clone, Copy)]
pub struct MetalFoam {
    /// Relative density ρ*/ρₛ (0–1).
    pub relative_density: f64,
    /// Parent solid Young's modulus (GPa).
    pub solid_modulus_gpa: f64,
    /// Parent solid yield strength (MPa).
    pub solid_yield_mpa: f64,
    /// Parent solid density (kg/m³).
    pub solid_density_kg_m3: f64,
    /// Cell topology: true = open-cell, false = closed-cell.
    pub open_cell: bool,
}

impl MetalFoam {
    /// Aluminium open-cell foam (e.g. ERG Duocel 6101-T6).
    pub fn aluminium_open_cell(relative_density: f64) -> Self {
        Self {
            relative_density,
            solid_modulus_gpa: 70.0,
            solid_yield_mpa: 270.0,
            solid_density_kg_m3: 2700.0,
            open_cell: true,
        }
    }

    /// Foam density (kg/m³).
    pub fn density_kg_m3(&self) -> f64 {
        self.relative_density * self.solid_density_kg_m3
    }

    /// Young's modulus (GPa) from Gibson–Ashby scaling.
    ///
    /// Open cell:   E* = C₁ · Eₛ · (ρ*/ρₛ)²
    /// Closed cell: E* = φ² · Eₛ · (ρ*/ρₛ)² + (1−φ) · Eₛ · (ρ*/ρₛ)
    ///                   (membrane term φ = 0.6 for typical closed-cell)
    pub fn youngs_modulus_gpa(&self) -> f64 {
        let r = self.relative_density;
        if self.open_cell {
            self.solid_modulus_gpa * r * r
        } else {
            let phi = 0.6f64;
            self.solid_modulus_gpa * (phi * phi * r * r + (1.0 - phi) * r)
        }
    }

    /// Plateau (yield) stress (MPa) — onset of plastic collapse.
    ///
    /// Open cell:   σ_y* = C₂ · σ_yₛ · (ρ*/ρₛ)^(3/2)
    /// Closed cell: σ_y* = 0.3 · σ_yₛ · ((φ · ρ*/ρₛ)^(1/2) + (1−φ)(ρ*/ρₛ))
    pub fn plateau_stress_mpa(&self) -> f64 {
        let r = self.relative_density;
        if self.open_cell {
            0.3 * self.solid_yield_mpa * r.powf(1.5)
        } else {
            let phi = 0.6f64;
            0.3 * self.solid_yield_mpa * ((phi * r).sqrt() + (1.0 - phi) * r)
        }
    }

    /// Densification strain ε_D (dimensionless) at which the foam has fully
    /// collapsed.
    pub fn densification_strain(&self) -> f64 {
        1.0 - 1.4 * self.relative_density
    }

    /// Energy absorption capacity (MJ/m³) up to densification, assuming a
    /// perfectly rectangular stress–strain plateau.
    pub fn energy_absorption_mj_m3(&self) -> f64 {
        self.plateau_stress_mpa() * self.densification_strain() * 1e-3
    }
}

// ---------------------------------------------------------------------------
// 8. Shape memory alloy (NiTi Nitinol) — superelastic model
// ---------------------------------------------------------------------------

/// Nitinol SMA material descriptor with superelastic behaviour.
///
/// Uses the simplified Brinson 1-D constitutive model.
#[derive(Debug, Clone)]
pub struct Nitinol {
    /// Service temperature (°C).
    pub temperature_c: f64,
    /// Austenite finish temperature A_f (°C).
    pub a_finish_c: f64,
    /// Martensite start temperature M_s (°C).
    pub m_start_c: f64,
    /// Austenite Young's modulus (GPa).
    pub ea_gpa: f64,
    /// Martensite Young's modulus (GPa).
    pub em_gpa: f64,
}

impl Nitinol {
    /// Biomedical / interventional-device grade Nitinol.
    pub fn biomedical_grade() -> Self {
        Self {
            temperature_c: 37.0, // body temperature
            a_finish_c: 10.0,
            m_start_c: -5.0,
            ea_gpa: 83.0,
            em_gpa: 28.0,
        }
    }

    /// Returns `true` when the alloy is in the superelastic (austenite) regime
    /// at the current temperature.
    pub fn is_superelastic(&self) -> bool {
        self.temperature_c > self.a_finish_c
    }

    /// Effective Young's modulus (GPa) using a rule-of-mixtures mixture
    /// of austenite and martensite fractions.
    ///
    /// ξ = martensite fraction (0 = fully austenite, 1 = fully martensite)
    pub fn effective_modulus_gpa(&self, martensite_fraction: f64) -> f64 {
        let xi = martensite_fraction.clamp(0.0, 1.0);
        self.ea_gpa * (1.0 - xi) + self.em_gpa * xi
    }

    /// Critical stress for onset of stress-induced martensite (MPa).
    ///
    /// σ_cr = C_A · (T − A_f) where C_A ≈ 8 MPa/°C for Nitinol.
    pub fn critical_sim_stress_mpa(&self) -> f64 {
        const C_A: f64 = 8.0; // MPa/°C
        let delta_t = (self.temperature_c - self.a_finish_c).max(0.0);
        200.0 + C_A * delta_t // 200 MPa is the lower plateau reference
    }

    /// Superelastic strain recovery (%) assuming the maximum transformation
    /// strain is 8 %.
    pub fn max_recoverable_strain_pct(&self) -> f64 {
        if self.is_superelastic() { 8.0 } else { 2.0 }
    }

    /// Density (kg/m³).
    pub fn density_kg_m3(&self) -> f64 {
        6450.0
    }
}

// ---------------------------------------------------------------------------
// 9. High-entropy alloy (HEA) mechanical model
// ---------------------------------------------------------------------------

/// Mechanical property model for a CoCrFeMnNi (Cantor) equiatomic HEA.
///
/// Reference: Gali & George, Intermetallics 2013; Otto et al. 2013.
#[derive(Debug, Clone)]
pub struct CantorhHea {
    /// Test/service temperature (°C).
    pub temperature_c: f64,
    /// Grain size (µm).
    pub grain_size_um: f64,
}

impl CantorhHea {
    /// Construct a Cantor HEA descriptor.
    pub fn new(temperature_c: f64, grain_size_um: f64) -> Self {
        Self {
            temperature_c,
            grain_size_um,
        }
    }

    /// Young's modulus (GPa).
    pub fn youngs_modulus_gpa(&self) -> f64 {
        let t = self.temperature_c;
        202.0 - 0.06 * t
    }

    /// Yield strength (MPa) including Hall–Petch grain-size strengthening.
    ///
    /// σ_y = σ₀ + k_HP / √d + thermal softening term
    pub fn yield_strength_mpa(&self) -> f64 {
        let sigma0 = 160.0; // MPa lattice friction stress
        let k_hp = 800.0; // MPa·µm^0.5 Hall-Petch coefficient
        let d = self.grain_size_um.max(1.0);
        let t = self.temperature_c;
        let hp = k_hp / d.sqrt();
        let thermal = (0.4 * t).max(0.0); // thermal softening ~0.4 MPa/°C
        (sigma0 + hp - thermal).max(50.0)
    }

    /// Ultimate tensile strength (MPa).
    pub fn uts_mpa(&self) -> f64 {
        self.yield_strength_mpa() * 1.5 // HEA exhibits high work hardening
    }

    /// Fracture toughness KIC (MPa·√m) — the Cantor alloy is notable for
    /// cryogenic toughening.
    pub fn fracture_toughness_mpa_sqrtm(&self) -> f64 {
        let t = self.temperature_c;
        if t < 0.0 {
            // Toughness increases at cryogenic temperatures
            217.0 - 0.5 * t // increases as T decreases below 0 °C
        } else {
            (217.0 - 0.1 * t).max(100.0)
        }
    }

    /// Stacking fault energy (mJ/m²) — low SFE in the Cantor alloy promotes
    /// twinning-induced plasticity.
    pub fn stacking_fault_energy_mj_m2(&self) -> f64 {
        let t = self.temperature_c;
        18.0 + 0.02 * t
    }

    /// Density (kg/m³).
    pub fn density_kg_m3(&self) -> f64 {
        8000.0
    }
}

// ---------------------------------------------------------------------------
// 10. Reentry vehicle thermal protection
// ---------------------------------------------------------------------------

/// Aerodynamic heating model for a blunt-nosed reentry vehicle.
///
/// Uses the Fay–Riddell stagnation-point heat flux correlation and a
/// one-dimensional heat conduction model for TPS sizing.
#[derive(Debug, Clone)]
pub struct ReentryVehicle {
    /// Nose radius (m).
    pub nose_radius_m: f64,
    /// Entry velocity (m/s).
    pub entry_velocity_m_s: f64,
    /// Freestream density (kg/m³).
    pub freestream_density_kg_m3: f64,
    /// TPS material type label (informational).
    pub tps_material: String,
    /// TPS thickness (mm).
    pub tps_thickness_mm: f64,
}

impl ReentryVehicle {
    /// Construct a reentry vehicle TPS descriptor.
    pub fn new(
        nose_radius_m: f64,
        entry_velocity_m_s: f64,
        freestream_density_kg_m3: f64,
        tps_material: &str,
        tps_thickness_mm: f64,
    ) -> Self {
        Self {
            nose_radius_m,
            entry_velocity_m_s,
            freestream_density_kg_m3,
            tps_material: tps_material.to_string(),
            tps_thickness_mm,
        }
    }

    /// Stagnation-point cold-wall heat flux (W/m²) using the Sutton–Graves
    /// correlation:
    ///
    /// q'' = k_sg · √(ρ/R_N) · V³
    ///
    /// where k_sg = 1.742×10⁻⁴ (SI units; ρ in kg/m³, R_N in m, V in m/s).
    pub fn stagnation_heat_flux_w_m2(&self) -> f64 {
        const K_SG: f64 = 1.742e-4;
        let rho = self.freestream_density_kg_m3;
        let r = self.nose_radius_m;
        let v = self.entry_velocity_m_s;
        K_SG * (rho / r).sqrt() * v.powi(3)
    }

    /// Peak surface temperature at the stagnation point (K) estimated from
    /// a radiation equilibrium condition:
    ///
    /// q'' = ε · σ · T⁴  →  T = (q'' / (ε·σ))^(1/4)
    pub fn radiative_equilibrium_temperature_k(&self, emissivity: f64) -> f64 {
        const SIGMA: f64 = 5.670_374_419e-8; // Stefan–Boltzmann constant, W/m²·K⁴
        let q = self.stagnation_heat_flux_w_m2();
        let eps = emissivity.clamp(0.1, 1.0);
        (q / (eps * SIGMA)).powf(0.25)
    }

    /// Integrated heat load (MJ/m²) for a given entry duration (s).
    ///
    /// Uses the parabolic profile q(t) = q_peak · 4t(t_entry-t)/t_entry².
    pub fn integrated_heat_load_mj_m2(&self, entry_duration_s: f64) -> f64 {
        let q_peak = self.stagnation_heat_flux_w_m2();
        // Integral of 4·q·t(T-t)/T² dt from 0 to T = (2/3)·q·T
        2.0 / 3.0 * q_peak * entry_duration_s / 1e6
    }

    /// Required TPS thickness (mm) to limit bondline temperature below
    /// `max_bondline_c` (°C), using a 1-D semi-infinite slab approximation.
    ///
    /// δ ≈ 2 · √(α · t_entry) · erfinv(1 − ΔT_allowed / ΔT_surface)
    /// (simplified as δ ≈ 1.5 · √(α · t_entry) for ΔT_allowed = 50 °C).
    pub fn required_tps_thickness_mm(
        &self,
        tps_thermal_diffusivity_m2_s: f64,
        entry_duration_s: f64,
    ) -> f64 {
        let alpha = tps_thermal_diffusivity_m2_s;
        let t = entry_duration_s;
        1.5 * (alpha * t).sqrt() * 1000.0 // convert m to mm
    }

    /// Vehicle ballistic coefficient (kg/m²): β = m / (C_D · A_ref).
    ///
    /// Higher β → deeper penetration before deceleration.
    pub fn ballistic_coefficient_kg_m2(&self, mass_kg: f64, cd: f64) -> f64 {
        let a_ref = PI * self.nose_radius_m * self.nose_radius_m;
        mass_kg / (cd * a_ref)
    }
}

// ---------------------------------------------------------------------------
// 11. Additional helpers — Titanium structural member sizing
// ---------------------------------------------------------------------------

/// Titanium structural beam sizing for aerospace applications.
///
/// Sizes a hollow circular tube to carry a given axial load with a safety
/// factor.
#[derive(Debug, Clone, Copy)]
pub struct TitaniumTube {
    /// Outer diameter (mm).
    pub outer_diameter_mm: f64,
    /// Wall thickness (mm).
    pub wall_thickness_mm: f64,
    /// Temperature (°C).
    pub temperature_c: f64,
}

impl TitaniumTube {
    /// Construct a titanium tube.
    pub fn new(outer_diameter_mm: f64, wall_thickness_mm: f64, temperature_c: f64) -> Self {
        Self {
            outer_diameter_mm,
            wall_thickness_mm,
            temperature_c,
        }
    }

    /// Cross-sectional area (mm²).
    pub fn cross_section_area_mm2(&self) -> f64 {
        let ro = self.outer_diameter_mm * 0.5;
        let ri = ro - self.wall_thickness_mm;
        PI * (ro * ro - ri * ri)
    }

    /// Second moment of area I (mm⁴).
    pub fn second_moment_mm4(&self) -> f64 {
        let ro = self.outer_diameter_mm * 0.5;
        let ri = ro - self.wall_thickness_mm;
        PI / 4.0 * (ro.powi(4) - ri.powi(4))
    }

    /// Euler buckling load (N) for a pin-ended column of length `l` (mm).
    pub fn euler_buckling_load_n(&self, length_mm: f64) -> f64 {
        let mat = Ti6Al4V::new(self.temperature_c);
        let e = mat.youngs_modulus_gpa() * 1e3; // MPa
        PI * PI * e * self.second_moment_mm4() / (length_mm * length_mm)
    }

    /// Axial load capacity (N) limited by yield.
    pub fn axial_yield_load_n(&self) -> f64 {
        let mat = Ti6Al4V::new(self.temperature_c);
        mat.yield_strength_mpa() * self.cross_section_area_mm2()
    }

    /// Margin of safety against yield under `applied_load_n` (positive = safe).
    pub fn margin_of_safety_yield(&self, applied_load_n: f64) -> f64 {
        self.axial_yield_load_n() / applied_load_n - 1.0
    }
}

// ---------------------------------------------------------------------------
// 12. Oxidation kinetics of superalloy bond coat
// ---------------------------------------------------------------------------

/// TGO (thermally-grown oxide) growth kinetics for MCrAlY bond coat.
///
/// Uses a parabolic oxidation rate law: h² = k_p · t, where h is oxide
/// thickness and k_p is the parabolic rate constant (µm²/h).
#[derive(Debug, Clone)]
pub struct BondCoatOxidation {
    /// Temperature (°C).
    pub temperature_c: f64,
    /// Al content of the bond coat (wt%).
    pub al_content_wt_pct: f64,
}

impl BondCoatOxidation {
    /// Construct an oxidation model.
    pub fn new(temperature_c: f64, al_content_wt_pct: f64) -> Self {
        Self {
            temperature_c,
            al_content_wt_pct,
        }
    }

    /// Parabolic rate constant kₚ (µm²/h) using an Arrhenius correlation.
    ///
    /// kₚ = A_0 · (Al/10)^0.5 · exp(−Q / RT)
    /// with A_0 = 2×10⁶ µm²/h, Q = 220 kJ/mol.
    pub fn parabolic_rate_constant_um2_h(&self) -> f64 {
        const A0: f64 = 2.0e6; // µm²/h
        const Q: f64 = 220_000.0; // J/mol
        const R: f64 = 8.314;
        let t_k = self.temperature_c + 273.15;
        let al_factor = (self.al_content_wt_pct / 10.0).sqrt().max(0.1);
        A0 * al_factor * (-Q / (R * t_k)).exp()
    }

    /// TGO thickness (µm) after `time_h` hours of oxidation.
    pub fn tgo_thickness_um(&self, time_h: f64) -> f64 {
        (self.parabolic_rate_constant_um2_h() * time_h).sqrt()
    }

    /// Critical TGO thickness (µm) above which spallation is likely (~7 µm
    /// for MCrAlY/7YSZ systems).
    pub fn critical_tgo_thickness_um() -> f64 {
        7.0
    }

    /// Life fraction consumed (0–1) at a given TGO thickness.
    pub fn life_fraction(&self, tgo_thickness_um: f64) -> f64 {
        (tgo_thickness_um / Self::critical_tgo_thickness_um())
            .powi(2)
            .min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Ti-6Al-4V ────────────────────────────────────────────────────────────

    #[test]
    fn test_ti64_youngs_modulus_at_rt() {
        let mat = Ti6Al4V::new(20.0);
        let e = mat.youngs_modulus_gpa();
        // At room temperature ~112–114 GPa
        assert!(e > 110.0 && e < 115.0, "E = {e}");
    }

    #[test]
    fn test_ti64_modulus_decreases_with_temperature() {
        let e_rt = Ti6Al4V::new(20.0).youngs_modulus_gpa();
        let e_hot = Ti6Al4V::new(600.0).youngs_modulus_gpa();
        assert!(e_hot < e_rt);
    }

    #[test]
    fn test_ti64_yield_strength_decreases() {
        let sy_rt = Ti6Al4V::new(20.0).yield_strength_mpa();
        let sy_hot = Ti6Al4V::new(500.0).yield_strength_mpa();
        assert!(sy_hot < sy_rt);
    }

    #[test]
    fn test_ti64_thermal_conductivity_increases() {
        let k_rt = Ti6Al4V::new(20.0).thermal_conductivity_w_mk();
        let k_hot = Ti6Al4V::new(500.0).thermal_conductivity_w_mk();
        assert!(k_hot > k_rt);
    }

    #[test]
    fn test_ti64_thermal_diffusivity_positive() {
        let alpha = Ti6Al4V::new(300.0).thermal_diffusivity_m2_s();
        assert!(alpha > 0.0);
    }

    // ── Inconel 718 ──────────────────────────────────────────────────────────

    #[test]
    fn test_inconel_creep_rate_increases_with_stress() {
        let mat = Inconel718::new(650.0);
        let cr_low = mat.norton_creep_rate_per_s(300.0);
        let cr_high = mat.norton_creep_rate_per_s(600.0);
        assert!(cr_high > cr_low);
    }

    #[test]
    fn test_inconel_creep_rate_increases_with_temperature() {
        let cr_cool = Inconel718::new(540.0).norton_creep_rate_per_s(500.0);
        let cr_hot = Inconel718::new(760.0).norton_creep_rate_per_s(500.0);
        assert!(cr_hot > cr_cool);
    }

    #[test]
    fn test_inconel_rupture_life_decreases_with_stress() {
        let mat = Inconel718::new(650.0);
        let life_low = mat.rupture_life_hours(400.0);
        let life_high = mat.rupture_life_hours(700.0);
        assert!(life_high < life_low);
    }

    #[test]
    fn test_inconel_density() {
        let mat = Inconel718::new(25.0);
        assert!((mat.density_kg_m3() - 8190.0).abs() < 1.0);
    }

    // ── CFRP / CLT ───────────────────────────────────────────────────────────

    #[test]
    fn test_cfrp_ply_reduced_stiffness_positive() {
        let ply = CfrpPly::im7_8552(0.0);
        let q = ply.reduced_stiffness_gpa();
        for &qi in &q {
            assert!(qi > 0.0, "Q component non-positive: {qi}");
        }
    }

    #[test]
    fn test_clt_quasi_isotropic_laminate() {
        // [0/45/-45/90]s quasi-isotropic laminate
        let angles = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0];
        let plies: Vec<CfrpPly> = angles.iter().map(|&a| CfrpPly::im7_8552(a)).collect();
        let res = clt_analysis(&plies);
        // Quasi-isotropic → Ex ≈ Ey within 5 %
        let diff = (res.ex_gpa - res.ey_gpa).abs() / res.ex_gpa;
        assert!(
            diff < 0.05,
            "Ex={:.2} Ey={:.2} diff={diff:.4}",
            res.ex_gpa,
            res.ey_gpa
        );
    }

    #[test]
    fn test_clt_zero_only_laminate_max_modulus() {
        let plies = vec![CfrpPly::im7_8552(0.0); 8];
        let unidirectional = clt_analysis(&plies);
        // Ex should be close to E1 = 161 GPa for all-0° laminate
        assert!(unidirectional.ex_gpa > 100.0);
    }

    #[test]
    fn test_clt_thickness_sum() {
        let plies: Vec<CfrpPly> = [0.0, 90.0, 0.0, 90.0]
            .iter()
            .map(|&a| CfrpPly::im7_8552(a))
            .collect();
        let res = clt_analysis(&plies);
        let expected = 4.0 * 0.125;
        assert!((res.total_thickness_mm - expected).abs() < 1e-9);
    }

    // ── SiC/SiC CMC ─────────────────────────────────────────────────────────

    #[test]
    fn test_sic_sic_modulus_range() {
        let cmc = SicSicCmc::new(1000.0, 0.45);
        let e = cmc.youngs_modulus_gpa();
        assert!(e > 100.0 && e < 220.0, "E = {e}");
    }

    #[test]
    fn test_sic_sic_uts_decreases_with_temperature() {
        let uts_low = SicSicCmc::new(800.0, 0.45).uts_mpa();
        let uts_high = SicSicCmc::new(1200.0, 0.45).uts_mpa();
        assert!(uts_high < uts_low);
    }

    #[test]
    fn test_sic_sic_max_use_temperature() {
        let cmc = SicSicCmc::new(25.0, 0.40);
        assert!((cmc.max_use_temperature_c() - 1200.0).abs() < 1.0);
    }

    // ── TBC ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_tbc_temperature_drop_proportional_to_flux() {
        let tbc = TbcYsz::new(120.0, 1300.0, 950.0);
        let dt1 = tbc.temperature_drop_c(1e6);
        let dt2 = tbc.temperature_drop_c(2e6);
        assert!((dt2 - 2.0 * dt1).abs() < 1e-6, "dt1={dt1} dt2={dt2}");
    }

    #[test]
    fn test_tbc_spallation_life_decreases_with_delta_t() {
        let tbc = TbcYsz::new(120.0, 1300.0, 950.0);
        let n_low = tbc.spallation_life_cycles(100.0);
        let n_high = tbc.spallation_life_cycles(200.0);
        assert!(n_high < n_low);
    }

    #[test]
    fn test_tbc_cte_mismatch_strain_sign() {
        let tbc = TbcYsz::new(120.0, 1300.0, 950.0);
        // Substrate has higher CTE → mismatch strain is positive (tensile in coating)
        let strain = tbc.cte_mismatch_strain(500.0);
        assert!(strain > 0.0);
    }

    // ── Ablative material ────────────────────────────────────────────────────

    #[test]
    fn test_ablative_char_fraction_below_onset() {
        let mat = AblativeMaterial::pica();
        assert!((mat.char_fraction(300.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_ablative_char_fraction_above_complete() {
        let mat = AblativeMaterial::pica();
        assert!((mat.char_fraction(1000.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ablative_char_fraction_midpoint() {
        let mat = AblativeMaterial::pica();
        let mid = (mat.pyrolysis_onset_c + mat.pyrolysis_complete_c) * 0.5;
        let beta = mat.char_fraction(mid);
        assert!((beta - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_ablative_recession_rate_increases_with_flux() {
        let mat = AblativeMaterial::pica();
        let r1 = mat.recession_rate_mm_s(1.0, 800.0);
        let r2 = mat.recession_rate_mm_s(5.0, 800.0);
        assert!(r2 > r1);
    }

    // ── Metal foam ───────────────────────────────────────────────────────────

    #[test]
    fn test_foam_modulus_increases_with_density() {
        let f1 = MetalFoam::aluminium_open_cell(0.05);
        let f2 = MetalFoam::aluminium_open_cell(0.10);
        assert!(f2.youngs_modulus_gpa() > f1.youngs_modulus_gpa());
    }

    #[test]
    fn test_foam_plateau_stress_positive() {
        let foam = MetalFoam::aluminium_open_cell(0.08);
        assert!(foam.plateau_stress_mpa() > 0.0);
    }

    #[test]
    fn test_foam_densification_strain_range() {
        let foam = MetalFoam::aluminium_open_cell(0.10);
        let eps_d = foam.densification_strain();
        assert!(eps_d > 0.0 && eps_d < 1.0, "eps_d = {eps_d}");
    }

    #[test]
    fn test_foam_energy_absorption_positive() {
        let foam = MetalFoam::aluminium_open_cell(0.10);
        assert!(foam.energy_absorption_mj_m3() > 0.0);
    }

    // ── Nitinol ──────────────────────────────────────────────────────────────

    #[test]
    fn test_nitinol_superelastic_at_body_temp() {
        let niti = Nitinol::biomedical_grade();
        assert!(niti.is_superelastic());
    }

    #[test]
    fn test_nitinol_effective_modulus_bounds() {
        let niti = Nitinol::biomedical_grade();
        let e_aust = niti.effective_modulus_gpa(0.0);
        let e_mart = niti.effective_modulus_gpa(1.0);
        assert!((e_aust - 83.0).abs() < 1e-9);
        assert!((e_mart - 28.0).abs() < 1e-9);
    }

    #[test]
    fn test_nitinol_critical_stress_increases_with_temperature() {
        let niti_hot = Nitinol {
            temperature_c: 60.0,
            ..Nitinol::biomedical_grade()
        };
        let niti_rt = Nitinol::biomedical_grade();
        assert!(niti_hot.critical_sim_stress_mpa() > niti_rt.critical_sim_stress_mpa());
    }

    // ── Cantor HEA ───────────────────────────────────────────────────────────

    #[test]
    fn test_hea_yield_hall_petch_effect() {
        let fine = CantorhHea::new(25.0, 10.0); // fine grain
        let coarse = CantorhHea::new(25.0, 100.0); // coarse grain
        assert!(fine.yield_strength_mpa() > coarse.yield_strength_mpa());
    }

    #[test]
    fn test_hea_cryogenic_toughening() {
        let cryo = CantorhHea::new(-196.0, 50.0);
        let rt = CantorhHea::new(25.0, 50.0);
        assert!(cryo.fracture_toughness_mpa_sqrtm() > rt.fracture_toughness_mpa_sqrtm());
    }

    #[test]
    fn test_hea_uts_greater_than_yield() {
        let hea = CantorhHea::new(25.0, 50.0);
        assert!(hea.uts_mpa() > hea.yield_strength_mpa());
    }

    // ── Reentry vehicle ──────────────────────────────────────────────────────

    #[test]
    fn test_reentry_heat_flux_positive() {
        let rv = ReentryVehicle::new(0.5, 7500.0, 0.001, "PICA", 80.0);
        assert!(rv.stagnation_heat_flux_w_m2() > 0.0);
    }

    #[test]
    fn test_reentry_equilibrium_temperature_increases_with_flux() {
        let rv_fast = ReentryVehicle::new(0.5, 8000.0, 0.001, "PICA", 80.0);
        let rv_slow = ReentryVehicle::new(0.5, 5000.0, 0.001, "PICA", 80.0);
        let t_fast = rv_fast.radiative_equilibrium_temperature_k(0.85);
        let t_slow = rv_slow.radiative_equilibrium_temperature_k(0.85);
        assert!(t_fast > t_slow);
    }

    #[test]
    fn test_reentry_integrated_heat_load_positive() {
        let rv = ReentryVehicle::new(0.5, 7500.0, 0.001, "PICA", 80.0);
        assert!(rv.integrated_heat_load_mj_m2(60.0) > 0.0);
    }

    #[test]
    fn test_required_tps_thickness_increases_with_time() {
        let rv = ReentryVehicle::new(0.5, 7500.0, 0.001, "PICA", 80.0);
        let alpha = 3.0e-7; // m²/s for PICA char
        let t1 = rv.required_tps_thickness_mm(alpha, 60.0);
        let t2 = rv.required_tps_thickness_mm(alpha, 120.0);
        assert!(t2 > t1);
    }

    // ── Titanium tube ────────────────────────────────────────────────────────

    #[test]
    fn test_ti_tube_cross_section_positive() {
        let tube = TitaniumTube::new(25.0, 1.5, 20.0);
        assert!(tube.cross_section_area_mm2() > 0.0);
    }

    #[test]
    fn test_ti_tube_euler_buckling_longer_is_weaker() {
        let tube = TitaniumTube::new(25.0, 1.5, 20.0);
        let p_short = tube.euler_buckling_load_n(200.0);
        let p_long = tube.euler_buckling_load_n(500.0);
        assert!(p_short > p_long);
    }

    #[test]
    fn test_ti_tube_margin_of_safety_positive_under_low_load() {
        let tube = TitaniumTube::new(25.0, 1.5, 20.0);
        let mos = tube.margin_of_safety_yield(100.0); // 100 N — very light load
        assert!(mos > 0.0, "MoS = {mos}");
    }

    // ── Bond coat oxidation ──────────────────────────────────────────────────

    #[test]
    fn test_tgo_growth_parabolic() {
        let ox = BondCoatOxidation::new(1050.0, 8.0);
        let h1 = ox.tgo_thickness_um(100.0);
        let h4 = ox.tgo_thickness_um(400.0);
        // Parabolic: h(4t) = √4 · h(t) = 2·h(t)
        assert!((h4 - 2.0 * h1).abs() < 1e-6, "h1={h1} h4={h4}");
    }

    #[test]
    fn test_tgo_life_fraction_at_critical() {
        let ox = BondCoatOxidation::new(1050.0, 8.0);
        let h_crit = BondCoatOxidation::critical_tgo_thickness_um();
        let lf = ox.life_fraction(h_crit);
        assert!((lf - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_tgo_rate_constant_increases_with_temperature() {
        let kp_low = BondCoatOxidation::new(950.0, 8.0).parabolic_rate_constant_um2_h();
        let kp_high = BondCoatOxidation::new(1100.0, 8.0).parabolic_rate_constant_um2_h();
        assert!(kp_high > kp_low);
    }
}
