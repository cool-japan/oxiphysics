// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geological and geomechanical material models: soils, rocks, sand,
//! permafrost, seabed sediments, and granular pressure models.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Soil Model — Mohr-Coulomb
// ─────────────────────────────────────────────────────────────────────────────

/// Soil model using the Mohr-Coulomb yield criterion.
///
/// Models elastic-perfectly-plastic soil behavior with friction and cohesion.
#[derive(Debug, Clone)]
pub struct SoilModel {
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
    /// Internal friction angle φ \[radians\].
    pub friction_angle: f64,
    /// Dilatancy angle ψ \[radians\].
    pub dilatancy_angle: f64,
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
}

impl SoilModel {
    /// Create a new Mohr-Coulomb soil model.
    pub fn new(
        cohesion: f64,
        friction_angle_deg: f64,
        dilatancy_angle_deg: f64,
        young_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            cohesion,
            friction_angle: friction_angle_deg.to_radians(),
            dilatancy_angle: dilatancy_angle_deg.to_radians(),
            young_modulus,
            poisson_ratio,
        }
    }

    /// Typical medium-density sand.
    pub fn medium_sand() -> Self {
        Self::new(0.0, 32.0, 5.0, 30e6, 0.3)
    }

    /// Typical stiff clay.
    pub fn stiff_clay() -> Self {
        Self::new(60e3, 22.0, 0.0, 80e6, 0.35)
    }

    /// Mohr-Coulomb yield function F(σ1, σ3).
    ///
    /// F = (σ1 - σ3) - 2*c*cos(φ) - (σ1 + σ3)*sin(φ)
    /// F < 0: elastic; F = 0: on yield surface.
    pub fn yield_function(&self, sigma1: f64, sigma3: f64) -> f64 {
        let phi = self.friction_angle;
        (sigma1 - sigma3) - 2.0 * self.cohesion * phi.cos() - (sigma1 + sigma3) * phi.sin()
    }

    /// Check if the stress state is yielding.
    pub fn is_yielding(&self, sigma1: f64, sigma3: f64) -> bool {
        self.yield_function(sigma1, sigma3) >= 0.0
    }

    /// Plastic multiplier Δλ for a given incremental plastic strain magnitude.
    ///
    /// Δλ = Δεp / (1 + sin(ψ))
    pub fn plastic_multiplier(&self, d_eps_p: f64) -> f64 {
        d_eps_p / (1.0 + self.dilatancy_angle.sin())
    }

    /// Critical state line (CSL) deviator stress q as a function of mean stress p \[Pa\].
    ///
    /// q = M_cs * p  where M_cs = 6*sin(φ) / (3 - sin(φ)) for compression.
    pub fn critical_state_line_q(&self, p: f64) -> f64 {
        let phi = self.friction_angle;
        let m_cs = 6.0 * phi.sin() / (3.0 - phi.sin());
        m_cs * p
    }

    /// Bulk modulus K = E / (3*(1-2ν)) \[Pa\].
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus G = E / (2*(1+ν)) \[Pa\].
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Uniaxial compressive strength (unconfined) \[Pa\].
    ///
    /// qu = 2*c*cos(φ) / (1 - sin(φ))
    pub fn ucs(&self) -> f64 {
        let phi = self.friction_angle;
        2.0 * self.cohesion * phi.cos() / (1.0 - phi.sin())
    }

    /// Ko (coefficient of earth pressure at rest) via Jaky's formula.
    ///
    /// Ko = 1 - sin(φ)
    pub fn k0_jaky(&self) -> f64 {
        1.0 - self.friction_angle.sin()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Modified Cam-Clay Model
// ─────────────────────────────────────────────────────────────────────────────

/// Modified Cam-Clay elastoplastic model for clays.
///
/// Captures volumetric hardening and critical state behavior.
#[derive(Debug, Clone)]
pub struct CamClayModel {
    /// Slope of CSL in q-p' space.
    pub m: f64,
    /// Compression index (slope of NCL in e-ln p' space).
    pub lambda: f64,
    /// Swelling index (slope of URL in e-ln p' space).
    pub kappa: f64,
    /// Specific volume at p'=1 kPa on CSL: N.
    pub n_value: f64,
    /// Initial void ratio.
    pub e0: f64,
    /// Initial preconsolidation pressure p'c0 \[Pa\].
    pub p_c0: f64,
}

impl CamClayModel {
    /// Create a new Modified Cam-Clay model.
    pub fn new(m: f64, lambda: f64, kappa: f64, n_value: f64, e0: f64, p_c0: f64) -> Self {
        Self {
            m,
            lambda,
            kappa,
            n_value,
            e0,
            p_c0,
        }
    }

    /// Typical normally consolidated kaolin clay.
    pub fn kaolin_ncc() -> Self {
        Self::new(0.92, 0.26, 0.05, 3.0, 1.5, 100e3)
    }

    /// Void ratio on CSL at mean stress p' \[Pa\].
    ///
    /// e_cs = N - λ * ln(p'/1 kPa)
    pub fn e_cs(&self, p: f64) -> f64 {
        self.n_value - self.lambda * (p / 1e3).ln()
    }

    /// Modified Cam-Clay yield surface value F(p, q, p0).
    ///
    /// F = q^2 / M^2 + p*(p - p0)  ; F < 0: elastic
    pub fn yield_surface_p0(&self, p: f64, q: f64) -> f64 {
        q * q / (self.m * self.m) + p * (p - self.p_c0)
    }

    /// Check if on or outside yield surface.
    pub fn is_yielding(&self, p: f64, q: f64) -> bool {
        self.yield_surface_p0(p, q) >= 0.0
    }

    /// Preconsolidation pressure from current void ratio.
    ///
    /// p0 = exp((N - e) / λ) \[Pa\] (×1 kPa reference)
    pub fn preconsolidation_pressure(&self, void_ratio: f64) -> f64 {
        let ln_p = (self.n_value - void_ratio) / self.lambda;
        1e3 * ln_p.exp()
    }

    /// Elastic volumetric strain increment for a change in mean stress dp at current p.
    ///
    /// dεv_e = (κ / (1 + e0)) * dp / p
    pub fn volumetric_strain_elastic(&self, dp: f64, p: f64) -> f64 {
        if p <= 0.0 {
            return 0.0;
        }
        (self.kappa / (1.0 + self.e0)) * dp / p
    }

    /// Plastic volumetric strain increment (hardening law).
    ///
    /// dεv_p = ((λ - κ) / (1 + e0)) * dp0 / p0
    pub fn volumetric_strain_plastic(&self, dp0: f64, p0: f64) -> f64 {
        if p0 <= 0.0 {
            return 0.0;
        }
        ((self.lambda - self.kappa) / (1.0 + self.e0)) * dp0 / p0
    }

    /// Critical state mean stress at current void ratio.
    pub fn critical_state_p(&self, void_ratio: f64) -> f64 {
        let exp_arg = (self.n_value - void_ratio) / self.lambda;
        1e3 * exp_arg.exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rock Model — Hoek-Brown Criterion
// ─────────────────────────────────────────────────────────────────────────────

/// Rock mass model using the generalized Hoek-Brown failure criterion.
#[derive(Debug, Clone)]
pub struct RockModel {
    /// Uniaxial compressive strength of intact rock UCS \[Pa\].
    pub ucs: f64,
    /// Intact rock parameter mi \[-\].
    pub mi: f64,
    /// Geological strength index GSI \[-\], 0–100.
    pub gsi: f64,
    /// Disturbance factor D \[-\], 0 (undisturbed) to 1 (heavily blasted).
    pub d: f64,
}

impl RockModel {
    /// Create a new Hoek-Brown rock mass model.
    pub fn new(ucs: f64, mi: f64, gsi: f64, d: f64) -> Self {
        Self { ucs, mi, gsi, d }
    }

    /// Typical granite (hard, good rock mass).
    pub fn granite() -> Self {
        Self::new(200e6, 33.0, 75.0, 0.0)
    }

    /// Typical shale (weak, moderate rock mass).
    pub fn shale() -> Self {
        Self::new(20e6, 8.0, 40.0, 0.2)
    }

    /// Reduced m_b parameter for rock mass.
    ///
    /// mb = mi * exp((GSI - 100) / (28 - 14D))
    pub fn mb_parameter(&self) -> f64 {
        self.mi * ((self.gsi - 100.0) / (28.0 - 14.0 * self.d)).exp()
    }

    /// Rock mass s parameter (intercept of HB criterion).
    ///
    /// s = exp((GSI - 100) / (9 - 3D))
    pub fn s_parameter(&self) -> f64 {
        ((self.gsi - 100.0) / (9.0 - 3.0 * self.d)).exp()
    }

    /// Rock mass a parameter (shape exponent).
    ///
    /// a = 0.5 + (1/6)*(exp(-GSI/15) - exp(-20/3))
    pub fn a_parameter(&self) -> f64 {
        0.5 + (1.0 / 6.0) * ((-self.gsi / 15.0).exp() - (-20.0_f64 / 3.0).exp())
    }

    /// Hoek-Brown yield check: returns true if (σ1, σ3) exceeds HB envelope.
    ///
    /// σ1 ≥ σ3 + UCS*(mb*σ3/UCS + s)^a  → yielding
    pub fn yield_hoek_brown(&self, sigma1: f64, sigma3: f64) -> bool {
        let mb = self.mb_parameter();
        let s = self.s_parameter();
        let a = self.a_parameter();
        let rhs = sigma3 + self.ucs * (mb * sigma3 / self.ucs + s).max(0.0).powf(a);
        sigma1 >= rhs
    }

    /// Uniaxial compressive strength of the rock mass \[Pa\].
    ///
    /// σ_cm = UCS * s^a
    pub fn uniaxial_compressive_strength_rock(&self) -> f64 {
        let s = self.s_parameter();
        let a = self.a_parameter();
        self.ucs * s.powf(a)
    }

    /// Tensile strength of rock mass \[Pa\].
    ///
    /// σ_t = -s * UCS / mb
    pub fn tensile_strength_rock(&self) -> f64 {
        let mb = self.mb_parameter();
        let s = self.s_parameter();
        -(s * self.ucs) / mb
    }

    /// Rock mass deformation modulus \[Pa\] (Hoek et al. 2002).
    ///
    /// Em = (1 - D/2)*sqrt(UCS/100)*10^((GSI-10)/40) \[GPa\] → Pa
    pub fn deformation_modulus(&self) -> f64 {
        let factor = 1.0 - self.d / 2.0;
        let ucs_mpa = self.ucs / 1e6;
        let em_gpa = factor * (ucs_mpa / 100.0).sqrt() * 10.0_f64.powf((self.gsi - 10.0) / 40.0);
        em_gpa * 1e9
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sand Model
// ─────────────────────────────────────────────────────────────────────────────

/// Sand model with relative density and dilatancy (Bolton's approach).
#[derive(Debug, Clone)]
pub struct SandModel {
    /// Relative density Dr \[-\], 0 (loose) to 1 (dense).
    pub relative_density: f64,
    /// Critical state friction angle φ_cs \[radians\].
    pub critical_state_friction_angle: f64,
    /// Mean particle diameter d50 \[m\].
    pub d50: f64,
    /// Specific gravity Gs \[-\].
    pub specific_gravity: f64,
}

impl SandModel {
    /// Create a new sand model.
    pub fn new(
        relative_density: f64,
        critical_state_friction_angle_deg: f64,
        d50: f64,
        specific_gravity: f64,
    ) -> Self {
        Self {
            relative_density,
            critical_state_friction_angle: critical_state_friction_angle_deg.to_radians(),
            d50,
            specific_gravity,
        }
    }

    /// Typical medium dense quartz sand.
    pub fn medium_dense_quartz() -> Self {
        Self::new(0.6, 32.0, 0.3e-3, 2.65)
    }

    /// Loose sand.
    pub fn loose_sand() -> Self {
        Self::new(0.3, 30.0, 0.2e-3, 2.65)
    }

    /// Peak friction angle \[radians\] using Bolton's simplified formula.
    ///
    /// φ_p = φ_cs + 5 * Ir  (radians, Ir = dilatancy index)
    pub fn peak_friction_angle(&self, dr: f64) -> f64 {
        let ir = self.dilatancy_index(dr, 100.0);
        self.critical_state_friction_angle + 5.0_f64.to_radians() * ir
    }

    /// Bolton dilatancy index Ir = Dr*(Q - ln(p/kPa)) - R.
    ///
    /// Q = 10 (quartz), R = 1.
    pub fn dilatancy_index(&self, dr: f64, _stress_level: f64) -> f64 {
        let q = 10.0;
        let r = 1.0;
        let p_kpa = 100.0_f64; // reference stress
        (dr * (q - (p_kpa).ln()) - r).max(0.0)
    }

    /// Bolton's dilatancy formula with explicit pressure.
    ///
    /// Ir = Dr*(Q - ln(p_kPa)) - R, clamped to \[0, ∞)
    pub fn bolton_dilatancy(&self, dr: f64, p_kpa: f64) -> f64 {
        let q = 10.0;
        let r = 1.0;
        if p_kpa <= 0.0 {
            return 0.0;
        }
        (dr * (q - p_kpa.ln()) - r).max(0.0)
    }

    /// Maximum void ratio (loosest packing) — approximate.
    pub fn e_max(&self) -> f64 {
        0.9 // typical for medium sand
    }

    /// Minimum void ratio (densest packing) — approximate.
    pub fn e_min(&self) -> f64 {
        0.6 // typical for medium sand
    }

    /// Current void ratio from relative density.
    pub fn current_void_ratio(&self) -> f64 {
        self.e_max() - self.relative_density * (self.e_max() - self.e_min())
    }

    /// Dry unit weight \[kN/m^3\].
    pub fn dry_unit_weight(&self) -> f64 {
        let e = self.current_void_ratio();
        self.specific_gravity * 9.81 / (1.0 + e)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Permafrost Model
// ─────────────────────────────────────────────────────────────────────────────

/// Permafrost material model with temperature-dependent properties.
///
/// Models unfrozen water content, ice saturation, and thermal/mechanical properties.
#[derive(Debug, Clone)]
pub struct PermafrostModel {
    /// Soil type identifier.
    pub soil_type: PermafrostSoilType,
    /// Total water content (mass fraction).
    pub total_water_content: f64,
    /// Soil dry density \[kg/m^3\].
    pub dry_density: f64,
    /// Mineral thermal conductivity λ_s \[W/(m·K)\].
    pub lambda_solid: f64,
}

/// Permafrost soil type classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PermafrostSoilType {
    /// Silty clay permafrost.
    SiltyClay,
    /// Sandy permafrost.
    Sand,
    /// Gravelly permafrost.
    Gravel,
}

impl PermafrostModel {
    /// Create a new permafrost model.
    pub fn new(
        soil_type: PermafrostSoilType,
        total_water_content: f64,
        dry_density: f64,
        lambda_solid: f64,
    ) -> Self {
        Self {
            soil_type,
            total_water_content,
            dry_density,
            lambda_solid,
        }
    }

    /// Typical Siberian silty permafrost.
    pub fn siberian_silt() -> Self {
        Self::new(PermafrostSoilType::SiltyClay, 0.35, 1400.0, 2.5)
    }

    /// Unfrozen water content θ_u(T) \[m^3/m^3\] using power law.
    ///
    /// θ_u = α * |T|^(-β) for T < 0°C; total water content for T ≥ 0°C.
    pub fn unfrozen_water_content(&self, t_celsius: f64) -> f64 {
        if t_celsius >= 0.0 {
            return self.total_water_content;
        }
        let (alpha, beta) = match self.soil_type {
            PermafrostSoilType::SiltyClay => (0.45, 0.6),
            PermafrostSoilType::Sand => (0.08, 0.5),
            PermafrostSoilType::Gravel => (0.05, 0.4),
        };
        let t_abs = (-t_celsius).max(1e-6);
        (alpha * t_abs.powf(-beta)).min(self.total_water_content)
    }

    /// Ice saturation S_i at temperature T \[°C\].
    pub fn ice_saturation(&self, t_celsius: f64) -> f64 {
        if t_celsius >= 0.0 {
            return 0.0;
        }
        let theta_u = self.unfrozen_water_content(t_celsius);
        let theta_i = (self.total_water_content - theta_u).max(0.0);
        if self.total_water_content > 0.0 {
            theta_i / self.total_water_content
        } else {
            0.0
        }
    }

    /// Bulk modulus of frozen soil \[Pa\] — increases with decreasing temperature.
    ///
    /// K_frozen ≈ K_unfrozen * (1 + 5*Si) where Si is ice saturation.
    pub fn bulk_modulus_frozen(&self, t_celsius: f64) -> f64 {
        let k_base = 1.5e8; // typical unfrozen bulk modulus Pa
        let si = self.ice_saturation(t_celsius);
        k_base * (1.0 + 5.0 * si)
    }

    /// Thermal conductivity λ \[W/(m·K)\] using geometric mean model.
    ///
    /// λ = λ_s^(1-n) * λ_w^(θ_u) * λ_i^(θ_i)
    pub fn thermal_conductivity(&self, theta_u: f64) -> f64 {
        let n = 0.4; // porosity
        let lambda_water = 0.56_f64;
        let lambda_ice = 2.2_f64;
        let theta_i = (n - theta_u).max(0.0);
        self.lambda_solid.powf(1.0 - n) * lambda_water.powf(theta_u) * lambda_ice.powf(theta_i)
    }

    /// Heat capacity \[J/(m^3·K)\] — weighted by volume fractions.
    pub fn volumetric_heat_capacity(&self, t_celsius: f64) -> f64 {
        let theta_u = self.unfrozen_water_content(t_celsius);
        let theta_i = (self.total_water_content - theta_u).max(0.0);
        let n = 0.4;
        // Solid: ~2.0e6, Water: 4.18e6, Ice: 1.9e6 J/(m^3 K)
        2.0e6 * (1.0 - n) + 4.18e6 * theta_u + 1.9e6 * theta_i
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Seabed Sediment Model
// ─────────────────────────────────────────────────────────────────────────────

/// Seabed sediment model for offshore geotechnical applications.
///
/// Includes undrained shear strength profiles and consolidation characteristics.
#[derive(Debug, Clone)]
pub struct SeabedSedimentModel {
    /// Undrained shear strength at mudline (z=0) \[Pa\].
    pub su0: f64,
    /// Strength gradient with depth \[Pa/m\].
    pub k_su: f64,
    /// Compression index Cc \[-\].
    pub compression_index_cc: f64,
    /// Recompression index Cr \[-\].
    pub recompression_index_cr: f64,
    /// Coefficient of consolidation Cv \[m^2/s\].
    pub cv: f64,
    /// Initial void ratio e0 \[-\].
    pub e0: f64,
    /// Saturated unit weight γ_sat \[kN/m^3\].
    pub gamma_sat: f64,
}

impl SeabedSedimentModel {
    /// Create a new seabed sediment model.
    pub fn new(
        su0: f64,
        k_su: f64,
        compression_index_cc: f64,
        recompression_index_cr: f64,
        cv: f64,
        e0: f64,
        gamma_sat: f64,
    ) -> Self {
        Self {
            su0,
            k_su,
            compression_index_cc,
            recompression_index_cr,
            cv,
            e0,
            gamma_sat,
        }
    }

    /// Typical soft marine clay (Gulf of Mexico).
    pub fn gulf_of_mexico_clay() -> Self {
        Self::new(5e3, 1.5e3, 0.7, 0.07, 2e-8, 2.8, 16.0)
    }

    /// Undrained shear strength Su at depth z \[m\] \[Pa\].
    ///
    /// Su(z) = su0 + k_su * z
    pub fn undrained_shear_strength(&self, depth_m: f64) -> f64 {
        self.su0 + self.k_su * depth_m
    }

    /// Primary consolidation settlement ΔS \[m\].
    ///
    /// ΔS = (Cc / (1 + e0)) * H * log10(σ_vf / σ_v0)
    pub fn settlement_primary(&self, sigma_v0: f64, sigma_vf: f64, h: f64, e0: f64) -> f64 {
        if sigma_v0 <= 0.0 || sigma_vf <= sigma_v0 {
            return 0.0;
        }
        (self.compression_index_cc / (1.0 + e0)) * h * (sigma_vf / sigma_v0).log10()
    }

    /// Time to reach a given degree of consolidation U (Terzaghi 1-D).
    ///
    /// T_v = Cv * t / H_dr^2, solved for t:
    /// t = T_v * H_dr^2 / Cv
    pub fn time_to_consolidation(&self, cv: f64, h_dr: f64, u: f64) -> f64 {
        // Time factor T_v from degree of consolidation U
        let t_v = if u <= 0.6 {
            PI * PI * u * u / 4.0
        } else {
            -0.933 * (1.0 - u).log10() - 0.085
        };
        t_v * h_dr * h_dr / cv
    }

    /// Overconsolidation ratio OCR at depth z given preconsolidation pressure.
    pub fn ocr(&self, sigma_v: f64, sigma_p: f64) -> f64 {
        if sigma_v <= 0.0 {
            return 1.0;
        }
        sigma_p / sigma_v
    }

    /// Pore pressure ratio ru = u / σ'v.
    pub fn pore_pressure_ratio(&self, u_excess: f64, sigma_v_eff: f64) -> f64 {
        if sigma_v_eff <= 0.0 {
            return 0.0;
        }
        u_excess / sigma_v_eff
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Granular Pressure Model — Hertz-Mindlin Contact
// ─────────────────────────────────────────────────────────────────────────────

/// Granular pressure model using Hertz-Mindlin contact theory.
///
/// Models pressure-dependent stiffness in granular assemblies (e.g., sand).
#[derive(Debug, Clone)]
pub struct GranularPressureModel {
    /// Shear modulus of grain material G \[Pa\].
    pub grain_shear_modulus: f64,
    /// Poisson's ratio of grain material ν.
    pub grain_poisson_ratio: f64,
    /// Mean grain radius R \[m\].
    pub mean_grain_radius: f64,
    /// Current void ratio e \[-\].
    pub void_ratio: f64,
}

impl GranularPressureModel {
    /// Create a new granular pressure model.
    pub fn new(
        grain_shear_modulus: f64,
        grain_poisson_ratio: f64,
        mean_grain_radius: f64,
        void_ratio: f64,
    ) -> Self {
        Self {
            grain_shear_modulus,
            grain_poisson_ratio,
            mean_grain_radius,
            void_ratio,
        }
    }

    /// Typical quartz sand grains.
    pub fn quartz_sand() -> Self {
        Self::new(44e9, 0.07, 0.15e-3, 0.65)
    }

    /// Hertz-Mindlin contact stiffnesses (Kn, Kt) for two identical spheres.
    ///
    /// Kn = (4/3) * G_eff * sqrt(R_eff * δ)  \[N/m\]
    /// Kt = 8 * G_eff * sqrt(R_eff * δ)       \[N/m\]
    ///
    /// Returns (Kn, Kt).
    pub fn hertz_mindlin_stiffness(&self, g: f64, nu: f64, r_eff: f64, sigma_n: f64) -> (f64, f64) {
        let g_eff = g / (2.0 * (1.0 - nu));
        // Normal overlap from Hertz theory: σ_n = Kn_HM / A * force
        // Contact stiffness from pressure: use Walton-Braun formulation
        let r_eff_m = r_eff;
        // δ from σ_n: F = (4/3)*G_eff*sqrt(R)*δ^(3/2) → δ = (3F/(4*G_eff*sqrt(R)))^(2/3)
        // F ≈ σ_n * A, A ≈ π*R^2 per grain
        let f_n = sigma_n * PI * r_eff_m * r_eff_m;
        let delta = (3.0 * f_n / (4.0 * g_eff * r_eff_m.sqrt()))
            .max(0.0)
            .powf(2.0 / 3.0);
        let kn = (4.0 / 3.0) * g_eff * (r_eff_m * delta).sqrt();
        let kt = 8.0 * g_eff * (r_eff_m * delta).sqrt();
        (kn, kt)
    }

    /// Coordination number Z from void ratio (empirical).
    ///
    /// Z = 10.726 - 11.469 * e (for e ∈ \[0.5, 0.9\])
    pub fn coordination_number(&self, void_ratio: f64) -> f64 {
        (10.726 - 11.469 * void_ratio).max(3.0)
    }

    /// Pack (granular assembly) elastic shear modulus \[Pa\].
    ///
    /// Uses Walton (1987) effective medium theory:
    /// G_pack = (5-4ν)/(5*(2-ν)) * (3*Z^2*(1-e)^2*G^2/(2*π^2*(1+e)^2))^(1/3) * p^(1/3)
    pub fn pack_stiffness(&self, g: f64, nu: f64, e: f64, sigma: f64) -> f64 {
        let z = self.coordination_number(e);
        let factor = (5.0 - 4.0 * nu) / (5.0 * (2.0 - nu));
        let inner =
            3.0 * z * z * (1.0 - e) * (1.0 - e) * g * g / (2.0 * PI * PI * (1.0 + e) * (1.0 + e));
        factor * inner.powf(1.0 / 3.0) * sigma.abs().powf(1.0 / 3.0)
    }

    /// Velocity of shear waves in granular pack \[m/s\].
    pub fn shear_wave_velocity(&self, rho_bulk: f64, sigma: f64) -> f64 {
        let g_pack = self.pack_stiffness(
            self.grain_shear_modulus,
            self.grain_poisson_ratio,
            self.void_ratio,
            sigma,
        );
        (g_pack / rho_bulk).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Loess Model
// ─────────────────────────────────────────────────────────────────────────────

/// Loess (wind-deposited silty soil) model with collapsible behavior.
#[derive(Debug, Clone)]
pub struct LoessModel {
    /// Natural void ratio.
    pub natural_void_ratio: f64,
    /// Liquid limit LL \[%\].
    pub liquid_limit: f64,
    /// Plastic limit PL \[%\].
    pub plastic_limit: f64,
    /// Natural water content w \[%\].
    pub water_content: f64,
    /// Saturated hydraulic conductivity \[m/s\].
    pub k_sat: f64,
}

impl LoessModel {
    /// Create a new loess model.
    pub fn new(
        natural_void_ratio: f64,
        liquid_limit: f64,
        plastic_limit: f64,
        water_content: f64,
        k_sat: f64,
    ) -> Self {
        Self {
            natural_void_ratio,
            liquid_limit,
            plastic_limit,
            water_content,
            k_sat,
        }
    }

    /// Typical Chinese Loess Plateau loess.
    pub fn chinese_loess() -> Self {
        Self::new(1.05, 28.0, 16.0, 12.0, 1e-6)
    }

    /// Plasticity index PI = LL - PL.
    pub fn plasticity_index(&self) -> f64 {
        self.liquid_limit - self.plastic_limit
    }

    /// Liquidity index LI = (w - PL) / PI.
    pub fn liquidity_index(&self) -> f64 {
        let pi = self.plasticity_index();
        if pi <= 0.0 {
            return 0.0;
        }
        (self.water_content - self.plastic_limit) / pi
    }

    /// Collapsibility coefficient δs (coefficient of subsidence).
    ///
    /// δs ≈ (e_p - e_s) / (1 + e0)  where e_p (natural) and e_s (saturated at same pressure).
    pub fn collapsibility_coefficient(&self, e_saturated: f64) -> f64 {
        (self.natural_void_ratio - e_saturated) / (1.0 + self.natural_void_ratio)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Rockfill Dam Material
// ─────────────────────────────────────────────────────────────────────────────

/// Rockfill and embankment dam core material model (Duncan-Chang hyperbolic).
#[derive(Debug, Clone)]
pub struct RockfillModel {
    /// Tangential modulus number K \[-\].
    pub k_modulus: f64,
    /// Modulus exponent n \[-\].
    pub n_exponent: f64,
    /// Failure ratio Rf \[-\] (≈ 0.7–0.9).
    pub rf: f64,
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
    /// Friction angle φ \[radians\].
    pub friction_angle: f64,
    /// Atmospheric pressure pa \[Pa\].
    pub pa: f64,
}

impl RockfillModel {
    /// Create a new rockfill model.
    pub fn new(
        k_modulus: f64,
        n_exponent: f64,
        rf: f64,
        cohesion: f64,
        friction_angle_deg: f64,
    ) -> Self {
        Self {
            k_modulus,
            n_exponent,
            rf,
            cohesion,
            friction_angle: friction_angle_deg.to_radians(),
            pa: 101325.0,
        }
    }

    /// Typical hard rockfill (quartzite).
    pub fn hard_quartzite_rockfill() -> Self {
        Self::new(800.0, 0.4, 0.8, 0.0, 40.0)
    }

    /// Initial tangent modulus Ei \[Pa\] (Duncan-Chang).
    ///
    /// Ei = K * pa * (σ3/pa)^n
    pub fn initial_tangent_modulus(&self, sigma3: f64) -> f64 {
        self.k_modulus * self.pa * (sigma3 / self.pa).abs().powf(self.n_exponent)
    }

    /// Deviatoric stress at failure (Mohr-Coulomb).
    ///
    /// (σ1 - σ3)_f = 2*(c*cos(φ) + σ3*sin(φ)) / (1 - sin(φ))
    pub fn deviatoric_stress_at_failure(&self, sigma3: f64) -> f64 {
        let phi = self.friction_angle;
        2.0 * (self.cohesion * phi.cos() + sigma3 * phi.sin()) / (1.0 - phi.sin())
    }

    /// Tangent modulus Et \[Pa\] from current stress ratio.
    ///
    /// Et = Ei * (1 - Rf * (σ1 - σ3) / (σ1 - σ3)_f)^2
    pub fn tangent_modulus(&self, sigma1: f64, sigma3: f64) -> f64 {
        let ei = self.initial_tangent_modulus(sigma3);
        let delta_f = self.deviatoric_stress_at_failure(sigma3);
        if delta_f <= 0.0 {
            return ei;
        }
        let delta = sigma1 - sigma3;
        let ratio = self.rf * delta / delta_f;
        let factor = (1.0 - ratio).max(0.0);
        ei * factor * factor
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── SoilModel ─────────────────────────────────────────────────────────

    #[test]
    fn test_soil_yield_function_elastic() {
        let soil = SoilModel::medium_sand();
        // Small stresses should be elastic
        let f = soil.yield_function(50e3, 100e3); // σ1 < σ3 (not realistic but tests sign)
        // With σ1 < σ3 the yield function should be negative
        assert!(f < 0.0, "Should be elastic when σ1 < σ3: {f}");
    }

    #[test]
    fn test_soil_yield_function_at_failure() {
        let soil = SoilModel::stiff_clay();
        // Compute σ1 that causes failure for σ3 = 100 kPa
        let sigma3 = 100e3;
        let phi = soil.friction_angle;
        let c = soil.cohesion;
        // (σ1-σ3) = 2c*cos(φ) + (σ1+σ3)*sin(φ) → solve for σ1
        // σ1*(1-sin(φ)) = 2c*cos(φ) + σ3*(1+sin(φ))
        let sigma1 = (2.0 * c * phi.cos() + sigma3 * (1.0 + phi.sin())) / (1.0 - phi.sin());
        let f = soil.yield_function(sigma1, sigma3);
        assert!(
            f.abs() < 1e-6 * sigma1,
            "Yield function should be ~0 at failure: {f}"
        );
    }

    #[test]
    fn test_soil_bulk_modulus() {
        let soil = SoilModel::medium_sand();
        let k = soil.bulk_modulus();
        assert!(k > 0.0, "Bulk modulus must be positive");
    }

    #[test]
    fn test_soil_shear_modulus() {
        let soil = SoilModel::medium_sand();
        let g = soil.shear_modulus();
        assert!(g > 0.0);
    }

    #[test]
    fn test_soil_ucs() {
        let soil = SoilModel::stiff_clay();
        let qu = soil.ucs();
        // For c=60kPa, φ=22°: qu = 2*60000*cos(22°)/(1-sin(22°)) ≈ 180 kPa
        assert!(qu > 100e3 && qu < 300e3, "UCS {qu} Pa out of range");
    }

    #[test]
    fn test_soil_k0() {
        let soil = SoilModel::medium_sand();
        let k0 = soil.k0_jaky();
        assert!(k0 > 0.0 && k0 < 1.0, "Ko must be in (0,1): {k0}");
    }

    #[test]
    fn test_soil_critical_state_q() {
        let soil = SoilModel::medium_sand();
        let q = soil.critical_state_line_q(200e3);
        assert!(q > 0.0);
    }

    #[test]
    fn test_soil_plastic_multiplier() {
        let soil = SoilModel::medium_sand();
        let lambda = soil.plastic_multiplier(0.01);
        assert!(lambda > 0.0 && lambda < 0.01);
    }

    // ── CamClayModel ──────────────────────────────────────────────────────

    #[test]
    fn test_camclay_e_cs() {
        let model = CamClayModel::kaolin_ncc();
        let e = model.e_cs(100e3); // p = 100 kPa
        assert!(e > 0.0, "Void ratio on CSL must be positive: {e}");
    }

    #[test]
    fn test_camclay_yield_surface_inside() {
        let model = CamClayModel::kaolin_ncc();
        // Inside yield surface (isotropic state well below preconsolidation)
        let f = model.yield_surface_p0(50e3, 0.0);
        assert!(f < 0.0, "Should be inside yield surface at p=50kPa<p0: {f}");
    }

    #[test]
    fn test_camclay_yield_surface_on() {
        let model = CamClayModel::kaolin_ncc();
        // On yield surface: q=0, p=p0
        let f = model.yield_surface_p0(model.p_c0, 0.0);
        assert!(
            f.abs() < 1e-6 * model.p_c0 * model.p_c0,
            "On yield surface: {f}"
        );
    }

    #[test]
    fn test_camclay_preconsolidation_pressure() {
        let model = CamClayModel::kaolin_ncc();
        let p0 = model.preconsolidation_pressure(model.e0);
        assert!(p0 > 0.0, "Preconsolidation pressure must be positive: {p0}");
    }

    #[test]
    fn test_camclay_elastic_volumetric_strain() {
        let model = CamClayModel::kaolin_ncc();
        let deps = model.volumetric_strain_elastic(10e3, 100e3);
        assert!(deps > 0.0, "Elastic strain must be positive");
    }

    #[test]
    fn test_camclay_plastic_strain() {
        let model = CamClayModel::kaolin_ncc();
        let deps = model.volumetric_strain_plastic(10e3, 100e3);
        assert!(deps > 0.0);
    }

    // ── RockModel ─────────────────────────────────────────────────────────

    #[test]
    fn test_rock_mb_granite() {
        let rock = RockModel::granite();
        let mb = rock.mb_parameter();
        assert!(mb > 0.0 && mb < rock.mi, "mb < mi: {mb}");
    }

    #[test]
    fn test_rock_s_parameter() {
        let rock = RockModel::granite();
        let s = rock.s_parameter();
        assert!(s > 0.0 && s <= 1.0, "s ∈ (0,1]: {s}");
    }

    #[test]
    fn test_rock_a_parameter() {
        let rock = RockModel::granite();
        let a = rock.a_parameter();
        assert!((0.5..0.65).contains(&a), "a ≈ 0.5 for intact rock: {a}");
    }

    #[test]
    fn test_rock_ucs() {
        let rock = RockModel::granite();
        let ucs_mass = rock.uniaxial_compressive_strength_rock();
        assert!(
            ucs_mass > 0.0 && ucs_mass <= rock.ucs,
            "Rock mass UCS ≤ intact UCS"
        );
    }

    #[test]
    fn test_rock_tensile_strength() {
        let rock = RockModel::granite();
        let t = rock.tensile_strength_rock();
        assert!(
            t < 0.0,
            "Tensile strength should be negative (tension): {t}"
        );
    }

    #[test]
    fn test_rock_hoek_brown_yield() {
        let rock = RockModel::granite();
        // At sigma3=0, the HB strength = ucs_mass; sigma1 = 2*ucs_mass should yield
        let sigma3 = 0.0;
        let ucs_mass = rock.uniaxial_compressive_strength_rock();
        let sigma1 = ucs_mass * 2.0;
        let yielding = rock.yield_hoek_brown(sigma1, sigma3);
        assert!(
            yielding,
            "Should yield at sigma1=2*ucs_mass, sigma3=0: ucs_mass={ucs_mass:.3e}"
        );
    }

    #[test]
    fn test_rock_deformation_modulus() {
        let rock = RockModel::granite();
        let em = rock.deformation_modulus();
        assert!(em > 0.0, "Deformation modulus must be positive");
    }

    #[test]
    fn test_rock_shale_params() {
        let rock = RockModel::shale();
        let mb = rock.mb_parameter();
        let mb_granite = RockModel::granite().mb_parameter();
        assert!(mb < mb_granite, "Shale has lower mb than granite");
    }

    // ── SandModel ─────────────────────────────────────────────────────────

    #[test]
    fn test_sand_peak_friction_angle() {
        let sand = SandModel::medium_dense_quartz();
        let phi_p = sand.peak_friction_angle(0.6);
        assert!(phi_p >= sand.critical_state_friction_angle, "Peak φ ≥ CS φ");
    }

    #[test]
    fn test_sand_dilatancy_index_positive() {
        let sand = SandModel::medium_dense_quartz();
        let ir = sand.dilatancy_index(0.7, 100.0);
        assert!(ir >= 0.0, "Dilatancy index non-negative");
    }

    #[test]
    fn test_sand_bolton_dilatancy() {
        let sand = SandModel::medium_dense_quartz();
        let ir = sand.bolton_dilatancy(0.6, 100.0);
        assert!(ir >= 0.0);
    }

    #[test]
    fn test_sand_void_ratio() {
        let sand = SandModel::medium_dense_quartz();
        let e = sand.current_void_ratio();
        assert!(
            e >= sand.e_min() && e <= sand.e_max(),
            "Void ratio out of range: {e}"
        );
    }

    #[test]
    fn test_sand_dry_unit_weight() {
        let sand = SandModel::medium_dense_quartz();
        let gd = sand.dry_unit_weight();
        assert!(
            gd > 0.0 && gd < 30.0,
            "Dry unit weight {gd} kN/m^3 out of range"
        );
    }

    #[test]
    fn test_sand_loose_has_lower_dr() {
        let loose = SandModel::loose_sand();
        let dense = SandModel::medium_dense_quartz();
        assert!(loose.relative_density < dense.relative_density);
    }

    // ── PermafrostModel ───────────────────────────────────────────────────

    #[test]
    fn test_permafrost_unfrozen_content_above_zero() {
        let pf = PermafrostModel::siberian_silt();
        let theta = pf.unfrozen_water_content(5.0);
        assert_eq!(
            theta, pf.total_water_content,
            "Above 0°C: all water unfrozen"
        );
    }

    #[test]
    fn test_permafrost_unfrozen_content_frozen() {
        let pf = PermafrostModel::siberian_silt();
        let theta = pf.unfrozen_water_content(-10.0);
        assert!(
            theta < pf.total_water_content,
            "Below 0°C: some water frozen"
        );
        assert!(theta >= 0.0);
    }

    #[test]
    fn test_permafrost_ice_saturation_above_zero() {
        let pf = PermafrostModel::siberian_silt();
        let si = pf.ice_saturation(5.0);
        assert_eq!(si, 0.0, "No ice above 0°C");
    }

    #[test]
    fn test_permafrost_ice_saturation_frozen() {
        let pf = PermafrostModel::siberian_silt();
        let si = pf.ice_saturation(-20.0);
        assert!(si > 0.0 && si <= 1.0, "Ice saturation in (0,1]: {si}");
    }

    #[test]
    fn test_permafrost_bulk_modulus_increases_with_cooling() {
        let pf = PermafrostModel::siberian_silt();
        let k_warm = pf.bulk_modulus_frozen(-1.0);
        let k_cold = pf.bulk_modulus_frozen(-20.0);
        assert!(k_cold > k_warm, "Frozen soil stiffens with cooling");
    }

    #[test]
    fn test_permafrost_thermal_conductivity() {
        let pf = PermafrostModel::siberian_silt();
        let lam = pf.thermal_conductivity(0.1);
        assert!(
            lam > 0.0 && lam < 10.0,
            "Conductivity {lam} W/(mK) out of range"
        );
    }

    #[test]
    fn test_permafrost_heat_capacity() {
        let pf = PermafrostModel::siberian_silt();
        let c = pf.volumetric_heat_capacity(-10.0);
        assert!(
            c > 1e5 && c < 5e6,
            "Heat capacity {c} J/(m^3 K) out of range"
        );
    }

    // ── SeabedSedimentModel ───────────────────────────────────────────────

    #[test]
    fn test_seabed_su_at_mudline() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let su = sed.undrained_shear_strength(0.0);
        assert_eq!(su, sed.su0, "At z=0, Su = su0");
    }

    #[test]
    fn test_seabed_su_increases_with_depth() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let su1 = sed.undrained_shear_strength(5.0);
        let su2 = sed.undrained_shear_strength(10.0);
        assert!(su2 > su1, "Su increases with depth");
    }

    #[test]
    fn test_seabed_settlement_positive() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let ds = sed.settlement_primary(50e3, 100e3, 10.0, sed.e0);
        assert!(ds > 0.0, "Settlement must be positive");
    }

    #[test]
    fn test_seabed_settlement_zero_when_no_loading() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let ds = sed.settlement_primary(50e3, 50e3, 10.0, sed.e0);
        assert_eq!(ds, 0.0, "No settlement without additional load");
    }

    #[test]
    fn test_seabed_time_consolidation() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let t = sed.time_to_consolidation(2e-8, 5.0, 0.5);
        assert!(t > 0.0, "Consolidation time must be positive");
    }

    #[test]
    fn test_seabed_ocr() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let ocr = sed.ocr(50e3, 100e3);
        assert!((ocr - 2.0).abs() < 1e-10, "OCR = 2 for sigma_p = 2*sigma_v");
    }

    // ── GranularPressureModel ─────────────────────────────────────────────

    #[test]
    fn test_granular_hertz_mindlin_stiffness() {
        let model = GranularPressureModel::quartz_sand();
        let (kn, kt) = model.hertz_mindlin_stiffness(44e9, 0.07, 0.15e-3, 100e3);
        assert!(kn > 0.0, "Normal stiffness must be positive: {kn}");
        assert!(kt > 0.0, "Tangential stiffness must be positive: {kt}");
    }

    #[test]
    fn test_granular_hertz_mindlin_kn_gt_kt_ratio() {
        let model = GranularPressureModel::quartz_sand();
        let (kn, kt) = model.hertz_mindlin_stiffness(44e9, 0.07, 0.15e-3, 100e3);
        // kt/kn = 6 from HM theory
        let ratio = kt / kn;
        assert!((ratio - 6.0).abs() < 1e-6, "kt/kn should be 6: {ratio}");
    }

    #[test]
    fn test_granular_coordination_number() {
        let model = GranularPressureModel::quartz_sand();
        let z = model.coordination_number(0.65);
        assert!(z > 3.0 && z < 12.0, "Coordination number {z} out of range");
    }

    #[test]
    fn test_granular_pack_stiffness_positive() {
        let model = GranularPressureModel::quartz_sand();
        let g_pack = model.pack_stiffness(44e9, 0.07, 0.65, 100e3);
        assert!(g_pack > 0.0, "Pack stiffness must be positive");
    }

    #[test]
    fn test_granular_pack_stiffness_pressure_dependent() {
        let model = GranularPressureModel::quartz_sand();
        let g1 = model.pack_stiffness(44e9, 0.07, 0.65, 100e3);
        let g2 = model.pack_stiffness(44e9, 0.07, 0.65, 400e3);
        assert!(g2 > g1, "Pack stiffness increases with pressure");
    }

    #[test]
    fn test_granular_shear_wave_velocity() {
        let model = GranularPressureModel::quartz_sand();
        let rho = 1600.0; // kg/m^3
        let vs = model.shear_wave_velocity(rho, 100e3);
        assert!(vs > 0.0 && vs < 2000.0, "Vs {vs} m/s out of range");
    }

    // ── LoessModel ────────────────────────────────────────────────────────

    #[test]
    fn test_loess_plasticity_index() {
        let loess = LoessModel::chinese_loess();
        let pi = loess.plasticity_index();
        assert!((pi - 12.0).abs() < 1e-10, "PI = LL - PL = 12");
    }

    #[test]
    fn test_loess_liquidity_index() {
        let loess = LoessModel::chinese_loess();
        let li = loess.liquidity_index();
        // LI = (12 - 16)/12 = -4/12 < 0 (stiff)
        assert!(li < 0.0, "LI should be negative for stiff loess");
    }

    #[test]
    fn test_loess_collapsibility() {
        let loess = LoessModel::chinese_loess();
        let delta_s = loess.collapsibility_coefficient(0.85);
        assert!(delta_s > 0.0, "Collapsibility must be positive");
    }

    // ── RockfillModel ─────────────────────────────────────────────────────

    #[test]
    fn test_rockfill_initial_modulus() {
        let rf = RockfillModel::hard_quartzite_rockfill();
        let ei = rf.initial_tangent_modulus(200e3);
        assert!(ei > 0.0, "Initial modulus must be positive");
    }

    #[test]
    fn test_rockfill_tangent_modulus_decreases() {
        let rf = RockfillModel::hard_quartzite_rockfill();
        let sigma3 = 200e3;
        let ei = rf.initial_tangent_modulus(sigma3);
        let delta_f = rf.deviatoric_stress_at_failure(sigma3);
        // At half failure stress
        let sigma1_half = sigma3 + 0.5 * delta_f;
        let et = rf.tangent_modulus(sigma1_half, sigma3);
        assert!(et < ei, "Tangent modulus decreases with shear");
    }

    #[test]
    fn test_rockfill_deviatoric_stress_at_failure_positive() {
        let rf = RockfillModel::hard_quartzite_rockfill();
        let dsf = rf.deviatoric_stress_at_failure(200e3);
        assert!(dsf > 0.0);
    }

    // ── Cross-model consistency tests ─────────────────────────────────────

    #[test]
    fn test_camclay_e_cs_decreases_with_pressure() {
        let model = CamClayModel::kaolin_ncc();
        let e1 = model.e_cs(100e3);
        let e2 = model.e_cs(1000e3);
        assert!(e2 < e1, "Void ratio on CSL decreases with pressure");
    }

    #[test]
    fn test_rock_granite_stronger_than_shale() {
        let granite = RockModel::granite();
        let shale = RockModel::shale();
        assert!(
            granite.uniaxial_compressive_strength_rock()
                > shale.uniaxial_compressive_strength_rock(),
            "Granite stronger than shale"
        );
    }

    #[test]
    fn test_soil_clay_has_cohesion() {
        let clay = SoilModel::stiff_clay();
        assert!(clay.cohesion > 0.0, "Clay has cohesion");
    }

    #[test]
    fn test_sand_has_no_cohesion() {
        let sand = SoilModel::medium_sand();
        assert_eq!(sand.cohesion, 0.0, "Sand has no cohesion");
    }

    #[test]
    fn test_seabed_pore_pressure_ratio() {
        let sed = SeabedSedimentModel::gulf_of_mexico_clay();
        let ru = sed.pore_pressure_ratio(20e3, 100e3);
        assert!((ru - 0.2).abs() < 1e-10);
    }
}
