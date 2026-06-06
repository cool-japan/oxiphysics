// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Additive manufacturing material models.
//!
//! This module covers the process-structure-property (PSP) chain for key AM
//! technologies:
//!
//! - **Powder bed fusion** (SLM / DMLS): [`PbfMaterial`], [`PbfProcessParams`],
//!   [`MeltPoolGeometry`], thermal gradient and residual stress computation.
//! - **Porosity effects**: [`PorosityModel`] — Ramakrishnan-Arunachalam,
//!   Mori-Tanaka, Coble-Kingery relationships.
//! - **Microstructure evolution**: [`GrainGrowthModel`], [`MartensiticModel`].
//! - **Binder jetting**: [`BinderJettingMaterial`], sintering shrinkage.
//! - **FDM polymers**: [`FdmPolymerMaterial`] — anisotropy, layer adhesion.
//! - **Support structures**: [`SupportMaterial`].
//! - **PSP linkage**: [`PspLinkage`] combining all above.
//! - **Scan strategy effects**: [`ScanStrategyEffect`].
//! - **Surface roughness**: [`AmSurfaceRoughness`].

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J/K).
pub const KB: f64 = 1.380_649e-23;

/// Gas constant (J/(mol·K)).
pub const R_GAS: f64 = 8.314_462;

// ---------------------------------------------------------------------------
// Common enums
// ---------------------------------------------------------------------------

/// AM process technology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmProcess {
    /// Selective laser melting (SLM) / laser powder bed fusion (LPBF).
    Slm,
    /// Direct metal laser sintering (DMLS).
    Dmls,
    /// Electron beam melting (EBM).
    Ebm,
    /// Binder jetting (BJ).
    BinderJetting,
    /// Fused deposition modelling (FDM) / fused filament fabrication (FFF).
    Fdm,
    /// Directed energy deposition (DED).
    Ded,
}

/// Metal alloy system commonly used in powder bed fusion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetalAlloy {
    /// Ti-6Al-4V (Grade 5 titanium).
    Ti6Al4V,
    /// 316L stainless steel.
    Steel316L,
    /// AlSi10Mg aluminium alloy.
    AlSi10Mg,
    /// IN718 nickel superalloy.
    In718,
    /// Maraging steel (1.2709).
    MaragingSteel,
    /// Cobalt-chromium (CoCr).
    CoCr,
    /// Copper (Cu).
    Copper,
}

// ---------------------------------------------------------------------------
// Powder Bed Fusion Material
// ---------------------------------------------------------------------------

/// Intrinsic material properties relevant to PBF processing.
#[derive(Debug, Clone)]
pub struct PbfMaterial {
    /// Alloy type.
    pub alloy: MetalAlloy,
    /// Density of fully dense material (kg/m³).
    pub density: f64,
    /// Thermal conductivity at room temperature (W/(m·K)).
    pub thermal_conductivity: f64,
    /// Specific heat capacity (J/(kg·K)).
    pub specific_heat: f64,
    /// Latent heat of fusion (J/kg).
    pub latent_heat_fusion: f64,
    /// Solidus temperature (K).
    pub solidus_temp: f64,
    /// Liquidus temperature (K).
    pub liquidus_temp: f64,
    /// Absorptivity for the process laser wavelength (0–1).
    pub absorptivity: f64,
    /// Young's modulus at room temperature (Pa).
    pub elastic_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Yield strength at room temperature (Pa).
    pub yield_strength: f64,
    /// Tensile strength at room temperature (Pa).
    pub tensile_strength: f64,
    /// Thermal expansion coefficient (1/K).
    pub thermal_expansion: f64,
    /// Powder bed packing fraction (0–1).
    pub packing_fraction: f64,
}

impl PbfMaterial {
    /// Pre-defined Ti-6Al-4V properties for SLM.
    pub fn ti6al4v_slm() -> Self {
        Self {
            alloy: MetalAlloy::Ti6Al4V,
            density: 4_430.0,
            thermal_conductivity: 6.7,
            specific_heat: 560.0,
            latent_heat_fusion: 286_000.0,
            solidus_temp: 1878.0,
            liquidus_temp: 1928.0,
            absorptivity: 0.35,
            elastic_modulus: 114e9,
            poisson_ratio: 0.342,
            yield_strength: 930e6,
            tensile_strength: 1_000e6,
            thermal_expansion: 8.6e-6,
            packing_fraction: 0.60,
        }
    }

    /// Pre-defined 316L stainless steel properties for SLM.
    pub fn steel_316l_slm() -> Self {
        Self {
            alloy: MetalAlloy::Steel316L,
            density: 7_990.0,
            thermal_conductivity: 16.3,
            specific_heat: 500.0,
            latent_heat_fusion: 272_000.0,
            solidus_temp: 1648.0,
            liquidus_temp: 1673.0,
            absorptivity: 0.40,
            elastic_modulus: 193e9,
            poisson_ratio: 0.290,
            yield_strength: 530e6,
            tensile_strength: 680e6,
            thermal_expansion: 16.0e-6,
            packing_fraction: 0.62,
        }
    }

    /// Pre-defined AlSi10Mg properties for SLM.
    pub fn alsi10mg_slm() -> Self {
        Self {
            alloy: MetalAlloy::AlSi10Mg,
            density: 2_680.0,
            thermal_conductivity: 130.0,
            specific_heat: 910.0,
            latent_heat_fusion: 396_000.0,
            solidus_temp: 833.0,
            liquidus_temp: 868.0,
            absorptivity: 0.09,
            elastic_modulus: 72e9,
            poisson_ratio: 0.330,
            yield_strength: 240e6,
            tensile_strength: 330e6,
            thermal_expansion: 21.0e-6,
            packing_fraction: 0.62,
        }
    }

    /// Thermal diffusivity α = k / (ρ c_p) (m²/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Melting range ΔT = T_liq − T_sol (K).
    pub fn melting_range(&self) -> f64 {
        (self.liquidus_temp - self.solidus_temp).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// PBF Process Parameters
// ---------------------------------------------------------------------------

/// Laser / electron beam process parameters for powder bed fusion.
#[derive(Debug, Clone)]
pub struct PbfProcessParams {
    /// Laser power (W).
    pub power: f64,
    /// Scan speed (m/s).
    pub scan_speed: f64,
    /// Hatch spacing (m).
    pub hatch_spacing: f64,
    /// Layer thickness (m).
    pub layer_thickness: f64,
    /// Laser spot radius (1/e² radius) (m).
    pub spot_radius: f64,
    /// Preheat temperature of the build plate (K).
    pub preheat_temp: f64,
    /// Scan strategy rotation angle between layers (degrees).
    pub rotation_angle_deg: f64,
}

impl PbfProcessParams {
    /// Volumetric energy density (J/m³): `E_v = P / (v * h * t)`.
    pub fn volumetric_energy_density(&self) -> f64 {
        let denom = self.scan_speed * self.hatch_spacing * self.layer_thickness;
        if denom < 1e-30 {
            f64::INFINITY
        } else {
            self.power / denom
        }
    }

    /// Linear energy density (J/m): `E_l = P / v`.
    pub fn linear_energy_density(&self) -> f64 {
        if self.scan_speed < 1e-12 {
            f64::INFINITY
        } else {
            self.power / self.scan_speed
        }
    }

    /// Interaction time of laser with powder (s): `t_int = 2*r / v`.
    pub fn interaction_time(&self) -> f64 {
        if self.scan_speed < 1e-12 {
            f64::INFINITY
        } else {
            2.0 * self.spot_radius / self.scan_speed
        }
    }
}

// ---------------------------------------------------------------------------
// Melt Pool Geometry (Rosenthal / Goldak approximation)
// ---------------------------------------------------------------------------

/// Estimated melt pool geometry from process parameters and material properties.
///
/// Based on the Eagar-Tsai analytical model (simplified).
#[derive(Debug, Clone)]
pub struct MeltPoolGeometry {
    /// Melt pool half-width (m).
    pub half_width: f64,
    /// Melt pool half-length (m) (in scan direction).
    pub half_length: f64,
    /// Melt pool depth (m).
    pub depth: f64,
}

impl MeltPoolGeometry {
    /// Estimate melt pool geometry using the simplified Eagar-Tsai model.
    ///
    /// # Arguments
    /// * `mat` — Material properties.
    /// * `params` — Process parameters.
    /// * `ambient_temp` — Ambient/preheat temperature (K).
    pub fn eagar_tsai(mat: &PbfMaterial, params: &PbfProcessParams, ambient_temp: f64) -> Self {
        let alpha = mat.thermal_diffusivity();
        let delta_t = mat.liquidus_temp - ambient_temp;
        let q = mat.absorptivity * params.power;
        let v = params.scan_speed;
        let sigma = params.spot_radius;
        // Characteristic length scales
        let r_char = ((2.0 * q * alpha) / (PI * mat.thermal_conductivity * v * delta_t))
            .sqrt()
            .max(sigma);
        let half_width = r_char;
        let half_length = r_char * (1.0 + v * r_char / (2.0 * alpha)).min(5.0);
        let depth = r_char * 0.5; // simplified: half the width
        Self {
            half_width,
            half_length,
            depth,
        }
    }

    /// Melt pool volume (m³), approximated as an ellipsoid.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.half_width * self.half_length * self.depth
    }

    /// Aspect ratio (length / width).
    pub fn aspect_ratio(&self) -> f64 {
        if self.half_width < 1e-20 {
            1.0
        } else {
            self.half_length / self.half_width
        }
    }
}

// ---------------------------------------------------------------------------
// Thermal Gradient and Cooling Rate
// ---------------------------------------------------------------------------

/// Thermal gradient and cooling rate estimates for AM.
#[derive(Debug, Clone)]
pub struct ThermalGradient {
    /// Thermal gradient magnitude at the melt pool boundary (K/m).
    pub gradient_magnitude: f64,
    /// Cooling rate at the melt pool boundary (K/s).
    pub cooling_rate: f64,
    /// Solidification rate (m/s).
    pub solidification_rate: f64,
}

impl ThermalGradient {
    /// Estimate thermal conditions from Rosenthal point-source solution.
    ///
    /// # Arguments
    /// * `mat` — Material properties.
    /// * `params` — Process parameters.
    /// * `ambient_temp` — Ambient temperature (K).
    pub fn from_rosenthal(mat: &PbfMaterial, params: &PbfProcessParams, ambient_temp: f64) -> Self {
        let alpha = mat.thermal_diffusivity();
        let v = params.scan_speed;
        let q = mat.absorptivity * params.power;
        let k = mat.thermal_conductivity;
        let r = params.spot_radius.max(1e-6);
        // Temperature at melt-pool edge (Rosenthal approximation)
        let t_peak = ambient_temp + q / (2.0 * PI * k * r) * (-v * r / (2.0 * alpha)).exp();
        let gradient_magnitude = (t_peak - ambient_temp) / r;
        let cooling_rate = gradient_magnitude * v;
        let solidification_rate = v;
        Self {
            gradient_magnitude,
            cooling_rate,
            solidification_rate,
        }
    }

    /// G·R product (K/s) — controls solidification microstructure.
    pub fn g_times_r(&self) -> f64 {
        self.gradient_magnitude * self.solidification_rate
    }

    /// G/R ratio (K·s/m²) — controls morphology (planar → cellular → dendritic).
    pub fn g_over_r(&self) -> f64 {
        if self.solidification_rate < 1e-12 {
            f64::INFINITY
        } else {
            self.gradient_magnitude / self.solidification_rate
        }
    }
}

// ---------------------------------------------------------------------------
// Residual Stress from Thermal Gradients
// ---------------------------------------------------------------------------

/// Residual stress model for PBF (temperature-gradient mechanism, TGM).
#[derive(Debug, Clone)]
pub struct ResidualStressModel {
    /// Elastic modulus (Pa).
    pub elastic_modulus: f64,
    /// Thermal expansion coefficient (1/K).
    pub thermal_expansion: f64,
    /// Yield strength at elevated temperature (Pa).
    pub yield_strength_hot: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}

impl ResidualStressModel {
    /// Build from a `PbfMaterial`.
    pub fn from_material(mat: &PbfMaterial) -> Self {
        Self {
            elastic_modulus: mat.elastic_modulus,
            thermal_expansion: mat.thermal_expansion,
            yield_strength_hot: mat.yield_strength * 0.5,
            poisson_ratio: mat.poisson_ratio,
        }
    }

    /// Estimate peak residual stress using the temperature-gradient mechanism.
    ///
    /// σ_res ≈ min(E α ΔT / (1 - ν), σ_y_hot)
    ///
    /// # Arguments
    /// * `delta_t` — Temperature difference across the heated layer (K).
    pub fn peak_residual_stress(&self, delta_t: f64) -> f64 {
        let elastic_stress =
            self.elastic_modulus * self.thermal_expansion * delta_t / (1.0 - self.poisson_ratio);
        elastic_stress.min(self.yield_strength_hot)
    }

    /// Estimate part distortion (curvature) from residual stress bilayer model.
    ///
    /// Uses Stoney's formula: κ = 6 σ t_f / (E_s t_s²)
    ///
    /// # Arguments
    /// * `sigma` — Residual stress in film (Pa).
    /// * `film_thickness` — Deposited layer thickness (m).
    /// * `substrate_thickness` — Substrate/part thickness (m).
    pub fn stoney_curvature(
        &self,
        sigma: f64,
        film_thickness: f64,
        substrate_thickness: f64,
    ) -> f64 {
        let es = self.elastic_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio);
        let ts2 = substrate_thickness * substrate_thickness;
        if ts2 < 1e-30 {
            0.0
        } else {
            6.0 * sigma * film_thickness / (es * ts2)
        }
    }

    /// Von Mises equivalent residual stress from in-plane biaxial state.
    ///
    /// For biaxial state σ_x = σ_y = σ:  σ_VM = σ.
    pub fn von_mises_biaxial(&self, sigma_inplane: f64) -> f64 {
        sigma_inplane.abs()
    }
}

// ---------------------------------------------------------------------------
// Porosity Effects on Mechanical Properties
// ---------------------------------------------------------------------------

/// Porosity model for relating relative density to effective properties.
#[derive(Debug, Clone)]
pub struct PorosityModel {
    /// Fully dense Young's modulus E₀ (Pa).
    pub e0: f64,
    /// Fully dense yield strength σ_y0 (Pa).
    pub yield_strength0: f64,
    /// Fully dense tensile strength σ_u0 (Pa).
    pub tensile_strength0: f64,
    /// Fully dense density ρ₀ (kg/m³).
    pub density0: f64,
}

impl PorosityModel {
    /// Effective Young's modulus using Ramakrishnan-Arunachalam model.
    ///
    /// E(p) = E₀ (1 - p)² / (1 + β p)  where β ≈ 2.
    pub fn elastic_modulus(&self, porosity: f64) -> f64 {
        let beta = 2.0;
        let p = porosity.clamp(0.0, 0.9999);
        self.e0 * (1.0 - p).powi(2) / (1.0 + beta * p)
    }

    /// Effective yield strength using power-law scaling.
    ///
    /// σ_y(p) = σ_y0 (1 - p)^n  with n ≈ 2.
    pub fn yield_strength(&self, porosity: f64) -> f64 {
        let n = 2.0;
        let p = porosity.clamp(0.0, 0.9999);
        self.yield_strength0 * (1.0 - p).powf(n)
    }

    /// Effective tensile strength.
    pub fn tensile_strength(&self, porosity: f64) -> f64 {
        let n = 1.8;
        let p = porosity.clamp(0.0, 0.9999);
        self.tensile_strength0 * (1.0 - p).powf(n)
    }

    /// Effective density.
    pub fn effective_density(&self, porosity: f64) -> f64 {
        self.density0 * (1.0 - porosity.clamp(0.0, 1.0))
    }

    /// Relative density (1 − porosity).
    pub fn relative_density(&self, porosity: f64) -> f64 {
        1.0 - porosity.clamp(0.0, 1.0)
    }

    /// Fatigue strength reduction factor due to porosity.
    ///
    /// Uses empirical Murakami √area model: Kt ≈ 1 + 2 √(πp/4)
    pub fn fatigue_strength_reduction(&self, porosity: f64) -> f64 {
        let p = porosity.clamp(0.0, 0.5);
        1.0 + 2.0 * (PI * p / 4.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Microstructure Evolution
// ---------------------------------------------------------------------------

/// Grain growth model for AM (isothermal / anisothermal).
///
/// Uses the parabolic grain growth law: D² − D₀² = K₀ t exp(−Q / RT).
#[derive(Debug, Clone)]
pub struct GrainGrowthModel {
    /// Pre-exponential factor K₀ (m²/s).
    pub k0: f64,
    /// Activation energy for grain boundary migration Q (J/mol).
    pub activation_energy: f64,
    /// Initial grain diameter D₀ (m).
    pub initial_diameter: f64,
}

impl GrainGrowthModel {
    /// Typical values for Ti-6Al-4V β-grain growth.
    pub fn ti6al4v_beta() -> Self {
        Self {
            k0: 7.0e-8,
            activation_energy: 175_000.0,
            initial_diameter: 50e-6,
        }
    }

    /// Grain diameter after isothermal annealing for time `t` at temperature `T` (K).
    pub fn grain_diameter_isothermal(&self, temp_k: f64, time_s: f64) -> f64 {
        let keff = self.k0 * (-self.activation_energy / (R_GAS * temp_k)).exp();
        let d2 = self.initial_diameter.powi(2) + keff * time_s;
        d2.max(0.0).sqrt()
    }

    /// Mean grain aspect ratio after directional solidification.
    ///
    /// In PBF columnar grains grow with aspect ratio AR ≈ G/R × const.
    pub fn columnar_aspect_ratio(&self, g_over_r: f64) -> f64 {
        // Empirical: AR ≈ 0.01 * (G/R) + 1; clamp to reasonable range
        (0.01 * g_over_r + 1.0).min(20.0)
    }
}

/// Martensitic transformation model for Ti alloys and steels in AM.
#[derive(Debug, Clone)]
pub struct MartensiticModel {
    /// Martensite start temperature M_s (K).
    pub ms_temp: f64,
    /// Martensite finish temperature M_f (K).
    pub mf_temp: f64,
    /// Maximum volume fraction of martensite (0–1).
    pub max_martensite_fraction: f64,
    /// Koistinen-Marburger coefficient b (K⁻¹).
    pub km_coefficient: f64,
}

impl MartensiticModel {
    /// Ti-6Al-4V α' martensite (hexagonal).
    pub fn ti6al4v_martensite() -> Self {
        Self {
            ms_temp: 875.0,
            mf_temp: 500.0,
            max_martensite_fraction: 1.0,
            km_coefficient: 0.011,
        }
    }

    /// 316L steel martensite (not typical at RT, Ms < 0).
    pub fn steel_316l_martensite() -> Self {
        Self {
            ms_temp: 233.0, // −40 °C
            mf_temp: 123.0, // −150 °C
            max_martensite_fraction: 0.30,
            km_coefficient: 0.033,
        }
    }

    /// Martensite volume fraction at temperature `T` (K), Koistinen-Marburger equation.
    ///
    /// f_m = f_max * (1 − exp(−b * (M_s − T))) for T < M_s; else 0.
    pub fn martensite_fraction(&self, temp_k: f64) -> f64 {
        if temp_k >= self.ms_temp {
            return 0.0;
        }
        let delta_t = self.ms_temp - temp_k;
        let f = self.max_martensite_fraction * (1.0 - (-self.km_coefficient * delta_t).exp());
        f.min(self.max_martensite_fraction)
    }

    /// Strength contribution from martensite (Hall-Petch-like).
    ///
    /// Δσ ≈ σ_α_prime * f_m  where σ_α_prime is the intrinsic martensite strength.
    pub fn strength_contribution(&self, temp_k: f64, martensite_strength: f64) -> f64 {
        self.martensite_fraction(temp_k) * martensite_strength
    }
}

// ---------------------------------------------------------------------------
// Binder Jetting Material
// ---------------------------------------------------------------------------

/// Material model for binder jetting (BJ) process.
#[derive(Debug, Clone)]
pub struct BinderJettingMaterial {
    /// Green body porosity after printing (0–1).
    pub green_porosity: f64,
    /// Target sintered porosity after sintering (0–1).
    pub sintered_porosity: f64,
    /// Linear shrinkage during sintering (0–1).
    pub linear_shrinkage: f64,
    /// Binder content (volume fraction).
    pub binder_fraction: f64,
    /// Sintering temperature (K).
    pub sintering_temp: f64,
    /// Sintering time (s).
    pub sintering_time: f64,
    /// Fully dense yield strength (Pa).
    pub dense_yield_strength: f64,
    /// Fully dense elastic modulus (Pa).
    pub dense_elastic_modulus: f64,
}

impl BinderJettingMaterial {
    /// Typical 316L stainless steel BJ parameters.
    pub fn steel_316l_bj() -> Self {
        Self {
            green_porosity: 0.40,
            sintered_porosity: 0.02,
            linear_shrinkage: 0.18,
            binder_fraction: 0.35,
            sintering_temp: 1593.0,
            sintering_time: 3600.0 * 6.0,
            dense_yield_strength: 530e6,
            dense_elastic_modulus: 193e9,
        }
    }

    /// Volumetric shrinkage from green to sintered state.
    pub fn volumetric_shrinkage(&self) -> f64 {
        1.0 - (1.0 - self.linear_shrinkage).powi(3)
    }

    /// Relative density of the sintered part.
    pub fn relative_density(&self) -> f64 {
        1.0 - self.sintered_porosity
    }

    /// Effective yield strength using power-law porosity correction.
    pub fn effective_yield_strength(&self) -> f64 {
        let p = self.sintered_porosity;
        self.dense_yield_strength * (1.0 - p).powf(2.0)
    }

    /// Effective elastic modulus using Ramakrishnan-Arunachalam model.
    pub fn effective_elastic_modulus(&self) -> f64 {
        let beta = 2.0;
        let p = self.sintered_porosity;
        self.dense_elastic_modulus * (1.0 - p).powi(2) / (1.0 + beta * p)
    }
}

// ---------------------------------------------------------------------------
// FDM Polymer Material
// ---------------------------------------------------------------------------

/// FDM/FFF polymer material with directional (anisotropic) properties.
#[derive(Debug, Clone)]
pub struct FdmPolymerMaterial {
    /// Polymer name/grade.
    pub name: String,
    /// Density (kg/m³).
    pub density: f64,
    /// In-layer (XY) tensile strength (Pa).
    pub tensile_strength_xy: f64,
    /// Through-layer (Z) tensile strength (Pa).
    pub tensile_strength_z: f64,
    /// In-layer Young's modulus (Pa).
    pub elastic_modulus_xy: f64,
    /// Through-layer Young's modulus (Pa).
    pub elastic_modulus_z: f64,
    /// In-layer elongation at break (0–1).
    pub elongation_xy: f64,
    /// Through-layer elongation at break (0–1).
    pub elongation_z: f64,
    /// Glass transition temperature (K).
    pub glass_transition_temp: f64,
    /// Layer adhesion strength (Pa) — inter-layer bond strength.
    pub layer_adhesion_strength: f64,
    /// Layer height (m).
    pub layer_height: f64,
    /// Raster width (m).
    pub raster_width: f64,
    /// Air gap between rasters (m, can be negative for overlap).
    pub air_gap: f64,
    /// Raster angle (degrees from X-axis).
    pub raster_angle_deg: f64,
}

impl FdmPolymerMaterial {
    /// Generic PLA-like material.
    pub fn pla_generic() -> Self {
        Self {
            name: "PLA".to_string(),
            density: 1_240.0,
            tensile_strength_xy: 60e6,
            tensile_strength_z: 35e6,
            elastic_modulus_xy: 3_500e6,
            elastic_modulus_z: 2_800e6,
            elongation_xy: 0.04,
            elongation_z: 0.02,
            glass_transition_temp: 333.0,
            layer_adhesion_strength: 30e6,
            layer_height: 0.2e-3,
            raster_width: 0.4e-3,
            air_gap: 0.0,
            raster_angle_deg: 45.0,
        }
    }

    /// Generic PEEK material (high-performance polymer).
    pub fn peek_generic() -> Self {
        Self {
            name: "PEEK".to_string(),
            density: 1_310.0,
            tensile_strength_xy: 100e6,
            tensile_strength_z: 60e6,
            elastic_modulus_xy: 4_000e6,
            elastic_modulus_z: 3_200e6,
            elongation_xy: 0.03,
            elongation_z: 0.02,
            glass_transition_temp: 416.0,
            layer_adhesion_strength: 50e6,
            layer_height: 0.15e-3,
            raster_width: 0.4e-3,
            air_gap: -0.05e-3,
            raster_angle_deg: 0.0,
        }
    }

    /// Anisotropy ratio: Z-strength / XY-strength.
    pub fn anisotropy_ratio(&self) -> f64 {
        if self.tensile_strength_xy < 1e-12 {
            1.0
        } else {
            self.tensile_strength_z / self.tensile_strength_xy
        }
    }

    /// Estimate effective tensile strength at a build angle `theta_deg`
    /// using a simple cosine interpolation.
    pub fn tensile_strength_at_angle(&self, theta_deg: f64) -> f64 {
        let theta = theta_deg.to_radians();
        let cos2 = theta.cos().powi(2);
        let sin2 = theta.sin().powi(2);
        cos2 * self.tensile_strength_xy + sin2 * self.tensile_strength_z
    }

    /// Estimated void volume fraction from air gap geometry.
    pub fn void_fraction(&self) -> f64 {
        if self.raster_width < 1e-12 || self.layer_height < 1e-12 {
            return 0.0;
        }
        let effective_gap = self.air_gap / self.raster_width;
        effective_gap.clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Support Structure Material
// ---------------------------------------------------------------------------

/// Material model for support structures in AM.
#[derive(Debug, Clone)]
pub struct SupportMaterial {
    /// Density scaling relative to bulk (0–1).
    pub density_fraction: f64,
    /// Elastic modulus fraction relative to bulk (0–1).
    pub modulus_fraction: f64,
    /// Critical force to remove support (N/m²).
    pub removal_force: f64,
    /// Support structure type.
    pub support_type: SupportType,
}

/// Type of AM support structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportType {
    /// Solid block support.
    Solid,
    /// Lattice / sparse support.
    Lattice,
    /// Cone-tree support.
    TreeLike,
    /// Water-soluble polymer support (FDM).
    Soluble,
    /// Powder-bed support (no structure needed in PBF).
    PowderBed,
}

impl SupportMaterial {
    /// Typical metal PBF block support.
    pub fn metal_block_support() -> Self {
        Self {
            density_fraction: 0.5,
            modulus_fraction: 0.4,
            removal_force: 1e6,
            support_type: SupportType::Solid,
        }
    }

    /// FDM soluble support (e.g. PVA).
    pub fn fdm_soluble() -> Self {
        Self {
            density_fraction: 0.8,
            modulus_fraction: 0.3,
            removal_force: 0.5e6,
            support_type: SupportType::Soluble,
        }
    }

    /// Estimate thermal resistance added by support block.
    ///
    /// R_thermal = h / (k_support * A)
    pub fn thermal_resistance(&self, height: f64, area: f64, bulk_conductivity: f64) -> f64 {
        let k = bulk_conductivity * self.modulus_fraction;
        if k < 1e-12 || area < 1e-20 {
            f64::INFINITY
        } else {
            height / (k * area)
        }
    }
}

// ---------------------------------------------------------------------------
// Process-Structure-Property (PSP) Linkage
// ---------------------------------------------------------------------------

/// PSP linkage: combines process parameters → structure descriptors → properties.
#[derive(Debug, Clone)]
pub struct PspLinkage {
    /// Material system.
    pub material: PbfMaterial,
    /// Process parameters.
    pub process: PbfProcessParams,
    /// Estimated porosity (0–1) from process window.
    pub estimated_porosity: f64,
}

impl PspLinkage {
    /// Construct PSP from material and process parameters.
    pub fn new(material: PbfMaterial, process: PbfProcessParams) -> Self {
        let porosity = estimate_porosity_from_energy_density(
            process.volumetric_energy_density(),
            material.density,
        );
        Self {
            material,
            process,
            estimated_porosity: porosity,
        }
    }

    /// Effective yield strength incorporating porosity.
    pub fn effective_yield_strength(&self) -> f64 {
        let p = self.estimated_porosity;
        self.material.yield_strength * (1.0 - p).powf(2.0)
    }

    /// Effective elastic modulus incorporating porosity.
    pub fn effective_elastic_modulus(&self) -> f64 {
        let beta = 2.0;
        let p = self.estimated_porosity;
        self.material.elastic_modulus * (1.0 - p).powi(2) / (1.0 + beta * p)
    }

    /// Thermal gradient estimate.
    pub fn thermal_gradient(&self, ambient_temp: f64) -> ThermalGradient {
        ThermalGradient::from_rosenthal(&self.material, &self.process, ambient_temp)
    }

    /// Melt pool geometry estimate.
    pub fn melt_pool(&self, ambient_temp: f64) -> MeltPoolGeometry {
        MeltPoolGeometry::eagar_tsai(&self.material, &self.process, ambient_temp)
    }

    /// Estimated relative density (1 − porosity).
    pub fn relative_density(&self) -> f64 {
        1.0 - self.estimated_porosity
    }
}

/// Estimate porosity from volumetric energy density using a sigmoid curve.
///
/// At optimal energy density (~60–80 J/mm³) porosity is minimal (~0.1%).
/// Too low → lack of fusion; too high → keyhole porosity.
pub fn estimate_porosity_from_energy_density(energy_density_j_per_m3: f64, _density: f64) -> f64 {
    // Convert to J/mm³
    let ev = energy_density_j_per_m3 * 1e-9;
    // Optimal energy density for most metals ≈ 70 J/mm³
    let ev_opt = 70.0_f64;
    let delta = (ev - ev_opt).abs() / ev_opt;
    // Porosity rises parabolically from ~0.1% at optimum
    let base = 0.001_f64;
    (base + 0.5 * delta * delta).min(0.30)
}

// ---------------------------------------------------------------------------
// Scan Strategy Effects
// ---------------------------------------------------------------------------

/// Effect of scan strategy on texture and residual stress.
#[derive(Debug, Clone)]
pub struct ScanStrategyEffect {
    /// Rotation angle between consecutive layers (degrees).
    pub rotation_angle_deg: f64,
    /// Estimated texture coefficient (0 = random, 1 = fully textured).
    pub texture_coefficient: f64,
    /// Residual stress scaling factor relative to unidirectional scanning.
    pub residual_stress_factor: f64,
    /// Estimated relative density achievable with this strategy.
    pub relative_density: f64,
}

impl ScanStrategyEffect {
    /// Effects for unidirectional (0° rotation) scanning.
    pub fn unidirectional() -> Self {
        Self {
            rotation_angle_deg: 0.0,
            texture_coefficient: 0.85,
            residual_stress_factor: 1.0,
            relative_density: 0.993,
        }
    }

    /// Effects for 67° rotation (common in SLM).
    pub fn rotating_67() -> Self {
        Self {
            rotation_angle_deg: 67.0,
            texture_coefficient: 0.35,
            residual_stress_factor: 0.65,
            relative_density: 0.997,
        }
    }

    /// Effects for 90° alternating scanning.
    pub fn alternating_90() -> Self {
        Self {
            rotation_angle_deg: 90.0,
            texture_coefficient: 0.50,
            residual_stress_factor: 0.75,
            relative_density: 0.995,
        }
    }

    /// Effects for island/checkerboard scanning.
    pub fn island() -> Self {
        Self {
            rotation_angle_deg: 90.0,
            texture_coefficient: 0.40,
            residual_stress_factor: 0.60,
            relative_density: 0.996,
        }
    }

    /// Estimate the number of unique scan orientations in N layers.
    pub fn unique_orientations_in_n_layers(&self, n_layers: usize) -> usize {
        if self.rotation_angle_deg < 1e-6 {
            return 1;
        }
        let period = (360.0 / self.rotation_angle_deg).round() as usize;
        period.min(n_layers)
    }
}

// ---------------------------------------------------------------------------
// Surface Roughness from AM
// ---------------------------------------------------------------------------

/// AM surface roughness model.
///
/// Covers staircase effect, spattering, and powder adhesion contributions.
#[derive(Debug, Clone)]
pub struct AmSurfaceRoughness {
    /// Layer thickness (m).
    pub layer_thickness: f64,
    /// Powder particle diameter D50 (m).
    pub powder_d50: f64,
    /// Build angle relative to horizontal (degrees, 0 = flat top, 90 = vertical wall).
    pub build_angle_deg: f64,
    /// Laser spot radius (m).
    pub spot_radius: f64,
}

impl AmSurfaceRoughness {
    /// Staircase roughness Ra from layer thickness and build angle.
    ///
    /// Ra_staircase ≈ (t / 4) |cos θ| / sin θ   for θ in (0°, 90°)
    pub fn ra_staircase(&self) -> f64 {
        let theta = self.build_angle_deg.to_radians();
        let sin_t = theta.sin().max(1e-6);
        let cos_t = theta.cos().abs();
        self.layer_thickness / 4.0 * cos_t / sin_t
    }

    /// Contribution from partially sintered/attached powder particles.
    ///
    /// Ra_powder ≈ D50 / 4
    pub fn ra_powder_adhesion(&self) -> f64 {
        self.powder_d50 / 4.0
    }

    /// Total estimated surface roughness Ra (m).
    pub fn ra_total(&self) -> f64 {
        (self.ra_staircase().powi(2) + self.ra_powder_adhesion().powi(2)).sqrt()
    }

    /// Roughness in micrometres.
    pub fn ra_total_um(&self) -> f64 {
        self.ra_total() * 1e6
    }

    /// Estimate fatigue life reduction factor due to surface roughness.
    ///
    /// Uses modified Goodman: Kf ≈ 1 + q * (Kt - 1) where q = notch sensitivity.
    pub fn fatigue_reduction_factor(&self) -> f64 {
        let ra_um = self.ra_total_um();
        // Empirical: Kf ≈ 1 + 0.06 * Ra_um for metals
        1.0 + 0.06 * ra_um
    }
}

// ---------------------------------------------------------------------------
// IN718 Specific Model
// ---------------------------------------------------------------------------

/// IN718 nickel superalloy properties for SLM/EBM.
#[derive(Debug, Clone)]
pub struct In718Properties {
    /// Gamma-prime (γ') volume fraction.
    pub gamma_prime_fraction: f64,
    /// Gamma-double-prime (γ'') volume fraction.
    pub gamma_double_prime_fraction: f64,
    /// Delta phase volume fraction.
    pub delta_fraction: f64,
    /// Creep coefficient (power-law).
    pub creep_coefficient: f64,
    /// Creep exponent n.
    pub creep_exponent: f64,
    /// Creep activation energy (J/mol).
    pub creep_activation_energy: f64,
}

impl In718Properties {
    /// Typical as-built SLM IN718 properties.
    pub fn as_built_slm() -> Self {
        Self {
            gamma_prime_fraction: 0.03,
            gamma_double_prime_fraction: 0.12,
            delta_fraction: 0.01,
            creep_coefficient: 2.5e-24,
            creep_exponent: 5.2,
            creep_activation_energy: 285_000.0,
        }
    }

    /// Typical heat-treated IN718 properties (solution annealed + aged).
    pub fn heat_treated() -> Self {
        Self {
            gamma_prime_fraction: 0.05,
            gamma_double_prime_fraction: 0.18,
            delta_fraction: 0.005,
            creep_coefficient: 1.5e-24,
            creep_exponent: 5.2,
            creep_activation_energy: 290_000.0,
        }
    }

    /// Estimate yield strength contribution from precipitates (Hall-Petch-like).
    ///
    /// Δσ ≈ M G b √ρ  simplified to ≈ C * (f_prime + f_doubleprime)^0.5
    pub fn precipitation_strengthening(&self) -> f64 {
        let c = 2000e6; // empirical constant
        c * (self.gamma_prime_fraction + self.gamma_double_prime_fraction).sqrt()
    }

    /// Minimum creep rate (s⁻¹) at stress σ and temperature T.
    ///
    /// Norton power law: ε̇ = A σⁿ exp(-Q / RT).
    pub fn creep_rate(&self, stress_pa: f64, temp_k: f64) -> f64 {
        self.creep_coefficient
            * stress_pa.powf(self.creep_exponent)
            * (-self.creep_activation_energy / (R_GAS * temp_k)).exp()
    }
}

// ---------------------------------------------------------------------------
// Energy Density Parameter Space
// ---------------------------------------------------------------------------

/// Process window analysis based on volumetric energy density.
#[derive(Debug, Clone)]
pub struct ProcessWindow {
    /// Minimum energy density for full melting (J/m³).
    pub ev_min: f64,
    /// Maximum energy density before keyhole porosity (J/m³).
    pub ev_max: f64,
    /// Optimal energy density for maximum density (J/m³).
    pub ev_optimal: f64,
}

impl ProcessWindow {
    /// Default process window for Ti-6Al-4V SLM (J/mm³ converted to J/m³).
    pub fn ti6al4v_slm() -> Self {
        Self {
            ev_min: 50e9,
            ev_max: 130e9,
            ev_optimal: 75e9,
        }
    }

    /// Default process window for 316L SLM.
    pub fn steel_316l_slm() -> Self {
        Self {
            ev_min: 45e9,
            ev_max: 120e9,
            ev_optimal: 70e9,
        }
    }

    /// Check whether given energy density is within the process window.
    pub fn is_in_window(&self, ev: f64) -> bool {
        ev >= self.ev_min && ev <= self.ev_max
    }

    /// Estimated relative density at energy density `ev`.
    pub fn estimated_relative_density(&self, ev: f64) -> f64 {
        if ev < self.ev_min {
            // Lack-of-fusion regime
            let frac = ev / self.ev_min;
            0.80 + 0.18 * frac
        } else if ev > self.ev_max {
            // Keyhole regime
            let excess = (ev - self.ev_max) / self.ev_max;
            0.999 - 0.15 * excess * excess
        } else {
            // In-window
            let x = (ev - self.ev_optimal).abs() / (self.ev_max - self.ev_min);
            0.999 - 0.05 * x * x
        }
        .clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Hardness Prediction
// ---------------------------------------------------------------------------

/// Vickers hardness (HV) from relative density and microstructure.
///
/// Simple linear model: HV = HV_dense * (1 − k_p * porosity)
pub fn vickers_hardness(hv_dense: f64, porosity: f64, k_p: f64) -> f64 {
    hv_dense * (1.0 - k_p * porosity.clamp(0.0, 1.0))
}

// ---------------------------------------------------------------------------
// Thermal Stress from Build Simulation
// ---------------------------------------------------------------------------

/// Layer-by-layer thermal stress accumulation.
///
/// Each deposited layer imposes a thermal strain `α ΔT` on the underlying body.
/// The cumulative residual stress is computed via a simple analytical laminate model.
pub fn cumulative_residual_stress(
    layer_temps: &[f64],
    ambient_temp: f64,
    thermal_expansion: f64,
    elastic_modulus: f64,
    poisson_ratio: f64,
) -> Vec<f64> {
    layer_temps
        .iter()
        .map(|&t| {
            let delta_t = (t - ambient_temp).abs();

            elastic_modulus * thermal_expansion * delta_t / (1.0 - poisson_ratio)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PbfMaterial ---

    #[test]
    fn test_ti6al4v_thermal_diffusivity_positive() {
        let m = PbfMaterial::ti6al4v_slm();
        assert!(m.thermal_diffusivity() > 0.0);
    }

    #[test]
    fn test_ti6al4v_melting_range_positive() {
        let m = PbfMaterial::ti6al4v_slm();
        assert!(m.melting_range() > 0.0);
    }

    #[test]
    fn test_steel_density_reasonable() {
        let m = PbfMaterial::steel_316l_slm();
        assert!(m.density > 7_000.0 && m.density < 9_000.0);
    }

    #[test]
    fn test_alsi10mg_absorptivity_range() {
        let m = PbfMaterial::alsi10mg_slm();
        assert!(m.absorptivity > 0.0 && m.absorptivity <= 1.0);
    }

    // --- PbfProcessParams ---

    #[test]
    fn test_volumetric_energy_density_positive() {
        let p = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.8,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let ev = p.volumetric_energy_density();
        assert!(ev > 0.0 && ev.is_finite());
    }

    #[test]
    fn test_linear_energy_density_consistent() {
        let p = PbfProcessParams {
            power: 200.0,
            scan_speed: 1.0,
            hatch_spacing: 100e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        assert!((p.linear_energy_density() - 200.0).abs() < 1e-6);
    }

    #[test]
    fn test_interaction_time_positive() {
        let p = PbfProcessParams {
            power: 200.0,
            scan_speed: 1.0,
            hatch_spacing: 100e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        assert!(p.interaction_time() > 0.0);
    }

    // --- MeltPoolGeometry ---

    #[test]
    fn test_melt_pool_volume_positive() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.5,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let mp = MeltPoolGeometry::eagar_tsai(&mat, &params, 300.0);
        assert!(mp.volume() > 0.0);
    }

    #[test]
    fn test_melt_pool_aspect_ratio_ge_1() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.5,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let mp = MeltPoolGeometry::eagar_tsai(&mat, &params, 300.0);
        assert!(mp.aspect_ratio() >= 1.0);
    }

    // --- ThermalGradient ---

    #[test]
    fn test_thermal_gradient_positive() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.5,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let tg = ThermalGradient::from_rosenthal(&mat, &params, 300.0);
        assert!(tg.gradient_magnitude > 0.0);
        assert!(tg.cooling_rate > 0.0);
    }

    #[test]
    fn test_g_times_r_positive() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.5,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let tg = ThermalGradient::from_rosenthal(&mat, &params, 300.0);
        assert!(tg.g_times_r() > 0.0);
    }

    // --- ResidualStressModel ---

    #[test]
    fn test_residual_stress_nonnegative() {
        let mat = PbfMaterial::ti6al4v_slm();
        let rs = ResidualStressModel::from_material(&mat);
        let sigma = rs.peak_residual_stress(500.0);
        assert!(sigma >= 0.0);
    }

    #[test]
    fn test_residual_stress_capped_at_yield() {
        let mat = PbfMaterial::ti6al4v_slm();
        let rs = ResidualStressModel::from_material(&mat);
        let sigma = rs.peak_residual_stress(10_000.0); // very large delta T
        assert!(sigma <= rs.yield_strength_hot + 1.0);
    }

    #[test]
    fn test_stoney_curvature_increases_with_stress() {
        let mat = PbfMaterial::ti6al4v_slm();
        let rs = ResidualStressModel::from_material(&mat);
        let k1 = rs.stoney_curvature(100e6, 30e-6, 10e-3);
        let k2 = rs.stoney_curvature(200e6, 30e-6, 10e-3);
        assert!(k2 > k1);
    }

    // --- PorosityModel ---

    #[test]
    fn test_porosity_zero_gives_full_properties() {
        let pm = PorosityModel {
            e0: 114e9,
            yield_strength0: 930e6,
            tensile_strength0: 1000e6,
            density0: 4430.0,
        };
        assert!((pm.elastic_modulus(0.0) - pm.e0).abs() < 1e-3 * pm.e0);
        assert!((pm.yield_strength(0.0) - pm.yield_strength0).abs() < 1e-6);
    }

    #[test]
    fn test_porosity_modulus_decreasing() {
        let pm = PorosityModel {
            e0: 114e9,
            yield_strength0: 930e6,
            tensile_strength0: 1000e6,
            density0: 4430.0,
        };
        let e1 = pm.elastic_modulus(0.01);
        let e2 = pm.elastic_modulus(0.05);
        assert!(e1 > e2);
    }

    #[test]
    fn test_fatigue_reduction_ge_1() {
        let pm = PorosityModel {
            e0: 114e9,
            yield_strength0: 930e6,
            tensile_strength0: 1000e6,
            density0: 4430.0,
        };
        assert!(pm.fatigue_strength_reduction(0.02) >= 1.0);
    }

    // --- GrainGrowthModel ---

    #[test]
    fn test_grain_growth_increases_with_time() {
        let gg = GrainGrowthModel::ti6al4v_beta();
        let d1 = gg.grain_diameter_isothermal(1200.0, 60.0);
        let d2 = gg.grain_diameter_isothermal(1200.0, 3600.0);
        assert!(d2 > d1);
    }

    #[test]
    fn test_grain_growth_initial_diameter_nonnegative() {
        let gg = GrainGrowthModel::ti6al4v_beta();
        let d = gg.grain_diameter_isothermal(300.0, 0.0);
        assert!(d >= 0.0);
    }

    // --- MartensiticModel ---

    #[test]
    fn test_martensite_below_ms() {
        let m = MartensiticModel::ti6al4v_martensite();
        let f = m.martensite_fraction(600.0);
        assert!(f > 0.0 && f <= 1.0);
    }

    #[test]
    fn test_martensite_above_ms_is_zero() {
        let m = MartensiticModel::ti6al4v_martensite();
        let f = m.martensite_fraction(1000.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_martensite_increases_cooling() {
        let m = MartensiticModel::ti6al4v_martensite();
        let f1 = m.martensite_fraction(800.0);
        let f2 = m.martensite_fraction(600.0);
        assert!(f2 > f1);
    }

    // --- BinderJettingMaterial ---

    #[test]
    fn test_bj_volumetric_shrinkage_positive() {
        let bj = BinderJettingMaterial::steel_316l_bj();
        assert!(bj.volumetric_shrinkage() > 0.0);
    }

    #[test]
    fn test_bj_relative_density_high() {
        let bj = BinderJettingMaterial::steel_316l_bj();
        assert!(bj.relative_density() > 0.95);
    }

    #[test]
    fn test_bj_effective_yield_strength_positive() {
        let bj = BinderJettingMaterial::steel_316l_bj();
        assert!(bj.effective_yield_strength() > 0.0);
    }

    // --- FdmPolymerMaterial ---

    #[test]
    fn test_fdm_anisotropy_ratio_range() {
        let fdm = FdmPolymerMaterial::pla_generic();
        let r = fdm.anisotropy_ratio();
        assert!(r > 0.0 && r <= 1.0);
    }

    #[test]
    fn test_fdm_tensile_at_0_equals_xy() {
        let fdm = FdmPolymerMaterial::pla_generic();
        let sigma = fdm.tensile_strength_at_angle(0.0);
        assert!((sigma - fdm.tensile_strength_xy).abs() < 1.0);
    }

    #[test]
    fn test_fdm_tensile_at_90_equals_z() {
        let fdm = FdmPolymerMaterial::pla_generic();
        let sigma = fdm.tensile_strength_at_angle(90.0);
        assert!((sigma - fdm.tensile_strength_z).abs() < 1.0);
    }

    // --- ScanStrategyEffect ---

    #[test]
    fn test_scan_67_lower_texture_than_unidirectional() {
        let uni = ScanStrategyEffect::unidirectional();
        let rot = ScanStrategyEffect::rotating_67();
        assert!(rot.texture_coefficient < uni.texture_coefficient);
    }

    #[test]
    fn test_scan_67_lower_stress_than_unidirectional() {
        let uni = ScanStrategyEffect::unidirectional();
        let rot = ScanStrategyEffect::rotating_67();
        assert!(rot.residual_stress_factor < uni.residual_stress_factor);
    }

    #[test]
    fn test_unique_orientations_unidirectional() {
        let uni = ScanStrategyEffect::unidirectional();
        assert_eq!(uni.unique_orientations_in_n_layers(100), 1);
    }

    // --- AmSurfaceRoughness ---

    #[test]
    fn test_surface_roughness_ra_positive() {
        let r = AmSurfaceRoughness {
            layer_thickness: 30e-6,
            powder_d50: 30e-6,
            build_angle_deg: 45.0,
            spot_radius: 35e-6,
        };
        assert!(r.ra_total() > 0.0);
    }

    #[test]
    fn test_surface_roughness_increases_with_layer_thickness() {
        let r1 = AmSurfaceRoughness {
            layer_thickness: 30e-6,
            powder_d50: 30e-6,
            build_angle_deg: 45.0,
            spot_radius: 35e-6,
        };
        let r2 = AmSurfaceRoughness {
            layer_thickness: 60e-6,
            powder_d50: 30e-6,
            build_angle_deg: 45.0,
            spot_radius: 35e-6,
        };
        assert!(r2.ra_total() > r1.ra_total());
    }

    #[test]
    fn test_surface_fatigue_reduction_ge_1() {
        let r = AmSurfaceRoughness {
            layer_thickness: 30e-6,
            powder_d50: 30e-6,
            build_angle_deg: 45.0,
            spot_radius: 35e-6,
        };
        assert!(r.fatigue_reduction_factor() >= 1.0);
    }

    // --- PSP linkage ---

    #[test]
    fn test_psp_relative_density_in_range() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.7,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let psp = PspLinkage::new(mat, params);
        let rd = psp.relative_density();
        assert!(rd > 0.0 && rd <= 1.0);
    }

    #[test]
    fn test_psp_effective_yield_positive() {
        let mat = PbfMaterial::ti6al4v_slm();
        let params = PbfProcessParams {
            power: 200.0,
            scan_speed: 0.7,
            hatch_spacing: 110e-6,
            layer_thickness: 30e-6,
            spot_radius: 35e-6,
            preheat_temp: 300.0,
            rotation_angle_deg: 67.0,
        };
        let psp = PspLinkage::new(mat, params);
        assert!(psp.effective_yield_strength() > 0.0);
    }

    // --- ProcessWindow ---

    #[test]
    fn test_process_window_in_range() {
        let pw = ProcessWindow::ti6al4v_slm();
        assert!(pw.is_in_window(pw.ev_optimal));
    }

    #[test]
    fn test_process_window_out_of_range_low() {
        let pw = ProcessWindow::ti6al4v_slm();
        assert!(!pw.is_in_window(pw.ev_min * 0.5));
    }

    #[test]
    fn test_relative_density_at_optimal_near_1() {
        let pw = ProcessWindow::ti6al4v_slm();
        let rd = pw.estimated_relative_density(pw.ev_optimal);
        assert!(rd > 0.99);
    }

    // --- In718 ---

    #[test]
    fn test_in718_precipitation_strengthening_positive() {
        let props = In718Properties::as_built_slm();
        assert!(props.precipitation_strengthening() > 0.0);
    }

    #[test]
    fn test_in718_creep_rate_positive() {
        let props = In718Properties::heat_treated();
        let cr = props.creep_rate(500e6, 923.0);
        assert!(cr > 0.0);
    }

    // --- Miscellaneous ---

    #[test]
    fn test_vickers_hardness_decreases_with_porosity() {
        let hv1 = vickers_hardness(350.0, 0.0, 3.0);
        let hv2 = vickers_hardness(350.0, 0.05, 3.0);
        assert!(hv1 > hv2);
    }

    #[test]
    fn test_cumulative_residual_stress_count() {
        let temps = vec![1000.0, 900.0, 800.0];
        let stresses = cumulative_residual_stress(&temps, 300.0, 8.6e-6, 114e9, 0.342);
        assert_eq!(stresses.len(), 3);
    }

    #[test]
    fn test_estimate_porosity_in_window_is_low() {
        let p = estimate_porosity_from_energy_density(70e9, 4430.0);
        assert!(p < 0.05);
    }

    #[test]
    fn test_support_thermal_resistance_finite() {
        let s = SupportMaterial::metal_block_support();
        let r = s.thermal_resistance(5e-3, 1e-4, 20.0);
        assert!(r.is_finite() && r > 0.0);
    }
}
