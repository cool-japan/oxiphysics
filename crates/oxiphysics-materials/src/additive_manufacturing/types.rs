//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Geometric parameters for a single AM layer.
#[derive(Debug, Clone)]
pub struct LayerModel {
    /// Layer thickness (m).
    pub layer_thickness: f64,
    /// Hatch spacing between adjacent scan tracks (m).
    pub hatch_spacing: f64,
    /// Nominal laser spot radius (m).
    pub spot_radius: f64,
    /// Scan speed (m/s).
    pub scan_speed: f64,
    /// Laser power (W).
    pub laser_power: f64,
    /// Scanning strategy applied to this layer.
    pub scan_strategy: ScanStrategy,
}
impl LayerModel {
    /// Construct a new [`LayerModel`].
    pub fn new(
        layer_thickness: f64,
        hatch_spacing: f64,
        spot_radius: f64,
        scan_speed: f64,
        laser_power: f64,
        scan_strategy: ScanStrategy,
    ) -> Self {
        Self {
            layer_thickness,
            hatch_spacing,
            spot_radius,
            scan_speed,
            laser_power,
            scan_strategy,
        }
    }
    /// Volumetric energy density (J/m³), often written as *E_v*.
    ///
    /// ```text
    /// E_v = P / (v · h · t)
    /// ```
    /// where *P* = power, *v* = scan speed, *h* = hatch spacing, *t* = layer thickness.
    pub fn volumetric_energy_density(&self) -> f64 {
        self.laser_power / (self.scan_speed * self.hatch_spacing * self.layer_thickness)
    }
    /// Linear energy density (J/m) — power divided by scan speed.
    pub fn linear_energy_density(&self) -> f64 {
        self.laser_power / self.scan_speed
    }
    /// Hatch overlap fraction (0–1).  Positive means tracks overlap.
    ///
    /// `overlap = 1 − hatch_spacing / (2 · spot_radius)`
    pub fn hatch_overlap(&self) -> f64 {
        1.0_f64 - self.hatch_spacing / (2.0_f64 * self.spot_radius)
    }
    /// Approximate number of scan tracks required to cover a square domain of
    /// side `width` (m).
    pub fn num_tracks(&self, width: f64) -> usize {
        ((width / self.hatch_spacing).ceil() as usize).max(1)
    }
}
/// Defect formation models for LPBF.
#[derive(Debug, Clone)]
pub struct DefectModels {
    /// Measured relative density (0–1).
    pub relative_density: f64,
    /// Volumetric energy density (J/m³).
    pub energy_density: f64,
    /// Threshold energy density below which lack-of-fusion pores form (J/m³).
    pub lof_threshold: f64,
    /// Threshold energy density above which keyhole pores form (J/m³).
    pub keyhole_threshold: f64,
    /// Ultimate tensile strength of the matrix material (Pa).
    pub uts: f64,
    /// Fracture toughness of the matrix (Pa·√m).
    pub fracture_toughness: f64,
    /// Characteristic defect size (m) — used for cracking criterion.
    pub defect_size: f64,
}
impl DefectModels {
    /// Construct a [`DefectModels`] struct.
    pub fn new(
        relative_density: f64,
        energy_density: f64,
        lof_threshold: f64,
        keyhole_threshold: f64,
        uts: f64,
        fracture_toughness: f64,
        defect_size: f64,
    ) -> Self {
        Self {
            relative_density,
            energy_density,
            lof_threshold,
            keyhole_threshold,
            uts,
            fracture_toughness,
            defect_size,
        }
    }
    /// True if the energy density is below the lack-of-fusion threshold.
    pub fn has_lack_of_fusion(&self) -> bool {
        self.energy_density < self.lof_threshold
    }
    /// True if the energy density is above the keyhole threshold.
    pub fn has_keyhole_porosity(&self) -> bool {
        self.energy_density > self.keyhole_threshold
    }
    /// Porosity fraction estimated from the relative density.
    pub fn porosity_fraction(&self) -> f64 {
        (1.0_f64 - self.relative_density).max(0.0_f64)
    }
    /// Cracking criterion using a stress-intensity-factor approach.
    ///
    /// Returns `true` when the mode-I SIF exceeds the fracture toughness:
    /// `K_I = σ_UTS · √(π · a) > K_IC`.
    pub fn cracking_criterion(&self) -> bool {
        let k_i = self.uts * (PI * self.defect_size).sqrt();
        k_i >= self.fracture_toughness
    }
    /// Quality flag: process is in the "sweet spot" (no LOF, no keyhole).
    pub fn in_process_window(&self) -> bool {
        !self.has_lack_of_fusion() && !self.has_keyhole_porosity()
    }
}
/// Overhang detection and support-structure estimation for AM parts.
///
/// Any surface facet inclined more than `max_overhang_angle` from the
/// horizontal is considered self-supporting; steeper facets require support.
#[derive(Debug, Clone)]
pub struct SupportStructure {
    /// Maximum self-supporting overhang angle from horizontal (radians).
    /// Typical LPBF value is ~45° (π/4).
    pub max_overhang_angle: f64,
    /// Material density used to compute support mass (kg/m³).
    pub material_density: f64,
    /// Support fill fraction relative to solid (0–1).
    /// Lattice supports typically use 0.05–0.20.
    pub support_fill_fraction: f64,
    /// Estimated support volume (m³) — populated by `estimate_support`.
    pub support_volume_m3: f64,
}
impl SupportStructure {
    /// Construct a [`SupportStructure`] model.
    ///
    /// # Arguments
    /// * `max_overhang_angle_deg` — maximum self-supporting angle in **degrees**.
    /// * `material_density` — kg/m³.
    /// * `support_fill_fraction` — lattice density relative to solid (0–1).
    pub fn new(
        max_overhang_angle_deg: f64,
        material_density: f64,
        support_fill_fraction: f64,
    ) -> Self {
        Self {
            max_overhang_angle: max_overhang_angle_deg.to_radians(),
            material_density,
            support_fill_fraction: support_fill_fraction.clamp(0.0_f64, 1.0_f64),
            support_volume_m3: 0.0_f64,
        }
    }
    /// Determine whether a facet with outward-normal vector `normal` (unit vector,
    /// z component positive upward) requires support.
    ///
    /// A facet requires support when it is downward-facing and the angle
    /// between its outward normal and the downward vertical is less than
    /// `max_overhang_angle` (i.e. the facet tilts more toward horizontal
    /// than the self-supporting limit).
    pub fn requires_support(&self, normal: [f64; 3]) -> bool {
        let nz = normal[2];
        if nz >= 0.0_f64 {
            return false;
        }
        let cos_a = (-nz).clamp(-1.0_f64, 1.0_f64);
        let angle_from_down = cos_a.acos();
        angle_from_down < self.max_overhang_angle
    }
    /// Estimate support volume from a bounding-box projection.
    ///
    /// Simplified model: support fills the volume between the part's lowest
    /// unsupported surface and the build plate, scaled by `support_fill_fraction`.
    ///
    /// # Arguments
    /// * `part_volume_m3` — volume of the part itself (m³).
    /// * `unsupported_fraction` — fraction of the part's projected area that
    ///   requires support (0–1), estimated from a geometry scan.
    /// * `part_height_m` — build height of the part (m).
    pub fn estimate_support(
        &mut self,
        part_volume_m3: f64,
        unsupported_fraction: f64,
        part_height_m: f64,
    ) {
        let projected_area = part_volume_m3 / part_height_m.max(1.0e-9_f64);
        let support_area = projected_area * unsupported_fraction.clamp(0.0_f64, 1.0_f64);
        self.support_volume_m3 =
            support_area * (0.5_f64 * part_height_m) * self.support_fill_fraction;
    }
    /// Estimated support mass (kg).
    pub fn support_mass_kg(&self) -> f64 {
        self.support_volume_m3 * self.material_density
    }
    /// Material waste fraction: support mass relative to (part + support) mass.
    ///
    /// # Arguments
    /// * `part_mass_kg` — mass of the finished part (kg).
    pub fn waste_fraction(&self, part_mass_kg: f64) -> f64 {
        let total = part_mass_kg + self.support_mass_kg();
        if total < 1.0e-15_f64 {
            return 0.0_f64;
        }
        self.support_mass_kg() / total
    }
}
/// Characteristic dimensions of the melt pool (all in metres).
#[derive(Debug, Clone, Copy)]
pub struct MeltPoolGeometry {
    /// Half-length of the melt pool along the scan direction (m).
    pub half_length: f64,
    /// Half-width of the melt pool (m).
    pub half_width: f64,
    /// Depth of the melt pool (m).
    pub depth: f64,
}
impl MeltPoolGeometry {
    /// Aspect ratio: depth / width.
    pub fn aspect_ratio(&self) -> f64 {
        self.depth / (2.0_f64 * self.half_width)
    }
    /// Volume of the melt pool modelled as a half-ellipsoid (m³).
    pub fn volume(&self) -> f64 {
        (2.0_f64 / 3.0_f64) * PI * self.half_length * self.half_width * self.depth
    }
}
/// Simplified residual stress build-up model for LPBF.
///
/// Uses the temperature-gradient mechanism (TGM) and cool-down shrinkage
/// to estimate in-plane longitudinal residual stress.
#[derive(Debug, Clone)]
pub struct ResidualStress {
    /// Young's modulus of the build material (Pa).
    pub youngs_modulus: f64,
    /// Coefficient of thermal expansion (1/K).
    pub cte: f64,
    /// Yield strength at room temperature (Pa).
    pub yield_strength: f64,
    /// Peak temperature during processing (K).
    pub peak_temperature: f64,
    /// Ambient (room) temperature (K).
    pub ambient_temp: f64,
}
impl ResidualStress {
    /// Construct a [`ResidualStress`] model.
    pub fn new(
        youngs_modulus: f64,
        cte: f64,
        yield_strength: f64,
        peak_temperature: f64,
        ambient_temp: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            cte,
            yield_strength,
            peak_temperature,
            ambient_temp,
        }
    }
    /// Estimated longitudinal residual stress (Pa) from the TGM.
    ///
    /// σ_res = min(E · α · ΔT, σ_yield) — capped at the yield strength.
    pub fn longitudinal_stress(&self) -> f64 {
        let delta_t = self.peak_temperature - self.ambient_temp;
        let elastic = self.youngs_modulus * self.cte * delta_t;
        elastic.min(self.yield_strength)
    }
    /// Stress ratio σ_res / σ_yield (dimensionless).
    pub fn stress_ratio(&self) -> f64 {
        self.longitudinal_stress() / self.yield_strength
    }
    /// Layer-by-layer accumulated stress assuming `n_layers` deposited.
    ///
    /// Uses a simplified accumulation model where each layer adds a fraction
    /// `layer_fraction` (0–1) of the single-layer residual stress.
    pub fn accumulated_stress(&self, n_layers: usize, layer_fraction: f64) -> f64 {
        let sigma_per_layer = self.longitudinal_stress() * layer_fraction;
        let total = sigma_per_layer * n_layers as f64;
        total.min(self.yield_strength)
    }
    /// Boolean flag: is the residual stress likely to cause delamination?
    ///
    /// Uses a simplified criterion: stress > 0.8 × σ_yield.
    pub fn delamination_risk(&self) -> bool {
        self.longitudinal_stress() > 0.8_f64 * self.yield_strength
    }
}
/// Functionally-graded material (FGM) deposition in multi-material AM.
///
/// Describes a composition transition between two constituent materials
/// along the build direction using a power-law profile.
#[derive(Debug, Clone)]
pub struct MultiMaterialAM {
    /// Young's modulus of material A (Pa).
    pub e_a: f64,
    /// Young's modulus of material B (Pa).
    pub e_b: f64,
    /// Density of material A (kg/m³).
    pub rho_a: f64,
    /// Density of material B (kg/m³).
    pub rho_b: f64,
    /// Thermal conductivity of material A (W/m·K).
    pub k_a: f64,
    /// Thermal conductivity of material B (W/m·K).
    pub k_b: f64,
    /// Power-law exponent *n* for the grading profile.
    /// n=0 → pure A; n→∞ → pure B; n=1 → linear.
    pub grading_exponent: f64,
    /// Total build height over which grading occurs (m).
    pub grading_height_m: f64,
}
impl MultiMaterialAM {
    /// Construct a [`MultiMaterialAM`] graded transition.
    pub fn new(
        e_a: f64,
        e_b: f64,
        rho_a: f64,
        rho_b: f64,
        k_a: f64,
        k_b: f64,
        grading_exponent: f64,
        grading_height_m: f64,
    ) -> Self {
        Self {
            e_a,
            e_b,
            rho_a,
            rho_b,
            k_a,
            k_b,
            grading_exponent,
            grading_height_m,
        }
    }
    /// Volume fraction of material B at height `z_m` (m from the bottom).
    ///
    /// `V_B(z) = (z / H)^n`
    pub fn volume_fraction_b(&self, z_m: f64) -> f64 {
        let xi = (z_m / self.grading_height_m.max(1.0e-15_f64)).clamp(0.0_f64, 1.0_f64);
        xi.powf(self.grading_exponent)
    }
    /// Effective Young's modulus at height `z_m` using the rule of mixtures.
    pub fn effective_modulus(&self, z_m: f64) -> f64 {
        let vb = self.volume_fraction_b(z_m);
        self.e_a * (1.0_f64 - vb) + self.e_b * vb
    }
    /// Effective density at height `z_m` (kg/m³).
    pub fn effective_density(&self, z_m: f64) -> f64 {
        let vb = self.volume_fraction_b(z_m);
        self.rho_a * (1.0_f64 - vb) + self.rho_b * vb
    }
    /// Effective thermal conductivity at height `z_m` (W/m·K).
    pub fn effective_conductivity(&self, z_m: f64) -> f64 {
        let vb = self.volume_fraction_b(z_m);
        self.k_a * (1.0_f64 - vb) + self.k_b * vb
    }
    /// Average (integrated) modulus over the full grading height.
    pub fn average_modulus(&self) -> f64 {
        self.e_a + (self.e_b - self.e_a) / (self.grading_exponent + 1.0_f64)
    }
}
/// Goldak double-ellipsoid heat source model.
///
/// Reference: Goldak et al. (1984) *Met. Trans. B*, 15, 299–305.
#[derive(Debug, Clone)]
pub struct GoldakHeatSource {
    /// Laser power (W).
    pub power: f64,
    /// Absorptivity (0–1).
    pub absorptivity: f64,
    /// Front ellipsoid semi-axis along scan direction (m).
    pub a_front: f64,
    /// Rear ellipsoid semi-axis along scan direction (m).
    pub a_rear: f64,
    /// Lateral semi-axis (m).
    pub b: f64,
    /// Depth semi-axis (m).
    pub c: f64,
    /// Fraction of heat deposited in the front ellipsoid.
    pub f_front: f64,
}
impl GoldakHeatSource {
    /// Construct a Goldak heat source with default fraction split (0.6 / 1.4).
    pub fn new(power: f64, absorptivity: f64, a_front: f64, a_rear: f64, b: f64, c: f64) -> Self {
        Self {
            power,
            absorptivity,
            a_front,
            a_rear,
            b,
            c,
            f_front: 0.6_f64,
        }
    }
    fn f_rear(&self) -> f64 {
        2.0_f64 - self.f_front
    }
    /// Volumetric heat flux (W/m³) at position `(x, y, z)` relative to the
    /// instantaneous heat source centre.  Positive *x* points in the scan
    /// direction, *z* points downward into the substrate.
    pub fn heat_flux(&self, x: f64, y: f64, z: f64) -> f64 {
        let q = self.absorptivity * self.power;
        let b2 = self.b * self.b;
        let c2 = self.c * self.c;
        let common = 6.0_f64 * f64::sqrt(3.0_f64) * q;
        if x >= 0.0_f64 {
            let a2 = self.a_front * self.a_front;
            let coeff =
                common * self.f_front / (PI * f64::sqrt(PI) * self.a_front * self.b * self.c);
            coeff * (-3.0_f64 * x * x / a2 - 3.0_f64 * y * y / b2 - 3.0_f64 * z * z / c2).exp()
        } else {
            let a2 = self.a_rear * self.a_rear;
            let coeff =
                common * self.f_rear() / (PI * f64::sqrt(PI) * self.a_rear * self.b * self.c);
            coeff * (-3.0_f64 * x * x / a2 - 3.0_f64 * y * y / b2 - 3.0_f64 * z * z / c2).exp()
        }
    }
    /// Peak heat flux at the source centre (W/m³).
    pub fn peak_heat_flux(&self) -> f64 {
        self.heat_flux(0.0_f64, 0.0_f64, 0.0_f64)
    }
}
/// Laser scanning strategy for each layer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanStrategy {
    /// Continuous raster scan lines parallel to the X-axis.
    Raster,
    /// Alternating 90° rotation between consecutive layers.
    Alternating90,
    /// 67° rotation between consecutive layers (reduces texture).
    Rotating67,
    /// Concentric inward spiral islands.
    Islands,
    /// Chessboard (checkerboard) pattern.
    Chessboard,
}
/// Anisotropic mechanical properties arising from the columnar-grain
/// microstructure of LPBF or DED parts.
#[derive(Debug, Clone)]
pub struct AnisotropicAM {
    /// Build direction.
    pub build_direction: BuildDirection,
    /// Young's modulus along the build direction (Pa).
    pub e_parallel: f64,
    /// Young's modulus transverse to the build direction (Pa).
    pub e_transverse: f64,
    /// Yield strength along the build direction (Pa).
    pub sigma_y_parallel: f64,
    /// Yield strength transverse to the build direction (Pa).
    pub sigma_y_transverse: f64,
    /// Texture coefficient (0 = random, 1 = perfectly aligned <001>).
    pub texture_coefficient: f64,
    /// Average columnar grain aspect ratio (length / width).
    pub grain_aspect_ratio: f64,
}
impl AnisotropicAM {
    /// Create an [`AnisotropicAM`] model.
    pub fn new(
        build_direction: BuildDirection,
        e_parallel: f64,
        e_transverse: f64,
        sigma_y_parallel: f64,
        sigma_y_transverse: f64,
        texture_coefficient: f64,
        grain_aspect_ratio: f64,
    ) -> Self {
        Self {
            build_direction,
            e_parallel,
            e_transverse,
            sigma_y_parallel,
            sigma_y_transverse,
            texture_coefficient,
            grain_aspect_ratio,
        }
    }
    /// Anisotropy index: (E_transverse − E_parallel) / E_parallel.
    pub fn anisotropy_index(&self) -> f64 {
        (self.e_transverse - self.e_parallel) / self.e_parallel
    }
    /// Effective modulus for a loading angle `theta` (radians) from the build
    /// direction, using a simple cosine interpolation.
    pub fn effective_modulus(&self, theta: f64) -> f64 {
        let ct = theta.cos();
        let st = theta.sin();
        1.0_f64 / (ct * ct / self.e_parallel + st * st / self.e_transverse)
    }
    /// Effective yield strength for loading angle `theta` (radians).
    pub fn effective_yield_strength(&self, theta: f64) -> f64 {
        let ct = theta.cos().abs();
        let st = theta.sin().abs();
        self.sigma_y_parallel * ct + self.sigma_y_transverse * st
    }
}
/// Directed Energy Deposition (DED) process model.
///
/// Captures clad geometry, dilution ratio and conditions favouring
/// epitaxial grain growth (substrate crystallographic continuity).
#[derive(Debug, Clone)]
pub struct DirectedEnergyDeposition {
    /// Laser power (W).
    pub laser_power: f64,
    /// Powder feed rate (kg/s).
    pub powder_feed_rate: f64,
    /// Scan speed (m/s).
    pub scan_speed: f64,
    /// Laser spot radius on the substrate (m).
    pub spot_radius: f64,
    /// Powder catchment efficiency (0–1).
    pub catchment_efficiency: f64,
    /// Substrate thermal conductivity (W/m·K).
    pub substrate_conductivity: f64,
    /// Melt pool absorptivity for the powder stream (0–1).
    pub absorptivity: f64,
}
impl DirectedEnergyDeposition {
    /// Construct a [`DirectedEnergyDeposition`] model.
    pub fn new(
        laser_power: f64,
        powder_feed_rate: f64,
        scan_speed: f64,
        spot_radius: f64,
        catchment_efficiency: f64,
        substrate_conductivity: f64,
        absorptivity: f64,
    ) -> Self {
        Self {
            laser_power,
            powder_feed_rate,
            scan_speed,
            spot_radius,
            catchment_efficiency,
            substrate_conductivity,
            absorptivity,
        }
    }
    /// Effective deposited power (W) after absorptivity losses.
    pub fn effective_power(&self) -> f64 {
        self.absorptivity * self.laser_power
    }
    /// Clad bead height (m) from a mass-balance estimate.
    ///
    /// `h ≈ ṁ_dep / (ρ_deposit × v × w)`
    /// Using a simplified density (7800 kg/m³ steel) and width ≈ 2·r.
    pub fn clad_height_m(&self, deposit_density: f64) -> f64 {
        let width = 2.0_f64 * self.spot_radius;
        let mass_flow = self.powder_feed_rate * self.catchment_efficiency;
        let volume_per_m = mass_flow / (deposit_density * self.scan_speed * width);
        volume_per_m.max(0.0_f64)
    }
    /// Dilution ratio D = (substrate melted volume) / (clad + substrate volume).
    ///
    /// Estimated from the energy balance: higher power → deeper penetration.
    pub fn dilution_ratio(&self) -> f64 {
        let e_s = self.effective_power() / (self.scan_speed * self.spot_radius).max(1.0e-15_f64);
        let e_ref = 5.0e6_f64;
        1.0_f64 - (-e_s / e_ref).exp()
    }
    /// True if conditions favour epitaxial grain growth (columnar solidification).
    ///
    /// Epitaxial growth occurs when the temperature gradient G is high and
    /// the solidification rate R is low (G/R > threshold).
    pub fn epitaxial_growth_likely(&self) -> bool {
        let g_over_r = self.effective_power()
            / (self.scan_speed.powi(2) * self.spot_radius * self.substrate_conductivity)
                .max(1.0e-30_f64);
        g_over_r > 1.0e6_f64
    }
    /// Powder utilisation efficiency (0–1): fraction of fed powder that is
    /// incorporated into the deposit.
    pub fn utilisation_efficiency(&self) -> f64 {
        self.catchment_efficiency.clamp(0.0_f64, 1.0_f64)
    }
}
/// Thermal history parameters for a melt-track or layer.
#[derive(Debug, Clone)]
pub struct ThermalHistory {
    /// Melt pool geometry.
    pub melt_pool: MeltPoolGeometry,
    /// Peak temperature reached in the melt pool (K).
    pub peak_temperature: f64,
    /// Solidification cooling rate (K/s).
    pub cooling_rate: f64,
    /// Thermal gradient at the solidification front (K/m).
    pub thermal_gradient: f64,
    /// Width of the heat-affected zone (HAZ) outside the melt pool (m).
    pub haz_width: f64,
    /// Whether recrystallisation is expected (temperature exceeded 0.5 × T_melt).
    pub recrystallised: bool,
}
impl ThermalHistory {
    /// Construct a [`ThermalHistory`] from process and material parameters using
    /// simplified analytical estimates.
    ///
    /// # Parameters
    /// - `layer`: layer model (power, speed, …)
    /// - `absorptivity`: laser absorptivity of the powder bed (0–1)
    /// - `thermal_conductivity`: bulk thermal conductivity (W/m·K)
    /// - `density`: material density (kg/m³)
    /// - `specific_heat`: specific heat capacity (J/kg·K)
    /// - `melting_point`: solidus temperature (K)
    /// - `ambient_temp`: substrate / build plate temperature (K)
    pub fn from_process(
        layer: &LayerModel,
        absorptivity: f64,
        thermal_conductivity: f64,
        density: f64,
        specific_heat: f64,
        melting_point: f64,
        ambient_temp: f64,
    ) -> Self {
        let effective_power = absorptivity * layer.laser_power;
        let _diffusivity = thermal_conductivity / (density * specific_heat);
        let r = layer.spot_radius.max(1.0e-6_f64);
        let peak_temperature =
            ambient_temp + effective_power / (2.0_f64 * PI * thermal_conductivity * r);
        let thermal_gradient = (peak_temperature - ambient_temp) / r;
        let cooling_rate = layer.scan_speed * thermal_gradient;
        let q_eff = effective_power;
        let half_length = q_eff
            / (2.0_f64 * PI * thermal_conductivity * (peak_temperature - ambient_temp))
                .max(1.0e-20_f64);
        let half_width = half_length * 0.6_f64;
        let depth = half_length * 0.4_f64;
        let haz_width = half_width * 1.5_f64;
        let recrystallised = peak_temperature > 0.5_f64 * melting_point;
        ThermalHistory {
            melt_pool: MeltPoolGeometry {
                half_length,
                half_width,
                depth,
            },
            peak_temperature,
            cooling_rate,
            thermal_gradient,
            haz_width,
            recrystallised,
        }
    }
    /// Solidification rate (m/s) — ratio of scan speed projected onto the
    /// liquid/solid interface.  Simplified: G_s = cooling_rate / thermal_gradient.
    pub fn solidification_rate(&self) -> f64 {
        if self.thermal_gradient.abs() < 1.0e-12_f64 {
            0.0_f64
        } else {
            self.cooling_rate / self.thermal_gradient
        }
    }
    /// Primary dendrite arm spacing (PDAS) in metres using the Hunt–Kurz
    /// power-law: λ₁ = A · (G · R)^(-0.25).
    ///
    /// Coefficient *A* ≈ 80 µm for typical Ti-6Al-4V (SI units).
    pub fn primary_dendrite_arm_spacing(&self) -> f64 {
        let a = 80.0e-6_f64;
        let gr = self.thermal_gradient * self.solidification_rate();
        if gr < 1.0e-20_f64 {
            a
        } else {
            a * gr.powf(-0.25_f64)
        }
    }
}
/// Binder jetting process model: powder spreading, binder saturation,
/// and sintering-induced shrinkage.
#[derive(Debug, Clone)]
pub struct BinderJetting {
    /// Nominal powder layer thickness (m).
    pub layer_thickness_m: f64,
    /// Powder bed packing density (0–1).
    pub packing_density: f64,
    /// Binder saturation ratio (volume of binder / volume of void).
    pub binder_saturation: f64,
    /// Powder particle mean diameter (m).
    pub particle_diameter_m: f64,
    /// Sintering shrinkage fraction (linear, 0–1).
    /// Typical: 0.15–0.22 for metals.
    pub sintering_shrinkage: f64,
    /// Green-part relative density (0–1) — after binding, before sinter.
    pub green_density: f64,
}
impl BinderJetting {
    /// Construct a [`BinderJetting`] model.
    pub fn new(
        layer_thickness_m: f64,
        packing_density: f64,
        binder_saturation: f64,
        particle_diameter_m: f64,
        sintering_shrinkage: f64,
    ) -> Self {
        let green_density = packing_density * binder_saturation.clamp(0.0_f64, 1.0_f64);
        Self {
            layer_thickness_m,
            packing_density: packing_density.clamp(0.0_f64, 1.0_f64),
            binder_saturation: binder_saturation.clamp(0.0_f64, 1.0_f64),
            particle_diameter_m,
            sintering_shrinkage: sintering_shrinkage.clamp(0.0_f64, 0.5_f64),
            green_density,
        }
    }
    /// Void fraction in the powder bed (1 − packing_density).
    pub fn void_fraction(&self) -> f64 {
        1.0_f64 - self.packing_density
    }
    /// Binder volume fraction relative to total layer volume.
    pub fn binder_volume_fraction(&self) -> f64 {
        self.void_fraction() * self.binder_saturation
    }
    /// Final linear dimension of a part after sintering.
    ///
    /// `L_final = L_green × (1 − sintering_shrinkage)`
    pub fn sintered_dimension(&self, green_dimension_m: f64) -> f64 {
        green_dimension_m * (1.0_f64 - self.sintering_shrinkage)
    }
    /// Volumetric shrinkage fraction after sintering.
    ///
    /// V_shrink / V_green = 1 − (1 − s)³
    pub fn volumetric_shrinkage(&self) -> f64 {
        let s = self.sintering_shrinkage;
        1.0_f64 - (1.0_f64 - s).powi(3)
    }
    /// Spreading speed limit estimate (m/s) based on Carman–Kozeny drag.
    ///
    /// Higher packing density or smaller particles lower the feasible
    /// recoater speed.
    pub fn max_spreading_speed_m_s(&self) -> f64 {
        let k_kozeny = 5.0_f64;
        let porosity = self.void_fraction().max(0.01_f64);
        let d = self.particle_diameter_m.max(1.0e-9_f64);
        porosity.powi(3) / (k_kozeny * (1.0_f64 - porosity).powi(2) * d)
    }
}
/// Steady-state Rosenthal temperature field for a moving point heat source.
///
/// Reference: Rosenthal (1946) *Trans. ASME*, 68, 849–866.
#[derive(Debug, Clone)]
pub struct RosenthalSolution {
    /// Laser power (W).
    pub power: f64,
    /// Absorptivity (0–1).
    pub absorptivity: f64,
    /// Scan speed (m/s).
    pub scan_speed: f64,
    /// Thermal conductivity (W/m·K).
    pub thermal_conductivity: f64,
    /// Thermal diffusivity (m²/s).
    pub diffusivity: f64,
    /// Ambient / preheat temperature (K).
    pub ambient_temp: f64,
}
impl RosenthalSolution {
    /// Construct a [`RosenthalSolution`].
    pub fn new(
        power: f64,
        absorptivity: f64,
        scan_speed: f64,
        thermal_conductivity: f64,
        diffusivity: f64,
        ambient_temp: f64,
    ) -> Self {
        Self {
            power,
            absorptivity,
            scan_speed,
            thermal_conductivity,
            diffusivity,
            ambient_temp,
        }
    }
    /// Temperature (K) at position `(xi, y, z)` in the moving heat-source
    /// frame.  `xi = x − v·t` is the distance behind the laser spot.
    ///
    /// Returns `f64::INFINITY` at the exact source centre.
    pub fn temperature(&self, xi: f64, y: f64, z: f64) -> f64 {
        let r = (xi * xi + y * y + z * z).sqrt();
        if r < 1.0e-15_f64 {
            return f64::INFINITY;
        }
        let q = self.absorptivity * self.power;
        let v = self.scan_speed;
        let alpha = self.diffusivity;
        let k = self.thermal_conductivity;
        let exponent = -v / (2.0_f64 * alpha) * (r + xi);
        self.ambient_temp + q / (2.0_f64 * PI * k * r) * exponent.exp()
    }
    /// Melt-pool width (m) at the surface (`z = 0`) given the melting point
    /// temperature `t_melt` (K).  Uses a bisection search.
    pub fn melt_pool_width(&self, t_melt: f64) -> f64 {
        let mut lo = 1.0e-9_f64;
        let mut hi = 1.0e-2_f64;
        for _ in 0..60 {
            let mid = 0.5_f64 * (lo + hi);
            let t = self.temperature(0.0_f64, mid, 0.0_f64);
            if t > t_melt {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        2.0_f64 * 0.5_f64 * (lo + hi)
    }
}
/// Principal build direction for an AM part.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuildDirection {
    /// Positive Z direction (standard upward build).
    PlusZ,
    /// Build tilted 45° from the Z axis toward X.
    Tilted45X,
    /// Horizontal build (powder-bed binder jetting or FDM on side).
    Horizontal,
    /// Custom direction stored as a unit vector index.
    Custom,
}
/// Process-window model mapping energy density to part quality metrics.
#[derive(Debug, Clone)]
pub struct ProcessWindow {
    /// Material density of the fully dense material (kg/m³).
    pub theoretical_density: f64,
    /// Lower bound energy density for good densification (J/m³).
    pub e_low: f64,
    /// Upper bound energy density before keyhole formation (J/m³).
    pub e_high: f64,
    /// Vickers hardness at theoretical density (HV).
    pub hardness_full_dense: f64,
    /// Fitted Hall–Petch coefficient (HV·m^0.5).
    pub hall_petch_k: f64,
    /// Grain size at optimal energy density (m).
    pub grain_size_opt: f64,
}
impl ProcessWindow {
    /// Construct a [`ProcessWindow`].
    pub fn new(
        theoretical_density: f64,
        e_low: f64,
        e_high: f64,
        hardness_full_dense: f64,
        hall_petch_k: f64,
        grain_size_opt: f64,
    ) -> Self {
        Self {
            theoretical_density,
            e_low,
            e_high,
            hardness_full_dense,
            hall_petch_k,
            grain_size_opt,
        }
    }
    /// Predicted relative density (0–1) as a function of `energy_density` (J/m³).
    ///
    /// Uses a sigmoidal model: `ρ_rel = 1 / (1 + exp(−k · (E − E_mid)))`.
    pub fn relative_density(&self, energy_density: f64) -> f64 {
        let e_mid = 0.5_f64 * (self.e_low + self.e_high);
        let k = 5.0_f64 / (self.e_high - self.e_low).max(1.0_f64);
        1.0_f64 / (1.0_f64 + (-k * (energy_density - e_mid)).exp())
    }
    /// Predicted Vickers hardness at a given `energy_density` (J/m³).
    ///
    /// Combines density-based degradation with a Hall–Petch grain-size term.
    pub fn hardness_prediction(&self, energy_density: f64, grain_size: f64) -> f64 {
        let rel = self.relative_density(energy_density);
        let hp = self.hall_petch_k / grain_size.sqrt();
        rel * (self.hardness_full_dense + hp - self.hall_petch_k / self.grain_size_opt.sqrt())
    }
    /// True if `energy_density` lies within the acceptable process window.
    pub fn is_in_window(&self, energy_density: f64) -> bool {
        energy_density >= self.e_low && energy_density <= self.e_high
    }
    /// Window width in J/m³.
    pub fn window_width(&self) -> f64 {
        (self.e_high - self.e_low).max(0.0_f64)
    }
}
/// Pre-deformation (negative distortion) to compensate for warpage in LPBF.
///
/// The method scales every point's displacement by a compensation factor so
/// that after springback the part converges to the nominal geometry.
#[derive(Debug, Clone)]
pub struct DistortionCompensation {
    /// Scaling factor applied to predicted displacements (typically −1 to −2).
    pub scale_factor: f64,
    /// Maximum allowable pre-deformation magnitude (m).
    pub max_deformation_m: f64,
    /// Convergence tolerance (m) — iterations stop when RMS update < tol.
    pub tolerance_m: f64,
    /// Number of compensation iterations performed.
    pub iteration_count: usize,
}
impl DistortionCompensation {
    /// Construct a [`DistortionCompensation`] model.
    ///
    /// `scale_factor` < 0 inverts the predicted distortion field.
    pub fn new(scale_factor: f64, max_deformation_m: f64, tolerance_m: f64) -> Self {
        Self {
            scale_factor,
            max_deformation_m,
            tolerance_m,
            iteration_count: 0,
        }
    }
    /// Apply compensation to a list of nodal displacement predictions.
    ///
    /// Each element of `displacements_m` is the predicted warpage for one
    /// node.  The function returns the compensated (pre-deformed) geometry
    /// offsets.
    pub fn apply(&self, displacements_m: &[f64]) -> Vec<f64> {
        displacements_m
            .iter()
            .map(|d| {
                let comp = self.scale_factor * d;
                comp.clamp(-self.max_deformation_m, self.max_deformation_m)
            })
            .collect()
    }
    /// Iterative loop: repeatedly apply compensation and check convergence.
    ///
    /// Returns the final compensated displacement field after at most
    /// `max_iter` iterations.
    pub fn iterate(&mut self, initial_displacements: &[f64], max_iter: usize) -> Vec<f64> {
        let mut current = initial_displacements.to_vec();
        self.iteration_count = 0;
        for _i in 0..max_iter {
            let next = self.apply(&current);
            let rms: f64 = next
                .iter()
                .zip(current.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                / (next.len().max(1) as f64);
            current = next;
            self.iteration_count += 1;
            if rms.sqrt() < self.tolerance_m {
                break;
            }
        }
        current
    }
    /// RMS magnitude of a displacement field (m).
    pub fn rms_displacement(displacements_m: &[f64]) -> f64 {
        if displacements_m.is_empty() {
            return 0.0_f64;
        }
        let sum_sq: f64 = displacements_m.iter().map(|d| d * d).sum();
        (sum_sq / displacements_m.len() as f64).sqrt()
    }
}
/// Scan-speed vs laser-power optimiser for minimal porosity and good surface.
///
/// Sweeps a discrete grid of (power, speed) pairs and scores each point
/// using a combined porosity-and-surface-roughness cost function.
#[derive(Debug, Clone)]
pub struct ProcessWindowOptimizer {
    /// Minimum laser power to evaluate (W).
    pub power_min: f64,
    /// Maximum laser power to evaluate (W).
    pub power_max: f64,
    /// Minimum scan speed to evaluate (m/s).
    pub speed_min: f64,
    /// Maximum scan speed to evaluate (m/s).
    pub speed_max: f64,
    /// Layer thickness (m) — fixed during optimisation.
    pub layer_thickness: f64,
    /// Hatch spacing (m) — fixed during optimisation.
    pub hatch_spacing: f64,
    /// Lower energy-density threshold for dense parts (J/m³).
    pub e_low: f64,
    /// Upper energy-density threshold before keyhole (J/m³).
    pub e_high: f64,
}
impl ProcessWindowOptimizer {
    /// Construct a [`ProcessWindowOptimizer`].
    pub fn new(
        power_min: f64,
        power_max: f64,
        speed_min: f64,
        speed_max: f64,
        layer_thickness: f64,
        hatch_spacing: f64,
        e_low: f64,
        e_high: f64,
    ) -> Self {
        Self {
            power_min,
            power_max,
            speed_min,
            speed_max,
            layer_thickness,
            hatch_spacing,
            e_low,
            e_high,
        }
    }
    /// Compute volumetric energy density for a (power, speed) pair (J/m³).
    pub fn energy_density(&self, power: f64, speed: f64) -> f64 {
        power / (speed * self.hatch_spacing * self.layer_thickness).max(1.0e-30_f64)
    }
    /// Porosity score for a given energy density (lower is better, 0 = optimal).
    ///
    /// Uses a quadratic penalty outside the \[e_low, e_high\] window.
    pub fn porosity_score(&self, ev: f64) -> f64 {
        if ev < self.e_low {
            ((self.e_low - ev) / self.e_low).powi(2)
        } else if ev > self.e_high {
            ((ev - self.e_high) / self.e_high).powi(2)
        } else {
            0.0_f64
        }
    }
    /// Surface roughness score for a given scan speed (lower is better).
    ///
    /// Higher scan speeds produce rougher surfaces (balling) above a threshold.
    pub fn roughness_score(&self, speed: f64) -> f64 {
        let v_ref = 0.8_f64 * self.speed_max;
        if speed > v_ref {
            ((speed - v_ref) / (self.speed_max - v_ref + 1.0e-30_f64)).powi(2)
        } else {
            0.0_f64
        }
    }
    /// Combined cost for a (power, speed) pair.
    pub fn cost(&self, power: f64, speed: f64) -> f64 {
        let ev = self.energy_density(power, speed);
        self.porosity_score(ev) + 0.5_f64 * self.roughness_score(speed)
    }
    /// Sweep `n_points × n_points` grid and return the (power, speed) with the
    /// lowest combined cost.
    pub fn optimise(&self, n_points: usize) -> (f64, f64) {
        let n = n_points.max(2);
        let mut best_cost = f64::INFINITY;
        let mut best_power = self.power_min;
        let mut best_speed = self.speed_min;
        for ip in 0..n {
            let p = self.power_min + ip as f64 * (self.power_max - self.power_min) / (n - 1) as f64;
            for iv in 0..n {
                let v =
                    self.speed_min + iv as f64 * (self.speed_max - self.speed_min) / (n - 1) as f64;
                let c = self.cost(p, v);
                if c < best_cost {
                    best_cost = c;
                    best_power = p;
                    best_speed = v;
                }
            }
        }
        (best_power, best_speed)
    }
}
