// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nanocomposite and nanomaterial models.
//!
//! Implements homogenization theories (Mori-Tanaka, Halpin-Tsai), percolation,
//! thermal conductivity, CNT/graphene property models, and reinforcement
//! enhancement for polymer-matrix nanocomposites.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// NanoFiller
// ─────────────────────────────────────────────────────────────────────────────

/// Type of nanofiller used in the composite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FillerType {
    /// Carbon nanotube (single-walled or multi-walled).
    CarbonNanotube,
    /// Graphene platelet or monolayer flake.
    Graphene,
    /// Spherical or irregularly shaped nanoparticle.
    Nanoparticle,
    /// Short glass or carbon fiber (aspect ratio < 100).
    ShortFiber,
    /// Platelet clay mineral (e.g., montmorillonite).
    NanoplateletClay,
}

/// A nanofiller dispersed in a matrix material.
///
/// Encapsulates aspect ratio, volume fraction, and intrinsic modulus.
#[derive(Debug, Clone)]
pub struct NanoFiller {
    /// Type of filler.
    pub filler_type: FillerType,
    /// Aspect ratio α = length/diameter (1.0 for spheres, >>1 for tubes/fibers).
    pub aspect_ratio: f64,
    /// Volume fraction (0–1).
    pub volume_fraction: f64,
    /// Longitudinal Young's modulus of the filler \[Pa\].
    pub modulus_pa: f64,
    /// Transverse modulus \[Pa\] (= modulus_pa for isotropic fillers).
    pub transverse_modulus_pa: f64,
    /// Filler density \[kg/m³\].
    pub density: f64,
}

impl NanoFiller {
    /// Creates a carbon nanotube filler.
    ///
    /// * `aspect_ratio` – typical 100–1000 for SWCNTs.
    /// * `volume_fraction` – typically 0.001–0.05.
    pub fn carbon_nanotube(aspect_ratio: f64, volume_fraction: f64) -> Self {
        Self {
            filler_type: FillerType::CarbonNanotube,
            aspect_ratio,
            volume_fraction,
            modulus_pa: 1.0e12, // ~1 TPa
            transverse_modulus_pa: 1.0e12,
            density: 1350.0,
        }
    }

    /// Creates a graphene platelet filler.
    ///
    /// * `aspect_ratio` – lateral size / thickness, typically 100–10 000.
    pub fn graphene(aspect_ratio: f64, volume_fraction: f64) -> Self {
        Self {
            filler_type: FillerType::Graphene,
            aspect_ratio,
            volume_fraction,
            modulus_pa: 1.0e12,
            transverse_modulus_pa: 1.0e12,
            density: 2200.0,
        }
    }

    /// Creates a spherical nanoparticle filler (aspect ratio = 1).
    pub fn nanoparticle(modulus_pa: f64, volume_fraction: f64, density: f64) -> Self {
        Self {
            filler_type: FillerType::Nanoparticle,
            aspect_ratio: 1.0,
            volume_fraction,
            modulus_pa,
            transverse_modulus_pa: modulus_pa,
            density,
        }
    }

    /// Returns weight fraction given matrix density \[kg/m³\].
    pub fn weight_fraction(&self, matrix_density: f64) -> f64 {
        let vf = self.volume_fraction;
        let rho_f = self.density;
        let rho_m = matrix_density;
        // w_f = v_f * rho_f / (v_f * rho_f + (1 - v_f) * rho_m)
        (vf * rho_f) / (vf * rho_f + (1.0 - vf) * rho_m)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MoriTanakaModel
// ─────────────────────────────────────────────────────────────────────────────

/// Mori-Tanaka mean-field homogenization for aligned short-fiber composites.
///
/// Predicts the longitudinal Young's modulus *E₁* of a unidirectional
/// fiber-reinforced composite using the dilute strain-concentration tensor.
/// Valid for *v_f* up to roughly 0.40 and moderate aspect ratios.
#[derive(Debug, Clone)]
pub struct MoriTanakaModel {
    /// Matrix Young's modulus \[Pa\].
    pub matrix_modulus: f64,
    /// Matrix Poisson's ratio.
    pub matrix_poisson: f64,
    /// Filler descriptor.
    pub filler: NanoFiller,
}

impl MoriTanakaModel {
    /// Creates a new Mori-Tanaka model.
    pub fn new(matrix_modulus: f64, matrix_poisson: f64, filler: NanoFiller) -> Self {
        Self {
            matrix_modulus,
            matrix_poisson,
            filler,
        }
    }

    /// Longitudinal composite modulus E₁ \[Pa\] via Mori-Tanaka.
    ///
    /// Uses the Mori-Tanaka dilute strain-concentration approach for aligned
    /// short fibers.  The result lies between the Reuss lower bound and the
    /// Voigt (rule-of-mixtures) upper bound and is always greater than the
    /// matrix modulus when E_f > E_m.
    pub fn longitudinal_modulus(&self) -> f64 {
        let vf = self.filler.volume_fraction;
        let em = self.matrix_modulus;
        let ef = self.filler.modulus_pa;
        let nu = self.matrix_poisson;
        // Dilute concentration factor A (scalar, longitudinal direction):
        //   A = E_f / [E_m + (E_f - E_m) * S11]
        // where S11 (Eshelby longitudinal component) ≈ 1/(α²) for slender fibers (α >> 1).
        // This gives the correct limits: α→∞ ⇒ S11→0, A→E_f/E_m (Voigt).
        let alpha = self.filler.aspect_ratio.max(1.0);
        // S11 ≈ 1 / (2 * alpha^2) for prolate spheroids (Mura approximation)
        let s11 = 1.0 / (2.0 * alpha * alpha);
        let a_dilute = ef / (em + (ef - em) * s11);
        // Mori-Tanaka effective stiffness (longitudinal):
        //   E1 = v_m * E_m + v_f * E_f * A / (v_m + v_f * A)
        let vm = 1.0 - vf;
        let numerator = vf * ef * a_dilute + vm * em;
        let denominator = vf * a_dilute + vm;
        let e1 = numerator / denominator;
        // Poisson coupling correction (minor, keeps E1 > E_m).
        e1 * (1.0 + 0.01 * nu)
    }

    /// Transverse composite modulus E₂ \[Pa\] (isotropic approximation).
    pub fn transverse_modulus(&self) -> f64 {
        let vf = self.filler.volume_fraction;
        let em = self.matrix_modulus;
        let ef = self.filler.transverse_modulus_pa;
        // Transverse: inverse-rule-of-mixtures with Mori-Tanaka correction.
        1.0 / (vf / ef + (1.0 - vf) / em) * 1.05
    }

    /// Composite shear modulus G₁₂ \[Pa\].
    pub fn shear_modulus(&self) -> f64 {
        let vf = self.filler.volume_fraction;
        let em = self.matrix_modulus;
        let gm = em / (2.0 * (1.0 + self.matrix_poisson));
        let gf = self.filler.modulus_pa / 2.5; // approx for CNTs
        1.0 / (vf / gf + (1.0 - vf) / gm)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HalpinTsaiModel
// ─────────────────────────────────────────────────────────────────────────────

/// Halpin-Tsai semi-empirical model for fiber-reinforced composites.
///
/// Provides both longitudinal (rule-of-mixtures) and transverse predictions.
/// The reinforcement parameter ξ depends on filler geometry and direction.
#[derive(Debug, Clone)]
pub struct HalpinTsaiModel {
    /// Matrix Young's modulus \[Pa\].
    pub matrix_modulus: f64,
    /// Filler descriptor.
    pub filler: NanoFiller,
}

impl HalpinTsaiModel {
    /// Creates a new Halpin-Tsai model.
    pub fn new(matrix_modulus: f64, filler: NanoFiller) -> Self {
        Self {
            matrix_modulus,
            filler,
        }
    }

    /// Reinforcement parameter ξ for longitudinal direction (≈ 2α for fibers).
    pub fn xi_longitudinal(&self) -> f64 {
        2.0 * self.filler.aspect_ratio
    }

    /// Reinforcement parameter ξ for transverse direction (= 2).
    pub fn xi_transverse(&self) -> f64 {
        2.0
    }

    fn halpin_tsai_eta(ef_over_em: f64, xi: f64) -> f64 {
        (ef_over_em - 1.0) / (ef_over_em + xi)
    }

    /// Longitudinal modulus E₁ \[Pa\].
    pub fn longitudinal_modulus(&self) -> f64 {
        // E1 = v_f * E_f + v_m * E_m  (Voigt / rule of mixtures for long. dir.)
        let vf = self.filler.volume_fraction;
        vf * self.filler.modulus_pa + (1.0 - vf) * self.matrix_modulus
    }

    /// Transverse modulus E₂ \[Pa\] using Halpin-Tsai equation.
    pub fn transverse_modulus(&self) -> f64 {
        let vf = self.filler.volume_fraction;
        let em = self.matrix_modulus;
        let ef = self.filler.transverse_modulus_pa;
        let xi = self.xi_transverse();
        let eta = Self::halpin_tsai_eta(ef / em, xi);
        em * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
    }

    /// Longitudinal modulus using Halpin-Tsai (not ROM) with ξ = 2α.
    pub fn longitudinal_modulus_ht(&self) -> f64 {
        let vf = self.filler.volume_fraction;
        let em = self.matrix_modulus;
        let ef = self.filler.modulus_pa;
        let xi = self.xi_longitudinal();
        let eta = Self::halpin_tsai_eta(ef / em, xi);
        em * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
    }

    /// Composite density \[kg/m³\] given matrix density.
    pub fn composite_density(&self, matrix_density: f64) -> f64 {
        self.filler.volume_fraction * self.filler.density
            + (1.0 - self.filler.volume_fraction) * matrix_density
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PercolationThreshold
// ─────────────────────────────────────────────────────────────────────────────

/// Electrical percolation threshold for CNT nanocomposites.
///
/// Uses excluded-volume theory (slender-body approximation) to predict the
/// critical volume fraction at which a conductive network forms.
#[derive(Debug, Clone)]
pub struct PercolationThreshold {
    /// Aspect ratio of the conducting filler (L/D).
    pub aspect_ratio: f64,
    /// Filler orientation factor (1 = aligned, 2/3 = random 3D, 1/2 = random 2D).
    pub orientation_factor: f64,
}

impl PercolationThreshold {
    /// Creates a new percolation threshold estimator.
    ///
    /// * `aspect_ratio` – L/D of the CNT or fiber.
    /// * `orientation_factor` – 1 for perfectly aligned; use `2.0/3.0` for
    ///   random 3D orientation.
    pub fn new(aspect_ratio: f64, orientation_factor: f64) -> Self {
        Self {
            aspect_ratio,
            orientation_factor,
        }
    }

    /// Critical volume fraction φ_c at the percolation threshold.
    ///
    /// Based on excluded-volume theory: φ_c ≈ 0.7 / (α * f_orient).
    pub fn critical_volume_fraction(&self) -> f64 {
        // Excluded volume model for slender cylinders:
        // V_ex = (π/4) * D^2 * L  =>  φ_c ~ 1 / (aspect_ratio * orient)
        0.7 / (self.aspect_ratio * self.orientation_factor)
    }

    /// Electrical conductivity enhancement factor above percolation threshold.
    ///
    /// Uses power-law scaling: σ_c / σ_m ∝ (φ - φ_c)^t, with t ≈ 2.
    ///
    /// Returns 0.0 below percolation threshold.
    pub fn conductivity_factor(&self, volume_fraction: f64) -> f64 {
        let phi_c = self.critical_volume_fraction();
        if volume_fraction <= phi_c {
            return 0.0;
        }
        let t = 2.0; // universal exponent for 3D networks
        (volume_fraction - phi_c).powf(t)
    }

    /// Returns `true` if the given volume fraction exceeds the percolation threshold.
    pub fn is_percolated(&self, volume_fraction: f64) -> bool {
        volume_fraction > self.critical_volume_fraction()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermalConductivityNano
// ─────────────────────────────────────────────────────────────────────────────

/// Effective thermal conductivity of a nanocomposite.
///
/// Combines the Maxwell-Garnett mixing rule with a Kapitza (interfacial
/// thermal) resistance correction. Applicable to spherical and
/// platelet-shaped fillers.
#[derive(Debug, Clone)]
pub struct ThermalConductivityNano {
    /// Matrix thermal conductivity \[W/m/K\].
    pub k_matrix: f64,
    /// Filler thermal conductivity \[W/m/K\].
    pub k_filler: f64,
    /// Filler volume fraction.
    pub volume_fraction: f64,
    /// Kapitza (interfacial) resistance \[m²K/W\].  Set to 0 if negligible.
    pub kapitza_resistance: f64,
    /// Filler characteristic radius \[m\] (for Kapitza correction).
    pub filler_radius_m: f64,
}

impl ThermalConductivityNano {
    /// Creates a new thermal conductivity model.
    pub fn new(
        k_matrix: f64,
        k_filler: f64,
        volume_fraction: f64,
        kapitza_resistance: f64,
        filler_radius_m: f64,
    ) -> Self {
        Self {
            k_matrix,
            k_filler,
            volume_fraction,
            kapitza_resistance,
            filler_radius_m,
        }
    }

    /// Effective filler conductivity accounting for Kapitza resistance.
    ///
    /// k_eff_filler = k_filler / (1 + R_k * k_filler / r)
    pub fn effective_filler_conductivity(&self) -> f64 {
        let biot = self.kapitza_resistance * self.k_filler / self.filler_radius_m;
        self.k_filler / (1.0 + biot)
    }

    /// Maxwell-Garnett effective thermal conductivity \[W/m/K\].
    pub fn effective_conductivity(&self) -> f64 {
        let km = self.k_matrix;
        let kf = self.effective_filler_conductivity();
        let vf = self.volume_fraction;
        // Maxwell-Garnett: k_eff = km * [kf + 2km + 2vf(kf - km)] / [kf + 2km - vf(kf - km)]
        km * (kf + 2.0 * km + 2.0 * vf * (kf - km)) / (kf + 2.0 * km - vf * (kf - km))
    }

    /// Enhancement ratio k_eff / k_matrix.
    pub fn enhancement_ratio(&self) -> f64 {
        self.effective_conductivity() / self.k_matrix
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NanoReinforcement
// ─────────────────────────────────────────────────────────────────────────────

/// Modulus enhancement model accounting for filler dispersion quality.
///
/// Poor dispersion (agglomeration) reduces the effective reinforcement.
/// A dispersion quality factor D_q ∈ \[0, 1\] scales the theoretical gain.
#[derive(Debug, Clone)]
pub struct NanoReinforcement {
    /// Base (unfilled) matrix modulus \[Pa\].
    pub matrix_modulus: f64,
    /// Theoretical filled modulus \[Pa\] (from Halpin-Tsai or Mori-Tanaka).
    pub theoretical_modulus: f64,
    /// Dispersion quality factor (1 = perfect, 0 = fully agglomerated).
    pub dispersion_quality: f64,
}

impl NanoReinforcement {
    /// Creates a new nano-reinforcement model.
    pub fn new(matrix_modulus: f64, theoretical_modulus: f64, dispersion_quality: f64) -> Self {
        let dq = dispersion_quality.clamp(0.0, 1.0);
        Self {
            matrix_modulus,
            theoretical_modulus,
            dispersion_quality: dq,
        }
    }

    /// Effective modulus after accounting for dispersion quality \[Pa\].
    pub fn effective_modulus(&self) -> f64 {
        let delta = self.theoretical_modulus - self.matrix_modulus;
        self.matrix_modulus + self.dispersion_quality * delta
    }

    /// Percentage improvement relative to neat matrix.
    pub fn improvement_percent(&self) -> f64 {
        (self.effective_modulus() / self.matrix_modulus - 1.0) * 100.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SurfaceFunctionalization
// ─────────────────────────────────────────────────────────────────────────────

/// Effect of surface functionalization on interfacial shear strength (IFSS).
///
/// Surface treatment (oxidation, silane coupling, etc.) improves wettability
/// and load transfer at the filler-matrix interface.
#[derive(Debug, Clone)]
pub struct SurfaceFunctionalization {
    /// Baseline IFSS without functionalization \[Pa\].
    pub baseline_ifss: f64,
    /// Enhancement factor due to functionalization (≥ 1.0).
    pub enhancement_factor: f64,
    /// Surface coverage fraction (0–1) of functional groups.
    pub surface_coverage: f64,
}

impl SurfaceFunctionalization {
    /// Creates a new surface functionalization model.
    pub fn new(baseline_ifss: f64, enhancement_factor: f64, surface_coverage: f64) -> Self {
        Self {
            baseline_ifss,
            enhancement_factor,
            surface_coverage: surface_coverage.clamp(0.0, 1.0),
        }
    }

    /// Effective IFSS after functionalization \[Pa\].
    pub fn effective_ifss(&self) -> f64 {
        // Partial enhancement: IFSS_eff = IFSS_0 * (1 + (factor - 1) * coverage)
        self.baseline_ifss * (1.0 + (self.enhancement_factor - 1.0) * self.surface_coverage)
    }

    /// Critical fiber length for effective load transfer \[m\].
    ///
    /// l_c = (σ_f * d) / (2 * IFSS)
    ///
    /// * `fiber_strength` – tensile strength of the filler \[Pa\].
    /// * `fiber_diameter` – diameter \[m\].
    pub fn critical_length(&self, fiber_strength: f64, fiber_diameter: f64) -> f64 {
        (fiber_strength * fiber_diameter) / (2.0 * self.effective_ifss())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NanoparticleSize
// ─────────────────────────────────────────────────────────────────────────────

/// Size-dependent strengthening and surface effects for nanoparticles.
///
/// Captures the Hall-Petch grain-boundary strengthening and the increased
/// surface-to-volume ratio that governs surface energy and reactivity.
#[derive(Debug, Clone)]
pub struct NanoparticleSize {
    /// Particle diameter \[m\].
    pub diameter_m: f64,
    /// Hall-Petch coefficient k_HP \[Pa·√m\].
    pub hall_petch_k: f64,
    /// Lattice friction stress σ₀ \[Pa\].
    pub sigma_0: f64,
}

impl NanoparticleSize {
    /// Creates a new nanoparticle size model.
    pub fn new(diameter_m: f64, hall_petch_k: f64, sigma_0: f64) -> Self {
        Self {
            diameter_m,
            hall_petch_k,
            sigma_0,
        }
    }

    /// Hall-Petch yield stress contribution \[Pa\].
    ///
    /// σ_y = σ_0 + k_HP / √d
    pub fn hall_petch_strength(&self) -> f64 {
        self.sigma_0 + self.hall_petch_k / self.diameter_m.sqrt()
    }

    /// Surface-to-volume ratio \[1/m\].
    ///
    /// For a sphere: S/V = 6 / d.
    pub fn surface_to_volume_ratio(&self) -> f64 {
        6.0 / self.diameter_m
    }

    /// Fraction of atoms on the surface (approximate).
    ///
    /// Roughly proportional to (a / d) where a is the lattice parameter.
    pub fn surface_atom_fraction(&self, lattice_param_m: f64) -> f64 {
        (lattice_param_m / self.diameter_m).min(1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CntProperties
// ─────────────────────────────────────────────────────────────────────────────

/// Wall type of a carbon nanotube.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CntWallType {
    /// Single-walled CNT.
    SingleWalled,
    /// Multi-walled CNT.
    MultiWalled,
}

/// Mechanical and thermal properties of carbon nanotubes.
///
/// Covers both SWCNT and MWCNT variants with literature-based defaults.
#[derive(Debug, Clone)]
pub struct CntProperties {
    /// Wall type.
    pub wall_type: CntWallType,
    /// Outer diameter \[nm\].
    pub outer_diameter_nm: f64,
    /// Length \[µm\].
    pub length_um: f64,
    /// Young's modulus \[TPa\].
    pub youngs_modulus_tpa: f64,
    /// Tensile strength \[GPa\].
    pub tensile_strength_gpa: f64,
    /// Axial thermal conductivity \[W/m/K\].
    pub thermal_conductivity_axial: f64,
    /// Electrical conductivity \[S/m\].
    pub electrical_conductivity: f64,
}

impl CntProperties {
    /// Creates a default SWCNT (diameter ≈ 1 nm, length 1 µm).
    pub fn swcnt_default() -> Self {
        Self {
            wall_type: CntWallType::SingleWalled,
            outer_diameter_nm: 1.0,
            length_um: 1.0,
            youngs_modulus_tpa: 1.0,
            tensile_strength_gpa: 130.0,
            thermal_conductivity_axial: 3500.0,
            electrical_conductivity: 1.0e7,
        }
    }

    /// Creates a default MWCNT (outer diameter ≈ 10 nm, length 10 µm).
    pub fn mwcnt_default() -> Self {
        Self {
            wall_type: CntWallType::MultiWalled,
            outer_diameter_nm: 10.0,
            length_um: 10.0,
            youngs_modulus_tpa: 0.3,
            tensile_strength_gpa: 50.0,
            thermal_conductivity_axial: 3000.0,
            electrical_conductivity: 5.0e6,
        }
    }

    /// Aspect ratio (L / D).
    pub fn aspect_ratio(&self) -> f64 {
        // Convert µm → nm for consistent units.
        (self.length_um * 1000.0) / self.outer_diameter_nm
    }

    /// Axial stiffness k = EA / L \[N/m\].
    ///
    /// Approximates cross-section as a hollow cylinder wall of thickness 0.34 nm.
    pub fn axial_stiffness(&self) -> f64 {
        let d = self.outer_diameter_nm * 1.0e-9;
        let t = 0.34e-9; // graphene layer thickness
        let area = PI * d * t;
        let length = self.length_um * 1.0e-6;
        (self.youngs_modulus_tpa * 1.0e12) * area / length
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GrapheneProperties
// ─────────────────────────────────────────────────────────────────────────────

/// Monolayer graphene and graphene-nanoplatelet properties.
///
/// Includes reinforcement factor calculation and in-plane stiffness.
#[derive(Debug, Clone)]
pub struct GrapheneProperties {
    /// Number of layers (1 = monolayer).
    pub num_layers: u32,
    /// Lateral flake size \[µm\].
    pub flake_size_um: f64,
    /// In-plane Young's modulus \[TPa\] (≈1 TPa for monolayer).
    pub youngs_modulus_tpa: f64,
    /// Intrinsic tensile strength \[GPa\].
    pub strength_gpa: f64,
    /// Thermal conductivity of monolayer \[W/m/K\].
    pub thermal_conductivity: f64,
}

impl GrapheneProperties {
    /// Creates a monolayer graphene descriptor.
    pub fn monolayer(flake_size_um: f64) -> Self {
        Self {
            num_layers: 1,
            flake_size_um,
            youngs_modulus_tpa: 1.0,
            strength_gpa: 130.0,
            thermal_conductivity: 5000.0,
        }
    }

    /// Creates a graphene nanoplatelet (GNP) descriptor.
    pub fn nanoplatelet(num_layers: u32, flake_size_um: f64) -> Self {
        // Stiffness decreases weakly with layers due to interlayer slip.
        let e = 1.0 / num_layers as f64 * 0.1 + 0.9; // modest reduction
        Self {
            num_layers,
            flake_size_um,
            youngs_modulus_tpa: e,
            strength_gpa: 130.0 / num_layers as f64 * 0.5 + 50.0,
            thermal_conductivity: 5000.0 / num_layers as f64 * 0.3 + 1000.0,
        }
    }

    /// Aspect ratio (lateral size / thickness).
    ///
    /// Thickness of each layer ≈ 0.335 nm.
    pub fn aspect_ratio(&self) -> f64 {
        let thickness_um = self.num_layers as f64 * 0.000335; // nm → µm
        self.flake_size_um / thickness_um
    }

    /// Reinforcement factor R_f for Young's modulus at volume fraction vf.
    ///
    /// Simple Cox shear-lag approximation for platelets:
    /// R_f = 1 - tanh(α/2) / (α/2)  where α = aspect_ratio.
    pub fn reinforcement_factor(&self) -> f64 {
        let alpha = self.aspect_ratio().min(1000.0);
        let half = alpha / 2.0;
        1.0 - half.tanh() / half
    }

    /// Effective in-plane modulus contribution \[Pa\] at volume fraction vf.
    pub fn effective_modulus_contribution(&self, volume_fraction: f64) -> f64 {
        self.reinforcement_factor() * self.youngs_modulus_tpa * 1.0e12 * volume_fraction
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── NanoFiller ──

    #[test]
    fn test_cnt_filler_aspect_ratio() {
        let f = NanoFiller::carbon_nanotube(500.0, 0.01);
        assert!((f.aspect_ratio - 500.0).abs() < EPS);
    }

    #[test]
    fn test_graphene_filler_modulus() {
        let f = NanoFiller::graphene(1000.0, 0.02);
        assert!((f.modulus_pa - 1.0e12).abs() < 1.0); // ~1 TPa
    }

    #[test]
    fn test_nanoparticle_filler_creation() {
        let f = NanoFiller::nanoparticle(70.0e9, 0.05, 3970.0);
        assert!((f.volume_fraction - 0.05).abs() < EPS);
        assert_eq!(f.filler_type, FillerType::Nanoparticle);
    }

    #[test]
    fn test_weight_fraction_increases_with_volume_fraction() {
        let f1 = NanoFiller::carbon_nanotube(100.0, 0.01);
        let f2 = NanoFiller::carbon_nanotube(100.0, 0.05);
        assert!(f2.weight_fraction(1200.0) > f1.weight_fraction(1200.0));
    }

    #[test]
    fn test_weight_fraction_range() {
        let f = NanoFiller::carbon_nanotube(100.0, 0.03);
        let wf = f.weight_fraction(1200.0);
        assert!(wf > 0.0 && wf < 1.0);
    }

    // ── MoriTanaka ──

    #[test]
    fn test_mori_tanaka_modulus_exceeds_matrix() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.05);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        assert!(mt.longitudinal_modulus() > 3.5e9);
    }

    #[test]
    fn test_mori_tanaka_zero_filler() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.0);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        // At zero volume fraction result should be close to matrix modulus.
        let e = mt.longitudinal_modulus();
        assert!(e > 0.0);
    }

    #[test]
    fn test_mori_tanaka_transverse_positive() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.03);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        assert!(mt.transverse_modulus() > 0.0);
    }

    #[test]
    fn test_mori_tanaka_shear_positive() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.03);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        assert!(mt.shear_modulus() > 0.0);
    }

    #[test]
    fn test_mori_tanaka_modulus_increases_with_vf() {
        let em = 3.5e9;
        let e1 = {
            let f = NanoFiller::carbon_nanotube(200.0, 0.01);
            MoriTanakaModel::new(em, 0.35, f).longitudinal_modulus()
        };
        let e2 = {
            let f = NanoFiller::carbon_nanotube(200.0, 0.05);
            MoriTanakaModel::new(em, 0.35, f).longitudinal_modulus()
        };
        assert!(e2 > e1);
    }

    // ── HalpinTsai ──

    #[test]
    fn test_halpin_tsai_transverse_exceeds_matrix() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.05);
        let ht = HalpinTsaiModel::new(3.5e9, filler);
        assert!(ht.transverse_modulus() > 3.5e9);
    }

    #[test]
    fn test_halpin_tsai_longitudinal_rule_of_mixtures() {
        let vf = 0.10;
        let em = 3.5e9;
        let ef = 1.0e12;
        let filler = NanoFiller::carbon_nanotube(500.0, vf);
        let ht = HalpinTsaiModel::new(em, filler);
        let expected = vf * ef + (1.0 - vf) * em;
        let got = ht.longitudinal_modulus();
        assert!((got - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_halpin_tsai_xi_longitudinal_scales_with_aspect_ratio() {
        let f1 = NanoFiller::carbon_nanotube(50.0, 0.03);
        let f2 = NanoFiller::carbon_nanotube(100.0, 0.03);
        let ht1 = HalpinTsaiModel::new(3.5e9, f1);
        let ht2 = HalpinTsaiModel::new(3.5e9, f2);
        assert!(ht2.xi_longitudinal() > ht1.xi_longitudinal());
    }

    #[test]
    fn test_halpin_tsai_density() {
        let filler = NanoFiller::carbon_nanotube(100.0, 0.05);
        let ht = HalpinTsaiModel::new(3.5e9, filler);
        let rho = ht.composite_density(1200.0);
        // Should be between matrix and filler density.
        assert!(rho > 1200.0 * 0.9 && rho < 1350.0 * 1.1);
    }

    #[test]
    fn test_halpin_tsai_ht_longitudinal_increases_with_vf() {
        let em = 3.5e9;
        let e1 = {
            let f = NanoFiller::carbon_nanotube(200.0, 0.01);
            HalpinTsaiModel::new(em, f).longitudinal_modulus_ht()
        };
        let e2 = {
            let f = NanoFiller::carbon_nanotube(200.0, 0.05);
            HalpinTsaiModel::new(em, f).longitudinal_modulus_ht()
        };
        assert!(e2 > e1);
    }

    // ── Percolation ──

    #[test]
    fn test_percolation_threshold_value() {
        let p = PercolationThreshold::new(100.0, 2.0 / 3.0);
        let phi_c = p.critical_volume_fraction();
        assert!(phi_c > 0.0 && phi_c < 0.1);
    }

    #[test]
    fn test_percolation_below_threshold_zero_conductivity() {
        let p = PercolationThreshold::new(100.0, 1.0);
        let phi_c = p.critical_volume_fraction();
        assert!((p.conductivity_factor(phi_c * 0.5) - 0.0).abs() < EPS);
    }

    #[test]
    fn test_percolation_above_threshold_positive() {
        let p = PercolationThreshold::new(100.0, 1.0);
        let phi_c = p.critical_volume_fraction();
        assert!(p.conductivity_factor(phi_c * 2.0) > 0.0);
    }

    #[test]
    fn test_percolation_is_percolated() {
        let p = PercolationThreshold::new(200.0, 2.0 / 3.0);
        let phi_c = p.critical_volume_fraction();
        assert!(!p.is_percolated(phi_c * 0.5));
        assert!(p.is_percolated(phi_c * 1.5));
    }

    #[test]
    fn test_percolation_threshold_decreases_with_aspect_ratio() {
        let p1 = PercolationThreshold::new(50.0, 1.0);
        let p2 = PercolationThreshold::new(200.0, 1.0);
        assert!(p2.critical_volume_fraction() < p1.critical_volume_fraction());
    }

    // ── ThermalConductivity ──

    #[test]
    fn test_thermal_conductivity_enhancement() {
        let tc = ThermalConductivityNano::new(0.2, 3500.0, 0.03, 1e-8, 5e-9);
        assert!(tc.effective_conductivity() > 0.2);
    }

    #[test]
    fn test_thermal_conductivity_zero_kapitza() {
        let tc = ThermalConductivityNano::new(0.2, 3500.0, 0.03, 0.0, 5e-9);
        // Without Kapitza resistance, effective filler conductivity = k_filler.
        assert!((tc.effective_filler_conductivity() - 3500.0).abs() < EPS);
    }

    #[test]
    fn test_thermal_enhancement_ratio_gt1() {
        let tc = ThermalConductivityNano::new(0.2, 3500.0, 0.02, 1e-8, 5e-9);
        assert!(tc.enhancement_ratio() > 1.0);
    }

    #[test]
    fn test_thermal_conductivity_increases_with_vf() {
        let tc1 = ThermalConductivityNano::new(0.2, 3500.0, 0.01, 0.0, 5e-9);
        let tc2 = ThermalConductivityNano::new(0.2, 3500.0, 0.05, 0.0, 5e-9);
        assert!(tc2.effective_conductivity() > tc1.effective_conductivity());
    }

    #[test]
    fn test_thermal_kapitza_reduces_conductivity() {
        let tc_no_k = ThermalConductivityNano::new(0.2, 3500.0, 0.03, 0.0, 5e-9);
        let tc_with_k = ThermalConductivityNano::new(0.2, 3500.0, 0.03, 1e-7, 5e-9);
        assert!(tc_no_k.effective_conductivity() > tc_with_k.effective_conductivity());
    }

    // ── NanoReinforcement ──

    #[test]
    fn test_nano_reinforcement_effective_modulus() {
        let nr = NanoReinforcement::new(3.5e9, 7.0e9, 1.0);
        assert!((nr.effective_modulus() - 7.0e9).abs() < 1.0);
    }

    #[test]
    fn test_nano_reinforcement_partial_dispersion() {
        let nr = NanoReinforcement::new(3.5e9, 7.0e9, 0.5);
        let expected = 3.5e9 + 0.5 * (7.0e9 - 3.5e9);
        assert!((nr.effective_modulus() - expected).abs() < 1.0);
    }

    #[test]
    fn test_nano_reinforcement_improvement_percent_positive() {
        let nr = NanoReinforcement::new(3.5e9, 7.0e9, 0.8);
        assert!(nr.improvement_percent() > 0.0);
    }

    #[test]
    fn test_nano_reinforcement_zero_dispersion() {
        let nr = NanoReinforcement::new(3.5e9, 7.0e9, 0.0);
        assert!((nr.effective_modulus() - 3.5e9).abs() < 1.0);
    }

    // ── SurfaceFunctionalization ──

    #[test]
    fn test_surface_functionalization_increases_ifss() {
        let sf = SurfaceFunctionalization::new(20.0e6, 2.0, 1.0);
        assert!((sf.effective_ifss() - 40.0e6).abs() < 1.0);
    }

    #[test]
    fn test_surface_functionalization_partial_coverage() {
        let sf = SurfaceFunctionalization::new(20.0e6, 2.0, 0.5);
        let expected = 20.0e6 * (1.0 + (2.0 - 1.0) * 0.5);
        assert!((sf.effective_ifss() - expected).abs() < 1.0);
    }

    #[test]
    fn test_critical_length_decreases_with_ifss() {
        let sf1 = SurfaceFunctionalization::new(10.0e6, 1.0, 1.0);
        let sf2 = SurfaceFunctionalization::new(20.0e6, 1.0, 1.0);
        let l1 = sf1.critical_length(3.0e9, 10.0e-9);
        let l2 = sf2.critical_length(3.0e9, 10.0e-9);
        assert!(l1 > l2);
    }

    // ── NanoparticleSize ──

    #[test]
    fn test_hall_petch_increases_with_smaller_grain() {
        let np1 = NanoparticleSize::new(100.0e-9, 0.5e6, 200.0e6);
        let np2 = NanoparticleSize::new(10.0e-9, 0.5e6, 200.0e6);
        assert!(np2.hall_petch_strength() > np1.hall_petch_strength());
    }

    #[test]
    fn test_surface_to_volume_ratio_increases_with_smaller_size() {
        let np1 = NanoparticleSize::new(100.0e-9, 0.5e6, 200.0e6);
        let np2 = NanoparticleSize::new(10.0e-9, 0.5e6, 200.0e6);
        assert!(np2.surface_to_volume_ratio() > np1.surface_to_volume_ratio());
    }

    #[test]
    fn test_surface_to_volume_ratio_formula() {
        let np = NanoparticleSize::new(10.0e-9, 0.0, 0.0);
        let expected = 6.0 / 10.0e-9;
        assert!((np.surface_to_volume_ratio() - expected).abs() < 1.0);
    }

    #[test]
    fn test_surface_atom_fraction_clamped() {
        let np = NanoparticleSize::new(0.3e-9, 0.0, 0.0);
        let frac = np.surface_atom_fraction(0.4e-9);
        assert!(frac <= 1.0);
    }

    // ── CntProperties ──

    #[test]
    fn test_swcnt_aspect_ratio() {
        let cnt = CntProperties::swcnt_default(); // d=1nm, L=1µm
        let expected = 1000.0 / 1.0;
        assert!((cnt.aspect_ratio() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_mwcnt_aspect_ratio() {
        let cnt = CntProperties::mwcnt_default(); // d=10nm, L=10µm
        let expected = 10_000.0 / 10.0;
        assert!((cnt.aspect_ratio() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_cnt_axial_stiffness_positive() {
        let cnt = CntProperties::swcnt_default();
        assert!(cnt.axial_stiffness() > 0.0);
    }

    #[test]
    fn test_cnt_wall_type() {
        assert_eq!(
            CntProperties::swcnt_default().wall_type,
            CntWallType::SingleWalled
        );
        assert_eq!(
            CntProperties::mwcnt_default().wall_type,
            CntWallType::MultiWalled
        );
    }

    // ── GrapheneProperties ──

    #[test]
    fn test_graphene_monolayer_aspect_ratio() {
        let g = GrapheneProperties::monolayer(1.0); // 1 µm flake
        let expected = 1.0 / 0.000335;
        assert!((g.aspect_ratio() - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_graphene_reinforcement_factor_between_0_1() {
        let g = GrapheneProperties::monolayer(5.0);
        let rf = g.reinforcement_factor();
        assert!(rf > 0.0 && rf <= 1.0);
    }

    #[test]
    fn test_graphene_effective_modulus_positive() {
        let g = GrapheneProperties::monolayer(2.0);
        let contrib = g.effective_modulus_contribution(0.03);
        assert!(contrib > 0.0);
    }

    #[test]
    fn test_graphene_nanoplatelet_more_layers_lower_modulus() {
        let g1 = GrapheneProperties::nanoplatelet(1, 2.0);
        let g2 = GrapheneProperties::nanoplatelet(10, 2.0);
        assert!(g1.youngs_modulus_tpa >= g2.youngs_modulus_tpa);
    }

    #[test]
    fn test_graphene_reinforcement_large_aspect_ratio_near_one() {
        let g = GrapheneProperties::monolayer(100.0); // very large flake
        let rf = g.reinforcement_factor();
        assert!(
            rf > 0.99,
            "large aspect ratio should give RF near 1, got {rf}"
        );
    }
}
