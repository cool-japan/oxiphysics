// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nanocomposite and nanomaterial models.
//!
//! Implements homogenization theories (Mori-Tanaka, Halpin-Tsai), percolation,
//! thermal conductivity, CNT/graphene property models, and reinforcement
//! enhancement for polymer-matrix nanocomposites.
//!
//! The Mori-Tanaka model is built on three public building blocks:
//! [`mandel`] (6×6 tensor algebra), [`eshelby_spheroid`] (Eshelby tensor of
//! prolate/oblate spheroids) and [`mori_tanaka`] (Benveniste's Mori-Tanaka
//! scheme and transversely isotropic constants).

use std::f64::consts::PI;

pub mod eshelby_spheroid;
pub mod mandel;
pub mod mori_tanaka;

use self::eshelby_spheroid::spheroid_eshelby_mandel;
use self::mandel::{Mat6, isotropic_stiffness};
use self::mori_tanaka::{
    TransverselyIsotropicModuli, mori_tanaka_stiffness, transversely_isotropic_moduli,
    transversely_isotropic_stiffness,
};
use crate::error::{Error, Result};

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

impl FillerType {
    /// Whether the filler is a platelet (oblate shape whose `aspect_ratio` is
    /// lateral size / thickness) rather than a fibre or particle (prolate or
    /// spherical shape whose `aspect_ratio` is length / diameter).
    fn is_platelet(self) -> bool {
        matches!(self, FillerType::Graphene | FillerType::NanoplateletClay)
    }
}

/// A nanofiller dispersed in a matrix material.
///
/// Encapsulates aspect ratio, volume fraction, and intrinsic modulus.
#[derive(Debug, Clone)]
pub struct NanoFiller {
    /// Type of filler.
    pub filler_type: FillerType,
    /// Aspect ratio α = length/diameter (1.0 for spheres, >>1 for tubes/fibers);
    /// lateral size / thickness for platelet fillers.
    pub aspect_ratio: f64,
    /// Volume fraction (0–1).
    pub volume_fraction: f64,
    /// Longitudinal Young's modulus of the filler \[Pa\]: along the fibre axis,
    /// or in the plane of a platelet.
    pub modulus_pa: f64,
    /// Transverse modulus \[Pa\] (= modulus_pa for isotropic fillers): perpendicular
    /// to the fibre axis, or through the thickness of a platelet.
    pub transverse_modulus_pa: f64,
    /// Filler density \[kg/m³\].
    pub density: f64,
    /// Poisson's ratio ν of the filler.  For an isotropic filler this is the usual
    /// Poisson's ratio; for a transversely isotropic filler it is the contraction
    /// ratio for loading along the stiffer principal direction (see
    /// [`mori_tanaka::transversely_isotropic_stiffness`]).
    pub poisson_ratio: f64,
    /// Shear modulus G \[Pa\] of the planes that contain the filler's axis of
    /// revolution (the fibre axis, or the platelet normal).  `None` selects the
    /// isotropic value E/(2(1 + ν)), which is only defined for an isotropic filler
    /// (`modulus_pa == transverse_modulus_pa`); an anisotropic filler must set it.
    pub shear_modulus_pa: Option<f64>,
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
            poisson_ratio: 0.28, // Lu, Phys. Rev. Lett. 79, 1297 (1997): E ≈ 1 TPa, ν ≈ 0.28
            shear_modulus_pa: None,
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
            poisson_ratio: 0.165, // in-plane graphite/graphene value (Blakslee et al., 1970)
            shear_modulus_pa: None,
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
            poisson_ratio: 0.25, // typical ceramic nanoparticle; override with `with_poisson_ratio`
            shear_modulus_pa: None,
        }
    }

    /// Returns the filler with the given Poisson's ratio.
    pub fn with_poisson_ratio(mut self, poisson_ratio: f64) -> Self {
        self.poisson_ratio = poisson_ratio;
        self
    }

    /// Returns the filler with the given shear modulus \[Pa\] of the planes
    /// containing its axis of revolution (needed for anisotropic fillers).
    pub fn with_shear_modulus(mut self, shear_modulus_pa: f64) -> Self {
        self.shear_modulus_pa = Some(shear_modulus_pa);
        self
    }

    /// Mandel stiffness of the filler in the frame whose direction 1 is the
    /// axis of revolution of the inclusion shape (the fibre axis, or the
    /// platelet normal).
    ///
    /// For fibres and particles `modulus_pa` acts along that axis and
    /// `transverse_modulus_pa` perpendicular to it; for platelets the roles are
    /// exchanged (in-plane modulus = `modulus_pa`).
    fn mandel_stiffness(&self) -> Result<Mat6> {
        let (e_axis, e_perp) = if self.filler_type.is_platelet() {
            (self.transverse_modulus_pa, self.modulus_pa)
        } else {
            (self.modulus_pa, self.transverse_modulus_pa)
        };
        let finite = e_axis.is_finite() && e_perp.is_finite();
        let isotropic =
            finite && (e_axis - e_perp).abs() <= 1.0e-9 * e_axis.abs().max(e_perp.abs());
        let g_axial = match self.shear_modulus_pa {
            Some(g) => g,
            None if isotropic => e_axis / (2.0 * (1.0 + self.poisson_ratio)),
            None if finite => return Err(Error::MissingProperty("shear_modulus_pa")),
            // Non-finite moduli are reported by the stiffness builder below.
            None => f64::NAN,
        };
        transversely_isotropic_stiffness(e_axis, e_perp, self.poisson_ratio, g_axial)
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

/// Mori-Tanaka mean-field homogenization for aligned spheroidal fillers.
///
/// Implements Benveniste's form of the Mori-Tanaka scheme in 6×6 Mandel
/// notation: the matrix is isotropic (`matrix_modulus`, `matrix_poisson`), the
/// filler is transversely isotropic about its axis of revolution (see
/// [`NanoFiller`]) and all fillers are aligned.  The inclusion shape is a
/// spheroid with the exact Eshelby tensor of
/// [`eshelby_spheroid::spheroid_eshelby_mandel`]:
///
/// * carbon nanotubes and short fibres are prolate spheroids with axis ratio
///   `aspect_ratio` (length/diameter);
/// * nanoparticles are spheres when `aspect_ratio = 1` (and prolate for a larger
///   value);
/// * graphene and clay platelets are oblate spheroids of axis ratio
///   `1/aspect_ratio` (thickness/lateral size), the platelet normal being the
///   axis of revolution.
///
/// The scheme reproduces the matrix at `v_f = 0` and the filler at `v_f = 1`,
/// stays between the Voigt and Reuss bounds, and is accurate from the dilute
/// limit up to moderate loadings (roughly `v_f ≲ 0.3`).
///
/// The directional getters follow the reinforcement direction of the filler:
///
/// | getter | fibres, particles | platelets |
/// |---|---|---|
/// | [`MoriTanakaModel::longitudinal_modulus`] | along the fibre axis | in the platelet plane |
/// | [`MoriTanakaModel::transverse_modulus`] | perpendicular to the fibre | through the thickness |
/// | [`MoriTanakaModel::shear_modulus`] | shear of the planes containing the fibre axis (G₁₂) | in-plane shear |
///
/// Getters return `NaN` when the model parameters are not physical (use
/// [`MoriTanakaModel::effective_moduli`] for the reason).
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

    /// Axis ratio (axis of revolution / equatorial radius) of the equivalent
    /// spheroid: `aspect_ratio` (at least 1) for fibres and particles, its
    /// reciprocal for platelets.
    fn spheroid_axis_ratio(&self) -> Result<f64> {
        let aspect = self.filler.aspect_ratio;
        if aspect.is_nan() {
            return Err(Error::ParameterOutOfRange {
                param: "aspect_ratio",
                value: aspect,
                min: 1.0,
                max: f64::INFINITY,
            });
        }
        let slenderness = aspect.max(1.0);
        Ok(if self.filler.filler_type.is_platelet() {
            1.0 / slenderness
        } else {
            slenderness
        })
    }

    /// Effective stiffness of the composite as a 6×6 Mandel matrix
    /// (`[11, 22, 33, 23, 13, 12]`, tensor shear scaled by √2) in the frame whose
    /// direction 1 is the filler's axis of revolution (the fibre axis, or the
    /// platelet normal).  The composite is transversely isotropic about that axis.
    ///
    /// # Errors
    /// [`Error::ParameterOutOfRange`] for a non-physical matrix modulus,
    /// matrix Poisson's ratio (`-1 < ν < 1/2` is required), volume fraction
    /// (`0 ≤ v_f ≤ 1`) or filler constants; [`Error::MissingProperty`] if an
    /// anisotropic filler has no `shear_modulus_pa`; [`Error::Singularity`] if
    /// the homogenization operators cannot be inverted.
    pub fn effective_stiffness(&self) -> Result<Mat6> {
        let e_matrix = self.matrix_modulus;
        let nu_matrix = self.matrix_poisson;
        if !(e_matrix.is_finite() && e_matrix > 0.0) {
            return Err(Error::ParameterOutOfRange {
                param: "matrix_modulus",
                value: e_matrix,
                min: 0.0,
                max: f64::INFINITY,
            });
        }
        if !(nu_matrix.is_finite() && nu_matrix > -1.0 && nu_matrix < 0.5) {
            return Err(Error::ParameterOutOfRange {
                param: "matrix_poisson",
                value: nu_matrix,
                min: -1.0,
                max: 0.5,
            });
        }
        let c_matrix = isotropic_stiffness(e_matrix, nu_matrix);
        let c_filler = self.filler.mandel_stiffness()?;
        let eshelby = spheroid_eshelby_mandel(nu_matrix, self.spheroid_axis_ratio()?)?;
        mori_tanaka_stiffness(&c_matrix, &c_filler, &eshelby, self.filler.volume_fraction)
    }

    /// Effective engineering constants of the composite, referred to the
    /// filler's axis of revolution (see [`TransverselyIsotropicModuli`]).
    ///
    /// # Errors
    /// As for [`MoriTanakaModel::effective_stiffness`].
    pub fn effective_moduli(&self) -> Result<TransverselyIsotropicModuli> {
        transversely_isotropic_moduli(&self.effective_stiffness()?)
    }

    /// (longitudinal modulus, transverse modulus, shear modulus) following the
    /// reinforcement direction of the filler.
    fn directional_moduli(&self) -> Result<(f64, f64, f64)> {
        let m = self.effective_moduli()?;
        Ok(if self.filler.filler_type.is_platelet() {
            (m.e_perp, m.e_axis, m.g_perp)
        } else {
            (m.e_axis, m.e_perp, m.g_axial)
        })
    }

    /// Longitudinal composite modulus E₁ \[Pa\]: along the fibre axis (in the
    /// platelet plane for platelet fillers).  `NaN` for non-physical inputs.
    pub fn longitudinal_modulus(&self) -> f64 {
        self.directional_moduli().map_or(f64::NAN, |d| d.0)
    }

    /// Transverse composite modulus E₂ \[Pa\]: perpendicular to the fibre axis
    /// (through the thickness for platelet fillers).  `NaN` for non-physical
    /// inputs.
    pub fn transverse_modulus(&self) -> f64 {
        self.directional_moduli().map_or(f64::NAN, |d| d.1)
    }

    /// Composite shear modulus G₁₂ \[Pa\]: shear of the planes containing the
    /// fibre axis (in-plane shear for platelet fillers).  `NaN` for
    /// non-physical inputs.
    pub fn shear_modulus(&self) -> f64 {
        self.directional_moduli().map_or(f64::NAN, |d| d.2)
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
        // At zero volume fraction the composite is the matrix: no empirical correction
        // factors (the former 1 + 0.01 nu and 1.05 factors broke this limit).
        let e = mt.longitudinal_modulus();
        assert!((e - 3.5e9).abs() / 3.5e9 < 1e-12, "E1(vf = 0) = {e}");
    }

    #[test]
    fn test_mori_tanaka_matrix_limit_is_exact_for_all_moduli() {
        let nu_m = 0.35;
        let gm = 3.5e9 / (2.0 * (1.0 + nu_m));
        for filler in [
            NanoFiller::carbon_nanotube(100.0, 0.0),
            NanoFiller::graphene(500.0, 0.0),
            NanoFiller::nanoparticle(70.0e9, 0.0, 3970.0),
        ] {
            let mt = MoriTanakaModel::new(3.5e9, nu_m, filler);
            assert!((mt.longitudinal_modulus() - 3.5e9).abs() / 3.5e9 < 1e-12);
            assert!((mt.transverse_modulus() - 3.5e9).abs() / 3.5e9 < 1e-12);
            assert!((mt.shear_modulus() - gm).abs() / gm < 1e-12);
        }
    }

    #[test]
    fn test_mori_tanaka_full_filler_limit_recovers_the_filler() {
        let nu_f = 0.28;
        let filler = NanoFiller::carbon_nanotube(100.0, 1.0);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        let gf = 1.0e12 / (2.0 * (1.0 + nu_f));
        assert!((mt.longitudinal_modulus() - 1.0e12).abs() / 1.0e12 < 1e-9);
        assert!((mt.transverse_modulus() - 1.0e12).abs() / 1.0e12 < 1e-9);
        assert!((mt.shear_modulus() - gf).abs() / gf < 1e-9);
    }

    #[test]
    fn test_mori_tanaka_moduli_stay_between_reuss_and_voigt() {
        let (em, nu_m) = (3.5e9, 0.35);
        for &alpha in &[10.0, 100.0, 1000.0] {
            for &vf in &[0.005, 0.02, 0.05, 0.1] {
                let filler = NanoFiller::carbon_nanotube(alpha, vf);
                let (ef, et) = (filler.modulus_pa, filler.transverse_modulus_pa);
                let mt = MoriTanakaModel::new(em, nu_m, filler);
                let voigt = vf * ef + (1.0 - vf) * em;
                let reuss = 1.0 / (vf / et + (1.0 - vf) / em);
                let e1 = mt.longitudinal_modulus();
                let e2 = mt.transverse_modulus();
                assert!(
                    e1 > em && e1 <= voigt,
                    "α = {alpha}, vf = {vf}: E1 = {e1}, Voigt = {voigt}"
                );
                assert!(
                    e2 >= reuss && e2 > em,
                    "α = {alpha}, vf = {vf}: E2 = {e2}, Reuss = {reuss}"
                );
            }
        }
    }

    #[test]
    fn test_mori_tanaka_modulus_increases_with_aspect_ratio() {
        let e1: Vec<f64> = [3.0, 10.0, 30.0, 100.0, 300.0]
            .iter()
            .map(|&alpha| {
                let filler = NanoFiller::carbon_nanotube(alpha, 0.03);
                MoriTanakaModel::new(3.5e9, 0.35, filler).longitudinal_modulus()
            })
            .collect();
        for w in e1.windows(2) {
            assert!(w[1] > w[0], "E1 must grow with the aspect ratio: {e1:?}");
        }
    }

    #[test]
    fn test_mori_tanaka_carbon_nanotube_matches_independent_reference() {
        // Reference: tensor-algebra Mori-Tanaka with a numerically integrated Eshelby tensor.
        let filler = NanoFiller::carbon_nanotube(100.0, 0.05).with_poisson_ratio(0.2);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        let m = mt.effective_moduli().unwrap();
        assert!((mt.longitudinal_modulus() - 42249265665.269).abs() / 42249265665.269 < 1e-9);
        assert!((mt.transverse_modulus() - 4260790429.5468).abs() / 4260790429.5468 < 1e-9);
        assert!((mt.shear_modulus() - 1431957505.0198).abs() / 1431957505.0198 < 1e-9);
        assert!((m.nu_axis_perp - 0.34287141917260).abs() < 1e-9);
        assert_eq!(m.e_axis, mt.longitudinal_modulus());
    }

    #[test]
    fn test_mori_tanaka_spherical_nanoparticles_are_isotropic() {
        // Hashin-Shtrikman/Mori-Tanaka closed form for spheres.
        let (em, nu_m) = (3.5e9, 0.35);
        let (ef, nu_f) = (70.0e9, 0.25);
        let vf = 0.1;
        let filler = NanoFiller::nanoparticle(ef, vf, 3970.0);
        let mt = MoriTanakaModel::new(em, nu_m, filler);
        let km = em / (3.0 * (1.0 - 2.0 * nu_m));
        let gm = em / (2.0 * (1.0 + nu_m));
        let kf = ef / (3.0 * (1.0 - 2.0 * nu_f));
        let gf = ef / (2.0 * (1.0 + nu_f));
        let k = km + vf / (1.0 / (kf - km) + 3.0 * (1.0 - vf) / (3.0 * km + 4.0 * gm));
        let g = gm
            + vf / (1.0 / (gf - gm)
                + 6.0 * (1.0 - vf) * (km + 2.0 * gm) / (5.0 * gm * (3.0 * km + 4.0 * gm)));
        let e = 9.0 * k * g / (3.0 * k + g);
        assert!((mt.longitudinal_modulus() - e).abs() / e < 1e-11);
        assert!((mt.transverse_modulus() - e).abs() / e < 1e-11);
        assert!((mt.shear_modulus() - g).abs() / g < 1e-11);
    }

    #[test]
    fn test_mori_tanaka_platelets_are_stiff_in_plane_and_soft_through_thickness() {
        let filler = NanoFiller::graphene(100.0, 0.02).with_poisson_ratio(0.2);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, filler);
        let in_plane = mt.longitudinal_modulus();
        let through_thickness = mt.transverse_modulus();
        // Reference values from the independent implementation (oblate spheroid, axis ratio 1/100).
        assert!((in_plane - 7460243731.7955).abs() / 7460243731.7955 < 1e-9);
        assert!((through_thickness - 4295082226.9394).abs() / 4295082226.9394 < 1e-9);
        assert!((mt.shear_modulus() - 2984072789.2369).abs() / 2984072789.2369 < 1e-9);
        assert!(in_plane > through_thickness && through_thickness > 3.5e9);
        // Bounds for the in-plane modulus.
        let voigt = 0.02 * 1.0e12 + 0.98 * 3.5e9;
        let reuss = 1.0 / (0.02 / 1.0e12 + 0.98 / 3.5e9);
        assert!(in_plane < voigt && through_thickness > reuss);
    }

    #[test]
    fn test_mori_tanaka_anisotropic_filler_needs_an_explicit_shear_modulus() {
        let mut fibre = NanoFiller::carbon_nanotube(50.0, 0.3);
        fibre.modulus_pa = 230.0e9;
        fibre.transverse_modulus_pa = 15.0e9;
        let fibre = fibre.with_poisson_ratio(0.2);
        let mt = MoriTanakaModel::new(3.5e9, 0.35, fibre.clone());
        assert!(matches!(
            mt.effective_moduli(),
            Err(Error::MissingProperty("shear_modulus_pa"))
        ));
        assert!(mt.longitudinal_modulus().is_nan());
        let mt = MoriTanakaModel::new(3.5e9, 0.35, fibre.with_shear_modulus(15.0e9));
        let m = mt.effective_moduli().unwrap();
        assert!((m.e_axis - 62202830144.875).abs() / 62202830144.875 < 1e-9);
        assert!((m.e_perp - 5446135641.9319).abs() / 5446135641.9319 < 1e-9);
        assert!((m.g_axial - 2172645476.3118).abs() / 2172645476.3118 < 1e-9);
        assert!((mt.shear_modulus() - m.g_axial).abs() < 1e-3);
    }

    #[test]
    fn test_mori_tanaka_invalid_parameters_are_reported() {
        let ok = NanoFiller::carbon_nanotube(100.0, 0.03);
        let bad_fraction = NanoFiller::carbon_nanotube(100.0, 1.5);
        let negative_fraction = NanoFiller::carbon_nanotube(100.0, -0.1);
        let nan_aspect = NanoFiller::carbon_nanotube(f64::NAN, 0.03);
        for mt in [
            MoriTanakaModel::new(3.5e9, 0.35, bad_fraction),
            MoriTanakaModel::new(3.5e9, 0.35, negative_fraction),
            MoriTanakaModel::new(3.5e9, 0.35, nan_aspect),
            MoriTanakaModel::new(-3.5e9, 0.35, ok.clone()),
            MoriTanakaModel::new(3.5e9, 0.5, ok.clone()),
            MoriTanakaModel::new(3.5e9, f64::NAN, ok.clone()),
            MoriTanakaModel::new(3.5e9, 0.35, ok.clone().with_poisson_ratio(0.9)),
        ] {
            assert!(
                matches!(
                    mt.effective_stiffness(),
                    Err(Error::ParameterOutOfRange { .. })
                ),
                "expected a parameter error for {mt:?}"
            );
            assert!(mt.longitudinal_modulus().is_nan());
            assert!(mt.transverse_modulus().is_nan());
            assert!(mt.shear_modulus().is_nan());
        }
        // A valid model is finite.
        let mt = MoriTanakaModel::new(3.5e9, 0.35, ok);
        assert!(mt.longitudinal_modulus().is_finite());
    }

    #[test]
    fn test_mori_tanaka_continuous_fibre_limit_matches_hill() {
        // Infinite aspect ratio: continuous aligned fibres.
        let (em, nu_m) = (3.5e9, 0.35);
        let vf = 0.3;
        let mut filler = NanoFiller::carbon_nanotube(f64::INFINITY, vf);
        filler.modulus_pa = 230.0e9;
        filler.transverse_modulus_pa = 230.0e9;
        let filler = filler.with_poisson_ratio(0.2);
        let mt = MoriTanakaModel::new(em, nu_m, filler);
        let (ef, nu_f) = (230.0e9, 0.2);
        let gm = em / (2.0 * (1.0 + nu_m));
        let gf = ef / (2.0 * (1.0 + nu_f));
        let km = em / (3.0 * (1.0 - 2.0 * nu_m)) + gm / 3.0;
        let kf = ef / (3.0 * (1.0 - 2.0 * nu_f)) + gf / 3.0;
        let hill = vf * ef
            + (1.0 - vf) * em
            + 4.0 * (nu_f - nu_m).powi(2) * vf * (1.0 - vf)
                / (vf / km + (1.0 - vf) / kf + 1.0 / gm);
        assert!((mt.longitudinal_modulus() - hill).abs() / hill < 1e-9);
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
