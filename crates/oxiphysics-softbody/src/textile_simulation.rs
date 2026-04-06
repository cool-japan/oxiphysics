// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Textile simulation: yarn mechanics, woven/knitted fabric, collision, and garment draping.
//!
//! # Overview
//!
//! This module provides a physically-based framework for simulating textile materials:
//!
//! * [`YarnModel`] — Twisted fiber bundle cross-section, yarn-yarn contact.
//! * [`WovenFabric`] — Plain/twill/satin weave geometry, interlacing pattern, unit cell.
//! * [`FabricMechanics`] — Tensile, shear, bending stiffness from classical lamination theory.
//! * [`KnittedFabric`] — Loop interlocking, loop length, stitch pattern, extensibility.
//! * [`FabricCollision`] — Self-collision with vertex-face detection, friction, layered fabric.
//! * [`GarmentDraping`] — Gravity draping simulation, seam constraints, fitting metrics.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use oxiphysics_core::math::{Real, Vec3};

// ─── YarnModel ────────────────────────────────────────────────────────────────

/// Model of a single yarn as a twisted bundle of fibers.
///
/// A yarn is characterised by its fiber bundle geometry (cross-section, twist),
/// elastic properties, and contact behaviour when yarns press against each other.
#[derive(Debug, Clone)]
pub struct YarnModel {
    /// Number of fibers in the bundle.
    pub num_fibers: usize,
    /// Radius of a single fiber \[m\].
    pub fiber_radius: Real,
    /// Twist angle per unit length \[rad/m\] (helix angle of the bundle).
    pub twist_per_length: Real,
    /// Young's modulus of the fiber material \[Pa\].
    pub fiber_young_modulus: Real,
    /// Shear modulus of the fiber material \[Pa\].
    pub fiber_shear_modulus: Real,
    /// Fiber packing fraction (0 < φ ≤ π/(2√3) for hexagonal close-packing).
    pub packing_fraction: Real,
    /// Yarn linear density \[kg/m\] (tex = g/km).
    pub linear_density: Real,
    /// Friction coefficient for yarn-yarn contact.
    pub friction_coefficient: Real,
}

impl YarnModel {
    /// Create a new yarn model.
    ///
    /// # Arguments
    /// * `num_fibers` — number of fibers in the twisted bundle
    /// * `fiber_radius` — individual fiber radius \[m\]
    /// * `twist_per_length` — twist angle per unit length \[rad/m\]
    /// * `young_modulus` — fiber Young's modulus \[Pa\]
    /// * `shear_modulus` — fiber shear modulus \[Pa\]
    /// * `linear_density` — yarn linear density \[kg/m\]
    /// * `friction_coefficient` — yarn-yarn friction coefficient
    pub fn new(
        num_fibers: usize,
        fiber_radius: Real,
        twist_per_length: Real,
        young_modulus: Real,
        shear_modulus: Real,
        linear_density: Real,
        friction_coefficient: Real,
    ) -> Self {
        let n = num_fibers as Real;
        let packing_fraction = (n * fiber_radius * fiber_radius * std::f64::consts::PI)
            / (yarn_cross_section_area(num_fibers, fiber_radius) * 1.1);
        Self {
            num_fibers,
            fiber_radius,
            twist_per_length,
            fiber_young_modulus: young_modulus,
            fiber_shear_modulus: shear_modulus,
            packing_fraction: packing_fraction.min(0.9),
            linear_density,
            friction_coefficient,
        }
    }

    /// Compute the effective yarn cross-section radius \[m\].
    ///
    /// Based on hexagonal packing of circular fibers.
    pub fn yarn_radius(&self) -> Real {
        let n = self.num_fibers as Real;
        self.fiber_radius * (n / self.packing_fraction).sqrt()
    }

    /// Compute the effective axial stiffness EA of the yarn \[N\].
    ///
    /// Uses the rule of mixtures: EA = φ · E_f · A_yarn
    pub fn axial_stiffness(&self) -> Real {
        let a_yarn = std::f64::consts::PI * self.yarn_radius().powi(2);
        self.packing_fraction * self.fiber_young_modulus * a_yarn
    }

    /// Compute the effective bending stiffness EI of the yarn \[N·m²\].
    ///
    /// EI = φ · E_f · I_yarn, where I_yarn = π r⁴ / 4
    pub fn bending_stiffness(&self) -> Real {
        let r = self.yarn_radius();
        let i_yarn = std::f64::consts::PI * r.powi(4) / 4.0;
        self.packing_fraction * self.fiber_young_modulus * i_yarn
    }

    /// Compute the helix angle of fibers in the twisted bundle \[rad\].
    ///
    /// θ = arctan(2π r_yarn · twist_per_length)
    pub fn helix_angle(&self) -> Real {
        (2.0 * std::f64::consts::PI * self.yarn_radius() * self.twist_per_length).atan()
    }

    /// Compute the yarn-to-yarn normal contact force using Hertzian contact theory.
    ///
    /// # Arguments
    /// * `overlap` — geometric overlap between yarn centers \[m\] (positive = penetrating)
    /// * `crossing_angle` — angle between the two yarn axes \[rad\]
    ///
    /// Returns the normal contact force \[N\].
    pub fn contact_force(&self, overlap: Real, crossing_angle: Real) -> Real {
        if overlap <= 0.0 {
            return 0.0;
        }
        // Hertz contact for two cylinders crossed at angle α:
        // F = k_contact * overlap^(3/2)
        let r = self.yarn_radius();
        let e_star = self.fiber_young_modulus / (2.0 * (1.0 - 0.3_f64.powi(2)));
        let r_eff = r / (2.0 * crossing_angle.sin().max(1e-6));
        let k_contact = (4.0 / 3.0) * e_star * r_eff.sqrt();
        k_contact * overlap.powf(1.5)
    }

    /// Compute friction force from yarn-yarn contact.
    ///
    /// Uses Coulomb's law: F_friction = μ · F_normal
    pub fn friction_force(&self, normal_force: Real) -> Real {
        self.friction_coefficient * normal_force.abs()
    }

    /// Estimate yarn tensile strength \[N\] based on fiber properties.
    ///
    /// Assumes all fibers share load proportional to the cosine of the helix angle.
    pub fn tensile_strength(&self, fiber_tensile_strength: Real) -> Real {
        let theta = self.helix_angle();
        let n = self.num_fibers as Real;
        n * fiber_tensile_strength * std::f64::consts::PI * self.fiber_radius.powi(2) * theta.cos()
    }

    /// Compute the strain energy density in the yarn for a given axial strain ε.
    pub fn strain_energy(&self, axial_strain: Real) -> Real {
        0.5 * self.axial_stiffness() * axial_strain.powi(2)
    }
}

/// Compute the approximate yarn cross-section area from fiber packing.
fn yarn_cross_section_area(num_fibers: usize, fiber_radius: Real) -> Real {
    // Approximate circular cross-section containing N hexagonally packed fibers
    let n = num_fibers as Real;
    std::f64::consts::PI * (fiber_radius * (n.sqrt() + 1.0)).powi(2)
}

// ─── WeavePattern ─────────────────────────────────────────────────────────────

/// Weave pattern types for woven fabric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeavePattern {
    /// Plain weave: 1-over-1-under interlacing.
    Plain,
    /// Twill weave: diagonal interlacing pattern.
    Twill(usize),
    /// Satin weave: long floats for smooth surface.
    Satin(usize),
    /// Basket weave: 2-over-2-under.
    Basket,
}

impl WeavePattern {
    /// Return the repeat size (unit cell width) in warp picks.
    pub fn repeat_size(&self) -> usize {
        match self {
            WeavePattern::Plain => 2,
            WeavePattern::Twill(n) => *n,
            WeavePattern::Satin(n) => *n,
            WeavePattern::Basket => 4,
        }
    }

    /// Compute whether the warp yarn is on top at position (i, j) in the weave matrix.
    ///
    /// Returns `true` if warp is over weft at the given intersection.
    pub fn warp_over_weft(&self, i: usize, j: usize) -> bool {
        match self {
            WeavePattern::Plain => (i + j).is_multiple_of(2),
            WeavePattern::Twill(n) => (i + j) % n < n / 2,
            WeavePattern::Satin(n) => {
                let shift = 2; // typical satin shift
                (i * shift) % n == j % n
            }
            WeavePattern::Basket => ((i / 2) + (j / 2)).is_multiple_of(2),
        }
    }

    /// Count the number of interlacings per unit cell.
    pub fn interlacings_per_unit_cell(&self) -> usize {
        let n = self.repeat_size();
        match self {
            WeavePattern::Plain => 2,
            WeavePattern::Twill(_) => n / 2,
            WeavePattern::Satin(_) => 1,
            WeavePattern::Basket => 2,
        }
    }
}

// ─── WovenFabric ─────────────────────────────────────────────────────────────

/// A woven fabric defined by its yarn properties and weave geometry.
///
/// Woven fabrics consist of interlaced warp (lengthwise) and weft (crosswise)
/// yarns. The unit cell describes one repeat of the weave pattern.
#[derive(Debug, Clone)]
pub struct WovenFabric {
    /// Warp yarn model (lengthwise yarns).
    pub warp_yarn: YarnModel,
    /// Weft yarn model (crosswise yarns).
    pub weft_yarn: YarnModel,
    /// Weave pattern (plain, twill, satin, etc.).
    pub pattern: WeavePattern,
    /// Warp yarn spacing (ends per unit length) \[1/m\].
    pub warp_spacing: Real,
    /// Weft yarn spacing (picks per unit length) \[1/m\].
    pub weft_spacing: Real,
    /// Fabric thickness \[m\].
    pub thickness: Real,
    /// Areal density \[kg/m²\].
    pub areal_density: Real,
}

impl WovenFabric {
    /// Create a new woven fabric definition.
    pub fn new(
        warp_yarn: YarnModel,
        weft_yarn: YarnModel,
        pattern: WeavePattern,
        warp_spacing: Real,
        weft_spacing: Real,
    ) -> Self {
        let thickness = 2.0 * (warp_yarn.yarn_radius() + weft_yarn.yarn_radius());
        let areal_density =
            warp_yarn.linear_density * warp_spacing + weft_yarn.linear_density * weft_spacing;
        Self {
            warp_yarn,
            weft_yarn,
            pattern,
            warp_spacing,
            weft_spacing,
            thickness,
            areal_density,
        }
    }

    /// Compute the unit cell dimensions (warp repeat × weft repeat) \[m × m\].
    pub fn unit_cell_size(&self) -> (Real, Real) {
        let n = self.pattern.repeat_size();
        (n as Real / self.warp_spacing, n as Real / self.weft_spacing)
    }

    /// Compute the weave crimp of warp yarns as a fraction (0..1).
    ///
    /// Crimp = (yarn path length - fabric length) / fabric length
    pub fn warp_crimp(&self) -> Real {
        let d_weft = 2.0 * self.weft_yarn.yarn_radius();
        let spacing = 1.0 / self.weft_spacing;
        let interlacings = self.pattern.interlacings_per_unit_cell() as Real;
        // Approximate: crimp ≈ π * d_weft / (2 * spacing) * interlacings
        std::f64::consts::PI * d_weft * interlacings
            / (2.0 * spacing * self.pattern.repeat_size() as Real)
    }

    /// Compute the weave crimp of weft yarns.
    pub fn weft_crimp(&self) -> Real {
        let d_warp = 2.0 * self.warp_yarn.yarn_radius();
        let spacing = 1.0 / self.warp_spacing;
        let interlacings = self.pattern.interlacings_per_unit_cell() as Real;
        std::f64::consts::PI * d_warp * interlacings
            / (2.0 * spacing * self.pattern.repeat_size() as Real)
    }

    /// Generate the interlacing pattern matrix for an n×n unit cell.
    ///
    /// Returns a 2D array of booleans: `true` = warp over weft.
    pub fn interlacing_matrix(&self, n: usize) -> Vec<Vec<bool>> {
        (0..n)
            .map(|i| (0..n).map(|j| self.pattern.warp_over_weft(i, j)).collect())
            .collect()
    }

    /// Compute the cover factor (fraction of fabric area covered by yarns).
    ///
    /// CF = d_warp * e_warp + d_weft * e_weft - d_warp * d_weft * e_warp * e_weft
    pub fn cover_factor(&self) -> Real {
        let d_warp = 2.0 * self.warp_yarn.yarn_radius();
        let d_weft = 2.0 * self.weft_yarn.yarn_radius();
        let e_warp = self.warp_spacing;
        let e_weft = self.weft_spacing;
        let cf = d_warp * e_warp + d_weft * e_weft - d_warp * d_weft * e_warp * e_weft;
        cf.min(1.0)
    }

    /// Compute the porosity of the fabric (fraction of open area).
    pub fn porosity(&self) -> Real {
        (1.0 - self.cover_factor()).max(0.0)
    }

    /// Estimate the tearing strength ratio (satin > twill > plain for equal yarns).
    pub fn relative_tearing_strength(&self) -> Real {
        match self.pattern {
            WeavePattern::Plain => 1.0,
            WeavePattern::Twill(n) => 1.0 + 0.1 * (n as Real - 2.0),
            WeavePattern::Satin(n) => 1.0 + 0.15 * (n as Real - 2.0),
            WeavePattern::Basket => 1.1,
        }
    }
}

// ─── FabricMechanics ─────────────────────────────────────────────────────────

/// Mechanical stiffness properties of a fabric computed from yarn-level properties.
///
/// Based on classical lamination theory adapted for textile composites.
#[derive(Debug, Clone)]
pub struct FabricMechanics {
    /// Tensile stiffness in warp direction \[N/m\].
    pub tensile_stiffness_warp: Real,
    /// Tensile stiffness in weft direction \[N/m\].
    pub tensile_stiffness_weft: Real,
    /// In-plane shear stiffness \[N/m\].
    pub shear_stiffness: Real,
    /// Bending stiffness in warp direction \[N·m\].
    pub bending_stiffness_warp: Real,
    /// Bending stiffness in weft direction \[N·m\].
    pub bending_stiffness_weft: Real,
    /// Poisson's ratio (warp-to-weft lateral contraction).
    pub poisson_ratio: Real,
    /// Fabric thickness \[m\].
    pub thickness: Real,
}

impl FabricMechanics {
    /// Compute fabric mechanics from a woven fabric definition.
    pub fn from_woven_fabric(fabric: &WovenFabric) -> Self {
        let t = fabric.thickness;

        // Tensile stiffness from yarn axial stiffness per unit width
        // K_t = EA_yarn * yarn_spacing (N/m per unit fabric width)
        let k_warp =
            fabric.warp_yarn.axial_stiffness() * fabric.warp_spacing * (1.0 - fabric.warp_crimp());
        let k_weft =
            fabric.weft_yarn.axial_stiffness() * fabric.weft_spacing * (1.0 - fabric.weft_crimp());

        // Shear stiffness: dominated by yarn-yarn friction and interlacing geometry
        let shear = 0.5
            * (k_warp * k_weft).sqrt()
            * fabric.pattern.interlacings_per_unit_cell() as Real
            * fabric.warp_yarn.friction_coefficient;

        // Bending stiffness from yarn bending stiffness EI per unit width
        let b_warp = fabric.warp_yarn.bending_stiffness() * fabric.warp_spacing;
        let b_weft = fabric.weft_yarn.bending_stiffness() * fabric.weft_spacing;

        // Poisson's ratio: determined by crimp interchange
        let nu = fabric.warp_crimp().min(0.5) * fabric.weft_crimp().min(0.5) * 2.0;

        Self {
            tensile_stiffness_warp: k_warp,
            tensile_stiffness_weft: k_weft,
            shear_stiffness: shear,
            bending_stiffness_warp: b_warp,
            bending_stiffness_weft: b_weft,
            poisson_ratio: nu,
            thickness: t,
        }
    }

    /// Compute the extensional stiffness matrix A \[N/m\] (3×3 in Voigt notation).
    ///
    /// Returns `[A11, A22, A12, A66]` (warp, weft, coupling, shear).
    pub fn extensional_stiffness(&self) -> [Real; 4] {
        let a11 = self.tensile_stiffness_warp;
        let a22 = self.tensile_stiffness_weft;
        let a12 = self.poisson_ratio * (a11 * a22).sqrt();
        let a66 = self.shear_stiffness;
        [a11, a22, a12, a66]
    }

    /// Compute the bending stiffness matrix D \[N·m\] (3×3 in Voigt notation).
    ///
    /// Returns `[D11, D22, D12, D66]`.
    pub fn bending_stiffness_matrix(&self) -> [Real; 4] {
        let d11 = self.bending_stiffness_warp;
        let d22 = self.bending_stiffness_weft;
        let d12 = self.poisson_ratio * (d11 * d22).sqrt();
        let d66 = self.shear_stiffness * self.thickness.powi(2) / 12.0;
        [d11, d22, d12, d66]
    }

    /// Compute the in-plane effective Young's modulus in warp direction \[Pa\].
    pub fn effective_modulus_warp(&self) -> Real {
        self.tensile_stiffness_warp / self.thickness
    }

    /// Compute the in-plane effective Young's modulus in weft direction \[Pa\].
    pub fn effective_modulus_weft(&self) -> Real {
        self.tensile_stiffness_weft / self.thickness
    }

    /// Compute the in-plane effective shear modulus \[Pa\].
    pub fn effective_shear_modulus(&self) -> Real {
        self.shear_stiffness / self.thickness
    }

    /// Compute the tensile force per unit width for given strains.
    ///
    /// # Arguments
    /// * `eps_warp` — warp direction strain
    /// * `eps_weft` — weft direction strain
    ///
    /// Returns `(N_warp, N_weft)` in N/m.
    pub fn tensile_force(&self, eps_warp: Real, eps_weft: Real) -> (Real, Real) {
        let [a11, a22, a12, _] = self.extensional_stiffness();
        let n_warp = a11 * eps_warp + a12 * eps_weft;
        let n_weft = a12 * eps_warp + a22 * eps_weft;
        (n_warp, n_weft)
    }

    /// Compute the bending moment per unit width for given curvatures.
    ///
    /// # Arguments
    /// * `kappa_warp` — curvature in warp direction \[1/m\]
    /// * `kappa_weft` — curvature in weft direction \[1/m\]
    ///
    /// Returns `(M_warp, M_weft)` in N.
    pub fn bending_moment(&self, kappa_warp: Real, kappa_weft: Real) -> (Real, Real) {
        let [d11, d22, d12, _] = self.bending_stiffness_matrix();
        let m_warp = d11 * kappa_warp + d12 * kappa_weft;
        let m_weft = d12 * kappa_warp + d22 * kappa_weft;
        (m_warp, m_weft)
    }

    /// Compute the shear force per unit width for a given shear angle γ \[rad\].
    pub fn shear_force(&self, shear_angle: Real) -> Real {
        self.shear_stiffness * shear_angle
    }

    /// Compute the energy density due to in-plane deformation \[J/m²\].
    pub fn strain_energy_density(&self, eps_warp: Real, eps_weft: Real, gamma: Real) -> Real {
        let [a11, a22, a12, a66] = self.extensional_stiffness();
        0.5 * (a11 * eps_warp.powi(2)
            + a22 * eps_weft.powi(2)
            + 2.0 * a12 * eps_warp * eps_weft
            + a66 * gamma.powi(2))
    }
}

// ─── StitchPattern ────────────────────────────────────────────────────────────

/// Stitch pattern for knitted fabrics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StitchPattern {
    /// Jersey (single) knit — one loop per row per column.
    Jersey,
    /// Rib stitch — alternating knit/purl columns.
    Rib,
    /// Purl stitch — all loops pulled toward the viewer.
    Purl,
    /// Cable stitch — groups of loops twisted periodically.
    Cable(usize),
}

impl StitchPattern {
    /// Return the ratio of width to height of a single loop for this stitch.
    pub fn loop_aspect_ratio(&self) -> Real {
        match self {
            StitchPattern::Jersey => 1.0,
            StitchPattern::Rib => 0.7,
            StitchPattern::Purl => 1.2,
            StitchPattern::Cable(n) => 1.0 / (*n as Real),
        }
    }

    /// Approximate extensibility ratio in the course direction.
    pub fn course_extensibility(&self) -> Real {
        match self {
            StitchPattern::Jersey => 1.6,
            StitchPattern::Rib => 2.0,
            StitchPattern::Purl => 1.4,
            StitchPattern::Cable(_) => 1.2,
        }
    }

    /// Approximate extensibility ratio in the wale direction.
    pub fn wale_extensibility(&self) -> Real {
        match self {
            StitchPattern::Jersey => 1.4,
            StitchPattern::Rib => 1.5,
            StitchPattern::Purl => 1.6,
            StitchPattern::Cable(_) => 1.1,
        }
    }
}

// ─── KnittedFabric ───────────────────────────────────────────────────────────

/// A knitted fabric made from interlocked yarn loops.
///
/// Knitted fabrics are characterised by their stitch pattern, loop geometry,
/// and the extensibility arising from loop straightening under tension.
#[derive(Debug, Clone)]
pub struct KnittedFabric {
    /// Yarn used to form the loops.
    pub yarn: YarnModel,
    /// Stitch pattern.
    pub stitch: StitchPattern,
    /// Loop length (yarn consumed per stitch) \[m\].
    pub loop_length: Real,
    /// Course density (loops per unit length in course direction) \[1/m\].
    pub course_density: Real,
    /// Wale density (loops per unit length in wale direction) \[1/m\].
    pub wale_density: Real,
    /// Areal density \[kg/m²\].
    pub areal_density: Real,
    /// Fabric thickness \[m\].
    pub thickness: Real,
}

impl KnittedFabric {
    /// Create a new knitted fabric model.
    pub fn new(
        yarn: YarnModel,
        stitch: StitchPattern,
        loop_length: Real,
        course_density: Real,
        wale_density: Real,
    ) -> Self {
        let areal_density = yarn.linear_density * loop_length * course_density * wale_density;
        let thickness = 4.0 * yarn.yarn_radius();
        Self {
            yarn,
            stitch,
            loop_length,
            course_density,
            wale_density,
            areal_density,
            thickness,
        }
    }

    /// Compute the stitch density (loops per unit area) \[1/m²\].
    pub fn stitch_density(&self) -> Real {
        self.course_density * self.wale_density
    }

    /// Compute the fabric cover factor.
    pub fn cover_factor(&self) -> Real {
        let d_yarn = 2.0 * self.yarn.yarn_radius();
        let cf = d_yarn * self.loop_length * self.stitch_density();
        cf.min(1.0)
    }

    /// Compute the loop shape factor (width-to-height ratio of a single loop).
    pub fn loop_shape_factor(&self) -> Real {
        let loop_width = 1.0 / self.wale_density;
        let loop_height = 1.0 / self.course_density;
        loop_width / loop_height
    }

    /// Estimate the tensile stiffness in the course direction \[N/m\].
    ///
    /// Under tension, loops straighten; stiffness is governed by loop geometry
    /// and yarn bending stiffness.
    pub fn course_tensile_stiffness(&self) -> Real {
        let r = self.yarn.yarn_radius();
        let ei = self.yarn.bending_stiffness();
        // Stiffness ~ EI / l^3 * extensibility factor
        let l = self.loop_length;
        let ext = self.stitch.course_extensibility();
        ei / (l * l * l) * self.course_density * ext * r
    }

    /// Estimate the tensile stiffness in the wale direction \[N/m\].
    pub fn wale_tensile_stiffness(&self) -> Real {
        let r = self.yarn.yarn_radius();
        let ei = self.yarn.bending_stiffness();
        let l = self.loop_length;
        let ext = self.stitch.wale_extensibility();
        ei / (l * l * l) * self.wale_density * ext * r
    }

    /// Compute the maximum extension in the course direction before yarn straightening.
    pub fn max_course_extension(&self) -> Real {
        self.stitch.course_extensibility() / self.course_density
    }

    /// Compute the maximum extension in the wale direction.
    pub fn max_wale_extension(&self) -> Real {
        self.stitch.wale_extensibility() / self.wale_density
    }

    /// Compute tensile force per unit width for a given course-direction strain ε.
    ///
    /// Uses a bilinear model: low stiffness in loop-straightening regime, high after.
    pub fn course_tensile_force(&self, strain: Real) -> Real {
        let ext = self.stitch.course_extensibility();
        let k = self.course_tensile_stiffness();
        if strain < (ext - 1.0) / ext {
            // Loop straightening regime
            k * 0.1 * strain
        } else {
            // Yarn stretching regime
            k * strain
        }
    }

    /// Compute the inter-loop contact force for a given compression.
    pub fn loop_contact_force(&self, compression: Real) -> Real {
        if compression <= 0.0 {
            return 0.0;
        }
        self.yarn
            .contact_force(compression, std::f64::consts::PI / 2.0)
    }
}

// ─── FabricCollision ─────────────────────────────────────────────────────────

/// A particle in the fabric mesh used for collision detection.
#[derive(Debug, Clone)]
pub struct FabricParticle {
    /// Current position \[m\].
    pub position: Vec3,
    /// Previous position (for PBD) \[m\].
    pub prev_position: Vec3,
    /// Velocity \[m/s\].
    pub velocity: Vec3,
    /// Inverse mass \[1/kg\] (0 = pinned).
    pub inv_mass: Real,
}

impl FabricParticle {
    /// Create a new fabric particle.
    pub fn new(position: Vec3, mass: Real) -> Self {
        let inv_mass = if mass > 1e-12 { 1.0 / mass } else { 0.0 };
        Self {
            position,
            prev_position: position,
            velocity: Vec3::zeros(),
            inv_mass,
        }
    }

    /// Create a pinned (static) fabric particle.
    pub fn pinned(position: Vec3) -> Self {
        Self {
            position,
            prev_position: position,
            velocity: Vec3::zeros(),
            inv_mass: 0.0,
        }
    }

    /// Check if the particle is pinned (immovable).
    pub fn is_pinned(&self) -> bool {
        self.inv_mass < 1e-12
    }
}

/// Self-collision detection and response for a fabric mesh.
///
/// Handles vertex-face (point-triangle) collision, including friction modelling
/// and support for layered (multi-layer) fabric stacks.
#[derive(Debug, Clone)]
pub struct FabricCollision {
    /// Fabric particles.
    pub particles: Vec<FabricParticle>,
    /// Triangle indices into the particle array.
    pub triangles: Vec<[usize; 3]>,
    /// Collision thickness (proxy radius for vertex-face tests) \[m\].
    pub collision_thickness: Real,
    /// Friction coefficient for fabric self-contact.
    pub friction: Real,
    /// Number of collision constraint iterations per step.
    pub collision_iterations: usize,
    /// Restitution coefficient.
    pub restitution: Real,
}

impl FabricCollision {
    /// Create a new fabric collision system.
    pub fn new(
        particles: Vec<FabricParticle>,
        triangles: Vec<[usize; 3]>,
        collision_thickness: Real,
        friction: Real,
    ) -> Self {
        Self {
            particles,
            triangles,
            collision_thickness,
            friction,
            collision_iterations: 5,
            restitution: 0.0,
        }
    }

    /// Test if a vertex is penetrating a triangle (vertex-face test).
    ///
    /// Returns `Some(penetration_depth, contact_normal)` if penetrating,
    /// `None` otherwise.
    pub fn vertex_face_test(&self, vertex_idx: usize, tri: [usize; 3]) -> Option<(Real, Vec3)> {
        // Skip if vertex is part of the triangle
        if tri.contains(&vertex_idx) {
            return None;
        }
        let p = self.particles[vertex_idx].position;
        let a = self.particles[tri[0]].position;
        let b = self.particles[tri[1]].position;
        let c = self.particles[tri[2]].position;

        let ab = b - a;
        let ac = c - a;
        let ap = p - a;

        let normal = ab.cross(&ac);
        let n_len = normal.norm();
        if n_len < 1e-10 {
            return None;
        }
        let n_hat = normal / n_len;

        let dist = ap.dot(&n_hat);
        let thickness = self.collision_thickness;

        if dist.abs() < thickness {
            // Check if projection is inside the triangle
            if point_in_triangle_barycentric(p, a, b, c) {
                let depth = thickness - dist.abs();
                let contact_normal = if dist >= 0.0 { n_hat } else { -n_hat };
                return Some((depth, contact_normal));
            }
        }
        None
    }

    /// Resolve vertex-face penetrations using position-based correction.
    pub fn resolve_self_collisions(&mut self) {
        for _ in 0..self.collision_iterations {
            let n_verts = self.particles.len();
            let tris = self.triangles.clone();
            for vi in 0..n_verts {
                for tri in &tris {
                    if let Some((depth, normal)) = self.vertex_face_test(vi, *tri) {
                        // Position correction: push vertex out
                        let w_v = self.particles[vi].inv_mass;
                        let w_tri: Real = tri.iter().map(|&i| self.particles[i].inv_mass).sum();
                        let total_w = w_v + w_tri / 3.0;
                        if total_w < 1e-12 {
                            continue;
                        }
                        let correction = depth / total_w;
                        if !self.particles[vi].is_pinned() {
                            self.particles[vi].position += normal * (correction * w_v);
                        }
                        for &ti in tri.iter() {
                            if !self.particles[ti].is_pinned() {
                                let w_i = self.particles[ti].inv_mass;
                                self.particles[ti].position -= normal * (correction * w_i / 3.0);
                            }
                        }
                    }
                }
            }
        }
    }

    /// Apply ground plane collision (floor at y = 0).
    pub fn apply_ground_collision(&mut self, floor_y: Real) {
        for p in &mut self.particles {
            if !p.is_pinned() && p.position.y < floor_y + self.collision_thickness {
                let penetration = floor_y + self.collision_thickness - p.position.y;
                p.position.y += penetration;
                // Apply friction to horizontal velocity
                let normal = Vec3::new(0.0, 1.0, 0.0);
                let v_n = p.velocity.dot(&normal);
                let v_t = p.velocity - normal * v_n;
                p.velocity = v_t * (1.0 - self.friction).max(0.0);
                if v_n < 0.0 {
                    p.velocity.y = -self.restitution * v_n;
                }
            }
        }
    }

    /// Apply sphere collision (for garment fitting over a body proxy).
    pub fn apply_sphere_collision(&mut self, center: Vec3, radius: Real) {
        for p in &mut self.particles {
            if !p.is_pinned() {
                let dp = p.position - center;
                let dist = dp.norm();
                let target = radius + self.collision_thickness;
                if dist < target && dist > 1e-10 {
                    p.position = center + dp / dist * target;
                }
            }
        }
    }

    /// Compute the average inter-particle separation (useful for diagnostics).
    pub fn average_particle_separation(&self) -> Real {
        if self.triangles.is_empty() {
            return 0.0;
        }
        let total: Real = self
            .triangles
            .iter()
            .map(|tri| {
                let a = self.particles[tri[0]].position;
                let b = self.particles[tri[1]].position;
                let c = self.particles[tri[2]].position;
                ((a - b).norm() + (b - c).norm() + (c - a).norm()) / 3.0
            })
            .sum();
        total / self.triangles.len() as Real
    }

    /// Count the number of active (non-pinned) particles.
    pub fn active_particle_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.is_pinned()).count()
    }
}

/// Test if point `p` projected onto the plane of triangle `(a, b, c)` lies inside
/// the triangle, using barycentric coordinates.
fn point_in_triangle_barycentric(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> bool {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;

    let d00 = ab.dot(&ab);
    let d01 = ab.dot(&ac);
    let d11 = ac.dot(&ac);
    let d20 = ap.dot(&ab);
    let d21 = ap.dot(&ac);

    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-12 {
        return false;
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;

    u >= 0.0 && v >= 0.0 && w >= 0.0
}

// ─── GarmentDraping ──────────────────────────────────────────────────────────

/// A seam constraint connecting two edge particles of adjacent fabric panels.
#[derive(Debug, Clone)]
pub struct SeamConstraint {
    /// Index of particle on panel A.
    pub particle_a: usize,
    /// Index of particle on panel B.
    pub particle_b: usize,
    /// Rest length of the seam (0 = coincident seam) \[m\].
    pub rest_length: Real,
    /// Stiffness of the seam constraint \[N/m\].
    pub stiffness: Real,
}

impl SeamConstraint {
    /// Create a new seam constraint.
    pub fn new(particle_a: usize, particle_b: usize, rest_length: Real, stiffness: Real) -> Self {
        Self {
            particle_a,
            particle_b,
            rest_length,
            stiffness,
        }
    }

    /// Compute the constraint violation (excess length) \[m\].
    pub fn violation(&self, particles: &[FabricParticle]) -> Real {
        let pa = particles[self.particle_a].position;
        let pb = particles[self.particle_b].position;
        (pa - pb).norm() - self.rest_length
    }

    /// Project the seam constraint using XPBD.
    pub fn project(&self, particles: &mut [FabricParticle], dt: Real) {
        let pa = particles[self.particle_a].position;
        let pb = particles[self.particle_b].position;
        let dp = pa - pb;
        let dist = dp.norm();
        if dist < 1e-12 {
            return;
        }
        let constraint = dist - self.rest_length;
        if constraint.abs() < 1e-12 {
            return;
        }
        let alpha = 1.0 / (self.stiffness * dt * dt);
        let wa = particles[self.particle_a].inv_mass;
        let wb = particles[self.particle_b].inv_mass;
        let total_w = wa + wb + alpha;
        if total_w < 1e-12 {
            return;
        }
        let lambda = -constraint / total_w;
        let dir = dp / dist;
        if !particles[self.particle_a].is_pinned() {
            particles[self.particle_a].position += dir * (lambda * wa);
        }
        if !particles[self.particle_b].is_pinned() {
            particles[self.particle_b].position -= dir * (lambda * wb);
        }
    }
}

/// Distance constraint for fabric stretch/compression.
#[derive(Debug, Clone)]
pub struct FabricDistanceConstraint {
    /// Index of first particle.
    pub idx_a: usize,
    /// Index of second particle.
    pub idx_b: usize,
    /// Rest length \[m\].
    pub rest_length: Real,
    /// Compliance (inverse stiffness) \[m/N\].
    pub compliance: Real,
}

impl FabricDistanceConstraint {
    /// Create a new distance constraint from two particles.
    pub fn new(idx_a: usize, idx_b: usize, rest_length: Real, compliance: Real) -> Self {
        Self {
            idx_a,
            idx_b,
            rest_length,
            compliance,
        }
    }

    /// Project the constraint (XPBD).
    pub fn project(&self, particles: &mut [FabricParticle], dt: Real) {
        let pa = particles[self.idx_a].position;
        let pb = particles[self.idx_b].position;
        let dp = pa - pb;
        let dist = dp.norm();
        if dist < 1e-12 {
            return;
        }
        let constraint = dist - self.rest_length;
        let alpha = self.compliance / (dt * dt);
        let wa = particles[self.idx_a].inv_mass;
        let wb = particles[self.idx_b].inv_mass;
        let total_w = wa + wb + alpha;
        if total_w < 1e-12 {
            return;
        }
        let lambda = -constraint / total_w;
        let dir = dp / dist;
        if !particles[self.idx_a].is_pinned() {
            particles[self.idx_a].position += dir * (lambda * wa);
        }
        if !particles[self.idx_b].is_pinned() {
            particles[self.idx_b].position -= dir * (lambda * wb);
        }
    }
}

/// Gravity draping simulation for garments.
///
/// Simulates a garment (or fabric panel) settling under gravity over a body
/// proxy, using PBD with distance constraints, seam constraints, and sphere
/// collision for the body.
#[derive(Debug, Clone)]
pub struct GarmentDraping {
    /// Fabric collision system (particles + triangles).
    pub fabric: FabricCollision,
    /// Distance constraints (stretch and shear resistance).
    pub distance_constraints: Vec<FabricDistanceConstraint>,
    /// Seam constraints connecting panel edges.
    pub seam_constraints: Vec<SeamConstraint>,
    /// Gravity vector \[m/s²\].
    pub gravity: Vec3,
    /// Body proxy spheres (center, radius) for collision.
    pub body_proxies: Vec<(Vec3, Real)>,
    /// Number of constraint solver iterations per time step.
    pub solver_iterations: usize,
    /// Areal density of the fabric \[kg/m²\].
    pub areal_density: Real,
    /// Damping factor (velocity attenuation per step).
    pub damping: Real,
}

impl GarmentDraping {
    /// Create a new garment draping simulation.
    pub fn new(
        particles: Vec<FabricParticle>,
        triangles: Vec<[usize; 3]>,
        areal_density: Real,
        collision_thickness: Real,
        gravity: Vec3,
    ) -> Self {
        let fabric = FabricCollision::new(particles, triangles, collision_thickness, 0.3);
        Self {
            fabric,
            distance_constraints: Vec::new(),
            seam_constraints: Vec::new(),
            gravity,
            body_proxies: Vec::new(),
            solver_iterations: 10,
            areal_density,
            damping: 0.99,
        }
    }

    /// Add a distance constraint between two particles.
    pub fn add_distance_constraint(&mut self, ia: usize, ib: usize, compliance: Real) {
        let pa = self.fabric.particles[ia].position;
        let pb = self.fabric.particles[ib].position;
        let rest = (pa - pb).norm();
        self.distance_constraints
            .push(FabricDistanceConstraint::new(ia, ib, rest, compliance));
    }

    /// Add a seam constraint between two particles.
    pub fn add_seam(&mut self, ia: usize, ib: usize, stiffness: Real) {
        let pa = self.fabric.particles[ia].position;
        let pb = self.fabric.particles[ib].position;
        let rest = (pa - pb).norm();
        self.seam_constraints
            .push(SeamConstraint::new(ia, ib, rest, stiffness));
    }

    /// Add a sphere body proxy for collision.
    pub fn add_body_sphere(&mut self, center: Vec3, radius: Real) {
        self.body_proxies.push((center, radius));
    }

    /// Perform one simulation step.
    pub fn step(&mut self, dt: Real) {
        let n = self.fabric.particles.len();

        // 1. Save previous positions, apply gravity
        for p in &mut self.fabric.particles {
            p.prev_position = p.position;
            if !p.is_pinned() {
                p.velocity += self.gravity * dt;
                p.velocity *= self.damping;
                p.position += p.velocity * dt;
            }
        }

        // 2. Constraint projection
        let dc_list: Vec<FabricDistanceConstraint> = self.distance_constraints.clone();
        let sc_list: Vec<SeamConstraint> = self.seam_constraints.clone();

        for _ in 0..self.solver_iterations {
            for dc in &dc_list {
                dc.project(&mut self.fabric.particles, dt);
            }
            for sc in &sc_list {
                sc.project(&mut self.fabric.particles, dt);
            }
        }

        // 3. Collision resolution
        self.fabric.resolve_self_collisions();
        self.fabric.apply_ground_collision(0.0);

        let proxies: Vec<(Vec3, Real)> = self.body_proxies.clone();
        for (center, radius) in &proxies {
            self.fabric.apply_sphere_collision(*center, *radius);
        }

        // 4. Update velocities from positions
        let inv_dt = if dt > 1e-12 { 1.0 / dt } else { 0.0 };
        for p in &mut self.fabric.particles {
            if !p.is_pinned() {
                p.velocity = (p.position - p.prev_position) * inv_dt;
            }
        }

        let _ = n;
    }

    /// Simulate `steps` steps of `dt` seconds each.
    pub fn simulate(&mut self, dt: Real, steps: usize) {
        for _ in 0..steps {
            self.step(dt);
        }
    }

    /// Compute the average fabric height (centre of mass Y coordinate) \[m\].
    pub fn average_height(&self) -> Real {
        if self.fabric.particles.is_empty() {
            return 0.0;
        }
        let sum: Real = self.fabric.particles.iter().map(|p| p.position.y).sum();
        sum / self.fabric.particles.len() as Real
    }

    /// Compute the maximum fabric height \[m\].
    pub fn max_height(&self) -> Real {
        self.fabric
            .particles
            .iter()
            .map(|p| p.position.y)
            .fold(Real::NEG_INFINITY, Real::max)
    }

    /// Compute the minimum fabric height \[m\].
    pub fn min_height(&self) -> Real {
        self.fabric
            .particles
            .iter()
            .map(|p| p.position.y)
            .fold(Real::INFINITY, Real::min)
    }

    /// Compute the fabric drape coefficient (ratio of draped to flat area).
    ///
    /// A value of 1.0 means perfectly flat; lower values indicate more draping.
    pub fn drape_coefficient(&self) -> Real {
        if self.fabric.triangles.is_empty() {
            return 1.0;
        }
        let projected_area: Real = self
            .fabric
            .triangles
            .iter()
            .map(|tri| {
                let a = self.fabric.particles[tri[0]].position;
                let b = self.fabric.particles[tri[1]].position;
                let c = self.fabric.particles[tri[2]].position;
                // Project onto XZ plane
                let a2 = Vec3::new(a.x, 0.0, a.z);
                let b2 = Vec3::new(b.x, 0.0, b.z);
                let c2 = Vec3::new(c.x, 0.0, c.z);
                (b2 - a2).cross(&(c2 - a2)).norm() * 0.5
            })
            .sum();

        let total_area: Real = self
            .fabric
            .triangles
            .iter()
            .map(|tri| {
                let a = self.fabric.particles[tri[0]].position;
                let b = self.fabric.particles[tri[1]].position;
                let c = self.fabric.particles[tri[2]].position;
                (b - a).cross(&(c - a)).norm() * 0.5
            })
            .sum();

        if total_area < 1e-12 {
            1.0
        } else {
            projected_area / total_area
        }
    }

    /// Compute the fitting metric: fraction of fabric particles within `tolerance`
    /// of a body sphere.
    pub fn fitting_metric(&self, tolerance: Real) -> Real {
        if self.fabric.particles.is_empty() || self.body_proxies.is_empty() {
            return 0.0;
        }
        let count = self
            .fabric
            .particles
            .iter()
            .filter(|p| {
                self.body_proxies.iter().any(|(center, radius)| {
                    let dist = (p.position - center).norm();
                    (dist - radius).abs() < tolerance
                })
            })
            .count();
        count as Real / self.fabric.particles.len() as Real
    }

    /// Count the total number of seam violations exceeding `tolerance`.
    pub fn seam_violation_count(&self, tolerance: Real) -> usize {
        self.seam_constraints
            .iter()
            .filter(|sc| sc.violation(&self.fabric.particles).abs() > tolerance)
            .count()
    }
}

// ─── Utility functions ────────────────────────────────────────────────────────

/// Create a rectangular fabric mesh (grid of particles and triangles).
///
/// # Arguments
/// * `nx`, `ny` — number of particles in X and Y directions
/// * `width`, `height` — physical dimensions \[m\]
/// * `mass_per_particle` — particle mass \[kg\]
///
/// Returns `(particles, triangles)`.
pub fn create_fabric_mesh(
    nx: usize,
    ny: usize,
    width: Real,
    height: Real,
    mass_per_particle: Real,
) -> (Vec<FabricParticle>, Vec<[usize; 3]>) {
    let dx = if nx > 1 {
        width / (nx - 1) as Real
    } else {
        0.0
    };
    let dy = if ny > 1 {
        height / (ny - 1) as Real
    } else {
        0.0
    };

    let mut particles = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let pos = Vec3::new(i as Real * dx, 0.0, j as Real * dy);
            particles.push(FabricParticle::new(pos, mass_per_particle));
        }
    }

    let mut triangles = Vec::new();
    for j in 0..ny.saturating_sub(1) {
        for i in 0..nx.saturating_sub(1) {
            let idx = |ii: usize, jj: usize| jj * nx + ii;
            triangles.push([idx(i, j), idx(i + 1, j), idx(i, j + 1)]);
            triangles.push([idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]);
        }
    }
    (particles, triangles)
}

/// Compute the total surface area of a fabric mesh \[m²\].
pub fn compute_fabric_area(particles: &[FabricParticle], triangles: &[[usize; 3]]) -> Real {
    triangles
        .iter()
        .map(|tri| {
            let a = particles[tri[0]].position;
            let b = particles[tri[1]].position;
            let c = particles[tri[2]].position;
            (b - a).cross(&(c - a)).norm() * 0.5
        })
        .sum()
}

/// Compute the centre of mass of a fabric mesh.
pub fn compute_fabric_com(particles: &[FabricParticle]) -> Vec3 {
    if particles.is_empty() {
        return Vec3::zeros();
    }
    let sum: Vec3 = particles
        .iter()
        .map(|p| p.position)
        .fold(Vec3::zeros(), |a, b| a + b);
    sum / particles.len() as Real
}

/// Compute the fabric wrinkle metric: standard deviation of triangle normals from
/// the mean normal direction.
///
/// A value near 0 indicates a flat fabric; higher values indicate more wrinkling.
pub fn wrinkle_metric(particles: &[FabricParticle], triangles: &[[usize; 3]]) -> Real {
    if triangles.is_empty() {
        return 0.0;
    }
    let normals: Vec<Vec3> = triangles
        .iter()
        .map(|tri| {
            let a = particles[tri[0]].position;
            let b = particles[tri[1]].position;
            let c = particles[tri[2]].position;
            let n = (b - a).cross(&(c - a));
            let len = n.norm();
            if len > 1e-10 { n / len } else { Vec3::zeros() }
        })
        .collect();

    let n = normals.len() as Real;
    let mean: Vec3 = normals.iter().fold(Vec3::zeros(), |a, b| a + b) / n;
    let variance: Real = normals
        .iter()
        .map(|nn| (nn - mean).norm_squared())
        .sum::<Real>()
        / n;
    variance.sqrt()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::math::Vec3;

    const EPS: Real = 1e-6;

    fn make_yarn() -> YarnModel {
        YarnModel::new(100, 5e-6, 1000.0, 1e9, 4e8, 1e-4, 0.3)
    }

    fn make_woven_fabric() -> WovenFabric {
        let yarn = make_yarn();
        WovenFabric::new(yarn.clone(), yarn, WeavePattern::Plain, 500.0, 500.0)
    }

    // ── YarnModel ───────────────────────────────────────────────────────────

    #[test]
    fn test_yarn_radius_positive() {
        let yarn = make_yarn();
        assert!(yarn.yarn_radius() > 0.0, "yarn radius should be positive");
    }

    #[test]
    fn test_yarn_axial_stiffness_positive() {
        let yarn = make_yarn();
        assert!(yarn.axial_stiffness() > 0.0);
    }

    #[test]
    fn test_yarn_bending_stiffness_positive() {
        let yarn = make_yarn();
        assert!(yarn.bending_stiffness() > 0.0);
    }

    #[test]
    fn test_yarn_helix_angle_positive() {
        let yarn = make_yarn();
        let theta = yarn.helix_angle();
        assert!(
            theta > 0.0,
            "helix angle should be positive with nonzero twist"
        );
    }

    #[test]
    fn test_yarn_contact_force_zero_no_overlap() {
        let yarn = make_yarn();
        let f = yarn.contact_force(0.0, std::f64::consts::PI / 4.0);
        assert!(f.abs() < EPS, "no overlap → zero contact force");
    }

    #[test]
    fn test_yarn_contact_force_increases_with_overlap() {
        let yarn = make_yarn();
        let f1 = yarn.contact_force(1e-6, std::f64::consts::PI / 4.0);
        let f2 = yarn.contact_force(2e-6, std::f64::consts::PI / 4.0);
        assert!(f2 > f1, "contact force should increase with overlap");
    }

    #[test]
    fn test_yarn_friction_positive_for_positive_normal() {
        let yarn = make_yarn();
        let ff = yarn.friction_force(10.0);
        assert!(ff > 0.0, "friction force should be positive");
    }

    #[test]
    fn test_yarn_tensile_strength_positive() {
        let yarn = make_yarn();
        let ts = yarn.tensile_strength(1e9);
        assert!(ts > 0.0, "tensile strength should be positive");
    }

    #[test]
    fn test_yarn_strain_energy_positive_nonzero_strain() {
        let yarn = make_yarn();
        let u = yarn.strain_energy(0.01);
        assert!(u > 0.0, "strain energy should be positive");
    }

    #[test]
    fn test_yarn_more_fibers_larger_stiffness() {
        let yarn_few = YarnModel::new(10, 5e-6, 1000.0, 1e9, 4e8, 1e-4, 0.3);
        let yarn_many = YarnModel::new(1000, 5e-6, 1000.0, 1e9, 4e8, 1e-4, 0.3);
        assert!(yarn_many.axial_stiffness() > yarn_few.axial_stiffness());
    }

    // ── WeavePattern ────────────────────────────────────────────────────────

    #[test]
    fn test_plain_weave_alternates() {
        let p = WeavePattern::Plain;
        // (0,0) → warp over weft; (0,1) → weft over warp
        assert!(p.warp_over_weft(0, 0));
        assert!(!p.warp_over_weft(0, 1));
        assert!(!p.warp_over_weft(1, 0));
    }

    #[test]
    fn test_twill_repeat_size() {
        let p = WeavePattern::Twill(4);
        assert_eq!(p.repeat_size(), 4);
    }

    #[test]
    fn test_satin_repeat_size() {
        let p = WeavePattern::Satin(8);
        assert_eq!(p.repeat_size(), 8);
    }

    #[test]
    fn test_plain_weave_interlacings() {
        let p = WeavePattern::Plain;
        assert_eq!(p.interlacings_per_unit_cell(), 2);
    }

    // ── WovenFabric ─────────────────────────────────────────────────────────

    #[test]
    fn test_woven_fabric_cover_factor_bounded() {
        let fabric = make_woven_fabric();
        let cf = fabric.cover_factor();
        assert!(
            (0.0..=1.0).contains(&cf),
            "cover factor should be in [0,1]: {cf}"
        );
    }

    #[test]
    fn test_woven_fabric_porosity_bounded() {
        let fabric = make_woven_fabric();
        let por = fabric.porosity();
        assert!(
            (0.0..=1.0).contains(&por),
            "porosity should be in [0,1]: {por}"
        );
    }

    #[test]
    fn test_woven_fabric_thickness_positive() {
        let fabric = make_woven_fabric();
        assert!(fabric.thickness > 0.0);
    }

    #[test]
    fn test_woven_fabric_areal_density_positive() {
        let fabric = make_woven_fabric();
        assert!(fabric.areal_density > 0.0);
    }

    #[test]
    fn test_woven_fabric_unit_cell_size() {
        let fabric = make_woven_fabric();
        let (uc_w, uc_h) = fabric.unit_cell_size();
        assert!(uc_w > 0.0 && uc_h > 0.0);
    }

    #[test]
    fn test_woven_fabric_interlacing_matrix_shape() {
        let fabric = make_woven_fabric();
        let mat = fabric.interlacing_matrix(4);
        assert_eq!(mat.len(), 4);
        assert_eq!(mat[0].len(), 4);
    }

    #[test]
    fn test_woven_fabric_warp_crimp_positive() {
        let fabric = make_woven_fabric();
        let crimp = fabric.warp_crimp();
        assert!(crimp >= 0.0, "warp crimp should be non-negative");
    }

    #[test]
    fn test_satin_higher_tearing_strength() {
        let yarn = make_yarn();
        let plain = WovenFabric::new(
            yarn.clone(),
            yarn.clone(),
            WeavePattern::Plain,
            500.0,
            500.0,
        );
        let satin = WovenFabric::new(yarn.clone(), yarn, WeavePattern::Satin(8), 500.0, 500.0);
        assert!(
            satin.relative_tearing_strength() > plain.relative_tearing_strength(),
            "satin should have higher relative tearing strength than plain"
        );
    }

    // ── FabricMechanics ─────────────────────────────────────────────────────

    #[test]
    fn test_fabric_mechanics_tensile_stiffness_positive() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        assert!(mech.tensile_stiffness_warp > 0.0);
        assert!(mech.tensile_stiffness_weft > 0.0);
    }

    #[test]
    fn test_fabric_mechanics_shear_stiffness_positive() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        assert!(mech.shear_stiffness > 0.0);
    }

    #[test]
    fn test_fabric_mechanics_tensile_force_positive_strain() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        let (nw, nf) = mech.tensile_force(0.01, 0.005);
        assert!(
            nw > 0.0 && nf > 0.0,
            "tensile forces should be positive for positive strains"
        );
    }

    #[test]
    fn test_fabric_mechanics_bending_moment() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        let (mw, mf) = mech.bending_moment(0.1, 0.05);
        assert!(
            mw > 0.0 && mf > 0.0,
            "bending moments should be positive for positive curvatures"
        );
    }

    #[test]
    fn test_fabric_mechanics_shear_force() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        let fs = mech.shear_force(0.1);
        assert!(
            fs > 0.0,
            "shear force should be positive for positive shear angle"
        );
    }

    #[test]
    fn test_fabric_mechanics_strain_energy_positive() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        let u = mech.strain_energy_density(0.01, 0.01, 0.05);
        assert!(u > 0.0);
    }

    #[test]
    fn test_fabric_mechanics_effective_modulus_positive() {
        let fabric = make_woven_fabric();
        let mech = FabricMechanics::from_woven_fabric(&fabric);
        assert!(mech.effective_modulus_warp() > 0.0);
        assert!(mech.effective_modulus_weft() > 0.0);
        assert!(mech.effective_shear_modulus() > 0.0);
    }

    // ── KnittedFabric ───────────────────────────────────────────────────────

    #[test]
    fn test_knitted_fabric_stitch_density_positive() {
        let yarn = make_yarn();
        let kf = KnittedFabric::new(yarn, StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        assert!(kf.stitch_density() > 0.0);
    }

    #[test]
    fn test_knitted_fabric_cover_factor_bounded() {
        let yarn = make_yarn();
        let kf = KnittedFabric::new(yarn, StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        let cf = kf.cover_factor();
        assert!((0.0..=1.0).contains(&cf), "cover factor={cf}");
    }

    #[test]
    fn test_knitted_fabric_extensibility() {
        let yarn = make_yarn();
        let jersey = KnittedFabric::new(yarn.clone(), StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        let rib = KnittedFabric::new(yarn, StitchPattern::Rib, 5e-3, 200.0, 200.0);
        // Rib has higher course extensibility than jersey
        assert!(rib.max_course_extension() > jersey.max_course_extension());
    }

    #[test]
    fn test_knitted_fabric_tensile_force_positive() {
        let yarn = make_yarn();
        let kf = KnittedFabric::new(yarn, StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        let f = kf.course_tensile_force(0.1);
        assert!(
            f > 0.0,
            "tensile force should be positive for positive strain"
        );
    }

    #[test]
    fn test_knitted_fabric_loop_contact_no_compression() {
        let yarn = make_yarn();
        let kf = KnittedFabric::new(yarn, StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        let f = kf.loop_contact_force(-0.01);
        assert!(f.abs() < EPS, "no contact force for negative overlap");
    }

    // ── FabricCollision ──────────────────────────────────────────────────────

    #[test]
    fn test_fabric_mesh_creation() {
        let (particles, triangles) = create_fabric_mesh(4, 4, 1.0, 1.0, 0.01);
        assert_eq!(particles.len(), 16);
        assert_eq!(triangles.len(), 18); // (nx-1)*(ny-1)*2 = 3*3*2
    }

    #[test]
    fn test_fabric_area_positive() {
        let (particles, triangles) = create_fabric_mesh(4, 4, 1.0, 1.0, 0.01);
        let area = compute_fabric_area(&particles, &triangles);
        assert!(area > 0.0, "fabric area should be positive");
        assert!(
            (area - 1.0).abs() < 0.01,
            "area of 1×1 mesh should be ≈ 1 m²: {area}"
        );
    }

    #[test]
    fn test_fabric_com_center() {
        let (particles, _) = create_fabric_mesh(3, 3, 2.0, 2.0, 0.01);
        let com = compute_fabric_com(&particles);
        assert!((com.x - 1.0).abs() < EPS, "com.x={}", com.x);
        assert!((com.z - 1.0).abs() < EPS, "com.z={}", com.z);
    }

    #[test]
    fn test_fabric_collision_sphere() {
        // Create mesh in XZ plane (Y=0), then move particles slightly into a sphere
        let (mut particles, triangles) = create_fabric_mesh(3, 3, 0.1, 0.1, 0.001);
        // Place the sphere below the fabric so particles are guaranteed to be outside
        // the sphere after collision; sphere center is at y = -1.0 well below
        let center = Vec3::new(0.05, -1.0, 0.05);
        let sphere_radius = 0.9; // sphere extends from y=-1.9 to y=-0.1
        // Move particle 4 (center, at 0.05,0,0.05) inside the sphere: place at y=-0.5
        particles[4].position = Vec3::new(0.05, -0.5, 0.05);
        let mut fc = FabricCollision::new(particles, triangles, 0.01, 0.3);
        fc.apply_sphere_collision(center, sphere_radius);
        // Particle 4 should now be pushed outside
        let p4 = &fc.particles[4];
        let dist = (p4.position - center).norm();
        let target = sphere_radius + fc.collision_thickness;
        assert!(
            dist >= target - EPS,
            "particle should be outside sphere after collision: dist={dist}, target={target}"
        );
    }

    #[test]
    fn test_fabric_ground_collision() {
        let mut fc = FabricCollision::new(
            vec![FabricParticle::new(Vec3::new(0.0, -0.1, 0.0), 0.01)],
            vec![],
            0.001,
            0.3,
        );
        fc.apply_ground_collision(0.0);
        assert!(
            fc.particles[0].position.y >= 0.0 - EPS,
            "particle should be above ground"
        );
    }

    #[test]
    fn test_vertex_face_no_self_collision() {
        let (particles, triangles) = create_fabric_mesh(3, 3, 1.0, 1.0, 0.01);
        let fc = FabricCollision::new(particles, triangles, 0.01, 0.3);
        // Vertices of the same triangle should not self-collide
        let result = fc.vertex_face_test(0, [0, 1, 3]);
        assert!(
            result.is_none(),
            "vertex in its own triangle should not collide"
        );
    }

    #[test]
    fn test_fabric_active_particle_count() {
        let particles = vec![
            FabricParticle::new(Vec3::zeros(), 0.01),
            FabricParticle::pinned(Vec3::new(1.0, 0.0, 0.0)),
            FabricParticle::new(Vec3::new(0.5, 0.0, 0.0), 0.01),
        ];
        let _ = particles[1].is_pinned();
        let fc = FabricCollision::new(particles, vec![], 0.01, 0.3);
        assert_eq!(fc.active_particle_count(), 2, "two non-pinned particles");
    }

    // ── GarmentDraping ──────────────────────────────────────────────────────

    #[test]
    fn test_garment_draping_fabric_falls() {
        // Create a mesh elevated above the ground (y = 5.0)
        let (mut particles, triangles) = create_fabric_mesh(3, 3, 1.0, 1.0, 0.01);
        // Move all particles up to y = 5.0
        for p in &mut particles {
            p.position.y = 5.0;
            p.prev_position.y = 5.0;
        }
        let initial_avg_height = 5.0_f64;
        let mut draping =
            GarmentDraping::new(particles, triangles, 0.2, 0.001, Vec3::new(0.0, -9.81, 0.0));

        // Add distance constraints
        let nx = 3;
        let ny = 3;
        for j in 0..ny {
            for i in 0..nx - 1 {
                draping.add_distance_constraint(j * nx + i, j * nx + i + 1, 1e-6);
            }
        }
        for j in 0..ny - 1 {
            for i in 0..nx {
                draping.add_distance_constraint(j * nx + i, (j + 1) * nx + i, 1e-6);
            }
        }

        draping.simulate(1.0 / 60.0, 60);
        let final_avg_height = draping.average_height();
        assert!(
            final_avg_height < initial_avg_height,
            "fabric should fall under gravity: initial={initial_avg_height}, final={final_avg_height}"
        );
    }

    #[test]
    fn test_garment_draping_pinned_corner_stays() {
        let (mut particles, triangles) = create_fabric_mesh(2, 2, 1.0, 1.0, 0.01);
        particles[0] = FabricParticle::pinned(Vec3::new(0.0, 5.0, 0.0));
        let pinned_pos = particles[0].position;
        let mut draping = GarmentDraping::new(
            particles,
            triangles,
            0.001,
            0.001,
            Vec3::new(0.0, -9.81, 0.0),
        );
        draping.simulate(1.0 / 60.0, 10);
        let pos = draping.fabric.particles[0].position;
        assert!(
            (pos - pinned_pos).norm() < EPS,
            "pinned particle should not move"
        );
    }

    #[test]
    fn test_seam_constraint_reduces_violation() {
        let mut particles = vec![
            FabricParticle::new(Vec3::new(0.0, 0.0, 0.0), 0.01),
            FabricParticle::new(Vec3::new(1.0, 0.0, 0.0), 0.01),
        ];
        let seam = SeamConstraint::new(0, 1, 0.0, 1000.0);
        let violation_before = seam.violation(&particles);
        for _ in 0..100 {
            seam.project(&mut particles, 1.0 / 60.0);
        }
        let violation_after = seam.violation(&particles);
        assert!(
            violation_after.abs() < violation_before.abs(),
            "seam constraint should reduce violation: before={violation_before}, after={violation_after}"
        );
    }

    #[test]
    fn test_distance_constraint_reduces_error() {
        let mut particles = vec![
            FabricParticle::new(Vec3::new(0.0, 0.0, 0.0), 0.01),
            FabricParticle::new(Vec3::new(2.0, 0.0, 0.0), 0.01),
        ];
        let dc = FabricDistanceConstraint::new(0, 1, 1.0, 1e-4);
        for _ in 0..100 {
            dc.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - 1.0).abs() < 0.01,
            "distance should converge to rest length: dist={dist}"
        );
    }

    #[test]
    fn test_wrinkle_metric_flat_fabric() {
        let (particles, triangles) = create_fabric_mesh(4, 4, 1.0, 1.0, 0.01);
        let metric = wrinkle_metric(&particles, &triangles);
        assert!(
            metric < 0.01,
            "flat fabric should have near-zero wrinkle metric: {metric}"
        );
    }

    #[test]
    fn test_garment_drape_coefficient_flat() {
        let (particles, triangles) = create_fabric_mesh(4, 4, 1.0, 1.0, 0.01);
        let draping =
            GarmentDraping::new(particles, triangles, 0.2, 0.001, Vec3::new(0.0, -9.81, 0.0));
        let dc = draping.drape_coefficient();
        // Flat initial fabric (Y=0 for all particles) → projected area = total area
        assert!(
            (dc - 1.0).abs() < 0.01,
            "flat fabric drape coefficient should be 1.0: {dc}"
        );
    }

    #[test]
    fn test_fabric_particle_pinned_zero_inv_mass() {
        let p = FabricParticle::pinned(Vec3::zeros());
        assert!(p.is_pinned(), "pinned particle should have is_pinned=true");
        assert!(
            p.inv_mass < EPS,
            "pinned particle should have zero inv_mass"
        );
    }

    #[test]
    fn test_point_in_triangle() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 0.0, 1.0);
        // Centroid should be inside
        let centroid = (a + b + c) / 3.0;
        assert!(point_in_triangle_barycentric(centroid, a, b, c));
        // Point far outside should not be inside
        let outside = Vec3::new(5.0, 0.0, 5.0);
        assert!(!point_in_triangle_barycentric(outside, a, b, c));
    }

    #[test]
    fn test_stitch_pattern_extensibility_order() {
        // Rib should have higher course extensibility than jersey
        assert!(
            StitchPattern::Rib.course_extensibility()
                > StitchPattern::Jersey.course_extensibility()
        );
    }

    #[test]
    fn test_knitted_rib_vs_jersey_stiffness() {
        let yarn = make_yarn();
        let jersey = KnittedFabric::new(yarn.clone(), StitchPattern::Jersey, 5e-3, 200.0, 200.0);
        let rib = KnittedFabric::new(yarn, StitchPattern::Rib, 5e-3, 200.0, 200.0);
        // Rib has higher extensibility → check the max extension values differ
        assert!(
            (rib.max_course_extension() - jersey.max_course_extension()).abs() > EPS,
            "rib and jersey should have different extensibilities"
        );
    }

    #[test]
    fn test_garment_fitting_metric_no_proxies() {
        let (particles, triangles) = create_fabric_mesh(3, 3, 1.0, 1.0, 0.01);
        let draping =
            GarmentDraping::new(particles, triangles, 0.2, 0.001, Vec3::new(0.0, -9.81, 0.0));
        let metric = draping.fitting_metric(0.05);
        assert!(
            (metric).abs() < EPS,
            "no body proxies → fitting metric should be 0"
        );
    }

    #[test]
    fn test_garment_seam_violation_count_initially_zero() {
        let (particles, triangles) = create_fabric_mesh(2, 2, 1.0, 1.0, 0.01);
        let mut draping =
            GarmentDraping::new(particles, triangles, 0.2, 0.001, Vec3::new(0.0, -9.81, 0.0));
        draping.add_seam(0, 1, 1000.0);
        // At t=0 the seam connects adjacent particles at their rest positions
        let violations = draping.seam_violation_count(0.001);
        assert_eq!(violations, 0, "no violations at rest positions");
    }
}
