// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biological materials module: soft tissues, bone, cartilage, vascular walls,
//! cell mechanics, hydrogels, and biomechanics analysis tools.
//!
//! All quantities use SI units (Pa, m, kg, N) unless otherwise stated.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper: 3×3 matrix ops (private)
// ---------------------------------------------------------------------------

#[cfg(test)]
fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn mat3_trace(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}

fn mat3_det(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn mat3_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Compute invariants (I1, I2, J) from deformation gradient F.
fn invariants(f: &[[f64; 3]; 3]) -> (f64, f64, f64) {
    let ft = mat3_transpose(f);
    let c = mat3_mul(&ft, f);
    let i1 = mat3_trace(&c);
    let c2 = mat3_mul(&c, &c);
    let i2 = 0.5 * (i1 * i1 - mat3_trace(&c2));
    let j = mat3_det(f);
    (i1, i2, j)
}

// ---------------------------------------------------------------------------
// SoftTissue — fiber-reinforced hyperelastic (Fung-type + collagen/elastin)
// ---------------------------------------------------------------------------

/// Fiber-reinforced hyperelastic soft tissue model based on the Fung
/// exponential strain energy function with embedded collagen and elastin
/// network contributions.
///
/// The strain energy density is:
/// `W = c/2 * (exp(b*(I1-3)) - 1) + k1/(2*k2) * (exp(k2*(I4-1)^2) - 1)`
/// where I1 is the first invariant and I4 = a·Ca is the fiber pseudo-invariant.
#[derive(Debug, Clone)]
pub struct SoftTissue {
    /// Matrix stiffness parameter c (Pa).
    pub c: f64,
    /// Matrix exponential coefficient b (dimensionless).
    pub b: f64,
    /// Collagen fiber stiffness k1 (Pa).
    pub k1: f64,
    /// Collagen exponential coefficient k2 (dimensionless).
    pub k2: f64,
    /// Elastin network modulus (Pa).
    pub elastin_modulus: f64,
    /// Collagen fiber volume fraction (0–1).
    pub collagen_fraction: f64,
    /// Elastin fiber volume fraction (0–1).
    pub elastin_fraction: f64,
    /// Fiber direction unit vector \[e0, e1, e2\].
    pub fiber_direction: [f64; 3],
    /// Tissue density (kg/m³).
    pub density: f64,
}

impl SoftTissue {
    /// Default properties for human skin dermis.
    pub fn skin_dermis() -> Self {
        Self {
            c: 1200.0,
            b: 12.0,
            k1: 1500.0,
            k2: 10.0,
            elastin_modulus: 300e3,
            collagen_fraction: 0.35,
            elastin_fraction: 0.04,
            fiber_direction: [1.0, 0.0, 0.0],
            density: 1100.0,
        }
    }

    /// Default properties for human myocardium.
    pub fn myocardium() -> Self {
        Self {
            c: 500.0,
            b: 8.0,
            k1: 2000.0,
            k2: 15.0,
            elastin_modulus: 200e3,
            collagen_fraction: 0.25,
            elastin_fraction: 0.03,
            fiber_direction: [1.0, 0.0, 0.0],
            density: 1050.0,
        }
    }

    /// Evaluate the Fung-type strain energy density (J/m³).
    ///
    /// `f` is the 3×3 deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (i1, _i2, _j) = invariants(f);
        let ft = mat3_transpose(f);
        let c_mat = mat3_mul(&ft, f);
        // Fiber pseudo-invariant I4 = a · C a
        let a = &self.fiber_direction;
        let ca = [
            c_mat[0][0] * a[0] + c_mat[0][1] * a[1] + c_mat[0][2] * a[2],
            c_mat[1][0] * a[0] + c_mat[1][1] * a[1] + c_mat[1][2] * a[2],
            c_mat[2][0] * a[0] + c_mat[2][1] * a[1] + c_mat[2][2] * a[2],
        ];
        let i4 = a[0] * ca[0] + a[1] * ca[1] + a[2] * ca[2];

        // Matrix (ground substance) term — Fung exponential
        let w_matrix = self.c / 2.0 * ((self.b * (i1 - 3.0)).exp() - 1.0);

        // Fiber term (only when stretched: I4 > 1)
        let w_fiber = if i4 > 1.0 {
            let xi = i4 - 1.0;
            self.k1 / (2.0 * self.k2) * ((self.k2 * xi * xi).exp() - 1.0)
        } else {
            0.0
        };

        // Elastin isotropic neo-Hookean contribution
        let w_elastin = self.elastin_modulus / 2.0 * (i1 - 3.0);

        w_matrix + w_fiber + w_elastin
    }

    /// Estimate the tangent modulus at small strains (Pa).
    pub fn small_strain_modulus(&self) -> f64 {
        // dW/dI1 at I1 = 3 contributes 2*(c*b + elastin) to Young's modulus
        let dw_di1 = self.c * self.b + self.elastin_modulus;
        4.0 * dw_di1
    }

    /// Collagen network toe-region end strain (dimensionless).
    ///
    /// Below this stretch ratio the collagen fibers are not yet load-bearing.
    pub fn toe_region_stretch(&self) -> f64 {
        1.0 + (1.0 / (2.0 * self.k2)).sqrt()
    }

    /// Compute the Cauchy stress in the fiber direction at a given uniaxial
    /// stretch `lambda` along the fiber axis.
    pub fn uniaxial_cauchy_stress(&self, lambda: f64) -> f64 {
        if lambda <= 0.0 {
            return 0.0;
        }
        // For incompressible uniaxial: F = diag(lambda, 1/sqrt(lambda), 1/sqrt(lambda))
        let lambda_perp = 1.0 / lambda.sqrt();
        let i1 = lambda * lambda + 2.0 * lambda_perp * lambda_perp;
        let i4 = lambda * lambda; // fiber aligned with stretch

        let dw_di1 = self.c * self.b * (self.b * (i1 - 3.0)).exp() + self.elastin_modulus;

        let dw_di4 = if i4 > 1.0 {
            let xi = i4 - 1.0;
            self.k1 * xi * (self.k2 * xi * xi).exp()
        } else {
            0.0
        };

        // Cauchy stress = 2*(lambda^2 dW/dI1 - lambda_perp^2 dW/dI1) + 2*lambda^2*dW/dI4

        2.0 * dw_di1 * (lambda * lambda - lambda_perp * lambda_perp)
            + 2.0 * lambda * lambda * dw_di4
    }
}

// ---------------------------------------------------------------------------
// BoneMaterial — cortical/cancellous, anisotropic elastic, Reuss–Voigt bounds
// ---------------------------------------------------------------------------

/// Anisotropic elastic model for bone (cortical and cancellous).
///
/// Cortical bone is treated as a transversely isotropic material with
/// independent axial (longitudinal) and transverse stiffnesses.
#[derive(Debug, Clone)]
pub struct BoneMaterial {
    /// Longitudinal (axial) Young's modulus E_L (Pa).  Cortical: ~17–25 GPa.
    pub e_longitudinal: f64,
    /// Transverse Young's modulus E_T (Pa).  Cortical: ~10–15 GPa.
    pub e_transverse: f64,
    /// Shear modulus G_LT (Pa).
    pub g_lt: f64,
    /// Poisson's ratio ν_LT (longitudinal–transverse).
    pub nu_lt: f64,
    /// Poisson's ratio ν_TT (transverse–transverse).
    pub nu_tt: f64,
    /// Bulk density (kg/m³).
    pub density: f64,
    /// Ultimate tensile strength (Pa).
    pub uts: f64,
    /// Ultimate compressive strength (Pa).
    pub ucs: f64,
    /// Volume fraction of mineralised matrix (0–1) for cancellous bone.
    pub volume_fraction: f64,
    /// Whether this represents cortical (`true`) or cancellous (`false`) bone.
    pub is_cortical: bool,
}

impl BoneMaterial {
    /// Human femoral cortical bone (mid-diaphysis).
    pub fn cortical_femur() -> Self {
        Self {
            e_longitudinal: 20.0e9,
            e_transverse: 12.0e9,
            g_lt: 4.5e9,
            nu_lt: 0.29,
            nu_tt: 0.63,
            density: 1900.0,
            uts: 120.0e6,
            ucs: 180.0e6,
            volume_fraction: 1.0,
            is_cortical: true,
        }
    }

    /// Vertebral cancellous bone (trabecular, ~15 % volume fraction).
    pub fn cancellous_vertebra() -> Self {
        Self {
            e_longitudinal: 0.3e9,
            e_transverse: 0.15e9,
            g_lt: 0.05e9,
            nu_lt: 0.25,
            nu_tt: 0.40,
            density: 300.0,
            uts: 5.0e6,
            ucs: 8.0e6,
            volume_fraction: 0.15,
            is_cortical: false,
        }
    }

    /// Voigt upper bound on effective longitudinal modulus for a two-phase
    /// composite (bone matrix + marrow/pore space).
    ///
    /// `e_phase2` is the modulus of the second phase (Pa).
    pub fn voigt_upper_bound(&self, e_phase2: f64) -> f64 {
        self.volume_fraction * self.e_longitudinal + (1.0 - self.volume_fraction) * e_phase2
    }

    /// Reuss lower bound on effective longitudinal modulus.
    ///
    /// `e_phase2` is the modulus of the second phase (Pa).
    pub fn reuss_lower_bound(&self, e_phase2: f64) -> f64 {
        let vf = self.volume_fraction;
        1.0 / (vf / self.e_longitudinal + (1.0 - vf) / e_phase2)
    }

    /// Hashin–Shtrikman arithmetic average of the two bounds.
    pub fn hashin_shtrikman_estimate(&self, e_phase2: f64) -> f64 {
        0.5 * (self.voigt_upper_bound(e_phase2) + self.reuss_lower_bound(e_phase2))
    }

    /// Morgan–Keaveny power-law density–modulus relationship.
    ///
    /// `ash_density_g_cm3` is the ash density in g/cm³.
    pub fn modulus_from_ash_density(ash_density_g_cm3: f64) -> f64 {
        // Keaveny & Hayes (1993) empirical fit for trabecular bone
        6.95e9 * ash_density_g_cm3.powf(1.49)
    }

    /// Anisotropy ratio (E_L / E_T).
    pub fn anisotropy_ratio(&self) -> f64 {
        self.e_longitudinal / self.e_transverse
    }

    /// Longitudinal wave speed (m/s).
    pub fn longitudinal_wave_speed(&self) -> f64 {
        (self.e_longitudinal / self.density).sqrt()
    }

    /// Transverse (shear) wave speed (m/s).
    pub fn shear_wave_speed(&self) -> f64 {
        (self.g_lt / self.density).sqrt()
    }

    /// Yield strain in compression (dimensionless, ~0.7–1 % for cortical).
    pub fn compressive_yield_strain(&self) -> f64 {
        self.ucs / self.e_longitudinal
    }
}

// ---------------------------------------------------------------------------
// CartilageModel — biphasic theory, creep, permeability
// ---------------------------------------------------------------------------

/// Biphasic (solid + fluid) model for articular cartilage.
///
/// Based on Mow et al. (1980) biphasic theory:
/// - Solid phase: linear elastic with Young's modulus `E_s` and Poisson's ratio `nu_s`.
/// - Fluid phase: interstitial water with hydraulic permeability `k`.
/// - Creep and stress relaxation arise from fluid exudation.
#[derive(Debug, Clone)]
pub struct CartilageModel {
    /// Solid-phase aggregate modulus H_A = E_s*(1-nu_s)/((1+nu_s)*(1-2*nu_s)) (Pa).
    pub aggregate_modulus: f64,
    /// Solid-phase Young's modulus E_s (Pa).  Typical: 0.5–2 MPa.
    pub e_solid: f64,
    /// Solid-phase Poisson's ratio ν_s.
    pub nu_solid: f64,
    /// Fluid-phase hydraulic permeability k (m⁴/N·s).  Typical: ~1e-15 m⁴/N·s.
    pub permeability: f64,
    /// Cartilage thickness (m).
    pub thickness: f64,
    /// Fluid volume fraction (porosity, 0–1).  Typical: 0.65–0.80.
    pub porosity: f64,
    /// Fixed charge density (mEq/mL) — for swelling pressure computations.
    pub fixed_charge_density: f64,
    /// Tissue density (kg/m³).
    pub density: f64,
}

impl CartilageModel {
    /// Typical human articular cartilage (femoral condyle).
    ///
    /// Alias for [`Self::human_knee_cartilage`].
    pub fn articular_cartilage() -> Self {
        Self::human_knee_cartilage()
    }

    /// Typical human articular cartilage (femoral condyle).
    pub fn human_knee_cartilage() -> Self {
        let e_solid = 0.8e6;
        let nu = 0.15;
        let ha = e_solid * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        Self {
            aggregate_modulus: ha,
            e_solid,
            nu_solid: nu,
            permeability: 1.0e-15,
            thickness: 3.0e-3,
            porosity: 0.75,
            fixed_charge_density: 0.1,
            density: 1100.0,
        }
    }

    /// Compute the aggregate modulus from E and ν.
    pub fn compute_aggregate_modulus(e: f64, nu: f64) -> f64 {
        e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Characteristic time constant for creep τ = h²/(H_A * k) (s).
    ///
    /// `h` is the cartilage thickness (m).
    pub fn creep_time_constant(&self) -> f64 {
        let h = self.thickness;
        h * h / (self.aggregate_modulus * self.permeability)
    }

    /// Biphasic creep compliance at time `t` (Pa⁻¹).
    ///
    /// Approximation using the first-term series solution from Mow et al. (1980).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let tau = self.creep_time_constant();
        let c_inf = 1.0 / self.aggregate_modulus;
        // First-term exponential approximation
        let c0 = 1.0 / (self.aggregate_modulus + 4.0 / 3.0 * self.e_solid);
        c_inf - (c_inf - c0) * (-t / tau).exp()
    }

    /// Donnan osmotic swelling pressure (Pa) at a given NaCl concentration.
    ///
    /// `c_ext` is the external NaCl concentration (mEq/mL).
    pub fn donnan_swelling_pressure(&self, c_ext: f64) -> f64 {
        // Simplified Donnan equilibrium: π = RT*(sqrt(cF^2/4 + c_ext^2) - c_ext)
        // using RT ≈ 2436 Pa·mL/mEq at 37°C
        let rt = 2436.0;
        let cf = self.fixed_charge_density;
        rt * ((cf * cf / 4.0 + c_ext * c_ext).sqrt() - c_ext)
    }

    /// Darcy velocity (fluid flux, m/s) under an applied pressure gradient `dp_dz` (Pa/m).
    pub fn darcy_fluid_flux(&self, dp_dz: f64) -> f64 {
        self.permeability * dp_dz
    }

    /// Strain-dependent permeability (Holmes–Mow model).
    ///
    /// `e` is the dilatational strain (volumetric, negative in compression).
    pub fn strain_dependent_permeability(&self, e: f64) -> f64 {
        // Holmes & Mow (1990): k(e) = k0 * exp(M * e)
        let m = 4.0; // empirical parameter
        self.permeability * (m * e).exp()
    }
}

// ---------------------------------------------------------------------------
// VascularWall — residual stress, opening angle, tri-layered structure
// ---------------------------------------------------------------------------

/// Model for arterial wall mechanics including residual stress, opening angle,
/// and three-layer (intima, media, adventitia) structure.
///
/// Based on the Holzapfel–Gasser–Ogden (HGO) model for arterial mechanics.
#[derive(Debug, Clone)]
pub struct VascularWall {
    /// Inner radius (unloaded, m).
    pub inner_radius: f64,
    /// Outer radius (unloaded, m).
    pub outer_radius: f64,
    /// Opening angle α (rad) — measures residual stress state.
    pub opening_angle: f64,
    /// Intima layer thickness (m).
    pub intima_thickness: f64,
    /// Media layer thickness (m).
    pub media_thickness: f64,
    /// Adventitia layer thickness (m).
    pub adventitia_thickness: f64,
    /// Isotropic neo-Hookean parameter c₁₀ (Pa) for the ground matrix.
    pub c10: f64,
    /// Fiber stiffness k1 (Pa) — HGO model.
    pub k1: f64,
    /// Fiber nonlinearity k2 (dimensionless) — HGO model.
    pub k2: f64,
    /// Mean fiber angle (rad) from the circumferential direction.
    pub fiber_angle: f64,
    /// Collagen fraction of wall composition (0–1).
    pub collagen_fraction: f64,
    /// Smooth muscle cell fraction (0–1).
    pub smooth_muscle_fraction: f64,
    /// Blood vessel density (kg/m³).
    pub density: f64,
}

impl VascularWall {
    /// Typical human aorta (thoracic).
    pub fn human_aorta() -> Self {
        Self {
            inner_radius: 12.5e-3,
            outer_radius: 15.0e-3,
            opening_angle: 100.0_f64.to_radians(),
            intima_thickness: 0.3e-3,
            media_thickness: 1.5e-3,
            adventitia_thickness: 0.7e-3,
            c10: 3000.0,
            k1: 2000.0,
            k2: 11.0,
            fiber_angle: 49.0_f64.to_radians(),
            collagen_fraction: 0.22,
            smooth_muscle_fraction: 0.15,
            density: 1050.0,
        }
    }

    /// Typical human coronary artery.
    pub fn coronary_artery() -> Self {
        Self {
            inner_radius: 1.5e-3,
            outer_radius: 2.2e-3,
            opening_angle: 115.0_f64.to_radians(),
            intima_thickness: 0.1e-3,
            media_thickness: 0.5e-3,
            adventitia_thickness: 0.1e-3,
            c10: 5000.0,
            k1: 3500.0,
            k2: 14.0,
            fiber_angle: 44.0_f64.to_radians(),
            collagen_fraction: 0.28,
            smooth_muscle_fraction: 0.20,
            density: 1050.0,
        }
    }

    /// Wall thickness (m).
    pub fn wall_thickness(&self) -> f64 {
        self.outer_radius - self.inner_radius
    }

    /// Circumferential residual stretch due to the opening angle.
    ///
    /// Returns the inner/outer residual stretch ratios at the given
    /// radial position `r` (m).
    pub fn residual_stretch_circumferential(&self, r: f64) -> f64 {
        // Fung & Liu (1989) formula for opening angle model
        let alpha = self.opening_angle;
        let ri = self.inner_radius;
        let ro = self.outer_radius;
        let r_mid = 0.5 * (ri + ro);
        // Approximate: lambda_theta = r / r_mid * (PI / (PI - alpha))
        (r / r_mid) * (PI / (PI - alpha))
    }

    /// HGO strain energy density (J/m³) at deformation gradient F.
    pub fn hgo_strain_energy(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (i1, _i2, j) = invariants(f);
        let ft = mat3_transpose(f);
        let c_mat = mat3_mul(&ft, f);

        // Fiber pseudo-invariant I4 (fiber aligned at ±fiber_angle from circumferential)
        let cos_a = self.fiber_angle.cos();
        let sin_a = self.fiber_angle.sin();
        let a1 = [cos_a, sin_a, 0.0];
        let i4 = a1[0] * (c_mat[0][0] * a1[0] + c_mat[0][1] * a1[1])
            + a1[1] * (c_mat[1][0] * a1[0] + c_mat[1][1] * a1[1]);

        // Penalty term for incompressibility
        let d = 1.0 / (2.0 * self.c10); // bulk modulus-like
        let w_vol = (j - 1.0).powi(2) / (2.0 * d);

        // Isochoric matrix
        let j_23 = j.powf(-2.0 / 3.0);
        let i1_bar = j_23 * i1;
        let w_iso = self.c10 * (i1_bar - 3.0);

        // Fiber (only when I4 > 1)
        let w_fiber = if i4 > 1.0 {
            let e = i4 - 1.0;
            self.k1 / (2.0 * self.k2) * ((self.k2 * e * e).exp() - 1.0)
        } else {
            0.0
        };

        w_vol + w_iso + w_fiber
    }

    /// Laplace tension (N/m) for a thin-walled cylinder at internal pressure p (Pa).
    pub fn laplace_tension(&self, internal_pressure: f64) -> f64 {
        internal_pressure * self.inner_radius
    }

    /// Circumferential wall stress (Pa) at mean wall radius under Laplace.
    pub fn circumferential_stress(&self, internal_pressure: f64) -> f64 {
        let h = self.wall_thickness();
        self.laplace_tension(internal_pressure) / h
    }

    /// Compliance (m/Pa) per unit length of vessel.
    ///
    /// Estimated from linearised pressure–radius relationship.
    pub fn vascular_compliance(&self, e_wall: f64) -> f64 {
        let r = 0.5 * (self.inner_radius + self.outer_radius);
        let h = self.wall_thickness();
        2.0 * PI * r * r * r / (e_wall * h)
    }
}

// ---------------------------------------------------------------------------
// CellMechanics — cortical tension, cytoskeletal stiffness, contractility
// ---------------------------------------------------------------------------

/// Mechanical model for a single eukaryotic cell including cortical tension,
/// cytoskeletal prestress, and active actomyosin contractility.
///
/// Based on the tensegrity model (Ingber, 1997) and Zaman et al. (2007).
#[derive(Debug, Clone)]
pub struct CellMechanics {
    /// Cortical tension (N/m) — surface tension of the cell cortex.
    pub cortical_tension: f64,
    /// Cell radius (m).
    pub radius: f64,
    /// Cytoskeletal Young's modulus (Pa).  Typical: 100–1000 Pa.
    pub cytoskeletal_modulus: f64,
    /// Active contractility stress (Pa) generated by actomyosin.
    pub active_stress: f64,
    /// Adhesion energy density (J/m²) — cell–substrate adhesion.
    pub adhesion_energy: f64,
    /// Cell density (kg/m³).
    pub density: f64,
    /// Number of actin stress fibers.
    pub num_stress_fibers: usize,
    /// Stress fiber prestress (Pa).
    pub prestress: f64,
}

impl CellMechanics {
    /// Typical adherent fibroblast.
    pub fn fibroblast() -> Self {
        Self {
            cortical_tension: 0.05e-3,
            radius: 10.0e-6,
            cytoskeletal_modulus: 400.0,
            active_stress: 1000.0,
            adhesion_energy: 1.0e-4,
            density: 1060.0,
            num_stress_fibers: 20,
            prestress: 500.0,
        }
    }

    /// Typical red blood cell.
    pub fn red_blood_cell() -> Self {
        Self {
            cortical_tension: 0.002e-3,
            radius: 4.0e-6,
            cytoskeletal_modulus: 6.0,
            active_stress: 0.0,
            adhesion_energy: 0.0,
            density: 1080.0,
            num_stress_fibers: 0,
            prestress: 0.0,
        }
    }

    /// Laplace pressure from cortical tension (Pa).
    pub fn laplace_pressure(&self) -> f64 {
        2.0 * self.cortical_tension / self.radius
    }

    /// Effective Young's modulus including active contractility (Pa).
    pub fn effective_modulus(&self) -> f64 {
        self.cytoskeletal_modulus + self.active_stress
    }

    /// Spreading area (m²) predicted from adhesion energy balance.
    ///
    /// Equates adhesion energy to cortical tension energy: A = π r² W / T_c.
    pub fn spreading_area(&self) -> f64 {
        let w = self.adhesion_energy;
        let tc = self.cortical_tension;
        if tc < 1e-20 {
            return 0.0;
        }
        PI * self.radius * self.radius * w / tc
    }

    /// Bending stiffness of the cell membrane (J), approximate.
    ///
    /// Uses the result from the Helfrich model: κ = 25 k_B T.
    pub fn membrane_bending_stiffness(&self) -> f64 {
        // k_B T at 37°C ≈ 4.28e-21 J
        let kbt = 4.28e-21_f64;
        25.0 * kbt
    }

    /// Cell volume (m³) assuming spherical shape.
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }

    /// Osmotic pressure (Pa) from a small volume change `dv/v`.
    ///
    /// Uses a linear osmotic modulus approximation K_osm ≈ E_cyto.
    pub fn osmotic_pressure(&self, volumetric_strain: f64) -> f64 {
        -self.cytoskeletal_modulus * volumetric_strain
    }
}

// ---------------------------------------------------------------------------
// HydrogelModel — swelling, crosslink density, Flory–Rehner theory
// ---------------------------------------------------------------------------

/// Hydrogel mechanical and swelling model based on Flory–Rehner theory.
///
/// Combines elastic (rubber-like) network energy with Flory–Huggins
/// polymer–solvent mixing free energy.
#[derive(Debug, Clone)]
pub struct HydrogelModel {
    /// Polymer volume fraction at synthesis (reference state, 0–1).
    pub phi0: f64,
    /// Equilibrium polymer volume fraction in swollen state (0–1).
    pub phi_eq: f64,
    /// Flory–Huggins interaction parameter χ (dimensionless).  Typical: 0.4–0.6.
    pub chi: f64,
    /// Crosslink density (mol/m³).
    pub crosslink_density: f64,
    /// Dry gel modulus (Pa).
    pub dry_modulus: f64,
    /// Solvent molar volume v_s (m³/mol).  Water: ~1.8e-5 m³/mol.
    pub solvent_molar_volume: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl HydrogelModel {
    /// PEG hydrogel typical parameters (10 % w/w, moderate crosslinking).
    pub fn peg_hydrogel() -> Self {
        Self {
            phi0: 0.10,
            phi_eq: 0.08,
            chi: 0.45,
            crosslink_density: 500.0,
            dry_modulus: 50.0e3,
            solvent_molar_volume: 1.8e-5,
            temperature: 310.0,
        }
    }

    /// Hyaluronic acid hydrogel.
    pub fn hyaluronic_acid() -> Self {
        Self {
            phi0: 0.02,
            phi_eq: 0.015,
            chi: 0.50,
            crosslink_density: 100.0,
            dry_modulus: 1.0e3,
            solvent_molar_volume: 1.8e-5,
            temperature: 310.0,
        }
    }

    /// Gas constant R (J/mol·K).
    const R: f64 = 8.314;

    /// Rubber-elastic network contribution to the swelling stress (Pa).
    ///
    /// Based on affine network theory: G_el = n_c * R * T * phi_eq^(1/3).
    pub fn network_modulus(&self) -> f64 {
        let n_c = self.crosslink_density;
        let phi = self.phi_eq;
        n_c * Self::R * self.temperature * phi.powf(1.0 / 3.0)
    }

    /// Flory–Huggins mixing chemical potential of the solvent (J/mol).
    ///
    /// Δμ_s = RT * \[ln(1-φ) + φ + χ φ²\]
    pub fn mixing_chemical_potential(&self) -> f64 {
        let phi = self.phi_eq;
        Self::R * self.temperature * ((1.0 - phi).ln() + phi + self.chi * phi * phi)
    }

    /// Osmotic (swelling) pressure Π (Pa) at equilibrium.
    ///
    /// Π = -R T / v_s * \[ln(1-φ) + φ + χ φ²\] + G_el*(φ^(1/3) - φ/2)
    pub fn osmotic_pressure(&self) -> f64 {
        let phi = self.phi_eq;
        let rtvs = Self::R * self.temperature / self.solvent_molar_volume;
        let pi_mix = -rtvs * ((1.0 - phi).ln() + phi + self.chi * phi * phi);
        let pi_el = self.network_modulus() * (phi.powf(1.0 / 3.0) - phi / 2.0);
        pi_mix + pi_el
    }

    /// Swelling ratio Q = V_swollen / V_dry = 1/φ_eq.
    pub fn swelling_ratio(&self) -> f64 {
        1.0 / self.phi_eq
    }

    /// Linear swelling strain ε = Q^(1/3) - 1.
    pub fn linear_swelling_strain(&self) -> f64 {
        self.swelling_ratio().powf(1.0 / 3.0) - 1.0
    }

    /// Gel modulus in the swollen state (Pa) via affine network theory.
    pub fn swollen_modulus(&self) -> f64 {
        self.network_modulus()
    }

    /// Diffusivity of solvent through the gel (m²/s) using free-volume theory.
    ///
    /// `d0` is the free-solution diffusivity (m²/s).
    pub fn effective_diffusivity(&self, d0: f64) -> f64 {
        // Obstruction model: D_eff = D0 * exp(-π r^2 L_c c / 4) simplified
        let phi = self.phi_eq;
        d0 * (-0.84 * phi / (1.0 - phi)).exp()
    }
}

// ---------------------------------------------------------------------------
// BiomechanicsAnalysis — stress shielding, failure envelope, fatigue life
// ---------------------------------------------------------------------------

/// Analysis tools for biomechanical systems: stress shielding in orthopaedic
/// implants, failure envelopes for bone/tissue, and fatigue life estimation.
#[derive(Debug, Clone)]
pub struct BiomechanicsAnalysis {
    /// Implant (or prosthesis) Young's modulus (Pa).
    pub implant_modulus: f64,
    /// Host bone Young's modulus (Pa).
    pub bone_modulus: f64,
    /// Cross-sectional area of implant (m²).
    pub implant_area: f64,
    /// Cross-sectional area of bone (m²).
    pub bone_area: f64,
    /// Applied axial load (N).
    pub applied_load: f64,
    /// Fatigue exponent m (material constant, typically 2–6 for bone).
    pub fatigue_exponent: f64,
    /// Reference fatigue stress σ₀ (Pa) at N₀ cycles.
    pub reference_fatigue_stress: f64,
    /// Reference cycle count N₀.
    pub reference_cycles: f64,
}

impl BiomechanicsAnalysis {
    /// Typical hip replacement scenario (titanium stem in femoral cortex).
    pub fn hip_replacement() -> Self {
        Self {
            implant_modulus: 110.0e9, // Ti-6Al-4V
            bone_modulus: 18.0e9,
            implant_area: 200.0e-6,
            bone_area: 400.0e-6,
            applied_load: 2500.0,
            fatigue_exponent: 6.0,
            reference_fatigue_stress: 60.0e6,
            reference_cycles: 1.0e6,
        }
    }

    /// Stress shielding factor (0–1): fraction of load carried by the implant.
    ///
    /// A value close to 1 means the implant carries almost all the load (severe
    /// stress shielding); 0 means no shielding.
    pub fn stress_shielding_factor(&self) -> f64 {
        let ea_implant = self.implant_modulus * self.implant_area;
        let ea_bone = self.bone_modulus * self.bone_area;
        ea_implant / (ea_implant + ea_bone)
    }

    /// Stress in the bone with the implant present (Pa).
    pub fn bone_stress_with_implant(&self) -> f64 {
        let f = self.stress_shielding_factor();
        (1.0 - f) * self.applied_load / self.bone_area
    }

    /// Stress in the bone without the implant (Pa).
    pub fn bone_stress_without_implant(&self) -> f64 {
        self.applied_load / self.bone_area
    }

    /// Tsai–Wu failure criterion for an orthotropic material.
    ///
    /// Returns the failure index (≥1 means failure).
    /// - `sigma11` axial stress (Pa)
    /// - `sigma22` transverse stress (Pa)
    /// - `tau12` shear stress (Pa)
    /// - `f1t` tensile strength in 1-direction (Pa)
    /// - `f1c` compressive strength in 1-direction (Pa)
    /// - `f2t` tensile strength in 2-direction (Pa)
    /// - `f2c` compressive strength in 2-direction (Pa)
    /// - `f12` shear strength (Pa)
    pub fn tsai_wu_failure_index(
        sigma11: f64,
        sigma22: f64,
        tau12: f64,
        f1t: f64,
        f1c: f64,
        f2t: f64,
        f2c: f64,
        f12: f64,
    ) -> f64 {
        let h1 = 1.0 / f1t - 1.0 / f1c;
        let h2 = 1.0 / f2t - 1.0 / f2c;
        let h11 = 1.0 / (f1t * f1c);
        let h22 = 1.0 / (f2t * f2c);
        let h66 = 1.0 / (f12 * f12);
        let h12 = -0.5 * (h11 * h22).sqrt(); // interactive term
        h1 * sigma11
            + h2 * sigma22
            + h11 * sigma11 * sigma11
            + h22 * sigma22 * sigma22
            + h66 * tau12 * tau12
            + 2.0 * h12 * sigma11 * sigma22
    }

    /// Von Mises equivalent stress (Pa).
    pub fn von_mises_stress(s11: f64, s22: f64, s33: f64, s12: f64, s23: f64, s13: f64) -> f64 {
        (0.5 * ((s11 - s22).powi(2)
            + (s22 - s33).powi(2)
            + (s33 - s11).powi(2)
            + 6.0 * (s12 * s12 + s23 * s23 + s13 * s13)))
            .sqrt()
    }

    /// Bone fatigue life (cycles) under cyclic stress amplitude `sigma_a` (Pa).
    ///
    /// Uses the power-law S-N relationship: N = N₀ * (σ₀/σ_a)^m.
    pub fn fatigue_life(&self, sigma_a: f64) -> f64 {
        if sigma_a <= 0.0 {
            return f64::INFINITY;
        }
        self.reference_cycles
            * (self.reference_fatigue_stress / sigma_a).powf(self.fatigue_exponent)
    }

    /// Miner's cumulative damage rule for variable-amplitude loading.
    ///
    /// Returns the total damage D (failure when D ≥ 1).
    /// `cycles_at_stress` is a slice of (n_applied, sigma_a) pairs.
    pub fn miners_damage(&self, cycles_at_stress: &[(f64, f64)]) -> f64 {
        cycles_at_stress
            .iter()
            .map(|(n, sigma)| n / self.fatigue_life(*sigma))
            .sum()
    }

    /// Maximum principal stress failure criterion.
    ///
    /// Returns `true` if any principal stress exceeds the material strength.
    pub fn max_principal_stress_failure(
        &self,
        sigma1: f64,
        sigma2: f64,
        sigma3: f64,
        tensile_strength: f64,
        compressive_strength: f64,
    ) -> bool {
        sigma1 > tensile_strength
            || sigma2 > tensile_strength
            || sigma3 > tensile_strength
            || sigma1 < -compressive_strength
            || sigma2 < -compressive_strength
            || sigma3 < -compressive_strength
    }

    /// Remodelling stimulus (dimensionless) — Huiskes et al. strain energy model.
    ///
    /// `sed` is the strain energy density (J/m³), `sed_ref` the reference SED.
    /// Returns positive for bone apposition, negative for resorption.
    pub fn remodelling_stimulus(sed: f64, sed_ref: f64, dead_zone: f64) -> f64 {
        let ratio = (sed - sed_ref) / sed_ref;
        if ratio.abs() < dead_zone {
            0.0
        } else {
            ratio - dead_zone.copysign(ratio)
        }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // --- SoftTissue ---

    #[test]
    fn test_soft_tissue_strain_energy_at_rest() {
        let st = SoftTissue::skin_dermis();
        let f_id = mat3_identity();
        let w = st.strain_energy_density(&f_id);
        // At identity, I1=3 so W_matrix ≈ 0 and W_fiber ≈ 0 (I4=1 → no fiber term)
        assert!(
            w.abs() < 1e-3,
            "strain energy at rest should be ~0, got {w}"
        );
    }

    #[test]
    fn test_soft_tissue_strain_energy_positive_under_stretch() {
        let st = SoftTissue::skin_dermis();
        let lambda = 1.2_f64;
        let lp = 1.0 / lambda.sqrt();
        let f = [[lambda, 0.0, 0.0], [0.0, lp, 0.0], [0.0, 0.0, lp]];
        let w = st.strain_energy_density(&f);
        assert!(
            w > 0.0,
            "strain energy should be positive under stretch, got {w}"
        );
    }

    #[test]
    fn test_soft_tissue_cauchy_stress_positive_stretch() {
        let st = SoftTissue::skin_dermis();
        let sigma = st.uniaxial_cauchy_stress(1.3);
        assert!(sigma > 0.0, "Cauchy stress should be positive, got {sigma}");
    }

    #[test]
    fn test_soft_tissue_cauchy_stress_at_one() {
        let st = SoftTissue::skin_dermis();
        let sigma = st.uniaxial_cauchy_stress(1.0);
        // At λ=1, (λ² - λ_perp²) = 0 so stress = 0
        assert!(
            sigma.abs() < 1e-3,
            "Cauchy stress at λ=1 should be ~0, got {sigma}"
        );
    }

    #[test]
    fn test_soft_tissue_small_strain_modulus_positive() {
        let st = SoftTissue::myocardium();
        let e = st.small_strain_modulus();
        assert!(e > 0.0, "small-strain modulus should be positive, got {e}");
    }

    #[test]
    fn test_soft_tissue_toe_region_stretch_greater_than_one() {
        let st = SoftTissue::skin_dermis();
        let toe = st.toe_region_stretch();
        assert!(toe > 1.0, "toe-region stretch should be >1, got {toe}");
    }

    #[test]
    fn test_soft_tissue_stress_increases_with_stretch() {
        let st = SoftTissue::skin_dermis();
        let s1 = st.uniaxial_cauchy_stress(1.1);
        let s2 = st.uniaxial_cauchy_stress(1.5);
        assert!(
            s2 > s1,
            "stress should increase with stretch: s1={s1}, s2={s2}"
        );
    }

    // --- BoneMaterial ---

    #[test]
    fn test_bone_anisotropy_ratio() {
        let bone = BoneMaterial::cortical_femur();
        let r = bone.anisotropy_ratio();
        assert!(
            (r - 20.0e9 / 12.0e9).abs() < TOL,
            "anisotropy ratio should be E_L/E_T"
        );
    }

    #[test]
    fn test_bone_voigt_reuss_bounds() {
        let bone = BoneMaterial::cancellous_vertebra();
        let e2 = 1.0e6; // marrow modulus ~1 MPa
        let voigt = bone.voigt_upper_bound(e2);
        let reuss = bone.reuss_lower_bound(e2);
        assert!(
            voigt >= reuss,
            "Voigt bound should be ≥ Reuss bound: {voigt} vs {reuss}"
        );
    }

    #[test]
    fn test_bone_hs_estimate_between_bounds() {
        let bone = BoneMaterial::cortical_femur();
        let e2 = 0.1e9;
        let v = bone.voigt_upper_bound(e2);
        let r = bone.reuss_lower_bound(e2);
        let hs = bone.hashin_shtrikman_estimate(e2);
        assert!(
            hs >= r && hs <= v,
            "HS estimate {hs} should be between Reuss {r} and Voigt {v}"
        );
    }

    #[test]
    fn test_bone_wave_speeds_positive() {
        let bone = BoneMaterial::cortical_femur();
        assert!(bone.longitudinal_wave_speed() > 0.0);
        assert!(bone.shear_wave_speed() > 0.0);
    }

    #[test]
    fn test_bone_compressive_yield_strain() {
        let bone = BoneMaterial::cortical_femur();
        let eps = bone.compressive_yield_strain();
        assert!(
            eps > 0.0 && eps < 0.02,
            "yield strain should be ~0.9 %, got {eps}"
        );
    }

    #[test]
    fn test_bone_ash_density_modulus_power_law() {
        let e = BoneMaterial::modulus_from_ash_density(0.8);
        assert!(e > 0.0, "modulus should be positive, got {e}");
    }

    // --- CartilageModel ---

    #[test]
    fn test_cartilage_aggregate_modulus_formula() {
        let e = 0.8e6;
        let nu = 0.15;
        let ha = CartilageModel::compute_aggregate_modulus(e, nu);
        let expected = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        assert!((ha - expected).abs() < 1.0, "aggregate modulus mismatch");
    }

    #[test]
    fn test_cartilage_creep_time_constant_positive() {
        let cart = CartilageModel::human_knee_cartilage();
        let tau = cart.creep_time_constant();
        assert!(
            tau > 0.0,
            "creep time constant should be positive, got {tau}"
        );
    }

    #[test]
    fn test_cartilage_creep_compliance_monotone() {
        let cart = CartilageModel::human_knee_cartilage();
        let c0 = cart.creep_compliance(0.0);
        let c100 = cart.creep_compliance(100.0);
        let cinf = cart.creep_compliance(1.0e12);
        assert!(c0 <= c100, "compliance should increase over time");
        assert!(c100 <= cinf, "compliance should approach equilibrium");
    }

    #[test]
    fn test_cartilage_darcy_flux() {
        let cart = CartilageModel::human_knee_cartilage();
        let flux = cart.darcy_fluid_flux(1.0e6);
        assert!(
            flux > 0.0,
            "Darcy flux should be positive for positive pressure gradient"
        );
    }

    #[test]
    fn test_cartilage_donnan_pressure_positive() {
        let cart = CartilageModel::human_knee_cartilage();
        let pi = cart.donnan_swelling_pressure(0.15);
        assert!(pi > 0.0, "Donnan pressure should be positive, got {pi}");
    }

    #[test]
    fn test_cartilage_strain_dependent_permeability() {
        let cart = CartilageModel::human_knee_cartilage();
        // Negative strain (compression) should reduce permeability
        let k_comp = cart.strain_dependent_permeability(-0.2);
        let k_tens = cart.strain_dependent_permeability(0.2);
        assert!(
            k_comp < cart.permeability,
            "permeability should decrease in compression"
        );
        assert!(
            k_tens > cart.permeability,
            "permeability should increase in tension"
        );
    }

    // --- VascularWall ---

    #[test]
    fn test_vascular_wall_thickness() {
        let aw = VascularWall::human_aorta();
        let h = aw.wall_thickness();
        assert!(
            (h - 2.5e-3).abs() < 1e-6,
            "aorta wall thickness should be ~2.5 mm, got {h}"
        );
    }

    #[test]
    fn test_vascular_laplace_tension() {
        let aw = VascularWall::human_aorta();
        let p = 13_332.0; // ~100 mmHg in Pa
        let t = aw.laplace_tension(p);
        assert!(t > 0.0, "Laplace tension should be positive, got {t}");
    }

    #[test]
    fn test_vascular_circumferential_stress() {
        let aw = VascularWall::human_aorta();
        let s = aw.circumferential_stress(13_332.0);
        assert!(
            s > 0.0,
            "circumferential stress should be positive, got {s}"
        );
    }

    #[test]
    fn test_vascular_hgo_energy_at_identity() {
        let aw = VascularWall::human_aorta();
        let f_id = mat3_identity();
        let w = aw.hgo_strain_energy(&f_id);
        // At identity, both isochoric and fiber terms should be small
        assert!(
            w.abs() < 1.0,
            "HGO energy at identity should be near 0, got {w}"
        );
    }

    #[test]
    fn test_vascular_hgo_energy_positive_under_stretch() {
        let aw = VascularWall::human_aorta();
        let lambda = 1.2_f64;
        let lp = 1.0 / lambda.sqrt();
        let f = [[lambda, 0.0, 0.0], [0.0, lp, 0.0], [0.0, 0.0, lp]];
        let w = aw.hgo_strain_energy(&f);
        assert!(
            w > 0.0,
            "HGO energy should be positive under stretch, got {w}"
        );
    }

    #[test]
    fn test_vascular_residual_stretch() {
        let aw = VascularWall::human_aorta();
        let r_mid = 0.5 * (aw.inner_radius + aw.outer_radius);
        let lam = aw.residual_stretch_circumferential(r_mid);
        assert!(lam > 0.0, "residual stretch should be positive, got {lam}");
    }

    // --- CellMechanics ---

    #[test]
    fn test_cell_laplace_pressure() {
        let cell = CellMechanics::fibroblast();
        let p = cell.laplace_pressure();
        assert!(p > 0.0, "Laplace pressure should be positive, got {p}");
    }

    #[test]
    fn test_cell_effective_modulus() {
        let cell = CellMechanics::fibroblast();
        let e = cell.effective_modulus();
        assert!(
            e > cell.cytoskeletal_modulus,
            "effective modulus should exceed cytoskeletal modulus"
        );
    }

    #[test]
    fn test_cell_volume_positive() {
        let cell = CellMechanics::fibroblast();
        assert!(cell.volume() > 0.0, "cell volume should be positive");
    }

    #[test]
    fn test_cell_osmotic_pressure_sign() {
        let cell = CellMechanics::fibroblast();
        // Compression (negative strain) → positive osmotic pressure
        let p = cell.osmotic_pressure(-0.1);
        assert!(
            p > 0.0,
            "osmotic pressure under compression should be positive"
        );
    }

    #[test]
    fn test_cell_membrane_bending_stiffness() {
        let cell = CellMechanics::red_blood_cell();
        let kappa = cell.membrane_bending_stiffness();
        assert!(
            kappa > 0.0,
            "bending stiffness should be positive, got {kappa}"
        );
    }

    // --- HydrogelModel ---

    #[test]
    fn test_hydrogel_swelling_ratio() {
        let gel = HydrogelModel::peg_hydrogel();
        let q = gel.swelling_ratio();
        assert!(
            q > 1.0,
            "swelling ratio should be >1 for swollen gel, got {q}"
        );
    }

    #[test]
    fn test_hydrogel_linear_swelling_strain() {
        let gel = HydrogelModel::peg_hydrogel();
        let eps = gel.linear_swelling_strain();
        assert!(
            eps > 0.0,
            "linear swelling strain should be positive, got {eps}"
        );
    }

    #[test]
    fn test_hydrogel_network_modulus_positive() {
        let gel = HydrogelModel::peg_hydrogel();
        let g = gel.network_modulus();
        assert!(g > 0.0, "network modulus should be positive, got {g}");
    }

    #[test]
    fn test_hydrogel_effective_diffusivity() {
        let gel = HydrogelModel::hyaluronic_acid();
        let d0 = 1.0e-9; // m²/s (typical protein in water)
        let deff = gel.effective_diffusivity(d0);
        assert!(
            deff > 0.0,
            "effective diffusivity should be positive, got {deff}"
        );
        assert!(
            deff < d0,
            "effective diffusivity should be less than free-solution value"
        );
    }

    // --- BiomechanicsAnalysis ---

    #[test]
    fn test_stress_shielding_factor_range() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let f = ba.stress_shielding_factor();
        assert!(
            f > 0.0 && f < 1.0,
            "stress shielding factor should be in (0,1), got {f}"
        );
    }

    #[test]
    fn test_bone_stress_with_implant_less_than_without() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let with_imp = ba.bone_stress_with_implant();
        let without = ba.bone_stress_without_implant();
        assert!(
            with_imp < without,
            "bone stress with implant should be less (stress shielding): {with_imp} vs {without}"
        );
    }

    #[test]
    fn test_fatigue_life_decreases_with_stress() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let n1 = ba.fatigue_life(50.0e6);
        let n2 = ba.fatigue_life(100.0e6);
        assert!(
            n1 > n2,
            "fatigue life should decrease with higher stress: N1={n1}, N2={n2}"
        );
    }

    #[test]
    fn test_fatigue_life_infinite_at_zero_stress() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let n = ba.fatigue_life(0.0);
        assert!(
            n.is_infinite(),
            "fatigue life should be infinite at zero stress"
        );
    }

    #[test]
    fn test_miners_damage_accumulates() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let loading = vec![(1.0e5, 60.0e6), (2.0e5, 80.0e6)];
        let d = ba.miners_damage(&loading);
        assert!(d > 0.0, "Miner's damage should be positive, got {d}");
    }

    #[test]
    fn test_von_mises_stress_positive() {
        let vm =
            BiomechanicsAnalysis::von_mises_stress(100.0e6, 50.0e6, 20.0e6, 10.0e6, 5.0e6, 3.0e6);
        assert!(vm > 0.0, "Von Mises stress should be positive, got {vm}");
    }

    #[test]
    fn test_von_mises_stress_hydrostatic_zero() {
        // Under hydrostatic (equal principal stresses, no shear), VM = 0
        let vm = BiomechanicsAnalysis::von_mises_stress(100.0e6, 100.0e6, 100.0e6, 0.0, 0.0, 0.0);
        assert!(
            vm.abs() < 1.0,
            "VM stress under hydrostatic loading should be 0, got {vm}"
        );
    }

    #[test]
    fn test_tsai_wu_failure_below_one_safe() {
        let fi = BiomechanicsAnalysis::tsai_wu_failure_index(
            10.0e6, 5.0e6, 2.0e6, 120.0e6, 180.0e6, 60.0e6, 80.0e6, 50.0e6,
        );
        assert!(
            fi < 1.0,
            "Tsai-Wu index should be <1 for safe loading, got {fi}"
        );
    }

    #[test]
    fn test_tsai_wu_failure_above_one_failed() {
        let fi = BiomechanicsAnalysis::tsai_wu_failure_index(
            200.0e6, 0.0, 0.0, 120.0e6, 180.0e6, 60.0e6, 80.0e6, 50.0e6,
        );
        assert!(
            fi > 1.0,
            "Tsai-Wu index should be >1 for overloaded case, got {fi}"
        );
    }

    #[test]
    fn test_remodelling_stimulus_zero_in_dead_zone() {
        let stimulus = BiomechanicsAnalysis::remodelling_stimulus(1.0, 1.0, 0.1);
        assert!(
            stimulus.abs() < TOL,
            "remodelling stimulus should be 0 in dead zone"
        );
    }

    #[test]
    fn test_remodelling_stimulus_positive_above_reference() {
        let stimulus = BiomechanicsAnalysis::remodelling_stimulus(2.0, 1.0, 0.05);
        assert!(
            stimulus > 0.0,
            "remodelling stimulus should be positive above reference"
        );
    }

    #[test]
    fn test_max_principal_stress_failure_trigger() {
        let failed = BiomechanicsAnalysis {
            implant_modulus: 0.0,
            bone_modulus: 0.0,
            implant_area: 0.0,
            bone_area: 1.0,
            applied_load: 0.0,
            fatigue_exponent: 6.0,
            reference_fatigue_stress: 60.0e6,
            reference_cycles: 1.0e6,
        }
        .max_principal_stress_failure(200.0e6, 0.0, 0.0, 120.0e6, 180.0e6);
        assert!(
            failed,
            "should detect failure when stress exceeds tensile strength"
        );
    }

    #[test]
    fn test_max_principal_stress_no_failure() {
        let ba = BiomechanicsAnalysis::hip_replacement();
        let safe = ba.max_principal_stress_failure(50.0e6, 30.0e6, 10.0e6, 120.0e6, 180.0e6);
        assert!(!safe, "should be safe when stresses are below strength");
    }
}
