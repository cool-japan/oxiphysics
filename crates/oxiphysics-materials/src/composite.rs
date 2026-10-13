// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite material micromechanics.
//!
//! Provides Voigt/Reuss bounds, Hashin-Shtrikman bounds, Halpin-Tsai fiber-reinforced
//! composite model, Classical Laminate Theory (CLT), and Mori-Tanaka particle
//! composite scheme.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper: 3×3 matrix operations
// ---------------------------------------------------------------------------

fn mat3_zero() -> [[f64; 3]; 3] {
    [[0.0; 3]; 3]
}

fn mat3_add_scaled(acc: &mut [[f64; 3]; 3], m: &[[f64; 3]; 3], scale: f64) {
    for i in 0..3 {
        for j in 0..3 {
            acc[i][j] += m[i][j] * scale;
        }
    }
}

/// Invert a 3×3 matrix. Returns None if singular.
fn mat3_inv(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-300 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}

// ---------------------------------------------------------------------------
// IsotropicConstituent
// ---------------------------------------------------------------------------

/// An isotropic material constituent for composite micromechanics.
#[derive(Debug, Clone, Copy)]
pub struct IsotropicConstituent {
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Volume fraction (0..1).
    pub volume_fraction: f64,
}

impl IsotropicConstituent {
    /// Create a new isotropic constituent.
    pub fn new(
        youngs_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        volume_fraction: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
            density,
            volume_fraction,
        }
    }

    /// Shear modulus G = E / (2*(1+nu)).
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus K = E / (3*(1-2*nu)).
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
}

// ---------------------------------------------------------------------------
// Voigt-Reuss bounds
// ---------------------------------------------------------------------------

/// Voigt upper bound: E_V = Σ V_i · E_i.
pub fn voigt_modulus(constituents: &[IsotropicConstituent]) -> f64 {
    constituents
        .iter()
        .map(|c| c.volume_fraction * c.youngs_modulus)
        .sum()
}

/// Reuss lower bound: 1/E_R = Σ V_i / E_i.
pub fn reuss_modulus(constituents: &[IsotropicConstituent]) -> f64 {
    let inv_e: f64 = constituents
        .iter()
        .map(|c| c.volume_fraction / c.youngs_modulus)
        .sum();
    1.0 / inv_e
}

/// Voigt rule-of-mixtures density: ρ = Σ V_i · ρ_i.
pub fn voigt_density(constituents: &[IsotropicConstituent]) -> f64 {
    constituents
        .iter()
        .map(|c| c.volume_fraction * c.density)
        .sum()
}

/// Voigt rule-of-mixtures Poisson's ratio: ν = Σ V_i · ν_i.
pub fn voigt_poisson(constituents: &[IsotropicConstituent]) -> f64 {
    constituents
        .iter()
        .map(|c| c.volume_fraction * c.poisson_ratio)
        .sum()
}

/// Hashin-Shtrikman upper bound on bulk modulus K_HS+.
///
/// Uses the phase with the higher shear modulus as the reference (stiffer phase).
pub fn hashin_shtrikman_upper(matrix: &IsotropicConstituent, fiber: &IsotropicConstituent) -> f64 {
    // Upper bound: reference material is the stiffer phase (higher K and G)
    let (ref_phase, other) = if fiber.bulk_modulus() >= matrix.bulk_modulus() {
        (fiber, matrix)
    } else {
        (matrix, fiber)
    };
    let k_ref = ref_phase.bulk_modulus();
    let g_ref = ref_phase.shear_modulus();
    let k_other = other.bulk_modulus();
    let v_other = other.volume_fraction;

    // K_HS+ = K_ref + V_other / (1/(K_other - K_ref) + 3*V_ref/(3*K_ref + 4*G_ref))
    let v_ref = 1.0 - v_other;
    let denom = 1.0 / (k_other - k_ref) + 3.0 * v_ref / (3.0 * k_ref + 4.0 * g_ref);
    k_ref + v_other / denom
}

/// Hashin-Shtrikman lower bound on bulk modulus K_HS-.
///
/// Uses the phase with the lower shear modulus as the reference (softer phase).
pub fn hashin_shtrikman_lower(matrix: &IsotropicConstituent, fiber: &IsotropicConstituent) -> f64 {
    // Lower bound: reference material is the softer phase (lower K and G)
    let (ref_phase, other) = if matrix.bulk_modulus() <= fiber.bulk_modulus() {
        (matrix, fiber)
    } else {
        (fiber, matrix)
    };
    let k_ref = ref_phase.bulk_modulus();
    let g_ref = ref_phase.shear_modulus();
    let k_other = other.bulk_modulus();
    let v_other = other.volume_fraction;

    let v_ref = 1.0 - v_other;
    let denom = 1.0 / (k_other - k_ref) + 3.0 * v_ref / (3.0 * k_ref + 4.0 * g_ref);
    k_ref + v_other / denom
}

// ---------------------------------------------------------------------------
// Halpin-Tsai
// ---------------------------------------------------------------------------

/// Halpin-Tsai model for short-fiber or continuous-fiber reinforced composites.
#[derive(Debug, Clone, Copy)]
pub struct HalpinTsai {
    /// Fiber constituent.
    pub fiber: IsotropicConstituent,
    /// Matrix constituent.
    pub matrix: IsotropicConstituent,
    /// Fiber aspect ratio l/d.
    pub aspect_ratio: f64,
}

impl HalpinTsai {
    /// Create a new Halpin-Tsai model.
    pub fn new(
        fiber: IsotropicConstituent,
        matrix: IsotropicConstituent,
        aspect_ratio: f64,
    ) -> Self {
        Self {
            fiber,
            matrix,
            aspect_ratio,
        }
    }

    /// Reinforcing factor ξ for longitudinal direction = 2 * aspect_ratio.
    pub fn xi_longitudinal(&self) -> f64 {
        2.0 * self.aspect_ratio
    }

    /// Reinforcing factor ξ for transverse direction = 2.0.
    pub fn xi_transverse(&self) -> f64 {
        2.0
    }

    /// Halpin-Tsai parameter η = (E_f/E_m - 1) / (E_f/E_m + ξ).
    pub fn eta(property_ratio: f64, xi: f64) -> f64 {
        (property_ratio - 1.0) / (property_ratio + xi)
    }

    /// Longitudinal modulus E_11 using Halpin-Tsai equation.
    pub fn longitudinal_modulus(&self) -> f64 {
        let v_f = self.fiber.volume_fraction;
        let xi = self.xi_longitudinal();
        let ratio = self.fiber.youngs_modulus / self.matrix.youngs_modulus;
        let eta = Self::eta(ratio, xi);
        self.matrix.youngs_modulus * (1.0 + xi * eta * v_f) / (1.0 - eta * v_f)
    }

    /// Transverse modulus E_22 using Halpin-Tsai equation.
    pub fn transverse_modulus(&self) -> f64 {
        let v_f = self.fiber.volume_fraction;
        let xi = self.xi_transverse();
        let ratio = self.fiber.youngs_modulus / self.matrix.youngs_modulus;
        let eta = Self::eta(ratio, xi);
        self.matrix.youngs_modulus * (1.0 + xi * eta * v_f) / (1.0 - eta * v_f)
    }

    /// In-plane shear modulus G_12 using Halpin-Tsai equation (ξ = 1 for shear).
    pub fn shear_modulus_12(&self) -> f64 {
        let v_f = self.fiber.volume_fraction;
        let xi = 1.0;
        let ratio = self.fiber.shear_modulus() / self.matrix.shear_modulus();
        let eta = Self::eta(ratio, xi);
        self.matrix.shear_modulus() * (1.0 + xi * eta * v_f) / (1.0 - eta * v_f)
    }

    /// Major Poisson ratio ν_12 = V_f·ν_f + V_m·ν_m (rule of mixtures).
    pub fn poisson_12(&self) -> f64 {
        let v_f = self.fiber.volume_fraction;
        let v_m = self.matrix.volume_fraction;
        v_f * self.fiber.poisson_ratio + v_m * self.matrix.poisson_ratio
    }

    /// Minor Poisson ratio ν_21 = ν_12 · E_22 / E_11.
    pub fn poisson_21(&self) -> f64 {
        self.poisson_12() * self.transverse_modulus() / self.longitudinal_modulus()
    }
}

// ---------------------------------------------------------------------------
// OrthotropicPly
// ---------------------------------------------------------------------------

/// A single orthotropic ply (unidirectional lamina) in plane stress.
#[derive(Debug, Clone, Copy)]
pub struct OrthotropicPly {
    /// Longitudinal Young's modulus E_11 (fiber direction).
    pub e11: f64,
    /// Transverse Young's modulus E_22.
    pub e22: f64,
    /// In-plane shear modulus G_12.
    pub g12: f64,
    /// Major Poisson ratio ν_12.
    pub nu12: f64,
    /// Ply thickness.
    pub thickness: f64,
}

impl OrthotropicPly {
    /// Create a new orthotropic ply.
    pub fn new(e11: f64, e22: f64, g12: f64, nu12: f64, thickness: f64) -> Self {
        Self {
            e11,
            e22,
            g12,
            nu12,
            thickness,
        }
    }

    /// Minor Poisson ratio ν_21 = ν_12 · E_22 / E_11.
    pub fn nu21(&self) -> f64 {
        self.nu12 * self.e22 / self.e11
    }

    /// Reduced stiffness matrix \[Q\] in the principal material directions (plane stress).
    ///
    /// ```text
    /// Q11 = E11 / (1 - nu12*nu21)
    /// Q22 = E22 / (1 - nu12*nu21)
    /// Q12 = nu12*E22 / (1 - nu12*nu21)
    /// Q66 = G12
    /// ```
    pub fn q_matrix(&self) -> [[f64; 3]; 3] {
        let nu21 = self.nu21();
        let denom = 1.0 - self.nu12 * nu21;
        let q11 = self.e11 / denom;
        let q22 = self.e22 / denom;
        let q12 = self.nu12 * self.e22 / denom;
        let q66 = self.g12;
        [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
    }

    /// Transformed (off-axis) reduced stiffness matrix \[Q_bar\] for a ply at `angle_deg`.
    pub fn q_bar_matrix(&self, angle_deg: f64) -> [[f64; 3]; 3] {
        let theta = angle_deg * PI / 180.0;
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let c2s2 = c2 * s2;

        let q = self.q_matrix();
        let q11 = q[0][0];
        let q12 = q[0][1];
        let q22 = q[1][1];
        let q66 = q[2][2];

        let qb11 = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * s4;
        let qb12 = (q11 + q22 - 4.0 * q66) * c2s2 + q12 * (c4 + s4);
        let qb22 = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * c4;
        let qb16 = (q11 - q12 - 2.0 * q66) * c2 * c * s - (q22 - q12 - 2.0 * q66) * s2 * c * s;
        let qb26 = (q11 - q12 - 2.0 * q66) * s2 * c * s - (q22 - q12 - 2.0 * q66) * c2 * c * s;
        let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2s2 + q66 * (c4 + s4);

        [[qb11, qb12, qb16], [qb12, qb22, qb26], [qb16, qb26, qb66]]
    }
}

// ---------------------------------------------------------------------------
// Laminate (Classical Laminate Theory)
// ---------------------------------------------------------------------------

/// A composite laminate described by Classical Laminate Theory (CLT).
///
/// Plies are ordered from bottom to top; mid-plane is at z = 0.
#[derive(Debug, Clone)]
pub struct Laminate {
    /// Plies as (ply, fiber angle in degrees).
    pub plies: Vec<(OrthotropicPly, f64)>,
}

impl Laminate {
    /// Create an empty laminate.
    pub fn new() -> Self {
        Self { plies: Vec::new() }
    }

    /// Append a ply with the given fiber angle.
    pub fn add_ply(&mut self, ply: OrthotropicPly, angle_deg: f64) {
        self.plies.push((ply, angle_deg));
    }

    /// Total laminate thickness h = Σ t_k.
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|(p, _)| p.thickness).sum()
    }

    /// Z-coordinates of ply interfaces measured from the mid-plane.
    fn ply_z_coords(&self) -> Vec<f64> {
        let h = self.total_thickness();
        let mut z = Vec::with_capacity(self.plies.len() + 1);
        let mut current = -h / 2.0;
        z.push(current);
        for (ply, _) in &self.plies {
            current += ply.thickness;
            z.push(current);
        }
        z
    }

    /// Extensional stiffness matrix A_ij = Σ Q_bar_ij * h_k.
    pub fn a_matrix(&self) -> [[f64; 3]; 3] {
        let mut a = mat3_zero();
        for (ply, angle) in &self.plies {
            mat3_add_scaled(&mut a, &ply.q_bar_matrix(*angle), ply.thickness);
        }
        a
    }

    /// Coupling stiffness matrix B_ij = 0.5 * Σ Q_bar_ij * (z_k² - z_{k-1}²).
    pub fn b_matrix(&self) -> [[f64; 3]; 3] {
        let z = self.ply_z_coords();
        let mut b = mat3_zero();
        for (k, (ply, angle)) in self.plies.iter().enumerate() {
            let scale = 0.5 * (z[k + 1].powi(2) - z[k].powi(2));
            mat3_add_scaled(&mut b, &ply.q_bar_matrix(*angle), scale);
        }
        b
    }

    /// Bending stiffness matrix D_ij = (1/3) * Σ Q_bar_ij * (z_k³ - z_{k-1}³).
    pub fn d_matrix(&self) -> [[f64; 3]; 3] {
        let z = self.ply_z_coords();
        let mut d = mat3_zero();
        for (k, (ply, angle)) in self.plies.iter().enumerate() {
            let scale = (z[k + 1].powi(3) - z[k].powi(3)) / 3.0;
            mat3_add_scaled(&mut d, &ply.q_bar_matrix(*angle), scale);
        }
        d
    }

    /// Effective x-direction modulus: 1 / (A⁻¹\[0\]\[0\] * h).
    pub fn effective_modulus_x(&self) -> f64 {
        let h = self.total_thickness();
        if let Some(a_inv) = mat3_inv(&self.a_matrix()) {
            1.0 / (a_inv[0][0] * h)
        } else {
            0.0
        }
    }

    /// Effective y-direction modulus: 1 / (A⁻¹\[1\]\[1\] * h).
    pub fn effective_modulus_y(&self) -> f64 {
        let h = self.total_thickness();
        if let Some(a_inv) = mat3_inv(&self.a_matrix()) {
            1.0 / (a_inv[1][1] * h)
        } else {
            0.0
        }
    }

    /// Effective in-plane shear modulus: 1 / (A⁻¹\[2\]\[2\] * h).
    pub fn effective_shear_modulus(&self) -> f64 {
        let h = self.total_thickness();
        if let Some(a_inv) = mat3_inv(&self.a_matrix()) {
            1.0 / (a_inv[2][2] * h)
        } else {
            0.0
        }
    }

    /// Returns true if the laminate is symmetric (B ≈ 0).
    pub fn is_symmetric(&self) -> bool {
        let b = self.b_matrix();
        let a = self.a_matrix();
        // Scale tolerance by largest A entry
        let a_norm = a[0][0].abs().max(a[1][1].abs()).max(a[2][2].abs());
        let tol = 1.0e-6 * a_norm;
        b.iter().flatten().all(|v| v.abs() < tol)
    }
}

impl Default for Laminate {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Mori-Tanaka
// ---------------------------------------------------------------------------

/// Mori-Tanaka effective medium scheme for particle-reinforced composites.
///
/// Assumes spherical inclusions in an isotropic matrix.
#[derive(Debug, Clone, Copy)]
pub struct MoriTanaka {
    /// Matrix phase.
    pub matrix: IsotropicConstituent,
    /// Inclusion (particle) phase.
    pub inclusion: IsotropicConstituent,
}

impl MoriTanaka {
    /// Create a new Mori-Tanaka model.
    pub fn new(matrix: IsotropicConstituent, inclusion: IsotropicConstituent) -> Self {
        Self { matrix, inclusion }
    }

    /// Effective bulk modulus K* via Mori-Tanaka.
    ///
    /// K* = K_m + V_i*(K_i - K_m) / (1 + (1 - V_i)*(K_i - K_m)/(K_m + 4*G_m/3))
    pub fn effective_bulk_modulus(&self) -> f64 {
        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        let k_i = self.inclusion.bulk_modulus();
        let v_i = self.inclusion.volume_fraction;

        let alpha = k_m + 4.0 * g_m / 3.0;
        k_m + v_i * (k_i - k_m) / (1.0 + (1.0 - v_i) * (k_i - k_m) / alpha)
    }

    /// Effective shear modulus G* via Mori-Tanaka.
    ///
    /// G* = G_m + V_i*(G_i - G_m) / (1 + (1 - V_i)*(G_i - G_m)/(G_m + zeta_m))
    /// where zeta_m = G_m*(9*K_m + 8*G_m)/(6*(K_m + 2*G_m)) is the Eshelby
    /// constraint parameter of a spherical inclusion, so that
    /// G_m + zeta_m = 5*G_m*(3*K_m + 4*G_m)/(6*(K_m + 2*G_m)).
    pub fn effective_shear_modulus(&self) -> f64 {
        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        let g_i = self.inclusion.shear_modulus();
        let v_i = self.inclusion.volume_fraction;

        let zeta_m = g_m * (9.0 * k_m + 8.0 * g_m) / (6.0 * (k_m + 2.0 * g_m));
        g_m + v_i * (g_i - g_m) / (1.0 + (1.0 - v_i) * (g_i - g_m) / (g_m + zeta_m))
    }

    /// Effective Young's modulus E* derived from K* and G*.
    ///
    /// E = 9*K*G / (3*K + G)
    pub fn effective_youngs_modulus(&self) -> f64 {
        let k = self.effective_bulk_modulus();
        let g = self.effective_shear_modulus();
        9.0 * k * g / (3.0 * k + g)
    }
}

// ---------------------------------------------------------------------------
// Ply strengths
// ---------------------------------------------------------------------------

/// Strength properties for a unidirectional ply.
#[derive(Debug, Clone, Copy)]
pub struct PlyStrength {
    /// Longitudinal tensile strength (Pa).
    pub xt: f64,
    /// Longitudinal compressive strength (Pa, positive value).
    pub xc: f64,
    /// Transverse tensile strength (Pa).
    pub yt: f64,
    /// Transverse compressive strength (Pa, positive value).
    pub yc: f64,
    /// In-plane shear strength (Pa).
    pub s12: f64,
}

impl PlyStrength {
    /// Create new ply strength values.
    pub fn new(xt: f64, xc: f64, yt: f64, yc: f64, s12: f64) -> Self {
        Self {
            xt,
            xc,
            yt,
            yc,
            s12,
        }
    }

    /// Typical CFRP ply strength (T300/914).
    pub fn cfrp_typical() -> Self {
        Self {
            xt: 1500.0e6,
            xc: 1200.0e6,
            yt: 50.0e6,
            yc: 250.0e6,
            s12: 70.0e6,
        }
    }

    /// Typical GFRP ply strength.
    pub fn gfrp_typical() -> Self {
        Self {
            xt: 1000.0e6,
            xc: 600.0e6,
            yt: 30.0e6,
            yc: 110.0e6,
            s12: 55.0e6,
        }
    }
}

// ---------------------------------------------------------------------------
// Maximum stress failure criterion
// ---------------------------------------------------------------------------

/// Result of a failure check on a single ply.
#[derive(Debug, Clone, Copy)]
pub struct PlyFailureResult {
    /// Fiber direction failure index (>= 1.0 means failure).
    pub fi_fiber: f64,
    /// Matrix direction failure index.
    pub fi_matrix: f64,
    /// Shear failure index.
    pub fi_shear: f64,
    /// Whether any mode has failed.
    pub failed: bool,
}

/// Maximum stress failure criterion.
///
/// Checks stress in principal material directions against allowable strengths.
pub fn max_stress_failure(sigma: [f64; 3], strength: &PlyStrength) -> PlyFailureResult {
    let s1 = sigma[0];
    let s2 = sigma[1];
    let s12 = sigma[2];

    let fi_fiber = if s1 >= 0.0 {
        s1 / strength.xt
    } else {
        (-s1) / strength.xc
    };

    let fi_matrix = if s2 >= 0.0 {
        s2 / strength.yt
    } else {
        (-s2) / strength.yc
    };

    let fi_shear = s12.abs() / strength.s12;

    PlyFailureResult {
        fi_fiber,
        fi_matrix,
        fi_shear,
        failed: fi_fiber >= 1.0 || fi_matrix >= 1.0 || fi_shear >= 1.0,
    }
}

// ---------------------------------------------------------------------------
// Tsai-Wu failure criterion
// ---------------------------------------------------------------------------

/// Tsai-Wu failure criterion.
///
/// Returns the failure index (>= 1.0 means failure).
pub fn tsai_wu_failure(sigma: [f64; 3], strength: &PlyStrength) -> f64 {
    let s1 = sigma[0];
    let s2 = sigma[1];
    let s12 = sigma[2];

    let f1 = 1.0 / strength.xt - 1.0 / strength.xc;
    let f2 = 1.0 / strength.yt - 1.0 / strength.yc;
    let f11 = 1.0 / (strength.xt * strength.xc);
    let f22 = 1.0 / (strength.yt * strength.yc);
    let f66 = 1.0 / (strength.s12 * strength.s12);
    // Interaction term: conservative estimate
    let f12 = -0.5 * (f11 * f22).sqrt();

    f1 * s1 + f2 * s2 + f11 * s1 * s1 + f22 * s2 * s2 + f66 * s12 * s12 + 2.0 * f12 * s1 * s2
}

// ---------------------------------------------------------------------------
// Hashin failure criterion
// ---------------------------------------------------------------------------

/// Hashin failure result with separate modes.
#[derive(Debug, Clone, Copy)]
pub struct HashinResult {
    /// Fiber tension failure index.
    pub fiber_tension: f64,
    /// Fiber compression failure index.
    pub fiber_compression: f64,
    /// Matrix tension failure index.
    pub matrix_tension: f64,
    /// Matrix compression failure index.
    pub matrix_compression: f64,
    /// Whether any mode has failed.
    pub failed: bool,
}

/// Hashin failure criterion (1980 formulation).
///
/// Provides distinct failure modes for fiber and matrix in tension/compression.
pub fn hashin_failure(sigma: [f64; 3], strength: &PlyStrength) -> HashinResult {
    let s1 = sigma[0];
    let s2 = sigma[1];
    let s12 = sigma[2];

    // Fiber tension (s1 >= 0)
    let fiber_tension = if s1 >= 0.0 {
        (s1 / strength.xt).powi(2) + (s12 / strength.s12).powi(2)
    } else {
        0.0
    };

    // Fiber compression (s1 < 0)
    let fiber_compression = if s1 < 0.0 {
        (-s1 / strength.xc).powi(2)
    } else {
        0.0
    };

    // Matrix tension (s2 >= 0)
    let matrix_tension = if s2 >= 0.0 {
        (s2 / strength.yt).powi(2) + (s12 / strength.s12).powi(2)
    } else {
        0.0
    };

    // Matrix compression (s2 < 0)
    let matrix_compression = if s2 < 0.0 {
        (s2 / (2.0 * strength.s12)).powi(2)
            + ((strength.yc / (2.0 * strength.s12)).powi(2) - 1.0) * s2 / strength.yc
            + (s12 / strength.s12).powi(2)
    } else {
        0.0
    };

    let failed = fiber_tension >= 1.0
        || fiber_compression >= 1.0
        || matrix_tension >= 1.0
        || matrix_compression >= 1.0;

    HashinResult {
        fiber_tension,
        fiber_compression,
        matrix_tension,
        matrix_compression,
        failed,
    }
}

// ---------------------------------------------------------------------------
// Puck failure criterion (plane stress)
// ---------------------------------------------------------------------------

/// Puck inter-fiber failure result.
#[derive(Debug, Clone, Copy)]
pub struct PuckResult {
    /// Inter-fiber failure effort (>= 1.0 means failure).
    pub iff_effort: f64,
    /// Fiber failure effort.
    pub ff_effort: f64,
    /// Whether any mode has failed.
    pub failed: bool,
}

/// Default Puck inclination parameter p⊥∥⁻ (dimensionless) of the
/// inter-fiber failure envelope under transverse compression.  The value 0.25
/// is the usual recommendation for carbon- and glass-fibre/epoxy plies
/// (Puck & Schürmann 2002, VDI 2014).
const PUCK_P_PERP_PAR_MINUS: f64 = 0.25;

/// Puck inter-fiber failure effort for transverse compression (`σ₂ < 0`) in
/// the plane-stress formulation of Puck & Schürmann.
///
/// Arguments are the transverse stress `sigma2` (negative), the in-plane
/// shear stress `tau21`, the transverse compressive strength `yc` (positive),
/// the in-plane shear strength `s21` and the inclination parameter
/// `p_perp_par_minus` (p⊥∥⁻).
///
/// The criterion distinguishes two fracture modes, separated by the stress
/// ratio `|σ₂/τ₂₁| = R⊥⊥ᴬ / τ₂₁c`:
///
/// * **Mode B** (shear dominated; the fracture plane is perpendicular to the
///   ply plane, θ_fp = 0°):
///   `f = (sqrt(τ₂₁² + (p⊥∥⁻·σ₂)²) + p⊥∥⁻·σ₂) / S₂₁`
/// * **Mode C** (compression dominated; oblique fracture plane, θ_fp > 0°):
///   `f = [(τ₂₁ / (2·(1 + p⊥⊥⁻)·S₂₁))² + (σ₂/Yc)²] · Yc / (−σ₂)`
///
/// with the fracture resistance of the action plane
/// `R⊥⊥ᴬ = S₂₁/(2·p⊥∥⁻) · (sqrt(1 + 2·p⊥∥⁻·Yc/S₂₁) − 1)`, the derived
/// parameter `p⊥⊥⁻ = p⊥∥⁻·R⊥⊥ᴬ/S₂₁` and the shear stress at the mode
/// transition `τ₂₁c = S₂₁·sqrt(1 + 2·p⊥⊥⁻)`.  `R⊥⊥ᴬ` is evaluated in the
/// algebraically identical, cancellation-free form
/// `Yc / (1 + sqrt(1 + 2·p⊥∥⁻·Yc/S₂₁))`, which stays finite for `p⊥∥⁻ → 0`.
///
/// The two branches meet continuously at the mode boundary, and pure
/// transverse compression `σ₂ = -Yc, τ₂₁ = 0` gives exactly `f = 1`.
pub(crate) fn puck_iff_transverse_compression(
    sigma2: f64,
    tau21: f64,
    yc: f64,
    s21: f64,
    p_perp_par_minus: f64,
) -> f64 {
    let p = p_perp_par_minus;
    let tau = tau21.abs();
    let compression = -sigma2;

    let r_perp_perp_a = yc / (1.0 + (1.0 + 2.0 * p * yc / s21).sqrt());
    let p_perp_perp_minus = p * r_perp_perp_a / s21;
    let tau_21c = s21 * (1.0 + 2.0 * p_perp_perp_minus).sqrt();

    // |σ₂/τ₂₁| <= R⊥⊥ᴬ/τ₂₁c written without a division so that τ₂₁ = 0
    // (pure compression) falls into Mode C.
    if compression * tau_21c <= r_perp_perp_a * tau {
        // Mode B
        let friction = p * sigma2;
        ((tau * tau + friction * friction).sqrt() + friction) / s21
    } else {
        // Mode C
        let shear_term = tau / (2.0 * (1.0 + p_perp_perp_minus) * s21);
        (shear_term * shear_term + (sigma2 / yc).powi(2)) * yc / compression
    }
}

/// Plane-stress Puck failure criterion.
///
/// Fiber failure (FF) is a maximum-stress check in the fiber direction.
/// Inter-fiber failure (IFF) follows the action-plane concept: transverse
/// tension (`σ₂ >= 0`, Mode A) uses the elliptical interaction of `σ₂/Yt`
/// and `τ₁₂/S₁₂`, while transverse compression (`σ₂ < 0`) uses the Puck
/// Mode B / Mode C formulation with the default inclination parameter
/// `p⊥∥⁻ = 0.25`, so that `σ₂ = -Yc` gives an effort of exactly 1.
pub fn puck_failure(sigma: [f64; 3], strength: &PlyStrength) -> PuckResult {
    let s1 = sigma[0];
    let s2 = sigma[1];
    let s12 = sigma[2];

    // Fiber failure
    let ff_effort = if s1 >= 0.0 {
        s1 / strength.xt
    } else {
        (-s1) / strength.xc
    };

    // Inter-fiber failure
    let iff_effort = if s2 >= 0.0 {
        // Mode A: transverse tension
        ((s2 / strength.yt).powi(2) + (s12 / strength.s12).powi(2)).sqrt()
    } else {
        // Mode B / Mode C: transverse compression
        puck_iff_transverse_compression(s2, s12, strength.yc, strength.s12, PUCK_P_PERP_PAR_MINUS)
    };

    PuckResult {
        iff_effort,
        ff_effort,
        failed: iff_effort >= 1.0 || ff_effort >= 1.0,
    }
}

// ---------------------------------------------------------------------------
// Ply-by-ply stress analysis
// ---------------------------------------------------------------------------

/// Stress state in a single ply (in material coordinates).
#[derive(Debug, Clone, Copy)]
pub struct PlyStress {
    /// Ply index.
    pub ply_index: usize,
    /// Fiber angle (degrees).
    pub angle_deg: f64,
    /// Stress in material coordinates \[sigma_1, sigma_2, tau_12\].
    pub stress: [f64; 3],
    /// Strain in material coordinates \[eps_1, eps_2, gamma_12\].
    pub strain: [f64; 3],
    /// Z-coordinate of ply midplane.
    pub z_mid: f64,
}

/// Transform global laminate strains to ply material coordinates.
///
/// `global_strain` is \[eps_x, eps_y, gamma_xy\] at the ply midplane.
/// `angle_deg` is the fiber orientation.
fn transform_strain_to_material(global_strain: [f64; 3], angle_deg: f64) -> [f64; 3] {
    let theta = angle_deg * PI / 180.0;
    let c = theta.cos();
    let s = theta.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;

    let ex = global_strain[0];
    let ey = global_strain[1];
    let gxy = global_strain[2];

    [
        c2 * ex + s2 * ey + cs * gxy,
        s2 * ex + c2 * ey - cs * gxy,
        -2.0 * cs * ex + 2.0 * cs * ey + (c2 - s2) * gxy,
    ]
}

/// Compute ply-by-ply stresses for a laminate under given midplane strains.
///
/// `midplane_strain` is \[eps_x0, eps_y0, gamma_xy0\].
/// `curvature` is \[kappa_x, kappa_y, kappa_xy\].
pub fn ply_by_ply_stress(
    laminate: &Laminate,
    midplane_strain: [f64; 3],
    curvature: [f64; 3],
) -> Vec<PlyStress> {
    let z_coords = laminate.ply_z_coords();
    let mut results = Vec::with_capacity(laminate.plies.len());

    for (k, (ply, angle)) in laminate.plies.iter().enumerate() {
        let z_mid = (z_coords[k] + z_coords[k + 1]) / 2.0;

        // Global strain at ply midplane: eps = eps0 + z * kappa
        let global_strain = [
            midplane_strain[0] + z_mid * curvature[0],
            midplane_strain[1] + z_mid * curvature[1],
            midplane_strain[2] + z_mid * curvature[2],
        ];

        // Transform to material coordinates
        let mat_strain = transform_strain_to_material(global_strain, *angle);

        // Compute stress using Q matrix: sigma = Q * eps (in material coords)
        let q = ply.q_matrix();
        let stress = [
            q[0][0] * mat_strain[0] + q[0][1] * mat_strain[1],
            q[0][1] * mat_strain[0] + q[1][1] * mat_strain[1],
            q[2][2] * mat_strain[2],
        ];

        results.push(PlyStress {
            ply_index: k,
            angle_deg: *angle,
            stress,
            strain: mat_strain,
            z_mid,
        });
    }

    results
}

// ---------------------------------------------------------------------------
// Progressive failure analysis
// ---------------------------------------------------------------------------

/// Result of progressive failure analysis.
#[derive(Debug, Clone)]
pub struct ProgressiveFailureResult {
    /// Load factor at first ply failure.
    pub first_ply_failure_load: f64,
    /// Index of the first ply to fail.
    pub first_failed_ply: usize,
    /// Load factor at last ply failure (ultimate).
    pub ultimate_load: f64,
    /// Number of plies that failed.
    pub n_failed_plies: usize,
}

/// Run progressive failure analysis using max-stress criterion.
///
/// Applies increasing load factor to the given strain until all plies fail.
/// Returns the load factor at first ply failure and ultimate failure.
pub fn progressive_failure_analysis(
    laminate: &Laminate,
    base_strain: [f64; 3],
    strength: &PlyStrength,
    max_load_factor: f64,
    n_steps: usize,
) -> ProgressiveFailureResult {
    let n_plies = laminate.plies.len();
    let mut failed = vec![false; n_plies];
    let mut first_failure_load = max_load_factor;
    let mut first_failed_ply = 0;
    let mut ultimate_load = max_load_factor;
    let mut n_failed = 0;
    let mut first_found = false;

    let curvature = [0.0, 0.0, 0.0];

    for step in 0..=n_steps {
        let lf = max_load_factor * (step as f64) / (n_steps as f64);
        let strain = [
            base_strain[0] * lf,
            base_strain[1] * lf,
            base_strain[2] * lf,
        ];

        let ply_stresses = ply_by_ply_stress(laminate, strain, curvature);

        let mut all_failed = true;
        for ps in &ply_stresses {
            if failed[ps.ply_index] {
                continue;
            }
            let result = max_stress_failure(ps.stress, strength);
            if result.failed {
                failed[ps.ply_index] = true;
                n_failed += 1;
                if !first_found {
                    first_failure_load = lf;
                    first_failed_ply = ps.ply_index;
                    first_found = true;
                }
            }
            if !failed[ps.ply_index] {
                all_failed = false;
            }
        }

        if all_failed && n_plies > 0 {
            ultimate_load = lf;
            break;
        }
    }

    ProgressiveFailureResult {
        first_ply_failure_load: first_failure_load,
        first_failed_ply,
        ultimate_load,
        n_failed_plies: n_failed,
    }
}

// ---------------------------------------------------------------------------
// Thermal residual stresses
// ---------------------------------------------------------------------------

/// Compute thermal residual strains in a ply due to temperature change.
///
/// `alpha` is \[alpha_1, alpha_2\] CTE in material coords.
/// `delta_t` is the temperature change.
pub fn thermal_strain(alpha: [f64; 2], delta_t: f64) -> [f64; 3] {
    [alpha[0] * delta_t, alpha[1] * delta_t, 0.0]
}

/// Compute thermal residual stresses in each ply of a laminate.
///
/// Assumes the laminate is constrained (cured as a whole) and computes
/// the residual stress due to CTE mismatch between plies.
///
/// `alpha` is \[alpha_1, alpha_2\] CTE for the ply material.
/// `delta_t` is the cure temperature change (e.g., from cure temp to room temp).
pub fn thermal_residual_stresses(
    laminate: &Laminate,
    alpha: [f64; 2],
    delta_t: f64,
) -> Vec<[f64; 3]> {
    // Compute effective laminate CTE (using A matrix)
    let a_mat = laminate.a_matrix();
    let a_inv = match mat3_inv(&a_mat) {
        Some(inv) => inv,
        None => return vec![[0.0; 3]; laminate.plies.len()],
    };

    // Compute thermal force resultant: N_T = sum(Q_bar * alpha_bar * dT * thickness)
    let mut n_thermal = [0.0; 3];
    for (ply, angle) in &laminate.plies {
        let theta = angle * PI / 180.0;
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;

        // Transform CTE to global coords
        let alpha_x = c2 * alpha[0] + s2 * alpha[1];
        let alpha_y = s2 * alpha[0] + c2 * alpha[1];
        let alpha_xy = 2.0 * cs * (alpha[0] - alpha[1]);

        let qbar = ply.q_bar_matrix(*angle);
        let t = ply.thickness;

        n_thermal[0] +=
            (qbar[0][0] * alpha_x + qbar[0][1] * alpha_y + qbar[0][2] * alpha_xy) * delta_t * t;
        n_thermal[1] +=
            (qbar[1][0] * alpha_x + qbar[1][1] * alpha_y + qbar[1][2] * alpha_xy) * delta_t * t;
        n_thermal[2] +=
            (qbar[2][0] * alpha_x + qbar[2][1] * alpha_y + qbar[2][2] * alpha_xy) * delta_t * t;
    }

    // Midplane strains due to thermal loads: eps0 = A^-1 * N_T
    let eps0 = [
        a_inv[0][0] * n_thermal[0] + a_inv[0][1] * n_thermal[1] + a_inv[0][2] * n_thermal[2],
        a_inv[1][0] * n_thermal[0] + a_inv[1][1] * n_thermal[1] + a_inv[1][2] * n_thermal[2],
        a_inv[2][0] * n_thermal[0] + a_inv[2][1] * n_thermal[1] + a_inv[2][2] * n_thermal[2],
    ];

    // Compute residual stress in each ply
    let mut results = Vec::with_capacity(laminate.plies.len());
    for (ply, angle) in &laminate.plies {
        let theta = angle * PI / 180.0;
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;

        // Free thermal strain in global coords
        let alpha_x = c2 * alpha[0] + s2 * alpha[1];
        let alpha_y = s2 * alpha[0] + c2 * alpha[1];
        let alpha_xy = 2.0 * cs * (alpha[0] - alpha[1]);

        // Mechanical strain = total strain - free thermal strain
        let mech_strain_global = [
            eps0[0] - alpha_x * delta_t,
            eps0[1] - alpha_y * delta_t,
            eps0[2] - alpha_xy * delta_t,
        ];

        // Transform to material coords
        let mat_strain = transform_strain_to_material(mech_strain_global, *angle);

        // Compute stress
        let q = ply.q_matrix();
        let stress = [
            q[0][0] * mat_strain[0] + q[0][1] * mat_strain[1],
            q[0][1] * mat_strain[0] + q[1][1] * mat_strain[1],
            q[2][2] * mat_strain[2],
        ];

        results.push(stress);
    }

    results
}

// ---------------------------------------------------------------------------
// Interlaminar stresses
// ---------------------------------------------------------------------------

mod interlaminar;

pub use interlaminar::interlaminar_shear_estimate;
pub(crate) use interlaminar::solve_linear_6;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-6;

    #[test]
    fn test_voigt_reuss_ordering() {
        let c1 = IsotropicConstituent::new(70.0e9, 0.33, 2700.0, 0.3);
        let c2 = IsotropicConstituent::new(400.0e9, 0.25, 3200.0, 0.7);
        let constituents = [c1, c2];
        let ev = voigt_modulus(&constituents);
        let er = reuss_modulus(&constituents);
        assert!(ev > er, "Voigt ({ev}) should be > Reuss ({er})");
    }

    #[test]
    fn test_voigt_density_weighted_average() {
        let c1 = IsotropicConstituent::new(70.0e9, 0.33, 2700.0, 0.4);
        let c2 = IsotropicConstituent::new(200.0e9, 0.30, 7800.0, 0.6);
        let constituents = [c1, c2];
        let rho = voigt_density(&constituents);
        let expected = 0.4 * 2700.0 + 0.6 * 7800.0;
        assert!(
            (rho - expected).abs() < TOL,
            "density = {rho} expected {expected}"
        );
    }

    #[test]
    fn test_halpin_tsai_longitudinal_between_fiber_matrix() {
        let fiber = IsotropicConstituent::new(230.0e9, 0.25, 1800.0, 0.6);
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 0.4);
        let ht = HalpinTsai::new(fiber, matrix, 50.0);
        let e11 = ht.longitudinal_modulus();
        assert!(
            e11 > matrix.youngs_modulus && e11 < fiber.youngs_modulus,
            "E11 = {e11} not between matrix ({}) and fiber ({})",
            matrix.youngs_modulus,
            fiber.youngs_modulus
        );
    }

    #[test]
    fn test_orthotropic_ply_q_matrix_symmetry() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let q = ply.q_matrix();
        assert!(
            (q[0][1] - q[1][0]).abs() < TOL,
            "Q[0][1]={} != Q[1][0]={}",
            q[0][1],
            q[1][0]
        );
    }

    #[test]
    fn test_symmetric_crossply_b_near_zero() {
        let t = 0.125e-3_f64;
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, t);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);
        lam.add_ply(ply, 0.0);
        assert!(
            lam.is_symmetric(),
            "Symmetric [0/90/0] laminate should have B ~ 0"
        );
    }

    #[test]
    fn test_laminate_single_0deg_effective_modulus_x() {
        let e11 = 140.0e9_f64;
        let e22 = 10.0e9_f64;
        let g12 = 5.0e9_f64;
        let nu12 = 0.3_f64;
        let t = 1.0e-3_f64;
        let ply = OrthotropicPly::new(e11, e22, g12, nu12, t);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        let ex = lam.effective_modulus_x();
        assert!(
            (ex - e11).abs() / e11 < 0.05,
            "effective_modulus_x = {ex}, E11 = {e11}"
        );
    }

    #[test]
    fn test_mori_tanaka_modulus_between_phases() {
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 0.7);
        let inclusion = IsotropicConstituent::new(70.0e9, 0.22, 2700.0, 0.3);
        let mt = MoriTanaka::new(matrix, inclusion);
        let e_eff = mt.effective_youngs_modulus();
        assert!(
            e_eff > matrix.youngs_modulus && e_eff < inclusion.youngs_modulus,
            "E_eff = {e_eff} not between matrix ({}) and inclusion ({})",
            matrix.youngs_modulus,
            inclusion.youngs_modulus
        );
    }

    /// Glass spheres in epoxy at V = 0.3: the shear modulus must agree with the
    /// independent closed form in `combination::mori_tanaka_homogenization`.
    #[test]
    fn test_mori_tanaka_shear_matches_combination_module() {
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 0.7);
        let inclusion = IsotropicConstituent::new(72.0e9, 0.22, 2500.0, 0.3);
        let mt = MoriTanaka::new(matrix, inclusion);
        let (k_ref, g_ref) = crate::combination::mori_tanaka_homogenization(
            0.3,
            inclusion.bulk_modulus(),
            inclusion.shear_modulus(),
            matrix.bulk_modulus(),
            matrix.shear_modulus(),
        );
        let g_eff = mt.effective_shear_modulus();
        assert!(
            (g_eff - g_ref).abs() / g_ref < 1e-9,
            "G* = {g_eff}, reference = {g_ref}"
        );
        let k_eff = mt.effective_bulk_modulus();
        assert!(
            (k_eff - k_ref).abs() / k_ref < 1e-9,
            "K* = {k_eff}, reference = {k_ref}"
        );
    }

    /// Mori-Tanaka limits: no inclusions returns the matrix, a pure inclusion
    /// phase returns the inclusion, and a dilute suspension follows the
    /// Eshelby dilute-limit slope dG/dV = (G_i - G_m)/(1 + (G_i - G_m)/(G_m + zeta_m)).
    #[test]
    fn test_mori_tanaka_shear_limits() {
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 1.0);
        let g_m = matrix.shear_modulus();
        let k_m = matrix.bulk_modulus();
        let inc = |v: f64| IsotropicConstituent::new(72.0e9, 0.22, 2500.0, v);

        let empty = MoriTanaka::new(matrix, inc(0.0)).effective_shear_modulus();
        assert!((empty - g_m).abs() / g_m < 1e-12, "V=0: {empty} vs {g_m}");

        let full = MoriTanaka::new(matrix, inc(1.0)).effective_shear_modulus();
        let g_i = inc(1.0).shear_modulus();
        assert!((full - g_i).abs() / g_i < 1e-12, "V=1: {full} vs {g_i}");

        let v = 1.0e-6;
        let dilute = MoriTanaka::new(matrix, inc(v)).effective_shear_modulus();
        let zeta_m = g_m * (9.0 * k_m + 8.0 * g_m) / (6.0 * (k_m + 2.0 * g_m));
        let slope = (g_i - g_m) / (1.0 + (g_i - g_m) / (g_m + zeta_m));
        let expected = g_m + v * slope;
        assert!(
            (dilute - expected).abs() / g_m < 1e-9,
            "dilute limit: {dilute} vs {expected}"
        );
    }

    /// Effective shear modulus lies between the Reuss and Voigt bounds and
    /// grows monotonically with the inclusion volume fraction.
    #[test]
    fn test_mori_tanaka_shear_bounded_and_monotone() {
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 1.0);
        let g_m = matrix.shear_modulus();
        let g_i = IsotropicConstituent::new(72.0e9, 0.22, 2500.0, 1.0).shear_modulus();
        let mut previous = g_m;
        for step in 1..=9 {
            let v = f64::from(step) / 10.0;
            let inclusion = IsotropicConstituent::new(72.0e9, 0.22, 2500.0, v);
            let g = MoriTanaka::new(matrix, inclusion).effective_shear_modulus();
            let voigt = (1.0 - v) * g_m + v * g_i;
            let reuss = 1.0 / ((1.0 - v) / g_m + v / g_i);
            assert!(
                g > reuss && g < voigt,
                "V={v}: Reuss {reuss} < G* {g} < Voigt {voigt}"
            );
            assert!(g > previous, "V={v}: not monotone ({g} <= {previous})");
            previous = g;
        }
    }

    #[test]
    fn test_halpin_tsai_poisson_12() {
        let fiber = IsotropicConstituent::new(230.0e9, 0.25, 1800.0, 0.6);
        let matrix = IsotropicConstituent::new(3.5e9, 0.35, 1200.0, 0.4);
        let ht = HalpinTsai::new(fiber, matrix, 50.0);
        let nu12 = ht.poisson_12();
        let expected = 0.6 * 0.25 + 0.4 * 0.35;
        assert!(
            (nu12 - expected).abs() < TOL,
            "nu_12 = {nu12} expected {expected}"
        );
    }

    // --- Max stress failure ---

    #[test]
    fn test_max_stress_no_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [100.0e6, 10.0e6, 5.0e6];
        let r = max_stress_failure(sigma, &s);
        assert!(!r.failed);
        assert!(r.fi_fiber < 1.0);
        assert!(r.fi_matrix < 1.0);
        assert!(r.fi_shear < 1.0);
    }

    #[test]
    fn test_max_stress_fiber_tension_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [2000.0e6, 0.0, 0.0]; // exceeds xt=1500 MPa
        let r = max_stress_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fi_fiber >= 1.0);
    }

    #[test]
    fn test_max_stress_fiber_compression_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [-1500.0e6, 0.0, 0.0]; // exceeds xc=1200 MPa
        let r = max_stress_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fi_fiber >= 1.0);
    }

    #[test]
    fn test_max_stress_matrix_tension_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 60.0e6, 0.0]; // exceeds yt=50 MPa
        let r = max_stress_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fi_matrix >= 1.0);
    }

    #[test]
    fn test_max_stress_shear_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 0.0, 80.0e6]; // exceeds s12=70 MPa
        let r = max_stress_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fi_shear >= 1.0);
    }

    // --- Tsai-Wu failure ---

    #[test]
    fn test_tsai_wu_no_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [100.0e6, 5.0e6, 5.0e6];
        let fi = tsai_wu_failure(sigma, &s);
        assert!(fi < 1.0, "Tsai-Wu FI = {fi} should be < 1");
    }

    #[test]
    fn test_tsai_wu_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 60.0e6, 0.0]; // matrix tension failure
        let fi = tsai_wu_failure(sigma, &s);
        assert!(fi >= 1.0, "Tsai-Wu FI = {fi} should be >= 1");
    }

    #[test]
    fn test_tsai_wu_zero_stress() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 0.0, 0.0];
        let fi = tsai_wu_failure(sigma, &s);
        assert!(
            fi.abs() < 1e-10,
            "zero stress should give zero FI, got {fi}"
        );
    }

    // --- Hashin failure ---

    #[test]
    fn test_hashin_no_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [100.0e6, 5.0e6, 5.0e6];
        let r = hashin_failure(sigma, &s);
        assert!(!r.failed);
    }

    #[test]
    fn test_hashin_fiber_tension() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [2000.0e6, 0.0, 0.0];
        let r = hashin_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fiber_tension >= 1.0);
    }

    #[test]
    fn test_hashin_fiber_compression() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [-1500.0e6, 0.0, 0.0];
        let r = hashin_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.fiber_compression >= 1.0);
    }

    #[test]
    fn test_hashin_matrix_tension() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 60.0e6, 0.0];
        let r = hashin_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.matrix_tension >= 1.0);
    }

    #[test]
    fn test_hashin_matrix_compression() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, -300.0e6, 0.0];
        let r = hashin_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.matrix_compression >= 1.0);
    }

    // --- Puck failure ---

    #[test]
    fn test_puck_no_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [100.0e6, 5.0e6, 5.0e6];
        let r = puck_failure(sigma, &s);
        assert!(!r.failed);
    }

    #[test]
    fn test_puck_fiber_failure() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [2000.0e6, 0.0, 0.0];
        let r = puck_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.ff_effort >= 1.0);
    }

    #[test]
    fn test_puck_iff_tension() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, 60.0e6, 0.0];
        let r = puck_failure(sigma, &s);
        assert!(r.failed);
        assert!(r.iff_effort >= 1.0);
    }

    #[test]
    fn test_puck_iff_compression() {
        let s = PlyStrength::cfrp_typical();
        let sigma = [0.0, -300.0e6, 50.0e6];
        let r = puck_failure(sigma, &s);
        // 300 MPa exceeds the 250 MPa transverse compressive strength, so the
        // Mode C effort must exceed 1 and flag failure.
        assert!(r.iff_effort > 1.0, "iff_effort = {}", r.iff_effort);
        assert!(r.failed);
    }

    /// Textbook transcription of the Puck plane-stress IFF compression
    /// branch, deliberately written with the un-rearranged formulas and an
    /// explicit stress-ratio test so it is independent of the implementation.
    fn puck_reference_compression(sigma2: f64, tau21: f64, yc: f64, s21: f64, p: f64) -> f64 {
        let r_a = s21 / (2.0 * p) * ((1.0 + 2.0 * p * yc / s21).sqrt() - 1.0);
        let p_pp = p * r_a / s21;
        let tau_c = s21 * (1.0 + 2.0 * p_pp).sqrt();
        if (sigma2 / tau21).abs() <= r_a / tau_c {
            ((tau21.powi(2) + (p * sigma2).powi(2)).sqrt() + p * sigma2) / s21
        } else {
            ((tau21 / (2.0 * (1.0 + p_pp) * s21)).powi(2) + (sigma2 / yc).powi(2)) * yc / (-sigma2)
        }
    }

    #[test]
    fn test_puck_iff_pure_transverse_compression_reaches_one_at_yc() {
        for s in [PlyStrength::cfrp_typical(), PlyStrength::gfrp_typical()] {
            let r = puck_failure([0.0, -s.yc, 0.0], &s);
            assert!(
                (r.iff_effort - 1.0).abs() < 1e-12,
                "iff_effort at sigma2 = -Yc must be 1, got {}",
                r.iff_effort
            );
            assert!(r.failed, "transverse compressive failure must be flagged");
        }
    }

    #[test]
    fn test_puck_iff_pure_transverse_compression_is_linear_in_stress() {
        let s = PlyStrength::cfrp_typical();
        for &fraction in &[0.1, 0.5, 0.9, 1.2, 2.0] {
            let r = puck_failure([0.0, -fraction * s.yc, 0.0], &s);
            assert!(
                (r.iff_effort - fraction).abs() < 1e-12,
                "fraction {fraction}: effort {}",
                r.iff_effort
            );
        }
    }

    #[test]
    fn test_puck_iff_compression_matches_textbook_reference_grid() {
        let s = PlyStrength::cfrp_typical();
        let sigmas = [
            -1.0e6, -10.0e6, -50.0e6, -100.0e6, -200.0e6, -250.0e6, -400.0e6,
        ];
        let taus = [5.0e6, 20.0e6, 45.0e6, 70.0e6, 90.0e6];
        for &sigma2 in &sigmas {
            for &tau in &taus {
                let got = puck_failure([0.0, sigma2, tau], &s).iff_effort;
                let want = puck_reference_compression(sigma2, tau, s.yc, s.s12, 0.25);
                assert!(
                    (got - want).abs() <= 1e-12 * want.abs().max(1.0),
                    "sigma2 {sigma2}, tau {tau}: got {got}, want {want}"
                );
            }
        }
    }

    #[test]
    fn test_puck_iff_modes_b_and_c_are_continuous_at_the_boundary() {
        let s = PlyStrength::cfrp_typical();
        let p = 0.25;
        let r_a = s.s12 / (2.0 * p) * ((1.0 + 2.0 * p * s.yc / s.s12).sqrt() - 1.0);
        let p_pp = p * r_a / s.s12;
        let tau_c = s.s12 * (1.0 + 2.0 * p_pp).sqrt();
        let tau = 40.0e6;
        let sigma_boundary = -(r_a / tau_c) * tau;

        // The two textbook branches agree exactly at the boundary.
        let mode_b =
            ((tau * tau + (p * sigma_boundary).powi(2)).sqrt() + p * sigma_boundary) / s.s12;
        let mode_c =
            ((tau / (2.0 * (1.0 + p_pp) * s.s12)).powi(2) + (sigma_boundary / s.yc).powi(2)) * s.yc
                / (-sigma_boundary);
        assert!(
            (mode_b - mode_c).abs() < 1e-12,
            "branch mismatch: B {mode_b}, C {mode_c}"
        );

        // The implementation is continuous across it.
        let below = puck_failure([0.0, sigma_boundary * (1.0 - 1.0e-9), tau], &s).iff_effort;
        let above = puck_failure([0.0, sigma_boundary * (1.0 + 1.0e-9), tau], &s).iff_effort;
        assert!(
            (below - above).abs() < 1e-7,
            "discontinuity at the B/C boundary: {below} vs {above}"
        );
        assert!((below - mode_b).abs() < 1e-7 && (above - mode_b).abs() < 1e-7);
    }

    #[test]
    fn test_puck_iff_transverse_compression_raises_shear_strength() {
        // Mode B: friction on the action plane lets the ply carry the full
        // in-plane shear strength plus a bit more when it is also compressed.
        let s = PlyStrength::cfrp_typical();
        let pure_shear = puck_failure([0.0, 0.0, s.s12], &s).iff_effort;
        let compressed = puck_failure([0.0, -20.0e6, s.s12], &s).iff_effort;
        assert!((pure_shear - 1.0).abs() < 1e-12);
        assert!(
            compressed < 1.0 && compressed > 0.8,
            "compressed shear effort {compressed}"
        );
        // Vanishing compression recovers the pure-shear limit.
        let limit = puck_failure([0.0, -1.0, s.s12], &s).iff_effort;
        assert!((limit - 1.0).abs() < 1e-6, "limit {limit}");
    }

    #[test]
    fn test_puck_iff_zero_friction_limit_is_finite() {
        // p = 0 must not produce 0/0 in R_perp_perp_A.
        let yc = 250.0e6;
        let s21 = 70.0e6;
        let shear = puck_iff_transverse_compression(-10.0e6, 35.0e6, yc, s21, 0.0);
        assert!(
            (shear - 0.5).abs() < 1e-12,
            "p = 0 Mode B reduces to tau/S21 = 0.5, got {shear}"
        );
        let compression = puck_iff_transverse_compression(-125.0e6, 0.0, yc, s21, 0.0);
        assert!((compression - 0.5).abs() < 1e-12, "got {compression}");
    }

    // --- Ply-by-ply stress ---

    #[test]
    fn test_ply_by_ply_stress_single_ply() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);

        let strain = [0.001, 0.0, 0.0];
        let curv = [0.0, 0.0, 0.0];
        let results = ply_by_ply_stress(&lam, strain, curv);

        assert_eq!(results.len(), 1);
        assert!(
            results[0].stress[0] > 0.0,
            "tensile stress in fiber direction"
        );
    }

    #[test]
    fn test_ply_by_ply_stress_symmetric_laminate() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);
        lam.add_ply(ply, 0.0);

        let strain = [0.001, 0.0, 0.0];
        let curv = [0.0, 0.0, 0.0];
        let results = ply_by_ply_stress(&lam, strain, curv);

        assert_eq!(results.len(), 3);
        // 0-degree plies should have fiber-direction tension
        assert!(results[0].stress[0] > 0.0);
        assert!(results[2].stress[0] > 0.0);
    }

    #[test]
    fn test_ply_by_ply_stress_with_curvature() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 0.0);

        let strain = [0.0, 0.0, 0.0];
        let curv = [10.0, 0.0, 0.0]; // bending
        let results = ply_by_ply_stress(&lam, strain, curv);

        // Top and bottom plies should have opposite signs of stress
        let s_top = results[1].stress[0];
        let s_bot = results[0].stress[0];
        assert!(
            s_top * s_bot < 0.0 || (s_top.abs() < 1e-3 && s_bot.abs() < 1e-3),
            "bending should cause opposite signs: top={s_top}, bot={s_bot}"
        );
    }

    // --- Progressive failure ---

    #[test]
    fn test_progressive_failure_returns_result() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);
        lam.add_ply(ply, 0.0);

        let s = PlyStrength::cfrp_typical();
        let base_strain = [0.01, 0.0, 0.0]; // 1% strain at full load

        let result = progressive_failure_analysis(&lam, base_strain, &s, 2.0, 100);
        assert!(result.first_ply_failure_load > 0.0);
        assert!(result.first_ply_failure_load <= result.ultimate_load);
        assert!(result.n_failed_plies > 0);
    }

    #[test]
    fn test_progressive_failure_90deg_fails_first() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0); // transverse ply, weaker in x-direction
        lam.add_ply(ply, 0.0);

        let s = PlyStrength::cfrp_typical();
        let base_strain = [0.01, 0.0, 0.0];

        let result = progressive_failure_analysis(&lam, base_strain, &s, 2.0, 200);
        // The 90-degree ply should fail first (its transverse strength is lower)
        assert_eq!(
            result.first_failed_ply, 1,
            "90-degree ply should fail first"
        );
    }

    // --- Thermal residual stresses ---

    #[test]
    fn test_thermal_strain() {
        let alpha = [1.0e-6, 30.0e-6]; // typical CFRP CTEs
        let dt = -150.0; // cooling from cure
        let eps = thermal_strain(alpha, dt);
        assert!((eps[0] - (-150.0e-6)).abs() < 1e-12);
        assert!((eps[1] - (-4500.0e-6)).abs() < 1e-12);
        assert!(eps[2].abs() < 1e-15);
    }

    #[test]
    fn test_thermal_residual_single_ply_zero() {
        // A single ply has no CTE mismatch with itself
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);

        let alpha = [1.0e-6, 30.0e-6];
        let residuals = thermal_residual_stresses(&lam, alpha, -150.0);
        assert_eq!(residuals.len(), 1);
        // Single unconstrained ply should have near-zero residual stress
        for s in &residuals[0] {
            assert!(s.abs() < 1.0, "single ply residual should be ~0, got {s}");
        }
    }

    #[test]
    fn test_thermal_residual_crossply_has_stress() {
        // [0/90] crossply should develop residual stresses from CTE mismatch
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);

        let alpha = [1.0e-6, 30.0e-6];
        let residuals = thermal_residual_stresses(&lam, alpha, -150.0);
        assert_eq!(residuals.len(), 2);
        // Should have non-zero residual stresses
        let has_nonzero = residuals
            .iter()
            .any(|s| s[0].abs() > 1.0 || s[1].abs() > 1.0);
        assert!(has_nonzero, "crossply should have residual stresses");
    }

    // --- Ply strength constructors ---

    #[test]
    fn test_ply_strength_cfrp() {
        let s = PlyStrength::cfrp_typical();
        assert!(s.xt > s.yt, "fiber strength should exceed transverse");
        assert!(
            s.xc > s.yc,
            "fiber compression should exceed transverse compression"
        );
    }

    #[test]
    fn test_ply_strength_gfrp() {
        let s = PlyStrength::gfrp_typical();
        assert!(s.xt > 0.0);
        assert!(s.s12 > 0.0);
    }

    #[test]
    fn test_ply_strength_new() {
        let s = PlyStrength::new(100.0, 80.0, 10.0, 50.0, 20.0);
        assert_eq!(s.xt, 100.0);
        assert_eq!(s.yc, 50.0);
    }

    // --- Transform strain ---

    #[test]
    fn test_transform_strain_zero_angle() {
        let strain = [0.001, 0.002, 0.0005];
        let transformed = transform_strain_to_material(strain, 0.0);
        for i in 0..3 {
            assert!((transformed[i] - strain[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_transform_strain_90deg() {
        let strain = [0.001, 0.002, 0.0];
        let transformed = transform_strain_to_material(strain, 90.0);
        // At 90 degrees, eps_1 and eps_2 should swap
        assert!((transformed[0] - 0.002).abs() < 1e-10);
        assert!((transformed[1] - 0.001).abs() < 1e-10);
    }
}
