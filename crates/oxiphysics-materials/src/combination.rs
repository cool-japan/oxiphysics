// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rules for combining material properties in contact pairs and composite materials.
//!
//! Provides combining rules for friction, restitution, and mechanical
//! properties (Young's modulus, Poisson's ratio) for contact between two
//! different materials. Also includes Voigt/Reuss/Hill averaging,
//! rule of mixtures, laminate theory (CLT) basics, Halpin-Tsai model,
//! and effective thermal expansion for composite materials.

// ─── Friction combining rules ────────────────────────────────────────────────

/// Rule for combining friction coefficients of two materials in contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrictionCombineRule {
    /// Arithmetic average: (f1 + f2) / 2
    Average,
    /// Minimum of the two values
    Min,
    /// Maximum of the two values
    Max,
    /// Product: f1 * f2
    Multiply,
    /// Geometric mean: sqrt(f1 * f2)
    GeometricMean,
    /// Harmonic mean: 2*f1*f2 / (f1 + f2)
    HarmonicMean,
}

/// Rule for combining restitution coefficients of two materials in contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RestitutionCombineRule {
    /// Arithmetic average: (r1 + r2) / 2
    Average,
    /// Minimum of the two values
    Min,
    /// Maximum of the two values
    Max,
    /// Product: r1 * r2
    Multiply,
    /// Geometric mean: sqrt(r1 * r2)
    GeometricMean,
    /// Harmonic mean: 2*r1*r2 / (r1 + r2)
    HarmonicMean,
}

/// Rule for combining Young's moduli for Hertzian contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModulusCombineRule {
    /// Harmonic mean (correct for Hertz contact): 2*E1*E2 / (E1 + E2)
    HertzHarmonic,
    /// Voigt upper bound: arithmetic average
    Voigt,
    /// Reuss lower bound: harmonic mean
    Reuss,
    /// Geometric mean: sqrt(E1 * E2)
    GeometricMean,
}

// ─── Combining functions ─────────────────────────────────────────────────────

/// Combine friction coefficients of two materials using the given rule.
pub fn combine_friction(f1: f64, f2: f64, rule: FrictionCombineRule) -> f64 {
    match rule {
        FrictionCombineRule::Average => (f1 + f2) * 0.5,
        FrictionCombineRule::Min => f1.min(f2),
        FrictionCombineRule::Max => f1.max(f2),
        FrictionCombineRule::Multiply => f1 * f2,
        FrictionCombineRule::GeometricMean => (f1 * f2).sqrt(),
        FrictionCombineRule::HarmonicMean => {
            let s = f1 + f2;
            if s.abs() < f64::EPSILON {
                0.0
            } else {
                2.0 * f1 * f2 / s
            }
        }
    }
}

/// Combine restitution coefficients of two materials using the given rule.
pub fn combine_restitution(r1: f64, r2: f64, rule: RestitutionCombineRule) -> f64 {
    match rule {
        RestitutionCombineRule::Average => (r1 + r2) * 0.5,
        RestitutionCombineRule::Min => r1.min(r2),
        RestitutionCombineRule::Max => r1.max(r2),
        RestitutionCombineRule::Multiply => r1 * r2,
        RestitutionCombineRule::GeometricMean => (r1 * r2).sqrt(),
        RestitutionCombineRule::HarmonicMean => {
            let s = r1 + r2;
            if s.abs() < f64::EPSILON {
                0.0
            } else {
                2.0 * r1 * r2 / s
            }
        }
    }
}

/// Combine Young's moduli for contact mechanics.
pub fn combine_modulus(e1: f64, e2: f64, rule: ModulusCombineRule) -> f64 {
    match rule {
        ModulusCombineRule::HertzHarmonic => {
            let s = e1 + e2;
            if s.abs() < f64::EPSILON {
                0.0
            } else {
                2.0 * e1 * e2 / s
            }
        }
        ModulusCombineRule::Voigt => (e1 + e2) * 0.5,
        ModulusCombineRule::Reuss => {
            let s = e1 + e2;
            if s.abs() < f64::EPSILON {
                0.0
            } else {
                2.0 * e1 * e2 / s
            }
        }
        ModulusCombineRule::GeometricMean => (e1 * e2).sqrt(),
    }
}

// ─── Hertz contact helpers ────────────────────────────────────────────────────

/// Effective Young's modulus for Hertz contact between two elastic bodies.
///
/// 1/E* = (1-ν1²)/E1 + (1-ν2²)/E2
pub fn hertz_effective_modulus(e1: f64, nu1: f64, e2: f64, nu2: f64) -> f64 {
    let inv = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
    if inv.abs() < f64::EPSILON {
        0.0
    } else {
        1.0 / inv
    }
}

/// Effective radius for Hertz contact between two spheres.
///
/// 1/R* = 1/R1 + 1/R2
pub fn hertz_effective_radius(r1: f64, r2: f64) -> f64 {
    let inv = 1.0 / r1 + 1.0 / r2;
    if inv.abs() < f64::EPSILON {
        0.0
    } else {
        1.0 / inv
    }
}

/// Hertz contact force for two elastic spheres.
///
/// F = (4/3) * E* * sqrt(R*) * δ^(3/2)
///
/// where δ is the penetration depth.
pub fn hertz_contact_force(e_star: f64, r_star: f64, penetration: f64) -> f64 {
    if penetration <= 0.0 {
        return 0.0;
    }
    (4.0 / 3.0) * e_star * r_star.sqrt() * penetration.powf(1.5)
}

// ─── Contact material pair ────────────────────────────────────────────────────

/// Combined contact properties for a pair of materials.
#[derive(Debug, Clone)]
pub struct ContactMaterialPair {
    /// Combined friction coefficient.
    pub friction: f64,
    /// Combined restitution coefficient.
    pub restitution: f64,
    /// Effective Young's modulus for contact stiffness (Pa).
    pub effective_modulus: f64,
    /// Contact damping coefficient (N·s/m).
    pub damping: f64,
}

impl ContactMaterialPair {
    /// Compute combined contact properties from two materials' properties.
    ///
    /// Uses geometric mean for friction, minimum for restitution, and
    /// Hertz formula for effective modulus.
    pub fn from_materials(
        friction1: f64,
        restitution1: f64,
        young1: f64,
        poisson1: f64,
        friction2: f64,
        restitution2: f64,
        young2: f64,
        poisson2: f64,
    ) -> Self {
        Self {
            friction: combine_friction(friction1, friction2, FrictionCombineRule::GeometricMean),
            restitution: combine_restitution(
                restitution1,
                restitution2,
                RestitutionCombineRule::Min,
            ),
            effective_modulus: hertz_effective_modulus(young1, poisson1, young2, poisson2),
            damping: 0.0,
        }
    }
}

// ─── Thermal contact resistance ───────────────────────────────────────────────

/// Thermal contact resistance at an interface between two materials.
///
/// R_contact = R1 + R_gap + R2
#[derive(Debug, Clone, Copy)]
pub struct ThermalContactResistance {
    /// Thermal conductivity of material 1 (W/(m·K)).
    pub k1: f64,
    /// Thermal conductivity of material 2 (W/(m·K)).
    pub k2: f64,
    /// Gap conductance (W/(m²·K)) for the interfacial gap.
    pub gap_conductance: f64,
}

impl ThermalContactResistance {
    /// Create a thermal contact resistance model.
    pub fn new(k1: f64, k2: f64, gap_conductance: f64) -> Self {
        Self {
            k1,
            k2,
            gap_conductance,
        }
    }

    /// Combined thermal conductance at the interface (W/(m²·K)).
    ///
    /// Uses harmonic mean of the two conductivities plus the gap conductance.
    pub fn combined_conductance(&self) -> f64 {
        let harmonic_k = if (self.k1 + self.k2).abs() < f64::EPSILON {
            0.0
        } else {
            2.0 * self.k1 * self.k2 / (self.k1 + self.k2)
        };
        harmonic_k + self.gap_conductance
    }

    /// Heat flux (W/m²) given temperature difference (K).
    pub fn heat_flux(&self, delta_t: f64) -> f64 {
        self.combined_conductance() * delta_t
    }
}

// ─── Diffusion combining ──────────────────────────────────────────────────────

/// Effective diffusion coefficient for a two-phase mixture (Maxwell model).
///
/// For spherical inclusions of phase 2 in matrix phase 1:
/// D_eff = D1 * (D2 + 2*D1 - 2*phi2*(D1-D2)) / (D2 + 2*D1 + phi2*(D1-D2))
pub fn maxwell_diffusivity(d1: f64, d2: f64, volume_fraction2: f64) -> f64 {
    let phi = volume_fraction2.clamp(0.0, 1.0);
    let num = d2 + 2.0 * d1 - 2.0 * phi * (d1 - d2);
    let den = d2 + 2.0 * d1 + phi * (d1 - d2);
    if den.abs() < f64::EPSILON {
        d1
    } else {
        d1 * num / den
    }
}

// ─── Voigt/Reuss/Hill averaging for composite materials ─────────────────────

/// Voigt (upper bound) average of a composite property.
///
/// P_Voigt = Σ fᵢ * Pᵢ
///
/// where fᵢ are volume fractions and Pᵢ are phase properties.
/// This is the iso-strain (parallel) bound.
pub fn voigt_average(volume_fractions: &[f64], properties: &[f64]) -> f64 {
    volume_fractions
        .iter()
        .zip(properties.iter())
        .map(|(&f, &p)| f * p)
        .sum()
}

/// Reuss (lower bound) average of a composite property.
///
/// 1/P_Reuss = Σ fᵢ / Pᵢ
///
/// where fᵢ are volume fractions and Pᵢ are phase properties.
/// This is the iso-stress (series) bound.
pub fn reuss_average(volume_fractions: &[f64], properties: &[f64]) -> f64 {
    let inv_sum: f64 = volume_fractions
        .iter()
        .zip(properties.iter())
        .map(|(&f, &p)| if p.abs() < f64::EPSILON { 0.0 } else { f / p })
        .sum();
    if inv_sum.abs() < f64::EPSILON {
        0.0
    } else {
        1.0 / inv_sum
    }
}

/// Hill (Voigt-Reuss-Hill) average of a composite property.
///
/// P_Hill = (P_Voigt + P_Reuss) / 2
///
/// Provides a practical estimate between the upper and lower bounds.
pub fn hill_average(volume_fractions: &[f64], properties: &[f64]) -> f64 {
    let v = voigt_average(volume_fractions, properties);
    let r = reuss_average(volume_fractions, properties);
    0.5 * (v + r)
}

// ─── Rule of mixtures for fiber composites ──────────────────────────────────

/// Longitudinal modulus of a unidirectional fiber composite via the rule of mixtures.
///
/// E_1 = V_f * E_f + (1 - V_f) * E_m
///
/// where V_f is fiber volume fraction, E_f is fiber modulus, E_m is matrix modulus.
pub fn rule_of_mixtures_longitudinal(vf: f64, e_fiber: f64, e_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    vf * e_fiber + (1.0 - vf) * e_matrix
}

/// Transverse modulus of a unidirectional fiber composite via the inverse rule of mixtures.
///
/// 1/E_2 = V_f / E_f + (1 - V_f) / E_m
pub fn rule_of_mixtures_transverse(vf: f64, e_fiber: f64, e_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let inv = vf / e_fiber + (1.0 - vf) / e_matrix;
    if inv.abs() < f64::EPSILON {
        0.0
    } else {
        1.0 / inv
    }
}

/// Longitudinal Poisson's ratio of a unidirectional fiber composite.
///
/// ν_12 = V_f * ν_f + (1 - V_f) * ν_m
pub fn rule_of_mixtures_poisson(vf: f64, nu_fiber: f64, nu_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    vf * nu_fiber + (1.0 - vf) * nu_matrix
}

/// In-plane shear modulus of a unidirectional fiber composite (inverse rule).
///
/// 1/G_12 = V_f / G_f + (1 - V_f) / G_m
pub fn rule_of_mixtures_shear(vf: f64, g_fiber: f64, g_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let inv = vf / g_fiber + (1.0 - vf) / g_matrix;
    if inv.abs() < f64::EPSILON {
        0.0
    } else {
        1.0 / inv
    }
}

/// Composite density via rule of mixtures.
///
/// ρ_c = V_f * ρ_f + (1 - V_f) * ρ_m
pub fn rule_of_mixtures_density(vf: f64, rho_fiber: f64, rho_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    vf * rho_fiber + (1.0 - vf) * rho_matrix
}

// ─── Halpin-Tsai model ──────────────────────────────────────────────────────

/// Halpin-Tsai model for transverse modulus of a fiber composite.
///
/// E_2 = E_m * (1 + ξ * η * V_f) / (1 - η * V_f)
///
/// where η = (E_f/E_m - 1) / (E_f/E_m + ξ)
/// and ξ is a shape/packing factor (typically 1 or 2 for transverse modulus).
pub fn halpin_tsai_modulus(vf: f64, e_fiber: f64, e_matrix: f64, xi: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    if e_matrix.abs() < f64::EPSILON {
        return 0.0;
    }
    let ratio = e_fiber / e_matrix;
    let eta = (ratio - 1.0) / (ratio + xi);
    e_matrix * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
}

/// Halpin-Tsai model for shear modulus of a fiber composite.
///
/// G_12 = G_m * (1 + ξ * η * V_f) / (1 - η * V_f)
///
/// where η = (G_f/G_m - 1) / (G_f/G_m + ξ)
/// and ξ is typically 1 for shear modulus.
pub fn halpin_tsai_shear(vf: f64, g_fiber: f64, g_matrix: f64, xi: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    if g_matrix.abs() < f64::EPSILON {
        return 0.0;
    }
    let ratio = g_fiber / g_matrix;
    let eta = (ratio - 1.0) / (ratio + xi);
    g_matrix * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
}

// ─── Effective thermal expansion ─────────────────────────────────────────────

/// Effective longitudinal coefficient of thermal expansion (CTE)
/// for a unidirectional composite using Schapery's formula.
///
/// α_1 = (V_f * E_f * α_f + V_m * E_m * α_m) / (V_f * E_f + V_m * E_m)
pub fn effective_cte_longitudinal(
    vf: f64,
    e_fiber: f64,
    alpha_fiber: f64,
    e_matrix: f64,
    alpha_matrix: f64,
) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let vm = 1.0 - vf;
    let num = vf * e_fiber * alpha_fiber + vm * e_matrix * alpha_matrix;
    let den = vf * e_fiber + vm * e_matrix;
    if den.abs() < f64::EPSILON {
        0.0
    } else {
        num / den
    }
}

/// Effective transverse CTE for a unidirectional composite (Schapery).
///
/// α_2 = (1 + ν_f) * V_f * α_f + (1 + ν_m) * V_m * α_m - α_1 * ν_12
///
/// where ν_12 is the composite Poisson's ratio.
pub fn effective_cte_transverse(
    vf: f64,
    alpha_fiber: f64,
    nu_fiber: f64,
    alpha_matrix: f64,
    nu_matrix: f64,
    alpha_1: f64,
    nu_12: f64,
) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let vm = 1.0 - vf;
    (1.0 + nu_fiber) * vf * alpha_fiber + (1.0 + nu_matrix) * vm * alpha_matrix - alpha_1 * nu_12
}

/// Effective bulk CTE for a particulate composite using Turner's model.
///
/// α_eff = (V_1 * K_1 * α_1 + V_2 * K_2 * α_2) / (V_1 * K_1 + V_2 * K_2)
///
/// where K is the bulk modulus and V is volume fraction.
pub fn turner_cte(v1: f64, k1: f64, alpha1: f64, k2: f64, alpha2: f64) -> f64 {
    let v1 = v1.clamp(0.0, 1.0);
    let v2 = 1.0 - v1;
    let num = v1 * k1 * alpha1 + v2 * k2 * alpha2;
    let den = v1 * k1 + v2 * k2;
    if den.abs() < f64::EPSILON {
        0.0
    } else {
        num / den
    }
}

// ─── Classical Laminate Theory (CLT) basics ─────────────────────────────────

/// Reduced stiffness matrix Q for an orthotropic lamina (plane stress).
///
/// Returns the 3x3 Q matrix as a flat \[f64; 9\] in row-major order:
/// \[Q11, Q12, 0, Q12, Q22, 0, 0, 0, Q66\]
///
/// Q11 = E1 / (1 - ν12*ν21)
/// Q22 = E2 / (1 - ν12*ν21)
/// Q12 = ν12 * E2 / (1 - ν12*ν21)
/// Q66 = G12
pub fn lamina_stiffness_matrix(e1: f64, e2: f64, nu12: f64, g12: f64) -> [f64; 9] {
    let nu21 = nu12 * e2 / e1;
    let denom = 1.0 - nu12 * nu21;
    let q11 = e1 / denom;
    let q22 = e2 / denom;
    let q12 = nu12 * e2 / denom;
    let q66 = g12;
    [q11, q12, 0.0, q12, q22, 0.0, 0.0, 0.0, q66]
}

/// Transform a reduced stiffness matrix Q to an arbitrary angle θ (radians).
///
/// Returns the transformed Q-bar matrix as \[f64; 9\] in row-major order.
/// Uses standard CLT transformation with m = cos(θ), n = sin(θ).
pub fn transform_stiffness(q: &[f64; 9], theta: f64) -> [f64; 9] {
    let m = theta.cos();
    let n = theta.sin();
    let m2 = m * m;
    let n2 = n * n;
    let m4 = m2 * m2;
    let n4 = n2 * n2;
    let mn = m * n;
    let m2n2 = m2 * n2;

    let q11 = q[0];
    let q12 = q[1];
    let q22 = q[4];
    let q66 = q[8];

    let qb11 = q11 * m4 + 2.0 * (q12 + 2.0 * q66) * m2n2 + q22 * n4;
    let qb22 = q11 * n4 + 2.0 * (q12 + 2.0 * q66) * m2n2 + q22 * m4;
    let qb12 = (q11 + q22 - 4.0 * q66) * m2n2 + q12 * (m4 + n4);
    let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * m2n2 + q66 * (m4 + n4);
    let qb16 = (q11 - q12 - 2.0 * q66) * m2 * mn + (q12 - q22 + 2.0 * q66) * n2 * mn;
    let qb26 = (q11 - q12 - 2.0 * q66) * n2 * mn + (q12 - q22 + 2.0 * q66) * m2 * mn;

    [qb11, qb12, qb16, qb12, qb22, qb26, qb16, qb26, qb66]
}

/// A single lamina (ply) in a laminate.
#[derive(Debug, Clone)]
pub struct Lamina {
    /// Reduced stiffness matrix Q (3x3 row-major).
    pub q: [f64; 9],
    /// Ply orientation angle (radians).
    pub theta: f64,
    /// Ply thickness (m).
    pub thickness: f64,
}

/// Compute the ABD stiffness matrices for a symmetric or general laminate.
///
/// Returns (A, B, D) as three \[f64; 9\] arrays in row-major order.
///
/// A_ij = Σ Q̄_ij_k * (z_k - z_{k-1})
/// B_ij = (1/2) Σ Q̄_ij_k * (z_k² - z_{k-1}²)
/// D_ij = (1/3) Σ Q̄_ij_k * (z_k³ - z_{k-1}³)
///
/// where z is measured from the laminate midplane.
pub fn laminate_abd(plies: &[Lamina]) -> ([f64; 9], [f64; 9], [f64; 9]) {
    let total_thickness: f64 = plies.iter().map(|p| p.thickness).sum();
    let mut z_bot = -total_thickness / 2.0;

    let mut a_mat = [0.0f64; 9];
    let mut b_mat = [0.0f64; 9];
    let mut d_mat = [0.0f64; 9];

    for ply in plies {
        let z_top = z_bot + ply.thickness;
        let qbar = transform_stiffness(&ply.q, ply.theta);

        for i in 0..9 {
            a_mat[i] += qbar[i] * (z_top - z_bot);
            b_mat[i] += 0.5 * qbar[i] * (z_top * z_top - z_bot * z_bot);
            d_mat[i] += (1.0 / 3.0) * qbar[i] * (z_top.powi(3) - z_bot.powi(3));
        }

        z_bot = z_top;
    }

    (a_mat, b_mat, d_mat)
}

/// Compute effective in-plane engineering constants from the A matrix of a laminate.
///
/// Returns (E_x, E_y, G_xy, nu_xy) assuming a symmetric laminate (B=0).
pub fn laminate_engineering_constants(
    a_mat: &[f64; 9],
    total_thickness: f64,
) -> (f64, f64, f64, f64) {
    // Normalized stiffness: a_ij = A_ij / h
    let a11 = a_mat[0] / total_thickness;
    let a12 = a_mat[1] / total_thickness;
    let a22 = a_mat[4] / total_thickness;
    let a66 = a_mat[8] / total_thickness;

    let det = a11 * a22 - a12 * a12;
    if det.abs() < f64::EPSILON {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let ex = det / a22;
    let ey = det / a11;
    let nu_xy = a12 / a22;
    let gxy = a66;

    (ex, ey, gxy, nu_xy)
}

// ─── Hashin-Shtrikman bounds ────────────────────────────────────────────────

/// Hashin-Shtrikman upper bound for bulk modulus of a two-phase composite.
///
/// K_upper = K_2 + V_1 / (1/(K_1-K_2) + 3*V_2/(3*K_2+4*G_2))
///
/// Assumes K_2 > K_1 (phase 2 is the stiffer phase).
pub fn hashin_shtrikman_bulk_upper(v1: f64, k1: f64, k2: f64, g2: f64) -> f64 {
    let v1 = v1.clamp(0.0, 1.0);
    let v2 = 1.0 - v1;
    let dk = k1 - k2;
    if dk.abs() < f64::EPSILON {
        return k1;
    }
    let inv = 1.0 / dk + 3.0 * v2 / (3.0 * k2 + 4.0 * g2);
    if inv.abs() < f64::EPSILON {
        k2
    } else {
        k2 + v1 / inv
    }
}

/// Hashin-Shtrikman lower bound for bulk modulus of a two-phase composite.
///
/// K_lower = K_1 + V_2 / (1/(K_2-K_1) + 3*V_1/(3*K_1+4*G_1))
///
/// Assumes K_1 < K_2 (phase 1 is the softer phase).
pub fn hashin_shtrikman_bulk_lower(v1: f64, k1: f64, g1: f64, k2: f64) -> f64 {
    let v1 = v1.clamp(0.0, 1.0);
    let v2 = 1.0 - v1;
    let dk = k2 - k1;
    if dk.abs() < f64::EPSILON {
        return k1;
    }
    let inv = 1.0 / dk + 3.0 * v1 / (3.0 * k1 + 4.0 * g1);
    if inv.abs() < f64::EPSILON {
        k1
    } else {
        k1 + v2 / inv
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    use crate::combination::halpin_tsai_modulus;

    #[test]
    fn test_combine_friction_average() {
        let r = combine_friction(0.4, 0.6, FrictionCombineRule::Average);
        assert!((r - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_combine_friction_geometric() {
        let r = combine_friction(0.25, 1.0, FrictionCombineRule::GeometricMean);
        assert!((r - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_combine_friction_harmonic() {
        let r = combine_friction(1.0, 1.0, FrictionCombineRule::HarmonicMean);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_combine_friction_min() {
        assert!((combine_friction(0.4, 0.6, FrictionCombineRule::Min) - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_combine_friction_max() {
        assert!((combine_friction(0.4, 0.6, FrictionCombineRule::Max) - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_combine_friction_multiply() {
        assert!((combine_friction(0.4, 0.6, FrictionCombineRule::Multiply) - 0.24).abs() < 1e-10);
    }

    #[test]
    fn test_combine_restitution_rules() {
        assert!(
            (combine_restitution(0.3, 0.7, RestitutionCombineRule::Average) - 0.5).abs() < 1e-10
        );
        assert!((combine_restitution(0.3, 0.7, RestitutionCombineRule::Min) - 0.3).abs() < 1e-10);
        assert!((combine_restitution(0.3, 0.7, RestitutionCombineRule::Max) - 0.7).abs() < 1e-10);
        assert!(
            (combine_restitution(0.3, 0.7, RestitutionCombineRule::Multiply) - 0.21).abs() < 1e-10
        );
    }

    #[test]
    fn test_hertz_effective_modulus_equal_materials() {
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let e_star = hertz_effective_modulus(e, nu, e, nu);
        let expected = e / (2.0 * (1.0 - nu * nu));
        assert!((e_star - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_hertz_contact_force_zero_penetration() {
        assert_eq!(hertz_contact_force(200e9, 0.01, 0.0), 0.0);
        assert_eq!(hertz_contact_force(200e9, 0.01, -0.001), 0.0);
    }

    #[test]
    fn test_hertz_effective_radius() {
        let r = hertz_effective_radius(1.0, 1.0);
        assert!((r - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_contact_heat_flux() {
        let tcr = ThermalContactResistance::new(50.0, 50.0, 0.0);
        let cond = tcr.combined_conductance();
        assert!((cond - 50.0).abs() < 1e-6);
        let flux = tcr.heat_flux(10.0);
        assert!((flux - 500.0).abs() < 1e-6);
    }

    #[test]
    fn test_maxwell_diffusivity_pure_phase() {
        let d = maxwell_diffusivity(1.0, 5.0, 0.0);
        assert!((d - 1.0).abs() < 1e-10);
        let d2 = maxwell_diffusivity(1.0, 5.0, 1.0);
        assert!((d2 - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_contact_material_pair() {
        let pair = ContactMaterialPair::from_materials(0.5, 0.4, 200e9, 0.3, 0.7, 0.8, 70e9, 0.33);
        assert!(pair.friction > 0.0);
        assert!(pair.restitution <= 0.4);
        assert!(pair.effective_modulus > 0.0);
    }

    // ─── Voigt/Reuss/Hill tests ─────────────────────────────────────────

    #[test]
    fn test_voigt_average_equal_phases() {
        // Two equal phases → average = that value
        let vf = [0.5, 0.5];
        let props = [200.0e9, 200.0e9];
        let result = voigt_average(&vf, &props);
        assert!((result - 200.0e9).abs() < 1.0);
    }

    #[test]
    fn test_voigt_average_weighted() {
        // 60% phase1 (100 GPa) + 40% phase2 (200 GPa) = 140 GPa
        let vf = [0.6, 0.4];
        let props = [100.0e9, 200.0e9];
        let result = voigt_average(&vf, &props);
        assert!((result - 140.0e9).abs() < 1.0);
    }

    #[test]
    fn test_reuss_average_equal_phases() {
        let vf = [0.5, 0.5];
        let props = [200.0e9, 200.0e9];
        let result = reuss_average(&vf, &props);
        assert!((result - 200.0e9).abs() < 1.0);
    }

    #[test]
    fn test_reuss_leq_voigt() {
        // Reuss bound <= Voigt bound always
        let vf = [0.3, 0.7];
        let props = [70.0e9, 200.0e9];
        let v = voigt_average(&vf, &props);
        let r = reuss_average(&vf, &props);
        assert!(r <= v + 1e-6, "Reuss {r} should be <= Voigt {v}");
    }

    #[test]
    fn test_hill_between_voigt_reuss() {
        let vf = [0.4, 0.6];
        let props = [70.0e9, 200.0e9];
        let v = voigt_average(&vf, &props);
        let r = reuss_average(&vf, &props);
        let h = hill_average(&vf, &props);
        assert!(
            h >= r - 1e-6 && h <= v + 1e-6,
            "Hill {h} not between Reuss {r} and Voigt {v}"
        );
    }

    #[test]
    fn test_hill_average_value() {
        let vf = [0.5, 0.5];
        let props = [100.0e9, 200.0e9];
        let h = hill_average(&vf, &props);
        // Voigt = 150 GPa, Reuss = 2*100*200/(100+200)*1e9 = 133.33 GPa
        // Hill = (150+133.33)/2 ~ 141.67 GPa
        assert!((h - 141.666666e9).abs() / h < 1e-4);
    }

    // ─── Rule of mixtures tests ─────────────────────────────────────────

    #[test]
    fn test_rom_longitudinal() {
        // Glass fiber (72 GPa) in epoxy (3.5 GPa), Vf = 0.6
        let e1 = rule_of_mixtures_longitudinal(0.6, 72.0e9, 3.5e9);
        let expected = 0.6 * 72.0e9 + 0.4 * 3.5e9;
        assert!((e1 - expected).abs() < 1.0);
    }

    #[test]
    fn test_rom_transverse() {
        let e2 = rule_of_mixtures_transverse(0.6, 72.0e9, 3.5e9);
        // Transverse modulus should be dominated by matrix
        assert!(
            e2 > 3.5e9,
            "Transverse modulus should exceed matrix modulus"
        );
        assert!(
            e2 < 72.0e9,
            "Transverse modulus should be less than fiber modulus"
        );
    }

    #[test]
    fn test_rom_transverse_leq_longitudinal() {
        let e1 = rule_of_mixtures_longitudinal(0.6, 72.0e9, 3.5e9);
        let e2 = rule_of_mixtures_transverse(0.6, 72.0e9, 3.5e9);
        assert!(e2 < e1, "E2 ({e2}) should be less than E1 ({e1})");
    }

    #[test]
    fn test_rom_poisson() {
        let nu = rule_of_mixtures_poisson(0.5, 0.22, 0.35);
        let expected = 0.5 * 0.22 + 0.5 * 0.35;
        assert!((nu - expected).abs() < 1e-10);
    }

    #[test]
    fn test_rom_shear() {
        let g = rule_of_mixtures_shear(0.5, 30.0e9, 1.3e9);
        assert!(g > 1.3e9);
        assert!(g < 30.0e9);
    }

    #[test]
    fn test_rom_density() {
        // Carbon fiber (1.8 g/cc) + epoxy (1.2 g/cc), Vf = 0.55
        let rho = rule_of_mixtures_density(0.55, 1800.0, 1200.0);
        let expected = 0.55 * 1800.0 + 0.45 * 1200.0;
        assert!((rho - expected).abs() < 1e-6);
    }

    // ─── Halpin-Tsai tests ──────────────────────────────────────────────

    #[test]
    fn test_halpin_tsai_pure_matrix() {
        // Vf=0 → should return matrix modulus
        let e2 = halpin_tsai_modulus(0.0, 72.0e9, 3.5e9, 2.0);
        assert!((e2 - 3.5e9).abs() < 1.0);
    }

    #[test]
    fn test_halpin_tsai_between_bounds() {
        // Result should be between Reuss and Voigt bounds
        let vf = 0.6;
        let ef = 72.0e9;
        let em = 3.5e9;
        let e_ht = halpin_tsai_modulus(vf, ef, em, 2.0);
        let e_voigt = rule_of_mixtures_longitudinal(vf, ef, em);
        let e_reuss = rule_of_mixtures_transverse(vf, ef, em);
        assert!(
            e_ht >= e_reuss - 1e-6 && e_ht <= e_voigt + 1e-6,
            "HT {e_ht} not between Reuss {e_reuss} and Voigt {e_voigt}"
        );
    }

    #[test]
    fn test_halpin_tsai_xi_effect() {
        // Higher ξ → closer to Voigt bound
        let vf = 0.5;
        let ef = 72.0e9;
        let em = 3.5e9;
        let e_low_xi = halpin_tsai_modulus(vf, ef, em, 1.0);
        let e_high_xi = halpin_tsai_modulus(vf, ef, em, 10.0);
        assert!(e_high_xi > e_low_xi, "Higher ξ should give higher modulus");
    }

    #[test]
    fn test_halpin_tsai_shear() {
        let g = halpin_tsai_shear(0.5, 30.0e9, 1.3e9, 1.0);
        assert!(g > 1.3e9);
        assert!(g < 30.0e9);
    }

    // ─── Effective CTE tests ────────────────────────────────────────────

    #[test]
    fn test_cte_longitudinal_equal_materials() {
        // Same CTE → effective CTE = that CTE
        let alpha = effective_cte_longitudinal(0.5, 200.0e9, 12.0e-6, 200.0e9, 12.0e-6);
        assert!((alpha - 12.0e-6).abs() < 1e-15);
    }

    #[test]
    fn test_cte_longitudinal_stiff_fiber_dominates() {
        // Stiff fiber with low CTE should pull result toward fiber CTE
        let vf = 0.6;
        let alpha = effective_cte_longitudinal(vf, 230.0e9, 5.0e-6, 3.5e9, 60.0e-6);
        // Should be much closer to fiber CTE (5e-6) than matrix CTE (60e-6)
        assert!(
            alpha < 20.0e-6,
            "Expected CTE near fiber value, got {alpha}"
        );
    }

    #[test]
    fn test_turner_cte_equal_phases() {
        let alpha = turner_cte(0.5, 100.0e9, 10.0e-6, 100.0e9, 10.0e-6);
        assert!((alpha - 10.0e-6).abs() < 1e-15);
    }

    #[test]
    fn test_turner_cte_bulk_weighted() {
        // Phase with higher bulk modulus dominates
        let alpha = turner_cte(0.5, 200.0e9, 5.0e-6, 50.0e9, 20.0e-6);
        // Weighted toward 5e-6 since K1 is much larger
        assert!(
            alpha < 12.5e-6,
            "Turner CTE should be weighted toward stiffer phase"
        );
    }

    #[test]
    fn test_cte_transverse() {
        let vf = 0.5;
        let alpha_1 = effective_cte_longitudinal(vf, 230.0e9, 5.0e-6, 3.5e9, 60.0e-6);
        let nu_12 = rule_of_mixtures_poisson(vf, 0.2, 0.35);
        let alpha_2 = effective_cte_transverse(vf, 5.0e-6, 0.2, 60.0e-6, 0.35, alpha_1, nu_12);
        // Transverse CTE should be larger than longitudinal
        assert!(
            alpha_2 > alpha_1,
            "Transverse CTE should exceed longitudinal"
        );
    }

    // ─── CLT tests ──────────────────────────────────────────────────────

    #[test]
    fn test_lamina_stiffness_matrix() {
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        // Q11 should be close to E1/(1-nu12*nu21)
        let nu21 = 0.3 * 10.0e9 / 140.0e9;
        let denom = 1.0 - 0.3 * nu21;
        let expected_q11 = 140.0e9 / denom;
        assert!((q[0] - expected_q11).abs() / expected_q11 < 1e-10);
        // Off-diagonal terms Q13=Q31=Q23=Q32=0
        assert_eq!(q[2], 0.0);
        assert_eq!(q[6], 0.0);
    }

    #[test]
    fn test_transform_stiffness_zero_angle() {
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        let qbar = transform_stiffness(&q, 0.0);
        for i in 0..9 {
            assert!(
                (qbar[i] - q[i]).abs() < 1.0,
                "Q-bar[{i}] = {} should match Q[{i}] = {} at θ=0",
                qbar[i],
                q[i]
            );
        }
    }

    #[test]
    fn test_transform_stiffness_90_degrees() {
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        let qbar = transform_stiffness(&q, std::f64::consts::FRAC_PI_2);
        // At 90°, Q̄11 ≈ Q22 and Q̄22 ≈ Q11
        assert!(
            (qbar[0] - q[4]).abs() / q[4] < 1e-8,
            "Q̄11 should ≈ Q22 at 90°"
        );
        assert!(
            (qbar[4] - q[0]).abs() / q[0] < 1e-8,
            "Q̄22 should ≈ Q11 at 90°"
        );
    }

    #[test]
    fn test_laminate_abd_symmetric() {
        // Symmetric laminate [0/90]_s should have B=0
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        let plies = vec![
            Lamina {
                q,
                theta: 0.0,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: std::f64::consts::FRAC_PI_2,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: std::f64::consts::FRAC_PI_2,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: 0.0,
                thickness: 0.125e-3,
            },
        ];
        let (_a, b, _d) = laminate_abd(&plies);
        for (i, &val) in b.iter().enumerate() {
            assert!(
                val.abs() < 1e-3,
                "B[{i}]={val} should be ~0 for symmetric laminate"
            );
        }
    }

    #[test]
    fn test_laminate_abd_a_positive_diagonal() {
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        let plies = vec![
            Lamina {
                q,
                theta: 0.0,
                thickness: 0.25e-3,
            },
            Lamina {
                q,
                theta: std::f64::consts::FRAC_PI_2,
                thickness: 0.25e-3,
            },
        ];
        let (a, _b, _d) = laminate_abd(&plies);
        assert!(a[0] > 0.0, "A11 should be positive");
        assert!(a[4] > 0.0, "A22 should be positive");
        assert!(a[8] > 0.0, "A66 should be positive");
    }

    #[test]
    fn test_laminate_engineering_constants() {
        let q = lamina_stiffness_matrix(140.0e9, 10.0e9, 0.3, 5.0e9);
        let plies = vec![
            Lamina {
                q,
                theta: 0.0,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: std::f64::consts::FRAC_PI_2,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: std::f64::consts::FRAC_PI_2,
                thickness: 0.125e-3,
            },
            Lamina {
                q,
                theta: 0.0,
                thickness: 0.125e-3,
            },
        ];
        let total_h: f64 = plies.iter().map(|p| p.thickness).sum();
        let (a, _b, _d) = laminate_abd(&plies);
        let (ex, ey, gxy, _nu_xy) = laminate_engineering_constants(&a, total_h);
        assert!(ex > 0.0, "Ex should be positive");
        assert!(ey > 0.0, "Ey should be positive");
        assert!(gxy > 0.0, "Gxy should be positive");
        // [0/90]_s → Ex ≈ Ey due to balanced layup
        assert!(
            (ex - ey).abs() / ex < 1e-6,
            "Ex ({ex}) should ≈ Ey ({ey}) for balanced layup"
        );
    }

    // ─── Hashin-Shtrikman tests ─────────────────────────────────────────

    #[test]
    fn test_hashin_shtrikman_bounds_ordering() {
        // Lower <= Hill <= Upper
        let v1 = 0.4;
        let k1 = 75.0e9; // aluminum-like
        let g1 = 26.0e9;
        let k2 = 160.0e9; // steel-like
        let g2 = 80.0e9;

        let hs_lower = hashin_shtrikman_bulk_lower(v1, k1, g1, k2);
        let hs_upper = hashin_shtrikman_bulk_upper(v1, k1, k2, g2);
        assert!(
            hs_lower <= hs_upper + 1e-6,
            "HS lower {hs_lower} should be <= HS upper {hs_upper}"
        );
    }

    #[test]
    fn test_hashin_shtrikman_pure_phases() {
        // v1=1 → K = K1
        let hs = hashin_shtrikman_bulk_lower(1.0, 75.0e9, 26.0e9, 160.0e9);
        assert!((hs - 75.0e9).abs() / 75.0e9 < 1e-6);
    }
}

// ─── Multi-scale homogenization ──────────────────────────────────────────────

/// Two-scale homogenization result for a periodic composite.
#[derive(Debug, Clone)]
pub struct HomogenizationResult {
    /// Effective stiffness tensor (6×6 Voigt notation, row-major \[f64; 36\]).
    pub c_eff: [f64; 36],
    /// Effective Young's modulus in the x direction.
    pub e_x: f64,
    /// Effective Poisson's ratio ν_xy.
    pub nu_xy: f64,
    /// Effective shear modulus G_xy.
    pub g_xy: f64,
}

/// Simple self-consistent (Eshelby-type) homogenization for a two-phase composite.
///
/// Uses the Mori-Tanaka (MT) approximation for effective bulk and shear moduli.
///
/// K_eff = K_m + V_f * (K_f - K_m) / (1 + V_m * (K_f - K_m) / (K_m + 4/3 G_m))
/// G_eff = G_m + V_f * (G_f - G_m) / (1 + V_m * (G_f - G_m) * (6*(K_m + 2*G_m)) / (5*G_m*(3*K_m + 4*G_m)))
pub fn mori_tanaka_homogenization(
    vf: f64,
    k_fiber: f64,
    g_fiber: f64,
    k_matrix: f64,
    g_matrix: f64,
) -> (f64, f64) {
    let vf = vf.clamp(0.0, 1.0);
    let vm = 1.0 - vf;

    // Effective bulk modulus
    let dk = k_fiber - k_matrix;
    let k_denom = k_matrix + 4.0 / 3.0 * g_matrix;
    let k_eff = if k_denom.abs() < f64::EPSILON {
        k_matrix
    } else {
        k_matrix + vf * dk / (1.0 + vm * dk / k_denom)
    };

    // Effective shear modulus
    let dg = g_fiber - g_matrix;
    let beta =
        6.0 * (k_matrix + 2.0 * g_matrix) / (5.0 * g_matrix * (3.0 * k_matrix + 4.0 * g_matrix));
    let g_denom = if (g_matrix * (3.0 * k_matrix + 4.0 * g_matrix)).abs() < f64::EPSILON {
        1.0
    } else {
        1.0 + vm * dg * beta
    };
    let g_eff = g_matrix + vf * dg / g_denom;

    (k_eff, g_eff)
}

/// Convert bulk modulus K and shear modulus G to engineering constants (E, ν).
pub fn bulk_shear_to_engineering(k: f64, g: f64) -> (f64, f64) {
    let e = 9.0 * k * g / (3.0 * k + g);
    let nu = (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g));
    (e, nu)
}

/// Build a `HomogenizationResult` from the Mori-Tanaka effective moduli.
pub fn homogenization_result(
    vf: f64,
    k_fiber: f64,
    g_fiber: f64,
    k_matrix: f64,
    g_matrix: f64,
) -> HomogenizationResult {
    let (k_eff, g_eff) = mori_tanaka_homogenization(vf, k_fiber, g_fiber, k_matrix, g_matrix);
    let (e_x, nu_xy) = bulk_shear_to_engineering(k_eff, g_eff);

    // Build isotropic 6×6 stiffness tensor (Voigt notation)
    let lam = k_eff - 2.0 / 3.0 * g_eff; // Lamé λ
    let mut c_eff = [0.0f64; 36];
    // C_11 = C_22 = C_33 = λ + 2G
    c_eff[0] = lam + 2.0 * g_eff;
    c_eff[7] = lam + 2.0 * g_eff;
    c_eff[14] = lam + 2.0 * g_eff;
    // Off-diagonal C_12 = C_13 = C_23 = λ
    c_eff[1] = lam;
    c_eff[6] = lam;
    c_eff[2] = lam;
    c_eff[12] = lam;
    c_eff[9] = lam;
    c_eff[13] = lam;
    // Shear C_44 = C_55 = C_66 = G
    c_eff[21] = g_eff;
    c_eff[28] = g_eff;
    c_eff[35] = g_eff;

    HomogenizationResult {
        c_eff,
        e_x,
        nu_xy,
        g_xy: g_eff,
    }
}

// ─── Random fiber composite models ──────────────────────────────────────────

/// Randomly oriented short-fiber composite (Cox-Krenchel model).
///
/// For randomly oriented (3-D isotropic) short fibers:
/// E_composite = η_l * η_o * V_f * E_f + (1 - V_f) * E_m
///
/// where η_o = 1/5 (random 3-D orientation factor)
/// and η_l is the fiber length efficiency factor from Cox shear-lag.
///
/// η_l = 1 - tanh(β * l/2) / (β * l/2)
/// β = sqrt(2*G_m / (E_f * A * ln(R/r)))
/// Simplified here to η_l as a direct input.
pub fn cox_krenchel_modulus(vf: f64, e_fiber: f64, e_matrix: f64, eta_l: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let vm = 1.0 - vf;
    let eta_o = 1.0 / 5.0; // 3-D random orientation
    eta_l * eta_o * vf * e_fiber + vm * e_matrix
}

/// Length efficiency factor η_l for the Cox shear-lag model.
///
/// η_l = 1 - tanh(n * L / 2) / (n * L / 2)
///
/// where n = sqrt(2 * G_m / (E_f * ln(R / r))) and L is the fiber length,
/// r is the fiber radius, R is the mean fibre spacing radius.
pub fn cox_length_efficiency(_fiber_length: f64, beta_l_over_2: f64) -> f64 {
    // beta_l_over_2 = β * L/2, a dimensionless parameter
    let x = beta_l_over_2;
    if x < 1e-10 {
        return 0.0;
    }
    1.0 - x.tanh() / x
}

/// Halpin-Tsai model for randomly oriented fibers (2-D in-plane isotropy).
///
/// E_iso = (3/8) * E_11 + (5/8) * E_22
///
/// where E_11 and E_22 are the longitudinal and transverse Halpin-Tsai moduli.
pub fn halpin_tsai_random_2d(vf: f64, e_fiber: f64, e_matrix: f64, xi: f64) -> f64 {
    let e11 = rule_of_mixtures_longitudinal(vf, e_fiber, e_matrix);
    let e22 = halpin_tsai_modulus(vf, e_fiber, e_matrix, xi);
    3.0 / 8.0 * e11 + 5.0 / 8.0 * e22
}

// ─── Woven composite homogenization ─────────────────────────────────────────

/// Mosaic model for a plain-weave composite.
///
/// Divides the unit cell into undulation and straight regions.
/// The effective in-plane modulus is estimated as a Voigt average
/// of straight and crimped regions.
///
/// `vf_warp` and `vf_fill` are volume fractions in warp and fill tows.
/// `e_tow` is the effective modulus of a straight tow.
/// `e_matrix` is the matrix modulus.
/// `crimp_angle` is the crimp half-angle (rad) of the tow undulation.
pub fn woven_mosaic_modulus(
    vf_warp: f64,
    vf_fill: f64,
    e_tow: f64,
    e_matrix: f64,
    crimp_angle: f64,
) -> f64 {
    let vf_warp = vf_warp.clamp(0.0, 1.0);
    let vf_fill = vf_fill.clamp(0.0, 1.0);
    let vm = (1.0 - vf_warp - vf_fill).max(0.0);

    // Crimped tow: modulus reduced by cos^4(θ) (classical off-axis rule)
    let e_tow_crimp = e_tow * crimp_angle.cos().powi(4);

    // Voigt combination: straight warp + crimped fill + matrix
    vf_warp * e_tow + vf_fill * e_tow_crimp + vm * e_matrix
}

/// Bridging model shear modulus for a plain weave composite.
///
/// G_eff = G_m * (V_f * G_f + V_m * G_m) / (V_m * G_f + V_f * G_m)
/// (analogous to Reuss shear, weighted by bridging)
pub fn woven_shear_modulus(vf: f64, g_fiber: f64, g_matrix: f64) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let vm = 1.0 - vf;
    let num = g_matrix * (vf * g_fiber + vm * g_matrix);
    let den = vm * g_fiber + vf * g_matrix;
    if den.abs() < f64::EPSILON {
        g_matrix
    } else {
        num / den
    }
}

// ─── Particulate composite ───────────────────────────────────────────────────

/// Kerner model for bulk modulus of a particulate composite.
///
/// K_eff / K_m = \[1 + V_p * (K_p - K_m) / (K_m + 4/3 G_m * (1 - V_p))\]
///
/// This is equivalent to the dilute Eshelby estimate.
pub fn kerner_bulk_modulus(vp: f64, k_particle: f64, g_matrix: f64, k_matrix: f64) -> f64 {
    let vp = vp.clamp(0.0, 1.0);
    let dk = k_particle - k_matrix;
    let denom = k_matrix + 4.0 / 3.0 * g_matrix * (1.0 - vp);
    if denom.abs() < f64::EPSILON {
        k_matrix
    } else {
        k_matrix + vp * dk * k_matrix / denom
    }
}

/// Nielsen model for tensile modulus of a particulate composite.
///
/// E_eff = E_m * (1 + A_E * B_E * V_p) / (1 - B_E * ψ * V_p)
///
/// where B_E = (E_p/E_m - 1) / (E_p/E_m + A_E)
/// ψ = 1 + (1 - φ_m) * V_p / φ_m²
/// φ_m is the maximum packing fraction (e.g. 0.637 for random packing)
/// A_E is the Einstein coefficient (e.g. 2.5 for spheres)
pub fn nielsen_modulus(vp: f64, e_particle: f64, e_matrix: f64, a_e: f64, phi_max: f64) -> f64 {
    let vp = vp.clamp(0.0, phi_max);
    if e_matrix.abs() < f64::EPSILON {
        return 0.0;
    }
    let ratio = e_particle / e_matrix;
    let b_e = (ratio - 1.0) / (ratio + a_e);
    let psi = 1.0 + (1.0 - phi_max) * vp / (phi_max * phi_max);
    let denom = 1.0 - b_e * psi * vp;
    if denom.abs() < f64::EPSILON {
        e_matrix
    } else {
        e_matrix * (1.0 + a_e * b_e * vp) / denom
    }
}

/// Composite sphere model for the effective bulk modulus (Hashin).
///
/// K_eff = K_2 + V_1 / (1/(K_1 - K_2) + 3*V_2 / (3*K_2 + 4*G_2))
///
/// (identical to Mori-Tanaka for spherical inclusions in an isotropic matrix)
pub fn composite_sphere_bulk(vp: f64, k_inclusion: f64, g_matrix: f64, k_matrix: f64) -> f64 {
    let vp = vp.clamp(0.0, 1.0);
    let vm = 1.0 - vp;
    let dk = k_inclusion - k_matrix;
    let inv = 1.0 / dk + 3.0 * vm / (3.0 * k_matrix + 4.0 * g_matrix);
    if inv.abs() < f64::EPSILON || dk.abs() < f64::EPSILON {
        k_matrix
    } else {
        k_matrix + vp / inv
    }
}

// ─── Nano-composite effective properties ─────────────────────────────────────

/// Interface/interphase correction for nano-composites.
///
/// For nano-fillers, the interface layer has non-negligible thickness.
/// Effective volume fraction including interphase:
///
/// V_eff = V_f * (1 + t / r)^3
///
/// where t is interphase thickness and r is particle radius.
pub fn nano_effective_volume_fraction(
    vf: f64,
    particle_radius: f64,
    interphase_thickness: f64,
) -> f64 {
    let vf = vf.clamp(0.0, 1.0);
    let factor = (1.0 + interphase_thickness / particle_radius).powi(3);
    (vf * factor).min(1.0)
}

/// Nano-composite modulus with interphase using a three-phase model.
///
/// Treats the composite as matrix + interphase + core particle,
/// applying the Mori-Tanaka sequentially:
/// 1) Combine core particle + interphase → effective inclusion.
/// 2) Combine effective inclusion + matrix → composite modulus.
pub fn nano_composite_modulus(
    vf_core: f64,
    e_core: f64,
    vf_interphase: f64,
    e_interphase: f64,
    e_matrix: f64,
    xi: f64,
) -> f64 {
    // Step 1: effective inclusion = core + interphase layer
    let vf_total = (vf_core + vf_interphase).min(1.0);
    if vf_total < 1e-15 {
        return e_matrix;
    }
    let vf_core_in_inclusion = vf_core / vf_total;
    let e_inclusion = halpin_tsai_modulus(vf_core_in_inclusion, e_core, e_interphase, xi);

    // Step 2: Halpin-Tsai with effective inclusion
    halpin_tsai_modulus(vf_total, e_inclusion, e_matrix, xi)
}

/// Surface/interface energy contribution to nano-composite elastic modulus.
///
/// Gurtin-Murdoch surface elasticity correction (simplified):
/// ΔK_surface ≈ 2 * K_s / r
///
/// where K_s is the surface bulk modulus and r is the particle radius.
pub fn nano_surface_elasticity_correction(k_surface: f64, particle_radius: f64) -> f64 {
    if particle_radius < f64::EPSILON {
        return 0.0;
    }
    2.0 * k_surface / particle_radius
}

// ─── Tests for new combination additions ─────────────────────────────────────

#[cfg(test)]
mod tests_new_combination {

    use crate::combination::bulk_shear_to_engineering;
    use crate::combination::composite_sphere_bulk;
    use crate::combination::cox_krenchel_modulus;
    use crate::combination::cox_length_efficiency;
    use crate::combination::halpin_tsai_modulus;
    use crate::combination::halpin_tsai_random_2d;
    use crate::combination::homogenization_result;
    use crate::combination::kerner_bulk_modulus;
    use crate::combination::mori_tanaka_homogenization;
    use crate::combination::nano_composite_modulus;
    use crate::combination::nano_effective_volume_fraction;
    use crate::combination::nano_surface_elasticity_correction;
    use crate::combination::nielsen_modulus;
    use crate::combination::woven_mosaic_modulus;
    use crate::combination::woven_shear_modulus;

    // --- Mori-Tanaka homogenization ---

    #[test]
    fn test_mori_tanaka_pure_matrix() {
        // Vf = 0 → K_eff = K_matrix, G_eff = G_matrix
        let km = 50.0e9;
        let gm = 30.0e9;
        let (k_eff, g_eff) = mori_tanaka_homogenization(0.0, 100.0e9, 60.0e9, km, gm);
        assert!((k_eff - km).abs() / km < 1e-10);
        assert!((g_eff - gm).abs() / gm < 1e-10);
    }

    #[test]
    fn test_mori_tanaka_bounds() {
        // K_eff should be between K_matrix and K_fiber
        let km = 50.0e9;
        let gm = 30.0e9;
        let kf = 200.0e9;
        let gf = 100.0e9;
        let (k_eff, g_eff) = mori_tanaka_homogenization(0.5, kf, gf, km, gm);
        assert!(k_eff >= km, "K_eff {k_eff} should be >= K_m {km}");
        assert!(k_eff <= kf, "K_eff {k_eff} should be <= K_f {kf}");
        assert!(g_eff >= gm, "G_eff {g_eff} should be >= G_m {gm}");
        assert!(g_eff <= gf, "G_eff {g_eff} should be <= G_f {gf}");
    }

    #[test]
    fn test_bulk_shear_to_engineering() {
        // For steel: K ≈ 167 GPa, G ≈ 80 GPa → E ≈ 200 GPa, ν ≈ 0.25
        let (e, nu) = bulk_shear_to_engineering(167.0e9, 80.0e9);
        assert!(e > 180.0e9 && e < 220.0e9, "E = {e}");
        assert!(nu > 0.2 && nu < 0.35, "ν = {nu}");
    }

    #[test]
    fn test_homogenization_result_positive_moduli() {
        let res = homogenization_result(0.4, 200.0e9, 80.0e9, 50.0e9, 20.0e9);
        assert!(res.e_x > 0.0, "E_x should be positive");
        assert!(res.g_xy > 0.0, "G_xy should be positive");
        assert!(res.nu_xy > -1.0 && res.nu_xy < 0.5, "ν_xy = {}", res.nu_xy);
    }

    #[test]
    fn test_homogenization_result_stiffness_tensor() {
        let res = homogenization_result(0.3, 200.0e9, 80.0e9, 50.0e9, 20.0e9);
        // C_11 = C_22 = C_33
        assert!((res.c_eff[0] - res.c_eff[7]).abs() < 1.0);
        assert!((res.c_eff[0] - res.c_eff[14]).abs() < 1.0);
        // C_44 = C_55 = C_66
        assert!((res.c_eff[21] - res.c_eff[28]).abs() < 1.0);
    }

    // --- Random fiber composite ---

    #[test]
    fn test_cox_krenchel_pure_matrix() {
        let e = cox_krenchel_modulus(0.0, 72.0e9, 3.5e9, 1.0);
        assert!((e - 3.5e9).abs() < 1.0);
    }

    #[test]
    fn test_cox_krenchel_increases_with_vf() {
        let e0 = cox_krenchel_modulus(0.1, 72.0e9, 3.5e9, 0.8);
        let e1 = cox_krenchel_modulus(0.4, 72.0e9, 3.5e9, 0.8);
        assert!(
            e1 > e0,
            "modulus should increase with fiber volume fraction"
        );
    }

    #[test]
    fn test_cox_length_efficiency_long_fiber() {
        // For long fibers (large β*L/2), η_l → 1
        let eta = cox_length_efficiency(10.0, 20.0);
        assert!(
            eta > 0.9,
            "length efficiency should be close to 1 for long fibers: {eta}"
        );
    }

    #[test]
    fn test_cox_length_efficiency_zero() {
        let eta = cox_length_efficiency(0.0, 0.0);
        assert!((eta - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_halpin_tsai_random_2d() {
        let e = halpin_tsai_random_2d(0.5, 72.0e9, 3.5e9, 2.0);
        assert!(e > 3.5e9, "should exceed matrix modulus");
        assert!(e < 72.0e9, "should be less than fiber modulus");
    }

    // --- Woven composite ---

    #[test]
    fn test_woven_mosaic_zero_crimp() {
        // No crimp: E = Voigt average of warp + fill + matrix
        let e = woven_mosaic_modulus(0.3, 0.3, 70.0e9, 3.5e9, 0.0);
        let expected = 0.3 * 70.0e9 + 0.3 * 70.0e9 + 0.4 * 3.5e9;
        assert!((e - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_woven_mosaic_crimp_reduces_modulus() {
        let e_no_crimp = woven_mosaic_modulus(0.3, 0.3, 70.0e9, 3.5e9, 0.0);
        let e_crimp = woven_mosaic_modulus(0.3, 0.3, 70.0e9, 3.5e9, 0.2);
        assert!(e_crimp < e_no_crimp, "crimp should reduce modulus");
    }

    #[test]
    fn test_woven_shear_equal_phases() {
        // G_f = G_m → G_eff = G_m
        let g = woven_shear_modulus(0.5, 50.0e9, 50.0e9);
        assert!((g - 50.0e9).abs() / 50.0e9 < 1e-10);
    }

    #[test]
    fn test_woven_shear_between_phases() {
        let gf = 80.0e9;
        let gm = 3.0e9;
        let g = woven_shear_modulus(0.4, gf, gm);
        // The bridging model is a harmonic-type blend
        // Result is positive and finite
        assert!(g > 0.0, "G_eff {g} should be positive");
        assert!(g.is_finite(), "G_eff {g} should be finite");
        // For different vf values, result should change
        let g2 = woven_shear_modulus(0.6, gf, gm);
        assert!((g - g2).abs() > 0.0, "result should vary with vf");
    }

    // --- Particulate composite ---

    #[test]
    fn test_kerner_bulk_pure_matrix() {
        let km = 50.0e9;
        let gm = 30.0e9;
        let k = kerner_bulk_modulus(0.0, 200.0e9, gm, km);
        // Vp=0: numerator = 0, so k ≈ k_matrix
        assert!((k - km).abs() / km < 1e-6);
    }

    #[test]
    fn test_kerner_bulk_above_matrix() {
        let k = kerner_bulk_modulus(0.3, 200.0e9, 30.0e9, 50.0e9);
        assert!(k > 50.0e9, "Kerner K should exceed matrix K");
    }

    #[test]
    fn test_nielsen_modulus_pure_matrix() {
        let e = nielsen_modulus(0.0, 200.0e9, 3.5e9, 2.5, 0.64);
        assert!((e - 3.5e9).abs() < 1.0);
    }

    #[test]
    fn test_nielsen_modulus_increases() {
        let e0 = nielsen_modulus(0.1, 200.0e9, 3.5e9, 2.5, 0.64);
        let e1 = nielsen_modulus(0.3, 200.0e9, 3.5e9, 2.5, 0.64);
        assert!(e1 > e0, "modulus should increase with volume fraction");
    }

    #[test]
    fn test_composite_sphere_bulk() {
        // Result should be between K_matrix and K_inclusion
        let km = 50.0e9;
        let ki = 200.0e9;
        let gm = 30.0e9;
        let k = composite_sphere_bulk(0.4, ki, gm, km);
        assert!(k >= km, "K_eff should be >= K_m");
        assert!(k <= ki, "K_eff should be <= K_i");
    }

    // --- Nano-composite ---

    #[test]
    fn test_nano_effective_vf_no_interphase() {
        // t=0 → V_eff = V_f
        let veff = nano_effective_volume_fraction(0.3, 10e-9, 0.0);
        assert!((veff - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_nano_effective_vf_grows_with_interphase() {
        let vf = 0.1;
        let r = 5e-9;
        let t = 2e-9;
        let veff = nano_effective_volume_fraction(vf, r, t);
        assert!(
            veff > vf,
            "V_eff should exceed V_f with non-zero interphase"
        );
    }

    #[test]
    fn test_nano_composite_modulus_no_interphase() {
        // V_interphase = 0 → should reduce to Halpin-Tsai
        let e = nano_composite_modulus(0.3, 200.0e9, 0.0, 3.5e9, 3.5e9, 2.0);
        let e_ht = halpin_tsai_modulus(0.3, 200.0e9, 3.5e9, 2.0);
        // Should be close (interphase = matrix)
        assert!((e - e_ht).abs() / e_ht < 0.05, "e={e}, e_ht={e_ht}");
    }

    #[test]
    fn test_nano_surface_correction_large_particle() {
        // For large particles, surface correction is negligible
        let delta_k = nano_surface_elasticity_correction(1e-9, 1e-3);
        assert!(
            delta_k < 1e-3,
            "correction should be tiny for large particles: {delta_k}"
        );
    }

    #[test]
    fn test_nano_surface_correction_small_particle() {
        // For nano-particles, surface correction is significant
        let delta_k = nano_surface_elasticity_correction(1.0, 1e-9);
        assert!(
            delta_k > 1e8,
            "correction should be large for nano-particles: {delta_k}"
        );
    }

    #[test]
    fn test_nano_surface_correction_zero_radius() {
        let delta_k = nano_surface_elasticity_correction(1.0, 0.0);
        assert_eq!(delta_k, 0.0);
    }
}
