// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended homogenization methods: asymptotic homogenization,
//! Voigt/Reuss scalar bounds, Hashin-Shtrikman scalar form,
//! fiber-reinforced composites, Classical Laminate Theory (CLT),
//! and additional composite moduli estimates.

#![allow(dead_code)]

use super::advanced::{mori_tanaka_bulk_modulus, mori_tanaka_shear_modulus};

// ============================================================================
// Asymptotic Homogenization (1D / scalar)
// ============================================================================

/// Asymptotic homogenization effective coefficient for a 1-D periodic cell.
///
/// For a 1-D cell Y = \[0, L\] with conductivity / stiffness a(y), the
/// effective coefficient is the harmonic mean:
///
///   a* = L / ∫₀ᴸ (1/a(y)) dy  (Reuss formula for 1-D)
///
/// # Arguments
/// * `a_vals`   – values of a(y) at uniformly spaced quadrature points
/// * `cell_len` – length L of the unit cell
#[allow(dead_code)]
pub fn asymptotic_homogenization_1d(a_vals: &[f64], cell_len: f64) -> f64 {
    if a_vals.is_empty() {
        return 0.0;
    }
    let n = a_vals.len() as f64;
    let dy = cell_len / n;
    let integral_inv: f64 = a_vals.iter().map(|&a| dy / a).sum();
    cell_len / integral_inv
}

/// Asymptotic homogenization effective coefficient via Voigt (upper) bound.
///
/// For a 1-D cell: a*_Voigt = (1/L) ∫₀ᴸ a(y) dy  (arithmetic mean).
#[allow(dead_code)]
pub fn asymptotic_homogenization_1d_voigt(a_vals: &[f64], cell_len: f64) -> f64 {
    if a_vals.is_empty() {
        return 0.0;
    }
    let n = a_vals.len() as f64;
    let dy = cell_len / n;
    let integral: f64 = a_vals.iter().map(|&a| a * dy).sum();
    integral / cell_len
}

/// Asymptotic homogenization of a 2-D composite via numerical integration
/// on a rectangular unit cell with two phases.
///
/// Returns the homogenized conductivity tensor (2×2 as \[f64; 4\] row-major)
/// using a simplified averaging approach (phase-averaged Voigt / Reuss).
///
/// # Arguments
/// * `phase_map` – n×n grid of phase index (0 or 1)
/// * `k0`        – conductivity of phase 0
/// * `k1`        – conductivity of phase 1
#[allow(dead_code)]
pub fn asymptotic_homogenization_2d_conductivity(
    phase_map: &[Vec<usize>],
    k0: f64,
    k1: f64,
) -> [f64; 4] {
    let n_rows = phase_map.len();
    if n_rows == 0 {
        return [0.0; 4];
    }
    let n_cols = phase_map[0].len();
    let total = (n_rows * n_cols) as f64;
    let f1: f64 = phase_map
        .iter()
        .flat_map(|row| row.iter())
        .filter(|&&p| p == 1)
        .count() as f64
        / total;
    let f0 = 1.0 - f1;

    // Simple Voigt average for x-direction, Reuss for y-direction
    let k_xx = f0 * k0 + f1 * k1;
    let k_yy = if (f0 / k0 + f1 / k1).abs() < 1e-60 {
        0.0
    } else {
        1.0 / (f0 / k0 + f1 / k1)
    };
    [k_xx, 0.0, 0.0, k_yy]
}

// ============================================================================
// Voigt / Reuss scalar bounds (extended)
// ============================================================================

/// Voigt (upper) bound on Young's modulus for N phases (rule of mixtures).
///
/// E_V = Σ f_i E_i
#[allow(dead_code)]
pub fn voigt_modulus(moduli: &[f64], fractions: &[f64]) -> f64 {
    moduli
        .iter()
        .zip(fractions.iter())
        .map(|(e, f)| e * f)
        .sum()
}

/// Reuss (lower) bound on Young's modulus for N phases.
///
/// 1/E_R = Σ f_i / E_i
#[allow(dead_code)]
pub fn reuss_modulus(moduli: &[f64], fractions: &[f64]) -> f64 {
    let inv_sum: f64 = moduli
        .iter()
        .zip(fractions.iter())
        .map(|(e, f)| f / e)
        .sum();
    if inv_sum.abs() < 1e-60 {
        0.0
    } else {
        1.0 / inv_sum
    }
}

/// Hill average Young's modulus: E_H = (E_V + E_R) / 2.
#[allow(dead_code)]
pub fn hill_modulus(moduli: &[f64], fractions: &[f64]) -> f64 {
    (voigt_modulus(moduli, fractions) + reuss_modulus(moduli, fractions)) / 2.0
}

/// Voigt shear modulus: G_V = Σ f_i G_i.
#[allow(dead_code)]
pub fn voigt_shear_modulus(shear_moduli: &[f64], fractions: &[f64]) -> f64 {
    shear_moduli
        .iter()
        .zip(fractions.iter())
        .map(|(g, f)| g * f)
        .sum()
}

/// Reuss shear modulus: 1/G_R = Σ f_i / G_i.
#[allow(dead_code)]
pub fn reuss_shear_modulus(shear_moduli: &[f64], fractions: &[f64]) -> f64 {
    let inv: f64 = shear_moduli
        .iter()
        .zip(fractions.iter())
        .map(|(g, f)| f / g)
        .sum();
    if inv.abs() < 1e-60 { 0.0 } else { 1.0 / inv }
}

// ============================================================================
// Hashin-Shtrikman bounds (scalar form)
// ============================================================================

/// Hashin-Shtrikman lower bound for bulk modulus (two-phase, K1 < K2, G1 < G2).
///
/// K_HS⁻ = K1 + f2 / (1/(K2-K1) + 3*f1/(3*K1+4*G1))
#[allow(dead_code)]
pub fn hs_bulk_lower_bound(k1: f64, g1: f64, k2: f64, f1: f64, f2: f64) -> f64 {
    let denom = 1.0 / (k2 - k1) + 3.0 * f1 / (3.0 * k1 + 4.0 * g1);
    if denom.abs() < 1e-60 {
        k1
    } else {
        k1 + f2 / denom
    }
}

/// Hashin-Shtrikman upper bound for bulk modulus.
///
/// K_HS⁺ = K2 + f1 / (1/(K1-K2) + 3*f2/(3*K2+4*G2))
#[allow(dead_code)]
pub fn hs_bulk_upper_bound(k2: f64, g2: f64, k1: f64, f1: f64, f2: f64) -> f64 {
    let denom = 1.0 / (k1 - k2) + 3.0 * f2 / (3.0 * k2 + 4.0 * g2);
    if denom.abs() < 1e-60 {
        k2
    } else {
        k2 + f1 / denom
    }
}

/// Hashin-Shtrikman lower bound for shear modulus (two-phase).
///
/// G_HS⁻ = G1 + f2 / (1/(G2-G1) + 6*f1*(K1+2*G1)/(5*G1*(3*K1+4*G1)))
#[allow(dead_code)]
pub fn hs_shear_lower_bound(k1: f64, g1: f64, g2: f64, f1: f64, f2: f64) -> f64 {
    if (g2 - g1).abs() < 1e-60 {
        return g1;
    }
    let alpha1 = 6.0 * f1 * (k1 + 2.0 * g1) / (5.0 * g1 * (3.0 * k1 + 4.0 * g1));
    let denom = 1.0 / (g2 - g1) + alpha1;
    if denom.abs() < 1e-60 {
        g1
    } else {
        g1 + f2 / denom
    }
}

/// Hashin-Shtrikman upper bound for shear modulus (two-phase).
///
/// G_HS⁺ = G2 + f1 / (1/(G1-G2) + 6*f2*(K2+2*G2)/(5*G2*(3*K2+4*G2)))
#[allow(dead_code)]
pub fn hs_shear_upper_bound(k2: f64, g2: f64, g1: f64, f1: f64, f2: f64) -> f64 {
    if (g1 - g2).abs() < 1e-60 {
        return g2;
    }
    let alpha2 = 6.0 * f2 * (k2 + 2.0 * g2) / (5.0 * g2 * (3.0 * k2 + 4.0 * g2));
    let denom = 1.0 / (g1 - g2) + alpha2;
    if denom.abs() < 1e-60 {
        g2
    } else {
        g2 + f1 / denom
    }
}

// ============================================================================
// Fiber-Reinforced Composites (transversely isotropic)
// ============================================================================

/// Engineering constants for a transversely isotropic fiber-reinforced composite.
#[allow(dead_code)]
pub struct FiberCompositeConstants {
    /// Longitudinal Young's modulus E₁ (fiber direction).
    pub e1: f64,
    /// Transverse Young's modulus E₂.
    pub e2: f64,
    /// In-plane Poisson's ratio ν₁₂.
    pub nu12: f64,
    /// Out-of-plane Poisson's ratio ν₂₃.
    pub nu23: f64,
    /// In-plane shear modulus G₁₂.
    pub g12: f64,
    /// Transverse shear modulus G₂₃.
    pub g23: f64,
}

/// Compute fiber-reinforced composite constants via the rule of mixtures
/// and Halpin-Tsai micromechanics.
///
/// # Arguments
/// * `e_f`  – fiber Young's modulus
/// * `e_m`  – matrix Young's modulus
/// * `nu_f` – fiber Poisson's ratio
/// * `nu_m` – matrix Poisson's ratio
/// * `g_f`  – fiber shear modulus
/// * `g_m`  – matrix shear modulus
/// * `vf`   – fiber volume fraction
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn fiber_composite_rom(
    e_f: f64,
    e_m: f64,
    nu_f: f64,
    nu_m: f64,
    g_f: f64,
    g_m: f64,
    vf: f64,
) -> FiberCompositeConstants {
    let vm = 1.0 - vf;

    // Longitudinal: rule of mixtures
    let e1 = vf * e_f + vm * e_m;

    // Transverse: Halpin-Tsai with ξ = 2 (circular fibers)
    let xi_e = 2.0_f64;
    let eta_e = (e_f / e_m - 1.0) / (e_f / e_m + xi_e);
    let e2 = e_m * (1.0 + xi_e * eta_e * vf) / (1.0 - eta_e * vf);

    // Poisson's ratios
    let nu12 = vf * nu_f + vm * nu_m;

    // ν₂₃ estimate via Christensen (simplified)
    let k_f = e_f / (3.0 * (1.0 - 2.0 * nu_f));
    let k_m = e_m / (3.0 * (1.0 - 2.0 * nu_m));
    let k12 = k_m + vf / (1.0 / (k_f - k_m) + vm / (k_m + g_m / 3.0));
    // E₂ = 4 k12 g23 / (k12 + g23 + k12 g23 nu12² / e1) – simplified
    let nu23_approx = e2 / (2.0 * (k12 + e2 / (1.0 + nu12)).max(1e-60)) - 1.0;
    let nu23 = nu23_approx.clamp(-0.5, 0.5);

    // In-plane shear: Halpin-Tsai with ξ = 1
    let xi_g = 1.0_f64;
    let eta_g = (g_f / g_m - 1.0) / (g_f / g_m + xi_g);
    let g12 = g_m * (1.0 + xi_g * eta_g * vf) / (1.0 - eta_g * vf);

    // Transverse shear: Halpin-Tsai with ξ = 1
    let eta_g23 = (g_f / g_m - 1.0) / (g_f / g_m + 1.0);
    let g23 = g_m * (1.0 + eta_g23 * vf) / (1.0 - eta_g23 * vf);

    let _ = k12; // used above
    FiberCompositeConstants {
        e1,
        e2,
        nu12,
        nu23,
        g12,
        g23,
    }
}

/// Longitudinal tensile strength of a fiber composite (rule of mixtures).
///
/// σ₁* = Vf * σ_fu + (1 - Vf) * σ_mu
///
/// where σ_fu is fiber strength and σ_mu is matrix stress at fiber failure strain.
#[allow(dead_code)]
pub fn fiber_strength_longitudinal(sigma_fu: f64, sigma_mu: f64, vf: f64) -> f64 {
    vf * sigma_fu + (1.0 - vf) * sigma_mu
}

/// Critical fiber volume fraction below which matrix failure dominates.
///
/// V_crit = (σ_mu - σ_mu_star) / (σ_fu + σ_mu - σ_mu_star)
#[allow(dead_code)]
pub fn critical_fiber_volume_fraction(sigma_mu: f64, sigma_mu_star: f64, sigma_fu: f64) -> f64 {
    let denom = sigma_fu + sigma_mu - sigma_mu_star;
    if denom.abs() < 1e-60 {
        0.0
    } else {
        (sigma_mu - sigma_mu_star) / denom
    }
}

// ============================================================================
// Rule of Mixtures (extended)
// ============================================================================

/// Rule of mixtures for density: ρ* = Σ f_i ρ_i.
#[allow(dead_code)]
pub fn rule_of_mixtures_density(densities: &[f64], fractions: &[f64]) -> f64 {
    densities
        .iter()
        .zip(fractions.iter())
        .map(|(rho, f)| rho * f)
        .sum()
}

/// Rule of mixtures for thermal conductivity (Voigt = parallel rule).
#[allow(dead_code)]
pub fn rule_of_mixtures_conductivity(lambdas: &[f64], fractions: &[f64]) -> f64 {
    lambdas
        .iter()
        .zip(fractions.iter())
        .map(|(l, f)| l * f)
        .sum()
}

/// Inverse rule of mixtures for compliance C_R = Σ f_i / E_i (Reuss).
#[allow(dead_code)]
pub fn inverse_rule_of_mixtures(moduli: &[f64], fractions: &[f64]) -> f64 {
    let s: f64 = moduli
        .iter()
        .zip(fractions.iter())
        .map(|(e, f)| f / e)
        .sum();
    if s.abs() < 1e-60 { 0.0 } else { 1.0 / s }
}

/// Rule of mixtures for CTE: α* = Σ f_i α_i E_i / E*.
///
/// Uses the Schapery model: α* = (Σ f_i α_i E_i) / (Σ f_i E_i).
#[allow(dead_code)]
pub fn rule_of_mixtures_cte(alphas: &[f64], moduli: &[f64], fractions: &[f64]) -> f64 {
    let num: f64 = alphas
        .iter()
        .zip(moduli.iter())
        .zip(fractions.iter())
        .map(|((a, e), f)| a * e * f)
        .sum();
    let den: f64 = moduli
        .iter()
        .zip(fractions.iter())
        .map(|(e, f)| e * f)
        .sum();
    if den.abs() < 1e-60 { 0.0 } else { num / den }
}

// ============================================================================
// Classical Laminate Theory (CLT)
// ============================================================================

/// A single lamina (ply) in a laminate stack.
#[allow(dead_code)]
pub struct LaminaLayer {
    /// Young's modulus in the fiber direction E₁.
    pub e1: f64,
    /// Young's modulus transverse to fibers E₂.
    pub e2: f64,
    /// In-plane shear modulus G₁₂.
    pub g12: f64,
    /// In-plane Poisson's ratio ν₁₂.
    pub nu12: f64,
    /// Ply thickness \[m\].
    pub thickness: f64,
    /// Fiber orientation angle θ \[radians\].
    pub theta: f64,
}

/// Compute the reduced stiffness matrix \[Q\] for a single lamina in its
/// principal material coordinate system.
///
/// Returns the 3×3 in-plane stiffness matrix Q as `[[f64; 3\]; 3]`.
#[allow(dead_code)]
pub fn lamina_reduced_stiffness(l: &LaminaLayer) -> [[f64; 3]; 3] {
    let nu21 = l.nu12 * l.e2 / l.e1;
    let denom = 1.0 - l.nu12 * nu21;
    if denom.abs() < 1e-60 {
        return [[0.0; 3]; 3];
    }
    let q11 = l.e1 / denom;
    let q12 = l.nu12 * l.e2 / denom;
    let q22 = l.e2 / denom;
    let q66 = l.g12;
    [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
}

/// Rotate the reduced stiffness matrix \[Q\] by angle θ to the laminate
/// coordinate system, producing the transformed stiffness Q̄.
///
/// Returns the 3×3 transformed reduced stiffness matrix.
#[allow(dead_code)]
pub fn transform_reduced_stiffness(q: &[[f64; 3]; 3], theta: f64) -> [[f64; 3]; 3] {
    let c = theta.cos();
    let s = theta.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;

    let q11 = q[0][0];
    let q12 = q[0][1];
    let q22 = q[1][1];
    let q66 = q[2][2];

    let q_bar_11 = q11 * c2 * c2 + 2.0 * (q12 + 2.0 * q66) * c2 * s2 + q22 * s2 * s2;
    let q_bar_12 = (q11 + q22 - 4.0 * q66) * c2 * s2 + q12 * (c2 * c2 + s2 * s2);
    let q_bar_22 = q11 * s2 * s2 + 2.0 * (q12 + 2.0 * q66) * c2 * s2 + q22 * c2 * c2;
    let q_bar_16 = (q11 - q12 - 2.0 * q66) * c2 * cs + (q12 - q22 + 2.0 * q66) * s2 * cs;
    let q_bar_26 = (q11 - q12 - 2.0 * q66) * cs * s2 + (q12 - q22 + 2.0 * q66) * cs * c2;
    let q_bar_66 =
        (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2 * s2 + q66 * (c2 * c2 + s2 * s2 - 2.0 * c2 * s2);

    [
        [q_bar_11, q_bar_12, q_bar_16],
        [q_bar_12, q_bar_22, q_bar_26],
        [q_bar_16, q_bar_26, q_bar_66],
    ]
}

/// Result of a Classical Laminate Theory analysis.
///
/// Contains the \[A\], \[B\], \[D\] matrices of the laminate.
#[allow(dead_code)]
pub struct CltResult {
    /// Extensional stiffness matrix A (3×3) \[N/m\].
    pub a_mat: [[f64; 3]; 3],
    /// Coupling stiffness matrix B (3×3) \[N\].
    pub b_mat: [[f64; 3]; 3],
    /// Bending stiffness matrix D (3×3) \[N·m\].
    pub d_mat: [[f64; 3]; 3],
    /// Total laminate thickness \[m\].
    pub total_thickness: f64,
    /// Effective in-plane Young's modulus Ex \[Pa\].
    pub ex: f64,
    /// Effective in-plane Young's modulus Ey \[Pa\].
    pub ey: f64,
    /// Effective in-plane shear modulus Gxy \[Pa\].
    pub gxy: f64,
}

/// Perform Classical Laminate Theory (CLT) analysis on a stack of plies.
///
/// Computes the A, B, D matrices and effective laminate engineering constants.
///
/// # Arguments
/// * `layers` – ordered list of lamina (from bottom z₁ to top z_n)
#[allow(dead_code)]
pub fn clt_analysis(layers: &[LaminaLayer]) -> CltResult {
    let total_thickness: f64 = layers.iter().map(|l| l.thickness).sum();

    // z-coordinates of ply interfaces (mid-plane at z = 0)
    let mut z_coords = vec![0.0_f64; layers.len() + 1];
    z_coords[0] = -total_thickness / 2.0;
    for (k, l) in layers.iter().enumerate() {
        z_coords[k + 1] = z_coords[k] + l.thickness;
    }

    let mut a_mat = [[0.0_f64; 3]; 3];
    let mut b_mat = [[0.0_f64; 3]; 3];
    let mut d_mat = [[0.0_f64; 3]; 3];

    for (k, layer) in layers.iter().enumerate() {
        let q = lamina_reduced_stiffness(layer);
        let q_bar = transform_reduced_stiffness(&q, layer.theta);

        let zk = z_coords[k];
        let zk1 = z_coords[k + 1];
        let dz1 = zk1 - zk;
        let dz2 = (zk1 * zk1 - zk * zk) / 2.0;
        let dz3 = (zk1 * zk1 * zk1 - zk * zk * zk) / 3.0;

        for i in 0..3 {
            for j in 0..3 {
                a_mat[i][j] += q_bar[i][j] * dz1;
                b_mat[i][j] += q_bar[i][j] * dz2;
                d_mat[i][j] += q_bar[i][j] * dz3;
            }
        }
    }

    // Effective laminate constants (for symmetric laminates B ≈ 0)
    let h = total_thickness;
    let ex = if a_mat[0][0].abs() > 1e-60 {
        (a_mat[0][0] * a_mat[1][1] - a_mat[0][1] * a_mat[0][1]) / (a_mat[1][1] * h)
    } else {
        0.0
    };
    let ey = if a_mat[1][1].abs() > 1e-60 {
        (a_mat[0][0] * a_mat[1][1] - a_mat[0][1] * a_mat[0][1]) / (a_mat[0][0] * h)
    } else {
        0.0
    };
    let gxy = a_mat[2][2] / h;

    CltResult {
        a_mat,
        b_mat,
        d_mat,
        total_thickness,
        ex,
        ey,
        gxy,
    }
}

/// Compute the in-plane Poisson's ratio ν_xy from the CLT A-matrix.
///
/// ν_xy = A12 / A22 for a symmetric laminate.
#[allow(dead_code)]
pub fn clt_poisson_xy(a: &[[f64; 3]; 3]) -> f64 {
    if a[1][1].abs() < 1e-60 {
        0.0
    } else {
        a[0][1] / a[1][1]
    }
}

/// Check if a laminate is symmetric (B ≈ 0) within a tolerance.
#[allow(dead_code)]
pub fn laminate_is_symmetric(b: &[[f64; 3]; 3], tol: f64) -> bool {
    b.iter().all(|row| row.iter().all(|&v| v.abs() < tol))
}

// ============================================================================
// Composite moduli (additional bounds and estimates)
// ============================================================================

/// Walpole's bounds for the effective bulk modulus (generalisation of H-S).
///
/// For a two-phase composite with phases ordered K1 < K2:
/// Returns (K_lower, K_upper).
#[allow(dead_code)]
pub fn walpole_bulk_bounds(k1: f64, g1: f64, k2: f64, g2: f64, f2: f64) -> (f64, f64) {
    let f1 = 1.0 - f2;
    let k_lower = hs_bulk_lower_bound(k1, g1, k2, f1, f2);
    let k_upper = hs_bulk_upper_bound(k2, g2, k1, f1, f2);
    (k_lower, k_upper)
}

/// Self-consistent estimate for bulk modulus of a two-phase composite.
///
/// Solves implicitly; here we use the explicit approximation:
/// K* ≈ (f1*K1*(3K2+4G2) + f2*K2*(3K1+4G1)) / (f1*(3K2+4G2) + f2*(3K1+4G1))
#[allow(dead_code)]
pub fn self_consistent_bulk(k1: f64, g1: f64, k2: f64, g2: f64, f1: f64, f2: f64) -> f64 {
    let w1 = 3.0 * k2 + 4.0 * g2;
    let w2 = 3.0 * k1 + 4.0 * g1;
    let num = f1 * k1 * w1 + f2 * k2 * w2;
    let den = f1 * w1 + f2 * w2;
    if den.abs() < 1e-60 { 0.0 } else { num / den }
}

/// Differential scheme (DS) estimate for effective bulk modulus.
///
/// Integrates ODE numerically from K_m to K* by adding inclusions
/// incrementally.  Uses n_steps Euler increments.
#[allow(dead_code)]
pub fn differential_scheme_bulk(k_m: f64, g_m: f64, k_i: f64, f_i: f64, n_steps: usize) -> f64 {
    if n_steps == 0 {
        return k_m;
    }
    let df = f_i / n_steps as f64;
    let mut k = k_m;
    let mut c = 0.0_f64; // current concentration
    for _ in 0..n_steps {
        let c_new = c + df;
        let f_frac = df / (1.0 - c);
        let p = 3.0 * k + 4.0 * g_m;
        let dk = if p.abs() < 1e-60 {
            0.0
        } else {
            (k_i - k) * f_frac * (3.0 * k + 4.0 * g_m) / p
        };
        k += dk;
        c = c_new;
    }
    k
}

/// Eshelby sphere inclusion factor for bulk modulus (κ).
///
/// κ = (K_i - K_m) / (K_i + 4/3 G_m)
#[allow(dead_code)]
pub fn eshelby_bulk_factor(k_i: f64, k_m: f64, g_m: f64) -> f64 {
    let denom = k_i + 4.0 / 3.0 * g_m;
    if denom.abs() < 1e-60 {
        0.0
    } else {
        (k_i - k_m) / denom
    }
}

/// Mori-Tanaka estimate for effective Young's modulus of a two-phase composite.
///
/// Uses the Mori-Tanaka bulk/shear to construct E* from Kfrom and G*.
#[allow(dead_code)]
pub fn mori_tanaka_young(k_m: f64, g_m: f64, k_i: f64, g_i: f64, f_i: f64) -> f64 {
    let k_star = mori_tanaka_bulk_modulus(k_m, g_m, k_i, f_i);
    let g_star = mori_tanaka_shear_modulus(k_m, g_m, g_i, f_i);
    9.0 * k_star * g_star / (3.0 * k_star + g_star)
}

/// Mori-Tanaka estimate for effective Poisson's ratio.
#[allow(dead_code)]
pub fn mori_tanaka_poisson(k_m: f64, g_m: f64, k_i: f64, g_i: f64, f_i: f64) -> f64 {
    let k_star = mori_tanaka_bulk_modulus(k_m, g_m, k_i, f_i);
    let g_star = mori_tanaka_shear_modulus(k_m, g_m, g_i, f_i);
    (3.0 * k_star - 2.0 * g_star) / (6.0 * k_star + 2.0 * g_star)
}

// ============================================================================
// Periodic Boundary Conditions helpers (CLT / composite)
// ============================================================================

/// Compute the average strain in a representative volume element (RVE).
///
/// Returns `[ε_xx, ε_yy, γ_xy]` from the area-weighted phase strains.
///
/// # Arguments
/// * `strains`   – strain vector per phase `[ε_xx, ε_yy, γ_xy]`
/// * `fractions` – volume fraction per phase
#[allow(dead_code)]
pub fn rve_average_strain(strains: &[[f64; 3]], fractions: &[f64]) -> [f64; 3] {
    let n = strains.len().min(fractions.len());
    let mut avg = [0.0_f64; 3];
    for i in 0..n {
        for k in 0..3 {
            avg[k] += fractions[i] * strains[i][k];
        }
    }
    avg
}

/// Compute the average stress in an RVE.
#[allow(dead_code)]
pub fn rve_average_stress(stresses: &[[f64; 3]], fractions: &[f64]) -> [f64; 3] {
    rve_average_strain(stresses, fractions)
}

/// Hill–Mandel energy condition check for in-plane problem.
///
/// Returns |<σ:ε> - <σ>:`ε`|.
#[allow(dead_code)]
pub fn hill_mandel_inplane(stresses: &[[f64; 3]], strains: &[[f64; 3]], fractions: &[f64]) -> f64 {
    let n = stresses.len().min(strains.len()).min(fractions.len());
    // Volume-averaged energy <σ:ε>
    let mut avg_energy = 0.0_f64;
    for i in 0..n {
        let se: f64 = stresses[i]
            .iter()
            .zip(strains[i].iter())
            .map(|(s, e)| s * e)
            .sum();
        avg_energy += fractions[i] * se;
    }
    // Energy of averages <σ>:<ε>
    let avg_s = rve_average_stress(stresses, fractions);
    let avg_e = rve_average_strain(strains, fractions);
    let energy_of_avg: f64 = avg_s.iter().zip(avg_e.iter()).map(|(s, e)| s * e).sum();

    (avg_energy - energy_of_avg).abs()
}

// ============================================================================
// Additional composite material functions
// ============================================================================

/// Puck's inter-fiber failure criterion (simplified, mode A).
///
/// Returns the failure index FI.  FI ≥ 1 implies failure.
///
/// FI = sqrt((τ_21/S21)² + (σ_2/R_T)²) + p_T * σ_2 / S21
///
/// where R_T is the transverse tensile strength, S21 is the in-plane shear
/// strength, and p_T is the inclination parameter for mode A (≈ 0.2–0.35).
#[allow(dead_code)]
pub fn puck_inter_fiber_failure(sigma_2: f64, tau_21: f64, r_t: f64, s21: f64, p_t: f64) -> f64 {
    if s21.abs() < 1e-60 || r_t.abs() < 1e-60 {
        return 0.0;
    }
    ((tau_21 / s21).powi(2) + (sigma_2 / r_t).powi(2)).sqrt() + p_t * sigma_2 / s21
}

/// Tsai-Hill failure criterion for an orthotropic ply.
///
/// Returns failure index.  Failure when FI ≥ 1.
///
/// FI = (σ₁/X)² - σ₁σ₂/X² + (σ₂/Y)² + (τ₁₂/S)²
#[allow(dead_code)]
pub fn tsai_hill_criterion(
    sigma_1: f64,
    sigma_2: f64,
    tau_12: f64,
    x_strength: f64,
    y_strength: f64,
    s_strength: f64,
) -> f64 {
    let t1 = (sigma_1 / x_strength).powi(2);
    let t2 = sigma_1 * sigma_2 / (x_strength * x_strength);
    let t3 = (sigma_2 / y_strength).powi(2);
    let t4 = (tau_12 / s_strength).powi(2);
    t1 - t2 + t3 + t4
}

/// Maximum stress failure criterion: returns true if any component exceeds its limit.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn max_stress_failure(
    sigma_1: f64,
    sigma_2: f64,
    tau_12: f64,
    x_t: f64,
    x_c: f64,
    y_t: f64,
    y_c: f64,
    s: f64,
) -> bool {
    sigma_1 > x_t || sigma_1 < -x_c || sigma_2 > y_t || sigma_2 < -y_c || tau_12.abs() > s
}

/// Compute the thermal residual stresses in a laminate after cool-down.
///
/// Simple approximation: σ_r = E* * (α_m - α_f) * ΔT for each phase.
#[allow(dead_code)]
pub fn thermal_residual_stress(e_eff: f64, alpha_m: f64, alpha_f: f64, delta_t: f64) -> f64 {
    e_eff * (alpha_m - alpha_f) * delta_t
}

#[cfg(test)]
mod tests_homogenization_extended {
    use super::*;

    // ── Asymptotic homogenization tests ──────────────────────────────────────

    #[test]
    fn test_asymptotic_1d_homogeneous() {
        // Uniform material: harmonic mean = value
        let a = vec![100.0_f64; 10];
        let ah = asymptotic_homogenization_1d(&a, 1.0);
        assert!((ah - 100.0).abs() < 1e-8, "homogeneous: ah = {ah}");
    }

    #[test]
    fn test_asymptotic_1d_voigt_homogeneous() {
        let a = vec![50.0_f64; 10];
        let av = asymptotic_homogenization_1d_voigt(&a, 1.0);
        assert!((av - 50.0).abs() < 1e-8, "homogeneous Voigt: av = {av}");
    }

    #[test]
    fn test_asymptotic_1d_voigt_ge_reuss() {
        let a = vec![100.0_f64, 200.0, 100.0, 200.0];
        let reuss = asymptotic_homogenization_1d(&a, 1.0);
        let voigt = asymptotic_homogenization_1d_voigt(&a, 1.0);
        assert!(voigt >= reuss, "Voigt={voigt} should be >= Reuss={reuss}");
    }

    #[test]
    fn test_asymptotic_1d_empty() {
        let result = asymptotic_homogenization_1d(&[], 1.0);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_asymptotic_2d_conductivity_homogeneous() {
        let n = 4;
        let map: Vec<Vec<usize>> = vec![vec![0; n]; n];
        let k = asymptotic_homogenization_2d_conductivity(&map, 1.0, 2.0);
        // f1 = 0: k_xx = k0 = 1.0, k_yy = k0 = 1.0
        assert!((k[0] - 1.0).abs() < 1e-10, "k_xx = {}", k[0]);
    }

    // ── Voigt / Reuss scalar bounds ───────────────────────────────────────────

    #[test]
    fn test_voigt_modulus_single_phase() {
        let e = voigt_modulus(&[200e9], &[1.0]);
        assert!((e - 200e9).abs() < 1.0, "single phase Voigt: {e}");
    }

    #[test]
    fn test_reuss_modulus_single_phase() {
        let e = reuss_modulus(&[200e9], &[1.0]);
        assert!((e - 200e9).abs() < 1.0, "single phase Reuss: {e}");
    }

    #[test]
    fn test_voigt_ge_reuss() {
        let moduli = vec![70e9_f64, 200e9];
        let fracs = vec![0.4_f64, 0.6];
        let ev = voigt_modulus(&moduli, &fracs);
        let er = reuss_modulus(&moduli, &fracs);
        assert!(ev >= er, "Voigt={ev} >= Reuss={er}");
    }

    #[test]
    fn test_hill_modulus_between() {
        let moduli = vec![70e9_f64, 200e9];
        let fracs = vec![0.4_f64, 0.6];
        let ev = voigt_modulus(&moduli, &fracs);
        let er = reuss_modulus(&moduli, &fracs);
        let eh = hill_modulus(&moduli, &fracs);
        assert!(
            eh >= er && eh <= ev,
            "Hill={eh} must be between Reuss={er} and Voigt={ev}"
        );
    }

    #[test]
    fn test_voigt_shear_two_phases() {
        let g = vec![26e9_f64, 77e9];
        let f = vec![0.5_f64, 0.5];
        let gv = voigt_shear_modulus(&g, &f);
        assert!((gv - 51.5e9).abs() / 51.5e9 < 1e-10, "gv = {gv}");
    }

    // ── Hashin-Shtrikman bounds ───────────────────────────────────────────────

    #[test]
    fn test_hs_bulk_lower_le_voigt() {
        let k1 = 50e9_f64;
        let g1 = 30e9_f64;
        let k2 = 150e9_f64;
        let f1 = 0.5_f64;
        let f2 = 0.5_f64;
        let k_hs_lower = hs_bulk_lower_bound(k1, g1, k2, f1, f2);
        let k_v = voigt_modulus(&[k1, k2], &[f1, f2]);
        assert!(k_hs_lower <= k_v + 1e-6, "HS lower <= Voigt");
    }

    #[test]
    fn test_hs_shear_lower_finite() {
        let k1 = 50e9_f64;
        let g1 = 30e9_f64;
        let g2 = 90e9_f64;
        let f1 = 0.4_f64;
        let f2 = 0.6_f64;
        let g_hs = hs_shear_lower_bound(k1, g1, g2, f1, f2);
        assert!(g_hs.is_finite() && g_hs > 0.0, "g_hs = {g_hs}");
    }

    #[test]
    fn test_hs_shear_upper_ge_lower() {
        let k1 = 50e9_f64;
        let g1 = 30e9_f64;
        let k2 = 200e9_f64;
        let g2 = 90e9_f64;
        let f1 = 0.5_f64;
        let f2 = 0.5_f64;
        let lower = hs_shear_lower_bound(k1, g1, g2, f1, f2);
        let upper = hs_shear_upper_bound(k2, g2, g1, f1, f2);
        assert!(upper >= lower - 1e-3, "HS+ >= HS-: {upper} vs {lower}");
    }

    // ── Fiber composite ROM tests ─────────────────────────────────────────────

    #[test]
    fn test_fiber_composite_rom_e1_rule_of_mixtures() {
        let vf = 0.5_f64;
        let e_f = 230e9_f64;
        let e_m = 3.5e9_f64;
        let c = fiber_composite_rom(e_f, e_m, 0.2, 0.35, 88e9, 1.3e9, vf);
        let e1_expected = vf * e_f + (1.0 - vf) * e_m;
        assert!(
            (c.e1 - e1_expected).abs() / e1_expected < 1e-10,
            "E1 rule of mixtures: {} vs {}",
            c.e1,
            e1_expected
        );
    }

    #[test]
    fn test_fiber_composite_rom_e2_finite() {
        let c = fiber_composite_rom(230e9, 3.5e9, 0.2, 0.35, 88e9, 1.3e9, 0.6);
        assert!(c.e2.is_finite() && c.e2 > 0.0, "E2 = {}", c.e2);
    }

    #[test]
    fn test_fiber_composite_rom_e1_ge_e2() {
        let c = fiber_composite_rom(230e9, 3.5e9, 0.2, 0.35, 88e9, 1.3e9, 0.6);
        assert!(c.e1 >= c.e2, "E1={} should be >= E2={}", c.e1, c.e2);
    }

    #[test]
    fn test_critical_fiber_volume_fraction_positive() {
        let v_crit = critical_fiber_volume_fraction(50e6, 35e6, 3500e6);
        assert!(v_crit >= 0.0 && v_crit.is_finite(), "v_crit = {v_crit}");
    }

    // ── Rule of mixtures tests ────────────────────────────────────────────────

    #[test]
    fn test_rule_of_mixtures_density_two_phase() {
        let rho = rule_of_mixtures_density(&[1000.0, 7800.0], &[0.8, 0.2]);
        let expected = 0.8 * 1000.0 + 0.2 * 7800.0;
        assert!((rho - expected).abs() < 1e-8, "density = {rho}");
    }

    #[test]
    fn test_rule_of_mixtures_cte_weighted() {
        let alphas = vec![12e-6_f64, 0.5e-6];
        let moduli = vec![70e9_f64, 400e9];
        let fracs = vec![0.5_f64, 0.5];
        let cte = rule_of_mixtures_cte(&alphas, &moduli, &fracs);
        assert!(cte.is_finite() && cte > 0.0, "CTE = {cte}");
    }

    #[test]
    fn test_inverse_rule_of_mixtures_equal_phases() {
        // Reuss: 1/E = f1/E1 + f2/E2
        let e = inverse_rule_of_mixtures(&[100.0, 100.0], &[0.5, 0.5]);
        assert!((e - 100.0).abs() < 1e-8, "uniform phases: e = {e}");
    }

    // ── CLT tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_lamina_reduced_stiffness_positive_diagonal() {
        let l = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let q = lamina_reduced_stiffness(&l);
        assert!(q[0][0] > 0.0, "Q11 > 0: {}", q[0][0]);
        assert!(q[1][1] > 0.0, "Q22 > 0: {}", q[1][1]);
        assert!(q[2][2] > 0.0, "Q66 > 0: {}", q[2][2]);
    }

    #[test]
    fn test_lamina_reduced_stiffness_symmetric() {
        let l = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let q = lamina_reduced_stiffness(&l);
        assert!((q[0][1] - q[1][0]).abs() < 1e-3, "Q should be symmetric");
    }

    #[test]
    fn test_transform_reduced_stiffness_zero_angle() {
        let l = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let q = lamina_reduced_stiffness(&l);
        let q_bar = transform_reduced_stiffness(&q, 0.0);
        assert!((q_bar[0][0] - q[0][0]).abs() < 1e-3, "0° transform: Q11");
        assert!((q_bar[1][1] - q[1][1]).abs() < 1e-3, "0° transform: Q22");
    }

    #[test]
    fn test_clt_single_ply_a_matrix_finite() {
        let ply = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let res = clt_analysis(&[ply]);
        for row in &res.a_mat {
            for &v in row {
                assert!(v.is_finite(), "A matrix must be finite: {v}");
            }
        }
    }

    #[test]
    fn test_clt_symmetric_laminate_b_near_zero() {
        // [0/90/90/0] symmetric laminate
        let ply0 = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let ply90 = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: std::f64::consts::PI / 2.0,
        };
        let layers = vec![
            LaminaLayer {
                e1: 140e9,
                e2: 10e9,
                g12: 5e9,
                nu12: 0.3,
                thickness: 0.001,
                theta: 0.0,
            },
            LaminaLayer {
                e1: 140e9,
                e2: 10e9,
                g12: 5e9,
                nu12: 0.3,
                thickness: 0.001,
                theta: std::f64::consts::PI / 2.0,
            },
            LaminaLayer {
                e1: 140e9,
                e2: 10e9,
                g12: 5e9,
                nu12: 0.3,
                thickness: 0.001,
                theta: std::f64::consts::PI / 2.0,
            },
            LaminaLayer {
                e1: 140e9,
                e2: 10e9,
                g12: 5e9,
                nu12: 0.3,
                thickness: 0.001,
                theta: 0.0,
            },
        ];
        let _ = (ply0, ply90); // suppress unused warning
        let res = clt_analysis(&layers);
        assert!(
            laminate_is_symmetric(&res.b_mat, 1.0),
            "B should be ~0 for symmetric laminate"
        );
    }

    #[test]
    fn test_clt_total_thickness() {
        let layers: Vec<LaminaLayer> = (0..5)
            .map(|_| LaminaLayer {
                e1: 100e9,
                e2: 8e9,
                g12: 4e9,
                nu12: 0.3,
                thickness: 0.002,
                theta: 0.0,
            })
            .collect();
        let res = clt_analysis(&layers);
        assert!(
            (res.total_thickness - 0.01).abs() < 1e-12,
            "total h = {}",
            res.total_thickness
        );
    }

    #[test]
    fn test_clt_effective_ex_positive() {
        let ply = LaminaLayer {
            e1: 140e9,
            e2: 10e9,
            g12: 5e9,
            nu12: 0.3,
            thickness: 0.001,
            theta: 0.0,
        };
        let res = clt_analysis(&[ply]);
        assert!(res.ex > 0.0, "Ex should be positive: {}", res.ex);
    }

    // ── Composite moduli tests ────────────────────────────────────────────────

    #[test]
    fn test_walpole_bulk_bounds_ordered() {
        let k1 = 50e9_f64;
        let g1 = 30e9_f64;
        let k2 = 200e9_f64;
        let g2 = 80e9_f64;
        let f2 = 0.4_f64;
        let (lo, hi) = walpole_bulk_bounds(k1, g1, k2, g2, f2);
        assert!(hi >= lo, "HS+ >= HS-: {hi} vs {lo}");
    }

    #[test]
    fn test_self_consistent_bulk_between_phases() {
        let k1 = 50e9_f64;
        let g1 = 30e9_f64;
        let k2 = 200e9_f64;
        let g2 = 80e9_f64;
        let k_sc = self_consistent_bulk(k1, g1, k2, g2, 0.5, 0.5);
        assert!(
            k_sc >= k1 && k_sc <= k2,
            "SC should be between phases: {k_sc}"
        );
    }

    #[test]
    fn test_differential_scheme_bulk_pure_matrix() {
        let k_m = 100e9_f64;
        let g_m = 40e9_f64;
        let k_ds = differential_scheme_bulk(k_m, g_m, 300e9, 0.0, 10);
        assert!((k_ds - k_m).abs() < 1e-6, "zero inclusion: k = {k_ds}");
    }

    #[test]
    fn test_mori_tanaka_young_between_phases() {
        let k_m = 50e9_f64;
        let g_m = 30e9_f64;
        let k_i = 300e9_f64;
        let g_i = 150e9_f64;
        let e_m = 9.0 * k_m * g_m / (3.0 * k_m + g_m);
        let e_i = 9.0 * k_i * g_i / (3.0 * k_i + g_i);
        let e_mt = mori_tanaka_young(k_m, g_m, k_i, g_i, 0.3);
        assert!(
            e_mt >= e_m && e_mt <= e_i,
            "MT E = {e_mt} vs [{e_m}, {e_i}]"
        );
    }

    #[test]
    fn test_mori_tanaka_poisson_bounded() {
        let nu = mori_tanaka_poisson(50e9, 30e9, 300e9, 150e9, 0.3);
        assert!(
            nu > 0.0 && nu < 0.5,
            "Poisson ratio must be in (0, 0.5): {nu}"
        );
    }

    // ── RVE / Hill-Mandel tests ───────────────────────────────────────────────

    #[test]
    fn test_rve_average_strain_single_phase() {
        let s = rve_average_strain(&[[1.0, 2.0, 0.5]], &[1.0]);
        assert!((s[0] - 1.0).abs() < 1e-10 && (s[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_rve_average_stress_two_phases() {
        let st = rve_average_stress(&[[100.0, 0.0, 0.0], [200.0, 0.0, 0.0]], &[0.5, 0.5]);
        assert!((st[0] - 150.0).abs() < 1e-8, "avg stress = {}", st[0]);
    }

    #[test]
    fn test_hill_mandel_inplane_uniform() {
        // Uniform σ and ε: Hill-Mandel error should be zero
        let stress = vec![[100.0_f64, 50.0, 10.0]; 2];
        let strain = vec![[0.001_f64, 0.0005, 0.0001]; 2];
        let fracs = vec![0.5_f64, 0.5];
        let err = hill_mandel_inplane(&stress, &strain, &fracs);
        assert!(err.abs() < 1e-15, "Hill-Mandel error = {err}");
    }

    // ── Failure criteria tests ────────────────────────────────────────────────

    #[test]
    fn test_tsai_hill_no_load_zero() {
        let fi = tsai_hill_criterion(0.0, 0.0, 0.0, 1500e6, 50e6, 70e6);
        assert_eq!(fi, 0.0);
    }

    #[test]
    fn test_tsai_hill_at_limit_is_one() {
        // σ₁ = X, σ₂ = 0, τ = 0: FI = 1
        let fi = tsai_hill_criterion(1500e6, 0.0, 0.0, 1500e6, 50e6, 70e6);
        assert!((fi - 1.0).abs() < 1e-10, "FI = {fi}");
    }

    #[test]
    fn test_max_stress_failure_no_failure() {
        let fail = max_stress_failure(100e6, 20e6, 10e6, 1500e6, 1200e6, 50e6, 150e6, 70e6);
        assert!(!fail, "Should not fail under low stresses");
    }

    #[test]
    fn test_max_stress_failure_triggers() {
        let fail = max_stress_failure(2000e6, 0.0, 0.0, 1500e6, 1200e6, 50e6, 150e6, 70e6);
        assert!(fail, "σ₁ > Xt should trigger failure");
    }

    #[test]
    fn test_thermal_residual_stress_sign() {
        // α_m > α_f and ΔT < 0 (cool-down): residual tension in fiber
        let sigma = thermal_residual_stress(70e9, 23e-6, 0.5e-6, -100.0);
        assert!(
            sigma < 0.0,
            "residual stress should be compressive: {sigma}"
        );
    }

    #[test]
    fn test_puck_inter_fiber_zero_stress() {
        let fi = puck_inter_fiber_failure(0.0, 0.0, 50e6, 70e6, 0.25);
        assert_eq!(fi, 0.0);
    }
}
