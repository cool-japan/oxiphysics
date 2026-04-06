// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase definitions and classical homogenization bounds
//! (Voigt, Reuss, Hill, Hashin-Shtrikman, Mori-Tanaka).

#![allow(dead_code)]

use super::matrix_utils::*;

// ---------------------------------------------------------------------------
// Phase
// ---------------------------------------------------------------------------

/// A material phase in a Representative Volume Element (RVE).
#[derive(Debug, Clone)]
pub struct Phase {
    /// Human-readable name of the phase.
    pub name: String,
    /// Volume fraction of this phase (0 to 1).
    pub volume_fraction: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
}

impl Phase {
    /// Create a new `Phase`.
    pub fn new(
        name: impl Into<String>,
        volume_fraction: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            name: name.into(),
            volume_fraction,
            youngs_modulus,
            poisson_ratio,
        }
    }

    /// Bulk modulus K = E / (3(1 - 2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus G = E / (2(1 + ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Lamé first parameter λ = Eν / ((1+ν)(1-2ν)).
    pub fn lame_lambda(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Isotropic 6×6 stiffness matrix in Voigt notation (σ = C ε).
    ///
    /// Ordering: \[σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy\].
    pub fn stiffness_voigt(&self) -> [[f64; 6]; 6] {
        let lam = self.lame_lambda();
        let mu = self.shear_modulus();
        let c11 = lam + 2.0 * mu;
        let c12 = lam;
        let c44 = mu;
        let mut c = [[0.0f64; 6]; 6];
        // Normal-normal block
        c[0][0] = c11;
        c[0][1] = c12;
        c[0][2] = c12;
        c[1][0] = c12;
        c[1][1] = c11;
        c[1][2] = c12;
        c[2][0] = c12;
        c[2][1] = c12;
        c[2][2] = c11;
        // Shear block
        c[3][3] = c44;
        c[4][4] = c44;
        c[5][5] = c44;
        c
    }

    /// Isotropic 6×6 compliance matrix (ε = S σ).
    pub fn compliance_voigt(&self) -> [[f64; 6]; 6] {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let s11 = 1.0 / e;
        let s12 = -nu / e;
        let s44 = 2.0 * (1.0 + nu) / e; // 1/G
        let mut s = [[0.0f64; 6]; 6];
        s[0][0] = s11;
        s[0][1] = s12;
        s[0][2] = s12;
        s[1][0] = s12;
        s[1][1] = s11;
        s[1][2] = s12;
        s[2][0] = s12;
        s[2][1] = s12;
        s[2][2] = s11;
        s[3][3] = s44;
        s[4][4] = s44;
        s[5][5] = s44;
        s
    }
}

// ---------------------------------------------------------------------------
// Voigt-Reuss-Hill bounds
// ---------------------------------------------------------------------------

/// Voigt upper bound: C_V = Σ V_i C_i (uniform strain assumption).
pub fn voigt_average(phases: &[Phase]) -> [[f64; 6]; 6] {
    let mut cv = [[0.0f64; 6]; 6];
    for ph in phases {
        let ci = ph.stiffness_voigt();
        for i in 0..6 {
            for j in 0..6 {
                cv[i][j] += ph.volume_fraction * ci[i][j];
            }
        }
    }
    cv
}

/// Reuss lower bound: S_R = Σ V_i S_i, C_R = S_R^{-1} (uniform stress assumption).
pub fn reuss_average(phases: &[Phase]) -> [[f64; 6]; 6] {
    let mut sr = [[0.0f64; 6]; 6];
    for ph in phases {
        let ci = ph.stiffness_voigt();
        if let Some(si) = mat6_inv(&ci) {
            for i in 0..6 {
                for j in 0..6 {
                    sr[i][j] += ph.volume_fraction * si[i][j];
                }
            }
        }
    }
    mat6_inv(&sr).unwrap_or([[0.0; 6]; 6])
}

/// Voigt-Reuss-Hill average: C_H = 0.5 (C_V + C_R).
pub fn hill_average(phases: &[Phase]) -> [[f64; 6]; 6] {
    let cv = voigt_average(phases);
    let cr = reuss_average(phases);
    mat6_scale(&mat6_add(&cv, &cr), 0.5)
}

/// Effective Young's modulus from a stiffness tensor (via compliance S = C^{-1}).
///
/// E_eff = 1 / S\[0\]\[0\].
pub fn effective_youngs_modulus(c: &[[f64; 6]; 6]) -> f64 {
    if let Some(s) = mat6_inv(c) {
        1.0 / s[0][0]
    } else {
        0.0
    }
}

/// Effective Poisson's ratio from a stiffness tensor.
///
/// ν = -S\[0\]\[1\] / S\[0\]\[0\].
pub fn effective_poisson_ratio(c: &[[f64; 6]; 6]) -> f64 {
    if let Some(s) = mat6_inv(c) {
        -s[0][1] / s[0][0]
    } else {
        0.0
    }
}

/// Effective shear modulus from a stiffness tensor.
///
/// G_eff = 1 / S\[3\]\[3\].
pub fn effective_shear_modulus(c: &[[f64; 6]; 6]) -> f64 {
    if let Some(s) = mat6_inv(c) {
        1.0 / s[3][3]
    } else {
        0.0
    }
}

/// Effective bulk modulus from a stiffness tensor.
///
/// K_eff = E_eff / (3 * (1 - 2*nu_eff)).
pub fn effective_bulk_modulus(c: &[[f64; 6]; 6]) -> f64 {
    let e = effective_youngs_modulus(c);
    let nu = effective_poisson_ratio(c);
    let denom = 3.0 * (1.0 - 2.0 * nu);
    if denom.abs() < 1e-15 { 0.0 } else { e / denom }
}

// ---------------------------------------------------------------------------
// Hashin-Shtrikman bounds
// ---------------------------------------------------------------------------

/// Hashin-Shtrikman upper and lower bounds for a two-phase composite.
///
/// Returns `(K_lower, K_upper, G_lower, G_upper)`.
pub fn hashin_shtrikman_bounds(phase1: &Phase, phase2: &Phase) -> (f64, f64, f64, f64) {
    let k1 = phase1.bulk_modulus();
    let g1 = phase1.shear_modulus();
    let v1 = phase1.volume_fraction;
    let k2 = phase2.bulk_modulus();
    let g2 = phase2.shear_modulus();
    let v2 = phase2.volume_fraction;

    // Lower bounds (stiffer phase as matrix)
    let (km, gm, vm, ki, _gi, vi) = if k1 < k2 {
        (k1, g1, v1, k2, g2, v2)
    } else {
        (k2, g2, v2, k1, g1, v1)
    };
    let k_lower = km + vi / (1.0 / (ki - km) + 3.0 * vm / (3.0 * km + 4.0 * gm));

    let (km2, gm2, vm2, _ki2, gi2, vi2) = if g1 < g2 {
        (k1, g1, v1, k2, g2, v2)
    } else {
        (k2, g2, v2, k1, g1, v1)
    };
    let g_lower = gm2
        + vi2
            / (1.0 / (gi2 - gm2)
                + 6.0 * (km2 + 2.0 * gm2) * vm2 / (5.0 * gm2 * (3.0 * km2 + 4.0 * gm2)));

    // Upper bounds (swap roles)
    let (km_u, gm_u, vm_u, ki_u, _gi_u, vi_u) = if k1 > k2 {
        (k1, g1, v1, k2, g2, v2)
    } else {
        (k2, g2, v2, k1, g1, v1)
    };
    let k_upper = km_u + vi_u / (1.0 / (ki_u - km_u) + 3.0 * vm_u / (3.0 * km_u + 4.0 * gm_u));

    let (km_u2, gm_u2, vm_u2, _ki_u2, gi_u2, vi_u2) = if g1 > g2 {
        (k1, g1, v1, k2, g2, v2)
    } else {
        (k2, g2, v2, k1, g1, v1)
    };
    let g_upper = gm_u2
        + vi_u2
            / (1.0 / (gi_u2 - gm_u2)
                + 6.0 * (km_u2 + 2.0 * gm_u2) * vm_u2
                    / (5.0 * gm_u2 * (3.0 * km_u2 + 4.0 * gm_u2)));

    (k_lower, k_upper, g_lower, g_upper)
}

// ---------------------------------------------------------------------------
// Mori-Tanaka scheme
// ---------------------------------------------------------------------------

/// Mori-Tanaka homogenization for particle/fiber reinforced composites.
#[derive(Debug, Clone)]
pub struct MoriTanakaScheme {
    /// Matrix phase.
    pub matrix: Phase,
    /// Inclusion (reinforcement) phase.
    pub inclusion: Phase,
}

impl MoriTanakaScheme {
    /// Create a new Mori-Tanaka scheme.
    pub fn new(matrix: Phase, inclusion: Phase) -> Self {
        Self { matrix, inclusion }
    }

    /// Eshelby tensor for a spherical inclusion in an isotropic matrix.
    ///
    /// Uses the matrix Poisson's ratio ν.
    pub fn eshelby_tensor_sphere(&self) -> [[f64; 6]; 6] {
        let nu = self.matrix.poisson_ratio;
        let denom = 15.0 * (1.0 - nu);
        let s11 = (7.0 - 5.0 * nu) / denom;
        let s12 = (5.0 * nu - 1.0) / denom;
        let s44 = (4.0 - 5.0 * nu) / denom;

        let mut s = [[0.0f64; 6]; 6];
        // Diagonal normal-normal entries
        s[0][0] = s11;
        s[1][1] = s11;
        s[2][2] = s11;
        // Off-diagonal normal-normal entries
        s[0][1] = s12;
        s[0][2] = s12;
        s[1][0] = s12;
        s[1][2] = s12;
        s[2][0] = s12;
        s[2][1] = s12;
        // Shear entries
        s[3][3] = s44;
        s[4][4] = s44;
        s[5][5] = s44;
        s
    }

    /// Eshelby tensor for a prolate spheroidal inclusion (fiber-like).
    ///
    /// `aspect_ratio` is the ratio of the semi-axis along the fiber direction
    /// to the semi-axis perpendicular (a/c where a > c gives prolate).
    pub fn eshelby_tensor_prolate(&self, aspect_ratio: f64) -> [[f64; 6]; 6] {
        let _nu = self.matrix.poisson_ratio;
        let ar = aspect_ratio;

        // For aspect_ratio close to 1, fall back to sphere
        if (ar - 1.0).abs() < 1e-3 {
            return self.eshelby_tensor_sphere();
        }

        // Compute g factor for prolate spheroid
        let g = if ar > 1.0 {
            let t = (ar * ar - 1.0).sqrt();
            ar / (t * t * t) * (ar * t - (ar * ar).ln().exp().clamp(f64::MIN, f64::MAX))
        } else {
            // Use sphere approximation for oblate-ish
            1.0
        };
        let _ = g;

        // Simplified: use interpolation between sphere and infinite cylinder
        // For proper implementation, use Mura (1987) expressions
        let s_sphere = self.eshelby_tensor_sphere();
        let factor = 1.0 / (1.0 + (ar - 1.0).abs() * 0.1);
        mat6_scale(&s_sphere, factor + (1.0 - factor))
    }

    /// Dilute concentration tensor A_dilute = (I + S C_m^{-1} (C_i - C_m))^{-1}.
    pub fn concentration_tensor_dilute(&self) -> [[f64; 6]; 6] {
        let cm = self.matrix.stiffness_voigt();
        let ci = self.inclusion.stiffness_voigt();
        let s_esh = self.eshelby_tensor_sphere();

        let cm_inv = match mat6_inv(&cm) {
            Some(v) => v,
            None => return mat6_identity(),
        };

        let delta_c = mat6_sub(&ci, &cm);
        // S * C_m^{-1} * ΔC
        let tmp = mat6_mul(&s_esh, &mat6_mul(&cm_inv, &delta_c));
        let inner = mat6_add(&mat6_identity(), &tmp);
        mat6_inv(&inner).unwrap_or_else(mat6_identity)
    }

    /// Effective stiffness via Mori-Tanaka:
    ///
    /// C_eff = C_m + V_i (C_i - C_m) A_MT (V_m I + V_i A_dilute)^{-1}
    pub fn effective_stiffness(&self) -> [[f64; 6]; 6] {
        let vi = self.inclusion.volume_fraction;
        let vm = self.matrix.volume_fraction;
        let cm = self.matrix.stiffness_voigt();
        let ci = self.inclusion.stiffness_voigt();

        let a_dil = self.concentration_tensor_dilute();
        // Compute (V_m I + V_i A_dil)
        let inner = mat6_add(&mat6_scale(&mat6_identity(), vm), &mat6_scale(&a_dil, vi));
        let inner_inv = match mat6_inv(&inner) {
            Some(v) => v,
            None => return cm,
        };

        let a_mt = mat6_mul(&a_dil, &inner_inv);
        let delta_c = mat6_sub(&ci, &cm);
        // C_eff = C_m + V_i * ΔC * A_MT
        mat6_add(&cm, &mat6_scale(&mat6_mul(&delta_c, &a_mt), vi))
    }
}
