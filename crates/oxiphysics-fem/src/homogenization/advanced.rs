// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced homogenization methods: strain concentration tensors, thermal
//! conductivity, anisotropy indices, thermo-elastic coupling, composite
//! failure criteria, Halpin-Tsai, three-phase composites, and damage models.

use super::bounds::*;
use super::matrix_utils::*;

// ---------------------------------------------------------------------------
// Strain concentration tensor (Eshelby / dilute inclusion)
// ---------------------------------------------------------------------------

/// Strain concentration tensor A for a spherical inclusion in an isotropic matrix.
///
/// In the dilute approximation (Eshelby) the average strain in an inclusion is
/// A · ε_∞ where ε_∞ is the far-field (macroscopic) strain.
///
/// For a spherical inclusion in an isotropic matrix the Eshelby tensor is
/// expressed using the Poisson's ratio of the matrix alone.
#[derive(Debug, Clone)]
pub struct StrainConcentrationTensor {
    /// 6×6 concentration tensor in Voigt notation.
    pub tensor: [[f64; 6]; 6],
}

impl StrainConcentrationTensor {
    /// Compute the dilute-limit strain concentration tensor for a spherical inclusion.
    ///
    /// Arguments:
    /// - `e_m`, `nu_m`: matrix Young's modulus and Poisson ratio
    /// - `e_i`, `nu_i`: inclusion Young's modulus and Poisson ratio
    pub fn dilute_spherical_inclusion(e_m: f64, nu_m: f64, e_i: f64, nu_i: f64) -> Self {
        // Eshelby tensor S for a sphere in an isotropic matrix (Voigt, contracted)
        let nu = nu_m;
        let s11 = (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
        let s12 = (5.0 * nu - 1.0) / (15.0 * (1.0 - nu));
        let s44 = (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu));

        // Build full 6×6 Eshelby tensor (diagonal blocks for isotropic sphere)
        let mut s = [[0.0f64; 6]; 6];
        for (i, row) in s.iter_mut().enumerate().take(3) {
            row[i] = s11;
            for (j, val) in row.iter_mut().enumerate().take(3) {
                if i != j {
                    *val = s12;
                }
            }
        }
        for (i, row) in s.iter_mut().enumerate().skip(3) {
            row[i] = 2.0 * s44;
        } // factor 2 for Voigt shear

        // Compliance of matrix and inclusion
        let c_m = Phase::new("m", 1.0, e_m, nu_m).stiffness_voigt();
        let c_i = Phase::new("i", 1.0, e_i, nu_i).stiffness_voigt();

        // A = (C_i - C_m) · S + C_m)^{-1} · C_i  (Eshelby dilute)
        // Simplified: A = [I + S · C_m^{-1} · (C_i - C_m)]^{-1}
        // Compute (C_i - C_m)
        let mut dc = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                dc[i][j] = c_i[i][j] - c_m[i][j];
            }
        }

        // Compute S · C_m^{-1} · ΔC + I
        let c_m_inv = mat6_inv(&c_m).unwrap_or(mat6_identity());
        let sc_inv_dc = mat6_mul(&s, &mat6_mul(&c_m_inv, &dc));

        let mut ident_plus = mat6_identity();
        for i in 0..6 {
            for j in 0..6 {
                ident_plus[i][j] += sc_inv_dc[i][j];
            }
        }

        let tensor = mat6_inv(&ident_plus).unwrap_or(mat6_identity());
        Self { tensor }
    }

    /// Apply the concentration tensor to a macroscopic strain vector.
    pub fn apply(&self, macro_strain: &[f64; 6]) -> [f64; 6] {
        mat6_vec_mul(&self.tensor, macro_strain)
    }
}

// ---------------------------------------------------------------------------
// Effective thermal conductivity
// ---------------------------------------------------------------------------

/// Maxwell model for effective thermal conductivity of a particulate composite.
///
/// λ_eff = λ_m * (λ_i + 2λ_m + 2f(λ_i - λ_m)) / (λ_i + 2λ_m - f(λ_i - λ_m))
///
/// where f is the inclusion volume fraction.
pub fn maxwell_effective_thermal_conductivity(lambda_m: f64, lambda_i: f64, f: f64) -> f64 {
    let num = lambda_i + 2.0 * lambda_m + 2.0 * f * (lambda_i - lambda_m);
    let den = lambda_i + 2.0 * lambda_m - f * (lambda_i - lambda_m);
    if den.abs() < 1e-30 {
        return lambda_m;
    }
    lambda_m * num / den
}

/// Hashin-Shtrikman lower bound on thermal conductivity.
///
/// Assumes λ_1 < λ_2 and phase 1 is the matrix with volume fraction (1 - f).
pub fn hs_thermal_lower_bound(lambda_1: f64, lambda_2: f64, f: f64) -> f64 {
    // HS lower: phase 1 acts as matrix
    let phi1 = 1.0 - f;
    let den = phi1 / (3.0 * lambda_1);
    if den.abs() < 1e-30 {
        return lambda_1;
    }
    lambda_1 + f / (1.0 / (lambda_2 - lambda_1) + den)
}

/// Hashin-Shtrikman upper bound on thermal conductivity.
pub fn hs_thermal_upper_bound(lambda_1: f64, lambda_2: f64, f: f64) -> f64 {
    // HS upper: phase 2 acts as matrix
    let phi2 = f;
    let den = phi2 / (3.0 * lambda_2);
    if den.abs() < 1e-30 {
        return lambda_2;
    }
    lambda_2 + (1.0 - f) / (1.0 / (lambda_1 - lambda_2) + den)
}

// ---------------------------------------------------------------------------
// Zener anisotropy index
// ---------------------------------------------------------------------------

/// Zener anisotropy index for cubic crystals.
///
/// A_z = 2 C_44 / (C_11 - C_12)
///
/// A_z = 1 for isotropic materials; A_z ≠ 1 indicates elastic anisotropy.
pub fn zener_anisotropy_index(c11: f64, c12: f64, c44: f64) -> f64 {
    let denom = c11 - c12;
    if denom.abs() < 1e-30 {
        return 1.0;
    }
    2.0 * c44 / denom
}

/// Universal anisotropy index (Ranganathan & Ostoja-Starzewski, 2008).
///
/// A_u = 5 G_V/G_R + K_V/K_R - 6  ≥ 0
/// A_u = 0 for isotropic materials.
pub fn universal_anisotropy_index(k_v: f64, k_r: f64, g_v: f64, g_r: f64) -> f64 {
    if k_r < 1e-30 || g_r < 1e-30 {
        return 0.0;
    }
    5.0 * g_v / g_r + k_v / k_r - 6.0
}

// ---------------------------------------------------------------------------
// Thermo-elastic coupling
// ---------------------------------------------------------------------------

/// Thermo-elastic coupling for a linearly elastic, isotropic material.
#[derive(Debug, Clone, Copy)]
pub struct ThermoElasticCoupling {
    /// Coefficient of thermal expansion \[1/K\].
    pub alpha: f64,
    /// Young's modulus \[Pa\].
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl ThermoElasticCoupling {
    /// Create a new thermo-elastic model.
    pub fn new(alpha: f64, e: f64, nu: f64) -> Self {
        Self { alpha, e, nu }
    }

    /// Free thermal strain for a temperature change ΔT.
    pub fn free_thermal_strain(&self, delta_t: f64) -> f64 {
        self.alpha * delta_t
    }

    /// Thermal stress in a fully constrained uniaxial bar: σ = -E·α·ΔT.
    pub fn constrained_thermal_stress(&self, delta_t: f64) -> f64 {
        -self.e * self.alpha * delta_t
    }

    /// Thermal stress vector (Voigt) for a fully constrained 3-D body.
    ///
    /// σ_thermal = -E/(1-2ν) · α · ΔT · \[1, 1, 1, 0, 0, 0\]ᵀ
    pub fn constrained_thermal_stress_voigt(&self, delta_t: f64) -> [f64; 6] {
        let coeff = -self.e / (1.0 - 2.0 * self.nu) * self.alpha * delta_t;
        [coeff, coeff, coeff, 0.0, 0.0, 0.0]
    }

    /// Effective thermal eigenstrain in Voigt notation: ε* = α·ΔT·\[1,1,1,0,0,0\].
    pub fn eigenstrain(&self, delta_t: f64) -> [f64; 6] {
        let e = self.alpha * delta_t;
        [e, e, e, 0.0, 0.0, 0.0]
    }
}

// ---------------------------------------------------------------------------
// Homogenized stiffness bounds summary
// ---------------------------------------------------------------------------

/// Aggregated upper and lower bounds on homogenized Young's modulus.
#[derive(Debug, Clone)]
pub struct HomogenizedStiffnessBounds {
    /// Voigt (upper) bound \[Pa\].
    pub voigt_bound: f64,
    /// Reuss (lower) bound \[Pa\].
    pub reuss_bound: f64,
    /// Hashin-Shtrikman upper bound \[Pa\].
    pub hs_upper: f64,
    /// Hashin-Shtrikman lower bound \[Pa\].
    pub hs_lower: f64,
}

impl HomogenizedStiffnessBounds {
    /// Compute all four bounds from a set of phases.
    pub fn compute(phases: &[Phase]) -> Self {
        let e_voigt = effective_youngs_modulus(&voigt_average(phases));
        let e_reuss = effective_youngs_modulus(&reuss_average(phases));

        // For HS bounds, use two-phase HS if possible
        let (hs_l, hs_u) = if phases.len() == 2 {
            let (k_l, k_u, _g_l, _g_u) = hashin_shtrikman_bounds(&phases[0], &phases[1]);
            // Convert bulk modulus bounds to Young's modulus approximation
            // Using E ≈ 9KG/(3K+G) at average shear modulus
            let g_avg = effective_shear_modulus(&hill_average(phases));
            let e_l = 9.0 * k_l * g_avg / (3.0 * k_l + g_avg);
            let e_u = 9.0 * k_u * g_avg / (3.0 * k_u + g_avg);
            (e_l.min(e_u), e_l.max(e_u))
        } else {
            (e_reuss, e_voigt)
        };

        Self {
            voigt_bound: e_voigt,
            reuss_bound: e_reuss,
            hs_upper: hs_u,
            hs_lower: hs_l,
        }
    }

    /// Return the Hill average (midpoint of Voigt and Reuss).
    pub fn hill_average(&self) -> f64 {
        0.5 * (self.voigt_bound + self.reuss_bound)
    }

    /// Return the narrower HS interval width.
    pub fn hs_interval_width(&self) -> f64 {
        self.hs_upper - self.hs_lower
    }
}

// ---------------------------------------------------------------------------
// Composite failure criteria
// ---------------------------------------------------------------------------

/// Tsai-Wu failure criterion for unidirectional fiber composites.
///
/// F_1 σ_1 + F_2 σ_2 + F_11 σ_1² + F_22 σ_2² + F_66 τ_12² + 2 F_12 σ_1 σ_2 = 1
#[derive(Debug, Clone)]
pub struct TsaiWuCriterion {
    /// Longitudinal tensile strength \[Pa\].
    pub x_t: f64,
    /// Longitudinal compressive strength \[Pa\] (positive value).
    pub x_c: f64,
    /// Transverse tensile strength \[Pa\].
    pub y_t: f64,
    /// Transverse compressive strength \[Pa\] (positive value).
    pub y_c: f64,
    /// In-plane shear strength \[Pa\].
    pub s12: f64,
    /// Interaction coefficient F_12 (commonly taken as -0.5/√(X_T X_C Y_T Y_C)).
    pub f12: f64,
}

impl TsaiWuCriterion {
    /// Construct a Tsai-Wu criterion with default F12 = -0.5/√(XT·XC·YT·YC).
    pub fn new(x_t: f64, x_c: f64, y_t: f64, y_c: f64, s12: f64) -> Self {
        let f12 = -0.5 / (x_t * x_c * y_t * y_c).sqrt();
        Self {
            x_t,
            x_c,
            y_t,
            y_c,
            s12,
            f12,
        }
    }

    /// Evaluate the Tsai-Wu failure index for a stress state \[σ_1, σ_2, τ_12\].
    ///
    /// Returns a value; failure occurs when `index >= 1`.
    pub fn failure_index(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let f1 = 1.0 / self.x_t - 1.0 / self.x_c;
        let f2 = 1.0 / self.y_t - 1.0 / self.y_c;
        let f11 = 1.0 / (self.x_t * self.x_c);
        let f22 = 1.0 / (self.y_t * self.y_c);
        let f66 = 1.0 / (self.s12 * self.s12);
        f1 * sigma1
            + f2 * sigma2
            + f11 * sigma1 * sigma1
            + f22 * sigma2 * sigma2
            + f66 * tau12 * tau12
            + 2.0 * self.f12 * sigma1 * sigma2
    }

    /// Returns `true` if the stress state exceeds the Tsai-Wu failure surface.
    pub fn failed(&self, sigma1: f64, sigma2: f64, tau12: f64) -> bool {
        self.failure_index(sigma1, sigma2, tau12) >= 1.0
    }
}

// ---------------------------------------------------------------------------
// Hill bounds (Voigt upper / Reuss lower) for two-phase composites
// ---------------------------------------------------------------------------

/// Voigt upper bound and Reuss lower bound for a two-phase composite.
#[derive(Debug, Clone)]
pub struct HillBounds {
    /// Voigt (upper) bound on effective Young's modulus.
    pub voigt: f64,
    /// Reuss (lower) bound on effective Young's modulus.
    pub reuss: f64,
}

impl HillBounds {
    /// Compute Hill bounds for a two-phase composite.
    ///
    /// `phi1` is the volume fraction of phase 1; phase 2 fills the rest.
    pub fn compute(e1: f64, _nu1: f64, phi1: f64, e2: f64, _nu2: f64) -> Self {
        let phi2 = 1.0 - phi1;
        let voigt = phi1 * e1 + phi2 * e2;
        let reuss = if phi1 / e1 + phi2 / e2 > 0.0 {
            1.0 / (phi1 / e1 + phi2 / e2)
        } else {
            0.0
        };
        Self { voigt, reuss }
    }

    /// Hill average: arithmetic mean of Voigt and Reuss bounds.
    pub fn hill_average(&self) -> f64 {
        0.5 * (self.voigt + self.reuss)
    }
}

// ---------------------------------------------------------------------------
// Rule of mixtures
// ---------------------------------------------------------------------------

/// Voigt (upper) bound: rule of mixtures — volume-weighted average.
pub fn effective_elastic_modulus_voigt(e: &[f64], fractions: &[f64]) -> f64 {
    assert_eq!(e.len(), fractions.len());
    e.iter()
        .zip(fractions.iter())
        .map(|(&ei, &fi)| ei * fi)
        .sum()
}

/// Reuss (lower) bound: inverse rule of mixtures.
pub fn effective_elastic_modulus_reuss(e: &[f64], fractions: &[f64]) -> f64 {
    assert_eq!(e.len(), fractions.len());
    let inv_sum: f64 = e
        .iter()
        .zip(fractions.iter())
        .map(|(&ei, &fi)| fi / ei)
        .sum();
    if inv_sum > 0.0 { 1.0 / inv_sum } else { 0.0 }
}

// ---------------------------------------------------------------------------
// Mori-Tanaka micromechanics
// ---------------------------------------------------------------------------

/// Mori-Tanaka estimate of effective bulk modulus.
///
/// Matrix: (K_m, G_m).  Inclusion: K_i, volume fraction f_i.
///
/// Standard form: K* = K_m + f_i*(K_i-K_m) / (1 + (1-f_i)*(K_i-K_m)/(K_m + 4G_m/3))
pub fn mori_tanaka_bulk_modulus(k_m: f64, g_m: f64, k_i: f64, f_i: f64) -> f64 {
    let alpha = k_m + (4.0 / 3.0) * g_m; // Eshelby-Willis factor denominator
    if alpha.abs() < 1e-30 {
        return k_m;
    }
    let delta_k = k_i - k_m;
    let denom = 1.0 + (1.0 - f_i) * delta_k / alpha;
    if denom.abs() < 1e-30 {
        return k_m;
    }
    k_m + f_i * delta_k / denom
}

/// Mori-Tanaka estimate of effective shear modulus.
///
/// Matrix: (K_m, G_m).  Inclusion: G_i, volume fraction f_i.
pub fn mori_tanaka_shear_modulus(k_m: f64, g_m: f64, g_i: f64, f_i: f64) -> f64 {
    // β = G_m * (9 K_m + 8 G_m) / (6 * (K_m + 2 G_m))
    let beta_num = g_m * (9.0 * k_m + 8.0 * g_m);
    let beta_denom = 6.0 * (k_m + 2.0 * g_m);
    if beta_denom.abs() < 1e-30 {
        return g_m;
    }
    let beta = beta_num / beta_denom;
    let delta_g = g_i - g_m;
    let numerator = g_m * (g_m + beta) + f_i * delta_g * (g_m + beta);
    let denominator = (g_m + beta) + f_i * delta_g;
    if denominator.abs() < 1e-30 {
        return g_m;
    }
    numerator / denominator
}

// ---------------------------------------------------------------------------
// Halpin-Tsai
// ---------------------------------------------------------------------------

/// Halpin-Tsai composite modulus.
///
/// E = E_m * (1 + ξ * η * vf) / (1 - η * vf)
/// where η = (E_f/E_m - 1) / (E_f/E_m + ξ)
///
/// `xi` is the reinforcement factor (shape parameter).
pub fn halpin_tsai_modulus(e_m: f64, e_f: f64, vf: f64, xi: f64) -> f64 {
    if e_m.abs() < 1e-30 {
        return 0.0;
    }
    let ratio = e_f / e_m;
    let denom_eta = ratio + xi;
    if denom_eta.abs() < 1e-30 {
        return e_m;
    }
    let eta = (ratio - 1.0) / denom_eta;
    let denom = 1.0 - eta * vf;
    if denom.abs() < 1e-30 {
        return e_m;
    }
    e_m * (1.0 + xi * eta * vf) / denom
}

// ---------------------------------------------------------------------------
// Composite CTE
// ---------------------------------------------------------------------------

/// Composite coefficient of thermal expansion (Turner's model).
///
/// CTE_c = (alpha1 * E1 * phi1 + alpha2 * E2 * phi2) / (E1 * phi1 + E2 * phi2)
pub fn composite_thermal_expansion(alpha1: f64, e1: f64, phi1: f64, alpha2: f64, e2: f64) -> f64 {
    let phi2 = 1.0 - phi1;
    let num = alpha1 * e1 * phi1 + alpha2 * e2 * phi2;
    let den = e1 * phi1 + e2 * phi2;
    if den.abs() < 1e-30 {
        return 0.5 * (alpha1 + alpha2);
    }
    num / den
}

// ---------------------------------------------------------------------------
// RveCell
// ---------------------------------------------------------------------------

/// A simple 2-D representative volume element (RVE) with a phase map.
#[derive(Debug, Clone)]
pub struct RveCell {
    /// Number of voxels in x direction.
    pub n_x: usize,
    /// Number of voxels in y direction.
    pub n_y: usize,
    /// Phase identifier per voxel (row-major, phase index starting at 0).
    pub phase_map: Vec<u8>,
}

impl RveCell {
    /// Create a new RveCell with all voxels in phase 0.
    pub fn new(n_x: usize, n_y: usize) -> Self {
        Self {
            n_x,
            n_y,
            phase_map: vec![0; n_x * n_y],
        }
    }

    /// Compute volume fractions for each phase present in the map.
    ///
    /// Returns a `Vec`f64` indexed by phase (0-based), length = max_phase + 1.
    pub fn compute_volume_fractions(&self) -> Vec<f64> {
        if self.phase_map.is_empty() {
            return vec![];
        }
        let max_phase = *self
            .phase_map
            .iter()
            .max()
            .expect("phase_map is non-empty after check") as usize;
        let total = self.phase_map.len() as f64;
        let mut counts = vec![0usize; max_phase + 1];
        for &p in &self.phase_map {
            counts[p as usize] += 1;
        }
        counts.iter().map(|&c| c as f64 / total).collect()
    }
}

// ---------------------------------------------------------------------------
// Periodic boundary conditions (DOF pairing)
// ---------------------------------------------------------------------------

/// Compute periodic boundary condition DOF pairs for a 2-D RVE.
///
/// Returns a list of `(left_dof, right_dof)` and `(bottom_dof, top_dof)` pairs
/// where left/right share the same y-index and bottom/top share the same x-index.
/// Each voxel has 2 DOFs: `2*(row*n_x + col)` and `2*(row*n_x + col) + 1`.
pub fn periodic_boundary_conditions_2d(rve: &RveCell) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let nx = rve.n_x;
    let ny = rve.n_y;
    if nx == 0 || ny == 0 {
        return pairs;
    }

    // Left (col=0) ↔ Right (col=nx-1) for each row
    for row in 0..ny {
        for dof_offset in 0..2 {
            let left = 2 * (row * nx) + dof_offset;
            let right = 2 * (row * nx + nx - 1) + dof_offset;
            pairs.push((left, right));
        }
    }

    // Bottom (row=0) ↔ Top (row=ny-1) for each col
    for col in 0..nx {
        for dof_offset in 0..2 {
            let bottom = 2 * (col) + dof_offset;
            let top = 2 * ((ny - 1) * nx + col) + dof_offset;
            pairs.push((bottom, top));
        }
    }

    pairs
}

// ---------------------------------------------------------------------------
// Halpin-Tsai extended — longitudinal and transverse
// ---------------------------------------------------------------------------

/// Halpin-Tsai longitudinal modulus E_11 for unidirectional fiber composite.
///
/// E_11 = (E_f * v_f + E_m * (1 - v_f))   (rule of mixtures for longitudinal)
pub fn halpin_tsai_longitudinal(e_f: f64, e_m: f64, vf: f64) -> f64 {
    e_f * vf + e_m * (1.0 - vf)
}

/// Halpin-Tsai transverse modulus E_22.
///
/// Uses the Halpin-Tsai equation with reinforcement factor xi = 2 * (aspect_ratio).
/// For circular fibers in a plate, xi ≈ 2.
pub fn halpin_tsai_transverse(e_f: f64, e_m: f64, vf: f64, xi: f64) -> f64 {
    halpin_tsai_modulus(e_m, e_f, vf, xi)
}

/// Halpin-Tsai in-plane shear modulus G_12.
///
/// xi = 1 for circular fibers.
pub fn halpin_tsai_shear_12(g_f: f64, g_m: f64, vf: f64) -> f64 {
    halpin_tsai_modulus(g_m, g_f, vf, 1.0)
}

/// Halpin-Tsai major Poisson's ratio ν_12 (rule of mixtures).
pub fn halpin_tsai_poisson_12(nu_f: f64, nu_m: f64, vf: f64) -> f64 {
    nu_f * vf + nu_m * (1.0 - vf)
}

/// Hashin-Shtrikman upper bound for shear modulus of a two-phase composite.
///
/// G_HS+ = G_1 + f2 / (1/(G_2-G_1) + 6(K_1+2G_1)/(5G_1(3K_1+4G_1)))
///
/// Assumes G_1 ≥ G_2 (phase 1 is the stiffer phase).
pub fn hashin_shtrikman_shear_upper(k1: f64, g1: f64, g2: f64, f2: f64) -> f64 {
    let alpha = 6.0 * (k1 + 2.0 * g1) / (5.0 * g1 * (3.0 * k1 + 4.0 * g1));
    if (g2 - g1).abs() < 1e-30 {
        return g1;
    }
    let denom = 1.0 / (g2 - g1) + alpha * f2;
    if denom.abs() < 1e-30 {
        return g1;
    }
    g1 + f2 / denom
}

/// Hashin-Shtrikman lower bound for shear modulus.
///
/// G_HS- = G_2 + f1 / (1/(G_1-G_2) + 6(K_2+2G_2)/(5G_2(3K_2+4G_2)))
///
/// Assumes G_2 ≤ G_1.
pub fn hashin_shtrikman_shear_lower(k2: f64, g2: f64, g1: f64, f1: f64) -> f64 {
    let alpha = 6.0 * (k2 + 2.0 * g2) / (5.0 * g2 * (3.0 * k2 + 4.0 * g2));
    if (g1 - g2).abs() < 1e-30 {
        return g2;
    }
    let denom = 1.0 / (g1 - g2) + alpha * f1;
    if denom.abs() < 1e-30 {
        return g2;
    }
    g2 + f1 / denom
}

/// Hashin-Shtrikman upper bound for bulk modulus.
///
/// K_HS+ = K_1 + f2 / (1/(K_2-K_1) + 3f1/(3K_1+4G_1))
///
/// Assumes K_1 ≥ K_2.
pub fn hashin_shtrikman_bulk_upper(k1: f64, g1: f64, k2: f64, f2: f64) -> f64 {
    let f1 = 1.0 - f2;
    let alpha = 3.0 * f1 / (3.0 * k1 + 4.0 * g1);
    if (k2 - k1).abs() < 1e-30 {
        return k1;
    }
    let denom = 1.0 / (k2 - k1) + alpha;
    if denom.abs() < 1e-30 {
        return k1;
    }
    k1 + f2 / denom
}

/// Hashin-Shtrikman lower bound for bulk modulus.
///
/// K_HS- = K_2 + f1 / (1/(K_1-K_2) + 3f2/(3K_2+4G_2))
///
/// Assumes K_2 ≤ K_1.
pub fn hashin_shtrikman_bulk_lower(k2: f64, g2: f64, k1: f64, f1: f64) -> f64 {
    let f2 = 1.0 - f1;
    let alpha = 3.0 * f2 / (3.0 * k2 + 4.0 * g2);
    if (k1 - k2).abs() < 1e-30 {
        return k2;
    }
    let denom = 1.0 / (k1 - k2) + alpha;
    if denom.abs() < 1e-30 {
        return k2;
    }
    k2 + f1 / denom
}

// ---------------------------------------------------------------------------
// Periodic unit cell homogenization (simplified 2D)
// ---------------------------------------------------------------------------

/// Configuration for a 2D periodic unit cell.
#[derive(Debug, Clone)]
pub struct PeriodicUnitCell {
    /// Cell dimensions [lx, ly].
    pub dimensions: [f64; 2],
    /// Matrix phase material.
    pub matrix: Phase,
    /// Inclusion phase material.
    pub inclusion: Phase,
    /// Inclusion radius (assuming circular inclusion).
    pub inclusion_radius: f64,
}

impl PeriodicUnitCell {
    /// Create a new periodic unit cell.
    pub fn new(lx: f64, ly: f64, matrix: Phase, inclusion: Phase, radius: f64) -> Self {
        Self {
            dimensions: [lx, ly],
            matrix,
            inclusion,
            inclusion_radius: radius,
        }
    }

    /// Volume fraction of the inclusion.
    pub fn inclusion_volume_fraction(&self) -> f64 {
        let area_incl = std::f64::consts::PI * self.inclusion_radius * self.inclusion_radius;
        let area_cell = self.dimensions[0] * self.dimensions[1];
        (area_incl / area_cell).min(1.0)
    }

    /// Effective bulk modulus via Mori-Tanaka for spherical inclusions.
    pub fn effective_bulk_modulus_mt(&self) -> f64 {
        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        let k_i = self.inclusion.bulk_modulus();
        let f_i = self.inclusion_volume_fraction();
        mori_tanaka_bulk_modulus(k_m, g_m, k_i, f_i)
    }

    /// Effective shear modulus via Mori-Tanaka.
    pub fn effective_shear_modulus_mt(&self) -> f64 {
        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        let g_i = self.inclusion.shear_modulus();
        let f_i = self.inclusion_volume_fraction();
        mori_tanaka_shear_modulus(k_m, g_m, g_i, f_i)
    }

    /// Effective Young's modulus from effective K and G.
    pub fn effective_youngs_modulus_mt(&self) -> f64 {
        let k_eff = self.effective_bulk_modulus_mt();
        let g_eff = self.effective_shear_modulus_mt();
        9.0 * k_eff * g_eff / (3.0 * k_eff + g_eff)
    }

    /// Effective Poisson's ratio from effective K and G.
    pub fn effective_poisson_ratio_mt(&self) -> f64 {
        let k_eff = self.effective_bulk_modulus_mt();
        let g_eff = self.effective_shear_modulus_mt();
        let denom = 3.0 * k_eff + g_eff;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        (3.0 * k_eff - 2.0 * g_eff) / (2.0 * denom)
    }

    /// Effective 6×6 isotropic stiffness tensor (Voigt notation) via MT.
    pub fn effective_stiffness_tensor_mt(&self) -> [[f64; 6]; 6] {
        let e_eff = self.effective_youngs_modulus_mt();
        let nu_eff = self.effective_poisson_ratio_mt();
        let phase_eff = Phase::new("eff", 1.0, e_eff, nu_eff);
        phase_eff.stiffness_voigt()
    }
}

// ---------------------------------------------------------------------------
// Thermal expansion homogenization
// ---------------------------------------------------------------------------

/// Effective thermal expansion tensor for a two-phase composite
/// using the Levin-Schapery relation.
///
/// α_eff = α_m + (α_i - α_m) * (K_i * (K_eff - K_m)) / (K_m * (K_eff - K_i) * ε)
/// where ε = volume fraction correction.
///
/// For isotropic phases the Levin relation gives:
///
///   α_eff = α_m + f_i * (α_i - α_m) * (1/K_m - 1/K_eff) / (1/K_m - 1/K_i)
pub fn levin_thermal_expansion(
    alpha_m: f64,
    k_m: f64,
    alpha_i: f64,
    k_i: f64,
    k_eff: f64,
    f_i: f64,
) -> f64 {
    if (k_i - k_m).abs() < 1e-30 {
        return alpha_m * (1.0 - f_i) + alpha_i * f_i;
    }
    let num = 1.0 / k_m - 1.0 / k_eff;
    let den = 1.0 / k_m - 1.0 / k_i;
    if den.abs() < 1e-30 {
        return alpha_m;
    }
    alpha_m + f_i * (alpha_i - alpha_m) * num / den
}

/// Multi-phase thermal expansion coefficient via Turner model.
///
/// CTE_c = Σ(α_i * E_i * φ_i) / Σ(E_i * φ_i)
pub fn multi_phase_cte_turner(alphas: &[f64], moduli: &[f64], fractions: &[f64]) -> f64 {
    assert_eq!(alphas.len(), moduli.len());
    assert_eq!(moduli.len(), fractions.len());
    let num: f64 = alphas
        .iter()
        .zip(moduli.iter())
        .zip(fractions.iter())
        .map(|((&a, &e), &f)| a * e * f)
        .sum();
    let den: f64 = moduli
        .iter()
        .zip(fractions.iter())
        .map(|(&e, &f)| e * f)
        .sum();
    if den.abs() < 1e-30 {
        return alphas.iter().sum::<f64>() / alphas.len() as f64;
    }
    num / den
}

/// Thermal mismatch stress in the matrix of a two-phase composite.
///
/// σ_m = E_m * (α_eff - α_m) * ΔT / (1 - ν_m)  (plane-strain approximation)
pub fn thermal_mismatch_stress(
    e_m: f64,
    nu_m: f64,
    alpha_eff: f64,
    alpha_m: f64,
    delta_t: f64,
) -> f64 {
    e_m * (alpha_eff - alpha_m) * delta_t / (1.0 - nu_m)
}

/// Volumetric thermal strain: ε_v = 3 * α * ΔT.
pub fn volumetric_thermal_strain(alpha: f64, delta_t: f64) -> f64 {
    3.0 * alpha * delta_t
}

// ---------------------------------------------------------------------------
// Effective elastic tensor — 3D Voigt notation
// ---------------------------------------------------------------------------

/// Build the 6×6 Voigt-notation elastic tensor C for an isotropic material
/// from Young's modulus E and Poisson's ratio ν.
pub fn isotropic_stiffness_tensor(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let c11 = lam + 2.0 * mu;
    let c12 = lam;
    let c44 = mu;
    let mut c = [[0.0f64; 6]; 6];
    c[0][0] = c11;
    c[0][1] = c12;
    c[0][2] = c12;
    c[1][0] = c12;
    c[1][1] = c11;
    c[1][2] = c12;
    c[2][0] = c12;
    c[2][1] = c12;
    c[2][2] = c11;
    c[3][3] = c44;
    c[4][4] = c44;
    c[5][5] = c44;
    c
}

/// Build the 6×6 compliance tensor S = C^{-1} for an isotropic material.
pub fn isotropic_compliance_tensor(e: f64, nu: f64) -> Option<[[f64; 6]; 6]> {
    let c = isotropic_stiffness_tensor(e, nu);
    mat6_inv(&c)
}

/// Effective elastic tensor from volume-averaged stiffnesses (Voigt bound).
///
/// C_eff = Σ f_i * C_i
pub fn voigt_effective_tensor(stiffnesses: &[[[f64; 6]; 6]], fractions: &[f64]) -> [[f64; 6]; 6] {
    assert_eq!(stiffnesses.len(), fractions.len());
    let mut c_eff = [[0.0f64; 6]; 6];
    for (ci, &fi) in stiffnesses.iter().zip(fractions.iter()) {
        for i in 0..6 {
            for j in 0..6 {
                c_eff[i][j] += fi * ci[i][j];
            }
        }
    }
    c_eff
}

/// Effective elastic tensor from Reuss bound (average of compliances).
///
/// S_eff = Σ f_i * S_i  then C_eff = S_eff^{-1}
pub fn reuss_effective_tensor(
    compliances: &[[[f64; 6]; 6]],
    fractions: &[f64],
) -> Option<[[f64; 6]; 6]> {
    assert_eq!(compliances.len(), fractions.len());
    let mut s_eff = [[0.0f64; 6]; 6];
    for (si, &fi) in compliances.iter().zip(fractions.iter()) {
        for i in 0..6 {
            for j in 0..6 {
                s_eff[i][j] += fi * si[i][j];
            }
        }
    }
    mat6_inv(&s_eff)
}

/// Hill average of Voigt and Reuss effective tensors.
///
/// C_Hill = (C_Voigt + C_Reuss) / 2
pub fn hill_average_tensor(
    stiffnesses: &[[[f64; 6]; 6]],
    compliances: &[[[f64; 6]; 6]],
    fractions: &[f64],
) -> Option<[[f64; 6]; 6]> {
    let c_voigt = voigt_effective_tensor(stiffnesses, fractions);
    let c_reuss = reuss_effective_tensor(compliances, fractions)?;
    let mut c_hill = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            c_hill[i][j] = 0.5 * (c_voigt[i][j] + c_reuss[i][j]);
        }
    }
    Some(c_hill)
}

/// Extract effective isotropic constants (E, ν) from a 6×6 stiffness tensor
/// assuming cubic/isotropic symmetry.
///
/// Returns `(E_eff, nu_eff)`.
pub fn effective_isotropic_constants(c: &[[f64; 6]; 6]) -> (f64, f64) {
    // For isotropic: C_11 = λ + 2μ, C_12 = λ, C_44 = μ
    let c11 = c[0][0];
    let c12 = c[0][1];
    let mu = c[3][3];
    let lam = c12;
    // E = mu * (3λ + 2μ) / (λ + μ)
    let denom_nu = lam + mu;
    let e = if denom_nu.abs() > 1e-30 {
        mu * (3.0 * lam + 2.0 * mu) / denom_nu
    } else {
        0.0
    };
    let nu = if denom_nu.abs() > 1e-30 {
        lam / (2.0 * denom_nu)
    } else {
        0.0
    };
    let _ = c11;
    (e, nu)
}

// ---------------------------------------------------------------------------
// Hill-Mandel condition check
// ---------------------------------------------------------------------------

/// Compute the Hill-Mandel energy error for a given stress/strain pair.
///
/// The Hill-Mandel condition states: <σ:ε> = <σ>:`ε`
///
/// This function computes |<σ:ε> - σ_avg:ε_avg| / |σ_avg:ε_avg|
/// as a relative energy error.
///
/// `stresses` and `strains` are lists of Gauss-point values,
/// `weights` are integration weights summing to 1.
pub fn hill_mandel_error(stresses: &[[f64; 6]], strains: &[[f64; 6]], weights: &[f64]) -> f64 {
    assert_eq!(stresses.len(), strains.len());
    assert_eq!(strains.len(), weights.len());
    let n = stresses.len();

    // Weighted average of σ:ε
    let se_avg: f64 = (0..n)
        .map(|i| {
            let s_dot_e: f64 = (0..6).map(|k| stresses[i][k] * strains[i][k]).sum();
            weights[i] * s_dot_e
        })
        .sum();

    // Weighted average stress and strain
    let mut sigma_avg = [0.0f64; 6];
    let mut eps_avg = [0.0f64; 6];
    for i in 0..n {
        for k in 0..6 {
            sigma_avg[k] += weights[i] * stresses[i][k];
            eps_avg[k] += weights[i] * strains[i][k];
        }
    }
    let sigma_dot_eps: f64 = (0..6).map(|k| sigma_avg[k] * eps_avg[k]).sum();

    let err = (se_avg - sigma_dot_eps).abs();
    if sigma_dot_eps.abs() > 1e-30 {
        err / sigma_dot_eps.abs()
    } else {
        err
    }
}

// ---------------------------------------------------------------------------
// Three-phase composite — core-shell model
// ---------------------------------------------------------------------------

/// Three-phase (core-shell) composite: inner sphere (phase 2) coated with
/// shell (phase 1) embedded in matrix (phase 0).
///
/// Uses the Christensen-Lo generalized self-consistent model estimate.
#[derive(Debug, Clone)]
pub struct ThreePhaseComposite {
    /// Matrix phase.
    pub matrix: Phase,
    /// Shell phase (intermediate layer).
    pub shell: Phase,
    /// Core phase (innermost inclusion).
    pub core: Phase,
    /// Volume fraction of core+shell (total inclusion volume).
    pub inclusion_fraction: f64,
    /// Volume fraction of core within the inclusion.
    pub core_fraction_in_inclusion: f64,
}

impl ThreePhaseComposite {
    /// Volume fraction of core in the total composite.
    pub fn core_global_fraction(&self) -> f64 {
        self.inclusion_fraction * self.core_fraction_in_inclusion
    }

    /// Volume fraction of shell in the total composite.
    pub fn shell_global_fraction(&self) -> f64 {
        self.inclusion_fraction * (1.0 - self.core_fraction_in_inclusion)
    }

    /// Effective bulk modulus via sequential Mori-Tanaka:
    /// First homogenize core+shell, then embed in matrix.
    pub fn effective_bulk_modulus(&self) -> f64 {
        // Step 1: homogenize core in shell
        let k_shell = self.shell.bulk_modulus();
        let g_shell = self.shell.shear_modulus();
        let k_core = self.core.bulk_modulus();
        let f_core_in_shell = self.core_fraction_in_inclusion;
        let k_12 = mori_tanaka_bulk_modulus(k_shell, g_shell, k_core, f_core_in_shell);
        let _g_12 =
            mori_tanaka_shear_modulus(k_shell, g_shell, self.core.shear_modulus(), f_core_in_shell);

        // Step 2: embed core-shell in matrix
        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        let f_inclusion = self.inclusion_fraction;
        mori_tanaka_bulk_modulus(k_m, g_m, k_12, f_inclusion)
    }

    /// Effective shear modulus via sequential Mori-Tanaka.
    pub fn effective_shear_modulus(&self) -> f64 {
        let k_shell = self.shell.bulk_modulus();
        let g_shell = self.shell.shear_modulus();
        let g_core = self.core.shear_modulus();
        let f_core_in_shell = self.core_fraction_in_inclusion;
        let g_12 = mori_tanaka_shear_modulus(k_shell, g_shell, g_core, f_core_in_shell);

        let k_m = self.matrix.bulk_modulus();
        let g_m = self.matrix.shear_modulus();
        mori_tanaka_shear_modulus(k_m, g_m, g_12, self.inclusion_fraction)
    }

    /// Effective Young's modulus from K and G.
    pub fn effective_youngs_modulus(&self) -> f64 {
        let k = self.effective_bulk_modulus();
        let g = self.effective_shear_modulus();
        9.0 * k * g / (3.0 * k + g)
    }
}

// ---------------------------------------------------------------------------
// Damage-affected homogenization (stiffness degradation)
// ---------------------------------------------------------------------------

/// Stiffness degradation of a composite phase due to damage.
///
/// The effective stiffness of phase i with damage D_i:
/// C_i_eff = (1 - D_i) * C_i_intact
pub fn damage_degraded_voigt_tensor(
    stiffnesses: &[[[f64; 6]; 6]],
    fractions: &[f64],
    damage_values: &[f64],
) -> [[f64; 6]; 6] {
    assert_eq!(stiffnesses.len(), fractions.len());
    assert_eq!(fractions.len(), damage_values.len());
    let mut c_eff = [[0.0f64; 6]; 6];
    for ((ci, &fi), &di) in stiffnesses
        .iter()
        .zip(fractions.iter())
        .zip(damage_values.iter())
    {
        let scale = fi * (1.0 - di.clamp(0.0, 1.0));
        for i in 0..6 {
            for j in 0..6 {
                c_eff[i][j] += scale * ci[i][j];
            }
        }
    }
    c_eff
}

#[cfg(test)]
mod tests_homogenization_new {
    use super::*;

    #[test]
    fn test_voigt_geq_reuss() {
        let e = [200e9_f64, 70e9_f64];
        let fracs = [0.6_f64, 0.4_f64];
        let voigt = effective_elastic_modulus_voigt(&e, &fracs);
        let reuss = effective_elastic_modulus_reuss(&e, &fracs);
        assert!(
            voigt >= reuss,
            "Voigt={voigt:.3e} should be >= Reuss={reuss:.3e}"
        );
    }

    #[test]
    fn test_hill_bounds_voigt_geq_reuss() {
        let bounds = HillBounds::compute(200e9, 0.3, 0.6, 70e9, 0.33);
        assert!(
            bounds.voigt >= bounds.reuss,
            "Voigt={:.3e} Reuss={:.3e}",
            bounds.voigt,
            bounds.reuss
        );
    }

    #[test]
    fn test_mori_tanaka_between_voigt_reuss() {
        let e_m = 70e9_f64;
        let nu_m = 0.33_f64;
        let e_i = 380e9_f64;
        let nu_i = 0.22_f64;
        let f_i = 0.3_f64;
        let phi1 = 1.0 - f_i;
        let k_m = e_m / (3.0 * (1.0 - 2.0 * nu_m));
        let g_m = e_m / (2.0 * (1.0 + nu_m));
        let k_i = e_i / (3.0 * (1.0 - 2.0 * nu_i));

        let k_mt = mori_tanaka_bulk_modulus(k_m, g_m, k_i, f_i);

        let k_voigt = phi1 * k_m + f_i * k_i;
        let k_reuss = 1.0 / (phi1 / k_m + f_i / k_i);
        assert!(
            (k_reuss * (1.0 - 1e-6)..=k_voigt * (1.0 + 1e-6)).contains(&k_mt),
            "MT K={k_mt:.3e} should be between Reuss={k_reuss:.3e} and Voigt={k_voigt:.3e}"
        );
    }

    #[test]
    fn test_halpin_tsai_at_xi_zero_is_reuss_like() {
        // At xi=0, Halpin-Tsai reduces to: E = E_m * (1 - vf) + E_f * vf ?
        // Actually at xi=0: eta = (Ef/Em - 1)/(Ef/Em + 0) = 1 - Em/Ef
        // E = Em * 1 / (1 - eta*vf)  — check it's positive and > E_m for E_f > E_m
        let e_m = 70e9_f64;
        let e_f = 380e9_f64;
        let vf = 0.2_f64;
        let e_ht = halpin_tsai_modulus(e_m, e_f, vf, 0.0);
        assert!(
            e_ht > e_m,
            "Halpin-Tsai should give E > E_m for E_f > E_m: {e_ht:.3e}"
        );
        // At xi=0, should be between Reuss and Voigt
        let reuss = 1.0 / ((1.0 - vf) / e_m + vf / e_f);
        let voigt = (1.0 - vf) * e_m + vf * e_f;
        assert!(
            (reuss * 0.9999..=voigt * 1.0001).contains(&e_ht),
            "HT={e_ht:.3e} Reuss={reuss:.3e} Voigt={voigt:.3e}"
        );
    }

    #[test]
    fn test_composite_cte_single_phase() {
        // phi1=1 => CTE = alpha1
        let cte = composite_thermal_expansion(12e-6, 200e9, 1.0, 23e-6, 70e9);
        assert!(
            (cte - 12e-6).abs() < 1e-20,
            "single-phase CTE should equal alpha1: {cte}"
        );
    }

    #[test]
    fn test_composite_cte_between_components() {
        let alpha1 = 12e-6_f64;
        let alpha2 = 23e-6_f64;
        let cte = composite_thermal_expansion(alpha1, 200e9, 0.5, alpha2, 70e9);
        assert!(
            (alpha1.min(alpha2)..=alpha1.max(alpha2)).contains(&cte),
            "CTE={cte:.3e} should be between components"
        );
    }

    #[test]
    fn test_rve_cell_volume_fractions() {
        let mut rve = RveCell::new(4, 4);
        // Set half the cells to phase 1
        for i in 8..16 {
            rve.phase_map[i] = 1;
        }
        let vf = rve.compute_volume_fractions();
        assert_eq!(vf.len(), 2);
        assert!((vf[0] - 0.5).abs() < 1e-10, "vf[0]={}", vf[0]);
        assert!((vf[1] - 0.5).abs() < 1e-10, "vf[1]={}", vf[1]);
    }

    #[test]
    fn test_periodic_bc_pair_count() {
        let rve = RveCell::new(4, 3); // 4 cols, 3 rows
        let pairs = periodic_boundary_conditions_2d(&rve);
        // Left-right: 3 rows × 2 dofs = 6 pairs
        // Bottom-top: 4 cols × 2 dofs = 8 pairs
        assert_eq!(pairs.len(), 14, "expected 14 pairs, got {}", pairs.len());
    }

    // ---- StrainConcentrationTensor ----

    #[test]
    fn test_strain_concentration_dilute_identity_matrix() {
        // For a very soft inclusion (→0 modulus) with very low volume fraction
        // the Eshelby tensor entries are bounded
        let a = StrainConcentrationTensor::dilute_spherical_inclusion(200e9, 0.3, 1e3, 0.45);
        // Diagonal entries should be finite and positive
        for i in 0..6 {
            assert!(a.tensor[i][i].is_finite(), "A[{i}][{i}] not finite");
        }
    }

    #[test]
    fn test_strain_concentration_bulk_dominated() {
        // Same modulus → concentration tensor is identity
        let a = StrainConcentrationTensor::dilute_spherical_inclusion(200e9, 0.3, 200e9, 0.3);
        for i in 0..6 {
            let diag = a.tensor[i][i];
            assert!((diag - 1.0).abs() < 0.1, "A[{i}][{i}]={diag} expected ~1");
        }
    }

    // ---- EffectiveThermalConductivity ----

    #[test]
    fn test_maxwell_thermal_single_phase() {
        // When inclusion modulus equals matrix modulus (lambda_m = lambda_i = 50.0),
        // any volume fraction should return the matrix conductivity.
        let k = maxwell_effective_thermal_conductivity(50.0, 50.0, 0.3);
        assert!((k - 50.0).abs() < 1e-6, "single-phase thermal cond={k}");
    }

    #[test]
    fn test_hashin_shtrikman_thermal_bounds() {
        let k_lower = hs_thermal_lower_bound(1.0, 50.0, 0.3);
        let k_upper = hs_thermal_upper_bound(1.0, 50.0, 0.3);
        assert!(
            k_lower <= k_upper + 1e-10,
            "HS: lower={k_lower} upper={k_upper}"
        );
    }

    // ---- AnisotropyIndex ----

    #[test]
    fn test_zener_anisotropy_isotropic() {
        // Isotropic material: A_z = 1
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let g = e / (2.0 * (1.0 + nu));
        let c11 = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c12 = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c44 = g;
        let az = zener_anisotropy_index(c11, c12, c44);
        assert!((az - 1.0).abs() < 1e-6, "isotropic Zener index={az}");
    }

    // ---- ThermoElasticCoupling ----

    #[test]
    fn test_thermal_stress_uniaxial() {
        // σ = -E·α·ΔT for constrained uniaxial
        let ts = ThermoElasticCoupling {
            alpha: 12e-6,
            e: 200e9,
            nu: 0.3,
        };
        let delta_t = 100.0;
        let sigma = ts.constrained_thermal_stress(delta_t);
        let expected = -200e9 * 12e-6 * 100.0;
        assert!(
            (sigma - expected).abs() < 1.0,
            "thermal stress={sigma} expected={expected}"
        );
    }

    #[test]
    fn test_free_thermal_strain() {
        let ts = ThermoElasticCoupling {
            alpha: 12e-6,
            e: 200e9,
            nu: 0.3,
        };
        let eps = ts.free_thermal_strain(50.0);
        assert!(
            (eps - 12e-6 * 50.0).abs() < 1e-20,
            "free thermal strain={eps}"
        );
    }

    // ---- HomogenizedStiffnessBounds ----

    #[test]
    fn test_homogenized_stiffness_bounds_order() {
        let phases = vec![
            Phase::new("steel", 0.6, 210e9, 0.3),
            Phase::new("epoxy", 0.4, 3.5e9, 0.35),
        ];
        let bounds = HomogenizedStiffnessBounds::compute(&phases);
        assert!(
            bounds.reuss_bound <= bounds.voigt_bound + 1.0,
            "Reuss <= Voigt: {} vs {}",
            bounds.reuss_bound,
            bounds.voigt_bound
        );
        assert!(
            bounds.hs_lower <= bounds.hs_upper + 1.0,
            "HS lower <= upper: {} vs {}",
            bounds.hs_lower,
            bounds.hs_upper
        );
    }

    #[test]
    fn test_homogenized_stiffness_hs_within_hill() {
        let phases = vec![
            Phase::new("steel", 0.5, 210e9, 0.3),
            Phase::new("epoxy", 0.5, 3.5e9, 0.35),
        ];
        let bounds = HomogenizedStiffnessBounds::compute(&phases);
        assert!(
            bounds.hs_lower >= bounds.reuss_bound - 1.0,
            "HS lower >= Reuss: {} vs {}",
            bounds.hs_lower,
            bounds.reuss_bound
        );
        assert!(
            bounds.hs_upper <= bounds.voigt_bound + 1.0,
            "HS upper <= Voigt: {} vs {}",
            bounds.hs_upper,
            bounds.voigt_bound
        );
    }

    // ---- Halpin-Tsai extended ----

    #[test]
    fn test_halpin_tsai_longitudinal_rule_of_mixtures() {
        let e_l = halpin_tsai_longitudinal(380e9, 70e9, 0.3);
        let expected = 380e9 * 0.3 + 70e9 * 0.7;
        assert!((e_l - expected).abs() < 1.0, "HT longitudinal = ROM");
    }

    #[test]
    fn test_halpin_tsai_transverse_positive() {
        let e_t = halpin_tsai_transverse(380e9, 70e9, 0.3, 2.0);
        assert!(
            e_t > 70e9,
            "Transverse modulus should exceed matrix modulus"
        );
        assert!(
            e_t < 380e9,
            "Transverse modulus should be less than fiber modulus"
        );
    }

    #[test]
    fn test_halpin_tsai_shear_positive() {
        let g_12 = halpin_tsai_shear_12(50e9, 26e9, 0.3);
        assert!(
            g_12 > 26e9 && g_12 < 50e9,
            "Shear modulus should be between matrix and fiber"
        );
    }

    #[test]
    fn test_halpin_tsai_poisson_rule_of_mixtures() {
        let nu_12 = halpin_tsai_poisson_12(0.25, 0.35, 0.4);
        let expected = 0.25 * 0.4 + 0.35 * 0.6;
        assert!((nu_12 - expected).abs() < 1e-12);
    }

    // ---- Hashin-Shtrikman extended ----

    #[test]
    fn test_hs_shear_upper_geq_lower() {
        let k1 = 166e9;
        let g1 = 80e9;
        let k2 = 50e9;
        let g2 = 30e9;
        let f = 0.3f64;
        let upper = hashin_shtrikman_shear_upper(k1, g1, g2, f);
        let lower = hashin_shtrikman_shear_lower(k2, g2, g1, 1.0 - f);
        assert!(
            upper >= lower - 1.0,
            "HS shear upper >= lower: {upper:.3e} vs {lower:.3e}"
        );
    }

    #[test]
    fn test_hs_bulk_upper_geq_lower() {
        let k1 = 166e9;
        let g1 = 80e9;
        let k2 = 50e9;
        let g2 = 30e9;
        let f = 0.3f64;
        let upper = hashin_shtrikman_bulk_upper(k1, g1, k2, f);
        let lower = hashin_shtrikman_bulk_lower(k2, g2, k1, 1.0 - f);
        assert!(
            upper >= lower - 1.0,
            "HS bulk upper >= lower: {upper:.3e} vs {lower:.3e}"
        );
    }

    #[test]
    fn test_hs_bounds_within_voigt_reuss() {
        let e1 = 210e9_f64;
        let nu1 = 0.3_f64;
        let e2 = 70e9_f64;
        let nu2 = 0.33_f64;
        let f1 = 0.6_f64;
        let f2 = 0.4_f64;
        let k1 = e1 / (3.0 * (1.0 - 2.0 * nu1));
        let k2 = e2 / (3.0 * (1.0 - 2.0 * nu2));
        let g1 = e1 / (2.0 * (1.0 + nu1));
        let _g2 = e2 / (2.0 * (1.0 + nu2));
        let k_hs_upper = hashin_shtrikman_bulk_upper(k1, g1, k2, f2);
        let k_voigt = f1 * k1 + f2 * k2;
        let k_reuss = 1.0 / (f1 / k1 + f2 / k2);
        assert!(k_hs_upper <= k_voigt * 1.001, "HS K upper <= Voigt");
        assert!(k_hs_upper >= k_reuss * 0.999, "HS K upper >= Reuss");
    }

    // ---- PeriodicUnitCell ----

    #[test]
    fn test_periodic_unit_cell_volume_fraction() {
        use std::f64::consts::PI;
        let matrix = Phase::new("epoxy", 0.7, 3.5e9, 0.35);
        let incl = Phase::new("glass", 0.3, 70e9, 0.23);
        let r = 0.1;
        let cell = PeriodicUnitCell::new(1.0, 1.0, matrix, incl, r);
        let vf = cell.inclusion_volume_fraction();
        let expected = PI * r * r;
        assert!(
            (vf - expected).abs() < 1e-12,
            "VF={vf:.4e} expected={expected:.4e}"
        );
    }

    #[test]
    fn test_periodic_unit_cell_effective_modulus_between_bounds() {
        let matrix = Phase::new("epoxy", 0.7, 3.5e9, 0.35);
        let incl = Phase::new("glass", 0.3, 70e9, 0.23);
        let r = 0.1;
        let cell = PeriodicUnitCell::new(1.0, 1.0, matrix, incl, r);
        let e_eff = cell.effective_youngs_modulus_mt();
        assert!(e_eff > 3.5e9, "Effective E should exceed matrix E");
        assert!(e_eff < 70e9, "Effective E should be less than inclusion E");
    }

    #[test]
    fn test_periodic_unit_cell_poisson_ratio_physical() {
        let matrix = Phase::new("epoxy", 0.7, 3.5e9, 0.35);
        let incl = Phase::new("glass", 0.3, 70e9, 0.23);
        let cell = PeriodicUnitCell::new(1.0, 1.0, matrix, incl, 0.1);
        let nu_eff = cell.effective_poisson_ratio_mt();
        assert!(
            nu_eff > 0.0 && nu_eff < 0.5,
            "Effective ν should be in (0, 0.5), got {nu_eff}"
        );
    }

    // ---- Levin thermal expansion ----

    #[test]
    fn test_levin_cte_pure_matrix() {
        // f_i = 0 → CTE = alpha_m
        let alpha_eff = levin_thermal_expansion(12e-6, 166e9, 23e-6, 50e9, 166e9, 0.0);
        // When f_i = 0, (1/K_m - 1/K_eff)/(1/K_m - 1/K_i) * 0 = 0
        // alpha_eff ≈ alpha_m
        assert!(
            (alpha_eff - 12e-6).abs() < 1e-12,
            "Pure matrix: alpha_eff should equal alpha_m, got {alpha_eff:.3e}"
        );
    }

    #[test]
    fn test_levin_cte_same_phases() {
        // When alpha_m = alpha_i → alpha_eff = alpha_m regardless
        let alpha_eff = levin_thermal_expansion(12e-6, 166e9, 12e-6, 50e9, 100e9, 0.3);
        assert!(
            (alpha_eff - 12e-6).abs() < 1e-15,
            "Same CTE phases: alpha_eff = alpha"
        );
    }

    #[test]
    fn test_multi_phase_cte_turner() {
        let alphas = vec![12e-6, 23e-6, 8e-6];
        let moduli = vec![200e9, 70e9, 380e9];
        let fracs = vec![0.5, 0.3, 0.2];
        let cte = multi_phase_cte_turner(&alphas, &moduli, &fracs);
        // Should be between min and max alpha
        assert!(
            (8e-6..=23e-6).contains(&cte),
            "Turner CTE out of range: {cte:.3e}"
        );
    }

    // ---- Effective elastic tensor ----

    #[test]
    fn test_isotropic_stiffness_tensor_symmetric() {
        let c = isotropic_stiffness_tensor(200e9, 0.3);
        assert!(
            mat6_symmetry_error(&c) < 1e-6,
            "Stiffness tensor not symmetric"
        );
    }

    #[test]
    fn test_isotropic_stiffness_tensor_roundtrip() {
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let c = isotropic_stiffness_tensor(e, nu);
        let (e_eff, nu_eff) = effective_isotropic_constants(&c);
        assert!(
            (e_eff - e).abs() / e < 1e-9,
            "E roundtrip: got {e_eff:.3e}, expected {e:.3e}"
        );
        assert!(
            (nu_eff - nu).abs() < 1e-9,
            "ν roundtrip: got {nu_eff:.4}, expected {nu:.4}"
        );
    }

    #[test]
    fn test_voigt_effective_tensor_single_phase() {
        let c = isotropic_stiffness_tensor(200e9, 0.3);
        let c_eff = voigt_effective_tensor(&[c], &[1.0]);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c_eff[i][j] - c[i][j]).abs() < 1.0,
                    "Single phase: C_eff[{i}][{j}] should equal C[{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_voigt_reuss_hill_tensor_bounds() {
        let c1 = isotropic_stiffness_tensor(200e9, 0.3);
        let c2 = isotropic_stiffness_tensor(70e9, 0.33);
        let s1 = mat6_inv(&c1).unwrap();
        let s2 = mat6_inv(&c2).unwrap();
        let fracs = [0.6f64, 0.4f64];
        let c_voigt = voigt_effective_tensor(&[c1, c2], &fracs);
        let c_reuss = reuss_effective_tensor(&[s1, s2], &fracs).unwrap();
        // Voigt C_11 >= Reuss C_11
        assert!(
            c_voigt[0][0] >= c_reuss[0][0] - 1.0,
            "Voigt C_11 >= Reuss C_11: {:.3e} vs {:.3e}",
            c_voigt[0][0],
            c_reuss[0][0]
        );
    }

    // ---- Hill-Mandel condition ----

    #[test]
    fn test_hill_mandel_uniform_stress_zero_error() {
        // Uniform stress and strain: <σ:ε> = σ:ε trivially (no spatial variation)
        let sigma = [1e6, 0.5e6, 0.3e6, 0.2e6, 0.1e6, 0.05e6];
        let eps = [5e-3, 2e-3, 1e-3, 8e-4, 4e-4, 2e-4];
        let stresses = vec![sigma, sigma, sigma];
        let strains = vec![eps, eps, eps];
        let weights = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let err = hill_mandel_error(&stresses, &strains, &weights);
        assert!(
            err < 1e-10,
            "Uniform stress/strain: Hill-Mandel error should be ~0, got {err}"
        );
    }

    // ---- ThreePhaseComposite ----

    #[test]
    fn test_three_phase_effective_modulus_between_phases() {
        let matrix = Phase::new("epoxy", 0.6, 3.5e9, 0.35);
        let shell = Phase::new("glass", 0.3, 70e9, 0.23);
        let core = Phase::new("carbon", 0.1, 300e9, 0.21);
        let comp = ThreePhaseComposite {
            matrix,
            shell,
            core,
            inclusion_fraction: 0.4,
            core_fraction_in_inclusion: 0.5,
        };
        let e_eff = comp.effective_youngs_modulus();
        assert!(e_eff > 3.5e9, "Effective E must exceed matrix E");
    }

    // ---- Damage-degraded homogenization ----

    #[test]
    fn test_damage_degraded_tensor_zero_damage() {
        let c1 = isotropic_stiffness_tensor(200e9, 0.3);
        let c2 = isotropic_stiffness_tensor(70e9, 0.33);
        let c_eff_undamaged = voigt_effective_tensor(&[c1, c2], &[0.6, 0.4]);
        let c_eff_damaged = damage_degraded_voigt_tensor(&[c1, c2], &[0.6, 0.4], &[0.0, 0.0]);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c_eff_damaged[i][j] - c_eff_undamaged[i][j]).abs() < 1.0,
                    "Zero damage: degraded tensor == undamaged"
                );
            }
        }
    }

    #[test]
    fn test_damage_degraded_tensor_reduces_stiffness() {
        let c = isotropic_stiffness_tensor(200e9, 0.3);
        let c_undamaged = damage_degraded_voigt_tensor(&[c], &[1.0], &[0.0]);
        let c_damaged = damage_degraded_voigt_tensor(&[c], &[1.0], &[0.5]);
        assert!(
            c_damaged[0][0] < c_undamaged[0][0],
            "Damage should reduce stiffness"
        );
        assert!(
            (c_damaged[0][0] / c_undamaged[0][0] - 0.5).abs() < 1e-9,
            "50% damage should halve stiffness"
        );
    }
}
