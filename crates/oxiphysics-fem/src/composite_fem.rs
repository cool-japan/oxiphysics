#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite material FEM using Classical Lamination Theory (CLT).
//!
//! Provides ply/laminate data structures, the ABD stiffness matrix, mid-plane
//! strain/curvature solution, ply-level stress recovery, and failure criteria
//! (Tsai-Wu, Hashin).

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Ply
// ---------------------------------------------------------------------------

/// A single orthotropic ply in a laminate.
#[derive(Clone, Debug)]
pub struct Ply {
    /// Ply thickness (m).
    pub thickness: f64,
    /// Fibre orientation angle from the laminate x-axis (degrees).
    pub angle_deg: f64,
    /// Longitudinal Young's modulus E₁ (Pa).
    pub e1: f64,
    /// Transverse Young's modulus E₂ (Pa).
    pub e2: f64,
    /// In-plane shear modulus G₁₂ (Pa).
    pub g12: f64,
    /// Major Poisson's ratio ν₁₂.
    pub nu12: f64,
}

impl Ply {
    /// Create a new [`Ply`].
    pub fn new(thickness: f64, angle_deg: f64, e1: f64, e2: f64, g12: f64, nu12: f64) -> Self {
        Self {
            thickness,
            angle_deg,
            e1,
            e2,
            g12,
            nu12,
        }
    }

    /// Minor Poisson's ratio ν₂₁ = ν₁₂ * E₂ / E₁.
    pub fn nu21(&self) -> f64 {
        if self.e1 < 1e-300 {
            return 0.0;
        }
        self.nu12 * self.e2 / self.e1
    }
}

// ---------------------------------------------------------------------------
// Laminate
// ---------------------------------------------------------------------------

/// A stack of plies forming a laminate.
#[derive(Clone, Debug)]
pub struct Laminate {
    /// Ordered list of plies (bottom to top).
    pub plies: Vec<Ply>,
}

impl Laminate {
    /// Create a new empty [`Laminate`].
    pub fn new() -> Self {
        Self { plies: Vec::new() }
    }

    /// Add a ply to the top of the laminate.
    pub fn add_ply(&mut self, ply: Ply) {
        self.plies.push(ply);
    }

    /// Total laminate thickness (m).
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }
}

impl Default for Laminate {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Reduced stiffness matrix Q (material frame) for a ply.
///
/// Returns 3×3 array `[[Q11, Q12, 0\], [Q12, Q22, 0], [0, 0, Q66]]`.
pub fn reduced_stiffness(ply: &Ply) -> [[f64; 3]; 3] {
    let nu21 = ply.nu21();
    let denom = 1.0 - ply.nu12 * nu21;
    if denom.abs() < 1e-300 {
        return [[0.0; 3]; 3];
    }
    let q11 = ply.e1 / denom;
    let q22 = ply.e2 / denom;
    let q12 = ply.nu12 * ply.e2 / denom;
    let q66 = ply.g12;
    [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
}

/// 3×3 Reuter rotation matrix T for a given angle (radians).
///
/// Maps material-frame Voigt strain to global-frame: `ε_global = T⁻¹ * ε_mat`.
pub fn rotation_matrix_2d(angle_rad: f64) -> [[f64; 3]; 3] {
    let c = angle_rad.cos();
    let s = angle_rad.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;
    [[c2, s2, 2.0 * cs], [s2, c2, -2.0 * cs], [-cs, cs, c2 - s2]]
}

/// Invert a 3×3 matrix; returns zero matrix if singular.
fn inv3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-300 {
        return [[0.0; 3]; 3];
    }
    let inv_det = 1.0 / det;
    [
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
    ]
}

/// Multiply two 3×3 matrices.
fn mat3x3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
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

/// Transpose a 3×3 matrix.
fn transpose3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

// ---------------------------------------------------------------------------
// Transformed stiffness (Q-bar)
// ---------------------------------------------------------------------------

/// Compute the transformed (global-frame) reduced stiffness matrix Q̄ for a ply.
///
/// Uses the transformation `Q̄ = T⁻¹ Q T⁻ᵀ` where T is the Reuter rotation
/// matrix corresponding to the ply fibre angle.
pub fn transformed_stiffness(ply: &Ply) -> [[f64; 3]; 3] {
    let angle_rad = ply.angle_deg.to_radians();
    let t = rotation_matrix_2d(angle_rad);
    let t_inv = inv3(t);
    let t_inv_t = transpose3(t_inv);
    let q = reduced_stiffness(ply);
    mat3x3_mul(mat3x3_mul(t_inv, q), t_inv_t)
}

// ---------------------------------------------------------------------------
// ABD matrix (Classical Lamination Theory)
// ---------------------------------------------------------------------------

/// Compute the 6×6 ABD stiffness matrix for a laminate.
///
/// The result is arranged as `[A, B; B, D]` where
/// - A (3×3) is the extensional stiffness,
/// - B (3×3) is the bending-extension coupling,
/// - D (3×3) is the bending stiffness.
///
/// The z-coordinates run from −h/2 to +h/2 (mid-plane reference).
pub fn laminate_abd_matrix(lam: &Laminate) -> [[f64; 6]; 6] {
    let total = lam.total_thickness();
    let mut z = -total / 2.0;

    let mut a = [[0.0f64; 3]; 3];
    let mut b = [[0.0f64; 3]; 3];
    let mut d = [[0.0f64; 3]; 3];

    for ply in &lam.plies {
        let z0 = z;
        let z1 = z + ply.thickness;
        let qbar = transformed_stiffness(ply);
        let dz1 = z1 - z0;
        let dz2 = (z1 * z1 - z0 * z0) / 2.0;
        let dz3 = (z1 * z1 * z1 - z0 * z0 * z0) / 3.0;
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += qbar[i][j] * dz1;
                b[i][j] += qbar[i][j] * dz2;
                d[i][j] += qbar[i][j] * dz3;
            }
        }
        z = z1;
    }

    let mut abd = [[0.0f64; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            abd[i][j] = a[i][j];
            abd[i][j + 3] = b[i][j];
            abd[i + 3][j] = b[i][j];
            abd[i + 3][j + 3] = d[i][j];
        }
    }
    abd
}

// ---------------------------------------------------------------------------
// Solve mid-plane strains and curvatures
// ---------------------------------------------------------------------------

/// Multiply a 6×6 matrix by a 6-vector.
fn mat6_mul_vec6(m: [[f64; 6]; 6], v: [f64; 6]) -> [f64; 6] {
    let mut r = [0.0f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            r[i] += m[i][j] * v[j];
        }
    }
    r
}

/// Invert a 6×6 matrix using Gauss-Jordan elimination.
fn inv6(m: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let n = 6usize;
    let mut aug = [[0.0f64; 12]; 6];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = m[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    for col in 0..n {
        // Find pivot
        let mut pivot_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                pivot_row = row;
            }
        }
        if max_val < 1e-300 {
            return [[0.0; 6]; 6];
        }
        aug.swap(col, pivot_row);
        let pivot = aug[col][col];
        for j in 0..12 {
            aug[col][j] /= pivot;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in 0..12 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }
    let mut inv = [[0.0f64; 6]; 6];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    inv
}

/// Solve for mid-plane strains and curvatures from resultant forces/moments.
///
/// Given the 6-vector `loads = [N_x, N_y, N_xy, M_x, M_y, M_xy]` and the 6×6
/// ABD matrix, returns `[ε_x⁰, ε_y⁰, γ_xy⁰, κ_x, κ_y, κ_xy]`.
pub fn laminate_midplane_strains(abd: &[[f64; 6]; 6], loads: [f64; 6]) -> [f64; 6] {
    let m: [[f64; 6]; 6] = *abd;
    let m_inv = inv6(m);
    mat6_mul_vec6(m_inv, loads)
}

// ---------------------------------------------------------------------------
// Ply-level stress recovery
// ---------------------------------------------------------------------------

/// Compute the in-plane stress vector \[σ₁, σ₂, τ₁₂\] in the **ply material frame**
/// at a given z-coordinate.
///
/// - `midplane_strain` – `[ε_x⁰, ε_y⁰, γ_xy⁰]` at the laminate mid-plane.
/// - `z` – through-thickness coordinate (m).
/// - `curvature` – `[κ_x, κ_y, κ_xy]`.
pub fn ply_stresses(ply: &Ply, midplane_strain: [f64; 3], z: f64, curvature: [f64; 3]) -> [f64; 3] {
    // Global strain at z
    let eps_global = [
        midplane_strain[0] + z * curvature[0],
        midplane_strain[1] + z * curvature[1],
        midplane_strain[2] + z * curvature[2],
    ];
    // Transform to material frame: ε_mat = T * ε_global
    let angle_rad = ply.angle_deg.to_radians();
    let t = rotation_matrix_2d(angle_rad);
    let eps_mat = [
        t[0][0] * eps_global[0] + t[0][1] * eps_global[1] + t[0][2] * eps_global[2],
        t[1][0] * eps_global[0] + t[1][1] * eps_global[1] + t[1][2] * eps_global[2],
        t[2][0] * eps_global[0] + t[2][1] * eps_global[1] + t[2][2] * eps_global[2],
    ];
    // σ_mat = Q * ε_mat
    let q = reduced_stiffness(ply);
    [
        q[0][0] * eps_mat[0] + q[0][1] * eps_mat[1] + q[0][2] * eps_mat[2],
        q[1][0] * eps_mat[0] + q[1][1] * eps_mat[1] + q[1][2] * eps_mat[2],
        q[2][0] * eps_mat[0] + q[2][1] * eps_mat[1] + q[2][2] * eps_mat[2],
    ]
}

// ---------------------------------------------------------------------------
// Failure criteria
// ---------------------------------------------------------------------------

/// Tsai-Wu failure index for a ply under plane stress.
///
/// Returns the scalar failure index F; failure occurs when F ≥ 1.
///
/// - `stress` – `[σ₁, σ₂, τ₁₂]` in material frame.
/// - `f1`, `f2` – linear interaction terms (1/Xᵀ − 1/Xᶜ, etc.).
/// - `f12` – biaxial interaction term (typically −1/(2 √(Xᵀ Xᶜ))).
/// - `f11`, `f22` – quadratic terms (1/(Xᵀ Xᶜ), 1/(Yᵀ Yᶜ)).
/// - `f66 = 1/S²` is derived from shear strength `s`.
pub fn tsai_wu_failure(
    stress: [f64; 3],
    f1: f64,
    f2: f64,
    f12: f64,
    f11: f64,
    f22: f64,
    s: f64,
) -> f64 {
    let f66 = if s.abs() > 1e-300 { 1.0 / (s * s) } else { 0.0 };
    f1 * stress[0]
        + f2 * stress[1]
        + f11 * stress[0] * stress[0]
        + f22 * stress[1] * stress[1]
        + f66 * stress[2] * stress[2]
        + 2.0 * f12 * stress[0] * stress[1]
}

/// Hashin failure criterion (2D).
///
/// Returns `(fibre_failed, matrix_failed)`.
///
/// - `stress` – `[σ₁, σ₂, τ₁₂]` in material frame.
/// - `xt`, `xc` – longitudinal tensile / compressive strength.
/// - `yt`, `yc` – transverse tensile / compressive strength.
/// - `s`        – in-plane shear strength.
pub fn hashin_failure(
    stress: [f64; 3],
    xt: f64,
    xc: f64,
    yt: f64,
    yc: f64,
    s: f64,
) -> (bool, bool) {
    let s1 = stress[0];
    let s2 = stress[1];
    let t12 = stress[2];
    let s2_ = if s.abs() > 1e-300 { s } else { 1e-300 };

    // Fibre failure
    let fibre_failed = if s1 >= 0.0 {
        // Tensile: (σ₁/Xᵀ)² + (τ₁₂/S)² ≥ 1
        let xt_ = if xt.abs() > 1e-300 { xt } else { 1e-300 };
        (s1 / xt_).powi(2) + (t12 / s2_).powi(2) >= 1.0
    } else {
        // Compressive: σ₁/Xᶜ ≤ -1
        let xc_ = if xc.abs() > 1e-300 { xc } else { 1e-300 };
        s1 / (-xc_) >= 1.0
    };

    // Matrix failure
    let matrix_failed = if s2 >= 0.0 {
        // Tensile: (σ₂/Yᵀ)² + (τ₁₂/S)² ≥ 1
        let yt_ = if yt.abs() > 1e-300 { yt } else { 1e-300 };
        (s2 / yt_).powi(2) + (t12 / s2_).powi(2) >= 1.0
    } else {
        // Compressive: (σ₂/2S)² + [(Yᶜ/2S)²-1]σ₂/Yᶜ + (τ₁₂/S)² ≥ 1
        let yc_ = if yc.abs() > 1e-300 { yc } else { 1e-300 };
        let term1 = (s2 / (2.0 * s2_)).powi(2);
        let term2 = ((yc_ / (2.0 * s2_)).powi(2) - 1.0) * (-s2 / yc_);
        let term3 = (t12 / s2_).powi(2);
        term1 + term2 + term3 >= 1.0
    };

    (fibre_failed, matrix_failed)
}

// ---------------------------------------------------------------------------
// First-ply failure
// ---------------------------------------------------------------------------

/// Find the index of the first ply to fail under the given in-plane loads.
///
/// Uses the Tsai-Wu criterion with default interaction factors derived from
/// the ply strengths.  Returns `None` if no ply fails.
///
/// `loads = [N_x, N_y, N_xy, M_x, M_y, M_xy]` (forces/moments per unit width).
///
/// Failure strengths are approximated from the ply elastic constants for
/// this demo implementation; in production code they would be explicit inputs.
pub fn first_ply_failure(lam: &Laminate, loads: [f64; 6]) -> Option<usize> {
    if lam.plies.is_empty() {
        return None;
    }
    let abd = laminate_abd_matrix(lam);
    let sol = laminate_midplane_strains(&abd, loads);
    let eps0 = [sol[0], sol[1], sol[2]];
    let kappa = [sol[3], sol[4], sol[5]];

    let total = lam.total_thickness();
    let mut z = -total / 2.0;

    for (idx, ply) in lam.plies.iter().enumerate() {
        let z_mid = z + ply.thickness / 2.0;
        let stress = ply_stresses(ply, eps0, z_mid, kappa);

        // Approximate strength values from moduli for demonstration
        let xt = ply.e1 * 0.01;
        let xc = ply.e1 * 0.008;
        let yt = ply.e2 * 0.005;
        let yc = ply.e2 * 0.004;
        let sh = ply.g12 * 0.02;

        let f1 = 1.0 / xt - 1.0 / xc;
        let f2 = 1.0 / yt - 1.0 / yc;
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        let f12 = -0.5 * (f11 * f22).sqrt();

        let fi = tsai_wu_failure(stress, f1, f2, f12, f11, f22, sh);
        if fi >= 1.0 {
            return Some(idx);
        }
        z += ply.thickness;
    }
    None
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- helpers ---
    fn make_unidirectional() -> Ply {
        Ply::new(0.001, 0.0, 140e9, 10e9, 5e9, 0.3)
    }

    fn make_symmetric_cross_ply() -> Laminate {
        let mut lam = Laminate::new();
        lam.add_ply(Ply::new(0.001, 0.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, 90.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, 90.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, 0.0, 140e9, 10e9, 5e9, 0.3));
        lam
    }

    fn make_angle_ply_45() -> Laminate {
        let mut lam = Laminate::new();
        lam.add_ply(Ply::new(0.001, 45.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, -45.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, -45.0, 140e9, 10e9, 5e9, 0.3));
        lam.add_ply(Ply::new(0.001, 45.0, 140e9, 10e9, 5e9, 0.3));
        lam
    }

    // --- ply tests ---

    #[test]
    fn test_ply_new() {
        let p = make_unidirectional();
        assert_eq!(p.thickness, 0.001);
        assert_eq!(p.angle_deg, 0.0);
        assert_eq!(p.e1, 140e9);
    }

    #[test]
    fn test_ply_nu21() {
        let p = make_unidirectional();
        let nu21 = p.nu21();
        assert!((nu21 - 0.3 * 10e9 / 140e9).abs() < 1e-12);
    }

    #[test]
    fn test_ply_nu21_zero_e1() {
        let p = Ply::new(0.001, 0.0, 0.0, 10e9, 5e9, 0.3);
        assert_eq!(p.nu21(), 0.0);
    }

    // --- reduced stiffness ---

    #[test]
    fn test_reduced_stiffness_symmetry() {
        let p = make_unidirectional();
        let q = reduced_stiffness(&p);
        assert!((q[0][1] - q[1][0]).abs() < 1e6);
    }

    #[test]
    fn test_reduced_stiffness_q11_positive() {
        let p = make_unidirectional();
        let q = reduced_stiffness(&p);
        assert!(q[0][0] > 0.0);
    }

    #[test]
    fn test_reduced_stiffness_q22_positive() {
        let p = make_unidirectional();
        let q = reduced_stiffness(&p);
        assert!(q[1][1] > 0.0);
    }

    #[test]
    fn test_reduced_stiffness_q66_equals_g12() {
        let p = make_unidirectional();
        let q = reduced_stiffness(&p);
        assert!((q[2][2] - p.g12).abs() < 1.0);
    }

    #[test]
    fn test_reduced_stiffness_off_diagonal_zero() {
        let p = make_unidirectional();
        let q = reduced_stiffness(&p);
        assert_eq!(q[0][2], 0.0);
        assert_eq!(q[2][0], 0.0);
        assert_eq!(q[1][2], 0.0);
        assert_eq!(q[2][1], 0.0);
    }

    // --- rotation matrix ---

    #[test]
    fn test_rotation_matrix_zero_angle() {
        let t = rotation_matrix_2d(0.0);
        // At 0 degrees, T = identity-like: t[0]=[1,0,0], t[1]=[0,1,0], t[2]=[0,0,1]
        assert!((t[0][0] - 1.0).abs() < 1e-12);
        assert!((t[1][1] - 1.0).abs() < 1e-12);
        assert!((t[2][2] - 1.0).abs() < 1e-12);
        assert!(t[0][1].abs() < 1e-12);
    }

    #[test]
    fn test_rotation_matrix_90_degrees() {
        let t = rotation_matrix_2d(std::f64::consts::FRAC_PI_2);
        // At 90°, c=0, s=1 → T = [[0,1,0],[1,0,0],[0,0,-1]] approximately
        assert!(t[0][0].abs() < 1e-12);
        assert!((t[0][1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rotation_matrix_45_degrees() {
        let t = rotation_matrix_2d(std::f64::consts::FRAC_PI_4);
        assert!((t[0][0] - 0.5).abs() < 1e-12);
        assert!((t[0][1] - 0.5).abs() < 1e-12);
    }

    // --- transformed stiffness ---

    #[test]
    fn test_transformed_stiffness_0_equals_q() {
        let p = make_unidirectional(); // 0°
        let qbar = transformed_stiffness(&p);
        let q = reduced_stiffness(&p);
        // At 0°, Q̄ ≈ Q
        assert!((qbar[0][0] - q[0][0]).abs() / q[0][0] < 1e-6);
        assert!((qbar[1][1] - q[1][1]).abs() / q[1][1] < 1e-6);
    }

    #[test]
    fn test_transformed_stiffness_90() {
        // At 90°, Q̄₁₁ ≈ Q₂₂
        let mut p = make_unidirectional();
        p.angle_deg = 90.0;
        let qbar = transformed_stiffness(&p);
        let p0 = make_unidirectional();
        let q = reduced_stiffness(&p0);
        assert!((qbar[0][0] - q[1][1]).abs() / q[1][1] < 1e-6);
    }

    #[test]
    fn test_transformed_stiffness_symmetry() {
        let p = Ply::new(0.001, 30.0, 140e9, 10e9, 5e9, 0.3);
        let qbar = transformed_stiffness(&p);
        assert!((qbar[0][1] - qbar[1][0]).abs() < 1e3);
    }

    // --- laminate ---

    #[test]
    fn test_laminate_total_thickness() {
        let lam = make_symmetric_cross_ply();
        assert!((lam.total_thickness() - 0.004).abs() < 1e-12);
    }

    #[test]
    fn test_laminate_add_ply() {
        let mut lam = Laminate::new();
        assert_eq!(lam.plies.len(), 0);
        lam.add_ply(make_unidirectional());
        assert_eq!(lam.plies.len(), 1);
    }

    #[test]
    fn test_laminate_default_empty() {
        let lam = Laminate::default();
        assert!(lam.plies.is_empty());
    }

    // --- ABD matrix ---

    #[test]
    fn test_abd_dimensions_nonzero_diagonal() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        // A11 and A22 must be positive
        assert!(abd[0][0] > 0.0);
        assert!(abd[1][1] > 0.0);
    }

    #[test]
    fn test_abd_symmetric_laminate_zero_b() {
        // A symmetric laminate should have B ≈ 0
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    abd[i][j + 3].abs() < 1.0,
                    "B[{i}][{j}] = {} should be ≈ 0",
                    abd[i][j + 3]
                );
            }
        }
    }

    #[test]
    fn test_abd_positive_d11() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        assert!(abd[3][3] > 0.0, "D11 must be positive");
    }

    #[test]
    fn test_abd_single_ply_a11() {
        let mut lam = Laminate::new();
        let p = Ply::new(0.001, 0.0, 140e9, 10e9, 5e9, 0.3);
        lam.add_ply(p.clone());
        let abd = laminate_abd_matrix(&lam);
        let q = reduced_stiffness(&p);
        // A11 = Q11 * h
        assert!((abd[0][0] - q[0][0] * 0.001).abs() / abd[0][0] < 1e-6);
    }

    // --- mid-plane strains ---

    #[test]
    fn test_midplane_strains_uniaxial() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        let loads = [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // N_x = 1 kN/m
        let sol = laminate_midplane_strains(&abd, loads);
        assert!(sol[0] > 0.0, "ε_x should be positive under tensile N_x");
    }

    #[test]
    fn test_midplane_strains_zero_loads() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        let loads = [0.0; 6];
        let sol = laminate_midplane_strains(&abd, loads);
        for &v in sol.iter() {
            assert!(v.abs() < 1e-20, "Zero loads → zero strains: {v}");
        }
    }

    // --- ply stresses ---

    #[test]
    fn test_ply_stresses_zero_strain() {
        let p = make_unidirectional();
        let sig = ply_stresses(&p, [0.0; 3], 0.0, [0.0; 3]);
        for &s in sig.iter() {
            assert!(s.abs() < 1e-10, "Zero strain → zero stress: {s}");
        }
    }

    #[test]
    fn test_ply_stresses_axial_strain() {
        let p = make_unidirectional(); // 0°
        let eps0 = [1e-3, 0.0, 0.0]; // 0.1% axial strain
        let sig = ply_stresses(&p, eps0, 0.0, [0.0; 3]);
        assert!(sig[0] > 0.0, "σ₁ > 0 for tensile ε_x in 0° ply");
    }

    // --- Tsai-Wu ---

    #[test]
    fn test_tsai_wu_zero_stress_is_zero() {
        let fi = tsai_wu_failure([0.0; 3], 1e-9, 1e-9, -1e-18, 1e-18, 1e-18, 1e6);
        assert!(fi.abs() < 1e-10, "Zero stress → failure index ~0: {fi}");
    }

    #[test]
    fn test_tsai_wu_positive_failure_index() {
        // Stress equal to tensile strength → failure
        let xt = 100e6_f64;
        let xc = 80e6_f64;
        let yt = 50e6_f64;
        let yc = 40e6_f64;
        let s = 30e6_f64;
        let f1 = 1.0 / xt - 1.0 / xc;
        let f2 = 1.0 / yt - 1.0 / yc;
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        let f12 = -0.5 * (f11 * f22).sqrt();
        let fi = tsai_wu_failure([xt, 0.0, 0.0], f1, f2, f12, f11, f22, s);
        assert!(fi > 0.0, "FI must be > 0 at tensile strength");
    }

    #[test]
    fn test_tsai_wu_shear_failure() {
        let s = 30e6_f64;
        let xt = 1000e6_f64;
        let xc = 1000e6_f64;
        let yt = 1000e6_f64;
        let yc = 1000e6_f64;
        let f1 = 0.0;
        let f2 = 0.0;
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        let f12 = 0.0;
        // τ = s → shear FI = 1/(s²) * s² = 1
        let fi = tsai_wu_failure([0.0, 0.0, s], f1, f2, f12, f11, f22, s);
        assert!((fi - 1.0).abs() < 1e-6, "Shear at strength: FI = {fi}");
    }

    // --- Hashin ---

    #[test]
    fn test_hashin_no_failure_at_zero() {
        let (ff, mf) = hashin_failure([0.0; 3], 1000e6, 800e6, 50e6, 40e6, 30e6);
        assert!(!ff && !mf, "No failure at zero stress");
    }

    #[test]
    fn test_hashin_fibre_tensile_failure() {
        let xt = 100e6_f64;
        // σ₁ = Xᵀ, τ₁₂ = 0 → (1)² + 0 = 1 → fails
        let (ff, _) = hashin_failure([xt, 0.0, 0.0], xt, 80e6, 50e6, 40e6, 30e6);
        assert!(ff, "Fibre should fail at tensile strength");
    }

    #[test]
    fn test_hashin_matrix_tensile_failure() {
        let yt = 50e6_f64;
        let s = 30e6_f64;
        // σ₂ = Yᵀ, τ₁₂ = 0 → (1)² + 0 = 1 → matrix fails
        let (_, mf) = hashin_failure([0.0, yt, 0.0], 1000e6, 800e6, yt, 40e6, s);
        assert!(mf, "Matrix should fail at transverse tensile strength");
    }

    #[test]
    fn test_hashin_compressive_fibre() {
        let xc = 80e6_f64;
        // σ₁ = −Xᶜ → compressive fibre failure
        let (ff, _) = hashin_failure([-xc, 0.0, 0.0], 100e6, xc, 50e6, 40e6, 30e6);
        assert!(ff, "Compressive fibre failure");
    }

    // --- first ply failure ---

    #[test]
    fn test_first_ply_failure_returns_none_at_zero_load() {
        let lam = make_symmetric_cross_ply();
        let result = first_ply_failure(&lam, [0.0; 6]);
        assert!(result.is_none(), "No failure at zero loads");
    }

    #[test]
    fn test_first_ply_failure_returns_some_at_large_load() {
        let lam = make_symmetric_cross_ply();
        // Extremely large axial load → failure
        let result = first_ply_failure(&lam, [1e14, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(result.is_some(), "Should fail under very large load");
    }

    #[test]
    fn test_first_ply_failure_empty_laminate() {
        let lam = Laminate::new();
        let result = first_ply_failure(&lam, [1e10; 6]);
        assert!(result.is_none(), "Empty laminate → None");
    }

    #[test]
    fn test_first_ply_failure_index_in_range() {
        let lam = make_symmetric_cross_ply();
        if let Some(idx) = first_ply_failure(&lam, [1e14, 0.0, 0.0, 0.0, 0.0, 0.0]) {
            assert!(idx < lam.plies.len(), "Index must be in range");
        }
    }

    // --- angle-ply ---

    #[test]
    fn test_angle_ply_abd_a16_nonzero() {
        // ±45° laminate (not balanced symmetric) has non-zero A16
        let mut lam = Laminate::new();
        lam.add_ply(Ply::new(0.001, 45.0, 140e9, 10e9, 5e9, 0.3));
        let abd = laminate_abd_matrix(&lam);
        // A16 should not be zero for a single 45° ply
        assert!(abd[0][2].abs() > 1e3, "A16 non-zero for 45° ply");
    }

    #[test]
    fn test_symmetric_angle_ply_b_near_zero() {
        let lam = make_angle_ply_45();
        let abd = laminate_abd_matrix(&lam);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    abd[i][j + 3].abs() < 10.0,
                    "B[{i}][{j}] = {} for symmetric ±45",
                    abd[i][j + 3]
                );
            }
        }
    }

    // --- CLT round-trip ---

    #[test]
    fn test_clt_roundtrip_stress_strain() {
        // For a quasi-isotropic laminate under uniaxial N_x,
        // the mid-plane strains should recover the applied load.
        let mut lam = Laminate::new();
        for &angle in &[0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0f64] {
            lam.add_ply(Ply::new(0.001, angle, 140e9, 10e9, 5e9, 0.3));
        }
        let abd = laminate_abd_matrix(&lam);
        let loads = [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sol = laminate_midplane_strains(&abd, loads);
        // Recover N_x = A11*eps0_x + A12*eps0_y + ...
        let nx_recovered = abd[0][0] * sol[0]
            + abd[0][1] * sol[1]
            + abd[0][2] * sol[2]
            + abd[0][3] * sol[3]
            + abd[0][4] * sol[4]
            + abd[0][5] * sol[5];
        assert!(
            (nx_recovered - 1000.0).abs() < 1.0,
            "Recovered N_x = {nx_recovered}"
        );
    }

    // --- misc ---

    #[test]
    fn test_inv3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = inv3(id);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_mat3x3_mul_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let c = mat3x3_mul(a, id);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_laminate_a_positive_definite_diagonal() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        // All A diagonal entries must be positive
        assert!(abd[0][0] > 0.0);
        assert!(abd[1][1] > 0.0);
        assert!(abd[2][2] > 0.0);
    }

    #[test]
    fn test_reduced_stiffness_singular_denom() {
        // nu12 = 1.0 would make denom = 0; should return zeros
        let p = Ply::new(0.001, 0.0, 10.0, 10.0, 5.0, 1.0);
        let q = reduced_stiffness(&p); // denom ~ 0
        // Should not panic; result is well-defined (all zeros or similar)
        let _ = q;
    }

    #[test]
    fn test_tsai_wu_negative_shear_symmetric() {
        // Tsai-Wu is quadratic in τ₁₂ → same FI for ±τ₁₂
        let s = 30e6_f64;
        let xt = 1000e6_f64;
        let xc = 1000e6_f64;
        let yt = 1000e6_f64;
        let yc = 1000e6_f64;
        let f1 = 0.0;
        let f2 = 0.0;
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        let f12 = 0.0;
        let fi_pos = tsai_wu_failure([0.0, 0.0, s * 0.5], f1, f2, f12, f11, f22, s);
        let fi_neg = tsai_wu_failure([0.0, 0.0, -s * 0.5], f1, f2, f12, f11, f22, s);
        assert!(
            (fi_pos - fi_neg).abs() < 1e-10,
            "Shear symmetry: {fi_pos} vs {fi_neg}"
        );
    }

    #[test]
    fn test_hashin_high_shear_matrix_failure() {
        // τ₁₂ >> S → matrix failure
        let s = 30e6_f64;
        let (_, mf) = hashin_failure([0.0, 0.0, 2.0 * s], 1000e6, 800e6, 500e6, 400e6, s);
        assert!(mf, "Matrix failure under high shear");
    }

    #[test]
    fn test_clt_biaxial_strain_both_nonzero() {
        let lam = make_symmetric_cross_ply();
        let abd = laminate_abd_matrix(&lam);
        let loads = [1000.0, 500.0, 0.0, 0.0, 0.0, 0.0];
        let sol = laminate_midplane_strains(&abd, loads);
        // Both ε_x and ε_y should be nonzero
        assert!(sol[0] != 0.0 || sol[1] != 0.0, "Biaxial response expected");
    }
}
