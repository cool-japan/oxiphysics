// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tsai-Wu ply failure analysis of [`CompositeLaminateClt`](super::CompositeLaminateClt)
//! laminates by classical lamination theory (CLT).
//!
//! For applied membrane resultants `N = [Nx, Ny, Nxy]` (N/mm) and zero applied
//! moments the mid-plane strains `ε⁰` and curvatures `κ` solve the
//! extension-bending equations
//!
//! ```text
//! [ A  B ] [ ε⁰ ]   [ N ]
//! [ B  D ] [ κ  ] = [ 0 ] ,
//!
//! (A, B, D)ᵢⱼ = Σₖ Q̄ᵢⱼ⁽ᵏ⁾ · (tₖ, (zₖ² − zₖ₋₁²)/2, (zₖ³ − zₖ₋₁³)/3)
//! ```
//!
//! with the transformed reduced stiffness `Q̄⁽ᵏ⁾` of every ply expressed in
//! MPa (the moduli are given in GPa) so that `A` has the unit N/mm.  For a
//! symmetric laminate `B = 0` and this reduces to `ε⁰ = A⁻¹ N`, `κ = 0`.  The
//! strain varies linearly through the thickness, `ε(z) = ε⁰ + z κ`, and the
//! stress of ply `k` in laminate axes is `σ = Q̄⁽ᵏ⁾ ε(z)`.  It is rotated into
//! the material axes of the ply before the Tsai-Wu criterion is evaluated.
//! Because the Tsai-Wu index is a convex function of the stress, its maximum
//! over a ply is attained on one of the two ply surfaces.

use crate::composite::solve_linear_6;

/// One ply as stored in `CompositeLaminateClt::plies`:
/// `(angle_deg, thickness_mm, e1_gpa, e2_gpa, g12_gpa, nu12)`.
pub(super) type Ply = (f64, f64, f64, f64, f64, f64);

/// Conversion of a stiffness in GPa to MPa (so that `stiffness · mm` is N/mm).
const GPA_TO_MPA: f64 = 1.0e3;

/// Why the ply stresses of a laminate could not be determined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum StressError {
    /// A ply, a load component or a strength is not a finite physical value.
    InvalidInput,
    /// The laminate cannot carry in-plane load: the extension-bending matrix
    /// is singular (no plies, or no ply thickness).
    Singular,
}

/// Stress state of one ply at both of its surfaces.
#[derive(Debug, Clone, Copy)]
pub(super) struct PlyStress {
    /// Cosine of the fibre angle.
    pub(super) cos: f64,
    /// Sine of the fibre angle.
    pub(super) sin: f64,
    /// Distance of the lower ply surface from the laminate mid-plane (mm).
    pub(super) z_bottom: f64,
    /// Distance of the upper ply surface from the laminate mid-plane (mm).
    pub(super) z_top: f64,
    /// `[σx, σy, τxy]` in laminate axes at the lower surface (MPa).
    pub(super) bottom: [f64; 3],
    /// `[σx, σy, τxy]` in laminate axes at the upper surface (MPa).
    pub(super) top: [f64; 3],
}

impl PlyStress {
    /// Rotate a laminate-axes stress `[σx, σy, τxy]` into the material axes
    /// of the ply, returning `[σ1, σ2, τ12]`.
    pub(super) fn to_material_axes(self, sigma: [f64; 3]) -> [f64; 3] {
        let (c, s) = (self.cos, self.sin);
        let [sx, sy, txy] = sigma;
        [
            c * c * sx + s * s * sy + 2.0 * c * s * txy,
            s * s * sx + c * c * sy - 2.0 * c * s * txy,
            -c * s * sx + c * s * sy + (c * c - s * s) * txy,
        ]
    }
}

/// Transformed reduced stiffness `Q̄` (3×3, MPa) of a ply in laminate axes,
/// or `None` when the ply data are not physical (non-finite values,
/// negative thickness, non-positive moduli, or `1 − ν12 ν21 <= 0`).
fn transformed_stiffness_mpa(ply: &Ply) -> Option<[[f64; 3]; 3]> {
    let (angle_deg, thickness, e1, e2, g12, nu12) = *ply;
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !(angle_deg.is_finite()
        && thickness.is_finite()
        && thickness >= 0.0
        && positive(e1)
        && positive(e2)
        && positive(g12)
        && nu12.is_finite())
    {
        return None;
    }
    let nu21 = nu12 * e2 / e1;
    let denominator = 1.0 - nu12 * nu21;
    if denominator.is_nan() || denominator <= 0.0 {
        return None;
    }
    let q11 = GPA_TO_MPA * e1 / denominator;
    let q22 = GPA_TO_MPA * e2 / denominator;
    let q12 = GPA_TO_MPA * nu12 * e2 / denominator;
    let q66 = GPA_TO_MPA * g12;

    let (s, c) = angle_deg.to_radians().sin_cos();
    let (c2, s2, cs) = (c * c, s * s, c * s);
    let (c4, s4) = (c2 * c2, s2 * s2);
    let q11b = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4;
    let q22b = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4;
    let q12b = (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (c4 + s4);
    let q66b = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (c4 + s4);
    // Q̄16 = (Q11 − Q12 − 2Q66) c³s − (Q22 − Q12 − 2Q66) c s³ and
    // Q̄26 = (Q11 − Q12 − 2Q66) c s³ − (Q22 − Q12 − 2Q66) c³s, with
    // c³s = c²·(cs) and c s³ = s²·(cs).
    let q16b = (q11 - q12 - 2.0 * q66) * c2 * cs - (q22 - q12 - 2.0 * q66) * s2 * cs;
    let q26b = (q11 - q12 - 2.0 * q66) * s2 * cs - (q22 - q12 - 2.0 * q66) * c2 * cs;
    Some([[q11b, q12b, q16b], [q12b, q22b, q26b], [q16b, q26b, q66b]])
}

/// Matrix–vector product of a 3×3 matrix with a 3-vector.
fn mat3_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Ply stresses of a laminate loaded by the in-plane resultants `resultants`
/// (`[Nx, Ny, Nxy]`, N/mm) with zero applied moments, see the module
/// documentation.
///
/// Returns one [`PlyStress`] per ply, bottom to top.
///
/// # Errors
///
/// [`StressError::InvalidInput`] for non-physical ply data or non-finite
/// resultants and [`StressError::Singular`] when the laminate has no
/// load-carrying thickness or its extension-bending matrix is singular.
pub(super) fn ply_stresses(
    plies: &[Ply],
    resultants: [f64; 3],
) -> Result<Vec<PlyStress>, StressError> {
    if resultants.iter().any(|n| !n.is_finite()) {
        return Err(StressError::InvalidInput);
    }
    let stiffness = plies
        .iter()
        .map(transformed_stiffness_mpa)
        .collect::<Option<Vec<_>>>()
        .ok_or(StressError::InvalidInput)?;
    let total_thickness: f64 = plies.iter().map(|ply| ply.1).sum();
    if !total_thickness.is_finite() || total_thickness <= 0.0 {
        return Err(StressError::Singular);
    }

    // Ply interfaces measured from the geometric mid-plane.
    let mut interfaces = Vec::with_capacity(plies.len() + 1);
    let mut z = -0.5 * total_thickness;
    interfaces.push(z);
    for ply in plies {
        z += ply.1;
        interfaces.push(z);
    }

    // [[A, B/h], [B/h, D/h²]] · [ε⁰, h·κ] = [N, 0]: scaling the coupling and
    // bending blocks by the laminate thickness h gives all entries the unit
    // N/mm, so the pivot test of the solver is meaningful for any thickness.
    let h = total_thickness;
    let mut system = [[0.0_f64; 6]; 6];
    for (q, bounds) in stiffness.iter().zip(interfaces.windows(2)) {
        let (z0, z1) = (bounds[0], bounds[1]);
        let dz = z1 - z0;
        // (z1² − z0²)/2 and (z1³ − z0³)/3, factored to avoid cancellation.
        let first_moment = 0.5 * dz * (z1 + z0);
        let second_moment = dz * (z1 * z1 + z1 * z0 + z0 * z0) / 3.0;
        for (i, q_row) in q.iter().enumerate() {
            for (j, &q_ij) in q_row.iter().enumerate() {
                system[i][j] += q_ij * dz;
                system[i][j + 3] += q_ij * first_moment / h;
                system[i + 3][j] += q_ij * first_moment / h;
                system[i + 3][j + 3] += q_ij * second_moment / (h * h);
            }
        }
    }
    let rhs = [resultants[0], resultants[1], resultants[2], 0.0, 0.0, 0.0];
    let solution = solve_linear_6(system, rhs).ok_or(StressError::Singular)?;
    let mid_strain = [solution[0], solution[1], solution[2]];
    let curvature = [solution[3] / h, solution[4] / h, solution[5] / h];

    let stress_at = |q: &[[f64; 3]; 3], z: f64| {
        let strain = [
            mid_strain[0] + z * curvature[0],
            mid_strain[1] + z * curvature[1],
            mid_strain[2] + z * curvature[2],
        ];
        mat3_vec(q, strain)
    };
    let mut stresses = Vec::with_capacity(plies.len());
    for ((ply, q), bounds) in plies
        .iter()
        .zip(stiffness.iter())
        .zip(interfaces.windows(2))
    {
        let (sin, cos) = ply.0.to_radians().sin_cos();
        stresses.push(PlyStress {
            cos,
            sin,
            z_bottom: bounds[0],
            z_top: bounds[1],
            bottom: stress_at(q, bounds[0]),
            top: stress_at(q, bounds[1]),
        });
    }
    Ok(stresses)
}

/// Tsai-Wu strength coefficients of a ply.
#[derive(Debug, Clone, Copy)]
struct TsaiWuCoefficients {
    f1: f64,
    f2: f64,
    f11: f64,
    f22: f64,
    f66: f64,
    f12: f64,
}

impl TsaiWuCoefficients {
    /// Coefficients from the positive strength magnitudes `xt`, `xc`
    /// (longitudinal tension/compression), `yt`, `yc` (transverse) and `s`
    /// (in-plane shear), with the interaction term
    /// `F12 = −½·√(F11·F22)`.  `None` unless all strengths are positive and
    /// finite.
    fn new(xt: f64, xc: f64, yt: f64, yc: f64, s: f64) -> Option<Self> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !(positive(xt) && positive(xc) && positive(yt) && positive(yc) && positive(s)) {
            return None;
        }
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        Some(Self {
            f1: 1.0 / xt - 1.0 / xc,
            f2: 1.0 / yt - 1.0 / yc,
            f11,
            f22,
            f66: 1.0 / (s * s),
            f12: -0.5 * (f11 * f22).sqrt(),
        })
    }

    /// Tsai-Wu index for the material-axes stress `[σ1, σ2, τ12]`:
    /// `F1σ1 + F2σ2 + F11σ1² + F22σ2² + F66τ12² + 2F12σ1σ2`.
    fn index(&self, material_stress: [f64; 3]) -> f64 {
        let [s1, s2, t12] = material_stress;
        self.f1 * s1
            + self.f2 * s2
            + self.f11 * s1 * s1
            + self.f22 * s2 * s2
            + self.f66 * t12 * t12
            + 2.0 * self.f12 * s1 * s2
    }
}

/// Maximum Tsai-Wu index over all plies (both surfaces of each) of the
/// laminate `plies` under the resultants `nx`, `ny`, `nxy` (N/mm).
///
/// Returns `f64::INFINITY` when the laminate cannot carry in-plane load
/// (no plies, zero thickness or singular extension-bending matrix) and
/// `f64::NAN` for non-physical ply data, non-finite loads or strengths that
/// are not positive and finite.
pub(super) fn tsai_wu_index(
    plies: &[Ply],
    nx: f64,
    ny: f64,
    nxy: f64,
    xt: f64,
    xc: f64,
    yt: f64,
    yc: f64,
    s: f64,
) -> f64 {
    let Some(coefficients) = TsaiWuCoefficients::new(xt, xc, yt, yc, s) else {
        return f64::NAN;
    };
    let stresses = match ply_stresses(plies, [nx, ny, nxy]) {
        Ok(stresses) => stresses,
        Err(StressError::InvalidInput) => return f64::NAN,
        Err(StressError::Singular) => return f64::INFINITY,
    };
    stresses
        .iter()
        .filter(|ply| ply.z_top > ply.z_bottom)
        .flat_map(|ply| {
            [ply.bottom, ply.top]
                .into_iter()
                .map(|sigma| coefficients.index(ply.to_material_axes(sigma)))
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aerospace::CompositeLaminateClt;

    // IM7/8552 ply properties (GPa) and strengths (MPa).
    const E1: f64 = 161.0;
    const E2: f64 = 11.4;
    const G12: f64 = 5.17;
    const NU12: f64 = 0.32;
    const XT: f64 = 2560.0;
    const XC: f64 = 1590.0;
    const YT: f64 = 73.0;
    const YC: f64 = 185.0;
    const S: f64 = 90.0;
    /// Ply thickness (mm).
    const T: f64 = 0.125;

    fn ply(angle: f64) -> Ply {
        (angle, T, E1, E2, G12, NU12)
    }

    fn laminate(angles: &[f64]) -> CompositeLaminateClt {
        CompositeLaminateClt::new(angles.iter().map(|&angle| ply(angle)).collect())
    }

    fn index(lam: &CompositeLaminateClt, n: [f64; 3]) -> f64 {
        lam.tsai_wu_index(n[0], n[1], n[2], XT, XC, YT, YC, S)
    }

    fn relative_error(value: f64, reference: f64) -> f64 {
        ((value - reference) / reference).abs()
    }

    /// Tsai-Wu index of a material-axes stress state, written out
    /// independently of [`TsaiWuCoefficients`].
    fn reference_index(s1: f64, s2: f64, t12: f64) -> f64 {
        let f11 = 1.0 / (XT * XC);
        let f22 = 1.0 / (YT * YC);
        let f12 = -0.5 * (f11 * f22).sqrt();
        (1.0 / XT - 1.0 / XC) * s1
            + (1.0 / YT - 1.0 / YC) * s2
            + f11 * s1 * s1
            + f22 * s2 * s2
            + t12 * t12 / (S * S)
            + 2.0 * f12 * s1 * s2
    }

    /// `Q̄ = T⁻¹ · Q · T⁻ᵀ` with the explicit stress-transformation matrix
    /// (engineering shear strain), in MPa.
    fn tensor_rotated_stiffness(angle_deg: f64) -> [[f64; 3]; 3] {
        let nu21 = NU12 * E2 / E1;
        let den = 1.0 - NU12 * nu21;
        let q = [
            [1.0e3 * E1 / den, 1.0e3 * NU12 * E2 / den, 0.0],
            [1.0e3 * NU12 * E2 / den, 1.0e3 * E2 / den, 0.0],
            [0.0, 0.0, 1.0e3 * G12],
        ];
        let (s, c) = angle_deg.to_radians().sin_cos();
        let t_inv = [
            [c * c, s * s, -2.0 * c * s],
            [s * s, c * c, 2.0 * c * s],
            [c * s, -c * s, c * c - s * s],
        ];
        let mut tmp = [[0.0_f64; 3]; 3];
        let mut out = [[0.0_f64; 3]; 3];
        for (row, tmp_row) in t_inv.iter().zip(tmp.iter_mut()) {
            for (j, entry) in tmp_row.iter_mut().enumerate() {
                *entry = row
                    .iter()
                    .zip(q.iter())
                    .map(|(a, q_row)| a * q_row[j])
                    .sum();
            }
        }
        for (tmp_row, out_row) in tmp.iter().zip(out.iter_mut()) {
            for (entry, t_row) in out_row.iter_mut().zip(t_inv.iter()) {
                *entry = tmp_row.iter().zip(t_row.iter()).map(|(a, b)| a * b).sum();
            }
        }
        out
    }

    #[test]
    fn transformed_stiffness_matches_the_tensor_rotation() {
        for &angle in &[0.0, 15.0, 30.0, 45.0, -45.0, 60.0, 90.0, 123.0, -170.0] {
            let got = transformed_stiffness_mpa(&ply(angle)).unwrap();
            let want = tensor_rotated_stiffness(angle);
            for (got_row, want_row) in got.iter().zip(want.iter()) {
                for (g, w) in got_row.iter().zip(want_row.iter()) {
                    assert!(
                        (g - w).abs() <= 1.0e-12 * w.abs().max(1.0e3),
                        "angle {angle}: {got:?} vs {want:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn single_zero_degree_ply_reaches_one_at_each_strength() {
        let lam = laminate(&[0.0]);
        let cases: [([f64; 3], &str); 7] = [
            ([XT * T, 0.0, 0.0], "longitudinal tension"),
            ([-XC * T, 0.0, 0.0], "longitudinal compression"),
            ([0.0, YT * T, 0.0], "transverse tension"),
            ([0.0, -YC * T, 0.0], "transverse compression"),
            ([0.0, 0.0, S * T], "positive shear"),
            ([0.0, 0.0, -S * T], "negative shear"),
            ([0.0, 0.0, 0.0], "no load"),
        ];
        for (n, label) in cases {
            let fi = index(&lam, n);
            let expected = if n == [0.0, 0.0, 0.0] { 0.0 } else { 1.0 };
            assert!(
                (fi - expected).abs() < 1.0e-12,
                "{label}: FI = {fi}, expected {expected}"
            );
        }
    }

    #[test]
    fn ninety_degree_ply_swaps_the_laminate_axes() {
        let lam = laminate(&[90.0]);
        let along_x = index(&lam, [YT * T, 0.0, 0.0]);
        let along_y = index(&lam, [0.0, XT * T, 0.0]);
        assert!((along_x - 1.0).abs() < 1.0e-12, "FI = {along_x}");
        assert!((along_y - 1.0).abs() < 1.0e-12, "FI = {along_y}");
    }

    #[test]
    fn plies_carry_the_applied_resultant_over_the_thickness() {
        // Regression: the stress used to be Q̄·N/t, i.e. ~4e5 "MPa" for a
        // 0° ply loaded to its tensile strength.
        let n = [XT * T, 0.0, 0.0];
        let stresses = ply_stresses(&[ply(0.0)], n).unwrap();
        assert_eq!(stresses.len(), 1);
        for surface in [stresses[0].bottom, stresses[0].top] {
            assert!(relative_error(surface[0], XT) < 1.0e-12, "{surface:?}");
            assert!(surface[1].abs() < 1.0e-9 && surface[2].abs() < 1.0e-9);
        }
        // For any single ply, whatever its stiffness, σ = N/t.
        let n = [310.0, -45.0, 62.0];
        for &angle in &[0.0, 17.0, 45.0, 90.0, -60.0] {
            let stresses = ply_stresses(&[ply(angle)], n).unwrap();
            for surface in [stresses[0].bottom, stresses[0].top] {
                for (s, applied) in surface.iter().zip(n.iter()) {
                    assert!(
                        (s - applied / T).abs() < 1.0e-9 * (applied / T).abs().max(1.0),
                        "angle {angle}: {surface:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn off_axis_single_ply_matches_the_rotated_applied_stress() {
        // σ = N/t in laminate axes, rotated into the material axes by hand.
        let n = [310.0, -45.0, 62.0];
        for &angle in &[15.0_f64, 30.0, 45.0, 60.0, -30.0, 120.0] {
            let (s, c) = angle.to_radians().sin_cos();
            let (sx, sy, txy) = (n[0] / T, n[1] / T, n[2] / T);
            let s1 = c * c * sx + s * s * sy + 2.0 * c * s * txy;
            let s2 = s * s * sx + c * c * sy - 2.0 * c * s * txy;
            let t12 = -c * s * sx + c * s * sy + (c * c - s * s) * txy;
            let expected = reference_index(s1, s2, t12);
            let got = index(&laminate(&[angle]), n);
            assert!(
                relative_error(got, expected) < 1.0e-10,
                "angle {angle}: {got} vs {expected}"
            );
        }
        // Reference values of an independent implementation.
        assert!(
            relative_error(
                index(&laminate(&[45.0]), [100.0, 0.0, 0.0]),
                34.17940229695736
            ) < 1.0e-10
        );
        assert!(relative_error(index(&laminate(&[30.0]), n), 120.21275618556996) < 1.0e-10);
        assert!(relative_error(index(&laminate(&[-30.0]), n), 316.72833944419136) < 1.0e-10);
        assert!(relative_error(index(&laminate(&[60.0]), n), 409.27146455378954) < 1.0e-10);
    }

    #[test]
    fn symmetric_cross_ply_matches_the_closed_form() {
        // [0/90/90/0]: A11 = A22 = (Q11 + Q22)/2 · h, A12 = Q12 · h, A16 = A26 = 0.
        let nu21 = NU12 * E2 / E1;
        let den = 1.0 - NU12 * nu21;
        let (q11, q22, q12) = (1.0e3 * E1 / den, 1.0e3 * E2 / den, 1.0e3 * NU12 * E2 / den);
        let h = 4.0 * T;
        let (a11, a12) = (0.5 * (q11 + q22) * h, q12 * h);
        let det = a11 * a11 - a12 * a12;
        let nx = 60.0;
        let (ex, ey) = (a11 * nx / det, -a12 * nx / det);
        // 0° plies: material axes coincide with the laminate axes.
        let fi_0 = reference_index(q11 * ex + q12 * ey, q12 * ex + q22 * ey, 0.0);
        // 90° plies: material axis 1 points along y.
        let fi_90 = reference_index(q12 * ex + q11 * ey, q22 * ex + q12 * ey, 0.0);
        let expected = fi_0.max(fi_90);
        let got = index(&laminate(&[0.0, 90.0, 90.0, 0.0]), [nx, 0.0, 0.0]);
        assert!(
            relative_error(got, expected) < 1.0e-10,
            "{got} vs closed form {expected}"
        );
        assert!(relative_error(got, 0.14962957252480227) < 1.0e-10);
        // The transversely loaded 90° plies are the critical ones here.
        assert!(fi_90 > fi_0);
    }

    #[test]
    fn symmetric_laminates_have_a_uniform_ply_stress() {
        let angles = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0];
        let stresses = ply_stresses(
            &angles.iter().map(|&a| ply(a)).collect::<Vec<_>>(),
            [150.0, 40.0, 30.0],
        )
        .unwrap();
        for ply_stress in &stresses {
            for (b, t) in ply_stress.bottom.iter().zip(ply_stress.top.iter()) {
                assert!((b - t).abs() < 1.0e-9 * b.abs().max(1.0), "{ply_stress:?}");
            }
        }
        let fi = index(&laminate(&angles), [150.0, 40.0, 30.0]);
        assert!(
            relative_error(fi, 0.2645395325342791) < 1.0e-10,
            "FI = {fi}"
        );
    }

    #[test]
    fn unsymmetric_laminates_include_the_bending_extension_coupling() {
        // [0/90] curves under pure membrane load; the outer surface of the 90°
        // ply carries far more stress than the mid-plane strain suggests
        // (a mid-plane-only estimate would give FI = 0.336 here).
        let fi = index(&laminate(&[0.0, 90.0]), [60.0, 0.0, 0.0]);
        assert!(
            relative_error(fi, 3.4583197396114618) < 1.0e-10,
            "FI = {fi}"
        );
        let fi = index(&laminate(&[0.0, 45.0, 90.0]), [80.0, 20.0, -15.0]);
        assert!(relative_error(fi, 6.160961270642342) < 1.0e-10, "FI = {fi}");
        let fi = index(&laminate(&[45.0, 90.0, 0.0]), [80.0, 20.0, -15.0]);
        assert!(relative_error(fi, 5.498682816492168) < 1.0e-10, "FI = {fi}");
        let fi = index(&laminate(&[30.0, -45.0, 10.0, 60.0]), [120.0, -35.0, 18.0]);
        assert!(relative_error(fi, 1.165245759068425) < 1.0e-10, "FI = {fi}");
        let fi = index(&laminate(&[0.0, 90.0, 0.0, 90.0]), [60.0, 0.0, 0.0]);
        assert!(
            relative_error(fi, 0.31997435586577405) < 1.0e-10,
            "FI = {fi}"
        );
    }

    #[test]
    fn ply_stresses_satisfy_force_and_moment_equilibrium() {
        let n = [120.0, -35.0, 18.0];
        for angles in [
            vec![0.0, 90.0],
            vec![0.0, 45.0, 90.0],
            vec![30.0, -45.0, 10.0, 60.0],
            vec![0.0, 90.0, 90.0, 0.0],
        ] {
            let plies: Vec<Ply> = angles.iter().map(|&a| ply(a)).collect();
            let stresses = ply_stresses(&plies, n).unwrap();
            let mut force = [0.0_f64; 3];
            let mut moment = [0.0_f64; 3];
            for st in &stresses {
                let dz = st.z_top - st.z_bottom;
                for (k, (b, t)) in st.bottom.iter().zip(st.top.iter()).enumerate() {
                    // Stress is linear in z across the ply: integrate exactly.
                    force[k] += 0.5 * (b + t) * dz;
                    moment[k] += dz / 6.0
                        * (b * (2.0 * st.z_bottom + st.z_top) + t * (st.z_bottom + 2.0 * st.z_top));
                }
            }
            for (k, ((f, m), applied)) in force.iter().zip(moment.iter()).zip(n.iter()).enumerate()
            {
                assert!(
                    (f - applied).abs() < 1.0e-9 * applied.abs(),
                    "{angles:?}: N[{k}] = {f} vs {applied}"
                );
                assert!(
                    m.abs() < 1.0e-9 * applied.abs() * T,
                    "{angles:?}: M[{k}] = {m}"
                );
            }
        }
    }

    #[test]
    fn mirrored_stacking_sequence_gives_the_same_index() {
        let n = [80.0, 20.0, -15.0];
        let forward = index(&laminate(&[0.0, 45.0, 90.0]), n);
        let mirrored = index(&laminate(&[90.0, 45.0, 0.0]), n);
        assert!(
            relative_error(forward, mirrored) < 1.0e-10,
            "{forward} vs {mirrored}"
        );
        // A different (non-mirrored) stacking order of an unsymmetric
        // laminate is a different laminate.
        let shifted = index(&laminate(&[45.0, 90.0, 0.0]), n);
        assert!((forward - shifted).abs() > 1.0e-2 * forward);
    }

    #[test]
    fn index_is_linear_plus_quadratic_in_the_load() {
        // For one ply the stress is linear in the load, hence
        // FI(λN) = λ·a + λ²·b, with a and b from the odd and even parts of
        // FI(±N).
        let single = laminate(&[35.0]);
        let n = [90.0, -25.0, 40.0];
        let minus_n = [-n[0], -n[1], -n[2]];
        let (plus, minus) = (index(&single, n), index(&single, minus_n));
        let linear = 0.5 * (plus - minus);
        let quadratic = 0.5 * (plus + minus);
        for &lambda in &[0.5, 2.0, 3.5] {
            let scaled = [lambda * n[0], lambda * n[1], lambda * n[2]];
            let got = index(&single, scaled);
            let want = lambda * linear + lambda * lambda * quadratic;
            assert!(
                (got - want).abs() < 1.0e-9 * want.abs().max(1.0),
                "lambda {lambda}: {got} vs {want}"
            );
        }
    }
    #[test]
    fn laminate_without_load_carrying_plies_is_singular() {
        assert_eq!(
            index(&CompositeLaminateClt::new(Vec::new()), [1.0, 0.0, 0.0]),
            f64::INFINITY
        );
        let massless = CompositeLaminateClt::new(vec![(0.0, 0.0, E1, E2, G12, NU12)]);
        assert_eq!(index(&massless, [1.0, 0.0, 0.0]), f64::INFINITY);
        assert_eq!(
            ply_stresses(&[], [1.0, 0.0, 0.0]).unwrap_err(),
            StressError::Singular
        );
    }

    #[test]
    fn non_physical_data_gives_nan() {
        let n = [10.0, 0.0, 0.0];
        let bad_modulus = CompositeLaminateClt::new(vec![(0.0, T, 0.0, E2, G12, NU12)]);
        assert!(index(&bad_modulus, n).is_nan());
        let bad_poisson = CompositeLaminateClt::new(vec![(0.0, T, E1, E2, G12, 30.0)]);
        assert!(index(&bad_poisson, n).is_nan());
        let bad_thickness = CompositeLaminateClt::new(vec![(0.0, -T, E1, E2, G12, NU12)]);
        assert!(index(&bad_thickness, n).is_nan());
        let lam = laminate(&[0.0]);
        assert!(index(&lam, [f64::NAN, 0.0, 0.0]).is_nan());
        assert!(
            lam.tsai_wu_index(1.0, 0.0, 0.0, 0.0, XC, YT, YC, S)
                .is_nan()
        );
        assert!(
            lam.tsai_wu_index(1.0, 0.0, 0.0, XT, XC, YT, YC, -S)
                .is_nan()
        );
        assert!(
            lam.tsai_wu_index(1.0, 0.0, 0.0, XT, f64::INFINITY, YT, YC, S)
                .is_nan()
        );
    }
}
