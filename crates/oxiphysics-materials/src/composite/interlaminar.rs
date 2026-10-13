// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Interlaminar shear stresses of laminates from through-thickness
//! equilibrium, with the small dense solver they need.

use super::Laminate;

/// Solve the 6×6 linear system `m · x = rhs` by Gaussian elimination with
/// partial pivoting.
///
/// Returns `None` when the matrix is numerically singular, i.e. a pivot is
/// smaller than `1e-13` times the largest entry of `m`.
pub(crate) fn solve_linear_6(mut m: [[f64; 6]; 6], mut rhs: [f64; 6]) -> Option<[f64; 6]> {
    let scale = m.iter().flatten().fold(0.0_f64, |acc, v| acc.max(v.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let tolerance = 1.0e-13 * scale;

    for col in 0..6 {
        let pivot_row = (col..6)
            .max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))
            .unwrap_or(col);
        let pivot_abs = m[pivot_row][col].abs();
        if pivot_abs.is_nan() || pivot_abs <= tolerance {
            return None;
        }
        m.swap(col, pivot_row);
        rhs.swap(col, pivot_row);

        let pivot_values = m[col];
        let pivot_rhs = rhs[col];
        let pivot = pivot_values[col];
        for (row_values, row_rhs) in m.iter_mut().zip(rhs.iter_mut()).skip(col + 1) {
            let factor = row_values[col] / pivot;
            if factor != 0.0 {
                for (target, source) in row_values.iter_mut().zip(pivot_values.iter()).skip(col) {
                    *target -= factor * source;
                }
                *row_rhs -= factor * pivot_rhs;
            }
        }
    }

    let mut x = [0.0_f64; 6];
    for row in (0..6).rev() {
        let tail: f64 = m[row]
            .iter()
            .zip(x.iter())
            .skip(row + 1)
            .map(|(coefficient, solved)| coefficient * solved)
            .sum();
        x[row] = (rhs[row] - tail) / m[row][row];
    }
    Some(x)
}

/// Dot product of a row of the transformed stiffness `Q̄` with a 3-vector.
fn qbar_row_dot(qbar: &[[f64; 3]; 3], row: usize, v: &[f64]) -> f64 {
    qbar[row][0] * v[0] + qbar[row][1] * v[1] + qbar[row][2] * v[2]
}

impl Laminate {
    /// Gradient of the generalised strains `[ε⁰, κ]` in one in-plane direction
    /// for a transverse shear force `shear` carried by the bending moment
    /// component `moment_row` (0 = Mₓ, 1 = M_y) while the membrane resultants
    /// stay constant.
    ///
    /// Solves the extension-bending equations
    /// `[[A, B], [B, D]] · [∂ε⁰, ∂κ] = [∂N, ∂M] = [0, shear·e_row]`.
    /// The system is non-dimensionalised with the laminate thickness `h`
    /// (unknowns `[∂ε⁰, h·∂κ]`, moment rows divided by `h`) so that the
    /// pivot test is meaningful for any thickness.
    fn strain_gradient_for_shear(
        &self,
        moment_row: usize,
        shear: f64,
    ) -> Option<([f64; 3], [f64; 3])> {
        let h = self.total_thickness();
        if !h.is_finite() || h <= 0.0 {
            return None;
        }
        let a = self.a_matrix();
        let b = self.b_matrix();
        let d = self.d_matrix();
        let mut system = [[0.0_f64; 6]; 6];
        for (i, (a_row, (b_row, d_row))) in a.iter().zip(b.iter().zip(d.iter())).enumerate() {
            for (j, ((&a_ij, &b_ij), &d_ij)) in
                a_row.iter().zip(b_row.iter()).zip(d_row.iter()).enumerate()
            {
                system[i][j] = a_ij;
                system[i][j + 3] = b_ij / h;
                system[i + 3][j] = b_ij / h;
                system[i + 3][j + 3] = d_ij / (h * h);
            }
        }
        let mut rhs = [0.0_f64; 6];
        rhs[3 + moment_row] = shear / h;
        let x = solve_linear_6(system, rhs)?;
        Some(([x[0], x[1], x[2]], [x[3] / h, x[4] / h, x[5] / h]))
    }
}

/// Interlaminar shear stresses `[τ_xz, τ_yz]` at every ply interface *and* at
/// the two outer surfaces (`n + 1` entries for `n` plies, bottom to top).
///
/// The outer entries are the traction-free boundary conditions: the bottom
/// value is exactly zero and the top value is zero up to round-off, which is
/// the equilibrium check `∫ ∂σ/∂x dz = ∂N/∂x = 0`.
fn interface_shear_profile(
    laminate: &Laminate,
    shear_resultant: [f64; 2],
) -> Result<Vec<[f64; 2]>, crate::Error> {
    let singular = || {
        crate::Error::singular(
            "interlaminar_shear_estimate",
            "the extension-bending (ABD) stiffness matrix is singular",
        )
    };
    // x-direction: dM_x/dx = Q_x;  y-direction: dM_y/dy = Q_y.
    let (eps_x, kappa_x) = laminate
        .strain_gradient_for_shear(0, shear_resultant[0])
        .ok_or_else(singular)?;
    let (eps_y, kappa_y) = laminate
        .strain_gradient_for_shear(1, shear_resultant[1])
        .ok_or_else(singular)?;

    let z = laminate.ply_z_coords();
    let mut profile = Vec::with_capacity(laminate.plies.len() + 1);
    let mut tau = [0.0_f64; 2];
    profile.push(tau);
    for (k, (ply, angle)) in laminate.plies.iter().enumerate() {
        let qbar = ply.q_bar_matrix(*angle);
        let dz = z[k + 1] - z[k];
        // (z_k+1² - z_k²) / 2 written as dz · z_mid to avoid cancellation.
        let dz_moment = dz * 0.5 * (z[k + 1] + z[k]);
        // ∂τ_xz/∂z = -∂σ_x/∂x  with  ∂σ_x/∂x = row 0 of Q̄ · (∂ε⁰ + z·∂κ)
        let d_sigma_x =
            qbar_row_dot(&qbar, 0, &eps_x) * dz + qbar_row_dot(&qbar, 0, &kappa_x) * dz_moment;
        // ∂τ_yz/∂z = -∂σ_y/∂y  with  ∂σ_y/∂y = row 1 of Q̄ · (∂ε⁰ + z·∂κ)
        let d_sigma_y =
            qbar_row_dot(&qbar, 1, &eps_y) * dz + qbar_row_dot(&qbar, 1, &kappa_y) * dz_moment;
        tau[0] -= d_sigma_x;
        tau[1] -= d_sigma_y;
        profile.push(tau);
    }
    Ok(profile)
}

/// Estimate the interlaminar shear stresses `[τ_xz, τ_yz]` at the interface
/// between ply `k` and ply `k + 1` for every interface of the laminate
/// (plies are ordered bottom to top, so the result has `n - 1` entries and is
/// empty for fewer than two plies).
///
/// `shear_resultant = [Q_x, Q_y]` is the transverse shear force per unit
/// width (N/m for SI inputs).  The result has the units of `Q / thickness`,
/// i.e. Pa for SI inputs.
///
/// The stresses follow from integrating the local equilibrium equations
/// through the thickness from the traction-free bottom surface,
///
/// ```text
/// ∂τ_xz/∂z = -∂σ_x/∂x ,    ∂τ_yz/∂z = -∂σ_y/∂y ,
/// ```
///
/// where the in-plane stress gradients follow from the ply stiffnesses and
/// the strain gradients of the laminate, `∂σ = Q̄_k·(∂ε⁰ + z·∂κ)`.  The strain
/// gradients are obtained from the extension-bending equations with constant
/// membrane resultants and a moment gradient equal to the shear force
/// (`∂N = 0`, `∂M_x/∂x = Q_x` for the x-direction, `∂M_y/∂y = Q_y` for the
/// y-direction), using the full 6×6 ABD matrix so unsymmetric laminates
/// (`B ≠ 0`) are handled as well.  For a symmetric laminate this reduces to
///
/// ```text
/// τ_xz(z_k) = -Σ_{j<=k} (row0(Q̄_j)·∂κ/∂x) · (z_j² - z_{j-1}²)/2 ,
/// ∂κ/∂x = D⁻¹·[Q_x, 0, 0]ᵀ ,
/// ```
///
/// and a homogeneous plate gives the classical parabola
/// `τ(z) = 3Q/(2h)·(1 - (2z/h)²)` with `τ_max = 1.5·Q/h`.  The estimate
/// neglects the in-plane shear-stress gradients `∂τ_xy` that would couple the
/// two directions, i.e. it is exact for cylindrical bending in either
/// direction.
///
/// # Errors
///
/// Returns [`crate::Error::Singularity`] when the laminate stiffness matrix
/// is singular (for example zero total thickness or vanishing ply moduli).
pub fn interlaminar_shear_estimate(
    laminate: &Laminate,
    shear_resultant: [f64; 2],
) -> Result<Vec<[f64; 2]>, crate::Error> {
    let n = laminate.plies.len();
    if n < 2 {
        return Ok(Vec::new());
    }
    let profile = interface_shear_profile(laminate, shear_resultant)?;
    // Drop the two traction-free outer surfaces; keep the n - 1 interfaces.
    Ok(profile[1..n].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composite::{OrthotropicPly, mat3_inv};

    /// Laminate built from `angles`, each ply subdivided into `sublayers`
    /// identical sub-plies (only to sample the through-thickness profile
    /// more finely; the mechanics are unchanged).
    fn laminate_from_angles(angles: &[f64], sublayers: usize) -> Laminate {
        let ply_thickness = 0.125e-3 / sublayers as f64;
        let mut lam = Laminate::new();
        for &angle in angles {
            for _ in 0..sublayers {
                lam.add_ply(
                    OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, ply_thickness),
                    angle,
                );
            }
        }
        lam
    }

    /// Independent closed-form reference for symmetric laminates (B = 0):
    /// dkappa/dx = D^-1 [Qx,0,0], dkappa/dy = D^-1 [0,Qy,0],
    /// tau(z_k) = -sum_{j<=k} row(Qbar_j) . dkappa * (z_j^2 - z_{j-1}^2)/2.
    fn reference_shear_symmetric(lam: &Laminate, qx: f64, qy: f64) -> Vec<[f64; 2]> {
        let d_inv = mat3_inv(&lam.d_matrix()).expect("D must be invertible");
        let dk_dx = [d_inv[0][0] * qx, d_inv[1][0] * qx, d_inv[2][0] * qx];
        let dk_dy = [d_inv[0][1] * qy, d_inv[1][1] * qy, d_inv[2][1] * qy];
        let z = lam.ply_z_coords();
        let mut acc = [0.0_f64; 2];
        let mut out = Vec::new();
        for (k, (ply, angle)) in lam.plies.iter().enumerate().take(lam.plies.len() - 1) {
            let qb = ply.q_bar_matrix(*angle);
            let moment = 0.5 * (z[k + 1].powi(2) - z[k].powi(2));
            acc[0] -= (qb[0][0] * dk_dx[0] + qb[0][1] * dk_dx[1] + qb[0][2] * dk_dx[2]) * moment;
            acc[1] -= (qb[1][0] * dk_dy[0] + qb[1][1] * dk_dy[1] + qb[1][2] * dk_dy[2]) * moment;
            out.push(acc);
        }
        out
    }

    #[test]
    fn test_interlaminar_shear_single_ply() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);

        let shear = interlaminar_shear_estimate(&lam, [100.0, 0.0]).unwrap();
        assert!(shear.is_empty(), "single ply has no interfaces");
    }

    #[test]
    fn test_interlaminar_shear_two_plies() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);

        let shear = interlaminar_shear_estimate(&lam, [100.0, 0.0]).unwrap();
        assert_eq!(shear.len(), 1);
    }

    #[test]
    fn test_interlaminar_shear_three_plies() {
        let ply = OrthotropicPly::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.125e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 90.0);
        lam.add_ply(ply, 0.0);

        let shear = interlaminar_shear_estimate(&lam, [100.0, 50.0]).unwrap();
        assert_eq!(shear.len(), 2);
    }

    /// For a homogeneous plate tau_xz is the classical parabola, with
    /// tau_xz(mid-plane) = 1.5 Qx / h.
    #[test]
    fn test_interlaminar_shear_homogeneous_plate_is_parabolic() {
        let lam = laminate_from_angles(&[0.0], 4);
        let h = lam.total_thickness();
        let qx = 1000.0;
        let shear = interlaminar_shear_estimate(&lam, [qx, 0.0]).unwrap();
        assert_eq!(shear.len(), 3);
        let tau_mid = 1.5 * qx / h;
        assert!(
            (shear[1][0] - tau_mid).abs() / tau_mid < 1e-9,
            "tau_xz(mid) = {}, expected {tau_mid}",
            shear[1][0]
        );
        // Interfaces at z = -h/4 and +h/4: 1.5 Q/h * (1 - (2z/h)^2) = 1.125 Q/h.
        let tau_quarter = 1.125 * qx / h;
        for k in [0, 2] {
            assert!(
                (shear[k][0] - tau_quarter).abs() / tau_quarter < 1e-9,
                "interface {k}: {} vs {tau_quarter}",
                shear[k][0]
            );
        }
        for interface in &shear {
            assert!(interface[1].abs() < 1e-9 * tau_mid, "tau_yz must vanish");
        }
    }

    /// The parabola holds for any single-material plate, also for an
    /// off-axis (fully anisotropic) ply, and for both shear directions.
    #[test]
    fn test_interlaminar_shear_homogeneous_off_axis_plate() {
        let lam = laminate_from_angles(&[30.0], 4);
        let h = lam.total_thickness();
        let (qx, qy) = (800.0, -500.0);
        let shear = interlaminar_shear_estimate(&lam, [qx, qy]).unwrap();
        assert!((shear[1][0] - 1.5 * qx / h).abs() < 1e-9 * (1.5 * qx / h).abs());
        assert!((shear[1][1] - 1.5 * qy / h).abs() < 1e-9 * (1.5 * qy / h).abs());
    }

    #[test]
    fn test_interlaminar_shear_matches_symmetric_closed_form() {
        let (qx, qy) = (1200.0, -700.0);
        for angles in [
            vec![0.0, 90.0, 90.0, 0.0],
            vec![0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0],
            vec![30.0, -30.0, -30.0, 30.0],
        ] {
            let lam = laminate_from_angles(&angles, 1);
            assert!(lam.is_symmetric());
            let got = interlaminar_shear_estimate(&lam, [qx, qy]).unwrap();
            let want = reference_shear_symmetric(&lam, qx, qy);
            assert_eq!(got.len(), want.len());
            let scale = 1.5 * qx.abs().max(qy.abs()) / lam.total_thickness();
            for (k, (g, w)) in got.iter().zip(want.iter()).enumerate() {
                for dir in 0..2 {
                    assert!(
                        (g[dir] - w[dir]).abs() < 1e-9 * scale,
                        "{angles:?} interface {k} dir {dir}: {} vs {}",
                        g[dir],
                        w[dir]
                    );
                }
            }
        }
    }

    /// Cross-ply stacking changes the distribution: the mid-plane value of a
    /// [0/90]s plate is no longer the homogeneous 1.5 Q/h, and it grows with
    /// the stiff 0-degree skins.
    #[test]
    fn test_interlaminar_shear_cross_ply_differs_from_homogeneous() {
        let lam = laminate_from_angles(&[0.0, 90.0, 90.0, 0.0], 1);
        let h = lam.total_thickness();
        let qx = 1000.0;
        let shear = interlaminar_shear_estimate(&lam, [qx, 0.0]).unwrap();
        let tau_mid = shear[1][0];
        assert!(
            (tau_mid - 1.5 * qx / h).abs() > 0.05 * (1.5 * qx / h),
            "[0/90]s must deviate from the homogeneous parabola: {tau_mid}"
        );
        assert!(tau_mid > 0.0);
    }

    /// Equilibrium of the through-thickness profile: both outer surfaces are
    /// traction free and the shear stress integrates to the applied shear
    /// force, for symmetric and unsymmetric (B != 0) laminates alike.
    #[test]
    fn test_interlaminar_shear_profile_satisfies_equilibrium() {
        let (qx, qy) = (1500.0, 900.0);
        for angles in [
            vec![0.0, 90.0, 90.0, 0.0],
            vec![0.0, 90.0],
            vec![0.0, 45.0, 90.0],
            vec![30.0, -45.0, 10.0, 60.0],
        ] {
            let lam = laminate_from_angles(&angles, 60);
            let profile = interface_shear_profile(&lam, [qx, qy]).unwrap();
            let n = lam.plies.len();
            assert_eq!(profile.len(), n + 1);
            assert_eq!(profile[0], [0.0, 0.0], "bottom surface is traction free");

            let peak = profile.iter().fold([0.0_f64; 2], |m, t| {
                [m[0].max(t[0].abs()), m[1].max(t[1].abs())]
            });
            assert!(
                profile[n][0].abs() < 1e-9 * peak[0] && profile[n][1].abs() < 1e-9 * peak[1],
                "{angles:?}: top surface not traction free: {:?}",
                profile[n]
            );

            let z = lam.ply_z_coords();
            let mut resultant = [0.0_f64; 2];
            for k in 0..n {
                let dz = z[k + 1] - z[k];
                resultant[0] += 0.5 * (profile[k][0] + profile[k + 1][0]) * dz;
                resultant[1] += 0.5 * (profile[k][1] + profile[k + 1][1]) * dz;
            }
            assert!(
                (resultant[0] - qx).abs() < 2e-3 * qx,
                "{angles:?}: integral of tau_xz = {}, Qx = {qx}",
                resultant[0]
            );
            assert!(
                (resultant[1] - qy).abs() < 2e-3 * qy,
                "{angles:?}: integral of tau_yz = {}, Qy = {qy}",
                resultant[1]
            );
        }
    }

    #[test]
    fn test_interlaminar_shear_is_linear_in_the_load() {
        let lam = laminate_from_angles(&[0.0, 45.0, -45.0, 90.0], 1);
        let base = interlaminar_shear_estimate(&lam, [100.0, 40.0]).unwrap();
        let doubled = interlaminar_shear_estimate(&lam, [200.0, 80.0]).unwrap();
        let none = interlaminar_shear_estimate(&lam, [0.0, 0.0]).unwrap();
        for ((b, d), z) in base.iter().zip(doubled.iter()).zip(none.iter()) {
            for dir in 0..2 {
                assert!((d[dir] - 2.0 * b[dir]).abs() < 1e-9 * b[dir].abs().max(1.0));
                assert_eq!(z[dir], 0.0, "no shear force means no interlaminar shear");
            }
        }
    }

    #[test]
    fn test_interlaminar_shear_ignores_in_plane_membrane_state() {
        // Regression for the old proxy: a pure membrane state used to return a
        // non-zero "shear".  Now only the transverse shear force drives it.
        let lam = laminate_from_angles(&[0.0, 90.0, 0.0], 1);
        let shear = interlaminar_shear_estimate(&lam, [0.0, 0.0]).unwrap();
        assert!(shear.iter().all(|t| t[0] == 0.0 && t[1] == 0.0));
    }

    #[test]
    fn test_interlaminar_shear_singular_laminate_is_an_error() {
        // Q = diag(1e9, 0, 0): rank-deficient in-plane stiffness.
        let ply = OrthotropicPly::new(1.0e9, 0.0, 0.0, 0.0, 1.0e-3);
        let mut lam = Laminate::new();
        lam.add_ply(ply, 0.0);
        lam.add_ply(ply, 0.0);
        let err = interlaminar_shear_estimate(&lam, [1.0, 1.0]).unwrap_err();
        assert!(
            matches!(err, crate::Error::Singularity { .. }),
            "expected a singularity error, got {err:?}"
        );
    }

    #[test]
    fn test_solve_linear_6_known_solution_and_singular_detection() {
        let mut m = [[0.0_f64; 6]; 6];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = 1.0 / (1.0 + (i as f64 - j as f64).abs()) + if i == j { 3.0 } else { 0.0 };
            }
        }
        let x_true = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0];
        let mut rhs = [0.0_f64; 6];
        for (i, r) in rhs.iter_mut().enumerate() {
            *r = m[i].iter().zip(x_true.iter()).map(|(a, b)| a * b).sum();
        }
        let x = solve_linear_6(m, rhs).expect("well-conditioned system");
        for (got, want) in x.iter().zip(x_true.iter()) {
            assert!((got - want).abs() < 1e-12, "x = {x:?}");
        }
        // Row pivoting: a permutation matrix has a zero on every diagonal
        // entry but is perfectly invertible.
        let mut permutation = [[0.0_f64; 6]; 6];
        for pair in 0..3 {
            permutation[2 * pair][2 * pair + 1] = 1.0;
            permutation[2 * pair + 1][2 * pair] = 1.0;
        }
        let permuted_rhs = [
            x_true[1], x_true[0], x_true[3], x_true[2], x_true[5], x_true[4],
        ];
        let y = solve_linear_6(permutation, permuted_rhs).expect("permutation is invertible");
        for (got, want) in y.iter().zip(x_true.iter()) {
            assert!((got - want).abs() < 1e-14, "y = {y:?}");
        }
        // Two identical rows make the matrix singular.
        let mut singular = m;
        singular[4] = singular[1];
        assert!(solve_linear_6(singular, rhs).is_none());
        assert!(solve_linear_6([[0.0; 6]; 6], rhs).is_none());
    }
}
