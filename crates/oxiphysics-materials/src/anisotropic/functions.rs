//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Ply;

/// Symmetric balanced laminate: computes the \[A\] membrane stiffness matrix.
///
/// `[A] = Σ Q̄_k · h_k`  (sum over all plies)
///
/// Returns the 3×3 A matrix as `[[A11, A12, A16], [A12, A22, A26], [A16, A26, A66]]`.
pub fn laminate_a_matrix(plies: &[Ply]) -> [[f64; 3]; 3] {
    let mut a = [[0.0_f64; 3]; 3];
    for ply in plies {
        let qb = ply.transformed_stiffness();
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += qb[i][j] * ply.thickness;
            }
        }
    }
    a
}
/// Effective in-plane Young's modulus `Ex` for a symmetric laminate.
///
/// Uses the inverse of `[A]`: `Ex = (A11·A22 − A12²) / (A22 · h)`.
pub fn laminate_effective_ex(plies: &[Ply]) -> f64 {
    let a = laminate_a_matrix(plies);
    let h: f64 = plies.iter().map(|p| p.thickness).sum();
    if h < f64::EPSILON {
        return 0.0;
    }
    let det = a[0][0] * a[1][1] - a[0][1] * a[0][1];
    det / (a[1][1] * h)
}
/// Compute a simplified pole figure intensity map.
///
/// A pole figure maps the distribution of a crystallographic direction
/// {hkl} in specimen coordinates.  For N crystal orientations (unit vectors),
/// this function accumulates intensities in a 2D (α, β) grid where
/// α = polar angle from z-axis (0..π/2) and β = azimuth (0..2π).
///
/// # Arguments
/// * `orientations` – list of unit vectors (pole directions in specimen frame).
/// * `n_alpha`      – number of polar angle bins.
/// * `n_beta`       – number of azimuth bins.
///
/// Returns a 2D histogram (n_alpha × n_beta) of pole intensities.
pub fn pole_figure(orientations: &[[f64; 3]], n_alpha: usize, n_beta: usize) -> Vec<Vec<u32>> {
    let mut pf = vec![vec![0u32; n_beta]; n_alpha];
    for &v in orientations {
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if norm < 1e-15 {
            continue;
        }
        let vn = [v[0] / norm, v[1] / norm, v[2] / norm];
        let alpha = vn[2].clamp(-1.0, 1.0).acos();
        let beta = vn[1].atan2(vn[0]) + std::f64::consts::PI;
        let ia = ((alpha / (std::f64::consts::PI / 2.0)) * n_alpha as f64) as usize;
        let ib = ((beta / (2.0 * std::f64::consts::PI)) * n_beta as f64) as usize;
        let ia = ia.min(n_alpha - 1);
        let ib = ib.min(n_beta - 1);
        pf[ia][ib] += 1;
    }
    pf
}
/// Compute the texture index of an orientation distribution function (ODF).
///
/// The texture index is defined as T = ∫ f(g)² dg, where f(g) is the ODF.
/// Here it is approximated as the L² norm squared of the normalised histogram
/// of orientation intensities.  A random texture has T = 1.
///
/// # Arguments
/// * `pole_fig` – 2D histogram (α × β bins) from [`pole_figure`].
///
/// Returns the texture index (≥ 1 for any texture, = 1 for random texture).
pub fn texture_index(pole_fig: &[Vec<u32>]) -> f64 {
    let total: u32 = pole_fig.iter().flat_map(|row| row.iter()).sum();
    if total == 0 {
        return 1.0;
    }
    let n_bins: usize = pole_fig.iter().map(|row| row.len()).sum();
    if n_bins == 0 {
        return 1.0;
    }
    let sum_sq: f64 = pole_fig
        .iter()
        .flat_map(|row| row.iter())
        .map(|&c| {
            let f = c as f64 / total as f64 * n_bins as f64;
            f * f
        })
        .sum();
    sum_sq / n_bins as f64
}
/// Coupling stiffness matrix \[B\] for a general laminate.
///
/// `[B] = 1/2 * Σ Q̄_k * (z_k² - z_{k-1}²)`
///
/// For symmetric laminates B = 0. The midplane is taken at z=0.
/// Ply positions are computed from the bottom of the stack.
///
/// Returns the 3×3 B matrix as `[[B11, B12, B16], [B12, B22, B26], [B16, B26, B66]]`.
pub fn laminate_b_matrix(plies: &[Ply]) -> [[f64; 3]; 3] {
    let total_thickness: f64 = plies.iter().map(|p| p.thickness).sum();
    let mut z_bottom = -total_thickness / 2.0;
    let mut b = [[0.0_f64; 3]; 3];
    for ply in plies {
        let z_top = z_bottom + ply.thickness;
        let qb = ply.transformed_stiffness();
        let factor = 0.5 * (z_top * z_top - z_bottom * z_bottom);
        for i in 0..3 {
            for j in 0..3 {
                b[i][j] += qb[i][j] * factor;
            }
        }
        z_bottom = z_top;
    }
    b
}
/// Bending stiffness matrix \[D\] for a general laminate.
///
/// `[D] = 1/3 * Σ Q̄_k * (z_k³ - z_{k-1}³)`
///
/// Returns the 3×3 D matrix.
pub fn laminate_d_matrix(plies: &[Ply]) -> [[f64; 3]; 3] {
    let total_thickness: f64 = plies.iter().map(|p| p.thickness).sum();
    let mut z_bottom = -total_thickness / 2.0;
    let mut d = [[0.0_f64; 3]; 3];
    for ply in plies {
        let z_top = z_bottom + ply.thickness;
        let qb = ply.transformed_stiffness();
        let factor = (z_top.powi(3) - z_bottom.powi(3)) / 3.0;
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] += qb[i][j] * factor;
            }
        }
        z_bottom = z_top;
    }
    d
}
/// Effective transverse modulus Ey for a symmetric laminate.
///
/// `Ey = (A11*A22 - A12²) / (A11 * h)`.
pub fn laminate_effective_ey(plies: &[Ply]) -> f64 {
    let a = laminate_a_matrix(plies);
    let h: f64 = plies.iter().map(|p| p.thickness).sum();
    if h < f64::EPSILON {
        return 0.0;
    }
    let det = a[0][0] * a[1][1] - a[0][1] * a[0][1];
    det / (a[0][0] * h)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::anisotropic::BraidedComposite;
    use crate::anisotropic::DiffusionTensor;
    use crate::anisotropic::FailureEnvelope;
    use crate::anisotropic::HashinFailureCriteria;
    use crate::anisotropic::HillYieldCriterion;
    use crate::anisotropic::LaRCFailureCriteria;
    use crate::anisotropic::OrthotropicMaterial;
    use crate::anisotropic::PuckFailureCriteria;
    use crate::anisotropic::SymmetryOperation;
    use crate::anisotropic::ThermalConductivityTensor;
    use crate::anisotropic::TransverselyIsotropic;
    use crate::anisotropic::WovenLamina;
    use crate::anisotropic::{MonoclinicCompliance, MonoclinicMaterial};
    /// Helper: multiply two 6×6 matrices.
    fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
        let mut c = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }
    #[test]
    fn test_orthotropic_compliance_diagonal() {
        let mat = OrthotropicMaterial::carbon_fiber_epoxy();
        let s = mat.compliance_matrix();
        let tol = 1.0e-20;
        assert!(
            (s[0][0] - 1.0 / mat.e1).abs() < tol,
            "S[0][0] should be 1/E1"
        );
        assert!(
            (s[1][1] - 1.0 / mat.e2).abs() < tol,
            "S[1][1] should be 1/E2"
        );
        assert!(
            (s[2][2] - 1.0 / mat.e3).abs() < tol,
            "S[2][2] should be 1/E3"
        );
        assert!(
            (s[3][3] - 1.0 / mat.g12).abs() < tol,
            "S[3][3] should be 1/G12"
        );
        assert!(
            (s[4][4] - 1.0 / mat.g23).abs() < tol,
            "S[4][4] should be 1/G23"
        );
        assert!(
            (s[5][5] - 1.0 / mat.g13).abs() < tol,
            "S[5][5] should be 1/G13"
        );
    }
    #[test]
    fn test_orthotropic_symmetry_carbon_fiber() {
        let mat = OrthotropicMaterial::carbon_fiber_epoxy();
        let nu21 = mat.nu12 * mat.e2 / mat.e1;
        assert!(
            (mat.nu12 / mat.e1 - nu21 / mat.e2).abs() < 1.0e-30,
            "nu12/E1 should equal nu21/E2"
        );
    }
    #[test]
    fn test_orthotropic_stiffness_inverts_compliance() {
        let mat = OrthotropicMaterial::carbon_fiber_epoxy();
        let d = mat.constitutive_matrix();
        let s = mat.compliance_matrix();
        let ds = mat6_mul(&d, &s);
        for (i, row) in ds.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1.0e-10,
                    "D*S[{i}][{j}] = {} but expected {expected}",
                    val
                );
            }
        }
    }
    #[test]
    fn test_transverse_isotropic_gt() {
        let mat = TransverselyIsotropic::new(140.0e9, 10.0e9, 0.3, 0.4, 5.0e9);
        let gt_expected = mat.et / (2.0 * (1.0 + mat.nut));
        assert!(
            (mat.gt() - gt_expected).abs() < 1.0e-6,
            "gt() = {} but expected {gt_expected}",
            mat.gt()
        );
    }
    #[test]
    fn test_transverse_to_orthotropic() {
        let ti = TransverselyIsotropic::new(140.0e9, 10.0e9, 0.3, 0.4, 5.0e9);
        let orth = ti.to_orthotropic();
        assert_eq!(
            orth.e2, orth.e3,
            "E2 should equal E3 for transversely isotropic"
        );
        assert_eq!(
            orth.g12, orth.g13,
            "G12 should equal G13 for transversely isotropic"
        );
        assert_eq!(
            orth.nu12, orth.nu13,
            "nu12 should equal nu13 for transversely isotropic"
        );
        assert!((orth.e2 - ti.et).abs() < 1.0e-10);
        assert!((orth.g12 - ti.ga).abs() < 1.0e-10);
        assert!((orth.g23 - ti.gt()).abs() < 1.0e-10);
    }
    #[test]
    fn test_isotropic_is_special_case() {
        use crate::elastic::LinearElastic;
        let e = 200.0e9;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let orth = OrthotropicMaterial::new(e, e, e, nu, nu, nu, g, g, g);
        let d_orth = orth.constitutive_matrix();
        let iso = LinearElastic::new(e, nu);
        let d_iso = iso.stress_strain_matrix_3d();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (d_orth[i][j] - d_iso[i][j]).abs() < 1.0e3,
                    "Normal block D_orth[{i}][{j}]={} vs D_iso[{i}][{j}]={}",
                    d_orth[i][j],
                    d_iso[i][j]
                );
            }
        }
        assert!((d_orth[3][3] - g).abs() < 1.0e3);
        assert!((d_orth[4][4] - g).abs() < 1.0e3);
        assert!((d_orth[5][5] - g).abs() < 1.0e3);
    }
    #[test]
    fn test_woven_balanced() {
        let lam = WovenLamina::balanced_plain_weave(0.5, 73.0e9, 3.5e9, 0.22, 0.35);
        assert!(
            (lam.e11 - lam.e22).abs() < 1.0e-6,
            "E11 should equal E22 for balanced weave: E11={}, E22={}",
            lam.e11,
            lam.e22
        );
    }
    #[test]
    fn test_compliance_positive_diagonal() {
        let mats = [
            OrthotropicMaterial::carbon_fiber_epoxy(),
            OrthotropicMaterial::douglas_fir(),
        ];
        for mat in &mats {
            let s = mat.compliance_matrix();
            for (i, row) in s.iter().enumerate() {
                assert!(
                    row[i] > 0.0,
                    "S[{i}][{i}] should be positive, got {}",
                    row[i]
                );
            }
        }
    }
    #[test]
    fn test_ply_reduced_stiffness_positive_diagonal() {
        let ply = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 0.0);
        let [q11, q22, _q12, q66] = ply.reduced_stiffness();
        assert!(q11 > 0.0, "Q11 should be positive");
        assert!(q22 > 0.0, "Q22 should be positive");
        assert!(q66 > 0.0, "Q66 should be positive");
    }
    #[test]
    fn test_ply_transformed_stiffness_0_deg_same_as_reduced() {
        let ply = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 0.0);
        let [q11, q22, q12, q66] = ply.reduced_stiffness();
        let qb = ply.transformed_stiffness();
        let tol = 1.0e3;
        assert!(
            (qb[0][0] - q11).abs() < tol,
            "Qb11={}, Q11={}",
            qb[0][0],
            q11
        );
        assert!(
            (qb[1][1] - q22).abs() < tol,
            "Qb22={}, Q22={}",
            qb[1][1],
            q22
        );
        assert!(
            (qb[0][1] - q12).abs() < tol,
            "Qb12={}, Q12={}",
            qb[0][1],
            q12
        );
        assert!(
            (qb[2][2] - q66).abs() < tol,
            "Qb66={}, Q66={}",
            qb[2][2],
            q66
        );
        assert!(
            qb[0][2].abs() < tol,
            "Qb16 should be 0 at 0°, got {}",
            qb[0][2]
        );
    }
    #[test]
    fn test_ply_transformed_stiffness_90_deg() {
        let ply0 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 0.0);
        let ply90 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 90.0);
        let qb0 = ply0.transformed_stiffness();
        let qb90 = ply90.transformed_stiffness();
        let tol = 1.0e3;
        assert!(
            (qb90[0][0] - qb0[1][1]).abs() < tol,
            "Qb11 at 90° = Qb22 at 0°"
        );
        assert!(
            (qb90[1][1] - qb0[0][0]).abs() < tol,
            "Qb22 at 90° = Qb11 at 0°"
        );
    }
    #[test]
    fn test_laminate_a_matrix_symmetric() {
        let ply0 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 0.0);
        let ply90 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 90.0);
        let plies = [ply0, ply90, ply90, ply0];
        let a = laminate_a_matrix(&plies);
        assert!((a[0][1] - a[1][0]).abs() < 1.0, "A should be symmetric");
        assert!((a[0][2] - a[2][0]).abs() < 1.0);
    }
    #[test]
    fn test_laminate_effective_ex_positive() {
        let ply = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.25e-3, 0.0);
        let plies = [ply, ply, ply, ply];
        let ex = laminate_effective_ex(&plies);
        assert!(ex > 0.0, "Effective Ex should be positive, got {ex}");
    }
    #[test]
    fn test_laminate_effective_ex_unidirectional() {
        let e1 = 140.0e9;
        let ply = Ply::new(e1, 10.0e9, 0.3, 5.0e9, 0.25e-3, 0.0);
        let plies = [ply, ply, ply, ply];
        let ex = laminate_effective_ex(&plies);
        assert!((ex - e1).abs() / e1 < 0.01, "Ex={ex}, E1={e1}");
    }
    #[test]
    fn test_monoclinic_e1_from_s11() {
        let e1 = 100.0e9;
        let mat = MonoclinicMaterial::new(MonoclinicCompliance {
            s11: 1.0 / e1,
            s12: -0.3 / e1,
            s13: -0.1 / e1,
            s16: 0.0,
            s22: 1.0 / 70.0e9,
            s23: -0.2 / 70.0e9,
            s26: 0.0,
            s33: 1.0 / 50.0e9,
            s36: 0.0,
            s44: 1.0 / 30.0e9,
            s45: 0.0,
            s55: 1.0 / 25.0e9,
            s66: 1.0 / 40.0e9,
        });
        assert!(
            (mat.e1() - e1).abs() < 1.0,
            "E1={}, expected {e1}",
            mat.e1()
        );
    }
    #[test]
    fn test_monoclinic_physical_validity() {
        let mat = MonoclinicMaterial::new(MonoclinicCompliance {
            s11: 1.0 / 100.0e9,
            s12: -0.3 / 100.0e9,
            s13: -0.1 / 100.0e9,
            s16: 0.0,
            s22: 1.0 / 70.0e9,
            s23: -0.2 / 70.0e9,
            s26: 0.0,
            s33: 1.0 / 50.0e9,
            s36: 0.0,
            s44: 1.0 / 30.0e9,
            s45: 0.0,
            s55: 1.0 / 25.0e9,
            s66: 1.0 / 40.0e9,
        });
        assert!(mat.is_physically_valid());
    }
    #[test]
    fn test_thermal_conductivity_orthotropic_heat_flux() {
        let tc = ThermalConductivityTensor::orthotropic(10.0, 5.0, 2.0);
        let grad_t = [1.0, 0.0, 0.0];
        let q = tc.heat_flux(grad_t);
        assert!(
            (q[0] + 10.0).abs() < 1e-12,
            "q_x should be -10, got {}",
            q[0]
        );
        assert!(q[1].abs() < 1e-12, "q_y should be 0");
        assert!(q[2].abs() < 1e-12, "q_z should be 0");
    }
    #[test]
    fn test_thermal_conductivity_effective_along_x() {
        let tc = ThermalConductivityTensor::orthotropic(10.0, 5.0, 2.0);
        let n = [1.0, 0.0, 0.0];
        let keff = tc.effective_conductivity(n);
        assert!(
            (keff - 10.0).abs() < 1e-12,
            "κ_eff along x should be 10, got {keff}"
        );
    }
    #[test]
    fn test_thermal_conductivity_mean() {
        let tc = ThermalConductivityTensor::orthotropic(10.0, 4.0, 1.0);
        assert!(
            (tc.mean_conductivity() - 5.0).abs() < 1e-12,
            "mean should be 5"
        );
    }
    #[test]
    fn test_diffusion_tensor_flux() {
        let dt = DiffusionTensor::orthotropic(1.0, 2.0, 3.0);
        let grad_c = [1.0, 0.0, 0.0];
        let j = dt.diffusive_flux(grad_c);
        assert!((j[0] + 1.0).abs() < 1e-12, "Jx = -Dxx*1 = -1");
        assert!(j[1].abs() < 1e-12);
    }
    #[test]
    fn test_diffusion_tensor_msd() {
        let dt = DiffusionTensor::orthotropic(2.0, 2.0, 2.0);
        let n = [1.0, 0.0, 0.0];
        let msd = dt.msd_along(n, 1.0);
        assert!((msd - 4.0).abs() < 1e-12, "MSD = 2*D*t = 4, got {msd}");
    }
    #[test]
    fn test_diffusion_tensor_mean() {
        let dt = DiffusionTensor::orthotropic(1.0, 2.0, 3.0);
        assert!((dt.mean_diffusivity() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_symmetry_identity_apply() {
        let op = SymmetryOperation::identity();
        let v = [1.0, 2.0, 3.0];
        let r = op.apply(v);
        assert!((r[0] - v[0]).abs() < 1e-12);
        assert!((r[1] - v[1]).abs() < 1e-12);
        assert!((r[2] - v[2]).abs() < 1e-12);
    }
    #[test]
    fn test_symmetry_inversion_det() {
        let op = SymmetryOperation::inversion();
        assert!((op.determinant() + 1.0).abs() < 1e-12, "inversion det = -1");
        assert!(!op.is_proper_rotation());
    }
    #[test]
    fn test_symmetry_c2z_is_proper() {
        let op = SymmetryOperation::c2_z();
        assert!(op.is_proper_rotation(), "C2z should be a proper rotation");
    }
    #[test]
    fn test_symmetry_compose_identity() {
        let id = SymmetryOperation::identity();
        let c2 = SymmetryOperation::c2_z();
        let composed = id.compose(&c2);
        let v = [1.0, 0.0, 0.0];
        let r1 = c2.apply(v);
        let r2 = composed.apply(v);
        for k in 0..3 {
            assert!((r1[k] - r2[k]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_symmetry_mirror_z_apply() {
        let op = SymmetryOperation::mirror_z();
        let v = [1.0, 2.0, 3.0];
        let r = op.apply(v);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 2.0).abs() < 1e-12);
        assert!((r[2] + 3.0).abs() < 1e-12, "z should be negated");
    }
    #[test]
    fn test_pole_figure_single_direction() {
        let orientations = vec![[0.0, 0.0, 1.0]];
        let pf = pole_figure(&orientations, 10, 10);
        let total: u32 = pf.iter().flat_map(|r| r.iter()).sum();
        assert_eq!(total, 1, "Single orientation should give 1 count total");
    }
    #[test]
    fn test_pole_figure_dimensions() {
        let pf = pole_figure(&[], 5, 8);
        assert_eq!(pf.len(), 5);
        for row in &pf {
            assert_eq!(row.len(), 8);
        }
    }
    #[test]
    fn test_texture_index_uniform() {
        let n = 4;
        let pf: Vec<Vec<u32>> = vec![vec![25u32; n]; n];
        let ti = texture_index(&pf);
        assert!(
            (ti - 1.0).abs() < 1e-10,
            "Uniform ODF should have texture index = 1, got {ti}"
        );
    }
    #[test]
    fn test_texture_index_single_peak() {
        let mut pf = vec![vec![0u32; 10]; 10];
        pf[0][0] = 100;
        let ti = texture_index(&pf);
        assert!(
            ti > 1.0,
            "Single-peak texture should have index > 1, got {ti}"
        );
    }
    #[test]
    fn test_hill_isotropic_von_mises() {
        let hill = HillYieldCriterion::isotropic(1.0);
        let yield_stress = 250.0e6_f64;
        let sigma = [yield_stress, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = hill.yield_function(sigma, yield_stress);
        assert!(f.abs() < 1.0, "Hill at yield: f={f}");
    }
    #[test]
    fn test_hill_yield_function_zero_stress() {
        let hill = HillYieldCriterion::from_r_values(2.0, 0.5, 2.0);
        let sigma = [0.0_f64; 6];
        let f = hill.yield_function(sigma, 200.0e6);
        assert!(f < 0.0, "zero stress should be inside yield surface, f={f}");
    }
    #[test]
    fn test_hill_yield_function_at_yield_x() {
        let hill = HillYieldCriterion::isotropic(1.0);
        let sy = 300.0e6_f64;
        let sigma = [sy, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = hill.yield_function(sigma, sy);
        assert!(f.abs() / (sy * sy) < 1e-8, "should be at yield, f={f}");
    }
    #[test]
    fn test_hill_r_values_asymmetry() {
        let hill = HillYieldCriterion::from_r_values(4.0, 0.5, 2.0);
        let sy = 200.0e6_f64;
        let sigma = [sy * 0.5, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = hill.yield_function(sigma, sy);
        assert!(f < 0.0, "sub-yield state should give f<0, f={f}");
    }
    #[test]
    fn test_hill_effective_stress_positive() {
        let hill = HillYieldCriterion::from_r_values(2.0, 2.0, 2.0);
        let sigma = [100.0e6, -50.0e6, 30.0e6, 20.0e6, 10.0e6, 5.0e6];
        let se = hill.effective_stress(sigma);
        assert!(se > 0.0, "effective stress should be positive, got {se}");
    }
    #[test]
    fn test_hashin_fibre_tension_at_threshold() {
        let xt = 1500.0e6_f64;
        let s12 = 70.0e6_f64;
        let hashin = HashinFailureCriteria {
            xt,
            xc: 1200.0e6,
            yt: 50.0e6,
            yc: 150.0e6,
            s12,
            s23: 40.0e6,
        };
        let fi = hashin.fibre_tension(xt, 0.0);
        assert!((fi - 1.0).abs() < 1e-10, "FI at threshold = 1, got {fi}");
    }
    #[test]
    fn test_hashin_fibre_tension_below_threshold() {
        let hashin = HashinFailureCriteria {
            xt: 1500.0e6,
            xc: 1200.0e6,
            yt: 50.0e6,
            yc: 150.0e6,
            s12: 70.0e6,
            s23: 40.0e6,
        };
        let fi = hashin.fibre_tension(500.0e6, 0.0);
        assert!(fi < 1.0, "sub-threshold FI should be < 1, got {fi}");
    }
    #[test]
    fn test_hashin_matrix_tension_at_threshold() {
        let yt = 60.0e6_f64;
        let s12 = 70.0e6_f64;
        let hashin = HashinFailureCriteria {
            xt: 1500.0e6,
            xc: 1200.0e6,
            yt,
            yc: 150.0e6,
            s12,
            s23: 40.0e6,
        };
        let fi = hashin.matrix_tension(yt, 0.0);
        assert!(
            (fi - 1.0).abs() < 1e-10,
            "matrix tension FI at Yt = 1, got {fi}"
        );
    }
    #[test]
    fn test_hashin_failure_index_max_mode() {
        let hashin = HashinFailureCriteria {
            xt: 1500.0e6,
            xc: 1200.0e6,
            yt: 50.0e6,
            yc: 150.0e6,
            s12: 70.0e6,
            s23: 40.0e6,
        };
        let (fi_max, _) = hashin.max_failure_index(100.0e6, -20.0e6, 10.0e6, 5.0e6);
        assert!(fi_max >= 0.0, "max FI should be non-negative, got {fi_max}");
    }
    #[test]
    fn test_puck_fibre_failure_below_threshold() {
        let puck = PuckFailureCriteria::carbon_epoxy_typical();
        let fi = puck.fibre_failure_index(500.0e6, 1500.0e6);
        assert!(fi < 1.0, "below threshold, FI={fi}");
    }
    #[test]
    fn test_puck_fibre_failure_at_threshold() {
        let puck = PuckFailureCriteria::carbon_epoxy_typical();
        let fi = puck.fibre_failure_index(1500.0e6, 1500.0e6);
        assert!((fi - 1.0).abs() < 1e-10, "at threshold FI=1, got {fi}");
    }
    #[test]
    fn test_puck_inter_fibre_failure_zero_stress() {
        let puck = PuckFailureCriteria::carbon_epoxy_typical();
        let fi = puck.inter_fibre_failure_index(0.0, 0.0, 50.0e6);
        assert!(fi < 1.0, "zero stress should be below threshold, FI={fi}");
    }
    #[test]
    fn test_larc_fibre_kinking_negative_axial() {
        let larc = LaRCFailureCriteria::default_carbon_epoxy();
        let fi = larc.fibre_kinking_index(-800.0e6, 50.0e6, 0.0);
        assert!(
            fi >= 0.0,
            "fibre kinking index should be non-negative, got {fi}"
        );
    }
    #[test]
    fn test_larc_matrix_cracking_below_threshold() {
        let larc = LaRCFailureCriteria::default_carbon_epoxy();
        let fi = larc.matrix_cracking_index(20.0e6, 5.0e6, 30.0e6);
        assert!(fi < 1.0, "below threshold, FI={fi}");
    }
    #[test]
    fn test_failure_envelope_biaxial_positive() {
        let env = FailureEnvelope::max_stress(1500.0e6, -1200.0e6, 50.0e6, -150.0e6);
        let fi = env.failure_index_biaxial(750.0e6, 25.0e6);
        assert!(fi < 1.0, "half-strength should be inside envelope, FI={fi}");
    }
    #[test]
    fn test_failure_envelope_at_limit() {
        let env = FailureEnvelope::max_stress(1000.0e6, -800.0e6, 60.0e6, -120.0e6);
        let fi = env.failure_index_biaxial(1000.0e6, 0.0);
        assert!((fi - 1.0).abs() < 1e-10, "at tensile limit, FI=1, got {fi}");
    }
    #[test]
    fn test_braided_composite_in_plane_modulus_positive() {
        let bc = BraidedComposite::triaxial_braid(0.6, 140.0e9, 3.5e9, 0.22, 0.35, 30.0);
        let e = bc.in_plane_modulus();
        assert!(e > 0.0, "in-plane modulus should be positive, E={e}");
    }
    #[test]
    fn test_braided_composite_shear_modulus_positive() {
        let bc = BraidedComposite::triaxial_braid(0.5, 73.0e9, 3.5e9, 0.22, 0.35, 45.0);
        let g = bc.effective_shear_modulus();
        assert!(g > 0.0, "shear modulus should be positive, G={g}");
    }
    #[test]
    fn test_braided_composite_axial_vs_45deg() {
        let bc0 = BraidedComposite::triaxial_braid(0.6, 140.0e9, 3.5e9, 0.22, 0.35, 20.0);
        let bc45 = BraidedComposite::triaxial_braid(0.6, 140.0e9, 3.5e9, 0.22, 0.35, 45.0);
        assert!(
            bc0.in_plane_modulus() > bc45.in_plane_modulus(),
            "lower braid angle should give higher modulus"
        );
    }
    #[test]
    fn test_laminate_b_matrix_symmetric_balanced() {
        let ply0 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 0.0);
        let ply90 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.125e-3, 90.0);
        let plies = [ply0, ply90, ply90, ply0];
        let b = laminate_b_matrix(&plies);
        for (i, row) in b.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < 1.0,
                    "B[{i}][{j}]={} should be ~0 for symmetric laminate",
                    val
                );
            }
        }
    }
    #[test]
    fn test_laminate_d_matrix_positive_diagonal() {
        let ply = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.25e-3, 0.0);
        let plies = [ply, ply, ply, ply];
        let d = laminate_d_matrix(&plies);
        for (i, row) in d.iter().enumerate() {
            assert!(row[i] > 0.0, "D[{i}][{i}]={} should be positive", row[i]);
        }
    }
    #[test]
    fn test_laminate_effective_ey_positive() {
        let ply = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, 0.25e-3, 90.0);
        let plies = [ply, ply, ply, ply];
        let ey = laminate_effective_ey(&plies);
        assert!(ey > 0.0, "Ey should be positive, got {ey}");
    }
    #[test]
    fn test_laminate_quasi_isotropic_a_matrix_nearly_isotropic() {
        let t = 0.125e-3_f64;
        let ply0 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, t, 0.0);
        let ply45 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, t, 45.0);
        let plyp45 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, t, -45.0);
        let ply90 = Ply::new(140.0e9, 10.0e9, 0.3, 5.0e9, t, 90.0);
        let plies = [ply0, ply45, plyp45, ply90, ply90, plyp45, ply45, ply0];
        let a = laminate_a_matrix(&plies);
        let ratio = a[0][0] / a[1][1];
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "A11/A22 = {ratio}, should be ~1 for quasi-isotropic"
        );
    }
    #[test]
    fn test_orthotropic_bulk_modulus_positive() {
        let mat = OrthotropicMaterial::carbon_fiber_epoxy();
        let k = mat.bulk_modulus_longitudinal();
        assert!(k > 0.0, "bulk modulus should be positive, K={k}");
    }
    #[test]
    fn test_orthotropic_symmetry_douglas_fir() {
        let mat = OrthotropicMaterial::douglas_fir();
        assert!(mat.check_symmetry(), "Douglas fir should satisfy symmetry");
    }
    #[test]
    fn test_orthotropic_constitutive_symmetric() {
        let mat = OrthotropicMaterial::carbon_fiber_epoxy();
        let d = mat.constitutive_matrix();
        for (i, row) in d.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - d[j][i]).abs() < 1.0,
                    "D[{i}][{j}]={} vs D[{j}][{i}]={}",
                    val,
                    d[j][i]
                );
            }
        }
    }
    #[test]
    fn test_transverse_isotropic_constitutive_dimensions() {
        let mat = TransverselyIsotropic::new(140.0e9, 10.0e9, 0.3, 0.4, 5.0e9);
        let d = mat.constitutive_matrix();
        for (i, row) in d.iter().enumerate() {
            assert!(row[i] > 0.0, "D[{i}][{i}] should be positive");
        }
    }
    #[test]
    fn test_braided_higher_vf_higher_stiffness() {
        let bc_lo = BraidedComposite::triaxial_braid(0.3, 73.0e9, 3.5e9, 0.22, 0.35, 30.0);
        let bc_hi = BraidedComposite::triaxial_braid(0.6, 73.0e9, 3.5e9, 0.22, 0.35, 30.0);
        assert!(
            bc_hi.in_plane_modulus() > bc_lo.in_plane_modulus(),
            "higher Vf should give higher stiffness"
        );
    }
    #[test]
    fn test_braided_in_plane_bounds() {
        let bc = BraidedComposite::triaxial_braid(0.5, 73.0e9, 3.5e9, 0.22, 0.35, 30.0);
        let e = bc.in_plane_modulus();
        assert!(e > 3.5e9 && e < 73.0e9, "E={e} should be between Em and Ef");
    }
}
