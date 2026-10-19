//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PlateParams;

/// Navier series first-term maximum deflection for a uniformly-loaded
/// simply-supported rectangular plate.
///
/// w_max = 0.00406 · q · a⁴ / D  (square plate)
pub fn max_deflection_uniform_load(load: f64, a: f64, _b: f64, params: &PlateParams) -> f64 {
    let d = params.flexural_rigidity();
    0.00406 * load * a * a * a * a / d
}
/// Natural frequency of simply-supported rectangular plate, mode (m, n).
///
/// f_mn = π²/2 · (m²/a² + n²/b²) · √(D / (ρ·h))
pub fn natural_frequency_mn(m: u32, n: u32, a: f64, b: f64, d: f64, rho: f64, h: f64) -> f64 {
    use std::f64::consts::PI;
    let mf = m as f64;
    let nf = n as f64;
    PI * PI / 2.0 * (mf * mf / (a * a) + nf * nf / (b * b)) * (d / (rho * h)).sqrt()
}
/// Maximum membrane stiffness: N_max = E·h / (1−ν²).
pub fn membrane_stiffness(e: f64, nu: f64, h: f64) -> f64 {
    e * h / (1.0 - nu * nu)
}
/// Classical external pressure at buckling of a long thin cylindrical shell
/// (Bresse's ring-collapse formula).
///
/// p_cr = E·h³ / (4·(1−ν²)·R³) = E·(h/R)³ / (4·(1−ν²))
///
/// This is the two-lobe buckling pressure of an infinitely long tube and the
/// lower bound for tubes of any length.  It is *not* the spherical-shell
/// (Zoelly) expression of [`spherical_shell_buckling`].  For tubes of finite
/// length see [`cylindrical_shell_buckling_finite`].
pub fn cylindrical_shell_buckling(e: f64, nu: f64, r: f64, h: f64) -> f64 {
    e * (h / r).powi(3) / (4.0 * (1.0 - nu * nu))
}
/// Critical external pressure of a finite-length thin cylindrical shell with
/// simply supported (or ring-stiffened) ends.
///
/// Windenburg–Trilling formula for intermediate-length tubes, bounded below by
/// the long-tube Bresse value of [`cylindrical_shell_buckling`]; `length` is the
/// unsupported length `L`.  Returns `NaN` for non-physical inputs or for tubes
/// too short for the formula to apply (`L/D − 0.45·√(h/D) ≤ 0`, `D = 2R`).
///
/// Argument order follows [`cylindrical_shell_buckling`]: `(e, nu, r, h, length)`.
pub fn cylindrical_shell_buckling_finite(e: f64, nu: f64, r: f64, h: f64, length: f64) -> f64 {
    crate::shell_fem::cylinder_buckling_pressure_finite(e, nu, h, r, length)
}
/// Classical external pressure at buckling of a thin spherical shell.
///
/// p_cr = 2·E / √(3·(1−ν²)) · (h/R)²
pub fn spherical_shell_buckling(e: f64, nu: f64, r: f64, h: f64) -> f64 {
    2.0 * e / (3.0 * (1.0 - nu * nu)).sqrt() * (h / r).powi(2)
}
/// Aspect ratio of a plate: a/b.
pub fn aspect_ratio(a: f64, b: f64) -> f64 {
    a / b
}
/// Center deflection of a simply-supported square plate under a concentrated load.
///
/// w = 0.01160 · P · a² / D
pub fn plate_concentrated_load_deflection(p: f64, a: f64, d: f64) -> f64 {
    0.01160 * p * a * a / d
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::shell::*;
    #[test]
    fn test_kirchhoff_flexural_rigidity() {
        let t = 0.01;
        let e = 210e9;
        let nu = 0.3;
        let plate = KirchhoffPlate::new(t, e, nu, 7800.0);
        let d = plate.flexural_rigidity();
        let expected = e * t * t * t / (12.0 * (1.0 - nu * nu));
        assert!(
            (d - expected).abs() / expected < 1e-12,
            "flexural rigidity: got {d}, expected {expected}"
        );
    }
    #[test]
    fn test_kirchhoff_constitutive_symmetric() {
        let plate = KirchhoffPlate::new(0.01, 210e9, 0.3, 7800.0);
        let dm = plate.constitutive_matrix();
        for (i, row) in dm.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - dm[j][i]).abs() < 1e-6,
                    "D_plate not symmetric at [{i},{j}]: {} vs {}",
                    val,
                    dm[j][i]
                );
            }
        }
    }
    #[test]
    fn test_kirchhoff_natural_frequency() {
        let plate = KirchhoffPlate::new(0.01, 210e9, 0.3, 7800.0);
        let omega = plate.natural_frequency_ss(1, 1, 1.0, 1.0);
        assert!(
            omega > 100.0 && omega < 5000.0,
            "natural frequency {omega} rad/s outside expected range [100, 5000]"
        );
    }
    #[test]
    fn test_mindlin_shear_stiffness() {
        let shell = MindlinShell::new(0.02, 200e9, 0.25);
        let g = shell.young_modulus / (2.0 * (1.0 + shell.poisson_ratio));
        let expected = (5.0 / 6.0) * g * shell.thickness;
        let ks = shell.shear_stiffness();
        assert!(
            (ks - expected).abs() / expected < 1e-12,
            "shear stiffness: got {ks}, expected {expected}"
        );
    }
    #[test]
    fn test_mindlin_is_thin() {
        let shell = MindlinShell::new(0.01, 200e9, 0.3);
        assert!(shell.is_thin(1.0), "t/L=0.01 should be thin");
        assert!(!shell.is_thin(0.05), "t/L=0.2 should not be thin");
    }
    #[test]
    fn test_membrane_area_unit_triangle() {
        let nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let area = MembraneTriangle::area(&nodes);
        assert!(
            (area - 0.5).abs() < 1e-14,
            "unit triangle area: got {area}, expected 0.5"
        );
    }
    #[test]
    fn test_membrane_stiffness_positive_definite() {
        let elem = MembraneTriangle::new(0.01, 200e9, 0.3);
        let nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let k = elem.stiffness_2d(&nodes);
        for (i, row) in k.iter().enumerate() {
            assert!(
                row[i] > 0.0,
                "stiffness diagonal k[{i}][{i}] = {} should be positive",
                row[i]
            );
        }
    }
    #[test]
    fn test_membrane_constitutive_symmetric() {
        let elem = MembraneTriangle::new(0.005, 70e9, 0.33);
        let d = elem.constitutive_matrix();
        let expected_d00 = elem.young_modulus / (1.0 - elem.poisson_ratio * elem.poisson_ratio);
        assert!(
            (d[0][0] - expected_d00).abs() / expected_d00 < 1e-12,
            "D[0][0]: got {}, expected {expected_d00}",
            d[0][0]
        );
        for (i, row) in d.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - d[j][i]).abs() < 1e-6,
                    "D not symmetric at [{i},{j}]: {} vs {}",
                    val,
                    d[j][i]
                );
            }
        }
    }
    #[test]
    fn test_plate_params_flexural_rigidity() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let expected = 200e9 * 0.01f64.powi(3) / (12.0 * (1.0 - 0.3f64 * 0.3));
        assert!((d - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_plate_params_rigidity_scales_cubic() {
        let p1 = PlateParams::new(200e9, 0.3, 0.01);
        let p2 = PlateParams::new(200e9, 0.3, 0.02);
        assert!((p2.flexural_rigidity() / p1.flexural_rigidity() - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_plate_params_bending_stiffness_symmetric() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let db = p.bending_stiffness_matrix();
        for (i, row) in db.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - db[j][i]).abs() < 1e-6,
                    "Db not symmetric at [{i},{j}]"
                );
            }
        }
    }
    #[test]
    fn test_plate_params_bending_stiffness_d11() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let db = p.bending_stiffness_matrix();
        assert!((db[0][0] - d).abs() / d < 1e-12);
    }
    #[test]
    fn test_plate_params_bending_stiffness_d12() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let db = p.bending_stiffness_matrix();
        assert!((db[0][1] - d * 0.3).abs() / (d * 0.3) < 1e-12);
    }
    #[test]
    fn test_plate_params_bending_stiffness_d33() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let db = p.bending_stiffness_matrix();
        let expected = d * (1.0 - 0.3) / 2.0;
        assert!((db[2][2] - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_rect_plate_node_positions() {
        let params = PlateParams::new(200e9, 0.3, 0.01);
        let elem = RectPlateElement::new(params, 0.5, 0.25);
        assert!((elem.node_positions[0][0] - (-0.5)).abs() < 1e-15);
        assert!((elem.node_positions[0][1] - (-0.25)).abs() < 1e-15);
    }
    #[test]
    fn test_rect_plate_stiffness_size() {
        let params = PlateParams::new(200e9, 0.3, 0.01);
        let elem = RectPlateElement::new(params, 0.5, 0.5);
        let k = elem.stiffness_matrix_12x12();
        assert_eq!(k.len(), 12);
        assert_eq!(k[0].len(), 12);
    }
    #[test]
    fn test_rect_plate_stiffness_positive_diagonal() {
        let params = PlateParams::new(200e9, 0.3, 0.01);
        let elem = RectPlateElement::new(params, 0.5, 0.5);
        let k = elem.stiffness_matrix_12x12();
        for (i, row) in k.iter().enumerate() {
            assert!(
                row[i] > 0.0,
                "k[{i}][{i}] should be positive, got {}",
                row[i]
            );
        }
    }
    #[test]
    fn test_rect_plate_mass_size() {
        let params = PlateParams::new(200e9, 0.3, 0.01);
        let elem = RectPlateElement::new(params, 0.5, 0.5);
        let m = elem.consistent_mass_matrix(7800.0);
        assert_eq!(m.len(), 12);
        assert_eq!(m[0].len(), 12);
    }
    #[test]
    fn test_rect_plate_mass_positive_diagonal() {
        let params = PlateParams::new(200e9, 0.3, 0.01);
        let elem = RectPlateElement::new(params, 0.5, 0.5);
        let m = elem.consistent_mass_matrix(7800.0);
        for (i, row) in m.iter().enumerate() {
            assert!(row[i] > 0.0, "m[{i}][{i}] should be positive");
        }
    }
    #[test]
    fn test_rect_plate_mass_scales_with_rho() {
        let params1 = PlateParams::new(200e9, 0.3, 0.01);
        let elem1 = RectPlateElement::new(params1, 0.5, 0.5);
        let params2 = PlateParams::new(200e9, 0.3, 0.01);
        let elem2 = RectPlateElement::new(params2, 0.5, 0.5);
        let m1 = elem1.consistent_mass_matrix(7800.0);
        let m2 = elem2.consistent_mass_matrix(15600.0);
        assert!((m2[0][0] / m1[0][0] - 2.0).abs() < 1e-12);
    }
    fn plate_element(a: f64, b: f64, nu: f64) -> RectPlateElement {
        RectPlateElement::new(PlateParams::new(200e9, nu, 0.01), a, b)
    }
    /// Nodal DOFs `(w, ∂w/∂y, −∂w/∂x)` of the field `w` with gradient `grad`.
    fn nodal_dofs(
        elem: &RectPlateElement,
        w: impl Fn(f64, f64) -> f64,
        grad: impl Fn(f64, f64) -> (f64, f64),
    ) -> Vec<f64> {
        elem.node_positions
            .iter()
            .flat_map(|&[x, y]| {
                let (w_x, w_y) = grad(x, y);
                [w(x, y), w_y, -w_x]
            })
            .collect()
    }
    fn mat_vec(k: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
        k.iter()
            .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
            .collect()
    }
    fn quad_form(k: &[Vec<f64>], v: &[f64]) -> f64 {
        mat_vec(k, v).iter().zip(v.iter()).map(|(a, b)| a * b).sum()
    }
    fn max_abs(k: &[Vec<f64>]) -> f64 {
        k.iter()
            .flat_map(|row| row.iter())
            .fold(0.0_f64, |m, v| m.max(v.abs()))
    }
    /// Eigenvalues (ascending) of a symmetric matrix by cyclic Jacobi rotations.
    fn symmetric_eigenvalues(matrix: &[Vec<f64>]) -> Vec<f64> {
        let n = matrix.len();
        let mut a: Vec<Vec<f64>> = matrix.to_vec();
        for _ in 0..100 {
            let off: f64 = (0..n)
                .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
                .map(|(i, j)| a[i][j] * a[i][j])
                .sum();
            let total: f64 = a.iter().flatten().map(|v| v * v).sum();
            if off <= 1e-30 * total {
                break;
            }
            for p in 0..n {
                for q in (p + 1)..n {
                    if a[p][q] == 0.0 {
                        continue;
                    }
                    let tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                    let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    for row in a.iter_mut() {
                        let (rp, rq) = (row[p], row[q]);
                        row[p] = c * rp - s * rq;
                        row[q] = s * rp + c * rq;
                    }
                    let row_p = a[p].clone();
                    let row_q = a[q].clone();
                    for k in 0..n {
                        a[p][k] = c * row_p[k] - s * row_q[k];
                        a[q][k] = s * row_p[k] + c * row_q[k];
                    }
                }
            }
        }
        let mut eigenvalues: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
        eigenvalues.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        eigenvalues
    }
    #[test]
    fn test_rect_plate_stiffness_is_symmetric_and_fully_coupled() {
        let k = plate_element(0.5, 0.25, 0.3).stiffness_matrix_12x12();
        let scale = max_abs(&k);
        let mut off_diagonal_nonzero = 0;
        for (i, row) in k.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                assert!(
                    (value - k[j][i]).abs() < 1e-12 * scale,
                    "asymmetric at ({i},{j})"
                );
                if i != j && value.abs() > 1e-9 * scale {
                    off_diagonal_nonzero += 1;
                }
            }
        }
        // A plate element couples every DOF to its neighbours: the matrix is
        // dense, not diagonal.
        assert!(
            off_diagonal_nonzero >= 100,
            "only {off_diagonal_nonzero} coupled pairs"
        );
    }
    #[test]
    fn test_rect_plate_stiffness_has_exactly_three_rigid_body_modes() {
        for (a, b, nu) in [
            (0.5, 0.5, 0.3),
            (0.5, 0.25, 0.3),
            (0.2, 0.9, 0.0),
            (1.0, 1.0, 0.45),
        ] {
            let k = plate_element(a, b, nu).stiffness_matrix_12x12();
            let scale = max_abs(&k);
            let scaled: Vec<Vec<f64>> = k
                .iter()
                .map(|row| row.iter().map(|v| v / scale).collect())
                .collect();
            let eig = symmetric_eigenvalues(&scaled);
            for (i, lambda) in eig.iter().take(3).enumerate() {
                assert!(
                    lambda.abs() < 1e-11,
                    "a={a} b={b}: eigenvalue {i} = {lambda}"
                );
            }
            assert!(
                eig[3] > 1e-4,
                "a={a} b={b}: fourth eigenvalue {} is not positive",
                eig[3]
            );
        }
    }
    #[test]
    fn test_rect_plate_rigid_motions_carry_no_force() {
        let elem = plate_element(0.5, 0.25, 0.3);
        let k = elem.stiffness_matrix_12x12();
        let scale = max_abs(&k);
        let modes = [
            nodal_dofs(&elem, |_, _| 1.0, |_, _| (0.0, 0.0)),
            nodal_dofs(&elem, |x, _| x, |_, _| (1.0, 0.0)),
            nodal_dofs(&elem, |_, y| y, |_, _| (0.0, 1.0)),
            nodal_dofs(&elem, |x, y| 2.0 - 3.0 * x + 5.0 * y, |_, _| (-3.0, 5.0)),
        ];
        for (m, mode) in modes.iter().enumerate() {
            let force = mat_vec(&k, mode);
            let norm_f = force.iter().map(|v| v * v).sum::<f64>().sqrt();
            let norm_d = mode.iter().map(|v| v * v).sum::<f64>().sqrt();
            assert!(
                norm_f < 1e-12 * scale * norm_d,
                "mode {m}: |K d| = {norm_f}"
            );
        }
    }
    #[test]
    fn test_rect_plate_strain_energy_of_constant_curvature_fields() {
        let (a, b, nu) = (0.5, 0.25, 0.3);
        let elem = plate_element(a, b, nu);
        let k = elem.stiffness_matrix_12x12();
        let d = elem.params.flexural_rigidity();
        let area = 4.0 * a * b;
        let energy = |dofs: &[f64]| 0.5 * quad_form(&k, dofs);
        // w = x²/2: κ = (1, 0, 0)
        let wxx = nodal_dofs(&elem, |x, _| 0.5 * x * x, |x, _| (x, 0.0));
        assert!((energy(&wxx) - 0.5 * d * area).abs() < 1e-12 * d * area);
        // w = y²/2: κ = (0, 1, 0)
        let wyy = nodal_dofs(&elem, |_, y| 0.5 * y * y, |_, y| (0.0, y));
        assert!((energy(&wyy) - 0.5 * d * area).abs() < 1e-12 * d * area);
        // w = (x² + y²)/2: κ = (1, 1, 0)  ->  ½ D (2 + 2ν) A
        let both = nodal_dofs(&elem, |x, y| 0.5 * (x * x + y * y), |x, y| (x, y));
        let expected = 0.5 * d * (2.0 + 2.0 * nu) * area;
        assert!((energy(&both) - expected).abs() < 1e-12 * d * area);
        // w = x y: κ = (0, 0, 2)  ->  ½ D (1 - ν)/2 · 4 · A
        let twist = nodal_dofs(&elem, |x, y| x * y, |x, y| (y, x));
        let expected = 0.5 * d * (1.0 - nu) / 2.0 * 4.0 * area;
        assert!((energy(&twist) - expected).abs() < 1e-12 * d * area);
    }
    #[test]
    fn test_rect_plate_stiffness_scales_with_element_size() {
        // Energy ~ D w²/s²: with scaled rotations (w, s θx, s θy) every entry
        // of the stiffness times s² is independent of the element size s.
        let kind_weight = |dof: usize, s: f64| if dof.is_multiple_of(3) { 1.0 } else { s };
        let (s1, s2) = (0.5, 1.3);
        let k1 = plate_element(s1, s1, 0.3).stiffness_matrix_12x12();
        let k2 = plate_element(s2, s2, 0.3).stiffness_matrix_12x12();
        let scale = max_abs(&k1) * s1 * s1;
        for i in 0..12 {
            for j in 0..12 {
                let m1 = k1[i][j] * s1 * s1 / (kind_weight(i, s1) * kind_weight(j, s1));
                let m2 = k2[i][j] * s2 * s2 / (kind_weight(i, s2) * kind_weight(j, s2));
                assert!(
                    (m1 - m2).abs() < 1e-12 * scale,
                    "entry ({i},{j}): {m1} vs {m2}"
                );
            }
        }
    }
    #[test]
    fn test_rect_plate_shape_functions_are_nodally_unit() {
        let elem = plate_element(0.5, 0.25, 0.3);
        let h = 1e-6;
        for (node, &[x, y]) in elem.node_positions.iter().enumerate() {
            let n = elem.shape_functions(x, y);
            let n_y_plus = elem.shape_functions(x, y + h);
            let n_y_minus = elem.shape_functions(x, y - h);
            let n_x_plus = elem.shape_functions(x + h, y);
            let n_x_minus = elem.shape_functions(x - h, y);
            for j in 0..12 {
                let value = if j == 3 * node { 1.0 } else { 0.0 };
                assert!(
                    (n[j] - value).abs() < 1e-12,
                    "N_{j}(node {node}) = {}",
                    n[j]
                );
                let d_dy = (n_y_plus[j] - n_y_minus[j]) / (2.0 * h);
                let minus_d_dx = -(n_x_plus[j] - n_x_minus[j]) / (2.0 * h);
                let theta_x = if j == 3 * node + 1 { 1.0 } else { 0.0 };
                let theta_y = if j == 3 * node + 2 { 1.0 } else { 0.0 };
                assert!(
                    (d_dy - theta_x).abs() < 1e-6,
                    "dN_{j}/dy at node {node} = {d_dy}"
                );
                assert!(
                    (minus_d_dx - theta_y).abs() < 1e-6,
                    "-dN_{j}/dx at node {node} = {minus_d_dx}"
                );
            }
        }
    }
    #[test]
    fn test_rect_plate_interpolates_every_polynomial_of_the_acm_space() {
        let elem = plate_element(0.5, 0.25, 0.3);
        let w = |x: f64, y: f64| {
            1.0 + 2.0 * x - y + 0.5 * x * x + x * y - 0.3 * y * y + 0.2 * x.powi(3)
                - 0.1 * x * x * y
                + 0.4 * x * y * y
                + 0.7 * y.powi(3)
                + 0.05 * x.powi(3) * y
                - 0.02 * x * y.powi(3)
        };
        let grad = |x: f64, y: f64| {
            (
                2.0 + x + y + 0.6 * x * x - 0.2 * x * y + 0.4 * y * y + 0.15 * x * x * y
                    - 0.02 * y.powi(3),
                -1.0 + x - 0.6 * y - 0.1 * x * x + 0.8 * x * y + 2.1 * y * y + 0.05 * x.powi(3)
                    - 0.06 * x * y * y,
            )
        };
        let dofs = nodal_dofs(&elem, w, grad);
        for &(x, y) in &[
            (0.0, 0.0),
            (0.2, -0.1),
            (-0.45, 0.2),
            (0.3, 0.24),
            (-0.5, -0.25),
        ] {
            let n = elem.shape_functions(x, y);
            let interpolated: f64 = n.iter().zip(dofs.iter()).map(|(a, b)| a * b).sum();
            assert!(
                (interpolated - w(x, y)).abs() < 1e-12,
                "({x}, {y}): {interpolated} vs {}",
                w(x, y)
            );
        }
    }
    #[test]
    fn test_rect_plate_consistent_mass_conserves_mass_and_inertia() {
        let (a, b, rho) = (0.5, 0.25, 7800.0);
        let elem = plate_element(a, b, 0.3);
        let m = elem.consistent_mass_matrix(rho);
        let h = elem.params.thickness;
        let scale = max_abs(&m);
        for (i, row) in m.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                assert!(
                    (value - m[j][i]).abs() < 1e-12 * scale,
                    "asymmetric at ({i},{j})"
                );
            }
        }
        let eig = symmetric_eigenvalues(
            &m.iter()
                .map(|r| r.iter().map(|v| v / scale).collect())
                .collect::<Vec<Vec<f64>>>(),
        );
        assert!(
            eig[0] > 1e-9,
            "mass matrix not positive definite: {}",
            eig[0]
        );
        let moment = |p: i32, q: i32| {
            let f = |len: f64, e: i32| {
                if e % 2 == 0 {
                    2.0 * len.powi(e + 1) / (e + 1) as f64
                } else {
                    0.0
                }
            };
            f(a, p) * f(b, q)
        };
        let fields: [(Vec<f64>, f64); 4] = [
            (
                nodal_dofs(&elem, |_, _| 1.0, |_, _| (0.0, 0.0)),
                rho * h * moment(0, 0),
            ),
            (
                nodal_dofs(&elem, |x, _| x, |_, _| (1.0, 0.0)),
                rho * h * moment(2, 0),
            ),
            (
                nodal_dofs(&elem, |_, y| y, |_, _| (0.0, 1.0)),
                rho * h * moment(0, 2),
            ),
            (
                nodal_dofs(&elem, |x, y| x * y, |x, y| (y, x)),
                rho * h * moment(2, 2),
            ),
        ];
        for (idx, (dofs, expected)) in fields.iter().enumerate() {
            let got = quad_form(&m, dofs);
            assert!(
                (got - expected).abs() < 1e-12 * expected,
                "field {idx}: {got} vs {expected}"
            );
        }
    }
    #[test]
    fn test_rect_plate_uniform_pressure_load_matches_numerical_integration() {
        let (a, b, q) = (0.5, 0.25, 3.0);
        let elem = plate_element(a, b, 0.3);
        let f = elem.uniform_pressure_load(q);
        let resultant: f64 = f.iter().step_by(3).sum();
        assert!(
            (resultant - q * 4.0 * a * b).abs() < 1e-13,
            "resultant = {resultant}"
        );
        // Independent check: 4-point Gauss-Legendre in each direction.
        let g = [
            (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
            (-0.339_981_043_584_856_26, 0.652_145_154_862_546_1),
            (0.339_981_043_584_856_26, 0.652_145_154_862_546_1),
            (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
        ];
        let mut integrated = [0.0_f64; 12];
        for &(xi, wx) in &g {
            for &(eta, wy) in &g {
                let n = elem.shape_functions(a * xi, b * eta);
                for (acc, nj) in integrated.iter_mut().zip(n.iter()) {
                    *acc += q * wx * wy * a * b * nj;
                }
            }
        }
        for j in 0..12 {
            assert!(
                (f[j] - integrated[j]).abs() < 1e-12,
                "load {j}: {} vs {}",
                f[j],
                integrated[j]
            );
        }
    }
    /// Centre deflection of a uniformly loaded unit square plate meshed with
    /// `nel × nel` ACM elements (`E = h = 1`, ν = 0.3, `q = 1`), simply
    /// supported (`w = 0` and zero tangential rotation on the edges) or clamped.
    fn square_plate_centre_deflection(nel: usize, clamped: bool) -> f64 {
        use crate::solvers::CholeskyDense;
        let params = PlateParams::new(1.0, 0.3, 1.0);
        let half = 0.5 / nel as f64;
        let elem = RectPlateElement::new(params, half, half);
        let ke = elem.stiffness_matrix_12x12();
        let fe = elem.uniform_pressure_load(1.0);
        let nn = nel + 1;
        let ndof = 3 * nn * nn;
        let node_id = |i: usize, j: usize| j * nn + i;
        let mut k = vec![0.0_f64; ndof * ndof];
        let mut f = vec![0.0_f64; ndof];
        for ey in 0..nel {
            for ex in 0..nel {
                let nodes = [
                    node_id(ex, ey),
                    node_id(ex + 1, ey),
                    node_id(ex + 1, ey + 1),
                    node_id(ex, ey + 1),
                ];
                let dofs: Vec<usize> = nodes
                    .iter()
                    .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
                    .collect();
                for r in 0..12 {
                    f[dofs[r]] += fe[r];
                    for c in 0..12 {
                        k[dofs[r] * ndof + dofs[c]] += ke[r][c];
                    }
                }
            }
        }
        let mut fixed = vec![false; ndof];
        for j in 0..nn {
            for i in 0..nn {
                let n = node_id(i, j);
                let on_x_edge = i == 0 || i == nel;
                let on_y_edge = j == 0 || j == nel;
                if on_x_edge || on_y_edge {
                    fixed[3 * n] = true;
                    if clamped {
                        fixed[3 * n + 1] = true;
                        fixed[3 * n + 2] = true;
                    } else {
                        if on_y_edge {
                            fixed[3 * n + 2] = true;
                        }
                        if on_x_edge {
                            fixed[3 * n + 1] = true;
                        }
                    }
                }
            }
        }
        let free: Vec<usize> = (0..ndof).filter(|&d| !fixed[d]).collect();
        let nf = free.len();
        let mut k_ff = vec![0.0_f64; nf * nf];
        let mut f_f = vec![0.0_f64; nf];
        for (r, &dr) in free.iter().enumerate() {
            f_f[r] = f[dr];
            for (c, &dc) in free.iter().enumerate() {
                k_ff[r * nf + c] = k[dr * ndof + dc];
            }
        }
        let chol = CholeskyDense::new(nf);
        let factor = chol.factorize(&k_ff).expect("constrained plate is SPD");
        let u = chol.solve(&factor, &f_f);
        let centre = 3 * node_id(nel / 2, nel / 2);
        let position = free
            .iter()
            .position(|&d| d == centre)
            .expect("centre w is free");
        u[position]
    }
    #[test]
    fn test_rect_plate_simply_supported_square_plate_matches_navier() {
        // w_max = 0.00406 q a⁴ / D  (Timoshenko), with D = E h³ / (12 (1 - ν²)).
        let d = PlateParams::new(1.0, 0.3, 1.0).flexural_rigidity();
        let reference = 0.00406 / d;
        let w8 = square_plate_centre_deflection(8, false);
        assert!(
            (w8 - reference).abs() < 0.03 * reference,
            "8x8 centre deflection {w8} vs reference {reference}"
        );
        let err4 = (square_plate_centre_deflection(4, false) - reference).abs();
        let err8 = (w8 - reference).abs();
        assert!(err8 < 0.4 * err4, "no mesh convergence: {err4} -> {err8}");
    }
    #[test]
    fn test_rect_plate_clamped_square_plate_matches_reference() {
        // w_max = 0.00126 q a⁴ / D for a clamped square plate.
        let d = PlateParams::new(1.0, 0.3, 1.0).flexural_rigidity();
        let reference = 0.00126 / d;
        let w8 = square_plate_centre_deflection(8, true);
        assert!(
            (w8 - reference).abs() < 0.05 * reference,
            "8x8 centre deflection {w8} vs reference {reference}"
        );
    }
    #[test]
    fn test_mindlin_params_shear_stiffness() {
        let p = MindlinParams::new(200e9, 0.3, 0.02, 5.0 / 6.0);
        let g = 200e9 / (2.0 * (1.0 + 0.3));
        let expected = (5.0 / 6.0) * g * 0.02;
        assert!((p.shear_stiffness() - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_mindlin_params_shear_scales_with_h() {
        let p1 = MindlinParams::new(200e9, 0.3, 0.01, 5.0 / 6.0);
        let p2 = MindlinParams::new(200e9, 0.3, 0.02, 5.0 / 6.0);
        assert!((p2.shear_stiffness() / p1.shear_stiffness() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_max_deflection_positive() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let w = max_deflection_uniform_load(1000.0, 1.0, 1.0, &p);
        assert!(w > 0.0);
    }
    #[test]
    fn test_max_deflection_scales_with_load() {
        let p1 = PlateParams::new(200e9, 0.3, 0.01);
        let p2 = PlateParams::new(200e9, 0.3, 0.01);
        let w1 = max_deflection_uniform_load(1000.0, 1.0, 1.0, &p1);
        let w2 = max_deflection_uniform_load(2000.0, 1.0, 1.0, &p2);
        assert!((w2 / w1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_max_deflection_formula() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let q = 500.0;
        let a = 2.0;
        let w = max_deflection_uniform_load(q, a, a, &p);
        let expected = 0.00406 * q * a * a * a * a / d;
        assert!((w - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_natural_frequency_mn_positive() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let f = natural_frequency_mn(1, 1, 1.0, 1.0, d, 7800.0, 0.01);
        assert!(f > 0.0);
    }
    #[test]
    fn test_natural_frequency_mn_increases_with_mode() {
        let p = PlateParams::new(200e9, 0.3, 0.01);
        let d = p.flexural_rigidity();
        let f11 = natural_frequency_mn(1, 1, 1.0, 1.0, d, 7800.0, 0.01);
        let f21 = natural_frequency_mn(2, 1, 1.0, 1.0, d, 7800.0, 0.01);
        assert!(f21 > f11);
    }
    #[test]
    fn test_shell_geometry_new() {
        let sg = ShellGeometry::new(1.0, 2.0, 0.01);
        assert!((sg.r1 - 1.0).abs() < 1e-15);
        assert!((sg.r2 - 2.0).abs() < 1e-15);
        assert!((sg.thickness - 0.01).abs() < 1e-15);
    }
    #[test]
    fn test_membrane_stiffness_formula() {
        let ns = membrane_stiffness(200e9, 0.3, 0.01);
        let expected = 200e9 * 0.01 / (1.0 - 0.3f64 * 0.3);
        assert!((ns - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_cylindrical_shell_buckling_positive() {
        let p = cylindrical_shell_buckling(200e9, 0.3, 1.0, 0.01);
        assert!(p > 0.0);
    }
    #[test]
    fn test_cylindrical_shell_buckling_formula() {
        let e = 200e9;
        let nu = 0.3;
        let r = 1.0;
        let h = 0.01;
        let p = cylindrical_shell_buckling(e, nu, r, h);
        let expected = e * h * h * h / (4.0 * (1.0 - nu * nu) * r * r * r);
        assert!((p - expected).abs() / expected < 1e-12);
        // t/R = 0.01, E = 200 GPa, nu = 0.3 gives 54.9 kPa, not the ~24 MPa of
        // the spherical-shell expression the old code evaluated.
        assert!((p - 54.9e3).abs() < 0.1e3, "p = {p}");
    }
    #[test]
    fn test_cylindrical_shell_buckling_agrees_with_cylindrical_shell_element() {
        let (e, nu, r, h) = (200e9, 0.3, 0.5, 0.01);
        let element = CylindricalShellElement::new(e, nu, 7800.0, h, 1.0, r);
        let from_element = element.critical_pressure_external();
        let from_fn = cylindrical_shell_buckling(e, nu, r, h);
        assert!((from_fn - from_element).abs() < 1e-12 * from_fn);
    }
    #[test]
    fn test_cylindrical_shell_buckling_finite_regimes() {
        let (e, nu, r, h) = (200e9, 0.3, 1.0, 0.01);
        let long_tube = cylindrical_shell_buckling(e, nu, r, h);
        let medium = cylindrical_shell_buckling_finite(e, nu, r, h, 8.0);
        let very_long = cylindrical_shell_buckling_finite(e, nu, r, h, 500.0);
        // Windenburg-Trilling for L/D = 4: 2.42 E (h/D)^2.5 / ((1-nu^2)^0.75 (L/D - 0.45 sqrt(h/D))).
        let t_over_d: f64 = h / (2.0 * r);
        let wt =
            2.42 * e * t_over_d.powf(2.5) / (0.91_f64.powf(0.75) * (4.0 - 0.45 * t_over_d.sqrt()));
        assert!((medium - wt).abs() < 1e-12 * wt, "medium = {medium}");
        assert!(medium > 3.0 * long_tube);
        assert!((very_long - long_tube).abs() < 1e-12 * long_tube);
        assert!(cylindrical_shell_buckling_finite(e, nu, r, h, 0.001).is_nan());
    }
    #[test]
    fn test_cylindrical_shell_buckling_scales_cubic() {
        let p1 = cylindrical_shell_buckling(200e9, 0.3, 1.0, 0.01);
        let p2 = cylindrical_shell_buckling(200e9, 0.3, 1.0, 0.02);
        assert!((p2 / p1 - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_spherical_shell_buckling_positive() {
        let p = spherical_shell_buckling(200e9, 0.3, 1.0, 0.01);
        assert!(p > 0.0);
    }
    #[test]
    fn test_spherical_shell_buckling_formula() {
        let e = 200e9;
        let nu = 0.3;
        let r = 1.0;
        let h = 0.01;
        let p = spherical_shell_buckling(e, nu, r, h);
        let expected = 2.0 * e / (3.0 * (1.0 - nu * nu)).sqrt() * (h / r).powi(2);
        assert!((p - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_spherical_shell_buckling_scales_quadratic() {
        let p1 = spherical_shell_buckling(200e9, 0.3, 1.0, 0.01);
        let p2 = spherical_shell_buckling(200e9, 0.3, 1.0, 0.02);
        assert!((p2 / p1 - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_aspect_ratio_square() {
        assert!((aspect_ratio(1.0, 1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_aspect_ratio_rectangle() {
        assert!((aspect_ratio(2.0, 1.0) - 2.0).abs() < 1e-15);
    }
    #[test]
    fn test_plate_concentrated_load_positive() {
        let w = plate_concentrated_load_deflection(1000.0, 1.0, 1e6);
        assert!(w > 0.0);
    }
    #[test]
    fn test_plate_concentrated_load_formula() {
        let p = 5000.0;
        let a = 2.0;
        let d = 1e7;
        let w = plate_concentrated_load_deflection(p, a, d);
        let expected = 0.01160 * p * a * a / d;
        assert!((w - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_laminate_single_layer_isotropic() {
        let ply = PlyLayer::new(200e9, 200e9, 77e9, 0.3, 0.01);
        let lam = CompositeLaminate::new(vec![ply], vec![0.0]);
        let abd = lam.abd_matrix();
        let e = 200e9;
        let nu = 0.3;
        let t = 0.01;
        let expected_a11 = e * t / (1.0 - nu * nu);
        assert!(
            (abd[0][0] - expected_a11).abs() / expected_a11 < 1e-6,
            "A11={} expected≈{}",
            abd[0][0],
            expected_a11
        );
    }
    #[test]
    fn test_laminate_abd_matrix_size() {
        let ply = PlyLayer::new(200e9, 20e9, 8e9, 0.3, 0.005);
        let lam = CompositeLaminate::new(vec![ply.clone(), ply], vec![0.0, 90.0_f64.to_radians()]);
        let abd = lam.abd_matrix();
        assert_eq!(abd.len(), 6);
        for row in &abd {
            assert_eq!(row.len(), 6);
        }
    }
    #[test]
    fn test_laminate_symmetric_cross_ply_b_zero() {
        let ply = PlyLayer::new(150e9, 10e9, 5e9, 0.3, 0.001);
        let thicknesses = vec![
            0.0,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI / 2.0,
            0.0,
        ];
        let lam = CompositeLaminate::new(
            vec![ply.clone(), ply.clone(), ply.clone(), ply],
            thicknesses,
        );
        let abd = lam.abd_matrix();
        for (i, row) in abd.iter().enumerate().skip(3) {
            for (j, &val) in row.iter().enumerate().take(3) {
                assert!(
                    val.abs() < 1.0,
                    "B[{i}][{j}]={} should be ~0 for symmetric laminate",
                    val
                );
            }
        }
    }
    #[test]
    fn test_laminate_in_plane_stiffness_positive_diagonal() {
        let ply = PlyLayer::new(200e9, 20e9, 8e9, 0.3, 0.003);
        let lam = CompositeLaminate::new(
            vec![ply.clone(), ply],
            vec![0.0, std::f64::consts::PI / 4.0],
        );
        let abd = lam.abd_matrix();
        assert!(
            abd[0][0] > 0.0 && abd[1][1] > 0.0 && abd[2][2] > 0.0,
            "A diagonal should be positive"
        );
    }
    #[test]
    fn test_ply_layer_q_matrix_isotropic_case() {
        let e = 200e9;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let ply = PlyLayer::new(e, e, g, nu, 0.01);
        let q = ply.q_matrix_global(0.0);
        let expected_q11 = e / (1.0 - nu * nu);
        assert!(
            (q[0][0] - expected_q11).abs() / expected_q11 < 1e-10,
            "Q11={} expected {}",
            q[0][0],
            expected_q11
        );
    }
    #[test]
    fn test_ply_layer_q_matrix_symmetry() {
        let ply = PlyLayer::new(150e9, 12e9, 5e9, 0.28, 0.005);
        let q = ply.q_matrix_global(std::f64::consts::PI / 6.0);
        for (i, row) in q.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - q[j][i]).abs() < 1.0,
                    "Q[{i}][{j}]={} vs Q[{j}][{i}]={}",
                    val,
                    q[j][i]
                );
            }
        }
    }
    #[test]
    fn test_sandwich_plate_flexural_rigidity() {
        let sw = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.02);
        let d = sw.effective_flexural_rigidity();
        assert!(d > 0.0, "Flexural rigidity should be positive");
    }
    #[test]
    fn test_sandwich_plate_shear_stiffness() {
        let sw = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.02);
        let ks = sw.transverse_shear_stiffness();
        assert!(ks > 0.0, "Shear stiffness should be positive");
    }
    #[test]
    fn test_sandwich_plate_total_thickness() {
        let sw = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.02);
        let t = sw.total_thickness();
        let expected = 2.0 * 0.001 + 0.02;
        assert!((t - expected).abs() < 1e-15);
    }
    #[test]
    fn test_sandwich_plate_face_governs_bending() {
        let sw1 = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.02);
        let sw2 = SandwichPlate::new(70e9, 0.33, 0.002, 1e6, 0.3, 0.02);
        assert!(sw2.effective_flexural_rigidity() > sw1.effective_flexural_rigidity());
    }
    #[test]
    fn test_sandwich_plate_core_shear() {
        let sw1 = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.01);
        let sw2 = SandwichPlate::new(70e9, 0.33, 0.001, 1e6, 0.3, 0.02);
        assert!(sw2.transverse_shear_stiffness() > sw1.transverse_shear_stiffness());
    }
    #[test]
    fn test_cylindrical_shell_natural_frequency() {
        let sh = CylindricalShellElement::new(200e9, 0.3, 7800.0, 0.01, 1.0, 0.5);
        let f = sh.ring_frequency();
        assert!(f > 0.0, "Ring frequency should be positive");
    }
    #[test]
    fn test_cylindrical_shell_critical_pressure() {
        let sh = CylindricalShellElement::new(200e9, 0.3, 7800.0, 0.01, 1.0, 0.5);
        let p = sh.critical_pressure_external();
        assert!(p > 0.0, "Critical pressure should be positive");
    }
    #[test]
    fn test_cylindrical_shell_membrane_hoop() {
        let sh = CylindricalShellElement::new(200e9, 0.3, 7800.0, 0.005, 1.0, 0.5);
        let p = 1e6;
        let sigma = sh.hoop_stress_internal_pressure(p);
        let expected = p * sh.radius / sh.thickness;
        assert!((sigma - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_cylindrical_shell_axial_stress() {
        let sh = CylindricalShellElement::new(200e9, 0.3, 7800.0, 0.005, 1.0, 0.5);
        let p = 1e6;
        let sigma = sh.axial_stress_closed_end(p);
        let expected = p * sh.radius / (2.0 * sh.thickness);
        assert!((sigma - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_cylindrical_shell_longitudinal_buckling() {
        let sh = CylindricalShellElement::new(200e9, 0.3, 7800.0, 0.01, 2.0, 0.5);
        let n_cr = sh.axial_buckling_load();
        assert!(n_cr > 0.0);
    }
    #[test]
    fn test_conical_shell_slant_length() {
        let cone = ConicalShell::new(200e9, 0.3, 7800.0, 0.01, 30.0_f64.to_radians(), 1.0);
        let s = cone.slant_length(0.5, 1.0);
        assert!(s > 0.0);
    }
    #[test]
    fn test_conical_shell_mean_radius() {
        let half_angle = 30.0_f64.to_radians();
        let cone = ConicalShell::new(200e9, 0.3, 7800.0, 0.01, half_angle, 0.5);
        let r_mean = cone.mean_radius_at(1.0);
        let expected = 0.5 + 1.0 * half_angle.sin();
        assert!((r_mean - expected).abs() < 1e-12);
    }
    #[test]
    fn test_conical_shell_meridional_stress() {
        let cone = ConicalShell::new(200e9, 0.3, 7800.0, 0.005, 30.0_f64.to_radians(), 0.5);
        let sigma = cone.meridional_stress_pressure(1e6, 0.5);
        assert!(sigma.is_finite() && sigma > 0.0);
    }
    #[test]
    fn test_shell_vibration_ring_mode() {
        let vib = ShellVibrationSolver::new(200e9, 0.3, 7800.0, 0.01, 0.5, 1.0);
        let omega = vib.ring_frequency_mode(2);
        assert!(omega > 0.0, "Ring mode frequency should be positive");
    }
    #[test]
    fn test_shell_vibration_frequency_increases_with_mode() {
        let vib = ShellVibrationSolver::new(200e9, 0.3, 7800.0, 0.01, 0.5, 1.0);
        let f2 = vib.ring_frequency_mode(2);
        let f3 = vib.ring_frequency_mode(3);
        let f4 = vib.ring_frequency_mode(4);
        assert!(f3 > f2 && f4 > f3, "f2={f2} f3={f3} f4={f4}");
    }
    #[test]
    fn test_shell_vibration_breathing_mode() {
        let vib = ShellVibrationSolver::new(200e9, 0.3, 7800.0, 0.01, 0.5, 1.0);
        let omega = vib.breathing_frequency();
        assert!(omega > 0.0);
    }
    #[test]
    fn test_flat_shell_stiffness_size() {
        let fs = FlatShellElement::new(200e9, 0.3, 7800.0, 0.01);
        let k = fs.stiffness_triangle(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        assert_eq!(k.len(), 18);
        for row in &k {
            assert_eq!(row.len(), 18);
        }
    }
    #[test]
    fn test_flat_shell_stiffness_nonzero_diagonal() {
        let fs = FlatShellElement::new(200e9, 0.3, 7800.0, 0.01);
        let k = fs.stiffness_triangle(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        let nonzero = k.iter().enumerate().any(|(i, row)| row[i].abs() > 0.0);
        assert!(nonzero, "Stiffness diagonal should have non-zero entries");
    }
    #[test]
    fn test_flat_shell_membrane_bending_decoupled() {
        let fs = FlatShellElement::new(70e9, 0.33, 2700.0, 0.002);
        let k = fs.stiffness_triangle(&[[0.0, 0.0], [0.5, 0.0], [0.0, 0.5]]);
        for row in &k {
            for &v in row {
                assert!(v.is_finite(), "Stiffness entry should be finite");
            }
        }
    }
    #[test]
    fn test_flat_shell_bending_stiffens_the_transverse_deflection_dofs() {
        // The previous "bending" matrix had only θ columns: every w column was
        // zero, so transverse deflection was a mechanism.
        let fs = FlatShellElement::new(200e9, 0.3, 7800.0, 0.01);
        let nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let kb = fs.bending_stiffness(&nodes);
        for w_dof in [0usize, 3, 6] {
            assert!(kb[w_dof][w_dof] > 0.0, "w DOF {w_dof} has no stiffness");
        }
        // A unit deflection of a single node, others held, needs a real force.
        let force_norm: f64 = (0..9).map(|r| kb[r][0] * kb[r][0]).sum::<f64>().sqrt();
        assert!(force_norm > 1.0e3, "reaction to w1 = 1: {force_norm}");
        // A pure rigid translation w = 1 still carries no force.
        let rigid = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        for row in &kb {
            let f: f64 = row.iter().zip(rigid.iter()).map(|(a, b)| a * b).sum();
            assert!(
                f.abs() < 1e-9 * kb[0][0].abs(),
                "rigid translation force {f}"
            );
        }
    }
    #[test]
    fn test_flat_shell_bending_has_rank_six() {
        let fs = FlatShellElement::new(70e9, 0.33, 2700.0, 0.004);
        let nodes = [[0.0, 0.0], [0.8, 0.1], [0.2, 0.6]];
        let kb = fs.bending_stiffness(&nodes);
        let rows: Vec<Vec<f64>> = kb.iter().map(|r| r.to_vec()).collect();
        let (eigenvalues, _) = crate::eigen::jacobi_eigen_dense(&rows, 9);
        let largest = eigenvalues[8];
        let null = eigenvalues
            .iter()
            .filter(|l| l.abs() < 1e-9 * largest)
            .count();
        let positive = eigenvalues.iter().filter(|&&l| l > 1e-4 * largest).count();
        assert_eq!((null, positive), (3, 6), "spectrum {eigenvalues:?}");
    }
    #[test]
    fn test_flat_shell_bending_block_matches_dkt_formulation() {
        let fs = FlatShellElement::new(200e9, 0.3, 7800.0, 0.01);
        let nodes = [[0.0, 0.0], [1.0, 0.2], [0.3, 0.9]];
        let plate = KirchhoffPlate::new(0.01, 200e9, 0.3, 7800.0);
        let expected =
            crate::shell::dkt::dkt_bending_stiffness(&nodes, &plate.constitutive_matrix());
        let direct = fs.bending_stiffness(&nodes);
        let full = fs.stiffness_triangle(&nodes);
        let bend_map = [2usize, 3, 4, 8, 9, 10, 14, 15, 16];
        for i in 0..9 {
            for j in 0..9 {
                assert_eq!(direct[i][j], expected[i][j]);
                assert_eq!(full[bend_map[i]][bend_map[j]], expected[i][j]);
            }
        }
        // The membrane DOFs (u, v) do not couple to the bending DOFs.
        for &m in &[0usize, 1, 6, 7, 12, 13] {
            for &b in &bend_map {
                assert_eq!(full[m][b], 0.0);
            }
        }
    }
    #[test]
    fn test_orthotropic_plate_rigidities() {
        let op = OrthotropicPlate::new(200e9, 20e9, 8e9, 0.3, 0.01);
        let (dx, dy, _dxy, dnu) = op.rigidities();
        assert!(
            dx > 0.0 && dy > 0.0 && dnu > 0.0,
            "All rigidities should be positive"
        );
    }
    #[test]
    fn test_orthotropic_plate_natural_frequency_ss() {
        let op = OrthotropicPlate::new(200e9, 20e9, 8e9, 0.3, 0.01);
        let omega = op.natural_frequency_ss(1, 1, 1.0, 1.0, 7800.0);
        assert!(omega > 0.0);
    }
    #[test]
    fn test_orthotropic_frequency_increases_with_mode() {
        let op = OrthotropicPlate::new(200e9, 20e9, 8e9, 0.3, 0.01);
        let f11 = op.natural_frequency_ss(1, 1, 1.0, 1.0, 7800.0);
        let f21 = op.natural_frequency_ss(2, 1, 1.0, 1.0, 7800.0);
        assert!(f21 > f11, "f21={f21} should be > f11={f11}");
    }
    #[test]
    fn test_orthotropic_plate_isotropic_limit() {
        let e = 200e9;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let op = OrthotropicPlate::new(e, e, g, nu, 0.01);
        let kp = KirchhoffPlate::new(0.01, e, nu, 7800.0);
        let f_op = op.natural_frequency_ss(1, 1, 1.0, 1.0, 7800.0);
        let f_kp = kp.natural_frequency_ss(1, 1, 1.0, 1.0);
        assert!((f_op - f_kp).abs() / f_kp < 0.01, "f_op={f_op} f_kp={f_kp}");
    }
    #[test]
    fn test_membrane_forces_output_length() {
        let sh = ShellElement::new(200e9, 0.3, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let disp = [0.0_f64; 6];
        let n_forces = sh.compute_membrane_forces(&nodes, &disp);
        assert_eq!(n_forces.len(), 3, "should return [Nx, Ny, Nxy]");
    }
    #[test]
    fn test_membrane_forces_zero_disp_zero_forces() {
        let sh = ShellElement::new(200e9, 0.3, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let disp = [0.0_f64; 6];
        let n_forces = sh.compute_membrane_forces(&nodes, &disp);
        for &f in &n_forces {
            assert!(f.abs() < 1e-6, "zero disp → zero membrane forces: {f}");
        }
    }
    #[test]
    fn test_membrane_forces_uniaxial_tension() {
        let e = 200e9;
        let nu = 0.3;
        let t = 0.01;
        let sh = ShellElement::new(e, nu, t);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let eps = 0.001;
        let disp = [0.0, 0.0, eps, 0.0, 0.0, 0.0_f64];
        let n_forces = sh.compute_membrane_forces(&nodes, &disp);
        assert!(
            n_forces[0].abs() > 0.0,
            "uniaxial: Nx should be non-zero: {}",
            n_forces[0]
        );
    }
    #[test]
    fn test_buckling_load_positive() {
        let sh = ShellElement::new(200e9, 0.3, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let n_cr = sh.compute_buckling_load(&nodes);
        assert!(n_cr > 0.0, "buckling load should be positive: {n_cr}");
    }
    #[test]
    fn test_buckling_load_thicker_is_higher() {
        let sh1 = ShellElement::new(200e9, 0.3, 0.005);
        let sh2 = ShellElement::new(200e9, 0.3, 0.010);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let n_cr1 = sh1.compute_buckling_load(&nodes);
        let n_cr2 = sh2.compute_buckling_load(&nodes);
        assert!(
            n_cr2 > n_cr1,
            "thicker shell → higher buckling load: {n_cr1} < {n_cr2}"
        );
    }
    #[test]
    fn test_buckling_load_stiffer_material_is_higher() {
        let sh_steel = ShellElement::new(200e9, 0.3, 0.01);
        let sh_alumin = ShellElement::new(70e9, 0.33, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let n_cr_s = sh_steel.compute_buckling_load(&nodes);
        let n_cr_a = sh_alumin.compute_buckling_load(&nodes);
        assert!(
            n_cr_s > n_cr_a,
            "steel (E=200GPa) → higher buckling than aluminium (E=70GPa): {n_cr_s} vs {n_cr_a}"
        );
    }
    #[test]
    fn test_large_deflection_correction_zero_deflection() {
        let sh = ShellElement::new(200e9, 0.3, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let disp = [0.0_f64; 9];
        let corr = sh.compute_large_deflection_correction(&nodes, &disp);
        for &v in &corr {
            assert!(
                v.abs() < 1e-10,
                "zero deflection → zero von Karman correction: {v}"
            );
        }
    }
    #[test]
    fn test_large_deflection_correction_scales_cubically() {
        let sh = ShellElement::new(200e9, 0.3, 0.01);
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let d1 = 1e-3_f64;
        let d2 = 2e-3_f64;
        let disp1 = [0.0_f64, 0.0, 0.0, d1, 0.0, 0.0, 0.0, 0.0, 0.0];
        let disp2 = [0.0_f64, 0.0, 0.0, d2, 0.0, 0.0, 0.0, 0.0, 0.0];
        let corr1 = sh.compute_large_deflection_correction(&nodes, &disp1);
        let corr2 = sh.compute_large_deflection_correction(&nodes, &disp2);
        let norm1: f64 = corr1.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm2: f64 = corr2.iter().map(|v| v * v).sum::<f64>().sqrt();
        let ratio = norm2 / norm1.max(1e-60);
        assert!(norm1 > 0.0, "correction should be non-zero: {norm1}");
        assert!(
            ratio > 4.0,
            "correction should grow super-linearly with deflection: ratio={ratio}"
        );
    }
    #[test]
    fn test_large_deflection_correction_finite() {
        let sh = ShellElement::new(70e9, 0.33, 0.002);
        let nodes = [[0.0_f64, 0.0], [0.5, 0.0], [0.0, 0.5]];
        let disp = [
            0.01_f64, 0.001, 0.001, 0.01, 0.001, 0.001, 0.01, 0.001, 0.001,
        ];
        let corr = sh.compute_large_deflection_correction(&nodes, &disp);
        for &v in &corr {
            assert!(v.is_finite(), "correction should be finite: {v}");
        }
    }
}
