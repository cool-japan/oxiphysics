//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::PlateParams;

/// Navier series first-term maximum deflection for a uniformly-loaded
/// simply-supported rectangular plate.
///
/// w_max = 0.00406 · q · a⁴ / D  (square plate)
#[allow(dead_code)]
pub fn max_deflection_uniform_load(load: f64, a: f64, _b: f64, params: &PlateParams) -> f64 {
    let d = params.flexural_rigidity();
    0.00406 * load * a * a * a * a / d
}
/// Natural frequency of simply-supported rectangular plate, mode (m, n).
///
/// f_mn = π²/2 · (m²/a² + n²/b²) · √(D / (ρ·h))
#[allow(dead_code)]
pub fn natural_frequency_mn(m: u32, n: u32, a: f64, b: f64, d: f64, rho: f64, h: f64) -> f64 {
    use std::f64::consts::PI;
    let mf = m as f64;
    let nf = n as f64;
    PI * PI / 2.0 * (mf * mf / (a * a) + nf * nf / (b * b)) * (d / (rho * h)).sqrt()
}
/// Maximum membrane stiffness: N_max = E·h / (1−ν²).
#[allow(dead_code)]
pub fn membrane_stiffness(e: f64, nu: f64, h: f64) -> f64 {
    e * h / (1.0 - nu * nu)
}
/// Classical external pressure at buckling of a thin cylindrical shell.
///
/// p_cr = 2·E·(h/R)³ / (3·(1−ν²))
#[allow(dead_code)]
pub fn cylindrical_shell_buckling(e: f64, nu: f64, r: f64, h: f64) -> f64 {
    2.0 * e * (h / r).powi(3) / (3.0 * (1.0 - nu * nu))
}
/// Classical external pressure at buckling of a thin spherical shell.
///
/// p_cr = 2·E / √(3·(1−ν²)) · (h/R)²
#[allow(dead_code)]
pub fn spherical_shell_buckling(e: f64, nu: f64, r: f64, h: f64) -> f64 {
    2.0 * e / (3.0 * (1.0 - nu * nu)).sqrt() * (h / r).powi(2)
}
/// Aspect ratio of a plate: a/b.
#[allow(dead_code)]
pub fn aspect_ratio(a: f64, b: f64) -> f64 {
    a / b
}
/// Center deflection of a simply-supported square plate under a concentrated load.
///
/// w = 0.01160 · P · a² / D
#[allow(dead_code)]
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
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (dm[i][j] - dm[j][i]).abs() < 1e-6,
                    "D_plate not symmetric at [{i},{j}]: {} vs {}",
                    dm[i][j],
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
        for i in 0..6 {
            assert!(
                k[i][i] > 0.0,
                "stiffness diagonal k[{i}][{i}] = {} should be positive",
                k[i][i]
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
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (d[i][j] - d[j][i]).abs() < 1e-6,
                    "D not symmetric at [{i},{j}]: {} vs {}",
                    d[i][j],
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
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (db[i][j] - db[j][i]).abs() < 1e-6,
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
        for i in 0..12 {
            assert!(
                k[i][i] > 0.0,
                "k[{i}][{i}] should be positive, got {}",
                k[i][i]
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
        for i in 0..12 {
            assert!(m[i][i] > 0.0, "m[{i}][{i}] should be positive");
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
        let expected = 2.0 * e * (h / r).powi(3) / (3.0 * (1.0 - nu * nu));
        assert!((p - expected).abs() / expected < 1e-12);
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
        for i in 3..6 {
            for j in 0..3 {
                assert!(
                    abd[i][j].abs() < 1.0,
                    "B[{i}][{j}]={} should be ~0 for symmetric laminate",
                    abd[i][j]
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
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (q[i][j] - q[j][i]).abs() < 1.0,
                    "Q[{i}][{j}]={} vs Q[{j}][{i}]={}",
                    q[i][j],
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
