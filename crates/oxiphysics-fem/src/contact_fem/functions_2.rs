//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::contact_fem::*;
    #[test]
    fn test_nts_gap_above_segment() {
        let slave = [1.0, 1.0, 0.0];
        let ma = [0.0, 0.0, 0.0];
        let mb = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let g = NodeToSegmentElement::gap(slave, ma, mb, normal);
        assert!(
            (g - 1.0).abs() < 1e-10,
            "Gap above segment should be 1.0: {g}"
        );
    }
    #[test]
    fn test_nts_gap_penetrating() {
        let slave = [1.0, -0.01, 0.0];
        let ma = [0.0, 0.0, 0.0];
        let mb = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let g = NodeToSegmentElement::gap(slave, ma, mb, normal);
        assert!(g < 0.0, "Gap for penetrating slave should be negative: {g}");
    }
    #[test]
    fn test_nts_penalty_force_above() {
        let elem = NodeToSegmentElement::new(1000.0, 500.0, 0.3);
        let slave = [1.0, 1.0, 0.0];
        let ma = [0.0, 0.0, 0.0];
        let mb = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let f = elem.penalty_force(slave, ma, mb, normal);
        assert!(f[1].abs() < 1e-10, "No force when above: {}", f[1]);
    }
    #[test]
    fn test_nts_penalty_force_penetrating() {
        let elem = NodeToSegmentElement::new(1000.0, 500.0, 0.3);
        let slave = [1.0, -0.01, 0.0];
        let ma = [0.0, 0.0, 0.0];
        let mb = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let f = elem.penalty_force(slave, ma, mb, normal);
        assert!(
            f[1] > 0.0,
            "Penalty force must be positive for penetration: {}",
            f[1]
        );
        assert!((f[1] - 10.0).abs() < 1e-9, "f = {} expected 10", f[1]);
    }
    #[test]
    fn test_nts_shape_functions_sum_to_one() {
        let slave = [0.7, 0.5, 0.0];
        let ma = [0.0, 0.0, 0.0];
        let mb = [1.0, 0.0, 0.0];
        let [n1, n2] = NodeToSegmentElement::shape_functions(slave, ma, mb);
        assert!(
            (n1 + n2 - 1.0).abs() < 1e-10,
            "Shape functions sum = {} != 1",
            n1 + n2
        );
    }
    #[test]
    fn test_nts_shape_functions_at_endpoints() {
        let ma = [0.0, 0.0, 0.0];
        let mb = [2.0, 0.0, 0.0];
        let [n1, n2] = NodeToSegmentElement::shape_functions([0.0, 0.0, 0.0], ma, mb);
        assert!((n1 - 1.0).abs() < 1e-10, "N_a at a = {n1}");
        assert!(n2.abs() < 1e-10, "N_b at a = {n2}");
        let [n1b, n2b] = NodeToSegmentElement::shape_functions([2.0, 0.0, 0.0], ma, mb);
        assert!(n1b.abs() < 1e-10, "N_a at b = {n1b}");
        assert!((n2b - 1.0).abs() < 1e-10, "N_b at b = {n2b}");
    }
    #[test]
    fn test_nts_contact_stiffness_in_contact() {
        let elem = NodeToSegmentElement::new(1000.0, 500.0, 0.3);
        let normal = [0.0, 1.0, 0.0];
        let ke = elem.contact_stiffness_matrix(normal, true);
        assert!((ke[1][1] - 1000.0).abs() < 1e-10, "k_yy = {}", ke[1][1]);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (ke[i][j] - ke[j][i]).abs() < 1e-10,
                    "Stiffness not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_nts_contact_stiffness_not_in_contact() {
        let elem = NodeToSegmentElement::new(1000.0, 500.0, 0.3);
        let normal = [0.0, 1.0, 0.0];
        let ke = elem.contact_stiffness_matrix(normal, false);
        for row in ke.iter() {
            for &v in row.iter() {
                assert_eq!(v, 0.0, "No stiffness when not in contact");
            }
        }
    }
    #[test]
    fn test_project_midpoint_of_segment() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let p = [1.0, 1.0, 0.0];
        let proj = project_point_to_segment(p, a, b);
        assert!((proj[0] - 1.0).abs() < 1e-10);
        assert!(proj[1].abs() < 1e-10);
    }
    #[test]
    fn test_project_beyond_end() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let p = [5.0, 1.0, 0.0];
        let proj = project_point_to_segment(p, a, b);
        assert!(
            (proj[0] - 1.0).abs() < 1e-10,
            "Should clamp to b: {}",
            proj[0]
        );
    }
    #[test]
    fn test_mortar_gap_integral_no_penetration() {
        let mortar = MortarContactElement {
            n_quad_pts: 3,
            penalty: 1000.0,
        };
        let slave_a = [0.0, 1.0, 0.0];
        let slave_b = [2.0, 1.0, 0.0];
        let master_a = [0.0, 0.0, 0.0];
        let master_b = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let g_int = mortar.mortar_gap_integral(slave_a, slave_b, master_a, master_b, normal);
        assert!(
            g_int > 0.0,
            "Integral should be positive for gap > 0: {g_int}"
        );
    }
    #[test]
    fn test_mortar_constraint_residual_no_contact() {
        let mortar = MortarContactElement {
            n_quad_pts: 2,
            penalty: 1000.0,
        };
        let slave_a = [0.0, 1.0, 0.0];
        let slave_b = [2.0, 1.0, 0.0];
        let master_a = [0.0, 0.0, 0.0];
        let master_b = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let r = mortar.constraint_residual(slave_a, slave_b, master_a, master_b, normal);
        assert!(r.abs() < 1e-10, "No contact → zero residual: {r}");
    }
    #[test]
    fn test_normal_pressure_recovery() {
        let csr = ContactStressRecovery { element_area: 1e-4 };
        let p = csr.normal_pressure(10.0);
        assert!((p - 100000.0).abs() < 0.1, "Pressure = {p}, expected 1e5");
    }
    #[test]
    fn test_shear_stress_recovery() {
        let csr = ContactStressRecovery { element_area: 1e-4 };
        let tau = csr.shear_stress(5.0);
        assert!((tau - 50000.0).abs() < 0.1, "Shear = {tau}, expected 5e4");
    }
    #[test]
    fn test_von_mises_contact_stress() {
        let csr = ContactStressRecovery { element_area: 1.0 };
        let vm = csr.von_mises_contact_stress(3.0, 4.0);
        let expected = (9.0_f64 + 48.0).sqrt();
        assert!(
            (vm - expected).abs() < 1e-10,
            "VM = {vm}, expected {expected}"
        );
    }
    #[test]
    fn test_max_shear_contact_stress() {
        let csr = ContactStressRecovery { element_area: 1.0 };
        let max_shear = csr.max_shear_stress(6.0, 4.0);
        assert!((max_shear - 5.0).abs() < 1e-10, "Max shear = {max_shear}");
    }
    #[test]
    fn test_nodal_stress_vector_direction() {
        let csr = ContactStressRecovery { element_area: 1.0 };
        let normal = [0.0, 0.0, 1.0];
        let sv = csr.nodal_stress_vector(10.0, normal);
        assert!(sv[0].abs() < 1e-10);
        assert!(sv[1].abs() < 1e-10);
        assert!((sv[2] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_traction_no_penetration() {
        let (t_n, t_t) = contact_traction(0.01, 0.0, 1000.0, 500.0, 0.3);
        assert_eq!(t_n, 0.0, "No traction for positive gap");
        assert_eq!(t_t, 0.0);
    }
    #[test]
    fn test_contact_traction_stick() {
        let (t_n, t_t) = contact_traction(-0.01, 0.001, 1000.0, 500.0, 0.3);
        assert!((t_n - 10.0).abs() < 1e-10, "t_n = {t_n}");
        assert!((t_t - 0.5).abs() < 1e-10, "t_t = {t_t}");
    }
    #[test]
    fn test_contact_traction_slip() {
        let (t_n, t_t) = contact_traction(-0.01, 1.0, 1000.0, 500.0, 0.3);
        assert!((t_n - 10.0).abs() < 1e-10, "t_n = {t_n}");
        assert!(
            (t_t - 3.0).abs() < 1e-10,
            "t_t = {t_t}, expected 3.0 (friction limit)"
        );
    }
    #[test]
    fn test_gap_to_quad_above_center() {
        let node = [0.5, 0.5, 1.0];
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [1.0, 1.0, 0.0];
        let v3 = [0.0, 1.0, 0.0];
        let (gap, normal) = gap_function_node_to_quad(node, v0, v1, v2, v3);
        assert!((gap - 1.0).abs() < 1e-10, "Gap should be 1.0: {gap}");
        assert!(
            (normal[2].abs() - 1.0).abs() < 1e-10,
            "Normal should be z: {}",
            normal[2]
        );
    }
    #[test]
    fn test_gap_to_quad_on_surface() {
        let node = [0.5, 0.5, 0.0];
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [1.0, 1.0, 0.0];
        let v3 = [0.0, 1.0, 0.0];
        let (gap, _) = gap_function_node_to_quad(node, v0, v1, v2, v3);
        assert!(gap.abs() < 1e-10, "On-surface gap should be 0: {gap}");
    }
    #[test]
    fn test_assemble_contact_stiffness_single_pair_in_contact() {
        let pairs = vec![(0, 1)];
        let gaps = vec![-0.01];
        let k = assemble_contact_stiffness(4, &pairs, &gaps, 1000.0);
        assert!((k[0] - 1000.0).abs() < 1e-10);
        assert!((k[5] - 1000.0).abs() < 1e-10);
        assert!((k[1] - (-1000.0)).abs() < 1e-10);
        assert!((k[4] - (-1000.0)).abs() < 1e-10);
    }
    #[test]
    fn test_assemble_contact_stiffness_no_contact() {
        let pairs = vec![(0, 1)];
        let gaps = vec![0.01];
        let k = assemble_contact_stiffness(4, &pairs, &gaps, 1000.0);
        for &v in k.iter() {
            assert_eq!(v, 0.0, "No stiffness for open gap");
        }
    }
    #[test]
    fn test_assemble_contact_force_penetrating() {
        let nodes = vec![2usize];
        let gaps = vec![-0.005];
        let f = assemble_contact_force(5, &nodes, &gaps, 1000.0);
        assert!(
            (f[2] - 5.0).abs() < 1e-10,
            "Contact force at node 2: {}",
            f[2]
        );
        assert_eq!(f[0], 0.0);
        assert_eq!(f[1], 0.0);
    }
    #[test]
    fn test_assemble_contact_force_no_contact() {
        let nodes = vec![0usize, 1];
        let gaps = vec![0.01, 0.02];
        let f = assemble_contact_force(3, &nodes, &gaps, 1000.0);
        for &v in f.iter() {
            assert_eq!(v, 0.0, "No force for open gaps");
        }
    }
    #[test]
    fn test_contact_pressure_field_at_node() {
        let contact_nodes = vec![[0.0, 0.0, 0.0_f64]];
        let nodal_forces = vec![10.0];
        let query = vec![[0.0, 0.0, 0.0_f64]];
        let field = contact_pressure_field(&query, &contact_nodes, &nodal_forces, 1e-4);
        assert!(field[0] > 0.0, "Pressure at contact node must be positive");
    }
    #[test]
    fn test_contact_pressure_field_empty() {
        let field = contact_pressure_field(&[[0.0, 0.0, 0.0]], &[], &[], 1e-4);
        assert_eq!(field[0], 0.0, "Empty contact nodes → zero pressure");
    }
    #[test]
    fn test_gauss_quad_1pt_integrates_constant() {
        let (pts, wts) = gauss_quad_1d(1);
        let integral: f64 = pts.iter().zip(wts.iter()).map(|(_, &w)| w * 1.0).sum();
        assert!(
            (integral - 2.0).abs() < 1e-10,
            "1-pt Gauss integral of 1: {integral}"
        );
    }
    #[test]
    fn test_gauss_quad_2pt_integrates_linear() {
        let (pts, wts) = gauss_quad_1d(2);
        let integral: f64 = pts.iter().zip(wts.iter()).map(|(&x, &w)| w * x).sum();
        assert!(
            integral.abs() < 1e-10,
            "Integral of x (odd) should be 0: {integral}"
        );
    }
    #[test]
    fn test_gauss_quad_3pt_integrates_quadratic() {
        let (pts, wts) = gauss_quad_1d(3);
        let integral: f64 = pts.iter().zip(wts.iter()).map(|(&x, &w)| w * x * x).sum();
        assert!(
            (integral - 2.0 / 3.0).abs() < 1e-10,
            "3-pt Gauss x² integral: {integral}"
        );
    }
    #[test]
    fn test_gauss_quad_weights_sum() {
        for n in [1, 2, 3, 4] {
            let (_, wts) = gauss_quad_1d(n);
            let s: f64 = wts.iter().sum();
            assert!((s - 2.0).abs() < 1e-6, "Gauss weights sum for n={n}: {s}");
        }
    }
}
