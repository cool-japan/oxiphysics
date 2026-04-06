//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use crate::vtu::types::*;
#[cfg(test)]
mod tests_vtu_expanded {
    use super::*;
    #[test]
    fn test_field_add() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let c = VtuFieldOps::add(&a, &b);
        assert_eq!(c, vec![5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_field_sub() {
        let a = vec![4.0, 5.0, 6.0];
        let b = vec![1.0, 2.0, 3.0];
        let c = VtuFieldOps::sub(&a, &b);
        assert_eq!(c, vec![3.0, 3.0, 3.0]);
    }
    #[test]
    fn test_field_mul() {
        let a = vec![2.0, 3.0];
        let b = vec![4.0, 5.0];
        let c = VtuFieldOps::mul(&a, &b);
        assert_eq!(c, vec![8.0, 15.0]);
    }
    #[test]
    fn test_field_scale() {
        let a = vec![1.0, 2.0, 3.0];
        let c = VtuFieldOps::scale(&a, 3.0);
        assert_eq!(c, vec![3.0, 6.0, 9.0]);
    }
    #[test]
    fn test_vector_magnitude() {
        let v = vec![3.0, 4.0, 0.0, 1.0, 0.0, 0.0];
        let m = VtuFieldOps::vector_magnitude(&v);
        assert_eq!(m.len(), 2);
        assert!((m[0] - 5.0).abs() < 1e-10, "mag={}", m[0]);
        assert!((m[1] - 1.0).abs() < 1e-10, "mag={}", m[1]);
    }
    #[test]
    fn test_normalize_vectors() {
        let v = vec![0.0, 3.0, 4.0];
        let n = VtuFieldOps::normalize_vectors(&v);
        assert_eq!(n.len(), 3);
        let mag: f64 = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_field_clamp() {
        let a = vec![-5.0, 0.5, 2.0];
        let c = VtuFieldOps::clamp(&a, 0.0, 1.0);
        assert_eq!(c, vec![0.0, 0.5, 1.0]);
    }
    #[test]
    fn test_dot_product() {
        let a = vec![1.0, 0.0, 0.0];
        let b = vec![0.0, 1.0, 0.0];
        let d = VtuFieldOps::dot_product(&a, &b);
        assert_eq!(d.len(), 1);
        assert!((d[0] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_cross_product() {
        let a = vec![1.0, 0.0, 0.0];
        let b = vec![0.0, 1.0, 0.0];
        let c = VtuFieldOps::cross_product(&a, &b);
        assert_eq!(c.len(), 3);
        assert!((c[0] - 0.0).abs() < 1e-10);
        assert!((c[1] - 0.0).abs() < 1e-10);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_field_stats() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (min_v, max_v, mean_v) = VtuFieldOps::stats(&a);
        assert!((min_v - 1.0).abs() < 1e-10);
        assert!((max_v - 5.0).abs() < 1e-10);
        assert!((mean_v - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_field_threshold() {
        let a = vec![0.5, 1.5, 2.5, 0.1];
        let t = VtuFieldOps::threshold(&a, 1.0);
        assert_eq!(t, vec![0.0, 1.0, 1.0, 0.0]);
    }
    #[test]
    fn test_cell_quality_triangle() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        w.add_triangle(0, 1, 2);
        let q = compute_cell_quality(&w);
        assert_eq!(q.len(), 1);
        assert!(q[0].aspect_ratio > 0.0 && q[0].aspect_ratio <= 1.0);
    }
    #[test]
    fn test_cell_quality_tet() {
        let mut w = VtuWriter::new();
        w.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        w.add_tet(0, 1, 2, 3);
        let q = compute_cell_quality(&w);
        assert_eq!(q.len(), 1);
        assert!(q[0].scaled_jacobian.abs() > 0.0);
    }
    #[test]
    fn test_cell_quality_multiple_cells() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 6]);
        w.add_triangle(0, 1, 2);
        w.add_triangle(3, 4, 5);
        let q = compute_cell_quality(&w);
        assert_eq!(q.len(), 2);
    }
    #[test]
    fn test_uniform_grid_to_vtu_1x1x1() {
        let w = uniform_grid_to_vtu(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        assert_eq!(w.num_points(), 8);
        assert_eq!(w.num_cells(), 1);
        assert_eq!(w.cell_types()[0], VTK_HEX);
    }
    #[test]
    fn test_uniform_grid_to_vtu_2x2x2() {
        let w = uniform_grid_to_vtu(2, 2, 2, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0);
        assert_eq!(w.num_points(), 27);
        assert_eq!(w.num_cells(), 8);
    }
    #[test]
    fn test_uniform_grid_bbox() {
        let w = uniform_grid_to_vtu(1, 1, 1, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
        let (min, max) = w.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-10);
        assert!((max[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_multi_piece_write() {
        let mut w0 = VtuWriter::new();
        w0.add_points(&[[0.0; 3]; 3]);
        w0.add_triangle(0, 1, 2);
        let mut w1 = VtuWriter::new();
        w1.add_points(&[[1.0; 3]; 4]);
        w1.add_quad(0, 1, 2, 3);
        let xml = VtuMultiPiece::write(&[&w0, &w1]);
        assert!(xml.contains("UnstructuredGrid"));
        assert_eq!(VtuMultiPiece::total_points(&[&w0, &w1]), 7);
        assert_eq!(VtuMultiPiece::total_cells(&[&w0, &w1]), 2);
    }
    #[test]
    fn test_point_to_cell_data_triangle() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 3]);
        w.add_triangle(0, 1, 2);
        let pd = vec![1.0, 2.0, 3.0];
        let cd = point_to_cell_data(&w, &pd);
        assert_eq!(cd.len(), 1);
        assert!((cd[0] - 2.0).abs() < 1e-10, "mean of [1,2,3] = 2");
    }
    #[test]
    fn test_cell_to_point_data_triangle() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 3]);
        w.add_triangle(0, 1, 2);
        let cd = vec![6.0];
        let pd = cell_to_point_data(&w, &cd);
        assert_eq!(pd.len(), 3);
        for &v in &pd {
            assert!((v - 6.0).abs() < 1e-10, "all nodes should get cell value");
        }
    }
    #[test]
    fn test_compute_bounding_box() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, 0.5, 0.0]];
        let (min, max) = compute_bounding_box(&pts).unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-10);
        assert!((max[1] - 2.0).abs() < 1e-10);
        assert!((max[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_compute_bounding_box_empty() {
        assert!(compute_bounding_box(&[]).is_none());
    }
    #[test]
    fn test_bbox_diagonal() {
        let d = bbox_diagonal([0.0; 3], [1.0, 1.0, 1.0]);
        assert!((d - f64::sqrt(3.0)).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_volume() {
        let v = bbox_volume([0.0; 3], [2.0, 3.0, 4.0]);
        assert!((v - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_threshold_cells() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 6]);
        w.add_triangle(0, 1, 2);
        w.add_triangle(3, 4, 5);
        let pd = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let above = threshold_cells(&w, &pd, 0.5);
        assert_eq!(above.len(), 1);
        assert_eq!(above[0], 1);
    }
    #[test]
    fn test_write_vtu_with_annotations() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 3]);
        w.add_triangle(0, 1, 2);
        let anns = vec![
            VtuAnnotation {
                key: "author".into(),
                value: "OxiPhysics".into(),
            },
            VtuAnnotation {
                key: "time".into(),
                value: "0.5".into(),
            },
        ];
        let xml = write_vtu_with_annotations(&w, &anns);
        assert!(xml.contains("<!-- VTU Annotations"));
        assert!(xml.contains("author: OxiPhysics"));
        assert!(xml.contains("time: 0.5"));
        assert!(xml.contains("<?xml version"));
    }
    #[test]
    fn test_higher_order_cell_num_nodes() {
        assert_eq!(higher_order_cell_num_nodes(VTK_QUADRATIC_TRIANGLE), 6);
        assert_eq!(higher_order_cell_num_nodes(VTK_QUADRATIC_QUAD), 8);
        assert_eq!(higher_order_cell_num_nodes(VTK_QUADRATIC_TET), 10);
        assert_eq!(higher_order_cell_num_nodes(VTK_QUADRATIC_HEX), 20);
        assert_eq!(higher_order_cell_num_nodes(VTK_TRIANGLE), 3);
    }
    #[test]
    fn test_estimate_gradient_uniform_field() {
        let w = uniform_grid_to_vtu(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let scalars = vec![5.0; w.num_points()];
        let grad = estimate_gradient(&w, &scalars);
        assert_eq!(grad.len(), w.num_points() * 3);
        for &g in &grad {
            assert!(
                g.abs() < 1e-10,
                "uniform field should have zero gradient, got {}",
                g
            );
        }
    }
    #[test]
    fn test_estimate_gradient_length() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 4]);
        w.add_tet(0, 1, 2, 3);
        let scalars = vec![1.0; 4];
        let grad = estimate_gradient(&w, &scalars);
        assert_eq!(grad.len(), 12);
    }
    #[test]
    fn test_write_cell_data_tensor_contains_name() {
        let tensors = vec![[1.0, 2.0, 3.0, 0.5, 0.3, 0.1]];
        let xml = VtuWriter::write_cell_data_tensor("stress", &tensors);
        assert!(xml.contains("Name=\"stress\""));
    }
    #[test]
    fn test_write_cell_data_tensor_six_components() {
        let tensors = vec![[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
        let xml = VtuWriter::write_cell_data_tensor("strain", &tensors);
        assert!(xml.contains("NumberOfComponents=\"6\""));
        assert!(xml.contains("1 0 0 0 0 0"));
    }
    #[test]
    fn test_write_cell_data_tensor_two_cells() {
        let tensors = vec![
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0, 0.1, 0.2, 0.3],
        ];
        let xml = VtuWriter::write_cell_data_tensor("T", &tensors);
        assert!(xml.contains("1 2 3 4 5 6 7 8 9"));
    }
    #[test]
    fn test_write_cell_data_tensor_empty() {
        let xml = VtuWriter::write_cell_data_tensor("empty_tensor", &[]);
        assert!(xml.contains("Name=\"empty_tensor\""));
        assert!(xml.contains("NumberOfComponents=\"6\""));
    }
    #[test]
    fn test_write_pvtu_parallel_xml_header() {
        let pvtu = VtuWriter::write_pvtu_parallel(&["part0.vtu", "part1.vtu"], &[], &[]);
        assert!(pvtu.starts_with("<?xml version=\"1.0\"?>"));
        assert!(pvtu.contains("PUnstructuredGrid"));
    }
    #[test]
    fn test_write_pvtu_parallel_piece_sources() {
        let pvtu = VtuWriter::write_pvtu_parallel(&["a.vtu", "b.vtu", "c.vtu"], &[], &[]);
        assert!(pvtu.contains("Source=\"a.vtu\""));
        assert!(pvtu.contains("Source=\"b.vtu\""));
        assert!(pvtu.contains("Source=\"c.vtu\""));
    }
    #[test]
    fn test_write_pvtu_parallel_point_data_names() {
        let pvtu = VtuWriter::write_pvtu_parallel(&["p0.vtu"], &["velocity", "pressure"], &[]);
        assert!(pvtu.contains("Name=\"velocity\""));
        assert!(pvtu.contains("Name=\"pressure\""));
        assert!(pvtu.contains("PPointData"));
    }
    #[test]
    fn test_write_pvtu_parallel_cell_data_names() {
        let pvtu = VtuWriter::write_pvtu_parallel(&["p0.vtu"], &[], &["stress", "strain"]);
        assert!(pvtu.contains("Name=\"stress\""));
        assert!(pvtu.contains("PCellData"));
    }
    #[test]
    fn test_read_cell_data_values() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 3]);
        w.add_triangle(0, 1, 2);
        w.add_cell_data_scalar("temperature", &[42.5]);
        let xml = w.to_xml();
        let vals = VtuReader::read_cell_data(&xml, "temperature").unwrap();
        assert_eq!(vals.len(), 1);
        assert!((vals[0] - 42.5).abs() < 1e-10);
    }
    #[test]
    fn test_read_cell_data_not_found() {
        let mut w = VtuWriter::new();
        w.add_points(&[[0.0; 3]; 3]);
        w.add_triangle(0, 1, 2);
        let xml = w.to_xml();
        let result = VtuReader::read_cell_data(&xml, "nonexistent");
        assert!(result.is_err());
    }
}
