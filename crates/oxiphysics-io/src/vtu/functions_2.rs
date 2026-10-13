//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_vtu_expanded {
    use crate::vtu::*;
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
        // Edges 1, sqrt(2), 1 and a worst corner angle of 45 degrees.
        assert!((q[0].aspect_ratio - 1.0 / 2.0_f64.sqrt()).abs() < 1e-14);
        assert!((q[0].scaled_jacobian - (2.0_f64 / 3.0).sqrt()).abs() < 1e-14);
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
        // Corner tetrahedron: sqrt(2) * 1 / max(1, 2, 2, 2) = sqrt(2) / 2.
        assert!((q[0].scaled_jacobian - 0.5 * 2.0_f64.sqrt()).abs() < 1e-14);
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
    fn quality_of(points: &[[f64; 3]], cell: &[usize], cell_type: u8) -> VtuCellQuality {
        let mut w = VtuWriter::new();
        w.add_points(points);
        w.add_cell(cell.to_vec(), cell_type);
        let q = compute_cell_quality(&w);
        assert_eq!(q.len(), 1);
        q[0].clone()
    }
    const UNIT_CUBE: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    const HEX_NODES: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];
    #[test]
    fn test_cell_quality_unit_cube_is_perfect() {
        // The old cyclic-node "edges" gave aspect 1/sqrt(2) and a scaled
        // Jacobian of 0 for a unit hexahedron.
        let q = quality_of(&UNIT_CUBE, &HEX_NODES, VTK_HEX);
        assert_eq!(q.cell_index, 0);
        assert_eq!(q.cell_type, VTK_HEX);
        assert!((q.min_edge - 1.0).abs() < 1e-14);
        assert!((q.max_edge - 1.0).abs() < 1e-14);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_stretched_cube_aspect_ratio() {
        // Orthogonal box 1 x 1 x 3: perfect Jacobian, edges 1 and 3.
        let mut pts = UNIT_CUBE;
        for p in pts.iter_mut() {
            p[2] *= 3.0;
        }
        let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
        assert!((q.min_edge - 1.0).abs() < 1e-14);
        assert!((q.max_edge - 3.0).abs() < 1e-14);
        assert!((q.aspect_ratio - 1.0 / 3.0).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_cell_quality_inverted_hex_is_negative() {
        let mut pts = UNIT_CUBE;
        for p in pts[4..].iter_mut() {
            p[2] = -1.0;
        }
        let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
        assert!((q.scaled_jacobian + 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_skewed_hex_uses_worst_corner() {
        // Pull the (1,1,1) corner out of the cube: the minimum over the eight
        // corner Jacobians is the reference value from an independent numpy
        // evaluation of the same definition.
        let mut pts = UNIT_CUBE;
        pts[6] = [1.8, 1.8, 1.0];
        let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
        assert!(
            (q.scaled_jacobian - 0.443_786_614_480_555_9).abs() < 1e-12,
            "{q:?}"
        );
    }
    #[test]
    fn test_cell_quality_is_scale_invariant_at_atomic_lengths() {
        // Edge products of 1e-30 must not be mistaken for degenerate cells.
        let mut pts = UNIT_CUBE;
        for p in pts.iter_mut() {
            for c in p.iter_mut() {
                *c *= 1.0e-10;
            }
        }
        let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-12, "{q:?}");
        assert!((q.max_edge - 1.0e-10).abs() < 1e-24);
    }
    #[test]
    fn test_cell_quality_voxel_matches_hex_with_reordered_nodes() {
        // VTK_VOXEL lists the nodes with x varying fastest.
        let voxel = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [0.0, 2.0, 2.0],
            [2.0, 2.0, 2.0],
        ];
        let q = quality_of(&voxel, &HEX_NODES, VTK_VOXEL);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.max_edge - 2.0).abs() < 1e-14);
    }
    #[test]
    fn test_cell_quality_unit_square_is_perfect() {
        let sq = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let q = quality_of(&sq, &[0, 1, 2, 3], VTK_QUAD);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_quad_in_arbitrary_plane() {
        // The same unit square lying in the yz-plane (counter-clockwise about +x).
        let sq = [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
        ];
        let q = quality_of(&sq, &[0, 1, 2, 3], VTK_QUAD);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_concave_quad_is_negative() {
        // Corner 2 is pushed into the quadrilateral: a reflex angle.
        let bad = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.2, 0.2, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let q = quality_of(&bad, &[0, 1, 2, 3], VTK_QUAD);
        assert!(
            (q.scaled_jacobian - (-0.882_352_941_176_470_3)).abs() < 1e-12,
            "{q:?}"
        );
    }
    #[test]
    fn test_cell_quality_equilateral_triangle_is_perfect() {
        let h = 3.0_f64.sqrt() / 2.0;
        let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]];
        let q = quality_of(&tri, &[0, 1, 2], VTK_TRIANGLE);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_right_triangle_values() {
        let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let q = quality_of(&tri, &[0, 1, 2], VTK_TRIANGLE);
        assert!((q.min_edge - 1.0).abs() < 1e-14);
        assert!((q.max_edge - 2.0_f64.sqrt()).abs() < 1e-14);
        assert!((q.aspect_ratio - 1.0 / 2.0_f64.sqrt()).abs() < 1e-14);
        // 2/sqrt(3) * sin(45 degrees)
        assert!((q.scaled_jacobian - (2.0_f64 / 3.0).sqrt()).abs() < 1e-14);
    }
    #[test]
    fn test_cell_quality_triangle_in_space_is_independent_of_orientation() {
        let h = 3.0_f64.sqrt() / 2.0;
        // Equilateral triangle rotated about the x-axis out of the xy-plane.
        let base = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]];
        let tilt = |p: [f64; 3]| [p[0], p[1] * 0.6 - p[2] * 0.8, p[1] * 0.8 + p[2] * 0.6];
        let tri = [tilt(base[0]), tilt(base[1]), tilt(base[2])];
        let q = quality_of(&tri, &[0, 1, 2], VTK_TRIANGLE);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-13, "{q:?}");
    }
    const REGULAR_TET: [[f64; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    #[test]
    fn test_cell_quality_regular_tet_is_perfect() {
        let q = quality_of(&REGULAR_TET, &[0, 1, 2, 3], VTK_TET);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.min_edge - 8.0_f64.sqrt()).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_inverted_tet_is_negative() {
        let q = quality_of(&REGULAR_TET, &[0, 2, 1, 3], VTK_TET);
        assert!((q.scaled_jacobian + 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_corner_tet_is_half_sqrt_two() {
        let corner = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let q = quality_of(&corner, &[0, 1, 2, 3], VTK_TET);
        assert!((q.aspect_ratio - 1.0 / 2.0_f64.sqrt()).abs() < 1e-14);
        assert!(
            (q.scaled_jacobian - 0.5 * 2.0_f64.sqrt()).abs() < 1e-14,
            "{q:?}"
        );
    }
    #[test]
    fn test_cell_quality_ideal_pyramid_is_perfect() {
        // Square base of side 1 with all eight edges equal.
        let apex = 1.0 / 2.0_f64.sqrt();
        let pyramid = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, apex],
        ];
        let q = quality_of(&pyramid, &[0, 1, 2, 3, 4], VTK_PYRAMID);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-13, "{q:?}");
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-13, "{q:?}");
        let mut flipped = pyramid;
        flipped[4][2] = -apex;
        let q = quality_of(&flipped, &[0, 1, 2, 3, 4], VTK_PYRAMID);
        assert!((q.scaled_jacobian + 1.0).abs() < 1e-13, "{q:?}");
    }
    #[test]
    fn test_cell_quality_wedge_in_vtk_orientation() {
        // VTK orders the base triangle so that its right-hand normal points
        // away from the top triangle.  Equilateral prism, all edges 1.
        let h = 3.0_f64.sqrt() / 2.0;
        let wedge = [
            [0.0, 0.0, 0.0],
            [0.5, h, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, h, 1.0],
            [1.0, 0.0, 1.0],
        ];
        let nodes = [0, 1, 2, 3, 4, 5];
        let q = quality_of(&wedge, &nodes, VTK_WEDGE);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-13);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-13, "{q:?}");
        // Mirroring the top triangle below the base inverts the element.
        let mut inverted = wedge;
        for p in inverted[3..].iter_mut() {
            p[2] = -1.0;
        }
        let q = quality_of(&inverted, &nodes, VTK_WEDGE);
        assert!((q.scaled_jacobian + 1.0).abs() < 1e-13, "{q:?}");
        // A right isosceles base: worst corner angle is 45 degrees.
        let right = [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
        ];
        let q = quality_of(&right, &nodes, VTK_WEDGE);
        assert!(
            (q.scaled_jacobian - (2.0_f64 / 3.0).sqrt()).abs() < 1e-13,
            "{q:?}"
        );
        assert!((q.aspect_ratio - 1.0 / 2.0_f64.sqrt()).abs() < 1e-13);
    }
    #[test]
    fn test_cell_quality_line_cell() {
        let q = quality_of(&[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]], &[0, 1], VTK_LINE);
        assert!((q.min_edge - 5.0).abs() < 1e-14);
        assert!((q.max_edge - 5.0).abs() < 1e-14);
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14);
        let zero = quality_of(&[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], &[0, 1], VTK_LINE);
        assert_eq!(zero.scaled_jacobian, 0.0);
        assert_eq!(zero.aspect_ratio, 0.0);
    }
    #[test]
    fn test_cell_quality_degenerate_tet_is_zero() {
        let q = quality_of(&[[1.0, 2.0, 3.0]; 4], &[0, 1, 2, 3], VTK_TET);
        assert_eq!(q.min_edge, 0.0);
        assert_eq!(q.max_edge, 0.0);
        assert_eq!(q.aspect_ratio, 0.0);
        assert_eq!(q.scaled_jacobian, 0.0);
    }
    #[test]
    fn test_cell_quality_flat_hex_has_zero_jacobian() {
        let mut pts = UNIT_CUBE;
        for p in pts.iter_mut() {
            p[2] = 0.0;
        }
        let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
        assert_eq!(q.scaled_jacobian, 0.0, "{q:?}");
    }
    #[test]
    fn test_cell_quality_quadratic_cells_use_corner_nodes() {
        // Quadratic tetrahedron: 4 corners followed by 6 mid-edge nodes.
        let mut pts = REGULAR_TET.to_vec();
        for &(a, b) in &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
            let mid = [
                0.5 * (REGULAR_TET[a][0] + REGULAR_TET[b][0]),
                0.5 * (REGULAR_TET[a][1] + REGULAR_TET[b][1]),
                0.5 * (REGULAR_TET[a][2] + REGULAR_TET[b][2]),
            ];
            pts.push(mid);
        }
        let nodes: Vec<usize> = (0..10).collect();
        let q = quality_of(&pts, &nodes, VTK_QUADRATIC_TET);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
        assert!((q.aspect_ratio - 1.0).abs() < 1e-14);
        assert_eq!(q.cell_type, VTK_QUADRATIC_TET);
    }
    #[test]
    fn test_cell_quality_unmeasurable_cells_report_zeros() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // Vertex and polygon cells are not evaluated.
        for (cell, ct) in [(vec![0], VTK_VERTEX), (vec![0, 1, 2], VTK_POLYGON)] {
            let q = quality_of(&pts, &cell, ct);
            assert_eq!(q.aspect_ratio, 0.0);
            assert_eq!(q.min_edge, 0.0);
            assert_eq!(q.max_edge, 0.0);
            assert_eq!(q.scaled_jacobian, 0.0);
            assert_eq!(q.cell_type, ct);
        }
        // Wrong node count for the type.
        let q = quality_of(&pts, &[0, 1], VTK_TRIANGLE);
        assert_eq!(q.scaled_jacobian, 0.0);
        // Connectivity entry outside the point list.
        let q = quality_of(&pts, &[0, 1, 7], VTK_TRIANGLE);
        assert_eq!(q.scaled_jacobian, 0.0);
        assert_eq!(q.max_edge, 0.0);
    }
    #[test]
    fn test_cell_quality_non_finite_coordinates_report_zeros() {
        // A NaN edge length used to leave `min_edge` at +infinity and the
        // Jacobian at +infinity; such a cell is simply unmeasurable.
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut pts = UNIT_CUBE;
            pts[6][1] = bad;
            let q = quality_of(&pts, &HEX_NODES, VTK_HEX);
            assert_eq!(q.min_edge, 0.0, "{bad}");
            assert_eq!(q.max_edge, 0.0, "{bad}");
            assert_eq!(q.aspect_ratio, 0.0, "{bad}");
            assert_eq!(q.scaled_jacobian, 0.0, "{bad}");
            assert_eq!(q.cell_type, VTK_HEX);
        }
        // Only the corner nodes of a quadratic cell are measured, so a bad
        // mid-edge node does not disqualify it.
        let mut pts = REGULAR_TET.to_vec();
        for &(a, b) in &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
            pts.push([
                0.5 * (REGULAR_TET[a][0] + REGULAR_TET[b][0]),
                0.5 * (REGULAR_TET[a][1] + REGULAR_TET[b][1]),
                0.5 * (REGULAR_TET[a][2] + REGULAR_TET[b][2]),
            ]);
        }
        pts[7][0] = f64::NAN;
        let nodes: Vec<usize> = (0..10).collect();
        let q = quality_of(&pts, &nodes, VTK_QUADRATIC_TET);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-14, "{q:?}");
    }
    #[test]
    fn test_cell_quality_mixed_mesh_keeps_cell_order() {
        let mut w = VtuWriter::new();
        w.add_points(&UNIT_CUBE);
        w.add_points(&[[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]]);
        w.add_hex(0, 1, 2, 3, 4, 5, 6, 7);
        w.add_triangle(8, 9, 10);
        w.add_cell(vec![0, 1], VTK_LINE);
        let q = compute_cell_quality(&w);
        assert_eq!(q.len(), 3);
        assert_eq!((q[0].cell_index, q[0].cell_type), (0, VTK_HEX));
        assert_eq!((q[1].cell_index, q[1].cell_type), (1, VTK_TRIANGLE));
        assert_eq!((q[2].cell_index, q[2].cell_type), (2, VTK_LINE));
        assert!((q[0].scaled_jacobian - 1.0).abs() < 1e-14);
        assert!((q[1].scaled_jacobian - (2.0_f64 / 3.0).sqrt()).abs() < 1e-14);
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
    fn linear_field(w: &VtuWriter, g: [f64; 3], offset: f64) -> Vec<f64> {
        w.points()
            .iter()
            .map(|p| g[0] * p[0] + g[1] * p[1] + g[2] * p[2] + offset)
            .collect()
    }
    fn assert_gradient_everywhere(grad: &[f64], expected: [f64; 3], tol: f64) {
        assert_eq!(grad.len() % 3, 0);
        for (node, g) in grad.chunks_exact(3).enumerate() {
            for k in 0..3 {
                assert!(
                    (g[k] - expected[k]).abs() < tol,
                    "node {node} component {k}: got {} expected {}",
                    g[k],
                    expected[k]
                );
            }
        }
    }
    #[test]
    fn test_estimate_gradient_unit_hex_reproduces_unit_slope() {
        // s = x on a unit hexahedron: the old edge-projection estimator gave
        // 1/3 at node 0.
        let w = uniform_grid_to_vtu(1, 1, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let scalars = linear_field(&w, [1.0, 0.0, 0.0], 0.0);
        let grad = estimate_gradient(&w, &scalars);
        assert_eq!(grad.len(), 24);
        assert_gradient_everywhere(&grad, [1.0, 0.0, 0.0], 1e-12);
    }
    #[test]
    fn test_estimate_gradient_recovers_linear_field_on_3x3x3_grid() {
        let w = uniform_grid_to_vtu(3, 3, 3, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        assert_eq!(w.num_points(), 64);
        let scalars = linear_field(&w, [2.0, -1.0, 3.0], 0.0);
        let grad = estimate_gradient(&w, &scalars);
        assert_gradient_everywhere(&grad, [2.0, -1.0, 3.0], 1e-10);
    }
    #[test]
    fn test_estimate_gradient_is_independent_of_offset_and_scale() {
        let w = uniform_grid_to_vtu(2, 3, 2, -5.0, 7.0, 100.0, 100.5, 0.0, 1.0e-3);
        let scalars = linear_field(&w, [0.5, -40.0, 1.0e3], 1.0e6);
        let grad = estimate_gradient(&w, &scalars);
        assert_gradient_everywhere(&grad, [0.5, -40.0, 1.0e3], 1e-6);
    }
    #[test]
    fn test_estimate_gradient_on_unstructured_tets() {
        // Two tetrahedra sharing a face, at irregular positions.
        let mut w = VtuWriter::new();
        w.add_points(&[
            [0.0, 0.0, 0.0],
            [1.3, 0.1, -0.2],
            [0.2, 0.9, 0.3],
            [0.4, 0.2, 1.1],
            [1.5, 1.4, 1.2],
        ]);
        w.add_tet(0, 1, 2, 3);
        w.add_tet(1, 2, 3, 4);
        let scalars = linear_field(&w, [-1.5, 0.25, 4.0], 2.0);
        let grad = estimate_gradient(&w, &scalars);
        assert_gradient_everywhere(&grad, [-1.5, 0.25, 4.0], 1e-12);
    }
    #[test]
    fn test_estimate_gradient_planar_mesh_gives_in_plane_gradient() {
        // A flat quadrilateral grid in z = 0: the normal matrix is singular and
        // the pseudo-inverse must keep the in-plane gradient, not zero it.
        let w = uniform_grid_to_vtu(3, 2, 1, 0.0, 3.0, 0.0, 2.0, 0.0, 0.0);
        let scalars = linear_field(&w, [3.0, 2.0, 0.0], 0.0);
        let grad = estimate_gradient(&w, &scalars);
        assert_gradient_everywhere(&grad, [3.0, 2.0, 0.0], 1e-10);
    }
    #[test]
    fn test_estimate_gradient_collinear_mesh() {
        let mut w = VtuWriter::new();
        w.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ]);
        w.add_cell(vec![0, 1], VTK_LINE);
        w.add_cell(vec![1, 2], VTK_LINE);
        w.add_cell(vec![2, 3], VTK_LINE);
        let scalars = linear_field(&w, [5.0, 0.0, 0.0], 0.0);
        let grad = estimate_gradient(&w, &scalars);
        assert_gradient_everywhere(&grad, [5.0, 0.0, 0.0], 1e-12);
    }
    #[test]
    fn test_estimate_gradient_exact_for_quadratic_on_symmetric_stencil() {
        // Interior node of a 2x2x2 grid sees a point-symmetric stencil, so the
        // even part of x^2 + y^2 cancels and the centre gradient is exact.
        let w = uniform_grid_to_vtu(2, 2, 2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let scalars: Vec<f64> = w
            .points()
            .iter()
            .map(|p| p[0] * p[0] + p[1] * p[1])
            .collect();
        let grad = estimate_gradient(&w, &scalars);
        let centre = w
            .points()
            .iter()
            .position(|p| p == &[0.5, 0.5, 0.5])
            .unwrap();
        assert!((grad[centre * 3] - 1.0).abs() < 1e-12);
        assert!((grad[centre * 3 + 1] - 1.0).abs() < 1e-12);
        assert!(grad[centre * 3 + 2].abs() < 1e-12);
    }
    #[test]
    fn test_estimate_gradient_isolated_and_unvalued_nodes_are_zero() {
        let mut w = VtuWriter::new();
        w.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [9.0, 9.0, 9.0],
        ]);
        w.add_triangle(0, 1, 2);
        // Node 3 belongs to no cell; scalars are also one entry short of the
        // number of points.
        let grad = estimate_gradient(&w, &[0.0, 1.0, 2.0]);
        assert_eq!(grad.len(), 12);
        assert_eq!(&grad[9..], &[0.0, 0.0, 0.0]);
        assert!((grad[0] - 1.0).abs() < 1e-12);
        assert!((grad[1] - 2.0).abs() < 1e-12);
        assert!(grad[2].abs() < 1e-12);
    }
    #[test]
    fn test_estimate_gradient_empty_mesh() {
        let w = VtuWriter::new();
        assert!(estimate_gradient(&w, &[]).is_empty());
    }
    #[test]
    fn test_estimate_gradient_ignores_coincident_neighbours() {
        let mut w = VtuWriter::new();
        w.add_points(&[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        w.add_cell(vec![0, 1, 2, 3, 4], VTK_POLY_VERTEX);
        let scalars = linear_field(&w, [1.0, 2.0, 3.0], 0.0);
        let grad = estimate_gradient(&w, &scalars);
        assert!(grad.iter().all(|g| g.is_finite()));
        assert!((grad[0] - 1.0).abs() < 1e-12);
        assert!((grad[1] - 2.0).abs() < 1e-12);
        assert!((grad[2] - 3.0).abs() < 1e-12);
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
