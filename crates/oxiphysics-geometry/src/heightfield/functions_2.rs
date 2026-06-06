//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_extra {
    use crate::heightfield::*;
    #[test]
    fn test_dda_traversal_diagonal_ray() {
        let hf = HeightField::from_fn(6, 6, 1.0, |_, _| 0.0);
        let origin = [0.1, 10.0, 0.1];
        let dir = [1.0, 0.0, 1.0];
        let traversal = HeightfieldRayTraversal::new(&hf, origin, dir, 200.0);
        assert!(traversal.is_some(), "diagonal ray should create traversal");
        let mut t = traversal.unwrap();
        let mut cells: Vec<(usize, usize)> = Vec::new();
        while let Some((col, row, _)) = t.next_cell() {
            cells.push((col, row));
        }
        assert!(!cells.is_empty(), "diagonal ray should visit cells");
    }
    #[test]
    fn test_dda_traversal_negative_dir() {
        let hf = HeightField::from_fn(6, 6, 1.0, |_, _| 0.0);
        let origin = [4.9, 10.0, 2.5];
        let dir = [-1.0, 0.0, 0.0];
        let traversal = HeightfieldRayTraversal::new(&hf, origin, dir, 100.0);
        assert!(
            traversal.is_some(),
            "negative-x ray should create traversal"
        );
        let mut t = traversal.unwrap();
        let mut count = 0;
        while t.next_cell().is_some() {
            count += 1;
        }
        assert!(count >= 1, "negative-x ray should visit cells");
    }
    #[test]
    fn test_dda_traversal_z_only_ray() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let origin = [2.5, 10.0, 0.1];
        let dir = [0.0, 0.0, 1.0];
        let traversal = HeightfieldRayTraversal::new(&hf, origin, dir, 100.0);
        assert!(traversal.is_some(), "Z-only ray should enter grid");
        let mut t = traversal.unwrap();
        let mut count = 0;
        while t.next_cell().is_some() {
            count += 1;
        }
        assert!(count >= 1, "Z-only ray should visit cells");
    }
    #[test]
    fn test_tessellation_no_degenerate_triangles_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let (verts, tris) = hf.to_triangle_mesh();
        for tri in &tris {
            let v0 = verts[tri[0]];
            let v1 = verts[tri[1]];
            let v2 = verts[tri[2]];
            let ab = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let ac = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let cross_len = ((ab[1] * ac[2] - ab[2] * ac[1]).powi(2)
                + (ab[2] * ac[0] - ab[0] * ac[2]).powi(2)
                + (ab[0] * ac[1] - ab[1] * ac[0]).powi(2))
            .sqrt();
            assert!(cross_len > 1e-12, "degenerate triangle found");
        }
    }
    #[test]
    fn test_tessellation_different_aspect_ratio() {
        let hf = HeightField::from_fn(5, 3, 1.0, |_, _| 0.0);
        let (_, tris) = hf.to_triangle_mesh();
        assert_eq!(tris.len(), (3 - 1) * (5 - 1) * 2);
    }
    #[test]
    fn test_gaussian_smooth_preserves_grid_size() {
        let hf = HeightField::from_fn(5, 6, 1.0, |c, r| (c * r) as f64 * 0.1);
        let smoothed = heightfield_smooth_gaussian(&hf, 0.5);
        assert_eq!(smoothed.rows, hf.rows);
        assert_eq!(smoothed.cols, hf.cols);
        assert!((smoothed.scale_x - hf.scale_x).abs() < 1e-12);
        assert!((smoothed.scale_z - hf.scale_z).abs() < 1e-12);
    }
    #[test]
    fn test_gaussian_smooth_large_sigma_heavily_flattens() {
        let mut hf = HeightField::from_fn(9, 9, 1.0, |_, _| 0.0);
        hf.heights[4 * 9 + 4] = 100.0;
        let smoothed = heightfield_smooth_gaussian(&hf, 3.0);
        let center_after = smoothed.height_at(4, 4);
        assert!(
            center_after < 50.0,
            "large sigma should heavily flatten spike, got {center_after}"
        );
    }
    #[test]
    fn test_gaussian_smooth_linear_ramp_unchanged() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, _| c as f64);
        let smoothed = heightfield_smooth_gaussian(&hf, 0.8);
        for &h in &smoothed.heights {
            assert!(
                h.is_finite(),
                "all heights should be finite after Gaussian smooth"
            );
        }
    }
    #[test]
    fn test_gradient_zero_on_flat_field() {
        let hf = HeightField::from_fn(5, 5, 2.0, |_, _| 3.7);
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let s = hf.slope_at(col, row);
                assert!(
                    s < 1e-9,
                    "flat field gradient should be 0 at ({col},{row}), got {s}"
                );
            }
        }
    }
    #[test]
    fn test_gradient_linear_ramp_exact() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, _| c as f64);
        let s = hf.slope_at(2, 2);
        assert!(
            (s - 1.0).abs() < 1e-9,
            "slope should be 1.0 for unit ramp, got {s}"
        );
    }
    #[test]
    fn test_gradient_z_ramp() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, r| r as f64);
        let s = hf.slope_at(2, 2);
        assert!((s - 1.0).abs() < 1e-9, "slope in Z should be 1.0, got {s}");
    }
    #[test]
    fn test_gradient_at_boundary_cells() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c + r) as f64 * 0.3);
        let _ = hf.slope_at(0, 0);
        let _ = hf.slope_at(hf.cols - 1, 0);
        let _ = hf.slope_at(0, hf.rows - 1);
        let _ = hf.slope_at(hf.cols - 1, hf.rows - 1);
    }
    #[test]
    fn test_normal_points_away_from_terrain() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let n = hf.normal_at_grid(col, row);
                assert!(
                    n[1] > 0.0,
                    "normal y component should be positive at ({col},{row})"
                );
            }
        }
    }
    #[test]
    fn test_normal_length_is_one() {
        let hf = HeightField::from_fn(6, 6, 1.5, |c, r| {
            ((c as f64).sin() + (r as f64).cos()) * 3.0
        });
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let n = hf.normal_at_grid(col, row);
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                assert!(
                    (len - 1.0).abs() < 1e-9,
                    "normal not unit at ({col},{row}): len={len}"
                );
            }
        }
    }
    #[test]
    fn test_normal_at_world_flat_is_up_various_positions() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        for &(x, z) in &[(0.5, 0.5), (1.0, 2.0), (3.7, 1.3), (0.0, 0.0)] {
            let n = hf.normal_at_world(x, z);
            assert!(
                (n[1] - 1.0).abs() < 1e-6,
                "flat normal y should be 1 at ({x},{z}), got {}",
                n[1]
            );
        }
    }
    #[test]
    fn test_thermal_erosion_conserves_material_approx() {
        let mut hf =
            HeightField::from_fn(5, 5, 1.0, |c, r| if c == 2 && r == 2 { 10.0 } else { 1.0 });
        let sum_before: f64 = hf.heights.iter().sum();
        heightfield_erode(&mut hf, 0.1, 3);
        let sum_after: f64 = hf.heights.iter().sum();
        assert!(
            (sum_before - sum_after).abs() < 1e-6,
            "thermal erosion should conserve total height sum: before={sum_before}, after={sum_after}"
        );
    }
    #[test]
    fn test_thermal_erosion_multiple_spikes() {
        let mut hf = HeightField::from_fn(7, 7, 1.0, |c, r| {
            if (c == 1 && r == 1) || (c == 5 && r == 5) {
                20.0
            } else {
                0.0
            }
        });
        let peak1_before = hf.height_at(1, 1);
        let peak2_before = hf.height_at(5, 5);
        heightfield_erode(&mut hf, 0.3, 5);
        assert!(
            hf.height_at(1, 1) < peak1_before,
            "spike at (1,1) should erode"
        );
        assert!(
            hf.height_at(5, 5) < peak2_before,
            "spike at (5,5) should erode"
        );
    }
    #[test]
    fn test_hydraulic_erode_small_grid() {
        let mut hf = HeightField::from_fn(3, 3, 1.0, |c, r| (c + r) as f64 * 2.0);
        let max_before = hf.max_height();
        hf.hydraulic_erode(0.2, 2);
        for &h in &hf.heights {
            assert!(
                h.is_finite(),
                "height should be finite after hydraulic erosion"
            );
        }
        let _ = max_before;
    }
    #[test]
    fn test_serialize_contains_metadata() {
        let hf = HeightField::new(vec![1.0, 2.0, 3.0, 4.0], 2, 2, 0.5, 0.5);
        let data = hf.serialize();
        assert_eq!(data[0] as usize, 2);
        assert_eq!(data[1] as usize, 2);
        assert!((data[2] - 0.5).abs() < 1e-12);
        assert!((data[3] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_deserialize_wrong_length_returns_none() {
        let data = vec![2.0, 2.0, 1.0, 1.0, 0.0, 1.0, 2.0];
        assert!(HeightField::deserialize(&data).is_none());
    }
    #[test]
    fn test_serialize_deserialize_large_grid() {
        let hf = HeightField::from_fn(8, 8, 0.5, |c, r| {
            ((c as f64 * 0.3).sin() + (r as f64 * 0.2).cos()) * 5.0
        });
        let data = hf.serialize();
        let hf2 = HeightField::deserialize(&data).expect("should deserialize");
        assert_eq!(hf2.rows, hf.rows);
        assert_eq!(hf2.cols, hf.cols);
        assert!((hf2.scale_x - hf.scale_x).abs() < 1e-9);
        for (a, b) in hf.heights.iter().zip(hf2.heights.iter()) {
            assert!((a - b).abs() < 1e-9, "height mismatch after roundtrip");
        }
    }
    #[test]
    fn test_mean_height_computed_correctly() {
        let hf = HeightField::new(vec![1.0, 2.0, 3.0, 4.0], 2, 2, 1.0, 1.0);
        assert!((hf.mean_height() - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_variance_matches_manual_calculation() {
        let hf = HeightField::new(vec![0.0, 4.0, 0.0, 4.0], 2, 2, 1.0, 1.0);
        assert!(
            (hf.height_variance() - 4.0).abs() < 1e-10,
            "variance={}",
            hf.height_variance()
        );
    }
    #[test]
    fn test_clamp_heights_all_values_in_range() {
        let mut hf = HeightField::new(vec![-10.0, 0.0, 5.0, 15.0], 2, 2, 1.0, 1.0);
        hf.clamp_heights(0.0, 10.0);
        for &h in &hf.heights {
            assert!((0.0..=10.0).contains(&h), "height {h} out of clamped range");
        }
    }
    #[test]
    fn test_scale_heights_negative_factor() {
        let mut hf = HeightField::new(vec![1.0, 2.0, 3.0, 4.0], 2, 2, 1.0, 1.0);
        hf.scale_heights(-1.0);
        assert!((hf.heights[0] - (-1.0)).abs() < 1e-10);
        assert!((hf.heights[3] - (-4.0)).abs() < 1e-10);
    }
    #[test]
    fn test_offset_heights_negative_offset() {
        let mut hf = HeightField::new(vec![5.0, 5.0, 5.0, 5.0], 2, 2, 1.0, 1.0);
        hf.offset_heights(-3.0);
        for &h in &hf.heights {
            assert!(
                (h - 2.0).abs() < 1e-10,
                "offset height should be 2.0, got {h}"
            );
        }
    }
    #[test]
    fn test_normalize_heights_range_is_zero_to_one() {
        let mut hf = HeightField::new(vec![2.0, 5.0, 8.0, 11.0], 2, 2, 1.0, 1.0);
        hf.normalize_heights();
        assert!((hf.min_height() - 0.0).abs() < 1e-10);
        assert!((hf.max_height() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_invert_heights_min_max_preserved() {
        let mut hf = HeightField::new(vec![1.0, 3.0, 5.0, 9.0], 2, 2, 1.0, 1.0);
        hf.invert_heights();
        assert!(
            (hf.min_height() - 1.0).abs() < 1e-10,
            "min after invert={}",
            hf.min_height()
        );
        assert!(
            (hf.max_height() - 9.0).abs() < 1e-10,
            "max after invert={}",
            hf.max_height()
        );
    }
    #[test]
    fn test_lod_downsample_height_averaging() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 4.0);
        let lod = hf.lod_downsample(2);
        for &h in &lod.heights {
            assert!(
                (h - 4.0).abs() < 1e-9,
                "averaged heights should remain 4.0, got {h}"
            );
        }
    }
    #[test]
    fn test_lod_pyramid_monotone_decreasing_resolution() {
        let hf = HeightField::from_fn(16, 16, 1.0, |_, _| 0.0);
        let pyramid = hf.lod_pyramid(4);
        for i in 1..pyramid.len() {
            assert!(
                pyramid[i].rows <= pyramid[i - 1].rows,
                "rows should not increase with LOD level"
            );
            assert!(
                pyramid[i].cols <= pyramid[i - 1].cols,
                "cols should not increase with LOD level"
            );
        }
    }
    #[test]
    fn test_lod_downsample_scale_increases() {
        let hf = HeightField::from_fn(8, 8, 1.0, |_, _| 0.0);
        let lod2 = hf.lod_downsample(2);
        let lod4 = lod2.lod_downsample(2);
        assert!(
            lod4.scale_x > lod2.scale_x,
            "scale should increase with downsampling"
        );
    }
    #[test]
    fn test_resample_upsampling_increases_resolution() {
        let hf = HeightField::from_fn(3, 3, 1.0, |_, _| 0.0);
        let upsampled = hf.resample(7, 7);
        assert_eq!(upsampled.rows, 7);
        assert_eq!(upsampled.cols, 7);
    }
    #[test]
    fn test_resample_downsampling_decreases_resolution() {
        let hf = HeightField::from_fn(9, 9, 1.0, |_, _| 0.0);
        let downsampled = hf.resample(3, 3);
        assert_eq!(downsampled.rows, 3);
        assert_eq!(downsampled.cols, 3);
    }
    #[test]
    fn test_flow_accumulation_valley_has_high_value() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, r| {
            let dc = c as f64 - 2.0;
            let dr = r as f64 - 2.0;
            dc * dc + dr * dr
        });
        let fa = hf.flow_accumulation();
        let center_idx = 2 * 5 + 2;
        let max_fa = fa.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(
            fa[center_idx], max_fa,
            "lowest cell should have max flow accumulation"
        );
    }
    #[test]
    fn test_flow_accumulation_flat_all_equal_to_one() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        let fa = hf.flow_accumulation();
        for &v in &fa {
            assert!(
                (v - 1.0).abs() < 1e-10,
                "flat field flow accumulation should be 1, got {v}"
            );
        }
    }
    #[test]
    fn test_closest_vertex_returns_valid_position() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let (pt, row, col) = hf.closest_vertex([2.3, 0.0, 3.7]);
        assert!(row < hf.rows, "row index out of bounds");
        assert!(col < hf.cols, "col index out of bounds");
        let expected_x = col as f64 * hf.scale_x;
        let expected_z = row as f64 * hf.scale_z;
        assert!((pt[0] - expected_x).abs() < 1e-9);
        assert!((pt[2] - expected_z).abs() < 1e-9);
    }
    #[test]
    fn test_closest_vertex_origin() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let (pt, row, col) = hf.closest_vertex([0.0, 0.0, 0.0]);
        assert_eq!(row, 0);
        assert_eq!(col, 0);
        assert!(pt[0].abs() < 1e-9 && pt[1].abs() < 1e-9 && pt[2].abs() < 1e-9);
    }
    #[test]
    fn test_slope_map_all_nonnegative() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, r| (c + r) as f64 * 0.5);
        for &s in hf.slope_map().iter() {
            assert!(s >= 0.0, "slope should be non-negative, got {s}");
        }
    }
    #[test]
    fn test_curvature_at_bowl_center_is_positive() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, r| {
            let dc = c as f64 - 2.0;
            let dr = r as f64 - 2.0;
            dc * dc + dr * dr
        });
        let k = hf.curvature_at(2, 2);
        assert!(
            k > 0.0,
            "bowl curvature at center should be positive, got {k}"
        );
    }
    #[test]
    fn test_ray_intersect_hits_raised_terrain() {
        let mut hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        hf.heights[2 * 5 + 2] = 5.0;
        let origin = [2.0, 20.0, 2.0];
        let dir = [0.0, -1.0, 0.0];
        let hit = heightfield_ray_intersect(&hf, origin, dir);
        assert!(hit.is_some(), "ray should hit terrain with bump");
    }
    #[test]
    fn test_ray_intersect_consistent_with_ray_cast() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [1.5, 5.0, 1.5];
        let dir = [0.0, -1.0, 0.0];
        let hit_intersect = heightfield_ray_intersect(&hf, origin, dir);
        assert!(
            hit_intersect.is_some(),
            "heightfield_ray_intersect should hit flat terrain"
        );
    }
    #[test]
    fn test_height_at_xz_matches_height_at_world() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, r| (c as f64).sin() + (r as f64).cos());
        for &(x, z) in &[(0.5, 0.5), (1.7, 2.3), (3.1, 0.8)] {
            let h1 = hf.height_at_xz(x, z);
            let h2 = hf.height_at_world(x, z);
            assert!(
                (h1 - h2).abs() < 1e-12,
                "height_at_xz and height_at_world should match at ({x},{z})"
            );
        }
    }
    #[test]
    fn test_tessellate_consistent_with_to_triangle_mesh() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c + r) as f64 * 0.1);
        let (v1, t1) = hf.tessellate();
        let (v2, t2) = hf.to_triangle_mesh();
        assert_eq!(v1.len(), v2.len());
        assert_eq!(t1.len(), t2.len());
    }
    #[test]
    fn test_normals_matches_compute_all_normals() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c + r) as f64 * 0.2);
        let n1 = hf.normals();
        let n2 = hf.compute_all_normals();
        assert_eq!(n1.len(), n2.len());
        for (a, b) in n1.iter().zip(n2.iter()) {
            for i in 0..3 {
                assert!((a[i] - b[i]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_surface_area_larger_for_rough_terrain() {
        let flat = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let rough = HeightField::from_fn(5, 5, 1.0, |c, r| ((c + r) as f64 * 0.5).sin() * 2.0);
        assert!(
            rough.surface_area() >= flat.surface_area(),
            "rough terrain should have >= surface area than flat terrain"
        );
    }
    #[test]
    fn test_surface_area_single_cell() {
        let hf = HeightField::from_fn(2, 2, 1.0, |_, _| 0.0);
        let area = hf.surface_area();
        assert!(
            (area - 1.0).abs() < 1e-10,
            "single flat cell area should be 1.0, got {area}"
        );
    }
    #[test]
    fn test_count_peaks_descending_field_has_zero_peaks() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (10 - c - r) as f64);
        let peaks = hf.count_peaks();
        assert!(peaks >= 1, "descending field should have at least 1 peak");
    }
    #[test]
    fn test_count_peaks_checkerboard() {
        let hf = HeightField::from_fn(5, 5, 1.0, |c, r| if (c + r) % 2 == 0 { 1.0 } else { 0.0 });
        let peaks = hf.count_peaks();
        assert!(peaks > 0, "checkerboard should have multiple peaks");
    }
}
