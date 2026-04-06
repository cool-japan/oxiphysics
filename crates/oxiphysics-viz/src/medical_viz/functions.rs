//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types_core::*;
use super::types_ext::*;
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_rgba_blend_transparent_background() {
        let c = Rgba::new(1.0, 0.0, 0.0, 1.0);
        let bg = Rgba::transparent();
        let result = c.blend_over(bg);
        assert!((result.r - 1.0).abs() < 1e-6);
        assert!((result.a - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_transfer_func_endpoints() {
        let tf = ColorTransferFunc::bone_preset();
        let c0 = tf.evaluate(0.0);
        let c1 = tf.evaluate(1.0);
        assert!(c0.a < 0.1, "c0 alpha={}", c0.a);
        assert!(c1.a > 0.9, "c1 alpha={}", c1.a);
    }
    #[test]
    fn test_transfer_func_monotone_alpha() {
        let tf = ColorTransferFunc::bone_preset();
        let a0 = tf.evaluate(0.0).a;
        let a1 = tf.evaluate(0.5).a;
        let a2 = tf.evaluate(1.0).a;
        assert!(a1 >= a0 - 1e-6);
        assert!(a2 >= a1 - 1e-6);
    }
    #[test]
    fn test_volume_data_get_set() {
        let mut v = VolumeData::new([4, 4, 4], [1.0; 3]);
        v.set(1, 2, 3, 0.75);
        assert!((v.get(1, 2, 3) - 0.75).abs() < 1e-6);
    }
    #[test]
    fn test_volume_data_sample_center() {
        let mut v = VolumeData::new([4, 4, 4], [1.0; 3]);
        for i in 0..v.data.len() {
            v.data[i] = 0.5;
        }
        let s = v.sample(0.5, 0.5, 0.5);
        assert!((s - 0.5).abs() < 1e-6, "s={s}");
    }
    #[test]
    fn test_volume_data_value_range() {
        let mut v = VolumeData::new([2, 2, 2], [1.0; 3]);
        v.data[0] = -1.0;
        v.data[7] = 1.0;
        let (mn, mx) = v.value_range();
        assert!((mn + 1.0).abs() < 1e-6);
        assert!((mx - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_volume_renderer_ray_casting() {
        let mut vol = VolumeData::new([4, 4, 4], [1.0; 3]);
        for v in &mut vol.data {
            *v = 1.0;
        }
        let tf = ColorTransferFunc::bone_preset();
        let renderer = VolumeRenderer::new(50);
        let color = renderer.cast_ray(&vol, &tf, [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]);
        assert!(color.a > 0.5, "alpha={}", color.a);
    }
    #[test]
    fn test_mip_renderer_max_value() {
        let mut vol = VolumeData::new([2, 2, 4], [1.0; 3]);
        vol.set(0, 0, 3, 0.9);
        let mip = MaximumIntensityProj::new(50);
        let v = mip.cast_ray(&vol, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert!(v > 0.5, "v={v}");
    }
    #[test]
    fn test_isosurface_all_inside() {
        let renderer = IsosurfaceRenderer::new(0.5);
        let corners = [1.0f32; 8];
        let positions = [[0.0f32; 3]; 8];
        let tris = renderer.marching_cube_triangles(corners, positions);
        assert!(tris.is_empty());
    }
    #[test]
    fn test_isosurface_mixed_produces_tris() {
        let renderer = IsosurfaceRenderer::new(0.5);
        let mut corners = [0.0f32; 8];
        corners[0] = 1.0;
        let positions: [[f32; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let tris = renderer.marching_cube_triangles(corners, positions);
        assert!(!tris.is_empty());
    }
    #[test]
    fn test_slice_plane_axial_extraction() {
        let mut vol = VolumeData::new([4, 4, 4], [1.0; 3]);
        for v in &mut vol.data {
            *v = 0.7;
        }
        let plane = SlicePlane::axial(0.5, 4, 4);
        let pixels = plane.extract(&vol);
        assert_eq!(pixels.len(), 16);
        for &p in &pixels {
            assert!((p - 0.7).abs() < 0.01, "p={p}");
        }
    }
    #[test]
    fn test_anatomy_atlas_find_by_label() {
        let mut atlas = AnatomyAtlas::new();
        atlas.add_region(AtlasRegion::new("bone", "Bone", Rgba::white(), 1));
        let r = atlas.find_by_label(1);
        assert!(r.is_some());
        assert_eq!(r.unwrap().name, "bone");
    }
    #[test]
    fn test_anatomy_atlas_toggle_visibility() {
        let mut atlas = AnatomyAtlas::new();
        atlas.add_region(AtlasRegion::new(
            "skin",
            "Skin",
            Rgba::new(0.9, 0.7, 0.6, 1.0),
            2,
        ));
        atlas.toggle_visibility("skin");
        let c = atlas.color_for_label(2);
        assert!(c.a < 1e-6);
    }
    #[test]
    fn test_measurement_distance() {
        let mut mt = MeasurementTool::new();
        let d = mt.measure_distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-10, "d={d}");
    }
    #[test]
    fn test_measurement_angle_90() {
        let mut mt = MeasurementTool::new();
        let a = mt.measure_angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((a - 90.0).abs() < 1e-6, "a={a}");
    }
    #[test]
    fn test_measurement_voi() {
        let mut mt = MeasurementTool::new();
        let v = mt.measure_voi([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        assert!((v - 24.0).abs() < 1e-10, "v={v}");
    }
    #[test]
    fn test_heatmap_color_at_center() {
        let data = vec![0.0f32, 0.5, 1.0];
        let hm = HeatmapOverlay::new(data, 0.8);
        let c = hm.color_at(1);
        assert!(c.a > 0.0);
    }
    #[test]
    fn test_heatmap_out_of_bounds() {
        let hm = HeatmapOverlay::new(vec![0.5], 1.0);
        let c = hm.color_at(999);
        assert!(c.a < 1e-6);
    }
    #[test]
    fn test_vessel_tree_segment_count() {
        let mut tree = VesselTree::new();
        tree.add_segment(VesselSegment::new([0.0; 3], [0.0, 1.0, 0.0], 0.05, 0.04));
        assert_eq!(tree.num_segments(), 1);
    }
    #[test]
    fn test_vessel_tree_dfs() {
        let mut tree = VesselTree::new();
        let root = tree.add_segment(VesselSegment::new([0.0; 3], [1.0, 0.0, 0.0], 0.1, 0.08));
        tree.add_child(
            root,
            VesselSegment::new([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.07, 0.06),
        );
        let order = tree.dfs_order();
        assert_eq!(order[0], 0);
        assert_eq!(order.len(), 2);
    }
    #[test]
    fn test_vessel_segment_length() {
        let seg = VesselSegment::new([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 0.05, 0.05);
        assert!((seg.length() - 5.0).abs() < 1e-5);
    }
    #[test]
    fn test_electro_physiology_activated_vertex() {
        let ep = ElectroPhysiology::new(vec![0.0, 50.0], 100.0);
        let c = ep.vertex_color(0);
        assert!(c.r > c.b, "should be reddish");
    }
    #[test]
    fn test_electro_physiology_time_advance() {
        let mut ep = ElectroPhysiology::new(vec![10.0], 100.0);
        ep.advance(20.0);
        assert!((ep.current_time - 20.0).abs() < 1e-6);
    }
    #[test]
    fn test_electro_physiology_isochrones() {
        let ep = ElectroPhysiology::new(vec![10.0, 10.5, 50.0], 100.0);
        let pairs = ep.isochrone_pairs(10.0, 1.0);
        assert!(!pairs.is_empty());
    }
    #[test]
    fn test_rgba_premultiply() {
        let c = Rgba::new(1.0, 1.0, 1.0, 0.5);
        let pm = c.pre_multiply_alpha();
        assert!((pm.r - 0.5).abs() < 1e-6);
        assert!((pm.a - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_soft_tissue_transfer_func() {
        let tf = ColorTransferFunc::soft_tissue_preset();
        let c = tf.evaluate(0.5);
        assert!(c.r > 0.0, "should be reddish at 0.5");
    }
    #[test]
    fn test_volume_gradient_near_edge() {
        let mut vol = VolumeData::new([8, 8, 8], [1.0; 3]);
        let [nx, ny, nz] = vol.dims;
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    vol.set(ix, iy, iz, ix as f32 / nx as f32);
                }
            }
        }
        let grad = vol.gradient(0.5, 0.5, 0.5, 0.05);
        assert!(grad[0].abs() > 1e-5, "gx={}", grad[0]);
    }
    #[test]
    fn test_heatmap_update_data() {
        let mut hm = HeatmapOverlay::new(vec![0.0, 1.0], 1.0);
        hm.update_data(vec![2.0, 5.0]);
        assert!((hm.min_val - 2.0).abs() < 1e-6);
        assert!((hm.max_val - 5.0).abs() < 1e-6);
    }
}
#[cfg(test)]
mod expanded_medical_tests {
    use super::*;
    #[test]
    fn test_incision_length() {
        let mut inc = IncisionPath::new();
        inc.add_point([0.0, 0.0, 0.0], 0.005, 0);
        inc.add_point([0.1, 0.0, 0.0], 0.005, 0);
        inc.add_point([0.2, 0.0, 0.0], 0.005, 0);
        assert!(
            (inc.total_length() - 0.2).abs() < 1e-10,
            "len={}",
            inc.total_length()
        );
    }
    #[test]
    fn test_incision_avg_depth() {
        let mut inc = IncisionPath::new();
        inc.add_point([0.0, 0.0, 0.0], 0.01, 0);
        inc.add_point([0.1, 0.0, 0.0], 0.03, 0);
        assert!((inc.average_depth() - 0.02).abs() < 1e-10);
    }
    #[test]
    fn test_trocar_volume() {
        let mut port = TrocarPort::new([0.0, 0.0, 0.0], 10.0, [0.0, 1.0, 0.0], "camera");
        port.depth_to_target = 150.0;
        let vol = port.working_cone_volume();
        assert!(vol > 0.0, "vol={vol}");
    }
    #[test]
    fn test_surgical_plan_no_collision() {
        let mut plan = SurgicalPlan::new(5.0);
        plan.add_port(TrocarPort::new([0.0, 0.0, 0.0], 10.0, [0.0, 1.0, 0.0], "A"));
        plan.add_port(TrocarPort::new(
            [50.0, 0.0, 0.0],
            10.0,
            [0.0, 1.0, 0.0],
            "B",
        ));
        assert!(
            plan.ports_collision_free(),
            "Well-separated ports should not collide"
        );
    }
    #[test]
    fn test_surgical_plan_collision() {
        let mut plan = SurgicalPlan::new(5.0);
        plan.add_port(TrocarPort::new([0.0, 0.0, 0.0], 12.0, [0.0, 1.0, 0.0], "A"));
        plan.add_port(TrocarPort::new([5.0, 0.0, 0.0], 12.0, [0.0, 1.0, 0.0], "B"));
        assert!(!plan.ports_collision_free(), "Close ports should collide");
    }
    #[test]
    fn test_vessel_length() {
        let seg = VesselSegmentExt::new(
            1,
            0,
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            2.0,
            1.5,
            VesselTypeExt::Artery,
        );
        assert!((seg.length() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_vessel_volume() {
        let seg = VesselSegmentExt::new(
            1,
            0,
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            2.0,
            2.0,
            VesselTypeExt::Vein,
        );
        let vol = seg.volume();
        let expected = std::f64::consts::PI * 4.0 * 10.0;
        assert!(
            (vol - expected).abs() / expected < 1e-3,
            "vol={vol}, expected={expected}"
        );
    }
    #[test]
    fn test_vessel_tree_length() {
        let mut tree = VesselTreeExt::new("liver");
        tree.add_segment(VesselSegmentExt::new(
            1,
            0,
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            1.0,
            1.0,
            VesselTypeExt::Artery,
        ));
        tree.add_segment(VesselSegmentExt::new(
            2,
            1,
            [5.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.8,
            0.8,
            VesselTypeExt::Artery,
        ));
        assert!((tree.total_length() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_tumor_volume() {
        let tumor = TumorRegion::spherical([0.0, 0.0, 0.0], 10.0);
        let vol = tumor.volume();
        let expected = (4.0 / 3.0) * std::f64::consts::PI * 1000.0;
        assert!((vol - expected).abs() / expected < 1e-6, "vol={vol}");
    }
    #[test]
    fn test_tumor_contains_inside() {
        let tumor = TumorRegion::spherical([0.0, 0.0, 0.0], 10.0);
        assert!(tumor.contains([1.0, 1.0, 1.0]));
        assert!(!tumor.contains([15.0, 0.0, 0.0]));
    }
    #[test]
    fn test_dose_set_get() {
        let mut dose = DoseDistribution::new([4, 4, 4], [1.0; 3], [0.0; 3]);
        dose.set(1, 2, 3, 50.0);
        assert!((dose.get(1, 2, 3) - 50.0).abs() < 0.01);
    }
    #[test]
    fn test_dose_max() {
        let mut dose = DoseDistribution::new([2, 2, 2], [1.0; 3], [0.0; 3]);
        dose.set(0, 0, 0, 60.0);
        dose.set(1, 1, 1, 45.0);
        assert!((dose.max_dose() - 60.0).abs() < 0.01);
    }
    #[test]
    fn test_dose_dvh() {
        let mut dose = DoseDistribution::new([4, 1, 1], [1.0; 3], [0.0; 3]);
        for i in 0..4 {
            dose.set(i, 0, 0, i as f32 * 20.0);
        }
        let dvh = dose.dvh(4);
        assert_eq!(dvh.len(), 4);
        assert!((dvh[0].1 - 1.0).abs() < 0.01, "dvh[0].1={}", dvh[0].1);
    }
    #[test]
    fn test_dose_v20() {
        let mut dose = DoseDistribution::new([4, 1, 1], [1.0; 3], [0.0; 3]);
        dose.set(0, 0, 0, 60.0);
        dose.set(1, 0, 0, 40.0);
        dose.set(2, 0, 0, 20.0);
        dose.set(3, 0, 0, 0.0);
        let v20 = dose.volume_at_dose(20.0);
        assert!((v20 - 0.75).abs() < 0.01, "v20={v20}");
    }
    #[test]
    fn test_dicom_to_hu() {
        let meta = DicomMetadata::default_ct();
        let hu = meta.to_hu(0);
        assert!((hu - (-1024.0)).abs() < 1e-6, "hu={hu}");
        let hu_water = meta.to_hu(1024);
        assert!(hu_water.abs() < 1e-6, "hu_water={hu_water}");
    }
    #[test]
    fn test_dicom_voxel_volume() {
        let meta = DicomMetadata::default_ct();
        assert!(meta.voxel_volume() > 0.0);
    }
    #[test]
    fn test_vessel_arteries_filter() {
        let mut tree = VesselTreeExt::new("kidney");
        tree.add_segment(VesselSegmentExt::new(
            1,
            0,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            1.0,
            VesselTypeExt::Artery,
        ));
        tree.add_segment(VesselSegmentExt::new(
            2,
            0,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            0.5,
            0.5,
            VesselTypeExt::Vein,
        ));
        assert_eq!(tree.arteries().len(), 1);
    }
    #[test]
    fn test_tumor_isosurface() {
        let tumor = TumorRegion::spherical([0.0, 0.0, 0.0], 5.0);
        let verts = tumor.isosurface_vertices(8, 16);
        assert!(!verts.is_empty(), "isosurface vertices should be generated");
    }
    #[test]
    fn test_dose_sample_at() {
        let mut dose = DoseDistribution::new([4, 4, 4], [1.0; 3], [0.0; 3]);
        dose.set(1, 1, 1, 55.0);
        let d = dose.sample_at([1.0, 1.0, 1.0]);
        assert!(d > 0.0, "Sample should be > 0 near voxel, d={d}");
    }
    #[test]
    fn test_surgical_incision_total() {
        let mut plan = SurgicalPlan::new(5.0);
        let mut inc = IncisionPath::new();
        inc.add_point([0.0, 0.0, 0.0], 0.005, 0);
        inc.add_point([0.05, 0.0, 0.0], 0.005, 0);
        plan.incisions.push(inc);
        assert!((plan.total_incision_length() - 0.05).abs() < 1e-10);
    }
}
/// Extract isodose contours from a 2D dose slice using marching squares (simplified).
///
/// Returns one `IsodoseContour` per dose level.
pub fn extract_isodose_contours(
    dose_slice: &[f32],
    width: usize,
    height: usize,
    dose_levels: &[f32],
    pixel_spacing: [f32; 2],
) -> Vec<IsodoseContour> {
    let mut contours: Vec<IsodoseContour> = dose_levels
        .iter()
        .map(|&dl| {
            let hue = dl
                / dose_levels
                    .iter()
                    .cloned()
                    .fold(f32::NEG_INFINITY, f32::max)
                    .max(1.0);
            let color = Rgba::new(hue, 1.0 - hue, 0.5, 1.0);
            IsodoseContour::new(dl, color)
        })
        .collect();
    if width == 0 || height == 0 {
        return contours;
    }
    let get = |x: usize, y: usize| -> f32 {
        if x >= width || y >= height {
            0.0
        } else {
            dose_slice[y * width + x]
        }
    };
    for (ci, contour) in contours.iter_mut().enumerate() {
        let level = dose_levels[ci];
        for cy in 0..(height.saturating_sub(1)) {
            for cx in 0..(width.saturating_sub(1)) {
                let c00 = get(cx, cy);
                let c10 = get(cx + 1, cy);
                let c01 = get(cx, cy + 1);
                let c11 = get(cx + 1, cy + 1);
                let interp = |a: f32, b: f32| -> f32 {
                    if (b - a).abs() < 1e-10 {
                        0.5
                    } else {
                        (level - a) / (b - a)
                    }
                };
                let mut cross_pts: Vec<[f32; 2]> = Vec::new();
                if (c00 >= level) != (c10 >= level) {
                    let t = interp(c00, c10);
                    cross_pts.push([
                        (cx as f32 + t) * pixel_spacing[0],
                        cy as f32 * pixel_spacing[1],
                    ]);
                }
                if (c10 >= level) != (c11 >= level) {
                    let t = interp(c10, c11);
                    cross_pts.push([
                        (cx as f32 + 1.0) * pixel_spacing[0],
                        (cy as f32 + t) * pixel_spacing[1],
                    ]);
                }
                if (c01 >= level) != (c11 >= level) {
                    let t = interp(c01, c11);
                    cross_pts.push([
                        (cx as f32 + t) * pixel_spacing[0],
                        (cy as f32 + 1.0) * pixel_spacing[1],
                    ]);
                }
                if (c00 >= level) != (c01 >= level) {
                    let t = interp(c00, c01);
                    cross_pts.push([
                        cx as f32 * pixel_spacing[0],
                        (cy as f32 + t) * pixel_spacing[1],
                    ]);
                }
                let mut i = 0;
                while i + 1 < cross_pts.len() {
                    contour.add_segment(cross_pts[i], cross_pts[i + 1]);
                    i += 2;
                }
            }
        }
    }
    contours
}
#[cfg(test)]
mod expanded_viz_tests {
    use super::*;
    #[test]
    fn test_ct_window_soft_tissue_clamp() {
        let w = CtWindow::soft_tissue();
        assert!((w.apply(-2000.0)).abs() < 1e-6, "below range should be 0");
        assert!(
            (w.apply(2000.0) - 1.0).abs() < 1e-6,
            "above range should be 1"
        );
    }
    #[test]
    fn test_ct_pixel_to_hu() {
        let hu = CtWindow::pixel_to_hu(0, 1.0, -1024.0);
        assert!((hu - (-1024.0)).abs() < 1e-4);
        let hu2 = CtWindow::pixel_to_hu(1024, 1.0, -1024.0);
        assert!(hu2.abs() < 1e-4);
    }
    #[test]
    fn test_ct_slice_intensity_at() {
        let hu = vec![0.0f32; 4];
        let mut hu2 = hu.clone();
        hu2[3] = 100.0;
        let slice = CtSlice::from_hu(hu2, 2, 2, 0, CtWindow::soft_tissue());
        assert!(slice.intensity_at(1, 1) > slice.intensity_at(0, 0) - 1e-6);
    }
    #[test]
    fn test_ct_slice_to_rgba_len() {
        let hu = vec![40.0f32; 9];
        let slice = CtSlice::from_hu(hu, 3, 3, 0, CtWindow::soft_tissue());
        let rgba = slice.to_rgba();
        assert_eq!(rgba.len(), 9);
    }
    #[test]
    fn test_mri_t1_csf_dark() {
        let csf = MriTissue::csf();
        let wm = MriTissue::white_matter();
        let s_csf = csf.t1_weighted_signal(500.0, 15.0);
        let s_wm = wm.t1_weighted_signal(500.0, 15.0);
        assert!(
            s_csf < s_wm,
            "CSF should be darker than WM on T1, csf={s_csf}, wm={s_wm}"
        );
    }
    #[test]
    fn test_mri_t2_csf_bright() {
        let csf = MriTissue::csf();
        let wm = MriTissue::white_matter();
        let s_csf = csf.t2_weighted_signal(100.0);
        let s_wm = wm.t2_weighted_signal(100.0);
        assert!(
            s_csf > s_wm,
            "CSF should be brighter than WM on T2, csf={s_csf}, wm={s_wm}"
        );
    }
    #[test]
    fn test_mri_signal_map_normalise() {
        let tissues = vec![
            MriTissue::white_matter(),
            MriTissue::grey_matter(),
            MriTissue::csf(),
        ];
        let labels = vec![0, 1, 2, 0, 1];
        let map = MriSignalMap::t1_weighted(&labels, &tissues, 5, 1, 500.0, 15.0);
        let norm = map.normalise();
        for &v in &norm {
            assert!((0.0..=1.0).contains(&v), "out of range: {v}");
        }
    }
    #[test]
    fn test_dicom_banner() {
        let info = DicomDisplayInfo::default_ct();
        let banner = info.banner();
        assert!(banner.contains("CT"), "banner={banner}");
    }
    #[test]
    fn test_dicom_is_axial() {
        let info = DicomDisplayInfo::default_ct();
        assert!(info.is_axial(), "default CT should be axial");
    }
    #[test]
    fn test_landmark_project_axial_inside() {
        let lm = AnatomicalLandmark::fiducial([10.0, 20.0, 50.0], "F1");
        let proj = lm.project_axial(50.0, 1.0);
        assert!(proj.is_some());
        let (x, y) = proj.unwrap();
        assert!((x - 10.0).abs() < 1e-10);
        assert!((y - 20.0).abs() < 1e-10);
    }
    #[test]
    fn test_landmark_project_axial_outside() {
        let lm = AnatomicalLandmark::fiducial([0.0, 0.0, 50.0], "F2");
        let proj = lm.project_axial(60.0, 1.0);
        assert!(proj.is_none());
    }
    #[test]
    fn test_landmark_set_within_radius() {
        let mut ls = LandmarkSet::new();
        ls.add(AnatomicalLandmark::fiducial([0.0, 0.0, 0.0], "A"));
        ls.add(AnatomicalLandmark::fiducial([1.0, 0.0, 0.0], "B"));
        ls.add(AnatomicalLandmark::fiducial([100.0, 0.0, 0.0], "C"));
        let near = ls.within_radius([0.0; 3], 5.0);
        assert_eq!(near.len(), 2);
    }
    #[test]
    fn test_landmark_count_by_type() {
        let mut ls = LandmarkSet::new();
        ls.add(AnatomicalLandmark::fiducial([0.0; 3], "F1"));
        ls.add(AnatomicalLandmark::anatomical([1.0, 0.0, 0.0], "A1"));
        ls.add(AnatomicalLandmark::fiducial([2.0, 0.0, 0.0], "F2"));
        assert_eq!(ls.count_by_type(LandmarkType::Fiducial), 2);
        assert_eq!(ls.count_by_type(LandmarkType::Anatomical), 1);
    }
    #[test]
    fn test_isodose_contour_length() {
        let mut contour = IsodoseContour::new(50.0, Rgba::white());
        contour.add_segment([0.0, 0.0], [3.0, 4.0]);
        assert!(
            (contour.total_length() - 5.0).abs() < 1e-5,
            "len={}",
            contour.total_length()
        );
    }
    #[test]
    fn test_extract_isodose_contours() {
        let dose: Vec<f32> = (0..16).map(|i| i as f32 * 4.0).collect();
        let contours = extract_isodose_contours(&dose, 4, 4, &[30.0, 45.0], [1.0, 1.0]);
        assert_eq!(contours.len(), 2);
    }
    #[test]
    fn test_mpr_standard_views_count() {
        let vol = VolumeData::new([8, 8, 8], [1.0; 3]);
        let views = MprView::standard_views(&vol, 0.5, 0.5, 0.5, 8);
        assert_eq!(views.len(), 3);
    }
    #[test]
    fn test_mpr_pixel_at() {
        let mut vol = VolumeData::new([8, 8, 8], [1.0; 3]);
        for v in &mut vol.data {
            *v = 0.5;
        }
        let views = MprView::standard_views(&vol, 0.5, 0.5, 0.5, 8);
        let p = views[0].pixel_at(4, 4);
        assert!((0.0..=1.0).contains(&p), "p={p}");
    }
    #[test]
    fn test_tumor_contour_area() {
        let mut c = TumorContour::new("PTV", Rgba::white(), 0);
        c.add_vertex([0.0, 0.0]);
        c.add_vertex([4.0, 0.0]);
        c.add_vertex([4.0, 3.0]);
        c.add_vertex([0.0, 3.0]);
        let area = c.area();
        assert!((area - 12.0).abs() < 1e-5, "area={area}");
    }
    #[test]
    fn test_tumor_contour_contains() {
        let mut c = TumorContour::new("GTV", Rgba::white(), 0);
        {
            static VERTS: [[f32; 2]; 8] = [
                [0.0, 5.0],
                [3.5, 3.5],
                [5.0, 0.0],
                [3.5, -3.5],
                [0.0, -5.0],
                [-3.5, -3.5],
                [-5.0, 0.0],
                [-3.5, 3.5],
            ];
            for v in VERTS {
                c.add_vertex(v);
            }
        }
        assert!(c.contains([0.0, 0.0]));
        assert!(!c.contains([10.0, 10.0]));
    }
    #[test]
    fn test_planar_measure_distance() {
        let mut tool = PlanarMeasurementTool::new();
        let d = tool.measure_distance_2d([0.0, 0.0], [3.0, 4.0], [1.0, 1.0]);
        assert!((d - 5.0).abs() < 1e-8, "d={d}");
    }
    #[test]
    fn test_pseudo_color_jet_low() {
        let c = PseudoColorMap::Jet.map(0.0);
        assert!(c.b > c.r, "Jet at 0 should be blue, b={}, r={}", c.b, c.r);
    }
    #[test]
    fn test_pseudo_color_hot_high() {
        let c = PseudoColorMap::Hot.map(1.0);
        assert!(
            c.r > 0.9 && c.g > 0.9 && c.b > 0.9,
            "r={} g={} b={}",
            c.r,
            c.g,
            c.b
        );
    }
    #[test]
    fn test_suv_bw_basic() {
        let params = SuvParameters::new(400e6, 70.0);
        let suv = params.suv_bw(10_000.0);
        assert!((suv - 1.75).abs() < 0.01, "suv={suv}");
    }
    #[test]
    fn test_dvh_curve_v_at_dose() {
        let mut dvh = DvhCurve::new("PTV", true, Rgba::white());
        let voxels: Vec<f64> = (0..100).map(|i| i as f64 * 0.6).collect();
        dvh.compute_from_voxels(&voxels, 50);
        let v0 = dvh.v_at_dose(0.0);
        assert!(v0 >= 0.99, "v0={v0}");
    }
    #[test]
    fn test_pet_suv_image_render_len() {
        let params = SuvParameters::new(370e6, 75.0);
        let activity: Vec<f64> = (0..16).map(|i| i as f64 * 1000.0).collect();
        let pet = PetSuvImage::from_activity(&activity, 4, 4, params, PseudoColorMap::Hot);
        let rendered = pet.render();
        assert_eq!(rendered.len(), 16);
    }
}
