//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
/// Apply window/level and rescale slope/intercept to a raw pixel buffer.
///
/// Returns an 8-bit display image (0–255).
pub fn apply_windowing(
    pixels: &[u16],
    rescale_slope: f64,
    rescale_intercept: f64,
    window_center: f64,
    window_width: f64,
) -> Vec<u8> {
    let low = window_center - window_width / 2.0;
    let high = window_center + window_width / 2.0;
    pixels
        .iter()
        .map(|&p| {
            let hu = p as f64 * rescale_slope + rescale_intercept;
            let norm = (hu - low) / (high - low);
            (norm.clamp(0.0, 1.0) * 255.0) as u8
        })
        .collect()
}
/// Trilinear-interpolation resampling of a 3D voxel grid.
///
/// `src` is a flat `f64` buffer of size `nx × ny × nz` (Z-major).
/// Returns a new buffer at the specified output dimensions.
pub fn resample_trilinear(src: &[f64], src_dims: [usize; 3], dst_dims: [usize; 3]) -> Vec<f64> {
    let [snx, sny, snz] = src_dims;
    let [dnx, dny, dnz] = dst_dims;
    let mut out = vec![0.0f64; dnx * dny * dnz];
    let sx = (snx as f64 - 1.0) / (dnx as f64 - 1.0).max(1.0);
    let sy = (sny as f64 - 1.0) / (dny as f64 - 1.0).max(1.0);
    let sz = (snz as f64 - 1.0) / (dnz as f64 - 1.0).max(1.0);
    let clamp_i = |v: isize| v.clamp(0, snx as isize - 1) as usize;
    let clamp_j = |v: isize| v.clamp(0, sny as isize - 1) as usize;
    let clamp_k = |v: isize| v.clamp(0, snz as isize - 1) as usize;
    for k in 0..dnz {
        for j in 0..dny {
            for i in 0..dnx {
                let fx = i as f64 * sx;
                let fy = j as f64 * sy;
                let fz = k as f64 * sz;
                let ix = fx.floor() as isize;
                let iy = fy.floor() as isize;
                let iz = fz.floor() as isize;
                let tx = fx - ix as f64;
                let ty = fy - iy as f64;
                let tz = fz - iz as f64;
                let lerp = |a: f64, b: f64, t: f64| a + t * (b - a);
                macro_rules! sv {
                    ($ii:expr, $jj:expr, $kk:expr) => {
                        src[clamp_k($kk) * snx * sny + clamp_j($jj) * snx + clamp_i($ii)]
                    };
                }
                let c000 = sv!(ix, iy, iz);
                let c100 = sv!(ix + 1, iy, iz);
                let c010 = sv!(ix, iy + 1, iz);
                let c110 = sv!(ix + 1, iy + 1, iz);
                let c001 = sv!(ix, iy, iz + 1);
                let c101 = sv!(ix + 1, iy, iz + 1);
                let c011 = sv!(ix, iy + 1, iz + 1);
                let c111 = sv!(ix + 1, iy + 1, iz + 1);
                let c_x0 = lerp(c000, c100, tx);
                let c_x1 = lerp(c010, c110, tx);
                let c_x2 = lerp(c001, c101, tx);
                let c_x3 = lerp(c011, c111, tx);
                let c_y0 = lerp(c_x0, c_x1, ty);
                let c_y1 = lerp(c_x2, c_x3, ty);
                out[k * dnx * dny + j * dnx + i] = lerp(c_y0, c_y1, tz);
            }
        }
    }
    out
}
/// Extract a 3D point cloud from a CT volume by intensity threshold.
///
/// Returns world-space coordinates for all voxels whose Hounsfield Unit value
/// is ≥ `threshold_hu`.
pub fn ct_point_cloud_threshold(
    pixels: &[u16],
    dims: [usize; 3],
    spacing: [f64; 3],
    origin: [f64; 3],
    rescale_slope: f64,
    rescale_intercept: f64,
    threshold_hu: f64,
) -> Vec<[f64; 3]> {
    let [nx, ny, nz] = dims;
    let mut points = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let idx = k * nx * ny + j * nx + i;
                let hu = pixels[idx] as f64 * rescale_slope + rescale_intercept;
                if hu >= threshold_hu {
                    points.push([
                        origin[0] + i as f64 * spacing[0],
                        origin[1] + j as f64 * spacing[1],
                        origin[2] + k as f64 * spacing[2],
                    ]);
                }
            }
        }
    }
    points
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::medical_io::types::*;
    #[test]
    fn test_dicom_tag_new() {
        let t = DicomTag::new(0x0028, 0x0010);
        assert_eq!(t.group, 0x0028);
        assert_eq!(t.element, 0x0010);
    }
    #[test]
    fn test_dicom_reader_set_get_tag() {
        let mut reader = DicomReader::new("test.dcm");
        reader.set_tag(DicomTag::patient_name(), "Smith^John");
        assert_eq!(reader.patient_name(), Some("Smith^John"));
    }
    #[test]
    fn test_dicom_reader_modality() {
        let mut reader = DicomReader::new("test.dcm");
        reader.set_tag(DicomTag::modality(), "CT");
        assert_eq!(reader.modality(), Some("CT"));
    }
    #[test]
    fn test_dicom_reader_load_synthetic() {
        let mut reader = DicomReader::new("test.dcm");
        reader.load_synthetic(8, 8);
        assert_eq!(reader.rows, 8);
        assert_eq!(reader.columns, 8);
        assert_eq!(reader.pixel_data.len(), 64);
    }
    #[test]
    fn test_dicom_reader_pixel() {
        let mut reader = DicomReader::new("test.dcm");
        reader.load_synthetic(4, 4);
        let p = reader.pixel(0, 0);
        assert!(p.is_some());
        let p_out = reader.pixel(10, 10);
        assert!(p_out.is_none());
    }
    #[test]
    fn test_dicom_image_data_new() {
        let img = DicomImageData::new(512, 512, 100);
        assert_eq!(img.rows, 512);
        assert_eq!(img.n_slices, 100);
    }
    #[test]
    fn test_dicom_image_data_to_hu() {
        let img = DicomImageData::new(10, 10, 1);
        let hu = img.to_hu(1024);
        assert!((hu - 0.0).abs() < 1.0);
    }
    #[test]
    fn test_dicom_image_data_window_level() {
        let img = DicomImageData::new(10, 10, 1);
        let wl = img.window_level(img.window_center);
        assert!(wl >= 127 && wl <= 128);
    }
    #[test]
    fn test_dicom_image_data_physical_size() {
        let img = DicomImageData::new(10, 20, 5);
        let size = img.physical_size();
        assert_eq!(size[0], 10.0);
        assert_eq!(size[1], 20.0);
        assert_eq!(size[2], 5.0);
    }
    #[test]
    fn test_dicom_image_data_voxel() {
        let mut img = DicomImageData::new(4, 4, 2);
        img.pixels[0] = 1000;
        assert_eq!(img.voxel(0, 0, 0), Some(1000));
        assert!(img.voxel(10, 10, 10).is_none());
    }
    #[test]
    fn test_nifti_reader_new() {
        let r = NiftiReader::new("test.nii");
        assert_eq!(r.dim[0], 3);
        assert_eq!(r.datatype, NiftiDtype::Float32);
    }
    #[test]
    fn test_nifti_reader_n_voxels() {
        let r = NiftiReader::new("test.nii");
        assert_eq!(r.n_voxels(), 64 * 64 * 32);
    }
    #[test]
    fn test_nifti_reader_init_data() {
        let mut r = NiftiReader::new("test.nii");
        r.init_data();
        assert_eq!(r.data.len(), r.n_voxels());
    }
    #[test]
    fn test_nifti_reader_voxel_to_world() {
        let r = NiftiReader::new("test.nii");
        let w = r.voxel_to_world(1.0, 2.0, 3.0);
        assert!((w[0] - 1.0).abs() < 1e-10);
        assert!((w[1] - 2.0).abs() < 1e-10);
        assert!((w[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_nifti_writer_build_header() {
        let writer = NiftiWriter::new("out.nii");
        let reader = NiftiReader::new("test.nii");
        let header = writer.build_header(&reader);
        assert_eq!(header.len(), 348);
        assert_eq!(header[344], b'n');
    }
    #[test]
    fn test_mha_mhd_format_new() {
        let m = MhaMhdFormat::new(true);
        assert_eq!(m.get("ObjectType"), Some("Image"));
        assert!(m.is_mha);
    }
    #[test]
    fn test_mha_mhd_format_dimensions() {
        let mut m = MhaMhdFormat::new(false);
        m.set("DimSize", "64 64 32");
        let dims = m.dimensions().unwrap();
        assert_eq!(dims, vec![64, 64, 32]);
    }
    #[test]
    fn test_mha_mhd_format_serialize() {
        let m = MhaMhdFormat::new(true);
        let s = m.serialize_header();
        assert!(s.contains("ElementDataFile = LOCAL"));
    }
    #[test]
    fn test_vtk_image_data_new() {
        let v = VtkImageData::new([10, 10, 10], [0.0; 3], [1.0; 3]);
        assert_eq!(v.n_voxels(), 1000);
    }
    #[test]
    fn test_vtk_image_data_index_to_world() {
        let v = VtkImageData::new([10, 10, 10], [1.0, 2.0, 3.0], [0.5, 0.5, 0.5]);
        let w = v.index_to_world(2, 4, 6);
        assert!((w[0] - 2.0).abs() < 1e-10);
        assert!((w[1] - 4.0).abs() < 1e-10);
        assert!((w[2] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_vtk_image_data_add_scalar() {
        let mut v = VtkImageData::new([4, 4, 4], [0.0; 3], [1.0; 3]);
        v.add_scalar("density", vec![1.0; 64]);
        assert!(v.scalar_arrays.contains_key("density"));
    }
    #[test]
    fn test_medical_mesh_new() {
        let m = MedicalMesh::new("bone");
        assert_eq!(m.name, "bone");
        assert_eq!(m.n_triangles(), 0);
    }
    #[test]
    fn test_medical_mesh_add_triangle() {
        let mut m = MedicalMesh::new("test");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(m.n_triangles(), 1);
        assert_eq!(m.vertices.len(), 3);
    }
    #[test]
    fn test_medical_mesh_surface_area() {
        let mut m = MedicalMesh::new("test");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let area = m.surface_area();
        assert!((area - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_landmark_set_new() {
        let ls = LandmarkSet::new("RAS");
        assert!(ls.landmarks.is_empty());
    }
    #[test]
    fn test_landmark_set_add_find() {
        let mut ls = LandmarkSet::new("RAS");
        ls.add(Landmark::new("AC", [0.0, 20.0, 0.0], 1.0));
        assert!(ls.find("AC").is_some());
        assert!(ls.find("PC").is_none());
    }
    #[test]
    fn test_landmark_set_centroid() {
        let mut ls = LandmarkSet::new("RAS");
        ls.add(Landmark::new("A", [0.0, 0.0, 0.0], 0.5));
        ls.add(Landmark::new("B", [2.0, 0.0, 0.0], 0.5));
        let c = ls.centroid().unwrap();
        assert!((c[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_segmentation_mask_new() {
        let seg = SegmentationMask::new([8, 8, 8], [1.0; 3]);
        assert_eq!(seg.count_label(0), 512);
    }
    #[test]
    fn test_segmentation_mask_set_get() {
        let mut seg = SegmentationMask::new([4, 4, 4], [1.0; 3]);
        seg.set_label(1, 2, 3, 2);
        assert_eq!(seg.get_label(1, 2, 3), Some(2));
    }
    #[test]
    fn test_segmentation_mask_volume() {
        let mut seg = SegmentationMask::new([4, 4, 4], [2.0; 3]);
        seg.set_label(0, 0, 0, 1);
        seg.set_label(1, 0, 0, 1);
        let vol = seg.volume_mm3(1);
        assert!((vol - 16.0).abs() < 1e-10);
    }
    #[test]
    fn test_segmentation_mask_bounding_box() {
        let mut seg = SegmentationMask::new([8, 8, 8], [1.0; 3]);
        seg.set_label(2, 3, 4, 1);
        let bb = seg.bounding_box(1).unwrap();
        assert_eq!(bb[0], 2);
        assert_eq!(bb[2], 3);
        assert_eq!(bb[4], 4);
    }
    #[test]
    fn test_dose_volume_new() {
        let dv = DoseVolume::new([10, 10, 10], [2.0; 3], 60.0);
        assert_eq!(dv.dose.len(), 1000);
        assert_eq!(dv.prescription_dose, 60.0);
    }
    #[test]
    fn test_dose_volume_set_get() {
        let mut dv = DoseVolume::new([4, 4, 4], [1.0; 3], 50.0);
        dv.set_dose(1, 2, 3, 45.0);
        assert!((dv.get_dose(1, 2, 3).unwrap() - 45.0).abs() < 1e-10);
    }
    #[test]
    fn test_dose_volume_max_dose() {
        let mut dv = DoseVolume::new([2, 2, 2], [1.0; 3], 60.0);
        dv.set_dose(0, 0, 0, 55.0);
        dv.set_dose(1, 1, 1, 62.0);
        assert!((dv.max_dose() - 62.0).abs() < 1e-10);
    }
    #[test]
    fn test_dose_volume_dvh() {
        let mut dv = DoseVolume::new([4, 4, 1], [1.0; 3], 60.0);
        for i in 0..16 {
            dv.dose[i] = i as f64 * 4.0;
        }
        let (bins, vf) = dv.dvh(10);
        assert_eq!(bins.len(), 10);
        assert_eq!(vf.len(), 10);
        assert!(vf[0] >= vf[9]);
    }
    #[test]
    fn test_dose_volume_v_dose() {
        let mut dv = DoseVolume::new([2, 2, 1], [1.0; 3], 60.0);
        dv.dose = vec![10.0, 20.0, 30.0, 40.0];
        let v30 = dv.v_dose(30.0);
        assert!((v30 - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_dose_volume_d_volume() {
        let mut dv = DoseVolume::new([2, 2, 1], [1.0; 3], 60.0);
        dv.dose = vec![10.0, 20.0, 30.0, 40.0];
        let d50 = dv.d_volume(0.5);
        assert!(d50 >= 20.0 && d50 <= 30.0);
    }
    #[test]
    fn test_dicom_vr_us() {
        let vr = DicomVr::Us(512);
        assert_eq!(vr.vr_code(), "US");
        assert_eq!(vr.as_us(), Some(512));
    }
    #[test]
    fn test_dicom_vr_ds() {
        let vr = DicomVr::Ds(3.125);
        assert_eq!(vr.vr_code(), "DS");
        assert!((vr.as_ds().unwrap() - 3.125).abs() < 1e-10);
    }
    #[test]
    fn test_dicom_vr_lo() {
        let vr = DicomVr::Lo("TestHospital".to_string());
        assert_eq!(vr.as_str(), Some("TestHospital"));
    }
    #[test]
    fn test_dicom_vr_ui() {
        let vr = DicomVr::Ui("1.2.840.10008.5.1.4.1.1.2".to_string());
        assert_eq!(vr.vr_code(), "UI");
    }
    #[test]
    fn test_dicom_dataset_set_get() {
        let mut ds = DicomDataSet::new();
        let tag = DicomTag::new(0x0028, 0x0010);
        ds.set(DicomElement::new(tag.clone(), DicomVr::Us(512)));
        assert_eq!(ds.get_us(&tag), Some(512));
    }
    #[test]
    fn test_dicom_dataset_roundtrip_bytes() {
        let mut ds = DicomDataSet::new();
        ds.set(DicomElement::new(
            DicomTag::new(0x0028, 0x0010),
            DicomVr::Us(256),
        ));
        ds.set(DicomElement::new(
            DicomTag::new(0x0008, 0x0060),
            DicomVr::Lo("CT".to_string()),
        ));
        let bytes = ds.to_bytes();
        let ds2 = DicomDataSet::parse_bytes(&bytes);
        assert_eq!(ds2.get_us(&DicomTag::new(0x0028, 0x0010)), Some(256));
        assert_eq!(ds2.get_str(&DicomTag::new(0x0008, 0x0060)), Some("CT"));
    }
    #[test]
    fn test_apply_windowing_basic() {
        let pixels = vec![1024u16, 2048, 0];
        let display = apply_windowing(&pixels, 1.0, -1024.0, 40.0, 400.0);
        assert_eq!(display.len(), 3);
        assert!(display[0] > 0 && display[0] < 255);
        assert_eq!(display[2], 0);
    }
    #[test]
    fn test_nifti_transform_identity() {
        let t = NiftiTransform::identity();
        let w = t.sform_to_world(1.0, 2.0, 3.0);
        assert!((w[0] - 1.0).abs() < 1e-10);
        assert!((w[1] - 2.0).abs() < 1e-10);
        assert!((w[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_nifti_qform_code_roundtrip() {
        let code = NiftiQformCode::Mni152;
        assert_eq!(
            NiftiQformCode::from_code(code.as_code()),
            NiftiQformCode::Mni152
        );
    }
    #[test]
    fn test_volume_reconstructor_add_slice() {
        let mut rec = VolumeReconstructor::new(4, 4, [1.0, 1.0], 2.0);
        let slice: Vec<u16> = (0..16).map(|i| i as u16).collect();
        rec.add_slice(&slice);
        rec.add_slice(&slice);
        assert_eq!(rec.n_slices, 2);
        assert_eq!(rec.voxel(0, 0, 0), Some(0));
        assert_eq!(rec.voxel(1, 2, 1), Some(6));
    }
    #[test]
    fn test_volume_reconstructor_physical_size() {
        let mut rec = VolumeReconstructor::new(10, 20, [0.5, 0.5], 1.5);
        rec.add_slice(&vec![0u16; 200]);
        let sz = rec.physical_size();
        assert!((sz[0] - 5.0).abs() < 1e-10);
        assert!((sz[1] - 10.0).abs() < 1e-10);
        assert!((sz[2] - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_resample_trilinear_identity() {
        let src: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let dims = [2usize, 2, 2];
        let out = resample_trilinear(&src, dims, dims);
        assert_eq!(out.len(), 8);
        for (a, b) in src.iter().zip(out.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn test_resample_trilinear_upscale() {
        let src = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let out = resample_trilinear(&src, [2, 2, 2], [3, 3, 3]);
        assert_eq!(out.len(), 27);
        assert!((out[0] - 0.0).abs() < 1e-8);
        assert!((out[26] - 7.0).abs() < 1e-8);
    }
    #[test]
    fn test_stl_exporter_round_trip() {
        let v0 = [0.0f32, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let normal = StlTriangle::compute_normal(v0, v1, v2);
        let tri = StlTriangle {
            normal,
            vertices: [v0, v1, v2],
        };
        let bytes = StlExporter::to_binary(&[tri]);
        let back = StlExporter::from_binary(&bytes);
        assert_eq!(back.len(), 1);
        assert!((back[0].vertices[1][0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_ct_point_cloud_threshold() {
        let pixels: Vec<u16> = vec![0, 500, 1000, 1500];
        let pts =
            ct_point_cloud_threshold(&pixels, [2, 2, 1], [1.0; 3], [0.0; 3], 1.0, -1024.0, 0.0);
        assert_eq!(pts.len(), 1);
    }
    #[test]
    fn test_nrrd_reader_parse_header() {
        let header = "type: float\ndimension: 3\nsizes: 64 64 32\nencoding: raw\n";
        let mut r = NrrdReader::new();
        r.parse_header(header);
        assert_eq!(r.sizes, vec![64, 64, 32]);
        assert_eq!(r.encoding, NrrdEncoding::Raw);
    }
    #[test]
    fn test_nrrd_reader_load_text() {
        let mut r = NrrdReader::new();
        r.parse_header("sizes: 2 2\nencoding: text\n");
        r.load_text("1.0 2.0 3.0 4.0");
        assert_eq!(r.data.len(), 4);
        assert!((r.data[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_minc_dimension_world_coord() {
        let dim = MincDimension::new("xspace", 10, 0.5, -2.0);
        assert!((dim.world_coord(0) - (-2.0)).abs() < 1e-10);
        assert!((dim.world_coord(4) - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_minc_volume_index3() {
        let dims = vec![
            MincDimension::new("xspace", 4, 1.0, 0.0),
            MincDimension::new("yspace", 4, 1.0, 0.0),
            MincDimension::new("zspace", 4, 1.0, 0.0),
        ];
        let vol = MincVolume::new(dims);
        assert_eq!(vol.index3(0, 0, 0), Some(0));
        assert_eq!(vol.index3(1, 0, 0), Some(16));
        assert!(vol.index3(5, 0, 0).is_none());
    }
    #[test]
    fn test_fiber_tract_arc_length() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let fiber = FiberTract::new(pts);
        assert!((fiber.arc_length() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_tract_vtk_writer_format() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let fiber = FiberTract::new(pts);
        let writer = TractVtkWriter::new("test_tracts");
        let vtk = writer.to_vtk_string(&[fiber]);
        assert!(vtk.contains("POLYDATA"));
        assert!(vtk.contains("LINES"));
        assert!(vtk.contains("POINTS 2"));
    }
}
