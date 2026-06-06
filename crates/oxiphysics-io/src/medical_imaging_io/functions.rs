//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Parses raw voxel bytes into f32 based on NIfTI datatype code.
pub(super) fn parse_nifti_voxels(data: &[u8], datatype: i16, num_voxels: usize) -> Vec<f32> {
    match datatype {
        2 => data.iter().take(num_voxels).map(|&b| b as f32).collect(),
        4 => {
            let count = num_voxels.min(data.len() / 2);
            (0..count)
                .map(|i| read_i16_le(data, i * 2) as f32)
                .collect()
        }
        16 => {
            let count = num_voxels.min(data.len() / 4);
            (0..count).map(|i| read_f32_le(data, i * 4)).collect()
        }
        64 => {
            let count = num_voxels.min(data.len() / 8);
            (0..count)
                .map(|i| {
                    let bytes: [u8; 8] = data[i * 8..i * 8 + 8].try_into().unwrap_or([0; 8]);
                    f64::from_le_bytes(bytes) as f32
                })
                .collect()
        }
        _ => vec![0.0; num_voxels],
    }
}
/// Reads a little-endian u16 from a byte slice.
pub(super) fn read_u16_le(data: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes([data[offset], data[offset + 1]])
}
/// Reads a little-endian u32 from a byte slice.
pub(super) fn read_u32_le(data: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ])
}
/// Reads a little-endian i16 from a byte slice.
pub(super) fn read_i16_le(data: &[u8], offset: usize) -> i16 {
    i16::from_le_bytes([data[offset], data[offset + 1]])
}
/// Reads a little-endian i32 from a byte slice.
pub(super) fn read_i32_le(data: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ])
}
/// Reads a little-endian f32 from a byte slice.
pub(super) fn read_f32_le(data: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ])
}
/// Converts flat index to (x, y, z) coordinates.
pub(super) fn index_to_coords(idx: usize, dims: [usize; 3]) -> (usize, usize, usize) {
    let xy = dims[0] * dims[1];
    let z = idx / xy;
    let rem = idx % xy;
    let y = rem / dims[0];
    let x = rem % dims[0];
    (x, y, z)
}
/// Converts (x, y, z) coordinates to a flat index.
pub(super) fn coords_to_index(x: usize, y: usize, z: usize, dims: [usize; 3]) -> usize {
    z * dims[0] * dims[1] + y * dims[0] + x
}
/// Builds a 1D Gaussian kernel with the given sigma and half-width.
pub(super) fn build_gaussian_kernel(sigma: f64, radius: usize) -> Vec<f32> {
    let size = 2 * radius + 1;
    let mut kernel = vec![0.0f32; size];
    let s2 = 2.0 * sigma * sigma;
    for (i, k) in kernel.iter_mut().enumerate() {
        let x = i as f64 - radius as f64;
        *k = (-x * x / s2).exp() as f32;
    }
    let sum: f32 = kernel.iter().sum();
    if sum > 0.0 {
        for k in &mut kernel {
            *k /= sum;
        }
    }
    kernel
}
/// Builds a synthetic DICOM byte buffer for testing.
///
/// Creates a minimal valid DICOM file with preamble, magic number,
/// and the specified data elements.
#[cfg(test)]
pub(super) fn build_test_dicom(elements: &[(u16, u16, &str, &[u8])]) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(&[0u8; 128]);
    buf.extend_from_slice(b"DICM");
    for &(group, element, vr, value) in elements {
        buf.extend_from_slice(&group.to_le_bytes());
        buf.extend_from_slice(&element.to_le_bytes());
        buf.extend_from_slice(vr.as_bytes());
        let padded_len = if value.len() % 2 == 1 {
            value.len() + 1
        } else {
            value.len()
        };
        buf.extend_from_slice(&(padded_len as u16).to_le_bytes());
        buf.extend_from_slice(value);
        if value.len() % 2 == 1 {
            buf.push(b' ');
        }
    }
    buf
}
/// Builds a synthetic NIfTI-1 byte buffer for testing.
#[cfg(test)]
pub(super) fn build_test_nifti(dims: [i16; 3], pixdim: [f32; 3], voxel_data: &[f32]) -> Vec<u8> {
    let mut buf = vec![0u8; 352];
    let hdr_size: i32 = 348;
    buf[0..4].copy_from_slice(&hdr_size.to_le_bytes());
    let ndim: i16 = 3;
    buf[40..42].copy_from_slice(&ndim.to_le_bytes());
    buf[42..44].copy_from_slice(&dims[0].to_le_bytes());
    buf[44..46].copy_from_slice(&dims[1].to_le_bytes());
    buf[46..48].copy_from_slice(&dims[2].to_le_bytes());
    let dt: i16 = 16;
    buf[70..72].copy_from_slice(&dt.to_le_bytes());
    let bp: i16 = 32;
    buf[72..74].copy_from_slice(&bp.to_le_bytes());
    let one_f: f32 = 1.0;
    buf[76..80].copy_from_slice(&one_f.to_le_bytes());
    buf[80..84].copy_from_slice(&pixdim[0].to_le_bytes());
    buf[84..88].copy_from_slice(&pixdim[1].to_le_bytes());
    buf[88..92].copy_from_slice(&pixdim[2].to_le_bytes());
    let vox_off: f32 = 352.0;
    buf[108..112].copy_from_slice(&vox_off.to_le_bytes());
    buf[112..116].copy_from_slice(&one_f.to_le_bytes());
    buf[344] = b'n';
    buf[345] = b'+';
    buf[346] = b'1';
    for v in voxel_data {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    buf
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::medical_imaging_io::types::*;
    #[test]
    fn test_voxel_volume_creation() {
        let vol = VoxelVolume::new([4, 5, 6], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        assert_eq!(vol.num_voxels(), 120);
        assert_eq!(vol.dims, [4, 5, 6]);
    }
    #[test]
    fn test_voxel_volume_get_set() {
        let mut vol = VoxelVolume::new([3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        vol.set(1, 2, 1, 42.0);
        assert_eq!(vol.get(1, 2, 1), Some(42.0));
        assert_eq!(vol.get(0, 0, 0), Some(0.0));
        assert_eq!(vol.get(3, 0, 0), None);
    }
    #[test]
    fn test_voxel_volume_from_data() {
        let data = vec![1.0; 8];
        let vol = VoxelVolume::from_data(data.clone(), [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(vol.is_some());
        let bad = VoxelVolume::from_data(data, [3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(bad.is_none());
    }
    #[test]
    fn test_voxel_volume_coords_roundtrip() {
        let vol = VoxelVolume::new([5, 6, 7], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        for z in 0..7 {
            for y in 0..6 {
                for x in 0..5 {
                    let idx = vol.index(x, y, z).unwrap();
                    let (rx, ry, rz) = vol.coords_from_index(idx).unwrap();
                    assert_eq!((x, y, z), (rx, ry, rz));
                }
            }
        }
    }
    #[test]
    fn test_voxel_to_world() {
        let vol = VoxelVolume::new([10, 10, 10], [0.5, 0.5, 2.0], [10.0, 20.0, 30.0]);
        let w = vol.voxel_to_world(4.0, 6.0, 2.0);
        assert!((w[0] - 12.0).abs() < 1e-10);
        assert!((w[1] - 23.0).abs() < 1e-10);
        assert!((w[2] - 34.0).abs() < 1e-10);
    }
    #[test]
    fn test_voxel_volume_stats() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        assert_eq!(vol.min_value(), 1.0);
        assert_eq!(vol.max_value(), 8.0);
        assert!((vol.mean_value() - 4.5).abs() < 1e-10);
    }
    #[test]
    fn test_dicom_tag_parsing() {
        let tag = DicomTag::new(0x0010, 0x0010, ValueRepresentation::PN);
        assert_eq!(tag.tag_string(), "(0010,0010)");
        assert!(!tag.is_private());
        assert!(!tag.is_group_length());
    }
    #[test]
    fn test_dicom_tag_private() {
        let tag = DicomTag::new(0x0009, 0x0001, ValueRepresentation::UN);
        assert!(tag.is_private());
    }
    #[test]
    fn test_vr_from_bytes() {
        assert_eq!(
            ValueRepresentation::from_bytes(b'P', b'N'),
            Some(ValueRepresentation::PN)
        );
        assert_eq!(
            ValueRepresentation::from_bytes(b'C', b'S'),
            Some(ValueRepresentation::CS)
        );
        assert_eq!(ValueRepresentation::from_bytes(b'X', b'X'), None);
    }
    #[test]
    fn test_vr_extended_length() {
        assert!(ValueRepresentation::OB.has_extended_length());
        assert!(ValueRepresentation::OW.has_extended_length());
        assert!(ValueRepresentation::SQ.has_extended_length());
        assert!(!ValueRepresentation::US.has_extended_length());
        assert!(!ValueRepresentation::CS.has_extended_length());
    }
    #[test]
    fn test_dicom_element_value_string() {
        let elem = DicomElement::new(
            DicomTag::new(0x0010, 0x0010, ValueRepresentation::PN),
            b"Doe^John \0".to_vec(),
        );
        assert_eq!(elem.value_as_string(), Some("Doe^John".to_string()));
    }
    #[test]
    fn test_dicom_element_value_u16() {
        let elem = DicomElement::new(
            DicomTag::new(0x0028, 0x0010, ValueRepresentation::US),
            512u16.to_le_bytes().to_vec(),
        );
        assert_eq!(elem.value_as_u16(), Some(512));
    }
    #[test]
    fn test_dicom_element_value_ds() {
        let elem = DicomElement::new(
            DicomTag::new(0x0018, 0x0050, ValueRepresentation::DS),
            b"2.500 ".to_vec(),
        );
        assert!((elem.value_as_ds().unwrap() - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_dicom_reader_parse() {
        let dicom = build_test_dicom(&[
            (0x0008, 0x0060, "CS", b"CT"),
            (0x0010, 0x0010, "PN", b"Smith^Alice"),
            (0x0028, 0x0010, "US", &512u16.to_le_bytes()),
            (0x0028, 0x0011, "US", &256u16.to_le_bytes()),
        ]);
        let mut reader = DicomReader::new();
        let header = reader.parse(&dicom).unwrap();
        assert_eq!(header.modality, "CT");
        assert_eq!(header.patient_name, "Smith^Alice");
        assert_eq!(header.rows, 512);
        assert_eq!(header.columns, 256);
        assert_eq!(reader.element_count(), 4);
    }
    #[test]
    fn test_dicom_reader_too_short() {
        let mut reader = DicomReader::new();
        assert!(reader.parse(&[0u8; 50]).is_err());
    }
    #[test]
    fn test_dicom_reader_bad_magic() {
        let mut data = vec![0u8; 136];
        data[128..132].copy_from_slice(b"XXXX");
        let mut reader = DicomReader::new();
        assert!(reader.parse(&data).is_err());
    }
    #[test]
    fn test_nifti_header_defaults() {
        let hdr = NiftiHeader::new();
        assert_eq!(hdr.sizeof_hdr, 348);
        assert_eq!(hdr.ndim(), 3);
        assert_eq!(hdr.volume_dims(), [1, 1, 1]);
        assert_eq!(hdr.magic, "n+1");
    }
    #[test]
    fn test_nifti_header_dim_field() {
        let mut hdr = NiftiHeader::new();
        hdr.dim = [3, 64, 64, 32, 0, 0, 0, 0];
        assert_eq!(hdr.volume_dims(), [64, 64, 32]);
        assert_eq!(hdr.num_voxels(), 64 * 64 * 32);
    }
    #[test]
    fn test_nifti_sform_affine() {
        let mut hdr = NiftiHeader::new();
        hdr.srow_x = [2.0, 0.0, 0.0, -100.0];
        hdr.srow_y = [0.0, 2.0, 0.0, -120.0];
        hdr.srow_z = [0.0, 0.0, 3.0, -50.0];
        let aff = hdr.sform_affine();
        assert!((aff[0][0] - 2.0).abs() < 1e-6);
        assert!((aff[2][2] - 3.0).abs() < 1e-6);
        assert!((aff[3][3] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_nifti_reader_parse() {
        let voxels: Vec<f32> = (0..27).map(|i| i as f32).collect();
        let buf = build_test_nifti([3, 3, 3], [1.0, 1.0, 1.0], &voxels);
        let mut reader = NiftiReader::new();
        reader.parse(&buf).unwrap();
        assert_eq!(reader.header.volume_dims(), [3, 3, 3]);
        assert_eq!(reader.voxels.len(), 27);
        assert!((reader.voxels[0] - 0.0).abs() < 1e-6);
        assert!((reader.voxels[26] - 26.0).abs() < 1e-6);
    }
    #[test]
    fn test_nifti_reader_to_volume() {
        let voxels: Vec<f32> = (0..8).map(|i| i as f32 * 10.0).collect();
        let buf = build_test_nifti([2, 2, 2], [0.5, 0.5, 2.0], &voxels);
        let mut reader = NiftiReader::new();
        reader.parse(&buf).unwrap();
        let vol = reader.to_volume();
        assert_eq!(vol.dims, [2, 2, 2]);
        assert!((vol.spacing[2] - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_flood_fill_single_component() {
        let mut vol = VoxelVolume::new([5, 5, 5], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        for z in 1..4 {
            for y in 1..4 {
                for x in 1..4 {
                    vol.set(x, y, z, 100.0);
                }
            }
        }
        let seg = Segmentation::from_threshold(&vol, 50.0);
        assert_eq!(seg.num_components, 1);
        assert_eq!(seg.component_size(1), 27);
    }
    #[test]
    fn test_flood_fill_two_components() {
        let mut vol = VoxelVolume::new([10, 1, 1], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        vol.set(0, 0, 0, 100.0);
        vol.set(1, 0, 0, 100.0);
        vol.set(5, 0, 0, 200.0);
        vol.set(6, 0, 0, 200.0);
        vol.set(7, 0, 0, 200.0);
        let seg = Segmentation::from_threshold(&vol, 50.0);
        assert_eq!(seg.num_components, 2);
        assert_eq!(seg.component_size(1), 2);
        assert_eq!(seg.component_size(2), 3);
    }
    #[test]
    fn test_label_statistics() {
        let mut vol = VoxelVolume::new([4, 4, 4], [2.0, 2.0, 2.0], [0.0, 0.0, 0.0]);
        for z in 0..2 {
            for y in 0..2 {
                for x in 0..2 {
                    vol.set(x, y, z, 100.0);
                }
            }
        }
        let seg = Segmentation::from_threshold(&vol, 50.0);
        let stats = seg.label_statistics(&vol);
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].voxel_count, 8);
        assert!((stats[0].volume_mm3 - 64.0).abs() < 1e-10);
        assert!((stats[0].mean_intensity - 100.0).abs() < 1e-10);
        assert!((stats[0].centroid[0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_trilinear_interpolation_at_corners() {
        let mut vol = VoxelVolume::new([3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        vol.set(0, 0, 0, 0.0);
        vol.set(1, 0, 0, 10.0);
        vol.set(0, 1, 0, 20.0);
        vol.set(1, 1, 0, 30.0);
        vol.set(0, 0, 1, 40.0);
        vol.set(1, 0, 1, 50.0);
        vol.set(0, 1, 1, 60.0);
        vol.set(1, 1, 1, 70.0);
        assert!(
            (VolumeProcessing::trilinear_interpolate(&vol, 0.0, 0.0, 0.0).unwrap() - 0.0).abs()
                < 1e-6
        );
        assert!(
            (VolumeProcessing::trilinear_interpolate(&vol, 1.0, 0.0, 0.0).unwrap() - 10.0).abs()
                < 1e-6
        );
        assert!(
            (VolumeProcessing::trilinear_interpolate(&vol, 0.5, 0.0, 0.0).unwrap() - 5.0).abs()
                < 1e-6
        );
    }
    #[test]
    fn test_trilinear_interpolation_midpoint() {
        let data = vec![10.0; 8];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let val = VolumeProcessing::trilinear_interpolate(&vol, 0.5, 0.5, 0.5).unwrap();
        assert!((val - 10.0).abs() < 1e-6);
    }
    #[test]
    fn test_trilinear_out_of_bounds() {
        let vol = VoxelVolume::new([3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        assert!(VolumeProcessing::trilinear_interpolate(&vol, -0.1, 0.0, 0.0).is_none());
        assert!(VolumeProcessing::trilinear_interpolate(&vol, 2.0, 0.0, 0.0).is_none());
    }
    #[test]
    fn test_hounsfield_unit_conversion() {
        let data = vec![1000.0, 2000.0, 0.0, 500.0];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 1], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let hu = VolumeProcessing::to_hounsfield_units(&vol, 1.0, -1024.0);
        assert!((hu.data[0] - (-24.0)).abs() < 1e-3);
        assert!((hu.data[1] - 976.0).abs() < 1e-3);
        assert!((hu.data[2] - (-1024.0)).abs() < 1e-3);
    }
    #[test]
    fn test_window_level() {
        let data = vec![-1000.0, 0.0, 40.0, 400.0];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 1], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let wl = VolumeProcessing::window_level(&vol, 40.0, 400.0);
        assert!((wl.data[0] - 0.0).abs() < 1e-6);
        assert!((wl.data[3] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_gaussian_smooth_uniform() {
        let data = vec![5.0; 27];
        let vol =
            VoxelVolume::from_data(data, [3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let smoothed = VolumeProcessing::gaussian_smooth(&vol, 1.0, 1);
        for &v in &smoothed.data {
            assert!((v - 5.0).abs() < 1e-4);
        }
    }
    #[test]
    fn test_histogram_equalize() {
        let data: Vec<f32> = (0..64).map(|i| i as f32).collect();
        let vol =
            VoxelVolume::from_data(data, [4, 4, 4], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let eq = VolumeProcessing::histogram_equalize(&vol, 64, 255.0);
        assert!(eq.min_value() >= -0.1);
        assert!(eq.max_value() <= 255.1);
    }
    #[test]
    fn test_threshold() {
        let data = vec![0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 25.0, 75.0];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let th = VolumeProcessing::threshold(&vol, 100.0);
        assert_eq!(th.data[0], 0.0);
        assert_eq!(th.data[2], 1.0);
        assert_eq!(th.data[3], 1.0);
    }
    #[test]
    fn test_histogram() {
        let data = vec![0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0];
        let vol =
            VoxelVolume::from_data(data, [2, 2, 2], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]).unwrap();
        let (centers, counts) = VolumeProcessing::histogram(&vol, 4);
        assert_eq!(centers.len(), 4);
        assert_eq!(counts.len(), 4);
        let total: u64 = counts.iter().sum();
        assert_eq!(total, 8);
    }
    #[test]
    fn test_gradient_magnitude() {
        let mut vol = VoxelVolume::new([5, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..5 {
                    vol.set(x, y, z, x as f32 * 10.0);
                }
            }
        }
        let grad = VolumeProcessing::gradient_magnitude(&vol);
        let val = grad.get(2, 1, 1).unwrap();
        assert!((val - 10.0).abs() < 1e-3);
    }
    #[test]
    fn test_display_formats() {
        let vol = VoxelVolume::new([64, 64, 32], [0.5, 0.5, 2.0], [0.0, 0.0, 0.0]);
        let s = format!("{}", vol);
        assert!(s.contains("64x64x32"));
        let hdr = NiftiHeader::new();
        let s = format!("{}", hdr);
        assert!(s.contains("n+1"));
        let tag = DicomTag::new(0x0010, 0x0010, ValueRepresentation::PN);
        let s = format!("{}", tag);
        assert!(s.contains("0010"));
    }
}
