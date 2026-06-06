// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Medical imaging format utilities (DICOM-inspired and NIfTI).
//!
//! Provides structures for DICOM tag/dataset handling, voxel volumes,
//! segmentation masks, NIfTI-1 style header I/O, and MRI phantom signal
//! simulation.

use std::collections::HashMap;
use std::io::{Read as IoRead, Write as IoWrite};

use crate::{Error, Result};

// ── DICOM helpers ─────────────────────────────────────────────────────────────

/// A single DICOM data element.
#[derive(Debug, Clone)]
pub struct DicomTag {
    /// Group number (first half of the tag address, e.g. `0x0008`).
    pub group: u16,
    /// Element number (second half, e.g. `0x0060`).
    pub element: u16,
    /// Two-character Value Representation code, e.g. `"CS"`, `"DS"`, `"UI"`.
    pub value_representation: String,
    /// Raw byte data of this element.
    pub data: Vec<u8>,
}

impl DicomTag {
    /// Create a new DICOM tag from raw bytes.
    pub fn new(group: u16, element: u16, vr: impl Into<String>, data: Vec<u8>) -> Self {
        Self {
            group,
            element,
            value_representation: vr.into(),
            data,
        }
    }

    /// Create a tag whose data is a UTF-8 string.
    pub fn from_str(group: u16, element: u16, vr: impl Into<String>, value: &str) -> Self {
        Self::new(group, element, vr, value.as_bytes().to_vec())
    }

    /// Create a tag whose data is a little-endian f64 (8 bytes).
    pub fn from_f64(group: u16, element: u16, vr: impl Into<String>, value: f64) -> Self {
        Self::new(group, element, vr, value.to_le_bytes().to_vec())
    }
}

/// A DICOM dataset — a collection of [`DicomTag`]s keyed by `(group, element)`.
#[derive(Debug, Clone, Default)]
pub struct DicomDataset {
    /// Internal tag store.
    pub tags: HashMap<(u16, u16), DicomTag>,
}

impl DicomDataset {
    /// Create an empty dataset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert or replace a tag.
    pub fn insert(&mut self, tag: DicomTag) {
        self.tags.insert((tag.group, tag.element), tag);
    }

    /// Retrieve the raw UTF-8 string value for the given `(group, element)`.
    ///
    /// Returns `None` if the tag is absent or the data is not valid UTF-8.
    pub fn get_string(&self, group: u16, elem: u16) -> Option<String> {
        let tag = self.tags.get(&(group, elem))?;
        String::from_utf8(tag.data.clone()).ok()
    }

    /// Retrieve a floating-point value for the given `(group, element)`.
    ///
    /// For VR codes `"FD"` or `"FL"` (binary float VRs) the data is
    /// interpreted as a little-endian f64 (8 bytes).  For all other VRs the
    /// raw bytes are decoded as a UTF-8 decimal string.  If the string parse
    /// fails and the data is exactly 8 bytes the binary fallback is tried.
    pub fn get_f64(&self, group: u16, elem: u16) -> Option<f64> {
        let tag = self.tags.get(&(group, elem))?;
        let vr = tag.value_representation.as_str();
        if vr == "FD" || vr == "FL" {
            // Binary floating-point VR
            if tag.data.len() == 8 {
                let bytes: [u8; 8] = tag.data[..8].try_into().ok()?;
                return Some(f64::from_le_bytes(bytes));
            }
        }
        // String-encoded decimal (DS, IS, …)
        if let Ok(s) = String::from_utf8(tag.data.clone())
            && let Ok(v) = s.trim().parse::<f64>()
        {
            return Some(v);
        }
        // Last-resort binary fallback
        if tag.data.len() == 8 {
            let bytes: [u8; 8] = tag.data[..8].try_into().ok()?;
            return Some(f64::from_le_bytes(bytes));
        }
        None
    }
}

// ── Voxel volume ──────────────────────────────────────────────────────────────

/// A 3-D voxel volume with uniform spacing, storing 16-bit signed pixel data.
#[derive(Debug, Clone)]
pub struct VoxelVolume {
    /// Dimensions `[width, height, depth]` in voxels.
    pub dimensions: [usize; 3],
    /// Physical size of one voxel in mm, `[dx, dy, dz]`.
    pub voxel_spacing: [f64; 3],
    /// Raw pixel values in row-major order (width-first).
    pub pixel_data: Vec<i16>,
}

impl VoxelVolume {
    /// Create a new voxel volume filled with zeros.
    pub fn new(dimensions: [usize; 3], voxel_spacing: [f64; 3]) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            dimensions,
            voxel_spacing,
            pixel_data: vec![0; n],
        }
    }

    /// Convert a stored pixel value to Hounsfield Units (HU).
    ///
    /// HU = `pixel * slope + intercept`.
    pub fn to_hounsfield(pixel: i16, slope: f64, intercept: f64) -> f64 {
        pixel as f64 * slope + intercept
    }

    /// Number of voxels in the volume.
    pub fn voxel_count(&self) -> usize {
        self.dimensions[0] * self.dimensions[1] * self.dimensions[2]
    }

    /// Physical volume in mm³.
    pub fn physical_volume_mm3(&self) -> f64 {
        self.voxel_count() as f64
            * self.voxel_spacing[0]
            * self.voxel_spacing[1]
            * self.voxel_spacing[2]
    }

    /// Get voxel at `(x, y, z)`.  Returns `None` if indices are out of range.
    pub fn get(&self, x: usize, y: usize, z: usize) -> Option<i16> {
        if x < self.dimensions[0] && y < self.dimensions[1] && z < self.dimensions[2] {
            Some(
                self.pixel_data
                    [z * self.dimensions[1] * self.dimensions[0] + y * self.dimensions[0] + x],
            )
        } else {
            None
        }
    }
}

// ── Segmentation ──────────────────────────────────────────────────────────────

/// A voxel-wise segmentation mask aligned to a [`VoxelVolume`].
#[derive(Debug, Clone)]
pub struct Segmentation {
    /// Label for each voxel in the same linear order as the parent volume.
    pub labels: Vec<u8>,
    /// Number of distinct classes (not counting background 0).
    pub n_classes: usize,
}

impl Segmentation {
    /// Create a segmentation mask of `n_voxels` voxels with `n_classes` classes.
    pub fn new(n_voxels: usize, n_classes: usize) -> Self {
        Self {
            labels: vec![0; n_voxels],
            n_classes,
        }
    }

    /// Count voxels with the given `label` and return the segmented volume in mm³.
    ///
    /// `spacing` is `[dx, dy, dz]` in mm.
    pub fn compute_volume(&self, label: u8, spacing: [f64; 3]) -> f64 {
        let count = self.labels.iter().filter(|&&l| l == label).count();
        count as f64 * spacing[0] * spacing[1] * spacing[2]
    }

    /// Fraction of all voxels assigned to `label`.
    pub fn label_fraction(&self, label: u8) -> f64 {
        if self.labels.is_empty() {
            return 0.0;
        }
        let count = self.labels.iter().filter(|&&l| l == label).count();
        count as f64 / self.labels.len() as f64
    }
}

// ── NIfTI header ──────────────────────────────────────────────────────────────

/// A minimal NIfTI-1 style header.
///
/// Only the fields needed for basic volume description are stored.
/// The on-disk format is a simple binary blob (not the full 348-byte NIfTI-1
/// standard) to keep the implementation self-contained.
#[derive(Debug, Clone)]
pub struct NiftiHeader {
    /// Dimensions of the volume: `dim[0]` is the rank (usually 3 or 4);
    /// `dim[1..=rank]` are the actual sizes.
    pub dim: [usize; 7],
    /// Voxel sizes along each dimension (mm or s for the time axis).
    pub pixdim: [f64; 7],
    /// NIfTI datatype code (e.g. 4 = INT16, 16 = FLOAT32).
    pub datatype: u16,
}

impl NiftiHeader {
    /// Create a default header for a 3-D INT16 volume.
    pub fn new_3d(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Self {
        Self {
            dim: [3, nx, ny, nz, 1, 1, 1],
            pixdim: [1.0, dx, dy, dz, 0.0, 0.0, 0.0],
            datatype: 4, // INT16
        }
    }

    /// Write the header to a binary file at `path`.
    ///
    /// Layout: 7×8 bytes (dims as u64 LE) + 7×8 bytes (pixdims as f64 LE)
    /// + 2 bytes (datatype as u16 LE) = 114 bytes total.
    pub fn write_header(&self, path: &str) -> Result<()> {
        let mut file = std::fs::File::create(path).map_err(Error::Io)?;
        for d in &self.dim {
            file.write_all(&(*d as u64).to_le_bytes())
                .map_err(Error::Io)?;
        }
        for p in &self.pixdim {
            file.write_all(&p.to_le_bytes()).map_err(Error::Io)?;
        }
        file.write_all(&self.datatype.to_le_bytes())
            .map_err(Error::Io)?;
        Ok(())
    }

    /// Read a header previously written by `write_header` from `path`.
    pub fn read_header(path: &str) -> Result<Self> {
        let mut file = std::fs::File::open(path).map_err(Error::Io)?;
        let mut buf = [0u8; 114];
        file.read_exact(&mut buf).map_err(Error::Io)?;
        let mut dim = [0usize; 7];
        for (i, d) in dim.iter_mut().enumerate() {
            let bytes: [u8; 8] = buf[i * 8..(i + 1) * 8]
                .try_into()
                .map_err(|_| Error::Parse("dim bytes".into()))?;
            *d = u64::from_le_bytes(bytes) as usize;
        }
        let offset = 7 * 8;
        let mut pixdim = [0f64; 7];
        for (i, p) in pixdim.iter_mut().enumerate() {
            let bytes: [u8; 8] = buf[offset + i * 8..offset + (i + 1) * 8]
                .try_into()
                .map_err(|_| Error::Parse("pixdim bytes".into()))?;
            *p = f64::from_le_bytes(bytes);
        }
        let dt_offset = offset + 7 * 8;
        let datatype = u16::from_le_bytes([buf[dt_offset], buf[dt_offset + 1]]);
        Ok(Self {
            dim,
            pixdim,
            datatype,
        })
    }
}

// ── MRI phantom ───────────────────────────────────────────────────────────────

/// Geometry type for an MRI phantom.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhantomGeometry {
    /// Spherical phantom with the given radius (mm).
    Sphere(f64),
    /// Cylindrical phantom with the given radius and height (mm).
    Cylinder {
        /// Radius of the cylinder (mm).
        radius: f64,
        /// Height of the cylinder (mm).
        height: f64,
    },
}

/// A simulated MRI phantom with known relaxation times.
///
/// Can simulate the spin-echo MRI signal for given TE/TR parameters.
#[derive(Debug, Clone)]
pub struct MriPhantom {
    /// Phantom geometry.
    pub geometry: PhantomGeometry,
    /// Longitudinal relaxation time T1 (ms).
    pub t1: f64,
    /// Transverse relaxation time T2 (ms).
    pub t2: f64,
    /// Proton density (a.u., typically normalised to water = 1.0).
    pub proton_density: f64,
}

impl MriPhantom {
    /// Create a new MRI phantom.
    pub fn new(geometry: PhantomGeometry, t1: f64, t2: f64, proton_density: f64) -> Self {
        Self {
            geometry,
            t1,
            t2,
            proton_density,
        }
    }

    /// Simulate the spin-echo MRI signal for echo time `te` and repetition
    /// time `tr` (both in ms).
    ///
    /// Uses the standard SE formula:
    /// `S = rho * (1 - exp(-TR/T1)) * exp(-TE/T2)`.
    pub fn simulate_signal(&self, te: f64, tr: f64) -> f64 {
        mri_signal_se(self.proton_density, self.t1, self.t2, tr, te)
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Map a Hounsfield Unit value to a material label string.
///
/// Thresholds are approximate clinical conventions:
/// - < −950 HU → `"air"`
/// - −950 to −100 HU → `"lung"`
/// - −100 to 20 HU → `"fat/soft_tissue"`
/// - 20 to 400 HU → `"soft_tissue/blood"`
/// - 400 to 1000 HU → `"bone"`
/// - ≥ 1000 HU → `"dense_bone/metal"`
pub fn hounsfield_to_material(hu: f64) -> &'static str {
    if hu < -950.0 {
        "air"
    } else if hu < -100.0 {
        "lung"
    } else if hu < 20.0 {
        "fat/soft_tissue"
    } else if hu < 400.0 {
        "soft_tissue/blood"
    } else if hu < 1000.0 {
        "bone"
    } else {
        "dense_bone/metal"
    }
}

/// Compute the spin-echo MRI signal.
///
/// `S = rho * (1 - exp(-TR / T1)) * exp(-TE / T2)`
///
/// All time parameters are in milliseconds.
pub fn mri_signal_se(rho: f64, t1: f64, t2: f64, tr: f64, te: f64) -> f64 {
    if t1 <= 0.0 || t2 <= 0.0 {
        return 0.0;
    }
    rho * (1.0 - (-tr / t1).exp()) * (-te / t2).exp()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── DicomTag ─────────────────────────────────────────────────────────

    #[test]
    fn test_dicom_tag_from_str() {
        let t = DicomTag::from_str(0x0008, 0x0060, "CS", "CT");
        assert_eq!(t.group, 0x0008);
        assert_eq!(t.element, 0x0060);
        assert_eq!(t.value_representation, "CS");
        assert_eq!(&t.data, b"CT");
    }

    #[test]
    fn test_dicom_tag_from_f64() {
        let t = DicomTag::from_f64(0x0028, 0x0030, "DS", 1.5);
        assert_eq!(t.data.len(), 8);
        let v = f64::from_le_bytes(t.data[..8].try_into().unwrap());
        assert!((v - 1.5).abs() < EPS);
    }

    // ── DicomDataset ─────────────────────────────────────────────────────

    #[test]
    fn test_dataset_get_string() {
        let mut ds = DicomDataset::new();
        ds.insert(DicomTag::from_str(0x0010, 0x0010, "PN", "Smith^John"));
        let name = ds.get_string(0x0010, 0x0010);
        assert_eq!(name, Some("Smith^John".to_string()));
    }

    #[test]
    fn test_dataset_get_string_missing() {
        let ds = DicomDataset::new();
        assert!(ds.get_string(0x0001, 0x0001).is_none());
    }

    #[test]
    fn test_dataset_get_f64_binary() {
        let mut ds = DicomDataset::new();
        ds.insert(DicomTag::from_f64(0x0028, 0x1053, "FD", 3.125));
        let v = ds.get_f64(0x0028, 0x1053).unwrap();
        assert!((v - 3.125).abs() < EPS);
    }

    #[test]
    fn test_dataset_get_f64_string() {
        let mut ds = DicomDataset::new();
        ds.insert(DicomTag::from_str(0x0028, 0x1052, "DS", "  42.5  "));
        let v = ds.get_f64(0x0028, 0x1052).unwrap();
        assert!((v - 42.5).abs() < EPS);
    }

    #[test]
    fn test_dataset_insert_overwrites() {
        let mut ds = DicomDataset::new();
        ds.insert(DicomTag::from_str(0x0010, 0x0010, "PN", "Old"));
        ds.insert(DicomTag::from_str(0x0010, 0x0010, "PN", "New"));
        assert_eq!(ds.get_string(0x0010, 0x0010), Some("New".to_string()));
    }

    // ── VoxelVolume ───────────────────────────────────────────────────────

    #[test]
    fn test_voxel_volume_count() {
        let v = VoxelVolume::new([4, 5, 6], [1.0; 3]);
        assert_eq!(v.voxel_count(), 120);
    }

    #[test]
    fn test_voxel_volume_physical_volume() {
        let v = VoxelVolume::new([10, 10, 10], [2.0, 2.0, 2.0]);
        assert!((v.physical_volume_mm3() - 8000.0).abs() < EPS);
    }

    #[test]
    fn test_to_hounsfield_water() {
        // CT water: pixel=0, slope=1, intercept=0 => HU=0
        let hu = VoxelVolume::to_hounsfield(0, 1.0, 0.0);
        assert!((hu).abs() < EPS);
    }

    #[test]
    fn test_to_hounsfield_bone() {
        let hu = VoxelVolume::to_hounsfield(700, 1.0, -1024.0);
        assert!((hu + 324.0).abs() < EPS);
    }

    #[test]
    fn test_voxel_get_in_bounds() {
        let v = VoxelVolume::new([3, 3, 3], [1.0; 3]);
        assert_eq!(v.get(0, 0, 0), Some(0));
    }

    #[test]
    fn test_voxel_get_out_of_bounds() {
        let v = VoxelVolume::new([3, 3, 3], [1.0; 3]);
        assert!(v.get(10, 0, 0).is_none());
    }

    // ── hounsfield_to_material ────────────────────────────────────────────

    #[test]
    fn test_hu_material_air() {
        assert_eq!(hounsfield_to_material(-1000.0), "air");
    }

    #[test]
    fn test_hu_material_lung() {
        assert_eq!(hounsfield_to_material(-500.0), "lung");
    }

    #[test]
    fn test_hu_material_fat() {
        assert_eq!(hounsfield_to_material(-50.0), "fat/soft_tissue");
    }

    #[test]
    fn test_hu_material_soft_tissue() {
        assert_eq!(hounsfield_to_material(50.0), "soft_tissue/blood");
    }

    #[test]
    fn test_hu_material_bone() {
        assert_eq!(hounsfield_to_material(700.0), "bone");
    }

    #[test]
    fn test_hu_material_dense_bone() {
        assert_eq!(hounsfield_to_material(1500.0), "dense_bone/metal");
    }

    // ── Segmentation ──────────────────────────────────────────────────────

    #[test]
    fn test_segmentation_volume_zero() {
        let seg = Segmentation::new(100, 3);
        let vol = seg.compute_volume(1, [1.0; 3]);
        assert!((vol).abs() < EPS);
    }

    #[test]
    fn test_segmentation_volume_all_labelled() {
        let mut seg = Segmentation::new(8, 1);
        seg.labels = vec![1; 8];
        let vol = seg.compute_volume(1, [2.0, 2.0, 2.0]);
        assert!((vol - 64.0).abs() < EPS);
    }

    #[test]
    fn test_segmentation_label_fraction() {
        let mut seg = Segmentation::new(10, 2);
        seg.labels[0] = 1;
        seg.labels[1] = 1;
        assert!((seg.label_fraction(1) - 0.2).abs() < EPS);
    }

    #[test]
    fn test_segmentation_empty() {
        let seg = Segmentation::new(0, 1);
        assert!((seg.label_fraction(1)).abs() < EPS);
    }

    // ── mri_signal_se ─────────────────────────────────────────────────────

    #[test]
    fn test_mri_signal_long_tr_short_te() {
        // Very long TR => (1 - exp(-TR/T1)) ≈ 1; very short TE => exp(-TE/T2) ≈ 1
        let s = mri_signal_se(1.0, 500.0, 100.0, 1e9, 0.0);
        assert!((s - 1.0).abs() < 1e-6, "signal should be ~rho: {s}");
    }

    #[test]
    fn test_mri_signal_zero_rho() {
        assert!((mri_signal_se(0.0, 500.0, 100.0, 1000.0, 10.0)).abs() < EPS);
    }

    #[test]
    fn test_mri_signal_invalid_t1() {
        assert!((mri_signal_se(1.0, 0.0, 100.0, 1000.0, 10.0)).abs() < EPS);
    }

    #[test]
    fn test_mri_signal_invalid_t2() {
        assert!((mri_signal_se(1.0, 500.0, 0.0, 1000.0, 10.0)).abs() < EPS);
    }

    #[test]
    fn test_mri_signal_t1_weighting() {
        // Higher T1 → lower signal for same TR
        let s_low = mri_signal_se(1.0, 300.0, 100.0, 600.0, 10.0);
        let s_high = mri_signal_se(1.0, 1500.0, 100.0, 600.0, 10.0);
        assert!(
            s_low > s_high,
            "lower T1 should give higher T1-weighted signal"
        );
    }

    #[test]
    fn test_mri_signal_t2_weighting() {
        // Longer TE → lower signal
        let s_short = mri_signal_se(1.0, 500.0, 80.0, 2000.0, 10.0);
        let s_long = mri_signal_se(1.0, 500.0, 80.0, 2000.0, 100.0);
        assert!(s_short > s_long, "short TE should give higher signal");
    }

    // ── MriPhantom ────────────────────────────────────────────────────────

    #[test]
    fn test_phantom_simulate_signal() {
        let p = MriPhantom::new(PhantomGeometry::Sphere(50.0), 800.0, 80.0, 1.0);
        let s = p.simulate_signal(10.0, 2000.0);
        assert!(s > 0.0 && s <= 1.0);
    }

    #[test]
    fn test_phantom_sphere_geometry() {
        let p = MriPhantom::new(PhantomGeometry::Sphere(25.0), 500.0, 60.0, 0.8);
        if let PhantomGeometry::Sphere(r) = p.geometry {
            assert!((r - 25.0).abs() < EPS);
        } else {
            panic!("expected sphere");
        }
    }

    #[test]
    fn test_phantom_cylinder_geometry() {
        let p = MriPhantom::new(
            PhantomGeometry::Cylinder {
                radius: 30.0,
                height: 100.0,
            },
            1000.0,
            100.0,
            1.0,
        );
        if let PhantomGeometry::Cylinder { radius, height } = p.geometry {
            assert!((radius - 30.0).abs() < EPS);
            assert!((height - 100.0).abs() < EPS);
        } else {
            panic!("expected cylinder");
        }
    }

    // ── NiftiHeader I/O ───────────────────────────────────────────────────

    #[test]
    fn test_nifti_roundtrip() {
        let path = std::env::temp_dir().join("test_nifti_header.bin");
        let hdr = NiftiHeader::new_3d(64, 128, 32, 0.5, 0.5, 1.0);
        hdr.write_header(path.to_str().unwrap_or("")).unwrap();
        let loaded = NiftiHeader::read_header(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(loaded.dim[0], 3);
        assert_eq!(loaded.dim[1], 64);
        assert_eq!(loaded.dim[2], 128);
        assert_eq!(loaded.dim[3], 32);
        assert!((loaded.pixdim[1] - 0.5).abs() < EPS);
        assert!((loaded.pixdim[3] - 1.0).abs() < EPS);
        assert_eq!(loaded.datatype, 4);
    }

    #[test]
    fn test_nifti_write_nonexistent_dir_fails() {
        let path = std::env::temp_dir()
            .join("nonexistent_dir_xyz")
            .join("header.bin");
        let hdr = NiftiHeader::new_3d(10, 10, 10, 1.0, 1.0, 1.0);
        assert!(hdr.write_header(path.to_str().unwrap_or("")).is_err());
    }

    #[test]
    fn test_nifti_read_nonexistent_fails() {
        let path = std::env::temp_dir().join("does_not_exist_nifti.bin");
        assert!(NiftiHeader::read_header(path.to_str().unwrap_or("")).is_err());
    }

    #[test]
    fn test_nifti_multiple_roundtrips() {
        for i in 0..3_u8 {
            let path = std::env::temp_dir().join(format!("test_nifti_{i}.bin"));
            let hdr = NiftiHeader::new_3d(10 + i as usize * 5, 20, 30, 1.0 + i as f64, 1.0, 1.0);
            hdr.write_header(path.to_str().unwrap_or("")).unwrap();
            let loaded = NiftiHeader::read_header(path.to_str().unwrap_or("")).unwrap();
            assert_eq!(loaded.dim[1], 10 + i as usize * 5);
            assert!((loaded.pixdim[1] - (1.0 + i as f64)).abs() < EPS);
        }
    }
}
