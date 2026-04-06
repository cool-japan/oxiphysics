//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use std::collections::VecDeque;

/// A DICOM tag identified by (group, element) pair.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DicomTag {
    /// Group number (e.g. 0x0010 for Patient).
    pub group: u16,
    /// Element number within the group.
    pub element: u16,
    /// Value Representation.
    pub vr: ValueRepresentation,
}
impl DicomTag {
    /// Creates a new DICOM tag.
    pub fn new(group: u16, element: u16, vr: ValueRepresentation) -> Self {
        Self { group, element, vr }
    }
    /// Returns true if this is a group-length tag (element == 0).
    pub fn is_group_length(&self) -> bool {
        self.element == 0
    }
    /// Returns true if this is a private tag (odd group number).
    pub fn is_private(&self) -> bool {
        self.group % 2 == 1
    }
    /// Formats the tag as `(GGGG,EEEE)`.
    pub fn tag_string(&self) -> String {
        format!("({:04X},{:04X})", self.group, self.element)
    }
}
/// A DICOM data element containing a tag and its value.
#[derive(Clone, Debug)]
pub struct DicomElement {
    /// The tag for this element.
    pub tag: DicomTag,
    /// The raw value bytes.
    pub value: Vec<u8>,
}
impl DicomElement {
    /// Creates a new DICOM element.
    pub fn new(tag: DicomTag, value: Vec<u8>) -> Self {
        Self { tag, value }
    }
    /// Interprets the value as a UTF-8 string, trimming trailing nulls/spaces.
    pub fn value_as_string(&self) -> Option<String> {
        let s = std::str::from_utf8(&self.value).ok()?;
        Some(s.trim_end_matches(['\0', ' ']).to_string())
    }
    /// Interprets the value as a little-endian u16.
    pub fn value_as_u16(&self) -> Option<u16> {
        if self.value.len() >= 2 {
            Some(u16::from_le_bytes([self.value[0], self.value[1]]))
        } else {
            None
        }
    }
    /// Interprets the value as a little-endian u32.
    pub fn value_as_u32(&self) -> Option<u32> {
        if self.value.len() >= 4 {
            Some(u32::from_le_bytes([
                self.value[0],
                self.value[1],
                self.value[2],
                self.value[3],
            ]))
        } else {
            None
        }
    }
    /// Interprets the value as a little-endian f64.
    pub fn value_as_f64(&self) -> Option<f64> {
        if self.value.len() >= 8 {
            let bytes: [u8; 8] = self.value[..8].try_into().ok()?;
            Some(f64::from_le_bytes(bytes))
        } else {
            None
        }
    }
    /// Interprets the value as a decimal string (DS VR) and returns f64.
    pub fn value_as_ds(&self) -> Option<f64> {
        let s = self.value_as_string()?;
        s.trim().parse::<f64>().ok()
    }
    /// Interprets the value as an integer string (IS VR) and returns i64.
    pub fn value_as_is(&self) -> Option<i64> {
        let s = self.value_as_string()?;
        s.trim().parse::<i64>().ok()
    }
}
/// Value Representation (VR) types used in DICOM data elements.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ValueRepresentation {
    /// Application Entity
    AE,
    /// Age String
    AS,
    /// Attribute Tag
    AT,
    /// Code String
    CS,
    /// Date
    DA,
    /// Decimal String
    DS,
    /// Date Time
    DT,
    /// Floating Point Double
    FD,
    /// Floating Point Single
    FL,
    /// Integer String
    IS,
    /// Long String
    LO,
    /// Long Text
    LT,
    /// Other Byte
    OB,
    /// Other Word
    OW,
    /// Person Name
    PN,
    /// Short String
    SH,
    /// Signed Long
    SL,
    /// Sequence of Items
    SQ,
    /// Signed Short
    SS,
    /// Short Text
    ST,
    /// Time
    TM,
    /// Unique Identifier (UID)
    UI,
    /// Unsigned Long
    UL,
    /// Unknown
    UN,
    /// Unsigned Short
    US,
    /// Unlimited Text
    UT,
}
impl ValueRepresentation {
    /// Parses a two-character VR code.
    pub fn from_bytes(b0: u8, b1: u8) -> Option<Self> {
        match (b0, b1) {
            (b'A', b'E') => Some(Self::AE),
            (b'A', b'S') => Some(Self::AS),
            (b'A', b'T') => Some(Self::AT),
            (b'C', b'S') => Some(Self::CS),
            (b'D', b'A') => Some(Self::DA),
            (b'D', b'S') => Some(Self::DS),
            (b'D', b'T') => Some(Self::DT),
            (b'F', b'D') => Some(Self::FD),
            (b'F', b'L') => Some(Self::FL),
            (b'I', b'S') => Some(Self::IS),
            (b'L', b'O') => Some(Self::LO),
            (b'L', b'T') => Some(Self::LT),
            (b'O', b'B') => Some(Self::OB),
            (b'O', b'W') => Some(Self::OW),
            (b'P', b'N') => Some(Self::PN),
            (b'S', b'H') => Some(Self::SH),
            (b'S', b'L') => Some(Self::SL),
            (b'S', b'Q') => Some(Self::SQ),
            (b'S', b'S') => Some(Self::SS),
            (b'S', b'T') => Some(Self::ST),
            (b'T', b'M') => Some(Self::TM),
            (b'U', b'I') => Some(Self::UI),
            (b'U', b'L') => Some(Self::UL),
            (b'U', b'N') => Some(Self::UN),
            (b'U', b'S') => Some(Self::US),
            (b'U', b'T') => Some(Self::UT),
            _ => None,
        }
    }
    /// Returns true if this VR uses a 4-byte length field (explicit VR).
    pub fn has_extended_length(&self) -> bool {
        matches!(
            self,
            Self::OB | Self::OW | Self::SQ | Self::UN | Self::UT | Self::LT
        )
    }
    /// Returns the two-character VR code as a string.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::AE => "AE",
            Self::AS => "AS",
            Self::AT => "AT",
            Self::CS => "CS",
            Self::DA => "DA",
            Self::DS => "DS",
            Self::DT => "DT",
            Self::FD => "FD",
            Self::FL => "FL",
            Self::IS => "IS",
            Self::LO => "LO",
            Self::LT => "LT",
            Self::OB => "OB",
            Self::OW => "OW",
            Self::PN => "PN",
            Self::SH => "SH",
            Self::SL => "SL",
            Self::SQ => "SQ",
            Self::SS => "SS",
            Self::ST => "ST",
            Self::TM => "TM",
            Self::UI => "UI",
            Self::UL => "UL",
            Self::UN => "UN",
            Self::US => "US",
            Self::UT => "UT",
        }
    }
}
/// NIfTI file reader that parses the 348-byte header and voxel data.
pub struct NiftiReader {
    /// The parsed header.
    pub header: NiftiHeader,
    /// The parsed voxel data as f32.
    pub voxels: Vec<f32>,
}
impl NiftiReader {
    /// Creates a new empty reader.
    pub fn new() -> Self {
        Self {
            header: NiftiHeader::new(),
            voxels: Vec::new(),
        }
    }
    /// Parses a NIfTI-1 (.nii) byte buffer.
    ///
    /// Reads the 348-byte header and interprets the voxel data based on
    /// the datatype field (supports float32=16, int16=4, uint8=2, float64=64).
    pub fn parse(&mut self, data: &[u8]) -> Result<(), String> {
        if data.len() < 348 {
            return Err("Data too short for NIfTI header".into());
        }
        let sizeof_hdr = read_i32_le(data, 0);
        if sizeof_hdr != 348 {
            return Err(format!("Invalid sizeof_hdr: {} (expected 348)", sizeof_hdr));
        }
        self.header.sizeof_hdr = sizeof_hdr;
        for i in 0..8 {
            self.header.dim[i] = read_i16_le(data, 40 + i * 2);
        }
        self.header.datatype = read_i16_le(data, 70);
        self.header.bitpix = read_i16_le(data, 72);
        for i in 0..8 {
            self.header.pixdim[i] = read_f32_le(data, 76 + i * 4);
        }
        self.header.vox_offset = read_f32_le(data, 108);
        self.header.scl_slope = read_f32_le(data, 112);
        self.header.scl_inter = read_f32_le(data, 116);
        self.header.qform_code = read_i16_le(data, 252);
        self.header.sform_code = read_i16_le(data, 254);
        self.header.quatern_b = read_f32_le(data, 256);
        self.header.quatern_c = read_f32_le(data, 260);
        self.header.quatern_d = read_f32_le(data, 264);
        self.header.qoffset_x = read_f32_le(data, 268);
        self.header.qoffset_y = read_f32_le(data, 272);
        self.header.qoffset_z = read_f32_le(data, 276);
        for i in 0..4 {
            self.header.srow_x[i] = read_f32_le(data, 280 + i * 4);
            self.header.srow_y[i] = read_f32_le(data, 296 + i * 4);
            self.header.srow_z[i] = read_f32_le(data, 312 + i * 4);
        }
        if data.len() >= 348 {
            let magic_bytes = &data[344..347.min(data.len())];
            self.header.magic = String::from_utf8_lossy(magic_bytes).to_string();
        }
        let offset = self.header.vox_offset as usize;
        if offset < data.len() {
            let voxel_data = &data[offset..];
            self.voxels =
                parse_nifti_voxels(voxel_data, self.header.datatype, self.header.num_voxels());
        }
        Ok(())
    }
    /// Converts parsed data into a `VoxelVolume`.
    pub fn to_volume(&self) -> VoxelVolume {
        let dims = self.header.volume_dims();
        let spacing = self.header.voxel_spacing();
        let origin = [
            self.header.qoffset_x as f64,
            self.header.qoffset_y as f64,
            self.header.qoffset_z as f64,
        ];
        let mut vol = VoxelVolume::new(dims, spacing, origin);
        let n = vol.num_voxels().min(self.voxels.len());
        vol.data[..n].copy_from_slice(&self.voxels[..n]);
        vol
    }
}
/// Volume processing operations for medical imaging data.
pub struct VolumeProcessing;
impl VolumeProcessing {
    /// Trilinear interpolation at fractional voxel coordinates.
    ///
    /// Returns `None` if the sample point is outside the volume bounds.
    pub fn trilinear_interpolate(volume: &VoxelVolume, x: f64, y: f64, z: f64) -> Option<f32> {
        let dims = volume.dims;
        if x < 0.0 || y < 0.0 || z < 0.0 {
            return None;
        }
        if x >= (dims[0] - 1) as f64 || y >= (dims[1] - 1) as f64 || z >= (dims[2] - 1) as f64 {
            return None;
        }
        let x0 = x.floor() as usize;
        let y0 = y.floor() as usize;
        let z0 = z.floor() as usize;
        let x1 = x0 + 1;
        let y1 = y0 + 1;
        let z1 = z0 + 1;
        let xd = (x - x0 as f64) as f32;
        let yd = (y - y0 as f64) as f32;
        let zd = (z - z0 as f64) as f32;
        let c000 = volume.get(x0, y0, z0)?;
        let c100 = volume.get(x1, y0, z0)?;
        let c010 = volume.get(x0, y1, z0)?;
        let c110 = volume.get(x1, y1, z0)?;
        let c001 = volume.get(x0, y0, z1)?;
        let c101 = volume.get(x1, y0, z1)?;
        let c011 = volume.get(x0, y1, z1)?;
        let c111 = volume.get(x1, y1, z1)?;
        let c00 = c000 * (1.0 - xd) + c100 * xd;
        let c01 = c001 * (1.0 - xd) + c101 * xd;
        let c10 = c010 * (1.0 - xd) + c110 * xd;
        let c11 = c011 * (1.0 - xd) + c111 * xd;
        let c0 = c00 * (1.0 - yd) + c10 * yd;
        let c1 = c01 * (1.0 - yd) + c11 * yd;
        Some(c0 * (1.0 - zd) + c1 * zd)
    }
    /// Applies 3D Gaussian smoothing to a volume.
    ///
    /// Uses separable convolution with the given sigma (in voxels) and
    /// kernel half-width `radius`.
    pub fn gaussian_smooth(volume: &VoxelVolume, sigma: f64, radius: usize) -> VoxelVolume {
        let dims = volume.dims;
        let n = dims[0] * dims[1] * dims[2];
        let kernel = build_gaussian_kernel(sigma, radius);
        let klen = kernel.len();
        let khalf = radius;
        let mut temp1 = vec![0.0f32; n];
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let mut sum = 0.0f32;
                    let mut wsum = 0.0f32;
                    for ki in 0..klen {
                        let sx = x as isize + ki as isize - khalf as isize;
                        if sx >= 0 && (sx as usize) < dims[0] {
                            let idx = coords_to_index(sx as usize, y, z, dims);
                            sum += volume.data[idx] * kernel[ki];
                            wsum += kernel[ki];
                        }
                    }
                    let idx = coords_to_index(x, y, z, dims);
                    temp1[idx] = if wsum > 0.0 { sum / wsum } else { 0.0 };
                }
            }
        }
        let mut temp2 = vec![0.0f32; n];
        for z in 0..dims[2] {
            for x in 0..dims[0] {
                for y in 0..dims[1] {
                    let mut sum = 0.0f32;
                    let mut wsum = 0.0f32;
                    for ki in 0..klen {
                        let sy = y as isize + ki as isize - khalf as isize;
                        if sy >= 0 && (sy as usize) < dims[1] {
                            let idx = coords_to_index(x, sy as usize, z, dims);
                            sum += temp1[idx] * kernel[ki];
                            wsum += kernel[ki];
                        }
                    }
                    let idx = coords_to_index(x, y, z, dims);
                    temp2[idx] = if wsum > 0.0 { sum / wsum } else { 0.0 };
                }
            }
        }
        let mut result = vec![0.0f32; n];
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                for z in 0..dims[2] {
                    let mut sum = 0.0f32;
                    let mut wsum = 0.0f32;
                    for ki in 0..klen {
                        let sz = z as isize + ki as isize - khalf as isize;
                        if sz >= 0 && (sz as usize) < dims[2] {
                            let idx = coords_to_index(x, y, sz as usize, dims);
                            sum += temp2[idx] * kernel[ki];
                            wsum += kernel[ki];
                        }
                    }
                    let idx = coords_to_index(x, y, z, dims);
                    result[idx] = if wsum > 0.0 { sum / wsum } else { 0.0 };
                }
            }
        }
        VoxelVolume {
            data: result,
            dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Performs histogram equalization on a volume.
    ///
    /// Maps voxel intensities to the range \[0, `max_output`\] using a
    /// cumulative histogram with `num_bins` bins.
    pub fn histogram_equalize(
        volume: &VoxelVolume,
        num_bins: usize,
        max_output: f32,
    ) -> VoxelVolume {
        let n = volume.data.len();
        if n == 0 {
            return volume.clone();
        }
        let vmin = volume.min_value();
        let vmax = volume.max_value();
        let range = vmax - vmin;
        if range <= f32::EPSILON {
            return volume.clone();
        }
        let mut histogram = vec![0u64; num_bins];
        for &v in &volume.data {
            let bin = ((v - vmin) / range * (num_bins - 1) as f32)
                .round()
                .clamp(0.0, (num_bins - 1) as f32) as usize;
            histogram[bin] += 1;
        }
        let mut cdf = vec![0u64; num_bins];
        cdf[0] = histogram[0];
        for i in 1..num_bins {
            cdf[i] = cdf[i - 1] + histogram[i];
        }
        let cdf_min = cdf.iter().copied().find(|&c| c > 0).unwrap_or(0);
        let denom = (n as u64).saturating_sub(cdf_min);
        let mut data = vec![0.0f32; n];
        if denom == 0 {
            return volume.clone();
        }
        for (i, &v) in volume.data.iter().enumerate() {
            let bin = ((v - vmin) / range * (num_bins - 1) as f32)
                .round()
                .clamp(0.0, (num_bins - 1) as f32) as usize;
            data[i] = ((cdf[bin] - cdf_min) as f64 / denom as f64 * max_output as f64) as f32;
        }
        VoxelVolume {
            data,
            dims: volume.dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Converts raw CT pixel values to Hounsfield Units (HU).
    ///
    /// HU = pixel_value * slope + intercept
    ///
    /// Standard HU values: air = -1000, water = 0, bone > 400.
    pub fn to_hounsfield_units(volume: &VoxelVolume, slope: f64, intercept: f64) -> VoxelVolume {
        let data: Vec<f32> = volume
            .data
            .iter()
            .map(|&v| (v as f64 * slope + intercept) as f32)
            .collect();
        VoxelVolume {
            data,
            dims: volume.dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Clamps voxel values to a Hounsfield Unit window.
    ///
    /// Values outside \[center - width/2, center + width/2\] are clamped
    /// and then rescaled to \[0, 1\].
    pub fn window_level(volume: &VoxelVolume, center: f64, width: f64) -> VoxelVolume {
        let low = (center - width / 2.0) as f32;
        let high = (center + width / 2.0) as f32;
        let range = high - low;
        let data: Vec<f32> = volume
            .data
            .iter()
            .map(|&v| {
                if range <= f32::EPSILON {
                    0.0
                } else {
                    ((v - low) / range).clamp(0.0, 1.0)
                }
            })
            .collect();
        VoxelVolume {
            data,
            dims: volume.dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Resamples a volume to a new set of dimensions using trilinear interpolation.
    pub fn resample(volume: &VoxelVolume, new_dims: [usize; 3]) -> VoxelVolume {
        let new_n = new_dims[0] * new_dims[1] * new_dims[2];
        let mut data = vec![0.0f32; new_n];
        let scale_x = if new_dims[0] > 1 {
            (volume.dims[0] - 1) as f64 / (new_dims[0] - 1) as f64
        } else {
            0.0
        };
        let scale_y = if new_dims[1] > 1 {
            (volume.dims[1] - 1) as f64 / (new_dims[1] - 1) as f64
        } else {
            0.0
        };
        let scale_z = if new_dims[2] > 1 {
            (volume.dims[2] - 1) as f64 / (new_dims[2] - 1) as f64
        } else {
            0.0
        };
        let new_spacing = [
            volume.spacing[0] * volume.dims[0] as f64 / new_dims[0].max(1) as f64,
            volume.spacing[1] * volume.dims[1] as f64 / new_dims[1].max(1) as f64,
            volume.spacing[2] * volume.dims[2] as f64 / new_dims[2].max(1) as f64,
        ];
        for nz in 0..new_dims[2] {
            for ny in 0..new_dims[1] {
                for nx in 0..new_dims[0] {
                    let sx = nx as f64 * scale_x;
                    let sy = ny as f64 * scale_y;
                    let sz = nz as f64 * scale_z;
                    let val =
                        Self::trilinear_interpolate(volume, sx, sy, sz).unwrap_or_else(|| {
                            let ix = sx.round().clamp(0.0, (volume.dims[0] - 1) as f64) as usize;
                            let iy = sy.round().clamp(0.0, (volume.dims[1] - 1) as f64) as usize;
                            let iz = sz.round().clamp(0.0, (volume.dims[2] - 1) as f64) as usize;
                            volume.get(ix, iy, iz).unwrap_or(0.0)
                        });
                    let idx = nz * new_dims[0] * new_dims[1] + ny * new_dims[0] + nx;
                    data[idx] = val;
                }
            }
        }
        VoxelVolume {
            data,
            dims: new_dims,
            spacing: new_spacing,
            origin: volume.origin,
        }
    }
    /// Computes a 3D gradient magnitude volume using central differences.
    pub fn gradient_magnitude(volume: &VoxelVolume) -> VoxelVolume {
        let dims = volume.dims;
        let n = dims[0] * dims[1] * dims[2];
        let mut data = vec![0.0f32; n];
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let gx = if x > 0 && x < dims[0] - 1 {
                        let v1 = volume.get(x + 1, y, z).unwrap_or(0.0);
                        let v0 = volume.get(x - 1, y, z).unwrap_or(0.0);
                        (v1 - v0) / (2.0 * volume.spacing[0] as f32)
                    } else {
                        0.0
                    };
                    let gy = if y > 0 && y < dims[1] - 1 {
                        let v1 = volume.get(x, y + 1, z).unwrap_or(0.0);
                        let v0 = volume.get(x, y - 1, z).unwrap_or(0.0);
                        (v1 - v0) / (2.0 * volume.spacing[1] as f32)
                    } else {
                        0.0
                    };
                    let gz = if z > 0 && z < dims[2] - 1 {
                        let v1 = volume.get(x, y, z + 1).unwrap_or(0.0);
                        let v0 = volume.get(x, y, z - 1).unwrap_or(0.0);
                        (v1 - v0) / (2.0 * volume.spacing[2] as f32)
                    } else {
                        0.0
                    };
                    let idx = coords_to_index(x, y, z, dims);
                    data[idx] = (gx * gx + gy * gy + gz * gz).sqrt();
                }
            }
        }
        VoxelVolume {
            data,
            dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Applies a binary threshold to a volume.
    ///
    /// Voxels >= threshold become 1.0, others become 0.0.
    pub fn threshold(volume: &VoxelVolume, threshold: f32) -> VoxelVolume {
        let data: Vec<f32> = volume
            .data
            .iter()
            .map(|&v| if v >= threshold { 1.0 } else { 0.0 })
            .collect();
        VoxelVolume {
            data,
            dims: volume.dims,
            spacing: volume.spacing,
            origin: volume.origin,
        }
    }
    /// Computes a histogram of voxel intensities.
    ///
    /// Returns `(bin_centers, counts)`.
    pub fn histogram(volume: &VoxelVolume, num_bins: usize) -> (Vec<f64>, Vec<u64>) {
        let vmin = volume.min_value() as f64;
        let vmax = volume.max_value() as f64;
        let range = vmax - vmin;
        if range <= f64::EPSILON || num_bins == 0 {
            return (vec![vmin], vec![volume.data.len() as u64]);
        }
        let bin_width = range / num_bins as f64;
        let centers: Vec<f64> = (0..num_bins)
            .map(|i| vmin + (i as f64 + 0.5) * bin_width)
            .collect();
        let mut counts = vec![0u64; num_bins];
        for &v in &volume.data {
            let bin = ((v as f64 - vmin) / range * (num_bins - 1) as f64)
                .round()
                .clamp(0.0, (num_bins - 1) as f64) as usize;
            counts[bin] += 1;
        }
        (centers, counts)
    }
}
/// A 3D voxel volume stored as a flat `Vec`f32`.
///
/// The data is stored in row-major order: index = z * (dimx * dimy) + y * dimx + x.
#[derive(Clone, Debug)]
pub struct VoxelVolume {
    /// Raw voxel intensity data, length = dims\[0\]*dims\[1\]*dims\[2\].
    pub data: Vec<f32>,
    /// Volume dimensions [x, y, z].
    pub dims: [usize; 3],
    /// Voxel spacing in mm [dx, dy, dz].
    pub spacing: [f64; 3],
    /// Volume origin in world coordinates [ox, oy, oz].
    pub origin: [f64; 3],
}
impl VoxelVolume {
    /// Creates a new `VoxelVolume` filled with zeros.
    ///
    /// # Arguments
    /// * `dims` - Dimensions [x, y, z]
    /// * `spacing` - Voxel spacing in mm
    /// * `origin` - World-space origin
    pub fn new(dims: [usize; 3], spacing: [f64; 3], origin: [f64; 3]) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        Self {
            data: vec![0.0; n],
            dims,
            spacing,
            origin,
        }
    }
    /// Creates a `VoxelVolume` from existing data.
    ///
    /// Returns `None` if `data.len() != dims\[0\]*dims\[1\]*dims\[2\]`.
    pub fn from_data(
        data: Vec<f32>,
        dims: [usize; 3],
        spacing: [f64; 3],
        origin: [f64; 3],
    ) -> Option<Self> {
        if data.len() != dims[0] * dims[1] * dims[2] {
            return None;
        }
        Some(Self {
            data,
            dims,
            spacing,
            origin,
        })
    }
    /// Returns the total number of voxels.
    pub fn num_voxels(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
    /// Converts (x, y, z) indices to a flat index.
    ///
    /// Returns `None` if out of bounds.
    pub fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        if x < self.dims[0] && y < self.dims[1] && z < self.dims[2] {
            Some(z * self.dims[0] * self.dims[1] + y * self.dims[0] + x)
        } else {
            None
        }
    }
    /// Gets the voxel value at (x, y, z).
    pub fn get(&self, x: usize, y: usize, z: usize) -> Option<f32> {
        self.index(x, y, z).map(|i| self.data[i])
    }
    /// Sets the voxel value at (x, y, z).
    pub fn set(&mut self, x: usize, y: usize, z: usize, value: f32) {
        if let Some(i) = self.index(x, y, z) {
            self.data[i] = value;
        }
    }
    /// Converts a flat index back to (x, y, z) coordinates.
    pub fn coords_from_index(&self, idx: usize) -> Option<(usize, usize, usize)> {
        if idx >= self.data.len() {
            return None;
        }
        let xy = self.dims[0] * self.dims[1];
        let z = idx / xy;
        let rem = idx % xy;
        let y = rem / self.dims[0];
        let x = rem % self.dims[0];
        Some((x, y, z))
    }
    /// Converts voxel coordinates to world coordinates using spacing and origin.
    pub fn voxel_to_world(&self, x: f64, y: f64, z: f64) -> [f64; 3] {
        [
            self.origin[0] + x * self.spacing[0],
            self.origin[1] + y * self.spacing[1],
            self.origin[2] + z * self.spacing[2],
        ]
    }
    /// Converts world coordinates to voxel coordinates.
    pub fn world_to_voxel(&self, wx: f64, wy: f64, wz: f64) -> [f64; 3] {
        [
            (wx - self.origin[0]) / self.spacing[0],
            (wy - self.origin[1]) / self.spacing[1],
            (wz - self.origin[2]) / self.spacing[2],
        ]
    }
    /// Returns the physical size of the volume in mm.
    pub fn physical_size(&self) -> [f64; 3] {
        [
            self.dims[0] as f64 * self.spacing[0],
            self.dims[1] as f64 * self.spacing[1],
            self.dims[2] as f64 * self.spacing[2],
        ]
    }
    /// Computes the minimum voxel value.
    pub fn min_value(&self) -> f32 {
        self.data.iter().copied().fold(f32::INFINITY, f32::min)
    }
    /// Computes the maximum voxel value.
    pub fn max_value(&self) -> f32 {
        self.data.iter().copied().fold(f32::NEG_INFINITY, f32::max)
    }
    /// Computes the mean voxel value.
    pub fn mean_value(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.data.iter().map(|&v| v as f64).sum();
        sum / self.data.len() as f64
    }
}
/// Parsed DICOM header information for a single image.
#[derive(Clone, Debug, Default)]
pub struct DicomHeader {
    /// Patient name (0010,0010).
    pub patient_name: String,
    /// Patient ID (0010,0020).
    pub patient_id: String,
    /// Study date (0008,0020) as YYYYMMDD string.
    pub study_date: String,
    /// Modality (0008,0060), e.g. "CT", "MR".
    pub modality: String,
    /// Slice thickness in mm (0018,0050).
    pub slice_thickness: f64,
    /// Pixel spacing [row, col] in mm (0028,0030).
    pub pixel_spacing: [f64; 2],
    /// Image rows (0028,0010).
    pub rows: u16,
    /// Image columns (0028,0011).
    pub columns: u16,
    /// Bits allocated (0028,0100).
    pub bits_allocated: u16,
    /// Bits stored (0028,0101).
    pub bits_stored: u16,
    /// Rescale intercept (0028,1052).
    pub rescale_intercept: f64,
    /// Rescale slope (0028,1053).
    pub rescale_slope: f64,
    /// Window center (0028,1050).
    pub window_center: f64,
    /// Window width (0028,1051).
    pub window_width: f64,
    /// Image position patient [x, y, z] (0020,0032).
    pub image_position: [f64; 3],
    /// Instance number (0020,0013).
    pub instance_number: i64,
}
impl DicomHeader {
    /// Creates a header with default values and rescale_slope = 1.
    pub fn new() -> Self {
        Self {
            rescale_slope: 1.0,
            ..Default::default()
        }
    }
}
/// Statistics for a single labeled region.
#[derive(Clone, Debug)]
pub struct LabelStats {
    /// The label id.
    pub label: u32,
    /// Number of voxels in this component.
    pub voxel_count: usize,
    /// Physical volume in mm^3.
    pub volume_mm3: f64,
    /// Centroid in voxel coordinates [x, y, z].
    pub centroid: [f64; 3],
    /// Bounding box minimum [x, y, z].
    pub bbox_min: [usize; 3],
    /// Bounding box maximum [x, y, z].
    pub bbox_max: [usize; 3],
    /// Mean intensity of the region (from the source volume).
    pub mean_intensity: f64,
}
/// NIfTI-1 header (348 bytes).
///
/// Stores dimensional information, voxel sizes, spatial transforms, and
/// data type metadata per the NIfTI-1 specification.
#[derive(Clone, Debug)]
pub struct NiftiHeader {
    /// Header size (must be 348 for NIfTI-1).
    pub sizeof_hdr: i32,
    /// Dimensions: dim\[0\] = ndim, dim\[1..ndim\] = sizes.
    pub dim: [i16; 8],
    /// Data type code (e.g., 16 = float32).
    pub datatype: i16,
    /// Bits per voxel.
    pub bitpix: i16,
    /// Voxel sizes: pixdim[1..ndim] = spacing in mm.
    pub pixdim: [f32; 8],
    /// Byte offset to voxel data.
    pub vox_offset: f32,
    /// Data scaling slope.
    pub scl_slope: f32,
    /// Data scaling intercept.
    pub scl_inter: f32,
    /// sform code (0=unknown, 1=scanner, 2=aligned, 3=talairach, 4=mni).
    pub sform_code: i16,
    /// qform code.
    pub qform_code: i16,
    /// sform affine row 0.
    pub srow_x: [f32; 4],
    /// sform affine row 1.
    pub srow_y: [f32; 4],
    /// sform affine row 2.
    pub srow_z: [f32; 4],
    /// Quaternion b parameter.
    pub quatern_b: f32,
    /// Quaternion c parameter.
    pub quatern_c: f32,
    /// Quaternion d parameter.
    pub quatern_d: f32,
    /// Quaternion x offset.
    pub qoffset_x: f32,
    /// Quaternion y offset.
    pub qoffset_y: f32,
    /// Quaternion z offset.
    pub qoffset_z: f32,
    /// Description string (up to 80 chars).
    pub descrip: String,
    /// Magic string ("n+1" for .nii, "ni1" for .hdr/.img pair).
    pub magic: String,
}
impl NiftiHeader {
    /// Creates a default NIfTI-1 header.
    pub fn new() -> Self {
        Self {
            sizeof_hdr: 348,
            dim: [3, 1, 1, 1, 0, 0, 0, 0],
            datatype: 16,
            bitpix: 32,
            pixdim: [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            vox_offset: 352.0,
            scl_slope: 1.0,
            scl_inter: 0.0,
            sform_code: 0,
            qform_code: 0,
            srow_x: [1.0, 0.0, 0.0, 0.0],
            srow_y: [0.0, 1.0, 0.0, 0.0],
            srow_z: [0.0, 0.0, 1.0, 0.0],
            quatern_b: 0.0,
            quatern_c: 0.0,
            quatern_d: 0.0,
            qoffset_x: 0.0,
            qoffset_y: 0.0,
            qoffset_z: 0.0,
            descrip: String::new(),
            magic: "n+1".to_string(),
        }
    }
    /// Returns the number of dimensions.
    pub fn ndim(&self) -> usize {
        self.dim[0].max(0) as usize
    }
    /// Returns the volume dimensions as (nx, ny, nz).
    ///
    /// For lower-dimensional data, missing dims default to 1.
    pub fn volume_dims(&self) -> [usize; 3] {
        let nd = self.ndim();
        [
            if nd >= 1 {
                self.dim[1].max(1) as usize
            } else {
                1
            },
            if nd >= 2 {
                self.dim[2].max(1) as usize
            } else {
                1
            },
            if nd >= 3 {
                self.dim[3].max(1) as usize
            } else {
                1
            },
        ]
    }
    /// Returns voxel spacing as [dx, dy, dz].
    pub fn voxel_spacing(&self) -> [f64; 3] {
        [
            self.pixdim[1] as f64,
            self.pixdim[2] as f64,
            self.pixdim[3] as f64,
        ]
    }
    /// Returns the total number of voxels.
    pub fn num_voxels(&self) -> usize {
        let d = self.volume_dims();
        d[0] * d[1] * d[2]
    }
    /// Returns the sform affine as a 4x4 matrix in row-major order.
    pub fn sform_affine(&self) -> [[f32; 4]; 4] {
        [self.srow_x, self.srow_y, self.srow_z, [0.0, 0.0, 0.0, 1.0]]
    }
    /// Computes the qform affine from quaternion parameters.
    pub fn qform_affine(&self) -> [[f64; 4]; 4] {
        let b = self.quatern_b as f64;
        let c = self.quatern_c as f64;
        let d = self.quatern_d as f64;
        let a_sq = 1.0 - b * b - c * c - d * d;
        let a = if a_sq > 0.0 { a_sq.sqrt() } else { 0.0 };
        let r11 = a * a + b * b - c * c - d * d;
        let r12 = 2.0 * (b * c - a * d);
        let r13 = 2.0 * (b * d + a * c);
        let r21 = 2.0 * (b * c + a * d);
        let r22 = a * a - b * b + c * c - d * d;
        let r23 = 2.0 * (c * d - a * b);
        let r31 = 2.0 * (b * d - a * c);
        let r32 = 2.0 * (c * d + a * b);
        let r33 = a * a - b * b - c * c + d * d;
        let dx = self.pixdim[1] as f64;
        let dy = self.pixdim[2] as f64;
        let dz = self.pixdim[3] as f64;
        let qfac = if self.pixdim[0] < 0.0 { -1.0 } else { 1.0 };
        [
            [r11 * dx, r12 * dy, r13 * dz * qfac, self.qoffset_x as f64],
            [r21 * dx, r22 * dy, r23 * dz * qfac, self.qoffset_y as f64],
            [r31 * dx, r32 * dy, r33 * dz * qfac, self.qoffset_z as f64],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}
/// DICOM file reader that parses the preamble, magic number, and data elements.
///
/// Supports Explicit VR Little Endian transfer syntax.
pub struct DicomReader {
    /// All parsed data elements.
    pub elements: Vec<DicomElement>,
    /// Pixel data bytes (7FE0,0010).
    pub pixel_data: Vec<u8>,
    /// Current read position in the byte stream.
    pub(super) _position: usize,
}
impl DicomReader {
    /// Creates a new empty `DicomReader`.
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
            pixel_data: Vec::new(),
            _position: 0,
        }
    }
    /// Parses a DICOM byte stream.
    ///
    /// The stream must start with a 128-byte preamble followed by "DICM".
    /// Data elements are read as Explicit VR Little Endian.
    pub fn parse(&mut self, data: &[u8]) -> Result<DicomHeader, String> {
        if data.len() < 132 {
            return Err("Data too short for DICOM preamble + magic".into());
        }
        if &data[128..132] != b"DICM" {
            return Err("Missing DICM magic number at offset 128".into());
        }
        self._position = 132;
        self.elements.clear();
        self.pixel_data.clear();
        let mut header = DicomHeader::new();
        while self._position + 8 <= data.len() {
            let group = read_u16_le(data, self._position);
            let element = read_u16_le(data, self._position + 2);
            if group == 0x7FE0 && element == 0x0010 {
                self._position += 4;
                if self._position + 2 > data.len() {
                    break;
                }
                let vr =
                    ValueRepresentation::from_bytes(data[self._position], data[self._position + 1]);
                self._position += 2;
                let length = if vr.is_none_or(|v| v.has_extended_length()) {
                    self._position += 2;
                    if self._position + 4 > data.len() {
                        break;
                    }
                    let len = read_u32_le(data, self._position) as usize;
                    self._position += 4;
                    len
                } else {
                    if self._position + 2 > data.len() {
                        break;
                    }
                    let len = read_u16_le(data, self._position) as usize;
                    self._position += 2;
                    len
                };
                let end = (self._position + length).min(data.len());
                self.pixel_data = data[self._position..end].to_vec();
                self._position = end;
                continue;
            }
            if self._position + 4 > data.len() {
                break;
            }
            let vr_opt =
                ValueRepresentation::from_bytes(data[self._position + 4], data[self._position + 5]);
            let vr = vr_opt.unwrap_or(ValueRepresentation::UN);
            self._position += 6;
            let length = if vr.has_extended_length() {
                self._position += 2;
                if self._position + 4 > data.len() {
                    break;
                }
                let len = read_u32_le(data, self._position) as usize;
                self._position += 4;
                len
            } else {
                if self._position + 2 > data.len() {
                    break;
                }
                let len = read_u16_le(data, self._position) as usize;
                self._position += 2;
                len
            };
            if length == 0xFFFFFFFF {
                break;
            }
            let end = (self._position + length).min(data.len());
            let value = data[self._position..end].to_vec();
            self._position = end;
            let tag = DicomTag::new(group, element, vr);
            let elem = DicomElement::new(tag, value);
            Self::populate_header(&elem, &mut header);
            self.elements.push(elem);
        }
        Ok(header)
    }
    /// Populates header fields from a parsed element.
    fn populate_header(elem: &DicomElement, header: &mut DicomHeader) {
        let g = elem.tag.group;
        let e = elem.tag.element;
        match (g, e) {
            (0x0010, 0x0010) => {
                if let Some(s) = elem.value_as_string() {
                    header.patient_name = s;
                }
            }
            (0x0010, 0x0020) => {
                if let Some(s) = elem.value_as_string() {
                    header.patient_id = s;
                }
            }
            (0x0008, 0x0020) => {
                if let Some(s) = elem.value_as_string() {
                    header.study_date = s;
                }
            }
            (0x0008, 0x0060) => {
                if let Some(s) = elem.value_as_string() {
                    header.modality = s;
                }
            }
            (0x0018, 0x0050) => {
                if let Some(v) = elem.value_as_ds() {
                    header.slice_thickness = v;
                }
            }
            (0x0028, 0x0030) => {
                if let Some(s) = elem.value_as_string() {
                    let parts: Vec<&str> = s.split('\\').collect();
                    if parts.len() >= 2
                        && let (Ok(r), Ok(c)) = (
                            parts[0].trim().parse::<f64>(),
                            parts[1].trim().parse::<f64>(),
                        )
                    {
                        header.pixel_spacing = [r, c];
                    }
                }
            }
            (0x0028, 0x0010) => {
                if let Some(v) = elem.value_as_u16() {
                    header.rows = v;
                }
            }
            (0x0028, 0x0011) => {
                if let Some(v) = elem.value_as_u16() {
                    header.columns = v;
                }
            }
            (0x0028, 0x0100) => {
                if let Some(v) = elem.value_as_u16() {
                    header.bits_allocated = v;
                }
            }
            (0x0028, 0x0101) => {
                if let Some(v) = elem.value_as_u16() {
                    header.bits_stored = v;
                }
            }
            (0x0028, 0x1052) => {
                if let Some(v) = elem.value_as_ds() {
                    header.rescale_intercept = v;
                }
            }
            (0x0028, 0x1053) => {
                if let Some(v) = elem.value_as_ds() {
                    header.rescale_slope = v;
                }
            }
            (0x0028, 0x1050) => {
                if let Some(v) = elem.value_as_ds() {
                    header.window_center = v;
                }
            }
            (0x0028, 0x1051) => {
                if let Some(v) = elem.value_as_ds() {
                    header.window_width = v;
                }
            }
            (0x0020, 0x0032) => {
                if let Some(s) = elem.value_as_string() {
                    let parts: Vec<&str> = s.split('\\').collect();
                    if parts.len() >= 3
                        && let (Ok(x), Ok(y), Ok(z)) = (
                            parts[0].trim().parse::<f64>(),
                            parts[1].trim().parse::<f64>(),
                            parts[2].trim().parse::<f64>(),
                        )
                    {
                        header.image_position = [x, y, z];
                    }
                }
            }
            (0x0020, 0x0013) => {
                if let Some(v) = elem.value_as_is() {
                    header.instance_number = v;
                }
            }
            _ => {}
        }
    }
    /// Finds an element by (group, element).
    pub fn find_element(&self, group: u16, element: u16) -> Option<&DicomElement> {
        self.elements
            .iter()
            .find(|e| e.tag.group == group && e.tag.element == element)
    }
    /// Returns the number of parsed data elements (excluding pixel data).
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
}
/// A binary segmentation mask with connected component labeling.
#[derive(Clone, Debug)]
pub struct Segmentation {
    /// Label map: 0 = background, >0 = component label.
    pub labels: Vec<u32>,
    /// Volume dimensions [x, y, z].
    pub dims: [usize; 3],
    /// Voxel spacing in mm.
    pub spacing: [f64; 3],
    /// Number of connected components found.
    pub num_components: u32,
}
impl Segmentation {
    /// Creates a segmentation from a binary mask volume using 3D flood fill.
    ///
    /// Any voxel with value >= `threshold` is considered foreground.
    /// Uses 6-connectivity (face-adjacent neighbors).
    pub fn from_threshold(volume: &VoxelVolume, threshold: f32) -> Self {
        let dims = volume.dims;
        let n = dims[0] * dims[1] * dims[2];
        let mut labels = vec![0u32; n];
        let mut current_label = 0u32;
        let foreground: Vec<bool> = volume.data.iter().map(|&v| v >= threshold).collect();
        for idx in 0..n {
            if !foreground[idx] || labels[idx] != 0 {
                continue;
            }
            current_label += 1;
            let mut queue = VecDeque::new();
            queue.push_back(idx);
            labels[idx] = current_label;
            while let Some(cidx) = queue.pop_front() {
                let (cx, cy, cz) = index_to_coords(cidx, dims);
                let neighbors: [(isize, isize, isize); 6] = [
                    (-1, 0, 0),
                    (1, 0, 0),
                    (0, -1, 0),
                    (0, 1, 0),
                    (0, 0, -1),
                    (0, 0, 1),
                ];
                for (dx, dy, dz) in &neighbors {
                    let nx = cx as isize + dx;
                    let ny = cy as isize + dy;
                    let nz = cz as isize + dz;
                    if nx < 0 || ny < 0 || nz < 0 {
                        continue;
                    }
                    let (nx, ny, nz) = (nx as usize, ny as usize, nz as usize);
                    if nx >= dims[0] || ny >= dims[1] || nz >= dims[2] {
                        continue;
                    }
                    let nidx = coords_to_index(nx, ny, nz, dims);
                    if foreground[nidx] && labels[nidx] == 0 {
                        labels[nidx] = current_label;
                        queue.push_back(nidx);
                    }
                }
            }
        }
        Self {
            labels,
            dims,
            spacing: volume.spacing,
            num_components: current_label,
        }
    }
    /// Creates a segmentation directly from a pre-computed label array.
    pub fn from_labels(labels: Vec<u32>, dims: [usize; 3], spacing: [f64; 3]) -> Option<Self> {
        if labels.len() != dims[0] * dims[1] * dims[2] {
            return None;
        }
        let max_label = labels.iter().copied().max().unwrap_or(0);
        Some(Self {
            labels,
            dims,
            spacing,
            num_components: max_label,
        })
    }
    /// Computes statistics for each labeled component.
    ///
    /// Uses the original volume for intensity statistics.
    pub fn label_statistics(&self, volume: &VoxelVolume) -> Vec<LabelStats> {
        if self.num_components == 0 {
            return Vec::new();
        }
        let voxel_vol = self.spacing[0] * self.spacing[1] * self.spacing[2];
        let nc = self.num_components as usize;
        let mut counts = vec![0usize; nc + 1];
        let mut sum_x = vec![0.0f64; nc + 1];
        let mut sum_y = vec![0.0f64; nc + 1];
        let mut sum_z = vec![0.0f64; nc + 1];
        let mut sum_intensity = vec![0.0f64; nc + 1];
        let mut bbox_min_x = vec![usize::MAX; nc + 1];
        let mut bbox_min_y = vec![usize::MAX; nc + 1];
        let mut bbox_min_z = vec![usize::MAX; nc + 1];
        let mut bbox_max_x = vec![0usize; nc + 1];
        let mut bbox_max_y = vec![0usize; nc + 1];
        let mut bbox_max_z = vec![0usize; nc + 1];
        for (idx, &lbl) in self.labels.iter().enumerate() {
            if lbl == 0 {
                continue;
            }
            let l = lbl as usize;
            let (x, y, z) = index_to_coords(idx, self.dims);
            counts[l] += 1;
            sum_x[l] += x as f64;
            sum_y[l] += y as f64;
            sum_z[l] += z as f64;
            if idx < volume.data.len() {
                sum_intensity[l] += volume.data[idx] as f64;
            }
            bbox_min_x[l] = bbox_min_x[l].min(x);
            bbox_min_y[l] = bbox_min_y[l].min(y);
            bbox_min_z[l] = bbox_min_z[l].min(z);
            bbox_max_x[l] = bbox_max_x[l].max(x);
            bbox_max_y[l] = bbox_max_y[l].max(y);
            bbox_max_z[l] = bbox_max_z[l].max(z);
        }
        let mut result = Vec::with_capacity(nc);
        for l in 1..=nc {
            if counts[l] == 0 {
                continue;
            }
            let c = counts[l] as f64;
            result.push(LabelStats {
                label: l as u32,
                voxel_count: counts[l],
                volume_mm3: counts[l] as f64 * voxel_vol,
                centroid: [sum_x[l] / c, sum_y[l] / c, sum_z[l] / c],
                bbox_min: [bbox_min_x[l], bbox_min_y[l], bbox_min_z[l]],
                bbox_max: [bbox_max_x[l], bbox_max_y[l], bbox_max_z[l]],
                mean_intensity: sum_intensity[l] / c,
            });
        }
        result
    }
    /// Returns the label at the given voxel coordinates.
    pub fn get_label(&self, x: usize, y: usize, z: usize) -> Option<u32> {
        if x < self.dims[0] && y < self.dims[1] && z < self.dims[2] {
            Some(self.labels[coords_to_index(x, y, z, self.dims)])
        } else {
            None
        }
    }
    /// Returns the voxel count for a given label.
    pub fn component_size(&self, label: u32) -> usize {
        self.labels.iter().filter(|&&l| l == label).count()
    }
    /// Extracts a binary mask for a single label as a `VoxelVolume`.
    pub fn extract_mask(&self, label: u32, spacing: [f64; 3]) -> VoxelVolume {
        let data: Vec<f32> = self
            .labels
            .iter()
            .map(|&l| if l == label { 1.0 } else { 0.0 })
            .collect();
        VoxelVolume {
            data,
            dims: self.dims,
            spacing,
            origin: [0.0, 0.0, 0.0],
        }
    }
}
