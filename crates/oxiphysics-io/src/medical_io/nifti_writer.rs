// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! NIfTI-1 single-file (`.nii` / `.nii.gz`) writer for [`NiftiReader`] volumes.
//!
//! The output is a complete NIfTI-1 file: the 348-byte header (dimensions,
//! datatype, voxel sizes, scaling, `sform` affine and the `n+1` magic), the
//! 4-byte extension flag, and the voxel data in little-endian order starting at
//! `vox_offset = 352`.

use super::types::{NiftiDtype, NiftiReader};
use crate::{Error, Result};

/// Size of the fixed NIfTI-1 header.
const HEADER_SIZE: usize = 348;
/// Offset of the voxel data: the header followed by the 4-byte extender.
const VOX_OFFSET: usize = 352;
/// Longest `descrip` text, leaving room for the terminating NUL.
const DESCRIP_CAPACITY: usize = 80;
/// `xyzt_units` value for "spatial units are millimetres".
const UNITS_MM: u8 = 2;
/// `sform_code` for "scanner anatomical" coordinates.
const SFORM_SCANNER_ANAT: i16 = 1;
/// Compression level used for `.nii.gz` output.
const GZIP_LEVEL: u8 = 6;

/// NIfTI-1 writer.
///
/// Voxel values of the [`NiftiReader`] are *real* values.  They are stored as
/// `raw = (value - scl_inter) / scl_slope`, so a reader that applies the NIfTI
/// rule `value = raw * scl_slope + scl_inter` recovers them.  A slope of zero
/// means "not scaled" in NIfTI, in which case the values are stored as they are.
#[derive(Clone, Debug)]
pub struct NiftiWriter {
    /// Output path used by [`NiftiWriter::write_to_file`].
    pub output_path: String,
    /// Data scaling: slope.
    pub scl_slope: f64,
    /// Data scaling: intercept.
    pub scl_inter: f64,
    /// Compress the output with gzip (`.nii.gz`).
    pub gzip: bool,
}

impl NiftiWriter {
    /// Create a new NIfTI writer.
    pub fn new(output_path: &str) -> Self {
        Self {
            output_path: output_path.to_string(),
            scl_slope: 1.0,
            scl_inter: 0.0,
            gzip: false,
        }
    }

    /// Build the 348-byte NIfTI-1 header describing `reader`.
    ///
    /// The header holds `dim`, `datatype`/`bitpix` (from the reader's
    /// [`NiftiDtype`]), `pixdim`, `vox_offset = 352`, the scaling slope and
    /// intercept, millimetre units, the description, `sform_code = 1` with the
    /// rows of the reader's affine, and the `n+1` magic.
    ///
    /// # Errors
    ///
    /// Fails when the dimensions cannot be stored (`dim[0]` outside `1..=7`, a
    /// size of zero or above `i16::MAX`) or when a header value (voxel size,
    /// affine entry, scaling) is not a finite number representable as `f32`.
    pub fn build_header(&self, reader: &NiftiReader) -> Result<Vec<u8>> {
        let (ndim, _) = checked_layout(reader)?;
        let (datatype, bitpix) = dtype_codes(&reader.datatype);
        let mut header = vec![0u8; HEADER_SIZE];
        put_i32(&mut header, 0, HEADER_SIZE as i32);
        put_i16(&mut header, 40, ndim as i16);
        for i in 1..8 {
            let size = if i <= ndim { reader.dim[i] as i16 } else { 1 };
            put_i16(&mut header, 40 + 2 * i, size);
        }
        put_i16(&mut header, 70, datatype);
        put_i16(&mut header, 72, bitpix);
        let qfac = if reader.pixdim[0] < 0.0 { -1.0 } else { 1.0 };
        put_f32(&mut header, 76, qfac);
        for i in 1..8 {
            put_f32(
                &mut header,
                76 + 4 * i,
                header_f32(reader.pixdim[i], "pixdim")?,
            );
        }
        put_f32(&mut header, 108, VOX_OFFSET as f32);
        let (slope, inter) = self.header_scaling()?;
        put_f32(&mut header, 112, slope as f32);
        put_f32(&mut header, 116, inter as f32);
        header[123] = UNITS_MM;
        let descrip = truncate_on_char_boundary(&reader.descrip, DESCRIP_CAPACITY - 1);
        header[148..148 + descrip.len()].copy_from_slice(descrip.as_bytes());
        put_i16(&mut header, 254, SFORM_SCANNER_ANAT);
        for (row, base) in [280usize, 296, 312].into_iter().enumerate() {
            for col in 0..4 {
                let value = header_f32(reader.affine[4 * row + col], "affine")?;
                put_f32(&mut header, base + 4 * col, value);
            }
        }
        header[344..348].copy_from_slice(b"n+1\0");
        Ok(header)
    }

    /// Serialise `reader` as a complete NIfTI-1 file in memory.
    ///
    /// The result is the header, the 4 extension bytes and the voxel data
    /// stored with the reader's datatype (little-endian, x fastest).  With
    /// [`NiftiWriter::gzip`] set the whole stream is gzip-compressed.
    ///
    /// # Errors
    ///
    /// Everything [`NiftiWriter::build_header`] rejects, plus: `reader.data`
    /// does not hold exactly the product of the dimensions, a value is NaN (or
    /// infinite) for an integer datatype, a scaled value does not fit the
    /// integer datatype, or gzip compression fails.
    pub fn write_to_buffer(&self, reader: &NiftiReader) -> Result<Vec<u8>> {
        let (_, n_voxels) = checked_layout(reader)?;
        if reader.data.len() != n_voxels {
            return Err(Error::General(format!(
                "NIfTI volume has {} voxel values but dim describes {n_voxels}",
                reader.data.len()
            )));
        }
        let mut buf = self.build_header(reader)?;
        buf.extend_from_slice(&[0u8; VOX_OFFSET - HEADER_SIZE]);
        self.encode_voxels(reader, &mut buf)?;
        if self.gzip {
            return oxiarc_deflate::gzip_compress(&buf, GZIP_LEVEL)
                .map_err(|e| Error::General(format!("gzip compression failed: {e}")));
        }
        Ok(buf)
    }

    /// Serialise `reader` (see [`NiftiWriter::write_to_buffer`]) and write the
    /// bytes to [`NiftiWriter::output_path`].
    ///
    /// # Errors
    ///
    /// Everything [`NiftiWriter::write_to_buffer`] rejects, plus file errors.
    pub fn write_to_file(&self, reader: &NiftiReader) -> Result<()> {
        let bytes = self.write_to_buffer(reader)?;
        std::fs::write(&self.output_path, bytes)?;
        Ok(())
    }

    /// Slope and intercept as stored in the header (rounded to `f32`).
    fn header_scaling(&self) -> Result<(f64, f64)> {
        let slope = header_f32(self.scl_slope, "scl_slope")?;
        let inter = header_f32(self.scl_inter, "scl_inter")?;
        Ok((f64::from(slope), f64::from(inter)))
    }

    /// Append the encoded voxel values of `reader` to `out`.
    fn encode_voxels(&self, reader: &NiftiReader, out: &mut Vec<u8>) -> Result<()> {
        let (mut slope, mut inter) = self.header_scaling()?;
        if slope == 0.0 {
            // NIfTI: a zero slope disables scaling altogether.
            slope = 1.0;
            inter = 0.0;
        }
        let raw = |v: f64| (v - inter) / slope;
        match reader.datatype {
            NiftiDtype::Float32 => {
                out.reserve(4 * reader.data.len());
                for &v in &reader.data {
                    out.extend_from_slice(&(raw(v) as f32).to_le_bytes());
                }
            }
            NiftiDtype::Float64 => {
                out.reserve(8 * reader.data.len());
                for &v in &reader.data {
                    out.extend_from_slice(&raw(v).to_le_bytes());
                }
            }
            NiftiDtype::Uint8 => {
                out.reserve(reader.data.len());
                for (idx, &v) in reader.data.iter().enumerate() {
                    out.push(integer_raw(raw(v), idx, "uint8", 0.0, f64::from(u8::MAX))? as u8);
                }
            }
            NiftiDtype::Int16 => {
                out.reserve(2 * reader.data.len());
                for (idx, &v) in reader.data.iter().enumerate() {
                    let r = integer_raw(
                        raw(v),
                        idx,
                        "int16",
                        f64::from(i16::MIN),
                        f64::from(i16::MAX),
                    )?;
                    out.extend_from_slice(&(r as i16).to_le_bytes());
                }
            }
            NiftiDtype::Int32 => {
                out.reserve(4 * reader.data.len());
                for (idx, &v) in reader.data.iter().enumerate() {
                    let r = integer_raw(
                        raw(v),
                        idx,
                        "int32",
                        f64::from(i32::MIN),
                        f64::from(i32::MAX),
                    )?;
                    out.extend_from_slice(&(r as i32).to_le_bytes());
                }
            }
        }
        Ok(())
    }
}

/// Validate the dimension fields of `reader`; returns `(ndim, n_voxels)`.
fn checked_layout(reader: &NiftiReader) -> Result<(usize, usize)> {
    let ndim = reader.dim[0];
    if !(1..=7).contains(&ndim) {
        return Err(Error::General(format!(
            "NIfTI-1 stores 1 to 7 dimensions, dim[0] is {ndim}"
        )));
    }
    let mut n_voxels = 1usize;
    for (axis, &size) in reader.dim.iter().enumerate().take(ndim + 1).skip(1) {
        if size == 0 || size > i16::MAX as usize {
            return Err(Error::General(format!(
                "NIfTI-1 dimension dim[{axis}] = {size} is outside 1..={}",
                i16::MAX
            )));
        }
        n_voxels = n_voxels
            .checked_mul(size)
            .ok_or_else(|| Error::General("NIfTI volume size overflows".to_string()))?;
    }
    Ok((ndim, n_voxels))
}

/// NIfTI `(datatype, bitpix)` codes for a voxel type.
fn dtype_codes(dtype: &NiftiDtype) -> (i16, i16) {
    match dtype {
        NiftiDtype::Uint8 => (2, 8),
        NiftiDtype::Int16 => (4, 16),
        NiftiDtype::Int32 => (8, 32),
        NiftiDtype::Float32 => (16, 32),
        NiftiDtype::Float64 => (64, 64),
    }
}

/// Convert a header value to `f32`, rejecting NaN, infinities and overflow.
fn header_f32(value: f64, what: &str) -> Result<f32> {
    let narrowed = value as f32;
    if value.is_finite() && narrowed.is_finite() {
        Ok(narrowed)
    } else {
        Err(Error::General(format!(
            "NIfTI header value {what} = {value} is not a finite f32"
        )))
    }
}

/// Round a scaled value to the nearest integer and check it fits `[lo, hi]`.
fn integer_raw(value: f64, index: usize, type_name: &str, lo: f64, hi: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::General(format!(
            "voxel {index} is {value} after scaling and cannot be stored as {type_name}"
        )));
    }
    let rounded = value.round();
    if !(lo..=hi).contains(&rounded) {
        return Err(Error::General(format!(
            "voxel {index} = {value} after scaling is outside the {type_name} range"
        )));
    }
    Ok(rounded)
}

/// Longest prefix of `text` that is at most `max_bytes` long and ends on a
/// character boundary.
fn truncate_on_char_boundary(text: &str, max_bytes: usize) -> &str {
    let mut end = text.len().min(max_bytes);
    while !text.is_char_boundary(end) {
        end -= 1;
    }
    &text[..end]
}

fn put_i16(buf: &mut [u8], offset: usize, value: i16) {
    buf[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_i32(buf: &mut [u8], offset: usize, value: i32) {
    buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_f32(buf: &mut [u8], offset: usize, value: f32) {
    buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::medical_imaging_io::NiftiReader as ParsedNifti;

    fn i16_at(b: &[u8], o: usize) -> i16 {
        i16::from_le_bytes([b[o], b[o + 1]])
    }

    fn i32_at(b: &[u8], o: usize) -> i32 {
        i32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])
    }

    fn f32_at(b: &[u8], o: usize) -> f32 {
        f32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])
    }

    /// A 4 x 3 x 2 volume with distinct voxel values and a non-trivial affine.
    fn sample_volume(datatype: NiftiDtype) -> NiftiReader {
        let mut r = NiftiReader::new("sample.nii");
        r.dim = [3, 4, 3, 2, 1, 1, 1, 1];
        r.pixdim = [1.0, 0.5, 0.75, 2.0, 1.0, 1.0, 1.0, 1.0];
        r.datatype = datatype;
        r.affine = [
            0.5, 0.0, 0.0, -10.0, //
            0.0, 0.75, 0.0, 20.0, //
            0.0, 0.0, 2.0, 30.0, //
            0.0, 0.0, 0.0, 1.0,
        ];
        r.descrip = "unit test volume".to_string();
        r.data = (0..24).map(|i| i as f64 * 3.0 - 6.0).collect();
        r
    }

    #[test]
    fn header_is_a_valid_nifti1_header() {
        let mut writer = NiftiWriter::new("out.nii");
        writer.scl_slope = 2.0;
        writer.scl_inter = 10.0;
        let reader = sample_volume(NiftiDtype::Float32);
        let h = writer.build_header(&reader).unwrap();
        assert_eq!(h.len(), 348);
        assert_eq!(i32_at(&h, 0), 348);
        assert_eq!(
            [40, 42, 44, 46, 48, 50, 52, 54].map(|o| i16_at(&h, o)),
            [3, 4, 3, 2, 1, 1, 1, 1]
        );
        assert_eq!(i16_at(&h, 70), 16, "datatype float32");
        assert_eq!(i16_at(&h, 72), 32, "bitpix");
        assert_eq!(f32_at(&h, 76), 1.0, "qfac");
        assert_eq!(f32_at(&h, 80), 0.5);
        assert_eq!(f32_at(&h, 84), 0.75);
        assert_eq!(f32_at(&h, 88), 2.0);
        assert_eq!(f32_at(&h, 108), 352.0, "vox_offset");
        assert_eq!(f32_at(&h, 112), 2.0, "scl_slope");
        assert_eq!(f32_at(&h, 116), 10.0, "scl_inter");
        assert_eq!(h[123], 2, "millimetres");
        assert_eq!(&h[148..164], b"unit test volume");
        assert_eq!(h[164], 0, "descrip is NUL terminated");
        assert_eq!(i16_at(&h, 252), 0, "qform_code");
        assert_eq!(i16_at(&h, 254), 1, "sform_code");
        assert_eq!(
            [280, 284, 288, 292].map(|o| f32_at(&h, o)),
            [0.5, 0.0, 0.0, -10.0]
        );
        assert_eq!(
            [296, 300, 304, 308].map(|o| f32_at(&h, o)),
            [0.0, 0.75, 0.0, 20.0]
        );
        assert_eq!(
            [312, 316, 320, 324].map(|o| f32_at(&h, o)),
            [0.0, 0.0, 2.0, 30.0]
        );
        assert_eq!(&h[344..348], b"n+1\0");
    }

    #[test]
    fn buffer_has_extension_bytes_and_data_at_offset_352() {
        let reader = sample_volume(NiftiDtype::Float32);
        let buf = NiftiWriter::new("o.nii").write_to_buffer(&reader).unwrap();
        assert_eq!(buf.len(), 352 + 24 * 4);
        assert_eq!(&buf[348..352], &[0, 0, 0, 0]);
        assert_eq!(f32_at(&buf, 352), -6.0);
        assert_eq!(f32_at(&buf, 352 + 4 * 23), 63.0);
    }

    #[test]
    fn float32_round_trips_through_the_nifti_parser() {
        let mut writer = NiftiWriter::new("o.nii");
        writer.scl_slope = 2.5;
        writer.scl_inter = -3.0;
        let reader = sample_volume(NiftiDtype::Float32);
        let buf = writer.write_to_buffer(&reader).unwrap();
        let mut parsed = ParsedNifti::new();
        parsed.parse(&buf).unwrap();
        let h = &parsed.header;
        assert_eq!(h.dim[..4], [3, 4, 3, 2]);
        assert_eq!(h.datatype, 16);
        assert_eq!(h.bitpix, 32);
        assert_eq!(h.vox_offset, 352.0);
        assert_eq!(h.pixdim[1..4], [0.5, 0.75, 2.0]);
        assert_eq!(h.sform_code, 1);
        assert_eq!(h.srow_x, [0.5, 0.0, 0.0, -10.0]);
        assert_eq!(h.srow_y, [0.0, 0.75, 0.0, 20.0]);
        assert_eq!(h.srow_z, [0.0, 0.0, 2.0, 30.0]);
        assert_eq!(h.magic, "n+1");
        assert_eq!(parsed.voxels.len(), 24);
        // The parser returns raw values; applying the header scaling gives the input back.
        for (raw, want) in parsed.voxels.iter().zip(&reader.data) {
            let real = f64::from(*raw) * f64::from(h.scl_slope) + f64::from(h.scl_inter);
            assert!((real - want).abs() < 1e-5, "{real} vs {want}");
        }
        // Scaling is applied in the NIfTI direction: raw = (value - inter) / slope.
        assert!((f64::from(parsed.voxels[0]) - (-6.0 + 3.0) / 2.5).abs() < 1e-6);
    }

    #[test]
    fn integer_and_double_datatypes_round_trip() {
        let writer = NiftiWriter::new("o.nii");
        for (dtype, code, bitpix) in [
            (NiftiDtype::Uint8, 2i16, 8i16),
            (NiftiDtype::Int16, 4, 16),
            (NiftiDtype::Float64, 64, 64),
        ] {
            let mut reader = sample_volume(dtype);
            // Values 0, 7, 14, ... fit every type and are exact in f32.
            reader.data = (0..24).map(|i| (7 * i) as f64).collect();
            let buf = writer.write_to_buffer(&reader).unwrap();
            let mut parsed = ParsedNifti::new();
            parsed.parse(&buf).unwrap();
            assert_eq!(parsed.header.datatype, code);
            assert_eq!(parsed.header.bitpix, bitpix);
            assert_eq!(parsed.voxels.len(), 24);
            for (got, want) in parsed.voxels.iter().zip(&reader.data) {
                assert_eq!(f64::from(*got), *want);
            }
        }
    }

    #[test]
    fn negative_int16_values_are_stored_signed() {
        let mut reader = sample_volume(NiftiDtype::Int16);
        reader.data = (0..24).map(|i| i as f64 * 100.0 - 1000.0).collect();
        let buf = NiftiWriter::new("o.nii").write_to_buffer(&reader).unwrap();
        let mut parsed = ParsedNifti::new();
        parsed.parse(&buf).unwrap();
        assert_eq!(parsed.voxels[0], -1000.0);
        assert_eq!(parsed.voxels[23], 1300.0);
    }

    #[test]
    fn int32_is_written_with_datatype_8() {
        let mut reader = sample_volume(NiftiDtype::Int32);
        reader.data = (0..24).map(|i| i as f64 * 1.0e6 - 5.0e6).collect();
        let buf = NiftiWriter::new("o.nii").write_to_buffer(&reader).unwrap();
        assert_eq!(i16_at(&buf, 70), 8);
        assert_eq!(i16_at(&buf, 72), 32);
        assert_eq!(buf.len(), 352 + 24 * 4);
        assert_eq!(i32_at(&buf, 352), -5_000_000);
        assert_eq!(i32_at(&buf, 352 + 4 * 23), 18_000_000);
    }

    #[test]
    fn integer_types_round_scaled_values_to_nearest() {
        let mut writer = NiftiWriter::new("o.nii");
        writer.scl_slope = 0.5;
        writer.scl_inter = 1.0;
        let mut reader = sample_volume(NiftiDtype::Int16);
        reader.dim = [1, 4, 1, 1, 1, 1, 1, 1];
        reader.data = vec![1.0, 1.5, 2.0, 11.2];
        let buf = writer.write_to_buffer(&reader).unwrap();
        // raw = (v - 1) / 0.5 = 0, 1, 2, 20.4 -> 20
        let raws: Vec<i16> = (0..4).map(|i| i16_at(&buf, 352 + 2 * i)).collect();
        assert_eq!(raws, vec![0, 1, 2, 20]);
    }

    #[test]
    fn zero_slope_means_unscaled() {
        let mut writer = NiftiWriter::new("o.nii");
        writer.scl_slope = 0.0;
        writer.scl_inter = 5.0;
        let reader = sample_volume(NiftiDtype::Float32);
        let buf = writer.write_to_buffer(&reader).unwrap();
        assert_eq!(f32_at(&buf, 112), 0.0);
        assert_eq!(f32_at(&buf, 352), -6.0);
        assert_eq!(f32_at(&buf, 352 + 4 * 5), 9.0);
    }

    #[test]
    fn out_of_range_and_non_finite_integer_voxels_are_errors() {
        let writer = NiftiWriter::new("o.nii");
        let mut reader = sample_volume(NiftiDtype::Uint8);
        reader.data = vec![0.0; 24];
        reader.data[5] = 300.0;
        let err = writer.write_to_buffer(&reader).unwrap_err().to_string();
        assert!(err.contains("voxel 5") && err.contains("uint8"), "{err}");
        reader.data[5] = -1.0;
        assert!(writer.write_to_buffer(&reader).is_err());
        reader.datatype = NiftiDtype::Int16;
        reader.data[5] = 40_000.0;
        assert!(writer.write_to_buffer(&reader).is_err());
        reader.data[5] = f64::NAN;
        assert!(writer.write_to_buffer(&reader).is_err());
        reader.datatype = NiftiDtype::Int32;
        reader.data[5] = f64::INFINITY;
        assert!(writer.write_to_buffer(&reader).is_err());
        // Floating-point storage keeps NaN.
        reader.datatype = NiftiDtype::Float32;
        reader.data[5] = f64::NAN;
        let buf = writer.write_to_buffer(&reader).unwrap();
        assert!(f32_at(&buf, 352 + 4 * 5).is_nan());
    }

    #[test]
    fn data_length_must_match_the_dimensions() {
        let writer = NiftiWriter::new("o.nii");
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.data.pop();
        let err = writer.write_to_buffer(&reader).unwrap_err().to_string();
        assert!(err.contains("23") && err.contains("24"), "{err}");
        reader.data.clear();
        assert!(writer.write_to_buffer(&reader).is_err());
        // The header alone does not need voxel data.
        assert_eq!(writer.build_header(&reader).unwrap().len(), 348);
    }

    #[test]
    fn invalid_dimensions_are_errors() {
        let writer = NiftiWriter::new("o.nii");
        for dim in [
            [0, 4, 3, 2, 1, 1, 1, 1],
            [8, 4, 3, 2, 1, 1, 1, 1],
            [3, 4, 0, 2, 1, 1, 1, 1],
            [3, 4, 3, 40_000, 1, 1, 1, 1],
        ] {
            let mut reader = sample_volume(NiftiDtype::Float32);
            reader.dim = dim;
            assert!(writer.build_header(&reader).is_err(), "{dim:?}");
            assert!(writer.write_to_buffer(&reader).is_err(), "{dim:?}");
        }
    }

    #[test]
    fn non_finite_header_values_are_errors() {
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.affine[3] = f64::NAN;
        assert!(NiftiWriter::new("o.nii").build_header(&reader).is_err());
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.pixdim[2] = f64::INFINITY;
        assert!(NiftiWriter::new("o.nii").build_header(&reader).is_err());
        let mut writer = NiftiWriter::new("o.nii");
        writer.scl_slope = f64::NAN;
        assert!(
            writer
                .build_header(&sample_volume(NiftiDtype::Float32))
                .is_err()
        );
        writer.scl_slope = 1.0;
        writer.scl_inter = 1.0e300;
        assert!(
            writer
                .build_header(&sample_volume(NiftiDtype::Float32))
                .is_err()
        );
    }

    #[test]
    fn unused_dimensions_are_normalised_to_one() {
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.dim = [2, 4, 6, 99, 77, 0, 0, 0];
        reader.data = vec![0.0; 24];
        let h = NiftiWriter::new("o.nii").build_header(&reader).unwrap();
        assert_eq!(
            [40, 42, 44, 46, 48, 50, 52, 54].map(|o| i16_at(&h, o)),
            [2, 4, 6, 1, 1, 1, 1, 1]
        );
    }

    #[test]
    fn four_dimensional_volumes_write_every_frame() {
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.dim = [4, 2, 2, 2, 3, 1, 1, 1];
        reader.data = (0..24).map(f64::from).collect();
        let writer = NiftiWriter::new("o.nii");
        let buf = writer.write_to_buffer(&reader).unwrap();
        assert_eq!(i16_at(&buf, 40), 4);
        assert_eq!(i16_at(&buf, 48), 3);
        assert_eq!(buf.len(), 352 + 24 * 4);
        assert_eq!(f32_at(&buf, 352 + 4 * 23), 23.0);
        reader.data.pop();
        assert!(writer.write_to_buffer(&reader).is_err());
    }

    #[test]
    fn negative_pixdim0_is_written_as_qfac_minus_one() {
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.pixdim[0] = -1.0;
        let h = NiftiWriter::new("o.nii").build_header(&reader).unwrap();
        assert_eq!(f32_at(&h, 76), -1.0);
    }

    #[test]
    fn long_multibyte_description_is_cut_on_a_character_boundary() {
        let mut reader = sample_volume(NiftiDtype::Float32);
        reader.descrip = "画".repeat(40);
        let h = NiftiWriter::new("o.nii").build_header(&reader).unwrap();
        let field = &h[148..228];
        let end = field.iter().position(|&b| b == 0).unwrap();
        assert_eq!(end, 78, "26 three-byte characters fit in 79 bytes");
        assert_eq!(std::str::from_utf8(&field[..end]).unwrap(), "画".repeat(26));
        assert!(field[end..].iter().all(|&b| b == 0));
    }

    #[test]
    fn gzip_output_decompresses_to_the_plain_stream() {
        let reader = sample_volume(NiftiDtype::Float32);
        let plain = NiftiWriter::new("o.nii").write_to_buffer(&reader).unwrap();
        let mut gz_writer = NiftiWriter::new("o.nii.gz");
        gz_writer.gzip = true;
        let gz = gz_writer.write_to_buffer(&reader).unwrap();
        assert_eq!(&gz[..2], &[0x1f, 0x8b], "gzip magic");
        assert_ne!(gz, plain);
        let restored = oxiarc_deflate::gzip_decompress(&gz).unwrap();
        assert_eq!(restored, plain);
        let mut parsed = ParsedNifti::new();
        parsed.parse(&restored).unwrap();
        assert_eq!(parsed.voxels.len(), 24);
    }

    #[test]
    fn write_to_file_stores_the_buffer() {
        let path = std::env::temp_dir().join(format!(
            "oxiphysics_io_nifti_writer_{}.nii",
            std::process::id()
        ));
        let writer = NiftiWriter::new(path.to_str().unwrap());
        let reader = sample_volume(NiftiDtype::Int16);
        writer.write_to_file(&reader).unwrap();
        let on_disk = std::fs::read(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(on_disk, writer.write_to_buffer(&reader).unwrap());
    }

    #[test]
    fn write_to_file_reports_io_errors() {
        let dir = std::env::temp_dir().join("oxiphysics_io_nifti_missing_dir");
        let path = dir.join("nested").join("out.nii");
        let writer = NiftiWriter::new(path.to_str().unwrap());
        let err = writer
            .write_to_file(&sample_volume(NiftiDtype::Float32))
            .unwrap_err();
        assert!(matches!(err, Error::Io(_)), "{err:?}");
    }
}
