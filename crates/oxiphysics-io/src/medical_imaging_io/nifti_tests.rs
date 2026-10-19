// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the NIfTI-1 reader: every scalar datatype decodes to its real
//! values, and anything the reader cannot decode is an error instead of a
//! volume of zeros.

use super::types::NiftiReader;

/// Small deterministic xorshift generator (the crate has no `rand` dependency).
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform index in `0..n` (`n` must not be zero).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Fields of a synthetic NIfTI-1 header.
struct Spec {
    ndim: i16,
    dims: [i16; 4],
    datatype: i16,
    bitpix: i16,
    vox_offset: f32,
}

impl Spec {
    fn new(dims: [i16; 3], datatype: i16, bitpix: i16) -> Self {
        Self {
            ndim: 3,
            dims: [dims[0], dims[1], dims[2], 1],
            datatype,
            bitpix,
            vox_offset: 352.0,
        }
    }

    /// A 348-byte header plus the 4-byte extender, followed by `payload`
    /// (padded with zeros up to `vox_offset` when that is larger than 352).
    fn build(&self, payload: &[u8]) -> Vec<u8> {
        let mut buf = vec![0u8; 352];
        buf[0..4].copy_from_slice(&348_i32.to_le_bytes());
        buf[40..42].copy_from_slice(&self.ndim.to_le_bytes());
        for (i, d) in self.dims.iter().enumerate() {
            buf[42 + 2 * i..44 + 2 * i].copy_from_slice(&d.to_le_bytes());
        }
        buf[70..72].copy_from_slice(&self.datatype.to_le_bytes());
        buf[72..74].copy_from_slice(&self.bitpix.to_le_bytes());
        for i in 0..4 {
            buf[76 + 4 * i..80 + 4 * i].copy_from_slice(&1.0_f32.to_le_bytes());
        }
        buf[108..112].copy_from_slice(&self.vox_offset.to_le_bytes());
        buf[112..116].copy_from_slice(&1.0_f32.to_le_bytes());
        buf[344..348].copy_from_slice(b"n+1\0");
        if self.vox_offset.is_finite() && self.vox_offset > 352.0 && self.vox_offset < 1.0e6 {
            buf.resize(self.vox_offset as usize, 0);
        }
        buf.extend_from_slice(payload);
        buf
    }
}

fn parse(buf: &[u8]) -> Result<NiftiReader, String> {
    let mut reader = NiftiReader::new();
    reader.parse(buf)?;
    Ok(reader)
}

fn error_of(buf: &[u8]) -> String {
    match parse(buf) {
        Ok(_) => panic!("expected the NIfTI parse to fail"),
        Err(message) => message,
    }
}

fn bytes_of<T: Copy, const N: usize>(values: &[T], to_bytes: impl Fn(T) -> [u8; N]) -> Vec<u8> {
    values.iter().flat_map(|&v| to_bytes(v)).collect()
}

fn assert_voxels(reader: &NiftiReader, want: &[f32]) {
    assert_eq!(reader.voxels, want);
}

// ---------------------------------------------------------------------------
// datatypes
// ---------------------------------------------------------------------------

#[test]
fn uint16_volume_decodes_to_its_values() {
    let values = [0u16, 1, 255, 256, 40_000, 65_535, 12, 7];
    let payload = bytes_of(&values, u16::to_le_bytes);
    let reader = parse(&Spec::new([2, 2, 2], 512, 16).build(&payload)).unwrap();
    assert_voxels(
        &reader,
        &[0.0, 1.0, 255.0, 256.0, 40_000.0, 65_535.0, 12.0, 7.0],
    );
}

#[test]
fn int32_volume_decodes_to_its_values() {
    let values = [-5_000_000_i32, 16_777_216, -1, 0, 123_456, i32::MIN, 3, 4];
    let payload = bytes_of(&values, i32::to_le_bytes);
    let reader = parse(&Spec::new([4, 2, 1], 8, 32).build(&payload)).unwrap();
    assert_voxels(
        &reader,
        &[
            -5_000_000.0,
            16_777_216.0,
            -1.0,
            0.0,
            123_456.0,
            -2_147_483_648.0,
            3.0,
            4.0,
        ],
    );
}

#[test]
fn every_scalar_datatype_decodes() {
    let n = 3;
    let dims = [3, 1, 1];
    // uint8
    let reader = parse(&Spec::new(dims, 2, 8).build(&[0, 128, 255])).unwrap();
    assert_voxels(&reader, &[0.0, 128.0, 255.0]);
    // int8
    let payload = bytes_of(&[-128_i8, -1, 127], i8::to_le_bytes);
    let reader = parse(&Spec::new(dims, 256, 8).build(&payload)).unwrap();
    assert_voxels(&reader, &[-128.0, -1.0, 127.0]);
    // int16
    let payload = bytes_of(&[-32_768_i16, -2, 32_767], i16::to_le_bytes);
    let reader = parse(&Spec::new(dims, 4, 16).build(&payload)).unwrap();
    assert_voxels(&reader, &[-32_768.0, -2.0, 32_767.0]);
    // uint32
    let payload = bytes_of(&[0_u32, 4_000_000_000, 7], u32::to_le_bytes);
    let reader = parse(&Spec::new(dims, 768, 32).build(&payload)).unwrap();
    assert_voxels(&reader, &[0.0, 4_000_000_000.0, 7.0]);
    // float32
    let payload = bytes_of(&[-1.5_f32, 0.0, 2.25], f32::to_le_bytes);
    let reader = parse(&Spec::new(dims, 16, 32).build(&payload)).unwrap();
    assert_voxels(&reader, &[-1.5, 0.0, 2.25]);
    // float64 is rounded to f32
    let payload = bytes_of(&[-1.5_f64, 0.1, 1.0e10], f64::to_le_bytes);
    let reader = parse(&Spec::new(dims, 64, 64).build(&payload)).unwrap();
    assert_voxels(&reader, &[-1.5, 0.1_f64 as f32, 1.0e10_f64 as f32]);
    // int64 / uint64 values that f32 holds exactly
    let payload = bytes_of(&[-(1_i64 << 52), 0, 1 << 40], i64::to_le_bytes);
    let reader = parse(&Spec::new(dims, 1024, 64).build(&payload)).unwrap();
    assert_voxels(&reader, &[-(2.0_f32.powi(52)), 0.0, 2.0_f32.powi(40)]);
    let payload = bytes_of(&[1_u64 << 63, 0, 3], u64::to_le_bytes);
    let reader = parse(&Spec::new(dims, 1280, 64).build(&payload)).unwrap();
    assert_voxels(&reader, &[2.0_f32.powi(63), 0.0, 3.0]);
    assert_eq!(reader.voxels.len(), n);
}

#[test]
fn unsupported_datatypes_are_errors_not_zero_volumes() {
    let payload = vec![0u8; 4096];
    for (datatype, name) in [
        (0, "unknown"),
        (1, "binary"),
        (32, "complex64"),
        (128, "rgb24"),
        (1536, "float128"),
        (1792, "complex128"),
        (2048, "complex256"),
        (2304, "rgba32"),
        (3, "not a NIfTI-1 datatype"),
        (9999, "not a NIfTI-1 datatype"),
        (-16, "not a NIfTI-1 datatype"),
    ] {
        let message = error_of(&Spec::new([2, 2, 2], datatype, 32).build(&payload));
        assert!(
            message.contains(&format!("unsupported NIfTI datatype {datatype}")),
            "{message}"
        );
        assert!(message.contains(name), "{message}");
    }
}

#[test]
fn bitpix_must_match_the_datatype() {
    let payload = vec![0u8; 64];
    let message = error_of(&Spec::new([2, 2, 2], 512, 8).build(&payload));
    assert!(message.contains("bitpix 8"), "{message}");
    let message = error_of(&Spec::new([2, 2, 2], 16, 64).build(&payload));
    assert!(message.contains("bitpix 64"), "{message}");
    assert!(parse(&Spec::new([2, 2, 2], 16, 32).build(&payload)).is_ok());
}

// ---------------------------------------------------------------------------
// sizes and offsets
// ---------------------------------------------------------------------------

#[test]
fn truncated_voxel_data_is_an_error() {
    let values: Vec<f32> = (0..8).map(|i| i as f32).collect();
    let payload = bytes_of(&values, f32::to_le_bytes);
    let spec = Spec::new([2, 2, 2], 16, 32);
    assert!(parse(&spec.build(&payload)).is_ok());
    let message = error_of(&spec.build(&payload[..payload.len() - 1]));
    assert!(message.contains("truncated"), "{message}");
    assert!(message.contains("32 bytes"), "{message}");
    // Bytes after the volume are ignored.
    let mut longer = payload.clone();
    longer.extend_from_slice(&[9; 100]);
    assert_eq!(parse(&spec.build(&longer)).unwrap().voxels.len(), 8);
    // A header with no data at all.
    assert!(parse(&spec.build(&[])).is_err());
}

#[test]
fn huge_dimensions_in_a_tiny_file_are_rejected_without_allocating() {
    for datatype_bitpix in [(2, 8), (16, 32), (64, 64)] {
        let spec = Spec::new([i16::MAX; 3], datatype_bitpix.0, datatype_bitpix.1);
        let message = error_of(&spec.build(&[0; 64]));
        assert!(message.contains("truncated"), "{message}");
    }
}

#[test]
fn vox_offset_must_be_a_whole_number_inside_the_file() {
    let payload = vec![0u8; 32];
    let mut spec = Spec::new([2, 2, 2], 16, 32);
    for bad in [
        f32::NAN,
        f32::INFINITY,
        f32::NEG_INFINITY,
        -352.0,
        -0.0,
        0.0,
        347.0,
        352.5,
        1.0e9,
    ] {
        spec.vox_offset = bad;
        let message = error_of(&spec.build(&payload));
        assert!(
            message.contains("vox_offset"),
            "vox_offset {bad}: {message}"
        );
    }
    // An offset beyond the data is the same error as truncated data.
    spec.vox_offset = 400.0;
    let short = spec.build(&payload)[..380].to_vec();
    assert!(error_of(&short).contains("vox_offset"));
}

#[test]
fn voxels_are_read_from_vox_offset() {
    let values = [10.0_f32, 20.0, 30.0, 40.0];
    let payload = bytes_of(&values, f32::to_le_bytes);
    for offset in [352.0, 368.0, 400.0] {
        let mut spec = Spec::new([2, 2, 1], 16, 32);
        spec.vox_offset = offset;
        let reader = parse(&spec.build(&payload)).unwrap();
        assert_voxels(&reader, &values);
        assert_eq!(reader.header.vox_offset, offset);
    }
}

#[test]
fn dimensions_are_validated() {
    let payload = vec![0u8; 256];
    for ndim in [0_i16, -1, 8, 100] {
        let mut spec = Spec::new([2, 2, 2], 16, 32);
        spec.ndim = ndim;
        let message = error_of(&spec.build(&payload));
        assert!(message.contains("dim[0]"), "{message}");
    }
    for bad in [[0, 2, 2], [2, 0, 2], [2, 2, 0], [-3, 2, 2]] {
        let message = error_of(&Spec::new(bad, 16, 32).build(&payload));
        assert!(message.contains("not positive"), "{message}");
    }
    // Sizes beyond dim[0] are not checked: they do not describe an axis.
    let mut spec = Spec::new([2, 2, 2], 16, 32);
    spec.ndim = 2;
    spec.dims = [2, 2, 0, 0];
    let reader = parse(&spec.build(&payload)).unwrap();
    assert_eq!(reader.header.volume_dims(), [2, 2, 1]);
    assert_eq!(reader.voxels.len(), 4);
}

#[test]
fn a_fourth_dimension_reads_only_the_first_volume() {
    let values: Vec<f32> = (0..24).map(|i| i as f32).collect();
    let payload = bytes_of(&values, f32::to_le_bytes);
    let mut spec = Spec::new([2, 2, 2], 16, 32);
    spec.ndim = 4;
    spec.dims = [2, 2, 2, 3];
    let reader = parse(&spec.build(&payload)).unwrap();
    assert_eq!(reader.header.dim[..5], [4, 2, 2, 2, 3]);
    assert_voxels(&reader, &values[..8]);
    // The data of the first volume must be there; later volumes are not read.
    let cut = &payload[..8 * 4];
    assert_eq!(parse(&spec.build(cut)).unwrap().voxels.len(), 8);
}

#[test]
fn byte_swapped_headers_are_reported_as_such() {
    let mut buf = Spec::new([2, 2, 2], 16, 32).build(&[0; 32]);
    buf[0..4].copy_from_slice(&348_i32.to_be_bytes());
    let message = error_of(&buf);
    assert!(message.contains("big-endian"), "{message}");
    buf[0..4].copy_from_slice(&540_i32.to_le_bytes());
    let message = error_of(&buf);
    assert!(message.contains("sizeof_hdr"), "{message}");
    assert!(parse(&buf[..347]).is_err());
}

// ---------------------------------------------------------------------------
// reader state
// ---------------------------------------------------------------------------

#[test]
fn a_failed_parse_leaves_the_reader_unchanged() {
    let values = [1.0_f32, 2.0, 3.0, 4.0];
    let good = Spec::new([2, 2, 1], 16, 32).build(&bytes_of(&values, f32::to_le_bytes));
    let mut reader = NiftiReader::new();
    reader.parse(&good).unwrap();
    let header_before = reader.header.clone();
    // A header that is acceptable until the datatype check, then truncated data.
    let bad = Spec::new([i16::MAX; 3], 16, 32).build(&[0; 16]);
    assert!(reader.parse(&bad).is_err());
    assert_eq!(reader.voxels, values);
    assert_eq!(reader.header.dim, header_before.dim);
    assert_eq!(reader.header.vox_offset, header_before.vox_offset);
    // The reader still converts to the volume that was read first.
    let volume = reader.to_volume();
    assert_eq!(volume.dims, [2, 2, 1]);
    assert_eq!(volume.data, values);
}

#[test]
fn to_volume_never_makes_up_voxels() {
    let values = [5.0_f32, 6.0, 7.0, 8.0];
    let reader =
        parse(&Spec::new([2, 2, 1], 16, 32).build(&bytes_of(&values, f32::to_le_bytes))).unwrap();
    let volume = reader.to_volume();
    assert_eq!(volume.dims, [2, 2, 1]);
    assert_eq!(volume.data, values);
    // A reader that parsed nothing gives an empty volume.
    let empty = NiftiReader::new().to_volume();
    assert_eq!(empty.dims, [0, 0, 0]);
    assert!(empty.data.is_empty());
    // A header edited to enormous dimensions does not make `to_volume` allocate.
    let mut edited = NiftiReader::new();
    edited.header.dim = [3, i16::MAX, i16::MAX, i16::MAX, 1, 1, 1, 1];
    let volume = edited.to_volume();
    assert_eq!(volume.dims, [0, 0, 0]);
    assert!(volume.data.is_empty());
}

#[test]
fn the_description_is_read_from_the_header() {
    let mut buf = Spec::new([1, 1, 1], 16, 32).build(&1.0_f32.to_le_bytes());
    buf[148..148 + 11].copy_from_slice(b"unit  test ");
    let reader = parse(&buf).unwrap();
    assert_eq!(reader.header.descrip, "unit  test");
    assert_eq!(reader.header.magic, "n+1");
}

#[test]
fn num_voxels_saturates_instead_of_overflowing() {
    let mut reader = NiftiReader::new();
    reader.header.dim = [3, i16::MAX, i16::MAX, i16::MAX, 1, 1, 1, 1];
    assert_eq!(
        reader.header.num_voxels(),
        i16::MAX as usize * i16::MAX as usize * i16::MAX as usize
    );
}

#[test]
fn corrupted_files_never_panic_and_never_decode_more_than_they_hold() {
    let values: Vec<f32> = (0..24).map(|i| i as f32).collect();
    let payload = bytes_of(&values, f32::to_le_bytes);
    let corpus = [
        Spec::new([4, 3, 2], 16, 32).build(&payload),
        Spec::new([4, 3, 2], 4, 16).build(&payload[..48]),
        Spec::new([4, 3, 2], 64, 64).build(&[payload.clone(), payload.clone()].concat()),
    ];
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    let (mut accepted, mut rejected) = (0_usize, 0_usize);
    for base in &corpus {
        for _ in 0..3000 {
            let mut bytes = base.clone();
            for _ in 0..1 + rng.below(4) {
                // Most of the damage lands in the header, where the sizes live.
                let limit = if rng.below(4) == 0 { bytes.len() } else { 352 };
                let at = rng.below(limit);
                bytes[at] = rng.next_u64() as u8;
            }
            if rng.below(5) == 0 {
                let keep = rng.below(bytes.len() + 1);
                bytes.truncate(keep);
            }
            let mut reader = NiftiReader::new();
            if reader.parse(&bytes).is_err() {
                rejected += 1;
                continue;
            }
            accepted += 1;
            let bytes_per_voxel = usize::try_from(reader.header.bitpix / 8).unwrap();
            let offset = reader.header.vox_offset as usize;
            assert_eq!(reader.voxels.len(), reader.header.num_voxels());
            assert!(
                offset + reader.voxels.len() * bytes_per_voxel <= bytes.len(),
                "{} voxels of {bytes_per_voxel} bytes from offset {offset} in {} bytes",
                reader.voxels.len(),
                bytes.len()
            );
            let volume = reader.to_volume();
            assert_eq!(volume.data.len(), reader.voxels.len());
        }
    }
    // The corruption reaches both outcomes, so the checks above are not vacuous.
    assert!(accepted > 1000, "{accepted} files were accepted");
    assert!(rejected > 1000, "{rejected} files were rejected");
}
