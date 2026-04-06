// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Custom binary format for physics field data.
//!
//! Provides a compact binary representation with a fixed header, run-length
//! encoding (RLE) for float sequences, CRC32 checksums, and integrity
//! verification.

// ── Header ────────────────────────────────────────────────────────────────────

/// Magic bytes that identify a physics binary field file.
pub const MAGIC: [u8; 4] = *b"OXPF";

/// Current format version.
pub const FORMAT_VERSION: u16 = 1;

/// Data type tag stored in the header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum FieldDataType {
    /// 64-bit IEEE 754 float.
    Float64 = 0,
    /// 32-bit IEEE 754 float.
    Float32 = 1,
    /// 32-bit signed integer.
    Int32 = 2,
}

impl FieldDataType {
    /// Attempt to parse a `u8` tag into a `FieldDataType`.
    pub fn from_u8(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::Float64),
            1 => Some(Self::Float32),
            2 => Some(Self::Int32),
            _ => None,
        }
    }
}

/// Header for a binary physics field record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BinaryHeader {
    /// Magic identifier (`OXPF`).
    pub magic: [u8; 4],
    /// Format version.
    pub version: u16,
    /// Element data type.
    pub data_type: FieldDataType,
    /// Number of dimensions (1, 2, or 3).
    pub ndim: u8,
    /// Dimension sizes (`[nx, ny, nz]`; unused dims are `1`).
    pub dims: [u32; 3],
    /// CRC32 checksum of the payload bytes.
    pub checksum: u32,
}

impl BinaryHeader {
    /// Serialised size in bytes.
    pub const SIZE: usize = 4 + 2 + 1 + 1 + 12 + 4; // 24 bytes

    /// Serialise the header to bytes.
    pub fn to_bytes(&self) -> [u8; Self::SIZE] {
        let mut buf = [0u8; Self::SIZE];
        buf[0..4].copy_from_slice(&self.magic);
        buf[4..6].copy_from_slice(&self.version.to_le_bytes());
        buf[6] = self.data_type as u8;
        buf[7] = self.ndim;
        buf[8..12].copy_from_slice(&self.dims[0].to_le_bytes());
        buf[12..16].copy_from_slice(&self.dims[1].to_le_bytes());
        buf[16..20].copy_from_slice(&self.dims[2].to_le_bytes());
        buf[20..24].copy_from_slice(&self.checksum.to_le_bytes());
        buf
    }

    /// Deserialise a header from a 24-byte slice.
    pub fn from_bytes(buf: &[u8]) -> Result<Self, String> {
        if buf.len() < Self::SIZE {
            return Err(format!(
                "Header too short: expected {} bytes, got {}",
                Self::SIZE,
                buf.len()
            ));
        }
        let magic: [u8; 4] = buf[0..4].try_into().expect("slice length must match");
        let version = u16::from_le_bytes(buf[4..6].try_into().expect("slice length must match"));
        let data_type = FieldDataType::from_u8(buf[6])
            .ok_or_else(|| format!("Unknown data type: {}", buf[6]))?;
        let ndim = buf[7];
        let d0 = u32::from_le_bytes(buf[8..12].try_into().expect("slice length must match"));
        let d1 = u32::from_le_bytes(buf[12..16].try_into().expect("slice length must match"));
        let d2 = u32::from_le_bytes(buf[16..20].try_into().expect("slice length must match"));
        let checksum = u32::from_le_bytes(buf[20..24].try_into().expect("slice length must match"));
        Ok(Self {
            magic,
            version,
            data_type,
            ndim,
            dims: [d0, d1, d2],
            checksum,
        })
    }

    /// Total number of elements (`dims[0] * dims[1] * dims[2]`).
    pub fn element_count(&self) -> usize {
        self.dims[0] as usize * self.dims[1] as usize * self.dims[2] as usize
    }
}

// ── Write / Read ──────────────────────────────────────────────────────────────

/// Encode a `f64` slice as raw little-endian bytes.
fn f64_slice_to_bytes(data: &[f64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(data.len() * 8);
    for &v in data {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    bytes
}

/// Decode a little-endian byte buffer into a `f64` Vec.
fn bytes_to_f64_slice(bytes: &[u8]) -> Result<Vec<f64>, String> {
    if !bytes.len().is_multiple_of(8) {
        return Err(format!(
            "Payload length {} is not a multiple of 8",
            bytes.len()
        ));
    }
    let mut out = Vec::with_capacity(bytes.len() / 8);
    for chunk in bytes.chunks_exact(8) {
        out.push(f64::from_le_bytes(
            chunk.try_into().expect("slice length must match"),
        ));
    }
    Ok(out)
}

/// Write a `f64` field with a binary header.
///
/// Returns the complete byte stream: `[header (24 bytes)] + [payload]`.
/// The header checksum is computed automatically.
///
/// * `data` – the field values in row-major order
/// * `dims` – `[nx, ny, nz]` (set unused dims to 1)
/// * `ndim` – number of active dimensions (1, 2, or 3)
pub fn write_binary_field(data: &[f64], dims: [u32; 3], ndim: u8) -> Vec<u8> {
    let payload = f64_slice_to_bytes(data);
    let crc = checksum_crc32(&payload);
    let header = BinaryHeader {
        magic: MAGIC,
        version: FORMAT_VERSION,
        data_type: FieldDataType::Float64,
        ndim,
        dims,
        checksum: crc,
    };
    let mut out = Vec::with_capacity(BinaryHeader::SIZE + payload.len());
    out.extend_from_slice(&header.to_bytes());
    out.extend_from_slice(&payload);
    out
}

/// Read and validate a binary field written by `write_binary_field`.
///
/// Returns `(header, data)` on success, or `Err` describing the problem.
pub fn read_binary_field(bytes: &[u8]) -> Result<(BinaryHeader, Vec<f64>), String> {
    if bytes.len() < BinaryHeader::SIZE {
        return Err("Data too short to contain a header".into());
    }
    let header = BinaryHeader::from_bytes(&bytes[..BinaryHeader::SIZE])?;

    if header.magic != MAGIC {
        return Err(format!("Bad magic: {:?}", header.magic));
    }
    if header.version != FORMAT_VERSION {
        return Err(format!("Unsupported version: {}", header.version));
    }

    let payload = &bytes[BinaryHeader::SIZE..];
    verify_integrity(payload, header.checksum)?;

    let data = bytes_to_f64_slice(payload)?;
    let expected = header.element_count();
    if data.len() != expected {
        return Err(format!(
            "Element count mismatch: header says {}, payload has {}",
            expected,
            data.len()
        ));
    }
    Ok((header, data))
}

// ── RLE compression ───────────────────────────────────────────────────────────

/// A run in the RLE encoding: `(value, count)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RleRun {
    /// The repeated value.
    pub value: f64,
    /// How many times it is repeated.
    pub count: usize,
}

/// Run-length encode a `f64` slice.
///
/// Adjacent values that are equal (bitwise) are collapsed into a single
/// `RleRun`. The encoding is lossless.
pub fn compress_rle(data: &[f64]) -> Vec<RleRun> {
    if data.is_empty() {
        return Vec::new();
    }
    let mut runs = Vec::new();
    let mut current_val = data[0];
    let mut count = 1usize;

    for &v in &data[1..] {
        if v.to_bits() == current_val.to_bits() {
            count += 1;
        } else {
            runs.push(RleRun {
                value: current_val,
                count,
            });
            current_val = v;
            count = 1;
        }
    }
    runs.push(RleRun {
        value: current_val,
        count,
    });
    runs
}

/// Decompress RLE runs back into a `f64` Vec.
pub fn decompress_rle(runs: &[RleRun]) -> Vec<f64> {
    let total: usize = runs.iter().map(|r| r.count).sum();
    let mut out = Vec::with_capacity(total);
    for run in runs {
        for _ in 0..run.count {
            out.push(run.value);
        }
    }
    out
}

// ── Checksum ──────────────────────────────────────────────────────────────────

/// Compute a CRC32 checksum of the given byte slice.
///
/// Uses the standard IEEE 802.3 polynomial (0xEDB88320, reflected).
pub fn checksum_crc32(data: &[u8]) -> u32 {
    // Build lookup table
    let table: [u32; 256] = {
        let mut t = [0u32; 256];
        for (i, entry) in t.iter_mut().enumerate() {
            let mut crc = i as u32;
            for _ in 0..8 {
                if crc & 1 != 0 {
                    crc = (crc >> 1) ^ 0xEDB8_8320;
                } else {
                    crc >>= 1;
                }
            }
            *entry = crc;
        }
        t
    };

    let mut crc = 0xFFFF_FFFFu32;
    for &byte in data {
        let idx = ((crc ^ byte as u32) & 0xFF) as usize;
        crc = (crc >> 8) ^ table[idx];
    }
    crc ^ 0xFFFF_FFFF
}

/// Verify that `data` matches the expected `checksum`.
///
/// Returns `Ok(())` on success or `Err` with a diagnostic message.
pub fn verify_integrity(data: &[u8], expected: u32) -> Result<(), String> {
    let actual = checksum_crc32(data);
    if actual == expected {
        Ok(())
    } else {
        Err(format!(
            "Checksum mismatch: expected 0x{:08X}, got 0x{:08X}",
            expected, actual
        ))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_data() -> Vec<f64> {
        (0..12).map(|i| i as f64 * 1.5).collect()
    }

    // Header serialisation
    #[test]
    fn test_header_roundtrip() {
        let header = BinaryHeader {
            magic: MAGIC,
            version: FORMAT_VERSION,
            data_type: FieldDataType::Float64,
            ndim: 1,
            dims: [12, 1, 1],
            checksum: 0xDEAD_BEEF,
        };
        let bytes = header.to_bytes();
        assert_eq!(bytes.len(), BinaryHeader::SIZE);
        let parsed = BinaryHeader::from_bytes(&bytes).unwrap();
        assert_eq!(parsed, header);
    }

    #[test]
    fn test_header_magic() {
        let header = BinaryHeader {
            magic: MAGIC,
            version: 1,
            data_type: FieldDataType::Float64,
            ndim: 1,
            dims: [1, 1, 1],
            checksum: 0,
        };
        let bytes = header.to_bytes();
        assert_eq!(&bytes[0..4], b"OXPF");
    }

    #[test]
    fn test_header_from_bytes_too_short() {
        assert!(BinaryHeader::from_bytes(&[0u8; 10]).is_err());
    }

    #[test]
    fn test_header_from_bytes_bad_data_type() {
        let mut bytes = [0u8; BinaryHeader::SIZE];
        bytes[6] = 99; // invalid data type
        assert!(BinaryHeader::from_bytes(&bytes).is_err());
    }

    #[test]
    fn test_header_element_count_1d() {
        let h = BinaryHeader {
            magic: MAGIC,
            version: 1,
            data_type: FieldDataType::Float64,
            ndim: 1,
            dims: [10, 1, 1],
            checksum: 0,
        };
        assert_eq!(h.element_count(), 10);
    }

    #[test]
    fn test_header_element_count_3d() {
        let h = BinaryHeader {
            magic: MAGIC,
            version: 1,
            data_type: FieldDataType::Float64,
            ndim: 3,
            dims: [4, 5, 6],
            checksum: 0,
        };
        assert_eq!(h.element_count(), 120);
    }

    // Write / read roundtrip
    #[test]
    fn test_write_read_roundtrip_1d() {
        let data = sample_data();
        let bytes = write_binary_field(&data, [12, 1, 1], 1);
        let (_hdr, recovered) = read_binary_field(&bytes).unwrap();
        assert_eq!(recovered.len(), data.len());
        for (a, b) in data.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn test_write_read_roundtrip_3d() {
        let data: Vec<f64> = (0..60).map(|i| i as f64).collect();
        let bytes = write_binary_field(&data, [3, 4, 5], 3);
        let (hdr, recovered) = read_binary_field(&bytes).unwrap();
        assert_eq!(hdr.dims, [3, 4, 5]);
        assert_eq!(recovered.len(), 60);
    }

    #[test]
    fn test_read_bad_magic() {
        let mut bytes = write_binary_field(&[1.0, 2.0], [2, 1, 1], 1);
        bytes[0] = b'X'; // corrupt magic
        assert!(read_binary_field(&bytes).is_err());
    }

    #[test]
    fn test_read_corrupted_payload() {
        let mut bytes = write_binary_field(&[1.0, 2.0], [2, 1, 1], 1);
        // Flip a byte in the payload
        let last = bytes.len() - 1;
        bytes[last] ^= 0xFF;
        assert!(read_binary_field(&bytes).is_err());
    }

    #[test]
    fn test_read_too_short() {
        assert!(read_binary_field(&[0u8; 5]).is_err());
    }

    #[test]
    fn test_write_read_empty() {
        let bytes = write_binary_field(&[], [0, 1, 1], 1);
        let (hdr, data) = read_binary_field(&bytes).unwrap();
        assert_eq!(hdr.dims[0], 0);
        assert!(data.is_empty());
    }

    // RLE compress / decompress
    #[test]
    fn test_rle_basic() {
        let data = vec![1.0, 1.0, 1.0, 2.0, 2.0, 3.0];
        let runs = compress_rle(&data);
        assert_eq!(runs.len(), 3);
        assert_eq!(runs[0].count, 3);
        assert_eq!(runs[1].count, 2);
        assert_eq!(runs[2].count, 1);
    }

    #[test]
    fn test_rle_empty() {
        assert!(compress_rle(&[]).is_empty());
    }

    #[test]
    fn test_rle_single_element() {
        let runs = compress_rle(&[42.0]);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].count, 1);
    }

    #[test]
    fn test_rle_no_repetition() {
        let data = vec![1.0, 2.0, 3.0];
        let runs = compress_rle(&data);
        assert_eq!(runs.len(), 3);
        for run in &runs {
            assert_eq!(run.count, 1);
        }
    }

    #[test]
    fn test_rle_roundtrip() {
        let data = vec![0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 3.0];
        let runs = compress_rle(&data);
        let recovered = decompress_rle(&runs);
        assert_eq!(recovered, data);
    }

    #[test]
    fn test_rle_all_same() {
        let data = vec![5.0; 100];
        let runs = compress_rle(&data);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].count, 100);
        let recovered = decompress_rle(&runs);
        assert_eq!(recovered, data);
    }

    #[test]
    fn test_rle_decompress_empty() {
        assert!(decompress_rle(&[]).is_empty());
    }

    // CRC32 checksum
    #[test]
    fn test_checksum_known_value() {
        // CRC32 of b"123456789" is 0xCBF43926 per IEEE specification
        let crc = checksum_crc32(b"123456789");
        assert_eq!(crc, 0xCBF4_3926);
    }

    #[test]
    fn test_checksum_empty() {
        let crc = checksum_crc32(b"");
        assert_eq!(crc, 0x0000_0000);
    }

    #[test]
    fn test_checksum_single_byte() {
        let c1 = checksum_crc32(b"A");
        let c2 = checksum_crc32(b"B");
        assert_ne!(c1, c2);
    }

    #[test]
    fn test_checksum_deterministic() {
        let data = b"physics_data_12345";
        assert_eq!(checksum_crc32(data), checksum_crc32(data));
    }

    // verify_integrity
    #[test]
    fn test_verify_integrity_ok() {
        let data = b"hello";
        let crc = checksum_crc32(data);
        assert!(verify_integrity(data, crc).is_ok());
    }

    #[test]
    fn test_verify_integrity_fail() {
        let data = b"hello";
        assert!(verify_integrity(data, 0xDEAD_BEEF).is_err());
    }

    // FieldDataType
    #[test]
    fn test_field_data_type_from_u8() {
        assert_eq!(FieldDataType::from_u8(0), Some(FieldDataType::Float64));
        assert_eq!(FieldDataType::from_u8(1), Some(FieldDataType::Float32));
        assert_eq!(FieldDataType::from_u8(2), Some(FieldDataType::Int32));
        assert_eq!(FieldDataType::from_u8(99), None);
    }

    // Integration: write, RLE-encode, compress, decompress, then verify
    #[test]
    fn test_write_read_with_rle_consistency() {
        let original = vec![0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
        let bytes = write_binary_field(&original, [10, 1, 1], 1);
        let (_hdr, recovered) = read_binary_field(&bytes).unwrap();

        // Apply RLE to the recovered data
        let runs = compress_rle(&recovered);
        let decompressed = decompress_rle(&runs);
        assert_eq!(decompressed, original);
    }
}
