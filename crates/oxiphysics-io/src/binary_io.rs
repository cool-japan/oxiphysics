// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Binary data I/O for OxiPhysics trajectory files.
//!
//! Provides a compact binary header, frame-level read/write of atomic
//! positions, little-endian byte-conversion helpers, and lossy quantized
//! compression for trajectory data.

use std::io::{Read, Write};

// ── Header ────────────────────────────────────────────────────────────────────

/// Binary file header for an OxiPhysics trajectory file.
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryHeader {
    /// Magic bytes identifying the file format (e.g. `b"OXIP"`).
    pub magic: [u8; 4],
    /// Format version number.
    pub version: u32,
    /// Number of atoms per frame.
    pub n_atoms: u32,
    /// Total number of frames in the file.
    pub n_frames: u32,
    /// Time step in picoseconds.
    pub dt: f64,
    /// Endianness flag: `0` = little-endian, `1` = big-endian.
    pub endian: u8,
}

/// Byte size of a serialised [`BinaryHeader`].
pub const HEADER_SIZE: usize = 4 + 4 + 4 + 4 + 8 + 1; // 25 bytes

/// Write a [`BinaryHeader`] to `path`, creating or truncating the file.
pub fn write_binary_header(path: &str, header: &BinaryHeader) -> Result<(), std::io::Error> {
    let mut file = std::fs::File::create(path)?;
    file.write_all(&header.magic)?;
    file.write_all(&header.version.to_le_bytes())?;
    file.write_all(&header.n_atoms.to_le_bytes())?;
    file.write_all(&header.n_frames.to_le_bytes())?;
    file.write_all(&f64_to_bytes_le(header.dt))?;
    file.write_all(&[header.endian])?;
    Ok(())
}

/// Read a [`BinaryHeader`] from `path`.
pub fn read_binary_header(path: &str) -> Result<BinaryHeader, std::io::Error> {
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)?;
    let mut buf4 = [0u8; 4];
    file.read_exact(&mut buf4)?;
    let version = u32::from_le_bytes(buf4);
    file.read_exact(&mut buf4)?;
    let n_atoms = u32::from_le_bytes(buf4);
    file.read_exact(&mut buf4)?;
    let n_frames = u32::from_le_bytes(buf4);
    let mut buf8 = [0u8; 8];
    file.read_exact(&mut buf8)?;
    let dt = f64_from_bytes_le(buf8);
    let mut endian_buf = [0u8; 1];
    file.read_exact(&mut endian_buf)?;
    Ok(BinaryHeader {
        magic,
        version,
        n_atoms,
        n_frames,
        dt,
        endian: endian_buf[0],
    })
}

// ── Frame I/O ─────────────────────────────────────────────────────────────────

/// Write one frame of atomic positions to an open file.
///
/// Each atom contributes 12 bytes (3 × f32 little-endian).
pub fn write_frame_binary(
    file: &mut std::fs::File,
    positions: &[[f32; 3]],
) -> Result<(), std::io::Error> {
    for pos in positions {
        file.write_all(&f32_to_bytes_le(pos[0]))?;
        file.write_all(&f32_to_bytes_le(pos[1]))?;
        file.write_all(&f32_to_bytes_le(pos[2]))?;
    }
    Ok(())
}

/// Read one frame of `n_atoms` atomic positions from an open file.
pub fn read_frame_binary(
    file: &mut std::fs::File,
    n_atoms: usize,
) -> Result<Vec<[f32; 3]>, std::io::Error> {
    let mut positions = Vec::with_capacity(n_atoms);
    let mut buf4 = [0u8; 4];
    for _ in 0..n_atoms {
        file.read_exact(&mut buf4)?;
        let x = f32_from_bytes_le(buf4);
        file.read_exact(&mut buf4)?;
        let y = f32_from_bytes_le(buf4);
        file.read_exact(&mut buf4)?;
        let z = f32_from_bytes_le(buf4);
        positions.push([x, y, z]);
    }
    Ok(positions)
}

// ── Byte-conversion helpers ───────────────────────────────────────────────────

/// Encode a `f64` value as 8 little-endian bytes.
pub fn f64_to_bytes_le(val: f64) -> [u8; 8] {
    val.to_bits().to_le_bytes()
}

/// Decode a `f64` value from 8 little-endian bytes.
pub fn f64_from_bytes_le(bytes: [u8; 8]) -> f64 {
    f64::from_bits(u64::from_le_bytes(bytes))
}

/// Encode a `f32` value as 4 little-endian bytes.
pub fn f32_to_bytes_le(val: f32) -> [u8; 4] {
    val.to_bits().to_le_bytes()
}

/// Decode a `f32` value from 4 little-endian bytes.
pub fn f32_from_bytes_le(bytes: [u8; 4]) -> f32 {
    f32::from_bits(u32::from_le_bytes(bytes))
}

// ── Quantized compression ─────────────────────────────────────────────────────

/// Lossy quantization of positions to `i16`.
///
/// Each coordinate is multiplied by `1.0 / scale` and clamped to
/// `[i16::MIN, i16::MAX]`.  Typical `scale` value: `0.001` for 1 pm
/// resolution with positions in Ångström.
pub fn compress_positions_quantized(positions: &[[f32; 3]], scale: f32) -> Vec<i16> {
    if scale == 0.0 {
        return vec![0i16; positions.len() * 3];
    }
    let inv = 1.0 / scale;
    let mut out = Vec::with_capacity(positions.len() * 3);
    for pos in positions {
        for &c in pos {
            let q = (c * inv).round().clamp(i16::MIN as f32, i16::MAX as f32) as i16;
            out.push(q);
        }
    }
    out
}

/// Reconstruct positions from quantized `i16` data.
///
/// The inverse of [`compress_positions_quantized`]: multiplies by `scale`.
pub fn decompress_positions_quantized(data: &[i16], scale: f32) -> Vec<[f32; 3]> {
    let n = data.len() / 3;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let x = data[i * 3] as f32 * scale;
        let y = data[i * 3 + 1] as f32 * scale;
        let z = data[i * 3 + 2] as f32 * scale;
        out.push([x, y, z]);
    }
    out
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_header() -> BinaryHeader {
        BinaryHeader {
            magic: *b"OXIP",
            version: 1,
            n_atoms: 100,
            n_frames: 50,
            dt: 0.002,
            endian: 0,
        }
    }

    #[test]
    fn test_f64_roundtrip() {
        let val = std::f64::consts::PI;
        let bytes = f64_to_bytes_le(val);
        let back = f64_from_bytes_le(bytes);
        assert!((back - val).abs() < 1e-15);
    }

    #[test]
    fn test_f64_zero() {
        let bytes = f64_to_bytes_le(0.0);
        let back = f64_from_bytes_le(bytes);
        assert!((back).abs() < 1e-15);
    }

    #[test]
    fn test_f64_negative() {
        let val = -1.23456789e-10;
        let bytes = f64_to_bytes_le(val);
        let back = f64_from_bytes_le(bytes);
        assert!((back - val).abs() < 1e-20);
    }

    #[test]
    fn test_f32_roundtrip() {
        let val = 3.15625_f32;
        let bytes = f32_to_bytes_le(val);
        let back = f32_from_bytes_le(bytes);
        assert!((back - val).abs() < 1e-6);
    }

    #[test]
    fn test_f32_negative() {
        let val = -0.001_f32;
        let bytes = f32_to_bytes_le(val);
        let back = f32_from_bytes_le(bytes);
        assert!((back - val).abs() < 1e-7);
    }

    #[test]
    fn test_write_read_header_roundtrip() {
        let path = "/tmp/test_oxiphysics_header.bin";
        let hdr = sample_header();
        write_binary_header(path, &hdr).unwrap();
        let hdr2 = read_binary_header(path).unwrap();
        assert_eq!(hdr2.magic, *b"OXIP");
        assert_eq!(hdr2.version, 1);
        assert_eq!(hdr2.n_atoms, 100);
        assert_eq!(hdr2.n_frames, 50);
        assert!((hdr2.dt - 0.002).abs() < 1e-12);
        assert_eq!(hdr2.endian, 0);
    }

    #[test]
    fn test_header_magic() {
        let path = "/tmp/test_oxiphysics_header_magic.bin";
        let hdr = BinaryHeader {
            magic: *b"TEST",
            ..sample_header()
        };
        write_binary_header(path, &hdr).unwrap();
        let hdr2 = read_binary_header(path).unwrap();
        assert_eq!(hdr2.magic, *b"TEST");
    }

    #[test]
    fn test_header_size() {
        assert_eq!(HEADER_SIZE, 25);
    }

    #[test]
    fn test_write_read_frame_roundtrip() {
        let path = "/tmp/test_oxiphysics_frame.bin";
        let positions = vec![[1.0_f32, 2.0, 3.0], [4.0, 5.0, 6.0], [-1.0, 0.0, 0.5]];
        {
            let mut file = std::fs::File::create(path).unwrap();
            write_frame_binary(&mut file, &positions).unwrap();
        }
        {
            let mut file = std::fs::File::open(path).unwrap();
            let back = read_frame_binary(&mut file, 3).unwrap();
            assert_eq!(back.len(), 3);
            assert!((back[0][0] - 1.0).abs() < 1e-6);
            assert!((back[2][2] - 0.5).abs() < 1e-6);
        }
    }

    #[test]
    fn test_write_read_frame_single_atom() {
        let path = "/tmp/test_oxiphysics_frame_single.bin";
        let positions = vec![[0.1_f32, 0.2, 0.3]];
        {
            let mut file = std::fs::File::create(path).unwrap();
            write_frame_binary(&mut file, &positions).unwrap();
        }
        {
            let mut file = std::fs::File::open(path).unwrap();
            let back = read_frame_binary(&mut file, 1).unwrap();
            assert!((back[0][1] - 0.2).abs() < 1e-6);
        }
    }

    #[test]
    fn test_write_read_empty_frame() {
        let path = "/tmp/test_oxiphysics_frame_empty.bin";
        {
            let mut file = std::fs::File::create(path).unwrap();
            write_frame_binary(&mut file, &[]).unwrap();
        }
        {
            let mut file = std::fs::File::open(path).unwrap();
            let back = read_frame_binary(&mut file, 0).unwrap();
            assert!(back.is_empty());
        }
    }

    #[test]
    fn test_compress_decompress_positions() {
        let positions = vec![[1.0_f32, 2.0, -1.5], [0.5, 0.0, 3.0]];
        let scale = 0.001;
        let compressed = compress_positions_quantized(&positions, scale);
        let decompressed = decompress_positions_quantized(&compressed, scale);
        assert_eq!(decompressed.len(), 2);
        assert!((decompressed[0][0] - 1.0).abs() < scale * 2.0);
        assert!((decompressed[1][2] - 3.0).abs() < scale * 2.0);
    }

    #[test]
    fn test_compress_empty() {
        let compressed = compress_positions_quantized(&[], 0.001);
        assert!(compressed.is_empty());
    }

    #[test]
    fn test_decompress_empty() {
        let decompressed = decompress_positions_quantized(&[], 0.001);
        assert!(decompressed.is_empty());
    }

    #[test]
    fn test_compress_zero_scale() {
        // scale=0 returns all zeros
        let positions = vec![[1.0_f32, 2.0, 3.0]];
        let compressed = compress_positions_quantized(&positions, 0.0);
        assert_eq!(compressed, vec![0i16, 0, 0]);
    }

    #[test]
    fn test_compress_large_positive() {
        // Value that would overflow i16: clamps to i16::MAX
        let positions = vec![[1e10_f32, 0.0, 0.0]];
        let compressed = compress_positions_quantized(&positions, 0.001);
        assert_eq!(compressed[0], i16::MAX);
    }

    #[test]
    fn test_compress_large_negative() {
        let positions = vec![[-1e10_f32, 0.0, 0.0]];
        let compressed = compress_positions_quantized(&positions, 0.001);
        assert_eq!(compressed[0], i16::MIN);
    }

    #[test]
    fn test_compress_count() {
        let positions = vec![[0.0_f32; 3]; 5];
        let compressed = compress_positions_quantized(&positions, 0.01);
        assert_eq!(compressed.len(), 15);
    }

    #[test]
    fn test_decompress_count() {
        let data = vec![0i16; 12];
        let out = decompress_positions_quantized(&data, 0.01);
        assert_eq!(out.len(), 4);
    }

    #[test]
    fn test_multiple_frames_sequential() {
        let path = "/tmp/test_oxiphysics_multiframe.bin";
        let frame1 = vec![[1.0_f32, 0.0, 0.0]];
        let frame2 = vec![[2.0_f32, 0.0, 0.0]];
        {
            let mut file = std::fs::File::create(path).unwrap();
            write_frame_binary(&mut file, &frame1).unwrap();
            write_frame_binary(&mut file, &frame2).unwrap();
        }
        {
            let mut file = std::fs::File::open(path).unwrap();
            let f1 = read_frame_binary(&mut file, 1).unwrap();
            let f2 = read_frame_binary(&mut file, 1).unwrap();
            assert!((f1[0][0] - 1.0).abs() < 1e-6);
            assert!((f2[0][0] - 2.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_header_dt_precision() {
        let path = "/tmp/test_oxiphysics_header_dt.bin";
        let hdr = BinaryHeader {
            dt: 1.23456789012345e-4,
            ..sample_header()
        };
        write_binary_header(path, &hdr).unwrap();
        let hdr2 = read_binary_header(path).unwrap();
        assert!((hdr2.dt - 1.23456789012345e-4).abs() < 1e-18);
    }

    #[test]
    fn test_header_large_n_atoms() {
        let path = "/tmp/test_oxiphysics_header_large.bin";
        let hdr = BinaryHeader {
            n_atoms: 1_000_000,
            ..sample_header()
        };
        write_binary_header(path, &hdr).unwrap();
        let hdr2 = read_binary_header(path).unwrap();
        assert_eq!(hdr2.n_atoms, 1_000_000);
    }

    #[test]
    fn test_quantize_negative_values() {
        let positions = vec![[-0.5_f32, -1.0, -2.0]];
        let scale = 0.001;
        let compressed = compress_positions_quantized(&positions, scale);
        let back = decompress_positions_quantized(&compressed, scale);
        assert!((back[0][0] - (-0.5)).abs() < scale * 2.0);
    }

    #[test]
    fn test_f32_bytes_all_zeros() {
        let bytes = f32_to_bytes_le(0.0);
        assert_eq!(bytes, [0, 0, 0, 0]);
    }

    #[test]
    fn test_f64_bytes_all_zeros() {
        let bytes = f64_to_bytes_le(0.0);
        assert_eq!(bytes, [0, 0, 0, 0, 0, 0, 0, 0]);
    }
}
