//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::io::{self, Write};

/// Magic number written at the start of every checkpoint file.
pub(super) const MAGIC: u32 = 0x4F58_4350;
/// Current format version.
pub(super) const FORMAT_VERSION: u32 = 1;
/// Tag byte identifying position records in a checkpoint file.
pub(super) const TAG_POSITIONS: u8 = 0x01;
/// Tag byte identifying velocity records in a checkpoint file.
pub(super) const TAG_VELOCITIES: u8 = 0x02;
/// Tag byte identifying named scalar array records.
pub(super) const TAG_SCALARS: u8 = 0x03;
/// Tag byte identifying named integer array records.
pub(super) const TAG_INTEGERS: u8 = 0x04;
/// Tag byte identifying the footer record.
pub(super) const TAG_FOOTER: u8 = 0xFF;
/// Size in bytes of the footer record: the [`TAG_FOOTER`] byte followed by the
/// little-endian CRC-32 of every byte that precedes it.
pub(super) const FOOTER_LEN: usize = 5;
/// Size in bytes of one stored three-component `f64` vector.
pub(super) const VEC3_BYTES: usize = 24;
/// Size in bytes of one stored `f64`.
const F64_BYTES: usize = 8;
/// Size in bytes of one stored `i32`.
const I32_BYTES: usize = 4;

fn truncated(what: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, what)
}

pub(super) fn corrupt(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Borrow `len` bytes of `data` starting at `*cursor` and advance the cursor
/// past them.
///
/// The cursor arithmetic is overflow-checked, so a hostile length prefix
/// yields an `UnexpectedEof` error (and leaves the cursor untouched) instead
/// of a wrapped cursor or a panic.
pub(super) fn take_bytes<'a>(
    data: &'a [u8],
    cursor: &mut usize,
    len: usize,
    what: &'static str,
) -> io::Result<&'a [u8]> {
    let end = cursor.checked_add(len).ok_or_else(|| truncated(what))?;
    let bytes = data.get(*cursor..end).ok_or_else(|| truncated(what))?;
    *cursor = end;
    Ok(bytes)
}
pub(super) fn read_u32(data: &[u8], cursor: &mut usize) -> io::Result<u32> {
    let bytes = take_bytes(data, cursor, 4, "read_u32")?;
    let array = <[u8; 4]>::try_from(bytes).map_err(|_| truncated("read_u32"))?;
    Ok(u32::from_le_bytes(array))
}
pub(super) fn read_u64(data: &[u8], cursor: &mut usize) -> io::Result<u64> {
    let bytes = take_bytes(data, cursor, 8, "read_u64")?;
    let array = <[u8; 8]>::try_from(bytes).map_err(|_| truncated("read_u64"))?;
    Ok(u64::from_le_bytes(array))
}
pub(super) fn read_f64(data: &[u8], cursor: &mut usize) -> io::Result<f64> {
    let bytes = take_bytes(data, cursor, 8, "read_f64")?;
    let array = <[u8; 8]>::try_from(bytes).map_err(|_| truncated("read_f64"))?;
    Ok(f64::from_le_bytes(array))
}

/// Table for the byte-wise CRC-32/ISO-HDLC algorithm (reflected polynomial
/// `0xEDB88320`).
const CRC32_TABLE: [u32; 256] = build_crc32_table();

const fn build_crc32_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    let mut index = 0usize;
    while index < 256 {
        let mut crc = index as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xEDB8_8320
            } else {
                crc >> 1
            };
            bit += 1;
        }
        table[index] = crc;
        index += 1;
    }
    table
}

/// Compute the CRC-32/ISO-HDLC checksum of `data`.
///
/// Uses the reflected polynomial `0xEDB88320` with an initial value and an
/// XOR-out of `0xFFFFFFFF` (the checksum used by zlib, PNG and Ethernet), so
/// `compute_checksum(b"123456789") == 0xCBF43926`.
pub fn compute_checksum(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in data {
        let index = ((crc ^ u32::from(byte)) & 0xFF) as usize;
        crc = (crc >> 8) ^ CRC32_TABLE[index];
    }
    crc ^ 0xFFFF_FFFF
}
pub(super) fn write_name(f: &mut impl Write, name: &str) -> io::Result<()> {
    let bytes = name.as_bytes();
    f.write_all(&(bytes.len() as u32).to_le_bytes())?;
    f.write_all(bytes)
}
pub(super) fn read_name(data: &[u8], cursor: &mut usize) -> io::Result<String> {
    let len = read_u32(data, cursor)? as usize;
    let bytes = take_bytes(data, cursor, len, "name truncated")?;
    std::str::from_utf8(bytes)
        .map(str::to_owned)
        .map_err(|e| corrupt(format!("UTF-8: {e}")))
}

/// Check that `data` begins with the checkpoint magic number.
///
/// Data too short to hold the magic number is a truncated file and is
/// reported like a missing footer.
pub(super) fn check_magic(data: &[u8]) -> io::Result<()> {
    let mut cursor = 0usize;
    let magic = read_u32(data, &mut cursor)
        .map_err(|_| corrupt("checkpoint not finalized or truncated"))?;
    if magic != MAGIC {
        return Err(corrupt("bad magic number"));
    }
    Ok(())
}

/// Verify the footer that `CheckpointWriter::finalize` appends and return the
/// bytes it protects (everything before the footer).
///
/// The footer is the last [`FOOTER_LEN`] bytes of a finalised checkpoint: the
/// [`TAG_FOOTER`] byte and the CRC-32 of the preceding bytes.  A missing footer
/// (an unfinished or truncated file) and a checksum that does not match (bit
/// rot, a torn write) are both reported as `InvalidData`.
pub(super) fn verify_footer(data: &[u8]) -> io::Result<&[u8]> {
    let not_finalized = || corrupt("checkpoint not finalized or truncated");
    let body_len = data
        .len()
        .checked_sub(FOOTER_LEN)
        .ok_or_else(not_finalized)?;
    let (body, footer) = data.split_at(body_len);
    let [tag, c0, c1, c2, c3] =
        <[u8; FOOTER_LEN]>::try_from(footer).map_err(|_| not_finalized())?;
    if tag != TAG_FOOTER {
        return Err(not_finalized());
    }
    let stored = u32::from_le_bytes([c0, c1, c2, c3]);
    let computed = compute_checksum(body);
    if stored != computed {
        return Err(corrupt(format!(
            "checkpoint checksum mismatch (stored {stored:#010x}, computed {computed:#010x})"
        )));
    }
    Ok(body)
}

/// Validate the file header at the start of `data` and split it up.
///
/// Returns the serialized metadata record and the offset of the first block
/// that follows it.  The magic number is checked and the metadata length is
/// bounded by the data actually present.
pub(super) fn split_header(data: &[u8]) -> io::Result<(&[u8], usize)> {
    let mut cursor = 0usize;
    if read_u32(data, &mut cursor)? != MAGIC {
        return Err(corrupt("bad magic number"));
    }
    let _version = read_u32(data, &mut cursor)?;
    let meta_len = read_u32(data, &mut cursor)? as usize;
    let meta = take_bytes(data, &mut cursor, meta_len, "metadata block truncated")?;
    Ok((meta, cursor))
}

/// Validate a block that claims to hold `count` elements of `elem_size` bytes
/// starting at `cursor` and return `(count, end)`: the element count as a
/// `usize` and the offset one past the block.
///
/// Every step is checked (`u64` to `usize`, `count * elem_size`,
/// `cursor + bytes`) and the block must lie entirely inside `data_len`.  A
/// corrupt count therefore produces an `InvalidData` error rather than a
/// wrapped cursor, an endless loop or an oversized allocation, and a caller
/// may safely reserve `count` elements.
pub(super) fn checked_block_end(
    data_len: usize,
    cursor: usize,
    count: u64,
    elem_size: usize,
) -> io::Result<(usize, usize)> {
    let count = usize::try_from(count)
        .map_err(|_| corrupt("block element count does not fit in memory"))?;
    let byte_len = count
        .checked_mul(elem_size)
        .ok_or_else(|| corrupt("block size overflows"))?;
    let end = cursor
        .checked_add(byte_len)
        .ok_or_else(|| corrupt("block extends past end of file"))?;
    if end > data_len {
        return Err(corrupt(format!(
            "block of {count} elements extends past end of file"
        )));
    }
    Ok((count, end))
}

/// One record of the block stream that follows the checkpoint header.
///
/// The payloads borrow the validated bytes of the file: positions and
/// velocities are packed `[f64; 3]` triples, scalar arrays packed `f64`
/// values and integer arrays packed `i32` values, all little-endian.
pub(super) enum Block<'a> {
    /// Particle positions.
    Positions(&'a [u8]),
    /// Particle velocities.
    Velocities(&'a [u8]),
    /// A named `f64` array.
    Scalars {
        /// Array name.
        name: String,
        /// Packed values.
        payload: &'a [u8],
    },
    /// A named `i32` array.
    Integers {
        /// Array name.
        name: String,
        /// Packed values.
        payload: &'a [u8],
    },
}

/// Read the record at `*cursor` and advance past it.
///
/// Returns `Ok(None)` at the footer or at the end of `data`.  Every count is
/// validated with [`checked_block_end`] before its payload is borrowed, so
/// all readers share one set of overflow and bounds guards.
pub(super) fn next_block<'a>(data: &'a [u8], cursor: &mut usize) -> io::Result<Option<Block<'a>>> {
    let Some(&tag) = data.get(*cursor) else {
        return Ok(None);
    };
    if tag == TAG_FOOTER {
        return Ok(None);
    }
    let mut at = *cursor + 1;
    let block = match tag {
        TAG_POSITIONS | TAG_VELOCITIES => {
            let count = read_u64(data, &mut at)?;
            let (_, end) = checked_block_end(data.len(), at, count, VEC3_BYTES)?;
            let payload = data
                .get(at..end)
                .ok_or_else(|| truncated("vector block truncated"))?;
            at = end;
            if tag == TAG_POSITIONS {
                Block::Positions(payload)
            } else {
                Block::Velocities(payload)
            }
        }
        TAG_SCALARS | TAG_INTEGERS => {
            let name = read_name(data, &mut at)?;
            let count = read_u64(data, &mut at)?;
            let elem_size = if tag == TAG_SCALARS {
                F64_BYTES
            } else {
                I32_BYTES
            };
            let (_, end) = checked_block_end(data.len(), at, count, elem_size)?;
            let payload = data
                .get(at..end)
                .ok_or_else(|| truncated("array block truncated"))?;
            at = end;
            if tag == TAG_SCALARS {
                Block::Scalars { name, payload }
            } else {
                Block::Integers { name, payload }
            }
        }
        other => return Err(corrupt(format!("unknown tag 0x{other:02X}"))),
    };
    *cursor = at;
    Ok(Some(block))
}

/// Decode a packed run of little-endian `[f64; 3]` triples.
pub(super) fn decode_vec3(payload: &[u8]) -> io::Result<Vec<[f64; 3]>> {
    let mut out = Vec::with_capacity(payload.len() / VEC3_BYTES);
    let mut cursor = 0usize;
    while cursor < payload.len() {
        let x = read_f64(payload, &mut cursor)?;
        let y = read_f64(payload, &mut cursor)?;
        let z = read_f64(payload, &mut cursor)?;
        out.push([x, y, z]);
    }
    Ok(out)
}

/// Decode a packed run of little-endian `f64` values.
pub(super) fn decode_f64_array(payload: &[u8]) -> io::Result<Vec<f64>> {
    let mut out = Vec::with_capacity(payload.len() / F64_BYTES);
    let mut cursor = 0usize;
    while cursor < payload.len() {
        out.push(read_f64(payload, &mut cursor)?);
    }
    Ok(out)
}

/// Decode a packed run of little-endian `i32` values.
pub(super) fn decode_i32_array(payload: &[u8]) -> io::Result<Vec<i32>> {
    let mut out = Vec::with_capacity(payload.len() / I32_BYTES);
    let mut cursor = 0usize;
    while cursor < payload.len() {
        let bytes = take_bytes(payload, &mut cursor, I32_BYTES, "read_i32")?;
        let array = <[u8; 4]>::try_from(bytes).map_err(|_| truncated("read_i32"))?;
        out.push(i32::from_le_bytes(array));
    }
    Ok(out)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::checkpoint_io::types::*;
    use std::collections::HashMap;
    use std::fs;
    use std::path::{Path, PathBuf};
    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "oxiphysics_ckpt_test_{}_{name}.bin",
            std::process::id()
        ));
        p
    }
    fn sample_meta() -> CheckpointMetadata {
        CheckpointMetadata::new(42, 0.001, 100, [0, 1, 0], "2026-01-01T00:00:00Z")
    }
    #[test]
    fn checksum_empty() {
        let c = compute_checksum(&[]);
        assert_eq!(c, 0x0000_0000);
    }
    #[test]
    fn checksum_deterministic() {
        let data = b"hello oxiphysics";
        assert_eq!(compute_checksum(data), compute_checksum(data));
    }
    #[test]
    fn checksum_differs_for_different_data() {
        let c1 = compute_checksum(b"aaa");
        let c2 = compute_checksum(b"bbb");
        assert_ne!(c1, c2);
    }
    #[test]
    fn checksum_single_byte() {
        let _ = compute_checksum(&[0xFF]);
    }
    #[test]
    fn metadata_roundtrip() {
        let meta = sample_meta();
        let bytes = meta.to_bytes();
        let restored = CheckpointMetadata::from_bytes(&bytes).unwrap();
        assert_eq!(meta, restored);
    }
    #[test]
    fn metadata_step_preserved() {
        let meta = CheckpointMetadata::new(9999, 3.125, 50, [1, 2, 3], "ts");
        let bytes = meta.to_bytes();
        let r = CheckpointMetadata::from_bytes(&bytes).unwrap();
        assert_eq!(r.step, 9999);
    }
    #[test]
    fn metadata_time_preserved() {
        let meta = CheckpointMetadata::new(0, 1.23456789, 10, [0, 0, 1], "ts");
        let bytes = meta.to_bytes();
        let r = CheckpointMetadata::from_bytes(&bytes).unwrap();
        assert!((r.time - 1.23456789).abs() < 1e-15);
    }
    #[test]
    fn metadata_version_preserved() {
        let meta = CheckpointMetadata::new(1, 0.0, 0, [3, 2, 1], "ts");
        let bytes = meta.to_bytes();
        let r = CheckpointMetadata::from_bytes(&bytes).unwrap();
        assert_eq!(r.crate_version, [3, 2, 1]);
    }
    #[test]
    fn metadata_created_at_preserved() {
        let meta = CheckpointMetadata::new(0, 0.0, 0, [0; 3], "hello world");
        let bytes = meta.to_bytes();
        let r = CheckpointMetadata::from_bytes(&bytes).unwrap();
        assert_eq!(r.created_at, "hello world");
    }
    #[test]
    fn write_read_header_only() {
        let path = tmp_path("header");
        let writer = CheckpointWriter::new(&path);
        let meta = sample_meta();
        writer.write_header(&meta).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let meta2 = reader.read_metadata().unwrap();
        assert_eq!(meta, meta2);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn write_read_positions() {
        let path = tmp_path("positions");
        let positions: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_positions(&positions).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let p2 = reader.read_positions().unwrap();
        assert_eq!(p2.len(), 2);
        assert!((p2[0][0] - 1.0).abs() < 1e-15);
        assert!((p2[1][2] - 6.0).abs() < 1e-15);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn write_read_velocities() {
        let path = tmp_path("velocities");
        let vel: Vec<[f64; 3]> = vec![[0.1, 0.2, 0.3]];
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_velocities(&vel).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let v2 = reader.read_velocities().unwrap();
        assert_eq!(v2.len(), 1);
        assert!((v2[0][1] - 0.2).abs() < 1e-15);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn write_read_scalars() {
        let path = tmp_path("scalars");
        let data = vec![1.1, 2.2, 3.3, 4.4];
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_scalars("temperature", &data).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let d2 = reader.read_scalars("temperature").unwrap();
        assert_eq!(d2.len(), 4);
        assert!((d2[2] - 3.3).abs() < 1e-14);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn write_read_integers() {
        let path = tmp_path("integers");
        let data = vec![10i32, 20, 30];
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_integers("types", &data).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let positions = reader.read_positions().unwrap();
        assert!(positions.is_empty());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn read_scalars_not_found_returns_error() {
        let path = tmp_path("scalars_missing");
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let result = reader.read_scalars("nonexistent");
        assert!(result.is_err());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn bad_magic_returns_error() {
        let path = tmp_path("bad_magic");
        fs::write(&path, [0xDE, 0xAD, 0xBE, 0xEF, 0x01, 0x00, 0x00, 0x00]).unwrap();
        let reader = CheckpointReader::new(&path);
        assert!(reader.read_metadata().is_err());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn restart_file_roundtrip() {
        let path = tmp_path("restart");
        let meta = sample_meta();
        let pos = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let vel = vec![[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]];
        let mut scalars = HashMap::new();
        scalars.insert("pressure".to_string(), vec![1.0, 2.0]);
        let rf = RestartFile::new(meta.clone(), pos.clone(), vel.clone(), scalars.clone());
        rf.save(&path).unwrap();
        let rf2 = RestartFile::load(&path).unwrap();
        assert_eq!(rf2.meta, meta);
        assert_eq!(rf2.positions.len(), 2);
        assert!((rf2.positions[0][0] - 1.0).abs() < 1e-15);
        assert_eq!(rf2.velocities.len(), 2);
        assert!(rf2.scalars.contains_key("pressure"));
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn restart_file_empty_particles() {
        let path = tmp_path("restart_empty");
        let meta = CheckpointMetadata::new(0, 0.0, 0, [0; 3], "ts");
        let rf = RestartFile::new(meta, vec![], vec![], HashMap::new());
        rf.save(&path).unwrap();
        let rf2 = RestartFile::load(&path).unwrap();
        assert!(rf2.positions.is_empty());
        assert!(rf2.velocities.is_empty());
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn restart_file_multiple_scalars() {
        let path = tmp_path("restart_multi");
        let meta = sample_meta();
        let mut scalars = HashMap::new();
        scalars.insert("a".to_string(), vec![1.0]);
        scalars.insert("b".to_string(), vec![2.0, 3.0]);
        let rf = RestartFile::new(meta, vec![], vec![], scalars);
        rf.save(&path).unwrap();
        let rf2 = RestartFile::load(&path).unwrap();
        assert!(rf2.scalars.contains_key("a"));
        assert!(rf2.scalars.contains_key("b"));
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn should_checkpoint_at_interval() {
        let mgr = CheckpointManager::new("/tmp", 5, 100);
        assert!(mgr.should_checkpoint(0));
        assert!(mgr.should_checkpoint(100));
        assert!(mgr.should_checkpoint(200));
        assert!(!mgr.should_checkpoint(1));
        assert!(!mgr.should_checkpoint(99));
    }
    #[test]
    fn should_checkpoint_zero_interval() {
        let tmpdir = std::env::temp_dir();
        let mgr = CheckpointManager::new(&tmpdir, 5, 0);
        assert!(!mgr.should_checkpoint(0));
        assert!(!mgr.should_checkpoint(100));
    }
    #[test]
    fn checkpoint_path_format() {
        let tmpdir = std::env::temp_dir();
        let mgr = CheckpointManager::new(tmpdir.join("sim"), 3, 50);
        let p = mgr.checkpoint_path(42);
        assert!(p.to_string_lossy().contains("checkpoint_0000000042.bin"));
    }
    #[test]
    fn list_checkpoints_empty_dir() {
        let dir = std::env::temp_dir().join("oxichk_empty_dir_test");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 5, 10);
        let list = mgr.list_checkpoints();
        assert!(list.is_empty());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn list_checkpoints_finds_files() {
        let dir = std::env::temp_dir().join("oxichk_list_test");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 5, 10);
        for step in [0u64, 10, 20] {
            let p = mgr.checkpoint_path(step);
            fs::write(&p, b"dummy").unwrap();
        }
        let list = mgr.list_checkpoints();
        assert_eq!(list.len(), 3);
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn latest_checkpoint_returns_last() {
        let dir = std::env::temp_dir().join("oxichk_latest_test");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 5, 10);
        for step in [0u64, 10, 20] {
            fs::write(mgr.checkpoint_path(step), b"x").unwrap();
        }
        let latest = mgr.latest_checkpoint().unwrap();
        assert!(latest.to_string_lossy().contains("0000000020"));
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn latest_checkpoint_empty_dir() {
        let dir = std::env::temp_dir().join("oxichk_latest_empty");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 5, 10);
        assert!(mgr.latest_checkpoint().is_none());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn prune_keeps_max_checkpoints() {
        let dir = std::env::temp_dir().join("oxichk_prune_test");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 3, 10);
        for step in [0u64, 10, 20, 30, 40] {
            fs::write(mgr.checkpoint_path(step), b"x").unwrap();
        }
        mgr.prune_old_checkpoints().unwrap();
        let list = mgr.list_checkpoints();
        assert_eq!(list.len(), 3);
        assert!(
            list.last()
                .unwrap()
                .to_string_lossy()
                .contains("0000000040")
        );
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn prune_no_op_when_under_limit() {
        let dir = std::env::temp_dir().join("oxichk_prune_noop");
        fs::create_dir_all(&dir).unwrap();
        let mgr = CheckpointManager::new(&dir, 10, 10);
        for step in [0u64, 10] {
            fs::write(mgr.checkpoint_path(step), b"x").unwrap();
        }
        mgr.prune_old_checkpoints().unwrap();
        assert_eq!(mgr.list_checkpoints().len(), 2);
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn writer_with_compress_flag() {
        let path = std::env::temp_dir().join("x.bin");
        let w = CheckpointWriter::new(path.to_str().unwrap_or("")).with_compress(true);
        assert!(w.compress);
    }
    #[test]
    fn full_write_read_pipeline() {
        let path = tmp_path("full_pipeline");
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let vel: Vec<[f64; 3]> = (0..10).map(|i| [0.0, i as f64, 0.0]).collect();
        let writer = CheckpointWriter::new(&path);
        let meta = CheckpointMetadata::new(100, 0.01, 10, [0, 1, 0], "2026-01-01");
        writer.write_header(&meta).unwrap();
        writer.write_positions(&pos).unwrap();
        writer.write_velocities(&vel).unwrap();
        writer.write_scalars("ke", &[0.5, 0.6, 0.7]).unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        let m2 = reader.read_metadata().unwrap();
        assert_eq!(m2.step, 100);
        let p2 = reader.read_positions().unwrap();
        assert_eq!(p2.len(), 10);
        let v2 = reader.read_velocities().unwrap();
        assert_eq!(v2.len(), 10);
        let ke2 = reader.read_scalars("ke").unwrap();
        assert_eq!(ke2.len(), 3);
        assert!((ke2[1] - 0.6).abs() < 1e-14);
        let _ = fs::remove_file(&path);
    }
    #[test]
    fn checkpoint_new_sets_checksum() {
        let state = vec![1u8, 2, 3, 4, 5];
        let ckpt = Checkpoint::new(1, 1000, 42, 0.042, state.clone());
        let expected = compute_checksum(&state);
        assert_eq!(ckpt.checksum, expected);
    }
    #[test]
    fn checkpoint_verify_valid() {
        let ckpt = Checkpoint::new(1, 999, 10, 1.0_f64, vec![0xAB, 0xCD]);
        assert!(ckpt.verify());
    }
    #[test]
    fn checkpoint_verify_detects_corruption() {
        let mut ckpt = Checkpoint::new(1, 999, 10, 1.0_f64, vec![1, 2, 3]);
        ckpt.state[0] ^= 0xFF;
        assert!(!ckpt.verify());
    }
    #[test]
    fn checkpoint_roundtrip_bytes() {
        let state: Vec<u8> = (0u8..20).collect();
        let ckpt = Checkpoint::new(2, 12345, 100, 3.125_f64, state);
        let bytes = ckpt.to_bytes();
        let ckpt2 = Checkpoint::from_bytes(&bytes).unwrap();
        assert_eq!(ckpt, ckpt2);
    }
    #[test]
    fn checkpoint_roundtrip_preserves_version() {
        let ckpt = Checkpoint::new(42, 0, 0, 0.0_f64, vec![]);
        let bytes = ckpt.to_bytes();
        let ckpt2 = Checkpoint::from_bytes(&bytes).unwrap();
        assert_eq!(ckpt2.version, 42);
    }
    #[test]
    fn checkpoint_roundtrip_preserves_sim_time() {
        let ckpt = Checkpoint::new(1, 0, 0, 9.99_f64, vec![]);
        let bytes = ckpt.to_bytes();
        let ckpt2 = Checkpoint::from_bytes(&bytes).unwrap();
        assert!((ckpt2.sim_time - 9.99_f64).abs() < 1e-14);
    }
    #[test]
    fn checkpoint_empty_state_roundtrip() {
        let ckpt = Checkpoint::new(1, 0, 0, 0.0_f64, vec![]);
        let bytes = ckpt.to_bytes();
        let ckpt2 = Checkpoint::from_bytes(&bytes).unwrap();
        assert!(ckpt2.state.is_empty());
    }
    #[test]
    fn checkpoint_compute_checksum_updates_field() {
        let mut ckpt = Checkpoint::new(1, 0, 0, 0.0_f64, vec![1, 2, 3]);
        ckpt.state.push(4);
        ckpt.compute_checksum();
        assert!(ckpt.verify());
    }
    #[test]
    fn diff_identical_states_has_no_edits() {
        let state = vec![0u8, 1, 2, 3];
        let diff = CheckpointDiff::compute(0, &state, 1, &state);
        assert_eq!(diff.diff_size(), 0);
    }
    #[test]
    fn diff_apply_reconstructs_target() {
        let base = vec![0u8, 0, 0, 0];
        let target = vec![0u8, 1, 0, 2];
        let diff = CheckpointDiff::compute(0, &base, 1, &target);
        let reconstructed = diff.apply(&base).unwrap();
        assert_eq!(reconstructed, target);
    }
    #[test]
    fn diff_detects_all_changes() {
        let base = vec![0u8; 8];
        let target = vec![1u8; 8];
        let diff = CheckpointDiff::compute(0, &base, 1, &target);
        assert_eq!(diff.diff_size(), 8);
    }
    #[test]
    fn diff_change_ratio_all_changed() {
        let base = vec![0u8; 10];
        let target = vec![1u8; 10];
        let diff = CheckpointDiff::compute(0, &base, 1, &target);
        let ratio = diff.change_ratio(10);
        assert!((ratio - 1.0_f64).abs() < 1e-10);
    }
    #[test]
    fn diff_change_ratio_none_changed() {
        let state = vec![5u8; 10];
        let diff = CheckpointDiff::compute(0, &state, 1, &state);
        let ratio = diff.change_ratio(10);
        assert!(ratio.abs() < 1e-10);
    }
    #[test]
    fn diff_partial_change() {
        let base = vec![0u8, 1, 2, 3, 4];
        let mut target = base.clone();
        target[2] = 99;
        let diff = CheckpointDiff::compute(0, &base, 1, &target);
        assert_eq!(diff.diff_size(), 1);
        let out = diff.apply(&base).unwrap();
        assert_eq!(out[2], 99);
    }
    #[test]
    fn diff_step_numbers_preserved() {
        let diff = CheckpointDiff::compute(10, &[0], 20, &[1]);
        assert_eq!(diff.base_step, 10);
        assert_eq!(diff.target_step, 20);
    }
    #[test]
    fn compressor_roundtrip_empty() {
        let c = CheckpointCompressor::new();
        let compressed = c.compress(&[]);
        let decompressed = c.decompress(&compressed).unwrap();
        assert!(decompressed.is_empty());
    }
    #[test]
    fn compressor_roundtrip_short() {
        let c = CheckpointCompressor::new();
        let input = b"hello world";
        let compressed = c.compress(input);
        let decompressed = c.decompress(&compressed).unwrap();
        assert_eq!(decompressed, input);
    }
    #[test]
    fn compressor_roundtrip_repeated_bytes() {
        let c = CheckpointCompressor::new();
        let input: Vec<u8> = vec![0xAB; 100];
        let compressed = c.compress(&input);
        let decompressed = c.decompress(&compressed).unwrap();
        assert_eq!(decompressed, input);
    }
    #[test]
    fn compressor_ratio_repeated_data_less_than_original() {
        let c = CheckpointCompressor::new();
        let pattern: Vec<u8> = vec![0xAA, 0xBB, 0xCC, 0xDD];
        let input: Vec<u8> = pattern.iter().cycle().take(256).cloned().collect();
        let compressed = c.compress(&input);
        let decompressed = c.decompress(&compressed).unwrap();
        assert_eq!(decompressed, input, "roundtrip must be lossless");
        let ratio = CheckpointCompressor::compression_ratio(input.len(), compressed.len());
        assert!(
            ratio < 2.0_f64,
            "compression ratio for repeated pattern should be reasonable, got {ratio}"
        );
    }
    #[test]
    fn compressor_ratio_zero_input() {
        let r = CheckpointCompressor::compression_ratio(0, 0);
        assert!((r - 1.0_f64).abs() < 1e-10);
    }
    #[test]
    fn compressor_roundtrip_random_like() {
        let c = CheckpointCompressor::new();
        let mut data = vec![0u8; 64];
        let mut x: u32 = 0xDEADBEEF;
        for b in &mut data {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            *b = (x >> 16) as u8;
        }
        let compressed = c.compress(&data);
        let decompressed = c.decompress(&compressed).unwrap();
        assert_eq!(decompressed, data);
    }
    #[test]
    fn compressor_roundtrip_long_literal() {
        let c = CheckpointCompressor::new();
        let input: Vec<u8> = (0u16..300)
            .map(|i| (i as u8).wrapping_add((i >> 8) as u8))
            .collect();
        let compressed = c.compress(&input);
        let decompressed = c.decompress(&compressed).unwrap();
        assert_eq!(decompressed, input);
    }
    #[test]
    fn catalog_empty_dir() {
        let dir = std::env::temp_dir().join("oxichk_catalog_empty");
        fs::create_dir_all(&dir).unwrap();
        let cat = CheckpointCatalog::scan(&dir);
        assert!(cat.is_empty());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_add_and_find() {
        let dir = std::env::temp_dir().join("oxichk_catalog_add");
        fs::create_dir_all(&dir).unwrap();
        let mut cat = CheckpointCatalog::scan(&dir);
        let ckpt = Checkpoint::new(1, 0, 5, 0.005_f64, vec![1, 2, 3]);
        cat.add(&ckpt).unwrap();
        assert!(!cat.is_empty());
        assert_eq!(cat.len(), 1);
        assert!(cat.path_for_step(5).is_some());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_load_step_roundtrip() {
        let dir = std::env::temp_dir().join("oxichk_catalog_load");
        fs::create_dir_all(&dir).unwrap();
        let mut cat = CheckpointCatalog::scan(&dir);
        let state = vec![10u8, 20, 30];
        let ckpt = Checkpoint::new(1, 1000, 7, 0.007_f64, state.clone());
        cat.add(&ckpt).unwrap();
        let loaded = cat.load_step(7).unwrap();
        assert_eq!(loaded.step, 7);
        assert_eq!(loaded.state, state);
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_steps_sorted() {
        let dir = std::env::temp_dir().join("oxichk_catalog_sorted");
        fs::create_dir_all(&dir).unwrap();
        let mut cat = CheckpointCatalog::scan(&dir);
        for step in [30u64, 10, 20] {
            cat.add(&Checkpoint::new(1, 0, step, 0.0_f64, vec![]))
                .unwrap();
        }
        let steps = cat.steps();
        assert_eq!(steps, vec![10, 20, 30]);
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_latest_and_earliest() {
        let dir = std::env::temp_dir().join("oxichk_catalog_latest");
        fs::create_dir_all(&dir).unwrap();
        let mut cat = CheckpointCatalog::scan(&dir);
        for step in [100u64, 50, 200] {
            cat.add(&Checkpoint::new(1, 0, step, 0.0_f64, vec![]))
                .unwrap();
        }
        let latest = cat.latest().unwrap();
        assert!(latest.to_string_lossy().contains("0000000200"));
        let earliest = cat.earliest().unwrap();
        assert!(earliest.to_string_lossy().contains("0000000050"));
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_remove_step() {
        let dir = std::env::temp_dir().join("oxichk_catalog_remove");
        fs::create_dir_all(&dir).unwrap();
        let mut cat = CheckpointCatalog::scan(&dir);
        cat.add(&Checkpoint::new(1, 0, 10, 0.0_f64, vec![]))
            .unwrap();
        cat.add(&Checkpoint::new(1, 0, 20, 0.0_f64, vec![]))
            .unwrap();
        cat.remove_step(10).unwrap();
        assert_eq!(cat.len(), 1);
        assert!(cat.path_for_step(10).is_none());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn catalog_not_found_step_returns_error() {
        let dir = std::env::temp_dir().join("oxichk_catalog_notfound");
        fs::create_dir_all(&dir).unwrap();
        let cat = CheckpointCatalog::scan(&dir);
        assert!(cat.load_step(999).is_err());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn file_writer_creates_bin_file() {
        let dir = std::env::temp_dir().join("oxichk_fw_test");
        fs::create_dir_all(&dir).unwrap();
        let writer = CheckpointFileWriter::new(&dir);
        let ckpt = Checkpoint::new(1, 0, 5, 0.005_f64, vec![9u8; 16]);
        let bin_path = writer.write(&ckpt).unwrap();
        assert!(bin_path.exists());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn file_writer_creates_json_sidecar() {
        let dir = std::env::temp_dir().join("oxichk_fw_json");
        fs::create_dir_all(&dir).unwrap();
        let writer = CheckpointFileWriter::new(&dir);
        let ckpt = Checkpoint::new(2, 500, 15, 0.015_f64, vec![1, 2, 3]);
        writer.write(&ckpt).unwrap();
        let json_path = dir.join("checkpoint_0000000015.json");
        assert!(json_path.exists());
        let json_data = fs::read_to_string(&json_path).unwrap();
        assert!(json_data.contains("\"step\":15"));
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn file_reader_validates_valid_checkpoint() {
        let dir = std::env::temp_dir().join("oxichk_fr_valid");
        fs::create_dir_all(&dir).unwrap();
        let writer = CheckpointFileWriter::new(&dir);
        let ckpt = Checkpoint::new(1, 0, 3, 0.003_f64, vec![7u8; 32]);
        let bin_path = writer.write(&ckpt).unwrap();
        let reader = CheckpointFileReader::new(&bin_path);
        let loaded = reader.read_and_validate().unwrap();
        assert_eq!(loaded.step, 3);
        assert!(loaded.verify());
        let _ = fs::remove_dir_all(&dir);
    }
    #[test]
    fn file_reader_detects_corrupted_checksum() {
        let dir = std::env::temp_dir().join("oxichk_fr_corrupt");
        fs::create_dir_all(&dir).unwrap();
        let writer = CheckpointFileWriter::new(&dir);
        let ckpt = Checkpoint::new(1, 0, 3, 0.003_f64, vec![7u8; 32]);
        let bin_path = writer.write(&ckpt).unwrap();
        let mut raw = fs::read(&bin_path).unwrap();
        if raw.len() > 32 {
            raw[32] ^= 0xFF;
        }
        fs::write(&bin_path, &raw).unwrap();
        let reader = CheckpointFileReader::new(&bin_path);
        assert!(reader.read_and_validate().is_err());
        let _ = fs::remove_dir_all(&dir);
    }

    // ── G2: Checkpoint compression ────────────────────────────────────────────

    #[test]
    fn compressed_checkpoint_roundtrip() {
        let path = tmp_path("compress_rt");
        let _ = fs::remove_file(&path);
        let meta = CheckpointMetadata::new(42, 0.042, 2, [0, 0, 0], "compress_test");
        let writer = CheckpointWriter::new(&path).with_compress(true);
        writer.write_header(&meta).unwrap();
        writer
            .write_positions(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            .unwrap();
        writer.write_scalars("energy", &[-7.5, -8.0]).unwrap();
        writer.finalize().unwrap();

        // File must start with zstd magic bytes after finalize with compress=true
        let raw = fs::read(&path).unwrap();
        assert_eq!(
            &raw[..4],
            &[0x28, 0xB5, 0x2F, 0xFD],
            "compressed file must start with zstd magic"
        );

        // Reader must transparently decompress and round-trip all data
        let reader = CheckpointReader::new(&path);
        let loaded_meta = reader.read_metadata().unwrap();
        assert_eq!(loaded_meta.step, 42);
        let positions = reader.read_positions().unwrap();
        assert_eq!(positions.len(), 2);
        assert!((positions[0][0] - 1.0).abs() < 1e-12);
        let energy = reader.read_scalars("energy").unwrap();
        assert_eq!(energy.len(), 2);
        assert!((energy[0] - (-7.5)).abs() < 1e-12);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn uncompressed_checkpoint_still_readable() {
        let path = tmp_path("no_compress_rt");
        let _ = fs::remove_file(&path);
        let meta = CheckpointMetadata::new(7, 0.007, 1, [0, 0, 0], "no_compress");
        let writer = CheckpointWriter::new(&path); // compress = false by default
        writer.write_header(&meta).unwrap();
        writer.write_positions(&[[0.0, 1.0, 0.0]]).unwrap();
        writer.finalize().unwrap();

        // Uncompressed file must NOT start with zstd magic
        let raw = fs::read(&path).unwrap();
        assert_ne!(
            &raw[..4],
            &[0x28, 0xB5, 0x2F, 0xFD],
            "uncompressed file must not start with zstd magic"
        );

        // Reader must still work correctly
        let reader = CheckpointReader::new(&path);
        let loaded_meta = reader.read_metadata().unwrap();
        assert_eq!(loaded_meta.step, 7);
        let positions = reader.read_positions().unwrap();
        assert_eq!(positions.len(), 1);
        let _ = fs::remove_file(&path);
    }

    // ── Integrity: footer CRC, magic, compressed restart, hostile counts ─────

    const ZSTD_FRAME_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];

    /// Serialise the file header (magic, version, metadata record).
    fn header_bytes(meta: &CheckpointMetadata) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&MAGIC.to_le_bytes());
        buf.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        let meta_bytes = meta.to_bytes();
        buf.extend_from_slice(&(meta_bytes.len() as u32).to_le_bytes());
        buf.extend_from_slice(&meta_bytes);
        buf
    }

    /// Append a valid footer: the footer tag plus the CRC-32 of `body`.
    fn seal(mut body: Vec<u8>) -> Vec<u8> {
        let csum = compute_checksum(&body);
        body.push(TAG_FOOTER);
        body.extend_from_slice(&csum.to_le_bytes());
        body
    }

    fn vec3_block(tag: u8, values: &[[f64; 3]]) -> Vec<u8> {
        let mut buf = vec![tag];
        buf.extend_from_slice(&(values.len() as u64).to_le_bytes());
        for v in values {
            for c in v {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        buf
    }

    fn scalar_block(name: &str, values: &[f64]) -> Vec<u8> {
        let mut buf = vec![TAG_SCALARS];
        buf.extend_from_slice(&(name.len() as u32).to_le_bytes());
        buf.extend_from_slice(name.as_bytes());
        buf.extend_from_slice(&(values.len() as u64).to_le_bytes());
        for v in values {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// A block header that claims `count` elements but carries no payload.
    fn hostile_header(tag: u8, name: &str, count: u64) -> Vec<u8> {
        let mut buf = vec![tag];
        if tag == TAG_SCALARS || tag == TAG_INTEGERS {
            buf.extend_from_slice(&(name.len() as u32).to_le_bytes());
            buf.extend_from_slice(name.as_bytes());
        }
        buf.extend_from_slice(&count.to_le_bytes());
        buf
    }

    /// Write `header + blocks` with a valid footer to `path`.
    fn write_sealed(path: &Path, blocks: &[Vec<u8>]) {
        let mut body = header_bytes(&sample_meta());
        for b in blocks {
            body.extend_from_slice(b);
        }
        fs::write(path, seal(body)).unwrap();
    }

    fn write_sample(path: &Path, compress: bool) {
        let writer = CheckpointWriter::new(path).with_compress(compress);
        writer.write_header(&sample_meta()).unwrap();
        writer
            .write_positions(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            .unwrap();
        writer
            .write_velocities(&[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
            .unwrap();
        writer.write_scalars("e", &[-7.5, 0.25, 1e-300]).unwrap();
        writer.write_integers("ids", &[3, 1, 4]).unwrap();
        writer.finalize().unwrap();
    }

    /// Offset of the first positions payload byte in a file written by
    /// [`write_sample`] without compression.
    fn positions_payload_offset() -> usize {
        header_bytes(&sample_meta()).len() + 1 + 8
    }

    fn assert_invalid<T>(result: io::Result<T>, label: &str) {
        match result {
            Ok(_) => panic!("{label}: a corrupt checkpoint was accepted"),
            Err(e) => assert_eq!(e.kind(), io::ErrorKind::InvalidData, "{label}: {e}"),
        }
    }

    /// Every reader that walks the block stream must reject the checkpoint at
    /// `path` with `InvalidData`.
    fn assert_block_readers_reject(path: &Path, scalar_name: &str) {
        let reader = CheckpointReader::new(path);
        assert_invalid(reader.read_positions(), "read_positions");
        assert_invalid(reader.read_velocities(), "read_velocities");
        assert_invalid(reader.read_scalars(scalar_name), "read_scalars");
        assert_invalid(reader.read_all(), "read_all");
        assert_invalid(RestartFile::load(path), "RestartFile::load");
    }

    /// Every public read path, including the header-only
    /// [`CheckpointReader::read_metadata`], must reject the checkpoint at
    /// `path` with `InvalidData`.
    fn assert_all_readers_reject(path: &Path, scalar_name: &str) {
        assert_invalid(CheckpointReader::new(path).read_metadata(), "read_metadata");
        assert_block_readers_reject(path, scalar_name);
    }

    #[test]
    fn checksum_matches_standard_check_value() {
        assert_eq!(compute_checksum(b"123456789"), 0xCBF4_3926);
    }

    #[test]
    fn checksum_matches_bitwise_reference() {
        fn bitwise(data: &[u8]) -> u32 {
            let mut crc: u32 = 0xFFFF_FFFF;
            for &byte in data {
                crc ^= u32::from(byte);
                for _ in 0..8 {
                    if crc & 1 != 0 {
                        crc = (crc >> 1) ^ 0xEDB8_8320;
                    } else {
                        crc >>= 1;
                    }
                }
            }
            crc ^ 0xFFFF_FFFF
        }
        let mut data = vec![0u8; 5000];
        let mut x: u32 = 0x1234_5678;
        for b in &mut data {
            x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            *b = (x >> 24) as u8;
        }
        for len in [0usize, 1, 2, 7, 255, 256, 257, 4096, 5000] {
            assert_eq!(
                compute_checksum(&data[..len]),
                bitwise(&data[..len]),
                "{len}"
            );
        }
    }

    #[test]
    fn restart_file_load_reads_compressed_checkpoint() {
        let path = tmp_path("restart_compressed");
        write_sample(&path, true);
        let raw = fs::read(&path).unwrap();
        assert_eq!(&raw[..4], &ZSTD_FRAME_MAGIC, "fixture must be compressed");

        let rf = RestartFile::load(&path).unwrap();
        assert_eq!(rf.meta, sample_meta());
        assert_eq!(rf.positions, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(rf.velocities, vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]);
        assert_eq!(rf.scalars.len(), 1);
        assert_eq!(rf.scalars["e"], vec![-7.5, 0.25, 1e-300]);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn read_all_matches_the_individual_readers() {
        for compress in [false, true] {
            let path = tmp_path(&format!("read_all_{compress}"));
            write_sample(&path, compress);
            let reader = CheckpointReader::new(&path);
            let all = reader.read_all().unwrap();
            assert_eq!(all.meta, reader.read_metadata().unwrap());
            assert_eq!(all.positions, reader.read_positions().unwrap());
            assert_eq!(all.velocities, reader.read_velocities().unwrap());
            assert_eq!(all.scalars["e"], reader.read_scalars("e").unwrap());
            let _ = fs::remove_file(&path);
        }
    }

    #[test]
    fn restart_file_save_load_roundtrip_with_integer_blocks_in_between() {
        let path = tmp_path("restart_with_integers");
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_positions(&[[1.0, 0.0, 0.0]]).unwrap();
        writer.write_integers("types", &[1, 2, 3]).unwrap();
        writer.write_scalars("rho", &[0.5, 1.5]).unwrap();
        writer.write_velocities(&[[0.0, 2.0, 0.0]]).unwrap();
        writer.finalize().unwrap();
        let rf = RestartFile::load(&path).unwrap();
        assert_eq!(rf.positions, vec![[1.0, 0.0, 0.0]]);
        assert_eq!(rf.velocities, vec![[0.0, 2.0, 0.0]]);
        assert_eq!(rf.scalars["rho"], vec![0.5, 1.5]);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn first_block_wins_when_a_name_or_kind_repeats() {
        let path = tmp_path("duplicates");
        write_sealed(
            &path,
            &[
                scalar_block("a", &[1.0]),
                scalar_block("a", &[2.0]),
                vec3_block(TAG_POSITIONS, &[[1.0, 1.0, 1.0]]),
                vec3_block(TAG_POSITIONS, &[[2.0, 2.0, 2.0]]),
            ],
        );
        let reader = CheckpointReader::new(&path);
        assert_eq!(reader.read_scalars("a").unwrap(), vec![1.0]);
        assert_eq!(reader.read_positions().unwrap(), vec![[1.0, 1.0, 1.0]]);
        let rf = reader.read_all().unwrap();
        assert_eq!(rf.scalars["a"], vec![1.0]);
        assert_eq!(rf.positions, vec![[1.0, 1.0, 1.0]]);
        assert!(rf.velocities.is_empty());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn finalising_twice_keeps_the_file_readable() {
        let path = tmp_path("finalize_twice");
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_positions(&[[1.0, 2.0, 3.0]]).unwrap();
        writer.finalize().unwrap();
        writer.finalize().unwrap();
        let reader = CheckpointReader::new(&path);
        assert_eq!(reader.read_positions().unwrap(), vec![[1.0, 2.0, 3.0]]);
        assert_eq!(reader.read_all().unwrap().positions.len(), 1);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn integer_arrays_round_trip_and_report_missing_names() {
        let path = tmp_path("integers_roundtrip");
        write_sample(&path, false);
        let reader = CheckpointReader::new(&path);
        assert_eq!(reader.read_integers("ids").unwrap(), vec![3, 1, 4]);
        let err = reader.read_integers("missing").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        let _ = fs::remove_file(&path);

        let negative = tmp_path("integers_negative");
        let writer = CheckpointWriter::new(&negative);
        writer.write_header(&sample_meta()).unwrap();
        writer
            .write_integers("signed", &[i32::MIN, -1, 0, i32::MAX])
            .unwrap();
        writer.finalize().unwrap();
        assert_eq!(
            CheckpointReader::new(&negative)
                .read_integers("signed")
                .unwrap(),
            vec![i32::MIN, -1, 0, i32::MAX]
        );
        let _ = fs::remove_file(&negative);
    }

    #[test]
    fn missing_scalar_in_a_valid_file_is_not_found() {
        let path = tmp_path("scalar_not_found");
        write_sample(&path, false);
        let err = CheckpointReader::new(&path)
            .read_scalars("nope")
            .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn flipped_payload_byte_is_detected_by_every_reader() {
        for compress in [false, true] {
            let path = tmp_path(&format!("bitrot_{compress}"));
            write_sample(&path, compress);
            let mut bytes = fs::read(&path).unwrap();
            if compress {
                bytes = oxiarc_zstd::decompress(&bytes).unwrap();
            }
            bytes[positions_payload_offset()] ^= 0x01;
            if compress {
                bytes = oxiarc_zstd::compress_with_level(&bytes, 3).unwrap();
            }
            fs::write(&path, &bytes).unwrap();
            assert_all_readers_reject(&path, "e");
            let _ = fs::remove_file(&path);
        }
    }

    #[test]
    fn checksum_mismatch_message_names_the_checksum() {
        let path = tmp_path("bitrot_message");
        write_sample(&path, false);
        let mut bytes = fs::read(&path).unwrap();
        let last = bytes.len() - 1;
        bytes[last] ^= 0x80;
        fs::write(&path, &bytes).unwrap();
        let err = CheckpointReader::new(&path).read_positions().unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("checksum mismatch"), "{err}");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn unfinalized_checkpoint_is_rejected() {
        let path = tmp_path("unfinalized");
        let writer = CheckpointWriter::new(&path);
        writer.write_header(&sample_meta()).unwrap();
        writer.write_positions(&[[1.0, 2.0, 3.0]]).unwrap();
        writer.write_scalars("e", &[1.0]).unwrap();
        assert_all_readers_reject(&path, "e");
        let err = CheckpointReader::new(&path).read_all().unwrap_err();
        assert!(err.to_string().contains("not finalized"), "{err}");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn truncated_checkpoint_is_rejected() {
        for cut in [1usize, FOOTER_LEN, FOOTER_LEN + 3] {
            let path = tmp_path(&format!("truncated_{cut}"));
            write_sample(&path, false);
            let bytes = fs::read(&path).unwrap();
            fs::write(&path, &bytes[..bytes.len() - cut]).unwrap();
            assert_all_readers_reject(&path, "e");
            let _ = fs::remove_file(&path);
        }
    }

    #[test]
    fn bad_magic_with_a_valid_footer_is_rejected_by_every_reader() {
        let path = tmp_path("bad_magic_sealed");
        let mut body = header_bytes(&sample_meta());
        body[..4].copy_from_slice(&0xDEAD_BEEFu32.to_le_bytes());
        body.extend_from_slice(&vec3_block(TAG_POSITIONS, &[[1.0, 2.0, 3.0]]));
        fs::write(&path, seal(body)).unwrap();
        assert_all_readers_reject(&path, "e");
        let err = CheckpointReader::new(&path).read_positions().unwrap_err();
        assert!(err.to_string().contains("bad magic"), "{err}");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn files_shorter_than_the_header_are_invalid_data() {
        for len in [0usize, 1, 3, 4, 8] {
            let path = tmp_path(&format!("tiny_{len}"));
            fs::write(&path, vec![0u8; len]).unwrap();
            assert_all_readers_reject(&path, "e");
            let _ = fs::remove_file(&path);
        }
    }

    #[test]
    fn a_foreign_file_is_reported_as_bad_magic_not_as_unfinalized() {
        let path = tmp_path("foreign_file");
        fs::write(&path, b"This is a text file, not an OxiPhysics checkpoint.").unwrap();
        for result in [
            CheckpointReader::new(&path).read_positions().map(|_| ()),
            CheckpointReader::new(&path).read_metadata().map(|_| ()),
            RestartFile::load(&path).map(|_| ()),
        ] {
            let err = result.unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::InvalidData);
            assert!(err.to_string().contains("bad magic"), "{err}");
        }
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn unknown_block_tag_is_invalid_data() {
        let path = tmp_path("unknown_tag");
        write_sealed(&path, &[vec![0x7E, 0, 0, 0, 0]]);
        assert_block_readers_reject(&path, "e");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn oversized_vector_count_is_rejected_without_panicking() {
        for tag in [TAG_POSITIONS, TAG_VELOCITIES] {
            for count in [u64::MAX, u64::MAX / 24 + 1, 1u64 << 40, 1_000_000] {
                let path = tmp_path(&format!("huge_vec_{tag}_{count}"));
                write_sealed(&path, &[hostile_header(tag, "", count)]);
                assert_block_readers_reject(&path, "e");
                let _ = fs::remove_file(&path);
            }
        }
    }

    #[test]
    fn wrapping_scalar_count_is_rejected_promptly() {
        // 2^61 - 2 elements of 8 bytes wrap to -16 modulo 2^64: with unchecked
        // arithmetic the cursor lands back on the block tag and the walker
        // loops forever.
        let path = tmp_path("wrapping_scalars");
        write_sealed(
            &path,
            &[hostile_header(TAG_SCALARS, "abc", (1u64 << 61) - 2)],
        );
        assert_block_readers_reject(&path, "abc");
        assert_block_readers_reject(&path, "other");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn oversized_integer_and_scalar_counts_are_rejected() {
        for tag in [TAG_INTEGERS, TAG_SCALARS] {
            for count in [u64::MAX, u64::MAX / 2 + 1, 1u64 << 62, 5_000_000_000] {
                let path = tmp_path(&format!("huge_arr_{tag}_{count}"));
                write_sealed(&path, &[hostile_header(tag, "n", count)]);
                assert_block_readers_reject(&path, "n");
                let _ = fs::remove_file(&path);
            }
        }
    }

    #[test]
    fn hostile_count_after_a_valid_block_is_rejected() {
        let path = tmp_path("hostile_after_valid");
        write_sealed(
            &path,
            &[
                vec3_block(TAG_VELOCITIES, &[[1.0, 2.0, 3.0]]),
                hostile_header(TAG_POSITIONS, "", u64::MAX),
            ],
        );
        let reader = CheckpointReader::new(&path);
        // The velocities block precedes the damage and stays readable.
        assert_eq!(reader.read_velocities().unwrap(), vec![[1.0, 2.0, 3.0]]);
        assert_invalid(reader.read_positions(), "read_positions");
        assert_invalid(reader.read_all(), "read_all");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn hostile_count_inside_a_compressed_checkpoint_is_rejected() {
        let path = tmp_path("hostile_compressed");
        let mut body = header_bytes(&sample_meta());
        body.extend_from_slice(&hostile_header(TAG_SCALARS, "abc", (1u64 << 61) - 2));
        let compressed = oxiarc_zstd::compress_with_level(&seal(body), 3).unwrap();
        fs::write(&path, &compressed).unwrap();
        assert_block_readers_reject(&path, "abc");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn oversized_metadata_length_is_rejected() {
        let path = tmp_path("huge_meta_len");
        let mut body = Vec::new();
        body.extend_from_slice(&MAGIC.to_le_bytes());
        body.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        body.extend_from_slice(&u32::MAX.to_le_bytes());
        fs::write(&path, seal(body)).unwrap();
        let reader = CheckpointReader::new(&path);
        assert!(reader.read_metadata().is_err());
        assert!(reader.read_positions().is_err());
        assert!(reader.read_all().is_err());
        let _ = fs::remove_file(&path);
    }

    /// Every public read path, `read_integers` included, must reject the
    /// checkpoint at `path` with `InvalidData`; a failure names `label`.
    fn assert_every_reader_rejects(path: &Path, label: &str) {
        let reader = CheckpointReader::new(path);
        assert_invalid(reader.read_metadata(), label);
        assert_invalid(reader.read_positions(), label);
        assert_invalid(reader.read_velocities(), label);
        assert_invalid(reader.read_scalars("e"), label);
        assert_invalid(reader.read_integers("ids"), label);
        assert_invalid(reader.read_all(), label);
        assert_invalid(RestartFile::load(path), label);
    }

    #[test]
    fn every_single_bit_flip_of_a_finalized_file_is_rejected() {
        // The magic check, the CRC-32 (which detects every single-bit error)
        // and the footer tag between them cover each byte of the file.
        let path = tmp_path("bit_flip_sweep");
        write_sample(&path, false);
        let pristine = fs::read(&path).unwrap();
        assert!(CheckpointReader::new(&path).read_all().is_ok());
        for byte in 0..pristine.len() {
            for bit in 0..8u32 {
                let mut damaged = pristine.clone();
                damaged[byte] ^= 1u8 << bit;
                fs::write(&path, &damaged).unwrap();
                assert_every_reader_rejects(&path, &format!("byte {byte}, bit {bit}"));
            }
        }
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn every_proper_prefix_of_a_finalized_file_is_rejected() {
        let path = tmp_path("prefix_sweep");
        write_sample(&path, false);
        let pristine = fs::read(&path).unwrap();
        for keep in 0..pristine.len() {
            fs::write(&path, &pristine[..keep]).unwrap();
            assert_every_reader_rejects(&path, &format!("prefix of {keep} bytes"));
        }
        fs::write(&path, &pristine).unwrap();
        assert!(CheckpointReader::new(&path).read_all().is_ok());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn block_counts_around_every_power_of_two_are_rejected() {
        // A non-empty block header with no payload behind it must fail for any
        // claimed count, whichever step of `count * size + cursor` would wrap.
        let path = tmp_path("power_of_two_counts");
        for tag in [TAG_POSITIONS, TAG_VELOCITIES, TAG_SCALARS, TAG_INTEGERS] {
            for name in ["", "abc"] {
                for shift in 0..64u32 {
                    let power = 1u64 << shift;
                    for count in [power - 1, power, power + 1] {
                        if count == 0 {
                            continue;
                        }
                        write_sealed(&path, &[hostile_header(tag, name, count)]);
                        assert_block_readers_reject(&path, name);
                    }
                }
            }
        }
        let _ = fs::remove_file(&path);
    }

    fn xorshift64(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn restart_file_round_trip_preserves_every_f64_bit_pattern() {
        // Random 64-bit patterns include NaNs with payloads, infinities,
        // subnormals and negative zero; none may be altered on the way.
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = || f64::from_bits(xorshift64(&mut state));
        let positions: Vec<[f64; 3]> = (0..1500).map(|_| [next(), next(), next()]).collect();
        let velocities: Vec<[f64; 3]> = (0..1500).map(|_| [next(), next(), next()]).collect();
        let noise: Vec<f64> = (0..700).map(|_| next()).collect();
        let special = vec![
            0.0,
            -0.0,
            f64::MIN_POSITIVE,
            5e-324,
            f64::MAX,
            f64::MIN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            f64::from_bits(0x7FF8_0000_0000_1234),
        ];
        let mut scalars = HashMap::new();
        scalars.insert("special".to_string(), special.clone());
        scalars.insert("noise".to_string(), noise.clone());
        let meta = CheckpointMetadata::new(11, 0.25, positions.len(), [0, 1, 3], "bit patterns");
        let path = tmp_path("bit_patterns");
        RestartFile::new(meta.clone(), positions.clone(), velocities.clone(), scalars)
            .save(&path)
            .unwrap();
        let loaded = RestartFile::load(&path).unwrap();
        assert_eq!(loaded.meta, meta);
        let bits3 =
            |v: &[[f64; 3]]| -> Vec<[u64; 3]> { v.iter().map(|p| p.map(f64::to_bits)).collect() };
        let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
        assert_eq!(bits3(&loaded.positions), bits3(&positions));
        assert_eq!(bits3(&loaded.velocities), bits3(&velocities));
        assert_eq!(loaded.scalars.len(), 2);
        assert_eq!(bits(&loaded.scalars["special"]), bits(&special));
        assert_eq!(bits(&loaded.scalars["noise"]), bits(&noise));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn take_bytes_rejects_cursor_overflow_without_moving() {
        let data = [1u8, 2, 3, 4];
        let mut cursor = usize::MAX - 1;
        assert!(take_bytes(&data, &mut cursor, 4, "t").is_err());
        assert_eq!(cursor, usize::MAX - 1);
        let mut cursor = 1usize;
        assert!(take_bytes(&data, &mut cursor, usize::MAX, "t").is_err());
        assert_eq!(cursor, 1);
        let mut cursor = 1usize;
        assert_eq!(take_bytes(&data, &mut cursor, 3, "t").unwrap(), &[2, 3, 4]);
        assert_eq!(cursor, 4);
    }

    #[test]
    fn scalar_readers_reject_out_of_range_cursors() {
        let data = [0u8; 16];
        // (start offset, whether a 4-byte / 8-byte read fits)
        let cases = [
            (0usize, true, true),
            (8, true, true),
            (12, true, false),
            (13, false, false),
            (16, false, false),
            (17, false, false),
            (usize::MAX - 7, false, false),
            (usize::MAX - 2, false, false),
            (usize::MAX, false, false),
        ];
        for (start, fits4, fits8) in cases {
            let mut c = start;
            assert_eq!(read_u32(&data, &mut c).is_ok(), fits4, "u32 at {start}");
            let mut c = start;
            assert_eq!(read_u64(&data, &mut c).is_ok(), fits8, "u64 at {start}");
            let mut c = start;
            assert_eq!(read_f64(&data, &mut c).is_ok(), fits8, "f64 at {start}");
        }
        let mut c = usize::MAX;
        assert!(read_u32(&data, &mut c).is_err());
        assert_eq!(c, usize::MAX, "a failed read must not move the cursor");
    }

    #[test]
    fn read_name_rejects_a_huge_length_prefix() {
        let mut data = u32::MAX.to_le_bytes().to_vec();
        data.extend_from_slice(b"abc");
        let mut cursor = 0usize;
        let err = read_name(&data, &mut cursor).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        let mut ok = 3u32.to_le_bytes().to_vec();
        ok.extend_from_slice(b"abc");
        let mut cursor = 0usize;
        assert_eq!(read_name(&ok, &mut cursor).unwrap(), "abc");
        assert_eq!(cursor, 7);
    }

    #[test]
    fn checked_block_end_guards_every_step() {
        assert_eq!(checked_block_end(100, 10, 3, 24).unwrap(), (3, 82));
        assert_eq!(checked_block_end(34, 10, 1, 24).unwrap(), (1, 34));
        assert_eq!(checked_block_end(34, 34, 0, 24).unwrap(), (0, 34));
        for (len, cursor, count, elem) in [
            (33usize, 10usize, 1u64, 24usize),
            (100, 0, u64::MAX, 8),
            (100, 0, (1u64 << 61) - 2, 8),
            (100, usize::MAX - 4, 1, 8),
            (usize::MAX, usize::MAX - 4, 1, 8),
            (100, 0, 1u64 << 63, 2),
        ] {
            let err = checked_block_end(len, cursor, count, elem).unwrap_err();
            assert_eq!(
                err.kind(),
                io::ErrorKind::InvalidData,
                "{len} {cursor} {count}"
            );
        }
    }

    #[test]
    fn verify_footer_returns_the_protected_body() {
        let sealed = seal(vec![1, 2, 3]);
        assert_eq!(verify_footer(&sealed).unwrap(), &[1, 2, 3]);
        assert!(verify_footer(&seal(Vec::new())).unwrap().is_empty());
        for short in [&sealed[..0], &sealed[..4]] {
            assert_eq!(
                verify_footer(short).unwrap_err().kind(),
                io::ErrorKind::InvalidData
            );
        }
        let mut wrong_tag = sealed.clone();
        let tag_at = wrong_tag.len() - FOOTER_LEN;
        wrong_tag[tag_at] = 0x00;
        assert_eq!(
            verify_footer(&wrong_tag).unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
        let mut wrong_body = sealed;
        wrong_body[0] ^= 0xFF;
        let err = verify_footer(&wrong_body).unwrap_err();
        assert!(err.to_string().contains("checksum mismatch"), "{err}");
    }

    #[test]
    fn checkpoint_from_bytes_rejects_a_huge_state_length() {
        let ckpt = Checkpoint::new(1, 2, 3, 0.5, vec![1, 2, 3]);
        let mut bytes = ckpt.to_bytes();
        // version(4) + timestamp(8) + step(8) + sim_time(8) + checksum(4)
        let len_at = 32;
        bytes[len_at..len_at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(Checkpoint::from_bytes(&bytes).is_err());
        bytes[len_at..len_at + 8].copy_from_slice(&4u64.to_le_bytes());
        assert!(Checkpoint::from_bytes(&bytes).is_err());
        bytes[len_at..len_at + 8].copy_from_slice(&3u64.to_le_bytes());
        assert_eq!(Checkpoint::from_bytes(&bytes).unwrap(), ckpt);
    }

    #[test]
    fn header_and_metadata_reject_huge_string_lengths() {
        let header = CheckpointHeader::new([0, 1, 0], 5, 6, "crate", 7);
        let mut bytes = header.to_bytes();
        // three u32 versions + timestamp + step
        let name_len_at = 12 + 8 + 8;
        bytes[name_len_at..name_len_at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(CheckpointHeader::from_bytes(&bytes).is_err());
        assert_eq!(
            CheckpointHeader::from_bytes(&header.to_bytes()).unwrap(),
            header
        );

        let meta = sample_meta();
        let mut bytes = meta.to_bytes();
        // step + time + n_particles + three u32 versions
        let ts_len_at = 8 + 8 + 8 + 12;
        bytes[ts_len_at..ts_len_at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(CheckpointMetadata::from_bytes(&bytes).is_err());
        assert_eq!(
            CheckpointMetadata::from_bytes(&meta.to_bytes()).unwrap(),
            meta
        );
    }

    #[test]
    fn simulation_state_from_bytes_rejects_a_huge_particle_count() {
        let mut bytes = u64::MAX.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 72]);
        assert!(SimulationState::from_bytes(&bytes).is_err());
        let mut bytes = 1_000_000u64.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 72]);
        assert!(SimulationState::from_bytes(&bytes).is_err());
    }

    #[test]
    fn simulation_state_bytes_roundtrip() {
        let mut state = SimulationState::new();
        state.positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        state.velocities = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        state.forces = vec![[-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0]];
        let restored = SimulationState::from_bytes(&state.to_bytes()).unwrap();
        assert_eq!(restored.positions, state.positions);
        assert_eq!(restored.velocities, state.velocities);
        assert_eq!(restored.forces, state.forces);
        let empty = SimulationState::from_bytes(&SimulationState::new().to_bytes()).unwrap();
        assert!(empty.is_empty());
    }
}
