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
pub(super) fn read_u32(data: &[u8], cursor: &mut usize) -> io::Result<u32> {
    if *cursor + 4 > data.len() {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "read_u32"));
    }
    let v = u32::from_le_bytes(
        data[*cursor..*cursor + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *cursor += 4;
    Ok(v)
}
pub(super) fn read_u64(data: &[u8], cursor: &mut usize) -> io::Result<u64> {
    if *cursor + 8 > data.len() {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "read_u64"));
    }
    let v = u64::from_le_bytes(
        data[*cursor..*cursor + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *cursor += 8;
    Ok(v)
}
pub(super) fn read_f64(data: &[u8], cursor: &mut usize) -> io::Result<f64> {
    if *cursor + 8 > data.len() {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "read_f64"));
    }
    let bits = u64::from_le_bytes(
        data[*cursor..*cursor + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *cursor += 8;
    Ok(f64::from_bits(bits))
}
/// Compute a simple CRC32-like checksum using polynomial XOR over the data.
///
/// Uses the standard CRC-32/ISO-HDLC polynomial (0xEDB88320) with an XOR-out
/// of 0xFFFFFFFF.
pub fn compute_checksum(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in data {
        crc ^= byte as u32;
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
pub(super) fn write_name(f: &mut impl Write, name: &str) -> io::Result<()> {
    let bytes = name.as_bytes();
    f.write_all(&(bytes.len() as u32).to_le_bytes())?;
    f.write_all(bytes)
}
pub(super) fn read_name(data: &[u8], cursor: &mut usize) -> io::Result<String> {
    let len = read_u32(data, cursor)? as usize;
    if *cursor + len > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "name truncated",
        ));
    }
    let s = String::from_utf8(data[*cursor..*cursor + len].to_vec())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("UTF-8: {e}")))?;
    *cursor += len;
    Ok(s)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::checkpoint_io::types::*;
    use std::collections::HashMap;
    use std::fs;
    use std::path::PathBuf;
    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("oxiphysics_ckpt_test_{name}.bin"));
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
        let mgr = CheckpointManager::new("/tmp", 5, 0);
        assert!(!mgr.should_checkpoint(0));
        assert!(!mgr.should_checkpoint(100));
    }
    #[test]
    fn checkpoint_path_format() {
        let mgr = CheckpointManager::new("/tmp/sim", 3, 50);
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
        let w = CheckpointWriter::new("/tmp/x.bin").with_compress(true);
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
}
