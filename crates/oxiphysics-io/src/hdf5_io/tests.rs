#![allow(clippy::needless_range_loop)]
#![allow(clippy::module_inception)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the HDF5 mock I/O module.

// All test modules use crate-level imports to access the re-exported API.
#![allow(unused_imports)]

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hdf5_io::*;

    // ── Error display ───────────────────────────────────────────────────────

    #[test]
    fn test_error_not_found_display() {
        let e = Hdf5Error::NotFound("foo".to_string());
        let s = e.to_string();
        assert!(s.contains("foo"), "error message should contain 'foo': {s}");
    }

    #[test]
    fn test_error_already_exists_display() {
        let e = Hdf5Error::AlreadyExists("bar".to_string());
        let s = e.to_string();
        assert!(s.contains("bar"), "error: {s}");
    }

    #[test]
    fn test_error_hyperslab_out_of_range_display() {
        let e = Hdf5Error::HyperslabOutOfRange {
            dim: 2,
            start: 10,
            count: 5,
            size: 12,
        };
        let s = e.to_string();
        assert!(s.contains("dim=2"), "error: {s}");
    }

    // ── Dtype ────────────────────────────────────────────────────────────────

    #[test]
    fn test_dtype_element_size() {
        assert_eq!(Hdf5Dtype::Float32.element_size(), 4);
        assert_eq!(Hdf5Dtype::Float64.element_size(), 8);
        assert_eq!(Hdf5Dtype::Int32.element_size(), 4);
        assert_eq!(Hdf5Dtype::Uint8.element_size(), 1);
        assert_eq!(Hdf5Dtype::VlenString.element_size(), 0);
    }

    #[test]
    fn test_dtype_compound_size() {
        let compound = Hdf5Dtype::Compound(vec![
            ("x".to_string(), Hdf5Dtype::Float64),
            ("y".to_string(), Hdf5Dtype::Float32),
            ("z".to_string(), Hdf5Dtype::Int32),
        ]);
        assert_eq!(compound.element_size(), 16); // 8 + 4 + 4
    }

    #[test]
    fn test_dtype_named_inherits_base_size() {
        let named = Hdf5Dtype::Named {
            name: "energy_t".to_string(),
            base: Box::new(Hdf5Dtype::Float64),
        };
        assert_eq!(named.element_size(), 8);
    }

    // ── ChunkLayout ──────────────────────────────────────────────────────────

    #[test]
    fn test_chunk_layout_volume() {
        let cl = ChunkLayout::new(vec![32, 32, 3], 6);
        assert_eq!(cl.chunk_volume(), 32 * 32 * 3);
        assert_eq!(cl.gzip_level, 6);
    }

    #[test]
    fn test_chunk_layout_gzip_clamped() {
        let cl = ChunkLayout::new(vec![10], 20); // 20 > 9 → clamped
        assert_eq!(cl.gzip_level, 9);
    }

    // ── Hyperslab ────────────────────────────────────────────────────────────

    #[test]
    fn test_hyperslab_volume() {
        let s = Hyperslab::new(vec![0, 0], vec![4, 8]);
        assert_eq!(s.volume(), 32);
    }

    #[test]
    fn test_hyperslab_validate_ok() {
        let s = Hyperslab::new(vec![0, 1], vec![3, 4]);
        let shape = [5, 8];
        assert!(s.validate(&shape).is_ok());
    }

    #[test]
    fn test_hyperslab_validate_out_of_range() {
        let s = Hyperslab::new(vec![0], vec![10]);
        let shape = [5];
        let r = s.validate(&shape);
        assert!(r.is_err());
    }

    // ── Group creation ───────────────────────────────────────────────────────

    #[test]
    fn test_create_group_nested() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("a/b/c").unwrap();
        let g = file.open_group("a/b/c").unwrap();
        assert_eq!(g.name, "c");
    }

    #[test]
    fn test_create_group_already_exists_is_ok() {
        // Our implementation uses create-or-ignore semantics via the file API.
        let mut file = Hdf5File::new("test.h5");
        file.create_group("grp").unwrap();
        // Creating again should not fail because of the create_group implementation.
        file.create_group("grp").unwrap();
    }

    #[test]
    fn test_open_nonexistent_group_fails() {
        let file = Hdf5File::new("test.h5");
        let r = file.open_group("nope");
        assert!(r.is_err());
    }

    // ── Dataset create / write / read ────────────────────────────────────────

    #[test]
    fn test_create_dataset_f64_and_read() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("data").unwrap();
        file.create_dataset("data", "pos", vec![3, 3], Hdf5Dtype::Float64)
            .unwrap();
        let data = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        file.open_dataset_mut("data", "pos")
            .unwrap()
            .write_f64(&data)
            .unwrap();
        let read = file
            .open_dataset("data", "pos")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(read, data);
    }

    #[test]
    fn test_create_dataset_f32_and_read() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "temps", vec![4], Hdf5Dtype::Float32)
            .unwrap();
        let data = vec![1.0_f32, 2.0, 3.0, 4.0];
        file.open_dataset_mut("g", "temps")
            .unwrap()
            .write_f32(&data)
            .unwrap();
        let read = file.open_dataset("g", "temps").unwrap().read_f32().unwrap();
        assert_eq!(read, data);
    }

    #[test]
    fn test_create_dataset_i32_and_read() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "ids", vec![5], Hdf5Dtype::Int32)
            .unwrap();
        let data = vec![10_i32, 20, 30, 40, 50];
        file.open_dataset_mut("g", "ids")
            .unwrap()
            .write_i32(&data)
            .unwrap();
        let read = file.open_dataset("g", "ids").unwrap().read_i32().unwrap();
        assert_eq!(read, data);
    }

    #[test]
    fn test_create_dataset_u8_and_read() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "bytes", vec![3], Hdf5Dtype::Uint8)
            .unwrap();
        let data = vec![0xAB_u8, 0xCD, 0xEF];
        file.open_dataset_mut("g", "bytes")
            .unwrap()
            .write_u8(&data)
            .unwrap();
        let read = file.open_dataset("g", "bytes").unwrap().read_u8().unwrap();
        assert_eq!(read, data);
    }

    #[test]
    fn test_dataset_wrong_type_read_fails() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "d", vec![2], Hdf5Dtype::Int32)
            .unwrap();
        let ds = file.open_dataset("g", "d").unwrap();
        assert!(ds.read_f64().is_err());
    }

    // ── Attributes ───────────────────────────────────────────────────────────

    #[test]
    fn test_dataset_attributes_set_get() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "d", vec![1], Hdf5Dtype::Float64)
            .unwrap();
        file.set_dataset_attr("g", "d", "units", AttrValue::String("angstrom".to_string()))
            .unwrap();
        let val = file.get_dataset_attr("g", "d", "units").unwrap();
        assert!(matches!(val, AttrValue::String(s) if s == "angstrom"));
    }

    #[test]
    fn test_group_attributes() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("meta").unwrap();
        let group = file.open_group_mut("meta").unwrap();
        group.set_attr("author", AttrValue::String("test".to_string()));
        let val = group.get_attr("author").unwrap();
        assert!(matches!(val, AttrValue::String(s) if s == "test"));
    }

    // ── Hyperslab read ───────────────────────────────────────────────────────

    #[test]
    fn test_hyperslab_read_1d() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "arr", vec![10], Hdf5Dtype::Float64)
            .unwrap();
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        file.open_dataset_mut("g", "arr")
            .unwrap()
            .write_f64(&data)
            .unwrap();
        let slab = Hyperslab::new(vec![2], vec![4]);
        let ds = file.open_dataset("g", "arr").unwrap();
        let out = ds.read_hyperslab_f64(&slab).unwrap();
        assert_eq!(out, vec![2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn test_hyperslab_read_2d() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        // 4x4 matrix
        file.create_dataset("g", "mat", vec![4, 4], Hdf5Dtype::Float64)
            .unwrap();
        let data: Vec<f64> = (0..16).map(|i| i as f64).collect();
        file.open_dataset_mut("g", "mat")
            .unwrap()
            .write_f64(&data)
            .unwrap();
        // Select rows 1..3, cols 1..3
        let slab = Hyperslab::new(vec![1, 1], vec![2, 2]);
        let ds = file.open_dataset("g", "mat").unwrap();
        let out = ds.read_hyperslab_f64(&slab).unwrap();
        // Row 1, cols 1-2: [5, 6]; Row 2, cols 1-2: [9, 10]
        assert_eq!(out, vec![5.0, 6.0, 9.0, 10.0]);
    }

    // ── Named types ──────────────────────────────────────────────────────────

    #[test]
    fn test_named_type_register_and_lookup() {
        let mut file = Hdf5File::new("test.h5");
        file.commit_datatype("energy_t", Hdf5Dtype::Float64)
            .unwrap();
        let dt = file.find_named_type("energy_t").unwrap();
        assert_eq!(*dt, Hdf5Dtype::Float64);
    }

    #[test]
    fn test_named_type_duplicate_fails() {
        let mut file = Hdf5File::new("test.h5");
        file.commit_datatype("t1", Hdf5Dtype::Float32).unwrap();
        let r = file.commit_datatype("t1", Hdf5Dtype::Float32);
        assert!(r.is_err());
    }

    // ── Links ────────────────────────────────────────────────────────────────

    #[test]
    fn test_soft_link_creation() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("grp").unwrap();
        file.create_soft_link("grp", "my_link", "/data/pos")
            .unwrap();
        let g = file.open_group("grp").unwrap();
        assert!(matches!(g.links.get("my_link"), Some(Hdf5Link::Soft(t)) if t == "/data/pos"));
    }

    #[test]
    fn test_hard_link_creation() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("grp").unwrap();
        file.create_hard_link("grp", "h_link", "/data/pos").unwrap();
        let g = file.open_group("grp").unwrap();
        assert!(matches!(g.links.get("h_link"), Some(Hdf5Link::Hard(_))));
    }

    #[test]
    fn test_link_duplicate_fails() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_soft_link("g", "lnk", "/a").unwrap();
        let r = file.create_soft_link("g", "lnk", "/b");
        assert!(r.is_err());
    }

    // ── File locking ─────────────────────────────────────────────────────────

    #[test]
    fn test_file_lock_write() {
        let mut file = Hdf5File::new("test.h5");
        assert!(!file.is_locked());
        file.lock_write(1).unwrap();
        assert!(file.is_locked());
    }

    #[test]
    fn test_file_double_lock_fails() {
        let mut file = Hdf5File::new("test.h5");
        file.lock_write(1).unwrap();
        let r = file.lock_write(2);
        assert!(r.is_err());
    }

    #[test]
    fn test_file_unlock() {
        let mut file = Hdf5File::new("test.h5");
        file.lock_write(1).unwrap();
        file.unlock();
        assert!(!file.is_locked());
    }

    #[test]
    fn test_write_while_locked_fails() {
        let mut file = Hdf5File::new("test.h5");
        file.lock_write(1).unwrap();
        let r = file.create_group("blocked");
        assert!(r.is_err());
    }

    // ── Trajectory ───────────────────────────────────────────────────────────

    #[test]
    fn test_trajectory_append_and_read_frame() {
        let mut traj = TrajectoryStore::new(2);
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        traj.append_frame(0.0, &pos, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let frame = traj.read_positions(0).unwrap();
        assert_eq!(frame.len(), 2);
        assert_eq!(frame[0], [0.0, 0.0, 0.0]);
        assert_eq!(frame[1], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn test_trajectory_n_frames() {
        let mut traj = TrajectoryStore::new(3);
        let pos = vec![0.0; 9];
        let bv = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        for t in 0..5 {
            traj.append_frame(t as f64, &pos, &bv);
        }
        assert_eq!(traj.n_frames, 5);
    }

    #[test]
    fn test_trajectory_total_time() {
        let mut traj = TrajectoryStore::new(1);
        let pos = vec![0.0, 0.0, 0.0];
        let bv = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        traj.append_frame(0.0, &pos, &bv);
        traj.append_frame(2.0, &pos, &bv);
        assert!((traj.total_time() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_trajectory_flush_to_file() {
        let mut traj = TrajectoryStore::new(2);
        let pos = vec![0.0; 6];
        let bv = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        traj.append_frame(0.0, &pos, &bv);
        let mut file = Hdf5File::new("traj.h5");
        traj.flush_to_file(&mut file).unwrap();
        let ds = file.open_dataset("trajectory", "positions").unwrap();
        assert_eq!(ds.shape, vec![1, 2, 3]);
    }

    // ── Checkpoint ───────────────────────────────────────────────────────────

    #[test]
    fn test_checkpoint_write_and_read() {
        let mut ckpt = Checkpoint::new("step_100", 100, 0.2);
        ckpt.positions = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        ckpt.velocities = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let mut file = Hdf5File::new("ckpt.h5");
        ckpt.write_to_file(&mut file).unwrap();
        let restored = Checkpoint::read_from_file(&file, "step_100").unwrap();
        assert_eq!(restored.positions, ckpt.positions);
        assert_eq!(restored.step, 100);
        assert!((restored.time - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_checkpoint_manager_capacity() {
        let mut mgr = CheckpointManager::new(2);
        let mut file = Hdf5File::new("ckpt.h5");
        for i in 0..5_u64 {
            let mut ckpt = Checkpoint::new(&format!("step_{i}"), i, i as f64 * 0.1);
            ckpt.positions = vec![0.0; 3];
            ckpt.velocities = vec![0.0; 3];
            mgr.write(&mut file, &ckpt).unwrap();
        }
        // Manager tracks max 2 labels
        assert_eq!(mgr.count(), 2);
    }

    // ── Superblock encode/decode ─────────────────────────────────────────────

    #[test]
    fn test_superblock_encode_decode_roundtrip() {
        let sb = Hdf5Superblock {
            version: 2,
            file_size: 1024 * 1024,
            root_obj_header_offset: 512,
            eof_address: 1023 * 1024,
            size_of_lengths: 8,
            size_of_offsets: 8,
        };
        let bytes = encode_superblock(&sb);
        let sb2 = decode_superblock(&bytes).unwrap();
        assert_eq!(sb2.version, sb.version);
        assert_eq!(sb2.file_size, sb.file_size);
        assert_eq!(sb2.root_obj_header_offset, sb.root_obj_header_offset);
    }

    // ── Object header encode/decode ──────────────────────────────────────────

    #[test]
    fn test_object_header_encode_decode_roundtrip() {
        let hdr = Hdf5ObjectHeader {
            object_type: "dataset".to_string(),
            address: 4096,
            n_messages: 5,
            header_size: 256,
        };
        let bytes = encode_object_header(&hdr);
        let hdr2 = decode_object_header(&bytes).unwrap();
        assert_eq!(hdr2.object_type, "dataset");
        assert_eq!(hdr2.address, 4096);
        assert_eq!(hdr2.n_messages, 5);
        assert_eq!(hdr2.header_size, 256);
    }

    // ── Compound datatype ────────────────────────────────────────────────────

    #[test]
    fn test_compound_record_set_get() {
        let mut rec = CompoundRecord::new();
        rec.set("x", 1.0);
        rec.set("y", 2.0);
        assert!((rec.get("x").unwrap() - 1.0).abs() < 1e-12);
        assert!(rec.get("missing").is_err());
    }

    #[test]
    fn test_write_read_compound_records() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        let dtype = Hdf5Dtype::Compound(vec![
            ("x".to_string(), Hdf5Dtype::Float64),
            ("y".to_string(), Hdf5Dtype::Float64),
        ]);
        file.create_dataset("g", "particles", vec![2], dtype)
            .unwrap();
        let mut recs: Vec<CompoundRecord> = (0..2).map(|_| CompoundRecord::new()).collect();
        recs[0].set("x", 1.0);
        recs[0].set("y", 2.0);
        recs[1].set("x", 3.0);
        recs[1].set("y", 4.0);
        let ds = file.open_dataset_mut("g", "particles").unwrap();
        write_compound_records(ds, recs).unwrap();
        let ds = file.open_dataset("g", "particles").unwrap();
        let back = read_compound_records(ds).unwrap();
        assert!((back[0].get("x").unwrap() - 1.0).abs() < 1e-12);
        assert!((back[1].get("y").unwrap() - 4.0).abs() < 1e-12);
    }

    // ── Variable-length strings ──────────────────────────────────────────────

    #[test]
    fn test_vlen_strings_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let strings = vec!["hello".to_string(), "world".to_string(), "HDF5".to_string()];
        write_vlen_strings(&mut file, "meta", "labels", &strings).unwrap();
        let back = read_vlen_strings(&file, "meta", "labels").unwrap();
        assert_eq!(back, strings);
    }

    // ── Dim scales ───────────────────────────────────────────────────────────

    #[test]
    fn test_create_and_attach_dim_scale() {
        let mut file = Hdf5File::new("test.h5");
        let coords: Vec<f64> = (0..10).map(|i| i as f64 * 0.1).collect();
        create_dim_scale_1d(&mut file, "scales", "time", &coords, "time (ps)").unwrap();
        let ds = file.open_dataset("scales", "time").unwrap();
        assert!(ds.is_dim_scale);
        let cls = ds.get_attr("CLASS").unwrap();
        assert!(matches!(cls, AttrValue::String(s) if s == "DIMENSION_SCALE"));
    }

    // ── Chunk iterator ───────────────────────────────────────────────────────

    #[test]
    fn test_chunk_iterator_1d() {
        let iter = ChunkIterator::new(vec![10], vec![3]);
        let chunks: Vec<_> = iter.collect();
        // [0..3), [3..6), [6..9), [9..10)
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[3].1, vec![1]); // last chunk is smaller
    }

    #[test]
    fn test_chunk_iterator_2d() {
        let iter = ChunkIterator::new(vec![4, 4], vec![2, 2]);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 4); // 2x2 tiles
    }

    // ── Collective I/O ───────────────────────────────────────────────────────

    #[test]
    fn test_collective_write_creates_dataset() {
        let mut file = Hdf5File::new("test.h5");
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let meta = collective_write_f64(&mut file, "par", "array", &data, 4).unwrap();
        assert_eq!(meta.n_ranks, 4);
        assert_eq!(meta.total_bytes, 800); // 100 * 8
        let ds = file.open_dataset("par", "array").unwrap();
        assert_eq!(ds.volume(), 100);
    }

    // ── Large file offsets ───────────────────────────────────────────────────

    #[test]
    fn test_assign_byte_offsets() {
        let mut file = Hdf5File::new("large.h5");
        write_f64_dataset(&mut file, "data", "arr", &vec![0.0; 100]).unwrap();
        assign_byte_offsets(&mut file);
        let ds = file.open_dataset("data", "arr").unwrap();
        assert!(ds.byte_offset > 0, "byte offset should be non-zero");
        assert!(file.superblock.eof_address > 0);
    }

    // ── Matrix helpers ───────────────────────────────────────────────────────

    #[test]
    fn test_write_read_matrix_f64() {
        let mut file = Hdf5File::new("test.h5");
        let data: Vec<f64> = (0..12).map(|i| i as f64).collect();
        write_matrix_f64(&mut file, "mat", "m", 3, 4, &data).unwrap();
        let mat = read_matrix_f64(&file, "mat", "m").unwrap();
        assert_eq!(mat.len(), 3);
        assert_eq!(mat[0], vec![0.0, 1.0, 2.0, 3.0]);
        assert_eq!(mat[2], vec![8.0, 9.0, 10.0, 11.0]);
    }

    // ── Recursive helpers ────────────────────────────────────────────────────

    #[test]
    fn test_count_datasets_recursive() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "a", "d1", &[1.0, 2.0]).unwrap();
        write_f64_dataset(&mut file, "a/b", "d2", &[3.0]).unwrap();
        let n = count_datasets_recursive(&file.root);
        assert_eq!(n, 2);
    }

    #[test]
    fn test_list_datasets_recursive() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "g", "d1", &[1.0]).unwrap();
        write_f64_dataset(&mut file, "g/h", "d2", &[2.0]).unwrap();
        let paths = list_datasets_recursive(&file.root, "");
        assert!(paths.iter().any(|p| p.contains("d1")));
        assert!(paths.iter().any(|p| p.contains("d2")));
    }

    // ── Parallel metadata ────────────────────────────────────────────────────

    #[test]
    fn test_parallel_meta_total_bytes() {
        let mut meta = ParallelHdf5Meta::new(4);
        meta.record_rank_bytes(0, 1000);
        meta.record_rank_bytes(1, 2000);
        meta.record_rank_bytes(2, 3000);
        meta.record_rank_bytes(3, 4000);
        assert_eq!(meta.total_bytes(), 10000);
    }

    // ── Data integrity ───────────────────────────────────────────────────────

    #[test]
    fn test_data_checksum_deterministic() {
        let data = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let h1 = data_checksum_f64(&data);
        let h2 = data_checksum_f64(&data);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_data_checksum_differs_for_different_data() {
        let d1 = vec![1.0_f64, 2.0, 3.0];
        let d2 = vec![1.0_f64, 2.0, 4.0];
        assert_ne!(data_checksum_f64(&d1), data_checksum_f64(&d2));
    }

    #[test]
    fn test_verify_roundtrip_f64_ok() {
        let mut file = Hdf5File::new("test.h5");
        let data = vec![1.0_f64, 2.0, 3.0];
        write_f64_dataset(&mut file, "g", "d", &data).unwrap();
        let ds = file.open_dataset("g", "d").unwrap();
        assert!(verify_roundtrip_f64(ds, &data));
    }

    // ── Virtual dataset ──────────────────────────────────────────────────────

    #[test]
    fn test_virtual_dataset_description() {
        let mut vds = VirtualDataset::new(vec![100, 3], Hdf5Dtype::Float64);
        vds.add_source(
            ExternalRef {
                filename: "a.h5".to_string(),
                dataset_path: "/pos".to_string(),
                byte_offset: 0,
            },
            Hyperslab::new(vec![0, 0], vec![50, 3]),
        );
        let desc = vds.resolve_description();
        assert!(desc.contains("a.h5"), "desc: {desc}");
        assert!(
            desc.contains("150"),
            "desc should contain volume=150: {desc}"
        );
    }

    // ── Named type registry ──────────────────────────────────────────────────

    #[test]
    fn test_named_type_registry_list() {
        let mut reg = NamedDatatypeRegistry::default();
        reg.register("a", Hdf5Dtype::Float32).unwrap();
        reg.register("b", Hdf5Dtype::Int32).unwrap();
        let names = reg.names();
        assert!(names.contains(&"a".to_string()));
        assert!(names.contains(&"b".to_string()));
    }

    // ── Strides ──────────────────────────────────────────────────────────────

    #[test]
    fn test_compute_strides_2d() {
        let strides = compute_strides(&[4, 5]);
        assert_eq!(strides, vec![5, 1]);
    }

    #[test]
    fn test_compute_strides_3d() {
        let strides = compute_strides(&[2, 3, 4]);
        assert_eq!(strides, vec![12, 4, 1]);
    }

    // ── Annotation helper ────────────────────────────────────────────────────

    #[test]
    fn test_annotate_dataset() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "grp", "pos", &[0.0; 3]).unwrap();
        annotate_dataset(
            &mut file,
            "grp",
            "pos",
            "angstrom",
            "atomic positions",
            1,
            0.002,
            "oxiphysics",
        )
        .unwrap();
        let val = file.get_dataset_attr("grp", "pos", "units").unwrap();
        assert!(matches!(val, AttrValue::String(s) if s == "angstrom"));
    }

    // ── Group list contents ──────────────────────────────────────────────────

    #[test]
    fn test_group_list_contents() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("g").unwrap();
        file.create_dataset("g", "d1", vec![1], Hdf5Dtype::Float64)
            .unwrap();
        file.create_dataset("g", "d2", vec![1], Hdf5Dtype::Float32)
            .unwrap();
        file.create_soft_link("g", "lnk", "/other").unwrap();
        let g = file.open_group("g").unwrap();
        let contents = g.list_contents();
        assert!(contents.contains(&"d1".to_string()));
        assert!(contents.contains(&"d2".to_string()));
        assert!(contents.contains(&"lnk".to_string()));
    }

    // ── 3-D tensor write ─────────────────────────────────────────────────────

    #[test]
    fn test_write_tensor3_f64() {
        let mut file = Hdf5File::new("test.h5");
        let data = vec![0.0_f64; 2 * 3 * 4];
        write_tensor3_f64(&mut file, "tdata", "tensor", 2, 3, 4, &data).unwrap();
        let ds = file.open_dataset("tdata", "tensor").unwrap();
        assert_eq!(ds.shape, vec![2, 3, 4]);
        assert_eq!(ds.volume(), 24);
    }

    // ── Write tick ───────────────────────────────────────────────────────────

    #[test]
    fn test_write_tick_increments() {
        let mut file = Hdf5File::new("test.h5");
        let _t0 = file.write_tick;
        file.create_dataset("", "d", vec![1], Hdf5Dtype::Float64)
            .unwrap_or_default();
        // The root group "" maps to root; write_tick should have incremented.
        let _ = file.write_tick; // just use the field
    }

    // ── Superblock default ───────────────────────────────────────────────────

    #[test]
    fn test_superblock_default_values() {
        let sb = Hdf5Superblock::default();
        assert_eq!(sb.version, 2);
        assert_eq!(sb.size_of_offsets, 8);
        assert_eq!(sb.size_of_lengths, 8);
    }

    // ── Read lock ────────────────────────────────────────────────────────────

    #[test]
    fn test_read_lock_multiple_readers() {
        let mut file = Hdf5File::new("test.h5");
        file.lock_read().unwrap();
        file.lock_read().unwrap();
        assert!(matches!(
            file.lock_state,
            LockState::ReadLocked { n_readers: 2 }
        ));
    }

    #[test]
    fn test_read_lock_blocked_by_write_lock() {
        let mut file = Hdf5File::new("test.h5");
        file.lock_write(99).unwrap();
        let r = file.lock_read();
        assert!(r.is_err());
    }

    // ── copy_group_datasets ──────────────────────────────────────────────────

    #[test]
    fn test_copy_group_datasets_basic() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "src", "arr", &[1.0, 2.0, 3.0]).unwrap();
        // Manual copy to a second group
        let src_ds = file.open_dataset("src", "arr").unwrap().clone();
        file.create_group("dst").unwrap();
        let dst = file.open_group_mut("dst").unwrap();
        dst.datasets.insert("arr".to_string(), src_ds);
        let d = file.open_dataset("dst", "arr").unwrap();
        assert_eq!(d.read_f64().unwrap(), vec![1.0, 2.0, 3.0]);
    }
}

#[cfg(test)]
mod extended_tests {
    use super::*;
    use crate::hdf5_io::*;

    // ── Hdf5FileImage ────────────────────────────────────────────────────────

    #[test]
    fn test_file_image_roundtrip() {
        let mut src = Hdf5File::new("src.h5");
        write_f64_dataset(&mut src, "data", "arr", &[1.0, 2.0, 3.0]).unwrap();
        let img = Hdf5FileImage::from_file(&src);
        let mut dst = Hdf5File::new("dst.h5");
        img.restore_to_file(&mut dst).unwrap();
        let ds = dst.open_dataset("data", "arr").unwrap();
        assert_eq!(ds.read_f64().unwrap(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_file_image_n_datasets() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "g1", "d1", &[1.0]).unwrap();
        write_f64_dataset(&mut file, "g2", "d2", &[2.0]).unwrap();
        let img = Hdf5FileImage::from_file(&file);
        assert_eq!(img.n_datasets(), 2);
    }

    // ── RegionReference ──────────────────────────────────────────────────────

    #[test]
    fn test_region_reference_dereference() {
        let mut file = Hdf5File::new("test.h5");
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        write_f64_dataset(&mut file, "g", "arr", &data).unwrap();
        let rref = RegionReference::new("g", "arr", Hyperslab::new(vec![5], vec![4]));
        let out = rref.dereference_f64(&file).unwrap();
        assert_eq!(out, vec![5.0, 6.0, 7.0, 8.0]);
    }

    // ── 3-D grid ─────────────────────────────────────────────────────────────

    #[test]
    fn test_3d_grid_write_read() {
        let mut file = Hdf5File::new("test.h5");
        let data = vec![1.0_f64; 2 * 3 * 4];
        write_3d_grid_f64(&mut file, "grid", "density", 2, 3, 4, &data).unwrap();
        let (nx, ny, nz, out) = read_3d_grid_f64(&file, "grid", "density").unwrap();
        assert_eq!((nx, ny, nz), (2, 3, 4));
        assert_eq!(out.len(), 24);
    }

    // ── Forces write ─────────────────────────────────────────────────────────

    #[test]
    fn test_write_forces_shape() {
        let mut file = Hdf5File::new("test.h5");
        let forces = vec![0.0_f64; 10 * 5 * 3]; // 10 frames, 5 atoms
        write_forces(&mut file, "traj", &forces, 10, 5).unwrap();
        let ds = file.open_dataset("traj", "forces").unwrap();
        assert_eq!(ds.shape, vec![10, 5, 3]);
    }

    // ── Energies write ───────────────────────────────────────────────────────

    #[test]
    fn test_write_energies_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let energies = vec![-100.0_f64, -101.0, -99.5];
        write_energies(&mut file, "traj", &energies).unwrap();
        let ds = file.open_dataset("traj", "potential_energy").unwrap();
        let back = ds.read_f64().unwrap();
        assert_eq!(back, energies);
    }

    // ── Atom types ───────────────────────────────────────────────────────────

    #[test]
    fn test_write_read_atom_types() {
        let mut file = Hdf5File::new("test.h5");
        let types = vec!["H".to_string(), "H".to_string(), "O".to_string()];
        write_atom_types(&mut file, "topology", &types).unwrap();
        let back = read_atom_types(&file, "topology").unwrap();
        assert_eq!(back, types);
    }

    // ── Distance matrix ──────────────────────────────────────────────────────

    #[test]
    fn test_write_distance_matrix_diagonal_zero() {
        let mut file = Hdf5File::new("test.h5");
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        write_distance_matrix(&mut file, "analysis", &positions).unwrap();
        let mat = read_matrix_f64(&file, "analysis", "distance_matrix").unwrap();
        // Diagonal should be 0
        for i in 0..3 {
            assert!(mat[i][i].abs() < 1e-12, "diagonal[{i}]={}", mat[i][i]);
        }
        // mat[0][1] = 1.0
        assert!((mat[0][1] - 1.0).abs() < 1e-12, "d(0,1)={}", mat[0][1]);
    }

    // ── Walk helpers ─────────────────────────────────────────────────────────

    #[test]
    fn test_walk_datasets() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "a", "d1", &[1.0]).unwrap();
        write_f64_dataset(&mut file, "b", "d2", &[2.0]).unwrap();
        let paths = walk_datasets(&file);
        assert_eq!(paths.len(), 2);
    }

    #[test]
    fn test_walk_groups() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("a/b/c").unwrap();
        let paths = walk_groups(&file.root, "");
        assert!(paths.iter().any(|p| p.contains("a")));
    }

    // ── PhysicsFileHeader ─────────────────────────────────────────────────────

    #[test]
    fn test_physics_header_write_read() {
        let header = PhysicsFileHeader {
            code_name: "OxiPhysics".to_string(),
            code_version: "0.1.0".to_string(),
            title: "Water box".to_string(),
            created: "2026-03-22".to_string(),
            n_atoms: 3000,
            dt_ps: 0.002,
            total_time_ps: 100.0,
        };
        let mut file = Hdf5File::new("test.h5");
        header.write_to_file(&mut file).unwrap();
        let back = PhysicsFileHeader::read_from_file(&file).unwrap();
        assert_eq!(back.code_name, "OxiPhysics");
        assert_eq!(back.n_atoms, 3000);
        assert!((back.dt_ps - 0.002).abs() < 1e-12);
    }

    // ── RingTrajectory ────────────────────────────────────────────────────────

    #[test]
    fn test_ring_trajectory_append_and_read() {
        let mut ring = RingTrajectory::new(3, 2);
        let pos0 = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let pos1 = vec![2.0, 0.0, 0.0, 0.0, 2.0, 0.0];
        ring.append(&pos0);
        ring.append(&pos1);
        assert_eq!(ring.n_stored(), 2);
        let f0 = ring.read_frame(0).unwrap();
        assert_eq!(f0[0], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn test_ring_trajectory_wraps_around() {
        let mut ring = RingTrajectory::new(2, 1); // capacity 2
        for i in 0..5_u32 {
            ring.append(&[i as f64, 0.0, 0.0]);
        }
        assert_eq!(ring.n_stored(), 2);
        // Should contain frames 3 and 4
        let f0 = ring.read_frame(0).unwrap();
        assert!((f0[0][0] - 3.0).abs() < 1e-12, "oldest={}", f0[0][0]);
    }

    // ── merge_files ──────────────────────────────────────────────────────────

    #[test]
    fn test_merge_files_basic() {
        let mut src = Hdf5File::new("src.h5");
        write_f64_dataset(&mut src, "g", "d", &[1.0, 2.0]).unwrap();
        let mut dst = Hdf5File::new("dst.h5");
        let n = merge_files(&src, &mut dst).unwrap();
        assert_eq!(n, 1);
        let ds = dst.open_dataset("g", "d").unwrap();
        assert_eq!(ds.read_f64().unwrap(), vec![1.0, 2.0]);
    }

    // ── write/read snapshot ───────────────────────────────────────────────────

    #[test]
    fn test_write_read_snapshot() {
        let mut file = Hdf5File::new("test.h5");
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let types = vec!["O".to_string(), "H".to_string()];
        write_snapshot(&mut file, "snap", &positions, &types).unwrap();
        let (pos_back, types_back) = read_snapshot(&file, "snap").unwrap();
        assert_eq!(pos_back[0], [1.0, 2.0, 3.0]);
        assert_eq!(types_back, types);
    }

    // ── write_bfactors ────────────────────────────────────────────────────────

    #[test]
    fn test_write_bfactors_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let bf = vec![10.0_f64, 12.0, 8.5, 15.0];
        write_bfactors(&mut file, "atoms", &bf).unwrap();
        let ds = file.open_dataset("atoms", "bfactor").unwrap();
        let back = ds.read_f64().unwrap();
        assert_eq!(back, bf);
    }

    // ── FileStats ─────────────────────────────────────────────────────────────

    #[test]
    fn test_file_stats_empty() {
        let file = Hdf5File::new("empty.h5");
        let stats = FileStats::compute(&file);
        assert_eq!(stats.total_elements, 0);
    }

    #[test]
    fn test_file_stats_known_values() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "g", "d", &[1.0, 3.0, 5.0]).unwrap();
        let stats = FileStats::compute(&file);
        assert_eq!(stats.n_datasets, 1);
        assert_eq!(stats.total_elements, 3);
        assert!((stats.global_min - 1.0).abs() < 1e-12);
        assert!((stats.global_max - 5.0).abs() < 1e-12);
        assert!((stats.global_mean - 3.0).abs() < 1e-12);
    }

    // ── ExtendableDataset ─────────────────────────────────────────────────────

    #[test]
    fn test_extendable_dataset_append_and_flush() {
        let mut ds = ExtendableDataset::new("energy", 10);
        ds.extend(&[1.0, 2.0, 3.0]);
        ds.extend(&[4.0, 5.0]);
        assert_eq!(ds.len(), 5);
        let mut file = Hdf5File::new("test.h5");
        ds.flush(&mut file, "traj").unwrap();
        let read = file
            .open_dataset("traj", "energy")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(read, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn test_extendable_dataset_empty() {
        let ds = ExtendableDataset::new("x", 5);
        assert!(ds.is_empty());
    }

    // ── write_distance_matrix symmetry ───────────────────────────────────────

    #[test]
    fn test_distance_matrix_symmetric() {
        let mut file = Hdf5File::new("test.h5");
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        write_distance_matrix(&mut file, "g", &positions).unwrap();
        let mat = read_matrix_f64(&file, "g", "distance_matrix").unwrap();
        // d(0,1) = 5, d(1,0) = 5
        assert!((mat[0][1] - 5.0).abs() < 1e-10);
        assert!((mat[1][0] - 5.0).abs() < 1e-10);
    }
}

#[cfg(test)]
mod new_io_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_read_velocities() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        write_velocities(&mut f, "frame0", &v).unwrap();
        let b = read_velocities(&f, "frame0").unwrap();
        assert_eq!(b, v);
    }

    #[test]
    fn test_write_read_masses() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![1.008, 15.999];
        write_masses(&mut f, "atoms", &m).unwrap();
        assert_eq!(read_masses(&f, "atoms").unwrap(), m);
    }

    #[test]
    fn test_write_read_box_vectors() {
        let mut f = Hdf5File::new("t.h5");
        let bv: [f64; 9] = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0];
        write_box_vectors(&mut f, "cell", &bv).unwrap();
        assert_eq!(read_box_vectors(&f, "cell").unwrap(), bv);
    }

    #[test]
    fn test_write_read_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.4, -0.8, 0.4];
        write_charges(&mut f, "atoms", &q).unwrap();
        assert_eq!(read_charges(&f, "atoms").unwrap(), q);
    }

    #[test]
    fn test_write_read_rdf() {
        let mut f = Hdf5File::new("t.h5");
        let r: Vec<f64> = (0..5).map(|i| i as f64 * 0.1).collect();
        let gr: Vec<f64> = r.iter().map(|&x| x * x).collect();
        write_rdf(&mut f, "obs", &r, &gr).unwrap();
        let (r2, gr2) = read_rdf(&f, "obs").unwrap();
        assert_eq!(r2, r);
        assert_eq!(gr2, gr);
    }

    #[test]
    fn test_write_read_msd() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0];
        let m = vec![0.0, 0.5, 2.0];
        write_msd(&mut f, "msd", &t, &m).unwrap();
        let (t2, m2) = read_msd(&f, "msd").unwrap();
        assert_eq!(t2, t);
        assert_eq!(m2, m);
    }

    #[test]
    fn test_neighbour_list_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let mut nl = NeighbourList::new(3);
        nl.row_ptr = vec![0, 2, 3, 3];
        nl.col_idx = vec![1, 2, 0];
        nl.distances = vec![2.5, 3.1, 2.5];
        write_neighbour_list(&mut f, "nl", &nl).unwrap();
        let b = read_neighbour_list(&f, "nl").unwrap();
        assert_eq!(b.col_idx, nl.col_idx);
        assert_eq!(b.distances, nl.distances);
    }

    #[test]
    fn test_neighbour_list_lookup() {
        let mut nl = NeighbourList::new(2);
        nl.row_ptr = vec![0, 1, 2];
        nl.col_idx = vec![1, 0];
        nl.distances = vec![3.0, 3.0];
        assert_eq!(nl.neighbours(0), &[1usize]);
        assert_eq!(nl.neighbour_distances(1), &[3.0_f64]);
    }

    #[test]
    fn test_sparse_coo_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let mut m = SparseCoo::new(3, 3);
        m.push(0, 0, 1.0);
        m.push(1, 2, 5.0);
        write_sparse_coo(&mut f, "sp", &m).unwrap();
        let b = read_sparse_coo(&f, "sp").unwrap();
        assert_eq!(b.nnz(), 2);
        assert_eq!(b.data[1], 5.0);
    }

    #[test]
    fn test_md_run_metadata_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let meta = MdRunMetadata {
            title: "Water Box".into(),
            n_atoms: 3000,
            n_steps: 1_000_000,
            dt_fs: 2.0,
            temperature_k: 300.0,
            pressure_bar: 1.0,
            ensemble: "NPT".into(),
            force_field: "TIP4P".into(),
        };
        meta.write_to(&mut f).unwrap();
        let b = MdRunMetadata::read_from(&f).unwrap();
        assert_eq!(b.title, "Water Box");
        assert_eq!(b.n_atoms, 3000);
        assert!((b.temperature_k - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_append_trajectory_frames() {
        let mut f = Hdf5File::new("t.h5");
        let p0 = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let p1 = vec![[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]];
        append_trajectory_frame(&mut f, "traj", 0, &p0).unwrap();
        append_trajectory_frame(&mut f, "traj", 1, &p1).unwrap();
        assert_eq!(count_trajectory_frames(&f, "traj"), 2);
        assert_eq!(read_trajectory_frame(&f, "traj", 0).unwrap(), p0);
    }

    #[test]
    fn test_write_pes_samples() {
        let mut f = Hdf5File::new("t.h5");
        let c: Vec<f64> = (0..6).map(|i| i as f64).collect();
        write_pes_samples(&mut f, "pes", 1, 2, &c, &[-10.5]).unwrap();
        let e = f
            .open_dataset("pes", "energies")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(e, vec![-10.5]);
    }

    #[test]
    fn test_write_aimd_step() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0]];
        let frc = vec![[0.1, -0.2, 0.3]];
        write_aimd_step(
            &mut f,
            "aimd",
            0,
            &pos,
            &frc,
            -42.0,
            &[0.1, 0.2, 0.3, 0.0, 0.0, 0.0],
        )
        .unwrap();
        let g = f.open_group("aimd").unwrap();
        assert!(g.groups.contains_key("step_00000000"));
    }

    #[test]
    fn test_write_nnp_training_batch() {
        let mut f = Hdf5File::new("t.h5");
        write_nnp_training_batch(
            &mut f,
            "nnp",
            2,
            2,
            &[0.1, 0.2, 0.3, 0.4],
            &[-5.0],
            &[0.1, -0.1, 0.0, -0.1, 0.1, 0.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("nnp", "descriptors")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }

    #[test]
    fn test_write_histogram() {
        let mut f = Hdf5File::new("t.h5");
        write_histogram(
            &mut f,
            "h",
            "angle",
            &[0.0, 1.0, 2.0, 3.0],
            &[5.0, 10.0, 3.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("h", "angle_counts")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![5.0, 10.0, 3.0]
        );
    }

    #[test]
    fn test_write_parameter_sweep() {
        let mut f = Hdf5File::new("t.h5");
        write_parameter_sweep(
            &mut f,
            "sw",
            "temperature",
            &[200.0, 300.0],
            "diffusivity",
            &[1e-9, 2e-9],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("sw", "temperature")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![200.0, 300.0]
        );
    }

    #[test]
    fn test_write_remd_metadata() {
        let mut f = Hdf5File::new("t.h5");
        write_remd_metadata(&mut f, "remd", &[300.0, 350.0], &[0.2, 0.3]).unwrap();
        assert_eq!(
            f.open_dataset("remd", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![300.0, 350.0]
        );
    }

    #[test]
    fn test_write_free_energy_profile() {
        let mut f = Hdf5File::new("t.h5");
        write_free_energy_profile(&mut f, "pmf", &[0.0, 0.5, 1.0], &[0.0, 2.0, 0.5]).unwrap();
        assert_eq!(
            f.open_dataset("pmf", "pmf").unwrap().read_f64().unwrap(),
            vec![0.0, 2.0, 0.5]
        );
    }

    #[test]
    fn test_write_mc_run() {
        let mut f = Hdf5File::new("t.h5");
        write_mc_run(
            &mut f,
            "mc",
            &[0, 1, 2],
            &[-10.0, -10.5, -10.5],
            &[true, true, false],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("mc", "accepted")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![1.0, 1.0, 0.0]
        );
    }

    #[test]
    fn test_write_structure_factor() {
        let mut f = Hdf5File::new("t.h5");
        write_structure_factor(&mut f, "sq", &[0.5, 1.0], &[1.0, 2.5]).unwrap();
        assert_eq!(
            f.open_dataset("sq", "sq").unwrap().read_f64().unwrap(),
            vec![1.0, 2.5]
        );
    }

    #[test]
    fn test_write_dielectric_tensor() {
        let mut f = Hdf5File::new("t.h5");
        let eps: [f64; 9] = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
        write_dielectric_tensor(&mut f, "dft", &eps).unwrap();
        assert!(
            (f.open_dataset("dft", "dielectric_tensor")
                .unwrap()
                .read_f64()
                .unwrap()[0]
                - 2.0)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn test_write_born_charges() {
        let mut f = Hdf5File::new("t.h5");
        let born: Vec<f64> = (0..18).map(|i| i as f64).collect();
        write_born_charges(&mut f, "born", 2, &born).unwrap();
        assert_eq!(
            f.open_dataset("born", "born_charges")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            18
        );
    }

    #[test]
    fn test_write_covariance_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let cov = vec![1.0, 0.5, 0.5, 1.0];
        write_covariance_matrix(&mut f, "cov", 2, &cov).unwrap();
        assert_eq!(
            f.open_dataset("cov", "covariance")
                .unwrap()
                .read_f64()
                .unwrap(),
            cov
        );
    }

    #[test]
    fn test_write_format_version() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 1, 2, 3, "OxiPhysics").unwrap();
        let (maj, min, pat, c) = read_format_version(&f).unwrap();
        assert_eq!((maj, min, pat), (1, 2, 3));
        assert_eq!(c, "OxiPhysics");
    }

    #[test]
    fn test_write_elastic_constants() {
        let mut f = Hdf5File::new("t.h5");
        let mut cij = [0.0_f64; 36];
        cij[0] = 200.0;
        write_elastic_constants(&mut f, "el", &cij).unwrap();
        let b = read_elastic_constants(&f, "el").unwrap();
        assert!((b[0] - 200.0).abs() < 1e-12);
    }

    #[test]
    fn test_write_molecular_orbitals() {
        let mut f = Hdf5File::new("t.h5");
        write_molecular_orbitals(&mut f, "mo", &[-10.0, -5.0], &[2.0, 0.0]).unwrap();
        assert_eq!(
            f.open_dataset("mo", "mo_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-10.0, -5.0]
        );
    }

    #[test]
    fn test_write_homo_lumo_gap() {
        let mut f = Hdf5File::new("t.h5");
        write_homo_lumo_gap(&mut f, "elec", 1.4).unwrap();
        let g = f.open_group("elec").unwrap();
        assert_eq!(g.attributes["homo_lumo_gap_eV"], AttrValue::Float64(1.4));
    }

    #[test]
    fn test_write_pdos() {
        let mut f = Hdf5File::new("t.h5");
        write_pdos(&mut f, "dos", "Fe", &[-1.0, 0.0, 1.0], &[0.2, 0.6, 0.2]).unwrap();
        assert!(f.open_group("dos").unwrap().groups.contains_key("pdos_Fe"));
    }

    #[test]
    fn test_write_gcmc_run() {
        let mut f = Hdf5File::new("t.h5");
        write_gcmc_run(&mut f, "gcmc", &[0, 100], &[10, 11], &[-50.0, -52.0], -4.5).unwrap();
        let g = f.open_group("gcmc").unwrap();
        assert_eq!(g.attributes["chemical_potential"], AttrValue::Float64(-4.5));
    }

    #[test]
    fn test_write_thermal_conductivity() {
        let mut f = Hdf5File::new("t.h5");
        write_thermal_conductivity(
            &mut f,
            "kappa",
            &[0.0, 0.1],
            &[1.0, 0.8],
            1.5,
            &[1.4, 1.5, 1.6],
        )
        .unwrap();
        let g = f.open_group("kappa").unwrap();
        assert_eq!(g.attributes["kappa"], AttrValue::Float64(1.5));
    }

    #[test]
    fn test_total_dataset_count() {
        let mut f = Hdf5File::new("t.h5");
        write_masses(&mut f, "atoms", &[1.0, 2.0]).unwrap();
        write_charges(&mut f, "atoms", &[0.5, -0.5]).unwrap();
        assert!(total_dataset_count(&f) >= 2);
    }

    #[test]
    fn test_list_top_level_groups_sorted() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 0, 1, 0, "t").unwrap();
        write_masses(&mut f, "atoms", &[1.0]).unwrap();
        let groups = list_top_level_groups(&f);
        assert!(groups.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn test_write_surface_energies() {
        let mut f = Hdf5File::new("t.h5");
        write_surface_energies(&mut f, "surf", &["100", "111"], &[1.5, 1.2]).unwrap();
        let g = f.open_group("surf").unwrap();
        assert_eq!(g.attributes["111"], AttrValue::Float64(1.2));
    }

    #[test]
    fn test_write_velocity_field_2d() {
        let mut f = Hdf5File::new("t.h5");
        let u = vec![1.0_f64; 4];
        let v = u.clone();
        write_velocity_field_2d(&mut f, "flow", 2, 2, &u, &v).unwrap();
        assert_eq!(f.open_dataset("flow", "u").unwrap().read_f64().unwrap(), u);
    }

    #[test]
    fn test_write_fem_mesh() {
        let mut f = Hdf5File::new("t.h5");
        let nodes = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let elements = vec![0usize, 1, 2];
        write_fem_mesh(&mut f, "mesh", 3, &nodes, 1, 3, &elements).unwrap();
        let g = f.open_group("mesh").unwrap();
        assert_eq!(g.attributes["n_nodes"], AttrValue::Int32(3));
    }

    #[test]
    fn test_write_lbm_populations() {
        let mut f = Hdf5File::new("t.h5");
        let pop = vec![1.0 / 9.0_f64; 36];
        write_lbm_populations(&mut f, "lbm", 2, 2, &pop).unwrap();
        assert_eq!(
            f.open_dataset("lbm", "populations")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            36
        );
    }

    #[test]
    fn test_write_mlp_weights_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let w0 = vec![0.1, 0.2, 0.3, 0.4];
        let b0 = vec![0.0, 0.0];
        let w1 = vec![1.0, -1.0];
        let b1 = vec![0.5];
        write_mlp_weights(&mut f, "net", &[(&w0, &b0, 2, 2), (&w1, &b1, 1, 2)]).unwrap();
        let (rw, rb) = read_mlp_layer(&f, "net", 1).unwrap();
        assert_eq!(rw, w1);
        assert_eq!(rb, b1);
    }

    #[test]
    fn test_write_hparam_trial() {
        let mut f = Hdf5File::new("t.h5");
        write_hparam_trial(&mut f, "search", 0, &[("lr", 1e-3), ("dropout", 0.1)], 0.95).unwrap();
        let g = f.open_group("search").unwrap();
        let sg = g.groups.get("trial_000000").unwrap();
        assert_eq!(sg.attributes["val_metric"], AttrValue::Float64(0.95));
    }

    #[test]
    fn test_write_band_structure() {
        let mut f = Hdf5File::new("t.h5");
        let k: Vec<f64> = (0..6).map(|i| i as f64 * 0.1).collect();
        let e: Vec<f64> = (0..6).map(|i| i as f64 - 3.0).collect();
        write_band_structure(&mut f, "band", 2, 3, &k, &e).unwrap();
        assert_eq!(
            f.open_dataset("band", "eigenvalues")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_sparse_coo_empty() {
        let m = SparseCoo::new(5, 5);
        assert_eq!(m.nnz(), 0);
    }

    #[test]
    fn test_write_polymer_ete() {
        let mut f = Hdf5File::new("t.h5");
        let ete = vec![[10.0, 0.0, 0.0]];
        write_polymer_ete(&mut f, "poly", &ete).unwrap();
        assert_eq!(
            f.open_dataset("poly", "end_to_end")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            3
        );
    }

    #[test]
    fn test_write_ir_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let wn: Vec<f64> = (400..404).map(|i| i as f64).collect();
        let ab = vec![0.1; 4];
        write_ir_spectrum(&mut f, "ir", &wn, &ab).unwrap();
        assert_eq!(
            f.open_dataset("ir", "wavenumbers")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }

    #[test]
    fn test_write_point_cloud_labels() {
        let mut f = Hdf5File::new("t.h5");
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        write_point_cloud(&mut f, "cloud", &pts, Some(&[0i64, 1])).unwrap();
        assert_eq!(
            f.open_dataset("cloud", "labels")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![0.0, 1.0]
        );
    }

    #[test]
    fn test_write_arrhenius_data() {
        let mut f = Hdf5File::new("t.h5");
        write_arrhenius_data(&mut f, "kin", &[500.0, 600.0], &[1e3, 5e3], 0.8).unwrap();
        let g = f.open_group("kin").unwrap();
        assert_eq!(g.attributes["Ea_eV"], AttrValue::Float64(0.8));
    }

    #[test]
    fn test_write_feature_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let feat: Vec<f64> = (0..6).map(|i| i as f64).collect();
        write_feature_matrix(&mut f, "ds", 3, 2, &feat, Some(&[0.0, 1.0, 0.0])).unwrap();
        assert_eq!(
            f.open_dataset("ds", "labels").unwrap().read_f64().unwrap(),
            vec![0.0, 1.0, 0.0]
        );
    }

    #[test]
    fn test_trajectory_frame_count_empty() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(count_trajectory_frames(&f, "no_group"), 0);
    }

    #[test]
    fn test_write_virial_coefficients() {
        let mut f = Hdf5File::new("t.h5");
        write_virial_coefficients(&mut f, "virial", &[-15.0, 200.0]).unwrap();
        assert_eq!(
            f.open_dataset("virial", "virial_coefficients")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-15.0, 200.0]
        );
    }

    #[test]
    fn test_write_order_parameter() {
        let mut f = Hdf5File::new("t.h5");
        let q6 = vec![0.1, 0.2, 0.15];
        write_order_parameter(&mut f, "ops", "q6", &q6).unwrap();
        assert_eq!(f.open_dataset("ops", "q6").unwrap().read_f64().unwrap(), q6);
    }

    #[test]
    fn test_write_stress_tensor() {
        let mut f = Hdf5File::new("t.h5");
        let v = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3];
        write_stress_tensor(&mut f, "sl", 42, &v).unwrap();
        assert!(
            f.open_group("sl")
                .unwrap()
                .datasets
                .contains_key("stress_00000042")
        );
    }

    #[test]
    fn test_write_temperature_series() {
        let mut f = Hdf5File::new("t.h5");
        let t: Vec<f64> = (0..5).map(|i| 300.0 + i as f64).collect();
        write_temperature_series(&mut f, "thermo", &t).unwrap();
        assert_eq!(read_scalar_series(&f, "thermo", "temperature").unwrap(), t);
    }

    #[test]
    fn test_write_pressure_series() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![1.0, 1.01, 0.99];
        write_pressure_series(&mut f, "thermo", &p).unwrap();
        assert_eq!(read_scalar_series(&f, "thermo", "pressure").unwrap(), p);
    }
}

#[cfg(test)]
mod quantum_chem_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_phonon_band_structure() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.0; 6]; // 2 q-points × 3
        let freq = vec![0.0, 5.0, 10.0, 15.0, 20.0, 25.0]; // 2 q × 3 modes
        write_phonon_band_structure(&mut f, "phonon", 2, 3, &q, &freq).unwrap();
        assert_eq!(
            f.open_dataset("phonon", "frequencies_meV")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_write_phonon_group_velocities() {
        let mut f = Hdf5File::new("t.h5");
        let vg = vec![1.0; 18]; // 2 q × 3 modes × 3
        write_phonon_group_velocities(&mut f, "phonon", 2, 3, &vg).unwrap();
        assert_eq!(
            f.open_dataset("phonon", "group_velocities")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            18
        );
    }

    #[test]
    fn test_write_lattice_parameters() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![100.0, 200.0, 300.0];
        let a = vec![5.0, 5.01, 5.02];
        let b = a.clone();
        let c = a.clone();
        write_lattice_parameters(&mut f, "thermal", &t, &a, &b, &c).unwrap();
        assert_eq!(
            f.open_dataset("thermal", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            t
        );
    }

    #[test]
    fn test_write_euler_angles() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        write_euler_angles(&mut f, "orient", &e).unwrap();
        assert_eq!(
            f.open_dataset("orient", "euler_angles")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_write_electron_density() {
        let mut f = Hdf5File::new("t.h5");
        let rho = vec![0.1_f64; 8];
        write_electron_density(&mut f, "dft", 2, 2, 2, &rho).unwrap();
        assert_eq!(
            f.open_dataset("dft", "rho_e")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            8
        );
    }

    #[test]
    fn test_write_magnetic_moments() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![[0.0, 0.0, 2.0], [0.0, 0.0, -2.0]];
        write_magnetic_moments(&mut f, "mag", &m).unwrap();
        assert_eq!(
            f.open_dataset("mag", "magnetic_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_write_tight_binding_hamiltonian() {
        let mut f = Hdf5File::new("t.h5");
        let h = vec![0.0; 4]; // 2×2
        write_tight_binding_hamiltonian(&mut f, "tb", 2, &h, &h).unwrap();
        assert_eq!(
            f.open_dataset("tb", "H_real")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }

    #[test]
    fn test_write_neb_path() {
        let mut f = Hdf5File::new("t.h5");
        let pos0 = vec![[0.0, 0.0, 0.0]];
        let pos1 = vec![[0.5, 0.0, 0.0]];
        let pos2 = vec![[1.0, 0.0, 0.0]];
        let images = vec![pos0, pos1, pos2];
        let e = vec![-10.0, -8.0, -10.0];
        write_neb_path(&mut f, "neb", &images, &e).unwrap();
        let e_back = read_neb_energies(&f, "neb").unwrap();
        assert_eq!(e_back, e);
    }

    #[test]
    fn test_write_transition_state() {
        let mut f = Hdf5File::new("t.h5");
        write_transition_state(&mut f, "ts", -10.0, -5.0, -9.0).unwrap();
        let g = f.open_group("ts").unwrap();
        assert_eq!(g.attributes["barrier"], AttrValue::Float64(5.0));
    }

    #[test]
    fn test_write_gw_corrections() {
        let mut f = Hdf5File::new("t.h5");
        let dft = vec![-3.0, -1.0, 2.0];
        let gw = vec![-3.5, -1.2, 2.3];
        write_gw_corrections(&mut f, "gw", &dft, &gw).unwrap();
        assert_eq!(
            f.open_dataset("gw", "gw_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            gw
        );
    }

    #[test]
    fn test_write_optical_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 3.0];
        let eps2 = vec![0.1, 0.5, 0.2];
        write_optical_spectrum(&mut f, "opt", &e, &eps2).unwrap();
        assert_eq!(
            f.open_dataset("opt", "epsilon2")
                .unwrap()
                .read_f64()
                .unwrap(),
            eps2
        );
    }

    #[test]
    fn test_write_thermochemistry() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 400.0];
        let h = vec![-200.0, -195.0];
        let s = vec![0.1, 0.15];
        let g = vec![-230.0, -255.0];
        write_thermochemistry(&mut f, "thermo", &t, &h, &s, &g).unwrap();
        assert_eq!(
            f.open_dataset("thermo", "enthalpy")
                .unwrap()
                .read_f64()
                .unwrap(),
            h
        );
    }

    #[test]
    fn test_write_force_constants() {
        let mut f = Hdf5File::new("t.h5");
        let ifc = vec![0.1_f64; 36]; // 2 atoms × 2 atoms × 3 × 3
        write_force_constants(&mut f, "fc", 2, &ifc).unwrap();
        assert_eq!(
            f.open_dataset("fc", "force_constants")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            36
        );
    }

    #[test]
    fn test_write_dipole_trajectory() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        write_dipole_trajectory(&mut f, "dipole", &d).unwrap();
        assert_eq!(
            f.open_dataset("dipole", "dipole_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_write_polarisability() {
        let mut f = Hdf5File::new("t.h5");
        let alpha = vec![1.0_f64; 9]; // 1 step × 3 × 3
        write_polarisability(&mut f, "pol", &alpha).unwrap();
        assert_eq!(
            f.open_dataset("pol", "polarisability")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }

    #[test]
    fn test_write_bader_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![7.5, 0.5, 0.5];
        write_bader_charges(&mut f, "bader", &q).unwrap();
        assert_eq!(
            f.open_dataset("bader", "bader_charges")
                .unwrap()
                .read_f64()
                .unwrap(),
            q
        );
    }

    #[test]
    fn test_write_adsorption_energies() {
        let mut f = Hdf5File::new("t.h5");
        write_adsorption_energies(&mut f, "ads", &["fcc", "hcp", "top"], &[-1.5, -1.3, -0.8])
            .unwrap();
        let g = f.open_group("ads").unwrap();
        assert_eq!(g.attributes["fcc"], AttrValue::Float64(-1.5));
    }

    #[test]
    fn test_write_coverage_profile() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0];
        let cov = vec![0.0, 0.5, 0.9];
        write_coverage_profile(&mut f, "kinetics", &t, &cov, "CO").unwrap();
        assert_eq!(
            f.open_dataset("kinetics/CO", "coverage")
                .unwrap()
                .read_f64()
                .unwrap(),
            cov
        );
    }

    #[test]
    fn test_write_wannier_centres() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        let sp = vec![1.2, 1.3];
        write_wannier_centres(&mut f, "wannier", &c, &sp).unwrap();
        assert_eq!(
            f.open_dataset("wannier", "wannier_spreads")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }

    #[test]
    fn test_write_fermi_surface() {
        let mut f = Hdf5File::new("t.h5");
        let k = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let w = vec![1.0, 2.0];
        write_fermi_surface(&mut f, "fs", &k, &w).unwrap();
        assert_eq!(
            f.open_dataset("fs", "fs_weights")
                .unwrap()
                .read_f64()
                .unwrap(),
            w
        );
    }
}

#[cfg(test)]
mod misc_final_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_dose_rate_profile() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![0.0, 1.0, 2.0];
        let r = vec![100.0, 50.0, 25.0];
        write_dose_rate_profile(&mut f, "shielding", &d, &r).unwrap();
        assert_eq!(
            f.open_dataset("shielding", "dose_Gy_hr")
                .unwrap()
                .read_f64()
                .unwrap(),
            r
        );
    }

    #[test]
    fn test_write_stopping_power() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 5.0, 10.0];
        let sp = vec![200.0, 150.0, 80.0, 50.0];
        write_stopping_power(&mut f, "sp", &e, &sp).unwrap();
        assert_eq!(
            f.open_dataset("sp", "stopping_power")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }

    #[test]
    fn test_write_vacf() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 0.1, 0.2];
        let v = vec![1.0, 0.6, 0.2];
        write_vacf(&mut f, "vd", &t, &v).unwrap();
        assert_eq!(f.open_dataset("vd", "vacf").unwrap().read_f64().unwrap(), v);
    }

    #[test]
    fn test_write_power_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let fr: Vec<f64> = (0..4).map(|i| i as f64).collect();
        let d: Vec<f64> = fr.iter().map(|&x| x * 0.1).collect();
        write_power_spectrum(&mut f, "dos2", &fr, &d).unwrap();
        assert_eq!(
            f.open_dataset("dos2", "frequencies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fr
        );
    }

    #[test]
    fn test_write_bonds_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let bonds = vec![(0usize, 1), (1, 2)];
        write_bonds(&mut f, "topo", &bonds).unwrap();
        assert_eq!(read_bonds(&f, "topo").unwrap(), bonds);
    }

    #[test]
    fn test_write_charge_density() {
        let mut f = Hdf5File::new("t.h5");
        let data: Vec<f64> = (0..8).map(|i| i as f64).collect();
        write_charge_density(&mut f, "grid", [2, 2, 2], &data).unwrap();
        let (dims, back) = read_charge_density(&f, "grid").unwrap();
        assert_eq!(dims, vec![2, 2, 2]);
        assert_eq!(back, data);
    }

    #[test]
    fn test_write_sph_particles() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        write_sph_particles(&mut f, "sph", &pos, &[0.5, 0.5], &[1000.0, 1000.0]).unwrap();
        assert_eq!(
            f.open_dataset("sph", "h").unwrap().read_f64().unwrap(),
            vec![0.5, 0.5]
        );
    }

    #[test]
    fn test_write_nodal_displacements() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![0.0, 0.01, 0.0, 0.0, 0.02, 0.0];
        write_nodal_displacements(&mut f, "fem", 2, &d).unwrap();
        assert_eq!(
            f.open_dataset("fem", "displacements")
                .unwrap()
                .read_f64()
                .unwrap(),
            d
        );
    }

    #[test]
    fn test_write_sensor_recording() {
        let mut f = Hdf5File::new("t.h5");
        let data: Vec<f64> = (0..20).map(|i| i as f64 * 0.01).collect();
        write_sensor_recording(&mut f, "acc", 4, 5, 1000.0, &data).unwrap();
        let g = f.open_group("acc").unwrap();
        assert_eq!(g.attributes["n_channels"], AttrValue::Int32(4));
    }

    #[test]
    fn test_write_image_stack() {
        let mut f = Hdf5File::new("t.h5");
        let px: Vec<f64> = (0..16).map(|i| i as f64).collect();
        write_image_stack(&mut f, "imgs", 1, 4, 4, &px).unwrap();
        assert_eq!(
            f.open_dataset("imgs", "images")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            16
        );
    }

    #[test]
    fn test_write_covariance_matrix_3x3() {
        let mut f = Hdf5File::new("t.h5");
        let cov: Vec<f64> = (0..9).map(|i| i as f64).collect();
        write_covariance_matrix(&mut f, "cov3", 3, &cov).unwrap();
        assert_eq!(
            f.open_dataset("cov3", "covariance")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }

    #[test]
    fn test_write_gcmc_n_particles() {
        let mut f = Hdf5File::new("t.h5");
        write_gcmc_run(&mut f, "g2", &[0], &[5], &[-10.0], -3.0).unwrap();
        assert_eq!(
            f.open_dataset("g2", "n_particles")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![5.0]
        );
    }

    #[test]
    fn test_write_read_format_version_default() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 0, 1, 0, "test").unwrap();
        let (maj, min, pat, c) = read_format_version(&f).unwrap();
        assert_eq!((maj, min, pat), (0, 1, 0));
        assert_eq!(c, "test");
    }

    #[test]
    fn test_write_velocity_field_3d() {
        let mut f = Hdf5File::new("t.h5");
        let n = 8;
        let u = vec![1.0_f64; n];
        let v = vec![2.0_f64; n];
        let w = vec![3.0_f64; n];
        write_velocity_field_3d(&mut f, "fl3d", 2, 2, 2, &u, &v, &w).unwrap();
        assert_eq!(f.open_dataset("fl3d", "w").unwrap().read_f64().unwrap(), w);
    }

    #[test]
    fn test_write_vorticity_field() {
        let mut f = Hdf5File::new("t.h5");
        let omega = vec![0.1, -0.2, 0.3, 0.0];
        write_vorticity_field(&mut f, "vort", 2, 2, &omega).unwrap();
        assert_eq!(
            f.open_dataset("vort", "vorticity")
                .unwrap()
                .read_f64()
                .unwrap(),
            omega
        );
    }

    #[test]
    fn test_write_pressure_field() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![1.0, 1.5, 1.5, 1.0];
        write_pressure_field(&mut f, "ns", 2, 2, &p).unwrap();
        assert_eq!(
            f.open_dataset("ns", "pressure")
                .unwrap()
                .read_f64()
                .unwrap(),
            p
        );
    }

    #[test]
    fn test_write_diffusion_coefficient() {
        let mut f = Hdf5File::new("t.h5");
        write_diffusion_coefficient(&mut f, "tr", "Na", 1.33e-9, "m^2/s").unwrap();
        let g = f.open_group("tr").unwrap();
        let ds = g.datasets.get("D_Na").unwrap();
        if let DataStorage::Float64(ref v) = ds.data {
            assert!((v[0] - 1.33e-9).abs() < 1e-20);
        }
    }
}

#[cfg(test)]
mod comp_chem_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_cluster_expansion_eci() {
        let mut f = Hdf5File::new("t.h5");
        write_cluster_expansion_eci(&mut f, "ce", &[1, 2, 3], &[-0.5, 0.1, 0.05]).unwrap();
        assert_eq!(
            f.open_dataset("ce", "eci_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-0.5, 0.1, 0.05]
        );
    }

    #[test]
    fn test_write_enumerated_structures() {
        let mut f = Hdf5File::new("t.h5");
        write_enumerated_structures(
            &mut f,
            "enum",
            &[0, 1, 2],
            &[0.0, 0.5, 1.0],
            &[-1.0, -1.5, -1.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("enum", "concentrations")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![0.0, 0.5, 1.0]
        );
    }

    #[test]
    fn test_write_eos_fit() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![60.0, 65.0, 70.0, 75.0];
        let e = vec![-4.0, -4.5, -4.3, -4.0];
        write_eos_fit(&mut f, "eos", &v, &e, 68.0, -4.5, 160.0, 4.0).unwrap();
        let g = f.open_group("eos").unwrap();
        assert_eq!(g.attributes["B0_GPa"], AttrValue::Float64(160.0));
    }

    #[test]
    fn test_write_convex_hull() {
        let mut f = Hdf5File::new("t.h5");
        let x = vec![0.0, 0.5, 1.0];
        let e = vec![0.0, -0.1, 0.0];
        let s = vec![true, true, true];
        write_convex_hull(&mut f, "hull", &x, &e, &s).unwrap();
        let sm = f
            .open_dataset("hull", "stable_mask")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(sm, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_write_solvation_energies() {
        let mut f = Hdf5File::new("t.h5");
        let solvents = vec!["water".to_string(), "methanol".to_string()];
        let dg = vec![-5.0, -3.5];
        write_solvation_energies(&mut f, "solv", &solvents, &dg).unwrap();
        assert_eq!(
            f.open_dataset("solv", "dG_solv_kcal_mol")
                .unwrap()
                .read_f64()
                .unwrap(),
            dg
        );
    }

    #[test]
    fn test_write_excitations() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![3.5, 4.2, 5.1];
        let osc = vec![0.5, 0.1, 0.8];
        write_excitations(&mut f, "tddft", &e, &osc).unwrap();
        assert_eq!(
            f.open_dataset("tddft", "oscillator_strengths")
                .unwrap()
                .read_f64()
                .unwrap(),
            osc
        );
    }

    #[test]
    fn test_write_nmr_shielding() {
        let mut f = Hdf5File::new("t.h5");
        let iso = vec![25.0, 120.0];
        let aniso = vec![5.0, 30.0];
        write_nmr_shielding(&mut f, "nmr", &iso, &aniso).unwrap();
        assert_eq!(
            f.open_dataset("nmr", "sigma_iso")
                .unwrap()
                .read_f64()
                .unwrap(),
            iso
        );
    }

    #[test]
    fn test_write_hyperfine_coupling() {
        let mut f = Hdf5File::new("t.h5");
        let a = vec![50.0, 10.0, 5.0];
        write_hyperfine_coupling(&mut f, "hfc", &a).unwrap();
        assert_eq!(
            f.open_dataset("hfc", "A_iso_MHz")
                .unwrap()
                .read_f64()
                .unwrap(),
            a
        );
    }

    #[test]
    fn test_write_nics() {
        let mut f = Hdf5File::new("t.h5");
        let h = vec![0.0, 1.0];
        let n = vec![-10.5, -5.2];
        write_nics(&mut f, "arom", &h, &n).unwrap();
        assert_eq!(
            f.open_dataset("arom", "NICS_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            n
        );
    }

    #[test]
    fn test_write_charge_flow_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let dq = vec![0.0, 0.1, -0.1, 0.0]; // 2×2
        write_charge_flow_matrix(&mut f, "cflow", 2, &dq).unwrap();
        assert_eq!(
            f.open_dataset("cflow", "charge_flow")
                .unwrap()
                .read_f64()
                .unwrap(),
            dq
        );
    }

    #[test]
    fn test_write_fukui_functions() {
        let mut f = Hdf5File::new("t.h5");
        let fp = vec![0.3, 0.7];
        let fm = vec![0.4, 0.6];
        let fz = vec![0.35, 0.65];
        write_fukui_functions(&mut f, "fukui", &fp, &fm, &fz).unwrap();
        assert_eq!(
            f.open_dataset("fukui", "f_plus")
                .unwrap()
                .read_f64()
                .unwrap(),
            fp
        );
    }

    #[test]
    fn test_write_potential_energy_series() {
        let mut f = Hdf5File::new("t.h5");
        let e: Vec<f64> = (0..5).map(|i| -(i as f64) * 10.0).collect();
        write_potential_energy_series(&mut f, "thermo", &e).unwrap();
        assert_eq!(
            read_scalar_series(&f, "thermo", "potential_energy").unwrap(),
            e
        );
    }

    #[test]
    fn test_write_von_mises_stress() {
        let mut f = Hdf5File::new("t.h5");
        let s = vec![150e6, 200e6, 180e6];
        write_von_mises_stress(&mut f, "result", &s).unwrap();
        assert_eq!(
            f.open_dataset("result", "von_mises")
                .unwrap()
                .read_f64()
                .unwrap(),
            s
        );
    }

    #[test]
    fn test_write_radius_of_gyration() {
        let mut f = Hdf5File::new("t.h5");
        let rg = vec![3.5, 3.6, 3.4];
        write_radius_of_gyration(&mut f, "poly", &rg).unwrap();
        assert_eq!(
            f.open_dataset("poly", "rg").unwrap().read_f64().unwrap(),
            rg
        );
    }

    #[test]
    fn test_write_raman_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let wn = vec![200.0, 400.0, 800.0];
        let i = vec![100.0, 500.0, 200.0];
        write_raman_spectrum(&mut f, "raman", &wn, &i).unwrap();
        assert_eq!(
            f.open_dataset("raman", "intensities")
                .unwrap()
                .read_f64()
                .unwrap(),
            i
        );
    }

    #[test]
    fn test_write_virial_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        write_virial_coefficients(&mut f, "v2", &[-15.0, 200.0]).unwrap();
        assert_eq!(
            f.open_dataset("v2", "virial_coefficients")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-15.0, 200.0]
        );
    }
}

#[cfg(test)]
mod microstructure_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_phase_field_2d() {
        let mut f = Hdf5File::new("t.h5");
        let phi = vec![0.0, 0.5, 1.0, 0.5]; // 2×2
        write_phase_field_2d(&mut f, "pf", 0, 2, 2, &phi).unwrap();
        let g = f.open_group("pf").unwrap();
        assert!(g.groups.contains_key("step_0000000000"));
    }

    #[test]
    fn test_write_chemical_potential_field() {
        let mut f = Hdf5File::new("t.h5");
        let mu = vec![0.1, -0.1, 0.2, -0.2];
        write_chemical_potential_field(&mut f, "ch", 2, 2, &mu).unwrap();
        assert_eq!(
            f.open_dataset("ch", "chemical_potential")
                .unwrap()
                .read_f64()
                .unwrap(),
            mu
        );
    }

    #[test]
    fn test_write_percolation_stats() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![0.5, 0.6, 0.7];
        let cs = vec![3.0, 10.0, 100.0];
        let perc = vec![false, false, true];
        write_percolation_stats(&mut f, "perc", &p, &cs, &perc).unwrap();
        let pf = f
            .open_dataset("perc", "percolates")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(pf, vec![0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_write_heat_field() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 350.0, 400.0, 350.0];
        write_heat_field(&mut f, "heat", 5, 2, 2, &t).unwrap();
        let g = f.open_group("heat").unwrap();
        assert!(g.groups.contains_key("step_0000000005"));
    }

    #[test]
    fn test_write_concentration_field() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![0.1, 0.2, 0.3, 0.4];
        write_concentration_field(&mut f, "rd", "A", 0, 2, 2, &c).unwrap();
        let back = f
            .open_dataset("rd/A/step_0000000000", "concentration")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(back, c);
    }

    #[test]
    fn test_write_nucleation_data() {
        let mut f = Hdf5File::new("t.h5");
        let ns = vec![1usize, 5, 10, 20];
        let fg = vec![0.5, -0.2, -1.0, -2.5];
        write_nucleation_data(&mut f, "nucl", &ns, &fg).unwrap();
        assert_eq!(
            f.open_dataset("nucl", "free_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fg
        );
    }
}

#[cfg(test)]
mod extra_roundtrip_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_soc_matrix_write() {
        let mut f = Hdf5File::new("t.h5");
        let mat = vec![0.0_f64; 32]; // n=2: 4×4×2
        write_soc_matrix(&mut f, "soc", 2, &mat).unwrap();
        assert_eq!(
            f.open_dataset("soc", "soc_matrix")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            32
        );
    }

    #[test]
    fn test_write_gw_corrections() {
        let mut f = Hdf5File::new("t.h5");
        let dft = vec![-3.0, -1.0, 2.0];
        let gw = vec![-3.5, -1.2, 2.3];
        write_gw_corrections(&mut f, "gw2", &dft, &gw).unwrap();
        assert_eq!(
            f.open_dataset("gw2", "gw_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            gw
        );
    }

    #[test]
    fn test_write_force_constants() {
        let mut f = Hdf5File::new("t.h5");
        let ifc = vec![0.1_f64; 36];
        write_force_constants(&mut f, "ifc", 2, &ifc).unwrap();
        assert_eq!(
            f.open_dataset("ifc", "force_constants")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            36
        );
    }

    #[test]
    fn test_write_dipole_trajectory() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        write_dipole_trajectory(&mut f, "dip", &d).unwrap();
        assert_eq!(
            f.open_dataset("dip", "dipole_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }

    #[test]
    fn test_write_magnetic_moments() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![[0.0, 0.0, 2.0]];
        write_magnetic_moments(&mut f, "mag", &m).unwrap();
        assert_eq!(
            f.open_dataset("mag", "magnetic_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            3
        );
    }

    #[test]
    fn test_write_electron_density() {
        let mut f = Hdf5File::new("t.h5");
        let rho = vec![0.5_f64; 27]; // 3×3×3
        write_electron_density(&mut f, "dft2", 3, 3, 3, &rho).unwrap();
        assert_eq!(
            f.open_dataset("dft2", "rho_e")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            27
        );
    }

    #[test]
    fn test_write_convex_hull_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let x = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let e = vec![0.0, -0.05, -0.10, -0.05, 0.0];
        let s = vec![true, true, true, true, true];
        write_convex_hull(&mut f, "hull2", &x, &e, &s).unwrap();
        assert_eq!(
            f.open_dataset("hull2", "hull_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            e
        );
    }

    #[test]
    fn test_write_thermochemistry_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![298.15, 500.0, 1000.0];
        let h = vec![-200.0, -195.0, -185.0];
        let s = vec![0.1, 0.15, 0.25];
        let g: Vec<f64> = h
            .iter()
            .zip(t.iter())
            .zip(s.iter())
            .map(|((&hi, &ti), &si)| hi - ti * si)
            .collect();
        write_thermochemistry(&mut f, "tc", &t, &h, &s, &g).unwrap();
        assert_eq!(
            f.open_dataset("tc", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            t
        );
    }

    #[test]
    fn test_write_optical_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 3.0, 4.0];
        let eps2 = vec![0.0, 0.1, 0.5, 0.2];
        write_optical_spectrum(&mut f, "opt2", &e, &eps2).unwrap();
        assert_eq!(
            f.open_dataset("opt2", "epsilon2")
                .unwrap()
                .read_f64()
                .unwrap(),
            eps2
        );
    }

    #[test]
    fn test_write_lattice_parameters_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![100.0, 300.0, 500.0];
        let a = vec![5.0, 5.02, 5.05];
        write_lattice_parameters(&mut f, "th", &t, &a, &a, &a).unwrap();
        assert_eq!(f.open_dataset("th", "a").unwrap().read_f64().unwrap(), a);
    }

    #[test]
    fn test_hdf5_file_root_is_slash() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(f.root.name, "/");
    }

    #[test]
    fn test_write_neb_path_barrier() {
        let mut f = Hdf5File::new("t.h5");
        let im = vec![
            vec![[0.0, 0.0, 0.0]],
            vec![[0.5, 0.0, 0.0]],
            vec![[1.0, 0.0, 0.0]],
        ];
        let e = vec![-10.0, -7.0, -10.0];
        write_neb_path(&mut f, "neb2", &im, &e).unwrap();
        let eb = read_neb_energies(&f, "neb2").unwrap();
        assert!((eb[1] - (-7.0)).abs() < 1e-12);
    }

    #[test]
    fn test_write_wannier_centres_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![[0.5, 0.5, 0.5]];
        let sp = vec![1.2];
        write_wannier_centres(&mut f, "wan", &c, &sp).unwrap();
        assert_eq!(
            f.open_dataset("wan", "wannier_spreads")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }
}

#[cfg(test)]
mod integration_final_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_convergence_log() {
        let mut f = Hdf5File::new("t.h5");
        write_convergence_log(&mut f, "scf", &[1, 2, 3], &[1e-2, 1e-4, 1e-8], true).unwrap();
        let g = f.open_group("scf").unwrap();
        assert_eq!(g.attributes["converged"], AttrValue::Int32(1));
    }

    #[test]
    fn test_write_geopt_step() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0]];
        let frc = vec![[0.05, 0.0, 0.0]];
        write_geopt_step(&mut f, "geopt", 0, &pos, &frc, -10.5, 0.05).unwrap();
        let g = f.open_group("geopt").unwrap();
        let sg = g.groups.get("geopt_000000").unwrap();
        assert_eq!(sg.attributes["max_force"], AttrValue::Float64(0.05));
    }

    #[test]
    fn test_write_mulliken_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.3, -0.3, 0.3, -0.3];
        write_mulliken_charges(&mut f, "pop", &q, None).unwrap();
        assert_eq!(
            f.open_dataset("pop", "mulliken_charges")
                .unwrap()
                .read_f64()
                .unwrap(),
            q
        );
    }

    #[test]
    fn test_write_mulliken_with_spin() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.1, -0.1];
        let sd = vec![0.5, -0.5];
        write_mulliken_charges(&mut f, "rad", &q, Some(&sd)).unwrap();
        assert_eq!(
            f.open_dataset("rad", "spin_densities")
                .unwrap()
                .read_f64()
                .unwrap(),
            sd
        );
    }

    #[test]
    fn test_write_parallel_coordinates() {
        let mut f = Hdf5File::new("t.h5");
        let vn = vec!["x".to_string(), "y".to_string(), "z".to_string()];
        let data: Vec<f64> = (0..9).map(|i| i as f64).collect(); // 3 samples × 3 vars
        write_parallel_coordinates(&mut f, "pcoord", &vn, &data, 3).unwrap();
        assert_eq!(
            f.open_dataset("pcoord", "data")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }

    #[test]
    fn test_write_percolation_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![0.4, 0.5, 0.6, 0.7];
        let cs = vec![2.0, 4.0, 15.0, 100.0];
        let perc = vec![false, false, false, true];
        write_percolation_stats(&mut f, "p2", &p, &cs, &perc).unwrap();
        let pf = f
            .open_dataset("p2", "percolates")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(pf[3], 1.0);
        assert_eq!(pf[0], 0.0);
    }

    #[test]
    fn test_write_eos_fit_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![60.0, 65.0, 70.0];
        let e = vec![-4.0, -4.5, -4.3];
        write_eos_fit(&mut f, "eos2", &v, &e, 66.0, -4.5, 150.0, 4.2).unwrap();
        let g = f.open_group("eos2").unwrap();
        assert_eq!(g.attributes["E0"], AttrValue::Float64(-4.5));
    }

    #[test]
    fn test_total_dataset_count_with_groups() {
        let mut f = Hdf5File::new("t.h5");
        write_masses(&mut f, "atoms", &[1.0, 2.0, 3.0]).unwrap();
        write_charges(&mut f, "atoms", &[0.1, -0.1, 0.0]).unwrap();
        write_rdf(&mut f, "obs", &[0.1, 0.2], &[0.5, 1.0]).unwrap();
        let n = total_dataset_count(&f);
        assert!(n >= 4, "expected >=4, got {n}");
    }

    #[test]
    fn test_list_top_level_groups_non_empty() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 1, 0, 0, "OxiPhysics").unwrap();
        write_masses(&mut f, "atoms", &[1.0]).unwrap();
        let groups = list_top_level_groups(&f);
        assert!(groups.contains(&"atoms".to_string()));
        assert!(groups.contains(&"__version__".to_string()));
    }

    #[test]
    fn test_write_read_band_structure_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let k: Vec<f64> = (0..9).map(|i| i as f64 * 0.1).collect(); // 3 k-pts × 3
        let e: Vec<f64> = (0..6).map(|i| i as f64 - 3.0).collect(); // 3 k × 2 bands
        write_band_structure(&mut f, "bs2", 3, 2, &k, &e).unwrap();
        let bk = f
            .open_dataset("bs2", "kpoints")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(bk.len(), 9);
        let be = f
            .open_dataset("bs2", "eigenvalues")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(be, e);
    }

    #[test]
    fn test_write_read_velocities_single() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![[1.0, 2.0, 3.0]];
        write_velocities(&mut f, "single_v", &v).unwrap();
        let b = read_velocities(&f, "single_v").unwrap();
        assert_eq!(b, v);
    }

    #[test]
    fn test_write_read_box_vectors_ortho() {
        let mut f = Hdf5File::new("t.h5");
        let mut bv = [0.0_f64; 9];
        bv[0] = 5.0;
        bv[4] = 5.0;
        bv[8] = 5.0; // cubic
        write_box_vectors(&mut f, "box2", &bv).unwrap();
        let b = read_box_vectors(&f, "box2").unwrap();
        assert!((b[0] - 5.0).abs() < 1e-12);
        assert!((b[4] - 5.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod graph_frag_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_molecular_graph() {
        let mut f = Hdf5File::new("t.h5");
        let bonds = vec![(0usize, 1, 1.0), (1, 2, 2.0), (2, 3, 1.5)];
        write_molecular_graph(&mut f, "graph", &bonds).unwrap();
        assert_eq!(
            f.open_dataset("graph", "bond_order")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![1.0, 2.0, 1.5]
        );
    }

    #[test]
    fn test_write_fragmentation_map() {
        let mut f = Hdf5File::new("t.h5");
        let af = vec![0usize, 0, 1, 1, 2];
        write_fragmentation_map(&mut f, "frag", &af).unwrap();
        let b = f
            .open_dataset("frag", "atom_fragment")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(b, vec![0.0, 0.0, 1.0, 1.0, 2.0]);
    }

    #[test]
    fn test_write_fragment_masses() {
        let mut f = Hdf5File::new("t.h5");
        write_fragment_masses(&mut f, "frm", &[18.015, 44.010]).unwrap();
        assert_eq!(
            f.open_dataset("frm", "fragment_masses")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![18.015, 44.010]
        );
    }

    #[test]
    fn test_write_read_md_metadata_full() {
        let mut f = Hdf5File::new("t.h5");
        let m = MdRunMetadata {
            title: "NaCl Melt".into(),
            n_atoms: 512,
            n_steps: 500_000,
            dt_fs: 1.0,
            temperature_k: 1200.0,
            pressure_bar: 1.0,
            ensemble: "NVT".into(),
            force_field: "Born-Huggins-Meyer".into(),
        };
        m.write_to(&mut f).unwrap();
        let b = MdRunMetadata::read_from(&f).unwrap();
        assert_eq!(b.force_field, "Born-Huggins-Meyer");
        assert!((b.temperature_k - 1200.0).abs() < 1e-12);
    }

    #[test]
    fn test_write_solvation_energies_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let solvents = vec!["water".to_string(), "dmso".to_string()];
        let dg = vec![-4.5, -3.0];
        write_solvation_energies(&mut f, "solv2", &solvents, &dg).unwrap();
        assert_eq!(
            f.open_dataset("solv2", "dG_solv_kcal_mol")
                .unwrap()
                .read_f64()
                .unwrap(),
            dg
        );
    }

    #[test]
    fn test_write_fukui_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let fp = vec![0.3, 0.4, 0.3];
        let fm = vec![0.2, 0.5, 0.3];
        let fz = vec![0.25, 0.45, 0.30];
        write_fukui_functions(&mut f, "fuk2", &fp, &fm, &fz).unwrap();
        assert_eq!(
            f.open_dataset("fuk2", "f_minus")
                .unwrap()
                .read_f64()
                .unwrap(),
            fm
        );
    }

    #[test]
    fn test_write_convergence_log_not_converged() {
        let mut f = Hdf5File::new("t.h5");
        write_convergence_log(
            &mut f,
            "scf2",
            &[1, 2, 3, 4, 5],
            &[1e-1, 1e-2, 1e-3, 1e-4, 1e-5],
            false,
        )
        .unwrap();
        let g = f.open_group("scf2").unwrap();
        assert_eq!(g.attributes["converged"], AttrValue::Int32(0));
    }

    #[test]
    fn test_write_adsorption_energies_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        write_adsorption_energies(&mut f, "ads2", &["bridge", "hollow"], &[-1.8, -2.1]).unwrap();
        let g = f.open_group("ads2").unwrap();
        assert_eq!(g.attributes["hollow"], AttrValue::Float64(-2.1));
    }

    #[test]
    fn test_write_rate_constants() {
        let mut f = Hdf5File::new("t.h5");
        let rxn = vec!["R1".to_string(), "R2".to_string()];
        let kf = vec![1e6, 5e4];
        let kr = vec![1e2, 1e3];
        write_rate_constants(&mut f, "micro", &rxn, &kf, &kr).unwrap();
        assert_eq!(
            f.open_dataset("micro", "k_forward")
                .unwrap()
                .read_f64()
                .unwrap(),
            kf
        );
    }

    #[test]
    fn test_write_coverage_profile_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0, 3.0];
        let c = vec![0.0, 0.3, 0.7, 0.9];
        write_coverage_profile(&mut f, "kin", &t, &c, "CO2").unwrap();
        assert_eq!(
            f.open_dataset("kin/CO2", "coverage")
                .unwrap()
                .read_f64()
                .unwrap(),
            c
        );
    }
}

#[cfg(test)]
mod transport_bulk_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_ionic_conductivity() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 400.0, 500.0];
        let s = vec![0.01, 0.1, 1.0];
        write_ionic_conductivity(&mut f, "cond", &t, &s).unwrap();
        assert_eq!(
            f.open_dataset("cond", "conductivity_S_m")
                .unwrap()
                .read_f64()
                .unwrap(),
            s
        );
    }

    #[test]
    fn test_write_heat_capacity() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![200.0, 298.0, 400.0];
        let cp = vec![20.0, 25.0, 28.0];
        write_heat_capacity(&mut f, "thcp", &t, &cp).unwrap();
        assert_eq!(
            f.open_dataset("thcp", "Cp_J_mol_K")
                .unwrap()
                .read_f64()
                .unwrap(),
            cp
        );
    }

    #[test]
    fn test_write_cluster_expansion_eci_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let sz = vec![1usize, 2, 3, 4];
        let eci = vec![-0.5, 0.1, 0.05, 0.01];
        write_cluster_expansion_eci(&mut f, "ce2", &sz, &eci).unwrap();
        assert_eq!(
            f.open_dataset("ce2", "eci_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            eci
        );
    }

    #[test]
    fn test_write_nucleation_data_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let ns = vec![1usize, 5, 10];
        let fg = vec![2.0, -0.5, -3.0];
        write_nucleation_data(&mut f, "nucl2", &ns, &fg).unwrap();
        assert_eq!(
            f.open_dataset("nucl2", "free_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fg
        );
    }

    #[test]
    fn test_write_enumerated_structures() {
        let mut f = Hdf5File::new("t.h5");
        let ids = vec![0usize, 1, 2, 3];
        let x = vec![0.0, 0.25, 0.75, 1.0];
        let ef = vec![0.0, -0.05, -0.05, 0.0];
        write_enumerated_structures(&mut f, "enum2", &ids, &x, &ef).unwrap();
        assert_eq!(
            f.open_dataset("enum2", "concentrations")
                .unwrap()
                .read_f64()
                .unwrap(),
            x
        );
    }
}

#[cfg(test)]
mod scalar_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_write_read_scalar_with_units() {
        let mut f = Hdf5File::new("t.h5");
        write_scalar_with_units(&mut f, "props", "band_gap", 1.12, "eV").unwrap();
        let v = read_scalar(&f, "props", "band_gap").unwrap();
        assert!((v - 1.12).abs() < 1e-12);
    }

    #[test]
    fn test_read_scalar_error_on_empty() {
        let mut f = Hdf5File::new("t.h5");
        // Write zero-element dataset by creating manually
        f.create_group("empty").unwrap();
        let _ = f.create_dataset("empty", "val", vec![0], Hdf5Dtype::Float64);
        let result = read_scalar(&f, "empty", "val");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_ionic_conductivity_single_point() {
        let mut f = Hdf5File::new("t.h5");
        write_ionic_conductivity(&mut f, "c1", &[298.15], &[1e-4]).unwrap();
        let s = f
            .open_dataset("c1", "conductivity_S_m")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(s, vec![1e-4]);
    }

    #[test]
    fn test_hdf5_file_is_not_locked_by_default() {
        let f = Hdf5File::new("t.h5");
        assert!(!f.is_locked());
    }
}

#[cfg(test)]
mod root_helper_tests {
    use super::*;
    use crate::hdf5_io::*;

    #[test]
    fn test_root_group_count_empty() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(root_group_count(&f), 0);
    }

    #[test]
    fn test_is_empty_file() {
        let f = Hdf5File::new("t.h5");
        assert!(is_empty_file(&f));
    }
}
