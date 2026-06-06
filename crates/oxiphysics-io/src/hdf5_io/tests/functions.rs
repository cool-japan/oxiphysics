//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::hdf5_io::*;
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
        assert_eq!(compound.element_size(), 16);
    }
    #[test]
    fn test_dtype_named_inherits_base_size() {
        let named = Hdf5Dtype::Named {
            name: "energy_t".to_string(),
            base: Box::new(Hdf5Dtype::Float64),
        };
        assert_eq!(named.element_size(), 8);
    }
    #[test]
    fn test_chunk_layout_volume() {
        let cl = ChunkLayout::new(vec![32, 32, 3], 6);
        assert_eq!(cl.chunk_volume(), 32 * 32 * 3);
        assert_eq!(cl.gzip_level, 6);
    }
    #[test]
    fn test_chunk_layout_gzip_clamped() {
        let cl = ChunkLayout::new(vec![10], 20);
        assert_eq!(cl.gzip_level, 9);
    }
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
    #[test]
    fn test_create_group_nested() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("a/b/c").unwrap();
        let g = file.open_group("a/b/c").unwrap();
        assert_eq!(g.name, "c");
    }
    #[test]
    fn test_create_group_already_exists_is_ok() {
        let mut file = Hdf5File::new("test.h5");
        file.create_group("grp").unwrap();
        file.create_group("grp").unwrap();
    }
    #[test]
    fn test_open_nonexistent_group_fails() {
        let file = Hdf5File::new("test.h5");
        let r = file.open_group("nope");
        assert!(r.is_err());
    }
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
        file.create_dataset("g", "mat", vec![4, 4], Hdf5Dtype::Float64)
            .unwrap();
        let data: Vec<f64> = (0..16).map(|i| i as f64).collect();
        file.open_dataset_mut("g", "mat")
            .unwrap()
            .write_f64(&data)
            .unwrap();
        let slab = Hyperslab::new(vec![1, 1], vec![2, 2]);
        let ds = file.open_dataset("g", "mat").unwrap();
        let out = ds.read_hyperslab_f64(&slab).unwrap();
        assert_eq!(out, vec![5.0, 6.0, 9.0, 10.0]);
    }
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
        assert_eq!(mgr.count(), 2);
    }
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
    #[test]
    fn test_vlen_strings_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let strings = vec!["hello".to_string(), "world".to_string(), "HDF5".to_string()];
        write_vlen_strings(&mut file, "meta", "labels", &strings).unwrap();
        let back = read_vlen_strings(&file, "meta", "labels").unwrap();
        assert_eq!(back, strings);
    }
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
    #[test]
    fn test_chunk_iterator_1d() {
        let iter = ChunkIterator::new(vec![10], vec![3]);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[3].1, vec![1]);
    }
    #[test]
    fn test_chunk_iterator_2d() {
        let iter = ChunkIterator::new(vec![4, 4], vec![2, 2]);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 4);
    }
    #[test]
    fn test_collective_write_creates_dataset() {
        let mut file = Hdf5File::new("test.h5");
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let meta = collective_write_f64(&mut file, "par", "array", &data, 4).unwrap();
        assert_eq!(meta.n_ranks, 4);
        assert_eq!(meta.total_bytes, 800);
        let ds = file.open_dataset("par", "array").unwrap();
        assert_eq!(ds.volume(), 100);
    }
    #[test]
    fn test_assign_byte_offsets() {
        let mut file = Hdf5File::new("large.h5");
        write_f64_dataset(&mut file, "data", "arr", &vec![0.0; 100]).unwrap();
        assign_byte_offsets(&mut file);
        let ds = file.open_dataset("data", "arr").unwrap();
        assert!(ds.byte_offset > 0, "byte offset should be non-zero");
        assert!(file.superblock.eof_address > 0);
    }
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
    #[test]
    fn test_parallel_meta_total_bytes() {
        let mut meta = ParallelHdf5Meta::new(4);
        meta.record_rank_bytes(0, 1000);
        meta.record_rank_bytes(1, 2000);
        meta.record_rank_bytes(2, 3000);
        meta.record_rank_bytes(3, 4000);
        assert_eq!(meta.total_bytes(), 10000);
    }
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
    #[test]
    fn test_named_type_registry_list() {
        let mut reg = NamedDatatypeRegistry::default();
        reg.register("a", Hdf5Dtype::Float32).unwrap();
        reg.register("b", Hdf5Dtype::Int32).unwrap();
        let names = reg.names();
        assert!(names.contains(&"a".to_string()));
        assert!(names.contains(&"b".to_string()));
    }
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
    #[test]
    fn test_write_tensor3_f64() {
        let mut file = Hdf5File::new("test.h5");
        let data = vec![0.0_f64; 2 * 3 * 4];
        write_tensor3_f64(&mut file, "tdata", "tensor", 2, 3, 4, &data).unwrap();
        let ds = file.open_dataset("tdata", "tensor").unwrap();
        assert_eq!(ds.shape, vec![2, 3, 4]);
        assert_eq!(ds.volume(), 24);
    }
    #[test]
    fn test_write_tick_increments() {
        let mut file = Hdf5File::new("test.h5");
        let _t0 = file.write_tick;
        file.create_dataset("", "d", vec![1], Hdf5Dtype::Float64)
            .unwrap_or_default();
        let _ = file.write_tick;
    }
    #[test]
    fn test_superblock_default_values() {
        let sb = Hdf5Superblock::default();
        assert_eq!(sb.version, 2);
        assert_eq!(sb.size_of_offsets, 8);
        assert_eq!(sb.size_of_lengths, 8);
    }
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
    #[test]
    fn test_copy_group_datasets_basic() {
        let mut file = Hdf5File::new("test.h5");
        write_f64_dataset(&mut file, "src", "arr", &[1.0, 2.0, 3.0]).unwrap();
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
    use crate::hdf5_io::*;
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
    #[test]
    fn test_region_reference_dereference() {
        let mut file = Hdf5File::new("test.h5");
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        write_f64_dataset(&mut file, "g", "arr", &data).unwrap();
        let rref = RegionReference::new("g", "arr", Hyperslab::new(vec![5], vec![4]));
        let out = rref.dereference_f64(&file).unwrap();
        assert_eq!(out, vec![5.0, 6.0, 7.0, 8.0]);
    }
    #[test]
    fn test_3d_grid_write_read() {
        let mut file = Hdf5File::new("test.h5");
        let data = vec![1.0_f64; 2 * 3 * 4];
        write_3d_grid_f64(&mut file, "grid", "density", 2, 3, 4, &data).unwrap();
        let (nx, ny, nz, out) = read_3d_grid_f64(&file, "grid", "density").unwrap();
        assert_eq!((nx, ny, nz), (2, 3, 4));
        assert_eq!(out.len(), 24);
    }
    #[test]
    fn test_write_forces_shape() {
        let mut file = Hdf5File::new("test.h5");
        let forces = vec![0.0_f64; 10 * 5 * 3];
        write_forces(&mut file, "traj", &forces, 10, 5).unwrap();
        let ds = file.open_dataset("traj", "forces").unwrap();
        assert_eq!(ds.shape, vec![10, 5, 3]);
    }
    #[test]
    fn test_write_energies_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let energies = vec![-100.0_f64, -101.0, -99.5];
        write_energies(&mut file, "traj", &energies).unwrap();
        let ds = file.open_dataset("traj", "potential_energy").unwrap();
        let back = ds.read_f64().unwrap();
        assert_eq!(back, energies);
    }
    #[test]
    fn test_write_read_atom_types() {
        let mut file = Hdf5File::new("test.h5");
        let types = vec!["H".to_string(), "H".to_string(), "O".to_string()];
        write_atom_types(&mut file, "topology", &types).unwrap();
        let back = read_atom_types(&file, "topology").unwrap();
        assert_eq!(back, types);
    }
    #[test]
    fn test_write_distance_matrix_diagonal_zero() {
        let mut file = Hdf5File::new("test.h5");
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        write_distance_matrix(&mut file, "analysis", &positions).unwrap();
        let mat = read_matrix_f64(&file, "analysis", "distance_matrix").unwrap();
        for (i, row) in mat.iter().enumerate().take(3) {
            assert!(row[i].abs() < 1e-12, "diagonal[{i}]={}", row[i]);
        }
        assert!((mat[0][1] - 1.0).abs() < 1e-12, "d(0,1)={}", mat[0][1]);
    }
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
        let mut ring = RingTrajectory::new(2, 1);
        for i in 0..5_u32 {
            ring.append(&[i as f64, 0.0, 0.0]);
        }
        assert_eq!(ring.n_stored(), 2);
        let f0 = ring.read_frame(0).unwrap();
        assert!((f0[0][0] - 3.0).abs() < 1e-12, "oldest={}", f0[0][0]);
    }
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
    #[test]
    fn test_write_bfactors_roundtrip() {
        let mut file = Hdf5File::new("test.h5");
        let bf = vec![10.0_f64, 12.0, 8.5, 15.0];
        write_bfactors(&mut file, "atoms", &bf).unwrap();
        let ds = file.open_dataset("atoms", "bfactor").unwrap();
        let back = ds.read_f64().unwrap();
        assert_eq!(back, bf);
    }
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
    #[test]
    fn test_distance_matrix_symmetric() {
        let mut file = Hdf5File::new("test.h5");
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        write_distance_matrix(&mut file, "g", &positions).unwrap();
        let mat = read_matrix_f64(&file, "g", "distance_matrix").unwrap();
        assert!((mat[0][1] - 5.0).abs() < 1e-10);
        assert!((mat[1][0] - 5.0).abs() < 1e-10);
    }
}
