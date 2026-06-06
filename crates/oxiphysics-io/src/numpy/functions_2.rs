//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::numpy::NpyArray;
    use crate::numpy::NpyDtype;
    use crate::numpy::NpyField;
    use crate::numpy::NpyMaskedArray;
    use crate::numpy::NpyRecordArray;
    use crate::numpy::NpySlice;
    use crate::numpy::NpzArchive;
    use crate::numpy::NpzWriter;
    use crate::numpy::*;
    /// write/read round-trip for f64: shape=\[3,2\], 6 values.
    #[test]
    fn test_roundtrip_f64() {
        let shape = vec![3usize, 2];
        let data: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let bytes = write_npy_f64(&shape, &data);
        let (got_shape, got_data) = read_npy_f64(&bytes).expect("read_npy_f64 failed");
        assert_eq!(got_shape, shape);
        assert_eq!(got_data.len(), data.len());
        for (a, b) in data.iter().zip(got_data.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    /// write/read round-trip for i32.
    #[test]
    fn test_roundtrip_i32() {
        let shape = vec![2usize, 3];
        let data: Vec<i32> = vec![-1, 0, 1, i32::MAX, i32::MIN, 42];
        let bytes = write_npy_i32(&shape, &data);
        let (got_shape, got_data) = read_npy_i32(&bytes).expect("read_npy_i32 failed");
        assert_eq!(got_shape, shape);
        assert_eq!(got_data, data);
    }
    /// Magic bytes check: `bytes[0..6] == b"\x93NUMPY"`.
    #[test]
    fn test_magic_bytes() {
        let shape = vec![2usize];
        let data = vec![1.0f64, 2.0f64];
        let bytes = write_npy_f64(&shape, &data);
        assert_eq!(&bytes[0..6], b"\x93NUMPY");
    }
    /// NpzWriter: add 2 arrays, serialize and deserialize, both recoverable.
    #[test]
    fn test_npz_roundtrip() {
        let mut writer = NpzWriter::new();
        let shape_f = vec![3usize, 2];
        let data_f: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let shape_i = vec![4usize];
        let data_i: Vec<i32> = vec![10, 20, 30, 40];
        writer.add_array_f64("matrix", &shape_f, &data_f);
        writer.add_array_i32("counts", &shape_i, &data_i);
        let bytes = writer.to_bytes();
        let recovered = NpzWriter::from_bytes(&bytes).expect("from_bytes failed");
        assert_eq!(recovered.files.len(), 2);
        let (s1, d1) = recovered
            .get_f64("matrix")
            .expect("matrix not found")
            .expect("read_npy_f64 failed");
        assert_eq!(s1, shape_f);
        assert_eq!(d1, data_f);
        let (s2, d2) = recovered
            .get_i32("counts")
            .expect("counts not found")
            .expect("read_npy_i32 failed");
        assert_eq!(s2, shape_i);
        assert_eq!(d2, data_i);
    }
    /// Shape encoding is correctly preserved.
    #[test]
    fn test_shape_encoding() {
        let shape = vec![5usize, 4, 3];
        let data = vec![0.0f64; 60];
        let bytes = write_npy_f64(&shape, &data);
        let (got_shape, _) = read_npy_f64(&bytes).expect("read_npy_f64 failed");
        assert_eq!(got_shape, shape);
    }
    /// NpyDtype::numpy_str returns expected strings.
    #[test]
    fn test_numpy_str() {
        assert_eq!(NpyDtype::Float64.numpy_str(), "<f8");
        assert_eq!(NpyDtype::Float32.numpy_str(), "<f4");
        assert_eq!(NpyDtype::Int32.numpy_str(), "<i4");
        assert_eq!(NpyDtype::Int64.numpy_str(), "<i8");
        assert_eq!(NpyDtype::Bool.numpy_str(), "?");
        assert_eq!(NpyDtype::Uint8.numpy_str(), "|u1");
    }
    /// 1-D array round-trip.
    #[test]
    fn test_1d_roundtrip() {
        let shape = vec![5usize];
        let data: Vec<f64> = vec![10.0, 20.0, 30.0, 40.0, 50.0];
        let bytes = write_npy_f64(&shape, &data);
        let (got_shape, got_data) = read_npy_f64(&bytes).expect("1d read failed");
        assert_eq!(got_shape, shape);
        assert_eq!(got_data, data);
    }
    #[test]
    fn test_roundtrip_f32() {
        let shape = vec![2usize, 3];
        let data: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let bytes = write_npy_f32(&shape, &data);
        let (got_shape, got_data) = read_npy_f32(&bytes).expect("read_npy_f32 failed");
        assert_eq!(got_shape, shape);
        assert_eq!(got_data, data);
    }
    #[test]
    fn test_roundtrip_i64() {
        let shape = vec![3usize];
        let data: Vec<i64> = vec![i64::MIN, 0, i64::MAX];
        let bytes = write_npy_i64(&shape, &data);
        let (got_shape, got_data) = read_npy_i64(&bytes).expect("read_npy_i64 failed");
        assert_eq!(got_shape, shape);
        assert_eq!(got_data, data);
    }
    #[test]
    fn test_element_size() {
        assert_eq!(NpyDtype::Float64.element_size(), 8);
        assert_eq!(NpyDtype::Float32.element_size(), 4);
        assert_eq!(NpyDtype::Int32.element_size(), 4);
        assert_eq!(NpyDtype::Int64.element_size(), 8);
        assert_eq!(NpyDtype::Bool.element_size(), 1);
        assert_eq!(NpyDtype::Uint8.element_size(), 1);
    }
    #[test]
    fn test_dtype_from_str() {
        assert_eq!(NpyDtype::from_numpy_str("<f8"), Ok(NpyDtype::Float64));
        assert_eq!(NpyDtype::from_numpy_str("<f4"), Ok(NpyDtype::Float32));
        assert_eq!(NpyDtype::from_numpy_str("<i4"), Ok(NpyDtype::Int32));
        assert!(NpyDtype::from_numpy_str("bad").is_err());
    }
    #[test]
    fn test_npy_array_validate() {
        let arr = NpyArray::from_f64(vec![2, 3], vec![1.0; 6]);
        assert!(arr.validate().is_ok());
        let bad = NpyArray::from_f64(vec![2, 3], vec![1.0; 5]);
        assert!(bad.validate().is_err());
    }
    #[test]
    fn test_npy_array_reshape() {
        let mut arr = NpyArray::from_f64(vec![2, 3], vec![1.0; 6]);
        assert!(arr.reshape(vec![3, 2]).is_ok());
        assert_eq!(arr.shape, vec![3, 2]);
        assert!(arr.reshape(vec![4, 2]).is_err());
    }
    #[test]
    fn test_npy_array_from_f32() {
        let arr = NpyArray::from_f32(vec![4], vec![1.0f32; 4]);
        assert_eq!(arr.dtype, NpyDtype::Float32);
        assert_eq!(arr.ndim(), 1);
        assert_eq!(arr.numel(), 4);
    }
    #[test]
    fn test_npy_array_from_i32() {
        let arr = NpyArray::from_i32(vec![2, 2], vec![1, 2, 3, 4]);
        assert_eq!(arr.dtype, NpyDtype::Int32);
        assert_eq!(arr.ndim(), 2);
    }
    #[test]
    fn test_validate_shape_ok() {
        assert!(validate_shape(&[2, 3], 6).is_ok());
        assert!(validate_shape(&[5], 5).is_ok());
    }
    #[test]
    fn test_validate_shape_err() {
        assert!(validate_shape(&[2, 3], 7).is_err());
    }
    #[test]
    fn test_flat_index() {
        assert_eq!(flat_index(&[1, 2], &[3, 4]).unwrap(), 6);
        assert_eq!(flat_index(&[0, 0], &[3, 4]).unwrap(), 0);
    }
    #[test]
    fn test_flat_index_out_of_range() {
        assert!(flat_index(&[3, 0], &[3, 4]).is_err());
    }
    #[test]
    fn test_unravel_index() {
        let indices = unravel_index(6, &[3, 4]).unwrap();
        assert_eq!(indices, vec![1, 2]);
    }
    #[test]
    fn test_unravel_flat_roundtrip() {
        let shape = vec![3, 4, 5];
        for flat in 0..60 {
            let indices = unravel_index(flat, &shape).unwrap();
            let recovered = flat_index(&indices, &shape).unwrap();
            assert_eq!(flat, recovered, "round-trip failed for flat={flat}");
        }
    }
    #[test]
    fn test_detect_dtype() {
        let bytes = write_npy_f64(&[2], &[1.0, 2.0]);
        assert_eq!(detect_npy_dtype(&bytes).unwrap(), NpyDtype::Float64);
        let bytes = write_npy_i32(&[3], &[1, 2, 3]);
        assert_eq!(detect_npy_dtype(&bytes).unwrap(), NpyDtype::Int32);
    }
    #[test]
    fn test_read_npy_shape() {
        let bytes = write_npy_f64(&[3, 4, 5], &vec![0.0; 60]);
        let shape = read_npy_shape(&bytes).unwrap();
        assert_eq!(shape, vec![3, 4, 5]);
    }
    #[test]
    fn test_npz_names_contains_remove() {
        let mut w = NpzWriter::new();
        w.add_array_f64("a", &[2], &[1.0, 2.0]);
        w.add_array_i32("b", &[3], &[1, 2, 3]);
        assert_eq!(w.len(), 2);
        assert!(w.contains("a"));
        assert!(!w.contains("c"));
        let names = w.names();
        assert!(names.contains(&"a"));
        assert!(names.contains(&"b"));
        assert!(w.remove("a"));
        assert_eq!(w.len(), 1);
        assert!(!w.contains("a"));
        assert!(!w.remove("a"));
    }
    #[test]
    fn test_npz_f32_i64() {
        let mut w = NpzWriter::new();
        w.add_array_f32("floats", &[3], &[1.0f32, 2.0, 3.0]);
        w.add_array_i64("longs", &[2], &[100i64, 200]);
        let bytes = w.to_bytes();
        let recovered = NpzWriter::from_bytes(&bytes).unwrap();
        let (shape, data) = recovered.get_f32("floats").unwrap().unwrap();
        assert_eq!(shape, vec![3]);
        assert_eq!(data, vec![1.0f32, 2.0, 3.0]);
        let (shape, data) = recovered.get_i64("longs").unwrap().unwrap();
        assert_eq!(shape, vec![2]);
        assert_eq!(data, vec![100i64, 200]);
    }
    #[test]
    fn test_scalar_array() {
        let shape: Vec<usize> = vec![];
        let data = vec![42.0_f64];
        let bytes = write_npy_f64(&shape, &data);
        let (got_shape, got_data) = read_npy_f64(&bytes).unwrap();
        assert!(got_shape.is_empty());
        assert_eq!(got_data.len(), 1);
        assert_eq!(got_data[0], 42.0);
    }
    #[test]
    fn test_wrong_dtype_error() {
        let bytes = write_npy_f64(&[2], &[1.0, 2.0]);
        assert!(read_npy_i32(&bytes).is_err());
    }
    #[test]
    fn test_truncated_data_error() {
        let bytes = write_npy_f64(&[10], &[0.0; 10]);
        let truncated = &bytes[..bytes.len() - 10];
        assert!(read_npy_f64(truncated).is_err());
    }
    #[test]
    fn test_bad_magic_error() {
        let mut bytes = write_npy_f64(&[2], &[1.0, 2.0]);
        bytes[0] = 0;
        assert!(read_npy_f64(&bytes).is_err());
    }
    #[test]
    fn test_npy_slice_row() {
        let data: Vec<f64> = (0..12).map(|i| i as f64).collect();
        let slice = NpySlice::new(&data, vec![3, 4]).unwrap();
        let row1 = slice.row(1).unwrap();
        assert_eq!(row1, &[4.0, 5.0, 6.0, 7.0]);
    }
    #[test]
    fn test_npy_slice_get() {
        let data: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let slice = NpySlice::new(&data, vec![2, 3]).unwrap();
        assert_eq!(slice.get(&[1, 2]).unwrap(), 5.0);
    }
    #[test]
    fn test_npy_slice_shape_mismatch() {
        let data = vec![1.0; 6];
        assert!(NpySlice::new(&data, vec![2, 4]).is_err());
    }
    #[test]
    fn test_masked_array_mean_valid() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let mask = vec![false, true, false, false];
        let ma = NpyMaskedArray::new(data, mask, vec![4], 1e20).unwrap();
        let mean = ma.mean_valid().unwrap();
        assert!((mean - 8.0 / 3.0).abs() < 1e-10, "mean={mean}");
    }
    #[test]
    fn test_masked_array_filled() {
        let data = vec![1.0, 2.0];
        let mask = vec![false, true];
        let ma = NpyMaskedArray::new(data, mask, vec![2], 999.0).unwrap();
        let filled = ma.filled();
        assert_eq!(filled[0], 1.0);
        assert_eq!(filled[1], 999.0);
    }
    #[test]
    fn test_masked_array_count_valid() {
        let ma = NpyMaskedArray::from_data(vec![1.0, 2.0, 3.0], vec![3]).unwrap();
        assert_eq!(ma.count_valid(), 3);
    }
    #[test]
    fn test_masked_array_mask_greater_than() {
        let mut ma = NpyMaskedArray::from_data(vec![1.0, 5.0, 2.0, 10.0], vec![4]).unwrap();
        ma.mask_greater_than(4.0);
        assert!(!ma.mask[0]);
        assert!(ma.mask[1]);
        assert!(!ma.mask[2]);
        assert!(ma.mask[3]);
        assert_eq!(ma.count_valid(), 2);
    }
    #[test]
    fn test_slice_mean() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(slice_mean(&data).unwrap(), 3.0);
        assert!(slice_mean(&[]).is_none());
    }
    #[test]
    fn test_slice_var() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let v = slice_var(&data).unwrap();
        assert!((v - 4.0).abs() < 1e-10, "var={v}");
    }
    #[test]
    fn test_slice_std() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s = slice_std(&data).unwrap();
        assert!((s - 2.0).abs() < 1e-10, "std={s}");
    }
    #[test]
    fn test_slice_min_max() {
        let data = vec![3.0, 1.0, 4.0, 1.5, 9.0, 2.6];
        let (min_v, min_i, max_v, max_i) = slice_min_max(&data).unwrap();
        assert_eq!(min_v, 1.0);
        assert_eq!(min_i, 1);
        assert_eq!(max_v, 9.0);
        assert_eq!(max_i, 4);
    }
    #[test]
    fn test_slice_percentile_median() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let median = slice_percentile(&data, 50.0).unwrap();
        assert!((median - 3.0).abs() < 1e-10, "median={median}");
    }
    #[test]
    fn test_slice_clip() {
        let data = vec![-1.0, 0.5, 2.0, 3.5];
        let clipped = slice_clip(&data, 0.0, 2.0);
        assert_eq!(clipped, vec![0.0, 0.5, 2.0, 2.0]);
    }
    #[test]
    fn test_slice_dot() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let dot = slice_dot(&a, &b).unwrap();
        assert!((dot - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_npz_archive_roundtrip() {
        let mut archive = NpzArchive::new();
        archive.insert("pos", NpyArray::from_f64(vec![3], vec![1.0, 2.0, 3.0]));
        archive.insert("idx", NpyArray::from_i32(vec![2], vec![10, 20]));
        let bytes = archive.to_bytes().unwrap();
        let recovered = NpzArchive::from_bytes(&bytes).unwrap();
        assert_eq!(recovered.len(), 2);
        let pos = recovered.get("pos").unwrap();
        assert_eq!(pos.dtype, NpyDtype::Float64);
        assert_eq!(pos.data_f64, vec![1.0, 2.0, 3.0]);
        let idx = recovered.get("idx").unwrap();
        assert_eq!(idx.data_i32, vec![10, 20]);
    }
    #[test]
    fn test_npz_archive_names_remove() {
        let mut archive = NpzArchive::new();
        archive.insert("a", NpyArray::from_f64(vec![1], vec![1.0]));
        archive.insert("b", NpyArray::from_f64(vec![1], vec![2.0]));
        assert!(archive.names().contains(&"a"));
        assert!(archive.remove("a"));
        assert!(!archive.names().contains(&"a"));
        assert_eq!(archive.len(), 1);
    }
    #[test]
    fn test_record_array_push_and_get() {
        let fields = vec![
            NpyField::scalar("x", NpyDtype::Float64),
            NpyField::scalar("y", NpyDtype::Float64),
            NpyField::scalar("mass", NpyDtype::Float64),
        ];
        let mut ra = NpyRecordArray::new(fields);
        ra.push_record(&[1.0, 2.0, 12.0]).unwrap();
        ra.push_record(&[3.0, 4.0, 16.0]).unwrap();
        assert_eq!(ra.n_records, 2);
        assert!((ra.get_scalar(1, "mass").unwrap() - 16.0).abs() < 1e-10);
    }
    #[test]
    fn test_record_array_column() {
        let fields = vec![
            NpyField::scalar("vx", NpyDtype::Float64),
            NpyField::scalar("vy", NpyDtype::Float64),
        ];
        let mut ra = NpyRecordArray::new(fields);
        ra.push_record(&[0.1, 0.2]).unwrap();
        ra.push_record(&[0.3, 0.4]).unwrap();
        let col = ra.column("vx").unwrap();
        assert_eq!(col, &[0.1, 0.3]);
    }
    #[test]
    fn test_linspace() {
        let v = linspace(0.0, 1.0, 5);
        assert_eq!(v.len(), 5);
        assert!((v[0] - 0.0).abs() < 1e-10);
        assert!((v[2] - 0.5).abs() < 1e-10);
        assert!((v[4] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_linspace_single() {
        let v = linspace(3.0, 3.0, 1);
        assert_eq!(v, vec![3.0]);
    }
    #[test]
    fn test_arange() {
        let v = arange(0.0, 1.0, 0.25).unwrap();
        assert_eq!(v.len(), 4);
        assert!((v[0] - 0.0).abs() < 1e-10);
        assert!((v[3] - 0.75).abs() < 1e-10);
    }
    #[test]
    fn test_arange_zero_step_error() {
        assert!(arange(0.0, 1.0, 0.0).is_err());
    }
    #[test]
    fn test_logspace() {
        let v = logspace(0.0, 2.0, 3);
        assert_eq!(v.len(), 3);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!((v[1] - 10.0).abs() < 1e-10);
        assert!((v[2] - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_transpose_2d() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let (t, shape) = transpose_2d(&data, &[2, 3]).unwrap();
        assert_eq!(shape, vec![3, 2]);
        assert_eq!(t[0], 1.0);
        assert_eq!(t[1], 4.0);
    }
    #[test]
    fn test_transpose_2d_square() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let (t, shape) = transpose_2d(&data, &[2, 2]).unwrap();
        assert_eq!(shape, vec![2, 2]);
        assert_eq!(t, vec![1.0, 3.0, 2.0, 4.0]);
    }
    #[test]
    fn test_matmul_2x3_3x2() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let b = vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let (c, shape) = matmul(&a, &[2, 3], &b, &[3, 2]).unwrap();
        assert_eq!(shape, vec![2, 2]);
        assert!((c[0] - 58.0).abs() < 1e-10);
        assert!((c[1] - 64.0).abs() < 1e-10);
        assert!((c[2] - 139.0).abs() < 1e-10);
        assert!((c[3] - 154.0).abs() < 1e-10);
    }
    #[test]
    fn test_matmul_identity() {
        let id = vec![1.0, 0.0, 0.0, 1.0];
        let a = vec![3.0, 4.0, 5.0, 6.0];
        let (c, _) = matmul(&id, &[2, 2], &a, &[2, 2]).unwrap();
        assert_eq!(c, a);
    }
    #[test]
    fn test_save_structured_magic() {
        let mut data_bytes = Vec::new();
        for v in &[1.0f64, 2.0f64, 3.0f64, 4.0f64] {
            data_bytes.extend_from_slice(&v.to_le_bytes());
        }
        let bytes =
            NpyArray::save_structured(&[("a", "<f8"), ("b", "<f8")], 2, &data_bytes).unwrap();
        assert_eq!(&bytes[0..6], NPY_MAGIC.as_ref());
    }
    #[test]
    fn test_save_structured_header_contains_field_names() {
        let data_bytes = vec![0u8; 8];
        let bytes = NpyArray::save_structured(&[("pressure", "<f8")], 1, &data_bytes).unwrap();
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let header = std::str::from_utf8(&bytes[10..10 + header_len]).unwrap();
        assert!(
            header.contains("pressure"),
            "header should contain field name 'pressure'"
        );
    }
    #[test]
    fn test_save_structured_empty_fields_error() {
        let result = NpyArray::save_structured(&[], 0, &[]);
        assert!(result.is_err());
    }
    #[test]
    fn test_save_structured_header_len_multiple_64() {
        let data_bytes = vec![0u8; 16];
        let bytes = NpyArray::save_structured(&[("x", "<f8")], 2, &data_bytes).unwrap();
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        assert_eq!(header_len % 64, 0, "header_len should be multiple of 64");
    }
    #[test]
    fn test_npz_archive_add_array_replaces() {
        let mut archive = NpzArchive::new();
        archive.add_array("v", NpyArray::from_f64(vec![2], vec![1.0, 2.0]));
        archive.add_array("v", NpyArray::from_f64(vec![3], vec![9.0, 8.0, 7.0]));
        assert_eq!(archive.len(), 1, "add_array should replace existing entry");
        assert_eq!(archive.get("v").unwrap().shape, vec![3]);
    }
    #[test]
    fn test_npz_archive_load_all_roundtrip() {
        let mut archive = NpzArchive::new();
        archive.add_array(
            "coords",
            NpyArray::from_f64(vec![2, 3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        );
        archive.add_array("ids", NpyArray::from_i32(vec![2], vec![0, 1]));
        let bytes = archive.to_bytes().unwrap();
        let loaded = NpzArchive::load_all(&bytes).unwrap();
        assert_eq!(loaded.len(), 2);
        let coords = loaded.get("coords").unwrap();
        assert_eq!(coords.shape, vec![2, 3]);
        assert!((coords.data_f64[5] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_npz_archive_iter() {
        let mut archive = NpzArchive::new();
        archive.add_array("a", NpyArray::from_f64(vec![1], vec![1.0]));
        archive.add_array("b", NpyArray::from_f64(vec![1], vec![2.0]));
        let names: Vec<&str> = archive.iter().map(|(n, _)| n).collect();
        assert!(names.contains(&"a"));
        assert!(names.contains(&"b"));
    }
    #[test]
    fn test_npz_archive_merge() {
        let mut a = NpzArchive::new();
        a.add_array("x", NpyArray::from_f64(vec![1], vec![1.0]));
        let mut b = NpzArchive::new();
        b.add_array("y", NpyArray::from_f64(vec![1], vec![2.0]));
        b.add_array("x", NpyArray::from_f64(vec![1], vec![99.0]));
        a.merge(b);
        assert_eq!(a.len(), 2);
        assert!((a.get("x").unwrap().data_f64[0] - 99.0).abs() < 1e-12);
    }
    #[test]
    fn test_npz_archive_total_elements() {
        let mut archive = NpzArchive::new();
        archive.add_array("a", NpyArray::from_f64(vec![3], vec![1.0, 2.0, 3.0]));
        archive.add_array("b", NpyArray::from_i32(vec![2, 2], vec![1, 2, 3, 4]));
        assert_eq!(archive.total_elements(), 7);
    }
}
