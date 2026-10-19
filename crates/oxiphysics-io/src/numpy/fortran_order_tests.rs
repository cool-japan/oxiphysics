// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reading `.npy` files that NumPy wrote in Fortran order, in other byte
//! orders or with newer header versions, and rejecting hostile headers.

use super::{
    NpyArray, NpyDtype, NpzArchive, NpzWriter, detect_npy_dtype, read_npy_f32, read_npy_f64,
    read_npy_i32, read_npy_i64, read_npy_shape, write_npy_f64,
};

/// `{'descr': ..., 'fortran_order': ..., 'shape': ..., }` as NumPy prints it.
fn header_dict(descr: &str, fortran_order: bool, shape: &str) -> String {
    let order = if fortran_order { "True" } else { "False" };
    format!("{{'descr': '{descr}', 'fortran_order': {order}, 'shape': {shape}, }}")
}

/// Assemble an `.npy` file of the given format version, padding the header with
/// spaces and a final newline so that the data start on a 64-byte boundary.
fn npy_file(version: u8, dict: &str, payload: &[u8]) -> Vec<u8> {
    let prefix = if version == 1 { 10 } else { 12 };
    let mut header = dict.as_bytes().to_vec();
    let unpadded = prefix + header.len() + 1;
    header.extend(std::iter::repeat_n(b' ', (64 - unpadded % 64) % 64));
    header.push(b'\n');
    let mut out = b"\x93NUMPY".to_vec();
    out.push(version);
    out.push(0);
    if version == 1 {
        out.extend_from_slice(&(header.len() as u16).to_le_bytes());
    } else {
        out.extend_from_slice(&(header.len() as u32).to_le_bytes());
    }
    out.extend_from_slice(&header);
    out.extend_from_slice(payload);
    out
}

fn f64_payload_le(values: &[f64]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

#[test]
fn read_npy_fortran_order_2x3_matches_c_order() {
    // The 2x3 matrix [[1, 2, 3], [4, 5, 6]] stored column by column.
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(2, 3)"),
        &f64_payload_le(&[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]),
    );
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert_eq!(shape, vec![2, 3]);
    assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn c_order_file_is_returned_unchanged() {
    let file = npy_file(
        1,
        &header_dict("<f8", false, "(2, 3)"),
        &f64_payload_le(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    );
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert_eq!(shape, vec![2, 3]);
    assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn transposed_array_saved_by_numpy_reads_back_as_the_transpose() {
    // `a = [[1, 2, 3], [4, 5, 6]]; np.save(f, a.T)`: a.T is F-contiguous, so
    // NumPy stores the bytes of `a` unchanged with fortran_order True and the
    // shape (3, 2).
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(3, 2)"),
        &f64_payload_le(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    );
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert_eq!(shape, vec![3, 2]);
    // a.T = [[1, 4], [2, 5], [3, 6]] in row-major order.
    assert_eq!(data, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
}

#[test]
fn fortran_order_three_axes_are_all_reordered() {
    let (n0, n1, n2) = (2usize, 3usize, 4usize);
    let value = |i: usize, j: usize, k: usize| (100 * i + 10 * j + k) as f64;
    // Fortran order: the first index varies fastest.
    let mut fortran = Vec::new();
    for k in 0..n2 {
        for j in 0..n1 {
            for i in 0..n0 {
                fortran.push(value(i, j, k));
            }
        }
    }
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(2, 3, 4)"),
        &f64_payload_le(&fortran),
    );
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert_eq!(shape, vec![2, 3, 4]);
    let mut expected = Vec::new();
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                expected.push(value(i, j, k));
            }
        }
    }
    assert_eq!(data, expected);
    assert_eq!(data[(3 + 2) * 4 + 1], 121.0);
}

#[test]
fn fortran_order_applies_to_every_element_type() {
    // [[1, 2], [3, 4]] stored column by column.
    let f32_payload: Vec<u8> = [1.0f32, 3.0, 2.0, 4.0]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let file = npy_file(1, &header_dict("<f4", true, "(2, 2)"), &f32_payload);
    assert_eq!(read_npy_f32(&file).unwrap().1, vec![1.0, 2.0, 3.0, 4.0]);

    let i32_payload: Vec<u8> = [1i32, 3, 2, 4]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let file = npy_file(1, &header_dict("<i4", true, "(2, 2)"), &i32_payload);
    assert_eq!(read_npy_i32(&file).unwrap().1, vec![1, 2, 3, 4]);

    let i64_payload: Vec<u8> = [1i64, 3, 2, 4]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let file = npy_file(1, &header_dict("<i8", true, "(2, 2)"), &i64_payload);
    assert_eq!(read_npy_i64(&file).unwrap().1, vec![1, 2, 3, 4]);
}

#[test]
fn fortran_order_is_the_identity_for_vectors_scalars_and_empty_arrays() {
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(4,)"),
        &f64_payload_le(&[4.0, 3.0, 2.0, 1.0]),
    );
    assert_eq!(read_npy_f64(&file).unwrap().1, vec![4.0, 3.0, 2.0, 1.0]);

    let file = npy_file(1, &header_dict("<f8", true, "()"), &f64_payload_le(&[7.5]));
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert!(shape.is_empty());
    assert_eq!(data, vec![7.5]);

    let file = npy_file(1, &header_dict("<f8", true, "(0, 3)"), &[]);
    let (shape, data) = read_npy_f64(&file).unwrap();
    assert_eq!(shape, vec![0, 3]);
    assert!(data.is_empty());
}

#[test]
fn shape_is_reported_unchanged_for_fortran_files() {
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(2, 3)"),
        &f64_payload_le(&[0.0; 6]),
    );
    assert_eq!(read_npy_shape(&file).unwrap(), vec![2, 3]);
    assert_eq!(detect_npy_dtype(&file).unwrap(), NpyDtype::Float64);
}

#[test]
fn missing_or_malformed_fortran_order_is_rejected() {
    let payload = f64_payload_le(&[1.0, 2.0]);
    let missing = npy_file(1, "{'descr': '<f8', 'shape': (2,), }", &payload);
    let err = read_npy_f64(&missing).unwrap_err();
    assert!(err.contains("fortran_order"), "{err}");

    let malformed = npy_file(
        1,
        "{'descr': '<f8', 'fortran_order': maybe, 'shape': (2,), }",
        &payload,
    );
    let err = read_npy_f64(&malformed).unwrap_err();
    assert!(err.contains("fortran_order"), "{err}");
}

#[test]
fn native_and_unprefixed_descriptors_are_accepted() {
    let payload = f64_payload_le(&[1.5, -2.5]);
    for descr in ["=f8", "f8"] {
        let file = npy_file(1, &header_dict(descr, false, "(2,)"), &payload);
        if cfg!(target_endian = "little") {
            assert_eq!(read_npy_f64(&file).unwrap().1, vec![1.5, -2.5], "{descr}");
        } else {
            assert!(read_npy_f64(&file).is_ok(), "{descr}");
        }
    }
}

#[test]
fn big_endian_arrays_are_decoded() {
    let payload: Vec<u8> = [1.5f64, -2.25, 1e300]
        .iter()
        .flat_map(|v| v.to_be_bytes())
        .collect();
    let file = npy_file(1, &header_dict(">f8", false, "(3,)"), &payload);
    assert_eq!(read_npy_f64(&file).unwrap().1, vec![1.5, -2.25, 1e300]);

    let payload: Vec<u8> = [-7i32, 258, i32::MAX]
        .iter()
        .flat_map(|v| v.to_be_bytes())
        .collect();
    let file = npy_file(1, &header_dict(">i4", true, "(3,)"), &payload);
    assert_eq!(read_npy_i32(&file).unwrap().1, vec![-7, 258, i32::MAX]);
}

#[test]
fn pipe_marker_is_only_valid_for_one_byte_types() {
    let payload = f64_payload_le(&[1.0]);
    let file = npy_file(1, &header_dict("|f8", false, "(1,)"), &payload);
    assert!(read_npy_f64(&file).is_err());
    assert!(NpyDtype::from_numpy_str("|f8").is_err());
    assert_eq!(NpyDtype::from_numpy_str("|u1"), Ok(NpyDtype::Uint8));
    assert_eq!(NpyDtype::from_numpy_str("|b1"), Ok(NpyDtype::Bool));
}

#[test]
fn dtype_strings_are_normalised() {
    assert_eq!(NpyDtype::from_numpy_str("<f8"), Ok(NpyDtype::Float64));
    assert_eq!(NpyDtype::from_numpy_str(">f8"), Ok(NpyDtype::Float64));
    assert_eq!(NpyDtype::from_numpy_str("=f4"), Ok(NpyDtype::Float32));
    assert_eq!(NpyDtype::from_numpy_str("<i4"), Ok(NpyDtype::Int32));
    assert_eq!(NpyDtype::from_numpy_str(">i8"), Ok(NpyDtype::Int64));
    assert_eq!(NpyDtype::from_numpy_str("?"), Ok(NpyDtype::Bool));
    for bad in [
        "",
        "<",
        "bad",
        "<f16",
        "<c8",
        "<f0",
        "<U5",
        "[('x', '<f8')]",
    ] {
        assert!(NpyDtype::from_numpy_str(bad).is_err(), "{bad:?}");
    }
}

#[test]
fn detect_dtype_handles_every_spelling_numpy_writes() {
    for (descr, expected) in [
        ("|u1", NpyDtype::Uint8),
        ("|b1", NpyDtype::Bool),
        ("<f4", NpyDtype::Float32),
        (">i8", NpyDtype::Int64),
    ] {
        let size = match descr {
            "|u1" | "|b1" => 1,
            "<f4" => 4,
            _ => 8,
        };
        let file = npy_file(1, &header_dict(descr, false, "(1,)"), &vec![0u8; size]);
        assert_eq!(detect_npy_dtype(&file).unwrap(), expected, "{descr}");
    }
}

#[test]
fn wrong_dtype_is_still_rejected() {
    let file = npy_file(
        1,
        &header_dict("<f8", false, "(1,)"),
        &f64_payload_le(&[1.0]),
    );
    assert!(read_npy_f32(&file).is_err());
    assert!(read_npy_i32(&file).is_err());
    assert!(read_npy_i64(&file).is_err());
    let err = read_npy_i64(&file).unwrap_err();
    assert!(err.contains("<i8") && err.contains("<f8"), "{err}");
}

#[test]
fn version_two_and_three_headers_are_read() {
    let payload = f64_payload_le(&[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
    for version in [2u8, 3] {
        let file = npy_file(version, &header_dict("<f8", true, "(2, 3)"), &payload);
        let (shape, data) = read_npy_f64(&file).unwrap();
        assert_eq!(shape, vec![2, 3], "version {version}");
        assert_eq!(
            data,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "version {version}"
        );
    }
    let mut unsupported = npy_file(1, &header_dict("<f8", false, "(1,)"), &payload);
    unsupported[6] = 4;
    assert!(read_npy_f64(&unsupported).is_err());
}

#[test]
fn hostile_shapes_are_errors_without_allocating() {
    for shape in [
        "(4294967296, 4294967296)",
        "(18446744073709551615, 2)",
        "(1000000000000,)",
        "(-1,)",
        "(2, x)",
    ] {
        let file = npy_file(1, &header_dict("<f8", true, shape), &f64_payload_le(&[1.0]));
        assert!(read_npy_f64(&file).is_err(), "{shape}");
        assert!(read_npy_i32(&file).is_err(), "{shape}");
    }
}

#[test]
fn damaged_headers_are_errors() {
    let good = write_npy_f64(&[2], &[1.0, 2.0]);
    for cut in [0usize, 5, 8, 9, 11, 40] {
        assert!(read_npy_f64(&good[..cut]).is_err(), "cut at {cut}");
    }
    let mut bad_len = good.clone();
    bad_len[8..10].copy_from_slice(&u16::MAX.to_le_bytes());
    assert!(read_npy_f64(&bad_len).is_err());
    let mut v2_huge = npy_file(2, &header_dict("<f8", false, "(1,)"), &[0u8; 8]);
    v2_huge[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(read_npy_f64(&v2_huge).is_err());
    let structured = npy_file(
        1,
        "{'descr': [('x', '<f8'), ('y', '<f8')], 'fortran_order': False, 'shape': (1,), }",
        &[0u8; 16],
    );
    let err = read_npy_f64(&structured).unwrap_err();
    assert!(err.contains("expected dtype"), "{err}");
}

#[test]
fn plain_round_trip_is_unchanged() {
    let data: Vec<f64> = (0..24).map(|i| f64::from(i) * 0.5).collect();
    let bytes = write_npy_f64(&[2, 3, 4], &data);
    let (shape, back) = read_npy_f64(&bytes).unwrap();
    assert_eq!(shape, vec![2, 3, 4]);
    assert_eq!(back, data);
}

#[test]
fn npz_archive_reads_fortran_ordered_members_in_c_order() {
    let mut writer = NpzWriter::new();
    let file = npy_file(
        1,
        &header_dict("<f8", true, "(2, 3)"),
        &f64_payload_le(&[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]),
    );
    writer.files.push(("m".to_string(), file));
    let archive = NpzArchive::from_bytes(&writer.to_bytes()).unwrap();
    let m: &NpyArray = archive.get("m").unwrap();
    assert_eq!(m.shape, vec![2, 3]);
    assert_eq!(m.data_f64, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
}

#[test]
fn npz_container_with_absurd_entry_count_is_an_error() {
    let mut bytes = u32::MAX.to_le_bytes().to_vec();
    bytes.extend_from_slice(&[0u8; 16]);
    assert!(NpzWriter::from_bytes(&bytes).is_err());
}
