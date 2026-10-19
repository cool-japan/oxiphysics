// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the NetCDF classic binary encoder/decoder, including golden
//! files produced by the reference netCDF-C library.

use super::*;
use crate::netcdf::types::{NetCdfFile, NetCdfVariable};
use crate::netcdf::{create_amber_trajectory, parse_ncf_text, write_ncf_text};

/// netCDF-C 4.x `NETCDF3_CLASSIC` file (written by netCDF4-python, fill off): record dim `time` (3), `x` = 3, `name_len` = 5; `title` = "odd", `version` = int [1, 2, 3]; `float temp(time, x)` (units "K"), `short counts(time)`, `char label(name_len)` = "hello", `byte flags(x)`, `double pi`.
const REF_CLASSIC: [u8; 412] = [
    0x43, 0x44, 0x46, 0x01, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x00, 0x00, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x04, 0x74, 0x69, 0x6d, 0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
    0x78, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x08, 0x6e, 0x61, 0x6d, 0x65,
    0x5f, 0x6c, 0x65, 0x6e, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x02,
    0x00, 0x00, 0x00, 0x05, 0x74, 0x69, 0x74, 0x6c, 0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02,
    0x00, 0x00, 0x00, 0x03, 0x6f, 0x64, 0x64, 0x00, 0x00, 0x00, 0x00, 0x07, 0x76, 0x65, 0x72, 0x73,
    0x69, 0x6f, 0x6e, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x00, 0x00, 0x05,
    0x00, 0x00, 0x00, 0x04, 0x74, 0x65, 0x6d, 0x70, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x05,
    0x75, 0x6e, 0x69, 0x74, 0x73, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x01,
    0x4b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x01, 0x6c,
    0x00, 0x00, 0x00, 0x06, 0x63, 0x6f, 0x75, 0x6e, 0x74, 0x73, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x01, 0x78, 0x00, 0x00, 0x00, 0x05, 0x6c, 0x61, 0x62, 0x65,
    0x6c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x01, 0x58,
    0x00, 0x00, 0x00, 0x05, 0x66, 0x6c, 0x61, 0x67, 0x73, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x01, 0x60, 0x00, 0x00, 0x00, 0x02, 0x70, 0x69, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06,
    0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x01, 0x64, 0x68, 0x65, 0x6c, 0x6c, 0x6f, 0x00, 0x00, 0x00,
    0xff, 0x00, 0x01, 0x00, 0x40, 0x09, 0x21, 0xfb, 0x54, 0x44, 0x2d, 0x18, 0x3f, 0xc0, 0x00, 0x00,
    0x40, 0x20, 0x00, 0x00, 0x40, 0x60, 0x00, 0x00, 0x00, 0x0a, 0x00, 0x00, 0x40, 0x90, 0x00, 0x00,
    0x40, 0xb0, 0x00, 0x00, 0x40, 0xd0, 0x00, 0x00, 0xff, 0xec, 0x00, 0x00, 0x40, 0xf0, 0x00, 0x00,
    0x41, 0x08, 0x00, 0x00, 0x41, 0x18, 0x00, 0x00, 0x00, 0x1e, 0x00, 0x00,
];

/// netCDF-C `NETCDF3_64BIT_OFFSET` file: lone record variable `short s(rec)` = [1, 2, 3], stored without per-record padding.
const REF_OFFSET64: [u8; 90] = [
    0x43, 0x44, 0x46, 0x02, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x03, 0x72, 0x65, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
    0x73, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x54, 0x00, 0x01, 0x00, 0x02, 0x00, 0x03,
];

/// netCDF-C `NETCDF3_64BIT_DATA` (CDF-5) file: `n` = 2, `int64 big(n)`, `ubyte ub(n)` with `ushort scale` = 65000.
const REF_DATA64: [u8; 240] = [
    0x43, 0x44, 0x46, 0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
    0x6e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x62, 0x69, 0x67, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xdc,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x75, 0x62, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0c,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05,
    0x73, 0x63, 0x61, 0x6c, 0x65, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0xfd, 0xe8, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xec, 0xff, 0xff, 0xff, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x07, 0xc8, 0xff, 0x00, 0x00,
];

/// The dataset stored in [`REF_CLASSIC`], built through the public API.
fn reference_classic_dataset() -> CdfDataset {
    let mut ds = CdfDataset::new();
    ds.add_unlimited_dimension("time", 3).unwrap();
    ds.add_dimension("x", 3).unwrap();
    ds.add_dimension("name_len", 5).unwrap();
    ds.set_text_attribute("title", "odd");
    ds.set_attribute(CdfAttribute::new("version", NcValues::Int(vec![1, 2, 3])));
    let t = ds
        .add_variable(
            "temp",
            &["time", "x"],
            NcValues::Float(vec![1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]),
        )
        .unwrap();
    ds.variables[t].set_attribute(CdfAttribute::text("units", "K"));
    ds.add_variable("counts", &["time"], NcValues::Short(vec![10, -20, 30]))
        .unwrap();
    ds.add_variable("label", &["name_len"], NcValues::text("hello"))
        .unwrap();
    ds.add_variable("flags", &["x"], NcValues::Byte(vec![-1, 0, 1]))
        .unwrap();
    ds.add_variable("pi", &[], NcValues::Double(vec![std::f64::consts::PI]))
        .unwrap();
    ds
}

fn be32(bytes: &[u8], at: usize) -> u32 {
    u32::from_be_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
}

#[test]
fn reads_reference_classic_file() {
    let ds = read_cdf_dataset(&REF_CLASSIC).unwrap();
    assert_eq!(ds.num_records(), 3);
    assert_eq!(ds.record_dimension(), Some(0));
    assert_eq!(ds.dimension("x").map(|d| d.len), Some(3));
    assert_eq!(ds.text_attribute("title").as_deref(), Some("odd"));
    assert_eq!(
        ds.attribute("version").map(|a| a.values.clone()),
        Some(NcValues::Int(vec![1, 2, 3]))
    );
    let temp = ds.variable("temp").unwrap();
    assert_eq!(temp.nc_type(), NcType::Float);
    assert_eq!(temp.text_attribute("units").as_deref(), Some("K"));
    assert_eq!(
        temp.values,
        NcValues::Float(vec![1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
    );
    assert!(ds.is_record_variable(temp));
    assert_eq!(
        ds.variable("counts").unwrap().values,
        NcValues::Short(vec![10, -20, 30])
    );
    assert_eq!(
        ds.variable("label").unwrap().values.as_text().as_deref(),
        Some("hello")
    );
    assert_eq!(
        ds.variable("flags").unwrap().values,
        NcValues::Byte(vec![-1, 0, 1])
    );
    assert_eq!(
        ds.variable("pi").unwrap().values,
        NcValues::Double(vec![std::f64::consts::PI])
    );
    assert_eq!(ds, reference_classic_dataset());
}

#[test]
fn writes_classic_file_identical_to_netcdf_c() {
    let bytes = write_cdf_dataset(&reference_classic_dataset(), CdfVersion::Classic).unwrap();
    assert_eq!(&bytes[..4], b"CDF\x01");
    assert_eq!(bytes.len(), REF_CLASSIC.len());
    assert_eq!(bytes, REF_CLASSIC.to_vec());
}

#[test]
fn reads_and_writes_lone_record_variable_without_padding() {
    let ds = read_cdf_dataset(&REF_OFFSET64).unwrap();
    assert_eq!(ds.num_records(), 3);
    assert_eq!(
        ds.variable("s").unwrap().values,
        NcValues::Short(vec![1, 2, 3])
    );
    let mut built = CdfDataset::new();
    built.add_unlimited_dimension("rec", 3).unwrap();
    built
        .add_variable("s", &["rec"], NcValues::Short(vec![1, 2, 3]))
        .unwrap();
    assert_eq!(built, ds);
    let bytes = write_cdf_dataset(&built, CdfVersion::Offset64).unwrap();
    assert_eq!(&bytes[..4], b"CDF\x02");
    // 84 header bytes + 3 records of 2 bytes each, with no padding between.
    assert_eq!(bytes.len(), 90);
    assert_eq!(bytes, REF_OFFSET64.to_vec());
}

#[test]
fn reads_and_writes_cdf5_extended_types() {
    let ds = read_cdf_dataset(&REF_DATA64).unwrap();
    assert_eq!(
        ds.variable("big").unwrap().values,
        NcValues::Int64(vec![-(1_i64 << 40), (1_i64 << 40) + 7])
    );
    let ub = ds.variable("ub").unwrap();
    assert_eq!(ub.values, NcValues::UByte(vec![200, 255]));
    assert_eq!(
        ub.attribute("scale").map(|a| a.values.clone()),
        Some(NcValues::UShort(vec![65000]))
    );
    let bytes = write_cdf_dataset(&ds, CdfVersion::Data64).unwrap();
    assert_eq!(bytes, REF_DATA64.to_vec());
    // The extended types cannot be stored in the 32-bit formats.
    assert!(write_cdf_dataset(&ds, CdfVersion::Offset64).is_err());
    assert!(write_cdf_dataset(&ds, CdfVersion::Classic).is_err());
}

#[test]
fn minimal_header_layout_matches_the_specification() {
    let mut ds = CdfDataset::new();
    ds.add_dimension("x", 2).unwrap();
    ds.add_variable("v", &["x"], NcValues::Double(vec![1.0, 2.0]))
        .unwrap();
    let classic = write_cdf_dataset(&ds, CdfVersion::Classic).unwrap();
    let mut expected: Vec<u8> = Vec::new();
    expected.extend_from_slice(b"CDF\x01");
    for word in [0u32, 0x0A, 1, 1] {
        expected.extend_from_slice(&word.to_be_bytes());
    }
    expected.extend_from_slice(b"x\0\0\0");
    for word in [2u32, 0, 0, 0x0B, 1, 1] {
        expected.extend_from_slice(&word.to_be_bytes());
    }
    expected.extend_from_slice(b"v\0\0\0");
    // ndims, dimid, ABSENT attributes, NC_DOUBLE, vsize, begin.
    for word in [1u32, 0, 0, 0, 6, 16, 80] {
        expected.extend_from_slice(&word.to_be_bytes());
    }
    expected.extend_from_slice(&1.0f64.to_be_bytes());
    expected.extend_from_slice(&2.0f64.to_be_bytes());
    assert_eq!(classic, expected);

    // CDF-2 widens only the begin field: same header plus four bytes.
    let offset64 = write_cdf_dataset(&ds, CdfVersion::Offset64).unwrap();
    assert_eq!(&offset64[..4], b"CDF\x02");
    assert_eq!(offset64.len(), classic.len() + 4);
    assert_eq!(&offset64[4..76], &classic[4..76]);
    assert_eq!(&offset64[76..84], &84u64.to_be_bytes());
    assert_eq!(read_cdf_dataset(&offset64).unwrap(), ds);
}

#[test]
fn every_truncated_prefix_is_rejected() {
    // (file, end of the last data value): only the trailing padding after
    // the last value may be missing.
    for (file, data_end) in [
        (&REF_CLASSIC[..], 410),
        (&REF_OFFSET64[..], 90),
        (&REF_DATA64[..], 238),
    ] {
        for len in 0..data_end {
            assert!(
                read_cdf_dataset(&file[..len]).is_err(),
                "prefix of {len} of {} bytes was accepted",
                file.len()
            );
        }
        for len in data_end..=file.len() {
            assert_eq!(
                read_cdf_dataset(&file[..len]).unwrap(),
                read_cdf_dataset(file).unwrap()
            );
        }
    }
}

#[test]
fn corrupted_headers_are_rejected_without_panicking() {
    let mut bad_magic = REF_CLASSIC;
    bad_magic[0] = b'X';
    assert!(read_cdf_dataset(&bad_magic).is_err());
    let mut bad_version = REF_CLASSIC;
    bad_version[3] = 3;
    assert!(read_cdf_dataset(&bad_version).is_err());
    // Every single-byte corruption must be handled gracefully.
    for at in 0..REF_CLASSIC.len() {
        for value in [0x00u8, 0x7F, 0xFF] {
            let mut bytes = REF_CLASSIC;
            bytes[at] = value;
            let _ = read_cdf_dataset(&bytes);
        }
    }
    // A huge element count must not trigger a huge allocation.
    let mut huge = REF_CLASSIC;
    huge[12..16].copy_from_slice(&0x7FFF_FFF0u32.to_be_bytes());
    assert!(read_cdf_dataset(&huge).is_err());
}

#[test]
fn streaming_record_count_is_derived_from_file_length() {
    let mut bytes = REF_CLASSIC;
    bytes[4..8].copy_from_slice(&u32::MAX.to_be_bytes());
    let ds = read_cdf_dataset(&bytes).unwrap();
    assert_eq!(ds.num_records(), 3);
    assert_eq!(ds, reference_classic_dataset());
}

#[test]
fn unknown_type_code_is_rejected() {
    let mut bytes = REF_CLASSIC;
    // nc_type of the `version` global attribute (NC_INT = 4) at offset 100.
    assert_eq!(be32(&bytes, 100), 4);
    bytes[100..104].copy_from_slice(&42u32.to_be_bytes());
    assert!(read_cdf_dataset(&bytes).is_err());
    // An extended type is not allowed in a CDF-1 file.
    bytes[100..104].copy_from_slice(&10u32.to_be_bytes());
    assert!(read_cdf_dataset(&bytes).is_err());
}

#[test]
fn writer_rejects_inconsistent_datasets() {
    let mut ds = CdfDataset::new();
    ds.add_dimension("x", 3).unwrap();
    assert!(
        ds.add_variable("v", &["x"], NcValues::Double(vec![1.0, 2.0]))
            .is_err()
    );
    assert!(
        ds.add_variable("v", &["y"], NcValues::Double(vec![1.0]))
            .is_err()
    );
    assert!(ds.add_dimension("x", 4).is_err());
    assert!(ds.add_dimension("zero", 0).is_err());
    assert!(ds.add_dimension("bad/name", 1).is_err());
    ds.add_unlimited_dimension("t", 2).unwrap();
    assert!(ds.add_unlimited_dimension("t2", 2).is_err());
    assert!(
        ds.add_variable("late", &["x", "t"], NcValues::Double(vec![0.0; 6]))
            .is_err()
    );
    // A dataset mutated into an invalid state is caught by the writer.
    ds.add_variable("ok", &["t", "x"], NcValues::Double(vec![0.0; 6]))
        .unwrap();
    ds.variables[0].values = NcValues::Double(vec![0.0; 5]);
    assert!(write_cdf_dataset(&ds, CdfVersion::Offset64).is_err());
}

#[test]
fn values_convert_from_f64_with_range_checks() {
    assert_eq!(
        NcValues::from_f64(NcType::Byte, &[-128.0, 127.0]).unwrap(),
        NcValues::Byte(vec![-128, 127])
    );
    assert!(NcValues::from_f64(NcType::Byte, &[128.0]).is_err());
    assert!(NcValues::from_f64(NcType::Int, &[1.5]).is_err());
    assert!(NcValues::from_f64(NcType::Short, &[f64::NAN]).is_err());
    assert!(NcValues::from_f64(NcType::Float, &[1e300]).is_err());
    assert_eq!(
        NcValues::from_f64(NcType::Char, &[120.0, 121.0, 122.0]).unwrap(),
        NcValues::text("xyz")
    );
    assert_eq!(
        NcValues::from_f64(NcType::Int64, &[-9_223_372_036_854_775_808.0]).unwrap(),
        NcValues::Int64(vec![i64::MIN])
    );
    assert!(NcValues::from_f64(NcType::Int64, &[9_223_372_036_854_775_808.0]).is_err());
    for ty in NcType::ALL {
        assert_eq!(NcType::from_code(ty.code()), Some(ty));
        assert_eq!(NcType::from_cdl_name(ty.cdl_name()), Some(ty));
    }
}

#[test]
fn netcdf_file_round_trips_through_binary() {
    let mut file = NetCdfFile::new();
    file.add_unlimited_dimension("step", 2);
    file.add_dimension("n", 3);
    file.add_attribute("source", "unit test");
    let mut energy = NetCdfVariable::new("energy", vec!["step".into()], vec![-1.25, 3.5]);
    energy.add_attribute("units", "kJ/mol");
    file.add_variable(energy);
    file.add_variable(
        NetCdfVariable::new(
            "ids",
            vec!["step".into(), "n".into()],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
        .with_nc_type(NcType::Int),
    );
    file.add_variable(NetCdfVariable::text("tag", vec!["n".into()], "abc"));
    for version in [
        CdfVersion::Classic,
        CdfVersion::Offset64,
        CdfVersion::Data64,
    ] {
        let bytes = write_cdf(&file, version).unwrap();
        assert_eq!(&bytes[..4], &version.magic());
        let back = read_cdf(&bytes).unwrap();
        assert_eq!(back.unlimited_dim.as_deref(), Some("step"));
        assert_eq!(back.get_dimension_size("step"), Some(2));
        assert_eq!(back.get_attribute("source"), Some("unit test"));
        let energy = back.get_variable("energy").unwrap();
        assert_eq!(energy.data, vec![-1.25, 3.5]);
        assert_eq!(energy.nc_type, NcType::Double);
        assert_eq!(energy.get_attribute("units"), Some("kJ/mol"));
        let ids = back.get_variable("ids").unwrap();
        assert_eq!(ids.nc_type, NcType::Int);
        assert_eq!(ids.data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(
            back.get_variable("tag").unwrap().as_text().as_deref(),
            Some("abc")
        );
    }
    // Fractional data cannot be stored in an integer variable.
    file.variables[1].data[0] = 0.5;
    assert!(write_cdf(&file, CdfVersion::Classic).is_err());
}

#[test]
fn numeric_attributes_are_rendered_as_text_in_netcdf_file() {
    let file = read_cdf(&REF_CLASSIC).unwrap();
    assert_eq!(file.get_attribute("version"), Some("1, 2, 3"));
    assert_eq!(file.get_attribute("title"), Some("odd"));
    assert_eq!(
        file.get_variable("flags").unwrap().data,
        vec![-1.0, 0.0, 1.0]
    );
    assert_eq!(
        file.get_variable("label").unwrap().as_text().as_deref(),
        Some("hello")
    );
}

#[test]
fn binary_file_helpers_use_the_filesystem() {
    let path =
        std::env::temp_dir().join(format!("oxiphysics_cdf_helpers_{}.nc", std::process::id()));
    let file = create_amber_trajectory(1, 1, vec![1.0, 2.0, 3.0], vec![0.5], None, None);
    write_cdf_file(&path, &file, CdfVersion::Offset64).unwrap();
    let raw = std::fs::read(&path).unwrap();
    assert_eq!(&raw[..4], b"CDF\x02");
    let back = read_cdf_file(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(
        back.get_variable("coordinates").unwrap().data,
        vec![1.0, 2.0, 3.0]
    );
    assert!(read_cdf_file(&path).is_err());
}

#[test]
fn cdl_text_keeps_types_char_data_and_record_count() {
    let file = create_amber_trajectory(
        1,
        2,
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![0.0, 1.0],
        Some(vec![10.0; 6]),
        Some(vec![90.0; 6]),
    );
    let cdl = write_ncf_text(&file);
    assert!(
        cdl.contains("\tfloat coordinates(frame, atom, spatial) ;"),
        "{cdl}"
    );
    assert!(cdl.contains("\tchar spatial(spatial) ;"), "{cdl}");
    assert!(cdl.contains("\tspatial = \"xyz\" ;"), "{cdl}");
    assert!(
        cdl.contains("\tcell_angular = \"alpha\", \"beta \", \"gamma\" ;"),
        "{cdl}"
    );
    let parsed = parse_ncf_text(&cdl).unwrap();
    assert_eq!(parsed.get_dimension_size("frame"), Some(2));
    for var in &file.variables {
        let back = parsed.get_variable(&var.name).unwrap();
        assert_eq!(back.nc_type, var.nc_type, "{}", var.name);
        assert_eq!(back.data, var.data, "{}", var.name);
        assert_eq!(back.dims, var.dims, "{}", var.name);
    }
    // The parsed text converts to binary exactly like the original.
    assert_eq!(
        write_cdf(&parsed, CdfVersion::Offset64).unwrap(),
        write_cdf(&file, CdfVersion::Offset64).unwrap()
    );
}
