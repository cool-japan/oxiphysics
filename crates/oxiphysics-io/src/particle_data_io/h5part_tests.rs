// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Corrupt or hostile H5Part streams must be rejected with an error; counts
//! read from the header must never drive an allocation.

use super::{H5partReader, H5partWriter, make_test_dataset};

const MAGIC: &[u8; 8] = b"H5PART\0\0";

fn stream(parts: &[&[u8]]) -> Vec<u8> {
    let mut out = MAGIC.to_vec();
    for part in parts {
        out.extend_from_slice(part);
    }
    out
}

#[test]
fn absurd_step_count_is_an_error_not_an_abort() {
    let bytes = stream(&[&u64::MAX.to_le_bytes()]);
    assert!(H5partReader::from_bytes(&bytes).is_err());
    let bytes = stream(&[&(1u64 << 62).to_le_bytes()]);
    assert!(H5partReader::from_bytes(&bytes).is_err());
}

#[test]
fn step_count_larger_than_the_stream_is_a_truncation_error() {
    let mut writer = H5partWriter::new();
    writer.write_step(make_test_dataset(3, 0.5));
    let mut bytes = writer.to_bytes();
    bytes[8..16].copy_from_slice(&1_000_000u64.to_le_bytes());
    let err = H5partReader::from_bytes(&bytes).unwrap_err();
    assert!(err.to_string().contains("truncated"), "{err}");
}

#[test]
fn absurd_particle_count_is_an_error_not_an_abort() {
    let bytes = stream(&[
        &1u64.to_le_bytes(),
        &7u64.to_le_bytes(),
        &0.25f64.to_le_bytes(),
        &u64::MAX.to_le_bytes(),
    ]);
    assert!(H5partReader::from_bytes(&bytes).is_err());
}

#[test]
fn zero_steps_is_an_empty_reader() {
    let reader = H5partReader::from_bytes(&stream(&[&0u64.to_le_bytes()])).unwrap();
    assert_eq!(reader.num_steps(), 0);
}

#[test]
fn honest_streams_still_round_trip() {
    let mut writer = H5partWriter::new();
    writer.write_step(make_test_dataset(5, 0.5));
    writer.write_step(make_test_dataset(2, 1.0));
    let reader = H5partReader::from_bytes(&writer.to_bytes()).unwrap();
    assert_eq!(reader.num_steps(), 2);
    assert_eq!(reader.step(0).unwrap().dataset.len(), 5);
    assert_eq!(reader.step(1).unwrap().dataset.len(), 2);
}
