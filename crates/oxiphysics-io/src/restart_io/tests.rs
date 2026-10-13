// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the restart/checkpoint module.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;

// ── Helpers ───────────────────────────────────────────────────────────────────

fn make_data(n: usize) -> RestartData {
    let meta = RestartMetadata::new("1.0", 1234567890, 42, 4.2, "oxiphysics-md", "unit test");
    let positions = (0..n).map(|i| [i as f64, i as f64 * 0.5, 0.0]).collect();
    let velocities = (0..n).map(|i| [0.1 * i as f64, 0.0, 0.0]).collect();
    let forces = (0..n).map(|_| [0.0, -9.81, 0.0]).collect();
    let masses = (0..n).map(|_| 1.0).collect();
    let types = (0..n).map(|i| (i % 3) as u32).collect();
    let box_matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
    RestartData {
        metadata: meta,
        positions,
        velocities,
        forces,
        masses,
        types,
        box_matrix,
        extra_scalars: Vec::new(),
        extra_vectors: Vec::new(),
    }
}

/// `n` particles with only `positions` (and the box) populated.
fn positions_only(n: usize) -> RestartData {
    let mut d = make_data(n);
    d.velocities.clear();
    d.forces.clear();
    d.masses.clear();
    d.types.clear();
    d
}

/// A snapshot that exercises every field, including the extras and awkward
/// metadata strings.
fn full_data() -> RestartData {
    let mut d = make_data(3);
    d.metadata.description =
        "He said \"hi\"\nline two\t\\ back\\slash # not a comment ".to_string();
    d.metadata.crate_name = "oxi physics\r\n".to_string();
    d.extra_scalars = vec![
        ("charge".to_string(), vec![0.1, -0.2, 0.3]),
        ("name with spaces\nand a newline".to_string(), vec![1e-300]),
        (String::new(), Vec::new()),
    ];
    d.extra_vectors = vec![
        ("spin".to_string(), vec![[0.0, 0.0, 1.0], [1.0, -0.0, 0.0]]),
        ("empty".to_string(), Vec::new()),
    ];
    d
}

fn all_formats() -> Vec<RestartFormat> {
    vec![
        RestartFormat::Binary,
        RestartFormat::Ascii,
        RestartFormat::Hdf5Like,
        RestartFormat::Json,
        RestartFormat::MessagePack,
    ]
}

/// A uniquely named directory under the OS temp dir, removed on drop.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "oxiphysics_restart_{tag}_{}_{id}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn bits3(v: &[[f64; 3]]) -> Vec<[u64; 3]> {
    v.iter()
        .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
        .collect()
}

fn bits1(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

/// Assert that two snapshots are identical down to the bit pattern of every
/// `f64` (so `-0.0` and `0.0` are distinguished).
fn assert_bit_exact(a: &RestartData, b: &RestartData) {
    assert_eq!(a.metadata.version, b.metadata.version);
    assert_eq!(a.metadata.timestamp, b.metadata.timestamp);
    assert_eq!(a.metadata.step, b.metadata.step);
    assert_eq!(a.metadata.time.to_bits(), b.metadata.time.to_bits());
    assert_eq!(a.metadata.crate_name, b.metadata.crate_name);
    assert_eq!(a.metadata.description, b.metadata.description);
    assert_eq!(bits3(&a.positions), bits3(&b.positions));
    assert_eq!(bits3(&a.velocities), bits3(&b.velocities));
    assert_eq!(bits3(&a.forces), bits3(&b.forces));
    assert_eq!(bits1(&a.masses), bits1(&b.masses));
    assert_eq!(a.types, b.types);
    assert_eq!(bits3(&a.box_matrix), bits3(&b.box_matrix));
    assert_eq!(a.extra_scalars.len(), b.extra_scalars.len());
    for ((na, va), (nb, vb)) in a.extra_scalars.iter().zip(&b.extra_scalars) {
        assert_eq!(na, nb);
        assert_eq!(bits1(va), bits1(vb));
    }
    assert_eq!(a.extra_vectors.len(), b.extra_vectors.len());
    for ((na, va), (nb, vb)) in a.extra_vectors.iter().zip(&b.extra_vectors) {
        assert_eq!(na, nb);
        assert_eq!(bits3(va), bits3(vb));
    }
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Finite `f64`s with uniformly random bit patterns (normals, subnormals, huge
/// and tiny exponents, both signs).
fn random_finite_f64s(count: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    let mut out = Vec::with_capacity(count);
    while out.len() < count {
        let v = f64::from_bits(splitmix64(&mut state));
        if v.is_finite() {
            out.push(v);
        }
    }
    out
}

/// A snapshot whose floating-point fields all hold random-bit-pattern values.
fn random_float_data(n: usize, seed: u64) -> RestartData {
    let mut values = random_finite_f64s(n * 10 + 20, seed).into_iter();
    let mut next = || values.next().unwrap();
    let mut d = make_data(n);
    for i in 0..n {
        d.positions[i] = [next(), next(), next()];
        d.velocities[i] = [next(), next(), next()];
        d.masses[i] = next();
    }
    for i in 0..n {
        d.forces[i] = [next(), next(), next()];
    }
    d.metadata.time = next();
    d.box_matrix = [
        [next(), next(), next()],
        [next(), next(), next()],
        [next(), next(), next()],
    ];
    d
}

fn file_roundtrip(format: RestartFormat, data: &RestartData, tag: &str) -> RestartData {
    let dir = TempDir::new(tag);
    let path = dir.path("restart.dat");
    RestartWriter::new(&path, format).write(data).unwrap();
    RestartReader::new(&path).read().unwrap()
}

/// Offset of the `N_PARTICLES` field in a binary blob written for `meta`.
fn count_offset(meta: &RestartMetadata) -> usize {
    8 + (4 + meta.version.len()) + 24 + (4 + meta.crate_name.len()) + (4 + meta.description.len())
}

/// The version-1 binary layout (no presence byte, all arrays always present),
/// kept here to prove the reader still accepts files written by old code.
fn encode_version_1(d: &RestartData) -> Vec<u8> {
    fn u32le(buf: &mut Vec<u8>, v: u32) {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    fn u64le(buf: &mut Vec<u8>, v: u64) {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    fn f64le(buf: &mut Vec<u8>, v: f64) {
        buf.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    fn s(buf: &mut Vec<u8>, v: &str) {
        u32le(buf, v.len() as u32);
        buf.extend_from_slice(v.as_bytes());
    }
    let mut b = Vec::new();
    b.extend_from_slice(b"OXRS");
    u32le(&mut b, 1);
    s(&mut b, &d.metadata.version);
    u64le(&mut b, d.metadata.timestamp);
    u64le(&mut b, d.metadata.step);
    f64le(&mut b, d.metadata.time);
    s(&mut b, &d.metadata.crate_name);
    s(&mut b, &d.metadata.description);
    u64le(&mut b, d.positions.len() as u64);
    for arr in [&d.positions, &d.velocities, &d.forces] {
        for p in arr.iter() {
            for &c in p {
                f64le(&mut b, c);
            }
        }
    }
    for &m in &d.masses {
        f64le(&mut b, m);
    }
    for &t in &d.types {
        u32le(&mut b, t);
    }
    for row in &d.box_matrix {
        for &c in row {
            f64le(&mut b, c);
        }
    }
    u32le(&mut b, d.extra_scalars.len() as u32);
    for (name, vals) in &d.extra_scalars {
        s(&mut b, name);
        u64le(&mut b, vals.len() as u64);
        for &v in vals {
            f64le(&mut b, v);
        }
    }
    u32le(&mut b, d.extra_vectors.len() as u32);
    for (name, vecs) in &d.extra_vectors {
        s(&mut b, name);
        u64le(&mut b, vecs.len() as u64);
        for v in vecs {
            for &c in v {
                f64le(&mut b, c);
            }
        }
    }
    b
}

// ── RestartMetadata / RestartData ─────────────────────────────────────────────

#[test]
fn test_metadata_new_fields() {
    let m = RestartMetadata::new("2.0", 100, 5, 0.5, "crate_x", "desc");
    assert_eq!(m.version, "2.0");
    assert_eq!(m.timestamp, 100);
    assert_eq!(m.step, 5);
    assert!((m.time - 0.5).abs() < 1e-12);
    assert_eq!(m.crate_name, "crate_x");
    assert_eq!(m.description, "desc");
}

#[test]
fn test_metadata_default_test() {
    let m = RestartMetadata::default_test();
    assert_eq!(m.version, "1.0");
    assert_eq!(m.step, 0);
}

#[test]
fn test_restart_data_n_particles() {
    let d = make_data(7);
    assert_eq!(d.n_particles(), 7);
}

#[test]
fn test_restart_data_empty() {
    let d = RestartData::empty(RestartMetadata::default_test());
    assert_eq!(d.n_particles(), 0);
    assert!(d.positions.is_empty());
}

#[test]
fn test_restart_data_single_particle_test() {
    let d = RestartData::single_particle_test();
    assert_eq!(d.n_particles(), 1);
    assert!((d.positions[0][0] - 1.0).abs() < 1e-12);
}

#[test]
fn validate_accepts_full_absent_and_empty_arrays() {
    assert!(make_data(5).validate().is_ok());
    assert!(positions_only(5).validate().is_ok());
    assert!(
        RestartData::empty(RestartMetadata::default_test())
            .validate()
            .is_ok()
    );
    let mut partial = make_data(4);
    partial.forces.clear();
    partial.types.clear();
    assert!(partial.validate().is_ok());
}

#[test]
fn validate_rejects_arrays_that_disagree_with_the_particle_count() {
    let mut d = make_data(4);
    d.masses.pop();
    let err = d.validate().unwrap_err();
    assert!(
        err.contains("masses") && err.contains("3") && err.contains("4"),
        "{err}"
    );

    let mut d = make_data(4);
    d.velocities.push([0.0; 3]);
    assert!(d.validate().unwrap_err().contains("velocities"));

    let mut d = make_data(4);
    d.types.truncate(1);
    assert!(d.validate().unwrap_err().contains("types"));

    let mut d = make_data(4);
    d.forces.truncate(2);
    assert!(d.validate().unwrap_err().contains("forces"));

    // Per-particle arrays without any positions cannot be stored under n = 0.
    let mut d = positions_only(0);
    d.masses = vec![1.0];
    assert!(d.validate().is_err());
}

// ── Binary round-trip ─────────────────────────────────────────────────────────

#[test]
fn test_binary_roundtrip_zero_particles() {
    let d = make_data(0);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let recovered = RestartReader::read_binary(&bytes).unwrap();
    assert_eq!(recovered.n_particles(), 0);
    assert_eq!(recovered.metadata.step, 42);
    assert_eq!(recovered, d);
}

#[test]
fn test_binary_roundtrip_positions() {
    let d = make_data(5);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    for i in 0..5 {
        assert!((r.positions[i][0] - d.positions[i][0]).abs() < 1e-12);
        assert!((r.positions[i][1] - d.positions[i][1]).abs() < 1e-12);
    }
    assert_eq!(r, d);
}

#[test]
fn test_binary_roundtrip_metadata() {
    let d = make_data(3);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    assert_eq!(r.metadata.version, "1.0");
    assert_eq!(r.metadata.timestamp, 1234567890);
    assert_eq!(r.metadata.step, 42);
    assert!((r.metadata.time - 4.2).abs() < 1e-10);
    assert_eq!(r.metadata.crate_name, "oxiphysics-md");
}

#[test]
fn test_binary_roundtrip_box_matrix() {
    let d = make_data(2);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    assert!((r.box_matrix[0][0] - 10.0).abs() < 1e-12);
    assert!((r.box_matrix[1][1] - 10.0).abs() < 1e-12);
    assert!((r.box_matrix[2][2] - 10.0).abs() < 1e-12);
    assert!(r.box_matrix[0][1].abs() < 1e-12);
}

#[test]
fn test_binary_roundtrip_types() {
    let d = make_data(6);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    for i in 0..6 {
        assert_eq!(r.types[i], (i % 3) as u32);
    }
}

#[test]
fn test_binary_roundtrip_extra_scalars() {
    let mut d = make_data(3);
    d.extra_scalars
        .push(("charge".into(), vec![0.1, -0.2, 0.3]));
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    assert_eq!(r.extra_scalars.len(), 1);
    assert_eq!(r.extra_scalars[0].0, "charge");
    assert!((r.extra_scalars[0].1[1] - (-0.2)).abs() < 1e-12);
}

#[test]
fn test_binary_roundtrip_extra_vectors() {
    let mut d = make_data(2);
    d.extra_vectors
        .push(("spin".into(), vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]));
    let bytes = RestartWriter::write_binary(&d).unwrap();
    let r = RestartReader::read_binary(&bytes).unwrap();
    assert_eq!(r.extra_vectors.len(), 1);
    assert_eq!(r.extra_vectors[0].0, "spin");
    assert!((r.extra_vectors[0].1[0][2] - 1.0).abs() < 1e-12);
}

#[test]
fn test_binary_magic_bytes() {
    let d = make_data(1);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    assert_eq!(&bytes[0..4], b"OXRS");
    assert_eq!(&RestartWriter::write_hdf5like(&d).unwrap()[0..4], b"OXH5");
    assert_eq!(&RestartWriter::write_msgpack(&d).unwrap()[0..4], b"MPK1");
}

#[test]
fn binary_roundtrip_is_bit_exact_for_random_floats_and_extras() {
    let mut d = random_float_data(40, 0xC0FFEE);
    d.extra_scalars = vec![("s".into(), random_finite_f64s(17, 5))];
    d.extra_vectors = vec![(
        "v".into(),
        random_finite_f64s(9, 6)
            .chunks(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect(),
    )];
    let r = RestartReader::read_binary(&RestartWriter::write_binary(&d).unwrap()).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn binary_roundtrip_preserves_nan_payloads_and_infinities() {
    let mut d = make_data(1);
    let nan_with_payload = f64::from_bits(0x7FF8_0000_DEAD_BEEF);
    d.positions[0] = [nan_with_payload, f64::INFINITY, f64::NEG_INFINITY];
    let r = RestartReader::read_binary(&RestartWriter::write_binary(&d).unwrap()).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn binary_roundtrip_every_tag_and_unicode_strings() {
    let mut d = full_data();
    d.metadata.description = "résumé — 日本語 🚀".to_string();
    for bytes in [
        RestartWriter::write_binary(&d).unwrap(),
        RestartWriter::write_hdf5like(&d).unwrap(),
        RestartWriter::write_msgpack(&d).unwrap(),
    ] {
        assert_bit_exact(&d, &RestartReader::read_binary(&bytes).unwrap());
    }
}

// F2437: arrays of different lengths under one particle count.

#[test]
fn restart_positions_only_roundtrips_in_every_format() {
    let d = positions_only(4);
    for format in all_formats() {
        let r = file_roundtrip(format.clone(), &d, "positions_only");
        assert_eq!(r, d, "format {format:?}");
        assert!(r.velocities.is_empty() && r.forces.is_empty());
        assert!(r.masses.is_empty() && r.types.is_empty());
    }
}

#[test]
fn restart_partial_arrays_roundtrip_in_every_format() {
    // velocities and types present, forces and masses absent.
    let mut d = make_data(5);
    d.forces.clear();
    d.masses.clear();
    for format in all_formats() {
        let r = file_roundtrip(format.clone(), &d, "partial");
        assert_eq!(r, d, "format {format:?}");
    }
}

#[test]
fn restart_inconsistent_arrays_are_rejected_by_every_encoder() {
    let mut d = make_data(4);
    d.masses.pop(); // masses of length n - 1
    assert!(RestartWriter::write_binary(&d).is_err());
    assert!(RestartWriter::write_hdf5like(&d).is_err());
    assert!(RestartWriter::write_msgpack(&d).is_err());
    assert!(RestartWriter::write_ascii(&d).is_err());
    assert!(RestartWriter::write_json(&d).is_err());
    for format in all_formats() {
        let writer = RestartWriter::new("never-created.dat", format);
        let err = writer.encode(&d).unwrap_err();
        assert!(err.contains("masses"), "{err}");
    }
}

#[test]
fn binary_presence_mask_controls_which_arrays_are_stored() {
    let n = 6;
    let full = RestartWriter::write_binary(&make_data(n)).unwrap();
    let slim = RestartWriter::write_binary(&positions_only(n)).unwrap();
    // velocities (24) + forces (24) + masses (8) + types (4) bytes per particle.
    assert_eq!(full.len() - slim.len(), n * (24 + 24 + 8 + 4));

    let meta = make_data(n).metadata;
    let mask_at = count_offset(&meta) + 8;
    assert_eq!(full[mask_at], 0b1111);
    assert_eq!(slim[mask_at], 0b0000);
    let mut velocities_and_masses = make_data(n);
    velocities_and_masses.forces.clear();
    velocities_and_masses.types.clear();
    let bytes = RestartWriter::write_binary(&velocities_and_masses).unwrap();
    assert_eq!(bytes[mask_at], 0b0101);
    let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    assert_eq!(version, 2);
}

#[test]
fn binary_reader_still_accepts_version_1_files() {
    let mut d = full_data();
    d.metadata.description = "legacy".to_string();
    let v1 = encode_version_1(&d);
    let r = RestartReader::read_binary(&v1).unwrap();
    assert_bit_exact(&d, &r);
    // And through the file-based reader, which sniffs the magic.
    let dir = TempDir::new("v1");
    let path = dir.path("old.bin");
    std::fs::write(&path, &v1).unwrap();
    assert_bit_exact(&d, &RestartReader::new(&path).read().unwrap());
}

#[test]
fn incremental_delta_with_absent_arrays_roundtrips() {
    let mut full = make_data(8);
    full.velocities.clear();
    full.types.clear();
    let mut inc = IncrementalRestart::new();
    inc.mark_changed(2);
    inc.mark_changed(5);
    inc.mark_changed(99); // out of range: skipped everywhere
    let delta = inc.extract_delta(&full);
    assert_eq!(delta.n_particles(), 2);
    assert!(delta.velocities.is_empty() && delta.types.is_empty());
    assert_eq!(delta.forces.len(), 2);
    assert_eq!(delta.masses.len(), 2);
    delta.validate().unwrap();
    for format in all_formats() {
        let r = file_roundtrip(format, &delta, "delta");
        assert_eq!(r, delta);
    }
}

// F2442: unchecked counts in the binary reader.

#[test]
fn restart_read_binary_huge_n_returns_err() {
    let d = make_data(0);
    let template = RestartWriter::write_binary(&d).unwrap();
    let at = count_offset(&d.metadata);
    for count in [
        u64::MAX,
        u64::MAX / 24,
        1u64 << 58,
        (1u64 << 58) + 7,
        1u64 << 40,
        5_000_000,
        1,
    ] {
        let mut bytes = template.clone();
        bytes[at..at + 8].copy_from_slice(&count.to_le_bytes());
        let err = RestartReader::read_binary(&bytes).unwrap_err();
        assert!(
            err.contains("particle arrays") || err.contains("EOF"),
            "count {count}: {err}"
        );
    }
}

#[test]
fn restart_read_binary_bad_magic_returns_err() {
    let bytes = RestartWriter::write_binary(&make_data(2)).unwrap();
    for bad in [b"XXXX", b"OXHS", b"\0\0\0\0", b"oxrs"] {
        let mut b = bytes.clone();
        b[0..4].copy_from_slice(bad);
        let err = RestartReader::read_binary(&b).unwrap_err();
        assert!(err.contains("magic"), "{err}");
    }
}

#[test]
fn binary_reader_rejects_unknown_versions_and_presence_bits() {
    let d = make_data(2);
    let template = RestartWriter::write_binary(&d).unwrap();
    for version in [0u32, 3, u32::MAX] {
        let mut b = template.clone();
        b[4..8].copy_from_slice(&version.to_le_bytes());
        let err = RestartReader::read_binary(&b).unwrap_err();
        assert!(err.contains("version"), "{err}");
    }
    let mut b = template.clone();
    b[count_offset(&d.metadata) + 8] = 0b1_0000;
    let err = RestartReader::read_binary(&b).unwrap_err();
    assert!(err.contains("presence"), "{err}");
}

#[test]
fn binary_reader_bounds_extra_array_counts() {
    let mut d = make_data(1);
    d.extra_scalars = vec![("a".into(), vec![1.0])];
    d.extra_vectors = vec![("b".into(), vec![[1.0, 2.0, 3.0]])];
    let template = RestartWriter::write_binary(&d).unwrap();
    // Layout tail: ... [n_scalar u32]["a"][len u64][1 f64][n_vec u32]["b"][len u64][3 f64]
    let name_bytes = 4 + 1; // u32 byte length + the one-byte name
    let tail_vec = 4 + name_bytes + 8 + 24; // n_vec + name + len + 3 values
    let tail_scalar = 4 + name_bytes + 8 + 8; // n_scalar + name + len + 1 value
    let n_scalar_at = template.len() - tail_vec - tail_scalar;
    let n_vec_at = template.len() - tail_vec;

    // Overwrite `field` (a little-endian integer) in a copy of the template and
    // return the reader's error.
    let corrupt = |at: usize, field: &[u8]| -> String {
        let mut b = template.clone();
        b[at..at + field.len()].copy_from_slice(field);
        RestartReader::read_binary(&b).unwrap_err()
    };

    let err = corrupt(n_scalar_at, &u32::MAX.to_le_bytes());
    assert!(err.contains("extra scalar arrays"), "{err}");
    let err = corrupt(n_vec_at, &u32::MAX.to_le_bytes());
    assert!(err.contains("extra vector arrays"), "{err}");

    // Per-array element counts.
    let scalar_len_at = n_scalar_at + 4 + name_bytes;
    let err = corrupt(scalar_len_at, &u64::MAX.to_le_bytes());
    assert!(err.contains("extra scalar values"), "{err}");
    let vec_len_at = n_vec_at + 4 + name_bytes;
    let err = corrupt(vec_len_at, &(1u64 << 59).to_le_bytes());
    assert!(err.contains("extra vector values"), "{err}");
    // An untouched copy still parses.
    assert_bit_exact(&d, &RestartReader::read_binary(&template).unwrap());
}

#[test]
fn binary_reader_rejects_every_proper_prefix() {
    let bytes = RestartWriter::write_binary(&full_data()).unwrap();
    for len in 0..bytes.len() {
        assert!(
            RestartReader::read_binary(&bytes[..len]).is_err(),
            "prefix of {len} bytes must not parse"
        );
    }
    assert!(RestartReader::read_binary(&bytes).is_ok());
}

#[test]
fn binary_reader_survives_every_single_byte_corruption() {
    // No byte pattern may panic or abort; Ok and Err are both acceptable.
    let bytes = RestartWriter::write_binary(&full_data()).unwrap();
    for i in 0..bytes.len() {
        for fill in [0x00u8, 0x7F, 0x80, 0xFF] {
            let mut b = bytes.clone();
            b[i] = fill;
            let _ = RestartReader::read_binary(&b);
        }
    }
}

#[test]
fn binary_reader_rejects_invalid_utf8_strings() {
    let d = make_data(0);
    let mut bytes = RestartWriter::write_binary(&d).unwrap();
    // The first metadata string ("1.0") starts right after magic + version +
    // its 4-byte length.
    bytes[12] = 0xFF;
    let err = RestartReader::read_binary(&bytes).unwrap_err();
    assert!(err.contains("UTF-8"), "{err}");
}

// ── ASCII ─────────────────────────────────────────────────────────────────────

#[test]
fn test_ascii_roundtrip_basic() {
    let d = make_data(4);
    let text = RestartWriter::write_ascii(&d).unwrap();
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_eq!(r.n_particles(), 4);
    assert_eq!(r.metadata.step, 42);
    assert_eq!(r, d);
}

#[test]
fn test_ascii_roundtrip_positions() {
    let d = make_data(3);
    let text = RestartWriter::write_ascii(&d).unwrap();
    let r = RestartReader::read_ascii(&text).unwrap();
    for i in 0..3 {
        assert!((r.positions[i][0] - d.positions[i][0]).abs() < 1e-12);
    }
}

#[test]
fn test_ascii_roundtrip_box() {
    let d = make_data(1);
    let text = RestartWriter::write_ascii(&d).unwrap();
    let r = RestartReader::read_ascii(&text).unwrap();
    assert!((r.box_matrix[0][0] - 10.0).abs() < 1e-12);
}

#[test]
fn test_ascii_contains_keywords() {
    let d = make_data(2);
    let text = RestartWriter::write_ascii(&d).unwrap();
    assert!(text.contains("BEGIN_POSITIONS"));
    assert!(text.contains("END_POSITIONS"));
    assert!(text.contains("VERSION"));
    assert!(text.contains("STEP"));
    assert!(text.contains("N_PARTICLES 2"));
}

// F2439: lossy 6-decimal formatting and silent parse failures.

#[test]
fn restart_ascii_roundtrip_bit_exact() {
    let mut d = make_data(2);
    d.positions = vec![
        [1e-12, 1e300, -0.0],
        [f64::MIN_POSITIVE, -1e-300, 0.1 + 0.2],
    ];
    d.velocities = vec![
        [5e-324, f64::MAX, -f64::MIN_POSITIVE],
        [1.0 / 3.0, -2.5e-9, 0.0],
    ];
    d.masses = vec![1.6605390666e-27, 9.1093837015e-31];
    d.metadata.time = 1e-15;
    d.box_matrix = [[2.5e-9, 0.0, -0.0], [0.0, 2.5e-9, 0.0], [0.0, 0.0, 2.5e-9]];
    let text = RestartWriter::write_ascii(&d).unwrap();
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_bit_exact(&d, &r);
    assert_eq!(r.positions[0][2].to_bits(), (-0.0f64).to_bits());
}

#[test]
fn restart_ascii_nanoscale_values_are_not_rounded_to_zero() {
    let mut d = make_data(1);
    d.positions = vec![[1e-9, 2e-9, 3e-9]];
    d.metadata.time = 1e-12;
    let text = RestartWriter::write_ascii(&d).unwrap();
    assert!(!text.contains("0.000000"), "{text}");
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_eq!(r.positions[0], [1e-9, 2e-9, 3e-9]);
    assert_eq!(r.metadata.time, 1e-12);
}

#[test]
fn ascii_roundtrip_is_bit_exact_for_random_floats() {
    for seed in [1u64, 2, 3] {
        let d = random_float_data(60, seed);
        let r = RestartReader::read_ascii(&RestartWriter::write_ascii(&d).unwrap()).unwrap();
        assert_bit_exact(&d, &r);
    }
}

#[test]
fn ascii_roundtrip_free_text_and_extras() {
    let d = full_data();
    let text = RestartWriter::write_ascii(&d).unwrap();
    // Every header value stays on one line.
    assert!(text.lines().any(|l| l.starts_with("DESC He said")));
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn ascii_roundtrip_nonfinite_values() {
    let mut d = make_data(1);
    d.positions[0] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    let r = RestartReader::read_ascii(&RestartWriter::write_ascii(&d).unwrap()).unwrap();
    assert!(r.positions[0][0].is_nan());
    assert_eq!(r.positions[0][1], f64::INFINITY);
    assert_eq!(r.positions[0][2], f64::NEG_INFINITY);
}

#[test]
fn ascii_reader_accepts_crlf_line_endings() {
    let d = full_data();
    let text = RestartWriter::write_ascii(&d)
        .unwrap()
        .replace('\n', "\r\n");
    // The escaped description never contains a raw CR, so CRLF is only ever a
    // line terminator here.
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn ascii_reader_accepts_files_from_the_six_decimal_writer() {
    let legacy = "# OxiPhysics restart file\nVERSION 1.0\nTIMESTAMP 7\nSTEP 3\nTIME 0.500000\nCRATE old\nDESC legacy file\nN_PARTICLES 1\nBEGIN_POSITIONS\n1.000000 2.000000 3.000000\nEND_POSITIONS\nBEGIN_VELOCITIES\n0.100000 0.200000 0.300000\nEND_VELOCITIES\nBEGIN_FORCES\n0.000000 -9.810000 0.000000\nEND_FORCES\nBEGIN_MASSES\n1.000000\nEND_MASSES\nBEGIN_TYPES\n0\nEND_TYPES\nBEGIN_BOX\n10.000000 0.000000 0.000000\n0.000000 10.000000 0.000000\n0.000000 0.000000 10.000000\nEND_BOX\n";
    let r = RestartReader::read_ascii(legacy).unwrap();
    assert_eq!(r.metadata.step, 3);
    assert_eq!(r.metadata.time, 0.5);
    assert_eq!(r.metadata.description, "legacy file");
    assert_eq!(r.positions, vec![[1.0, 2.0, 3.0]]);
    assert_eq!(r.forces, vec![[0.0, -9.81, 0.0]]);
    assert_eq!(r.box_matrix[2][2], 10.0);
}

/// A minimal valid ASCII file with `positions` lines substituted.
fn ascii_with(body: &str) -> String {
    format!("N_PARTICLES 1\n{body}")
}

#[test]
fn ascii_reader_reports_bad_numbers_with_the_line_number() {
    let text = ascii_with("BEGIN_POSITIONS\n1.0 2.0 oops\nEND_POSITIONS\n");
    let err = RestartReader::read_ascii(&text).unwrap_err();
    assert!(err.contains("line 3") && err.contains("oops"), "{err}");

    let text = "N_PARTICLES 1\nSTEP twelve\n";
    let err = RestartReader::read_ascii(text).unwrap_err();
    assert!(err.contains("line 2") && err.contains("STEP"), "{err}");

    let text = "N_PARTICLES many\n";
    assert!(
        RestartReader::read_ascii(text)
            .unwrap_err()
            .contains("line 1")
    );
    let text = "TIME 1.0.0\nN_PARTICLES 0\n";
    assert!(
        RestartReader::read_ascii(text)
            .unwrap_err()
            .contains("TIME")
    );
    let text = "TIMESTAMP -4\nN_PARTICLES 0\n";
    assert!(
        RestartReader::read_ascii(text)
            .unwrap_err()
            .contains("TIMESTAMP")
    );
}

#[test]
fn ascii_reader_rejects_unparsable_masses_and_types_instead_of_skipping_them() {
    let d = make_data(3);
    let text = RestartWriter::write_ascii(&d).unwrap();
    // The old reader dropped the bad mass, leaving 2 masses for 3 particles.
    let bad_mass = text.replacen("BEGIN_MASSES\n1.0\n", "BEGIN_MASSES\nnot-a-number\n", 1);
    assert_ne!(bad_mass, text);
    let err = RestartReader::read_ascii(&bad_mass).unwrap_err();
    assert!(
        err.contains("mass") && err.contains("not-a-number"),
        "{err}"
    );

    let bad_type = text.replacen("BEGIN_TYPES\n0\n", "BEGIN_TYPES\n-1\n", 1);
    assert_ne!(bad_type, text);
    let err = RestartReader::read_ascii(&bad_type).unwrap_err();
    assert!(err.contains("type"), "{err}");
}

#[test]
fn ascii_reader_rejects_wrong_column_counts() {
    for row in ["1.0 2.0", "1 2 3 4", "1.0"] {
        let text = ascii_with(&format!("BEGIN_POSITIONS\n{row}\nEND_POSITIONS\n"));
        let err = RestartReader::read_ascii(&text).unwrap_err();
        assert!(err.contains("exactly 3 values"), "{row}: {err}");
    }
    let text =
        ascii_with("BEGIN_POSITIONS\n1 2 3\nEND_POSITIONS\nBEGIN_MASSES\n1.0 2.0\nEND_MASSES\n");
    assert!(
        RestartReader::read_ascii(&text)
            .unwrap_err()
            .contains("exactly 1 value")
    );
}

#[test]
fn ascii_reader_rejects_arrays_that_disagree_with_n_particles() {
    let d = make_data(3);
    let text = RestartWriter::write_ascii(&d).unwrap();

    // Truncated positions (one line removed).
    let truncated = text.replacen("2.0 1.0 0.0\n", "", 1);
    assert_ne!(truncated, text);
    let err = RestartReader::read_ascii(&truncated).unwrap_err();
    assert!(
        err.contains("N_PARTICLES is 3") && err.contains("POSITIONS"),
        "{err}"
    );

    // Masses one short.
    let short_masses = text.replacen("BEGIN_MASSES\n1.0\n", "BEGIN_MASSES\n", 1);
    assert_ne!(short_masses, text);
    let err = RestartReader::read_ascii(&short_masses).unwrap_err();
    assert!(
        err.contains("MASSES") && err.contains("expected 0 or 3"),
        "{err}"
    );

    // Declared count larger than the data.
    let wrong_n = text.replacen("N_PARTICLES 3", "N_PARTICLES 4", 1);
    assert!(RestartReader::read_ascii(&wrong_n).is_err());
}

#[test]
fn ascii_reader_structural_errors() {
    let missing_n = "BEGIN_POSITIONS\nEND_POSITIONS\n";
    assert!(
        RestartReader::read_ascii(missing_n)
            .unwrap_err()
            .contains("missing N_PARTICLES")
    );

    let unterminated = ascii_with("BEGIN_POSITIONS\n1 2 3\n");
    let err = RestartReader::read_ascii(&unterminated).unwrap_err();
    assert!(
        err.contains("line 2") && err.contains("END_POSITIONS"),
        "{err}"
    );

    let wrong_end = ascii_with("BEGIN_POSITIONS\n1 2 3\nEND_VELOCITIES\n");
    let err = RestartReader::read_ascii(&wrong_end).unwrap_err();
    assert!(
        err.contains("line 4") && err.contains("unexpected"),
        "{err}"
    );

    let nested = ascii_with("BEGIN_POSITIONS\nBEGIN_VELOCITIES\n");
    assert!(
        RestartReader::read_ascii(&nested)
            .unwrap_err()
            .contains("unexpected")
    );

    let unknown = ascii_with("BOGUS 1\n");
    let err = RestartReader::read_ascii(&unknown).unwrap_err();
    assert!(
        err.contains("line 2") && err.contains("unrecognised"),
        "{err}"
    );

    let stray_end = ascii_with("END_POSITIONS\n");
    assert!(
        RestartReader::read_ascii(&stray_end)
            .unwrap_err()
            .contains("unrecognised")
    );

    let duplicate = ascii_with(
        "BEGIN_POSITIONS\n1 2 3\nEND_POSITIONS\nBEGIN_POSITIONS\n4 5 6\nEND_POSITIONS\n",
    );
    let err = RestartReader::read_ascii(&duplicate).unwrap_err();
    assert!(err.contains("duplicate"), "{err}");

    let short_box =
        ascii_with("BEGIN_POSITIONS\n1 2 3\nEND_POSITIONS\nBEGIN_BOX\n1 0 0\n0 1 0\nEND_BOX\n");
    let err = RestartReader::read_ascii(&short_box).unwrap_err();
    assert!(err.contains("BOX") && err.contains("3 rows"), "{err}");
}

#[test]
fn ascii_reader_ignores_a_utf8_byte_order_mark() {
    // Editors on some platforms prepend a BOM; it must not turn the first line
    // into an "unrecognised line" error.
    let d = full_data();
    let text = format!("\u{feff}{}", RestartWriter::write_ascii(&d).unwrap());
    assert_bit_exact(&d, &RestartReader::read_ascii(&text).unwrap());

    let dir = TempDir::new("ascii_bom");
    let path = dir.path("restart.txt");
    std::fs::write(&path, text.as_bytes()).unwrap();
    assert_bit_exact(&d, &RestartReader::new(&path).read().unwrap());
}

#[test]
fn ascii_reader_survives_every_single_byte_corruption() {
    // No byte pattern may panic; Ok and Err are both acceptable.
    let text = RestartWriter::write_ascii(&full_data()).unwrap();
    for i in 0..text.len() {
        for fill in [b'\n', b' ', b'#', b'-', b'e', b'9', b'\\', b'\0'] {
            let mut bytes = text.clone().into_bytes();
            bytes[i] = fill;
            if let Ok(mutated) = String::from_utf8(bytes) {
                let _ = RestartReader::read_ascii(&mutated);
            }
        }
    }
}

#[test]
fn ascii_reader_defaults_absent_optional_sections() {
    let text = ascii_with("BEGIN_POSITIONS\n1 2 3\nEND_POSITIONS\n");
    let r = RestartReader::read_ascii(&text).unwrap();
    assert_eq!(r.positions, vec![[1.0, 2.0, 3.0]]);
    assert!(r.velocities.is_empty() && r.forces.is_empty());
    assert!(r.masses.is_empty() && r.types.is_empty());
    assert_eq!(r.box_matrix, [[0.0; 3]; 3]);
    assert_eq!(r.metadata.version, "");
}

// ── JSON ──────────────────────────────────────────────────────────────────────

#[test]
fn test_json_contains_particles_key() {
    let d = make_data(3);
    let j = RestartWriter::write_json(&d).unwrap();
    assert!(j.contains("n_particles"));
    assert!(j.contains("positions"));
    assert!(j.contains("velocities"));
}

#[test]
fn test_json_step_present() {
    let d = make_data(1);
    let j = RestartWriter::write_json(&d).unwrap();
    assert!(j.contains("\"step\": 42"));
}

#[test]
fn json_writer_emits_every_field_as_valid_json() {
    let d = full_data();
    let text = RestartWriter::write_json(&d).unwrap();
    let doc: serde_json::Value = serde_json::from_str(&text).unwrap();
    let obj = doc.as_object().unwrap();
    for key in [
        "format",
        "format_version",
        "version",
        "timestamp",
        "step",
        "time",
        "crate",
        "description",
        "n_particles",
        "positions",
        "velocities",
        "forces",
        "masses",
        "types",
        "box_matrix",
        "extra_scalars",
        "extra_vectors",
    ] {
        assert!(obj.contains_key(key), "missing member `{key}`");
    }
    assert_eq!(obj["format"], "oxiphysics-restart");
    assert_eq!(obj["format_version"], 2);
    assert_eq!(obj["n_particles"], 3);
    assert_eq!(obj["forces"].as_array().unwrap().len(), 3);
    assert_eq!(obj["box_matrix"].as_array().unwrap().len(), 3);
    assert_eq!(obj["extra_scalars"].as_array().unwrap().len(), 3);
    assert_eq!(obj["extra_vectors"][0]["name"], "spin");
    // Quotes, newlines, tabs and backslashes in the description are escaped,
    // so the document parses and the string comes back verbatim.
    assert_eq!(obj["description"], d.metadata.description.as_str());
}

#[test]
fn json_roundtrip_every_field() {
    let d = full_data();
    let r = RestartReader::read_json(&RestartWriter::write_json(&d).unwrap()).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn json_roundtrip_is_bit_exact_for_random_floats() {
    for seed in [11u64, 12, 13] {
        let d = random_float_data(60, seed);
        let r = RestartReader::read_json(&RestartWriter::write_json(&d).unwrap()).unwrap();
        assert_bit_exact(&d, &r);
    }
}

#[test]
fn restart_json_roundtrip_bit_exact_tricky_values() {
    let mut d = make_data(2);
    d.positions = vec![[1e-12, 1e300, -0.0], [f64::MIN_POSITIVE, 5e-324, f64::MAX]];
    d.masses = vec![0.1 + 0.2, 1.0 / 3.0];
    d.metadata.time = 2.2250738585072014e-308;
    let r = RestartReader::read_json(&RestartWriter::write_json(&d).unwrap()).unwrap();
    assert_bit_exact(&d, &r);
}

#[test]
fn json_roundtrip_nonfinite_values() {
    let mut d = make_data(1);
    d.positions[0] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    d.metadata.time = f64::NAN;
    let text = RestartWriter::write_json(&d).unwrap();
    serde_json::from_str::<serde_json::Value>(&text).unwrap();
    let r = RestartReader::read_json(&text).unwrap();
    assert!(r.positions[0][0].is_nan() && r.metadata.time.is_nan());
    assert_eq!(r.positions[0][1], f64::INFINITY);
    assert_eq!(r.positions[0][2], f64::NEG_INFINITY);
}

#[test]
fn json_reader_accepts_documents_from_the_six_decimal_writer() {
    // The pre-fix writer: no forces, no box, no extras, 6 decimals.
    let legacy = "{\n  \"version\": \"1.0\",\n  \"timestamp\": 5,\n  \"step\": 9,\n  \"time\": 0.250000,\n  \"crate\": \"c\",\n  \"description\": \"old\",\n  \"n_particles\": 2,\n  \"positions\": [[1.000000,2.000000,3.000000],[4.000000,5.000000,6.000000]],\n  \"velocities\": [[0.000000,0.000000,0.000000],[1.000000,1.000000,1.000000]],\n  \"masses\": [1.000000,2.000000],\n  \"types\": [0,1]\n}\n";
    let r = RestartReader::read_json(legacy).unwrap();
    assert_eq!(r.n_particles(), 2);
    assert_eq!(r.metadata.step, 9);
    assert_eq!(r.masses, vec![1.0, 2.0]);
    assert_eq!(r.types, vec![0, 1]);
    assert!(r.forces.is_empty());
    assert_eq!(r.box_matrix, [[0.0; 3]; 3]);
}

#[test]
fn json_reader_ignores_unknown_members_and_accepts_a_bom() {
    let text = "\u{feff}{\"unknown\": [1, 2, {\"x\": null}], \"positions\": [[1,2,3]], \"n_particles\": 1}";
    let r = RestartReader::read_json(text).unwrap();
    assert_eq!(r.positions, vec![[1.0, 2.0, 3.0]]);
}

#[test]
fn json_reader_rejects_malformed_documents() {
    let cases: [(&str, &str); 18] = [
        ("[]", "top-level"),
        ("42", "top-level"),
        ("{", "restart JSON"),
        ("{\"positions\": [[1,2]]}", "exactly 3"),
        ("{\"positions\": [[1,2,\"x\"]]}", "positions[0][2]"),
        ("{\"positions\": [1]}", "array of 3"),
        ("{\"positions\": 5}", "must be an array"),
        (
            "{\"positions\": [[1,2,3]], \"n_particles\": 2}",
            "n_particles",
        ),
        ("{\"positions\": [[1,2,3]], \"masses\": [1, 2]}", "masses"),
        ("{\"positions\": [[1,2,3]], \"types\": [-1]}", "types[0]"),
        (
            "{\"positions\": [[1,2,3]], \"types\": [4294967296]}",
            "types[0]",
        ),
        ("{\"step\": \"x\"}", "step"),
        ("{\"step\": 1.5}", "step"),
        ("{\"format\": \"other\"}", "format"),
        ("{\"format_version\": 3}", "format_version"),
        ("{\"box_matrix\": [[1,0,0]]}", "box_matrix"),
        ("{\"extra_scalars\": [{\"values\": []}]}", "name"),
        (
            "{\"extra_vectors\": [{\"name\": \"a\", \"values\": [[1,2]]}]}",
            "extra_vectors[0].values[0]",
        ),
    ];
    for (text, needle) in cases {
        let err = RestartReader::read_json(text).unwrap_err();
        assert!(
            err.contains(needle),
            "{text}: expected `{needle}` in `{err}`"
        );
    }
}

#[test]
fn json_reader_survives_every_single_byte_corruption() {
    // No byte pattern may panic; Ok and Err are both acceptable.
    let text = RestartWriter::write_json(&full_data()).unwrap();
    for i in 0..text.len() {
        for fill in [b'{', b'}', b'[', b']', b'"', b',', b':', b'9', b'-', b'\\'] {
            let mut bytes = text.clone().into_bytes();
            bytes[i] = fill;
            if let Ok(mutated) = String::from_utf8(bytes) {
                let _ = RestartReader::read_json(&mutated);
            }
        }
    }
}

#[test]
fn json_reader_rejects_deeply_nested_input_without_overflowing_the_stack() {
    let nested = format!("{{\"x\": {}{}}}", "[".repeat(100_000), "]".repeat(100_000));
    assert!(RestartReader::read_json(&nested).is_err());
    assert!(RestartReader::read_json(&"[".repeat(1_000_000)).is_err());
}

// ── detect_format / file-based routing (F2441) ────────────────────────────────

#[test]
fn test_detect_format_binary() {
    let d = make_data(1);
    let bytes = RestartWriter::write_binary(&d).unwrap();
    assert_eq!(RestartReader::detect_format(&bytes), RestartFormat::Binary);
}

#[test]
fn test_detect_format_ascii() {
    let text = "# OxiPhysics restart\nVERSION 1.0\n";
    assert_eq!(
        RestartReader::detect_format(text.as_bytes()),
        RestartFormat::Ascii
    );
}

#[test]
fn test_detect_format_json() {
    let j = "{\"n_particles\": 0}";
    assert_eq!(
        RestartReader::detect_format(j.as_bytes()),
        RestartFormat::Json
    );
}

#[test]
fn detect_format_tags_and_edge_cases() {
    let d = make_data(1);
    assert_eq!(
        RestartReader::detect_format(&RestartWriter::write_hdf5like(&d).unwrap()),
        RestartFormat::Hdf5Like
    );
    assert_eq!(
        RestartReader::detect_format(&RestartWriter::write_msgpack(&d).unwrap()),
        RestartFormat::MessagePack
    );
    // JSON with leading whitespace, a byte-order mark, or fewer than 4 bytes.
    for json in ["  \n\t{\"a\":1}", "\u{feff}{\"a\":1}", "{}", " {"] {
        assert_eq!(
            RestartReader::detect_format(json.as_bytes()),
            RestartFormat::Json,
            "{json:?}"
        );
    }
    assert_eq!(RestartReader::detect_format(b""), RestartFormat::Ascii);
    assert_eq!(RestartReader::detect_format(b"ab"), RestartFormat::Ascii);
    assert_eq!(
        RestartReader::detect_format(b"STEP 4"),
        RestartFormat::Ascii
    );
}

#[test]
fn restart_every_format_roundtrips_single_particle_test_through_a_file() {
    let d = RestartData::single_particle_test();
    for format in all_formats() {
        let r = file_roundtrip(format.clone(), &d, "single");
        assert_eq!(r, d, "format {format:?}");
    }
}

#[test]
fn restart_every_format_roundtrips_every_field_through_a_file() {
    let d = full_data();
    for format in all_formats() {
        let r = file_roundtrip(format.clone(), &d, "full");
        assert_bit_exact(&d, &r);
    }
}

#[test]
fn restart_json_file_is_read_back_by_the_file_reader() {
    // Before the fix `read` routed `{` to the binary decoder, so every JSON
    // restart file failed to load.
    let dir = TempDir::new("jsonfile");
    let path = dir.path("restart.json");
    let d = make_data(5);
    RestartWriter::new(&path, RestartFormat::Json)
        .write(&d)
        .unwrap();
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(RestartReader::detect_format(&bytes), RestartFormat::Json);
    assert_eq!(RestartReader::new(&path).read().unwrap(), d);
}

#[test]
fn reader_reports_missing_files_and_garbage() {
    let dir = TempDir::new("garbage");
    let missing = dir.path("nope.dat");
    let err = RestartReader::new(&missing).read().unwrap_err();
    assert!(err.contains("failed to read"), "{err}");

    let garbage = dir.path("garbage.dat");
    std::fs::write(&garbage, [0xFFu8, 0xFE, 0x00, 0x01, 0x02]).unwrap();
    // Not a known magic, not JSON: treated as text, which is not valid UTF-8.
    assert!(RestartReader::new(&garbage).read().is_err());

    let empty = dir.path("empty.dat");
    std::fs::write(&empty, b"").unwrap();
    assert!(RestartReader::new(&empty).read().is_err());
}

// ── RestartWriter::write ──────────────────────────────────────────────────────

#[test]
fn write_overwrites_existing_files_and_leaves_no_temporary() {
    let dir = TempDir::new("overwrite");
    let path = dir.path("restart.bin");
    let writer = RestartWriter::new(&path, RestartFormat::Binary);
    writer.write(&make_data(3)).unwrap();
    writer.write(&make_data(5)).unwrap();
    assert_eq!(RestartReader::new(&path).read().unwrap(), make_data(5));
    let leftovers: Vec<_> = std::fs::read_dir(&dir.0)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(leftovers, vec!["restart.bin".to_string()]);
}

#[test]
fn write_rejects_invalid_data_before_touching_the_file() {
    let dir = TempDir::new("invalid");
    let path = dir.path("restart.bin");
    let writer = RestartWriter::new(&path, RestartFormat::Binary);
    writer.write(&make_data(3)).unwrap();
    let before = std::fs::read(&path).unwrap();

    let mut bad = make_data(3);
    bad.masses.pop();
    let err = writer.write(&bad).unwrap_err();
    assert!(err.contains("masses"), "{err}");
    assert_eq!(
        std::fs::read(&path).unwrap(),
        before,
        "good file must survive"
    );
}

#[test]
fn write_reports_io_errors() {
    let dir = TempDir::new("ioerr");
    // Parent directory does not exist.
    let missing_dir = dir.path("no/such/dir/restart.bin");
    let err = RestartWriter::new(&missing_dir, RestartFormat::Binary)
        .write(&make_data(1))
        .unwrap_err();
    assert!(err.contains("failed to write"), "{err}");

    // The target is an existing directory: the rename fails and the temporary
    // file is cleaned up.
    let target = dir.path("target_is_a_dir");
    std::fs::create_dir(&target).unwrap();
    let err = RestartWriter::new(&target, RestartFormat::Binary)
        .write(&make_data(1))
        .unwrap_err();
    assert!(err.contains("failed to write"), "{err}");
    assert!(!std::path::Path::new(&format!("{target}.tmp")).exists());
}

// ── RestartValidator ──────────────────────────────────────────────────────────

#[test]
fn test_validator_sum_consistent() {
    let bytes = b"hello world";
    let sum = RestartValidator::checksum_sum(bytes);
    assert_eq!(sum, 1116);
    assert!(RestartValidator::verify_sum(bytes, sum));
}

#[test]
fn test_validator_xor_consistent() {
    let bytes = b"test data";
    let xor = RestartValidator::checksum_xor(bytes);
    assert!(RestartValidator::verify_xor(bytes, xor));
    assert_eq!(RestartValidator::checksum_xor(&[0b1010, 0b0110]), 0b1100);
}

#[test]
fn test_validator_sum_detects_corruption() {
    let bytes = b"original";
    let sum = RestartValidator::checksum_sum(bytes);
    let corrupt = b"0riginal";
    assert!(!RestartValidator::verify_sum(corrupt, sum));
}

#[test]
fn test_validator_xor_empty() {
    assert_eq!(RestartValidator::checksum_xor(b""), 0);
    assert_eq!(RestartValidator::checksum_sum(b""), 0);
}

#[test]
fn validator_detects_every_single_bit_flip_in_a_restart_blob() {
    let bytes = RestartWriter::write_binary(&make_data(3)).unwrap();
    let sum = RestartValidator::checksum_sum(&bytes);
    let xor = RestartValidator::checksum_xor(&bytes);
    assert!(RestartValidator::verify_sum(&bytes, sum));
    assert!(RestartValidator::verify_xor(&bytes, xor));
    for i in 0..bytes.len() {
        for bit in 0..8 {
            let mut flipped = bytes.clone();
            flipped[i] ^= 1 << bit;
            assert!(
                !RestartValidator::verify_sum(&flipped, sum),
                "sum checksum missed byte {i} bit {bit}"
            );
            assert!(
                !RestartValidator::verify_xor(&flipped, xor),
                "xor checksum missed byte {i} bit {bit}"
            );
        }
    }
}

// ── IncrementalRestart ────────────────────────────────────────────────────────

#[test]
fn test_incremental_mark_and_extract() {
    let full = make_data(5);
    let mut inc = IncrementalRestart::new();
    inc.mark_changed(1);
    inc.mark_changed(3);
    let delta = inc.extract_delta(&full);
    assert_eq!(delta.n_particles(), 2);
    assert_eq!(delta.positions, vec![full.positions[1], full.positions[3]]);
    assert_eq!(
        delta.velocities,
        vec![full.velocities[1], full.velocities[3]]
    );
    assert_eq!(delta.forces, vec![full.forces[1], full.forces[3]]);
    assert_eq!(delta.masses, vec![full.masses[1], full.masses[3]]);
    assert_eq!(delta.types, vec![1, 0]);
    assert_eq!(delta.metadata, full.metadata);
    assert_eq!(delta.box_matrix, full.box_matrix);
    delta.validate().unwrap();
}

#[test]
fn test_incremental_reset_clears_changes() {
    let mut inc = IncrementalRestart::new();
    inc.mark_changed(0);
    inc.mark_changed(2);
    inc.reset();
    assert!(inc.changed.is_empty());
    assert_eq!(inc.extract_delta(&make_data(3)).n_particles(), 0);
}

#[test]
fn test_incremental_no_duplicates() {
    let mut inc = IncrementalRestart::new();
    inc.mark_changed(0);
    inc.mark_changed(0);
    inc.mark_changed(0);
    assert_eq!(inc.changed.len(), 1);
}

#[test]
fn incremental_delta_follows_marking_order_and_drops_extras() {
    let mut full = make_data(4);
    full.extra_scalars = vec![("charge".to_string(), vec![1.0; 4])];
    full.extra_vectors = vec![("spin".to_string(), vec![[0.0, 0.0, 1.0]; 4])];
    let mut inc = IncrementalRestart::new();
    inc.mark_changed(3);
    inc.mark_changed(0);
    let delta = inc.extract_delta(&full);
    assert_eq!(delta.positions, vec![full.positions[3], full.positions[0]]);
    assert!(delta.extra_scalars.is_empty() && delta.extra_vectors.is_empty());
}

#[test]
fn incremental_delta_of_nothing_is_a_valid_empty_snapshot() {
    let delta = IncrementalRestart::new().extract_delta(&make_data(3));
    assert_eq!(delta.n_particles(), 0);
    delta.validate().unwrap();
    for format in all_formats() {
        let r = file_roundtrip(format.clone(), &delta, "empty_delta");
        assert_eq!(r, delta, "format {format:?}");
    }
}

// ── Conversion utilities ──────────────────────────────────────────────────────

#[test]
fn test_restart_to_xyz_header_line_count() {
    let d = make_data(4);
    let xyz = restart_to_xyz(&d);
    let lines: Vec<&str> = xyz.lines().collect();
    assert_eq!(
        lines,
        vec![
            "4",
            "Restart step=42 time=4.200000",
            "H 0.000000 0.000000 0.000000",
            "C 1.000000 0.500000 0.000000",
            "N 2.000000 1.000000 0.000000",
            "H 3.000000 1.500000 0.000000",
        ]
    );
}

#[test]
fn restart_to_xyz_uses_x_for_unknown_and_missing_types() {
    let mut d = make_data(3);
    d.types = vec![3, 7, 1];
    let xyz = restart_to_xyz(&d);
    let symbols: Vec<&str> = xyz
        .lines()
        .skip(2)
        .filter_map(|l| l.split_whitespace().next())
        .collect();
    assert_eq!(symbols, vec!["O", "X", "C"]);

    let bare = restart_to_xyz(&positions_only(2));
    let symbols: Vec<&str> = bare
        .lines()
        .skip(2)
        .filter_map(|l| l.split_whitespace().next())
        .collect();
    assert_eq!(symbols, vec!["X", "X"]);
}

#[test]
fn test_restart_to_lammps_dump_has_timestep() {
    let d = make_data(2);
    let dump = restart_to_lammps_dump(&d);
    let lines: Vec<&str> = dump.lines().collect();
    assert_eq!(
        lines,
        vec![
            "ITEM: TIMESTEP",
            "42",
            "ITEM: NUMBER OF ATOMS",
            "2",
            "ITEM: BOX BOUNDS pp pp pp",
            "0.0 10.000000",
            "0.0 10.000000",
            "0.0 10.000000",
            "ITEM: ATOMS id type x y z vx vy vz fx fy fz mass",
            "1 0 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 -9.810000 0.000000 1.000000",
            "2 1 1.000000 0.500000 0.000000 0.100000 0.000000 0.000000 0.000000 -9.810000 0.000000 1.000000",
        ]
    );
}

#[test]
fn test_restart_to_lammps_dump_atom_count() {
    let d = make_data(3);
    let dump = restart_to_lammps_dump(&d);
    let atom_lines: Vec<&str> = dump
        .lines()
        .skip_while(|l| !l.starts_with("ITEM: ATOMS"))
        .skip(1)
        .collect();
    assert_eq!(atom_lines.len(), 3);
}

#[test]
fn restart_to_lammps_dump_defaults_absent_arrays() {
    // A positions-only snapshot (e.g. a delta): zero velocity and force, unit
    // mass, type 0.
    let dump = restart_to_lammps_dump(&positions_only(1));
    let atom = dump.lines().last().unwrap();
    assert_eq!(
        atom,
        "1 0 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000"
    );
}

// ── CheckpointManager ─────────────────────────────────────────────────────────

/// Path `CheckpointManager` uses for `step` under `base`.
fn checkpoint_file(base: &str, step: u64) -> String {
    Path::new(base)
        .join(format!("checkpoint_{step:010}.bin"))
        .to_string_lossy()
        .into_owned()
}

/// Replace the checkpoint file at `path` with a non-empty directory so that
/// deleting it fails with something other than "not found".
fn block_removal(path: &str) {
    std::fs::remove_file(path).unwrap();
    std::fs::create_dir(path).unwrap();
    std::fs::write(Path::new(path).join("occupant"), b"x").unwrap();
}

/// Sorted names of the entries of `dir`.
fn dir_entries(dir: &str) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn test_checkpoint_manager_list_empty() {
    let dir = TempDir::new("ckpt_empty");
    let base = dir.path("ckpts");
    let mgr = CheckpointManager::new(&base, 3);
    assert!(mgr.list_checkpoints().is_empty());
    // Nothing is created before the first save.
    assert!(!Path::new(&base).exists());
    assert_eq!(mgr.load_latest().unwrap(), None);
}

#[test]
fn test_checkpoint_manager_prune_keeps_max() {
    let dir = TempDir::new("ckpt_prune");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 2);
    let d = make_data(1);
    let p1 = mgr.save_checkpoint(&d, 1).unwrap();
    let p2 = mgr.save_checkpoint(&d, 2).unwrap();
    let p3 = mgr.save_checkpoint(&d, 3).unwrap();
    assert_eq!(p1, checkpoint_file(&base, 1));
    assert_eq!(p3, checkpoint_file(&base, 3));
    assert_eq!(
        mgr.list_checkpoints(),
        vec![(2, p2.clone()), (3, p3.clone())]
    );
    // The evicted checkpoint is really gone from disk, the others are not.
    assert!(!Path::new(&p1).exists());
    assert!(Path::new(&p2).exists() && Path::new(&p3).exists());
}

#[test]
fn checkpoint_rolling_window_over_many_saves() {
    let dir = TempDir::new("ckpt_rolling");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 3);
    for step in 1..=6 {
        mgr.save_checkpoint(&make_data(step as usize), step)
            .unwrap();
    }
    let steps: Vec<u64> = mgr.list_checkpoints().iter().map(|(s, _)| *s).collect();
    assert_eq!(steps, vec![4, 5, 6]);
    assert_eq!(
        dir_entries(&base),
        vec![
            "checkpoint_0000000004.bin".to_string(),
            "checkpoint_0000000005.bin".to_string(),
            "checkpoint_0000000006.bin".to_string(),
        ]
    );
    assert_eq!(mgr.load_latest().unwrap(), Some(make_data(6)));
}

#[test]
fn checkpoint_manager_keeps_at_least_one_checkpoint() {
    assert_eq!(CheckpointManager::new("unused", 0).max_checkpoints, 1);

    // `max_checkpoints` is a public field; zeroing it afterwards must not make
    // a save delete the checkpoint it has just written.
    let dir = TempDir::new("ckpt_zero_max");
    let mut mgr = CheckpointManager::new(&dir.path("ckpts"), 2);
    mgr.max_checkpoints = 0;
    let first = mgr.save_checkpoint(&make_data(1), 1).unwrap();
    let second = mgr.save_checkpoint(&make_data(2), 2).unwrap();
    assert!(!Path::new(&first).exists());
    assert!(Path::new(&second).exists());
    assert_eq!(mgr.list_checkpoints(), vec![(2, second)]);
    assert_eq!(mgr.load_latest().unwrap(), Some(make_data(2)));
}

#[test]
fn checkpoint_load_latest_returns_the_most_recent_save() {
    let dir = TempDir::new("ckpt_latest");
    let mut mgr = CheckpointManager::new(&dir.path("ckpts"), 3);
    let first = make_data(2);
    let second = make_data(4);
    mgr.save_checkpoint(&first, 10).unwrap();
    mgr.save_checkpoint(&second, 20).unwrap();
    assert_eq!(mgr.load_latest().unwrap(), Some(second));
}

#[test]
fn checkpoint_saving_a_step_twice_replaces_the_file_and_keeps_one_entry() {
    let dir = TempDir::new("ckpt_resave");
    let mut mgr = CheckpointManager::new(&dir.path("ckpts"), 3);
    let p = mgr.save_checkpoint(&make_data(1), 5).unwrap();
    let second = make_data(2);
    let q = mgr.save_checkpoint(&second, 5).unwrap();
    assert_eq!(p, q);
    assert_eq!(mgr.list_checkpoints(), vec![(5, p)]);
    assert_eq!(mgr.load_latest().unwrap(), Some(second));
}

// F2444: save_checkpoint swallowed write errors, recorded the failed checkpoint
// and then pruned the older good ones.

#[test]
fn checkpoint_failed_save_keeps_older_checkpoints() {
    let dir = TempDir::new("ckpt_fail_keep");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 1);
    let good = make_data(3);
    let good_path = mgr.save_checkpoint(&good, 1).unwrap();

    // A directory squatting on step 2's file name makes the final rename fail.
    let blocked = checkpoint_file(&base, 2);
    std::fs::create_dir(&blocked).unwrap();
    let err = mgr.save_checkpoint(&make_data(3), 2).unwrap_err();
    assert!(err.contains("failed to write"), "{err}");

    // With max_checkpoints = 1 the old code would now have deleted `good_path`.
    assert_eq!(mgr.list_checkpoints(), vec![(1, good_path.clone())]);
    assert!(Path::new(&good_path).exists());
    assert_eq!(mgr.load_latest().unwrap(), Some(good));
    // The failed attempt left no temporary file behind: only the good
    // checkpoint and the blocking directory exist.
    assert_eq!(
        dir_entries(&base),
        vec![
            "checkpoint_0000000001.bin".to_string(),
            "checkpoint_0000000002.bin".to_string(),
        ]
    );
}

#[test]
fn checkpoint_base_dir_that_is_a_regular_file_is_an_error() {
    let dir = TempDir::new("ckpt_file_base");
    let not_a_dir = dir.path("not_a_dir");
    std::fs::write(&not_a_dir, b"occupied").unwrap();
    let mut mgr = CheckpointManager::new(&not_a_dir, 3);
    let before = mgr.list_checkpoints();
    let err = mgr.save_checkpoint(&make_data(1), 1).unwrap_err();
    assert!(
        err.contains("failed to create checkpoint directory"),
        "{err}"
    );
    assert_eq!(mgr.list_checkpoints(), before);
    assert!(mgr.list_checkpoints().is_empty());
    assert_eq!(mgr.load_latest().unwrap(), None);
    assert_eq!(std::fs::read(&not_a_dir).unwrap(), b"occupied");
}

#[test]
fn checkpoint_unusable_base_dir_after_a_good_save_keeps_the_good_one() {
    let dir = TempDir::new("ckpt_swap_base");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 1);
    let good = make_data(2);
    let good_path = mgr.save_checkpoint(&good, 1).unwrap();

    let not_a_dir = dir.path("not_a_dir");
    std::fs::write(&not_a_dir, b"x").unwrap();
    mgr.base_dir = not_a_dir;
    assert!(mgr.save_checkpoint(&make_data(2), 2).is_err());
    assert_eq!(mgr.list_checkpoints(), vec![(1, good_path.clone())]);
    assert!(Path::new(&good_path).exists());
    assert_eq!(mgr.load_latest().unwrap(), Some(good));
}

#[test]
fn checkpoint_inconsistent_data_is_rejected_without_side_effects() {
    let dir = TempDir::new("ckpt_invalid");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 1);
    let good_path = mgr.save_checkpoint(&make_data(2), 1).unwrap();

    let mut bad = make_data(2);
    bad.masses.pop();
    let err = mgr.save_checkpoint(&bad, 2).unwrap_err();
    assert!(err.contains("masses"), "{err}");
    assert_eq!(mgr.list_checkpoints(), vec![(1, good_path.clone())]);
    assert!(Path::new(&good_path).exists());
    assert!(!Path::new(&checkpoint_file(&base, 2)).exists());
}

#[test]
fn checkpoint_load_latest_reports_corrupt_and_missing_files() {
    let dir = TempDir::new("ckpt_load_err");
    let mut mgr = CheckpointManager::new(&dir.path("ckpts"), 3);
    let path = mgr.save_checkpoint(&make_data(3), 1).unwrap();

    // The old `load_latest` turned every failure into `None`.
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&path, &bytes[..bytes.len() / 2]).unwrap();
    assert!(mgr.load_latest().is_err());

    std::fs::remove_file(&path).unwrap();
    let err = mgr.load_latest().unwrap_err();
    assert!(err.contains("failed to read"), "{err}");
}

#[test]
fn checkpoint_prune_reports_undeletable_files_and_retries_later() {
    let dir = TempDir::new("ckpt_prune_err");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 2);
    let d = make_data(1);
    let p1 = mgr.save_checkpoint(&d, 1).unwrap();
    let p2 = mgr.save_checkpoint(&d, 2).unwrap();

    block_removal(&p1);
    mgr.max_checkpoints = 1;
    let err = mgr.prune_old_checkpoints().unwrap_err();
    assert!(err.contains("failed to remove checkpoint"), "{err}");
    // The stuck entry is still listed, so a later call retries it.
    assert_eq!(
        mgr.list_checkpoints(),
        vec![(1, p1.clone()), (2, p2.clone())]
    );

    // Once the obstacle is gone the file counts as deleted.
    std::fs::remove_dir_all(&p1).unwrap();
    mgr.prune_old_checkpoints().unwrap();
    assert_eq!(mgr.list_checkpoints(), vec![(2, p2)]);
}

#[test]
fn checkpoint_save_reports_prune_failure_but_keeps_the_new_checkpoint() {
    let dir = TempDir::new("ckpt_save_prune_err");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 1);
    let p1 = mgr.save_checkpoint(&make_data(1), 1).unwrap();
    block_removal(&p1);

    let newest = make_data(2);
    let err = mgr.save_checkpoint(&newest, 2).unwrap_err();
    assert!(
        err.contains("was written but pruning old checkpoints failed"),
        "{err}"
    );
    let p2 = checkpoint_file(&base, 2);
    assert_eq!(mgr.list_checkpoints(), vec![(1, p1), (2, p2)]);
    assert_eq!(mgr.load_latest().unwrap(), Some(newest));
}

#[test]
fn checkpoint_prune_treats_already_deleted_files_as_deleted() {
    let dir = TempDir::new("ckpt_prune_gone");
    let base = dir.path("ckpts");
    let mut mgr = CheckpointManager::new(&base, 2);
    let d = make_data(1);
    let p1 = mgr.save_checkpoint(&d, 1).unwrap();
    let p2 = mgr.save_checkpoint(&d, 2).unwrap();
    std::fs::remove_file(&p1).unwrap();
    mgr.max_checkpoints = 1;
    mgr.prune_old_checkpoints().unwrap();
    assert_eq!(mgr.list_checkpoints(), vec![(2, p2)]);
}
