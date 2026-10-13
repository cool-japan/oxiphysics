// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for [`NetcdfTrajectoryBuilder::write_cdl`]: every frame must reach the
//! `data:` section and inconsistent frames must be rejected.

use crate::netcdf::{NetcdfReader, NetcdfTrajectoryBuilder, TrajectoryFrame};

fn frame(index: usize, n_atoms: usize) -> TrajectoryFrame {
    let base = 100.0 * index as f64;
    TrajectoryFrame {
        time_ps: 0.5 * index as f64,
        positions: (0..n_atoms)
            .map(|a| {
                let v = base + 10.0 * a as f64;
                [v, v + 0.25, v + 0.5]
            })
            .collect(),
        velocities: None,
        box_lengths: None,
    }
}

fn builder(n_frames: usize, n_atoms: usize) -> NetcdfTrajectoryBuilder {
    let mut b = NetcdfTrajectoryBuilder::new().with_title("three frames");
    for i in 0..n_frames {
        b.add_frame(frame(i, n_atoms));
    }
    b
}

fn cdl_of(b: &NetcdfTrajectoryBuilder) -> String {
    let mut buf = Vec::new();
    b.write_cdl(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

/// Values of the `name = ... ;` statement inside the `data:` section.
fn data_values(cdl: &str, name: &str) -> Option<Vec<f64>> {
    let data = cdl.split("data:\n").nth(1)?;
    let prefix = format!("\t{name} = ");
    let line = data.lines().find(|l| l.starts_with(&prefix))?;
    let body = line.strip_prefix(&prefix)?.trim_end_matches(';').trim();
    Some(body.split(',').map(|v| v.trim().parse().unwrap()).collect())
}

#[test]
fn every_frame_reaches_the_data_section() {
    let b = builder(3, 2);
    let cdl = cdl_of(&b);
    assert!(cdl.contains("data:"), "no data section:\n{cdl}");
    let coords = data_values(&cdl, "coordinates").unwrap();
    assert_eq!(coords.len(), 18, "3 frames x 2 atoms x 3 components");
    let expected: Vec<f64> = b
        .frames
        .iter()
        .flat_map(|f| f.positions.iter().flatten().copied())
        .collect();
    assert_eq!(coords, expected);
    assert_eq!(data_values(&cdl, "time").unwrap(), vec![0.0, 0.5, 1.0]);
}

#[test]
fn declarations_match_the_written_data() {
    let cdl = cdl_of(&builder(3, 2));
    assert!(
        cdl.contains("\tframe = UNLIMITED ; // (3 currently)"),
        "{cdl}"
    );
    assert!(cdl.contains("\tatom = 2 ;"));
    assert!(cdl.contains("\tspatial = 3 ;"));
    assert!(
        !cdl.contains("label"),
        "unused label dimension must be gone"
    );
    assert!(cdl.contains("double coordinates(frame, atom, spatial)"));
    assert!(cdl.contains("double time(frame)"));
    assert!(!cdl.contains("cell_lengths"), "no frame has a box:\n{cdl}");
    assert!(
        !cdl.contains("velocities"),
        "no frame has velocities:\n{cdl}"
    );
}

#[test]
fn round_trips_through_the_cdl_reader() {
    let b = builder(3, 2);
    let reader = NetcdfReader::from_cdl(&cdl_of(&b)).unwrap();
    assert_eq!(reader.unlimited_dim.as_deref(), Some("frame"));
    assert_eq!(reader.get_dimension("atom"), Some(2));
    assert_eq!(reader.get_dimension("spatial"), Some(3));
    let expected: Vec<f64> = b
        .frames
        .iter()
        .flat_map(|f| f.positions.iter().flatten().copied())
        .collect();
    assert_eq!(reader.get_data("coordinates").unwrap(), expected.as_slice());
    assert_eq!(reader.get_data("time").unwrap(), &[0.0, 0.5, 1.0]);
    let coords = reader.get_variable("coordinates").unwrap();
    assert_eq!(coords.dimensions, vec!["frame", "atom", "spatial"]);
    assert_eq!(coords.units, "angstrom");
}

#[test]
fn cell_lengths_are_written_only_when_every_frame_has_a_box() {
    let mut b = builder(3, 1);
    for (i, f) in b.frames.iter_mut().enumerate() {
        f.box_lengths = Some([10.0 + i as f64, 20.0, 30.0]);
    }
    let cdl = cdl_of(&b);
    assert!(cdl.contains("double cell_lengths(frame, spatial)"));
    assert_eq!(
        data_values(&cdl, "cell_lengths").unwrap(),
        vec![10.0, 20.0, 30.0, 11.0, 20.0, 30.0, 12.0, 20.0, 30.0]
    );
    // One frame without a box: the variable disappears entirely.
    b.frames[1].box_lengths = None;
    let partial = cdl_of(&b);
    assert!(!partial.contains("cell_lengths"), "{partial}");
}

#[test]
fn velocities_are_written_only_when_every_frame_has_them() {
    let mut b = builder(2, 2);
    for f in b.frames.iter_mut() {
        f.velocities = Some(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    }
    let cdl = cdl_of(&b);
    assert!(cdl.contains("double velocities(frame, atom, spatial)"));
    assert!(cdl.contains("velocities:units = \"nanometer/picosecond\""));
    assert_eq!(
        data_values(&cdl, "velocities").unwrap(),
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    );
    b.frames[0].velocities = None;
    assert!(!cdl_of(&b).contains("velocities"));
}

#[test]
fn nanometre_units_label_positions_and_box() {
    let mut b = builder(1, 1).in_nanometres();
    b.frames[0].box_lengths = Some([1.0, 2.0, 3.0]);
    let cdl = cdl_of(&b);
    assert!(cdl.contains("coordinates:units = \"nanometer\""));
    assert!(cdl.contains("cell_lengths:units = \"nanometer\""));
    assert!(!cdl.contains("angstrom"));
}

#[test]
fn ragged_frames_are_rejected() {
    let mut b = builder(2, 3);
    b.add_frame(frame(2, 2));
    let mut buf = Vec::new();
    let err = b.write_cdl(&mut buf).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("frame 2"), "{err}");
}

#[test]
fn mismatched_velocity_count_is_rejected() {
    let mut b = builder(2, 3);
    b.frames[1].velocities = Some(vec![[0.0; 3]; 2]);
    let mut buf = Vec::new();
    let err = b.write_cdl(&mut buf).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn empty_builder_writes_a_header_without_data() {
    let cdl = cdl_of(&NetcdfTrajectoryBuilder::new().with_title("empty"));
    assert!(cdl.contains("(0 currently)"));
    assert!(!cdl.contains("data:"));
    let reader = NetcdfReader::from_cdl(&cdl).unwrap();
    assert!(reader.get_data("coordinates").is_none_or(|d| d.is_empty()));
}

#[test]
fn frames_without_atoms_do_not_emit_empty_statements() {
    let cdl = cdl_of(&builder(2, 0));
    assert!(data_values(&cdl, "coordinates").is_none(), "{cdl}");
    assert_eq!(data_values(&cdl, "time").unwrap(), vec![0.0, 0.5]);
    assert!(NetcdfReader::from_cdl(&cdl).is_ok());
}

#[test]
fn title_and_application_are_escaped() {
    let b = NetcdfTrajectoryBuilder::new()
        .with_title("a \"quoted\" back\\slash\nsecond line")
        .with_application("tab\there");
    let cdl = cdl_of(&b);
    assert!(
        cdl.contains("\t:title = \"a \\\"quoted\\\" back\\\\slash\\nsecond line\" ;"),
        "{cdl}"
    );
    assert!(cdl.contains("\t:application = \"tab\\there\" ;"), "{cdl}");
    assert_eq!(cdl.lines().filter(|l| l.contains("second line")).count(), 1);
}

#[test]
fn values_round_trip_exactly() {
    let mut b = NetcdfTrajectoryBuilder::new();
    b.add_frame(TrajectoryFrame {
        time_ps: 0.1,
        positions: vec![[1.0 / 3.0, -2.5e-7, 1.0e10]],
        velocities: None,
        box_lengths: None,
    });
    let reader = NetcdfReader::from_cdl(&cdl_of(&b)).unwrap();
    assert_eq!(
        reader.get_data("coordinates").unwrap(),
        &[1.0 / 3.0, -2.5e-7, 1.0e10]
    );
    assert_eq!(reader.get_data("time").unwrap(), &[0.1]);
}
