// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the frame writer and the per-atom accessors on frames whose arrays
//! disagree with `n_atoms`: they must report an error (or skip the missing
//! entries) instead of panicking.

use super::{
    LammpsDumpFrame, LammpsDumpReader, LammpsDumpWriter, extract_forces, extract_stress_tensor,
    mean_force_magnitude,
};
use std::io::Cursor;

/// A frame of three atoms with velocities and one custom property.
fn full_frame() -> LammpsDumpFrame {
    LammpsDumpFrame {
        timestep: 10,
        n_atoms: 3,
        box_bounds: [[0.0, 10.0]; 3],
        atom_ids: vec![1, 2, 3],
        types: vec![1, 2, 1],
        positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
        velocities: Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]),
        tilt_factors: None,
        custom_properties: vec![("energy".to_string(), vec![-1.5, 2.5, 3.5])],
    }
}

fn write(frame: &LammpsDumpFrame) -> (Result<(), String>, Vec<u8>) {
    let mut buf = Vec::new();
    let result = LammpsDumpWriter::write_frame(&mut buf, frame);
    (result, buf)
}

fn read(bytes: &[u8]) -> LammpsDumpFrame {
    LammpsDumpReader::read_frame(Cursor::new(bytes))
        .expect("the written text holds a frame")
        .expect("the written frame parses")
}

#[test]
fn a_consistent_frame_round_trips_through_the_writer() {
    let frame = full_frame();
    let (result, bytes) = write(&frame);
    result.unwrap();
    let parsed = read(&bytes);
    assert_eq!(parsed.atom_ids, frame.atom_ids);
    assert_eq!(parsed.types, frame.types);
    assert_eq!(parsed.positions, frame.positions);
    assert_eq!(parsed.velocities, frame.velocities);
    assert_eq!(
        parsed.get_custom_property("energy"),
        Some(&[-1.5, 2.5, 3.5][..])
    );
}

#[test]
fn a_short_array_is_an_error_naming_it_and_nothing_is_written() {
    type Shorten = fn(&mut LammpsDumpFrame);
    let cases: [(&str, Shorten); 5] = [
        ("atom_ids", |f| {
            f.atom_ids.pop();
        }),
        ("types", |f| {
            f.types.pop();
        }),
        ("positions", |f| {
            f.positions.pop();
        }),
        ("velocities", |f| {
            if let Some(v) = f.velocities.as_mut() {
                v.pop();
            }
        }),
        ("energy", |f| {
            if let Some((_, values)) = f.custom_properties.first_mut() {
                values.pop();
            }
        }),
    ];
    for (name, shorten) in cases {
        let mut frame = full_frame();
        shorten(&mut frame);
        let (result, bytes) = write(&frame);
        let message = result.expect_err(name);
        assert!(
            message.contains(name) && message.contains('3') && message.contains('2'),
            "{name}: {message}"
        );
        assert!(
            bytes.is_empty(),
            "{name}: the frame is rejected before any output"
        );
    }
}

#[test]
fn a_frame_that_declares_more_atoms_than_it_holds_is_an_error() {
    let mut frame = full_frame();
    frame.n_atoms = 4;
    let (result, bytes) = write(&frame);
    let message = result.unwrap_err();
    assert!(message.contains("4 atoms"), "{message}");
    assert!(bytes.is_empty());
}

#[test]
fn surplus_entries_beyond_n_atoms_are_not_written() {
    let mut frame = full_frame();
    frame.n_atoms = 2;
    let (result, bytes) = write(&frame);
    result.unwrap();
    let parsed = read(&bytes);
    assert_eq!(parsed.n_atoms, 2);
    assert_eq!(parsed.atom_ids, vec![1, 2]);
    assert_eq!(parsed.positions, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    assert_eq!(
        parsed.velocities,
        Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    );
    assert_eq!(parsed.get_custom_property("energy"), Some(&[-1.5, 2.5][..]));
}

#[test]
fn a_frame_without_atoms_is_written_and_read_back() {
    let frame = LammpsDumpFrame {
        timestep: 0,
        n_atoms: 0,
        box_bounds: [[0.0, 1.0]; 3],
        atom_ids: Vec::new(),
        types: Vec::new(),
        positions: Vec::new(),
        velocities: None,
        tilt_factors: None,
        custom_properties: Vec::new(),
    };
    let (result, bytes) = write(&frame);
    result.unwrap();
    assert_eq!(read(&bytes).n_atoms, 0);
}

/// Six stress components for `n` atoms, each with the given length.
fn stress_frame(n_atoms: usize, lengths: [usize; 6]) -> LammpsDumpFrame {
    let names = ["sxx", "syy", "szz", "sxy", "sxz", "syz"];
    LammpsDumpFrame {
        timestep: 0,
        n_atoms,
        box_bounds: [[0.0, 10.0]; 3],
        atom_ids: (1..=n_atoms as u32).collect(),
        types: vec![1; n_atoms],
        positions: vec![[0.0; 3]; n_atoms],
        velocities: None,
        tilt_factors: None,
        custom_properties: names
            .iter()
            .zip(lengths)
            .enumerate()
            .map(|(k, (name, len))| {
                let values = (0..len).map(|i| (10 * k + i) as f64).collect();
                ((*name).to_string(), values)
            })
            .collect(),
    }
}

#[test]
fn stress_components_are_returned_per_atom_in_voigt_order() {
    let tensors = extract_stress_tensor(&stress_frame(2, [2; 6])).unwrap();
    assert_eq!(
        tensors,
        vec![
            [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
            [1.0, 11.0, 21.0, 31.0, 41.0, 51.0]
        ]
    );
}

#[test]
fn a_stress_component_shorter_than_n_atoms_gives_none_instead_of_a_panic() {
    for short in 0..6 {
        let mut lengths = [3; 6];
        lengths[short] = 2;
        assert!(
            extract_stress_tensor(&stress_frame(3, lengths)).is_none(),
            "component {short} is missing an atom"
        );
    }
}

#[test]
fn a_frame_without_atoms_has_an_empty_stress_table() {
    assert_eq!(
        extract_stress_tensor(&stress_frame(0, [0; 6])),
        Some(Vec::new())
    );
}

/// Force components `fx fy fz` for `n_atoms` atoms, each with the given length.
fn force_frame(n_atoms: usize, lengths: [usize; 3]) -> LammpsDumpFrame {
    let names = ["fx", "fy", "fz"];
    LammpsDumpFrame {
        timestep: 0,
        n_atoms,
        box_bounds: [[0.0, 10.0]; 3],
        atom_ids: (1..=n_atoms as u32).collect(),
        types: vec![1; n_atoms],
        positions: vec![[0.0; 3]; n_atoms],
        velocities: None,
        tilt_factors: None,
        custom_properties: names
            .iter()
            .zip(lengths)
            .enumerate()
            .map(|(k, (name, len))| {
                let values = (0..len).map(|i| (10 * k + i) as f64).collect();
                ((*name).to_string(), values)
            })
            .collect(),
    }
}

#[test]
fn forces_are_returned_per_atom() {
    let frame = force_frame(2, [2; 3]);
    assert_eq!(
        extract_forces(&frame),
        Some(vec![[0.0, 10.0, 20.0], [1.0, 11.0, 21.0]])
    );
    // |(0, 10, 20)| = sqrt(500) and |(1, 11, 21)| = sqrt(563).
    let want = (500.0_f64.sqrt() + 563.0_f64.sqrt()) / 2.0;
    let got = mean_force_magnitude(&frame).unwrap();
    assert!((got - want).abs() < 1e-12, "{got} vs {want}");
}

#[test]
fn a_force_component_shorter_than_n_atoms_gives_none_instead_of_a_panic() {
    for short in 0..3 {
        let mut lengths = [3; 3];
        lengths[short] = 2;
        let frame = force_frame(3, lengths);
        assert!(
            extract_forces(&frame).is_none(),
            "component {short} is missing an atom"
        );
        assert!(mean_force_magnitude(&frame).is_none());
    }
}

#[test]
fn a_frame_without_atoms_has_an_empty_force_table() {
    let frame = force_frame(0, [0; 3]);
    assert_eq!(extract_forces(&frame), Some(Vec::new()));
    assert_eq!(mean_force_magnitude(&frame), Some(0.0));
}

#[test]
fn positions_by_type_skips_types_without_a_position() {
    let mut frame = full_frame();
    frame.types = vec![1, 2, 1, 1, 1];
    assert_eq!(
        frame.positions_by_type(1),
        vec![[1.0, 2.0, 3.0], [7.0, 8.0, 9.0]]
    );
    assert_eq!(frame.positions_by_type(2), vec![[4.0, 5.0, 6.0]]);
    assert!(frame.positions_by_type(9).is_empty());
}
