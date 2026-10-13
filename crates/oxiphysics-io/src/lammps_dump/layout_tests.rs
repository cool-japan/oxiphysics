// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for reading the `ITEM: ATOMS` table by column name.

use super::{LammpsDumpFrame, LammpsDumpReader};
use std::io::Cursor;

/// An orthogonal dump frame with the box `x in [-5, 5]`, `y in [0, 20]`,
/// `z in [10, 40]`.
fn dump_text(columns: &str, rows: &[&str]) -> String {
    let mut text = format!(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{}\nITEM: BOX BOUNDS pp pp pp\n-5 5\n0 20\n10 40\nITEM: ATOMS {columns}\n",
        rows.len()
    );
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

fn read(text: &str) -> Result<LammpsDumpFrame, String> {
    LammpsDumpReader::read_frame(Cursor::new(text.as_bytes())).expect("a frame")
}

fn assert_close3(got: [f64; 3], want: [f64; 3]) {
    for k in 0..3 {
        assert!((got[k] - want[k]).abs() < 1e-12, "{got:?} vs {want:?}");
    }
}

#[test]
fn lammps_dump_scaled_coords_unscaled() {
    // `dump atom` writes scaled coordinates: x = xlo + xs * lx, ...
    let frame = read(&dump_text(
        "id type xs ys zs",
        &["1 1 0.5 0.5 0.5", "2 2 0.0 1.0 0.25"],
    ))
    .unwrap();
    assert_eq!(frame.atom_ids, vec![1, 2]);
    assert_eq!(frame.types, vec![1, 2]);
    assert_close3(frame.positions[0], [0.0, 10.0, 25.0]);
    assert_close3(frame.positions[1], [-5.0, 20.0, 17.5]);
    assert!(frame.custom_properties.is_empty());
    assert!(frame.velocities.is_none());
}

#[test]
fn scaled_unwrapped_coordinates_are_unscaled_too() {
    let frame = read(&dump_text("id type xsu ysu zsu", &["7 3 1.5 -0.5 2.0"])).unwrap();
    assert_close3(frame.positions[0], [10.0, -10.0, 70.0]);
}

#[test]
fn unwrapped_coordinates_are_used_as_they_are() {
    let frame = read(&dump_text("id type xu yu zu", &["1 1 -12.5 33.0 1.0"])).unwrap();
    assert_close3(frame.positions[0], [-12.5, 33.0, 1.0]);
}

#[test]
fn triclinic_scaled_coordinates_use_the_true_box() {
    // True box: xlo,ylo,zlo = 1,2,3; lx,ly,lz = 10,8,6; xy,xz,yz = 2,1,0.5.
    // LAMMPS writes the bounding box: x in [1, 11 + 3], y in [2, 10 + 0.5].
    let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS xy xz yz pp pp pp\n1.0 14.0 2.0\n2.0 10.5 1.0\n3.0 9.0 0.5\nITEM: ATOMS id type xs ys zs\n1 1 0.0 0.0 0.0\n2 1 1.0 1.0 1.0\n3 1 0.5 0.25 0.5\n";
    let frame = read(text).unwrap();
    assert_eq!(frame.tilt_factors, Some([2.0, 1.0, 0.5]));
    // The origin corner and the opposite corner of the cell.
    assert_close3(frame.positions[0], [1.0, 2.0, 3.0]);
    assert_close3(frame.positions[1], [14.0, 10.5, 9.0]);
    // x = 1 + 0.5*10 + 0.25*2 + 0.5*1, y = 2 + 0.25*8 + 0.5*0.5, z = 3 + 0.5*6.
    assert_close3(frame.positions[2], [7.0, 4.25, 6.0]);
}

#[test]
fn lammps_dump_mol_column_order() {
    // 'id mol type x y z' used to read mol as the type and type as x.
    let frame = read(&dump_text(
        "id mol type x y z",
        &["1 7 2 1.0 2.0 3.0", "2 7 1 4.0 5.0 6.0"],
    ))
    .unwrap();
    assert_eq!(frame.atom_ids, vec![1, 2]);
    assert_eq!(frame.types, vec![2, 1]);
    assert_close3(frame.positions[0], [1.0, 2.0, 3.0]);
    assert_close3(frame.positions[1], [4.0, 5.0, 6.0]);
    assert_eq!(frame.get_custom_property("mol"), Some(&[7.0, 7.0][..]));
}

#[test]
fn columns_may_appear_in_any_order() {
    let frame = read(&dump_text(
        "q z type y id x c_pe",
        &["-0.5 3.0 2 2.0 9 1.0 -4.25"],
    ))
    .unwrap();
    assert_eq!(frame.atom_ids, vec![9]);
    assert_eq!(frame.types, vec![2]);
    assert_close3(frame.positions[0], [1.0, 2.0, 3.0]);
    let names: Vec<&str> = frame
        .custom_properties
        .iter()
        .map(|(n, _)| n.as_str())
        .collect();
    assert_eq!(names, vec!["q", "c_pe"]);
    assert_eq!(frame.get_custom_property("q"), Some(&[-0.5][..]));
    assert_eq!(frame.get_custom_property("c_pe"), Some(&[-4.25][..]));
}

#[test]
fn velocities_are_found_by_name() {
    let frame = read(&dump_text(
        "id type vx vy vz x y z",
        &["1 1 0.1 0.2 0.3 1.0 2.0 3.0"],
    ))
    .unwrap();
    let v = frame.velocities.expect("velocities");
    assert_close3(v[0], [0.1, 0.2, 0.3]);
    assert_close3(frame.positions[0], [1.0, 2.0, 3.0]);
    assert!(frame.custom_properties.is_empty());
}

#[test]
fn incomplete_velocity_triple_stays_a_custom_column() {
    let frame = read(&dump_text("id type x y z vx", &["1 1 1.0 2.0 3.0 0.5"])).unwrap();
    assert!(frame.velocities.is_none());
    assert_eq!(frame.get_custom_property("vx"), Some(&[0.5][..]));
}

#[test]
fn plain_position_columns_win_over_other_families() {
    let frame = read(&dump_text(
        "id type x y z xs ys zs",
        &["1 1 1.0 2.0 3.0 0.5 0.5 0.5"],
    ))
    .unwrap();
    assert_close3(frame.positions[0], [1.0, 2.0, 3.0]);
    assert_eq!(frame.get_custom_property("xs"), Some(&[0.5][..]));
}

#[test]
fn text_element_column_is_skipped() {
    let frame = read(&dump_text(
        "id type element x y z",
        &["1 1 Si 1.0 2.0 3.0", "2 2 O 4.0 5.0 6.0"],
    ))
    .unwrap();
    assert_eq!(frame.types, vec![1, 2]);
    assert_close3(frame.positions[1], [4.0, 5.0, 6.0]);
    assert!(frame.custom_properties.is_empty());
}

#[test]
fn missing_required_columns_are_errors() {
    for (columns, missing) in [
        ("type x y z", "'id'"),
        ("id x y z", "'type'"),
        ("id type x y", "position"),
        ("id type xs ys", "position"),
        ("id type x ys zs", "position"),
        ("id type", "position"),
        ("", "'id'"),
    ] {
        let err = read(&dump_text(columns, &[])).unwrap_err();
        assert!(err.contains(missing), "{columns:?}: {err}");
    }
}

#[test]
fn rows_shorter_than_the_header_are_errors() {
    let err = read(&dump_text(
        "id type x y z c_pe",
        &["1 1 1.0 2.0 3.0 0.5", "2 1 1.0 2.0 3.0"],
    ))
    .unwrap_err();
    assert!(
        err.contains("atom line 1") && err.contains("5 columns"),
        "{err}"
    );
}

#[test]
fn unparsable_values_are_errors_and_name_the_column() {
    let err = read(&dump_text("id type x y z c_pe", &["1 1 1.0 2.0 3.0 oops"])).unwrap_err();
    assert!(err.contains("c_pe") && err.contains("oops"), "{err}");
    let err = read(&dump_text("id type x y z", &["1 1 1.0 two 3.0"])).unwrap_err();
    assert!(err.contains("two"), "{err}");
    let err = read(&dump_text("id type x y z", &["1 x 1.0 2.0 3.0"])).unwrap_err();
    assert!(err.contains("type"), "{err}");
    let err = read(&dump_text(
        "id type x y z vx vy vz",
        &["1 1 1 2 3 0 0 nope"],
    ))
    .unwrap_err();
    assert!(err.contains("vz"), "{err}");
}

#[test]
fn huge_atom_count_does_not_preallocate_or_hang() {
    let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n4000000000000\nITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\nITEM: ATOMS id type x y z\n1 1 0 0 0\n";
    let err = read(text).unwrap_err();
    assert!(err.contains("unexpected EOF"), "{err}");
}

#[test]
fn empty_frame_with_a_valid_header_is_fine() {
    let frame = read(&dump_text("id type x y z", &[])).unwrap();
    assert_eq!(frame.n_atoms, 0);
    assert!(frame.positions.is_empty() && frame.atom_ids.is_empty());
}

#[test]
fn writer_output_with_custom_columns_round_trips() {
    use super::LammpsDumpWriter;
    let frame = LammpsDumpFrame {
        timestep: 12,
        n_atoms: 2,
        box_bounds: [[0.0, 5.0]; 3],
        atom_ids: vec![4, 9],
        types: vec![1, 2],
        positions: vec![[1.0, 2.0, 3.0], [4.0, 0.5, 0.25]],
        velocities: Some(vec![[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]]),
        tilt_factors: None,
        custom_properties: vec![
            ("fx".to_string(), vec![1.0, -1.0]),
            ("mol".to_string(), vec![3.0, 3.0]),
        ],
    };
    let mut buf = Vec::new();
    LammpsDumpWriter::write_frame(&mut buf, &frame).unwrap();
    let back = read(&String::from_utf8(buf).unwrap()).unwrap();
    assert_eq!(back.atom_ids, frame.atom_ids);
    assert_eq!(back.types, frame.types);
    assert_eq!(back.positions, frame.positions);
    assert_eq!(back.velocities, frame.velocities);
    assert_eq!(back.custom_properties, frame.custom_properties);
}

#[test]
fn consecutive_frames_with_different_layouts_are_read_independently() {
    let first = dump_text("id type x y z", &["1 1 1.0 2.0 3.0"]);
    let second = dump_text("id mol type xs ys zs", &["1 4 2 0.5 0.5 0.5"]);
    let frames = LammpsDumpReader::read_all_frames(&format!("{first}{second}")).unwrap();
    assert_eq!(frames.len(), 2);
    assert_close3(frames[0].positions[0], [1.0, 2.0, 3.0]);
    assert_eq!(frames[1].types, vec![2]);
    assert_close3(frames[1].positions[0], [0.0, 10.0, 25.0]);
    assert_eq!(frames[1].get_custom_property("mol"), Some(&[4.0][..]));
}
