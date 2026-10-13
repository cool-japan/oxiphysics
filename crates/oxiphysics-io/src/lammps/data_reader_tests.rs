// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for [`LammpsDataReader`]: atom styles, section handling and the
//! consistency checks that replace the old silent fallbacks.

use super::data_reader::LammpsDataReader;
use super::types::{LammpsAtomStyle, LammpsDataSectionBuilder, LammpsDataWriter};
use crate::Error;

const POSITIONS: [[f64; 3]; 3] = [[1.0, 2.0, 3.0], [-4.5, 0.25, 6.0], [7.0, -8.0, 9.5]];

fn assert_close3(got: [f64; 3], want: [f64; 3]) {
    for k in 0..3 {
        assert!((got[k] - want[k]).abs() < 1e-9, "{got:?} vs {want:?}");
    }
}

fn header(n_atoms: usize) -> String {
    let mut s = LammpsDataWriter::write_header(n_atoms, 0, [-5.0, -6.0, -7.0], [5.0, 6.0, 7.0]);
    s.push_str("2 atom types\n");
    s
}

fn masses_section() -> String {
    let mut b = LammpsDataSectionBuilder::new();
    b.add_mass(1, 12.011);
    b.add_mass(2, 15.999);
    format!("\n{}", b.build())
}

fn error_text(result: crate::Result<LammpsDataReader>) -> String {
    match result {
        Ok(_) => panic!("expected a parse error"),
        Err(Error::Parse(message)) => message,
        Err(other) => panic!("expected Error::Parse, got {other:?}"),
    }
}

#[test]
fn atomic_style_round_trips() {
    let text = LammpsDataWriter::write_complete(
        &POSITIONS,
        &[12.0, 16.0, 12.0],
        &[1, 2, 1],
        [-5.0, -6.0, -7.0],
        [5.0, 6.0, 7.0],
    );
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Atomic));
    assert_eq!(r.types(), &[1, 2, 1]);
    assert_eq!(r.atom_ids(), &[1, 2, 3]);
    assert_eq!(r.masses(), &[12.0, 16.0, 12.0]);
    for (got, want) in r.positions().iter().zip(POSITIONS) {
        assert_close3(*got, want);
    }
    assert!(r.charges().is_empty() && r.molecule_ids().is_empty());
    assert_eq!(r.box_bounds(), ([-5.0, -6.0, -7.0], [5.0, 6.0, 7.0]));
}

#[test]
fn charge_style_round_trips() {
    let charges = [0.5, -0.25, 0.0];
    let types = [1u32, 2, 2];
    let mut text = header(3);
    text.push_str(&masses_section());
    text.push_str(&LammpsDataWriter::write_atoms_charge(
        &POSITIONS,
        &[12.011, 15.999, 15.999],
        &charges,
        &types,
    ));
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Charge));
    assert_eq!(r.types(), &types);
    assert_eq!(r.charges(), &charges);
    assert!(r.molecule_ids().is_empty());
    for (got, want) in r.positions().iter().zip(POSITIONS) {
        assert_close3(*got, want);
    }
    assert!((r.masses()[0] - 12.011).abs() < 1e-9);
    assert!((r.masses()[1] - 15.999).abs() < 1e-9);
}

#[test]
fn full_style_round_trips() {
    let charges = [-0.82, 0.41, 0.41];
    let types = [1u32, 2, 2];
    let mols = [7u32, 7, 8];
    let mut text = header(3);
    text.push_str(&masses_section());
    text.push_str(&LammpsDataWriter::write_atoms_full(
        &POSITIONS, &charges, &types, &mols,
    ));
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Full));
    assert_eq!(r.types(), &types);
    assert_eq!(r.charges(), &charges);
    assert_eq!(r.molecule_ids(), &[7, 7, 8]);
    for (got, want) in r.positions().iter().zip(POSITIONS) {
        assert_close3(*got, want);
    }
}

#[test]
fn molecular_style_round_trips() {
    let types = [2u32, 1, 2];
    let mols = [1u32, 1, 2];
    let mut text = header(3);
    text.push_str(&masses_section());
    text.push_str(&LammpsDataWriter::write_atoms_molecular(
        &POSITIONS, &types, &mols,
    ));
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Molecular));
    assert_eq!(r.types(), &types);
    assert_eq!(r.molecule_ids(), &[1, 1, 2]);
    assert!(r.charges().is_empty());
    assert_close3(r.positions()[1], POSITIONS[1]);
    assert!((r.masses()[0] - 15.999).abs() < 1e-9);
}

#[test]
fn masses_followed_by_coefficient_sections_parse() {
    // The builder emits Masses, Pair Coeffs, Bond Coeffs and Angle Coeffs; the old
    // reader treated "Pair" as a mass row and failed the whole file.
    let mut b = LammpsDataSectionBuilder::new();
    b.add_mass(1, 12.011);
    b.add_mass(2, 15.999);
    b.add_pair_coeff(1, 1, 0.086, 3.4);
    b.add_pair_coeff(2, 2, 0.2, 3.1);
    b.add_bond_coeff(1, 350.0, 1.5);
    b.add_angle_coeff(1, 50.0, 109.5);
    let mut text = header(3);
    text.push('\n');
    text.push_str(&b.build());
    text.push_str(&LammpsDataWriter::write_atoms_atomic(
        &POSITIONS,
        &[12.011, 15.999, 12.011],
        &[1, 2, 1],
    ));
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.types(), &[1, 2, 1]);
    assert_eq!(r.masses().len(), 3);
    assert!((r.masses()[0] - 12.011).abs() < 1e-9);
    assert!((r.masses()[1] - 15.999).abs() < 1e-9);
    assert!((r.masses()[2] - 12.011).abs() < 1e-9);
    assert_close3(r.positions()[2], POSITIONS[2]);
}

#[test]
fn coefficient_rows_are_never_read_as_atoms() {
    // Coefficient and topology sections *after* Atoms contain rows that look like
    // atom rows; they must not be appended to the atom list.
    let mut text = header(2);
    text.push_str(&masses_section());
    text.push_str(&LammpsDataWriter::write_atoms_atomic(
        &POSITIONS[..2],
        &[12.0, 16.0],
        &[1, 2],
    ));
    text.push_str(
        "\nAngles\n\n1 1 1 2 1\n2 1 2 1 2\n\nDihedrals\n\n1 1 1 2 3 4\n\nImpropers\n\n1 1 1 2 3 4\n\
         \nPair Coeffs # lj/cut\n\n1 0.1 3.4\n2 0.2 3.1\n\nCS-Info\n\n1 1 1 1 1 2.5\n",
    );
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.positions().len(), 2);
    assert_eq!(r.types(), &[1, 2]);
    assert!(r.bonds().is_empty());
}

#[test]
fn write_data_style_file_with_velocities_and_bonds() {
    let text = "LAMMPS data file via write_data, version 29 Sep 2021, timestep = 0\n\
\n\
3 atoms\n\
2 atom types\n\
2 bonds\n\
1 bond types\n\
\n\
-5.0 5.0 xlo xhi\n\
-5.5 5.5 ylo yhi\n\
-6.0 6.0 zlo zhi\n\
\n\
Masses\n\
\n\
1 15.9994 # O\n\
2 1.008 # H\n\
\n\
Pair Coeffs # lj/cut/coul/long\n\
\n\
1 0.1553 3.166\n\
2 0.0 0.0\n\
\n\
Bond Coeffs # harmonic\n\
\n\
1 450.0 1.0\n\
\n\
Atoms # full\n\
\n\
1 1 1 -0.82 0.0 0.0 0.0 0 0 0\n\
3 1 2 0.41 -0.757 0.586 0.0 0 0 0\n\
2 1 2 0.41 0.757 0.586 0.0 1 0 -1\n\
\n\
Velocities\n\
\n\
2 0.4 0.5 0.6\n\
1 0.1 0.2 0.3\n\
3 -0.1 -0.2 -0.3\n\
\n\
Bonds\n\
\n\
1 1 1 2\n\
2 1 1 3\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Full));
    assert_eq!(r.atom_ids(), &[1, 3, 2]);
    assert_eq!(r.types(), &[1, 2, 2]);
    assert_eq!(r.molecule_ids(), &[1, 1, 1]);
    assert_eq!(r.charges(), &[-0.82, 0.41, 0.41]);
    assert_eq!(r.masses(), &[15.9994, 1.008, 1.008]);
    assert_close3(r.positions()[1], [-0.757, 0.586, 0.0]);
    assert_close3(r.positions()[2], [0.757, 0.586, 0.0]);
    assert_eq!(r.box_bounds(), ([-5.0, -5.5, -6.0], [5.0, 5.5, 6.0]));
    assert_eq!(r.tilt_factors(), None);
    // Velocities are matched by atom id, not by file order.
    let v = r.velocities().unwrap();
    assert_close3(v[0], [0.1, 0.2, 0.3]);
    assert_close3(v[1], [-0.1, -0.2, -0.3]);
    assert_close3(v[2], [0.4, 0.5, 0.6]);
    assert_eq!(r.bonds(), &[(1, 1, 2), (1, 1, 3)]);
}

#[test]
fn bonds_written_by_the_writer_are_read_back() {
    let mut text = LammpsDataWriter::write_header(3, 2, [0.0; 3], [10.0; 3]);
    text.push_str(&LammpsDataWriter::write_atoms_atomic(
        &POSITIONS,
        &[1.0; 3],
        &[1, 1, 1],
    ));
    text.push_str(&LammpsDataWriter::write_bonds(&[(1, 1, 2), (2, 2, 3)]));
    let r = LammpsDataReader::from_data_str(&text).unwrap();
    assert_eq!(r.bonds(), &[(1, 1, 2), (2, 2, 3)]);
}

#[test]
fn box_comes_from_the_file() {
    let text = "title\n\n1 atoms\n\n-1.5 2.5 xlo xhi\n10 30 ylo yhi\n0.125 0.25 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.box_bounds(), ([-1.5, 10.0, 0.125], [2.5, 30.0, 0.25]));
}

#[test]
fn triclinic_tilt_is_reported() {
    let text = "title\n\n1 atoms\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n1.5 -2.0 0.25 xy xz yz\n\nAtoms # atomic\n\n1 1 0 0 0\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.tilt_factors(), Some([1.5, -2.0, 0.25]));
    assert_eq!(r.box_bounds(), ([0.0; 3], [10.0; 3]));
}

#[test]
fn missing_box_lines_are_errors_not_defaults() {
    let atoms = "\nAtoms # atomic\n\n1 1 0 0 0\n";
    let all = ["0 10 xlo xhi\n", "0 10 ylo yhi\n", "0 10 zlo zhi\n"];
    for missing in 0..3 {
        let mut text = String::from("title\n\n1 atoms\n\n");
        for (k, line) in all.iter().enumerate() {
            if k != missing {
                text.push_str(line);
            }
        }
        text.push_str(atoms);
        let message = error_text(LammpsDataReader::from_data_str(&text));
        let axis = ["x", "y", "z"][missing];
        assert!(
            message.contains(&format!("{axis}lo {axis}hi")),
            "missing {axis}: {message}"
        );
    }
    // No header at all.
    let message = error_text(LammpsDataReader::from_data_str(
        "title\n\nAtoms # atomic\n\n1 1 0 0 0\n",
    ));
    assert!(message.contains("xlo xhi"), "{message}");
}

#[test]
fn invalid_box_bounds_are_errors() {
    for bad in ["5 5 xlo xhi", "6 5 xlo xhi", "nan 5 xlo xhi", "a 5 xlo xhi"] {
        let text = format!("title\n\n{bad}\n0 1 ylo yhi\n0 1 zlo zhi\n");
        let message = error_text(LammpsDataReader::from_data_str(&text));
        assert!(message.contains("line 3"), "{bad}: {message}");
    }
}

#[test]
fn title_line_is_never_interpreted() {
    let text = "Atoms of water\n\n1 atoms\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\nAtoms # atomic\n\n1 1 0.5 0.5 0.5\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.positions(), &[[0.5, 0.5, 0.5]]);
    // A title that looks like a header line must not set the atom count either.
    let text = "7 atoms\n\n1 atoms\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\nAtoms # atomic\n\n1 1 0.5 0.5 0.5\n";
    assert_eq!(LammpsDataReader::from_data_str(text).unwrap().types(), &[1]);
}

#[test]
fn atom_count_must_match_the_header() {
    let text = "title\n\n3 atoms\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n2 1 1 1 1\n";
    let message = error_text(LammpsDataReader::from_data_str(text));
    assert!(
        message.contains("3 atoms") && message.contains("holds 2"),
        "{message}"
    );
}

#[test]
fn bond_count_must_match_the_header() {
    let text = "title\n\n1 atoms\n2 bonds\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n\nBonds\n\n1 1 1 1\n";
    let message = error_text(LammpsDataReader::from_data_str(text));
    assert!(
        message.contains("2 bonds") && message.contains("holds 1"),
        "{message}"
    );
}

#[test]
fn atom_style_hint_selects_the_columns() {
    let row = "1 1 2 0.5 1.0 2.0 3.0\n";
    for (hint, style, ok) in [
        ("full", LammpsAtomStyle::Full, true),
        ("bond", LammpsAtomStyle::Molecular, false),
        ("angle", LammpsAtomStyle::Molecular, false),
    ] {
        let text = format!(
            "t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # {hint}\n\n{row}"
        );
        let result = LammpsDataReader::from_data_str(&text);
        if ok {
            assert_eq!(result.unwrap().atom_style(), Some(style));
        } else {
            // 'id mol type x y z' reads 0.5 as x, so the row is valid but the
            // positions differ from the full-style reading.
            let r = result.unwrap();
            assert_eq!(r.atom_style(), Some(style));
            assert_eq!(r.positions(), &[[0.5, 1.0, 2.0]]);
        }
    }
}

#[test]
fn no_style_comment_recognises_atomic_layouts_only() {
    let make = |row: &str| {
        format!("t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms\n\n{row}\n")
    };
    let r = LammpsDataReader::from_data_str(&make("1 2 1.0 2.0 3.0")).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Atomic));
    assert_eq!(r.types(), &[2]);
    let r = LammpsDataReader::from_data_str(&make("1 2 1.0 2.0 3.0 0 0 -1")).unwrap();
    assert_eq!(r.positions(), &[[1.0, 2.0, 3.0]]);
    for ambiguous in ["1 2 0.5 1.0 2.0 3.0", "1 1 2 0.5 1.0 2.0 3.0"] {
        let message = error_text(LammpsDataReader::from_data_str(&make(ambiguous)));
        assert!(message.contains("infer the atom style"), "{message}");
    }
    // The caller can resolve the ambiguity.
    let r = LammpsDataReader::from_data_str_with_style(
        &make("1 2 0.5 1.0 2.0 3.0"),
        LammpsAtomStyle::Charge,
    )
    .unwrap();
    assert_eq!(r.charges(), &[0.5]);
    assert_eq!(r.positions(), &[[1.0, 2.0, 3.0]]);
}

#[test]
fn explicit_style_overrides_the_comment() {
    let text = "t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # atomic\n\n1 3 0.5 1.0 2.0 3.0\n";
    let r = LammpsDataReader::from_data_str_with_style(text, LammpsAtomStyle::Charge).unwrap();
    assert_eq!(r.atom_style(), Some(LammpsAtomStyle::Charge));
    assert_eq!(r.types(), &[3]);
    assert_eq!(r.charges(), &[0.5]);
}

#[test]
fn unsupported_styles_are_rejected() {
    for style in ["sphere", "dipole", "hybrid sphere dipole", "ellipsoid"] {
        let text = format!(
            "t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # {style}\n\n1 1 1.0 1.0 0.5 0.5 0.5\n"
        );
        let message = error_text(LammpsDataReader::from_data_str(&text));
        assert!(message.contains("not supported"), "{style}: {message}");
        assert!(message.contains("line 9"), "{style}: {message}");
    }
}

#[test]
fn short_atom_rows_are_errors() {
    let text = "t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # full\n\n1 1 1 0.5 1.0 2.0\n";
    let message = error_text(LammpsDataReader::from_data_str(text));
    assert!(message.contains("at least 7 columns"), "{message}");
    assert!(message.contains("line 11"), "{message}");
}

#[test]
fn bad_numbers_and_types_are_errors() {
    let make = |row: &str| {
        format!(
            "t\n\n1 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # atomic\n\n{row}\n"
        )
    };
    for row in ["1 1 x 0 0", "1 0 0 0 0", "1 -1 0 0 0", "q 1 0 0 0"] {
        assert!(
            LammpsDataReader::from_data_str(&make(row)).is_err(),
            "row {row:?} must be rejected"
        );
    }
}

#[test]
fn atom_types_must_not_exceed_the_declared_count() {
    let text = "t\n\n1 atoms\n2 atom types\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # atomic\n\n1 3 0 0 0\n";
    let message = error_text(LammpsDataReader::from_data_str(text));
    assert!(message.contains("atom type 3"), "{message}");
}

#[test]
fn masses_must_be_positive_and_complete() {
    let make = |masses: &str| {
        format!(
            "t\n\n1 atoms\n2 atom types\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nMasses\n\n{masses}\nAtoms # atomic\n\n1 1 0 0 0\n"
        )
    };
    assert!(LammpsDataReader::from_data_str(&make("1 12.0\n2 16.0\n")).is_ok());
    assert!(LammpsDataReader::from_data_str(&make("1 0.0\n2 16.0\n")).is_err());
    assert!(LammpsDataReader::from_data_str(&make("1 -3.0\n2 16.0\n")).is_err());
    assert!(LammpsDataReader::from_data_str(&make("0 3.0\n2 16.0\n")).is_err());
    assert!(LammpsDataReader::from_data_str(&make("1 abc\n2 16.0\n")).is_err());
    assert!(LammpsDataReader::from_data_str(&make("1\n2 16.0\n")).is_err());
    let message = error_text(LammpsDataReader::from_data_str(&make("1 12.0\n")));
    assert!(
        message.contains("2 atom types") && message.contains("1 entries"),
        "{message}"
    );
}

#[test]
fn types_without_a_masses_entry_get_the_documented_placeholder() {
    let text = "t\n\n2 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nMasses\n\n1 4.0\n\nAtoms # atomic\n\n1 1 0 0 0\n2 2 1 1 1\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.masses(), &[4.0, 1.0]);
}

#[test]
fn velocities_must_cover_every_atom_exactly_once() {
    let make = |velocities: &str| {
        format!(
            "t\n\n2 atoms\n\n0 9 xlo xhi\n0 9 ylo yhi\n0 9 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n2 1 1 1 1\n\nVelocities\n\n{velocities}"
        )
    };
    assert!(LammpsDataReader::from_data_str(&make("1 0 0 0\n2 1 1 1\n")).is_ok());
    for bad in [
        "1 0 0 0\n",
        "1 0 0 0\n1 1 1 1\n",
        "1 0 0 0\n9 1 1 1\n",
        "1 0 0\n2 1 1 1\n",
    ] {
        assert!(
            LammpsDataReader::from_data_str(&make(bad)).is_err(),
            "velocities {bad:?} must be rejected"
        );
    }
}

#[test]
fn velocities_are_absent_without_the_section() {
    let text =
        LammpsDataWriter::write_complete(&POSITIONS, &[1.0; 3], &[1, 1, 1], [0.0; 3], [9.0; 3]);
    assert!(
        LammpsDataReader::from_data_str(&text)
            .unwrap()
            .velocities()
            .is_none()
    );
}

#[test]
fn inline_comments_blank_runs_and_crlf_are_handled() {
    let text = "title # ignored\r\n\r\n\r\n2 atoms # two\r\n1 atom types\r\n\r\n0 9 xlo xhi # x\r\n0 9 ylo yhi\r\n0 9 zlo zhi\r\n\r\n\r\nMasses # label\r\n\r\n1 2.5 # C\r\n\r\nAtoms # atomic # extra words\r\n\r\n1 1 0.5 0.5 0.5 # first\r\n\r\n2 1 1.5 1.5 1.5\r\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert_eq!(r.positions(), &[[0.5, 0.5, 0.5], [1.5, 1.5, 1.5]]);
    assert_eq!(r.masses(), &[2.5, 2.5]);
}

#[test]
fn from_str_delegates_to_the_reader() {
    let text =
        LammpsDataWriter::write_complete(&POSITIONS, &[1.0; 3], &[1, 1, 1], [0.0; 3], [9.0; 3]);
    let r: LammpsDataReader = text.parse().unwrap();
    assert_eq!(r.positions().len(), 3);
    assert!("not a data file".parse::<LammpsDataReader>().is_err());
}

#[test]
fn file_without_atoms_section_is_an_empty_system() {
    let text = "t\n\n0 atoms\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n";
    let r = LammpsDataReader::from_data_str(text).unwrap();
    assert!(r.positions().is_empty() && r.atom_style().is_none());
}

#[test]
fn huge_declared_counts_do_not_allocate() {
    let text = "t\n\n999999999999 atoms\n999999999999 bonds\n\n0 1 xlo xhi\n0 1 ylo yhi\n0 1 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n";
    let message = error_text(LammpsDataReader::from_data_str(text));
    assert!(message.contains("999999999999"), "{message}");
}
