// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the CHGCAR reader: grid indexing, truncated input and untrusted
//! header sizes.

use super::functions::lattice_volume;
use super::types::{ChgcarFile, MolVizError, PoscarCoordMode, PoscarFile};

/// POSCAR block of a single H atom in a 2 A cubic cell (volume 8 A^3).
const POSCAR_BLOCK: &str =
    "Test cell\n1.0\n2.0 0.0 0.0\n0.0 2.0 0.0\n0.0 0.0 2.0\nH\n1\nDirect\n0.0 0.0 0.0\n";

/// Format `values` five per line, the layout VASP uses.
fn grid_text(values: &[f64]) -> String {
    let mut out = String::new();
    for chunk in values.chunks(5) {
        let row: Vec<String> = chunk.iter().map(|v| format!("{v}")).collect();
        out.push(' ');
        out.push_str(&row.join(" "));
        out.push('\n');
    }
    out
}

fn ramp(n: usize) -> Vec<f64> {
    (0..n).map(|i| i as f64).collect()
}

fn chgcar_text(dims: &str, values: &[f64]) -> String {
    format!("{POSCAR_BLOCK}\n{dims}\n{}", grid_text(values))
}

#[test]
fn density_uses_x_fastest_index() {
    let text = format!(
        "{}augmentation occupancies   1   1\n 0.1 0.2\n",
        chgcar_text("2 3 4", &ramp(24))
    );
    let chg = ChgcarFile::parse(&text).unwrap();
    assert_eq!(chg.grid_dims, [2, 3, 4]);
    assert_eq!(chg.charge.len(), 24);
    let vol = lattice_volume(chg.poscar.lattice);
    assert!((vol - 8.0).abs() < 1e-12);
    // ix is the fastest index, then iy, then iz.
    assert!((chg.density(1, 0, 0).unwrap() * vol - 1.0).abs() < 1e-12);
    assert!((chg.density(0, 1, 0).unwrap() * vol - 2.0).abs() < 1e-12);
    assert!((chg.density(0, 0, 1).unwrap() * vol - 6.0).abs() < 1e-12);
    for iz in 0..4 {
        for iy in 0..3 {
            for ix in 0..2 {
                let want = (ix + 2 * (iy + 3 * iz)) as f64;
                let got = chg.density(ix, iy, iz).unwrap() * vol;
                assert!(
                    (got - want).abs() < 1e-12,
                    "({ix},{iy},{iz}): {got} vs {want}"
                );
            }
        }
    }
}

#[test]
fn density_outside_grid_is_none() {
    let chg = ChgcarFile::parse(&chgcar_text("2 3 4", &ramp(24))).unwrap();
    assert!(chg.density(2, 0, 0).is_none());
    assert!(chg.density(0, 3, 0).is_none());
    assert!(chg.density(0, 0, 4).is_none());
    assert!(chg.density(usize::MAX, usize::MAX, usize::MAX).is_none());
    assert!(chg.density(1, 2, 3).is_some());
}

#[test]
fn density_of_degenerate_cell_is_none() {
    let mut chg = ChgcarFile::parse(&chgcar_text("1 1 2", &[1.0, 2.0])).unwrap();
    chg.poscar.lattice = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(chg.density(0, 0, 0).is_none());
}

#[test]
fn poscar_only_input_is_an_error_not_a_panic() {
    let cases = [
        POSCAR_BLOCK.to_string(),
        format!("{POSCAR_BLOCK}\n"),
        format!("{POSCAR_BLOCK}\n\n\n"),
        String::new(),
        "just a comment\n".to_string(),
    ];
    for text in &cases {
        match ChgcarFile::parse(text) {
            Err(MolVizError::MissingSection { .. }) => {}
            other => panic!("expected MissingSection for {text:?}, got {other:?}"),
        }
    }
}

#[test]
fn truncated_grid_is_a_dimension_mismatch() {
    let err = ChgcarFile::parse(&chgcar_text("2 3 4", &ramp(23))).unwrap_err();
    match err {
        MolVizError::DimensionMismatch { message } => {
            assert!(
                message.contains("24") && message.contains("23"),
                "{message}"
            );
        }
        other => panic!("expected DimensionMismatch, got {other:?}"),
    }
}

#[test]
fn grid_followed_by_blank_line_and_no_augmentation_block_is_truncated() {
    let text = format!("{}\n", chgcar_text("2 3 4", &ramp(10)));
    assert!(matches!(
        ChgcarFile::parse(&text),
        Err(MolVizError::DimensionMismatch { .. })
    ));
}

#[test]
fn reads_exactly_the_expected_values_and_ignores_following_blocks() {
    // A spin-polarised file repeats "nx ny nz" and a second grid right after
    // the first one, without an augmentation block in between.
    let mut text = chgcar_text("2 3 4", &ramp(24));
    text.push_str("   2   3   4\n");
    text.push_str(&grid_text(&[100.0, 101.0, 102.0]));
    let chg = ChgcarFile::parse(&text).unwrap();
    assert_eq!(chg.charge, ramp(24));
}

#[test]
fn surplus_values_on_the_last_line_are_rejected() {
    let mut values = ramp(24);
    values.push(24.0);
    let err = ChgcarFile::parse(&chgcar_text("2 3 4", &values)).unwrap_err();
    assert!(
        matches!(err, MolVizError::DimensionMismatch { .. }),
        "{err:?}"
    );
}

#[test]
fn malformed_dimension_lines_are_rejected() {
    for dims in ["2 3", "2 x 4", "2 3 4 5", "-1 3 4", "2.5 3 4"] {
        let err = ChgcarFile::parse(&chgcar_text(dims, &ramp(24))).unwrap_err();
        assert!(
            matches!(err, MolVizError::ParseError { .. }),
            "{dims}: {err:?}"
        );
    }
    let err = ChgcarFile::parse(&chgcar_text("0 3 4", &ramp(24))).unwrap_err();
    assert!(
        matches!(err, MolVizError::DimensionMismatch { .. }),
        "{err:?}"
    );
}

#[test]
fn hostile_grid_sizes_are_errors_not_allocations() {
    for dims in [
        "18446744073709551615 2 2",
        "4294967296 4294967296 1",
        "100000 100000 100000",
    ] {
        let err = ChgcarFile::parse(&chgcar_text(dims, &ramp(8))).unwrap_err();
        assert!(
            matches!(err, MolVizError::DimensionMismatch { .. }),
            "{dims}: {err:?}"
        );
    }
}

#[test]
fn non_numeric_density_value_is_a_parse_error() {
    let text = format!("{POSCAR_BLOCK}\n1 1 2\n 0.5 abc\n");
    assert!(matches!(
        ChgcarFile::parse(&text),
        Err(MolVizError::ParseError { .. })
    ));
}

#[test]
fn fortran_three_digit_exponents_are_repaired() {
    let text = format!("{POSCAR_BLOCK}\n1 1 3\n 0.5 0.12500000000-100 -0.25000000000+101\n");
    let chg = ChgcarFile::parse(&text).unwrap();
    assert!((chg.charge[0] - 0.5).abs() < 1e-15);
    assert!((chg.charge[1] / 1.25e-101 - 1.0).abs() < 1e-12);
    assert!((chg.charge[2] / -2.5e100 - 1.0).abs() < 1e-12);
}

#[test]
fn vasp4_poscar_without_species_line_is_supported() {
    let text = "Old style\n1.0\n3.0 0.0 0.0\n0.0 3.0 0.0\n0.0 0.0 3.0\n2\nCartesian\n0.0 0.0 0.0\n1.5 1.5 1.5\n\n1 1 2\n 0.25 0.75\n";
    let chg = ChgcarFile::parse(text).unwrap();
    assert!(chg.poscar.species.is_empty());
    assert_eq!(chg.poscar.counts, vec![2]);
    assert_eq!(chg.poscar.coord_mode, PoscarCoordMode::Cartesian);
    assert_eq!(chg.grid_dims, [1, 1, 2]);
    assert_eq!(chg.charge, vec![0.25, 0.75]);
    let vol = lattice_volume(chg.poscar.lattice);
    assert!((chg.density(0, 0, 1).unwrap() * vol - 0.75).abs() < 1e-12);
}

#[test]
fn selective_dynamics_block_does_not_shift_the_grid() {
    let text = "Selective\n1.0\n2.0 0.0 0.0\n0.0 2.0 0.0\n0.0 0.0 2.0\nFe O\n1 1\nSelective dynamics\nDirect\n0.0 0.0 0.0 T T T\n0.5 0.5 0.5 F F F\n\n1 2 1\n 7.0 9.0\n";
    let chg = ChgcarFile::parse(text).unwrap();
    assert!(chg.poscar.selective_dynamics);
    assert_eq!(chg.poscar.n_atoms(), 2);
    assert_eq!(chg.grid_dims, [1, 2, 1]);
    assert_eq!(chg.charge, vec![7.0, 9.0]);
}

#[test]
fn poscar_block_is_preserved() {
    let chg = ChgcarFile::parse(&chgcar_text("1 1 1", &[3.0])).unwrap();
    let direct = PoscarFile::parse(POSCAR_BLOCK).unwrap();
    assert_eq!(chg.poscar.comment, direct.comment);
    assert_eq!(chg.poscar.positions, direct.positions);
    assert_eq!(chg.poscar.counts, direct.counts);
}
