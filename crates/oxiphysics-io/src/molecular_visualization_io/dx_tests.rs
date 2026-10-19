// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the DX volumetric grid reader and writer, with files laid out the
//! way APBS writes them (gridconnections is object 2, the data array object 3).

use super::types::{DxGrid, MolVizError};

/// Format `v` like C's `%.6e` (two-digit signed exponent), as APBS does.
fn c_exp(v: f64) -> String {
    let text = format!("{v:.6e}");
    let (mantissa, exponent) = text.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exponent.abs())
}

/// An APBS-style file for a `counts` grid whose data block holds `values`
/// while the array header declares `declared_items` items.
fn apbs_file(counts: [usize; 3], declared_items: usize, values: &[f64]) -> String {
    let [nx, ny, nz] = counts;
    let mut s = String::new();
    s.push_str("# Data from APBS 1.5\n#\n# POTENTIAL (kT/e)\n#\n");
    s.push_str(&format!(
        "object 1 class gridpositions counts {nx} {ny} {nz}\n"
    ));
    s.push_str("origin -1.000000e+01 -5.000000e+00 2.500000e+00\n");
    s.push_str("delta 5.000000e-01 0.000000e+00 0.000000e+00\n");
    s.push_str("delta 0.000000e+00 5.000000e-01 0.000000e+00\n");
    s.push_str("delta 0.000000e+00 0.000000e+00 5.000000e-01\n");
    s.push_str(&format!(
        "object 2 class gridconnections counts {nx} {ny} {nz}\n"
    ));
    s.push_str(&format!(
        "object 3 class array type double rank 0 items {declared_items} data follows\n"
    ));
    for row in values.chunks(3) {
        let cells: Vec<String> = row.iter().map(|v| c_exp(*v)).collect();
        s.push_str(&cells.join(" "));
        s.push('\n');
    }
    s.push_str("attribute \"dep\" string \"positions\"\n");
    s.push_str("object \"regular positions regular connections\" class field\n");
    s.push_str("component \"positions\" value 1\n");
    s.push_str("component \"connections\" value 2\n");
    s.push_str("component \"data\" value 3\n");
    s
}

fn ramp(n: usize) -> Vec<f64> {
    (0..n).map(|i| i as f64).collect()
}

fn expect_dimension_mismatch(text: &str) -> String {
    match DxGrid::parse(text) {
        Err(MolVizError::DimensionMismatch { message }) => message,
        other => panic!("expected DimensionMismatch, got {other:?}"),
    }
}

#[test]
fn apbs_file_loads_its_data_from_object_three() {
    let text = apbs_file([2, 3, 4], 24, &ramp(24));
    let grid = DxGrid::parse(&text).unwrap();
    assert_eq!((grid.nx, grid.ny, grid.nz), (2, 3, 4));
    assert_eq!(grid.data, ramp(24), "values must not be zeros");
    assert_eq!(grid.origin, [-10.0, -5.0, 2.5]);
    assert_eq!(grid.delta_x, [0.5, 0.0, 0.0]);
    assert_eq!(grid.delta_y, [0.0, 0.5, 0.0]);
    assert_eq!(grid.delta_z, [0.0, 0.0, 0.5]);
    assert!(grid.units.is_empty(), "no comment mentions units");
}

#[test]
fn windows_line_endings_are_accepted() {
    let text = apbs_file([2, 3, 4], 24, &ramp(24)).replace('\n', "\r\n");
    let grid = DxGrid::parse(&text).unwrap();
    assert_eq!(grid.data, ramp(24));
    assert_eq!(grid.origin, [-10.0, -5.0, 2.5]);
    assert_eq!(grid.delta_z, [0.0, 0.0, 0.5]);
}

#[test]
fn grid_indexing_keeps_z_fastest() {
    let grid = DxGrid::parse(&apbs_file([2, 3, 4], 24, &ramp(24))).unwrap();
    // Flat index = (ix * ny + iy) * nz + iz, the last index varying fastest.
    assert_eq!(grid.get(0, 0, 1), 1.0);
    assert_eq!(grid.get(0, 1, 0), 4.0);
    assert_eq!(grid.get(1, 0, 0), 12.0);
    assert_eq!(grid.get(1, 2, 3), 23.0);
    let c = grid.coord(1, 2, 3);
    assert_eq!(c, [-9.5, -4.0, 4.0]);
}

#[test]
fn legacy_layout_with_the_array_as_object_two_is_still_read() {
    let text = "object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class array type double rank 0 items 8 data follows\n1 2 3\n4 5 6\n7 8\nattribute \"dep\" string \"positions\"\n";
    let grid = DxGrid::parse(text).unwrap();
    assert_eq!(grid.data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
}

#[test]
fn quoted_object_names_and_an_end_marker_are_handled() {
    let text = "object \"positions\" class gridpositions counts 1 1 3\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject \"connections\" class gridconnections counts 1 1 3\nobject \"densities\" class array type float category real rank 0 items 3 data follows\n1.5 -2.5e-1 +3\nend\n";
    let grid = DxGrid::parse(text).unwrap();
    assert_eq!(grid.data, vec![1.5, -0.25, 3.0]);
}

#[test]
fn comments_and_blank_lines_inside_the_data_block_are_skipped() {
    let text = "object 1 class gridpositions counts 1 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class gridconnections counts 1 2 2\nobject 3 class array type double rank 0 items 4 data follows\n1 2\n\n# a remark\n3 4\nattribute \"dep\" string \"positions\"\n";
    let grid = DxGrid::parse(text).unwrap();
    assert_eq!(grid.data, vec![1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn the_last_units_comment_is_kept() {
    let text = apbs_file([1, 1, 3], 3, &[1.0, 2.0, 3.0])
        .replace("# POTENTIAL (kT/e)", "# potential in units of kT/e");
    let grid = DxGrid::parse(&text).unwrap();
    assert_eq!(grid.units, "# potential in units of kT/e");
}

#[test]
fn only_the_first_data_array_is_used() {
    let mut text = apbs_file([1, 1, 3], 3, &[1.0, 2.0, 3.0]);
    text.push_str("object 4 class array type double rank 0 items 3 data follows\n9 9 9\n");
    let grid = DxGrid::parse(&text).unwrap();
    assert_eq!(grid.data, vec![1.0, 2.0, 3.0]);
}

#[test]
fn round_trip_reproduces_every_field_exactly() {
    let mut grid = DxGrid::empty(2, 3, 2);
    grid.origin = [-1.0 / 3.0, 2.5e-9, 1.0e10];
    grid.delta_x = [0.1, 0.0, 0.0];
    grid.delta_y = [0.0, 0.2, 1.0 / 7.0];
    grid.delta_z = [1e-300, 0.0, 0.3];
    grid.units = "# potential units kT/e".to_string();
    grid.data = vec![
        1.0 / 3.0,
        -2.5e-7,
        6.02214076e23,
        f64::MIN_POSITIVE,
        -0.0,
        -1.0,
        1e-12,
        123456.789012345,
        f64::MAX,
        -0.1,
        2.0f64.powi(-30),
        std::f64::consts::PI,
    ];
    let text = grid.to_dx_string();
    let back = DxGrid::parse(&text).unwrap();
    assert_eq!((back.nx, back.ny, back.nz), (2, 3, 2));
    assert_eq!(back.origin, grid.origin);
    assert_eq!(back.delta_x, grid.delta_x);
    assert_eq!(back.delta_y, grid.delta_y);
    assert_eq!(back.delta_z, grid.delta_z);
    assert_eq!(back.units, grid.units);
    assert_eq!(back.data.len(), grid.data.len());
    for (a, b) in back.data.iter().zip(&grid.data) {
        assert_eq!(a.to_bits(), b.to_bits(), "{a} vs {b}");
    }
    // Writing what was read gives the same text again.
    assert_eq!(back.to_dx_string(), text);
}

#[test]
fn writer_emits_gridconnections_and_a_complete_field() {
    let mut grid = DxGrid::empty(2, 3, 4);
    grid.data = ramp(24);
    let text = grid.to_dx_string();
    assert!(text.contains("object 1 class gridpositions counts 2 3 4\n"));
    assert!(text.contains("object 2 class gridconnections counts 2 3 4\n"));
    assert!(text.contains("object 3 class array type double rank 0 items 24 data follows\n"));
    assert!(
        !text.contains("object 2 class array"),
        "the array must not be object 2"
    );
    assert!(text.contains("attribute \"dep\" string \"positions\"\n"));
    assert!(text.contains("class field\n"));
    assert!(text.contains("component \"positions\" value 1\n"));
    assert!(text.contains("component \"connections\" value 2\n"));
    assert!(text.contains("component \"data\" value 3\n"));
    // Three values per line and the objects appear in order.
    let first = text.find("object 1 ").unwrap();
    let second = text.find("object 2 ").unwrap();
    let third = text.find("object 3 ").unwrap();
    assert!(first < second && second < third);
    let data_lines = text
        .lines()
        .skip_while(|l| !l.contains("data follows"))
        .skip(1)
        .take_while(|l| !l.starts_with("attribute"))
        .count();
    assert_eq!(data_lines, 8);
}

#[test]
fn writer_output_of_a_ragged_data_vector_is_rejected_when_read_back() {
    let mut grid = DxGrid::empty(2, 2, 2);
    grid.data.pop();
    let message = expect_dimension_mismatch(&grid.to_dx_string());
    assert!(message.contains("7") && message.contains("8"), "{message}");
}

#[test]
fn short_data_block_is_an_error_not_zero_padding() {
    let text = apbs_file([2, 3, 4], 24, &ramp(23));
    let message = expect_dimension_mismatch(&text);
    assert!(
        message.contains("24") && message.contains("23"),
        "{message}"
    );
}

#[test]
fn surplus_values_are_an_error() {
    let text = apbs_file([2, 3, 4], 24, &ramp(25));
    expect_dimension_mismatch(&text);
}

#[test]
fn declared_item_count_must_match_the_grid() {
    // The data block has the right size but the header claims 8 items.
    let text = apbs_file([2, 3, 4], 8, &ramp(24));
    let message = expect_dimension_mismatch(&text);
    assert!(message.contains("8") && message.contains("24"), "{message}");
}

#[test]
fn zero_and_overflowing_counts_are_errors() {
    for counts in [[0, 3, 4], [2, 0, 4], [2, 3, 0]] {
        expect_dimension_mismatch(&apbs_file(counts, 0, &[]));
    }
    let hostile = apbs_file([usize::MAX, 2, 2], 4, &ramp(4));
    expect_dimension_mismatch(&hostile);
    let hostile = apbs_file([usize::MAX / 2 + 1, 2, 1], 4, &ramp(4));
    expect_dimension_mismatch(&hostile);
}

#[test]
fn two_dimensional_counts_are_rejected() {
    let text = "object 1 class gridpositions counts 4 4\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\nobject 2 class array type double rank 0 items 16 data follows\n";
    expect_dimension_mismatch(text);
}

#[test]
fn missing_sections_are_reported() {
    let no_positions = "object 2 class array type double rank 0 items 1 data follows\n1.0\n";
    assert!(matches!(
        DxGrid::parse(no_positions),
        Err(MolVizError::MissingSection { .. })
    ));
    assert!(matches!(
        DxGrid::parse(""),
        Err(MolVizError::MissingSection { .. })
    ));
    let no_array = "object 1 class gridpositions counts 1 1 1\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class gridconnections counts 1 1 1\n";
    match DxGrid::parse(no_array) {
        Err(MolVizError::MissingSection { section }) => assert!(section.contains("data")),
        other => panic!("expected MissingSection, got {other:?}"),
    }
    let no_origin = "object 1 class gridpositions counts 1 1 1\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 class array type double rank 0 items 1 data follows\n1\n";
    match DxGrid::parse(no_origin) {
        Err(MolVizError::MissingSection { section }) => assert!(section.contains("origin")),
        other => panic!("expected MissingSection, got {other:?}"),
    }
    let no_delta = "object 1 class gridpositions counts 1 1 1\norigin 0 0 0\nobject 3 class array type double rank 0 items 1 data follows\n1\n";
    match DxGrid::parse(no_delta) {
        Err(MolVizError::MissingSection { section }) => assert!(section.contains("delta")),
        other => panic!("expected MissingSection, got {other:?}"),
    }
}

#[test]
fn malformed_numbers_are_parse_errors() {
    let bad_value = apbs_file([1, 1, 3], 3, &[1.0, 2.0, 3.0]).replace("2.000000e+00", "two");
    assert!(matches!(
        DxGrid::parse(&bad_value),
        Err(MolVizError::ParseError { .. })
    ));
    let bad_origin = apbs_file([1, 1, 3], 3, &ramp(3)).replace("-1.000000e+01", "x");
    assert!(matches!(
        DxGrid::parse(&bad_origin),
        Err(MolVizError::ParseError { .. })
    ));
    let bad_delta = apbs_file([1, 1, 3], 3, &ramp(3))
        .replace("delta 5.000000e-01 0.000000e+00 0.000000e+00", "delta 1 0");
    assert!(matches!(
        DxGrid::parse(&bad_delta),
        Err(MolVizError::ParseError { .. })
    ));
    let bad_count =
        apbs_file([1, 1, 3], 3, &ramp(3)).replace("counts 1 1 3\norigin", "counts 1 one 3\norigin");
    assert!(matches!(
        DxGrid::parse(&bad_count),
        Err(MolVizError::ParseError { .. })
    ));
}

#[test]
fn binary_and_vector_arrays_are_rejected() {
    let binary = apbs_file([1, 1, 3], 3, &ramp(3)).replace(
        "rank 0 items 3 data follows",
        "rank 0 items 3 msb binary data follows",
    );
    assert!(matches!(
        DxGrid::parse(&binary),
        Err(MolVizError::ParseError { .. })
    ));
    let vectors = apbs_file([1, 1, 3], 3, &ramp(3)).replace(
        "rank 0 items 3 data follows",
        "rank 1 shape 3 items 3 data follows",
    );
    assert!(matches!(
        DxGrid::parse(&vectors),
        Err(MolVizError::ParseError { .. })
    ));
}

#[test]
fn origin_and_delta_lines_outside_the_gridpositions_object_are_ignored() {
    let mut text = apbs_file([1, 1, 3], 3, &ramp(3));
    text.push_str("origin 99 99 99\ndelta 9 9 9\n");
    let grid = DxGrid::parse(&text).unwrap();
    assert_eq!(grid.origin, [-10.0, -5.0, 2.5]);
    assert_eq!(grid.delta_x, [0.5, 0.0, 0.0]);
}
